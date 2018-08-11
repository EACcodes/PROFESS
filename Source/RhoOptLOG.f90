MODULE RhoOptLOG
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoOptLOG
!     |_SUBROUTINE LogTruncatedNewton
!     |_FUNCTION ConstrainLog
!     |_FUNCTION Propogate
!     |_SUBROUTINE LogLineSearch
!     |_FUNCTION PrecondForLog
!
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
! REFERENCES: 
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE RhoOptimizers
  USE OUTPUT, ONLY: outputUnit
  USE CellInfo, ONLY: cell
  USE MathFunctions, ONLY: volume
  USE CellInfo, ONLY: m3G
  USE Timer, ONLY: TimerStart
  USE Timer, ONLY: TimerStop
  USE Timer, ONLY: stopwatch 
  USE Report, ONLY: MinimizerReportFooter ! Subroutine that prints out the results.
  USE Report, ONLY: MinimizerReportHeader ! Prints out the header for the minimization report
  USE Report, ONLY: MinimizerReportSteps  ! Subroutine that prints individual steps
  USE Output, ONLY: WrtOut

  IMPLICIT NONE

CONTAINS

SUBROUTINE LogTruncatedNewton(rhoR, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This is a truncated Newton minimizer, where the equality constraint is
! handled by introducing chi=ln[rho] and minimizing with respect to chi.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   03/31/2008  Subroutine adapted from RhoTNWithLogLine below.  (VLL)
!
!------------------------------------------------------------------------------

  USE KEDF_WGCkernel, ONLY: firstOrderWGC

  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:,:) :: &
    rhoR           ! The spin-density in real space.

  REAL(kind=DP), DIMENSION(:) :: &
    energy         ! the energy table.

                    !>> INTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(SIZE(rhoR,1),SIZE(rhoR,2), SIZE(rhoR,3), &
                           SIZE(rhoR,4)) :: &
    bestDir, &     ! The best approximation to the Newton step as computed by LCG.
    lowDir, &      ! The direction for the line search as it was when res.res was lowest.
    dir, &         ! The running CG direction for the linear CG algrorithm.
    res, &         ! The residual for the linear CG.
    why, &         ! The preconditioned residual for the LCG.
    hessDir, &     ! The product of the local Hessian with the CG direction.
    trialRho, &    ! Tentative density.
    currentPot     ! The total potential for the current step density.

  REAL(kind=DP) :: &
    totWeight, &   ! The integral of the weighing function for convergence.
    numpoints, &   ! The total number of gridpoints in the real space array.
    nele, &        ! The number of electrons, not normalized.
    normPot, &     ! The standard deviation of the potential.
    crit1, &       ! normPot subpart
    crit2, &       ! normPot subpart
    tolPot, &      ! The maximum potential fluctuation below which CG stops.
    tolCG, &       ! Exit criterion for the linear CG algorithm.
    lambdachi, &   ! Lagrange multiplier estimate for constraining the potential.
    eps, &         ! Small step for the LCG Hessian-vector product computation.
    aCG, &         ! Running alpha in the linear CG algorithm.
    reswhy, &      ! Sum of the residual times the preconditioned direction.
    resres, &      ! Sum of the residual squared.
    oldreswhy, &   ! Storage for the previous value of whyres.
    refresres, &   ! Storage for the first value of resres.
    lowresres, &   ! The lowest value of the res.res inner product so far.
    denwolf, &     ! Subpart of wolfe below.
    wolfe, &       ! slope of the line search a dt=0, useful for Wolfe conditions.
    dt             ! current time step for the line minimization.

  REAL(kind=DP), DIMENSION(SIZE(energy)) :: &
    lineEnergy, &  ! The best energy along a given line search.
    startE         ! The best energy at the beginning of a given mu-loop.

  INTEGER :: &
    eCalc, &       ! Counter for the total number of energy computations.
    potCalc, &     ! Running number of potential computations in the current LCG.
    totalPot, &    ! Counter for the total number of potential computations.
    linesteps, &   ! number of energy calculations along a given line search.
    lineOutcome, & ! The reason the line minimization stopped.
    linecare, &    ! How careful the linesearch should be.
    lcgoutcome, &  ! The reason the current LCG stopped.
    tnOutcome, &   ! The reason the Newton algorithm stopped.
    precondType, & ! The kind of preconditioner used, in any.
    i, &           ! Linear CG step counter.
    j, &           ! Newton step counter.
    maxi           ! Maximum number of linear CG steps in a single Newton step.

  LOGICAL :: &
    badDir, &      ! TRUE if the search direction has to be inverted to be descent.
    bailedout      ! Flag to know if the current LCG ended in failure.

  TYPE(stopwatch):: &
    watch          ! Timer

                    !>> INITIALIZATION <<!
  tolPot = 1.0E-6_DP
  tolCG = 1.0E-1_DP
  eps = 1.0E-7_DP

  tnOutcome = -1
  watch = TimerStart()

                    !>> FUNCTION BODY <<!

  IF(rankGlobal==0) THEN
    WRITE(*,*) " Use Truncated Newton method for Log(rho)"
  ENDIF

  ! Report pagination.
  CALL MinimizerReportHeader
  totalPot = 0

  IF (usePreconditioner) THEN
    precondType = 1
  ELSE
    precondType = 0
  END IF

#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(SUM(rhoR), nele, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
  numpoints = REAL(SIZE(rhoR,1)*SIZE(rhoR,2)*m3G*SIZE(rhoR,4),kind=DP)
#else
  nele = SUM(rhoR)
  numpoints = REAL(SIZE(rhoR),kind=DP)
#endif
  CALL CalculatePotentialPlus(rhoR, .TRUE., currentPot, startE)
  ecalc = 1



  ! --------------------------------------
  ! start Newton Loop
  ! --------------------------------------
  WRITE(message, '(A)')'Iter    totEnergy(Ha)'
  CALL WrtOut(6,message)
  DO j=1, maxIter ! This is the Newton loop, each pass is a Newton step.
    potCalc = 0
#ifdef __USE_PARALLEL 
    CALL MPI_ALLREDUCE(SUM(currentPot*rhoR**2), lambdachi, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr) 
    IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***" 
    CALL MPI_ALLREDUCE(SUM(rhoR**2), totWeight, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr) 
    IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***" 
#else 
    lambdachi = SUM(currentPot*rhoR**2) 
    totWeight = SUM(rhoR**2) 
#endif 
    lambdachi = lambdachi / totWeight
    res = -rhoR * (currentPot - lambdachi)
#ifdef __USE_PARALLEL
    CALL MPI_ALLREDUCE(SUM(res**2), crit1, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
    CALL MPI_ALLREDUCE(SUM(res), crit2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
#else
    crit1 = SUM(res**2)
    crit2 = SUM(res)
#endif
    crit2 = crit2 / numpoints
    normPot = SQRT((crit1/numpoints - crit2**2))
    ! The "norm" of the potential is actually its standard deviation over the cell.
    ! If the standard deviation of the potential is low enough, we stop.
    IF (normPot<tolPot) THEN
      CALL MinimizerReportSteps(j, startE, normPot, 0, 0, .FALSE., 1, 0.0_DP)
      tnOutcome = 0
      EXIT
    ELSE IF (normPot<100._DP*tolPot.AND.firstOrderWGC==1) THEN
      firstOrderWGC = -1
      CALL CalculatePotentialPlus(rhoR, .TRUE., currentPot, startE)
      ecalc = ecalc + 1
      CYCLE
    END IF
    IF(j==1) CALL MinimizerReportSteps(0, startE, normPot, 1, 1, .FALSE., 1, 0.0_DP)


    ! Ready to obtain our best Newton step approximation.
    ! Store the reference for the exit criterion.
    refresres = crit1

    ! Initialize the LCG
    bestDir = 0._DP
    lowDir = res
    lowresres = refresres
    dir = 0._DP
    reswhy = 1._DP
    bailedout = .FALSE.
    lcgoutcome = -1

    IF (j<5) THEN
      maxi = 5
    ELSE
      maxi = 100
    END IF

    !-----------------------------------------
    ! start line search loop
    !-----------------------------------------
    DO i=1, maxi

      ! Preconditioning.
      why = PrecondForLog(rhoR, res, precondType)
      oldreswhy = reswhy
#ifdef __USE_PARALLEL
      CALL MPI_ALLREDUCE(SUM(res*why), reswhy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
      IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
#else
      reswhy = SUM(res*why)
#endif
      dir = reswhy/oldreswhy * dir + why

      ! Compute the Hessian-vector product by difference.
      trialRho = rhoR * EXP(eps * dir)
      CALL CalculatePotentialPlus(trialRho, .TRUE., hessDir, startE)

      hessDir = (hessDir - currentPot) / eps + (currentPot - lambdachi) * dir
      hessDir = hessDir * rhoR

#ifdef __USE_PARALLEL
      CALL MPI_ALLREDUCE(SUM(hessDir * dir), aCG, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
      IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
#else
      aCG = SUM(hessDir * dir)
#endif
      aCG = reswhy / aCG
      IF (aCG<0) THEN
        lcgoutcome = 2
        bailedout = .TRUE.
        IF(i==1) bestDir = res
        EXIT
      END IF

      bestDir = bestDir + aCG * dir
      res = res - aCG * hessDir

#ifdef __USE_PARALLEL
      IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
      CALL MPI_ALLREDUCE(SUM(res**2), resres, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
      IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
#else
      resres = SUM(res**2)
#endif
      IF (resres<lowresres) THEN
        lowresres = resres
        lowDir = bestDir
      END IF
      ! Have we made enough progress that we can exit yet?
      IF ((resres/refresres)<tolCG) THEN
        ! We have beaten our convergence criterion. Get out.
        lcgoutcome = 0
        EXIT
      END IF
    END DO ! Linear CG loop.


    ! We have exited the LCG loop. lcgoutcome tells us why:
    ! lcgoutcome = 0 means LCG reached convergence.
    ! lcgoutcome = 2 means the LCG found a non-positive definite Hessian.
    ! lcgoutcome =-1 means LCG ran for maxi steps without converging.
    ! bailedout = true is redundant with lcgoutcome = 2.
    linecare = ABS(lcgoutcome)
    potCalc = i
    totalPot = totalPot + potCalc

    ! Whatever the outcome, we will use the best direction for the line search.
    ! We need to make sure that this is a descent direction.
#ifdef __USE_PARALLEL
    CALL MPI_ALLREDUCE(SUM(rhoR*bestDir), wolfe, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
    CALL MPI_ALLREDUCE(SUM(rhoR**2), denwolf, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
#else
    wolfe = SUM(rhoR*bestDir)
    denwolf = SUM(rhoR**2)
#endif
    wolfe = wolfe / denwolf
#ifdef __USE_PARALLEL
    CALL MPI_ALLREDUCE(SUM(currentPot*rhoR*(bestDir-rhoR*wolfe)), wolfe, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
#else
    wolfe = SUM(currentPot*rhoR*(bestDir-rhoR*wolfe))
#endif
    IF (wolfe>0._DP) THEN
      ! Our best direction is NOT a descent direction! Panic.
      IF (rankGlobal==0) WRITE(*,*) "Panic in LOGTRUNCATEDNEWTON! SEND HELP!"
      bestDir = -bestDir
      wolfe = -wolfe
      badDir = .TRUE.
    ELSE
      badDir = .FALSE.
    END IF

    ! Proceed with the line search. First try the predicted minimum.
    dt = 1._DP
    CALL LogLineSearch(rhoR, bestDir, startE(1), wolfe, dt, trialRho, currentPot, lineEnergy, lineOutcome, lineSteps, linecare)
    eCalc = eCalc + linesteps
    CALL MinimizerReportSteps(j, lineEnergy, normPot, potCalc, lineSteps, badDir, lcgoutcome, dt)
    SELECT CASE (lineOutcome)
    CASE (0) ! Minimum found in line.
      rhoR = trialRho
      startE = lineEnergy
    CASE (2) ! Line search has failed.
      ! We are stuck. Exit, and say we had a problem.
      tnOutcome = 2
      EXIT
    CASE DEFAULT ! Ran out of TN steps.
      tnOutcome = 1
      EXIT
    END SELECT

    WRITE(message,'(I5,ES20.12)') j, startE(1)
    CALL WrtOut(6, message)
  END DO ! Newton loop.
  energy = startE
  CALL MinimizerReportFooter(j, tnOutcome, energy, eCalc, totalPot, TimerStop(watch))
  
  ! mohan add 2012-12-20
  WRITE(message,'(A,Es19.12,A,Es11.4,A,F10.1)') & 
    & ' Final total energy: ', energy(1), &
    & ' Ha, volume=', volume(cell%cellReal), ' bohr^3' 

  CALL WrtOut(6,message)

  RETURN

END SUBROUTINE LogTruncatedNewton


FUNCTION ConstrainLog(nRho, oRho, numE, dt)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This subroutine takes the proposed step from the line search and constrains
! it so that the corresponding density integrates to the number of electrons.
! The constraining is done iteratively.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   04/21/2008  Subroutine created. (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:,:) :: &
    oRho, &      ! The spin-density that we start from.
    nRho         ! The density that needs constrained.

  REAL(kind=DP) :: &
    numE, &      ! The number of electrons, raw.
    dt, &        ! The step size used, so we make sure it's large enough.
    ConstrainLog ! The constant that makes it all happen.

                    !>> INTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(SIZE(oRho,1),SIZE(oRho,2),SIZE(oRho,3),SIZE(oRho,4)) :: &
    array        ! Temporary array to save on exponentiation costs.

  REAL(kind=DP) :: &
    res, &       ! The difference between the integrated constrained rho and N
    dres, &      ! The first derivative of res with repect to mu.
    mu, &        ! The constraining constant to be determined iteratively.
    dmu          ! Proposed change in mu.

                    !>> INITIALIZATION <<!
  mu = 0._DP
                    !>> FUNCTION BODY <<!
  IF (dt<1.0E-10_DP) THEN
    ! dt is too small for the iterative process. This direct calculation may still
    ! return garbage due to dramatic cancellation. In general, small dt is bad news.
#ifdef __USE_PARALLEL
    CALL MPI_ALLREDUCE(SUM(nRho), res, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP &
      "***PROBLEMS WITH MPI_ALLREDUCE IN CONSTRAINLOG***"
    CALL MPI_ALLREDUCE(SUM(oRho), dres, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP &
      "***PROBLEMS WITH MPI_ALLREDUCE IN CONSTRAINLOG***"
#else
    res = SUM(nRho)
    dres = SUM(oRho)
#endif
    mu = (res-numE) / dres
  ELSE
    DO
      array = nRho * EXP(-mu*oRho)
#ifdef __USE_PARALLEL
    CALL MPI_ALLREDUCE(SUM(array), res, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP &
      "***PROBLEMS WITH MPI_ALLREDUCE IN CONSTRAINLOG***"
    CALL MPI_ALLREDUCE(SUM(array*oRho), dres, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) STOP &
      "***PROBLEMS WITH MPI_ALLREDUCE IN CONSTRAINLOG***"
#else
    res = SUM(array)
    dres = SUM(array*oRho)
#endif
      res = res - numE
      dmu = res / dres
      mu = mu + dmu
      IF (ABS(dmu/mu)<1.0E-10_DP .OR. ABS(res/numE)<1.0E-10_DP) EXIT
    END DO
  END IF

  ConstrainLog = mu

END FUNCTION ConstrainLog




FUNCTION Propagate(dir, den)
!------------------------------------------------------------------------------ 
! DESCRIPTION: 
! This functions transforms a search direction for the density into a search
! direction for its log, while keeping the direction identical for an infinitesimal
! step.
! 
! GLOBAL/MODULE VARIABLES CHANGED: 
! 
! CONDITIONS AND ASSUMPTIONS: 
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
! 
! REFERENCES: 
! 
!------------------------------------------------------------------------------
! REVISION LOG: 
!   01/09/2008  Function created.  (VLL)
! 
!------------------------------------------------------------------------------
  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: &
    dir, &       ! The search direction to be tranformed.
    den          ! The reference electron density.

  REAL(kind=DP), DIMENSION(0:SIZE(dir,1)-1, 0:SIZE(dir,2)-1, &
                           0:SIZE(dir,3)-1, SIZE(dir,4)) :: &
    Propagate    ! The modified direction.

                    !>> INTERNAL VARIABLES <<! 
  REAL(kind=DP) :: &
    cutoff, &    ! Lowest value propagator can take.
    lowbound     ! Lowest value of propagator in this case.

                    !>> INITIALIZATION <<! 
cutoff = -2._DP
                    !>> FUNCTION BODY <<! 
  Propagate = 0._DP
  WHERE (dir>-den)
    propagate = LOG(1._DP+dir/den)
!ELSEWHERE
!propagate = dir/den
!END WHERE

  END WHERE
  lowbound = MINVAL(Propagate)
  lowbound = MIN(lowbound, -1._DP)
  WHERE (dir<(-den))
    propagate = lowbound
  END WHERE
  WHERE (propagate<cutoff)
    propagate = cutoff
  END WHERE

END FUNCTION Propagate


SUBROUTINE LogLineSearch(rho, line, phi0, dphi0, step, newrho, newpot, energy, outcome, totalsteps, care)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This subroutine takes the current log(rho) and the search direction to find
! the step that minimizes the energy along that line. It returns it, as well
! as the resulting energy.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   03/26/2008  Subroutine adapted from LineSearch.  (VLL)
!
!------------------------------------------------------------------------------
  USE CellInfo, ONLY : cell
  USE MATHFUNCTIONS, ONLY : volume

  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: &
    rho, &       ! The spin-density in real space.
    line         ! The minimization direction.

  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(OUT) :: &
    newrho, &    ! Density at line minimum.
    newpot       ! Potential at line minimum.

  REAL(kind=DP), INTENT(INOUT) :: &
    step, &      ! The step that yields the energy minimum along the line.
    phi0, &      ! Initial energy.
    dphi0        ! Initial slope to be normalized.

  REAL(kind=DP), DIMENSION(:), INTENT(OUT) :: &
    energy       ! the energy table.

  INTEGER, INTENT(IN) :: &
    care         ! How careful we should be, related to how the LCG ended.

  INTEGER, INTENT(OUT) :: &
    outcome, &   ! Reason for exiting this algorithm
    totalsteps   ! Total number of line minimization steps

                    !>> INTERNAL VARIABLES <<!
  REAL(kind=DP) :: &
    a1, a2, a3, pt1, pt2, slp1, slp2, num, denom, dV, debugstart, debugstep, &
    delta, tol1, tol6, cap, linemu, slpavg, &
    nele         ! The number of electrons, not normalized.

  REAL(kind=DP) :: &
    c1, &        ! Wolfe condition on value parameter.
    c2           ! Wolfe condition on slope parameter.

  REAL(kind=DP), DIMENSION(3) :: &
    amin, &      ! The four predicted extrema position in the parabola fits.
    curve        ! 'a' in phi(x)= a(x-a1)^2 + b(x-a1) + c

  LOGICAL :: &
    debugmode, careful

  REAL(kind=DP) :: &
    numpoints    ! Total number of gridpoints.
    
  INTEGER :: &
    maxk, &      ! Maximum number of line search steps.
    k            ! Line step counter.

                    !>> INITIALIZATION <<!
  ! DEBUG
  debugmode = .FALSE.
  debugstart = 0.0_DP
  debugstep = 0.02_DP

  IF (care>0) THEN
    careful = .TRUE.
  ELSE
    careful = .FALSE.
  END IF

  c1 = 1.0E-4_DP
  c2 = 0.9_DP
  outcome = 2
  tol1 = 0.8_DP
  tol6 = 0.7_DP
  cap = 1.0E20_DP
  maxk = 100
#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(SUM(rho), nele, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, &
                     mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN RHOTNWITHLOGLINE***"
  numpoints = REAL(SIZE(rho,1)*SIZE(rho,2)*m3G*SIZE(rho,4),kind=DP)
#else
  nele = SUM(rho)
  numpoints = REAL(SIZE(rho),kind=DP)
#endif

  dV = cell%dV

  a1 = 0._DP
  a2 = step
  pt1 = phi0
  slp1 = dphi0 * dV
  dphi0 = slp1
  totalsteps = 0

!  IF(rankGlobal==0) THEN
!  WRITE(*,*) "Starting a new line."
!  WRITE(*,*) "Step                   Energy                  Slope"
!  WRITE(*,*) a1, pt1, slp1
!  END IF

                    !>> FUNCTION BODY <<!
  DO k=1, maxk
    ! This line is for debugging purposes only.
    IF (debugmode) a2 = debugstart + k * debugstep

    newrho = rho * EXP(a2 * line)
    IF (MINVAL(newrho)<=0._DP.OR.MAXVAL(newrho)>HUGE(0._DP)) THEN
      a2 = 0.5 * (a1+a2)
      CYCLE
    END IF
    linemu = ConstrainLog(newrho, a2*rho, nele, a2)
    IF (linemu.NE.linemu) THEN
      a2 = 0.5 * (a1+a2)
      CYCLE
    END IF
    newrho = newrho * EXP(-a2*linemu*rho)
    IF (MINVAL(newrho)<=0._DP.OR.MAXVAL(newrho)>HUGE(0._DP)) THEN
      a2 = 0.5 * (a1+a2)
      CYCLE
    END IF
    CALL CalculatePotentialPlus(newrho, .TRUE., newpot, energy)
    IF (energy(1).NE.energy(1)) THEN
      a2 = 0.5 * (a1+a2)
      CYCLE
    END IF
    totalsteps = totalsteps + 1
    pt2 = energy(1)

    ! Now compute the slope at a2.
#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(SUM(newrho*line), num, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN LINESEARCH***"
  CALL MPI_ALLREDUCE(SUM(newrho*rho), denom, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN LINESEARCH***"
#else
  num = SUM(newrho*line)
  denom = SUM(newrho*rho)
#endif
    num = num / denom
#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(SUM(newrho*newpot*(line-num*rho)), slp2, 1, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN LINESEARCH***"
#else
  slp2 = SUM(newrho*newpot*(line-num*rho))
#endif
    slp2 = slp2 * dV
!  IF(rankGlobal==0) WRITE(*,*) a2, pt2, slp2

    IF (debugmode) CYCLE

    ! Now we have our two points with slopes.
    ! Check the strong Wolfe conditions
    IF (pt2<=phi0+c1*a2*dphi0 .AND. ABS(slp2)<=c2*ABS(dphi0)) THEN
      ! Wolfe conditions satisfied. We're done here.
      outcome = 0
      EXIT
    END IF

    ! Based on this information, estimate the shape of the energy profile between a1 and a2.
    delta = a2 - a1
    ! Parabola 1 is made of pt1, pt2 and slp1
    curve(1) = (pt2 - pt1)/delta**2 - slp1/delta
    amin(1) = -0.5_DP * slp1 / curve(1)
    ! Parabola 2 is made of slp1 and slp2
    curve(2) = 0.5_DP * (slp2 - slp1)/delta
    amin(2) = -0.5_DP * slp1 / curve(2)
    ! Parabola 3 is made of pt1, pt2 and slp2
    curve(3) = slp2/delta + (pt1 - pt2)/delta**2
    amin(3) = (0.5_DP * slp2 - (pt2 - pt1)/delta) / curve(3)
    slpavg = (pt2 - pt1)/delta

    IF (pt2>pt1) THEN
      ! Energy increased, this is a bracket case.
      cap = MIN(cap, a2)
      a2 = a1 + amin(1)
      CYCLE

    ELSE IF (slp2>0._DP) THEN
      ! We passed a min, another bracket case.
      cap = MIN(cap, a2)
      a2 = a1 + amin(2)
      CYCLE

    ELSE IF(pt2>pt1+delta*slpavg) THEN
      ! parabola 1 predicts a min, para 2 may as well.
      a3 = a1+amin(1)
      IF (amin(2)>0._DP.AND.curve(2)>0._DP) a3 = MIN(a3, a1+amin(2))
      IF (a3>=cap) a3 = MIN(a2+10._DP*delta,0.5_DP*(a2+cap))
      ! And if that's still not good just look halfway between the bounds.
      IF (a3>=cap) a3 = 0.5_DP*(a1+cap)
      IF (a3>a2) THEN
        a1 = a2
        pt1 = pt2
        slp1 = slp2
      END IF
      a2 = a3

    ELSE IF (slp2>slp1) THEN
      ! parabola 2 predicts a minimum, let's see...
      a3 = a1 + amin(2)
      IF (a3>=cap) a3 = MIN(a2+10._DP*delta,0.5_DP*(a2+cap))
      ! And if that's still not good just look halfway between the bounds.
      IF (a3>=cap) a3 = 0.5_DP*(a1+cap)
      IF (a3>a2) THEN
        a1 = a2
        pt1 = pt2
        slp1 = slp2
      END IF
      a2 = a3

    ELSE
      ! Most probleamtic case, can't use parabolae 1 or 2 at all.
      IF (amin(3)>0._DP.AND.curve(3)>0._DP) THEN
        ! Parabola 3 might work. go for it.
        a3 = a1 + amin(3)
        IF (a3>=cap) a3 = MIN(a2+10._DP*delta,0.5_DP*(a2+cap))
        ! And if that's still not good just look halfway between the bounds.
        IF (a3>=cap) a3 = 0.5_DP*(a1+cap)
      IF (a3>a2) THEN
        a1 = a2
        pt1 = pt2
        slp1 = slp2
      END IF
        a2 = a3
        careful = .TRUE.
      ELSE
        ! There should be no minimum between a1 and a2.
        a3 = a2 + 10._DP * delta
        IF (a3>cap) a3 = 0.5_DP * (cap+a2)
        a1 = a2
        pt1 = pt2
        slp1 = slp2
        a2 = a3
        careful = .TRUE.
      END IF
    END IF

    ! Make sure the exit criteria are current.
    IF (careful) THEN
      c1 = 0.2_DP
      c2 = 0.2_DP
      IF (care==2) c2 = 0.001_DP
    END IF

  END DO

  ! We are out of the line search loop because:
  ! outcome = 0 : We found a stationary point.
  ! outcome = 2 : We ran out of step.
  step = a2

END SUBROUTINE LogLineSearch


FUNCTION PrecondForLog(rho, p, which)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the preconditioned direction for the TN algorithm with rho.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! We are assuming that the most ill-conditioned energy term is the von Weizsacker.
! If there is no von Weizsacker (mu=0) this function will not do anything.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/29/2007 Function created.  (VLL)
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: k1G, k2G, k3G
  USE CellInfo, ONLY: n1G, n2G, n3G
  USE CellInfo, ONLY: numSpin
  USE PlaneWave, ONLY: qTable
  USE KEDF_VW, ONLY: mu
  USE Fourier, ONLY: FFT

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: p
  ! residual in the LCG.
  !
  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  ! Electron density in real space.
  !
  REAL(kind=DP), DIMENSION(n1G, n2G, n3G, numSpin) :: PrecondForLog
  ! The direction after application on the preconditioner.
  !
  INTEGER :: which
  ! Type of preconditioning that needs doing.

                          !>> INTERNAL VARIABLES <<!

  COMPLEX(kind=DP), DIMENSION(k1G, k2G, k3G, numSpin) :: recipArray
  ! Array for reciprocal space business.
  !

                           !>> INITIALIZATION <<!
                           !>> FUNCTION BODY <<!

  IF (mu == 0._DP) THEN
    ! No von Weizsacker term. We cannot precondition.
    PrecondForLog = p
  ELSE
    SELECT CASE(which)
    CASE(1)
      ! Von Weizsacker solver for log. Usually good.
      recipArray = FFT(p / SQRT(rho))
      recipArray = recipArray / SPREAD(qTable, DIM=4, NCOPIES=SIZE(rho,4))**2
      recipArray(1,1,1,1) = (0._DP,0._DP)
      PrecondForLog = (4._DP / mu) * FFT(recipArray)
      precondForLog = precondForLog / SQRT(rho)

    CASE DEFAULT
      ! Not accounted for. Do nothing.
      precondForLog = p

    END SELECT
  END IF
END FUNCTION PrecondForLog


END MODULE RhoOptLOG

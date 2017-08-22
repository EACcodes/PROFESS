MODULE RhoOptSTN
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoOptSTN
!     |_SUBROUTINE SqrtNewtonMinimization
!       |_FUNCTION Constrain
!       |_FUNCTION NewtonDirection
!       |_FUNCTION Step
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
! 2013 mohan separate this file from big file.
!
!------------------------------------------------------------------------------

  USE RhoOptimizers
  USE OUTPUT, ONLY: outputUnit
  USE Output, ONLY: WrtOut
  USE CellInfo, ONLY: cell
  USE MPI_Functions, ONLY: rankGlobal
  USE Timer, ONLY: TimerStart
  USE Timer, ONLY: TimerStop
  USE Timer, ONLY: stopwatch 

  USE Report, ONLY: MinimizerReportFooter ! Subroutine that prints out the results.
  USE Report, ONLY: MinimizerReportHeader ! Prints out the header for the minimization report
  USE Report, ONLY: MinimizerReportSteps  ! Subroutine that prints individual steps

  IMPLICIT NONE

  REAL(KIND=DP) :: sumRho = 0.0_DP

CONTAINS


SUBROUTINE SqrtNewtonMinimization(energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine minimizes the energy w.r.t. to the sqrt of the electron
!   density using the Truncated Newton method.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE KEDF_WGCkernel, ONLY: firstOrderWGC
  USE Sys, ONLY: rhoR 
  USE SYS, ONLY: interior

  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:), INTENT(OUT) :: energy 
  ! where the final energies are stored   
  !

                    !>> INTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(SIZE(energy)) :: tenEnergiesAgo

  REAL(KIND=DP), DIMENSION(LBOUND(rhoR,1):UBOUND(rhoR,1), &
                           LBOUND(rhoR,2):UBOUND(rhoR,2), &
                           LBOUND(rhoR,3):UBOUND(rhoR,3), &
                           SIZE(rhoR,4)) :: &
    dirNext, &    ! The next direction that we would like to search in,
                  ! can be either the conjugate or steepest descent direction.
    pot, &        ! The total potential for the current step density.
    tempRA

  REAL(KIND=DP) :: &
#ifdef __USE_PARALLEL
    mpiLocalSum, &
    mpiGlobalSum, &
#endif
    potentialNorm, & ! The SQRT of the sum of the potential
    lambda

  INTEGER :: &
#ifdef __USE_PARALLEL
    mpiLocalIntSum, &
#endif
    i, &          ! Minimization step counter
    newtonIter, & ! Number of CGMIN loops in the newton step
    success, &    ! whether our line minimization succeeded
                  !  =1 if it did succed
                  !  =2 if it exited due to a warning
                  !  =3 if it exited because of error in line minimizer
    numPoints, &
    allocatestatus

  LOGICAL :: &
    restart, &    ! whether the present step restarted from the steepest
                  ! descent vector if we were doing conjugate gradient.
    cheat

  REAL(KIND=DP), ALLOCATABLE :: tempRho(:,:,:,:)


  ! Variables for the line search algorithm:
  CHARACTER(len=60) :: task ! Status of line minimizer 

  REAL(KIND=DP) :: &
    stp, &
    ftol = 1.E-4_DP, &
    gtol = 0.6_DP, &
    xtol = 1.E-12_DP, &
    stpmin = 1.E-20_DP, &
    stpmax = 100._DP, &
    e(3)

  INTEGER,DIMENSION(2) :: isave
  REAL(KIND=DP), DIMENSION(13) :: dsave
  TYPE(stopwatch):: watch        ! Timer
  EXTERNAL DCSRCH


                      !>> INITIALIZATION <<!

  CALL WrtOut(6,' RhoOptimizing ... (SqrtNewtonMinimization in RhoOptSTN.f90)')

  watch = TimerStart()
  ALLOCATE(tempRho(LBOUND(rhoR,1):UBOUND(rhoR,1), &
                   LBOUND(rhoR,2):UBOUND(rhoR,2), &
                   LBOUND(rhoR,3):UBOUND(rhoR,3), &
                   SIZE(rhoR,4)), STAT=allocatestatus)
  IF (allocatestatus>0) THEN
    WRITE(*,*) "ALLOCATE DID NOT WORK WITH ERROR TYPE", allocatestatus
    STOP
  END IF

  ! In order for the method to work with taking points out, the tolerance
  ! must be sufficiently small.  That way, when we do "GetActive", the 
  ! energy will lower.
  rhoOutcome = -1
  success = 1
  restart = .FALSE.
  cheat = .TRUE. .AND. cheatPot

                           ! Cheat the inner CG loop by not actually computing
                           ! some expensive potentials

  ! This just sets the sum of the density so that functions like Step can 
  ! use it.
#ifdef __USE_PARALLEL
  mpiLocalSum = SUM(rhoR, interior)
  CALL MPI_ALLREDUCE(mpiLocalSum, sumRho, 1, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpierr /= MPI_SUCCESS) THEN
    CALL WrtOut(6, "***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
    STOP
  ENDIF
#else
  sumRho = SUM(rhoR, interior)
#endif

  ! First compute the starting energy, and count it as an evaluation
  CALL MinimizerReportHeader

  ! Calculate the initial energy
  CALL CalculatePotentialPlus(rhoR, .TRUE., pot, energy)
  numPoints = SIZE(rhoR)

#ifdef __USE_PARALLEL
  mpiLocalIntSum = numPoints
  CALL MPI_ALLREDUCE(mpiLocalIntSum, numPoints, 1, MPI_INTEGER, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) THEN
    CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
    STOP
  ENDIF
#endif

  ! We're going to do something messy here to conserve memory.  We're going to
  ! use the storage of rhoR to store the square root of rho temporarily
  ! for the purposes of this calculation.
  rhoR= SQRT(rhoR)

  pot = Constrain(pot, rhoR, interior, lambda,size(pot,1),size(pot,2),size(pot,3),size(pot,4))
  numEnergyTotal = 1
  numPotentialTotal = 1

  ! Calculate the norm of the potential
  potentialNorm = SUM(pot**2)

#ifdef __USE_PARALLEL
  mpiLocalSum = potentialNorm
  CALL MPI_ALLREDUCE(mpiLocalSum, potentialNorm, 1, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpierr /= MPI_SUCCESS) THEN
    CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
    STOP
  ENDIF
#endif

  potentialNorm = SQRT(potentialNorm/REAL(numPoints,KIND=DP))

  i = 0
  CALL MinimizerReportSteps(i, energy, potentialNorm, 0, 0, .FALSE., 1, 0._DP)

  ! Set up e such that it isn't considered converged in the first step
  e(2) = energy(1) + 2._DP*tole
  e(1) = energy(1) + 4._DP*tole


                    !>> FUNCTION BODY <<!

  WRITE(message, *) '(RhoOptimizers): NEW method'
  CALL WrtOut(6, message)
  WRITE(message, '(A)') ' Iter    totEnergy(Ha)    potNorm     deltaE' 
  CALL WrtOut(6, message)

  ! Loop until we find a minimum!
  DO

    e(3) = e(2); 
    e(2) = e(1); 
    e(1) = energy(1);

    IF (ABS(e(1)-e(2))<tole .AND. ABS(e(1)-e(3))<tole) THEN
      rhoOutcome = 0
    ENDIF



    IF (rhoOutcome==0) THEN
      EXIT
!    ! If constrained potential is tiny, we are close to a critical point, stop.
!    IF (potentialNorm < tol) THEN
!      rhoOutcome = 0
!      EXIT        
    ELSE IF (potentialNorm < 100._DP*pot_tol .AND. firstOrderWGC==1) THEN 
      firstOrderWGC = -1
      CALL CalculatePotentialPlus(rhoR, .TRUE., pot, energy)
      numEnergyTotal = numEnergyTotal + 1
      numPotentialTotal = numPotentialTotal + 1
      CYCLE
    END IF


    ! Here we exit if we took more than maxIter iterations.  
    IF (i > maxIter) THEN 
      rhoOutcome = 1
      EXIT
    END IF

    i = i + 1

    !--------------------------------------------------------------------------------------------
    ! If we are already doing steepest descent and we have a catastropic error
    ! in the line search, there's not much more we can do
    IF((success == 3 .OR. success ==2) .AND. restart .OR. &
       (i > 20 .AND. MOD(i,10)==8 .AND. &
        ABS((tenEnergiesAgo(1)-energy(1))/energy(1)) < 1.E-8_DP)) THEN
      WRITE(*,*) "FATAL ERROR: CATASTROPHIC FAILURE IN LINE MINIMIZER."
      rhoOutcome = 2
      EXIT

    ! Switch to steepest descent if there is an error using truncated newton
    ELSE IF ((success == 3 .OR. success == 2) .AND. .NOT. cheat) THEN
      restart = .TRUE.
      dirNext = - pot

    ! Normal case: Truncated newton    
    ELSE

      ! Let's not terminate the cheating too early, otherwise we miss out
      ! on the benefit!
      IF (success /= 1 .AND. cheat .AND. i > 3) THEN
        cheat = .FALSE.
        WRITE(*,*) "Turning off cheat for finding Newton direction at iteration ", i

      END IF

      restart = .FALSE.

      ! Here we cheat the thing by not calcualating WT and WGC potential terms
      ! in this loop. Still seems to work though, amazingly. It costs another
      ! array to do this. We could also not carry around an extra array 
      ! and then recompute it with the right value at the end. We can 
      ! eliminate this array if it becomes a problem.
      IF(cheat) THEN
        CALL CalculatePotentialPlus(rhoR**2, .TRUE., tempRA)
        numPotentialTotal = numPotentialTotal + 1
      ELSE
        tempRA = pot+lambda*rhoR
      END IF


      dirNext = NewtonDirection(pot, tempRA, lambda, newtonIter)

      numPotentialTotal = numPotentialTotal + newtonIter

    END IF
    !--------------------------------------------------------------------------------------------

    ! Now that we have our next direction in dirNext, go down the line and 
    ! minimize.  
    numEnergyLineMin = 0
    task = 'START'
    stp = 1.0_DP    ! Current initial estimate of satisfactory step
    DO

      !oldtheta = theta ! (Shin, mohan add 10-01-12)
    
#ifdef __USE_PARALLEL
      mpiLocalSum = SUM(pot*dirNext, interior)
      CALL MPI_ALLREDUCE(mpiLocalSum, mpiGlobalSum, 1, MPI_REAL8, MPI_SUM, &
                         MPI_COMM_WORLD, mpiErr)
      IF (mpierr /= MPI_SUCCESS) THEN
        CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
        STOP
      ENDIF

      CALL Dcsrch(stp, energy(1), mpiGlobalSum,&
                  ftol, gtol, xtol, task, stpmin, stpmax, &
                  isave, dsave)
#else
      CALL Dcsrch(stp, energy(1), SUM(pot*dirNext, interior),&
                  ftol, gtol, xtol, task, stpmin, stpmax, &
                  isave, dsave)
#endif

      IF(task .EQ. 'CONVERGENCE') THEN
        EXIT
      ELSE IF(.NOT. (task .EQ. 'FG')) THEN 
        stp = 0._DP
      END IF

      IF(MOD(i,10)==9) tenEnergiesAgo=energy

      ! Update energy and potential
      tempRho = Step(rhoR, stp, dirNext, interior)
      CALL CalculatePotentialPlus(tempRho, .TRUE., pot, energy)
      tempRA = SQRT(tempRho)
      pot = Constrain(pot, tempRA, interior, lambda, size(pot,1), size(pot,2), size(pot,3), size(pot,4) ) 
      numEnergyLineMin = numEnergyLineMin+1

      IF(.NOT. (task .EQ. 'FG')) EXIT

    END DO

    ! What is the status of the line minimizer?
    IF(task(1:5) .EQ. 'ERROR') THEN
      success = 3
      WRITE(*,*) i, ": ", task
    ELSE IF(.NOT. (task .EQ. 'CONVERGENCE')) THEN
      success = 2
      WRITE(*,*) i, ": ", task
    ELSE
      success = 1
      rhoR = tempRA
    END IF

    numEnergyTotal = numEnergyTotal + numEnergyLineMin
    numPotentialTotal = numPotentialTotal + numEnergyLineMin

    ! Calculate the norm of the potential
    potentialNorm = SUM(pot**2)

#ifdef __USE_PARALLEL
    mpiLocalSum = potentialNorm
    CALL MPI_ALLREDUCE(mpiLocalSum, potentialNorm, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpierr /= MPI_SUCCESS) THEN
      CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
      STOP
    ENDIF
#endif

    potentialNorm = SQRT(potentialNorm/REAL(numPoints,KIND=DP))

    CALL MinimizerReportSteps(i, energy, potentialNorm, newtonIter, numEnergyLineMin, restart, success, stp)

    ! output information
    IF(rankGlobal==0) THEN
      WRITE(*, '(I3,ES20.12,ES12.4,ES12.4)') i , energy(1), potentialNorm, e(1)-e(2)
    ENDIF


  END DO

  CALL MinimizerReportFooter(i, rhoOutcome, energy, numEnergyTotal, numPotentialTotal, TimerStop(watch))

  ! Make sure to end up with the actual rho, not the square root of rho
  rhoR = rhoR**2

  ! mohan add 2012-12-20
  WRITE(message,'(A,Es19.12,A,Es11.4,A,F10.1)') & 
    & ' Final total energy: ', energy(1), &
    & ' Ha, volume=', cell%vol, ' bohr^3' 

  CALL WrtOut(6,message)

  DEALLOCATE(tempRho)

  RETURN

CONTAINS


FUNCTION Constrain(step, sqrtRho, mask, mu, dimX, dimY, dimZ, nspin)
!----------------------------------------------------------------------------
! DESCRIPTION:
! This function takes a potential and essentially makes the average of
! its components zero by subtracting an offset.  It's essentially a lagrange
! multiplier that makes the potential norm-conserving.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------

  IMPLICIT NONE
                         !>> EXTERNAL VARIABLES <<!

  INTEGER, INTENT(IN) :: dimX, dimY, dimZ, nspin

  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin), INTENT(IN) :: step
  ! The step we want to make norm-conserving
  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin), INTENT(IN) :: sqrtRho

  LOGICAL, DIMENSION(dimX,dimY,dimZ,nspin), INTENT(IN) :: mask

  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: Constrain
  ! Our answer

  REAL(KIND=DP), INTENT(OUT) :: mu
  ! The lagrange multiplier

                     !>> INTERNAL VARIABLES <<!

#ifdef __USE_PARALLEL
  REAL(KIND=DP), DIMENSION(2) :: mpiLocalSum2, mpiGlobalSum2
#endif

                      !>> INITIALIZATION <<!
                      !>> FUNCTION BODY <<! 
    
#ifdef __USE_PARALLEL
  mpiLocalSum2(1) = SUM(step)
  mpiLocalSum2(2) = SUM(sqrtRho)
  CALL MPI_ALLREDUCE(mpiLocalSum2, mpiGlobalSum2, 2, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpierr /= MPI_SUCCESS) THEN
    CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
    STOP
  ENDIF

  mu = mpiGlobalSum2(1)/mpiGlobalSum2(2)
#else
  mu = SUM(step)/SUM(sqrtRho)
#endif

  Constrain = step - sqrtRho * mu

  RETURN

END FUNCTION Constrain


FUNCTION NewtonDirection(b, b2, lambda, iter) RESULT(x)
!----------------------------------------------------------------------------
! DESCRIPTION:
!   We want to solve the linear equation Ax = b where A is the Hessian of 
!   the energy with respect to rhoR, and b is the negative of the 
!   first gradient of the energy w.r.t. rhoR.
!
!   Note that this is equivalent to minimizing the equation:
!     f(x) = 1/2 xTAx - bTx
!     The gradient to this equation is Ax - b.  The "steepest descent"
!     direction would be negative this, or b - Ax (which is also the residual
!     to the linear equation)
!
! CONDITIONS AND ASSUMPTIONS:
!  For dirichlet boundary conditions, b must be 0 along the boundaries.
!  This function works sort of fortuitously in dirichlet boundary conditions.
!  Technically we should enforce for every array that the boundary is 0.
!  However, this seems like a lot of housekeeping.  It works out that if this
!  is true of b, the code preserves this to be true for the other arrays
!  as well.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------

  USE Hartree, ONLY : SqrtPreconditioner

  IMPLICIT NONE



                         !>> EXTERNAL VARIABLES <<!
    

  REAL(KIND=DP), DIMENSION(LBOUND(rhoR,1):, LBOUND(rhoR,2):, &
                           LBOUND(rhoR,3):, :), INTENT(IN) :: &
      b, &             ! The constrained gradient dE/dx - 2*lamdba*rhoR 
                       !   with 0's
                       !   on the boundary.
      b2

  REAL(KIND=DP), INTENT(IN) :: lambda
  ! The lagrange multiplier
  !
  INTEGER, INTENT(OUT) :: iter

  REAL(KIND=DP), DIMENSION(LBOUND(b,1):UBOUND(b,1),LBOUND(b,2):UBOUND(b,2), &
                           LBOUND(b,3):UBOUND(b,3),SIZE(b,4)) :: x
  ! The newton direction
    
                         !>> INTERNAL VARIABLES <<!
    
  REAL(KIND=DP), DIMENSION(LBOUND(b,1):UBOUND(b,1),LBOUND(b,2):UBOUND(b,2), &
                           LBOUND(b,3):UBOUND(b,3),SIZE(b,4)) :: &
      p, &             ! The next conjugate gradient direction
      Ap, &            ! The result of applying A (the Hessian) to p
      y                ! The solution of My = r, where M is the preconditioner

  REAL(KIND=DP), DIMENSION(LBOUND(b,1):UBOUND(b,1),LBOUND(b,2):UBOUND(b,2), &
                           LBOUND(b,3):UBOUND(b,3),SIZE(b,4)), TARGET :: r
  ! residual, r = b - Ax, also the gradient direction

  REAL(KIND=DP), DIMENSION(:,:,:), POINTER :: r3d

  REAL(KIND=DP) :: &
#ifdef __USE_PARALLEL
      mpiLocalSum, &
#endif
      criteria, &
      alpha, &
      pAp, &
      ryLast, &
      ry, &
      rr, &
      rrLast, &
      epsilon        ! The size of the finite difference step


                           !>> INITIALIZATION <<!
    r3d=>r(:,:,:,1)

    ry = 0

                      !>> FUNCTION BODY <<! 

    ! Start out with initial guess of 0 for the solution
    x = 0._DP

    ! If this is the case, then r = b - Ax
    r = -b!JPotentialCLone(-b)

#ifdef __USE_PARALLEL
    mpiLocalSum = SUM(r**2, interior)
    CALL MPI_ALLREDUCE(mpiLocalSum, rr, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) THEN
      CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
      STOP
    ENDIF
#else
    rr = SUM(r**2, interior)
#endif
    rrLast = rr*2._DP
    criteria = rr
    iter = 1
    DO

    ! Again, the preconditioner: don't work yet.
    ! Solve the preconditioning equation. We got this preconditoner because
    ! Nocedal and Wright said that one possible preconditioner could be 
    ! a simplified form of the main equation. The Laplacian is probably
    ! the most elliptical term, so we will use that.

    ! CAUTION: For optimization purposes, the variable y is used later
    ! as a temporary variable since the value is not "cumulative" ... that is,
    ! later y's do not depend on the value of previous y's. If this is 
    ! not true, we will have to compensate.
      IF(calcReal) THEN
      ELSE
        ! Here is the preconditioner
        IF(usePreconditioner) THEN
          y = SqrtPreconditioner(r3d)
        ELSE
          y = r
        END IF
      END IF
      
      ryLast = ry
#ifdef __USE_PARALLEL
      mpiLocalSum = SUM(r*y, interior)
      CALL MPI_ALLREDUCE(mpiLocalSum, ry, 1, MPI_REAL8, MPI_SUM, &
                         MPI_COMM_WORLD, mpiErr)
      IF (mpiErr/=MPI_SUCCESS) THEN
        CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
        STOP
      ENDIF
#else
      ry = SUM(r*y, interior)
#endif      
      IF(iter == 1) THEN
        ! Calculate the next direction, p.
        p = y
      ELSE
        ! Conjugate gradient part
        p = y + ry/ryLast * p 
      END IF

      ! Find Ap, the result of applying the hessian A to vector p
      ! Do this by a 1st order finite differences approximation.

      ! Picking a good epsilon is sort of a tricky excercise.
      ! A straight epsilon = 1.E-9 seems to work pretty well, but here we have
      ! tried to tailor it to be dependent on the ratio of rhoR and p.
      ! We find the minimum over all processors, so use mpi_allreduce
      ! mpiLocalsum is actually the local minimum, not the sum, but I didn't
      ! feel like creating another variable
#ifdef __USE_PARALLEL
      mpiLocalSum = 1.E-6_DP * MINVAL(ABS(rhoR/p), interior)
      CALL MPI_ALLREDUCE(mpiLocalSum, epsilon, 1, MPI_REAL8, MPI_MIN, &
                         MPI_COMM_WORLD, mpiErr)
      IF (mpiErr/=MPI_SUCCESS) STOP &
        "**PROBLEMS WITH FINDING GLOBAL MIN IN NEWTONDIRECTION**"
#else
      epsilon = 1.E-6_DP * MINVAL(ABS(rhoR/p), interior)
#endif
      IF(calcReal) THEN
      ELSE
!          CALL CalculatePotentialSqrt((rhoR + epsilon * p)**2, Ap, cheat)
          CALL CalculatePotentialPlus((rhoR + epsilon * p)**2, .TRUE., Ap)

      END IF

      ! Here we attempt to calculate the effect of the operator A on the vector
      ! p.  Here, A = Hessian - 2*lambda, since the gradient we have passed in
      ! is modified by the lagrange multiplier.
      ! We need to adjust b a bit because b is really dE/dx - 2 * lambda * x
      Ap = (Ap - b2)/epsilon - lambda*p

      IF(calcReal) THEN
      END IF

#ifdef __USE_PARALLEL
      mpiLocalSum = SUM(p*Ap, interior)
      CALL MPI_ALLREDUCE(mpiLocalSum, pAp, 1, MPI_REAL8, MPI_SUM, &
                         MPI_COMM_WORLD, mpiErr)
      IF (mpiErr/=MPI_SUCCESS) THEN
        CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
        STOP
      ENDIF
#else
      pAp = SUM(p*Ap, interior)
#endif
      ! Check for positive-definiteness here
      IF(pAp < 0._DP) THEN 
        IF(iter == 1) THEN
          IF(calcReal) THEN
          ELSE
            x=-(x-ry/pAp*p)
          END IF
        END IF
        EXIT
      END IF

      ! This is the exact distance we need to move in direction p to get to
      ! the minimum of energy (analytical "line-search")
      alpha = ry / pAp
      x = x + alpha * p
      ! Really, r = b - Ax (or Axnew)  But, we don't have Ax.  
      ! However, we do have Ap, and if we apply that to the equation 
      ! xnew = xold + alpha*p, we get Axnew = Axold + alpha*Ap.
      ! Thus, r = b - Ax_old - alpha*Ap.  But rold = b - Ax_old.  So,
      ! rold = rnew - alpha*Ap
      r = r - alpha * Ap
      
      rrLast = rr
#ifdef __USE_PARALLEL
      mpiLocalSum = SUM(r**2, interior)
      CALL MPI_ALLREDUCE(mpiLocalSum, rr, 1, MPI_REAL8, MPI_SUM, &
                         MPI_COMM_WORLD, mpiErr)
      IF (mpiErr/=MPI_SUCCESS) THEN
        CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
        STOP
      ENDIF
#else
      rr = SUM(r**2, interior)
#endif
      ! Exit criteria. The residual has gone down by 10%, or we have past
      ! 50 iterations, or the residual has stopped changing.
      IF(rr < .1_DP*criteria .OR. iter > 50 &
        .OR. (ABS(rr-rrLast)/rr < .01_DP .AND. iter>9)) EXIT
      iter = iter + 1

    END DO

END FUNCTION NewtonDirection


FUNCTION Step(sqrtRho, time, dir, mask)
!----------------------------------------------------------------------------
! DESCRIPTION:
!   This function takes a step in direction dir from density rho given the 
!   timestep time.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Move to Calculator in the future?
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------

    IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!
    REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: &
      sqrtRho, &
      dir
 
    LOGICAL, DIMENSION(:,:,:,:), INTENT(IN) :: &
      mask

    REAL(KIND=DP), INTENT(IN) :: &
      time

    REAL(KIND=DP), DIMENSION(SIZE(sqrtRho,1), SIZE(sqrtRho,2), &
                             SIZE(sqrtRho,3), SIZE(sqrtRho,4)) :: &
      Step

                         !>> INTERNAL VARIABLES <<!
    REAL(KIND=DP) :: &
      normFactor

    INTEGER :: &
      i, j, k

                         !>> INITIALIZAITION <<!

    Step = 0._DP
    normFactor = 0._DP

                         !>> FUNCTION BODY <<!

    ! Doing it this way with the loop is faster (via an speed test w/ full
    ! optimization) than the easier-to-read way.
    DO k = 1, SIZE(Step,3)
      DO j = 1, SIZE(Step,2)
        DO i = 1, SIZE(Step,1)
          IF(mask(i,j,k,1)) THEN

            Step(i,j,k,1) = sqrtRho(i,j,k,1) + time * dir(i,j,k,1)

            ! One thing you can do is to divide rho by 10 if after the step is
            ! taken we have a negative value. This is to discourage "bouncing"
            ! of density that has already been set to 0. But this seems to
            ! sometimes accelerate deterioration of the WGC. Feel free to try it
            ! IF(Step(i,j,k,1) < 0._DP) Step(i,j,k,1) = Step(i,j,k,1)/10._DP

            Step(i,j,k,1) = Step(i,j,k,1)**2
            normFactor = normFactor + Step(i,j,k,1)

          ELSE
            Step(i,j,k,1) = sqrtRho(i,j,k,1)**2
          END IF
        END DO
      END DO
    END DO

#ifdef __USE_PARALLEL
    mpiLocalSum = normFactor
    CALL MPI_ALLREDUCE(mpiLocalSum, normFactor, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) THEN
      CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***")
      STOP
    ENDIF
#endif

    normFactor = sumRho/normFactor
    WHERE(interior) Step = Step * normFactor

END FUNCTION Step

END SUBROUTINE SqrtNewtonMinimization

END MODULE RhoOptSTN

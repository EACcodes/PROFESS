MODULE RhoOptN
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoOptN
!     |_SUBROUTINE ProjNormRhoMinimization (mainly used)
!       |_SUBROUTINE DensityOptimization
!       |_SUBROUTINE getNextDirection
!       |_SUBROUTINE lineSearch
!       |_SUBROUTINE StepInformation 
!       |_Subroutine SpinThetaSearch
!       |_FUNCTION thetaCGDirection
!     |_SUBROUTINE CheckExit
!     |_SUBROUTINE DenOptReport
!
! DESCRIPTION:
!   This module contains electron minimization method which conserving
!   the electron number during linear search.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
! 1.  Press, W.H., et. al.  Numerical Recipies: The Art of Scientific
!     Computing.  Cambridge University Press, Cambridge, 1986.
! 2.  Gill, P.E., et.al.  Practical Optimization.  Academic Press, San Diego, 
!     1986 (12th printing in 2000).
!------------------------------------------------------------------------------
                              !<< GLOBAL >>
  USE RhoOptimizers
  USE OutputFiles, ONLY: outputUnit
  USE Output, ONLY: QUIT
  USE Output, ONLY: WrtOut
  USE cellInfo, ONLY: kinetic ! Steven's KEDF
  USE CellInfo, ONLY: numSpin
  USE CONSTANTS, ONLY: bohr
  USE RhoLineSearch, ONLY: gradOftheta
  USE RhoLineSearch, ONLY: getNextGradient
  USE RhoLineSearch, ONLY: checkGradient
  USE RhoLineSearch, ONLY: LineSearchTN
  USE RhoDirCG, ONLY: InnerProd

  USE Timer, ONLY: TimerStart
  USE Timer, ONLY: TimerStop
  USE Timer, ONLY: stopwatch 

  USE Report, ONLY: MinimizerReportFooter ! Subroutine that prints out the results.
  USE Report, ONLY: MinimizerReportHeader ! Prints out the header for the minimization report
  USE Report, ONLY: MinimizerReportSteps  ! Subroutine that prints individual steps

  IMPLICIT NONE

  INTEGER :: infoUnit=6

  ! dimension of the FFT box in this subroutine.
  INTEGER :: dimX, dimY, dimZ, dimXYZ


CONTAINS


SUBROUTINE ProjNormRhoMinimization(energy, dir_type)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is a subroutine which used either Truncated Newton (TN) or 
!   Conjugate gradient (CG) method or BFGS method. 
! 
!   Details:
!
!     (1) After calculated the Newton/CG/BFGS directions, 
!         this direction will be projected to make it perpendicular
!         to the current sqrt(rho) vector, then be normalized to the total 
!         electron numbers in the system. 
!
!         Then we optimize the mixing factor: theta,
!         and the new density is formulated as: 
!            sqrt(rho)_new = sqrt(rho)_old*cos(theta) + 
!                      projected_and_normalized_Newton_dir*sin(theta)
!
!         Apparently, the norm of sqrt(rho)_new is automatically conserved 
!         with any theta. And theta should only be allowed to vary from [0,pi/2]
!
!     (2) In the Newton equations: 
!         Hessian*p = g is solved by 1st order finite difference, still 
!         as in Greg Ho's code. But the Hessian*p now is calculated with a
!         corrected formula.
!   
!     (3) Chemical potential (i.e. the Lagrangian multiplier) is calculated
!         with the correct formula (the one in Greg's code is not correct, but
!         works)
!            mu = SUM(sqrt(rho)*sqrt_pot*dV)/2/toteleNum
!
!     (4) For spin-polarized case, there will be two theta values. Using two line searches
!         seperately is not the most efficient way. Instead, a small CG solver is used
!         to obtain the optimal theta vector (theta1,theta2) to minimize the energy.
!
!   *********************************************************************
!   **This subroutine only works for periodic boundry condition now.*****
!   *********************************************************************
!
! REFERENCES:
!   Jiang Hong and Weitao Yang, J. Chem. Phys. 2004 (chemical potential formula)
!   M.C. Payne, M.P. Teter, D.C. Allenet al. Rev. Mod. Phys. 1992
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/02/2004  Subroutine adapted from direct minimization.  (GSH)
!   06/01/2008  (1) Projects the truncated Newton dir to conserve the total
!                   electron number, therefore we get rid of rescaling the
!                   electron number every time.  "Rescaling" is not a proper
!                   procedure, the electron number should be conserved in the
!                   way implemented in this subroutine. 
!               (2) Removed all the cheat option
!               (3) Removed all code related to the dirichlet boundry
!               (4) NewtonDirection subroutine: correctly compute Ap now
!                                                               (Chen Huang)
!   08/29/08    the finite difference in NewtonDir routine is corrected, 
!               previous one is not efficient once having vacuum. (Chen Huang)
!   09/08/08    Add CG direction, now can do CG (Chen Huang)
!   10/07/11    Add spin-polarized case (Steven)
!   2014        Reorganize the code structure (Mohan Chen)
!------------------------------------------------------------------------------

  USE MathFunctions, ONLY: MinMaxVal  ! this one takes parallel into consideration
  USE Output, ONLY: QUIT
  USE CellInfo, ONLY: cell
  USE Sys, ONLY: rhoR, rho0 ! electronic density in real space
  USE KEDF_EvW, ONLY: kloc, aloc, tolk
  USE KEDF_WGCD, ONLY: InitializeWGCD   ! Initialize Steven's WGCD method
  USE KEDF_WGCD, ONLY: ComputeFrMatrix  ! Compute the Fr function in Steven's WGCD method
  USE RhoDirBFGS, ONLY: InitializeBFGS
  USE RhoDirBFGS, ONLY: CleanBFGS
  USE RhoDirBFGS, ONLY: UpdateBFGS     
  USE RhoDirBFGS, ONLY: ChemicalPotential
  USE KEDF_DenDec, ONLY: do_den_dec
  USE KEDF_DenDec, ONLY: potDD 

  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(OUT) :: energy(:) 
  ! where the final energies are stored   
  !
  INTEGER, INTENT(INOUT) :: dir_type        
  ! type of direction
  ! 1: we do truncated newton
  ! 2: we do conjugate gradient
  ! 3: we do lbfgs

                    !>> INTERNAL VARIABLES <<!
                    
  REAL(KIND=DP), DIMENSION(0:SIZE(rhoR,1)-1, 0:SIZE(rhoR,2)-1, &
                           0:SIZE(rhoR,3)-1, SIZE(rhoR,4)) :: &
    phi,  &       ! phi^2 = rho, can be negative
    tempPhi, &    ! temporary storage of phi
    dirNext, &    ! The next direction that we would like to search in,
    dirNext_old, &! previous dirNext
    dEdPhi, &     ! The d(Energy)/d(Phi)
    dLdPhi ,&     ! The d(Lagrangian)/d(Phi)
    dLdPhi_old    ! previous dLdPhi

  REAL(KIND=DP) :: oldr = 0.0_DP ! for SC-KEDF
  REAL(KIND=DP) :: olda = 0.0_DP ! for SC-KEDF

  REAL(KIND=DP) :: &
    totalpotentialNorm, &
    potentialNorm(Size(rhoR,4)), &
    numPoints,   &
    mu0(size(rhoR,4)), &            ! chemical potenital
    totEleNum(SIZE(rhoR,4)),     & ! total number of electrons in the system
    dV, &            ! volume element on the real grid point
    e(3)             ! to store the last 3 successive total energy
    
  INTEGER :: &
    doingSD, &    
    newtonIter, & ! Number of CGMIN loops in the newton step
    success, &    ! whether our line minimization succeeded
                  !  =1 if it did succed
                  !  =-1 if it exited because of error in line minimizer
                  !  =-2 if it exited due to a warning
    NTdirFlag

  LOGICAL :: restart  ! whether the present step restarted from the steepest

  ! Variables for the line search algorithm:
  CHARACTER(LEN=60) :: date, time

  REAL(KIND=DP) :: potNorm
  REAL(KIND=DP) :: theta(SIZE(rhoR,4))        ! the mixing factor
  REAL(KIND=DP) :: gradth(SIZE(rhoR,4))       ! dL/d(theta), L = E - mu (SUM(rho)-toteleNum), the Lagrangian
  REAL(KIND=DP) :: ftol = 1.E-4_DP            ! new energy should be smaller than the older one.
  REAL(KIND=DP) :: gtol = 0.2_DP              ! 0.1 is for CG which needs a more exact line-search,
                                              ! here we use 0.5. gtol=0.9 is too inaccurate
                                              ! gtol MUST > ftol to have wolfe condition always satisified
  REAL(KIND=DP),parameter :: xtol = 1.E-12_DP
    
  INTEGER :: isp

  INTEGER :: magIter
  ! magnetization iteration number

  REAL(KIND=DP) :: Ntot
  REAL(KIND=DP) :: magnetization
  REAL(KIND=DP) :: ee1 = 0.0_DP
  REAL(KIND=DP) :: ee2 = 0.0_DP
  REAL(KIND=DP) :: m1, m2, de, ee3, m3, de1, b

  TYPE(stopwatch):: watch          ! Timer
  LOGICAL :: switchStp, switchStp2 !
  LOGICAL :: exit_flag = .FALSE.   ! (Shin, mohan add 10-01-12)
  LOGICAL :: fixedN

  ! check if new iteration will be applied ( extra density iteration ) 
  LOGICAL :: new_iteration = .FALSE.
  LOGICAL :: useWGCD = .FALSE.
  INTEGER :: iter = 0
  ! density optimization iteration counter
  !
  INTEGER :: imag ! iteration flag for magnetic loop
  INTEGER :: iter_extra = 0 ! index of iteratoin for extra density optimization.
  LOGICAL :: grad_exit_flag ! test is the gradient is ok. 

  INTEGER :: total ! mohan testing

  EXTERNAL DCSRCH

  REAL(KIND=DP) :: sc_max_density = 0.0_DP


                      !>> INITIALIZATION <<!

  CALL Title("RhoOptimizers::ProjNormRhoMinimization")

  dimX  = SIZE(rhoR,1)
  dimY  = SIZE(rhoR,2)
  dimZ  = SIZE(rhoR,3)
  dimXYZ = dimX * dimY * dimZ

  WRITE(outputUnit, '(/A)') " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           DENSITY OPTIMIZATION"
 
  watch = TimerStart()

  ! autosetting 1
  IF(numSpin == 1)THEN
     fixedN = .TRUE.
  ELSE
     IF(fixedmag)THEN
        fixedN = .TRUE.
     ELSE
        fixedN = .FALSE.
     ENDIF
  ENDIF

  ! Set Wolfe Conditions for convergence
  ! criteria of line search
  SELECT CASE(dir_type)
    CASE(1) !NTN
      ftol = 1e-4
      gtol = 2e-1
    CASE(2) !NCG
      ftol = 1e-4
      gtol = 1e-2
    CASE(3) !NBF
      ftol = 1e-4
      gtol = 2e-1
      CALL InitializeBFGS(rhoR)
  END SELECT
 

  ! initialize 2
  dirNext     = 0.D0
  dirNext_old = 0.D0
  dLdPhi      = 0.D0
  dLdPhi_old  = 0.D0
  dEdPhi      = 0.D0

  ! initialize 3
  rhoOutcome = -1
  doingSD = -1       ! we will do sqrt truncated Newton first
  restart = .FALSE.
  
  ! STEVEN's NEW KEDF !
  IF(kinetic == 12) THEN
    CALL InitializeWGCD( dimX, dimY, dimZ, kinetic, pot_tol )
    useWGCD = .TRUE.
  ENDIF

  ! calculate the number of points.
  numPoints =  REAL(dimX, KIND=DP) &
              *REAL(dimY, KIND=DP) &
              *REAL(dimZ, KIND=DP) ! mohan fix bug 2013-07-23

  CALL ReduceRealLevel1(numPoints)
  !WRITE(outputUnit,*) "numPoints after reduce = ",numPoints


  ! compute the total electron in cell
  dV = cell%vol/numPoints
  DO isp = 1, numSpin
    totEleNum(isp) = SUM(rhoR(:,:,:,isp))*cell%dV
    CALL ReduceRealLevel1(totEleNum(isp))
    WRITE(outputUnit,*) "Electron number for spin ", isp, " is ",totEleNum(isp)
  ENDDO
  
  ! total electron number
  Ntot = SUM(totEleNum)

  WRITE(outputUnit,'(A,F20.6,A)') " Cell Volume      : " , cell%vol*bohr**3 , " Angstrom^3"
  WRITE(outputUnit,'(A,I20)')     " Ion number       : " , cell%numIon
  WRITE(outputUnit,'(A,F20.6)')   " Electron number  : " , Ntot 
  WRITE(outputUnit,'(A,F20.6,A)') " Ion density      : " , cell%numIon/cell%vol/bohr**3 , " Ion/Angstrom^3"
  WRITE(outputUnit,'(A,F20.6,A)') " Electron density : " , Ntot/cell%vol/bohr**3 , " Electron/Angstrom^3"
  
  ! magnetization
  IF (numSpin == 2) THEN
    magnetization = totEleNum(1)-totEleNum(2)
  ENDIF
  
  de = 1.d0
  magIter = 0


  iter_extra = 0

! WGCD loop for Fr matrix self-consistent
! or loop for EvW KEDF
10 CONTINUE

   iter_extra = iter_extra + 1
!  WRITE(outputUnit,'(A,I5,A,I5)') " Extra Eleconic Iteration Index", iter_extra, "/", niter_extra

  IF( numSpin == 2) THEN
    WRITE(message, '(A)')'Miter       totEnergy(Ha)      dE(Ha)      TotMag     Ele(UP)   Ele(DOWN) Eiter' 
    CALL WrtOut(infoUnit, message)
  ENDIF

  DO imag=1,nmag ! Loop for magnetic moment

    IF(imag>1) THEN
      WRITE(outputUnit,*) "New loop for magnetic moment, fixedN=",fixedN
    ENDIF

    ! (1.1) Initialize
    rhoOutcome = -1
    success = 1
    doingSD = -1       ! we will do sqrt truncated Newton first
    restart = .FALSE.
    theta(:) = 2.d-1     ! the starting value for mixing factor

    CALL MinimizerReportHeader

    ! (1.2) initialize our phi
    DO isp = 1, numSpin
        phi(:,:,:,isp) = SQRT(ABS(rhoR(:,:,:,isp)))
    ENDDO

    ! (1.3) calculate the potential
    CALL CalculatePotentialPlus(phi(:,:,:,:)**2,.TRUE.,dEdPhi(:,:,:,:),energy) ! Update for Steven's KEDF
    dEdPhi = dEdPhi * SIGN(1._DP,phi)    ! Everything still positive at this point
    

    ! (1.4) calculate the chemical potential.
    CALL ChemicalPotential(dEdPhi,phi,totEleNum,mu0)

    DO isp = 1, numSpin
       dLdPhi(:,:,:,isp) = dEdPhi(:,:,:,isp) - mu0(isp)*2._DP*phi(:,:,:,isp)
    ENDDO
  
    numEnergyTotal = 1
    numPotentialTotal = 1

    ! (1.5) Calculate the norm of the potential
    ! if potentialNorm very close to zero, then we are almost
    ! done the density optimization
    DO isp=1, numSpin
      potentialNorm(isp) = SUM(dLdPhi(:,:,:,isp)**2)
      CALL ReduceRealLevel1(potentialNorm(isp))
    ENDDO
 
    IF (numSpin == 1) THEN
       totalpotentialNorm = potentialNorm(1)
       totalpotentialNorm = SQRT(totalpotentialNorm/numPoints)
    ELSE
       totalpotentialNorm = SUM(potentialNorm)
       totalpotentialNorm = SQRT(totalpotentialNorm/numPoints/2.d0)
    ENDIF
    CALL MinimizerReportSteps(0, energy, totalpotentialNorm, 0, 0, .FALSE., 1, 0._DP)
    
    ! (1.6) backup histroy energies, used for stop criteria
    e(1) = energy(1)
    e(2) = 0._DP
    e(3) = 0._DP
 
    CALL DenOptReport(dir_type)

    ! (1.7) Optimize the density.
    CALL DensityOptimization()

    ! For Isaac's KEDF
    ! only if the density is converged, then we consider to do
    ! extra iteration.
    IF (exit_flag .AND. kloc>0._DP) THEN ! (Shin)

      sc_max_density=MinMaxVal(tempphi(:,:,:,1)**2, 'max')
      aloc = kloc*((sc_max_density-rho0)/rho0)**2

      WRITE(message,*) "SC KEDF: Error between new and old density = ", ABS(oldr-sc_max_density)
      CALL WrtOut(6, message)

      WRITE(message,*) "SC KEDF: [ (rho_max-rho0)/rho0 ]^2 = ", ((sc_max_density-rho0)/rho0)**2
      CALL WrtOut(6, message)

      WRITE(message,*) "SC KEDF: Error between new and old A = ", ABS(olda-aloc)
      CALL WrtOut(6, message)

      IF (ABS(olda-aloc)>tolk) THEN 
        WRITE(outputUnit,*) " A Parameter in EvW KEDF is ", aloc
        WRITE(outputUnit,*) " K Parameter in EvW KEDF is ", kloc
        WRITE(message,*) "new A=", aloc, " K=", kloc
        CALL WrtOut(6, message)

        oldr = sc_max_density 
        olda = aloc

        IF( iter_extra < niter_extra ) THEN
          exit_flag = .FALSE.
          goto 10
        ENDIF
      END IF

    ENDIF

    !If fixed-Mag, then exit. if not, do tuned mag
    IF(fixedN .EQV. .TRUE.) THEN
      magIter = magIter + 1

      IF(numSpin==2) THEN
        WRITE(message, '(I5,ES20.12,ES12.4,ES12.4,ES12.4,ES12.4,I5)') magIter, e(1), e(1)-ee1, &
        ABS(magnetization), totEleNum(1), totEleNum(2), iter
        CALL WrtOut(6, message)
      ENDIF

      EXIT
    ELSE
      !WRITE(message,51) totEleNum(1), totEleNum(2), totEleNum(1)-totEleNum(2)
      magIter = magIter + 1
      WRITE(message, '(I5,ES20.12,ES12.4,ES12.4,ES12.4,ES12.4,I5)') magIter, e(1), e(1)-ee1, &
      ABS(magnetization), totEleNum(1), totEleNum(2), iter
     ! CALL WrtOut(6, message)
 
      ! difference of energy is less than a small value,
      ! finish the iteration of magnetization,
      ! otherwise, go on.
      IF(ABS(de) .LT. tolm) GOTO 30
      ee3 = ee2
      ee2 = ee1
      ee1 = e(1)
      e(2:3) = 0.d0
      m3 = m2
      m2 = m1
      m1 = magnetization
      de = (ee1-ee2)/(m1-m2)

      IF(magIter == 1) magnetization = m1 + 0.1
      IF(magIter == 2) magnetization = m1 - 0.2
      IF(magIter == 3) magnetization = m1 - de/ee1
      IF(magIter .GT. 3) THEN
        b=(m1-m2)/(m3-m2)
        de1=2.d0*(ee1-(1-b)*ee2 - b*ee3)/((1.d0-1.d0/b)*(m1-m2)**2)
        ! maximal iteration number: 60
        IF(magIter .GT. 60) GOTO 30
        magnetization = m1 - de/de1
      ENDIF

      rhoR(:,:,:,1) = rhoR(:,:,:,1)/totEleNum(1)
      rhoR(:,:,:,2) = rhoR(:,:,:,2)/totEleNum(2)
      totEleNum(1) = 0.5d0*(Ntot + magnetization)
      totEleNum(2) = 0.5d0*(Ntot - magnetization)
      rhoR(:,:,:,1) = rhoR(:,:,:,1)*totEleNum(1)
      rhoR(:,:,:,2) = rhoR(:,:,:,2)*totEleNum(2)
    ENDIF

  
  ENDDO ! End spin Loop

30 CONTINUE  

  IF ( rankGlobal==0 ) THEN
    CALL MinimizerReportFooter( & 
       iter, rhoOutcome, energy, numEnergyTotal, &
       numPotentialTotal, TimerStop(watch))
  ENDIF

  ! Make sure to end up with the actual rho, not the square root of rho
  ! get the new charge density.
  rhoR = phi**2

  ! calculate the number of points.
  total = SIZE(rhoR,1) * SIZE(rhoR,2) * SIZE(rhoR,3)

  IF(useWGCD .EQV. .TRUE.) THEN
    ! this is WGCD formalism
    CALL ComputeFrMatrix(rhoR, dimX, dimY, dimZ, numSpin, new_iteration, kinetic, pot_tol)
    IF(new_iteration .EQV. .TRUE.) THEN
      IF( iter_extra .LE. niter_extra ) THEN
        GOTO 10
      ENDIF
    ENDIF
  ENDIF


  ! check the total electron number
  DO isp = 1, numSpin
     totEleNum(isp) =  SUM(rhoR(:,:,:,isp))*dV
     CALL ReduceRealLevel1(totEleNum(isp))
!     WRITE(*,*) "mohan test, totEleNum=", totEleNum(isp)
  ENDDO
!  Ntot = Sum(totEleNum)
  
  ! write out the message
  WRITE(message,'(A,Es19.12,A,Es11.4,A,F10.1)') & 
    & ' Final total energy: ', energy(1), &
    & ' Ha, volume=', cell%vol, ' bohr^3, totQ=', Ntot 
  CALL WrtOut(infoUnit,message)


  ! Print the header
  CALL DATE_AND_TIME(date, time)  
  WRITE(message,'(13a)') " Current Time : " // date(5:6) // "/" // &
    date(7:8) // "/" // date(1:4) // " at " // & 
    time(1:2)// ":" // time(3:4) // ":" // time(5:6) // " "
  CALL WrtOut(infoUnit,message);
  CALL WrtOut(infoUnit,'')


  ! if we density decomposition, we need to save the 
  ! potential here
  ! note, we need dE/dRho, not dE/dPhi
  ! only do spin=1 right now
  IF(do_den_dec==1) THEN
!    potDD(:,:,:) = dEdPhi(:,:,:,1) / (2.0d0*SQRT(rhoR(:,:,:,1)))
    potDD(:,:,:) = dLdPhi(:,:,:,1) / (2.0d0*SQRT(rhoR(:,:,:,1)))
  ENDIF


  RETURN

CONTAINS


SUBROUTINE DensityOptimization()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Get the next direction through different method.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
 
  IMPLICIT NONE
  
  !! >> INITIALIZE << !!
  CALL Title("RhoOptN::DensityOptimization")

  iter = 1  ! the counter for electronic iteration number

  !! >> FUNCTION << !!
  DO ! Level 2: electronic iteration
    switchStp  = .False.
    switchStp2 = .False.

    CALL getNextDirection()

    DO 
      CALL getNextGradient(totEleNum, phi, dEdPhi, theta, gradth, dirNext)
      CALL checkGradient(switchStp, switchStp2, gradth, grad_exit_flag, rhoOutcome, dirNext, dLdPhi)
      IF(grad_exit_flag) EXIT
    ENDDO 

    ! save the direction
    dirNext_old = dirNext

    IF (rhoOutcome==1) EXIT

    ! LINE SEARCH based on theta
    numEnergyLineMin = 0 
    CALL LineSearchTN(iter, theta, energy, gradth, ftol, gtol, xtol, success, &
    rhoOutcome, tempPhi, phi, dirNext, dEdPhi)

    ! Output step information --------------------------------
    CALL StepInformation(iter,mu0,dEdPhi(:,:,:,:),potNorm,tempphi(:,:,:,:),theta,NTdirFlag, & 
                         totEleNum,energy,numEnergyLineMin, restart,success,newtonIter,dir_type)

    
    !----------- CHECK THE CONVERGENCE ---------------
    CALL CheckExit(maxIter,iter,conv_check,e,energy(1),potNorm,rhoOutcome,exit_flag)    

    !update information for BFGS method.
    IF(dir_type==3) THEN
      IF(numSpin==2) THEN
        DO isp =1, Numspin
          tempPhi(:,:,:,isp) = phi(:,:,:,isp)*COS(theta(isp)) + dirNext(:,:,:,isp)*SIN(theta(isp))
        ENDDO
      ENDIF
      CALL UpdateBFGS(phi,tempPhi,dEdPhi,dLdPhi,totEleNum,dV)
    ENDIF

    !update the phi, i.e., the density.
    IF(numSpin == 1) THEN
      phi=tempPhi
    ELSE
      DO isp = 1, numSpin
        Phi(:,:,:,isp) = phi(:,:,:,isp)*COS(theta(isp)) + dirNext(:,:,:,isp)*SIN(theta(isp))
      ENDDO
    ENDIF

    ! update chemical potential and the dL/dPhi with new phi
    ! dEdPhi is already up to date from the line search loop above
    CALL ChemicalPotential(dEdPhi,phi,totEleNum,mu0)

    DO isp=1,numSpin
      dLdPhi(:,:,:,isp)  = dEdPhi(:,:,:,isp) - 2._DP*mu0(isp)*phi(:,:,:,isp)
    ENDDO
        
    !check if converged, if yes, exit the fixed-Mag optimization loop
    IF (exit_flag .EQV. .TRUE.) EXIT

    iter = iter + 1
    ! WRITE(outputUnit,*) " Electronic Iteration Number ", iter

    ! total energy line search number
    numEnergyTotal = numEnergyTotal + numEnergyLineMin
    numPotentialTotal = numPotentialTotal + numEnergyLineMin

  END DO ! End Density Optimization Loop

END SUBROUTINE DensityOptimization


SUBROUTINE getNextDirection()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Get the next direction through different method.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
 
  USE RhoDirCG, ONLY: CGDirection
  USE RhoDirNew, ONLY: NewtonDirection
  USE RhoDirBFGS, ONLY : BFGSDirection

  IMPLICIT NONE 

  CALL Title("RhoOptN::getNextDirection")

  SELECT CASE(dir_type)
    CASE(1) 
    ! compute the Newton direction
    DO isp=1,numSpin
      CALL NewtonDirection(phi,dLdPhi,mu0(isp),newtonIter,NTdirFlag,isp,dimX,dimY,dimZ,numSpin,dirNext(:,:,:,isp))
      ! numPotentialTotal = numPotentialTotal + newtonIter
    ENDDO

    CASE(2)
      ! compute CG direction. if mixed method, CG is first used
      CALL CGDirection(iter,dimX,dimY,dimZ,numSpin, &
                       dLdPhi,dLdPhi_old,dirNext,dirNext_old)

    CASE(3)
      ! compute BFGS direction
      ! dirNext(:,:,:,:) = BFGSDirection(energy(1),phi,dLdPhi,dLdPhi_old,dV)
      CALL BFGSDirection(energy(1),phi,dLdPhi,dLdPhi_old,dV,dirNext)

  END SELECT

  RETURN

END SUBROUTINE getNextDirection
 

  
SUBROUTINE StepInformation(i,mu0,dEdPhi,potentialNorm,phi,theta,NTdirFlag, & 
                           totEleNum,energy,numEnergyLineMin, & 
                           restart,success,newtonIter,do_CG)
!----------------------------------------------------------------------------
! DESCRIPTION:
! Simply output the step information for the
! density optimization .
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014-07-24 Mohan Chen
!------------------------------------------------------------------------------

  USE KEDF_HC10, ONLY: nref

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: do_CG
  REAL(KIND=DP), INTENT(OUT) :: potentialNorm
  REAL(KIND=DP), INTENT(IN) :: mu0(:)
  REAL(KIND=DP), INTENT(IN) :: theta(:)
  REAL(KIND=DP), INTENT(IN) :: totEleNum(:)
  REAL(KIND=DP), INTENT(IN) :: energy(9)
  REAL(KIND=DP), INTENT(IN) :: dEdPhi(:,:,:,:)  ! dEnergy/d(phi)
  REAL(KIND=DP), INTENT(IN) :: phi(:,:,:,:)

  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: numEnergyLineMin 
  INTEGER, INTENT(IN) :: success
  INTEGER, INTENT(IN) :: newtonIter
  INTEGER, INTENT(IN) :: NTdirFlag

  LOGICAL :: restart
     
  REAL(KIND=DP), DIMENSION(numSpin) :: totEleNum2
  ! a temporay var to store total electron number
  !
  REAL(KIND=DP), DIMENSION(numSpin) :: minPot, maxPot
  REAL(KIND=DP), DIMENSION(numSpin) :: minDen, maxDen

  INTEGER :: isp

  DO isp = 1, numSpin
    ! For output purpose only
     maxPot(isp) = MinMaxVal(dEdPhi(:,:,:,isp),'max')  ! MinMaxVal function returns the
     minPot(isp) = MinMaxVal(dEdPhi(:,:,:,isp),'min')  ! max or min of its input, 
     maxDen(isp) = MinMaxVal(tempphi(:,:,:,isp)**2, 'max')  ! automatically taking 
     minDen(isp) = MinMaxVal(tempphi(:,:,:,isp)**2, 'min')  ! parallel into consideration
  ENDDO
    
  IF(numSpin == 1) THEN
     potentialNorm = SUM((dEdPhi(:,:,:,1)-2._DP*tempphi(:,:,:,1)*mu0(1))**2)
  ELSE
    potentialNorm = SUM((dEdPhi(:,:,:,1)-2._DP*tempphi(:,:,:,1)*mu0(1))**2)
    potentialNorm = potentialNorm + SUM((dEdPhi(:,:,:,2)-2._DP*tempphi(:,:,:,2)*mu0(2))**2)
  ENDIF

  CALL ReduceRealLevel1(potentialNorm)
    
  IF(numSpin == 1) THEN
      potentialNorm = SQRT(potentialNorm/numPoints)
  ELSE
      potentialNorm = SQRT(potentialNorm/numPoints/2.d0)
  ENDIF

  DO isp =1, numSpin
    CALL MinimizerReportSteps( & 
    i, energy, potentialNorm, newtonIter, &
    numEnergyLineMin,restart,success,theta(isp))
  ENDDO
    
  DO isp = 1,numSpin
    totEleNum2(isp) = InnerProd(tempphi(:,:,:,isp),tempphi(:,:,:,isp))
  ENDDO
    
  IF(.FALSE.) then
      ! silence compiler about unused phi. I will keep the dummy, as it may be useful one day for
      ! debugging (JMD 2013-10-23)
    WRITE(*,*) "DEBUG: Phi dimension ", SIZE(phi,1), SIZE(phi,2), SIZE(phi,3)
  ENDIF

  IF(numSpin == 2) THEN
    infoUnit = outputUnit
  ENDIF

  IF(numSpin == 1) then
    IF (kinetic == 11) THEN 
    ! Interpolating the kernels
      IF (do_CG==1) THEN 
        WRITE(message, & 
              '(I3, Es20.10, Es10.2, Es10.2, I3, Es11.3, Es10.2,A,Es9.2,  Es10.2,A,Es9.2,  a,I4)')&
              i,energy(1),theta(1),potentialNorm,NTdirFlag,(totEleNum2(1)-totEleNum(1)), & 
              minDen(1),'/',maxDen(1),minPot(1),'/',maxPot(1),' ',nref
      ELSE
        WRITE(message, & 
              '(I3, Es20.10, Es10.2, Es10.2, Es11.3, Es10.2,A,Es9.2,  Es10.2,A,Es9.2,  a,I4)')&
              i,energy(1),theta(1),potentialNorm,(totEleNum2(1)-totEleNum(1)), & 
              minDen(1),'/',maxDen(1),minPot(1),'/',maxPot(1),' ',nref
        ENDIF
    ELSE
      IF (do_CG==1) THEN 
         WRITE(message, & 
              '(I3, Es20.10, Es10.2, Es10.2, I3, Es11.3, Es10.2,A,Es9.2,  Es10.2,A,Es9.2)')&
              i,energy(1),theta(1),potentialNorm,NTdirFlag,(totEleNum2(1)-totEleNum(1)), & 
              minDen(1),'/',maxDen(1),minPot(1),'/',maxPot(1)
      ELSE
        WRITE(message, & 
              '(I3, Es20.10, Es10.2, Es10.2, Es11.3, Es10.2,A,Es9.2,  Es10.2,A,Es9.2)')&
              i,energy(1),theta(1),potentialNorm,(totEleNum2(1)-totEleNum(1)), & 
              minDen(1),'/',maxDen(1),minPot(1),'/',maxPot(1)
        ENDIF
      ENDIF
      CALL WrtOut(infoUnit,message)
    ELSE
      Do isp =1, numSpin 
        IF (kinetic == 11) THEN 
          ! Interpolating the kernels
            IF (do_CG==1) THEN 
                 WRITE(message, & 
                 '(I3, Es20.10, Es10.2, Es10.2, I3, Es11.3, Es10.2,A,Es9.2,  Es10.2,A,Es9.2,  a,I4)')&
                 i,energy(1),theta(isp),potentialNorm,NTdirFlag,(totEleNum2(isp)-totEleNum(isp)), & 
                 minDen(isp),'/',maxDen(isp),minPot(isp),'/',maxPot(isp),' ',nref
            ELSE
                 WRITE(message, & 
                 '(I3, Es20.10, Es10.2, Es10.2, Es11.3, Es10.2,A,Es9.2,  Es10.2,A,Es9.2,  a,I4)')&
                 i,energy(1),theta(isp),potentialNorm,(totEleNum2(isp)-totEleNum(isp)), & 
                 minDen(isp),'/',maxDen(isp),minPot(isp),'/',maxPot(isp),' ',nref 

            ENDIF
        ELSE
            IF (do_CG==1) THEN 
                 WRITE(message, & 
                 '(I3, Es20.10, Es10.2, Es10.2, I3, Es11.3, Es10.2,A,Es9.2,  Es10.2,A,Es9.2)')&
                 i,energy(1),theta(isp),potentialNorm,NTdirFlag,(totEleNum2(isp)-totEleNum(isp)), & 
                 minDen(isp),'/',maxDen(isp),minPot(isp),'/',maxPot(isp)
            ELSE
                 WRITE(message, & 
                 '(I3, Es20.10, Es10.2, Es10.2, Es11.3, Es10.2,A,Es9.2,  Es10.2,A,Es9.2)')&
                 i,energy(1),theta(isp),potentialNorm,(totEleNum2(isp)-totEleNum(isp)), & 
                 minDen(isp),'/',maxDen(isp),minPot(isp),'/',maxPot(isp)
            ENDIF
        ENDIF
        CALL WrtOut(infoUnit,message)
      ENDDO
    ENDIF

  IF(numSpin == 2) THEN
    infoUnit = 6
  ENDIF

  RETURN 

END SUBROUTINE StepInformation


END SUBROUTINE ProjNormRhoMinimization


SUBROUTINE CheckExit(maxIter,i,conv_check,e,e_now,potNorm,rhoOutcome,exit_flag)
!------------------------------------------------------------------------------
! DESCRIPTION:
!  This subroutine checks if rho optimization is done.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014-07-24 Mohan Chen
!-------------------------------------------------------------------

  USE CellInfo, ONLY: cell
    
  IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!

  INTEGER, INTENT(IN) :: maxIter
  ! maximum allowed iteration
  !
  INTEGER, INTENT(IN) :: i
  ! current # of iteration
  !
  CHARACTER(len=500), INTENT(IN) :: conv_check   
  ! use pot or energy to check the convergence
  !
  REAL(KIND=DP), INTENT(IN) :: e_now
  ! current total energy
  !
  REAL(KIND=DP), INTENT(IN) :: potNorm       
  ! potential norm
  !
  REAL(KIND=DP), INTENT(INOUT) :: e(3)
  ! store energies
  !
  INTEGER, INTENT(INOUT) :: rhoOutcome 
  ! status for rho optimization
  !
  LOGICAL, INTENT(OUT) :: exit_flag
  ! true if exit
  !

                           !>> INTERNAL VARIABLES <<!

  LOGICAL :: eConv, pConv

                           !>> INITIALIZATION <<!

  exit_flag = .FALSE.
  pConv = .FALSE.
  eConv = .FALSE.

  IF(numSpin == 2) THEN
    infoUnit = outputUnit
  ENDIF

                           !>> FUNCTION BODY <<!

  ! Check potential convergence criterion met
  IF ( potNorm < pot_tol ) THEN
    pConv = .TRUE.
  ENDIF

  ! Check energy convergence criterion met
  e(3) = e(2); 
  e(2) = e(1); 
  e(1) = e_now;
  IF (i>=3 .AND. ABS(e(1)-e(2))<tole  & 
           .AND. ABS(e(1)-e(3))<tole) THEN
     eConv = .TRUE.
  ENDIF

  ! See if desired convergence criteria met


  IF (conv_check(1:3) == "ENE") THEN
    IF (eConv) THEN 
       rhoOutcome = 0
       WRITE(message,'(A,Es12.5,A)')  ' Total energy converged in the last three iteration within:', tole, ' Ha'
       CALL WrtOut(infoUnit,message)
       exit_flag = .TRUE.
    ENDIF
  ENDIF
    
  IF (conv_check(1:3) == "POT") THEN
    IF (pConv) THEN 
      rhoOutcome = 0
      WRITE(message,'(A,Es12.5,A)')  ' Potential norm is converged to ', pot_tol, ', density optimization finished.'
      CALL WrtOut(infoUnit,message)
      exit_flag = .TRUE.
    ENDIF
  ENDIF

  IF (conv_check(1:4) == "BOTH") THEN
    IF (pConv .or. eConv) THEN 
      rhoOutcome = 0
      WRITE(message,'(A,Es12.5,A,E12.5,A)')  ' Electronic iteration converged. potTol=', pot_tol, ' eTol=', tole, ' done!'
      CALL WrtOut(infoUnit,message)
      exit_flag = .TRUE.
    ENDIF
  ENDIF


  ! Check for failure in optimization
  if (rhoOutcome == 1) then 
    WRITE(message,'(A)') ' Density optimization failed.'
    CALL WrtOut(infoUnit,message)
    exit_flag = .TRUE.
  endif

  ! Here we exit if we took more than maxIter iterations.  
  IF (i >= maxIter) THEN 
    rhoOutcome = 1
    exit_flag = .TRUE.
  ENDIF

  IF(numSpin==2) THEN
    infoUnit = 6
  ENDIF

  RETURN

END SUBROUTINE  CheckExit


SUBROUTINE DenOptReport(dir_type)
!------------------------------------------------------------------------------
! DESCRIPTION:
! print out the density optimization information for spin=1 case,
! for spin=2 case, we didn't show this detail.
! 
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014-07-24 Mohan Chen
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell
  USE RhoDirCG, ONLY: cg_alg
 
  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  INTEGER, INTENT(IN) ::dir_type 

  IF(numSpin == 2) THEN
    infoUnit = outputUnit
  ENDIF

  !! >> FUNCTION << !!
  SELECT CASE(dir_type)

    CASE(1) ! NTN 
      WRITE(message, *) " "
      CALL WrtOut(infoUnit,message)
      WRITE(message, *) '(RhoOptimizers) Use NTN (Truncated Newton) minimization method for density.'
      CALL WrtOut(infoUnit,message) 
      WRITE(message, '(A)')' Iter    totEnergy(Ha)     Theta   & 
            &potNorm  New    ele#       min/max(den)      min/max(dE/dPhi)'
      CALL WrtOut(infoUnit, message)

    CASE(2) ! NCG
      WRITE(message, *) " "
      CALL WrtOut(infoUnit,message)
      WRITE(message, *) "(RhoOptimizers) Use NCG (Conjugate Gradient)", cg_alg(1:2), &
      " minimization method for density."
      CALL WrtOut(infoUnit,message) 
      WRITE(message, '(A)')' Iter    totEnergy(Ha)     Theta   & 
            &potNorm       ele_ch       min/max(den)      min/max(dE/dPhi)'
      CALL WrtOut(infoUnit, message)

    CASE(3) ! NBF 
      WRITE(message, *) " "
      CALL WrtOut(infoUnit,message)
      WRITE(message, *) '(RhoOptimizers) Use NBF (BFGS) minimization method for density.'
      CALL WrtOut(infoUnit,message) 
      WRITE(message, '(A)')' Iter    totEnergy(Ha)     Theta   & 
         &potNorm       ele_ch       min/max(den)      min/max(dE/dPhi)'
      CALL WrtOut(infoUnit, message)

  END SELECT

  !! >> Finish << !!

  IF(numSpin == 2) THEN
    infoUnit = 6 
  ENDIF

  RETURN

END SUBROUTINE DenOptReport


END MODULE RhoOptN

MODULE IonOptBFGS
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE IonOptBFGS
!     |_SUBROUTINE IonLBFGS 
!
! DESCRIPTION:
!   This module contains optimization procedures used to optimize ion positions
!   given the forces.  
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
  USE CONSTANTS, ONLY : DP 
  USE CellInfo, ONLY: cell
  USE OUTPUT, ONLY : WrtOut
  USE OutputFiles, ONLY : outputUnit
  USE REPORT, ONLY: GeometryMinimizerReportHeader
  USE REPORT, ONLY: GeometryMinimizerReportFooter
  USE REPORT, ONLY: GeometryMinimizerReportSteps
  USE IonOptimizers, ONLY: forceCutoff
  USE IonOptimizers, ONLY: watch
  USE IonOptimizers, ONLY: watch2
  USE TIMER, ONLY : TimerStart
  USE TIMER, ONLY : TimerStop

  IMPLICIT NONE

CONTAINS


SUBROUTINE IonLBFGS(RhoOptimizer,rho,energy,forces,frozenIon)
!------------------------------------------------------------------------------
! DESCRIPTION:
!  The ion relaxation is done by limited-memory BFGS method.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/20/2009  Created (Chen Huang)
!------------------------------------------------------------------------------

  USE MathFunctions, ONLY: Inverse
  USE MPI_Functions, ONLY: message

  IMPLICIT NONE
                          !>> EXTERNAL VARIABLES <<!

  EXTERNAL RhoOptimizer 
  ! Density optimiation routine
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho 
  ! charge density in real space
  !
  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy
  ! energy array
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: forces
  ! First index ion index
  ! Second index is x/y/z, 
  ! 3rd index: kinds of forces 1 for total force
  !
  LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: frozenIon
  !
  !

                     ! >>>>>> INTERNAL VARIABLES <<<<<<<< !

  INTEGER :: iter = 0
  ! index for interation
  !
  INTEGER :: nMax = 1000
  !
  !
  INTEGER :: ia
  !
  !
  LOGICAL :: relaxDone
  !
  !
  REAL(KIND=DP) :: cart(cell%numIon,3) 
  ! cartesian coordinates for atoms
  !
  REAL(KIND=DP) :: cart2frac(3,3)   
  ! for coordinate conversion
  !
  REAL(KIND=DP) :: frac2cart(3,3)
  ! for coordinate conversion
  !
  REAL(KIND=DP) :: maxForce
  ! maximal force
  !
  
  !!L-BFGS related 

  INTEGER  :: IPRINT(2),ICALL,IFLAG, N, M, LP, MP
  !
  REAL(KIND=DP) :: GTOL,STPMIN,STPMAX,EPS,XTOL, F
  !
  REAL(KIND=DP) :: G(3*cell%numIon),X(3*cell%numIon),DIAG(3*cell%numIon)
  !
  REAL(KIND=DP), ALLOCATABLE :: W(:) 
  ! working array for BFGS 
  !
  LOGICAL  :: DIAGCO
  ! diagonal hessian matrix
  !
  INTEGER :: step = 1    
  ! ionic step counter
  !
  EXTERNAL LB2
  !
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
  !

                         !!>> FUNCTIONS <<!!

  WRITE(message,*) "(Ion-relax): BFGS to be used."
  CALL WrtOut(6,message)

  watch2 = TimerStart()

  CALL GeometryMinimizerReportHeader
 
  relaxDone = .false.
  N = 3*cell%numIon             ! number of variables.
  M = 5                         ! number of corrections used in the BFGS update
  ALLOCATE( W(N*(2*M+1)+2*M) )  ! Working array

  ! BFGS parameters
  IPRINT(1)= -1  ! no output
  IPRINT(2)= 0
  DIAGCO= .FALSE. ! We do not wish to provide the diagonal matrices Hk0
  EPS  = 1.0D-5   ! We exit by monitoring force by ourselves
  XTOL = 1.0D-4   ! tolerance for atoms' movement in the line search
  ICALL= 0
  IFLAG= 0

  ! Matrices for coordinate conversion
  frac2cart = cell%cellReal
  cart2frac = Inverse(cell%cellReal)
  
  ! Initialize X (the coordinates in one-dimension)
  DO ia=1,cell%numIon
    cart(ia,:) = MATMUL(frac2cart,cell%ionTable(ia)%coord)
  ENDDO
  X = RESHAPE(cart,(/3*cell%numIon/))

  watch = TimerStart()
  
  ! Update F (total energy) and G ( -forces )
  CALL ComputeEnergyAndGradient(relaxDone)
  step = step + 1

  CALL WrtOut(6,"============================================= ")
  WRITE(message,'(a,I4,a,ES24.16,a,ES24.16,a)') & 
   "(Ion-Relax)BFGS Iter=", iter, & 
   ", totEnergy=", F, " (Ha), maxForce=", maxForce, " (Ha/bohr)"
  CALL WrtOut(6,message)
  CALL WrtOut(6,"============================================= ")
  CALL WrtOut(6,"")


  DO ! ========== Main loop ===========

    WRITE(outputUnit,*) " ------ ION (BFGS METHOD) STEP ", step, " ------"

    CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG) 
    
    ! Update coordinates
    cart = RESHAPE(X,(/cell%numIon,3/))
    DO ia=1,cell%numIon
      cell%ionTable(ia)%coord(:) = MATMUL(cart2frac,cart(ia,:))
    ENDDO

    watch = TimerStart()

    ! Update F (total energy) and G ( -forces )
    CALL ComputeEnergyAndGradient(relaxDone)

    ! Complete one BFGS Iteration
    IF ( IFLAG==0 ) THEN
      iter = iter + 1
      CALL WrtOut(outputUnit,"============================================= ")
      WRITE(message,'(a,I4,a,ES24.16,a,ES24.16,a)') & 
        "(Ion-Relax)BFGS Iter=", iter, & 
        ", totEnergy=", F, " (Ha), maxForce=", maxForce, " (Ha/bohr)"
      CALL WrtOut(outputUnit,message)
      CALL WrtOut(outputUnit,"============================================= ")
      CALL WrtOut(outputUnit,"")
    ENDIF

    ! Force drops below tolerence?
    IF ( relaxDone ) THEN
      WRITE(message,'(a,ES12.4,a,ES12.4,a)') & 
        "(Ion-Relax): Max Force=", maxForce, " < ",forceCutoff, & 
        " ,ion-relax is successful. Exit BFGS."
      CALL WrtOut(6,message)
      EXIT
    ENDIF

    IF ( IFLAG < 0 ) THEN 
      WRITE(message,*) & 
        "(Ion-Relax): BFGS negative FLAG, maxForce=", & 
        maxforce, ". Exit BFGS."
      CALL WrtOut(6,message) 
      EXIT
    ENDIF

    ICALL = ICALL + 1
    ! We allow at most nMax 
    ! evaluations of total energy and gradient
    IF (ICALL > nMax ) THEN
      WRITE(message,*) "(Ion-Relax) BFGS,  maxForce=",maxForce, & 
        ". Exceed maximum evaluation of total energy and gradient, exit BFGS."
      CALL WrtOut(6,message) 
      EXIT
    ENDIF

    step = step + 1

  ENDDO !End of main BFGS loop

  DEALLOCATE(W)

  CALL  GeometryMinimizerReportFooter(TimerStop(watch2))

  RETURN

CONTAINS


SUBROUTINE ComputeEnergyAndGradient(relaxDone)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Compute the total energy and forces then update F and G
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  USE RefreshIons, ONLY: RefreshIonTerms ! Subroutine that does some set up
  USE CalForces, ONLY : CalculateForces

  IMPLICIT NONE

  !! >> EXTERNAL PARAMETERS << !!

  LOGICAL :: relaxDone

  CALL RefreshIonTerms  ! Call each time ion positions change
  CALL RhoOptimizer
  CALL CalculateForces(rho, forces)
  IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP

  ! The total energy
  F = energy(1)
  ! The gradient
  G = -RESHAPE(forces(:,:,1),(/N/))

  maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))
  relaxDone = .false.
  relaxDone = maxForce < forceCutoff

  ! Print out information
  ! mohan add 2014-02-19
  CALL GeometryMinimizerReportSteps(step, TimerStop(watch), &
                                    forces, maxForce, &
                                    0.D0, .FALSE. , &
                                    MAXVAL(ABS(forces(:,:,1))), .TRUE.)

  RETURN

END SUBROUTINE ComputeEnergyAndGradient

END SUBROUTINE IonLBFGS

END MODULE IonOptBFGS

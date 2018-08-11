MODULE IonOptCG
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE IonOptCG
!     |_SUBROUTINE GradientOptimization (CG)
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
  USE OUTPUT, ONLY : WrtOut
  USE OutputFiles, ONLY : outputUnit
  USE Timer, ONLY : TimerStart
  USE Timer, ONLY : TimerStop 
  USE MPI_Functions

  IMPLICIT NONE
                     
CONTAINS


SUBROUTINE GradientOptimization(Optimizer, rho, energy, forces, frozenIon)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/10/2006  Created (GSH)
!   09/11/2008  Added wrapping ions back into periodic box (L. Hung)
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell

  USE MathFunctions, ONLY : Inverse

  USE IonOptimizers, ONLY: maxIonStep
  USE IonOptimizers, ONLY: watch
  USE IonOptimizers, ONLY: watch2
  USE IonOptimizers, ONLY: forceCutoff

  USE Report, ONLY: GeometryMinimizerReportHeader
  USE Report, ONLY: GeometryMinimizerReportFooter
  USE Report, ONLY: GeometryMinimizerReportSteps

  USE RefreshIons, ONLY: RefreshIonTerms
  USE CalForces, ONLY : CalculateForces          ! Calculate the forces.    

  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!

  EXTERNAL Optimizer 
  ! A subroutine that is called to optimize the electron 
  ! density relative to the ion positions
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho 
  ! The electron density, real space
  !
  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy
  !
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: forces
  ! Total forces.  Final index is 1 for total force,
  ! First is ion number, second direction (1,2,3 for x,y,z)
  !
  LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: frozenIon
  !

                        !>> INTERNAL VARIABLES <<!
  INTEGER :: i
  ! Counter to parse through the ion table.
  !
  INTEGER :: numEnergyLineMin
  !
  INTEGER :: step
  ! step number we are on
  !
  REAL(KIND=DP), DIMENSION(SIZE(cell%ionTable), 3) :: origCoord
  ! The last coordinates
  !
  REAL(KIND=DP), DIMENSION(SIZE(cell%ionTable), 3) :: fracStep
  ! step in fractional coordinates
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1), SIZE(forces,2)):: g, dirNext
  !
  REAL(KIND=DP) :: gamma
  REAL(KIND=DP) :: maxForce
  REAL(KIND=DP) :: fp
  REAL(KIND=DP) :: stp
  REAL(KIND=DP) :: ftol = 0.01_DP
  REAL(KIND=DP) :: gtol = 0.1_DP
  REAL(KIND=DP) :: xtol = 1.E-10_DP
  REAL(KIND=DP) :: stpmin = 1.E-100_DP
  REAL(KIND=DP) :: stpmax = 1000._DP
  REAL(KIND=DP) :: fpLast
  REAL(KIND=DP) :: checkpt1
  REAL(KIND=DP) :: checkpt2
  CHARACTER(len=60) :: task
  INTEGER,DIMENSION(2) :: isave
  REAL(KIND=DP), DIMENSION(13) :: dsave
  EXTERNAL DCSRCH
  REAL(KIND=DP) :: cart2frac (3,3)

                        !>> INITIALIZATION <<!

  watch2 = TimerStart()
  forces = 0._DP
  stp = 0._DP
  cart2frac = Inverse(cell%cellreal)

  CALL GeometryMinimizerReportHeader

  step = 1
  watch = TimerStart()

  ! Optimize initial density and calculate forces
  CALL Optimizer
  CALL CalculateForces(rho, forces)

  IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP
  maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)) )

  ! If we will leave the optForces loop immediately, need some initial output
  IF (maxForce < forceCutoff) THEN
    CALL GeometryMinimizerReportSteps(step, TimerStop(watch), &
      forces, maxForce, stp, .FALSE., MAXVAL(ABS(forces(:,:,1))), .TRUE.)
  ENDIF


                        !>> FUNCTION BODY <<!

  ! Loop to minimize the force
  optForces: DO

    WRITE(outputUnit,*) " ------ ION (GRADIENT METHOD) STEP ", step, " ------"

    IF (maxForce < forceCutoff) EXIT

    IF(step == 1) THEN
      dirNext =  forces(:,:,1)
      fp = SUM(-forces(:,:,1)*dirNext)
      stp = 5._DP
    ELSE
      gamma = MAX(SUM(forces(:,:,1)*(forces(:,:,1)-g))/SUM(g*g), 0.0_DP)
      dirNext =  forces(:,:,1) + gamma * dirNext
      fpLast = fp
      fp = SUM(-forces(:,:,1)*dirNext)
      stp = MAX(MIN(stp* fpLast / fp, stp*2._DP),stp*.25_DP)
    END IF

    g = forces(:,:,1)

    DO i = 1, 3
      origCoord(:,i) = cell%ionTable%coord(i)
    END DO

    ! Do the linesearch
    numEnergyLineMin = 0
    task = 'START'
    linesearch: DO
      
      CALL Dcsrch(stp,energy(1),SUM(-forces(:,:,1)*dirNext),&
                  ftol,gtol,xtol,task,stpmin,stpmax, &
                  isave,dsave)

      IF(task .EQ. 'FG') THEN
        CALL GeometryMinimizerReportSteps(step, TimerStop(watch), &
                                          forces, maxForce, &
                                          stp, .FALSE., &
                                          MAXVAL(ABS(forces(:,:,1))), .FALSE.)

        watch = TimerStart()
      ELSE
        EXIT
      END IF

      ! This is an experimental method to exit if things have gone to heck
      IF(MOD(numEnergyLineMin,5) == 4) THEN
        checkpt2 = checkpt1
        checkpt1 = energy(1)
        IF(numEnergyLineMin > 10 &
           .AND. ABS(checkpt2 - checkpt1)/SIZE(cell%ionTable) < 1.E-9_DP) THEN
          ! Attempt restart from steepest descent
          stp = 1._DP
          g = forces(:,:,1) ! Should force gamma to 0
          EXIT
        END IF
      END IF

      ! Takes a step 'stp'.  Here we actually change the ion positions.
      ! We could just loop over ions, but we suspect that vectorizing is 
      ! faster over the longer vector.  We don't right now wrap ions back
      ! into the cell that have fallen off the edge.
      ! We keep track of the last position we came from.
      do i = 1,cell%numIon
        fracStep(i,:) = MATMUL(cart2frac, stp*dirNext(i,:))
      enddo
      DO i = 1, 3
        cell%ionTable%coord(i) = MODULO(origCoord(:,i) + fracStep(:,i), 1._DP)
      END DO

      CALL RefreshIonTerms ! This must be called every time
                                            ! the ion positions are changed
      CALL Optimizer                ! Optimize the density
      CALL CalculateForces(rho, forces)
      IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP
      maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))
      numEnergyLineMin = numEnergyLineMin+1

    END DO linesearch

    CALL GeometryMinimizerReportSteps(step, TimerStop(watch), &
                                      forces, maxForce, &
                                      stp, .FALSE., &
                                      MAXVAL(ABS(forces(:,:,1))), .TRUE.)
    watch = TimerStart()

    step = step + 1

  END DO optForces

  CALL GeometryMinimizerReportFooter(TimerStop(watch2))

  RETURN

END SUBROUTINE GradientOptimization

END MODULE IonOptCG

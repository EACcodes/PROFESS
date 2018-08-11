MODULE IonOptQui
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE IonOptQui
!     |_SUBROUTINE QuickMinOptimization (QUI)
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
  USE TIMER, ONLY : TimerStart
  USE TIMER, ONLY : TimerStop 
  USE OUTPUT, ONLY : WrtOut
  USE OutputFiles, ONLY : outputUnit
  USE MPI_Functions

  IMPLICIT NONE

CONTAINS


SUBROUTINE QuickMinOptimization(Optimizer, rho, energy, forces, frozenIon)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine implements the QUICKMIN (projected Velocity Verlet) 
!   geometry optimization algorithm.  It assumes for now that all masses are 1.
!
!   The quickmin optimization procedure basically is a straight dynamics
!   model where you start with an initial (small) timestep, then accelerate,
!   keeping track of the velocities.  When the dot product of the force that
!   you started with and the current force becomes zero, then you back up
!   one step and start over.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   I'm not sure that I got all the criteria right ... I'm using the sum of the
!   dot products of all the directions as the criteria as to whether a step is 
!   bad.  This might not be (probably is not) the best thing to do.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!     Long time ago: Created by Greg Ho
!     June-2-2008  : Make the time_step_(n+1) = 2._DP*time_step_n during the 
!                    ion relaxation, this will accelerate the convergence a lot.
!                    If the time step is too big, it will be divided by 5.
!                    (Chen Huang)
!     09/11/08: Ions now wrap into periodic box (L. Hung)
!------------------------------------------------------------------------------
 
  USE IonOptimizers, ONLY: maxIonStep
  USE IonOptimizers, ONLY: timeStep
  USE IonOptimizers, ONLY: forceCutoff
  USE IonOptimizers, ONLY: watch, watch2

  USE CellInfo, ONLY: cell

  USE MathFunctions, ONLY : Inverse

  USE Report, ONLY : GeometryMinimizerReportHeader
  USE Report, ONLY : GeometryMinimizerReportFooter
  USE Report, ONLY : GeometryMinimizerReportSteps

  USE CalForces, ONLY : CalculateForces          ! Calculate the forces.    
  USE RefreshIons, ONLY: RefreshIonTerms

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
  ! differnt components of energy
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: forces
  ! Total forces.  Final index is 1 for total force,
  ! First is ion number, second direction (1,2,3 for x,y,z)
  !
  LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: frozenIon
  ! flag to determine whether atoms are frozen
  ! 

                        !>> INTERNAL VARIABLES <<!


  INTEGER :: i
  ! Counter to parse through the ion table.
  !
  INTEGER :: step
  ! step number we are on
  !
  REAL(KIND=DP), DIMENSION(cell%numIon, 3) :: velocity
  ! velocity of ions  
  !
  REAL(KIND=DP), DIMENSION(cell%numIon, 3) :: lastCoord
  ! The last coordinates
  !
  REAL(KIND=DP), DIMENSION(cell%numIon, 3) ::  minCoord
  !
  REAL(KIND=DP), DIMENSION(cell%numIon, 3) ::  fracStep
  ! step in fractional coordinates
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1), SIZE(forces,2), SIZE(forces,3)):: initForces
  ! The force in the direction we're traversing
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1), SIZE(forces,2), SIZE(forces,3)):: prevForces
  ! The previous forces
  !
  REAL(KIND=DP), DIMENSION(cell%numIon) :: mass
  !
  REAL(KIND=DP) :: maxForce
  !
  REAL(KIND=DP) :: minMaxForce
  !
  REAL(KIND=DP) :: bestForce
  !
  REAL(KIND=DP), SAVE :: lastTimeStep = -1._DP
  !
  LOGICAL :: overStep
  !
  REAL(KIND=DP), PARAMETER :: rhoOptTolMax = 100._DP
  ! Maximum factor changing tolerance for rho
  !
  REAL(KIND=DP), PARAMETER :: rhoOptTolMin = 1._DP      
  ! Minimum factor changing tolerance for rho
  !
  REAL(KIND=DP), DIMENSION(SIZE(energy)) :: oldEnergy
  !
  REAL(KIND=DP), DIMENSION(SIZE(energy)) :: testEnergy 
  ! Energies of line search at test point (totals at location 1)
  !
  REAL(KIND=DP) :: cart2frac(3,3) 
  !
  ! REAL(KIND=DP) :: rhoOptimizerTolAdjust = 1._DP 
  !

                        !>> INITIALIZATION <<!

  IF(rankGlobal == 0) THEN 
    WRITE(*,*) 'Quickmin Ionic Optimization Starting.'
  ENDIF

  watch2 = TimerStart()

  IF(lastTimeStep > 0._DP) timeStep = lastTimeStep
  overStep = .FALSE.
  forces = 0._DP
  mass = 1._DP
  
  cart2frac = Inverse(cell%cellReal)
  DO i = 1, 3
    minCoord(:,i) = cell%ionTable%coord(i)
  END DO

  CALL GeometryMinimizerReportHeader


                        !>> FUNCTION BODY <<!
  step = 1
  velocity = 0._DP   
  
  watch = TimerStart()

  WRITE(outputUnit,'(/A)')   " --------------------------------------------------"
  WRITE(outputUnit,'(A,I5,A,I5)') "  IONIC OPTIMIZATION (QUI METHOD), step=", step, "/", maxIonStep
  WRITE(outputUnit,'(A)')    " --------------------------------------------------"
  

  CALL Optimizer                ! Optimize density

  oldEnergy = energy


  CALL CalculateForces(rho, forces)  ! Compute initial forces

  IF(ALLOCATED(frozenIon)) THEN
    WHERE(frozenIon) forces(:,:,1)=0._DP
  ENDIF

  ! 1 here in forces stand for total force
  maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))
  bestForce = maxForce
  minMaxForce = maxForce
  initForces = forces

  ! Loop to minimize the force
  DO
    
    ! Print output for step
    WRITE(outputUnit,*)
    WRITE(outputUnit,'(A,I3,A,Es11.4,A,Es11.4,A)') &
        ' Quickmin: After step:',step, &
        ' maxForce=',maxForce,' (Force Tolerence=',forceCutoff,')'

    ! Print out information
    CALL GeometryMinimizerReportSteps(step, TimerStop(watch), &
                                      forces, maxForce, timeStep, &
                                      overStep, SUM(SQRT(velocity**2)), &
                                      .NOT.overStep, minCoord)

    ! it's better to put it just behind the GeometryMinimizerReportSteps
    ! functions, otherwise if we forget to add 1 to the step, it will
    ! have some mistake on opening geometry file next time.
    ! mohan note 2013-07-24
    step = step + 1

    IF( step > maxIonStep ) THEN
      WRITE(outputUnit,*) "Reach the allowed maximal ion steps = ", step
      EXIT
    ENDIF


    ! Exit if the magnitude of the forces are small enough
    IF (maxForce < forceCutoff) THEN
      WRITE(outputUnit,'(A)') " Force tolerence has been achieved."
      EXIT
    ENDIF

    watch = TimerStart()

    WRITE(outputUnit,'(/A)')   " --------------------------------------------------"
    WRITE(outputUnit,'(A,I5,A,I5)') "  IONIC OPTIMIZATION (QUI METHOD), step=", step, "/", maxIonStep
    WRITE(outputUnit,'(A)')    " --------------------------------------------------"

    WRITE(*,'(/A)')   " --------------------------------------------------"
    WRITE(*,'(A,I5,A,I5)') "  IONIC OPTIMIZATION (QUI METHOD), step=", step, "/", maxIonStep
    WRITE(*,'(A)')    " --------------------------------------------------"
  

    ! **EXPERIMENTAL** If we have taken 10 steps w/o any progress,
    !  decrease the timestep.
    IF(MOD(step,10) == 9) THEN
      IF (maxForce > minMaxForce*3._DP) THEN
        timeStep = timeStep * .5_DP
      ELSE
        minMaxForce = maxForce
      END IF
    END IF

    ! D = Do + Vot + 0.5 * a * t^2
    DO i = 1,size(cell%ionTable,1)
      fracStep(i,:) = MATMUL(cart2frac, timeStep * velocity(i,:) &
                      + 0.5_DP * timeStep**2 * forces(i,:,1))
    ENDDO

    ! Here we change the ion positions.
    ! We could just loop over ions, but we suspect that vectorizing is 
    ! faster over the longer vector.  We don't right now wrap ions back
    ! into the cell that have fallen off the edge.
    ! We save the last position we came from.
    DO i = 1, 3
      lastCoord(:,i) = cell%ionTable%coord(i)
      cell%ionTable%coord(i) = MODULO(cell%ionTable%coord(i) + fracStep(:,i), 1._DP)
    END DO

    ! Save forces
    prevForces = forces

    ! This is REALLY UGLY.
    ! Here, we adjust the tolerance of the density minimizer (fast convg?)
!    rhoOptimizerTolAdjust = MIN(MAX(rhoOptTolMin, maxForce/tol/5._DP), &
!                                rhoOptTolMax)

    ! Compute new forces and energy
    CALL RefreshIonTerms ! This must be called every time the 
                                          ! ion positions are changed
    CALL Optimizer            ! Optimize the density

    testEnergy=energy

    
    CALL CalculateForces(rho, forces)


    IF(ALLOCATED(frozenIon)) THEN
      WHERE(frozenIon) forces(:,:,1)=0._DP
    ENDIF

    ! calculate the maximal force
    maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))
    IF(maxForce < bestForce) THEN
      DO i = 1, 3
        minCoord(:,i) = cell%ionTable%coord(i)
      END DO
      bestForce = maxForce
    END IF

    ! Prepare for next iteration, depending on whether rho optimization was
    ! unsuccessful, whether the most recent step has made ions "overStep",
    ! or whether we proceed with standard quickMin.

    ! Change timestep and/or velocity parameters if "overstepped"
    IF (SUM(forces(:,:,1) * velocity) < 0._DP .OR. &
             testEnergy(1) > oldEnergy(1)) THEN

      ! Initial step must be valid.
      ! Also, no steps should allow energy to increase.
      IF(step == 1 .OR. testEnergy(1) > oldEnergy(1)) THEN

        overStep = .TRUE.
        timeStep = timeStep / 5._DP
        velocity = 0._DP ! Reset velocity if cutting timestep is not enough

        IF (rankGlobal==0) THEN
          WRITE(outputUnit,*) " Quickmin: Note, too big time step (over step), I have decreased it."
          WRITE(outputUnit,*) " Time step is ", timeStep
        ENDIF

        lastTimeStep = timeStep

        ! Restore the previous ion positions and forces
        forces = prevForces
        DO i=1,3
          cell%ionTable%coord(i) = lastCoord(:,i)
        END DO
        CYCLE

      ! If only a slight overstep (energy still decreases), make velocity 0 and
      ! proceed
      ELSE 
        overStep = .FALSE.
        oldEnergy = testEnergy
        velocity = 0._DP
      END IF

   ! Update the velocity as in standard quickmin for the next iteration
    ELSE
      ! Experimental - timestep acceleration
      timeStep = timeStep * 2._DP ! we increase the timestep here to accelerate
                                  ! the convergence 
      oldEnergy = testEnergy
      overStep = .FALSE.
      ! V = Vo + at projected in the directon of the forces
      velocity = SUM(forces(:,:,1) * velocity) * forces(:,:,1) / &
                 SUM(forces(:,:,1)**2) + forces(:,:,1) * timeStep
    ENDIF


  END DO

  CALL GeometryMinimizerReportSteps(step, TimerStop(watch), &
                                    forces, maxForce, timeStep, &
                                    overStep, SUM(SQRT(velocity**2)), &
                                    .NOT.overStep, minCoord)

  CALL GeometryMinimizerReportFooter(TimerStop(watch2))

  RETURN

END SUBROUTINE QuickMinOptimization


END MODULE IonOptQui

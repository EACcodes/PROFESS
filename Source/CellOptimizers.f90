MODULE CellOptimizers
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE CellOptimizers
!     |_SUBROUTINE NoOptimization
!     |_SUBROUTINE SteepestDecent
!       |_SUBROUTINE WriteLogFile
!
! DESCRIPTION:
!   This module contains optimization procedures used to optimize a cell
!   given the stress.  These subroutines are passed the function that it needs
!   to call to relax the ionic configuration and electron density after each
!   change of the lattice vector.  
!
!   The only module this functon knows about is Calculator, which is the 
!   module that returns stress terms given cell parameters.
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

  USE CONSTANTS, ONLY: DP 
  !
  USE CellInfo, ONLY: cellStruct 
  ! A cell type
  !
  USE RefreshCell, ONLY: RefreshCellSetup
  !
  USE CalStress, ONLY : CalculateStress  
  ! The subroutine that computes the stress.
  !
  USE RefreshIons, ONLY: RescaleDensity  
  ! Subroutine that maintains the integral of rho
  ! constant at the number of electrons.
  !
  USE OutputFiles, ONLY: outputUnit

  IMPLICIT NONE

  INTEGER :: maxCellStep = 100
  ! maximal cell relaxation steps
  !
  REAL(KIND=DP) :: tols = 5.0E-7_DP
  ! tolerence for stress
  ! should be 100 times smaller than the tol of force which 
  ! is generally 5e-5 a.u.
  !

CONTAINS


SUBROUTINE NoOptimization(Optimizer, rhoR, energy, stress, calcStress)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Does nothing.  Simply minimizes the geometry of the given cell dimentions,
!   using Calculator (which is a subroutine that is passed in). 
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
                        !>> EXTERNAL VARIABLES <<!

  EXTERNAL Optimizer 
  ! A subroutine that is called to optimize the ions 
  ! positions relative to the cell dimensions.
  !
  LOGICAL, INTENT(IN) :: calcStress       
  ! Whether or not to calculate the stress
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rhoR 
  ! The density of the system
  !
  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy     
  ! The energy of the system 
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(OUT) :: stress  
  ! The final answer, the stress
  !

                        !>> FUNCTION BODY <<!

  CALL Optimizer
  ! Calculate the stress if that was requested
  IF (calcStress) THEN
    CALL CalculateStress(rhoR, energy, stress)
  ENDIF

  RETURN

END SUBROUTINE NoOptimization


SUBROUTINE SteepestDecent(Optimizer, rhoR, energy, stress, axis)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Tries to minimize the stress by changing the lattice vectors. At the moment
!   the process is somewhat arbitrary and may or may not work well, but it is
!   simple.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Find optimal values of tol and dt, or create a variable timestep method.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE MathFunctions, ONLY: Inverse
  USE Output, ONLY: WrtOut
  USE OUTPUT, ONLY: PrintStress
  USE Output, ONLY: PrintGeometry
  USE Output, ONLY: celOutputFreq
  USE CellInfo, ONLY: RefreshLattice
  USE CellInfo, ONLY: cell

  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!

  EXTERNAL Optimizer  
  ! A subroutine that is called to optimize the ions positions.
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: rhoR  
  ! The real-space electron density
  !
  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy 
  ! The energies of the system 
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(OUT) :: stress            
  ! The final answer, the stress
  !
  INTEGER, INTENT(IN) :: axis 
  ! Direction of stress optimization
  ! 1 = x-direction only; 2 = y-direction only
  ! 3 = z-direction only; otherwise optimize in all dims

                        !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: tmpStress(3,3)
  !
  REAL(KIND=DP) :: oldEnergy
  ! Total energy at the previous step
  !
  REAL(KIND=DP) :: dt
  ! The arbitrary step to minimize the stress
  !
  REAL(KIND=DP), DIMENSION(3,3) :: newCell           
  ! Transformed cell including the rotation from stress. 
  !
  INTEGER :: i,j,k               
  ! Loop counters
  !
  CHARACTER (len=500) :: message
  !
  LOGICAL :: cellRelaxFlag = .TRUE.
  ! used for printing geometry
  !


                        !>> INITIALIZATION <<!
  dt = 10._DP
  i = 0
  oldEnergy = 0._DP

                        !>> FUNCTION BODY <<!

  SELECT CASE(axis)
  CASE(1) 
    CALL WrtOut(6,"Source/CellOptimizers.f90 >>> Only X direction is to be relaxed")
  CASE(2) 
    CALL WrtOut(6,"Source/CellOptimizers.f90 >>> Only Y direction is to be relaxed")
  CASE(3) 
    CALL WrtOut(6,"Source/CellOptimizers.f90 >>> Only Z direction is to be relaxed")
  CASE DEFAULT
    CALL WrtOut(6,"Source/CellOptimizers.f90 >>> Cell is to be FULLY relaxed.")
  ENDSELECT

  ! Loop to relax the cell
  DO
     
    i = i + 1
    IF( i > maxCellStep ) THEN
      WRITE(outputUnit,*) "Reach the allowed maximal cell steps = ", i 
      EXIT
    ENDIF

    ! print the geometry for cell relaxation
    ! default: celOutputFreq == -1
    IF((celOutputFreq>0 .AND. MOD(i,celOutputFreq)==0)) THEN
      CALL PrintGeometry(cellRelaxFlag, i)
    ENDIF

    CALL Optimizer

    CALL CalculateStress(rhoR, energy, stress)! Compute stress as needed.

    CALL PrintStress(stress)

    CALL WriteLogFile


    IF (axis == 1 .or. axis == 2 .or. axis == 3) THEN 
      tmpStress = stress
      stress = 0.d0
      stress(axis,axis) = tmpStress(axis,axis)
    ENDIF

    ! Is the stress minimized?
    IF (axis>0 .AND. axis<4) THEN
      IF (ABS(stress(axis,axis))<tols) EXIT
    ELSE
      IF (MAXVAL(ABS(stress))<tols) EXIT
    END IF

    ! Here is the minimization scheme in itself.
    IF (energy(1)>oldEnergy) THEN
      dt = dt/ 4._DP
    ELSE
      dt = dt * 1.2_DP
    END IF
    oldEnergy = energy(1)
    
    ! This may rotate the cell.
    newCell = cell%cellReal - cell%vol * dt * MATMUL(stress, TRANSPOSE(Inverse(cell%cellReal)))

    !Get rid of most changes if uniaxial deformation
    IF (axis>0 .AND. axis<4) THEN
      DO j = 1,3
        DO k = 1,3
          IF (k/=axis .OR. j/=axis) newCell(k,j) = cell%cellReal(k,j)
        END DO
      END DO
    END IF

    ! Refresh the cell, including the length and angles of cell
    CALL RefreshLattice(newCell)
    ! Update the properties that depends on the cell size
    CALL RefreshCellSetup
    ! Rescale the electron density
    CALL RescaleDensity(rhoR)

  END DO 

  CONTAINS

  SUBROUTINE WriteLogFile
  !-------------------------------------------
  ! DESCRIPTION:
  ! A simple subroutine to output messages 
  ! regarding the cell optimization process
  !
  ! CONDITIONS AND ASSUMPTIONS:
  !
  ! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
  !
  !-------------------------------------------
  ! REVISION LOG:
  !
  !-------------------------------------------

    CALL WrtOut(6,'')
    WRITE(message,'(A,i3,A,i3,A,G12.4,A,G16.8,A)') & 
       '; Cell optimization at step', i, ' / ', maxCellStep, & 
       ', max(stress):',maxval(abs(stress)), & 
       ', total energy:', energy(1), ';'

    WRITE(outputUnit,'(A,i3,A,G12.4,A,G16.8,A)') & 
       '; Cell optimization at step', i, & 
       ', max(stress):',maxval(abs(stress)), & 
       ', total energy:', energy(1), ';'

    CALL WrtOut(6,message)
    CALL WrtOut(6,'stress (Ha/bohr^3):')
    WRITE(message,'(A,3es12.4)') ' ', stress(1,:); CALL WrtOut(6,message)
    WRITE(message,'(A,3es12.4)') ' ', stress(2,:); CALL WrtOut(6,message)
    WRITE(message,'(A,3es12.4)') ' ', stress(3,:); CALL WrtOut(6,message)
    CALL WrtOut(6,'lattice (bohr):')
    WRITE(message,'(A,3es12.4)') ' ', cell%cellReal(:,1); CALL WrtOut(6,message)
    WRITE(message,'(A,3es12.4)') ' ', cell%cellReal(:,2); CALL WrtOut(6,message)
    WRITE(message,'(A,3es12.4)') ' ', cell%cellReal(:,3); CALL WrtOut(6,message)
    CALL WrtOut(6,'')
    CALL WrtOut(6,'')

    RETURN

  END SUBROUTINE WriteLogFile

END SUBROUTINE SteepestDecent


END MODULE CellOptimizers

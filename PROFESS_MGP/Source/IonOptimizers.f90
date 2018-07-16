MODULE IonOptimizers
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE IonOptimizers
!     |_SUBROUTINE NoOptimization
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
  USE Timer, ONLY : stopwatch

  IMPLICIT NONE

  REAL(KIND=DP) :: timeStep = .4_DP
  ! For QUI ion optimization method
  !
  REAL(KIND=DP) :: forceCutoff = 5e-5 
  ! According to Dr. Carter 1.E-2 ev/A is enough. 
  ! 1.E-2 eV/A is 1.94469057410952468E-4 hartree/bohr.
  !                  
  TYPE(stopwatch) :: watch, watch2 
  ! Timer objects
  !
  INTEGER :: cg_type = 2 
  ! 1=>  Polak.Ribieare (more robust)
  ! 2=>  Hager and Zhang
  ! 3=>  Dai and Yuan
  !
  INTEGER :: maxIonStep = 200 
  ! maximal allowed ion relaxation steps
  !

CONTAINS


SUBROUTINE NoOptimization(Optimizer, rho, forces, calcForces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This program leaves ion positioning unchanged.  It optimizes only electron
!   density (given the cell dimensions) using calcForces, a subroutine that is
!   passed in.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  
  USE CalForces, ONLY: CalculateForces

  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!

  EXTERNAL Optimizer 
  ! A subroutine that is called to optimize the electron 
  ! density relative to the ion positions.
  !
  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho 
  ! The electron density, real space
  !
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: forces  
  ! Total forces.  Final index is 1 for total force,
  ! First is ion number, second direction (1,2,3 for x,y,z)
  !
  LOGICAL :: calcForces 
  ! Whether to calculate the forces
  !

  CALL Optimizer
  
  IF(calcForces) THEN
    CALL CalculateForces(rho, forces)
  ENDIF

  RETURN

END SUBROUTINE NoOptimization

END MODULE IonOptimizers

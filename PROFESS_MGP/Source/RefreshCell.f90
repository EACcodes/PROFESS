MODULE RefreshCell 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RefreshCell
!     |_SUBROUTINE RefreshCellSetup
!
! DESCRIPTION:
! Refresh the cell setup,  
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
                              !<< GLOBAL >>
  USE CONSTANTS, ONLY: DP
 
  IMPLICIT NONE

CONTAINS

SUBROUTINE RefreshCellSetup
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is the routine that has to be run after each time the cell 
!   dimensions are changed. 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell
  USE CellInfo, ONLY: RefreshLattice

  USE EWALD, ONLY : EwaldSetup 
  USE CellInfo, ONLY: m1G, m2G, m3G
  USE IonElectronSpline, ONLY: iiSpline

  USE SetupKEDF, ONLY : KEDFRefresh    
  USE CellInfo, ONLY: kinetic

  USE RefreshIons, ONLY: RefreshIonTerms



  IMPLICIT NONE

                     ! >> FUNCTION BODY <<!

  ! Refresh cell volume, dV, lattice vector angles
  CALL RefreshLattice(cell%cellReal)

  ! Setting up for the Ion-Ion term calculators for Energy and forces, etc.
  CALL EwaldSetup(cell%cellReal, cell%ionTable, cell%elementTable, iiSpline, m1G, m2G, m3G)

  ! refresh parameters in KEDF if needed
  CALL KEDFRefresh(kinetic)

  ! reset density if needed 
  ! recalculate ion-ion energy 
  ! recalculate ion-electron potential
  CALL RefreshIonTerms()

  RETURN

END SUBROUTINE RefreshCellSetup


END MODULE RefreshCell

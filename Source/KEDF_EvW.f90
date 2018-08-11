MODULE KEDF_EvW
!----------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_EvW
!     |_SUBROUTINE Cal_EVT
!     |_SUBROUTINE Cal_EVC 
!
! DESCRIPTION:
! EvW Kinetic energy density functional (EvW-KEDF).
! A KEDF that for semiconductors. EvW-KEDF enhances the von Weizsacker KEDF 
! term. The enhancement factor is related to localized electron density.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
! J. Chem. Phys. 140, 18A531 (2014) Shin and E. A. Carter 
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012 ALGORITHM WRITTEN BY SHIN
! 2014 MODULE RECREATED BY MOHAN CHEN
!------------------------------------------------------------------------------

                              !<< GLOBAL >>

  USE Constants, ONLY: DP
  USE MPI_Functions

  IMPLICIT NONE

  REAL(KIND=DP) :: kloc = -100.0_DP   
  ! Parameters for EvW KEDF (Shin, mohan add 09-26-12)
  !
  REAL(KIND=DP) :: aloc = -100.0_DP   
  ! Parameters for EvW KEDF (Shin, mohan add 09-26-12)
  !
  REAL(KIND=DP) :: tolk = 1.0e-8 
  ! the difference of aloc between two extra steps 
  ! tolerence for extra loop of iterations to converge aloc
  ! for EvW KEDF (Shin, mohan add 10-01-12)
  !

CONTAINS


SUBROUTINE Cal_EVT(potential, rho, calcEnergy, locETable, optSqrt)
!------------------------------------------------------------------------------
! DESCRIPTION:
! EvW KEDF based on Want-Teter KEDF.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 
!------------------------------------------------------------------------------

  USE KEDF_TF, ONLY: TFPotential, TFEnergy
  USE KEDF_VW, ONLY: VWPotentialSqrtPlus 
  USE CellInfo, ONLY: cell, numSpin
  USE CellInfo, ONLY: n1G, n2G, n3G
  USE SYS, ONLY: bvac ! use vacuum treatment or not
  USE KEDF_WT, ONLY: WTPotentialPlus

  IMPLICIT NONE

  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(OUT) :: potential
  ! TF potential, last dimension is spin
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  ! charge density on real space
  ! 
  LOGICAL, INTENT(IN) :: calcEnergy
  ! calculate TF energy or not
  ! 
  REAL(KIND=DP), DIMENSION(9), INTENT(OUT) :: locETable 
  ! contains each energy term
  !
  LOGICAL, INTENT(IN) :: optSqrt  
  ! Take derivative relative to sqrt(rho) (instead of rho)
  !
  !! >> LOCAL VARIABLE << !!
  !
  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:) :: tempPot
  ! temp array
  !
  INTEGER :: allocateStatus           
  ! Checks that tables are allocated alright.
  !

  IF(numSpin==2) THEN
    WRITE(message,*) "EvW only works for non spin-polarized now"
    CALL Error(6, message)
  ENDIF

  AlLOCATE(tempPot(n1G, n2G, n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(message,*)'Error allocating tempPot array. Leaving.'
    CALL Error(6, message) 
  END IF

  potential = potential + SPREAD(TFPotential(rho(:,:,:,1)), 4, numSpin)*(1._DP-aloc/2._DP)
  IF (calcEnergy) THEN
    locETable(7) = TFEnergy(rho(:,:,:,1))*(1._DP-aloc/2._DP)
  ENDIF

  ! Compute the Wang-Teter term
  CALL WTPotentialPlus(rho(:,:,:,1), tempPot, calcEnergy, locETable(9), bvac)
  potential = potential + SPREAD(tempPot, 4, numSpin)*(1._DP-aloc/2._DP)

  IF (calcEnergy) THEN
    locETable(9)=locETable(9)*(1._DP-aloc/2._DP)
  ENDIF

  ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
  IF (optSqrt) THEN
    potential = 2._DP * SQRT(rho) * potential
  ENDIF

  ! Von Weiszacker (vW)
  CALL VWPotentialSqrtPlus(SQRT(rho(:,:,:,1)), tempPot, calcEnergy, locETable(8))

  IF(calcEnergy) THEN
    locETable(8) = locEtable(8)*(1._DP+aloc) 
  ENDIF

  ! Chain rule: dE/d(rho) = dE/d(sqrt(rho)) / (2*sqrt(rho))
  IF (.NOT. optSqrt) THEN
    tempPot = tempPot/(2._DP*SQRT(rho(:,:,:,1)))
  ENDIF
             
  potential = potential + SPREAD(tempPot, 4, numSpin)*(1.0_DP+aloc)


  DEALLOCATE(tempPot)

  RETURN

END SUBROUTINE Cal_EVT



SUBROUTINE Cal_EVC(potential, rho, calcEnergy, locETable, optSqrt)
!------------------------------------------------------------------------------
! DESCRIPTION:
! EvW KEDF based on WGC KEDF.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE KEDF_TF, ONLY: TFPotential, TFEnergy
  USE KEDF_VW, ONLY: VWPotentialSqrtPlus 
  USE CellInfo, ONLY: cell, numSpin
  USE CellInfo, ONLY: n1G, n2G, n3G
  USE SYS, ONLY: bvac ! use vacuum treatment or not
  USE KEDF_WGC, ONLY: WGCPotentialPlus

  IMPLICIT NONE

  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(OUT) :: potential
  ! TF potential, last dimension is spin
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  ! charge density on real space
  ! 
  LOGICAL, INTENT(IN) :: calcEnergy
  ! calculate TF energy or not
  ! 
  REAL(KIND=DP), DIMENSION(9), INTENT(OUT) :: locETable 
  ! contains each energy term
  !
  LOGICAL, INTENT(IN) :: optSqrt  
  ! Take derivative relative to sqrt(rho) (instead of rho)
  !
  !! >> LOCAL VARIABLE << !!
  !
  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:) :: tempPot
  ! temp array
  !
  INTEGER :: allocateStatus           
  ! Checks that tables are allocated alright.
  !

  IF(numSpin==2) THEN
    WRITE(message,*) "EvW KEDF only works for non spin-polarized now"
    CALL Error(6, message)
  ENDIF

  AlLOCATE(tempPot(n1G, n2G, n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(message,*)'Error allocating tempPot array. Leaving.'
    CALL Error(6, message) 
  END IF

  potential = potential + SPREAD(TFPotential(rho(:,:,:,1)), 4, numSpin)*(1._DP-aloc/2._DP)

  IF (calcEnergy) THEN
    locETable(7) = TFEnergy(rho(:,:,:,1))*(1._DP-aloc/2._DP)
  ENDIF

  CALL WGCPotentialPlus(rho(:,:,:,1), tempPot, calcEnergy, locETable(9), bvac)

  IF (calcEnergy) THEN
    locETable(9) = locEtable(9)*(1._DP-aloc/2._DP)
  ENDIF

  potential = potential + SPREAD(tempPot, 4, numSpin)*(1._DP-aloc/2._DP)

  ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
  IF (optSqrt) THEN
    potential = 2._DP * SQRT(rho) * potential
  ENDIF

  CALL VWPotentialSqrtPlus(SQRT(rho(:,:,:,1)), tempPot, calcEnergy, locETable(8))
  
  IF (calcEnergy) THEN
    locETable(8) = locEtable(8)*(1._DP+aloc)
  ENDIF

  ! Chain rule: dE/d(rho) = dE/d(sqrt(rho)) / (2*sqrt(rho))
  IF (.NOT. optSqrt) tempPot = tempPot/(2._DP*SQRT(rho(:,:,:,1)))

  potential = potential + SPREAD(tempPot, 4, numSpin)*(1._DP+aloc) 

  DEALLOCATE(tempPot)

  RETURN

END SUBROUTINE Cal_EVC


END MODULE KEDF_EvW

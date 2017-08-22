MODULE RefreshIons 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Calculator
!     |_SUBROUTINE RescaleDensity
!     |_SUBROUTINE RefreshIonTerms
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
                              !<< GLOBAL >>
  USE CONSTANTS, ONLY: DP
  USE MPI_Functions

  IMPLICIT NONE

  LOGICAL :: trashPreviousRho = .FALSE.   
  ! Reset electron density to uniform every step
  !

CONTAINS


SUBROUTINE RescaleDensity(rhoR)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine ensures that the density still integrates to the number of 
!   electrons after the cell vectors have been altered, as a change in volume
!   will change the total number number of electrons if the density is kept.
!   This function is necessary for stress optimization.
!
! CONDITIONS AND ASSUMPTIONS:  Only rescales properly for non-padded array
!   i.e. periodic boundaries.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell, numSpin ! cell information in it
  USE CellInfo, ONLY: m1G, m2G, m3G, m123G ! global real space FFT dimension
  USE Sys, ONLY: magmom  ! magnetic moments

  IMPLICIT NONE

                  !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: rhoR 
  ! Electronic density in real-space
  !

                  !>> INTERNAL VARIABLES <<! 

  INTEGER :: numEle                  
  ! number of electrons in the system
  !
  REAL(KIND=DP) :: sumrho = 0.0_DP
  ! sum of total rho
  ! 
  REAL(KIND=DP) :: sumrho0(2)
  ! sum rho for each spin
  !
  REAL(KIND=DP) :: scalingFactor       
  ! used to keep the number of electrons constant.
  !
  INTEGER :: isp
  ! spin index
  !

                  !>> INITIALIZATION <<!  

  CALL TITLE("RefreshIon::RescaleDensity")

  ! The cell contains most of the info we need. 
  numEle = NINT(SUM(cell%elementTable%chargeTot))

                  ! >> FUNCTION BODY <<!

  ! Now, we normalize the density to the number of electrons.
  IF(numspin == 1) THEN
    sumrho = SUM(rhoR) * cell%dv
    CALL ReduceRealLevel1(sumrho)
  ELSE IF(numspin == 2) Then
    DO isp = 1,2
      rhoR(:,:,:,isp) = 0.5d0*(rhoR(:,:,:,isp) - ((-1)**isp) * magmom/ cell%vol )
      sumrho0(isp) = SUM(rhoR(:,:,:,isp)) * cell%dv
      CALL ReduceRealLevel1( sumrho0(isp) )
      WRITE(outputUnit,*) "Electron number for spin ", isp, " is ", sumrho0(isp)
    ENDDO
    sumrho = sumrho0(1) + sumrho0(2)
  ENDIF
    
  scalingFactor = REAL(numEle,kind=DP)/sumrho

  !WRITE(outputUnit,*) "sumrho=", sumrho
  !WRITE(outputUnit,*) "cell%dv=", cell%dv
  !WRITE(outputUnit,*) "numEle=", numEle
  !WRITE(outputUnit,*) "scalingFactor=", scalingFactor

  rhoR = rhoR * scalingFactor

  WRITE(outputUnit,'(A,F10.6,A,I10)') " Rescale the total density from ", sumrho * cell%dV, " to ", numEle 

!! DEBUG
!!  WRITE(message,*) " Rescale density, totQ=", numEle
!!  CALL WrtOut(6,message)
!! END DEBUG

  RETURN

END SUBROUTINE RescaleDensity


SUBROUTINE RefreshIonTerms
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine will perform all the necessary mopping up procedures
!   to execute whenever the ion positions are changed.
!   Currently, this includes:
!     1.  Recalculating the ion-ion energy
!     2.  Recalculating the ionic potential.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE SYS, ONLY: rhoR

  USE NearestDistance, ONLY: CheckNearestDistanceAtoms
  USE Output, ONLY : QUIT

  USE Ewald, ONLY: IonIonEnergy
  USE Ewald, ONLY: EwaldEnergy
  USE IonElectronSpline, ONLY: iiSpline
  USE CellInfo, ONLY: cell

  USE CellInfo, ONLY: n1G, n2G, n3G
  USE CellInfo, ONLY: k1G, k2G, k3G
  USE IonElectronSpline, ONLY: ionPotReal
  USE IonElectronSpline, ONLY: ieSpline
  USE IonElectron, ONLY: IonElectronPotentialRecip ! Function that returns the
  USE IonElectronSpline, ONLY: IonElectronPotRecipSpline
  USE Fourier_NEW, ONLY: FFT_NEW, FFT_STD_STATE

  USE KEDF_DenDec, ONLY: do_den_dec
  USE KEDF_DenDec, ONLY: SetupCoreDensity

  IMPLICIT NONE

                  !>> INTERNAL VARIABLES <<! 

  REAL(KIND=DP) :: minDistance
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ionPotRecip
  INTEGER :: allocateStatus

                  !>> INITIALIZATION <<!

  CALL Title("Calculator::RefreshIonTerms")
  CALL StartClock('IonRefresh')

  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           REFRESH ION POSITIONS"

  ! Reset electron density to uniform if specified
  IF (trashPreviousRho .EQV. .TRUE.) THEN
    ! set to uniform electron gas
    rhoR=SUM(cell%elementTable%chargeTot)/cell%vol
    CALL RescaleDensity(rhoR)
  END IF

  ! Check how close the nearest neighbors are
  IF (CheckNearestDistanceAtoms(minDistance)<0) THEN 
    message="  CheckNearestDistanceAtoms < 0"
    CALL QUIT(message)
  ENDIF

                  !>> FUNCTION BODY <<!

  ! 1) Recalculate the ion-ion energy
  ! Periodic boundary condition
  ionIonEnergy = EwaldEnergy(cell%cellReal, cell%ionTable, cell%elementTable, iiSpline)

  ! we don't need to do this for the first time,
  ! we use DEALLOCATE --> ALLOCATE for general usage of this subroutine,
  IF(ALLOCATED(ionPotReal)) DEALLOCATE(ionPotReal)

  ALLOCATE(ionPotReal(n1G, n2G, n3G),stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*)'Error allocating ionPotReal. Leaving.'
    STOP
  END IF

  ALLOCATE(ionPotRecip(k1G, k2G, k3G),stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*)'Error allocating ionPotRecip. Leaving.'
    STOP
  END IF

  ! 2) Recalculate the ionic potential
  IF (ieSpline .EQV. .TRUE.) THEN
    ionPotRecip = IonElectronPotRecipSpline(cell%ionTable, cell%elementTable)
  ELSE
    ionPotRecip = IonElectronPotentialRecip(cell%cellReal, cell%ionTable, cell%elementTable)
  END IF

  CALL FFT_NEW(FFT_STD_STATE, ionPotRecip, ionPotReal)

  DEALLOCATE(ionPotRecip)

  IF(do_den_dec == 1) THEN
    ! 3) for density decomposition method
    ! this is used to renew the density for the next ionic step,
    ! we should have a better plan to do this,
    ! noted by mohan 2013-09-03
    CALL SetupCoreDensity
    WRITE(outputUnit,*) "Make the core density one is in RefreshIonPositions"
  ENDIF

  CALL StopClock('IonRefresh')

  RETURN

END SUBROUTINE RefreshIonTerms


END MODULE RefreshIons

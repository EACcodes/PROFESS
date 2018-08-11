MODULE CalStress
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE CalStress
!     |_SUBROUTINE CalculateStress
!
! DESCRIPTION:
! Calculate the stress for the system, including KEDF stress,
! Hartree stress, exchange-correlation stress (LDA or PBE),
! ion-electron stress (spline or not), ion-ion stress (spline or not)
!
! CONDITIONS AND ASSUMPTIONS:
! Here is just the driver for the calculating each part of stress.
! Each stress is still calculated from each seperate module.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! before 2014: exist in Calculator.f90 in PROFESS 2.0
! 2014: mohan cretate CalStress module
!------------------------------------------------------------------------------

                              !<< GLOBAL >>

  USE Constants, ONLY: DP 
  USE Constants, ONLY: auToGPa 
  USE IonElectronSpline, ONLY: iiSpline
  USE IonElectronSpline, ONLY: ieSpline
  USE CellInfo, ONLY: cell, numSpin
  USE Output, ONLY: QUIT
  USE CellInfo, ONLY: kinetic
  USE MPI_Functions, ONLY: ReduceRealLevel1 
  USE CellInfo, ONLY: k1G, k2G, k3G
  USE OutputFiles

  IMPLICIT NONE

CONTAINS


SUBROUTINE CalculateStress(rhoR, energy, stress)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine computes the stress tensor for the current cell, ion table and
!   density. It presupposes that the electronic energy has been computed for 
!   the current density which, in current state of affairs, is always the case.
!
! CONDITIONS AND ASSUMPTIONS:
!   We're assuming that the modules we're using that actually calculate the 
!   various terms return ONLY UPPER TRIANGULAR parts of stress tensors that
!   must be symmetrized at the end. So at the end, we symmetrize the stress
!   tensor.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! 
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/05/2004  Subroutine created.  (Vincent Ligneres)
!   02/25/2004  Revised to make the individual routines return a full tensor.
!               (GSH)
!------------------------------------------------------------------------------

  USE IonElectron, ONLY: IonElectronStress ! The ion-electron stress function
  USE IonElectronSpline, ONLY: IonElectronStressSpline ! The ion-electron stress

  USE Hartree, ONLY: JStress     ! The coulombic stress function

  USE XC_PBE, ONLY: PBEStress    ! The PBE ex-corr stress function
  USE XC_LDA, ONLY: exchangeCorrelation ! type of exchange correlation
  USE XC_LDA, ONLY: LDAStress      ! The LDA exchange-correlation stress function 
  USE XC_LDA, ONLY: LSDAStress     ! The SPIN LDA exchange-correlation stress function

  USE KEDF_TF, ONLY: TFStress      ! The Thomas-Fermi stress function
  USE KEDF_VW, ONLY: VWStress      ! The Von Weisaker stress function
  USE KEDF_WT, ONLY: WTStress      ! The Wang-Teter stress function
  USE KEDF_WGC, ONLY: WGCStress    ! The Wang-Govind-Carter stress function
  USE KEDF_HC10, ONLY: IntStress             ! The HC Stress function

  USE Ewald, ONLY : EwaldStress

  USE CellInfo, ONLY: m123G

  USE Fourier_NEW, ONLY: FFT_STD_STATE
  USE Fourier_NEW, ONLY: FFT_NEW
  USE Fourier, ONLY: FFT

  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN), TARGET :: rhoR 
  ! The electron density in real-space, spin depend.
  !
  REAL(KIND=DP), DIMENSION(:) :: energy 
  ! A table of energies that are CORRECT for the given system configuration
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(OUT) :: stress 
  ! The 3x3 2-tensor, the answer
  !
 
  !>> INTERNAL VARIABLES <<!

  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G,SIZE(rhoR,4)) :: rhoRecip_SI
  ! Electron density in recip space, spin-independ.
  !
  INTEGER :: a, b    
  ! Stress tensor components
  !
  INTEGER :: isp
  ! spin
  !
  REAL(KIND=DP) :: tmp1, tmp2

                    !>> INITIALIZATION <<!
  CALL TITLE("CalStress")
  CALL StartClock("CalStress")

  stress = 0._DP

  ! Store the density independent real and reciprocal space densities locally

  DO isp = 1, numSpin
    CALL FFT_NEW(FFT_STD_STATE, rhoR(:,:,:,isp), rhoRecip_SI(:,:,:,isp))
  ENDDO

  tmp1 = REAL(numSpin,KIND=dp)
  tmp2 = 1.d0/tmp1


                    !>> FUNCTION BODY <<!
  
  !============== KINETIC STRESS ===============
  SELECT CASE(kinetic)

    CASE(1) 
      ! Thomas-Fermi (TF)
      stress = stress + TFStress(cell%Vol, energy(7))

    CASE(2) 
      ! Von Weiszacker (vW)
       DO isp = 1, numSpin
          stress = stress + tmp2*VWStress(FFT(SQRT(tmp1*rhoR(:,:,:,isp))), m123G)
       ENDDO

    CASE(3) 
      ! TF + vW
      stress = stress + TFStress(cell%vol, energy(7))
       DO isp = 1, numSpin
          stress = stress &
               + tmp2*VWStress(FFT(SQRT(tmp1*rhoR(:,:,:,isp))), m123G)
       ENDDO

    CASE(4)
      ! TF + vW + Wang-Teter(WT)
       stress = stress + TFStress(cell%vol, energy(7))
       DO isp = 1, numSpin
          stress = stress &
               + tmp2*VWStress(FFT(SQRT(tmp1*rhoR(:,:,:,isp))), m123G) &
               + tmp2*WTStress(tmp1*rhoR(:,:,:,isp), energy(9))
       ENDDO

    CASE(5)
      ! TF + vW + Wang-Govind-Carter(WGC)
      stress = stress + TFStress(cell%vol, energy(7))
      DO isp = 1, numSpin
        stress = stress + tmp2 * VWStress(FFT(SQRT(tmp1*rhoR(:,:,:,isp))), m123G) 
        stress = stress + tmp2 * WGCStress(rhoR(:,:,:,isp))
      ENDDO

    CASE(11)
      ! TF + vW + HC
       stress = stress + TFStress(cell%vol, energy(7))
       DO isp = 1, numSpin
          stress = stress &
               + tmp2*VWStress(FFT(SQRT(tmp1*rhoR(:,:,:,isp))), m123G) &
               + tmp2*intStress(cell%vol, tmp1*rhoR(:,:,:,isp))
       ENDDO

    CASE DEFAULT                    ! This case should never occur.

      message=" The KEDF selected is not implemented for stress calculation. Leaving."
      WRITE(errorUnit,*) " kinetic=", kinetic
      WRITE(errorUnit,*) message
      CALL QUIT(message)

  END SELECT

  !-----------
  ! FOR TEST 
  !-----------
  !WRITE(*,*) " stress for Kinetic (Gpa): "
  !WRITE(*,'(3F12.6)') Stress(1,:)*auToGPa
  !WRITE(*,'(3F12.6)') Stress(2,:)*auToGPa
  !WRITE(*,'(3F12.6)') Stress(3,:)*auToGPa

  
  !------------------------------
  ! Add the ion-electron stress
  !------------------------------
  IF(ieSpline) THEN ! B-spline approximation for ion-electron stress
     stress = stress + &
          IonElectronStressSpline(SUM(rhoRecip_SI,4), cell%ionTable, &
             cell%elementTable, energy(3), SIZE(rhoR,1))
  ELSE ! Exact calculation
     stress = stress + &
          IonElectronStress(SUM(rhoRecip_SI,4), cell%ionTable, cell%elementTable, &
             cell%cellReal, energy(3))
  END IF

  ! Add the coulombic stress
  stress = stress + JStress(cell%vol, SUM(rhoRecip_SI,4), energy(4))

  !-----------------------------------------
  ! Add the exchange-correlation stress
  !-----------------------------------------
  SELECT CASE (exchangeCorrelation)
    CASE(1) ! LDA 
      IF(numSpin == 1) THEN
        stress = stress + LDAStress(cell%vol, rhoR, energy(5))
      ELSEIF(numSpin == 2) THEN
        stress = stress + LSDAStress(cell%vol, rhoR, energy(5))
      ENDIF
    CASE(2) ! PBE
      stress = stress + PBEStress(rhoR(:,:,:,1), rhoRecip_SI(:,:,:,1))
    CASE DEFAULT ! This case should never occur.
    message="The XCF selected is not implemented. Leaving."
    CALL QUIT(message)
  END SELECT

  ! The previous stresses are different on each local processor
  ! We now sum the stresses to get the total stress (apart from ion-ion
  ! contributions, which will be added by calling EwaldStress)
  CALL ReduceRealLevel1(stress, SIZE(stress,1), SIZE(stress,2) )


  ! Add the ion-ion (Ewald) stress
  stress = stress + EwaldStress(cell%cellReal, cell%ionTable, cell%elementTable, iiSpline)

  !------------------------------
  ! Symmetrize the stress tensor
  !------------------------------
  DO b=1,3
    DO a=b+1,3
      stress(a,b) = stress(b,a) 
    END DO
  END DO

  CALL StopClock("CalStress")
  RETURN

END SUBROUTINE CalculateStress

END MODULE CalStress


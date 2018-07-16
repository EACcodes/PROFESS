MODULE KEDF_MGP
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_MGP
!     |_SUBROUTINE MGPPotentialPlus 
!     |_SUBROUTINE MGPPotPlus 
!   
! DESCRIPTION:
! Mi-Genova-Pavanello KEDF (2018)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
! REFERENCES: 
! Mi, Wenhui, Alessandro Genova, and Michele Pavanello. 
!"Nonlocal kinetic energy functionals by functional integration." 
!The Journal of chemical physics 148.18 (2018): 184107. 
!------------------------------------------------------------------------------
! REVISION LOG:
!------------------------------------------------------------------------------

  USE Constants,      ONLY: DP
  USE OutputFiles,    ONLY: outputUnit
  USE KEDF_WTkernel,  ONLY: keKernel
  USE KEDF_WGC

  IMPLICIT NONE

CONTAINS

SUBROUTINE MGPPotentialPlus(rho, potential, calcEnergy, energy,vacCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
! An interface to decide the calculation of MGP kinetic potential
! and the energy.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!------------------------------------------------------------------------------

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN)  :: rho        ! Electron density in real space (spin-independent)
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential  ! The Mi-Genova-Pavanello contribution to the kinetic energy.
  LOGICAL, INTENT(IN) :: calcEnergy
  LOGICAL, INTENT(IN) :: vacCutoff
  REAL(KIND=DP), INTENT(OUT) :: energy

  CALL MGPPotPlus(rho, potential, calcEnergy, energy, vacCutoff)

  RETURN

END SUBROUTINE MGPPotentialPlus


SUBROUTINE MGPPotPlus(rho, potential, calcEnergy, energy,vacCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine calculates the Mi-Genova-Pavanello kinetic potential 
!
! GLOBAL/MODULE VARIABLES CHANGED:
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

  IMPLICIT NONE


                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho        ! Electron density in real space (spin-independent)
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential ! The Mi-Genova-Pavanello contribution to the kinetic energy.
  LOGICAL, INTENT(IN) :: calcEnergy
  LOGICAL, INTENT(IN) :: vacCutoff
  REAL(KIND=DP), INTENT(OUT) :: energy

                     !>> INTERNAL VARIABLES <<!

! auxiliary vector in G space
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),      SIZE(rho,2),      SIZE(rho,3))       :: cutf,rhofs 

                       !>> INITIALIZATION <<!

  IF( ALLOCATED(keKernel) .EQV. .FALSE. ) THEN
    WRITE(*, *) "keKernel is not allocated."
    !CALL Error(6, message)
    stop
  ENDIF

  IF (SIZE(keKernel,4) .ne. 3) THEN
    WRITE(*, *) "keKernel is not allocated properly."
   ! CALL Error(6, message)
   stop
  ENDIF

  IF (vacCutoff .EQV. .TRUE.) THEN
    cutf = CutoffFunc(rho)
  END IF

  rhofs = 0.0
  rhofs = rho**(5._DP/6._DP) 
  energy = 0._DP

 ! Calculate the potential
  IF (vacCutoff) THEN
    potential = FFT(FFT(cutf*rhofs) * keKernel(:,:,:,1))
  ELSE
    potential = FFT(FFT(rhofs) * keKernel(:,:,:,1))
  END IF
  potential = potential * rho**(-1.0d0/6.0d0)
!=================================================================
  IF (calcEnergy) then
    energy = (3._DP/5._DP)*SUM(rho * potential)
  ENDIF
!=============================================================
END SUBROUTINE MGPPotPlus

END MODULE KEDF_MGP

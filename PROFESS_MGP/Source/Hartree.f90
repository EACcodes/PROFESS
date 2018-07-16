MODULE Hartree 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Hartree 
!     |_FUNCTION SqrtPreconditioner
!     |_FUNCTION JPotentialPlus
!     |_FUNCTION JStress
!
! DESCRIPTION:
!   Calculates the potential contibutions to the forces, stress, and energies.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2003-12-15  File Created by consolidating existing energy files (GSH) 
!   2003-12-18  Changed SUBROUTINE to FUNCTIONS
!   2004-01-07  Changed qTable and qMask and most instances of rhoR to 
!               spin-independent forms.
!   2004-02-05  Added nine stress functions. (VLL)
!   2007-12-14  Added Choly-Kaxiras methods for ion-electron terms.  (LH)
!   2013-12-09  Remove some redundant calculations about qVectors. (Mohan)
!
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!
  USE CONSTANTS, ONLY: DP
  USE CONSTANTS, ONLY: PI
  USE PlaneWave, ONLY: qTable
  USE PlaneWave, ONLY: qMask
  USE MathFunctions, ONLY: Vecmul
  USE PlaneWave, ONLY: WrapQ
  USE Fourier, ONLY: FFT ! The Fast Fourier Transform routine.

  IMPLICIT NONE
  
CONTAINS


FUNCTION SqrtPreconditioner(rhoRecip)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the coulombic potential in real-space.  The only thing special here
!   is that the q=0 term must be set to 0 before the final FFT.  Otherwise it 
!   will blow up.  Luckily this term is cancelled by the coulombic part of the
!   ion-electron term at g=0.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! The dimension of SqrtPreconditioner is dangerous to start from 0.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: k1G, k2G, k3G, k3Goff

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(0:,0:,0:), INTENT(IN) :: rhoRecip  ! Electron density in reciprocal space (spin-neutral)
  REAL(KIND=DP), DIMENSION(0:SIZE(rhoRecip,1)-1, 0:SIZE(rhoRecip,2)-1, &
                           0:SIZE(rhoRecip,3)-1, 1) :: &
    SqrtPreconditioner           ! The Coulombic potential.


                         !>> INTERNAL VARIABLES <<!
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: jPotentialRecip ! Electron density in reciprocal space (spin-neutral)

                          !>> INITIALIZATION <<!
                           !>> FUNCTION BODY <<!

   ! Assign the reciprocal space potential to a temporary array.  Sorry, I
   ! thought about it a lot but couldn't figure out how to do it without one.
   jPotentialRecip = FFT(rhoRecip) / qTable**2

   ! Set the q=0 term to zero.
   IF (k3Goff==0) THEN
     jPotentialRecip(1,1,1) = (0._DP,0._DP)
   ENDIF

   ! FFT to real space.
   SqrtPreconditioner(:,:,:,1) = FFT(jPotentialRecip)

END FUNCTION SqrtPreconditioner


SUBROUTINE JPotentialPlus(rhoReal, potential, calcEnergy, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the coulombic potential in real-space.  The only thing special here
!   is that the q=0 term must be set to 0 before the final FFT.  Otherwise it 
!   will blow up.  Luckily this term is cancelled by the coulombic part of the
!   ion-electron term at g=0.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   09/05/2008  Function Created following JPotential and JEnergy (L. Hung)
!
!------------------------------------------------------------------------------
  USE CellInfo, ONLY: k1G, k2G, k3G, k3Goff

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoReal        ! Spin neutral electron density in real space
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential      ! Calculated Hartree potential
  LOGICAL, INTENT(IN) :: calcEnergy     ! Whether to calculate energy
  REAL(KIND=DP), INTENT(OUT) :: energy

                         !>> INTERNAL VARIABLES <<!

  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: jPotentialRecip

                          !>> INITIALIZATION <<!

!   CALL StartClock('JPotPlus');

   energy = 0._DP
   jPotentialRecip = FFT(rhoReal) ! Store rhoRecip temporarily to save memory


  IF (k3Goff==0) THEN      ! qTable(1,1,1)=0, so we cannot divide jPotentialRecip
    qTable(1,1,1) = 1._DP! by it.  So we give qtable(1,1,1) some nonzero number
                         ! temporarily and restore it later.
  ENDIF

                           !>> FUNCTION BODY <<!

   ! Assign the reciprocal space potential to a temporary array.  Sorry, I
   ! thought about it a lot but couldn't figure out how to do it without one.
   IF(SIZE(jPotentialRecip,1) == SIZE(qTable,1) .AND. &
      SIZE(jPotentialRecip,2) == SIZE(qTable,2) .AND. &
      SIZE(jPotentialRecip,3) == SIZE(qTable,3)) THEN
      
      jPotentialRecip = (4._DP * pi) * jPotentialRecip / qTable**2

   ELSE
       jPotentialRecip = (4._DP * pi) * jPotentialRecip / &
                   WrapQ(qTable, SIZE(jPotentialRecip,1), SIZE(jPotentialRecip,2), &
                         SIZE(jPotentialRecip,3))**2
   END IF

   IF (k3Goff==0) THEN
     jPotentialRecip(1,1,1) = (0._DP,0._DP) ! Set q=0 term to 0
     qTable(1,1,1) = 0._DP   ! Done with qTable, set q(1,1,1) back to 0._DP
   END IF

   ! FFT to real space.
   potential = FFT(jPotentialRecip)

   ! Calculate energy
   IF (calcEnergy) THEN
     energy = SUM(potential*rhoReal)/2._DP
   ENDIF

!   CALL StopClock('JPotPlus');

END SUBROUTINE JPotentialPlus



FUNCTION JStress(cellVol, rhoQ, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the coulombic stress component 
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
!   02/05/2004  Function created.  (Vincent Ligneres)
!
!------------------------------------------------------------------------------
  USE Constants, ONLY: auToGPa
  USE PlaneWave, ONLY: qVectors
  USE CellInfo, ONLY: k1G, k2G, k3G
  USE CellInfo, ONLY: m3G ! equals global Z dimensions of recip FFT

  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN) :: cellVol ! The volume of the cell
  REAL(KIND=DP), INTENT(IN) :: energy  ! The electron-electron energy
  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoQ  ! Electron density in reciprocal space (spin-neutral)
  REAL(KIND=DP), DIMENSION(3,3) :: JStress              ! The final answer

  !>> INTERNAL VARIABLES <<!
  INTEGER :: a,b ! The two indices that define the stress component
  INTEGER :: ix, i2, i3    ! counter for parsing q over the recip. sphere
  REAL(KIND=DP) :: temp    ! Temporary variable for efficiency

                           !>> INITIALIZATION <<!

  ! We copy the dimensions of our grid from qMask, to ensure consistancy.
  JStress = 0._DP
                           !>> FUNCTION BODY <<!
  DO i3=1, k3G 
    DO i2=1, k2G 
      DO ix=1, k1G
        IF(qMask(ix,i2,i3))THEN
          ! Now compute the product of the components of interest.
          temp = CONJG(rhoQ(ix,i2,i3))*rhoQ(ix,i2,i3)/qTable(ix,i2,i3)**4
          DO a=1,3
            DO b=a,3          
                JStress(a,b) = JStress(a,b) + temp * qVectors(ix,i2,i3,a) * qVectors(ix,i2,i3,b) 
            END DO
          END DO
        END IF
      END DO
    END DO
  END DO

  JStress = JStress * 8._DP * pi

  DO a = 1,3
#ifdef __USE_PARALLEL
    ! m3G equals the third dimension of reciprocal space (total)
    JStress(a,a) = JStress(a,a) - energy / cellVol &
                                   *REAL(k3G,kind=DP)/REAL(m3G,kind=DP)
#else
    JStress(a,a) = JStress(a,a) - energy / cellVol
#endif
  END DO

!  WRITE(*,*) " stress for Hartree (Gpa): "
!  WRITE(*,'(3F12.6)') JStress(1,:)*auToGPa
!  WRITE(*,'(3F12.6)') JStress(2,:)*auToGPa
!  WRITE(*,'(3F12.6)') JStress(3,:)*auToGPa

END FUNCTION JStress

END MODULE Hartree 

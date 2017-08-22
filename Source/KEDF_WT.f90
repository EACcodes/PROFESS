MODULE KEDF_WT
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_WT
!     |_SUBROUTINE WTPotentialPlus 
!     |_SUBROUTINE WTPotPlusSB 
!     |_SUBROUTINE WTPotPlus 
!     |_FUNCTION WTEnergy
!     |_FUNCTION WTPotential
!     |_FUNCTION WTStress 
!   
! DESCRIPTION:
! Wang-Teter KEDF introduces the non-local term in KEDF,
! A very good KEDF for light metals.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
! REFERENCES: 
! Phys. Rev. B. 45, 13196 (1992) Lin-Wang Wang and Michael P. Teter
! "Kinetic-energy functional of the electron density"
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014 MODULE CREATED BY MOHAN CHEN, BASED ON PROFESS 2.0.
!------------------------------------------------------------------------------

  USE Constants,    ONLY : DP
  USE OutputFiles,  ONLY : outputUnit
  USE KEDF_WGC

  IMPLICIT NONE

CONTAINS

SUBROUTINE WTPotentialPlus(rho, potential, calcEnergy, energy, vacCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   An interface to decide the calculation of Wang-Teter kinetic potential
! and the energy.
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

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho        ! Electron density in real space (spin-independent)
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential ! The Wang-Teter contribution to the kinetic energy.
  LOGICAL, INTENT(IN) :: calcEnergy
  LOGICAL, INTENT(IN) :: vacCutoff
  REAL(KIND=DP), INTENT(OUT) :: energy

  CALL WTPotPlus(rho, potential, calcEnergy, energy, vacCutoff)

  RETURN

END SUBROUTINE WTPotentialPlus


SUBROUTINE WTPotPlus(rho, potential, calcEnergy, energy, vacCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine calculates the Wang-Teter kinetic potential based on the real-
!   space electron density. This is a generalized version of the WT functional
!   that admits variable exponents alpha and beta on rho on each side of the
!   kernel. The functional derivative for this is NOT in the paper, it is the
!   sum of the derivatives with respect to rho(r) and rho(r').
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Wang Y.A. Govind N. and Carter E.A., Phys. Rev. B 60(24), 1999, 16350
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/14/2004  File created.  (Vincent Ligneres)
!   08/15/2006  Changed to get rid of rho in denominator (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE


                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho        ! Electron density in real space (spin-independent)
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential ! The Wang-Teter contribution to the kinetic energy.
  LOGICAL, INTENT(IN) :: calcEnergy
  LOGICAL, INTENT(IN) :: vacCutoff
  REAL(KIND=DP), INTENT(OUT) :: energy

                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: cutf ! Cutoff function - used only if vacCutoff

  !-------------------------------- TEST -------------------------------------------!
  ! REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: pot1 
  ! REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: pot2 
  ! INTEGER :: i,j,k
  !-------------------------------- TEST -------------------------------------------!

                       !>> INITIALIZATION <<!

  IF( ALLOCATED(keKernel) .EQV. .FALSE. ) THEN
    WRITE(message, *) "keKernel is not allocated."
    CALL Error(6, message)
  ENDIF

  energy = 0._DP

  IF (vacCutoff .EQV. .TRUE.) THEN
    cutf = CutoffFunc(rho)
  END IF

                       !>> FUNCTION BODY <<!

  ! Kernel is always SIZE 1 in the 4th dimension, we just do this so that it
  ! will be compatible with the WGC kernel

  ! Temporary storage
  IF (vacCutoff) THEN
    potential = FFT(FFT(cutf*rho**beta) * keKernel(:,:,:,1))
  ELSE
    potential = FFT(FFT(rho**beta) * keKernel(:,:,:,1))
  END IF

  ! Calculate energy
  IF (calcEnergy) energy = cTF * SUM(rho**alpha * potential)

  ! Now calculate actual potential
  IF (vacCutoff) THEN
    potential = cTF * &
          (alpha * rho**(alpha-1._DP) * potential +&
           (beta * cutf * rho**(beta-1._DP) + CutoffFuncD(rho) * rho**beta) &
             * FFT(FFT(rho**alpha) * keKernel(:,:,:,1)))
  ELSE
    potential = cTF * &
          (alpha * rho**(alpha-1._DP) * potential +&
           beta * rho**(beta-1._DP) * FFT(FFT(rho**alpha) * keKernel(:,:,:,1)))

    !-------------------------------- TEST -------------------------------------------!
!    WRITE(outputUnit,*) "Output Wang-Teter Potential on i=9, j=9"
!    pot1 = cTF * (alpha * rho**(alpha-1._DP) * potential )
!    pot2 = cTF * (beta * rho**(beta-1._DP) * FFT(FFT(rho**alpha) * keKernel(:,:,:,1)))
!    DO i=1, SIZE(rho,1)
!      DO j=1, SIZE(rho,2)
!        DO k=1, SIZE(rho,3)
!          IF(i==9 .AND. j==9) THEN
!            WRITE(outputUnit,'(1I5,4ES20.2)') &
!            k,pot1(i,j,k),pot2(i,j,k),potential(i,j,k), rho(i,j,k)
!          ENDIF
!        ENDDO
!      ENDDO
!    ENDDO
    !-------------------------------- TEST -------------------------------------------!
  END IF
END SUBROUTINE WTPotPlus


FUNCTION WTEnergy(rhoR_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine calculates the Wang-Teter kinetic energy based on the real-
!   space electron density. This is a generalized version of the WT functional
!   that admits variable exponents alpha and beta on rho on each side of the
!   kernel.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Wang Y.A. Govind N. and Carter E.A., Phys. Rev. B 60(24), 1999, 16350
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/04/2003  File created.  (Vincent Ligneres)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP) :: WTEnergy ! The Wang-Teter contribution to the kinetic energy.
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoR_SI ! Electron density in real space (spin-independent)

  !>> INTERNAL VARIABLES <<!
  !>> INITIALIZATION <<!
  !>> FUNCTION BODY <<!

  ! Kernel is always SIZE 1 in the 4th dimension, we just do this so that it
  ! will be compatible with the WGC kernel
  !  WTEnergy = 2.0_DP * cTF * cellVol * &
  !    SUM( FFT(rhoR_SI**beta) * CONJG(FFT(rhoR_SI**alpha)) * keKernel(:,:,:,1), &
  !        qMask)

  ! This is an alternate method of computing the WT energy.
  WTEnergy = cTF * SUM(rhoR_SI**alpha * FFT(FFT(rhoR_SI**beta) * keKernel(:,:,:,1)))

END FUNCTION WTEnergy


FUNCTION WTPotential(rho)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine calculates the Wang-Teter kinetic potential based on the real-
!   space electron density. This is a generalized version of the WT functional
!   that admits variable exponents alpha and beta on rho on each side of the
!   kernel. The functional derivative for this is NOT in the paper, it is the
!   sum of the derivatives with respect to rho(r) and rho(r').
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Wang Y.A. Govind N. and Carter E.A., Phys. Rev. B 60(24), 1999, 16350
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/14/2004  File created.  (Vincent Ligneres)
!   08/15/2006  Changed to get rid of rho in denominator (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    rho                   ! Electron density in real space (spin-independent)

  REAL(KIND=DP), DIMENSION(SIZE(rho,1), SIZE(rho,2), SIZE(rho,3)) :: &
    WTPotential           ! The Wang-Teter contribution to the kinetic energy.

                     !>> INTERNAL VARIABLES <<!
                       !>> INITIALIZATION <<!
                       !>> FUNCTION BODY <<!
  ! Kernel is always SIZE 1 in the 4th dimension, we just do this so that it
  ! will be compatible with the WGC kernel
  WTPotential = cTF * &
        (alpha * rho**(alpha-1._DP) * FFT(FFT(rho**beta) * keKernel(:,:,:,1)) +&
         beta * rho**(beta-1._DP) * FFT(FFT(rho**alpha) * keKernel(:,:,:,1)))

END FUNCTION WTPotential


FUNCTION WTStress(rhoReal, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the Wang-Teter portion of the stress 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE MathFunctions, ONLY: GPrime   ! Derivative of the Lindhard G function.
  USE CellInfo, ONLY: k1G, k2G, k3G, k3Goff
  USE CellInfo, ONLY: n3G, m3G
  USE CellInfo, ONLY: cell
  USE Constants, ONLY: PI
  USE PlaneWave, ONLY: qVectors
  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(3,3) :: WTStress
  ! the answer
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoReal  
  ! Spin-independant density in real space.
  !
  REAL(KIND=DP), INTENT(IN) :: energy 
  ! The WT energy from the energy table.
  !

                         !>> INTERNAL VARIABLES <<!
  INTEGER :: a, b   
  ! counters
  !
  REAL(KIND=DP) :: eta 
  ! Dimensionless momentum (norm of q over twice the Fermi vector)
  !
  REAL(KIND=DP) :: diff
  ! Temporary storage for the Lindhard differential
  !
  REAL(KIND=DP) :: coef
  ! Extra diagonal component of the stress if rho0 fluctuates.
  !
  REAL(KIND=DP) :: mult
  ! Diagonal element multiplier based on wether rho0 is held.
  !
  REAL(KIND=DP) :: kF
  ! Fermi wave-vector
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FrhoA
  ! FFT of rho^alpha.
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FrhoB
  ! FFT of rho^beta
  !
  INTEGER :: ix, i2, i3
  ! counter for parsing q over the recip sphere
  !
  REAL(KIND=DP), DIMENSION(3) :: qPoint 
  ! cartesian coordinates of the q vector
  !

                           !>> INITIALIZATION <<!

  WTStress = 0._DP

  kF = (3 * pi**2 * rho0)**(1._DP/3._DP)
  ! technically we should FFT rho^a - rho0^a but this only shifts rho^a in
  ! real space and therefore only changes the q=0 term. Since we don't include
  ! the q=0 term in the reciprocal space sum, it doesn't matter.

  FrhoA = FFT(rhoReal**alpha)

  FrhoB = CONJG(FFT(rhoReal**beta))

  ! The kernel is a function of rho0 and eta. eta is a function of rho0. If
  ! rho0 is N/V then that volume dependance translates into a stress component.
  ! Otherwise rho0 is simply a constant.
  IF (hold0) THEN
    coef = 0._DP
    mult = 1._DP - alpha - beta
  ELSE
    coef = -1._DP / 3._DP
    mult = -2._DP / 3._DP
  END IF

                           !>> FUNCTION BODY <<!

  ! Loop through all q vectors in the half-sphere. 
  DO i3 = 1, k3G
    DO i2 = 1, k2G
      DO ix = 1, k1G 
        ! Go further only if this point matters (within the KE cutoff)
        IF (qMask(ix,i2,i3)) THEN
          ! Calculate the qPoint cartesian coordinates
          qPoint = qVectors(ix,i2,i3,:) 

          eta = qTable(ix,i2,i3) / (2._DP * kF)
          ! Then calculate the differential of the Lindhard function at this q.
          diff = - eta * GPrime(eta, mu)
          ! Overwrite diff to include the two rho terms
          diff = diff * FrhoA(ix,i2,i3) * FrhoB(ix,i2,i3)
          ! and now for the 6 stress terms
          DO a = 1, 3
            DO b = a, 3
              WTStress(a,b) = WTStress(a,b) + diff * &
                (qPoint(a)*qPoint(b) / qTable(ix,i2,i3)**2 + coef*(b-a-1)*(b-a-2)/2._DP)
            END DO ! b
          END DO ! a
        END IF ! in sphere
      END DO ! ix
    END DO ! i2
  END DO ! i3

  ! Now multiply by the prefactor
  WTStress = WTStress * pi**2 / (alpha * beta * kF * rho0**(alpha+beta-2))

  ! Finally, add the diagonal term
  DO a = 1, 3
    WTStress(a,a) = WTStress(a,a) + mult * energy / cell%vol & 
                    * REAL(n3G,kind=DP)/REAL(m3G,kind=DP)
  END DO

  RETURN

END FUNCTION WTStress

END MODULE KEDF_WT

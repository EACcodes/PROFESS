MODULE KEDF_Q 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_Q 
!     |_SUBROUTINE CalLHQ
!     |_FUNCTION LQEnergy (LQ Functional (Jeng Da))
!     |_FUNCTION LQPotential
!     |_FUNCTION HQEnergy (HQ Functional (Jeng Da)) 
!     |_FUNCTION HQPotential
!
! DESCRIPTION:
!   Calculates the kinetic forces, stress, and energies.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Orbital-Free Density Functional Theory: Kinetic Potentials and Ab Initio
!       Local Pseudopotentials, Jeng-Da Chai and John D. Weeks
!       Phys. Rev. B 75, 205122 (2007)
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014 Module created by Mohan.
!
!------------------------------------------------------------------------------

  USE PlaneWave, ONLY: qVectors
  USE PlaneWave, ONLY: qTable
  USE Fourier, ONLY: FFT
  USE CONSTANTS, ONLY: DP
  USE CONSTANTS, ONLY: PI
  USE CONSTANTS, ONLY: IMAG
  USE KEDF_TF, ONLY: lambda
  USE KEDF_TF, ONLY: cTF
  USE KEDF_VW, ONLY: mu
  USE SYS, ONLY: rho0

  IMPLICIT NONE

CONTAINS

SUBROUTINE CalLHQ(potential, rho, calcEnergy, optSqrt, eTF, eVW, eQ, LQflag)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE KEDF_TF, ONLY: TFPotential
  USE KEDF_TF, ONLY: TFEnergy
  USE KEDF_VW, ONLY: VWPotentialSqrtPlus
  USE CellInfo, ONLY: n1G, n2G, n3G

  IMPLICIT NONE

  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G), INTENT(OUT) :: potential
  ! potential we want for the answer
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G), INTENT(IN) :: rho
  ! rho for optimization
  !
  LOGICAL :: calcEnergy
  ! calculate the energies for KEDF or not
  !
  LOGICAL :: optSqrt
  ! optimize sqrt(rho) if optSqrt=.TRUE.
  !
  REAL(KIND=DP), INTENT(OUT) :: eTF
  ! KEDF energy for Thomas-Fermi term
  !
  REAL(KIND=DP), INTENT(OUT) :: eVW
  ! KEDF energy for VW term
  !
  REAL(KIND=DP), INTENT(OUT) :: eQ
  ! KEDF energy for KEDF Q term
  !
  LOGICAL :: LQflag
  ! If true, use LQ KEDF
  ! If false, use HQ KEDF
  !
  !! >> LOCAL VARIABLSE << !!
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tempPotential
  !
  INTEGER :: allocateStatus
  ! error flag 
  !

  
  !! >> INITIALIZATION << !!

  ALLOCATE(tempPotential(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*)'Error allocating tempPotential in KEDF_CAT. Leaving.'
    STOP
  END IF

  ! Compute the TF term. 
  potential = potential + TFPotential(rho)
  IF (calcEnergy) THEN
    eTF = TFEnergy(rho)
  ENDIF
 
  ! Lindhard-like term.
  IF(LQflag .EQV. .TRUE.) THEN
    potential = potential + LQPotential(rho)
  ELSE
    potential = potential + HQPotential(rho)
  ENDIF

  IF (calcEnergy) THEN
    IF(LQflag .EQV. .TRUE.) THEN
      eQ = LQEnergy(rho)
    ELSE
      eQ = HQEnergy(rho)
    ENDIF
  ENDIF

  ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
  IF (optSqrt) THEN
    potential = 2._DP * SQRT(rho) * potential
  ENDIF

  ! Compute vW term. 
  CALL VWPotentialSqrtPlus(SQRT(rho), tempPotential, calcEnergy, eVW)

  ! Chain rule: dE/d(rho) = dE/d(sqrt(rho)) / (2*sqrt(rho))
  IF (.NOT. optSqrt) THEN
    tempPotential = tempPotential/(2._DP*SQRT(rho))
  ENDIF

  potential = potential + tempPotential

  DEALLOCATE(tempPotential)

  RETURN

END SUBROUTINE CalLHQ


FUNCTION LQEnergy(rho) 
!------------------------------------------------------------------------------ 
! DESCRIPTION: 
! 
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG: 
!   06/22/2005  Function created.  (Vincent Ligneres)
!------------------------------------------------------------------------------
  USE CellInfo, ONLY: cell, m3G, n3Goff ! The cell vectors, ions, etc...
  USE MATHFUNCTIONS, ONLY : LindG                 ! The Lindhard function

  IMPLICIT NONE 

                      !>> EXTERNAL VARIABLES <<!
 
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    rho                   ! Electron density in real space (spin-independent)
 
  REAL(KIND=DP) :: &
    LQEnergy
 
                     !>> INTERNAL VARIABLES <<!
  COMPLEX(KIND=DP), DIMENSION(SIZE(qTable,1), SIZE(qTable,2), SIZE(qTable,3))::&
    workArray             ! recip. space array to store intermediate results.
 
  REAL(KIND=DP), DIMENSION(SIZE(rho,1), SIZE(rho,2), SIZE(rho,3)) :: &
    rArray, &             ! Projection of the r vectors along one component 
                          ! (x, y, z) 
    gradRho               ! The Wang-Teter contribution to the kinetic energy. 

  REAL(KIND=DP) :: &
    tkf                   ! 2 * the rho-independent fermi wave-vector.
 
  INTEGER :: ix, iy, iz, i2, i3, ia            ! counters
 
                       !>> INITIALIZATION <<!
  tkf = 2._DP * (3._DP * pi**2 * rho0)**(1._DP/3._DP)


                       !>> FUNCTION BODY <<! 
  workArray = FFT(rho)
  gradRho=0._DP
  DO ia=1, 3
    DO iz=1, SIZE(rho,3)
      DO iy=1, SIZE(rho,2)
        DO ix=1, SIZE(rho,1)
          rArray(ix,iy,iz) = (ix-1000) *cell%cellReal(1,ia) &
                              /REAL(SIZE(rho,1),KIND=DP) &
                           + (iy-1000) *cell%cellReal(2,ia) &
                              /REAL(SIZE(rho,2),KIND=DP) &
                           + (iz+n3Goff-1000) *cell%cellReal(3,ia) &
                              /REAL(m3G,KIND=DP) 
        END DO ! ix
      END DO ! iy
    END DO ! iz
    gradRho = gradRho + FFT(imag * qVectors(:,:,:,ia) * workArray) * rArray
  END DO ! ia
  gradRho = 3._DP * rho - gradRho
  workArray=FFT(SQRT(rho))
  DO i3=1, SIZE(qTable,3) 
    DO i2=1, SIZE(qTable,2) 
      DO ix=1, SIZE(qTable,1) 
        workArray(ix,i2,i3)=LindG(qTable(ix,i2,i3)/tkf,lambda,mu) * &
                            workArray(ix,i2,i3)
      END DO ! ix 
    END DO ! i2
  END DO !i3
  gradRho = FFT(workArray) * gradRho
  gradRho = gradRho / SQRT(rho)
  LQEnergy = SUM(gradRho) 
  LQEnergy = LQEnergy * tkf**2 / 12._DP

  RETURN
 
END FUNCTION LQEnergy


FUNCTION LQPotential(rho) 
!------------------------------------------------------------------------------
! DESCRIPTION: 
! 
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG: 
!   06/22/2005  Function created.  (Vincent Ligneres) 
!------------------------------------------------------------------------------ 
  USE MATHFUNCTIONS, ONLY : LindG
  ! The Lindhard function
  !

  IMPLICIT NONE 
 
                      !>> EXTERNAL VARIABLES <<! 
 
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho
  ! Electron density in real space (spin-independent) 
  !
  REAL(KIND=DP), DIMENSION(SIZE(rho,1), SIZE(rho,2), SIZE(rho,3)) :: LQPotential
  ! The Wang-Teter contribution to the kinetic energy. 
  !
 
                     !>> INTERNAL VARIABLES <<! 

  COMPLEX(KIND=DP), DIMENSION(SIZE(qTable,1), SIZE(qTable,2), SIZE(qTable,3)):: workArray
  ! recip. space array to store intermediate results.
  !
  REAL(KIND=DP) :: tkf
  ! 2 * the rho-independent fermi wave-vector.
  !
  INTEGER :: ix, i2, i3
  ! counters
  !

                       !>> INITIALIZATION <<! 

  tkf = 2._DP * (3._DP * pi**2 * rho0)**(1._DP/3._DP)
  workArray = FFT(SQRT(rho))
                       !>> FUNCTION BODY <<! 
  DO i3= 1, SIZE(qTable,3)
    DO i2=1, SIZE(qTable,2)
      DO ix=1, SIZE(qTable,1)
        workArray(ix,i2,i3)=LindG(qTable(ix,i2,i3)/tkf,lambda,mu) * &
                            workArray(ix,i2,i3)
      END DO !ix
    END DO ! i2
  END DO ! i3
  LQPotential = FFT(workArray)/SQRT(rho)
  LQPotential = LQPotential * tkf**2 / 6._DP

  RETURN
 
END FUNCTION LQPotential



FUNCTION HQEnergy(rho) 
!------------------------------------------------------------------------------ 
! DESCRIPTION: 
! 
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! 
!------------------------------------------------------------------------------ 
! REVISION LOG:
!   06/22/2005  Function created.  (Jeng-Da Chai)
!------------------------------------------------------------------------------ 
  USE CellInfo, ONLY :  cell, n3Goff     ! The cell vectors, ions, etc... 
  USE MATHFUNCTIONS, ONLY :  LindG
  ! The Lindhard function 
 
  IMPLICIT NONE 
 
                      !>> EXTERNAL VARIABLES <<! 
 
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho
  ! Electron density in real space (spin-independent) 
  !
  REAL(KIND=DP) :: HQEnergy 
 
                     !>> INTERNAL VARIABLES <<! 

  COMPLEX(KIND=DP), DIMENSION(SIZE(qTable,1), SIZE(qTable,2), SIZE(qTable,3)):: workArray
  ! recip. space array to store intermediate results. 
  !
  REAL(KIND=DP), DIMENSION(SIZE(rho,1), SIZE(rho,2), SIZE(rho,3)) :: rArray
  ! Projection of the r vectors along one component (x, y, z) 
  !
  REAL(KIND=DP), DIMENSION(SIZE(rho,1), SIZE(rho,2), SIZE(rho,3)) :: gradRho
  ! The Wang-Teter contribution to the kinetic energy. 
  !
  REAL(KIND=DP) :: tkf
  ! 2 * the rho-independent fermi wave-vector. 
  !
  INTEGER :: ix, iy, iz, i2, i3, ia            ! counters 
 
                       !>> INITIALIZATION <<! 
  tkf = 2._DP * (3._DP * pi**2 * rho0)**(1._DP/3._DP) 

                       !>> FUNCTION BODY <<! 
  workArray = FFT(rho) 
  gradRho=0._DP 
  DO ia=1, 3 
    DO iz=1, SIZE(rho,3)
      DO iy=1, SIZE(rho,2)
        DO ix=1, SIZE(rho,1) 
          rArray(ix,iy,iz) = ix * cell%cellReal(1,ia) + &
                             iy * cell%cellReal(2,ia) + &
                             (iz+n3Goff) * cell%cellReal(3,ia)
        END DO ! ix 
      END DO ! iy
    END DO ! iz 
    gradRho = gradRho + FFT(imag * qVectors(:,:,:,ia) * workArray) * rArray 
  END DO ! ia
  gradRho = 3._DP * rho - gradRho 
  workArray=FFT(cTF*(5._DP/3._DP)*rho**(2._DP/3._DP)) 
  DO i3= 1, SIZE(qTable,3) 
    DO i2=1, SIZE(qTable,2) 
      DO ix=1, SIZE(qTable,1)
        workArray(ix,i2,i3)=LindG(qTable(ix,i2,i3)/tkf,lambda,mu) * &
                            workArray(ix,i2,i3)
      END DO ! ix 
    END DO ! i2
  END DO !i3
  gradRho = FFT(workArray) * gradRho
  HQEnergy = SUM(gradRho) 
  HQEnergy = HQEnergy / 2._DP
 
END FUNCTION HQEnergy 
 

FUNCTION HQPotential(rho) 
!------------------------------------------------------------------------------
! DESCRIPTION: 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! 
!------------------------------------------------------------------------------ 
! REVISION LOG: 
!   06/22/2005  Function created.  (Vincent Ligneres)
!------------------------------------------------------------------------------
  USE MATHFUNCTIONS, ONLY : LindG
  ! The Lindhard function 
  !
 
  IMPLICIT NONE 
 
                      !>> EXTERNAL VARIABLES <<! 
 
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho
  ! Electron density in real space (spin-independent) 
  !
  REAL(KIND=DP), DIMENSION(SIZE(rho,1), SIZE(rho,2), SIZE(rho,3)) :: HQPotential
  ! The Wang-Teter contribution to the kinetic energy. 
  !
 
                     !>> INTERNAL VARIABLES <<! 

  COMPLEX(KIND=DP), DIMENSION(SIZE(qTable,1), SIZE(qTable,2), SIZE(qTable,3)) :: workArray
  ! recip. space array to store intermediate results. 
  !
  REAL(KIND=DP) :: tkf
  ! 2 * the rho-independent fermi wave-vector. 
  !
  INTEGER :: ix, i2, i3
  ! counters 
  !
 
                       !>> INITIALIZATION <<! 

  tkf = 2._DP * (3._DP * pi**2 * rho0)**(1._DP/3._DP) 
  workArray = FFT(cTF*(5._DP/3._DP)*rho**(2._DP/3._DP)) 
                       !>> FUNCTION BODY <<! 
  DO i3= 1, SIZE(qTable,3) 
    DO i2=1, SIZE(qTable,2) 
      DO ix=1, SIZE(qTable,1) 
        workArray(ix,i2,i3)=LindG(qTable(ix,i2,i3)/tkf,lambda,mu) * &
                            workArray(ix,i2,i3) 
      END DO !ix
    END DO ! i2
  END DO ! i3
  HQPotential = FFT(workArray)

  RETURN
 
END FUNCTION HQPotential 

END MODULE KEDF_Q

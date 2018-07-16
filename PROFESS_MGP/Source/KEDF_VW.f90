MODULE KEDF_VW 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_VW 
!     |_SUBROUTINE CalVW
!     |_FUNCTION VWEnergy
!     |_FUNCTION VWPotential
!     |_FUNCTION VWPotentialSqrt
!     |_FUNCTION VWPotentialSqrtPlus
!     |_FUNCTION VWStress
!
! DESCRIPTION:
!   Calculates the kinetic forces, stress, and energies of
!   von Weiszacker KEDF
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/15/2003  File Created by consolidating existing energy files (GSH) 
!   12/18/2003  Changed SUBROUTINE to FUNCTIONS
!   01/07/2004  Changed qTable and qMask and most instances of rhoR to 
!               spin-independent forms.
!   02/05/2004  Added nine stress functions. (VLL)
!   02/25/2004  Reorganized into different files (GSH)
!   03/23/2005  Added energy and potential for hybrid TF (VLL)
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!
  USE Constants, ONLY: DP
  USE Constants, ONLY: IMAG
  USE Fourier, ONLY: FFT
  USE PlaneWave, ONLY: qTable
  USE PlaneWave, ONLY: qVectors
  USE OutputFiles, ONLY: outputUnit
  USE OutputFiles, ONLY: errorUnit

  IMPLICIT NONE

  REAL(KIND=DP) :: mu = -100.0_DP     
  ! Multiplier for the VW contribution.
  !

CONTAINS


SUBROUTINE CalVW(potential, rho, calcEnergy, optSqrt, eVW)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the VW potential and energy.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell, numSpin
  USE CellInfo, ONLY: n1G, n2G, n3G

  IMPLICIT NONE
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(OUT) :: potential
  ! TF potential, last dimension is spin
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  ! charge density on real space
  ! 
  LOGICAL, INTENT(IN) :: calcEnergy
  ! calculate TF energy or not
  ! 
  LOGICAL, INTENT(IN) :: optSqrt
  ! optimize sqrt(rho) in the whole algorithm or not
  !
  REAL(KIND=DP), INTENT(OUT) :: eVW
  ! VW term energy
  !
  !! >> LOCAL VARIABLES << !!
  !
  INTEGER :: isp
  ! spin index
  !
  REAL(KIND=DP) :: etmp
  ! tmp energy for VW energy
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tempPotential
  !
  INTEGER :: allocateStatus
  !

  !! >> INITIALIZATION << !!

  ALLOCATE(tempPotential(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*)'Error allocating tempPotential in KEDF_VW. Leaving.'
    STOP
  END IF

  !! >> FUNCTION << !!

  SELECT CASE(numSpin)
    CASE(1)

      ! Von Weiszacker (vW)
      CALL VWPotentialSqrtPlus(SQRT(rho(:,:,:,1)), tempPotential, calcEnergy, eVW)

      ! Chain rule: dE/d(rho) = dE/d(sqrt(rho)) / (2*sqrt(rho))
      IF (.NOT. optSqrt) THEN
        tempPotential = tempPotential/(2._DP*SQRT(rho(:,:,:,1)))
      ENDIF
             
      potential = potential + SPREAD(tempPotential, 4, numSpin)

    CASE(2)
      
      eVW = 0.0_DP

      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      DO isp = 1,2
        ! Von Weiszacker (vW)
!       potential(:,:,:,isp) = potential(:,:,:,isp) +  VWPotential(2.d0*rho(:,:,:,isp))
!       IF (calcEnergy) locETable(8) = locETable(8) +  0.5d0*VWEnergy(2.d0*rho(:,:,:,isp))

        CALL VWPotentialSqrtPlus(SQRT(2.d0*rho(:,:,:,isp)), tempPotential, calcEnergy, etmp)

        tempPotential = tempPotential/SQRT(2.d0)

        IF (.NOT. optSqrt) THEN
          tempPotential = tempPotential/(2._DP*SQRT(rho(:,:,:,isp)))
        ENDIF

        potential(:,:,:,isp) = potential(:,:,:,isp) + tempPotential
                
        IF (calcEnergy) THEN
          eVW = eVW +  etmp*0.5
        ENDIF

      ENDDO
  END SELECT

  RETURN

END SUBROUTINE CalVW


FUNCTION VWEnergy(sqrtRhoR_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the von Weizsacker kinetic energy based on real-
!   space electron density and the table of norms of q-vectors in recip. space.
!   This function uses FFT's.
!
!   Here is an alternate definition of the VW energy:
!
!   VWEnergy = mu/8._DP * SUM((FFT(imag*qVectors(:,:,:,1)*FFT(rhoR_SI))**2 + &
!              FFT(imag*qVectors(:,:,:,2)*FFT(rhoR_SI))**2 + &
!              FFT(imag*qVectors(:,:,:,3)*FFT(rhoR_SI))**2)/rhoR_SI) * &
!              cellVol / SIZE(rhoR_SI)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/24/2003  File created.  (VLL)
!   07/31/2006  Changed to take square root of density as input (GSH)
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP) :: VWEnergy
  ! The VW energy.
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: sqrtRhoR_SI
  ! Sqrt of Electron density in real space 
  ! (spin-independent)

  VWEnergy = 0.5_DP * mu * SUM(sqrtRhoR_SI*FFT(FFT(sqrtRhoR_SI)*qTable**2))

  RETURN

END FUNCTION VWEnergy


FUNCTION VWPotential(rhoR_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   1/12/2003 Function created (GSH)
!------------------------------------------------------------------------------

  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoR_SI
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(SIZE(rhoR_SI, 1),SIZE(rhoR_SI, 2), SIZE(rhoR_SI, 3)) :: VWPotential
  ! The VW potential.
  !
  !>> LOCAL VARIABLES <<!
  !
  ! Store SQRT(rhoR_SI) first to avoid having to do the calculation twice.
  VWPotential = SQRT(rhoR_SI)
  VWPotential = mu * 0.5_DP * FFT( FFT(VWPotential) * qTable**2 ) / VWPotential

  RETURN

END FUNCTION VWPotential


FUNCTION VWPotentialSqrt(rhoR_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the Von-Weizsacker potential in real-space.  Gotten by taking the 
!   functional derivative of the formula used for the VWEnergy.
!
!   I am puzzled.  By taking the functional derivitive of the expression used
!   for the VW Energy, I get the following formula:
!
!    VWPotential = 0.25_DP * &
!                  (SQRT(rhoR_SI) * FFT( FFT(1/SQRT(rhoR_SI)) * qTable**2) + &
!                   1/SQRT(rhoR_SI) * FFT( FFT(SQRT(rhoR_SI)) * qTable**2))
!
!   However, Alex's code doesn't use this formula... instead, it uses a simpler
!   formula (the one that is implemented here).  It seems to give the same
!   answer, so I'll use it.   But I'd be more comfortable if someone could show
!   that they are the same.
!
!   UPDATE:  Vincent and I figured out why this formula works.  Apperantly,
!            the two forms of the VW energy (sqrt version and 1/p version)
!            are not equal ... they are equal in their kinetic energy
!            SUM, but not equal in the kinetic energy density...they are off
!            by a factor that becomes 0 if you sum over the entire grid.
!            Anyway, the derivitives of the two forms are also different.
!            The derivitive we're using here is NOT consistant with the
!            energy expression here (we're using the derivitive of the
!            1/p form) but that's okay as long as we're using periodic
!            cells. Otherwise, we need to beware (as in embedding). Actually,
!            in that case, the potential expression here is probably correct,
!            its the energy we'd need to change.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Figure out why this formula works?
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   1/12/2003 Function created (GSH)
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoR_SI 
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(SIZE(rhoR_SI, 1),SIZE(rhoR_SI, 2), SIZE(rhoR_SI, 3)) :: VWPotentialSqrt           
  ! The VW potential.
  !
  VWPotentialSqrt = mu * FFT( FFT(SQRT(rhoR_SI)) * qTable**2)

  RETURN
  
END FUNCTION VWPotentialSqrt


SUBROUTINE VWPotentialSqrtPlus(sqrtRho_SI, potentialSqrt, calcEnergy, energy, weightingfunc, dV, IfCalForce)    
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the Von-Weizsacker potential in real-space and energy.  Gotten by
!   taking the functional derivative of the formula used for the VWEnergy.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/12/2003 Function created (GSH)
!   09/08/2008 Combined energy and potential subroutine (LH)
!------------------------------------------------------------------------------

  USE Fourier_NEW, ONLY: FFT_NEW, FFT_STD_STATE
  USE CellInfo, ONLY: k1G, k2G, k3G
  USE MathFunctions, ONLY : GlobalVal3d

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: sqrtRho_SI      
  ! Sqrt electron density in real space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potentialSqrt           
  ! The VW potential.
  !
  LOGICAL, INTENT(IN) :: calcEnergy
  !
  REAL(KIND=DP), INTENT(OUT) :: energy
  !
  REAL(KIND=DP), INTENT(IN) , OPTIONAL:: dV
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: weightingfunc
  !
  LOGICAL, INTENT(IN), OPTIONAL :: IfCalForce
  !

                          !>> INTERNAL VARIABLES <<!
  !                        
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G)::trafo
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: trafoX
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: Rho, DeltaRho2, Delta2Rho, deltaRhox, deltaRhoy, deltaRhoz, tmpX
  !
  INTEGER::status,lX,uX,lY,uY,lZ,uZ
  !
                           !>> INITIALIZATION <<!
  energy = 0._DP
  CALL StartClock('VWPotSqrPlus')
                           !>> FUNCTION BODY <<!

!  WRITE(*,*) "ifCalForce=", ifCalForce

  if(.NOT.present(ifCalForce) .OR. .NOT.ifCalForce) THEN

    ! could not avoid this energy difference due to differnt FFTW we use.
    CALL FFT_NEW(FFT_STD_STATE,sqrtRho_SI,trafo)
    trafo = trafo * qTable * qTable
    CALL FFT_NEW(FFT_STD_STATE,trafo,potentialSqrt)
  
    IF (calcEnergy) THEN
      energy = 0.5_DP * SUM(sqrtRho_SI*potentialSqrt)
    ENDIF

    energy = energy * mu
    potentialSqrt = potentialSqrt * mu

  ELSE
                   
  ENDIF

  CALL StopClock('VWPotSqrPlus')
  
  ! JMD: Fortran standard defines automatic deallocation at the end of the variables scope
  
END SUBROUTINE VWPotentialSqrtPlus


FUNCTION VWStress(sqrtRhoQ, sizeRealRho)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the von-weisaker stress terms.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! qA and qB should be deleted
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/05/2004  Function created.  (Vincent Ligneres)
!   09/03/2008  Changed to derivative of d(sqrt(rho))/dr to remove factor of
!               1/rho  (Linda Hung)
!------------------------------------------------------------------------------
  USE MATHFUNCTIONS, ONLY : Vecmul
  USE CellInfo, ONLY: k1G, k2G, k3G

  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3,3) :: VWStress
  ! The answer
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: sqrtRhoQ      
  ! Spin-independant density in recip. space.
  !
  INTEGER, INTENT(IN) :: sizeRealRho   
  ! Total number of gridpoints in 3D real array for density
  !

                         !>> INTERNAL VARIABLES <<!
  INTEGER :: a, b   ! counters 
  !
  REAL(KIND=DP), DIMENSION(k1G,k2G,k3G) :: qA
  ! The a component of every half-sphere q.
  !
  REAL(KIND=DP), DIMENSION(k1G,k2G,k3G) :: qB
  ! The b component
  !
  INTEGER :: ix, i2, i3    
  ! counter for parsing q over the recip sphere
  !

                           !>> INITIALIZATION <<!
  VWStress = 0._DP
 
                           !>> FUNCTION BODY <<!

  ! Loop through all q vectors in the half-sphere. 
  DO a = 1, 3
  DO b = a, 3

  DO i3 = 1, k3G
    DO i2 = 1, k2G
      DO ix = 1, k1G 
        qA(ix,i2,i3) = qVectors(ix,i2,i3,a)
        qB(ix,i2,i3) = qVectors(ix,i2,i3,b)
      END DO
    END DO
  END DO

  VWStress(a,b) = -mu  &
                 * SUM(FFT(imag * qA *sqrtRhoQ) * FFT(imag * qB * sqrtRhoQ))&
                 / REAL(sizeRealRho, KIND=DP)

! This old implementation is unstable in vacuum
!#ifdef __USE_PARALLEL
!  VWStress(a,b) = mu * (-0.25_DP) &
!                 * SUM(FFT(imag * qA *rhoQ) * FFT(imag * qB * rhoQ) / rhoReal)&
!                 / REAL(SIZE(rhoReal,1)*totY*recipDim2, KIND=DP)
!#else
!  VWStress(a,b) = mu * (-0.25_DP) &
!                 * SUM(FFT(imag * qA *rhoQ) * FFT(imag * qB * rhoQ) / rhoReal)&
!                 / REAL(SIZE(rhoReal), KIND=DP)
!#endif

  END DO
  END DO

  RETURN

END FUNCTION VWStress


END MODULE KEDF_VW 

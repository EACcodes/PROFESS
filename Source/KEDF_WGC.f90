MODULE KEDF_WGC
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_WGC
!     |_FUNCTION WGCPotentialPlus 
!     |_FUNCTION CutoffFunc
!     |_FUNCTION CutoffFuncD
!     |_FUNCTION WGCStress
!
! DESCRIPTION:
! The non-local KEDF kernel contains a double-density-dependent kernel.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
! REFERENCES: 
! Phys. Rev. B. 60, 16350 (1999) Yan Alexander Wang, Niranjan Govind
! and Emily A. Carter
! "Orbital-free kinetic-energy density functionals with a density-dependent
!  kernel"
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014 MODULE CREATED BY MOHAN CHEN, BASED ON PROFESS 2.0
!------------------------------------------------------------------------------

  USE Constants, ONLY: DP
  USE OutputFiles, ONLY: outputUnit
  USE KEDF_TF, ONLY: cTF
  USE KEDF_WTkernel, ONLY: alpha, beta
  USE SetupFFT
  USE Fourier_NEW
  USE Output, ONLY: WrtOut
  USE Fourier, ONLY: FFT
  USE PlaneWave, ONLY: qTable
  USE PlaneWave, ONLY: qMask
  USE CellInfo, ONLY: n1G, n2G, n3G
  USE CellInfo, ONLY: k1G, k2G, k3G
  USE KEDF_WGCkernel, ONLY: WGCT, firstOrderWGC, gamma 
  USE KEDF_WTkernel, ONLY: keKernel 
  USE Sys, ONLY: rho0, rhoS, hold0, holdS

  IMPLICIT NONE

  REAL(KIND=DP) :: rhoV = 1.E-6    ! for cutoff function used in WGCvac KEDF
  REAL(KIND=DP) :: rhoStep = 1.E-6 ! for cutoff function used in WGCvac KEDF

CONTAINS


SUBROUTINE WGCPotentialPlus(rho, potential, calcEnergy, energy, vacCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine calculates the Wang-Govind-Cearter potential based on the
!   real-space electron density. This is a general version of WGC functional
!   that can take practically any alpha, beta, gamma combination without prior
!   requirement. It is based on Stuart's code, and my improvements upon it.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/03/2003  File adapted from WGCPot_cheap.  (Vincent Ligneres)
!   10/22/2013  FFTW3 new interface, removal of temporaries (JMD)
!------------------------------------------------------------------------------

  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho ! Electron density in real space (spin-independent)
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential ! The WGC contribution to the potential.

  LOGICAL, INTENT(IN) :: calcEnergy   ! Whether to calculate energy
  LOGICAL, INTENT(IN) :: vacCutoff    ! Whether to use the cutoff function to regulate
                                      ! behavior of functional in vacuum

  REAL(KIND=DP), INTENT(OUT) :: energy ! WGC term energy

                      !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: &
    cutf, &     ! Cutoff function - used only if vacCutoff
    cutfp, &    ! Derivative of cutoff function - used only if vacCutoff
    rhoX, &     ! Storage for rho^alpha and rho^beta
    theta, &    ! Rho - rho*, expansion term - used only with WGCT=2
    workRA, &   ! Array for intermediate computations
    workRA2     ! Array for intermediate computations - only when WGCT=2

  ! ------------------------------- TEST -------------------------------------
  !REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: &
  !  pot1, &     ! (alpha+1) * rho1^{alpha} * <W2|rho2^{beta}>
  !  pot2, &     ! alpha * rho1^{alpha-1} * <W1+W2*(rho2-2*rhos)|rho2^{beta}>
  !  pot3, &     ! (beta+1) * rho2^{beta} *  <rho1^{alpha}|W2>
  !  pot4        ! beta * rho2^{beta-1} * <W1+W2*(rho1-2*rhos)|rho1^{alpha}>
  !
  !INTEGER :: i, j, k
  ! ------------------------------- TEST -------------------------------------


  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: &
    FrhoX, &    ! Fourier transform of rho^alpha or beta.
    FtrhoX, &   ! F-transform of theta* rho^alpha or beta - only when WGCT=2
    trafo0

                       !>> INITIALIZATION <<!

  CALL StartClock('WGCPotPlus')


  IF (vacCutoff) THEN
    cutf = CutoffFunc(rho)
    cutfp = CutoffFuncD(rho)
  END IF
  
                       !>> FUNCTION BODY <<!

  IF (firstOrderWGC==1) THEN  ! First order Taylor expanded WGC

    !--------------------------------------------------------------
    ! DERIVATIVE OF WGC ENERGY WITH RESPECT TO RHO^BETA(r^{PRIME}) 
    !--------------------------------------------------------------
    rhoX = rho**alpha                  ! rho^alpha in real space. 
    CALL FFT_NEW(FFT_STD_STATE,rhoX,FrhoX)                  ! rho^alpha in G space. FFT 1
    workRA = (rho - 2._DP*rhoS) * rhoX ! first order expansion in real space multiplied by rho^alpha 

    ! rhoB may be regulated in vacuum
    IF (vacCutoff) THEN
      rhoX = rho**beta * cutf
    ELSE
      rhoX = rho**beta
    END IF

    ! FFT 2; FFT 3; FFT4
    CALL FFT_NEW(FFT_STD_STATE,workRA,trafo0)
    CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,1)*FrhoX+keKernel(:,:,:,2)*trafo0,potential)
    CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,2)*FrhoX,workRA)

    ! Calculate energy in real space.
    ! rhoX = rho**beta
    IF (calcEnergy) THEN
      energy = cTF * SUM(rhoX * (potential + rho*workRA))
    ENDIF

    ! Calculate potential (part 1)
    IF (vacCutoff) THEN
      potential = (beta*rho**(beta-1._DP)*cutf + rho**beta*cutfp) * potential &
                + ((beta+1._DP)*rhoX + cutfp*rho**(beta+1._DP)) * workRA

      ! ------------------------------- TEST -------------------------------------
      ! pot4 = beta*rho**(beta-1._DP) * potential
      ! pot3 = (beta+1._DP)*rhoX * workRA
      ! ------------------------------- TEST -------------------------------------
    ELSE

      ! ------------------------------- TEST -------------------------------------
      ! pot4 = beta*rho**(beta-1._DP) * potential
      ! pot3 = (beta+1._DP)*rhoX * workRA
      ! ------------------------------- TEST -------------------------------------

      potential = beta*rho**(beta-1._DP) * potential &
                  + (beta+1._DP)*rhoX * workRA
    END IF


    !-----------------------------------------------------------------
    ! CHANGE TO THE DERIVATIVE OF E(WGC) WITH RESPECT TO RHO^ALPHA(r)
    !-----------------------------------------------------------------
    ! Setup more reused variables (1 FFT)
    ! rhoX = rho**beta
    CALL FFT_NEW(FFT_STD_STATE,rhoX,FrhoX)                  ! FrhoB
    workRA = (rho - 2._DP*rhoS) * rhoX

    ! Temporary storage (3 FFTs)
    CALL FFT_NEW(FFT_STD_STATE,workRA,trafo0)
    CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,1) * FrhoX + keKernel(:,:,:,2) * trafo0,rhoX)
    CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,2) * FrhoX, workRA)! Save memory by temp storage here
         ! above comment ignores the temporaries created (JMD)

    ! Second half of potential
    potential = potential +  alpha * rho**(alpha-1._DP) * rhoX &
                + (alpha+1._DP) * rho**alpha * workRA

    ! ------------------------------- TEST -------------------------------------
!    pot2 = alpha * rho**(alpha-1._DP) * rhoX
!    pot1 = (alpha+1._DP) * rho**alpha * workRA
!    WRITE(outputUnit,*) " potential of WGC"
!    DO i=1, SIZE(rho,1)
!      DO j=1, SIZE(rho,2)
!        DO k=1, SIZE(rho,3)
!          IF(i==9 .AND. j==9) THEN
!           WRITE(outputUnit,'(I3,6ES20.2)') k, &
!             cTF*pot1(i,j,k),cTF*pot2(i,j,k),cTF*pot3(i,j,k),cTF*pot4(i,j,k), &
!             cTF*potential(i,j,k), rho(i,j,k)
!          ENDIF
!          WRITE(outputUnit,'(3I3,4F10.5,2ES10.2,1F10.5)') &
!            i,j,k,pot1(i,j,k),pot2(i,j,k),pot3(i,j,k),pot4(i,j,k),&
!            rho(i,j,k),rhoS,&
      !      pot1(i,j,k)+pot2(i,j,k)+pot3(i,j,k)+pot4(i,j,k), &
!            potential(i,j,k)
      !    END IF
!        END DO
!      END DO
!    END DO
    ! ------------------------------- TEST -------------------------------------

  ELSE IF(firstOrderWGC==-1) THEN ! Second order Taylor expanded WGC

    ! Setup reused variables (2 FFTs)
    rhoX = rho**alpha        ! rhoA
    CALL FFT_NEW(FFT_STD_STATE,rhoX,FrhoX)        ! FrhoA
    theta = rho - rhoS
    CALL FFT_NEW(FFT_STD_STATE,theta*rhoX,FTrhoX) ! FtrhoA

    ! Temporary storage and intermediate computations (4 FFTs)
    rhoX = theta*theta * rhoX/2._DP
    IF (WGCT == -3) THEN 
      CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,1) * FrhoX + keKernel(:,:,:,2) * FtrhoX,workRA)
    ELSE
      CAlL FFT_NEW(FFT_STD_STATE,rhoX,trafo0)
      CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,1) * FrhoX + keKernel(:,:,:,2) * FtrhoX + & 
          keKernel(:,:,:,3) * trafo0,workRA)
    ENDIF
    CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,2) * FrhoX + keKernel(:,:,:,4) * FtrhoX,workRA2)

    IF ( WGCT == -3 ) THEN
      potential = 0.d0
    ELSE
      CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,3) * FrhoX,potential)
    ENDIF

    ! rhoB may be regulated in vacuum
    IF (vacCutoff) THEN
      rhoX = rho**beta * cutf
    ELSE
      rhoX = rho**beta
    END IF

    ! Calculate energy
    IF (calcEnergy) THEN
      energy = cTF * SUM(rhoX*(workRA + theta*workRA2 + theta*theta/2._DP*potential))
    ENDIF


    ! Sum first half of the potential
    IF (vacCutoff) THEN
      potential = rho**(beta-1._DP) * &
                    ((beta*cutf + rho*cutfp) &
                      * (workRA - rhoS*workRA2 + rhoS*rhoS/2._DP*potential) + &
                     ((beta+1._DP)*rho*cutf + rho**2*cutfp) &
                      * (workRA2 - rhoS*potential) + &
                     ((beta+2._DP)*rho**2*cutf+rho**3*cutfp) * potential/2._DP)
    ELSE
      potential = rho**(beta-1._DP) * &
                    (beta*workRA + (rho+beta*theta)*workRA2 &
                     + theta*(rho+beta*theta/2._DP)*potential)

    END IF

    ! Setup more reused variables (2 FFTs)
    CALL FFT_NEW(FFT_STD_STATE,rhoX,FrhoX)
    CALL FFT_NEW(FFT_STD_STATE,theta * rhoX,FtrhoX)

    ! More intermediate computations (4 FFTs)
    rhoX = theta*theta * rhoX/2._DP
    IF (WGCT==-3) THEN
       CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,1) * FrhoX + keKernel(:,:,:,2) * FtrhoX,workRA)
    ELSE
       CALL FFT_NEW(FFT_STD_STATE,rhoX,trafo0)
       CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,1) * FrhoX + keKernel(:,:,:,2) * FtrhoX &
                  + keKernel(:,:,:,3) * trafo0,workRA)
    ENDIF
    CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,2) * FrhoX + keKernel(:,:,:,4) * FtrhoX,workRA2)

    IF (WGCT==-3) THEN 
      rhoX = 0.d0
    ELSE
      CALL FFT_NEW(FFT_STD_STATE,keKernel(:,:,:,3) * FrhoX,rhoX)
    ENDIF

    ! Calculate second half of potential
    potential = potential + rho**(alpha-1._DP) * &
                  (alpha*workRA + (rho+alpha*theta)*workRA2 &
                   + theta*(rho+alpha*theta/2._DP)*rhoX)

  ELSE
    CALL WrtOut(6, 'WGCHessianForRho: unknown firstOrderWGC value, STOP')
    STOP
  END IF

  potential = potential * cTF

  CALL StopClock('WGCPotPlus')

END SUBROUTINE WGCPotentialPlus


FUNCTION CutoffFuncOld(rho)
!---------------------------------------------------------------
! DESCRIPTION: 
! compute the CutoffFunction which is used for WGCvac (
! short for WGC + vacuum), As we know WGC potentil divergence
! as density -> 0 due to the fact that beta<1
! with this cutoff function which go to zero very fast when 
! rho->0, we can avoid the divegence of WGC pot in vacuum
!
!   with CutFunc(r) = ( Exp[rho/rhos] - 1 ) / ( Exp[rho/rhos] + Exp[rhov/rhos] )
! 
!   CutFunc(r->inf) --> 1
!   CutFUnc(r->0  ) ~ r + o(r^2)
!
!   rhov and rhos are two parameters to be fixed, rhov is called the 
!   vacuum density default value is 1e-5
!   rhos = 1e-6,
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
  
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=dp), DIMENSION(n1G, n2G, n3G) :: CutoffFuncOld
  REAL(KIND=dp), DIMENSION(n1G, n2G, n3G) :: power 
  REAL(KIND=dp), DIMENSION(n1G, n2G, n3G) :: const 

  ! >> FUNCTION BEGINS << !
  const = exp(rhov/rhoStep)
  power = rho/rhoStep
  
  WHERE (power<1E-2) 
    power = rho/(rhoStep+const*rhoStep) &   ! the cutoff function is taylor expanded to 2nd order
       + (-1._DP+const)*rho**2/(2._DP*(1._DP+const)**2*rhoStep**2)
  ELSE WHERE (power>100._DP) 
    power = exp(100._DP)  ! rho is much bigger
  ELSE WHERE (power<=100._DP .AND. power>=1.E-2)
    power = exp(power)   
  END WHERE

  CutoffFuncOld = (power-1._DP)/(power+const)


!! DEBUG
!  WRITE(*,*) 'rhoStep, rhoV=>',rhoStep,rhoV
!  WRITE(*,*) 'exponent, rho(1,1,1)=>',rho(1,1,1)/rhoStep, rho(1,1,1)
!  WRITE(*,*) exp(rho(1,1,1)/rhoStep)
!  stop
!! ENDDEBUG
  
END FUNCTION CutoffFuncOld


FUNCTION CutoffFunc(rho)
!---------------------------------------------------------------
! compute the CutoffFunction which is used for WGCvac (
! short for WGC + vacuum), As we know WGC potentil divergence
! as density -> 0 due to the fact that beta<1
! with this cutoff function which go to zero very fast when 
! rho->0, we can avoid the divegence of WGC pot in vacuum
!
!   with CutFunc(r) = ( Exp[rho/rhos] - 1 ) / ( Exp[rho/rhos] + Exp[rhov/rhos] )
! 
!   CutFunc(r->inf) --> 1
!   CutFUnc(r->0  ) ~ r + o(r^2)
!
!   rhov and rhos are two parameters to be fixed, rhov is called the 
!   vacuum density default value is 1e-5
!   rhos = 1e-6,
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 06/22/2013 Mohan Chen fix a bug in the cutoff function 
!------------------------------------------------------------------------------
  IMPLICIT NONE
  
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho ! Electron density in real space (spin-independent)
  REAL(KIND=dp), DIMENSION(n1G, n2G, n3G) ::  &
    CutoffFunc, & 
    power, & 
    const

  ! >> FUNCTION BEGINS << !
  const = exp(rhov/rhoStep)
  power = rho/rhoStep
  
  WHERE (power<1E-2) 
    CutoffFunc = rho/(rhoStep+const*rhoStep) &   ! the cutoff function is taylor expanded to 2nd order
       + (-1._DP+const)*rho**2/(2._DP*(1._DP+const)**2*rhoStep**2)
  ELSE WHERE (power>100._DP) 
    power = exp(100._DP)  ! rho is much bigger
    CutoffFunc = (power-1._DP)/(power+const)
  ELSE WHERE (power<=100._DP .AND. power>=1.E-2)
    power = exp(power)   
    CutoffFunc = (power-1._DP)/(power+const)
  END WHERE

END FUNCTION CutoffFunc



FUNCTION CutoffFuncD(rho)
!-----------------------------------------------------------------------------
! Compute the d(CutoffFunction)/d(rho), used for WGCvac (
! short for WGC + vacuum), As we know WGC potentil divergence
! as density -> 0 due to the fact that beta<1
! with this cutoff function which go to zero very fast when 
! rho->0, we can avoid the divegence of WGC pot in vacuum
!
!   with CutFunc(r) = (1+Exp[-n0/D])(1/1+Exp[-(rho-n0)/D])-1/(1+Exp[n0/D])
!
!   CutFunc(r->inf) --> 1
!   CutFUnc(r->0  ) ~ r + o(r^2)
!
!   n0 and D are two parameters to be fixed, n0 is called the 
!   vacuum density default value is 1e-5
!   D = 1e-6,
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
!   Aug/24/2008   Created (Chen Huang)
!------------------------------------------------------------------------------
  IMPLICIT NONE
  
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: CutoffFuncD
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: power
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: const 

  ! >> FUNCTION BEGINS << !
  power = rho/rhoStep
  const = exp(rhov/rhoStep)
  
  WHERE (power<1e-2)     ! rho is very small
    CutoffFuncD = 1._DP/(rhoStep+const*rhoStep) &   ! Tylor expand the formula
     + (-1._DP+const)*rho/((1._DP+const)**2*rhoStep**2)
  ELSE WHERE (power>100._DP)
    power = exp(100._DP)
    CutoffFuncD = power*(1._DP+const)/((power+ const)**2*rhoStep)
  ELSE WHERE (power>=1.E-2 .AND. power<=100._DP)
    power = exp(power)
    CutoffFuncD = power*(1._DP+const)/((power+ const)**2*rhoStep)
  ENDWHERE
  
END FUNCTION CutoffFuncD


FUNCTION CutoffFuncDold(rho)
!-----------------------------------------------------------------------------
! DESCRIPTION: 
! Compute the d(CutoffFunction)/d(rho), used for WGCvac (
! short for WGC + vacuum), As we know WGC potentil divergence
! as density -> 0 due to the fact that beta<1
! with this cutoff function which go to zero very fast when 
! rho->0, we can avoid the divegence of WGC pot in vacuum
!
!   with CutFunc(r) = (1+Exp[-n0/D])(1/1+Exp[-(rho-n0)/D])-1/(1+Exp[n0/D])
!
!   CutFunc(r->inf) --> 1
!   CutFUnc(r->0  ) ~ r + o(r^2)
!
!   n0 and D are two parameters to be fixed, n0 is called the 
!   vacuum density default value is 1e-5
!   D = 1e-6,
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
  
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: CutoffFuncDold
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: power
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: const 

  ! >> FUNCTION BEGINS << !
  power = rho/rhoStep
  const = exp(rhov/rhoStep)
  
  WHERE (power<1e-2)     ! rho is very small
    power = 1._DP/(rhoStep+const*rhoStep) &   ! Tylor expand the formula
     + (-1._DP+const)*rho/((1._DP+const)**2*rhoStep**2)
  ELSE WHERE (power>100._DP)
    power = exp(100._DP)
  ELSE WHERE (power>=1.E-2 .AND. power<=100._DP)
    power = exp(power)
  ENDWHERE
  
  CutoffFuncDold = power*(1._DP+const)/((power+ const)**2*rhoStep)
  
END FUNCTION CutoffFuncDold


FUNCTION WGCStress(rhoReal)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Here's the function that computes the WGC portion of the stress and its
! contribution to the total stress tensor. It is calculated a lot less often
! than the energy or potential, so one routine is enough. Either way, the
! computation is a pain anyways.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REFERENCES:
! Robin Hayes Thesis. 
!------------------------------------------------------------------------------

  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu
  USE Fourier_NEW, ONLY: FFT_NEW
  USE Fourier_NEW, ONLY: FFT_STD_STATE
  USE PlaneWave, ONLY: qVectors
  USE MathFunctions, ONLY: Vecmul
  USE MathFunctions, ONLY: LindG ! Returns the Lindhard G function.
  USE MathFunctions, ONLY: GPrime ! Returns the first derivative of G.
  USE Constants, ONLY: PI 

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3,3) :: WGCStress
  ! The answer
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoReal
  ! Spin-independant density in real space.
  !

                         !>> INTERNAL VARIABLES <<!
  INTEGER :: a, b   
  ! counters
  !
  REAL(KIND=DP) :: c00, c01, c02, c11
  ! See Robin's group meeting for the expressions.
  !
  REAL(KIND=DP) :: c001, c002, c003, c011, c012, c021, c111
  ! Subparameters
  !
  REAL(KIND=DP) :: d00, d01, d02, d11
  ! Idem.
  !
  REAL(KIND=DP) :: diff
  ! comes from differentiating eta* when rho* is not held.
  !
  REAL(KIND=DP) :: eta
  ! Dimensionless momentum (norm of q over twice the Fermi vector)
  !
  REAL(KIND=DP) :: diago
  ! Temporary storage for the diagonal term of the stress.
  !
  REAL(KIND=DP) :: offdia
  ! Same for the off-diagonal contribution.
  !
  REAL(KIND=DP) :: ki
  ! Multiplicative constant before G(eta) and G'(eta).
  !
  REAL(KIND=DP) :: kF
  ! Fermi wave-vector based on rhoS instead of rho0 (as in WT)
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FrhoA
  ! FFT of rho^alpha
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FrhoB
  ! FFT of rho^beta
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FrhoA1
  ! FFT of rho^(alpha+1)
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FrhoB1
  ! FFT of rho^(beta+1)
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FrhoA2
  ! FFT of rho^(alpha+2)
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FrhoB2
  ! FFT of rho^(beta+2)
  !

  INTEGER :: ix, i2, i3  
  ! counter for parsing q over the recip sphere
  !
  REAL(KIND=DP), DIMENSION(3) :: q  
  ! cartesian coordinates of the q vector
  !
  CHARACTER(len=500) :: message
  !

                           !>> INITIALIZATION <<!

  CALL StartClock('WGCStress')

  WGCStress = 0._DP
  kF = (3._DP * pi**2 * rhoS)**(1._DP/3._DP)


  ! use new fft
  CALL FFT_NEW(FFT_STD_STATE, rhoReal**alpha, FrhoA)
  FrhoA = CONJG(FrhoA)
  CALL FFT_NEW(FFT_STD_STATE, rhoReal**beta, FrhoB)
  CALL FFT_NEW(FFT_STD_STATE, rhoReal**(alpha+1), FrhoA1)
  FrhoA1 = CONJG(FrhoA1)
  CALL FFT_NEW(FFT_STD_STATE, rhoReal**(beta+1), FrhoB1)
  CALL FFT_NEW(FFT_STD_STATE, rhoReal**(alpha+2), FrhoA2)
  FrhoA2 = CONJG(FrhoA2)
  CALL FFT_NEW(FFT_STD_STATE, rhoReal**(beta+2), FrhoB2)

! old fft
!  FrhoA = CONJG(FFT(rhoReal**alpha))
!  FrhoB = FFT(rhoReal**beta)
!  FrhoA1 = CONJG(FFT(rhoReal**(alpha+1)))
!  FrhoB1 = FFT(rhoReal**(beta+1))
!  FrhoA2 = CONJG(FFT(rhoReal**(alpha+2)))
!  FrhoB2 = FFT(rhoReal**(beta+2))

  ki = 20._DP * rho0**(5._DP/3._DP - (alpha + beta))

  ! Here's the tricky part: rho0 and rho* are functions of the volume unless
  ! they are fixed, and therefore give rise to stress terms. But only if they
  ! are not held constant, so we need to distinguish the cases.
  IF (hold0) THEN
    c001 = 1._DP - alpha - beta
    c002 = 2._DP * rhoS * (alpha + beta - 1._DP)
    c003 = rhoS**2 * (1._DP - alpha - beta)
    c011 = -alpha - beta
    c012 = rhoS * (alpha + beta)
    c021 = -(alpha + beta + 1._DP) / 2._DP
    c111 = -(alpha + beta + 1._DP)
  ELSE
    c001 = -2._DP / 3._DP
    c002 = rhoS * 4._DP / 3._DP
    c003 = -rhoS**2 * 2._DP / 3._DP
    c011 = -5._DP / 3._DP
    c012 = rhoS * 5._DP / 3._DP
    c021 = -4._DP / 3._DP
    c111 = -8._DP / 3._DP
  END IF
  IF (holdS) THEN
    diff = 0._DP
  ELSE
    c011 = c011 + 1._DP
    c012 = c012 - rhoS
    c021 = c021 + 1._DP
    c111 = c111 + 2._DP
    diff = -1._DP / 3._DP
  END IF


                           !>> FUNCTION BODY <<!

  ! Loop through all q vectors in the half-sphere. 
  DO i3 = 1, k3G 
    DO i2 = 1, k2G
      DO ix = 1, k1G
        ! Go further only if this point matters (within the KE cutoff)
        IF (qMask(ix,i2,i3)) THEN
          ! Calculate the qPoint cartesian coordinates given by this mVector.
          q(:) = qVectors(ix,i2,i3,:) 
          eta = qTable(ix,i2,i3) / (2._DP * kF)

          IF (WGCT == 2) THEN 
            ! First compute the diagonal term
            c00 = c001 * keKernel(ix,i2,i3,1) + c002 * keKernel(ix,i2,i3,2) + &
                  c003 * (keKernel(ix,i2,i3,3) + keKernel(ix,i2,i3,4))
            c01 = c011 * keKernel(ix,i2,i3,2) + &
                  c012 * (keKernel(ix,i2,i3,3) + keKernel(ix,i2,i3,4))
            c02 = c021 * keKernel(ix,i2,i3,3)
            c11 = c111 * keKernel(ix,i2,i3,4)
          ELSEIF(WGCT==-3) THEN 
            c00 = c001 * keKernel(ix,i2,i3,1) + c002 * keKernel(ix,i2,i3,2) + &
                  c003 * (keKernel(ix,i2,i3,4))
            c01 = c011 * keKernel(ix,i2,i3,2) + &
                  c012 * (keKernel(ix,i2,i3,4))
            c02 = c021 * 0.d0          
            c11 = c111 * keKernel(ix,i2,i3,4)
          ELSE 
            WRITE(message,*)"(FunctionalKinetic): WGC Stress is not implemented for WGCT = ", & 
              WGCT, " code stop!" 
            CALL WrtOut(6,message)
            STOP
          ENDIF
            
          diago = FrhoA(ix,i2,i3) * c00 * FrhoB(ix,i2,i3) &
                + FrhoA1(ix,i2,i3) * c01 * FrhoB(ix,i2,i3) &
                + FrhoA(ix,i2,i3) * c01 * FrhoB1(ix,i2,i3) &
                + FrhoA2(ix,i2,i3) * c02 * FrhoB(ix,i2,i3) &
                + FrhoA(ix,i2,i3) * c02 * FrhoB2(ix,i2,i3) &
                + FrhoA1(ix,i2,i3) * c11 * FrhoB1(ix,i2,i3)

          ! Now do the same for the off-diagonal term
          d11 = -ki * eta / (36._DP * rhoS**2) * GPrime(eta, mu) &
             - keKernel(ix,i2,i3,4) * 6._DP * (alpha + beta) &
             - keKernel(ix,i2,i3,2) * (gamma * (alpha + beta) + 6._DP * alpha &
                                     * beta) / rhoS
          IF (WGCT==2)  THEN                            
            d02 = -ki * eta / (72._DP * rhoS**2) * GPrime(eta, mu) &
               - keKernel(ix,i2,i3,4) * (3._DP * (alpha+beta+1._DP) - gamma) &
               - keKernel(ix,i2,i3,2) * (3._DP * gamma * (alpha+beta+1._DP) - &
               gamma**2 + 18._DP * alpha * beta) / (6._DP * rhoS)
          ELSE IF( WGCT==-3) THEN 
            d02 = 0.d0
          ELSE 
            WRITE(message,*)"(FunctionalKinetic): WGC Stress is not implemented for WGCT = ", & 
              WGCT, " code stop!" 
            CALL WrtOut(6,message)
            STOP
          ENDIF

             
          ! Temporary store -dF[w'AB]/deta* in d01 to use in d00. 
          d01 = keKernel(ix,i2,i3,4) * 6._DP*rhoS + keKernel(ix,i2,i3,2)*gamma 
          
          d00 = 6._DP * rhoS * keKernel(ix,i2,i3,2) - 2 * rhoS * d01 &
                + rhoS**2 * (2._DP * d02 + d11)
                
          d01 = -rhoS * (d11 + 2._DP * d02) + d01
          
          offdia = FrhoA(ix,i2,i3) * d00 * FrhoB(ix,i2,i3) &
                 + FrhoA1(ix,i2,i3) * d01 * FrhoB(ix,i2,i3) &
                 + FrhoA(ix,i2,i3) * d01 * FrhoB1(ix,i2,i3) &
                 + FrhoA2(ix,i2,i3) * d02 * FrhoB(ix,i2,i3) &
                 + FrhoA(ix,i2,i3) * d02 * FrhoB2(ix,i2,i3) &
                 + FrhoA1(ix,i2,i3) * d11 * FrhoB1(ix,i2,i3)
          ! and now for the 6 stress terms
          DO a = 1, 3
            WGCStress(a,a) = WGCStress(a,a) + diago + offdia * diff
            DO b = a, 3
              WGCStress(a,b) = WGCStress(a,b) + offdia * &
                q(a) * q(b) / qTable(ix,i2,i3)**2
            END DO ! b
          END DO ! a
        END IF ! in sphere
      END DO ! ix
    END DO ! i2
  END DO ! i3

  ! Now multiply by the prefactor
  ! 2 accounts for the double counting of plane waves
  WGCStress = WGCStress * 2._DP * cTF

  CALL StopClock('WGCStress')

  RETURN

END FUNCTION WGCStress

END MODULE KEDF_WGC

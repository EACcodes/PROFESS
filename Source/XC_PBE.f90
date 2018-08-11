MODULE XC_PBE 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE XC_PBE 
!     |_SUBROUTINE PBEPot
!     |_FUNCTION PBESpin_LibXC
!     |_FUNCTION PBE_LibXC
!     |_FUNCTION PBEStress
!
! DESCRIPTION:
!   Calculates the potential contibutions to the forces, stress, and energies.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
!   [2] libxc and algorithm in ABINIT rhohxc,drivexc,libxc_functional, xcpot, xcmult
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/15/2003  File Created by consolidating existing energy files (GSH) 
!   12/18/2003  Changed SUBROUTINE to FUNCTIONS
!   01/07/2004  Changed qTable and qMask and most instances of rhoR to 
!               spin-independent forms.
!   02/05/2004  Added nine stress functions. (VLL)
!   12/14/2007  Added Choly-Kaxiras methods for ion-electron terms.  (LH)
!   2013-08-05 Redesigned by Mohan Chen
!
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!
  USE CONSTANTS, ONLY: DP  
  USE CONSTANTS, ONLY: PI
  USE CONSTANTS, ONLY: IMAG 
  USE CONSTANTS, ONLY: TINY
  USE PlaneWave, ONLY: qVectors
  USE FOURIER, ONLY : FFT ! The Fast Fourier Transform routine.
  USE Fourier_NEW
  USE SetupFFT
  USE CellInfo, ONLY: k1G, k2G, k3G

  IMPLICIT NONE

  REAL(KIND=DP) :: pbeCutoff = 0.d0  
  ! used in PBEPot in make it 
  ! stable in low density region, we should not use it in principle
  ! I do not know why PBE potenital is not stable

CONTAINS

SUBROUTINE PBEPot(rhoReal_SI, potential, calcEnergy, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function computes the exchange-correlation potential in the GGA
!   (PBE) based on the real-space electron density.  This version works for
!   NO SPIN POLARIZATION ONLY.
!
!   Complete calculation to handle spin-polarized cases.  We tried here to 
!   minimize the number of temporary arrays necessary at the expense of excess
!   calculation.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   This thing uses a whole lot of temporary arrays
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/16/2003   Function created.  (VLL)
!   10/24/2008   Totally rewritten by Chen Huang, remove the numerically 
!                instabilty of the previous code when rho -> zero, this should
!                be the final code for PBE xc functional.
!   10/22/2013   Not final. Reoptimized w.r.t. memory usage, not finished (JMD)
!
!------------------------------------------------------------------------------

  USE MathFunctions, ONLY: MinMaxVal                

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoReal_SI
  ! Electron density in real space, spin INDEPENDENT
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential
  ! The XC potential
  !
  LOGICAL, INTENT(IN) :: calcEnergy
  !
  REAL(KIND=DP), INTENT(OUT) :: energy
  !

                     !>> INTERNAL VARIABLES <<! 

  REAL(KIND=DP), parameter :: &
    kappa = 0.8040_DP, &  ! The constant used in exchange, equation 10
    beta = 0.06672455060314922_DP, &
    pi2 = pi**2, &
    mu = beta * (pi2 / 3._DP)  , &
    a = 0.0310907_DP, &   ! Parameters in the uniform e-density correlation e.
    a1 = 0.21370_DP, &
    b1 = 7.5957_DP, &
    b2 = 3.5876_DP, &
    b3 = 1.6382_DP, &
    b4 = 0.49294_DP, & 
    onethird = 1._DP/3._DP, & 
    twothird = 2._DP/3._DP , & 
    fourthird = 4._DP/3._DP, & 
    seventhird = 7._DP/3._DP,   &
    eightthird = 8._DP/3._DP, &
    threesecond = 3._DP/2._DP

  REAL(KIND=DP) :: &
                          ! The constant part of exunif in equation 10
                          ! of [1].
    coeff,&
    gamma

  REAL(KIND=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), &
                           SIZE(rhoReal_SI,3)) :: &
    ecunif,  &
    decunif_dRho, &
    s2,      & 
    tmp1, tmp0, &
    dHdt2,   & 
    tmpRho, &           ! Everything will be done on this tmpRho
    Rho13, &            ! tmpRho^{1/3}
    Rho43               ! tmpRho^{4/3}

  REAL(KIND=DP) :: &
    bigB, &
    rs


  REAL(KIND=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), &
                           SIZE(rhoReal_SI,3), 3) :: &
    exGrad, &
    gradRho

  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: rhoRecip,trafo1,trafo2

  INTEGER :: i
  REAL(KIND=DP) :: coeff1, coeff2, coeff3, coeff4, coeff5
  REAL(KIND=DP) :: coeff6, coeff7
  REAL(KIND=DP) :: pre1, pre2, pre3, pre4, pre5, pre6, pre7, pre8
  INTEGER :: n1, n2, n3
  INTEGER :: ir1, ir2, ir3
  REAL(KIND=DP) :: ff, gg
  REAL(KIND=DP) :: rs_half
  REAL(KIND=DP) :: kf_without_rho ! K_{F}=(3pi^{2}rho)^{1/3},this is the coefficient without rho
  REAL(KIND=DP) :: ecut2, grad2
  REAL(KIND=DP) :: Eec
  REAL(KIND=DP) :: tmpX,tmpXX

                      !>> INITIALIZATION <<!
    
    ! If density becomes very small, PBE subroutine will not be stable
    ! (generate NaN, and Infinity). Therefore we make a hard cutoff here, 
    ! only aims to make this subroutine correcly behave
    ! we use 1e-20 as the bottom line for our density, which will affect
    ! nothing in normal cases.

    
    ! Start the clock
    CALL StartClock('PBEPot')
!    CALL StartClock('PBEtest1')

    ! (1.1) Get the three dimensions of charge density in real space
    n1 = SIZE(rhoReal_SI,1)
    n2 = SIZE(rhoReal_SI,2)
    n3 = SIZE(rhoReal_SI,3)
!    write(*,*) "n1=",n1,"n2=",n2,"n3=",n3

    ! (1.2) Set the energy and 3D potential array to zero.
    ! potential(n1,nr,n3)
    energy = 0._DP
    potential = 0._DP

    ! (1.3) Get the temporary charge density, and eliminate the
    ! very small terms.
    tmpRho = rhoReal_SI
    where ( rhoReal_SI < 1e-20 ) 
      tmpRho = 1e-20
    endwhere
  
    ! (1.4) Transform the charge from real space to G space.
    ! rhoRecip is charge density in reciprocal space,
    ! the dimension of charge is less than (n1,n2,n3) in G space
    CALL FFT_NEW(FFT_STD_STATE,tmpRho,rhoRecip)

    ! (1.5) see equation(6) in ref.1
    gamma = (1._DP - LOG(2._DP)) / pi2 

    ! (1.6) First calculate the Rho13=rho^{1/3}
    ! and the Rho43=rho^{4/3}
    ! both will be used many times in the following
    Rho13 = tmpRho**onethird
    Rho43 = tmpRho*Rho13

    ! (1.7) Calculate the derivative of charge in G space, and then FFT back to real space,
    ! 3 FFTs are needed here.
    ! TODO I am really upset about the imag here
    trafo1 = imag*qVectors(:,:,:,1)*rhoRecip
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,1))
    trafo1 = imag*qVectors(:,:,:,2)*rhoRecip
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,2))
    trafo1 = imag*qVectors(:,:,:,3)*rhoRecip
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,3))

    ! (1.8) the s(r)^2 parameter in PBE paper, eqn(9) from ref.1
    ! because s = |nabla_rho(r)| / 2K_{F}rho(r)
    ! so s(r)^2 = |nabla_rho(r)|^2 / 4K_{F}^{2}rho(r)^{2}
    !           = |nabla_rho(r)|^2 / 4*pre1*rho(r)^{8/3}   
    kf_without_rho=4._DP * (3._DP*pi2)**twothird
    ! In order to compare the |nabla_rho(r)|^2 and the energy cutoff
    ecut2=pbeCutoff*pbeCutoff
    DO ir1=1, n1
      DO ir2=1, n2
        DO ir3=1, n3
          grad2=gradRho(ir1,ir2,ir3,1)*gradRho(ir1,ir2,ir3,1) &
               +gradRho(ir1,ir2,ir3,2)*gradRho(ir1,ir2,ir3,2) &
               +gradRho(ir1,ir2,ir3,3)*gradRho(ir1,ir2,ir3,3)
          IF( grad2 < ecut2 ) THEN
            gradRho(ir1,ir2,ir3,:)=0.0_DP
          ELSE
            s2(ir1,ir2,ir3)=grad2 / ( kf_without_rho * ( Rho43(ir1,ir2,ir3)*Rho43(ir1,ir2,ir3)))
          END IF
        END DO
      END DO
    END DO
    
!    CALL StopClock('PBEtest1')
!    CALL StartClock('PBEtest2')

                      ! >> FUNCTION BODY <<!

    !-------------------------------------------------
    ! Compute the exhange potential first
    !-------------------------------------------------

    ! (2.1) the Fx for exchange part, eqn(14) from ref.1
    ! tmp1 = Fx
    ! tmp0 = dFx/d(s^2)
    pre1=1.0_DP+kappa
    pre2=mu/kappa
    DO ir1=1, n1
      DO ir2=1, n2
        DO ir3=1, n3
          pre3=1.0_DP+pre2*s2(ir1,ir2,ir3)
          tmp1(ir1,ir2,ir3) = pre1-kappa/pre3
          tmp0(ir1,ir2,ir3) = mu/(pre3*pre3)  
        END DO
      END DO
    END DO
    
    ! (2.2) Calculate the Exchange Potential and correlation potential
    coeff = -3._DP/(4._DP*pi)*(3._DP*pi2)**onethird
    coeff1 = coeff * fourthird
    coeff2 = coeff *(-eightthird)
    potential = potential + coeff1 * Rho13 * tmp1
    potential = potential + coeff2 * Rho13 * tmp0 * s2
 
    ! (2.3) Calculate exchange energy
    IF (calcEnergy) energy = SUM(Rho43 * coeff * tmp1)

!    CALL StopClock('PBEtest2')

    ! (2.4) calculate the correlation potential
    ! the potential due to rho(r)*epsion^unif part
    pre1=2.0_DP*a
    pre2=b1*0.5_DP
    pre3=threesecond*b3
    pre4=b4*2._DP
    pre5=(3._DP/4._DP/pi)**onethird

    DO ir1 = 1, n1
      DO ir2 = 1, n2
        DO ir3 = 1, n3

           rs = pre5 / Rho13(ir1,ir2,ir3)
           rs_half = sqrt(rs)

           tmpX = pre1*(b1*rs_half + rs*(b2+b3*rs_half+b4*rs))
           tmpXX = pre1*(pre2/rs_half + b2 + pre3*rs_half + pre4*rs)

           coeff1 = LOG(1._DP+1._DP/tmpX)!tmp3(ir1,ir2,ir3))
           coeff2 = 1._DP+a1*rs

           ! the ec^unif in PBE (prl) paper
           ecunif(ir1,ir2,ir3) = -pre1*coeff1*coeff2

           decunif_dRho(ir1,ir2,ir3) = pre1*(-a1*coeff1 + coeff2/(tmpX*(1+tmpX)) *tmpXX) &
                * ( -rs * onethird / tmpRho(ir1,ir2,ir3) )
        
        END DO
      END DO
    END DO

!    CALL StopClock('PBEtest2')
!    CALL StartClock('PBEtest3')

    ! (3.1) Calculate the potential and the energy
    potential = potential + ecunif + tmpRho*decunif_dRho
    if (calcEnergy) energy = energy + SUM(tmpRho*ecunif)
    
    ! (3.1) Now to compute the potential due to 
    ! rho*H(r,xi,t) in PBE correlation part
    
    ! we convert the s^2 to t^2 now, 
    ! the t^2 is in equation(7) in ref 1
    coeff7 = pi*(3._DP*pi2)**onethird/4._DP
    s2 = s2*coeff7*Rho13
    
    ! d(t^2)/d(rho)
    coeff7 = -7._DP/3._DP
    tmp1 = coeff7*s2/tmpRho
    
    pre5=gamma/beta
    pre6=-1.d0/beta
    pre7=beta/gamma

    DO ir1 = 1, n1
      DO ir2 = 1, n2
        DO ir3 = 1, n3
    
         ! The capital A in eqn(8) in Ref[1] will blow up (become Infinity)
         ! wherever density becomes very small (large vacuum case),
         ! because small rho -> small r_s -> ecunif->0 -> A->1/0.
         ! To void this blows up of capital A, we define bigB as
         ! the inverse of A (=1/A), and B will will always be finite in low density
         ! regions
          Eec = EXP(-ecunif(ir1,ir2,ir3)/gamma)
          bigB = pre5*(Eec-1._DP)

          coeff1 = bigB * (bigB + s2(ir1, ir2, ir3))
          coeff4 = 1.0/(coeff1 + s2(ir1, ir2, ir3)*s2(ir1, ir2, ir3))
          tmpXX = tmp1(ir1,ir2,ir3)

          ! tmp1 will be used to calculate the energy
          ! s2 has been conveted to t2 in Eq(7) of PBE (PRL paper)
          tmp1(ir1,ir2,ir3) = coeff1 * coeff4
          coeff3 = tmp1(ir1,ir2,ir3) * coeff4

          ff = 2.d0*bigB+s2(ir1,ir2,ir3)

          tmpX = bigB*coeff4-(bigB+2.d0*s2(ir1,ir2,ir3))*coeff3

          coeff5 = 1.0_DP + (pre7*s2(ir1,ir2,ir3)*tmp1(ir1,ir2,ir3))

          gg = beta / coeff5

          ! dH / dt2
          dHdt2(ir1,ir2,ir3) = gg*(tmpX*s2(ir1,ir2,ir3)+tmp1(ir1,ir2,ir3))

          ! d_H /d_BigB, re-use tmp2 here to save memory
          tmpX = gg*s2(ir1,ir2,ir3)*ff * ( coeff4 - coeff3 )

          potential(ir1, ir2, ir3) = potential( ir1, ir2,ir3)+gamma*LOG(coeff5) + &
             tmpRho(ir1,ir2,ir3)*(dHdt2(ir1,ir2,ir3)*tmpXX+tmpX* &
             ( pre6 * Eec * decunif_dRho(ir1,ir2,ir3) ))

        END DO
      ENDDO
    ENDDO
        
!    CALL StopClock('PBEtest3')
!    CALL StartClock('PBEtest4')
    
    ! compute energy for the 2nd part in correlation
    coeff1 = beta/gamma
    if (calcEnergy) energy = energy + SUM(tmpRho*gamma*LOG(1._DP+coeff1*s2*tmp1))

    ! prepare the dEx/dGradRho for later usage
    coeff6 = (-3._DP/(4._DP*pi)*(3._DP*pi**2)**onethird)/(2._DP*(3._DP*pi**2)**twothird)

    ! dFx/d(gradRho)
    !d/d(gradient of Rho) part
    pre8 = (pi/3._DP)**onethird/8._DP
    do i = 1, 3 
      ! rho * d(H)/d(t^2) * d(t^2)/d(gradRho)
      ! and then add up the XC part
      exGrad(:,:,:,i) = ( coeff6 * tmp0 + dHdt2 * pre8 ) * gradRho(:,:,:,i) / Rho43
    enddo 

    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,1),trafo1)
    trafo1 = trafo1*imag*qVectors(:,:,:,1) ! TODO REMOVE THE IMAG (STORE QVECTORS WITH IT)
    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,2),trafo2)
    trafo1 = trafo1 + trafo2*imag*qVectors(:,:,:,2)
    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,3),trafo2)
    trafo1 = trafo1 + trafo2*imag*qVectors(:,:,:,3)
    CALL FFT_NEW(FFT_STD_STATE,trafo1,tmp0)

    ! Final part of the potential
    potential = potential - tmp0
!    CALL StopClock('PBEtest4')

    CALL StopClock('PBEPot')

END SUBROUTINE PBEPot


SUBROUTINE PBE_LibXC(rho, potential, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is a general version of PBE using libxc and work for both spin
!   unpolarized and polarized cases. The algo is the same as ABINIT libxc PBE
!   subroutine.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!  02/28/2014 : Created By Steven Xia
!  03/04/2014 : Ported to new FFT, general optimizations (JMD)
!------------------------------------------------------------------------------
  USE CellInfo, ONLY: numSpin 
  ! Number of spins in calculation (1 or 2) 
 
  USE xc_f90_types_m
  USE xc_f90_lib_m
  USE libxc_funcs_m

  IMPLICIT NONE

  REAL(KIND=DP),DIMENSION(:,:,:,:),INTENT(IN):: rho
  REAL(KIND=DP),DIMENSION(:,:,:,:),INTENT(OUT):: potential
  REAL(KIND=DP),INTENT(OUT):: energy

  !working variables
  !gradup2 = gradup * gradup
  !sigma(:,:,:,1) = gradup*gradup
  !sigma(2)=gradup*graddn
  !simga(3)=graddn*graddn
  !vsigma = dE/dsigma
  !depsxc(:,:,:,1) = dE/dup
  !depsxc(:,:,:,2) = dE/ddn
  !depsxc(:,:,:,3) = 1/|gradup|*dE/d|gradup|
  !depsxc(:,:,:,4) = 1/|graddn|*dE/d|graddn|
  !depsxc(:,:,:,5) = 1/|gradtot|*dE/d|gradtot|
  !rhonow(:,:,:,1,1) = up
  ! if numspin = 1
  !rhonow(:,:,:,1,2) = gradtot/|gradtot|*dE/d|gradtot|
  ! if numspin = 2
  !rhonow(:,:,:,1,2) = (gradup/|gradup|*dE/d|gradup| + gradtot/|gradtot| * dE/d|gradtot|)xdir
  !rhonow(:,:,:,1,3) = (gradup/|gradup|*dE/d|gradup| + gradtot/|gradtot| * dE/d|gradtot|)ydir
  !rhonow(:,:,:,1,4) = (gradup/|gradup|*dE/d|gradup| + gradtot/|gradtot| * dE/d|gradtot|)zdir
  !rhonow(:,:,:,2,:) is for spin_dn and the fifth dimension is the same as described above

  !rhotmp,vctmp,vxtmp,sigmatmp,vxsigmatmp,vcsigmatmp,extmp,ectmp are temp variables for libxc calling
  
  ! local variables (-> please at some point move them to heap [JMD])
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: rhoRecip,trafo1,trafo2
  REAL(KIND=DP),DIMENSION(size(rho,1),size(rho,2),size(rho,3)):: rhotot,up,dn,excpbe
  REAL(KIND=DP),DIMENSION(size(rho,1),size(rho,2),size(rho,3),3):: gradup,graddn,gradtot,sigma
  REAL(KIND=DP),DIMENSION(size(rho,1),size(rho,2),size(rho,3),5):: depsxc
  REAL(KIND=DP),DIMENSION(size(rho,1),size(rho,2),size(rho,3),size(rho,4),4):: rhonow     
  REAL(KIND=DP)::rhotmp(size(rho,4)),vctmp(size(rho,4)),vxtmp(size(rho,4))
  REAL(KIND=DP)::sigmatmp(3),vxsigmatmp(3),vcsigmatmp(3)
  REAL(KIND=DP)::extmp,ectmp

  INTEGER:: i,j,k,ispin
  
  TYPE(xc_f90_pointer_t) :: x_func, c_func
  TYPE(xc_f90_pointer_t) :: x_info, c_info

  !! >> FUNCTION << !!
     
  !   write(*,*) "USING NEW SPIN LIBXC!!!!"

  !first get parameter needed for libxc
  If(numSpin==1) then
      rhotot = rho(:,:,:,1)
  Else
      up = rho(:,:,:,1)
      dn = rho(:,:,:,2)
      rhotot = up + dn
  Endif
  
  !compute gradtot
  CALL FFT_NEW(FFT_STD_STATE,rhotot,rhoRecip)
  
  ! TODO I am really upset about the imag here (JMD)
  trafo1 = imag*qVectors(:,:,:,1)*rhoRecip
  CALL FFT_NEW(FFT_STD_STATE,trafo1,gradtot(:,:,:,1))
  trafo1 = imag*qVectors(:,:,:,2)*rhoRecip
  CALL FFT_NEW(FFT_STD_STATE,trafo1,gradtot(:,:,:,2))
  trafo1 = imag*qVectors(:,:,:,3)*rhoRecip
  CALL FFT_NEW(FFT_STD_STATE,trafo1,gradtot(:,:,:,3))
  
  If(numSpin==1) then
      sigma(:,:,:,1) = gradtot(:,:,:,1)*gradtot(:,:,:,1) &
                      +gradtot(:,:,:,2)*gradtot(:,:,:,2) &
                      +gradtot(:,:,:,3)*gradtot(:,:,:,3)
      sigma(:,:,:,2) = 0.d0
      sigma(:,:,:,3) = 0.d0
  Else
      ! compute gradup, gradup2
      CALL FFT_NEW(FFT_STD_STATE,up,rhoRecip)
      trafo1 = imag*qVectors(:,:,:,1)*rhoRecip
      CALL FFT_NEW(FFT_STD_STATE,trafo1,gradup(:,:,:,1))
      trafo1 = imag*qVectors(:,:,:,2)*rhoRecip
      CALL FFT_NEW(FFT_STD_STATE,trafo1,gradup(:,:,:,2))
      trafo1 = imag*qVectors(:,:,:,3)*rhoRecip
      CALL FFT_NEW(FFT_STD_STATE,trafo1,gradup(:,:,:,3))
      
      sigma(:,:,:,1) = gradup(:,:,:,1)*gradup(:,:,:,1) &
                      +gradup(:,:,:,2)*gradup(:,:,:,2) &
                      +gradup(:,:,:,3)*gradup(:,:,:,3)
      ! compute graddn, graddn2
      CALL FFT_NEW(FFT_STD_STATE,dn,rhoRecip)
      trafo1 = imag*qVectors(:,:,:,1)*rhoRecip
      CALL FFT_NEW(FFT_STD_STATE,trafo1,graddn(:,:,:,1))
      trafo1 = imag*qVectors(:,:,:,2)*rhoRecip
      CALL FFT_NEW(FFT_STD_STATE,trafo1,graddn(:,:,:,2))
      trafo1 = imag*qVectors(:,:,:,3)*rhoRecip
      CALL FFT_NEW(FFT_STD_STATE,trafo1,graddn(:,:,:,3))
      
      sigma(:,:,:,3) = graddn(:,:,:,1)*graddn(:,:,:,1) &
                      +graddn(:,:,:,2)*graddn(:,:,:,2) &
                      +graddn(:,:,:,3)*graddn(:,:,:,3)
      
      !compute gradup*graddn
      sigma(:,:,:,2) = gradup(:,:,:,1)*graddn(:,:,:,1) &
                      +gradup(:,:,:,2)*graddn(:,:,:,2) &
                      +gradup(:,:,:,3)*graddn(:,:,:,3)
                      
  Endif
  
  !use libxc to get dE/drho_up/dn and dE/dsigma
  call xc_f90_func_init(x_func, x_info, XC_GGA_X_PBE, numspin)
  call xc_f90_func_init(c_func, c_info, XC_GGA_C_PBE, numspin)
  
  ! JMD: it may be possible to get more speed out of this by reordering the loops.
  do k = 1,size(rho,3)
      do j = 1,size(rho,2)
          do i = 1,size(rho,1)
              rhotmp = rho(i,j,k,:)
              sigmatmp = sigma(i,j,k,:)
              call xc_f90_gga_exc_vxc(x_func,1,rhotmp(1),sigmatmp(1),extmp,vxtmp(1),vxsigmatmp(1))
              call xc_f90_gga_exc_vxc(c_func,1,rhotmp(1),sigmatmp(1),ectmp,vctmp(1),vcsigmatmp(1))
              
              excpbe(i,j,k) =  extmp + ectmp
              depsxc(i,j,k,1) = vxtmp(1) + vctmp(1)
         
              If(numSpin==1) then
                  depsxc(i,j,k,5) =  2.d0 * (vcsigmatmp(1) + vxsigmatmp(1))
              Else
                  depsxc(i,j,k,2) = vxtmp(2) + vctmp(2)
                  depsxc(i,j,k,3) =  2.d0*(vxsigmatmp(1)+vcsigmatmp(1)) - vxsigmatmp(2) - vcsigmatmp(2)
                  depsxc(i,j,k,4) =  2.d0*(vxsigmatmp(3)+vcsigmatmp(3)) - vxsigmatmp(2) - vcsigmatmp(2)
                  depsxc(i,j,k,5) =  vxsigmatmp(2) + vcsigmatmp(2)
              Endif
          enddo
      enddo
  enddo
  
  ! now already get depsxc, calculate rhonow
  If(numspin==1) then
      rhonow(:,:,:,1,1) = rhotot
      rhonow(:,:,:,1,2) = depsxc(:,:,:,5) * gradtot(:,:,:,1)
      rhonow(:,:,:,1,3) = depsxc(:,:,:,5) * gradtot(:,:,:,2)
      rhonow(:,:,:,1,4) = depsxc(:,:,:,5) * gradtot(:,:,:,3)   
  Else
      rhonow(:,:,:,1,1) = up(:,:,:)
      rhonow(:,:,:,2,1) = dn(:,:,:)
      
      rhonow(:,:,:,1,2) = gradup(:,:,:,1)*depsxc(:,:,:,3)+gradtot(:,:,:,1)*depsxc(:,:,:,5)
      rhonow(:,:,:,1,3) = gradup(:,:,:,2)*depsxc(:,:,:,3)+gradtot(:,:,:,2)*depsxc(:,:,:,5)
      rhonow(:,:,:,1,4) = gradup(:,:,:,3)*depsxc(:,:,:,3)+gradtot(:,:,:,3)*depsxc(:,:,:,5)


      rhonow(:,:,:,2,2) = graddn(:,:,:,1)*depsxc(:,:,:,4)+gradtot(:,:,:,1)*depsxc(:,:,:,5)
      rhonow(:,:,:,2,3) = graddn(:,:,:,2)*depsxc(:,:,:,4)+gradtot(:,:,:,2)*depsxc(:,:,:,5)
      rhonow(:,:,:,2,4) = graddn(:,:,:,3)*depsxc(:,:,:,4)+gradtot(:,:,:,3)*depsxc(:,:,:,5)
  Endif

  ! now already get rhonow, do divergence calculation 
  DO ispin =1,size(rho,4)
      
      CALL FFT_NEW(FFT_STD_STATE,rhonow(:,:,:,ispin,2),trafo1)
      CALL FFT_NEW(FFT_STD_STATE,rhonow(:,:,:,ispin,3),trafo2)
      trafo1 = imag*qVectors(:,:,:,1)*trafo1 + &
               imag*qVectors(:,:,:,2)*trafo2
      CALL FFT_NEW(FFT_STD_STATE,rhonow(:,:,:,ispin,4),trafo2)
      trafo1 = trafo1 + imag*qVectors(:,:,:,3)*trafo2
      
      !trafo1 = imag*qVectors(:,:,:,1)*trafo1 + &
      !         imag*qVectors(:,:,:,2)*trafo2 + &   
      !         imag*qVectors(:,:,:,3)*trafo3
               
      CALL FFT_NEW(FFT_STD_STATE,trafo1,potential(:,:,:,ispin))
      
      potential(:,:,:,ispin) = -potential(:,:,:,ispin) + depsxc(:,:,:,ispin) ! now vscpbe already contain the gradient derivative term
                                                        ! add the dE/drho_updn term
                                  
  ENDDO
  
  !now energy
  energy = SUM(rhotot*excpbe)
       
  ! debug:
  !write(*,*) "New PBE max potential is: ", maxval(potential)
  !call PBEPot(rho(:,:,:,1),potential(:,:,:,1),.true.,energy)
  !write(*,*) "old PBE max potential is: ", maxval(potential)

  RETURN

END SUBROUTINE PBE_LibXC


FUNCTION PBEStress(rhoReal_SI, rhoRecip_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the exchange-correlation stress component 
!   specified by a and b in the local density approximation.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Further optimization is almost definately possible
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/05/2004  Function created.  (Vincent Ligneres)
!
!------------------------------------------------------------------------------
  USE Constants, ONLY: auToGPa

  USE CellInfo, ONLY : m123G
  USE CellInfo, ONLY : cell, numSpin

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoReal_SI
  ! Electron density in real space, spin INDEPENDENT
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoRecip_SI
  ! Electron density in real space, spin INDEPENDENT
  !
  REAL(KIND=DP), DIMENSION(3,3) :: PBEStress
  ! The final answer


                      !>> INTERNAL VARIABLES <<!
  REAL(KIND=DP), PARAMETER :: &
                          ! The constant part of exunif in equation 10
                          ! of [1].
    kappa = 0.8040_DP, &  ! The constant used in exchange, equation 10
    beta = 0.06672455060314922_DP, &
    mu = beta * (pi**2 / 3._DP), &
    a = 0.0310907_DP, &   ! Parameters in the uniform e-density correlation e.
    a1 = 0.21370_DP, &
    b1 = 7.5957_DP, &
    b2 = 3.5876_DP, &
    b3 = 1.6382_DP, &
    b4 = 0.49294_DP

  INTEGER :: &
    i, j   ! The two indices that define the stress component we're looking at.

  REAL(KIND=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), &
                           SIZE(rhoReal_SI,3)) :: &
    ecunif, &
    decdr, &
    t_2, &                ! In the exchange, this is s^2, 
                          !  in the correlation, t^2
    rtrs, &               ! SQRT(rs)
    modGradRho, &
    temp1, temp2

  REAL(KIND=DP) :: &
    cst, &
    cKf, &
    ceXUnif, &
    gamma, &
    delta, &              ! this is not in the paper, but we define it thus
    cKs, &                ! Thomas-fermi screening wavenumber
    tempReal

  REAL(KIND=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), &
                           SIZE(rhoReal_SI,3), 3) :: &
    gradRho

                           !>> INITIALIZATION <<!

  CALL StartClock("PBEStress")

  cst = (3._DP * pi**2 / 16._DP)**(2._DP/3._DP)
  cKf = (3._DP * pi**2)**(1._DP/3._DP)
  ceXUnif = -0.75_DP * (cKf / pi)
  gamma = (1._DP - LOG(2._DP)) / pi**2
  delta = beta / gamma     ! this is not in the paper, but we define it thus
  cKs = (4._DP * cKf/ pi)**(1._DP/2._DP) ! Thomas-fermi screening wavenumber

  gradRho(:,:,:,1) = -FFT(imag*qVectors(:,:,:,1)*rhoRecip_SI)
  gradRho(:,:,:,2) = -FFT(imag*qVectors(:,:,:,2)*rhoRecip_SI)
  gradRho(:,:,:,3) = -FFT(imag*qVectors(:,:,:,3)*rhoRecip_SI)

  modGradRho = SQRT(gradRho(:,:,:,1)**2 + gradRho(:,:,:,2)**2 + gradRho(:,:,:,3)**2)

                           !>> FUNCTION BODY <<!
  ! Do we have a spin-neutral or a spin-polarized density?
  SELECT CASE (numSpin)

  ! Spin-neutral case.
  CASE(1)

    ! First calculate the exchange stress
    ! This is s_2
    t_2 = (modGradRho/(2._DP * cKf * rhoReal_SI**(4._DP/3._DP)))**2

    ! Temporary values to assist computation
    temp1 = (1._DP + mu * t_2 / kappa)
    temp2 = -ceXUnif*mu / (2._DP * cKf**2 * rhoReal_SI**(4._DP/3._DP) &
                             * temp1**2)
    tempReal = SUM(ceXunif * rhoReal_SI**(4._DP/3._DP) / 3._DP &
           * ((2._DP * mu * t_2 / temp1**2) - (1._DP + kappa - kappa / temp1)))

    DO i = 1, 3
      DO j = i, 3
        PBEStress(i,j) = SUM(temp2 *  gradRho(:,:,:,i) * gradRho(:,:,:,j))
        IF (i==j) THEN
          PBEStress(i,j) = PBEStress(i,j) + tempReal
        END IF
      END DO
    END DO

    ! Now calculate the correlation stress
    rtrs = (3._DP / (4._DP * pi * rhoReal_SI))**(1._DP/6._DP)
    t_2 = cst * t_2 / rtrs**2
    temp1 = 2._DP * a * rtrs*(b1+rtrs*(b2+rtrs*(b3+rtrs*b4)))
    temp2 = 1._DP + a1 * rtrs**2
 
    ecunif = -2._DP * a * temp2 * LOG(1._DP + (1._DP / temp1))
    decdr = 2._DP * a**2 / 3._DP * (temp2 &
            * rtrs*(b1+rtrs*(2._DP*b2+rtrs*(3._DP*b3+rtrs*4._DP*b4))) &
            / ((1._DP + temp1)*temp1) &
            - a1*rtrs**2/a * LOG(1._DP + 1._DP / temp1))

    temp1 = delta/(EXP(-ecunif/gamma) - 1._DP)*t_2       ! At^2
    temp2 = 1._DP + temp1 +temp1**2                      ! 1 + At^2 + (At)^2
    temp2 = temp2**2 + delta * t_2*(1._DP + temp1)*temp2

    tempReal = SUM(rhoReal_SI &
                * (decdr + (beta*t_2*(1._DP + 2._DP * temp1)/3._DP &
                    - temp1**3 * (2._DP + temp1) * EXP(-ecunif/gamma) * decdr)&
                /temp2))

    temp2 = - (1._DP + 2._DP * temp1) * beta / (2._DP * cKs**2 &
            * rhoReal_SI**(4._DP/3._DP) * temp2)

 
    DO i = 1, 3
      DO j = i, 3
        PBEStress(i,j) = PBEStress(i,j) &
                         + SUM(temp2 *  gradRho(:,:,:,i) * gradRho(:,:,:,j))
        IF (i==j) THEN
          PBEStress(i,j) = PBEStress(i,j) + tempReal
        END IF
      END DO
    END DO

    ! This is the 1 / Volume * dV part
    PBEStress = PBEStress / REAL(m123G, KIND=DP)

  CASE DEFAULT
    WRITE(*,*) 'PBEStress error: Can only handle one spin.'
    STOP

  END SELECT


!  WRITE(*,*) " stress for PBE (Gpa): "
!  WRITE(*,'(3F12.6)') PBEStress(1,:)*auToGPa
!  WRITE(*,'(3F12.6)') PBEStress(2,:)*auToGPa
!  WRITE(*,'(3F12.6)') PBEStress(3,:)*auToGPa

  CALL StopClock("PBEStress")

  RETURN

END FUNCTION PBEStress

END MODULE XC_PBE

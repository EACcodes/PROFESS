MODULE XC_LDA 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE XC_LDA 
!     |_FUNCTION LDAEnergy                  (Local Density Approx. exch.-corr.)
!     |_FUNCTION LDAPot
!     |_SUBROUTINE LSDAPotPZ
!     |_SUBROUTINE LSDAPotPW92
!     |_SUBROUTINE spn_gcor
!     |_FUNCTION LSDAStress
!     |_FUNCTION LDAStress
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1]  Perdew J.P. and Zunger A., Phys. Rev. B 23(10), 1981, 5048-79
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/15/2003  File Created by consolidating existing energy files (GSH) 
!   12/18/2003  Changed SUBROUTINE to FUNCTIONS
!   01/07/2004  Changed qTable and qMask and most instances of rhoR to 
!               spin-independent forms.
!   02/05/2004  Added nine stress functions. (VLL)
!   12/14/2007  Added Choly-Kaxiras methods for ion-electron terms.  (LH)
!   2013-08-05 Redesigned by mohan
!
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!


  USE CONSTANTS, ONLY: DP
  USE CONSTANTS, ONLY: PI

  IMPLICIT NONE

  INTEGER :: exchangeCorrelation=1
  ! Specifies the XC functional: 
  !  (1) Local Density Approximation (default)
  !  (2) Generalized Gradient Approx (PBE version)
  !
  INTEGER :: lsda = 1
  ! Specifies the LSDA XC functional used for spin-polarzied calculation
  !  (1) Perdew-Zunger parametrization (default)
  !  (2) Perdew-Wang parametrization

CONTAINS

FUNCTION LDAEnergy(rhoReal)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the exchange-correlation energy in Local Density
!   Approximation (LDA) based on the real-space electron density.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Complete calculation to handle spin-polarized cases.
!   I could also change this to use LDAPointEnergy, for cleanness
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/01/2003  File created.  (VLL)
!   06/19/2006  Optimized (GSH)
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP) :: LDAEnergy
  ! The XC energy.                     
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rhoReal
  ! Electron density in real space, spin DEPENDENT
  !

                     !>> INTERNAL VARIABLES <<! 

  REAL(KIND=DP), PARAMETER :: &
    ft = 4.0_DP / 3.0_DP, &              ! four thirds
    ot = 1.0_DP/3.0_DP, &                ! one third
    mot = -1.0_DP/3.0_DP, &              ! negative (minus) one third
    cX = -0.73855876638202234_DP, &      ! -0.75_DP * (3._DP/pi)**ot    
    cC = 0.62035049089940009_DP, &       ! (7.5E-1_DP / pi)**ot
    fzz = 1.709921d0, &                     ! For PW92 LSDA
    gamma = 0.5198421d0, &
    !ax = -0.7385588d0, &
    one =1.d0, &
    two = 2.d0, &
    three = 3.d0, &
    four = 4.d0

  REAL(KIND=DP) :: &
    exch, &                              ! Running exchange energy
    rs, &                                ! (3/(4.pi.rho))^1/3
    corr, &                              ! Running correlation energy
    zet, &                             ! Spin-polarization
    eC, &                                 ! exchange and correlation energy
    exa, &
    exb, &
    eu, f,z4,eurs,dd, &  ! work variables
    ep,eprs,alfrsm,alfm,exc       ! work variables

  REAL(KIND=DP), DIMENSION(2), PARAMETER :: &
    a = (/3.11E-2_DP,1.555E-2_DP/), &      ! Parameter A for rs < 1 in [1]
    b = (/-4.8E-2_DP,-2.69E-2_DP/), &      ! Parameter B for rs < 1 in [1]
    c = (/2E-3_DP,7E-4_DP/), &             ! Parameter C (Ceperley-Alder) [1]
    d = (/-1.16E-2_DP,-4.8E-3_DP/), &      ! Parameter D (C-A) from [1]
    g = (/-1.423E-1_DP,-8.43E-2_DP/), &    ! Parameter gamma for rs > 1 [1]
    b1 = (/1.0529_DP,1.3981_DP/), &        ! Parameter beta1 for rs > 1 [1]
    b2 = (/3.334E-1_DP,2.611E-1_DP/)       ! Parameter beta2 for rs > 1 [1]

  INTEGER :: &
    ns, &                                ! Shorthand notation for num. spins
    ix, iy, iz                           ! Dummy counters

                           !>> INITIALIZATION <<!   
                           !>> FUNCTION BODY <<!

  ns=SIZE(rhoReal,4)

  ! Do we have a spin-neutral or a spin-polarized density?
  SELECT CASE (ns)

  ! Spin-neutral case.
  CASE(1)
    ! Exchange energy.
    exch = cX * SUM(rhoReal**ft)

    ! Correlation energy.
    corr = 0._DP
    DO iz=1, SIZE(rhoReal,3)
      DO iy=1, SIZE(rhoReal,2)
        DO ix=1, SIZE(rhoReal,1)
          ! Calculate r subscript s [1] in atomic units.
          rs = cC * rhoReal(ix,iy,iz,1)**mot 
          IF (rs<1) THEN
            eC = LOG(rs)*(a(1) + c(1) * rs) + b(1) + d(1) * rs
          ELSE
            eC = g(1) / (1 + b1(1) * SQRT(rs) + b2(1) * rs)
          END IF
          corr = corr + rhoReal(ix,iy,iz,1) * eC
        END DO
      END DO
    END DO
    LDAEnergy= exch + corr

  CASE(2)
     exc = 0.d0
     ! Perdew-Wang 1992 XC functional
     DO iz=1, SIZE(rhoReal,3)
        DO iy=1, SIZE(rhoReal,2)
           DO ix=1, SIZE(rhoReal,1)
              ! exchange
              dd=two* rhoReal(ix,iy,iz,1)
              exa = cx*dd**ot
              dd=two* rhoReal(ix,iy,iz,2)
              exb = cx*dd**ot
              ! local correlation
              dd=rhoReal(ix,iy,iz,1)+rhoReal(ix,iy,iz,2)
              !rs=(0.75d0/pi/dd)**ot
!              rs = cC * dd**mot
              rs=(0.75d0/(pi*dd))**ot
              call spn_gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0, &
                   1.6382d0,0.49294d0,1.00d0,rs,eu,eurs)
              call spn_gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0, &
                   3.3662d0,0.62517d0,1.00D0,rs,ep,eprs)
              call spn_gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0, &
                   0.88026d0,0.49671d0,1.00d0,rs,alfm,alfrsm)
              zet=(rhoReal(ix,iy,iz,1) - rhoReal(ix,iy,iz,2))/dd
              f = ((one+zet)**ft+(one-zet)**ft-two)/gamma
              z4 = zet**4
              ec = eu*(one-f*z4)+ep*f*z4-alfm*f*(one-z4)/fzz
!              ec=two*(rhoReal(ix,iy,iz,1)*exa+rhoReal(ix,iy,iz,2)*exb+ec*dd)
              ec=(rhoReal(ix,iy,iz,1)*exa+rhoReal(ix,iy,iz,2)*exb+ec*dd)
              exc = exc + ec
           END DO
        END DO
     END DO

     LDAEnergy= exc
     
  CASE DEFAULT
  
     write(*,*) "ERROR: Unknown case in LDAEnergy. Case is: ",ns
     flush(6)
     stop


  END SELECT

END FUNCTION LDAEnergy


FUNCTION LDAPot(rho)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function computes the exchange-correlation potential in the Local 
!   Density Approximation (LDA) based on the real-space electron density.
!
! CONDITIONS AND ASSUMPTIONS:
!   We can use LDAPointPot for moudularity.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Complete calculation to handle spin-polarized cases.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/16/2003  Function created.  (VLL)
!   06/19/2006  Optimized (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  ! Electron density in real space, spin DEPENDENT
  !
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2), &
                           SIZE(rho,3),SIZE(rho,4)) :: LDAPot
 ! The XC potential.

                     !>> INTERNAL VARIABLES <<! 


  REAL(KIND=DP), PARAMETER :: &
    mot = -1._DP/3._DP, &                ! minus one-third.
    ot = 1.0_DP/3.0_DP, &                ! one third
    cX = -0.98474502184269641_DP, &      ! Exchange constant:
                                         !  (4/3) * -0.75_DP * (3._DP/pi)**ot 
    ft = 4.0_DP / 3.0_DP, &              ! four thirds
    cC = 0.62035049089940009_DP, &       ! Correlation constant
    cC1 = -5.83666666666666666E-2_DP, &  ! (b(1) - a(1) * ot)
    cC2 = 1.33333333333333333E-3_DP, &   ! 2*ot * c(1)
    cC3 = -8.4E-3_DP, &                  ! ot * (2*d(1) - c(1)) 
    cC4 = -0.174798948333333333_DP, &    ! g(1) * (7._DP/6._DP) * b1(1)
    cC5 = -6.325709333333333333E-2_DP, & ! g(1) * 4._DP * ot * b2(1)
    fzz = 1.709921d0, &                     ! For PW92 LSDA
    gamma = 0.5198421d0, &
    ax = -0.7385588d0, &
    one =1.d0, &
    two = 2.d0, &
    three = 3.d0, &
    four = 4.d0

  REAL(KIND=DP) :: rs                    ! (3/(4.pi.rho))^1/3
  REAL(KIND=DP) :: sqrtRS
  REAL(KIND=DP), DIMENSION(2), PARAMETER :: &
    a = (/3.11E-2_DP,1.555E-2_DP/), &      ! Parameter A for rs < 1 in [1]
!    b = (/-4.8E-2_DP,-2.69E-2_DP/), &      ! Parameter B for rs < 1 in [1]
!    c = (/2E-3_DP,7E-4_DP/), &             ! Parameter C (Ceperley-Alder) [1]
!    d = (/-1.16E-2_DP,-4.8E-3_DP/), &      ! Parameter D (C-A) from [1]
    g = (/-1.423E-1_DP,-8.43E-2_DP/), &    ! Parameter gamma for rs > 1 [1]
    b1 = (/1.0529_DP,1.3981_DP/), &        ! Parameter beta1 for rs > 1 [1]
    b2 = (/3.334E-1_DP,2.611E-1_DP/)       ! Parameter beta2 for rs > 1 [1]

  INTEGER :: ns                           ! Shorthand notation for num. spins
  INTEGER :: ix, iy, iz                   ! Dummy counters

                           !>> FUNCTION BODY <<!
  ns=SIZE(rho,4)

  ! Do we have a spin-neutral or a spin-polarized density?
  SELECT CASE (ns)

  ! Spin-neutral case.
  CASE(1)
    ! Exchange potential. (will multiply by cX later)
    LDAPot = rho**ot

    ! Correlation potential.
    DO iz=1, SIZE(rho,3)
      DO iy=1, SIZE(rho,2)
        DO ix=1, SIZE(rho,1)

            ! Calculate r subscript s [1] in atomic units.
            rs = cC / LDAPot(ix,iy,iz,1)
            IF (rs<1) THEN

              LDAPot(ix,iy,iz,1) = cX * LDAPot(ix,iy,iz,1) + &
                LOG(rs) * (a(1) + cC2 * rs) + cC1 + cC3*rs

            ELSE
              sqrtRS = SQRT(rs)
              LDAPot(ix,iy,iz,1) = cX * LDAPot(ix,iy,iz,1) + &
                (g(1) + cC4 * sqrtRS + cC5 * rs) / &
                (1._DP + b1(1) * sqrtRS + b2(1) * rs)**2
            END IF

        END DO
      END DO
    END DO

  CASE(2)
    
    WRITE(*,*) "ERROR: ec would be used unintialized to compute a quantity in LDAPot (case 2)."
    WRITE(*,*) "ERROR: Obviously, this will yield garbage. STOP"
    FLUSH(6)
    STOP
    
  END SELECT

  RETURN

END FUNCTION LDAPot


SUBROUTINE LSDAPotPW92(rho, LDAPotential,  LDAEnergy)
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

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  ! Electron density in real space, spin DEPENDENT
  !
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2), &
       SIZE(rho,3),SIZE(rho,4)), INTENT(OUT) :: LDAPotential
  ! The XC potential.
  !
  REAL(KIND=DP), INTENT(OUT) :: LDAEnergy

                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), PARAMETER :: &
    p75vpi = 0.75d0/pi, &
    ax = -0.7385588d0, &
    fzz = 1.709921d0, &        
    gamma = 0.5198421d0, &
    one =1.d0, &
    two = 2.d0, &
    three = 3.d0, &
    four = 4.d0

  REAL(KIND=DP) :: &
    rs, &                              ! (3/(4.pi.rho))^1/3
    zet, &                             ! Spin-polarization
    ex(2), vx(2),exc, &                           ! exchange energy
    ec, vc(2), &                              ! correlation energy
    eu,f,z4,ecrs,eurs,eczet,comm,fz,d, &  ! work variables
    ep,eprs,alfrsm,alfm,ac2,third         ! work variables

  INTEGER :: ix, iy, iz
  ! Dummy counters

  CALL StartClock('LSDAPW92')

  third = one/three
  ac2 = four/three
  exc= 0.d0
  DO iz=1, SIZE(rho,3)
     DO iy=1, SIZE(rho,2)
        DO ix=1, SIZE(rho,1)
           ! exchange
           d=two*rho(ix,iy,iz,1)
           ex(1)= ax*d**third
           vx(1)=ex(1)*ac2
           d=two*rho(ix,iy,iz,2)
           ex(2)= ax*d**third
           vx(2)=ex(2)*ac2
           ! local correlation
           d=sum(rho(ix,iy,iz,:))
           rs=(p75vpi/d)**third
           call spn_gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0, &
                1.6382d0,0.49294d0,1.00d0,rs,eu,eurs)
           call spn_gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0, &
                3.3662d0,0.62517d0,1.00D0,rs,ep,eprs)
           call spn_gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0, &
                0.88026d0,0.49671d0,1.00d0,rs,alfm,alfrsm)
           zet=(rho(ix,iy,iz,1)-rho(ix,iy,iz,2))/d
           f = ((one+zet)**ac2+(one-zet)**ac2-two)/gamma
           z4 = zet**4
           ec = eu*(one-f*z4)+ep*f*z4-alfm*f*(one-z4)/fzz
           ecrs = eurs*(one-f*z4)+eprs*f*z4-alfrsm*f*(one-z4)/fzz
           fz = ac2*((one+zet)**third-(one-zet)**third)/gamma
           eczet = four*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu &
                -(one-z4)*alfm/fzz)
           comm = ec -rs*ecrs/three-zet*eczet
           vc(1) = comm + eczet
           vc(2) = comm - eczet

!           vxc(i,1)=two*(vx(1)+vc(1))
!           vxc(i,2)=two*(vx(2)+vc(2))
           ec=(rho(ix,iy,iz,1)*ex(1)+rho(ix,iy,iz,2)*ex(2)+ec*d)
           exc = exc + ec
           LDAPotential(ix,iy,iz,1)=(vx(1)+vc(1)) 
           LDAPotential(ix,iy,iz,2)=(vx(2)+vc(2)) 
        ENDDO
     ENDDO
  ENDDO

  LDAEnergy = exc

  CALL StopClock('LSDAPW92')

  RETURN

END SUBROUTINE LSDAPotPW92


SUBROUTINE LSDAPotPZ(rho, LDAPotential,  LDAEnergy)
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
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  ! Electron density in real space, spin DEPENDENT
  !
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2), &
       SIZE(rho,3),SIZE(rho,4)), INTENT(OUT) :: LDAPotential
  ! The XC potential.
  !
  REAL(KIND=DP), INTENT(OUT) :: LDAEnergy

                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), PARAMETER :: &
    p75vpi = 0.75d0/pi, &
    ax = -0.7385588d0, &
    fzz = 1.709921d0, &
    gamma = 0.5198421d0, &
    one =1.d0, &
    two = 2.d0, &
    three = 3.d0, &
    four = 4.d0

  REAL(KIND=DP) :: &
    vxc(2), &
    vs, &
    lrs, &
    es, &
    srs, &
    rs, &                              ! (3/(4.pi.rho))^1/3
    zet, &                             ! Spin-polarization
    ex(2), vx(2),exc, &                           ! exchange energy
    ec, vc(2), &                              ! correlation energy
    eu,f,fz,d, &  ! work variables
    ac2,third         ! work variables
  
  ! High density (rs<1) contants, LDA
  REAL(KIND=DP), PARAMETER :: ca_c1 = 0.0622d0
  REAL(KIND=DP), PARAMETER :: ca_c2 = 0.0960d0
  REAL(KIND=DP), PARAMETER :: ca_c3 = 0.0040d0
  REAL(KIND=DP), PARAMETER :: ca_c4 = 0.0232d0
  REAL(KIND=DP), PARAMETER :: ca_c5 = 0.0192d0

  REAL(KIND=DP), PARAMETER :: ca_as =  0.0311d0
  REAL(KIND=DP), PARAMETER :: ca_bs = -0.0538d0
  REAL(KIND=DP), PARAMETER :: ca_cs =  0.0014d0
  REAL(KIND=DP), PARAMETER :: ca_ds = -0.0096d0
  ! Low density (rs>1) constants, LDA
  REAL(KIND=DP), PARAMETER :: ca_g = -0.2846d0
  REAL(KIND=DP), PARAMETER :: ca_b1 = 1.0529d0
  REAL(KIND=DP), PARAMETER :: ca_b2 = 0.3334d0
  REAL(KIND=DP), PARAMETER :: ca_gs = -0.1686d0
  REAL(KIND=DP), PARAMETER :: ca_b1s = 1.3981d0
  REAL(KIND=DP), PARAMETER :: ca_b2s = 0.2611d0

  INTEGER :: ix, iy, iz
  ! Dummy counters

  CALL StartClock('LSDAPZ')

  third = one/three
  ac2 = four/three
  exc= 0.d0
  DO iz=1, SIZE(rho,3)
     DO iy=1, SIZE(rho,2)
        DO ix=1, SIZE(rho,1)
           ! exchange
           d=two*rho(ix,iy,iz,1)
           ex(1)= ax*d**third
           vx(1)=ex(1)*ac2
           d=two*rho(ix,iy,iz,2)
           ex(2)= ax*d**third
           vx(2)=ex(2)*ac2
           ! local correlation
           d=sum(rho(ix,iy,iz,:))
           rs=(p75vpi/d)**third

           zet=(rho(ix,iy,iz,1)-rho(ix,iy,iz,2))/d
           f = ((one+zet)**ac2+(one-zet)**ac2-two)/gamma
           fz = ac2*((one+zet)**third-(one-zet)**third)/gamma
           if (rs < one) then
           ! correlation: high density
              lrs = log(rs)
              eu = ca_c1 * lrs - ca_c2 + (ca_c3 * lrs - ca_c4) * rs  
              vc(1) = eu - (ca_c1 + (ca_c3 * lrs + ca_c3 - ca_c4) * rs) /3.d0
              es = ca_as * lrs + ca_bs + (ca_cs * lrs + ca_ds) * rs
              vs = es - (ca_as + (ca_cs * lrs + ca_cs + ca_ds) * rs) /3.d0
           else
              ! correlation: low density
              srs = sqrt(rs)
              eu = ca_g / (1.d0 + ca_b1 * srs + ca_b2 * rs)
              vc(1) = eu * eu * (1.d0 + 7.0d0 * ca_b1 * srs / 6.d0 + &
                   4.d0/3.d0 * ca_b2 * rs) / ca_g
              es = ca_gs / (1.d0 + ca_b1s * srs + ca_b2s * rs)
              vs = es * es * (1.d0 + 7.0d0 * ca_b1s * srs / 6.d0 + &
                   4.d0/3.d0 * ca_b2s * rs) / ca_gs
           endif
           vc(2) = vc(1) + f * (vs - vc(1)) - (es - eu) * &
                (1.d0 + zet) * fz
           vc(1) = vc(1) + f * (vs - vc(1)) + (es - eu) * &
                (1.d0 - zet) * fz
           
           ec = eu + f * (es - eu)
           
           vxc(1)=two*vx(1)+vc(1)
           vxc(2)=two*vx(2)+vc(2)
           ec=two*(rho(ix,iy,iz,1)*ex(1)+rho(ix,iy,iz,2)*ex(2))+ec*d
           exc = exc + ec
           LDAPotential(ix,iy,iz,1)= 0.5d0*vxc(1)
           LDAPotential(ix,iy,iz,2)= 0.5d0*vxc(2)
        enddo
     enddo
  enddo
  LDAEnergy = 0.5d0*exc

  CALL StopClock('LSDAPZ')
END SUBROUTINE LSDAPotPZ


SUBROUTINE spn_gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This subroutine computes the local correlation energy and
! potential for the Perdew-Wang exchange-correlation scheme.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  IMPLICIT NONE
  !
  ! Input/Output variables:
  !
  REAL(KIND=DP), INTENT(IN) :: a,a1,b1,b2,b3,b4,p,rs
  REAL(KIND=DP), INTENT(OUT) :: gg,ggrs
  !
  ! Work variables:
  !
  REAL(KIND=DP) :: p1,q0,rsp,q1,q2,q3,rs12,rs32,two, one,three
  !---------------------------------------------------------------
  one = 1.d0
  two = 2.d0
  three = 3.d0
  p1 = p + 1.d0
  q0 = -two*a*(one+a1*rs)
  rs12 = sqrt(rs)
  rs32 = rs12**3
  rsp = rs**p
  q1 = two*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
  q2 = log(one+one/q1)
  gg = q0*q2
  q3 = a*(b1/rs12+two*b2+three*b3*rs12+two*b4*p1*rsp)
  ggrs = -two*a*a1*q2-q0*q3/(q1**2+q1)

  RETURN

END SUBROUTINE spn_gcor


FUNCTION LSDAStress(cellVol, rhoR, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the exchange-correlation stress component
!   specified by a and b in the local density approximation.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: m123G
  ! m123G is global FFT grid size 
  USE CellInfo, ONLY: m3G
  ! m3G is global FFT grid size for Z
  USE CellInfo, ONLY: n3G
  ! n3G is local FFT grid size for Z

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(3,3) :: LSDAStress
  ! The final answer
  !
  REAL(KIND=DP), INTENT(IN) :: cellVol
  ! The volume of the cell
  !
  REAL(KIND=DP), INTENT(IN) :: energy
  ! The LDA energy
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rhoR
  ! Spin-independent electron density
  !

                      !>> INTERNAL VARIABLES <<!
  INTEGER :: a
  ! Index of stress component
  !
  REAL(KIND=DP) :: LDAEnergy
  !
  REAL(KIND=DP) :: LDAPotential(size(rhoR,1),size(rhoR,2),size(rhoR,3),2)
  !
  REAL(KIND=DP) :: factor1, factor2


                           !>> INITIALIZATION <<!

  ! n3G: local FFT dimension for Z
  ! m3G: global FFT dimension for Z
  factor1 = REAL(n3G,KIND=DP) / REAL(m3G, KIND=DP)
  factor2 = 1.0_DP / REAL(m123G, KIND=DP)

  LSDAStress = 0._DP
  CALL LSDAPotPZ(rhoR, LDAPotential,  LDAEnergy)
                           !>> FUNCTION BODY <<!
  DO a = 1, 3
     LSDAStress(a,a) = energy / cellVol * factor1 &
          - (SUM(rhoR(:,:,:,1)*LDAPotential(:,:,:,1))+ &
             SUM(rhoR(:,:,:,2)*LDAPotential(:,:,:,2))) * factor2
  END DO

END FUNCTION LSDAStress


FUNCTION LDAStress(cellVol, rhoR, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the exchange-correlation stress component 
!   specified by a and b in the local density approximation.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/05/2004  Function created.  (Vincent Ligneres)
!
!------------------------------------------------------------------------------

  USE Constants, ONLY: auToGPa

  USE CellInfo, ONLY: m123G
  ! m123G is global FFT grid size 
  USE CellInfo, ONLY: m3G
  ! m3G is global FFT grid size for Z
  USE CellInfo, ONLY: n3G
  ! n3G is local FFT grid size for Z

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3,3) :: LDAStress
  ! The final answer
  !
  REAL(KIND=DP), INTENT(IN) :: cellVol
  ! The volume of the cell
  !
  REAL(KIND=DP), INTENT(IN) :: energy
  ! The LDA energy
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rhoR
  ! Spin-independent electron density
  !

                      !>> INTERNAL VARIABLES <<!

  INTEGER :: a
  ! Index of stress component
  !
  REAL(KIND=DP) :: factor1, factor2
  !

                           !>> INITIALIZATION <<!

  LDAStress = 0._DP

  ! n3G: local FFT dimension for Z
  ! m3G: global FFT dimension for Z
  factor1 = REAL(n3G,KIND=DP) / REAL(m3G, KIND=DP)
  factor2 = 1.0_DP / REAL(m123G, KIND=DP)

                           !>> FUNCTION BODY <<!
  DO a = 1, 3
    LDAStress(a,a) = energy / cellVol * factor1 - SUM(rhoR*LDAPot(rhoR)) * factor2
  END DO

!  WRITE(*,*) " stress for LDA (Gpa): "
!  WRITE(*,'(3F12.6)') LDAStress(1,:)*auToGPa
!  WRITE(*,'(3F12.6)') LDAStress(2,:)*auToGPa
!  WRITE(*,'(3F12.6)') LDAStress(3,:)*auToGPa

END FUNCTION LDAStress

END MODULE XC_LDA

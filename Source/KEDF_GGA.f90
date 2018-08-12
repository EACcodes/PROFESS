MODULE KEDF_GGA
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_GGA
!     |_SUBROUTINE GGAPotentialPlus 
!     |_SUBROUTINE vWGTF
!------------------------------------------------------------------------------

  USE Constants, ONLY: DP, pi, imag
  USE Fourier, ONLY: FFT
  USE OutputFiles, ONLY: outputUnit

  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu

  USE PlaneWave, ONLY: qTable, qVectors, qMask
  USE CellInfo, ONLY: k1G, k2G, k3G

  IMPLICIT NONE
  
  Integer :: GGA_functional=1        ! model number for GGA functionals
  Integer :: model=1                 ! model number for vWGTF KEDFs
  Real(DP) :: CP=0.d0       ! penalty cooeffecient for the vW term, mostly zero except for TFvW models

  ! For vW+GTF model interpolation, Steven added
  Integer :: numint=200
  REAL(KIND=DP), DIMENSION(200) :: &
    d, &               ! rho/rho0 grid to be interpolated.
    ELFvsd, ELFvsdDD  ! ELF values

  REAL(KIND=DP) :: LKTa0 = 1.0

CONTAINS

SUBROUTINE GGAPotentialPlus(rho, potential, calcEnergy, energy, functionals) 
!------------------------------------------------------------------------------
! This is a subroutine calculate all kinds of GGA functionals. See reviews in 
! ref(1): JCTC 5, 3161 (2009), 
! ref(2): JCP 127, 144109 (2007)
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
!   July/24/2013   Created (Steven)
!------------------------------------------------------------------------------

  IMPLICIT NONE

                                !>> External variables <<!
  INTEGER, INTENT(in) :: &
    functionals          ! the flag telling which GGA funtional we are using
  
  LOGICAL, INTENT(in) :: &
    calcEnergy          ! if energy is calculated
  
  REAL(kind=DP), Dimension(:,:,:), INTENT(in) :: &
    rho                 ! total electron density
  
  REAL(kind=DP), INTENT(out) :: &
    energy              ! kinetic energy

  REAL(kind=DP), Dimension(:,:,:), INTENT(out) :: &
    potential           ! dE/drho, i.e., kinetic energy potential

                
                                !>> Internal variables <<!
  Integer :: &
    i, j, k
  
  Real(KIND=DP), parameter :: &
    cTF = 2.87123400018819_DP, &
    cs = 2.d0 * (3.d0*pi**2)**(1.d0/3.d0),& ! difference between two different s definitions
    twothird = 2.d0/3.d0, &
    fourthird = 4.d0/3.d0, &
    fivethird = 5.d0/3.d0

    ! Constants in functionals:
  REAL(KIND=DP) :: LLPb, LLPc, &
    OL1d, OL2D, &
    PWA1, PWA2, PWA3, PWA4, PWA, PWB1, &
    LCA1, LCA2, LCA3, LCA4, LCA, LCB1, &
    TWk, TWmu, &
    PBE2k, PBE2mu, &
    E00A1, E00A2, E00A3, E00A4, &
    P92A1, P92A2, P92A3, P92A4, &
    PW86A1, PW86A2, PW86A3, &
    DKA1, DKA2, DKA3,DKB1, DKB2, DKB3, &
    B86AC1, B86AC2, &
    B86BC1, B86BC2, &
    Thakb, Thakc, Thakd, &
    LG94A2, LG94A4, LG94A6, LG94A8, &
    LG94A10, LG94A12, LG94b, &
    DK87A, DK87B, DK87C, &
    LKTa, &
    ! penalty energy and its multiplication
    PenaltyE
  
  REAL(KIND=DP), Dimension(size(rho,1),size(rho,2),size(rho,3)) :: &
    temprho, rho23, rho83, rho43, &
    s, s2, x,&          ! reduced gradient, x in DK functional
    dsdrho, &           ! ds/drho
    tf, &               ! thomas fermi energy density
    dtfdrho, &          ! dtf/drho
    F, &                ! Enhancemence factor
    dFds, &             ! dF/ds
    PenaltyP            ! penalty potential term

  REAL(KIND=DP), Dimension(size(rho,1),size(rho,2),size(rho,3),3) :: &
    grad, &             ! density gradient
    dsdgrad, &          ! ds/dgradrho
    group               ! for calculating gradient use

  COMPLEX(kind=DP), DIMENSION(k1G, k2G, k3G) :: rhoRecip

                                !>> Initialization <<!
  tf = CTF * rho**(5.0_DP/3.0_DP)
  LLPb = 0.0044188*2.d0**(2.d0/3.d0)
  LLPc = 0.0253*2.d0**(1.d0/3.d0)
  OL1d = 0.00187*2.d0**(1.d0/3.d0)
  OL2D = 0.0245*2.d0**(1.d0/3.d0)
  PWA1 = 0.19645/cs
  PWA2 = 0.2743/cs**2
  PWA3 = 0.1508/cs**2
  PWA4 = 100.00/cs**2
  PWA = 7.7956/cs
  PWB1 = 0.004/cs**4
  LCA1 = 0.093907/cs
  LCA2 = 0.26608/cs**2
  LCA3 = 0.0809615/cs**2
  LCA4 = 100.00/cs**2
  LCA = 76.320/cs
  LCB1 = 0.57767e-4/cs**4
  TWk = 0.8438
  TWmu = 0.2319/cs
  PBE2k = 6.9031
  PBE2mu = 2.0309/cs**2
  LKTa = LKTa0/cs
  E00A1 = 135.d0
  E00A2 = 28.d0/cs**2
  E00A3 = 5.d0/cs**4
  E00A4 = 3.d0/cs**2
  P92A1 = 1.d0
  P92A2 = 88.3960d0/cs**2
  P92A3 = 16.3683/cs**4
  P92A4 = 88.2108/cs**2
  PW86A1 = 1.296d0/cs**2
  PW86A2 = 14.d0/cs**4
  PW86A3 = 0.2d0/cs**6
  DKA1 = 0.95
  DKA2 = 14.28111
  DKA3 = -19.5762
  DKB1 = -0.05
  DKB2 = 9.99802
  DKB3 = 2.96085
  B86AC1 = 0.0039*2.d0**(2.d0/3.d0)
  B86AC2 = 0.004*2.d0**(2.d0/3.d0)
  B86BC1 = 0.00403*2.d0**(2.d0/3.d0)
  B86BC2 = 0.007*2.d0**(2.d0/3.d0)
  Thakb = 0.0055*2.d0**(2.d0/3.d0)
  Thakc = 0.0253*2.d0**(1.d0/3.d0)
  Thakd = 0.072*2.d0**(1.d0/3.d0)
  LG94A2 = (1e-8+0.1234)/0.024974/cs**2
  LG94A4 = 29.790/cs**4
  LG94A6 = 22.417/cs**6
  LG94A8 = 12.119/cs**8
  LG94A10 = 1570.1/cs**10
  LG94A12 = 55.944/cs**12
  LG94b = 0.024974
  DK87A = 7.d0/324.d0/(18.d0*pi**4)**(1.d0/3.d0) * 2.d0**(2.d0/3.d0)
  DK87B = 0.861504*2.d0**(1.d0/3.d0)
  DK87C = 0.044286*2.d0**(2.d0/3.d0)
  
  ! remove super small density points
  tempRho = rho
  where ( rho < 1e-20 )
    tempRho = 1e-20
  endwhere

  Rho23 = tempRho**(2.d0/3.d0)
  Rho43 = Rho23*Rho23
  Rho83 = Rho43*Rho43

  ! calculate gradient of rho
  rhoRecip = FFT(temprho)
  grad(:,:,:,1) = FFT(imag*qVectors(:,:,:,1)*rhoRecip)
  grad(:,:,:,2) = FFT(imag*qVectors(:,:,:,2)*rhoRecip)
  grad(:,:,:,3) = FFT(imag*qVectors(:,:,:,3)*rhoRecip)

  ! calculate s
  s = sqrt( grad(:,:,:,1)**2 + grad(:,:,:,2)**2 + grad(:,:,:,3)**2 ) / rho43
  s2 = s*s

  ! calculate dtf/drho, ds/drho, and ds/dgrad
  dtfdrho = fivethird * CTF * rho23
  
  ! here ds/drho actually equal to ds/drho*s to cancel out 1/s in dF/ds
  dsdrho = -fourthird * s2 / temprho
  
  ! sometimes s=0, we do not divide s here, and we substract s from df/ds for
  ! every functional
  dsdgrad(:,:,:,1) = grad(:,:,:,1) / rho83
  dsdgrad(:,:,:,2) = grad(:,:,:,2) / rho83
  dsdgrad(:,:,:,3) = grad(:,:,:,3) / rho83
 

                               !>> Funtional Body <<!
  ! each dF/ds actually is dF/ds/s
  Select case(functionals)
    CASE(0)
      write(*,*) "Error! No GGA functional is specified"
      STOP
    
    CASE(1) ! TF only
      F = 1.d0
      dFds = 0.d0
    
    CASE(2) ! vW only
      F = s**2/CTF/8.d0
      dFds = 1.d0/CTF/4.d0
    
    CASE(3) ! aTF + bvW
      F = lambda + mu * s**2/CTF/8.d0
      dFds = mu/CTF/4.d0
    
    CASE(4) ! LLP functional ( #? in ref(2) )
      F = 1.d0 + LLPb * s**2 / ( 1.d0 + LLPc * s * asinh(2**(1.d0/3.d0)*s) )
      dFds = - LLPb*s*(2**(1.d0/3.d0)*LLPc*s/Sqrt(1.d0+2**(2.d0/3.d0)*s**2) + LLPc*ASinh(2**(1.d0/3.d0)*s))/&
               (1.d0+LLPc*s*ASinh(2**(1.d0/3.d0)*s))**2 &
             + 2.d0*LLPb/(1.d0+LLPc*s*ASinh(2**(1.d0/3.d0)*s))
    
    CASE(5) ! OL2 functional
      F = 1.d0 + s**2/72.d0/CTF + OL2D*s/(1.d0 + 2.d0**(5.d0/3.d0)*s)
      dFds = s/36.d0/CTF - 2.d0**(5.d0/3.d0)*OL2D*s/(1.d0+2**(5.d0/3.d0)*s)**2 + OL2D/(1.d0+2.d0**(5.d0/3.d0)*s)
      do i=1,size(dFds,1)
          do j=1,size(dFds,2)
              do k=1,size(dFds,3)
                  if(s(i,j,k)>1e-10) then
                      dFds(i,j,k) = dFds(i,j,k)/s(i,j,k)
                  endif
              enddo
          enddo
      enddo

    CASE(6) ! PW91k functional
      F = (1.d0 + (PWA2-PWA3*Exp(-PWA4*s**2))*s**2 + PWA1*s*ASinh(PWA*s)) / (1.d0 + PWB1*s**4 + PWA1*s*ASinh(PWA*s) )
      dfds = - ( 4.d0*PWB1*s**3+PWA*PWA1*s/Sqrt(1.d0+PWA**2*s**2) + PWA1*ASinh(PWA*s) ) * &
               ( 1.d0+(PWA2-PWA3*Exp(-PWA4*s**2))*s**2+PWA1*s*ASinh(PWA*s) ) / &
               ( 1.d0 + PWB1*s**4+PWA1*s*ASinh(PWA*s) )**2 &
             + ( 2.d0*(PWA2-PWA3*Exp(-PWA4*s**2))*s+2.d0*PWA3*PWA4*Exp(-PWA4*s**2)*s**3 &
             +PWA*PWA1*s/Sqrt(1.d0+PWA**2*s**2)+PWA1*ASinh(PWA*s) ) / &
               ( 1.d0+PWB1*s**4+PWA1*s*ASinh(PWA*s))
      do i=1,size(dFds,1)
          do j=1,size(dFds,2)
              do k=1,size(dFds,3)
                  if(s(i,j,k)>1e-10) then
                      dFds(i,j,k) = dFds(i,j,k)/s(i,j,k)
                  endif
              enddo
          enddo
      enddo
    
    CASE(7) ! TW02
      F = 1.d0 + TWk - TWk/(1.d0+TWmu/TWk*s2)
      dFds = 2.d0*TWmu/(1.d0+TWmu/TWk*s2)**2
    
    CASE(8) ! PBE2 by Trickey's
      F = 1.d0 + PBE2k - PBE2k/(1.d0+PBE2mu/PBE2k*s2)
      dFds = 2.d0*PBE2mu/(1.d0+PBE2mu/PBE2k*s2)**2
    
    CASE(9) ! LC94 functional
      F = (1.d0 + (LCA2-LCA3*Exp(-LCA4*s**2))*s**2 + LCA1*s*ASinh(LCA*s)) / (1.d0 + LCB1*s**4 + LCA1*s*ASinh(LCA*s) )
      dfds = - ( 4.d0*LCB1*s**3+LCA*LCA1*s/Sqrt(1.d0+LCA**2*s**2) + LCA1*ASinh(LCA*s) ) * &
               ( 1.d0+(LCA2-LCA3*Exp(-LCA4*s**2))*s**2+LCA1*s*ASinh(LCA*s) ) / &
               ( 1.d0 + LCB1*s**4+LCA1*s*ASinh(LCA*s) )**2 &
             + ( 2.d0*(LCA2-LCA3*Exp(-LCA4*s**2))*s+2.d0*LCA3*LCA4*Exp(-LCA4*s**2)*s**3 &
             +LCA*LCA1*s/Sqrt(1.d0+LCA**2*s**2)+LCA1*ASinh(LCA*s) ) / &
               ( 1.d0+LCB1*s**4+LCA1*s*ASinh(LCA*s))
      do i=1,size(dFds,1)
          do j=1,size(dFds,2)
              do k=1,size(dFds,3)
                  if(s(i,j,k)>1e-10) then
                      dFds(i,j,k) = dFds(i,j,k)/s(i,j,k)
                  endif
              enddo
          enddo
      enddo

    Case(10)   ! E00 functional
      F = (E00A1 + E00A2*s**2 + E00A3*s**4) / (E00A1 + E00A4*s**2)
      dfds = (2.d0*E00A2+4.d0*E00A3*s**2)/(E00A1 + E00A4*s**2) - &
             2.d0*E00A4*(E00A1+E00A2*s**2+E00A3*s**4) / (E00A1 + E00A4*s**2)**2
    
    Case(11)   ! P92 functional
      F = (P92A1 + P92A2*s**2 + P92A3*s**4) / (P92A1 + P92A4*s**2)
      dfds = (2.d0*P92A2+4.d0*P92A3*s**2)/(P92A1 + P92A4*s**2) - &
             2.d0*P92A4*(P92A1+P92A2*s**2+P92A3*s**4) / (P92A1 + P92A4*s**2)**2

    Case(12)   ! PW86 functional
      F = (1.d0 + PW86A1*s**2 + PW86A2*s**4 + PW86A3*s**6) ** (1.d0/15.d0)
      dfds = (2.d0*PW86A1+4.d0*PW86A2*s**2+6.d0*PW86A3*s**4) &
           / 15.d0 / (1.d0+PW86A1*s**2+PW86A2*s**4+PW86A3*s**6)**(14.d0/15.d0)
    
    Case(13)   ! DK functional
      x = s**2/72.d0/CTF
      F = (9.d0*DKB3*x**4+DKA3*x**3+DKA2*x**2+DKA1*x+1.d0) / &
          (DKB3*x**3+DKB2*x**2+DKB1*x+1.d0)
      dFds = (DKA1+2.d0*DKA2*x+3.d0*DKA3*x**2+36.d0*DKB3*x**3) / &
             (1.d0+DKB1*x+DKB2*x**2+DKB3*x**3) - &
             (DKB1+2.d0*DKB2*x+3.d0*DKB3*x**2)*(1.d0+DKA1*x+DKA2*x**2+DKA3*x**3+9.d0*DKB3*x**4) / &
             (1.d0+DKB1*x+DKB2*x**2+DKB3*x**3)**2
      dFds = dfds/36.d0/CTF

    Case(14)   ! OL1 functional
      F = 1.d0 + s**2/72.d0/CTF + OL1d*s
      dFds = s/36.d0/CTF + OL1d
      do i=1,size(dFds,1)
          do j=1,size(dFds,2)
              do k=1,size(dFds,3)
                  if(s(i,j,k)>1e-10) then
                      dFds(i,j,k) = dFds(i,j,k)/s(i,j,k)
                  endif
              enddo
          enddo
      enddo
   
    Case(15)   ! B86A functional
      F = 1.d0 + B86AC1*s**2/(1.d0+B86AC2*s**2)
      dFds = -2.d0*B86AC1*B86AC2*s**2/(1.d0+B86AC2*s**2)**2 + &
              2.d0*B86AC1/(1.d0+B86AC2*s**2)
    
    Case(16)   ! B86B functional
      F = 1.d0 + B86BC1*s**2/(1.d0+B86BC2*s**2)**(4.d0/5.d0)
      dFds = - 8.d0*B86BC1*B86BC2*s**2/5.d0/(1.d0+B86BC2*s**2)**(9.d0/5.d0) + &
               2.d0*B86BC1/(1.d0+B86BC2*s**2)**(4.d0/5.d0)

    Case(17)   ! Thak functional
      F = 1.d0 + Thakb * s**2 / (1.d0 + Thakc * s * asinh(2**(1.d0/3.d0)*s)) - ThakD*s/(1.d0+2**(5.d0/3.d0)*s)
      dFds = - Thakb*s**2*(2**(1.d0/3.d0)*Thakc*s/Sqrt(1.d0+2**(2.d0/3.d0)*s**2) + Thakc*ASinh(2**(1.d0/3.d0)*s))/&
               (1.d0+Thakc*s*ASinh(2**(1.d0/3.d0)*s))**2 &
             + 2.d0*Thakb*s/(1.d0+Thakc*s*ASinh(2**(1.d0/3.d0)*s)) +  2.d0**(5.d0/3.d0)*ThakD*s/(1.d0+2**(5.d0/3.d0)*s)**2 &
             - ThakD/(1.d0+2**(5.d0/3.d0)*s)
      do i=1,size(dFds,1)
          do j=1,size(dFds,2)
              do k=1,size(dFds,3)
                  if(s(i,j,k)>1e-10) then
                      dFds(i,j,k) = dFds(i,j,k)/s(i,j,k)
                  endif
              enddo
          enddo
      enddo

    Case(18)   ! DK87
      F = 1.d0 + DK87A*s**2*(1.d0+DK87B*s)/(1.d0+DK87C*s**2)
      dFds = - 2.d0*DK87A*DK87C*s**2*(1.d0+DK87B*s)/(1.d0+DK87C*s**2)**2 + &
               DK87A*DK87B*s/(1.d0+DK87C*s**2) + 2.d0*DK87A*(1.d0+DK87B*s)/(1.d0+DK87C*s**2)
    
    Case(19)   ! LKT
        ! code contributed by Kai Luo (University of Florida)
        ! reference:
        ! K. Luo, V.V. Karasiev, S.B. Trickey, Phys. Rev. B 98 (2018): 041111.

         where (s<1e-5)
           ! (10/3 - a^2) + (5 a^4 s^2)/6 - (61 a^6 s^4)/120
           F=1.0+(5.d0/3.d0/(cs**2)-LKTa**2/2.0)*s2+5.d0*LKTa**4*s2**2/24.d0
           dFds = 10.d0/3.d0/(cs**2) - LKTa**2 + 5*LKTa**4*(s2)/6.d0 - 61*LKTa**6*(s2)**2/120.d0
         elsewhere
           F=5.d0/3.d0*(s2/cs**2) + 1.d0/cosh(LKTa*sqrt(s2))
           dFds=10.d0/3.d0/(cs**2) - LKTa/cosh(LKTa*sqrt(s2))*tanh(LKTa*sqrt(s2))/sqrt(s2)
         endwhere

    Case(-99)   ! return zero
      F = 0.d0
      dFds = 0.d0

  End Select
   
  ! set CP accorading the vW coefficient, because the way we do vW here is not stable
  if(functionals==2) then
    CP = 1.d0
  elseif(functionals==3) then
    CP = mu
  else
    CP = 0.d0
  endif

  ! calculate energy if needed
  If(calcEnergy) then
    energy = Sum(tf*F)
    PenaltyE = CP * ( 0.5_DP * SUM(sqrt(temprho) * FFT( FFT( sqrt(temprho) ) * qTable**2)) - &
                    SUM( (grad(:,:,:,1)**2 + grad(:,:,:,2)**2 + grad(:,:,:,3)**2)/temprho )/8.d0 )
    energy = energy + PenaltyE
  Endif
  
  ! calculate potential
  ! the term in the divergence calculation
  group(:,:,:,1) = tf * dfds * dsdgrad(:,:,:,1)
  group(:,:,:,2) = tf * dfds * dsdgrad(:,:,:,2)
  group(:,:,:,3) = tf * dfds * dsdgrad(:,:,:,3)
  
  ! get divergence term first
  potential = FFT( imag*qVectors(:,:,:,1)*FFT( group(:,:,:,1) ) + &
                   imag*qVectors(:,:,:,2)*FFT( group(:,:,:,2) ) + &
                   imag*qVectors(:,:,:,3)*FFT( group(:,:,:,3) ) )

  ! final potential
  potential = dtfdrho * F + tf * dFds * dsdrho - potential

  PenaltyP = FFT( imag*qVectors(:,:,:,1)*FFT( grad(:,:,:,1)/temprho/4.d0 )  + &
                  imag*qVectors(:,:,:,2)*FFT( grad(:,:,:,2)/temprho/4.d0 )  + &
                  imag*qVectors(:,:,:,3)*FFT( grad(:,:,:,3)/temprho/4.d0 )  )
  PenaltyP = -( grad(:,:,:,1)**2.d0 + grad(:,:,:,2)**2.d0 +  grad(:,:,:,3)**2.d0 ) / temprho**2.d0 / 8.d0 - PenaltyP

  PenaltyP = CP * ( FFT( FFT(sqrt(temprho)) * qTable**2)/2.d0/sqrt(temprho) - PenaltyP )
  
  potential = potential + PenaltyP

END Subroutine GGAPotentialPlus


SUBROUTINE vWGTF(rho, potential, calcEnergy, energy, model) 
!------------------------------------------------------------------------------
! This is a subroutine for vW + GTF, where G is a functional of rho
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   J. Xia and E. A. Carter, submitted to PRB
!------------------------------------------------------------------------------
! REVISION LOG:
!   Sept/05/2013   Created (Steven)
!------------------------------------------------------------------------------
  USE MathSplines, ONLY: spline_cubic_set
  USE MathSplines, ONLY: spline_cubic_val
  USE Sys, ONLY: rho0

  IMPLICIT NONE

                                !>> External variables <<!
  LOGICAL, INTENT(in) :: &
    calcEnergy          ! if energy is calculated
  
  Integer, INTENT(in) :: &
    model

  REAL(kind=DP), Dimension(:,:,:), INTENT(in) :: &
    rho                 ! total electron density
  
  REAL(kind=DP), INTENT(out) :: &
    energy              ! kinetic energy

  REAL(kind=DP), Dimension(:,:,:), INTENT(out) :: &
    potential           ! dE/drho, i.e., kinetic energy potential

                                !>> Internal variables <<!
  Integer :: &
    i, j, k
  
  Real(DP), parameter :: &
    cTF = 2.87123400018819_DP, &
    twothird = 2.d0/3.d0, &
    fivethird = 5.d0/3.d0, &
    ! model 1 parameters:
    expo = -1.29941357_DP, &
    co = exp(-0.01088208_DP), &
    ! model 2 parameters:
    ga = 5.70010115157, &
    gb = 0.256298586338

  Real(DP) :: &
    yval, ypval, yppval ! for model 3 interpolation use
  
  REAL(DP), Dimension(size(rho,1),size(rho,2),size(rho,3)) :: &
    rho53, rho23, sqrtrho, Drhob, &
    DlessRho, &  ! rho/rho0
    G, ELF, &         ! the enhancemence factor for TF term
    dGdrho
  
!>>>>>>>>>>>>>>>>>>  Initialization <<<<<<<<<<<<<<<<<<<<<<<
    energy = 0.d0
    potential = 0.d0
    rho23 = rho**twothird
    rho53 = rho*rho23
    sqrtrho = sqrt(rho)
    DlessRho = rho/rho0
    Drhob = DlessRho**gb

!>>>>>>>>>>>>>>>>>>  functional body <<<<<<<<<<<<<<<<<<<<<<
    ! first vW term
    potential = FFT( FFT(sqrtrho) * qTable**2)
    IF (calcEnergy) energy = 0.5_DP * SUM(sqrtRho*potential)
    potential = potential/(2._DP*sqrtrho)

    ! G term
    Select case(model)
    Case(1)
        G = co * DlessRho**expo
        dGdrho = co*expo * DlessRho**(expo-1.d0) / rho0
    Case(2)
        ELF = (tanh(ga*DRhob - ga) + 1.d0)/2.d0
        G = sqrt(1.d0/ELF - 1.d0)
        dGdrho = -0.5d0/ELF/ELF/G
        dGdrho = dGdrho*0.5d0*ga*gb*Drhob/DlessRho/cosh(ga-ga*DRhob)**2
        dGdrho = dGdrho/rho0
    Case(3)  ! interpolate ELF vs d
        do i = lbound(G,1),ubound(G,1)
            do j = lbound(G,2),ubound(G,2)
                do k = lbound(G,3),ubound(G,3)
                    if(DlessRho(i,j,k)<2.0 .AND. DlessRho(i,j,k)>0.0) then
                        call spline_cubic_val(numint, d(1:numint), ELFvsd(1:numint),ELFvsdDD(1:numint), &
                                              DlessRho(i,j,k), yval, ypval, yppval )
                        ELF(i,j,k) = yval
                        dGdrho(i,j,k) = ypval/rho0  ! actually it is dELF/drho here
                    ELSE
                        ELF(i,j,k) = (tanh(ga*DRhob(i,j,k) - ga) + 1.d0)/2.d0
                        dGdrho(i,j,k) = 0.5d0*ga*gb*Drhob(i,j,k)/DlessRho(i,j,k)/ &
                                        cosh(ga-ga*DRhob(i,j,k))**2
                        dGdrho(i,j,k) = dGdrho(i,j,k)/rho0
                    Endif
                enddo
            enddo
        enddo
        G = sqrt(1.d0/ELF - 1.d0)
        dGdrho = (-0.5d0/ELF/ELF/G)*dGdrho
    End select
    
    If(CalcEnergy) then
        energy = energy + Sum(CTF*rho53*G)
    Endif

    potential = potential + CTF*fivethird*rho23*G + CTF*rho53*dGdrho

End Subroutine vWGTF


END MODULE KEDF_GGA













 



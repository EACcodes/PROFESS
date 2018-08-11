MODULE KEDF_HC10 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_HC10 
!     |_FUNCTION intStress 
!     |_FUNCTION intEnergy
!     |_FUNCTION intPot 
!     |_SUBROUTINE MakeXi
!     |_SUBROUTINE MakeHC10Kernel
!     |_FUNCTION H10
!     |_FUNCTION H01
!     |_FUNCTION H11
!     |_FUNCTION H01
!
! DESCRIPTION:
!   Calculates the kinetic forces, stress, and energies.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
! "Nonlocal orbital-free kinetic energy density functional for semiconductors"
! Phys. Rev. B 81, 045206 (2010) Chen Huang and Emily A. Carter
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Written by Chen Huang.
!   2013-08-05  Recreated by Mohan Chen
!
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!
  USE CONSTANTS, ONLY: DP
  USE CONSTANTS, ONLY: PI
  USE CONSTANTS, ONLY: IMAG
  USE CellInfo, ONLY: k1G, k2G, k3G, k123G ! recip space array dimension (for each proc.)
  USE CellInfo, ONLY: n1G, n2G, n3G, n123G ! real space array dimension (for each proc.)
  USE Sys, ONLY: rhoS ! useful for binning
  USE KEDF_WTkernel, ONLY: alpha, beta ! alpha and beta for non-local KEDF
  USE PlaneWave, ONLY: qTable ! plane wave norm
  USE PlaneWave, ONLY: qMask ! be careful about double counting
  USE PlaneWave, ONLY: qVectors ! plane wave vector
  USE FOURIER, ONLY : FFT  ! The Fast Fourier Transform routine.
  USE MPI_Functions
  USE OUTPUT, ONLY: WrtOut
  USE KEDF_TF, ONLY: cTF

  IMPLICIT NONE

  INTEGER  :: nref = 0       
  ! number of bins used for interpolation 
  !
  REAL(KIND=DP) :: cutrho=0.d0
  !
  INTEGER :: interp_method = 1
  ! 1: algrithmic sequence, 2: geometric sequence
  !
  INTEGER :: nband  = -1
  ! How many intervals used in the spline?
  ! -1: means the code will increase the number of intervals 
  ! depending on xi automatically
  !
  REAL(KIND=DP) :: hc_lambda_val = 0.01_DP
  ! the value will be set to hc_lambda arrary
  !
  REAL(KIND=DP) :: refRatio = 1.15_DP
  ! for set bins
  !
  REAL(KIND=DP), ALLOCATABLE :: refxi(:)  
  ! a set of xi values for interpolation
  !
  REAL(KIND=DP), ALLOCATABLE :: hc_lambda(:,:,:)  
  ! lambda array for Huang-Carter KEDF
  !
  LOGICAL :: have_RK = .false. 
  ! flag for if the Runge-Kutta method has been called
  !
  REAL(KIND=DP), EXTERNAL :: DDOT
  ! if we use blas
  !

CONTAINS


FUNCTION intStress(vol, rho)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! Here's the function that computes the HC KEDF portion of the stress and its
! contribution to the total stress tensor. It is calculated a lot less often
! than the energy or potential, so one routine is enough. Either way, the
! computation is a pain anyways.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   09/06/2011  Function created.  (Steven Xia)
!
!------------------------------------------------------------------------------

  USE MathFunctions, ONLY: GPrime
  ! Returns the first derivative of G.
  !
  USE MathFunctions, ONLY: stepfun
  !
  USE CellInfo, ONLY: cell 
  !
  USE PlaneWave, ONLY: qVectors
  !

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(3,3) :: intStress
  ! The answer we want
  !
  REAL(KIND=DP), INTENT(IN) :: vol
  ! Volume of the cell in real space.
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho   
  ! Spin-independant density in real space.
  !

                         !>> INTERNAL VARIABLES <<!
  INTEGER :: i, j, k, ii, jj, a, b   
  ! counters
  !
  INTEGER :: upperpower, downpower, allocatestat
  !
  REAL(KIND=DP) :: C = 8.d0*3.d0*pi**2*cTF
  !
  REAL(KIND=DP) :: maxxi, minxi, midxi, h, tempSum
  !
  REAL(KIND=DP), ALLOCATABLE :: refxi(:)
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: rhoA
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: rhoB
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: kF
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: xi
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: ss
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: mask
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: t
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: t2
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: h00p
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: h01p
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: h10p
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: h11p
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: temp
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: temp2
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G,3) :: gradrho
  ! gradient of rho
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FFTb
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: kernel0, kernel1
  ! w(g) and F{dw/dxi}*xi, for xi(i)
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: kernel0a, kernel1a
  ! w(g) and F{dw/dxi}*xi, for xi(i+1)
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: dwdeta, d2wdeta2
  ! those are dw(eta)/deta and d2w(eta)/deta2, for xi(i)
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: dwdetaa, d2wdeta2a     
  ! those are dw(eta)/deta and d2w(eta)/deta2, for xi(i+1)
  !
  REAL(KIND=DP), DIMENSION(k1G, k2G, k3G) :: eta, etaa, Flind, Flinda  
  ! dimensionless qvector and its corresponding Lind Function
  !

  !! >> INITIALIZATION << !!
  intStress = 0._DP

  ! determine dV for real space integration
  RhoA = rho**alpha
  RhoB = rho**beta
  gradRho(:,:,:,1) = FFT(imag*qVectors(:,:,:,1)*FFT(rho))
  gradRho(:,:,:,2) = FFT(imag*qVectors(:,:,:,2)*FFT(rho))
  gradRho(:,:,:,3) = FFT(imag*qVectors(:,:,:,3)*FFT(rho))
  FFTb = FFT(RhoB)

  ! first calculate kf, xi and ss according 
  kF = (3.d0*pi**2*rho) ** 0.333333333333_dp
  CALL MakeXi(rho,xi,ss)

  ! Now set the refxi array in order for intepolation
  maxxi=maxval(xi)
  minxi=minval(xi)

#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(maxxi,maxxi,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpiErr)
  CALL MPI_ALLREDUCE(minxi,minxi,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpiErr)
#endif

  ! make refxi array
  !==================
  midXi = (3.d0*pi**2.d0*rhoS)**(1.d0/3.d0) 
  upperPower = ceiling(LOG(maxXi/midXi) / LOG(refRatio))
  downPower  = floor  (LOG(minXi/midXi) / LOG(refRatio))
  nref = upperpower-downpower+1
  ALLOCATE (refxi(nref), stat = allocatestat)
  IF(allocatestat /= 0) then
    WRITE(*,*) 'error in allocating rhobnd() in dividedensity (), stop'
    STOP
  ENDIF

  DO ii = nref,1,-1
    refxi(ii) = midxi * refRatio**(downPower + ii-1)
  ENDDO
  ! setting refxi array done
  
  ! loop over bins
  !===============
  DO jj =1,nref-1 
    ! get w(g) and F{dw(xi*deltar)/dxi}*xi, at refxi(jj) and refxi(jj+1)
    CALL MakeHC10Kernel(refxi(jj),kernel0,kernel1)
    CALL MakeHC10Kernel(refxi(jj+1),kernel0a,kernel1a)
    
    ! calculate eta and Flind for both refxi(jj) and refxi(jj+1), because the refxi is changing
    ! the eta and etaa needs to updated during each bin calculation
    eta = qTable/2.d0/refxi(jj)
    ! zero eta would be problematic, we do the trick the same as in the int_ode subroutine
    WHERE (eta==0.d0)
      eta = 1e-3
    ENDWHERE

    WHERE (eta==1.d0) 
      FLind = 2.d0
    ELSEWHERE (eta==0.0d0) 
      FLind = 1.d0
    ELSEWHERE
      FLind = 1.d0/(1.d0/2.d0+(1.d0-eta**2)/(4.d0*eta)*LOG(ABS((1.d0+eta)/(1.d0-eta))))
    ENDWHERE
    
    etaa = qTable/2.d0/refxi(jj+1)
    where (etaa==0.d0)
      etaa = 1e-3
    endwhere
    WHERE (etaa==1.d0) 
      FLinda = 2.d0
    ELSEWHERE (etaa==0.0d0) 
      FLinda = 1.d0
    ELSEWHERE
      FLinda = 1.d0/(1.d0/2.d0+(1.d0-etaa**2)/(4.d0*etaa)*LOG(ABS((1.d0+etaa)/(1.d0-etaa))))
    ENDWHERE

    ! calculate h, t, mask, h10p and others  
    h = refXi(jj+1) - refXi(jj)
    t = (xi - refXi(jj)) / h
    mask = stepfun(xi-refXi(jj)) * stepfun(refXi(jj+1)-xi)
    t2 = t**2
    h00p = 6.d0*t2 - 6.d0*t
    h10p = 3.d0*t2 - 4.d0*t + 1.d0
    h01p = -6.d0*t2 + 6.d0*t
    h11p = 3.d0*t2 - 2.d0*t

    ! Loop over stress matrix
    DO a = 1, 3
      DO b = a, 3
         ! first take care of the diag part which involves delta(ab) term
         if(a==b) then
            ! first part, actually is proportional to the HC energy
            temp = h00(t) * FFT(kernel0(:,:,:)*FFTb) + &
                   h10(t) * h * FFT(kernel1(:,:,:)/refxi(jj)*FFTb) + &
                   h01(t) * FFT(kernel0a(:,:,:)*FFTb) + &
                   h11(t) * h * FFT(kernel1a(:,:,:)/refxi(jj+1)*FFTb)
            
            tempSum = SUM(RhoA*mask*temp)
            intStress(a,b) = intStress(a,b) + C * (1.d0-alpha-beta) * tempSum/vol*cell%dV

            ! then the second part
            temp = 1.d0/h * h00p * FFT(kernel0(:,:,:)*FFTb) + &
                   h10p * FFT(kernel1(:,:,:)/refxi(jj)*FFTb) + &
                   1.d0/h * h01p * FFT(kernel0a(:,:,:)*FFTb) + &
                   h11p * FFT(kernel1a(:,:,:)/refXi(jj+1)*FFTb)
            tempSum = SUM(rhoA*mask*1.d0/3.d0*kF*(hc_lambda*ss-1.d0)*temp)
            intStress(a,b) = intStress(a,b) + C * tempSum/vol*cell%dV
         endif  ! the delta(ab) term added

         ! Now the third term in the equation, we need set dwdeta, and d2wdeta2
         dwdeta = 5.d0/3.d0*(Flind-3.d0*eta*eta-1) - (5.d0-3.d0*beta)*beta*kernel0*8.d0*refxi(jj)**3
         dwdeta = -dwdeta/beta/eta
         dwdetaa = 5.d0/3.d0*(Flinda-3.d0*etaa*etaa-1) - (5.d0-3.d0*beta)*beta*kernel0a*8.d0*refxi(jj+1)**3
         dwdetaa = -dwdetaa/beta/etaa
         do i=lbound(d2wdeta2,1),ubound(d2wdeta2,1)
            do j=lbound(d2wdeta2,2),ubound(d2wdeta2,2)
               do k=lbound(d2wdeta2,3),ubound(d2wdeta2,3)
                  ! this is derived from the kernel equation
                  d2wdeta2(i,j,k) = 5.d0/3.d0*GPrime(eta(i,j,k),1.d0) + beta*dwdeta(i,j,k) - &
                                    (5.d0-3.d0*beta)*beta*dwdeta(i,j,k)
                  d2wdeta2(i,j,k) = -d2wdeta2(i,j,k)/beta/eta(i,j,k)
                  d2wdeta2a(i,j,k) = 5.d0/3.d0*GPrime(etaa(i,j,k),1.d0) + beta*dwdetaa(i,j,k) - &
                                    (5.d0-3.d0*beta)*beta*dwdetaa(i,j,k)
                  d2wdeta2a(i,j,k) = -d2wdeta2a(i,j,k)/beta/etaa(i,j,k)
               enddo
            enddo
         enddo
         ! dwdeta and d2wdeta2 done
         ! calculate the third term now:
         temp2 = h00(t) * FFT(-FFTb*dwdeta/16.d0/(refxi(jj)**4)*qVectors(:,:,:,a)*qVectors(:,:,:,b)/(qTable+1e-18)) + &
                 h * h10(t) * FFT(-FFTb*(1.d0/16.d0/refxi(jj)**5*(-4.d0*dwdeta-eta*d2wdeta2))* &
                                                 qVectors(:,:,:,a)*qVectors(:,:,:,b)/(qTable+1e-18)) + &
                 h01(t) * FFT(-FFTb*dwdetaa/16.d0/(refxi(jj+1)**4)*qVectors(:,:,:,a)*qVectors(:,:,:,b)/(qTable+1e-18)) + &
                 h * h11(t) * FFT(-FFTb*(1.d0/16.d0/refxi(jj+1)**5*(-4.d0*dwdetaa-eta*d2wdeta2a))* &
                                                    qVectors(:,:,:,a)*qVectors(:,:,:,b)/(qTable+1e-18))

         tempSum = SUM(rhoA*mask*temp2)
         intStress(a,b) = intStress(a,b) + C * tempSum/vol*cell%dV

         ! The final term
         if(a/=b) then  ! the temp term has not already been calculated 
            temp = 1.d0/h * h00p * FFT(kernel0(:,:,:)*FFTb) + &
                   h10p * FFT(kernel1(:,:,:)/refxi(jj)*FFTb) + &
                   1.d0/h * h01p * FFT(kernel0a(:,:,:)*FFTb) + &
                   h11p * FFT(kernel1a(:,:,:)/refXi(jj+1)*FFTb)
         endif
         
         tempSum = SUM(rhoA*mask*(3.d0*pi**2)**0.333333333333_dp*hc_lambda/rho**(7.d0/3.d0)* &
                       (-2.d0*gradrho(:,:,:,a)*gradrho(:,:,:,b)) * temp)
         intStress(a,b) = intStress(a,b) + C * tempSum/vol*cell%dV

      END DO ! b
    END DO ! a
  END DO ! jj

  RETURN
  
END FUNCTION intStress


FUNCTION intEnergy(rho) 
!------------------------------------------------------------------------------
! DESCRIPTION:
! Total energy of this single density dependent kernel is computed first by
! splining the kernel in real space. see details on Wiki
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Nov/10/2009   Created (Chen Huang)
!------------------------------------------------------------------------------
  USE MathFunctions, ONLY : stepfun

  IMPLICIT NONE

  !>>> EXTERNAL VARS <<<<!
  REAL(KIND=DP) :: rho(:,:,:)
  !
  REAL(KIND=DP) :: intEnergy
  !
  
  !>>> INTERNAL VARS <<<<!
  INTEGER :: ii, jj
  !
  INTEGER :: upperpower, downpower
  !
  INTEGER :: allocatestat
  !
  REAL(KIND=DP) :: c = 8.d0*3.d0*pi**2
  !
  REAL(KIND=DP) :: h
  !
  REAL(KIND=DP) :: maxxi, minxi, midxi
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: rhoA
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: xi
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: ss
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: rhoMask
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: t
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: t2
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: t3
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: h00
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: h01
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: h10
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: h11

  REAL(KIND=DP), ALLOCATABLE :: refxi(:)

  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FFTx
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: kernel0,kernel1
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: kernel0a, kernel1a  
  ! Fourier transform of rho for temp stroage.
  !
 
  CALL Title("KEDF_HC10::intEnergy")
  CALL StartClock('intEnergy')

  CALL MakeXi(rho,xi,ss)
  
  maxxi=maxval(xi)
  minxi=minval(xi)

#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(maxxi,maxxi,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpiErr)
  CALL MPI_ALLREDUCE(minxi,minxi,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpiErr)
#endif

  ! make refxi array
  !==================
  midXi = (3.d0*pi**2.d0*rhoS)**(1.d0/3.d0) 
  upperPower = CEILING(LOG(maxXi/midXi) / LOG(refRatio))
  downPower  = FLOOR  (LOG(minXi/midXi) / LOG(refRatio))
  nref = upperpower-downpower+1

  ALLOCATE (refxi(nref), stat = allocatestat)
  IF (allocatestat /= 0) THEN
    WRITE(*,*) 'error in allocating rhobnd() in dividedensity (), stop'
    STOP
  ENDIF

  DO ii = nref,1,-1
    refxi(ii) = midxi * refRatio**(downPower + ii-1)
  ENDDO

!  CALL MPI_BARRIER(MPI_COMM_WORLD,err)
  
  rhoA = rho**alpha
  FFTx = 0.d0

  ! loop over bins
  !===============
  ! WRITE(*,*) "nref in intEnergy is ", nref
  DO jj =1,nref-1 

    CALL StartClock('hc1')

    CALL MakeHC10Kernel(refxi(jj),kernel0,kernel1)
    CALL MakeHC10Kernel(refxi(jj+1),kernel0a,kernel1a)

    CALL StopClock('hc1')
    CALL StartClock('hc2')

    h = refXi(jj+1) - refXi(jj)
    t = (xi - refXi(jj)) / h

!new ones
    rhoMask = rhoA
    IF(cutrho==0.d0)THEN
      WHERE( (xi .LT. refxi(jj)) .OR. (xi .GE. refxi(jj+1)))
        rhoMask = 0.d0
      ENDWHERE
    ELSE
      WHERE( (xi .LT. refxi(jj)) .OR. (xi .GE. refxi(jj+1)) .OR. (rho .LT. cutrho))
        rhoMask = 0.d0
      ENDWHERE
    ENDIF

!old ones
!    IF (cutrho == 0.d0) THEN   
!      mask = stepfun(xi-refXi(jj)) * stepfun(refXi(jj+1)-xi)
!    ELSE
!      mask = stepfun(rho-cutrho) * stepfun(xi-refXi(jj)) * stepfun(refXi(jj+1)-xi)
!    ENDIF
!    rhoMask = rhoA*mask

    t2 = t * t
    t3 = t2 * t

    h00 = 2.d0*t3-3.d0*t2+1.d0
    h10 = t3 - 2.d0*t2 + t
    h01 = 1.d0 - h00
    h11 = t3 - t2

    FFTx = FFTx &
       + kernel0  * FFT( rhoMask*h00 ) &
       + kernel0a * FFT( rhoMask*h01 )                  &
       + h * ( kernel1 /refXi(jj) * FFT( rhoMask*h10 )  &
       + kernel1a /refXi(jj+1) * FFT( rhoMask*h11 )  )

    CALL StopClock('hc2')
  ENDDO 
   
   intEnergy = cTF*C*SUM(rho**beta * FFT(FFTx))
!   WRITE(*,*) "intEnergy=", intEnergy

!  intEnergy = cTF*C*SUM(rho**beta * stepfun(rho-cutrho) * FFT(FFTx)) !this is for consistent vaccum cut
!  CALL WrtOut(6,'(exit) intEnergy ')
!  PRINT *,'leave intEnergy. intEnergy =', intEnergy
!  STOP

  deallocate (refxi)
  
  if(maxval(abs(rho))<1e-20) then
      write(*,*) 'no this spin channel electron'
      intEnergy=0.d0
  endif

  CALL StopClock('intEnergy')

END FUNCTION intEnergy


FUNCTION intPot(rho, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Total energy of this single density dependent kernel is computed first by
! splining the kernel in real space. This is the subroutine computing 
! the potenital associated with it. See details on Wiki.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Nov/10/2009   Created (Chen Huang)
!------------------------------------------------------------------------------
  USE MATHFUNCTIONS, ONLY : stepfun 

  IMPLICIT NONE

  !>>> EXTERNAL VARS <<<<!
  REAL(KIND=DP) :: rho(n1G,n2G,n3G)
  REAL(KIND=DP) :: intpot(n1G,n2G,n3G)
  REAL(KIND=DP) :: energy 
  
  !>>> INTERNAL VARS <<<<!
  INTEGER :: ii, jj
  INTEGER :: upperpower
  INTEGER :: downpower
  INTEGER :: allocatestat

  REAL(KIND=DP) :: c = 8.d0*3.d0*pi**2
  REAL(KIND=DP) :: h = 0.d0
  REAL(KIND=DP) :: maxxi = 0.d0    ! small number
  REAL(KIND=DP) :: minxi = 1000.d0 ! random large number
  REAL(KIND=DP) :: midxi

  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: rhoA
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: rho83
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: xi
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: ss
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: mask
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: kF
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: intpot1
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: potTmp
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: pxiprho

  ! mohan add
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: fkf0
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: fkf0_saved
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: fk1rf
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: fk1rf_saved

  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: fkfa
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: fkarf

  ! in G space
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kref1
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kref1_saved
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kref2

  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: t
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: t2
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: t3

  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: h00
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: h01
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: h10
  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: h11

  REAL(KIND=DP), DIMENSION(n1G,n2G,n3G) :: rhoMask


  REAL(KIND=DP), ALLOCATABLE :: refxi(:)  ! nodes for interpolation
  REAL(KIND=DP) :: gradRho(n1G, n2G, n3G, 3)

  ! in G space
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kernel0
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kernel1
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kernel0a
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kernel1a

  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kernel0_saved
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: kernel1_saved

  ! both in G space
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: FFTb
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: FFTx

  REAL(KIND=DP) :: factor

  CALL Title("KEDF_HC10::intPot")

  !>>>>>>>>>>>>> FUNCTION BEGINS <<<<<<<<<<<!
 
  CALL StartClock('intPot')

  ! grad(rho)
  FFTb = FFT(rho)
  gradRho(:,:,:,1) = FFT(imag*qVectors(:,:,:,1)*FFTb)
  gradRho(:,:,:,2) = FFT(imag*qVectors(:,:,:,2)*FFTb)
  gradRho(:,:,:,3) = FFT(imag*qVectors(:,:,:,3)*FFTb)

  !CALL MakeXi(rho,xi,ss)
  rho83 = rho**(8.d0/3.d0)
  CALL MakeXiG(rho,xi,ss,gradRho,rho83)
  
  maxxi=maxval(xi)
  minxi=minval(xi)

#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(maxxi,maxxi,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpiErr)
  CALL MPI_ALLREDUCE(minxi,minxi,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpiErr)
#endif
 
  !WRITE(outputUnit,*) "Maximal and minimal of Xi and number of bin: ", maxxi, minxi, nref

  !------------------
  ! make refxi array
  !------------------
  midXi = (3.d0*pi**2.d0*rhoS)**(1.d0/3.d0) 
  upperPower = ceiling(LOG(maxXi/midXi) / LOG(refRatio))
  downPower  = floor  (LOG(minXi/midXi) / LOG(refRatio))
  nref = upperpower-downpower+1

  ALLOCATE (refxi(nref), stat = allocatestat)
  IF (allocatestat /= 0) then
    WRITE(*,*) 'error in allocating rhobnd() in dividedensity (), stop'
    STOP
  ENDIF

  DO ii = nref,1,-1
    refxi(ii) = midxi * refRatio**(downPower + ii-1)
  ENDDO

  kF = (3.d0*pi**2*rho)**(1.d0/3.d0) 
     
  intpot1 = 0.d0
  rhoA = rho**alpha 
  FFTb = FFT(rho**beta)
  FFTx = 0.d0
  potTmp = 0.d0
    
  !-----------------
  ! loop over bins
  !-----------------
  DO jj = 1,nref-1
    Call StartClock('intPot1')

    !---------------------------------------------------------------
    ! mohan found in fact each time we only need
    ! to calculate two new kernels in G space
    IF(jj > 1) THEN
        
        CALL MakeHC10Kernel(refXi(jj+1),kernel0a,kernel1a)


! mohan update, using BLAS routines
! be careful the dimensions of kernel is k123,
! and they are all complex values.
        CALL DCOPY(n123G, fkf0_saved, 1, fkf0, 1)
        CALL ZCOPY(k123G, kref1_saved, 1, kref1, 1)
        CALL DCOPY(n123G, fk1rf_saved, 1, fk1rf, 1)
        CALL ZCOPY(k123G, kernel0_saved, 1, kernel0, 1)
        CALL ZCOPY(k123G, kernel1_saved, 1, kernel1, 1)

! old ones, has been replaced by blas
!        fkf0 = fkf0_saved
!        kref1 = kref1_saved
!        fk1rf = fk1rf_saved
!        kernel0 = kernel0_saved
!        kernel1 = kernel1_saved
         
        ! new ones
        fkfa  = FFT(kernel0a(:,:,:) * FFTb)
        kref2 = kernel1a/refXi(jj+1)
        fkarf = FFT(kref2 * FFTb)

    ELSE IF(jj == 1) THEN

        CALL MakeHC10Kernel(refXi(jj),kernel0,kernel1)
        CALL MakeHC10Kernel(refXi(jj+1),kernel0a,kernel1a)

        ! G space to real space
        fkf0  = FFT(kernel0(:,:,:) * FFTb) ! don't keep this for next iteration
        kref1 = kernel1/refXi(jj)          ! don't keep this for next iteration
        fk1rf = FFT(kref1 * FFTb)          ! don't keep this for next iteration

        fkfa  = FFT(kernel0a(:,:,:) * FFTb)  ! keep this
        kref2 = kernel1a/refXi(jj+1)         ! keep this
        fkarf = FFT(kref2 * FFTb)            ! keep this

    ENDIF

    IF( jj .NE. nref-1) THEN

! use blas, updated by mohan
        CALL DCOPY(n123G, fkfa, 1, fkf0_saved, 1)
        CALL ZCOPY(k123G, kref2, 1, kref1_saved, 1)
        CALL DCOPY(n123G, fkarf, 1, fk1rf_saved, 1)
        CALL ZCOPY(k123G, kernel0a, 1, kernel0_saved, 1)
        CALL ZCOPY(k123G, kernel1a, 1, kernel1_saved, 1)

! old ones
!        fkf0_saved = fkfa
!        kref1_saved = kref2
!        fk1rf_saved = fkarf
!        kernel0_saved = kernel0a ! keep kernel for next iteration
!        kernel1_saved = kernel1a ! keep kernel for next teration
    ENDIF
    !---------------------------------------------------------------


    Call StopClock('intPot1')
    Call StartClock('intPot2')

    h = refXi(jj+1) - refXi(jj)
    t = (xi - refXi(jj)) / h
    t2 = t * t
    t3 = t * t2


!new ones in intPot()
!different from intEnergy because mask
!functions are used in two places
    mask = 1.d0
    IF(cutrho==0.d0)THEN
      WHERE( (xi .LT. refxi(jj)) .OR. (xi .GE. refxi(jj+1)))
        mask = 0.d0
      ENDWHERE
    ELSE
      WHERE( (xi .LT. refxi(jj)) .OR. (xi .GE. refxi(jj+1)) .OR. (rho .LT. cutrho))
        mask = 0.d0
      ENDWHERE
    ENDIF
    rhoMask = rhoA * mask

!    IF (cutrho == 0.d0) THEN
!      mask = stepfun(xi-refXi(jj)) * stepfun(refXi(jj+1)-xi)
!    ELSE
!      mask = stepfun(rho-cutrho) * stepfun(xi-refXi(jj)) * stepfun(refXi(jj+1)-xi)
!    ENDIF
!    rhoMask = rhoA * mask

    !============ PART I ==========

    h00 = 2.d0*t3-3.d0*t2+1.d0
    h10 = t3 - 2.d0*t2 + t
    h01 = 1.d0 - h00
    h11 = t3 - t2

    intpot1 = intpot1 +  mask * & 
       ( h00 * fkf0 + &        ! h00
         h01 * fkfa + &        ! h01
         h * ( h10 * fk1rf + h11 * fkarf) )

    Call StopClock('intPot2')
    Call StartClock('intPot3')

    !============ PART II ============

    ! FFTx is in G space
    ! here we calculate kernel(G)*rho^alpha(G)
    ! in order to use FFT, we set each zeta
    ! as a constant.
    FFTx = FFTx &
      + kernel0  * FFT( rhoMask * h00  )       &
      + kernel0a * FFT( rhoMask * h01  )       &
      + h * ( kref1 * FFT( rhoMask * h10  ) + kref2 * FFT( rhoMask * h11  ) )

    Call StopClock('intPot3')
    Call StartClock('intPot4')

    !=========== PART III ==========
    ! here h00 is the derivative of h00
    ! in order to save the memory
    ! we use array h00 instead of allocate a new one
    ! d (h(t))/ dt
    h00 = 6.d0*t2 - 6.d0*t
    h10 = 3.d0*t2 - 4.d0*t + 1.d0
    h01 = -h00
    !h11 = 3.d0*t2 - 2.d0*t
    h11 = h10 + 2.d0*t - 1.d0

!---- USE BLAS ------
!   CALL DSCAL(dimXYZ, 1.d0/h, h00, 1)
!( h00 * fkf0 + & ! used with DSCAL 

    potTmp = potTmp +  mask * & 
        ( (h00 * fkf0 + h01 * fkfa)/h + & 
          h10 * fk1rf + h11 * fkarf)

    Call StopClock('intPot4')

  ENDDO ! jj

  ! due to \partial (xi) / \partial rho
  pxiprho = kF/(3.d0*rho)*(1.d0-7.d0*hc_lambda*ss)
  potTmp = rhoA * potTmp
 
  ! due to  \partial (xi) / \partial grad(rho)
  ! note now xi is the derivative with respect to grad(rho)
  t3 = ( kF*hc_lambda*2.d0/rho83 ) * potTmp

  ! t is used as tmporary array here
  ! t2 and t3 is also used as temporay array here
  t = FFT(FFTx) 
  t2 = rho**(beta -1) * t

  !--------------------------------------
  ! Final potential due to this term
  ! intpot consists four parts
  !--------------------------------------
  intPot = &
  !---------
  ! part 1
  !---------
   alpha * rho**(alpha-1) * intpot1  &
  !---------
  ! part 2
  !---------
  + beta * t2 &
  !---------
  ! part 3
  !---------
  + pxiprho * potTmp &
  !---------
  ! part 4 
  !---------
  - ( & 
      FFT(imag*qVectors(:,:,:,1)*FFT(gradRho(:,:,:,1)*t3)) + &
      FFT(imag*qVectors(:,:,:,2)*FFT(gradRho(:,:,:,2)*t3)) + &
      FFT(imag*qVectors(:,:,:,3)*FFT(gradRho(:,:,:,3)*t3)))

  factor = cTF * c
  intPot = intPot * factor 

  !----------------------------
  ! PLEASE KEEP THE OLD ONE,
  ! OTHERWISE, SOMETIMES THE
  ! BLAS ROUTINES GIVE WRONG
  ! ANSWERS.
  !----------------------------
  ! old
  !energy = factor * SUM ( rho * t2 )
  !WRITE(*,*) " energy from intPot=", energy

  !-------------------
  ! new : using BLAS
  !-------------------
  energy = DDOT(n123G, rho, 1, t2, 1)
  energy = energy * factor

!  WRITE(*,*) " energy from intPot=", energy

  IF(maxval(abs(rho))<1e-20) THEN
    WRITE(*,*) 'no this spin channel electron'
    intpot=0.d0
  ENDIF

!! DEBUG
!!  CALL WrtOut(6,'(exit) intPot ')
!  print *,'max min intpot = ',maxval(intpot),minval(intpot)
!!  stop

  DEALLOCATE(refxi)

  CALL StopClock('intPot')

END FUNCTION intPot


SUBROUTINE MakeXiG(rho, xi, ss, gradRho, rho83)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Make xi(i,j,k) and ss(:,:,:) according to rho(i,j,k) 
! with formula xi = kF * (1+lambda*s^2)
! where s^2 = |grad_rho|^2/rho^(8/3)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!-------------------------------------------------------------------------------
! REVISION LOG:
! Created by Chen Huang (2008-Mar-16)
! updated by Mohan Chen (2013-08-06)
!--------------------------------------------------------------------------------
  IMPLICIT NONE
  
  ! >>> EXTERNAL VARS <<<<!
  REAL(KIND=DP), INTENT(IN)  :: rho(:,:,:)
  REAL(KIND=DP), INTENT(IN)  :: rho83(:,:,:)
  REAL(KIND=DP), INTENT(IN)  :: gradRho(:,:,:,:)
  REAL(KIND=DP), INTENT(OUT) :: xi(:,:,:)
  REAL(KIND=DP), INTENT(OUT) :: ss(:,:,:)

  ! >>> function body <<< !
  !=============================
  CALL StartClock('MakeXiG')
  
  ss = ( gradRho(:,:,:,1)**2 + & ! ss = |grad_rho|^2/rho^(8/3)
         gradRho(:,:,:,2)**2 + &
         gradRho(:,:,:,3)**2 ) / rho83
  xi = (3._DP * PI**2 * rho)**(1._DP/3._DP) * (1._DP + hc_lambda * ss)
  
  ! limit xi 
  WHERE(xi > 1000.d0)
    xi = 1000.d0
  ENDWHERE

  WHERE(xi < 0.001d0)
    xi = 0.001 
  ENDWHERE 

  CALL StopClock('MakeXiG')

  RETURN

END SUBROUTINE MakeXiG


SUBROUTINE MakeXi(rho, xi, ss)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Make xi(i,j,k) and ss(:,:,:) according to rho(i,j,k) 
! with formula xi = kF * (1+lambda*s^2)
! where s^2 = |grad_rho|^2/rho^(8/3)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!-------------------------------------------------------------------------------
! REVISION LOG:
! Created by Chen Huang (2008-Mar-16)
!--------------------------------------------------------------------------------
  IMPLICIT NONE
  
  ! >>> EXTERNAL VARS <<<<!
  REAL(KIND=DP), INTENT(IN)  :: rho(:,:,:)
  REAL(KIND=DP), INTENT(OUT) :: xi(:,:,:)
  REAL(KIND=DP), INTENT(OUT) :: ss(:,:,:)

  ! >>> local vars <<< 
  COMPLEX(KIND=DP), DIMENSION(SIZE(qMask,1),SIZE(qMask,2),SIZE(qMask,3)) :: &
    FFTrhoA                ! Fourier transform of rho for temp stroage.

  ! >>> function body <<< !
  !=============================
  CALL StartClock('MakeXi')
  
  FFTrhoA = FFT(rho)
  ss = ( FFT(imag*qVectors(:,:,:,1)*FFTrhoA)**2 + & ! ss = |grad_rho|^2/rho^(8/3)
         FFT(imag*qVectors(:,:,:,2)*FFTrhoA)**2 + &
         FFT(imag*qVectors(:,:,:,3)*FFTrhoA)**2 ) / rho**(8._DP/3._DP)
  xi = (3._DP * pi**2 * rho)**(1._DP/3._DP) * (1._DP + hc_lambda * ss)
  
  ! limit xi 
  WHERE(xi > 1000.d0)
    xi = 1000.d0
  ENDWHERE

  WHERE(xi < 0.001d0)
    xi = 0.001 
  ENDWHERE 

  CALL StopClock('MakeXi')

END SUBROUTINE MakeXi


SUBROUTINE MakeHC10Kernel(xi,k0,k1)
!------------------------------------------------------------------------------
! DESCRIPTION: 
!   Used for kernel interpolation. Compute the kernels: k0,k1 in q-space
!
!   k0 = FFT[w(r)]
!   k1 = FFT[dw(r)/dr * r]
!
! CONDITIONS AND ASSUMPTIONS
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Nov/10/2009: Chen Huang (Created)
!
!------------------------------------------------------------------------------
  USE IntKernelODE, ONLY: &
    ode_beta, &
    wInf, &    
    int_eta=>eta,&
    int_w=>w, &
    int_wp=>wp, &
    makeKernel_ode1, &
    intKernelODE_Clean => Clean

  IMPLICIT NONE

  !>>>>>> external vars <<<<<<<<!
  REAL(KIND=DP), INTENT(IN) :: xi

  COMPLEX(KIND=DP), INTENT(OUT) :: k0(:,:,:), k1(:,:,:)

  !>>>>>> internal vars <<<<<<<<!
 
  INTEGER :: qx, qy ,qz
  INTEGER :: rindex, num_eta  

  REAL(KIND=DP) :: etaStep 
  REAL(KIND=DP) :: eta
  REAL(KIND=DP) :: tmpw, tmpwp

  CHARACTER (len=500) :: message

  CALL StartClock('MakeHC10Kernel')
  
  !>>>>>>>>>>> function body <<<<<<<<<<<<<<!

!  OPEN (unit=1981, file='kernelw.out', status = 'unknown', form='formatted', action='write')
!  OPEN (unit=1982, file='kernelwp.out', status = 'unknown', form='formatted', action='write')

  ! Only solve for the kernel once
  IF( have_RK .eqv. .FALSE. ) THEN
    have_RK = .TRUE.
    WRITE(message,*) 'Making Huang-Carter 10 KEDF Kernel ...  '
    CALL WrtOut(6,message)
    ode_beta = beta
    wInf = -8.d0/3.d0/((5.d0-3.d0*ode_beta)*ode_beta)
    CALL makeKernel_ode1
    int_w(1) = 0.d0
    int_w    = int_w
    int_wp   = int_eta * int_wp   ! make int_wp => w'*eta
  END IF
    
  num_eta = size(int_eta)
  etaStep = int_eta(2)-int_eta(1) ! uniform eta 1-D mesh

  ! loop over q-space filling the kernels 
  DO qz = 1, SIZE(qTable,3)
    DO qy = 1, size(qTable,2)
      DO qx = 1, size(qTable,1)
    
        eta = qTable(qx, qy, qz) / (2._DP * xi)
        rindex = floor(eta/etaStep) + 1 

        if (rindex>=num_eta)  then 
          tmpw = int_w(num_eta)
          tmpwp = int_wp(num_eta)
        else
          tmpw  = (int_w(rindex)  + int_w(rindex+1))/2.d0
          tmpwp = (int_wp(rindex) + int_wp(rindex+1))/2.d0
        endif
          
        k0(qx,qy,qz) = tmpw/(8.d0*xi**3)
        k1(qx,qy,qz) = (-tmpwp-3.d0*tmpw)/(8.d0*xi**3)  ! FFT[ dw/d\Delta(r) * \Delta(r) ]

      ENDDO
    ENDDO
  ENDDO

  CALL StopClock('MakeHC10Kernel')

  RETURN

END SUBROUTINE MakeHC10Kernel


FUNCTION h00(t)
!---------------------------------------------------
! hermite function basis h00
! Used for cubic Hermite spline
!---------------------------------------------------
  IMPLICIT NONE

  REAL(KIND=dp), intent(in) :: t(:,:,:)
  REAL(KIND=dp) :: h00(size(t,1), size(t,2), size(t,3))

  h00 = 2._DP * t**3 - 3._DP*t**2 + 1._DP
  return

END FUNCTION h00

FUNCTION h10(t)
!---------------------------------------------------
! hermite function basis h10
! Used for cubic Hermite spline
!---------------------------------------------------
  IMPLICIT NONE

  REAL(KIND=dp), intent(in) :: t(:,:,:)
  REAL(KIND=dp) :: h10(size(t,1), size(t,2), size(t,3))

  h10 = t**3 - 2._DP*t**2 + t
  return

END FUNCTION h10

FUNCTION h01(t)
!---------------------------------------------------
! hermite function basis h10
! Used for cubic Hermite spline
!---------------------------------------------------
  IMPLICIT NONE

  REAL(KIND=dp), intent(in) :: t(:,:,:)
  REAL(KIND=dp) :: h01(size(t,1), size(t,2), size(t,3))

  h01 = -2._DP*t**3 + 3._DP*t**2 
  return 

END FUNCTION h01

FUNCTION h11(t)
!---------------------------------------------------
! hermite function basis h10
! Used for cubic Hermite spline
!---------------------------------------------------
  IMPLICIT NONE

  REAL(KIND=dp), intent(in) :: t(:,:,:)
  REAL(KIND=dp) :: h11(size(t,1), size(t,2), size(t,3))

  h11 = t**3 - t**2 
  return 

END FUNCTION h11


END MODULE KEDF_HC10 

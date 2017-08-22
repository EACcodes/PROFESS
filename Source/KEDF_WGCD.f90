MODULE KEDF_WGCD 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_WGCD 
!     |_SUBROUTINE InitializeWGCD
!     |_SUBROUTINE ComputeFrMatrix
!     |_SUBROUTINE DecomposeDensityKEDF
!
! DESCRIPTION:
!   This modules provides variables and functions that are needed in WGCD. 
!
! REFERENCES:
!   Junchao Xia, Emily A. Carter Phys. Rev. B 86, 235109 (2012) 
!   "Density-decomposed orbital-free density functional theory for 
!    covalently-bonded molecules and materials"
!------------------------------------------------------------------------------
! REVISION LOG:
!   04/15/2011 Steven write the WGCD code. 
!   01/31/2012 Mohan rewrite the WGCD code.
!------------------------------------------------------------------------------
                              !<< GLOBAL >>

  USE Constants, ONLY: DP, PI, IMAG
  USE OutputFiles
  USE MPI_Functions
  USE MathFunctions, ONLY: stepfun
  USE MathSplines, ONLY: spline_cubic_set, spline_cubic_val    
  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu
  USE KEDF_WTkernel, ONLY: keKernel
  USE KEDF_WTkernel, ONLY: fillWT    ! fill the WT kernel first! 
  USE KEDF_WTkernel, ONLY: alpha, beta
  USE KEDF_WGCkernel, ONLY: FillWGC
  USE KEDF_WGC, ONLY : rhov
  USE KEDF_WGC, ONLY: WGCPotentialPlus 
  USE SYS, ONLY: rho0, rhoS,  bvac 

  IMPLICIT NONE

  INTEGER :: WGCDiter = 0
  LOGICAL :: WGCDflag = .FALSE.
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: Fr
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: oldFr
  REAL(KIND=DP) :: oldalpha, oldbeta, yval, ypval, yppval, & ! interpolation use
                   sumrho, sumsize, DFr, maxFr, minFr
  INTEGER :: intn, ke1, ke2, ke3
  INTEGER :: filestatus

  REAL(KIND=DP), DIMENSION(151) :: t 
  ! (Steven add) Data of the scaling function to be interpolated in WGCD model.
  !
  REAL(KIND=DP), DIMENSION(151) :: scalefunt 
  REAL(KIND=DP), DIMENSION(151) :: scalefunDD 

  REAL(KIND=DP) :: rhoc = 6.84e-3 
  ! Can be changed from the input
  !
  REAL(KIND=DP) :: scfc = 1e-4
  ! F(r) self consistent convergence criteria. (Steven, mohan add 12-12-12)
  !
  REAL(KIND=DP) :: shiftm = 0.d0  
  ! Parameters for decomposition. (Steven, mohan add 12-12-12)
  !
  REAL(KIND=DP) :: mrhos = 1.d0   
  ! Parameters for decomposition. (Steven, mohan add 12-12-12) 
  !

CONTAINS


SUBROUTINE InitializeWGCD(dimX, dimY, dimZ, kinetic, pot_tol)
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
    
  INTEGER, INTENT(IN) :: dimX, dimY, dimZ
  INTEGER, INTENT(INOUT) :: kinetic 
  ! type of KEDF
  !
  REAL(KIND=DP), INTENT(INOUT) :: pot_tol 
  ! potential tolerance for electron optimization
  !

  CALL StartClock('InitWGCD')

  IF(ALLOCATED(Fr)) THEN
    !WRITE(outputUnit,*) "Has already allocate Fr and oldFr"
  ELSE
    ALLOCATE( Fr(dimX, dimY, dimZ) )
    ALLOCATE( oldFr(dimX, dimY, dimZ) )
  ENDIF

  WGCDiter = 0
  WGCDflag = .FALSE.
  Fr = 1.d0  
    
  ! for first iterations (guesses), use a very large pot_tol
  pot_tol = pot_tol*1e-10 + 1e-2

  !in ion relaxation, each time it enters rhoOptimizer, set KEDF to be WT
  WGCDflag = .TRUE.

  IF(rankGlobal==0) THEN
    WRITE(*,*) ' WGCD: first iteration, set KEDF as WT'
  ENDIF
  kinetic = 4

  oldalpha = alpha
  oldbeta = beta
  alpha = 5.0_DP / 6.0_DP
  beta = 5.0_DP / 6.0_DP

  ke1 = SIZE(kekernel,1)
  ke2 = SIZE(kekernel,2)
  ke3 = SIZE(kekernel,3)

  DEALLOCATE(kekernel)
  ALLOCATE(kekernel(ke1,ke2,ke3,1))

  CALL FillWT()

  CALL StopClock('InitWGCD')

  RETURN
      
END SUBROUTINE InitializeWGCD


SUBROUTINE ComputeFrMatrix(rho, dimX, dimY, dimZ, nspin, new_iteration, kinetic, pot_tol)
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

  ! >> EXTERNAL VARIABLES << !

  INTEGER, INTENT(IN) :: dimX, dimY, dimZ, nspin
  ! dimension of rho
  !
  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin), INTENT(IN) :: rho
  ! charge density on FFT box.
  !
  LOGICAL, INTENT(INOUT) :: new_iteration
  ! if new_iteration == .TRUE. goto 10 in RhoOptN.f90.
  !
  INTEGER, INTENT(INOUT) :: kinetic 
  ! type of KEDF
  !
  REAL(KIND=DP), INTENT(INOUT) :: pot_tol 
  ! potential tolerance for electron optimization
  !

  ! >> INNER VARIABLES << !

  INTEGER :: ii, jj, kk

  ! >> INITIALIZE << !

  CALL StartClock('Fr')

  new_iteration = .FALSE.
  
  ! >> FUNCTION << !


  IF(WGCDflag .EQV. .TRUE.) THEN
  ! interpolate the scaling function to compute Fr matrix

    WRITE(*,*) 'WGCDiter=',WGCDiter
    intn = 151
    IF(WGCDiter==0) THEN
      ! now set everything ready for WGC
      IF(rankGlobal==0) THEN
        WRITE(*,*) 'WGCD: First iteration done, set all thing to WGC'
      ENDIF
      kinetic = 12
      beta = oldbeta
      alpha = oldalpha
      ke1 = SIZE(keKernel,1)
      ke2 = SIZE(keKernel,2)
      ke3 = SIZE(keKernel,3)
      DEALLOCATE(keKernel)
      ALLOCATE(keKernel(ke1,ke2,ke3,4))
    ENDIF
 
    WGCDiter = WGCDiter + 1
    IF(rankGlobal==0) THEN
      WRITE(*,*) 'WGCD: ITERATION ',WGCDiter, ': calculating Fr matrix'
    ENDIF

    oldFr = Fr

    IF(bvac .EQV. .TRUE.) THEN
      rho0 = SUM(rho(:,:,:,1)*Fr*stepfun(rho(:,:,:,1)*Fr-rhoc)) &
            /SUM(rho(:,:,:,1)/rho(:,:,:,1)*stepfun(rho(:,:,:,1)*Fr-rhoc))
    ELSE
      rho0 = SUM(rho(:,:,:,1)*Fr)/SIZE(rho(:,:,:,1))
    ENDIF

#ifdef __USE_PARALLEL
    IF(bvac .EQV. .TRUE.) THEN
        sumrho = SUM(rho(:,:,:,1)*Fr*stepfun(rho(:,:,:,1)*Fr-rhoc))
        sumsize = SUM(rho(:,:,:,1)/rho(:,:,:,1)*stepfun(rho(:,:,:,1)*Fr-rhoc))
    ELSE
        sumrho =  SUM(rho(:,:,:,1)*Fr)
        sumsize = SIZE(rho(:,:,:,1))
    ENDIF
    CALL ReduceRealLevel1(sumrho)
    CALL ReduceRealLevel1(sumsize)
    rho0=sumrho/sumsize
#endif

    DO kk=lbound(rho,3),ubound(rho,3)
       DO jj=lbound(rho,2),ubound(rho,2)
          DO ii=lbound(rho,1),ubound(rho,1)
             IF(rho(ii,jj,kk,1)/rho0-shiftm<0) THEN
               Fr(ii,jj,kk) = 1.d0
             ELSE IF(rho(ii,jj,kk,1)/rho0-shiftm>10.d0) THEN
               Fr(ii,jj,kk) = (1.2 + 0.00556*((1.d0 - 0.9**((rho(ii,jj,kk,1)/rho0-shiftm-3.4)/0.1))/0.1 - 1.d0)) &
                            /(rho(ii,jj,kk,1)/rho0-shiftm)
             ELSE
               CALL spline_cubic_val(intn, t(1:intn), scalefunt(1:intn),scalefunDD(1:intn), &
                                     rho(ii,jj,kk,1)/rho0-shiftm,  yval, ypval, yppval )
               Fr(ii,jj,kk) = yval
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF(bvac .EQV. .TRUE.) THEN
        rho0 = Sum(rho(:,:,:,1)*Fr*stepfun(rho(:,:,:,1)*Fr-rhoc)) &
              /sum(rho(:,:,:,1)/rho(:,:,:,1)*stepfun(rho(:,:,:,1)*Fr-rhoc))
        rhoS = maxval(rho(:,:,:,1)*Fr)/1.5d0
    ELSE
        rho0 = Sum(rho(:,:,:,1)*Fr)/size(rho(:,:,:,1))
        rhoS = rho0
    ENDIF

#ifdef __USE_PARALLEL
    IF(bvac .EQV. .TRUE.) THEN
        sumrho = SUM(rho(:,:,:,1)*Fr*stepfun(rho(:,:,:,1)*Fr-rhoc))
        sumsize = SUM(rho(:,:,:,1)/rho(:,:,:,1)*stepfun(rho(:,:,:,1)*Fr-rhoc))
    ELSE
        sumrho =  Sum(rho(:,:,:,1)*Fr)
        sumsize = SIZE(rho(:,:,:,1))
    ENDIF

    CALL ReduceRealLevel1(sumrho)
    CALL ReduceRealLevel1(sumsize)
    rho0=sumrho/sumsize

    IF(bvac .EQV. .TRUE.) THEN
      rhoS=rho0/1.5d0
    ELSE
      rhoS=rho0
    ENDIF
#endif

    CALL FillWGC()

    DFr = MAXVAL(ABS(Fr-oldFr))
#ifdef __USE_PARALLEL
    CALL MPI_ALLREDUCE(DFr,DFr,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
#endif

    maxFr = MAXVAL(Fr)
    minFr = MINVAL(Fr)
#ifdef __USE_PARALLEL
    CALL MPI_ALLREDUCE(maxFr,maxFr,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
    CALL MPI_ALLREDUCE(minFr,minFr,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif
    
    ! if Dfr is less than some critical value, set to normal pot tolerance
    IF( DFr<1.0e-4 .AND. pot_tol>0.99e-2 ) THEN
      pot_tol = (pot_tol-1e-2) * 1e10
    ENDIF

    IF( DFr<scfc .AND. WGCDiter>1) THEN
      IF(rankGlobal==0) THEN
         WRITE(*,*) 'WGCD: Scaling function F(r) now is fully self-consistent. Job finished!'
         WRITE(*,*) 'WGCD: Output density info'
         WRITE(*,*) 'WGCD: rho0 of rhodel is: ', rho0
         WRITE(*,*) 'WGCD: rhoS in WGC is: ', rhoS
         WRITE(*,*) 'WGCD: maximum in F(r): ', maxFr
         WRITE(*,*) 'WGCD: minimum in F(r): ', minFr
      ENDIF
    ELSE 
      IF(rankGlobal==0) THEN
        WRITE(*,*) 'WGCD: Decomposing: Density is rescaled. Now round back'
        WRITE(*,*) 'WGCD: Maxium F(r) difference is: ',DFr
        WRITE(*,*) ' '
      ENDIF

      IF(WGCDiter>=30 .and. DFr>0.01 .or. WGCDiter>70) THEN
        IF(rankGlobal==0) THEN
           WRITE(*,*) 'WGCD: Self-consistent F(r) calculation fails. Program stops!'
        ENDIF
      ELSE
        new_iteration = .true.
        RETURN
      ENDIF
    ENDIF
  ENDIF ! kinetic=12

  CALL StopClock('Fr')

  RETURN
    
END SUBROUTINE ComputeFrMatrix


SUBROUTINE DecomposeDensityKEDF(rho,Fr,energy,potential,vaccutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Density decomposition of the total density to rho1 and rho2
!   rho1 = rho * F(rho/rho0) and rho2 = rho - rho1
!   Thus, rho1 is a smooth density, and rho2 is localized density in bonding area
!
!   Then T[rho] = (T1[rho]-T1[rho1]) + T2[rho1]
!   T1 is chosen as lambda*TF + mu*vW
!   T2 is chosen as TF + vW + WGC
!   But here WGC kernel is remade with new Rho0 and RhoS calculated with rho1
!
!   First calculate WGC then get reasonable density
!   so that we know where is bonding region. We then calculate F(r) matrix and
!   keep it constant in the second round optimization. WGC kernel and rho0, rhoS are
!   also updated just once after the first round optimization is done.
!
!   Afterwards, F(r) matrix is re-evaluated. And loop goes on until F(r) self-consistent
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!------------------------------------------------------------------------------
! REVISION LOG:
!   3/26/2003  created  (steven)
!------------------------------------------------------------------------------
  USE PlaneWave, ONLY: qTable
  USE FOURIER, ONLY : FFT

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: Fr
  ! scale function
  !
  REAL(KIND=DP), INTENT(OUT) :: energy
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential
  !
  LOGICAL, INTENT(IN) :: vaccutoff
  !

                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), PARAMETER :: cTF = 2.87123400018819_DP

  REAL(KIND=DP) :: ekin, e1 

  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: rho1
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: rho23
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: rho123
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: sqrtrho
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: sqrtrho1
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: kinPot
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: temp2 
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: potTF 

  REAL(KIND=DP) :: sumTF = 0.d0
  REAL(KIND=DP) :: sumVW = 0.d0

                       !>> INITIALIZATION <<!

  CALL StartClock('WGCD')

  e1 = 0.d0
  temp2 = 0.d0
  ekin = 0.d0
  kinPot = 0.d0
  sqrtrho=sqrt(rho)
  rho23=rho**(2.d0/3.d0)
  ! first decompose density
  rho1 = rho * Fr
  rho123=rho1**(2.d0/3.d0)
  sqrtrho1=sqrt(rho1)
  

                        !! >> FUNCTION << !!

  !-------------------------------------------------------------
  ! first part of the interaction KEDF involving global density
  ! lambda is the coefficient for TF KEDF
  ! while mu is the coefficient for VW KEDF
  ! calculate T1[rho] term
  ! TF
  ekin = ekin + lambda * cTF*SUM(rho*rho23)
  kinPot = kinPot + lambda * 5.d0/3.d0 * cTF * rho23
  ! vW
  temp2 = FFT( FFT(sqrtrho) * qTable**2)
  ekin = ekin + mu * 0.5_DP * SUM(sqrtrho*temp2) 
  kinPot = kinPot + mu * temp2/(2._DP*sqrtrho)

  !-------------------------------------------------------------
  ! seoncd part of the interaction KEDF involving localized density
  ! Now calculate -T1[rho1]
  ! TF
  sumTF = cTF*SUM(rho1*rho123)
  ekin = ekin - lambda * sumTF 
  potTF = Fr * 5.d0/3.d0 * cTF * rho123 
  kinPot = kinPot - lambda * potTF 
  ! vW
  temp2 = FFT( FFT(sqrtrho1) * qTable**2)
  sumVW = SUM(sqrtrho1*temp2)
  temp2 = temp2 * Fr / (2._DP*sqrtrho1)
  ekin = ekin - mu * 0.5_DP * sumVW 
  kinPot = kinPot - mu * temp2 

  ! the last term: T2[rho1]
  ! TF
  ekin = ekin + sumTF 
  kinPot = kinPot + potTF 
  ! vW
  ekin = ekin + 0.5_DP * sumVW 
  kinPot = kinPot + temp2

  ! finally WGC, rho0, rhoS, and WGC kernel are updated in the RhoOptimization.F90 subroutine
  ! simply just use rho1 instead of rho
  CALL WGCPotentialPlus(rho1, temp2, .TRUE., e1, vacCutoff)
  ekin = ekin + e1
  kinPot = kinPot + temp2 * Fr

  energy = ekin
  potential = kinPot

  CALL StopClock('WGCD')

  RETURN

END SUBROUTINE DecomposeDensityKEDF

END MODULE KEDF_WGCD

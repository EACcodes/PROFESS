MODULE SetupKEDF
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE SetupKEDF
!     |_SUBROUTINE SetupFunctional 
!
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   1. Watson, S.C. and Carter, E.A.  "Linear-Scaling Parallel Algorithms for
!      the First-Principles Treatment of Metals."  Computer Physics 
!      Communications, 128 (2000) 67-92
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 12/15/03  File created as a spliter off of the disorganized DataTypes.
!             Also added commentary.  (GSH)
! 01/08/03  Took out spin from qTable and qMask. (GSH)
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP
  USE MPI_Functions
  USE Output, ONLY: WrtOut
  USE Output, ONLY: QUIT
  USE OutputFiles, ONLY: outputUnit

  IMPLICIT NONE


CONTAINS

SUBROUTINE SetupFunctional
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine takes the number of points on a real space grid, and 
!   initializes qTable, qMask, etc. and sets the FUNCTIONAL module to start
!   working.  It must be run once per program execution.
!
!   In the event that the WT or WGC kinetic energy
!   functionals are used, it will also allocate and fill the kernel table. Once
!   again, the allocation and filling parts are kept separate to allow updating
!   of the kernel when the cell size is altered.
!
! REFERENCES:
!   1.  Ashcroft and Mermin, Solid State Physics.  Harcourt College Publishers,
!       Fort Worth, 1976.
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/30/2003  File created.  (Vincent Ligneres)
!   11/19/2003  Checking validity on Al_fcc (4.03A) Grid size does not match
!               that of Stuart's code. I get 16 he has 24 for 1200eV cut. (VLL)
!   11/20/2003  Corrected formulas numX, Y, Z. Tested on Al_fcc, it works. &
!               (VLL)
!   11/22/2003  Cosmetic Changes (Greg Ho)
!   12/05/2003  Added the kernel table and the FillKernel call. (VLL)
!   12/12/2003  Removed call to FillQTable, actually removed the whole 
!               procedure.  Replaced with a few lines of code. (GSH)
!   01/08/2004  Renamed this procedure "SetupGridObjects" (GSH)
!   03/23/2005  Added initialization for the Hybrid KEDF. (VLL)
!   09/21/2006  Fixed indices to allow parallelization (LH)
!------------------------------------------------------------------------------

  USE SYS, ONLY: rho0, rhoS, hold0, holdS

  USE CellInfo, ONLY: kinetic
  USE CellInfo, ONLY: cell, numSpin
  USE CellInfo, ONLY: k1G, k2G, k3G, k3Goff
  USE CellInfo, ONLY: m1G, m3G, n3Goff
  USE CellInfo, ONLY: n1G, n2G, n3G

  USE IonElectronSpline, ONLY: iiSpline 
  USE IonElectronSpline, ONLY: ieSpline 
  USE CBSpline, ONLY: FillB
  USE CBSpline, ONLY: splineOrder
  USE CBSpline, ONLY: bSpline1
  USE CBSpline, ONLY: bSpline2
  USE CBSpline, ONLY: bSpline3

  USE PlaneWave, ONLY: qMask
  USE PlaneWave, ONLY: qVectors
  USE PlaneWave, ONLY: qTable

  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu
  USE KEDF_WTkernel, ONLY: alpha, beta
  USE KEDF_WTkernel, ONLY: keKernel 
  USE KEDF_WGCkernel, ONLY: firstOrderWGC
  USE KEDF_WGCkernel, ONLY: gamma
  USE KEDF_WGCkernel, ONLY: alpha5
  USE KEDF_WGCkernel, ONLY: beta5
  USE KEDF_HC10, ONLY: hc_lambda
  USE KEDF_HC10, ONLY: hc_lambda_val
  USE KEDF_WGCD, ONLY: t
  USE KEDF_WGCD, ONLY: scalefunt
  USE KEDF_WGCD, ONLY: scalefunDD
  USE KEDF_EvW, ONLY: kloc, aloc
  USE KEDF_GGA, ONLY: numint, d, ELFvsd, ELFvsdDD

  USE MathSplines, ONLY: spline_cubic_set ! For Steven's KEDF.
  USE MathSplines, ONLY: spline_cubic_val

  IMPLICIT NONE

                       !>> INTERNAL VARIABLES <<! 
  INTEGER :: fileStatus          
  ! Checks that tables are allocated alright.
  !
  INTEGER :: ii
  ! For Steven's KEDF.
  !
  REAL(kind=DP) :: tmp
  ! another temp var
  !
  CHARACTER(len=500) :: message
  !

                           !>> INITIALIZATION <<!    

  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           SETUP KEDF"
  WRITE(outputUnit,*) "kinetic energy density functional number is ", kinetic
  WRITE(outputUnit,*) "( 1=TF, 2=VW, 3=TF+VW, 4=WT,  5=WGC,  "
  WRITE(outputUnit,*) "  7=LQ, 8=HQ, 10=CAT,  11=HC, 12=WGCD,"
  WRITE(outputUnit,*) "  15=GGA, 16=WGCD+GGA, 17=EVT, 18=EVC )"

                           !>> FUNCTION BODY <<!

  ! Allocate a bunch of memory
  ! Periodic boundary condition, 
  ! 'q' is plane wave
  ! qTable is for storing |q|,
  ! qVectors is for each q (q_x, q_y, q_z)

  WRITE(outputUnit,'(A,3I5)') " qTable dimension: ", k1G, k2G, k3G

  ! because we only use Gamma k point, the number of plane waves
  ! can be deduced to almost half of the original set of plane waves.
  ALLOCATE(qTable(k1G, k2G, k3G), stat=fileStatus)
  IF (fileStatus/=0) THEN
    WRITE(*,*)'Error allocating the q-norm table. Leaving.'
    STOP
  END IF

  ALLOCATE(qVectors(k1G, k2G, k3G, 3), stat=fileStatus)
  IF (fileStatus/=0) THEN
    WRITE(*,*)'Error allocating the q-vector table. Leaving.'
    STOP
  END IF

  ALLOCATE(qMask(k1G, k2G, k3G), stat=fileStatus)
  IF (fileStatus/=0) THEN
    WRITE(*,*)'Error allocating the q-in/out table. Leaving.'
    STOP
  END IF

  ! If at least one of particle-mesh Ewald or Choly-Kaxiras ion-electron
  ! calculations are being used
  IF (ieSpline .OR. iiSpline) THEN

    ALLOCATE(bSpline1(k1G), stat=fileStatus)
    IF (fileStatus/=0) STOP &
      'Error allocating the b-spline table 1. Leaving.'

    ALLOCATE(bSpline2(k2G), stat=fileStatus)
    IF (fileStatus/=0) STOP &
      'Error allocating the b-spline table 2. Leaving.'

    ALLOCATE(bSpline3(k3G), stat=fileStatus)
    IF (fileStatus/=0) STOP &
      'Error allocating the b-spline table 3. Leaving.'

    ! Check that spline order is positive and even - change if necessary
    IF (splineOrder<=0) THEN
      splineOrder = 10
      WRITE(*,*) "Spline order must be positive and even: changing to 10"
    ELSE IF (MOD(splineOrder,2)==1) THEN
      splineOrder = splineOrder + 1
      WRITE(*,*) "Spline order must be even - changed to", splineOrder
    END IF

    CALL FillB(bSpline1,bSpline2,bSpline3,splineOrder,m1G,m3G,n3Goff)

  END IF


  ! If we use a non-trivial KEDF we have to initialize the parameters
  SELECT CASE (kinetic)

!------------------------------------------------------------------------------
  CASE(1) ! TF
    ! This is technically lambda * TF. We do nothing.
    IF (lambda<-99.0) lambda = 1.0_DP
    mu = 0._DP

    WRITE(outputUnit,'(A,F10.4)') " Coefficient for TF KEDF is ", lambda
!------------------------------------------------------------------------------
  CASE(2) ! VW
    ! Technically, it's actually mu * VW.
    IF (mu<-99._DP) mu = 1.0_DP
    lambda = 0._DP

    WRITE(outputUnit,'(A,F10.4)') " Coefficient for VW KEDF is ", mu
!------------------------------------------------------------------------------
  CASE(3) ! lambda * TF + mu * VW
    ! Do nothing. Default is TF + 1/9 vW.
    IF (lambda<-99) lambda = 1.0_DP
    IF (mu<-99) mu = 1.0_DP/9._DP

    WRITE(outputUnit,'(A,F10.4)') " Coefficient for TF KEDF is ", lambda
    WRITE(outputUnit,'(A,F10.4)') " Coefficient for VW KEDF is ", mu
!------------------------------------------------------------------------------
  CASE(4, 17) ! 4 is Wang-Teter, standard and custom response
              ! 17 is EvW based on Wang-Teter
    ! This is now lambda TF + mu vW + WT(alpha,beta).
    IF (lambda<-99) lambda = 1.0_DP 
    IF (mu<-99) mu = 1.0_DP

    WRITE(outputUnit,'(A,F10.4)') " Coefficient for TF KEDF is ", lambda
    WRITE(outputUnit,'(A,F10.4)') " Coefficient for VW KEDF is ", mu

    ! At this point, if alpha, beta and rho0 have not been set - we have to do
    ! it assuming standard WT calculation.
    IF (alpha<0._DP) alpha = 5._DP/6._DP
    IF (beta<0._DP) beta = alpha
    IF (rho0>0._DP) THEN
      hold0 = .TRUE. ! This should be redundant but one cannot be too careful.
    ELSE
      rho0 = cell%numEle/cell%vol
    END IF

    WRITE(outputUnit,'(A,F10.4)') " Coefficient for alpha in WT KEDF is ", alpha
    WRITE(outputUnit,'(A,F10.4)') " Coefficient for beta  in WT KEDF is ", beta 
    WRITE(outputUnit,'(A,F10.4)') " rho0 (Electron/Bohr^3) in WGC KEDF is ", rho0 
    WRITE(outputUnit,*) "hold0 = ", hold0 

    ! Periodic Boundary Condition
    ! the pre-FFT array, before padding, is saved for convolutions.
    ALLOCATE(keKernel(k1G, k2G, k3G,1), stat=fileStatus)

    IF (fileStatus/=0) THEN
      WRITE(message,*)'Error allocating the WT kernel table. Leaving.'
      CALL Error(6, message)
    END IF 

    IF (kloc<0._DP) kloc = 0._DP       ! for SC KEDF (Shin, mohan add 09-26-12)
    IF (aloc<0._DP) aloc = 0._DP       ! for SC KEDF (Shin, mohan add 09-26-12)
!------------------------------------------------------------------------------
  CASE(5, 18) ! Wang-Govind-Carter mohan add 2013-07-31
  ! This is lamba * TF + mu * vW + WGC(alpha, beta, gamma).
    IF (lambda<-99) lambda = 1.0_DP
    IF (mu<-99) mu = 1.0_DP 
    IF (firstOrderWGC==-100) THEN
      CALL WrtOut(6,'You have not specify the keyword of WGCT in your .inpt file')
      CALL WrtOut(6,'WGCT is a mandatory keyword, see manual for details')
      CALL WrtOut(6,'code stoped.')
      STOP
    ENDIF

    WRITE(outputUnit,'(A,F10.4)') " Coefficient for TF KEDF is ", lambda
    WRITE(outputUnit,'(A,F10.4)') " Coefficient for VW KEDF is ", mu

    ! keKernel is a very important array for storage of KEDF kernels.
    ! For small box, the kernel should be allocated in each small box
    ! with buffer.
    ALLOCATE(keKernel(k1G, k2G, k3G, 4), stat=fileStatus)
    IF (fileStatus/=0) THEN
      WRITE(message,*)'Error allocating the WGC kernel table. Leaving.'
      CALL Error(6, message)
    END IF

    IF (fileStatus/=0) THEN
      WRITE(message,*)'Error allocating the WGC kernel table. Leaving.'
      CALL Error(6, message)
    END IF 
!

!!
    ! Now we set alpha, beta, gamma and rho0 if they have not been set yet.
    ! Using the "ideal" Alex values.
    IF (alpha<0._DP) alpha = (5.0_DP + SQRT(5.0_DP)) / 6.0_DP
    IF (beta<0._DP) beta = (5.0_DP - SQRT(5.0_DP)) / 6.0_DP
    IF (gamma<0._DP) gamma = 2.7_DP
    IF (alpha5<0._DP) alpha5 = 1.0_DP  ! for wgct -5
    IF (beta5<0._DP)  beta5  = 1.0_DP  ! for wgct -5
    IF (alpha < beta) THEN
      tmp = alpha
      alpha = beta
      beta = tmp
      CALL WrtOut(6,'alpha < beta for WGC KEDF, I interchanged these two.')
    ENDIF

    WRITE(outputUnit,'(A,F10.4)') " Coefficient for alpha in WGC KEDF is ", alpha 
    WRITE(outputUnit,'(A,F10.4)') " Coefficient for beta  in WGC KEDF is ", beta 
    WRITE(outputUnit,'(A,F10.4)') " Coefficient for gamma in WGC KEDF is ", gamma
      
    IF (rho0>0._DP) THEN
      hold0 = .TRUE.
    ELSE
      rho0 = cell%numEle/cell%vol 
    END IF
    IF (rhoS>0._DP) THEN 
      holdS = .TRUE.
    ELSE
      IF (holdS .EQV. .TRUE.) THEN
        rhoS = rho0
      ELSE
        rhoS = cell%numEle/cell%vol 
      END IF
    END IF

    WRITE(outputUnit,'(A,F10.4)') " rho0 (Electron/Bohr^3) in WGC KEDF is ", rho0 
    WRITE(outputUnit,'(A,F10.4)') " rhoS (Electron/Bohr^3) in WGC KEDF is ", rhoS 
    WRITE(outputUnit,*) "hold0 = ", hold0 
    WRITE(outputUnit,*) "holdS = ", holdS 

    IF (kloc<0._DP) kloc = 0._DP       ! for SC KEDF (Shin, mohan add 09-26-12)
    IF (aloc<0._DP) aloc = 0._DP       ! for SC KEDF (Shin, mohan add 09-26-12)
!------------------------------------------------------------------------------
  CASE(7) ! LQ Functional (Jeng-Da Chai) 
    lambda = 1._DP ! Making sure we get 1*TF+1*VW+LQ. 
    mu = 1._DP 
    ! At this point, if rho0 has not been set we have to do 
    ! it assuming standard LQ calculation. 
    IF (rho0>0._DP) THEN 
      hold0 = .TRUE. ! This should be redundant but one cannot be too careful.
    ELSE 
      rho0 = cell%numEle/cell%vol 
    END IF 
!------------------------------------------------------------------------------
  CASE(8) ! HQ Functional (Jeng-Da Chai)
    lambda = 1._DP ! Making sure we get 1*TF+1*VW+HQ.
    mu = 1._DP
    ! At this point, if rho0 has not been set we have to do
    ! it assuming standard LQ calculation.
    IF (rho0>0._DP) THEN
      hold0 = .TRUE. ! This should be redundant but one cannot be too careful.
    ELSE 
      rho0 = cell%numEle/cell%vol 
    END IF 
!------------------------------------------------------------------------------
  CASE(10) ! CAT KEDF
    ALLOCATE(keKernel(k1G, k2G, k3G, 4), stat=fileStatus)
    IF (fileStatus/=0) THEN
      WRITE(*,*)'Error allocating the WGC kernel table. Leaving.'
      STOP
    END IF

    IF (rho0>0._DP) THEN
      hold0 = .TRUE.
    ELSE
      rho0 = cell%numEle/cell%vol 
    END IF
    IF (rhoS>0._DP) THEN 
      holdS = .TRUE.
    ELSE
      IF (holdS) THEN
        rhoS = rho0
      ELSE
        rhoS = cell%numEle/cell%vol 
      END IF
    END IF

!------------------------------------------------------------------------------
  CASE(11) ! Huang-Carter KEDF 
    ! This is TF + VW + Huang-Carter KEDF (2011) PRB .
    ALLOCATE(hc_lambda(n1G, n2G, n3G), stat=fileStatus)
    IF (fileStatus/=0) THEN 
      PRINT *, 'Error in allocating hc_lambda, code =', fileStatus, ' STOP' ; STOP
    ENDIF

    hc_lambda = hc_lambda_val
    lambda = 1.0_DP  ! coeff for TF
    mu     = 1.0_DP  ! coeff for VW

    
    ! the default is setting rhoS=-1.0, meaning the average density
    IF(rhoS==-1.d0) then
      rhoS = cell%numEle/cell%vol 
      WRITE(outputUnit,*) "rhoS=",rhoS
    ENDIF

    IF(rho0==-1.d0) then
      rho0 = cell%numEle/cell%vol 
      WRITE(outputUnit,*) "rho0=",rho0
    ENDIF
   
    IF(rhoS<0.d0 .AND. rho0<0.d0) THEN 
      PRINT *,'For HC10, either rho0 or rhoS needs to be given! error, stop' 
      WRITE(*,*) " rhoS=", rhoS
      WRITE(*,*) " rho0=", rho0
      STOP
    ENDIF
   
    IF (rhoS<0.d0) rhoS = rho0
    IF (rho0<0.d0) rho0 = rhoS

    IF (abs(rho0-rhoS) > 1e-10) then
      PRINT *,'rho0 must equal to rhoS! error, stop'
      STOP
    ENDIF

    IF (beta < 0.d0) then 
      PRINT *,' beta < 0.0 error stop'
    ELSE
      alpha = 8.d0/3.d0 - beta
    ENDIF
    
    WRITE(outputUnit,'(A,F12.6)') " coefficients for TF KEDF is ", lambda
    WRITE(outputUnit,'(A,F12.6)') " coefficients for VW KEDF is ", mu
    WRITE(outputUnit,'(A,F12.6)') " Lambda in HC10 is set to ", hc_lambda_val
    WRITE(outputUnit,'(A,F12.6)') " ALPHA  in HC10 is set to ", alpha
    WRITE(outputUnit,'(A,F12.6)') " BETA   in HC10 is set to ", beta
    WRITE(outputUnit,'(A,F12.6)') " Sum of alpha and beta is ", alpha+beta

!------------------------------------------------------------------------------
  CASE(12) ! Wang-Govind-Carter-based decomposition
    IF (lambda<-99) lambda = 1.0_DP
    IF (mu<-99) mu = 1.0_DP

    IF (firstOrderWGC==-100) THEN
      CALL WrtOut(6,'You have not specify the keyword of WGCT in your .inpt file')
      CALL WrtOut(6,'WGCT is a mandatory keyword, see manual for details')
      CALL WrtOut(6,'code stoped.')
      STOP
    ENDIF

    ! See comment in CASE(4) WT for why kernels saved differently depending
    ! on boundary conditions.
    ! Periodic boundary condition
    ALLOCATE(keKernel(k1G, k2G, k3G, 4), stat=fileStatus)

    IF (fileStatus/=0) THEN
      WRITE(*,*)'Error allocating the WGC kernel table. Leaving.'
      STOP
    END IF

    ! Now we set alpha, beta, gamma and rho0 if they have not been set yet.
    ! Using the "ideal" Alex values.
    IF (alpha<0._DP) alpha = (5.0_DP + SQRT(5.0_DP)) / 6.0_DP
    IF (beta<0._DP) beta = (5.0_DP - SQRT(5.0_DP)) / 6.0_DP
    IF (gamma<0._DP) gamma = 2.7_DP
    IF (alpha5<0._DP) alpha5 = 1.0_DP  ! for wgct -5
    IF (beta5<0._DP)  beta5  = 1.0_DP  ! for wgct -5
    IF (alpha < beta) THEN
      tmp = alpha
      alpha = beta
      beta = tmp
      CALL WrtOut(6,'alpha < beta for WGC KEDF, I interchanged these two.')
    ENDIF

    IF (rho0>0._DP) THEN
      hold0 = .TRUE.
    ELSE
      rho0 = cell%numEle/cell%vol 
    END IF
    IF (rhoS>0._DP) THEN
      holdS = .TRUE.
    ELSE
      IF (holdS .EQV. .TRUE.) THEN
        rhoS = rho0
      ELSE
        rhoS = cell%numEle/cell%vol 
      END IF
    END IF

    ! set up F(rho/rho0)
    DO ii=1,151
        t(ii) = -5.d0 + 0.1*(ii-1)
    ENDDO
    scalefunt = &
    ! this is the first tial, in the excel "scalefun" sheet 2, in red.
                (/1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, 1.000000000000, &
                  0.999411764680, 0.999411764536, 0.999411763575, 0.999411757168, 0.999411714453, &
                  0.999411429688, 0.999409531250, 0.999396875000, 0.999312500000, 0.998750000000, &
                  0.995000000000, 0.970000000000, 0.893750000000, 0.829230769231, 0.773928571429, &
                  0.726000000000, 0.684062500000, 0.647058823529, 0.614197530864, 0.584795321637, &
                  0.558333333333, 0.534391534392, 0.512626262626, 0.492753623188, 0.474537037037, &
                  0.457777777778, 0.442307692308, 0.427983539095, 0.414682539683, 0.402298850575, &
                  0.390740740741, 0.379928315412, 0.369791666667, 0.360269360269, 0.351307189542, &
                  0.342857142857, 0.334723333333, 0.326893945946, 0.319358115789, 0.312104809231, &
                  0.305122970100, 0.298401632283, 0.291930005506, 0.285697539724, 0.279693972439, &
                  0.273909362413, 0.268334112559, 0.262958984254, 0.257775104874, 0.252773970011, &
                  0.247947441550, 0.243287742544, 0.238787449669, 0.234439483858, 0.230237099631, &
                  0.226173873492, 0.222243691747, 0.218440737966, 0.214759480304, 0.211194658845, &
                  0.207741273078, 0.204394569610, 0.201150030187, 0.198003360071, 0.194950476812, &
                  0.191987499452, 0.189110738151, 0.186316684271, 0.183602000905, 0.180963513846, &
                  0.178398202998, 0.175903194209, 0.173475751528, 0.171113269849, 0.168813267961, &
                  0.166573381949, 0.164391358968, 0.162265051343, 0.160192411001, 0.158171484206, &
                  0.156200406588, 0.154277398448, 0.152400760328, 0.150568868822, 0.148780172631, &
                  0.147033188834, 0.145326499370, 0.143658747715, 0.142028635751, 0.140434920803, &
                  0.138876412848, 0.137351971876, 0.135860505398, 0.134400966096, 0.132972349599, &
                  0.131573692379, 0.130204069775, 0.128862594109, 0.127548412916, 0.126260707264, &
                  0.124998690172 /)
    ! Interpolate first
    call spline_cubic_set ( 151, t(1:151), scalefunt(1:151), 1, 0._DP, 0, 0._DP, scalefunDD(1:151) )
  
!------------------------------------------------------------------------------
  CASE(15) ! GGA
    IF (lambda<-99.0) lambda = 1.0_DP
    IF (mu<-99._DP) mu = 1.0_DP
    
!------------------------------------------------------------------------------
  CASE(16) ! GGA + WGCD
    if(rho0==-1.d0) then
      rho0 = cell%numEle/cell%vol 
    ENDIF
    
    DO ii=1,numint
        d(ii) = 0.d0 + 0.01d0*(ii-1.0) 
    ENDDO
    
    ELFvsd = &
                (/0.00005000, 0.00021000, 0.00037127, 0.00054472, 0.00073318, & 
                  0.00093765, 0.00115868, 0.00139664, 0.00165186, 0.00192465, & 
                  0.00236011, 0.00309081, 0.00412453, 0.00520969, 0.00608653, & 
                  0.00804831, 0.00934243, 0.01138334, 0.01366375, 0.01558003, & 
                  0.01800965, 0.02043927, 0.02286888, 0.02529850, 0.02772812, & 
                  0.02988452, 0.03126635, 0.03314065, 0.03671603, 0.04029141, & 
                  0.04386678, 0.04744216, 0.05101753, 0.05459291, 0.05816829, & 
                  0.06174366, 0.06531904, 0.06957662, 0.07408569, 0.07884882, & 
                  0.08360986, 0.08869948, 0.09387569, 0.09905189, 0.10422809, & 
                  0.10940430, 0.11458050, 0.11975670, 0.12493291, 0.13010911, & 
                  0.13575203, 0.14156364, 0.14777025, 0.15409835, 0.16042646, & 
                  0.16675457, 0.17308268, 0.17941079, 0.18573890, 0.19254076, & 
                  0.19900258, 0.20546441, 0.21192623, 0.21838806, 0.22489239, & 
                  0.23158703, 0.23835919, 0.24628165, 0.25467238, 0.26306311, & 
                  0.27145384, 0.27984457, 0.28823530, 0.29603890, 0.30357624, & 
                  0.31111359, 0.31865093, 0.32618828, 0.33372562, 0.34126297, & 
                  0.34880031, 0.35633766, 0.36387501, 0.37141235, 0.37887925, & 
                  0.38634356, 0.39414054, 0.40193752, 0.40973450, 0.41753148, & 
                  0.42532847, 0.43312545, 0.43966414, 0.44604622, 0.45242831, & 
                  0.45938698, 0.46679639, 0.47420580, 0.48161520, 0.48902461, & 
                  0.49643402, 0.50384343, 0.51125283, 0.51866224, 0.52607165, & 
                  0.53348106, 0.54089046, 0.54829987, 0.55570928, 0.56311868, & 
                  0.57052809, 0.57793750, 0.58713237, 0.59276530, 0.59823903, & 
                  0.60428059, 0.61121037, 0.61413646, 0.61795201, 0.62262063, & 
                  0.63018503, 0.63820852, 0.64623202, 0.65214957, 0.65693808, & 
                  0.66172660, 0.66722826, 0.67265272, 0.67800015, 0.68327078, & 
                  0.68846488, 0.69358277, 0.69862482, 0.70359142, 0.70848303, & 
                  0.71330013, 0.71804324, 0.72271290, 0.72730971, 0.73183426, & 
                  0.73628720, 0.74066918, 0.74498089, 0.74922303, 0.75339632, & 
                  0.75750151, 0.76153934, 0.76551058, 0.76941602, 0.77325645, & 
                  0.77703267, 0.78074548, 0.78439571, 0.78798418, 0.79151172, & 
                  0.79497915, 0.79838731, 0.80173704, 0.80502917, 0.80826454, & 
                  0.81144399, 0.81456834, 0.81763844, 0.82065510, 0.82361917, & 
                  0.82653145, 0.82939276, 0.83220392, 0.83496573, 0.83767899, & 
                  0.84034450, 0.84296303, 0.84553538, 0.84806231, 0.85054458, & 
                  0.85298296, 0.85537819, 0.85773101, 0.86004215, 0.86231233, & 
                  0.86454227, 0.86673267, 0.86888422, 0.87099761, 0.87307352, & 
                  0.87511261, 0.87711554, 0.87908296, 0.88101551, 0.88291381, & 
                  0.88477848, 0.88661014, 0.88840938, 0.89017680, 0.89191297, & 
                  0.89361848, 0.89529387, 0.89693971, 0.89855653, 0.90014488/) 
    ! Interpolate first
    call spline_cubic_set ( numint, d(1:numint), ELFvsd(1:numint), 1, 0._DP, 0, 0._DP, ELFvsdDD(1:numint) )

!------------------------------------------------------------------------------
    CASE(17) ! MGP KEDF 
      ALLOCATE(keKernel(k1G, k2G, k3G,3), stat=fileStatus)
      IF (fileStatus/=0) THEN
         WRITE(message,*)'Error allocating the MGP kernel table. Leaving.'
         CALL Error(6, message)
      END IF
      IF (lambda<-99) lambda = 1.0_DP
      IF (mu<-99) mu = 1.0_DP
      WRITE(6,'(A,F12.6)') " Lambda for TF KEDF in MGP is ", lambda
      WRITE(6,'(A,F12.6)') " Mu     for vW KEDF in MGP is ", mu
      IF (rho0>0._DP) THEN
         hold0 = .TRUE.
      ELSE
         rho0 = cell%numEle/cell%vol
      END IF
      IF (gamma<0._DP) gamma = 1._DP
      WRITE(6,'(A,F12.6)') " gamma in MGP KEDF is ", gamma

      WRITE(6,'(A,F10.4)') " Original rho0 (Electron/Bohr^3) is ", rho0

      IF (.NOT. holdS) THEN
         rhoS=0._DP
      END IF

!------------------------------------------------------------------------------
  CASE DEFAULT ! By default we do nothing.
    
    WRITE(message,*) "SetupFunctional: Please check your type of KEDF, ", kinetic 
    CALL QUIT(message)

  END SELECT 

  RETURN

END SUBROUTINE SetupFunctional


SUBROUTINE KEDFRefresh(kinetic)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is intended to be run to refresh this module and make it 
!   ready for use each time the lattice vector is changed.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   We need to also refresh the value for rho0 here.
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell
  USE PlaneWave, ONLY: FillQTable
  USE KEDF_WTkernel, ONLY: FillWT
  USE KEDF_WGCkernel, ONLY: FillWGC
<<<<<<< HEAD
=======
  USE KEDF_MGPkernel, ONLY: FillMGP
>>>>>>> 86ef2760f4a74e31322f3edcbdbd6c06a2415ee6
  USE KEDF_CAT, ONLY: FillCAT
  USE Sys, ONLY: hold0, holdS, rho0, rhoS
  USE KEDF_WGCD, ONLY: mrhos

  IMPLICIT NONE
                     !>> EXTERNAL VARIABLES <<!

  INTEGER, INTENT(IN) :: kinetic
  REAL(KIND=DP) :: tmpRho


                     !>> INTERNAL VARIABLES <<! 
                     !>> INITIALIZATION <<!
  CALL Title("SetupKEDF::KEDFRefresh")
  CALL StartClock('KEDFRefresh')

  tmpRho = cell%numEle/cell%vol

                       !>> FUNCTION BODY <<!
  ! Run FillQTable with the correct reciprocal space lattice vectors
  ! pereorid boundary condition
  CALL FillQTable(cell%cellRecip)

  SELECT CASE(kinetic)
    ! 4: WT, 17: EVT 
    CASE(4,17)
      IF (.NOT.hold0) rho0 = tmpRho 
      CALL FillWT()
    ! 5: WGC, 12: WGCD, 18: EVC
    CASE(5,12,18)
      IF (.NOT.hold0) rho0 = tmpRho 
      IF (.NOT.holdS) rhoS = mrhos * tmpRho 
      CALL FillWGC()
<<<<<<< HEAD
=======
    ! 19: MGP
    CASE(19)
     IF (.NOT.hold0) rho0 = tmpRho
     CALL FillMGP()
>>>>>>> 86ef2760f4a74e31322f3edcbdbd6c06a2415ee6
    ! 10: CAT
    CASE(10)
      IF (.NOT.hold0) rho0 = tmpRho 
      IF (.NOT.holdS) rhoS = tmpRho 
      CALL FillCAT()
  END SELECT

  CALL StopClock('KEDFRefresh')

  RETURN

END SUBROUTINE KEDFRefresh

END MODULE SetupKEDF

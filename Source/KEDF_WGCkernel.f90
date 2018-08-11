MODULE KEDF_WGCkernel
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_WGCkernel 
!     |_SUBROUTINE FillWGC
!     |_SUBROUTINE FillWGC_ReciprocalSpace
!     |_SUBROUTINE FillWGC_RestartFrom_ReciprocalSpace 
!
! DESCRIPTION:
!   Setup the WGC KEDF kernels. 
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
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP, pi 
  USE FOURIER, ONLY: FFT
  USE PlaneWave, ONLY: qTable
  USE KEDF_WTkernel, ONLY: alpha, beta
  USE KEDF_WTkernel, ONLY: keKernel
  USE SYS, ONLY: rho0, rhoS

  USE MPI_Functions
  USE OutputFiles
  USE OUTPUT, ONLY: WrtOut

  IMPLICIT NONE

                    !>> INTERNAL VARIABLES <<!  

  REAL(KIND=DP) :: gamma = -1.0_DP    
  ! Averaging parameter for the WGC functional.
  !
  INTEGER :: firstOrderWGC= -1
  ! default value is -1 corresponding to WGCT==2 
  !
  INTEGER :: WGCT = 2             
  ! by default WGC will be expanded to the 2nd order
  !
  REAL(KIND=DP), ALLOCATABLE :: nls_wpp(:)   
  ! used for kernel interploation
  !
  REAL(KIND=DP), ALLOCATABLE :: nls_w1pp(:)  
  ! lifetime is during the code lifetime
  !
  REAL(KIND=DP), ALLOCATABLE :: nls_w2pp(:)
  !

  REAL(KIND=DP) :: alpha5 = -100.0_DP 
  ! Coefficient to diagonal terms of the kernel in WGC
  !
  REAL(KIND=DP) :: beta5 = -100.0_DP  
  ! Coefficient to diagonal terms of the kernel in WGC
  !
  INTEGER :: i, j, k, ix, i2, i3, in, it, status
  !
  REAL(KIND=DP) :: eta, tkFstar
  !
  LOGICAL :: have_RK = .false. 
  ! flag for if the Runge-Kutta method has been called
  !

CONTAINS


SUBROUTINE FillWGC
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
  
  ! >> INITIALIZATION << !!
  WRITE(outputUnit, '(/A)') " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')  "                           INITIALIZING WGC KERNEL"

  ! tkFstar = 2 k_{F}
  tkFstar = 2._DP * (3._DP * pi**2 * rhoS)**(1._DP/3._DP)
  WRITE(outputUnit,*) "Fermi vector : ", tkFstar/2.0_DP
  WRITE(outputUnit,*) "Reference density : ", rhoS

  !-------------------------------------------------------
  ! Periodic boundary conditions
  CALL FillWGC_ReciprocalSpace()  

  RETURN

END SUBROUTINE FillWGC


SUBROUTINE FillWGC_ReciprocalSpace()
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

  USE IntKernelODE, ONLY: &  ! this module used rksuite to integrate kernel's ODE equ.a
    KEDFType,  & 
    winf,  &                   ! kernel at infinit, starting value for RK method
    int_eta => eta , &         ! kernel is defined on eta
    int_w  => w,  &            ! the kernel
    int_w1 => w1, &            ! d kernel/d eta
    int_w2 => w2 ,&            ! d^2 kernel /d eta^2
    ode_alpha, &               ! the alpha in David's paper
    ode_beta,  &               ! the beta
    ode_gamma, &               ! the gamma
    makeKernel2, &             ! the driver to do ODE integration
    IntKernelODEClean => Clean ! clean the all the temperory vars

  USE MathSplines, ONLY: &
    spline_cubic_set, &
    spline_cubic_val

  USE OUTPUT, ONLY : outputKernel

  IMPLICIT NONE

  REAL(KIND=DP) :: ypval, yppval
  REAL(KIND=DP) :: zero = 0.0 ! is needed in spline_cubic_set, we can not use 0.0,
  INTEGER :: iq
! because it is real. 
!  ! debug
!  ! the following code is used to solve WGC differential equation with RK method
!  ! here we commnet it off and still use the original code by Vincent

    ! outputKernel can be set in the input file.
    IF(outputKernel .AND. rankGlobal==0) THEN
      OPEN (unit=1981, file='KEDF-kernel-G.dat', status = 'unknown', form='formatted', action='write')
      WRITE(outputUnit,*) "This is the KEDF kernels with parameters: "
      WRITE(1981,'(A)') "WGC Kernels in Gspace: eta, w(eta), w'(eta), w''(eta) "
      WRITE(1981,'(ES20.12, A)') tkFstar, " Fermi Vector"
      WRITE(1981,'(ES20.12, A)') rhoS, " Reference Density"
    ENDIF

    ! RK method is only used in the first time.
    IF (.not. have_RK) THEN 

      WRITE(outputUnit,*) "Make WGC kernels."

      !! Setting parameters for doing Runge-Kutta method
      ode_alpha = alpha
      ode_beta  = beta 
      ode_gamma = gamma
      wInf      = -1.6d0 * 20.d0/ (36.d0*alpha*beta)  ! w(eta) at eta=infinity
                                                    ! from eq(26) in WGC paper (1999,PRB)
      KEDFtype  = 2  ! 1: for CAT KEDF,
                     ! 2: for WGC
      WRITE(message,*) " "; CALL WrtOut(6,message)
      WRITE(message,*) '(FillWGC) Make WGC kernel';   CALL WrtOut(6,message)
      WRITE(message,*) '(FillWGC) Alpha                         : ', ode_alpha;   CALL WrtOut(6,message)
      WRITE(message,*) '(FillWGC) Beta                          : ', ode_beta;    CALL WrtOut(6,message)
      WRITE(message,*) '(FillWGC) Gamma                         : ', ode_gamma;   CALL WrtOut(6,message)
    
      ! Perform Runge-Kutta method 
      CALL makeKernel2
      int_w1 = int_eta*int_w1        ! not as before, makeKernel will give int_w1 => w'
      int_w2 = int_eta**2 * int_w2   !                int_w2 =>  w''
    
!      WRITE(message,*) ' Number of eta(:) = ', SIZE(int_eta)
!      CALL WrtOut(6,message)
      CALL WrtOut(6,message)

      ALLOCATE(nls_wpp(SIZE(int_eta)), nls_w1pp(SIZE(int_eta)), nls_w2pp(SIZE(int_eta)))

      CALL spline_cubic_set ( SIZE(int_eta), int_eta, int_w , 0,zero,0,zero, nls_wpp )
      CALL spline_cubic_set ( SIZE(int_eta), int_eta, int_w1, 0,zero,0,zero, nls_w1pp)
      CALL spline_cubic_set ( SIZE(int_eta), int_eta, int_w2, 0,zero,0,zero, nls_w2pp)
      

      ! Print out the one dimension kernel in G space,
      IF(outputKernel .AND. rankGlobal==0) THEN
          WRITE(1981,'(I8,A)') SIZE(int_eta), " Size of eta."
          DO iq=1, SIZE(int_eta)
            WRITE(1981,'(4ES20.12)') int_eta(iq), int_w(iq), int_w1(iq), int_w2(iq)
          ENDDO
      ENDIF


      have_RK = .TRUE.
    ENDIF  !! have_RK 

    ! And now we can compute the kernel and its derivatives at every q-point.

!   CALL WrtOut(6,' (FillWGC) Interpolating WGC kernels')
    flush(outputUnit)
    keKernel = 0._DP

    WRITE(outputUnit, *) "Dimension of KEDF Kernel : ", &
      SIZE(keKernel,1), SIZE(keKernel,2), SIZE(keKernel,3)

!    WRITE(outputUnit, *) "size of int_eta = ", SIZE(int_eta)

    DO i3=1, SIZE(keKernel,3)
      DO i2=1, SIZE(keKernel,2)
        DO ix=1, SIZE(keKernel,1)

          ! eta = |q| / (2*kF)
          eta = qTable(ix,i2,i3) / tkFstar

          ! Get the kernel value, 1st and 2nd derivitives as a function of
          ! eta.
          ! Temporarily, here keKernel(:,:,:,1) <= w(eta), 
          !                   keKernel(:,:,:,2) <= eta * w'(eta)
          !                   keKernel(:,:,:,3) <= eta**2 * w'(eta)

          CALL spline_cubic_val ( SIZE(int_eta), int_eta, int_w,   nls_wpp, eta, kekernel(ix,i2,i3,1), ypval, yppval )
          CALL spline_cubic_val ( SIZE(int_eta), int_eta, int_w1, nls_w1pp, eta, kekernel(ix,i2,i3,2), ypval, yppval )
          CALL spline_cubic_val ( SIZE(int_eta), int_eta, int_w2, nls_w2pp, eta, kekernel(ix,i2,i3,3), ypval, yppval )

          !IF(outputKernel .AND. rankGlobal==0) THEN
          !  WRITE(1981, '(4ES20.12)') eta,kekernel(ix,i2,i3,1),kekernel(ix,i2,i3,2),kekernel(ix,i2,i3,3)
          !ENDIF

          keKernel(ix,i2,i3,4) = (keKernel(ix,i2,i3,3) &
            + (1+gamma)*keKernel(ix,i2,i3,2)) &
            / (36._DP * rhoS**2)
          keKernel(ix,i2,i3,3) = (keKernel(ix,i2,i3,3) &
            + (7-gamma)*keKernel(ix,i2,i3,2)) &
            / (36._DP * rhoS**2)
          keKernel(ix,i2,i3,2) = -keKernel(ix,i2,i3,2) / (6._DP * rhoS)

        END DO ! ix
      END DO ! i2
    ENDDO ! i3
    
    !Post-Process keKernel based on input keyword WGCT
    IF (WGCT==0) THEN
        !only 0th order WGC kernel will be used
        keKernel(:,:,:,2:4) = 0._DP
    ELSE IF (WGCT==1) THEN
        keKernel(:,:,:,3:4) = 0._DP
    ELSE IF (WGCT==2) THEN
    ELSE IF (WGCT==-3) THEN
        keKernel(:,:,:,3) = 0._DP  ! ddW/ddRho(r) = 0
    ELSE IF (WGCT==-4) THEN
        keKernel(:,:,:,4) = 0._DP  ! ddW/dRho(r)dRho(r')=0
    ELSE IF (WGCT==-5) THEN        ! for WGCT -5, alpha5*ddW/ddRho(r) + beta5*ddW/dRho(r)dRho(r')
        keKernel(:,:,:,3) = alpha5*keKernel(:,:,:,3)
        keKernel(:,:,:,4) = beta5*keKernel(:,:,:,4)
    ELSE 
        CALL WrtOut(6, ' (FillWGC) WGCT is set to a value which is not defined, error, stop code')
        STOP
    ENDIF

    IF(outputKernel .AND. rankGlobal==0) THEN
     CLOSE(1981)
     outputKernel=.FALSE. ! only print once
    ENDIF

    RETURN

END SUBROUTINE FillWGC_ReciprocalSpace


SUBROUTINE FillWGC_RestartFrom_ReciprocalSpace(kernelFile, kernel, qNorm)
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

  USE MATHSPLINES, ONLY : spline_cubic_val

  IMPLICIT NONE

  ! >> INPUT VARIABLES << !
  CHARACTER(LEN=*), INTENT(IN) :: kernelFile
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: qNorm
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: kernel

  ! >> LOCAL VARIABLES << !
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: eta_in
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: w0_in, w1_in, w2_in
  INTEGER :: nq
  REAL(KIND=DP) :: ypval, yppval
  INTEGER :: i3, i2, ix
  INTEGER :: iq
  INTEGER :: inputUnit
  INTEGER :: fileStatus

  CHARACTER(LEN=5) :: KEDF_type
  REAL(KIND=DP) :: Fermiv
  REAL(KIND=DP) :: RhoS_in

  WRITE(outputUnit,'(A)') " Read in WGC kernels."

  ! open the kernel file
  OPEN(UNIT=inputUnit, ACTION="read", BLANK="null", &
       FILE=kernelFile, FORM="formatted",&
       IOSTAT=fileStatus, PAD="no", STATUS="old")
  IF (fileStatus/=0) THEN
    WRITE(errorUnit,*)'Could not open ',TRIM(kernelFile),'.'
    WRITE(errorUnit,*)'Please make sure this density file does exist.'
    STOP
  ELSE IF(fileStatus==0) THEN

    READ (inputUnit,*) KEDF_type
    READ (inputUnit,*) Fermiv
    READ (inputUnit,*) RhoS_in
    READ (inputUnit,*) nq

    WRITE(outputUnit,'(A,A)')       " Kernel file name  : ", kernelFile
    WRITE(outputUnit,'(A,A)')       " KEDF type         : ", KEDF_type
    WRITE(outputUnit,'(A,ES20.12)') " Fermi vector      : ", Fermiv
    WRITE(outputUnit,'(A,ES20.12)') " Reference density : ", RhoS_in
    WRITE(outputUnit,'(A,I8)')      " Mesh points       : ", nq

  END IF

  ! read in kernels.
  ALLOCATE(eta_in(nq))
  ALLOCATE(w0_in(nq))
  ALLOCATE(w1_in(nq))
  ALLOCATE(w2_in(nq))

  DO iq=1, nq
    READ (inputUnit,*) eta_in(iq), w0_in(iq), w1_in(iq), w2_in(iq)
  ENDDO


  !! >>> FUNCTION << !!!
  DO i3=1, SIZE(kernel,3)
    DO i2=1, SIZE(kernel,2)
      DO ix=1, SIZE(kernel,1)

        ! eta = |q| / (2*kF)
        eta = qNorm(ix,i2,i3) / Fermiv

        ! Get the kernel value, 1st and 2nd derivitives as a function of
        ! eta.
        ! Temporarily, here kernel(:,:,:,1) <= w(eta), 
        !                   kernel(:,:,:,2) <= eta * w'(eta)
        !                   kernel(:,:,:,3) <= eta**2 * w'(eta)

        CALL spline_cubic_val ( nq, eta_in, w0_in,  nls_wpp, eta, kernel(ix,i2,i3,1), ypval, yppval )
        CALL spline_cubic_val ( nq, eta_in, w1_in, nls_w1pp, eta, kernel(ix,i2,i3,2), ypval, yppval )
        CALL spline_cubic_val ( nq, eta_in, w2_in, nls_w2pp, eta, kernel(ix,i2,i3,3), ypval, yppval )

        !IF(outputKernel .AND. rankGlobal==0) THEN
        !  WRITE(1981, '(4ES20.12)') eta,kernel(ix,i2,i3,1),kernel(ix,i2,i3,2),kernel(ix,i2,i3,3)
        !ENDIF

        kernel(ix,i2,i3,4) = (kernel(ix,i2,i3,3) &
            + (1+gamma)*kernel(ix,i2,i3,2)) &
            / (36._DP * rhoS**2)
        kernel(ix,i2,i3,3) = (kernel(ix,i2,i3,3) &
            + (7-gamma)*kernel(ix,i2,i3,2)) &
            / (36._DP * rhoS**2)
        kernel(ix,i2,i3,2) = -kernel(ix,i2,i3,2) / (6._DP * rhoS)

      END DO ! ix
    END DO ! i2
  ENDDO ! i3



    !Post-Process kernel based on input keyword WGCT
    IF (WGCT==0) THEN
        !only 0th order WGC kernel will be used
        kernel(:,:,:,2:4) = 0._DP
    ELSE IF (WGCT==1) THEN
        kernel(:,:,:,3:4) = 0._DP
    ELSE IF (WGCT==2) THEN
    ELSE IF (WGCT==-3) THEN
        kernel(:,:,:,3) = 0._DP  ! ddW/ddRho(r) = 0
    ELSE IF (WGCT==-4) THEN
        kernel(:,:,:,4) = 0._DP  ! ddW/dRho(r)dRho(r')=0
    ELSE IF (WGCT==-5) THEN        ! for WGCT -5, alpha5*ddW/ddRho(r) + beta5*ddW/dRho(r)dRho(r')
        kernel(:,:,:,3) = alpha5*kernel(:,:,:,3)
        kernel(:,:,:,4) = beta5*kernel(:,:,:,4)
    ELSE 
        CALL WrtOut(6, ' (FillWGC) WGCT is set to a value which is not defined, error, stop code')
        STOP
    ENDIF
    
  CLOSE (inputUnit, STATUS="keep")

  RETURN

END SUBROUTINE FillWGC_RestartFrom_ReciprocalSpace


END MODULE KEDF_WGCkernel

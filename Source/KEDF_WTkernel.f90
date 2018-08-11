MODULE KEDF_WTkernel 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_WTkernel 
!     |_SUBROUTINE FillWT
!     |_SUBROUTINE FillWT_ReciprocalSpace
!     |_SUBROUTINE FillWT_RestartFrom_GSpace
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
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP, PI 
  USE FOURIER, ONLY: FFT
  USE SYS, ONLY: rho0
  USE SYS, ONLY: rhoS
  USE MathFunctions, ONLY: LindG
  USE PlaneWave, ONLY: qTable 
  USE MPI_Functions
  USE OutputFiles
  USE OUTPUT, ONLY: WrtOut

  IMPLICIT NONE

                    !>> INTERNAL VARIABLES <<!  

  REAL(KIND=DP) :: alpha = -100.0_DP  
  ! Exponent to rho left of the kernel in WT and WGC.
  !
  REAL(KIND=DP) :: beta = -100.0_DP   
  ! Exponent to rho on the right side of the kernel in WT and WGC.
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: keKernel 

  ! When WT and WGC w/ periodic b.c., stored in recip space.

  ! rinc is emperically determined. 0.1 seems to be
  ! converged w.r.t. 0.01 for a single atom in vacuum by 
  ! .004 meV/atom. This is no guarantee for other systems, though.
  REAL(KIND=DP), PARAMETER :: rinc = 0.1_DP
  REAL(KIND=DP), PARAMETER :: maxr4K2 = 500._DP        
  ! Past this number, we will use the analytical formula for K_II.
  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: wtKernelr
  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: wtKernelFr

  INTEGER :: &
    allocateStatus, &
    ix, i2, i3, &          ! Dummy counters
    i,j,k, &
    noPadX, noPadY, noPadZ, &  ! Size of unpadded array
    size1DKernel             ! Size of the wt/wgc 1D kernel

  REAL(KIND=DP) :: ft = 5.0_DP / 3.0_DP  ! five thirds
  REAL(KIND=DP) :: coef                  ! Lindhard function multiplier.
  REAL(KIND=DP) :: tkF                   ! 2 * kF [1] (10)
  REAL(KIND=DP) :: x2, y2, z2 ! ???
  REAL(KIND=DP) :: scaling    ! ???
  REAL(KIND=DP) :: qMax
  REAL(KIND=DP) :: lowerLim
  REAL(KIND=DP) :: upperLim
  REAL(KIND=DP) :: addValue
  REAL(KIND=DP) :: rNorm
  REAL(KIND=DP) :: k2Value

CONTAINS



SUBROUTINE FillWT()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine initializes the kernel for the WT kinetic energy
!   functional. It uses the table of norms qTable and the average electron 
!   density rho0 to output a 4-D table that contains the corresponding kernel 
!   at every q-point. This subroutine should be executed every time the cell 
!   shape is altered, after the qTable is reset.
!
!   After some trial and error, I arbitrarily decided that 1000 is an
!   appropriate place to cut off the WT kernel in real space.
!
!   Note that for periodic boundary conditions, the kernel is saved in
!   reciprocal space, to avoid having to recompute it when calculating energy.
!   For Dirichlet boundary conditions, however, the kernel is saved in real
!   space before padding - this cuts down on the amount of memory needed.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Wang Y.A. Govind N. and Carter E.A., Phys. Rev. B 60(24), 1999, 16350
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/05/2003  Created (Vincent Ligneres)
!   12/1/2005   Modified to do kernel in real space
!
!------------------------------------------------------------------------------

  IMPLICIT NONE


                      !>> INITIALIZATION <<!

  coef = 5._DP/(9._DP*alpha*beta*rho0**(alpha+beta-ft))    ! See [1] eq. (16)
  !WRITE(outputUnit,*) " WT coef, alpha=", alpha
  !WRITE(outputUnit,*) " WT coef, beta=", beta
  !WRITE(outputUnit,*) " WT coef, ft=", ft
  !WRITE(outputUnit,*) " WT coef, rho0=", rho0

  ! two * Kf
  tkF = 2._DP * (3._DP * rho0 * pi**2)**(1._DP/3._DP) 
!  coef2 = -tKf**2 * 0.13714285714285712_DP          ! Coefficient for 
                                                    ! limit as as q->in

                       !>> FUNCTION BODY << 

  ! This is the usual case with the wt kernel in reciprocal space
  ! Periodic boundary condition
  CALL FillWT_ReciprocalSpace()

  RETURN

END SUBROUTINE FillWT


SUBROUTINE FillWT_ReciprocalSpace()
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

  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu
  USE OUTPUT, ONLY: outputKernel

  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!  

  INTEGER :: n_eta 
  ! how many eta points I would like to print
  !
  INTEGER :: ie    
  ! the counter
  !
  REAL(KIND=DP) :: d_eta 
  ! delta eta
  !

  ! outputKernel can be set in the input file.
  ! only the main processor is allowed to write.
  IF(outputKernel .AND. rankGlobal==0) THEN

    OPEN (unit=1985, file='KEDF-kernel-G.dat', status = 'unknown', form='formatted', action='write')
    WRITE(outputUnit,*) "Print the G space Wang-Teter KEDF kernels into KEDF-kernel-G.dat "
    WRITE(1985,'(A)') "4"
    WRITE(1985,'(ES22.15, A)') tkF
    WRITE(1985,'(ES22.15, A)') rhoS

    ! set the value here, generally this is large enough for
    ! restart calculations
    n_eta = 50000
    d_eta = 0.001

    WRITE(1985,'(I8, A)') n_eta

    DO ie=1, n_eta
      WRITE(1985,'(2ES20.12)') ie*d_eta, LindG(ie*d_eta,lambda,mu)
    ENDDO

    CLOSE(1985)
    outputKernel=.FALSE.

  ENDIF


  WRITE(outputUnit,*) "Fill the 3D Wang-Teter kernel in G space"
  ! fill the Wang-Teter kernel
  DO i3=1, SIZE(keKernel,3)
    DO i2=1, SIZE(keKernel,2)
      DO ix=1, SIZE(keKernel,1)
        keKernel(ix,i2,i3,1) = LindG(qTable(ix,i2,i3)/tkF,lambda,mu) * coef 
      END DO ! ix
    END DO ! i2
  ENDDO ! i3

  RETURN

END SUBROUTINE FillWT_ReciprocalSpace



SUBROUTINE FillWT_RestartFrom_GSpace(kernelFile, kernel, qNorm)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 05/28/2013 Created by Mohan Chen
!------------------------------------------------------------------------------

  !USE MATHFUNCTIONS, ONLY : spline_cubic_val

  IMPLICIT NONE

  ! >> INPUT VARIABLES << !
  CHARACTER(LEN=*), INTENT(IN) :: kernelFile ! name of the kernel file
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: qNorm ! norm of q vectors
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: kernel ! the target kernel
  !array we would like to fill

  ! >> LOCAL VARIABLES << !
  !---------------
  ! read in info 
  !---------------
  INTEGER :: KEDF_type ! the type of this KEDF kernel file
  REAL(KIND=DP) :: kF_in         ! read in fermi vector
  REAL(KIND=DP) :: RhoS_in      ! read in reference density
  REAL(KIND=DP) :: maxEta_in     ! maximal eta read in
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: eta_in ! read in eta
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: w0_in ! read in Wang-Teter kernel 
  INTEGER :: nq_in ! read in mesh points for Wang-Teter kernel

  REAL(KIND=DP) :: maxG         ! maximal G
  REAL(KIND=DP) :: maxEta       ! maximal eta

  !---------------
  ! file related 
  !---------------
  INTEGER :: inputUnit  ! kernel file pointer 
  INTEGER :: fileStatus ! check the status of kernel file

  !-------------
  ! counters
  !-------------
  INTEGER :: iq         ! used for read in kernel file
  INTEGER :: i3, i2, ix ! used with 'kernel'
  REAL(KIND=DP) :: eta  ! eta is calculated from qNorm

  !---------------
  ! interpolation
  !---------------
  INTEGER :: ind1, ind2, ind_mid ! for new interpolation
  REAL(KIND=DP) :: fac1, fac2    ! for new interpolation
  ! REAL(KIND=DP) :: ypval, yppval ! for old interpolation method


  CALL Title("KEDF_WTkernel:FillWT_RestartFrom_GSpace")
  WRITE(outputUnit,'(A)') " Read in WT kernel."

  !-------------------------
  ! open the kernel file
  !-------------------------
  OPEN(UNIT=inputUnit, ACTION="read", BLANK="null", &
       FILE=kernelFile, FORM="formatted",&
       IOSTAT=fileStatus, PAD="no", STATUS="old")

  !---------------------------
  ! if the file doesn't exist
  !---------------------------
  IF (fileStatus/=0) THEN

    WRITE(errorUnit,*)'Could not open ',TRIM(kernelFile),'.'
    WRITE(errorUnit,*)'Please make sure this density file does exist.'
    STOP

  !---------------------------
  ! else if the file exists
  !---------------------------
  ELSE IF(fileStatus==0) THEN

    READ (inputUnit,*) KEDF_type  ! Type should be 'WT'. 
    READ (inputUnit,*) kF_in      ! Fermi vector value.
    READ (inputUnit,*) RhoS_in    ! Read in reference density.
    READ (inputUnit,*) nq_in      ! Number of 1D q points.

    ! be careful of the format, that might cause segment fault!!!
    WRITE(outputUnit,'(A,A)')       " kernel file name  : ", kernelFile
    WRITE(outputUnit,'(A,I5)')      " KEDF type         : ", KEDF_type
    WRITE(outputUnit,'(A,ES20.12)') " Fermi vector      : ", kF_in
    WRITE(outputUnit,'(A,ES20.12)') " Reference density : ", RhoS_in
    WRITE(outputUnit,'(A,I8)')      " Mesh points       : ", nq_in

  END IF

 
  ! allocate the 1D kernel
  ALLOCATE(eta_in(1:nq_in))
  ALLOCATE(w0_in(1:nq_in))

  ! read in one dimension eta=q/2kF, and w0_in is the value of kernel
    DO iq=1, nq_in
    READ (inputUnit,*) eta_in(iq), w0_in(iq)
    ! then this is the real linear response with proper factor
    w0_in(iq) = w0_in(iq) * coef
    ! in fact, read in is G, not G/kF_in
    ! WRITE(outputUnit,*) eta_in(iq), " ", w0_in(iq)
  ENDDO

  ! save the maximal eta value
  maxEta_in = eta_in(nq_in)


  IF( tkF .NE. kF_in ) THEN
    WRITE(outputUnit,*) "Fermi vector in this cell is ",tkF
    WRITE(outputUnit,*) "Fermi vector from kernel file is ",kF_in
  ENDIF


  !! >>> FUNCTION << !!
  maxG = 0.D0
  maxEta = 0.D0
  WRITE(outputUnit,*) "First qNorm is ", qNorm(1,1,1)
  WRITE(outputUnit,*) "size of Wang-Teter kernel in this run is ", SIZE(kernel,1), SIZE(kernel,2), SIZE(kernel,3)
  DO i3=1, SIZE(kernel,3)
    DO i2=1, SIZE(kernel,2)
      DO ix=1, SIZE(kernel,1)

        ! eta = |q| / (2*kF)
        eta = qNorm(ix,i2,i3) / tkF
        maxG = max(qNorm(ix,i2,i3), maxG)
        maxEta = max(eta, maxEta)

        !WRITE(outputUnit,*) "Kernel Old = ", kernel(ix,i2,i3,1)
        !----------------------------------------------------
        ! do interpolation here !!!
        ! if the eta is too small, we should use
        ! the first value in the kernel file directly
        !----------------------------------------------------
        IF( eta <= eta_in(1) ) THEN
          kernel(ix,i2,i3,1) = w0_in(1)
        ELSE IF( eta > maxEta_in ) THEN
          kernel(ix,i2,i3,1) = w0_in(nq_in)
        ELSE
          ind1 = 1
          ind2 = nq_in
          DO WHILE (ind1 < ind2-1)
            ind_mid = (ind1 + ind2)/2
            IF(eta > eta_in(ind_mid) ) THEN
              ind1 = ind_mid
            ELSE
              ind2 = ind_mid
            ENDIF
          ENDDO
          fac1 = (eta_in(ind2)-eta)/(eta_in(ind2)-eta_in(ind1))
          fac2 = (eta-eta_in(ind1))/(eta_in(ind2)-eta_in(ind1))
          kernel(ix,i2,i3,1) = fac1*w0_in(ind1) + fac2*w0_in(ind2)
        ENDIF

!---------------------------------------------
! old spline method, not accurate enough
! mohan 2013-07-17
!---------------------------------------------
!        CALL spline_cubic_val ( nq_in, eta_in, w0_in,  w1_in, eta, kernel(ix,i2,i3,1), ypval, yppval )
        !WRITE(outputUnit,*) "Kernel New = ", kernel(ix,i2,i3,1)


!---------------------------------------
! for test
!---------------------------------------
        !IF( DABS( kernel(ix, i2, i3, 1) ) > 1.0e-5 ) THEN
        !WRITE(outputUnit,'(3I5,ES20.12)') ix,i2,i3,kernel(ix,i2,i3,1)
        !ENDIF


      END DO ! ix
    END DO ! i2
  ENDDO ! i3


  WRITE(outputUnit,'(A,ES20.12)') " Max |G| in this calculation is ", maxG
  WRITE(outputUnit,'(A,ES20.12)') " Max eta used from kernel is ", maxEta

  IF(maxEta_in .LT. maxEta) THEN
    WRITE(outputUnit,*) " WARNING! Please increase the maximal eta value in KEDF kernel file."
    STOP
  ENDIF

  DEALLOCATE(eta_in)
  DEALLOCATE(w0_in)
  CLOSE (inputUnit) 

  RETURN

END SUBROUTINE FillWT_RestartFrom_GSpace


END MODULE KEDF_WTkernel 

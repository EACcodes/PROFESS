MODULE Fourier_new
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!    MODULE Fourier_new
!       |_SUBROUTINE PlanFFT
!       |_SUBROUTINE GetFFTDims
!       |_SUBROUTINE GetFFTComplexDims
!       |_INTERFACE FFT
!         |_FUNCTION ForwardFFT_4D (Private)
!         |_FUNCTION BackFFT_4D (Private)
!         |_FUNCTION ForwardFFT_3D (Private)
!         |_FUNCTION BackFFT_3D (Private)
!       |_SUBROUTINE CleanFFT
!
! DESCRIPTION:
!   This module interfaces with the Fastest Fourier Transform in the World 
!   (FFTW) public library v3.X to provide Fourier transform facilities 
!   for our quantities. Each Fourier transform has to be planned for first 
!   (PlanFFT) then executed as many times as necessary (FFT() ) and finally 
!   cleaned up (free up the memory) using CleanFFT.
!
! REFERENCES:
!   You are encouraged to consult the FFTW3 manual online at 
!   http://www.fftw.org
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!   12/15/2003  Changed INTEGER*8 to INTEGER(KIND=8) to make the compiler happy
!               Also reformatted a bunch of stuff, made blurbs (GSH)
!   08/03/2013 Finally porting this over to FFTW3, just 14 years after the last
!              commit to FFTW2 (JMD)
!   08/18/2013 Abstracting things, removing global state, also I removed the
!              individual revision logs as they were redundant with this (JMD)
!   08/26/2013 Finalizing the serial version (1.0) of this. Still some issues
!              outstanding that I'd like to address though. (JMD)
!   09/13/2013 Change the normalization (make it more stable against floating
!              point errors, start working on wisdom import/export.
!              Initial branch off the serial code for a parallel code. This
!              will merge in stuff from the old parallel interface if needed.
!              (JMD)
!   10/11/2013 Preparing for merge into Profess using a shim. Please note that
!              I think this is almost the worst possible option (only second to
!              leaving things as they are now) but this is how things are. (JMD) 
!------------------------------------------------------------------------------
                             ! << GLOBAL >> !
USE CONSTANTS, ONLY: DP      ! A double precision number
USE,INTRINSIC::iso_c_binding

IMPLICIT NONE

include 'fftw3-mpi.f03'
#ifdef __USE_PARALLEL
INCLUDE 'mpif.h'
#endif

INTEGER,PARAMETER :: &
     FFTW3_ONEPLAN_ESTIMATE = 1, &
     FFTW3_ONEPLAN_MEASURE = 2, &
     FFTW3_ONEPLAN_PATIENT = 3


TYPE FFT_CONFIG

   INTEGER :: fftAlgo        ! the mode for the FFT

   INTEGER(C_INTPTR_T):: &
     totalDimX, &            ! dimX
     totalDimY, &            ! dimY, total
     totalDimZ, &            ! dimZ, total
     localDimZ, &         ! the local part of the z dimension
     localDimZOff

   INTEGER :: &
     iCountFFT = 0        ! FFT counter for reporting purposes.
     
#ifdef FFT_NORM_EXACT
   INTEGER(KIND=DP) :: normConst
#else
   REAL(kind=DP) :: normConst
#endif
   
   TYPE(C_PTR)::fftw3GlobFWD,fftw3GlobBWD
   TYPE(C_PTR)::cdata
   REAL(C_DOUBLE),POINTER :: in(:,:,:) ! a pointer into dat representing a real scratch field
   COMPLEX(C_DOUBLE_COMPLEX), POINTER:: out(:,:,:) ! pointer for complex scratch
   
END TYPE FFT_CONFIG

PRIVATE :: &
  ForwardFFT_4D, &        ! Use FFT
  BackFFT_4D, &           ! Use FFT
  ForwardFFT_3D, &        ! Use FFT
  BackFFT_3D              ! Use FFT

! This interface picks the right transform to perform based on the nature of
! the incomming array: if it's real the FFT is done forward, if complex the 
! back transform is done. All the calls in OFDFT should be of this type: 
! FFT(config, f).
INTERFACE FFT_NEW
  MODULE PROCEDURE ForwardFFT_4D
  MODULE PROCEDURE BackFFT_4D
  MODULE PROCEDURE ForwardFFT_3D
  MODULE PROCEDURE BackFFT_3D
END INTERFACE

TYPE(FFT_CONFIG),SAVE::FFT_STD_STATE

CONTAINS

SUBROUTINE PlanFFT_NEW(MPI_COMM,config,mode,dimX,dimY,dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is the initialization procedure that sets up the FFT and its derived
!   data type for a FFT in the above specified dimensions.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(IN) :: &
    dimX, dimY, dimZ, &     ! The dimensions of the cell to be FFT'd
    mode
  TYPE(FFT_CONFIG),INTENT(INOUT)::config  
  INTEGER,INTENT(IN):: MPI_COMM                    ! communicator
  
  INTEGER(C_SIZE_T)::totalLocalSize
  INTEGER(C_INT)::planDepth
  INTEGER(C_INTPTR_T)::x,y,z

                       !>> INTERNAL VARIABLES <<! 
#ifndef FFT_NORM_EXACT
  INTEGER(KIND=8)::divNorm
#endif
                         !>> INITIALIZATION <<! 
                         
  config%fftAlgo = mode 

  config%totalDimX = dimX
  config%totalDimY = dimY
  config%totalDimZ = dimZ
  x = dimX
  y = dimY
  z = dimZ
  
#ifdef FFT_NORM_EXACT
  config%normConst = INT(dimX,kind=DP)*INT(dimY,kind=DP)*INT(dimZ,kind=DP)
#else
  divNorm = INT(dimX,KIND=8)*INT(dimY,KIND=8)*INT(dimZ,KIND=8)
  config%normConst = (REAL(1,KIND=DP))/ (REAL(divNorm,KIND=DP))
#endif
  

#ifdef DEBUG_FFT_INTER
  write(*,*) "DEBUG: Sizes before calling into FFTW3"
  write(*,*) "DEBUG: Total dim X ",config%totalDimX
  write(*,*) "DEBUG: Total dim Y ",config%totalDimY
  write(*,*) "DEBUG: Total dim Z ",config%totalDimZ
  write(*,*) "DEBUG: Norming: ",config%normConst
  write(*,*) "DEBUG: Enering FFT local call with ",z,y,x/2+1
  flush(6)
#endif 
 
  ! figure out the local size
  ! please note that FFTW sometimes needs extra memory and fiddles around with
  ! the dimensions in that case. So we hand a copy over.
  totalLocalSize = fftw_mpi_local_size_3d(z, y, x/2+1, &
                         MPI_COMM, config%localDimZ,config%localDimZOff)
                                                  
  ! allcate our scratch space (in-place FFT!)
  config%cdata = fftw_alloc_complex(totalLocalSize)
  ! alias a little around... ;-)
  call c_f_pointer(config%cdata, config%in, [2*(config%totalDimX/2+1),config%totalDimY,config%localDimZ])
!  call c_f_pointer(config%cdata, config%in,[config%totalDimX,config%totalDimY,config%localDimZ])
  call c_f_pointer(config%cdata, config%out, [config%totalDimX/2+1,config%totalDimY,config%localDimZ])

#ifdef DEBUG_FFT_INTER
  write(*,*) "DEBUG: Total local memory for FFT: ",totalLocalSize
  write(*,*) "DEBUG: Total dim X ",config%totalDimX
  write(*,*) "DEBUG: Total dim Y ",config%totalDimY
  write(*,*) "DEBUG: Total dim Z ",config%totalDimZ
  write(*,*) "DEBUG: Local dim Z ",config%localDimZ
  write(*,*) "DEBUG: Local offset Z ",config%localDimZOff
  write(*,*) "DEBUG: Size X real scratch pointer ",size(config%in,1)
  write(*,*) "DEBUG: Size Y real scratch pointer ",size(config%in,2)
  write(*,*) "DEBUG: Size Z real scratch pointer ",size(config%in,3)
  write(*,*) "DEBUG: Size X complex scratch pointer ",size(config%out,1)
  write(*,*) "DEBUG: Size Y complex scratch pointer ",size(config%out,2)
  write(*,*) "DEBUG: Size Z complex scratch pointer ",size(config%out,3)
#endif
       
  ! plan here
  if(mode .eq. FFTW3_ONEPLAN_ESTIMATE) then
    planDepth = FFTW_ESTIMATE
  else if (mode .eq. FFTW3_ONEPLAN_MEASURE) then
    planDepth = FFTW_MEASURE
  else if (mode .eq. FFTW3_ONEPLAN_PATIENT) then
    planDepth = FFTW_PATIENT
  else
    write(*,*) "ERROR: Unknown mode! Stopping.",mode
  endif
  config%fftw3GlobFWD = fftw_mpi_plan_dft_r2c_3d(config%totalDimZ,config%totalDimY,&
        config%totalDimX,config%in,config%out,MPI_COMM,planDepth)
          
  config%fftw3GlobBWD = fftw_mpi_plan_dft_c2r_3d(config%totalDimZ,config%totalDimY,&
        config%totalDimX,config%out,config%in,MPI_COMM,planDepth)
       
#ifdef DEBUG_FFT_INTER
   write(*,*) "DEBUG: Setup stage complete, FWD and BWD plan created!"
   flush(6)
#endif

END SUBROUTINE PlanFFT_NEW

SUBROUTINE ImportWisdom(ioNum)

  IMPLICIT NONE
  INTEGER, INTENT(IN)::ioNum
  ! TODO

END SUBROUTINE ImportWisdom

SUBROUTINE ExportWisdom(ioNum)

  IMPLICIT NONE
  INTEGER, INTENT(IN)::ioNum
  !TODO

END SUBROUTINE ExportWisdom

SUBROUTINE GetFFTDims_NEW(config,dimX,dimY,dimZ,zoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (real-space part)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(OUT) :: &
    dimX, dimY, dimZ, zoff           ! The dimensions of the cell to be FFT'd
    
  TYPE(FFT_CONFIG),INTENT(IN) :: config

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  dimX = config%totalDimX
  dimY = config%totalDimY
  dimZ = config%localDimZ
  zoff = config%localDimZOff

END SUBROUTINE GetFFTDims_NEW


SUBROUTINE GetFFTComplexDims_NEW(config,dimX,dimY,locDimZ,locDimZOffset)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (reciprocal space part)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(OUT) :: &
    dimX, dimY, locDimZ, &            ! The dimensions of the cell to be FFT'd
    locDimZOffset    

  TYPE(FFT_CONFIG),INTENT(IN) :: config

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  dimX = config%totalDimX/2 + 1
  dimY = config%totalDimY
  locDimZ = config%localDimZ
  locDimZOffset = config%localDimZOff

END SUBROUTINE GetFFTComplexDims_NEW

SUBROUTINE ForwardFFT_4D(config,array,transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Real Space --> G space
!   This function is not called directly from the OFDFT code. Use the FFT 
!   interface instead. It performs the transformation of a real 4-dimensional 
!   array into its complex 4-dimensional transform. The first dimension is 
!   halved.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  TYPE(FFT_CONFIG),INTENT(INOUT)::config
  REAL(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to transform

  COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), &
                              SIZE(array,3), SIZE(array,4)) :: &
    transform         ! The answer

                       !>> INTERNAL VARIABLES <<! 
  INTEGER :: &
    is                ! Counter for spin

                         !>> INITIALIZATION <<!
  CALL StartClock('ForwFFT_4D')
                         !>> FUNCTION BODY <<!

  DO is=1, SIZE(array,4)
    CALL ForwardFFT_3D(config,array(:,:,:,is),transform(:,:,:,is))
  END DO !is

  CALL StopClock('ForwFFT_4D')

END SUBROUTINE ForwardFFT_4D

SUBROUTINE BackFFT_4D(config,array,transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of 
!   a complex function over the half-box in reciprocal space back to real 
!   space. It acts on 4-dimensional arrays, the fourth dimension being spin.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!

  TYPE(FFT_CONFIG),INTENT(INOUT)::config
  COMPLEX(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to be back FFT'd

  REAL(kind=DP), DIMENSION(config%totalDimX,config%totalDimY,config%localDimZ, &
                           SIZE(array,4)) :: &
    transform         ! The answer
 
                       !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    is                ! Counter for spin

                        !>> INITIALIZATION <<!   
  CALL StartClock('BackFFT_4D')
                        !>> FUNCTION BODY <<!

  DO is=1, SIZE(array,4)
    CALL BackFFT_3D(config,array(:,:,:,is),transform(:,:,:,is))
  END DO 

  CALL StopClock('BackFFT_4D')

END SUBROUTINE BackFFT_4D

SUBROUTINE ForwardFFT_3D(config,array,transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code. Use the FFT 
!   interface instead. It performs the transformation of a real 3-dimensional 
!   array into its complex 3-dimensional transform. The first dimension is 
!   halved.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!
  TYPE(FFT_CONFIG), INTENT(INOUT) :: config
  REAL(kind=DP), DIMENSION(:,:,:) :: &
    array             ! The array to transform

  COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), &
                              SIZE(array,3)) :: &
    transform         ! The answer
  
                         !>> INITIALIZATION <<! 
#ifdef PROFESS
  CALL StartClock('ForwFFT_3D')
#endif
                         !>> FUNCTION BODY <<!
  config%in(1:config%totalDimX,:,:) = array(1:config%totalDimX,:,:)
  CALL fftw_mpi_execute_dft_r2c(config%fftw3GlobFWD,config%in,config%out)
  
  config%iCountFFT = config%iCountFFT + 1

  ! The forward transform needs to be renormalized afterwards.
#ifdef DEBUG_FFT_INTER
  write(*,*)"DEBUG: unnormalized FWD transformed"
  write(*,*) config%out(:,:,:)
  if(config%out(1,1,1) /= config%out(1,1,1)) then
    write(*,*) "ERROR: NaN case in config%out. Stop."
    flush(6)
    stop
  endif
#endif

#ifdef FFT_NORM_EXACT
  transform(:,:,:) = config%out(:,:,:) / config%normConst
#else
  transform(:,:,:) = config%out(:,:,:) * config%normConst
#endif

#ifdef DEBUG_FFT_INTER
  write(*,*)"DEBUG: normalized FWD transformed"
  write(*,*) transform(:,:,:)
  if(transform(1,1,1) /= transform(1,1,1)) then
    write(*,*) "ERROR: NaN case in transformed. Stop."
    flush(6)
    stop
  endif
#endif

#ifdef PROFESS
  CALL StopClock('ForwFFT_3D')
#endif

END SUBROUTINE ForwardFFT_3D


SUBROUTINE BackFFT_3D(config,array,transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of a 
!   complex function over the half-box in reciprocal space back to real 
!   space. It acts on 3-dimensional arrays.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  TYPE(FFT_CONFIG), INTENT(INOUT) :: config
  COMPLEX(kind=DP), DIMENSION(:,:,:) :: &
    array             ! The array to be back FFT'd

  REAL(kind=DP), DIMENSION(config%totalDimX,config%totalDimY, &
                           config%localDimZ):: &
    transform         ! The answer
 
                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<!
#ifdef PROFESS
  CALL StartClock('BackFFT_3D')
#endif
                         !>> FUNCTION BODY <<!

  config%out(:,:,:) = array(:,:,:)
  CALL fftw_mpi_execute_dft_c2r(config%fftw3GlobBWD,config%out,config%in)
  transform(1:config%totalDimX,:,:) = config%in(1:config%totalDimX,:,:)

  config%iCountFFT = config%iCountFFT + 1

#ifdef PROFESS
  CALL StopClock('BackFFT_3D')
#endif

END SUBROUTINE BackFFT_3D


SUBROUTINE CleanFFT_NEW(config)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine is called at the end of the run to free the memory 
!   associated with the fftw3Plan.
!------------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(FFT_CONFIG), INTENT(IN) :: config
    
  if(config%fftAlgo .eq. FFTW3_ONEPLAN_ESTIMATE  &
     .or. config%fftAlgo .eq. FFTW3_ONEPLAN_MEASURE  &
     .or. config%fftAlgo .eq. FFTW3_ONEPLAN_PATIENT) then
    
     call fftw_destroy_plan(config%fftw3GlobFWD)
     call fftw_destroy_plan(config%fftw3GlobBWD)
  endif
  
  call fftw_free(config%cdata)
  
END SUBROUTINE CleanFFT_NEW

END MODULE Fourier_new

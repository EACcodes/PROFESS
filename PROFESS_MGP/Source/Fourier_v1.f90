MODULE Fourier_NEW
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!    MODULE Fourier_NEW
!       |_SUBROUTINE PlanFFT_NEW
!       |_SUBROUTINE GetFFTDims_NEW
!       |_SUBROUTINE GetFFTComplexDims_NEW
!       |_INTERFACE FFT_NEW
!         |_FUNCTION ForwardFFT_4D (Private)
!         |_FUNCTION BackFFT_4D (Private)
!         |_FUNCTION ForwardFFT_3D (Private)
!         |_FUNCTION BackFFT_3D (Private)
!       |_SUBROUTINE CleanFFT_NEW
!
! DESCRIPTION:
!   This module interfaces with the Fastest Fourier Transform in the World 
!   (FFTW) public library v3.X to provide Fourier transform facilities 
!   for our quantities. Each Fourier transform has to be fftw3Planned for first 
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
!              point errors, start working on wisdom import/export (JMD)
!   10/11/2013 Use fast normalization, fix up for Profess import (JMD)
!   10/15/2013 Start pushing this into some more sane way of dealing with
!              memory...
!------------------------------------------------------------------------------
                             ! << GLOBAL >> !
USE CONSTANTS, ONLY: DP      ! A double precision number

IMPLICIT NONE

INCLUDE 'fftw3.f'

INTEGER,PARAMETER :: &
     FFTW3_ESTIMATE = 0, &
     FFTW3_ONEPLAN_ESTIMATE = 1, &
     FFTW3_ONEPLAN_MEASURE = 2, &
     FFTW3_ONEPLAN_PATIENT = 3


TYPE FFT_CONFIG

   INTEGER :: fftAlgo        ! the mode for the 

   INTEGER :: &
     totalDimX, &            ! dimX
     totalDimY, &            ! dimY, total
     totalDimZ, &            ! dimZ, total
     localDimZ, &         ! the local part of the z dimension
     localDimZOff

   INTEGER :: &
     iCountFFT = 0, &        ! FFT counter for reporting purposes.
     offset                  ! parity of reelRA's X-size.
     
#ifdef FFT_NORM_EXACT
   INTEGER(kind=DP) :: normConst
#else
   REAL(kind=DP) :: normConst
#endif
   
   INTEGER(kind=8)::fftw3GlobFWD,fftw3GlobBWD
   
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

SUBROUTINE PlanFFT_NEW(config,mode,dimX,dimY,dimZ)
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
  
  REAL(kind=DP),DIMENSION(:,:,:),ALLOCATABLE::tmpArrREAL
  COMPLEX(kind=DP),DIMENSION(:,:,:),ALLOCATABLE::tmpArrCMPL
  INTEGER::status

                       !>> INTERNAL VARIABLES <<! 
#ifndef FFT_NORM_EXACT
  INTEGER(KIND=8)::divNorm
#endif
                         !>> INITIALIZATION <<! 
                         
  config%fftAlgo = mode 

  config%totalDimX = dimX
  config%totalDimY = dimY
  config%totalDimZ = dimZ
  config%localDimZ = dimZ
  config%localDimZOff = 0
  
#ifdef FFT_NORM_EXACT
  config%normConst = INT(dimX,kind=DP)*INT(dimY,kind=DP)*INT(dimZ,kind=DP)
#else
  divNorm = INT(dimX,KIND=8)*INT(dimY,KIND=8)*INT(dimZ,KIND=8)
  config%normConst = (REAL(1,KIND=DP))/ (REAL(divNorm,KIND=DP))
#endif
  

  config%offset = MOD(dimX, 2)
  
  if(mode .eq. FFTW3_ONEPLAN_ESTIMATE  &
     .or. mode .eq. FFTW3_ONEPLAN_MEASURE  &
     .or. mode .eq. FFTW3_ONEPLAN_PATIENT) then
     
     ! plan here
     allocate(tmpArrCMPL(dimX,dimY,dimZ),tmpArrREAL(dimX,dimY,dimZ),stat=status)
     if(status .ne. 0) then
       write(*,*) "Failure to allocate temporary arrays for oneshot interface."
       stop
     endif
     
     if(mode .eq. FFTW3_ONEPLAN_ESTIMATE) then
     
       CALL dfftw_plan_dft_r2c_3d(config%fftw3GlobFWD,config%totalDimX,config%totalDimY,&
          config%totalDimZ,tmpArrREAL,tmpArrCMPL,IOR(FFTW_ESTIMATE,FFTW_UNALIGNED))
          
       CALL dfftw_plan_dft_c2r_3d(config%fftw3GlobBWD,config%totalDimX,config%totalDimY,&
          config%totalDimZ,tmpArrCMPL,tmpArrREAL,IOR(FFTW_ESTIMATE,FFTW_UNALIGNED))
          
     else if(mode .eq. FFTW3_ONEPLAN_MEASURE) then
     
       CALL dfftw_plan_dft_r2c_3d(config%fftw3GlobFWD,config%totalDimX,config%totalDimY,&
          config%totalDimZ,tmpArrREAL,tmpArrCMPL,IOR(FFTW_MEASURE,FFTW_UNALIGNED))
          
       CALL dfftw_plan_dft_c2r_3d(config%fftw3GlobBWD,config%totalDimX,config%totalDimY,&
          config%totalDimZ,tmpArrCMPL,tmpArrREAL,IOR(FFTW_MEASURE,FFTW_UNALIGNED))
     
     else if(mode .eq. FFTW3_ONEPLAN_PATIENT) then
     
       CALL dfftw_plan_dft_r2c_3d(config%fftw3GlobFWD,config%totalDimX,config%totalDimY,&
          config%totalDimZ,tmpArrREAL,tmpArrCMPL,IOR(FFTW_PATIENT,FFTW_UNALIGNED))
          
       CALL dfftw_plan_dft_c2r_3d(config%fftw3GlobBWD,config%totalDimX,config%totalDimY,&
          config%totalDimZ,tmpArrCMPL,tmpArrREAL,IOR(FFTW_PATIENT,FFTW_UNALIGNED))
     
     endif
     
#ifdef DEBUG_FFT_INTER
     write(*,*) "DEBUG: FFTW3 forward plan"
     CALL dfftw_print_fftw3Plan(config%fftw3GlobFWD)
     write(*,*) "DEBUG: FFTW3 backward plan"
     CALL dfftw_print_fftw3Plan(config%fftw3GlobBWD)
#endif
     
  endif

END SUBROUTINE PlanFFT_NEW

SUBROUTINE ImportWisdom(ioNum)

  IMPLICIT NONE
  INTEGER, INTENT(IN)::ioNum
  ! TODO
  write(*,*) "INFO: Implementation for ImportWisdom not complete. Assigned I/O unit ",ioNum

END SUBROUTINE ImportWisdom

SUBROUTINE ExportWisdom(ioNum)

  IMPLICIT NONE
  INTEGER, INTENT(IN)::ioNum
  !TODO
  write(*,*) "INFO: Implementation for ExportWisdom not complete. Assigned I/O unit ",ioNum

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

  ! The x-size of the returned array is computed from the reciprocal-space
  ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
  ! k+1 in reciprocal space. We solve the problem by storing the parity.
  REAL(kind=DP), DIMENSION(2*(SIZE(array,1)-1)+config%offset, SIZE(array,2), &
                           SIZE(array,3), SIZE(array,4)) :: &
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
  
  INTEGER(kind=8)::fftw3Plan

                         !>> INITIALIZATION <<! 
  CALL StartClock('ForwFFT_3D')
                         !>> FUNCTION BODY <<!
                         
  if(config%fftAlgo .eq. FFTW3_ONEPLAN_ESTIMATE  &
     .or. config%fftAlgo .eq. FFTW3_ONEPLAN_MEASURE  &
     .or. config%fftAlgo .eq. FFTW3_ONEPLAN_PATIENT) then
    
    CALL dfftw_execute_dft_r2c(config%fftw3GlobFWD,array,transform)
     
  else

    ! FIRST Plan the FFT
    ! A couple of notes:
    !   * Planning a "known" FFT is fast according to the FFTW3 manual
    !   * we are not using wisdom, because... potato...
    !   * I'd like to optionally do FFTW_PATIENT in the future, but this means we 
    !     need to copy the arrays as they are overridden at the fftw3Planning stage
    CALL dfftw_plan_dft_r2c_3d(fftw3Plan,config%totalDimX,config%totalDimY, &
       config%totalDimZ,array,transform,FFTW_ESTIMATE)
       
#ifdef DEBUG_FFT_INTER
    write(*,*) "DEBUG: Local FFTW3 forward plan:"
    call dfftw_print_fftw3Plan(fftw3Plan)
#endif
       
    ! SECOND Execute the fftw3Plan
    CALL dfftw_execute_dft_r2c(fftw3Plan,array,transform)
    ! THIRD Destroy the fftw3Plan
    CALL dfftw_destroy_plan(fftw3Plan)
  
  endif
  
  config%iCountFFT = config%iCountFFT + 1

  ! The forward transform needs to be renormalized afterwards.
#ifdef FFT_NORM_EXACT
  transform = transform / config%normConst
#else
  transform = transform * config%normConst
#endif
  
  CALL StopClock('ForwFFT_3D')

END SUBROUTINE ForwardFFT_3D


SUBROUTINE BackFFT_3D(config,array,transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of a 
!   complex function over the half-box in reciprocal space back to real 
!   space. It acts on 3-dimensional arrays.
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!   08/03/2013  File altered to support FFTW3 (JMD)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  TYPE(FFT_CONFIG), INTENT(INOUT) :: config
  COMPLEX(kind=DP), DIMENSION(:,:,:) :: &
    array             ! The array to be back FFT'd

  ! The x-size of the returned array is computed from the reciprocal-space
  ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
  ! k+1 in reciprocal space. We solve the problem by assuming an odd real size
  REAL(kind=DP), DIMENSION(2*(SIZE(array,1)-1)+config%offset, SIZE(array,2), &
                           SIZE(array,3)) :: &
    transform         ! The answer
    
  INTEGER(kind=8)::fftw3Plan
 
                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<!
  CALL StartClock('BackFFT_3D')
                         !>> FUNCTION BODY <<!

  if(config%fftAlgo .eq. FFTW3_ONEPLAN_ESTIMATE  &
     .or. config%fftAlgo .eq. FFTW3_ONEPLAN_MEASURE  &
     .or. config%fftAlgo .eq. FFTW3_ONEPLAN_PATIENT) then
    
    CALL dfftw_execute_dft_c2r(config%fftw3GlobBWD,array,transform)
     
  else
     
    ! FIRST Plan the FFT
    ! A couple of notes:
    !   * Planning a "known" FFT is fast according to the FFTW3 manual
    !   * we are not using wisdom, because... potato...
    !   * I'd like to optionally do FFTW_PATIENT in the future, but this means we 
    !     need to copy the arrays as they are overridden at the fftw3Planning stage
    CALL dfftw_plan_dft_c2r_3d(fftw3Plan,config%totalDimX,config%totalDimY, &
       config%totalDimZ,array,transform,FFTW_ESTIMATE)
       
#ifdef DEBUG_FFT_INTER
    write(*,*) "DEBUG: Local FFTW3 backward plan:"
    call dfftw_print_fftw3Plan(fftw3Plan)
#endif
       
    ! SECOND Execute the fftw3Plan
    CALL dfftw_execute_dft_c2r(fftw3Plan,array,transform)
    ! THIRD Destroy the fftw3Plan
    CALL dfftw_destroy_plan(fftw3Plan)
  
  endif

  config%iCountFFT = config%iCountFFT + 1

  CALL StopClock('BackFFT_3D')

END SUBROUTINE BackFFT_3D


SUBROUTINE CleanFFT_NEW(config)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine is called at the end of the run to free the memory 
!   associated with the fftw3Plan in case of the guru interface being used for
!   planning..
!------------------------------------------------------------------------------
  IMPLICIT NONE
  TYPE(FFT_CONFIG), INTENT(IN) :: config
    
  if(config%fftAlgo .eq. FFTW3_ONEPLAN_ESTIMATE  &
     .or. config%fftAlgo .eq. FFTW3_ONEPLAN_MEASURE  &
     .or. config%fftAlgo .eq. FFTW3_ONEPLAN_PATIENT) then
    
     call dfftw_destroy_plan(config%fftw3GlobFWD)
     call dfftw_destroy_plan(config%fftw3GlobBWD)
  endif
  
END SUBROUTINE CleanFFT_NEW

END MODULE Fourier_new

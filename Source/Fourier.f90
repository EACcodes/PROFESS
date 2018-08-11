 MODULE Fourier
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!    MODULE Fourier
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
!   This is an ugly shim to make the new, state-free interface work properly
!   with the client code.
!
! CONDITIONS AND ASSUMPTIONS:
!   Nobody will touch for a long time although it should be... 
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Delete this and use the new, state-free interface directly.
!
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/11/2013 File created. (JMD)
!------------------------------------------------------------------------------
                             ! << GLOBAL >> !
USE CONSTANTS, ONLY: DP      ! A double precision number
USE FOURIER_NEW

IMPLICIT NONE

PRIVATE :: &
  BackFFT_4D, &           ! Use FFT
  ForwardFFT_4D, &        ! Use FFT
  BackFFT_3D, &           ! Use FFT
  ForwardFFT_3D           ! Use FFT

INTEGER, SAVE :: offset
INTEGER, SAVE :: iCountFFT

! This interface picks the right transform to perform based on the nature of
! the incomming array: if it's real the FFT is done forward, if complex the 
! back transform is done. All the calls in OFDFT should be of this type: 
! FFT(f).
INTERFACE FFT
  MODULE PROCEDURE ForwardFFT_4D
  MODULE PROCEDURE BackFFT_4D
  MODULE PROCEDURE ForwardFFT_3D
  MODULE PROCEDURE BackFFT_3D
END INTERFACE

CONTAINS

#ifdef __USE_PARALLEL
SUBROUTINE PlanFFT(MPI_COMM,dimX,dimY,dimZ)
#else
SUBROUTINE PlanFFT(dimX,dimY,dimZ)
#endif
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is the initialization procedure that first gets the system name as is
!   called as an argument to OFDFT, and turns it into the various input file
!   names.  Then, it calls all the programs necessary to set variables to 
!   default values, then reads the geometry file to get all the variables sets
!   to the correct values.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(IN) :: dimX, dimY, dimZ ! The dimensions of the cell to be FFT'd
#ifdef __USE_PARALLEL
  INTEGER, INTENT(IN) :: MPI_COMM         ! MPI communicator
#endif

#ifdef __USE_PARALLEL
  CALL PlanFFT_NEW(MPI_COMM,FFT_STD_STATE,FFTW3_ONEPLAN_MEASURE,dimX,dimY,dimZ)
#else
  CALL PlanFFT_NEW(FFT_STD_STATE,FFTW3_ONEPLAN_MEASURE,dimX,dimY,dimZ)
#endif

  offset = MOD(dimX, 2)
  iCountFFT = 0

END SUBROUTINE PlanFFT

SUBROUTINE GetFFTDims(dimX,dimY,locDimZ,locDimZOffset)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (real-space part)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(OUT) :: &
    dimX, dimY           ! The dimensions of the cell to be FFT'd
  INTEGER, INTENT(OUT) :: &
    locDimZ, locDimZOffset  ! Local slice info

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

   CALL GetFFTDims_NEW(FFT_STD_STATE,dimX,dimY,locDimZ,locDimZOffset)


END SUBROUTINE GetFFTDims

SUBROUTINE GetFFTComplexDims(dimX,dimY,locDimZ,locDimZOffset)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (reciprocal space part)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(OUT) :: &
    dimX, dimY            ! The dimensions of the cell to be FFT'd
  INTEGER, INTENT(OUT) :: locDimZOffset, locDimZ

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  CALL GetFFTComplexDims_NEW(FFT_STD_STATE,dimX,dimY,locDimZ,locDimZOffset)

END SUBROUTINE GetFFTComplexDims

! MARKED FOR DEPRECATION!!!
FUNCTION ForwardFFT_4D(array) RESULT(transform)
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
    CALL FFT_NEW(FFT_STD_STATE,array(:,:,:,is),transform(:,:,:,is))
  END DO !is

  CALL StopClock('ForwFFT_4D')

END FUNCTION ForwardFFT_4D

! MARKED FOR DEPRECATION!!!
FUNCTION BackFFT_4D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of 
!   a complex function over the half-box in reciprocal space back to real 
!   space. It acts on 4-dimensional arrays, the fourth dimension being spin.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!

  COMPLEX(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to be back FFT'd

  ! The x-size of the returned array is computed from the reciprocal-space
  ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
  ! k+1 in reciprocal space. We solve the problem by storing the parity.
  REAL(kind=DP), DIMENSION(FFT_STD_STATE%totalDimX,FFT_STD_STATE%totalDimY, &
                           FFT_STD_STATE%localDimZ, SIZE(array,4)) :: &
    transform         ! The answer
 
                       !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    is                ! Counter for spin

                        !>> INITIALIZATION <<!   
  CALL StartClock('BackFFT_4D')
                        !>> FUNCTION BODY <<!

  DO is=1, SIZE(array,4)
    CALL FFT_NEW(FFT_STD_STATE,array(:,:,:,is),transform(:,:,:,is))
  END DO 

  CALL StopClock('BackFFT_4D')

END FUNCTION BackFFT_4D


FUNCTION ForwardFFT_3D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code. Use the FFT 
!   interface instead. It performs the transformation of a real 3-dimensional 
!   array into its complex 3-dimensional transform. The first dimension is 
!   halved.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:) :: &
    array             ! The array to transform

  COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), &
                              SIZE(array,3)) :: &
    transform         ! The answer

  CALL FFT_NEW(FFT_STD_STATE,array,transform)
  iCountFFT = iCountFFT + 1

END FUNCTION ForwardFFT_3D


FUNCTION BackFFT_3D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of a 
!   complex function over the half-box in reciprocal space back to real 
!   space. It acts on 3-dimensional arrays.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!

  COMPLEX(kind=DP), DIMENSION(:,:,:) :: &
    array             ! The array to be back FFT'd

  ! The x-size of the returned array is computed from the reciprocal-space
  ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
  ! k+1 in reciprocal space. We solve the problem by assuming an odd real size
  REAL(kind=DP), DIMENSION(FFT_STD_STATE%totalDimX,FFT_STD_STATE%totalDimY, &
                           FFT_STD_STATE%localDimZ) :: &
    transform         ! The answer
 
                       !>> INTERNAL VARIABLES <<! 
                       !>> INTERNAL VARIABLES <<! 

  CALL FFT_NEW(FFT_STD_STATE,array,transform)
  iCountFFT = iCountFFT + 1

END FUNCTION BackFFT_3D


SUBROUTINE CleanFFT
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine is called at the end of the run to free the memory 
!   associated with the plan.
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!  
                      !>> INTERNAL VARIABLES <<!   
                       !>> INITIALIZATION <<!   
                        !>> FUNCTION BODY <<!

  CALL CleanFFT_NEW(FFT_STD_STATE)

END SUBROUTINE CleanFFT

END MODULE Fourier

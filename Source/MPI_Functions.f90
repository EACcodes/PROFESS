MODULE MPI_Functions 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE MPI_Functions
!     |_SUBROUTINE InitializeMPI 
!     |_SUBROUTINE QuitMPI 
!     |_SUBROUTINE BcastReal
!     |_SUBROUTINE BcastReal_Dim2
!     |_SUBROUTINE BcastCharacter
!     |_INTERFACE ReduceReal
!
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/31/2013  created by Mohan Chen
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP
  USE OutputFiles

  IMPLICIT NONE

#ifdef __USE_PARALLEL
  INCLUDE 'mpif.h'
#endif
  
  INTEGER :: rankGlobal = 0 ! processor rank in total processors
  INTEGER :: sizeGlobal = 1 ! processor size in total 

  !----------------------------
  ! NEW FEATURE IN PROFESS 3.0
  !----------------------------
  INTEGER :: MPI_GFFT_WORLD   ! Communicator name for global FFT.
  INTEGER :: rankGFFT = 0     ! Processor rank for global FFT.
  INTEGER :: sizeGFFT = 1     ! Processor size for global FFT.

  INTEGER :: indGFFT = 0      ! Group index for global FFT.
                              ! There may be more than 1 group,
                              ! but we will only use the first group 
                              ! (indGFFT = 0 ) for global FFT calculations.

  INTEGER :: numProcGFFT = -1 ! The read in processor number for global FFT.
                              ! If this value is -1, the user don't
                              ! specify the number of processors for global
                              ! FFT, then the sizeGFFT = sizeGlobal.
                              ! However, if numProcGFFT > 0, 
                              ! The processor number for global FFT
                              ! will be set to numProcGFFT
  INTEGER :: mpiErr
  INTEGER :: warning = 0      ! close warning


INTERFACE ReduceRealLevel1
  MODULE PROCEDURE ReduceReal_single 
  MODULE PROCEDURE ReduceReal_array
  MODULE PROCEDURE ReduceRealLevel1_2Darray
  MODULE PROCEDURE ReduceRealLevel1_3Darray
END INTERFACE


CONTAINS


SUBROUTINE InitializeMPI()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Get the rank and size for each processor for two communicators.
!   First one is the MPI_COMM_WORLD, for all processors.
!   The second one is MPI_GFFT_WORLD, which is for the processors
!   that will use global FFT.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 06/21/2013 Created (Mohan Chen)
!------------------------------------------------------------------------------
 
  IMPLICIT NONE

#ifdef __USE_PARALLEL
  CALL MPI_INIT(mpiErr) ! Begin the MPI abilities of this program
  IF (mpiErr/=MPI_SUCCESS) STOP "**PROBLEMS INITIALIZING MPI**"

  ! (1) setup communicator for all processors
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rankGlobal, mpiErr)
  CALL MPI_Comm_SIZE(MPI_COMM_WORLD, sizeGlobal, mpiErr);
  IF (rankGlobal==0) WRITE(*,*) "Parallel version, total Processor Number is ", sizeGlobal
  IF (mpiErr/=MPI_SUCCESS) STOP "**PROBLEMS FINDING RANK OF PROCESS**"

  ! (2) new communicator for global FFT group
  ! note that we will only use indGFFT = 0 group to do global FFT.
  IF( numProcGFFT == -1 ) THEN
    indGFFT = 0
    rankGFFT = rankGlobal
    sizeGFFT = sizeGlobal
    MPI_GFFT_WORLD=MPI_COMM_WORLD ! mohan fix bug 2013-07-12
  ELSE IF( numProcGFFT .GT. 0 ) THEN
    indGFFT = rankGlobal/numProcGFFT
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, indGFFT, rankGlobal, MPI_GFFT_WORLD, mpiErr)
    IF(mpiErr /= MPI_SUCCESS) STOP "***Error in calling MPI_COMM_SPLIT***"
    CALL MPI_COMM_RANK(MPI_GFFT_WORLD, rankGFFT, mpiErr)
    IF(mpiErr /= MPI_SUCCESS) STOP "***Error in calling MPI_COMM_RANK***"
    CALL MPI_Comm_SIZE(MPI_GFFT_WORLD, sizeGFFT, mpiErr);
    IF(mpiErr /= MPI_SUCCESS) STOP "***Error in calling MPI_COMM_RANK***"
  ENDIF


  IF( rankGlobal == 0 ) THEN
    WRITE(*,*) "Parallel version of PROFESS"
  ENDIF
#else
  WRITE(*,*) "Serial version of PROFESS"
#endif

  ! outRank is in OutputFiles
  outRank = rankGlobal

END SUBROUTINE InitializeMPI



SUBROUTINE QuitMPI 
!------------------------------------------------------------------------------
! DESCRIPTION:
! finish MPI
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
#ifdef __USE_PARALLEL
  CALL MPI_FINALIZE(mpiErr) ! End the MPI abilities of the program
  IF (mpiErr/=MPI_SUCCESS) STOP "**PROBLEMS FINALIZING MPI**"
#endif

END SUBROUTINE QuitMPI



SUBROUTINE BcastReal(value,valueLen)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Bcast the vaue to all the processors
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
  INTEGER, INTENT(IN) :: valueLen
  REAL(KIND=DP), DIMENSION(valueLen), INTENT(INOUT) :: value

#ifdef __USE_PARALLEL
  CALL MPI_BCAST(value, valueLen, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
#elif DEBUG_PROFESS
  IF( warning .GT. 0 ) THEN
    WRITE(*,*) "ERROR: You are in a MPI function during serial execution. Nothing here for ",valueLen,value
  ENDIF
#endif

END SUBROUTINE BcastReal



SUBROUTINE BcastReal_Dim2(value,valueLen)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Bcast 2D array to all the processors 
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
  INTEGER, INTENT(IN) :: valueLen
  REAL(KIND=DP), DIMENSION(3,valueLen), INTENT(INOUT) :: value

#ifdef __USE_PARALLEL
  CALL MPI_BCAST(value, 3*valueLen, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
#elif DEBUG_PROFESS
  IF( warning .GT. 0 ) THEN
    WRITE(*,*) "ERROR: You are in a MPI function during serial execution. Nothing here for ",valueLen,value
  ENDIF
#endif

END SUBROUTINE BcastReal_Dim2



SUBROUTINE BcastCharacter(word,wordLen)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Bcast the chacracter to all the processors on 
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
  INTEGER, INTENT(IN) :: wordLen
  CHARACTER, DIMENSION(wordLen), INTENT(INOUT) :: word

#ifdef __USE_PARALLEL
  CALL MPI_BCAST(word, wordLen, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) STOP "***PROBLEMS BROADCASTING SYSTEMNAME***"
#elif DEBUG_PROFESS
  IF( warning .GT. 0 ) THEN
    WRITE(*,*) "ERROR: You are in a MPI function during serial execution. Nothing here for ",wordLen,word
  ENDIF
#endif

END SUBROUTINE BcastCharacter


SUBROUTINE ReduceReal_single(value)
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
   
  REAL(KIND=DP), INTENT(INOUT) :: value

  ! >> LOCAL VARIABLES << !
#ifdef __USE_PARALLEL
  REAL(KIND=DP) :: sumValue

  ! summarlize the total points.
  CALL MPI_ALLREDUCE(value, sumValue, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)

  IF (mpiErr /= MPI_SUCCESS) THEN
    WRITE(errorUnit,*) 'Reduce error happens in ReduceRealLevel1()' 
    STOP
  ENDIF
  value = sumValue

#elif DEBUG_PROFESS
  IF( warning .GT. 0 ) THEN
    WRITE(*,*) "ERROR: You are in a MPI function during serial execution. Nothing here for ",value
  ENDIF
#endif

END SUBROUTINE ReduceReal_single
  


SUBROUTINE ReduceReal_array(value, nvalue)
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
  INTEGER, INTENT(IN) :: nvalue
  REAL(KIND=DP), DIMENSION(nvalue), INTENT(INOUT) :: value
#ifdef __USE_PARALLEL
  ! >> LOCAL VARIABLES << !
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: sumValue

  ! summarlize the total points.
  ALLOCATE(sumValue(nvalue))

  CALL MPI_ALLREDUCE(value, sumValue, nvalue, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) THEN
    WRITE(errorUnit,*) 'Reduce error happens in ReduceRealLevel1()' 
    STOP
  ENDIF
  value = sumValue
  DEALLOCATE(sumValue)
#elif DEBUG_PROFESS
  IF( warning .GT. 0 ) THEN
    WRITE(*,*) "ERROR: You are in a MPI function during serial execution. Nothing here for ",nvalue,value
  ENDIF
#endif
END SUBROUTINE ReduceReal_array



SUBROUTINE ReduceRealLevel1_2Darray(value, nx, ny)
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
  INTEGER, INTENT(IN) :: nx, ny
  REAL(DP), DIMENSION(nx,ny), INTENT(INOUT) :: value
  !! >> LOCAL VARIABLES << !!
#ifdef __USE_PARALLEL
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sumValue
  INTEGER :: nxy
  !! >> INITIALIZE << !!
  ! summarlize the total points.
  nxy = nx*ny
  ALLOCATE(sumValue(nx,ny))
  CALL MPI_ALLREDUCE(value, sumValue, nxy, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) THEN
    WRITE(errorUnit,*) 'Reduce error happens in ReduceRealLevel1()' 
    STOP
  ENDIF
  value = sumValue
  DEALLOCATE(sumValue)
#elif DEBUG_PROFESS
  IF( warning .GT. 0 ) THEN
    WRITE(*,*) "ERROR: You are in a MPI function during serial execution. Nothing here for ",value,nx,ny
  ENDIF
#endif
END SUBROUTINE ReduceRealLevel1_2Darray
  


SUBROUTINE ReduceRealLevel1_3Darray(value, nx, ny, nz)
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
  INTEGER, INTENT(IN) :: nx, ny, nz
  REAL(DP), DIMENSION(nx,ny,nz), INTENT(INOUT) :: value
  !! >> LOCAL VARIABLES << !!
#ifdef __USE_PARALLEL
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: sumValue
  INTEGER :: nxyz
  !! >> INITIALIZE << !!
  ! summarlize the total points.
  nxyz = nx*ny*nz
  ALLOCATE(sumValue(nx,ny,nz))
  CALL MPI_ALLREDUCE(value, sumValue, nxyz, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) THEN
    WRITE(errorUnit,*) 'Reduce error happens in ReduceRealLevel1()' 
    STOP
  ENDIF
  value = sumValue
  DEALLOCATE(sumValue)
#elif DEBUG_PROFESS
  IF( warning .GT. 0 ) THEN
    WRITE(*,*) "ERROR: You are in a MPI function during serial execution. Nothing here for ",value,nx,ny,nz
  ENDIF
#endif
END SUBROUTINE ReduceRealLevel1_3Darray


SUBROUTINE ReduceInteger_single(value)
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

  INTEGER, INTENT(INOUT) :: value

  ! >> LOCAL VARIABLES << !
#ifdef __USE_PARALLEL
  INTEGER :: sumValue

  ! summarlize the total points.
  CALL MPI_ALLREDUCE(value, sumValue, 1, MPI_INTEGER, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) THEN
    WRITE(errorUnit,*) 'Reduce error happens in ReduceRealLevel1()'
    STOP
  ENDIF
  value = sumValue
#elif DEBUG_PROFESS
  IF( warning .GT. 0 ) THEN
    WRITE(*,*) "ERROR: You are in a MPI function during serial execution. Nothing here for ",value
  ENDIF
#endif

END SUBROUTINE ReduceInteger_single



END MODULE MPI_Functions

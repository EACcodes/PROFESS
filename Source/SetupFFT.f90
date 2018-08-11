MODULE SetupFFT
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE SetupFFT 
!     |_SUBROUTINE InitializeFFT
!     |_SUBROUTINE FFTcounterReport
!     |_SUBROUTINE SizeSystem
!       |_FUNCTION PadFFTdims
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
!   02/25/2004 Pieced together for organizational purposes. (GSH)
!   03/10/2004 Deleted numEle (the number of electrons) (GSH)
!   08/22/2013 MOHAN CHEN update
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY : DP ! Double Precision
  USE MPI_Functions

  IMPLICIT NONE

  INTEGER :: dimType = 1 
  ! Type of dimensions we want in the cell.
  !   1 = the lowest multiples of 2, 3, 5, 7
  !   2 = the lowest multiples of 2, 3, 5, 7 but at 
  !       least one occurance of 2 (even number)
  !   3 = the lowest multiples of 3, 5, 7 (odd number)
  !   4 = strictly powers of 2 if periodic
  !       powers of 2 + 1 if dirichlet
  !   5 = strictly powers of 2 plus one instance of 
  !       3, 5, or 7
  !

CONTAINS


SUBROUTINE InitializeFFT(energyCutoff, gridSpacing, totX, totY, totZ, locZ, locZOff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine takes the cell from ReadGeometry to compute the table of
!   q-vectors. It then calls the subroutine FillQTable, that is separate
!   because it needs to be called every time we alter the cell shape. In event
!   that the WT or WGC kinetic energy functionals are used, it will also
!   allocate and fill the kernel table. Once again, the allocation and filling
!   parts are kept separate to allow updating of the kernel when the cell size
!   is altered.
!
! CONDITIONS AND ASSUMPTIONS:
!   It is assumed that the real-space grid is odd along each of its
!   dimensions. If you want to change this, you need to read the warning note
!   in the initialization part or expose yourself to considerable trouble.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   1. Ashcroft and Mermin, Solid State Physics.  Harcourt College Publishers,
!      Fort Worth, 1976.
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 08/22/2013 MOHAN CHEN
! 10/23/2013 Update to FFTW3 interface (JMD)
!------------------------------------------------------------------------------
#ifdef __USE_PARALLEL
  USE FOURIER, ONLY: GetFFTComplexDims, GetFFTDims
#endif
  USE FOURIER_NEW
  USE OutputFiles, ONLY : outputUnit ! Add by Mohan
  USE CellInfo, ONLY : cell

  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  REAL(KIND=DP), INTENT(IN) :: gridSpacing
  REAL(kind=DP), INTENT(IN) :: energyCutoff        ! The kinetic energy cutoff
  INTEGER, INTENT(INOUT):: totX, totY, totZ, locZ, locZOff

  ! >> LOCAL VARIABLES << !!

  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           SETUP FFT"

  IF(energyCutoff > 0._DP) THEN
    ! we determine the FFT size according to 'sizeGFFT',
    ! 'sizeGFFT' is the number of processors that used for
    ! global FFT calculations, it doesn't need to be the 
    ! total numbe of processors that you use. In fact,
    ! it can be smaller than it. This is useful when
    ! we want to scale the calculations up to thousands of
    ! processors.
    CALL SizeSystem(energyCutoff, sizeGFFT, totX, totY, totZ)
#ifdef DEBUG_FFT_INTER
    WRITE(*,*) "DEBUG: SizeSystem computed FFT dimensions ",totX,totY,totZ
#endif
  ELSE
    totX = NINT(gridSpacing*cell%lengthA)
    totY = NINT(gridSpacing*cell%lengthB)
    totZ = NINT(gridSpacing*cell%lengthC)

#ifdef DEBUG_FFT_INTER
    WRITE(*,*) "DEBUG: Computed FFT dimensions ",totX,totY,totZ
    WRITE(*,*) "DEBUG: gridspaceing ",gridSpacing
    WRITE(*,*) "DEBUG: cell dimensions: ",cell%lengthA,cell%lengthB,cell%lengthC
#endif

    WRITE(*,*) "WARNING: GRIDPOINT DENSITY HAVE BEEN SPECIFIED MANUALLY."
    WRITE(*,*) "  THIS MAY BE NON-OPTIMAL FOR THE FFT SPEED. BE CAREFUL UNTIL"
    WRITE(*,*) "  THIS FEATURE IS IDIOT-PROOFED."

  END IF


  WRITE(outputUnit, '(A,3I5)') " Global FFT Dimension : ", totX, totY, totZ
  IF(rankGlobal==0) THEN
    WRITE(6, '(A,3I5)') " (Setup System) Global FFT Dims.         : ", totX, totY, totZ
  ENDIF


   ! Set up FFTs for interior points
   ! Periodic boundary condition
#ifdef __USE_PARALLEL
    ! MEASURE WILL INTRODUCE RANDOM SMALL ERRORS AT FIRST. (mohan note)
    !CALL PlanFFT_NEW(MPI_GFFT_WORLD,FFT_STD_STATE,FFTW3_ONEPLAN_MEASURE,totX, totY, totZ)
    CALL PlanFFT_NEW(MPI_GFFT_WORLD,FFT_STD_STATE,FFTW3_ONEPLAN_ESTIMATE,totX, totY, totZ)
#else
    ! MEASURE WILL INTRODUCE RANDOM SMALL ERRORS AT FIRST. (mohan note)
    !CALL PlanFFT_NEW(FFT_STD_STATE,FFTW3_ONEPLAN_MEASURE,totX, totY, totZ)
    CALL PlanFFT_NEW(FFT_STD_STATE,FFTW3_ONEPLAN_ESTIMATE,totX, totY, totZ)
#endif
 
#ifdef __USE_PARALLEL

  WRITE(outputUnit,*) "Rank/Size in Global Communicator     : ", rankGlobal, "/", sizeGlobal
  WRITE(outputUnit,*) "Rank/Size in Global FFT Communicator : ", rankGFFT, "/", sizeGFFT
!  WRITE(outputUnit,*) "Processor number for Global FFT  : ", numProcGFFT
!  WRITE(outputUnit,*) "Index in Global FFT Communicator : ", indGFFT

  ! For parallel OFDFT, find what information is on local processors
  CALL GetFFTComplexDims(totX,totY,locZ,locZOff)
  CALL GetFFTDims(totX,totY,locZ,locZOff)
  
  WRITE(outputUnit,"(A,I5)") " Parallel FFT, totX    = ",totX
  WRITE(outputUnit,"(A,I5)") " Parallel FFT, totY    = ",totY
  WRITE(outputUnit,"(A,I5)") " Parallel FFT, locZ    = ",locZ
  WRITE(outputUnit,"(A,I5)") " Parallel FFT, locZoff = ",locZoff

#else
  locZ = totZ
  locZOff = 0 ! mohan add
#endif

#ifdef __USE_PARALLEL

  ! Double check that locZOff ranks are the same as MPI_COMM_WORLD ranks.
  ! This should be true for FFTW, according to FFTW people, but isn't
  ! explicitly stated in the documentation. This check MUST BE TRUE for
  ! interprocessor communications to proceed correctly
  IF (locZOff/=rankGFFT*locZ) STOP &
    "**FFTW Z-DIM RANKINGS DIFFERENT FROM EXPECTED: PARALLEL WON'T WORK**"

  ! Reset locZOff if using Dirichlet boundary conditions to allow optimal
  ! communications between padded and non-padded arrays in real space.
  !
  ! Example: 5 processors (P0-P4), with 10 gridpoints in z-direction for FFTs,
  !          5 for typical unpadded arrays.  Data is N0, N1, ..., N9, with
  !          N5-N9 the padding.  For FFTs, P0 has N0,N1; P1 has N2,N3; etc.
  !          When unpadded, P0 has N0; P1 has N2; P2 has N4; P3 has N3;
  !          P4 has N1.
  !
  ! Variables totX, totY, totZ, locZ is reset for non-padded arrays, with each
  ! processor keeping only half the number of gridpoints needed for FFTs.
  ! The processors with ranks 0 to (sizeGFFT-1)/2 maintain the same offsets,
  ! but offsets for the remainder of the processors are reset.
  !
  
#else

#endif

END SUBROUTINE InitializeFFT



SUBROUTINE SizeSystem(energyCutoff, numProc, numX, numY, numZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine takes the energy cutoff assigned in ReadOptions and the cell
!   from ReadGeometry to compute the size of the real-space (and equivalently,
!   the reciprocal-space) grid.
!
! CONDITIONS AND ASSUMPTIONS:
!   It is assumed that the real-space grid is odd along each of its
!   dimensions.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   08/24/2006 Separated from SetupSystem so that PlanFFT will know the total
!              size of rhoR and can allocate local rhoR arrays (Linda Hung)
!
!------------------------------------------------------------------------------
  USE CellInfo, ONLY : cell
  USE CONSTANTS, ONLY : PI

  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN) :: energyCutoff  ! The kinetic energy cutoff
  INTEGER, INTENT(IN) :: numProc             ! Number of processors used in this run
  INTEGER, INTENT(OUT) :: numX, numY, numZ   ! Number of points on real space (and also
                                             ! reciprocal space, as it turns out) axes.

                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: qMax                ! Maximum value of |q| based on energyCutoff
  REAL(KIND=DP) :: m1Max, m2Max, m3Max ! Maximum coefficients for b1, b2, b3.

                           !>> INITIALIZATION <<!

                           !>> FUNCTION BODY <<!

  ! Sets the cutoff value of q.  E = hbar**2 * qMax**2 / 2me, but hbar and me
  ! are 1. (Eq. 2.7 in Ref 1)
  qMax = SQRT(energyCutoff * 2._DP)

  ! Since q = 2*pi*m / L, mMax (the number of times a unit 2*pi / L goes
  ! into qMax) can be found by : mMax = qMax / (2*pi/L).  However, we multiply
  ! qMax here by 2 and we shift the index m so that it starts from 0 (for the
  ! purposes of the FFT).  Originally, -mMin = mMax.  So then we really
  ! get mMax = qMax * 2 / (2 * pi / L) ==> mMax = qMax * L / pi (Eq. 2.16)

  m1Max = qMax * cell%lengthA / PI
  m2Max = qMax * cell%lengthB / PI
  m3Max = qMax * cell%lengthC / PI

  numX = PadFFTdims(m1Max,  .FALSE.)
  !!!GREG: I temporarily made this "TRUE" because otherwise some of my
  !!!      parallel jobs will have different x and y grid #'s when
  !!!      the cell sizes are the same.  This is temporarily inconvenient
  !!!      for me.
!  numX = PadFFTdims(m1Max, bound(1), .TRUE.)
  numY = PadFFTdims(m2Max, .TRUE.)
  numZ = PadFFTdims(m3Max, .TRUE.)
  ! Changes the real number miMax into integer that is a multiple of 3, 5, and
  ! 7 only, for FFT.  Numbers odd by construction, so no info lost.
  ! If TRUE on third input, Pad is a multiple of the number of processors.
  ! WARNING on changing these formulas: Real-space grid dimensions are
  ! assumed to be odd. If you tweak these formulas to get even grid sizes in
  ! either x or z, you are asking for trouble. For z even, the x=0 plane
  ! does not contain all the information necessary to rebuild the full sphere.
  ! FFTW seems to be getting around this problem somehow, possibly by storing
  ! the necessary information in a space that would otherwise be occupied by
  ! redundant data. It may or may not cause trouble, but try it at your own
  ! risk. In x, the back FFT routines (3D AND 4D) reconstruct the real-space
  ! dimension from the reciprocal space length. The ambiguity is lifted by
  ! assuming the real-space length is odd. You need to change BackFFT_4D and
  ! BackFFT_3D in the FOURIER module if you want to have numX even or there
  ! WILL be problems. (VLL)

  CONTAINS

  FUNCTION PadFFTdims(R, procMultiple)
  !----------------------------------------------------------------------------
  ! DESCRIPTION:
  ! This function converts a real, positive number R into an integer of equal
  ! or greater value that is a multiple of 3, 5 and 7 only.  It also may
  ! require that the number be a multiple of the number of processors.
  !
  ! REVISION LOG:
  !   12/25/2003 This file was created by VLL.  I'm just moving it here. (GSH)
  !   09/14/2006 Pad can be a multiple of the number of processors, if needed
  !              for better FFTW parallelization (LH)
  !
  !----------------------------------------------------------------------------
    IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

    INTEGER :: PadFFTdims                  ! The integer output by the function.
    REAL(KIND=DP), INTENT(IN) :: R         ! The real number in input.
    LOGICAL, INTENT(IN) :: procMultiple    ! If TRUE, Pad is a multiple of numProc

                       !>> INTERNAL VARIABLES <<!

    INTEGER :: P2, P3, P5, P7     ! Minimum powers of 2, 3, 5, and 7 that
                                  ! top R/REAL(numProc).
    INTEGER :: minimum            ! Current smallest combination that top R.
    INTEGER :: start_i2, end_i2   ! The start and end value of i2
    INTEGER :: i2, i3, i5, i7     ! Dummy counters
    INTEGER :: try                ! Tryout value
    REAL(KIND=DP) :: logR         ! Natural logarithm of R
    REAL(KIND=DP) :: procR        ! Minimum grid size per processor

    REAL(kind=DP), PARAMETER :: eps = 1.E-14_DP ! small number for integer case

                       !>> INITIALIZATION <<!

    ! If we want Pad to be a multiple of the number of processors, we proceed
    ! the same as usual, except with R replaced with R/numProc to find out the
    ! appropriate powers of 3,5,7, etc.  The value of numProc is multiplied
    ! back in the end.
    IF (procMultiple) THEN
      procR = R/REAL(numProc,kind=DP)
    ELSE
      procR = R
    END IF

    logR = LOG(R)

                       !>> FUNCTION BODY <<!
    IF(dimType == 1 .OR. dimType == 2 .OR. dimType == 3) THEN

      P2 = INT(logR/LOG(2.0_DP)-eps+1)
      P3 = INT(logR/LOG(3.0_DP)-eps+1) ! This formula gives correct PN value
      P5 = INT(logR/LOG(5.0_DP)-eps+1) ! even if R/numProc=N^PN, thanks to eps.
      P7 = INT(logR/LOG(7.0_DP)-eps+1)

      minimum = 7 ** P7                ! Sets the starting point.

      IF(dimType == 1) THEN
        start_i2 = 0
        end_i2 = P2
      ELSE IF(dimType == 2) THEN
        start_i2 = 1
        end_i2 = P2
      ELSE IF(dimType == 3) THEN
        start_i2 = 0
        end_i2 = 0
      END IF

      ! This is less efficient than above, but it's fast anyway.
      DO i7 = 0, P7
        DO i5 = 0, P5
          DO i3 = 0, P3
            DO i2 = start_i2, end_i2
              try = 7**i7 * 5**i5 * 3**i3 * 2**i2
              IF(try>=procR) THEN
                IF (try < minimum) minimum = try
                EXIT
              END IF

            END DO
          END DO
        END DO
      END DO

      PadFFTdims = minimum

    ELSE IF(dimType == 4) THEN
      PadFFTdims = 2 ** (INT(LOG(procR)/LOG(2._DP) - eps + 1))

    ELSE IF(dimType == 5) THEN
      PadFFTdims = MIN(2 ** (INT(LOG(procR/1._DP)/LOG(2._DP) - eps + 1)) * 1, &
                2 ** (INT(LOG(procR/3._DP)/LOG(2._DP) - eps + 1)) * 3, &
                2 ** (INT(LOG(procR/5._DP)/LOG(2._DP) - eps + 1)) * 5, &
                2 ** (INT(LOG(procR/7._DP)/LOG(2._DP) - eps + 1)) * 7)

    ELSE
      STOP "**Invalid DIMENSION option**"

    END IF

    IF (procMultiple) THEN
      PadFFTdims = PadFFTdims*numProc
    ENDIF

  END FUNCTION PadFFTdims

END SUBROUTINE SizeSystem

END MODULE SetupFFT

MODULE Output
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Output
!     |_SUBROUTINE PrintForces
!     |_SUBROUTINE PrintStress
!     |_SUBROUTINE PrintLattice
!     |_SUBROUTINE PrintIons
!     |_SUBROUTINE PrintPseudo
!     |_SUBROUTINE PrintGeometry
!     |_SUBROUTINE PrintDensity
!     |_SUBROUTINE PrintPotential
!     |_INTERFACE Dump
!       |_SUBROUTINE DumpReal4D
!       |_SUBROUTINE DumpReal3D
!       |_SUBROUTINE DumpComplex
!       |_SUBROUTINE DumpLogical
!     |_SUBROUTINE PrintEwald
!     |_SUBROUTINE PrintTransitionState
!     |_SUBROUTINE CleanOutput
!     |_SUBROUTINE WrtOut
!     |_SUBROUTINE QUIT
!     |_SUBROUTINE PrintDensityClearly
!     |_SUBROUTINE TITLE
!
! DESCRIPTION:
!   This is a module that contains all the pretty printing routines that print
!   to any and all of the OFDFT output files.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/15/03  File created (GSH)
!
!------------------------------------------------------------------------------
                              !<< GLOBAL >>!

  USE CONSTANTS, ONLY : &
  DP, &                    ! Shortcut for "double precision"
  hartreeToeV, &           ! conversion from hartrees to eV
  bohr, &                  ! an atomic unit of length
  boltzmann, &             ! Boltzmann's constant
  auToGPa                  ! Conversion to gigapascals

  USE SYS, ONLY : &
  grids, &
  forceIon, &       ! The forces on the ions
  rhoR, &           ! Electron density in real space
  potReal, &        ! Potential in real space
  stress, &         ! The stress on the cell
  frozenIon, &      ! whether to optimize ion position
  energy            ! The energy of the system

  USE CellInfo, ONLY :  cell, m3G, n3G, n3Goff
  ! cell: Contains data about the cell, lattice positions ion positions

  USE MATHFUNCTIONS, ONLY : Volume
  ! Gets the volume of the cell

  USE TIMER, ONLY : &
  TimerStart, &
  TimerStop, &
  stopwatch

  USE OUTPUTFILES, ONLY : &
  outputUnit, &
  errorUnit

  IMPLICIT NONE

  ! The following are the parameters for all the output files.
  ! Note:  We've used 5 for the input files, so don't use 5.
  INTEGER, PARAMETER :: &
  densityUnit = 11, &      ! density files
  potentialUnit = 12, &    ! potential files
  outputGeomUnit = 13, &   ! geometry/stress files
  forceUnit = 14, &        ! force files
  transUnit = 15, &        ! transition state output file
  junkUnit = 16, &         ! junk/dumping output for debugging
  lineLength = 78          ! As a matter of convention, we have 
                           ! 78 characters in a line.
                                                          
  ! What should be printed out
  LOGICAL :: &
  outputOptionGeom = .FALSE., &
  outputFinalForces = .FALSE., &   ! Print out forces, controlled by user 
  outputEwald = .FALSE., &            ! Print out information regarding the 
                                      ! Ewaldsummation
  outputFinalPotential = .FALSE. , &  ! Print out the final electronic 
                                      ! potential.
  outputIntermediateDensity = .FALSE., & ! Prints the intermedate density
  outputFinalDensity = .FALSE., &     ! Print out the FINAL electron density
  outputKernel = .FALSE., &           ! Print out the kernels for KEDF.
  outputFinalGeometry = .FALSE., &    ! Print out the FINAL electron density
  outputFinalStress = .FALSE., &      ! Print out the FINAL stress
  outputForcesSeparate = .FALSE., &   ! Separate forces into another file
  outputTransitionState = .FALSE. 

  REAL(KIND=DP) :: maxMB=2000
  ! Maximum size of density output file (default 2GB)

  INTEGER, DIMENSION(:), POINTER :: outputIonTable
  ! The ion table, unsorted, for output
  !

  INTEGER :: &
  geoOutputFreq = -1,       &    ! Ion coords output every geoOutputFreq steps
  celOutputFreq = -1,       &    ! Ion coords output every celOutputFreq steps
  exchangeCorrelationFunct, &    ! Which exchange correlation functional
                                 ! is being used
  outputRhoMethod, &             ! the METHOD of density minimization
  outputIonMethod, &             ! the METHOD of ion minimization
  outputMinimizeGeometry = 0, &  ! The level of detail for minimizing  
                                 ! the geometry.  Default is low.
  outputMinimizeDensity = 0, &   ! The level of detail for minimizing the 
                                 ! density.  Default is low (0)
  outputRank = 0                 ! Rank of process

  CHARACTER(LEN=30) :: &
  outputSystemName               ! The raw argument on the command line 

  REAL(KIND=DP), DIMENSION(SIZE(energy)) :: &
  energyTime = 0.0_DP, &         ! the time to calculate a single energy
  singlePointTime = 0.0_DP, &    ! total time spent in this calculation
                                 ! for this ionic structure
  totalEnergyTime = 0.0_DP, &       ! total time spent for each energy in this
  totalPotentialTime = 0.0_DP       ! total time spent for each energy in this
                                 ! calculation 

  LOGICAL :: &
  jumpToMax, &                   ! We jumped to the max value
  restarted                   ! Did we restart the minimization?

  TYPE(stopwatch) :: &
  rhoWatch                       ! Timer for electronic minimization

! This dumps out an array into a readable form for pasting directly into a
! spreadsheet, such as excel, for debugging
INTERFACE Dump
  MODULE PROCEDURE DumpReal3D
  MODULE PROCEDURE DumpReal4D
  MODULE PROCEDURE DumpComplex
  MODULE PROCEDURE DumpLogical
END INTERFACE


CONTAINS


SUBROUTINE PrintForces(forces, systemName, ionStep)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out the forces into a file.  The file is named 
!   by the argument given, "fileName" followed by a number tag assigned 
!   sequentially every time the subroutine is run.
!  
! CONDITIONS AND ASSUMPTIONS: Called if outputFinalForces=.TRUE. or
!   outputMinimizeGeometry >=3
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/20/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:) :: forces 
  ! Total forces
  !
  CHARACTER(LEN=*) :: systemName         
  ! The raw argument on the command line 
  !
  INTEGER, OPTIONAL, INTENT(in) :: ionStep  
  ! the current step # of ion relax
  !

                       !>> INTERNAL VARIABLES <<!

  ! Because this is initialized in the declaration, the value will be 
  ! preserved from one execution of the program to another. 
  INTEGER :: &
    i, &
    fileStatus, &    ! Status integer to see if the file opening succeeded
    fileUnit

  CHARACTER(len=50) :: &
    tmp, &
    forceFile

  CHARACTER :: &
    border

                       !>> INITIALIZATION <<!

  IF (PRESENT(ionStep)) THEN 
    WRITE(tmp,*) ionStep
    forceFile = TRIM(systemName) // "." // TRIM(ADJUSTL(tmp)) // ".force"
  ELSE
    forceFile = TRIM(systemName) // ".force.out"
  ENDIF


  IF(outputForcesSeparate) THEN
    fileUnit = forceUnit  
    border = " "
    OPEN(unit=forceUnit, access="sequential", action="write", blank="null", &
         file = trim(forceFile), form="formatted", &
         iostat=fileStatus, position="rewind", status="replace")
    IF (fileStatus/=0) THEN
      WRITE(*,*)'The OUTPUT force file ', trim(forceFile), &
      ' had a problem opening.'
      STOP
    END IF  
  ELSE
    fileUnit = outputUnit
    border = ":"
  END IF

  ! Format Descriptors
  11 FORMAT(A1, 1X, A2, 1X, I5, 1X, 3(ES18.10e2, 1X), 9X,  A1)
  12 FORMAT(A1, A76, A1)

                       !>> FUNCTION BODY <<!

  ! Output the forces, remember to convert
  WRITE(fileUnit, 12) " "
  WRITE(fileUnit, 12) border, REPEAT("-", 76), border
  WRITE(fileUnit, '(A)') " TOTAL-FORCE (eV/A)"
  WRITE(fileUnit, 12) border, REPEAT("-", 76), border
  WRITE(fileUnit, 12) border, &
                      REPEAT(" ", 6) // "x" // REPEAT(" ", 18) // "y" // &
                      REPEAT(" ", 18)// "z" // REPEAT(" ", 22), &
                      border

  ! print the force for each atom
  DO i = 1, SIZE(cell%ionTable)
    WRITE(fileUnit, 11) &
      border, &
      cell%elementTable(cell%ionTable(outputIonTable(i))%elementID)%elementName, &
      i, &
      forces(outputIonTable(i),:)*hartreeToeV/bohr, &
      border
  END DO
  WRITE(fileUnit, 12) border, REPEAT(".", 76), border
  CLOSE(forceUnit)


  ! in a.u. Hartree/Bohr
  ! mohan add 2014-09-18
  WRITE(fileUnit, 12) " "
  WRITE(fileUnit, 12) border, REPEAT("-", 76), border
  WRITE(fileUnit, '(A)') " TOTAL-FORCE (Hartree/Bohr)"
  WRITE(fileUnit, 12) border, REPEAT("-", 76), border
  WRITE(fileUnit, 12) border, &
                      REPEAT(" ", 6) // "x" // REPEAT(" ", 18) // "y" // &
                      REPEAT(" ", 18)// "z" // REPEAT(" ", 22), &
                      border

  ! print the force for each atom
  DO i = 1, SIZE(cell%ionTable)
    WRITE(fileUnit, 11) &
      border, &
      cell%elementTable(cell%ionTable(outputIonTable(i))%elementID)%elementName, &
      i, &
      forces(outputIonTable(i),:), &
      border
  END DO
  WRITE(fileUnit, 12) border, REPEAT(".", 76), border
  CLOSE(forceUnit)


END SUBROUTINE PrintForces


SUBROUTINE PrintStress(stress)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints the stress tensor, along with the lattice vectors.
!
! CONDITIONS AND ASSUMPTIONS: Only call when outputFinalStress=.TRUE.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:) :: stress ! the stress tensor
                       !>> INTERNAL VARIABLES <<!
  INTEGER :: i
  CHARACTER :: border
  CHARACTER(len=30) :: comp

                       !>> INITIALIZATION <<!

  border = ";"

  ! Format Descriptors
  11 FORMAT(A1, 1X, A1, 1X, 3(ES18.10e2, 1X), 16X,  A1)
  12 FORMAT(A1, A76, A1)

                       !>> FUNCTION BODY <<!

  ! Output the lattice vectors, convert to A.
  WRITE(outputUnit, *) " "
  WRITE(outputUnit, *) "------------------------------------------------------------------------------"
  WRITE(outputUnit, *) "                        LATTICE VECTORS (Angstrom)"
  DO i = 1, 3
    IF(i == 1) comp = "LAT 1"
    IF(i == 2) comp = "LAT 2"
    IF(i == 3) comp = "LAT 3"
    WRITE(outputUnit, '(1X, A8, 1X, 3(F15.6, 1X))') comp(1:8),cell%cellReal(:,i) * bohr
  END DO

  ! Output the stresses, remember to convert
  WRITE(outputUnit, *) "------------------------------------------------------------------------------"
  WRITE(outputUnit, *) "                               STRESS (GPa)"
  DO i = 1, 3
    IF(i == 1) comp = "ST_GPa 1"
    IF(i == 2) comp = "ST_GPa 2"
    IF(i == 3) comp = "ST_GPa 3"
    WRITE(outputUnit, '(1X, A8, 1X, 3(F15.6, 1X))') comp(1:8), & 
       stress(i,:) * auToGPa 
  END DO
  
! mohan comment this out on 2014-05-02
! mohan comment back 2014-09-18, in order to be consistent with manual
  WRITE(outputUnit, *) " "
  WRITE(outputUnit, *) "--------------------------------------------------"
  WRITE(outputUnit, *) "STRESS (a.u.)"
  WRITE(outputUnit, *) "--------------------------------------------------"
  DO i = 1, 3
    IF(i == 1) comp = "ST_AU  1"
    IF(i == 2) comp = "ST_AU  2"
    IF(i == 3) comp = "ST_AU  3"
    WRITE(outputUnit, '(1X, A8, 1X, 3(ES18.10, 1X))') comp(1:8), stress(i,:)
  END DO

END SUBROUTINE PrintStress


SUBROUTINE PrintLattice()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function prints out the new lattice in the output geometry file
!
! CONDITIONS AND ASSUMPTIONS: Only call when outputFinalStress=.TRUE.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/20/2003  Created (Greg Ho)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

                       !>> INTERNAL VARIABLES <<!

  ! Because this is initialized in the declaration, the value will be 
  ! preserved from one execution of the program to another. 
  INTEGER :: &
    i                ! counter

                       !>> INITIALIZATION <<!

  ! Format Descriptors
  10 FORMAT(A)
  11 FORMAT(3F20.15)

                       !>> FUNCTION BODY <<!
  WRITE(outputGeomUnit, 10) "%BLOCK LATTICE_CART"
  DO i=1, 3
    WRITE(outputGeomUnit, 11) cell%cellReal(:,i) * bohr
  END DO
  WRITE(outputGeomUnit, 10) "%END BLOCK LATTICE_CART"

END SUBROUTINE PrintLattice


SUBROUTINE PrintIons(positions)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function prints out the new ion positions in the output geometry file 
!
! CONDITIONS AND ASSUMPTIONS: Printed only when outputOptionGeom=.TRUE.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/20/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:), OPTIONAL :: positions

                       !>> INTERNAL VARIABLES <<!

  ! Because this is initialized in the declaration, the value will be 
  ! preserved from one execution of the program to another. 
  INTEGER :: i ! counter

                       !>> INITIALIZATION <<!

  ! Format Descriptors
  10 FORMAT(A)
  11 FORMAT(A6, 3F20.15)

                       !>> FUNCTION BODY <<!

  WRITE(outputGeomUnit, 10) "%BLOCK POSITIONS_FRAC"
  DO i=1, SIZE(cell%ionTable)
    IF(PRESENT(positions)) THEN
      WRITE(outputGeomUnit, 11) &
        cell%elementTable(cell%ionTable(outputIonTable(i))%elementID)%elementName, &
        positions(i,:)
    ELSE

      WRITE(outputGeomUnit, 11) &
        cell%elementTable(cell%ionTable(outputIonTable(i))%elementID)%elementName, &
        cell%ionTable(outputIonTable(i))%coord
    END IF
  END DO
  WRITE(outputGeomUnit, 10) "%END BLOCK POSITIONS_FRAC"

END SUBROUTINE PrintIons


SUBROUTINE PrintPseudo
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function prints out the pseudopotential name in the output geometry
!   file.
!
! CONDITIONS AND ASSUMPTIONS: Printed only when outputOptionGeom=.TRUE.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/20/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!
                       !>> INTERNAL VARIABLES <<!
  INTEGER :: i

                       !>> INITIALIZATION <<!

  ! Format Descriptors
  10 FORMAT(A)
  11 FORMAT(4X, A, 1X, A)

                       !>> FUNCTION BODY <<!
  WRITE(outputGeomUnit, 10) "%BLOCK SPECIES_POT"
  DO i=1, SIZE(cell%elementTable) - 1
    WRITE(outputGeomUnit, 11) &
      cell%elementTable(i)%elementName, &
      cell%elementTable(i)%psp%pspFileName
  END DO
  WRITE(outputGeomUnit, 10) "%END BLOCK SPECIES_POT"

END SUBROUTINE PrintPseudo


SUBROUTINE PrintOptimizeIons
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function prints out the list of ions with 1 if their positions were
!   optimized in a given dimension, and 0 if positions were fixed.
!
! CONDITIONS AND ASSUMPTIONS: Printed only when outputOptionGeom=.TRUE.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   08/27/2008  Created (Linda Hung)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                          !>> EXTERNAL VARIABLES <<!
                          !>> INTERNAL VARIABLES <<!
  INTEGER :: i,j   ! counters

  INTEGER, DIMENSION(3) :: writeFrozenIon

                          !>> INITIALIZATION <<!

  10 FORMAT(A)

                          !>> FUNCTION BODY << !

  ! Only print out block of a 'frozenion' array exists
  IF (ALLOCATED(frozenIon)) THEN

    WRITE(outputGeomUnit, 10) "%BLOCK ION_OPTIMIZATION"

    DO i=1, SIZE(cell%ionTable)
      DO j=1,3
        IF (frozenIon(i,j)) THEN
          writeFrozenIon(j) = 0   ! ion not to be optimized
        ELSE
          writeFrozenIon(j) = 1   ! ion to be optimized
        END IF
      END DO

      WRITE(outputGeomUnit, *) writeFrozenIon
    END DO

    WRITE(outputGeomUnit, 10) "%END BLOCK ION_OPTIMIZATION"

  END IF

END SUBROUTINE PrintOptimizeIons


SUBROUTINE PrintGeometry(cellRelaxFlag, Step, positions)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function prints out the entire geometry (*.*.ion) file, including
!   cell lattice, ion positions, and pseudopotential name.
!
! CONDITIONS AND ASSUMPTIONS: Printed only when outputOptionGeom=.TRUE.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   06/23/2008 Created (Linda Hung)
!   06/20/2014 Add function for printing geometry for cell relaxation (mohan)
!------------------------------------------------------------------------------
  IMPLICIT NONE

                           !>> EXTERNAL VARIABLES <<!

  LOGICAL, OPTIONAL :: cellRelaxFlag
  ! whether this printting geometry is after ion relaxation
  ! or after cell relaxation
  !
  INTEGER, INTENT(IN), OPTIONAL:: Step  
  ! the current step # of ion/cell relax
  ! 
  REAL(KIND=DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: positions
  ! fractional atomic coordinates 
  !

                           !>> INTERNAL VARIABLES <<!

  INTEGER :: fileStatus
  !
  CHARACTER(LEN=50) :: tmp
  !
  CHARACTER(LEN=50) :: geomFile
  !
                           !>> INITIALIZATION <<!

                           !>> FUNCTION BODY <<!


  IF (PRESENT(cellRelaxFlag)) THEN
    IF (PRESENT(Step)) THEN 
      WRITE(tmp,*) Step
      IF(cellRelaxFlag .eqv. .TRUE.) THEN
        geomFile = TRIM(outputSystemName) // ".cel." // TRIM(ADJUSTL(tmp)) // ".geom" 
      ELSE
        geomFile = TRIM(outputSystemName) // ".ion." // TRIM(ADJUSTL(tmp)) // ".geom"
      ENDIF
    ENDIF
  ELSE
    geomFile = TRIM(outputSystemName) // ".final.geom"
  ENDIF


  ! open the file,
  OPEN(unit=outputGeomUnit, access="sequential", action="write", &
       blank="null", delim="none", file=TRIM(geomFile), form="formatted", &
       iostat=fileStatus, pad="no", position="rewind", status="replace")

  ! in case we can't open the file
  IF (fileStatus/=0) THEN
    WRITE(*,*)'The file ',TRIM(geomFile),' had a problem on opening.'
    WRITE(*,*)'And I cannot even begin to imagine what it might be.'
    STOP
  ELSE
    ! print the lattice first in geometry file.
    CALL PrintLattice  
    ! print the atom positions in geometry file.
    IF(PRESENT(positions)) THEN 
      CALL PrintIons(positions)
    ELSE
      CALL PrintIons
    END IF
    ! print the atom species pseudopotential file name
    CALL PrintPseudo
    CALL PrintOptimizeIons
  END IF

END SUBROUTINE PrintGeometry


SUBROUTINE Print3DArray(array, fileName, arrayUnit, recOffset, scaling)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out a 3D array into the direct formatted file
!   associated with outputArrayUnit.  It will write to records starting from
!   record number recordOffset+1.  Processors write to the file one at a time -
!   sequentially, not simultaneously.
!  
! CONDITIONS AND ASSUMPTIONS: Only pass in portions of arrays that you want to
!   have written out.  For example, do not include boundary values when
!   printing array with Dirichlet boundary conditions.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   06/23/2008  Array printing made more modular.  See PrintDensity, 
!               PrintPotential, Dump__ for log notes on previous versions.
!               (Linda Hung)
!
!------------------------------------------------------------------------------
  USE MPI_Functions

  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(0:,0:,0:), INTENT(IN) :: array
  ! 3D real space array to be printed out
  !
  CHARACTER(LEN=*), INTENT(IN) :: fileName
  ! Name of output file
  !
  INTEGER, INTENT(IN) :: arrayUnit
  ! Number associated with file
  !
  INTEGER, INTENT(IN) :: recOffset    
  ! Value of 4th dimension (for 3d array, this is 0) so that
  ! records from previous 4th dimension don't get overwritten
  !
  REAL(KIND=DP), INTENT(IN) :: scaling
  ! For scaling indices of printed out array
  ! Input negative values for no scaling
  !

                       !>> INTERNAL VARIABLES <<!

  ! Because this is initialized in the declaration, the value will be 
  ! preserved from one execution of the program to another. 
  INTEGER :: &
#ifdef __USE_PARALLEL
    numProc, &
    nextRank, prevRank, &
#else
#endif
    startFlag, &
    ix,iy,iz, &       ! Counters for printing out
    indX, indY, indZ, &     ! Counters for array
    sizex, sizey, sizez, &  ! Size of non-padded array in each direction
    printx, printy, printz,&! Size of printed out array in each dimension
    fileStatus              ! Status to see if the file/mpi calls succeeded

#ifdef __USE_PARALLEL
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: &
    mpiStatus
#endif
  REAL(kind=dp) :: &
    fracX, fracY, fracZ, &  ! Place in array - not an integer
    diffX, diffY, diffZ, &  ! Difference frac? - INT(frac?)
    point1, point2          ! Used in interpolation

                              !>> INITIALIZATION <<!

  ! Format Descriptors
  10 FORMAT(A14 I12.1)
  11 FORMAT(E26.20e2)

  ! Only use non-padded array
  sizex=SIZE(array,1)
  sizey=SIZE(array,2)
#ifdef __USE_PARALLEL
  sizez=m3G
#else
  sizez=SIZE(array,3)
#endif

  ! Size of printed out array
  IF (scaling > 0._DP) THEN
    printx = INT(sizex*scaling)
    printy = INT(sizey*scaling)
    printz = INT(sizez*scaling)
  ELSE
    printx = sizex
    printy = sizey
    printz = sizez
  END IF

  ! Set up neighboring processor information
#ifdef __USE_PARALLEL
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numProc, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) STOP "**MPI_COMM_SIZE ERROR IN PRINTDENSITY**"
  prevRank = outputRank-1
  nextRank = outputRank+1
  IF (outputRank==0) prevRank=numProc-1
  IF (outputRank==numProc-1) nextRank=0
#endif

  startFlag = 0 ! Changes to 1 before a processor begins output

                          !>> FUNCTION BODY <<!

  ! First (or in serial case, only) processor opens file with replacement
  IF (outputRank==0) THEN
    startFlag = 1
    OPEN(UNIT=arrayUnit, ACCESS="direct", ACTION="write", BLANK="null", &
         file = trim(fileName), form="formatted", recl=27, status="replace",&
         iostat=fileStatus)
    IF (fileStatus/=0) THEN
      WRITE(*,*)'The OUTPUT array file ', trim(fileName), &
      ' had a problem opening.'
      STOP
    END IF

    ! Writes information about file length in a header (avoid repetitious
    ! header writing if 4d array calls for 4th dimension >1)
    IF (recOffset==0) THEN
      WRITE(arrayUnit, 10, REC=1) "x-dimension:", printx
      WRITE(arrayUnit, 10, REC=2) "y-dimension:", printy
      WRITE(arrayUnit, 10, REC=3) "z-dimension:", printz
    END IF
    WRITE(arrayUnit, 10, REC=4) "# of spins: ", recOffset + 1

  END IF

  zIndex: DO iz = 0, printz-1

    ! Check which z-index is needed
    IF (scaling > 0._DP) THEN
      fracZ = (REAL(sizez,kind=DP) / REAL(printz,kind=DP) &
               * (1._DP + 2._DP*REAL(iz,kind=DP)) - 1._DP) * 0.5_DP
      indZ = INT(fracZ)
      diffZ = fracZ - REAL(indZ,kind=DP)
      IF (indZ+1<n3Goff .OR. startFlag==2) CYCLE
    ELSE
      indZ = iz
      IF (indZ<n3Goff .OR. startFlag==2) CYCLE
    END IF

#ifdef __USE_PARALLEL
    ! If indZ is on a new processor - must set up write access
    IF (indZ>=n3Goff .AND. outputRank>0 .AND. startFlag==0) THEN
      CALL MPI_RECV(startFlag, 1, MPI_INTEGER, prevRank, 9, &
                    MPI_COMM_WORLD, mpiStatus, mpiErr)
      IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_RECV***"

      OPEN(UNIT=arrayUnit, ACCESS="direct", ACTION="write", BLANK="null", &
           file = trim(fileName), form="formatted", recl=27, status="old",&
           iostat=fileStatus)

      IF (fileStatus/=0) THEN
        WRITE(*,*)'The OUTPUT array file ', trim(fileName), &
        ' had a problem opening.'
        STOP
      END IF
    END IF

    ! Close file and allow next processor to access array output if this
    ! process is done
    IF (indZ>=n3Goff+SIZE(array,3)) THEN
      CLOSE(arrayUnit)
      IF ((outputRank/=numProc-1)  .OR. &
          (outputRank/=numProc/2)) THEN
        CALL MPI_SEND(startFlag, 1, MPI_INTEGER, nextRank, 9, &
                      MPI_COMM_WORLD, mpiErr)
        IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS SENDING MPI STUFF**"
      END IF
      startFlag = 2
      CYCLE
    END IF
#endif

    yIndex: DO iy = 0, printy-1
      IF (scaling > 0._DP) THEN
        fracY = (REAL(sizey,kind=DP) / REAL(printy,kind=DP) &
                 * (1._DP + 2._DP*REAL(iy,kind=DP)) - 1._DP) * 0.5_DP
        indY = INT(fracY)
        diffY = fracY - REAL(indY,kind=DP)
      END IF

      xIndex: DO ix = 0, printx-1

        IF (scaling < 0._DP) THEN ! Print out exactly what we have
          WRITE(arrayUnit, 11, &
                REC = recOffset*printx*printy*printz + iz*printx*printy &
                       + iy*printx + ix + 5) &
            array(ix,iy,iz-n3Goff)

        ELSE ! Get interpolated values on a smaller grid
          fracX = (REAL(sizex,kind=DP)/REAL(printx,kind=DP) &
                   * (1._DP + 2._DP*REAL(ix,kind=DP)) - 1._DP) * 0.5_DP
          indX = INT(fracX)
          diffX = fracX - REAL(indX,kind=DP)

          IF (indZ>=n3Goff .AND. indZ<n3Goff+SIZE(array,3)) &
            point1 = (array(indX,indY,indZ-n3Goff)*(1-diffX) &
                      + array(indX+1,indY,indZ-n3Goff)*diffX) &
                     *(1-diffY) &
                   + (array(indX,indY+1,indZ-n3Goff)*(1-diffX) &
                      + array(indX+1,indY+1,indZ-n3Goff)*diffX) &
                     *diffY

          indZ = indZ + 1

          IF (indZ>=n3Goff .AND. indZ<n3Goff+SIZE(array,3)) &
            point2 = (array(indX,indY,indZ-n3Goff)*(1-diffX) &
                        + array(indX+1,indY,indZ-n3Goff)*diffX) &
                     *(1-diffY) &
                   + (array(indX,indY+1,indZ-n3Goff)*(1-diffX) &
                      + array(indX+1,indY+1,indZ-n3Goff)*diffX) &
                     *diffY

          indZ = indZ - 1

#ifdef __USE_PARALLEL
          ! Communicate point2 value if indZ+1 is on a different processor
          ! from indZ.  Higher processor returns to beginning of loop.
          IF (indZ+1==n3Goff) THEN
            CALL MPI_SEND(point2, 1, MPI_REAL8, prevRank, 11, &
                          MPI_COMM_WORLD, mpiErr)
            IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS SENDING MPI STUFF**"
            CYCLE
          ELSE IF (indZ==n3Goff+SIZE(array,3)-1) THEN
            CALL MPI_RECV(point2, 1, MPI_REAL8, nextRank, 11, &
                            MPI_COMM_WORLD, mpiStatus, mpiErr)
            IF (mpiErr/=MPI_SUCCESS) STOP "***PROBLEMS WITH MPI_RECV***"
          END IF
#endif
          WRITE(arrayUnit, 11, &
                REC = recOffset*printx*printy*printz + iz*printx*printy &
                       + iy*printx + ix + 5) &
            point1*(1-diffZ) + point2*diffZ

        END IF
      END DO xIndex
    END DO yIndex
  END DO zIndex

  ! Final process closes file
#ifdef __USE_PARALLEL
  IF ((outputRank==numProc-1) .OR. &
      (outputRank==numProc/2)) &
    CLOSE(arrayUnit)
#else
  CLOSE(arrayUnit)
#endif

END SUBROUTINE Print3DArray


SUBROUTINE PrintDensity(density)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out the density into a file.  The file is named 
!   by the argument given, "systemName" followed by a number tag assigned 
!   sequentially every time the subroutine is run.
!  
! CONDITIONS AND ASSUMPTIONS: All processes with non-padded density should call
!   this since distributed density must be sent to head node to print out.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Will need to be upgraded if spin>1 is used
!   DOES THIS HAVE TO BE IN SI UNITS?
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/20/2003  Created (Greg Ho)
!   02/07/2007  Parallelized (LH)
!   01/29/2008  Writes to formatted direct output file (LH)
!   06/23/2008  3D array part separated for access by *.pot, *.den, and *.junk
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:,:) :: density
  ! Electron Density in real space to be printed
  !

                       !>> INTERNAL VARIABLES <<!

  INTEGER :: iSpin
  ! spin index
  !
  REAL(KIND=DP) :: fileSize  
  ! Expected density file size in MB (no scaling)
  !
  REAL(KIND=DP) :: scaling   
  ! Factor to give desired rho file size
  !
  CHARACTER(len=50) :: densityFile
  !

                       !>> INITIALIZATION <<!

  densityFile = TRIM(outputSystemName) // ".den"
  WRITE(outputUnit,*) "Print density file to ",densityFile


#ifdef __USE_PARALLEL
  fileSize = (REAL(SIZE(density,1),kind=DP) * REAL(SIZE(density,2),kind=DP) &
              * REAL(m3G,kind=DP) * REAL(SIZE(density,4),kind=DP) &
              + 4._DP) * 27._DP/1048576._DP
#else
  fileSize = (REAL(SIZE(density,1),kind=DP) * REAL(SIZE(density,2),kind=DP) &
              * REAL(SIZE(density,3),kind=DP) * REAL(SIZE(density,4),kind=DP) &
              + 4._DP) * 27._DP/1048576._DP
#endif

  IF (fileSize>maxMB) THEN
    WRITE(*,*) "WARNING: DENSITY FILE LARGER THAN ", maxMB, " MB."
    WRITE(*,*) "NOW RESIZING FROM ", fileSize, "MB TO SATISFY MAX MB." 
    scaling = (maxMB/fileSize)**(1._DP/3._DP)
  ELSE
    scaling = -1._DP ! Do not scale array
  END IF

                       !>> FUNCTION BODY <<!

  ! Everybody writes output to appropriate record number
  DO iSpin = 1, SIZE(density,4)
    CALL Print3DArray(density(:,:,:,iSpin), densityFile, densityUnit, iSpin-1, scaling)
  END DO ! spin index

  RETURN

END SUBROUTINE PrintDensity


SUBROUTINE PrintPotential(potential)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out the electronic potential into a file.  The file 
!   is named by the argument given, "fileName" followed by a number assigned 
!   sequentially every time the subroutine is run.  THE OUTPUT IS IN UNITS
!   eV * A^3!!
!  
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   ***SHOULD UNDERGO UNIT CONVERSION***
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/07/2007  Parallelized (LH)
!   01/30/2008  Writes to formatted direct file (LH)
!   06/23/2008  3D array part separated for access by *.pot, *.den, and *.junk
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(0:,0:,0:,:), INTENT(IN) :: &
    potential        ! Real space electronic potential on current proccesor

                       !>> INTERNAL VARIABLES <<!

  INTEGER :: iSpin

  REAL(kind=dp) :: &
    fileSize, &             ! Expected potential file size in MB (no scaling)
    scaling                 ! Factor to give desired rho file size

  CHARACTER(len=50) :: &
    potentialFile

                       !>> INITIALIZATION <<!

  potentialFile = TRIM(outputSystemName) // ".pot" 

#ifdef __USE_PARALLEL
  fileSize = (REAL(SIZE(potential,1),kind=DP)*REAL(SIZE(potential,2),kind=DP) &
              *REAL(m3G,kind=DP)*REAL(SIZE(potential,4),kind=DP) &
              + 4._DP) * 27._DP/1048576._DP
#else
  fileSize = (REAL(SIZE(potential,1),kind=DP)*REAL(SIZE(potential,2),kind=DP) &
              *REAL(SIZE(potential,3),kind=DP)*REAL(SIZE(potential,4),kind=DP) &
              + 4._DP) * 27._DP/1048576._DP
#endif

  IF (fileSize>maxMB) THEN
    WRITE(*,*) "WARNING: POTENTIAL FILE LARGER THAN ", maxMB, " MB."
    WRITE(*,*) "NOW RESIZING FROM ", fileSize, "MB TO SATISFY MAX MB." 
    scaling = (maxMB/fileSize)**(1._DP/3._DP)
  ELSE
    scaling = -1._DP ! Array not scaled
  END IF

                       !>> FUNCTION BODY <<!

  ! Everybody writes output to appropriate record number
  DO iSpin = 1, SIZE(potential,4)
    CALL Print3DArray(potential(:,:,:,iSpin), potentialFile, potentialUnit, &
                      iSpin-1, scaling)
  END DO ! spin index

END SUBROUTINE PrintPotential


SUBROUTINE DumpReal4D(array, systemName, increment)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out a grid quantity (density, potentials) into a 
!   file. 
!  
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   06/23/2008  3D array part separated for access by *.pot, *.den, and *.junk
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: &
    array

  CHARACTER(len=*), INTENT(IN) :: &
    systemName         ! The raw argument on the command line 

  LOGICAL, INTENT(IN) :: &
    increment

                       !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    count=0, &
    iSpin

  CHARACTER(len=30) :: &
    dumpNum
  
  CHARACTER(len=50) :: &
    dumpFile
                       !>> INITIALIZATION <<!

  IF(increment) count = count + 1
  WRITE(dumpNum, *) count
  dumpFile = TRIM(systemName) // "." // TRIM(ADJUSTL(dumpNum)) // ".junk"

                       !>> FUNCTION BODY <<!

  ! Everybody writes output to appropriate record number
  DO iSpin = 1, SIZE(array,4)
    CALL Print3DArray(array(:,:,:,iSpin), dumpFile, junkUnit, iSpin-1, -1._DP)
  END DO ! spin index

END SUBROUTINE DumpReal4D


SUBROUTINE DumpReal3D(array, systemName, increment)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out a grid quantity (density, potentials) into a 
!   file.  
!  
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   06/23/2008  3D array part separated for access by *.pot, *.den, and *.junk
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    array

  CHARACTER(len=*), INTENT(IN) :: &
    systemName         ! The raw argument on the command line 

  LOGICAL, INTENT(IN) :: &
    increment

                       !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    count=0

  CHARACTER(len=30) :: &
    dumpNum
  
  CHARACTER(len=50) :: &
    dumpFile
                       !>> INITIALIZATION <<!

  IF(increment) count = count + 1
  WRITE(dumpNum, *) count
  dumpFile = TRIM(systemName) // "." // TRIM(ADJUSTL(dumpNum)) // ".junk"

                       !>> FUNCTION BODY <<!

  ! Everybody writes output to appropriate record number
  CALL Print3DArray(array(:,:,:), dumpFile, junkUnit, 0, -1._DP)

END SUBROUTINE DumpReal3D


SUBROUTINE DumpComplex(grid, systemName, preserve)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out a complex grid quantity into a file.  
!  
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Should extend to 4-D someday, current not parallelized.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  COMPLEX(kind=DP), DIMENSION(:,:,:) :: &
    grid

  CHARACTER(len=*) :: &
    systemName         ! The raw argument on the command line 

                       !>> INTERNAL VARIABLES <<!

  ! Because this is initialized in the declaration, the value will be 
  ! preserved from one call to the routine to the next. 
  INTEGER :: &
    i, j, k, &
    sizex, sizey, &
    sizez, &
    count = 0, &     ! The number of times this subroutine has been executed
    fileStatus       ! Status integer to see if the file opening succeeded

  LOGICAL, OPTIONAL :: &
    preserve 

  CHARACTER(len=30) :: &
    potentialNum
  
  CHARACTER(len=50) :: &
    potentialFile

                       !>> INITIALIZATION <<!
  count = count + 1
  WRITE(potentialNum, *) count

  IF(PRESENT(preserve)) THEN
    IF(preserve) THEN
      potentialFile = TRIM(systemName) // "dump" // &
        TRIM(ADJUSTL(potentialNum)) // ".den"
    ELSE
      potentialFile = TRIM(systemName) // "dump" // ".den"
    END IF
  ELSE
    potentialFile = TRIM(systemName) // "dump" // ".den"
  END IF

  sizex = SIZE(grid, 1)
  sizey = SIZE(grid, 2)
  sizez = SIZE(grid, 3)

  ! Format Descriptors
  11 FORMAT(E26.20e2)

                       !>> FUNCTION BODY <<!

  OPEN(unit=junkUnit, access="sequential", action="write", blank="null", &
       file = trim(potentialFile), form="formatted", &
       iostat=fileStatus, position="rewind", status="replace")
  IF (fileStatus/=0) THEN
    WRITE(*,*)'The OUTPUT potential file ', trim(potentialFile), &
    ' had a problem opening.'
    STOP
  END IF  

  DO k = 1, sizez
    DO j = 1, sizey

        WRITE(junkUnit,'(1000(E20.10e2,1X))') (ABS(grid(i,j,k)), i=1,sizex)

    END DO
    WRITE(junkUnit,*) " "
  END DO

  CLOSE(junkUnit)

END SUBROUTINE DumpComplex


SUBROUTINE DumpLogical(grid, systemName, preserve)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out a logical grid quantity into a file.  
!  
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Should extend to 4-D someday, current not parallelized.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  LOGICAL, DIMENSION(:,:,:) :: &
    grid

  CHARACTER(len=*) :: &
    systemName         ! The raw argument on the command line 

                       !>> INTERNAL VARIABLES <<!

  ! Because this is initialized in the declaration, the value will be 
  ! preserved from one call to the routine to the next. 
  INTEGER :: &
    i, j, k, &
    sizex, sizey, &
    sizez, &
    count = 0, &     ! The number of times this subroutine has been executed
    fileStatus       ! Status integer to see if the file opening succeeded

  LOGICAL, OPTIONAL :: &
    preserve 

  CHARACTER(len=30) :: &
    potentialNum
  
  CHARACTER(len=50) :: &
    potentialFile

                       !>> INITIALIZATION <<!
  count = count + 1
  WRITE(potentialNum, *) count

  IF(PRESENT(preserve)) THEN
    IF(preserve) THEN
      potentialFile = TRIM(systemName) // "dump" // &
        TRIM(ADJUSTL(potentialNum)) // ".den"
    ELSE
      potentialFile = TRIM(systemName) // "dump" // ".den"
    END IF
  ELSE
    potentialFile = TRIM(systemName) // "dump" // ".den"
  END IF

  sizex = SIZE(grid, 1)
  sizey = SIZE(grid, 2)
  sizez = SIZE(grid, 3)

                       !>> FUNCTION BODY <<!

  OPEN(unit=junkUnit, access="sequential", action="write", blank="null", &
       file = trim(potentialFile), form="formatted", &
       iostat=fileStatus, position="rewind", status="replace")
  IF (fileStatus/=0) THEN
    WRITE(*,*)'The OUTPUT potential file ', trim(potentialFile), &
    ' had a problem opening.'
    STOP
  END IF  

  DO k = 1, sizez
    DO j = 1, sizey
        WRITE(junkUnit,'(1000(L1,1X))') (grid(i,j,k), i=1,sizex)
    END DO
    WRITE(junkUnit,*) " "
  END DO

  CLOSE(junkUnit)

END SUBROUTINE DumpLogical


SUBROUTINE PrintEwald(ewaldSum, eta, realCutoff, recipCutoff, realError, &
                      recipError, ewaldComponents, ewaldTimes)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out information related to the Ewald summation.
!  
! CONDITIONS AND ASSUMPTIONS: Only called if outputEwald=.TRUE.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/10/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), INTENT(IN) :: &
    ewaldSum, &       ! The final answer
    eta, &            ! Parameter which controls how much of the calculation
                      ! is done in real space vs reciprocal space.  
    realCutoff, &     ! The optimum real-space cutoff radius (calculated)
    recipCutoff, &    ! The optimum reciprocal-space cutoff radius (calculated)
    realError, &      ! The maximum error in real space
    recipError        ! The maximum error in reciprocal space

  REAL(kind=DP), DIMENSION(:), INTENT(IN) :: &
    ewaldComponents   ! Array that holds all the values from the various
                      ! components

  REAL(kind=DP), DIMENSION(:), INTENT(IN) :: &
    ewaldTimes        ! Array holding the execution time of each procedure
  
                       !>> INTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(SIZE(ewaldComponents)) :: &
    ewaldEnergiesAU   ! Ewald energies in AU
 
                       !>> INITIALIZATION <<!
  ewaldEnergiesAU = ewaldComponents * hartreeToeV

  ! Format descriptors
  10 FORMAT(A25, 3X, F15.12, 1X, A3)
  11 FORMAT(A28, 1X, ES18.12, 1X, A2)
  13 FORMAT(A22, 4X, F20.12, 1X, "eV", 3X, "time: ", E8.3, 1X, "s")


                        !>> FUNCTION BODY <<!
  IF(outputEwald) THEN
    WRITE(outputUnit,*) " "
    WRITE(outputUnit,*) "      ********THE EWALD SUMMATION********"
    WRITE(outputUnit,10) "eta:                     ", eta / bohr, "1/A"
    WRITE(outputUnit,10) "Real-Space Cutoff:       ", realCutoff * bohr, "A"
    WRITE(outputUnit,10) "Reciprocal-Space Cutoff: ", recipCutoff / bohr, &
      "1/A"  
    WRITE(outputUnit,*) " "

    WRITE(outputUnit,11) "Max Real-Space Error:       ", &
                         realError * hartreeToeV, "eV" 
    WRITE(outputUnit,11) "Max Reciprocal-Space Error: ", &
                         recipError * hartreeToeV, "eV"
    WRITE(outputUnit,*) " "
 
    WRITE(outputUnit,13) "  Real-Space Term:      ", ewaldEnergiesAU(1), &
                         ewaldTimes(1)
    WRITE(outputUnit,13) "  Recip-Space Term:     ", ewaldEnergiesAU(3), &
                         ewaldTimes(3)
    WRITE(outputUnit,13) "  Real-Space Self Term: ", ewaldEnergiesAU(2), &
                         ewaldTimes(2)
    WRITE(outputUnit,13) "  Ave. Potential Term:  ", ewaldEnergiesAU(4), &
                         ewaldTimes(4)

    WRITE(outputUnit,*) &
      "+______________________________________________________________"
    WRITE(outputUnit,13) "Total Ewald Energy:        ", &
      ewaldSum * hartreeToeV, SUM(ewaldTimes)

    WRITE(outputUnit,*) "       ******************************************  "
  END IF
END SUBROUTINE PrintEwald




SUBROUTINE PrintTransitionState
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out the file necessary to be output to interface
!   with CINEB.f90, the transition state finder.
!
! CONDITIONS AND ASSUMPTIONS: Call only if outputTransitionState = .TRUE.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
    

                     !>> EXTERNAL VARIABLES <<!
                    !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    i
                       !>> INITIALIZATION <<!

  ! Format Descriptors
  10 FORMAT(E25.15)
  11 FORMAT(3E25.15)
                       !>> FUNCTION BODY <<!

  ! Write out the energy
  WRITE(transUnit, 10) energy(1)*hartreeToeV
  WRITE(transUnit, *)

  ! Write out the forces
  DO i = 1, SIZE(forceIon, 1)
    WRITE(transUnit, 11) forceIon(i,:,1)*hartreeToeV/bohr
  END DO
  WRITE(transUnit, *)

  ! Write out the stress
  WRITE(transUnit, 11) stress* auToGPa
  WRITE(transUnit, *)

  ! Write out the fractional ion coordinates
  DO i=1, SIZE(cell%ionTable)
    WRITE(transUnit, 11) cell%ionTable(outputIonTable(i))%coord
  END DO

END SUBROUTINE PrintTransitionState


SUBROUTINE WrtOut(unitOut, msg, nextLine)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Output messages to the unit (6 is the screen), this subroutine take care
!   of the parallel: only the master node writes message.
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

                         !>> EXTERNAL VARIABLES <<!

   CHARACTER (LEN=*), INTENT(IN) :: msg       ! Message to output
   INTEGER, INTENT(IN)  :: unitOut            ! Unit number for writing out
   LOGICAL, INTENT(IN), OPTIONAL :: nextLine  ! Advance to next line?

                         !>> FUNCTION BODY <<!

   IF (outputRank==0) THEN
     IF (PRESENT(nextLine)) THEN
       IF (nextLine) THEN
         WRITE(unitOut,'(A)') TRIM(msg)
       ELSE
         WRITE(unitOut,'(A)',ADVANCE="NO") TRIM(msg)
       ENDIF
     ELSE
       WRITE(unitOut,'(A)') TRIM(msg)
     ENDIF
     flush(unitOut)
   ENDIF

END SUBROUTINE WrtOut


SUBROUTINE CleanOutput
!------------------------------------------------------------------------------
! DESCRIPTION:
!  we clean all the allocated arrays in this module
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Created by Chen Huang May/21/2009
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

  DEALLOCATE(outputIonTable)
  
  return

END SUBROUTINE CleanOutput

! mohan add 22-10-12
SUBROUTINE QUIT(message)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012 Created by Mohan Chen
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=500),OPTIONAL :: message
  !
#ifdef __USE_PARALLEL
  INTEGER :: mpiErr   ! Error status from initializing, finalizing MPI
  CALL MPI_FINALIZE(mpiErr) ! End the MPI abilities of the program  
#endif
  !
  CALL WrtOut(6," ")
  CALL WrtOut(6," -------------------------------- WARNING ------------------------------------")
  CALL WrtOut(6," PROFESS Stop for the following reason:")
  IF(PRESENT(message)) THEN
    CALL WrtOut(6, message)
  ENDIF
  CALL WrtOut(6," Find out the reason and restart calculation again.")
  CALL WrtOut(6," -----------------------------------------------------------------------------")

  !
  CALL StopClock('PROFESS')

  ! pritn the time information into outputUnit
  IF(outputRank==0) THEN
    CALL printClock(' ',outputUnit)
  ENDIF

  STOP

END SUBROUTINE Quit


SUBROUTINE PrintDensityClearly(rhoR, iter)
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
#ifdef __USE_PARALLEL
  INCLUDE 'mpif.h'
#endif
  !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: rhoR ! Electronic density in real-space
  INTEGER, OPTIONAL :: iter

  !
  !>> INTERNAL VARIABLES <<!
  !
  INTEGER :: ii, jj, kk
  INTEGER :: NumSpin
#ifdef __USE_PARALLEL
  INTEGER :: myRank, mySize, mpiErr
  INTEGER :: zn, isp
#endif
  CHARACTER(LEN=100) :: message
  CHARACTER(LEN=100) :: message2
  CHARACTER(LEN=100) :: message3
  

  NumSpin=SIZE(rhoR,4) 

  IF(PRESENT(iter)) THEN
    WRITE(message,'(I8)') iter
    WRITE(message,'(A,A,A)') 'den.',TRIM(ADJUSTL(message)),'.dat'

    WRITE(message2,'(I8)') iter
    WRITE(message2,'(A,A,A)') 'den2.',TRIM(ADJUSTL(message2)),'.dat'

    WRITE(message3,'(I8)') iter
    WRITE(message3,'(A,A,A)') 'den_diff.',TRIM(ADJUSTL(message3)),'.dat'
  ELSE
    WRITE(message,'(A)') "den.dat"
    WRITE(message2,'(A)') "den2.dat"
    WRITE(message3,'(A)') "den_diff.dat"
  ENDIF

#ifndef __USE_PARALLEL

  OPEN(unit=2000,FILE=message,ACTION='WRITE')
  WRITE(2000,*) 'TITLE="OFdensity"'
  WRITE(2000,*) 'VARIABLES = "X"  "Y"  "Z"  "OFDENSITY"'
  WRITE(2000,*) 'ZONE I= ', size(rhoR,1),' J=  ', size(rhoR,2),' K= ',size(rhoR,3),' F=POINT'
  DO kk=lbound(rhoR,3),ubound(rhoR,3)
    DO jj=lbound(rhoR,2),ubound(rhoR,2)
      DO ii=lbound(rhoR,1),ubound(rhoR,1)
        WRITE(2000,'(3I4,ES20.12)') ii,jj,kk,rhoR(ii,jj,kk,1)
      ENDDO
    ENDDO
  ENDDO
  CLOSE(2000)

#else

  
  !! DEBUG DUMP DENSITY 
  CALL MPI_Comm_size(MPI_COMM_WORLD, mysize, mpiErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpiErr)

  ! write titles for tecplot
  IF(myrank == 0) then
      OPEN(unit=1999,FILE=message,ACTION='WRITE')
      WRITE(1999,*) 'TITLE = "OFdensity"'
      WRITE(1999,*) 'VARIABLES = "X"  "Y"  "Z"  "OFDENSITY"'
      WRITE(1999,*) 'ZONE I= ', size(rhoR,1),' J=  ', size(rhoR,2),' K= ',mysize*size(rhoR,3),' F=POINT'
      CLOSE(1999)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD, mpiErr)

 
  zn = 0
  DO isp = 0, mysize-1
    IF(myrank == isp) THEN
      OPEN(unit=1999,FILE=message,ACTION='WRITE',POSITION='APPEND')
      DO kk=lbound(rhoR,3),ubound(rhoR,3)
        DO jj=lbound(rhoR,2),ubound(rhoR,2)
          DO ii=lbound(rhoR,1),ubound(rhoR,1)
            WRITE(1999,'(3I4,ES20.12)') ii,jj,kk+zn,rhoR(ii,jj,kk,1)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(1999)
  
      IF(NumSpin == 2) THEN
        OPEN(unit=1995,FILE=message2,ACTION='WRITE',POSITION='APPEND')
          DO kk=lbound(rhoR,3),ubound(rhoR,3)
            DO jj=lbound(rhoR,2),ubound(rhoR,2)
              DO ii=lbound(rhoR,1),ubound(rhoR,1)
                WRITE(1995,'(3I4,ES20.12)') ii,jj,kk+zn,rhoR(ii,jj,kk,2)
              ENDDO
            ENDDO
           ENDDO
        CLOSE(1995)
           
        OPEN(unit=1993,FILE=message3,ACTION='WRITE',POSITION='APPEND')
        DO kk=lbound(rhoR,3),ubound(rhoR,3)
          DO jj=lbound(rhoR,2),ubound(rhoR,2)
            DO ii=lbound(rhoR,1),ubound(rhoR,1)
              WRITE(1993,'(3I4,ES20.12)') ii,jj,kk+zn,(rhoR(ii,jj,kk,1)-rhoR(ii,jj,kk,2))/(rhoR(ii,jj,kk,1)+rhoR(ii,jj,kk,2))
            ENDDO
          ENDDO
        ENDDO
        CLOSE(1993)
      ENDIF
      zn = zn + ubound(rhoR,3)+1
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpiErr)
    CALL MPI_BCAST(zn,1,mpi_integer,isp,MPI_COMM_WORLD,mpiErr)
  ENDDO

!  WRITE(message,'('' (NE2) Max Rho            : '',F10.8)') MAXVAL(rhoR)
!  CALL WrtOut(6,message)
!  WRITE(message,'('' (NE2) Min Rho            : '',F10.8)') MINVAL(rhoR)
!  CALL WrtOut(6,message)
  !! END DEBUG

#endif   


END SUBROUTINE PrintDensityClearly

END MODULE Output


! Use to debug, print out the title, 
SUBROUTINE TITLE(message)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012 Created by Mohan Chen
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP
  USE OutputFiles, ONLY: outputUnit
  USE MPI_Functions, ONLY: rankGlobal
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: message

! If you want to print TITLE, screen the 'RETURN'
  RETURN
  IF(rankGlobal==0) THEN
    WRITE(6,'(" ==> ",A)') message
    WRITE(outputUnit,'(" ==> ",A)') message
  ENDIF

END SUBROUTINE TITLE

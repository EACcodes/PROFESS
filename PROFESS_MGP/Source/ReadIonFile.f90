MODULE ReadIonFile
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE ReadIonFile 
!     |_SUBROUTINE ReadGeometry
!       |_FUNCTION Locate
!     |_SUBROUTINE ReadPseudo
!     |_SUBROUTINE ReadDensity
!
! DESCRIPTION:
!   Subroutines that read input files from the outside world, and sets up
!   the module SYSTEM to reflect what is in the files.  This module is 
!   responsible for files such as the .ion files, the pseudopotential files,
!   and the density files.
!
!   The reason this is a separate module from InitializeInputs is that the
!   .inpt file that ReadInputFile handles is very specific to our OFDFT 
!   code.  The other input files, however, are not ... in fact, we have tried
!   to model them after the CASTEP input files of the same type.  If desired,
!   these modules could be useful for other programs.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Perhaps PASS things like cell, rhoR, and bound instead of using SYSTEM.
! 
! REFERENCES:
!   [1] M.C. Payne, M.P. Teter, D.C. Allan, T.A. Arias, and J.D. Joannopoulos,
!         Reviews of Modern Physics, Vol. 64, No. 4, p. 1045, October 1992.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/15/2003  Consolidated all of these functions in to one module (Greg Ho)
!   12/18/2003  Commenting, formatting, commenting, formatting! (GSH)
!
!------------------------------------------------------------------------------

  USE Constants, ONLY: DP
  USE Constants, ONLY: systemNameLen ! The allowed length of the system name
  USE Constants, ONLY: pi            ! 3.1415...
  USE Constants, ONLY: bohr          ! The bohr radius
  USE Constants, ONLY: hartreeToeV   ! hartree to eV

  USE Output, ONLY: WrtOut
  USE Output, ONLY: QUIT
  USE OutputFiles, ONLY: errorUnit
  USE OutputFiles, ONLY: outputUnit
  USE OutputFiles, ONLY: Error

  USE CellInfo, ONLY : cell   ! the cell, contains lattice vectors and ions
  USE CellInfo, ONLY : RefreshLattice ! refresh the cell information

  USE Sys, ONLY: frozenIon
  
  USE MathFunctions, ONLY : Volume
  ! Function that gets the volume of a cell
  
  USE MathSplines, ONLY: spline_cubic_set
  ! used for cubic spline, this routine will set up 
  ! spline for each pseudopotential

  USE AtomicDensity, ONLY : atomicDensityFile
  USE KEDF_DenDec, ONLY: atomicCoreFile ! mohan add 2013-07-31

  IMPLICIT NONE

  INTEGER, PARAMETER :: inputUnit = 5
  ! Make an input unit varable for input files.  

  CHARACTER(len=systemNameLen + 4), DIMENSION(:), ALLOCATABLE :: &
  pseudoFile         ! Another temporary array.
                     ! This one stores the names of the psp files
                     ! to be read by the ReadPseudo routine. 
                     ! Alloctated in ReadGeometry and deallocated in 
                     ! ReadPseudo

CONTAINS


SUBROUTINE ReadGeometry(geometryFile, defaultPseudo)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine opens and reads the input file that describes the system.
!   It sets the box size and shape, number and position of ions in fractional
!   coordinates.
!
!   This procedure is written to read input files in CASTEP format, but most
!   CASTEP keywords will be ignored. Only the cell size and ionic positions are
!   read, and units are assumed to be Angstroms.
!   WARNING: Specifying a different unit will cause this parser to bail out. 
!   This is not a bug, it's a feature, as this program will assume distances 
!   to be in Angstroms and would output garbage if they were not.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   The ions are sorted by a simple sort.  If this turns out to be slow, 
!   which it shouldn't, Quicksort could be implemented.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/15/2003  File created.  (VLL)
!   10/21/2003  The routine reads cell size and ions properly (VLL)
!   10/22/2003  Pseudopotential file names have been added (VLL)
!   11/13/2003  Added loop to sort the ions by type (VLL)
!   11/21/2003  Updated modules used and prettified.  Added loop to calculate
!               ionTypeIndexTable (GSH)
!   12/15/2003  Changed many comments (GSH)
!   01/16/2003  Changed so that internally the units are in atomic units. (GSH)
!   03/05/2003  Added a small hack in here so that we can output the ions in
!               the order they were read in.
!   05/11/2005  Fixed the bug that it crashes with only one ion.
!
!------------------------------------------------------------------------------
  USE CONSTANTS, ONLY: AtomicMass
  USE CellInfo, ONLY: ion
  USE MathFunctions, ONLY: Cross
  USE MathFunctions, ONLY: Vecmul
  USE MathFunctions, ONLY: Inverse
  USE Output, ONLY: outputIonTable
  ! Table of ion pointers, unsorted, for output

  IMPLICIT NONE
  

                           !>> EXTERNAL VARIABLES <<!

  CHARACTER(len=*), INTENT(IN) :: geometryFile
  CHARACTER(len=*), INTENT(IN) :: defaultPseudo

                           !>> INTERNAL VARIABLES <<! 

  CHARACTER(len=6) :: firstArg
  ! First argument on the ion line, i.e. ion type.
  !
  CHARACTER(len=34) :: command
  ! Temporary storage for the current line.
  !
  CHARACTER(len=6), DIMENSION(:), ALLOCATABLE :: ionElement
  ! Temporary array to store the various kinds of ions declared in the 
  ! geometry file. Used to call the corresponding pseudo-potential files.
  !
  INTEGER :: &
    numFrozen, &
    i,j, &                     ! dummy counter
    numIonType, &              ! Number of types of ions in the system
    numIon, &                  ! Number of ions in the system
    fileStatus, &              ! Checks that file operations are going smoothly
    lineCounter, &             ! Counts the total number of lines in the file.
    tempElement                ! Temporary element for initializing 
                               ! ionTypeIndexTable

  INTEGER, DIMENSION(3) :: readFrozenIon
  ! Used for reading input on whether ion position is fixed

  LOGICAL :: &
    foundLattice = .FALSE., &  ! Control flag for lattice parameters
    foundIons = .FALSE., &     ! Control flag for ions
    ionExist, &                ! Flag for counting atom types in the system.
    ionCoordInCart = .FALSE., &! Flag for projecting ionic coordinates.
    ionCoordInBohr = .FALSE.   ! Flag for projecting ionic coordinates.

  REAL(KIND=DP) :: &
    a, b, c, &                 ! Cell lattice dimensions. 
    alpha, beta, gamma, & 
    tmp_vector(3)

  REAL(KIND=DP), DIMENSION(3,3) :: converter
  ! Matrix for ion projection.
  !
  TYPE(ion), DIMENSION(:), ALLOCATABLE :: ionStorage
  ! Temporary table for ions
  !
  TYPE(ion) :: transfer
  ! dummy ion for switching table elements.
  ! we make it a pointer because of a compiler 
  ! bug in ifc7.
  !
  CHARACTER(len=500) :: message
  !

                          !>> INITIALIZATION <<!   

  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           READ GEOMETRY"

                           !>> FUNCTION BODY <<!
  ! I'm going to do something that's kind of ugly and wasteful: I'll count the
  ! lines in the geometry file and allocate the ionstorage table based on that.
  ! Then I'll count allocated lines in ionstorage and allocate the ion table
  ! (Permanent one, declared in OFDFT) based on that. Why? Because I have no
  ! way of knowing in advance how many ions are in the file, and I can't fill a
  ! table that has not been allocated. Unless we're talking millions of atoms,
  ! the whole business will take negligible time anyways, so it's no big deal.
  ! Also, the ionElement table is now allocated this way too.  Heehee.
  ! If anyone has a better way, feel free to modify this routine.
  WRITE(message,*) " "
  CALL WrtOut(6,message)
  WRITE(message,*) '(Read Geometry) from file ',geometryFile
  CALL WrtOut(6,message)


  OPEN(unit=inputUnit, access="sequential", action="read", blank="null", &
       delim="none", file=geometryFile, form="formatted", iostat=fileStatus, &
       pad="no", position="rewind",  status="old")

  IF (fileStatus/=0) THEN
    WRITE(message,*)' Could not open ',TRIM(geometryFile),'.'
    CALL QUIT(message)
  END IF

  lineCounter = 0
  DO
    ! One should note that before compiler PATHSCALE 2.5, we used here the
    ! format '(A)'. However, this stopped working. So now we use *, which 
    ! means that WHITESPACES are not read. So if you were expecting here
    ! the number of lines, you will be disappointed -- linecounter here only
    ! contains the number of NON-BLANK lines! This is okay for our purposes,
    ! however, so we won't lose too much sleep over it.
    READ(inputUnit, *, iostat=fileStatus) command
    ! End of file reached
    IF (fileStatus<0) EXIT 
    lineCounter =lineCounter + 1
  END DO

  ! Reset the input file.
  REWIND inputUnit


  ! Allocate the ion table
  ALLOCATE(ionStorage(lineCounter), stat=fileStatus)
  IF (filestatus/=0) THEN
    WRITE(message,*)'Error allocating the temporary ion table. Leaving.'
    CALL WrtOut(6,message)
    STOP
  END IF
  ionStorage%elementID = 1


  ! Then allocate the table in which we'll keep the ion elements.  This is
  ! very wasteful, but whatever.  It'll be deallocated at the end of this
  ! subroutine so it doesn't matter.  The other alternative would be to 
  ! actually count the number of ions here.
  ALLOCATE(ionElement(lineCounter), stat=fileStatus)
  IF (filestatus/=0) THEN
    WRITE(errorUnit,*)'Error allocating the temporary ion label. Leaving.'
    STOP
  END IF 
  ionElement = ''

  DO
    ! Read the next line.
    READ(inputUnit,*,iostat=fileStatus) firstArg
    
    IF (fileStatus<0) THEN
      EXIT ! We have encountered the end of the file, time to leave.
    END IF

    IF(TRIM(firstArg(1:1))=='#') CYCLE
    BACKSPACE inputUnit
    READ(inputUnit,*,iostat=fileStatus) firstArg, command


    IF (firstArg=="%BLOCK") THEN
      SELECT CASE (TRIM(command))

!--------------------------------------------------------------------------------------
! INPUT THE CELL
!--------------------------------------------------------------------------------------
        CASE ("LATTICE_CART")
          foundLattice = .TRUE.
          READ(inputUnit,*,iostat=fileStatus) cell%cellReal
          IF (fileStatus/=0) THEN
            WRITE(errorUnit,*)'Error reading the CART cell size. Leaving.'
            STOP
          END IF

!--------------------------------------------------------------------------------------
! INPUT THE ANGLE IS ANOTHER OPTION
!--------------------------------------------------------------------------------------
        CASE ("LATTICE_ABC")
          foundLattice = .TRUE.
          READ(inputUnit,*,iostat=fileStatus) a, b, c, alpha, beta, gamma
          IF (fileStatus/=0) THEN
            WRITE(errorUnit,*)'Error reading the ABC cell size. Leaving.'
            STOP
          END IF

          cell%lengthA = a
          cell%lengthB = b
          cell%lengthC = c
          cell%angleAlpha = alpha  ! unit is angles
          cell%angleBeta = beta
          cell%angleGamma = gamma

          cell%cellReal = 0._DP
          alpha = alpha * pi / 180._DP       ! Converting angles to radians.
          beta = beta * pi / 180._DP
          gamma = gamma * pi /180._DP
          cell%cellReal(1,1) = a              ! First vector is along x-axis.
          cell%cellReal(1,2) = b * COS(gamma) ! Second vector is in x-y plane.
          cell%cellReal(2,2) = b * SIN(gamma) ! Projecting on y-axis.
          cell%cellReal(1,3) = c * COS(beta)  ! Projecting on x-axis.
          ! These last two are based on vector b dot vector c = b c cos(alpha) 
          ! and norm conservation, respectively. I checked them, they work.
          cell%cellReal(2,3) = c * (COS(alpha) - COS(beta) * COS(gamma)) / &
                               SIN(gamma)
          cell%cellReal(3,3) = &
            SQRT(c**2-cell%cellReal(1,3)**2-cell%cellReal(2,3)**2)


!--------------------------------------------------------------------------------------
! INPUT THE ATOM COORDINATES
!--------------------------------------------------------------------------------------
        CASE ("POSITIONS_FRAC", "POSITIONS_CART", "POSITIONS_BOHR")
          IF (foundIons) EXIT
          foundIons = .TRUE.

          ! If the ion coordinates are absolute, we will need to convert them 
          ! to the fractional system later on. This marker will remind us.
          IF (TRIM(command(1:13))=="POSITIONS_ABS" .OR. &
              TRIM(command(1:14))=="POSITIONS_CART") &
              ionCoordInCart = .TRUE.

          IF (TRIM(command(1:13))=="POSITIONS_BOHR" ) & 
              ionCoordInBohr = .TRUE.
           
          numIon = 0
          numIonType = 0
          DO
          ! Assuming this new ion is of unknown type, we read its type and 
          ! we check against all other known types.  If we know it, we 
          ! reference it.
            READ(inputUnit,*) firstArg
            
            ! a comment line
            if (firstArg(1:1) == '#' ) cycle
            if (firstArg(1:4) == '%END' ) exit

            ! This line had a valid ion on it, keep going.
            numIon = numIon + 1


            ! ionExist = .FALSE. means it's a new type of ion
            ionExist = .FALSE.                   
            DO i=1, numIonType
              IF (firstArg==ionElement(i)) THEN  
                ionStorage(numIon)%elementID = i 
                ionExist = .TRUE.
                EXIT
              END IF
            END DO

            ! Now, if we have not found it then it's really new. So we put it 
            ! in the array.
            IF (.NOT.ionExist) THEN              
              numIonType = numIonType + 1
              ionElement(numIonType)=firstArg
              ionStorage(numIon)%elementID = numIonType
            END IF

            ! Reading the ion type made us advance a line.
            BACKSPACE inputUnit 
            ! Re-reading, this time coords also
            READ(inputUnit,*, iostat=fileStatus) firstArg,& 
                 ionStorage(numIon)%coord     

            ! We've reached the end of the list of ions.  In fact, that last 
            ! one wasn't even an ion.
            IF (fileStatus/=0) THEN 
              numIon = numIon - 1   
              EXIT
            END IF

          END DO

          ! So now we've found all the ions and they are recorded, in ion 
          ! format, in our temporary storage table ionStorage. We also know 
          ! we have numion of them.  Their element matches the index of the 
          ! ionelement table, for psp purposes.  Now, allocate the ion table
          ALLOCATE(cell%ionTable(numIon), stat=fileStatus) 
          IF (filestatus/=0) THEN                     
            WRITE(errorUnit,*)'Error allocating the main ion table. Leaving.'
            STOP
          END IF

          ALLOCATE(outputIonTable(numIon), stat=fileStatus)
          IF (filestatus/=0) THEN                     
            WRITE(errorUnit,*)'Error allocating the output ion table. Leaving.'
            STOP
          END IF

          ! We transfer our storage table and free the memory it was hoarding.
          cell%ionTable = ionStorage(1:numIon) 
          DEALLOCATE(ionStorage) 

          ! Initialize the array of pseudopotential filenames so that each 
          ! entry has the pseudopotential suffix inside.
          ALLOCATE(pseudoFile(numIonType), stat=fileStatus)
          ! for atomic density
          ALLOCATE(atomicDensityFile(numIonType))
          ! for atomic core density; mohan add 2013-07-31
          ALLOCATE(atomicCoreFile(numIonType))

          ! set the number of ion types and the total ion number
          ! in the cell
          cell%numIonType = numIonType
          cell%numIon = numIon
          !WRITE(*,*) "numIonType=", numIonType
          !WRITE(*,*) "numIon=", numIon

          IF (filestatus/=0) THEN                     
            WRITE(errorUnit,*)'Error allocating pseudofile array. Leaving.'
            STOP
          END IF
          pseudoFile = defaultPseudo

        ! This case uses numIonType and therefore depends on the previous case 
        ! being executed before. That means the ion definition has to precede 
        ! the pseudopotential file list.
        CASE ("SPECIES_POT")

          DO
            ! command is the psp filename
            READ(inputUnit,*) firstArg, command 

            if (firstArg(1:1) == '#' ) cycle
            IF (firstArg(1:4)=="%END") EXIT          
            
            ! If we haven't reached the end, we're going to have to look at 
            ! our table of elements and match the file.  Note that if there is 
            ! no match we do nothing at all.
            DO i=1, numIonType                  
              IF (firstArg==ionElement(i)) THEN 
                pseudoFile(i) = command     
                EXIT                        
              END IF                        
            END DO
          END DO

        ! the input core density file for each element
        CASE ("SPECIES_CORE")
          DO
            ! command is the psp filename
            READ(inputUnit,*) firstArg, command 

            if (firstArg(1:1) == '#' ) cycle
            IF (firstArg(1:4)=="%END") EXIT          
            
            ! If we haven't reached the end, we're going to have to look at 
            ! our table of elements and match the file.  Note that if there is 
            ! no match we do nothing at all.
            DO i=1, numIonType                  
              IF (firstArg==ionElement(i)) THEN 
                atomicCoreFile(i) = command     
                EXIT                        
              END IF                        
            END DO
          END DO

        ! the pseudopotential file list.
        CASE ("SPECIES_RHOA")
          DO
            READ(inputUnit,*) firstArg, command
            IF (firstArg(1:1) == '#' ) CYCLE
            IF (firstArg(1:4)=="%END") EXIT
            ! If we haven't reached the end, we're going to have to look at
            ! our table of elements and match the file.  Note that if there is
            ! no match we do nothing at all.
            DO i=1, numIonType
              IF (firstArg==ionElement(i)) THEN
                atomicDensityFile(i) = command
                EXIT
              END IF
            END DO
          END DO

        ! This (optional) case assumes that the total number of ions is known.
        ! If a dimension for an ion has value 0, it is not optimized (frozen).
        ! Value of 1 (or actually all nonzero values), ion is optimized.
        ! In the absence of this block, all ion positions will be optimized.
        CASE ("ION_OPTIMIZATION")
          ALLOCATE(frozenIon(numIon,3))
          frozenIon = .FALSE.
          numFrozen = 0
          DO i=1, numIon
            READ(inputUnit,*, iostat=fileStatus) readFrozenIon
            IF (fileStatus/=0) THEN
              CALL WrtOut(6,'Error reading ION_OPTIMIZATION. Leaving.')
              STOP
            END IF
            WHERE (readFrozenIon==0) frozenIon(i,:)=.TRUE.
            IF (MAXVAL(readFrozenIon)==0 ) THEN 
              numFrozen = numFrozen + 1
            ENDIF
          END DO
          WRITE(message,'(a,I10)') & 
            " (Read Geometry) Frozen atom number      : ",numFrozen
          CALL WrtOut(6,message) 

        CASE DEFAULT
          WRITE(message,*) 'Invalid BLOCK keyword in ion geometry file! key =', firstArg
          CALL WrtOut(6,message)
          STOP

      END SELECT
    END IF
  END DO

  WRITE(message,'('' (geometry) Atom Species                 : '', I10)') numIonType
  CALL WrtOut(6,message)
  WRITE(message,'('' (geometry) Atom Number                  : '', I10)') numIon
  CALL WrtOut(6,message)


  IF (.NOT.foundLattice) THEN
    WRITE(errorUnit,*)'Could not find lattice parameter. Leaving.'
    STOP
  END IF
  IF (.NOT.foundIons) THEN
    WRITE(errorUnit,*)'Could not find any ionic coordinates. Leaving'
    STOP
  END IF

  ! Cleaning: We may have overcounted the number of types by reading a "%END"
  ! string as an ion type. If so, that type would be the last one and we need
  ! to remove it from numIonType.
  IF (ionElement(numIonType)(1:4)=="%END") numIonType = numIonType - 1
  
  ! If the coordinates are absolutes, we convert to fractional by projection.
  ! ATTENTION:  We must do this BEFORE we convert the cell lattice vectors
  !             to bohr!
  IF (ionCoordInCart) THEN
    converter=Inverse(cell%cellReal)
    DO i=1, numIon
      ! This is technically a projection along the cell axes. Cute, eh?
      cell%ionTable(i) % coord = Vecmul(converter,cell%ionTable(i)%coord)
    END DO 
  END IF
 
  !! We have coordinates in bohr
  IF (ionCoordInBohr) THEN
    converter=Inverse(cell%cellReal)
    DO i=1, numIon
      ! This is technically a projection along the cell axes. Cute, eh?
      cell%ionTable(i) % coord = Vecmul(converter,cell%ionTable(i)%coord * bohr )
    END DO 
  END IF

  ! steven modified, in case that ion position is out of the primitive cell, move all images
  ! into the cell:
  DO i=1, numIon
    cell%ionTable(i)%coord(1) = cell%ionTable(i)%coord(1) - FLOOR(cell%ionTable(i)%coord(1))
    cell%ionTable(i)%coord(2) = cell%ionTable(i)%coord(2) - FLOOR(cell%ionTable(i)%coord(2))
    cell%ionTable(i)%coord(3) = cell%ionTable(i)%coord(3) - FLOOR(cell%ionTable(i)%coord(3))
  EndDo


  !-------------------------------------------------------------------------------------------------
  ! Now, we want the lattice vectors in ATOMIC units but they were given in
  ! Angstroms.  Convert them here.
  CALL RefreshLattice( cell%cellReal / bohr )

  !-------------------------------------------------------------------------------------------------



  ! Initialize the ion table (For output)
  DO i=1, numIon
    outputIonTable(i) = i 
  END DO

  ! Now we need to make sure that the atoms are sorted by type.
  ! This is sort of the simplest sorting algorithms known.  If this turns
  ! out to be slow, we can implement quicksort
  ! Only sort if we have more than one ion!
  IF(numIon > 1) THEN
    i=2
    DO
      IF (cell%ionTable(i)%elementID < cell%ionTable(i-1)%elementID) THEN
        ! swap the ions
        transfer = cell%ionTable(i)
        cell%ionTable(i) = cell%ionTable(i-1)
        cell%ionTable(i-1) = transfer

        ! Store the swap information for printing later.  We need to swap the 
        ! ORIGINAL positions to get the result we want.
        j = outputIonTable(Locate(i, outputIonTable))
        outputIonTable(Locate(i, outputIonTable)) = &
          outputIonTable(Locate(i-1, outputIonTable))     
        outputIonTable(Locate(i-1, outputIonTable)) = j     
        
        i=i-1
        IF (i<2) i=2
      ELSE
        i=i+1
      END IF
      IF (i>numIon) EXIT
    END DO
  END IF

  ! Now we allocate and initialize elementTable, which contains the index of
  ! ionTable each different type of ion starts at, the atomic name,
  ! pseudopotential, charge, and mass.  Assume that ionTable is sorted.
  ALLOCATE(cell%elementTable(numIonType + 1), stat=fileStatus) 
  IF (filestatus/=0) THEN                     
    WRITE(errorUnit,*)'Error allocating the ion type number table.  Leaving.'
    STOP
  END IF

  ! Start out with a non-element
  tempElement = 0
  j = 1
  DO i = 1, numIon
    IF (cell%ionTable(i)%elementID .NE. tempElement) THEN
      cell%elementTable(j)%firstIonID = i
      j = j + 1      
      tempElement = cell%ionTable(i)%elementID
    ENDIF
  END DO

  ! Assign the last element of elementTable%firstIonID to numIon+1, so that
  ! we can use this array to loop over all ions of a given type
  cell%elementTable(numIonType + 1)%firstIonID = numIon + 1

  ! Save element names
  DO i=1, numIonType
    cell%elementTable(i)%elementName = ionElement(i)
    ! Input the first three characters of ionElement, update by mohan 09-27-12
    cell%elementTable(i)%mass = AtomicMass(ionElement(i)(1:3))
  END DO

  ! At this point, any psp file name that has not yet been overriden by the
  ! geometry file needs to be changed from the default (just a suffix) to a
  ! full name using the type of ion specified by the user.
  DO i=1, numIonType
    IF (pseudoFile(i)==defaultPseudo) THEN
      pseudoFile(i) = TRIM(ionElement(i))//defaultPseudo
    END IF
  END DO

  CLOSE (inputUnit, status="keep") ! Closing the file, we're done.

  DEALLOCATE(ionElement)

  !==============================================
  ! Test if cell is in right-hand coordinate
  tmp_vector = Cross(cell%cellReal(:,1),cell%cellReal(:,2))
  IF ( dot_product(tmp_vector, cell%cellReal(:,3)) < 0.d0 ) THEN
    CALL WrtOut(6,"Box is not in right-hand coordinate, Change your cell definition in .ion file. code stop")
    STOP
  ENDIF

  ! At this point, if everything went fine, we know:
  ! numIon and numIonType, plus their associated tables ionTable and 
  ! ionElement. The coordinates of all the ions in fractional form.
  ! The complete geometry of the cell, contained in cellReal, as three vectors.


  CONTAINS
  
  FUNCTION Locate(x, a)
  !----------------------------------------------------------------------------
  ! DESCRIPTION:
  ! This internal subroutine returns the position of the first occurance of 
  ! integer x in integer array a. 
  !
  ! CONDITIONS AND ASSUMPTIONS: 
  ! x exists in a, and a starts with an index of 1
  !
  ! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
  !
  !----------------------------------------------------------------------------
  ! REVISION LOG:
  !   03/4/2004  Created (GSH)
  !----------------------------------------------------------------------------
  IMPLICIT NONE 
                 !>> EXTERNAL VARIABLES <<!
  INTEGER, DIMENSION(:), INTENT(IN) :: a
  ! The array
  !
  INTEGER, INTENT(IN) :: x
  ! The integer
  !
  INTEGER :: Locate
  ! The answer
  !

                 !>> INTERNAL VARIABLES <<!
  INTEGER :: i
  ! counter
  !

                 !>> INITIALIZATION <<!
                 !>> FUNCTION BODY <<!
  DO i = 1, SIZE(a)
    IF(a(i) == x) THEN
      Locate = i
      EXIT
    END IF
  END DO

  END FUNCTION Locate

END SUBROUTINE ReadGeometry


SUBROUTINE ReadPseudo
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine reads all the pseudopotential files necessary to run the 
!   system and fills out the corresponding pseudopotential table and 
!   parameters.
!
!   It also figures out the charges on each ion from the pseudopotential and
!   initializes it in ionTable.  To do this, we consider the following from
!   [1].  All pseudopotentials die off to Z/r at large r.  Fourier transforming
!   this, we get that at low G, the pseudopotential is dominated by a
!   4*pi*Z/G**2 term.  The value of the pseudopotential at this G consists
!   of a coulombic part (which follows the above expression) and a 
!   non-coulombic part, which is very close to the value of the pseudopotential
!   at G=0 (since there is no coulombic term at G=0, by construction).  
!   Thus, we take the value of the pseudopotential at the smallest nonzero
!   G that we have, and we subtract the G=0 term.  Call this V.  Then
!   Z = - V * G**2 / 4 / pi (see [1])
!
! CONDITIONS AND ASSUMPTIONS:
!   We're using the NEW CASSTEP format for .recpot files.  Pretty much, the
!   pseudopotential file contains a big comment block at the beginning that
!   ends with END COMMENT.  The next line is 2 numbers that correspond to
!   the version number of the PseudoPotential, which we ignore.  The next line
!   corresponds to the max G value (in A^-1) and then the pseudopotential
!   values, in ev*AU^3/A^3 are displayed, three to a line.  The first one is 
!   of course the G=0 value, and the last one is the G=maxG value.  WE 
!   ASSUME THAT THE VALUES ARE EVENLY SPACED IN G, AND ALSO THAT TOTAL NUMBER
!   OF VALUES IS A NUMBER DIVISIBLE BY 3 (so that the last line will have 3
!   values).
!
!   Also, since these are local pseudopotentials, we assume that everything
!   is in one channel...i.e., everything after the comments all constitute
!   one channel all the way to the end.
!
!   There's also a one-line terminator "1000" at the end of the file.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/22/2003  File created.  (VLL)
!   10/23/2003  Routine completed and tested for one ion type. Works. (VLL)
!   11/22/2003  Cosmetic changes, added some modules (GSH)
!   12/02/2003  I changed the Pseudopotential format to the new CASSTEP format.
!               Also, the ionic charge is now calculated.  Actually, I did a
!               whole bunch of stuff.  Hard to explain everything.(GSH)
!   12/03/2003  Removed "PAD=NO" from the OPEN parameters.  Didn't work
!               otherwise.(GSH)
!   02/25/2004  Changed so that it doesn't need pseudotable or 
!               allpseudopotvalues anymore.
!   05/13/2005  Added .realpot type support.  .realpot format is similar except
!               its in logarithmic form and for now the charge is included in
!               the file.
!------------------------------------------------------------------------------
  USE CellInfo, ONLY : pseudoPot ! A pseudopotential type

  IMPLICIT NONE
                         !>> EXTERNAL VARIABLES <<!

                         !>> INTERNAL VARIABLES <<! 
  CHARACTER(len=80) :: line ! a temporary storage for a line of the file
              
  INTEGER :: &
    numIonType, &     ! The number of types of ions in the system
    ionCharge, &      ! The charge on an ion
    i, j, ii, &       ! dummy counter.
    allocateStatus, & ! Control integer.
    fileStatus, &     ! Control integer. Checks if file opening goes right.
    readStatus        ! Control integer, Checks whether reads go right

  TYPE(pseudoPot), POINTER :: psp               ! pointer to a temporary pseudopotential

  INTEGER :: position
  INTEGER :: tempInt

  REAL(kind=DP) :: tempReal
  CHARACTER (len=500) :: message

                          !>> INITIALIZATION <<!


  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           READ PSEUDOPOTENTIAL"

  numIonType = SIZE(cell%elementTable) - 1

  cell%elementTable%charge = 0.d0
  cell%elementTable%chargeTot = 0.d0

                          !>> FUNCTION BODY <<!
  ! (1) Now, for each type of ion, initialize the pseudopotentials
  DO i = 1, numIonType

    ! (1.1) Allocate a new memory location for our temporary pseudopotenial pointer
    ALLOCATE(psp)

    ! (1.2) Open the pseudopotential file
    OPEN(UNIT=inputUnit, ACCESS="sequential", ACTION="read", BLANK="null", &
         DELIM="none", FILE=pseudoFile(i), FORM="formatted",& 
         IOSTAT=fileStatus, POSITION="rewind",  STATUS="old")
  
    IF (fileStatus/=0) THEN
      WRITE(message,*)' Could not open ',TRIM(pseudoFile(i)),'.'
      CALL Quit(message)
    END IF
  
    ! (1.3) Store the file name we got it from 
    psp%pspFileName = TRIM(pseudoFile(i))
   
    ! Figure out what kind of file it is
    position = INDEX(psp%pspFileName, ".", back=.TRUE.)

    IF(TRIM(psp%pspFileName(position:)) == ".recpot") THEN
      psp%type = 1     ! type .recpot
    ELSE
      STOP "***PSP FILE TYPE NOT SUPPORTED***"
    END IF

    ! (1.4) 
    IF(psp%type == 1 .OR. psp%type == 2) THEN

      ! (1.4.1) Read until we find an END COMMENT statement
      line = ''
      DO WHILE(INDEX(line, 'END COMMENT') == 0)
        READ(inputUnit,'(A)', IOSTAT=readStatus) line
        IF(readStatus < 0) THEN
          WRITE(message,*) " Could not find 'END COMMENT' line in ", psp%pspFileName
          CALL QUIT(message)
        ENDIF
      END DO

      ! (1.4.2) Also skip the version information
      READ(inputUnit, '(A)') line


      ! (1.4.3)
      IF(psp%type == 1) THEN

        ! Now read in the max G value (in A^-1)
        READ(inputUnit, *) psp%maxG

        ! Convert maxG to AU^-1
        psp%maxG = psp%maxG * bohr

      ELSE IF(psp%type == 2) THEN

        READ(inputUnit, *) psp%maxR

      ELSE
        STOP "***PSP FILE TYPE NOT SUPPORTED***"
      END IF

    ! Now just parse through the file and figure out how many lines there are
    ! WE ASSUME HERE THAT THERE ARE 3 ENTRIES PER LINE, INCLUDING THE LAST ONE!
    ! We also assume that since these are LOCAL pseudopotentials, there is only
    ! one channel (ie, everything after the comments constitute one file and
    ! so we can read to the end.)
      psp%numPoints = 0
      DO
        READ(inputUnit, *, IOSTAT=readStatus) line
        IF(readStatus < 0) EXIT
        ! Every line of psp file contains 3 values!
        psp%numPoints = psp%numPoints + 3
      END DO
      ! since the last line is a terminator, we read 3 points too far.
      psp%numPoints = psp%numPoints - 3

      ! Check for a proper terminator
      IF(line(1:4) /= "1000") STOP "***Pseudopotential lacks the proper &
                                    &terminator***"

      ! Allocate the pseudopotential values for this pseudopotential
      ALLOCATE(psp%potValues(psp%numPoints), stat=allocateStatus)

      IF (allocateStatus/=0) THEN
        WRITE(errorUnit,*)'Error allocating space for values of pseudopotential.'
        STOP
      END IF

      ! Now that we know how much memory we needed to allocate, rewind the file
      ! and go back to the beginning.  Actually read the values this time.
      REWIND(inputUnit) 
    
      ! Read until we find an END COMMENT statement
      line = ''
      DO WHILE(INDEX(line, 'END COMMENT') == 0)
        READ(inputUnit,'(a)') line
      END DO

      ! Also skip the version information and max G value, which we already
      ! read in.
      DO j=1,2
        READ(inputUnit, '(a)') line
      END DO

      ! Now read in the rest of the pseudopotential (terminator, of course, is
      ! discarded)
      READ (inputUnit, *, IOSTAT=readStatus) psp%potValues
      IF(readStatus < 0) THEN
        WRITE(message,*) "Last line in the .recpot file &
        should contain the full three values"
        CALL Error(6, message)
      ENDIF

      If(psp%potValues(size(psp%potValues))==1000) then
          write(*,*) 'Warning: 2nd last line of psp fine misses one number! Corrected to 0.00 here!'
          psp%potValues(size(psp%potValues)) = 0.d0
      Endif

      IF(psp%type == 1) THEN

        ! Adjust to be in atomic units
        psp%potValues = psp%potValues / (bohr**3 * hartreeToeV)

        ! Now figure out the charge of the ion based on the pseudopotential
        ! for all ions of this type. (see equation in description) [1]
        ionCharge = NINT( &
          - (psp%potValues(2) - psp%potValues(1)) * &
          (psp%maxG / REAL(psp%numPoints-1,DP))**2 / &
          (4.0_DP * pi))

        WRITE(message,'(A,I10)')' (ReadPseudo) ionCharge                  : ',ionCharge ! (Shin, mohan add 10-01-12)
        CALL WrtOut(6,message) !


        !
        ! To spline our V(q) + 4*pi*Charge/q**2, Yes! we are NOT going to spline v(q) itself
        ! Simply because, there is NO way for us to spline psp%potValues(1) (which is finite)
        ! and psp%potValues(2) which is going to -infinity
        !
        ! We are goint to store the 2nd derivative at each q-points to potDD array
        ! (defined psp data strucutre) after doing spline
        ! these 2nd derivatives will later be used by us to interpolate thoese pseudopotential
        !

        !    WRITE(*,'(A,i2,A,i3)')  'Doing Spline for pseudopotential of ion type:',&
        !                            i, ',  of ionCharge=', ionCharge

        ! Allocate the pseudopotential values for this pseudopotential
        ALLOCATE(psp%potDD(psp%numPoints), stat=allocateStatus)
        IF (allocateStatus/=0) THEN
          WRITE(*,*)'Error allocating space for 2nd derivatives of pseudopotential.'
          STOP
        END IF

    
        ! preparing knots for spline
        ALLOCATE(psp%t(psp%numPoints), stat=allocateStatus)
        IF (allocateStatus/=0) THEN
          WRITE(*,*)'Error allocating space for t(:). knots, in ReadInputFiles.f90, stop'
          STOP
        END IF
        DO ii = 1, psp%numPoints 
          psp%t(ii) = REAL(ii,kind=DP)   ! preparing knots for spline, we 
        ENDDO                      ! use index of each q-point as knots

        !preparing  vqS(:)
        ALLOCATE(psp%vqS(psp%numPoints), stat=allocateStatus)
        IF (allocateStatus/=0) THEN
          WRITE(*,*)'Error allocating space for vqS in ReadInputFiles.f90, stop'
          STOP
        END IF
        psp%vqS(1) = psp%potValues(1)
        psp%vqS(2:psp%numPoints) = psp%potValues(2:psp%numPoints) + &
                               4._DP * pi * ionCharge / & 
                               ( (psp%t(2:psp%numPoints)-1._DP)/REAL(psp%numPoints-1, kind=DP)*psp%maxG )**2._DP
        ! set up the spline
        CALL spline_cubic_set ( psp%numPoints, psp%t, psp%vqS, 1, 0._DP, 1, 0._DP, psp%potDD )

!    WRITE(*,*) psp%potDD
    
!    WRITE(*,'(A,i2,A,i3)')  ' --- Splined pseudopotential of ion type:',&
!                            i, ',  of ionCharge=', ionCharge
    ! end of our spline ------------------------:

      ELSE IF(psp%type == 2) THEN
        ! The real space PSP is already in atomic units

        ! Now figure out the charge of the ion based on the pseudopotential
        ! for all ions of this type. (see equation in description) [1]
        ionCharge = -NINT(psp%potvalues(psp%numPoints) * psp%maxR)

        ! this is the pseudopotential format for ABINIT
      END IF

    ELSE IF(psp%type == 3) THEN
      line = ''

      ! skip first line (comment line)
      READ(inputUnit,'(a)') line

      ! Third line contains the charge
      READ(inputUnit,*) tempReal, tempReal
      ionCharge = NINT(tempReal)

      ! Fourth line includes the # of points
      READ(inputUnit,*) psp%numPoints,  psp%numPoints,  psp%numPoints,  &
                        psp%numPoints,psp%numPoints  

      ALLOCATE(psp%potValues(psp%numPoints), stat=allocateStatus)

      ! skip four lines
      DO j=1,4
        READ(inputUnit,'(a)') line
      END DO

      ! Get maxR
      READ (inputUnit,*) tempInt, tempReal, psp%potValues(1) 
      READ (inputUnit,*) tempInt, tempReal, psp%potValues(2)

      psp%maxR = tempReal*REAL(psp%numPoints-1,kind=DP)

      ! now read in the pseudopotential values
      DO j=3, psp%numPoints
        READ (inputUnit,*) tempInt, tempReal, psp%potValues(j) 
      END DO

    END IF

    ! Assign this charge to all the ions of this type in ionTable.
    ! Also create a pointer to the corresponding pseudopotential.
    psp%charge = ionCharge
    cell%elementTable(i)%charge = ionCharge
    cell%elementTable(i)%chargeTot = ionCharge * &
                                     REAL(cell%elementTable(i+1)%firstIonID - &
                                      cell%elementTable(i)%firstIonID, kind=DP)

    cell%elementTable(i)%psp => psp

    CLOSE (inputUnit) 
    NULLIFY(psp)

  END DO


  cell%numEle = SUM(cell%elementTable%chargeTot)
  WRITE(outputUnit,*) "Total electron number in ReadPseudo is ", cell%numEle
  
  ! We won't need the pseudopotential file names anymore.
  DEALLOCATE(pseudoFile)
 

END SUBROUTINE ReadPseudo


SUBROUTINE ReadDensity(densityFile, xRepeat,yRepeat,zRepeat)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine reads the density as it is stored, in real space, in the
!   specified file, and stores it in the density table provided. If the size
!   of the file grid matches what has, at this point, been determined to be the
!   density grid size in real space, good. Otherwise, values are recalculated
!   by interpolation to match the current grid size.
!   Density file format is: m1G, m2G, m3G, numSpin, followed by the 
!   density at each grid point (without any mention of coordinates.)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/28/2003  File created.  (VLL)
!   10/29/2003  Density read from file and interpolated if necessary (VLL)
!   11/13/2003  Adapted the extrapolation part to 4-dimensional grid (VLL)
!   11/22/2003  Cosmetic changes, added a module (GSH)
!   03/10/2004  Removed numEle from SYSTEM
!   01/30/2008  Input file is now formatted, direct access (LH)
!
!------------------------------------------------------------------------------

  USE MPI_Functions
  USE SYS, ONLY: rhoR ! The density in real space
  USE SYS, ONLY: interior
  USE CellInfo, ONLY: m3G, n3Goff

  IMPLICIT NONE

                     !>> EXTERNAL VARIABLES <<!
  CHARACTER(LEN=*), INTENT(IN) :: densityFile
  INTEGER, INTENT(IN) :: &
    xRepeat, &   ! Default is 1 - reads in density file only once
    yRepeat, &   ! Larger values indicate repetitions in either the x, y, or z
    zRepeat      ! dimension

                     !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    ix,iy,iz,iSpin, &   ! dummy counters.
    indX, indY, indZ, & ! Indexes of the file grid that are directly below
                        ! a given current grid point in that grid. 
    fileStatus, &       ! Control integer. Checks whether files exist.
    fileSpin, &         ! The number of spin-densities contained in the file.
    sX, sY, sZ, &       ! Shifts in indexes for points close to boundaries.
    xLen, yLen, zLen, & ! Number of points in the file in x, y and z resp
    numX, numY, numZ, & ! Number of points in the rhoR in the x, y, z dirs.
    numS, &             ! Number of spins in the system (1 or 2)
    numEle              ! Number of electrons in the system

  REAL(KIND=DP) :: &
#ifdef __USE_PARALLEL
    tempLocInt, &       ! Value of integral on local processor
#endif
    integral, &         ! integral of rho over the interpolation.
    fracX, fracY, &     ! Fractional coordinates of the
      fracZ, &          ! current real sp. grid points in the file grid.
    diffX, diffY, &     ! Difference between fracX and indX, etc.
      diffZ, &
    data1, data2, &     ! Values from density file used in interpolation
    data3, data4, &
    point1, point2      ! intermediate interpolation values

  CHARACTER(LEN=14) :: comment ! Comment about file length - junked

                           !>> INITIALIZATION <<!
  
  rhoR = 0._DP

  numX = SIZE(rhoR, 1)
  numY = SIZE(rhoR, 2)
  numZ = SIZE(rhoR, 3)

  numS = SIZE(rhoR, 4)
  numEle = NINT(SUM(cell%elementTable%chargeTot))

  ! Format Descriptors
  10 FORMAT(A14 I12.1)
  11 FORMAT(E26.20e2)

  rhoR = 0._DP
  integral = 0._DP

                           !>> FUNCTION BODY <<!

  OPEN(unit=inputUnit, access="direct", action="read", blank="null", &
       delim="none", file=densityFile, form="formatted",& 
       iostat=fileStatus, pad="no", status="old", recl=27)
  IF (fileStatus/=0) THEN
    WRITE(errorUnit,*)'Could not open ',TRIM(densityFile),'.'
    WRITE(errorUnit,*)'Please make sure this density file does exist.'
    WRITE(*,*)'Please make sure this density file does exist.'
    WRITE(outputUnit,*)'Please make sure this density file does exist.'
    STOP
  ELSE IF(fileStatus==0) THEN 
    WRITE(outputUnit,'(A,A)') " Read in the denisty file from ",densityFile
  END IF

  ! Read header for density file length
  READ (inputUnit, 10, REC=1, iostat=fileStatus) comment, xLen
  IF (fileStatus/=0) THEN
    WRITE(errorUnit,*)'Density file ',TRIM(densityFile),' is not properly ',&
      'formatted for x. Leaving.'
    STOP
  END IF
  READ (inputUnit, 10, REC=2, iostat=fileStatus) comment, yLen
  IF (fileStatus/=0) THEN
    WRITE(errorUnit,*)'Density file ',TRIM(densityFile),' is not properly ',&
      'formatted for y. Leaving.'
    STOP
  END IF
  READ (inputUnit, 10, REC=3, iostat=fileStatus) comment, zLen
  IF (fileStatus/=0) THEN
    WRITE(errorUnit,*)'Density file ',TRIM(densityFile),' is not properly ',&
      'formatted for z. Leaving.'
    STOP
  END IF
  READ (inputUnit, 10, REC=4, iostat=fileStatus) comment, fileSpin
  IF (fileStatus/=0) THEN
    WRITE(errorUnit,*)'Density file ',TRIM(densityFile),' is not properly ',&
      'formatted for spin. Leaving.'
    STOP
  END IF

  IF (numS /= fileSpin) THEN
    WRITE(errorUnit,*)'The current job has either too many or too few ',&
      'spin-densities compared with ',TRIM(densityFile), '. Leaving.'
    STOP
  END IF

  IF (xLen*xRepeat==numX.AND.yLen*yRepeat==numY.AND.zLen*zRepeat==m3G) THEN


    ! Check the grid sizes, and if they match, copy. 
    ! We also need to renormalize, in case the incoming density was in 
    ! units other than A^-3 (often, rho will be expressed in bohr^-3.)
    ! To do this, we sum the density over the whole grid, and turn that into
    ! an integral over the cell volume with the ratio.
    DO iSpin = 1, numS
      DO iz = 0, numZ-1
        indZ = MOD((n3Goff+iz),zLen)
        DO iy = 0, numY-1
          indY = MOD(iy,yLen)
          DO ix = 0, numX-1
            indX = MOD(ix,xLen)
            READ(inputUnit, 11, &
                 REC = (5 + indX + xLen*(indY + yLen*(indZ + zLen*(iSpin-1)))))&
                 rhoR(ix,iy,iz,iSpin)
          END DO
        END DO
      END DO
    END DO

  ELSE
  ! Now, here comes the problem: there is no guarantee that the file grid 
  ! matches the current real space grid. And if they don't, we'll have to 
  ! interpolate.
    DO iSpin = 1, numS
      DO iz = 0, numZ-1
        fracZ = (REAL(zLen*zRepeat,kind=DP) / REAL(m3G,kind=DP) &
                 * (1._DP + 2._DP*REAL(iz+n3Goff,kind=DP)) - 1._DP) * 0.5_DP
        fracZ = MODULO(fracZ,REAL(zLen,kind=DP))
        indZ = INT(fracZ)
        diffZ = fracZ - REAL(indZ,kind=DP)
        IF (indZ==zLen-1) THEN
          sZ=zLen ! Shift for periodic wrapping if necessary
        ELSE
          sZ = 0                    ! By default there should be no shift.
        END IF
        DO iy = 0, numY-1
          fracY = (REAL(yLen*yRepeat,kind=DP) / REAL(numY,kind=DP) &
                   * (1._DP + 2._DP*REAL(iy,kind=DP)) - 1._DP) * 0.5_DP
          fracY = MODULO(fracY,REAL(yLen,kind=DP))
          indY = INT(fracY)
          diffY = fracY - REAL(indY,kind=DP)
          IF (indY==yLen-1) THEN
            sY=yLen
          ELSE
            sY = 0
          END IF
          DO ix = 0, numX-1
            fracX = (REAL(xLen*xRepeat,kind=DP) / REAL(numX,kind=DP) &
                     * (1._DP + 2._DP*REAL(ix,kind=DP)) - 1._DP) * 0.5_DP
            fracX = MODULO(fracX,REAL(xLen,kind=DP))
            indX = INT(fracX)
            diffX = fracX - REAL(indX,kind=DP)
            IF (indX==xLen-1) THEN
              sX=xLen
            ELSE
              sX = 0
            END IF
         ! Here is how this ugly bunch works: imagine the cell. Now put 
         ! the real space grid from the file in the cell. Now superimpose 
         ! the real space grid for this calculation to it. Each point of 
         ! this grid can be viewed as being inside a cube formed by eight 
         ! points of the file grid (a pen and paper may help.)  We get 
         ! the density at these eight points from the file and calculate 
         ! rho at the ninth by doing a weighed average. Since the formula 
         ! was a bit long, I split it into two averages over faces, then 
         ! averaged along the z direction. sX, sY and sZ are used for the 
         ! upper x, y and z boundaries. Take the example of sX: if fracX+1 
         ! is between xLen and xLen+1, then the two densities we want
         ! to average are at xlen and... 1! So we need to remove xLen 
         ! from the index we look at to ensure we are considering the 
         ! right point. Same for sY and sZ.
            READ(inputUnit, 11, REC=(iSpin-1)*xLen*yLen*zLen+indZ*xLen*yLen &
                                     +indY*xLen+indX+5) data1
            READ(inputUnit, 11, REC=(iSpin-1)*xLen*yLen*zLen+indZ*xLen*yLen &
                                     +indY*xLen+indX-sX+6) data2
            READ(inputUnit, 11, REC=(iSpin-1)*xLen*yLen*zLen+indZ*xLen*yLen &
                                     +(indY+1-sY)*xLen+indX+5) data3
            READ(inputUnit, 11, REC=(iSpin-1)*xLen*yLen*zLen+indZ*xLen*yLen &
                                     +(indY+1-sY)*xLen+indX-sX+6) data4
            point1 = (data1*(1-diffX) + data2*diffX) * (1-diffY) &
                   + (data3*(1-diffX) + data4*diffX) * diffY
            indZ = indZ + 1 - sZ
            READ(inputUnit, 11, REC=(iSpin-1)*xLen*yLen*zLen+indZ*xLen*yLen &
                                     +indY*xLen+indX+5) data1
            READ(inputUnit, 11, REC=(iSpin-1)*xLen*yLen*zLen+indZ*xLen*yLen &
                                     +indY*xLen+indX-sX+6) data2
            READ(inputUnit, 11, REC=(iSpin-1)*xLen*yLen*zLen+indZ*xLen*yLen &
                                     +(indY+1-sY)*xLen+indX+5) data3
            READ(inputUnit, 11, REC=(iSpin-1)*xLen*yLen*zLen+indZ*xLen*yLen &
                                     +(indY+1-sY)*xLen+indX-sX+6) data4
            point2 = (data1*(1-diffX) + data2*diffX) * (1-diffY) &
                   + (data3*(1-diffX) + data4*diffX) * diffY
            indZ = indZ - 1 + sZ
            rhoR(ix,iy,iz,iSpin) = point1*(1-diffZ) + point2*diffZ
          END DO
        END DO
      END DO
    END DO

  END IF

  integral = SUM(rhoR, interior)

#ifdef __USE_PARALLEL
  tempLocInt = integral
  CALL MPI_ALLREDUCE(tempLocInt, integral, 1, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "**PROBLEMS WITH MPI_ALLREDUCE IN READDENSITY**"
#endif

  ! The reason we did the integral is that now we need to renormalize rho.
  integral = REAL(numEle,kind=DP) * REAL(numX,kind=DP) * REAL(numY,kind=DP) &
             * REAL(m3G, kind=DP) * REAL(numS, kind=DP) &
             / (integral * cell%vol)

  WHERE(interior)
    rhoR = rhoR * integral
  END WHERE

  ! Closing the file, we're done.
  CLOSE (inputUnit, status="keep")

END SUBROUTINE ReadDensity


END MODULE ReadIonFile

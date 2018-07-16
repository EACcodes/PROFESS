PROGRAM Convert
!------------------------------------------------------------------------------
! STRUCTURE OF MAIN PROGRAM:
!   PROGRAM Convert
!
! DESCRIPTION:
!   This program converts files between different formats.  Right now, it only
!   convert density files:
!     (1) Old density file format with each density values output on a
!         separate line, and a single header line containing file length info
!     (2) New density file format (direct formatted) with all information
!         appearing on a the same line.  Header information takes up 4 records.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2008/02/06  File created (Linda Hung)
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

                           !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    fileStatus, &
    ix, iy, iz, iSpin, &
    xLen, yLen, zLen, &
    fileSpin

  INTEGER, PARAMETER :: &
    inputUnit = 5, &
    outputUnit = 6

  CHARACTER(len=80) :: &
    inputFile, &        ! Name of original file
    outputFile          ! Name of file in new format

  CHARACTER(len=26) :: &
    temp                ! Temporary storage of variable being transferred

  CHARACTER(len=14) :: &
    comment             ! Junked part of filetype 2 header

  CHARACTER(len=3) :: &
    inputType, &
    outputType

                           !>> INITIALIZATION <<!

  CALL GETARG(1,inputFile)
  CALL GETARG(2,inputType)
  CALL GETARG(3,outputFile)
  CALL GETARG(4,outputType)

  IF (outputFile=="") THEN
    WRITE(*,*) ''
    WRITE(*,*) 'You need to specify appropriate parameters on the command line.'
    WRITE(*,*) 'Format should be:'
    WRITE(*,*) ''
    WRITE(*,*) '  RhoConvert inputFilename inputFiletype outputFilename &
               &outputFiletype'
    WRITE(*,*) ''
    WRITE(*,*) 'Currently implemented filetypes are "dir", "seq", and "tec" &
               &(direct, sequential, tecplot).'
    WRITE(*,*) 'Program cannot run.  LEAVING.'
    WRITE(*,*) ''
    STOP
  END IF

  ! Open input file as appropriate for type
  SELECT CASE (inputType)
  
  CASE("seq")
    OPEN(unit=inputUnit, access="sequential", action="read", blank="null", &
        delim="none", file=TRIM(inputFile), form="formatted", &
        iostat=fileStatus, pad="no", position="rewind", status="old")

  CASE("dir")
    OPEN(unit=inputUnit, access="direct", action="read", blank="null", &
        delim="none", file=TRIM(inputFile), form="formatted", &
        iostat=fileStatus, pad="no", status="old", recl=27)

  CASE DEFAULT
    WRITE(*,*) 'I do not recognize input type "', inputType, '".'
    WRITE(*,*) 'Known input types are "dir" and "seq".  Leaving.'
    STOP

  END SELECT

  IF (fileStatus/=0) THEN
    WRITE(*,*) 'The file ', TRIM(inputFile), ' had a problem on opening.'
    WRITE(*,*) 'My uneducated guess: this input file does not exist.'
    WRITE(*,*) fileStatus
    STOP
  END IF

  ! Open output file as appropriate for type
  SELECT CASE (outputType)

  CASE ("tec")
    OPEN(unit=outputUnit, access="sequential", action="write", blank="null", &
        delim="none", file=TRIM(outputFile), form="formatted", &
        iostat=fileStatus, pad="no", status="new")

  CASE ("seq")
    OPEN(unit=outputUnit, access="sequential", action="write", blank="null", &
        delim="none", file=TRIM(outputFile), form="formatted", &
        iostat=fileStatus, pad="no", status="new")

  CASE ("dir")
    OPEN(unit=outputUnit, access="direct", action="write", blank="null", &
        delim="none", file=TRIM(outputFile), form="formatted", &
        iostat=fileStatus, pad="no", status="new", recl=27)

  CASE DEFAULT
    WRITE(*,*) 'I do not recognize output type "', outputType, '".'
    WRITE(*,*) 'Known output types are "dir", "seq", and "tec".  Leaving.'
    STOP

  END SELECT

  IF (fileStatus/=0) THEN
    WRITE(*,*) 'The file ', TRIM(outputFile), ' had a problem opening.'
    WRITE(*,*) 'My uneducated guess: a file with this name already exists.'
    STOP
  END IF

  ! Format descriptors
  10 FORMAT(A14 I12.1)
  11 FORMAT(A26)
  12 FORMAT(3(I12.1) A3 A26)

                           !>> FUNCTION BODY <<!
  ! Read Header
  SELECT CASE (inputType)

  CASE("seq")
    READ(inputUnit, *, iostat=fileStatus) xLen, yLen, zLen, fileSpin
    IF (fileStatus/=0) THEN
      WRITE(*,*) 'File ', TRIM(inputFile), ' is not in sequential format.'
      WRITE(*,*) 'Leaving.'
    STOP
  END IF

  CASE("dir")
    READ (inputUnit, 10, REC=1, iostat=fileStatus) comment, xLen
    IF (fileStatus/=0) THEN
      WRITE(*,*) 'Input file ',TRIM(inputFile),' is not in direct format.  &
                 &Leaving.'
      STOP
    END IF
    READ (inputUnit, 10, REC=2, iostat=fileStatus) comment, yLen
    IF (fileStatus/=0) THEN
      WRITE(*,*) 'Input file ',TRIM(inputFile),' is not in direct format.  &
                 &Leaving.'
      STOP
    END IF
    READ (inputUnit, 10, REC=3, iostat=fileStatus) comment, zLen
    IF (fileStatus/=0) THEN
      WRITE(*,*) 'Input file ',TRIM(inputFile),' is not in direct format.  &
                 &Leaving.'
      STOP
    END IF
    READ (inputUnit, 10, REC=4, iostat=fileStatus) comment, fileSpin
    IF (fileStatus/=0) THEN
      WRITE(*,*) 'Input file ',TRIM(inputFile),' is not in direct format.  &
                 &Leaving.'
      STOP
    END IF

  END SELECT

  ! Write header
  SELECT CASE (outputType)
  CASE("tec")
    IF (fileSpin/=1) STOP "ERROR: CARTESIAN DOESN'T WORK FOR SPIN>1!"
    WRITE(outputUnit,*) 'variables = "X", "Y", "Z", "RHO"'
    WRITE(outputUnit,*) "zone i=", xLen, ", j=", yLen, ", k=", zLen, &
                                 ", DATAPACKING=POINT"
  CASE("seq")
    WRITE(outputUnit, '(4I7.1)') xLen, yLen, zLen, fileSpin
  CASE("dir")
    WRITE(outputUnit, 10, REC=1) "x-dimension:", xLen
    WRITE(outputUnit, 10, REC=2) "y-dimension:", yLen
    WRITE(outputUnit, 10, REC=3) "z-dimension:", zLen
    WRITE(outputUnit, 10, REC=4) "# of spins: ", fileSpin
  END SELECT

  ! Transcribe values
  DO iSpin = 0, fileSpin-1
    DO iz = 0, zLen-1
      DO iy = 0, yLen-1
        DO ix = 0, xLen-1

          ! Read a value (in character format)
          SELECT CASE (inputType)
          CASE("seq")
            READ(inputUnit,*) temp
          CASE("dir")
            READ(inputUnit, 11, &
              REC = iSpin*xLen*yLen*zLen + iz*xLen*yLen + iy*xLen + ix + 5)&
              temp
          END SELECT

          ! Write a value (in character format)
          SELECT CASE (outputType)
          CASE("tec")
            WRITE(outputUnit, 12) ix, iy, iz, "   ", temp
          CASE("seq")
            WRITE(outputUnit, 11) temp
          CASE("dir")
            WRITE(outputUnit, 11, &
              REC = iSpin*xLen*yLen*zLen + iz*xLen*yLen + iy*xLen + ix + 5)&
              temp
          END SELECT

        END DO
      END DO
    END DO
  END DO

END PROGRAM Convert

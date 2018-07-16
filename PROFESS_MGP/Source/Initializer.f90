MODULE Initializer
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Initializer
!     |_SUBROUTINE InitPROFESS 
!     |_SUBROUTINE SetupAllModules
!     |_SUBROUTINE CleanInitializeInputs
!     |_SUBROUTINE CleanCalculator 
!
! DESCRIPTION:
!   Reads the input file (specifically the .inpt file) from the outside world,
!   and sets all the options, arrays, etc. that are specified to make this 
!   program run!  Also calls ReadIonFile to read the other input files, and
!   asks other modules to initialize themselves.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   All the modules should have a SETUP and a CLEAN funct. to initialize their
!   data, instead of how we USE all of them here, which is yucky.
!
!   Read in and initialization should be seperated.
!
! REFERENCES: 
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE ReadInputFile 

  IMPLICIT NONE

CONTAINS


SUBROUTINE InitPROFESS(systemName)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This procedure runs all the procedures necessary in this module to 
!   initialize the inputs.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!------------------------------------------------------------------------------

  USE SYS, ONLY: rhoR
  USE CellInfo, ONLY: numSpin
  USE Output, ONLY: printDensityClearly
  USE Output, ONLY: transUnit
  USE AtomicDensity, ONLY: GenerateAtomicDensity 

  USE KEDF_WTkernel, ONLY : keKernel ! KEDF kernel in G space
  USE PlaneWave, ONLY : qTable       ! Norm of q
  USE MPI_Functions

  USE KEDF_DenDec, ONLY: do_den_dec
  USE KEDF_DenDec, ONLY: SubstrateCoreDensity
  USE KEDF_WTkernel, ONLY: FillWT_RestartFrom_GSpace
  USE KEDF_WGCkernel, ONLY: FillWGC_RestartFrom_ReciprocalSpace

  USE ReadIonFile, ONLY: ReadGeometry 
  USE ReadIonFile, ONLY: ReadPseudo
  USE ReadIonFile, ONLY: ReadDensity

  IMPLICIT NONE


                    !>> EXTERNAL VARIABLES <<!

  CHARACTER(LEN=*), INTENT(IN) :: systemName
  ! The base name of the system
  !
  
#ifdef __USE_PARALLEL
  !-----------------------------
  ! FOR PARALLEL DEBUGGING
  !-----------------------------
  CHARACTER(LEN=500) :: tempName
#endif


                    !>> INTERNAL VARIABLES <<!
  INTEGER :: fileStatus         ! to see if the file opening worked

                      !>> INITIALIZATION <<!
                      !>> FUNCTION BODY <<!

  ! From the system name, infer all the names of all the input files
  ! (except for the pseudopotential files, which must be read in).
  inputFile = TRIM(systemName)    // defaultInput
  geometryFile = TRIM(systemName) // defaultGeometry
  densityFile = TRIM(systemName)  // defaultDensity
  outputFile = TRIM(systemName)   // defaultOutput
  errorFile = TRIM(systemName)    // defaultError

  ! First read all the options.  This also might overwrite the error and 
  ! output file names.

  CALL StartClock('InitPROFESS')

  CALL ReadOptions

#ifdef __USE_PARALLEL

  ! Two possibilities: if output all the files, 
  ! or just output the rankGlobal==0 file
  IF (output_all_files==1 .OR. & 
    (output_all_files==0 .AND. rankGlobal==0) ) THEN 

    !-----------------------------
    ! FOR PARALLEL DEBUGGING
    !-----------------------------
    IF(output_all_files==1) THEN
      WRITE(tempName,'(I8)') rankGlobal+1
      WRITE(outputFile, '(a,a)') TRIM(outputFile), TRIM(ADJUSTL(tempName))
      WRITE(errorFile, '(a,a)') TRIM(errorFile), TRIM(ADJUSTL(tempName))
      outputFile=TRIM(ADJUSTL(outputFile))
      errorFile=TRIM(ADJUSTL(errorFile))
    ENDIF

#endif
 

  ! Open the ERROR file
  OPEN(UNIT=errorUnit, ACCESS="sequential", ACTION="write", BLANK="null", &
       DELIM="none", FILE=TRIM(errorFile), FORM="formatted",& 
       IOSTAT=fileStatus, PAD="no", POSITION="rewind", STATUS="replace")

  IF (fileStatus/=0) THEN
    WRITE(*,*) "****", TRIM(errorFile), "*****", errorFile, "*****"
    WRITE(*,*)'The file ',TRIM(errorFile),' had a problem on opening.'
    WRITE(*,*)'And I cannot even begin to imagine what it might be.'
    STOP
  END IF

  ! Open the OUTPUT file
  OPEN(UNIT=outputUnit, ACCESS="sequential", ACTION="write", BLANK="null", &
       DELIM="none", FILE=TRIM(outputFile), FORM="formatted",& 
       IOSTAT=fileStatus, PAD="no", POSITION="rewind", STATUS="replace")

  IF (fileStatus/=0) THEN
    WRITE(errorUnit,*)'The file ',TRIM(outputFile),' had a problem opening.'
    WRITE(errorUnit,*)'And I cannot even begin to imagine what it might be.'
    STOP
  END IF

  WRITE(outputUnit,*) "******************************************************************************" 
  WRITE(outputUnit,*) "                            WELCOME TO PROFESS                   "
  WRITE(outputUnit,*) "           (PRinceton Orbital-Free Electronic Structure Software) "
  WRITE(outputUnit,*) "******************************************************************************" 

  CALL CheckOptions

#ifdef __USE_PARALLEL
  END IF
#endif

  ! Now call the routine to set default values of all the modules we're
  ! using, and call the subroutines to read the input files, and set up
  ! the Q-table and the FFT algorithms.  We should now have everything in
  ! place to run the actual program.
  CALL ReadGeometry(geometryFile, defaultPseudo)
  CALL ReadPseudo                     ! The number of electron in the cell
                                      ! is calculated here!
  ! This fully initializes the modules in the system to do our calculation.
  ! (includes stuff like setting up the FFT, etc)

  CALL SetupAllModules

  ! In case we need to use something like an FFT procedure for the density,
  ! we don't read in the density until here.
  ! True => uniform starting density.
  IF( fileDensityFlag .AND. atomicDensityFlag ) THEN
    WRITE(message,*) "Be clear the input charge density is from file or atomic density."
    CALL Error(6, message)
  ENDIF

  IF( fileDensityFlag ) THEN
     CALL ReadDensity(densityFile, xRepeat, yRepeat, zRepeat)
     IF(do_den_dec==1) THEN
       CALL SubstrateCoreDensity(rhoR)
      ENDIF
     ! mohan add, read in atomic density
  ENDIF

  IF( atomicDensityFlag ) THEN
    CALL GenerateAtomicDensity(rhoR,SIZE(rhoR,1),SIZE(rhoR,2),SIZE(rhoR,3),SIZE(rhoR,4))

    ! test
    !DO i=1, SIZE(rhoR,1)
    !  DO j=1, SIZE(rhoR,2)
    !    DO k=1, SIZE(rhoR,3)
    !      WRITE(outputUnit,*) rhoR(i,j,k,1)
    !    END DO
    !  END DO
    !END DO
  ENDIF

  !CALL PrintDensityClearly(rhoR,999)

  IF (rankGlobal == 0) THEN
    ! Open the OUTPUTTRANSITIONFILE if necessary
    IF(outputTransitionState)  THEN
      transFile = TRIM(systemName) // defaultTrans
      OPEN(unit=transUnit, access="sequential", action="write", &
        blank="null", delim="none", file=TRIM(transFile), &
        form="formatted", iostat=fileStatus, pad="no", position="rewind",  &
        status="replace")
      IF (fileStatus/=0) THEN
        WRITE(*,*)'The file ',TRIM(transFile),' had a problem on opening.'
        WRITE(*,*)'And I cannot even begin to imagine what it might be.'
        STOP
      END IF 
    END IF
  END IF

  CALL StopClock('InitPROFESS')

END SUBROUTINE InitPROFESS 


SUBROUTINE SetupAllModules
!----------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine sets up various modules with the parameters that are 
!   available to this module.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Instead of using Calculator module, we should set up the Optimizer module.
!----------------------------------------------------------------------------
! REVISION LOG:
!
!--------------------------------------------------------------------------

  USE SetupKEDF, ONLY : SetupFunctional
  USE KEDF_DenDec, ONLY: AllocateCoreDensity 
  USE RefreshCell, ONLY: RefreshCellSetup
  USE FOURIER, ONLY : PlanFFT ! Procedure to set up various FFT routines
  USE SYS, ONLY: SetupSystem  ! Subroutine to set up the system

  USE CellInfo, ONLY : SetupCellFFTdims

  IMPLICIT NONE

                 !>> FUNCTION BODY <<!

  ! Set up some of the modules with the variables available here
  CALL SetupSystem(energyCutoff, kinetic)

  CALL SetupCellFFTdims()

  CALL SetupFunctional()

  ! This is only necessary when density decomposition 
  ! method is used.
  ! please note we need to allocate space for potDD
  ! in this subroutine, and this requires that the
  ! plane wave related parameters must be initialized
  ! first, which is done in SetupFunctional
  CALL AllocateCoreDensity() ! mohan add 2013-07-30

  ! First three input arguments are dimesnions of local reciprocal space array
  CALL RefreshCellSetup()

  RETURN

END SUBROUTINE SetupAllModules


SUBROUTINE CleanInitializeInputs
!----------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine cleans up some of the stuff initialized in this module
!   in the beginning.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Instead of using Calculator module, we should set up the Optimizer 
!   module.
!--------------------------------------------------------------------------
! REVISION LOG:
!
!--------------------------------------------------------------------------
  USE SYS, ONLY : CleanSystem            ! Procedure to set up various FFT routines
  USE Output, ONLY : CleanOutput         ! clean the allocated vars in output module

  IMPLICIT NONE

  CALL CleanSystem
  CALL CleanCalculator
  CALL CleanOutput

  RETURN

END SUBROUTINE CleanInitializeInputs


SUBROUTINE CleanCalculator
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine just does some cleanup like deallocates memory, etc.  It also
!   calls the cleanup routines for all the modules that it uses to calculate
!   the individual components.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2014-06-04 Mohan 
!------------------------------------------------------------------------------

  USE PlaneWave, ONLY: qTable
  USE PlaneWave, ONLY: qVectors
  USE PlaneWave, ONLY: qMask
  USE IntKernelODE, ONLY: IntKernelODEClean => clean 
  ! clean the all the temperory vars
  USE KEDF_WTkernel, ONLY: keKernel
  USE KEDF_WGCkernel, ONLY: nls_wpp, nls_w1pp, nls_w2pp
  USE CBSpline, ONLY: bSpline1
  USE CBSpline, ONLY: bSpline2
  USE CBSpline, ONLY: bSpline3
  USE IonElectronSpline, ONLY: ionPotReal

  IMPLICIT NONE

                     !>> EXTERNAL VARIABLES <<!
                     !>> INTERNAL VARIABLES <<! 
                     !>> INITIALIZATION <<!

                       !>> FUNCTION BODY <<!
  IF(ALLOCATED(qTable)) DEALLOCATE(qTable)
  IF(ALLOCATED(qVectors)) DEALLOCATE(qVectors)
  IF(ALLOCATED(qMask)) DEALLOCATE(qMask)
  IF(ALLOCATED(bSpline1)) DEALLOCATE(bSpline1)
  IF(ALLOCATED(bSpline2)) DEALLOCATE(bSpline2)
  IF(ALLOCATED(bSpline3)) DEALLOCATE(bSpline3)
  IF(ALLOCATED(keKernel)) DEALLOCATE(keKernel)
  IF(ALLOCATED(nls_wpp)) DEALLOCATE(nls_wpp)
  IF(ALLOCATED(nls_w1pp)) DEALLOCATE(nls_w1pp)
  IF(ALLOCATED(nls_w2pp)) DEALLOCATE(nls_w2pp)
  IF(ALLOCATED(ionPotReal)) DEALLOCATE(ionPotReal)

  CALL IntKernelODEClean

  RETURN

END SUBROUTINE CleanCalculator


END MODULE Initializer 


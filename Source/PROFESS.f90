PROGRAM OFDFT
!------------------------------------------------------------------------------
! STRUCTURE OF MAIN PROGRAM:
!   |_ PROGRAM OFDFT
!
! DESCRIPTION:
!   This is the MAIN program that wraps everything else into its protective
!   shell.  It simply is the skeleton of the program.  It runs OFDFT!
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
!------------------------------------------------------------------------------
! REFERENCES: 
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY : DP                    ! Double precision
  USE CONSTANTS, ONLY : systemNameLen         ! Max # of chars allowed for the system name

  USE MPI_Functions, ONLY: rankGlobal, sizeGlobal
  USE MPI_Functions, ONLY: InitializeMPI, QuitMPI
  USE MPI_Functions, ONLY: BcastCharacter

  USE OutputFiles, ONLY : outputUnit
  USE Output, ONLY : outputRank

  USE REPORT, ONLY : FinalizeOutput           ! Subroutine that prints information on exit
  USE REPORT, ONLY : ReportHeader             ! Subroutine to print the header files

  !-------------------------------------------
  ! STEP 1: initialize PROFESS
  !-------------------------------------------
  USE Initializer, ONLY : InitPROFESS           ! Routine that initalizes all our input files, 
  USE Initializer, ONLY : CleanInitializeInputs ! Cleanup routines to be run on exit

  !-------------------------------------------
  ! STEP 2: optimize cell/ions/atoms 
  !-------------------------------------------
  USE Optimizer, ONLY : OptimizeCell          ! The function to run to optimize the cell
  USE Optimizer, ONLY : mdType                ! default = -1, no MD
  USE Optimizer, ONLY : DoMolecularDynamics   ! Do MolecularDynamics
    
  IMPLICIT NONE

                        !>> INTERNAL VARIABLES <<!

  CHARACTER(LEN=systemNameLen) :: systemName  ! The raw argument on the command line

                         !>> INITIALIZATION <<!

  ! Start the clock 
  CALL InitClocks() 
  CALL StartClock('PROFESS')

  ! Initialize the rankGlobal and sizeGlobal
  CALL InitializeMPI()
  outputRank=rankGlobal

                          !>> FUNCTION BODY <<!

  ! Read the argument on the command line
  ! Non-root nodes have problems reading the line, so root reads and broadcasts
  IF (rankGlobal==0) THEN
    PRINT *, "========================================================"
    PRINT *, "                   WELCOME TO PROFESS                   "
    PRINT *, " (PRinceton Orbital-Free Electronic Structure Software) "
    PRINT *, "========================================================"
    PRINT *, " "
    CALL GETARG(1,systemName) 
    IF (systemName=="") THEN
      WRITE (*,*) 'You need to specify the file on the command line. Leaving.'
      STOP
    END IF
  END IF


  CALL BcastCharacter( systemName, systemNameLen )

  IF (rankGlobal==0) THEN
    WRITE(*,*) ' >> Calling Initialize() ...'
  ENDIF

  !---------------------------
  ! 1) READ THE PARAMETERS
  !---------------------------
  CALL InitPROFESS(systemName)
  CALL ReportHeader(systemName, sizeGlobal)

  !---------------------------
  ! 2) MAIN ROUTINES 
  !---------------------------
  IF (mdType == -1) THEN
    IF (rankGlobal==0) THEN 
      PRINT *, ' >> Calling OptimizeCell() ...'
    ENDIF
    CALL OptimizeCell    
  ELSE
    IF (rankGlobal==0) PRINT *, ' >> Calling DoMolecularDynamics() ...'
    CALL DoMolecularDynamics
  ENDIF

  !---------------------------
  ! 3) FINISH ALL 
  !---------------------------
  CALL FinalizeOutput
  CALL CleanInitializeInputs 
  CALL QuitMPI()

END PROGRAM OFDFT

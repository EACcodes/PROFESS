MODULE Report
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Report
!     |_SUBROUTINE ReportHeader
!     |_SUBROUTINE ReportAnswers
!     |_SUBROUTINE MinimizerReportHeader
!     |_SUBROUTINE MinimizerReportSteps
!     |_SUBROUTINE MinimizerReportFooter
!     |_SUBROUTINE GeometryMinimizerReportHeader
!     |_SUBROUTINE GeometryMinimizerReportSteps
!     |_SUBROUTINE GeometryMinimizerReportFooter
!     |_SUBROUTINE FinalizeOutput
!
! DESCRIPTION:
!   This is a module that contains all the pretty printing routines that print
!   to any and all of the OFDFT output files.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/12  File created (GSH)
!------------------------------------------------------------------------------
                              !<< GLOBAL >>!


  USE OutputFiles, ONLY : outputUnit
  USE MPI_Functions, ONLY : rankGlobal
  USE OUTPUT
  USE CellInfo, ONLY: m3G

  REAL(KIND=DP), SAVE, PRIVATE :: timeCumul
  ! Cumulative time spent in this minimization.
  !

CONTAINS


SUBROUTINE ReportHeader(systemName, numProc)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out the header information in the output file when
!   the run first begins.
!  
! CONDITIONS AND ASSUMPTIONS: This should only be called by the head node 0, 
!   since that's the only node with knowledge of the output files.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE PlaneWave, ONLY: energyCutoff

  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  CHARACTER(LEN=*), INTENT(IN) :: systemName  
  ! The raw argument on the command line     
  !
  INTEGER, INTENT(IN) :: numProc              
  ! Total number of processes
  !

                       !>> INTERNAL VARIABLES <<!
  
  CHARACTER(LEN=10) :: date
  ! Stores the date
  !
  CHARACTER(LEN=10) :: time
  ! Stores the time
  !

                       !>> INITIALIZATION <<!

  outputSystemName = systemName

  ! Format descriptors
  11 FORMAT(A)
  12 FORMAT("*", A75, " *")
  13 FORMAT("*", A30, I6, 40X, "*")
  14 FORMAT("*", A30, 3F12.5, 10X, "*")
  15 FORMAT("*", A39, F14.3, " eV", 20X, "*")
  16 FORMAT("*", A39, I10, 27X, "*")
  17 FORMAT("*", A53, F21.15, 2X, "*")
  18 FORMAT("*", A72, I2, 2X, "*")


                       !>> FUNCTION BODY <<!

  headnode: IF (outputRank==0) THEN
  CALL DATE_AND_TIME(date, time)  

  ! Print the header
  WRITE(*,*) " "
  WRITE(*, 11) REPEAT("*", lineLength)

  WRITE(*, 11) REPEAT("*", (lineLength - 44)/2) // & 
                        "   ORBITAL-FREE DENSITY FUNCTIONAL THEORY   " // & 
                        REPEAT("*", (lineLength - 44)/2)
  WRITE(*, 11) REPEAT("*", lineLength)
  WRITE(*, 12) "Run started on: " // date(5:6) // "/" // &
                        date(7:8) // "/" // date(1:4) // " at " // & 
                        time(1:2)// ":" // time(3:4) // ":" // time(5:6) // " "
  WRITE(*, 12) " System Name: " // systemName // &
                        REPEAT(" ", MAX(lineLength - LEN(systemName),0))
  WRITE(*, 12) " "
  WRITE(*, 12) REPEAT(" ", (linelength - 24)/2) // &
                        "<<< CELL INFORMATION >>>" // &
                        REPEAT(" ", (lineLength - 24)/2)
  WRITE(*, 12) " "
  WRITE(*, 13) "Total Number of Ions: ", SIZE(cell%ionTable)
  WRITE(*, 13) "Number of Types of Ions: ", &
                        SIZE(cell%elementTable)-1

  WRITE(*, 14) "1st Cell Lattice Vector (A):", cell%cellReal(:,1)*bohr
  WRITE(*, 14) "2nd Cell Lattice Vector (A):", cell%cellReal(:,2)*bohr
  WRITE(*, 14) "3rd Cell Lattice Vector (A):", cell%cellReal(:,3)*bohr
  WRITE(*, 12) " "
  WRITE(*, 12) "Boundary Conditions for 1st Lattice Dir.: Periodic" 
  WRITE(*, 12) "Boundary Conditions for 2nd Lattice Dir.: Periodic" 
  WRITE(*, 12) "Boundary Conditions for 3rd Lattice Dir.: Periodic" 
  WRITE(*, 12) " "

  WRITE(*, 12) REPEAT(" ", (linelength - 30)/2) // &
                        "<<< SIMULATION INFORMATION >>>" // &
                        REPEAT(" ", (lineLength - 30)/2)
  WRITE(*, 12) " "
  IF (energyCutoff > 0._DP) THEN
    WRITE(*, 15) "Kinetic Energy Cutoff: ", energyCutoff * hartreeToeV
  END IF
  WRITE(*, 16) "# Spins: ", SIZE(rhoR, 4)
  WRITE(*, 16) "# Electrons in system: ", NINT(SUM(cell%elementTable%chargeTot))

  WRITE(*, 16) "# Gridpoints Along 1st Lattice Dir.: ", SIZE(rhoR, 1)
  WRITE(*, 16) "# Gridpoints Along 2nd Lattice Dir.: ", SIZE(rhoR, 2)

#ifdef __USE_PARALLEL
  WRITE(*, 16) "# Gridpoints Along 3rd Lattice Dir.: ", m3G 
#else
  WRITE(*, 16) "# Gridpoints Along 3rd Lattice Dir.: ", SIZE(rhoR,3)
#endif

  WRITE(*, 16) "# Processors: ", numProc
  WRITE(*, 12) " "

  SELECT CASE (outputRhoMethod)
    CASE(0)
      WRITE(*,12) "Electronic Density Minimization Algorithm: None"
    CASE(1)
      WRITE(*,12) "Electronic Density Minimization Algorithm: NTN (default)"
    CASE(2)
      WRITE(*,12) "Electronic Density Minimization Algorithm: NCG"
    CASE(3)
      WRITE(*,12) "Electronic Density Minimization Algorithm: NBF"
    CASE(4)
      WRITE(*,12) "Electronic Density Minimization Algorithm: STN"
    CASE(5)
      WRITE(*,12) "Electronic Density Minimization Algorithm: SCG"
    CASE(6)
      WRITE(*,12) "Electronic Density Minimization Algorithm: Hybrid SCG + STN"
    CASE(7)
      WRITE(*,12) "Electronic Density Minimization Algorithm: LOG"
  END SELECT

  SELECT CASE (outputIonMethod)
    CASE(0)
      WRITE(*,12) "Ion Geometry Optimization Algorithm: None"
    CASE(2)
      WRITE(*,12) "Ion Geometry Optimization Algorithm: Quickmin Minimization"
    CASE(3)
      WRITE(*,12) "Ion Geometry Optimization Algorithm: Conjugate Gradient Minimization"
    CASE(4)
      WRITE(*,12) "Ion Geometry Optimization Algorithm: Conjugate Gradient Minimization 2"
    CASE(5)
      WRITE(*,12) "Ion Geometry Optimization Algorithm: BFGS Minimization"
    CASE(6)
      WRITE(*,12) "Molecular Dynamics Algorithm: Nose-Hoover Thermostat (NVT)"
    CASE(7)
      WRITE(*,12) "Molecular Dynamics Algorithm: Parrinello-Rahman (NPT)"
    CASE(8)
      WRITE(*,12) "Molecular Dynamics Algorithm: NVE"
    CASE DEFAULT
      WRITE(*,12) "Ion Geometry Optimization Algorithm: Unknown"
  END SELECT

  WRITE(*, 11) REPEAT("*", lineLength)

  END IF headnode

  RETURN

END SUBROUTINE ReportHeader



SUBROUTINE ReportAnswers(energy, title)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine takes an array of energies and prints it out in a pretty
!   format.
!
! CONDITIONS AND ASSUMPTIONS: Prints for head node (final energy), if
!   outputMinimizeDensity>=1, outputMinimizeGeometry>=5
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/10/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: kinetic

  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!
  ! Table of energies. From left to right: total energy, kinetic, external, 
  ! coulombic, exchange-correlation, ion-ion, Thomas-Fermi, von Weiszacker and 
  ! third term Wang-Teter, WGC, ...)
  REAL(kind=DP), DIMENSION(:), INTENT(IN) :: energy          
  !
  CHARACTER(len=*) , INTENT(IN) :: title
  ! The title of this output, max 60 characters
  !

                    !>> INTERNAL VARIABLES <<!

  INTEGER :: titleLength
  ! The length of the title, in the # of characters
  !
  INTEGER :: titleOffset     
  ! 1 if there is the titleLen is an odd number, 0 if it's 
  !
  REAL(kind=DP), DIMENSION(9) :: eVenergy = 0.0_DP
  ! final energies expressed in eV   
  !


                       !>> INITIALIZATION <<!
  titleLength = LEN(title)
  titleOffset = MOD(titleLength, 2)
  eVenergy = energy * hartreeToeV

  ! Format descriptors
  10 FORMAT("# ", A24, " =", 1X, ES20.12, " eV", 12X, F10.3, " s", 1X, "#")
  11 FORMAT(A)
  12 FORMAT("#", A76, "#")
  13 FORMAT("# ", A24, " =", 1X, ES20.12, " eV", 25X, "#")
                      
                       !>> FUNCTION BODY <<!
  ! Print the header
  
  WRITE(outputUnit, 11) REPEAT("#", lineLength)

  WRITE(outputUnit,11) &
    "# " // REPEAT(" ", (lineLength-3-titleLength)/2) // title // &
    REPEAT(" ", (lineLength-3-titleLength-titleOffset)/2-7) // "time    #"

  WRITE(outputUnit, 11) REPEAT("#", lineLength)

  ! Print out kinetic energies
  IF (kinetic /= 2) &
    WRITE(outputUnit,10) "Thomas-Fermi Kin. Energy",  eVenergy(7), energyTime(7)
  IF (kinetic /= 1) &
    WRITE(outputUnit,10) "Von-Weizsacker Kin. Energy",  eVenergy(8), energyTime(8)
  IF (kinetic == 4) &
    WRITE(outputUnit,10) "Wang-Teter Kinetic Energy",  eVenergy(9), energyTime(9)
  IF (kinetic == 5) &
    WRITE(outputUnit,10) "WGC Kinetic Energy",  eVenergy(9), energyTime(9)
  IF (kinetic == 7) & 
    WRITE(outputUnit,10) "LQ Kinetic Energy",  eVenergy(9), energyTime(9) 
  IF (kinetic == 8) & 
    WRITE(outputUnit,10) "HQ Kinetic Energy",  eVenergy(9), energyTime(9) 
  WRITE(outputUnit, 12) " "

  ! Print out Exchange Correlation Energy
  IF (exchangeCorrelationFunct == 1) &
    WRITE(outputUnit,10) "LDA Exch-Corr Energy", eVenergy(5), &
                         energyTime(5)

  IF(exchangeCorrelationFunct == 2) &
    WRITE(outputUnit,10) "GGA(PBE) Exch-Cor Energy", eVenergy(5), &
                         energyTime(5)

  ! Print out Coulombic Energy
  WRITE(outputUnit, 10) "Coulombic Energy", eVenergy(4), energyTime(4)

  ! Print out the Ion-Electron Energy
  WRITE(outputUnit, 10) "Ion-Electron Energy", eVenergy(3), energyTime(3)

  ! Print out the Ion-Ion Energy
  WRITE(outputUnit, 10) "Ion-Ion Energy", eVenergy(6), energyTime(6)

  WRITE(outputUnit, 12) " "
  ! Print out total energies
  WRITE(outputUnit,13) "TOTAL KINETIC ENERGY",  eVenergy(2)
  WRITE(outputUnit,13) "TOTAL POTENTIAL ENERGY", &
                       eVenergy(5) + eVenergy(4) + eVenergy(3) + eVenergy(6)

  WRITE(outputUnit,13) "TOTAL ENERGY", eVenergy(1)

  ! Print the footer
  WRITE(outputUnit, 11) REPEAT("#", lineLength)

  ! Flush the write buffer
  flush(outputUnit)

END SUBROUTINE ReportAnswers



SUBROUTINE MinimizerReportHeader
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Prints out the header of the density minimization report
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
                      !>> INTERNAL VARIABLES <<!
                      !>> INITIALIZATION <<!
  timeCumul = 0._DP
  rhoWatch = TimerStart()
  jumpToMax = .FALSE.
  restarted = .FALSE.

  ! Format descriptors
  11 FORMAT(A)

                       !>> FUNCTION BODY <<!

! outputMinimizeDensity should be always be 0 for non-head nodes
  IF(outputMinimizeDensity >= 1) THEN

    ! Print the header
    WRITE(outputUnit,*) 
    WRITE(outputUnit,11) "|" // REPEAT("-", lineLength-2) // "|"
 
    SELECT CASE (outputRhoMethod)

      ! SQR
      CASE(1)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-55)/2) // &
          "Sqrt Conjugate Gradient Electronic Minimization Report" // &
          REPEAT(" ", (linelength-56)/2) // "|"
      ! NE2
      CASE(2)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-58)/2) // &
          "Sqrt Truncated Newton (v2) Electronic Minimization Report" // &
          REPEAT(" ", (linelength-59)/2) // "|"
      ! NEW
      CASE(3)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-53)/2) // &
          "Sqrt Truncated Newton Electronic Minimization Report" // &
          REPEAT(" ", (linelength-54)/2) // "|"
      ! HYB
      CASE(4)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
          "Hybrid CG/NT Electronic Minimization Report" // &
          REPEAT(" ", (linelength-45)/2) // "|"
      ! BFGS
      CASE(5)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
          "BFGS         Electronic Minimization Report" // &
          REPEAT(" ", (linelength-45)/2) // "|"
      ! NCG
      CASE(6)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
          "NCG          Electronic Minimization Report" // &
          REPEAT(" ", (linelength-45)/2) // "|"

      CASE DEFAULT
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-46)/2) // &
          "Unknown Method Electronic Minimization Report" // &
          REPEAT(" ", (linelength-47)/2) // "|"

    END SELECT

    WRITE(outputUnit,11) "+" // REPEAT("-", lineLength-2) // "+"
    WRITE(outputUnit,11) "| Step         Energy        Pot. Norm" // &
      " #LCG /#Line  line     #FFT   Time (s) |"
    WRITE(outputUnit,11) "| number        (eV)         (in code)" // &
      " Steps/Steps  step    (cumul) (cumul)  |"
    WRITE(outputUnit,11) "+" // REPEAT("-", lineLength-2) // "+"

  END IF

END SUBROUTINE MinimizerReportHeader


SUBROUTINE MinimizerReportSteps(step, energy, potentialNorm, numEnergyLine, &
                                numEnergyBrack, restart, success, duration)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Prints out the individual step data for the density minimization report.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/23/2004  Subroutine Created (GSH)
!------------------------------------------------------------------------------

  USE FOURIER, ONLY : iCountFFT      
  ! Cumulative FFT counter for tracking.
  !

  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy
  ! Table of energies. From left to right: total energy, kinetic, external, 
  ! coulombic, exchange-correlation, ion-ion, Thomas-Fermi, 
  ! von Weiszacker and third term Wang-Teter, WGC, ...)

  INTEGER, INTENT(IN) :: &
    step, &          ! The step number of the minimizer
    numEnergyLine, & ! The number of times the energy was evaluated in the
                     ! total line minimization part
    numEnergyBrack, &! The numberof times the energy was evaluated for 
                     ! the bracketing step.
    success          ! 1 if the line minimizer from the last step 
                     ! successfully found a minima.
                     ! 2 if the line minimizer went to max timestep
                     ! 3 if it it discovered increasing energy at small step

  LOGICAL, INTENT(IN) :: restart
  ! .true. if the last step restarted the conjugate
  ! gradient with the steepest descent vector.

  REAL(KIND=DP), INTENT(IN) :: duration
  ! Total time spent in the minimizer subroutine.
  !
  REAL(KIND=DP), INTENT(IN) :: potentialNorm    
  ! The norm of the potential
  !

                    !>> INTERNAL VARIABLES <<!

  CHARACTER(LEN=2) :: note

                       !>> INITIALIZATION <<!

  ! Format descriptors
  11 FORMAT(A)
  16 FORMAT("|", A2, I4, 1X, ES20.12, 1X, ES11.5, 1X, I3, "/", I3, 2X, ES9.3, &
            1X, I6, 1X, ES9.3, " |")
  note = "  "

                       !>> FUNCTION BODY <<!

  IF (outputMinimizeDensity >= 2) THEN

    IF(restart) THEN
      note(1:1) = "+"
      restarted = .TRUE.
    END IF

    IF(success == 2) THEN
      jumpToMax = .TRUE.
      note(2:2) = "*"
    END IF

    IF(success == 4) THEN
      jumpToMax = .TRUE.
      note(2:2) = "^"
    END IF

    timeCumul = timeCumul + TimerStop(rhoWatch)
    rhoWatch = TimerStart()
    WRITE(outputUnit, 16) note, step, energy(1) * hartreeToeV, &
      potentialNorm, numEnergyLine, numEnergyBrack, duration, &
      iCountFFT, timeCumul

    ! Flush the write buffer
    flush(outputUnit)

    IF(outputMinimizeDensity >= 3) THEN
      CALL ReportAnswers(energy, "Step energy report")
      WRITE(outputUnit, 11) " "
    END IF

  END IF

  IF (outputIntermediateDensity) THEN
    CALL PrintDensity(rhoR)
  END IF

  RETURN

END SUBROUTINE MinimizerReportSteps


SUBROUTINE MinimizerReportFooter(steps, cause, energy, numEnergyEvals, numPotentialEvals, duration)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine gets called when an electronic energy minimization is done
!   and it prints out status information about what happened.
!  
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/03/2004  Subroutine created (Vincent Ligneres)
!   02/23/2004  Modified to include more information (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy
  ! Table of energies. From left to right: total energy, kinetic, external, 
  ! coulombic, exchange-correlation, ion-ion, Thomas-Fermi, von Weiszacker & 
  ! third term Wang-Teter, WGC, ...)

  INTEGER, INTENT(IN) :: &
    steps, &          ! Number of steps the minimizer ran for.
    cause, &          ! The reason why the minimizer stopped.
    numEnergyEvals, & ! Total number of times the energy has been evaluated
    numPotentialEvals ! Total number of times potential has been evaluated

  REAL(KIND=DP), INTENT(IN) :: duration     
  ! Total time spent in the minimizer subroutine.

                    !>> INTERNAL VARIABLES <<!
                       !>> INITIALIZATION <<!

  ! Format descriptors
  10 FORMAT("| ", A75, "|")
  11 FORMAT(A)
  13 FORMAT("| ", A22, 8X, I9, 36X, "|")
  14 FORMAT("| ", A22, 12X, F9.3, " s", 30X, "|")
  15 FORMAT("| ", A29, 1X, I9, 36X, "|")
  16 FORMAT("| ", A29, 8X, ES8.1, 30X, "|")

                       !>> FUNCTION BODY <<!
  IF (outputMinimizeDensity >= 2 .AND. jumpToMax) THEN
    WRITE(outputUnit, 10) " "
    WRITE(outputUnit, 10) "    * = No minimum found in line minimization ... &
                               &went to max timestep "
  END IF

  IF (outputMinimizeDensity >= 2 .AND. restarted) THEN
    WRITE(outputUnit, 10) " "
    WRITE(outputUnit, 10) "    + = Triggered a restart from Steepest Descent "
  END IF

  IF (outputMinimizeDensity >= 1) THEN
    WRITE(outputUnit, 10) " "  
    SELECT CASE (cause)
      CASE (0)
        WRITE(outputUnit,10) "Stationary point found!  "
      CASE (1)
        WRITE(outputUnit,10) "Maximum number of steps exceeded!  "
      CASE (2)
        WRITE(outputUnit,10) "** Major method breakdown occured! **"
      CASE DEFAULT
        WRITE(outputUnit,10) "I don't know what happened but it can't be good.."
    END SELECT

    WRITE(outputUnit,13) "Total number of steps:", steps
    WRITE(outputUnit,15) "Total number of energy evals:", numEnergyEvals
    WRITE(outputUnit,15) "Total number of potential evals:", numPotentialEvals
    WRITE(outputUnit,16) "Machine precision for energy:", &
                  (NEAREST(energy(1),energy(1)*2._DP)-energy(1))*hartreeToeV

    WRITE(outputUnit,14) "Minimization time:", duration

    WRITE(outputUnit,11) "+" // REPEAT("-", lineLength-2) // "+"

    IF(outputMinimizeDensity >= 2) THEN
      CALL ReportAnswers(energy, "Intermediate energy report")
      WRITE(outputUnit, 11) " "
    END IF

    IF(cause == 2) THEN
      CALL ReportAnswers(energy, "FAILURE energy report")
      WRITE(outputUnit,10) "**Major method breakdown occured in DENSITY MIN!**"
      WRITE(outputUnit,10) "**Major method breakdown occured in DENSITY MIN!**"
      WRITE(outputUnit,10) "**Major method breakdown occured in DENSITY MIN!**"
    END IF

    IF(cause == 3) THEN
      WRITE(outputUnit,10) "** EXITED BECAUSE OF ENERGY CRITERIA! **"
    END IF

  END IF

END SUBROUTINE MinimizerReportFooter



SUBROUTINE GeometryMinimizerReportHeader(extraInfo)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Prints out the header for the geometry minimization report.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/23/2004  Created (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

   CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: extraInfo      
   ! Line of information added to the header

                      !>> INTERNAL VARIABLES <<!
                      !>> INITIALIZATION <<!

  ! Format descriptors
  11 FORMAT(A)

                       !>> FUNCTION BODY <<!
  IF(outputMinimizeGeometry >= 1) THEN

    ! Print the header
    WRITE(outputUnit,*)
    WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"
    
    ! IonMethod 0 (NON): No optimization.
    ! IonMethod 2 (QUI): Quickmin optimization method.
    ! IonMethod 3 (CON): CG method.
    ! IonMethod 4 (CG2): CG method version 2. 
    ! IonMethod 5 (BFG): BFGS method.
    SELECT CASE(outputIonMethod)
      CASE(2)
        WRITE(outputUnit,11) "[" // REPEAT(" ", (lineLength-42)/2) // &
          "Quickmin Ion Geometry Minimization Report" // &
          REPEAT(" ", (linelength-43)/2) // "]"
      CASE(3)
        WRITE(outputUnit,11) "[" // REPEAT(" ", (lineLength-52)/2) // &
          "Conjugate Gradient Ion Geometry Minimization Report" // &
          REPEAT(" ", (linelength-53)/2) // "]"
      CASE(4)
        WRITE(outputUnit,11) "[" // REPEAT(" ", (lineLength-57)/2) // &
          "Conjugate Gradient (v2) Ion Geometry Minimization Report" // &
          REPEAT(" ", (linelength-58)/2) // "]"
      CASE(5)
        WRITE(outputUnit,11) "[" // REPEAT(" ", (lineLength-38)/2) // &
          "BFGS Ion Geometry Minimization Report" // &
          REPEAT(" ", (linelength-39)/2) // "]"
      CASE DEFAULT
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-40)/2) // &
          "Unknown Method Ionic Minimization Report" // &
          REPEAT(" ", (linelength-41)/2) // "|"

    END SELECT

    IF (PRESENT(extraInfo)) WRITE(outputUnit,'(A)') TRIM(extraInfo)

    WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"

    IF(outputMinimizeGeometry >= 2) THEN
      IF (outputIonMethod > 6) THEN
        WRITE(outputUnit,11) "[ Step     MaxForce       MaxStep  " // &
          "    Hamiltonian   Temperature       Time  ]" 
        WRITE(outputUnit,11) "[   #       (eV/A)          (A)    " // &
          "       (eV)            (K)          (s)   ]"
        WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"
      ELSE
        WRITE(outputUnit,11) "[ Step     MaxForce     Time Step  " // &
          "Velocity   Overstep?                Time  ]" 
        WRITE(outputUnit,11) "[   #       (eV/A)                 " // &
          "                                     (s)  ]"
        WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"
      END IF
    END IF

  END IF

END SUBROUTINE GeometryMinimizerReportHeader


SUBROUTINE GeometryMinimizerReportSteps(step, duration, forces, maxForce, &
                                        stepSize, overStep, velocityMag, &
                                        stepDone, positions)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Prints out the per-step information for the geometry minimization report
!
! outputMinimizeGeometry
! 3 : Print total force.
! 4 : plus ion-ion and ion-electron forces.
! 5 : plus energy information
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/23/2004  Subroutine Created (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(IN) :: step 
  ! The step number of the minimizer
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: forces 
  ! forces for this step
  !
  REAL(KIND=DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: positions 
  ! atom positions for this step
  !
  REAL(KIND=DP), INTENT(IN) :: duration 
  ! Total time spent in the minimizer subroutine.
  !
  REAL(KIND=DP), INTENT(IN) :: maxForce 
  ! Maximum value of the forces
  !
  REAL(KIND=DP), INTENT(IN) :: stepSize
  !
  REAL(KIND=DP), INTENT(IN) :: velocityMag
  !
  LOGICAL, INTENT(IN) :: overStep 
  ! means the energy is higher than previous step
  !
  LOGICAL, INTENT(IN) :: stepDone 
  ! means this step decrease the energy
  !
  LOGICAL :: cellRelaxFlag = .FALSE.
  ! this is not a report for cell relaxation
  !

                    !>> INTERNAL VARIABLES <<!
                       !>> INITIALIZATION <<!
  ! Format descriptors
  11 FORMAT(A)
  16 FORMAT("[ ", I4, 3X, ES12.6, 3X, ES9.3, 3X, ES9.3, 3X, L1, &
            18X, F9.3, 1X, "]")

                       !>> FUNCTION BODY <<!

  IF(outputMinimizeGeometry >= 2) THEN
      !----------------------------------------
      ! EXPLANATION FOR EACH PRINTOUT PARAMETER
      !----------------------------------------
      ! IONIC STEP
      ! MAX FORCCE (eV/Angstrom)
      ! STEP SIZE (Angstrom)
      ! VELOCITY
      ! OVE STEP (flag)
      ! TIME (Second)
      WRITE(outputUnit, 16) step, maxForce*hartreeToeV/bohr, &
                            stepSize, velocityMag*hartreeToeV/bohr, &
                            overStep, duration
  END IF


  IF(outputMinimizeGeometry >= 3 .AND. stepDone) THEN
    CALL PrintForces(forces(:,:,1), " TOTAL_FORCE", step) ! Total forces
    IF(outputMinimizeGeometry >= 4) THEN
      CALL PrintForces(forces(:,:,2), " ION-ION FORCE", step)
      CALL PrintForces(forces(:,:,3), " ION-ELECTRON FORCE", step)
    END IF
    WRITE(outputUnit, 11) " "
  ENDIF
  
  IF(outputMinimizeGeometry >= 5) THEN
    CALL ReportAnswers(energy, "Energy Report")
    WRITE(outputUnit, 11) " "
  END IF

  ! Save ion positions / geometry
  IF(PRESENT(positions) .AND. geoOutputFreq>0 .AND. &
     (step==1 .OR. MOD(step,geoOutputFreq)==0)) THEN 
    IF(rankGlobal==0) THEN ! mohan add 2013-07-24
!      WRITE(outputUnit,*) " print the geometry, step is ", step
      CALL PrintGeometry(cellRelaxFlag, step, positions)
    ENDIF
  ENDIF

  ! Flush the write buffer
  flush(outputUnit)

  RETURN

END SUBROUTINE GeometryMinimizerReportSteps


SUBROUTINE GeometryMinimizerReportFooter(duration)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine gets called when a geometry minimization is done
!   and it prints out status information about what happened.
!  
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/03/2004  Subroutine created (Vincent Ligneres)
!   02/23/2004  Modified to include more information (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), INTENT(IN) :: duration
  ! Total time spent in the minimizer subroutine.
  !

                    !>> INTERNAL VARIABLES <<!
                       !>> INITIALIZATION <<!

  ! Format descriptors
  10 FORMAT("[ ", A75, "]")
  11 FORMAT(A)
  14 FORMAT("[ ", A22, 12X, F9.3, " s", 30X, "]")

                       !>> FUNCTION BODY <<!

  IF (outputMinimizeGeometry >= 1) THEN
    WRITE(outputUnit, 10) " "  
    WRITE(outputUnit,14) "Minimization time:", duration
    WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"
  END IF

  RETURN

END SUBROUTINE GeometryMinimizerReportFooter



SUBROUTINE FinalizeOutput
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Finalizes the output at the end of the OFDFT program.
!
! CONDITIONS AND ASSUMPTIONS: Must be called by at least all processors
!   containing non-padded density
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE MPI_Functions
  USE CellInfo, ONLY: cell, m123G
  USE KEDF_DenDec, ONLY: do_den_dec
  USE KEDF_DenDec, ONLY: AddCoreDensity

  IMPLICIT NONE
                      !>> INTERNAL VARIABLES <<!
  
  CHARACTER(LEN=10) :: date          
  ! Stores the date
  !
  CHARACTER(LEN=10) :: time               
  ! Stores the time
  !

                       !>> INITIALIZATION <<!

                        !>> FUNCTION BODY <<!

  ! Print out the final density
  IF (outputFinalDensity) THEN
    ! add the core density into final density for 
    ! density decomposition scheme.
    IF(do_den_dec==1) THEN
      CALL AddCoreDensity(rhoR)
    ENDIF
    !WRITE(outputUnit,*) "Print out final density with charge ", SUM(rhoR)*cell%dV 
    CALL PrintDensity(rhoR)
  END IF

  IF (outputRank==0) THEN
    ! Print our final energy
    CALL ReportAnswers(energy, "FINAL ENERGIES")

    ! Print out the final forces (but has already been done if outputMinGeom>=3)
    IF (outputFinalForces .AND. outputMinimizeGeometry<3) THEN
      CALL PrintForces(forceIon(:,:,1), outputSystemName)
    ENDIF

    ! Print out the final stresses
    IF (outputFinalStress .EQV. .TRUE.) THEN 
      CALL PrintStress(stress)
    ENDIF

    ! Print out the final geometry - this may print a repeat file, if we were
    ! asking for geometry output during the run, but we'll not worry about that
    ! for now

    IF (outputFinalGeometry .EQV. .TRUE.) THEN
      CALL PrintGeometry
    ENDIF

    CALL StopClock('PROFESS')

    ! whole clocks
    CALL PrintClock(' ',outputUnit) 
    ! some chosen clocks
    WRITE(outputUnit,*) "------------------------------- TOTAL ---------------------------------------"
    CALL PrintClockWith('PROFESS',outputUnit)
    CALL PrintClockWith('Initialize',outputUnit)
    CALL PrintClockWith('MD',outputUnit)
    WRITE(outputUnit,*) "-----------------------------------------------------------------------------"

    WRITE(6,*) "--------------------- END OF PROFESS, HAVE A GREAT DAY ! --------------------"
    CALL PrintClockWith('PROFESS',6)
    CALL PrintClockWith('Initialize',6)
    !CALL PrintClockWith('ForwFFT_3D',6)
    !CALL PrintClockWith('BackFFT_3D',6)
    WRITE(6,*) "-----------------------------------------------------------------------------"


    CALL DATE_AND_TIME(date, time)

    WRITE(outputUnit, '(/ A)') "Run completed on: " // date(5:6) // "/" // &
                        date(7:8) // "/" // date(1:4) // " at " // & 
                        time(1:2)// ":" // time(3:4) // ":" // time(5:6) // " "
    WRITE(outputUnit,*) " "

    CLOSE(outputUnit)
    CLOSE(errorUnit)

    IF (outputOptionGeom) CLOSE(outputGeomUnit)

    IF(outputTransitionState) CALL PrintTransitionState

    flush(6)

  END IF

#ifdef __USE_PARALLEL
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
  IF (mpiErr /= MPI_SUCCESS) STOP "***MPI_BARRIER PROBLEM IN FINAL OUTPUT***"
#endif

  RETURN

END SUBROUTINE FinalizeOutput

END MODULE Report


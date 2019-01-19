MODULE ReadInputFile
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE ReadInputFile
!     |_SUBROUTINE ReadOptions
!     |_SUBROUTINE CheckOptions
!
! DESCRIPTION:
!   Reads the input file (specifically the .inpt file) from the outside world,
!   and sets all the options, arrays, etc. that are specified to make this 
!   program run! 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   All the modules should have a SETUP and a CLEAN funct. to initialize their
!   data, instead of how we USE all of them here, which is yucky.
!
!   Read in and initialization should be seperated.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

                         !>> GLOBAL <<!

  USE Constants, ONLY: DP                ! Double precision
  USE Constants, ONLY: pi                ! 3.1415...
  USE Constants, ONLY: auPressure        ! convert pressure from Pa to a.u.
  USE Constants, ONLY: boltzmann         ! Boltzmann constant
  USE Constants, ONLY: fundamentalTime   ! Atomic time unit, in seconds
  USE Constants, ONLY: systemNameLen     ! The number of characters allowed in the system name
  USE Constants, ONLY: bohr              ! The bohr radius
  USE Constants, ONLY: hartreeToeV       ! energy unit from Hartree to eV

  USE OutputFiles
  USE Output, ONLY: WrtOut

  USE ReadIonFile, ONLY: inputUnit  ! The unit used to read input files

!!!!! BELOW THIS LINE ARE PARAMETERS DIRECTLY MODIFIED BY THIS MODULE !!!!!

  USE CellInfo, ONLY : cell ! The cell, contains lattice vectors and ions.
  USE CellInfo, ONLY: kinetic ! The kinetic energy functionals used
  USE CellInfo, ONLY: usingGridPack
  USE CellInfo, ONLY: numSpin ! Number of spins in calculation (1 or 2) 
  USE Sys, ONLY: magmom  ! Magnetic moment
  USE Sys, ONLY: rho0
  USE Sys, ONLY: rhoS
  USE Sys, ONLY: hold0
  USE Sys, ONLY: holdS
  USE Sys, ONLY: bvac
  USE Sys, ONLY: gridSpacing

  USE Sys, ONLY: LumgpExp, LumgpFactor  ! Used by MGP KEDF
  USE PlaneWave, ONLY: energyCutoff ! The kinetic energy cutoff
  USE SetupFFT, ONLY: dimType  ! Specifies how to calculate the grid size  

  USE Ewald, ONLY: errorTolerance ! The error tolerance for the Ewald
  USE Ewald, ONLY: maxRealCutoff ! The maximum real space cutoff for Ewald
  USE Ewald, ONLY: etaIncrement ! The stepsize for eta
  USE Ewald, ONLY: recipCutoffIncrement! The stepsize for the recip. cutoff

  USE IonElectronSpline, ONLY: iiSpline ! Use cardinal b-spline for ion-ion terms
  USE IonElectronSpline, ONLY: ieSpline ! Use cardinal b-spline for ion-electron terms
  USE NearestDistance, ONLY: nnDist ! User defined nearest distance between atoms
  USE NearestDistance, ONLY: checkNNdist_bin ! True: space is binned to treat more
  USE RefreshIons, ONLY: trashPreviousRho

  USE CBSpline, ONLY: splineOrder    

  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu

  USE KEDF_WTkernel, ONLY: alpha, beta
  USE KEDF_WGCkernel, ONLY: firstOrderWGC
  USE KEDF_WGCkernel, ONLY: WGCT
  USE KEDF_WGCkernel, ONLY: gamma
  USE KEDF_WGCkernel, ONLY: alpha5
  USE KEDF_WGCkernel, ONLY: beta5
  USE KEDF_WGC, ONLY : rhoV, rhoStep            ! set the rhoV and rhoStep used for WGV

  USE KEDF_CAT, ONLY: cat_alpha    ! parameters for CAT KEDF
  USE KEDF_CAT, ONLY: cat_beta     ! parameters for CAT KEDF 
  USE KEDF_CAT, ONLY: cat_gamma    ! parameters for CAT KEDF

  USE KEDF_HC10, ONLY: refRatio      ! ratio for interpolating kernels 
  USE KEDF_HC10, ONLY: hc_lambda_val ! the value that hc_lambda array will take
  USE KEDF_HC10, ONLY: cutrho        ! for rho smaller than cutrho, set it to zero

  USE KEDF_DenDec, ONLY: do_den_dec ! Do transition metal.
  USE KEDF_DenDec, ONLY: aTF        ! Coefficient for TF used for transition metal.
  USE KEDF_DenDec, ONLY: bVW        ! Coefficient for VW used for transition metal.

  USE KEDF_WGCD, ONLY: rhoc
  USE KEDF_WGCD, ONLY: shiftm
  USE KEDF_WGCD, ONLY: scfc
  USE KEDF_WGCD, ONLY: mrhoS

  USE KEDF_EvW, ONLY: kloc, aloc, tolk

  USE KEDF_GGA, ONLY: model

  USE XC_LDA, ONLY: exchangeCorrelation ! The exchange correlation functionals used
  USE XC_PBE, ONLY: pbeCutoff
  USE KEDF_GGA, ONLY: GGA_functional
  USE KEDF_GGA, ONLY : CP               ! penalty coefficient to help GGA convergence

  USE Optimizer, ONLY: calOption
  USE Optimizer, ONLY: rhoMethod
  USE Optimizer, ONLY: ionMethod
  USE Optimizer, ONLY: cellRelax
  USE Optimizer, ONLY: mdType

  USE RhoOptimizers, ONLY: maxIter
  USE RhoOptimizers, ONLY: niter_extra
  USE RhoOptimizers, ONLY: nmag !maximal magnetici iteration steps.
  USE RhoOptimizers, ONLY: tolm !tolerence for magnetic moment
  USE RhoOptimizers, ONLY: conv_check ! check energy or potential
  USE RhoOptimizers, ONLY: pot_tol
  USE RhoOptimizers, ONLY: tole
  USE RhoOptimizers, ONLY: usePreconditioner
  USE RhoOptimizers, ONLY: cheatPot
  USE RhoOptimizers, ONLY: fixedmag 

  USE RhoDirCG, ONLY: cg_alg
  USE RhoDirBFGS, ONLY : MBFGS ! parameter for lbfgs method

  USE IonOptimizers, ONLY: forceCutoff
  USE IonOptimizers, ONLY: cg_type
  USE IonOptimizers, ONLY: timeStep
  USE IonOptimizers, ONLY: maxIonStep

  USE CellOptimizers, ONLY: maxCellStep
  USE CellOptimizers, ONLY:tols ! tolerence for stress

  USE MolecularDynamics, ONLY: md_output_path ! output path of md files.
  USE MolecularDynamics, ONLY: nResn
  USE MolecularDynamics, ONLY: nYosh
  USE MolecularDynamics, ONLY: temperature
  USE MolecularDynamics, ONLY: timeTot
  USE MolecularDynamics, ONLY: dt
  USE MolecularDynamics, ONLY: tau_thermo 
  USE MolecularDynamics, ONLY: tau_baro
  USE MolecularDynamics, ONLY: Qmass         ! thermostat's mass
  USE MolecularDynamics, ONLY: dump_md_freq  ! how frequent to dump MD information,
  USE MolecularDynamics, ONLY: msd_startTime 
  USE MolecularDynamics, ONLY: rstMD         ! whether or not to restart MD calculation
  USE MolecularDynamics, ONLY: velRescale    ! whether to rescale the initial velocities
  USE MolecularDynamics, ONLY: fixTemperature
  USE MolecularDynamics, ONLY: doMSDAtom     ! whether print out MSD for each atom

  USE NPT, ONLY : constr_type ! how to constrain box during NPT
  USE NPT, ONLY : extPres     ! external pressure 
  USE NPT, ONLY : bMass       ! barastat's mass
  
! These are duplicates of various parameters in the OUTPUT module. The reason
! we have to have a copy in OUTPUT is that OUTPUT is compiled before 
! everything else and so can't references a lot of the modules that these
! reside in.  Thus, we have a copy specifically for output that is used only
! in printing out the parameters.


  USE OUTPUT, ONLY : &
  outputFinalForces, & ! Print out forces at the end of the calc in a separatefile 
  outputFinalDensity, &    ! Print the final electron density
  outputFinalPotential, &  ! Print the final electronic potential
  outputFinalGeometry, &   ! Print the final coordinates of ions
  outputFinalStress, &     ! Whether to output the final stress
  outputEwald, &           ! Print ewald summation data
  outputKernel, &          ! Print the kernels
  outputIntermediateDensity, & ! Print the intermediate density
  outputMinimizeDensity, & ! How much detail to print for the density min.
  outputMinimizeGeometry, &! How much detail to print for the geometry min.
  outputForcesSeparate, &  ! Whether to seperate forces into a separate file.
  maxMB, &
  exchangeCorrelationFunct, &
  outputRhoMethod,&        ! copy of method for density optimizations
  outputIonMethod, &       ! copy of method for minimizing ions
  outputOptionGeom, &      ! did we do an ion optimization?
  outputTransitionState, & ! Output for transition state
  geoOutputFreq, &         ! dump geometry every geoOutputFreq ion-opt steps
  celOutputFreq            ! dump geometry every celOutputFreq cel-opt steps

  USE MPI_Functions
  USE MathFunctions, ONLY: UpperCase


  IMPLICIT NONE

  CHARACTER(LEN=7), PARAMETER :: defaultPseudo = ".recpot"  ! Default psp suffix.
  CHARACTER(LEN=6), PARAMETER :: defaultTrans = ".trans"    ! Default suffix for transistion state files
  CHARACTER(LEN=5), PARAMETER :: defaultInput = ".inpt"     ! Default suffix for input files
  CHARACTER(LEN=4), PARAMETER :: defaultGeometry = ".ion"   ! Default suffix for the geometry file
  CHARACTER(LEN=4), PARAMETER :: defaultDensity = ".den"    ! " " for the density file
  CHARACTER(LEN=4), PARAMETER :: defaultOutput = ".out"     ! " " for the output file
  CHARACTER(LEN=4), PARAMETER :: defaultError = ".err"      ! " " for the error file

  CHARACTER(LEN=systemNameLen + 6) :: transFile          ! The name of the transition state file

  CHARACTER(LEN=systemNameLen + 5) :: &
  inputFile, &       ! The name of the input file
  geometryFile, &    ! The name of the file with the lattice info.
                     ! Set by SetDefaults, possibly modified by 
                     ! ReadOptions and used by ReadGeometry.
  densityFile, &     ! The initial electron density file name.
  outputFile, &      ! The output file name.
  errorFile          ! The error file name.

  ! Private variables
  LOGICAL :: fileDensityFlag = .FALSE. ! Is there a density file? true => yes, false => no.
  LOGICAL :: fileKernelFlag = .FALSE. !  Is there a kernel file? true => yes, false => no.  
  LOGICAL :: atomicDensityFlag = .FALSE. ! do we use atomic density? true = > yes, false => no.

  INTEGER :: &
  coulombOption = 1, &     ! Copy of the coulombOption in Multigrid
  xRepeat = 1, &           ! Number of times to repeat density file in x-dim
  yRepeat = 1, &           ! Number of times to repeat density file in y-dim
  zRepeat = 1              ! Number of times to repeat density file in z-dim

  INTEGER :: output_all_files=0 ! Add by mohan, 2013-01-20.

CONTAINS


SUBROUTINE ReadOptions
!------------------------------------------------------------------------------
! DESCRIPTION:
! This routine takes the input file, which has the name of the argument on the
! command line plus a '.inpt' suffix, and reads it to assign all the options
! for the current calculation. Read the manual for the exact use and role of
! each option. The values set here override those from SetDefaults.  This 
! subroutine also converts things from SI to AU.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! All right, there are several ways this file can be "pruned down" from its
! monstrously large size :)  One is just to write a simple subroutine
! "GetOption" that does a lot of the reading work for you ...
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!--------------------------------------------------------------------------

  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!

                       !>> INTERNAL VARIABLES <<!
  INTEGER :: tempInt
  ! Temporary integer
  !
  INTEGER :: fileStatus
  ! Control integer. Checks whether files exist.
  !
  CHARACTER(LEN=20) :: keyword, option 
  ! Character strings to read options.
  !
  REAL(kind=DP) :: tempReal 
  ! Actual value of a numerical option.
  !
  
                      !>> INITIALIZATION <<!


                      !>> FUNCTION BODY <<!

  ! First, open the input file with all the options in it
  OPEN(unit=inputUnit, access="sequential", action="read", blank="null", &
       delim="none", file=TRIM(inputFile), form="formatted",& 
       iostat=fileStatus, pad="no", position="rewind",  status="old")
  
  IF (fileStatus/=0) THEN
    WRITE(*,*)'The file ', TRIM(inputFile), ' had a problem on opening. filestatue=',fileStatus
    WRITE(*,*)'My uneducated guess: this input file does not exist.'
    STOP
  END IF

  ! default spin is 1
  numSpin = 1

  ! Now parse through the input file one line at a time.
  DO
    fileStatus = 0

    ! Reading first word or eof.
    READ(inputUnit,*, IOSTAT=fileStatus) keyword 

    ! If this is the end of the file, exit
    IF (fileStatus<0) EXIT 

    ! Put the keyword read in all uppercase to compare.
    CALL Uppercase(keyword)

    BACKSPACE inputUnit

    ! And now we're looking for the keyword.  Ignore all but the first four
    ! letters.
    SELECT CASE (TRIM(keyword(1:4)))

      !-----------------------------------------------------------------------------
      ! 01) Specifying the plane wave energy cutoff, in eV. 
      CASE("ECUT") 
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, energyCutoff
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'Ecut format error. Format is: Ecut [REAL VALUE (eV)]'
          CALL Error(6, message)
        END IF
        IF (energyCutoff<=0) THEN
          WRITE(message,*) 'Specified ecut value must be positive: ', energyCutoff
          CALL Error(6, message)
        END IF 
        WRITE(message,'('' (input) Energy Cutoff                   : '',E10.2, '' eV'')') energyCutoff
        CALL WrtOut(6, message)
        ! Adjust the KE cutoff to appropriate units.
        energyCutoff = energyCutoff / hartreeToeV

      !-----------------------------------------------------------------------------
      ! 02) Alternatively, we can specify a grid density (in inverse angstroms)
      CASE("GDEN") 
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, gridSpacing
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'GDEN format error. Format is: GDEN [REAL VALUE]'
          CALL Error(6, message)
        END IF
        IF (gridSpacing<=0) THEN
          WRITE(message,*) 'Specified GDEN value must be positive: ', gridSpacing
          CALL Error(6, message)
        END IF 
        ! Adjust the grid density to Bohr^-1 
        gridSpacing = gridSpacing * bohr

      !-----------------------------------------------------------------------------
      ! 03) Specifying how the dimension should be calculated
      CASE("DIME") 
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(*,*) 'DIME format error. Format is: DIME [OPTION]'
          STOP
        END IF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:3))) 
        CASE("ALL") ! Both even and odd gridsizes
          dimType = 1
        CASE("EVE") ! Only even gridsizes
          dimType = 2
        CASE("ODD") ! Only odd gridsizes
          dimType = 3
        CASE("TWO") ! Strictly powers of two
          dimType = 4
        CASE DEFAULT
          WRITE(message,*) 'Warning: Encountered unknown DIME argument ', TRIM(option) 
          CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------
      ! 04) Setting the job as a minimization.
      CASE("MINI") 
        mdType=-1 ! no MD for default
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(*,*) 'MINI format error. Format is: MINI [OPTION]'
          STOP
        END IF

        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:3)))
        CASE("RHO", "DEN") ! Minimize only the electron density.
          calOption(1) = .TRUE. ! Minimize rho.
        CASE("ION", "GEO") ! Minimize the ion coordinates in the cell, and rho.
          calOption(1) = .TRUE. ! Minimize rho.
          calOption(2) = .TRUE. ! Minimize the forces.
          calOption(4) = .TRUE. ! Calculate the forces.
        CASE("CEL", "STR") ! Relax the cell dimensions, ion positions, and rho.
          calOption(1:5) = .TRUE. ! Minimize rho, forces, cell Calculate F, stress.
        CASE("NVT") ! Perform molecular dynamics using canonical ensemble (Nose-Hoover)
          CALL WrtOut(6, " (input) Perform Molecular Dynamics      : NVT")
          mdType =1
          ionMethod=6
        CASE("NPT") ! Perform molecular dynamics using isothermal-isobaric
                    ! ensemble (Parinello-Rahman)
          CALL WrtOut(6, " (input) Perform Molecular Dynamics      : NPT")
          mdType =2
          ionMethod=7
        CASE("NVE")
          CALL WrtOut(6, " (input) Perform Molecular Dynamics      : NVE")
          mdType =3
          ionMethod=8
        CASE DEFAULT
          WRITE(message,*) 'Warning: Encountered unknown MINI argument ', TRIM(option)
          CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------
      ! 05) allowed maximal number for each calculations
      CASE("MAXI")
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, option, tempInt 
        IF (fileStatus/=0) THEN
          WRITE(*,*) 'MAXI format error. Format is: MAXI [OPTION] [INT]'
          STOP
        END IF

        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:3)))
        CASE("DEN","RHO") ! max electronic steps
          maxIter=tempInt
        CASE("ION", "GEO") ! max ion steps 
          maxIonStep=tempInt 
        CASE("CEL", "STR") ! max ion steps 
          maxCellStep=tempInt 
        CASE("EXT") ! extra electronic optimization steps
          niter_extra=tempInt
        CASE("MAG") ! max magnetic steps
          nmag=tempInt
        CASE DEFAULT
          WRITE(message,*) 'Warning: Encountered unknown MAXI argument ', TRIM(option)
          CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------
      ! 06) specifying the energy minimization method.
      CASE("METH") 
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(*,*) 'METH format error. Format is: METH [OPTION]'
          STOP
        END IF

        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:3)))
          CASE("NON") ! No minimization for density
            rhoMethod = 0
          !---------------------------------------------------
          CASE("NTN") ! This is sqrt TN method that conserves electron number
                      ! systematically (not scaling, after each line search
                      ! like the method=3 below).
                      ! This should be used as the default optimizer for rho.
            rhoMethod = 1
          CASE("NCG") ! This is conjugate gradient method that conserves electron number
                      ! systematically (not scaling).
            WRITE(message,*) "(input) Energy Minimization Method      : NCG"
            CALL WrtOut(6,message)
            rhoMethod = 2
          CASE("NBF")   ! BFGS density optimizer   (
            WRITE(message,*) "(input) Energy Minimization Method      : BFGS"
            CALL WrtOut(6,message)
            rhoMethod = 3
          !---------------------------------------------------
          CASE("STN") ! Select the sqrt truncated newton minimization
            rhoMethod = 4
          CASE("SCG") ! Select the SQRT Conjugate Gradient
            rhoMethod = 5
          CASE("HYB") ! Select conjugate gradient/newton hybrid
            rhoMethod = 6
          !---------------------------------------------------
          CASE("LOG") ! Select conjugate gradient/newton hybrid
            rhoMethod = 7
          !---------------------------------------------------
          CASE("ION", "GEO", "COO", "FOR") ! Select the ionic relaxation method 
            IF (rankGlobal==0) outputFinalGeometry = .TRUE.
            BACKSPACE inputUnit
            READ(inputUnit,*, IOSTAT=fileStatus) keyword, option, option
            IF (fileStatus/=0) THEN
              WRITE(*,*) 'METH ION format error. Format is: METH ION [OPTION]'
              STOP
            END IF

            CALL UPPERCASE(option)
            SELECT CASE (TRIM(option(1:3)))
              CASE("NON") 
                ionMethod = 0
              ! QuickMin
              CASE ("QUI")
                ionMethod = 2
              ! Conjugate gradient
              CASE("CON")
                ionMethod = 3
              ! Conjugate gradient(Chen Huang's revised)
              CASE("CG2") 
                ionMethod = 4
              ! BFGS(Steven add this)
              CASE("BFG") 
                ionMethod = 5
              CASE DEFAULT
                WRITE(message,*) 'Warning: Encountered unknown METH ION argument ', TRIM(option)
                CALL Error(6, message)
            END SELECT
          CASE("CEL")
            IF (rankGlobal==0) outputFinalGeometry = .TRUE.
            BACKSPACE inputUnit
            READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, tempInt
            IF (fileStatus/=0) THEN
              WRITE(message,*) 'METH CEL format error. Format is: METH CEL [INT]'
              CALL Error(6, message)
            END IF
            SELECT CASE(tempInt)
              CASE (1)
                cellRelax = tempInt
              CASE (2)
                cellRelax = tempInt
              CASE DEFAULT
                WRITE(message,*) 'Warning: Encountered unknown METH CEL argument (only 1, 2 and 3): ', tempInt
                CALL Error(6, message)
            ENDSELECT
          CASE DEFAULT
            WRITE(message,*) 'Warning: Encountered unknown METH argument ', TRIM(option)
            CALL Error(6, message)
          END SELECT

      !-----------------------------------------------------------------------------
      ! 07)  select which CG algorithm to use with METH NCG 
      CASE("DENC")  
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, tempInt
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'DENC format error. Format is: DENC [1 or 2]'
          CALL Error(6,message)
        END IF
        SELECT CASE(tempInt)
        CASE (1)
          cg_alg = "PR"
          CALL WrtOut(6," Polak-Ribiere CG algorithm for density optimization")
        CASE (2)
          cg_alg = "HZ"
          CALL WrtOut(6," Hager-Zhang CG algorithm for density optimization")
        CASE DEFAULT 
          WRITE(message,*) 'Warning: Encountered unknown DENC argument (only 1 or 2)', tempInt
          CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------
      ! 08) select which CG algorithm to use with METH ION CG2 
      CASE("IONC") !! the cg_type to be used in 
                   !! nonlinear conjugate gradient ion relaxation (CG2)
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, tempInt
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'IONC format error. Format is: IONC [1,2 or 3]'
          CALL Error(6,message)
        END IF
        SELECT CASE(tempInt)
        CASE (1)
          cg_type = 1
          CALL WrtOut(6," Polak-Ribiere CG algorithm for ionic optimization")
        CASE (2)
          cg_type = 2 
          CALL WrtOut(6," Hager-Zhang CG algorithm for ionic optimization")
        CASE (3)
          cg_type = 3 
          CALL WrtOut(6," Dai-Yuan CG algorithm for ionic optimization")
        CASE DEFAULT 
          WRITE(message,*) 'Warning: Encountered unknown IONC argument (only 1, 2 and 3)', tempInt
          CALL Error(6, message)
        END SELECT

      !------------------------------------------------------------------------------
      ! the number of storage steps in L-bfgs method
      CASE("MBFG") 
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, tempInt
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'MBFG format error. Format is: MBFG [INT]'
          CALL Error(6,message)
        END IF
        MBFGS = tempInt 
        WRITE(message,*) '(input) M_BFGS (used for BFGS density optimizer) is ', MBFGS
        CALL Wrtout(6,message)

      !-----------------------------------------------------------------------------  
      ! calculate force / stress
      CASE("CALC") 
        READ(inputUnit,*) keyword, option
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:3)))
        CASE("FOR") ! Calculate the forces.
          calOption(4) = .TRUE. 
          IF (rankGlobal==0) outputFinalForces = .TRUE.
        CASE("STR") ! Calculate the stress tensor.
          calOption(5) = .TRUE. 
          IF (rankGlobal==0) outputFinalStress = .TRUE.
        CASE DEFAULT
          WRITE(*,*) 'Warning: Encountered unknown CALCULATE argument ',&
            TRIM(option), '. Skipping.'
        END SELECT

      !-----------------------------------------------------------------------------
      ! 09) convergence based on energy or potential 
      CASE("CONV")
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, tempInt
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'CONV format error. Format is: CONV [1,2 or 3]'
          CALL Error(6,message)
        END IF
        SELECT CASE ( tempInt )
        CASE(1)
          conv_check(1:4) = "ENE " 
        CASE(2)
          conv_check(1:4) = "POT "
        CASE(3)
          conv_check(1:4) = "BOTH"
        CASE DEFAULT
          WRITE(message,*) 'Warning: Encountered unknown CONV argument (only 1, 2 and 3)', tempInt
          CALL Error(6, message)
        ENDSELECT
        WRITE(message,'(a,a)') " (input) Density Convergence is based on : ", conv_check(1:4)
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------
      ! 10) The stop criteria for total energy convergence in eV. 
      CASE("TOLE") 
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, tole
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'TOLE format error. Format is: TOLE [REAL VALUE]'
          CALL Error(6,message)
        END IF
        IF (tole<=0) THEN
          WRITE(message,*) 'Specified TOLE value must be positive.'
          CALL Error(6,message)
        END IF 
        WRITE(message,'('' (input) Tolerence of energy (Ha)        : '',E10.2)') tole
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------
      ! 11) tolerence for norm of potential.
      CASE("TOLP")
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, pot_tol
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'TOLP format error. Format is: TOLP [REAL VALUE]'
          CALL Error(6, message)
        END IF
        IF (pot_tol<=0) THEN
          WRITE(message,*) 'Specified TOLP value must be positive.'
          CALL Error(6, message)
        END IF 
        WRITE(message,'('' (input) Tolerence of Potential          : '',E10.2)') pot_tol
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------
      ! 12) The stop criteria for stress (cell relaxation)
      CASE("TOLS")
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, tols
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'TOLS format error. Format is: TOLS [REAL VALUE]'
          CALL Error(6, message)
        ENDIF
        IF (tols<=0) THEN
          WRITE(message,*) 'Specified TOLS value must be positive.'
          CALL Error(6, message)
        END IF
        WRITE(message,'('' (input) Tolerence of stress            : '',E10.2)') tols

      !-----------------------------------------------------------------------------
      ! 13) The stop criteria for force
      CASE("TOLF") 
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, forceCutoff 
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'TOLF format error. Format is: TOLF [REAL VALUE]'
          CALL Error(6, message)
        ENDIF
        IF (forceCutoff<=0) THEN
          WRITE(message,*) 'Specified TOLF value must be positive: ', forceCutoff
          CALL Error(6, message)
        END IF
        WRITE(message,'(a,ES12.4,a)') '(input) force cutoff is set to ', forceCutoff, ' Ha/bohr'
        CALL WrtOut(6,message)


      !-----------------------------------------------------------------------------
      ! 14) The stop criteria for magnetic moment
      CASE("TOLM")
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, tolm
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'TOLM format error. Format is: TOLM [REAL VALUE]'
          CALL Error(6, message)
        END IF
        IF(tolm<=0) THEN
          WRITE(message,*) 'Specified TOLM value must be positive.'
          CALL Error(6, message)
        ENDIF
        WRITE(message,'(a,ES12.4,a)') '(input) magnetic tolerence is set to ', tolm 
        CALL WrtOut(6,message)


      !-----------------------------------------------------------------------------
      ! 15) Sets timestep for ion minimization
      CASE("TIME") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, timeStep
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'TIME format error. Format is: TIME [REAL VALUE]'
          CALL Error(6,message)
        END IF

      !-----------------------------------------------------------------------------
      ! 16) use precondition or not 
      CASE("PREC") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'PREC format error. Format is: PREC [ON/OFF]'
          CALL Error(6,message)
        END IF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:2)))
        CASE("ON")
          usePreconditioner = .TRUE.
          WRITE(message,'(A)') ' (input) Precondition is turning on (only for METH NEW/MIX).'
          CALL WrtOut(6,message)
        CASE("OF")
          usePreconditioner = .FALSE.
        CASE DEFAULT
          WRITE(message,*) 'Warning: Encountered unknown PREC argument (ON/OFF): ', option
          CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------
      ! 17) Reset electron density to uniform for each ion optimization step
      CASE("RHOU")
        READ(inputUnit,*) keyword
        trashPreviousRho=.TRUE.
        WRITE(message,'(A)') ' (input) Use uniform electron gas for each new geometry.' 
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------
      ! 18) Reset electron density to uniform for each ion optimization step
      ! Specifying the number of spin-densities.
      CASE("SPIN") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, numSpin
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'Spin format error. Format is: Spin [1 or 2]'
          CALL Error(6, message)
        END IF
        IF (numSpin<=0) THEN
          WRITE(message,*) 'Number of spin-densities must be positive.'
          CALL Error(6, message)
        END IF
        IF (numSpin>2) THEN
          WRITE(message,*) 'Cannot handle more than two spins.'
          CALL Error(6, message)
        END IF

      !-----------------------------------------------------------------------------
      ! 19) Initial value of the total magnetic moment (in units of Mu_B)
      CASE("MAGM")
        READ(inputUnit,*, IOSTAT=fileStatus) keyword, magmom
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'MAGM format error. Format is: MAGM [REAL VALUE]'
          CALL WrtOut(6,message)
        END IF
        WRITE(message,'(A,F10.4)') ' (input) Set initial magnetization to be ', magmom  
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------
      ! 20) fix the magnetization or not. 
      CASE("FIXM")
        READ(inputUnit,*) keyword
        fixedmag = .TRUE.

      !-----------------------------------------------------------------------------
      ! 21) Specifying the kinetic energy functional to be used.
      CASE("KINE") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'KINE format error. Format is: KINE [OPTION]'
          CALL Error(6,message)
        END IF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:5)))
        CASE("THO", "TF") ! Thomas-Fermi(TF) functional. 
          kinetic = 1
          WRITE(message, *) '(input) TF KEDF to be used.'
        CASE("VON", "VW") ! Von Weiszacker(vW) functional.
          kinetic = 2
          WRITE(message, *) '(input) vW KEDF to be used.'
        CASE("TF+")       ! Functional alpha * TF + beta * vW.
          kinetic = 3
          WRITE(message, *) '(input) alpha * T + beta * vW KEDF to be used.'
        CASE("WT")        ! TF + vW + Wang-Teter functional.
          kinetic = 4
          WRITE(message, *) '(input) WT KEDF to be used.'
        CASE("WTV")       ! TF + vW + Wang-Teter functional for use in vacuum
          kinetic = 4
          WRITE(message, *) '(input) WT + vacuum is activated.'
          bvac = .TRUE.   ! this is the only flag tell PROFESS that we are regulating vacuum
        CASE("WGC")       ! TF + vW + Wang-Govind-Carter functional.
          kinetic = 5
          WRITE(message, *) '(input) WGC KEDF to be used.'
        CASE("WGV")         ! TF + vW + Wang-Govind-Carter functional.
          kinetic = 5       ! still 5, but the subroutine called by Calculator will
                            ! be WGCvacEnergy and WGCvacPot
          WRITE(message, *) '(input) WGC + vacuum is activated.'
          bvac = .TRUE.     ! this is the only flag tell PROFESS that we are using WGCvac not WGC
        CASE("LQ")        ! LQ Functional (Jeng Da)
          kinetic = 7
          WRITE(message, *) '(input) LQ KEDF to be used.'
        CASE("HQ")        ! HQ functional (Jeng-Da Chai)
          kinetic = 8
          WRITE(message, *) '(input) HQ KEDF to be used.'
        CASE("CAT")       ! CAT KEDF, see the reference in FunctinalKinetic.F90
          kinetic = 10
          WRITE(message, *) '(input) CAT KEDF to be used.'
        CASE("CAV")         ! TF + vW + Wang-Govind-Carter functional.
          kinetic = 10 
          WRITE(message, *) '(input) CAT + vacuum is activated.'
          bvac = .TRUE.     ! this is the only flag tell PROFESS that we are using CATvac
        CASE("HC")         ! Huang-Carter KEDF (2011) PRB
          kinetic = 11
          WRITE(message, *) '(input) Huang-Carter 10 KEDF to be used'
        CASE("DEC")         ! Density Decompose method, WGCD
          kinetic = 12
          WRITE(message, *) '(input) WGC-Decomposition KEDF to be used'
          bvac = .false.
        CASE("DEV")       ! Density Decompose method, WGCD, with vacuum
          kinetic = 12
          WRITE(message, *) '(input) WGC-Decomposition KEDF to be used'
          bvac = .true.
        CASE("GGA")       ! many GGA functionals
          kinetic = 15
          WRITE(message, *) '(input) GGA-type KEDF to be used'
          ! backward to read specific GGA functional
          BACKSPACE inputUnit
          READ(inputUnit, *, IOSTAT=fileStatus) keyword, message, option
          IF (fileStatus/=0) THEN
            WRITE(message,*) ' PARA GGA format error. Format is: PARA GGA [OPTION]'
            CALL Error(6, message)
          ELSE
            CALL Uppercase(option)
            SELECT CASE ( TRIM(option(1:4)) )
            CASE("TF")
              GGA_functional = 1
            CASE("VW")
              GGA_functional = 2
            CASE("TF+")
              GGA_functional = 3
            CASE("LLP")
              GGA_functional = 4
            CASE("OL2")
              GGA_functional = 5
            CASE("PW91")
              GGA_functional = 6
            CASE("TW")
              GGA_functional = 7
            CASE("PBE2")
              GGA_functional = 8
            CASE("LC94")
              GGA_functional = 9
            CASE("E00")
              GGA_functional = 10
            CASE("P92")
              GGA_functional = 11
            CASE("PW86")
              GGA_functional = 12
            CASE("DK")
              GGA_functional = 13
            CASE("OL1")
              GGA_functional = 14
            CASE("B86A")
              GGA_functional = 15
            CASE("B86B")
              GGA_functional = 16
            CASE("THAK")
              GGA_functional = 17
            CASE("DK87")
              GGA_functional = 18
            End Select
          ENDIF
        CASE("VW+") ! Density Decomposition using vW+G*TF KEDF 
          kinetic = 16
          WRITE(message,*) '(input) Density Decomposition using vW+G*TF KEDF'
        CASE("EVT") ! EvW using WT
          kinetic = 17
          WRITE(message,*) '(input) EvW KEDF based on WT KEDF' 
        CASE("EVC") ! EvW using WGC
          kinetic =18
          WRITE(message,*) '(input) EvW KEDF based on WGC KEDF' 
        CASE("MGP")        ! TF + vW + MGP functional.
          kinetic = 19
          WRITE(message, *) '(input) MGP KEDF to be used.'
        CASE DEFAULT
          WRITE(message,*) 'Warning: Encountered unknown KINE argument ', TRIM(option)
          CALL Error(6,message) 
        END SELECT
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------
      ! 22) Set the kinetic energy functional parameters.     
      CASE("PARA")
        READ(inputUnit, *, IOSTAT=fileStatus) keyword, option, tempReal
        IF (fileStatus/=0) THEN
          WRITE(message,*) ' PARA format error. Format is: PARA [OPTION] [REAL VALUE]'
          CALL Error(6, message)
        ENDIF
        CALL Uppercase(option)

        SELECT CASE (TRIM(option(1:4)))
          !------------------------------------------------------------------------------
          CASE("LAMB", "LAM") ! parameter for TF
            lambda = tempReal
          CASE("MU") ! parameter for VW
            mu = tempReal
          CASE("ALPH", "ALP") ! parameter for nonlocal KEDF
            alpha = tempReal
          CASE("BETA", "BET") ! parameter for nonlocal KEDF
            beta = tempReal
          CASE("GAMM", "GAM") ! parameter for nonlocal KEDF
            gamma = tempReal
          !------------------------------------------------------------------------------
          CASE("CATA") ! alpha for CAT KEDF
            ! cat_alpha = tempReal
            WRITE(message,*) '(input) cat_alpha must be zero, because we need the vW KEDF to', & 
                               ' be the lower bound of the KEDF!'
            CALL Error(6, message) 
          CASE("CATB") ! beta for CAT KEDF 
            ! cat_beta = tempReal
            WRITE(message,*) '(input) cat_beta must be 2/3, this is the optimum value'
            CALL Error(6,message)
          CASE("CATG") ! gamma for CAT KEDF 
            cat_gamma = tempReal
            Write(message,*) '(input) Gamma in CAT is set to be ',cat_gamma
            CALL WrtOut(6,message)
          !------------------------------------------------------------------------------
          CASE("RH0", "R0", "RHO0") ! Setting the value of rho0.
            rho0 = tempReal*bohr**3
            hold0 = .TRUE.
          CASE("RHS", "RS", "RHOS") ! Setting rho* to its value
            rhoS = tempReal*bohr**3
            holdS = .TRUE.
          CASE("MRHO")   ! multiplier of rhoS, mostly used in Steven's KEDF
            mrhos = tempReal
            WRITE(message,*) '(input) MRHOS is ', tempReal
            CALL WrtOut(6, message)
          CASE("LUME", "LE", "LUE") ! Setting the value of Lumgp exponent.
            LumgpExp = tempReal
            WRITE(message,*) '(input) LumgpExp set to ', LumgpExp
            CALL WrtOut(6, message)
          CASE("LUMF", "LF", "LUF") ! Setting the value of Lumgp factor.
            LumgpFactor = tempReal
            WRITE(message,*) '(input) LumgpFactor set to ', LumgpFactor
            CALL WrtOut(6, message)
          !------------------------------------------------------------------------------
          CASE("PBEC") ! Setting rho* to its value
            pbeCutoff = tempReal
            WRITE(message,*) '(input) PBE cutoff is set to ', pbeCutoff
            CALL WrtOut(6, message)
          !------------------------------------------------------------------------------
          CASE("AL5") ! for WGCT -5
            alpha5 = tempReal
          CASE("BE5") ! for WGCT -5
            beta5 = tempReal
          !------------------------------------------------------------------------------
          ! 2010 HUANG-CARTER KEDF
          CASE("RR")   ! ratio beteen bins for interpolating kernel
            refRatio = tempReal**(1._DP/3._DP)
            WRITE(message,'(2(a,ES12.4))') ' (input) ratio for interpolating kernels:', tempReal,', refRatio=', refRatio
            CALL WrtOut(6, message)
          CASE("HCLA") 
            hc_lambda_val = tempReal
            WRITE(message,*) '(input) hc_lambda is going to be set to = ', hc_lambda_val
            CALL Wrtout(6,message)
          CASE("CUTR") 
            cutrho = tempReal
            WRITE(message,*) '(input) density in HC KEDF smaller than ', cutrho, " wil set to 0"
            CALL Wrtout(6,message)
          !------------------------------------------------------------------------------
          ! 2012 CHEN HUANG's decomposition scheme.
          CASE("ATF")  ! Ratio of Thomas Fermi KEDF for transition metal
            aTF = tempReal
            WRITE(message,'(2(A,ES12.4))') ' (input) aTF = ', tempReal
            CALL WrtOut(6, message)
          CASE("BVW")  ! Ratio of Von Weizacker KEDF for transition metal
            bVW = tempReal
            WRITE(message,'(2(A,ES12.4))') ' (input) bVW = ', tempReal
            CALL WrtOut(6, message)
          !------------------------------------------------------------------------------
          ! 2012 JUNCHAO XIA's decomposition scheme.
          CASE("M")   ! Steven's KEDF
            shiftm = tempReal
            WRITE(message,*) '(input) M in WGCD is ', tempReal
            CALL WrtOut(6, message)
          CASE("SCFC")   ! Steven's KEDF
            scfc = tempReal
            WRITE(message,*) '(input) self-consistency for Fr is set to be ', tempReal
            CALL WrtOut(6, message)
          CASE("RHOC") ! Steven's KEDF 
            rhoc = tempReal
            WRITE(message,*) '(input) RHOC is set to be: ', tempReal
          !------------------------------------------------------------------------------
          ! 2014 Issac's KEDF.
          CASE("KLOC") ! for EvW KEDF (Shin, mohan add 09-26-12)
            kloc = tempReal
            WRITE(message,'(2(A,ES12.4))') ' (input) kloc in EvW KEDF is ', tempReal
            CALL WrtOut(6, message)
          CASE("ALOC") ! for EvW KEDF (Shin, mohan add 09-26-12)
            aloc = tempReal
            WRITE(message,'(2(A,ES12.4))') ' (input) aloc in EvW KEDF is ', tempReal
            CALL WrtOut(6, message)
          CASE("TOLK") ! for EvW KEDF
            tolk = tempReal
            WRITE(message,'(2(A,ES12.4))') ' (input) tolk in EvW KEDF is ', tempReal, ' (a.u.)'
            CALL WrtOut(6, message)
          !------------------------------------------------------------------------------
          ! GGA functionals
          CASE("CP")   
            CP = tempReal
            WRITE(message,*) '(input) GGA penalty coeff to be ', tempReal
            CALL WrtOut(6, message)
          CASE("VWMO")   
            model = NINT(tempReal) ! jmd: unhappy about this change but the whole parser seems rather duct-taped
            WRITE(message,*) '(input) vW+G*TF KEDF uses model ', tempReal
            CALL WrtOut(6, message)
          !------------------------------------------------------------------------------
          CASE DEFAULT ! No parameter by that name found.
            WRITE(message,*) 'Warning: Encountered unknown PARA argument ', TRIM(option)
            CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 23) orders of WGC
      CASE("WGCT") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, WGCT
        WRITE(message,'('' (input) Kinetic Energy Functional       : WGC '',I1)') WGCT
        CALL WrtOut(6, message)
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'WGCT format error. Format is: WGCT [INTEGER]'
          CALL Error(6, message)
        END IF
        IF (WGCT == 0) THEN  ! only use the 0th order of WGC
          firstOrderWGC = 1
        ELSE IF (WGCT == 1) THEN  ! use the 1st order of WGC
          firstOrderWGC = 1
        ELSE IF (WGCT == 2) THEN  ! use all the 2nd order of WGC
          firstOrderWGC = -1
        ELSE IF (WGCT == -3) THEN ! ddw/drho(r)^2 term will be set to zero in FillWGC routine
          firstOrderWGC = -1 
        ELSE IF (WGCT == -4) THEN ! ddW/dRho(r)dRho(r') term will be set to zero in FillWGC routine
          firstOrderWGC = -1
        ELSE IF (WGCT == -5) THEN ! use 2 coefficients for diagonal and non-diagonal terms in the 2nd order WGC
          firstOrderWGC = -1
        ELSE 
          CALL WrtOut(6,'WGCT has been assigned undefined value, error, stop')
          STOP
        ENDIF

      !-----------------------------------------------------------------------------  
      ! 24) the density at which the vacuum cutoff function quickly goes to zero
      CASE("RHOV") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, rhoV
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'RHOV format error. Format is: RHOV [REAL VALUE]'
          CALL Error(6, message)
        END IF
        WRITE(message,*) '(input) rhoV for vacuum is set to ', rhoV
        CALL WrtOut(6, message)

      !-----------------------------------------------------------------------------  
      ! 25) how fast the vacuum is going to zero 
      CASE("RHST") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, rhoStep
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'RHST format error. Format is: RHST [REAL VALUE]'
          CALL Error(6, message)
        END IF
        WRITE(message,*) '(input) rhoStep is set to ', rhoStep
        CALL WrtOut(6, message)

      !-----------------------------------------------------------------------------  
      ! 26) Hold rho0 or rhoS at their strating values even if the cell reshapes.
      CASE("HOLD")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'HOLD format error. Format is: HOLD [OPTION].'
          CALL Error(6, message)
        ENDIF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:4)))
          CASE("RH0", "R0", "RHO0") ! Setting the value of rho0.
            hold0 = .TRUE.
          CASE("RHS", "RS", "RHOS") ! Setting rho* to its value
            holdS = .TRUE.
          CASE DEFAULT
            WRITE(message,*) 'Warning: Encountered unknown HOLD argument ', TRIM(option)
            CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 27) Specifying the exchange-correlation functional.
      CASE("EXCH") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'EXCH format error. Format is: EXCH [LDA or PBE].'
          CALL Error(6, message)
        ENDIF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:3)))
        CASE("LDA") ! Local Density Approximation (default.)
          CALL WrtOut(6, ' (input) Exchange-Correlation Functional : LDA')
          exchangeCorrelation = 1
        CASE("PBE") ! Generalized Gradient Approximation, PBE version
          CALL WrtOut(6, ' (input) Exchange-Correlation Functional : PBE')
          exchangeCorrelation = 2
        CASE DEFAULT
          WRITE(message,*) 'Warning: Unknown EXCH option ',TRIM(option)
          CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 28) Specifying the Ewald parameters.
      CASE("EWAL") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'EWAL format error. Format is: EWAL [OPTION] [REAL VALUE].'
          CALL Error(6, message)
        ENDIF
        CALL Uppercase(option)
        ! Rewind so we can read it again!
        BACKSPACE inputUnit
        SELECT CASE (TRIM(option(1:3)))
          CASE("TOL") ! Ewald Error Tolerance
            READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, errorTolerance
            IF (fileStatus/=0) THEN
              WRITE(message,*) 'EWAL TOL format error. Format is: EWAL TOL [REAL VALUE].'
              CALL Error(6, message)
            END IF
            IF (errorTolerance <= 0) THEN
              WRITE(message,*) 'Specified EWAL TOL value must be positive: ',errorTolerance
              CALL Error(6, message) 
            END IF 
            errorTolerance = errorTolerance / hartreeToeV

          CASE("MAX") ! Ewald Max Real Cutoff
            READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, maxRealCutoff
            IF (fileStatus/=0) THEN
              WRITE(message, *) 'EWAL MAX format error. Format is: EWAL MAX [REAL VALUE].'
              CALL Error(6, message)
            END IF
            IF (maxRealCutoff <= 0) THEN
              WRITE(message,*) 'Specified EWAL MAX value must be positive: ', maxRealCutoff
              CALL Error(6, message)
            END IF
            maxRealCutoff = maxRealCutoff / bohr 

          CASE("ETA") ! Ewald Eta Increment
            READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, etaIncrement
            IF (fileStatus/=0) THEN
              WRITE(message,*) 'EWAL ETA format error. Format is: EWAL ETA [REAL VALUE].'
              CALL Error(6, message)
            END IF
            IF (etaIncrement <= 0) THEN
              WRITE(message,*) 'Specified EWAL ETA value must be positive: ', etaIncrement
              CALL Error(6, message)
            END IF 
            etaIncrement = etaIncrement * BOHR

          CASE("REC") ! Ewald Recip Cutoff Increment
            READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, recipCutoffIncrement
            IF (fileStatus/=0) THEN
              WRITE(message,*) 'EWAL REC format error. Format is: EWAL REC [REAL VALUE].'
              CALL Error(6, message)
            END IF
            IF (recipCutoffIncrement <= 0) THEN
              WRITE(message,*) 'Specified EWAL REC value must be positive:', recipCutoffIncrement
              CALL Error(6, message)
            END IF
            recipCutoffIncrement = recipCutoffIncrement * BOHR

          CASE DEFAULT
            WRITE(message,*) 'Warning: Unknown EWAL parameter ',TRIM(option)
            CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 29) Specifying whether to use cardinal b-spline approximations
      CASE ("SPLI")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, tempInt
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'SPLI format error. Format is: SPLI [OPTION] [INT].'
          CALL Error(6, message)
        END IF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:2)))
        CASE("AL") ! B-splines for ion-ion and ion-electron terms
          ieSpline = .TRUE.
          iiSpline = .TRUE.
          splineOrder = tempInt
          CALL WrtOut(6, ' (input) spline algorithm for ion-ion and ion-electron terms')
        CASE("IE") ! B-splines for ion-electron only (Choly-Kaxiras)
          ieSpline = .TRUE.
          iiSpline = .FALSE.
          splineOrder = tempInt
          CALL WrtOut(6, ' (input) spline algorithm for ion-electron term')
        CASE("II") ! B-splines for ion-ion only (particle-mesh Ewald)
          ieSpline = .FALSE.
          iiSpline = .TRUE.
          splineOrder = tempInt
          CALL WrtOut(6, ' (input) spline algorithm for ion-ion term')
        CASE("NO") ! Do not use approximation
          ieSpline = .FALSE.
          iiSpline = .FALSE.
        CASE DEFAULT
          WRITE(message,*) 'Warning: Unknown SPLI option:',TRIM(option)
          CALL Error(6, message)
          STOP
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 30) MD PARAMETER: md output path, mohan add 2013-01-17
      CASE("MDOP")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'MDOP format error. Format is: MDOP [PATH].'
          CALL Error(6, message)
        END IF
        md_output_path = TRIM(option) 
        
      !-----------------------------------------------------------------------------  
      ! 31) MD PARAMETER: Flag to restart Molecular Dynamics calculations.
      ! rstMD must be an integer
      ! 1) if rstMD < 0, for example, rstMD = -4, then the MD will restart
      ! from file ion.4.dat and vel.4.dat, this function is useful when
      ! you want to get consistent information compared to previous MD
      ! runs.
      ! 2) if rstMD == 0, we restart a new MD.
      ! 3) if rstMD == 1, restart from MD from files ion.restar and 
      ! vel.restart
      ! 4) if rstMD > 1, error.
      !------------------------------------------------------------------
      CASE("RSTM")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, rstMD
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'RSTM format error. Format is: RSTM [INT].'
          CALL Error(6, message)
        END IF
        IF (rstMD < 0)  THEN 
          CALL WrtOut(6, " (input) Restart MD from selected step.")
        ELSE IF(rstMD == 0) THEN
          CALL WrtOut(6," (input) Restart MD option is not used.")
        ELSE IF(rstMD == 1) THEN
          CALL WrtOut(6," (input) Restart MD from ion.restart and vel.restart")
        ELSE IF(rstMD > 1) THEN
          WRITE(message,*) 'Specified RSTM value must be <= 1.', rstMD 
          CALL Error(6, message)
        ENDIF

      !-----------------------------------------------------------------------------
      ! 32) MD PARAMETERS: NVE Molecular Dynamics
      CASE("NVEM")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, tempReal
        IF (fileStatus/=0) THEN
          WRITE(*,*) ' NVEM format error. Format is: NVEM [OPTION] [REAL VALUE].'
        ENDIF
        IF(tempReal < 0._DP) THEN
          WRITE(message,*) 'NVEM parameter must be positive: ', tempReal
          CALL Error(6, message)
        END IF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:4)))
          CASE("TEMP") ! Temperature
            temperature = tempReal * BOLTZMANN ! convert from [K] to a.u. and combine with boltzmann,
          CASE("TIME") ! Total time
            timeTot = tempReal / fundamentalTime  ! convert to a.u.
          CASE("DTIM") ! Time increment, in second
            dt = tempReal / fundamentalTime       ! convert to a.u.
          CASE DEFAULT
            WRITE(message,*) 'Warning: Unknown NVEM option: ',TRIM(option)
            CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 33) MD PARAMETERS: NVT thermostat
      CASE("NOSE") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, tempReal
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'NOSE format error. Format is: NOSE [OPTION] [REAL VALUE].'
          CALL Error(6, message)
        ELSE
        IF(tempReal < 0._DP) THEN
          WRITE(message,*) 'NOSE parameter must be positive: ', tempReal
          CALL Error(6, message)
        END IF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:4)))
          CASE("TEMP") ! Temperature
            WRITE(message,'('' (input) Temperature in NVT              : '',F10.2)') tempReal
            temperature = tempReal * BOLTZMANN ! convert to temperature in a.u.
            CALL WrtOut(6, message)
          CASE("TIME") ! Total time              
            WRITE(message,'('' (input) Total time in NVT               : '',E10.2)') tempReal
            CALL WrtOut(6, message)
            timeTot = tempReal / fundamentalTime ! convert to time in a.u.
          CASE("DTIM") ! Time interval
            WRITE(message,'('' (input) Time interval in NVT            : '',E10.2)') tempReal
            CALL WrtOut(6, message)
            dt = tempReal / fundamentalTime   ! convert to time in a.u.              
          CASE("QMAS") ! Q mass (in atomic units, mass * length**2)
            Qmass = tempReal 
          CASE("NRES") ! integer for nResn
            nResn = NINT(tempReal) 
            IF(nResn .GE. 1 ) THEN
              WRITE(message,'('' (input) NRESN         in NVT            : '',I3)') nResn
              CALL WrtOut(6,message)
            ELSE
              WRITE(message,*) "NOSE NRES can only be integer larger than 1."
              CALL Error(6, message)
            ENDIF
          CASE("NYOS") ! integer for nYosh
            nYosh = NINT(tempReal) 
            IF(nYosh==1 .OR. nYosh==3 .OR. nYosh==5) THEN
              WRITE(message,'('' (input) NYOSH         in NVT            : '',I3)') nYosh
              CALL WrtOut(6,message)
            ELSE
              WRITE(message,*) "NOSE NYOS can only be 1, 3 or 5."
              CALL Error(6, message)
            ENDIF
          CASE DEFAULT
            WRITE(message,*) 'Warning: Unknown NOSE option: ',TRIM(option)
            CALL Error(6, message)
        END SELECT
      END IF

      !-----------------------------------------------------------------------------  
      ! 34) MD PARAMETERS: NPT full cell relaxation
      CASE("NPTM") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option, tempReal
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'NPTM format error. Format is:NPTM [OPTION] [REAL VALUE].'
          CALL Error(6, message)
        ENDIF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:4)))
          CASE("TEMP") ! Temperature
            temperature = tempReal * boltzmann ! convert from [K] to a.u. and combine with boltzmann,
                                             ! no need for bolzmann in future
          CASE("TIME") ! Total time
            timeTot = tempReal / fundamentalTime  ! convert to a.u.
          CASE("DTIM") ! Time increment, in second
            dt = tempReal / fundamentalTime       ! convert to a.u.
          CASE("QMAS") ! thermostat mass (in atomic units, mass * length**2)
            Qmass = tempReal 
          CASE("BMAS") ! barostat mass (in atomic units, mass * length**2)
            bMass = tempReal 
          CASE("EXTP") ! exteranl pressure
            extPres = tempReal*auPressure   ! convert from Pa to a.u. pressure
          CASE("NRES") ! integer for nResn
            nResn = NINT(tempReal) ! jmd: unhappy about this change but the whole parser seems rather duct-taped
            IF(nResn .GT. 1) THEN
              WRITE(message,'('' (input) NRESN         in NPT            : '',I3)') nResn
              CALL WrtOut(6,message)
            ELSE
              WRITE(message,*) "NRESN can only be integer larger than 1."
              CALL Error(6, message) 
            ENDIF
          CASE("NYOS") ! integer for nYosh
            nYosh = NINT(tempReal)  ! jmd: unhappy about this change but the whole parser seems rather duct-taped
            IF(nYosh==1 .OR. nYosh==3 .OR. nYosh==5) THEN
              WRITE(message,'('' (input) NYOSH         in NPT            : '',I3)') nYosh
              CALL WrtOut(6,message)
            ELSE
              WRITE(message,*) "NYOSH can only be 1, 3 or 5."
              CALL Error(6,message)
            ENDIF
          CASE("TAUT") ! period of thermo stat vibration
            tau_thermo = tempReal/fundamentalTime
          CASE("TAUB") ! vibration period of barostat
            tau_baro   = tempReal/fundamentalTime
          CASE("CONS") ! vibration period of barostat
            constr_type = NINT(tempReal)  ! jmd: unhappy about this change but the whole parser seems rather duct-taped
          CASE DEFAULT
            WRITE(message,*) 'Warning: Unknown NPTM option: ',TRIM(option)
            CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 35) MD PARAMETER: rescale the initial read in velocity
      CASE("VELR")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword
        velRescale = .TRUE.
        WRITE(message,'(A)') ' (input) Rescale the velocities to target temperature when restarting MD.' 
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------  
      ! 36) MD PARAMETER: Mean Square Displacements are calculated beyond this
      ! step
      CASE("MSDI")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, tempInt
        msd_startTime = tempInt 
        WRITE(message,'(A,I5)') ' (input) Print out the MSD for each atom after step ', tempInt
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------  
      ! 37) MD PARAMETER: Mean Square Displacements are printed for each atom
      CASE("MSDA")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword
        doMSDAtom = .TRUE.
        WRITE(message,'(A)') ' (input) Print out the MSD for each atom during MD.'
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------  
      ! 38) Specifying the name of the geometry file.
      CASE("GEOM") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, geometryFile
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'GEOM format error. Format is: GEOM [FILE].'
          CALL Error(6, message)
        ENDIF
        WRITE(message,*) '(input) Read in geometry file from: ',geometryFile 
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------  
      ! 39) Specifying the name of the output file.
      CASE("OUTP") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, outputFile
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'OUTP format error. Format is: OUTP [FILE].'
          CALL Error(6, message)
        ENDIF
        WRITE(message,*) '(input) Output file is: ', outputFile 
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------  
      ! 40) Specifying the name of the error file.
      CASE("ERRO") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, errorFile
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'ERRO format error. Format is: ERRO [FILE].'
          CALL Error(6, message)
        ENDIF
        WRITE(message,*) '(input) Error file is: ', errorFile 
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------  
      ! 41) Specifying the name of the input density file.
      CASE("RHOF") 
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, densityFile
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'RHOF format error. Format is: RHOF [FILE].'
          CALL Error(6, message)
        ENDIF
        fileDensityFlag=.TRUE.
        WRITE(message,*) ' (input) Read in density file from: ',densityFile 
        CALL WrtOut(6,message)

      !-----------------------------------------------------------------------------  
      ! 42) Whether to use atomic density.
      CASE("RHOA")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword
        atomicDensityFlag=.TRUE.

      !-----------------------------------------------------------------------------  
      ! 43) Specifying whether there is a density file.
      CASE("OLDD", "OLDR") 
        READ(inputUnit, *) keyword
        fileDensityFlag=.TRUE.

      !-----------------------------------------------------------------------------  
      ! 44) Specify the number of time the input density file is read in
      CASE("RHOR")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, option
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'RHOR format error. Format is: RHOR [OPTION] [INT].'
          CALL Error(6, message)
        ENDIF
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:1)))
        CASE("X") ! 
            BACKSPACE inputUnit
            READ(inputUnit,*) keyword, option, xRepeat
        CASE("Y") ! 
            BACKSPACE inputUnit
            READ(inputUnit,*) keyword, option, yRepeat
        CASE("Z") ! 
            BACKSPACE inputUnit
            READ(inputUnit,*) keyword, option, zRepeat
        CASE DEFAULT
          WRITE(message,*) 'Warning: Encountered unknown RHOR argument ', TRIM(option)
          CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 45) Specifying the maximum size of the output density file
      CASE("RHOM")
        READ(inputUnit,*,IOSTAT=fileStatus) keyword, maxMB
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'RHOM format error. Format is: RHOM [REAL VALUE].'
          CALL Error(6, message)
        ENDIF

      !-----------------------------------------------------------------------------  
      ! 46) What do we want to print out for this calculation?
      CASE("PRIN")
        READ(inputUnit,*) keyword, option
        CALL Uppercase(option)

        SELECT CASE (TRIM(option(1:3)))
          CASE("MIN") ! Print various minimizer details
            BACKSPACE inputUnit
            READ(inputUnit, *) keyword, option, option, tempInt
            CALL UPPERCASE(option)
            IF (rankGlobal == 0) THEN
              SELECT CASE(TRIM(option(1:3)))
                CASE("RHO","DEN") ! We're talking about the density minimizer
                  outputMinimizeDensity = tempInt
                CASE("ION","GEO","COO","FOR") ! We're talking about the geometry minimizer
                  outputMinimizeGeometry = tempInt
                CASE DEFAULT
                  WRITE(message,*) 'Warning: Unknown PRIN MIN option: ',TRIM(option)
                  CALL Error(6, message)
              END SELECT
            END IF
          ! mohan add 2014-06-27 
          CASE("CEL") ! frequency for geometry output during cell relax
            IF(rankGlobal==0) THEN
              BACKSPACE inputUnit
              READ(inputUnit, *) keyword, option, tempInt
              celOutputFreq = tempInt ! Output geometry every tempInt steps
            ENDIF
          CASE("EWA") ! Print out parameters relating to the Ewald Summation
                      ! every time it is run
            IF (rankGlobal==0) outputEwald = .TRUE.
          CASE("STR") ! Print out the stress whenever computed
            IF (rankGlobal==0) outputFinalStress = .TRUE.
            ! Make sure the stress gets calculated.
            calOption(5) = .TRUE.
          CASE("RES") ! Prints restart information, such as the density at
                      ! each iteration, in case the job fails
            outputIntermediateDensity = .TRUE.
          CASE("RHO","DEN") ! Print out the electron density at the end of the calculation
            outputFinalDensity = .TRUE.
          CASE("ION","GEO","COO") ! frequency for geemetry output during ion relax
            IF (rankGlobal==0) THEN
              outputFinalGeometry = .TRUE.
              BACKSPACE inputUnit
              READ(inputUnit, *) keyword, option, tempInt
              geoOutputFreq = tempInt ! Output geometry every tempInt steps
            ENDIF
          CASE("FOR") ! Print out the total forces in a seperate file
            IF (rankGlobal==0) outputFinalForces = .TRUE.
            ! Ensure forces are computed so we don't report ghosts.
            calOption(4) = .TRUE.
            BACKSPACE inputUnit
            READ(inputUnit, *) keyword, option, option
            CALL UPPERCASE(option)
            IF (rankGlobal==0) THEN
              SELECT CASE (TRIM(option(1:3)))
                CASE("SEP")
                  outputForcesSeparate = .TRUE.
                CASE("SAM")
                  outputForcesSeparate = .FALSE.
                CASE DEFAULT
                  WRITE(message,*) 'Warning: Unknown PRINT FORCES option ', TRIM(option)
                  CALL Error(6, message)
              END SELECT
            END IF
          CASE("MDI") ! MD PARAMETERS: output ion molecular dynamics information for restarting MD every dump_md_freq steps
            BACKSPACE inputUnit
            READ(inputUnit,*) keyword, option, dump_md_freq
            WRITE(message,*) "(input) dump MD information at frqeuency: ", dump_md_freq
            CALL WrtOut(6,message)
        
          ! test only, not shown in PROFESS 3 manual
          CASE("KER") ! Print out the kernel after the kernels have been generated
            outputKernel = .TRUE.
        !  CASE("POT") ! Print out the electronic potential at the end of the calculation
        !    outputFinalPotential = .TRUE.
          CASE DEFAULT
            WRITE(message,*) 'Warning: Unknown PRIN option ',TRIM(option)
            CALL Error(6, message)
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 47) Transistion state finder.  When the transition state finder is
      ! integrated into this code, this part will CHANGE.
      CASE("TRAN") 
        READ(inputUnit,*) keyword, option
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:2)))          
        CASE("ON")
          IF (rankGlobal==0) outputTransitionState = .TRUE.
        CASE("OF")
          IF (rankGlobal==0) outputTransitionState = .FALSE.
        END SELECT

      !-----------------------------------------------------------------------------  
      ! 48) 
      CASE("NNDI") 
        READ(inputUnit,*, iostat=fileStatus) keyword, nnDist
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'Error in reading nnDist keyword.'
          CALL Error(6, message)
        END IF
        IF(nnDist<0) THEN
          WRITE(message,*) 'Turning off computation of nearest neighbor distance.'
          CALL WrtOut(6,message)
        ELSE
          nnDist = nnDist/bohr
          WRITE(message,'('' (input) Threshold of atom distance      : '',F16.8,'' Bohr'')') nnDist 
          CALL WrtOut(6,message)
        ENDIF

      !-----------------------------------------------------------------------------  
      ! 49) 
      CASE("NNBI")
        ! Compute nearest atom distance using binning method?
        READ(inputUnit,*, iostat=fileStatus) keyword, tempReal
        IF (fileStatus/=0) THEN
          WRITE(message,*) 'Error in reading NNBI keyword.'
          CALL Error(6, message)
        END IF
        IF(tempReal<0) THEN
          WRITE(message,*) '(Input) NNDist compuated without binning.'
          checkNNDist_bin = .false.
          CALL WrtOut(6,message)
        ELSE
          WRITE(message,*) '(Input) NNDist compuated with binning.'
          checkNNDist_bin = .true.
          CALL WrtOut(6,message)
        ENDIF


      !-----------------------------------------------------------------------------  
      ! 51) Use Density Decomposition Method 
      CASE("UDDM")
        READ(inputUnit, *)
        do_den_dec=1
        CALL WrtOut(6," (input) Use Density Decomposition Method (DDM)")

      !-----------------------------------------------------------------------------  
      ! BELOW WILL NOT BE WRITTEN INTO MANUAL
      ! FOR TESTS ONLY
      !-----------------------------------------------------------------------------  
      ! only used in RhoOptSTN
      CASE("CHEA") 
        READ(inputUnit,*) keyword, option
        CALL Uppercase(option)
        SELECT CASE (TRIM(option(1:2)))
        CASE("ON")
          cheatPot = .TRUE.
        CASE("OF")
          cheatPot = .FALSE.
        END SELECT

      ! MD PARAMETER: If we change temperature during MD simulation?
      CASE("CHGT") 
        READ(inputUnit, *, iostat=fileStatus) keyword, tempInt
        IF (tempInt > 0)  THEN 
          fixTemperature = tempInt
          CALL WrtOut(6," (input) Temperature will be changed during the MD")
        ELSE
          fixTemperature = tempInt
          CALL WrtOut(6," (input) Temperature will NOT be changed during the MD")
        ENDIF

      ! DEBUG PARAMETERS: Output file generated by each core.
      CASE("ALLF")
        READ(inputUnit, *, iostat=fileStatus) keyword, tempInt
        IF(tempInt<0 .OR. tempInt>1) THEN
          WRITE(message,*) 'ALLF can only be chosen 0 or 1'
          CALL WrtOut(6,message)
          STOP
        ENDIF
        output_all_files=tempInt
        CALL WrtOut(6," (input) Output all files.")

      ! *************Insert new keywords here.*****************

      CASE DEFAULT

        ! Comment line. Ignore.
        IF(keyword(1:1) == "#") THEN
          READ(inputUnit,*, iostat=fileStatus) keyword

        ELSE
          WRITE(*,*) 'Warning: Unknown keyword ', TRIM(keyword), '. code stopped.'
          STOP
          READ(inputUnit,*) keyword
        END IF
    END SELECT
  END DO

  ! Closing the file, we're done with it.
  CLOSE (inputUnit, status="keep") 


  ! These variables are needed for OUTPUT.  Initialize them here.  But
  ! later, these may be initialized elsewhere.
  IF (rankGlobal==0) THEN
    exchangeCorrelationFunct = exchangeCorrelation
    outputRhoMethod = rhoMethod
    outputIonMethod = ionMethod
    outputOptionGeom = outputOptionGeom .OR. calOption(2) .OR. calOption(4)
  END IF

  RETURN


END SUBROUTINE ReadOptions


SUBROUTINE CheckOptions
!------------------------------------------------------------------------------
! DESCRIPTION:
! This routine checks all the read in parameters to see whether they are
! in a reasonable range.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/04/2012  File created.  (mohan)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

  CHARACTER(LEN=500) :: message(10)
  INTEGER :: num

  num = 0

  ! Check to see if some mandatory parameters were specified.
  IF (energyCutoff<=0._DP .AND. gridSpacing <=0._DP) THEN
    IF (rankGlobal == 0) WRITE(*,*) 'The cutoff energy or gridpt density must be specified in the ',&
      'input file, but the ECUT or GDEN keyword could not be located. Leaving.'
    STOP
  END IF

  IF (energyCutoff>0._DP .AND. gridSpacing>0._DP) THEN
    IF (rankGlobal==0) WRITE(*,*) 'The cutoff energy and gridpt density cannot both be specified! Leaving.'
    STOP
  END IF

  IF ( (bvac .EQV. .TRUE.) .AND. &
       (rho0<0._DP .OR. rhoS<0._DP)) THEN
    IF (rankGlobal == 0) &
      WRITE(*,*) &
        "If there's vacuum, parameters rho0 and rhoS must be specified in inpt!"
    STOP
  END IF

  
  IF( mdType .EQ. 1 ) THEN

    IF( temperature < 0.0_DP ) THEN
      num = num + 1
      WRITE(message(num),*) " Temperature (e.g. 'NOSE TM 300.0') should be > 0, TM=",temperature
    END IF

    IF( timeTot < 0.0_DP ) THEN 
      num = num + 1
      WRITE(message(num),*) " Total time (e.g. 'NOSE TI 1.0e-12') should be > 0, TI=",timeTot
    END IF

    IF( dt < 0.0_DP ) THEN 
      num = num + 1
      WRITE(message(num),*) " Delta time (e.g. 'NOSE DT 0.1e-15') should be > 0, DT=",dt
    END IF

    IF( Qmass < 0.0_DP ) THEN 
      num = num + 1
      WRITE(message(num),*) " Qmass (e.g. 'NOSE QM 10.0') should be > 0, QM=",Qmass
    END IF
  
  ENDIF


  IF ( nnDist < 0 ) THEN
    WRITE(outputUnit,'(a)') " nnDist < 0, skipped nearest neighbour checking."
  ENDIF

  ! If there are any errors, print the message and stop.
  IF ( num .GE. 1 ) THEN
    CALL Error(6,message,num) 
  ENDIF

END SUBROUTINE CheckOptions


END MODULE ReadInputFile 

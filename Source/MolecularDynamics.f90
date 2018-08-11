MODULE MolecularDynamics
!----------------------- -------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE MolecularDynamics
!     |_SUBROUTINE RestartMD
!     |_SUBROUTINE RemoveMovementOfCenterOfMass
!     |_SUBROUTINE GetAtomKE
!     |_SUBROUTINE MakeIntCoeff
!     |_SUBROUTINE InitRandomSeed
!     |_SUBROUTINE InitVelocity
!     |_SUBROUTINE MonitorMeanSquareDisplacement
!     |_SUBROUTINE OutputIonVelocityAndStates
!     |_SUBROUTINE OutputMDGeom
!     |_SUBROUTINE ReadNewTemp
!     |_SUBROUTINE AdjustCenterOfMass
!     |_SUBROUTINE CalFreedom
!
! DESCRIPTION:
!   This module contains molecular dynamics subroutines.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!  Ref[1]: Martyna et. al., "Explicit reversible integrators for extended system
!           dynamics," Mol. Phys. 87 (1996), 1117.
!   Allen and Tildesley, "Computer Simulation of Liquids" (predictor-corrector)
!   Frenkel and Smit, "Understanding Molecular Simulation".
!
!------------------------------------------------------------------------------
! REVISION LOG:
! Written by Linda Hung, Chen Huang, and Mohan Chen
!------------------------------------------------------------------------------

  USE Constants, ONLY: fundamentalTime
  USE Constants, ONLY: PI
  USE Constants, ONLY: DP
  USE Constants, ONLY: BOHR
  USE Constants, ONLY: boltzmann
  USE Constants, ONLY: hartreeToeV 
  
  USE MathFunctions, ONLY : Inverse
  USE MathFunctions, ONLY : Vecmul

  USE Output, ONLY: WrtOut
  USE OutputFiles

  USE MPI_Functions, ONLY : rankGlobal
  USE CellInfo, ONLY : Cell

  USE Timer, ONLY: TimerStart
  USE Timer, ONLY: TimerStop
  USE Timer, ONLY: stopwatch

  USE SYS, ONLY: rhoR
  ! real space density

  IMPLICIT NONE

  REAL(KIND=DP) :: tau_thermo = 0.d0
  ! The thermostat fluctuation period, for Al fcc at 300K ~ 1000fs
  !
  REAL(KIND=DP) :: tau_baro   = 0.d0
  ! The barostat   fluctuation period, for Al fcc at 300K ~ 500 fs
  !
  REAL(KIND=DP) :: timeTot    = -1.d0
  ! Total time of simulation
  !
  REAL(KIND=DP) :: Qmass = -1._DP
  ! Inertia of extended system variable
  !
  REAL(KIND=DP) :: dt = -1._DP
  ! Time increment (hbar/E_hartree)
  !
  REAL(KIND=DP) :: temperature = -1.d0
  ! Temperature (in Hartree, 1 Hartree ~ 3E5 K)
  !
  LOGICAL :: doMSDAtom = .FALSE.
  ! compute the mean square displacement for each atom.
  !
  REAL(KIND=DP) :: msd_startTime = 10
  ! beyond this time, msd is going to be computed
  !
  REAL(KIND=DP) :: msd = -1.D0
  !
  REAL(KIND=DP) :: diffuCoeff = 0.D0
  !

  REAL(KIND=DP), ALLOCATABLE :: cartNoWrap(:,:)  
  ! cartensian coordinates of atoms, *not* wrapped
  !
  REAL(KIND=DP), ALLOCATABLE :: msd_coordst(:,:)
  REAL(KIND=DP), ALLOCATABLE :: msd_coords0(:,:)
 
  INTEGER :: nResn=3   
  ! Coarseness of Trotter factorization of Nose-Hoover thermostat
  !
  INTEGER :: nYosh=3   
  ! Order of thermostat scheme corresponding with "intCoeff"
  !
  INTEGER :: output_stress_period = 1 
  ! The period to output stress
  !
  INTEGER :: dump_md_freq = 1         
  ! The period to dump MD information for monitoring and restarting MD
  !
  INTEGER :: fixTemperature = 1       
  ! The period to read in new temperature during MD run.
  !
  CHARACTER(LEN=500) :: md_output_path = '.' 
  ! output directory of md files: .ion .vel
  !
  INTEGER :: rstMD = 0  ! 1 : restart MD, vel.restart and ion.restart files will be read
                        ! 0 : not restart from ion and vel files
                        ! -n: restart from vel.n and ion.n in directory
                        ! "md_output_path"

  INTEGER :: startStep=0 
  ! the first step number of MD
  !
  LOGICAL :: velRescale = .false. 
  ! whether to rescale velocities at the first step
  !
  REAL(KIND=DP) :: freedom
  ! number of freedoms
  !
  TYPE(stopwatch):: watch, watch2
  ! timing variables
  !

CONTAINS


SUBROUTINE RestartMD(md_type,vel,xLogS,vLogS,vBoxG)
!---------------------------------------------------------------------------
! DESCRIPTION:
!   Read in MD restart informations from
!   (1) ion.restart (box and ion position info) and 
!   (2) vel.restart (velocities of atoms and box's and thermostats, as well
!                    as thermostat's positions)
!   Note:
!    two file names, ion.restart vel.restart are fixed, cannot be changed
!    and the file vel.restart is just the file dumped from last MD run, you 
!    cannot modify it. (simply to rename the vel.dat to vel.restart)
!    
!    You can add lines of comments into the ion.restart file, but
!    you cannot mess up the blocks which contain ion position and lattice
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   1) We should enable the restart MD simulations with fractional coordinates
! of atoms.
!   2) Not sure how to consider fixing the mass center, do we need to consdier
! thermostat?
!   3) When restart, we want to read in potential information and fixed atom
! information.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2014-02-20 Mohan add warnings for reading ion file. 
!------------------------------------------------------------------------------
  USE CONSTANTS, ONLY: BOHR
  USE CONSTANTS, ONLY : BOLTZMANN
  USE CellInfo, ONLY: RefreshLattice
  USE MathFunctions, ONLY : Inverse
  USE SYS, ONLY: frozenIon

  IMPLICIT NONE
  
  !! >> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(IN) :: md_type ! 1: NVT, 2: NPT, 3: NVE 

  REAL(KIND=DP), INTENT(OUT) :: xLogS ! thermostat's position
  REAL(KIND=DP), INTENT(OUT) :: vLogS ! thermostat's speed
  REAL(KIND=DP), INTENT(OUT) :: vel(3,cell%numIon)  ! velocity of each atom, unit is a.u.
  REAL(KIND=DP), OPTIONAL :: vBoxG(3,3)        ! barostat's velocity

  INTEGER :: & 
    ii, &
    state = 0, &
    ion_file = 110,  & 
    vel_file = 119
  
  REAL(KIND=DP) :: & 
    ionPos(3),  &
    invMat(3,3)

  CHARACTER(LEN=500) :: & 
    message, &
    eleName, &
    line, line1, line2

  REAL(KIND=DP) :: ke
  !
  INTEGER :: openStatus1 
  ! check if the files exist
  !
  INTEGER :: openStatus2
  !
  INTEGER :: jj
  ! integer for x,y,z index
  !
  CHARACTER(LEN=50) :: filename_ion
  CHARACTER(LEN=50) :: filename_vel
  CHARACTER(LEN=30) :: tmp_step

  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           RESTART MOLECULAR DYNAMICS"


  !---------------------------------------
  ! Initialize the name of restart files.
  !---------------------------------------
  ! mohan add 2013-07-11
  ! rstMD is the input parameter,
  ! if equals 1, restart from current directory,
  ! the ionic coordinates are from ion.restart
  ! while velocities are from vel.restart
  IF(rstMD==1) THEN
    filename_ion="ion.restart"
    filename_vel="vel.restart"
  ! If the value is negative, start the calculation
  ! from selected step (inside the md_output_path),
  ! this function is usually used to keep on doing
  ! MD simulations from broken point.
  ELSE IF(rstMD<0) THEN
    WRITE(tmp_step, *) abs(rstMD) 
    filename_ion = TRIM(md_output_path) // "/ion." // TRIM(ADJUSTL(tmp_step)) // ".dat"
    filename_vel = TRIM(md_output_path) // "/vel." // TRIM(ADJUSTL(tmp_step)) // ".dat"
  ENDIF

  WRITE(outputUnit,*) "Restart MD, geometry from file ", filename_ion
  WRITE(outputUnit,*) "Restart MD, velocity from file ", filename_vel


  !---------------------------------------
  ! open the ionic file and velocity file
  !---------------------------------------
  OPEN(UNIT=ion_file, FILE=filename_ion, ACTION="read", IOSTAT=openStatus1)
  OPEN(UNIT=vel_file, FILE=filename_vel, ACTION="read", IOSTAT=openStatus2)
  IF(openStatus1.NE.0) THEN
    WRITE(message,*) filename_ion, " file doesn't exist."
    CALL Error(outputUnit, message)
  ENDIF
  IF(openStatus2.NE.0) THEN
    WRITE(message,*) filename_vel, " file doesn't exist."
    CALL Error(outputUnit, message)
  ENDIF

  ! >>> FUNCTION <<<!


  !----------------------------------
  ! begin to read in ion coordinates
  !----------------------------------
  state = 0
  ! For ion_file
  CALL WrtOut(6," (RestartMD) Reading positions ...")
  DO WHILE ( .TRUE. )

    !------------------------------------------------------
    ! read in the lattice vectors of the cell,
    ! the unit of lattice vectors are in Ansstrom,
    ! transfer to Bohr here.
    !------------------------------------------------------
    IF (state==1) THEN 
      READ(ion_file,*) cell%cellReal(:,1)
      READ(ion_file,*) cell%cellReal(:,2)
      READ(ion_file,*) cell%cellReal(:,3)
      ! read in cell is in Angstrom unit,
      ! translate it to Bohr unit
      CALL RefreshLattice(cell%cellReal/BOHR)
      ! Inverse matrix is used to transfer the ion
      ! coordinates to direct coordinates 'cell%ionTable(ii)%coord'
      invMat = INVERSE(cell%cellReal)
      state = 0
      CALL WrtOut(6," (RestartMD) Finish BLOCK LATTICE_CART")
      CYCLE
    ENDIF

    !--------------------------------------------------------------
    ! restart from Cartesian coordinates (unit of atom coordinates
    ! is Angstrom) 
    !--------------------------------------------------------------
    IF (state==21 .OR. state==22) THEN 
      DO ii=1,cell%numIon
        READ(ion_file,*) eleName, ionPos(:)
        IF(state==21) THEN
          ! translate the 'ionPos' from Angstrom into Bohr.
          cell%ionTable(ii)%coord(1:3) = MATMUL(invMat,ionPos/BOHR)
        ELSE IF(state==22) THEN
          ! fractional coordinates
          cell%ionTable(ii)%coord(1:3) = ionPos
        ENDIF
      ENDDO

      ! in case that ion position is out of the primitive cell, move all images
      ! into the cell:
      DO ii=1, cell%numIon
        cell%ionTable(ii)%coord(1) = cell%ionTable(ii)%coord(1) - FLOOR(cell%ionTable(ii)%coord(1))
        cell%ionTable(ii)%coord(2) = cell%ionTable(ii)%coord(2) - FLOOR(cell%ionTable(ii)%coord(2))
        cell%ionTable(ii)%coord(3) = cell%ionTable(ii)%coord(3) - FLOOR(cell%ionTable(ii)%coord(3))
      ENDDO

      state = 0
      CALL WrtOut(6," (RestartMD) Finish BLOCK POSITIONS_CART.")
      EXIT
    ENDIF

    READ(ion_file,*) line

    IF( line(1:6) .EQ. "%BLOCK") THEN 
      BACKSPACE (ion_file)
      READ(ion_file,*) line1,line2
      IF (line2(1:12) == "LATTICE_CART") THEN 
        CALL WRTOUT(6," (RestartMD) Reading LATTICE_CART block ...")
        state = 1   ! Reading %BLOCK LATTICE_CART block
        CYCLE
      ELSE IF( line2(1:14) == "POSITIONS_CART" ) THEN 
        CALL WRTOUT(6," (RestartMD) Reading POSITIONS_CART block ...")
        state = 21   ! Reading %BLOCK POSITIONS_CART block
        CYCLE
      ELSE IF( line2(1:14) == "POSITIONS_FRAC" ) THEN 
        CALL WRTOUT(6," (RestartMD) Reading POSITIONS_CART block ...")
        state = 22   ! Reading %BLOCK POSITIONS_CART block
        CYCLE
      ELSE 
        WRITE(message,*) "Restart MD: only supports&
        LATTICE_CART and POSITIONS_CART in ion file."
        CALL Error(6, message)
      ENDIF
      CYCLE
    ENDIF

  ENDDO ! for ion_file

  CALL WRTOUT(6," (RestartMD) Finish reading ion positions.")
  CALL WRTOUT(6," ")



  !----------------------------------------------
  ! Read in vel.restart file
  !----------------------------------------------

  !----------
  ! NVT Case 
  !----------
  IF (md_type == 1)  THEN
    CALL WrtOut(6," (RestartMD) Reading velocities (NVT) ...")
    READ(vel_file,*) eleName  ! discard the first line
    READ(vel_file,*) xLogS    ! does not influence the center of mass of the system
    READ(vel_file,*) vLogS    ! does not influence the center of mass of the system
    WRITE(message,'(A,ES15.6)') ' Read in thermostat position xLogS ', xLogS
    CALL WrtOut(6,message)
    WRITE(message,'(A,ES15.6)') ' Read in thermostat speed    vLogS ', vLogS
    CALL WrtOut(6,message)
    READ(vel_file,*) eleName  ! discard barostat velocity
    READ(vel_file,*) eleName  ! discard barostat velocity
    READ(vel_file,*) eleName  ! discard barostat velocity
    READ(vel_file,*) eleName  ! discard this first line (string)

    ! read in velocity, unit is ??
    DO ii=1,cell%numIon
      READ(vel_file,*) vel(1:3,ii)
    ENDDO
    CALL WRTOUT(6," (RestartMD) Finish reading velocities (NVT).")
  ENDIF

  !----------
  ! NPT Case 
  !----------
  IF (md_type == 2)  THEN
    CALL WrtOut(6," (RestartMD) Reading velocities (NPT)...")
    READ(vel_file,*) eleName  ! discard the first line
    READ(vel_file,*) xLogS    ! does not influence the center of mass of the system
    READ(vel_file,*) vLogS    ! does not influence the center of mass of the system
    WRITE(message,'(A,ES15.6)') ' Read in thermostat position xLogS ', xLogS
    CALL WrtOut(6,message)
    WRITE(message,'(A,ES15.6)') ' Read in thermostat speed    vLogS ', vLogS
    CALL WrtOut(6,message)
    READ(vel_file,*) vBoxG(1,1), vBoxG(1,2), vBoxG(1,3)    ! vBoxG(1,:)
    READ(vel_file,*) vBoxG(2,1), vBoxG(2,2), vBoxG(2,3)    ! vboxG(2,:)
    READ(vel_file,*) vBoxG(3,1), vBoxG(3,2), vBoxG(3,3)    ! vBoxG(3,:)
    ! read in velocity for each atom.
    READ(vel_file,*) eleName   ! discard first line
    DO ii=1,cell%numIon
      READ(vel_file,*) vel(1:3,ii)
    ENDDO
  ENDIF


  !----------
  ! NVE Case 
  !----------
  IF (md_type == 3)  THEN
    CALL WrtOut(6," (RestartMD) Reading velocities (NVE) ...")
    READ(vel_file,*) eleName  ! discard MD step 
    READ(vel_file,*) eleName  ! discard thermostat position
    READ(vel_file,*) eleName  ! discard terhmostat velocity
    READ(vel_file,*) eleName  ! discard barostat velocity
    READ(vel_file,*) eleName  ! discard barostat velocity
    READ(vel_file,*) eleName  ! discard barostat velocity
    ! read in velocity for each atom.
    READ(vel_file,*) eleName  ! discard first line 
    DO ii=1,cell%numIon
      READ(vel_file,*) vel(1:3,ii)
    ENDDO
    CALL WrtOut(6," (RestartMD) Finish reading velocities (NVE).")
  ENDIF


  ! mohan add 2014-07-29
  IF( ALLOCATED(frozenIon) ) THEN
    CALL WrtOut(6," (RestartMD) Set velocities to zero for frozen ions.")
    DO ii=1, cell%numIon
      DO jj=1, 3
        IF(frozenIon(ii,jj) .EQV. .TRUE.) THEN
          vel(jj, ii) = 0.0_DP
        ENDIF
      ENDDO
    ENDDO
  ENDIF

! mohan read in vBoxG on 2013-07-11
  ! We set zeros to all the 
  ! thermostat and barostat related
!  IF (PRESENT(vBOXG)) vBoxG = 0.D0
!  xLogS = 0.d0
!  vLogS = 0.d0

  ! add md_type=1 2014-11-14 by mohan
  IF (md_type ==1 .or. md_type == 2 .or. md_type == 3) THEN
    CALL RemoveMovementOfCenterOfMass(vel, .TRUE.)
    ! can not do this right now because of Ewald setup or things like that.
    !CALL AdjustCenterOfMass(.TRUE.) 
  ENDIF

!------------------------------------------------------------------------------------
! mohan add 2013-07-11
!  ! Rescale all velocities to start to exact temperature desired
  CALL GetAtomKE(vel,ke)

  WRITE(outputunit,*) "Initial temperature from velocities is ", &
    ke/(1.5_dp*REAL( cell%numIon,kind=DP)) /BOLTZMANN, " Kelvin"

  
  IF(velRescale .eqv. .TRUE.) THEN
    WRITE(outputUnit,*) "scale the velocities to desired starting temperature ", temperature/BOLTZMANN
    vel = vel * SQRT(temperature*1.5_dp * REAL(cell%numIon,kind=DP) / ke)
  ENDIF
!------------------------------------------------------------------------------------

  CLOSE(ion_file)
  CLOSE(vel_file)

  CALL WrtOut(6," (RestartMD) Finish")

  RETURN

END SUBROUTINE RestartMD


SUBROUTINE GetAtomKE(vel,ke)
!---------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the classical kinetic energy of a group of atoms.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------

  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN)   :: vel(3,size(cell%ionTable,1))  ! Velocities of atoms
  REAL(KIND=DP), INTENT(OUT)  :: ke ! kinetic energy we want.

  !>> INTERIOR VARIABLES <<!
  REAL(KIND=DP) :: mass
  INTEGER  :: ion
  
  !>> INITIALIZATION <<!
  !>> FUNCTION BODY <<!
  ke = 0._DP

  ! kinetic energy = \sum_{i} 1/2*m_{i}*(vx^2+vy^2+vz^2)
  DO ion = 1, SIZE(cell%ionTable)
    mass = cell%elementTable(cell%ionTable(ion)%elementID)%mass
    ke = ke + 0.5_dp * mass * SUM(vel(:,ion)**2)
!   WRITE(outputUnit,*) vel(1,ion), vel(2,ion), vel(3,ion)
  END DO

  RETURN

END SUBROUTINE GetAtomKE


SUBROUTINE MakeIntCoeff(nYosh,intCoeff)
!---------------------------------------------------------------------------
! DESCRIPTION:
!  Initialize the Yoshida/Suzuki coefficient
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nYosh ! approximation order
  REAL(KIND=DP),INTENT(OUT) :: intCoeff(nYosh)

  CALL Title("MolecularDynamics::MakeIntCoeff")

  SELECT CASE (nYosh)
    CASE(1)
      intCoeff(1) = 1._DP
    CASE(3)
      intCoeff = 1._DP/(2._DP-2._DP**(1._DP/3._DP))
      intCoeff(2) = 1._DP - 2._DP*intCoeff(2)
    CASE(5)
      intCoeff = 1._DP/(4._DP-4._DP**(1._DP/3._DP))
      intCoeff(3) = 1._DP - 4._DP*intCoeff(3)
    CASE DEFAULT
      STOP "ERROR: NYOSH MUST BE 1, 3, or 5.  STOPPING."
  END SELECT

  RETURN
END SUBROUTINE MakeIntCoeff


SUBROUTINE InitRandomSeed
!---------------------------------------------------------------------------
! DESCRIPTION:
!  Initialize the random seed used for random number generator.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------
  IMPLICIT NONE

  !! >>> local vars <<< !! 
  INTEGER :: seedSize, & 
             errorFlag, & 
             clock, i 
  INTEGER, ALLOCATABLE:: randomSeed(:)
  
  CALL RANDOM_SEED(SIZE=seedSize)
  ALLOCATE(randomSeed(seedSize), stat=errorFlag)
  CALL SYSTEM_CLOCK(COUNT=clock)
  randomSeed = clock + 37*(/ (i-1, i=1,seedSize) /)
  CALL RANDOM_SEED(PUT=randomSeed)

  DEALLOCATE(randomSeed)
  RETURN
END SUBROUTINE InitRandomSeed


SUBROUTINE InitVelocity(temperature,v)
!----------------------------------------------------------------------------
! DESCRIPTION:
!   This function initializes the velocities of atom before 
!   doing MD simulation.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  USE SYS, ONLY: frozenIon
  USE MPI_Functions, ONLY: rankGlobal
  USE MPI_Functions, ONLY: BcastReal_Dim2

  IMPLICIT NONE
  
  !! >> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN)  :: temperature ! temperature to be simulated
  REAL(KIND=DP), INTENT(OUT) :: v(3,cell%numIon) ! velocity of each atom, in a.u.

  !! >> LOCAL VARIABLES << !!
  INTEGER  :: i, ii, nfrozen
  REAL(KIND=DP) :: & 
     rand(2) ,         &
     vel1, vel2, mass, & 
     KE                           ! to store kinetic energy

  CALL Title("MolecularDynamics::InitVelocity")

  ! initialize here as it otherwise causes a runtime error with check enabled on
  ! ifort (JMD, 2013)
  nfrozen = 0

  ! note from mohan 2013-03
  ! get the random wave functions from processor 1!!
  ! otherwise the velocity on each node will be different!
  ! then the geometry will be different at each MD step
  ! the results are totally wrong.
  IF(rankGlobal == 0) THEN

    ! Assign random numbers to velocities
    DO i = 1, cell%numIon*3/2 + 1
      CALL RANDOM_NUMBER(rand)
      vel1 = SQRT(-2._DP*LOG(rand(1)))*COS(2*pi*rand(2))
      vel2 = SQRT(-2._DP*LOG(rand(1)))*SIN(2*pi*rand(2))
      IF ((2*i-2)/3 < cell%numIon) v(MOD(2*i-2,3)+1,(2*i-2)/3+1) = vel1
      IF ((2*i-1)/3 < cell%numIon) v(MOD(2*i-1,3)+1,(2*i-1)/3+1) = vel2
    END DO

    ! Scale velocity according to mass of atoms
    DO i = 1,cell%numIon
      mass = cell%elementTable(cell%ionTable(i)%elementID)%mass
      v(:,i) = v(:,i)/mass
    END DO
  
    CALL RemoveMovementOfCenterOfMass(v, .TRUE.)
  
    ! Do not move atoms which are frozen
    nfrozen = 0
    IF (ALLOCATED(frozenIon)) THEN 
      DO ii=1,cell%numIon
        IF (frozenIon(ii,1) .EQV. .TRUE.) v(1,ii) = 0.d0
        IF (frozenIon(ii,2) .EQV. .TRUE.) v(2,ii) = 0.d0
        IF (frozenIon(ii,3) .EQV. .TRUE.) v(3,ii) = 0.d0
        IF (frozenIon(ii,1) .EQV. .TRUE.) nfrozen = nfrozen + 1
      ENDDO
    ENDIF

    ! Rescale all velocities to start to exact temperature desired
    CALL  GetAtomKE(v,KE)
    WRITE(outputUnit,*) "Kinetic energy from random distribution of velocities = ", KE
    if(abs(KE) < 1.d-20) then
      write(*,*) "ERROR: Kinetic energy in init is almost zero ", KE
      write(*,*) "       Stopping calculation as MD is pointless"
      flush(6)
      stop
    endif

    v = v * SQRT(temperature*1.5_dp * REAL( cell%numIon-nfrozen,kind=DP) / KE)

  ENDIF

  ! mohan fix bug 2013-03
  ! distribute the random velocity from first processor to other
  ! processors
  CALL BcastReal_Dim2( v , cell%numIon )

  ! mohan add 2013-04-17
  ! mohan fix 2014-05-02
  CALL  GetAtomKE(v,KE)
  WRITE(outputUnit,*) "In order to start from the input temperature, kinetic energy should be = ",KE

!  WRITE(outputUnit,*) " init velocity:"
!  DO i=1, cell%numIon
!    WRITE(outputUnit,*) v(1,i), v(2,i), v(3,i)
!  ENDDO

!  ENDIF
    
  RETURN
END SUBROUTINE InitVelocity


SUBROUTINE RemoveMovementOfCenterOfMass(vel, adjust)
!---------------------------------------------------------------------------
! DESCRIPTION:
!    Compute the velocity for the center of mass, and remove it
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
  !! >> EXTERNAL VARIABLES << !!
  REAL(KIND=DP), INTENT(INOUT) :: vel(3,cell%numIon) ! velocties of each atom.

  !! >> LOCAL VARIABLES << !!
  INTEGER  :: ii, id
  REAL(KIND=DP) :: avgvel(3)
  REAL(KIND=DP) :: avel(3)
  REAL(KIND=DP) :: mass
  REAL(KIND=DP) :: totMass
  LOGICAL :: adjust

  !! >> FUNCTION << !!
  avgvel(:) = 0.d0
  totMass = 0.d0

  !---------------------------------------------------------
  ! sum up the total momentum of all the atoms
  ! 'i' is the index of atom
  ! in the system. (Px,Py,Pz) =\sum m_{i}*(vx,vy,vz)_{i}
  !---------------------------------------------------------
  DO ii=1,cell%numion
    id = cell%ionTable(ii)%elementID
    mass = cell%elementTable(id)%mass
    avgvel = avgvel + mass * vel(:,ii)
    totMass = totMass + mass
  ENDDO

  ! unit of 'avgvel' is Bohr/a.u.
  avgvel = avgvel / totMass

  avel = avgvel * bohr ! change distance unit form Bohr to Angstrom
  avel = avel / dt * 1.0e-15 ! change time unit from a.u. to picosecond

  WRITE(outputUnit,'(A,3F12.6,A)') " Velocity of Mass Center  : " &
  , avel(1), avel(2), avel(3), " (Angstrom/fs)"

  ! mohan add 2014-05-02
  IF(adjust) THEN
    WRITE(outputUnit,*) "Conservation of total momentum, substrate velocity: "
    WRITE(outputUnit,'(A,3F10.6)') " Average velocity : ", avgvel(1), avgvel(2), avgvel(3)
    DO ii=1,cell%numIon
      vel(:,ii) = vel(:,ii) - avgvel(:)
    ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE RemoveMovementOfCenterOfMass


SUBROUTINE MonitorMeanSquareDisplacement(step,msd,diffuCoeff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!  Calculate the Mean Square Displcacement of particles
!  Compute the mean square displacement of atoms according to formula:
!    MSD = <|r(t)-r(0)|^2>
!  <...> denotes the averaging over all the atoms
!
!  For solids, MSD saturates to a finte value, for liquids, 
!  MSD grows linearly with time.
! 
!  Also compute the diffusion coefficient D:  
!   D = MSD/(6*t)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  USE Output, ONLY:  WrtOut
  USE OutputFiles, ONLY: outputUnit
  USE Constants, ONLY: fundamentalTime
      
  IMPLICIT NONE
      
  ! >> EXTERNAL VARIABLES << !
  INTEGER :: step
  REAL(KIND=DP), INTENT(OUT) :: msd, diffuCoeff

  ! >> LOCAL VARIABLES << !
  REAL(KIND=DP) :: box(3,3)
  INTEGER :: ii 
  REAL(KIND=DP) :: timePassed, diff(3), tmp

  ! INITIALIZATION
  box = cell%cellReal

  ! mohan add 2013-07-11
  IF(msd_startTime>ABS(rstMD)) THEN
    msd_startTime=msd_startTime
  ELSE IF(rstMD<0) THEN
    msd_startTime = ABS(rstMD)+1
  ELSE
    msd_startTime = 1
  ENDIF


  ! >> FUNCTION << !
    
  IF (step == msd_startTime) THEN
    ALLOCATE(msd_coords0(3,cell%numIon))
    ALLOCATE(msd_coordst(3,cell%numIon))
    DO ii=1,cell%numIon
      msd_coords0(:,ii) = cartNoWrap(:,ii)
    ENDDO
  ENDIF

  IF (step > msd_startTime) THEN
    DO ii=1,cell%numIon
      msd_coordst(:,ii) = cartNoWrap(:,ii)
    ENDDO
    !-----------------------------------
    ! compute msd and diffuCoeff  
    msd = 0._DP

    IF(doMSDAtom .EQV. .TRUE.) THEN
      WRITE(outputUnit,*) "MSD for each atom (Bohr^2):"
    ENDIF 

    DO ii=1,cell%numIon  ! loop over ions
      diff(:) = msd_coords0(:,ii) - msd_coordst(:,ii)
      ! the mean square displacement has unit of Bohr^2,
      ! so we should not use 'SQRT' here.
      ! mohan note 2013-05-01
      !tmp = SQRT( SUM(diff**2) )
      tmp = SUM(diff**2)
      msd = msd + tmp

      ! mohan add 2013-04-22
      IF(doMSDAtom .EQV. .TRUE.) THEN
        WRITE(outputUnit,'(1X,I8,1X,ES14.6)') ii, tmp
      ENDIF
    ENDDO

    msd = msd / REAL(cell%numIon,KIND=DP) * BOHR * BOHR
    timePassed = dt * ( step - msd_startTime )
    diffuCoeff = msd / (6._DP*timePassed) / (fundamentalTime / 1.0e-12)

    WRITE(outputUnit,'(A,F10.4,A)') " Mean Square Displacement : ", msd, " Angstrom^2"
    WRITE(outputUnit,'(A,F10.4,A)') " Diffusion Coefficient    : ", diffuCoeff, " Angstrom^2/ps"

!--------------------------------------------------
! NOTE !
! to convert the diffuCoeff to A^2/picosecond, 
! D * bohr * bohr / ( fundamentalTime / 1.0e-12 )
! where fundamentalTime is 2.418884326505E-17_DP
!--------------------------------------------------
!!DEBUG          
!          DO ii = 1,cell%numIon
!            WRITE(message,'(a,I6,a,ES12.4)') & 
!              ' ion:',ii,' ion_sd:', & 
!              SQRT(SUM((msd_coordst(:,ii)-msd_coords0(:,ii))**2))
!            CALL WRTOUT(6,Message)
!            WRITE(message,'(a,I6,a,3ES12.4)') & 
!              ' ion:',ii,' ion_fracCoord:', cell%ionTable(ii)%coord(:)
!              CALL WRTOUT(6,Message)
!          ENDDO ! loop over ions
!!~DEBUG

    ENDIF ! if (step > msd_startTime)


    RETURN

END SUBROUTINE MonitorMeanSquareDisplacement


!SUBROUTINE MonitorRDF(step,cell,m,s,g)
!!--------------------------------------------------------------------------
!! DESCRIPTION:
!!  Compute radial distribution function (g), 
!!  and structure factor (s)
!!----------------------------------------------------------------------------
!  use DataTypes, ONLY: & 
!    cellStruct
!  implicit none
!  integer,intent(in) :: &  
!    step, & 
!    m
!  type(cellStruct), intent(in) :: & 
!    cell
!  real(dp),intent(out) :: & 
!    s(m), &   ! Structure factor
!    g(m)      ! Radial distribution function
!  
!  
!END SUBROUTINE MonitorRDF


SUBROUTINE OutputIonVelocityAndStates(iter, vel, md_type, xLogS, vLogS, vBoxG)
!------------------------------------------------------------------------------
! DESCRIPTION:
!  Output the 
!   (1) ion's velocity 
!   (2) thermostat's xLogS, and vLogs (for NVT, NPT)
!   (3) barostat's vbox (box's velocity)  (NPT only)
!  for restart purpose
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!    Dec/29/2008: Created by Chen Huang
!------------------------------------------------------------------------------
  USE Output, ONLY : Wrtout

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: iter, md_type !number of ions, iteration number, NVT/NPT
  REAL(DP), INTENT(IN) :: vel(3,cell%numIon) 
  REAL(KIND=DP), INTENT(IN) :: xLogS
  REAL(KIND=DP), INTENT(IN) :: vLogS
  REAL(DP), INTENT(IN), OPTIONAL :: vBoxG(3,3)
  
  ! >>>>>>>>>>>> local vars <<<<<<<<<<<!
  LOGICAL :: pass
  INTEGER :: i
  CHARACTER(LEN=100) :: iteration_number
  CHARACTER(LEN=500) :: message

  pass = .FALSE.

  ! only rankGlobal=0 will output the velocity.
  IF (rankGlobal/=0) RETURN


  IF (dump_md_freq==1 .OR. & 
      iter==1 .OR. & 
      ( dump_md_freq > 1 .and. MOD(iter,dump_md_freq)==0 ) ) & 
    pass = .TRUE.
  IF (.NOT. pass) RETURN


  WRITE(iteration_number,'(I8)') iter
  WRITE(message,'(a,a,a,a)') TRIM(md_output_path),'/vel.',TRIM(ADJUSTL(iteration_number)),'.dat'
  ! open the file: 'message'
  OPEN (unit=1988, file=message, action='write')
  WRITE(message,'(a,I8)') 'MD STEP=',iter
  CALL Wrtout(1988,message)

  !----- thermo-stat's velocity and positions
  WRITE(message,'(ES20.12, a)') xLogS, " <= thermostat position : xLogS (a.u.)"; CALL WrtOut(1988,message)
  WRITE(message,'(ES20.12, a)') vLogS, " <= thermostat velocity : vLogS (a.u.)"; CALL WrtOut(1988,message)
  SELECT CASE(md_type)
  CASE (1) ! NVT
    WRITE(message,'(A)') " 0.0 0.0 0.0 <= barostat velocity 1"; CALL WrtOut(1988,message)
    WRITE(message,'(A)') " 0.0 0.0 0.0 <= barostat velocity 2"; CALL WrtOut(1988,message)
    WRITE(message,'(A)') " 0.0 0.0 0.0 <= barostat velocity 3"; CALL WrtOut(1988,message)
    ! add addiontional outputs if needed
  CASE (2) ! NPT
    WRITE(message,'(3ES20.12, a)') vBoxG(1,:), " <= barostat velocity : vBoxG(1,:)"; CALL WrtOut(1988,message)
    WRITE(message,'(3ES20.12, a)') vBoxG(2,:), " <= barostat velocity : vBoxG(2,:)"; CALL WrtOut(1988,message)
    WRITE(message,'(3ES20.12, a)') vBoxG(3,:), " <= barostat velocity : vBoxG(3,:)"; CALL WrtOut(1988,message)
  CASE (3) ! NVE
    WRITE(message,'(A)') " 0.0 0.0 0.0 <= barostat velocity 1"; CALL WrtOut(1988,message)
    WRITE(message,'(A)') " 0.0 0.0 0.0 <= barostat velocity 2"; CALL WrtOut(1988,message)
    WRITE(message,'(A)') " 0.0 0.0 0.0 <= barostat velocity 3"; CALL WrtOut(1988,message)
  END SELECT
  
  ! output velocities
  CALL WrtOut(1988, "ION VELOCITIES (a.u.):")
  DO i = 1,cell%numIon
    !WRITE(message,'(3(a,Es12.4))') & 
    WRITE(message,'(3(a,Es16.8))') & 
      "  ", vel(1,i),"  ",vel(2,i),"  ",vel(3,i)
    CALL WrtOut(1988,message)
  ENDDO



  CLOSE(1988)
  
END SUBROUTINE OutputIonVelocityAndStates


SUBROUTINE OutputMDGeom(iter)
!---------------------------------------------------------------------------
! DESCRIPTION:
!  Output the ion's position and cell's structure for monitoring and MD
!  restart.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!---------------------------------------------------------------------------
! REVISION LOG:
!    Dec/29/2008: Created by Chen Huang
!----------------------------------------------------------------------------
  USE Output, ONLY : Wrtout
  USE Constants, ONLY : bohr

  IMPLICIT NONE
  LOGICAL :: pass
  INTEGER  :: i,iter,uf
  CHARACTER(LEN=100) :: iteration_number
  CHARACTER(LEN=500) :: message, eleName
   
  ! >> FUNCTION << !
  pass = .FALSE.
  IF (rankGlobal/=0) RETURN
  IF (dump_md_freq==1 .OR. & 
      iter==1 .OR. & 
      ( dump_md_freq > 1 .AND. MOD(iter,dump_md_freq)==0 ) ) & 
    pass = .TRUE.
  IF (.not. pass) RETURN

  uf = 1999

  ! the final file name is the md_geometry plus the iteraction number
  WRITE(iteration_number,'(I8)') iter
  WRITE(message,'(a,a,a,a)') TRIM(md_output_path),'/ion.',TRIM(ADJUSTL(iteration_number)),'.dat'

  !-------------------------------
  ! open the file: 'ion.1.dat'
  !-------------------------------
  OPEN (UNIT=uf, FILE=message, ACTION='write')
  ! translate the ionic coordinates from Angstrom to Bohr
  WRITE(message,'(a)') '%BLOCK LATTICE_CART'
  CALL WrtOut(uf,message)
  WRITE(message,'(3(a,Es20.12))') '  ',cell%cellReal(1,1)*bohr,'  ',cell%cellReal(2,1)*bohr,'  ',cell%cellReal(3,1)*bohr
  CALL WrtOut(uf,message)
  WRITE(message,'(3(a,Es20.12))') '  ',cell%cellReal(1,2)*bohr,'  ',cell%cellReal(2,2)*bohr,'  ',cell%cellReal(3,2)*bohr
  CALL WrtOut(uf,message)
  WRITE(message,'(3(a,Es20.12))') '  ',cell%cellReal(1,3)*bohr,'  ',cell%cellReal(2,3)*bohr,'  ',cell%cellReal(3,3)*bohr
  CALL WrtOut(uf,message)
  WRITE(message,'(a)') '%ENDBLOCK LATTICE_CART'
  CALL WrtOut(uf,message)
  
  !-------------------------------
  ! Print the atomic coordinates
  !-------------------------------
  WRITE(message,'(a)') '%BLOCK POSITIONS_CART'
  CALL WrtOut(uf,message)
  DO i = 1, cell%numIon
    WRITE(eleName ,'(a)') & 
      cell%elementTable(cell%ionTable(i)%elementID)%elementName
    ! translate the ionic coordinates from Angstrom to Bohr
    WRITE(message,'(a,a,3(a,Es16.8))') & 
      "  ",TRIM(ADJUSTL(eleName)), & 
      "  ",cartNoWrap(1,i)*bohr, & ! unit of cartNorap is 'Bohr' 
      "  ",cartNoWrap(2,i)*bohr, &  
      "  ",cartNoWrap(3,i)*bohr
    CALL WrtOut(uf,message)
  ENDDO
  WRITE(message,'(A)') '%END BLOCK POSITIONS_CART'
  CALL WrtOut(uf,message)

  !-------------------------------
  ! Print the local potential.
  !-------------------------------
  WRITE(message,'(A)') '%BLOCK SPECIES_POT'
  CALL WrtOut(uf, message)
  DO i = 1, cell%numIonType
    WRITE(eleName, '(A)') cell%elementTable(i)%elementName
    WRITE(message, '(A,A,A)') TRIM(eleName), ' ', cell%elementTable(i)%psp%pspFileName
   CALL WrtOut(uf, message)
  ENDDO
  WRITE(message,'(A)') '%END BLOCK SPECIES_POT'
  CALL WrtOut(uf, message)

  CLOSE(uf)

  RETURN
  
END SUBROUTINE OutputMDGeom



SUBROUTINE ReadNewTemp( temp, step )
!----------------------------------------------------------------------------
! DESCRIPTION:
!  If fixTemperature == 0, then this subroutine will be skipped
!  otherwise, we are going to read a new temperature from the disk.
!  You only need to create a file named temperature.dat, and then 
!  change the temperature in it. 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------
 
  USE Constants, ONLY: BOLTZMANN

  IMPLICIT NONE

  INTEGER :: step, fs
  REAL(DP) :: temp, in_temp
  CHARACTER (LEN=500) :: message

  IF ( fixTemperature > 0 ) THEN

    ! Change the temperature every 'fixTemperature' steps.
    IF ( MOD(step, fixTemperature) == 0 ) THEN
    
      ! Read in new temperature from file.
      OPEN(unit=111, file='ChangeTemp.dat', action='READ', iostat=fs)
      IF (fs/=0) THEN
        CALL WRTOUT(6,"ERROR IN OPENING ChangeTemp.dat, CODE STOP!")
        STOP
      ENDIF
      READ(111,*) in_temp
      CLOSE(111)

      ! Renew information.
      in_temp =  in_temp * BOLTZMANN
      IF ( ABS(in_temp-temp) > 1e-6 ) THEN
        WRITE(message,'(a,F12.4,a,F12.4)') & 
          "(ReadNewTemp): Read in new temp:", in_temp/BOLTZMANN, & 
          " previous temp:", temp/BOLTZMANN
        CALL WrtOut(6,message)
        temp = in_temp
      ELSE
        WRITE(message,'(a,F12.4,a,F12.4,a)') & 
          "(ReadNewTemp): new temp:", in_temp/BOLTZMANN & 
          , " previous temp:", temp/BOLTZMANN, ". No change of temp."
        CALL WrtOut(6,message)
      ENDIF

    ENDIF
  ENDIF

  RETURN 
END SUBROUTINE ReadNewTemp


SUBROUTINE AdjustCenterOfMass(shiftCM)
!---------------------------------------------------------------------------
! DESCRIPTION:
!  Compute the center of mass of the atoms and adjust the center of mass to 
!  the origin (0,0,0).
!  On exit, cell%ionTable%coord will be updated (if shiftCM=true), to have 
!  the center of mass at (0,0,0)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------
  USE MathFunctions, ONLY: Inverse

  IMPLICIT NONE
  REAL(KIND=DP) :: cx,cy,cz
  LOGICAL, OPTIONAL :: shiftCM
  INTEGER :: ii  
  REAL(KIND=DP) :: mass 
  REAL(KIND=DP) :: totMass 
  REAL(KIND=DP) :: invCellReal(3,3)

  ! >>> FUNCTION BEGINS <<<!

  cx = 0.d0
  cy = 0.d0
  cz = 0.d0
  totMass = 0.d0
  DO ii=1, cell%numIon
    mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
    totMass = totMass + mass
    cx = cx + cartNoWrap(1,ii) * mass
    cy = cy + cartNoWrap(2,ii) * mass
    cz = cz + cartNoWrap(3,ii) * mass
  ENDDO
  cx = cx / totMass
  cy = cy / totMass
  cz = cz / totMass
  
  ! Adjust the fractional coords if shiftCM is TRUE
  ! Notice: we are in fact tracking cartNoWrap not cell%ionTable%coord
  ! the later one is only used for OFDFT computations
  IF (PRESENT(shiftCM) .AND. shiftCM .EQV. .true.) THEN
     invCellReal = Inverse(cell%cellReal)
     DO ii=1, cell%numIon 
       cartNoWrap(1,ii) = cartNoWrap(1,ii) - cx
       cartNoWrap(2,ii) = cartNoWrap(2,ii) - cy
       cartNoWrap(3,ii) = cartNoWrap(3,ii) - cz
       cell%ionTable(ii)%coord = MODULO(MATMUL(invCellReal,cartNoWrap(:,ii)), 1._DP)
     ENDDO
  ENDIF

  WRITE(outputUnit,'(A,3F10.3,A)') " Center of Mass           :" &
  , cx*bohr, cy*bohr, cz*bohr &
  , " (Angstrom * a.u.)"

  RETURN
END SUBROUTINE AdjustCenterOfMass


SUBROUTINE CalFreedom()
!---------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the total freedom number in the system.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
! Created by Mohan 2014-07-21.
!----------------------------------------------------------------------------

  USE SYS, ONLY: frozenIon

  IMPLICIT NONE

  !! >> INTERNAL VARIABLES << !!

  INTEGER :: nfrozen
  ! number of frozen freedom
  !
  INTEGER :: ii
  ! index
  !

  !! >> INITIALIZATION << !!
  freedom = 0
  nfrozen = 0

  !! >> FUNCTION << !!

  ! Calculate # of frozen atoms 
  nfrozen = 0
  IF (ALLOCATED(frozenIon)) THEN 
    DO ii=1,cell%numIon
      IF (frozenIon(ii,1) .EQV. .TRUE.) nfrozen = nfrozen + 1
      IF (frozenIon(ii,2) .EQV. .TRUE.) nfrozen = nfrozen + 1
      IF (frozenIon(ii,3) .EQV. .TRUE.) nfrozen = nfrozen + 1
    ENDDO
  ENDIF

  freedom = 3*REAL(cell%numIon,KIND=DP) - REAL(nfrozen, KIND=DP)

  WRITE(outputUnit,'(A,F10.2)') " Freedom             : ",freedom
  WRITE(outputUnit,'(A,I10)')   " Frozen  freedom     : ",nfrozen
  
  RETURN

END SUBROUTINE CalFreedom


END MODULE MolecularDynamics

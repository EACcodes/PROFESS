MODULE NVE
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE MolecularDynamics
!     |_SUBROUTINE RunNVE
!     |_SUBROUTINE RescaleVelocity
!     |_FUNCTION Conserved
!
! DESCRIPTIONS:
!   This is used to carry out molecular dynamics for NVE ensemble.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   Martyna et. al., "Explicit reversible integrators for extended system
!      dynamics," Mol. Phys. 87 (1996), 1117.
!   Allen and Tildesley, "Computer Simulation of Liquids" (predictor-corrector)
!   Frenkel and Smit, "Understanding Molecular Simulation".
!------------------------------------------------------------------------------
! REVISION LOG:
! 2013 Created by Mohan Chen
!------------------------------------------------------------------------------

  USE MolecularDynamics

  IMPLICIT NONE

CONTAINS


SUBROUTINE RunNVE(Optimizer, energy, forces)
!------------------------------------------------------------------------------
! DESCRIPTION: 
!   The driver routine for NVE ensemble OFDFT molecular dynamics.
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

  USE Report, ONLY: GeometryMinimizerReportHeader
  USE Report, ONLY: GeometryMinimizerReportSteps
  USE Report, ONLY: GeometryMinimizerReportFooter

  USE CalForces, ONLY : CalculateForces

  USE RefreshCell, ONLY: RefreshCellSetup
  USE RefreshIons, ONLY: RescaleDensity
  USE RefreshIons, ONLY: RefreshIonTerms

  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!


  EXTERNAL Optimizer 
  ! A subroutine that is called to optimize the electron 
  ! density relative to the ion positions
  !
  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy
  ! Energy from electronic optimization
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: forces
  ! Total forces.  Final index is 1 for total force,
  ! First is ion number, second direction (1,2,3 for x,y,z)
  !

                        !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    ii, &
    i, &                       ! Counters to parse through the ion table.
    step, &                    ! step number we are on
    numIon, &                  ! total number of atoms
    nStep, &                   ! Total number of steps
    errorFlag

  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: velocity
  ! Velocity of ions (atomic units)
  !
  REAL(KIND=DP), DIMENSION(3) :: fracStep
  ! Step (due to velocity) in fractional coordinates
  !
  REAL(KIND=DP) :: &
    mass,      &   ! atom mass, temperay var
    maxStep, &     ! Largest movement by an ion
    xLogs, &
    vLogs, &
    hamiltonian, &
    conservedE, &
    maxForce, &
    msd, &                  !  mean square displacement
    diffuCoeff, &           ! diffusion coefficient
    twiceKE, &     ! Kinetic energy x 2
    oldEtot        ! old energy

  CHARACTER(len=80) :: parameterInfo
  ! Line containing values for Qmass, dt, temperature for output to header
  !
  REAL(KIND=DP) :: tempNow
  !
  INTEGER :: md_type = 3  
  ! NVE = 3,


                        !>> INITIALIZATION <<!

  CALL Title("MolecularDynamics::NVE")
  CALL WrtOut(6," (NVE) Start.")

  watch2 = TimerStart()
  watch = TimerStart()
  numIon = SIZE(cell%ionTable)
  
  ! Calculate freedom of the system 
  CALL CalFreedom()

  ! Calculate the number of total MD steps.
  nStep = CEILING(timeTot/dt)
  forces = 0._DP

  ! Allocate arrays for velocities, cartesian coordinates.
  ALLOCATE(velocity(3,numIon), stat=errorFlag)
  ALLOCATE(cartNoWrap(3,numIon), stat=errorFlag)
  
  ! Set up extra output to ion optimizer / MD header
  WRITE(outputUnit,'(a,ES20.12,a)') "Time interval       : ", dt*fundamentalTime*10._DP**15, " (fs)"
  WRITE(outputUnit,'(a,ES20.12,a)') "Target temperature  : ", temperature/boltzmann, " (K)"
  WRITE(outputUnit,'(a,I8)') "Number of NVE steps : ", nStep

  IF (rankGlobal==0) CALL GeometryMinimizerReportHeader(parameterInfo)

  IF ( rstMD == 0) THEN 
    CALL InitRandomSeed
    CALL InitVelocity(temperature,velocity)
    vLogs = 0._DP
    xLogs = 0._DP
  ELSE IF ( (rstMD == 1) .OR. (rstMD<0) ) THEN
    CALL RestartMD(md_type,velocity,xLogS,vLogS)
    CALL RefreshCellSetup
    CALL RescaleDensity(rhoR)
  ENDIF
  CALL GetAtomKE(velocity,twiceKE)
  twiceKE = twiceKE * 2._DP

  tempNow = twiceKE/freedom/boltzmann

  IF (rankGlobal==0) WRITE(*,*) " start temperature = ", tempNow


  ! mohan add 2013-07-11
  IF(rstMD<0) THEN 
    startStep=ABS(rstMD)
  ELSE
    startStep=0
  ENDIF


  CALL Optimizer                        
  ! Optimize the density

  CALL CalculateForces(rhoR,forces)
  ! calculate the ground state forces

  IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP
  maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))

  conservedE = Conserved(twiceKE/2._DP, energy(1))
  oldEtot=energy(1)

  WRITE(outputUnit,'(A)') "Kin/Ext/Col/Exc/Ewald/TF/vW/nonlocalKEDF/"
  WRITE(outputUnit,'(ES20.10, ES20.10, ES20.10)') energy(2), energy(3), energy(4)
  WRITE(outputUnit,'(ES20.10, ES20.10, ES20.10)') energy(5), energy(6), energy(7)
  WRITE(outputUnit,'(ES20.10, ES20.10)') energy(8), energy(9)

  ! Print intial output
  CALL GeometryMinimizerReportSteps(startStep, TimerStop(watch), forces, maxForce, &
                                    0._DP, .TRUE., hamiltonian, .TRUE.)

  DO i=1,numIon
    cartNoWrap(:,i) = MATMUL(cell%cellReal,cell%ionTable(i)%coord)
  ENDDO


  WRITE(message,'(1X,A,7X,A,11X,A,9X,A,4X,A)') "NVE_STEP ", "SystemEnergy", "Conserved", "DeltaE", "Temperature"
  CALL WrtOut(6,message)

  ! Output the message to the screen.
  WRITE(message,'(I8,ES20.10, ES20.10, ES15.5, ES15.5)') &
      startStep, energy(1), conservedE, energy(1)-oldEtot, twiceKE/freedom/boltzmann
  CALL WrtOut(6,message)
  oldEtot=energy(1)


                        !>> FUNCTION BODY <<!

  ! Loop for MD step 
  DO step = startStep+1, nStep

    WRITE(outputUnit,'(''------------------------------------------------------------------------------'')')
    WRITE(outputUnit,'(A,I6)') "MD(NVE) STEP ", step

    ! mohan add 2014-12-01
    CALL RemoveMovementOfCenterOfMass(velocity, .TRUE.)
    CALL AdjustCenterOfMass(.TRUE.) 

    ! Calculate the Mean-Square-Displacement.
    CALL MonitorMeanSquareDisplacement(step,msd,diffuCoeff)
  
    ! (1) 1st step of Verlet-Velocity
    ! New velocity are obtained
    DO ii=1,size(cell%ionTable,1)        
      mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
      velocity(:,ii) = velocity(:,ii) + forces(ii,:,1)/mass*dt/2._DP
    ENDDO
    
    ! (2) 2nd step of Verlet-Velocity 
    ! Update the Non-Wrapped cartesion coordinate
    DO i=1,numIon
      cartNoWrap(:,i) = velocity(:,i)*dt + cartNoWrap(:,i)
    ENDDO
    
    ! Calculate the maximal velocities.
    ! The step in fractional coordinates
    maxStep = 0._DP
    DO i = 1, numIon
      maxStep = MAX(maxStep, SUM(velocity(:,i)**2))
      fracStep = Vecmul(Inverse(cell%cellReal), velocity(:,i)*dt)
      cell%ionTable(i)%coord = MODULO(cell%ionTable(i)%coord + fracStep, 1._DP)
    END DO
    maxStep = SQRT(maxStep)*dt

    ! (3.1) OFDFT optimizer: in order to calculate the force
    ! Refresh the ionpositions, including the Ewald term
    ! and the local pseudopotentials term.
    CALL RefreshIonTerms
    ! This must be called every time the ion positions are changed

    ! (3.2)
    ! Get the new charge density in the new ion positions.
    CALL Optimizer                        ! Optimize the density

    ! (3.3) Calculate the forces. 
    CALL CalculateForces(rhoR, forces)
    IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP
    maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))

    ! (3.4) 3rd step to Update velocity
    DO ii=1,size(cell%ionTable,1)        
      mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
      velocity(:,ii) = velocity(:,ii) + forces(ii,:,1)/mass*dt/2._DP
    ENDDO

    ! calculate the kinetic energy of ions
    CALL GetAtomKE(velocity,twiceKE)
    twiceKE = 2._DP * twiceKE

    !------------------------------------------------------
    ! this part should be not necessary in NVE
    ! mohan 2013-05-09
    ! calculate the conserved quantity during MD 
    !hamiltonian = Conserved(twiceKE/2._DP, energy(1))
    !newKE = 0.5*twiceKE-(hamiltonian-conservedE) 
    !newKE = 2.0*newKE;

    ! updates velocity in order to adjust temperature
    !CALL RescaleVelocity(twiceKE, newKE, numIon)
    !------------------------------------------------------

    ! calculate the conserved quantity during MD 
    hamiltonian = Conserved(twiceKE/2._DP, energy(1))

! kinetic, external, coulombic, exchange-correlation,
! ion-ion, Thomas-Fermi, von Weiszacker and
! third term Wang-Teter, WGC,

    WRITE(outputUnit,'(A)') "Kin/Ext/Col/Exc/Ewald/TF/vW/nonlocalKEDF/"
    WRITE(outputUnit,'(ES20.10, ES20.10, ES20.10)') energy(2), energy(3), energy(4)
    WRITE(outputUnit,'(ES20.10, ES20.10, ES20.10)') energy(5), energy(6), energy(7)
    WRITE(outputUnit,'(ES20.10, ES20.10)') energy(8), energy(9)

    ! For monitoring MD or restarting MD
    CALL OutputMDGeom(step)
    CALL OutputIonVelocityAndStates(step,velocity, md_type,xLogS,vLogS)

    ! Print output
    CALL GeometryMinimizerReportSteps(step, TimerStop(watch), &
                                      forces, maxForce, maxStep, &
                                      .TRUE., hamiltonian, .TRUE.)
                                      

    watch = TimerStart()
    
    ! Output the message to the screen.
    WRITE(message,'(I8,ES20.10, ES20.10, ES15.5, ES15.5)') &
      step, energy(1), hamiltonian, energy(1)-oldEtot, twiceKE/freedom/boltzmann
    CALL WrtOut(6,message)
    oldEtot=energy(1)

  END DO
  !----------------------------------------
  ! End of NVE MD loop
  !----------------------------------------

  CALL GeometryMinimizerReportFooter(TimerStop(watch2))

  CALL WrtOut(6,"(NVE): finished.")
  
  DEALLOCATE(cartNoWrap)
  DEALLOCATE(velocity)
  IF (ALLOCATED(msd_coords0)) DEALLOCATE(msd_coords0)
  IF (ALLOCATED(msd_coordst)) DEALLOCATE(msd_coordst)

  RETURN

CONTAINS


SUBROUTINE RescaleVelocity(KE, newKE, numIon)
!---------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2013 Created by Mohan Chen
!----------------------------------------------------------------------------
  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN) :: KE, newKE
  INTEGER, INTENT(IN) :: numIon

                        !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: scaling
  ! scaling apply to velocities
  !
  REAL(KIND=DP) :: tempNow
  ! temperature for now 
  !
  REAL(KIND=DP) :: tempTarget
  ! target temperature
  !


                        !>> INITIALIZATION <<!

  tempNow = KE/freedom/boltzmann
  tempTarget = newKE/freedom/boltzmann

  WRITE(outputUnit,'(a,ES20.12,a)') "Temperature now     : ", tempNow, " (K)"
  WRITE(outputUnit,'(a,ES20.12,a)') "New temperature  : ", temperature/boltzmann, " (K)"

  scaling = SQRT(tempTarget/tempNow)
                        !>> FUNCTION BODY <<!

  WRITE(outputUnit,*) " scaling factor of velocity : ", scaling
  velocity = velocity*scaling
  twiceKE = twiceKE*scaling**2

  RETURN

END SUBROUTINE RescaleVelocity


FUNCTION Conserved(KE, PE)
!---------------------------------------------------------------------------
! DESCRIPTION:
! This function calculates the conserved quantity for the NVE system. 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2013 Created by Mohan Chen
!----------------------------------------------------------------------------

  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN) :: KE
  ! Kinetic energy of particles (due to their velocity)
  !
  REAL(KIND=DP), INTENT(IN) :: PE
  ! Potential energy of particles (energy calculated using OFDFT)
  !
  REAL(KIND=DP) :: Conserved 
  ! The conserved quantity
  !

                        !>> INITIALIZATION <<!

                        !>> FUNCTION BODY <<!

   Conserved = KE + PE 
                    
   WRITE(outputUnit,'(a,ES25.12,a)') "NVE Conservation     : ", Conserved," (Hartree)"
   WRITE(outputUnit,'(a,ES25.12,a)') "NVE Temperature      : ", twiceKE/freedom/boltzmann," (K)"
   WRITE(outputUnit,'(a,ES25.12,a)') "NVE Kinetic energy   : ", KE," (Hartree)"
   WRITE(outputUnit,'(a,ES25.12,a)') "NVE Potential energy : ", PE," (Hartree)"

   RETURN

END FUNCTION Conserved


END SUBROUTINE RunNVE


END MODULE NVE

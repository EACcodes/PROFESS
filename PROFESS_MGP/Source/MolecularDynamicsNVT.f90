MODULE NVT
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE NVT 
!     |_SUBROUTINE RunNVT
!     |_SUBROUTINE NHIntegrator
!     |_SUBROUTINE NHhamiltonian
!
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!------------------------------------------------------------------------------
! REVISION LOG:
! First written by Linda Hung
! 2012-2013 Reconstructed by Mohan
!------------------------------------------------------------------------------

  USE MolecularDynamics

  IMPLICIT NONE

CONTAINS


SUBROUTINE RunNVT(Optimizer, energy, forces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine integrates motion for molecular dynamics using the
!   Nose-Hoover thermostat coupled with a predictor-corrector integrator.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!     10/20/2008 : Adapted from subroutines in IonOptimizers and symplectic
!                  algorithm from Martyna et. al. with help from Frenkel/Smit.
!
!------------------------------------------------------------------------------

  USE SYS, ONLY: frozenIon
  USE OUTPUT, ONLY: printStress

  USE REPORT, ONLY: GeometryMinimizerReportHeader 
  USE REPORT, ONLY: GeometryMinimizerReportSteps
  USE REPORT, ONLY: GeometryMinimizerReportFooter

  USE RefreshCell, ONLY: RefreshCellSetup
  USE CalStress, ONLY : CalculateStress
  USE CalForces, ONLY: CalculateForces

  USE RefreshIons, ONLY: RefreshIonTerms
  USE RefreshIons, ONLY: RescaleDensity

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
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: intCoeff
  ! coefficient for integrations
  !
  REAL(KIND=DP), DIMENSION(3) :: fracStep 
  ! Step (due to velocity) in fractional coordinates
  !
  REAL(KIND=DP) :: &
    stress(3,3),   &
    mass,          &  ! atom mass, temperay var
    maxStep,       &  ! Largest movement by an ion
    xLogS,         &  ! position of thermostat
    vLogS,         &  ! velocity of thermostat
    hamiltonian,   &  ! conserved energy
    maxForce,      &
    msd,           &  ! mean square displacement
    diffuCoeff,    &  ! diffusion coefficient
    twiceKE,       &  ! Kinetic energy x 2
    oldEtot           ! old energy

  CHARACTER(len=80) :: parameterInfo      ! Line containing values for Qmass, dt, temperature
                                          ! for output to header


                        !>> INITIALIZATION <<!
  CALL Title("MolecularDynamics::NVT")
  CALL WrtOut(6," (NVT) Start.")

  watch2 = TimerStart()
  watch = TimerStart()
  numIon = cell%numIon 

  ! calculate the freedom
  CALL CalFreedom()

  ! Calculate the number of total MD steps.
  nStep = CEILING(timeTot/dt)
  forces = 0._DP

  ! Allocate arrays for velocities, cartesian coordinates.
  ALLOCATE(velocity(3,numIon), stat=errorFlag)
  ALLOCATE(cartNoWrap(3,numIon), stat=errorFlag)
  ALLOCATE(intCoeff(nYosh), stat=errorFlag)
  
  CALL MakeIntCoeff(nYosh,intCoeff)

  ! Set up extra output to ion optimizer / MD header
  WRITE(outputUnit,'(A,ES20.12,A)') " Qmass for NVT       : ", Qmass," (a.u.)"
  WRITE(outputUnit,'(A,ES20.12,A)') " Time interval       : ", dt*fundamentalTime*10._DP**15, " (fs)"
  WRITE(outputUnit,'(A,ES20.12,A)') " Target temperature  : ", temperature/boltzmann, " (K)"
  WRITE(outputUnit,'(A,I8)') " Number of NVT steps : ", nStep

  IF (rankGlobal==0) THEN
    CALL GeometryMinimizerReportHeader(parameterInfo)
  ENDIF

  flush(6)

  IF ( rstMD == 0 ) THEN 
    WRITE(outputUnit,*) "Initialize the velocities from random."
    CALL InitRandomSeed
    CALL InitVelocity(temperature,velocity)
    vLogS = 0._DP
    xLogS = 0._DP
  ELSE IF ( (rstMD==1) .OR. (rstMD<0) ) THEN
    WRITE(outputUnit,*) "Initialize the ion positions and velocities from files."
    CALL RestartMD(1, velocity, xLogS, vLogS)
    CALL RefreshCellSetup
    CALL RescaleDensity(rhoR)
  ENDIF


        !>> FUNCTION BODY <<!


  ! mohan add 2013-07-11
  IF(rstMD<0) THEN 
    startStep=ABS(rstMD)
  ELSE
    startStep=0
  ENDIF

  ! get the kinetic energy
  CALL GetAtomKE(velocity,twiceKE)
  twiceKE = twiceKE * 2._DP

  ! Set up forces 
  CALL Optimizer                        ! Optimize the density
  CALL CalculateForces(rhoR, forces)
  IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP
  maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))

  hamiltonian = NHhamiltonian(twiceKE/2._DP, energy, xLogS, vLogS)
  oldEtot=energy(1)

  ! Print intial output
  CALL GeometryMinimizerReportSteps(startStep, TimerStop(watch), forces, maxForce, &
                                    0._DP, .TRUE., hamiltonian, .TRUE.)

  DO i=1,numIon
    cartNoWrap(:,i) = MATMUL(cell%cellReal,cell%ionTable(i)%coord)
  ENDDO


  !-----------------------------------------------
  ! big loop
  !-----------------------------------------------
  WRITE(message,'(1X,A,7X,A,11X,A,9X,A,4X,A)') "MD_STEP ", "SystemEnergy", "Conserved", "DeltaE", "Temperature"
  CALL WrtOut(6,message)

  WRITE(message,'(I8,ES20.10, ES20.10, ES15.5, ES15.5)') &
      startStep, energy(1), hamiltonian, energy(1)-oldEtot, &
      twiceKE/freedom/boltzmann
  CALL WrtOut(6,message)
 
  oldEtot=energy(1)

 
  ! Loop for MD step 
  DO step = startStep+1, nStep

    WRITE(outputUnit, '(//A)') " ******************************************************************************"
    WRITE(outputUnit,'(A,I8)') "      Nose-Hoover Molecular Dynamics (NVT) STEP ", step
    WRITE(outputUnit, '(A)')   " ******************************************************************************"

    CALL RemoveMovementOfCenterOfMass(velocity, .TRUE.)
    CALL AdjustCenterOfMass(.TRUE.) 

    ! (1) Calculate the Mean-Square-Displacement.
    CALL MonitorMeanSquareDisplacement(step,msd,diffuCoeff)
  
    ! (2) First themorstat step: updates velocity, twiceKE, xLogS, vLogS
    CALL NHIntegrator
    
    ! (3) New velocity obtained (Verlet-like step)
    DO ii=1,numIon        
      mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
      velocity(:,ii) = velocity(:,ii) + forces(ii,:,1)/mass*dt/2._DP
    ENDDO
    
    ! (4) Update the Non-Wrapped cartesion coordinates
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

    ! (5) OFDFT optimizer
    ! Refresh the ionpositions, including the Ewald term
    ! and the local pseudopotentials term.
    CALL RefreshIonTerms ! This must be called every time the 
                                    ! ion positions are changed

    ! Get the new charge density in the new ion positions.
    CALL Optimizer                        ! Optimize the density
    CALL CalculateStress(rhoR, energy, stress)
    CALL PrintStress(stress)

    ! (6) Calculate the forces. 
    CALL CalculateForces(rhoR, forces)
    IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP
    maxForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))

    ! (7) 2nd to Update velocity (Verlet-like step)
    DO ii=1,size(cell%ionTable,1)        
      mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
      velocity(:,ii) = velocity(:,ii) + forces(ii,:,1)/mass*dt/2._DP
    ENDDO


    ! (8)
    CALL GetAtomKE(velocity,twiceKE)
    twiceKE = 2._DP * twiceKE


    
    ! (9) Second thermostat step: updates velocity, twiceKE, xLogS, vLogS
    CALL NHIntegrator

    ! For monitoring MD or restarting MD   
    CALL OutputMDGeom(step)
    CALL OutputIonVelocityAndStates(step,velocity, 1,xLogS,vLogS)


    ! (10) Conserved quantity during MD 
    hamiltonian = NHhamiltonian(twiceKE/2._DP, energy, xLogS, vLogS)

! kinetic, external, coulombic, exchange-correlation,
! ion-ion, Thomas-Fermi, von Weiszacker and
! third term Wang-Teter, WGC,



    ! Print output
    CALL GeometryMinimizerReportSteps(step, TimerStop(watch), &
                                      forces, maxForce, maxStep, &
                                      .TRUE., hamiltonian, .TRUE.)
                                      

    watch = TimerStart()
    
    ! Change temperature if needed
    ! CALL ReadNewTemp( temperature, step )

    ! Output the message to the screen.
    WRITE(message,'(I8,ES20.10, ES20.10, ES15.5, ES15.5)') &
      step, energy(1), hamiltonian, energy(1)-oldEtot, twiceKE/freedom/boltzmann
    CALL WrtOut(6,message)
    oldEtot=energy(1)

  END DO
  !----------------------------------------
  ! End of NVT MD loop
  !----------------------------------------

  CALL GeometryMinimizerReportFooter(TimerStop(watch2))

  CALL WrtOut(6,"(NVT): finished.")
  
  DEALLOCATE(cartNoWrap)
  DEALLOCATE(velocity)
  DEALLOCATE(intCoeff)
  IF (ALLOCATED(msd_coords0)) DEALLOCATE(msd_coords0)
  IF (ALLOCATED(msd_coordst)) DEALLOCATE(msd_coordst)

CONTAINS


SUBROUTINE NHIntegrator
!---------------------------------------------------------------------------
! DESCRIPTION:
!   This function propagates the Nose-Hoover thermostat extended-system
!   variables using multiple step technique.
!   Note that only thermostat variables are updated, including
!   thermostat positions (xLogS),
!   thermostat velocities (vLogS),
!   thermostat acceleration (gLogS),
!   but not including any variables for atoms. 
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
                        !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: scaling, gLogS
  INTEGER :: iResn, iYosh
  REAL(KIND=DP) :: wdt2, wdt4

                        !>> INITIALIZATION <<!

  scaling = 1._DP

  ! calculate the force of thermostat
  gLogS = (twiceKE - freedom*temperature) / Qmass

                        !>> FUNCTION BODY <<!

!  WRITE(outputUnit,*) " gLogS=",gLogS
!  WRITE(outputUnit,*) " twiceKE=",twiceKE
!  WRITE(outputUnit,*) " temperature=",temperature
  
  DO iResn = 1, nResn
    DO iYosh = 1, nYosh
      wdt2 = intCoeff(iYosh)*dt/(2.d0*REAL(nResn,KIND=dp))
      wdt4 = intCoeff(iYosh)*dt/(4.d0*REAL(nResn,KIND=dp))

      ! calculate the first half of integration.
      vLogS = vLogS + gLogS*wdt4
      scaling = scaling * EXP(-wdt2 * vLogS)

      ! In some cases, the expotnential part will be NaN,
      ! Try to give warning at this case, added by mohan
      IF( check_isnan(scaling) ) THEN
        WRITE(outputUnit,*) " iResn=", iResn, " iYosh=", iYosh
        WRITE(outputUnit,*) " gLogS=", gLogS
        WRITE(outputUnit,*) " vLogS=", vLogS
        WRITE(outputUnit,'("Scaling in NHIntegrator is NaN !")')
        STOP
      ENDIF

      gLogS = (scaling**2*twiceKE - freedom*temperature) / Qmass
      ! Calculate xLogS only as a check of code
      xLogS = xLogS + vLogS*wdt2

      ! calculate the second half of integration
      vLogS = vLogS + gLogS*wdt4
    END DO
  END DO

  WRITE(outputUnit,*) " "
  WRITE(outputUnit,*) "After Nose-Hoover integrator: Scaling for velocity : ", scaling
  velocity = velocity*scaling
  twiceKE = twiceKE*scaling**2

  RETURN

END SUBROUTINE NHIntegrator


FUNCTION NHhamiltonian(KE, PE, xLogS, vLogS)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This function calculates the conserved quantity for the Nose-Hoover system.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014-05-02 Mohan update the output information
!------------------------------------------------------------------------------
 
  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: PE ! Energy from electronic optimization
  REAL(KIND=DP), INTENT(IN) :: KE ! Kinetic energy of particles (due to their velocity)
  REAL(KIND=DP), INTENT(IN) :: xLogS ! "position" of constraint variable
  REAL(KIND=DP), INTENT(IN) :: vLogS ! "velocity" of constraint variable
  REAL(KIND=DP) :: NHhamiltonian     ! The conserved quantity

  !>> LOCAL VARIABLES <<!
  REAL(KIND=DP) :: thermoKin
  REAL(KIND=DP) :: thermoPot

  !>> FUNCTION BODY <<!

  thermoKin = 0.5_DP*vLogS**2*Qmass
  thermoPot = freedom*temperature*xLogS

  NHhamiltonian = KE + PE(1) + thermoKin + thermoPot

  WRITE(outputUnit,'(A)') " "
  WRITE(outputUnit, *) "------------------------------------------------------------------------------"
  WRITE(outputUnit, *) "                             ENERGY REPORT "
  WRITE(outputUnit,'(A)') " Nose-Hoover Thermostat Details (Energy unit is eV)"
  WRITE(outputUnit,'(A,F20.6)') " NVT Temperature (K)  : ", twiceKE/freedom/boltzmann
  WRITE(outputUnit,'(A,F20.6)') " NVT Conserved Energy : ", NHhamiltonian * hartreeToeV
  WRITE(outputUnit,'(A,F20.6)') " Particle   Kinetic   Energy  : ", KE * hartreeToeV  
  WRITE(outputUnit,'(A,F20.6)') " Particle   Potential Energy  : ", PE(1) * hartreeToeV
  WRITE(outputUnit,'(A,F20.6)') " Thermostat Kinetic   Energy  : ", thermoKin * hartreeToeV
  WRITE(outputUnit,'(A,F20.6)') " Thermostat Potential Energy  : ", thermoPot * hartreeToeV

  WRITE(outputUnit,'(A)') " "
  WRITE(outputUnit,'(A)') " Potential Energy Details (Energy unit is eV): "
  WRITE(outputUnit,'(A,F20.6)') " Total    Energy : ", PE(1)*hartreeToeV 
  WRITE(outputUnit,'(A,F20.6)') " Kinetic  Energy : ", PE(2)*hartreeToeV 
  WRITE(outputUnit,'(A,F20.6)') "      TF  Term   : ", PE(7)*hartreeToeV 
  WRITE(outputUnit,'(A,F20.6)') "      vW  Term   : ", PE(8)*hartreeToeV 
  WRITE(outputUnit,'(A,F20.6)') "      NL  Term   : ", PE(9)*hartreeToeV 
  WRITE(outputUnit,'(A,F20.6)') " Ion-Elec Energy : ", PE(3)*hartreeToeV 
  WRITE(outputUnit,'(A,F20.6)') " Coulomb  Energy : ", PE(4)*hartreeToeV 
  WRITE(outputUnit,'(A,F20.6)') " Exc      Energy : ", PE(5)*hartreeToeV 
  WRITE(outputUnit,'(A,F20.6)') " Ewald    Energy : ", PE(6)*hartreeToeV 

  RETURN

END FUNCTION NHhamiltonian

END SUBROUTINE RunNVT


LOGICAL FUNCTION check_isnan(a) 
!------------------------------------------------------------------------------
! DESCRIPTION:
! check if a value is 'NAN'
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  REAL(KIND=DP) a 

  IF (a.NE.a) THEN 
    check_isnan = .true. 
  ELSE 
    check_isnan = .false. 
  END IF 

  RETURN

END FUNCTION check_isnan


END MODULE NVT

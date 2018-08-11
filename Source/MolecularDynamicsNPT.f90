MODULE NPT 
!----------------------- -------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE NPT 
!     |_SUBROUTINE NPTFullCell
!     |_SUBROUTINE MVV
!     |_SUBROUTINE NPTFullCell_Integrator
!     |_SUBROUTINE FrozeAtom
!     |_SUBROUTINE NPTFullCell_Conserved
!     |_SUBROUTINE GetAtomKETensor
!     |_SUBROUTINE GetBoxKE
!     |_SUBROUTINE BoxConstraint
!     |_SUBROUTINE EstimateQMass
!     |_SUBROUTINE EstimateBMass
!     |_SUBROUTINE CheckNPTFullCellParam
!     |_SUBROUTINE AdjustCenterOfMass
!     |_SUBROUTINE OutputStress
!     |_SUBROUTINE StressTensorKE
!
! DESCRIPTION:
!   This module contains isotermal-isobaric (NPT) 
!   molecular dynamics subroutines.
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
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2008    Chen Huang create
!   2013-04 Mohan Chen update
!------------------------------------------------------------------------------

  USE MolecularDynamics
  USE CellInfo, ONLY: cell 

  IMPLICIT NONE

  REAL(KIND=DP) :: extPres = -1.d0  
  ! the external pressure
  !
  REAL(KIND=DP) :: effTemp          
  ! effective temperature
  !
  REAL(KIND=DP) :: atomKE           
  ! ion's kinetic energy
  !
  REAL(KIND=DP) :: boxKE            
  ! box vector's kinetic energy
  !
  REAL(KIND=DP) :: conserved        
  ! the quantity need to be conserved in NPT
  !
  REAL(KIND=DP) :: bMass = -1._DP   
  ! mass of barostat 
  !
  INTEGER  :: constr_type=-1 
  ! How to apply constrain onto box change in NPT simualtion.
  ! -1: means no constraint
  !  1: the initial orthognol box is constained to be orthognol
  !  you can add your own constraints later
  !  2: the angle between a/b/c will be fixed during NPT, NOT coded yet!
  !  3: X direction will be fixed during NPT, box needs to be orthogonal
  !  4: Y direction will be fixed during NPT, box needs to be orthogonal
  !  5: Z direction will be fixed during NPT, box needs to be orthogonal

CONTAINS


SUBROUTINE NPTFullCell(RhoOptimizer)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine gives you iosthermal-isobaric ensemble.
!   (the pressure is iostropic, with given pressure on the cell) 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! Dec/21/2008 Created by Chen Huang
!
!------------------------------------------------------------------------------

  USE OUTPUT, ONLY : WrtOut
  USE SYS, ONLY : energy
  USE CalStress, ONLY : CalculateStress
  USE RefreshCell, ONLY: RefreshCellSetup
  USE RefreshIons, ONLY: RescaleDensity  
      
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!

  EXTERNAL RhoOptimizer ! A subroutine that is called to optimize the electron 
                        ! density relative to the ion positions

                        !>> INTERNAL VARIABLES <<!
  INTEGER :: ii 
  INTEGER :: step
  INTEGER :: numIon
  INTEGER :: nStep ! number of total MD steps
  CHARACTER(len=500) :: message

  !! --------------------------- KEY PROPERTIES ----------------------------------------!!
  REAL(KIND=DP) :: velocity(3,cell%numIon) ! Velocity of ions (atomic units)
  REAL(KIND=DP), ALLOCATABLE :: intCoeff(:)
  REAL(KIND=DP) :: keStress(3,3) ! 1st part in the virial stress tensor (ionic kinetic part)
  REAL(KIND=DP) :: stress(3,3)   ! 2nd part in the virial stress, the total OFDFT stress tensor
  REAL(KIND=DP) :: sysDim        ! dimension of the system
  REAL(KIND=DP) :: xLogS         ! position of thermostat
  REAL(KIND=DP) :: vLogS         ! velocity of thermostat
  REAL(KIND=DP) :: vBoxG(3,3)    ! barostat's velocity
  REAL(KIND=DP) :: msd           !  mean square displacement
  REAL(KIND=DP) :: diffuCoeff    ! diffusion coefficient
  REAL(KIND=DP) :: oldEtot = 0.0_DP


  CALL StartClock('NPT')

                        !>> INITIALIZATION <<!

  !! check key paramteres
  CALL CheckNPTFullCellParam()

  nStep = CEILING(timeTot/dt)      ! total step for MD
  numIon = cell%numIon             ! number of atoms in the cell
  sysDim = 3._dp                   ! dimension of the system

  ! Calculate the freedom
  CALL CalFreedom()

  ALLOCATE(intCoeff(nYosh))
  
  CALL MakeIntCoeff(nYosh,intCoeff)
  CALL InitRandomSeed 

  IF (rstMD==0)  THEN 
    !A fresh new MD: Do not restart MD
    CALL InitVelocity(temperature,velocity) 
    ! Initialize thermostat, and barostat
    xLogS = 0.0_DP        ! position of thermostat
    vLogS = 0.0_DP        ! velocity of thermostat
    vBoxG = 0.0_DP        ! box velocity
  ELSE IF ( (rstMD==1) .OR. (rstMD<0) ) THEN
    ! If we restart simulation from previous MD run
    ! RestartMD will read in ion coords and ion velocities from previous run
    CALL RestartMD(2,velocity,xLogS,vLogS,vBoxG)
    CALL RefreshCellSetup
    CALL RescaleDensity(rhoR)
  ENDIF ! rstMD

  ALLOCATE(cartNoWrap(3,numIon))  ! Store cartesion coords, 
                                  ! no wrapping for Periodic boundary condition
  ! Initialize cartNoWrap
  DO ii=1,numIon
    cartNoWrap(:,ii) = MATMUL(cell%cellReal,cell%ionTable(ii)%coord(:))
  ENDDO


                        !>> FUNCTION BODY <<!
  
  ! VERY IMPORTANT: We move the center of mass to the origin (0,0,0)
  ! In NPT, if center of mass is not at (0,0,0), the center of mass
  ! cannot be conserved due to coupling between barostat and atoms
  CALL AdjustCenterOfMass(.TRUE.)
  ! mohan add 2014-04-30, fix a bug.
  ! when the atom coordinates are shifted, the ion-electron
  ! potential should be regeranted.
  CALL RefreshCellSetup

  CALL BoxConstraint(vBoxG,constr_type)

  ! Compute the total H which is supposed to be conserved

  ! mohan add 2013-07-11
  IF(rstMD<0) THEN 
    startStep=ABS(rstMD)
  ELSE
    startStep=0 ! the first preparation step is called '0' step.
  ENDIF

  !------------------------------------------------------
  ! First calculation 
  !------------------------------------------------------
  CALL RhoOptimizer
  CALL GetAtomKE(velocity,atomKE)
  CALL GetBoxKE(bMass,vBoxG,boxKE)
  CALL NPTFullCell_Conserved(energy,xLogS,vLogS,sysDim)


  !------------------------------------------------------
  ! NPT main loop
  !------------------------------------------------------
  WRITE(message,'(1X,A,7X,A,11X,A,9X,A,4X,A,4X,A)') "MD_STEP ", "SystemEnergy", "Conserved", "DeltaE", "Temperature", "Volume"
  CALL WrtOut(6,message)
  WRITE(message,'(I8,ES20.10, ES20.10, ES15.5, ES15.5, ES15.5)') &
      startStep, energy(1), conserved, energy(1)-oldEtot, effTemp, cell%vol 
  CALL WrtOut(6,message)
  oldEtot=energy(1)


  CALL CalculateStress(rhoR, energy, stress)
  CALL StressTensorKE(velocity,keStress)
  CALL OutputStress( 0, stress, keStress )

  !--------------------------------------------------
  ! read in from startStep, here is the new step
  !--------------------------------------------------
  DO step = startStep+1, nStep

    WRITE(outputUnit,'(//////A)') " --------------------------------------------------"
    WRITE(outputUnit,'(A,I6)') " Molecular Dynamics (NPT) STEP ", step
    WRITE(outputUnit,'(A)') " --------------------------------------------------"


    !---------------------------------------------
    ! Monitor the mean square displacement
    ! CALL MonitorMeanSquareDisplacement(step,msd,diffuCoeff)

    
    !----------------------------------------------------------------
    ! 1st NHCP step: to update, Eq(40)  0->dt/2
    ! (1) v(t/2) (2)xLogS (3)vLogS (4)box velocity
    CALL NPTFullCell_Integrator(energy,intCoeff,xLogS,vLogS,vBoxG,velocity,stress)

    !-----------------------------------------------
    ! Modified velocity-Verlet step, to update
    ! (1) cell (2)atom velocities (3)positions of atoms           
    CALL MVV(RhoOptimizer,dt,velocity,vBoxG)

    !--------------------------------------------------------------------
    ! Second NHCP step: to update < Eq(40) >
    ! (1) v(t) (2)xLogS (3)vLogS (4)box velocity (vBoxG)
    CALL NPTFullCell_Integrator(energy,intCoeff, xLogS,vLogS,vBoxG,velocity,stress)
    
    ! Check conserved quantity
    CALL GetAtomKE(velocity,atomKE)
    CALL GetBoxKE(bMass,vBoxG,boxKE)
    CALL NPTFullCell_Conserved(energy,xLogS,vLogS,sysDim)
    CALL OutputIonVelocityAndStates(step,velocity,2,xLogS,vLogS,vBoxG)
    CALL OutputMDGeom( step ) 
    CALL StressTensorKE(velocity,keStress)
    CALL OutputStress( step, stress, keStress )
          
    ! Change temperature if needed
    ! CALL ReadNewTemp( temperature, step )

    ! Output the message to the screen.
    WRITE(message,'(I8,ES20.10, ES20.10, ES15.5, ES15.5, ES15.5)') &
      step, energy(1), conserved, energy(1)-oldEtot, effTemp, cell%vol 
    CALL WrtOut(6,message)

    oldEtot=energy(1)

  END DO
  ! End of propagating the NPT system
  !------------------------------------------------------------

  DEALLOCATE(intCoeff)
  DEALLOCATE(cartNoWrap)
  IF (ALLOCATED(msd_coords0)) DEALLOCATE(msd_coords0)
  IF (ALLOCATED(msd_coordst)) DEALLOCATE(msd_coordst)

  CALL StopClock('NPT')

  RETURN
END SUBROUTINE NPTFullCell


SUBROUTINE MVV(RhoOptimizer,dt,vel,vBoxG)
!---------------------------------------------------------------------------
! DESCRIPTION:
! This subroutine implements the modified-velocity-Verlet (MVV) to update:
!  - atoms' positions 
!  - atoms' velocities
!  - box position (barostat position)  
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!   Dec/29/2008: Created by Chen Huang
!   Apr.09 2013: Updated by Mohan Chen
!----------------------------------------------------------------------------
   USE CalForces, ONLY : CalculateForces
   USE RefreshCell, ONLY: RefreshCellSetup
   USE RefreshIons, ONLY: RescaleDensity
   USE MathFunctions, ONLY: Symmetrize
   USE MathFunctions, ONLY: Diag
   USE SYS, ONLY: forceIon
     
   IMPLICIT NONE
   
   !!>> EXTERNAL VARS <<!
   EXTERNAL RhoOptimizer  ! call it to optimize density

   REAL(KIND=DP), INTENT(INOUT) :: vel(3,cell%numIon) ! atom velocity
   REAL(KIND=DP), INTENT(IN) :: vBoxG(3,3)       ! barostat velocity
   REAL(KIND=DP) :: dt                           ! the time step to perform this NHCP 

   !! >> LOCAL VARS <<<!
   INTEGER  :: ii
   REAL(KIND=DP) :: dt2, mass
   REAL(KIND=DP) :: vTemp1(3,3)
   REAL(KIND=DP) :: vTemp2(3,3)
   REAL(KIND=DP) :: invBox(3,3)
   REAL(KIND=DP) :: ie_mat(3,3)
   REAL(KIND=DP) :: is_mat(3,3)
   REAL(KIND=DP) :: poly,arg2
   REAL(KIND=DP) :: veigv(3,3)
   REAL(KIND=DP) :: eig(3)

                                         ! >> FUNCTION BEGINS <<<!
   dt2 = dt/2._DP

   
   !----------------------------------------------------
   ! Update atom's velocities   v(0)->v(dt/2)
   ! Density coming in is already relaxed

   CALL CalculateForces(rhoR,forceIon)
   CALL FrozeAtom(forceIon(:,:,1),2)
   do ii=1,cell%numIon
      mass =  cell%elementTable(cell%ionTable(ii)%elementID)%mass
      vel(:,ii)=vel(:,ii) + dt2 * forceIon(ii,:,1)/mass
   enddo
   CALL FrozeAtom(vel,1)

   
   !------------------------------------------------------
   ! Update particle position  x(t=0)->x(t)

   ! Making propogator
   CALL Diag(3,vBoxG,eig,veigv)
   ie_mat = 0.d0
   is_mat = 0.d0
   DO ii=1,3
     ie_mat(ii,ii) = EXP(dt*eig(ii))     ! the Ie matrix in equ(50)
     ! make sinh[x]/x
     arg2 = (eig(ii)*dt2)**2
     poly = 1.d0 + arg2/6.d0 +   &
                   arg2**2/120.d0 + &
                   arg2**3/5040.d0 +  &
                   arg2**4/362880.d0 ! Taylor expansion of Sinh[x]/x
     is_mat(ii,ii) = EXP(dt2*eig(ii))*poly
   enddo
   vTemp1 = MATMUL(veigv,MATMUL(ie_mat,TRANSPOSE(veigv)))  ! 1st coeff in Eq(50)
   vTemp2 = MATMUL(veigv,MATMUL(is_mat,TRANSPOSE(veigv)))  ! 2nd coeff in Eq(50)   
   
   ! Update atom cartesion coords, x(0)->x(t)
   DO ii=1,cell%numIon
     cartNoWrap(:,ii) = MATMUL(vTemp1,cartNoWrap(:,ii))  & 
                        +  dt*MATMUL(vTemp2,vel(:,ii))
   ENDDO

   ! Monitor, center of mass position
   CALL AdjustCenterOfMass(.FALSE.)

   ! Monitor, forceIon
   WRITE(outputUnit,'(a,ES12.4,a)') " (MVV): max_force=", & 
     MAXVAL(SQRT(SUM(forceIon(:,:,1)**2,2))),' a.u.'

   ! Update the box, h(0)->h(t)
   cell%cellReal = MATMUL( vTemp1,cell%cellReal )
   CALL Symmetrize(3,cell%cellReal)  
                           ! symetrize the box, prevent accumulated numerical
                           ! in principle, box is always a symmetric matrix, so does vBoxG.
                           ! Because, qtStress is always a symmetric matrix, no net torque on box

   ! Refresh cell%cellReal AND cell%ionTable(ii)%coord(:) to make OFDFT calculations
   ! here MODULO is used for cell%ionTable(ii)%coord(:) to make Ewald more fast,
   ! however not neccessary. cartNoWrap must be used to trap the real position of atoms!
   invBox = Inverse(cell%cellReal)
   DO ii = 1, cell%numIon
     cell%ionTable(ii)%coord(:) = MODULO(MATMUL(invBox,cartNoWrap(:,ii)),1._DP)
   ENDDO
   CALL RefreshCellSetup
   CALL RescaleDensity(rhoR)
    
   !------------------------------------------------------------
   ! Update atoms' velocities: v(t/2)->v(t)
   !
   CALL RhoOptimizer




   CALL CalculateForces(rhoR,forceIon)
   CALL FrozeAtom(forceIon(:,:,1),2)
   DO ii=1, cell%numIon
     mass = cell%elementTable(cell%ionTable(ii)%elementID)%mass
     vel(:,ii)=vel(:,ii) + dt2*forceIon(ii,:,1)/mass
   ENDDO
   CALL FrozeAtom(vel,1)

   ! Udpated v(0)->v(t)
   !         x(0)->x(t)
   !         h(0)->h(t)
   
   RETURN
END SUBROUTINE MVV 


SUBROUTINE NPTFullCell_Integrator(energy,intCoeff,xLogS,vLogS,vBoxG,velocity,qtStress)
!---------------------------------------------------------------------------
! DESCRIPTION:
!   This function integrate system undergoing Parrinello-Rahman-Hoover NPT
!   dynamics. A system of N particles and the box variables is assumed to 
!   be coupled to a single Nose-Hoover thermostat.
!  
!   This is the integrator of NHCP, it updates:
!      xLogS    :    position of thermostat
!      vLogS    :    speed of thermostat
!      vBoxG    :    speed of barostat
!      velocity :    velocities of atoms
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
! 
!----------------------------------------------------------------------------

  USE MathFunctions, ONLY: Diag
  USE SYS, ONLY : frozenIon
  USE CalStress, ONLY : CalculateStress

  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN) :: energy(9) 
  ! incoming OFDFT energy with current system
  !
  REAL(KIND=DP), INTENT(OUT) :: qtStress(3,3) 
  ! internel stress
  !
  REAL(KIND=DP), INTENT(OUT) :: xLogS
  ! position of thermstat
  !
  REAL(KIND=DP), INTENT(OUT) :: vLogS
  ! speed of thermostat
  !
  REAL(KIND=DP), INTENT(OUT) :: vBoxG(3,3)
  ! speed of barostat
  !
  REAL(KIND=DP), INTENT(OUT) :: velocity(:,:)
  ! atom's velocities
  !
  REAL(KIND=DP), INTENT(IN) :: intCoeff(:)           
  ! weight for integration
  !
                      !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: pInt(3,3)  
  ! internal stress
  !
  REAL(KIND=DP) :: gLogS
  ! force on thermostat , eq(47)
  !
  REAL(KIND=DP) :: gBoxG(3,3)
  ! force on barostat   , eq(46)
  !
  REAL(KIND=DP) :: AA
  ! temperary var
  !
  REAL(KIND=DP) :: AtomKETensor(3,3)
  ! kinetic energy tensor, 1st term in Eq(46)
  !
  REAL(KIND=DP) :: boxKE
  REAL(KIND=DP) :: vTemp(3,3)

  REAL(KIND=DP) :: ie_mat(3,3)
  ! tmp array : veig * wdt2
  !
  REAL(KIND=DP) :: veigv(3,3) 
  ! eigenvalue of vTemp
  !
  REAL(KIND=DP) :: veig(3)
  ! eigenvectores of vTemp
  !
  REAL(KIND=DP) :: trvG
  ! trace of vg
  !
  REAL(KIND=DP) :: wdt2, wdt4, wdt8
  ! weight for integration
  !
  INTEGER :: numIon
  ! number of atoms
  !
  INTEGER :: iResn
  ! index for Resn
  !
  INTEGER :: iYosh
  ! index for Yosh 
  !
  INTEGER :: ii
  ! index for diagonal term in 'gBoxG'
  !

                        !>> INITIALIZATION <<!
     
  numIon = cell%numIon
  
  CALL GetAtomKE(velocity,atomKE)
  CALL GetAtomKETensor(velocity,atomKETensor)
  CALL GetBoxKE(bMass,vBoxG,boxKE)
  
  !--------------------------------------------------------
  ! force for thermostat
  gLogS = (2.0_DP*atomKE+2.0_DP*boxKE-(freedom+3.0_DP**2)*temperature)/Qmass ! eq(47)

  CALL CalculateStress(rhoR, energy, qtStress)
  !--------------------------------------------------------
  ! force for barostat
  pInt = -qtStress
  gBoxG = (2._dp*AtomKETensor + cell%Vol*pInt)/BMass ! eq(46)
  DO ii=1,3
    gBoxG(ii,ii) = gBoxG(ii,ii) + (2._dp*AtomKE/freedom - cell%Vol*extPres)/BMass
  ENDDO
  CALL BoxConstraint(gBoxG,constr_type)
   
                        !>> FUNCTION BODY <<!

  ! Start the multiple time step procedure
  DO iResn = 1, nResn
    DO iYosh = 1, nYosh
      
      !---------------------------------------------------
      ! prepare weight for integration
      wdt2 = intCoeff(iYosh)*dt/(2.d0*REAL(nResn,KIND=dp))
      wdt4 = intCoeff(iYosh)*dt/(4.d0*REAL(nResn,KIND=dp))
      wdt8 = intCoeff(iYosh)*dt/(8.d0*REAL(nResn,KIND=dp))
      ! print *,"wdt2,wdt4,wdt8:",wdt2,wdt4,wdt8
      
      !-------------------------------------------------------------------
      ! Update the thermostat velocity, and box velocity
      vLogS = vLogS + gLogS*wdt4        
      AA = EXP(-wdt8*vLogS)
      vBoxG  = vBoxG*AA*AA + wdt4*gBoxG*AA

      !-------------------------------------------------------
      ! Update the particle velocities 
      ! (1) apply the half of the thermo part, this is different from the paper
      ! why do we do this?
      ! mohan 2013-05-08
      ! AA = EXP(-wdt4*vLogS)
      ! DO ii = 1,numIon
      !  velocity(:,ii) = velocity(:,ii) * AA
      ! ENDDO
      ! CALL FrozeAtom(velocity,1)
      
      ! (2) Apply the baro part
      trvG = (vBoxG(1,1)+vBoxG(2,2)+vBoxG(3,3))/freedom
      vTemp = vBoxG
      DO ii = 1,3
        ! Chen Huang's version
        ! vTemp(ii,ii) = vTemp(ii,ii) + trvG !+ vLogS
        ! Mohan's version, original Martyna's version,
        ! but be careful, the velocity is not updated foreward and
        ! afterwards!!!
        vTemp(ii,ii) = vTemp(ii,ii) + trvG + vLogS
      ENDDO        
      CALL Diag(3,vTemp,vEig,vEigV)  !eigen value and vectores of vTemp
      ie_mat = 0.d0
      DO ii =1,3
        ie_mat(ii,ii) = EXP(-veig(ii)*wdt2)
      ENDDO

      vTemp = MATMUL(veigv,MATMUL(ie_mat,TRANSPOSE(veigv)))

      ! update the velocities of atoms
      DO ii=1,numIon
         velocity(:,ii) = MATMUL(vTemp,velocity(:,ii))
      ENDDO
      CALL FrozeAtom(velocity,1)

      ! why do we do this?
      ! mohan 2013-05-08
      ! (3) Apply another half of the thermo part to atom velocities
      ! AA = exp(-wdt4*vLogS)
      ! DO ii = 1,numIon
      !  velocity(:,ii) = velocity(:,ii) * AA
      ! ENDDO
      ! CALL FrozeAtom(velocity,1)
      ! end of updating atom's velocities

      
      !--------------------------------------------------------
      ! Update the thermo positions (xLogS)
      xLogS = xLogS + vLogS*wdt2
      
      !---------------------------------------------------------
      ! Update the box velocity (vLogG)
      ! the forces on Box, using Eq(46) in ref[1]
      call GetAtomKE(velocity,atomKE)
      call GetAtomKETensor(velocity,atomKETensor)
      gBoxG = (2._DP*AtomKETensor + cell%Vol*pInt)/bMass

      ! only diagonalization part of 
      ! Barostat force matrix (3*3)
      ! BMass : weight of Barostat (Wg in Eq(46) )
      DO ii=1,3
        gBoxG(ii,ii) = gBoxG(ii,ii) + & 
           (2._DP*atomKE/freedom - cell%Vol*extPres)/BMass
      ENDDO
      CALL BoxConstraint(gBoxG,constr_type)

      ! Update the barostat velocity (vBoxG)
      AA = EXP(-wdt8*vLogS)
      vBoxG = vBoxG * AA * AA + wdt4*gBoxG * AA
      CALL GetBoxKE(bMass,vBoxG,boxKE)

      !-----------------------------------------------
      ! Update the thermostat velocity (vLogS)
      ! Update forces on thermostat (gLogS)
      gLogS = (2._DP*atomKE+2._DP*boxKE-(freedom+3._DP**2)*temperature)/Qmass ! eq(47)
      vLogS = vLogS + gLogS * wdt4

    END DO
  END DO
  
  ! Now we have updated:
  !      xLogS    :    position of thermostat
  !      vLogS    :    speed of thermostat
  !      vBoxG    :    speed of barostat
  !      velocity :    velocities of atoms
  RETURN 

END SUBROUTINE NPTFullCell_Integrator


SUBROUTINE FrozeAtom(dataIn,flag)
!---------------------------------------------------------------------------
! DESCRIPTION:
!   Set velocity to be zero, if some ions are frozened
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------

  USE SYS, ONLY : frozenIon

  IMPLICIT NONE
  
  INTEGER :: flag, ii
  !
  REAL(KIND=DP) :: dataIn(:,:) 
  ! 1st index: x/y/z, 2nd index: ion
  !

  IF (ALLOCATED(frozenIon)) THEN 
    SELECT CASE (flag)
    CASE (1) ! Incoming is velocity
      DO ii=1, SIZE(dataIn,2)
        IF( frozenIon(ii,1) ) dataIn(1,ii)=0.D0
        IF( frozenIon(ii,2) ) dataIn(2,ii)=0.D0
        IF( frozenIon(ii,3) ) dataIn(3,ii)=0.D0
      ENDDO
    CASE (2) ! Incoming is forces
      DO ii=1, SIZE(dataIn,1)
        IF( frozenIon(ii,1) ) dataIn(ii,1)=0.D0
        IF( frozenIon(ii,2) ) dataIn(ii,2)=0.D0
        IF( frozenIon(ii,3) ) dataIn(ii,3)=0.D0
      ENDDO
    CASE DEFAULT
      PRINT *, "(FrozenAtom) flag is not defined! code stop."
      STOP
    END SELECT
  END IF

  RETURN
END SUBROUTINE FrozeAtom


SUBROUTINE NPTFullCell_Conserved(PE, xLogS,vLogS,d)
!---------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the conserved quantity for the Nose-Hoover
!   system.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------

  USE Constants, ONLY: hartreeToeV
  USE Constants, ONLY: boltzmann
  USE MathFunctions, ONLY : Det

  IMPLICIT NONE
  
  !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(:) :: PE
  ! potential energy
  !
  REAL(KIND=DP) :: thermoKin
  ! kinetic energy of thermostat
  !
  REAL(KIND=DP) :: thermoPot
  ! potential energy of thermostat
  !
  REAL(KIND=DP), INTENT(IN) :: d
  ! system's dimension, 
  !
  REAL(KIND=DP), INTENT(IN) :: xLogS
  ! thermostat position
  !
  REAL(KIND=DP), INTENT(IN) :: vLogS
  ! thermostat velocity
  !



  !>> INITIALIZATION <<!
  !>> FUNCTION BODY <<!

  ! The conserved hamiltonian should written like this
  ! the sum of the kinetic energy of ions : KE
  ! the potential energy : PE
  ! the multiply of external pressure (extPres) and 
  ! the volume of the cell (cellVol) : P*V
  ! boxKE ???
  ! the kinetic energy of bath : KE(bath)
  ! N*T*X ???
  thermoKin = 0.5_DP*vLogS**2*Qmass
  thermoPot = (freedom+d**2)*temperature*xLogS

  conserved = & 
      atomKE + PE(1) & 
    + extPres*cell%vol        &   ! related to volume change
    + boxKE                  &   ! related to the volume change
    + 0.5_DP*vLogS**2*Qmass  &   ! related to the bath
    + (freedom+d**2)*temperature*xLogS  ! related to the bath

  ! The effective temperature for ions.
  effTemp = atomKE / (freedom/2.0_DP) / boltzmann
  
  !-------------------------------
  ! Write to Log file
  !-------------------------------

  WRITE(outputUnit,'(A)') " "
  WRITE(outputUnit, *) "------------------------------------------------------------------------------"
  WRITE(outputUnit, *) "                             ENERGY REPORT "
  WRITE(outputUnit,'(A)') " NPT Details (Energy unit is eV)"
  WRITE(outputUnit,'(A,F20.6)') " Temperature (K)      : ", effTemp 
  WRITE(outputUnit,'(A,F20.6)') " NPT Conserved Energy : ", conserved * hartreeToeV
  WRITE(outputUnit,'(A,F20.6)') " Particle   Kinetic   Energy  : ", atomKE * hartreeToeV  
  WRITE(outputUnit,'(A,F20.6)') " Particle   Potential Energy  : ", PE(1) * hartreeToeV
  WRITE(outputUnit,'(A,F20.6)') " Cell Kinetic Energy          : ", boxKE
  WRITE(outputUnit,'(A,F20.6)') " Cell Volume (Bohr^3)         : ", cell%vol
  WRITE(outputUnit,'(A,F20.6)') " Ext Pressure * Vol   Energy  : ", extPres*cell%vol*hartreeToeV
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
  
END SUBROUTINE NPTFullCell_Conserved


SUBROUTINE GetAtomKETensor(vel,keTensor)
!---------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the classical kinetic energy 
!   tensor of a group of atoms, used in NPTFullCell subroutine.
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

  REAL(KIND=DP),INTENT(IN) :: vel(3,int(cell%vol))  
  ! Velocities of atoms
  !
  REAL(KIND=DP),INTENT(OUT) :: keTensor(3,3)
  ! Tensor of kinetic term

                        !>> INTERIOR VARIABLES <<!

  INTEGER :: ion, ii , jj
  ! index for ion and tensor index
  !
  REAL(KIND=DP) :: ionMass
  ! tmp mass for each element
  !
                        !>> FUNCTION BODY <<!

  ! Reset the keTensor
  keTensor = 0.0_DP
  DO ii=1,3
    DO jj=1,3
      DO ion = 1, cell%numIon 
         ionMass = cell%elementTable(cell%ionTable(ion)%elementID)%mass
         keTensor(ii,jj) = keTensor(ii,jj) + & 
           0.5_DP * ionMass * vel(ii,ion) * vel(jj,ion)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE GetAtomKETensor


SUBROUTINE GetBoxKE(bMass,vBoxG,boxKE)
!----------------------------------------------------------------------------
! DESCRIPTION:
!  Compute the kinetic energy of the Box, in NPTFullCell simulation
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
 
 REAL(KIND=DP), INTENT(IN) :: bMass       
 ! box (barostat) 's mass
 !
 REAL(KIND=DP), INTENT(IN) :: vBoxG(3,3)  
 ! box's velocity
 !
 REAL(KIND=DP), INTENT(OUT) :: boxKE      
 ! box's Kinetic energy
 !

 ! >> FUNCTION BEGINS <<!
 !! boxKE = 1/2 * bMass * Tr[vBoxG^T*vBoxG] 
 
 boxKE = 0.50_DP * bMass * SUM(vBoxG**2) ! Tr[Vg*Vg] eq(47)
  
 RETURN      
END SUBROUTINE GetBoxKE


SUBROUTINE BoxConstraint(force_box,constr_type)
!----------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine apply constrains onto the box which is a 3*3 matrix
!   with the column vectors to the lattice vectors. For example, if we do
!   not want to keep box to be a cube durint the MD, we can set the 
!   offdiagonal terms to be zero.
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

  !!! >> EXTERNAL VARIABLES << !!!

  REAL(KIND=DP), INTENT(INOUT) :: force_box(3,3) 
  ! the force on box
  !
  INTEGER :: constr_type, ii, jj
  !

  !----- LOCAL VARS ----------
!  real(dp) :: alp,bet,gam   ! Angles between (c,b),(a,c),(b,a)
!  real(dp) :: norm_a, norm_b, norm_c
!  real(dp) :: a(3),b(3),c(3)

  !!! >> FUNCTION BEGINS << !!!
  IF (constr_type < 0) THEN
    ! we do not apply any constraint on the box
    ! box is free to change both size & shape during NPT DM
    RETURN
  ENDIF


  SELECT CASE(constr_type)
  CASE(1) ! keep the box to be orthogonal
    force_box(1,2) = 0.d0
    force_box(1,3) = 0.d0
    force_box(2,1) = 0.d0
    force_box(2,3) = 0.d0
    force_box(3,1) = 0.d0
    force_box(3,2) = 0.d0
    
  CASE(2) ! angle between a,b,c will be preserved
    PRINT *, 'Have not implemented yet, stop'
    STOP
  
  CASE(3) ! keep the box fixed in X direction
    force_box(1,1) = 0.d0
    ! remove rotation of the cell
    DO ii=1,3
      DO jj=1,3
        IF (ii/=jj) & 
          force_box(ii,jj) = 0.d0
      ENDDO
    ENDDO

  CASE(4) ! keep the box fixed in Y direction
    force_box(2,2) = 0.d0
    ! remove rotation of the cell
    DO ii=1,3
      DO jj=1,3
        IF (ii/=jj) & 
          force_box(ii,jj) = 0.d0
      ENDDO
    ENDDO

  CASE(5) ! keep the box fixed in Z direction
    force_box(3,3) = 0.d0
    ! remove rotation of the cell
    DO ii=1,3
      DO jj=1,3
        IF (ii/=jj) & 
          force_box(ii,jj) = 0.d0
      ENDDO
    ENDDO

  CASE(6) ! keep the box fixed in X and Y directions
    force_box(1,1) = 0.d0
    force_box(2,2) = 0.d0
    DO ii=1,3
      DO jj=1,3
        IF (ii/=jj) & 
          force_box(ii,jj) = 0.d0
      ENDDO
    ENDDO

  CASE(7) ! keep the box fixed in X and Z directions
    force_box(1,1) = 0.d0
    force_box(3,3) = 0.d0
    DO ii=1,3
      DO jj=1,3
        IF (ii/=jj) & 
          force_box(ii,jj) = 0.d0
      ENDDO
    ENDDO

  CASE(8) ! keep the box fixed in Y and Z directions
    force_box(2,2) = 0.d0
    force_box(3,3) = 0.d0
    DO ii=1,3
      DO jj=1,3
        IF (ii/=jj) & 
          force_box(ii,jj) = 0.d0
      ENDDO
    ENDDO
    

  CASE DEFAULT
    PRINT *, 'Source/MolecularDynamics.f90 >>> constr_type = invalid value', & 
             'code stop'
    STOP
  END SELECT

  RETURN
END SUBROUTINE BoxConstraint



SUBROUTINE EstimateQMass(Nf,temp,tau,Qmass)
!---------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine estimate the Qmass used for the thermostat,
!    QMass = N_freedom * kT * tau_Qmass^2
!   where tau_QMass is the period for themostat fluctuation
!
! REFERENCES:
!   The PINY_MD code
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

   INTEGER :: Nf                
   ! number of freedom
   !
   REAL(KIND=DP), INTENT(IN) :: temp 
   ! temperature, already multiplied with boltzmann constatn
   !
   REAL(KIND=DP), INTENT(IN) :: tau  
   ! period in a.u
   !
   REAL(KIND=DP), INTENT(OUT):: QMass
   ! final answer

   QMass = REAL(Nf,KIND=dp)*temp*tau*tau

   RETURN

END SUBROUTINE EstimateQMass


SUBROUTINE EstimateBMass(Nf,temp,tau,Bmass)
!---------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine estimate the Bmass used for the thermostat,
!
!    BMass = (N_freedom + d) * kT * tau_Bmass^2
!
!   where tau_BMass is the period for barostat fluctuation
!   d = 3 for 3-dimention system
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

  INTEGER :: Nf 
  ! number of freedom
  !
  REAL(KIND=DP), INTENT(IN) :: temp 
  ! temperature, already multiplied with boltzmann constatn
  !
  REAL(KIND=DP), INTENT(IN) :: tau  
  ! period in a.u
  !
  REAL(KIND=DP), INTENT(OUT):: BMass
  ! the final answer

  BMass = REAL(Nf+3,KIND=dp)*temp*tau*tau

  RETURN

END SUBROUTINE EstimateBMass


SUBROUTINE CheckNPTFullCellParam()
!---------------------------------------------------------------------------
! DESCRIPTION:
!   Check the input parameters before doing NPT_fullCell simulation
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!----------------------------------------------------------------------------

  USE Constants, ONLY: auPressure
  USE Constants, ONLY: boltzmann

  IMPLICIT NONE
  
  ! check initialized variables.
  IF (nResn<0 .OR. &
      nYosh<0 .OR. & 
      extPres     == -1.0_DP .OR. &
      dt          == -1.0_DP .OR. &
      timeTot     == -1.0_DP .OR. &
      temperature == -1.0_DP .OR. &
      (bMass == -1.0_DP .AND. tau_baro   == 0.0_DP ) .OR.  &
      (Qmass == -1.0_DP .AND. tau_thermo == 0.0_DP ) ) THEN
    WRITE(outputUnit,*) "Source/MolecularDynamics.f90 >>> " ,  &
     "one key parameter for NPT is not initialized, error, code stop"
    PRINT *,"nResn,nYosh,extPres,dt,timeTot,temperature,bMass,tau_baro,Qmass,tau_thermo"
    PRINT *,nResn,nYosh,extPres,dt,timeTot,temperature,bMass,tau_baro,Qmass,tau_thermo
    STOP
  ENDIF

  IF ( tau_thermo > 0.0_DP ) THEN
    CALL EstimateBMass(cell%numIon*3,temperature,tau_thermo,Qmass)
    WRITE(outputUnit,'(a,f9.1,a,a,ES8.1,a)') &
     " tau_thermo:",tau_thermo*fundamentalTime/1.e-12," ps", & 
     " => QMass (est.): ",QMass, " a.u."
  ENDIF

  IF ( tau_baro > 0.0_DP ) THEN
    CALL EstimateQMass(cell%numIon*3,temperature,tau_baro,Bmass)
    WRITE(outputUnit,'(a,f9.1,a,a,ES8.1,a)') &
     " tau_baro  :",tau_baro*fundamentalTime/1.e-12," ps", & 
     " => BMass (est.): ",BMass, " a.u."
  ENDIF
  
  WRITE(outputUnit,'(/a)') " --------------------------------------------------"
  WRITE(outputUnit,'(a)')  " NPT(MD) PARAMETERS"
  WRITE(outputUnit,'(a)')  " --------------------------------------------------"
  WRITE(outputUnit,'(a,ES12.4,a)')  " QMass   = ", QMass, " a.u."
  WRITE(outputUnit,'(a,ES12.4,a)')  " BMass   = ", BMass, " a.u."
  WRITE(outputUnit,'(a,ES12.4,a)')  " tau_T   = ", tau_thermo*fundamentalTime, " (sec)"
  WRITE(outputUnit,'(a,ES12.4,a)')  " tau_B   = ", tau_baro*fundamentalTime, " (sec)"
  WRITE(outputUnit,'(a,ES12.4,a)')  " totTime = ", timeTot*fundamentalTime, " (sec)"
  WRITE(outputUnit,'(a,ES12.4,a)')  " dt      = ", dt*fundamentalTime     , " (sec)"
  WRITE(outputUnit,'(a,ES12.4,a)')  " temp    = ", temperature/boltzmann,   " (K)"
  WRITE(outputUnit,'(a,ES12.4,a)')  " extPres = ", extPres/auPressure,      " (Pa)"
  WRITE(outputUnit,'(a,2I2)')       " nResn, nYosh       = ", nResn, nYosh
  WRITE(outputUnit,'(a, I2)')       " NPT constrain type = ", constr_type

END SUBROUTINE CheckNPTFullCellParam


SUBROUTINE OutputStress(iter,stress_ofdft,stress_ke)
!---------------------------------------------------------------------------
! DESCRIPTION:
!  Output the stress every output_stress_period steps.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------
  USE Output, ONLY: PrintStress
  USE CONSTANTS, ONLY: auToGPa ! Conversion to gigapascals
  USE CONSTANTS, ONLY: Bohr

  IMPLICIT NONE

  INTEGER :: iter, i
  ! index
  !
  REAL(KIND=DP) :: stress_ofdft(3,3) 
  ! stress from OFDFT calculations
  !
  REAL(KIND=DP) :: stress_ke(3,3) 
  ! stress from kinetic term
  !
  CHARACTER(LEN=30) :: comp
  ! 

  IF (output_stress_period == 1 .OR. mod(iter, output_stress_period)==1 ) THEN 
    WRITE(outputUnit, *) " "
    WRITE(outputUnit, *) "--------------------------------------------------"
    WRITE(outputUnit, *) "LATTICE INFORMATION (Angstrom)"
    WRITE(outputUnit, *) "--------------------------------------------------"
    DO i = 1, 3
      IF(i == 1) comp = "LA_A   1"
      IF(i == 2) comp = "LA_A   2"
      IF(i == 3) comp = "LA_A   3"
      WRITE(outputUnit, '(1X, A8, 1X, 3(ES18.10, 1X))') comp(1:8),cell%cellReal(:,i) * bohr
    END DO


    WRITE(outputUnit, *) " "
    WRITE(outputUnit, *) "--------------------------------------------------"
    WRITE(outputUnit, *) "Barostat STRESS (GPa)"
    WRITE(outputUnit, *) "--------------------------------------------------"
    DO i = 1, 3
      IF(i == 1) comp = "ST_GPa 1"
      IF(i == 2) comp = "ST_GPa 2"
      IF(i == 3) comp = "ST_GPa 3"
      WRITE(outputUnit, '(1X, A8, 1X, 3(ES18.10, 1X))') comp(1:8), & 
         stress_ke(i,:) * auToGPa 
    END DO
  

    WRITE(outputUnit, *) " "
    WRITE(outputUnit, *) "--------------------------------------------------"
    WRITE(outputUnit, *) "STRESS (GPa)"
    WRITE(outputUnit, *) "--------------------------------------------------"
    DO i = 1, 3
      IF(i == 1) comp = "ST_GPa 1"
      IF(i == 2) comp = "ST_GPa 2"
      IF(i == 3) comp = "ST_GPa 3"
      WRITE(outputUnit, '(1X, A8, 1X, 3(ES18.10, 1X))') comp(1:8), & 
         stress_ofdft(i,:) * auToGPa 
    END DO


    WRITE(outputUnit, *) " "
    WRITE(outputUnit, *) "--------------------------------------------------"
    WRITE(outputUnit, *) "TOTAL STRESS (GPa)"
    WRITE(outputUnit, *) "--------------------------------------------------"
    DO i = 1, 3
      IF(i == 1) comp = "ST_GPa 1"
      IF(i == 2) comp = "ST_GPa 2"
      IF(i == 3) comp = "ST_GPa 3"
      WRITE(outputUnit, '(1X, A8, 1X, 3(ES18.10, 1X))') comp(1:8), & 
         (stress_ofdft(i,:)+stress_ke(i,:)) * auToGPa 
    END DO

    !CALL PrintStress(stress_ofdft)
    RETURN 
  ENDIF

END SUBROUTINE OutputStress


SUBROUTINE StressTensorKE(vel,keStress)
!-----------------------------------------------------------------
! DESCRIPTION:
! The stress tensor part due to the kinetic energy
! part, the total stress tensor by virial theorem is 
!  sigma_tot = sigma_ke + sigma_pot
! the sigma_pot is just the OFDFT stress tensor
! which contains ion-ion, ion-ele, KEDF, XC, Hatree ...
!
! According to wikipeida, and is very easy to prove:
! No minus sign, and consistent with Robin's thesis (ask.
! Emily for Robin's thesis on OFDFT stress@page 166)
!
!   sigma_ke_ij = 1/Volume * sigma_i (mass_i*vel_i*vel_j)
!  
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!
!-----------------------------------------------------------------

  IMPLICIT NONE
  
  REAL(KIND=DP), INTENT(IN) :: vel(3,cell%numIon)         
  ! velocity of ions
  !
  REAL(KIND=DP), INTENT(out) :: keStress(3,3)
  ! stress from kinetic term
  
  REAL(KIND=DP):: mass
  ! tmp mass
  !
  INTEGER :: i,j,a,id
  ! index
  !
  
  keStress = 0.d0
  DO a = 1, cell%numIon 
    id = cell%ionTable(a)%elementID
    mass = cell%elementTable(id)%mass
    DO i = 1,3
      DO j = 1,3         
        keStress(i,j) = keStress(i,j) + mass*vel(i,a)*vel(j,a)
      ENDDO
    ENDDO
  ENDDO

  keStress = keStress / cell%vol

  RETURN

END SUBROUTINE StressTensorKE


END MODULE NPT

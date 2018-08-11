PROGRAM CINEB
!------------------------------------------------------------------------------
! STRUCTURE OF MAIN PROGRAM:
!   PROGRAM CINEB
!
! INPUT SUBROUTINE
!     |_SUBROUTINE ReadCinebInputFile
!     |_SUBROUTINE ReadCinebParamFile
!
! CALCULATION SUBROUTINE
!     |_SUBROUTINE TransitionStateVelocityVerlot
!     |_SUBROUTINE TotalCalculation
!     |_SUBROUTINE TrueCalculation
!     |_SUBROUTINE GetSpringConstant
!     |_SUBROUTINE GetTotalForce
!     |_SUBROUTINE HenkelTangent
!     |_SUBROUTINE HenkelForce
!     |_SUBROUTINE SplineTangent
!     |_SUBROUTINE SplineForce
!
! MATHEMATICAL SUBROUTINE
!     |_SUBROUTINE Spline
!     |_SUBROUTINE Spline2
!     |_SUBROUTINE Recenter
!     |_SUBROUTINE Rotate
!     |_FUNCTION CartToFract
!     |_FUNCTION FractToCart
!     |_SUBROUTINE GaussianElimination
!
! OUTPUT SUBROUTINE
!     |_SUBROUTINE OutputWriteHeader
!     |_SUBROUTINE OutputWrite
!     |_SUBROUTINE OutputWriteFooter
!     |_SUBROUTINE RestartPosition
!
! DESCRIPTION:
!   This is a Climbing Image Nudged Elastic Band code that determines 
!   transistion states using the ORBITAL FREE DFT code.  It moves the ions
!   around by itself, and performs single point calculations by calling OFDFT.
!
! CONDITIONS AND ASSUMPTIONS:
!   If you have trouble with this, there are a few things you may want to do:
!     1.  Play with the climbing image.  Perhaps turn it on only after
!         convergence of the NEB.
!     2.  String a band between the two states that are of highest energy,
!         rerun the calculation
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Someone may want to integrate this directly into the OFDFT code someday,
!   probably when it has been parallelized.
!
!   Things to do:
!     Clean up Code:
!       - Comment all variables
!       - Capitalize properly
!       - Restrict lines to about 80 chars
!       - Fill in headers of functions
!       - Put external variables in right place, label with (INTENT)
!     Implement variable spring constants?  spring contant depends linearly
!       on the energies?
!     Figure out whether the "improved" tangent method is really better?
!     Can someone please make it so that you don't have to specify # of atoms
!       in the unit cell????!!  Or the number of atom types or the number of
!       images, for that matter.  And make a prettier parameter file?
!     Recenter and rotate -- do they work?
!     Eliminiate the need for an ORBITALFREE.inp file by merging with 
!       CINEB.param  (Speaking of which, CINEB.inp is so terrible right now)
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2/10/2005 Kyle Caspersen originally wrote this program for the OLD OFDFT
!             code.  Greg Ho is now adapting and integrating it into newOFDFT, 
!             as well as prettifying / commenting it.
!   4/19/2006 Merged with CINEBMOD
!
!------------------------------------------------------------------------------
  USE CONSTANTS, ONLY : &
    DP, &                  ! Double precision 
    pi, &                  ! You know, pi
    AtomicMass             ! The atomic mass table for atoms

  USE MATHFUNCTIONS, ONLY : &
    Inverse, &             ! Performs the inverse of a 3x3 matrix
    Volume, &              ! A function to find the volume
    Cross                  ! Performs a cross product 

  IMPLICIT NONE

                          !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    numAtoms, &            ! The number of atoms in the cell
    numAtomTypes, &        ! The number of types of atoms in cell
    numImages, &           ! The number of images to use in the NEB
    nudgeType, &           ! value : -1  ->  NO BAND
                           !       :  0  ->  ELASTIC BAND
                           !       :  1  ->  NEB ("Improved" tangent)
                           !       :  2  ->  CINEB ("Improved tangent")
                           !       :  3  ->  NEB (Original spline tangent)
                           !       :  4  ->  CINEB (Original spline tangent)
                           !       :  5  ->  NEB (New spline tangent)
                           !       :  6  ->  CINEB (New spline tangent)
    initialConditions, &   ! The initial conditions:  
                           !       1-> use initial conditions straight up
                           !       2-> relax first and last images with 
                           !           QUICKMIN 
    transOutputUnit = 22, &       ! The "simplified" output file
    transImagesOutputUnit = 23, & ! The detailed output file with per image 
                                  ! data
    transImagesOutputForcesUnit = 24, &
    openStatus

  REAL(KIND=DP) :: &
    forceCut            ! The criteria for the force to be converged on a
                           ! relaxed band.

  REAL(KIND=DP), DIMENSION(3,3) :: &
    cell                   ! the cell lattice parameters

  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    cartPosition           ! the ion positions in cartesian coordinates

  CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: &
    atomicSymbol           ! the atomic symbol

  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: &
    atomicPSP              ! the pseudopotential file

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: &
    frozenAtoms

  CHARACTER(LEN=1) :: &
    posType               ! "F" if we're doing fractional coordiates, "C" for
                           ! cartesian

  CHARACTER(LEN=3) :: &
    numProc

                         !>> INITIALIZATION <<!
                          !>> FUNCTION BODY <<!


  CALL GETARG(1, numProc)
  IF(numProc == "") numProc = "1  "

  ! Open the output file
  OPEN(unit=transOutputUnit, position="rewind", status="replace", & 
       action = "write", blank="null", file="cineb.out", iostat=openStatus)

  ! Read input files
  CALL ReadCinebInputFile
  CALL ReadCinebParamFile

  ! Write out the initial ion coordinates of all the images.
  !CALL WriteRestartInformation(cartPosition, cell, atomicPSP, atomicSymbol)

  ! Use minimization scheme to find the minimum energy path
  CALL OutputWriteHeader(cell)



  CALL TransitionStateVelocityVerlot

  STOP ! mohan test
  CALL OutputWriteFooter

  CLOSE(transOutputUnit)
!  CLOSE(transImagesOutputUnit)
!  CLOSE(transImagesOutputForcesUnit)

CONTAINS


SUBROUTINE ReadCinebInputFile
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Reads the .inp file
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!
                            !>> INTERNAL VARIABLES <<!
  CHARACTER(LEN=90) :: &
    chjunk              ! Junk character

  INTEGER :: &
    ni, &
    m,n, &
    readStatus, &
    oldni               ! Temp variable to store the number of the last
                        !  image

                            !>> INITIALIZATION <<!
                             !>> FUNCTION BODY <<!

  ! This is a very simplistic input routine we may want to change?
  !WRITE(*,*) "WARNING: AFTER SWITCHING TO NEW PATHSCALE COMPILER, THE FOLLOWING line is suspected not to work. Please check and remove this comment."

  OPEN(UNIT=1,FILE="CINEB.inp")
  READ(1,'(A)') chjunk

  ! Cell lattice parameters
  READ(1,*) cell(1,1),cell(2,1),cell(3,1)
  READ(1,*) cell(1,2),cell(2,2),cell(3,2)
  READ(1,*) cell(1,3),cell(2,3),cell(3,3)

  WRITE(*,*) cell

  ! Number of atoms and number of atom types
  READ(1,*) chjunk,numAtoms
  READ(1,*) chjunk,numatomtypes

  ! Number of images
  READ(1,*) chjunk,numImages

  WRITE(*,*) "atom type: ", numatomtypes
  WRITE(*,*) "image num: ", numImages

  ALLOCATE(cartPosition(numImages,numAtoms,3))
  ALLOCATE(atomicSymbol(numAtoms))
  ALLOCATE(atomicPSP(numatomtypes))
  ALLOCATE(frozenAtoms(numAtoms,3))

  DO m = 1, numatomtypes
    READ(1,*) chjunk, atomicPSP(m)
    atomicPSP(m) = chjunk(1:2)//" "//atomicPSP(m)
    WRITE(*,*) atomicPSP(m)
  END DO

  ! Read whether we're doing FRACTIONAL or CARTESIAN coordinates
  READ(1,'(A1)') posType
  WRITE(*,*) "postype = ", posType


  ! Read the first image
  READ(1,*,iostat=readStatus) chjunk,ni
  IF(readStatus == 1) THEN
    WRITE(*,*) "***ERROR IN CINEB.inp on image 1 header ***"
    STOP
  END IF
  oldni=ni

  ! Read all the ions
  DO m=1,numAtoms
    READ(1,*,iostat=readStatus) atomicSymbol(m), cartPosition(ni,m,1), &
              cartPosition(ni,m,2), cartPosition(ni,m,3)
    IF(readStatus == 1) THEN
      WRITE(*,*) "***ERROR IN CINEB.inp on image 1, atom ", m
      STOP
    END IF
  END DO


  ! Now read all the other images.  I'm not sure why this is a seperate loop
  ! that's different from the main loop, except that the atomic symbol, the
  ! mass and charge are JUNKED.

  ! Note that this loop fails if we don't specify the last image.
  DO n=2,numImages

    READ(1,*,iostat=readStatus) chjunk,ni
    IF(readStatus == 1) THEN
      WRITE(*,*) "***ERROR IN CINEB.inp on image ", n," header ***"
      STOP
    END IF
    DO m=1,numAtoms
      READ(1,*,iostat=readStatus) &
        chjunk, cartPosition(ni,m,1), cartPosition(ni,m,2), &
                cartPosition(ni,m,3)
      IF(readStatus == 1) THEN
        WRITE(*,*) "***ERROR IN CINEB.inp on image ",n, ", atom ", m
        STOP
      END IF
    END DO

    ! If we're skipping images (we're not specifying ALL images from 1 to 
    ! numImages) then linearly interpolate all the images in between.
    IF (ni-oldni>1) THEN
      DO m=oldni+1,ni-1
        cartPosition(m,:,:) = cartPosition(oldni,:,:) + &
          ((cartPosition(ni,:,:) - cartPosition(oldni,:,:)) &
          * (REAL(m-oldni,KIND=DP)/REAL(ni-oldni, kind=DP)))
      END DO
    END IF
    
    oldni=ni
    IF (ni>=numImages) EXIT

  END DO

  ! Read frozen atoms
  READ(1,*,iostat=readStatus) chjunk
  WRITE(*,*) chjunk
  IF(readStatus == 1) THEN
    WRITE(*,*) "***ERROR IN CINEB.inp on frozen atoms header ***"
    STOP
  END IF
  DO m=1,numAtoms
    READ(1,*,iostat=readStatus) frozenAtoms(m,1), frozenAtoms(m,2), frozenAtoms(m,3)
    IF(readStatus == 1) THEN
      WRITE(*,*) "***ERROR IN CINEB.inp on frozen atoms ", m
      STOP
    END IF
  END DO

  ! If the ion positions are in FRACTIONAL coordinates:
  IF (posType=="F") THEN
    cartPosition = FractToCart(cell, cartPosition)

  ! Ion positions are in cartesian coordinates.
  ELSE IF (posType=="C") THEN

  ELSE
    PRINT *,"Attention you have an error with the position type ... SO FIX IT."
  END IF

  CLOSE(1)

  RETURN

END SUBROUTINE ReadCinebInputFile


SUBROUTINE ReadCinebParamFile
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Reads the .param file
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Haha, you're kidding, right?  This function is the epitome of laziness.
!   It's really, really sad. 
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!
                            !>> INTERNAL VARIABLES <<!
  CHARACTER(LEN=90) :: label         ! Junk character
    
                            !>> INITIALIZATION <<!
                             !>> FUNCTION BODY <<!

  OPEN(UNIT=777,FILE="CINEB.param")
  READ(777,*) label, nudgeType
  READ(777,*) label, initialConditions
  READ(777,*) label, forceCut
  WRITE(*,*) "NUDGE_TYPE=", nudgeType
  WRITE(*,*) "INIT_COND=", initialConditions
  WRITE(*,*) "FORCE_CUT=", forceCut
  CLOSE(777)

END SUBROUTINE ReadCinebParamFile


SUBROUTINE TransitionStateVelocityVerlot()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Does a QUICKMIN search to try to minimize the total force (sum of true
!   force and spring force) in the band.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/18/2005 Changed to obey correct projectile motion laws, optimized, 
!             commented (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!
                           !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    ni, na, m, &          ! Counters
    linecounts, &         ! number of steps in this linesearch
    numiterations, &      ! total number of linesearches done
    numCalcs              ! number of times we have ran OFDFT

  REAL(kind=DP) :: &
    maxForce, &           ! The maximum force for any one ion within all the
                          !   images of the system
    lastMaxForce, &       ! The last maximum force
    aveForce, &
    lastAveForce, &
    minAveForce, &
    minMaxForce, &        ! The minimum max force so far
    dotP, &               ! Dot product of the initial force in linesearch
                          !   and the current force
    atomLength, &         ! this is approximately the length in one direction 
                          !   for one atom, assuming atoms are about evenly 
                          !   distributed.  Used for checking the initial
                          !   timestep.
    timeStep              ! our timestep

  REAL(kind=DP), DIMENSION(numImages) :: &
    springConstant        ! our spring constant

  REAL(kind=DP), DIMENSION(SIZE(cartPosition,1), SIZE(cartPosition,2), &
                           SIZE(cartPosition,3)) :: &
    initialForce, &       ! The array of forces for the initial step in the
                          !   linesearch
    totalForce, &         ! The array of forces for the current state
    forceOld, &           ! The array of forces for the previous step in the 
                          !   linesearch
    cartOld, &            ! The atom positions (cartesian) for previous step  
                          !   in the linesearch
    velocity, &           ! The velocity
    springForce, &        ! The spring force
    trueForce, &          ! The true force
    stressTensor,&        ! The stress tensor
    positionTangent, &    ! The tangent in the CINEB
    massPosition          ! Center of mass positions of the ions

  REAL(kind=DP), DIMENSION(numImages) :: &
    energy                ! Storage for the energy at each state

                            !>> INITIALIZATION <<!

  !timeStep = .49_DP          ! This is our initial timestep
  timeStep = .10_DP          ! This is our initial timestep
  numCalcs = 0
  lineCounts = 0

  ! If we change the cell lattice vectors, atomLength will have to be
  ! refreshed every loop
  atomLength=(Volume(cell)/REAL(numAtoms,kind=DP))**(1._DP/3._DP)

  numiterations=0     ! Number of total linesearches
 
  ! Calculate the initial forces    
  CALL TotalCalculation(numCalcs, cell, cartPosition, atomicSymbol, &
                        atomicPSP, trueForce, springForce, totalForce, &
                        energy, stressTensor, springConstant)

  ! Writes the initial starting data to the output file
  CALL OutputWrite(cartPosition, massPosition, energy, trueForce, &
                   springForce, totalForce, stressTensor, numCalcs, &
                   numiterations, 0, .FALSE., 0._DP)


                             !>> FUNCTION BODY <<!

  ! For convergence, here is a steepest-descent linesearch algorithm.
  DO

    numiterations=numiterations+1
    ! Define the initial state before the line search
    initialForce = totalForce
 
    ! Make sure our time increment is not too big.  We don't want the
    ! atoms flying all over the place.  Here, we make sure that the timestep
    ! doesn't cause any one atom to move more than 10% the average length 
    ! of an atom in the cell, assuming bulk densities.
    ! Obviously, this is heuristic, and may need to be changed to be more
    ! strict for, say, vacuum systems, where this is a poor approximation.

    ! Figure out the maximum force on one ion in any of the images
    maxForce = SQRT(MAXVAL(SUM(totalForce**2,3)))
    aveForce = SQRT(SUM(totalForce**2)/numImages/numAtoms)

    IF(numIterations > 1) THEN
      minMaxForce = MIN(maxForce, minMaxForce)
      minAveForce = MIN(aveForce, minAveForce)
    ELSE
      minMaxForce = 2._DP*maxForce
      minAveForce = 2._DP*aveForce
    END IF
    IF(numIterations == 1) THEN
      lastMaxForce = maxForce * 2._DP
      lastAveForce = aveForce * 2._DP
    END IF

    ! Limit the timestep to any one atom traveling 10% of the space
    ! in the cell it is likely confined to
    IF (.5_DP *maxForce*timeStep**2 > 1._DP) THEN
      timeStep = SQRT(2._DP*maxForce)
    ENDIF


    ! Begin linesearch
    linecounts=0
    velocity = 0._DP
    DO

      linecounts=linecounts+1

      ! Remember our current position
      cartOld=cartPosition
      forceOld=totalForce

      ! Move the ions .... d = d0 + vt + .5 a t^2
      ! Project the velocity in the direction of the force
      cartPosition(1:numImages,:,:) = &
        cartPosition(1:numImages,:,:) + &
        velocity(1:numImages,:,:) * timeStep + &
        .5_DP * totalForce(1:numImages,:,:) * timeStep**2

      ! Update the current instantaneous velocity with previous force
      ! .... v = vo + at(1/2)
      velocity(1:numImages,:,:) = &
        velocity(1:numImages,:,:) + &
        .5_DP * totalForce(1:numImages,:,:) * timeStep

      ! Calculate the new forces and energies
      CALL TotalCalculation(numCalcs, cell, cartPosition, atomicSymbol, &
                            atomicPSP, trueForce, springForce, totalForce, &
                            energy, stressTensor, springConstant)

      ! Calculate the dot product and the maximum forces
      dotp = SUM(initialForce(1:numImages,:,:)*totalForce(1:numImages,:,:))
      lastMaxForce = maxForce
      lastAveForce = aveForce
      maxForce = SQRT(MAXVAL(SUM(totalForce**2,3)))
      aveForce = SQRT(SUM(totalForce**2)/numImages/numAtoms)
      minAveForce = MIN(aveForce, minAveForce)

      ! Update the current instantaneous velocity with this force 
      ! .... v = vo + at(1/2)
      velocity(1:numImages,:,:) = &
        velocity(1:numImages,:,:) + &
        .5_DP * totalForce(1:numImages,:,:) * timeStep
      
      ! Project the velocity only in the direction of the forces
      velocity(1:numImages,:,:) = &     
        SUM(velocity(1:numImages,:,:) * totalForce(1:numImages,:,:)) &
        / SUM(totalForce(1:numImages,:,:)**2) * &
        totalForce(1:numImages,:,:)

      ! Write the output files the new forces, etc.
      CALL OutputWrite(cartPosition, massPosition, energy, trueForce, &
                       springForce, totalForce, stressTensor, numCalcs, &
                       numiterations, linecounts, dotp < 0._DP, timeStep)

      IF (ABS(lastMaxForce-maxForce)<1.E-6_DP) THEN
        WRITE(*,'(A,F15.10,A,F15.10,A)') '# diffForce < 1E-6 (lastMaxForce=', lastMaxForce, ', maxForce=', maxForce, ') => EXIT'
        timeStep=.49_DP
        velocity=0._DP
        linecounts=0
        EXIT
      END IF

      IF (lineCounts == 1 .AND. ((dotp < 0._DP .AND. maxForce > forceCut)))THEN
        WRITE(*,'(A)') '# lineCounts == 1 and dotP < 0 => CYCLE'
        IF(aveForce > lastAveForce .OR. aveForce>minAveForce) THEN
          cartPosition=cartOld
          totalForce=forceOld
        END IF
        timeStep=timeStep/2._DP
        velocity=0._DP
        linecounts=0
        CYCLE
      END IF

      IF (lineCounts > 1 .AND. maxForce > forceCut .AND. &
        (aveForce>minAveForce*1.5_DP .OR. aveForce>lastAveForce*1.25_DP &
          .OR. dotp < 0._DP)) THEN
        WRITE(*,'(A)') '# lineCounts > 1 and dotP < 0 => EXIT'
        IF(aveForce > minAveForce) THEN
          cartPosition=cartOld
          totalForce=forceOld
        END IF
        IF(lineCounts == 2) timeStep = timeStep *2._DP/3._DP
        EXIT
      END IF

      WRITE(*,'(A)') '# Good productive step => CONTINUE'
      ! This was a good, productive step, so we write out the restart 
      ! information to the restart file
      !CALL WriteRestartInformation(cartPosition, cell, atomicPSP, atomicSymbol)

      IF(timeStep < 1.E-6_DP) STOP "***MIN TIMESTEP REACHED***"

      WRITE(*,*) "maxForce=", maxForce
      WRITE(*,*) "forceCut=", forceCut

      ! Exit condition
      IF (maxForce <= forceCut) THEN
        WRITE(*,'(A,F15.10,A,F15.10,A)') '# maxForce <= forceCut (maxForce=', maxForce, ', forceCut=', forceCut, ') => EXIT'
        EXIT
      END IF

      ! If the search is going too slowly, speed it up a bit.
      IF(lineCounts > 3) THEN
        WRITE(*,'(A,F15.10,A,F15.10,A)') '# Too slow search (oldTimeStep=', timeStep, ', newTimeStep=', timeStep*2._DP, ') => BOOST'
        timeStep = timeStep*2._DP
      END IF
    END DO
    
    ! Exit completely if we found our solution!
    IF (maxForce <= forceCut) THEN
      WRITE(*,'(A,F15.10,A,F15.10,A)') '# maxForce <= forceCut (maxForce=', maxForce, ', forceCut=', forceCut, ') => DONE'
      EXIT
    END IF

  END DO 

END SUBROUTINE TransitionStateVelocityVerlot


SUBROUTINE TotalCalculation(numCalcs, cell, cartPosition, atomicSymbol, &
                            atomicPSP, trueForce, springForce, totalForce, &
                            energy, stressTensor, springConstant)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Finds all the forces (true force, spring force, and total) including the
!   effects of nudging and the climbing image, on the band of images.  Also
!   find the energy and the stress.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!

  INTEGER :: numCalcs
  REAL(KIND=DP), DIMENSION(:) :: springConstant
  REAL(KIND=DP), DIMENSION(3,3) :: cell
  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3) :: &
    trueForce, &
    springForce, &
    totalForce,&
    cartPosition

  REAL(KIND=DP), DIMENSION(numImages) :: energy
  CHARACTER(LEN=3), DIMENSION(numImages) :: atomicSymbol
  CHARACTER(LEN=*), DIMENSION(:) :: atomicPSP


                           !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3) :: massPosition
  REAL(KIND=DP), DIMENSION(numImages,3,3) :: stressTensor


                            !>> INITIALIZATION <<!

                             !>> FUNCTION BODY <<!

  CALL TrueCalculation(numcalcs, cell, cartPosition, atomicSymbol, &
                       atomicPSP, energy, trueForce, stressTensor)
  numCalcs = numCalcs + 1

!  WRITE(*,*) "CHECK: maxE=", MAXVAL(energy), ", E1=", energy(1)

  CALL Recenter(cell, cartPosition, atomicSymbol, massPosition)

  CALL GetSpringConstant(energy, trueForce, massPosition, springConstant)

  CALL GetTotalForce(cartPosition,trueForce,springForce,energy,springConstant, &
                     totalForce)

END SUBROUTINE TotalCalculation


SUBROUTINE TrueCalculation(numcalcs, cell, cartPosition, atomicSymbol, &
                           atomicPSP, energy, trueForce, stressTensor)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!

  INTEGER, INTENT(IN) :: numcalcs
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cell
  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3), INTENT(INOUT) :: cartPosition
  REAL(KIND=DP), DIMENSION(numImages,3,3), INTENT(OUT) :: stressTensor
  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3) :: trueForce
  REAL(KIND=DP), DIMENSION(numImages) :: energy
  CHARACTER(LEN=3), DIMENSION(numAtoms) :: atomicSymbol
  CHARACTER(LEN=80), DIMENSION(:) :: atomicPSP

                            !>> INTERNAL VARIABLES <<!

  CHARACTER(LEN=80) :: chjunk
  CHARACTER(80) :: line
  CHARACTER(LEN=9) :: imageName
  CHARACTER(LEN=19) :: dirImageName
  CHARACTER(LEN=4) :: sCalc, sCalc2
  INTEGER :: &
    ni, i, &
    stat, fs, &
    start_im, &             ! Start image
    end_im, &               ! End image
    sum, &                  ! sum of all the monitors
    status                  ! status of the monitor

  REAL(KIND=DP), DIMENSION(numAtoms,3) :: fp
  CHARACTER(LEN=4) :: label
  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3) :: fractPosition
  REAL(kind=DP), DIMENSION(3) :: tmpforce

                            !>> INITIALIZATION <<!
  fractPosition = CartToFract(cell, cartPosition)

                             !>> FUNCTION BODY <<!

  ! If this is the first calculation, we will do the first and last 
  ! images, as well as set up the working directory.
  IF (numcalcs == 0) THEN
    start_im=1
    end_im=numImages
    CALL SYSTEM("rm -rf IMAGE00*")
  ELSE
    start_im=2
    end_im=numImages-1
  END IF

  ! We first set up all the OFDFT files necessary
  DO ni = start_im, end_im

    WRITE(imageName,'(A5, I4.4)') "IMAGE", ni
    WRITE(dirImageName,'(A19)') imageName // "/" // imageName
    IF (numcalcs == 0) CALL SYSTEM("mkdir " // imageName)

    ! First write out the input files
    OPEN(UNIT=11, FILE="ORBITALFREE.inp", POSITION="REWIND", ACTION="READ")
    OPEN(UNIT=1, FILE=dirImageName//".inpt", POSITION="REWIND", &
         ACTION="WRITE", STATUS="REPLACE")

    DO
      READ(11,'(A80)',IOSTAT=stat) chjunk
      IF (stat/=0) EXIT
      WRITE(1,'(A80)') chjunk
    END DO
    CLOSE(11)

    IF(numCalcs == 0) THEN
      IF(initialConditions == 2 .AND. (ni == 1 .OR. ni == numImages)) THEN
        WRITE(1,'(A)') "MINI ION"
        WRITE(1,'(A)') "METHOD ION CG2"
        WRITE(1,'(A)') "TOLF 1e-4"
      END IF
    END IF

    CLOSE(1)

    ! Now write out the ion files
    OPEN(UNIT=1, FILE=dirImageName//".ion", POSITION="REWIND", &
         ACTION="WRITE",STATUS="REPLACE")

    WRITE(1,'(A)') "%BLOCK LATTICE_CART"
    WRITE(1,'(2X, 3F20.15)') cell
    WRITE(1,'(A)') "%END BLOCK LATTICE_CART"

    WRITE(1,'(A)') "%BLOCK POSITIONS_FRAC"
    DO i=1,numAtoms
      WRITE(1,'(2X, A2, 3F20.15)') atomicSymbol(i), fractPosition(ni,i,:)
    END DO
    WRITE(1,'(A)') "%END BLOCK POSITIONS_FRAC"

    WRITE(1,'(A)') "%BLOCK SPECIES_POT"
    DO i = 1, SIZE(atomicPSP)
      WRITE(1,'(2X, A)') atomicPSP(i)
    END DO
    WRITE(1,'(A)') "%END BLOCK SPECIES_POT"

    WRITE(1,'(A)') "%BLOCK ION_OPTIMIZATION"
    DO i=1,numAtoms
      WRITE(1,'(3I10)') frozenAtoms(i,:)
    END DO
    WRITE(1,'(A)') "%END BLOCK ION_OPTIMIZATION"

    CLOSE(1)

  END DO

  ! Now we jump through some hoops to submit the jobs to the queue
  ! We first set up all the OFDFT files necessary
  OPEN(UNIT=1, FILE="executeOFDFT.s", POSITION="REWIND", ACTION="WRITE", STATUS="REPLACE")
  WRITE(1,'(A)') "#!/bin/csh"
  DO ni = start_im, end_im
    WRITE(imageName,'(A5, I4.4)') "IMAGE", ni
    WRITE(sCalc,'(I4.4)') numcalcs
    IF (numcalcs < 1) THEN
      WRITE(sCalc2,'(I4.4)') numcalcs
    ELSE
      WRITE(sCalc2,'(I4.4)') numcalcs-1
    END IF
    WRITE(1,*) "cd " // imageName

    IF(numcalcs>0) THEN ! mohan add 2014-10-07
      WRITE(1,*) "mv -f " // imageName // ".trans " // imageName // "." // sCalc2 // ".trans"
      WRITE(1,*) "mv -f " // imageName // ".force.out " // imageName // "." // sCalc2 // ".force.out"
    ENDIF

    WRITE(1,*) "rm -rf " // imageName // ".force.out "
    WRITE(1,*) "cp -f " // imageName // ".ion " // imageName // "." // sCalc // ".ion"
    WRITE(1,*) "cp ../*.recpot ."
    !WRITE(1,*) "cp ../*.den ."
    WRITE(1,*) "subofdftp_cineb " // imageName
!
! PLEASE WRITE DOWN THE PATH OF PROFESS EXEUCTABLE
!
    WRITE(1,*) "../PROFESS " // imageName // " > " // imageName // ".log"
    WRITE(1,*) "cd .."
  END DO
  CLOSE(1)
  CALL SYSTEM("chmod u+x executeOFDFT.s;./executeOFDFT.s")

  ! Now we monitor the running OFDFT jobs until they are all done.
  DO
    OPEN(UNIT=11,FILE="./MONITOR.dat")
    WRITE(11,'(A5,A10)') "IMAGE","STATUS"
    sum = 0
    DO ni=1,numImages
      ! Check for the .trans file
      WRITE(imageName,'(A5, I4.4)') "IMAGE", ni
      WRITE(dirImageName,'(A19)') imageName // "/" // imageName
      OPEN(UNIT=1,ACTION="READ",STATUS="OLD",FILE=dirImageName//".force.out",IOSTAT=stat)
      CLOSE(1)
      ! File was successfully opened, exists!
      IF (stat==0) THEN
        WRITE(11,'(I5,A10)') ni,"DONE"
      ELSE
        WRITE(11,'(I5,A10)') ni,"RUNNING"
      END IF
      sum=sum+stat
    END DO
    CLOSE(11)
    IF (sum==0) EXIT
    WRITE(*,*) "wait 15 seconds to harvest data from .trans files"
    CALL SLEEP(15) 
  END DO

  ! Now we harvest the data output from the .trans files
  DO ni=1,numImages
    WRITE(imageName,'(A5, I4.4)') "IMAGE", ni
    WRITE(dirImageName,'(A19)') imageName // "/" // imageName
    ! Open the interface file, get the energy, forces, and stress
    OPEN(UNIT=1, ACTION="READ",STATUS="OLD",FILE=dirImageName // ".trans")
    ! Read energy
    READ(1,*) energy(ni)
    WRITE(*,*) "CHECK: ", energy(ni)
    ! Read force
    DO i = 1, numAtoms
      READ(1,*) tmpforce(1), tmpforce(2), tmpforce(3)
      IF ( frozenAtoms(i,1)==0 ) tmpforce(1) = 0._DP
      IF ( frozenAtoms(i,2)==0 ) tmpforce(2) = 0._DP
      IF ( frozenAtoms(i,3)==0 ) tmpforce(3) = 0._DP
!      WRITE(*,*) "CHECK: tmpforce=(", tmpforce(1), ", ", tmpforce(2), ", ", tmpforce(3), ")"
      trueForce(ni,i,:) = tmpforce(:)
    END DO

    ! Note that this is actually the transpose, but oh well.
    READ(1,*) stressTensor(ni,:,:)

    ! If this is the first or last image, and we did a geometry optimization
    ! Read the relaxed ion positions in
    IF(initialConditions == 2 .AND. (ni == 1 .OR. ni == numImages)) THEN

      ! Convert to fractional coordinates first, the whole thing ... this is a 
      ! bit messy and more work than required.  You can fix it if you want.
      cartPosition = CartToFract(cell, cartPosition)

      ! Read in the new relaxed fractional positions
      DO i = 1, numAtoms
        READ(1,*) cartPosition(ni,i,:)
      END DO

      ! Convert back to cartesian coordinates
      cartPosition = FractToCart(cell, cartPosition)
          
    END IF

    CLOSE(1)

    ! Copy the density file to the next start density file
    !IF(ni>1 .AND. ni<numImages) THEN
    !  CALL SYSTEM("mv " // dirImageName // ".0.den " // imageName // "/" // "startDensity.den")
    !END IF

  END DO

END SUBROUTINE TrueCalculation


SUBROUTINE GetSpringConstant(energy, trueForce, massPosition, springConstant)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Finds an appropriate spring constant to use in the climbing-image
!   nudged elastic band
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   We may want to implement variable spring constants, as described in
!   JCP 113, pg 9901   
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    trueForce, &           ! The true force
    massPosition           ! the center of mass positions of the ions

  REAL(kind=DP), DIMENSION(:), INTENT(OUT) :: &
    springConstant         ! The spring constant we're trying to get

  REAL(kind=DP), DIMENSION(:), INTENT(IN) :: &
    energy

                            !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    ni                     ! counter for # of images

  REAL(kind=DP) :: &
    emax, &
    eref, &
    ei, &
    kmax, &
    dk, &
    maxForce, &           ! Max true force over all ions in each image
    minDis                ! Max displacement over all ions in each images

                            !>> INITIALIZATION <<!

  ! Initialize the spring constant
  springConstant=1._DP

                             !>> FUNCTION BODY <<!

  ! Figure out the maximum true force (over all ions) over all the images.
  maxForce=0._DP
  DO ni=1, SIZE(trueForce,1)
    maxForce = MAX(SUM(trueForce(ni,:,:)**2), maxForce)
  END DO
  maxForce = SQRT(maxForce)

  ! Figure out the maximum displacement (over all ions) on any image.
  minDis=100000000000._DP             ! Initialized to large number
  DO ni=1, SIZE(massPosition,1)-1
    minDis = MIN(SUM((massPosition(ni,:,:) &
                           - massPosition(ni+1,:,:))**2), minDis)
  END DO
  minDis = SQRT(minDis)

  ! Calculate the spring constant
  kmax = maxForce / minDis
  dk = kmax / 2._DP

  emax = MAXVAL(energy)
  eref = MAX(energy(1), energy(SIZE(energy)))


  DO ni=1, SIZE(springConstant)-1
    ei = MAX(energy(ni), energy(ni+1))
    IF(ei > eref) THEN
      springConstant(ni) = kmax - dk*(emax - ei)/(emax-eref)
    ELSE
      springConstant(ni) = kmax - dk
    END IF
  END DO


!  WRITE(*,*) springConstant
!  WRITE(*,*)
END SUBROUTINE GetSpringConstant


SUBROUTINE GetTotalForce(positions, trueForces, springForces, energies, springConstant, & 
                         totalForces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Finds the total force depending on the type of run (nudgeType).
!
! GLOBAL/MODULE VARIABLES CHANGES:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                         !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    positions                 ! Positions of ions

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(INOUT) :: &
    trueForces                ! Forces on images/ions due to energy. OUT: perp. component only

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: &
    energies                  ! Energies of the images

  REAL(kind=DP), INTENT(IN), DIMENSION(:) :: &
    springConstant

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: &
    springForces, &
    totalForces

  REAL(kind=DP), DIMENSION(3) :: &
    tmpforce

                         !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(SIZE(positions,1),SIZE(positions,2), &
                           SIZE(positions,3)) :: &
    tangents                  ! Tangents of each image/ion

  INTEGER :: &
    numImages, &
    i, ni
  
                         !>> INITIALIZATION <<!

  numImages = SIZE(positions,1)
  springForces = 0._DP
  totalForces = 0._DP

                         !>> FUNCTION BODY <<!

  ! NO BAND
  IF (nudgeType==-1) THEN
    totalForces = trueForces

  ! NEB or CINEB using the "improved tangent method
  ELSE IF ((nudgeType == 1).OR.(nudgeType == 2)) THEN
    CALL HenkelTangent(positions, energies, tangents)
    CALL HenkelForce(positions, tangents, trueForces, springForces, energies, &
                     springConstant, totalForces)

  ! (PLAIN) ELASTIC BAND
  ELSE
    springForces(2:numImages-1,:,:) = springConstant(1) &
                                       * (  positions(3:numImages,:,:) &
                                          + positions(1:numImages-2,:,:) &
                                          - 2._DP*positions(2:numImages-1,:,:))
    IF (nudgeType==0) THEN
      totalForces = trueForces + springForces

    ! NEB or CINEB using a spline fit for tangent
    ELSE ! nudgeTypes 3,4,5,6
      CALL SplineTangent(positions, tangents)
      CALL SplineForce(tangents, trueForces, springForces, energies, &
                       totalForces)
    END IF
  END IF

  ! Ensure that the forces on the initial and final images are 0
  totalForces(1,:,:) = 0._DP
  totalForces(numImages,:,:) = 0._DP

  DO ni=1,numImages
    DO i = 1, numAtoms
      tmpforce(:) = totalForces(ni,i,:)
      IF ( frozenAtoms(i,1)==0 ) tmpforce(1) = 0._DP
      IF ( frozenAtoms(i,2)==0 ) tmpforce(2) = 0._DP
      IF ( frozenAtoms(i,3)==0 ) tmpforce(3) = 0._DP
!      WRITE(*,*) "CHECK: tmpforce=(", tmpforce(1), ", ", tmpforce(2), ", ", tmpforce(3), "), numAtoms=", numAtoms, ", numImages=", numImages
      totalForces(ni,i,:) = tmpforce(:)
    END DO
  END DO

END SUBROUTINE GetTotalForce


SUBROUTINE HenkelTangent(positions, energies, tangents)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine calculates the tangent defined by Henkelman et. al. at each
!   image.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1]  J. Chem. Phys. Vol 113, No. 22, 8 December 2000, Pg. 9978.
!        "Improved tangent estimate in the nudged elastic band method for 
!         minimum energy paths and saddle points."
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2/16/2005  Subroutine Created (Greg)
!
!------------------------------------------------------------------------------  
  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    positions                         ! The ion positions in the lattice

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: &
    energies                          ! Energies of the images

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: &
    tangents                          ! Normalized tangents


                         !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    i, &                              ! a counter
    numImages, &                      ! number of images
    numIons                           ! number of ions

  REAL(kind=DP), DIMENSION(SIZE(positions,2), SIZE(positions,3)) :: &
    t_plus, t_minus                   ! See reference for definitions

  REAL(kind=DP) :: &
    v_max, v_min, &                   ! See reference for definitions
    norm                              ! Un-normalized length of tangent vector

                         !>> INITIALIZATION <<!
  numImages = SIZE(positions, 1)
  numIons = SIZE(positions, 2)
  tangents = 0._DP

! Tangent defined as 0 at 1st and last images
! Calculation for other tangents below.

                         !>> FUNCTION BODY <<!
  DO i = 2, numImages-1

    ! Calculate t+ (and t-), secant between i and i+1 (and i-1) images
    t_plus = positions(i+1,:,:) - positions(i,:,:)
    t_minus = positions(i,:,:) - positions(i-1,:,:)

    ! Set the secant toward higher energy adjacent image as the "tangent"
    IF (energies(i+1) > energies(i) .AND. energies(i) > energies(i-1)) THEN
      tangents(i,:,:) = t_plus

    ELSE IF(energies(i+1) < energies(i) .AND. energies(i) < energies(i-1)) THEN
      tangents(i,:,:) = t_minus

    ELSE ! If image i at max/min, use weighted average to determine "tangent"
      ! Calculate v_max and v_min (weighting depends on energy)
      v_max = MAX(ABS(energies(i+1) - energies(i)), &
                  ABS(energies(i-1) - energies(i)))
      v_min = MIN(ABS(energies(i+1) - energies(i)), &
                  ABS(energies(i-1) - energies(i)))
      IF(energies(i+1)>energies(i-1)) THEN
        tangents(i,:,:) = t_plus*v_max + t_minus*v_min
      ELSE
        tangents(i,:,:) = t_plus*v_min + t_minus*v_max
      END IF

    END IF

    ! Normalize the tangent of image
    norm = SQRT(SUM(tangents(i,:,:)**2))
    tangents(i,:,:) = tangents(i,:,:) / norm

  END DO

END SUBROUTINE HenkelTangent


SUBROUTINE HenkelForce(positions, tangents, trueForces, springForces, energies, &
                       springConstant, totalForces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   CINEB and NEB total forces on images are calculated using the "improved
!   tangent method by Henkelman et al.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Implement variable length spring constants
!
! REFERENCES:
!   [1]  J. Chem. Phys. Vol 113, No. 22, 8 December 2000, Pg. 9978.
!        "Improved tangent estimate in the nudged elastic band method for 
!         minimum energy paths and saddle points."
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2/16/2005  Subroutine Created (Greg)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                         !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    positions, &                      ! Position of each ion
    tangents                          ! Normalized tangents

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: &
    energies                          ! Energies of the images

  REAL(kind=DP), DIMENSION(:), INTENT(IN) :: &
    springConstant

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: &
    totalForces, &                    ! The force on (ions in) each image
    springForces

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(INOUT) :: &
    trueForces                        ! Force due to energy of image, when out,
                                      ! the perperdicular direction


                         !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(SIZE(positions,1),SIZE(positions,2), &
                           SIZE(positions,3)) :: &
    springForceParallel, &            ! Spring force in tangent direction
    trueForceParallel                 ! True force in tangent direction

  REAL(kind=DP), DIMENSION(SIZE(positions,2), SIZE(positions,3)) :: &
    r_upper, r_lower

  INTEGER :: &
    numImages, &
    i
  
  INTEGER, DIMENSION(1) :: &
    highImage                         ! Index of the highest-energy image

                         !>> INITIALIZATION <<!

  springForceParallel = 0._DP
  trueForceParallel = 0._DP
  numImages = SIZE(positions,1)


                         !>> FUNCTION BODY <<!

  ! Equation 12 in the above reference
  DO i = 2, numImages - 1
    r_upper = positions(i+1,:,:) - positions(i,:,:)
    r_lower = positions(i,:,:) - positions(i-1,:,:)

    springForceParallel(i,:,:) = tangents(i,:,:) * (SQRT(SUM(r_upper**2))*springConstant(i) - SQRT(SUM(r_lower**2)) * springConstant(i-1))

    trueForceParallel(i,:,:) = tangents(i,:,:) * &
                               SUM(trueForces(i,:,:) * tangents(i,:,:))
  END DO

  totalForces = (trueForces - trueForceParallel) + springForceParallel

  ! The highest image climbs in energy for CINEB
  IF (nudgeType == 2) THEN
    highImage = MAXLOC(energies(2:numImages-1)) ! Gives location in subarray
    highImage = highImage + 1                   ! Location in energies array
    totalForces(highImage,:,:) = trueForces(highImage,:,:) &
                                - 2._DP * trueForceParallel(highImage,:,:)
  END IF

  springForces = springForceParallel
  trueForces = trueForces - trueForceParallel
  
END SUBROUTINE HenkelForce


SUBROUTINE SplineTangent(positions, tangents)
!------------------------------------------------------------------------------
! DESCRIPTION: Calculates the tangents of ions at each image used in
!              spline-fitting.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------  
  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    positions                ! Cartesian positions of ions

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: &
    tangents                 ! Normed tangents corresponding to ions/images

                         !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    i, j, ni, &              ! Counters
    numImages, &             ! Number of images
    numAtoms                 ! Number of ions per image

  REAL(KIND=DP) :: &
    spacing, &               ! The distance between two images on the band
    a, b, c, &               ! Spline coefficients (ax^3+bx^2+cx+d)
    x, & 
    norm                     ! Norm of tangents per image

  REAL(KIND=DP), DIMENSION(SIZE(positions,1)) :: &
    x_values, &              ! The position of the mth image on the band
                             ! (distance term)
    y_values                 ! Position of an ion corresponding to image

  REAL(KIND=DP), DIMENSION(4*(SIZE(positions,1)-1)) :: &
    spline_fit

                         !>> INITIALIZATION <<!

  numImages = SIZE(positions,1)
  numAtoms = SIZE(positions,2)
  x_values(1)=0._DP

                         !>> FUNCTION BODY <<

  ! Initialize the x-values (distance of the images from one another on band)
  DO ni=2,numImages
    spacing = SQRT(SUM((positions(ni,:,:)-positions(ni-1,:,:))**2))
    x_values(ni) = x_values(ni-1) + spacing
  END DO
  
  DO i=1,numAtoms
    DO j=1,3
      
      y_values(:)=positions(:,i,j)

      spline_fit=0._DP
      IF (nudgeType > 4) THEN ! nudgeType 5 or 6
        CALL Spline2(x_values, y_values, spline_fit)
      ELSE ! nudgeType 3 or 4
        CALL Spline(x_values, y_values, spline_fit)
      END IF
      DO ni=1,numImages-1
        x=x_values(ni)
        a=spline_fit(4*ni-3)     ! x^3 coefficient of spline fit
        b=spline_fit(4*ni-2)     ! x^2 coefficient of spline fit
        c=spline_fit(4*ni-1)     ! x coefficient of spline fit
        tangents(ni,i,j)=(3._DP*a*(x**2))+(2._DP*b*x)+c
      END DO

      ! Don't forget the last image's tangents
      x=x_values(numImages)
      tangents(numImages,i,j)=(3._DP*a*(x**2))+(2._DP*b*x)+c

    END DO
  END DO

  ! Normalize the tangents within each image
  DO ni=1,numImages
    norm = SQRT(SUM(tangents(ni,:,:)**2))
    tangents(ni,:,:) = tangents(ni,:,:) / norm
  END DO

END SUBROUTINE SplineTangent



SUBROUTINE SplineForce(tangents, trueForces, springForces, energies, & 
                       totalForces)
!------------------------------------------------------------------------------
! DESCRIPTION: Gives the total force for CINEB and NEB using the spline
!              tangents.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                            !>> EXTERNAL VARIABLES <<!


  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    tangents, &          ! Normalized tangents
    trueForces, &        ! Force from energy of image
    springForces         ! Force from "spring"

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: &
    energies             ! Energy of each image

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: &
    totalForces          ! Overall force acting on image


                            !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    i, j, ni, &          ! Counters
    numImages, &
    numAtoms

  INTEGER,DIMENSION(1) :: &
    highImage            ! Index of highest-energy image

  REAL(KIND=DP), DIMENSION(SIZE(tangents,1),SIZE(tangents,2), &
                           SIZE(tangents,3)) :: &
    trueparallel, &      ! Component of true force parallel to tangent
    springparallel       ! Component of spring force parallel to tangent


                            !>> INITIALIZATION <<!

  numImages = SIZE(tangents,1)
  numAtoms = SIZE(tangents,2)

                            !>> FUNCTION BODY <<!

  DO ni=1,numImages

    ! Find projections of true and spring forces parallel to tangen
    trueparallel(ni,:,:)=trueForces(ni,:,:)*tangents(ni,:,:)
    springparallel(ni,:,:)=springForces(ni,:,:)*tangents(ni,:,:)

  END DO

  ! Nudge!
  totalForces = (trueForces - trueparallel) + springparallel

  ! Loop executes for CINEB
  IF ((nudgeType==4).OR.(nudgeType==6)) THEN
    highImage = MAXLOC(energies(2:numImages-1)) ! Gives index in subarray
    highImage = highImage + 1                   ! Index in energies array
    totalForces(highImage,:,:) = trueForces(highImage,:,:) &
                               -2._DP*trueparallel(highImage,:,:)
  END IF

END SUBROUTINE SplineForce



SUBROUTINE Spline(xvalues,yvalues,spfit)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: &
    xvalues, &
    yvalues

  REAL(KIND=DP), DIMENSION(:), INTENT(OUT) :: &
    spfit


                         !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    numvalues, &
    pos, &
    nv, &
    i

  REAL(KIND=DP), DIMENSION(SIZE(spfit)+1,SIZE(spfit)) :: &
    fitlat

  REAL(KIND=DP), DIMENSION(1,SIZE(spfit)) :: &
    gauss_fit

  REAL(KIND=DP), DIMENSION(SIZE(spfit),SIZE(spfit)) :: &
    gauss_lat

                         !>> INITIALIZATION <<!
  numvalues = SIZE(xvalues)

                         !>> FUNCTION BODY <<!


  !-- construct the fit matrix --
  fitlat=0._DP
  pos=0
  DO nv=1,numvalues-1
    i=(4*(nv-1))+1

    pos=pos+1
    fitlat(i,pos)=xvalues(nv)**3
    fitlat(i+1,pos)=xvalues(nv)**2
    fitlat(i+2,pos)=xvalues(nv)
    fitlat(i+3,pos)=1._DP
    fitlat(SIZE(spfit)+1,pos)=yvalues(nv)

    pos=pos+1
    fitlat(i,pos)=xvalues(nv+1)**3
    fitlat(i+1,pos)=xvalues(nv+1)**2
    fitlat(i+2,pos)=xvalues(nv+1)
    fitlat(i+3,pos)=1._DP
    fitlat(SIZE(spfit)+1,pos)=yvalues(nv+1)
  END DO

  DO nv=2,numvalues-1
    i=(4*(nv-2))+1

    pos=pos+1
    fitlat(i,pos)=3._DP*xvalues(nv)**2
    fitlat(i+1,pos)=2._DP*xvalues(nv)
    fitlat(i+2,pos)=1._DP
    fitlat(i+4,pos)=-3._DP*xvalues(nv)**2
    fitlat(i+5,pos)=-2._DP*xvalues(nv)
    fitlat(i+6,pos)=-1._DP

    pos=pos+1
    fitlat(i,pos)=6._DP*xvalues(nv)
    fitlat(i+1,pos)=2._DP
    fitlat(i+4,pos)=-6._DP*xvalues(nv)
    fitlat(i+5,pos)=-2._DP
  END DO

  pos=pos+1
  fitlat(1,pos)=1._DP
  pos=pos+1
  fitlat(4*(numvalues-2)+1,pos)=1._DP
  !-----------------------------

  gauss_lat(:,:)=fitlat(1:SIZE(spfit),:)

  DO i=1,SIZE(spfit)
    gauss_fit(1,i)=fitlat(SIZE(spfit)+1,i)
  END DO

  CALL GaussElimination(SIZE(spfit),SIZE(spfit),gauss_lat,1,gauss_fit)

  DO i=1,SIZE(spfit)
    spfit(i)=gauss_fit(1,i)
  END DO

END SUBROUTINE Spline


SUBROUTINE Spline2(x,y,spfit)
!-----------------------------------------------------------------------------
! DESCRIPTION: Gives the coefficients for a cubic spline fit going through
!              data points (x,y) with natural boundary conditions.  Fit is:
!                S(k) = a(k)(x-x_k)^3 + b(k)(x-x_k)^2 + c(k)(x-x_k) + d(k)
!              Solves using Crout factorization of a tri-diagonal system.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPTROVEMENTS:
!
! REFERENCES: Burden and Faires, Numerical Analysis.
!
!------------------------------------------------------------------------------
! REVISION LOG: (4/24/06) Created by Linda Hung to replace previous spline-
!               fitting subroutine.
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                           !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: &
    x, &
    y

  REAL(KIND=DP), DIMENSION(:), INTENT(OUT) :: &
    spfit


                           !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    k, &
    numvalues

  REAL(KIND=DP), DIMENSION(SIZE(x)-1) :: &
    h, &        !Difference between two adjacent x-values
    alpha, &
    mu, &
    a, &
    c

  REAL(KIND=DP), DIMENSION(SIZE(x)) :: &
    L, &
    z, &
    b


                           !>> INITIALIZATION <<!

  numvalues = SIZE(x)

                           !>> FUNCTION BODY <<!

  DO k=1,numvalues-1
    h(k) = x(k+1)-x(k)
  END DO

  DO k = 2,numvalues-1
    alpha(k) = 3._DP*(y(k+1)-y(k))/h(k) - 3._DP*(y(k)-y(k-1))/h(k-1)
  END DO

  ! Solve a tri-diagonal linear system using Crout factorization

  L(1)=1._DP
  mu(1)=-1._DP
  z(1)=0._DP

  DO k=2,numvalues-1
    L(k) = 2._DP*(x(k+1)-x(k-1)) - h(k-1)*mu(k-1)
    mu(k) = h(k)/L(k)
    z(k) = (alpha(k)-h(k-1)*z(k-1))/L(k)
  END DO

  L(numvalues) = -1._DP-mu(numvalues-1)
  z(numvalues) = z(numvalues-1)/L(numvalues)
  b(numvalues) = z(numvalues)

  ! Obtain values for spline coefficients and put into spfit vector

  DO k=numvalues-1,1,-1
    b(k)=z(k)-mu(k)*b(k+1)
    c(k)=(y(k+1)-y(k))/h(k) - h(k)*(b(k+1)+2._DP*b(k))/3._DP
    a(k)=(b(k+1)-b(k))/(3._DP*h(k))
  END DO

  DO k=1,numvalues-1
    spfit(4*k-3)=a(k)
    spfit(4*k-2)=b(k)
    spfit(4*k-1)=c(k)
    spfit(4*k)=y(k)
  END DO

END SUBROUTINE Spline2


SUBROUTINE Recenter(cell, cartPosition, atomicSymbol, massPosition)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine outputs massPosition, or the position in center-of-mass
!   coordinates.  We need to use these to prevent net translation/rotation
!   of the images.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   Czerminski and Elber, "Reaction path study of conformational transitions
!     in Flexible systems:  Applications to peptides."  J Chem Phys 92 (9), 
!     1 May 1990.  Pg 5580 (pay attention to eq. 3)
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2/17/2005 Created by Greg
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!


  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cell
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: cartPosition
  REAL(KIND=DP), DIMENSION(SIZE(cartPosition,1), SIZE(cartPosition, 2), &
                           SIZE(cartPosition,3)), INTENT(OUT) :: massPosition
  CHARACTER(LEN=3), DIMENSION(SIZE(cartPosition,2)), INTENT(IN) :: atomicSymbol

                          !>> INTERNAL VARIABLES <<! 
  INTEGER :: &
    i, &
    na, &
    m, &
    ni, &
    numAtoms, &
    numImages

  REAL(KIND=DP), DIMENSION(3,3) :: invcell
  REAL(kind=DP), DIMENSION(3) :: &
    center, &             
    theta        ! The angles we need to rotate by

  REAL(KIND=DP) :: totalMass
  REAL(kind=DP), DIMENSION(SIZE(cartPosition,2)) :: mass

                            !>> INITIALIZATION <<!
  IF(SIZE(cartPosition,3) /= 3) STOP "**WRONG SIZE FOR CARTPOSITION, 3**"

  numImages = SIZE(cartPosition,1)   
  numAtoms = SIZE(cartPosition,2)

  DO i = 1, numAtoms
    mass(i) = AtomicMass(atomicSymbol(i))
  END DO

  totalMass=SUM(mass)

                            ! >> FUNCTION BODY <<!

  DO ni=1,numImages

    ! Determine center of mass
    center = 0._DP
    DO na=1,numAtoms
      center = center + mass(na) * cartPosition(ni, na, :)
    END DO

    center = center / totalMass

    ! Offset everything to center of mass coordinates
    DO na=1,numAtoms
      massPosition(ni,na,:) = cartPosition(ni,na,:) - center
    END DO

    ! Now take care of rotation.
    theta = 0._DP
    DO na = 1, numAtoms
      theta = theta + (mass(na) * Cross(massPosition(ni, na, :), &
                                        massPosition(1, na, :)))
    END DO

    theta = theta / totalMass

    ! Offset everything to the right rotation
    DO na=1, numAtoms
!      CALL Rotate(massPosition(ni,na,:), -theta)
    END DO

  END DO

END SUBROUTINE Recenter


SUBROUTINE Rotate(vector, theta)
!------------------------------------------------------------------------------
! DESCRIPTION: 
!   Rotates a vector by theta(1) around the x axis, theta(2) around the y axis,
!   theta(3) around the z axis.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2/17/2005 Created by Greg 
!
!------------------------------------------------------------------------------
  USE MATHFUNCTIONS, ONLY : &
    Vecmul ! Multiplies a matrix and a vector

  IMPLICIT NONE
                            !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(3), INTENT(INOUT) :: &
    vector

  REAL(kind=DP), DIMENSION(3) :: &
    theta

                            !>> INTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(3,3) :: &
    x_rot, y_rot, z_rot

                             !>> INITIALIZATION <<!
                             !>> FUNCTION BODY <<!

  x_rot(1,1) = 1._DP
  x_rot(1,2) = 0._DP
  x_rot(1,3) = 0._DP
  x_rot(2,1) = 0._DP
  x_rot(2,2) = COS(theta(1))
  x_rot(2,3) = -SIN(theta(1))
  x_rot(3,1) = 0._DP
  x_rot(3,2) = SIN(theta(1))
  x_rot(3,3) = COS(theta(1))

  y_rot(1,1) = COS(theta(2))
  y_rot(1,2) = 0._DP
  y_rot(1,3) = SIN(theta(2))
  y_rot(2,1) = 0._DP
  y_rot(2,2) = 1._DP
  y_rot(2,3) = 0._DP
  y_rot(3,1) = -SIN(theta(2))
  y_rot(3,2) = 0._DP
  y_rot(3,3) = COS(theta(2))

  z_rot(1,1) = COS(theta(3))
  z_rot(1,2) = -SIN(theta(3))
  z_rot(1,3) = 0._DP
  z_rot(2,1) = SIN(theta(3))
  z_rot(2,2) = COS(theta(3))
  z_rot(2,3) = 0._DP
  z_rot(3,1) = 0._DP
  z_rot(3,2) = 0._DP
  z_rot(3,3) = 1._DP

  vector = Vecmul(x_rot, Vecmul(y_rot, Vecmul(z_rot, vector)))

END SUBROUTINE Rotate


FUNCTION CartToFract(cell,cartPosition)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                            !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: &
    cell

  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3), INTENT(IN) :: &
    cartPosition

  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3) :: &
    CartToFract

                            !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    na,m,ni

  REAL(KIND=DP), DIMENSION(3,3) :: &
    invCell

                             !>> INITIALIZATION <<!
  invCell = Inverse(cell)

                             !>> FUNCTION BODY <<!

  DO ni=1,numImages
    DO na=1,numAtoms
      DO m=1,3
        CartToFract(ni,na,m)=(invCell(1,m)*cartPosition(ni,na,1)) +&
                             (invCell(2,m)*cartPosition(ni,na,2)) +&
                             (invCell(3,m)*cartPosition(ni,na,3))
      END DO
    END DO
  END DO

END FUNCTION CartToFract


FUNCTION FractToCart(cell, fractPosition)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function converts the positions of the ions in all the images from
!   fractional to cartesian coordinates.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/21/06 Converted to a function from a subroutine.
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                            !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: &
    cell

  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3), INTENT(IN) :: &
    fractPosition

  REAL(KIND=DP), DIMENSION(numImages, numAtoms, 3) :: &
    FractToCart

                            !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    na,m,ni

                             !>> INITIALIZATION <<!
                             !>> FUNCTION BODY <<!
  DO ni=1,numImages
    DO na=1,numAtoms
      DO m=1,3
        FractToCart(ni,na,m)=(cell(1,m)*fractPosition(ni,na,1))+&
                             (cell(2,m)*fractPosition(ni,na,2))+&
                             (cell(3,m)*fractPosition(ni,na,3))
      END DO
    END DO
  END DO

END FUNCTION FractToCart


SUBROUTINE GaussElimination(NR,NCLHS,LEFT,NCRHS,RIGHT)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gaussian elimination with partial pivoting
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!
  INTEGER :: &
    NR, &
    NCLHS, &
    NCRHS

  REAL(KIND=DP), DIMENSION(NCLHS,NR) :: &
    LEFT

  REAL(KIND=DP), DIMENSION(NCRHS,NR) :: &
    RIGHT

                         !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    i, j, m, &
    ind_piv

  INTEGER, DIMENSION(NR) :: &
    IPIV

  REAL(KIND=DP) :: &
    max_piv, &
    minv,maxv

  REAL(KIND=DP), DIMENSION(NCLHS+NCRHS,NR) :: &
    worklat

  REAL(KIND=DP), DIMENSION(NCLHS+NCRHS) :: &
    work2

                         !>> INITIALIZATION <<!
                         !>> FUNCTION BODY <<!

  worklat(1:NCLHS,:)=LEFT
  worklat(NCLHS+1:NCLHS+NCRHS,:)=RIGHT

  DO i=1,NR
    IPIV(i)=i
  END DO

  DO i=1,NR

    !----- PIVOTING SECTION -----
    DO j=i,NR
      minv=MINVAL(worklat(i:NCLHS,j))
      maxv=MAXVAL(worklat(i:NCLHS,j))
      minv=ABS(minv)
      maxv=ABS(maxv)
      IF (minv>maxv) THEN
        worklat(:,j)=worklat(:,j)/minv
      ELSE
        worklat(:,j)=worklat(:,j)/maxv
      END IF
      IF (worklat(i,j)<0._DP) THEN
        worklat(:,j)=-1._DP*worklat(:,j)
      END IF
    END DO

    max_piv=0._DP
    ind_piv=j

    DO j=i,NR
      IF (worklat(i,j)>max_piv) THEN
        max_piv=worklat(i,j)
        ind_piv=j
      END IF
    END DO

    work2=worklat(:,i)
    worklat(:,i)=worklat(:,ind_piv)
    worklat(:,ind_piv)=work2(:)
    m=IPIV(i)
    IPIV(i)=ind_piv
    IPIV(ind_piv)=m
    !----------------------------

    !----- ELIMINATION SECTION -----
    worklat(:,i)=worklat(:,i)/worklat(i,i)
    DO j=1,NR
      IF (j/=i) THEN
        IF (worklat(i,j)/=0._DP) THEN
          worklat(:,j)=worklat(:,j)/worklat(i,j)
          worklat(:,j)=worklat(:,j)-worklat(:,i)
        END IF
      END IF
    END DO
    !-------------------------------

  END DO

  DO i=1,NR
    worklat(:,i)=worklat(:,i)/worklat(i,i)
  END DO

  m=0
  IF (m==1) THEN
    DO j=1,NR
      IF (IPIV(j)/=j) THEN

        DO m=1,NR
          IF (IPIV(m)==j) THEN
            ind_piv=m
          END IF
        END DO

        work2=worklat(:,j)
        worklat(:,j)=worklat(:,ind_piv)
        worklat(:,ind_piv)=work2(:)
        m=IPIV(j)
        IPIV(j)=IPIV(ind_piv)
        IPIV(ind_piv)=j

      END IF
    END DO
  END IF

  RIGHT=worklat(NCLHS+1:NCLHS+NCRHS,:)

END SUBROUTINE GaussElimination


SUBROUTINE OutputWriteHeader(cell)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Writes to the output file for the transition state finder
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/18/2005 Commented and cleaned up the output file (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                           !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3,3) :: &
    cell              ! the cell vectors

                           !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    n                 ! counter

                            !>> INITIALIZATION <<!
                             !>> FUNCTION BODY <<!

  WRITE(transOutputUnit,'(A5,A45)') '','------------------- CELL --------------------'
  DO n=1,3
    WRITE(transOutputUnit,'(A5,3F15.10)') '',cell(1,n),cell(2,n),cell(3,n)
  END DO
  WRITE(transOutputUnit,'(A5,A45)') '','---------------------------------------------'
  WRITE(transOutputUnit,*) ''

  WRITE(transOutputUnit,'(A)') "                     ENERGY            ENERGY           MAXIMUM            NORM                     "
  WRITE(transOutputUnit,'(A)') " CALC LINE STEP   BARRIER 1 (eV)    BARRIER 2 (eV)  TOTAL FORCE (eV/A)  OF FORCES(eV/A)     TIMESTEP"
  WRITE(transOutputUnit,'(A)') "----------------------------------------------------------------------------------------------------"
  WRITE(*, '(A)') "   RUNS LNSR STEP   E_BARRIER 1(eV)   E_BARRIER 2(eV)   MAX FORCE(eV/A)   FORCE NORM(eV/A)  TIMESTEP CHECK"

!  WRITE(transImagesOutputUnit,'(A17,A45)') '','------------------- CELL --------------------'
!  DO n=1,3
!    WRITE(transImagesOutputUnit,'(A17,3F15.10)') '',cell(1,n),cell(2,n),cell(3,n)
!  END DO
!  WRITE(transImagesOutputUnit,'(A17,A45)') '','---------------------------------------------'
!  WRITE(transImagesOutputUnit,*) ''

  ! Flush the write buffer
  flush(transOutputUnit)

END SUBROUTINE OutputWriteHeader


SUBROUTINE OutputWrite(cartPosition, massPosition, energy, trueForce, &
                       springForce, totalForce, stressTensor, numcalcs, &
                       numiterations, linecounts, badStep, timeStep)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Writes to the output file for the transition state finder
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/18/2005 Commented and cleaned up the output file (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                           !>> EXTERNAL VARIABLES <<!
  INTEGER :: &
    numCalcs, &
    numIterations, &
    lineCounts

  REAL(kind=DP), INTENT(IN) :: &
    timeStep

  REAL(KIND=DP), DIMENSION(numImages) :: &
    energy

  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3) :: &
    trueForce, &
    totalForce, &
    springForce, &
    cartPosition, &
    massPosition

  REAL(KIND=DP), DIMENSION(numImages,3,3) :: &
    stressTensor

  LOGICAL, INTENT(IN) :: &
    badStep

                           !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    ni,na,m,n

  REAL(KIND=DP) :: &
    dmag

  REAL(KIND=DP), DIMENSION(numImages) :: &
    tmag, &
    smag, &
    totmag, &
    stress_mag

  CHARACTER(len=3) :: &
    stepStatus

  REAL(kind=DP), DIMENSION(numImages,numAtoms,3) :: &
    fractPosition
 
                            !>> INITIALIZATION <<!
  tmag=0._DP
  smag=0._DP
  totmag=0._DP
  stress_mag=0._DP
  IF(badStep) THEN
    stepStatus = "bad"
  ELSE
    stepStatus = "   "
  END IF
  fractPosition = CartToFract(cell, cartPosition)

                             !>> FUNCTION BODY <<!
  DO ni=1,numImages
    DO na=1,numAtoms
      DO m=1,3
        tmag(ni)=tmag(ni)+(trueForce(ni,na,m)**2)
        smag(ni)=smag(ni)+(springForce(ni,na,m)**2)
        totmag(ni)=totmag(ni)+(totalForce(ni,na,m)**2)
      END DO
    END DO
    DO m=1,3
      DO n=1,3
        stress_mag(ni)=stress_mag(ni)+(stressTensor(ni,m,n)**2)
      END DO
    END DO
  END DO
  stress_mag=SQRT(stress_mag)
  tmag=SQRT(tmag)
  smag=SQRT(smag)
  totmag=SQRT(totmag)

  WRITE(transOutputUnit,'(A2, 3I5, 4ES18.9, 1X, E10.5, 1X, A3, A2)') "* ", numcalcs, numiterations, linecounts, MAXVAL(energy)-energy(1), MAXVAL(energy)-energy(numImages), SQRT(MAXVAL(SUM(totalForce**2,3))), SQRT(SUM(totalForce**2)/numImages/numAtoms), timeStep, stepStatus, " *"
  WRITE(*, '(A2, 3I5, 4ES18.9, 1X, E10.5, 1X, A3)') "* ", numcalcs, numiterations, linecounts, MAXVAL(energy)-energy(1), MAXVAL(energy)-energy(numImages), SQRT(MAXVAL(SUM(totalForce**2,3))), SQRT(SUM(totalForce**2)/numImages/numAtoms), timeStep, stepStatus

  WRITE(transOutputUnit,*) ''
  WRITE(transOutputUnit,'(A)') "   | IMAGE      ENERGY               TRUE_F              SPRING_F           TOTAL_F           INTERVAL        |"
  DO ni=1,numImages
    dmag = 0._DP
    IF(ni /= numImages) THEN
      DO na=1,numAtoms
        DO m=1,3
          dmag=dmag+((cartPosition(ni,na,m)-cartPosition(ni+1,na,m))**2)
        END DO
      END DO
    END IF
    dmag = SQRT(dmag)
    WRITE(transOutputUnit,'(A,I4,1X,5ES20.12,A)') "   |",ni,energy(ni),tmag(ni),smag(ni),totmag(ni),dmag," |"
  END DO
  WRITE(transOutputUnit,*) ''

  ! Flush the write buffer
  flush(transOutputUnit)

END SUBROUTINE OutputWrite



SUBROUTINE OutputWriteFooter()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Writes to the output file for the transition state finder
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/18/2005 Commented and cleaned up the output file (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                           !>> EXTERNAL VARIABLES <<!

                           !>> INTERNAL VARIABLES <<!
                            !>> INITIALIZATION <<!
                             !>> FUNCTION BODY <<!

  WRITE(transOutputUnit,'(A)') "-------------------------------------------------------------------"

  ! Flush the write buffer
  flush(transOutputUnit)

END SUBROUTINE OutputWriteFooter


SUBROUTINE WriteRestartInformation(cartPosition, cell, atomicPSP, &
                                   atomicSymbol)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine writes out the ion coordinates in all the images so we
!   can keep track of them as the run goes ... also, for restarting if the
!   job crashes.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: &
    cell

  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3), INTENT(IN) :: & 
    cartPosition

  CHARACTER(LEN=3), DIMENSION(numAtoms), INTENT(IN) :: &
    atomicSymbol

  CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: &
    atomicPSP

                            !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    ni, na

  REAL(KIND=DP), DIMENSION(numImages,numAtoms,3) :: & 
    fractPosition


                            !>> INITIALIZATION <<!
  fractPosition = CartToFract(cell, cartPosition)

                             !>> FUNCTION BODY <<!


  OPEN(UNIT=1,FILE="restart.out")

  WRITE(1,'(A4)') "CELL"
  WRITE(1,'(3F20.15)') cell(1,1),cell(2,1),cell(3,1)
  WRITE(1,'(3F20.15)') cell(1,2),cell(2,2),cell(3,2)
  WRITE(1,'(3F20.15)') cell(1,3),cell(2,3),cell(3,3)
  WRITE(1,'(A9,I10)') "numAtoms=",numAtoms
  WRITE(1,'(A12,I10)') "numatomtype=",numatomtypes
  WRITE(1,'(A9,I10)') "numimage=",numImages
  WRITE(1,'(A)') atomicPSP
  WRITE(1,'(A10)') "FRACTIONAL"

  DO ni=1,numImages
    WRITE(1,'(A5,I5)') "IMAGE",ni
    DO na=1,numAtoms
      WRITE(1,'(A2,3F20.15)') atomicSymbol(na), fractPosition(ni,na,1), &
                              fractPosition(ni,na,2), fractPosition(ni,na,3)
    END DO
  END DO

  CLOSE(1)

END SUBROUTINE WriteRestartInformation


END PROGRAM CINEB

MODULE RhoOptSCG
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoOptSCG 
!     |_SUBROUTINE SqrtGradientMinimization
!       |_FUNCTION Constrain
!       |_SUBROUTINE BrentLineMin
!       |_SUBROUTINE MinBracket
!       |_FUNCTION Step
!------------------------------------------------------------------------------

  USE Constants, ONLY: tiny
  USE RhoOptimizers
  USE KEDF_WGCkernel, ONLY: firstOrderWGC
  USE OUTPUT, ONLY: outputUnit
  USE Output, ONLY: WrtOut
  USE CellInfo, ONLY: cell
  USE MathFunctions, ONLY: volume
  USE Sys, ONLY: rhoR
  USE Sys, ONLY: interior
  USE Sys, ONLY: gridPack
  USE Timer, ONLY: TimerStart
  USE Timer, ONLY: TimerStop
  USE Timer, ONLY: stopwatch 

  USE Report, ONLY: MinimizerReportFooter ! Subroutine that prints out the results.
  USE Report, ONLY: MinimizerReportHeader ! Prints out the header for the minimization report
  USE Report, ONLY: MinimizerReportSteps  ! Subroutine that prints individual steps

  IMPLICIT NONE

  REAL(KIND=DP) :: sumRho = 0.0_DP
  REAL(KIND=DP) :: tiny2 = tiny*0.999_DP ! A number slightly below tiny

CONTAINS


SUBROUTINE SqrtGradientMinimization(grid, energy, cgmin, maxRhoIter)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine minimizes the energy with respect to the SQRT of the 
!   electron density  using the conjugate gradient / steepest descent 
!   methodology
!------------------------------------------------------------------------------
  USE OUTPUT, ONLY : DUMP
  USE MathFunctions, ONLY: & 
    MinMaxVal                  

  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!

  TYPE(gridPack), INTENT(INOUT) :: grid
  ! The pack of grids that we have in our system, contains the density, etc.
  !
  REAL(kind=DP), DIMENSION(:), INTENT(OUT) :: energy
  ! where the final energies are stored   
  !
  LOGICAL, INTENT(IN) :: cgmin          
  ! .true. for conjugate gradient
  ! .false. for steepest descent
  !
  INTEGER :: maxRhoIter
  
                    !>> INTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(LBOUND(rhoR,1):UBOUND(rhoR,1), &
                           LBOUND(rhoR,2):UBOUND(rhoR,2), &
                           LBOUND(rhoR,3):UBOUND(rhoR,3), &
                           SIZE(rhoR,4)) :: &
    g, &          ! The current gradient after applying constraints
    dirNext, &    ! The next direction that we would like to search in,
                  ! can be either the conjugate or steepest descent direction.
    pot           ! The total potential for the current step density.
    
  REAL(kind=DP) :: &
#ifdef __USE_PARALLEL
    mpiLocalSum, &! Sum on local processor before MPI call
#endif
    gamma, &      ! The parameter by which we scale the previous direction in
                  ! the conjugate-gradient scheme to find the new direction.
                  ! See 10.6.5 in Ref 1.  This parameter is beta in Ref 2.
    lastStep, &   ! A storage variable for the step size in the last direction
                  ! returned by our line minimizer.  Needed for  
                  ! preconditioning.
    potentialNorm, &     ! The SQRT of the sum of the potential
    lastPotentialNorm, & ! Last value of the SQRT of sum of potential
    lambda, &
    minDen, maxDen, &
    minPot, maxPot

#ifdef __USE_PARALLEL
  REAL(kind=DP), DIMENSION(2) :: &
    mpiLocalSum2, &! 2-dimensional array version of mpiLocalSum
    mpiGlobalSum2  ! 2-d array version of mpiGlobalSum
#endif

  INTEGER :: &
#ifdef __USE_PARALLEL
    mpiLocalIntSum, &   ! Local sum before MPI call for integer values
    mpiGlobalIntSum, &  ! Sum after MPI summing call
#endif
    iter, &       ! Minimization step counter
    success       ! whether our line minimization succeeded
                  !  =1 if it did
                  !  =2 if it didn't, went to max timestep
                  !  =3 if it exited because energy is increasing at small step

  LOGICAL :: &
    restart, &    ! whether the present step restarted from the steepest
                  ! descent vector if we were doing conjugate gradient.
    quitNextFail  ! whether the last iteration was a restart or not

  REAL(kind=DP), DIMENSION(9) :: &
    tenEnergiesAgo

  TYPE(stopwatch):: &
    watch        ! Timer

  character (len=500) :: &
    message

                      !>> INITIALIZATION <<!

  watch=TimerStart()
  rhoOutcome = -1
  success = 1
  quitNextFail = .FALSE.  

  lastStep = 0._DP

  ! This just sets the sum of the density so that functions like Step can 
  ! use it.
#ifdef __USE_PARALLEL
  mpiLocalSum = SUM(rhoR, interior)
  CALL MPI_ALLREDUCE(mpiLocalSum, sumRho, 1, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) THEN
    CALL WrtOut(6,"SqrtGradientMinimization: PROBLEMS WITH MPI_ALLREDUCE, code STOP!!")
    STOP
  ENDIF
#else
  sumRho = SUM(rhoR, interior)
#endif

  ! First compute the starting energy, and count it as an evaluation
  CALL MinimizerReportHeader

  ! Calculate the initial energy and potential
  CALL CalculatePotentialPlus(rhoR, .TRUE., pot, energy)

  ! We're going to do something messy here to conserve memory.  We're going to
  ! use the storage of rhoR to store the square root of rho  temporarily
  ! for the purposes of this calculation.
  rhoR= SQRT(rhoR)

  ! Make sure variable storing last energies are greater than the energy
  tenEnergiesAgo = energy*2._DP
  pot = Constrain(pot, rhoR, interior, lambda,size(pot,1),size(pot,2),size(pot,3),size(pot,4))
  numEnergyTotal = 1
  numPotentialTotal = 1

                    !>> FUNCTION BODY <<!

  ! Loop until we find a minimum!
  iter = 0

  WRITE(message, '(A)')'Iter    totEnergy(Ha)     potNorm  suc       min/max(den)         min/max(pot)'
  CALL WrtOut(6, message)

  DO

    restart = .FALSE.

    ! Here we exit if we took more than maxRhoIter iterations.  
    IF (iter > maxRhoIter) THEN 
      rhoOutcome = 1 
      EXIT
    END IF

    lastPotentialNorm = potentialNorm
    ! Calculate (using MPI) SUM(pot**2,interior)/COUNT(interior)
#ifdef __USE_PARALLEL
    mpiLocalSum = SUM(pot**2, interior)
    CALL MPI_ALLREDUCE(mpiLocalSum, potentialNorm, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr /= MPI_SUCCESS) THEN 
      CALL WrtOut(6,"m PROBLEMS WITH MPI_ALLREDUCE IN SqrtGradientMinimization***")
      STOP
    ENDIF

    mpiLocalIntSum = COUNT(interior)
    CALL MPI_ALLREDUCE(mpiLocalIntSum, mpiGlobalIntSum, 1, MPI_INTEGER, &
                       MPI_SUM, MPI_COMM_WORLD, mpiErr)
    IF (mpiErr /= MPI_SUCCESS) THEN
      CALL WrtOut(6, "m PROBLEMS WITH MPI_ALLREDUCE IN SqrtGradientMinimization***")
      STOP
    ENDIF

    potentialNorm = SQRT(potentialNorm/REAL(mpiGlobalIntSum,kind=DP))
#else
    potentialNorm = SQRT(SUM(pot**2, interior) &
                           / REAL(COUNT(interior),kind=DP))
#endif

    ! Try to detect a WGC blowup ...
    IF(iter > 1 .AND. potentialNorm/lastPotentialNorm > 1.E3_DP) THEN
      WRITE(*,*) "WGC BLOWUP DETECTED ... COMPENSATING ... "
    END IF

    ! If the constrained potential is tiny we are close to a critical point, 
    ! so stop.
    IF (potentialNorm <  pot_tol) THEN
      rhoOutcome = 0
      iter = iter + 1
      CALL MinimizerReportSteps(iter, energy, potentialNorm, 1, 0, &
                                restart, 1, lastStep)
      EXIT
    ELSE IF (potentialNorm < 100._DP*pot_tol .AND. firstOrderWGC==1) THEN
      firstOrderWGC = -1
      IF(calcReal) THEN
      ELSE ! periodic case
!        CALL CalculatePotentialSqrt(rhoR**2, pot)
        CALL CalculatePotentialPlus(rhoR**2, .TRUE., pot, energy)

      END IF

      numEnergyTotal = numEnergyTotal + 1
      numPotentialTotal = numPotentialTotal + 1
      CYCLE
    END IF

    ! Add a second convergence criteria based on the energy
    ! We will consider the system converged if the energy after 
    ! 10 iterations is within 10x the machine precision as the energy
    ! before 10 iterations.
    IF (MOD(iter,10) == 1) THEN
      IF( iter > 10 .AND. ABS(tenEnergiesAgo(1) - energy(1)) &
           < 10._DP*(NEAREST(energy(1),energy(1)*2._DP)-energy(1))) THEN
        rhoOutcome = 0
        iter = iter + 1
        CALL MinimizerReportSteps(iter, energy, potentialNorm, 1, 0, &
                                  restart, 1, lastStep)
        EXIT 
      END IF
      tenEnergiesAgo = energy
    END IF

    ! If this is the first loop, then use the steepest descent vector no matter
    ! what.  Otherwise, use either the steepest descent or conjugate gradient
    ! direction, depending on what we told it to do.
    IF(iter == 0) THEN
      CALL MinimizerReportSteps(iter, energy, potentialNorm, 1, 0, &
                                restart, 1, 0._DP)
      dirNext = - pot
    ELSE
      ! If we're not doing CGMIN, set the conjugate direction to be same as 
      ! the steepest descent direction.
      IF(.NOT. cgmin) THEN
        dirNext = -pot

        IF(success == 3) THEN
            CALL WrtOut(6,"FATAL ERROR: INCREASING ENERGY CLOSE TO THE INITIAL POINT.")
            rhoOutcome = 2
            EXIT
        END IF

      ! Here, find the conjugate gradient.
      ELSE 

        ! If the line minimizer failed ...
        IF(success == 3) THEN

          ! If we've done this one time already, tried steepest descent,
          ! and things still failed, give up
          IF(quitNextFail) THEN
            CALL WrtOut(6,"FATAL ERROR: INCREASING ENERGY CLOSE TO INITIAL POINT.")
            rhoOutcome = 2
            EXIT
          END IF

          ! Here, restart from steepest descent
          quitNextFail = .TRUE.
          restart = .TRUE.
          dirNext = - pot

        ELSE
          quitNextFail = .FALSE.

          ! Here is the polak-ribere method for determining the next direction.
          ! If gamma is less than 0, then set it to zero (steepest-descent 
          ! direction) to guarantee convergence of the polak-ribiere method.
          ! Use MPI to get gamma = SUM(pot*(pot-g))/SUM(g*g) globally
#ifdef __USE_PARALLEL
          mpiLocalSum2(1) = SUM(pot*(pot-g),interior)
          mpiLocalSum2(2) = SUM(g*g,interior)
          CALL MPI_ALLREDUCE(mpiLocalSum2, mpiGlobalSum2, 2, MPI_REAL8, &
                             MPI_SUM, MPI_COMM_WORLD, mpiErr)
          IF (mpiErr /= MPI_SUCCESS) STOP &
            "***PROBLEMS WITH MPI_ALLREDUCE IN SqrtGradientMinimization***"

          gamma = MAX(mpiGlobalSum2(1)/mpiGlobalSum2(2), 0._DP)
#else
          gamma = MAX(SUM(pot*(pot-g))/SUM(g*g), 0._DP)
#endif
          dirNext = - pot + gamma * dirNext

        END IF
      END IF
    END IF

    ! Now that we have our next direction in dirNext, go down the line and 
    ! minimize.  

    CALL BrentLineMin(&
      grid, dirNext, &    ! the density and the minimization direction.
      lambda, &
      energy,&            ! line minimizer puts energy here at min point
      lastStep, &         ! the timestep taken to the minimum point
      success &           ! .true. if the minimization was successful
    )

    ! update some counters.
    iter = iter + 1
    numEnergyTotal = numEnergyTotal + numEnergyLineMin + numEnergyBracket

    ! Update G
    g = pot

    ! Calculate the new potential
    IF(calcReal) THEN
    ELSE
!      CALL CalculatePotentialSqrt(rhoR**2, pot)
      CALL CalculatePotentialPlus(rhoR**2, .TRUE., pot, energy)

    END IF
    numPotentialTotal = numPotentialTotal + 1

    pot = Constrain(pot, rhoR, interior, lambda,size(pot,1),size(pot,2),size(pot,3),size(pot,4)) 

    CALL MinimizerReportSteps(iter, energy, potentialNorm, &
      numEnergyLineMin + numEnergyBracket, numEnergyBracket, &
      restart, success, lastStep)

    maxDen = MinMaxVal(rhoR(:,:,:,1)**2, 'max')
    minDen = MinMaxVal(rhoR(:,:,:,1)**2, 'min')
    maxPot = MinMaxVal(pot(:,:,:,1)**2, 'max')
    minPot = MinMaxVal(pot(:,:,:,1)**2, 'min')

    WRITE(message, & 
      '(I3, Es20.10, Es11.2, I3, Es11.2,A,Es9.2, Es11.2,A,Es9.2)')&
      iter,energy(1),potentialNorm,success,minDen,'/',maxDen, minPot,'/',maxPot
    CALL WrtOut(6,message)

  END DO

  CALL MinimizerReportFooter(iter, rhoOutcome, energy, numEnergyTotal, &
                             numPotentialTotal, TimerStop(watch))


  ! Make sure to end up with the actual rho, not the square root of rho
  rhoR = rhoR**2

! mohan add 2012-12-20
  WRITE(message,'(A,Es19.12,A,Es11.4,A,F10.1)') & 
    & ' Final total energy: ', energy(1), &
    & ' Ha, volume=', volume(cell%cellReal), ' bohr^3' 

  CALL WrtOut(6,message)



  CONTAINS

  FUNCTION Constrain(step, rho, mask, mu, dimX, dimY, dimZ, nspin)
  !----------------------------------------------------------------------------
  ! DESCRIPTION:
  ! This function takes a potential and essentially makes the average of
  ! its components zero by subtracting an offset.  It's essentially a lagrange
  ! multiplier that makes the potential norm-conserving.
  !----------------------------------------------------------------------------
    IMPLICIT NONE
                         !>> EXTERNAL VARIABLES <<!

    INTEGER, INTENT(IN) :: dimX, dimY, dimZ, nspin

    REAL(kind=DP), DIMENSION(dimX,dimY,dimZ,nspin), INTENT(IN) :: &
      step, &                   ! The step we want to make norm-conserving
      rho

    LOGICAL, DIMENSION(dimX,dimY,dimZ,nspin), INTENT(IN) :: &
      mask

    REAL(kind=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: & 
      Constrain                ! Our answer

    REAL(kind=DP), INTENT(OUT) :: &
      mu                       ! The lagrange multiplier

                     !>> INTERNAL VARIABLES <<!
                      !>> INITIALIZATION <<!
                      !>> FUNCTION BODY <<! 


  ! Use MPI calls to calculate mu = SUM(step,mask)/(2*SUM(rho,mask)) globally
#ifdef __USE_PARALLEL
    mpiLocalSum2(1) = SUM(step,mask)
    mpiLocalSum2(2) = SUM(rho,mask)
    CALL MPI_ALLREDUCE(mpiLocalSum2, mpiGlobalSum2, 2, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr /= MPI_SUCCESS) STOP &
      "***PROBLEMS WITH MPI_ALLREDUCE IN Constrain***"

    mu = mpiGlobalSum2(1)/mpiGlobalSum2(2)/2._DP
#else
    mu = SUM(step,mask)/SUM(rho,mask)/2._DP
#endif

    WHERE(mask)
      Constrain = step - 2._DP * rho * mu
    ELSEWHERE
      Constrain = 0.0_DP
    END WHERE

  END FUNCTION Constrain


  SUBROUTINE BrentLineMin(grid, unitDirection, lambda, finalEnergy, timeB, &
                          success)
  !----------------------------------------------------------------------------
  ! DESCRIPTION:
  !   This is a line-minimization routine that runs off in the direction of
  !   unitDirection (assumes that unitDirection is norm-conserving) and CHANGES
  !   the stored density in GridObjects to the density at which the energy is 
  !   minimized in that direction.  It also CHANGES the stored energy to the 
  !   energy corresponded to by the current density.
  !   It also returns "finalDir", which is the potential (norm-conserving) 
  !   associated with the final electron density.  It also returns timeB, or
  !   the timestep required to get to the minimum.  It also returns the number 
  !   of energy evaluations it took for both the line minimization part and the 
  !   bracketing part, as well as if the line minimization was successful.
  !
  !   "runType" is an argument that dictates some of the parameters in this line
  !   minimization.  0 is for just the golden section search without parabolic
  !   interpolation, 1 is for the standard Brent line minimization with the
  !   parabolic interpolation.
  !
  !   The exit condition is one of the trickiest things in this whole thing.
  !   Right now, it exits if the distance between the min point and one of the
  !   bracketing intervals is less than a certain percent of the min point 
  !   away from the min point.
  !
  ! GLOBAL/MODULE VARIABLES CHANGED:
  !   GridObjects.f90:
  !     rhoR
  !   Output.f90
  !     energy
  !
  ! CONDITIONS AND ASSUMPTIONS:
  !   BE FOREWARNED:  Right now, it exits when the next parabolic fit that is
  !   found is close to a point we've already calculated.  Truth be told, this
  !   is a LOUSY idea.  It only works for very smooth curves.
  !
  ! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
  !
  ! REFERENCES:
  ! 1.  Press, William H., et. al.  Numerical Recipies in Fortran: The Art
  !     of Scientific Computing.  Cambridge University Press, Cambridge, 1992.
  !     pg. 399 and on.
  !----------------------------------------------------------------------------
  ! REVISION LOG:
  !   01/25/2004  Subroutine created.  (GSH)
  !
  !----------------------------------------------------------------------------
    IMPLICIT NONE

                             !>> EXTERNAL VARIABLES <<!

    TYPE(gridPack), INTENT(INOUT) :: &
      grid

    REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: &
      unitDirection             ! The direction that we want to move in.

    REAL(kind=DP), DIMENSION(:), INTENT(INOUT) :: &
      finalEnergy               ! Coming in, this is the energy at the initial
                                ! starting density.  Going out, this is the 
                                ! final energy at the minimized point, x

    REAL(kind=DP), INTENT(IN) :: &
      lambda

    REAL(kind=DP), INTENT(INOUT) :: &
      timeB                     ! The timestep to the minimum.

    INTEGER, INTENT(OUT) :: &
      success                   ! whether we succeeded in finding a minimum
                                ! to our specifications.

  
                            !>> INTERNAL VARIABLES <<!

    REAL(kind=DP), DIMENSION(SIZE(finalEnergy)):: &
      energyA, energyC, &       ! The energies at the above times (timeA, etc..)
      energyB, energyB1, &
      energyB2, &
      energyGuess               ! Storage variable for the energy corresponding
                                ! to the guessed time.  

    REAL(kind=DP), DIMENSION(LBOUND(rhoR,1):UBOUND(rhoR,1), & 
                             LBOUND(rhoR,2):UBOUND(rhoR,2), &
                             LBOUND(rhoR,3):UBOUND(rhoR,3), &
                             SIZE(rhoR,4)) :: &
      tmpPot

    INTEGER :: &
      iter                      ! A counter to keep track of how many points
                                ! we've sampled

    REAL(kind=DP) :: &
      timeA, timeC, &   ! Bracketing interval of the minimum we're looking at
      timeB1, timeB2, & ! The 2nd min point so far and the 3rd min point so far,
                        ! repectively.
      mid, &            ! the midpoint between timeA and timeB
      timeGuess, &      ! our trial time to evaluate
      guess, &           ! the current proposed step
      parabolicMin, &   ! the step that our parabolic mimization would lead to
      temp1, temp2, &   ! temporary variables for parabolic calculation
      stepB4, &         ! the previous step
      stepB4Last, &     ! the step before the previous step
      ratio             ! our ratio for successive steps

    TYPE(gridPack) :: &
      tempGrid


                            !>> INITIALIZATION <<!

    ! Get 1 - the golden ratio
    ratio = (3._DP - SQRT(5._DP))/2._DP 
    numEnergyLineMin = 0

                            !>> FUNCTION BODY <<!

    ! First bracket the minima.  Find timeA, timeB, timeC and their energies 
    ! such that timeA < timeB < timeC, energyB < energyA, energyB < energyC.
    ! return SUCCESS = 1 if it was able to do so.
    CALL MinBracket(grid, unitDirection, lambda, finalEnergy, timeA, timeB, &
                    timeC, energyA, energyB, energyC, success)

    ! Initialize stepB4 and stepB4Last so that the first step will try a 
    ! parabolic fit.
    stepB4 = (timeC - timeA)*2.1_DP   
    stepB4Last = stepB4  

    ! If we successfully bracketed the minima, then ...
    IF(success == 1) THEN

      ! Initialize B1 and B2 for the first parabolic fit.  Assign B1 to the
      ! lower of A or C, and B2 to the other one.
      IF(energyA(1) < energyC(1)) THEN
        timeB1 = timeA
        energyB1 = energyA

        timeB2 = timeC
        energyB2 = energyC
      ELSE 
        timeB1 = timeC
        energyB1 = energyC

        timeB2 = timeA
        energyB2 = energyA
      END IF

      ! Start going down the line!
      iter = 1
      DO

        ! If the maximum number of iterations is exceeded, stop.
        IF(iter >= lineMaxIter) EXIT

        ! This is the normal exit condition.  If the bracketed interval is 
        ! a certain percentage of timeB, then we're done.
        IF((timeC - timeA)/timeB < lineTolPercent) EXIT

        ! Set the midpoint 
        mid = 0.5_DP * (timeA + timeC)

        ! Here is the golden ratio. We will default to this.
        IF (timeB > mid) THEN
          guess = ratio * (timeA - timeB)
        ELSE
          guess = ratio * (timeC - timeB)         
        END IF

        ! Now do the parabolic fit 
        IF(lineRunType == 1) THEN

          temp1 = (timeB - timeB1) * (energyB(1) - energyB2(1))
          temp2 = (timeB - timeB2) * (energyB(1) - energyB1(1)) 
          parabolicMin = 0.5_DP * &
            ((temp2 * (timeB - timeB2)) - (temp1 * (timeB - timeB1))) / &
             (temp1 - temp2)

          ! Test to see if the fit is acceptable.  Two conditions:
          ! 1.  Must fall between the interval timeA and timeC, and
          ! 2.  Must be less than 1/2 the step before last.
          IF(ABS(parabolicMin) < ABS(.5_DP * stepB4Last) .AND. &
            timeB + parabolicMin > timeA + tiny2 .AND. &
            timeB + parabolicMin < timeC - tiny2) THEN
   
            ! This is to safeguard the guess from being too small.
            IF(ABS(parabolicMin) < tiny2) parabolicMin = &
              SIGN(tiny2, parabolicMin)
            guess = parabolicMin
            stepB4Last = stepB4

          ELSE
            ! After every golden ratio step, make sure we attempt a parabolic.
            stepB4Last = MAX(timeA - timeB, timeC - timeB)
          END IF

        ELSE
          ! This is so that after every golden ratio step, we'll do a parabolic.
          stepB4Last = MAX(timeA - timeB, timeC - timeB)

        END IF

        ! Remember the last step
        stepB4 = guess

        ! This is the time that we're going to guess.
        timeGuess = timeB + guess

        ! Our function evaluation, increment the # of times we've evaluated
        ! the energy.
        IF(calcReal) THEN
        ELSE 
          CALL CalculatePotentialPlus(Step(rhoR, timeGuess, unitDirection, interior), & 
                .TRUE., tmpPot, energyGuess) 

        END IF
        numEnergyLineMin = numEnergyLineMin + 1

        ! Now we do housekeeping.

        ! First possibility is that the guessed energy is lower than at B.
        ! Then, Great!
        IF (energyGuess(1) < energyB(1)) THEN
          ! replace either a or b with x ...
          IF(timeGuess >= timeB) THEN
            timeA = timeB
            energyA = energyB
          ELSE 
            timeC = timeB
            energyC = energyB
          END IF
   
          timeB2 = timeB1
          energyB2 = energyB1 
    
          timeB1 = timeB
          energyB1 = energyB
 
          timeB = timeGuess
          energyB = energyGuess

        ! Otherwise, still shrink our borders
        ELSE
          IF(timeGuess >= timeB) THEN
            timeC = timeGuess
            energyC = energyGuess
          ELSE
            timeA = timeGuess
            energyA = energyGuess  
          END IF

          ! Check if its the second lowest energy found
          IF (energyGuess(1) < energyB1(1)) THEN
            timeB2 = timeB1
            energyB2 = energyB1
   
            timeB1 = timeGuess
            energyB1 = energyGuess

          ! Check if its the third lowest energy found
          ELSE IF(energyGuess(1) < energyB2(1)) THEN
            timeB2 = timeGuess
            energyB2 = energyGuess
          END IF

        END IF

        iter = iter + 1
      END DO
   
    ENDIF

    ! Finally, change the stored density and energy for this density.
    rhoR = SQRT(Step(rhoR, timeB, unitDirection, interior))
    finalEnergy = energyB

  END SUBROUTINE BrentLineMin


  SUBROUTINE MinBracket(grid, direction, lambda, currentE, timeA, timeB, &
                        timeC, energyA, energyB, energyC, success)
  !----------------------------------------------------------------------------
  ! DESCRIPTION:
  !   This subroutine simply seeks to find three points that point to a 
  !   minima.  The criteria for this is that it seeks to find points a, b, and c
  !   such that a < b < c and fa > fb and fc > fb.
  !
  !   Upon failure, it returns the minimum point found in timeB with its 
  !   corresponding energy in ENERGYB and success = .FALSE.
  !
  ! CONDITIONS AND ASSUMPTIONS:
  !   We assume that DIRECTION gives the proper direction to go downhill in 
  !   energy.
  ! 
  ! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
  !   This is a pretty heuristic algorithm that may be improved.
  !   This is also not foolproof.
  !----------------------------------------------------------------------------
    IMPLICIT NONE
                           !>> EXTERNAL VARIABLES <<!

    TYPE(gridPack), INTENT(IN) :: &
      grid

    REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: &
      direction         ! The gradient direction that we're headed in. 
                        ! Normalized.

    REAL(kind=DP), INTENT(OUT) :: &
      timeA, &          ! The three timesteps that brackets the minima.
      timeB, &          ! If this subroutine failed, then timeB contains
      timeC             ! the timestep that gives the minimum energy found.
  
    REAL(kind=DP), DIMENSION(:), INTENT(IN) :: & 
      currentE          ! The energy of the current density without the timestep

    REAL(kind=DP), DIMENSION(:), INTENT(OUT) :: &
      energyA, &        ! Contains the energy at points B and C upon success.
      energyB, &        ! Upon failure, the lowest energy found is returned
      energyC           ! in energyB.

    REAL(kind=DP), INTENT(IN) :: &
      lambda

    INTEGER, INTENT(OUT) :: &
      success           ! Whether we succeeded in our attempt to bracket the 
                        ! minima.


                               !>> INTERNAL VARIABLES <<!

    REAL(kind=DP), PARAMETER :: & 
      initFrac = 0.1_DP, & ! Right now this is set to be optimized for a first
                             ! step that is very large (to the max)
      optRatio = 10._DP      ! The optimum ratio found so far for determining
                             ! successive steps with high accuracy

    REAL(kind=DP), DIMENSION(SIZE(energyA)) :: &
      energyGuess            ! The energy at our temporary point

    REAL(kind=DP), DIMENSION(LBOUND(rhoR,1):UBOUND(rhoR,1), & 
                             LBOUND(rhoR,2):UBOUND(rhoR,2), &
                             LBOUND(rhoR,3):UBOUND(rhoR,3), &
                             SIZE(rhoR,4)) :: &
      tmpPot
  
    REAL(kind=DP) :: &
      gold, &             ! The golden ratio
      timeGuess, &        ! The timestep that we're going to next
      parabolicMin, &     ! The result of our parabolic extrapolation for the 
                          ! next point.
      temp1, temp2, &     ! temporary variables used only for the parabolic 
                          ! extrapolation.
      init, &             ! The initial value of timeB that we'll start guessing
                          ! with.  We base it on the successful value of timeB
                          ! from the last step.
      gridMinimum

    REAL(kind=DP), SAVE :: &
      lastTime = -1._DP   ! The last value of timeB.  Used to determine init.

    TYPE(gridPack) :: &
      tempGrid

                           !>> INITIALIZATION <<!
                           
    ! silence compiler warning  (JMD)
    if(.false.) then
       write(*,*) "INFO: Lambda is ", lambda
    endif

    numEnergyBracket = 0

    ! This is the golden ratio.
    gold = (3._DP - SQRT(5._DP))/2._DP

    ! Here we determine what our timeB will start with based on the successful
    ! value of timeB from the previous iteration.
    IF(lastTime < 0._DP) THEN
      init = 1._DP
    ELSE
      init = MAX(lastTime, 5.E-3_DP)
    END IF

    ! Here are our initial three points.
    timeA = 0._DP
    energyA = currentE

    timeB = timeA + init
    timeC = timeA + init * optRatio

    CALL CalculatePotentialPlus(Step(rhoR, timeB, direction, interior), &
         .TRUE., tmpPot, energyB) 
    CALL CalculatePotentialPlus(Step(rhoR, timeC, direction, interior), &
         .TRUE., tmpPot, energyC) 


    numEnergyBracket = numEnergyBracket + 2

                            !>> FUNCTION BODY <<! 
    DO

      ! This is our exit condition.  Mission accomplished!
      IF(energyA(1) > energyB(1) .AND. energyC(1) > energyB(1)) THEN
        success = 1
        lastTime = timeB*2._DP
        EXIT

      ! If energyB is less than energyA, and we're NOT going downhill...
      ELSE IF(energyB(1) >= energyA(1)) THEN
#ifdef __USE_PARALLEL
        CALL MPI_ALLREDUCE(MINVAL(rhoR, interior), gridMinimum, 1, &
                           MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpiErr)
        IF (mpiErr /= MPI_SUCCESS) STOP &
          "***PROBLEMS WITH MPI_ALLREDUCE IN MinBracket***"
#else
        gridMinimum = MINVAL(rhoR, interior)
#endif

        IF ((timeB - timeA)/gridMinimum< 1.E-12) THEN
          success = 3
          timeB = 0._DP
          lastTime = -1._DP

          CALL CalculatePotentialPlus(Step(rhoR, timeB, direction, interior), &
                .TRUE., tmpPot, energyB)

          numEnergyBracket = numEnergyBracket + 1

          EXIT

        ELSE
          ! We just try a smaller value of timeB until we do get something in 
          ! the downhill direction.
          timeC = timeB
          energyC = energyB
          timeB = timeA + (timeC - timeA) / optRatio

          CALL CalculatePotentialPlus(Step(rhoR, timeB, direction, interior), &
                .TRUE., tmpPot, energyB) 

          numEnergyBracket = numEnergyBracket + 1

          CYCLE
        END IF

      ! The only other possiblity is if energyA > energyB > energyC. Then we
      ! just try to go out to further and further energy C.
      ELSE IF(energyC(1) <= energyB(1)) THEN
        ! Do a parabolic interpolation of the minima
        temp1 = (timeB - timeA) * (energyB(1) - energyC(1))
        temp2 = (timeB - timeC) * (energyB(1) - energyA(1))
        parabolicMin = timeB - &
          0.5_DP * (temp1 * (timeB - timeA) - temp2 * (timeB - timeC)) / &
                   (temp1 - temp2)

        ! Make sure that the parabolic interpolation doesn't give a point 
        ! that's too close.
        IF(parabolicMin > timeC + tiny2) THEN
          timeGuess = MAX(parabolicMin, timeC * optRatio)

        ! Otherwise, just go out to a larger version of timeC.
        ELSE
          timeGuess = timeC * optRatio
        END IF

        ! Do our function evaluation.
        CALL CalculatePotentialPlus(Step(rhoR, timeGuess, direction, interior), &
            .TRUE., tmpPot, energyGuess)

        numEnergyBracket = numEnergyBracket + 1

        ! ...And housekeep.
        timeA = timeB
        energyA = energyB
        timeB = timeC
        energyB = energyC

        timeC = timeGuess
        energyC = energyGuess

        CYCLE

      END IF

    END DO

  END SUBROUTINE MinBracket


  FUNCTION Step(sqrtRho, time, dir, mask)
  !----------------------------------------------------------------------------
  ! DESCRIPTION:
  !   This function takes a step in direction dir from density rho given the 
  !   timestep time.
  !----------------------------------------------------------------------------
    IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!
    REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: &
      sqrtRho, &
      dir
 
    REAL(kind=DP), INTENT(IN) :: &
      time

    LOGICAL, DIMENSION(:,:,:,:), INTENT(IN) :: &
      mask

    REAL(kind=DP), DIMENSION(SIZE(sqrtRho,1), SIZE(sqrtRho,2), &
                             SIZE(sqrtRho,3), SIZE(sqrtRho,4)) :: &
      Step

                         !>> INTERNAL VARIABLES <<!
    REAL(kind=DP) :: &
      normFactor

    INTEGER :: &
      i, j, k
  
                         !>> INITIALIZAITION <<!

    IF(SIZE(Step,4) /=1) STOP "***STEP CAN'T HANDLE MULTIPLE SPINS ****"

    Step = 0._DP         ! Removes built up numerical noise on boundaries
    normFactor = 0._DP   ! Normalizes electron number

                       !>> FUNCTION BODY <<!
    ! Doing it this way with the loop is faster (via an speed test w/ full
    ! optimization) than the easier-to-read way.
    DO k = 1, SIZE(Step,3)
      DO j = 1, SIZE(Step,2)
        DO i = 1, SIZE(Step,1)
          IF(mask(i,j,k,1)) THEN

            Step(i,j,k,1) = sqrtRho(i,j,k,1) + time * dir(i,j,k,1)

            ! One thing you can do is to divide rho by 10 if after the step is
            ! taken we have a negative value. This is to discourage "bouncing"
            ! of density that has already been set to 0. But this seems to
            ! sometimes accelerate deterioration of the WGC. Feel free to try it
            ! IF(Step(i,j,k,1) < 0._DP) Step(i,j,k,1) = Step(i,j,k,1)/10._DP

            Step(i,j,k,1) = Step(i,j,k,1)**2
            normFactor = normFactor + Step(i,j,k,1)

          ELSE
             Step(i,j,k,1) = sqrtRho(i,j,k,1)**2

          END IF
        END DO
      END DO
    END DO

#ifdef __USE_PARALLEL
    mpiLocalSum = normFactor
    CALL MPI_ALLREDUCE(mpiLocalSum, normFactor, 1, MPI_REAL8, MPI_SUM, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr/=MPI_SUCCESS) THEN
      CALL WrtOut(6,"***PROBLEMS WITH MPI_ALLREDUCE IN SQRTNEWTONMINIMIZATION***,STOP")
      STOP
    ENDIF
#endif

    normFactor = sumRho/normFactor

    WHERE(mask) Step = Step * normFactor

  END FUNCTION Step



END SUBROUTINE SqrtGradientMinimization

END MODULE RhoOptSCG

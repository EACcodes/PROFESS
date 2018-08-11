MODULE RhoLineSearch
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoLineSearch
!       |_SUBROUTINE getNextGradient
!       |_FUNCTION gradOftheta
!       |_SUBROUTINE checkGradient
!       |_SUBROUTINE lineSearch
!       |_Subroutine SpinThetaSearch
!       |_FUNCTION thetaCGDirection
!
! DESCRIPTION:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
! 1.  Press, W.H., et. al.  Numerical Recipies: The Art of Scientific
!     Computing.  Cambridge University Press, Cambridge, 1986.
! 2.  Gill, P.E., et.al.  Practical Optimization.  Academic Press, San Diego, 
!     1986 (12th printing in 2000).
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

                              !<< GLOBAL >>

  USE Constants, ONLY: DP
  USE MPI_Functions
  USE RhoDirCG, ONLY: InnerProd

  IMPLICIT NONE

  INTEGER :: maxLine = 200  
  ! the max allowed line search number

CONTAINS


SUBROUTINE getNextGradient(totEleNum, phi, dEdPhi, theta, gradth, dirNext)
!------------------------------------------------------------------------------
! DESCRIPTION:
! project & normalize the Newton direction, dirNext
! Also initialize theta such that tempPhi ~ phi + NewtonDirection.
! This is the theta that will be found by the line search if we are
! near the minimum.  The expression for theta is
!      theta = |dir_perp| / |phi + Newtondirection|
! where dir_perp is the part of the newton direction that is
! perpendicular to phi.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: numSpin

  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: totEleNum 
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: phi 
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: dEdPhi

  REAL(KIND=DP), DIMENSION(:), INTENT(OUT) :: theta 
  REAL(KIND=DP), DIMENSION(:), INTENT(OUT) :: gradth 
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(OUT) :: dirNext

  !! >> INTERNAL VARIABLES << !!

  REAL(KIND=DP) :: normDirPerp(2) 
  ! norm of new direction, [max spin dimension]
  !
  REAL(KIND=DP) :: normPhiDir(2)  
  ! norm of phi + direction
  !
  INTEGER :: isp = 0
  ! spin index

  CALL Title("getNextGradient")
  CALL StartClock('NextGradient')

  IF (numSpin == 1) THEN

    DO isp =1,numSpin

      ! (1) get the norm of (phi + direction_in)
      normPhiDir(isp) = SQRT(InnerProd(phi(:,:,:,isp)+dirNext(:,:,:,isp),phi(:,:,:,isp)+dirNext(:,:,:,isp)))

      ! (2) orthogonalize the direction_in to phi
      dirNext(:,:,:,isp) = dirNext(:,:,:,isp) - phi(:,:,:,isp)*InnerProd(phi(:,:,:,isp),dirNext(:,:,:,isp)) / totEleNum(isp)
      
      ! (3) calculate the norm of direction_new
      normDirPerp(isp) = SQRT(InnerProd(dirNext(:,:,:,isp),dirNext(:,:,:,isp)))

      ! (4) Here's the initial theta, the smaller the theta is,
      ! the more stable the iteration is.
      theta(isp) = MIN(theta(isp), normDirPerp(isp)/normPhiDir(isp))

      ! (5) Here's the part of the Newton direction perpendicular to phi, normalized
      dirNext(:,:,:,isp) = dirNext(:,:,:,isp) * SQRT(totEleNum(isp))/normDirPerp(isp)

      ! (6) get the new gradient
      gradth(isp) = gradOftheta(0.D0, dEdPhi(:,:,:,isp), phi(:,:,:,isp), dirNext(:,:,:,isp))
      ! not sure why 0.D0 is used here not theta, mohan note 2014-08-03,
      !WRITE(*,*) "theta=", theta(isp)
      !gradth(isp) = gradOftheta(theta(isp), dEdPhi(:,:,:,isp), phi(:,:,:,isp), dirNext(:,:,:,isp))

    ENDDO ! End isp

  ELSE ! spin =2

    DO isp = 1,numSpin

      normPhiDir(isp) = SQRT(InnerProd(phi(:,:,:,isp)+dirNext(:,:,:,isp), phi(:,:,:,isp)+dirNext(:,:,:,isp)))

      dirNext(:,:,:,isp) = dirNext(:,:,:,isp) - phi(:,:,:,isp)*InnerProd(phi(:,:,:,isp),dirNext(:,:,:,isp)) / totEleNum(isp)

      normDirPerp(isp) = SQRT(InnerProd(dirNext(:,:,:,isp),dirNext(:,:,:,isp)))

      ! Here's the initial theta
      theta(isp) = 0.d0

      ! Here's the part of the Newton direction perpendicular to phi, normalized
      dirNext(:,:,:,isp) = dirNext(:,:,:,isp) * SQRT(totEleNum(isp))/normDirPerp(isp)

      gradth(isp) = gradOftheta(theta(isp), dEdPhi(:,:,:,isp), phi(:,:,:,isp), dirNext(:,:,:,isp))

    ENDDO ! End isp

  ENDIF ! End numSpin
  CALL StopClock('NextGradient')

END SUBROUTINE getNextGradient


FUNCTION gradOftheta(theta,dEdPhi,phi,dirNext)
!----------------------------------------------------------------------------
! DESCRIPTION:
!
! This is a surboutine computes the dE/d(theta) (not dL/dtheta)
! at a given theta, because the number of electron will be conserved
! by tuning theta, so only dE/d(theta) is needed now, dL/dTheta is not
! needed, and will be equal to dE/dTheta in principle.
!
! Easy to get:
!    newPhi = sqrtRho*cos(theta) + dirNext*sin(theta)
!    dE/d(theta) = dE/d(sqrtRho)*d(sqrtRho)/d(rho)*d(Rho)/d(newPhi)*d(newPhi)/d(Theta)
!                = sqrtPotential * 1/(2*sqrt(rho))*2*newPhi*(-phi*sin(theta)+dirNext*cos(theta)
!                = sqrtPotential * (newPhi/sqrt(rho)) * (-phi*sin(theta)+dirNext*cos(theta)
!
!  since newPhi^2 = rho
!
!  so newPhi/sqrt(rho) = sign(newPhi), we get:
!    dE/d(theta) = sqrtPotential * sign(newPhi) * [-phi*sin(theta)+dirNext*cos(theta) ]
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! June-04-2008, Created by Chen Huang
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell

  IMPLICIT NONE

  REAL(KIND=DP), INTENT(IN) :: theta             
  ! the mixing parameter between old phi and new gradient
  !
  REAL(KIND=DP), INTENT(IN) :: dEdPhi(:,:,:)  
  ! dE/d(newPhi), newPhi is the new phi at new theta
  !
  REAL(KIND=DP), INTENT(IN) :: phi(:,:,:)        
  ! sqrt rho (newest one), no spin is considered now
  !
  REAL(KIND=DP), INTENT(IN) :: dirNext(:,:,:)    
  ! the direction that to go
  !
  REAL(KIND=DP):: gradOftheta

  ! >> function begins <<!

  gradOftheta = cell%dV*SUM( dEdPhi*(-phi*SIN(theta)+dirNext*COS(theta)) )

  ! reduce to a unit value.
  CALL ReduceRealLevel1(gradOftheta)

  RETURN

END FUNCTION gradOftheta


SUBROUTINE checkGradient(switchStp, switchStp2, gradth, grad_exit_flag, rhoOutcome, dirNext, dLdPhi)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: numSpin
  USE Output, ONLY: wrtOut
  USE Output, ONLY: QUIT

  IMPLICIT NONE

  LOGICAL, INTENT(OUT) :: switchStp, switchStp2
  !
  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: gradth
  ! gradient of theta 
  !
  LOGICAL, INTENT(INOUT) :: grad_exit_flag

  INTEGER, INTENT(OUT) :: rhoOutcome 
  ! 0: everything is fine
  ! 1: gradient is too large
  ! 2: other mistakes
  !
  REAL(KIND=DP) :: gradThreshould = 1.0e5
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: dirNext
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: dLdPhi
  !

  CALL Title("checkGradient")
  grad_exit_flag = .false.

  IF(numSpin==1) THEN
    IF(gradth(1)>gradThreshould) THEN !mohan add 10-28-12
      WRITE(message,*) "The gradient is too large, the value is ", gradth(1)
      CALL QUIT(message)
    ELSE IF(gradth(1)>=0._DP) THEN
      IF (switchStp .EQV. .TRUE.) THEN
        CALL WrtOut(6,'With the steep decent dir , I still get the gradth>0, I might have reached the minimum.')
        rhoOutcome = 1
        grad_exit_flag = .TRUE.
        RETURN
      ELSE
        CALL WrtOut(6,'gradth>0, replace Newton dir with the steepest decent dir')
        switchStp = .TRUE.
        dirNext(:,:,:,1) = -dLdPhi(:,:,:,1)
        WRITE(message,*) 'gradth due to newton dir => ', gradth(1)
        CALL WrtOut(6,message)
      ENDIF
    ELSE  ! check gradth
      grad_exit_flag = .TRUE.
      RETURN
    ENDIF

  ELSE ! Else numSpin
    IF(gradth(1)>gradThreshould) THEN !mohan add 10-28-12
      WRITE(message,'("The gradient 1 is too large, so PROFESS decide to stop")')
      CALL QUIT(message)
    ELSE IF(gradth(1)>=0._DP) THEN
      IF (switchStp .EQV. .TRUE. .AND. switchStp2 .EQV. .TRUE.) THEN
        CALL WrtOut(6,'With the steep decent dir I still get the gradth>0, I might have reached the minimum.')
        rhoOutcome = 1
        grad_exit_flag = .TRUE.
        RETURN
      ELSE
        CALL WrtOut(6,'gradth>0, I replaced Newton dir with the steepest decent dir')
        switchStp = .TRUE.
        dirNext(:,:,:,1) = -dLdPhi(:,:,:,1)
        WRITE(message,*) 'gradient due to newton dir => ', gradth(1)
        CALL WrtOut(6,message)
      ENDIF ! End switchStp
    ENDIF ! End gradth(1)

    IF(gradth(2)>gradThreshould) THEN !mohan add 10-28-12
      WRITE(message,'("The gradient 2 is too large, so PROFESS decide to stop")')
      CALL QUIT(message)
    ELSE IF(gradth(2)>=0._DP) THEN
      IF (switchStp .eqv. .TRUE. .AND. switchStp2 .eqv. .TRUE.) THEN
        CALL WrtOut(6,'with the steep decent dir , I still get the gradth>0, I might have reached the minimum.')
        rhoOutcome = 1
        grad_exit_flag = .true.
        RETURN
      ELSE
        CALL WrtOut(6,"gradth>0, I replaced Newton dir with the steepest decent dir")
        switchStp2 = .TRUE.
        dirNext(:,:,:,2) = -dLdPhi(:,:,:,2)
        WRITE(message,*) "gradient due to newton dir => ", gradth(2)
        CALL WrtOut(6,message)
      ENDIF ! End switchStp
    ENDIF ! End gradth(2)

    IF(gradth(1)< 0.0_DP .AND. gradth(2)<0.0_DP) THEN
      grad_exit_flag = .true.
      RETURN
    ENDIF ! End gradth(1)
  ENDIF ! End numSpin

  RETURN

END SUBROUTINE checkGradient


SUBROUTINE LineSearchTN(iter, theta, energy, gradth, ftol, gtol, xtol, success, rhoOutcome, &
 tempPhi, phi, dirNext, dEdPhi)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   do the line search using "Dcsrch"
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE Constants, ONLY: pi
  USE cellInfo, ONLY: numSpin
  USE RhoOptimizers, ONLY: numEnergyLineMin
  USE RhoOptimizers, ONLY: maxIter
  USE CalPotPlus, ONLY: CalculatePotentialPlus
  USE Output, ONLY: WrtOut
  USE Output, ONLY: QUIT 

  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!

  INTEGER :: iter
  REAL(KIND=DP), DIMENSION(:) :: theta
  REAL(KIND=DP), DIMENSION(:) :: energy
  REAL(KIND=DP), DIMENSION(:) :: gradth
  REAL(KIND=DP) :: ftol
  REAL(KIND=DP) :: gtol
  REAL(KIND=DP) :: xtol
  INTEGER :: success
  INTEGER :: rhoOutcome
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: tempPhi
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: phi
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: dirNext
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: dEdPhi


  !! >> LOCAL VARIABLES << !!

  INTEGER,DIMENSION(2,numSpin) :: isave
  ! variable used in Dcsrch
  !
  REAL(KIND=DP), DIMENSION(13,numSpin) :: dsave
  ! variable used in Dcsrch
  !
  CHARACTER(LEN=60) :: task1
  ! Status of line minimizer
  !
  REAL(KIND=DP) :: stpmin, stpmax

  CALL Title("LineSearch")

  stpmin = 0.0_DP ! lower bound for theta
  stpmax = pi     ! upper bound for theta

  IF(numSpin==1) THEN !spin-unpolarized case
    task1 = 'START'
    DO
      ! Do line search based on Wolfe conditions
      CALL Dcsrch(theta(1),energy(1),gradth(1),ftol,gtol,xtol,task1,stpmin,stpmax,isave,dsave)
      ! What's the message from line minimization ?
      IF(task1(1:4) .EQ. 'CONV')  THEN
        success = 1
        EXIT

      ELSE IF(task1(1:5) .EQ. 'ERROR') THEN
        success = -1
        WRITE(message,*) 'sqrtNewtonMinimization2: ERROR: Dcsrch()=>', task1
        CALL WrtOut(6,message)
        EXIT

      ELSE IF(task1(1:4) .EQ. 'WARN') THEN
        success = -2
        WRITE(message,*) 'sqrtNewtonMinimization2: WARNING: &
              & Dcsrch()=>', task1,' gradth:',gradth
        !!DEBUG
!        IF (task(10:13) .EQ. 'ROUN' .AND. i==2 ) THEN
!          CALL DebugLineSearch(tempRA, dirNext, dV, totEleNum, 1)
!        endif
        !!ENDDEBUG
        CALL WrtOut(6,message)
        EXIT  ! we ignore warning here, continue our optimization
      ELSE IF(task1(1:2) .EQ. 'FG') THEN ! need more line search, continue
        !Compute trial phi with new theta
        tempPhi = phi*COS(theta(1)) + dirNext*SIN(theta(1))
        CALL CalculatePotentialPlus(tempPhi**2,.TRUE.,dEdPhi,energy)
        dEdPhi = dEdPhi * sign(1._DP, tempPhi)

        gradth = gradOftheta(theta(1),dEdPhi(:,:,:,1),phi(:,:,:,1),dirNext(:,:,:,1))

        ! Count line minimization #
        numEnergyLineMin = numEnergyLineMin+1
        IF(numEnergyLineMin > maxLine) THEN
          CALL WrtOut(6,'sqrtNewtonMinimization2: WARNING:  &
                & exceed the max iteration number for line minimization.')
          CALL WrtOut(6,'skip rest line minimization steps, exit &
                & the line minimization loop.')
          EXIT
        ENDIF

      END IF

    END DO ! task1
  ELSE  !spin-polarized case
    !---------------------------------------------------------------------------------------------------
    !--instead of LINE SEARCH based on theta, do cg optimization based on vector theta -----------------
    !---------------------------------------------------------------------------------------------------
    CALL SpinThetaSearch(theta,energy,gradth,ftol,gtol,xtol,task1,isave,dsave,success, &
    tempPhi, phi, dirNext, dEdPhi)
    ! update the theta and gradth expression here.
    ! mohan 09/24/12
  ENDIF
! If we get an error message from the SpinThetaSearch  minimizer, we will not update
  ! rhoR, and we just restart our loop from the top
  IF (success == -1) THEN ! .OR. success(2) == -1 ) THEN
    IF (iter>maxIter) THEN
      CALL WrtOut(6, "sqrtNewtonMinimization2: Get ERROR message from line minimization,")
      CALL WrtOut(6, "and reached the maximum allowed iteration number,")
      CALL WrtOut(6, "Sever error in Rho Optimization! You should not trust the final results.")
      rhoOutcome = 1
      message="Fail rho optimizations"
      CALL QUIT(message)
    ENDIF
  END IF

  RETURN

END SUBROUTINE LineSearchTN


SUBROUTINE SpinThetaSearch(theta,energy,gradth,ftol,gtol,xtol,task1,isave,dsave,success,&
tempPhi, phi, dirNext, dEdPhi)
!-----------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012 Junchao Xia
! 2014-07-24 Mohan Chen
!------------------------------------------------------------------------------

  USE Constants, ONLY: pi
  USE CalPotPlus, ONLY: CalculatePotentialPlus
  USE CellInfo, ONLY: numSpin
  USE Output, ONLY: WrtOut

  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(INOUT) :: theta(numSpin)
  REAL(KIND=DP), INTENT(INOUT) :: gradth(numSpin)
  REAL(KIND=DP), INTENT(INOUT) :: energy(:)

  REAL(KIND=DP), INTENT(IN) :: ftol,gtol,xtol

  INTEGER, INTENT(INOUT) :: isave(2),success

  REAL(KIND=DP), INTENT(INOUT), DIMENSION(13):: dsave

  CHARACTER(LEN=60), INTENT(INOUT) :: task1

  REAL(KIND=DP), DIMENSION(:,:,:,:) :: tempPhi
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: Phi
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: dirNext
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: dEdPhi

                       !>> INTERNAL VARIABLES <<!
  INTEGER:: j,k,numThetaLineMin

  REAL(KIND=DP) :: alpha
  REAL(KIND=DP) :: dEdalpha
  REAL(KIND=DP) :: thetadir(numSpin)
  REAL(KIND=DP) :: absthetadir(numSpin)
  REAL(KIND=DP) :: thetadir_old(numSpin)
  REAL(KIND=DP) :: gradth_old(numSpin)
  REAL(KIND=DP) :: maxthetadir
  REAL(KIND=DP) :: alphatol, temtheta(2)

                          !>> FUNCTION BODY <<!

  alpha=0.0
  thetadir = 0.d0
  thetadir_old = 0.d0
  gradth_old = 0.d0
  alphatol = 1e-4
  j = 1
  numThetaLineMin = 1

  !write(message,*) 'Entering Thetasearch subroutine'
  !call wrtout(6,message)

  DO k=1,2
     tempPhi(:,:,:,k) = phi(:,:,:,k)*COS(theta(k)) + dirNext(:,:,:,k)*SIN(theta(k))
  ENDDO

  CALL CalculatePotentialPlus(tempPhi**2,.TRUE.,dEdPhi,energy)
  dEdPhi = dEdPhi * sign(1._DP, tempPhi)

  DO
    !debug
    !write(message,*) j,' th line search started'
    !call wrtout(6,message)
    !write(message,*) 'theta vector = (', theta(1),',',theta(2),')'
    !call wrtout(6,message)

    thetadir = thetaCGDirection(j,gradth_old,gradth,thetadir_old)
    !write(message,*) 'thetadir is (',thetadir(1),',',thetadir(2),')'
    !call wrtout(6,message)

    DO k=1,2
       gradth(k) = gradOftheta(theta(k),dEdPhi(:,:,:,k),phi(:,:,:,k), &
                            dirNext(:,:,:,k))
    enddo

    !write(message,*) 'gradth is (',gradth(1),',',gradth(2),')'
    !call wrtout(6,message)
    do
       dEdalpha = gradth(1)*thetadir(1)+gradth(2)*thetadir(2)
       if(dEdalpha>=0.d0) then
          thetadir = -gradth
       else
          exit
       endif
    enddo

    !write(message,*) 'alpha = ',  alpha, ' Energy = ', energy(1), ' Potential = ', dEdalpha
    !call wrtout(6,message)

    thetadir_old = thetadir
    absthetadir(1)=abs(thetadir(1))
    absthetadir(2)=abs(thetadir(2))
    maxthetadir = maxval(absthetadir)

    alpha = min (0.1, 0.1*pi/maxthetadir)
    !------------------------------------------------------------
    !--------------- LINE SEARCH based on alpha -----------------
    !------------------------------------------------------------
    task1 = 'START'

!      IF(abs(dEdalpha)<1e-3) then
!         CALL WrtOut(6,'Assume reached minimum')
!         success = -4

     numThetaLineMin = 0
!      else
       DO ! start newton minimization
       ! Do line search based on Wolfe conditions
        CALL Dcsrch(alpha,energy(1),dEdalpha,ftol,gtol,xtol,task1,0.d0,pi/maxthetadir, &
                   isave,dsave)
        ! What's the message from line minimization ?
        IF(task1(1:4) .EQ. 'CONV')  THEN
          success = 2
          EXIT
        ELSE IF(task1(1:5) .EQ. 'ERROR') THEN
          success = -1
          WRITE(message,*) 'sqrtNewtonMinimization2: &
          & ERROR: Dcsrch()=>', task1
          CALL WrtOut(6,message)
          EXIT
        ELSE IF(task1(1:4) .EQ. 'WARN') THEN
          success = -2
          WRITE(message,*) 'sqrtNewtonMinimization2: WARNING: &
          & Dcsrch()=>', task1,' dEdalpha:',dEdalpha
!!DEBUG
!         IF (task(10:13) .EQ. 'ROUN' .AND. i==2 ) THEN
!           CALL DebugLineSearch(tempRA, dirNext, dV, totEleNum, 1)
!         endif
!!ENDDEBUG
          CALL WrtOut(6,message)
          EXIT  ! we ignore warning here, continue our optimization

        ELSE IF(task1(1:2) .EQ. 'FG') THEN ! need more line search, continue

         !write(message,*) 'FG result from Line search'
         !call wrtout(6,message)

         !Compute trial phi with new alpha and theta
         temtheta = theta + alpha*thetadir

         !write(message,*) 'theta vector is updated = (', temtheta(1),',',temtheta(2),')'
         !call wrtout(6,message)

         do k=1,2
           tempPhi(:,:,:,k) = phi(:,:,:,k)*COS(temtheta(k)) + dirNext(:,:,:,k)*SIN(temtheta(k))
         enddo

         CALL CalculatePotentialPlus(tempPhi**2,.TRUE.,dEdPhi,energy)
         dEdPhi = dEdPhi * sign(1._DP, tempPhi)
         do k=1,2
           gradth(k) = gradOftheta(temtheta(k),dEdPhi(:,:,:,k),phi(:,:,:,k), &
                              dirNext(:,:,:,k))
         enddo
         dEdalpha = gradth(1)*thetadir(1)+gradth(2)*thetadir(2)
         !write(message,*) 'alpha = ',  alpha, ' Energy = ', energy(1), ' Potential = ', dEdalpha
         !call wrtout(6,message)
         !alpha = alpha + 0.1
         !if(alpha >5) exit
         absthetadir(1)=abs(thetadir(1))
         absthetadir(2)=abs(thetadir(2))
         maxthetadir = maxval(absthetadir)

         ! Count line minimization #
         numThetaLineMin = numThetaLineMin+1

!         IF(abs(alpha*dEdalpha)<1e-8) then
!           CALL WrtOut(6,'Assume reached minimum')
!           success = -4
!           exit
!         ENDIF

         IF(numThetaLineMin > 10) THEN
           success = -1
           CALL WrtOut(6,'sqrtNewtonMinimization2: WARNING:  &
            & exceed the max iteration number for line minimization.')
           CALL WrtOut(6,'skip rest line minimization steps, exit &
            & the line minimization loop.')
           EXIT
         ENDIF

        END IF

       END DO
!      ENDIF
    !-------------------------------------------------------------------
    !------- END OF LINE SEARCH ----------------------------------------
    !-------------------------------------------------------------------
     theta = theta + alpha * thetadir
     if(success == -1) exit
     if(success == -4) exit
     if(SQRT(gradth(1)*gradth(1) + gradth(2)*gradth(2))<alphatol) exit
     j = j + 1
     if(j>2) exit
  enddo

  RETURN

END SUBROUTINE SpinThetaSearch


FUNCTION thetaCGDirection(j,gradth_old,gradth,thetadir_old)
!----------------------------------------------------------------------------
! DESCRIPTION:
! Compute the conjugate gradient direction for vector theta
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014-07-24 Mohan Chen
!----------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(KIND=DP), DIMENSION(2) :: thetaCGDirection
  REAL(KIND=DP), DIMENSION(2), INTENT(INOUT) :: gradth_old
  REAL(KIND=DP), DIMENSION(2), INTENT(IN) :: gradth, thetadir_old
  INTEGER, INTENT(IN) :: j

  !! >> INTERNAL VARS << !!
  REAL(KIND=DP) :: cg_gamma

  IF (j==1) THEN
    cg_gamma = 0._DP
  ELSE
    ! Polak and Ribire formula
    cg_gamma = MAX(0._DP, &
    ThetaInnerProd(gradth,gradth-gradth_old) / &
    ThetaInnerProd(gradth_old,gradth_old))
  ENDIF
  gradth_old = gradth
  ThetaCGDirection = - gradth + cg_gamma * thetadir_old

  RETURN

END FUNCTION ThetaCGDirection


FUNCTION ThetaInnerProd(theta1,theta2)
!---------------------------------------------------------------
! Description:
! Compute the inner product of two vectors v1 and v2.
! v1 and v2 are three dimensional v1(x,y,z),v2(x,y,z)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014-07-24 Mohan Chen
!---------------------------------------------------------------
  IMPLICIT NONE

  REAL(KIND=DP), INTENT(in) :: theta1(:), theta2(:)

  REAL(KIND=DP) :: ThetaInnerProd

  !! >> FUNCTIONS << !!

  ThetaInnerProd = SUM(theta1*theta2)

  RETURN

END FUNCTION ThetaInnerProd


END MODULE RhoLineSearch

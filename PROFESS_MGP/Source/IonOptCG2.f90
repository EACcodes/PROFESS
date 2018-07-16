MODULE IonOptCG2
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE IonOptCG2
!     |_SUBROUTINE ConjugateGradient (CG2)
!
! DESCRIPTION:
!   This module contains optimization procedures used to optimize ion positions
!   given the forces.  
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

  USE CONSTANTS, ONLY : DP 

  USE Output, ONLY : WrtOut
  USE OutputFiles, ONLY : outputUnit

  USE MPI_Functions

  USE Timer, ONLY : TimerStart
  USE Timer, ONLY : TimerStop 

  USE CellInfo, ONLY: cell

  IMPLICIT NONE

CONTAINS


SUBROUTINE ConjugateGradient(RhoOptimizer, rho, energy, forces, frozenIon)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine uses steepest descent or conjugate gradient (Polak-Ribere)
!   to find the minimum energy ionic configuration. The code used a home-made
!   line search subroutine, based on strong wolfe condition. The code can exit
!   during the outer CG loop or over the inner line search, as long as the
!   maximum forces drops below the criterion. 
!  Our line search subroutine is more robust than dcserch.f(More and Thuente Alg.)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   7/20/2006  Created (Linda Hung)
!   Aug/19/2008,  (1) dcserch.f (More and Thuente Alg.) is used for line search
!                 (2) line search is revised by (Chen Huang)
!   09/11/2008  Added wrapping back into periodic box (L. Hung)
!   06/08/2008  A new home-made line search (LnSrch) is used, 
!               Add Wolfe condition check in the line search. (Chen Huang)
!
!------------------------------------------------------------------------------

  USE MathFunctions, ONLY: Inverse  

  USE LnSrch_MOD, ONLY: lnsrch 
  USE LnSrch_MOD, ONLY: work => points      !! only for monitor purpose, should no change

  USE Report, ONLY: GeometryMinimizerReportHeader
  USE Report, ONLY: GeometryMinimizerReportFooter
  USE Report, ONLY: GeometryMinimizerReportSteps

  USE IonOptimizers, ONLY: cg_type
  USE IonOptimizers, ONLY: maxIonStep
  USE IonOptimizers, ONLY: watch
  USE IonOptimizers, ONLY: watch2

  USE NearestDistance, ONLY: CheckNearestDistanceAtoms
  USE NearestDistance, ONLY: nnDist

  USE RefreshIons, ONLY: RefreshIonTerms
  USE CalForces, ONLY : CalculateForces          ! Calculate the forces.    

  IMPLICIT NONE

                          !>> EXTERNAL VARIABLES <<!

  EXTERNAL RhoOptimizer 
  ! A subroutine that is called to optimize the electron
  ! density relative to the ion positions
  !
  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho 
  ! Electron density, realspace
  !
  REAL(kind=DP), DIMENSION(:), INTENT(IN) :: energy
  !
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: forces           
  ! First index ion index, second index is x/y/z, 
  ! 3rd index: kinds of forces 1 for total force
  !
  LOGICAL, ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: frozenIon
  ! flag to indicate whether ions are frozen or not
  !

                          !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(cell%numIon,3) :: fracCoord
  ! fractional coordinates of atoms
  !
  REAL(KIND=DP), DIMENSION(cell%numIon,3) :: oldCoordCart
  ! Ion positions at beginning of line search
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: bracket_coord
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: delta_coord
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: dir_old
  ! old cg direction
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: dir
  ! current cg direction
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: oldgrad
  ! Vector in direction of old forces
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: newgrad
  ! Vector in direction of new forces
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: y_k
  ! newGrad - oldGrad
  !
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: tmp1
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: tmp2
  !
  REAL(KIND=DP) :: cg_eta
  !
  REAL(KIND=DP) :: maxForce
  ! Maximum component of force among all ions
  !
  REAL(KIND=DP) :: gam
  ! Gamma, for conjugate gradient
  !
  REAL(KIND=DP) :: alpha
  ! step in line search
  !
  REAL(KIND=DP) :: norm_dk
  REAL(KIND=DP) :: norm_gk
  REAL(KIND=DP) :: norm_yk
  REAL(KIND=DP) :: gradAlp     
  ! dE/dAlpha
  !
  INTEGER, PARAMETER :: maxit = 50 
  ! Maximum number of iterations for midpoint and endpoint search
  !
  INTEGER :: nLine
  ! line search counter
  !
  INTEGER :: maxLine
  ! max Line Search 
  !
  INTEGER :: ia
  !
  INTEGER :: step 
  !
  CHARACTER(len=70) :: task 
  ! Status of line minimizer 
  !
  CHARACTER(len=500) :: message
  INTEGER  :: checkNNdist
  REAL(KIND=DP) :: minDistance
  REAL(KIND=DP) :: cart2frac(3,3) 
  REAL(KIND=DP) :: frac2cart(3,3)
  LOGICAL :: restart_cg
  LOGICAL :: relax_done

  !----------------------------
  ! For line search subroutine 
  !----------------------------
  REAL(KIND=DP) :: disTol
  !
  REAL(KIND=DP) :: ftol = 1e-4  
  ! new energy should be smaller than the older one.
  !
  REAL(KIND=DP) :: gtol = 1e-1
  ! 5e-2 is for CG which needs a more exact line-search,
  !
  REAL(KIND=DP) :: xtol = 1.E-12_DP
  !
  REAL(KIND=DP) :: stpmin = 0._DP
  ! lower bound for theta
  !
  REAL(KIND=DP) :: stpmax = 0._dp     
  ! upper bound for theta
  !

                            !>> INITIALIZATION <<!

  CALL WrtOut(6,' Ion relaxation based on conjugate gradient, starting ...')
  SELECT CASE(cg_type)
  CASE(1)
    CALL WrtOut(6," Polak-Ribieare CG for ion relax")
  CASE(2)
    CALL WrtOut(6," Hager-Zhang CG for ion relax")
  CASE(3)
    CALL WrtOut(6," Dai-Yuan CG for ion relax")
  CASE DEFAULT
    CALL WrtOut(6," cg algrithm is not defined! code stop!")
    STOP
  END SELECT
  watch = TimerStart()
  watch2 = TimerStart()

  CALL GeometryMinimizerReportHeader
  oldGrad  = 0.d0
  newgrad  = 0.d0  ! Initialize to any nonzero number
  dir_old  = 0.d0
  dir      = 0.d0  ! Initialize to zero 
  maxLine  = 10
  
  ! Matrices for coordinate conversion
  frac2cart = cell%cellReal
  cart2frac = Inverse(cell%cellReal)

                      

                            !>> FUNCTION BODY <<!

  step = 0 ! counter for CG loops
  WRITE(outputUnit,'(/A)')   " --------------------------------------------------"
  WRITE(outputUnit,'(A,I5,A,I5)') "  IONIC OPTIMIZATION (CG2 method), step=", step, "/", maxIonStep
  WRITE(outputUnit,'(A)')    " --------------------------------------------------"

! Already set up ion positions during RefreshLatticeVectors
! Don't need to separately call RefreshIonTerms
! LINDA FIX THIS COMMENT
!  CALL RefreshIonTerms(cell, grids)  ! Call each time ion positions change
  CALL RhoOptimizer
  CALL CalculateForces(rho, forces)  ! Forces are in cartesion coord
  IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP

  !---------------------------------------------------------------------
  !---     CG LOOP     -------------------------------------------------
  !---------------------------------------------------------------------
  alpha  = 1.d0
  restart_cg = .false.
  relax_done = .false.
  CGloop: DO

    step = step+1

    IF( step > maxIonStep ) THEN
      WRITE(outputUnit,*) "Reach the allowed maximal ion steps = ", step
      EXIT
    ENDIF

    !-------- Check point (1) -------------
    CALL CheckExit(maxForce, relax_done)

    CALL StepInfo(step,maxForce,alpha,energy(1))

    IF (relax_done) EXIT

    ! Convert to cartesion coordinates
    DO ia=1,cell%numIon
      oldCoordCart(ia,:) = MATMUL(frac2cart,cell%ionTable(ia)%coord(:)) 
      fracCoord(ia,:) = cell%ionTable(ia)%coord(:)
    END DO
    CALL GeometryMinimizerReportSteps(& 
      step, TimerStop(watch), forces, maxForce, & 
      alpha, .FALSE., REAL(nLine,kind=DP),&
      .TRUE.,fracCoord )

    watch = TimerStart()

    WRITE(outputUnit,'(/A)')   " --------------------------------------------------"
    WRITE(outputUnit,'(A,I5,A,I5)') "  IONIC OPTIMIZATION (CG2 method), step=", step, "/", maxIonStep
    WRITE(outputUnit,'(A)')    " --------------------------------------------------"


    newGrad = forces(:,:,1)
    IF (step == 1 .or. restart_cg .eqv. .TRUE. .OR. MOD(step,30)==0) THEN
      gam = 0.d0
      restart_cg = .FALSE.
      IF (step/=1) CALL WRTOUT(6,"======= ION OPTIMIZATION: CG Restarted =========")
    ELSE
      SELECT CASE ( cg_type )
      CASE(1)
         ! More robust, Polak.Ribiere CG 
         gam = SUM((newGrad-oldGrad)*newgrad)/SUM(oldGrad**2) 
         gam = MAX(0._DP, gam)              
      CASE(2) 
         ! Reference: A Survey of nonlinear conjugate gradient methods, 
         !            William W. Hager and Hongchao Zhang  (2005)
         y_k = - newGrad + oldGrad
         norm_yk = SQRT(SUM(y_k**2))
         tmp1 = y_k - 2.d0*dir_old*norm_yk**2/SUM(dir_old*y_k)
         tmp2 = -newGrad / SUM(dir_old*y_k)
         gam = SUM( tmp1 * tmp2 )

         cg_eta = 0.01_DP
         norm_dk = SQRT(SUM(dir_old**2))
         norm_gk = SQRT(SUM(oldGrad**2))
         gam = MAX(gam, -1._DP/(norm_dk*MIN(cg_eta,norm_gk)))
      CASE(3) 
         y_k = - newGrad + oldGrad
         gam = SUM(newGrad**2) / SUM(dir_old * y_k)
         gam = MAX(0.d0, gam)
      END SELECT
    ENDIF
     
    ! Make new CG dir, and backup old CG dir and old gradient
    dir = newGrad + gam*dir_old
    dir_old = dir
    oldGrad = newGrad
    
    CALL gradOfAlpha(forces,dir,cell%numIon,gradAlp)

    stpMin = 0.d0
    IF ( nnDist > 0) THEN      
      checkNNdist = CheckNearestDistanceAtoms(minDistance)
      stpMax = ABS( (minDistance - nnDist) / MAXVAL(dir) )
    ELSE 
      ! we have disabled the CheckNearestDistanceAtoms() function
      ! we assume we at most move 10 bohr
      stpMax = ABS( 1.d0 / MAXVAL(dir)) 
    ENDIF

    alpha = MIN( alpha, stpMax )    ! initial alpha

    WRITE(message,'(a,Es12.4,a,Es12.4,a,Es12.4,a,Es12.4)') & 
      "For linesearch: stpMax = ", stpMax, " stpMin=", stpMin, & 
      " gradAlp=", gradAlp, " alpha_0=", alpha
    CALL WrtOut(6,message)
    
    !--------------------------------------------------
    !--- LINE SEARCH ----------------------------------
    !--------------------------------------------------
    task = "START"
    nLine = 0
    LINESEARCH: DO
       
       CALL LnSrch(alpha,energy(1),gradAlp,task,ftol,gtol,xtol,stpmin,stpmax)

       CALL WrtOut(outputUnit,"---------------------------------------------------------------------------")
       WRITE(message,'(3(a,ES12.4),a,I2)') " max_force=",maxForce," dE/d(alpha) =", gradAlp, "-> NEXT_Alp:=", alpha, &
         " LineSteps#:", nLine
       CALL WrtOut(outputUnit,message)
       WRITE(message,'(3(a,ES16.8))') " LEFT : alpha = ", work(1,1), " E=", work(1,2), " dEdAlp=", work(1,3)
       CALL WrtOut(outputUnit,message)
       WRITE(message,'(3(a,ES16.8))') " RIGHT: alpha = ", work(2,1), " E=", work(2,2), " dEdAlp=", work(2,3)
       CALL WrtOut(outputUnit,message)
       CALL WrtOut(outputUnit,"---------------------------------------------------------------------------")
        
       !===========================================
       ! Total energy change are close to zero
       IF ( nLine>=4 .and. ABS(work(1,2)-work(2,2))< 1.e-6 ) THEN 
         CALL WrtOut(outputUnit,"Ion-Optimization/CG2: energy changes < 1e-6 Ha (in 4 line searches). exit line search.")
         EXIT
       ENDIF

       !===========================================
       ! Displacement of atoms are close to zero
       IF (work(1,1) > 0.d0) THEN
         disTol = 1e-5  ! in bohr
         bracket_coord = (work(1,1)-work(2,1)) * dir
         IF ( MAXVAL(SQRT(SUM(bracket_coord**2,2))) < disTol ) THEN
!           restart_cg = .true.
!           alpha = 1.d0
           WRITE(message,'(a,ES12.4,a)') & 
            'Ion-Optimization/CG2: Maximum atom displacment is less than ', & 
             MAXVAL(SQRT(SUM(bracket_coord**2,2))), & 
             " bohr. Exit line search. " ! CG to be restarted."
           CALL WrtOut(6,message)
           EXIT
         ENDIF
       ENDIF
       
       !===============================================
       ! What's the message from line minimization ?
       IF(task(1:4) .EQ. 'CONV')  THEN
         WRITE(message,'(a,I2,a)') & 
           'Ion-Optimization/CG2: Line search is successful! ==> ', nLine, " steps."
         CALL WrtOut(6, message)
         EXIT

       ELSE IF(task(1:5) .EQ. 'ERROR') THEN
         WRITE(message,*) 'CG of ion optimization has ERROR: LnSrch()=>', task
         restart_cg = .true.
         alpha = 1.d0 ! Reset alpha
         CALL WrtOut(6,message)
         CALL WrtOut(6,"Ion-optimization: CG to be restarted.")
         EXIT

       ELSE IF( (task(1:2) .EQ. 'FG') .OR. (task(1:4) .EQ. 'WARN')) THEN 
         IF( task(1:4) .EQ. 'WARN') THEN 
           WRITE(message,*) 'NCG: WARNING: LnSrch()=>', task,' gradAlp:',gradAlp
           CALL WrtOut(6,message)
         ENDIF
         
         WRITE(message,*) " Try alpha=", alpha
         CALL WrtOut(6,message)
         ! need more line search, continue
         ! Compute trial phi with new theta
         DO ia=1,size(cell%iontable,1)
           !mohan test
           !WRITE(outputUnit,'(I5,A,3ES25.6)') ia, " move ", alpha*dir(ia,1), alpha*dir(ia,2), alpha*dir(ia,3)  
           delta_coord(ia,:) = alpha*dir(ia,:)
           cell%ionTable(ia)%coord(:) = MODULO( &
             MATMUL(cart2frac, & 
                    oldCoordCart(ia,:) + delta_coord(ia,:)) , 1._DP)
         ENDDO

         !=======================================================
         ! Compute gradOfAlp at new trial ion coordinates

         CALL RefreshIonTerms
         !DO ia=1,size(cell%iontable,1)
         !  write(*,*) cell%ionTable(ia)%coord(1)*cell%cellReal(1,1) + &
         !             cell%ionTable(ia)%coord(2)*cell%cellReal(1,2) + &
         !             cell%ionTable(ia)%coord(3)*cell%cellReal(1,3), &
         !             cell%ionTable(ia)%coord(1)*cell%cellReal(2,1) + &
         !             cell%ionTable(ia)%coord(2)*cell%cellReal(2,2) + &
         !             cell%ionTable(ia)%coord(3)*cell%cellReal(2,3), &
         !             cell%ionTable(ia)%coord(1)*cell%cellReal(3,1) + &
         !             cell%ionTable(ia)%coord(2)*cell%cellReal(3,2) + &
         !             cell%ionTable(ia)%coord(3)*cell%cellReal(3,3)
         ! Enddo

         CALL RhoOptimizer   ! Optimize the density to get total energy
         CALL CalculateForces(rho, forces)

         IF(ALLOCATED(frozenIon)) WHERE(frozenIon) forces(:,:,1)=0._DP
         CALL gradOfAlpha(forces,dir,size(cell%ionTable,1),gradAlp)

         !---- Check Point (2)-------
         CALL CheckExit(maxForce,relax_done)
         IF (relax_done) THEN 
           CALL WrtOut(6,"Exit during line search.")
           EXIT
         ENDIF
         

         ! Count line minimization #
         nLine = nLine + 1
         IF(nLine > maxLine) THEN
           restart_cg = .true.
           alpha = 1.d0
           CALL WrtOut(6,'Source/IonOptimizers.f90 ==> CG: WARNING:  &
                       & Exceeded the max iteration number for line search.&
                       & CG to be restarted')
           EXIT
         ENDIF
         
       END IF ! Test of Task: IF (task)

    END DO linesearch    
    !--- END OF LINE SEARCH ------------------------------
    !-----------------------------------------------------

  END DO CGloop

  CALL GeometryMinimizerReportSteps(& 
    step, TimerStop(watch), forces, maxForce, & 
    alpha, .FALSE., REAL(nLine,kind=DP),&
    .TRUE.)

  CALL GeometryMinimizerReportFooter(TimerStop(watch2))

  RETURN

CONTAINS


SUBROUTINE CheckExit(maximumForce, relaxDone)
!-------------------------------------------------------------------
! DESCRIPTION:
! Return true if the max(forces) has dropped below forceCutoff
! also print out the step information upon each call
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!-------------------------------------------------------------------
! REVISION LOG:
! by Chen Huang  (Aug/2009)
!-------------------------------------------------------------------

  USE IonOptimizers, ONLY: forceCutoff

  IMPLICIT NONE
  REAL(KIND=DP) :: maximumForce
  LOGICAL :: relaxDone
    
  maximumForce = MAXVAL(SQRT(SUM(forces(:,:,1)**2,2)))
  relaxDone = maximumForce < forceCutoff
  IF (relaxDone) THEN 
    CALL WrtOut(6,'(CG2): force drops below the criterion successfully, exit.')
    CALL WrtOut(6,'')
  ENDIF
  RETURN

END SUBROUTINE CheckExit


END SUBROUTINE ConjugateGradient


SUBROUTINE StepInfo(step, maxForce, stepSizeAB, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Print CG loop information
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Dec/2008  Created (Chen Huang)
!------------------------------------------------------------------------------

  USE RhoOptimizers, ONLY:  tole

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: step
  REAL(KIND=DP), INTENT(IN) :: maxForce
  REAL(KIND=DP), INTENT(IN) :: stepSizeAB, energy
  CHARACTER (LEN=500) :: message
  
  CALL WrtOut(6,'') 
  CALL WrtOut(6,'===================================================================================================')
  WRITE(message, '(A,I5,A,G24.16,A,Es12.4,A)') & 
              "(CG2): CG iter: ", step,' totEnergy=', energy, '(Ha)   maxForce=', maxForce, ' Ha/bohr'
  CALL WrtOut(6,message)
  WRITE(message, '(A,Es10.2,A,Es10.2,A)') & 
              "(CG2):              step  =", stepSizeAb, '      (bohr) tolE    =',tole, '   Ha'
  CALL WrtOut(6,message)
  CALL WrtOut(6,'===================================================================================================')
  CALL WrtOut(6,'') 

  RETURN

END SUBROUTINE StepInfo


SUBROUTINE gradOfAlpha(forces,gradient,numIons,gradAlp)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Compute the dE/d(alpha), i.e. the gradient of alpha w.r.t. total 
!   energy, used for line search in CG
!   for formula: 
!       dE/dAlpha = SUM(-Forces_i \dot gradient_i ) 
!
!   since gradient is in the same direction as forces
!   in the above formula will give a value < 0, and this is correct
!   since Alpha > 0 ( the step size in line search >0 )
!   Will decide if to restart the CG, gam will be set to 
!   zero if CG restart is needed
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Dec/2008  Created (Chen Huang)
!------------------------------------------------------------------------------
  IMPLICIT NONE

  REAL(KIND=DP), INTENT(IN) :: forces(:,:,:)  
  ! the forces at the given alpha
  !
  REAL(KIND=DP), INTENT(IN) :: gradient(SIZE(forces,1),SIZE(forces,2))  
  ! the dir for line search
  !
  INTEGER, INTENT(in) :: numIons  
  ! number of ions
  !
  REAL(KIND=DP), INTENT(OUT) :: gradAlp
  !
  INTEGER :: ia 
  !

  !>> FUNCTION << !

  gradAlp = 0._DP
  DO ia = 1, numIons
    gradAlp = gradAlp - dot_product(forces(ia,:,1),gradient(ia,:))
  ENDDO

  RETURN

END SUBROUTINE gradOfAlpha

END MODULE IonOptCG2

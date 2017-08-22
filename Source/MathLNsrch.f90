MODULE LnSrch_MOD
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE LnSrch_MOD
!     |_SUBROUTINE LnSrch 
!       |_FUNCTION newAlp 
!
! DESCRIPTION:
!  This is a line search module used for line search, written by Chen Huang
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
!   Aug/2009  Created (Chen Huang)
!   2014-07-01 Updated by Mohan
!------------------------------------------------------------------------------
  USE Constants, ONLY: DP
  USE MPI_Functions, ONLY: rankGlobal

  REAL(KIND=DP) :: dsave(13)
  REAL(KIND=DP) :: points(2,3)

CONTAINS

SUBROUTINE LnSrch(alp,f,g,task,c1,c2,xtol,stpMin,stpMax)
!------------------------------------------------------------------------------
! DESCRIPTION:
!  A home-made line search by Chen Huang, based on linear interpolation 
!  the g = df/d(alpha), in order to find the minimum.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/2009  Created (Chen Huang)
!------------------------------------------------------------------------------
  IMPLICIT NONE
  
  REAL(KIND=DP), INTENT(IN) :: f
  ! the function value at alp

  REAL(KIND=DP), INTENT(IN) :: g
  ! d(f)/d(alp) at alp
  !
  REAL(KIND=DP), INTENT(IN) :: stpMax, stpMin, xtol
  !
  REAL(KIND=DP) :: alp
  ! line search var
  !
  CHARACTER(LEN=70) :: task
  !
  REAL(KIND=DP) :: c1,c2
  !

  ! >> Internal vars << !!
  LOGICAL :: cond1, cond2
  !

  !! >> FUNCTION BEGINS <<!!
  IF(rankGlobal==0) THEN
    WRITE(6,*) "INFO: Current LnSrch does not respect xtol. Assigned value is: ",xtol
  ENDIF

  !-----------------------------
  ! Line search start
  !-----------------------------
  IF (task(1:5) == "START") THEN
    IF ( g>0.d0) THEN 
      task = "ERROR: Initial G>0"
      RETURN
    ENDIF
    points      = -1.d0
    points(1,1) = 0.d0
    points(1,2) = f
    points(1,3) = g
    alp         = alp
    dsave(1)    = f
    dsave(2)    = g
    points(2,1) = -1.d0
    task(1:2)   = "FG"
    RETURN
  ENDIF

  
  !--------------------------------------
  ! Test of Wolfe condition
  !--------------------------------------
  IF (task(1:2) == "FG" .and. points(2,1)>0.d0 ) THEN 
    cond1 = .false.
    cond2 = .false.
    !(1) sufficient decent
    cond1 = ( f .LE.  dsave(1) + c1*alp*dsave(2) )
    !(2) strong wolfe condition
    cond2 = ( ABS(g) .LE. c2*ABS(dsave(2)) )
    IF ( cond1 .and. cond2 ) THEN
      task = "CONV"
      RETURN
    ENDIF
  ENDIF


  !-----------------------------------------
  ! Comes the second point
  !-----------------------------------------
  IF ( points(2,1) < 0.d0) THEN 
    IF ( g<points(1,3)) THEN 
      IF ( f<points(1,2)) THEN 
        alp = alp*5.d0
      ELSE
        alp =(points(1,1) + alp)/2.d0
      ENDIF
      task = "WARNING: Need a better second point."
      RETURN 
    ENDIF
    points(2,1) = alp
    points(2,2) = f
    points(2,3) = g
    alp = newAlp(points(1,1),points(1,3),points(2,1),points(2,3),stpMin,stpMax)
    IF (alp < 0.d0 .OR. alp<points(1,1) )THEN
      IF( points(2,2) < points(1,2) )THEN 
        task = "WARNING: Interpolation error, code:-1, energy decreases, step further"
        alp = 5.d0 * points(2,1) 
      ELSE
        task = "ERROR: Interpolation error, code:-1"
      ENDIF
      RETURN
    ENDIF
    task = "FG"
    RETURN
  ENDIF
  
  !------------------------
  ! g = dE/d(alp) > 0 
  !------------------------
  IF (g > 0) THEN 
    IF ( points(2,3) > 0.d0 ) THEN
      IF (points(2,1) > alp) THEN 
        !! normal case
        ! Make right limit more tight
        points(2,1) = alp
        points(2,2) = f
        points(2,3) = g
        alp = newAlp(points(1,1),points(1,3),points(2,1),points(2,3),stpMin,stpMax)
        IF ( alp<points(1,1) ) THEN 
          !! wrong case
          alp = (points(1,1) + points(2,1))/2.d0
          task = "WARNING: Order flips, code:0, energy decreases, step further"
        ELSE
          !! normal case
          task = "FG"
        ENDIF
        RETURN
      ELSE  ! alp > points(2,1)
        IF ( f<points(2,2) ) THEN 
          task = "WARNINIG: Order flipped, code:1"
          alp = alp * 5.d0
        ELSE
          alp = (points(1,1) + points(2,1))/2.d0
          task = "WARNING: Order flipped, code:2"
        ENDIF
        RETURN
      ENDIF
    ENDIF

    IF ( points(2,3) < 0.d0 ) THEN
      !! make the brackets more tight
      points(1,1) = points(2,1)
      points(1,2) = points(2,2)
      points(1,3) = points(2,3)
      points(2,1) = alp
      points(2,2) = f
      points(2,3) = g
      alp = newAlp(points(1,1),points(1,3),points(2,1),points(2,3),stpMin,stpMax)
      IF (alp<points(1,1) .OR. alp>points(2,1) ) THEN 
        !! wrong case, 
        IF ( points(1,2) > points(2,2) ) THEN 
          alp = points(2,1) * 5.d0
          task = "WARNING: Interpolation error, code:1, energy decreases, step further"
        ELSE
          alp = (points(1,1) + points(2,1) ) /2.d0
          task = "WARNING: Interpolation error, code:1"
        ENDIF
      ELSE
        !! normal case
        task = "FG"
      ENDIF
      RETURN
    ENDIF
  ENDIF !! if (g>0)
   
  !-------------------
  ! g < 0 
  !-------------------
  IF ( g < 0.d0 ) THEN
    
    !! if alp is calculated to be wrong, due to in accurate g
    IF ( points(2,3) > 0.d0 ) THEN 
      !! update LEFT point
      points(1,1) = alp
      points(1,2) = f
      points(1,3) = g        
      alp = newAlp(points(1,1),points(1,3),points(2,1),points(2,3),stpMin,stpMax)
      IF (alp<points(1,1) .OR. alp>points(2,1))THEN
        !! wrong case
        IF ( points(1,2) > points(2,2) ) THEN 
          task = "WARNING: Interpolation, code:2, f decreases, step further"
          alp = 5.d0 * points(2,1) 
        ELSE
          alp = ( points(1,1) + points(2,1) ) /2.d0
          task = "WARNING: Interpolation, code:4, brackets the minimum by half."
        ENDIF
      ELSE
        !! normal case
        task = "FG"
      ENDIF
      RETURN
    ENDIF
    
    IF ( points(2,3) < 0.d0 ) THEN
      points(1,1) = points(2,1)
      points(1,2) = points(2,2)
      points(1,3) = points(2,3)
      points(2,1) = alp
      points(2,2) = f
      points(2,3) = g
      alp = newAlp(points(1,1),points(1,3),points(2,1),points(2,3),stpMin,stpMax)
      IF ( alp<points(2,1) )THEN
        !! Wrong case
        IF ( points(1,2) > points(2,2) ) THEN
          alp = points(2,1) * 5.d0
          task = "WARNING: Interpolation error: code:3, f decrease, step further"
        ELSE
          !! FATAL ERROR 
          task = "ERROR: Interpolation error: code:3"
        ENDIF
      ELSE
        !! Normal case: alp > points(2,1)
        task = "FG"
      ENDIF
      RETURN
    ENDIF

  ENDIF !! if (g < 0) 

  CONTAINS 

  FUNCTION newAlp(x1,y1,x2,y2,stpMin,stpMax)
  !------------------------------------------------------------------------------
  ! DESCRIPTION:
  !  find the newX to make newY = 0
  !  based on x_i and y_i
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
    REAL(KIND=DP) :: x1,x2,y1,y2
    ! two points for finding the y(x0)=0
    !
    REAL(KIND=DP) :: newAlp, stpMin, stpMax
    !
    REAL(KIND=DP) :: slope,b, alp_t

    slope = (y1-y2)/(x1-x2)
    b = y1 - x1*slope
    alp_t = -b/slope
    
    IF (alp_t>=stpMax .OR. alp_t<=stpMin) THEN 
      ! Only enlarge alp by 5 times
      newAlp = -1.d0
    ELSE
      newAlp = alp_t
    ENDIF

    RETURN
  END FUNCTION newAlp 

ENDSUBROUTINE LnSrch


END MODULE LnSrch_MOD

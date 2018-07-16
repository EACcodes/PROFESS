MODULE LocalPseudoPot
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE LocalPseudoPot
!     |_FUNCTION PseudoPotLookup
!       |_FUNCTION PseudoPotLookupReal
!       |_FUNCTION PseudoPotLookupRecip
!     |_FUNCTION PseudoPotDiffLookup
!
! DESCRIPTION:
!   This module contains universal data types for use in Atomic calculations.
!   For instance, It defines type for a pseudopotential, an ion, etc.  It also
!   defines a few functions to look up values given a numerical
!   pseudopotential.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!    
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2003-11-21  File Created (GSH)
!   2003-12-15  Added the Pseudopotential structures to this module (GSH)
!   2003-03-01  Really, really changed.  A lot. (GSH)
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP
  USE CONSTANTS, ONLY: PI
  USE MathSplines, ONLY: spline_cubic_val ! used for interpolate the cubic spline
  USE CellInfo

                              !>> GLOBAL <<!
  IMPLICIT NONE

CONTAINS


FUNCTION PseudoPotLookup(psp, lookUp)  
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function looks up a value in a numerically described pseudopotential
!   by some form of interpolation.  It can handle both real and reciprocal
!   space.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   7/11/2003  Created (Greg Ho)
!------------------------------------------------------------------------------
  IMPLICIT NONE


                            !>> EXTERNAL VARIABLES <<!
  TYPE(pseudoPot), INTENT(IN) :: psp
  ! The pseudopotntial passed in
  !
  REAL(kind=DP), INTENT(IN) :: lookup
  ! The norm of the q-vector (or r) to be looked up
  !
  REAL(kind=DP) :: PseudoPotLookup   
  ! The answer
  !

                            !>> INTERNAL VARIABLES <<!
                            !>> INITIALIZATION <<!
                             ! >> FUNCTION BODY <<!

  ! First for reciprocal space
  IF(psp%type==1) THEN
    PseudoPotLookup = PspRecipCubicSpline(psp, lookup)

  ! Then for real space
  ELSE IF(psp%type==2) THEN
    WRITE(*,*) "ERROR: PSEUDOPOTLOOKUP FOR REAL DISABLED"
    flush(6)
    stop
!    PseudoPotLookup = PseudoPotLookupReal(psp, lookup)

  ELSE
    STOP "PSP ERROR"

  END IF

  CONTAINS


  FUNCTION PspRecipCubicSpline(psp, qNorm)
  !----------------------------------------------------------------------------
  ! DESCRIPTION:
  !   This function takes a pseudopotential (let's call it V(q) 
  !   object and magnitude of a q-vector
  !   and looks up the value of the pseudopotential at that q-vector.  If the
  !   value is out of range, it returns 0. It will subtract Coulombic part (-4*piZ/q^2)
  !   from v(q), and then spline the non-Coulombic part. 
  !
  !   By doing this, we have a finite
  !   value at q=0, so we can spline the region near q=0. Otherwise,
  !   becaue V(q=0) = -infiity, we can NOT spline it. 
  !
  ! CONDITIONS AND ASSUMPTIONS:
  !
  ! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
  !
  !----------------------------------------------------------------------------
  ! REVISION LOG:
  !    Nov-28-2007, Created by Chen Huang
  !----------------------------------------------------------------------------

  IMPLICIT NONE
  TYPE(pseudoPot), INTENT(IN) :: psp
  ! The pseudopotntial passed in
  !
  REAL(KIND=DP), INTENT(IN) :: qNorm
  ! qNorm of each point in q-mesh, size = Nqx*Nqy*Nqz
  !
  REAL(KIND=DP) :: PspRecipCubicSpline
  !
  REAL(KIND=DP)  :: yval, &                  ! the wanted value at qNorm
                    ypp_left, &              ! the 1st derivative at left nearest point of qNorm
                    ypp_right,&              ! the 2nd derivaive at right nearest point of qNorm
                    y_left, &                ! vqS at left nearest point of qNorm
                    y_right,&                ! vqS at right nearest point of qNorm
                    pos,&                    ! the position of qNorm
                    qSpacing ,&              ! interval between nearest two q points in psp file
                    dt
             
  INTEGER :: left, &              ! left nearest point to qNorm
             right                ! right nearest point to qNorm
    
  !> function body <!
  qSpacing = psp%maxG / REAL((psp%numPoints - 1),KIND=DP)
 
  ! qNorm cannot be negative ...
  IF (qNorm<0._DP) THEN
    WRITE(*,*) 'In PspLookup(), qNorms < 0, code stop!'
    STOP
  ENDIF

  ! check some spectial cases
  IF (qNorm == 0._DP) THEN
    PspRecipCubicSpline = psp%potValues(1)
    RETURN
  ENDIF
          
  ! if qNorm bigger than defined maxG in psp, return 0._DP
  IF (qNorm >= psp%maxG) THEN
    PspRecipCubicSpline = 0._DP
    RETURN
  ENDIF

  ! read the following words to see why I do this check here
  ! since NaN is not equal to anyting in FORTRAN
  ! we use the following to see if 4*pi/qNorm is too big to case a NaN
  IF (-4._DP*pi/qNorm**2_DP .NE. -4_DP*pi/qNorm**2_DP ) THEN
    WRITE(*,*)&     
     ' There is another case which should be considered: very large bulk. ', &
     '  because larger system will give denser q points and will put q points ', &
     ' very closer to q=0 and this will make -4*pi*Z/q**2 to -infinity',&
     ' If computer find -4*pi*Z/q**2 is too big, it will generate NaN',&
     ' But currently, in our group, we have not meet a system large enough to ',&
     ' cause CPU to generate NaN, if you meet such kind of problem, do something!, code STOP!!'
    STOP
  ENDIF

  ! Not meet any special cases above, OK! let's do interpolation
  pos = qNorm/qSpacing + 1._DP

  ! evaluate the interpolation from cubic spline
  ! CALL spline_cubic_val ( psp%numPoints, psp%t, psp%vqS, psp%potDD, pos, yval, ypval, yppval )
  ! we do the interpolation directly here
  ! I just graped the code from spline_cubic_val subroutine here
  ! instead of calling it, this will be faster ?

  left = FLOOR(pos)
  right = left + 1
  !
  !  Evaluate the polynomial.
  !
  dt = pos - REAL(left, KIND=DP)
!
  y_left = psp%vqS(left)
  y_right = psp%vqS(right)
  ypp_left = psp%potDD(left)
  ypp_right = psp%potDD(right)

  yval = y_left + dt * ( y_right - y_left &
         - ( ypp_right / 6.0_DP + ypp_left / 3.0_DP ) &
         + dt * ( 0.5_DP * ypp_left &
         + dt * ( ( ypp_right - ypp_left ) / 6.0_DP ) ) )

  ! return value from spline_cubic_va is yval
  ! add back -4*pi*Z/q**2 term
  PspRecipCubicSpline = yval - 4._DP*pi*REAL(psp%charge, KIND=DP) / qNorm**2._DP
    
  RETURN
    
  END FUNCTION PspRecipCubicSpline

END FUNCTION PseudoPotLookup


FUNCTION PseudoPotDiffLookup(psp, qNorm)  
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function takes a pseudopotential object and magnitude of a q-vector
!   and computes the value of that pseudopotential's first derivative at that
!   q-vector. If the value is out of range, it returns 0. Otherwise it uses
!   a simple interpolation to estimate the numerical value returned.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   When q becomes very, very small, this method will fail: first, the q=0
!   value is so out of place that the point of index 2 cannot be reliably
!   differentiated. Second, should the point of index 1 have to be 
!   differentiated, the formula ends up out of the table without safeguards.
!   This can only happen with very large cells though.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/06/2004  Created based on PseudoPotLookup above.(Vincent Ligneres)
!   07/19/2004  Replaced the step functions with a linear interpolation of
!               the derivatives at each gridpoint, based on finite diff.
!               with the points immediately above and below. (VLL)
!   12/17/2007  Corrected by Chen Huang. here we use a standard spline procedure to 
!               calculate the derivative of V(q)-4*pi*Z/q**2.
!------------------------------------------------------------------------------
  IMPLICIT NONE

                            !>> EXTERNAL VARIABLES <<!

  TYPE(pseudoPot), INTENT(IN) :: psp
  ! The pseudopotential passed in
  !
  REAL(KIND=DP), INTENT(IN) :: qNorm
  ! The norm of the q-vector to be looked up
  !
  REAL(KIND=DP) :: PseudoPotDiffLookup  
  ! The answer!
  !

                            !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: & ! mohan add (KIND=DP) on 02-14-2013
    pspSpacing, &    ! The spacing between values of the pseudopotential.
    rawIndex, &      ! The projected index (may be fractional) at which we 
    yval, & 
    ypval, &
    yppval

                            !>> INITIALIZATION <<!

  IF(psp%type == 2) STOP "***PSP ERROR***"

  ! Calculate pspSpacing and rawIndex.  Add 1 to rawIndex because of q=0
  pspSpacing = psp%maxG / REAL(psp%numPoints - 1,KIND=DP)
  rawIndex = qNorm / pspSpacing + 1._DP

                             ! >> FUNCTION BODY <<!

  ! First take care of the possiblity that our q is out of range.
  IF (qNorm > psp%maxG) THEN
    PseudoPotDiffLookup = 0.0_DP
  ! Otherwise, do a pseudo-cubic-spline interpolation.
  ELSE 
     call spline_cubic_val( psp%numpoints, psp%t, psp%vqS, psp%potDD, rawIndex, yval, ypval, yppval )
     PseudoPotDiffLookup = ypval/pspSpacing + 8.0_DP*pi*psp%charge/qNorm**3
  END IF  
    
  RETURN
END FUNCTION PseudoPotDiffLookup


END MODULE LocalPseudoPot

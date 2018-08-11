MODULE CBSpline 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE CBSpline 
!     |_SUBROUTINE GetCardinalBSpline 
!     |_SUBROUTINE FillB
!     |_SUBROUTINE BSplineProduct
!
! DESCRIPTION:
! Used in spline method for ion-ion and ion-electron interactions
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES: 
!   Choly and Kaxiras, PRB 67, 155101
!   Essman, Perera, Berkowitz, Darden, Lee, and Pdersen, J Chem Phys 103, 8577
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/08/2007  Created (Linda Hung)
!------------------------------------------------------------------------------

  USE Constants, ONLY: DP
  USE Constants, ONLY: PI
  USE Constants, ONLY: IMAG
  USE OutputFiles

  IMPLICIT NONE

  COMPLEX(KIND=DP), DIMENSION(:), ALLOCATABLE :: bSpline1
  COMPLEX(KIND=DP), DIMENSION(:), ALLOCATABLE :: bSpline2
  COMPLEX(KIND=DP), DIMENSION(:), ALLOCATABLE :: bSpline3

  INTEGER :: splineOrder = 10       
  ! Order of cardinal b-splines (for II, IE terms)
  !

CONTAINS

SUBROUTINE GetCardinalBSpline(spline, order, x)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Solve for M_n(x) for integer values of x, 0<x<n, where M_n is the
!   nth order Cardinal B-spline.  Use iterative method with
!    M_n(x) = x/(n-1)*M_(n-1)(x) + (n-x)/(n-1)*M_(n-1)(x-1)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/08/2007  Created (Linda Hung)
!------------------------------------------------------------------------------
  IMPLICIT NONE
  
                          !>> EXTERNAL VARIABLES <<!

  INTEGER, INTENT(IN) :: order 
  ! order of the spline
  !
  REAL(kind=DP), INTENT(IN) :: x
  ! Value in interval (0,1] for which spline
  ! [0, M_n(x), M_n(x+1), ..., M_n(x-1)] will be generated
  !
  REAL(kind=DP), DIMENSION(0:), INTENT(OUT) :: spline
  ! Values outside of range 1:order will be 0

                          !>> INTERNAL VARIABLES <<!

  INTEGER :: n, k   
  ! Counters
  !

                          !>> INITIALIZATION <<!

  ! Initialize spline values to 0
  spline = 0._DP
  ! Initialize iterative method with order 2 spline
  spline(1) = x
  spline(2) = 1._DP - x

                     !>> FUNCTION BODY <<!

  DO n=3,order
    DO k=0,n-1
      spline(n-k) = ((x+REAL(n-k-1,kind=DP)) * spline(n-k) &
                     + (REAL(k+1,kind=DP)-x) * spline(n-k-1)) &
                    / REAL(n-1,kind=DP)
    END DO
  END DO

  RETURN

END SUBROUTINE GetCardinalBSpline



SUBROUTINE FillB(b1, b2, b3, order, totX, tot3, locOffset)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine is during module setup.  It is called only once, assuming 
!   fixed recip grid size/spacing and b-spline order
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/05/2007  Subroutine created.  (Linda Hung)
!   11/24/2008  Only vectors saved, instead of entire 3D array, to cut down
!               on memory usage (Linda Hung)
!------------------------------------------------------------------------------

  IMPLICIT NONE

                          !>> EXTERNAL VARIABLES <<!

  COMPLEX(KIND=DP), DIMENSION(:), INTENT(OUT) :: b1, b2, b3
  !
  INTEGER, INTENT(IN) :: order
  ! Order of Euler exponential spline
  !
  INTEGER, INTENT(IN) :: totX
  ! Number of gridpts for x-direction
  !
  INTEGER, INTENT(IN) :: tot3
  ! Tot # b-gridpts in 3rd recip space dimension
  !
  INTEGER, INTENT(IN) :: locOffset
  ! Offset value from 3rd recip space dimension
  !

                          !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(0:order) :: M 
  ! Spline array
  !
  COMPLEX(KIND=DP) :: denom
  INTEGER :: sizeX, recipDim2, recipDim3
  INTEGER :: ix, i2, i3, mValue, k

                          !>> INITIALIZATION <<!

  ! Obtain M_n(1), M_n(2), ... M_n(n-1)
  CALL GetCardinalBSpline(M, order, 1._DP)

  sizeX = SIZE(b1)
  recipDim2 = SIZE(b2)
  recipDim3 = SIZE(b3)
                          !>> FUNCTION BODY <<!

  DO ix = 1, sizeX
    denom = 0._DP
    mValue = ix - 1
    DO k = 0, order-2
      denom = denom + M(k+1)*EXP(2._DP*pi*imag*REAL(mValue,KIND=DP) &
                                 *REAL(k,KIND=DP)/REAL(totX,KIND=DP))
    END DO
    b1(ix) = EXP(2._DP*pi*imag*REAL(order-1,KIND=DP) &
                 *REAL(mValue,KIND=DP)/REAL(totX,KIND=DP)) &
              / denom
  END DO

  DO i2 = 1, recipDim2
    denom = 0._DP
    mValue = i2 - 1
    DO k = 0, order-2
      denom = denom + M(k+1)*EXP(2._DP*pi*imag*REAL(mValue,KIND=DP) &
                                 *REAL(k,KIND=DP)/REAL(recipDim2,KIND=DP))
    END DO
    b2(i2) = EXP(2._DP*pi*imag*REAL(order-1,KIND=DP) &
                 *REAL(mValue,KIND=DP)/REAL(recipDim2,KIND=DP)) &
              / denom
  END DO

  DO i3 = 1, recipDim3
    denom = 0._DP
    mValue = locOffset + i3 - 1
    DO k = 0, order-2
      denom = denom + M(k+1)*EXP(2._DP*pi*imag*REAL(mValue,KIND=DP) &
                                 *REAL(k,KIND=DP)/REAL(tot3,KIND=DP))
    END DO
    b3(i3) = EXP(2._DP*pi*imag*REAL(order-1,KIND=DP)&
                 *REAL(mValue,KIND=DP)/REAL(tot3,KIND=DP)) &
              / denom
  END DO

  RETURN

END SUBROUTINE FillB


SUBROUTINE BSplineProduct(array, results, squareProduct)
!------------------------------------------------------------------------------
! DESCRIPTION: 
!   Muliplies an array with b-splines, which were previously set up.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/08/2007  Created (Linda Hung)
!------------------------------------------------------------------------------

  IMPLICIT NONE

                          !>> EXTERNAL VARIABLES <<!

  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: array
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: results
  !
  LOGICAL, INTENT(IN), OPTIONAL :: squareProduct  
  ! Multiply by |bSpline|^2, instead of bSpline?
  !

                          !>> INTERNAL VARIABLES <<!

  INTEGER :: i, j, k   
  ! Counters
  !
  LOGICAL :: square
  !

                            !>> INITIALIZATION <<!

  square = .FALSE.

  IF(.NOT. ALLOCATED(bSpline1)) THEN
    message="SPLINE METHOD: Not allocated bSpline1 yet."
    CALL ERROR(6,message)
  ENDIF

  IF(.NOT. ALLOCATED(bSpline2)) THEN
    message="SPLINE METHOD: Not allocated bSpline2 yet."
    CALL ERROR(6,message)
  ENDIF

  IF(.NOT. ALLOCATED(bSpline3)) THEN
    message="SPLINE METHOD: Not allocated bSpline3 yet."
    CALL ERROR(6,message)
  ENDIF

                            !>> FUNCTION BODY <<!

  IF (PRESENT(squareProduct)) square = squareProduct

  IF (square) THEN
    DO i = 1, SIZE(array,1)
      DO j = 1, SIZE(array,2)
        DO k = 1, SIZE(array,3)
          results(i,j,k) = array(i,j,k) * &
                                  bSpline1(i) * CONJG(bSpline1(i)) * &
                                  bSpline2(j) * CONJG(bSpline2(j)) * &
                                  bSpline3(k) * CONJG(bSpline3(k))
        END DO
      END DO
    END DO

  ELSE
    DO i = 1, SIZE(array,1)
      DO j = 1, SIZE(array,2)
        DO k = 1, SIZE(array,3)
          results(i,j,k) = array(i,j,k) &
                                * bSpline1(i) * bSpline2(j) * bSpline3(k)
        END DO
      END DO
    END DO
  END IF

  RETURN

END SUBROUTINE BSplineProduct


END MODULE CBSpline

MODULE NMS
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE NMS 
!     |_SUBROUTINE pchez
!     |_SUBROUTINE pchfd
!     |_SUBROUTINE pchev
!     |_SUBROUTINE pchim
!     |_SUBROUTINE chfdv 
!     |_SUBROUTINE pchsp
!       |_FUNCTION pchdf
!     |_FUNCTION pchst
!
! DESCRIPTION:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
!------------------------------------------------------------------------------
! REVISION LOG:
! Created by Chen Huang and revised by Johannes and Mohan Chen
!------------------------------------------------------------------------------

  USE Constants, ONLY: DP

  PUBLIC :: pchez,pchev ! Profess only needs these two

CONTAINS


SUBROUTINE pchez ( n, x, f, d, spline, wk, lwk, ierr )
!*****************************************************************************80
!
!! PCHEZ carries out easy to use spline or cubic Hermite interpolation.
!
!  Discussion:
!
!    This routine sets derivatives for spline (two continuous derivatives)
!    or Hermite cubic (one continuous derivative) interpolation.
!    Spline interpolation is smoother, but may not "look" right if the
!    data contains both "steep" and "flat" sections.  Hermite cubics
!    can produce a "visually pleasing" and monotone interpolant to
!    monotone data.
!
!    This routine is an easy to use driver for the PCHIP routines.
!    Various boundary conditions are set to default values by PCHEZ.
!    Many other choices are available in the subroutines PCHIC,
!    PCHIM and PCHSP.
!
!    Use PCHEV to evaluate the resulting function and its derivative.
!
!    If SPLINE is TRUE, the interpolating spline satisfies the default
!    "not-a-knot" boundary condition, with a continuous third derivative
!    at X(2) and X(N-1).
!
!    If SPLINE is FALSE, the interpolating Hermite cubic will be monotone
!    if the input data is monotone.  Boundary conditions are computed from
!    the derivative of a local quadratic unless this alters monotonicity.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!    Carl deBoor,
!    A Practical Guide to Splines, Chapter IV,
!    Springer-Verlag,
!    New York, 1978.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  
!    N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent
!    variable values.
!
!    Input, real F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Output, real D(N), the derivative values at the data points.
!
!    Input, logical SPLINE, specifies if the interpolant is to be a spline
!    with two continuous derivatives (SPLINE is TRUE), or a Hermite cubic
!    interpolant with one continuous derivative (SPLINE is FALSE).
!
!    Workspace, real WK(LWK), required only if SPLINE is TRUE.
!
!    Input, integer LWK, the length of the work array WK, which 
!    must be at least 2*N.  However, WK is not needed if SPLINE is FALSE,
!    and in this case LWK is arbitrary.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    positive, can only occur when SPLINE is FALSE,  means that there were
!      IERR switches in the direction of monotonicity.  When SPLINE is
!      FALSE, PCHEZ guarantees that if the input data is monotone, the
!      interpolant will be too.  This warning is to alert you to the fact
!      that the input data was not monotone.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -7, if LWK is less than 2*N and SPLINE is TRUE.
!
  IMPLICIT NONE

  INTEGER :: lwk
  INTEGER :: n

  REAL(KIND=DP) ::  d(n)
  REAL(KIND=DP) ::  f(n)
  INTEGER , SAVE, DIMENSION ( 2 ) :: ic = (/ 0, 0 /)
  INTEGER :: ierr
  INTEGER , parameter :: incfd = 1
  LOGICAL :: spline
  REAL(KIND=DP) ::  vc(2)
  REAL(KIND=DP) ::  wk(lwk)
  REAL(KIND=DP) ::  x(n)

  if ( spline ) then
    CALL pchsp ( ic, vc, n, x, f, d, incfd, wk, lwk, ierr )
  else
    CALL pchim ( n, x, f, d, incfd, ierr )
  end if

  return
END SUBROUTINE pchez


SUBROUTINE pchfd ( n, x, f, d, incfd, skip, ne, xe, fe, de, ierr )
!*****************************************************************************80
!
!! PCHFD evaluates a piecewise cubic Hermite function.
!
!  Discsussion:
!
!    PCHFD evaluates a piecewise cubic Hermite function and its first
!    derivative at an array of points.  PCHFD may be used by itself
!    for Hermite interpolation, or as an evaluator for PCHIM
!    or PCHIC.
!
!    PCHFD evaluates the cubic Hermite function and its first derivative
!    at the points XE.
!
!    If only function values are required, use PCHFE instead.
!
!    To provide compatibility with PCHIM and PCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!    Most of the coding between the call to CHFDV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFDV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of PCHFD that assumes a strictly
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points less than
!    X(N-1), followed by points greater than X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XE(J)
!    implies
!      X(I) <= XE(K).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Modified:
!
!    13 August 2005
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  
!    N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent
!    variable values.
!
!    Input, real F(INCFD,N), the function values.
!    F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, real D(INCFD,N), the derivative values.
!    D(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Input, integer INCFD, increment between successive values in 
!    F and D.
!
!    Input/output, logical SKIP, controls whether data validity checks
!    should be made.  Setting the input value to FALSE will skip the checks.
!    On output with 0 <= IERR, SKIP will be set to TRUE.
!
!    Input, integer NE, the number of evaluation points.
!
!    Input, real XE(NE), points at which the function is
!    to be evaluated.
!
!    Output, real FE(NE), the values of the cubic Hermite
!    function at XE.
!
!    Output, real DE(NE), the derivative of the cubic
!    Hermite function at XE.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if NE < 1.
!    -5, if an error has occurred in the lower-level routine CHFDV.
!
  IMPLICIT NONE

  INTEGER ::  incfd
  INTEGER ::  n
  INTEGER ::  ne

  REAL(KIND=DP) ::  d(incfd,n)
  REAL(KIND=DP) ::  de(ne)
  REAL(KIND=DP) ::  f(incfd,n)
  REAL(KIND=DP) ::  fe(ne)
  INTEGER ::  i
  INTEGER ::  ierc
  INTEGER ::  ierr
  INTEGER ::  ir
  INTEGER ::  j
  INTEGER ::  j_first
  INTEGER ::  j_new
  INTEGER ::  j_save
  INTEGER ::  next(2)
  INTEGER ::  nj
  LOGICAL :: skip
  REAL(KIND=DP) ::  x(n)
  REAL(KIND=DP) ::  xe(ne)
!
!  Check arguments.
!
  if ( .not. skip ) then

    if ( n < 2 ) then
      ierr = -1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFD - Fatal error!'
      write ( *, '(a)' ) '  Number of data points less than 2.'
      return
    end if

    if ( incfd < 1 ) then
      ierr = -2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHFD - Fatal error!'
      write ( *, '(a)' ) '  Increment less than 1.'
      return
    end if

    do i = 2, n
      if ( x(i) <= x(i-1) ) then
        ierr = -3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFD - Fatal error!'
        write ( *, '(a)' ) '  X array not strictly increasing.'
        return
      end if
    end do

  end if

  if ( ne < 1 ) then
    ierr = -4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHFD - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points less than 1.'
    return
  end if

  ierr = 0
  skip = .true.
!
!  Loop over intervals.
!  The interval index is IL = IR - 1.
!  The interval is X(IL) <= X < X(IR).
!
  j_first = 1
  ir = 2

  do
!
!  Skip out of loop if have processed all evaluation points.
!
    if ( ne < j_first ) then
      exit
    end if
!
!  Locate all points in interval.
!
    j_save = ne + 1

    do j = j_first, ne
      if ( x(ir) <= xe(j) ) then
        j_save = j
        if ( ir == n ) then
          j_save = ne + 1
        end if
        exit
      end if
    end do
!
!  Have located first point beyond interval.
!
    j = j_save

    nj = j - j_first
!
!  Skip evaluation if no points in interval.
!
    if ( nj /= 0 ) then
!
!  Evaluate cubic at XE(J_FIRST:J-1).
!
      call chfdv ( x(ir-1), x(ir), f(1,ir-1), f(1,ir), d(1,ir-1), d(1,ir), &
        nj, xe(j_first), fe(j_first), de(j_first), next, ierc )

      if ( ierc < 0 ) then
        ierr = -5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PCHFD - Fatal error!'
        write ( *, '(a)' ) '  Error return from CHFDV.'
        return
      end if
!
!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
!
      if ( next(2) /= 0 ) then

        if ( ir < n ) then
          ierr = -5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PCHFD - Fatal error!'
          write ( *, '(a)' ) '  IR < N.'
          return
        end if
!
!  These are actually extrapolation points.
!
        ierr = ierr + next(2)

      end if
!
!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
!
      if ( next(1) /= 0 ) then
!
!  These are actually extrapolation points.
!
        if ( ir <= 2 ) then
          ierr = ierr + next(1)
!
!  XE is not ordered relative to X, so must adjust evaluation interval.
!
!  First, locate first point to left of X(IR-1).
!
        else

          j_new = -1

          do i = j_first, j-1
            if ( xe(i) < x(ir-1) ) then
              j_new = i
              exit
            end if
          end do

          if ( j_new == -1 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PCHFD - Fatal error!'
            write ( *, '(a)' ) '  Could not bracket the data point.'
            return
          end if
!
!  Reset J.  This will be the new J_FIRST.
!
          j = j_new
!
!  Now find out how far to back up in the X array.
!
          do i = 1, ir-1
            if ( xe(j) < x(i) ) then
              exit
            end if
          end do
!
!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
          ir = max ( 1, i-1 )

        end if

      end if

      j_first = j

    end if

    ir = ir + 1

    if ( n < ir ) then
      exit
    end if

  end do

  return
END SUBROUTINE pchfd


SUBROUTINE pchev ( n, x, f, d, nval, xval, fval, dval, ierr )
!*****************************************************************************80
!
!! PCHEV evaluates a piecewise cubic Hermite or spline function.
!
!  Discussion:
!
!    PCHEV evaluates the function and first derivative of a piecewise
!    cubic Hermite or spline function at an array of points XVAL.
!
!    The evaluation will be most efficient if the elements of XVAL are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XVAL(J)
!    implies
!      X(I) <= XVAL(K).
!
!    If any of the XVAL are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch,
!    Piecewise Cubic Hermite Interpolation Package, Final Specifications,
!    Lawrence Livermore National Laboratory,
!    Computer Documentation UCID-30194, August 1982.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  
!    N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent
!    variable values.
!
!    Input, real F(N), the function values.  F(I) is the value
!    corresponding to X(I).
!
!    Input, real D(N), the derivative values.  D(i) is the value
!    corresponding to X(I).
!
!    Input, integer NVAL, the number of points at which the 
!    functions are to be evaluated.
!
!    Input, real XVAL(NVAL), the points at which the functions
!    are to be evaluated.
!
!    Output, real FVAL(NVAL), the values of the cubic Hermite
!    function at XVAL.
!
!    Output, real DVAL(NVAL), the derivatives of the cubic
!    Hermite function at XVAL.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    positive, means that extrapolation was performed at IERR points.
!    -1, if N < 2.
!    -3, if the X array is not strictly increasing.
!    -4, if NVAL < 1.
!    -5, if an error has occurred in CHFDV.
!
  IMPLICIT NONE

  INTEGER :: n
  INTEGER :: nval

  REAL(KIND=DP) :: d(n)
  REAL(KIND=DP) :: dval(nval)
  REAL(KIND=DP) :: f(n)
  REAL(KIND=DP) :: fval(nval)
  INTEGER :: ierr
  INTEGER , SAVE :: incfd = 1
  LOGICAL, SAVE :: skip = .true.
  REAL(KIND=DP) :: x(n)
  REAL(KIND=DP) :: xval(nval)

  call pchfd ( n, x, f, d, incfd, skip, nval, xval, fval, dval, ierr )

  return
END SUBROUTINE pchev


SUBROUTINE pchim ( n, x, f, d, incfd, ierr )
!*****************************************************************************80
!
!! PCHIM sets derivatives for a piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    The routine set derivatives needed to determine a monotone piecewise
!    cubic Hermite interpolant to given data.
!
!    The interpolant will have an extremum at each point where
!    monotonicity switches direction.  See PCHIC if user control is desired
!    over boundary or switch conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!  Parameters:
!
!    Input, integer N, the number of data points.  
!    N must be at least 2.
!
!    Input, real X(N), the strictly increasing independent
!    variable values.
!
!    Input, real F(INCFD,N), dependent variable values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!    PCHIM is designed for monotonic data, but it will work for any F-array.
!    It will force extrema at points where monotonicity switches direction.
!    If some other treatment of switch points is desired, PCHIC should be
!    used instead.
!
!    Output, real D(INCFD,N), the derivative values at the
!    data points.  If the data are monotonic, these values will determine
!    a monotone cubic Hermite function.  The value corresponding to X(I)
!    is stored in D(1+(I-1)*INCFD).
!
!    Input, integer INCFD, the increment between successive
!    values in F and D.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    positive, means IERR switches in the direction of monotonicity
!    were detected.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!
  IMPLICIT NONE

  INTEGER :: incfd
  INTEGER :: n

  REAL(KIND=DP) :: d(incfd,n)
  REAL(KIND=DP) :: del1
  REAL(KIND=DP) :: del2
  REAL(KIND=DP) :: dmax
  REAL(KIND=DP) :: dmin
  REAL(KIND=DP) :: drat1
  REAL(KIND=DP) :: drat2
  REAL(KIND=DP) :: dsave
  REAL(KIND=DP) :: f(incfd,n)
  REAL(KIND=DP) :: h1
  REAL(KIND=DP) :: h2
  REAL(KIND=DP) :: hsum
  REAL(KIND=DP) :: hsumt3
  INTEGER :: i
  INTEGER :: ierr
  INTEGER :: nless1
  REAL(KIND=DP) :: temp
  REAL(KIND=DP) :: w1
  REAL(KIND=DP) :: w2
  REAL(KIND=DP) :: x(n)
!
!  Check the arguments.
!
  if ( n < 2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHIM - Fatal error!'
    write ( *, '(a)' ) '  Number of data points less than 2.'
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PCHIM - Fatal error!'
    write ( *, '(a)' ) '  Increment less than 1.'
    return
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCHIM - Fatal error!'
      write ( *, '(a)' ) '  X array not strictly increasing.'
      return
    end if
  end do

  ierr = 0
  nless1 = n - 1
  h1 = x(2) - x(1)
  del1 = ( f(1,2) - f(1,1) ) / h1
  dsave = del1
!
!  Special case N=2, use linear interpolation.
!
  if ( n == 2 ) then
    d(1,1) = del1
    d(1,n) = del1
    return
  end if
!
!  Normal case, 3 <= N.
!
  h2 = x(3) - x(2)
  del2 = ( f(1,3) - f(1,2) ) / h2
!
!  Set D(1) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  hsum = h1 + h2
  w1 = ( h1 + hsum ) / hsum
  w2 = -h1 / hsum
  d(1,1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,1), del1 ) <= 0.0D+00 ) then

    d(1,1) = 0.0D+00
!
!  Need do this check only if monotonicity switches.
!
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then

     dmax = 3.0D+00 * del1

     if ( abs ( dmax ) < abs ( d(1,1) ) ) then
       d(1,1) = dmax
     end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( 2 < i ) then
      h1 = h2
      h2 = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = ( f(1,i+1) - f(1,i) ) / h2
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0D+00

    temp = pchst ( del1, del2 )

    if ( temp < 0.0D+00 ) then

      ierr = ierr + 1
      dsave = del2
!
!  Count number of changes in direction of monotonicity.
!
    else if ( temp == 0.0D+00 ) then

      if ( del2 /= 0.0D+00 ) then
        if ( pchst ( dsave, del2 ) < 0.0D+00 ) then
          ierr = ierr + 1
        end if
        dsave = del2
      end if
!
!  Use Brodlie modification of Butland formula.
!
    else

      hsumt3 = 3.0D+00 * hsum
      w1 = ( hsum + h1 ) / hsumt3
      w2 = ( hsum + h2 ) / hsumt3
      dmax = max ( abs ( del1 ), abs ( del2 ) )
      dmin = min ( abs ( del1 ), abs ( del2 ) )
      drat1 = del1 / dmax
      drat2 = del2 / dmax
      d(1,i) = dmin / ( w1 * drat1 + w2 * drat2 )

    end if

  end do
!
!  Set D(N) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  w1 = -h2 / hsum
  w2 = ( h2 + hsum ) / hsum
  d(1,n) = w1 * del1 + w2 * del2

  if ( pchst ( d(1,n), del2 ) <= 0.0D+00 ) then
    d(1,n) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0D+00 * del2

    if ( abs ( dmax ) < abs ( d(1,n) ) ) then
      d(1,n) = dmax
    end if

  end if

  return
END SUBROUTINE pchim


SUBROUTINE chfdv ( x1, x2, f1, f2, d1, d2, ne, xe, fe, de, next, ierr )

!*****************************************************************************80
!
!! CHFDV evaluates a cubic polynomial and its derivative given in Hermite form.
!
!  Discussion:
!
!    CHFDV evaluates a cubic polynomial and its first derivative.
!    The cubic polynomial is given in Hermite form.  The evaluation
!    is carried out at an array of points.
!
!    This routine was designed for use by PCHFD, but it may also be
!    useful directly as an evaluator for a piecewise cubic Hermite
!    function in applications, such as graphing, where the interval
!    is known in advance.
!
!    If only function values are required, use CHFEV instead.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real X1, X2, the endpoints of the interval of
!    definition of  the cubic.  X1 and X2 must be distinct.
!
!    Input, real F1, F2, the values of the function at X1 and
!    X2, respectively.
!
!    Input, real D1, D2, the derivative values at the ends
!     of the interval.
!
!    Input, integer NE, the number of evaluation points.
!
!    Input, real XE(NE), the points at which the functions are to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in next.
!
!    Output, real FE(NE), DE(NE), the values of the cubic
!    function and its derivative at the points XE(*).
!
!    Output, integer NEXT(2), indicates the number of 
!    extrapolation points:
!    NEXT(1) = number of evaluation points to left of interval.
!    NEXT(2) = number of evaluation points to right of interval.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!
  IMPLICIT NONE

  INTEGER :: ne

  REAL(KIND=DP) :: c2
  REAL(KIND=DP) :: c2t2
  REAL(KIND=DP) :: c3
  REAL(KIND=DP) :: c3t3
  REAL(KIND=DP) :: d1
  REAL(KIND=DP) :: d2
  REAL(KIND=DP) :: de(ne)
  REAL(KIND=DP) :: del1
  REAL(KIND=DP) :: del2
  REAL(KIND=DP) :: delta
  REAL(KIND=DP) :: f1
  REAL(KIND=DP) :: f2
  REAL(KIND=DP) :: fe(ne)
  REAL(KIND=DP) :: h
  INTEGER :: i
  INTEGER :: ierr
  INTEGER :: next(2)
  REAL(KIND=DP) :: x
  REAL(KIND=DP) :: x1
  REAL(KIND=DP) :: x2
  REAL(KIND=DP) :: xe(ne)
  REAL(KIND=DP) :: xma
  REAL(KIND=DP) :: xmi
!
!  Check arguments.
!
  if ( ne < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The number of evaluation points was less than 1.'
    stop
  end if

  h = x2 - x1

  if ( h == 0.0D+00 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFDV - Fatal error!'
    write ( *, '(a)' ) '  The interval endpoints are equal.'
    return
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0D+00, h )
  xma = max ( 0.0D+00, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h

  c2 = -( del1 + del1 + del2 )
  c2t2 = c2 + c2
  c3 = ( del1 + del2 ) / h
  c3t3 = c3 + c3 + c3
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
    de(i) = d1 + x * ( c2t2 + x * c3t3 )
!
!  Count extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( xma < x ) then
      next(2) = next(2) + 1
    end if

  end do

  return
END SUBROUTINE chfdv


SUBROUTINE pchsp ( ic, vc, n, x, f, d, incfd, wk, nwk, ierr )
!*****************************************************************************80
!
!! PCHSP sets derivatives needed for Hermite cubic spline interpolant.
!
!  Description:
!
!    PCHSP sets derivatives needed to determine the Hermite representation
!    of the cubic spline interpolant to given data, with specified boundary
!    conditions.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F and D arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by PCHFE or PCHFD.
!
!    This is a modified version of Carl de Boor's cubic spline routine CUBSPL.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer-Verlag (new york, 1978).
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer IC(2), specifies desired boundary conditions:
!    IC(1) = IBEG, desired condition at beginning of data.
!    0, to set D(1) so that the third derivative is continuous at X(2).
!      This is the "not a knot" condition provided by de Boor's cubic spline
!      routine CUBSPL, and is the default boundary condition here.
!    1, if first derivative at X(1) is given in VC(1).
!    2, if second derivative at X(1) is given in VC(1).
!    3, to use the 3-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 3.
!    4, to use the 4-point difference formula for D(1).
!      Reverts to the default boundary condition if N < 4.
!    For the "natural" boundary condition, use ibeg=2 and vc(1)=0.
!    IC(2) = IEND, desired condition at end of data.
!    IEND may take on the same values as IBEG, but applied to derivative at
!    X(N).  In case IEND = 1 or 2, the value is given in VC(2).
!
!    Input, real VC(2), specifies desired boundary values,
!    as indicated above.  VC(1) need be set only if IC(1) = 1 or 2.
!    VC(2) need be set only if IC(2) = 1 or 2.
!
!    Input, integer N, the number of data points.  N must be
!    at least 2.
!
!    Input, real X(N), the strictly increasing independent
!    variable values.
!
!    Input, real F(INCFD,N), the dependent values to be
!    interpolated.  F(1+(I-1)*INCFD) is the value corresponding to X(I).
!
!    Output, real D(INCFD,N), the derivative values at the
!    data points.  These values will determine the cubic spline interpolant
!    with the requested boundary conditions.  The value corresponding to
!    X(I) is stored in D(1+(I-1)*INCFD).
!
!    Input, integer INCFD, increment between successive values
!    in F and D.
!
!    Workspace, real WK(NWK).
!
!    Input, integer NWK, the size of WK, which must be 
!    at least 2 * N.
!
!    Output, integer IERR, error flag.
!    0, no errors.
!    -1, if N < 2.
!    -2, if INCFD < 1.
!    -3, if the X array is not strictly increasing.
!    -4, if IBEG < 0 or 4 < IBEG.
!    -5, if IEND < 0 or 4 < IEND.
!    -6, if both of the above are true.
!    -7, if NWK is too small.
!    -8, in case of trouble solving the linear system
!        for the interior derivative values.
!
  IMPLICIT NONE

  INTEGER :: incfd
  INTEGER :: n

  REAL(KIND=DP) :: d(incfd,n)
  REAL(KIND=DP) :: f(incfd,n)
  REAL(KIND=DP) :: g
  INTEGER :: ibeg
  INTEGER :: ic(2)
  INTEGER :: iend
  INTEGER :: ierr
  INTEGER :: index
  INTEGER :: j
  INTEGER :: nwk
  REAL(KIND=DP) :: stemp(3)
  REAL(KIND=DP) :: vc(2)
  REAL(KIND=DP) :: wk(2,n)
  REAL(KIND=DP) :: x(n)
  REAL(KIND=DP) :: xtemp(4)

  if ( n < 2 ) then
    ierr = -1
    write(*,*) 'pchsp -- number of data points less than two', ierr
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    write(*,*) 'pchsp -- increment less than one', ierr
    return
  end if

  do j = 2, n
    if ( x(j) <= x(j-1) ) then
      ierr = -3
      write(*,*) 'pchsp -- x-array not strictly increasing', ierr
      return
    end if
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0

  if ( ibeg < 0 .or. 4 < ibeg ) then
    ierr = ierr - 1
  end if

  if ( iend < 0 .or. 4 < iend ) then
    ierr = ierr - 2
  end if

  if ( ierr < 0 ) then
    go to 5004
  end if
!
!  Function definition is ok -- go on.
!
  if ( nwk < 2 * n ) then
    go to 5007
  end if
!
!  Compute first differences of X sequence and store in wk(1,.). also,
!  compute first divided difference of data and store in wk(2,.).
!
  do j = 2, n
    wk(1,j) = x(j) - x(j-1)
    wk(2,j) = ( f(1,j) - f(1,j-1) ) / wk(1,j)
  end do
!
!  Set to default boundary conditions if N is too small.
!
  if ( n < ibeg ) then
    ibeg = 0
  end if

  if ( n < iend ) then
    iend = 0
  end if
!
!  Set up for boundary conditions.
!
  if ( ibeg == 1 .or. ibeg == 2 ) then
     d(1,1) = vc(1)
  else if ( 2 < ibeg ) then
!
!  Pick up first IBEG points, in reverse order.
!
     do j = 1, ibeg
       index = ibeg - j + 1
       xtemp(j) = x(index)
       if ( j < ibeg ) then
         stemp(j) = wk(2,index)
       end if
     end do

     d(1,1) = pchdf ( ibeg, xtemp, stemp, ierr )
     if ( ierr /= 0 ) then
       go to 5009
     end if

     ibeg = 1
  end if

  if ( iend == 1 .or. iend == 2 ) then
     d(1,n) = vc(2)
  else if ( 2 < iend ) then
!
!  Pick up last IEND points.
!
     do j = 1, iend
       index = n - iend + j
       xtemp(j) = x(index)
       if ( j < iend ) then
         stemp(j) = wk(2,index+1)
       end if
     end do

     d(1,n) = pchdf ( iend, xtemp, stemp, ierr )

     if ( ierr /= 0 ) then
       go to 5009
     end if

     iend = 1

  end if
!
!  Begin coding from cubspl
!
!  A tridiagonal linear system for the unknown slopes S(1:N) of
!  F at X(1:N) is generated and then solved by Gauss elimination,
!  with s(j) ending up in d(1,j), all j.
!  wk(1,.) and wk(2,.) are used for temporary storage.
!
!  Construct first equation from first boundary condition, of the form
!    wk(2,1) * s(1) + wk(1,1) * s(2) = D(1,1)
!
  if ( ibeg == 0 ) then
!
!  No condition at left end and N = 2.
!
     if ( n == 2 ) then
        wk(2,1) = 1.0D+00
        wk(1,1) = 1.0D+00
        d(1,1) = 2.0D+00 * wk(2,2)
!
!  Not-a-knot condition at left end and 2 < N.
!
     else
        wk(2,1) = wk(1,3)
        wk(1,1) = wk(1,2) + wk(1,3)
        d(1,1) =(( wk(1,2) + 2.0D+00 * wk(1,1) ) * wk(2,2) * wk(1,3) &
                             + wk(1,2)**2 * wk(2,3) ) / wk(1,1)
     end if
  else if ( ibeg == 1 ) then
!
!  Slope prescribed at left end.
!
     wk(2,1) = 1.0D+00
     wk(1,1) = 0.0D+00
  else
!
!  Second derivative prescribed at left end.
!
     wk(2,1) = 2.0D+00
     wk(1,1) = 1.0D+00
     d(1,1) = 3.0D+00 * wk(2,2) - 0.5D+00 * wk(1,2) * d(1,1)
  end if
!
!  If there are interior knots, generate the corresponding equations and
!  carry out the forward pass of Gauss elimination, after which the J-th
!  equation reads
!
!    wk(2,j) * s(j) + wk(1,j) * s(j+1) = d(1,j).
!
  if ( 1 < n-1 ) then
    do j = 2, n-1
        if ( wk(2,j-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -wk(1,j+1) / wk(2,j-1)
        d(1,j) = g * d(1,j-1) + 3.0D+00 &
          * ( wk(1,j) * wk(2,j+1) + wk(1,j+1) * wk(2,j) )
        wk(2,j) = g * wk(1,j-1) + 2.0D+00 * ( wk(1,j) + wk(1,j+1) )
    end do
  end if
!
!  Construct last equation from second boundary condition, of the form
!
!    (-g * wk(2,n-1)) * s(n-1) + wk(2,n) * s(n) = d(1,n)
!
!  If slope is prescribed at right end, one can go directly to back-
!  substitution, since arrays happen to be set up just right for it
!  at this point.
!
  if ( iend == 1 ) then
    go to 30
  end if

  if ( iend == 0 ) then
     if ( n == 2 .and. ibeg == 0 ) then
!
!  Not-a-knot at right endpoint and at left endpoint and N = 2.
!
        d(1,2) = wk(2,2)
        go to 30
     else if ( n == 2 .or. ( n == 3 .and. ibeg == 0 ) ) then
!
!  Either ( N = 3 and not-a-knot also at left) or (N=2 and *not*
!  not-a-knot at left end point).
!
        d(1,n) = 2.0D+00 * wk(2,n)
        wk(2,n) = 1.0D+00
        if ( wk(2,n-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -1.0D+00 / wk(2,n-1)
     else
!
!  Not-a-knot and 3 <= N, and either 3 < N or also not-a-
!  knot at left end point.
!
        g = wk(1,n-1) + wk(1,n)
!
!  Do not need to check following denominators (x-differences).
!
        d(1,n) = ( ( wk(1,n) + 2.0D+00 * g ) * wk(2,n) * wk(1,n-1) &
          + wk(1,n)**2 * ( f(1,n-1) - f(1,n-2) ) / wk(1,n-1) ) / g
        if ( wk(2,n-1) == 0.0D+00 ) then
          go to 5008
        end if
        g = -g / wk(2,n-1)
        wk(2,n) = wk(1,n-1)
     end if
  else
!
!  Second derivative prescribed at right endpoint.
!
     d(1,n) = 3.0D+00 *wk(2,n) + 0.5D+00 * wk(1,n) * d(1,n)
     wk(2,n) = 2.0D+00
     if ( wk(2,n-1) == 0.0D+00 ) then
       go to 5008
     end if
     g = -1.0D+00 / wk(2,n-1)
  end if
!
!  Complete forward pass of Gauss elimination.
!
  wk(2,n) = g * wk(1,n-1) + wk(2,n)

  if ( wk(2,n) == 0.0D+00 ) then
    go to 5008
  end if

  d(1,n) = ( g * d(1,n-1) + d(1,n) ) / wk(2,n)
!
!  Carry out back substitution.
!
   30 continue

  do j = n-1, 1, -1
    if ( wk(2,j) == 0.0D+00 ) then
      go to 5008
    end if
    d(1,j) = ( d(1,j) - wk(1,j) * d(1,j+1) ) / wk(2,j)
  end do

  return
!
!  error returns.
!
 5004 continue
!
!  ic out of range return.
!
  ierr = ierr - 3
  write(*,*) 'pchsp -- ic out of range', ierr
  return

 5007 continue
!
!  nwk too small return.
!
  ierr = -7
  write(*,*) 'pchsp -- work array too small', ierr
  return

 5008 continue
!
!  singular system.
!  theoretically, this can only occur if successive x-values
!  are equal, which should already have been caught (ierr=-3).
!
  ierr = -8
  write(*,*) 'pchsp -- singular linear system', ierr
  return
!
 5009 continue
!
!  error return from pchdf.
!  this case should never occur.
!
  ierr = -9
  write(*,*) 'pchsp -- error return from pchdf', ierr

  return

CONTAINS


FUNCTION pchdf ( k, x, s, ierr )
!*****************************************************************************80
!
!! PCHDF approximates a derivative using divided differences of data.
!
!  Discussion:
!
!    The routine uses a divided difference formulation to compute a K-point
!    approximation to the derivative at X(K) based on the data in X and S.
!
!    It is called by PCHCE and PCHSP to compute 3 and 4 point boundary
!    derivative approximations.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Pages 10-16,
!    Springer-Verlag, 1978.
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer K, is the order of the desired derivative 
!    approximation.  K must be at least 3.
!
!    Input, real X(K), contains the K values of the independent
!    variable.  X need not be ordered, but the values must be distinct.
!
!    Input/output, real S(K-1).  On input, the associated slope
!    values:
!      S(I) = ( F(I+1)-F(I))/(X(I+1)-X(I))
!    On output, S is overwritten.
!
!    Output, integer IERR, error flag.
!    0, no error.
!    -1, if K < 2.
!
!    Output, real PCHDF, the desired derivative approximation if
!    IERR=0 or to zero if IERR=-1.
!
  IMPLICIT NONE
  INTEGER :: k
  INTEGER :: i
  INTEGER :: ierr
  INTEGER :: j
  REAL(KIND=DP) :: pchdf
  REAL(KIND=DP) :: s(k-1)
  REAL(KIND=DP) :: value
  REAL(KIND=DP) :: x(k)
!
!  Check for legal value of K.
!
  if ( k < 3 ) then
    ierr = -1
    write(*,*) 'pchdf -- k less than three', ierr
    pchdf = 0.0D+00
    return
  end if
!
!  Compute coefficients of interpolating polynomial.
!
  do j = 2, k-1
    do i = 1, k-j
      s(i) = ( s(i+1) - s(i) ) / ( x(i+j) - x(i) )
    end do
  end do
!
!  Evaluate the derivative at X(K).
!
  value = s(1)

  do i = 2, k-1
    value = s(i) + value * ( x(k) - x(i) )
  end do

  ierr = 0
  pchdf = value

  return
END FUNCTION pchdf

END SUBROUTINE pchsp


FUNCTION pchst ( arg1, arg2 )

!*****************************************************************************80
!
!! PCHST: PCHIP sign-testing routine.
!
!  Discussion:
!
!    This routine essentially computes the sign of ARG1 * ARG2.
!
!    The object is to do this without multiplying ARG1 * ARG2, to avoid
!    possible over/underflow problems.
!
!  Author:
!
!    Fred Fritsch,
!    Mathematics and Statistics Division,
!    Lawrence Livermore National Laboratory.
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ARG1, ARG2, two values to check.
!
!    Output, real PCHST,
!    -1.0, if ARG1 and ARG2 are of opposite sign.
!     0.0, if either argument is zero.
!    +1.0, if ARG1 and ARG2 are of the same sign.
!
  IMPLICIT NONE

  REAL(KIND=DP) :: arg1
  REAL(KIND=DP) :: arg2
  REAL(KIND=DP) :: pchst

  pchst = sign ( 1.0D+00, arg1 ) * sign ( 1.0D+00, arg2 )

  if ( arg1 == 0.0D+00 .or. arg2 == 0.0D+00 ) then
    pchst = 0.0D+00
  end if

  RETURN
END FUNCTION pchst

END MODULE NMS

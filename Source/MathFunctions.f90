MODULE MathFunctions
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE MathFunctions
!     |_INTERFACE Cross 
!        |_FUNCTION CrossProductComplex
!        |_FUNCTION CrossProductReal
!     |_INTERFACE Norm
!        |_FUNCTION NormComplex     
!        |_FUNCTION NormReal
!     |_INTERFACE Det
!        |_FUNCTION DetComplex     
!        |_FUNCTION DetReal
!     |_INTERFACE Inverse
!        |_FUNCTION InverseComplex     
!        |_FUNCTION InverseReal  
!     |_FUNCTION LindG
!     |_FUNCTION GPrime
!     |_FUNCTION Volume
!
! DESCRIPTION:
!   This module includes some useful math formulas and equations.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!    
! REFERENCES:
!   Press, et. al.  Numerical Recipies, The Art of Scientific Computing.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/24/2003  File Created (Greg Ho)
!
!------------------------------------------------------------------------------
                                 !>> GLOBAL <<!  
  USE CONSTANTS, ONLY : &
  DP, &   ! Shortcut for "double precision"
  pi      ! come on, PI

  USE TIMER, ONLY : &
  stopWatch, &
  TimerStart, &
  TimerStop

  IMPLICIT NONE

  PRIVATE :: &
  CrossProductReal, &    ! use .X.
  CrossProductComplex, & ! use .X.
  DetReal, &             ! use Det
  DetComplex, &          ! use Det
  NormReal, &            ! use Norm 
  NormComplex, &         ! use Norm
  InverseReal, &         ! use Inverse
  InverseComplex         ! use Inverse

                            !>> INTERNAL VARIABLES <<!

  ! Cross product
  INTERFACE Cross 
    MODULE PROCEDURE CrossProductComplex
    MODULE PROCEDURE CrossProductReal
  END INTERFACE

  ! Gets the norm of a 3-vector
  INTERFACE Norm
    MODULE PROCEDURE NormReal
    MODULE PROCEDURE NormComplex
  END INTERFACE

  ! Gets the determinant of a matrix
  INTERFACE Det
    MODULE PROCEDURE DetReal
    MODULE PROCEDURE DetComplex
  END INTERFACE

  ! Inverts a matrix
  INTERFACE Inverse
    MODULE PROCEDURE InverseReal
    MODULE PROCEDURE InverseComplex
  END INTERFACE

CONTAINS


SUBROUTINE Uppercase(text)
!----------------------------------------------------------------------------
! DESCRIPTION:
!   This INTERNAL subroutine converts mixed case text into full uppercase 
!   text.  
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!   12/15/03 Moved this file originally in Datatypes written by VL to here
!            (GSH)
!--------------------------------------------------------------------------

  IMPLICIT NONE

                   !>> EXTERNAL VARIABLES <<!

  CHARACTER(len=*) :: text
  ! the string to modify.
  !
  INTEGER :: i
  ! dummy index for counting the characters in text.

                   !>> FUNCTION BODY <<!

    DO i=1, LEN(text)
      IF (text(i:i)>="a".AND.text(i:i)<="z") THEN
        text(i:i) = ACHAR(IACHAR(text(i:i))-IACHAR("a")+IACHAR("A"))
      END IF
    END DO

END SUBROUTINE Uppercase


FUNCTION CrossProductComplex(a,b)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Simply takes the cross product of 2 complex 3x3 vectors.
!
! CONDITIONS AND ASSUMPTIONS:
!   This should never be called directly, but the interface .CROSS. should be
!   used.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/30/2003  Function created.  (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                     !>> EXTERNAL VARIABLES <<!
  COMPLEX(KIND=DP), INTENT(IN), DIMENSION(3) :: a, b
  ! 2 vectors to be cross-producted
  !
  COMPLEX(KIND=DP), DIMENSION(3) :: CrossProductComplex    
  ! The answer
  !

                       !>> FUNCTION BODY <<!

  CrossProductComplex(1) = a(2)*b(3) - a(3)*b(2)
  CrossProductComplex(2) = a(3)*b(1) - a(1)*b(3)
  CrossProductComplex(3) = a(1)*b(2) - a(2)*b(1)

END FUNCTION CrossProductComplex


FUNCTION CrossProductReal(a,b)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Simply takes the cross product of 2 real 3x3 vectors.
!
! CONDITIONS AND ASSUMPTIONS:
!   This should never be called directly, but the interface .CROSS. should be
!   used.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/30/2003  Function created.  (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                     !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN), DIMENSION(3) :: a, b
  ! Two vectors to be cross-producted
  !
  REAL(KIND=DP), DIMENSION(3) :: CrossProductReal 
  ! The answer
  !

                       !>> FUNCTION BODY <<!
  CrossProductReal(1) = a(2)*b(3) - a(3)*b(2)
  CrossProductReal(2) = a(3)*b(1) - a(1)*b(3)
  CrossProductReal(3) = a(1)*b(2) - a(2)*b(1)

END FUNCTION CrossProductReal


FUNCTION NormComplex(a)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gives the norm of a complex 3-Vector
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  Function created.  (Greg Ho)
!   10/15/2013  move to intrinsic, as it should be!, and silence ifort (JMD)
!
!------------------------------------------------------------------------------

  IMPLICIT NONE
                     !>> EXTERNAL VARIABLES <<!

  COMPLEX(KIND=DP), INTENT(IN), DIMENSION(3) :: a
  ! The vector 
  !
  COMPLEX(KIND=DP) :: NormComplex          
  ! the answer
  !

  NormComplex = SQRT(dot_product(a,a))

END FUNCTION NormComplex


FUNCTION NormReal(a)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gives the norm of a real 3-Vector
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/30/2003  Function created.  (Greg Ho)
!   10/15/2013  move to intrinsic, as it should be!, and silence ifort (JMD)
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

                     !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN), DIMENSION(3) :: a
  ! The vector
  !
  REAL(KIND=DP) :: NormReal
  ! The answer
  !

  NormReal = SQRT(dot_product(a,a))

END FUNCTION NormReal


FUNCTION DetComplex(M) 
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gives the determinant of a 3x3 matrix (complex)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/21/2003  Function created.  (Vincent Ligneres)
!   11/20/2003  Changed name to DetComplex (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!  
  COMPLEX(KIND=DP), DIMENSION(3,3), INTENT(IN) :: M
  ! The matrix
  !
  COMPLEX(KIND=DP) :: DetComplex
  ! The answer
  !

  DetComplex = M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))&
               -M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1))&
               +M(1,3)*(M(2,1)*M(3,2)-M(3,1)*M(2,2))

END FUNCTION DetComplex


FUNCTION DetReal(M) 
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gives the determinant of a 3x3 matrix (real)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/21/2003  Function created.  (Vincent Ligneres)
!   11/20/2003  Changed name to DetReal (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!  
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: M
  ! The matrix
  !
  REAL(KIND=DP) :: DetReal
  ! The answer
  !

  DetReal = M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))&
               -M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1))&
               +M(1,3)*(M(2,1)*M(3,2)-M(3,1)*M(2,2))

END FUNCTION DetReal


FUNCTION InverseComplex(M) 
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function returns the inverse of a complex 3x3 matrix. It is the user's
!   responsibilty to make sure the matrix inverted is not singular. I do not
!   test for it here. You've been warned.
!
! CONDITIONS AND ASSUMPTIONS:
!   1. The matrix is not singular!!!
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  Created(Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!  
  COMPLEX(KIND=DP), DIMENSION(3,3), INTENT(IN) :: M
  ! the matrix
  !
  COMPLEX(KIND=DP), DIMENSION(3,3) :: InverseComplex    
  ! the answer
  !

                    !>> INTERNAL VARIABLES <<!  
  
  COMPLEX(KIND=DP) :: d
  ! d is the determinant of M.
  !
  INTEGER :: i, j
  ! Dummy indexes
  !

                     !>> INITIALIZATION <<!
  d = Det(M)
            
                     !>> FUNCTION BODY << 
  DO i=1, 3
    DO j=1,3
      InverseComplex(i,j)=(-1)**(REAL(i,KIND=DP)+REAL(j,KIND=DP))/d
    END DO 
  END DO
  InverseComplex(1,1)=InverseComplex(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))
  InverseComplex(2,1)=InverseComplex(2,1)*(M(1,2)*M(3,3)-M(1,3)*M(3,2))
  InverseComplex(3,1)=InverseComplex(3,1)*(M(1,2)*M(2,3)-M(2,2)*M(1,3))
  InverseComplex(1,2)=InverseComplex(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1))
  InverseComplex(2,2)=InverseComplex(2,2)*(M(1,1)*M(3,3)-M(1,3)*M(3,1))
  InverseComplex(3,2)=InverseComplex(3,2)*(M(1,1)*M(2,3)-M(2,1)*M(1,3))
  InverseComplex(1,3)=InverseComplex(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
  InverseComplex(2,3)=InverseComplex(2,3)*(M(1,1)*M(3,2)-M(3,1)*M(1,2))
  InverseComplex(3,3)=InverseComplex(3,3)*(M(1,1)*M(2,2)-M(2,1)*M(1,2))
  InverseComplex = TRANSPOSE(InverseComplex)

END FUNCTION InverseComplex


FUNCTION InverseReal(M) 
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function returns the inverse of a real 3x3 matrix. It is the user's
!   responsibilty to make sure the matrix inverted is not singular. I do not
!   test for it here. You've been warned.
!
! CONDITIONS AND ASSUMPTIONS:
!   1. The matrix is not singular!!!
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/21/2003  Created (Vincent Ligneres)
!   11/20/2003  Changed to be in MathFunctions Module, and cosmetically
!               modified (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!  
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: M
  ! the matrix
  !
  REAL(KIND=DP), DIMENSION(3,3) :: InverseReal
  ! the answer
  !

                    !>> INTERNAL VARIABLES <<!  
  
  REAL(KIND=DP) :: d
  ! d is the determinant of M.
  !
  INTEGER :: i, j
  ! Dummy indexes
  !

                     !>> INITIALIZATION <<!
  d = Det(M)
            
                     !>> FUNCTION BODY << 

  DO i=1, 3
    DO j=1,3
      InverseReal(i,j)=(-1)**(REAL(i,KIND=DP)+REAL(j,KIND=DP))/d
    END DO 
  END DO
  InverseReal(1,1)=InverseReal(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))
  InverseReal(2,1)=InverseReal(2,1)*(M(1,2)*M(3,3)-M(1,3)*M(3,2))
  InverseReal(3,1)=InverseReal(3,1)*(M(1,2)*M(2,3)-M(2,2)*M(1,3))
  InverseReal(1,2)=InverseReal(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1))
  InverseReal(2,2)=InverseReal(2,2)*(M(1,1)*M(3,3)-M(1,3)*M(3,1))
  InverseReal(3,2)=InverseReal(3,2)*(M(1,1)*M(2,3)-M(2,1)*M(1,3))
  InverseReal(1,3)=InverseReal(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
  InverseReal(2,3)=InverseReal(2,3)*(M(1,1)*M(3,2)-M(3,1)*M(1,2))
  InverseReal(3,3)=InverseReal(3,3)*(M(1,1)*M(2,2)-M(2,1)*M(1,2))
  InverseReal = TRANSPOSE(InverseReal)

END FUNCTION InverseReal


FUNCTION LindG(eta,lambda,mu) 
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the q-dependent part of WT kernel, as described
!   in [1] eq. (13) (and (11).)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Wang Y.A. Govind N. and Carter E.A., Phys. Rev. B 60(24), 1999, 16350
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/04/2003  Created (Vincent Ligneres)
!   10/17/2005  Modified to handle l * TF + mu * vW + WT-type functionals (VLL)
!   12/1/2005   Modified to do taylor expansion`at high eta (> 3.65) (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!  

  REAL(KIND=DP), INTENT(IN) :: eta
  ! the point at which the function gets evaluated.
  !
  REAL(KIND=DP), INTENT(IN) :: lambda
  ! the TF multiplier for compensating.
  !
  REAL(KIND=DP), INTENT(IN) :: mu
  ! the vW multiplier
  !

  REAL(KIND=DP) :: LindG
  ! The Lindhard G function as described in [1]
  !

                    !>> INTERNAL VARIABLES <<!  
  REAL(KIND=DP) :: eta2
  ! eta ** 2
  !
  REAL(KIND=DP) :: invEta2
  ! 1/eta**2
  !

                     !>> INITIALIZATION <<!
                     !>> FUNCTION BODY << 

  IF (eta<0._DP) THEN
    LindG = 0._DP

  ! Limit for small eta
  ELSE IF (eta < 1E-10_DP) THEN
    LindG = 1._DP - lambda + eta**2 * (1._DP / 3._DP - 3._DP * mu)

  ! Around the singularity
  ELSE IF (ABS(eta-1._DP) < 1E-10_DP) THEN
    LindG = 2._DP - lambda - 3._DP * mu + 20._DP * (eta-1._DP)

  ! Taylor expansion for high eta
  ELSE IF (eta > 3.65_DP) THEN ! we determined empircally that 3.65 was a 
                               ! good crossover point to the taylor expansion
    eta2 = eta**2
    invEta2 = 1._DP/eta2
    LindG = 3._DP*(1._DP-mu)*eta2 &
            -lambda-0.6_DP &
            + invEta2 * (-0.13714285714285712_DP &
            + invEta2 * (-6.39999999999999875E-2_DP &
            + invEta2 * (-3.77825602968460128E-2_DP &
            + invEta2 * (-2.51824061652633074E-2_DP &
            + invEta2 * (-1.80879839616166146E-2_DP &
            + invEta2 * (-1.36715733124818332E-2_DP &
            + invEta2 * (-1.07236045520990083E-2_DP &
            + invEta2 * (-8.65192783339199453E-3_DP & 
            + invEta2 * (-7.1372762502456763E-3_DP & 
            + invEta2 * (-5.9945117538835746E-3_DP & 
            + invEta2 * (-5.10997527675418131E-3_DP & 
            + invEta2 * (-4.41060829979912465E-3_DP & 
            + invEta2 * (-3.84763737842981233E-3_DP & 
            + invEta2 * (-3.38745061493813488E-3_DP & 
            + invEta2 * (-3.00624946457977689E-3_DP)))))))))))))))

  ELSE
    LindG = 1._DP / (0.5_DP + 0.25_DP * (1._DP-eta**2) * LOG((1._DP + eta)/&
             abs(1._DP-eta))/eta) - 3._DP * mu * eta**2 - lambda
  END IF

END FUNCTION LindG


FUNCTION GPrime(eta, mu) 
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the first derivative of the G lindhard function 
!   described in [1] eq. (13) (and (11).) Useful for WGC stress.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Wang Y.A. Govind N. and Carter E.A., Phys. Rev. B 60(24), 1999, 16350
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   05/13/2004  Created (Vincent Ligneres)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!  

  REAL(KIND=DP), INTENT(IN) :: eta
  ! the point at which the function gets evaluated.
  !
  REAL(KIND=DP), INTENT(IN) ::  mu
  ! vW multiplier.
  !
  REAL(KIND=DP) :: GPrime
  ! The first derivative of the function G from [1]
  !

                    !>> INTERNAL VARIABLES <<!  
                     !>> INITIALIZATION <<!
                     !>> FUNCTION BODY << 
  IF (eta<0._DP) THEN
    GPrime = 0._DP
  ELSE IF (eta < 1E-10_DP) THEN
    GPrime = 2._DP * eta * (1._DP / 3._DP - 3._DP * mu)
  ELSE IF (ABS(eta-1._DP) < 1E-10_DP) THEN
    GPrime = 40._DP
  ELSE
    GPrime = ((eta**2 + 1._DP) * 0.25_DP / eta**2 * LOG(ABS((1._DP + eta)/ &
             (1._DP-eta))) - 0.5_DP/eta) / (0.5_DP + 0.25_DP * (1._DP-eta**2) &
             * LOG((1._DP + eta) / ABS(1._DP-eta))/eta)**2 - 6._DP * eta * mu
  END IF

END FUNCTION GPrime


FUNCTION Volume(M)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the volume of a cell, when the cell is expressed
!   as a 3x3 matrix of its unit vectors (as is cellReal.)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!    
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/21/2003  File Created (Vincent Ligneres)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: M
  ! The 3x3 matrix of 3-D  unit vectors

  REAL(KIND=DP) :: Volume
  ! The answer we've all been waiting for
                
                          !>> INTERNAL VARIABLES <<! 
                            !>> INITIALIZATION <<!   
                           ! >> FUNCTION BODY <<!
  Volume = abs(Det(M))

END FUNCTION Volume


FUNCTION Vecmul(matrix,vector)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Multiplies a 3x3 matrix with a vector
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
  
                     !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: &
    matrix

  REAL(KIND=DP), DIMENSION(3), INTENT(IN) :: vector

  REAL(KIND=DP), DIMENSION(3) :: Vecmul

  Vecmul(1) = matrix(1,1)*vector(1)+matrix(1,2)*vector(2)+matrix(1,3)*vector(3)
  Vecmul(2) = matrix(2,1)*vector(1)+matrix(2,2)*vector(2)+matrix(2,3)*vector(3)
  Vecmul(3) = matrix(3,1)*vector(1)+matrix(3,2)*vector(2)+matrix(3,3)*vector(3)

END FUNCTION Vecmul


FUNCTION MinMaxVal(data, flag)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   take the data as input, the shape for the data is of (:,:,:), then get the 
!   maximum or minimum value depending on flag = 'MIN' or 'MAX'
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Sep/11/2008  Chen Huang
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

#ifdef __USE_PARALLEL
  INCLUDE 'mpif.h'
#endif

  REAL(KIND=DP), INTENT(IN) :: DATA(:,:,:)
  CHARACTER(LEN=*) :: flag
  real(KIND=DP) :: MinMaxVal, tmp
#ifdef __USE_PARALLEL
  real(KIND=DP) :: mpiGlobalVal
  integer ::  mpiErr
#endif
  CHARACTER(LEN=500) :: msg

  msg(1:3) = flag(1:3)
  CALL Uppercase(msg)

  if (msg(1:3)=='MIN') then
      tmp = MINVAL(data)
#ifdef __USE_PARALLEL
      CALL MPI_ALLREDUCE(tmp, mpiGlobalVal, 1, MPI_REAL8, MPI_MIN, &
                         MPI_COMM_WORLD, mpiErr)
      if (mpiErr /= MPI_SUCCESS) then
        STOP '(MinMaxval:) error in this subroutint due to MPI, code stop.'
      endif
      tmp = mpiGlobalVal
#endif
      MinMaxVal = tmp
  else if (msg(1:3)=='MAX') then
      tmp = MAXVAL(data)
#ifdef __USE_PARALLEL
      CALL MPI_ALLREDUCE(tmp, mpiGlobalVal, 1, MPI_REAL8, MPI_MAX, &
                         MPI_COMM_WORLD, mpiErr)
      if (mpiErr /= MPI_SUCCESS) then
        STOP '(MinMaxval:) error in this subroutint due to MPI, code stop.'
      endif
      tmp = mpiGlobalVal
#endif
      MinMaxVal = tmp
  else
      write(*,*) "ERROR: Called MinMaxVal w/o a proper message. Msg is ",msg
      flush(6)
      stop
  endif

  return


END FUNCTION MinMaxVal


SUBROUTINE Metric3d(cellReal, rmet, gmet)
!----------------------------------------------------------------------------
!  DESCRIPTION:
!    Compute the metric tensor of the cell.
!   input : cellreal = |a1  b1  c1|
!                      |a2  b2  c2|
!                      |a3  b3  c3|
!    a b c are the lattice vectors
!   
!    rmet = \sum_k ( cellreal_{k,i}*cellreal_{k,j} )  
!    gmet = \sum_k ( cellrecip_{k,i}*cellrecip_{k,j} )
! 
!   where cellrecip = inverse(transpose(cellreal))
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!   Dec/2008 by Chen Huang
!--------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=DP), INTENT(IN) :: & 
    cellReal(3,3)
  REAL(KIND=DP), INTENT(OUT) ::  & 
    rmet(3,3), & 
    gmet(3,3)

  REAL(KIND=DP) :: cellRecip(3,3)
  INTEGER :: ii 

  cellRecip = Inverse(Transpose(cellReal)) * 2.d0 * pi

  !Compute real space metrics
   do ii=1,3
     rmet(ii,:)=cellReal(1,ii)*cellReal(1,:)+&
                cellReal(2,ii)*cellReal(2,:)+&
                cellReal(3,ii)*cellReal(3,:)
   end do

  !Compute reciprocal space metrics
   do ii=1,3
     gmet(ii,:)=cellRecip(1,ii)*cellRecip(1,:)+&
                cellRecip(2,ii)*cellRecip(2,:)+&
                cellRecip(3,ii)*cellRecip(3,:)
   end do

END SUBROUTINE Metric3d


SUBROUTINE Diag(n,a,eig,eigv)
!----------------------------------------------------------------------------
!  DESCRIPTION:
!    Compute the eigen vector and eigen values of matrix a,
!    eigen value is stored eig(n)
!    eigen vectors are stored as eig(:,1/2/3.../n), COLUMN-WISE
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!   Dec/2008 by Chen Huang
!--------------------------------------------------------------------------
  implicit none
  integer :: n
  real(dp), intent(in) :: & 
    a(n,n)                    ! input matrix n*n
  real(dp), intent(out) :: & 
    eig(n), &                    ! output eigen values
    eigv(n,n)                    ! output eigen vectors, stored column-wise
  
  !! >> local vars
  real(dp) :: tmpA(n*(n+1)/2)
  real(dp) :: work(3*n)         ! working array, must of 3*n length
  integer  :: info,i,j

  ! first we need to zip A into tmpA, this is for 'U'
  ! based on tip from http://www.netlib.org/lapack/double/dspev.f
  do j=1,n
    do i=1,j
      tmpA(i+(j-1)*j/2) = A(i,j)
    enddo
  enddo

  ! DSPEV is a lapack subroutine, you need to link Lapack lib to use it
  call DSPEV('V','U',n,tmpA,eig,eigV,n,work,info) 
  if (info/=0) then
    print *,"Source/MathFunctions.f90>>Diag>> stop, error, info/=0, DSPEV() exit with problem!"
    stop
  endif

  return
END SUBROUTINE Diag


SUBROUTINE Symmetrize(n,a)
!----------------------------------------------------------------------------
!  DESCRIPTION:
!   Symmtrize incoming matrix a.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!   Dec/2008 by Chen Huang
!--------------------------------------------------------------------------
  implicit none 
  integer,intent(in) :: n
  real(dp),intent(inout) :: a(n,n)
  real(dp) :: b(n,n)
  integer :: ii,jj
  
  do ii=1,n
    do jj=1,n
      if (ii/=jj) then
        b(ii,jj) = (a(ii,jj)+a(jj,ii))/2._DP
      else
        b(ii,jj) = a(ii,jj)
      endif
    enddo
  enddo

  a = b

  return

END SUBROUTINE Symmetrize


FUNCTION stepfun(array)
!------------------------------------------------------------------------------
! Array is an input array, 
! stepfun(i,j,k)=1 if array(i,j,k)>=0, otherwise stepfun(i,j,k)=0
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!-------------------------------------------------------------------------------
! Created by Chen Huang (2008-Mar-16)
!--------------------------------------------------------------------------------
  IMPLICIT NONE

  REAL(KIND=DP), intent(in) :: array(:,:,:)
  REAL(KIND=DP) :: stepfun(size(array,1), size(array,2), size(array,3))

  stepfun = 1._DP
  WHERE (array < 0._DP) 
    stepfun = 0._DP
  END WHERE

  RETURN

END FUNCTION stepfun


FUNCTION GlobalVal3d(data, flag, useData)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   take the data as input, the shape for the data is of (:,:,:), then get the 
!   maximum or minimum value depending on flag = 'MIN' or 'MAX' or 'SUM'
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Sep/11/2008  Chen Huang
!   07/27/2010: Integrated sum option, added option of which points to count
!               (L. Hung)
!------------------------------------------------------------------------------
  IMPLICIT NONE

#ifdef __USE_PARALLEL
  INCLUDE 'mpif.h'
#endif

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    data        ! 3D array of real numbers

  CHARACTER(LEN=*), INTENT(IN) :: &
    flag        ! Find global minimum, maximum, or sum?

  LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: &
    useData     ! Whether to include data point in computing min/max/sum

  REAL(KIND=DP) :: &
    globalVal3d ! Result

                       !>> INTERNAL VARIABLES <<!

  CHARACTER(LEN=3) :: &
    msg

#ifdef __USE_PARALLEL
  INTEGER :: mpiErr
#endif

                           !>> INITIALIZATION <<!

  CALL StartClock('GlobalVal3d')

  msg = flag(1:3)
  CALL Uppercase(msg)

                           !>> FUNCTION BODY <<!

  IF(PRESENT(useData)) THEN
    IF (SIZE(useData,1)/=SIZE(data,1) .OR. &
        SIZE(useData,2)/=SIZE(data,2) .OR. &
        SIZE(useData,3)/=SIZE(data,3)) &
      STOP "GlobalVal3d: Mask and data arrays must have same dimension. Stop."


    SELECT CASE (msg)
#ifdef __USE_PARALLEL
      CASE ('MIN')
        CALL MPI_ALLREDUCE(MINVAL(data,useData),globalVal3d, 1, MPI_REAL8, &
                           MPI_MIN, MPI_COMM_WORLD, mpiErr)
      CASE ('MAX')
        CALL MPI_ALLREDUCE(MAXVAL(data,useData),globalVal3d, 1, MPI_REAL8, &
                           MPI_MAX, MPI_COMM_WORLD, mpiErr)
      CASE('SUM')
        CALL MPI_ALLREDUCE(SUM(data,useData),globalVal3d, 1, MPI_REAL8, &
                           MPI_SUM, MPI_COMM_WORLD, mpiErr)
#else
      CASE ('MIN')
        globalVal3d = MINVAL(data,useData)
      CASE ('MAX')
        globalVal3d = MAXVAL(data,useData)
      CASE('SUM')
        globalVal3d = SUM(data,useData)
#endif
    CASE DEFAULT
      STOP "GlobalVal3d: flag undefined, code stop."
    END SELECT

  ELSE ! There was no logical useDataIn array - check all points
    SELECT CASE (msg)
#ifdef __USE_PARALLEL
      CASE ('MIN')
        CALL MPI_ALLREDUCE(MINVAL(data),globalVal3d, 1, MPI_REAL8, MPI_MIN, &
                           MPI_COMM_WORLD, mpiErr)
      CASE ('MAX')
        CALL MPI_ALLREDUCE(MAXVAL(data),globalVal3d, 1, MPI_REAL8, MPI_MAX, &
                           MPI_COMM_WORLD, mpiErr)
      CASE('SUM')
        CALL MPI_ALLREDUCE(SUM(data),globalVal3d, 1, MPI_REAL8, MPI_SUM, &
                           MPI_COMM_WORLD, mpiErr)
#else
      CASE ('MIN')
        globalVal3d = MINVAL(data)
      CASE ('MAX')
        globalVal3d = MAXVAL(data)
      CASE('SUM')
        globalVal3d = SUM(data)
#endif
    CASE DEFAULT
      STOP "GlobalVal3d: flag undefined, code stop."
    END SELECT
  END IF

! Finally check that MPI call went okay
#ifdef __USE_PARALLEL
  IF(mpiErr /= MPI_SUCCESS) THEN
    STOP "GlobalVal3d: error in MPI_ALLREDUCE, code stop."
  ENDIF
#endif

  CALL StopClock('GlobalVal3d')

END FUNCTION GlobalVal3d

END MODULE MathFunctions

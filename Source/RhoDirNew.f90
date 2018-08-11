MODULE RhoDirNew
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoDirNew
!     |_SUBROUTINE NewtonDirection 
!
! DESCRIPTION:
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

CONTAINS
 

SUBROUTINE NewtonDirection(phi, b, mu0, iter, flag, isp, dimX, dimY, dimZ, nspin, dirNext)
!----------------------------------------------------------------------------
! DESCRIPTION:
!   We want to solve the linear equation Ax = -b where A is the Hessian of
!   the energy with respect to grid%rho, and b is the gradient of the energy
!   w.r.t. grid%rho.
!
! CONDITIONS AND ASSUMPTIONS:
!    only works for periodic boundary condition now
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!----------------------------------------------------------------------------
! REVISION LOG:
!   4/1/2004 Function created (GSH)
!   6/2/2008 Revised by Chen Huang
!            (1) no need for old incoming parameters: b2, tempGrid
!            (2) the expression for Ap is corrected
!----------------------------------------------------------------------------

    USE CONSTANTS, ONLY: DP, machPrec
    USE Hartree, ONLY : SqrtPreconditioner
    USE MPI_Functions, ONLY : ReduceRealLevel1
    USE CalPotPlus, ONLY: CalculatePotentialPlus 

    IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

    INTEGER, INTENT(IN) :: dimX, dimY, dimZ, nspin
    REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin), INTENT(IN) :: &
      b, &             ! The gradient: dL/dPhi = dE/dPhi - 2*mu*Phi
      phi              ! current wavefunction phi

    REAL(KIND=DP), INTENT(IN) :: mu0
    INTEGER, INTENT(IN) :: isp
    INTEGER, INTENT(INOUT) :: flag
    INTEGER, INTENT(OUT) :: iter
    REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ), INTENT(INOUT) :: dirNext ! The newton direction

                         !>> INTERNAL VARIABLES <<!

    REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: newPhi ! temperary phi
    REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: Ap ! The result of applying A (the Hessian) to p 
    REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ), TARGET :: &
      r, &             ! residual, r = b - Ax, also the gradient direction
      p, &             ! The next conjugate gradient direction
      y                ! The solution of My = r, where M is the preconditioner

    REAL(KIND=DP), DIMENSION(:,:,:), POINTER :: r3d

    REAL(KIND=DP) :: &
      mpiGlobalSum1, &
      mpiGlobalSum2, &
      criteria, &
      alpha, &
      pAp, &
      ryLast, &
      ry, &
      rr, &
      rrLast, &
      epsilon        ! The size of the finite difference step

    CALL Title("NewtonDirection") 
      
                           !>> INITIALIZATION <<!
    flag = -1
    r3d=>r

    ! Preconditioner related
    ! not implmented yet.
    ry = 0

                      !>> FUNCTION BODY <<!

    ! Start out with initial guess of 0 for the solution
    dirNext = 0._DP

    ! If this is the case, then r = -b - Ax
    r = -b(:,:,:,isp)


    rr = SUM(r**2)
    CALL ReduceRealLevel1(rr)

    rrLast = rr*2._DP
    criteria = rr

    iter = 1
    newPhi = phi
    !---------------------------------------------
    ! the inner CG loop starts here --------------
    !---------------------------------------------
    DO

      y = r  ! if no precond. r = -b ,so y is initilized to -b

      ryLast = ry

      ry = SUM(r*y)
      CALL ReduceRealLevel1(ry)

      IF(iter == 1) THEN
        ! Calculate the next direction, p.
        p = y
      ELSE
        ! Conjugate gradient part
        p = y + ry/ryLast * p
      END IF
      ! Find Ap, the result of applying the hessian A to vector p
      ! Do this by a 1st order finite differences approximation.

      ! Picking a good epsilon is sort of a tricky excercise.
      ! A straight epsilon = 1.E-9 seems to work pretty well, but here we have
      ! tried to tailor it to be dependent on the ratio of grid%rho and p.
      ! We find the minimum over all processors, so use mpi_allreduce
      ! mpiLocalsum is actually the local minimum, not the sum, but I didn't
      ! feel like creating another variable
      mpiGlobalSum1 = SUM(phi(:,:,:,isp)**2)
      CALL ReduceRealLevel1( mpiGlobalSum1 )
      mpiGlobalSum2 = SUM(p**2)
      CALL ReduceRealLevel1( mpiGlobalSum2 )
      epsilon = 2._DP*sqrt(machPrec)*(1._DP+SQRT(mpiGlobalSum1))/SQRT(mpiGlobalSum2)

      ! Calculate Ap with 1st order finite difference
      newPhi(:,:,:,isp) = phi(:,:,:,isp) + epsilon*p                     ! make the trial phi

      CALL CalculatePotentialPlus(newPhi**2, .TRUE., Ap) ! For Steven's KEDF

      Ap(:,:,:,isp) = Ap(:,:,:,isp)*sign(1._DP,newPhi(:,:,:,isp))                   ! convert to d(energy)/d(phi)
      Ap(:,:,:,isp) = Ap(:,:,:,isp) - mu0*2._DP*newPhi(:,:,:,isp)                   ! Now, we have done Ap at (rho+p*epsilon)
      Ap(:,:,:,isp) = (Ap(:,:,:,isp) - b(:,:,:,isp))/epsilon                        ! 1st order finite difference

      pAp = SUM(p*Ap(:,:,:,isp))
      CALL ReduceRealLevel1( pAp )

      ! Check for positive-definiteness here
      IF(pAp < 0._DP) THEN
        IF(iter == 1) THEN
          dirNext=-(dirNext-ry/pAp*p)
        END IF
        flag = -2
        EXIT
      END IF

      ! This is the exact distance we need to move in direction p to get to
      ! the minimum of energy (analytical "line-search")
      alpha = ry / pAp
      dirNext = dirNext + alpha * p
      ! Really, r = -b - Ax (or Axnew)  But, we don't have Ax.
      ! However, we do have Ap, and if we apply that to the equation
      ! xnew = xold + alpha*p, we get Axnew = Axold + alpha*Ap.
      ! Thus, r = -b - Ax_old - alpha*Ap.  But rold = -b - Ax_old.  So,
      ! rnew = rold - alpha*Ap
      r = r - alpha * Ap(:,:,:,isp)

      rrLast = rr
      rr = SUM(r**2)
      CALL ReduceRealLevel1(rr)

      ! Exit criteria. The residual has gone down by 10%, or we have past
      ! 50 iterations, or the residual has stopped changing.
      IF(rr < .1_DP*criteria) THEN
      ! the residual stops changing
!        WRITE(*,'(A)') 'the residue converged in NewtonDir function.'
        flag = 0
        EXIT
      ENDIF
      IF(iter > 50) THEN
!        WRITE(*,'(A)') '50 iteration steps has been done in NewtonDir function, exit.'
        flag = 1
        EXIT
      ENDIF
      IF((ABS(rr-rrLast)/rr < .01_DP .AND. iter>9)) THEN
!        WRITE(*,'(A)') 'the residual has dropped by less 10% from last iteration, exit the NewtonDir function'
        flag = 2
        EXIT
      ENDIF

      iter = iter + 1

   END DO
   
   RETURN

END SUBROUTINE NewtonDirection

END MODULE RhoDirNew

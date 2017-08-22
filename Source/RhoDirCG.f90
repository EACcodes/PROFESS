MODULE RhoDirCG
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoDirCG
!     |_SUBROUTINE CGDirection 
!     |_FUNCTION InnerProd
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

  USE Constants, ONLY: DP
  USE MPI_Functions

  IMPLICIT NONE

  CHARACTER(len=500) :: cg_alg = 'HZ'


CONTAINS


SUBROUTINE CGDirection(i,dimX,dimY,dimZ,nspin, &
                       dLdPhi,dLdPhi_old,dirNext,dirNext_old)
!----------------------------------------------------------------------------
! DESCRIPTION:
!  Compute the conjugate gradient direction for density optimization
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!----------------------------------------------------------------------------

  USE CellInfo, ONLY: cell

  IMPLICIT NONE

  !!>> EXTERNAL VARIABLES <<!!

  INTEGER , INTENT(IN) :: i   ! current CG steps
  INTEGER, INTENT(IN) :: dimX
  INTEGER, INTENT(IN) :: dimY
  INTEGER, INTENT(IN) :: dimZ
  INTEGER, INTENT(IN) :: nspin
  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: dLdPhi  
  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: dLdPhi_old 
  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: dirNext_old
  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: dirNext
  
  !!>> LOCAL VARIABLES << !!

  INTEGER :: isp
  REAL(KIND=DP), DIMENSION(nspin) :: cg_gamma 
  REAL(KIND=DP), DIMENSION(nspin) :: cg_eta
  REAL(KIND=DP), DIMENSION(nspin) :: norm_dk
  REAL(KIND=DP), DIMENSION(nspin) :: norm_y
  REAL(KIND=DP), DIMENSION(nspin) :: norm_gk
  REAL(KIND=DP), DIMENSION(nspin) :: dTy
  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin) :: y_k
             
  DO isp=1,nspin
    IF (i==1 .OR. MOD(i,100)==1) THEN  ! restart cg every 100 cg steps 
                                       ! to remove the accumulated cg errors
      cg_gamma(isp) = 0._DP
    ELSE
      SELECT CASE(cg_alg(1:2))
      CASE ('PR')
        ! Polak and Ribire formula
        cg_gamma(isp) = MAX(0._DP, & 
          InnerProd(dLdPhi(:,:,:,isp),dLdPhi(:,:,:,isp)-dLdPhi_old(:,:,:,isp)) / & 
          InnerProd(dLdPhi_old(:,:,:,isp),dLdPhi_old(:,:,:,isp)))
      CASE ('HZ')
        ! HAGER. & ZHANG. formula (2005) ACM
         y_k(:,:,:,isp) = dLdPhi(:,:,:,isp)-dLdPhi_old(:,:,:,isp)
         dTy(isp) = InnerProd(dirNext_old(:,:,:,isp),y_k(:,:,:,isp))
         norm_y(isp) = sqrt(InnerProd(y_k(:,:,:,isp),y_k(:,:,:,isp))/cell%dV)
         cg_gamma(isp) = 1._DP / dTy(isp) * &
             InnerProd(y_k(:,:,:,isp)-2._DP*dirNext_old(:,:,:,isp)*norm_y(isp)**2/dTy(isp), dLdPhi(:,:,:,isp))
         
         ! Hager and Zhang SIAM J. Optim. 2005. Eq(1.5)(1.6)
         cg_eta(isp) = 0.01_DP
         norm_dk(isp) = SQRT(InnerProd(dirNext_old(:,:,:,isp),dirNext_old(:,:,:,isp))/cell%dV)
         norm_gk(isp) = SQRT(InnerProd(dLdPhi_old(:,:,:,isp), dLdPhi_old(:,:,:,isp)) /cell%dV)
         cg_gamma(isp) = MAX(cg_gamma(isp), -1._DP/(norm_dk(isp)*MIN(cg_eta(isp),norm_gk(isp))))
      END SELECT
    ENDIF
    dLdPhi_old(:,:,:,isp) = dLdPhi(:,:,:,isp)
    dirNext(:,:,:,isp) = - dLdPhi(:,:,:,isp) + cg_gamma(isp) * dirNext_old(:,:,:,isp)
  ENDDO

  RETURN

END SUBROUTINE CGDirection


FUNCTION InnerProd(v1,v2)
!---------------------------------------------------------------
! DESCRIPTION:
! Compute the inner product of two vectors v1 and v2.
! v1 and v2 are three dimensional v1(x,y,z),v2(x,y,z)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! Created by Chen Huang, June-02-2008
!---------------------------------------------------------------

  USE CellInfo, ONLY: cell

  IMPLICIT NONE

  REAL(KIND=DP), INTENT(in) :: v1(:,:,:)
  REAL(KIND=DP), INTENT(in) :: v2(:,:,:)
  REAL(kind=DP) :: InnerProd

  !! >> FUNCTIONS << !!

  InnerProd = SUM(v1*v2)

  CALL ReduceRealLevel1(InnerProd)

  InnerProd = InnerProd * cell%dV

  RETURN

END FUNCTION InnerProd


END MODULE RhoDirCG

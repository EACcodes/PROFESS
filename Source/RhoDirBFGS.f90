MODULE RhoDirBFGS
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoDirBFGS
!     |_SUBROUTINE InitializeBFGS
!     |_SUBROUTINE CleanBFGS
!     |_FUNCTION BFGSDirection
!     |_SUBROUTINE UpdateBFGS
!     |_SUBROUTINE ChemicalPotential
!
! DESCRIPTION:
!   This modules provides BFGS algorithm to RhoTN.f90.
!
! REFERENCES:
! 1.  p779, Math. of Comp. 35,151, 1980
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   04/15/2011 Steven write the BFGS code. 
!   01/30/2013 Mohan put the BFGS into a new file RhoDirBFGS.f90
!------------------------------------------------------------------------------
                              !<< GLOBAL >>
  USE CONSTANTS, ONLY: DP
  USE OutputFiles
  USE MPI_Functions
  USE RhoDirCG, ONLY: InnerProd

  IMPLICIT NONE
  !
  INTEGER :: MBFGS = 5 
  ! the number of steps the L-BFGS method stores
  ! 
  INTEGER :: BFGSPT
  !
  INTEGER :: BFGSiter
  !
  REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: BFGSrho, BFGSalpha                              
  ! dimension of the two arrays are (MBFGS, NSPIN)
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: BFGSs, BFGSy
  ! work array for lbfgs subroutine.
  ! dimension are the same as rho

CONTAINS


SUBROUTINE InitializeBFGS(rho)
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

  IMPLICIT NONE

  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  INTEGER :: dimX, dimY, dimZ
  INTEGER :: nspin

  dimX  = SIZE(rho, 1)
  dimY  = SIZE(rho, 2)
  dimZ  = SIZE(rho, 3)
  nspin = SIZE(rho, 4)

  ! work array for lbfgs subroutin, scalar rho and alpha, see Math. of Comp. 35, 151, 1980
  IF( ALLOCATED( BFGSrho ) ) THEN
    !WRITE(outputUnit,*) "Have allocate BFGS."
  ELSE
    ALLOCATE( BFGSrho( MBFGS, nspin ) )
    ALLOCATE( BFGSalpha( MBFGS, nspin ) )
    ALLOCATE( BFGSs( MBFGS, dimX, dimY, dimZ, nspin ) )
    ALLOCATE( BFGSy( MBFGS, dimX, dimY, dimZ, nspin ) )
  ENDIF

  ! initialize 1 
  ! work variable for bfgs method
  BFGSrho = 0.d0
  BFGSalpha = 0.d0
  BFGSs = 0.d0
  BFGSy = 0.d0
  BFGSPT = 0
  BFGSiter = 0

  RETURN

END SUBROUTINE InitializeBFGS


SUBROUTINE CleanBFGS()
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

  IMPLICIT NONE
    
  DEALLOCATE( BFGSrho )
  DEALLOCATE( BFGSalpha )
  DEALLOCATE( BFGSs )
  DEALLOCATE( BFGSy )

  RETURN

END SUBROUTINE CleanBFGS


SUBROUTINE BFGSDirection(energy, phi,dLdPhi,dLdPhi_old,dV,dirNext)
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

  IMPLICIT NONE
  REAL(KIND=DP), INTENT(IN):: energy, dV
  REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:) :: phi, dLdPhi
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: dLdPhi_old
  REAL(KIND=DP),DIMENSION(SIZE(phi,1),SIZE(phi,2),SIZE(phi,3),SIZE(phi,4)) :: dirNext  
  ! Final result

  !! >> LOCAL VARIABLES << !!
  !work variables
  REAL(KIND=DP),DIMENSION(size(phi,1),size(phi,2),size(phi,3),size(phi,4)) :: Temp,Diag   
  !temp is a variable for temporary use and loop. Diag is H0, but only the diaganal part
  !
  REAL(KIND=DP),DIMENSION(size(phi,4)) :: Ys,YY,SQ,YR,Beta       
  !y*s and y**2 and s*Q, y*r
  !
  INTEGER :: NPT, CP, Bound, isp, i, NS   
  ! NPT is BFGSpoint-1, CP is for looping, Bound is for looping bound, NS=NumSpin
  !

  !! >> INITIALIZE << !!

  NS = size(phi,4)  
  Temp = 0.d0
  dirNext = 0.d0
  Bound = 0
  NPT = BFGSPT - 1
  IF(NPT == -1) NPT = MBFGS - 1
    
  !BFGSiter plus
  BFGSiter = BFGSiter + 1
  !if this is the first entry, record the old dLdPhi, then the direction is the minus dLdPhi, reture
  IF(BFGSiter == 1) THEN
    dirNext= -dLdPhi
    dLdPhi_old = dLdPhi
    RETURN
  ! if not the first entry, compute the direction described in paper
  ELSE
    Bound = BFGSiter-1
    IF(BFGSiter>MBFGS) THEN
      Bound = MBFGS
    ENDIF
        
    !calculate y*s,y*y and H0
    DO isp = 1,NS
      YS(isp) = InnerProd(BFGSs(NPT+1,:,:,:,isp),BFGSy(NPT+1,:,:,:,isp))
      YY(isp) = InnerProd(BFGSy(NPT+1,:,:,:,isp),BFGSy(NPT+1,:,:,:,isp))
      Diag(:,:,:,isp) = YS(isp)/YY(isp)
    ENDDO
        
    CP = BFGSPT
    IF(CP == 0) CP=MBFGS

    BFGSrho(CP,:) = 1/YS(:)  !compute BFGSrho
    Temp = -dLdPhi      !Qbound = Giter
    CP = BFGSPT

    ! now compute step 3 in H*G formulat, see p779, Math. of Comp. 35,151, 1980
    ! first compute Q0
    DO i=1,Bound
      CP = CP -1
        IF(CP == -1) CP = MBFGS - 1
        DO isp = 1,NS
          SQ(isp) = InnerProd(BFGSs(CP+1,:,:,:,isp),Temp(:,:,:,isp))  !S*Q
          BFGSalpha(CP+1,isp) = BFGSrho(CP+1,isp)*SQ(isp)            !compute alpha
          temp(:,:,:,isp) = temp(:,:,:,isp) - BFGSalpha(CP+1,isp)*BFGSy(CP+1,:,:,:,isp) !compute qi
        ENDDO  
    ENDDO ! finally get Q0

    temp = temp*Diag    !compute r0=H0*q0

    !now compute ri, from r0 to get r1,r2....finally riter
    DO i=1,Bound
      DO isp = 1,NS
        YR(isp) = InnerProd(BFGSy(CP+1,:,:,:,isp),temp(:,:,:,isp)) !YR=y*r
        beta(isp) = YR(isp)*BFGSrho(CP+1,isp) !beta=rho*y*r
        beta(isp) = BFGSalpha(CP+1,isp) - beta(isp) !compute alpha-beta
        temp(:,:,:,isp) = temp(:,:,:,isp) + beta(isp)*BFGSs(CP+1,:,:,:,isp)
      ENDDO
      CP = CP + 1
      IF (CP == MBFGS) CP = 0
    ENDDO !now we get the next direction

    dirNext = temp
    dLdPhi_old = dLdPhi

  ENDIF ! END BFGS ITER

  RETURN
    
END SUBROUTINE BFGSDirection    
  

SUBROUTINE UpdateBFGS(phi, tempPhi, dEdPhi, dLdPhi, totEleNum, dV)
!------------------------------------------------------------------------------
! DESCRIPTION: 
! This is a subroutine for lbfgs method. since we don't use the 
! line search in lbfgs subroutine, we should tell the W work array
! our new line search direction.
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
 
  REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:) :: phi
  REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:) :: tempPhi
  REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:) :: dEdPhi
  REAL(KIND=DP), INTENT(IN), DIMENSION(:,:,:,:) :: dLdPhi

  REAL(KIND=DP), INTENT(IN) :: totEleNum(:)
  ! Total number of electrons in the system
  !
  REAL(KIND=DP), INTENT(IN) :: dV
  ! Volume element on the real grid point
  !
    
  !work variables
  REAL(KIND=DP), DIMENSION(size(phi,1),size(phi,2),size(phi,3),size(phi,4)) :: tempdLdPhi !tempdLdphi

  REAL(KIND=DP), DIMENSION(size(phi,4)) :: mu 

  INTEGER :: isp,NS


  !! >> INITIALIZE << !!
  CALL Title("RhoOptMethodBFGS::UpdateBFGS")
  CALL StartClock('UpdateBFGS')
  tempdLdPhi = 0.d0
  mu = 0.d0
  NS = size(phi,4)
    
  !! >> FUNCTION << !!
  ! update chemical potential and the dL/dPhi with new phi
  CALL ChemicalPotential(dEdPhi,tempPhi,totEleNum,mu)

  DO isp=1,size(phi,4)
    tempdLdPhi(:,:,:,isp)  = dEdPhi(:,:,:,isp) - 2._DP*mu(isp)*tempPhi(:,:,:,isp)
  ENDDO
 
  !now update BFGSs and BFGSy, etc
  DO isp = 1, NS
    BFGSs(BFGSPT+1,:,:,:,isp) = tempPhi(:,:,:,isp)-phi(:,:,:,isp)
    BFGSy(BFGSPT+1,:,:,:,isp) = tempdLdPhi(:,:,:,isp)-dLdPhi(:,:,:,isp)
  ENDDO
    
  BFGSPT = BFGSPT + 1
  IF(BFGSPT == MBFGS) BFGSPT = 0

  CALL StopClock('UpdateBFGS')

  RETURN

END SUBROUTINE UpdateBFGS


SUBROUTINE ChemicalPotential(dEdPhi,phi,totEle,mu0)
!-------------------------------------------------------------------------------
! DESCRIPTION:
! This subroutine updates the chemical potential, also called
! the Lagrangian multiplier in the formula:
!      mu0 = \int(dEdPhi*phi)/(2*N_e).
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! June-1-2008: Created By Chen Huang 
!
!-------------------------------------------------------------------------------

  USE CONSTANTS, ONLY : DP
  USE CellInfo, ONLY: cell, numSpin
    
  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!

  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: dEdPhi
  ! d(energy)/d(phi)
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: phi
  ! the wavefunction phi
  !
  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: totEle
  ! total electron for each spin channel 
  !
  REAL(KIND=DP), DIMENSION(:), INTENT(INOUT) :: mu0
  !

  !! >> INTERNAL VARIABLES << !!

  INTEGER :: isp

  !! >> FUNCTION << !!

  DO isp = 1, numSpin
    mu0(isp)= cell%dV * SUM(dEdPhi(:,:,:,isp)*phi(:,:,:,isp))/(2._DP*totEle(isp))
    CALL ReduceRealLevel1(mu0(isp))
  ENDDO

  RETURN

END SUBROUTINE ChemicalPotential


END MODULE RhoDirBFGS

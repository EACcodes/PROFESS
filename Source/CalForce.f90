MODULE CalForces
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE CalForces
!     |_SUBROUTINE CalculateForces
!     |_SUBROUTINE ForcesII_IE
!     |_SUBROUTINE ForcesDenDec
!
! DESCRIPTION:
! Calculate forces, including Ion-Ion force (Ewald term),
! Ion-Electron force (local pseudopotential term), 
! and should include other terms that related to local charge.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! force should be allocated in a better array
! force should be gathered more nicely
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014-06-04 Mohan revised.
!------------------------------------------------------------------------------
                              !<< GLOBAL >>
  USE CONSTANTS, ONLY : DP
  USE CellInfo, ONLY: cell
  USE Fourier, ONLY: FFT
  USE MPI_Functions

  IMPLICIT NONE

CONTAINS

SUBROUTINE CalculateForces(rhoR, forces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine calculates the forces in a system given the density and the
!   ionic configuration.  The answer is stored in FORCES.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!  11/19/2003 Created (Greg Ho)
!  2013-08-18 Created by Mohan Chen
!------------------------------------------------------------------------------

  USE KEDF_DenDec, ONLY: do_den_dec

  IMPLICIT NONE

  !! >> External variables << !!
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: rhoR  ! Electronic density in real-space
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(INOUT) :: forces  ! forces

  CALL TITLE("CalForces::CalculateForces")
  CALL StartClock("CalForce")

  IF(do_den_dec .EQ. 1) THEN
    CALL ForcesDenDec(rhoR, forces)
  ELSE
    CALL ForcesII_IE(rhoR, forces)
  ENDIF
  
  CALL StopClock("CalForce")

  RETURN

END SUBROUTINE CalculateForces


SUBROUTINE ForcesDenDec(rhoR, forces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Calculate force in density decomposition method.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   2013-08-18 Created by Mohan Chen 
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell, numSpin
  USE CellInfo, ONLY: n1G, n2G, n3G
  USE KEDF_DenDec, ONLY: core_den
  USE KEDF_DenDec, ONLY: PulayForceDenDec 
  USE KEDF_DenDec, ONLY: ForceDenDec 

  IMPLICIT NONE

  !! >> External variables << !!
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: rhoR  ! Electronic density in real-space
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(INOUT) :: forces  ! total forces

  !! >> Internal variables << !!
  REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: rhoTotal ! for total density
  ! decomposition method

  INTEGER :: is ! index for spin

  !! >> Initialization << !!
  ALLOCATE(rhoTotal(n1G,n2G,n3G,numSpin)) 

  !! >> Function << !!
  
  ! sum up the total density
  DO is=1, numSpin
    rhoTotal(:,:,:,is) = rhoR(:,:,:,is) + core_den(:,:,:)
  ENDDO

  CALL ForcesII_IE(rhoTotal, forces)  ! Compute initial forces


  ! The second part is the Pulay forces
  ! \int grad(core_den) * (V_local + V_Hxc + a*V_TF + b*V_vW) dr
  !CALL PulayForceDenDec(forces(:,:,1)) ! 1 means total force
  CALL ForceDenDec(forces(:,:,1)) ! 1 means total force

  DEALLOCATE(rhoTotal)

  RETURN

END SUBROUTINE ForcesDenDec



SUBROUTINE ForcesII_IE(rhoR, forces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine calculates the forces in a system given the density and the
!   ionic configuration.  The answer is stored in FORCES.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/19/2003  Created (Greg Ho)
!------------------------------------------------------------------------------

  USE Ewald, ONLY : EwaldForces 
  USE IonElectron, ONLY: IonElectronForces
  USE IonElectronSpline, ONLY: IonElectronForcesSpline
  USE IonElectronSpline, ONLY: iiSpline, ieSpline
  USE OUTPUT, ONLY: WrtOut

  IMPLICIT NONE

                  !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: forces ! Total forces
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rhoR ! The electron density in real space

                  !>> INTERNAL VARIABLES <<!

#ifdef __USE_PARALLEL
  REAL(KIND=DP), DIMENSION(SIZE(forces,1),SIZE(forces,2)) :: tempLocForces ! Array to store forces before MPI summation
#endif

                  !>> INITIALIZATION <<!
  forces = 0._DP
                  ! >> FUNCTION BODY <<!

  ! The first component of forces is the total forces.  The second component
  ! is that ion-ion, the third is the ion-electron. 
  ! We store the individual components because we hope there aren't really that many ions so it 
  ! shouldn't be limiting in terms of memory, but its good to be able to have
  ! this sort of information.
  forces(:,:,2) = EwaldForces(cell%cellReal, cell%ionTable, cell%elementTable, iiSpline) 

  IF(ieSpline) THEN  ! IE spline approx
    forces(:,:,3) = IonElectronForcesSpline(FFT(SUM(rhoR,4)), cell%ionTable, cell%elementTable, cell%cellReal)
  ELSE               ! IE exact
    forces(:,:,3) = IonElectronForces(FFT(SUM(rhoR,4)), cell%ionTable, cell%elementTable, cell%cellReal)
  END IF

  ! If we're running in parallel, ion-electron interactions from separate
  ! processors still need to be added to get total forces for all ions.
#ifdef __USE_PARALLEL
  tempLocForces = forces(:,:,3)
  CALL MPI_ALLREDUCE(tempLocForces, forces(:,:,3), SIZE(tempLocForces), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) STOP "**PROBLEMS WITH MPI_ALLREDUCE IN CALCULATEFORCES***"
#else
#endif


!!!!!!!!!!!!!!!! BEGIN DEBUGGING OUTPUT
!  ! Conserve momentum in total system (remove average force in each dimension)
!  DO i_dim=1,3
!    avg = SUM(forces(:,i_dim,2))/REAL(SIZE(cell%ionTable),kind=DP)
!    forces(:,i_dim,2) = forces(:,i_dim,2) - avg
!    WRITE(*,*) i_dim, "dim avg ewald force was", avg
!  END DO
!  ! Conserve momentum in total system (remove average force in each dimension)
!  DO i_dim=1,3
!    avg = SUM(forces(:,i_dim,3))/REAL(SIZE(cell%ionTable),kind=DP)
!    forces(:,i_dim,3) = forces(:,i_dim,3) - avg
!    WRITE(*,*) i_dim, "dim avg i-e force was", avg
!  END DO
!
!WRITE(*,*) "ewald forces:"
!DO i_dim=1,SIZE(forces,1)
!WRITE(*,*) forces(i_dim,:,2)
!END DO
!WRITE(*,*) "i-e forces:"
!DO i_dim=1,SIZE(forces,1)
!WRITE(*,*) forces(i_dim,:,3)
!END DO
!CALL F77flush()
!!!!!!!!!!!!!!!! END DEBUGGING OUTPUT

!-----------------------------------------------------------
! forces(:,:,2) is ion-ion forces
! forces(:,:,3) is ion-electron forces
!------------------------------------------------------------
  forces(:,:,1) = forces(:,:,2) + forces(:,:,3)


  ! Conserve momentum in total system (remove average force in each dimension)

  ! *hartreeToeV/bohr
!  WRITE(outputUnit,*) " In order to conserve the momentum in total system."
!  DO i_dim=1,3
!    avg = SUM(forces(:,i_dim,1))/REAL(cell%numIon)
!    forces(:,i_dim,1) = forces(:,i_dim,1) - avg
!    WRITE(outputUnit,*) "average force along dimension ", i_dim, " is ", avg
!  END DO

  RETURN

END SUBROUTINE ForcesII_IE

END MODULE CalForces

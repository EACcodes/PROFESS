MODULE KEDF_TF 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_TF 
!     |_SUBROUTINE CalTF
!     |_FUNCTION TFEnergy
!     |_FUNCTION TFPotentialSqrt
!     |_FUNCTION TFStress
!
! DESCRIPTION:
!   Calculates the Thomas-Fermi kinetic forces, stress, and energies.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/15/2003  File Created by consolidating existing energy files (GSH) 
!   12/18/2003  Changed SUBROUTINE to FUNCTIONS
!   01/07/2004  Changed qTable and qMask and most instances of rhoR to 
!               spin-independent forms.
!   02/05/2004  Added nine stress functions. (VLL)
!   02/25/2004  Reorganized into different files (GSH)
!   03/23/2005  Added energy and potential for hybrid TF (VLL)
!   11/23/2013  Mohan split TF from FunctionalKinetic
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!

  USE Constants, ONLY: DP       ! shortcut for Double Precision
  USE PlaneWave, ONLY: qTable   ! Table of the norms of the q-vectors.
  USE PlaneWave, ONLY: qVectors ! Table of actual q vectors.
  USE PlaneWave, ONLY: qMask    ! Array of boolean that marks the recip. space 

  IMPLICIT NONE

  REAL(KIND=DP) :: lambda = -100.0_DP 
  ! Multiplier for the TF contribution.
  !
  REAL(KIND=DP), PARAMETER :: cTF = 2.87123400018819_DP               
  ! Coefficients for Thosmas-Fermi term
  !

CONTAINS


SUBROUTINE CalTF(potential, rho, calcEnergy, eTF)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the Thomas-Fermi potential and energy.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell, numSpin

  IMPLICIT NONE
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(OUT) :: potential
  ! TF potential, last dimension is spin
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho
  ! charge density on real space
  ! 
  LOGICAL, INTENT(IN) :: calcEnergy
  ! calculate TF energy or not
  ! 
  REAL(KIND=DP), INTENT(OUT) :: eTF
  ! Thomas-Fermi term energy
  !
  !! >> LOCAL VARIABLES << !!
  !
  INTEGER :: isp
  ! spin index
  !

  SELECT CASE(numSpin)
    CASE(1)
      ! Thomas-Fermi (TF)
      potential = potential + SPREAD(TFPotential(rho(:,:,:,1)), 4, numSpin)
      IF (calcEnergy) THEN
        eTF = TFEnergy(rho(:,:,:,1))
      ENDIF

    CASE(2)
      
      eTF = 0.0_DP
      ! Thomas-Fermi (TF)
      DO isp = 1, 2
        potential(:,:,:,isp) = potential(:,:,:,isp) + TFPotential(2.d0*rho(:,:,:,isp))
        IF (calcEnergy) THEN
          eTF = eTF + 0.5d0 * (TFEnergy(2.d0*rho(:,:,:,isp)))
        ENDIF
      ENDDO

    END SELECT
        
  RETURN

END SUBROUTINE CalTF


FUNCTION TFEnergy(rhoR_SI, weightingfunc, IfCalForce)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the Thomas-Fermi kinetic energy based on the real-
!   space electron density.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/14/2003  File created.  (VLL)
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP) :: TFEnergy             
  ! The TF energy. 
  !
  REAL(KIND=DP),DIMENSION(:,:,:),INTENT(IN) ::rhoR_SI 
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=DP),DIMENSION(:,:,:),INTENT(IN),OPTIONAL::weightingfunc
  !
  LOGICAL, INTENT(IN),OPTIONAL:: IfCalForce
  
                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP),PARAMETER :: ft = 5.0_DP / 3.0_DP 
  ! five thirds
  !

                           !>> FUNCTION BODY <<!

 IF(.NOT.PRESENT(IfCalForce).OR..NOT.IfCalForce) THEN 
    TFEnergy = lambda * cTF * SUM(rhoR_SI**ft) !*weightingfunc**frord) 
 ELSE 
    TFEnergy = lambda * cTF * SUM(rhoR_SI**ft*weightingfunc)
 END IF
 
 RETURN
 
END FUNCTION TFEnergy 


FUNCTION TFPointEnergy(rhoR_SI, i, j, k)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the Thomas-Fermi kinetic energy at a single grid
!   point.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/10/2004 (GSH)
!------------------------------------------------------------------------------
  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  !
  REAL(kind=DP) :: TFPointEnergy
  ! The TF energy at the gridpoint
  !
  INTEGER, INTENT(IN) :: i, j, k
  ! The location of hte gridpoint in the grid
  !
  REAL(kind=DP), DIMENSION(0:,0:,0:), INTENT(IN) :: rhoR_SI
  ! Electron density in real space (spin-independent)
  !
  ! >> INTERNAL VARIABLES <<!
  !
  REAL(kind=DP), PARAMETER :: ft = 5.0_DP / 3.0_DP 
  ! five thirds
  !

                           !>> INITIALIZATION <<!
                           !>> FUNCTION BODY <<!
  TFPointEnergy = lambda * cTF * rhoR_SI(i,j,k)**ft

  RETURN

END FUNCTION TFPointEnergy


FUNCTION TFPotential(rhoR_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the Thomas-Fermi potential in real-space
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   1/11/2003  Function created.  (GSH)
!------------------------------------------------------------------------------
  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoR_SI
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(SIZE(rhoR_SI, 1), SIZE(rhoR_SI, 2), SIZE(rhoR_SI, 3)) :: TFPotential
  ! The TF potential.
  !
  !>> INTERNAL VARIABLES <<!
  !
  REAL(KIND=DP), PARAMETER :: tt = 2.0_DP / 3.0_DP 
  ! two thirds
  !

                           !>> INITIALIZATION <<!
                           !>> FUNCTION BODY <<!

  TFPotential = (5._DP / 3._DP * cTF * lambda) * rhoR_SI**tt

  RETURN

END FUNCTION TFPotential


FUNCTION TFPotentialSqrt(rhoR_SI, weightingfunc, IfCalForce)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the Thomas-Fermi potential in real-space
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   1/11/2003  Function created.  (GSH)
!------------------------------------------------------------------------------
  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoR_SI 
  ! Electron density in real space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN),OPTIONAL:: weightingfunc
  !
  LOGICAL, INTENT(IN),OPTIONAL:: IfCalForce
  !
  REAL(KIND=DP), DIMENSION(SIZE(rhoR_SI, 1), SIZE(rhoR_SI, 2), SIZE(rhoR_SI, 3)) :: TFPotentialSqrt 
  ! The TF potential wrt the sqrt of the density
  !
  !>> INTERNAL VARIABLES <<!
  !
  REAL(KIND=DP), PARAMETER :: ss = 7.0_DP / 6.0_DP

                           !>> INITIALIZATION <<!
                           !>> FUNCTION BODY <<! 
  IF(.NOT.present(IfCalForce).OR..NOT.IfCalForce) THEN
    TFPotentialSqrt = (10._DP / 3._DP * cTF * lambda)*rhoR_SI**ss
  ELSE 
    TFPotentialSqrt = weightingfunc*(10._DP / 3._DP * cTF * lambda)*rhoR_SI**ss 
  END IF

  RETURN

END FUNCTION TFPotentialSqrt


FUNCTION TFStress(cellVol, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the Thomas-Fermi stress on local process
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/05/2004  Function created.  (Vincent Ligneres)
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: n3G, m3G

  IMPLICIT NONE
  !
  !>> EXTERNAL VARIABLES <<!
  !
  REAL(KIND=DP), INTENT(IN) :: cellVol  
  ! The volume of the cell
  !
  REAL(KIND=DP), INTENT(IN) :: energy   
  ! The TF energy
  !
  REAL(KIND=DP), DIMENSION(3,3) :: TFStress    
  ! The 3x3 tensor
  !
  !>> INTERNAL VARIABLES <<!
  !
  INTEGER :: a  
  ! counter
  !
  REAL(KIND=DP) :: temp        
  ! temporary storage
  !

  TFStress = 0._DP
  temp = -2._DP * energy / (3._DP * cellVol)
#ifdef __USE_PARALLEL
  temp = temp * REAL(n3G,KIND=DP)/REAL(m3G,KIND=DP)
#endif

  ! Put temp along the diagonal.
  DO a = 1, 3
    TFStress(a,a) = temp
  END DO

  RETURN

END FUNCTION TFStress


END MODULE KEDF_TF

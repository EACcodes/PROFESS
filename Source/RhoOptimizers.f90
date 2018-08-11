MODULE RhoOptimizers
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE RhoOptimizers
!     |_SUBROUTINE NoMinimization
!     |_SUBROUTINE LineSearch
!     |_SUBROUTINE RescaleRho
!
! DESCRIPTION:
!   This module contains the various electron minimization schemes we've come 
!   up with so far. 
!
!   It transforms the electron density (rho) into the lowest energy rho it can
!   find, hopefully the global energy minimum, at least a local min. 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! Some of this is heuristic and can be improved.  One way to do so for a 
! particular type of system is to make use of the GRAPHLINE function to plot
! various line searches.
!
! REFERENCES:
! 1.  Press, W.H., et. al.  Numerical Recipies: The Art of Scientific
!     Computing.  Cambridge University Press, Cambridge, 1986.
! 2.  Gill, P.E., et.al.  Practical Optimization.  Academic Press, San Diego, 
!     1986 (12th printing in 2000).
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/16/2004  File created.  (VLL)
!   02/23/2004  Added CGMIN, SDMIN (influenced by VLL's original codes).
!     Optimized CGMIN, added restarts and preconditioners.
!
!------------------------------------------------------------------------------
                              !<< GLOBAL >>
  USE CONSTANTS, ONLY : DP ! Shortcut for Double Precision
  USE CalPotPlus, ONLY: CalculatePotentialPlus
  USE MPI_Functions

  IMPLICIT NONE

  ! Variables potentially set by the initializer
  INTEGER :: &
  rhoOutcome = -1, &     ! Success of rho optimizer
                         ! 0 if outcome successful
                         ! 1 if exceeded max # iterations without finding min
                         ! 2 if other failure (e.g wrong line search direction)
  maxIter = 500,  &    ! Maximum steps take in rho minimizer
  lineMaxIter = 5, &     ! The maximum iterations for the brent line minimizer
  lineRunType = 1        ! 1 for parabolic interpolation, 0 for no parabolic
                         ! interpolation

  REAL(kind=DP) :: &
  pot_tol = 1.E-5, &      ! Tolerance for the norm of the potential
                          ! 5e-5 is too high to see dislocation dissociation. 
  tole    = 1.E-6, &         ! Default is 1E-5 eV
  lineTolPercent = .025_DP, &! The tolerance for the Brent line minizer
  fracMinPot = 1.5_DP      ! In Greg's alternative minimizer, the ratio of 
                           ! conjugate gradient of the current step to the last
                           ! step that will cause restart in steepest descent.

  LOGICAL :: calcReal = .FALSE.     
  ! Calculating all terms in real space?
  !
  LOGICAL :: cheatPot = .FALSE.     
  ! Initially skip WT/WGC terms to save time?
  !
  LOGICAL :: usePreconditioner = .FALSE.
  !
  INTEGER :: numEnergyLineMin
  ! The # of energy calculations it took to find the 
  ! minima during the linesearch
  !
  INTEGER :: numEnergyBracket
  ! The # of energy calculations it took to initially bracket the minima.
  !
  INTEGER :: numEnergyTotal
  ! The total number of energy calculations performed
  ! Throughout the minimization.
  !
  INTEGER :: numPotentialTotal  
  ! Total number of potential calculations performed
  !
  LOGICAL :: fixedmag = .FALSE.
  ! .true. for constrained spin calculations
  !
  CHARACTER(len=500) :: conv_check = "ENE"

  LOGICAL (KIND=DP) :: outDen 
  !output density or not, mohan add 10-12-12,
  !if this set to true, will reduce the efficiency of parallel
  !
  INTEGER :: niter_extra = 50   
  ! mohan add, maximal electronic steps and for extra electronic iteration.
  !
  INTEGER :: nmag = 100 
  ! max magnetic iteration
  !
  REAL(KIND=DP) :: tolm = 5.0e-6 
  ! default stop criteria for magnetic moment.  
  !

CONTAINS


SUBROUTINE NoMinimization(energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine computes the electronic energy associated with the current
!   density but does no minimization whatsoever. It is intended mostly for
!   testing purposes.
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
!   01/22/2004  Subroutine created.  (Greg Ho)
!
!------------------------------------------------------------------------------

  USE Sys, ONLY: rhoR
  USE CellInfo, ONLY: n1G, n2G, n3G, numSpin

  IMPLICIT NONE

  REAL(KIND=DP), DIMENSION(:), INTENT(OUT) :: energy
  ! where the final energies are stored   
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G, numSpin) :: tmpPot

  CALL CalculatePotentialPlus(rhoR, .TRUE., tmpPot, energy)

  RETURN

END SUBROUTINE NoMinimization


END MODULE RhoOptimizers

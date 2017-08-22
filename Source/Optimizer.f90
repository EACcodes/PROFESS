MODULE Optimizer
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Optimizer
!     |_OptimizeCell
!     |_OptimizeIons
!     |_OptimizeRho
!
! DESCRIPTION:
!   This module contains procedures that choose which optimizer to use from
!   the various optimization libraries available to optimize the cell, ions, 
!   and density in the OFDFT program. It has knowledge of the entire system.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
!------------------------------------------------------------------------------
! REFERENCES: 
!
!------------------------------------------------------------------------------

         !! >> GLOBAL << !!

  USE CONSTANTS, ONLY: DP        ! A double precision number

  USE OUTPUT, ONLY : errorUnit   ! The error file
  USE OUTPUT, ONLY : WrtOut

  USE SYS, ONLY: rhoR      ! The electron density in real space
  USE SYS, ONLY: frozenIon ! Logical array for whether ions are fixed (no optimization)
  USE SYS, ONLY: potReal   ! The potential in real-space
  USE SYS, ONLY: forceIon  ! The forces
  USE SYS, ONLY: stress    ! The stress
  USE SYS, ONLY: energy
  USE SYS, ONLY: grids     ! Multigrid hierarchy

  IMPLICIT NONE

                                 !>> GLOBAL <<
           ! ***** PARAMETERS FOR WHAT TO DO OVERALL *******
  LOGICAL, DIMENSION(5) :: calOption = & 
              (/.TRUE., &   ! Table of calculation options. The elements 
               .FALSE., &  ! correspond to whether following should be done:
               .FALSE., &  !  (1) Minimize the electron density. (default on)
               .FALSE., &  !  (2) Minimize ionic positions inside unit cell. (default off)
               .FALSE. /)  !  (3) Relax the cell, minimize stress tensor. (default off)
                           !  (4) Calculate the forces on the ions. (default off)
                           !  (5) Calculate the stress tensor. (default off)

  ! ***** PARAMETERS FOR CELL LATTICE MINIMIZATION *******
  INTEGER :: cellRelax = -1 ! Cell lattice minimization dimensions
                            ! Default: no minimization or all directions
                            ! (1) Optimize only in the x-direction
                            ! (2) Optimize only in the y-direction
                            ! (3) Optimize only in the z-direction

  ! ***** PARAMETERS FOR ION POSITION MINIMIZATION *******
  INTEGER :: ionMethod = -1 ! algorithm used for ion position minimization
                            !  (0) No minimization. (default) 
                            !  (2) Quickmin
                            !  (3) Legacy CG
                            !  (4) Conjugate gradient 2
                            !  (5) BFGS

  ! ***** PARAMETERS FOR ELECTRONIC DENSITY MINIMIZATION *******
  INTEGER :: rhoMethod = 1 ! algorithm used for electron density minimization:
                           !  (0) No minimization 
                           !  (1) Sqrt Newton minimization conserving #
                           !  (2) CG conserving #
                           !  (3) BFGS conserving #
                           !  (4) normal Sqrt Newton minimization
                           !  (5) normal Conjugate Gradient Minimization
                           !  (6) Conjugate gradient/Newton hybrid
                           !  (7) Log Truncated Newton Minimization

  !***** PARAMETERS FOR MD simulation *******
  INTEGER :: mdType = -1   ! -1:  no MD
                           !  1:  NVT, nose - hoover
                           !  2:  NPT full cell fluctuation, MD
                           !  3:  NVE

CONTAINS


SUBROUTINE DoMolecularDynamics
!------------------------------------------------------------------------------
! DESCRIPTION:
!  This subroutine does MD simulation. It passes the handler for density optimizer
!  to the MD subroutine.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE MOLECULARDYNAMICS, ONLY : md_output_path
  USE NVT, ONLY : RunNVT
  USE NPT, ONLY : NPTFullCell
  USE NVE, ONLY : RunNVE
  USE Output, ONLY : outputUnit
  USE RhoOptN, ONLY : infoUnit
  USE MPI_Functions, ONLY : rankGlobal

  IMPLICIT NONE

  CHARACTER(LEN=100) :: command

  ! >> INITIALIZE << !
  CALL StartClock('MD')

  ! >> FUNCTION << !

  ! Choose which type of output of information
  ! If infoUnit = 6, output them on the screen,
  ! otherwise infoUnit = outputUnit, output them into .out file.
  infoUnit = outputUnit

  ! found the directory that contains mdfiles.
  command='test -e ' // TRIM(md_output_path) // ' || mkdir ' // TRIM(md_output_path)
  IF(rankGlobal==0) CALL SYSTEM(command)
 
  SELECT CASE (mdType)
    CASE(1) 
    ! MD with Nose-Hoover thermostat (NVT)
    CALL RunNVT(OptimizeRho, energy, forceIon)
    CASE(2) 
    ! MD with NPT
    CALL NPTFullCell(OptimizeRho)
    CASE(3)
    ! MD with NVE
    CALL RunNVE(OptimizeRho, energy, forceIon)
    CASE (-1)
    PRINT *, "Optimizer : md = -1, invalid value for MD, error, code stop"
    STOP
  CASE DEFAULT
    PRINT *, "Optimizer : invalid value for key md found, stop!"
    PRINT *, "Optimizer : md: ", mdType
    STOP
  END SELECT

  CALL StopClock('MD')

  RETURN

END SUBROUTINE DoMolecularDynamics


SUBROUTINE OptimizeCell
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This the main routine of the program, the one that runs all the actual
!   computations. At the point this routine starts, we are done reading all the
!   information we need about the system and the job to perform, and on our way
!   to actually doing it. This routine relaxes the cell vectors so stress is
!   minimized. For this, it requires at each step that the ion positions be
!   optimized and that the stress be calculated. If we do not wish to change 
!   cell lattice vectors, the ion optimization is still run once.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE Output, ONLY: outputFinalGeometry
  ! print out the geometry at the end.
  !
  USE CellOptimizers, ONLY: NoOptimization 
  ! No Optimization
  !
  USE CellOptimizers, ONLY: SteepestDecent 
  ! Optimization by steepest decent
  !

  IMPLICIT NONE

  !>> FUNCTION BODY <<!

  ! 3 means do cell relaxation
  IF(calOption(3)) THEN
  
    CALL SteepestDecent(OptimizeIons,rhoR,energy,stress,cellRelax)
  ELSE
    ! 5 means only calculate the stress
    IF(calOption(5)) THEN
      CALL NoOptimization(OptimizeIons, rhoR, energy, stress, .TRUE.)
    ! else we don't do cell relaxtion and we don't calculate the stress
    ELSE
      CALL NoOptimization(OptimizeIons, rhoR, energy, stress, .FALSE.)
    END IF
  END IF

  RETURN

END SUBROUTINE OptimizeCell

 
SUBROUTINE OptimizeIons
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function will compute the real-space self-term contribution to the
!   Ewald sum (Eqn. 5 of [2]).  Uses units, of course, from [3].  
!   It does virtually nothing else.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE IonOptimizers, ONLY : NoOptimization  ! No optimization is done
  USE IonOptQui, ONLY: QuickMinOptimization ! Implements quickmin
  USE IonOptCG, ONLY: GradientOptimization  ! Greg's CG optimization
  USE IonOptCG2, ONLY: ConjugateGradient    ! Revised conjugate gradient (Chen Huang) (CG2)
  USE IonOptBFGS, ONLY: IonLBFGS            ! limited-memory BFGS method

  IMPLICIT NONE

  !>> FUNCTION BODY <<!

  ! Pick which minimization routine to run
  IF(calOption(2)) THEN

    IF (ionMethod == -1) THEN
      CALL WrtOut(6,'(OptimizeIons): you have not specified the method for ion relaxation, code stop')
      STOP
    ENDIF

    SELECT CASE(ionMethod)

      CASE(0) ! no minimization
        IF(calOption(4)) THEN
          CALL NoOptimization(OptimizeRho, rhoR, forceIon, .TRUE.)
        ELSE
          CALL NoOptimization(OptimizeRho, rhoR, forceIon, .FALSE.)
        ENDIF

      CASE(2) ! Velocity Verlet (quick minimization)
        CALL QuickMinOptimization(OptimizeRho, rhoR, energy, forceIon, frozenIon)

      CASE(3) ! Conjugate gradient (old)
        CALL GradientOptimization(OptimizeRho, rhoR, energy, forceIon, frozenIon)

      CASE(4) ! Better conjugate gardient (CG2)
        CALL ConjugateGradient(OptimizeRho, rhoR, energy, forceIon, frozenIon)

      CASE(5) ! BFGS (a qausi-newton method)
        CALL IonLBFGS(OptimizeRho, rhoR, energy, forceIon, frozenIon)

      ! This should normally not occur.
      CASE DEFAULT   
        WRITE(errorUnit,*) &
          '(OptimizeIons): The ion method selected was not implemented. STOP!'
        STOP

      END SELECT
  ELSE
    IF(calOption(4)) THEN
      CALL NoOptimization(OptimizeRho, rhoR, forceIon, .TRUE.)
    ELSE
      CALL NoOptimization(OptimizeRho, rhoR, forceIon, .FALSE.)
    END IF
  END IF

END SUBROUTINE OptimizeIons


SUBROUTINE OptimizeRho
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine minimizes the total energy with respect to electron density
!   according to the chosen method. It computes the initial energy, selects the
!   appropriate minimization algorithm and runs it.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
  USE RhoOptimizers, ONLY: NoMinimization ! Method #0: No minimization
  USE RhoOptimizers, ONLY: maxIter
  USE RhoOptSCG, ONLY: SqrtGradientMinimization    ! Method #1, part of 4: Conjugate gradient / steepest descent
  USE RhoOptN,  ONLY: ProjNormRhoMinimization     ! Method #2, 5, 6: TN with electron number conserved
  USE RhoOptSTN, ONLY: SqrtNewtonMinimization      ! Method 3, part of 4: NEW
  USE RhoOptLOG, ONLY: LogTruncatedNewton          ! Method 7

  IMPLICIT NONE

                     !>> EXTERNAL VARIABLES <<!

                     !>> INTERNAL VARIABLES <<!

  INTEGER :: dir_type 
  ! type of directions for Truncated-Newton
  !
  LOGICAL :: cgmin
  ! whether we use cg method in SqrtGradientMinmization
  !
  INTEGER :: firstIterMax
  ! For hybrid method, the number of iterations for first
  ! method.
  !
                       !>> INITIALIZATION <<!

                        !>> FUNCTION BODY <<!

  IF(calOption(1)) THEN

    if (rhoMethod==-1) then 
      call WrtOut(6,'(OptimizeRho): you have not specified the method for density relaxation, code stop')
      stop
    endif

    ! Select the minimization method for the density.
    SELECT CASE(rhoMethod) 

      CASE(0)
        CALL NoMinimization(energy)

      !-------------------------------------------------------------------

      CASE(1) ! NTN
        ! This truncated newton method conserve electron # automatically,
        ! not by scaling which is case for method = 4 below
        IF(ALLOCATED(potReal)) DEALLOCATE(potReal)
        dir_type = 1
        CALL ProjNormRhoMinimization(energy, dir_type)

      CASE(2) ! NCG
        ! This CG method conserve electron # automatically,
        IF(ALLOCATED(potReal)) DEALLOCATE(potReal)
        dir_type = 2
        CALL ProjNormRhoMinimization(energy, dir_type)

      CASE(3) ! NBF 
        dir_type = 3
        CALL ProjNormRhoMinimization(energy, dir_type)

      !-------------------------------------------------------------------

      CASE(4) ! STN
        ! This is like, so incredibly ugly.  But saves memory.
        IF(ALLOCATED(potReal)) DEALLOCATE(potReal)
        CALL SqrtNewtonMinimization(energy)

      CASE(5) ! SCG
        ! This is like, so incredibly ugly.  But saves memory.
        IF(ALLOCATED(potReal)) DEALLOCATE(potReal)
        cgmin = .TRUE.
        CALL SqrtGradientMinimization(grids(0), energy, cgmin, maxIter) 

      CASE(6) ! SCG+STN 
        IF(ALLOCATED(potReal)) DEALLOCATE(potReal)
        firstIterMax = 15
        cgmin = .TRUE.
        CALL SqrtGradientMinimization(grids(0), energy, cgmin, firstIterMax)
        CALL SqrtNewtonMinimization(energy)
        
      CASE(7) ! LOG
        CALL LogTruncatedNewton(rhoR, energy)

      CASE(8)
      !  CALL RhoTNWithLogLine(rhoR, energy)

      ! This should normally not occur.
      CASE DEFAULT   
        WRITE(errorUnit,*) 'The method selected was not implemented. Leaving.'
        STOP

    END SELECT
  ELSE
    CALL NoMinimization(energy)
  END IF

END SUBROUTINE OptimizeRho


END MODULE Optimizer

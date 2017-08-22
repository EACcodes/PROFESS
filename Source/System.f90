MODULE SYS
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE System
!     |_SUBROUTINE SetupSystem
!     |_SUBROUTINE CleanSystem
!     |_SUBROUTINE GetMachinePrecision
!     |_SUBROUTINE CalculateElectronNumber
!
! DESCRIPTION:
!   This is the module that describes all the important parameters of our
!   system at any given time, such as the density, the potential, the ion
!   positions, which gridpoints constitute the boundaries of the cell, the 
!   stress, etc.  This information is limited to be known only by a SELECT
!   FEW modules (not ALL).  As of this writing, the priviliged modules are
!   InitialzeInputs, Optimizer, and Output.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! 1) xmin, xmax, ymin, ymax, zmin, zmax for rhoR should be updated,
! it's very dangerous to use xmin=0, ymin=0. zmin=0. Because usually
! Fortran starts from 1.
!
! 2) to replace these variables:
! G space FFT dimensions:
! totX/2+1 -> k1G
! totZ -> k2G
! locY -> k3G
! locYOff -> k3Goff
! Real space FFT dimensions in local:
! totX -> n1G
! totY -> n2G
! locZ -> n3G
! Real space FFT dimension in global:
! totX -> m1G
! totY -> m2G
! totZ -> m3G
! locZoff -> n3Goff
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2004-02-25 Pieced together for organizational purposes. (GSH)
! 2004-03-10 Deleted numEle (the number of electrons) (GSH)
! 2013-12-06 Mohan update.
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY : DP                  ! Double Precision
  USE CONSTANTS, ONLY : machprec
  USE MPI_Functions

                              !>>GLOBAL<<!
  IMPLICIT NONE

#ifdef __USE_PARALLEL
! When running in parallel, arrays used in electronic calculations use
! domain decomposition.  Arrays in real space are sliced in the z-direction,
! with dimensions (totX,totY,locZ) with offset locZOff.  Arrays in reciprocal
! space have are sliced in the same way.
  INTEGER :: numIonLoc  ! Number of ions on local processor (used in Ewald)
  INTEGER :: numIonInit ! Ion number offset (used in Ewald)
#endif

  REAL(KIND=DP), DIMENSION(9) :: energy = 0.0_DP        
  ! Table of energies. From left to right: 1) total energy,
  ! 2) kinetic, 3) external, 4) coulombic, 5) exchange-correlation,
  ! 6) ion-ion, 7) Thomas-Fermi, 8) von Weiszacker and 
  ! 9) third term Wang-Teter, WGC, ...)

  REAL(KIND=DP) :: magmom = 0.0_DP 
  ! Magnetic moment, for initialization of density
  !
 
  !!! RELATED TO ELECTRONIC OPTIMIZATION !!!
  
  REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: rhoR 
  ! Electron Density in real space.
  !
  LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: interior 
  ! The interior of the system (not boundary)
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: potReal 
  ! The real-space potential
  !

  !!! RELATED TO IONIC/CELL OPTIMIZATION !!!

  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: frozenIon
  ! Whether each ion is fixed in a given dimension
  ! 1st index is ion, second is direction x, y, or z.
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: forceIon
  ! Forces on ions.  First index is the ion and second
  ! index is the direction x, y, or z.  Final index
  ! is 1 if we're talking about the total force, 
  ! 2 for ion-ion and 3 for ion-electron.
  !
  REAL(KIND=DP), DIMENSION(3,3) ::  stress  
  ! The stress

  !------------------------------------------
  ! (1) A Type that combines all the grids 
  ! necessary to describe the system.
  !------------------------------------------
  TYPE :: gridPack

    ! for AMD
    REAL(KIND=DP), DIMENSION(:,:,:), POINTER :: rhoCore
    REAL(KIND=DP), DIMENSION(3,3) :: &
    section              ! The coarse grid can be split up into three
                         ! sections along the cell vectors.  The first section
                         ! will be the section before the ULTIMATE fine grid.
                         ! The second section is the section of the fine grid,
                         ! and the third is the section after the fine grid.
                         ! Here is the fraction that each section takes up of
                         ! the whole grid, in each of the 3 directions.

    ! end for amd

  END TYPE gridPack

  TYPE(gridPack), DIMENSION(:), ALLOCATABLE :: grids

  REAL(KIND=DP) :: gridSpacing = -1.0_DP  
  ! Density in 1/A of grid points (in 1 dimension)
  !
  LOGICAL :: bvac =.FALSE. 
  ! no vacuum regulation (using cutoff) in WGC, WT, and CAT
  !
  REAL(KIND=DP) :: rho0 = -1.0_DP     
  ! Average or user-specified value of the density for WT and WGC.
  !
  REAL(KIND=DP) :: rhoS = -1.0_DP     
  ! Value around which the kernel is Taylor-expanded in WGC (rho*)
  !
  LOGICAL :: hold0=.FALSE. 
  ! True means rho0 does not change when the cell varies.
  !
  LOGICAL :: holdS=.FALSE. 
  ! False means rhoS is always numele/Volume(cellReal).
  !
  
CONTAINS



SUBROUTINE SetupSystem(energyCutoff, kinetic)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine takes the cell from ReadGeometry to compute the table of 
!   q-vectors. It then calls the subroutine FillQTable, that is separate
!   because it needs to be called every time we alter the cell shape. In event
!   that the WT or WGC kinetic energy functionals are used, it will also
!   allocate and fill the kernel table. Once again, the allocation and filling
!   parts are kept separate to allow updating of the kernel when the cell size
!   is altered.
!
! CONDITIONS AND ASSUMPTIONS:
!   It is assumed that the real-space grid is odd along each of its 
!   dimensions. If you want to change this, you need to read the warning note
!   in the initialization part or expose yourself to considerable trouble.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   1. Ashcroft and Mermin, Solid State Physics.  Harcourt College Publishers,
!      Fort Worth, 1976.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/30/2003  File created.  (Vincent Ligneres)
!   11/19/2003  Checking validity on Al_fcc (4.03A) Grid size does not match
!               that of Stuart's code. I get 16 he has 24 for 1200eV cut. (VLL)
!   11/20/2003  Corrected formulas numX, Y, Z. Tested on Al_fcc, it works. &
!               (VLL)
!   11/22/2003  Cosmetic Changes (Greg Ho)
!   12/05/2003  Added the kernel table and the FillKernel call. (VLL)
!   12/12/2003  Removed call to FillQTable, actually removed the whole 
!               procedure.  Replaced with a few lines of code.
!   01/08/2003  Renamed this procedure "SetupGridObjects"
!   08/23/2013  Update this subroutine by MOHAN CHEN 
!
!------------------------------------------------------------------------------
  USE OutputFiles, ONLY : outputUnit ! Add by Mohan
  USE CellInfo, ONLY : m1G, m2G, m3G, n1G, n2G, n3G, n3Goff ! global FFT dimensions, add by Mohan
  USE CellInfo, ONLY : cell, numSpin
  USE SetupFFT, ONLY: InitializeFFT

  IMPLICIT NONE

                     !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(IN) :: kinetic      ! Parameter depending on choice of KEDF
  REAL(KIND=DP), INTENT(IN) :: energyCutoff ! The kinetic energy cutoff

                     !>> INTERNAL VARIABLES <<! 

#ifdef __USE_PARALLEL
  INTEGER :: numExtra    ! MOD(numIons,sizeGFFT)
  INTEGER :: numIons
#endif
  INTEGER :: xmin, xmax, ymin, ymax, zmin, zmax ! index for density
  INTEGER :: allocateStatus           ! Checks that tables are allocated alright.
  INTEGER :: isp                      ! spin index
  INTEGER :: kePotDim                 ! dimensions of kePot
  REAL(KIND=DP) :: startElectronNumber = 0.D0
  REAL(KIND=DP) :: sumrho = 0.0_DP

  ! >>>>>>>>>>>> FUNCTION BEGINS <<<<<<<<<<<<<<<<!

  CALL Title("System::SetupSystem")

  machprec = GetMachinePrecision()
  WRITE(outputUnit,*) '(System) Machine Precision              : ', machprec

  CALL InitializeFFT(energyCutoff, gridSpacing, m1G, m2G, m3G, n3G, n3Goff)

  n1G = m1G
  n2G = m2G

  !----------------------------------------------------------------------------------------------

#ifdef __USE_PARALLEL

  ! Setup some ion parallelization information
  numIons = SIZE(cell%ionTable)
  numIonLoc = numIons/sizeGFFT
  numExtra = MOD(numIons,sizeGFFT)
  IF (rankGFFT<numExtra) THEN
     numIonLoc = numIonLoc+1
     numIonInit = rankGFFT*numIonLoc + 1
  ELSE IF (rankGFFT >= numIons) THEN
     numIonLoc = 0
  ELSE
     numIonInit = numExtra + rankGFFT*numIonLoc + 1
  ENDIF
  
#endif

  ! for periodic boundary condition
  xmin = 0
  xmax = n1G-1
  ymin = 0
  ymax = n2G-1
  zmin = 0
  zmax = n3G-1
  
  ! For kePot, if we are using WT, we want kePot to have 2 dimensions
  ! in the last dimension. If WGC, then 3
  IF(kinetic == 4) THEN
    kePotDim = 2
  ELSE IF(kinetic == 5 .or. kinetic==10 .or. kinetic==11 ) THEN  ! WGC or CAT KEDF or Huang-Carter KEDF
    kePotDim = 3
  ELSE
    kePotDim = -1
  END IF

                           !>> FUNCTION BODY <<!

  ! Allocate memory for rhoR.
  ALLOCATE(rhoR(xmin:xmax, ymin:ymax, zmin:zmax, numSpin), stat=allocateStatus)

  WRITE(outputUnit, *) "Dimension of rho for this processor : "
  WRITE(outputUnit, *) "xmin, xmax : ",xmin, xmax
  WRITE(outputUnit, *) "ymin, ymax : ",ymin, ymax
  WRITE(outputUnit, *) "zmin, zmax : ",zmin, zmax

  IF (allocateStatus/=0) THEN
    WRITE(*,*)'Error allocating the real sp. density table. Leaving.'
    STOP
  END IF


  ! Allocate memory for the boundary

  ALLOCATE(interior(xmin:xmax, ymin:ymax, zmin:zmax, numSpin), &
       stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*)'Error allocating the interior table. Leaving.'
    STOP
  END IF

  ! Allocate memory for the potentials
  ALLOCATE(potReal(xmin:xmax, ymin:ymax, zmin:zmax, numSpin), &
           stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*)'Error allocating the real sp. potential table. Leaving.'
    STOP
  END IF

  ! Allocate the force table
  !ALLOCATE(forceIon(SIZE(cell%ionTable),3,3), stat=allocateStatus) 
  ! mohan changes the last dimension to 6, because 4,5,6 dimensions
  ! are required by AMD algorithms
  ALLOCATE(forceIon(SIZE(cell%ionTable),3,6), stat=allocateStatus) 
  IF (allocateStatus/=0) THEN                    
    WRITE(*,*)'Error allocating the ion forces table. Leaving.'
    STOP
  END IF

  ! Initialize the interior array


  ! Periodic boundary conditions.  All points valid
  interior(:,:,:,:) = .TRUE.

  ! Initialize RhoR to a uniform starting density
  startElectronNumber = NINT(SUM(cell%elementTable%chargeTot))

  IF( cell%vol .EQ. 0.0_DP ) THEN
    WRITE(outputUnit,*) "The cell volume is 0, check the lattice vectors."
    STOP
  ENDIF

  rhoR = REAL(startElectronNumber,KIND=DP) / cell%vol

  WRITE(outputUnit,*) "Total electron number = ",SUM(cell%elementTable%chargeTot)
  WRITE(outputUnit,*) "Uniform rhoR = ", rhoR(xmin,ymin,zmin,1)

  ! initialize of the density according to magnetic moment
  IF(numspin == 2) Then
     DO isp = 1,2
        rhoR(:,:,:,isp) = 0.5d0*(rhoR(:,:,:,isp) - ((-1)**isp) * magmom/ cell%vol )
        sumrho = SUM(rhoR(:,:,:,isp)) * cell%vol / m1G / m2G / m3G 
        CALL ReduceRealLevel1(sumrho)
        WRITE(outputUnit,*) "Electron number for spin ", isp, " is ", sumrho
     ENDDO
  ENDIF

  ! Assign the minimum parameters required for conjugate gradient.  This
  ! code is partially duplicated in SetupGridPack, which is inelegant.
  ALLOCATE(grids(0:0))

 RETURN

END SUBROUTINE SetupSystem


SUBROUTINE CleanSystem(ionKeepOpt)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine just does some cleanup like deallocates memory, etc.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2013-12-06 Mohan update
!
!------------------------------------------------------------------------------
  USE Fourier, ONLY : CleanFFT
  USE CellInfo, ONLY : cell

  IMPLICIT NONE

                     !>> EXTERNAL VARIABLES <<!

  LOGICAL, INTENT(IN), OPTIONAL :: ionKeepOpt

                     !>> INTERNAL VARIABLES <<! 
  INTEGER :: i  ! dummy counter
  INTEGER :: status
  LOGICAL :: ionKeep

                     !>> INITIALIZATION <<!
  
  ionKeep = .FALSE.
  IF(PRESENT(ionKeepOpt)) THEN
    IF(ionKeepOpt .eqv. .TRUE.) ionKeep=.TRUE.
  END IF
                       !>> FUNCTION BODY <<!

  ! Deallocate the density
  DEALLOCATE(rhoR)
  DEALLOCATE(interior)
  IF (ALLOCATED(potReal)) DEALLOCATE(potReal)
  DEALLOCATE(forceIon)

  ! Deallocate the pseudopotentials if requested (also by default)
  IF (.NOT. ionKeep) THEN
    DO i = 1, SIZE(cell%elementTable) - 1
      DEALLOCATE(cell%elementTable(i)%psp%potValues)
      DEALLOCATE(cell%elementTable(i)%psp%vqS)
      DEALLOCATE(cell%elementTable(i)%psp%potDD)
      DEALLOCATE(cell%elementTable(i)%psp%t)
      DEALLOCATE(cell%elementTable(i)%psp)
    END DO
  END IF
  DEALLOCATE(cell%ionTable)
  DEALLOCATE(cell%elementTable)

  ! Deallocate frozen ion array
  IF (ALLOCATED(frozenIon)) DEALLOCATE(frozenIon)


  DEALLOCATE(grids,STAT=status)
  IF(status/=0) THEN
    WRITE(outputUnit,*) "Problem deallocating grids in CleanSystem"
    STOP
  ENDIF
 
  CALL CleanFFT

  RETURN

END SUBROUTINE CleanSystem


FUNCTION GetMachinePrecision() RESULT(eps)
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

  REAL(KIND=DP) :: eps
  
  eps = 1._DP
  eps = nearest(eps,1._DP)-eps

  RETURN
 
END FUNCTION GetMachinePrecision


SUBROUTINE CalculateElectronNumber(cellLoc,rho,dimX,dimY,dimZ,dimTotal,totEleNum)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the electron number of the system, used in small box technique
! test.
! 
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! Move this function to somewhere else.
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2013-12-06 Mohan update.
!
!------------------------------------------------------------------------------
  
  USE CellInfo
  USE MathFunctions, ONLY : Volume

  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  TYPE(cellStruct), INTENT(IN) :: cellLoc
  INTEGER, INTENT(IN) :: dimX, dimY, dimZ ! dimension of rho
  INTEGER, INTENT(IN) :: dimTotal  ! total dimension of global FFT, used to 
                                   ! calculate dV
  REAL(KIND=DP), INTENT(INOUT) :: totEleNum ! total electron number
  REAL(KIND=DP), DIMENSION(dimX, dimY, dimZ), INTENT(IN) :: rho 
                                   ! density, may be from small box

  !! >> INTERNAL VARIABLES << !!
  REAL(KIND=DP) :: cellVol
  REAL(KIND=DP) :: dV

  !! >> FUNCTION
  cellVol = Volume(cellLoc%cellReal) ! get volume of the cell
  dV = cellVol / dimTotal         ! get delta volume
  totEleNum = SUM(rho(:,:,:))*dV  ! get total electron number

  !WRITE(outputUnit,*) "dV=",dV

  CALL ReduceRealLevel1(totEleNum)

  RETURN

END SUBROUTINE CalculateElectronNumber 


SUBROUTINE CleanSystemForESPDOpt
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine just does some cleanup like deallocates memory, etc.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!   Everything in this module.
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
  USE FOURIER, ONLY : CleanFFT
  USE CellInfo, ONLY : CellStruct

  IMPLICIT NONE

                     !>> INTERNAL VARIABLES <<!
  INTEGER :: i
  ! dummy counter
  !
  INTEGER :: status

                     !>> INITIALIZATION <<!

                       !>> FUNCTION BODY <<!
  ! Deallocate the density
  DEALLOCATE(rhoR)
  write(*,*) "INFO: rhocore, auxidm and fixingvariables must be deallocated, not complete implementation. will probably fail." ! jmd
  DEALLOCATE(interior)
  IF (ALLOCATED(frozenIon))DEALLOCATE(frozenIon)

  DEALLOCATE(grids,STAT=status)
  IF(status/=0) write(errorUnit,*) "Problem deallocating grids in CleanSystem"

  CALL CleanFFT

  RETURN

END SUBROUTINE CleanSystemForESPDOpt 


END MODULE SYS

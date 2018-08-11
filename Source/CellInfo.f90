MODULE CellInfo 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE CellInfo
!     |_SUBROUTINE SetupCellFFTdims
!     |_SUBROUTINE RefreshLattice
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
! 06/2013 Created by MOHAN CHEN
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP, systemNameLen
  USE OUTPUTFILES, ONLY: outputUnit

  IMPLICIT NONE

  !! >> REAL SPACE PARAMETERS < !!
  INTEGER :: n1G, n2G, n3G ! local fft dimension
  INTEGER :: n123G         ! n1G * n2G * n3G
  INTEGER :: m1G, m2G, m3G ! global fft dimension (the whole cell)
  INTEGER :: m123G         ! m1G * m2G * m3G
  INTEGER :: n3Goff        ! offset in this processor

  !! >> RECIRPOCAL SPACE PARAMETERS << !!
  INTEGER :: k1G, k2G, k3G ! G-space fft dimension
  INTEGER :: k123G         ! k1G * k2G * k3G
  INTEGER :: k3Goff        ! offset in this processor

  INTEGER :: numSpin = 1
  ! default value is 1, another option is 2
  !


  !------------------------------------------
  ! (1) A Type for a Pseudopotential
  !------------------------------------------
  TYPE :: pseudoPot
  !
  INTEGER :: numPoints
  ! numPoints: The total number of points we have in the pseudopotential.
  ! 
  INTEGER :: charge       
  ! charge: the charge with this pseudopotential
  !
  REAL(KIND=DP) :: maxG=0._DP 
  ! maxG: The maximum G that this pseudopotential goes to, in A^-1 (if recpot)
  ! 
  REAL(KIND=DP) :: maxR=0._DP 
  ! maxR: The maximum r that this pseudopotential goes to, in A (if realpot)
  !
  ! REAL(KIND=DP) :: r2=1.E-6_DP
  ! r2: (realpot only) Second r value in the pseudopot
  !
  REAL(KIND=DP), DIMENSION(:), POINTER :: potValues 
  ! potValues: An array of the values of the pseudopotential.
  ! In reciprocal space it is in evenly spaced increments, 
  ! running from G values of 0 to maxG, and in real space it is on a logarithmic grid.
  !
  REAL(KIND=DP), DIMENSION(:), POINTER :: potDD
  ! potDD : this array contains 2nd derivative at each q ponit
  ! This array is used by spline 
  !
  REAL(KIND=DP), DIMENSION(:), POINTER :: vqS 
  ! vqS : this array contains v(q) + 4*pi*Charge/q**2
  ! This array is used by spline 
  !
  REAL(KIND=DP), DIMENSION(:), POINTER :: t 
  ! t: this array contains knots for this pseudopotential, t(i)=i
  ! we use index of each q-point as knots in spline
  !
  CHARACTER(systemNameLen+7) :: pspFileName           
  ! pspFileName: The name of the pseudopotential file
  !
  INTEGER :: type=1 
  ! type: The type of pseudopotential.
  ! 1 = .recpot  2 = .realpot
  !
  END TYPE pseudoPot



  !------------------------------------------
  ! (2) A Type for element
  !------------------------------------------
  TYPE :: element
  CHARACTER(LEN=2) :: elementName 
  ! elementName: The name of the ion as listed in the ion file.     
  !
  REAL(KIND=DP) :: charge         
  ! chage: Z-core electrons.  The electronic charge per atom.
  !
  REAL(KIND=DP) :: chargeTot      
  ! chargeTot: Total electronic charge for this type of atom
  !
  REAL(KIND=DP) :: mass           
  ! mass: Mass of ion
  !
  TYPE(pseudoPot), POINTER :: psp 
  ! psp: The pseudopotential associated with the ion.
  !
  INTEGER :: firstIonID           
  ! firstIonID: index of the FIRST ion of this type in iontable
  !
  END TYPE element



  !------------------------------------------
  ! (3) A Type for an ion
  !------------------------------------------
  TYPE :: ion
  !
  INTEGER :: elementID
  ! elementID: Not chemical, but in the list of encountered types.
  !
  REAL(kind=DP), DIMENSION(3) :: coord=0._DP           
  ! coord: position of the ion In fractional coordinates.
  !
  END TYPE ion


  !------------------------------------------
  ! (4) Type for a cell (contains ions, etc)
  !------------------------------------------
  TYPE :: cellStruct
  !
  REAL(KIND=DP), DIMENSION(3,3) :: cellReal=0._DP 
  ! cellReal : Real-space unit-cell vectors (A, B, C) (Bohr)
  ! Note: read in cellReal is in the unit of Angstrom,
  ! but in that function cellReal is also transformed to Bohr unit.
  !
  REAL(KIND=DP), DIMENSION(3,3) :: cellRecip=0._DP
  ! cellRecip: Reciprocal-space unit cell vectors.
  ! cellRecip = TRANSPOSE(Inverse(cellReal)) * 2._DP * PI 
  !
  REAL(KIND=DP) :: lengthA=0._DP 
  REAL(KIND=DP) :: lengthB=0._DP
  REAL(KIND=DP) :: lengthC=0._DP
  ! cell Length for vector A, B, C (Bohr)
  !
  REAL(KIND=DP) :: angleAlpha=0._DP
  REAL(KIND=DP) :: angleBeta=0._DP
  REAL(KIND=DP) :: angleGamma=0._DP
  ! cell angleAlpha is the angle between A and B
  !
  REAL(KIND=DP) :: vol=0._DP
  ! cell volume (Bohr^3)
  !
  REAL(KIND=DP) :: dv=0._DP
  ! delta volume of real space grid
  !
  INTEGER :: numIonType = 0
  ! number of ion types
  ! 
  INTEGER :: numIon = 0
  ! number of total ions in the cell
  !
  REAL(KIND=DP) :: numEle = 0.0_DP
  ! number of total electron number
  !
  TYPE(ion), DIMENSION(:), POINTER :: ionTable 
  ! ionTable: Type + location of each ion, sorted by ion type
  ! We can get the number of ions (formerly numIon) by
  ! the size of this array.
  !
  TYPE(element), DIMENSION(:), POINTER :: elementTable
  ! elementTable: This array actually has one more element than 
  ! numIonType, because the last entry here numIon+1
  ! to make it easy to use this to loop through stuff
  !
  REAL(kind=DP) :: totNumEle
  
  END TYPE cellStruct



  ! Our cell, includes the cell lattice vectors and the ions
  TYPE(cellStruct), SAVE :: cell 

  LOGICAL :: usingGridPack = .FALSE.   ! Are we using grid packs for Multigrid?

  ! move somewhere else in the future
  LOGICAL :: storeRealIonPot = .TRUE. ! Store the real-space ionic potential (required
                                      ! for multigrid runs)

  INTEGER :: kinetic = 5 ! default is WGC KEDF


CONTAINS


SUBROUTINE SetupCellFFTdims
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
 

  USE FOURIER, ONLY : GetFFTDims        ! global FFT in real space
  USE FOURIER, ONLY : GetFFTComplexDims ! global FFT in G space

  IMPLICIT NONE

  CALL Title('CellInfo::SetupCellFFTdims')

  ! real space FFT dimensions
  CALL GetFFTDims(n1G, n2G, n3G, n3Goff)
  ! G space FFT dimensions
  CALL GetFFTComplexDims(k1G, k2G, k3G, k3Goff)

  ! m1G, m2G, m3G should be set in SetupSystem() 
  IF((m1G .LE. 0) .OR. (m2G .LE. 0) .OR. (m3G .LE. 0)) THEN
    WRITE(outputUnit,*) " check the value of m1G, m2G, m3G."
    STOP
  ENDIF
 
  n123G = n1G * n2G * n3G
  m123G = m1G * m2G * m3G
  k123G = k1G * k2G * k3G

  WRITE(outputUnit,'(A,3I8)') " (GB INFO) FFT size      : ", n1G, n2G, n3G 
  WRITE(outputUnit,'(A,I8)')  " (GB INFO) FFT offset    : ", n3Goff
  WRITE(outputUnit,'(A,3I8)') " (GB INFO) FFT G size    : ", k1G, k2G, k3G 
  WRITE(outputUnit,'(A,I8)')  " (GB INFO) FFT G offset  : ", k3Goff

  RETURN

END SUBROUTINE SetupCellFFTdims


SUBROUTINE RefreshLattice(newCell)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Refresh cell information, like lattices related variables. 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! The cell%dV is not set to a value if global grid size 'm123G' is 0,
! so, be very careful about this, don't use cell%dV until it is initialized
! we will optimize this in the future.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY : PI          ! 3.14159...
  USE MathFunctions, ONLY : Volume  ! calculate the volume of lattice
  USE MathFunctions, ONLY : Inverse ! inverse 3*3 matrix
  USE OutputFiles, ONLY : outputUnit
  USE MPI_Functions, ONLY: message, Error

  REAL(KIND=DP), DIMENSION(3,3) :: newCell ! Transformed cell including the rotation from stress. 

  !>> LOCAL VARIABLES <<!
  REAL(KIND=DP) :: a, b, c, ab, ac, bc

  CALL Title('CellInfo::RefreshLattice')
  
  ! real spac evectors
  cell%cellReal = newCell

  ! calculate the reciprocal space vectors
  cell%cellRecip = TRANSPOSE(Inverse(cell%cellReal)) * 2._DP * PI

 
  ! volume of real space cell
  cell%vol=Volume(cell%cellReal)

  IF(m123G > 0) THEN
    cell%dv=cell%vol/REAL(m123G, KIND=DP)
    WRITE(outputUnit,*) " dV of each real space grid is now ", cell%dV
  ENDIF

  ! calculate the length of the three new lattice vectors (unit is Bohr!)
  ! be careful of the correct reading order in fortran.
  ! mohan fix bug 2013-08-04
  a = SQRT(cell%cellReal(1,1)**2 + cell%cellReal(2,1)**2 + cell%cellReal(3,1)**2)
  b = SQRT(cell%cellReal(1,2)**2 + cell%cellReal(2,2)**2 + cell%cellReal(3,2)**2)
  c = SQRT(cell%cellReal(1,3)**2 + cell%cellReal(2,3)**2 + cell%cellReal(3,3)**2)

  ab = cell%cellReal(1,1)*cell%cellReal(1,2)+&
       cell%cellReal(2,1)*cell%cellReal(2,2)+&
       cell%cellReal(3,1)*cell%cellReal(3,2)
  bc = cell%cellReal(1,3)*cell%cellReal(1,2)+&
       cell%cellReal(2,3)*cell%cellReal(2,2)+&
       cell%cellReal(3,3)*cell%cellReal(3,2)
  ac = cell%cellReal(1,1)*cell%cellReal(1,3)+&
       cell%cellReal(2,1)*cell%cellReal(2,3)+&
       cell%cellReal(3,1)*cell%cellReal(3,3)

  ! calculate the length of the three new lattice vectors (unit is Bohr!)
  cell%lengthA = a
  cell%lengthB = b
  cell%lengthC = c

  ! calculate the angle between the new three lattice vectors.
  cell%angleAlpha = ACOS(bc/b/c)/PI*180.0_DP
  cell%angleBeta  = ACOS(ac/a/c)/PI*180.0_DP
  cell%angleGamma = ACOS(ab/a/b)/PI*180.0_DP

  ! print out the cell
  WRITE(outputUnit,'(A,3(ES12.3))') " Cell length (Bohr) ", cell%lengthA, cell%lengthB, cell%lengthC
  WRITE(outputUnit,'(A,3(ES12.3))') " Angle     (Degree) ", cell%angleAlpha, cell%angleBeta, cell%angleGamma
  WRITE(outputUnit,'(A,ES12.3)')    " Volume    (Bohr^3) ", cell%vol 

  ! Volume should be positive
  IF( cell%vol < 0._DP) THEN
    WRITE(*,*) ' cell volume < 0, check your lattice parameter'
    STOP
  ENDIF

END SUBROUTINE RefreshLattice


END MODULE CellInfo

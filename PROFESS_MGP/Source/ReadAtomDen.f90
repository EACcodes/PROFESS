MODULE AtomicDensity
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE AtomicDensity
!     |_SUBROUTINE GenerateAtomicDensity ! General routine.
!       |_SUBROUTINE ReadAtomicDensity ! read in atomic density in G space.
!     |_FUNCTION ADenLookup ! do interpolation to get 3D rho(G), private function
!     |_FUNCTION ADenRecip  ! get 3D rho(G), private function
!     |_SUBROUTINE ADenReal ! calculate atomic density rho(r), private function
!
! DESCRIPTION:
!   This modules read in atomic charge density and put it in the density
!   array.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/18/2012 Mohan created. 
!------------------------------------------------------------------------------

  !! << GLOBAL >> !!

  USE CONSTANTS, ONLY: DP, PI, bohr
  USE OutputFiles, ONLY: outputUnit
  USE PlaneWave, ONLY: qMask
  USE CellInfo, ONLY : cell
  USE CellInfo, ONLY: k1G, k2G, k3G

  IMPLICIT NONE

  TYPE :: atomDen 
  ! the type of atomic density
  !
    INTEGER :: nG 
    ! the total number of points we have in the .arho file.
    !
    REAL(KIND=DP) :: maxG 
    ! max G value in the file.
    !
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: rhoG 
    ! read in 1D rho in G space
    !
    REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: normG 
    ! norm of G vectors
    !
  END TYPE atomDen

  CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: atomicDensityFile 
  ! be allocated in ReadInputFiles.f90
  ! save the atomic density files be allocated in ReadInputFiles.f90
  !
  INTEGER :: atomTypeNum = 0 
  ! number of atom types
  !

CONTAINS


SUBROUTINE GenerateAtomicDensity(rho,dimX,dimY,dimZ,nspin)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Read in atomic density in real space, and transforms it into 1 dimensional
! G space array. Then generate the 3D recip array and transform it back
! to real space to get the final atomic charge density.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! This subroutine is called in function Initialize in InitializeInputs.f90
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012 Created by Mohan Chen.
!------------------------------------------------------------------------------

  USE MathFunctions, ONLY: Volume ! be removed! 

  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  INTEGER :: dimX, dimY, dimZ, nspin
  REAL(KIND=DP), DIMENSION(dimX,dimY,dimZ,nspin), INTENT(INOUT) :: rho
  TYPE(atomDen), DIMENSION(:), ALLOCATABLE :: aden
 
  !! >> INTERNAL VARIABLES << !!
  INTEGER :: is

  !! >> INITIALIZE << !!
  CALL Title("AtomicDensity::GenerateAtomicDensity")
  WRITE(outputUnit,'(A)') " Read in atomic charge density."

  ! how many types of atoms
  atomTypeNum = SIZE(atomicDensityFile)

  ! allocate atomic density type.
  ALLOCATE(aden(atomTypeNum))

  !! >> FUNCTION << !!

  ! read in atom density rho(G) in one dimension.
  CALL ReadAtomicDensity

  ! get the starting density as atomic density
  ! for each spin.
  DO is=1, nspin
    WRITE(outputUnit,*) "Generate atomic density for spin ",is
    CALL ADenReal(rho(:,:,:,is),aden)
    WRITE(outputUnit,*) "Minimal denisty value is ", MINVAL(rho(:,:,:,1))
    WRITE(outputUnit,*) "Maximal denisty value is ", MAXVAL(rho(:,:,:,1))
  END DO

  !------------------------------------------------------------
  ! calculate the number of points.
!  total = dimX * dimY * dimZ 
  ! compute the total electron in cell
!  dV = Volume(cell%cellreal)/total
!  DO is = 1, nspin
!    totEleNum = SUM(rho(:,:,:,is))*dV
!    WRITE(*,*) "dV=",dV
!    WRITE(*,*) "output total electrons number : ", totEleNum
!  ENDDO
  !------------------------------------------------------------
 

  !! >> FINISH << !!
!  DO itype = 1, atomTypeNum
!    DEALLOCATE(aden(itype)%normG)
!    DEALLOCATE(aden(itype)%rhoG)
!  ENDDO
!  DEALLOCATE(aden)
  
  RETURN

  CONTAINS

SUBROUTINE ReadAtomicDensity
!------------------------------------------------------------------------------
! DESCRIPTION:
! Read atomic density from files.
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
  INTEGER :: openStatus
  INTEGER :: afile=33 
  ! atomic file index
  INTEGER :: ig
  INTEGER :: itype 
  ! index

  !! >> INITIALIZE << !!
  CALL Title("AtomicDensity::ReadAtomicDensity")

  !! >> FUNCTION << !!
  DO itype = 1, atomTypeNum
    WRITE(outputUnit,*) "atomType=",itype," FileName=", TRIM(atomicDensityFile(itype))

    ! open the file
    OPEN(UNIT=afile,FILE=TRIM(atomicDensityFile(itype)),ACTION="read",IOSTAT=openStatus)

    IF(openStatus.NE.0) THEN
      WRITE(outputUnit,*) "File ",TRIM(atomicDensityFile(itype)), " doesn't exist"
      STOP
    ENDIF

    ! number of G points.
    READ(afile,*) aden(itype)%nG
    WRITE(outputUnit,*) "gridPoints=",aden(itype)%nG


    ! allocate g and rho(g)
    ALLOCATE(aden(itype)%normG(aden(itype)%nG))
    ALLOCATE(aden(itype)%rhoG(aden(itype)%nG))
    ! read in g and rho(g)
    DO ig = 1, aden(itype)%nG
      READ(afile,*) aden(itype)%normG(ig), aden(itype)%rhoG(ig)
      !WRITE(outputUnit,*), aden(itype)%normG(ig), aden(itype)%rhoG(ig)
    END DO

    DO ig = 1, aden(itype)%nG
      aden(itype)%normG(ig) = aden(itype)%normG(ig) * bohr
    END DO

    aden(itype)%maxG = aden(itype)%normG(aden(itype)%nG)

    CLOSE(afile)

  END DO
END SUBROUTINE ReadAtomicDensity

END SUBROUTINE GenerateAtomicDensity


FUNCTION ADenLookup(aden, qNorm)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Get the atomic density for each |q| (norm of plane wave) 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! 
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012 Created by Mohan Chen.
!------------------------------------------------------------------------------
  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  TYPE(atomDen), INTENT(IN) :: aden
  REAL(KIND=DP) :: qNorm
  REAL(KIND=DP) :: ADenLookup

  !! >> INTERNAL VARIABLES << !!
  REAL(KIND=DP) :: qSpacing
  REAL(KIND=DP) :: x0, x1, x2, x3
  REAL(KIND=DP) :: pos
  INTEGER :: iq

  !! >> INITIALIZE << !!
!  CALL Title("AtomicDensity::AtomicDensityLookup")

  !! >> FUNCTION << !!
  qSpacing = aden%normG(2) 

  !WRITE(outputUnit,*) "qspacing=",qspacing,"qnorm=",qnorm,"maxG=",aden%maxG

  ! qNorm cannot be negative ...
  IF (qNorm<0._DP) THEN
    WRITE(*,*) 'In AtomicDensityLookup(), qNorms < 0, code stop!'
    STOP
  ENDIF

  ! check some spectial cases
  IF (qNorm == 0._DP) THEN
    ADenLookup=aden%rhoG(1)
    RETURN
  ENDIF

  ! if qNorm bigger than defined maxG in psp, return 0._DP
  IF (qNorm >= aden%maxG) THEN
    ADenLookup = 0._DP
    RETURN
  ENDIF


  ! read the following words to see why I do this check here
  ! since NaN is not equal to anyting in FORTRAN
  ! we use the following to see if 4*i/qNorm is too big to case a NaN
  IF (-4._DP*PI/qNorm**2_DP .NE. -4_DP*PI/qNorm**2_DP ) THEN
    WRITE(*,*)&
     ' There is another case which should be considered: very large bulk. ', &
     '  because larger system will give denser q points and will put q points ', &
     ' very closer to q=0 and this will make -4*pi*Z/q**2 to -infinity',&
     ' If computer find -4*pi*Z/q**2 is too big, it will generate NaN',&
     ' But currently, in our group, we have not meet a system large enough to ',&
     ' cause CPU to generate NaN, if you meet such kind of problem, do something!, code STOP!!'

     STOP
  ENDIF


  ! Not meet any special cases above, OK! let's do interpolation
  pos = qNorm/qSpacing + 1._DP
  
  IF ( pos > aden%nG-3 ) THEN
    WRITE(*,*) " Not enough G points for interpolations."
    WRITE(*,*) " qNorm=",qNorm
    WRITE(*,*) " pos=",pos," nG=",aden%nG-3
    STOP
  ENDIF

  iq = FLOOR(pos)

  x0 = pos - iq
  x1 = 1.0 - x0
  x2 = 2.0 - x0
  x3 = 3.0 - x0

  ADenLookup = x1*x2*(aden%rhoG(iq)*x3+aden%rhoG(iq+3)*x0)/6.0+ &
    x0*x3*(aden%rhoG(iq+1)*x2-aden%rhoG(iq+2)*x1)/2.0

  !WRITE(outputUnit,*) "qspacing=",qspacing,"qnorm=",qnorm,"iq=",iq," value=",ADenLookup

  RETURN

END FUNCTION ADenLookup


FUNCTION ADenRecip(ionTable, elementTable, aden)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Fill the 3D recip array of atomic density.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012 Created by Mohan Chen.
!------------------------------------------------------------------------------
 
  USE MATHFUNCTIONS, ONLY : Norm
  USE CellInfo, ONLY : ion, element
  USE PlaneWave, ONLY: CCStructureFactor
  USE PlaneWave, ONLY: qVectors
  USE CellInfo, ONLY: k1G, k2G, k3G

  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!

  TYPE(ion), DIMENSION(:), INTENT(IN) :: ionTable
  ! Type + location of each ion, sorted by ion type
  ! We can get the number of ions (formerly numIon) by
  ! the size of this array.
  !
  TYPE(element), DIMENSION(:), INTENT(IN) :: elementTable
  ! index of the FIRST ion of each type in iontable
  ! This array actually has one more element than
  ! numIonType, because the last entry here numIon+1
  ! to make it easy to use this to loop through stuff
  !
  TYPE(atomDen), DIMENSION(:), INTENT(IN) :: aden  
  ! atomic density information 
  !
  !! >> INTERNAL VARIABLES << !!
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: ADenRecip
  !
  INTEGER :: ix, i2, i3 
  ! counter for parsing q over the recip. sphere
  !
  INTEGER :: i_type 
  ! Counters to parse through all types of ions
  !
  REAL(KIND=DP) :: value 
  ! atomic charge value.
  !
  REAL(KIND=DP) :: qNorm 
  ! norm of q vectors. 
  !
  REAL(KIND=DP), DIMENSION(3) :: qPoint 
  ! q point value
  !

  !! >> INITIALIZE << !!
  CALL StartClock('ADenRec')
  CALL Title("AtomicDensity::AtomicDensityRecip")
  WRITE(outputUnit,'(A)') " Generating Atomic Charge Density in Reciprocal Space"

  ADenRecip = 0.0

  !! >> FUNCTION << !!
  ! Also, we index over the size of the qTable to be sure we have the same
  ! number of q-vectors.
  DO i3 = 1, k3G
    DO i2 = 1, k2G
      DO ix = 1, k1G

        ! Calculate the qPoint cartesian coordinates from this mVector
        qPoint(:) = qVectors(ix,i2,i3,:) 
        qNorm = Norm(qPoint)

        ! Loop through all types of ions
        DO i_type = 1, SIZE(aden)

          ! Pick out the pseudopotential value for this ion type at q
          value = ADenLookup(aden(i_type), qNorm)

          ADenRecip(ix,i2,i3) = &
          ADenRecip(ix,i2,i3) + value &
            * CCStructureFactor(ionTable(elementTable(i_type)%firstIonID:&
                                         elementTable(i_type+1)%firstIonID-1), &
                                cell%cellReal, &
                                -qPoint)
        END DO
      END DO
    END DO 
  ENDDO

! TEST
!  DO ix=1,sizeX
!    DO i2=1,recipDim2
!      DO i3=1,recipDim3
!        WRITE(outputUnit,*) ADenRecip(ix,i2,i3)
!      END DO
!    END DO
!  END DO

  ADenRecip = ADenRecip / cell%vol

  CALL StopClock('ADenRec')

  RETURN

END FUNCTION ADenRecip


SUBROUTINE ADenReal(rho, aden)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the atomic density in real space (3D). 
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012 Created by Mohan Chen
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell
  USE FOURIER_NEW 

  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(INOUT) :: rho
  TYPE(atomDen), DIMENSION(:), INTENT(IN) :: aden

  !! >> INITIALIZE << !!
  CALL Title("AtomicDensity::ADenReal")

  rho = 0.d0
  CALL FFT_NEW(FFT_STD_STATE,ADenRecip(cell%ionTable, &
  cell%elementTable, aden),rho)

  RETURN

END SUBROUTINE ADenReal


END MODULE AtomicDensity

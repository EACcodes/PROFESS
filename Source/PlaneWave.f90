MODULE PlaneWave 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE PlaneWave
!     |_SUBROUTINE CCStructureFactor 
!     |_SUBROUTINE FillQTable
!     |_FUNCTION WrapQ
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
!   01/31/2013  this module is reconstructed by Mohan Chen
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP
  USE OutputFiles, ONLY: outputUnit

  IMPLICIT NONE

  REAL(KIND=DP) :: energyCutoff = -1 
  ! Recip. space plane waves cutoff in eV.
  !

!------------------------------------------------------------------------------
! The following tables are indexed by the b-vectors (the reciprocal-space
! lattice vectors).  (minus one offset as of 12/15/2003)
! Thus, Table(1,2,3, <spin>) references the tabulated value at the reciprocal
! space vector b1*(1-1) + b2*(2-1) + b3*(3-1).
! 
! The following tables are all half of a cube, to make use of the fact 
! that in reciprocal space, you only need to store half of the points, 
! because F(-x) = F(x*) when we have real objects.  (see ref. 1 above).  
! By convention, the first dimension is the one that we will store only half
! the necessary values.
!
! However, to satisfy one of the requirements of the Fast Fourier Transform, 
! the indexing is sort of strange.  Both positive and negative values of the
! indices are allowed...however, the FFT only allows positive indices.  Thus,
! we will make it so that all our dimensions START at 0, and value stored 
! will also correspond to 0 repitions of that lattice vector.  The index
! will then rise towards the maximum repititions allowed...but then, while
! the INDEX will continue going up to more positive numbers,  value stored 
! at that index will actually WRAP AROUND to the value at the most negative 
! index possible, and up until the value right before 0.
!------------------------------------------------------------------------------

! So, here's an example.  Lets say for our second index, we want to store 
! values from -3 to 3 reps of our lattice vector.  Then:
!
!   Index:                      1 2 3 4  5  6  7 
!   # Reps of Lattice vector:   0 1 2 3 -3 -2 -1

  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: qTable 
  ! A table of the norms of q-vectors.
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: qVectors
  ! Here are the actual q vectors. 
  !
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: qMask          
  ! .true. if the positive half-sphere is within the kinetic
  ! energy cutoff, .false. if otherwise.  EXCLUDES THE Q=0
  ! TERM (always .false.) which is used in ion-electron calcs!

  CONTAINS

FUNCTION CCStructureFactor(singleTypeIonTable, cellReal, qPoint)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This function calculates the complex conjugate of the structure factor for
! a given ion type, at a given q-point.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/18/2007  Function separated from ion-electron functions.  (Linda Hung)
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: IMAG
  USE CellInfo, ONLY: ion
  USE MathFunctions, ONLY: Vecmul

  IMPLICIT NONE

                          !>> EXTERNAL VARIABLES <<!

  TYPE(ion), DIMENSION(:), INTENT(IN) :: singleTypeIonTable    
  ! Type + location of all ions of a given type
  ! We can get the number of ions (formerly numIon) by
  ! the size of this array.
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cellReal              
  ! The real-space cell lattice vectors
  !
  REAL(KIND=DP), DIMENSION(3), INTENT(IN) :: qPoint                
  ! The location of the structure factor array, in
  ! reciprocal space
  !
  COMPLEX(KIND=DP) :: CCStructureFactor     
  ! Complex conjugate of the structure factor at qPoint

                          !>> INTERNAL VARIABLES <<!

  INTEGER :: i_ion
  ! Counter for ions
  !
  REAL(KIND=DP), DIMENSION(3) :: r 
  ! Absolute coordinates of an ion in real-space.

                          !>> INITIALIZATION <<!
  CCStructureFactor = (0.0_DP, 0.0_DP)

                          !>> FUNCTION BODY <<!

  ! Loops through all ions of a given type
  DO i_ion = 1, SIZE(singleTypeIonTable)

    ! Pick out the real coordinates of this ion
    r = Vecmul(cellReal, singleTypeIonTable(i_ion)%coord)

    ! Complex conjugate of structure factor at qPoint
    CCStructureFactor = CCStructureFactor + EXP(IMAG*CMPLX(DOT_PRODUCT(qPoint,r),0.0_DP))

  END DO

  RETURN

END FUNCTION CCStructureFactor


SUBROUTINE FillQTable(cellRecip)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine updates the table of norms that tells energy calculators in
!   reciprocal space whether to include a given q-point in the sum or leave it
!   out (based on whether it is within the energyCutoff sphere in recip. sp.)
!   It should be called every time the cell dimensions are altered.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Right now, this subroutine works perfectly for periodic boundary conditions.
!   However, qVectors (and therefore, qTable) is not correct for Dirichlet
!   boundaries.  We can still do calculations for some stresses (like WT and
!   IonElectron, I hope!) by setting qMask for the padded area as .FALSE.
!   I need to check if this actually makes sense though.  Calculations for 
!   PBE, LQ, HQ, WGC, and Hybrid energy functionals will require that qVectors
!   be properly padded, or other methods must be used.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/07/2003  File created.  (Vincent Ligneres)
!   11/22/2003  Cosmetic modifications (Greg Ho)
!   11/24/2003  Changing the table from booleans (in/out) to reals (norm.) 
!               (VLL)
!   12/10/2003  Cosmetic modifications (Greg Ho)
!   12/12/2003  Changed this subroutine to simplify it and use ApplyToQ.
!               (Greg Ho)
!   01/16/2003  Uhhh...I think that there might have been a bug where we 
!               set qMask to true at too many points (we missed a sqrt)
!   12/17/2007  Set up qMask to remove small bug - 'works' now for even arrays
!               (Linda Hung)
!
!------------------------------------------------------------------------------
  USE FOURIER, ONLY : offset

  USE MathFunctions, ONLY: Vecmul ! A fn. multiplies 3x3 matrix with 3d vector
  USE MathFunctions, ONLY: Norm   ! Returns the norm of a real vector

  USE CellInfo, ONLY: k1G, k2G, k3G, k3Goff ! G space FFT dimensions
  USE CellInfo, ONLY: m1G, m2G, m3G ! real space dimensions (Y)

  IMPLICIT NONE
                          !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:), INTENT(IN) :: cellRecip            
  ! the reciprocal space lattice vectors
  !
                         !>> INTERNAL VARIABLES <<! 
                         
  INTEGER :: i1, i2, i3                
  ! dummy counters parse over entire q-grid 
  !
  REAL(kind=DP), DIMENSION(3) :: mVector 
  ! GRIDPOINT coordinates in reciprocal space.
  ! according to the reciprocal lattice vectors.
  !
  REAL(KIND=DP), DIMENSION(3) :: qPoint  
  ! cartesian coordinates in reciprocal space.
  !

                           !>> INITIALIZATION <<!   
 
  CALL Title("PlaneWave:FillQTable")

  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           SETUP PLANE WAVE "

  qMask = .FALSE.
  qVectors = 1._DP
  qTable = 1._DP
      

                            !>> FUNCTION BODY <<!

  ! Now we will fill the qTable, a table of the norms of the q-vectors.
  ! Compute the matri1, put the smaller loop on the outside,
  ! for longer loops, better parallelization.
  ! Also, we index over the size of the qTable to be sure we have the same
  ! number of q-vectors.
  DO i3 = 1, k3G
    mVector(3) = i3 + k3Goff - 1
    ! Wrap around if index is greater than 1/2 the table size
    IF (mVector(3)>m3G/2) mVector(3)=mVector(3)-m3G
    DO i2 = 1, k2G
      mVector(2) = i2 - 1
      ! Wrap around if index is greater than 1/2 the table size
      IF (mVector(2)>m2G/2) mVector(2)=mVector(2)-m2G
      DO i1 = 1, k1G
        mVector(1) = i1 - 1

          ! Calculate the qPoint cartesian coordinates given by this mVector.
          qPoint = Vecmul(cellRecip, mVector)

          ! Assign this point in the qTable and qVectors.
          qVectors(i1, i2, i3, :) = qPoint(:)
          qTable(i1, i2, i3) = Norm(qPoint)

          ! Now to initialize qMask, an array the same size as qTable.
          ! For every entry in qTable where q is below the cutoff energy,
          ! make it true.  Otherwise, make it false.

          ! In FFTW, real to complex transforms of size n save complex values
          ! with a smaller x-dimension (n/2+1).  This is because the FFT of
          ! any real array H of dimension n has the Hermitian symmetry that
          !                         H_k=H_(n-k)*
          ! That is, H_0, H_1, ..., H_(n/2) are saved in memory.

          ! We set up qMask to be FALSE at the half disk avoid duplication of
          ! data.  This involves properly initializing values at k=0; if n is
          ! even, k=n/2 values are also tricky.

          ! Note that this qMask skips values at:
          !   (1,1,numZ/2+1)                if MOD(numZ,2)=0,
          !   (1,m3G/2+1,numZ/2+1)          if MOD(m3G,2)=MOD(numZ,2)=0,
          !   (numX/2+1,1,numZ/2+1)         if MOD(numX,2)=MOD(numZ,2)=0,
          !   (numX/2+1, m3G/2+1, numZ/2+1) if all dimensions are even.
          ! These points actually should be included, but shouldn't be doubled
          ! (as all others are), which causes difficulties in coding.  So think
          ! of their omission as using a KE cutoff just large enough that those
          ! points aren't relevant.

#ifdef __USE_PARALLEL
          IF ((offset==1 .AND. i1>=2) &
            .OR. (offset==0 .AND. i1>=2 .AND. i1<k1G) &
            .OR. (i1==1 .AND. i3+k3Goff>=2 .AND. i3+k3Goff<=(m3G-1)/2+1) &
            .OR. (i1==1 .AND. i3+k3Goff==1 .AND. i2 >=2 &
                 .AND. i2<=(k2G-1)/2+1) &
            .OR. (i1==1 .AND. MOD(m3G,2)==0 .AND. i3+k3Goff==m3G/2+1 &
                 .AND. i2 >=2 .AND. i2<=(k2G-1)/2+1) &
            .OR. (offset==0 .AND. i1==k1G .AND. i3+k3Goff>=2 &
                 .AND. i3+k3Goff<=(m3G-1)/2+1)&
            .OR. (offset==0 .AND. i1==k1G .AND. i3+k3Goff==1 .AND. i2 >=2 &
                 .AND. i2<=(k2G-1)/2+1) &
            .OR. (offset==0 .AND. i1==k1G .AND. MOD(m3G,2)==0 &
                 .AND. i3+k3Goff==m3G/2+1 .AND. i2 >=2 &
                 .AND. i2<=(k2G-1)/2+1)) THEN
#else
          IF ((offset==1 .AND. i1>=2) &
            .OR. (offset==0 .AND. i1>=2 .AND. i1<k1G) &
            .OR. (i1==1 .AND. i2>=2 .AND. i2<=(k2G-1)/2+1) &
            .OR. (i1==1 .AND. i2==1 .AND. i3 >=2 .AND. i3<=(k3G-1)/2+1) &
            .OR. (i1==1 .AND. MOD(k2G,2)==0 .AND. i2==k2G/2+1 &
                 .AND. i3 >=2 .AND. i3<=(k3G-1)/2+1) &
            .OR. (offset==0 .AND. i1==k1G .AND. i2>=2 &
                 .AND. i2<=(k2G-1)/2+1) &
            .OR. (offset==0 .AND. i1==k1G .AND. i2==1 .AND. i3 >=2 &
                 .AND. i3<=(k3G-1)/2+1) &
            .OR. (offset==0 .AND. i1==k1G .AND. MOD(k2G,2)==0 &
                 .AND. i2==k2G/2+1 .AND. i3 >=2 &
                 .AND. i3<=(k3G-1)/2+1)) &
          THEN
#endif
            qMask(i1, i2, i3) = .TRUE.
          ELSE
            qMask(i1, i2, i3) = .FALSE.
          END IF

      END DO
    END DO
  END DO

END SUBROUTINE FillQTable


FUNCTION WrapQ(table, x, y, z)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is simply a small helper function.
!   This function takes an array, presumably one that is smaller than the size
!     of the qtable (like the density in a coarser grid as used in the 
!     Multigrid code) and returns a qTable that is appropriate for that size,
!     if possible.
!   CAUTION:  THIS DOES NOT DO ANYTHING WITH QMASK...IF YOU WERE PLANNING ON
!             USING QMASK, YOU WILL HAVE TO WRITE A PROGRAM TO DO THIS YOURSELF
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
                    !>> EXTERNAL VARIABLES <<!  s

  INTEGER, INTENT(IN) :: x, y, z
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: table
  !
  REAL(KIND=DP), DIMENSION(x,y,z) :: WrapQ
  !

                    !>> INTERNAL VARIABLES <<!  
  INTEGER :: &
    i, j, k, & ! Dummy counters
    sizex, sizey, sizez

                      !>> INITIALIZATION <<!
  sizex = SIZE(table,1)
  sizey = SIZE(table,2)
  sizez = SIZE(table,3)

                       !>> FUNCTION BODY << 
  IF(x == sizex .AND. &
     y == sizey .AND. &
     z == sizez) THEN
    WrapQ = table
  ELSE IF(x > sizex .OR. &
          y > sizey .OR. &
          z > sizez) THEN
    STOP "***CANNOT HANDLE THIS QTABLE***"
  ELSE

    DO i = 1, x
      DO j = 1, y/2+1
        DO k = 1, z/2+1
          WrapQ(i,j,k) = table(i,j,k)
        END DO
      END DO
    END DO

    
    DO i = 1, x
      DO j = 1, y/2+1
        DO k = z/2+2, z
          WrapQ(i,j,k) = table(i,j,z+2-k)
        END DO
      END DO
    END DO

    DO i = 1, x
      DO j = y/2+2, y
        DO k = 1, z/2+1
          WrapQ(i,j,k) = table(i,y+2-j,k)
        END DO
      END DO
    END DO

    DO i = 1, x
      DO j = y/2+2, y
        DO k = z/2+2, z
          WrapQ(i,j,k) = table(i,y+2-j,z+2-k)
        END DO
      END DO
    END DO

  END IF

  RETURN

END FUNCTION WrapQ


END MODULE PlaneWave

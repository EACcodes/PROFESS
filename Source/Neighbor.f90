MODULE NearestDistance
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE NearestDistance
!     |_FUNCTION CheckNearestDistanceAtoms 
!     |_FUNCTION NNDist_noBin 
!     |_FUNCTION NNDist_noBin
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
!   08/18/2013 Seperate this from Calculator.f90 (mohan)
!
!------------------------------------------------------------------------------

                              !<< GLOBAL >>

  USE CONSTANTS, ONLY: DP, PI
  USE TIMER, ONLY : TimerStart, TimerStop, stopWatch

  IMPLICIT NONE

  LOGICAL :: checkNNDist_bin = .TRUE.  
  ! True: Nearest distance between atoms are computed with binning method
  !
  REAL(kind=DP):: nnDist = 2._DP  
  ! This is the shortest distance (bohr) that is
  ! permitted for nearest neighbor atoms.  Here 4 bohr
  ! is the default because it is the shortest bond
  ! between two Ga atoms in their alpha-Ga structure
  ! I think this should be the lower bound for any
  ! other materials in general cases, but of course
  ! you can specify this number to smaller value.

CONTAINS


FUNCTION  CheckNearestDistanceAtoms(minDistance)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Driver, calls 
!     (1) NNDist_nobin or 
!     (2) NNDist_bin
!   depends on checkNNDist_method_bin
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!    Created by Chen Huang Sep/10/2008
!    Add binning function  Chen Huang Dec/2009 
!------------------------------------------------------------------------------

  IMPLICIT NONE
  REAL(KIND=DP) :: minDistance
  INTEGER :: CheckNearestDistanceAtoms

  CALL Title("Calculator::CheckNearestDistanceAtoms")


  IF (checkNNDist_bin) THEN 
    CheckNearestDistanceAtoms = NNDist_bin(minDistance)
  ELSE
    CheckNearestDistanceAtoms = NNDist_nobin(minDistance)
  ENDIF

  RETURN
END FUNCTION  CheckNearestDistanceAtoms


FUNCTION NNDist_bin (minDistance)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   this subroutine directly computed the distance between two nearest atoms
!   in the system, to prevent the unphysical case that two atoms are two close
!   to each other
!   This subroutin first divides the system into bins, and then 
!   calculate the nearest neighour distance based on bins
!   to dramatically reduce the computational time
!   It shoule be used for large system
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!    Created by Chen Huang Sep/10/2008
!    Add binning function  Chen Huang Dec/2009 
!    bug fixed for binning Steven Jan/2011
!    bug fix for binning   Kaili Jiang Jan/2019
!------------------------------------------------------------------------------

  USE CellInfo, ONLY : cell
  USE MATHFUNCTIONS, ONLY: Norm, Cross
  USE OUTPUT, only : WrtOut

  IMPLICIT NONE
  REAL(KIND=DP) :: minDistance
  INTEGER :: NNDist_bin

  !>> internal vars <<!
    
  INTEGER  ::  ia, ib, ic, i, j, k
  CHARACTER(LEN=500) :: msg
    
  TYPE(stopWatch) :: watch1
  CHARACTER (len=500) :: message
  REAL(KIND=DP) :: & 
    sort(3,SIZE(cell%ionTable,1)), &
    aa,bb,cc, & 
    angle_ab, angle_cc, tmp_vector(3), &
    intervalx,intervaly,intervalz,& 
    tmp, &
    tmpCoord1(3), &
    tmpCoord(3), &
    tmpDiff(3), &
    tmpX(SIZE(cell%ionTable,1)),&
    tmpY(SIZE(cell%ionTable,1)),&
    tmpZ(SIZE(cell%ionTable,1))

  INTEGER :: & 
    y, x, &
    tmp_ptr, &
    next_ptr, &
    in,jn,kn, &
    bi,bj,bk, &
    nions,  &
    nbx,nby,nbz, & 
    idx(SIZE(cell%ionTable,1)), & 
    idy(SIZE(cell%ionTable,1)), & 
    idz(SIZE(cell%ionTable,1))

  LOGICAL :: & 
    xout, yout, zout, & 
    xneg, yneg, zneg

  INTEGER, ALLOCATABLE :: & 
    ptr(:,:,:) , &
    reg(:,:,:) ,&
    counter(:,:,:)
 
  LOGICAL :: test ! whether to output the bin information, added by mohan


  !>> FUNCTION BEGINS <<!

  Call Title("Calculator::NNDist_bin")

  !================================================
  ! Some special cases to skip
  IF (nnDist <= 0 ) THEN 
! Has declare this information in the subroutine checkOption
!   CALL WrtOut(6," nnDist < 0, then Skipped NNDistCheck (bin).")
    minDistance = 1e6
    NNDist_bin = 1
    RETURN
  ENDIF

  IF (SIZE(cell%ionTable,1) == 1) then
    WRITE(msg,'(A)')' (NNDist_bin) Only one atom in the box, Skip.'
    CALL WrtOut (6,msg)
    NNDist_bin = 1
    RETURN
  ENDIF

  !==============================================
  ! Define bins
  aa = Norm(cell%cellReal(:,1))
  bb = Norm(cell%cellReal(:,2))
  cc = Norm(cell%cellReal(:,3))

  intervalx = 20.d0            ! Spacing of bins along vector a, in bohr

  angle_ab = ACOS(ABS(DOT_PRODUCT(cell%CellReal(:,1), cell%CellReal(:,2))/(aa*bb)))
  intervaly = intervalx / SIN(angle_ab)

  tmp_vector = Cross(cell%cellReal(:,1),cell%cellReal(:,2))
  angle_cc = ACOS(ABS(DOT_PRODUCT(tmp_vector, cell%CellReal(:,3))/(cc*Norm(tmp_vector))))
  intervalz = intervalx / COS(angle_cc)

  nbx = ceiling(aa/intervalx)  ! # of bins in X direction
  nby = ceiling(bb/intervaly)  ! # of bins in Y direction
  nbz = ceiling(cc/intervalz)  ! # of bins in Z direction
  
  WRITe(message, *) " "
  CALL WrtOut(6,message)

  test = .TRUE.
  IF(test) THEN

    WRITE(message,'(a,F10.4,a,F10.4,a,F10.4,a)') &
      " (NNDist_bin) Cell Size                  : ", aa, " ", bb, " ", cc, " (Bohr)"
    CALL WrtOut(6,message)
    WRITE(message,'(a,F10.4,a)') " (NNDist_bin) Angles <a,b>               : ", angle_ab*180.d0/PI, " Degree"
    CALL WrtOut(6,message)
    WRITE(message,'(a,F10.4,a)') " (NNDist_bin) Angles <c,a-b plane>       : ", angle_cc*180_dp/PI, " Degree"
    CALL WrtOut(6,message)
    WRITE(message,'(a,3F10.4,a)')& 
      " (NNDist_bin) Bin Spacing                : ",& 
      intervalx,intervaly,intervalz, " (Bohr)"
    CALL WrtOut(6,message)
    WRITE(message,'(a,3I6)') & 
      " (NNDist_bin) Bin Number                 : ", nbx,nby,nbz
    CALL WrtOut(6,message)
  ENDIF

  ! get number of ions
  nions = SIZE(cell%ionTable,1)
  
  ALLOCATE(reg(nbx,nby,nbz))
  ALLOCATE(ptr(nbx,nby,nbz))
  ALLOCATE(counter(nbx,nby,nbz))

  reg = 0
  DO i=1, nions
    ! Fractional Coordinates
    tmpCoord = cell%ionTable(i)%coord
    
    ! Wrap ion fractoinal coordinates into cell box
    tmpCoord(1) = MODULO( tmpCoord(1),1._DP )
    tmpCoord(2) = MODULO( tmpCoord(2),1._DP )
    tmpCoord(3) = MODULO( tmpCoord(3),1._DP )

    ! Get the Cartesian coordinates of each atom
    tmpCoord = MATMUL(cell%cellReal, tmpCoord)

    tmpX(i) = tmpCoord(1)
    tmpY(i) = tmpCoord(2)
    tmpZ(i) = tmpCoord(3)
   
    ! Steven: the binning is kind of problematic, it should not use X,Y,Z coord.
    ! Compute which bin this atom is in?
    idx(i) = ceiling(cell%ionTable(i)%coord(1) * aa/intervalx)
    idy(i) = ceiling(cell%ionTable(i)%coord(2) * bb/intervaly)
    idz(i) = ceiling(cell%ionTable(i)%coord(3) * cc/intervalz)
!    write(*,*) "idx,y,z=",idx(i),idy(i),idz(i)
    ! Mohan: the values of idx,idy,idz may be problematic, the values might be <0,
    ! so I try to fix it.
! old version 10-26-2012
!    IF ( idx(i) == 0 )  idx(i) = 1
!    IF ( idy(i) == 0 )  idy(i) = 1
!    IF ( idz(i) == 0 )  idz(i) = 1
! new version providd by mohan.
    IF ( idx(i) <= 0 )  idx(i) = 1
    IF ( idy(i) <= 0 )  idy(i) = 1
    IF ( idz(i) <= 0 )  idz(i) = 1
    
    reg(idx(i),idy(i),idz(i)) = reg(idx(i),idy(i),idz(i)) + 1
  ENDDO

  !===================================================
  ! ptr: the starting index of each in sort() array
  Do k=1,nbz
    DO j=1,nby
      DO i=1,nbx
        IF (i==1 .and. j==1 .and. k==1)  THEN 
          ptr(1,1,1) = 1
        ELSE
          ptr(i,j,k) = next_ptr
        ENDIF
        next_ptr = ptr(i,j,k) + reg(i,j,k)
        !write(*,'(a,3I,a,I,a,I)') "ptr=",i,j,k,"  =>", ptr(i,j,k)," reg=",reg(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  !==============================================================
  ! Sort: atoms arranged according to their bins.
  counter = 0
  DO i=1,nions
    bi = idx(i)
    bj = idy(i)
    bk = idz(i)
    tmp_ptr = ptr(bi,bj,bk)
     !print *,tmp_ptr, " bi/bj/bk=", bi,bj,bk
     !write(*,*) 'nions: ',i
     !write(*,*) 'bin info: ',bi,' ',bj,' ',bk
    sort(1,counter(bi,bj,bk)+tmp_ptr) = tmpX(i)
    sort(2,counter(bi,bj,bk)+tmp_ptr) = tmpY(i)
    sort(3,counter(bi,bj,bk)+tmp_ptr) = tmpZ(i)
    counter(bi,bj,bk) = counter(bi,bj,bk) + 1
  ENDDO

  watch1 = TimerStart()
  
  !========================================================
  ! Check nearest distance
  minDistance = 1e6  ! initialized to be a big value
  
  ! Loop over bins
  DO ic = 1,nbz
    DO ib = 1,nby
      DO ia =1,nbx
        
        ! Loop over neighbour bins and itself
        DO k = -1,1 
          DO j = -1,1
            DO i = -1,1
              
              ! Neighbour bin's index
              xneg = .false.
              yneg = .false.
              zneg = .false.
              xout = .false.
              yout = .false.
              zout = .false.

              in = ia + i
              jn = ib + j
              kn = ic + k

              IF (in==nbx+1) THEN; xout = .TRUE.; in = 1; ENDIF
              IF (jn==nby+1) THEN; yout = .TRUE.; jn = 1; ENDIF
              IF (kn==nbz+1) THEN; zout = .TRUE.; kn = 1; ENDIF
              IF (in==0) THEN; xneg = .TRUE.; in = nbx; ENDIF
              IF (jn==0) THEN; yneg = .TRUE.; jn = nby; ENDIF
              IF (kn==0) THEN; zneg = .TRUE.; kn = nbz; ENDIF 

              ! Loop over atoms in bins
!              write(*,'(a,3I2,a,I2,a,I2)') & 
!                "ia ib ic    =",ia,ib,ic," reg1=",reg(ia,ib,ic), " ptr1=",ptr(ia,ib,ic)
!              write(*,'(a,3I2,a,I2,a,I2)') & 
!                "in,jn,kn    =",in,jn,kn," reg2=",reg(in,jn,kn), " ptr2=",ptr(in,jn,kn)

              DO y=0,reg(ia,ib,ic)-1
                DO x=0,reg(in,jn,kn)-1

!                  print *, " y x = ", x, y
                  IF( i==0 .and. j==0 .and. k==0 .and. y==x ) CYCLE

                  tmpCoord1 = sort(:,ptr(ia,ib,ic)+y)
                  tmpCoord  = sort(:,ptr(in,jn,kn)+x)

                  IF ( xout ) THEN
                    tmpCoord(1) = tmpCoord(1) + cell%cellReal(1,1)
                    tmpCoord(2) = tmpCoord(2) + cell%cellReal(2,1)
                    tmpCoord(3) = tmpCoord(3) + cell%cellReal(3,1)
                  ENDIF

                  IF ( yout ) THEN 
                    tmpCoord(1) = tmpCoord(1) + cell%cellReal(1,2)
                    tmpCoord(2) = tmpCoord(2) + cell%cellReal(2,2)
                    tmpCoord(3) = tmpCoord(3) + cell%cellReal(3,2)
                  ENDIF

                  IF ( zout ) THEN
                    tmpCoord(1) = tmpCoord(1) + cell%cellReal(1,3)
                    tmpCoord(2) = tmpCoord(2) + cell%cellReal(2,3)
                    tmpCoord(3) = tmpCoord(3) + cell%cellReal(3,3)
                  ENDIF

                  IF ( xneg ) THEN
                    tmpCoord(1) = tmpCoord(1) - cell%cellReal(1,1)
                    tmpCoord(2) = tmpCoord(2) - cell%cellReal(2,1)
                    tmpCoord(3) = tmpCoord(3) - cell%cellReal(3,1)
                  ENDIF

                  IF ( yneg ) THEN
                    tmpCoord(1) = tmpCoord(1) - cell%cellReal(1,2)
                    tmpCoord(2) = tmpCoord(2) - cell%cellReal(2,2)
                    tmpCoord(3) = tmpCoord(3) - cell%cellReal(3,2)
                  ENDIF

                  IF ( zneg ) THEN
                    tmpCoord(1) = tmpCoord(1) - cell%cellReal(1,3)
                    tmpCoord(2) = tmpCoord(2) - cell%cellReal(2,3)
                    tmpCoord(3) = tmpCoord(3) - cell%cellReal(3,3)
                  ENDIF

                  tmpDiff = tmpCoord1-tmpCoord
                  tmp = Norm(tmpDiff)
                  IF (tmp<minDistance) THEN
                     minDistance = tmp
                  ENDIF

                ENDDO
              ENDDO !Loop over atoms in each bin

            END DO !Loop Over neighbour bins
          END DO
        END DO

      ENDDO !Loop over bins
    ENDDO
  ENDDO

  IF (minDistance < nnDist) THEN
    WRITE(msg, '(A,Es12.4,A,Es12.4,A)') & 
        ' (NNDist_bin): there are two atoms have a minDistance of',& 
        minDistance, ' which is SMALLER than the allowed threshold:', nnDist, &
        ' , the code stop!'
    CALL WrtOut(6,msg)
    NNDist_bin = -1
  ELSE
    WRITE(msg,'(A,Es12.4,A,Es12.4,A)') & 
         ' (NNDist_bin) Min Distance in Cell       : ', & 
         minDistance,' > ',nnDist,' (Bohr)'
    CALL WrtOut(6,msg)
    NNDist_bin  =  1
  ENDIF
  WRITE(message,'(a,Es12.4,a)') &  
    ' (NNDist_bin) Time Cost                  : ', TimerStop(watch1), " sec."
  CALL WrtOut(6,message)

  DEALLOCATE(reg)
  DEALLOCATE(ptr)
  DEALLOCATE(counter)

  RETURN
END FUNCTION  NNDist_bin


FUNCTION NNDist_noBin(minDistance)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   this subroutine directly computed the distance between two nearest atoms
!   in the system, to prevent the unphysical case that two atoms are two close
!   to each other
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!    created by Chen Huang Sep/10/2008
!------------------------------------------------------------------------------

  USE CellInfo, ONLY : cell
  USE MATHFUNCTIONS, ONLY: Norm
  USE OUTPUT, only : WrtOut

  IMPLICIT NONE
  real(dp):: minDistance
  integer :: NNDist_noBin
  INTEGER :: atom1, atom2

  !>> internal vars <<!
    
  INTEGER  ::  ia, ib, i, j, k
  REAL(KIND=DP) :: tmp
  character(len=500) :: msg
  REAL(kind=DP), DIMENSION(3) :: image
  TYPE(stopWatch) :: watch1
  CHARACTER (len=500) ::  message
  REAL(KIND=DP) :: tmpMatrix(3)
  

  !>> FUNCTION BEGINS <<!
  
  IF (nnDist <= 0 ) THEN 
    CALL WrtOut(6,"nnDist < 0, then Skipped NNDistCheck.")
    minDistance = 1e6
    NNDist_noBin = 1
    RETURN

  ELSE
    watch1 = TimerStart()
    CALL WrtOut(6, "Check NNdist...")
  
    if (size(cell%ionTable,1) == 1) then
       write(msg,*) &
         '(CheckNNDist): Only one atom in the box, no need to', & 
         ' check distance between nearest atoms, skip'
       call WrtOut (6,msg)
       NNDist_noBin = 1
       return
    endif

    minDistance = 1e6  ! initialized to be a big value
    DO ia = 1,size(cell%ionTable,1)
      DO ib = ia+1,size(cell%ionTable,1)
        DO i = -1,1
          image(1) = REAL(i, kind=DP)
          DO j = -1,1
            image(2) = REAL(j, kind=DP)
            DO k = -1,1
              image(3) = REAL(k, kind=DP)
              tmpMatrix = MATMUL(cell%cellReal, (cell%ionTable(ia)%coord-cell%ionTable(ib)%coord+ image))
              tmp = Norm(tmpMatrix)
              IF (tmp<minDistance) THEN 
                minDistance = tmp
                atom1 = ia
                atom2 = ib
              ENDIF

            END DO
          END DO
        END DO
      ENDDO
    ENDDO

    if (minDistance < nnDist) then
      WRITE(msg, '(A,Es12.4,A,Es12.4,A,I6,a,I6)') & 
        '(CheckNNDist-NoBin): there are two atoms have a minDistance of',& 
        minDistance, ' which is SMALLER than the allowed threshold:', nnDist, &
        ' , the code stop! atom1=',atom1," atom2=",atom2
      CALL WrtOut(6,msg)
      NNDist_noBin  =  -1
    ELSE
      write(msg,'(A,Es12.4,A,Es12.4,A)') & 
         '(CheckNNDist-NoBin): minDistance=', & 
         minDistance,' > ',nnDist,' (bohr) OK!'
      call WrtOut(6,msg)
      NNDist_noBin  =  1
    endif
    WRITE(message,'(a,Es12.4)') &  
      'Exit NNDist-NoBin: Time Cost: ', TimerStop(watch1)
    CALL WrtOut(6,message)
    return

  ENDIF  ! IF  ( nndist )
END FUNCTION  NNDist_noBin


END MODULE NearestDistance

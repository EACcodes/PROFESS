MODULE IonElectronSpline 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE IonElectronSpline 
!     |_FUNCTION IonElectronPotRecipSpline
!     |_FUNCTION IonElectronForcesSpline
!     |_FUNCTION IonElectronStressSpline
!     |_SUBROUTINE FillQIonTable
!     |_SUBROUTINE CalculateSplineForces
!
! DESCRIPTION:
!   Calculates the potential contibutions to the forces, stress, and energies.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   Choly and Kaxiras, Physical Review B 67, 155101.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/10/2007  Created (Linda Hung) 
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!
  USE CONSTANTS, ONLY: DP
  USE CONSTANTS, ONLY: IMAG
  USE CONSTANTS, ONLY: PI

  USE MathFunctions, ONLY: Norm 
  USE MathFunctions, ONLY: Vecmul 
  USE MathFunctions, ONLY: Inverse 
  USE MathFunctions, ONLY: Volume

  USE CellInfo, ONLY: ion
  USE CellInfo, ONLY: element
  USE CellInfo, ONLY: n1G, n2G, n3G, n3Goff
  USE CellInfo, ONLY: m1G, m2G, m3G, m123G
  USE CellInfo, ONLY: k1G, k2G, k3G, k3Goff

  USE LocalPseudoPot, ONLY: PseudoPotLookup 
  USE LocalPseudoPot, ONLY: PseudoPotDiffLookup

  USE PlaneWave, ONLY: qVectors
  USE PlaneWave, ONLY: qTable
  USE PlaneWave, ONLY: qMask

  USE OutputFiles, ONLY: outputUnit

  USE CBSpline, ONLY: GetCardinalBSpline
  USE CBSpline, ONLY: BSplineProduct
  USE CBSpline, ONLY: splineOrder

  USE TIMER, ONLY : TimerStart, TimerStop, stopWatch
  USE Output, ONLY: WrtOut

  USE Fourier_NEW, ONLY: FFT_STD_STATE
  USE Fourier_NEW, ONLY: FFT_NEW

  IMPLICIT NONE

  ! in the future this term will be moved into Ewald
  LOGICAL :: iiSpline = .TRUE.  
  ! TRUE = use cardinal b-spline for ion-ion terms
  !
  LOGICAL :: ieSpline = .TRUE.  
  ! TRUE = use cardinal b-splines for ion-electron terms
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: ionPotReal 
  ! The ion-electron potential in real space.
  !
  REAL(KIND=DP) :: ieTime
  ! used to calculate the time
  !
  TYPE(stopwatch):: watch 
  ! Timer
  !
  CHARACTER (LEN=100) :: message
  !

CONTAINS

FUNCTION IonElectronPotRecipSpline(ionTable, elementTable)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the Ion-Electron potential in reciprocal space.  The Euler spline
!   is used (with cardinal B splines) to allow a O(NlnN) scaling algorithm.
!   The structure of this function differs from that of 
!   IonElectronPotentialRecip in that complex conjugate of the structure
!   factor is calculated for the entire array before multiplication with the
!   pseudopotential value.  This takes more memory, but is necessary for the
!   appropriate scaling.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE PlaneWave, ONLY: qVectors
  USE CellInfo, ONLY: cell

  IMPLICIT NONE
                            !>> EXTERNAL VARIABLES <<!

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
  COMPLEX(KIND=DP), DIMENSION(k1G,k2G,k3G) :: IonElectronPotRecipSpline
  ! The final answer
  !

                            !>> INTERNAL VARIABLES <<!
  INTEGER :: ix, i2, i3
  ! dummy counters to pares over the entire q-grid 
  !
  INTEGER :: i_type
  ! Counters to parse through all types of ions
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: qIonTable
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: structureFactor
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tmpRecip
  !
  REAL(KIND=DP) :: qNorm
  !
  INTEGER :: allocateStatus

                             !>> INITIALIZATION <<!

  WRITE(outputUnit,*) " "
  WRITE(outputUnit,'(A)') " Calculating Ion-Electron Potential (spline) ... "
  flush(outputUnit)
  watch=TimerStart()

  IonElectronPotRecipSpline = 0._DP

  ALLOCATE(qIonTable(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating qIonTable array for IonElectronPotRecipSpline"
  STOP
  END IF

  ALLOCATE(structureFactor(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating structureFactor array for IonElectronPotRecipSpline"
  STOP
  END IF

  ALLOCATE(tmpRecip(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating tmpRecip array for IonElectronPotRecipSpline"
  STOP
  END IF



                          !>> FUNCTION BODY <<!

  DO i_type = 1, SIZE(elementTable)-1

    ! Obtain array of structure factors
    CALL FillQIonTable(qIonTable, ionTable(elementTable(i_type)%firstIonID: &
                                   (elementTable(i_type+1)%firstIonID-1)), &
                       splineOrder, m1G, m2G, m3G, n3Goff)

    ! Sorry...bad notation, since the structure factor is actually
    ! given by CONJG(bSpline*FFT(qIonTable))
    CALL FFT_NEW(FFT_STD_STATE, qIonTable, structureFactor)

    ! Now iterate over gridpoints to multiply with pseudopotential value at q
    DO i3 = 1, k3G
      DO i2 = 1, k2G
        DO ix = 1, k1G 
          qNorm = Norm(qVectors(ix,i2,i3,:))

          ! Put ion-type-specific information together
          IonElectronPotRecipSpline(ix,i2,i3) = &
            IonElectronPotRecipSpline(ix,i2,i3) &
          + PseudoPotLookup(elementTable(i_type)%psp, qNorm) &
            * structureFactor(ix,i2,i3)

        END DO

      END DO
    END DO

  END DO

  CALL BSplineProduct(IonElectronPotRecipSpline, tmpRecip)
  IonElectronPotRecipSpline = CONJG(tmpRecip) / cell%vol * REAL( m123G, KIND=DP )

  ! calculating the time to do ion-electron force"
  ieTime = TimerStop(watch)
  WRITE(message,'(A, F10.3, A)') " complete.", ieTime, " s"
  CALL WrtOut(outputUnit,message)

  DEALLOCATE( qIonTable )
  DEALLOCATE( structureFactor )
  DEALLOCATE( tmpRecip )

  RETURN

END FUNCTION IonElectronPotRecipSpline


FUNCTION IonElectronForcesSpline(rhoRecip, ionTable, elementTable, cellReal)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   A procedure that computes the ion-electron forces, by taking the 
!   negative derivative of the spline approximation to the structure factor.
!   We use the useful fact for cardinal B-splines that
!     d/dx(M_n(x)) = M_(n-1)(x) - M_(n-1)(x-1)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Will be very inefficient in parallel if there ion types are not
!   homogeneous throughout the simulation domain.
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE MATHFUNCTIONS, ONLY: Inverse
  USE CONSTANTS, ONLY: hartreeToeV

  IMPLICIT NONE

                            !>> EXTERNAL VARIABLES <<!

  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoRecip
  ! Electron density in recip space (spin-independent)
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cellReal
  ! The real-space cell lattice vectors
  !
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
 
  REAL(KIND=DP), DIMENSION(SIZE(ionTable),3) :: IonElectronForcesSpline   
  ! the answer

                            !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3,3) :: cellInv
  ! Inverse of cellReal
  !
  INTEGER :: ix, i2, i3, i_type
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: Factor
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: pspRecip
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tmpRecip
  !
  INTEGER :: allocateStatus
  ! check the allocated disk


                            !>> INITIALIZATION <<!

  ! title and time
  WRITE(outputUnit,*) " "
  WRITE(outputUnit,'(A)', ADVANCE="NO") " Calculating Ion-Electron Forces (spline) ... "
  flush(outputUnit)
  watch=TimerStart()

  ALLOCATE( Factor(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating Factor array for IonElectronForcesSpline"
  STOP
  END IF

  ALLOCATE( pspRecip(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating pspRecip array for IonElectronForcesSpline"
  STOP
  END IF

  ALLOCATE( tmpRecip(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating tmpRecip array for IonElectronForcesSpline"
  STOP
  END IF

  IonElectronForcesSpline = 0._DP
  cellInv = TRANSPOSE(Inverse(cellReal))

                            !>> FUNCTION BODY <<!

  DO i_type = 1, SIZE(elementTable)-1
    pspRecip = 0._DP
    ! Calculate PSP array for the ion type
    DO i3 = 1, k3G
      DO i2 = 1, k2G 
        DO ix = 1, k1G 
          pspRecip(ix,i2,i3) = &
            PseudoPotLookup(elementTable(i_type)%psp, qTable(ix,i2,i3))
        END DO
      END DO
    END DO

    ! Remove the pspRecip(1,1,1) term
    IF(n3Goff==0) THEN
      pspRecip(1,1,1) = 0
    ENDIF

    ! Take FFT to find array not dependant on ion position
    ! old one
    CALL BSplineProduct(rhoRecip*pspRecip, tmpRecip)
    CALL FFT_NEW( FFT_STD_STATE, tmpRecip, Factor)

    CALL CalculateSplineForces(&
           IonElectronForcesSpline(elementTable(i_type)%firstIonID: &
                                   elementTable(i_type+1)%firstIonID-1,:), &
           ionTable(elementTable(i_type)%firstIonID: &
                    elementTable(i_type+1)%firstIonID-1), &
           Factor, splineOrder, cellInv, &
           n3Goff, m3G)
           
  END DO

  ! calculating the time to do ion-electron force"
  ieTime = TimerStop(watch)
  WRITE(message,'(A, F10.3, A)') " complete.", ieTime, " s"
  CALL WrtOut(outputUnit,message)

  DEALLOCATE( Factor )
  DEALLOCATE( pspRecip )
  DEALLOCATE( tmpRecip )

  RETURN

END FUNCTION IonElectronForcesSpline



FUNCTION IonElectronStressSpline(rhoRecip, ionTable, elementTable, energy, numX)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This function calculates the ion-electron stress component specified by a 
! and b.  Analogous to Choly-Kaxiras's energy calculation above.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE PlaneWave, ONLY: qVectors
  USE CellInfo, ONLY: cell

  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoRecip
  ! Electron density in reciprocal space (spin-neutral)
  !
  TYPE(ion), DIMENSION(:), INTENT(IN) :: ionTable
  ! Type + location of each ion, sorted by ion type
  ! We can get the number of ions (formerly numIon) by 
  ! the size of this array.
  !
  TYPE(element), DIMENSION(:), INTENT(IN) :: elementTable
  ! Index of the first ion of each type in iontable
  ! This array actually has one more element than
  ! numIonType, because the last entry here numIon+1
  ! to make it easy to use this to loop through stuff
  !
  REAL(KIND=DP), INTENT(IN) :: energy
  ! The ion-electron energy
  !
  INTEGER, INTENT(IN) :: numX
  ! Tot number of gridpoints in x-direction (realspace)
  !
  REAL(KIND=DP), DIMENSION(3,3) :: IonElectronStressSpline    
  ! The final answer
  !

                          !>> INTERNAL VARIABLES <<!

  INTEGER :: a, b
  ! The two indices that define the stress component
  !
  INTEGER :: ix, i2, i3
  ! counter for parsing q over the recip. sphere
  !
  INTEGER :: i_type
  ! counter for the type of ion
  !
  REAL(KIND=DP), DIMENSION(3) :: qPoint
  ! cartesian coordinates of the q vector
  !
  REAL(KIND=DP) :: qNorm
  ! The norm of q
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: qIonTable
  !
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: structureFactor       
  ! Actually the complex conjugate of the structure
  ! factor, with b-spline approximation
  !
  COMPLEX(KIND=DP) :: pspPart               
  ! the pseudopotential part of the expression.
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: qInG
  !
  !
  INTEGER :: allocateStatus
  ! check allocation disk

  
                          !>> INITIALIZATION <<!

  ! title and time
  WRITE(outputUnit,*) " "
  WRITE(outputUnit,'(A)', ADVANCE="NO") " Calculating Ion-Electron Stress (spline) ... "
  flush(outputUnit)
  watch=TimerStart()

  IonElectronStressSpline = 0._DP

  ALLOCATE( qIonTable(n1G, n2G, n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating qIonTable array for IonElectronStressSpline"
  STOP
  END IF

  ALLOCATE( structureFactor(k1G, k2G, k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating structureFactor array for IonElectronStressSpline"
  STOP
  END IF

  ALLOCATE( qInG(k1G, k2G, k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating qInG array for IonElectronStressSpline"
  STOP
  END IF

                          !>> FUNCTION BODY <<!
 
  ! Loops through all ions of a given type
  DO i_type = 1, cell%numIonType 

    ! Calculate the complex conjugate of the structure factor (approx)
    CALL FillQIonTable(qIonTable, ionTable(elementTable(i_type)%firstIonID: &
                                   (elementTable(i_type+1)%firstIonID-1)), &
                       splineOrder, numX, m2G, m3G, n3Goff)

    ! This is actually the complex conjugate of the structure factor
    ! Sorry for the bad notation!
    CALL FFT_NEW( FFT_STD_STATE, qIonTable, qInG)
    CALL BSplineProduct(qInG, structureFactor)

    ! Now iterate through gridpoints
    DO i3 = 1, k3G
!    WRITE(*,*) "ie i3=", i3
      DO i2 = 1, k2G
        DO ix = 1, k1G 
          ! Go further only if this point matters (within the KE cutoff)
          IF(qMask(ix, i2, i3)) THEN

            ! Find distance qNorm
            qPoint(:) = qVectors(ix,i2,i3,:)
            qNorm = Norm(qPoint)

            pspPart = PseudoPotDiffLookup(elementTable(i_type)%psp, qNorm)

            ! Now multiply by the components of interest.
            DO a = 1, 3
              DO b = a, 3
                IonElectronStressSpline(a,b) = IonElectronStressSpline(a,b) &
                  + qPoint(a) * qPoint(b) / qNorm * &
                    rhoRecip(ix,i2,i3) * structureFactor(ix,i2,i3) * pspPart
              END DO
            END DO

          END IF

        END DO
      END DO
    END DO

  END DO

  IonElectronStressSpline = -2._DP * IonElectronStressSpline * REAL(m123G, KIND=DP)

  DO a = 1, 3
    IonElectronStressSpline(a,a) = IonElectronStressSpline(a,a) - energy &
                                   * REAL( n3G, KIND=DP ) / REAL( m3G, KIND=DP)
  END DO

  IonElectronStressSpline = IonElectronStressSpline / cell%vol

  ! calculating the time to do ion-electron force"
  ieTime = TimerStop(watch)

  WRITE(message,'(A, F10.3, A)') " complete.", ieTime, " s"
  CALL WrtOut(outputUnit,message)

  DEALLOCATE( qIonTable )
  DEALLOCATE( structureFactor )
  DEALLOCATE( qInG )

  RETURN

END FUNCTION IonElectronStressSpline


SUBROUTINE FillQIonTable(qIonTable, singleTypeIonTable, splineOrder, &
                         totX, totY, totZ, locZOff)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This function calculates an array which, when FFTed and muliplied by bSpline,
! gives an approximation for CONJG(structure factor) of a given ion type.
! Approximation is through the method of cardinal B-splines.
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
                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(0:,0:,0:) :: qIonTable 
  ! An array containing part of structure factor appx
  !
  TYPE(ion), DIMENSION(:), INTENT(IN) :: singleTypeIonTable    
  ! Type + location of all ions of a given type
  ! We can get the number of ions of one type by 
  ! the size of this array.
  !
  INTEGER, INTENT(IN) :: splineOrder
  ! Even integer for order of spline appx
  !
  INTEGER, INTENT(IN) :: totX
  ! Total number gridpoints in x-dim (realspace)
  !
  INTEGER, INTENT(IN) :: totY
  ! Total number gridpoints in y-dimension
  !
  INTEGER, INTENT(IN) :: totZ
  ! Total number gridpoints in z-dimension
  !
  INTEGER, INTENT(IN) :: locZOff
  ! Offset in z-dimension
  !

                       !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(0:splineOrder) :: M1, M2, M3
  !
  !
  REAL(KIND=DP), DIMENSION(3) :: gridCoord
  !
  !
  INTEGER :: i_ion
  !
  !
  INTEGER :: ix, i2, i3
  !
  !
  INTEGER :: l1, l2, l3

                       !>> INITIALIZATION <<!

  qIonTable = 0._DP

                       !>> FUNCTION BODY <<!

  ! Calculate ion-based q-table for set ion/psp type
  DO i_ion = 1, SIZE(singleTypeIonTable)

    ! Identify closest gridpoint to ion position on electron grid
    gridCoord(1) = singleTypeIonTable(i_ion)%coord(1)*REAL(totX,KIND=DP)
    gridCoord(2) = singleTypeIonTable(i_ion)%coord(2)*REAL(totY,KIND=DP)
    gridCoord(3) = singleTypeIonTable(i_ion)%coord(3)*REAL(totZ,KIND=DP)

    ! Obtain M_n(x), M_n(x+1), ... M_n(n-1)
    CALL GetCardinalBSpline(M1,splineOrder,gridCoord(1)-FLOOR(gridCoord(1)))
    CALL GetCardinalBSpline(M2,splineOrder,gridCoord(2)-FLOOR(gridCoord(2)))
    CALL GetCardinalBSpline(M3,splineOrder,gridCoord(3)-FLOOR(gridCoord(3)))

    DO i3=0,splineOrder-1
      ! Coords would be mod gridCoord(3)-i3 if calculating for structure factor
      ! But here we are looking at the complex conjugate of the s.f.
      l3 = MODULO(-FLOOR(gridCoord(3))+i3, totZ)
#ifdef __USE_PARALLEL
      ! Screen whether it affects current processor (z-dim, parallel only)
      IF (l3<locZOff .OR. l3>=locZOff+SIZE(qIonTable,3)) CYCLE
#endif
      DO i2=0,splineOrder-1
        l2 = MODULO(-FLOOR(gridCoord(2))+i2, totY)
        DO ix=0,splineOrder-1
          l1 = MODULO(-FLOOR(gridCoord(1))+ix, totX)
          ! Add value if not screened out
          qIonTable(l1,l2,l3-locZOff) = qIonTable(l1,l2,l3-locZOff) &
                                      + M1(ix+1)*M2(i2+1)*M3(i3+1)
        END DO
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE FillQIonTable



SUBROUTINE CalculateSplineForces(splineForces, ionTable, Factor, splineOrder, &
                                 cellInv, locZOff, numZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This function calculates the forces using the cardinal b-spline approximation
! by taking analytical derivatives of the qIonTable array.
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

                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:) :: splineForces 
  ! An array containing part of structure factor appx
  !
  TYPE(ion), DIMENSION(:), INTENT(IN) :: ionTable
  ! Information about ions
  !
  REAL(KIND=DP), DIMENSION(0:,0:,0:), INTENT(IN) :: Factor
  !
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cellInv
  !
  !
  INTEGER, INTENT(IN) :: splineOrder
  ! Even integer for order of spline appx
  !
  INTEGER, INTENT(IN) :: locZOff
  ! Z-dimension offset (0 for serial)
  !
  INTEGER, INTENT(IN) :: numZ
  ! Total number of gridpoints in 3rd dimension
  !

                       !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(0:splineOrder) :: &
    M1, M2, M3, &
    M1forDeriv, M2forDeriv, M3forDeriv

  REAL(KIND=DP), DIMENSION(3) :: &
    qIonDeriv, &
    gridCoord

  INTEGER :: &
    i_ion, &
    numX, numY, &
    ix, iy, iz, &
    l1, l2, l3

                       !>> INITIALIZATION <<!

  splineForces = 0._DP
  numX = SIZE(Factor,1)
  numY = SIZE(Factor,2)

                       !>> FUNCTION BODY <<!

  DO i_ion = 1, SIZE(ionTable)

    ! Identify closest gridpoint to ion position on electron grid
    gridCoord(1) = ionTable(i_ion)%coord(1) * REAL(numX, KIND=DP)
    gridCoord(2) = ionTable(i_ion)%coord(2) * REAL(numY, KIND=DP)
    gridCoord(3) = ionTable(i_ion)%coord(3) * REAL(numZ, KIND=DP)

    CALL GetCardinalBSpline(M1,splineOrder,gridCoord(1)-FLOOR(gridCoord(1)))
    CALL GetCardinalBSpline(M2,splineOrder,gridCoord(2)-FLOOR(gridCoord(2)))
    CALL GetCardinalBSpline(M3,splineOrder,gridCoord(3)-FLOOR(gridCoord(3)))

    CALL GetCardinalBSpline(M1forDeriv, splineOrder-1, gridCoord(1)-FLOOR(gridCoord(1)))
    CALL GetCardinalBSpline(M2forDeriv, splineOrder-1, gridCoord(2)-FLOOR(gridCoord(2)))
    CALL GetCardinalBSpline(M3forDeriv, splineOrder-1, gridCoord(3)-FLOOR(gridCoord(3)))

    DO iz=0,splineOrder-1
      l3 = MODULO(FLOOR(gridCoord(3))-iz, numZ)

#ifdef __USE_PARALLEL
      ! Screen whether it affects current processor (z-dim, parallel only)
      IF (l3<locZOff .OR. l3>=locZOff+SIZE(Factor,3)) CYCLE
#endif
      DO iy=0,splineOrder-1
        l2 = MODULO(FLOOR(gridCoord(2))-iy, numY)

        DO ix=0,splineOrder-1
          l1 = MODULO(FLOOR(gridCoord(1))-ix, numX)

          qIonDeriv(1) = M2(iy+1) * M3(iz+1) * REAL(numX,KIND=DP) &
                         * (M1forDeriv(ix+1)-M1forDeriv(ix))
          qIonDeriv(2) = M1(ix+1) * M3(iz+1) * REAL(numY,KIND=DP) &
                         * (M2forDeriv(iy+1)-M2forDeriv(iy))
          qIonDeriv(3) = M1(ix+1) * M2(iy+1) * REAL(numZ,KIND=DP) &
                         * (M3forDeriv(iz+1)-M3forDeriv(iz))

          splineForces(i_ion,:) = &
            splineForces(i_ion,:) - Vecmul(cellInv,qIonDeriv) * Factor(l1,l2,l3-locZOff)

        END DO
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE CalculateSplineForces


END MODULE IonElectronSpline 

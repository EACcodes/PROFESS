MODULE IonElectron 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE IonElectron 
!     |_FUNCTION IonElectronEnergy          (Ion-electron, or external)
!     |_FUNCTION IonElectronEnergyReal
!     |_FUNCTION IonElectronPotentialRecip
!     |_FUNCTION IonElectronForces
!     |_FUNCTION IonElectronStress
!
! DESCRIPTION:
!   Calculates the potential contibutions to the forces, stress, and energies.
!   12/14/2007  Added Choly-Kaxiras methods for ion-electron terms.  (LH)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Watson, Stuart.  "Structural Relaxation at Defects by Ab Initio 
!       Molecular Dynamics."  (D. Phil. Thesis).  Trinity College, Oxford, 
!       1996.  Pg. 13.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
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
  USE LocalPseudoPot, ONLY: PseudoPotLookup 
  USE LocalPseudoPot, ONLY: PseudoPotDiffLookup
  USE PlaneWave, ONLY: qTable
  USE PlaneWave, ONLY: qMask
  USE PlaneWave, ONLY: qVectors
  USE OutputFiles, ONLY: outputUnit
  USE IonElectronSpline, ONLY: IonElectronPotRecipSpline

  IMPLICIT NONE

CONTAINS


FUNCTION IonElectronEnergy(IonElectPotRecip, cellVol, rhoRecip)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   A procedure that computes the ion-electron energy, given the electron 
!   density in real space.  It simply executes the equation found in 
!   ref. 1 (2.42).  Basically it just multiplies
!   the density in recip. space with the potential and sums up for all g.
!   There is only one tricky thing about this procedure ... the q=0 term must
!   be included and calculated! qCutMask takes the q=0 term out for every other
!   calculation, but for this one it must be left in.  
!   We also must double every single q to account for only using
!   the half-sphere in reciprocal space, and we CANNOT double the value at 
!   q=0.  
!
!   Presently, this procedure only does this calculation for reciprocal space.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Right now everything is done in recip space, but things can readily be
!   implemented for real-space, also.
! 
!   Also, this thing copies the FFT of the density into local memory, which
!   is probably a waste of memory.
!    
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/10/2003  Created (Greg Ho)
!    1/08/2003  Changed this procedure to not ever modify the global qMask.
!    1/30/2003  Robin convinced me that the correct formula was the real
!               part of the complex conjugate of the potential time the density
!               as implied by Ben Jesson's Thesis. (because we use +iGr for 
!               the structure factor instead of -iGr. (Greg Ho)
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

                            !>> EXTERNAL VARIABLES <<!

  COMPLEX(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: IonElectPotRecip
  ! The ion-electron potential in recip space 
  !
  COMPLEX(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoRecip
  ! Electron density in recip space (spin-independent)
  !
  REAL(kind=DP) :: cellVol
  ! The cell volume
  !
  REAL(kind=DP) :: IonElectronEnergy     
  ! The final answer - the external energy
  !

                            ! >> FUNCTION BODY <<!
  ! Calculate the ion-electron term. 
  ! We also add in the value at q=0 (qMask is assumed to have a value
  ! .false. at q=0, cuz in every other place we use it it's to be false.
  ! Except here.  Joy.) 
  ! We multiply by 2 to account for the half-sphere at every point except
  ! q=0.

  IonElectronEnergy = cellVol * 2.0_DP * REAL(SUM(IonElectPotRecip * CONJG(rhoRecip), qMask))
  IonElectronEnergy = IonElectronEnergy + cellVol * REAL(IonElectPotRecip(1,1,1)*CONJG(rhoRecip(1,1,1)))

  RETURN

END FUNCTION IonElectronEnergy


FUNCTION IonElectronEnergyReal(rhoR_SI, pot)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the coulombic energy at a given gridpoint
!   based on the real-space electron density and a potential.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/24/2003  File created.  (VLL)
!   01/12/2003  Changed the format of the formula so that instead of 
!               multiplying by the complex conjugate, it squares the real
!               and complex parts and adds them.
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoR_SI
  ! Electron density in real space (spin-independent)
  !
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: pot
  !
  !
  REAL(kind=DP) :: IonElectronEnergyReal 
  ! The coulombic electronic energy component.
  !

  IonElectronEnergyReal = SUM(rhoR_SI * pot)

  RETURN

END FUNCTION IonElectronEnergyReal


FUNCTION IonElectronPotentialRecip(cellReal, ionTable, elementTable)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Returns the Ion-Electron potential in reciprocal space.  
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell
  USE CellInfo, ONLY: k1G, k2G, k3G
  USE CellInfo, ONLY: k3Goff, m3G
  USE MPI_Functions, ONLY: rankGlobal 
  USE PlaneWave, ONLY: CCStructureFactor

  IMPLICIT NONE
                            !>> EXTERNAL VARIABLES <<!
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
  !
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: IonElectronPotentialRecip
  ! The final answer
  !

                            !>> INTERNAL VARIABLES <<!
  INTEGER :: ix, i2, i3
  ! dummy counters to pares over the entire q-grid 
  !
  INTEGER :: i_type
  ! Counters to parse through all types of ions
  !
  INTEGER :: charWritten

  REAL(KIND=DP) :: pseudoPotVal
  ! Intermediate storage for a pseudopotential value that was looked up.
  !
  REAL(KIND=DP), DIMENSION(3) :: qPoint 
  ! cartesian coordinates in reciprocal space.
  !
  REAL(KIND=DP) :: qNorm 
  ! The norm of the q-vector
  !
         
                            !>> INITIALIZATION <<!

  CALL StartClock('IonEleRecip')

  WRITE(outputUnit,'(A)') " Calculating Ion-Electron Potential "

  ! We get the dimensions of the q-grid
  IonElectronPotentialRecip = (0.0_DP, 0.0_DP)
  charWritten = 0

                             ! >> FUNCTION BODY <<!

  ! Also, we index over the size of the qTable to be sure we have the same
  ! number of q-vectors.
  DO i3 = 1, k3G
    DO i2 = 1, k2G 
      DO ix = 1, k1G 

        qPoint = -1.0_DP * qVectors(ix, i2, i3, :) 
        qNorm = Norm(qPoint)

        ! Loop through all types of ions
        DO i_type = 1, SIZE(elementTable)-1

          ! Pick out the pseudopotential value for this ion type at q
          pseudoPotVal = &
            PseudoPotLookup(elementTable(i_type)%psp, qNorm)

          ! Multiply the the structure factor (see Ben Jesson's thesis) with
          ! the pseudopotential in reciprocal space.
          !
          ! You may ask: Why the negative sign in front of qPoint? CCStructure
          ! factor returns e^iqr (complex conjugate of structure factor),
          ! but we want the actual structure factor here, e^-iqr. The structure
          ! structure factor has a negative sign because of how our FFT is
          ! defined: a forward transform has a negative sign. In the forces and
          ! stress, however, we use the complex conjugate of the structure
          ! factor. Why? Because when one calculates the ENERGY, 
          !
          !          Eie = SUM(rho(r)*Pot(r))
          !              = SUM(rho(g)*Pot(-g))
          !              = SUM(rho(-g)*Pot(g))
          ! 
          ! If we use the first expression, in real space, then 
          ! we want Pot(r) = FFT(Pot(g)), which is what we're calculating here.
          ! However, the forces and stress in reciprocal space can sometimes
          ! be obtained more conveniently by taking the derivatives of the
          ! second and third expressions instead (though of course one is free
          ! to take the first expression as well, if you wish.) In some of our
          ! stress and force subroutines, we chose the second (the one that
          ! uses Pot(-g)) because it results in a more compact form
          ! (fewer floating - signs).

          IonElectronPotentialRecip(ix,i2,i3) = &
          IonElectronPotentialRecip(ix,i2,i3) + pseudoPotVal &
            * CCStructureFactor(ionTable(elementTable(i_type)%firstIonID:&
                                         elementTable(i_type+1)%firstIonID-1), &
                                cellReal,  qPoint)


        END DO
      END DO
    END DO
  END DO

  IonElectronPotentialRecip = IonElectronPotentialRecip / cell%vol
  
  CALL StopClock('IonEleRecip')

  RETURN

END FUNCTION IonElectronPotentialRecip


FUNCTION IonElectronForces(rhoRecip, ionTable, elementTable, cellReal)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   A procedure that computes the ion-electron forces, by taking the 
!   negative derivative of the IonElectronEnergy with respect to the position
!   of the ions.
! 
!   Um...by all rights (at least by my derivation) there should be a volume
!   term in here.  Appernatly I can only get the right magnitude for the
!   forces if I don't include the volume term.  Its probably cancelled 
!   somwhere else ... we've been doing funny things with the volume for the 
!   ion-electron stuff all day now. :P
!
!   The only thing tricky is the - signs and the i's.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!    
!------------------------------------------------------------------------------
! REVISION LOG:
!   2/4/2004  Created (Greg Ho)
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: Cell
  USE CellInfo, ONLY: k1G, k2G, k3G

  IMPLICIT NONE

                            !>> EXTERNAL VARIABLES <<!

  ! Electron density in recip space (spin-independent)
  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoRecip 

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
  !
  REAL(KIND=DP), DIMENSION(SIZE(ionTable),3) :: IonElectronForces 
  ! The answer
  !

                           !>> INTERNAL VARIABLES <<!
  INTEGER :: i_ion, i_type 
  ! counter for the number of ions                            
  !
  INTEGER :: ix, i2, i3
  ! counter over reciprocal space tables
  !
  REAL(KIND=DP), DIMENSION(3, SIZE(ionTable)) :: r 
  ! ion coordinates
  !
  REAL(KIND=DP) :: vq

                           !>> INITIALIZATION <<!

  CALL StartClock('IonEleForce')

  WRITE(outputUnit,'(A)',ADVANCE="NO") " Calculating Ion-Electron Forces "
  WRITE(outputUnit, *)

  IonElectronForces = 0._DP

  ! Solve for the distances of the ions beforehand
  DO i_ion = 1, cell%numIon 
    r(:,i_ion) = Vecmul(cellReal,ionTable(i_ion)%coord)
  END DO

                            ! >> FUNCTION BODY <<!

  ! Loop over gridpoints .. this is faster

  DO i3 = 1, k3G 
    DO i2 = 1, k2G 
      DO ix = 1, k1G 

        IF(qMask(ix, i2, i3)) THEN

          ! Loop over all types of ions
          DO i_type = 1, cell%numIonType 

            vq = PseudoPotLookup(elementTable(i_type)%psp, qTable(ix, i2, i3))

            ! Loops through all ions of that type
            DO i_ion = elementTable(i_type)%firstIonID, &
                       elementTable(i_type+1)%firstIonID-1

                ! You may ask: Why no negative sign in front of qPoint? See
                ! comment in IonElectronPotentialRecip.

                IonElectronForces(i_ion,:) = IonElectronForces(i_ion,:) + &
                   qVectors(ix,i2,i3,:) * vq * &
                   AIMAG(rhoRecip(ix, i2, i3) &
                   * EXP(imag * DOT_PRODUCT(qVectors(ix, i2, i3, :), r(:,i_ion))))

            END DO

          END DO
        END IF

      END DO
    END DO
  END DO

  ! Multiply in factor of 2 due to symmetry
  IonElectronForces = IonElectronForces * 2._DP

  CALL StopClock('IonEleForce')

  RETURN

END FUNCTION IonElectronForces


FUNCTION IonElectronStress(rhoRecip, ionTable, elementTable, cellReal, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This function calculates the ion-electron stress component specified by a 
! and b.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
! We can further increase the efficiency by first store r(:) instead
! of calculating it inside CCStructureFactor
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/06/2004 Function created.  (Vincent Ligneres)
!   05/12/2014 Mohan Update.
!
!------------------------------------------------------------------------------
  USE Constants, ONLY: auToGPa 
  USE CellInfo, ONLY: Cell
  USE CellInfo, ONLY: k1G, k2G, k3G, n3G, m3G
  USE PlaneWave, ONLY: CCStructureFactor
  USE PlaneWave, ONLY: qVectors

  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3,3) :: IonElectronStress 
  ! The final answer
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoRecip 
  ! Electron density in reciprocal space (spin-neutral)
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cellReal
  !
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

                         !>> INTERNAL VARIABLES <<!


  INTEGER :: ix, i2, i3
  ! counter for parsing q over the recip. sphere
  !
  INTEGER :: i_type
  ! counter for the type of ion
  !
  INTEGER :: a, b
  ! The two indices that define the stress component
  !

  REAL(KIND=DP), DIMENSION(3) :: qPoint 
  ! cartesian coordinates of the q vector
  !
  COMPLEX(KIND=DP) :: pspPart 
  ! the pseudopotential part of the expression.
  !
  REAL(KIND=DP) :: qNorm  
  ! The norm of q (plane wave, G)
  !

  
                           !>> INITIALIZATION <<!

  WRITE(outputUnit,'(A)', ADVANCE="NO") " Calculating Ion-Electron Stress "
  CALL StartClock('IonEleStress')
  IonElectronStress = 0._DP

                           !>> FUNCTION BODY <<!

  ! Loop through all q vectors in the half-sphere.
  DO i3 = 1, k3G
    DO i2 = 1, k2G
      DO ix = 1, k1G 
        ! Go further only if this point matters (within the KE cutoff)
        IF(qMask(ix, i2, i3)) THEN

          qPoint(:) = qVectors(ix, i2, i3, :)
          qNorm = Norm(qPoint)

          pspPart = (0._DP, 0._DP)

          ! Loops through all ions of a given type. Trust me, qPoint is 
          ! positive here.
          DO i_type = 1, cell%numIonType 

            ! You may ask: Why no negative sign in front of qPoint? See
            ! comment in IonElectronPotentialRecip.
            pspPart = pspPart+PseudoPotDiffLookup(elementTable(i_type)%psp, qNorm)&
                  * CCStructureFactor(ionTable(elementTable(i_type)%firstIonID:&
                                      elementTable(i_type+1)%firstIonID-1), &
                                      cellReal, qPoint)
          END DO

          ! Finally, multiply by the components of interest.
          DO a = 1, 3
            DO b = a, 3
              IonElectronStress(a,b) = IonElectronStress(a,b) &
                + rhoRecip(ix, i2, i3) * pspPart * qPoint(a)*qPoint(b)/qNorm
            END DO
          END DO

        END IF
      END DO
    END DO
  END DO

  ! two accounts for the other part of plane wave
  IonElectronStress = -2._DP * IonElectronStress

  DO a = 1, 3
    IonElectronStress(a,a) = IonElectronStress(a,a) - energy &
                               *REAL(n3G,kind=DP)/REAL(m3G,kind=DP)
  END DO

  IonElectronStress = IonElectronStress / cell%Vol

!  WRITE(*,*) " stress for IonElectron (Gpa): "
!  WRITE(*,'(3F12.6)') IonElectronStress(1,:)*auToGPa
!  WRITE(*,'(3F12.6)') IonElectronStress(2,:)*auToGPa
!  WRITE(*,'(3F12.6)') IonElectronStress(3,:)*auToGPa

  CALL StopClock('IonEleStress')

  RETURN

END FUNCTION IonElectronStress


END MODULE IonElectron

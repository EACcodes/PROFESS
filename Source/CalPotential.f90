MODULE CalPotPlus
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE CalPotPlus
!     |_SUBROUTINE CalculatePotentialPlus 
!       |_SUBROUTINE CalculateIonPotPlus
!       |_SUBROUTINE CalculateHartreePotPlus
!       |_SUBROUTINE CalculateExcPotPlus
!       |_SUBROUTINE CalculateKEDFPotPlus
!       |_SUBROUTINE CalculateEnergy 
! 
! DESCRIPTION:
!   this module calculate the direction=dE/dphi, where phi = sqrt(rho), 
!   this direction is then used to generate the next direction.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 08-24-2013 Updated by MOHAN CHEN
!------------------------------------------------------------------------------

  USE Constants, ONLY: DP
  USE OUTPUTFILES, ONLY : outputUnit
  USE OUTPUTFILES, ONLY : errorUnit
  USE CellInfo, ONLY: n1G, n2G, n3G ! local fft dims in real space
  USE CellInfo, ONLY: k1G, k2G, k3G ! local fft dims in G space
  USE CellInfo, ONLY: cell
  USE OUTPUT, ONLY : QUIT

  IMPLICIT NONE

  REAL(KIND=DP) :: energyHC

CONTAINS

SUBROUTINE CalculatePotentialPlus(rho, optSqrt, potential, eTable)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Calculates the potential in real-space (same dims as rho),
!   stores in POTENTIAL.  Potential will be dE/d(sqrt(rho)) if optSqrt=.TRUE., 
!   dE/d(rho) if optSqrt=.FALSE.  This calculates energy if eTable is present.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE IonElectron, ONLY: IonElectronEnergyReal
  USE KEDF_DenDec, ONLY: do_den_dec
  USE KEDF_DenDec, ONLY: core_den
  USE KEDF_DenDec, ONLY: potDenDec
  USE CellInfo, ONLY: Cell, numSpin
  USE OUTPUT, ONLY : WrtOut

  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:,:), INTENT(IN) :: rho 
  ! The electron density.
  !
  LOGICAL, INTENT(IN) :: optSqrt  
  ! Take derivative relative to sqrt(rho) (instead of rho)
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G,numSpin), INTENT(OUT) :: potential
  ! Final answer potential (spin-dependent)
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G, numSpin) :: work
  !
  REAL(KIND=DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: eTable 
  ! Energies

                    !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: rhoReal_SI  
  ! Electron density in real space, spin independent
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: tempPotential
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G) :: del_rho 
  ! Delocalized electron density.
  ! For density decomposition method 
  !
  INTEGER :: isp  
  ! spin index
  LOGICAL :: calcEnergy
  ! calculate the energy for each term or not
  !
  REAL(KIND=DP), DIMENSION(9) :: locETable 
  ! Table of energies within this subroutine
  !
  REAL(KIND=DP) :: tmp, etmp           
  ! Volume element
  !

                    !>> INITIALIZATION <<!

  CALL TITLE("CalPotPlus::CalculatePotentialPlus")

  ALLOCATE(rhoReal_SI(n1G, n2G, n3G))

  ! Calculate energy only if eTable is present
  calcEnergy = .FALSE.
  IF(PRESENT(eTable)) THEN
    calcEnergy = .TRUE.
    locETable = 0._DP
  END IF

  rhoReal_SI = 0.d0
  DO isp = 1, numSpin  
     rhoReal_SI =  rhoReal_SI + rho(:,:,:,isp)
  ENDDO

                    !>> FUNCTION BODY <<!
  IF (do_den_dec==1) THEN
    del_rho = rhoReal_SI
    rhoReal_SI = rhoReal_SI + core_den ! Becomes the total density.
  ENDIF

  ! (1) ionic potential
  CALL CalculateIonPotPlus()

  IF (calcEnergy) THEN
    locETable(3) = IonElectronEnergyReal(rhoReal_SI, potential(:,:,:,1)) 
  ENDIF

  ! (2) Hartree potential
  CALL CalculateHartreePotPlus()

  ! Update the nonmagnetic parts of the potential to minority spin channel
  IF (numSpin == 2) THEN
    potential(:,:,:,2) = potential(:,:,:,1)
  ENDIF

  ! (3) Exchange correlation potential
  CALL CalculateExcPotPlus()

  ! (4) KEDF potential
  IF(do_den_dec==1) THEN
    ! Here we can not use total density, 
    ! only delocalized density is used here.
    ! because delocalized density is the optimized variable.
    CALL CalculateKEDFPotPlus(del_rho)
    CALL potDenDec(rhoReal_SI, del_rho, potential, locETable, calcEnergy)
  ELSE
    ! rhoReal_SI is the total density
    ! because the total density is the optimized variable.
    CALL CalculateKEDFPotPlus(rhoReal_SI)
  ENDIF

  ! energy term
  IF (calcEnergy) THEN
    CALL CalculateEnergy
  END IF
    
  IF(ALLOCATED(rhoReal_SI)) THEN
    DEALLOCATE(rhoReal_SI)
  ENDIF

  RETURN

CONTAINS


SUBROUTINE CalculateIonPotPlus
!------------------------------------------------------------------------------
! DESCRIPTION:
! Do FFT to transform Local pseudopotential (including structure factor)
! from G space to real space.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE IonElectronSpline, ONLY: ionPotReal

  IMPLICIT NONE

  CALL StartClock('Pot_Ion')

  !========================
  ! Ion-electron terms
  !========================
  ! mohan update this on 2014-06-30, remove the required FFT  
  !CALL FFT_NEW(FFT_STD_STATE, ionPotRecip, potential(:,:,:,1))
  potential(:,:,:,1) = ionPotReal(:,:,:)

  CALL StopClock('Pot_Ion')

  RETURN

END SUBROUTINE CalculateIonPotPlus


SUBROUTINE CalculateHartreePotPlus
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the Hartree potential.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
  USE Hartree, ONLY: JPotentialPlus  ! The Coulombic potential.

  IMPLICIT NONE

  CALL StartClock('Pot_Har')

  !========================
  ! Hartree terms
  !========================
  ! Do it in G space and FFT.
  CALL JPotentialPlus(rhoReal_SI, tempPotential, calcEnergy, locETable(4))
  potential(:,:,:,1) = potential(:,:,:,1) + tempPotential

  CALL StopClock('Pot_Har')

  RETURN

END SUBROUTINE CalculateHartreePotPlus


SUBROUTINE CalculateExcPotPlus
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the exchange-correlation energy Exc[rho(r)] and 
! potential Vxc[rho(r)] for a given charge density
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE XC_LDA, ONLY: exchangeCorrelation ! type of exchange correlation
  USE XC_LDA, ONLY: LSDA        ! type of LSDA. 1:PZ 2:PW92
  USE XC_LDA, ONLY: LDAEnergy   ! Exc from LDA using Perdew-Zunger
  USE XC_LDA, ONLY: LDAPot      ! Vxc from LDA using Perdew-Zunger
  USE XC_LDA, ONLY: LSDAPotPW92 ! Vxc from LSDA using PW92
  USE XC_LDA, ONLY: LSDAPotPZ  ! Vxc from LSDA using Perdew-Zunger
  USE XC_PBE, ONLY: PBEPot     ! GGA ex-corr pot.
  USE XC_PBE, ONLY: PBE_LibXC  ! PBE for spin-polarized
 
  IMPLICIT NONE

  CALL StartClock('Pot_Exc')

  !===========================
  !Exchange-correlation terms
  !===========================
  SELECT CASE (exchangeCorrelation)
    CASE(1) ! Local Density Approximation
      IF(numSpin == 1) THEN
        potential = potential + LDAPot(rho)
        IF (calcEnergy) locETable(5) = LDAEnergy(rho)
      ELSE  ! spin-polarized
        SELECT CASE(lsda)
          CASE(1)
            CALL LSDAPotPZ(rho,work,tmp)
          CASE(2)
            CALL LSDAPotPW92(rho,work,tmp)
        END SELECT
        potential = potential + work
        IF(calcEnergy) locETable(5) = tmp
      ENDIF

    CASE(2) ! Generalized Gradient Approximation (PBE version)
      IF(numSpin==1) THEN
        CALL PBEPot(rhoReal_SI,tempPotential,calcEnergy,locETable(5))
        potential(:,:,:,1) = potential(:,:,:,1) + tempPotential
      ELSE
        CALL PBE_LibXC(rho,work,tmp)
        potential = potential + work
        IF(calcEnergy) locETable(5) = tmp
       END IF

    CASE DEFAULT ! This case should never occur.
      CALL WrtOut(6, 'The XCF selected is not implemented. STOP')
      STOP
  END SELECT

  CALL StopClock('Pot_Exc')

END SUBROUTINE CalculateExcPotPlus


SUBROUTINE CalculateKEDFPotPlus(rhoOpt)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the KEDF potential. We highly suggest that if one wants to develope
! new KEDF, he/she should add a new file named "KEDF_XX.f90" and implement
! all the details in that file.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: kinetic
  USE KEDF_TF, ONLY: CalTF
  USE KEDF_VW, ONLY: CalVW
  USE KEDF_WT, ONLY : WTPotentialPlus
  USE KEDF_WGC, ONLY : WGCPotentialPlus   
  USE KEDF_MGP, ONLY : MGPPotentialPlus
  USE KEDF_Q, ONLY: CalLHQ
  USE KEDF_CAT, ONLY: CalCAT
  USE KEDF_HC10, ONLY: intPot 
  USE KEDF_WGCD, ONLY : DecomposeDensityKEDF 
  USE KEDF_DenDec, ONLY: potDD 
  USE KEDF_GGA, ONLY: GGA_functional
  USE KEDF_GGA, ONLY : model, GGAPotentialPlus, vWGTF
  USE KEDF_WGCD, ONLY: Fr
  USE KEDF_EvW, ONLY: Cal_EVC, Cal_EVT
  USE SYS, ONLY: bvac 
  !consider better algorithm for vacuum or not
  !
  IMPLICIT NONE

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoOpt
  CHARACTER(LEN=500) :: message
  INTEGER :: ii, jj, kk

  CALL Title("CalPotPlus::CalculateKEDFPotPlus") 
  CALL StartClock("Pot_KEDF")
  
  !=======================
  ! Kinetic energy terms 
  !=======================
  SELECT CASE(kinetic)

    !-----------------------------------------
    ! Thomas Fermi 
    !-----------------------------------------
    CASE(1)
      ! TF potential will add to previous potential,
      ! and then if we optimize sqrt(rho), we need to
      ! multiply 2sqrt(rho) after this.
      CALL CalTF(potential, rho, calcEnergy, locETable(7))
      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF

    !-----------------------------------------
    ! VW 
    !-----------------------------------------
    CASE(2)
      ! VW potential will calculate dE/d(sqrt(rho)) directly,
      ! so we multiply 2sqrt(rho) first.
      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF
      ! calculate VW potential and add it to potential. 
      CALL CalVW(potential, rho, calcEnergy, optSqrt, locETable(8))
    
    !-----------------------------------------
    ! VW + TF
    !-----------------------------------------
    CASE(3)    
      CALL CalTF(potential, rho, calcEnergy, locETable(7))
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF
      CALL CalVW(potential, rho, calcEnergy, optSqrt, locETable(8))

    !-----------------------------------------
    ! Wang-Teter 
    !-----------------------------------------
    CASE(4)
      CALL CalTF(potential, rho, calcEnergy, locETable(7))
      SELECT CASE(numSpin)
      CASE(1)
        ! Compute the Wang-Teter term
        CALL WTPotentialPlus(rhoOpt, tempPotential, calcEnergy, locETable(9), bvac)
        potential = potential + SPREAD(tempPotential, 4, numSpin)
      CASE(2)
        DO isp = 1, 2
          CALL WTPotentialPlus(2.d0*rho(:,:,:,isp), tempPotential, calcEnergy, etmp, bvac)
          potential(:,:,:,isp) = potential(:,:,:,isp) + tempPotential
          IF (calcEnergy) THEN
            locETable(9) = locETable(9) + 0.5d0 * etmp
          ENDIF
        ENDDO
      END SELECT
      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF
      CALL CalVW(potential, rho, calcEnergy, optSqrt, locETable(8))

    !-----------------------------------------
    ! WGC KEDF
    !-----------------------------------------
    CASE(5)
      CALL CalTF(potential, rho, calcEnergy, locETable(7))
      SELECT CASE(numSpin)
        CASE(1)
          CALL WGCPotentialPlus(rhoOpt, tempPotential, calcEnergy, locETable(9), bvac)
          potential = potential + SPREAD(tempPotential, 4, numSpin)
        CASE(2)
          DO isp = 1,2
            CALL WGCPotentialPlus(2.d0*rho(:,:,:,isp), tempPotential, calcEnergy, etmp, bvac)
            potential(:,:,:,isp) = potential(:,:,:,isp) + tempPotential
            locETable(9) = locETable(9) + 0.5d0 * etmp
          ENDDO
      END SELECT
      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF
      CALL CalVW(potential, rho, calcEnergy, optSqrt, locETable(8))
     
    !-----------------------------------------
    ! LQ-KEDF
    !-----------------------------------------
    CASE(7) 
      SELECT CASE(numSpin)
        CASE(1)
          ! the last parameter '.TRUE.' indicates LQ KEDF
          CALL CalLHQ(potential(:,:,:,1), rhoOpt, calcEnergy, optSqrt, &
                     locETable(7), locETable(8), locETable(9), .TRUE.)      
        CASE(2)
          STOP "not implement spin-polarized for LQ KEDF, stop."
      END SELECT
      
    !-----------------------------------------
    ! HQ - KEDF
    !-----------------------------------------
    CASE(8) 
      SELECT CASE(numSpin)
        CASE(1)
          ! the last parameter '.FALSE.' indicates HQ KEDF
          CALL CalLHQ(potential(:,:,:,1), rhoOpt, calcEnergy, optSqrt, &
                     locETable(7), locETable(8), locETable(9), .FALSE.)      
        CASE(2)
          STOP "not implement spin-polarized for HQ KEDF, stop."
      END SELECT
                   
    !-----------------------------------------
    ! CAT KEDF
    !-----------------------------------------
    CASE(10)
      SELECT CASE(numSpin)
        CASE(1)
          CALL CalCAT(potential(:,:,:,1), rhoOpt, calcEnergy, optSqrt, bvac, &
                      locETable(7), locETable(8), locETable(9))
        CASE(2)
          STOP "not implement spin-polarized for CAT KEDF, stop."
      END SELECT
      
    !-----------------------------------------
    ! Huang-Carter KEDF 
    !-----------------------------------------
    CASE(11) 
      CALL CalTF(potential, rho, calcEnergy, locETable(7))
      SELECT CASE(numSpin)
        CASE(1)    
           ! Compute the interpolation term
           tempPotential = intPot(rhoOpt, energyHC)
           IF(calcEnergy) THEN
             locETable(9)  = energyHC 
           ENDIF
           potential = potential + SPREAD(tempPotential, 4, numSpin)
        CASE(2)
          DO isp=1,numSpin
            potential(:,:,:,isp) = potential(:,:,:,isp)+intPot(2.d0*rho(:,:,:,isp), energyHC)
            ! mohan modify 2014-07-05,
            IF(calcEnergy) THEN
              locETable(9) = locETable(9) + energyHC * 0.5
            ENDIF
          ENDDO
          ! original formula to calculate the energy
          !If(calcEnergy) THEN
          !  locETable(9) = intEnergy(2.d0*rho(:,:,:,1))*0.5 + intEnergy(2.d0*rho(:,:,:,2))*0.5
          !ENDIF
        END SELECT
        ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
        IF (optSqrt) THEN
          potential = 2._DP * SQRT(rho) * potential
        ENDIF
        CALL CalVW(potential, rho, calcEnergy, optSqrt, locETable(8))

    !-----------------------------------------
    ! WGCD KEDF 
    !-----------------------------------------
    CASE (12)  ! Decompose the density
      ! only for spin-unpolarized calculation now
      CALL DecomposeDensityKEDF(rhoOpt, Fr, etmp, temppotential, bvac)
      potential = potential + SPREAD(temppotential, 4, numSpin)
      IF (calcEnergy) then
         locETable(7) = 0.d0
         locETable(8) = 0.d0
         locETable(9) = etmp
      ENDIF
      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF

    !-----------------------------------------
    ! GGA semilocal functionals 
    !-----------------------------------------
    CASE (15)  ! GGA semilocal functionals
      ! only for spin-unpolarized calculation now
      call GGAPotentialPlus(rhoOpt, temppotential, .true., etmp, GGA_functional)
      potential = potential + SPREAD(temppotential, 4, numSpin)

      IF (calcEnergy) then
         locETable(7) = etmp
         locETable(8) = 0.d0
         locETable(9) = 0.d0
      ENDIF
      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF

    !-----------------------------------------
    ! vW + G*TF KEDF
    !-----------------------------------------
    CASE (16)
      ! only for spin-unpolarized calculation now
      CALL vWGTF(rhoOpt, temppotential, .true., etmp, model)
      potential = potential + SPREAD(temppotential, 4, numSpin)

      IF (calcEnergy) then
         locETable(7) = etmp
         locETable(8) = 0.d0
         locETable(9) = 0.d0
      ENDIF
      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF

    !-----------------------------------------
    ! EvW KEDF with WT 
    !-----------------------------------------
    CASE (17)
      CALL Cal_EVT(potential, rho, calcEnergy, locETable, optSqrt)

    !-----------------------------------------
    ! EvW KEDF with WGC
    !-----------------------------------------
    CASE (18)
      CALL Cal_EVC(potential, rho, calcEnergy, locETable, optSqrt)
  
    !----------------------------------------
    ! Mi-Genova-Pavanello (MGP) 
    !---------------------------------------
    CASE (19)
      CALL CalTF(potential, rho, calcEnergy, locETable(7))
      SELECT CASE(numSpin)
      CASE(1)
        ! Compute the nonlocal term
        tempPotential = 0._DP
        CALL MGPPotentialPlus(rhoOpt, tempPotential, calcEnergy, locETable(9), bvac)
        potential = potential + SPREAD(tempPotential, 4, numSpin)
      CASE(2)
        message="MGP not implemented for open-shell systems YET. STOP."
        WRITE(errorUnit,'(A)') message
        CALL QUIT(message);
      END SELECT

      ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
      IF (optSqrt) THEN
        potential = 2._DP * SQRT(rho) * potential
      ENDIF
      CALL CalVW(potential, rho, calcEnergy, optSqrt, locETable(8))

    CASE DEFAULT  ! This case should never occur.
      message=" CalculatePotentialPlus: The KEDF selected is not implemented. STOP."
      WRITE(errorUnit,'(A)') message 
      CALL QUIT(message);
      STOP

  END SELECT

  CALL StopClock("Pot_KEDF")

  RETURN

END SUBROUTINE CalculateKEDFPotPlus


SUBROUTINE CalculateEnergy
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the energies.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

    USE MPI_Functions, ONLY: ReduceRealLevel1
    USE Ewald, ONLY: ionIonEnergy
 
    IMPLICIT NONE

    ! Finish up energy computations
    ! Sum up KE terms
    locETable(2) = SUM(locETable(7:9))
    ! Multiply all terms by rho array volume element
    eTable = locETable*cell%dV

    CALL ReduceRealLevel1( eTable, SIZE(eTable) )

    ! Add the ion-ion energy
    eTable(6) = ionIonEnergy

    ! Now, sum up the total energy
    eTable(1) = SUM(eTable(2:6))

    RETURN
END SUBROUTINE CalculateEnergy


END SUBROUTINE CalculatePotentialPlus


END MODULE CalPotPlus

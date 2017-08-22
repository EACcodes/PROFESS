MODULE KEDF_DenDec 
!----------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_DenDec 
!     |_SUBROUTINE SetupCoreDensity
!     |_SUBROUTINE AllocateCoreDensity
!     |_SUBROUTINE MakeCoreDensity 
!     |_SUBROUTINE ForceWithCoreDensity
!     |_SUBROUTINE SubstrateCoreDensity
!     |_SUBROUTINE AddCoreDensity
!     |_SUBROUTINE PulayForceDenDec
!     |_SUBROUTINE PotDenDec
!
! DESCRIPTION:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
! Phys. Rev. B 85, 045126 (2012)
! "Toward an orbital-free density functional theory of transition metals
! based on an electron density decomposition"
! Chen Huang and Emily A. Carter
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2012       Density decomposition code is first written by Chen Huang.
! 2013-08-19 Mohan add structure relaxation for density decomposition
!            method and also WGC for the delocalized density KEDF.
!------------------------------------------------------------------------------
                              !<< GLOBAL >>

  USE CONSTANTS, ONLY: DP, PI
  USE MATHFUNCTIONS, ONLY: Volume  ! Returns the volume of a 3x3 matrix
  USE MATHFUNCTIONS, ONLY: Inverse ! Returns the inverse of a 3x3 matrix 
  USE MATHFUNCTIONS, ONLY: Vecmul  ! Vector multiply...
  USE MATHFUNCTIONS, ONLY: Norm    ! Calculate the norm
  USE OutputFiles, ONLY: outputUnit 
  USE MPI_Functions

  IMPLICIT NONE

  INTEGER :: do_den_dec=0 
  ! 1: do density decomposition 0: no
  !
  INTEGER :: ncore=0 
  ! number of 1D core density array points 
  !
  REAL(KIND=DP), ALLOCATABLE :: recip_core(:,:,:) 
  ! reciprocal space 1D core density
  !
  REAL(KIND=DP), ALLOCATABLE :: core_den(:,:,:)  
  ! The localized electron density.
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: core_den_recip(:,:,:)  
  ! The localized electron density in reciprocal space
  !
  REAL(KIND=DP), ALLOCATABLE :: potDD(:,:,:)  
  ! The potential used to calculate the Pulay force 
  !
  REAL(KIND=DP) :: aTF=0.8d0
  !
  REAL(KIND=DP) :: bVW=0.2d0
  !
  REAL(KIND=DP) :: tot_charge 
  ! If transition metal are used, this is the total charge, for Ag this is 19.

  ! be allocated in ReadInputFiles.f90
  CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: atomicCoreFile 
  ! save the atomic density files be allocated in ReadInputFiles.f90
  !
CONTAINS


SUBROUTINE SetupCoreDensity
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Setup the core densities.
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
! 2013-08-19 Mohan add 
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell, m123G
  ! cell information
  !
  USE Sys, ONLY: rhoR 
  ! electronic density in real space
  !

  IMPLICIT NONE

  !! >> Internal variables << !!
  REAL(KIND=DP) :: tot_core_charge ! 
  REAL(KIND=DP) :: tot_del_charge
  REAL(KIND=DP) :: tot_charge
  REAL(KIND=DP) :: scaling
  REAL(KIND=DP) :: sum_core ! sum the core charge

  INTEGER :: numEle


  !! >> Initialize << !!
  IF(do_den_dec .NE. 1) THEN
    RETURN
  ENDIF

  CALL TITLE("KEDF_DenDec::SetupCoreDensity")

  numEle = NINT(SUM(cell%elementTable%chargeTot))

  !! >> Function << !!
  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "              SETUP CORE DENSIY IN DENSITY DECOMPOSITION METHOD "
  CALL MakeCoreDensity()
  WRITE(outputUnit,*) "Setup Calculator: Reset initial electron density ..."

  sum_core=SUM(core_den)
  CALL ReduceRealLevel1( sum_core )

  tot_core_charge=sum_core*cell%vol/REAL(m123G, KIND=DP)
  tot_del_charge = numEle - tot_core_charge

  tot_charge=SUM(rhoR)
  CALL ReduceRealLevel1( tot_charge )
  tot_charge=tot_charge*cell%vol/REAL(m123G, KIND=DP)
  scaling=tot_del_charge / tot_charge 
  rhoR = rhoR * scaling
  
  WRITE(outputUnit,*) 'Total core charge from core density file ',tot_core_charge
  WRITE(outputUnit,*) 'Scaling for delocalized density is ', scaling

! need a reduce here
! WRITE(outputUnit,*) 'Total charge after scaling: ', SUM(rhoR)*cell%vol/m123G

  RETURN
END SUBROUTINE SetupCoreDensity




SUBROUTINE AllocateCoreDensity()
!------------------------------------------------------------------------------
! DESCRIPTION:
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
!   The algorithm part is created by Chen Huang
!   12/10/2012  Mohan Chen
!
!------------------------------------------------------------------------------

  USE SYS, ONLY: rhoR      ! to allocate the array of core desity
  USE PlaneWave, ONLY: qMask

  ! Allocate memeory for core density
  IF (do_den_dec==1) THEN
    !WRITE(outputUnit,*) " the dimension for core density in G space is", &
    !  SIZE(qMask,1), SIZE(qMask,2), SIZE(qMask,3)
    
    IF( SIZE(qMask,1) .EQ. 0 ) STOP "Check 1st dim of qMask"
    IF( SIZE(qMask,2) .EQ. 0 ) STOP "Check 2nd dim of qMask"
    IF( SIZE(qMask,3) .EQ. 0 ) STOP "Check 3rd dim of qMask"
    IF( SIZE(rhoR,1) .EQ. 0 ) STOP "Check 1st dim of rhoR"
    IF( SIZE(rhoR,2) .EQ. 0 ) STOP "Check 2nd dim of rhoR"
    IF( SIZE(rhoR,3) .EQ. 0 ) STOP "Check 3rd dim of rhoR"


    ! core density
    ALLOCATE(core_den(SIZE(rhoR,1),SIZE(rhoR,2),SIZE(rhoR,3)))
    ! potential used to calculate the pulay force term
    ALLOCATE(potDD(SIZE(rhoR,1),SIZE(rhoR,2),SIZE(rhoR,3)))
    ! core density in reciprocal space
    ALLOCATE(core_den_recip(SIZE(qMask,1), SIZE(qMask,2), SIZE(qMask,3)))
  ENDIF

  RETURN

END SUBROUTINE AllocateCoreDensity



SUBROUTINE MakeCoreDensity
!------------------------------------------------------------------------------
! DESCRIPTION:
!   make up the core density, first we load the core density (in recipcal space!
!   from the disk, then we
!   Returns the Ion-Electron potential in reciprocal space.
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
!   The algorithm part is created by Chen Huang
!   12/10/2012  Mohan Chen
!
!------------------------------------------------------------------------------
  USE CellInfo
  USE MPI_Functions, ONLY: rankGlobal
  USE PlaneWave, ONLY: qMask
  USE PlaneWave, ONLY: qVectors

  USE FOURIER, ONLY: FFT ! The Fast Fourier Transform routine.
  USE PlaneWave, ONLY: CCStructureFactor
  USE NMS, ONLY: pchez,pchev

  IMPLICIT NONE

                            !>> EXTERNAL VARIABLES <<!

                            !>> INTERNAL VARIABLES <<!
  INTEGER :: ix, i2, i3
  ! dummy counters to pares over the entire q-grid
  !
  INTEGER :: i_type               
  ! Counters to parse through all types of ions
  !
  INTEGER :: charWritten
  !
  REAL(KIND=DP),DIMENSION(1) :: pseudoPotVal
  ! Intermediate storage for a pseudopotential value that was looked up.
  !
  REAL(KIND=DP),DIMENSION(1):: qNorm               
  ! The norm of the q-vector
  !
  REAL(KIND=DP), DIMENSION(3) :: qPoint
  ! cartesian coordinates in reciprocal space.
  !
  INTEGER :: ii, ierr, mpiErr
  REAL(KIND=DP), ALLOCATABLE :: rr(:,:) 
  ! tmp reciprocal space 1D core density
  !
  REAL(KIND=DP), ALLOCATABLE :: ss(:) 
  ! tmp also used for spline
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: tmp_recip(:,:,:)
  REAL(KIND=DP) :: q_space
  REAL(KIND=DP) :: sigma_wk(1)
  REAL(KIND=DP),DIMENSION(1) :: ypval
  REAL(KIND=DP) :: core_den_max
  REAL(KIND=DP) :: core_den_min
  REAL(KIND=DP) :: core_den_sum
  REAL(KIND=DP) :: cpu_grid_size

                            !>> INITIALIZATION <<!
  CALL StartClock('MakeCoreDen')

  core_den_recip = (0.0_DP, 0.0_DP)
  charWritten = 0

!  WRITE(*,*) "numIonType=", cell%numIonType

  ! we need to furthre improve this, after the first
  ! time, ncore > 0 and we will never allocate reciop_core again
  IF( ncore == 0 ) THEN

    ! Loop through all types of ions, find ncore
    DO i_type = 1, cell%numIonType 
      !! read the file containing the core density (in recipcal space)
      ! atomicCoreDesity has been allocated and read in ReadInputFiles.f90 

      IF( atomicCoreFile(i_type) == "none" ) THEN
        WRITE(outputUnit,'(/A,I5)') " No core density for type", i_type
        CYCLE
      ELSE
        WRITE(outputUnit,'(A,A)') " Open the core density file ", atomicCoreFile(i_type)
      ENDIF

      OPEN(FILE=atomicCoreFile(i_type), UNIT=111, ACTION='read')
      READ(111,*) ncore,q_space
      WRITE(outputUnit,*) "core density, nq=", ncore, " dq=", q_space
      CLOSE(111)
    ENDDO
    ! allocate the core density for each element
    ALLOCATE(recip_core(cell%numIonType, ncore,2))
  ENDIF

  ! Loop through all types of ions
  DO i_type = 1, SIZE(cell%elementTable)-1
    !! read the file containing the core density (in recipcal space)
    ! atomicCoreDesity has been allocated and read in ReadInputFiles.f90 

    IF( atomicCoreFile(i_type) == "none" ) THEN
      WRITE(outputUnit,'(/A,I5)') " No core density for type", i_type
      CYCLE
    ELSE
      WRITE(outputUnit,'(A,A)') " Open the core density file ", atomicCoreFile(i_type)
    ENDIF

    OPEN(FILE=atomicCoreFile(i_type), UNIT=111, ACTION='read')
    READ(111,*) ncore,q_space

    ! read in the reciprocal space core density
    DO ii=1,ncore
      READ(111,*) recip_core(i_type,ii,:)
    ENDDO

    ! if we don't use this, the efficiency
    ! is extermely low, I don't know why!
    ALLOCATE(rr(ncore,2))
    ALLOCATE(ss(ncore))
    rr(:,:) = recip_core(i_type,:,:)

    CLOSE(111)


    ! prepare the spline
    CALL pchez(ncore,rr(:,1),rr(:,2),ss(:),.false.,sigma_wk,1,ierr)

    ! >> FUNCTION BODY <<!

    ! Also, we index over the size of the qTable to be sure we have the same
    ! number of q-vectors.
    DO i3 = 1, k3G
      DO i2 = 1, k2G
        DO ix = 1, k1G
          ! Calculate the qPoint cartesian coordinates from this mVector
          qPoint(:) = qVectors(ix,i2,i3,:) 
          qNorm(1) = Norm(qPoint)

          ! Pick out the pseudopotential value for this ion type at q
          ! pseudoPotVal = PseudoPotLookup(elementTable(i_type)%psp, qNorm)

          !-----------------------------------------
          ! hermite interpolation (more accurate)
          !-----------------------------------------
          IF ( qNorm(1) > rr(ncore-10,1) ) THEN
            pseudoPotVal(1) = 0.d0
          ELSE
            !PCHEV evaluates a piecewise cubic Hermite or spline function.
            CALL pchev(ncore,rr(:,1),rr(:,2),ss,1,qNorm,pseudoPotVal,ypval,ierr)
            ! too slow!!
            !call pchev(ncore,recip_core(i_type,:,1),recip_core(i_type,:,2),ss(:),1,qNorm,pseudoPotVal,ypval,ierr)
          ENDIF
          
          core_den_recip(ix,i2,i3) = core_den_recip(ix,i2,i3) + pseudoPotVal(1) &
                               * CCStructureFactor( &
                                cell%ionTable( & 
                                 cell%elementTable(i_type)%firstIonID: &
                                 cell%elementTable(i_type+1)%firstIonID-1), &
                                cell%cellReal, &
                                -qPoint)

        END DO  !end sizeX
      END DO ! end i2
    END DO ! end i3

    DEALLOCATE( rr )
    DEALLOCATE( ss )

  END DO ! end itype

  core_den_recip = core_den_recip / cell%vol 

  ! this part is very very very important !!!
  ! because the input array into FFT will be changed after this subroutine.
  ! so, in order to keep core_den_recip, I need to put core_den_recip into
  ! another array to keep core_den_recip unaffected.
  ! mohan 2013-09-08
  ALLOCATE( tmp_recip(k1G, k2G, k3G) )
  tmp_recip = core_den_recip
  core_den = FFT(tmp_recip)
  DEALLOCATE(tmp_recip)

  !
  ! output information about core density in the space, and sum
  !
  core_den_max  = maxval(core_den)
  core_den_min  = minval(core_den)
  core_den_sum  = sum(core_den)
  cpu_grid_size = dble(size(core_den))
#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(core_den_max,core_den_max,  1,MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, mpiErr)
  CALL MPI_ALLREDUCE(core_den_min,core_den_min,  1,MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpiErr)
  CALL MPI_ALLREDUCE(core_den_sum,core_den_sum,  1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  CALL MPI_ALLREDUCE(cpu_grid_size,cpu_grid_size,1,MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
#endif
  WRITE(outputUnit,*) "Total core density = ", &
  core_den_sum*volume(cell%cellreal)/cpu_grid_size

  WRITE(outputUnit,*) "Minimal core density = ",core_den_min
  WRITE(outputUnit,*) "Maximal core density = ",core_den_max

  !! make the core_density positive everywhere
  WHERE (core_den<0.d0) core_den = 0.d0
  WRITE(outputUnit,*) "Core density is set to >= 0 everywhere"
  WRITE(outputUnit,*) "Setup core density finished."

  CALL StopClock('MakeCoreDen')

  RETURN

END SUBROUTINE MakeCoreDensity


SUBROUTINE SubstrateCoreDensity(rhoR)
!------------------------------------------------------------------------------
! DESCRIPTION:
! After reading the density, we need to substrate the core density,
! and the remaining delocalized density can be used to do optimization
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
! 2013-09-03 Mohan add 
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell, m123G, numSpin

  IMPLICIT NONE

  !! >> External variables << !!
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: rhoR  ! Electronic density in real-space

  !! >> Internal variables << !!
  INTEGER :: is

  !! >> Initialization << !!

  WRITE(outputUnit,*) "Substrate the core density."
  !WRITE(outputUnit,*) "charge from global rho is", SUM(rhoR)*cell%vol/m123G
  !WRITE(outputUnit,*) "charge from core rho is ", SUM(core_den)*cell%vol/m123G

  DO is=1, numSpin
    rhoR(:,:,:,is) = rhoR(:,:,:,is) - core_den(:,:,:)
  ENDDO

  RETURN
END SUBROUTINE SubstrateCoreDensity


SUBROUTINE AddCoreDensity(rhoR)
!------------------------------------------------------------------------------
! DESCRIPTION:
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
! 2013-09-03 Mohan add 
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: cell, m123G, numSpin

  IMPLICIT NONE

  !! >> External variables << !!
  REAL(KIND=DP), DIMENSION(:,:,:,:) :: rhoR  ! Electronic density in real-space

  !! >> Internal variables << !!
  INTEGER :: is

  !! >> Initialization << !!

  WRITE(outputUnit,*) "Add the core density."
!  WRITE(outputUnit,*) "charge from global rho is", SUM(rhoR)*cell%vol/m123G
!  WRITE(outputUnit,*) "charge from core rho is ", SUM(core_den)*cell%vol/m123G

  DO is=1, numSpin
    rhoR(:,:,:,is) = rhoR(:,:,:,is) + core_den(:,:,:)
  ENDDO

  RETURN
END SUBROUTINE AddCoreDensity



SUBROUTINE PulayForceDenDec(forces)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the Pulay force contributed by the ion coordinate dependent 
! core density in density decomposition scheme.
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
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: imag
  USE CONSTANTS, ONLY: hartreeToeV
  USE CONSTANTS, ONLY: bohr
  USE Fourier, ONLY: FFT 
  USE PlaneWave, ONLY: qVectors, qMask
  USE CellInfo, ONLY: Cell, m123G
  USE SYS, ONLY: rhoR
  
  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  REAL(KIND=DP), DIMENSION(:,:), INTENT(INOUT) :: forces  ! total forces

  !! >> INTERNAL VARIABLES << !!
  REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: core_den_grad
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tmp
  INTEGER :: i_type, i_ion ! index for atom element and atoms
  REAL(KIND=DP) :: Fx, Fy, Fz
  REAL(KIND=DP) :: dV

  ! used to double check the Pulay forces.
  COMPLEX(KIND=DP), ALLOCATABLE :: potDDrecip(:,:,:)
  REAL(KIND=DP) :: FxG, FyG, FzG
!  INTEGER :: i,j,k
  INTEGER :: dimX, dimY, dimZ


  !! >> INITIALIZATION << !!
  CALL StartClock("PulayFDD")
  WRITE(outputUnit,*) "Calculate the Pulay force term in density decomposition method"


  dimX = SIZE(core_den,1)
  dimY = SIZE(core_den,2)
  dimZ = SIZE(core_den,3)

  ! allocate the array for gradient of core density 
  ALLOCATE(core_den_grad(SIZE(core_den,1), SIZE(core_den,2), SIZE(core_den,3),3))
  ALLOCATE(tmp(SIZE(core_den,1), SIZE(core_den,2), SIZE(core_den,3)))


  core_den_grad = 0.D0
  tmp=0.D0

  WRITE(outputUnit,*) "size of core_den_grad is ", SIZE(core_den,1), &
  SIZE(core_den,2), SIZE(core_den,3)
  
  dV = cell%vol / m123G

  ALLOCATE(potDDrecip(SIZE(qMask,1), SIZE(qMask,2), SIZE(qMask,3)))

  ! do FFT to get reciprocal space potential
  potDDrecip = FFT(potDD)

  ! Loop through all types of ions
  DO i_type = 1, SIZE(cell%elementTable)-1
    !! read the file containing the core density (in recipcal space)
    ! atomicCoreDesity has been allocated and read in ReadInputFiles.f90 

    IF( atomicCoreFile(i_type) == "none" ) THEN
      CYCLE
    ELSE

      DO i_ion = cell%elementTable(i_type)%firstIonID, &
                 cell%elementTable(i_type+1)%firstIonID-1



        ! calculate the gradient of each atom's core density
        ! from G space to real space
        ! BE CAREFUL! the core_den_recip will not be changed here
        ! for calculating gradient, I should not use conjugate of CORE DENSITY(G)
        core_den_grad(:,:,:,1)=FFT(imag*qVectors(:,:,:,1)*core_den_recip) 
        core_den_grad(:,:,:,2)=FFT(imag*qVectors(:,:,:,2)*core_den_recip)
        core_den_grad(:,:,:,3)=FFT(imag*qVectors(:,:,:,3)*core_den_recip)

!        WHERE(core_den(:,:,:) < 1.E-7)
!          core_den_grad(:,:,:,1) = 0.d0
!          core_den_grad(:,:,:,2) = 0.d0
!          core_den_grad(:,:,:,3) = 0.d0
!        ENDWHERE

        ! force is F=dE/dR
        Fx = SUM(potDD(:,:,:)*core_den_grad(:,:,:,1)) * dV
        Fy = SUM(potDD(:,:,:)*core_den_grad(:,:,:,2)) * dV
        Fz = SUM(potDD(:,:,:)*core_den_grad(:,:,:,3)) * dV

!        Fx = 0.d0
!        Fy = 0.d0
!        Fz = 0.d0
!        DO i=1, dimX
!          DO j=1, dimY
!            DO k=1, dimZ
!              Fx = Fx + (potDD(i,j,k)) * core_den_grad(i,j,k,1) * dV
!              Fy = Fy + (potDD(i,j,k)) * core_den_grad(i,j,k,2) * dV
!              Fz = Fz + (potDD(i,j,k)) * core_den_grad(i,j,k,3) * dV
!            ENDDO
!          ENDDO
!        ENDDO

        !WRITE(*,*) " dV=", dV
        !WRITE(*,*) " cellvol=", cell%vol
     
        !IonElectronEnergy = cellVol * 2.0_DP * &
        !  REAL(SUM(IonElectPotRecip * CONJG(rhoRecip), qMask))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! multiply the results with the cell volume because in core_den_recip we
        ! already divide by volume, this is just make what core_den_recip should
        ! be in reciprocal space.
        ! 2.0 accounts for another part of plane wave
        ! we don't calculate q=(0,0,0) here, that's why we need qMask here
        FxG = REAL(SUM(CONJG(potDDrecip(:,:,:))*imag*qVectors(:,:,:,1)*(core_den_recip),qMask)) * cell%vol * 2.d0
        FyG = REAL(SUM(CONJG(potDDrecip(:,:,:))*imag*qVectors(:,:,:,2)*(core_den_recip),qMask)) * cell%vol * 2.d0
        FzG = REAL(SUM(CONJG(potDDrecip(:,:,:))*imag*qVectors(:,:,:,3)*(core_den_recip),qMask)) * cell%vol * 2.d0


        CALL ReduceRealLevel1(Fx)
        CALL ReduceRealLevel1(Fy)
        CALL ReduceRealLevel1(Fz)

        CALL ReduceRealLevel1(FxG)
        CALL ReduceRealLevel1(FyG)
        CALL ReduceRealLevel1(FzG)

        !WRITE(outputUnit,*) " Calculate the pulay force for atom ", i_ion, " in type ", i_type
        !WRITE(outputUnit,*) " Force X on H (eV/A) = ", Fx*hartreeToeV/bohr,FxG*hartreeToeV/bohr
        !WRITE(outputUnit,*) " Force Y on H (eV/A) = ", Fy*hartreeToeV/bohr,FyG*hartreeToeV/bohr
        !WRITE(outputUnit,*) " Force Z on H (eV/A) = ", Fz*hartreeToeV/bohr,FzG*hartreeToeV/bohr

        WRITE(*,*) " Calculate the pulay force for atom ", i_ion, " in type ", i_type
        WRITE(*,*) " Force X on H (eV/A) = ", Fx*hartreeToeV/bohr,FxG*hartreeToeV/bohr
        WRITE(*,*) " Force Y on H (eV/A) = ", Fy*hartreeToeV/bohr,FyG*hartreeToeV/bohr
        WRITE(*,*) " Force Z on H (eV/A) = ", Fz*hartreeToeV/bohr,FzG*hartreeToeV/bohr

        ! add to the total force
        forces(i_ion,1) = forces(i_ion,1) + Fx
        forces(i_ion,2) = forces(i_ion,2) + Fy
        forces(i_ion,3) = forces(i_ion,3) + Fz

      ENDDO

    ENDIF
  ENDDO

!  WRITE(outputUnit,*) " core density and its gradient, ssss"  
!  DO i=24, 24
!    DO j=24, 24
!      DO k=1, 160
!        WRITE(outputUnit,'(I5,4(F10.6))') k, core_den(i,j,k), core_den_grad(i,j,k,1), &
!        core_den_grad(i,j,k,2), core_den_grad(i,j,k,3)
!      ENDDO
!    ENDDO
!  ENDDO
      
 

  ! get the total potential

  ! calculate the force


  DEALLOCATE(core_den_grad)
  DEALLOCATE(potDDrecip)
  CALL StopClock("PulayFDD")

END SUBROUTINE PulayForceDenDec



SUBROUTINE ForceDenDec(forces)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the Pulay force contributed by the ion coordinate dependent 
! core density in density decomposition scheme.
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
!
!------------------------------------------------------------------------------
  USE CellInfo, ONLY: cell
  USE PlaneWave, ONLY: qMask, qTable, qVectors
  USE PlaneWave, ONLY: CCStructureFactor
  USE FOURIER, ONLY: FFT ! The Fast Fourier Transform routine.
  USE CONSTANTS, ONLY: imag
  USE CONSTANTS, ONLY: hartreeToeV
  USE CONSTANTS, ONLY: bohr
  USE NMS, ONLY: pchev


  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  REAL(KIND=DP), DIMENSION(:,:), INTENT(INOUT) :: forces  ! total forces

  !! >> INTERNAL VARIABLES << !!
  INTEGER :: i3, i2, ix
  INTEGER :: ierr
  INTEGER :: sizeX, recipDim2, recipDim3! the size of the q-grid in the 3 dimensions
  INTEGER :: i_type, i_ion
  REAL(KIND=DP), DIMENSION(3) :: qPoint ! cartesian coordinates in reciprocal space.
  REAL(KIND=DP), DIMENSION(1) :: qNorm
  REAL(KIND=DP), DIMENSION(1) :: coreDenVal ! core density value in G space
  COMPLEX(KIND=DP) :: coreDenR  ! core density value in G multiply by structure factor 
  REAL(KIND=DP), DIMENSION(1) :: ypval ! used in interpolation
  REAL(KIND=DP), ALLOCATABLE :: ddforce(:,:) ! force contributed by decomposed density
  COMPLEX(KIND=DP), ALLOCATABLE :: potDDrecip(:,:,:) ! the potential that related to core density
  REAL(KIND=DP), ALLOCATABLE :: sigma_ypp(:) ! used for interpolation to get core density
  REAL(KIND=DP), ALLOCATABLE :: rr(:,:) ! used for interpolation to get core density

  CALL StartClock('ForceDenDec')

  !! >> Initialize << !!
  sizeX = SIZE(qMask,1)
  recipDim2 = SIZE(qMask,2)
  recipDim3 = SIZE(qMask,3)

  ALLOCATE( ddforce(cell%numIon, 3 ) ) 
  ALLOCATE( potDDrecip(sizeX, recipDim2, recipDim3) )
  ALLOCATE( sigma_ypp(ncore) )
  ALLOCATE( rr(ncore, 2) )

  ! do FFT to get reciprocal space potential
  potDDrecip = FFT(potDD)

  ddforce = 0.d0

  !! >> Function << !!

  ! Loop over all types of ions
  DO i_type = 1, SIZE(cell%elementTable)-1 
    ! only calculate those terms that has decomposed density
    IF( atomicCoreFile(i_type) == "none" ) THEN
      CYCLE
    ENDIF

    rr(:,:) = recip_core(i_type,:,:)
    
    DO i3 = 1, recipDim3
      DO i2 = 1, recipDim2
        DO ix = 1, sizex

          IF(qMask(ix, i2, i3)) THEN

            ! Calculate the qPoint cartesian coordinates from this mVector
            qPoint(:) = qVectors(ix, i2, i3, :) 
            qNorm(1) = Norm(qPoint)

            !-------------------------------------------------
            ! hermite interpolation (more accurate)
            ! recip_core: the read in core density
            ! ncore: number of points in G space core density
            !-------------------------------------------------
            !IF ( qNorm > recip_core(i_type,ncore-10,1) ) THEN
            IF ( qNorm(1) > rr(ncore-10,1) ) THEN
              coreDenVal(1) = 0.d0
            ELSE
              !PCHEV evaluates a piecewise cubic Hermite or spline function.
              !CALL pchev(ncore,recip_core(i_type,:,1),recip_core(i_type,:,2),sigma_ypp(:),1,qNorm,coreDenVal,ypval,ierr)
             
              ! use this one
              CALL pchev(ncore,rr(:,1),rr(:,2),sigma_ypp,1,qNorm,coreDenVal,ypval,ierr)
            ENDIF
  
            ! Loops through all ions of that type
            DO i_ion = cell%elementTable(i_type)%firstIonID, &
                       cell%elementTable(i_type+1)%firstIonID-1

              coreDenR = coreDenVal(1) * CCStructureFactor(cell%ionTable(i_ion: i_ion), cell%cellReal, -qPoint)

              ddforce(i_ion,1) = ddforce(i_ion,1) + REAL(CONJG(potDDrecip(ix,i2,i3))*imag*qPoint(1)*coreDenR) * 2.d0
              ddforce(i_ion,2) = ddforce(i_ion,2) + REAL(CONJG(potDDrecip(ix,i2,i3))*imag*qPoint(2)*coreDenR) * 2.d0
              ddforce(i_ion,3) = ddforce(i_ion,3) + REAL(CONJG(potDDrecip(ix,i2,i3))*imag*qPoint(3)*coreDenR) * 2.d0

            END DO ! i_ion
          END IF ! qmask
        END DO ! ix
      END DO ! i2
    END DO ! i3
  END DO ! i_type

  CALL ReduceReal_array(ddforce, cell%numIon*3)

  !! print out the force
  WRITE(outputUnit,*) " Forces arise from frozen decomposed density."
  DO i_ion=1, cell%numIon
    forces(i_ion,:) = forces(i_ion,:) + ddforce(i_ion,:)
    WRITE(outputUnit,'(I5,3(A, ES15.8))') i_ion, " ", ddforce(i_ion,1)*hartreeToeV/bohr, " ", &
     ddforce(i_ion,2)*hartreeToeV/bohr, " ", ddforce(i_ion,3)*hartreeToeV/bohr
  ENDDO

  DEALLOCATE( ddforce )
  DEALLOCATE( potDDrecip )
  DEALLOCATE( sigma_ypp )
  DEALLOCATE( rr )

  CALL StopClock('ForceDenDec')

  RETURN

END SUBROUTINE ForceDenDec



SUBROUTINE PotDenDec(rhoReal_SI, del_rho, potential, locETable, calcEnergy)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Calculate the potential due to density decomposition method.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 2014-06-04 Mohan
!------------------------------------------------------------------------------

  USE KEDF_TF, ONLY: TFPotential
  USE KEDF_TF, ONLY: TFEnergy
  USE KEDF_VW, ONLY: VWPotentialSqrtPlus 

  IMPLICIT NONE

  !! >> EXTERNAL VARIABLES << !!
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT):: potential 
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN):: rhoReal_SI 
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN):: del_rho
  REAL(KIND=DP), DIMENSION(9), INTENT(INOUT) :: locETable ! Table of energies within this subroutine
  LOGICAL :: calcEnergy

  !! >> INTERNAL VARIABLES << !!
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tmpPot1
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tmpPot2
  REAL(KIND=DP) :: tmpEnergy

  !! >> INITIALIZATION << !!

  ALLOCATE(tmpPot1(SIZE(potential,1), SIZE(potential,2), SIZE(potential,3)))
  ALLOCATE(tmpPot2(SIZE(potential,1), SIZE(potential,2), SIZE(potential,3)))

  !! >> FUNCTION << !!
  ! for correction term (a*TF+b*vW)
  tmpPot1 = aTF * ( TFPotential(rhoReal_SI)-TFPotential(del_rho) )
        
  IF (calcEnergy) THEN
    locETable(7) = locETable(7) + aTF * (TFEnergy(rhoReal_SI)-TFEnergy(del_rho))
  ENDIF

  tmpPot1 = 2._DP * SQRT(del_rho) * tmpPot1
        
  ! 2) Compute vW term.
        
  ! for correction term we use KEDF: (a*TF+b*vW)
  CALL VWPotentialSqrtPlus(SQRT(rhoReal_SI), tmpPot2, calcEnergy, tmpEnergy) 
  tmpPot1 = tmpPot1 + bVW* tmpPot2*sqrt(del_rho)/sqrt(rhoReal_SI)
  locETable(8) = locETable(8) + bVW*tmpEnergy
       
  CALL VWPotentialSqrtPlus(SQRT(del_rho), tmpPot2, calcEnergy, tmpEnergy)
  tmpPot1 = tmpPot1 - bVW*tmpPot2
  locETable(8) = locETable(8) - bVW*tmpEnergy

  ! put the interaction potential into the total potential
  Potential(:,:,:,1) = Potential(:,:,:,1) + tmpPot1

  !! >> FINILIZE << !!
  DEALLOCATE(tmpPot1)
  DEALLOCATE(tmpPot2)

  RETURN

END SUBROUTINE PotDenDec


END MODULE KEDF_DenDec

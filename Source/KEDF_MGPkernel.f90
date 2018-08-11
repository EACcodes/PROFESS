MODULE KEDF_MGPkernel 
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_MGPkernel 
!     |_SUBROUTINE FillMGP
!     |_SUBROUTINE FillMGP_ReciprocalSpace
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
!
!------------------------------------------------------------------------------

  USE CONSTANTS,     ONLY: DP, PI 
  USE FOURIER,       ONLY: FFT
  USE SYS,           ONLY: rho0, rhoS, LumgpExp, LumgpFactor
  USE MathFunctions, ONLY: LindG
  USE PlaneWave,     ONLY: qTable 
  USE MPI_Functions
  USE OutputFiles
  USE OUTPUT,        ONLY: WrtOut, outputUnit
  USE KEDF_WTkernel, ONLY: keKernel
 !USE KEDF_WGCkernel, ONLY: gamma
  USE KEDF_TF,       ONLY: lambda
  USE KEDF_VW,       ONLY: mu
  USE CellInfo,      ONLY: cell

  IMPLICIT NONE

  REAL(KIND=DP) :: tkF, cMGP             

CONTAINS



SUBROUTINE FillMGP()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine initializes the kernel for the MGP kinetic energy

! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
! JCP 148, 184107 (2018)
!------------------------------------------------------------------------------
! REVISION LOG:
!------------------------------------------------------------------------------

  IMPLICIT NONE

                      !>> INITIALIZATION <<!

  cMGP = pi**2/(3._DP * pi**2)**(1._DP/3._DP) 

                       !>> FUNCTION BODY << 

  ! This is the usual case with the wt kernel in reciprocal space
  ! Periodic boundary condition
 CALL FillMGP_ReciprocalSpace()

  RETURN

END SUBROUTINE FillMGP


SUBROUTINE FillMGP_ReciprocalSpace()
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

  INTEGER :: i_var, MaxPoints    
  INTEGER :: i3,i2,ix    
  REAL(KIND=DP) :: t_var, deltat, kertmp, q, q2, myrho
  REAL(KIND=DP) :: dt=0.00001 
  WRITE(outputUnit,'(a)') " Fill the 3D Mi-Genova-Pavanello kernel in G space"
 ! fill the Mi-Genova-Pavanello kernel
  myrho=rho0
  WRITE(6,'(a,f10.4)') ' RHO0 in MGP kernel: ', myrho 

  tkF = 2._DP * (3._DP * myrho * pi**2)**(1._DP/3._DP) 
 ! Also can make this a input variable
  MaxPoints = 500

  DO i3=1, SIZE(keKernel,3)
    DO i2=1, SIZE(keKernel,2)
      DO ix=1, SIZE(keKernel,1)

        t_var  = 1._DP/(MaxPoints)
        deltat = 1._DP/(MaxPoints)
        kertmp = 0._DP

        q2=qTable(ix,i2,i3)**2
        q=qTable(ix,i2,i3)
     !====================Hypercorrelation term =================================================   
       hypercorrelation_: do i_var = 1,MaxPoints
           kertmp = kertmp +  0.5_DP*((LindG(q/(tkF*(t_var+dt)**(1._DP/3._DP)),-0.6_DP,1._DP)&
              -LindG(q/(tkF*(t_var-dt)**(1._DP/3._DP)),-0.6_DP,1._DP))/dt)*t_var**(5.0_DP/6.0_DP)

           t_var = t_var + deltat
       end do hypercorrelation_       

       KeKernel(ix,i2,i3,2) = -1.2_DP*kertmp*deltat
     !=========================================================================================== 

     !========Kinetic Electron term ===============================================
          Kekernel(ix,i2,i3,3) = 4*pi*erf(q)**2*LumgpFactor*exp(-q2*LumgpEXP)/q2/CMGP
     !=============================================================================
     
     !==========Sum of Three terms ============================================================================ 
          KeKernel(ix,i2,i3,1) = 1.2_DP*LindG(q/tkF,1.0_DP,1._DP) + KeKernel(ix,i2,i3,3) + Kekernel(ix,i2,i3,2)
     !=========================================================================================================
         ! It already does, but let's impose correct q = 0 limit
          IF ( q2 < 1.0e-8_DP ) KeKernel(ix,i2,i3,:) = 0._DP
 
         END DO ! ix
       END DO ! i2
     ENDDO ! i3
   keKernel = keKernel * cMGP
  
  RETURN

END SUBROUTINE FillMGP_ReciprocalSpace

END MODULE KEDF_MGPkernel 

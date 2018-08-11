MODULE KEDF_CAT
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_CAT  
!     |_SUBROUTINE CalCAT
!     |_FUNCTION CATEnergy 
!     |_FUNCTION CATPot 
!     |_FUNCTION CATPotentialPlus 
!     |_FUNCTION CATrhot 
!     |_FUNCTION CATvacEnergy 
!     |_FUNCTION CATvacPot
!     |_SUBROUTINE FillCAT
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES: 
! David Garcia-Aldea and J.E. Alvarellos's Phys. Rev. A 76, 052504 (2007) 
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE Constants, ONLY: DP
  USE OutputFiles, ONLY: outputUnit
  USE KEDF_WGC
  USE Output, ONLY: WrtOut

  IMPLICIT NONE

  ! I have found that cat_beta = 2/3, cat_alpha=0 will make w_inf to be single-value=>
  ! w_inf = -1/3, the only thing you can tune is cat_gamma now
  ! and the optimum cat_gamma for alumimum diamod and fcc is 1.4
  ! then I tested the surface energies of aluminum fcc(111/100/110)
  ! = 943/1156/1177 (mj/m^2), nice!
  ! in fact, in David Garcia_aldea and JE Alvarellos paper, PRA, (2007)
  ! they also found that beta =2/3 is the best value, and when beta=2/3
  ! the CAT becomes same as the WGC KEDF.
  REAL(KIND=DP) :: cat_alpha=0.0_DP
  ! the lowest bound is vW term
  !
  REAL(KIND=DP) :: cat_beta =2._DP/3._DP
  ! this is determined by w(inf) has only one value
  !
  REAL(KIND=DP) :: cat_gamma=1.4_DP           
  ! the only parameter that can be tuned
  !

CONTAINS

SUBROUTINE CalCAT(potential, rho, calcEnergy, optSqrt, bvac, eTF, eVW, eNL)
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

  USE KEDF_TF, ONLY: TFPotential
  USE KEDF_TF, ONLY: TFEnergy
  USE KEDF_VW, ONLY: VWPotentialSqrtPlus
  USE CellInfo, ONLY: n1G, n2G, n3G
  USE CellInfo, ONLY: cell 
  
  IMPLICIT NONE

  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G), INTENT(OUT) :: potential
  ! potential we want for the answer
  !
  REAL(KIND=DP), DIMENSION(n1G, n2G, n3G), INTENT(IN) :: rho
  ! rho for optimization
  !
  LOGICAL :: calcEnergy
  ! calculate the energies for KEDF or not
  !
  LOGICAL :: optSqrt
  ! optimize sqrt(rho) if optSqrt=.TRUE.
  !
  LOGICAL :: bvac
  ! use vacuum corrections or not
  !
  REAL(KIND=DP), INTENT(OUT) :: eTF
  ! KEDF energy for Thomas-Fermi term
  !
  REAL(KIND=DP), INTENT(OUT) :: eVW
  ! KEDF energy for VW term
  !
  REAL(KIND=DP), INTENT(OUT) :: eNL
  ! KEDF energy for non-local term
  !
  !! >> LOCAL VARIABLSE << !!
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tempPotential
  !
  INTEGER :: allocateStatus
  !

  !! >> INITIALIZATION << !!

  ALLOCATE(tempPotential(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*)'Error allocating tempPotential in KEDF_CAT. Leaving.'
    STOP
  END IF


  ! Compute the TF term.
  potential = potential - cat_alpha* TFPotential(rho)

  IF (calcEnergy) THEN
    eTF = -cat_alpha*TFEnergy(rho)
  ENDIF

  ! Compute the CAT term
  IF (bvac) THEN
    potential = potential + (1._DP+cat_alpha) * CATvacPot(rho)
    IF (calcEnergy) THEN
      eNL = (1._DP+cat_alpha)*CATvacEnergy(rho)
    ENDIF
  ELSE
    potential = potential + (1._DP+cat_alpha) * CATPot(rho)
    IF (calcEnergy) THEN
      eNL = (1._DP+cat_alpha)*CATEnergy(rho)
    ENDIF
  ENDIF

!---------------------------------------------------------------------------------------------
! The potentialplus subroutine is still being debugged
!      CALL CATPotentialPlus(rho, tempPotential, calcEnergy, locETable(9), bvac)
!      potential = potential  + (1._DP + cat_alpha) * SPREAD(tempPotential, 4, cell%numSpin)
!---------------------------------------------------------------------------------------------

  ! Chain rule: dE/d(sqrt(rho)) = dE/d(rho) * 2*sqrt(rho)
  IF (optSqrt) THEN
    potential = 2._DP * SQRT(rho) * potential
  ENDIF

  ! Compute vW term.
  CALL VWPotentialSqrtPlus(SQRT(rho), tempPotential, calcEnergy, eVW)

  ! Chain rule: dE/d(rho) = dE/d(sqrt(rho)) / (2*sqrt(rho))
  IF (.NOT. optSqrt) THEN
    tempPotential = tempPotential/(2._DP*SQRT(rho))
  ENDIF

  potential = potential + tempPotential 

  DEALLOCATE(tempPotential)

  RETURN

END SUBROUTINE CalCAT


FUNCTION CATEnergy(rho) 
!------------------------------------------------------------------------------
! DESCRIPTION:
! This is CAT KEDF defined in the equation (24) in 
! David Garcia-Aldea and J.E. Alvarellos's  paper (PRB)2007
! The kernel is expaned to the 1st order
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP) :: CATEnergy
  REAL(KIND=DP) :: rho(:,:,:)

                     !>> INTERNAL VARIABLES <<!
  REAL(KIND=DP) :: rhot(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3))
  REAL(KIND=DP) :: third = 1._DP/3._DP
  
  
  !!!! FORTRAN: it is illegal to raise a negative number with real value power ....
  !!!! so we squre rhot first to make it positive and then take the 1/(3*beta) exponent

  rhot = CATrhot(rho)
  rhot = rhot**2
  CATEnergy = SUM(rho * rhot**(third/CAT_beta))
  CATEnergy = CATEnergy*cTF

  RETURN
END FUNCTION CATEnergy


FUNCTION CATPot(rho) 
!------------------------------------------------------------------------------
! DESCRIPTION:
! This is CAT KEDF defined in the equation (24) in 
! David Garcia-Aldea and J.E. Alvarellos's  paper (PRB)2007
! The kernel is expaned to the 1st order
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP) :: rho(:,:,:)
  REAL(KIND=DP) :: CATPot(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3))

                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: rhoA(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3))
  REAL(KIND=DP) :: rhot(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3))
  REAL(KIND=DP) :: rhoCATbetaMinusOne(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3))
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FFTtmp
  REAL(KIND=DP) :: third = 1._dp/3._dp
  
  !!!!!!FORTRAN: it is illegal to raise a negative number with real value power ....
  
  rhot = CATrhot(rho)
  rhoA = rho*(rhot**2)**(third/CAT_beta)/rhot*2._DP/3._DP/CAT_beta
  FFTtmp = FFT(rhoA)
  rhoCATbetaMinusOne = rho**(CAT_beta-1)
  
  !0:
  CATPot = (rhot**2)**(third/CAT_beta)
  !1:
  CATPot = CATpot+CAT_beta*rhoCATbetaMinusOne*& 
                     ( FFT(FFTtmp*keKernel(:,:,:,1) + FFT(rhoA*(rho-rhoS))*keKernel(:,:,:,2)) &
                       + (rho-rhoS)*FFT(FFTtmp*keKernel(:,:,:,2)) )
                
  !2:       
  CATPot = CATPot+rho**CAT_beta*FFT(FFTtmp*keKernel(:,:,:,2)) & 
                 +rhoA*FFT(FFT(rho**CAT_beta)*keKernel(:,:,:,2))

  ! the off-diagonal term in the CAT
  IF (WGCT==-3) THEN
    CATpot = CATpot + & 
      (rho**cat_beta * (cat_beta+1) - rhoS * rho**(cat_beta-1) * cat_beta) * FFT(FFT((rho-rhoS)*rhoA)*keKernel(:,:,:,4)) + & 
      rhoA * FFT(FFT(rho**cat_beta*(rho-rhoS))*keKernel(:,:,:,4))
  ENDIF

  CATPot = CATPot*cTF
  
  RETURN
END FUNCTION CATPot


SUBROUTINE CATPotentialPlus(rho, potential, calcEnergy, energy, vacCutoff) 
!------------------------------------------------------------------------------
! DESCRIPTION: 
! This is CAT KEDF defined in the equation (24) in 
! David Garcia-Aldea and J.E. Alvarellos's  paper (PRB)2007
! The kernel is expaned to the 1st order
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!------------------------------------------------------------------------------

  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: rho
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential
  LOGICAL, INTENT(IN) :: calcEnergy, vacCutoff
  REAL(KIND=DP), INTENT(OUT) :: energy

                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: cutf
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: cutfp
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: rhoX
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: rhot

  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FFTtmp

  REAL(KIND=DP), PARAMETER :: third = 1._dp/3._dp

                         !>> INITIALIZATION <<! (5 FFTs)

  ! Temporary storage
  IF (vacCutoff) THEN
    cutf =  CutoffFunc(rho)
    cutfp = CutoffFuncD(rho)
    rhoX = rho**CAT_beta*cutf
  ELSE
    rhoX = rho**CAT_beta
  END IF
  FFTtmp = FFT(rhoX)
  potential = FFT(FFTtmp*keKernel(:,:,:,2))

  rhot = FFT(FFTtmp*keKernel(:,:,:,1) &
           + FFT(rhoX*(rho-2._DP*rhoS))*keKernel(:,:,:,2) &
           + rho*potential)

  !Note: in FORTRAN, cannot raise a negative number to real value power
  
  rhoX = rho*(rhot**2)**(third/CAT_beta)/rhot*2._DP/3._DP/CAT_beta

                         !>> FUNCTION BODY <<!
  
  ! Calculate potential (part 1)
  rhot = rhot**2
  potential = rhoX*potential + (rhot)**(third/CAT_beta)

  ! Calculate energy
  IF (calcEnergy) energy = cTF * SUM(rho * rhot**(third/CAT_beta))

  ! Temporary storage
  FFTtmp = FFT(rhoX)
  rhot = FFT(FFTtmp*keKernel(:,:,:,2))

  ! Calculate potential (part 2)
  IF (vacCutoff) THEN
    potential = potential &
              + (CAT_beta*rho**(CAT_beta-1)*cutf + rho**CAT_beta*cutfp) &
                * (FFT(FFTtmp*keKernel(:,:,:,1) &
                     + FFT(rhoX*(rho-rhoS))*keKernel(:,:,:,2)) &
                  + (rho-rhoS)*rhot)
    potential = potential + rho**CAT_beta*cutf*rhot
  ELSE
    potential = potential &
              + CAT_beta*rho**(CAT_beta-1) * &
                  (FFT(FFTtmp*keKernel(:,:,:,1) &
                     + FFT(rhoX*(rho-rhoS))*keKernel(:,:,:,2)) &
                  + (rho-rhoS+rho**CAT_beta)*rhot)
    potential = potential + rho**CAT_beta * rhot
  END IF

  potential = potential*cTF

  RETURN
  
END SUBROUTINE CATPotentialPlus


FUNCTION CATrhot(rho) 
!------------------------------------------------------------------------------
! DESCRIPTION:
! This is CAT KEDF defined in the equation (24) in 
! David Garcia-Aldea and J.E. Alvarellos's  paper (PRB)2007
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!------------------------------------------------------------------------------

  IMPLICIT NONE
  REAL(KIND=DP) :: rho(:,:,:)
  REAL(KIND=DP) :: CATrhot(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) ! the rho_tilde
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FFTTmp
  
  ! has been optimized to reduce the # of FFT calls
  FFTtmp = FFT(rho**CAT_beta)
  CATrhot = FFT(FFTtmp*keKernel(:,:,:,1)+FFT(rho**CAT_beta*(rho-rhoS))*keKernel(:,:,:,2)) &
                    + (rho-rhoS)*FFT(FFTtmp*keKernel(:,:,:,2))

  IF (WGCT==-3) THEN
    CATrhot = CATrhot + & 
      (rho-rhoS) * FFT( FFT(rho**cat_beta*(rho-rhoS))*keKernel(:,:,:,4) )
  ENDIF

  RETURN

END FUNCTION CATrhot


FUNCTION CATvacEnergy(rho)
!------------------------------------------------------------------------------
! DESCRIPTION:
! This is a modified CAT KEDF which works for vacuum
! modified from the original CAT KEDF defined in the equation (24) in 
! David Garcia-Aldea and J.E. Alvarellos's  paper (PRB)2007
!
! The kernel is expaned to the 1st order.
!
! Because beta < 1, (you cannot set beta >1 that will give you worse 
! bulk properties for Aluminum), so the KEDF potential divergence
! when density < 1e-5 in vacuum region
!
! So the rho^beta in CAT KEDF is now multiplied by a cutoff function:
!   CutFunc(r) = (1+Exp[-n0/D])(1/1+Exp[-(rho-n0)/D])-1/(1+Exp[n0/D])
!
!   CutFunc(r->inf) --> 1
!   CutFUnc(r->0  ) ~ r + o(r^2)
!
!   n0 and D are two parameters to be fixed, n0 is called the 
!   vacuum density default value is 1e-5
!   D = 1e-6,
! 
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP) :: CATvacEnergy
  REAL(KIND=DP) :: rho(:,:,:)

                     !>> INTERNAL VARIABLES <<!
  REAL(KIND=DP) :: rhot(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3))
  real(KIND=dp) :: third = 1._DP/3._DP
  
 
  !!!! FORTRAN: it is illegal to raise a negative number with real value power ....
  !!!! so we squre rhot first to make it positive and then take the 1/(3*beta) exponent

  rhot = CATvacrhot(rho)
  rhot = rhot**2
  CATvacEnergy = SUM(rho * rhot**(third/CAT_beta))
  CATvacEnergy = CATvacEnergy*cTF

  RETURN
END FUNCTION CATvacEnergy


FUNCTION CATvacPot(rho) 
!------------------------------------------------------------------------------
! DESCRIPTION:
! This is CAT KEDF defined in the equation (24) in 
! David Garcia-Aldea and J.E. Alvarellos's  paper (PRB)2007
! The kernel is expaned to the 1st order
!
! Because beta < 1, (you cannot set beta >1 that will give you worse 
! bulk properties for Aluminum), so the KEDF potential divergence
! when density < 1e-5 in vacuum region
!
! So the rho^beta in CAT KEDF is now multiplied by a cutoff function:
!   CutFunc(r) = (1+Exp[-n0/D])(1/1+Exp[-(rho-n0)/D])-1/(1+Exp[n0/D])
!
!   CutFunc(r->inf) --> 1
!   CutFUnc(r->0  ) ~ r + o(r^2)
!
!   n0 and D are two parameters to be fixed, n0 is called the 
!   vacuum density default value is 1e-5
!   D = 1e-6,
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!------------------------------------------------------------------------------

  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP) :: rho(:,:,:)
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: CATvacPot 
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: cutf 
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: cutfp 

                     !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: rhoA(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3))
  REAL(KIND=Dp) :: rhot(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3))
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FFTtmp
  REAL(KIND=DP) :: third = 1._dp/3._dp

  !!!!!!FORTRAN: it is illegal to raise a negative number with real value power ....

  cutf = CutoffFunc(rho)
  cutfp = CutoffFuncD(rho)
  
  rhot = CATvacRhot(rho)
  !0:
  CATvacPot = (rhot**2)**(third/CAT_beta)

  !1: 
  rhoA = rho*(rhot**2)**(third/CAT_beta)/rhot*2._DP/3._DP/CAT_beta
    
  FFTtmp = FFT(rhoA)
  CATvacPot = CATvacPot + & 
       (CAT_beta*rho**(CAT_beta-1)*cutf+rho**CAT_beta*cutfp) & 
           * (FFT(kekernel(:,:,:,1)*FFTtmp+kekernel(:,:,:,2)*FFT((rho-rhoS)*rhoA)) & 
          +(rho-rhoS)*FFT(kekernel(:,:,:,2)*FFTtmp))

  CATvacPot = CATvacPot & 
     + rhoA*FFT(FFT(rho**CAT_beta*cutf)*keKernel(:,:,:,2)) & 
     + rho**CAT_beta*cutf*FFT(FFTtmp*keKernel(:,:,:,2))  


  CATvacPot = CATvacPot*cTF
  
  RETURN
END FUNCTION CATvacPot


FUNCTION CATvacRhot(rho) 
!------------------------------------------------------------------------------
! DESCRIPTION:
! This is CAT KEDF defined in the equation (24) in 
! David Garcia-Aldea and J.E. Alvarellos's  paper (PRB)2007
!
! Because beta < 1, (you cannot set beta >1 that will give you worse 
! bulk properties for Aluminum), so the KEDF potential divergence
! when density < 1e-5 in vacuum region
!
! So the rho^beta in CAT KEDF is now multiplied by a cutoff function:
!   CutFunc(r) = (1+Exp[-n0/D])(1/1+Exp[-(rho-n0)/D])-1/(1+Exp[n0/D])
!
!   CutFunc(r->inf) --> 1
!   CutFUnc(r->0  ) ~ r + o(r^2)
!
!   n0 and D are two parameters to be fixed, n0 is called the 
!   vacuum density default value is 1e-5, and D = 1e-6
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   Aug/24/2008   Created (Chen Huang)
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=DP) :: rho(:,:,:)
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: CATVacRhot
  REAL(KIND=DP), DIMENSION(SIZE(rho,1),SIZE(rho,2),SIZE(rho,3)) :: cutff
  COMPLEX(KIND=DP), DIMENSION(k1G, k2G, k3G) :: FFTTmp
  
  ! has been optimized to reduce the # of FFT calls
  cutff =  CutoffFunc(rho)
  FFTtmp = FFT(rho**CAT_beta*cutff)
  CATvacRhot = FFT(   FFTtmp*keKernel(:,:,:,1) & 
                    + FFT(rho**CAT_beta*(rho-rhoS)*cutff)*keKernel(:,:,:,2)) &
            + (rho-rhoS)*FFT(FFTtmp*keKernel(:,:,:,2))

  RETURN

END FUNCTION CATvacRhot


SUBROUTINE FillCAT
!------------------------------------------------------------------------------
! DESCRIPTION: 
! prepare the kinetic energy kernel of CAT 
! (David Garcia-Aldea and J.E. Alvarellos's  paper (PRB)2007)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE KEDF_TF, ONLY: lambda
  !
  USE KEDF_VW, ONLY: mu
  !
  USE IntKernelODE, ONLY: KEDFType
  ! this module used rksuite to integrate kernel's ODE equ.a
  !
  USE IntKernelODE, ONLY: winf
  ! kernel at infinit, starting value for RK method
  !
  USE IntKernelODE, ONLY: int_eta => eta 
  ! kernel is defined on eta
  !
  USE IntKernelODE, ONLY: int_w  => w
  ! the kernel
  !
  USE IntKernelODE, ONLY: int_w1 => w1
  ! d kernel/d eta
  !
  USE IntKernelODE, ONLY: int_w2 => w2
  ! d^2 kernel /d eta^2
  !
  USE IntKernelODE, ONLY: ode_alpha
  ! the alpha in David's paper
  !
  USE IntKernelODE, ONLY: ode_beta
  ! the beta
  !
  USE IntKernelODE, ONLY: ode_gamma
  ! the gamma
  !
  USE IntKernelODE, ONLY: makeKernel2
  ! the driver to do ODE integration
  !
  USE IntKernelODE, ONLY: IntKernelODEClean => Clean 
  ! clean the all the temperory vars

  USE MathSplines, ONLY: spline_cubic_set
  USE MathSplines, ONLY: spline_cubic_val

  USE Constants, ONLY: PI
    
  IMPLICIT NONE
  
  REAL(KIND=DP) :: ypval, yppval
  REAL(KIND=DP) :: eta, tkFStar
  INTEGER :: ix,i2,i3
  REAL(KIND=DP), ALLOCATABLE :: nls_wpp(:), nls_w1pp(:), nls_w2pp(:)
  REAL(KIND=DP) :: zero = 0.0

!! DEBUG
 !   OPEN (unit=1981, file='kernelw.out', status = 'unknown', form='formatted', action='write')
 !   OPEN (unit=1982, file='kernelwp.out', status = 'unknown', form='formatted', action='write')
 !   OPEN (unit=1983, file='kernelwpp.out', status = 'unknown', form='formatted', action='write')
!! ENDDEBUG
    
!   KEDFType = 1 ! for CAT KEDF
!   CAT_Beta  = 0.6666666666666d0
!   CAT_alpha = 0.0
!   CAT_gamma = 1.4d0
!   ode_beta =  CAT_Beta
!   ode_alpha = CAT_alpha
!   ode_gamma = CAT_gamma
!   lambda    = 1.d0  ! coeff for TF
!   mu        = 1.d0        ! coeff for vW
!   wInf      = -0.3333333333d0

   KEDFType = 1 ! for CAT KEDF
   CAT_Beta  = 1.d0
   CAT_alpha = 0.6d0
   !CAT_gamma = 1.4d0
   ode_beta =  CAT_Beta
   ode_alpha = CAT_alpha
   ode_gamma = CAT_gamma
   lambda    = 1.d0        ! coeff for TF, just set to 1 , we take care of it automatically
   mu        = 1.d0        ! coeff for vW, same as for lambda
   wInf      = 0.d0

!   a = 2._DP*(2._DP-3._DP*cat_beta)
!   b = 12._DP
!   c = 6._DP*(cat_beta-1)-10._DP/(1+cat_alpha)*(-0.6_DP+cat_alpha)
!   wInf = (-b+sqrt(b**2-4._DP*a*c))/(2._DP*a) 

! I have found that cat_beta = 2/3, cat_alpha=0 will make w_inf to be single-value=>
! w_inf = -1/3, the only thing you can tune is cat_gamma now
! and the optimum cat_gamma for alumimum diamod and fcc is 1.4
! then I tested the surface energies of aluminum fcc(111/100/110)
! = 943/1156/1177 (mj/m^2), nice!
! in fact, in David Garcia_aldea and JE Alvarellos paper, PRA, (2007)
! they also found that beta =2/3 is the best value, and when beta=2/3
! the CAT becomes same and the WGC KEDF.
    
  CALL WrtOut(6,'')
  CALL WrtOut(6,' Using RK method to solve CAT kernel differential equation...')
  CALL makeKernel2
  CALL WrtOut(6,' Finish solving CAT kernel differential equation.')
  CALL WrtOut(6,'')

  int_w1 = int_eta*int_w1        ! not as before, makeKernel will give int_w1 => w'
  int_w2 = int_eta**2 * int_w2   !                int_w2 =>  w''
  ALLOCATE(nls_wpp(size(int_eta)), nls_w1pp(size(int_eta)), nls_w2pp(size(int_eta)))
    
  CALL spline_cubic_set ( size(int_eta), int_eta, int_w , 0,zero,0,zero, nls_wpp )
  CALL spline_cubic_set ( size(int_eta), int_eta, int_w1, 0,zero,0,zero, nls_w1pp)
  CALL spline_cubic_set ( size(int_eta), int_eta, int_w2, 0,zero,0,zero, nls_w2pp)

  flush(6)

  ! And now we can compute the kernel and its derivatives at every q-point.
  tkFstar = 2._DP * (3._DP * pi**2 * rhoS)**(1._DP/3._DP)
  keKernel = 0._DP
  DO i3=1, SIZE(keKernel,3)
    DO i2=1, SIZE(keKernel,2)
      DO ix=1, SIZE(keKernel,1)
        eta = qTable(ix,i2,i3) / tkFstar

!! DEBUG
          ! Get the kernel value, 1st and 2nd derivitives as a function of
          ! eta.
          ! Temporarily, here keKernel(:,:,:,1) <= w(eta), 
          !                   keKernel(:,:,:,2) <= eta * w'(eta)
          !                   keKernel(:,:,:,3) <= eta**2 * w'(eta)
!! ENDDEBUG

        CALL spline_cubic_val ( size(int_eta), int_eta, int_w,   nls_wpp, eta, kekernel(ix,i2,i3,1), ypval, yppval )
        CALL spline_cubic_val ( size(int_eta), int_eta, int_w1, nls_w1pp, eta, kekernel(ix,i2,i3,2), ypval, yppval )
        CALL spline_cubic_val ( size(int_eta), int_eta, int_w2, nls_w2pp, eta, kekernel(ix,i2,i3,3), ypval, yppval )
!! DEBUG
!
! file operation only work for sequential code...
!          write(1981, *)eta,kekernel(ix,i2,i3,1)
!          write(1982, *)eta,kekernel(ix,i2,i3,2)
!          write(1983, *)eta,kekernel(ix,i2,i3,3)
!! ENDDEBUG

        keKernel(ix,i2,i3,4) = (keKernel(ix,i2,i3,3) &
          + (1+cat_gamma)*keKernel(ix,i2,i3,2)) &
          / (36._DP * rhoS**2)
        keKernel(ix,i2,i3,3) = (keKernel(ix,i2,i3,3) &
          + (7-cat_gamma)*keKernel(ix,i2,i3,2)) &
          / (36._DP * rhoS**2)
        keKernel(ix,i2,i3,2) = -keKernel(ix,i2,i3,2) / (6._DP * rhoS)
      END DO ! ix
    END DO ! i2
  ENDDO ! i3
    
!! DEBUG
!    CLOSE(1981)
!    CLOSE(1982)
!    CLOSE(1983)
!! ENDDEBUG

  CALL IntKernelODEClean
  DEALLOCATE(nls_wpp, nls_w1pp, nls_w2pp)

  RETURN
 
END SUBROUTINE FillCAT


END MODULE KEDF_CAT

MODULE IntKernelODE
!------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE IntKernelODE
!     |_SUBROUTINE makeKernel_ode1 
!     |_SUBROUTINE makeKernel2
!     |_SUBROUTINE WGC
!     |_SUBROUTINE CAT
!     |_SUBROUTINE SingleDep
! 
! DESCRIPTION:
! This module takes the kernel's ODE and do the RK integration with
! the rkSUITE package, the one named rksuite.f in this folder
!
! The makeKernel2() subroutine take a 2nd order ODE to integrate
!                  if your new kernel ode is 1st order eqn. you need to 
!                  write your own  makeKernel subroutine, let's name it
!                  makeKernel1(). makeKernel2() will integrate the 2nd
!                  order ode from infinity to zero. It will obtain w(eta)
!                  w(eta)', w(eta)'' at a sequence of eta points.
!
! Do not modify makeKernel2() unless you really find something wrong
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
! 
!------------------------------------------------------------------------------
! REVISION LOG:
! 2008-02 created by Chen Huang
!-----------------------------------------------------------------------

  USE Constants, ONLY: DP
  USE Output, ONLY: Wrtout 

  REAL(KIND=DP) :: ode_alpha  
  ! this is just the parameter alpha defined in Alvarellos PRB 2007
  !
  REAL(KIND=DP) :: ode_beta   
  ! this is just the parameter beta defined in Alvarellos PRB 2007
  !
  REAL(KIND=DP) :: ode_gamma  
  ! the gamma in the above paper
  !
  REAL(KIND=DP) :: wInf
  !
  INTEGER :: KEDFtype
  !
  REAL(KIND=DP), ALLOCATABLE ::  eta(:), w(:), w1(:), w2(:), wp(:)

CONTAINS


SUBROUTINE makeKernel_ode1
!------------------------------------------------------------------------------
! DESCRIPTION:
! This subroutine calls the user-defined kernel ODE, 
! ode1 => means this is subroutine for 1st order ode problem.
!
! It calls the rksuite.f code to solve the kernels ODE.
!
! The ODE is solved by integrating the kernel ODE from +inf to 0
! If you have a new Kernel ODE, please just plug off the SingleDep() routine
! and plug in your own defined ODE function. 
!
! It has been shown that integration from zero to +inf to be NOT 
! stable.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! Feb 2008 by Chen Huang
!------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: NEQ, LENWRK, METHOD
  !
  PARAMETER (NEQ=1,LENWRK=32*NEQ,METHOD=2)
  !
  REAL(KIND=DP) :: ZERO, ONE, TWO, FOUR
  !
  PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0)
  !
  REAL(KIND=DP) :: HNEXT, HSTART, T, TEND,  &
                   TOL, TSTART, TWANT, WASTE
  !
  INTEGER :: L, STPCST, STPSOK, TOTF, UFLAG
  !
  LOGICAL :: ERRASS, MESAGE
  !
  REAL(KIND=DP):: THRES(NEQ), WORK(LENWRK), Y(NEQ), YMAX(NEQ), &
                  YP(NEQ), YSTART(NEQ)
  !
  EXTERNAL :: SETUP, STAT, UT
  !
  INTEGER :: i, t_num 
  !
  REAL(KIND=DP) :: t_start = 0.E0
  !
  REAL(KIND=DP) :: t_end = 50.E0
  !
  REAL(KIND=DP) :: t_step = 1.E-3
  !
  CHARACTER(len=500) :: message
  !

  ! We use a uniform grid
  t_num = FLOOR(t_end/t_step) + 1

  ! make eta(:)               
  ALLOCATE(eta(t_num),w(t_num), wp(t_num))
  eta = -1.d0
  w = -500.d0
  wp = -500.d0

  DO i=1, t_num
    eta(i) = t_start + REAL(i-1,kind=8)*t_step
  ENDDO
    
  WRITE(message,*) '(IntKernelODE): Solving 1st order ode of kernel now ...' 
  CALL WrtOut(6,message)
  WRITE(message,*) '                # of eta points = ', t_num
  CALL WrtOut(6,message)
  WRITE(message,*) '                eta(min) =',t_start
  CALL WrtOut(6,message)
  WRITE(message,*) '                eta(max) =',t_end
  CALL WrtOut(6,message)
  WRITE(message,*) '                eta_step =',t_step

!
!  Set the initial conditions.  Note that TEND is taken well past
!  the last output point, TLAST.  When this is possible, and it
!  usually is, it is good practice.
!
   TSTART = eta(size(eta))
   YSTART(1) = wInf
   TEND = eta(2)*0.9d0  ! eta(1) = 0._DP, so we will not do this point, we start from eta(2)
!
! Initialize output.
!   
  WRITE (message,'(a,Es12.4)') '                 ode_beta = ',ode_beta
  CALL WrtOut(6,message)
  WRITE (message,'(a,I10,a)') '                 METHOD = ', METHOD, ' 2=>RK(4,5), 3=>RK(7,8)'
  CALL WrtOut(6,message)
  WRITE (message,'(a,Es12.4)') '                 w[+inf] = ',YStart(1)
  CALL WrtOut(6,message)
!
!  Set error control parameters.
!
  TOL = 5.0D-5
  DO L = 1, NEQ
    THRES(L) = 1.0D-10
  ENDDO
!
!  Call the setup routine. Because messages are requested, MESAGE = .TRUE.,
!  there is no need later to test values of flags and print out explanations.
!  In this variant no error assessment is done, so ERRASS is set .FALSE..
!  By setting HSTART to zero, the code is told to find a starting (initial)
!  step size automatically .
!
  MESAGE = .FALSE.
  ERRASS = .FALSE.
  HSTART = ZERO
  CALL SETUP_RK(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,'Usual Task', &
                ERRASS,HSTART,WORK,LENWRK,MESAGE)
!
!  Compute answers at NOUT equally spaced output points. It is good
!  practice to code the calculation of TWANT so that the last value
!  is exactly TLAST.
!
  w(size(eta))= YSTART(1)
  wp(size(eta))= 0.d0

  DO L = size(eta)-1, 2, -1
    TWANT = eta(L)
20  CONTINUE 
    CALL UT_RK(SingleDep,TWANT,T,Y,YP,YMAX,WORK,UFLAG)

    !
    ! make a choice based on UFlag
    !
    SELECT CASE (UFLAG)
      CASE (1) ! success
        w(L)  = Y(1)
        wp(L) = YP(1)
!       WRITE (*,*) eta(L), w(L), w1(L), w2(L)
      CASE (2); GOTO 20
      CASE (3); GOTO 20
      CASE (4); WRITE(*,*)'stiff problem, stop!!'; EXIT
      CASE (5); WRITE(*,*)'too demanding for accuracy, stop'; EXIT
      CASE (6); WRITE(*,*)'global error assessement may not be reliable &
                        & beyond current point, I need to stop!'; EXIT
      CASE (911); WRITE(*,*) 'UFlag=911, stop!!'; EXIT
    END SELECT
!
  ENDDO ! end of loop over eta(t_num:2)

  ! make up w,wp at i=1 without doing integration
  w(1) = w(2)
  wp(1) = wp(2)
!
!  YMAX(L) is the largest magnitude computed for the solution component
!  Y(L) in the course of the integration from TSTART to the last T.  It
!  is used to decide whether THRES(L) is reasonable and to select a new
!  THRES(L) if it is not.
!
  IF (L>1) THEN
    WRITE (*,'(A/)') '             YMAX(L) '
    DO L = 1, NEQ
      WRITE (*,'(13X,1PE9.2)')    YMAX(L)
    ENDDO
    STOP
  ENDIF

!
! check if w, w1, w2 are all assigned during ode solving
!
  DO L=1,size(eta)
    IF (w(L)==-500.d0 .or. wp(L)==-500.d0) THEN
      WRITE(*,*) 'error in final values of w.',w(L), wp(L)
      STOP
    ENDIF
  ENDDO
!
!  The subroutine STAT is used to obtain some information about the progress
!  of the integration. TOTF is the total number of calls to F made so far
!  in the integration; it is a machine-independent measure of work.  At present
!  the integration is finished, so the value printed out refers to the overall
!  cost of the integration.
!
  CALL RK_STATISTICS(TOTF,STPCST,WASTE,STPSOK,HNEXT)
  WRITE (message,'(a,I10)') &
   '                 The cost of the integration in evaluations of F is', TOTF
  CALL WrtOut(6,message)
  WRITE (message,'(a,ES12.4)')& 
    '                 (makeKernel_ode1): w(1)=', w(1)
  CALL WrtOut(6,message)

!     open (unit=110,action='write')
!     do L = 1,size(w,1)
!       write (110,*) eta(L),w(L)
!     enddo
!     close(110)
!     stop
!
END SUBROUTINE makeKernel_ode1


SUBROUTINE  makeKernel2
!----------------------------------------------------------------------------
! DESCRIPTION:
! This subroutine calls the user-defined kernel ODE, i.e. NLS_TF() function.
! and calls the rksuite.f code to solve the kernels ODE.
! '2' means this is subroutine for 2nd order ode problem.
! this ODE is solved by integrating the kernel ODE from +inf to 0
! If you have a new Kernel ODE, please just plug off the NLS_TF() and plug in
! your own defined ODE function. Ofcoure, this makeKernel2() will still 
! integrate from +inf to 0. This has been shown to be the only way to 
! integrate this kernel ode. Integration from zero to +inf has been show by 
! many papers to be not stable and cannot be used.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! Feb 2008 by Chen Huang
!------------------------------------------------------------------------------

  USE Output, ONLY: WrtOut

  IMPLICIT NONE

  INTEGER :: NEQ, LENWRK, METHOD
  !
  PARAMETER (NEQ=2,LENWRK=32*NEQ,METHOD=2)
  !
  REAL(KIND=DP) :: ZERO, ONE, TWO, FOUR
  !
  PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0)
  !
  !     .. Local Scalars ..
  REAL(KIND=DP) :: HNEXT, HSTART, T, TEND,  &
                   TOL, TSTART, TWANT, WASTE
  !
  INTEGER :: L, STPCST, STPSOK, TOTF, UFLAG
  !
  LOGICAL :: ERRASS, MESAGE
  !
  !     .. Local Arrays ..
  REAL(KIND=DP) :: THRES(NEQ), WORK(LENWRK), Y(NEQ), YMAX(NEQ), &
                   YP(NEQ), YSTART(NEQ)
  !
  !     .. External Subroutines ..
  EXTERNAL :: SETUP, STAT, UT
  !
  !     .. vars of chen
  REAL(KIND=DP) :: t_start(2), t_end(2), t_step(2), maxT
  !
  INTEGER :: i, npoints(2)
  !
  CHARACTER(len=500) :: message
  !

! initialize each time-span
!      WRITE(message,'('' (FillWGC) Enter makeKernel2'')')
!      CALL WrtOut(6,message)


  maxT = 100.d0

  t_start(1) = 1.d-3; 
  t_end(1) = 3.d0; 
  t_step(1) = 1d-3;
  
  t_start(2) = t_end(1) + t_step(1)
  t_end(2) = maxT
  t_step(2) = t_step(1)*10.d0
      

! make eta(:)
  npoints(:) = FLOOR((t_end(:)-t_start(:))/t_step(:)) - 1
!    CALL WrtOut(6,'    Solving 2nd order ode of kernel with rksuite' )
!    WRITE(message,*) '   Number of eta points : ', SUM(npoints)
!    CALL WrtOut(6,message)
!    WRITE(message,*) '   Max of eta           : ',maxT
!    CALL WrtOut(6,message)
!    WRITE(message,*) '   Min of eta           : ',t_start(1)
!    CALL WrtOut(6,message)
               
  ALLOCATE(eta(SUM(npoints)))
  ALLOCATE(w(size(eta)), w1(size(eta)), w2(size(eta)))
  eta = -1.d0
  w = -500.d0
  w1 = -500.d0
  w2 = -500.d0

  DO i=1, npoints(1)
    eta(i) = t_start(1) + REAL(i-1,kind=8)*t_step(1)
  ENDDO
  DO i=npoints(1)+1, size(eta)
    eta(i) = t_start(2) + REAL(i-npoints(1)-1,kind=8)*t_step(2)
  ENDDO

!
!  Set the initial conditions.  Note that TEND is taken well past
!  the last output point, TLAST.  When this is possible, and it
!  usually is, it is good practice.
!
  TSTART = eta(size(eta))
  YSTART(1) = wInf
  YSTART(2) = 0.0d0
  TEND = eta(1)*0.9d0
!
! Initialize output.
!   
!    WRITE (message,*) '   ode_alpha =',ode_alpha
!    CALL WrtOut(6,message)
!    WRITE (message,*) '   ode_beta  =',ode_beta
!    CALL WrtOut(6,message)
!    WRITE (message,*) '   ode_gamma =',ode_gamma
!    CALL WrtOut(6,message)
!    WRITE (message,*) '   Method               : ', METHOD, ' (2=>RK(4,5), 3=>RK(7,8))'
!    CALL WrtOut(6,message)
!    WRITE (message,*) '   w_inf                : ', YSTART(1)
!    CALL WrtOut(6,message)
!
!  Set error control parameters.
!
  TOL = 5.0D-5
  DO L = 1, NEQ
    THRES(L) = 1.0D-10
  ENDDO
!
!  Call the setup routine. Because messages are requested, MESAGE = .TRUE.,
!  there is no need later to test values of flags and print out explanations.
!  In this variant no error assessment is done, so ERRASS is set .FALSE..
!  By setting HSTART to zero, the code is told to find a starting (initial)
!  step size automatically .
!
  MESAGE = .FALSE.
  ERRASS = .FALSE.
  HSTART = ZERO
  CALL SETUP_RK(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,'Usual Task', &
                ERRASS,HSTART,WORK,LENWRK,MESAGE)
!
!  Compute answers at NOUT equally spaced output points. It is good
!  practice to code the calculation of TWANT so that the last value
!  is exactly TLAST.
!
  w(size(eta))=0.d0
  w1(size(eta))=0.d0
  w2(size(eta))=0.d0

  DO L = size(eta)-1, 1,-1
     TWANT = eta(L)
20 CONTINUE
     SELECT CASE(KEDFtype)
       CASE(1)
         call UT_RK(CAT,TWANT,T,Y,YP,YMAX,WORK,UFLAG)
       CASE(2)
         call UT_RK(WGC,TWANT,T,Y,YP,YMAX,WORK,UFLAG)
     ENDSELECT
     !
     ! make a choice based on UFlag
     !
     SELECT CASE (UFLAG)
       CASE (1) ! success
         w(L) = Y(1)
         w1(L) = Y(2)
         w2(L) = YP(2)
       CASE (2); GOTO 20
       CASE (3); GOTO 20
       CASE (4); CALL WrtOut(6,'stiff problem, stop!!'); EXIT
       CASE (5); CALL WrtOut(6,'too demanding for accuracy, stop'); EXIT
       CASE (6); CALL WrtOut(6,'global error assessement may not be reliable beyond current point, I need to stop!'); EXIT
       CASE (911); CALL WrtOut(6, 'UFlag=911, stop!!'); EXIT
     END SELECT
!
   ENDDO ! end of loop over eta(i)
!
!  YMAX(L) is the largest magnitude computed for the solution component
!  Y(L) in the course of the integration from TSTART to the last T.  It
!  is used to decide whether THRES(L) is reasonable and to select a new
!  THRES(L) if it is not.
!
   IF (L>1) THEN
     WRITE (message,*) '         YMAX(L) '
     CALL WrtOut(6,message)
     DO L = 1, NEQ
       WRITE (message,*)    YMAX(L)
       CALL WrtOut(6,message)
     ENDDO
     STOP
   ENDIF

!
! check if w, w1, w2 are all assigned during ode solving
!
   DO L=1,size(eta)
     IF (w(L)==-500.d0 .or. w1(L)==-500.d0 .or. w2(L)==-500.d0) THEN
       WRITE(message,*) 'error in final values of w.',w(L), w1(L),w2(L)
       CALL WrtOut(6,message)
       STOP
     ENDIF
   ENDDO
!
!  The subroutine STAT is used to obtain some information about the progress
!  of the integration. TOTF is the total number of calls to F made so far
!  in the integration; it is a machine-independent measure of work.  At present
!  the integration is finished, so the value printed out refers to the overall
!  cost of the integration.
!
   CALL RK_STATISTICS(TOTF,STPCST,WASTE,STPSOK,HNEXT)
!     WRITE (message,*) '   Cost of evaluating F : ', TOTF
!     CALL WrtOut(6,message)
!     WRITE (message,*) '   w(1)                 : ', w(1)
!     CALL WrtOut(6,message)

!DEBUG
!     OPEN(1999,file='ode_kernel',status='replace')
!     DO L=1,size(eta)
!       WRITE(1999,*)eta(L), w(L)
!     ENDDO
!     CLOSE(1999)
!     STOP
!ENDDEBUG.......
!
!    WRITE(message,'('' (FillWGC) Leave makeKernel2'')')
!    CALL WrtOut(6,message)
      
END SUBROUTINE makeKernel2


SUBROUTINE WGC(T,Y,YP)
!------------------------------------------------------------------------------
! DESCRIPTION:
! here is your kernel ODE, you need to give y, and y' upon given t
! i.e. here y is just the w(eta) and t is just the eta
! this NLS_TF ODE is the equation (26) in David Garcia-Aldea and 
! J.E. Alvarellos's  paper (PRB)2007
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  !     .. Scalar Arguments ..
  REAL(KIND=DP) :: T, eta
  !     .. Array Arguments ..
  REAL(KIND=DP) :: Y(*), YP(*), FLind, beta, gamma
  !     .. Executable Statements ..
  !     YP(1) =  Y(2)

  eta = T
  beta  = ode_beta
  gamma = ode_gamma

  ! check eta
  IF (eta<0.d0) THEN
    WRITE(*,*) 'in NLS_TF(), eta<0, stop!'
    STOP
  END IF
      

  YP(1) = Y(2)
  
  IF (eta==1.d0) THEN
    FLind = 2.d0
    WRITE(*,*) 'encounter the case eta=1'
  ELSE IF (eta==0.0d0) THEN 
    FLind = 1.d0
    WRITE(*,*) 'encounter the case eta=0'
  ELSE 
    FLind = 1.d0/(1.d0/2.d0+(1.d0-eta**2)/(4.d0*eta)*LOG(ABS((1.d0+eta)/(1.d0-eta))))
  ENDIF

! WGC formula, Alex Wang, Govind Carter, PRB 1999
  YP(2) = 20.d0 * (FLind- 3.d0*eta**2 -1.d0) - &
         ( gamma+1.d0-6.d0*5.d0/3.d0)*eta*Y(2) - &
         36.d0 * (5.d0/3.d0-beta) * beta * Y(1)

  YP(2) = YP(2)/eta**2

  RETURN

END SUBROUTINE WGC


SUBROUTINE CAT(T,Y,YP)
!------------------------------------------------------------------------
! DESCRIPTION: 
! here is your kernel ODE, you need to give y, and y' upon given t
! i.e. here y is just the w(eta) and t is just the eta
! this NLS_TF ODE is the equation (15) in David Garcia-Aldea and 
! J.E. Alvarellos's  paper (PRB)2007
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------

  USE Output, only: WrtOut

  IMPLICIT NONE
    
  REAL(KIND=DP) :: T, eta
  !
  REAL(KIND=DP) :: Y(*), YP(*), FLind, alpha, beta, gamma
  !

  REAL(KIND=DP) :: lhs
  !
  CHARACTER (LEN=500) :: message
  !

  YP(1) = Y(2)
      
  eta = T

  ! check eta
  IF (eta<0.d0) THEN
    WRITE(*,*) 'in CAT(), eta<0, stop!'
    STOP
  END IF


  beta  = ode_beta
  gamma = ode_gamma
  alpha = ode_alpha

  
  IF (eta==1.d0) THEN
    FLind = 2.d0
    WRITE(message,*) 'encounter the case eta=1'
    CALL WrtOut(6,message)
  ELSE IF (eta==0.0d0) THEN 
    FLind = 1.d0
    WRITE(message,*) 'encounter the case eta=0'
    CALL WrtOut(6,message)
  ELSE 
    FLind = 1.d0/(1.d0/2.d0+(1.d0-eta**2)/(4.d0*eta)*LOG(ABS((1.d0+eta)/(1.d0-eta))))
  ENDIF

! David Garcia-Aldea and J.E. Alvarellos's  paper (24) (PRB)2007
      ! left handside of the (24)
  lhs = 10.d0*(Flind-3.d0*eta**2+alpha)/(1.d0+alpha)
  
  YP(2) = lhs - (  12.d0*Y(1) + 2.d0*(2.d0-3.d0*beta)*Y(1)**2 + & 
                   1.d0/36.d0*2.d0/beta*(2.d0/beta-3.d0)*(Y(2)*eta)**2 &
                   -2.d0*(2.d0/3.d0/beta-1.d0)*Y(1)*Y(2)*eta + 6.d0*(beta-1.d0)  &
                   +(-2.d0/beta-2.d0+(1.d0+gamma)/(3.d0*beta))*Y(2)*eta   &
                 )
  YP(2) = YP(2) *(3.d0*beta)/eta**2

  RETURN

END SUBROUTINE CAT


SUBROUTINE SingleDep(T,Y,YP)
!------------------------------------------------------------------------
! DESCRIPTION:
! This is the kernel solver for teh single denesity dependent kernel 
! by Huang and Carter, submitted (2009). Ask professor for this KEDF which 
! is designed for semiconductors.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! nov/10/2009: Chen Huang (Created)
!------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(KIND=DP) :: T, eta
  REAL(KIND=DP) :: Y(*), YP(*), FLind,beta

  eta = T
  beta  = ode_beta
  
  ! check eta < 0 ?
  IF (eta<0.d0) THEN
    WRITE(*,*) 'in NLS_TF(), eta<0, stop!'
    STOP
  END IF

  IF (eta==1.d0) THEN
    FLind = 2.d0
    WRITE(*,*) 'encounter the case eta=1'
  ELSE IF (eta==0.0d0) THEN 
    FLind = 1.d0
    WRITE(*,*) 'encounter the case eta=0'
  ELSE 
    FLind = 1.d0/(1.d0/2.d0+(1.d0-eta**2)/(4.d0*eta)*LOG(ABS((1.d0+eta)/(1.d0-eta))))
  ENDIF

! kenerl depends on r only => single density dependent kernel
!      YP(1) = (FLind - 3.d0*eta**2 - 1.d0)*5.d0/9.d0 - (5.d0-3.d0*beta)*beta*Y(1)/3.d0
!      YP(1) = YP(1) / ( - beta/3.d0*eta)
  YP(1) = (FLind - 3.d0*eta**2 - 1.d0)*5.d0/3.d0 - (5.d0-3.d0*beta)*beta*Y(1)
  YP(1) = YP(1) / (-beta*eta)

  RETURN
END SUBROUTINE SingleDep


SUBROUTINE Clean
!------------------------------------------------------------------------
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

  IF (ALLOCATED(eta)) DEALLOCATE(eta)
  IF (ALLOCATED(w))   DEALLOCATE(w)
  IF (ALLOCATED(w1))  DEALLOCATE(w1)
  IF (ALLOCATED(w2))  DEALLOCATE(w2)

  RETURN

END SUBROUTINE Clean


END MODULE IntKernelODE

MODULE Clock
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Clock
!     |_FUNCTION TimerStart
!     |_FUNCTION TimerStop
!   SUBROUTINE InitClocks
!   SUBROUTINE StartClock
!   SUBROUTINE StopClock 
!   SUBROUTINE PrintClock
!   SUBROUTINE PrintClockWith
!   MODULE Timer
!     |_FUNCTION TimerStart
!     |_FUNCTION TimerStop
!
! DESCRIPTION:
!   This module holds clocks that can be used to get execution time for
!   each subroutine/function.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/05/2012  File Created (Mohan Chen)
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP

  IMPLICIT NONE
                            
  ! Definition for a type STOPWATCH
  TYPE :: stopwatch
    INTEGER(8) :: & 
      tStart, &        ! the starting time
      tStop, &         ! the stopping time
      tResolution      ! the resolution we want on our time
  END TYPE stopwatch

  SAVE
  ! max number of clocks allowed
  INTEGER, PARAMETER :: maxclock = 200
  !
  ! record wall time (wall_time)
  REAL(KIND=DP)     :: wall_time(maxclock)
  CHARACTER(LEN=30) :: clock_label(maxclock)
  INTEGER           :: called(maxclock)
  TYPE(stopwatch)   :: timer(maxclock) 
  INTEGER :: nclock = 0

END MODULE Clock


SUBROUTINE InitClocks()
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
!
!------------------------------------------------------------------------------
  USE CONSTANTS, ONLY: DP
  USE Clock, ONLY: called
  USE Clock, ONLY: maxclock
  USE Clock, ONLY: clock_label
  USE Clock, ONLY: wall_time
  USE Clock, ONLY: timer

  IMPLICIT NONE

  INTEGER :: i

                           !>> FUNCTION BODY <<!
  !
  DO i = 1, maxclock
     !
     called(i)      = 0
     wall_time(i)   = 0.0_DP
     clock_label(i) = ' '
     timer(i)%tStart= 0 
     timer(i)%tStop = 0
     timer(i)%tResolution = 1000000
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE InitClocks


SUBROUTINE StartClock( label )
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
  USE CONSTANTS, ONLY: DP
  USE Clock,    ONLY : nclock, clock_label, maxclock, timer
  !
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  CHARACTER(LEN=*) :: label
                      !>> INTERNAL VARIABLES <<!
  CHARACTER(LEN=30) :: label_
  INTEGER :: i
                      !>> FUNCTION BODY <<!
  !
  ! ... prevent trouble if label is longer than 30 characters
  !
  label_ = TRIM ( label )
  !
  DO i = 1, nclock
    IF ( clock_label(i) == label_ ) THEN
        CALL SYSTEM_CLOCK(timer(i)%tStart, timer(i)%tResolution)
        RETURN
     END IF
  END DO
  !
  ! ... clock not found : add new clock for given label
  !
  IF ( nclock == maxclock ) THEN
     WRITE( *, '("start_clock: Too many clocks! call ignored")' )
  ELSE
     nclock                = nclock + 1
     clock_label(nclock)   = label_
     CALL SYSTEM_CLOCK(timer(i)%tStart, timer(i)%tResolution)
  END IF
  RETURN
END SUBROUTINE StartClock



SUBROUTINE StopClock( label )
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

  Use Constants, ONLY: DP
  USE Clock, ONLY: nclock
  USE Clock, ONLY: wall_time
  USE Clock, ONLY: clock_label
  USE Clock, ONLY: called
  USE Clock, ONLY: timer
  !
  IMPLICIT NONE
                         !>> EXTERNAL VARIABLES <<!
  CHARACTER(LEN=*) :: label
                         !>> INTERNAL VARIABLES <<!
  CHARACTER(LEN=30) :: label_
  INTEGER :: i
                         !>> FUNCTION BODY <<!
  !
  ! ... prevent trouble if label is longer than 30 characters
  !
  label_ = TRIM ( label )

  !if(label=="PROFESS") then
  !write(*,*) " nclock=", nclock
  !endif
  !return

  !
  DO i = 1, nclock
    IF ( clock_label(i) == label_ ) THEN
        !
        ! ... found previously defined clock : check if properly initialised,
        ! ... add elapsed time, increase the counter of calls
        !
        CALL SYSTEM_CLOCK(timer(i)%tStop)
        wall_time(i)  = wall_time(i) + REAL(timer(i)%tStop - timer(i)%tStart, KIND=DP)/ REAL(timer(i)%tResolution) 

        ! for test
!        write(*,*) " timer::tStop=", timer(i)%tStop
!        write(*,*) " timer::tStart=", timer(i)%tStart

        called(i)    = called(i) + 1
        RETURN
     END IF
  END DO
  RETURN
END SUBROUTINE StopClock
 


SUBROUTINE PrintClock( label, outUnit )
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
  USE Clock, ONLY : nclock, clock_label, wall_time, called

  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!
  CHARACTER(LEN=*) :: label
                      !>> INTERNAL VARIABLES <<!
  INTEGER          :: i
  INTEGER, INTENT(IN) :: outUnit
                      !>> FUNCTION BODY <<!



  WRITE(outUnit,*) " "
  WRITE(outUnit,*) "----------------------------- TIME INFORMATION ------------------------------"

  IF ( label == ' ' ) THEN
     !
     DO i = 1, nclock
        !
        IF(called(i)>0) THEN
          WRITE(outUnit,'(4X,A12," : ",F9.2,"s WALL ",F9.2,"% (" , I8 , " calls)" , F9.2 , "s/call" )' ) & 
            clock_label(i),wall_time(i),wall_time(i)/wall_time(1)*100,called(i),wall_time(i)/called(i)
        ENDIF
        !
     END DO
     !
  END IF
  WRITE(outUnit,*) "-----------------------------------------------------------------------------"

  RETURN
  
END SUBROUTINE PrintClock


SUBROUTINE PrintClockWith(label, outUnit)
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
  USE Clock,    ONLY : nclock, clock_label, wall_time, called
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: label
  INTEGER, INTENT(IN) :: outUnit

  INTEGER :: i

  DO i = 1, nclock
  !
    IF(clock_label(i)==label) THEN
      WRITE(outUnit,'(4X,A12," : ",F9.2,"s WALL ",F9.2,"% (" , I8 , " calls)" , F9.2 , "s/call" )' ) & 
      clock_label(i),wall_time(i),wall_time(i)/wall_time(1)*100,called(i),wall_time(i)/called(i)
    ENDIF
  !
  END DO
 
END SUBROUTINE PrintClockWith



MODULE Timer
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Timer
!     |_FUNCTION TimerStart
!     |_FUNCTION TimerStop
!
! DESCRIPTION:
!   This module holds a timer that can be used to get execution time.
!   This exploits the primitive function SYSTEM_CLOCK.  It wraps around every
!   24 hours or so, so if we see a negative time, that's why.  It returns
!   the wall-clock elapsed time, not the CPU time.
!
!   Please use caution when using this timer.  This is NOT meant to be used
!   in fast loops where it will be executed many times in a small period of
!   time.  It WILL add significant overhead this way.
!
!   Rather, we intend this only to be used to time things that take a 
!   a significant amount of time and contain many (useful) calculations in
!   between ... for instance, if we are doing a geometry optimization, it 
!   might be useful to know (when monitoring the calculation) how long each 
!   density minimization in the loop takes, so we can optimize the geometry
!   minization code, and also time our coffee breaks.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!    
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/24/2003  File Created (Greg Ho)
!
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!
  USE CONSTANTS, ONLY: DP                    ! Abbreviation for double precision

  IMPLICIT NONE
                            
  ! Definition for a type STOPWATCH
  TYPE :: stopwatch
  INTEGER(8) :: & 
    timeStart, &        ! the starting time
    timeStop, &         ! the stopping time
    timeResolution      ! the resolution we want on our time
  END TYPE stopwatch

 
CONTAINS

FUNCTION TimerStart()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function returns a new timer
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!    
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/21/2003  File Created (Vincent Ligneres)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  TYPE(stopwatch) :: &
    TimerStart                  ! the new stopwatch

                          !>> INTERNAL VARIABLES <<! 
                            !>> INITIALIZATION <<!   

  CALL SYSTEM_CLOCK(TimerStart%timeStart, TimerStart%timeResolution)
                           ! >> FUNCTION BODY <<!
  

END FUNCTION TimerStart


FUNCTION TimerStop(watch)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function stops the timer and returns the time in seconds.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!    
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/21/2003  File Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  TYPE(stopwatch), INTENT(INOUT) :: &
    watch                          ! the input timer

  REAL(kind=dp) :: &
    TimerStop                      ! stop the timer

                          !>> INTERNAL VARIABLES <<! 
                            !>> INITIALIZATION <<!   
                           ! >> FUNCTION BODY <<!
  CALL SYSTEM_CLOCK(watch%timeStop)

  TimerStop = REAL(watch%timeStop - watch%timeStart, KIND=dp) &
              / REAL(watch%timeResolution)

END FUNCTION TimerStop


END MODULE Timer

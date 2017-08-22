MODULE OutputFiles
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE OutputFiles
!     |_INTERFACE Error
!     |_SUBROUTINE OneError
!     |_SUBROUTINE ManyError
!
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/20/13  File created (Mohan)
!------------------------------------------------------------------------------
                              !<< GLOBAL >>!


  IMPLICIT NONE

  INTEGER, PARAMETER :: outputUnit = 8        ! output files
  INTEGER, PARAMETER :: errorUnit = 9         ! error files

  CHARACTER(LEN=500) :: message

  INTEGER :: outRank = 0

  ! Thits Error output error messages for character msg or 
  ! charcter array msg(:)
  INTERFACE Error
    MODULE PROCEDURE OneError
    MODULE PROCEDURE ManyError
  END INTERFACE

CONTAINS


SUBROUTINE OneError(unitOut, msg)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Output messages to the unit (6 is the screen), this subroutine receives
!   a warning sign and quit the program.
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

   CHARACTER (LEN=*), INTENT(IN) :: msg    
   ! Message to output
   !
   INTEGER, INTENT(IN)  :: unitOut         
   ! Unit number for writing out
   !

                         !>> FUNCTION BODY <<!

   IF (outRank==0) THEN
       WRITE(unitOut,'(" " )')
       WRITE(unitOut,'(" ---------------------------------- ERROR ------------------------------------" )')
       WRITE(unitOut,'(" " )')
       WRITE(unitOut,'(" ",A)') TRIM(msg)
       WRITE(unitOut,'(" " )')
       WRITE(unitOut,'(" -----------------------------------------------------------------------------" )')
       WRITE(unitOut,'(" QUIT PROFESS.")')

       WRITE(errorUnit,'(" " )')
       WRITE(errorUnit,'(" ---------------------------------- ERROR ------------------------------------" )')
       WRITE(errorUnit,'(" " )')
       WRITE(errorUnit,'(" ",A)') TRIM(msg)
       WRITE(errorUnit,'(" " )')
       WRITE(errorUnit,'(" -----------------------------------------------------------------------------" )')
       WRITE(errorUnit,'(" QUIT PROFESS.")')
   ENDIF
   flush(unitOut)

   STOP

END SUBROUTINE OneError


SUBROUTINE ManyError(unitOut, msg, num)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Output messages to the unit (6 is the screen), this subroutine receives
!   a warning sign and quit the program.
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

   CHARACTER (LEN=*), INTENT(IN) :: msg(:)    
   ! Message to output
   !
   INTEGER, INTENT(IN)  :: unitOut            
   ! Unit number for writing out
   !
   INTEGER, OPTIONAL :: num                   
   ! How many messages need to be printed
   !
   INTEGER :: n                               
   ! Number of messages actually printed
   !
   INTEGER :: i                               
   ! Used for cicle
   !

   IF( PRESENT(num) ) then
      n = num
   ELSE
      n = 1
   END IF

                         !>> FUNCTION BODY <<!
   IF (outRank==0) THEN
       WRITE(unitOut,'(" " )')
       WRITE(unitOut,'(" ------------------------------------- ERROR ---------------------------------------" )')
       WRITE(unitOut,'(" " )')

       DO i = 1, n
       WRITE(unitOut,'(" ",I3,A)') i, TRIM(msg(i))
       ENDDO

       WRITE(unitOut,'(" " )')
       WRITE(unitOut,'(" -----------------------------------------------------------------------------------" )')
       WRITE(unitOut,'(" QUIT PROFESS.")')
   ENDIF
   flush(unitOut)

   STOP

END SUBROUTINE ManyError


END MODULE OutputFiles

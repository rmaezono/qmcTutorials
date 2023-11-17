SUBROUTINE fatal(instring)
!------------------------------------------------------------------!
! Code has hit a bit of a problem so print error message and stop. !
!------------------------------------------------------------------!
 CHARACTER(*),INTENT(in) :: instring
 write(*,fmt="(a)") instring
 stop
END SUBROUTINE fatal

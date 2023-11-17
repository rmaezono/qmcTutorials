MODULE run_control
 USE format_utils, ONLY : wordwrap
 IMPLICIT NONE


CONTAINS


 SUBROUTINE errstop(subroutine,error)
!-------------------------------------------------------!
! Write out routine name and error message then stop.   !
!-------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: subroutine,error
 write(6,*)
 write(6,*)'ERROR : '//trim(subroutine)
 call wordwrap(trim(error))
 write(6,*)
 write(6,'(1x,78(''-''))')
 write(6,*)
 stop
 END SUBROUTINE errstop


 SUBROUTINE errstop_master(subroutine,error)
!-------------------------------------------------------!
! Write out routine name and error message then stop.   !
!-------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: subroutine,error
 write(6,*)
 write(6,*)'ERROR : '//trim(subroutine)
 call wordwrap(trim(error))
 write(6,*)
 write(6,'(1x,78(''-''))')
 write(6,*)
 stop
 END SUBROUTINE errstop_master


 SUBROUTINE errwarn(routine,warning)
!--------------------------!
! Print a warning message. !
!--------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: routine,warning
 call wordwrap('Warning: ['//trim(routine)//'] '//trim(warning))
 write(6,*)
 END SUBROUTINE errwarn


 SUBROUTINE qmc_stop
 IMPLICIT NONE
 stop
 END SUBROUTINE qmc_stop


 SUBROUTINE skip(iunit,nskip)
!---------------------------------------!
! Skip records in a free format file.   !
!---------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: iunit,nskip
 INTEGER iskip
 do iskip=1,nskip
  read(iunit,fmt=*,err=1,end=1)
 enddo
 return
1 call errstop('SKIP','Error reading file with the skip routine.')
 END SUBROUTINE skip


END MODULE run_control

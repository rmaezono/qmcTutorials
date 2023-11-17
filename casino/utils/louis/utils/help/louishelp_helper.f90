PROGRAM casinohelp_helper
 USE esdf
 IMPLICIT NONE
 CHARACTER(128) arg1,arg2
 INTEGER ierr

 read(5,*,iostat=ierr)arg1
 if(ierr/=0)then
  write(6,*)'CASINOHELP_HELPER: Cannot read from standard input.'
  stop
 endif
 read(5,*,iostat=ierr)arg2
 if(ierr/=0)then
  call help_system(trim(arg1))
 else
  call help_system(trim(arg1),trim(arg2))
 endif

END PROGRAM casinohelp_helper

MODULE run_control

 USE dsp
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


 SUBROUTINE errstop_quiet(subroutine,error)
!-------------------------------------------------------------------------!
! Write out routine name and error message then stop without whingeing    !
! about typewriters. Often won't work on parallel machines!               !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: subroutine,error
 write(6,1)subroutine,error
1 format(/1x,'ERROR : ',a,/1x,a/)
 write(6,'(1x,78(''-''))')
 stop
 END SUBROUTINE errstop_quiet


 SUBROUTINE all_stop
 IMPLICIT NONE
 stop
 END SUBROUTINE all_stop


 SUBROUTINE wordwrap(text,unit_in,linelength_in)
!-------------------------------------------------------------------------!
! This subroutine prints out the contents of the character string 'text', !
! ensuring that line breaks only occur at space characters.  The output   !
! is written to unit unit_in if this parameter is supplied; otherwise the !
! output is written to unit 6.  The maximum length of each line is given  !
! by linelength_in if this is supplied; otherwise the default line length !
! is 79 characters.                                                       !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in),OPTIONAL :: unit_in,linelength_in
 CHARACTER(*),INTENT(in) :: text
 CHARACTER(260) :: temp
 INTEGER i,unit,lentext,startpoint,stoppoint,lastpos,linelength
 if(present(unit_in))then
  unit=unit_in
 else
  unit=7
 endif ! unit supplied.
 lentext=len(trim(text))
 if(lentext<1)then
  write(unit,*)
  return
 endif ! No text
 if(present(unit_in))then
  unit=unit_in
 else
  unit=6
 endif ! unit supplied.
 if(present(linelength_in))then
  if(linelength_in>=2)then
   linelength=linelength_in
  else
   linelength=2
  endif ! sensible line-length supplied.
 else
  linelength=79
 endif ! linelength present.
 startpoint=1
 i=0
 do
  i=i+1
  stoppoint=startpoint+linelength-1
  if(stoppoint<=lentext)then
   lastpos=index(trim(text(startpoint:stoppoint))," ",.true.)
   if(lastpos>0)stoppoint=startpoint+lastpos-1
  else
   stoppoint=lentext
  endif ! stoppoint <= length of text
  if(i==1)then
! Allow the user to indent the first line, if (s)he wishes.
   temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
   write(unit,*)trim(temp)
  else
   temp=text(startpoint:stoppoint) ! or else pathscale f90 fails to compile
   write(unit,*)trim(adjustl(temp))
  endif ! i=1
  if(stoppoint==lentext)then
   exit
  else
   startpoint=stoppoint+1
  endif ! Finished text?
 enddo ! Loop over lines.
 END SUBROUTINE wordwrap


END MODULE run_control

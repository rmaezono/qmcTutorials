module error_mod

contains

 subroutine error(err_no,i1,i2)
 implicit none
 integer,intent(in) :: err_no,i1,i2

 select case (err_no)
 case (0)
  write(*,*) 'No error.'
 case(1)
  write(*,'(a)')       'extractcsf: max_targ is too small'
  write(*,'(a,1x,i6)') ' should be:',i1
  write(*,'(a,1x,i6)') ' currently:',i2
 case(2)
  write(*,'(a)')       'extractcsf: max_subshell is too small'
  write(*,'(a,1x,i6)') ' should be:',i1
  write(*,'(a,1x,i6)') ' currently:',i2
 case(3)
  write(*,'(a)')       'extractcsf: max_csf is too small'
  write(*,'(a,1x,i6)') ' should be:',i1
  write(*,'(a,1x,i6)') ' currently:',i2
 case(4)
  write(*,'(a)')       'extractcsf: max_wk is too small'
  write(*,'(a,1x,i6)') ' should be:',i1
  write(*,'(a,1x,i6)') ' currently:',i2
 case(5)
  write(*,'(a)')       '  subshell: maxic is too small'
  write(*,'(a,1x,i6)') ' should be:',i1
  write(*,'(a,1x,i6)') ' currently:',i2
 case(6)
  write(*,'(a)')       '  subshell: max_targ is too small'
  write(*,'(a,1x,i6)') ' should be:',i1
  write(*,'(a,1x,i6)') ' currently:',i2
 case(7)
  write(*,'(a)')       '  subshell: max_det is too small'
  write(*,'(a,1x,i6)') ' should be:',i1
  write(*,'(a,1x,i6)') ' currently:',i2
 case(21)
  write(*,'(a)')       'extractcsf: Sum ^2 Ylm should=1'
  write(*,'(a)')       '             it doesnt.       '
 case(22)
  write(*,'(a)')       'extractcsf: Sum ^2 Xlm should=1'
  write(*,'(a)')       '             it doesnt.       '
 end select

 stop
 
 end subroutine error

end module error_mod


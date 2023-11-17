module file_utils
implicit none
private
public open_units,skip
contains
subroutine open_units(io_no,ierr)
use store, only : open_unit
implicit none
integer,intent(out) :: io_no,ierr
ierr=0
do io_no=10,99
if(.not.open_unit(io_no))exit
enddo
open_unit(io_no)=.true.
if(io_no==99)ierr=1
end subroutine open_units
subroutine skip(iunit,nskip)
use run_control, only : errstop
implicit none
integer,intent(in) :: iunit,nskip
integer xyzzyaaaa3
do xyzzyaaaa3=1,nskip
read(iunit,fmt=*,err=1,end=1)
enddo
return
1 call errstop('SKIP','Error reading file with the skip routine.')
end subroutine skip
end module file_utils

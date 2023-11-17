module singleton
implicit none
private
public fft,fftn,fftkind
integer,parameter :: fftkind=kind(0.d0)
real(fftkind),parameter :: xyzzyaaaa1=0.86602540378443865_fftkind,xyzz&
&yaaab1=0.30901699437494742_fftkind,xyzzyaaac1=0.95105651629515357_fft&
&kind,xyzzyaaad1=3.14159265358979323_fftkind
interface fft
module procedure xyzzyaaae1
module procedure xyzzyaaaf1
module procedure xyzzyaaag1
module procedure xyzzyaaah1
module procedure xyzzyaaai1
module procedure xyzzyaaaj1
module procedure xyzzyaaak1
end interface
contains
function xyzzyaaae1(array,dim,inv,stat) result(ft)
implicit none
complex(fftkind),dimension(:),intent(in) :: array
integer,dimension(:),intent(in),optional :: dim
logical,intent(in),optional :: inv
integer,intent(out),optional :: stat
complex(fftkind),dimension(size(array,1)) :: ft
ft=array
call fftn(ft,shape(array),inv=inv,stat=stat)
end function xyzzyaaae1
function xyzzyaaaf1(array,dim,inv,stat) result(ft)
implicit none
complex(fftkind),dimension(:,:),intent(in) :: array
integer,dimension(:),intent(in),optional :: dim
logical,intent(in),optional :: inv
integer,intent(out),optional :: stat
complex(fftkind),dimension(size(array,1),size(array,2)) :: ft
ft=array
call fftn(ft,shape(array),dim,inv,stat)
end function xyzzyaaaf1
function xyzzyaaag1(array,dim,inv,stat) result(ft)
implicit none
complex(fftkind),dimension(:,:,:),intent(in) :: array
integer,dimension(:),intent(in),optional :: dim
logical,intent(in),optional :: inv
integer,intent(out),optional :: stat
complex(fftkind),dimension(size(array,1),size(array,2),size(array,3)):&
&: ft
ft=array
call fftn(ft,shape(array),dim,inv,stat)
end function xyzzyaaag1
function xyzzyaaah1(array,dim,inv,stat) result(ft)
implicit none
complex(fftkind),dimension(:,:,:,:),intent(in) :: array
integer,dimension(:),intent(in),optional :: dim
logical,intent(in),optional :: inv
integer,intent(out),optional :: stat
complex(fftkind),dimension(size(array,1),size(array,2),size(array,3),s&
&ize(array,4)):: ft
ft=array
call fftn(ft,shape(array),dim,inv,stat)
end function xyzzyaaah1
function xyzzyaaai1(array,dim,inv,stat) result(ft)
implicit none
complex(fftkind),dimension(:,:,:,:,:),intent(in) :: array
integer,dimension(:),intent(in),optional :: dim
logical,intent(in),optional :: inv
integer,intent(out),optional :: stat
complex(fftkind),dimension(size(array,1),size(array,2),size(array,3),s&
&ize(array,4),size(array,5)):: ft
ft=array
call fftn(ft,shape(array),dim,inv,stat)
end function xyzzyaaai1
function xyzzyaaaj1(array,dim,inv,stat) result(ft)
implicit none
complex(fftkind),dimension(:,:,:,:,:,:),intent(in) :: array
integer,dimension(:),intent(in),optional :: dim
logical,intent(in),optional :: inv
integer,intent(out),optional :: stat
complex(fftkind),dimension(size(array,1),size(array,2),size(array,3),s&
&ize(array,4),size(array,5),size(array,6)) :: ft
ft=array
call fftn(ft,shape(array),dim,inv,stat)
end function xyzzyaaaj1
function xyzzyaaak1(array,dim,inv,stat) result(ft)
complex(fftkind),dimension(:,:,:,:,:,:,:),intent(in) :: array
integer,dimension(:),intent(in),optional :: dim
logical,intent(in),optional :: inv
integer,intent(out),optional :: stat
complex(fftkind),dimension(size(array,1),size(array,2),size(array,3),s&
&ize(array,4),size(array,5),size(array,6),size(array,7)):: ft
ft=array
call fftn(ft,shape(array),dim,inv,stat)
end function xyzzyaaak1
subroutine fftn(array,shapea,dim,inv,stat)
implicit none
complex(fftkind),dimension(*),intent(inout) :: array
integer,dimension(:),intent(in) :: shapea
integer,dimension(:),intent(in),optional :: dim
logical,intent(in),optional :: inv
integer,intent(out),optional :: stat
integer xyzzyaaaa16(size(shapea))
integer xyzzyaaab16,xyzzyaaac16,xyzzyaaad16
real(fftkind) xyzzyaaae16
logical xyzzyaaaf16
if(present(inv))then
xyzzyaaaf16=inv
else
xyzzyaaaf16=.false.
endif
if(present(dim))then
xyzzyaaac16=min(size(dim),size(xyzzyaaaa16))
xyzzyaaaa16(1:xyzzyaaac16)=dim(1:xyzzyaaac16)
else
xyzzyaaac16=size(xyzzyaaaa16)
do xyzzyaaab16=1,xyzzyaaac16
xyzzyaaaa16(xyzzyaaab16)=xyzzyaaab16
enddo
endif
xyzzyaaad16=product(shapea)
xyzzyaaae16=sqrt(1.0_fftkind/product(shapea(xyzzyaaaa16(1:xyzzyaaac16)&
&)))
array(1:xyzzyaaad16)=array(1:xyzzyaaad16)*xyzzyaaae16
do xyzzyaaab16=1,xyzzyaaac16
call xyzzyaaal1(array,xyzzyaaad16,shapea(xyzzyaaaa16(xyzzyaaab16)),pro&
&duct(shapea(1:xyzzyaaaa16(xyzzyaaab16))),xyzzyaaaf16,stat)
if(present(stat))then
if(stat/=0)return
endif
enddo
end subroutine fftn
subroutine xyzzyaaal1(array,ntotal,npass,nspan,inv,stat)
implicit none
integer,intent(in) :: ntotal,npass,nspan
complex(fftkind),dimension(*),intent(inout) :: array
logical,intent(in) :: inv
integer,intent(out),optional :: stat
integer xyzzyaaaa17(bit_size(0))
complex(fftkind),dimension(:),allocatable :: xyzzyaaab17
real(fftkind),dimension(:),allocatable :: xyzzyaaac17,xyzzyaaad17
integer,dimension(:),allocatable :: xyzzyaaae17
integer xyzzyaaaf17,xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,xyzzyaaaj17,xy&
&zzyaaak17,xyzzyaaal17,xyzzyaaam17,xyzzyaaan17,xyzzyaaao17,xyzzyaaap17&
&,xyzzyaaaq17,xyzzyaaar17,xyzzyaaas17,xyzzyaaat17,xyzzyaaau17,xyzzyaaa&
&v17,xyzzyaaaw17,xyzzyaaax17,xyzzyaaay17
real(fftkind) xyzzyaaaz17,xyzzyaaba17,xyzzyaabb17,xyzzyaabc17,xyzzyaab&
&d17,xyzzyaabe17,xyzzyaabf17,xyzzyaabg17,xyzzyaabh17,xyzzyaabi17,xyzzy&
&aabj17,xyzzyaabk17,xyzzyaabl17,xyzzyaabm17
complex(fftkind) xyzzyaabn17,xyzzyaabo17,xyzzyaabp17,xyzzyaabq17,xyzzy&
&aabr17,xyzzyaabs17,xyzzyaabt17
if(npass<=1)return
xyzzyaaba17=xyzzyaaab1
if(inv)then
xyzzyaabb17=xyzzyaaac1
xyzzyaaaz17=xyzzyaaaa1
xyzzyaabc17=xyzzyaaad1
else
xyzzyaabb17=-xyzzyaaac1
xyzzyaaaz17=-xyzzyaaaa1
xyzzyaabc17=-xyzzyaaad1
endif
xyzzyaaav17=ntotal
xyzzyaaau17=nspan
xyzzyaaag17=xyzzyaaau17
xyzzyaaat17=xyzzyaaav17-1
xyzzyaaaj17=xyzzyaaau17/npass
xyzzyaabd17=xyzzyaabc17*xyzzyaaaj17
xyzzyaabc17=xyzzyaabc17*2.0_fftkind
call xyzzyaabu17
xyzzyaaaw17=maxval(xyzzyaaaa17(1:xyzzyaaax17))
if(xyzzyaaax17-ishft(xyzzyaaas17,1)>0)then
xyzzyaaay17=max(xyzzyaaax17+1,product(xyzzyaaaa17(xyzzyaaas17+1:xyzzya&
&aax17-xyzzyaaas17))-1)
else
xyzzyaaay17=xyzzyaaax17+1
endif
if(present(stat))then
stat=0
allocate(xyzzyaaab17(xyzzyaaaw17),xyzzyaaac17(xyzzyaaaw17),xyzzyaaad17&
&(xyzzyaaaw17),stat=stat)
if(stat/=0)return
call xyzzyaabv17
deallocate(xyzzyaaac17,xyzzyaaad17,stat=stat)
if(stat/=0)return
allocate(xyzzyaaae17(xyzzyaaay17),stat=stat)
if(stat/=0)return
call xyzzyaabw17
deallocate(xyzzyaaae17,xyzzyaaab17,stat=stat)
if(stat/=0)return
else
allocate(xyzzyaaab17(xyzzyaaaw17),xyzzyaaac17(xyzzyaaaw17),xyzzyaaad17&
&(xyzzyaaaw17))
call xyzzyaabv17
deallocate(xyzzyaaac17,xyzzyaaad17)
allocate(xyzzyaaae17(xyzzyaaay17))
call xyzzyaabw17
deallocate(xyzzyaaae17,xyzzyaaab17)
endif
contains
subroutine xyzzyaabu17
implicit none
xyzzyaaax17=0
xyzzyaaam17=npass
do while(mod(xyzzyaaam17,16)==0)
xyzzyaaax17=xyzzyaaax17+1
xyzzyaaaa17(xyzzyaaax17)=4
xyzzyaaam17=xyzzyaaam17/16
enddo
xyzzyaaai17=3
xyzzyaaal17=9
do
do while(mod(xyzzyaaam17,xyzzyaaal17)==0)
xyzzyaaax17=xyzzyaaax17+1
xyzzyaaaa17(xyzzyaaax17)=xyzzyaaai17
xyzzyaaam17=xyzzyaaam17/xyzzyaaal17
enddo
xyzzyaaai17=xyzzyaaai17+2
xyzzyaaal17=xyzzyaaai17*xyzzyaaai17
if(xyzzyaaal17>xyzzyaaam17)exit
enddo
if(xyzzyaaam17<=4)then
xyzzyaaas17=xyzzyaaax17
xyzzyaaaa17(xyzzyaaax17+1)=xyzzyaaam17
if(xyzzyaaam17/=1)xyzzyaaax17=xyzzyaaax17+1
else
if(xyzzyaaam17-ishft(xyzzyaaam17/4,2)==0)then
xyzzyaaax17=xyzzyaaax17+1
xyzzyaaaa17(xyzzyaaax17)=2
xyzzyaaam17=xyzzyaaam17/4
endif
xyzzyaaas17=xyzzyaaax17
xyzzyaaai17=2
do
if(mod(xyzzyaaam17,xyzzyaaai17)==0)then
xyzzyaaax17=xyzzyaaax17+1
xyzzyaaaa17(xyzzyaaax17)=xyzzyaaai17
xyzzyaaam17=xyzzyaaam17/xyzzyaaai17
endif
xyzzyaaai17=ishft((xyzzyaaai17+1)/2,1)+1
if(xyzzyaaai17>xyzzyaaam17)exit
enddo
endif
if(xyzzyaaas17>0)then
xyzzyaaai17=xyzzyaaas17
do
xyzzyaaax17=xyzzyaaax17+1
xyzzyaaaa17(xyzzyaaax17)=xyzzyaaaa17(xyzzyaaai17)
xyzzyaaai17=xyzzyaaai17-1
if(xyzzyaaai17==0)exit
enddo
endif
end subroutine xyzzyaabu17
subroutine xyzzyaabv17
implicit none
xyzzyaaaf17=0
xyzzyaaak17=0
do
xyzzyaabm17=xyzzyaabd17/xyzzyaaag17
xyzzyaabh17=sin(xyzzyaabm17)
xyzzyaabh17=2.0_fftkind*xyzzyaabh17*xyzzyaabh17
xyzzyaabm17=sin(xyzzyaabm17+xyzzyaabm17)
xyzzyaaar17=1
xyzzyaaaf17=xyzzyaaaf17+1
select case(xyzzyaaaa17(xyzzyaaaf17))
case(2)
xyzzyaaag17=xyzzyaaag17/2
xyzzyaaan17=xyzzyaaag17+2
do
do
xyzzyaaao17=xyzzyaaar17+xyzzyaaag17
xyzzyaabp17=array(xyzzyaaao17)
array(xyzzyaaao17)=array(xyzzyaaar17)-xyzzyaabp17
array(xyzzyaaar17)=array(xyzzyaaar17)+xyzzyaabp17
xyzzyaaar17=xyzzyaaao17+xyzzyaaag17
if(xyzzyaaar17>xyzzyaaat17)exit
enddo
xyzzyaaar17=xyzzyaaar17-xyzzyaaat17
if(xyzzyaaar17>xyzzyaaaj17)exit
enddo
if(xyzzyaaar17>xyzzyaaag17)return
do
xyzzyaabe17=1.0_fftkind-xyzzyaabh17
xyzzyaabj17=xyzzyaabm17
do
do
do
xyzzyaaao17=xyzzyaaar17+xyzzyaaag17
xyzzyaabp17=array(xyzzyaaar17)-array(xyzzyaaao17)
array(xyzzyaaar17)=array(xyzzyaaar17)+array(xyzzyaaao17)
array(xyzzyaaao17)=xyzzyaabp17*cmplx(xyzzyaabe17,xyzzyaabj17,kind=fftk&
&ind)
xyzzyaaar17=xyzzyaaao17+xyzzyaaag17
if(xyzzyaaar17>=xyzzyaaav17)exit
enddo
xyzzyaaao17=xyzzyaaar17-xyzzyaaav17
xyzzyaabe17=-xyzzyaabe17
xyzzyaaar17=xyzzyaaan17-xyzzyaaao17
if(xyzzyaaar17<=xyzzyaaao17)exit
enddo
xyzzyaabi17=xyzzyaabe17-(xyzzyaabh17*xyzzyaabe17+xyzzyaabm17*xyzzyaabj&
&17)
xyzzyaabj17=xyzzyaabm17*xyzzyaabe17-xyzzyaabh17*xyzzyaabj17+xyzzyaabj1&
&7
xyzzyaabe17=2.0_fftkind-(xyzzyaabi17*xyzzyaabi17+xyzzyaabj17*xyzzyaabj&
&17)
xyzzyaabj17=xyzzyaabj17*xyzzyaabe17
xyzzyaabe17=xyzzyaabe17*xyzzyaabi17
xyzzyaaar17=xyzzyaaar17+xyzzyaaaj17
if(xyzzyaaar17>=xyzzyaaao17)exit
enddo
xyzzyaaan17=xyzzyaaan17+1+1
xyzzyaaar17=(xyzzyaaan17-xyzzyaaag17)/2+xyzzyaaaj17
if(xyzzyaaar17>xyzzyaaaj17+xyzzyaaaj17)exit
enddo
case(4)
xyzzyaaah17=xyzzyaaag17
xyzzyaaag17=xyzzyaaag17/4
do
xyzzyaabe17=1.0_fftkind
xyzzyaabj17=0.0_fftkind
do
do
xyzzyaaan17=xyzzyaaar17+xyzzyaaag17
xyzzyaaao17=xyzzyaaan17+xyzzyaaag17
xyzzyaaap17=xyzzyaaao17+xyzzyaaag17
xyzzyaabs17=array(xyzzyaaar17)+array(xyzzyaaao17)
xyzzyaabt17=array(xyzzyaaar17)-array(xyzzyaaao17)
xyzzyaabq17=array(xyzzyaaan17)+array(xyzzyaaap17)
xyzzyaabr17=array(xyzzyaaan17)-array(xyzzyaaap17)
array(xyzzyaaar17)=xyzzyaabs17+xyzzyaabq17
xyzzyaabq17=xyzzyaabs17-xyzzyaabq17
if(inv)then
xyzzyaabs17=xyzzyaabt17+cmplx(-aimag(xyzzyaabr17),real(xyzzyaabr17),ki&
&nd=fftkind)
xyzzyaabt17=xyzzyaabt17+cmplx(aimag(xyzzyaabr17),-real(xyzzyaabr17),ki&
&nd=fftkind)
else
xyzzyaabs17=xyzzyaabt17+cmplx(aimag(xyzzyaabr17),-real(xyzzyaabr17),ki&
&nd=fftkind)
xyzzyaabt17=xyzzyaabt17+cmplx(-aimag(xyzzyaabr17),real(xyzzyaabr17),ki&
&nd=fftkind)
endif
if(xyzzyaabj17==0.0_fftkind)then
array(xyzzyaaan17)=xyzzyaabs17
array(xyzzyaaao17)=xyzzyaabq17
array(xyzzyaaap17)=xyzzyaabt17
else
array(xyzzyaaan17)=xyzzyaabs17*cmplx(xyzzyaabe17,xyzzyaabj17,kind=fftk&
&ind)
array(xyzzyaaao17)=xyzzyaabq17*cmplx(xyzzyaabf17,xyzzyaabk17,kind=fftk&
&ind)
array(xyzzyaaap17)=xyzzyaabt17*cmplx(xyzzyaabg17,xyzzyaabl17,kind=fftk&
&ind)
endif
xyzzyaaar17=xyzzyaaap17+xyzzyaaag17
if(xyzzyaaar17>xyzzyaaav17)exit
enddo
xyzzyaabf17=xyzzyaabe17-(xyzzyaabh17*xyzzyaabe17+xyzzyaabm17*xyzzyaabj&
&17)
xyzzyaabj17=xyzzyaabm17*xyzzyaabe17-xyzzyaabh17*xyzzyaabj17+xyzzyaabj1&
&7
xyzzyaabe17=2.0_fftkind-(xyzzyaabf17*xyzzyaabf17+xyzzyaabj17*xyzzyaabj&
&17)
xyzzyaabj17=xyzzyaabj17*xyzzyaabe17
xyzzyaabe17=xyzzyaabe17*xyzzyaabf17
xyzzyaabf17=xyzzyaabe17*xyzzyaabe17-xyzzyaabj17*xyzzyaabj17
xyzzyaabk17=2.0_fftkind*xyzzyaabe17*xyzzyaabj17
xyzzyaabg17=xyzzyaabf17*xyzzyaabe17-xyzzyaabk17*xyzzyaabj17
xyzzyaabl17=xyzzyaabf17*xyzzyaabj17+xyzzyaabk17*xyzzyaabe17
xyzzyaaar17=xyzzyaaar17-xyzzyaaav17+xyzzyaaaj17
if(xyzzyaaar17>xyzzyaaag17)exit
enddo
xyzzyaaar17=xyzzyaaar17-xyzzyaaag17+1
if(xyzzyaaar17>xyzzyaaaj17)exit
enddo
if(xyzzyaaag17==xyzzyaaaj17)return
case default
xyzzyaaam17=xyzzyaaaa17(xyzzyaaaf17)
xyzzyaaah17=xyzzyaaag17
xyzzyaaag17=xyzzyaaag17/xyzzyaaam17
select case(xyzzyaaam17)
case(3)
do
do
xyzzyaaan17=xyzzyaaar17+xyzzyaaag17
xyzzyaaao17=xyzzyaaan17+xyzzyaaag17
xyzzyaabp17=array(xyzzyaaar17)
xyzzyaabo17=array(xyzzyaaan17)+array(xyzzyaaao17)
array(xyzzyaaar17)=xyzzyaabp17+xyzzyaabo17
xyzzyaabp17=xyzzyaabp17-0.5_fftkind*xyzzyaabo17
xyzzyaabo17=(array(xyzzyaaan17)-array(xyzzyaaao17))*xyzzyaaaz17
array(xyzzyaaan17)=xyzzyaabp17+cmplx(-aimag(xyzzyaabo17),real(xyzzyaab&
&o17),kind=fftkind)
array(xyzzyaaao17)=xyzzyaabp17+cmplx(aimag(xyzzyaabo17),-real(xyzzyaab&
&o17),kind=fftkind)
xyzzyaaar17=xyzzyaaao17+xyzzyaaag17
if(xyzzyaaar17>=xyzzyaaat17)exit
enddo
xyzzyaaar17=xyzzyaaar17-xyzzyaaat17
if(xyzzyaaar17>xyzzyaaag17)exit
enddo
case(5)
xyzzyaabf17=xyzzyaaba17*xyzzyaaba17-xyzzyaabb17*xyzzyaabb17
xyzzyaabk17=2.0_fftkind*xyzzyaaba17*xyzzyaabb17
do
do
xyzzyaaan17=xyzzyaaar17+xyzzyaaag17
xyzzyaaao17=xyzzyaaan17+xyzzyaaag17
xyzzyaaap17=xyzzyaaao17+xyzzyaaag17
xyzzyaaaq17=xyzzyaaap17+xyzzyaaag17
xyzzyaabs17=array(xyzzyaaan17)+array(xyzzyaaaq17)
xyzzyaabt17=array(xyzzyaaan17)-array(xyzzyaaaq17)
xyzzyaabq17=array(xyzzyaaao17)+array(xyzzyaaap17)
xyzzyaabr17=array(xyzzyaaao17)-array(xyzzyaaap17)
xyzzyaabn17=array(xyzzyaaar17)
array(xyzzyaaar17)=xyzzyaabn17+xyzzyaabs17+xyzzyaabq17
xyzzyaabp17=xyzzyaabs17*xyzzyaaba17+xyzzyaabq17*xyzzyaabf17+xyzzyaabn1&
&7
xyzzyaabo17=xyzzyaabt17*xyzzyaabb17+xyzzyaabr17*xyzzyaabk17
array(xyzzyaaan17)=xyzzyaabp17+cmplx(-aimag(xyzzyaabo17),real(xyzzyaab&
&o17),kind=fftkind)
array(xyzzyaaaq17)=xyzzyaabp17+cmplx(aimag(xyzzyaabo17),-real(xyzzyaab&
&o17),kind=fftkind)
xyzzyaabp17=xyzzyaabs17*xyzzyaabf17+xyzzyaabq17*xyzzyaaba17+xyzzyaabn1&
&7
xyzzyaabo17=xyzzyaabt17*xyzzyaabk17-xyzzyaabr17*xyzzyaabb17
array(xyzzyaaao17)=xyzzyaabp17+cmplx(-aimag(xyzzyaabo17),real(xyzzyaab&
&o17),kind=fftkind)
array(xyzzyaaap17)=xyzzyaabp17+cmplx(aimag(xyzzyaabo17),-real(xyzzyaab&
&o17),kind=fftkind)
xyzzyaaar17=xyzzyaaaq17+xyzzyaaag17
if(xyzzyaaar17>=xyzzyaaat17)exit
enddo
xyzzyaaar17=xyzzyaaar17-xyzzyaaat17
if(xyzzyaaar17>xyzzyaaag17)exit
enddo
case default
if(xyzzyaaam17/=xyzzyaaak17)then
xyzzyaaak17=xyzzyaaam17
xyzzyaabj17=xyzzyaabc17/xyzzyaaam17
xyzzyaabe17=cos(xyzzyaabj17)
xyzzyaabj17=sin(xyzzyaabj17)
xyzzyaaad17(xyzzyaaak17)=1.0_fftkind
xyzzyaaac17(xyzzyaaak17)=0.0_fftkind
xyzzyaaai17=1
do
xyzzyaaad17(xyzzyaaai17)=xyzzyaaad17(xyzzyaaam17)*xyzzyaabe17+xyzzyaaa&
&c17(xyzzyaaam17)*xyzzyaabj17
xyzzyaaac17(xyzzyaaai17)=xyzzyaaad17(xyzzyaaam17)*xyzzyaabj17-xyzzyaaa&
&c17(xyzzyaaam17)*xyzzyaabe17
xyzzyaaam17=xyzzyaaam17-1
xyzzyaaad17(xyzzyaaam17)=xyzzyaaad17(xyzzyaaai17)
xyzzyaaac17(xyzzyaaam17)=-xyzzyaaac17(xyzzyaaai17)
xyzzyaaai17=xyzzyaaai17+1
if(xyzzyaaai17>=xyzzyaaam17)exit
enddo
endif
do
do
xyzzyaaan17=xyzzyaaar17
xyzzyaaao17=xyzzyaaar17+xyzzyaaah17
xyzzyaabn17=array(xyzzyaaar17)
xyzzyaabp17=xyzzyaabn17
xyzzyaaai17=1
xyzzyaaan17=xyzzyaaan17+xyzzyaaag17
do
xyzzyaaao17=xyzzyaaao17-xyzzyaaag17
xyzzyaaai17=xyzzyaaai17+1
xyzzyaaab17(xyzzyaaai17)=array(xyzzyaaan17)+array(xyzzyaaao17)
xyzzyaabp17=xyzzyaabp17+xyzzyaaab17(xyzzyaaai17)
xyzzyaaai17=xyzzyaaai17+1
xyzzyaaab17(xyzzyaaai17)=array(xyzzyaaan17)-array(xyzzyaaao17)
xyzzyaaan17=xyzzyaaan17+xyzzyaaag17
if(xyzzyaaan17>=xyzzyaaao17)exit
enddo
array(xyzzyaaar17)=xyzzyaabp17
xyzzyaaan17=xyzzyaaar17
xyzzyaaao17=xyzzyaaar17+xyzzyaaah17
xyzzyaaai17=1
do
xyzzyaaan17=xyzzyaaan17+xyzzyaaag17
xyzzyaaao17=xyzzyaaao17-xyzzyaaag17
xyzzyaaal17=xyzzyaaai17
xyzzyaabp17=xyzzyaabn17
xyzzyaabo17=(0.0_fftkind,0.0_fftkind)
xyzzyaaam17=1
do
xyzzyaaam17=xyzzyaaam17+1
xyzzyaabp17=xyzzyaabp17+xyzzyaaab17(xyzzyaaam17)*xyzzyaaad17(xyzzyaaal&
&17)
xyzzyaaam17=xyzzyaaam17+1
xyzzyaabo17=xyzzyaabo17+xyzzyaaab17(xyzzyaaam17)*xyzzyaaac17(xyzzyaaal&
&17)
xyzzyaaal17=xyzzyaaal17+xyzzyaaai17
if(xyzzyaaal17>xyzzyaaak17)xyzzyaaal17=xyzzyaaal17-xyzzyaaak17
if(xyzzyaaam17>=xyzzyaaak17)exit
enddo
xyzzyaaam17=xyzzyaaak17-xyzzyaaai17
array(xyzzyaaan17)=xyzzyaabp17+cmplx(-aimag(xyzzyaabo17),real(xyzzyaab&
&o17),kind=fftkind)
array(xyzzyaaao17)=xyzzyaabp17+cmplx(aimag(xyzzyaabo17),-real(xyzzyaab&
&o17),kind=fftkind)
xyzzyaaai17=xyzzyaaai17+1
if(xyzzyaaai17>=xyzzyaaam17)exit
enddo
xyzzyaaar17=xyzzyaaar17+xyzzyaaah17
if(xyzzyaaar17>xyzzyaaat17)exit
enddo
xyzzyaaar17=xyzzyaaar17-xyzzyaaat17
if(xyzzyaaar17>xyzzyaaag17)exit
enddo
end select
if(xyzzyaaaf17==xyzzyaaax17)return
xyzzyaaar17=xyzzyaaaj17+1
do
xyzzyaabf17=1.0_fftkind-xyzzyaabh17
xyzzyaabj17=xyzzyaabm17
do
xyzzyaabe17=xyzzyaabf17
xyzzyaabk17=xyzzyaabj17
xyzzyaaar17=xyzzyaaar17+xyzzyaaag17
do
do
array(xyzzyaaar17)=cmplx(xyzzyaabf17,xyzzyaabk17,kind=fftkind)*array(x&
&yzzyaaar17)
xyzzyaaar17=xyzzyaaar17+xyzzyaaah17
if(xyzzyaaar17>xyzzyaaav17)exit
enddo
xyzzyaabi17=xyzzyaabj17*xyzzyaabk17
xyzzyaabk17=xyzzyaabj17*xyzzyaabf17+xyzzyaabe17*xyzzyaabk17
xyzzyaabf17=xyzzyaabe17*xyzzyaabf17-xyzzyaabi17
xyzzyaaar17=xyzzyaaar17-xyzzyaaav17+xyzzyaaag17
if(xyzzyaaar17>xyzzyaaah17)exit
enddo
xyzzyaabf17=xyzzyaabe17-(xyzzyaabh17*xyzzyaabe17+xyzzyaabm17*xyzzyaabj&
&17)
xyzzyaabj17=xyzzyaabj17+xyzzyaabm17*xyzzyaabe17-xyzzyaabh17*xyzzyaabj1&
&7
xyzzyaabe17=2.0_fftkind-(xyzzyaabf17*xyzzyaabf17+xyzzyaabj17*xyzzyaabj&
&17)
xyzzyaabj17=xyzzyaabj17*xyzzyaabe17
xyzzyaabf17=xyzzyaabf17*xyzzyaabe17
xyzzyaaar17=xyzzyaaar17-xyzzyaaah17+xyzzyaaaj17
if(xyzzyaaar17>xyzzyaaag17)exit
enddo
xyzzyaaar17=xyzzyaaar17-xyzzyaaag17+xyzzyaaaj17+1
if(xyzzyaaar17>xyzzyaaaj17+xyzzyaaaj17)exit
enddo
end select
enddo
end subroutine xyzzyaabv17
subroutine xyzzyaabw17
implicit none
xyzzyaaae17(1)=xyzzyaaau17
if(xyzzyaaas17>0)then
xyzzyaaam17=xyzzyaaas17+xyzzyaaas17+1
if(xyzzyaaax17<xyzzyaaam17)xyzzyaaam17=xyzzyaaam17-1
xyzzyaaai17=1
xyzzyaaae17(xyzzyaaam17+1)=xyzzyaaaj17
do
xyzzyaaae17(xyzzyaaai17+1)=xyzzyaaae17(xyzzyaaai17)/xyzzyaaaa17(xyzzya&
&aai17)
xyzzyaaae17(xyzzyaaam17)=xyzzyaaae17(xyzzyaaam17+1)*xyzzyaaaa17(xyzzya&
&aai17)
xyzzyaaai17=xyzzyaaai17+1
xyzzyaaam17=xyzzyaaam17-1
if(xyzzyaaai17>=xyzzyaaam17)exit
enddo
xyzzyaaap17=xyzzyaaae17(xyzzyaaam17+1)
xyzzyaaag17=xyzzyaaae17(2)
xyzzyaaar17=xyzzyaaaj17+1
xyzzyaaao17=xyzzyaaag17+1
xyzzyaaai17=1
if(npass/=ntotal)then
permute_multi: do
do
do
xyzzyaaam17=xyzzyaaar17+xyzzyaaaj17
do
xyzzyaabp17=array(xyzzyaaar17)
array(xyzzyaaar17)=array(xyzzyaaao17)
array(xyzzyaaao17)=xyzzyaabp17
xyzzyaaar17=xyzzyaaar17+1
xyzzyaaao17=xyzzyaaao17+1
if(xyzzyaaar17>=xyzzyaaam17)exit
enddo
xyzzyaaar17=xyzzyaaar17+xyzzyaaau17-xyzzyaaaj17
xyzzyaaao17=xyzzyaaao17+xyzzyaaau17-xyzzyaaaj17
if(xyzzyaaar17>=xyzzyaaav17)exit
enddo
xyzzyaaar17=xyzzyaaar17-xyzzyaaav17+xyzzyaaaj17
xyzzyaaao17=xyzzyaaao17-xyzzyaaav17+xyzzyaaag17
if(xyzzyaaao17>=xyzzyaaau17)exit
enddo
do
do
xyzzyaaao17=xyzzyaaao17-xyzzyaaae17(xyzzyaaai17)
xyzzyaaai17=xyzzyaaai17+1
xyzzyaaao17=xyzzyaaae17(xyzzyaaai17+1)+xyzzyaaao17
if(xyzzyaaao17<=xyzzyaaae17(xyzzyaaai17))exit
enddo
xyzzyaaai17=1
do
if(xyzzyaaar17<xyzzyaaao17)cycle permute_multi
xyzzyaaar17=xyzzyaaar17+xyzzyaaaj17
xyzzyaaao17=xyzzyaaao17+xyzzyaaag17
if(xyzzyaaao17>=xyzzyaaau17)exit
enddo
if(xyzzyaaar17>=xyzzyaaau17)exit
enddo
exit
enddo permute_multi
else
permute_single: do
do
xyzzyaabp17=array(xyzzyaaar17)
array(xyzzyaaar17)=array(xyzzyaaao17)
array(xyzzyaaao17)=xyzzyaabp17
xyzzyaaar17=xyzzyaaar17+1
xyzzyaaao17=xyzzyaaao17+xyzzyaaag17
if(xyzzyaaao17>=xyzzyaaau17)exit
enddo
do
do
xyzzyaaao17=xyzzyaaao17-xyzzyaaae17(xyzzyaaai17)
xyzzyaaai17=xyzzyaaai17+1
xyzzyaaao17=xyzzyaaae17(xyzzyaaai17+1)+xyzzyaaao17
if(xyzzyaaao17<=xyzzyaaae17(xyzzyaaai17))exit
enddo
xyzzyaaai17=1
do
if(xyzzyaaar17<xyzzyaaao17)cycle permute_single
xyzzyaaar17=xyzzyaaar17+1
xyzzyaaao17=xyzzyaaao17+xyzzyaaag17
if(xyzzyaaao17>=xyzzyaaau17)exit
enddo
if(xyzzyaaar17>=xyzzyaaau17)exit
enddo
exit
enddo permute_single
endif
xyzzyaaaj17=xyzzyaaap17
endif
if(ishft(xyzzyaaas17,1)+1>=xyzzyaaax17)return
xyzzyaaah17=xyzzyaaae17(xyzzyaaas17+1)
xyzzyaaai17=xyzzyaaax17-xyzzyaaas17
xyzzyaaaa17(xyzzyaaai17+1)=1
do
xyzzyaaaa17(xyzzyaaai17)=xyzzyaaaa17(xyzzyaaai17)*xyzzyaaaa17(xyzzyaaa&
&i17+1)
xyzzyaaai17=xyzzyaaai17-1
if(xyzzyaaai17==xyzzyaaas17)exit
enddo
xyzzyaaas17=xyzzyaaas17+1
xyzzyaaat17=xyzzyaaaa17(xyzzyaaas17)-1
xyzzyaaai17=0
xyzzyaaal17=0
do
xyzzyaaam17=xyzzyaaas17+1
xyzzyaaao17=xyzzyaaaa17(xyzzyaaas17)
xyzzyaaar17=xyzzyaaaa17(xyzzyaaam17)
xyzzyaaai17=xyzzyaaai17+1
if(xyzzyaaai17>xyzzyaaat17)exit
xyzzyaaal17=xyzzyaaal17+xyzzyaaar17
do while(xyzzyaaal17>=xyzzyaaao17)
xyzzyaaal17=xyzzyaaal17-xyzzyaaao17
xyzzyaaao17=xyzzyaaar17
xyzzyaaam17=xyzzyaaam17+1
xyzzyaaar17=xyzzyaaaa17(xyzzyaaam17)
xyzzyaaal17=xyzzyaaal17+xyzzyaaar17
enddo
xyzzyaaae17(xyzzyaaai17)=xyzzyaaal17
enddo
xyzzyaaai17=0
do
do
xyzzyaaai17=xyzzyaaai17+1
xyzzyaaar17=xyzzyaaae17(xyzzyaaai17)
if(xyzzyaaar17>=0)exit
enddo
if(xyzzyaaar17/=xyzzyaaai17)then
do
xyzzyaaam17=xyzzyaaar17
xyzzyaaar17=xyzzyaaae17(xyzzyaaam17)
xyzzyaaae17(xyzzyaaam17)=-xyzzyaaar17
if(xyzzyaaar17==xyzzyaaai17)exit
enddo
xyzzyaaap17=xyzzyaaar17
else
xyzzyaaae17(xyzzyaaai17)=-xyzzyaaai17
if(xyzzyaaai17==xyzzyaaat17)exit
endif
enddo
do
xyzzyaaai17=xyzzyaaap17+1
xyzzyaaav17=xyzzyaaav17-xyzzyaaah17
xyzzyaaaf17=xyzzyaaav17-1+1
if(xyzzyaaav17<0)exit
do
do
xyzzyaaai17=xyzzyaaai17-1
if(xyzzyaaae17(xyzzyaaai17)>=0)exit
enddo
xyzzyaaal17=xyzzyaaaj17
do
xyzzyaaag17=xyzzyaaal17
if(xyzzyaaal17>xyzzyaaaw17)xyzzyaaag17=xyzzyaaaw17
xyzzyaaal17=xyzzyaaal17-xyzzyaaag17
xyzzyaaam17=xyzzyaaae17(xyzzyaaai17)
xyzzyaaar17=xyzzyaaaj17*xyzzyaaam17+xyzzyaaaf17+xyzzyaaal17
xyzzyaaan17=xyzzyaaar17+xyzzyaaag17
xyzzyaaao17=0
do
xyzzyaaao17=xyzzyaaao17+1
xyzzyaaab17(xyzzyaaao17)=array(xyzzyaaan17)
xyzzyaaan17=xyzzyaaan17-1
if(xyzzyaaan17==xyzzyaaar17)exit
enddo
do
xyzzyaaan17=xyzzyaaar17+xyzzyaaag17
xyzzyaaao17=xyzzyaaan17-xyzzyaaaj17*(xyzzyaaam17+xyzzyaaae17(xyzzyaaam&
&17))
xyzzyaaam17=-xyzzyaaae17(xyzzyaaam17)
do
array(xyzzyaaan17)=array(xyzzyaaao17)
xyzzyaaan17=xyzzyaaan17-1
xyzzyaaao17=xyzzyaaao17-1
if(xyzzyaaan17==xyzzyaaar17)exit
enddo
xyzzyaaar17=xyzzyaaao17
if(xyzzyaaam17==xyzzyaaai17)exit
enddo
xyzzyaaan17=xyzzyaaar17+xyzzyaaag17
xyzzyaaao17=0
do
xyzzyaaao17=xyzzyaaao17+1
array(xyzzyaaan17)=xyzzyaaab17(xyzzyaaao17)
xyzzyaaan17=xyzzyaaan17-1
if(xyzzyaaan17==xyzzyaaar17)exit
enddo
if(xyzzyaaal17==0)exit
enddo
if(xyzzyaaai17==1)exit
enddo
enddo
end subroutine xyzzyaabw17
end subroutine xyzzyaaal1
end module singleton

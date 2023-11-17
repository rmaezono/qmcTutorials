module format_utils
implicit none
private
public wout,wordwrap,capitalize,global_time_heading,i2s,i2s64,r2s,r2ss&
&,r2s2,d2s,d2s0,z2s,l2s,log2int,int2log,write_list_int,byte2human,time&
&2human,labelled_list,display_param,switch_case,allow_slave_write,prin&
&t_centred_line,get_field,pad_int
public r2s_length
public arg_separator
integer,parameter :: r2s_length=80
logical :: allow_slave_write=.true.
type arg_separator
logical,pointer :: dummy=>null()
end type arg_separator
interface wout
module procedure xyzzyaaad1,xyzzyaaae1,xyzzyaaaf1,xyzzyaaag1,xyzzyaaah&
&1,xyzzyaaai1,xyzzyaaaj1,xyzzyaaak1,xyzzyaaal1,xyzzyaaam1,xyzzyaaan1,x&
&yzzyaaao1,xyzzyaaap1,xyzzyaaaq1,xyzzyaaar1,xyzzyaaas1,xyzzyaaat1,xyzz&
&yaaau1
end interface wout
contains
logical function xyzzyaaaa1(io)
use store, only : wout_inhibit_node,o
integer,intent(in) :: io
xyzzyaaaa1=allow_slave_write.or..not.wout_inhibit_node.or.io/=o
end function xyzzyaaaa1
logical function xyzzyaaab1(io)
use store, only : o,output_file
integer,intent(in) :: io
logical xyzzyaaaa4
if(io==o.and.o>6)then
inquire(unit=o,opened=xyzzyaaaa4)
if(.not.xyzzyaaaa4)then
open(unit=o,file=trim(output_file),status='unknown',position='append')
xyzzyaaab1=.false.
else
xyzzyaaab1=.true.
endif
else
xyzzyaaab1=.true.
endif
end function xyzzyaaab1
subroutine xyzzyaaac1(io,was_open)
use store, only : o
integer,intent(in) :: io
logical,intent(in) :: was_open
if(.not.was_open.and.io==o.and.o>6)close(io)
end subroutine xyzzyaaac1
subroutine xyzzyaaad1(unused_arg_separator,unit,fmt)
use store, only : o
integer,intent(in),optional :: unit
character(*),intent(in),optional :: fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa6
logical xyzzyaaab6
xyzzyaaaa6=o
if(present(unit))xyzzyaaaa6=unit
if(.not.xyzzyaaaa1(xyzzyaaaa6))return
xyzzyaaab6=xyzzyaaab1(xyzzyaaaa6)
if(present(fmt))then
write(xyzzyaaaa6,fmt)
else
write(xyzzyaaaa6,'(a)')''
endif
call xyzzyaaac1(xyzzyaaaa6,xyzzyaaab6)
end subroutine xyzzyaaad1
subroutine xyzzyaaae1(c,unused_arg_separator,unit,fmt)
use store, only : o
integer,intent(in),optional :: unit
character(*),intent(in) :: c
character(*),intent(in),optional :: fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa7
logical xyzzyaaab7
xyzzyaaaa7=o
if(present(unit))xyzzyaaaa7=unit
if(.not.xyzzyaaaa1(xyzzyaaaa7))return
xyzzyaaab7=xyzzyaaab1(xyzzyaaaa7)
if(present(fmt))then
write(xyzzyaaaa7,fmt)trim(c)
else
write(xyzzyaaaa7,'(1x,a)')trim(c)
endif
call xyzzyaaac1(xyzzyaaaa7,xyzzyaaab7)
end subroutine xyzzyaaae1
subroutine xyzzyaaaf1(c1,c2,unused_arg_separator,unit,fmt)
use store, only : o
integer,intent(in),optional :: unit
character(*),intent(in) :: c1,c2
character(*),intent(in),optional :: fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa8
logical xyzzyaaab8
xyzzyaaaa8=o
if(present(unit))xyzzyaaaa8=unit
if(.not.xyzzyaaaa1(xyzzyaaaa8))return
xyzzyaaab8=xyzzyaaab1(xyzzyaaaa8)
if(present(fmt))then
write(xyzzyaaaa8,fmt)trim(c1),trim(c2)
else
write(xyzzyaaaa8,'(1x,a)')trim(c1)//trim(c2)
endif
call xyzzyaaac1(xyzzyaaaa8,xyzzyaaab8)
end subroutine xyzzyaaaf1
subroutine xyzzyaaag1(c,i,unused_arg_separator,unit,fmt,rfmt,adjust)
use dsp
use store, only : o
integer,intent(in) :: i
integer,intent(in),optional :: unit
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa9
logical xyzzyaaab9,xyzzyaaac9
character(128) tmpr
xyzzyaaaa9=o
if(present(unit))xyzzyaaaa9=unit
if(.not.xyzzyaaaa1(xyzzyaaaa9))return
xyzzyaaab9=xyzzyaaab1(xyzzyaaaa9)
xyzzyaaac9=.false.
if(present(adjust))xyzzyaaac9=adjust
if(present(fmt))then
write(xyzzyaaaa9,fmt)c,i
else
if(present(rfmt))then
write(tmpr,rfmt)i
if(xyzzyaaac9)tmpr=adjustl(tmpr)
else
write(tmpr,*)i
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa9,'(1x,a)')c//trim(tmpr)
endif
call xyzzyaaac1(xyzzyaaaa9,xyzzyaaab9)
end subroutine xyzzyaaag1
subroutine xyzzyaaah1(c1,i,c2,unused_arg_separator,unit,fmt,rfmt,adjus&
&t)
use dsp
use store, only : o
integer,intent(in) :: i
integer,intent(in),optional :: unit
logical,intent(in),optional :: adjust
character(*),intent(in) :: c1,c2
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa10
logical xyzzyaaab10,xyzzyaaac10
character(128) tmpr
xyzzyaaaa10=o
if(present(unit))xyzzyaaaa10=unit
if(.not.xyzzyaaaa1(xyzzyaaaa10))return
xyzzyaaab10=xyzzyaaab1(xyzzyaaaa10)
xyzzyaaac10=.false.
if(present(adjust))xyzzyaaac10=adjust
if(present(fmt))then
write(xyzzyaaaa10,fmt)c1,i,c2
else
if(present(rfmt))then
write(tmpr,rfmt)i
if(xyzzyaaac10)tmpr=adjustl(tmpr)
else
write(tmpr,*)i
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa10,'(1x,a)')c1//trim(tmpr)//trim(c2)
endif
call xyzzyaaac1(xyzzyaaaa10,xyzzyaaab10)
end subroutine xyzzyaaah1
subroutine xyzzyaaai1(c,ii,unused_arg_separator,unit,fmt,rfmt,adjust,r&
&sep)
use dsp
use store, only : o
integer,intent(in) :: ii(:)
integer,intent(in),optional :: unit
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt,rsep
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11
integer,parameter :: xyzzyaaae11=8
logical xyzzyaaaf11,xyzzyaaag11
character(128) tmpr,tmpr2
xyzzyaaad11=o
if(present(unit))xyzzyaaad11=unit
if(.not.xyzzyaaaa1(xyzzyaaad11))return
xyzzyaaaf11=xyzzyaaab1(xyzzyaaad11)
xyzzyaaag11=.false.
if(present(adjust))xyzzyaaag11=adjust
if(present(fmt))then
write(xyzzyaaad11,fmt)c,ii
else
xyzzyaaaa11=size(ii,1)
do xyzzyaaab11=1,xyzzyaaaa11,xyzzyaaae11
if(present(rfmt))then
write(tmpr,rfmt)ii(xyzzyaaab11)
do xyzzyaaac11=1,min(xyzzyaaae11-1,xyzzyaaaa11-xyzzyaaab11)
write(tmpr2,rfmt)ii(xyzzyaaab11+xyzzyaaac11)
if(xyzzyaaag11)tmpr2=adjustl(tmpr2)
if(present(rsep))then
tmpr=trim(tmpr)//rsep//trim(tmpr2)
else
if(xyzzyaaag11)then
tmpr=trim(tmpr)//' '//trim(tmpr2)
else
tmpr=trim(tmpr)//trim(tmpr2)
endif
endif
enddo
else
write(tmpr,*)ii(xyzzyaaab11:min(xyzzyaaaa11,xyzzyaaab11+xyzzyaaae11-1)&
&)
tmpr=adjustl(tmpr)
endif
if(xyzzyaaab11==1)then
write(xyzzyaaad11,'(1x,a)')c//trim(tmpr)
else
write(xyzzyaaad11,'(1x,a)')trim(tmpr)
endif
enddo
endif
call xyzzyaaac1(xyzzyaaad11,xyzzyaaaf11)
end subroutine xyzzyaaai1
subroutine xyzzyaaaj1(c,l,unused_arg_separator,unit,fmt,rfmt,adjust)
use dsp
use store, only : o
integer,intent(in),optional :: unit
logical,intent(in) :: l
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa12
logical xyzzyaaab12,xyzzyaaac12
character(128) tmpr
xyzzyaaaa12=o
if(present(unit))xyzzyaaaa12=unit
if(.not.xyzzyaaaa1(xyzzyaaaa12))return
xyzzyaaab12=xyzzyaaab1(xyzzyaaaa12)
xyzzyaaac12=.false.
if(present(adjust))xyzzyaaac12=adjust
if(present(fmt))then
write(xyzzyaaaa12,fmt)c,l
else
if(present(rfmt))then
write(tmpr,rfmt)l
if(xyzzyaaac12)tmpr=adjustl(tmpr)
else
write(tmpr,*)l
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa12,'(1x,a)')c//trim(tmpr)
endif
call xyzzyaaac1(xyzzyaaaa12,xyzzyaaab12)
end subroutine xyzzyaaaj1
subroutine xyzzyaaak1(c1,l,c2,unused_arg_separator,unit,fmt,rfmt,adjus&
&t)
use dsp
use store, only : o
integer,intent(in),optional :: unit
logical,intent(in) :: l
logical,intent(in),optional :: adjust
character(*),intent(in) :: c1,c2
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa13
logical xyzzyaaab13,xyzzyaaac13
character(128) tmpr
xyzzyaaaa13=o
if(present(unit))xyzzyaaaa13=unit
if(.not.xyzzyaaaa1(xyzzyaaaa13))return
xyzzyaaab13=xyzzyaaab1(xyzzyaaaa13)
xyzzyaaac13=.false.
if(present(adjust))xyzzyaaac13=adjust
if(present(fmt))then
write(xyzzyaaaa13,fmt)c1,l,c2
else
if(present(rfmt))then
write(tmpr,rfmt)l
if(xyzzyaaac13)tmpr=adjustl(tmpr)
else
write(tmpr,*)l
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa13,'(1x,a)')c1//trim(tmpr)//trim(c2)
endif
call xyzzyaaac1(xyzzyaaaa13,xyzzyaaab13)
end subroutine xyzzyaaak1
subroutine xyzzyaaal1(c,ll,unused_arg_separator,unit,fmt,rfmt,adjust,r&
&sep)
use dsp
use store, only : o
integer,intent(in),optional :: unit
logical,intent(in) :: ll(:)
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt,rsep
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14,xyzzyaaad14
integer,parameter :: xyzzyaaae14=20
logical xyzzyaaaf14,xyzzyaaag14
character(128) tmpr,tmpr2
xyzzyaaad14=o
if(present(unit))xyzzyaaad14=unit
if(.not.xyzzyaaaa1(xyzzyaaad14))return
xyzzyaaaf14=xyzzyaaab1(xyzzyaaad14)
xyzzyaaag14=.false.
if(present(adjust))xyzzyaaag14=adjust
if(present(fmt))then
write(xyzzyaaad14,fmt)c,ll
else
xyzzyaaaa14=size(ll,1)
do xyzzyaaab14=1,xyzzyaaaa14,xyzzyaaae14
if(present(rfmt))then
write(tmpr,rfmt)ll(xyzzyaaab14)
do xyzzyaaac14=1,min(xyzzyaaae14-1,xyzzyaaaa14-xyzzyaaab14)
write(tmpr2,rfmt)ll(xyzzyaaab14+xyzzyaaac14)
if(xyzzyaaag14)tmpr2=adjustl(tmpr2)
if(present(rsep))then
tmpr=trim(tmpr)//rsep//trim(tmpr2)
else
if(xyzzyaaag14)then
tmpr=trim(tmpr)//' '//trim(tmpr2)
else
tmpr=trim(tmpr)//trim(tmpr2)
endif
endif
enddo
else
write(tmpr,*)ll(xyzzyaaab14:min(xyzzyaaaa14,xyzzyaaab14+xyzzyaaae14-1)&
&)
tmpr=adjustl(tmpr)
endif
if(xyzzyaaab14==1)then
write(xyzzyaaad14,'(1x,a)')c//trim(tmpr)
else
write(xyzzyaaad14,'(1x,a)')trim(tmpr)
endif
enddo
endif
call xyzzyaaac1(xyzzyaaad14,xyzzyaaaf14)
end subroutine xyzzyaaal1
subroutine xyzzyaaam1(c,r,unused_arg_separator,unit,fmt,rfmt,adjust)
use dsp
use store, only : o
integer,intent(in),optional :: unit
real(sp),intent(in) :: r
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa15
logical xyzzyaaab15,xyzzyaaac15
character(128) tmpr
xyzzyaaaa15=o
if(present(unit))xyzzyaaaa15=unit
if(.not.xyzzyaaaa1(xyzzyaaaa15))return
xyzzyaaab15=xyzzyaaab1(xyzzyaaaa15)
xyzzyaaac15=.false.
if(present(adjust))xyzzyaaac15=adjust
if(present(fmt))then
write(xyzzyaaaa15,fmt)c,r
else
if(present(rfmt))then
write(tmpr,rfmt)r
if(xyzzyaaac15)tmpr=adjustl(tmpr)
else
write(tmpr,*)r
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa15,'(1x,a)')c//trim(tmpr)
endif
call xyzzyaaac1(xyzzyaaaa15,xyzzyaaab15)
end subroutine xyzzyaaam1
subroutine xyzzyaaan1(c1,r,c2,unused_arg_separator,unit,fmt,rfmt,adjus&
&t)
use dsp
use store, only : o
integer,intent(in),optional :: unit
real(sp),intent(in) :: r
logical,intent(in),optional :: adjust
character(*),intent(in) :: c1,c2
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa16
logical xyzzyaaab16,xyzzyaaac16
character(128) tmpr
xyzzyaaaa16=o
if(present(unit))xyzzyaaaa16=unit
if(.not.xyzzyaaaa1(xyzzyaaaa16))return
xyzzyaaab16=xyzzyaaab1(xyzzyaaaa16)
xyzzyaaac16=.false.
if(present(adjust))xyzzyaaac16=adjust
if(present(fmt))then
write(xyzzyaaaa16,fmt)c1,r,c2
else
if(present(rfmt))then
write(tmpr,rfmt)r
if(xyzzyaaac16)tmpr=adjustl(tmpr)
else
write(tmpr,*)r
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa16,'(1x,a)')c1//trim(tmpr)//trim(c2)
endif
call xyzzyaaac1(xyzzyaaaa16,xyzzyaaab16)
end subroutine xyzzyaaan1
subroutine xyzzyaaao1(c,rr,unused_arg_separator,unit,fmt,rfmt,adjust,r&
&sep)
use dsp
use store, only : o
integer,intent(in),optional :: unit
real(sp),intent(in) :: rr(:)
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt,rsep
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17
integer,parameter :: xyzzyaaae17=6
logical xyzzyaaaf17,xyzzyaaag17
character(128) tmpr,tmpr2
xyzzyaaad17=o
if(present(unit))xyzzyaaad17=unit
if(.not.xyzzyaaaa1(xyzzyaaad17))return
xyzzyaaaf17=xyzzyaaab1(xyzzyaaad17)
xyzzyaaag17=.false.
if(present(adjust))xyzzyaaag17=adjust
if(present(fmt))then
write(xyzzyaaad17,fmt)c,rr
else
xyzzyaaaa17=size(rr,1)
do xyzzyaaab17=1,xyzzyaaaa17,xyzzyaaae17
if(present(rfmt))then
write(tmpr,rfmt)rr(xyzzyaaab17)
do xyzzyaaac17=1,min(xyzzyaaae17-1,xyzzyaaaa17-xyzzyaaab17)
write(tmpr2,rfmt)rr(xyzzyaaab17+xyzzyaaac17)
if(xyzzyaaag17)tmpr2=adjustl(tmpr2)
if(present(rsep))then
tmpr=trim(tmpr)//rsep//trim(tmpr2)
else
if(xyzzyaaag17)then
tmpr=trim(tmpr)//' '//trim(tmpr2)
else
tmpr=trim(tmpr)//trim(tmpr2)
endif
endif
enddo
else
write(tmpr,*)rr(xyzzyaaab17:min(xyzzyaaaa17,xyzzyaaab17+xyzzyaaae17-1)&
&)
tmpr=adjustl(tmpr)
endif
if(xyzzyaaab17==1)then
write(xyzzyaaad17,'(1x,a)')c//trim(tmpr)
else
write(xyzzyaaad17,'(1x,a)')trim(tmpr)
endif
enddo
endif
call xyzzyaaac1(xyzzyaaad17,xyzzyaaaf17)
end subroutine xyzzyaaao1
subroutine xyzzyaaap1(c,d,unused_arg_separator,unit,fmt,rfmt,adjust)
use dsp
use store, only : o
integer,intent(in),optional :: unit
real(dp),intent(in) :: d
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa18
logical xyzzyaaab18,xyzzyaaac18
character(128) tmpr
xyzzyaaaa18=o
if(present(unit))xyzzyaaaa18=unit
if(.not.xyzzyaaaa1(xyzzyaaaa18))return
xyzzyaaab18=xyzzyaaab1(xyzzyaaaa18)
xyzzyaaac18=.false.
if(present(adjust))xyzzyaaac18=adjust
if(present(fmt))then
write(xyzzyaaaa18,fmt)c,d
else
if(present(rfmt))then
write(tmpr,rfmt)d
if(xyzzyaaac18)tmpr=adjustl(tmpr)
else
write(tmpr,*)d
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa18,'(1x,a)')c//trim(tmpr)
endif
call xyzzyaaac1(xyzzyaaaa18,xyzzyaaab18)
end subroutine xyzzyaaap1
subroutine xyzzyaaaq1(c1,d,c2,unused_arg_separator,unit,fmt,rfmt,adjus&
&t)
use dsp
use store, only : o
integer,intent(in),optional :: unit
real(dp),intent(in) :: d
logical,intent(in),optional :: adjust
character(*),intent(in) :: c1,c2
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa19
logical xyzzyaaab19,xyzzyaaac19
character(128) tmpr
xyzzyaaaa19=o
if(present(unit))xyzzyaaaa19=unit
if(.not.xyzzyaaaa1(xyzzyaaaa19))return
xyzzyaaab19=xyzzyaaab1(xyzzyaaaa19)
xyzzyaaac19=.false.
if(present(adjust))xyzzyaaac19=adjust
if(present(fmt))then
write(xyzzyaaaa19,fmt)c1,d,c2
else
if(present(rfmt))then
write(tmpr,rfmt)d
if(xyzzyaaac19)tmpr=adjustl(tmpr)
else
write(tmpr,*)d
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa19,'(1x,a)')c1//trim(tmpr)//trim(c2)
endif
call xyzzyaaac1(xyzzyaaaa19,xyzzyaaab19)
end subroutine xyzzyaaaq1
subroutine xyzzyaaar1(c,dd,unused_arg_separator,unit,fmt,rfmt,adjust,r&
&sep)
use dsp
use store, only : o
integer,intent(in),optional :: unit
real(dp),intent(in) :: dd(:)
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt,rsep
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa20,xyzzyaaab20,xyzzyaaac20,xyzzyaaad20
integer,parameter :: xyzzyaaae20=3
logical xyzzyaaaf20,xyzzyaaag20
character(128) tmpr,tmpr2
xyzzyaaad20=o
if(present(unit))xyzzyaaad20=unit
if(.not.xyzzyaaaa1(xyzzyaaad20))return
xyzzyaaaf20=xyzzyaaab1(xyzzyaaad20)
xyzzyaaag20=.false.
if(present(adjust))xyzzyaaag20=adjust
if(present(fmt))then
write(xyzzyaaad20,fmt)c,dd
else
xyzzyaaaa20=size(dd,1)
do xyzzyaaab20=1,xyzzyaaaa20,xyzzyaaae20
if(present(rfmt))then
write(tmpr,rfmt)dd(xyzzyaaab20)
do xyzzyaaac20=1,min(xyzzyaaae20-1,xyzzyaaaa20-xyzzyaaab20)
write(tmpr2,rfmt)dd(xyzzyaaab20+xyzzyaaac20)
if(xyzzyaaag20)tmpr2=adjustl(tmpr2)
if(present(rsep))then
tmpr=trim(tmpr)//rsep//trim(tmpr2)
else
if(xyzzyaaag20)then
tmpr=trim(tmpr)//' '//trim(tmpr2)
else
tmpr=trim(tmpr)//trim(tmpr2)
endif
endif
enddo
else
write(tmpr,*)dd(xyzzyaaab20:min(xyzzyaaaa20,xyzzyaaab20+xyzzyaaae20-1)&
&)
tmpr=adjustl(tmpr)
endif
if(xyzzyaaab20==1)then
write(xyzzyaaad20,'(1x,a)')c//trim(tmpr)
else
write(xyzzyaaad20,'(1x,a)')trim(tmpr)
endif
enddo
endif
call xyzzyaaac1(xyzzyaaad20,xyzzyaaaf20)
end subroutine xyzzyaaar1
subroutine xyzzyaaas1(c,z,unused_arg_separator,unit,fmt,rfmt,adjust)
use dsp
use store, only : o
integer,intent(in),optional :: unit
complex(dp),intent(in) :: z
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa21
logical xyzzyaaab21,xyzzyaaac21
character(128) tmpr
xyzzyaaaa21=o
if(present(unit))xyzzyaaaa21=unit
if(.not.xyzzyaaaa1(xyzzyaaaa21))return
xyzzyaaab21=xyzzyaaab1(xyzzyaaaa21)
xyzzyaaac21=.false.
if(present(adjust))xyzzyaaac21=adjust
if(present(fmt))then
write(xyzzyaaaa21,fmt)c,z
else
if(present(rfmt))then
write(tmpr,rfmt)z
if(xyzzyaaac21)tmpr=adjustl(tmpr)
else
write(tmpr,*)z
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa21,'(1x,a)')c//trim(tmpr)
endif
call xyzzyaaac1(xyzzyaaaa21,xyzzyaaab21)
end subroutine xyzzyaaas1
subroutine xyzzyaaat1(c1,c2,z,unused_arg_separator,unit,fmt,rfmt,adjus&
&t)
use dsp
use store, only : o
integer,intent(in),optional :: unit
complex(dp),intent(in) :: z
logical,intent(in),optional :: adjust
character(*),intent(in) :: c1,c2
character(*),intent(in),optional :: rfmt,fmt
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa22
logical xyzzyaaab22,xyzzyaaac22
character(128) tmpr
xyzzyaaaa22=o
if(present(unit))xyzzyaaaa22=unit
if(.not.xyzzyaaaa1(xyzzyaaaa22))return
xyzzyaaab22=xyzzyaaab1(xyzzyaaaa22)
xyzzyaaac22=.false.
if(present(adjust))xyzzyaaac22=adjust
if(present(fmt))then
write(xyzzyaaaa22,fmt)c1,z,c2
else
if(present(rfmt))then
write(tmpr,rfmt)z
if(xyzzyaaac22)tmpr=adjustl(tmpr)
else
write(tmpr,*)z
tmpr=adjustl(tmpr)
endif
write(xyzzyaaaa22,'(1x,a)')c1//trim(tmpr)//trim(c2)
endif
call xyzzyaaac1(xyzzyaaaa22,xyzzyaaab22)
end subroutine xyzzyaaat1
subroutine xyzzyaaau1(c,zz,unused_arg_separator,unit,fmt,rfmt,adjust,r&
&sep)
use dsp
use store, only : o
integer,intent(in),optional :: unit
complex(dp),intent(in) :: zz(:)
logical,intent(in),optional :: adjust
character(*),intent(in) :: c
character(*),intent(in),optional :: rfmt,fmt,rsep
type(arg_separator),intent(in),optional :: unused_arg_separator
integer xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyzzyaaad23
integer,parameter :: xyzzyaaae23=1
logical xyzzyaaaf23,xyzzyaaag23
character(128) tmpr,tmpr2
xyzzyaaad23=o
if(present(unit))xyzzyaaad23=unit
if(.not.xyzzyaaaa1(xyzzyaaad23))return
xyzzyaaaf23=xyzzyaaab1(xyzzyaaad23)
xyzzyaaag23=.false.
if(present(adjust))xyzzyaaag23=adjust
if(present(fmt))then
write(xyzzyaaad23,fmt)c,zz
else
xyzzyaaaa23=size(zz,1)
do xyzzyaaab23=1,xyzzyaaaa23,xyzzyaaae23
if(present(rfmt))then
write(tmpr,rfmt)zz(xyzzyaaab23)
do xyzzyaaac23=1,min(xyzzyaaae23-1,xyzzyaaaa23-xyzzyaaab23)
write(tmpr2,rfmt)zz(xyzzyaaab23+xyzzyaaac23)
if(xyzzyaaag23)tmpr2=adjustl(tmpr2)
if(present(rsep))then
tmpr=trim(tmpr)//rsep//trim(tmpr2)
else
if(xyzzyaaag23)then
tmpr=trim(tmpr)//' '//trim(tmpr2)
else
tmpr=trim(tmpr)//trim(tmpr2)
endif
endif
enddo
else
write(tmpr,*)zz(xyzzyaaab23:min(xyzzyaaaa23,xyzzyaaab23+xyzzyaaae23-1)&
&)
tmpr=adjustl(tmpr)
endif
if(xyzzyaaab23==1)then
write(xyzzyaaad23,'(1x,a)')c//trim(tmpr)
else
write(xyzzyaaad23,'(1x,a)')trim(tmpr)
endif
enddo
endif
call xyzzyaaac1(xyzzyaaad23,xyzzyaaaf23)
end subroutine xyzzyaaau1
subroutine wordwrap(text,unit_in,linelength_in)
use store, only : o
implicit none
integer,intent(in),optional :: unit_in,linelength_in
character(*),intent(in) :: text
integer xyzzyaaaa24,unit,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaae&
&24,xyzzyaaaf24
character(260) :: xyzzyaaag24
if(present(unit_in))then
unit=unit_in
else
unit=o
endif
xyzzyaaab24=len(trim(text))
if(xyzzyaaab24<1)then
write(unit,*)
return
endif
if(present(linelength_in))then
if(linelength_in>=2)then
xyzzyaaaf24=linelength_in
else
xyzzyaaaf24=2
endif
else
xyzzyaaaf24=79
endif
xyzzyaaac24=1
xyzzyaaaa24=0
do
xyzzyaaaa24=xyzzyaaaa24+1
xyzzyaaad24=xyzzyaaac24+xyzzyaaaf24-1
if(xyzzyaaad24<=xyzzyaaab24)then
xyzzyaaae24=index(trim(text(xyzzyaaac24:xyzzyaaad24))," ",.true.)
if(xyzzyaaae24>0)xyzzyaaad24=xyzzyaaac24+xyzzyaaae24-1
else
xyzzyaaad24=xyzzyaaab24
endif
if(xyzzyaaaa24==1)then
xyzzyaaag24=text(xyzzyaaac24:xyzzyaaad24)
call wout(trim(xyzzyaaag24),unit=unit)
else
xyzzyaaag24=text(xyzzyaaac24:xyzzyaaad24)
call wout(trim(adjustl(xyzzyaaag24)),unit=unit)
endif
if(xyzzyaaad24==xyzzyaaab24)then
exit
else
xyzzyaaac24=xyzzyaaad24+1
endif
enddo
end subroutine wordwrap
subroutine capitalize(string,decapitalize_in)
implicit none
logical,intent(in),optional :: decapitalize_in
character(*),intent(inout) :: string
integer xyzzyaaaa25,xyzzyaaab25
integer,parameter :: xyzzyaaac25=ichar('a'),xyzzyaaad25=ichar('z'),xyz&
&zyaaae25=ichar('A'),xyzzyaaaf25=ichar('Z'),xyzzyaaag25=xyzzyaaae25-xy&
&zzyaaac25
logical xyzzyaaah25
if(present(decapitalize_in))then
xyzzyaaah25=decapitalize_in
else
xyzzyaaah25=.false.
endif
if(xyzzyaaah25)then
do xyzzyaaaa25=1,len_trim(string)
xyzzyaaab25=ichar(string(xyzzyaaaa25:xyzzyaaaa25))
if(xyzzyaaab25>=xyzzyaaae25.and.xyzzyaaab25<=xyzzyaaaf25) string(xyzzy&
&aaaa25:xyzzyaaaa25)=achar(xyzzyaaab25-xyzzyaaag25)
enddo
else
do xyzzyaaaa25=1,len_trim(string)
xyzzyaaab25=ichar(string(xyzzyaaaa25:xyzzyaaaa25))
if(xyzzyaaab25>=xyzzyaaac25.and.xyzzyaaab25<=xyzzyaaad25) string(xyzzy&
&aaaa25:xyzzyaaaa25)=achar(xyzzyaaab25+xyzzyaaag25)
enddo
endif
end subroutine capitalize
function switch_case(string,to_lower) result(switchcase)
implicit none
logical,intent(in),optional :: to_lower
character(*),intent(in) :: string
integer xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26
integer,parameter :: xyzzyaaad26(2)=(/ichar('a'),ichar('A')/),xyzzyaaa&
&e26(2)=(/ichar('z'),ichar('Z')/),                     xyzzyaaaf26(2)=&
&(/ichar('A')-ichar('a'),ichar('a')-ichar('A')/)
character(len(string)) switchcase
switchcase=string
xyzzyaaac26=1
if(present(to_lower))then
if(to_lower)xyzzyaaac26=2
endif
do xyzzyaaaa26=1,len(switchcase)
xyzzyaaab26=ichar(switchcase(xyzzyaaaa26:xyzzyaaaa26))
if(xyzzyaaab26>=xyzzyaaad26(xyzzyaaac26).and.xyzzyaaab26<=xyzzyaaae26(&
&xyzzyaaac26))then
xyzzyaaab26=xyzzyaaab26+xyzzyaaaf26(xyzzyaaac26)
switchcase(xyzzyaaaa26:xyzzyaaaa26)=char(xyzzyaaab26)
endif
enddo
end function switch_case
subroutine global_time_heading(beginning)
implicit none
logical,intent(in) :: beginning
character(20) chdate,chtime
chdate(:)=' '
chtime(:)=' '
call date_and_time(chdate,chtime)
if(beginning)then
call wout('Started '//chdate(1:4)//'/'//chdate(5:6)//'/'//chdate(7:8)/&
&/' '//chtime(1:2)//':'//chtime(3:4)//':'//chtime(5:))
else
call wout('Ends '//chdate(1:4)//'/'//chdate(5:6)//'/'//chdate(7:8)//' &
&'//chtime(1:2)//':'//chtime(3:4)//':'//chtime(5:))
endif
end subroutine global_time_heading
character(20) function i2s(n)
implicit none
integer,intent(in) :: n
integer xyzzyaaaa28,xyzzyaaab28
character tmp,sign
if(n==0)then
i2s='0'
return
endif
sign=' '
if(n<0)sign='-'
do xyzzyaaaa28=1,len(i2s)
i2s(xyzzyaaaa28:xyzzyaaaa28)=' '
enddo
xyzzyaaaa28=abs(n)
do xyzzyaaab28=1,len(i2s)
if(xyzzyaaaa28==0)exit
i2s(xyzzyaaab28:xyzzyaaab28)=achar(ichar('0')+mod(xyzzyaaaa28,10))
xyzzyaaaa28=xyzzyaaaa28/10
enddo
xyzzyaaaa28=1
xyzzyaaab28=len_trim(i2s)
do
if(xyzzyaaaa28>=xyzzyaaab28)exit
tmp=i2s(xyzzyaaab28:xyzzyaaab28)
i2s(xyzzyaaab28:xyzzyaaab28)=i2s(xyzzyaaaa28:xyzzyaaaa28)
i2s(xyzzyaaaa28:xyzzyaaaa28)=tmp
xyzzyaaaa28=xyzzyaaaa28+1
xyzzyaaab28=xyzzyaaab28-1
enddo
i2s=trim(sign)//i2s
end function i2s
character(20) function i2s64(n)
use dsp
integer(i64),intent(in) :: n
integer(i64) xyzzyaaaa29,xyzzyaaab29
character tmp,sign
if(n==0)then
i2s64='0'
return
endif
sign=' '
if(n<0)sign='-'
do xyzzyaaaa29=1,len(i2s64)
i2s64(xyzzyaaaa29:xyzzyaaaa29)=' '
enddo
xyzzyaaaa29=abs(n)
do xyzzyaaab29=1,len(i2s64)
if(xyzzyaaaa29==0)exit
i2s64(xyzzyaaab29:xyzzyaaab29)=achar(ichar('0')+mod(xyzzyaaaa29,10_i64&
&))
xyzzyaaaa29=xyzzyaaaa29/10
enddo
xyzzyaaaa29=1
xyzzyaaab29=len_trim(i2s64)
do
if(xyzzyaaaa29>=xyzzyaaab29)exit
tmp=i2s64(xyzzyaaab29:xyzzyaaab29)
i2s64(xyzzyaaab29:xyzzyaaab29)=i2s64(xyzzyaaaa29:xyzzyaaaa29)
i2s64(xyzzyaaaa29:xyzzyaaaa29)=tmp
xyzzyaaaa29=xyzzyaaaa29+1
xyzzyaaab29=xyzzyaaab29-1
enddo
i2s64=trim(sign)//i2s64
end function i2s64
character(r2s_length) function r2s(r,real_format)
use dsp
implicit none
real(dp),intent(in) :: r
character(*),intent(in) :: real_format
write(r2s,real_format)r
r2s=adjustl(r2s)
end function r2s
character(r2s_length) function r2ss(r,real_format)
use dsp
implicit none
real(sp),intent(in) :: r
character(*),intent(in) :: real_format
write(r2ss,real_format)r
r2ss=adjustl(r2ss)
end function r2ss
character(r2s_length) function r2s2(r,real_format)
use dsp
implicit none
real(dp),intent(in) :: r
character(*),intent(in) :: real_format
write(r2s2,real_format)r
if(r<0.d0)then
r2s2=adjustl(r2s2)
else
r2s2=' '//adjustl(r2s2)
endif
end function r2s2
character(r2s_length) function z2s(z,real_format)
use dsp
implicit none
complex(dp),intent(in) :: z
character(*),intent(in) :: real_format
character(r2s_length) :: xyzzyaaaa33,xyzzyaaab33
write(xyzzyaaaa33,real_format)dble(z)
xyzzyaaaa33=adjustl(xyzzyaaaa33)
write(xyzzyaaab33,real_format)aimag(z)
xyzzyaaab33=adjustl(xyzzyaaab33)
z2s='('//trim(xyzzyaaaa33)//','//trim(xyzzyaaab33)//')'
end function z2s
character(72) function d2s(x,prec_in,exp_char_in,space_for_pos_in)
use dsp
implicit none
integer,intent(in),optional :: prec_in
real(dp),intent(in) :: x
logical,intent(in),optional :: space_for_pos_in
character(1),intent(in),optional :: exp_char_in
integer,parameter :: xyzzyaaaa34=50
integer xyzzyaaab34(xyzzyaaaa34+1),xyzzyaaac34,xyzzyaaad34,xyzzyaaae34&
&,xyzzyaaaf34,xyzzyaaag34
real(dp) xyzzyaaah34,xyzzyaaai34
logical xyzzyaaaj34,xyzzyaaak34
character(1) exp_char,sign_char
if(present(prec_in))then
if(prec_in<=0)then
xyzzyaaac34=precision(1.d0)+1
else
xyzzyaaac34=prec_in
endif
else
xyzzyaaac34=precision(1.d0)+1
endif
if(xyzzyaaac34>xyzzyaaaa34)xyzzyaaac34=xyzzyaaaa34
if(present(exp_char_in))then
exp_char=exp_char_in
else
exp_char='E'
endif
if(present(space_for_pos_in))then
xyzzyaaaj34=space_for_pos_in
else
xyzzyaaaj34=.false.
endif
sign_char=' '
if(x<0)then
sign_char='-'
xyzzyaaak34=.true.
xyzzyaaah34=-x
else
xyzzyaaak34=xyzzyaaaj34
xyzzyaaah34=x
endif
if(xyzzyaaah34/=0.d0)then
xyzzyaaag34=floor(log10(xyzzyaaah34))
else
xyzzyaaag34=0
endif
xyzzyaaae34=xyzzyaaag34
xyzzyaaai34=10.d0**xyzzyaaae34
do xyzzyaaad34=1,xyzzyaaac34
xyzzyaaab34(xyzzyaaad34)=floor(xyzzyaaah34/xyzzyaaai34)
xyzzyaaah34=xyzzyaaah34-dble(xyzzyaaab34(xyzzyaaad34))*xyzzyaaai34
if(xyzzyaaah34<0.d0)xyzzyaaah34=0.d0
xyzzyaaae34=xyzzyaaae34-1
xyzzyaaai34=xyzzyaaai34*0.1d0
enddo
xyzzyaaab34(xyzzyaaac34+1)=floor(xyzzyaaah34/xyzzyaaai34)
if(any(xyzzyaaab34(1:xyzzyaaac34+1)>9).or.any(xyzzyaaab34(1:xyzzyaaac3&
&4+1)<0))then
d2s='NaN'
return
endif
do
if(xyzzyaaab34(xyzzyaaac34+1)>=5)then
xyzzyaaaf34=1
xyzzyaaab34(xyzzyaaac34+1)=0
else
exit
endif
do xyzzyaaad34=xyzzyaaac34,1,-1
xyzzyaaab34(xyzzyaaad34)=xyzzyaaab34(xyzzyaaad34)+xyzzyaaaf34
xyzzyaaaf34=xyzzyaaab34(xyzzyaaad34)/10
xyzzyaaab34(xyzzyaaad34)=mod(xyzzyaaab34(xyzzyaaad34),10)
enddo
if(xyzzyaaaf34>0)then
do xyzzyaaad34=xyzzyaaac34+1,2,-1
xyzzyaaab34(xyzzyaaad34)=xyzzyaaab34(xyzzyaaad34-1)
enddo
xyzzyaaab34(1)=xyzzyaaaf34
xyzzyaaag34=xyzzyaaag34+1
endif
enddo
do xyzzyaaad34=xyzzyaaac34,2,-1
if(xyzzyaaab34(xyzzyaaad34)==0)then
xyzzyaaac34=xyzzyaaac34-1
else
exit
endif
enddo
d2s=char(48+xyzzyaaab34(1))//'.'
if(xyzzyaaak34)d2s=sign_char//trim(d2s)
do xyzzyaaad34=2,xyzzyaaac34
d2s=trim(d2s)//char(48+xyzzyaaab34(xyzzyaaad34))
enddo
if(xyzzyaaag34>0)then
d2s=trim(d2s)//exp_char//'+'//trim(i2s(xyzzyaaag34))
elseif(xyzzyaaag34<0)then
d2s=trim(d2s)//exp_char//trim(i2s(xyzzyaaag34))
endif
end function d2s
character(72) function d2s0(x,prec_in,exp_char_in,space_for_pos_in)
use dsp
implicit none
real(dp),intent(in) :: x
integer,intent(in),optional :: prec_in
character(1),intent(in),optional :: exp_char_in
logical,intent(in),optional :: space_for_pos_in
integer,parameter :: xyzzyaaaa35=50
integer xyzzyaaab35(xyzzyaaaa35+1),xyzzyaaac35,xyzzyaaad35,xyzzyaaae35&
&,xyzzyaaaf35,xyzzyaaag35
real(dp) xyzzyaaah35,xyzzyaaai35
character(1) exp_char,sign_char
logical xyzzyaaaj35,xyzzyaaak35
if(present(prec_in))then
if(prec_in<=0)then
xyzzyaaac35=precision(1.d0)+1
else
xyzzyaaac35=prec_in
endif
else
xyzzyaaac35=precision(1.d0)+1
endif
if(xyzzyaaac35>xyzzyaaaa35)xyzzyaaac35=xyzzyaaaa35
if(present(exp_char_in))then
exp_char=exp_char_in
else
exp_char='E'
endif
if(present(space_for_pos_in))then
xyzzyaaaj35=space_for_pos_in
else
xyzzyaaaj35=.false.
endif
sign_char=' '
if(x<0)then
sign_char='-'
xyzzyaaak35=.true.
xyzzyaaah35=-x
else
xyzzyaaak35=xyzzyaaaj35
xyzzyaaah35=x
endif
if(xyzzyaaah35/=0.d0)then
xyzzyaaag35=floor(log10(xyzzyaaah35))
else
xyzzyaaag35=0
endif
xyzzyaaae35=xyzzyaaag35
xyzzyaaai35=10.d0**xyzzyaaae35
do xyzzyaaad35=1,xyzzyaaac35
xyzzyaaab35(xyzzyaaad35)=floor(xyzzyaaah35/xyzzyaaai35)
xyzzyaaah35=xyzzyaaah35-dble(xyzzyaaab35(xyzzyaaad35))*xyzzyaaai35
if(xyzzyaaah35<0.d0)xyzzyaaah35=0.d0
xyzzyaaae35=xyzzyaaae35-1
xyzzyaaai35=xyzzyaaai35*0.1d0
enddo
xyzzyaaab35(xyzzyaaac35+1)=floor(xyzzyaaah35/xyzzyaaai35)
if(any(xyzzyaaab35(1:xyzzyaaac35+1)>9).or.any(xyzzyaaab35(1:xyzzyaaac3&
&5+1)<0))then
d2s0='NaN'
return
endif
do
if(xyzzyaaab35(xyzzyaaac35+1)>=5)then
xyzzyaaaf35=1
xyzzyaaab35(xyzzyaaac35+1)=0
else
exit
endif
do xyzzyaaad35=xyzzyaaac35,1,-1
xyzzyaaab35(xyzzyaaad35)=xyzzyaaab35(xyzzyaaad35)+xyzzyaaaf35
xyzzyaaaf35=xyzzyaaab35(xyzzyaaad35)/10
xyzzyaaab35(xyzzyaaad35)=mod(xyzzyaaab35(xyzzyaaad35),10)
enddo
if(xyzzyaaaf35>0)then
do xyzzyaaad35=xyzzyaaac35+1,2,-1
xyzzyaaab35(xyzzyaaad35)=xyzzyaaab35(xyzzyaaad35-1)
enddo
xyzzyaaab35(1)=xyzzyaaaf35
xyzzyaaag35=xyzzyaaag35+1
endif
enddo
d2s0=char(48+xyzzyaaab35(1))//'.'
if(xyzzyaaak35)d2s0=sign_char//trim(d2s0)
do xyzzyaaad35=2,xyzzyaaac35
d2s0=trim(d2s0)//char(48+xyzzyaaab35(xyzzyaaad35))
enddo
if(xyzzyaaag35>0)then
d2s0=trim(d2s0)//exp_char//'+'//trim(i2s(xyzzyaaag35))
elseif(xyzzyaaag35<0)then
d2s0=trim(d2s0)//exp_char//trim(i2s(xyzzyaaag35))
endif
end function d2s0
character(1) function l2s(l)
implicit none
logical,intent(in) :: l
if(l)then
l2s='T'
else
l2s='F'
endif
end function l2s
integer function log2int(log_in)
implicit none
logical,intent(in) :: log_in
if(log_in)then
log2int=1
else
log2int=0
endif
end function log2int
logical function int2log(int_in)
implicit none
integer,intent(in) :: int_in
int2log=(int_in/=0)
end function int2log
subroutine write_list_int(n,list_item,cols,int_length,space_length,wri&
&te_to)
use store,only : o
implicit none
integer,intent(in) :: n,list_item(1:n),int_length,space_length,cols
integer,intent(in),optional :: write_to
integer xyzzyaaaa39,xyzzyaaab39,xyzzyaaac39,xyzzyaaad39
character(2) int_string
character(80) wline
character(1000) space_string,normal_row,last_row
xyzzyaaad39=o
if(present(write_to))xyzzyaaad39=write_to
if(cols<=0.or.cols>=150.or.n<=0.or.int_length<=0.or.int_length>9.or.sp&
&ace_length<=0)then
call wout('N/A',unit=xyzzyaaad39)
return
endif
int_string='i'//char(48+int_length)
space_string='"'
do xyzzyaaaa39=1,space_length
space_string=' '//trim(space_string)
enddo
space_string='"'//trim(space_string)
xyzzyaaab39=n/cols+1
xyzzyaaac39=mod(n,cols)
if(xyzzyaaab39>1)then
normal_row='('
if(cols>1)then
do xyzzyaaaa39=1,cols-1
normal_row=trim(normal_row)//trim(space_string)//','//int_string//','
enddo
endif
normal_row=trim(normal_row)//trim(space_string)//','//int_string//')'
do xyzzyaaaa39=1,xyzzyaaab39-1
write(wline,fmt=trim(normal_row))list_item(1+(xyzzyaaaa39-1)*cols:xyzz&
&yaaaa39*cols)
call wout(wline,unit=xyzzyaaad39,fmt='(a)')
enddo
endif
if(xyzzyaaac39>=1)then
last_row='('
if(xyzzyaaac39>1)then
do xyzzyaaaa39=1,xyzzyaaac39
last_row=trim(last_row)//trim(space_string)//','//int_string//','
enddo
endif
last_row=trim(last_row)//trim(space_string)//','//int_string//')'
write(wline,trim(last_row))list_item(n-xyzzyaaac39+1:n)
call wout(wline,unit=xyzzyaaad39,fmt='(a)')
endif
end subroutine write_list_int
subroutine labelled_list(rdata,label,flag,comment)
use dsp
implicit none
real(dp),intent(in) :: rdata(:)
logical,intent(in),optional :: flag(:)
character(2),intent(in) :: label(:)
character(*),intent(in),optional :: comment
integer xyzzyaaaa40,xyzzyaaab40
integer,parameter :: xyzzyaaac40=4
logical xyzzyaaad40
character(2) char2
character(80) line
xyzzyaaab40=size(rdata,1)
if(xyzzyaaab40<1)return
line=''
char2=label(1)
xyzzyaaad40=.false.
do xyzzyaaaa40=1,xyzzyaaab40
if(xyzzyaaaa40>1)then
char2='  '
if(label(xyzzyaaaa40)/=label(xyzzyaaaa40-1))char2=label(xyzzyaaaa40)
endif
if(present(flag))xyzzyaaad40=flag(xyzzyaaaa40)
if(xyzzyaaad40)then
line=trim(line)//' '//char2//' '//trim(r2s2(rdata(xyzzyaaaa40),'(es15.&
&7)'))//'*'
else
line=trim(line)//' '//char2//' '//trim(r2s2(rdata(xyzzyaaaa40),'(es16.&
&8)'))
endif
if(mod(xyzzyaaaa40,xyzzyaaac40)==0)then
call wout(line,fmt='(a)')
line=''
endif
enddo
if(len_trim(line)>0)call wout(line,fmt='(a)')
if(present(flag).and.present(comment))then
if(any(flag))call wout('   [*] : '//trim(adjustl(comment)))
endif
end subroutine labelled_list
character(r2s_length) function byte2human(nbyte,binary_units)
use dsp
implicit none
real(dp),intent(in) :: nbyte
logical,intent(in),optional :: binary_units
integer xyzzyaaaa41,xyzzyaaab41,xyzzyaaac41
real(dp) xyzzyaaad41,xyzzyaaae41
logical xyzzyaaaf41
character(40) mem_string,suffix
if(nbyte==0.d0)then
byte2human='0.00 B'
return
endif
xyzzyaaaf41=.true.
if(present(binary_units))xyzzyaaaf41=binary_units
xyzzyaaae41=1024.d0
if(.not.xyzzyaaaf41)xyzzyaaae41=1.d3
xyzzyaaaa41=0
xyzzyaaad41=nbyte
do while(xyzzyaaad41>=1.d3)
xyzzyaaad41=xyzzyaaad41/xyzzyaaae41
xyzzyaaaa41=xyzzyaaaa41+3
enddo
xyzzyaaab41=int(xyzzyaaad41)
if(xyzzyaaab41>99)then
mem_string=trim(i2s(xyzzyaaab41))
elseif(xyzzyaaab41>9)then
xyzzyaaab41=int(xyzzyaaad41*1.d1)
mem_string=trim(i2s(xyzzyaaab41/10))//'.'//trim(i2s(mod(xyzzyaaab41,10&
&)))
elseif(xyzzyaaab41>=1)then
xyzzyaaab41=int(xyzzyaaad41*1.d2)
xyzzyaaac41=mod(xyzzyaaab41,100)
mem_string=trim(i2s(xyzzyaaab41/100))//'.'
if(xyzzyaaac41<10)then
mem_string=trim(mem_string)//'0'//trim(i2s(xyzzyaaac41))
else
mem_string=trim(mem_string)//trim(i2s(xyzzyaaac41))
endif
else
xyzzyaaab41=int(xyzzyaaad41*1.d3)
mem_string='0.'//trim(i2s(xyzzyaaab41))
endif
if(xyzzyaaaf41)then
select case(xyzzyaaaa41)
case(0)
suffix='B'
case(3)
suffix='KiB'
case(6)
suffix='MiB'
case(9)
suffix='GiB'
case(12)
suffix='TiB'
case(15)
suffix='PiB'
case(18)
suffix='EiB'
case(21)
suffix='ZiB'
case(24)
suffix='YiB'
case default
mem_string=trim(mem_string)//'E'//trim(i2s(xyzzyaaaa41-24))
suffix='YiB'
end select
else
select case(xyzzyaaaa41)
case(0)
suffix='B'
case(3)
suffix='kB'
case(6)
suffix='MB'
case(9)
suffix='GB'
case(12)
suffix='TB'
case(15)
suffix='PB'
case(18)
suffix='EB'
case(21)
suffix='ZB'
case(24)
suffix='YB'
case default
mem_string=trim(mem_string)//'E'//trim(i2s(xyzzyaaaa41-24))
suffix='YB'
end select
endif
byte2human=trim(mem_string)//' '//trim(suffix)
end function byte2human
character(r2s_length) function time2human(time_in_seconds)
use dsp
implicit none
real(sp),intent(in) :: time_in_seconds
integer(i64) xyzzyaaaa42
integer,parameter :: xyzzyaaab42=4,xyzzyaaac42=2
integer xyzzyaaad42(xyzzyaaab42),xyzzyaaae42,xyzzyaaaf42
integer,parameter :: xyzzyaaag42(xyzzyaaab42-1)=(/60,60,24/)
character,parameter :: xyzzyaaah42(xyzzyaaab42)=(/'s','m','h','d'/)
xyzzyaaaa42=int(time_in_seconds+sign(1.0,time_in_seconds)*0.5,i64)
xyzzyaaad42=0
do xyzzyaaae42=1,xyzzyaaab42-1
xyzzyaaad42(xyzzyaaae42)=int(mod(xyzzyaaaa42,int(xyzzyaaag42(xyzzyaaae&
&42),i64)))
xyzzyaaaa42=xyzzyaaaa42/int(xyzzyaaag42(xyzzyaaae42),i64)
if(xyzzyaaaa42==0_i64)exit
enddo
xyzzyaaaf42=xyzzyaaae42
if(xyzzyaaaa42/=0_i64)then
xyzzyaaaf42=xyzzyaaab42
xyzzyaaad42(xyzzyaaaf42)=int(xyzzyaaaa42)
endif
time2human=trim(i2s(xyzzyaaad42(xyzzyaaaf42)))//xyzzyaaah42(xyzzyaaaf4&
&2)
do xyzzyaaae42=xyzzyaaaf42-1,max(xyzzyaaaf42-xyzzyaaac42+1,1),-1
if(xyzzyaaad42(xyzzyaaae42)==0)cycle
time2human=trim(time2human)//trim(i2s(xyzzyaaad42(xyzzyaaae42)))//xyzz&
&yaaah42(xyzzyaaae42)
enddo
end function time2human
subroutine display_param(value,optimizable,param_name,fmtstring_in,ind&
&ent_in,comment_in,write_to_in)
use dsp
use store,only : o
implicit none
integer,intent(in) :: optimizable
integer,intent(in),optional :: indent_in,write_to_in
real(dp),intent(in) :: value
character(*),intent(in) :: param_name
character(*),intent(in),optional :: fmtstring_in,comment_in
integer xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43
character(30) string2
character(34) string1
character(80) tmpr,fmtstring,string3
xyzzyaaab43=o
if(present(write_to_in))xyzzyaaab43=write_to_in
xyzzyaaaa43=3
if(present(indent_in))xyzzyaaaa43=indent_in
if(xyzzyaaaa43<0.or.xyzzyaaaa43>10)then
call wout('Bug: INDENT should be between 0 and 10.',unit=xyzzyaaab43)
return
endif
if(present(fmtstring_in))then
if(len_trim(adjustl(fmtstring_in))>len(fmtstring))then
call wout('Bug: format string is too long.',unit=xyzzyaaab43)
return
endif
fmtstring=trim(adjustl(fmtstring_in))
else
fmtstring='(es23.15)'
endif
if(present(comment_in))then
if(len_trim(adjustl(comment_in))>len(string3))then
call wout('Comment is too long.',unit=xyzzyaaab43)
return
endif
string3=trim(adjustl(comment_in))
else
string3=''
endif
if(len_trim(adjustl(param_name))+xyzzyaaaa43+15>len(string1))then
call wout('Parameter name is too long.',unit=xyzzyaaab43)
return
endif
string1=trim(adjustl(param_name))
do xyzzyaaac43=1,xyzzyaaaa43-1
string1=' '//trim(string1)
enddo
tmpr=r2s(value,trim(fmtstring))
if(optimizable==1)then
string1(len(string1)-12:len(string1))='(optimizable)'
elseif(optimizable==0)then
string1(len(string1)-6:len(string1))='(fixed)'
elseif(optimizable==-1)then
string1(len(string1)-11:len(string1))='(determined)'
endif
if(value<0.d0)then
string2=' : '//trim(tmpr)
else
string2=' :  '//trim(tmpr)
endif
call wout(trim(string1//string2//string3),unit=xyzzyaaab43)
end subroutine display_param
subroutine print_centred_line(line,width)
implicit none
integer,intent(in) :: width
character(*),intent(in) :: line
integer xyzzyaaaa44,xyzzyaaab44
xyzzyaaaa44=len(line)
xyzzyaaab44=max(0,(width-xyzzyaaaa44)/2)
call wout(repeat(' ',xyzzyaaab44)//trim(line),fmt='(a)')
end subroutine print_centred_line
function get_field(n,line) result(field)
implicit none
integer,intent(in) :: n
character(*),intent(in) :: line
integer xyzzyaaaa45,xyzzyaaab45,xyzzyaaac45
logical xyzzyaaad45
character(len(line)) field,tline
if(n==0)then
field=''
return
endif
xyzzyaaac45=abs(n)
tline=trim(adjustl(line))
xyzzyaaad45=(n<0)
do xyzzyaaaa45=1,xyzzyaaac45-1
xyzzyaaab45=scan(trim(adjustl(tline)),' ',xyzzyaaad45)
if(xyzzyaaab45==0)then
field=''
return
endif
if(xyzzyaaad45)then
tline=trim(adjustl(tline(1:xyzzyaaab45-1)))
else
tline=trim(adjustl(tline(xyzzyaaab45+1:)))
endif
enddo
xyzzyaaab45=scan(trim(adjustl(tline)),' ',xyzzyaaad45)
if(xyzzyaaab45==0)then
field=trim(adjustl(tline))
elseif(xyzzyaaad45)then
field=trim(adjustl(tline(xyzzyaaab45+1:)))
else
field=trim(adjustl(tline(1:xyzzyaaab45-1)))
endif
end function get_field
function pad_int(n,nmax,pad_with) result(output)
implicit none
integer,intent(in) :: n,nmax
character,optional,intent(in) :: pad_with
integer xyzzyaaaa46
character(20) output
output=i2s(n)
xyzzyaaaa46=len_trim(i2s(nmax))-len_trim(output)
if(present(pad_with))then
output=repeat(pad_with,xyzzyaaaa46)//trim(output)
else
output=repeat(' ',xyzzyaaaa46)//trim(output)
endif
end function pad_int
end module format_utils

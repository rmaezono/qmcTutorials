module run_control
use dsp
use format_utils,only : wout,i2s,wordwrap,time2human
use parallel,    only : am_master,my_node
use store,       only : max_cpu_time,max_real_time
implicit none
private
public timer,timer_start,timer_end,tcputime,walltime,exceeds_time_limi&
&t,errstop,errstop_master,errstop2,check_alloc,errwarn,qmc_stop,loop_t&
&ime_estimate,time_report,traceback,errwarn_silent,get_total_time,tota&
&l_time_previous
integer :: xyzzyaaaa1=0
real(sp) xyzzyaaab1(3),xyzzyaaac1,total_time_previous
real(sp),parameter :: xyzzyaaad1=0.001
integer,parameter :: xyzzyaaae1=20
logical,parameter :: xyzzyaaaf1=.false.
logical,parameter :: xyzzyaaag1=.false.
character(30) :: xyzzyaaah1(1:xyzzyaaae1)='(unknown)'
type timed_routine
real(sp) time(3),start(3),rtime(3)
character(30) name
logical error,collapse
integer(i64) hits
type(timed_routine),pointer :: first_sublevel,next,parent
end type timed_routine
type(timed_routine),pointer :: xyzzyaaai1,xyzzyaaaj1
contains
subroutine timer(label,activate,collapse)
use parallel, only : use_timer,am_master
implicit none
logical,intent(in) :: activate
logical,intent(in),optional :: collapse
character(*),intent(in) :: label
real(sp) xyzzyaaaa2(3),xyzzyaaab2(3)
logical err
type(timed_routine),pointer :: xyzzyaaac2
if(xyzzyaaag1.and.am_master)then
if(activate)call wout('+'//repeat(' ',xyzzyaaaa1)//trim(label),fmt='(a&
&)')
endif
if(.not.use_timer)then
if(activate)then
xyzzyaaaa1=min(xyzzyaaaa1+1,xyzzyaaae1)
xyzzyaaah1(xyzzyaaaa1)=trim(label)
else
xyzzyaaaa1=max(xyzzyaaaa1-1,0)
endif
return
endif
if(activate)then
err=.false.
if(xyzzyaaaa1>=xyzzyaaae1)return
xyzzyaaaa1=xyzzyaaaa1+1
xyzzyaaah1(xyzzyaaaa1)=trim(label)
xyzzyaaac2=>xyzzyaaaj1%first_sublevel
if(associated(xyzzyaaac2,xyzzyaaaj1))then
call xyzzyaaak1(xyzzyaaaj1,xyzzyaaac2,label,.true.,err)
else
do
if(trim(xyzzyaaac2%name)==label)exit
if(associated(xyzzyaaac2%next,xyzzyaaaj1))then
call xyzzyaaak1(xyzzyaaaj1,xyzzyaaac2,label,.false.,err)
exit
endif
xyzzyaaac2=>xyzzyaaac2%next
enddo
endif
if(err)return
xyzzyaaaj1=>xyzzyaaac2
xyzzyaaaj1%hits=xyzzyaaaj1%hits+1_i64
if(present(collapse))xyzzyaaaj1%collapse=xyzzyaaaj1%collapse.or.collap&
&se
call xyzzyaaam1(xyzzyaaaa2)
xyzzyaaaj1%start(:)=xyzzyaaaa2(:)
else
call xyzzyaaam1(xyzzyaaaa2)
if(trim(xyzzyaaaj1%name)==label)then
xyzzyaaab2(:)=xyzzyaaaa2(:)-xyzzyaaaj1%start(:)
xyzzyaaaj1%time(:)=xyzzyaaaj1%time(:)+xyzzyaaab2(:)
xyzzyaaaj1=>xyzzyaaaj1%parent
xyzzyaaaa1=xyzzyaaaa1-1
else
xyzzyaaac2=>xyzzyaaaj1
do
xyzzyaaac2=>xyzzyaaac2%parent
if(associated(xyzzyaaac2,xyzzyaaai1))return
if(trim(xyzzyaaac2%name)==label)exit
enddo
do
xyzzyaaab2(:)=xyzzyaaaa2(:)-xyzzyaaaj1%start(:)
xyzzyaaaj1%time(:)=xyzzyaaaj1%time(:)+xyzzyaaab2(:)
if(associated(xyzzyaaaj1,xyzzyaaac2))exit
xyzzyaaaj1%error=.true.
xyzzyaaaj1=>xyzzyaaaj1%parent
xyzzyaaaa1=xyzzyaaaa1-1
enddo
endif
endif
end subroutine timer
subroutine xyzzyaaak1(my_parent,new_level,label,first,err)
implicit none
logical,intent(in) :: first
logical,intent(out) :: err
character(*),intent(in) :: label
integer xyzzyaaaa3
type(timed_routine),pointer :: my_parent,new_level,xyzzyaaab3
err=.false.
allocate(xyzzyaaab3,stat=xyzzyaaaa3)
if(xyzzyaaaa3/=0)then
err=.true.
return
endif
if(first)then
my_parent%first_sublevel=>xyzzyaaab3
else
new_level%next=>xyzzyaaab3
endif
new_level=>xyzzyaaab3
new_level%name=label
new_level%time=0.0
new_level%rtime=0.0
new_level%start=0.0
new_level%hits=0
new_level%error=.false.
new_level%collapse=.false.
new_level%first_sublevel=>new_level
new_level%next=>my_parent
new_level%parent=>my_parent
end subroutine xyzzyaaak1
subroutine timer_start
use parallel, only : am_master
implicit none
if(.not.am_master)return
allocate(xyzzyaaai1)
xyzzyaaaj1=>xyzzyaaai1
xyzzyaaai1%name='LEVEL_ZERO'
xyzzyaaai1%time=0.0
xyzzyaaai1%rtime=0.0
xyzzyaaai1%start=0.0
xyzzyaaai1%hits=1
xyzzyaaai1%error=.false.
xyzzyaaai1%collapse=.false.
xyzzyaaai1%first_sublevel=>xyzzyaaai1
xyzzyaaai1%next=>xyzzyaaai1
xyzzyaaai1%parent=>xyzzyaaai1
call xyzzyaaam1(xyzzyaaab1)
xyzzyaaac1=walltime()
total_time_previous=0.0
end subroutine timer_start
subroutine timer_end
use parallel, only : use_timer
implicit none
integer xyzzyaaaa5
integer(i64) hits
real(sp) xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaaaf5(3),xyz&
&zyaaag5
logical xyzzyaaah5
character(40) label
character(80) line
type(timed_routine),pointer :: xyzzyaaai5,xyzzyaaaj5
line=' '
xyzzyaaab5=tcputime()-xyzzyaaab1(3)
xyzzyaaac5=walltime()-xyzzyaaac1
call wout()
call wout('Total CASINO CPU time  : : : ',dble(xyzzyaaab5),rfmt='(f13.&
&4)')
call wout('Total CASINO real time : : : ',dble(xyzzyaaac5),rfmt='(f13.&
&4)')
call wout()
if(total_time_previous/=0.0)then
call wout('Cumulative CPU time    : : : ',dble(xyzzyaaab5+total_time_p&
&revious),rfmt='(f13.4)')
call wout('(summed over restarts)')
call wout()
endif
if(use_timer)then
xyzzyaaaa5=0
xyzzyaaai5=>xyzzyaaai1
xyzzyaaah5=.false.
xyzzyaaai5%time(3)=xyzzyaaab5
call wout('Timing information:')
call wout('----------------------------------CPU---Sys/CPU--CPU/brn--C&
&PU/tot--Hits--')
mainprint: do
if(associated(xyzzyaaai5%first_sublevel,xyzzyaaai5).or.(xyzzyaaai5%col&
&lapse.and..not.xyzzyaaaf1))then
if(associated(xyzzyaaai5%next,xyzzyaaai5%parent))then
do
if(associated(xyzzyaaai5%parent,xyzzyaaai5))exit mainprint
xyzzyaaai5=>xyzzyaaai5%parent
xyzzyaaaa5=xyzzyaaaa5-1
if(xyzzyaaai5%rtime(3)>0.)then
label=repeat('-',xyzzyaaaa5)//'[rest]'
xyzzyaaag5=0.
xyzzyaaae5=0.
xyzzyaaad5=0.
if(xyzzyaaai5%rtime(3)>0.)xyzzyaaag5=min(1.,max(0.,xyzzyaaai5%rtime(2)&
&/xyzzyaaai5%rtime(3)))
if(xyzzyaaab5>0.)xyzzyaaae5=min(1.,max(0.,xyzzyaaai5%rtime(3)/xyzzyaaa&
&b5))
if(xyzzyaaai5%time(3)>0.)xyzzyaaad5=min(1.,max(0.,xyzzyaaai5%rtime(3)/&
&xyzzyaaai5%time(3)))
call xyzzyaaal1(label,xyzzyaaag5,max(0.0,xyzzyaaai5%rtime(3)),xyzzyaaa&
&d5,xyzzyaaae5,-1_i64,line)
call wout(trim(line),fmt='(a)')
endif
if(associated(xyzzyaaai5%next,xyzzyaaai5%parent))cycle
xyzzyaaai5=>xyzzyaaai5%next
exit
enddo
else
xyzzyaaai5=>xyzzyaaai5%next
endif
else
xyzzyaaaf5(:)=0.
xyzzyaaaj5=>xyzzyaaai5%first_sublevel
do
xyzzyaaaf5=xyzzyaaaf5+xyzzyaaaj5%time
if(associated(xyzzyaaaj5%next,xyzzyaaaj5%parent))exit
xyzzyaaaj5=>xyzzyaaaj5%next
enddo
xyzzyaaai5%rtime=xyzzyaaai5%time-xyzzyaaaf5
xyzzyaaai5=>xyzzyaaai5%first_sublevel
xyzzyaaaa5=xyzzyaaaa5+1
endif
xyzzyaaag5=0.
xyzzyaaae5=0.
xyzzyaaad5=0.
if(xyzzyaaai5%time(3)>0.)xyzzyaaag5=min(1.,max(0.,xyzzyaaai5%time(2)/x&
&yzzyaaai5%time(3)))
if(xyzzyaaab5>0.)xyzzyaaae5=min(1.,max(0.,xyzzyaaai5%time(3)/xyzzyaaab&
&5))
xyzzyaaaj5=>xyzzyaaai5%parent
if(xyzzyaaaj5%time(3)>0.)xyzzyaaad5=min(1.,max(0.,xyzzyaaai5%time(3)/x&
&yzzyaaaj5%time(3)))
if(xyzzyaaai5%collapse.and.xyzzyaaaa5>1)then
label=repeat('-',xyzzyaaaa5-2)//'+'//trim(xyzzyaaai5%name)
else
label=repeat('-',xyzzyaaaa5-1)//trim(xyzzyaaai5%name)
endif
if(xyzzyaaai5%error)then
label=trim(label)//'*'
xyzzyaaah5=.true.
endif
hits=xyzzyaaai5%hits
if(xyzzyaaab5>0.)xyzzyaaae5=min(1.,max(0.,xyzzyaaai5%time(3)/xyzzyaaab&
&5))
call xyzzyaaal1(label,xyzzyaaag5,max(0.0,xyzzyaaai5%time(3)),xyzzyaaad&
&5,xyzzyaaae5,hits,line)
call wout(trim(line),fmt='(a)')
enddo mainprint
call wout(repeat('-',73))
call wout('      CPU: CPU time (seconds)     Sys/CPU: system-to-CPU ti&
&me ratio')
call wout('   CPU/brn: CPU time (% of branch)      CPU/tot: CPU time (&
&% of total)')
if(xyzzyaaah5)call wout('             *: routines without correct time&
&r deactivation.')
call wout(repeat('=',73))
else
call wout()
call wout('Subroutine timers deactivated (use TIMING_INFO input keywor&
&d)')
call wout()
call wout(repeat('=',73))
endif
end subroutine timer_end
subroutine xyzzyaaal1(label,s_ratio,c_time,c_ratio_branch,c_ratio_tot,&
&hits,line)
implicit none
integer(i64),intent(in) :: hits
real(sp),intent(in) :: s_ratio,c_time,c_ratio_branch,c_ratio_tot
character(40),intent(in) :: label
character(80),intent(out) :: line
integer xyzzyaaaa6,xyzzyaaab6
integer(i64) xyzzyaaac6
real(sp) xyzzyaaad6
character(40) hits_string,suffix
if(hits>=0_i64)then
xyzzyaaac6=hits
xyzzyaaaa6=0
do
if(xyzzyaaac6>=1000_i64)then
xyzzyaaac6=xyzzyaaac6/1000_i64
xyzzyaaaa6=xyzzyaaaa6+3
else
exit
endif
enddo
if(xyzzyaaac6<10.and.xyzzyaaaa6>0)then
xyzzyaaad6=real(hits)*10.**(-xyzzyaaaa6)
if(nint(xyzzyaaad6)>9)then
hits_string='10'
else
xyzzyaaab6=int(hits/10_i64**(xyzzyaaaa6-1))
hits_string=trim(i2s(xyzzyaaab6/10))//'.'//trim(i2s(mod(xyzzyaaab6,10)&
&))
endif
else
hits_string=trim(i2s(int(xyzzyaaac6)))
endif
select case(xyzzyaaaa6)
case(0)
suffix=''
case(3)
suffix='k'
case(6)
suffix='M'
case(9)
suffix='G'
case(12)
suffix='T'
case default
suffix='E'//trim(i2s(xyzzyaaaa6))
end select
else
hits_string=''
suffix=''
endif
write(line,"(1x,a,t30,f10.2,t42,f6.2,'%',t51,f6.2,'%',t60,f6.2,'%',t69&
&,a3,a)")trim(label),c_time,100.*s_ratio,100.*c_ratio_branch,100.*c_ra&
&tio_tot,trim(hits_string),trim(suffix)
end subroutine xyzzyaaal1
logical function exceeds_time_limit(first_block)
use parallel
use store, only : max_cpu_time,max_real_time
implicit none
logical,intent(in) :: first_block
integer,save :: xyzzyaaaa7
real(sp) xyzzyaaab7(3),xyzzyaaac7
real(sp),save :: xyzzyaaad7,xyzzyaaae7,xyzzyaaaf7,xyzzyaaag7
logical,save :: xyzzyaaah7=.true.,xyzzyaaai7,xyzzyaaaj7,xyzzyaaak7
exceeds_time_limit=.false.
if(xyzzyaaah7)then
xyzzyaaai7=max_cpu_time>0.
xyzzyaaaj7=max_real_time>0.
xyzzyaaak7=xyzzyaaai7.or.xyzzyaaaj7
xyzzyaaah7=.false.
endif
if(am_master)then
if(xyzzyaaak7)then
if(first_block)then
xyzzyaaaa7=0
else
xyzzyaaaa7=xyzzyaaaa7+1
endif
endif
if(xyzzyaaai7)then
call xyzzyaaam1(xyzzyaaab7)
if(first_block)then
xyzzyaaae7=xyzzyaaab7(3)
xyzzyaaad7=0.
else
xyzzyaaad7=(xyzzyaaab7(3)-xyzzyaaae7)/real(xyzzyaaaa7,sp)
endif
if(xyzzyaaab7(3)-xyzzyaaab1(3)+xyzzyaaad7*1.2>max_cpu_time)exceeds_tim&
&e_limit=.true.
endif
if(xyzzyaaaj7)then
xyzzyaaac7=walltime()
if(first_block)then
xyzzyaaag7=xyzzyaaac7
xyzzyaaaf7=0.
else
xyzzyaaaf7=(xyzzyaaac7-xyzzyaaag7)/real(xyzzyaaaa7,sp)
endif
if(xyzzyaaac7-xyzzyaaac1+xyzzyaaaf7*1.2>max_real_time)exceeds_time_lim&
&it=.true.
endif
endif
if(xyzzyaaak7)then
call mpi_bcast(exceeds_time_limit,1,mpi_logical,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'Broadcasting exceeds_time_limit in exceeds_time_&
&limit')
endif
end function exceeds_time_limit
subroutine xyzzyaaam1(time)
implicit none
real(sp),intent(out) :: time(3)
real(sp) user_and_system_time(2)
interface
real(kind(1.0)) function etime(user_and_system_time)
real(kind(1.0)),intent(inout) :: user_and_system_time(2)
end function etime
end interface
time(3)=etime(user_and_system_time)
time(1)=user_and_system_time(1)
time(2)=user_and_system_time(2)
end subroutine xyzzyaaam1
real(sp) function tcputime()
implicit none
real(sp) user_and_system_time(2)
interface
real(kind(1.0)) function etime(user_and_system_time)
real(kind(1.0)),intent(inout) :: user_and_system_time(2)
end function etime
end interface
tcputime=etime(user_and_system_time)
end function tcputime
real(sp) function walltime()
implicit none
integer(i64) xyzzyaaaa12
integer(i64),save :: xyzzyaaab12
logical,save :: xyzzyaaac12=.true.
call clock_long(xyzzyaaaa12)
if(xyzzyaaac12)then
xyzzyaaab12=xyzzyaaaa12
xyzzyaaac12=.false.
endif
walltime=real(xyzzyaaaa12-xyzzyaaab12,sp)*xyzzyaaad1
end function walltime
subroutine get_total_time(t)
real(sp),intent(out) :: t
t=tcputime()-xyzzyaaab1(3)
t=t+total_time_previous
end subroutine get_total_time
subroutine clock_long(t)
implicit none
integer(i64),intent(out) :: t
integer xyzzyaaaa14,xyzzyaaab14(8)
t=0_i64
call date_and_time(values=xyzzyaaab14)
do xyzzyaaaa14=2000,xyzzyaaab14(1)-1
t=t+days_upto_month(xyzzyaaaa14,12)
enddo
t=t+days_upto_month(xyzzyaaab14(1),xyzzyaaab14(2)-1)
t=t+xyzzyaaab14(3)-1_i64
t=t*24_i64
t=t+xyzzyaaab14(5)
t=t*60_i64
t=t+xyzzyaaab14(6)
t=t*60_i64
t=t+xyzzyaaab14(7)
t=t*1000_i64
t=t+xyzzyaaab14(8)
end subroutine clock_long
integer(i64) function days_upto_month(year,month)
implicit none
integer,intent(in) :: year,month
integer(i64),parameter :: xyzzyaaaa15(12)=(/31,59,90,120,151,181,212,2&
&43,273,304,334,365/)
days_upto_month=0
if(month<=0.or.month>=13)return
days_upto_month=xyzzyaaaa15(month)
if(month>1.and.((mod(year,4)==0.and.mod(year,100)/=0).or.mod(year,400)&
&==0))days_upto_month=days_upto_month+1_i64
end function days_upto_month
subroutine loop_time_estimate(description,iter,num_iters,report_wallti&
&me,no_newline)
implicit none
integer,intent(in),optional :: iter,num_iters
logical,intent(in),optional :: report_walltime,no_newline
character(*),intent(in) :: description
real(sp) xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16,xyzzyaaae16,x&
&yzzyaaaf16
real(sp),save :: xyzzyaaag16=0.,xyzzyaaah16=0.,xyzzyaaai16=0.,xyzzyaaa&
&j16=0.
real(sp),parameter :: xyzzyaaak16=60.
logical xyzzyaaal16,xyzzyaaam16
logical,save :: xyzzyaaan16=.false.
if(.not.am_master)return
xyzzyaaal16=.false.
if(present(report_walltime))xyzzyaaal16=report_walltime
xyzzyaaam16=.false.
if(present(no_newline))xyzzyaaam16=no_newline
if(.not.present(iter))then
if(xyzzyaaal16)then
xyzzyaaac16=walltime()-xyzzyaaai16
call wout(trim(description)//' [total walltime: '//trim(time2human(xyz&
&zyaaac16))//']')
else
xyzzyaaaa16=tcputime()-xyzzyaaag16
call wout(trim(description)//' [total CPU time: '//trim(time2human(xyz&
&zyaaaa16))//']')
endif
if(.not.xyzzyaaam16)call wout()
elseif(iter==1)then
xyzzyaaag16=tcputime()
xyzzyaaai16=walltime()
xyzzyaaan16=.false.
call wout(trim(description))
elseif(iter==2)then
xyzzyaaah16=tcputime()
xyzzyaaaj16=walltime()
elseif(.not.xyzzyaaan16)then
xyzzyaaae16=tcputime()
xyzzyaaaa16=xyzzyaaae16-xyzzyaaah16
xyzzyaaaf16=walltime()
xyzzyaaac16=xyzzyaaaf16-xyzzyaaaj16
if(xyzzyaaaa16>=xyzzyaaak16)then
xyzzyaaab16=real(num_iters-iter+1)*xyzzyaaaa16/real(iter-2)
xyzzyaaad16=real(num_iters-iter+1)*xyzzyaaac16/real(iter-2)
if(xyzzyaaal16)then
call wout(' [walltime: '//trim(time2human(xyzzyaaaa16))//' elapsed, '/&
&/trim(time2human(xyzzyaaab16))//' remaining]')
else
call wout(' [CPU time: '//trim(time2human(xyzzyaaaa16))//' elapsed, '/&
&/trim(time2human(xyzzyaaab16))//' remaining]')
endif
xyzzyaaan16=.true.
if(max_cpu_time>0.)then
if(max_cpu_time+xyzzyaaab1(3)<xyzzyaaae16+xyzzyaaab16)call errstop('LO&
&OP_TIME_ESTIMATE','Estimated CPU time exceeds available time.')
endif
if(max_real_time>0.)then
if(max_real_time+xyzzyaaac1<xyzzyaaaf16+xyzzyaaad16)call errstop('LOOP&
&_TIME_ESTIMATE','Estimated walltime exceeds available time.')
endif
endif
endif
end subroutine loop_time_estimate
subroutine time_report(description,is_start,report_walltime,no_newline&
&)
implicit none
logical,intent(in),optional :: is_start,report_walltime,no_newline
character(*),intent(in) :: description
real(sp) xyzzyaaaa17,xyzzyaaab17
real(sp),save :: xyzzyaaac17=0.,xyzzyaaad17=0.
logical xyzzyaaae17,xyzzyaaaf17,xyzzyaaag17
xyzzyaaae17=.false.
if(present(is_start))xyzzyaaae17=is_start
xyzzyaaaf17=.false.
if(present(report_walltime))xyzzyaaaf17=report_walltime
xyzzyaaag17=.false.
if(present(no_newline))xyzzyaaag17=no_newline
if(xyzzyaaae17)then
xyzzyaaac17=tcputime()
xyzzyaaad17=walltime()
call wout(trim(description))
else
if(xyzzyaaaf17)then
xyzzyaaab17=walltime()-xyzzyaaad17
call wout(trim(description)//' [total walltime: '//trim(time2human(xyz&
&zyaaab17))//']')
else
xyzzyaaaa17=tcputime()-xyzzyaaac17
call wout(trim(description)//' [total CPU time: '//trim(time2human(xyz&
&zyaaaa17))//']')
endif
if(.not.xyzzyaaag17)call wout()
endif
end subroutine time_report
subroutine traceback
use parallel,only : use_timer
implicit none
integer xyzzyaaaa18
type(timed_routine),pointer :: xyzzyaaab18
call wout('CASINO internal traceback:')
if(.not.use_timer)then
if(xyzzyaaaa1>0)then
call wout(' Problem detected at '//trim(xyzzyaaah1(xyzzyaaaa1)))
else
call wout(' Problem detected at an unknown routine')
endif
do xyzzyaaaa18=xyzzyaaaa1-1,1,-1
call wout(' Called from '//trim(xyzzyaaah1(xyzzyaaaa18)))
enddo
else
call wout(' Problem detected at '//trim(xyzzyaaaj1%name))
xyzzyaaab18=>xyzzyaaaj1
do
if(associated(xyzzyaaab18%parent,xyzzyaaai1))exit
xyzzyaaab18=>xyzzyaaab18%parent
call wout(' Called from '//trim(xyzzyaaab18%name))
enddo
endif
call wout(' Called from MAIN')
end subroutine traceback
subroutine errstop(routine,error)
implicit none
character(*),intent(in) :: routine,error
call wout()
call wout('ERROR : '//trim(routine))
call wordwrap(trim(error))
call wout()
call traceback
call wout()
call wout(repeat('-',78))
call wout()
call qmc_abort
end subroutine errstop
subroutine check_alloc(ialloc,routine,symbol,suggestion)
implicit none
integer,intent(in) :: ialloc
character(*),intent(in) :: routine,symbol
character(*),intent(in),optional :: suggestion
if(ialloc==0)return
call wout()
if(len_trim(symbol)==0)then
call wout('ERROR : Allocation problem in '//trim(routine)//' on node #&
&'//trim(i2s(my_node))//'.')
else
call wout('ERROR : Allocation problem in '//trim(routine)//' ('//trim(&
&symbol)//') on node #'//trim(i2s(my_node))//'.')
endif
call wout()
if(present(suggestion))then
call wout('SUGGESTION : '//trim(suggestion))
call wout()
endif
call traceback
call wout()
call wout(repeat('-',78))
call wout()
call qmc_abort
end subroutine check_alloc
subroutine errstop_master(routine,error)
use parallel
implicit none
character(*),intent(in) :: routine,error
if(am_master)then
call wout()
call wout('ERROR : '//trim(routine))
call wordwrap(trim(error))
call wout()
call traceback
call wout()
call wout(repeat('-',78))
if(nnodes>1)then
call wout()
call wout('Waiting for all nodes to exit.')
endif
endif
call qmc_barrier
if(nnodes>1.and.am_master)then
call wout()
call wout('Done.')
endif
call qmc_stop
end subroutine errstop_master
subroutine errstop2(routine,error,int1)
implicit none
integer,intent(in) :: int1
character(*),intent(in) :: routine,error
call wout()
call wout('ERROR : '//routine)
call wout(error//' '//trim(i2s(int1)))
call wout()
call traceback
call wout()
call wout(repeat('-',78))
call wout()
call qmc_abort
end subroutine errstop2
subroutine errwarn(routine,warning)
use parallel
implicit none
character(*),intent(in) :: routine,warning
if(am_master)then
call wordwrap('Warning: ['//trim(routine)//'] '//trim(warning))
call wout()
endif
end subroutine errwarn
subroutine errwarn_silent(routine,warning)
use parallel
implicit none
character(*),intent(in) :: routine,warning
if(am_master)then
call wordwrap('Warning : ['//trim(routine)//'] '//trim(warning))
call wout()
endif
end subroutine errwarn_silent
subroutine qmc_abort
use parallel
implicit none
call mpi_abort(mpi_comm_world,-1,ierror)
stop
end subroutine qmc_abort
subroutine qmc_stop
use parallel
implicit none
call mpi_finalize(ierror)
stop
end subroutine qmc_stop
end module run_control

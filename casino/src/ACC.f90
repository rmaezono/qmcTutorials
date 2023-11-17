module slaarnacc
use dsp
use parallel
use format_utils,only : wout,i2s
use run_control, only : errstop,errstop_master,errwarn,timer,check_all&
&oc
use store,       only : open_unit
implicit none
private
public ranx,ranx_gaussian,initialize_random,get_random_state,put_rando&
&m_state,ranx_max,ranx_gauss_max,ranx_buffer,ranx_gauss_buffer,ranx_in&
&dx,ranx_gauss_indx,random_seed_kw
integer,save :: ranx_indx,ranx_gauss_indx,xyzzyaaaa1,xyzzyaaab1,xyzzya&
&aac1
integer,parameter :: ranx_max=63,ranx_gauss_max=20
real ranx_buffer(ranx_max),ranx_gauss_buffer(ranx_gauss_max)
character(20) random_seed_kw
integer,parameter :: xyzzyaaad1=4,xyzzyaaae1=3,xyzzyaaaf1=2**24,xyzzya&
&aag1=2147483563
integer :: xyzzyaaah1(0:xyzzyaaad1)=(/0,24,73,199,365/),xyzzyaaai1=100&
&0000000,xyzzyaaaj1=24,xyzzyaaak1=10,xyzzyaaal1=0,xyzzyaaam1=0,xyzzyaa&
&an1=0,xyzzyaaao1
integer,save :: next(24),xyzzyaaap1=xyzzyaaae1,xyzzyaaaq1,xyzzyaaar1
real,save :: xyzzyaaas1(24),xyzzyaaat1=0.,xyzzyaaau1,xyzzyaaav1
real,parameter :: xyzzyaaaw1=4096.
logical,save :: xyzzyaaax1=.true.
contains
subroutine initialize_random(ranluxlevel_in,ranprint_in,ranlog_in)
implicit none
integer,intent(in) :: ranluxlevel_in,ranprint_in,ranlog_in
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2
real(dp) xyzzyaaag2
logical,save :: xyzzyaaah2=.true.
integer,parameter :: xyzzyaaai2=16807,xyzzyaaaj2=2147483647,xyzzyaaak2&
&=127773,xyzzyaaal2=2836,mask=123459786
xyzzyaaac1=ranluxlevel_in
xyzzyaaaa1=ranprint_in
xyzzyaaab1=ranlog_in
if(.not.xyzzyaaah2.and.trim(random_seed_kw)/='timer_reset')call errsto&
&p_master('INITIALIZE_RANDOM','This subroutine should only be called o&
&nce per run unless RANDOM_SEED is "timer_reset".')
select case(trim(random_seed_kw))
case('standard')
xyzzyaaao1=314159265
case('timer','timer_reset')
if(am_master)call system_clock(xyzzyaaao1,xyzzyaaad2,xyzzyaaae2)
call mpi_bcast(xyzzyaaao1,1,mpi_integer,0,mpi_comm_world,xyzzyaaaf2)
call checkmpi(ierror,'Broadcasting main_random_seed in INITIALIZE_RAND&
&OM')
case default
read(random_seed_kw,*,iostat=xyzzyaaaf2)xyzzyaaao1
if(xyzzyaaaf2/=0)call errstop_master('SET_INPUT_PARAMETERS','Problem w&
&ith RANDOM_SEED value.')
random_seed_kw='integer'
end select
xyzzyaaac2=xyzzyaaao1
xyzzyaaac2=ieor(xyzzyaaac2,mask)
do xyzzyaaaa2=1,my_node
xyzzyaaab2=xyzzyaaac2/xyzzyaaak2
xyzzyaaac2=xyzzyaaai2*(xyzzyaaac2-xyzzyaaab2*xyzzyaaak2)-xyzzyaaal2*xy&
&zzyaaab2
if(xyzzyaaac2<0)xyzzyaaac2=xyzzyaaac2+xyzzyaaaj2
enddo
xyzzyaaac2=ieor(xyzzyaaac2,mask)
call xyzzyaabd1(xyzzyaaac1,xyzzyaaac2,0,0)
ranx_indx=0
ranx_gauss_indx=0
do xyzzyaaaa2=1,24
xyzzyaaag2=ranx()
enddo
if(nnodes<100)then
call xyzzyaaay1
else
if(xyzzyaaaa1>0)call xyzzyaaay1
endif
xyzzyaaah2=.false.
end subroutine initialize_random
subroutine xyzzyaaay1
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3
real(dp),allocatable :: xyzzyaaaf3(:),xyzzyaaag3(:,:)
xyzzyaaad3=1
if(xyzzyaaaa1>1)xyzzyaaad3=xyzzyaaaa1
allocate(xyzzyaaaf3(xyzzyaaad3),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'RANDOM_CHECK_LOG','initial_random_values'&
&)
if(am_master)then
allocate(xyzzyaaag3(xyzzyaaad3,nnodes),stat=xyzzyaaac3)
else
allocate(xyzzyaaag3(1,1),stat=xyzzyaaac3)
endif
call check_alloc(xyzzyaaac3,'RANDOM_CHECK_LOG','random_values')
do xyzzyaaaa3=1,xyzzyaaad3
xyzzyaaaf3(xyzzyaaaa3)=ranx()
enddo
if(xyzzyaaad3==1)ranx_indx=ranx_indx+1
call mpi_gather(xyzzyaaaf3,xyzzyaaad3,mpi_double_precision,xyzzyaaag3,&
&xyzzyaaad3,mpi_double_precision,0,mpi_comm_world,ierror)
if(am_master)then
do xyzzyaaaa3=1,nnodes-1
do xyzzyaaab3=xyzzyaaaa3+1,nnodes
if(xyzzyaaag3(1,xyzzyaaaa3)==xyzzyaaag3(1,xyzzyaaab3))call errwarn('RA&
&NDOM_CHECK_LOG','First random numbers on nodes '//trim(i2s(xyzzyaaaa3&
&))//' and '//trim(i2s(xyzzyaaab3))//' are equal. Continuing regardles&
&s.')
enddo
enddo
if(xyzzyaaaa1>0)then
open(xyzzyaaab1,file='random.log',status='unknown',iostat=xyzzyaaae3)
if(xyzzyaaae3/=0)call errstop('RANDOM_CHECK_LOG','Problem opening rand&
&om.log')
if(nnodes==1)then
do xyzzyaaab3=1,xyzzyaaaa1
write(xyzzyaaab1,*)xyzzyaaag3(xyzzyaaab3,1)
enddo
else
write(xyzzyaaab1,'(1x,a)')'Node ; #ran'
do xyzzyaaaa3=1,nnodes
write(xyzzyaaab1,*)xyzzyaaaa3,xyzzyaaaa1
write(xyzzyaaab1,'((5(1x,f14.10)))')(xyzzyaaag3(xyzzyaaab3,xyzzyaaaa3)&
&,xyzzyaaab3=1,xyzzyaaaa1)
enddo
close(xyzzyaaab1)
open_unit(xyzzyaaab1)=.false.
endif
endif
endif
deallocate(xyzzyaaaf3,xyzzyaaag3)
end subroutine xyzzyaaay1
subroutine get_random_state(isdext)
implicit none
integer,intent(out) :: isdext(25)
integer xyzzyaaaa4
do xyzzyaaaa4=1,24
isdext(xyzzyaaaa4)=int(xyzzyaaas1(xyzzyaaaa4)*xyzzyaaaw1*xyzzyaaaw1)
enddo
isdext(25)=xyzzyaaaj1+100*xyzzyaaak1+10000*xyzzyaaal1+1000000*xyzzyaaa&
&p1
if(xyzzyaaat1>0.)isdext(25)=-isdext(25)
end subroutine get_random_state
subroutine put_random_state(isdext)
implicit none
integer,intent(in) :: isdext(25)
if(any(isdext/=0))then
call xyzzyaaaz1(isdext)
else
call xyzzyaaba1
endif
end subroutine put_random_state
subroutine xyzzyaaaz1(isdext)
implicit none
integer,intent(in) :: isdext(25)
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6
integer,parameter :: xyzzyaaae6=16807,xyzzyaaaf6=2147483647,xyzzyaaag6&
&=127773,xyzzyaaah6=2836,mask=123459786
if(xyzzyaaax1)then
call wout('Proper results only with initialization from 25 integers ob&
&tained with RLUXOUT.')
xyzzyaaax1=.false.
endif
xyzzyaaau1=1.
do xyzzyaaaa6=1,24
next(xyzzyaaaa6)=xyzzyaaaa6-1
xyzzyaaau1=xyzzyaaau1*0.5
enddo
next(1)=24
xyzzyaaav1=xyzzyaaau1*4096.
do xyzzyaaaa6=1,24
xyzzyaaas1(xyzzyaaaa6)=real(isdext(xyzzyaaaa6))*xyzzyaaau1
enddo
xyzzyaaat1=0.
if(isdext(25)<0)xyzzyaaat1=xyzzyaaau1
xyzzyaaab6=abs(isdext(25))
xyzzyaaaj1=mod(xyzzyaaab6,100)
xyzzyaaab6=xyzzyaaab6/100
xyzzyaaak1=mod(xyzzyaaab6,100)
xyzzyaaab6=xyzzyaaab6/100
xyzzyaaal1=mod(xyzzyaaab6,100)
xyzzyaaab6=xyzzyaaab6/100
xyzzyaaap1=xyzzyaaab6
if(xyzzyaaap1<=xyzzyaaad1)then
xyzzyaaaq1=xyzzyaaah1(xyzzyaaap1)
elseif(xyzzyaaap1>=24)then
xyzzyaaaq1=xyzzyaaap1-24
else
xyzzyaaaq1=xyzzyaaah1(xyzzyaaad1)
xyzzyaaap1=xyzzyaaad1
endif
if(xyzzyaaap1/=xyzzyaaac1)then
if(am_master)then
call wout('RANLUXLEVEL changed in input from '//trim(i2s(xyzzyaaap1))/&
&/' to '//trim(i2s(xyzzyaaac1))//'.')
if(xyzzyaaap1>xyzzyaaac1)then
call wout('Reinitializing the random number generator with lower quali&
&ty stream.')
else
call wout('Reinitializing the random number generator with higher qual&
&ity stream.')
endif
call wout()
endif
xyzzyaaad6=xyzzyaaao1
xyzzyaaad6=ieor(xyzzyaaad6,mask)
do xyzzyaaaa6=1,my_node
xyzzyaaac6=xyzzyaaad6/xyzzyaaag6
xyzzyaaad6=xyzzyaaae6*(xyzzyaaad6-xyzzyaaac6*xyzzyaaag6)-xyzzyaaah6*xy&
&zzyaaac6
if(xyzzyaaad6<0)xyzzyaaad6=xyzzyaaad6+xyzzyaaaf6
enddo
xyzzyaaad6=ieor(xyzzyaaad6,mask)
call xyzzyaabd1(xyzzyaaac1,xyzzyaaad6,0,0)
elseif(am_master)then
call wout('Random number generator reset to state in config.in.')
call wout()
endif
end subroutine xyzzyaaaz1
subroutine xyzzyaaba1
implicit none
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7
integer,parameter :: xyzzyaaad7=16807,xyzzyaaae7=2147483647,xyzzyaaaf7&
&=127773,xyzzyaaag7=2836,mask=123459786
real(dp) xyzzyaaah7
xyzzyaaac7=xyzzyaaao1
xyzzyaaac7=ieor(xyzzyaaac7,mask)
do xyzzyaaaa7=1,my_node
xyzzyaaab7=xyzzyaaac7/xyzzyaaaf7
xyzzyaaac7=xyzzyaaad7*(xyzzyaaac7-xyzzyaaab7*xyzzyaaaf7)-xyzzyaaag7*xy&
&zzyaaab7
if(xyzzyaaac7<0)xyzzyaaac7=xyzzyaaac7+xyzzyaaae7
enddo
xyzzyaaac7=ieor(xyzzyaaac7,mask)
call xyzzyaabd1(xyzzyaaac1,xyzzyaaac7,0,0)
do xyzzyaaaa7=1,24
xyzzyaaah7=ranx()
enddo
end subroutine xyzzyaaba1
real(dp) function ranx()
implicit none
if(ranx_indx==0)call xyzzyaabc1
ranx=dble(ranx_buffer(ranx_indx))
ranx_indx=ranx_indx-1
end function ranx
real(dp) function ranx_gaussian(sqrt_variance)
implicit none
real(dp),intent(in) :: sqrt_variance
if(ranx_gauss_indx==0)call xyzzyaabb1
ranx_gaussian=sqrt_variance*dble(ranx_gauss_buffer(ranx_gauss_indx))
ranx_gauss_indx=ranx_gauss_indx-1
end function ranx_gaussian
subroutine xyzzyaabb1
implicit none
integer xyzzyaaaa10
real xyzzyaaab10,xyzzyaaac10,xyzzyaaad10,xyzzyaaae10
do xyzzyaaaa10=1,ranx_gauss_max,2
do
xyzzyaaab10=real(ranx())
xyzzyaaab10=xyzzyaaab10+xyzzyaaab10-1.
xyzzyaaac10=real(ranx())
xyzzyaaac10=xyzzyaaac10+xyzzyaaac10-1.
xyzzyaaad10=xyzzyaaab10*xyzzyaaab10+xyzzyaaac10*xyzzyaaac10
if(xyzzyaaad10<=1..and.xyzzyaaad10>0.)exit
enddo
xyzzyaaae10=log(xyzzyaaad10)/xyzzyaaad10
xyzzyaaae10=sqrt(-xyzzyaaae10-xyzzyaaae10)
ranx_gauss_buffer(xyzzyaaaa10)=xyzzyaaac10*xyzzyaaae10
if(xyzzyaaaa10<ranx_gauss_max)ranx_gauss_buffer(xyzzyaaaa10+1)=xyzzyaa&
&ab10*xyzzyaaae10
enddo
ranx_gauss_indx=ranx_gauss_max
end subroutine xyzzyaabb1
subroutine xyzzyaabc1
implicit none
integer xyzzyaaaa11,xyzzyaaab11
real xyzzyaaac11
if(xyzzyaaax1)call errstop('RANLUX','RANLUX called without initializat&
&ion.')
call timer('RANDOM_NUMBER',.true.)
do xyzzyaaaa11=1,ranx_max
xyzzyaaac11=xyzzyaaas1(xyzzyaaak1)-xyzzyaaas1(xyzzyaaaj1)-xyzzyaaat1
if(xyzzyaaac11<0.)then
xyzzyaaac11=xyzzyaaac11+1.0
xyzzyaaat1=xyzzyaaau1
else
xyzzyaaat1=0.
endif
xyzzyaaas1(xyzzyaaaj1)=xyzzyaaac11
xyzzyaaaj1=next(xyzzyaaaj1)
xyzzyaaak1=next(xyzzyaaak1)
ranx_buffer(xyzzyaaaa11)=xyzzyaaac11
if(xyzzyaaac11<xyzzyaaav1)then
ranx_buffer(xyzzyaaaa11)=ranx_buffer(xyzzyaaaa11)+xyzzyaaau1*xyzzyaaas&
&1(xyzzyaaak1)
if(ranx_buffer(xyzzyaaaa11)==0.)ranx_buffer(xyzzyaaaa11)=xyzzyaaau1*xy&
&zzyaaau1
endif
xyzzyaaal1=xyzzyaaal1+1
if(xyzzyaaal1==24)then
xyzzyaaal1=0
xyzzyaaam1=xyzzyaaam1+xyzzyaaaq1
do xyzzyaaab11=1,xyzzyaaaq1
xyzzyaaac11=xyzzyaaas1(xyzzyaaak1)-xyzzyaaas1(xyzzyaaaj1)-xyzzyaaat1
if(xyzzyaaac11<0.)then
xyzzyaaac11=xyzzyaaac11+1.0
xyzzyaaat1=xyzzyaaau1
else
xyzzyaaat1=0.
endif
xyzzyaaas1(xyzzyaaaj1)=xyzzyaaac11
xyzzyaaaj1=next(xyzzyaaaj1)
xyzzyaaak1=next(xyzzyaaak1)
enddo
endif
enddo
xyzzyaaam1=xyzzyaaam1+ranx_max
if(xyzzyaaam1>=xyzzyaaai1)then
xyzzyaaan1=xyzzyaaan1+1
xyzzyaaam1=xyzzyaaam1-xyzzyaaai1
endif
ranx_indx=ranx_max
call timer('RANDOM_NUMBER',.false.)
end subroutine xyzzyaabc1
subroutine xyzzyaabd1(lux,ins,k1,k2)
implicit none
integer,intent(in) :: lux,ins,k1,k2
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12(24),xyzzyaaad12,xyzzyaaae1&
&2,xyzzyaaaf12,xyzzyaaag12,xyzzyaaah12,xyzzyaaai12
real xyzzyaaaj12
xyzzyaaap1=lux
xyzzyaaaq1=xyzzyaaah1(xyzzyaaap1)
xyzzyaaal1=0
if(am_master)then
call wout('Initialize random number generator')
call wout('==================================')
call wout('Generator                                 :  RANLUX')
call wout('RANLUX luxury level                       :  '//trim(i2s(xy&
&zzyaaap1)))
xyzzyaaai12=xyzzyaaaq1+24
call wout('p value                                   :  '//trim(i2s(xy&
&zzyaaai12)))
endif
if(ins<0)call errstop('RANLUX_INIT','Illegal initialization : negative&
& input seed.')
if(ins>0)then
xyzzyaaar1=ins
else
xyzzyaaar1=xyzzyaaao1
endif
if(nnodes==1)then
call wout('Value of random seed                      :  '//trim(i2s(xy&
&zzyaaar1)))
elseif(am_master)then
call wout('Value of random seed on node 0            :  '//trim(i2s(xy&
&zzyaaar1)))
endif
if(am_master)then
call wout()
select case(trim(random_seed_kw))
case('standard')
call wout('Initialized from standard internal seed.')
case('timer','timer_reset')
call wout('Initialized with seed derived from timer.')
case('integer')
call wout('Initialized from integer value of random_seed in input.')
case default
call errstop('RANLUX_INIT','Unknown value of random_seed keyword.')
end select
call wout()
endif
xyzzyaaax1=.false.
xyzzyaaau1=1.
do xyzzyaaaa12=1,24
xyzzyaaau1=xyzzyaaau1*0.5
xyzzyaaae12=xyzzyaaar1/53668
xyzzyaaar1=40014*(xyzzyaaar1-xyzzyaaae12*53668)-xyzzyaaae12*12211
if(xyzzyaaar1<0)xyzzyaaar1=xyzzyaaar1+xyzzyaaag1
xyzzyaaac12(xyzzyaaaa12)=mod(xyzzyaaar1,xyzzyaaaf1)
enddo
xyzzyaaav1=xyzzyaaau1*4096.
do xyzzyaaaa12=1,24
xyzzyaaas1(xyzzyaaaa12)=real(xyzzyaaac12(xyzzyaaaa12))*xyzzyaaau1
next(xyzzyaaaa12)=xyzzyaaaa12-1
enddo
next(1)=24
xyzzyaaaj1=24
xyzzyaaak1=10
xyzzyaaat1=0.
if(xyzzyaaas1(24)==0.)xyzzyaaat1=xyzzyaaau1
xyzzyaaam1=k1
xyzzyaaan1=k2
if(k1+k2/=0)then
do xyzzyaaab12=1,k2+1
xyzzyaaaf12=xyzzyaaai1
if(xyzzyaaab12==k2+1)xyzzyaaaf12=k1
do xyzzyaaad12=1,xyzzyaaaf12
xyzzyaaaj12=xyzzyaaas1(xyzzyaaak1)-xyzzyaaas1(xyzzyaaaj1)-xyzzyaaat1
if(xyzzyaaaj12<0.)then
xyzzyaaaj12=xyzzyaaaj12+1.0
xyzzyaaat1=xyzzyaaau1
else
xyzzyaaat1=0.
endif
xyzzyaaas1(xyzzyaaaj1)=xyzzyaaaj12
xyzzyaaaj1=next(xyzzyaaaj1)
xyzzyaaak1=next(xyzzyaaak1)
enddo
enddo
xyzzyaaal1=mod(xyzzyaaam1,xyzzyaaaq1+24)
if(xyzzyaaan1>0)then
xyzzyaaag12=mod(xyzzyaaai1,xyzzyaaaq1+24)
xyzzyaaah12=xyzzyaaan1*xyzzyaaag12+xyzzyaaal1
xyzzyaaal1=mod(xyzzyaaah12,xyzzyaaaq1+24)
endif
if(xyzzyaaal1>23)then
call wout(' Error in restarting with RLUXGO:')
call wout(' The values '//trim(i2s(ins))//', '//trim(i2s(k1))//', '//t&
&rim(i2s(k2))//' cannot occur at luxury level '//trim(i2s(xyzzyaaap1))&
&)
xyzzyaaal1=0
endif
endif
end subroutine xyzzyaabd1
end module slaarnacc

module slaarnacd
use dsp
use run_control, only : check_alloc,errstop,errwarn
use store,       only : popstats
implicit none
private
public init_reblock,finish_reblock,reblock_add,get_reblocked_ave,reblo&
&ck_dump,popstats_dump,reblock_nstep,save_reblock_data,restore_reblock&
&_data
type reblock_step_data
integer :: block_length=1,steps_until_flush=0
real(dp) :: sum_w2_closed=0.d0,sum_w_open=0.d0
real(dp),pointer :: sum_o2_closed(:)=>null(),sum_o_open(:)=>null()
type(reblock_step_data),pointer :: next=>null()
type(reblock_step_data),pointer :: prev=>null()
real(dp),pointer :: stderr(:)=>null(),dstderr(:)=>null()
end type reblock_step_data
integer reblock_nstep,xyzzyaaaa1
real(dp) xyzzyaaab1
real(dp),allocatable :: xyzzyaaac1(:)
real(dp) xyzzyaaad1
real(dp),allocatable :: xyzzyaaae1(:)
real(dp),allocatable :: xyzzyaaaf1(:)
real(dp),allocatable :: xyzzyaaag1(:)
real(dp),allocatable :: xyzzyaaah1(:)
type(reblock_step_data),pointer :: xyzzyaaai1=>null()
type(reblock_step_data),pointer :: xyzzyaaaj1=>null()
real(dp),allocatable :: xyzzyaaak1(:)
integer,allocatable :: xyzzyaaal1(:)
real(dp),allocatable :: xyzzyaaam1(:),xyzzyaaan1(:)
real(dp),allocatable :: xyzzyaaao1(:)
real(dp),allocatable :: xyzzyaaap1(:)
real(dp),allocatable :: xyzzyaaaq1(:)
real(dp),allocatable :: xyzzyaaar1(:)
real(dp),allocatable :: xyzzyaaas1(:)
real(dp),allocatable :: xyzzyaaat1(:)
contains
subroutine xyzzyaaau1(pt,block_length)
type(reblock_step_data),pointer :: pt
integer,intent(in) :: block_length
integer xyzzyaaaa2
allocate(pt,stat=xyzzyaaaa2)
call check_alloc(xyzzyaaaa2,'NEW_REBLOCK','')
pt%block_length=block_length
pt%steps_until_flush=block_length
pt%sum_w2_closed=0.d0
pt%sum_w_open=0.d0
allocate(pt%sum_o2_closed(xyzzyaaaa1),pt%sum_o_open(xyzzyaaaa1),pt%std&
&err(xyzzyaaaa1),pt%dstderr(xyzzyaaaa1),stat=xyzzyaaaa2)
call check_alloc(xyzzyaaaa2,'NEW_REBLOCK','sum_o')
pt%sum_o2_closed(:)=0.d0
pt%sum_o_open(:)=0.d0
pt%stderr=0.d0
pt%dstderr=0.d0
nullify(pt%next)
nullify(pt%prev)
end subroutine
subroutine init_reblock(nobs)
implicit none
integer,intent(in) :: nobs
integer xyzzyaaaa3
xyzzyaaaa1=nobs
call xyzzyaaau1(xyzzyaaai1,1)
xyzzyaaaj1=>xyzzyaaai1
reblock_nstep=0
xyzzyaaab1=0.d0
allocate(xyzzyaaac1(xyzzyaaaa1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'INIT_REBLOCK','reblock_sum_o')
xyzzyaaac1(:)=0.d0
allocate(xyzzyaaak1(xyzzyaaaa1),xyzzyaaal1(xyzzyaaaa1),xyzzyaaam1(xyzz&
&yaaaa1),xyzzyaaan1(xyzzyaaaa1),xyzzyaaao1(xyzzyaaaa1),stat=xyzzyaaaa3&
&)
call check_alloc(xyzzyaaaa3,'INIT_REBLOCK','reblock_best')
if(popstats)then
allocate(xyzzyaaae1(xyzzyaaaa1),xyzzyaaaf1(xyzzyaaaa1),xyzzyaaag1(xyzz&
&yaaaa1),xyzzyaaah1(xyzzyaaaa1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'INIT_REBLOCK','reblock_sum_ow2')
xyzzyaaad1=0.d0
xyzzyaaae1(:)=0.d0
xyzzyaaaf1(:)=0.d0
xyzzyaaag1(:)=0.d0
xyzzyaaah1(:)=0.d0
allocate(xyzzyaaap1(xyzzyaaaa1),xyzzyaaas1(xyzzyaaaa1),xyzzyaaaq1(xyzz&
&yaaaa1),xyzzyaaar1(xyzzyaaaa1),xyzzyaaat1(xyzzyaaaa1),stat=xyzzyaaaa3&
&)
call check_alloc(xyzzyaaaa3,'INIT_REBLOCK','reblock_ave_w2')
endif
end subroutine init_reblock
subroutine save_reblock_data
use slaarnaaf
implicit none
integer xyzzyaaaa4,xyzzyaaab4
type(reblock_step_data),pointer :: xyzzyaaac4
reblock_nbs_config=0
xyzzyaaac4=>xyzzyaaai1
do while(associated(xyzzyaaac4))
reblock_nbs_config=reblock_nbs_config+1
xyzzyaaac4=>xyzzyaaac4%next
enddo
reblock_nstep_config=reblock_nstep
reblock_nobs_config=xyzzyaaaa1
reblock_sum_w_config=xyzzyaaab1
if(allocated(reblock_block_length_config))deallocate(reblock_block_len&
&gth_config,reblock_s_u_f_config,reblock_sum_w2_closed_config,reblock_&
&sum_w_open_config,reblock_sum_o2_closed_config,reblock_sum_o_open_con&
&fig,reblock_sum_o_config)
allocate(reblock_block_length_config(reblock_nbs_config),reblock_s_u_f&
&_config(reblock_nbs_config),reblock_sum_w2_closed_config(reblock_nbs_&
&config),reblock_sum_w_open_config(reblock_nbs_config),reblock_sum_o2_&
&closed_config(xyzzyaaaa1,reblock_nbs_config),reblock_sum_o_open_confi&
&g(xyzzyaaaa1,reblock_nbs_config),reblock_sum_o_config(xyzzyaaaa1),sta&
&t=xyzzyaaab4)
call check_alloc(xyzzyaaab4,'SAVE_REBLOCK_DATA','reblock_block_length_&
&config')
xyzzyaaaa4=1
xyzzyaaac4=>xyzzyaaai1
do while(associated(xyzzyaaac4))
reblock_block_length_config(xyzzyaaaa4)=xyzzyaaac4%block_length
reblock_s_u_f_config(xyzzyaaaa4)=xyzzyaaac4%steps_until_flush
reblock_sum_w2_closed_config(xyzzyaaaa4)=xyzzyaaac4%sum_w2_closed
reblock_sum_w_open_config(xyzzyaaaa4)=xyzzyaaac4%sum_w_open
reblock_sum_o2_closed_config(:,xyzzyaaaa4)=xyzzyaaac4%sum_o2_closed(:)
reblock_sum_o_open_config(:,xyzzyaaaa4)=xyzzyaaac4%sum_o_open(:)
xyzzyaaaa4=xyzzyaaaa4+1
xyzzyaaac4=>xyzzyaaac4%next
enddo
reblock_sum_o_config(:)=xyzzyaaac1(:)
popstats_config=popstats
if(popstats)then
reblock_sum_w2_config=xyzzyaaad1
if(allocated(reblock_sum_ow0_config))deallocate(reblock_sum_ow0_config&
&,reblock_sum_ow2_config,reblock_sum_o2w0_config,reblock_sum_o2w2_conf&
&ig)
allocate(reblock_sum_ow0_config(xyzzyaaaa1),reblock_sum_ow2_config(xyz&
&zyaaaa1),reblock_sum_o2w0_config(xyzzyaaaa1),reblock_sum_o2w2_config(&
&xyzzyaaaa1),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,'SAVE_REBLOCK_DATA','reblock_sum_ow0_confi&
&g')
reblock_sum_ow0_config(:)=xyzzyaaae1(:)
reblock_sum_ow2_config(:)=xyzzyaaaf1(:)
reblock_sum_o2w0_config(:)=xyzzyaaag1(:)
reblock_sum_o2w2_config(:)=xyzzyaaah1(:)
endif
end subroutine save_reblock_data
subroutine restore_reblock_data
use slaarnaaf
implicit none
integer xyzzyaaaa5
type(reblock_step_data),pointer :: xyzzyaaab5
if(associated(xyzzyaaai1))call finish_reblock
call init_reblock(reblock_nobs_config)
do xyzzyaaaa5=1,reblock_nbs_config
if(xyzzyaaaa5>1)then
call xyzzyaaau1(xyzzyaaab5,reblock_block_length_config(xyzzyaaaa5))
xyzzyaaab5%prev=>xyzzyaaaj1
xyzzyaaaj1%next=>xyzzyaaab5
xyzzyaaaj1=>xyzzyaaab5
endif
xyzzyaaaj1%steps_until_flush=reblock_s_u_f_config(xyzzyaaaa5)
xyzzyaaaj1%sum_w2_closed=reblock_sum_w2_closed_config(xyzzyaaaa5)
xyzzyaaaj1%sum_w_open=reblock_sum_w_open_config(xyzzyaaaa5)
xyzzyaaaj1%sum_o2_closed(:)=reblock_sum_o2_closed_config(:,xyzzyaaaa5)
xyzzyaaaj1%sum_o_open(:)=reblock_sum_o_open_config(:,xyzzyaaaa5)
enddo
reblock_nstep=reblock_nstep_config
xyzzyaaab1=reblock_sum_w_config
xyzzyaaac1(:)=reblock_sum_o_config(:)
if(popstats)then
if(.not.popstats_config)call errstop("RESTORE_REBLOCK_DATA",'could not&
& find expected POPSTATS in config file')
xyzzyaaad1=reblock_sum_w2_config
xyzzyaaae1(:)=reblock_sum_ow0_config(:)
xyzzyaaaf1(:)=reblock_sum_ow2_config(:)
xyzzyaaag1(:)=reblock_sum_o2w0_config(:)
xyzzyaaah1(:)=reblock_sum_o2w2_config(:)
endif
end subroutine restore_reblock_data
recursive subroutine xyzzyaaav1(xyzzyaaab6,sum_o,sum_w,xyzzyaaaa6)
integer,intent(in) :: xyzzyaaaa6
real(dp),intent(in) :: sum_o(xyzzyaaaa1),sum_w
type(reblock_step_data),pointer :: xyzzyaaab6
xyzzyaaab6%sum_o_open(:)=xyzzyaaab6%sum_o_open(:)+sum_o(:)
xyzzyaaab6%sum_w_open=xyzzyaaab6%sum_w_open+sum_w
xyzzyaaab6%steps_until_flush=xyzzyaaab6%steps_until_flush-xyzzyaaaa6
if(xyzzyaaab6%steps_until_flush==0)then
if(.not.associated(xyzzyaaab6%next))then
call xyzzyaaau1(xyzzyaaab6%next,xyzzyaaab6%block_length*2)
xyzzyaaab6%next%prev=>xyzzyaaab6
xyzzyaaaj1=>xyzzyaaab6%next
endif
call xyzzyaaav1(xyzzyaaab6%next,xyzzyaaab6%sum_o_open,xyzzyaaab6%sum_w&
&_open,xyzzyaaab6%block_length)
xyzzyaaab6%steps_until_flush=xyzzyaaab6%block_length
xyzzyaaab6%sum_o2_closed(:)=xyzzyaaab6%sum_o2_closed(:)+(xyzzyaaab6%su&
&m_o_open(:)**2)/xyzzyaaab6%sum_w_open
xyzzyaaab6%sum_w2_closed=xyzzyaaab6%sum_w2_closed+xyzzyaaab6%sum_w_ope&
&n**2
xyzzyaaab6%sum_o_open(:)=0.d0
xyzzyaaab6%sum_w_open=0.d0
endif
end subroutine xyzzyaaav1
subroutine reblock_add(olocal,w)
implicit none
real(dp),intent(in) :: olocal(xyzzyaaaa1),w
reblock_nstep=reblock_nstep+1
xyzzyaaab1=xyzzyaaab1+w
xyzzyaaac1(:)=xyzzyaaac1(:)+olocal(:)
if(popstats)then
xyzzyaaad1=xyzzyaaad1+w**2
xyzzyaaae1(:)=xyzzyaaae1(:)+olocal(:)/w
xyzzyaaaf1(:)=xyzzyaaaf1(:)+w*olocal(:)
xyzzyaaag1(:)=xyzzyaaag1(:)+(olocal(:)/w)**2
xyzzyaaah1(:)=xyzzyaaah1(:)+olocal(:)**2
endif
call xyzzyaaav1(xyzzyaaai1,olocal,w,1)
end subroutine reblock_add
subroutine get_reblocked_ave(o_ave,o_stderr,o_dstderr)
implicit none
real(dp),intent(out) :: o_ave(xyzzyaaaa1),o_stderr(xyzzyaaaa1),o_dstde&
&rr(xyzzyaaaa1)
call xyzzyaaaw1
o_ave(:)=xyzzyaaak1(:)
o_stderr(:)=xyzzyaaam1(:)
o_dstderr(:)=xyzzyaaan1(:)
end subroutine get_reblocked_ave
subroutine xyzzyaaaw1
implicit none
integer xyzzyaaaa9
real(dp) xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,n,sum_w_open
real(dp),dimension(xyzzyaaaa1) :: xyzzyaaae9,xyzzyaaaf9,xyzzyaaag9,xyz&
&zyaaah9,s,ds,sum_o_open,xyzzyaaai9
type(reblock_step_data),pointer :: xyzzyaaaj9
xyzzyaaak1(:)=0.d0
xyzzyaaai9(:)=reblock_nstep
xyzzyaaal1(:)=0
xyzzyaaam1(:)=xyzzyaaai1%stderr(:)
xyzzyaaan1(:)=0.d0
xyzzyaaao1(:)=1.d0
if(xyzzyaaab1<=0.d0.or.reblock_nstep==0)return
xyzzyaaak1(:)=xyzzyaaac1(:)/xyzzyaaab1
if(popstats)then
xyzzyaaae9(:)=xyzzyaaai1%sum_o2_closed(:)
xyzzyaaaf9(:)=0.d0
if(xyzzyaaai1%sum_w_open>0.d0)xyzzyaaaf9(:)=xyzzyaaai1%sum_o_open(:)**&
&2/xyzzyaaai1%sum_w_open
xyzzyaaap1(:)=(xyzzyaaae9(:)+xyzzyaaaf9(:))/xyzzyaaab1
xyzzyaaaq1(:)=xyzzyaaae1(:)/real(reblock_nstep,dp)
xyzzyaaas1(:)=xyzzyaaaf1(:)/xyzzyaaad1
xyzzyaaar1(:)=xyzzyaaag1(:)/real(reblock_nstep,dp)
xyzzyaaat1(:)=xyzzyaaah1(:)/xyzzyaaad1
endif
xyzzyaaaa9=0
sum_w_open=0.d0
sum_o_open(:)=0.d0
xyzzyaaaj9=>xyzzyaaai1
do while(associated(xyzzyaaaj9))
xyzzyaaaj9%stderr(:)=0.d0
xyzzyaaaj9%dstderr(:)=0.d0
if(xyzzyaaaj9%block_length*2>reblock_nstep)exit
sum_w_open=sum_w_open+xyzzyaaaj9%sum_w_open
sum_o_open(:)=sum_o_open(:)+xyzzyaaaj9%sum_o_open(:)
xyzzyaaae9(:)=xyzzyaaaj9%sum_o2_closed(:)
xyzzyaaaf9(:)=0.d0
if(sum_w_open>0.d0)xyzzyaaaf9(:)=sum_o_open(:)**2/sum_w_open
xyzzyaaag9(:)=(xyzzyaaae9(:)+xyzzyaaaf9(:))/xyzzyaaab1
xyzzyaaab9=xyzzyaaaj9%sum_w2_closed
xyzzyaaac9=sum_w_open**2
xyzzyaaad9=(xyzzyaaab9+xyzzyaaac9)/xyzzyaaab1**2
if(xyzzyaaad9<=0.or.1.d0-xyzzyaaad9<=0.d0)exit
xyzzyaaah9(:)=abs(xyzzyaaag9(:)-xyzzyaaak1(:)**2)/(1.d0-xyzzyaaad9)
n=real(reblock_nstep,dp)/real(xyzzyaaaj9%block_length,dp)
s(:)=sqrt(max(xyzzyaaah9(:),0.d0)/n)
ds(:)=s(:)/sqrt(n+n-2.d0)
xyzzyaaaj9%stderr(:)=s(:)
xyzzyaaaj9%dstderr(:)=ds(:)
xyzzyaaaa9=xyzzyaaaj9%block_length
xyzzyaaaj9=>xyzzyaaaj9%next
enddo
xyzzyaaaj9=>xyzzyaaai1
do while(associated(xyzzyaaaj9))
if(xyzzyaaaj9%block_length>xyzzyaaaa9)exit
where(xyzzyaaai1%stderr(:)>0.d0)
xyzzyaaai9(:)=(xyzzyaaaj9%stderr(:)/xyzzyaaai1%stderr(:))**2
endwhere
where(xyzzyaaal1(:)==0.and.real(xyzzyaaaj9%block_length,dp)**3 >= 2*re&
&al(reblock_nstep,dp) * xyzzyaaai9(:)**2)
xyzzyaaal1(:)=xyzzyaaaj9%block_length
xyzzyaaam1(:)=xyzzyaaaj9%stderr(:)
xyzzyaaan1(:)=xyzzyaaaj9%dstderr(:)
xyzzyaaao1(:)=xyzzyaaai9(:)
endwhere
xyzzyaaaj9=>xyzzyaaaj9%next
enddo
end subroutine xyzzyaaaw1
subroutine reblock_dump(idx,ident,unit)
use format_utils, only : wout
use store,        only : o
implicit none
integer,intent(in) :: idx
integer,intent(in),optional :: unit
character(*),intent(in) :: ident
integer xyzzyaaaa10
real(dp) xyzzyaaab10,xyzzyaaac10,xyzzyaaad10
character(80) tmpr
type(reblock_step_data),pointer :: xyzzyaaae10
if(idx<=0)return
xyzzyaaaa10=o
if(present(unit))xyzzyaaaa10=unit
call wout('Dumping reblock data for '//ident//':',unit=xyzzyaaaa10)
write(tmpr,'(a,f20.12,a,f20.12)')'     mean:',xyzzyaaak1(idx),' +/- ',&
&xyzzyaaam1(idx)
call wout(tmpr,unit=xyzzyaaaa10)
write(tmpr,'(a,f20.12,a,f20.12)')'   stderr:',xyzzyaaam1(idx),' +/- ',&
&xyzzyaaan1(idx)
call wout(tmpr,unit=xyzzyaaaa10)
xyzzyaaab10=xyzzyaaai1%stderr(idx)
if(xyzzyaaab10>0.d0)then
xyzzyaaac10=xyzzyaaam1(idx)/xyzzyaaab10
xyzzyaaad10=xyzzyaaan1(idx)/xyzzyaaab10
write(tmpr,'(a,f20.12,a,f20.12)')'   errfac:',xyzzyaaac10,' +/- ',xyzz&
&yaaad10
call wout(tmpr,unit=xyzzyaaaa10)
write(tmpr,'(a,f20.12,a,f20.12)')'   N_corr:',xyzzyaaac10**2,' +/- ',x&
&yzzyaaad10*2*xyzzyaaac10
call wout(tmpr,unit=xyzzyaaaa10)
endif
if(reblock_nstep>1)then
call wout('  ------------------------------------------------------',u&
&nit=xyzzyaaaa10)
call wout('   Block len      Std error   Err in error',unit=xyzzyaaaa1&
&0)
xyzzyaaae10=>xyzzyaaai1
do while(associated(xyzzyaaae10))
if(xyzzyaaae10%block_length==xyzzyaaal1(idx))then
write(tmpr,'(3x,i9,1x,es14.6,1x,es14.6,2x,a)')xyzzyaaae10%block_length&
&,xyzzyaaae10%stderr(idx),xyzzyaaae10%dstderr(idx),'*** BEST ***'
call wout(tmpr,unit=xyzzyaaaa10)
elseif(reblock_nstep/xyzzyaaae10%block_length>=2)then
write(tmpr,'(3x,i9,1x,es14.6,1x,es14.6)')xyzzyaaae10%block_length,xyzz&
&yaaae10%stderr(idx),xyzzyaaae10%dstderr(idx)
call wout(tmpr,unit=xyzzyaaaa10)
endif
xyzzyaaae10=>xyzzyaaae10%next
enddo
call wout('  ------------------------------------------------------',u&
&nit=xyzzyaaaa10)
nullify(xyzzyaaae10)
endif
end subroutine reblock_dump
subroutine popstats_dump(dtdmc,targ_wt,unit)
use slaarnaag,    only : pi
use format_utils, only : wout
use slaarnabg,     only : isperiodic,model_system,npcells
use slaarnabj,      only : tagh_energy,tagh_esqr,tagh_nconf
use store,        only : o
implicit none
integer,intent(in),optional :: unit
real(dp),intent(in) :: dtdmc,targ_wt
integer xyzzyaaaa11
real(dp) xyzzyaaab11,xyzzyaaac11,xyzzyaaad11,xyzzyaaae11,xyzzyaaaf11,x&
&yzzyaaag11,xyzzyaaah11,xyzzyaaai11,xyzzyaaaj11,xyzzyaaak11,xyzzyaaal1&
&1,xyzzyaaam11,xyzzyaaan11,xyzzyaaao11,xyzzyaaap11,xyzzyaaaq11,xyzzyaa&
&ar11
character(80) tmpr
type(reblock_step_data),pointer :: xyzzyaaas11
if(.not.popstats)return
xyzzyaaaa11=o
if(present(unit))xyzzyaaaa11=unit
xyzzyaaas11=>xyzzyaaai1
do while(associated(xyzzyaaas11))
if(xyzzyaaas11%block_length==xyzzyaaal1(tagh_energy))exit
xyzzyaaas11=>xyzzyaaas11%next
enddo
if(.not.associated(xyzzyaaas11))return
xyzzyaaab11=(xyzzyaaas11%stderr(tagh_energy)/xyzzyaaai1%stderr(tagh_en&
&ergy))
xyzzyaaac11=(xyzzyaaas11%dstderr(tagh_energy)/xyzzyaaai1%stderr(tagh_e&
&nergy))
xyzzyaaad11=xyzzyaaab11**2
xyzzyaaae11=xyzzyaaac11*2*xyzzyaaab11
xyzzyaaaf11=dtdmc*xyzzyaaad11
xyzzyaaag11=dtdmc*xyzzyaaae11
xyzzyaaai11=real(reblock_nstep,dp)
xyzzyaaaj11=xyzzyaaab1**2/xyzzyaaad1
xyzzyaaak11=xyzzyaaaq1(tagh_nconf)
xyzzyaaal11=sqrt(abs((xyzzyaaar1(tagh_nconf)-xyzzyaaak11**2))/xyzzyaaa&
&i11)
xyzzyaaah11=xyzzyaaak1(tagh_esqr)-xyzzyaaak1(tagh_energy)**2
if(xyzzyaaah11<=0.d0)return
xyzzyaaam11=xyzzyaaah11/(xyzzyaaah11-(xyzzyaaas1(tagh_esqr)-xyzzyaaat1&
&(tagh_energy)))
xyzzyaaan11=xyzzyaaam11*xyzzyaaaj11/xyzzyaaad11
if(isperiodic.and..not.model_system)xyzzyaaah11=xyzzyaaah11*real(npcel&
&ls,dp)**2
xyzzyaaao11=xyzzyaaah11*sqrt(2/abs(xyzzyaaan11-1))
xyzzyaaap11=exp(sqrt(xyzzyaaah11)*xyzzyaaaf11/sqrt(2*pi))
xyzzyaaaq11=(xyzzyaaao11*xyzzyaaaf11/sqrt(8*pi*xyzzyaaah11)+xyzzyaaag1&
&1*sqrt(xyzzyaaah11/(2*pi)))*xyzzyaaap11
xyzzyaaar11=(xyzzyaaai11*xyzzyaaaq1(tagh_nconf))/(xyzzyaaaj11*xyzzyaaa&
&m11)
call wout('Analysis of statistical efficiency -- see PRB 81, 035119 (2&
&010).',unit=xyzzyaaaa11)
call wout('-----------------------------------------------------------&
&-------',unit=xyzzyaaaa11)
write(tmpr,'(f21.12," +/- ",f21.12)')xyzzyaaad11,xyzzyaaae11
call wout('Int corr length (steps)      = '//trim(tmpr),unit=xyzzyaaaa&
&11)
call wout('DMC time step (au)           = ',dtdmc,rfmt='(f21.12)',unit&
&=xyzzyaaaa11)
write(tmpr,'(f21.12," +/- ",f21.12)')xyzzyaaaf11,xyzzyaaag11
call wout('Int correlation time (au)    = '//trim(tmpr),unit=xyzzyaaaa&
&11)
write(tmpr,'(f21.12," +/- ",f21.12)')xyzzyaaah11,xyzzyaaao11
call wout('Var of loc en (au / simcell) = '//trim(tmpr),unit=xyzzyaaaa&
&11)
call wout('Std dev of local energy      = ',sqrt(xyzzyaaah11),rfmt='(f&
&21.12)',unit=xyzzyaaaa11)
call wout('Number of steps of accum data= ',xyzzyaaai11,rfmt='(f21.12)&
&',unit=xyzzyaaaa11)
call wout('Effective number of steps    = ',xyzzyaaaj11,rfmt='(f21.12)&
&',unit=xyzzyaaaa11)
call wout('Target weight                = ',targ_wt,rfmt='(f21.12)',un&
&it=xyzzyaaaa11)
write(tmpr,'(f21.12," +/- ",f21.12)')xyzzyaaak11,xyzzyaaal11
call wout('Average population           = '//trim(tmpr),unit=xyzzyaaaa&
&11)
call wout('Effective population         = ',xyzzyaaam11,rfmt='(f21.12)&
&',unit=xyzzyaaaa11)
write(tmpr,'(f21.12," +/- ",f21.12)')xyzzyaaap11,xyzzyaaaq11
call wout('Stat inefficiency (est)      = '//trim(tmpr),unit=xyzzyaaaa&
&11)
call wout('Stat inefficiency (measured) = ',xyzzyaaar11,rfmt='(f21.12)&
&',unit=xyzzyaaaa11)
call wout(unit=xyzzyaaaa11)
if(xyzzyaaap11>=2.d0+xyzzyaaaq11.or.xyzzyaaar11>=2.d0+xyzzyaaaq11)call&
& errwarn('POPSTATS','Significant inefficiency due to population corre&
&lation. Be sure to understand the implications.')
end subroutine popstats_dump
subroutine finish_reblock
implicit none
type(reblock_step_data),pointer :: xyzzyaaaa12,xyzzyaaab12
xyzzyaaaa12=>xyzzyaaai1
do while(associated(xyzzyaaaa12))
xyzzyaaab12=>xyzzyaaaa12
xyzzyaaaa12=>xyzzyaaaa12%next
deallocate(xyzzyaaab12%sum_o2_closed)
deallocate(xyzzyaaab12%sum_o_open)
deallocate(xyzzyaaab12%stderr)
deallocate(xyzzyaaab12%dstderr)
deallocate(xyzzyaaab12)
enddo
nullify(xyzzyaaai1,xyzzyaaaj1)
deallocate(xyzzyaaac1)
if(popstats)then
deallocate(xyzzyaaae1,xyzzyaaaf1,xyzzyaaag1,xyzzyaaah1)
deallocate(xyzzyaaap1,xyzzyaaaq1,xyzzyaaas1,xyzzyaaar1,xyzzyaaat1)
endif
xyzzyaaaa1=0
deallocate(xyzzyaaak1,xyzzyaaal1,xyzzyaaam1,xyzzyaaan1,xyzzyaaao1)
end subroutine finish_reblock
end module slaarnacd

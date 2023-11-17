module slaarnacp
use dsp
implicit none
private
public vmc_main,equilibration,vmc_ionjump,corper_default_vmc,corper_de&
&fault_opt
real(dp) vmc_ionjump
integer :: corper_default_opt=15
integer,parameter :: corper_default_vmc=3
logical,parameter :: xyzzyaaaa1=.false.
logical,parameter :: xyzzyaaab1=.true.
logical,parameter :: xyzzyaaac1=.false.
logical,parameter :: xyzzyaaad1=.false.
logical,parameter :: xyzzyaaae1=.false.
logical,parameter :: xyzzyaaaf1=.false.
logical,parameter :: xyzzyaaag1=.false.
character(80) tmpr,tmpr2,tmpr3,tmpr4,tmpr5
contains
subroutine equilibration(is0,is1,rele,sele,dt_array,sqrt_dt_array,dt_s&
&hift_array,opt_dtvmc_array,nequil,corper)
use slaarnaaa
use parallel
use store
use slaarnacs
use slaarnaam,  only : eval_local_energy
use slaarnaaq,        only : particle_is_fixed,pair_corr,pair_corr_sph&
&,pcf_rfix,pcfs_rfix,fixed_particle
use slaarnaan, only : vector_difference
use format_utils,  only : wout,i2s
use slaarnabg,      only : dimensionality
use slaarnabt,     only : parabolic_min,correlation_time,reaverage,dco&
&py,correlation_time_alt
use slaarnaca,         only : have_ppots
use slaarnacc,only : ranx
use run_control,   only : errstop,timer,loop_time_estimate,check_alloc&
&,tcputime
implicit none
integer,intent(in) :: is0,is1,nequil,opt_dtvmc_array(no_difftypes)
integer,intent(inout) :: sele(netot),corper
real(dp),intent(inout) :: dt_array(no_difftypes),sqrt_dt_array(no_diff&
&types),dt_shift_array(no_difftypes),rele(3,netot)
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaal2&
&(no_difftypes),xyzzyaaam2,xyzzyaaan2,xyzzyaaao2,xyzzyaaap2,xyzzyaaaq2
integer,parameter :: xyzzyaaar2=200
integer,allocatable :: xyzzyaaas2(:)
real(sp) xyzzyaaat2
real(dp) xyzzyaaau2(3),xyzzyaaav2,prob,xyzzyaaaw2(no_difftypes),xyzzya&
&aax2(no_difftypes),xyzzyaaay2(no_difftypes),xyzzyaaaz2(no_difftypes),&
&xyzzyaaba2(no_difftypes),xyzzyaabb2(no_difftypes),xyzzyaabc2(no_difft&
&ypes),xyzzyaabd2(no_difftypes),xyzzyaabe2(no_difftypes),xyzzyaabf2,xy&
&zzyaabg2(no_difftypes),xyzzyaabh2,xyzzyaabi2,xyzzyaabj2,xyzzyaabk2(no&
&_difftypes),xyzzyaabl2,xyzzyaabm2,xyzzyaabn2,xyzzyaabo2,xyzzyaabp2,xy&
&zzyaabq2,xyzzyaabr2,xyzzyaabs2(3),xyzzyaabt2,xyzzyaabu2,xyzzyaabv2,xy&
&zzyaabw2,xyzzyaabx2(no_difftypes),xyzzyaaby2(no_difftypes),xyzzyaabz2&
&(no_difftypes),xyzzyaaca2(no_difftypes),xyzzyaacb2(no_difftypes),xyzz&
&yaacc2(no_difftypes),xyzzyaacd2
real(dp),parameter :: xyzzyaace2=0.1d0
real(dp),allocatable :: xyzzyaacf2(:,:),xyzzyaacg2(:),xyzzyaach2(:),xy&
&zzyaaci2(:),xyzzyaacj2(:)
logical xyzzyaack2,isnan,isinf,xyzzyaacl2,xyzzyaacm2,xyzzyaacn2,xyzzya&
&aco2
if(nequil==0)return
call timer('EQUILIBRATION',.true.,collapse=.true.)
if(particle_is_fixed)then
if(pair_corr)then
rele(1:3,fixed_particle)=pcf_rfix(1:3)
elseif(pair_corr_sph)then
rele(1:3,fixed_particle)=pcfs_rfix(1:3)
else
call errstop('EQUILIBRATION','Particle fixed in inappropriate circumst&
&ances.')
endif
endif
call define_config(is0,rele,sele)
if(vmc_cfg_by_cfg)then
allocate(xyzzyaacf2(3,netot),xyzzyaaas2(netot),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'EQUILIBRATION','cfg-by-cfg arrays')
endif
if(xyzzyaaab1)then
xyzzyaaap2=nsampling_levels
else
xyzzyaaap2=1
endif
xyzzyaaak2=no_difftypes
if(vmc_cfg_by_cfg)xyzzyaaak2=1
xyzzyaacd2=1.d0
if(vmc_cfg_by_cfg)xyzzyaacd2=dble(netot)
xyzzyaacm2=any(opt_dtvmc_array(1:xyzzyaaak2)/=0)
xyzzyaacl2=any(opt_dtvmc_array(1:xyzzyaaak2)==2)
xyzzyaaan2=10
if(xyzzyaacl2)xyzzyaaan2=20
xyzzyaaag2=max(1,nequil/(xyzzyaaan2+1))
xyzzyaaaw2=0.d0
xyzzyaaaz2=0.d0
xyzzyaaax2=0.d0
xyzzyaaba2=0.d0
xyzzyaabc2=0.d0
xyzzyaaay2=0.d0
xyzzyaabb2=0.d0
xyzzyaabx2=dt_array
xyzzyaabz2=0.d0
xyzzyaaby2=.5d0
xyzzyaaca2=0.d0
xyzzyaabd2=0.d0
xyzzyaacc2=0.d0
xyzzyaaam2=0
xyzzyaaal2=0
xyzzyaacn2=.false.
xyzzyaaad2=1
if(noncoll_spin)xyzzyaaad2=2
if(use_altsamp)call complextosimple(is0,is1)
xyzzyaaaa2=0
xyzzyaaao2=0
if(xyzzyaaan2==0.and.corper==0)xyzzyaaat2=tcputime()
do
xyzzyaaaa2=xyzzyaaaa2+1
if(xyzzyaaaa2<=nequil)then
if(nequil==1)then
call loop_time_estimate('Running VMC equilibration ('//trim(i2s(nequil&
&))//' move).',xyzzyaaaa2,nequil)
else
call loop_time_estimate('Running VMC equilibration ('//trim(i2s(nequil&
&))//' moves).',xyzzyaaaa2,nequil)
endif
if(xyzzyaacm2.and.xyzzyaaaa2>1.and.mod(xyzzyaaaa2,xyzzyaaag2)==1.and.x&
&yzzyaaam2<xyzzyaaan2)then
xyzzyaaam2=xyzzyaaam2+1
call mpi_reduce(xyzzyaaca2,xyzzyaacb2,no_difftypes,mpi_double_precisio&
&n,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing acc_count in equilibration.')
call mpi_reduce(xyzzyaabd2,xyzzyaabe2,no_difftypes,mpi_double_precisio&
&n,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing sum_dr2 in equilibration.')
if(am_master)then
xyzzyaacc2(1:xyzzyaaak2)=xyzzyaacc2(1:xyzzyaaak2)*dble(nnodes)
where(xyzzyaacc2(1:xyzzyaaak2)/=0.d0)
xyzzyaabk2(1:xyzzyaaak2)=xyzzyaacb2(1:xyzzyaaak2)/xyzzyaacc2(1:xyzzyaa&
&ak2)
xyzzyaabg2(1:xyzzyaaak2)=xyzzyaabe2(1:xyzzyaaak2)/(xyzzyaacd2*xyzzyaac&
&c2(1:xyzzyaaak2))
elsewhere
xyzzyaabk2(1:xyzzyaaak2)=0.d0
xyzzyaabg2(1:xyzzyaaak2)=0.d0
endwhere
endif
xyzzyaaca2=0.d0
xyzzyaabd2=0.d0
xyzzyaacc2=0.d0
if(am_master)then
if(xyzzyaaam2==1)call wout(' Performing time-step optimization.')
call xyzzyaacp2
if(xyzzyaaam2==xyzzyaaan2)then
if(vmc_cfg_by_cfg)then
dt_array(1)=xyzzyaabx2(1)
do xyzzyaaaf2=2,no_difftypes
dt_array(xyzzyaaaf2)=dt_array(1)*difftype_mass(1)/difftype_mass(xyzzya&
&aaf2)
enddo
else
where(opt_dtvmc_array(1:xyzzyaaak2)>0)dt_array(1:xyzzyaaak2)=xyzzyaabx&
&2(1:xyzzyaaak2)
endif
if(no_difftypes==1.or.vmc_cfg_by_cfg)then
call wout(' Optimized DTVMC:',dt_array(1),rfmt='(es12.4)')
else
do xyzzyaaah2=1,no_difftypes
if(opt_dtvmc_array(xyzzyaaah2)==1)then
call wout(' DTVMC #'//trim(i2s(xyzzyaaah2))//': ',dt_array(xyzzyaaah2)&
&,' (optimized)',rfmt='(es12.4)')
else
call wout(' DTVMC #'//trim(i2s(xyzzyaaah2))//': ',dt_array(xyzzyaaah2)&
&,' (fixed)',rfmt='(es12.4)')
endif
enddo
endif
endif
endif
call mpi_bcast(dt_array,no_difftypes,mpi_double_precision,0,mpi_comm_w&
&orld,ierror)
call checkmpi(ierror,'broadcasting dt_array in equilibration.')
sqrt_dt_array(1:no_difftypes)=sqrt(dt_array(1:no_difftypes))
if(xyzzyaaam2==xyzzyaaan2.and.corper==0)xyzzyaaat2=tcputime()
elseif(.not.xyzzyaacm2.and.xyzzyaaaa2==1.and.corper==0)then
xyzzyaaat2=tcputime()
endif
xyzzyaacn2=.false.
elseif(xyzzyaaaa2==nequil+1)then
if(nequil>0)call loop_time_estimate('Done.')
if(corper==0)then
xyzzyaabl2=dble(tcputime()-xyzzyaaat2)
call mpi_reduce(xyzzyaaca2,xyzzyaacb2,no_difftypes,mpi_double_precisio&
&n,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing acc_count in equilibration.')
call mpi_reduce(xyzzyaabl2,xyzzyaabp2,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing time_move in equilibration.')
if(am_master)then
xyzzyaacc2(1:xyzzyaaak2)=xyzzyaacc2(1:xyzzyaaak2)*dble(nnodes)
xyzzyaabk2(1:xyzzyaaak2)=xyzzyaacb2(1:xyzzyaaak2)/xyzzyaacc2(1:xyzzyaa&
&ak2)
if(xyzzyaacm2)then
xyzzyaabl2=xyzzyaabp2/(dble(nequil-xyzzyaaag2*xyzzyaaan2)*dble(nnodes)&
&)
else
xyzzyaabl2=xyzzyaabp2/(dble(nequil)*dble(nnodes))
endif
if(.not.vmc_cfg_by_cfg)then
xyzzyaabt2=1.d0
do xyzzyaaah2=1,no_difftypes
xyzzyaabt2=xyzzyaabt2*(1.d0-xyzzyaabk2(xyzzyaaah2))**sum(nele(:),mask=&
&which_difftype==xyzzyaaah2)
enddo
else
xyzzyaabt2=1.d0-xyzzyaabk2(1)
endif
xyzzyaabt2=max(1.d0-xyzzyaabt2,0.05d0)
xyzzyaaao2=nint(max(xyzzyaace2*dble(nequil),dble(xyzzyaaar2)/xyzzyaabt&
&2))
endif
call mpi_bcast(xyzzyaaao2,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting nequil_add in equilibration.')
allocate(xyzzyaacg2(xyzzyaaao2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'EQUILIBRATION','e_array')
xyzzyaaca2=0.d0
xyzzyaacc2=0.d0
xyzzyaacg2=0.d0
if(use_altsamp)then
allocate(xyzzyaach2(xyzzyaaao2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'EQUILIBRATION','w_array')
xyzzyaach2=0.d0
endif
xyzzyaaat2=tcputime()
else
exit
endif
endif
if(xyzzyaaaa2>nequil.and.xyzzyaaaa2<=nequil+xyzzyaaao2)then
call loop_time_estimate('Finding optimal inner loop length ('//trim(i2&
&s(xyzzyaaao2))//' additional moves).',xyzzyaaaa2-nequil,xyzzyaaao2)
xyzzyaacn2=.true.
xyzzyaacl2=.false.
xyzzyaacm2=.false.
elseif(xyzzyaaaa2==nequil+xyzzyaaao2+1)then
xyzzyaabm2=dble(tcputime()-xyzzyaaat2)
if(use_altsamp)then
call correlation_time_alt(xyzzyaacg2(:),xyzzyaach2(:),xyzzyaabn2,xyzzy&
&aabo2)
else
call correlation_time(xyzzyaacg2(:),xyzzyaabn2,xyzzyaabo2)
endif
deallocate(xyzzyaacg2)
if(use_altsamp)deallocate(xyzzyaach2)
call mpi_reduce(xyzzyaabm2,xyzzyaabq2,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing time_energy in equilibration.')
call mpi_reduce(xyzzyaaca2,xyzzyaacb2,no_difftypes,mpi_double_precisio&
&n,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing acc_count in equilibration.')
if(am_master)then
allocate(xyzzyaaci2(nnodes),xyzzyaacj2(nnodes),stat=xyzzyaaac2)
else
allocate(xyzzyaaci2(1),xyzzyaacj2(1),stat=xyzzyaaac2)
endif
call check_alloc(xyzzyaaac2,'EQUILIBRATION','corr_tau_array/corr_tau_e&
&rr_array')
call mpi_gather(xyzzyaabn2,1,mpi_double_precision,xyzzyaaci2,1,mpi_dou&
&ble_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering corr_tau in equilibration.')
call mpi_gather(xyzzyaabo2,1,mpi_double_precision,xyzzyaacj2,1,mpi_dou&
&ble_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering corr_tau_err in equilibration.')
if(am_master)then
xyzzyaabm2=abs(xyzzyaabq2/(dble(xyzzyaaao2)*dble(nnodes))-xyzzyaabl2)
xyzzyaacc2(1:xyzzyaaak2)=xyzzyaacc2(1:xyzzyaaak2)*dble(nnodes)
xyzzyaabk2(1:xyzzyaaak2)=xyzzyaacb2(1:xyzzyaaak2)/xyzzyaacc2(1:xyzzyaa&
&ak2)
if(any(xyzzyaacj2<0.d0).or.xyzzyaabl2==0.d0.or.xyzzyaabm2==0.d0)then
xyzzyaabn2=-1.d0
xyzzyaabo2=-1.d0
corper=1
else
call reaverage(nnodes,xyzzyaaao2,xyzzyaaci2,xyzzyaacj2,1.d0,xyzzyaabn2&
&,xyzzyaabo2)
call xyzzyaacq2(xyzzyaabn2,xyzzyaabm2,xyzzyaabl2,xyzzyaabk2,vmc_cfg_by&
&_cfg,corper)
endif
if(xyzzyaaag1)then
call wout(' Correlation time : ',(/xyzzyaabn2,xyzzyaabo2/),rfmt='(es12&
&.4)',rsep=' +- ')
call wout(' Acceptance ratios:',xyzzyaabk2(1:xyzzyaaak2),rfmt='(es12.4&
&)')
call wout(' Time per move    : ',xyzzyaabl2,rfmt='(es12.4)')
call wout(' Time per energy  : ',xyzzyaabm2,rfmt='(es12.4)')
endif
call wout(' Optimized vmc_decorr_period: '//trim(i2s(corper)))
endif
deallocate(xyzzyaaci2,xyzzyaacj2)
call mpi_bcast(corper,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting corper in equilibration.')
call loop_time_estimate('Done.')
exit
endif
if(.not.vmc_cfg_by_cfg)then
do xyzzyaaab2=1,netot
if((particle_is_fixed.and.xyzzyaaab2==fixed_particle).or.xyzzyaaab2==o&
&n_top_jj)cycle
if(xyzzyaaad1)then
call wout('PARTICLE #'//trim(i2s(xyzzyaaab2))//':')
call wout(' R_OLD: ',rele(1:dimensionality,xyzzyaaab2))
if(xyzzyaaad2>1)call wout(' S_OLD: '//trim(i2s(sele(xyzzyaaab2))))
endif
xyzzyaaah2=which_difftype(which_spin(xyzzyaaab2))
do xyzzyaaae2=1,xyzzyaaad2
call xyzzyaaak1(xyzzyaaab2,is0,is1,dt_array(xyzzyaaah2),sqrt_dt_array(&
&xyzzyaaah2),dt_shift_array(xyzzyaaah2),xyzzyaaae2,xyzzyaaau2,xyzzyaaa&
&j2,xyzzyaabr2)
call vector_difference(rele(1,xyzzyaaab2),xyzzyaaau2,xyzzyaabs2)
xyzzyaabf2=sum(xyzzyaabs2**2)
if(xyzzyaaad1)then
if(xyzzyaaae2==1)then
call wout(' R_NEW: ',xyzzyaaau2(1:dimensionality))
else
call wout(' S_NEW: '//trim(i2s(xyzzyaaaj2)))
endif
endif
if(altsamp==1.and.have_ppots)then
if(grid_angles_valid(is0))grid_angles_scr(:,:,:,is1)=grid_angles_scr(:&
&,:,:,is0)
grid_angles_valid(is1)=grid_angles_valid(is0)
endif
do xyzzyaaai2=1,xyzzyaaap2
if(xyzzyaaab1)then
xyzzyaaaq2=xyzzyaaai2
else
xyzzyaaaq2=0
endif
xyzzyaaco2=any(rele(:,xyzzyaaab2)/=xyzzyaaau2).or.sele(xyzzyaaab2)/=xy&
&zzyaaaj2
if(xyzzyaaco2)then
if(.not.use_altsamp)then
call wfn_ratio(is0,is1,xyzzyaaaq2,relprob=xyzzyaabh2,isnan=isnan,isinf&
&=isinf)
else
call wfn_altratio(is0,is1,xyzzyaaaq2,relprob=xyzzyaabh2,isnan=isnan,is&
&inf=isinf)
endif
if(xyzzyaaai2==1)xyzzyaabh2=xyzzyaabh2*xyzzyaabr2
prob=min(1.d0,xyzzyaabh2)
if(xyzzyaaad1)then
if(isnan)then
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaaai2))//&
&': zero [was Not a Number]')
elseif(isinf)then
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaaai2))//&
&': zero [it diverged]')
else
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaaai2))//&
&': ',prob)
endif
endif
if(isnan.or.isinf)prob=0.d0
else
prob=1.d0
if(xyzzyaaad1)call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(&
&xyzzyaaai2))//': unity [no-change move]')
endif
xyzzyaaav2=0.d0
if(prob<1.d0)xyzzyaaav2=ranx()
xyzzyaack2=xyzzyaaav2<prob
if(.not.xyzzyaack2)exit
enddo
if(xyzzyaaae2==1)xyzzyaacc2(xyzzyaaah2)=xyzzyaacc2(xyzzyaaah2)+1.d0
if(xyzzyaack2)then
if(xyzzyaaco2)then
call accept_move(is0,is1)
if(use_altsamp)call alt_accept(is0,is1)
rele(:,xyzzyaaab2)=xyzzyaaau2(:)
sele(xyzzyaaab2)=xyzzyaaaj2
if(xyzzyaaab2==on_top_ii)then
rele(:,on_top_jj)=xyzzyaaau2(:)
sele(on_top_jj)=xyzzyaaaj2
endif
endif
if(xyzzyaaae2==1)then
xyzzyaaca2(xyzzyaaah2)=xyzzyaaca2(xyzzyaaah2)+1.d0
xyzzyaabd2(xyzzyaaah2)=xyzzyaabd2(xyzzyaaah2)+xyzzyaabf2
if(xyzzyaaab2==on_top_ii)xyzzyaabd2(which_difftype(which_spin(on_top_j&
&j)))=xyzzyaabd2(which_difftype(which_spin(on_top_jj)))+xyzzyaabf2
endif
if(xyzzyaaad1)call wout(' MOVE ACCEPTED')
else
if(use_altsamp)call alt_reject(is1)
if(xyzzyaaad1)call wout(' MOVE REJECTED')
endif
enddo
enddo
if(xyzzyaacn2)then
if(.not.use_altsamp)then
call eval_local_energy(is0,etot=xyzzyaacg2(xyzzyaaaa2-nequil))
else
call alt_halflogpdf(is0,xyzzyaabt2)
call simpletocomplex(is0,is1)
call eval_local_energy(is0,etot=xyzzyaacg2(xyzzyaaaa2-nequil))
call weight_config_vmc(is0,xyzzyaabu2)
xyzzyaacg2(xyzzyaaaa2-nequil)=xyzzyaabu2*xyzzyaacg2(xyzzyaaaa2-nequil)
xyzzyaach2(xyzzyaaaa2-nequil)=xyzzyaabu2
if(altsamp==1.and.have_ppots)grid_angles_valid(is0)=.false.
call complextosimple(is0,is1)
endif
endif
else
if(xyzzyaaad1)then
do xyzzyaaab2=1,netot
call wout('R_'//trim(i2s(xyzzyaaab2))//'_OLD: ',rele(1:dimensionality,&
&xyzzyaaab2))
if(xyzzyaaad2>1)call wout('S_OLD: '//trim(i2s(sele(xyzzyaaab2))))
enddo
endif
do xyzzyaaae2=1,xyzzyaaad2
call xyzzyaaal1(is0,is1,dt_array,sqrt_dt_array,dt_shift_array,xyzzyaaa&
&e2,xyzzyaacf2,xyzzyaaas2,xyzzyaabr2)
xyzzyaabf2=0.d0
do xyzzyaaab2=1,netot
call vector_difference(rele(1,xyzzyaaab2),xyzzyaacf2(1,xyzzyaaab2),xyz&
&zyaabs2)
xyzzyaabf2=xyzzyaabf2+sum(xyzzyaabs2**2)
enddo
if(xyzzyaaad1)then
do xyzzyaaab2=1,netot
if(xyzzyaaae2==1)then
call wout('R_'//trim(i2s(xyzzyaaab2))//'_NEW: ',xyzzyaacf2(1:dimension&
&ality,xyzzyaaab2))
else
call wout('S_'//trim(i2s(xyzzyaaab2))//'_NEW: '//trim(i2s(sele(xyzzyaa&
&ab2))))
endif
enddo
endif
if(altsamp==1.and.have_ppots)then
if(grid_angles_valid(is0))grid_angles_scr(:,:,:,is1)=grid_angles_scr(:&
&,:,:,is0)
grid_angles_valid(is1)=grid_angles_valid(is0)
endif
xyzzyaabi2=1.d0
xyzzyaack2=.true.
do xyzzyaaai2=1,xyzzyaaap2
if(xyzzyaaab1)then
xyzzyaaaq2=xyzzyaaai2
else
xyzzyaaaq2=0
endif
if(.not.use_altsamp)then
call wfn_ratio(is0,is1,xyzzyaaaq2,relprob=xyzzyaabh2,isnan=isnan,isinf&
&=isinf)
else
call wfn_altratio(is0,is1,xyzzyaaaq2,relprob=xyzzyaabh2,isnan=isnan,is&
&inf=isinf)
endif
if(xyzzyaaai2==1)xyzzyaabh2=xyzzyaabh2*xyzzyaabr2
prob=min(1.d0,xyzzyaabh2)
if(xyzzyaaad1)then
if(isnan)then
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaaai2))//&
&': zero [was Not a Number]')
elseif(isinf)then
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaaai2))//&
&': zero [it diverged]')
else
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaaai2))//&
&': ',prob)
endif
endif
if(isnan.or.isinf)then
prob=0.d0
xyzzyaabh2=0.d0
endif
xyzzyaabi2=xyzzyaabi2*xyzzyaabh2
if(.not.xyzzyaack2)cycle
xyzzyaaav2=0.d0
if(prob<1.d0)xyzzyaaav2=ranx()
xyzzyaack2=xyzzyaaav2<prob
if(.not.xyzzyaack2)then
if(xyzzyaabi2>0.d0)cycle
exit
endif
enddo
xyzzyaabi2=min(xyzzyaabi2,1.d0)
prob=xyzzyaabi2
xyzzyaabj2=1.d0-xyzzyaabi2
if(xyzzyaaad2==2)then
prob=0.5d0*prob
xyzzyaabj2=0.5d0*xyzzyaabj2
endif
xyzzyaabd2(1)=xyzzyaabd2(1)+xyzzyaabf2*prob
if(xyzzyaacn2.and.xyzzyaaac1)then
if(.not.use_altsamp)then
call eval_local_energy(is0,etot=xyzzyaabt2)
call eval_local_energy(is1,etot=xyzzyaabu2)
xyzzyaacg2(xyzzyaaaa2-nequil)=xyzzyaacg2(xyzzyaaaa2-nequil)+xyzzyaabt2&
&*xyzzyaabj2+xyzzyaabu2*prob
else
call alt_halflogpdf(is0,xyzzyaabv2)
call alt_halflogpdf(is1,xyzzyaabw2)
call simpletocomplex(is0,is1)
call eval_local_energy(is0,etot=xyzzyaabt2)
call eval_local_energy(is1,etot=xyzzyaabu2)
call weight_config_vmc(is0,xyzzyaabv2)
call weight_config_vmc(is1,xyzzyaabw2)
xyzzyaacg2(xyzzyaaaa2-nequil)=xyzzyaacg2(xyzzyaaaa2-nequil)+xyzzyaabt2&
&*xyzzyaabv2*xyzzyaabj2+xyzzyaabu2*xyzzyaabv2*prob
xyzzyaach2(xyzzyaaaa2-nequil)=xyzzyaach2(xyzzyaaaa2-nequil)+xyzzyaabv2&
&*xyzzyaabj2+xyzzyaabv2*prob
if(simplepdf==1.and.isitcomplex)then
call empty_scratch_wfn(is0)
call empty_scratch_wfn(is1)
endif
if(altsamp==1.and.have_ppots)grid_angles_valid(is0)=.false.
call complextosimple(is0,is1)
endif
endif
if(xyzzyaaae2==1)xyzzyaacc2(1)=xyzzyaacc2(1)+1.d0
if(xyzzyaack2)then
call accept_move(is0,is1)
if(use_altsamp)call alt_accept(is0,is1)
call dcopy(three_netot,xyzzyaacf2(1,1),1,rele(1,1),1)
sele(:)=xyzzyaaas2(:)
if(xyzzyaaae2==1)xyzzyaaca2(1)=xyzzyaaca2(1)+1.d0
if(xyzzyaaad1)call wout('MOVE ACCEPTED')
else
if(use_altsamp)call alt_reject(is1)
if(xyzzyaaad1)call wout('MOVE REJECTED')
endif
enddo
if(xyzzyaacn2.and..not.xyzzyaaac1)then
if(.not.use_altsamp)then
call eval_local_energy(is0,etot=xyzzyaacg2(xyzzyaaaa2-nequil))
else
call alt_halflogpdf(is0,xyzzyaabv2)
call simpletocomplex(is0,is1)
call eval_local_energy(is0,etot=xyzzyaabt2)
call weight_config_vmc(is0,xyzzyaabv2)
xyzzyaacg2(xyzzyaaaa2-nequil)=xyzzyaabt2*xyzzyaabv2
xyzzyaach2(xyzzyaaaa2-nequil)=xyzzyaabv2
if(altsamp==1.and.have_ppots)grid_angles_valid(is0)=.false.
call complextosimple(is0,is1)
endif
endif
endif
enddo
if(use_altsamp)call simpletocomplex
call timer('EQUILIBRATION',.false.)
if(vmc_cfg_by_cfg)deallocate(xyzzyaacf2,xyzzyaaas2)
contains
subroutine xyzzyaacp2
use slaarnabg, only : isperiodic,wigner_seitz_radius
implicit none
integer xyzzyaaaa3
real(dp) xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaaaf3,xyzzya&
&aag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3,xyzzyaaal3,xyzzyaaam&
&3,xyzzyaaan3,xyzzyaaao3
logical xyzzyaaap3
real(dp),parameter :: xyzzyaaaq3=0.199995d0,xyzzyaaar3=0.019995d0
real(dp),parameter :: xyzzyaaas3=0.1d0
real(dp),parameter :: xyzzyaaat3=4.d0,xyzzyaaau3=1.15d0
real(dp),parameter :: xyzzyaaav3=1.1d0,xyzzyaaaw3=1.01d0
real(dp),parameter :: xyzzyaaax3=10.d0
if(xyzzyaaaf1)call wout(' Iteration '//trim(i2s(xyzzyaaam2)))
xyzzyaaan3=dble(xyzzyaaam2)
xyzzyaaao3=dble(xyzzyaaan2)
xyzzyaaal3=(xyzzyaaan3*(xyzzyaaao3*xyzzyaaau3-xyzzyaaat3)-xyzzyaaao3*(&
&xyzzyaaau3-xyzzyaaat3))/(xyzzyaaan3*(xyzzyaaao3-1.d0))
xyzzyaaam3=(xyzzyaaan3*(xyzzyaaao3*xyzzyaaaw3-xyzzyaaav3)-xyzzyaaao3*(&
&xyzzyaaaw3-xyzzyaaav3))/(xyzzyaaan3*(xyzzyaaao3-1.d0))
do xyzzyaaaf2=1,xyzzyaaak2
if(opt_dtvmc_array(xyzzyaaaf2)==0)cycle
xyzzyaaab3=dt_array(xyzzyaaaf2)
xyzzyaaac3=xyzzyaaaw2(xyzzyaaaf2)
xyzzyaaad3=xyzzyaaax2(xyzzyaaaf2)
xyzzyaaae3=xyzzyaaay2(xyzzyaaaf2)
xyzzyaaaj3=xyzzyaabk2(xyzzyaaaf2)-.5d0
xyzzyaaak3=xyzzyaabc2(xyzzyaaaf2)
xyzzyaaaf3=xyzzyaabg2(xyzzyaaaf2)
xyzzyaaag3=xyzzyaaaz2(xyzzyaaaf2)
xyzzyaaah3=xyzzyaaba2(xyzzyaaaf2)
xyzzyaaai3=xyzzyaabb2(xyzzyaaaf2)
xyzzyaaaa3=xyzzyaaal2(xyzzyaaaf2)
if(xyzzyaaaf1)then
call wout('  For family #'//trim(i2s(xyzzyaaaf2))//':')
call wout('   Current dt = ',xyzzyaaab3,rfmt='(es12.4)')
call wout('   Measured D = ',xyzzyaaaf3,rfmt='(es12.4)')
call wout('   Measured a = ',(xyzzyaaaj3+.5d0)*100.d0,' %',rfmt='(f5.1&
&)')
endif
select case(xyzzyaaaa3)
case(0,1)
if(xyzzyaaaf1)call wout('   Target is to achieve a ~ 50%')
if(abs(xyzzyaaby2(xyzzyaaaf2))>abs(xyzzyaaaj3))then
if(xyzzyaaaf1)call wout('   Flagging as best dt so far')
xyzzyaabx2(xyzzyaaaf2)=xyzzyaaab3
xyzzyaaby2(xyzzyaaaf2)=xyzzyaaaj3
xyzzyaabz2(xyzzyaaaf2)=xyzzyaaaf3
endif
if(abs(xyzzyaaaj3)>=xyzzyaaaq3)then
if(xyzzyaaaf1)call wout('   a is far from 50%')
xyzzyaaad3=xyzzyaaab3
if(xyzzyaaaj3>0.d0)then
if(xyzzyaaaf1)call wout('   Next dt is an increase from current')
xyzzyaaab3=xyzzyaaab3*xyzzyaaal3
else
if(xyzzyaaaf1)call wout('   Next dt is a decrease from current')
xyzzyaaab3=xyzzyaaab3/xyzzyaaal3
endif
xyzzyaaak3=xyzzyaaaj3
xyzzyaaaa3=1
elseif(xyzzyaaaa3==0)then
xyzzyaaad3=xyzzyaaab3
if(xyzzyaaaj3>0.d0)then
if(xyzzyaaaf1)call wout('   Next dt is an increase from current')
xyzzyaaab3=xyzzyaaab3*2.d0
else
if(xyzzyaaaf1)call wout('   Next dt is a decrease from current')
xyzzyaaab3=xyzzyaaab3*.5d0
endif
xyzzyaaak3=xyzzyaaaj3
xyzzyaaaa3=1
else
if(abs(xyzzyaaaj3-xyzzyaaak3)>xyzzyaaar3)then
xyzzyaaae3=xyzzyaaab3
xyzzyaaab3=xyzzyaaab3*exp(-xyzzyaaaj3*(log(xyzzyaaab3)-log(xyzzyaaad3)&
&)/(xyzzyaaaj3-xyzzyaaak3))
if(xyzzyaaaf1)then
if(xyzzyaaab3>max(xyzzyaaad3,xyzzyaaae3).or.xyzzyaaab3<min(xyzzyaaad3,&
&xyzzyaaae3))then
call wout('   Next dt is a linear extrapolation of current and previou&
&s values')
else
call wout('   Next dt is a linear interpolation of current and previou&
&s values')
endif
endif
xyzzyaaad3=xyzzyaaae3
xyzzyaaak3=xyzzyaaaj3
else
xyzzyaaad3=xyzzyaaab3
if(xyzzyaaaj3>=0.d0)then
if(xyzzyaaaf1)call wout('   Next dt is a small increase from current')
xyzzyaaab3=xyzzyaaab3*xyzzyaaam3
else
if(xyzzyaaaf1)call wout('   Next dt is a small decrease from current')
xyzzyaaab3=xyzzyaaab3/xyzzyaaam3
endif
xyzzyaaak3=xyzzyaaaj3
endif
if(opt_dtvmc_array(xyzzyaaaf2)==2.and.abs(xyzzyaaaj3)<xyzzyaaas3)then
if(xyzzyaaaf1)call wout('   Switching to maximizing D in next iteratio&
&n')
xyzzyaaaa3=2
endif
endif
case default
if(xyzzyaaaf1)call wout('   Target is to maximize D')
if(abs(xyzzyaabz2(xyzzyaaaf2))<abs(xyzzyaaaf3))then
if(xyzzyaaaf1)call wout('   Flagging as best dt so far')
xyzzyaabx2(xyzzyaaaf2)=xyzzyaaab3
xyzzyaaby2(xyzzyaaaf2)=xyzzyaaaj3
xyzzyaabz2(xyzzyaaaf2)=xyzzyaaaf3
endif
select case(xyzzyaaaa3)
case(2)
if(xyzzyaaaf1)call wout('   Next dt is an increase from current (popul&
&ating bracket points)')
xyzzyaaad3=xyzzyaaab3
xyzzyaaah3=xyzzyaaaf3
xyzzyaaab3=xyzzyaaab3*xyzzyaaal3
xyzzyaaaa3=3
case(3)
if(xyzzyaaaf3>xyzzyaaah3.eqv.xyzzyaaab3>xyzzyaaad3)then
if(xyzzyaaaf1)call wout('   Next dt is an increase from current (popul&
&ating bracket points)')
xyzzyaaac3=xyzzyaaad3
xyzzyaaag3=xyzzyaaah3
xyzzyaaad3=xyzzyaaab3
xyzzyaaah3=xyzzyaaaf3
xyzzyaaab3=xyzzyaaad3*xyzzyaaal3
xyzzyaaaa3=4
else
if(xyzzyaaaf1)call wout('   Next dt is a decrease from previous (popul&
&ating bracket points)')
xyzzyaaae3=xyzzyaaab3
xyzzyaaai3=xyzzyaaaf3
xyzzyaaab3=xyzzyaaad3/xyzzyaaal3
xyzzyaaaa3=5
endif
case default
if(xyzzyaaaa3==4)then
xyzzyaaae3=xyzzyaaab3
xyzzyaaai3=xyzzyaaaf3
xyzzyaaaa3=6
elseif(xyzzyaaaa3==5)then
xyzzyaaac3=xyzzyaaab3
xyzzyaaag3=xyzzyaaaf3
xyzzyaaaa3=6
else
if(xyzzyaaab3<xyzzyaaac3)then
xyzzyaaae3=xyzzyaaad3
xyzzyaaai3=xyzzyaaah3
xyzzyaaad3=xyzzyaaac3
xyzzyaaah3=xyzzyaaag3
xyzzyaaac3=xyzzyaaab3
xyzzyaaag3=xyzzyaaaf3
elseif(xyzzyaaab3<xyzzyaaad3)then
if(xyzzyaaaf3>xyzzyaaah3)then
xyzzyaaae3=xyzzyaaad3
xyzzyaaai3=xyzzyaaah3
xyzzyaaad3=xyzzyaaab3
xyzzyaaah3=xyzzyaaaf3
else
xyzzyaaac3=xyzzyaaab3
xyzzyaaag3=xyzzyaaaf3
endif
elseif(xyzzyaaab3>xyzzyaaae3)then
xyzzyaaac3=xyzzyaaad3
xyzzyaaag3=xyzzyaaah3
xyzzyaaad3=xyzzyaaae3
xyzzyaaah3=xyzzyaaai3
xyzzyaaae3=xyzzyaaab3
xyzzyaaai3=xyzzyaaaf3
elseif(xyzzyaaab3>xyzzyaaad3)then
if(xyzzyaaaf3>xyzzyaaah3)then
xyzzyaaac3=xyzzyaaad3
xyzzyaaag3=xyzzyaaah3
xyzzyaaad3=xyzzyaaab3
xyzzyaaah3=xyzzyaaaf3
else
xyzzyaaae3=xyzzyaaab3
xyzzyaaai3=xyzzyaaaf3
endif
endif
endif
if(xyzzyaaag3<xyzzyaaah3)then
if(xyzzyaaah3<xyzzyaaai3)then
if(xyzzyaaaf1)call wout('   Next dt is an increase from right bracket &
&point (bracketing)')
xyzzyaaab3=xyzzyaaae3*xyzzyaaal3
else
if(xyzzyaaaf1)call wout('   Next dt is a parabolic interpolation')
call parabolic_min(log(xyzzyaaac3),log(xyzzyaaad3),log(xyzzyaaae3),xyz&
&zyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaab3,xyzzyaaaf3,xyzzyaaap3)
if(xyzzyaaap3)then
xyzzyaaab3=xyzzyaaad3
else
xyzzyaaab3=exp(xyzzyaaab3)
endif
if(xyzzyaaae3-xyzzyaaac3<0.1d0*xyzzyaaad3)then
if(xyzzyaaaf1)call wout('   Will restart maximization from there')
xyzzyaaaa3=2
endif
endif
else
if(xyzzyaaah3<xyzzyaaai3)then
if(xyzzyaaai3>xyzzyaaag3)then
if(xyzzyaaaf1)call wout('   Next dt is an increase from right bracket &
&point (escape from minimum)')
xyzzyaaab3=xyzzyaaae3*xyzzyaaal3
xyzzyaaac3=xyzzyaaad3
xyzzyaaag3=xyzzyaaah3
xyzzyaaad3=xyzzyaaae3
xyzzyaaah3=xyzzyaaai3
xyzzyaaaa3=4
else
if(xyzzyaaaf1)call wout('   Next dt is a decrease from left bracket po&
&int (escape from minimum)')
xyzzyaaab3=xyzzyaaac3/xyzzyaaal3
xyzzyaaae3=xyzzyaaad3
xyzzyaaai3=xyzzyaaah3
xyzzyaaad3=xyzzyaaac3
xyzzyaaah3=xyzzyaaag3
xyzzyaaaa3=5
endif
else
if(xyzzyaaaf1)call wout('   Next dt is a decrease from left bracket po&
&int (bracketing)')
xyzzyaaab3=xyzzyaaac3/xyzzyaaal3
endif
endif
end select
end select
if(isperiodic)then
if(xyzzyaaab3>xyzzyaaax3*wigner_seitz_radius**2)then
if(xyzzyaaaf1)call wout('   dt over limit, setting to limit')
xyzzyaaab3=xyzzyaaax3*wigner_seitz_radius**2
endif
endif
dt_array(xyzzyaaaf2)=xyzzyaaab3
xyzzyaaaw2(xyzzyaaaf2)=xyzzyaaac3
xyzzyaaax2(xyzzyaaaf2)=xyzzyaaad3
xyzzyaaay2(xyzzyaaaf2)=xyzzyaaae3
xyzzyaabc2(xyzzyaaaf2)=xyzzyaaak3
xyzzyaaaz2(xyzzyaaaf2)=xyzzyaaag3
xyzzyaaba2(xyzzyaaaf2)=xyzzyaaah3
xyzzyaabb2(xyzzyaaaf2)=xyzzyaaai3
xyzzyaaal2(xyzzyaaaf2)=xyzzyaaaa3
enddo
if(vmc_cfg_by_cfg)then
do xyzzyaaaf2=1,no_difftypes
dt_array(xyzzyaaaf2)=xyzzyaaab3*difftype_mass(1)/difftype_mass(xyzzyaa&
&af2)
enddo
endif
end subroutine xyzzyaacp2
subroutine xyzzyaacq2(xyzzyaabn2,xyzzyaabm2,xyzzyaabl2,xyzzyaabk2,vmc_&
&cfg_by_cfg,corper)
implicit none
logical,intent(in) :: vmc_cfg_by_cfg
real(dp),intent(in) :: xyzzyaabn2,xyzzyaabm2,xyzzyaabl2,xyzzyaabk2(:)
integer,intent(out) :: corper
integer xyzzyaaaa4,xyzzyaaab4
real(dp) xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaaaf4
if(.not.vmc_cfg_by_cfg)then
xyzzyaaaf4=1.d0
do xyzzyaaaa4=1,no_difftypes
xyzzyaaaf4=xyzzyaaaf4*(1.d0-xyzzyaabk2(xyzzyaaaa4))**sum(nele(:),mask=&
&which_difftype==xyzzyaaaa4)
enddo
else
xyzzyaaaf4=1.d0-xyzzyaabk2(1)
endif
if(xyzzyaaag1)call wout(' Prob. of not moving = ',xyzzyaaaf4)
corper=1
xyzzyaaad4=xyzzyaacr2(xyzzyaabn2,xyzzyaabm2,xyzzyaabl2,xyzzyaaaf4,vmc_&
&cfg_by_cfg,corper)
if(xyzzyaaag1)call wout(' At p='//trim(i2s(corper))//' Eff(p)=',xyzzya&
&aad4)
if(xyzzyaaad4==0.d0)return
xyzzyaaae4=xyzzyaaad4
xyzzyaaab4=1
do
corper=corper+1
xyzzyaaac4=xyzzyaacr2(xyzzyaabn2,xyzzyaabm2,xyzzyaabl2,xyzzyaaaf4,vmc_&
&cfg_by_cfg,corper)
if(xyzzyaaag1)call wout(' At p='//trim(i2s(corper))//' Eff(p)=',xyzzya&
&aac4)
if(xyzzyaaac4>xyzzyaaae4)then
xyzzyaaae4=xyzzyaaac4
xyzzyaaab4=corper
elseif(xyzzyaaac4<=xyzzyaaad4.and.corper>=100)then
corper=xyzzyaaab4
exit
endif
xyzzyaaad4=xyzzyaaac4
enddo
end subroutine xyzzyaacq2
real(dp) function xyzzyaacr2(xyzzyaabn2,xyzzyaabm2,xyzzyaabl2,prob_nom&
&ove,vmc_cfg_by_cfg,corper)
implicit none
integer,intent(in) :: corper
real(dp),intent(in) :: xyzzyaabn2,xyzzyaabm2,xyzzyaabl2,prob_nomove
logical,intent(in) :: vmc_cfg_by_cfg
real(dp) xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5
xyzzyaaae5=max(1.d0,xyzzyaabn2)
xyzzyaaaa5=dble(corper)*xyzzyaabl2
if(vmc_cfg_by_cfg.and.xyzzyaaac1)then
xyzzyaaab5=(2.d0-prob_nomove**(corper-1))*xyzzyaabm2
else
xyzzyaaab5=(1.d0-prob_nomove**corper)*xyzzyaabm2/(1.d0-prob_nomove)
endif
xyzzyaaac5=1.d0+2.d0*((xyzzyaaae5-1.d0)**corper)/((xyzzyaaae5+1)**corp&
&er-(xyzzyaaae5-1)**corper)
xyzzyaaad5=(xyzzyaaaa5+xyzzyaaab5)*xyzzyaaac5
if(xyzzyaaad5==0.d0)then
xyzzyaacr2=0.d0
else
xyzzyaacr2=1.d0/xyzzyaaad5
endif
end function xyzzyaacr2
end subroutine equilibration
subroutine xyzzyaaak1(ii,is,js,dt,sqrt_dt,dt_shift,iattempt,rnew,snew,&
&db_factor)
use slaarnach
use slaarnabg,      only : nitot,isperiodic,dimensionality
use slaarnacc,only : ranx
use store,         only : on_top_ii,on_top_jj
use slaarnacs,     only : define_config_oneelec,define_config_twoelec
implicit none
integer,intent(in) :: ii,is,js,iattempt
integer,intent(out) :: snew
real(dp),intent(in) :: dt,sqrt_dt,dt_shift
real(dp),intent(out) :: rnew(3),db_factor
integer xyzzyaaaa6
real(dp) xyzzyaaab6(3),xyzzyaaac6,xyzzyaaad6
logical xyzzyaaae6,xyzzyaaaf6
logical,save :: xyzzyaaag6=.false.,xyzzyaaah6=.false.,xyzzyaaai6=.true&
&.
if(xyzzyaaai6)then
if(nitot>1.and..not.isperiodic.and.vmc_ionjump>0.d0)xyzzyaaag6=.true.
if(xyzzyaaaa1.and.nitot>0)xyzzyaaah6=.true.
xyzzyaaai6=.false.
endif
xyzzyaaaf6=.false.
xyzzyaaaa6=1
xyzzyaaab6(:)=0.d0
select case(iattempt)
case(1)
xyzzyaaae6=.false.
if(xyzzyaaag6)xyzzyaaae6=ranx()<vmc_ionjump
if(xyzzyaaae6)then
call xyzzyaaas1(rele_scr(:,ii,is),xyzzyaaab6)
else
if(xyzzyaaah6)then
call get_eivecs(is)
call xyzzyaaaq1(dimensionality,dt,sqrt_dt,eivecs_scr(1,1,ii,is),xyzzya&
&aab6,xyzzyaaac6)
xyzzyaaaf6=.true.
else
call xyzzyaaam1(dimensionality,sqrt_dt,dt_shift,xyzzyaaab6)
xyzzyaaac6=1.d0
endif
endif
case(2)
call xyzzyaaat1(xyzzyaaaa6)
end select
rnew=rele_scr(:,ii,is)+xyzzyaaab6(:)
if(xyzzyaaaa6==-1)then
snew=3-sele_scr(ii,is)
else
snew=sele_scr(ii,is)
endif
if(ii/=on_top_ii)then
call define_config_oneelec(ii,is,js,rnew,snew)
else
call define_config_twoelec(ii,on_top_jj,is,js,rnew,snew,rnew,snew)
endif
db_factor=1.d0
if(xyzzyaaaf6)then
if(xyzzyaaac6==0.d0)then
db_factor=0.d0
else
xyzzyaaab6=-xyzzyaaab6
if(xyzzyaaah6)then
call get_eivecs1_ch(ii,js)
call xyzzyaaar1(dimensionality,dt,eivecs1_chscr(1,1,js),xyzzyaaab6,xyz&
&zyaaad6)
else
xyzzyaaad6=1.d0
endif
if(xyzzyaaad6==0.d0)then
db_factor=0.d0
else
db_factor=xyzzyaaad6/xyzzyaaac6
endif
endif
endif
end subroutine xyzzyaaak1
subroutine xyzzyaaal1(is,js,dt_array,sqrt_dt_array,dt_shift_array,iatt&
&empt,rele_new,sele_new,db_factor)
use slaarnach
use slaarnaaq,        only : particle_is_fixed,fixed_particle
use slaarnabg,      only : nitot,isperiodic,dimensionality
use slaarnacc,only : ranx
use store,         only : netot,which_difftype,no_difftypes,which_spin&
&,on_top_ii,on_top_jj
use slaarnacs,     only : define_config
implicit none
integer,intent(in) :: is,js,iattempt
integer,intent(out) :: sele_new(netot)
real(dp),intent(in) :: dt_array(no_difftypes),sqrt_dt_array(no_difftyp&
&es),dt_shift_array(no_difftypes)
real(dp),intent(out) :: rele_new(3,netot),db_factor
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7
real(dp) xyzzyaaad7(3),xyzzyaaae7(netot),xyzzyaaaf7
logical xyzzyaaag7,xyzzyaaah7(netot),xyzzyaaai7
logical,save :: xyzzyaaaj7=.false.,xyzzyaaak7=.false.,xyzzyaaal7=.true&
&.
if(xyzzyaaal7)then
if(nitot>1.and..not.isperiodic.and.vmc_ionjump>0.d0)xyzzyaaaj7=.true.
if(xyzzyaaaa1.and.nitot>0)xyzzyaaak7=.true.
xyzzyaaal7=.false.
endif
xyzzyaaah7(:)=.false.
xyzzyaaai7=.false.
if(xyzzyaaak7)call get_eivecs(is)
do xyzzyaaaa7=1,netot
xyzzyaaac7=1
xyzzyaaad7(:)=0.d0
if((particle_is_fixed.and.fixed_particle==xyzzyaaaa7).or.xyzzyaaaa7==o&
&n_top_jj)cycle
select case(iattempt)
case(1)
xyzzyaaag7=.false.
if(xyzzyaaaj7)xyzzyaaag7=ranx()<vmc_ionjump
if(xyzzyaaag7)then
call xyzzyaaas1(rele_scr(:,xyzzyaaaa7,is),xyzzyaaad7)
else
xyzzyaaab7=which_difftype(which_spin(xyzzyaaaa7))
if(xyzzyaaak7)then
call xyzzyaaaq1(dimensionality,dt_array(xyzzyaaab7),sqrt_dt_array(xyzz&
&yaaab7),eivecs_scr(1,1,xyzzyaaaa7,is),xyzzyaaad7,xyzzyaaae7(xyzzyaaaa&
&7))
xyzzyaaah7(xyzzyaaaa7)=.true.
xyzzyaaai7=.true.
else
call xyzzyaaam1(dimensionality,sqrt_dt_array(xyzzyaaab7),dt_shift_arra&
&y(xyzzyaaab7),xyzzyaaad7)
xyzzyaaae7(xyzzyaaaa7)=1.d0
endif
endif
case(2)
call xyzzyaaat1(xyzzyaaac7)
end select
rele_new(:,xyzzyaaaa7)=rele_scr(:,xyzzyaaaa7,is)+xyzzyaaad7(:)
if(xyzzyaaac7==-1)then
sele_new(xyzzyaaaa7)=3-sele_scr(xyzzyaaaa7,is)
else
sele_new(xyzzyaaaa7)=sele_scr(xyzzyaaaa7,is)
endif
if(xyzzyaaaa7==on_top_ii)then
rele_new(:,on_top_jj)=rele_new(:,xyzzyaaaa7)
sele_new(on_top_jj)=sele_new(xyzzyaaaa7)
endif
enddo
call define_config(js,rele_new,sele_new)
db_factor=1.d0
if(xyzzyaaai7)then
if(xyzzyaaak7)call get_eivecs(js)
do xyzzyaaaa7=1,netot
if(.not.xyzzyaaah7(xyzzyaaaa7))cycle
if(xyzzyaaae7(xyzzyaaaa7)==0.d0)then
db_factor=0.d0
exit
endif
if(xyzzyaaak7)then
xyzzyaaad7=rele_scr(1:3,xyzzyaaaa7,is)-rele_scr(1:3,xyzzyaaaa7,js)
xyzzyaaab7=which_difftype(which_spin(xyzzyaaaa7))
call xyzzyaaar1(dimensionality,dt_array(xyzzyaaab7),eivecs_scr(1,1,xyz&
&zyaaaa7,js),xyzzyaaad7,xyzzyaaaf7)
else
xyzzyaaaf7=1.d0
endif
if(xyzzyaaaf7==0.d0)then
db_factor=0.d0
exit
else
db_factor=db_factor*xyzzyaaaf7/xyzzyaaae7(xyzzyaaaa7)
endif
enddo
endif
end subroutine xyzzyaaal1
subroutine xyzzyaaam1(d,sqrt_dt,dt_shift,dr)
use slaarnacc,only : ranx_gaussian
implicit none
integer,intent(in) :: d
real(dp),intent(in) :: sqrt_dt,dt_shift
real(dp),intent(inout) :: dr(*)
integer xyzzyaaaa8
real(dp) xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaaaf8
if(dt_shift<=0.d0)then
do xyzzyaaaa8=1,d
dr(xyzzyaaaa8)=ranx_gaussian(sqrt_dt)
enddo
else
xyzzyaaab8=ranx_gaussian(sqrt_dt)+sqrt_dt*dt_shift
select case(d)
case(1)
dr(1)=xyzzyaaab8
case(2)
xyzzyaaaf8=0.d0
do while(xyzzyaaaf8==0.d0)
xyzzyaaac8=ranx_gaussian(1.d0)
xyzzyaaad8=ranx_gaussian(1.d0)
xyzzyaaaf8=sqrt(xyzzyaaac8*xyzzyaaac8+xyzzyaaad8*xyzzyaaad8)
enddo
xyzzyaaac8=xyzzyaaac8/xyzzyaaaf8
xyzzyaaad8=xyzzyaaad8/xyzzyaaaf8
dr(1)=xyzzyaaab8*xyzzyaaac8
dr(2)=xyzzyaaab8*xyzzyaaad8
case(3)
xyzzyaaaf8=0.d0
do while(xyzzyaaaf8==0.d0)
xyzzyaaac8=ranx_gaussian(1.d0)
xyzzyaaad8=ranx_gaussian(1.d0)
xyzzyaaae8=ranx_gaussian(1.d0)
xyzzyaaaf8=sqrt(xyzzyaaac8*xyzzyaaac8+xyzzyaaad8*xyzzyaaad8+xyzzyaaae8&
&*xyzzyaaae8)
enddo
xyzzyaaac8=xyzzyaaac8/xyzzyaaaf8
xyzzyaaad8=xyzzyaaad8/xyzzyaaaf8
xyzzyaaae8=xyzzyaaae8/xyzzyaaaf8
dr(1)=xyzzyaaab8*xyzzyaaac8
dr(2)=xyzzyaaab8*xyzzyaaad8
dr(3)=xyzzyaaab8*xyzzyaaae8
end select
endif
end subroutine xyzzyaaam1
subroutine xyzzyaaan1(d,dt,dr,tprob)
use slaarnaag,only : twopi
implicit none
integer,intent(in) :: d
real(dp),intent(in) :: dt,dr(*)
real(dp),intent(out) :: tprob
real(dp) xyzzyaaaa9
xyzzyaaaa9=sum((dr(1:d))**2)
tprob=exp(-.5d0*xyzzyaaaa9/dt)/(twopi*dt)**(.5d0*dble(d))
end subroutine xyzzyaaan1
subroutine xyzzyaaao1(r,dr)
use slaarnaag,     only : two_over_root_pi
use slaarnacc,only : ranx
implicit none
real(dp),intent(in) :: r
real(dp),intent(inout) :: dr
dr=two_over_root_pi*r*sqrt(-log(1.d0-ranx()))-r
end subroutine xyzzyaaao1
subroutine xyzzyaaap1(r,dr,tprob)
use slaarnaag, only : pi_over_four
implicit none
real(dp),intent(in) :: r,dr
real(dp),intent(out) :: tprob
real(dp) xyzzyaaaa11,xyzzyaaab11
xyzzyaaab11=abs(r+dr)
xyzzyaaaa11=pi_over_four/(r*r)
tprob=(xyzzyaaaa11+xyzzyaaaa11)*xyzzyaaab11*exp(-xyzzyaaaa11*xyzzyaaab&
&11*xyzzyaaab11)
end subroutine xyzzyaaap1
subroutine xyzzyaaaq1(d,dt,sqrt_dt,eivecs1,dr,tprob)
use slaarnabg, only : nitot
implicit none
integer,intent(in) :: d
real(dp),intent(in) :: dt,sqrt_dt,eivecs1(4,nitot)
real(dp),intent(inout) :: dr(3),tprob
integer xyzzyaaaa12,xyzzyaaab12
real(dp) xyzzyaaac12,xyzzyaaad12(3,3),xyzzyaaae12(3),xyzzyaaaf12
xyzzyaaaa12=minloc(eivecs1(4,1:nitot),1)
xyzzyaaac12=eivecs1(4,xyzzyaaaa12)
if(xyzzyaaac12==0.d0)then
call xyzzyaaam1(d,sqrt_dt,-1.d0,dr(1))
call xyzzyaaan1(d,dt,dr(1),tprob)
return
endif
xyzzyaaad12(1:d,1)=eivecs1(1:d,xyzzyaaaa12)/xyzzyaaac12
if(d==2)then
xyzzyaaad12(1:2,2)=(/-xyzzyaaad12(2,1),xyzzyaaad12(1,1)/)
else
if(any(xyzzyaaad12(1:2,1)/=0.d0))then
xyzzyaaad12(1:3,2)=(/-xyzzyaaad12(2,1),xyzzyaaad12(1,1),0.d0/)
else
xyzzyaaad12(1:3,2)=(/0.d0,-xyzzyaaad12(3,1),xyzzyaaad12(2,1)/)
endif
xyzzyaaad12(1:3,2)=xyzzyaaad12(1:3,2)/sqrt(sum(xyzzyaaad12(1:3,2)**2))
xyzzyaaad12(1,3)=xyzzyaaad12(2,1)*xyzzyaaad12(3,2)-xyzzyaaad12(3,1)*xy&
&zzyaaad12(2,2)
xyzzyaaad12(2,3)=xyzzyaaad12(3,1)*xyzzyaaad12(1,2)-xyzzyaaad12(1,1)*xy&
&zzyaaad12(3,2)
xyzzyaaad12(3,3)=xyzzyaaad12(1,1)*xyzzyaaad12(2,2)-xyzzyaaad12(2,1)*xy&
&zzyaaad12(1,2)
endif
call xyzzyaaao1(xyzzyaaac12,xyzzyaaae12(1))
call xyzzyaaap1(xyzzyaaac12,xyzzyaaae12(1),tprob)
if(d>1)then
call xyzzyaaam1(d-1,sqrt_dt,-1.d0,xyzzyaaae12(2))
call xyzzyaaan1(d-1,dt,xyzzyaaae12(2),xyzzyaaaf12)
tprob=tprob*xyzzyaaaf12
endif
dr(1:d)=0.d0
do xyzzyaaab12=1,d
dr(1:d)=dr(1:d)+xyzzyaaae12(xyzzyaaab12)*xyzzyaaad12(1:d,xyzzyaaab12)
enddo
end subroutine xyzzyaaaq1
subroutine xyzzyaaar1(d,dt,eivecs1,dr,tprob)
use slaarnabg, only : nitot
implicit none
integer,intent(in) :: d
real(dp),intent(in) :: dt,dr(3),eivecs1(4,nitot)
real(dp),intent(out) :: tprob
integer xyzzyaaaa13
real(dp) xyzzyaaab13,xyzzyaaac13(3),xyzzyaaad13(3),xyzzyaaae13
xyzzyaaaa13=minloc(eivecs1(4,1:nitot),1)
xyzzyaaab13=eivecs1(4,xyzzyaaaa13)
if(xyzzyaaab13==0.d0)then
call xyzzyaaan1(d,dt,dr(1),tprob)
return
endif
xyzzyaaac13(1:d)=eivecs1(1:d,xyzzyaaaa13)/xyzzyaaab13
xyzzyaaad13(1)=sum(dr(1:d)*xyzzyaaac13(1:d))
if(d>1)then
xyzzyaaad13(2)=sqrt(sum((dr(1:d)-xyzzyaaac13(1:d)*xyzzyaaad13(1))**2))
xyzzyaaad13(3)=0.d0
endif
call xyzzyaaap1(xyzzyaaab13,xyzzyaaad13(1),tprob)
if(d>1)then
call xyzzyaaan1(d-1,dt,xyzzyaaad13(2),xyzzyaaae13)
tprob=tprob*xyzzyaaae13
endif
end subroutine xyzzyaaar1
subroutine xyzzyaaas1(rold,dr)
use slaarnabg,      only : nitot,rion
use slaarnacc,only : ranx
real(dp),intent(in) :: rold(3)
real(dp),intent(out) :: dr(3)
integer xyzzyaaaa14,xyzzyaaab14
real(dp) xyzzyaaac14(3),xyzzyaaad14(3)
xyzzyaaaa14=int(ranx()*dble(nitot))+1
xyzzyaaab14=int(ranx()*dble(nitot-1))+1
if(xyzzyaaab14>=xyzzyaaaa14)xyzzyaaab14=xyzzyaaab14+1
xyzzyaaac14(:)=rion(:,xyzzyaaaa14)
xyzzyaaad14(:)=rion(:,xyzzyaaab14)
if(sum((rold(:)-xyzzyaaac14(:))**2)<sum((rold(:)-xyzzyaaad14(:))**2))t&
&hen
dr(:)=xyzzyaaad14(:)-xyzzyaaac14(:)
else
dr(:)=xyzzyaaac14(:)-xyzzyaaad14(:)
endif
end subroutine xyzzyaaas1
subroutine xyzzyaaat1(dspin)
use slaarnacc,only : ranx
implicit none
integer,intent(out) :: dspin
if(ranx()<0.5d0)then
dspin=-1
else
dspin=1
endif
end subroutine xyzzyaaat1
subroutine xyzzyaaau1(iblock,nmove,nblock,no_data,ave,std_err,asratio_&
&lev,corr_tau,corr_tau_err,vartemp_tot,vartemp_err,diff_const,r2maxn,f&
&accspinflipn,elapsed_block_time,nmove_total)
use slaarnaam
use store
use format_utils,only : wout,r2s,i2s
use slaarnabk,only : mpc_hartree
use slaarnaca,       only : have_ppots,have_veep
use slaarnaaq,      only : int_sf,eval_int_sf,eval_contact_den
use slaarnabg,    only : isperiodic,npcells,periodicity,model_system,n&
&itot,ignore_ionic_interactions
use parallel,    only : nnodes
use slaarnace,  only : relativistic,eval_maspol,eval_masvel,eval_darwi&
&nen,eval_darwinee,eval_retard
implicit none
integer,intent(in) :: no_data,iblock,nmove,nblock,nmove_total
real(sp),intent(in) :: elapsed_block_time
real(dp),intent(in) :: asratio_lev(:,:),corr_tau,corr_tau_err,vartemp_&
&tot,vartemp_err,diff_const(no_difftypes),r2maxn,faccspinflipn
real(dp),intent(inout) :: ave(no_data),std_err(no_data)
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16
real(dp) xyzzyaaae16,xyzzyaaaf16,xyzzyaaag16,xyzzyaaah16,xyzzyaaai16,x&
&yzzyaaaj16,xyzzyaaak16
character(11) corr_tau_str,corr_tau_err_str,diff_const_str,eff_str
if(esupercell)then
ave(1:no_data)=ave(1:no_data)*real(npcells,dp)
std_err(1:no_data)=std_err(1:no_data)*real(npcells,dp)
xyzzyaaae16=constant_energy
if(finite_size_corr.and.isvmc)then
xyzzyaaah16=ke_corr
xyzzyaaai16=xc_corr
endif
xyzzyaaak16=real(npcells,dp)
elseif(isperiodic.and.model_system)then
xyzzyaaae16=constant_energy*inv_netot
if(finite_size_corr.and.isvmc)then
xyzzyaaah16=ke_corr*inv_netot
xyzzyaaai16=xc_corr*inv_netot
endif
xyzzyaaak16=1.d0
else
xyzzyaaae16=constant_energy/real(npcells,dp)
if(finite_size_corr.and.isvmc)then
xyzzyaaah16=ke_corr/real(npcells,dp)
xyzzyaaai16=xc_corr/real(npcells,dp)
endif
xyzzyaaak16=1.d0
endif
xyzzyaaac16=size(asratio_lev,1)
xyzzyaaaa16=size(asratio_lev,2)
call wout()
call wout(repeat('=',73))
call wout('In block : '//trim(i2s(iblock)))
call wout()
corr_tau_str=adjustl(r2s(corr_tau,'(es12.4)'))
corr_tau_err_str=adjustl(r2s(corr_tau_err,'(es12.4)'))
if(vartemp_tot>0.d0.and.elapsed_block_time>0.)then
xyzzyaaaj16=dble(nmove_total)/(dble(elapsed_block_time)*vartemp_tot*ma&
&x(corr_tau,1.d0))
eff_str=adjustl(r2s(xyzzyaaaj16,'(es12.4)'))
else
eff_str='Infinity'
endif
if(xyzzyaaaa16==1)then
diff_const_str=r2s(diff_const(1),'(es12.4)')
if(xyzzyaaac16==1)then
call wout('Acceptance ratio         (%)  = ',asratio_lev(1,1)*100.d0,r&
&fmt='(f8.4)')
if(noncoll_spin)call wout('Acc. rat. for spin flips (%)  = ',faccspinf&
&lipn*100.d0,rfmt='(f8.4)')
call wout('Diffusion constant  (Bohr^2)  = '//trim(diff_const_str))
call wout('Correlation time     (steps)  = '//trim(corr_tau_str)//' +-&
& '//trim(corr_tau_err_str))
call wout('Efficiency      (au^-2 s^-1)  = '//trim(eff_str))
if(nnodes>1)then
call wout('No. of VMC steps per process  = '//i2s(nmove))
else
call wout('Number of VMC steps           = '//i2s(nmove))
endif
else
call wout('Acceptance ratio <level 1>        (%)  = ',asratio_lev(1,1)&
&*100.d0,rfmt='(f8.4)')
do xyzzyaaad16=2,xyzzyaaac16
call wout('Acceptance ratio <levels 1-'//trim(i2s(xyzzyaaad16))//'>   &
&  (%)  = ',asratio_lev(xyzzyaaad16,1)*100.d0,rfmt='(f8.4)')
enddo
if(noncoll_spin)call wout('Acceptance ratio for spin flips (%)    = ',&
&faccspinflipn*100.d0,rfmt='(f8.4)')
call wout('Diffusion constant           (Bohr^2)  = '//trim(diff_const&
&_str))
call wout('Correlation time              (steps)  = '//trim(corr_tau_s&
&tr)//' +- '//trim(corr_tau_err_str))
call wout('Efficiency               (au^-2 s^-1)  = '//trim(eff_str))
if(nnodes>1)then
call wout('No. of VMC steps per process           = '//i2s(nmove))
else
call wout('Number of VMC steps                    = '//i2s(nmove))
endif
endif
else
if(xyzzyaaac16==1)then
do xyzzyaaab16=1,xyzzyaaaa16
call wout('Acceptance ratio #'//trim(i2s(xyzzyaaab16))//'         (%) &
& = ',asratio_lev(1,xyzzyaaab16)*100.d0,rfmt='(f8.4)')
enddo
if(noncoll_spin)call wout('Acc. ratio for spin flips   (%)  = ',faccsp&
&inflipn*100.d0,rfmt='(f8.4)')
do xyzzyaaab16=1,xyzzyaaaa16
diff_const_str=r2s(diff_const(xyzzyaaab16),'(es12.4)')
call wout('Diffusion constant #'//trim(i2s(xyzzyaaab16))//'  (Bohr^2) &
& = '//trim(diff_const_str))
enddo
call wout('Correlation time        (steps)  = '//trim(corr_tau_str)//'&
& +- '//trim(corr_tau_err_str))
call wout('Efficiency         (au^-2 s^-1)  = '//trim(eff_str))
if(nnodes>1)then
call wout('No. of VMC steps per process     = '//i2s(nmove))
else
call wout('Number of VMC steps              = '//i2s(nmove))
endif
else
do xyzzyaaab16=1,xyzzyaaaa16
call wout('Acceptance ratio #'//trim(i2s(xyzzyaaab16))//' <level 1>   &
&     (%)  = ',asratio_lev(1,xyzzyaaab16)*100.d0,rfmt='(f8.4)')
enddo
do xyzzyaaad16=2,xyzzyaaac16
do xyzzyaaab16=1,xyzzyaaaa16
call wout('Acceptance ratio #'//trim(i2s(xyzzyaaab16))//' <levels 1-'/&
&/trim(i2s(xyzzyaaad16))//'>     (%)  = ',asratio_lev(xyzzyaaad16,xyzz&
&yaaab16)*100.d0,rfmt='(f8.4)')
enddo
enddo
if(noncoll_spin)call wout('Acceptance ratio for spin flips      (%)  =&
& ',faccspinflipn*100.d0,rfmt='(f8.4)')
do xyzzyaaab16=1,xyzzyaaaa16
diff_const_str=r2s(diff_const(xyzzyaaab16),'(es12.4)')
call wout('Diffusion constant #'//trim(i2s(xyzzyaaab16))//'           &
&(Bohr^2)  = '//trim(diff_const_str))
enddo
call wout('Correlation time                 (steps)  = '//trim(corr_ta&
&u_str)//' +- '//trim(corr_tau_err_str))
call wout('Efficiency                  (au^-2 s^-1)  = '//trim(eff_str&
&))
if(nnodes>1)then
call wout('No. of VMC steps per process              = '//i2s(nmove))
else
call wout('Number of VMC steps                       = '//i2s(nmove))
endif
endif
endif
call wout()
if(isperiodic)then
if(model_system)then
call wout(' Block average energies (au per particle)')
else
if(esupercell.or.npcells==1)then
call wout(' Block average energies (au per simulation cell)')
else
call wout(' Block average energies (au per primitive cell)')
endif
endif
else
call wout(' Block average energies (au)')
endif
call wout()
select case(trim(interaction))
case('mpc')
call xyzzyaaav1('Total energy',ave(i_eloc_mpc),std_err(i_eloc_mpc))
case('ewald_mpc')
call xyzzyaaav1('Total energy with Ewald',ave(i_eloc_def),std_err(i_el&
&oc_def))
call xyzzyaaav1('Total energy with MPC',ave(i_eloc_mpc),std_err(i_eloc&
&_mpc))
case('mpc_ewald')
call xyzzyaaav1('Total energy with MPC',ave(i_eloc_mpc),std_err(i_eloc&
&_mpc))
call xyzzyaaav1('Total energy with Ewald',ave(i_eloc_def),std_err(i_el&
&oc_def))
case default
call xyzzyaaav1('Total energy',ave(i_eloc_def),std_err(i_eloc_def))
end select
if(relativistic)then
call xyzzyaaav1('Total energy (relativistic)',ave(i_eloc_def)+ave(i_re&
&ltot),sqrt(std_err(i_eloc_def)**2+std_err(i_reltot)**2))
endif
call xyzzyaaav1('Kinetic energy KEI (used in Total)',ave(i_kei),std_er&
&r(i_kei))
call xyzzyaaav1('Kinetic energy TI',ave(i_ti),std_err(i_ti))
call xyzzyaaav1('Kinetic energy FISQ',ave(i_fisq),std_err(i_fisq))
call xyzzyaaav1('Potential energy',ave(i_pei),std_err(i_pei))
select case(trim(interaction))
case('coulomb')
call xyzzyaaav1('e-e interaction',ave(i_pote),std_err(i_pote))
case('ewald')
call xyzzyaaav1('Ewald e-e interaction',ave(i_pote),std_err(i_pote))
case('mpc')
call xyzzyaaav1('MPC e-e interaction',ave(i_mpc),std_err(i_mpc))
call xyzzyaaav1('Short range part of MPC',ave(i_short),std_err(i_short&
&))
call xyzzyaaav1('Long range part of MPC',ave(i_long),std_err(i_long))
if(hartree_xc.and.periodicity==3)then
call xyzzyaaav1('MPC estimate of Hartree energy',mpc_hartree*xyzzyaaak&
&16)
call xyzzyaaav1('MPC estimate of XC energy',ave(i_mpc)-mpc_hartree*xyz&
&zyaaak16)
endif
case('ewald_mpc')
call xyzzyaaav1('Ewald e-e interaction',ave(i_pote),std_err(i_pote))
call xyzzyaaav1('MPC e-e interaction',ave(i_mpc),std_err(i_mpc))
call xyzzyaaav1('Short range part of MPC',ave(i_short),std_err(i_short&
&))
call xyzzyaaav1('Long range part of MPC',ave(i_long),std_err(i_long))
if(hartree_xc.and.periodicity==3)then
call xyzzyaaav1('MPC estimate of Hartree energy',mpc_hartree*xyzzyaaak&
&16)
call xyzzyaaav1('MPC estimate of XC energy',ave(i_mpc)-mpc_hartree*xyz&
&zyaaak16)
endif
case('mpc_ewald')
call xyzzyaaav1('MPC e-e interaction',ave(i_mpc),std_err(i_mpc))
call xyzzyaaav1('Short range part of MPC',ave(i_short),std_err(i_short&
&))
call xyzzyaaav1('Long range part of MPC',ave(i_long),std_err(i_long))
if(hartree_xc.and.periodicity==3)then
call xyzzyaaav1('MPC estimate of Hartree energy',mpc_hartree*xyzzyaaak&
&16)
call xyzzyaaav1('MPC estimate of XC energy',ave(i_mpc)-mpc_hartree*xyz&
&zyaaak16)
endif
call xyzzyaaav1('Ewald e-e interaction',ave(i_pote),std_err(i_pote))
case('manual')
call xyzzyaaav1('Manual e-e interaction',ave(i_pote),std_err(i_pote))
end select
if(nitot/=0.and..not.ignore_ionic_interactions)then
if(isperiodic)then
if(have_ppots)then
if(use_expot)then
call xyzzyaaav1('Ewald e-i int. (local) + ext. pot.',ave(i_potil),std_&
&err(i_potil))
else
call xyzzyaaav1('Ewald e-i interaction (local)',ave(i_potil),std_err(i&
&_potil))
endif
call xyzzyaaav1('Ewald e-i interaction (non-local)',ave(i_potinl),std_&
&err(i_potinl))
else
if(use_expot)then
call xyzzyaaav1('Ewald e-n int. + external pot.',ave(i_potil),std_err(&
&i_potil))
else
call xyzzyaaav1('Ewald e-n interaction',ave(i_potil),std_err(i_potil))
endif
endif
else
if(have_ppots)then
if(use_expot)then
call xyzzyaaav1('e-i int. (local) + external pot.',ave(i_potil),std_er&
&r(i_potil))
else
call xyzzyaaav1('e-i interaction (local)',ave(i_potil),std_err(i_potil&
&))
endif
call xyzzyaaav1('e-i interaction (non-local)',ave(i_potinl),std_err(i_&
&potinl))
else
if(use_expot)then
call xyzzyaaav1('e-n interaction + external pot.',ave(i_potil),std_err&
&(i_potil))
else
call xyzzyaaav1('e-n interaction',ave(i_potil),std_err(i_potil))
endif
endif
endif
else
if(use_expot)call xyzzyaaav1('External potential',ave(i_potil),std_err&
&(i_potil))
endif
if(have_veep)then
if(isperiodic)then
call xyzzyaaav1('Core polarization energy',ave(i_ecpp_tot),std_err(i_e&
&cpp_tot))
else
call xyzzyaaav1('Total core polarization energy',ave(i_ecpp_tot),std_e&
&rr(i_ecpp_tot))
call xyzzyaaav1('Core polarization (e-i term)',ave(i_vcpp_ei),std_err(&
&i_vcpp_ei))
call xyzzyaaav1('Core polarization (e term)',ave(i_vcpp_e),std_err(i_v&
&cpp_e))
call xyzzyaaav1('Core polarization (e-e term)',ave(i_vcpp_ee),std_err(&
&i_vcpp_ee))
endif
endif
if(relativistic)then
call wout(' Relativistic terms:')
if(eval_maspol)call xyzzyaaav1('Mass polarization',ave(i_emasspol),std&
&_err(i_emasspol),blank=.false.)
if(eval_masvel)call xyzzyaaav1('Mass velocity term',ave(i_emassvel),st&
&d_err(i_emassvel),blank=.false.)
if(eval_darwinen)call xyzzyaaav1('Electron-nucleus Darwin term',ave(i_&
&edarwin_en),std_err(i_edarwin_en),blank=.false.)
if(eval_darwinee)call xyzzyaaav1('Electron-electron Darwin term',ave(i&
&_edarwin_ee),std_err(i_edarwin_ee),blank=.false.)
if(eval_retard)call xyzzyaaav1('Retardation term',ave(i_eretard),std_e&
&rr(i_eretard),blank=.false.)
call xyzzyaaav1('Total relativistic correction',ave(i_reltot),std_err(&
&i_reltot))
endif
if(constant_energy/=0.d0)call xyzzyaaav1('Constant energy contribution&
&s',xyzzyaaae16)
if(isvmc.and.finite_size_corr)then
if(nblock>1)call wout(' Finite size correction data (using block energ&
&ies):')
if(chi_squared_sf>0.d0)call xyzzyaaav1('Chi squared in structure facto&
&r fit',chi_squared_sf,blank=.false.,au=.false.)
call xyzzyaaav1('Total XC correction (dXC)',xyzzyaaai16,blank=.false.)
call xyzzyaaav1('Total KE correction (dKE)',xyzzyaaah16,blank=.false.)
select case(trim(interaction))
case('mpc')
call xyzzyaaav1('Total energy with MPC',ave(i_eloc_mpc),blank=.false.)
call xyzzyaaav1('Total energy with MPC + dKE',ave(i_eloc_mpc)+xyzzyaaa&
&h16,blank=.false.)
case('ewald_mpc')
call xyzzyaaav1('Total energy with Ewald',ave(i_eloc_def),blank=.false&
&.)
call xyzzyaaav1('Total energy with Ewald + dXC+dKE',ave(i_eloc_def)+xy&
&zzyaaah16+xyzzyaaai16,blank=.false.)
call xyzzyaaav1('Total energy with Ewald + dXC',ave(i_eloc_def)+xyzzya&
&aai16,blank=.false.)
call xyzzyaaav1('Total energy with Ewald + dKE',ave(i_eloc_def)+xyzzya&
&aah16,blank=.false.)
call xyzzyaaav1('Total energy with MPC',ave(i_eloc_mpc),blank=.false.)
call xyzzyaaav1('Total energy with MPC + dKE',ave(i_eloc_mpc)+xyzzyaaa&
&h16,blank=.false.)
case('mpc_ewald')
call xyzzyaaav1('Total energy with MPC',ave(i_eloc_mpc),blank=.false.)
call xyzzyaaav1('Total energy with MPC + dKE',ave(i_eloc_mpc)+xyzzyaaa&
&h16,blank=.false.)
call xyzzyaaav1('Total energy with Ewald',ave(i_eloc_def),blank=.false&
&.)
call xyzzyaaav1('Total energy with Ewald + dXC+dKE',ave(i_eloc_def)+xy&
&zzyaaah16+xyzzyaaai16,blank=.false.)
call xyzzyaaav1('Total energy with Ewald + dXC',ave(i_eloc_def)+xyzzya&
&aai16,blank=.false.)
call xyzzyaaav1('Total energy with Ewald + dKE',ave(i_eloc_def)+xyzzya&
&aah16,blank=.false.)
case default
call xyzzyaaav1('Total energy with Ewald',ave(i_eloc_def),blank=.false&
&.)
call xyzzyaaav1('Total energy with Ewald + dXC+dKE',ave(i_eloc_def)+xy&
&zzyaaah16+xyzzyaaai16,blank=.false.)
call xyzzyaaav1('Total energy with Ewald + dXC',ave(i_eloc_def)+xyzzya&
&aai16,blank=.false.)
call xyzzyaaav1('Total energy with Ewald + dKE',ave(i_eloc_def)+xyzzya&
&aah16,blank=.false.)
end select
if(iblock/=nblock)call wout(' The dXC accuracy improves with accumulat&
&ing structure factor.')
call wout()
endif
if(int_sf.and.trim(runtype)=='vmc')then
call eval_int_sf(xyzzyaaaf16,xyzzyaaag16)
if(esupercell)then
xyzzyaaaf16=xyzzyaaaf16*real(npcells,dp)
xyzzyaaag16=xyzzyaaag16*real(npcells,dp)
endif
call wout()
call wout(' Electron-electron interaction energy from structure factor&
&:')
if(model_system)then
call xyzzyaaav1('Hartree energy',xyzzyaaaf16*inv_netot,blank=.false.)
call xyzzyaaav1('Exchange-correlation energy',xyzzyaaag16*inv_netot,bl&
&ank=.false.)
call xyzzyaaav1('Total e-e interaction energy',(xyzzyaaaf16+xyzzyaaag1&
&6)*inv_netot,blank=.false.)
else
call xyzzyaaav1('Hartree energy',xyzzyaaaf16,blank=.false.)
call xyzzyaaav1('Exchange-correlation energy',xyzzyaaag16,blank=.false&
&.)
call xyzzyaaav1('Total e-e interaction energy',xyzzyaaaf16+xyzzyaaag16&
&,blank=.false.)
endif
call wout()
endif
if(nmove>1)then
if(isperiodic)then
call xyzzyaaav1('Var. of local energy per sim. cell',vartemp_tot,varte&
&mp_err)
else
call xyzzyaaav1('Variance of local energy',vartemp_tot,vartemp_err)
endif
endif
if(eval_contact_den)then
call wout(' Contact density between electrons and a positron:')
call xyzzyaaav1(' w',ave(i_contact_den),std_err(i_contact_den),blank=.&
&false.,au=.false.)
call wout()
endif
if(.not.isperiodic)then
if(nitot==0)then
call xyzzyaaav1('Maximum distance from CoM',sqrt(r2maxn))
else
call xyzzyaaav1('Maximum distance from origin',sqrt(r2maxn))
endif
endif
end subroutine xyzzyaaau1
subroutine xyzzyaaav1(label,datum,error,blank,au)
use format_utils, only : wout
implicit none
real(dp),intent(in) :: datum
real(dp),intent(in),optional :: error
logical,intent(in),optional :: blank,au
character(*),intent(in) :: label
logical xyzzyaaaa17,xyzzyaaab17
character(80) wline
xyzzyaaaa17=.true.
if(present(blank))xyzzyaaaa17=blank
xyzzyaaab17=.true.
if(present(au))xyzzyaaab17=au
if(xyzzyaaab17)then
write(wline,'(1x,a,t37,a,t47,f21.12)')label,'(au) =',datum
else
write(wline,'(1x,a,t37,a,t47,f21.12)')label,'     =',datum
endif
call wout(wline)
if(present(error))then
if(error>0.d0)then
write(wline,'(1x,a,t40,a,t47,f21.12)')'Standard error','+/-',error
call wout(wline)
endif
endif
if(xyzzyaaaa17)call wout()
end subroutine xyzzyaaav1
subroutine vmc_main(dt_array_in,dt_shift_array_in,opt_dtvmc_array,nequ&
&il,corper,nmove_in,nblock_in,nwrcon,nvmcave,writeout_vmc_hist,vmc_twi&
&st_av,nequil_ta,vmc_energy,vmc_stderr,aborted,initial_rele,initial_re&
&le_set)
use slaarnaaa
use slaarnaaf
use slaarnaam
use slaarnabg
use parallel
use slaarnacd
use slaarnach
use store
use slaarnacs
use slaarnaab, only : forces_scratch_request,setup_forces_accum,finish&
&_forces_accum,eval_local_forces,tagh_forces,nfterms,nfcomps
use slaarnaad,      only : plot_backflow,backflow_stats
use slaarnaan, only : make_movie,points,vector_difference
use slaarnaaq,        only : expvals,expval_scratch_request,setup_expv&
&al_accum,finish_expval_accum,accumulate_expvals,write_expval,particle&
&_is_fixed,pair_corr,pair_corr_sph,pcf_rfix,pcfs_rfix,fixed_particle,f&
&inite_size_corr_xc
use file_utils,    only : open_units
use format_utils,  only : wout,i2s,r2s,wordwrap
use slaarnaas,     only : mc_twist_offset,eval_finite_hf_energies
use slaarnabk,  only : interaction_mpc_use
use slaarnabt,     only : correlation_time,reaverage,dcopy,correlation&
&_time_alt
use slaarnaca,         only : have_ppots
use slaarnacc,only : ranx,get_random_state,put_random_state
use run_control,   only : errstop,errwarn,timer,exceeds_time_limit,tcp&
&utime,errstop_master,check_alloc
use slaarnacr,     only : check_kinetic
implicit none
integer,intent(in) :: nblock_in,nequil,nmove_in,nwrcon,nvmcave,nequil_&
&ta
real(dp),intent(in) :: initial_rele(3,netot),dt_array_in(no_difftypes)&
&,dt_shift_array_in(no_difftypes)
real(dp),intent(out) :: vmc_energy,vmc_stderr
integer,intent(inout) :: opt_dtvmc_array(no_difftypes),corper
logical,intent(in) :: writeout_vmc_hist,vmc_twist_av,initial_rele_set(&
&netot)
logical,intent(out) :: aborted
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18
real(dp) xyzzyaaae18,xyzzyaaaf18,xyzzyaaag18,xyzzyaaah18
character(80) wline
integer xyzzyaaai18,xyzzyaaaj18,xyzzyaaak18,xyzzyaaal18,xyzzyaaam18,xy&
&zzyaaan18,xyzzyaaao18,xyzzyaaap18,xyzzyaaaq18,xyzzyaaar18,xyzzyaaas18&
&,xyzzyaaat18,xyzzyaaau18,xyzzyaaav18,xyzzyaaaw18,xyzzyaaax18,xyzzyaaa&
&y18,xyzzyaaaz18
real(dp) xyzzyaaba18,xyzzyaabb18,xyzzyaabc18,xyzzyaabd18,xyzzyaabe18,x&
&yzzyaabf18,xyzzyaabg18,xyzzyaabh18(no_difftypes)
integer xyzzyaabi18(no_difftypes),xyzzyaabj18,xyzzyaabk18,xyzzyaabl18,&
&xyzzyaabm18,xyzzyaabn18,xyzzyaabo18
integer,allocatable :: sele(:),xyzzyaabp18(:)
real(dp) xyzzyaabq18(3),xyzzyaabr18,prob,xyzzyaabs18,xyzzyaabt18,xyzzy&
&aabu18(no_difftypes),xyzzyaabv18,xyzzyaabw18,xyzzyaabx18,xyzzyaaby18(&
&3),xyzzyaabz18,xyzzyaaca18(no_difftypes),xyzzyaacb18(no_difftypes),xy&
&zzyaacc18,xyzzyaacd18(no_difftypes),xyzzyaace18(no_difftypes),xyzzyaa&
&cf18,xyzzyaacg18,xyzzyaach18,xyzzyaaci18
real(dp),allocatable :: rele(:,:),xyzzyaacj18(:,:),xyzzyaack18(:,:),xy&
&zzyaacl18(:,:)
logical xyzzyaacm18,xyzzyaacn18,xyzzyaaco18
logical,save :: xyzzyaacp18=.false.,xyzzyaacq18=.false.,xyzzyaacr18=.f&
&alse.,xyzzyaacs18=.false.
integer xyzzyaact18,xyzzyaacu18,xyzzyaacv18,xyzzyaacw18
integer,save :: xyzzyaacx18
real(dp) xyzzyaacy18,etot,xyzzyaacz18,xyzzyaada18,xyzzyaadb18,xyzzyaad&
&c18,xyzzyaadd18,xyzzyaade18,xyzzyaadf18,xyzzyaadg18,xyzzyaadh18,xyzzy&
&aadi18,xyzzyaadj18,xyzzyaadk18,xyzzyaadl18,xyzzyaadm18,xyzzyaadn18,xy&
&zzyaado18,xyzzyaadp18,xyzzyaadq18,xyzzyaadr18,xyzzyaads18,xyzzyaadt18&
&,xyzzyaadu18,xyzzyaadv18,xyzzyaadw18(5)
real(dp),allocatable :: xyzzyaadx18(:),xyzzyaady18(:),xyzzyaadz18(:),x&
&yzzyaaea18(:),xyzzyaaeb18(:),xyzzyaaec18(:,:),xyzzyaaed18(:,:),xyzzya&
&aee18(:),xyzzyaaef18(:,:),xyzzyaaeg18(:),xyzzyaaeh18(:),xyzzyaaei18(:&
&),xyzzyaaej18(:),xyzzyaaek18(:),xyzzyaael18(:),xyzzyaaem18(:),xyzzyaa&
&en18(:),xyzzyaaeo18(:),xyzzyaaep18(:)
logical xyzzyaaeq18,xyzzyaaer18,isnan,isinf
integer xyzzyaaes18,xyzzyaaet18,xyzzyaaeu18,xyzzyaaev18
real(dp) xyzzyaaew18,xyzzyaaex18,xyzzyaaey18,xyzzyaaez18
real(sp) xyzzyaafa18,xyzzyaafb18,xyzzyaafc18
aborted=.false.
call timer('VMC',.true.)
call timer('SETUP',.true.,collapse=.true.)
if(am_master)then
call wout()
call wout('BEGIN VMC CALCULATION')
call wout('=====================')
call wout()
endif
call xyzzyaaff18
if(newrun)then
call points(rele,sele,initial_rele,initial_rele_set,printout=.true.)
call equilibration(xyzzyaact18,xyzzyaacu18,rele,sele,xyzzyaabu18,xyzzy&
&aace18,xyzzyaacd18,opt_dtvmc_array,nequil,corper)
if(corper==0)corper=1
xyzzyaaat18=0
endif
if(any(xyzzyaabu18<=0.d0))call errstop_master('VMC_MAIN','VMC time ste&
&p is zero.')
call define_config(xyzzyaact18,rele,sele)
if(am_master.and.newrun)then
if(use_backflow)call plot_backflow(rele)
call timer('WFN_CHECK',.true.,collapse=.true.)
call check_kinetic(xyzzyaact18,xyzzyaacw18)
call timer('WFN_CHECK',.false.)
endif
if(use_backflow)call backflow_stats(xyzzyaaae18,xyzzyaaaf18)
if(am_master)call wout('Starting VMC.')
call xyzzyaafg18
call timer('SETUP',.false.)
block: do xyzzyaaak18=1,xyzzyaaaw18
if(.not.xyzzyaaak18==1.and.xyzzyaaak18==xyzzyaaai18.and.xyzzyaaay18/=0&
&)then
xyzzyaaaj18=xyzzyaaay18
xyzzyaaax18=xyzzyaaay18
call xyzzyaafj18(.false.,.false.)
endif
if(exceeds_time_limit(xyzzyaaak18==1))then
if(am_master)then
call wout()
call errwarn('VMC','Time limit exceeded or about to be exceeded. Emerg&
&ency stop.')
if(am_master.and.chkpoint_level==-1.and.xyzzyaacp18)then
call wout('Writing config.out file for restart despite CHECKPOINT==-1.&
&')
call wout()
endif
call wout('CONTINUATION INFO:')
if(isopt_vmc.or.isvmc_opt)then
if(opt_cycle>opt_cycles)then
call wout(' Suggested action: roll back and continue')
call wout(' Set RUNTYPE = vmc')
if(old_input)then
call wout(' Set NWRCON = 0')
else
call wout(' Set VMC_NCONFIG_WRITE = 0')
endif
call wout(' Set WRITEOUT_VMC_HIST = F')
call wout(' Move correlation.out.'//trim(i2s(opt_cycle-1))//' to corre&
&lation.data')
else
if(opt_cycle==1)then
call wout(' Suggested action: restart with more time')
call wout(' PROBLEM: first VMC run won''t fit in time slot!')
call wout(' Probably an optimization cycle wouldn''t fit either.')
call wout(' You must re-run using a longer time slot.')
else
call wout(' Suggested action: roll back and continue')
if(isopt_vmc)call wout(' Set RUNTYPE = vmc_opt')
call wout(' Set OPT_CYCLES = '//trim(i2s(opt_cycles-opt_cycle+1)))
call wout(' Move correlation.out.'//trim(i2s(opt_cycle-1))//' to corre&
&lation.data')
endif
endif
else
call wout(' Suggested action: continue run directly')
if(old_input)then
if(use_blocktime)call wout(' Set NMOVE = '//trim(i2s(xyzzyaaaj18)))
call wout(' Set NBLOCK = '//trim(i2s(xyzzyaaai18-xyzzyaaak18+1)))
if(nwrcon>0)call wout(' Set NWRCON = '//trim(i2s(nwrcon-xyzzyaaan18)))
else
if(.not.use_blocktime)then
call wout(' Set VMC_NBLOCK = '//trim(i2s(xyzzyaaai18-xyzzyaaak18+1)))
call wout(' Set VMC_NSTEP = '//trim(i2s((xyzzyaaai18-xyzzyaaak18+1)*xy&
&zzyaaaj18*nnodes)))
else
if(xyzzyaaak18==1)then
call wout(' Set VMC_NSTEP = '//trim(i2s(xyzzyaaaj18*nnodes)))
else
call wout(' Set VMC_NSTEP = '//trim(i2s(((xyzzyaaai18-xyzzyaaak18)*xyz&
&zyaaaj18+xyzzyaaay18)*nnodes)))
endif
endif
if(nwrcon>0)call wout(' Set VMC_NCONFIG_WRITE = '//trim(i2s(nnodes*(nw&
&rcon-xyzzyaaan18))))
endif
call wout(' Set NEWRUN = F')
call wout(' Move config.out to config.in')
call wout()
call wout(' NB: The runqmc options --continue/--auto-continue do this &
&automatically.')
endif
call wout()
endif
if(xyzzyaacp18)then
xyzzyaacq18=.true.
call xyzzyaafk18
endif
aborted=.true.
exit
endif
xyzzyaafa18=tcputime()
if(vmc_twist_av)then
call mc_twist_offset
call eval_finite_hf_energies(xyzzyaadu18,xyzzyaadv18)
if(nnodes>1)then
call mpi_gather(k_offset,3,mpi_double_precision,xyzzyaaef18,3,mpi_doub&
&le_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Gathering k-vector offsets.')
call mpi_reduce(xyzzyaadu18,xyzzyaads18,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing hf_ke in vmc_main.')
call mpi_reduce(xyzzyaadv18,xyzzyaadt18,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing hf_ex in vmc_main.')
xyzzyaads18=xyzzyaads18/dble(nnodes)
xyzzyaadt18=xyzzyaadt18/dble(nnodes)
else
xyzzyaads18=xyzzyaadu18
xyzzyaadt18=xyzzyaadv18
endif
if(am_master)then
call wout()
call wout()
call wout('===========================================================&
&=================')
if(nnodes==1)then
call wout('NEW K-VECTOR OFFSET: ')
call wout('    ',k_offset(1:3))
else
call wout('PROCESSOR                         NEW K-VECTOR OFFSETS')
do xyzzyaaad18=0,nnodes-1
write(wline,'(i6,4x,3(1x,es21.13))')xyzzyaaad18,xyzzyaaef18(1:3,xyzzya&
&aad18)
call wout(wline)
enddo
endif
call wout('===========================================================&
&=================')
call wout()
endif
xyzzyaabi18=0
call equilibration(xyzzyaact18,xyzzyaacu18,rele,sele,xyzzyaabu18,xyzzy&
&aace18,xyzzyaacd18,xyzzyaabi18,nequil_ta,corper)
call define_config(xyzzyaact18,rele,sele)
endif
xyzzyaaal18=0
xyzzyaacl18=0.d0
xyzzyaabs18=0.d0
xyzzyaaec18=0.d0
xyzzyaacf18=0.d0
xyzzyaaca18=0.d0
if(use_altsamp)call complextosimple(xyzzyaact18,xyzzyaacu18)
do xyzzyaaao18=1,xyzzyaaax18
do xyzzyaaap18=1,xyzzyaaaq18
if(xyzzyaaal18==huge(1))call errstop('VMC_MAIN','Integer overflow: ist&
&ep.')
xyzzyaaal18=xyzzyaaal18+1
xyzzyaaeq18=(mod(xyzzyaaap18,corper)==0)
xyzzyaaer18=.false.
if(xyzzyaaeq18.and.xyzzyaaan18<nwrcon)then
xyzzyaaav18=xyzzyaaav18-1
xyzzyaaer18=(xyzzyaaav18==0)
endif
if(.not.vmc_cfg_by_cfg)then
call timer('RANDOM_WALK',.true.)
do xyzzyaabm18=1,netot
if((particle_is_fixed.and.xyzzyaabm18==fixed_particle).or.xyzzyaabm18=&
&=on_top_jj)cycle
if(xyzzyaaad1)then
call wout('PARTICLE #'//trim(i2s(xyzzyaabm18))//':')
call wout(' R_OLD: ',rele(1:dimensionality,xyzzyaabm18))
if(xyzzyaaar18>1)call wout(' S_OLD: '//trim(i2s(sele(xyzzyaabm18))))
endif
xyzzyaabl18=which_difftype(which_spin(xyzzyaabm18))
do xyzzyaaas18=1,xyzzyaaar18
call xyzzyaaak1(xyzzyaabm18,xyzzyaact18,xyzzyaacu18,xyzzyaabu18(xyzzya&
&abl18),xyzzyaace18(xyzzyaabl18),xyzzyaacd18(xyzzyaabl18),xyzzyaaas18,&
&xyzzyaabq18,xyzzyaabj18,xyzzyaacc18)
call vector_difference(rele(1,xyzzyaabm18),xyzzyaabq18,xyzzyaaby18)
xyzzyaabz18=sum(xyzzyaaby18**2)
if(xyzzyaaad1)then
if(xyzzyaaas18==1)then
call wout(' R_NEW: ',xyzzyaabq18(1:dimensionality))
else
call wout(' S_NEW: '//trim(i2s(xyzzyaabj18)))
endif
endif
if(altsamp==1.and.have_ppots)then
if(grid_angles_valid(xyzzyaact18))grid_angles_scr(:,:,:,xyzzyaacu18)=g&
&rid_angles_scr(:,:,:,xyzzyaact18)
grid_angles_valid(xyzzyaacu18)=grid_angles_valid(xyzzyaact18)
endif
do xyzzyaabk18=1,xyzzyaabn18
if(xyzzyaaab1)then
xyzzyaabo18=xyzzyaabk18
else
xyzzyaabo18=0
endif
xyzzyaacn18=any(rele(:,xyzzyaabm18)/=xyzzyaabq18).or.sele(xyzzyaabm18)&
&/=xyzzyaabj18
if(xyzzyaacn18)then
if(.not.use_altsamp)then
call wfn_ratio(xyzzyaact18,xyzzyaacu18,xyzzyaabo18,relprob=xyzzyaabv18&
&,isnan=isnan,isinf=isinf)
else
call wfn_altratio(xyzzyaact18,xyzzyaacu18,xyzzyaabo18,relprob=xyzzyaab&
&v18,isnan=isnan,isinf=isinf)
endif
if(xyzzyaabk18==1)xyzzyaabv18=xyzzyaabv18*xyzzyaacc18
prob=min(1.d0,xyzzyaabv18)
if(xyzzyaaad1)then
if(isnan)then
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaabk18))/&
&/': zero [was Not a Number]')
elseif(isinf)then
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaabk18))/&
&/': zero [it diverged]')
else
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaabk18))/&
&/': ',prob)
endif
endif
if(isnan.or.isinf)prob=0.d0
else
prob=1.d0
if(xyzzyaaad1)call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(&
&xyzzyaabk18))//': unity [no-change move]')
endif
if(xyzzyaaad1)call wout(' ACCEPTANCE PROBABILITY: ',prob)
xyzzyaabr18=0.d0
if(prob<1.d0)xyzzyaabr18=ranx()
xyzzyaacm18=xyzzyaabr18<prob
if(.not.xyzzyaacm18)exit
if(xyzzyaaas18==1)xyzzyaacl18(xyzzyaabk18,xyzzyaabl18)=xyzzyaacl18(xyz&
&zyaabk18,xyzzyaabl18)+1.d0
enddo
if(xyzzyaacm18)then
if(xyzzyaacn18)then
call accept_move(xyzzyaact18,xyzzyaacu18)
if(use_altsamp)call alt_accept(xyzzyaact18,xyzzyaacu18)
rele(:,xyzzyaabm18)=xyzzyaabq18
sele(xyzzyaabm18)=xyzzyaabj18
if(xyzzyaabm18==on_top_ii)then
rele(:,on_top_jj)=xyzzyaabq18
sele(on_top_jj)=xyzzyaabj18
endif
endif
if(xyzzyaaas18==1)then
call xyzzyaafw18(xyzzyaabm18,rele,xyzzyaabs18)
xyzzyaaca18(xyzzyaabl18)=xyzzyaaca18(xyzzyaabl18)+xyzzyaabz18
if(xyzzyaabm18==on_top_ii)xyzzyaaca18(which_difftype(which_spin(on_top&
&_jj)))=xyzzyaaca18(which_difftype(which_spin(on_top_jj)))+xyzzyaabz18
else
xyzzyaacf18=xyzzyaacf18+1.d0
endif
if(xyzzyaaad1)call wout(' MOVE ACCEPTED')
else
if(use_altsamp)call alt_reject(xyzzyaacu18)
if(xyzzyaaad1)call wout(' MOVE REJECTED')
endif
enddo
enddo
call timer('RANDOM_WALK',.false.)
if(xyzzyaaeq18)then
if(use_altsamp)then
call alt_halflogpdf(xyzzyaact18,xyzzyaaey18)
call simpletocomplex(xyzzyaact18,xyzzyaacu18)
endif
call eval_local_energy(xyzzyaact18,etot=etot,ecomps=xyzzyaaee18(1:xyzz&
&yaaet18),isnan=isnan,isinf=isinf)
if(xyzzyaaae1)then
if(isnan)then
call wout('ENERGY: Not a Number')
elseif(isinf)then
call wout('ENERGY: Diverges')
else
call wout('ENERGY: ',etot)
endif
endif
if(isnan.or.isinf)call errstop('VMC_MAIN','Floating-point exception re&
&ported when calculating energy. Don''t know what to do, so stopping. &
&Please file a bug report.')
if(expvals)call accumulate_expvals(xyzzyaact18,1.d0,.false.)
if(forces)call eval_local_forces(xyzzyaact18,etot,xyzzyaaee18(i_kei),x&
&yzzyaaee18(i_potinl),0.d0,xyzzyaaee18(xyzzyaaeu18:xyzzyaaev18))
if(use_altsamp)then
call weight_config_vmc(xyzzyaact18,xyzzyaaew18)
xyzzyaaee18(i_wght)=1.d0
xyzzyaaec18(1:xyzzyaaes18,xyzzyaaao18)=xyzzyaaec18(1:xyzzyaaes18,xyzzy&
&aaao18)+xyzzyaaew18*xyzzyaaee18(1:xyzzyaaes18)
xyzzyaadx18((xyzzyaaal18-1)/corper+1)=xyzzyaaew18*etot
xyzzyaaeg18((xyzzyaaal18-1)/corper+1)=xyzzyaaew18
if(altsamp==1.and.have_ppots)grid_angles_valid(xyzzyaact18)=.false.
else
xyzzyaaec18(1:xyzzyaaes18,xyzzyaaao18)=xyzzyaaec18(1:xyzzyaaes18,xyzzy&
&aaao18)+xyzzyaaee18(1:xyzzyaaes18)
xyzzyaadx18((xyzzyaaal18-1)/corper+1)=etot
endif
endif
else
call timer('RANDOM_WALK',.true.)
if(xyzzyaaad1)then
do xyzzyaabm18=1,netot
call wout('R_'//trim(i2s(xyzzyaabm18))//'_OLD: ',rele(1:dimensionality&
&,xyzzyaabm18))
if(xyzzyaaar18>1)call wout('S_OLD: '//trim(i2s(sele(xyzzyaabm18))))
enddo
endif
do xyzzyaaas18=1,xyzzyaaar18
call xyzzyaaal1(xyzzyaact18,xyzzyaacu18,xyzzyaabu18,xyzzyaace18,xyzzya&
&acd18,xyzzyaaas18,xyzzyaacj18,xyzzyaabp18,xyzzyaacc18)
xyzzyaabz18=0.d0
do xyzzyaabm18=1,netot
call vector_difference(rele(1,xyzzyaabm18),xyzzyaacj18(1,xyzzyaabm18),&
&xyzzyaaby18)
xyzzyaabz18=xyzzyaabz18+sum(xyzzyaaby18**2)
enddo
if(xyzzyaaad1)then
do xyzzyaabm18=1,netot
if(xyzzyaaas18==1)then
call wout('R_'//trim(i2s(xyzzyaabm18))//'_NEW: ',xyzzyaacj18(1:dimensi&
&onality,xyzzyaabm18))
else
call wout('S_'//trim(i2s(xyzzyaabm18))//'_NEW: '//trim(i2s(sele(xyzzya&
&abm18))))
endif
enddo
endif
if(altsamp==1.and.have_ppots)then
if(grid_angles_valid(xyzzyaact18))grid_angles_scr(:,:,:,xyzzyaacu18)=g&
&rid_angles_scr(:,:,:,xyzzyaact18)
grid_angles_valid(xyzzyaacu18)=grid_angles_valid(xyzzyaact18)
endif
xyzzyaabw18=1.d0
xyzzyaacm18=.true.
do xyzzyaabk18=1,xyzzyaabn18
if(xyzzyaaab1)then
xyzzyaabo18=xyzzyaabk18
else
xyzzyaabo18=0
endif
if(.not.use_altsamp)then
call wfn_ratio(xyzzyaact18,xyzzyaacu18,xyzzyaabo18,relprob=xyzzyaabv18&
&,isnan=isnan,isinf=isinf)
else
call wfn_altratio(xyzzyaact18,xyzzyaacu18,xyzzyaabo18,relprob=xyzzyaab&
&v18,isnan=isnan,isinf=isinf)
endif
if(xyzzyaabk18==1)xyzzyaabv18=xyzzyaabv18*xyzzyaacc18
prob=min(1.d0,xyzzyaabv18)
if(xyzzyaaad1)then
if(isnan)then
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaabk18))/&
&/': zero [was Not a Number]')
elseif(isinf)then
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaabk18))/&
&/': zero [it diverged]')
else
call wout(' ACCEPTANCE PROBABILITY AT LEVEL '//trim(i2s(xyzzyaabk18))/&
&/': ',prob)
endif
endif
if(isnan.or.isinf)then
prob=0.d0
xyzzyaabv18=0.d0
endif
xyzzyaabw18=xyzzyaabw18*xyzzyaabv18
if(.not.xyzzyaacm18)cycle
xyzzyaabr18=0.d0
if(prob<1.d0)xyzzyaabr18=ranx()
xyzzyaacm18=xyzzyaabr18<prob
if(.not.xyzzyaacm18)then
if(xyzzyaaeq18.and.xyzzyaabw18>0.d0)cycle
exit
endif
if(xyzzyaaas18==1)xyzzyaacl18(xyzzyaabk18,1)=xyzzyaacl18(xyzzyaabk18,1&
&)+1.d0
enddo
if(xyzzyaaeq18.and.xyzzyaaac1)then
xyzzyaabw18=min(xyzzyaabw18,1.d0)
prob=xyzzyaabw18
xyzzyaabx18=1.d0-xyzzyaabw18
if(xyzzyaaar18==2)then
prob=.5d0*prob
xyzzyaabx18=.5d0*xyzzyaabx18
endif
if(use_altsamp)then
call alt_halflogpdf(xyzzyaact18,xyzzyaaey18)
call alt_halflogpdf(xyzzyaacu18,xyzzyaaez18)
call simpletocomplex(xyzzyaact18,xyzzyaacu18)
endif
call eval_local_energy(xyzzyaact18,etot=xyzzyaadq18,ecomps=xyzzyaaee18&
&(1:xyzzyaaet18),isnan=isnan,isinf=isinf)
if(xyzzyaaae1)then
if(isnan)then
call wout('OLD ENERGY: Not a Number')
elseif(isinf)then
call wout('OLD ENERGY: Diverges')
else
call wout('OLD ENERGY: ',xyzzyaadq18)
endif
endif
if(isnan.or.isinf)then
if(.not.xyzzyaacm18)call errstop('VMC_MAIN','Floating-point exception &
&reported when calculating old energy. Don''t know what to do, so stop&
&ping. Please file a bug report.')
prob=1.d0
xyzzyaabx18=0.d0
endif
if(expvals)call accumulate_expvals(xyzzyaact18,xyzzyaabx18,.false.)
if(forces)call eval_local_forces(xyzzyaact18,xyzzyaadq18,xyzzyaaee18(i&
&_kei),xyzzyaaee18(i_potinl),0.d0,xyzzyaaee18(xyzzyaaeu18:xyzzyaaev18)&
&)
if(use_altsamp)then
call weight_config_vmc(xyzzyaact18,xyzzyaaew18)
xyzzyaaee18(i_wght)=1.d0
xyzzyaaec18(1:xyzzyaaes18,xyzzyaaao18)=xyzzyaaec18(1:xyzzyaaes18,xyzzy&
&aaao18)+xyzzyaaew18*xyzzyaabx18*xyzzyaaee18(1:xyzzyaaes18)
if(altsamp==1.and.have_ppots)grid_angles_valid(xyzzyaact18)=.false.
else
xyzzyaaec18(1:xyzzyaaes18,xyzzyaaao18)=xyzzyaaec18(1:xyzzyaaes18,xyzzy&
&aaao18)+xyzzyaabx18*xyzzyaaee18(1:xyzzyaaes18)
endif
call eval_local_energy(xyzzyaacu18,etot=xyzzyaadr18,ecomps=xyzzyaaee18&
&(1:xyzzyaaet18),isnan=isnan,isinf=isinf)
if(xyzzyaaae1)then
if(isnan)then
call wout('NEW ENERGY: Not a Number')
elseif(isinf)then
call wout('NEW ENERGY: Diverges')
else
call wout('NEW ENERGY: ',xyzzyaadr18)
endif
endif
if(isnan.or.isinf)then
if(xyzzyaacm18)call errstop('VMC_MAIN','Floating-point exception repor&
&ted when calculating new energy. Don''t know what to do, so stopping.&
& Please file a bug report.')
prob=0.d0
xyzzyaabx18=1.d0
endif
if(expvals)call accumulate_expvals(xyzzyaacu18,prob,.false.)
if(forces)call eval_local_forces(xyzzyaacu18,xyzzyaadr18,xyzzyaaee18(i&
&_kei),xyzzyaaee18(i_potinl),0.d0,xyzzyaaee18(xyzzyaaeu18:xyzzyaaev18)&
&)
if(use_altsamp)then
call weight_config_vmc(xyzzyaacu18,xyzzyaaex18)
xyzzyaaee18(i_wght)=1.d0
xyzzyaaec18(1:xyzzyaaes18,xyzzyaaao18)=xyzzyaaec18(1:xyzzyaaes18,xyzzy&
&aaao18)+xyzzyaaex18*prob*xyzzyaaee18(1:xyzzyaaes18)
else
xyzzyaaec18(1:xyzzyaaes18,xyzzyaaao18)=xyzzyaaec18(1:xyzzyaaes18,xyzzy&
&aaao18)+prob*xyzzyaaee18(1:xyzzyaaes18)
endif
if(xyzzyaaas18==1)then
if(use_altsamp)then
xyzzyaadx18((xyzzyaaal18-1)/corper+1)=(1.d0-xyzzyaabw18)*xyzzyaaew18*x&
&yzzyaadq18+xyzzyaabw18*xyzzyaaex18*xyzzyaadr18
xyzzyaaeg18((xyzzyaaal18-1)/corper+1)=(1.d0-xyzzyaabw18)*xyzzyaaew18+x&
&yzzyaabw18*xyzzyaaex18
if(simplepdf==1.and.isitcomplex)then
call empty_scratch_wfn(xyzzyaact18)
call empty_scratch_wfn(xyzzyaacu18)
endif
if(altsamp==1.and.have_ppots)grid_angles_valid(xyzzyaact18)=.false.
else
xyzzyaadx18((xyzzyaaal18-1)/corper+1)=xyzzyaadq18*(1.d0-xyzzyaabw18)+x&
&yzzyaadr18*xyzzyaabw18
endif
endif
endif
if(xyzzyaacm18)then
call accept_move(xyzzyaact18,xyzzyaacu18)
if(use_altsamp)call alt_accept(xyzzyaact18,xyzzyaacu18)
call dcopy(three_netot,xyzzyaacj18(1,1),1,rele(1,1),1)
sele(:)=xyzzyaabp18(:)
if(xyzzyaaas18==1)then
call xyzzyaafw18(0,rele,xyzzyaabs18)
xyzzyaaca18(1)=xyzzyaaca18(1)+xyzzyaabz18
else
xyzzyaacf18=xyzzyaacf18+1.d0
endif
if(xyzzyaaad1)call wout('MOVE ACCEPTED')
else
if(use_altsamp)call alt_reject(xyzzyaacu18)
if(xyzzyaaad1)call wout('MOVE REJECTED')
endif
enddo
call timer('RANDOM_WALK',.false.)
if(xyzzyaaeq18.and..not.xyzzyaaac1)then
if(use_altsamp)then
call alt_halflogpdf(xyzzyaact18,xyzzyaaey18)
call simpletocomplex(xyzzyaact18,xyzzyaacu18)
endif
call eval_local_energy(xyzzyaact18,etot=etot,ecomps=xyzzyaaee18(1:xyzz&
&yaaet18),isnan=isnan,isinf=isinf)
if(xyzzyaaae1)then
if(isnan)then
call wout('ENERGY: Not a Number')
elseif(isinf)then
call wout('ENERGY: Diverges')
else
call wout('ENERGY: ',etot)
endif
endif
if(isnan.or.isinf)call errstop('VMC_MAIN','Floating-point exception re&
&ported when calculating energy. Don''t know what to do, so stopping. &
&Please file a bug report.')
if(expvals)call accumulate_expvals(xyzzyaact18,1.d0,.false.)
if(forces)call eval_local_forces(xyzzyaact18,etot,xyzzyaaee18(i_kei),x&
&yzzyaaee18(i_potinl),0.d0,xyzzyaaee18(xyzzyaaeu18:xyzzyaaev18))
if(use_altsamp)then
call weight_config_vmc(xyzzyaact18,xyzzyaaew18)
xyzzyaaee18(i_wght)=1.d0
xyzzyaaec18(1:xyzzyaaes18,xyzzyaaao18)=xyzzyaaec18(1:xyzzyaaes18,xyzzy&
&aaao18)+xyzzyaaew18*xyzzyaaee18(1:xyzzyaaes18)
xyzzyaadx18((xyzzyaaal18-1)/corper+1)=xyzzyaaew18*etot
xyzzyaaeg18((xyzzyaaal18-1)/corper+1)=xyzzyaaew18
if(altsamp==1.and.have_ppots)grid_angles_valid(xyzzyaact18)=.false.
else
xyzzyaaec18(1:xyzzyaaes18,xyzzyaaao18)=xyzzyaaec18(1:xyzzyaaes18,xyzzy&
&aaao18)+xyzzyaaee18(1:xyzzyaaes18)
xyzzyaadx18((xyzzyaaal18-1)/corper+1)=etot
endif
endif
endif
if(xyzzyaaer18)then
call add_config(rele=rele,sele=sele)
if(use_altsamp)call add_config(modify=.true.,logp=xyzzyaaey18)
call add_config_energy_items(xyzzyaact18)
call add_config_wfn_items(xyzzyaact18)
if(vmc_twist_av)then
xyzzyaadw18(1:3)=k_offset(1:3)
xyzzyaadw18(4:5)=(/xyzzyaadu18,xyzzyaadv18/)
call add_config(modify=.true.,twist=xyzzyaadw18)
endif
xyzzyaaan18=xyzzyaaan18+1
xyzzyaaav18=xyzzyaaau18
endif
if(use_altsamp.and.xyzzyaaeq18)call complextosimple(xyzzyaact18,xyzzya&
&acu18)
if(makemovie.and.my_node==movienode.and.mod(xyzzyaaal18,movieplot)==0)&
&call make_movie(rele)
enddo
if(xyzzyaaak18==1)then
if(use_blocktime.and.(am_master.or.(am_slave.and..not.xyzzyaacs18)))th&
&en
call xyzzyaafi18
endif
if(xyzzyaaao18>=xyzzyaaaj18)then
allocate(xyzzyaaed18(xyzzyaaes18,xyzzyaaaj18),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'VMC_MAIN','ave_nodes')
if(nwrcon>0.and.xyzzyaaao18>xyzzyaaaj18)then
xyzzyaaac18=xyzzyaaan18-xyzzyaaaj18/xyzzyaaau18
if(xyzzyaaac18>0)then
xyzzyaaan18=xyzzyaaan18-xyzzyaaac18
call delete_config_vmc(xyzzyaaac18)
endif
endif
if(xyzzyaaai18==1)xyzzyaaco18=.true.
exit
endif
else
if(xyzzyaaak18==xyzzyaaai18)xyzzyaaco18=.true.
endif
enddo
if(use_altsamp)call simpletocomplex(xyzzyaact18,xyzzyaacu18)
call xyzzyaafm18
xyzzyaaae18=xyzzyaabf18*xyzzyaacy18*xyzzyaabg18
xyzzyaaec18(:,1:xyzzyaaaj18)=xyzzyaaec18(:,1:xyzzyaaaj18)*xyzzyaaae18
if(use_altsamp)then
xyzzyaaec18(i_wght,1:xyzzyaaaj18)=xyzzyaaec18(i_wght,1:xyzzyaaaj18)/xy&
&zzyaacy18
xyzzyaaae18=sum(xyzzyaaec18(i_wght,1:xyzzyaaaj18))
call mpi_reduce(xyzzyaaae18,xyzzyaaaf18,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
call mpi_bcast(xyzzyaaaf18,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
endif
call mpi_reduce(xyzzyaaec18,xyzzyaaed18,xyzzyaaes18*xyzzyaaaj18,mpi_do&
&uble_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Summing sum1 over nodes in VMC.')
if(use_altsamp)then
if(interaction_mpc_use)then
call correlation_time_alt(xyzzyaaec18(i_eloc_mpc,1:xyzzyaaaj18),xyzzya&
&aec18(i_wght,1:xyzzyaaaj18),xyzzyaacz18,xyzzyaada18)
else
call correlation_time_alt(xyzzyaaec18(i_eloc_def,1:xyzzyaaaj18),xyzzya&
&aec18(i_wght,1:xyzzyaaaj18),xyzzyaacz18,xyzzyaada18)
endif
else
if(interaction_mpc_use)then
call correlation_time(xyzzyaaec18(i_eloc_mpc,1:xyzzyaaaj18),xyzzyaacz1&
&8,xyzzyaada18)
else
call correlation_time(xyzzyaaec18(i_eloc_def,1:xyzzyaaaj18),xyzzyaacz1&
&8,xyzzyaada18)
endif
endif
if(am_master)then
allocate(xyzzyaaeh18(nnodes),xyzzyaaei18(nnodes),stat=xyzzyaaaa18)
else
allocate(xyzzyaaeh18(1),xyzzyaaei18(1),stat=xyzzyaaaa18)
endif
call mpi_gather(xyzzyaacz18,1,mpi_double_precision,xyzzyaaeh18,1,mpi_d&
&ouble_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering corr_tau in VMC')
call mpi_gather(xyzzyaada18,1,mpi_double_precision,xyzzyaaei18,1,mpi_d&
&ouble_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering corr_tau_err in VMC')
if(am_master)then
if(use_altsamp)then
xyzzyaaae18=xyzzyaaba18/sum(xyzzyaaed18(i_wght,1:xyzzyaaaj18))
do xyzzyaaac18=1,xyzzyaaes18
xyzzyaadz18(xyzzyaaac18)=sum(xyzzyaaed18(xyzzyaaac18,1:xyzzyaaaj18))*x&
&yzzyaabd18*xyzzyaaae18
xyzzyaaea18(xyzzyaaac18)=sum((xyzzyaaed18(xyzzyaaac18,1:xyzzyaaaj18)-x&
&yzzyaaed18(i_wght,1:xyzzyaaaj18)*xyzzyaadz18(xyzzyaaac18))**2)*xyzzya&
&abd18
xyzzyaaaf18=xyzzyaaea18(xyzzyaaac18)*xyzzyaabe18
xyzzyaaeb18(xyzzyaaac18)=sqrt(max(0.d0,xyzzyaaaf18))*xyzzyaaae18
enddo
do xyzzyaaao18=1,xyzzyaaaj18
do xyzzyaaac18=1,xyzzyaaes18
if(xyzzyaaac18==i_wght)cycle
xyzzyaaed18(xyzzyaaac18,xyzzyaaao18)=xyzzyaaed18(xyzzyaaac18,xyzzyaaao&
&18)*xyzzyaaae18
enddo
enddo
else
do xyzzyaaac18=1,xyzzyaaes18
xyzzyaadz18(xyzzyaaac18)=sum(xyzzyaaed18(xyzzyaaac18,1:xyzzyaaaj18))*x&
&yzzyaabd18
xyzzyaaea18(xyzzyaaac18)=sum((xyzzyaaed18(xyzzyaaac18,1:xyzzyaaaj18)-x&
&yzzyaadz18(xyzzyaaac18))**2)*xyzzyaabd18
xyzzyaaaf18=xyzzyaaea18(xyzzyaaac18)*xyzzyaabe18
xyzzyaaeb18(xyzzyaaac18)=sqrt(max(0.d0,xyzzyaaaf18))
enddo
endif
if(any(xyzzyaaei18<0.d0))then
xyzzyaacz18=-1.d0
xyzzyaada18=-1.d0
else
call reaverage(nnodes,xyzzyaaaj18,xyzzyaaeh18,xyzzyaaei18,1.d0,xyzzyaa&
&cz18,xyzzyaada18)
endif
if(writeout_vmc_hist.and..not.particle_is_fixed)call xyzzyaafn18
xyzzyaaat18=xyzzyaaat18+xyzzyaaaj18
if(use_altsamp)then
xyzzyaaae18=xyzzyaaba18/sum(xyzzyaaed18(i_wght,1:xyzzyaaaj18))
do xyzzyaaac18=1,xyzzyaaaj18
call reblock_add(xyzzyaaed18(:,xyzzyaaac18)-xyzzyaaed18(i_wght,xyzzyaa&
&ac18)*xyzzyaadz18(:)*xyzzyaaae18,1.d0)
enddo
else
do xyzzyaaac18=1,xyzzyaaaj18
call reblock_add(xyzzyaaed18(:,xyzzyaaac18),1.d0)
enddo
endif
endif
deallocate(xyzzyaaeh18,xyzzyaaei18)
if(use_altsamp)then
call xyzzyaafp18
else
call xyzzyaafo18
endif
if(expvals)call write_expval
if(finite_size_corr.and.am_master.and.isvmc)call finite_size_corr_xc
if(am_master)then
xyzzyaafb18=tcputime()-xyzzyaafa18
call xyzzyaaau1(xyzzyaaak18,xyzzyaaaj18,xyzzyaaai18,xyzzyaaes18,xyzzya&
&adz18,xyzzyaaeb18,xyzzyaack18,xyzzyaacz18,xyzzyaada18,xyzzyaadb18,xyz&
&zyaadc18,xyzzyaacb18,xyzzyaabt18,xyzzyaacg18,xyzzyaafb18,xyzzyaaaj18*&
&nnodes)
endif
call xyzzyaafu18
if(use_backflow.and.simplepdf==0)then
call backflow_stats(xyzzyaaae18,xyzzyaaaf18)
xyzzyaaag18=xyzzyaaaf18
if(vmc_cfg_by_cfg)xyzzyaaag18=xyzzyaaae18
call mpi_reduce(xyzzyaaag18,xyzzyaaah18,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
if(am_master)then
xyzzyaaah18=xyzzyaaah18*xyzzyaabf18
tmpr=r2s(xyzzyaaah18,'(f8.4)')
if(.not.vmc_cfg_by_cfg)then
call wout(' Particles affected per move       (%) :  '//trim(tmpr))
else
call wout(' Particles within backflow range   (%) :  '//trim(tmpr))
endif
call wout()
endif
endif
call xyzzyaafk18
if(am_master)then
call wout()
call wout('Time taken in block    : : : ',dble(xyzzyaafb18),rfmt='(f13&
&.4)')
call wout()
endif
xyzzyaacp18=.true.
if(xyzzyaaco18)exit block
enddo block
if(am_master)call wout(repeat('=',73))
if(am_master)call xyzzyaafv18
call end_config_accumulation(xyzzyaacq18)
call xyzzyaafr18
call timer('VMC',.false.)
contains
subroutine xyzzyaafd18
implicit none
call open_units(xyzzyaacx18,xyzzyaaab18)
if(xyzzyaaab18/=0)call errstop('VMC','Unable to find free i/o unit (hi&
&st).')
end subroutine xyzzyaafd18
subroutine xyzzyaafe18
implicit none
open_unit(xyzzyaacx18)=.false.
end subroutine xyzzyaafe18
subroutine xyzzyaaff18
implicit none
integer xyzzyaaaa21,xyzzyaaab21,xyzzyaaac21
logical xyzzyaaad21(4)
logical,allocatable :: xyzzyaaae21(:)
character(20) config_dum(0),extra_dum(0),extra_item(4)
character(20),dimension(:),allocatable :: xyzzyaaaf21
call xyzzyaafd18
allocate(rele(3,netot),sele(netot),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','basic arrays')
rele=0.d0
sele=0
if(vmc_cfg_by_cfg)then
allocate(xyzzyaacj18(3,netot),xyzzyaabp18(netot),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','cfg-by-cfg arrays')
xyzzyaacj18=0.d0
xyzzyaabp18=0
endif
xyzzyaabu18=dt_array_in
xyzzyaace18=sqrt(xyzzyaabu18)
xyzzyaacd18=dt_shift_array_in
xyzzyaaaa21=7
if(noncoll_spin)xyzzyaaaa21=xyzzyaaaa21+1
if(use_altsamp)xyzzyaaaa21=xyzzyaaaa21+1
if(vmc_twist_av)xyzzyaaaa21=xyzzyaaaa21+1
allocate(xyzzyaaaf21(xyzzyaaaa21),xyzzyaaae21(xyzzyaaaa21),stat=xyzzya&
&aaa18)
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','config_item')
xyzzyaaab21=7
xyzzyaaaf21(1:7)=(/'RELE           ','LOGDET         ','ETOT          &
& ','FIDET          ',  'LAPDET         ','LOCAL_POTENTIAL','NLTOT    &
&      '/)
if(noncoll_spin)then
xyzzyaaab21=xyzzyaaab21+1
xyzzyaaaf21(xyzzyaaab21)='SELE'
endif
if(use_altsamp)then
xyzzyaaab21=xyzzyaaab21+1
xyzzyaaaf21(xyzzyaaab21)='LOGP'
endif
if(vmc_twist_av)then
xyzzyaaab21=xyzzyaaab21+1
xyzzyaaaf21(xyzzyaaab21)='TWIST'
endif
extra_item(1:4)=(/'VMC_SAVED_STATE ','RANDOM          ','FINAL_VMC_RES&
&ULT','REBLOCK_DATA    '/)
xyzzyaaac21=0
if(nwrcon>0)xyzzyaaac21=-1
call load_configs(xyzzyaaac21,'VMC_OPT',config_dum,xyzzyaaaf21,xyzzyaa&
&ae21,extra_dum,extra_item,xyzzyaaad21)
deallocate(xyzzyaaae21)
if(xyzzyaaad21(2))call put_random_state(random_state_config)
call xyzzyaaft18(xyzzyaaad21(3))
if(newrun)then
if(xyzzyaaac21>0.or.any(xyzzyaaad21))then
if(xyzzyaaad21(1))call errwarn('INITIALIZE_VMC','Discarding existing c&
&heckpoint data as NEWRUN=T.')
call dismantle_configs
endif
else
if(.not.xyzzyaaad21(1))call errstop_master('INITIALIZE_VMC','NEWRUN is&
& .false. but no previous VMC state has been loaded. Check that the co&
&nfig file corresponds to a previous VMC run.')
if(am_master)then
call wout('Starting from saved VMC state.')
call wout()
endif
call xyzzyaafl18
endif
xyzzyaaac21=0
if(nwrcon>0)xyzzyaaac21=nwrcon
call init_config_accumulation('VMC',xyzzyaaac21,xyzzyaaaf21,extra_item&
&)
deallocate(xyzzyaaaf21)
if(use_blocktime)then
stop_method='nstep'
if(trim(adjustl(stop_method))=='nstep')then
xyzzyaaaw18=huge(1)
xyzzyaaai18=1
xyzzyaaaj18=nblock_in*nmove_in
xyzzyaaax18=xyzzyaaaj18
else
xyzzyaaaw18=huge(1)
xyzzyaaax18=huge(1)
xyzzyaaai18=1
if(netot<20)then
xyzzyaaaj18=1000000
elseif(netot<100)then
xyzzyaaaj18=100000
elseif(netot<1000)then
xyzzyaaaj18=50000
else
xyzzyaaaj18=10000
endif
xyzzyaafc18=0.
endif
else
xyzzyaaaj18=nmove_in
xyzzyaaai18=nblock_in
xyzzyaaaw18=nblock_in
xyzzyaaax18=nmove_in
endif
xyzzyaaco18=.false.
xyzzyaaay18=xyzzyaaaj18
allocate(xyzzyaadx18(xyzzyaaaj18*nvmcave),stat=xyzzyaaaa18)
if(use_altsamp)allocate(xyzzyaaeg18(xyzzyaaaj18*nvmcave),stat=xyzzyaaa&
&a18)
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','variance arrays (1).')
if(am_master)then
if(use_blocktime)then
allocate(xyzzyaady18(1),stat=xyzzyaaaa18)
else
allocate(xyzzyaady18(nblock_in),stat=xyzzyaaaa18)
endif
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','variance arrays (2).')
endif
xyzzyaact18=0
xyzzyaacu18=0
xyzzyaacv18=0
xyzzyaacw18=0
call scratch_protect(xyzzyaact18)
if(vmc_cfg_by_cfg)call scratch_protect(xyzzyaacu18)
call scratch_request(kinetic=xyzzyaact18)
call scratch_request(ratio1_from=xyzzyaact18,ratio1_to=xyzzyaacw18)
if(vmc_cfg_by_cfg)then
call scratch_request(ratiocfg_from=xyzzyaact18,ratiocfg_to=xyzzyaacu18&
&)
else
call scratch_request(ratio1_from=xyzzyaact18,ratio1_to=xyzzyaacu18)
if(on_top_ii/=0)call scratch_request(ratio2_from=xyzzyaact18,ratio2_to&
&=xyzzyaacu18)
endif
call energy_scratch_request(xyzzyaact18)
if(use_altsamp)then
call scratch_request(wfn_detail=xyzzyaact18)
call scratch_request(wfn_detail=xyzzyaacu18)
call scratch_request(drift=xyzzyaact18)
call scratch_request(drift=xyzzyaacu18)
call scratch_request(kinetic=xyzzyaacu18)
call energy_scratch_request(xyzzyaacu18)
endif
if(altsamp==1.and.have_ppots)call scratch_request(ratio1_from=xyzzyaac&
&u18,ratio1_to=xyzzyaacv18)
if(expvals)call expval_scratch_request(xyzzyaact18)
if(forces)call forces_scratch_request(xyzzyaact18)
if(vmc_cfg_by_cfg.and.xyzzyaaac1)then
call energy_scratch_request(xyzzyaacu18)
if(expvals)call expval_scratch_request(xyzzyaacu18)
if(forces)call forces_scratch_request(xyzzyaacu18)
endif
call setup_scratch
call which_scratch(xyzzyaact18)
call which_scratch(xyzzyaacu18)
if(use_altsamp.and.altsamp==1.and.have_ppots)call which_scratch(xyzzya&
&acv18)
call which_scratch(xyzzyaacw18)
call setup_wfn_utils
call setup_energy_utils
if(use_altsamp)call setup_alt_utils
if(expvals)call setup_expval_accum
if(forces)call setup_forces_accum
call xyzzyaafh18
call timer('SETUP',.false.)
if(use_altsamp)scr_nl_2=xyzzyaacv18
if(xyzzyaaab1)then
xyzzyaabn18=nsampling_levels
else
xyzzyaabn18=1
endif
if(.not.vmc_cfg_by_cfg)then
allocate(xyzzyaacl18(xyzzyaabn18,no_difftypes),xyzzyaack18(xyzzyaabn18&
&,no_difftypes),stat=xyzzyaaaa18)
else
allocate(xyzzyaacl18(xyzzyaabn18,1),xyzzyaack18(xyzzyaabn18,1),stat=xy&
&zzyaaaa18)
endif
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','iaccept_lev')
if(am_master)then
if(newrun)then
call init_reblock(xyzzyaaes18)
else
if(xyzzyaaad21(1).and.xyzzyaaad21(4))then
call restore_reblock_data
else
call errwarn('INITIALIZE_VMC','VMC reblock data not in config.in file &
&- probably produced with version of CASINO earlier than 2.13.333. On-&
&the-fly reblocked error bars through restart will be incorrect - use &
&correlation time estimate instead.')
call init_reblock(xyzzyaaes18)
endif
endif
endif
if(vmc_twist_av.and.nnodes>1)then
if(am_master)then
allocate(xyzzyaaef18(3,0:nnodes-1),stat=xyzzyaaaa18)
else
allocate(xyzzyaaef18(1,1),stat=xyzzyaaaa18)
endif
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','k_offset_all')
endif
end subroutine xyzzyaaff18
subroutine xyzzyaafg18
use slaarnaag, only : max_rep_int
implicit none
integer xyzzyaaaa22
xyzzyaack18=0.d0
xyzzyaacb18=0.d0
xyzzyaabb18=real(nnodes,dp)
xyzzyaaba18=real(xyzzyaaaj18,dp)
xyzzyaabc18=real(nvmcave,dp)
xyzzyaaan18=0
if(nwrcon>0)then
if((isvmc.or.isopt_vmc.or.isvmc_opt).and.trim(adjustl(stop_method))=='&
&target_error')then
xyzzyaaau18=1
else
if(real(xyzzyaaaj18,dp)*real(xyzzyaaai18,dp)*real(nvmcave,dp)>max_rep_&
&int)call errstop ('INPUT_SETUP','Integer overflow : nblock*nmove*nvmc&
&ave')
xyzzyaaau18=(xyzzyaaaj18*xyzzyaaai18*nvmcave)/nwrcon
endif
if(xyzzyaaau18==0)call errstop('INPUT_SETUP','Not enough moves to writ&
&e all configs.')
if(isvmc_opt.or.isopt_vmc.or.isvmc_dmc)then
if(xyzzyaaau18>corper)then
corper=min(corper_default_vmc,corper)
else
corper=corper/xyzzyaaau18
if(mod(corper_default_opt,xyzzyaaau18)>0)corper=corper+1
endif
endif
else
if(isvmc_opt.or.isopt_vmc)corper=corper_default_vmc
endif
xyzzyaaav18=xyzzyaaau18
xyzzyaaaq18=nvmcave*corper
xyzzyaabf18=1.d0/xyzzyaabb18
xyzzyaabd18=1.d0/xyzzyaaba18
xyzzyaabg18=1.d0/xyzzyaabc18
if(xyzzyaaaj18>1)then
xyzzyaabe18=1.d0/(xyzzyaaba18-1.d0)
else
xyzzyaabe18=1.d0
endif
if(model_system.and.isperiodic)then
xyzzyaacy18=inv_netot
else
xyzzyaacy18=1.d0/real(npcells,dp)
endif
xyzzyaaam18=xyzzyaaaj18*corper*nvmcave
xyzzyaach18=1.d0/(real(xyzzyaaam18,dp)*xyzzyaabb18)
xyzzyaaci18=xyzzyaach18/real(netot,dp)
xyzzyaabh18=0.d0
do xyzzyaaaa22=1,nspin
xyzzyaabh18(which_difftype(xyzzyaaaa22))=xyzzyaabh18(which_difftype(xy&
&zzyaaaa22))+real(nele(xyzzyaaaa22),dp)
enddo
where(xyzzyaabh18>0.d0)xyzzyaabh18=xyzzyaach18/xyzzyaabh18
xyzzyaaar18=1
if(noncoll_spin)xyzzyaaar18=2
end subroutine xyzzyaafg18
subroutine xyzzyaafh18
implicit none
xyzzyaaes18=n_ecomp
xyzzyaaet18=xyzzyaaes18
if(forces)then
xyzzyaaeu18=xyzzyaaes18+1
xyzzyaaev18=xyzzyaaes18+nfcomps
xyzzyaaes18=xyzzyaaev18
endif
allocate(xyzzyaaec18(xyzzyaaes18,xyzzyaaaj18),xyzzyaaee18(xyzzyaaes18)&
&,xyzzyaadz18(xyzzyaaes18),xyzzyaaea18(xyzzyaaes18),xyzzyaaeb18(xyzzya&
&aes18),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'SETUP_ACCUMULATION_ARRAYS','basic')
end subroutine xyzzyaafh18
subroutine xyzzyaafi18
use slaarnaag, only : max_rep_int
implicit none
xyzzyaafc18=tcputime()-xyzzyaafa18
if(xyzzyaafc18>block_time)then
if(am_master)xyzzyaaaj18=xyzzyaaao18
call mpi_bcast(xyzzyaaaj18,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting istop in check_blocktime')
xyzzyaaaz18=xyzzyaaaj18
xyzzyaacs18=.true.
if(trim(adjustl(stop_method))=='nstep')then
xyzzyaaai18=xyzzyaaax18/xyzzyaaaj18
xyzzyaaay18=mod(xyzzyaaax18,xyzzyaaaj18)
if(xyzzyaaay18>0)xyzzyaaai18=xyzzyaaai18+1
else
endif
if(xyzzyaaaj18<xyzzyaaax18)call xyzzyaafj18(.true.,.true.)
xyzzyaaax18=xyzzyaaaj18
else
if(xyzzyaaao18==xyzzyaaaj18)then
if(trim(adjustl(stop_method))/='nstep')then
xyzzyaaaz18=xyzzyaaaj18
if(am_master)then
xyzzyaaba18=1.1d0*(xyzzyaaba18/(real(xyzzyaafc18/block_time,dp)))
if(xyzzyaaba18>max_rep_int)then
call errstop('CHECK_BLOCKTIME','Encountered integer overflow when incr&
&easing nmove.')
else
xyzzyaaaj18=int(xyzzyaaba18)
endif
endif
call mpi_bcast(xyzzyaaaj18,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nmove in check_blocktime <2>')
call xyzzyaafj18(.true.,.false.)
endif
endif
endif
end subroutine xyzzyaafi18
subroutine xyzzyaafj18(resize_nmove_arrays,know_nblock)
implicit none
logical,intent(in) :: resize_nmove_arrays,know_nblock
real(dp),allocatable :: xyzzyaaaa25(:,:),xyzzyaaab25(:)
xyzzyaaba18=real(xyzzyaaaj18,dp)
xyzzyaabd18=1.d0/xyzzyaaba18
if(xyzzyaaaj18>1)xyzzyaabe18=1.d0/(xyzzyaaba18-1.d0)
if(xyzzyaaaj18==0)xyzzyaabe18=1.d0
xyzzyaaam18=xyzzyaaaj18*corper*nvmcave
xyzzyaach18=1.d0/(real(xyzzyaaam18,dp)*xyzzyaabb18)
xyzzyaaci18=xyzzyaach18/real(netot,dp)
xyzzyaabh18=0.d0
do xyzzyaaac18=1,nspin
xyzzyaabh18(which_difftype(xyzzyaaac18))=xyzzyaabh18(which_difftype(xy&
&zzyaaac18))+real(nele(xyzzyaaac18),dp)
enddo
where(xyzzyaabh18>0.d0)xyzzyaabh18=xyzzyaach18/xyzzyaabh18
if(resize_nmove_arrays)then
allocate(xyzzyaaaa25(xyzzyaaes18,xyzzyaaaz18),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','blocktemp array')
xyzzyaaaa25(:,1:xyzzyaaaz18)=xyzzyaaec18(:,1:xyzzyaaaz18)
deallocate(xyzzyaaec18)
allocate(xyzzyaaec18(xyzzyaaes18,xyzzyaaaj18),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','sum1 array')
xyzzyaaec18(:,1:xyzzyaaaz18)=xyzzyaaaa25(:,1:xyzzyaaaz18)
deallocate(xyzzyaaaa25)
allocate(xyzzyaaab25(xyzzyaaaz18*nvmcave),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','blocktemp2 array')
xyzzyaaab25(1:xyzzyaaaz18*nvmcave)=xyzzyaadx18(1:xyzzyaaaz18*nvmcave)
deallocate(xyzzyaadx18)
allocate(xyzzyaadx18(xyzzyaaaj18*nvmcave),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','energy_array')
xyzzyaadx18(1:xyzzyaaaz18*nvmcave)=xyzzyaaab25(1:xyzzyaaaz18*nvmcave)
if(use_altsamp)then
xyzzyaaab25(1:xyzzyaaaz18*nvmcave)=xyzzyaaeg18(1:xyzzyaaaz18*nvmcave)
deallocate(xyzzyaaeg18)
allocate(xyzzyaaeg18(xyzzyaaaj18*nvmcave),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','weight_array')
xyzzyaaeg18(1:xyzzyaaaz18*nvmcave)=xyzzyaaab25(1:xyzzyaaaz18*nvmcave)
endif
deallocate(xyzzyaaab25)
endif
if(know_nblock)then
xyzzyaaae18=xyzzyaaeo18(1)
deallocate(xyzzyaaeo18)
allocate(xyzzyaaeo18(xyzzyaaai18))
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','nmove_block')
xyzzyaaeo18(1)=xyzzyaaae18
xyzzyaaaf18=xyzzyaaep18(1)
deallocate(xyzzyaaep18)
allocate(xyzzyaaep18(xyzzyaaai18))
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','nnm1_block')
xyzzyaaep18(1)=xyzzyaaaf18
xyzzyaaag18=xyzzyaaej18(1)
deallocate(xyzzyaaej18)
allocate(xyzzyaaej18(xyzzyaaai18))
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','E_block')
xyzzyaaej18(1)=xyzzyaaag18
xyzzyaaah18=xyzzyaaek18(1)
deallocate(xyzzyaaek18)
allocate(xyzzyaaek18(xyzzyaaai18))
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','dE_block')
xyzzyaaek18(1)=xyzzyaaah18
xyzzyaaae18=xyzzyaael18(1)
deallocate(xyzzyaael18)
allocate(xyzzyaael18(xyzzyaaai18))
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','tau_block')
xyzzyaael18(1)=xyzzyaaae18
if(xyzzyaacr18)then
xyzzyaaaf18=xyzzyaaem18(1)
deallocate(xyzzyaaem18)
allocate(xyzzyaaem18(xyzzyaaai18))
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','E_block2')
xyzzyaaem18(1)=xyzzyaaaf18
xyzzyaaag18=xyzzyaaen18(1)
deallocate(xyzzyaaen18)
allocate(xyzzyaaen18(xyzzyaaai18))
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','dE_block2')
xyzzyaaen18(1)=xyzzyaaag18
endif
if(am_master)then
xyzzyaaah18=xyzzyaady18(1)
deallocate(xyzzyaady18)
allocate(xyzzyaady18(xyzzyaaai18))
call check_alloc(xyzzyaaaa18,'RESET_NMOVE_DEPS','variance_array')
xyzzyaady18(1)=xyzzyaaah18
endif
endif
end subroutine xyzzyaafj18
subroutine xyzzyaafk18
implicit none
logical,save :: xyzzyaaaa26=.false.
rele_vmc_config=rele
if(noncoll_spin)sele_vmc_config=sele
dtvmc_array_config=xyzzyaabu18
vmc_steps_config=xyzzyaaat18
call get_random_state(random_state_config)
if((xyzzyaaak18==xyzzyaaai18.or.xyzzyaacq18).and.rng_restart_safe)xyzz&
&yaaaa26=.true.
if(am_master)call save_reblock_data
select case(chkpoint_level)
case(-1)
if(xyzzyaacq18)call write_configs(xyzzyaaaa26)
case(0)
if(xyzzyaacq18)then
call write_configs(xyzzyaaaa26)
else
if(xyzzyaaak18==xyzzyaaai18.and.(isvmc.or.nwrcon==0))call write_config&
&s(xyzzyaaaa26)
endif
case(1)
if(xyzzyaaak18==xyzzyaaai18.or.isvmc)call write_configs(xyzzyaaaa26)
case(2)
call write_configs(xyzzyaaaa26)
end select
end subroutine xyzzyaafk18
subroutine xyzzyaafl18
implicit none
rele=rele_vmc_config
if(all(rele==0.d0))then
call errwarn('VMC_LOAD_STATE','Config file does not contain a saved st&
&ate for node '//trim(i2s(my_node))//'. Using an unequilibrated config&
&uration.')
call points(rele,sele,initial_rele,initial_rele_set)
else
if(noncoll_spin)sele=sele_vmc_config
where(opt_dtvmc_array>0)
xyzzyaabu18=dtvmc_array_config
xyzzyaace18=sqrt(dtvmc_array_config)
opt_dtvmc_array=-1
endwhere
endif
xyzzyaaat18=vmc_steps_config
if(makemovie.and.my_node==movienode)call make_movie(rele)
if(particle_is_fixed)call xyzzyaafs18
end subroutine xyzzyaafl18
subroutine xyzzyaafm18
implicit none
integer xyzzyaaaa28
real(dp) xyzzyaaab28(no_difftypes),xyzzyaaac28(xyzzyaabn18,no_difftype&
&s)
xyzzyaaac28(:,:)=0.d0
if(.not.vmc_cfg_by_cfg)then
call mpi_reduce(xyzzyaacl18,xyzzyaaac28,xyzzyaabn18*no_difftypes,mpi_d&
&ouble_precision,mpi_sum,0,mpi_comm_world,ierror)
do xyzzyaaaa28=1,xyzzyaabn18
xyzzyaack18(xyzzyaaaa28,:)=xyzzyaaac28(xyzzyaaaa28,:)*xyzzyaabh18(:)
enddo
call mpi_reduce(xyzzyaaca18,xyzzyaaab28,no_difftypes,mpi_double_precis&
&ion,mpi_sum,0,mpi_comm_world,ierror)
xyzzyaacb18(:)=xyzzyaaca18(:)*xyzzyaabh18(:)
else
call mpi_reduce(xyzzyaacl18,xyzzyaaac28,xyzzyaabn18,mpi_double_precisi&
&on,mpi_sum,0,mpi_comm_world,ierror)
xyzzyaack18(:,1)=xyzzyaaac28(:,1)*xyzzyaach18
call mpi_reduce(xyzzyaaca18,xyzzyaaab28,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
xyzzyaacb18(1)=xyzzyaaca18(1)*xyzzyaaci18
endif
if(.not.isperiodic)then
call mpi_reduce(xyzzyaabs18,xyzzyaabt18,1,mpi_double_precision,mpi_max&
&,0,mpi_comm_world,ierror)
else
xyzzyaabt18=0.d0
endif
if(noncoll_spin)then
call mpi_reduce(xyzzyaacf18,xyzzyaacg18,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
if(vmc_cfg_by_cfg)then
xyzzyaacg18=xyzzyaacg18*xyzzyaach18
else
xyzzyaacg18=xyzzyaacg18*xyzzyaaci18
endif
else
xyzzyaacg18=0.d0
endif
end subroutine xyzzyaafm18
subroutine xyzzyaafn18
use slaarnabj
implicit none
integer xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaae29,xy&
&zzyaaaf29
character(30) anumber
character(640) writeme
open(xyzzyaacx18,file='vmc.hist',status='unknown',position='append',io&
&stat=xyzzyaaab18)
if(xyzzyaaab18/=0)call errstop('VMC','Error opening vmc.hist file.')
do xyzzyaaao18=1,xyzzyaaaj18
writeme=trim(i2s(xyzzyaaat18+xyzzyaaao18))
do xyzzyaaaa29=2,no_cols_qmc
if(xyzzyaaaa29==tagh_energy)then
if(interaction_mpc_use)then
write(anumber,*)xyzzyaaed18(i_eloc_mpc,xyzzyaaao18)
else
write(anumber,*)xyzzyaaed18(i_eloc_def,xyzzyaaao18)
endif
elseif(xyzzyaaaa29==tagh_etotalt)then
if(interaction_mpc_use)then
write(anumber,*)xyzzyaaed18(i_eloc_def,xyzzyaaao18)
else
write(anumber,*)xyzzyaaed18(i_eloc_mpc,xyzzyaaao18)
endif
elseif(xyzzyaaaa29==tagh_esqr)then
write(anumber,*)xyzzyaaed18(i_esqr,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_k)then
write(anumber,*)xyzzyaaed18(i_kei,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_t)then
write(anumber,*)xyzzyaaed18(i_ti,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_fisq)then
write(anumber,*)xyzzyaaed18(i_fisq,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_ewald)then
write(anumber,*)xyzzyaaed18(i_pote,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_local)then
write(anumber,*)xyzzyaaed18(i_potil,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_nonlocal)then
write(anumber,*)xyzzyaaed18(i_potinl,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_short)then
write(anumber,*)xyzzyaaed18(i_short,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_long)then
write(anumber,*)xyzzyaaed18(i_long,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_cppei)then
if(isperiodic)then
write(anumber,*)xyzzyaaed18(i_ecpp_tot,xyzzyaaao18)
else
write(anumber,*)xyzzyaaed18(i_vcpp_ei,xyzzyaaao18)
endif
elseif(xyzzyaaaa29==tagh_cppe)then
write(anumber,*)xyzzyaaed18(i_vcpp_e,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_cppee)then
write(anumber,*)xyzzyaaed18(i_vcpp_ee,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_masspol)then
write(anumber,*)xyzzyaaed18(i_emasspol,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_massvel)then
write(anumber,*)xyzzyaaed18(i_emassvel,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_darwinen)then
write(anumber,*)xyzzyaaed18(i_edarwin_en,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_darwinee)then
write(anumber,*)xyzzyaaed18(i_edarwin_ee,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_retard)then
write(anumber,*)xyzzyaaed18(i_eretard,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_dipole1)then
write(anumber,*)xyzzyaaed18(i_dipole1,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_dipole2)then
write(anumber,*)xyzzyaaed18(i_dipole2,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_dipole3)then
write(anumber,*)xyzzyaaed18(i_dipole3,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_dipole_sq)then
write(anumber,*)xyzzyaaed18(i_dipole_sq,xyzzyaaao18)
elseif(forces)then
xyzzyaaaf29=xyzzyaaeu18-1
do xyzzyaaad29=1,nitot_forces
do xyzzyaaae29=1,naxis_forces
do xyzzyaaab29=1,nfterms
xyzzyaaaf29=xyzzyaaaf29+1
if(xyzzyaaaa29==tagh_forces(xyzzyaaab29,xyzzyaaae29,xyzzyaaad29))write&
&(anumber,*)xyzzyaaed18(xyzzyaaaf29,xyzzyaaao18)
enddo
enddo
enddo
elseif(xyzzyaaaa29==tagh_contact_den)then
write(anumber,*)xyzzyaaed18(i_contact_den,xyzzyaaao18)
elseif(xyzzyaaaa29==tagh_hf_ke)then
write(anumber,*)xyzzyaads18
elseif(xyzzyaaaa29==tagh_hf_ex)then
write(anumber,*)xyzzyaadt18
else
call errstop('WRITE_VMC_HIST','Bug.')
endif
writeme=trim(writeme)//' '//trim(anumber)
xyzzyaaac29=modulo(xyzzyaaaa29,25)
if(xyzzyaaac29==0)then
write(xyzzyaacx18,'(a)')trim(writeme)
writeme=''
endif
enddo
if(.not.(xyzzyaaac29==0))write(xyzzyaacx18,'(a)')trim(writeme)
enddo
close(xyzzyaacx18)
end subroutine xyzzyaafn18
subroutine xyzzyaafo18
implicit none
real(dp) xyzzyaaaa30,xyzzyaaab30,xyzzyaaac30,xyzzyaaad30,xyzzyaaae30,x&
&yzzyaaaf30
integer xyzzyaaag30
xyzzyaaac30=0.d0
do xyzzyaaag30=1,nvmcave
xyzzyaaaa30=sum(xyzzyaadx18(xyzzyaaag30:nvmcave*xyzzyaaaj18:nvmcave))
xyzzyaaab30=sum(xyzzyaadx18(xyzzyaaag30:nvmcave*xyzzyaaaj18:nvmcave)**&
&2)
xyzzyaaac30=xyzzyaaac30+(xyzzyaaab30-xyzzyaaaa30**2*xyzzyaabd18)*xyzzy&
&aabe18
enddo
call mpi_reduce(xyzzyaaac30,xyzzyaadb18,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
if(am_master)then
xyzzyaadb18=xyzzyaadb18*xyzzyaabg18*xyzzyaabf18
xyzzyaady18(xyzzyaaak18)=max(xyzzyaadb18,0.d0)
endif
if(am_master)then
if(interaction_mpc_use)then
xyzzyaaaf30=xyzzyaadz18(i_eloc_mpc)/xyzzyaacy18
else
xyzzyaaaf30=xyzzyaadz18(i_eloc_def)/xyzzyaacy18
endif
else
xyzzyaaaf30=0.d0
endif
call mpi_bcast(xyzzyaaaf30,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting etemp in vmc')
xyzzyaaaa30=0.d0
xyzzyaaab30=0.d0
do xyzzyaaal18=1,xyzzyaaaj18
xyzzyaaad30=sum((xyzzyaadx18((xyzzyaaal18-1)*nvmcave+1:xyzzyaaal18*nvm&
&cave)-xyzzyaaaf30)**2)*xyzzyaabg18
xyzzyaaaa30=xyzzyaaaa30+xyzzyaaad30
xyzzyaaab30=xyzzyaaab30+xyzzyaaad30**2
enddo
xyzzyaaae30=(xyzzyaaab30*xyzzyaabd18-(xyzzyaaaa30*xyzzyaabd18)**2)*xyz&
&zyaabe18
call mpi_reduce(xyzzyaaae30,xyzzyaadc18,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
if(am_master)xyzzyaadc18=sqrt(max(xyzzyaadc18,0.d0))*xyzzyaabf18
end subroutine xyzzyaafo18
subroutine xyzzyaafp18
implicit none
integer xyzzyaaaa31
real(dp) xyzzyaaab31,xyzzyaaac31,xyzzyaaad31,xyzzyaaae31,xyzzyaaaf31,x&
&yzzyaaag31,xyzzyaaah31,xyzzyaaai31,xyzzyaaaj31,xyzzyaaak31,xyzzyaaal3&
&1,xyzzyaaam31,xyzzyaaan31
xyzzyaaae31=0.d0
xyzzyaaaf31=0.d0
do xyzzyaaaa31=1,nvmcave
xyzzyaaab31=sum(xyzzyaadx18(xyzzyaaaa31::nvmcave))/sum(xyzzyaaeg18(xyz&
&zyaaaa31::nvmcave))
xyzzyaaac31=sum((xyzzyaadx18(xyzzyaaaa31::nvmcave)-xyzzyaaeg18(xyzzyaa&
&aa31::nvmcave)*xyzzyaaab31)**2)
xyzzyaaad31=sum(xyzzyaaeg18(xyzzyaaaa31::nvmcave))
xyzzyaaae31=xyzzyaaae31+xyzzyaaac31
xyzzyaaaf31=xyzzyaaaf31+xyzzyaaad31
enddo
xyzzyaaae31=xyzzyaaae31*xyzzyaabe18
xyzzyaaaf31=xyzzyaaaf31*xyzzyaabd18
call mpi_reduce(xyzzyaaae31,xyzzyaaag31,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
call mpi_reduce(xyzzyaaaf31,xyzzyaaah31,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
if(am_master)then
xyzzyaaag31=xyzzyaaag31*xyzzyaabg18*xyzzyaabf18
xyzzyaaah31=xyzzyaaah31*xyzzyaabg18*xyzzyaabf18
xyzzyaadb18=xyzzyaaag31/xyzzyaaah31**2
xyzzyaady18(xyzzyaaak18)=max(xyzzyaadb18,0.d0)
endif
if(am_master)then
if(interaction_mpc_use)then
xyzzyaaam31=xyzzyaadz18(i_eloc_mpc)/xyzzyaacy18
else
xyzzyaaam31=xyzzyaadz18(i_eloc_def)/xyzzyaacy18
endif
else
xyzzyaaam31=0.d0
endif
call mpi_bcast(xyzzyaaam31,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting etemp in vmc')
xyzzyaaac31=sum(xyzzyaaeg18(1:xyzzyaaaj18*nvmcave))
call mpi_reduce(xyzzyaaac31,xyzzyaaan31,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
xyzzyaaan31=xyzzyaaan31*xyzzyaabd18*xyzzyaabg18*xyzzyaabf18
call mpi_bcast(xyzzyaaan31,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting wtemp in vmc')
xyzzyaaai31=0.d0
xyzzyaaaj31=0.d0
do xyzzyaaal18=1,xyzzyaaaj18
xyzzyaaak31=sum((xyzzyaadx18((xyzzyaaal18-1)*nvmcave+1:xyzzyaaal18*nvm&
&cave)-xyzzyaaeg18((xyzzyaaal18-1)*nvmcave+1:xyzzyaaal18*nvmcave)*xyzz&
&yaaam31)**2)
xyzzyaaak31=xyzzyaaak31/xyzzyaaan31**2
xyzzyaaak31=xyzzyaaak31*xyzzyaabg18
xyzzyaaai31=xyzzyaaai31+xyzzyaaak31
xyzzyaaaj31=xyzzyaaaj31+xyzzyaaak31**2
enddo
xyzzyaaal31=(xyzzyaaaj31*xyzzyaabd18-(xyzzyaaai31*xyzzyaabd18)**2)*xyz&
&zyaabe18
call mpi_reduce(xyzzyaaal31,xyzzyaadc18,1,mpi_double_precision,mpi_sum&
&,0,mpi_comm_world,ierror)
if(am_master)xyzzyaadc18=sqrt(max(xyzzyaadc18,0.d0))*xyzzyaabf18
end subroutine xyzzyaafp18
subroutine xyzzyaafq18
implicit none
real(dp) xyzzyaaaa32,xyzzyaaab32,xyzzyaaac32
if(am_master)then
xyzzyaaac32=sum(xyzzyaaeo18(1:xyzzyaaai18))
xyzzyaaaa32=max(dot_product(xyzzyaaeo18(1:xyzzyaaai18),xyzzyaady18(1:x&
&yzzyaaai18))/xyzzyaaac32,0.d0)
xyzzyaaab32=0.d0
do xyzzyaaak18=1,xyzzyaaai18
xyzzyaaab32=xyzzyaaab32+(xyzzyaady18(xyzzyaaak18)-xyzzyaaaa32)**2
enddo
if(xyzzyaaai18>1)then
xyzzyaaab32=sqrt(xyzzyaaab32/real(xyzzyaaai18*(xyzzyaaai18-1),dp))
else
xyzzyaaab32=-1.d0
endif
if(xyzzyaaaj18>1)then
call wout()
tmpr=r2s(xyzzyaaaa32,'(f21.12)')
if(xyzzyaaab32>0.d0)then
tmpr2=r2s(xyzzyaaab32,'(f21.12)')
call wout(' Sample variance of E_L (au^2/sim.cell) : '//trim(tmpr)//' &
&+- '//trim(tmpr2))
else
call wout(' Sample variance of E_L (au^2/sim.cell) : '//trim(tmpr))
endif
call wout()
call wout(repeat('=',73))
call wout()
endif
endif
end subroutine xyzzyaafq18
subroutine xyzzyaafr18
implicit none
deallocate(rele,sele)
if(vmc_cfg_by_cfg)deallocate(xyzzyaacj18,xyzzyaabp18)
if(am_master)deallocate(xyzzyaady18)
deallocate(xyzzyaaec18,xyzzyaaee18,xyzzyaaed18,xyzzyaadz18,xyzzyaaea18&
&,xyzzyaaeb18,xyzzyaadx18)
if(use_altsamp)deallocate(xyzzyaaeg18)
deallocate(xyzzyaacl18)
if(forces)call finish_forces_accum
if(expvals)call finish_expval_accum
call xyzzyaafe18
if(am_master)call finish_reblock
call finish_energy_utils
call finish_wfn_utils
if(use_altsamp)call finish_alt_utils
call finish_scratch
if(vmc_twist_av.and.nnodes>1)deallocate(xyzzyaaef18)
end subroutine xyzzyaafr18
subroutine xyzzyaafs18
implicit none
logical xyzzyaaaa34
xyzzyaaaa34=.false.
if(pair_corr)then
xyzzyaaaa34=any(rele(:,fixed_particle)/=pcf_rfix(:))
elseif(pair_corr_sph)then
xyzzyaaaa34=any(rele(:,fixed_particle)/=pcfs_rfix(:))
else
call errstop('CHECK_PCF_RFIX','Called in inappropriate circumstances.'&
&)
endif
if(xyzzyaaaa34)then
call wordwrap('Coordinates of fixed particle read from config.in are n&
&ot what was specified in input/expval.data.')
call errstop('CHECK_PCF_RFIX','Quitting.')
endif
end subroutine xyzzyaafs18
subroutine xyzzyaaft18(have_data)
implicit none
logical,intent(in) :: have_data
if(trim(interaction)=='ewald_mpc'.or.trim(interaction)=='mpc_ewald')xy&
&zzyaacr18=.true.
if(use_blocktime)then
allocate(xyzzyaaeo18(1),xyzzyaaep18(1),xyzzyaaej18(1),xyzzyaaek18(1),x&
&yzzyaael18(1),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','E_block..')
if(xyzzyaacr18)then
allocate(xyzzyaaem18(1),xyzzyaaen18(1),stat=xyzzyaaaa18)
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','E_block2..')
endif
else
allocate(xyzzyaaeo18(nblock_in),xyzzyaaep18(nblock_in),xyzzyaaej18(nbl&
&ock_in),xyzzyaaek18(nblock_in),xyzzyaael18(nblock_in),stat=xyzzyaaaa1&
&8)
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','E_block.. <2>')
if(xyzzyaacr18)then
allocate(xyzzyaaem18(nblock_in),xyzzyaaen18(nblock_in),stat=xyzzyaaaa1&
&8)
call check_alloc(xyzzyaaaa18,'INITIALIZE_VMC','E_block2.. <2>')
endif
endif
if(have_data)then
xyzzyaadd18=final_vmce_config
xyzzyaade18=final_vmcde_config
xyzzyaadf18=final_vmcdeu_config
xyzzyaadk18=final_vmcvar_config
xyzzyaadl18=final_vmcmove_config
xyzzyaadg18=final_vmctau_config
if(xyzzyaacr18)then
xyzzyaadh18=final_vmce2_config
xyzzyaadi18=final_vmcde2_config
xyzzyaadj18=final_vmcde2u_config
xyzzyaadi18=final_vmcde2_config
xyzzyaadj18=final_vmcde2u_config
endif
else
xyzzyaadl18=0.d0
endif
end subroutine xyzzyaaft18
subroutine xyzzyaafu18
implicit none
real(dp) xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36,xyzzyaaad36,xyzzyaaae36,x&
&yzzyaaaf36,xyzzyaaag36,xyzzyaaah36,xyzzyaaai36,xyzzyaaaj36,xyzzyaaak3&
&6,xyzzyaaal36,xyzzyaaam36,xyzzyaaan36,xyzzyaaao36,xyzzyaaap36,xyzzyaa&
&aq36,xyzzyaaar36,xyzzyaaas36
if(am_master)then
if(esupercell)then
xyzzyaaas36=1.d0
else
xyzzyaaas36=1.d0/xyzzyaacy18
endif
if(xyzzyaacr18)then
if(interaction_mpc_use)then
xyzzyaaej18(xyzzyaaak18)=xyzzyaadz18(i_eloc_mpc)*xyzzyaaas36
xyzzyaaek18(xyzzyaaak18)=xyzzyaaeb18(i_eloc_mpc)*xyzzyaaas36
xyzzyaaem18(xyzzyaaak18)=xyzzyaadz18(i_eloc_def)*xyzzyaaas36
xyzzyaaen18(xyzzyaaak18)=xyzzyaaeb18(i_eloc_def)*xyzzyaaas36
else
xyzzyaaej18(xyzzyaaak18)=xyzzyaadz18(i_eloc_def)*xyzzyaaas36
xyzzyaaek18(xyzzyaaak18)=xyzzyaaeb18(i_eloc_def)*xyzzyaaas36
xyzzyaaem18(xyzzyaaak18)=xyzzyaadz18(i_eloc_mpc)*xyzzyaaas36
xyzzyaaen18(xyzzyaaak18)=xyzzyaaeb18(i_eloc_mpc)*xyzzyaaas36
endif
else
if(interaction_mpc_use)then
xyzzyaaej18(xyzzyaaak18)=xyzzyaadz18(i_eloc_mpc)*xyzzyaaas36
xyzzyaaek18(xyzzyaaak18)=xyzzyaaeb18(i_eloc_mpc)*xyzzyaaas36
else
xyzzyaaej18(xyzzyaaak18)=xyzzyaadz18(i_eloc_def)*xyzzyaaas36
xyzzyaaek18(xyzzyaaak18)=xyzzyaaeb18(i_eloc_def)*xyzzyaaas36
endif
endif
xyzzyaael18(xyzzyaaak18)=xyzzyaacz18
xyzzyaaeo18(xyzzyaaak18)=real(xyzzyaaaj18,dp)*xyzzyaabb18
xyzzyaaep18(xyzzyaaak18)=xyzzyaaeo18(xyzzyaaak18)*(xyzzyaaeo18(xyzzyaa&
&ak18)-1.d0)
xyzzyaaan36=sum(xyzzyaaeo18(1:xyzzyaaak18))
xyzzyaaaa36=1.d0/xyzzyaaan36
if(xyzzyaaak18==1)then
xyzzyaaae36=xyzzyaael18(1)
if(xyzzyaaae36<0.d0)xyzzyaaae36=1.d0
vmc_energy=xyzzyaaej18(1)
xyzzyaadm18=xyzzyaaek18(1)
vmc_stderr=sqrt(xyzzyaaae36)*xyzzyaadm18
else
xyzzyaaab36=dot_product(xyzzyaaeo18(1:xyzzyaaak18),xyzzyaaej18(1:xyzzy&
&aaak18))*xyzzyaaaa36
xyzzyaaac36=dot_product(xyzzyaaeo18(1:xyzzyaaak18),xyzzyaaej18(1:xyzzy&
&aaak18)**2)*xyzzyaaaa36
xyzzyaaad36=dot_product(xyzzyaaep18(1:xyzzyaaak18),xyzzyaaek18(1:xyzzy&
&aaak18)**2)*xyzzyaaaa36
xyzzyaaae36=dot_product(xyzzyaaeo18(1:xyzzyaaak18),xyzzyaael18(1:xyzzy&
&aaak18))*xyzzyaaaa36
if(xyzzyaaae36<0.d0)xyzzyaaae36=1.d0
vmc_energy=xyzzyaaab36
xyzzyaadm18=(xyzzyaaac36-xyzzyaaab36*xyzzyaaab36+xyzzyaaad36)/(xyzzyaa&
&an36-1.d0)
xyzzyaadm18=sqrt(max(xyzzyaadm18,0.d0))
vmc_stderr=sqrt(xyzzyaaae36)*xyzzyaadm18
endif
xyzzyaaao36=dot_product(xyzzyaaeo18(1:xyzzyaaak18),xyzzyaady18(1:xyzzy&
&aaak18))*xyzzyaaaa36
if(xyzzyaacr18)then
if(xyzzyaaak18>1)then
xyzzyaaab36=dot_product(xyzzyaaeo18(1:xyzzyaaak18),xyzzyaaem18(1:xyzzy&
&aaak18))*xyzzyaaaa36
xyzzyaaac36=dot_product(xyzzyaaeo18(1:xyzzyaaak18),xyzzyaaem18(1:xyzzy&
&aaak18)**2)*xyzzyaaaa36
xyzzyaaad36=dot_product(xyzzyaaep18(1:xyzzyaaak18),xyzzyaaen18(1:xyzzy&
&aaak18)**2)*xyzzyaaaa36
xyzzyaadn18=xyzzyaaab36
xyzzyaado18=(xyzzyaaac36-xyzzyaaab36*xyzzyaaab36+xyzzyaaad36)/(xyzzyaa&
&an36-1.d0)
xyzzyaado18=sqrt(max(xyzzyaado18,0.d0))
xyzzyaadp18=sqrt(xyzzyaaae36)*xyzzyaado18
else
xyzzyaadn18=xyzzyaaem18(1)
xyzzyaado18=xyzzyaaen18(1)
xyzzyaadp18=sqrt(xyzzyaaae36)*xyzzyaado18
endif
endif
if(xyzzyaadl18>0.d0)then
xyzzyaaaf36=xyzzyaaan36
xyzzyaaag36=xyzzyaadl18
xyzzyaaah36=vmc_energy
xyzzyaaai36=xyzzyaadd18
xyzzyaaaj36=xyzzyaadm18**2
xyzzyaaak36=xyzzyaadf18**2
xyzzyaaal36=xyzzyaaao36
xyzzyaaam36=xyzzyaadk18
xyzzyaaan36=xyzzyaaaf36+xyzzyaaag36
xyzzyaaab36=(xyzzyaaah36*xyzzyaaaf36+xyzzyaaai36*xyzzyaaag36)/xyzzyaaa&
&n36
xyzzyaaae36=(xyzzyaaae36*xyzzyaaaf36+xyzzyaadg18*xyzzyaaag36)/xyzzyaaa&
&n36
xyzzyaaap36=xyzzyaaah36*xyzzyaaah36*xyzzyaaaf36+xyzzyaaai36*xyzzyaaai3&
&6*xyzzyaaag36-xyzzyaaab36*xyzzyaaab36*xyzzyaaan36
xyzzyaaaq36=xyzzyaaaf36*(xyzzyaaaf36-1.d0)*xyzzyaaaj36+xyzzyaaag36*(xy&
&zzyaaag36-1.d0)*xyzzyaaak36
xyzzyaaar36=(xyzzyaaaf36-1.d0)*xyzzyaaal36+(xyzzyaaag36-1.d0)*xyzzyaaa&
&m36
vmc_energy=xyzzyaaab36
xyzzyaadm18=sqrt((xyzzyaaaq36+xyzzyaaap36)/(xyzzyaaan36*(xyzzyaaan36-1&
&.d0)))
xyzzyaaao36=(xyzzyaaar36+xyzzyaaap36)/(xyzzyaaan36-1.d0)
vmc_stderr=sqrt(xyzzyaaae36)*xyzzyaadm18
if(xyzzyaacr18)then
xyzzyaaah36=xyzzyaadn18
xyzzyaaai36=xyzzyaadh18
xyzzyaaaj36=xyzzyaado18**2
xyzzyaaak36=xyzzyaadj18**2
xyzzyaaab36=(xyzzyaaah36*xyzzyaaaf36+xyzzyaaai36*xyzzyaaag36)/xyzzyaaa&
&n36
xyzzyaaap36=xyzzyaaah36*xyzzyaaah36*xyzzyaaaf36+xyzzyaaai36*xyzzyaaai3&
&6*xyzzyaaag36-xyzzyaaab36*xyzzyaaab36*xyzzyaaan36
xyzzyaaaq36=xyzzyaaaf36*(xyzzyaaaf36-1.d0)*xyzzyaaaj36+xyzzyaaag36*(xy&
&zzyaaag36-1.d0)*xyzzyaaak36
xyzzyaadn18=xyzzyaaab36
xyzzyaado18=sqrt((xyzzyaaaq36+xyzzyaaap36)/(xyzzyaaan36*(xyzzyaaan36-1&
&.d0)))
xyzzyaadp18=sqrt(xyzzyaaae36)*xyzzyaado18
endif
endif
endif
call mpi_bcast(vmc_energy,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting vmc_energy in VMC.')
call mpi_bcast(vmc_stderr,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting vmc_stderr in VMC.')
call mpi_bcast(xyzzyaadm18,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting vmc_stderru in VMC.')
call mpi_bcast(xyzzyaaao36,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting vmc_variance in VMC.')
call mpi_bcast(xyzzyaaae36,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting tc in VMC.')
call mpi_bcast(xyzzyaaan36,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting totmove in VMC.')
if(xyzzyaacr18)then
call mpi_bcast(xyzzyaadn18,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting vmc_energy2 in VMC.')
call mpi_bcast(xyzzyaadp18,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting vmc_stderr2 in VMC.')
call mpi_bcast(xyzzyaado18,1,mpi_double_precision,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting vmc_stderr2u in VMC.')
endif
final_vmce_config=vmc_energy
final_vmcde_config=vmc_stderr
final_vmcdeu_config=xyzzyaadm18
final_vmcvar_config=xyzzyaaao36
final_vmctau_config=xyzzyaaae36
final_vmcmove_config=xyzzyaaan36
if(xyzzyaacr18)then
final_vmce2_config=xyzzyaadn18
final_vmcde2_config=xyzzyaadp18
final_vmcde2u_config=xyzzyaado18
endif
end subroutine xyzzyaafu18
subroutine xyzzyaafv18
implicit none
integer xyzzyaaaa37,xyzzyaaab37,xyzzyaaac37
real(dp) xyzzyaaad37(xyzzyaaes18),xyzzyaaae37(xyzzyaaes18),xyzzyaaaf37&
&(xyzzyaaes18)
if(.not.esupercell)then
vmc_energy=vmc_energy*xyzzyaacy18
xyzzyaadm18=xyzzyaadm18*xyzzyaacy18
vmc_stderr=vmc_stderr*xyzzyaacy18
if(xyzzyaacr18)then
xyzzyaadn18=xyzzyaadn18*xyzzyaacy18
xyzzyaado18=xyzzyaado18*xyzzyaacy18
xyzzyaadp18=xyzzyaadp18*xyzzyaacy18
endif
endif
if(abs(vmc_energy)<10.d0)then
xyzzyaaac37=1
elseif(abs(vmc_energy)<100.d0)then
xyzzyaaac37=2
elseif(abs(vmc_energy)<1000.d0)then
xyzzyaaac37=3
else
xyzzyaaac37=4
endif
call wout('FINAL RESULT:')
call wout()
if(xyzzyaadl18>0.d0)then
call wout('Restarted calculation. Include data from previous runs.')
call wout()
endif
call wout(' VMC energy (au)    Standard error      Correction for seri&
&al correlation')
call wout()
if(xyzzyaacr18)then
if(interaction_mpc_use)then
call wout(' MPC interaction')
else
call wout(' Ewald interaction')
endif
endif
tmpr=r2s(vmc_energy,'(f21.12)')
tmpr2=r2s(xyzzyaadm18,'(f21.12)')
call wout(trim(tmpr)//' +/- '//trim(tmpr2)//'      No correction')
tmpr3=r2s(vmc_stderr,'(f21.12)')
if(trim(tmpr2)==trim(tmpr3))then
select case(xyzzyaaac37)
case (1)
call wout(' Insufficient data                      Correlation time me&
&thod')
case (2)
call wout(' Insufficient data                       Correlation time m&
&ethod')
case (3)
call wout(' Insufficient data                        Correlation time &
&method')
case default
call wout(' Insufficient data                         Correlation time&
& method')
end select
else
call wout(trim(tmpr)//' +/- '//trim(tmpr3)//'      Correlation time me&
&thod')
endif
call get_reblocked_ave(xyzzyaaad37,xyzzyaaae37,xyzzyaaaf37)
if(esupercell.and..not.model_system.and.npcells>1)then
xyzzyaaad37(:)=xyzzyaaad37(:)*real(npcells,dp)
xyzzyaaae37(:)=xyzzyaaae37(:)*real(npcells,dp)
xyzzyaaaf37(:)=xyzzyaaaf37(:)*real(npcells,dp)
endif
if(.not.use_altsamp)then
if(xyzzyaacr18)then
if(interaction_mpc_use)then
xyzzyaaaa37=i_eloc_mpc
tmpr=r2s(xyzzyaaad37(xyzzyaaaa37),'(f21.12)')
tmpr2=r2s(xyzzyaaae37(xyzzyaaaa37),'(f21.12)')
xyzzyaaab37=i_eloc_def
tmpr3=r2s(xyzzyaaad37(xyzzyaaab37),'(f21.12)')
tmpr4=r2s(xyzzyaaae37(xyzzyaaab37),'(f21.12)')
else
xyzzyaaaa37=i_eloc_def
tmpr=r2s(xyzzyaaad37(xyzzyaaaa37),'(f21.12)')
tmpr2=r2s(xyzzyaaae37(xyzzyaaaa37),'(f21.12)')
xyzzyaaab37=i_eloc_mpc
tmpr3=r2s(xyzzyaaad37(xyzzyaaab37),'(f21.12)')
tmpr4=r2s(xyzzyaaae37(xyzzyaaab37),'(f21.12)')
endif
else
if(interaction_mpc_use)then
xyzzyaaaa37=i_eloc_mpc
tmpr=r2s(xyzzyaaad37(xyzzyaaaa37),'(f21.12)')
tmpr2=r2s(xyzzyaaae37(xyzzyaaaa37),'(f21.12)')
else
xyzzyaaaa37=i_eloc_def
tmpr=r2s(xyzzyaaad37(xyzzyaaaa37),'(f21.12)')
tmpr2=r2s(xyzzyaaae37(xyzzyaaaa37),'(f21.12)')
endif
endif
else
if(xyzzyaacr18)then
if(interaction_mpc_use)then
xyzzyaaaa37=i_eloc_mpc
tmpr=r2s(vmc_energy,'(f21.12)')
tmpr2=r2s(xyzzyaaae37(xyzzyaaaa37),'(f21.12)')
xyzzyaaab37=i_eloc_def
tmpr3=r2s(xyzzyaadn18,'(f21.12)')
tmpr4=r2s(xyzzyaaae37(xyzzyaaab37),'(f21.12)')
else
xyzzyaaaa37=i_eloc_def
tmpr=r2s(vmc_energy,'(f21.12)')
tmpr2=r2s(xyzzyaaae37(xyzzyaaaa37),'(f21.12)')
xyzzyaaab37=i_eloc_mpc
tmpr3=r2s(xyzzyaadn18,'(f21.12)')
tmpr4=r2s(xyzzyaaae37(xyzzyaaab37),'(f21.12)')
endif
else
if(interaction_mpc_use)then
xyzzyaaaa37=i_eloc_mpc
tmpr=r2s(vmc_energy,'(f21.12)')
tmpr2=r2s(xyzzyaaae37(xyzzyaaaa37),'(f21.12)')
else
xyzzyaaaa37=i_eloc_def
tmpr=r2s(vmc_energy,'(f21.12)')
tmpr2=r2s(xyzzyaaae37(xyzzyaaaa37),'(f21.12)')
endif
endif
endif
if(xyzzyaaaf37(xyzzyaaaa37)==0.d0)then
select case(xyzzyaaac37)
case (1)
call wout(' Insufficient data                      On-the-fly reblocki&
&ng method')
case (2)
call wout(' Insufficient data                       On-the-fly reblock&
&ing method')
case (3)
call wout(' Insufficient data                        On-the-fly rebloc&
&king method')
case default
call wout(' Insufficient data                         On-the-fly reblo&
&cking method')
end select
else
call wout(trim(tmpr)//' +/- '//trim(tmpr2)//'      On-the-fly reblocki&
&ng method')
endif
if(xyzzyaacr18)then
call wout()
if(interaction_mpc_use)then
call wout(' Ewald interaction')
else
call wout(' MPC interaction')
endif
tmpr=r2s(xyzzyaadn18,'(f21.12)')
tmpr2=r2s(xyzzyaado18,'(f21.12)')
call wout(trim(tmpr)//' +/- '//trim(tmpr2)//'      No correction')
tmpr5=r2s(xyzzyaadp18,'(f21.12)')
if(trim(tmpr2)==trim(tmpr5))then
select case(xyzzyaaac37)
case (1)
call wout(' Insufficient data                      Correlation time me&
&thod')
case (2)
call wout(' Insufficient data                       Correlation time m&
&ethod')
case (3)
call wout(' Insufficient data                        Correlation time &
&method')
case default
call wout(' Insufficient data                         Correlation time&
& method')
end select
else
call wout(trim(tmpr)//' +/- '//trim(tmpr5)//'      Correlation time me&
&thod')
endif
if(xyzzyaaaf37(xyzzyaaab37)==0.d0)then
select case(xyzzyaaac37)
case (1)
call wout(' Insufficient data                      On-the-fly reblocki&
&ng method')
case (2)
call wout(' Insufficient data                       On-the-fly reblock&
&ing method')
case (3)
call wout(' Insufficient data                        On-the-fly rebloc&
&king method')
case default
call wout(' Insufficient data                         On-the-fly reblo&
&cking method')
end select
else
call wout(trim(tmpr3)//' +/- '//trim(tmpr4)//'      On-the-fly reblock&
&ing method')
endif
endif
if(xyzzyaaaf37(xyzzyaaaa37)>0.1d0*xyzzyaaae37(xyzzyaaaa37))then
call wout()
call wordwrap('Bad reblock convergence - probably not enough data samp&
&les.')
call wout()
if(xyzzyaacr18)then
if(interaction_mpc_use)then
call reblock_dump(xyzzyaaaa37,"energy (MPC)")
else
call reblock_dump(xyzzyaaaa37,"energy (Ewald)")
endif
else
call reblock_dump(xyzzyaaaa37,"energy")
endif
else
if(xyzzyaacr18)then
if(xyzzyaaaf37(xyzzyaaab37)>0.1d0*xyzzyaaae37(xyzzyaaab37))then
call wout()
call wordwrap('Bad reblock convergence - probably not enough data samp&
&les.')
call wout()
if(interaction_mpc_use)then
call reblock_dump(xyzzyaaab37,"energy (Ewald)")
else
call reblock_dump(xyzzyaaab37,"energy (MPC)")
endif
endif
endif
endif
call xyzzyaafq18
end subroutine xyzzyaafv18
subroutine xyzzyaafw18(xyzzyaabm18,rele,xyzzyaabs18)
implicit none
integer,intent(in) :: xyzzyaabm18
real(dp),intent(in) :: rele(3,netot)
real(dp),intent(inout) :: xyzzyaabs18
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38
real(dp) xyzzyaaad38(3),xyzzyaaae38,xyzzyaaaf38
logical xyzzyaaag38
if(isperiodic)return
xyzzyaaad38=0.d0
if(nitot==0)then
xyzzyaaag38=any(pinfmass)
xyzzyaaae38=0.d0
xyzzyaaac38=0
do xyzzyaaaa38=1,nspin
if(xyzzyaaag38)then
if(.not.pinfmass(xyzzyaaaa38))then
xyzzyaaac38=xyzzyaaac38+nele(xyzzyaaaa38)
cycle
endif
xyzzyaaaf38=1.d0
else
xyzzyaaaf38=pmass(xyzzyaaaa38)
endif
xyzzyaaae38=xyzzyaaae38+nele(xyzzyaaaa38)*xyzzyaaaf38
do xyzzyaaab38=1,nele(xyzzyaaaa38)
xyzzyaaac38=xyzzyaaac38+1
xyzzyaaad38(1:dimensionality)=xyzzyaaad38(1:dimensionality)+xyzzyaaaf3&
&8*rele(1:dimensionality,xyzzyaaac38)
enddo
enddo
xyzzyaaad38(1:dimensionality)=xyzzyaaad38(1:dimensionality)/xyzzyaaae3&
&8
endif
if(xyzzyaabm18/=0)then
xyzzyaabs18=max(xyzzyaabs18,sum((rele(1:dimensionality,xyzzyaabm18)-xy&
&zzyaaad38(1:dimensionality))**2))
else
do xyzzyaaac38=1,netot
xyzzyaabs18=max(xyzzyaabs18,sum((rele(1:dimensionality,xyzzyaaac38)-xy&
&zzyaaad38(1:dimensionality))**2))
enddo
endif
end subroutine xyzzyaafw18
end subroutine vmc_main
end module slaarnacp

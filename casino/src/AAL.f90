module slaarnaal
use slaarnaaa
use slaarnaaf
use dsp
use slaarnaam
use slaarnabg
use parallel
use slaarnach
use store
use slaarnacs
use casl,          only : write_casl
use slaarnaag,     only : maximum_exp_arg,minimum_exp_arg,pi
use slaarnaah,   only : write_correlation_header
use slaarnaas,     only : mc_twist_offset,ppmcta_partial_sum,ppmcta_fi&
&t
use format_utils,  only : wout,i2s,l2s,r2s,byte2human,labelled_list,wo&
&rdwrap
use slaarnabt,     only : eigenproblem,dmatmul,dscal,make_representabl&
&e,inverse,ddot,parabolic_min,dcopy,daxpy
use slaarnaca,         only : have_ppots
use slaarnacc,only : put_random_state,get_random_state
use run_control,   only : errstop,timer,loop_time_estimate,time_report&
&,check_alloc,errstop_master,errwarn
implicit none
private
public emin_main,emin_xi_value,emin_min_energy,emin_mine_present,emin_&
&auto_varmin,emin_auto_varmin_present
real(dp) emin_xi_value,emin_min_energy
logical emin_mine_present,emin_auto_varmin,emin_auto_varmin_present
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1
real(dp) xyzzyaaad1
integer xyzzyaaae1,xyzzyaaaf1
real(dp),allocatable :: xyzzyaaag1(:),xyzzyaaah1(:)
logical,allocatable :: xyzzyaaai1(:),xyzzyaaaj1(:),xyzzyaaak1(:),xyzzy&
&aaal1(:)
character(2),allocatable :: xyzzyaaam1(:)
real(dp),allocatable :: xyzzyaaan1(:,:),xyzzyaaao1(:,:),xyzzyaaap1(:,:&
&),  xyzzyaaaq1(:,:),xyzzyaaar1(:),xyzzyaaas1(:),xyzzyaaat1(:),xyzzyaa&
&au1(:,:),xyzzyaaav1(:,:)
complex(dp),allocatable :: xyzzyaaaw1(:),xyzzyaaax1(:,:)
integer xyzzyaaay1
logical xyzzyaaaz1,xyzzyaaba1
real(dp) xyzzyaabb1,xyzzyaabc1,xyzzyaabd1
character(80) unit_string
real(dp) xyzzyaabe1,xyzzyaabf1
integer,parameter :: nsignal=10,xyzzyaabg1=1,xyzzyaabh1=2,xyzzyaabi1=3&
&,   xyzzyaabj1=4,xyzzyaabk1=5,           xyzzyaabl1=6,xyzzyaabm1=7,  &
&   xyzzyaabn1=8,xyzzyaabo1=9,xyzzyaabp1=10
character(80),parameter :: signal_text(nsignal)=(/'could not solve eig&
&enproblem',  'no valid eigenvalues        ',  'eigenvector divide by &
&zero  ',  'energy below threshold      ',  'individual weights too la&
&rge',  'weight variance too large   ',  'zero total weight           &
&',  'local energy is Not a Number',  'local energy diverges       ', &
& 'wave function is zero       '/)
logical xyzzyaabq1
logical xyzzyaabr1
logical xyzzyaabs1
real(dp) xyzzyaabt1,xyzzyaabu1
real(dp),allocatable :: xyzzyaabv1(:),xyzzyaabw1(:)
contains
subroutine emin_main()
implicit none
real(dp) xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzya&
&aaf2
call timer('EMIN',.true.)
call xyzzyaabx1(xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2&
&,xyzzyaaaf2)
call xyzzyaaca1(xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2&
&,xyzzyaaaf2)
call xyzzyaaby1
call timer('EMIN',.false.)
end subroutine emin_main
subroutine xyzzyaabx1(energy0,energy0i,errorbar0,errorbar0i,variance0,&
&variance0i)
implicit none
real(dp),intent(out) :: energy0,energy0i,errorbar0,errorbar0i,variance&
&0,variance0i
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3
real(dp) xyzzyaaae3,xyzzyaaaf3,xyzzyaaag3
real(dp),allocatable :: xyzzyaaah3(:)
logical,allocatable :: xyzzyaaai3(:),xyzzyaaaj3(:),xyzzyaaak3(:),xyzzy&
&aaal3(:),xyzzyaaam3(:)
character(2),allocatable :: xyzzyaaan3(:)
if(am_master)then
call wout('Energy minimization configuration')
call wout('=================================')
call wout()
endif
if(isperiodic)then
if(model_system)then
xyzzyaabf1=dble(netot)
unit_string='a.u./particle'
else
xyzzyaabf1=dble(npcells)
unit_string='a.u./unit cell'
endif
else
xyzzyaabf1=dble(npcells)
unit_string='a.u.'
endif
xyzzyaabe1=1.d0/xyzzyaabf1
if(emin_xi_value<0.d0.or.emin_xi_value>1.d0)call errstop_master('SETUP&
&_EMIN','The xi parameter value is invalid.')
xyzzyaabq1=.false.
xyzzyaabr1=.true.
call setup_wfn_params(xyzzyaaac3)
allocate(xyzzyaaah3(xyzzyaaac3),xyzzyaaai1(xyzzyaaac3),       xyzzyaaa&
&k3(xyzzyaaac3),xyzzyaaai3(xyzzyaaac3),   xyzzyaaal3(xyzzyaaac3),xyzzy&
&aaaj3(xyzzyaaac3),    xyzzyaaam3(xyzzyaaac3),xyzzyaaan3(xyzzyaaac3),s&
&tat=xyzzyaaad3)
call check_alloc(xyzzyaaad3,'SETUP_EMIN','1')
call get_params(xyzzyaaah3,xyzzyaaak3,xyzzyaaai3,xyzzyaaal3,xyzzyaaaj3&
&,xyzzyaaam3,xyzzyaaan3)
xyzzyaaai1(:)=.false.
if(fix_cutoffs)xyzzyaaai1(:)=xyzzyaaai1(:).or.xyzzyaaak3(:)
xyzzyaaae1=count(.not.xyzzyaaai1)
if(xyzzyaaae1<1)call errstop_master('SETUP_EMIN','No parameters to opt&
&imize!')
allocate(xyzzyaaag1(xyzzyaaae1),xyzzyaaaj1(xyzzyaaae1),xyzzyaaak1(xyzz&
&yaaae1),xyzzyaaam1(xyzzyaaae1),xyzzyaaal1(xyzzyaaae1),xyzzyaaah1(xyzz&
&yaaae1),stat=xyzzyaaad3)
call check_alloc(xyzzyaaad3,'SETUP_EMIN','2')
xyzzyaaaa3=0
do xyzzyaaab3=1,xyzzyaaac3
if(xyzzyaaai1(xyzzyaaab3))cycle
xyzzyaaaa3=xyzzyaaaa3+1
xyzzyaaag1(xyzzyaaaa3)=xyzzyaaah3(xyzzyaaab3)
xyzzyaaaj1(xyzzyaaaa3)=xyzzyaaak3(xyzzyaaab3)
xyzzyaaak1(xyzzyaaaa3)=xyzzyaaal3(xyzzyaaab3)
xyzzyaaal1(xyzzyaaaa3)=xyzzyaaam3(xyzzyaaab3)
xyzzyaaam1(xyzzyaaaa3)=xyzzyaaan3(xyzzyaaab3)
xyzzyaaah1(xyzzyaaaa3)=1.d0
enddo
if(am_master)then
call wout('Energy minimization internal setup:')
call wout(repeat('-',51))
if(xyzzyaabq1)then
call wout('Optimize                       :  variance')
call wout('Target function                :  variance')
else
call wout('Optimize                       :  energy')
call wout('Target function                :  energy + 3*error')
call wout('xi parameter for semiorthog.   :  ',emin_xi_value,rfmt='(f1&
&2.6)',adjust=.true.)
endif
if(xyzzyaabr1)then
call wout('Weights in corr. sampling      :  yes')
else
call wout('Weights in corr. sampling      :  no')
endif
call wout('Matrix regularization          :  normalized basis')
call wout('H matrix manipulation          :  on')
if(have_ppots)call wout('Fix E_NL in correl. sampling   :  '//trim(l2s&
&(opt_fixnl)))
call wout(repeat('-',51))
call wout()
call wout('There are '//trim(i2s(xyzzyaaac3))//' optimizable parameter&
&s.')
if(xyzzyaaae1/=xyzzyaaac3)then
if(count(xyzzyaaak3)==1)then
call wout(trim(i2s(count(xyzzyaaak3)))//' cut-off parameter fixed in t&
&his cycle.')
else
call wout(trim(i2s(count(xyzzyaaak3)))//' cut-off parameters fixed in &
&this cycle.')
endif
else
call wout('Will optimize all of them.')
endif
call wout()
endif
deallocate(xyzzyaaah3,xyzzyaaak3,xyzzyaaai3,xyzzyaaal3,xyzzyaaaj3,xyzz&
&yaaam3,xyzzyaaan3)
xyzzyaaay1=0
call scratch_protect(xyzzyaaay1)
call energy_scratch_request(xyzzyaaay1)
call scratch_request(wfn_detail=xyzzyaaay1)
call ederiv_scratch_request(xyzzyaaay1,xyzzyaaai1,wfn_deriv=.true.)
call setup_scratch
call which_scratch(xyzzyaaay1)
call setup_wfn_utils
if(use_altsamp)call setup_alt_utils
call setup_energy_utils
call xyzzyaabz1(energy0,energy0i,errorbar0,errorbar0i,variance0,varian&
&ce0i)
xyzzyaaaf1=xyzzyaaae1+1
allocate(xyzzyaaan1(xyzzyaaaf1,xyzzyaaaf1),xyzzyaaao1(xyzzyaaaf1,xyzzy&
&aaaf1),xyzzyaaap1(xyzzyaaaf1,xyzzyaaaf1),xyzzyaaaq1(xyzzyaaaf1,xyzzya&
&aaf1),xyzzyaaax1(xyzzyaaaf1,xyzzyaaaf1),stat=xyzzyaaad3)
call check_alloc(xyzzyaaad3,'SETUP_EMIN','global')
if(complex_wf)then
allocate(xyzzyaaau1(xyzzyaaaf1,xyzzyaaaf1),xyzzyaaav1(xyzzyaaaf1,xyzzy&
&aaaf1),stat=xyzzyaaad3)
call check_alloc(xyzzyaaad3,'SETUP_EMIN','complex')
endif
if(am_master)then
xyzzyaaae3=dble(size(rele_config))
xyzzyaaag3=2.d0*dble(size(xyzzyaaaw1))+dble(size(xyzzyaaar1))+dble(siz&
&e(xyzzyaaat1))
xyzzyaaaf3=dble(size(xyzzyaaan1))+dble(size(xyzzyaaao1))+dble(size(xyz&
&zyaaap1))+dble(size(xyzzyaaaq1))+2.d0*dble(size(xyzzyaaax1))
if(complex_wf)then
xyzzyaaag3=xyzzyaaag3+dble(size(xyzzyaaas1))
xyzzyaaaf3=xyzzyaaaf3+dble(size(xyzzyaaau1))+dble(size(xyzzyaaav1))
endif
call wout('Optimization workspace:')
call wout(repeat('-',51))
call wout('No. of variable parameters (P) :  ',trim(i2s(xyzzyaaae1)),f&
&mt='(1x,a,a18)')
call wout('No. of configurations (C)      :  ',trim(i2s(xyzzyaaab1)),f&
&mt='(1x,a,a18)')
call wout(repeat('-',51))
call wout('Configuration storage          :  ',trim(byte2human(8.d0*xy&
&zzyaaae3)),fmt='(1x,a,a18)')
call wout('Vectors of size C              :  ',trim(byte2human(8.d0*xy&
&zzyaaag3)),fmt='(1x,a,a18)')
call wout('Matrices of size P^2           :  ',trim(byte2human(8.d0*xy&
&zzyaaaf3)),fmt='(1x,a,a18)')
call wout(repeat('-',51))
if(nnodes==1)then
call wout('Total memory required          :  ',trim(byte2human(8.d0*(x&
&yzzyaaae3+xyzzyaaag3+xyzzyaaaf3))),fmt='(1x,a,a18)')
else
call wout('Total memory required per node :  ',trim(byte2human(8.d0*(x&
&yzzyaaae3+xyzzyaaag3+xyzzyaaaf3))),fmt='(1x,a,a18)')
endif
call wout(repeat('-',51))
call wout()
endif
end subroutine xyzzyaabx1
subroutine xyzzyaaby1
implicit none
logical :: xyzzyaaaa4=.true.
character(20) dum_config_item(0),extra_item(1)
if(.not.rng_restart_safe)xyzzyaaaa4=.false.
call dismantle_configs
extra_item(1)='RANDOM'
call init_config_accumulation('OPT',0,dum_config_item,extra_item)
call get_random_state(random_state_config)
call end_config_accumulation(.false.)
select case(chkpoint_level)
case(-1)
continue
case(0)
if(isopt.or.((isvmc_opt.or.isopt_vmc).and.opt_cycle==opt_cycles))call &
&write_configs(xyzzyaaaa4)
case(1,2)
call write_configs(xyzzyaaaa4)
end select
call finish_storage_energy
call finish_storage
deallocate(xyzzyaaaw1)
deallocate(xyzzyaaar1,xyzzyaaat1)
if(complex_wf)deallocate(xyzzyaaas1)
if(xyzzyaabs1)deallocate(xyzzyaabv1,xyzzyaabw1)
call finish_energy_utils
call finish_wfn_utils
if(use_altsamp)call finish_alt_utils
call finish_scratch
xyzzyaaay1=0
call finish_wfn_params
deallocate(xyzzyaaag1,xyzzyaaaj1,xyzzyaaak1,xyzzyaaal1,xyzzyaaam1,xyzz&
&yaaai1,xyzzyaaah1)
deallocate(xyzzyaaan1,xyzzyaaao1,xyzzyaaap1,xyzzyaaaq1,xyzzyaaax1)
if(complex_wf)deallocate(xyzzyaaau1,xyzzyaaav1)
end subroutine xyzzyaaby1
subroutine xyzzyaabz1(energy0,energy0i,errorbar0,errorbar0i,variance0,&
&variance0i)
implicit none
real(dp),intent(out) :: energy0,energy0i,errorbar0,errorbar0i,variance&
&0,variance0i
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5,xyzzyaaag5,xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5
real(dp) xyzzyaaak5,xyzzyaaal5,xyzzyaaam5,xyzzyaaan5,xyzzyaaao5,xyzzya&
&aap5,xyzzyaaaq5,xyzzyaaar5,xyzzyaaas5,xyzzyaaat5,xyzzyaaau5(9),xyzzya&
&aav5(9)
logical xyzzyaaaw5(1),xyzzyaaax5(1),isnan,isinf,xyzzyaaay5,nan,inf,xyz&
&zyaaaz5
character(20) opt_item(1),req_extra(0),opt_extra(1)
character(20),allocatable :: xyzzyaaba5(:)
character(80) tmpr
xyzzyaaaa1=-1
xyzzyaaab5=1
if(noncoll_spin)xyzzyaaab5=xyzzyaaab5+1
if(use_altsamp)xyzzyaaab5=xyzzyaaab5+1
allocate(xyzzyaaba5(xyzzyaaab5),stat=xyzzyaaad5)
xyzzyaaac5=1
xyzzyaaba5(1:1)=(/'RELE'/)
if(noncoll_spin)then
xyzzyaaac5=xyzzyaaac5+1
xyzzyaaba5(xyzzyaaac5)='SELE'
endif
if(use_altsamp)then
xyzzyaaac5=xyzzyaaac5+1
xyzzyaaba5(xyzzyaaac5)='LOGP'
endif
opt_item=(/'TWIST'/)
opt_extra=(/'RANDOM'/)
call load_configs(xyzzyaaaa1,'VMC',xyzzyaaba5,opt_item,xyzzyaaaw5,req_&
&extra,opt_extra,xyzzyaaax5)
deallocate(xyzzyaaba5)
xyzzyaabs1=xyzzyaaaw5(1)
if(xyzzyaaax5(1))call put_random_state(random_state_config)
xyzzyaaaz1=virtual_nnodes>1
xyzzyaaab1=xyzzyaaaa1*nnodes
xyzzyaaac1=xyzzyaaab1
xyzzyaaba1=am_master.and.(.not.xyzzyaaaz1.or.virtual_node==0)
if(xyzzyaaaz1)then
xyzzyaaac1=virtual_nconfig
if(virtual_nconfig<xyzzyaaaa1)call errstop_master('EMIN_READ_CONFIGS',&
&'VIRTUAL_NCONFIG < number of configs actually present.')
endif
if(xyzzyaaac1<xyzzyaaae1)call errstop_master('EMIN_READ_CONFIGS','Fewe&
&r configs than parameters in config file.')
if(xyzzyaaac1<2)call errstop_master('EMIN_READ_CONFIGS','Cannot perfor&
&m EMIN with less than two configs.')
xyzzyaaad1=dble(xyzzyaaac1)
xyzzyaabc1=0.05d0*xyzzyaaad1
xyzzyaabd1=0.3d0**2
allocate(xyzzyaaaw1(xyzzyaaaa1),stat=xyzzyaaad5)
call check_alloc(xyzzyaaad5,'EMIN_READ_CONFIGS','logwfn_original')
xyzzyaaaw1=cmplx(0.d0,0.d0,dp)
if(xyzzyaaba1)then
allocate(xyzzyaaar1(xyzzyaaac1),xyzzyaaat1(xyzzyaaac1),stat=xyzzyaaad5&
&)
elseif(am_master)then
allocate(xyzzyaaar1(xyzzyaaab1),xyzzyaaat1(xyzzyaaab1),stat=xyzzyaaad5&
&)
else
allocate(xyzzyaaar1(xyzzyaaaa1),xyzzyaaat1(xyzzyaaaa1),stat=xyzzyaaad5&
&)
endif
call check_alloc(xyzzyaaad5,'EMIN_READ_CONFIGS','cs')
xyzzyaaar1=0.d0
xyzzyaaat1=0.d0
if(complex_wf)then
if(xyzzyaaba1)then
allocate(xyzzyaaas1(xyzzyaaac1),stat=xyzzyaaad5)
elseif(am_master)then
allocate(xyzzyaaas1(xyzzyaaab1),stat=xyzzyaaad5)
else
allocate(xyzzyaaas1(xyzzyaaaa1),stat=xyzzyaaad5)
endif
call check_alloc(xyzzyaaad5,'EMIN_READ_CONFIGS','csim')
xyzzyaaas1=0.d0
endif
call setup_storage(xyzzyaaaa1,xyzzyaaai1)
call setup_storage_energy(xyzzyaaaa1)
isnan=.false.
isinf=.false.
xyzzyaaay5=.false.
do xyzzyaaaa5=1,xyzzyaaaa1
call load_from_storage(xyzzyaaay1,xyzzyaaaa5)
call load_from_storage_energy(xyzzyaaay1,xyzzyaaaa5)
if(xyzzyaabs1)call mc_twist_offset(twist_config(1:3,xyzzyaaaa5))
call eval_local_energy(xyzzyaaay1,etot=xyzzyaaar1(xyzzyaaaa5),etot_ima&
&g=xyzzyaaal5,fix_nl_grid=opt_fixnl,isnan=nan,isinf=inf)
isnan=isnan.or.nan
isinf=isinf.or.inf
if(complex_wf)xyzzyaaas1(xyzzyaaaa5)=xyzzyaaal5
call wfn_logval(xyzzyaaay1,xyzzyaaaw1(xyzzyaaaa5),xyzzyaaaz5)
xyzzyaaay5=xyzzyaaay5.or.xyzzyaaaz5
call save_to_storage(xyzzyaaay1,xyzzyaaaa5)
call save_to_storage_energy(xyzzyaaay1,xyzzyaaaa5)
if(xyzzyaabr1.and.use_altsamp)xyzzyaaat1(xyzzyaaaa5)=2.d0*dble(xyzzyaa&
&aw1(xyzzyaaaa5)-cmplx(logp_config(xyzzyaaaa5),0.d0,dp))
enddo
if(isnan.or.isinf.or.xyzzyaaay5)call errstop('EMIN_READ_CONFIGS','Floa&
&ting-point exception reported when calculating initial energies. Don'&
&'t know what to do, so stopping. Please file a bug report.')
if(xyzzyaabs1)then
allocate(xyzzyaabv1(xyzzyaaaa1),xyzzyaabw1(xyzzyaaaa1),stat=xyzzyaaad5&
&)
call check_alloc(xyzzyaaad5,'EMIN_READ_CONFIGS','opt_hfke')
xyzzyaabv1(1:xyzzyaaaa1)=twist_config(4,1:xyzzyaaaa1)*dble(netot)
xyzzyaabw1(1:xyzzyaaaa1)=twist_config(5,1:xyzzyaaaa1)*dble(netot)
call ppmcta_partial_sum(xyzzyaaaa1,xyzzyaaar1(1:xyzzyaaaa1),xyzzyaabv1&
&,xyzzyaabw1,xyzzyaaau5)
call mpi_reduce(xyzzyaaau5,xyzzyaaav5,9,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
if(xyzzyaaaz1.and.am_master)then
if(xyzzyaaba1)then
do xyzzyaaae5=1,virtual_nnodes-1
write(6,*)xyzzyaaae5
read(5,*,iostat=ierror)(xyzzyaaau5(xyzzyaaaj5),xyzzyaaaj5=1,9)
if(ierror/=0)call errstop('EMIN_READ_CONFIGS','VP error (twist_partial&
&_sum).')
xyzzyaaav5=xyzzyaaav5+xyzzyaaau5
enddo
else
do xyzzyaaae5=1,virtual_nnodes-1
read(5,*,iostat=ierror)xyzzyaaaf5
if(ierror/=0)call errstop('SETUP_VARMIN','VP error (jnode).')
if(xyzzyaaaf5/=virtual_node)cycle
write(6,*)(xyzzyaaav5(xyzzyaaaj5),xyzzyaaaj5=1,9)
enddo
endif
endif
if(xyzzyaaba1)call ppmcta_fit(xyzzyaaav5,xyzzyaabt1,xyzzyaabu1)
if(xyzzyaaaz1.and.am_master)then
if(xyzzyaaba1)then
write(6,*)xyzzyaabt1,xyzzyaabu1
else
read(5,*,iostat=ierror)xyzzyaabt1,xyzzyaabu1
if(ierror/=0)call errstop('SETUP_VARMIN','VP error (twist_a).')
endif
endif
call mpi_bcast(xyzzyaabt1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcasting twist_a in emin_read_cfgs')
call mpi_bcast(xyzzyaabu1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcasting twist_b in emin_read_cfgs')
if(am_master)then
call wout('MCTA post-processing fit parameters for this run:')
call wout(' a = ',xyzzyaabt1)
call wout(' b = ',xyzzyaabu1)
call wout()
endif
xyzzyaaar1(1:xyzzyaaaa1)=xyzzyaaar1(1:xyzzyaaaa1)+xyzzyaabt1*xyzzyaabv&
&1(1:xyzzyaaaa1)+xyzzyaabu1*xyzzyaabw1(1:xyzzyaaaa1)
endif
energy0=0.d0
energy0i=0.d0
errorbar0=0.d0
errorbar0i=0.d0
variance0=0.d0
variance0i=0.d0
if(.not.use_altsamp)then
xyzzyaaam5=sum(xyzzyaaar1(1:xyzzyaaaa1))
xyzzyaaan5=sum(xyzzyaaar1(1:xyzzyaaaa1)**2)
if(complex_wf)then
xyzzyaaao5=sum(xyzzyaaas1(1:xyzzyaaaa1))
xyzzyaaap5=sum(xyzzyaaas1(1:xyzzyaaaa1)**2)
endif
call mpi_reduce(xyzzyaaam5,xyzzyaaaq5,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing sum_E in emin_read_cfgs')
call mpi_reduce(xyzzyaaan5,xyzzyaaar5,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing sum_E2 in emin_read_cfgs')
if(complex_wf)then
call mpi_reduce(xyzzyaaao5,xyzzyaaas5,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing sum_Ei in emin_read_cfgs')
call mpi_reduce(xyzzyaaap5,xyzzyaaat5,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'reducing sum_Ei2 in emin_read_cfgs')
endif
if(xyzzyaaaz1.and.am_master)then
call time_report('VP communication: gathering energies',.true.)
if(xyzzyaaba1)then
do xyzzyaaae5=1,virtual_nnodes-1
write(6,*)xyzzyaaae5
read(5,*,iostat=xyzzyaaag5)xyzzyaaam5,xyzzyaaan5
if(xyzzyaaag5/=0)call errstop('EMIN_READ_CONFIG','VP error (sum_E).')
xyzzyaaaq5=xyzzyaaaq5+xyzzyaaam5
xyzzyaaar5=xyzzyaaar5+xyzzyaaan5
if(complex_wf)then
read(5,*,iostat=xyzzyaaag5)xyzzyaaao5,xyzzyaaap5
if(xyzzyaaag5/=0)call errstop('EMIN_READ_CONFIGS','VP error (sum_Ei).'&
&)
xyzzyaaas5=xyzzyaaas5+xyzzyaaao5
xyzzyaaat5=xyzzyaaat5+xyzzyaaap5
endif
enddo
else
do xyzzyaaae5=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaag5)xyzzyaaaf5
if(xyzzyaaag5/=0)call errstop('EMIN_READ_CONFIGS','VP error (jnode).')
if(xyzzyaaaf5/=virtual_node)cycle
write(6,*)xyzzyaaaq5,xyzzyaaar5
if(complex_wf)write(6,*)xyzzyaaas5,xyzzyaaat5
enddo
endif
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
if(xyzzyaaba1)then
energy0=xyzzyaaaq5/xyzzyaaad1
variance0=xyzzyaaar5/(xyzzyaaad1-1.d0)-(xyzzyaaad1/(xyzzyaaad1-1.d0))*&
&energy0**2
variance0=max(0.d0,variance0)
errorbar0=sqrt(variance0/xyzzyaaad1)
if(complex_wf)then
energy0i=xyzzyaaas5/xyzzyaaad1
variance0i=xyzzyaaat5/(xyzzyaaad1-1.d0)-(xyzzyaaad1/(xyzzyaaad1-1.d0))&
&*energy0i**2
variance0i=max(0.d0,variance0i)
errorbar0i=sqrt(variance0i/xyzzyaaad1)
endif
endif
if(xyzzyaaaz1.and.am_master)then
call time_report('VP communication: broadcasting initial energies',.tr&
&ue.)
if(xyzzyaaba1)then
write(6,*)energy0,errorbar0,variance0
if(complex_wf)write(6,*)energy0i,errorbar0i,variance0i
else
read(5,*,iostat=xyzzyaaag5)energy0,errorbar0,variance0
if(xyzzyaaag5/=0)call errstop('EMIN_READ_CONFIGS','VP error (energy0).&
&')
if(complex_wf)then
read(5,*,iostat=xyzzyaaag5)energy0i,errorbar0i,variance0i
if(xyzzyaaag5/=0)call errstop('EMIN_READ_CONFIGS','VP error (energy0i)&
&.')
endif
endif
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(energy0,1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting energy0 in emin_read_cfgs')
call mpi_bcast(errorbar0,1,mpi_double_precision,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'broadcasting errorbar0 in emin_read_cfgs')
call mpi_bcast(variance0,1,mpi_double_precision,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'broadcasting variance0 in emin_read_cfgs')
if(complex_wf)then
call mpi_bcast(energy0i,1,mpi_double_precision,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'broadcasting energy0i in emin_read_cfgs')
call mpi_bcast(errorbar0i,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcasting errorbar0i in emin_read_cfgs')
call mpi_bcast(variance0i,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcasting variance0i in emin_read_cfgs')
endif
else
call mpi_reduce(isnan,nan,1,mpi_logical,mpi_lor,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'reducing isNaN is emin_read_configs')
call mpi_reduce(isinf,inf,1,mpi_logical,mpi_lor,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'reducing isInf is emin_read_configs')
call mpi_reduce(xyzzyaaay5,xyzzyaaaz5,1,mpi_logical,mpi_lor,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'reducing isZero is emin_read_configs')
if(am_master)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaar1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_energies in emin_read_configs')
if(complex_wf)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaas1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_energies_imag in emin_read_configs'&
&)
endif
if(xyzzyaabr1)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaat1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_weights in emin_read_configs')
endif
else
call mpi_gather(xyzzyaaar1,xyzzyaaaa1,mpi_double_precision,xyzzyaaar1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_energies in emin_read_configs')
if(complex_wf)then
call mpi_gather(xyzzyaaas1,xyzzyaaaa1,mpi_double_precision,xyzzyaaas1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_energies_imag in emin_read_confgis'&
&)
endif
if(xyzzyaabr1)then
call mpi_gather(xyzzyaaat1,xyzzyaaaa1,mpi_double_precision,xyzzyaaat1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_weights in emin_read_configs')
endif
endif
if(xyzzyaaaz1.and.am_master)then
call time_report(' VP communication: gathering energies',.true.)
if(xyzzyaaba1)then
xyzzyaaah5=xyzzyaaab1
do xyzzyaaae5=1,virtual_nnodes-1
write(6,*)xyzzyaaae5
read(5,*,iostat=ierror)nan,inf,xyzzyaaaz5
if(ierror/=0)call errstop('EMIN_READ_CONFIGS','VP error (NaN).')
if(nan)isnan=.true.
if(inf)isinf=.true.
if(xyzzyaaaz5)xyzzyaaay5=.true.
read(5,*,iostat=ierror)xyzzyaaai5
if(ierror/=0)call errstop('EMIN_READ_CONFIGS','VP error (kcfg).')
read(5,*,iostat=ierror)(xyzzyaaar1(xyzzyaaaa5),xyzzyaaaa5=xyzzyaaah5+1&
&,xyzzyaaah5+xyzzyaaai5)
if(ierror/=0)call errstop('EMIN_READ_CONFIGS','VP error (cs_energies).&
&')
if(complex_wf)then
read(5,*,iostat=ierror)(xyzzyaaas1(xyzzyaaaa5),xyzzyaaaa5=xyzzyaaah5+1&
&,xyzzyaaah5+xyzzyaaai5)
if(ierror/=0)call errstop('EMIN_READ_CONFIGS','VP error (cs_energies_i&
&mag).')
endif
if(xyzzyaabr1)then
read(5,*,iostat=ierror)(xyzzyaaat1(xyzzyaaaa5),xyzzyaaaa5=xyzzyaaah5+1&
&,xyzzyaaah5+xyzzyaaai5)
if(ierror/=0)call errstop('EMIN_READ_CONFIGS','VP error (cs_weights).'&
&)
endif
xyzzyaaah5=xyzzyaaah5+xyzzyaaai5
enddo
if(xyzzyaaah5/=virtual_nconfig)call errstop('EMIN_READ_CONFIGS','VP er&
&ror (config count '//trim(i2s(xyzzyaaah5))//'/='//trim(i2s(virtual_nc&
&onfig))//').')
else
do xyzzyaaae5=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaag5)xyzzyaaaf5
if(xyzzyaaag5/=0)call errstop('EMIN_READ_CONFIGS','VP error (jnode).')
if(xyzzyaaaf5/=virtual_node)cycle
write(6,*)isnan,isinf,xyzzyaaay5
write(6,*)xyzzyaaab1
write(6,*)(xyzzyaaar1(xyzzyaaaa5),xyzzyaaaa5=1,xyzzyaaab1)
if(complex_wf)write(6,*)(xyzzyaaas1(xyzzyaaaa5),xyzzyaaaa5=1,xyzzyaaab&
&1)
if(xyzzyaabr1)write(6,*)(xyzzyaaat1(xyzzyaaaa5),xyzzyaaaa5=1,xyzzyaaab&
&1)
enddo
endif
call time_report(' Done.',report_walltime=.true.,no_newline=.true.)
endif
if(xyzzyaaba1)then
xyzzyaaal5=maxval(xyzzyaaat1)
xyzzyaaat1=xyzzyaaat1-xyzzyaaal5
do xyzzyaaaa5=1,xyzzyaaac1
if(xyzzyaaat1(xyzzyaaaa5)>maximum_exp_arg*0.8d0)then
xyzzyaaat1(xyzzyaaaa5)=0.d0
else
xyzzyaaat1(xyzzyaaaa5)=exp(xyzzyaaat1(xyzzyaaaa5))
endif
enddo
xyzzyaaal5=sum(xyzzyaaat1)
call dscal(xyzzyaaac1,xyzzyaaad1/xyzzyaaal5,xyzzyaaat1(1),1)
xyzzyaaal5=sum((xyzzyaaat1(:)-1.d0)**2)/(xyzzyaaad1*(xyzzyaaad1-1.d0))
energy0=sum(xyzzyaaat1(:)*xyzzyaaar1(:))/xyzzyaaad1
variance0=sum((xyzzyaaat1(:)*(xyzzyaaar1(:)-energy0))**2)/(xyzzyaaad1-&
&1.d0)
variance0=max(0.d0,variance0)
errorbar0=sqrt(variance0/xyzzyaaad1)
if(complex_wf)then
energy0i=sum(xyzzyaaat1(:)*xyzzyaaas1(:))/xyzzyaaad1
variance0i=sum((xyzzyaaat1(:)*(xyzzyaaas1(:)-energy0i))**2)/(xyzzyaaad&
&1-1.d0)
variance0i=max(0.d0,variance0i)
errorbar0i=sqrt(variance0i/xyzzyaaad1)
endif
endif
if(xyzzyaaaz1.and.am_master)then
call time_report(' VP communication: broadcasting energy',.true.)
if(xyzzyaaba1)then
write(6,*)energy0,errorbar0,variance0
if(complex_wf)write(6,*)energy0i,errorbar0i,variance0i
else
read(5,*,iostat=xyzzyaaag5)energy0,errorbar0,variance0
if(xyzzyaaag5/=0)call errstop('EMIN_READ_CONFIGS','VP error (energy).'&
&)
if(complex_wf)then
read(5,*,iostat=xyzzyaaag5)energy0i,errorbar0i,variance0i
if(xyzzyaaag5/=0)call errstop('EMIN_READ_CONFIGS','VP error (energyi).&
&')
endif
endif
call time_report(' Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(energy0,1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting energy in emin_read_configs.')
call mpi_bcast(errorbar0,1,mpi_double_precision,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'broadcasting errorbar in emin_read_configs.')
call mpi_bcast(variance0,1,mpi_double_precision,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'broadcasting variance in emin_read_configs.')
if(complex_wf)then
call mpi_bcast(energy0i,1,mpi_double_precision,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'broadcasting energyi in emin_read_configs.')
call mpi_bcast(errorbar0i,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcasting errorbari in emin_read_configs.')
call mpi_bcast(variance0i,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcasting variancei in emin_read_configs.')
endif
endif
xyzzyaaak5=energy0-4.d0*sqrt(variance0)
if(emin_mine_present)then
xyzzyaabb1=emin_min_energy*xyzzyaabf1
else
xyzzyaabb1=xyzzyaaak5
endif
if(am_master)then
if(.not.xyzzyaaaz1)then
call wout('Number of nodes                                       : '//&
&trim(i2s(nnodes)))
if(nnodes>1)call wout('Number of configurations per node              &
&       : '//trim(i2s(xyzzyaaaa1)))
call wout('Total number of configurations                        : '//&
&trim(i2s(xyzzyaaac1)))
else
if(xyzzyaaba1)then
call wout('VP optimization: this is the VP master node.')
else
call wout('VP optimization: this is VP slave node #'//trim(i2s(virtual&
&_node))//'.')
endif
if(nnodes>1)then
call wout('Number of MPI nodes on this VP node                  : '//t&
&rim(i2s(nnodes)))
call wout('Number of configurations per MPI node on this VP node: '//t&
&rim(i2s(xyzzyaaaa1)))
endif
call wout('Number of configurations on this VP node              : '//&
&trim(i2s(xyzzyaaab1)))
call wout('Number of VP nodes                                    : '//&
&trim(i2s(virtual_nnodes)))
call wout('Total number of configurations                        : '//&
&trim(i2s(xyzzyaaac1)))
endif
call wout()
write(tmpr,*)xyzzyaabe1*xyzzyaaak5
call wordwrap('Stored VMC result suggests minimum energy of '//trim(ad&
&justl(tmpr))//' '//trim(unit_string))
if(emin_mine_present)then
write(tmpr,*)xyzzyaabe1*xyzzyaabb1
call wordwrap('User-supplied minimum of '//trim(adjustl(tmpr))//' '//t&
&rim(unit_string)//' will be enforced instead.')
else
call wout('This minimum will be enforced.')
endif
call wout()
endif
end subroutine xyzzyaabz1
subroutine xyzzyaaca1(energy,energyi,errorbar,errorbari,variance,varia&
&ncei)
implicit none
real(dp),intent(inout) :: energy,energyi,errorbar,errorbari,variance,v&
&ariancei
character(80) tmpr,tmpr2,tmpr3,tmpr4
character(512) errmsg
integer xyzzyaaaa6,xyzzyaaab6
real(dp) xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaaaf6,xyzzyaaag6,xyzzya&
&aah6,params_print(xyzzyaaae1),xyzzyaaai6(xyzzyaaaf1),xyzzyaaaj6,xyzzy&
&aaak6,criterion,xyzzyaaal6
real(dp) xyzzyaaam6,xyzzyaaan6
if(am_master)then
call wout('Optimization start')
call wout('==================')
call wout()
if(xyzzyaaaz1.and..not.xyzzyaaba1)then
call wout('STARTING VIRTUAL SLAVE OPTIMIZATION')
call wout('===================================')
call wout()
endif
endif
call xyzzyaacp1(energy,energyi,errorbar,errorbari,variance,variancei,x&
&yzzyaaaj6)
call put_params(xyzzyaaag1,xyzzyaaai1,.false.,params_print=params_prin&
&t)
if(xyzzyaaba1.and.opt_info>1)call xyzzyaacb1(0,energy,energyi,errorbar&
&,errorbari,variance,variancei,params_print)
do xyzzyaaaa6=1,opt_maxiter
xyzzyaaak6=xyzzyaaaj6
xyzzyaaac6=energy
xyzzyaaad6=energyi
xyzzyaaae6=errorbar
xyzzyaaaf6=errorbari
xyzzyaaag6=variance
xyzzyaaah6=variancei
energy=0.d0
energyi=0.d0
errorbar=-1.d0
errorbari=0.d0
variance=-1.d0
variancei=0.d0
call xyzzyaacc1(xyzzyaaac6,xyzzyaaad6)
call xyzzyaacf1(xyzzyaaai6,xyzzyaaal6,xyzzyaaab6)
call qmc_barrier
if(xyzzyaaab6/=0)then
if(am_master)then
call wordwrap('Exiting optimization after failing in iteration #'//tri&
&m(i2s(xyzzyaaaa6))//'.')
call wout()
endif
exit
endif
call xyzzyaaci1(xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaaaf6,xyzzyaaag6&
&,xyzzyaaah6,xyzzyaaai6,energy,energyi,errorbar,errorbari,variance,var&
&iancei)
xyzzyaaag1(:)=xyzzyaaag1(:)+xyzzyaaai6(2:xyzzyaaaf1)*xyzzyaaah1(:)
call put_params(xyzzyaaag1,xyzzyaaai1,.false.,params_print=params_prin&
&t)
if(am_master)then
if(xyzzyaaba1.and.opt_info>1)call xyzzyaacb1(xyzzyaaaa6,energy,energyi&
&,errorbar,errorbari,variance,variancei,params_print)
call write_correlation_header(.true.,.true.)
call update_wfn_casl
call write_casl(':parameters.casl','parameters.'//trim(i2s(max(1,opt_c&
&ycle)))//'.casl',errmsg)
if(len_trim(errmsg)>0)call errstop_master('EMIN_DRIVER',trim(errmsg))
endif
call xyzzyaacp1(energy,energyi,errorbar,errorbari,variance,variancei,x&
&yzzyaaaj6,criterion=criterion)
if(xyzzyaaak6-xyzzyaaaj6<=criterion)then
if(xyzzyaaba1.and.opt_info>1)then
call wout('Criterion for convergence satisfied.')
call wout()
endif
exit
endif
if(xyzzyaaaa6==opt_maxiter.and.xyzzyaaba1.and.opt_info>1)then
call wout('Exceeded maximum number of iterations.')
call wout()
endif
enddo
if(use_altsamp.and.altsamp==1)then
xyzzyaaam6=energy
xyzzyaaan6=2.d0*errorbar
if(am_master)then
call wout()
if(xyzzyaaam6>vmc_optimum_e0.or.xyzzyaaan6>vmc_optimum_ew.or.opt_cycle&
&==1)then
call wout(' Unmodified PDF parameters.')
xyzzyaaam6=vmc_optimum_e0
xyzzyaaan6=vmc_optimum_ew
else
call wout(' Modified PDF parameters.')
endif
if(simplepdf==1)xyzzyaaan6=vmc_optimum_ew
tmpr=r2s(vmc_optimum_e0,'(f10.4)')
tmpr2=r2s(vmc_optimum_ew,'(f10.4)')
tmpr3=r2s(xyzzyaaam6,'(f10.4)')
tmpr4=r2s(xyzzyaaan6,'(f10.4)')
call wout('(e0,ew)_old -> (e0,ew)_new = ('//trim(tmpr)//','//trim(tmpr&
&2)//') -> ('//trim(tmpr3)//','//trim(tmpr4)//')')
call wout()
vmc_optimum_e0=xyzzyaaam6
vmc_optimum_ew=xyzzyaaan6
endif
call mpi_bcast(vmc_optimum_e0,1,mpi_double_precision,0,mpi_comm_world,&
&ierror)
call mpi_bcast(vmc_optimum_ew,1,mpi_double_precision,0,mpi_comm_world,&
&ierror)
endif
if(xyzzyaaaz1.and..not.xyzzyaaba1)then
call wout('VP slave optimization ended.')
call wout()
endif
end subroutine xyzzyaaca1
subroutine xyzzyaacb1(iter,energy,energyi,errorbar,errorbari,variance,&
&variancei,params_print)
implicit none
integer,intent(in) :: iter
real(dp),intent(in) :: energy,energyi,errorbar,errorbari,variance,vari&
&ancei,params_print(xyzzyaaae1)
integer xyzzyaaaa7
character(80) tmpr,char_80
xyzzyaaaa7=34
if(complex_wf)xyzzyaaaa7=25
if(.not.am_master)return
call wout('Optimization monitor :')
call wout('----------------------')
if(iter==0)then
call wout('Start of minimization process')
call wout('Parameters:')
call labelled_list(params_print,xyzzyaaam1,flag=xyzzyaaaj1,comment='sh&
&allow parameters')
else
call wout('After iteration : '//trim(i2s(iter)))
call wout('Parameters:')
call labelled_list(params_print,xyzzyaaam1)
endif
call wout()
if(errorbar>=0.d0.or.variance<0.d0)then
char_80='Energy ('//trim(unit_string)//')'
char_80=repeat(' ',xyzzyaaaa7-len_trim(char_80))//trim(char_80)
write(tmpr,*)energy*xyzzyaabe1
char_80=trim(char_80)//' : '//trim(adjustl(tmpr))
if(complex_wf)then
write(tmpr,*)abs(energyi*xyzzyaabe1)
if(energyi<0.d0)then
char_80=trim(char_80)//' - i '//trim(adjustl(tmpr))
else
char_80=trim(char_80)//' + i '//trim(adjustl(tmpr))
endif
endif
call wout(trim(char_80))
endif
if(errorbar>=0.d0)then
char_80='Error ('//trim(unit_string)//')'
char_80=repeat(' ',xyzzyaaaa7-len_trim(char_80))//trim(char_80)
write(tmpr,*)errorbar*xyzzyaabe1
char_80=trim(char_80)//' : '//trim(adjustl(tmpr))
if(complex_wf)then
write(tmpr,*)errorbari*xyzzyaabe1
char_80=trim(char_80)//' + i '//trim(adjustl(tmpr))
endif
call wout(trim(char_80))
endif
if(variance>=0.d0)then
char_80='Variance (a.u.)'
char_80=repeat(' ',xyzzyaaaa7-len_trim(char_80))//trim(char_80)
write(tmpr,*)variance
char_80=trim(char_80)//' : '//trim(adjustl(tmpr))
if(complex_wf)then
write(tmpr,*)variancei
char_80=trim(char_80)//' + i '//trim(adjustl(tmpr))
endif
call wout(trim(char_80))
endif
call wout()
end subroutine xyzzyaacb1
subroutine xyzzyaacc1(energy0,energy0i)
implicit none
real(dp),intent(in) :: energy0,energy0i
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8
real(dp) xyzzyaaad8,ei0,xyzzyaaae8(xyzzyaaae1),xyzzyaaaf8(xyzzyaaaf1),&
&xyzzyaaag8(xyzzyaaaf1),xyzzyaaah8(xyzzyaaaf1),xyzzyaaai8(xyzzyaaaf1),&
&xyzzyaaaj8(xyzzyaaaf1),xyzzyaaak8(xyzzyaaaf1),xyzzyaaal8,xyzzyaaam8,x&
&yzzyaaan8,xyzzyaaao8,xyzzyaaap8,xyzzyaaaq8
complex(dp) xyzzyaaar8
logical xyzzyaaas8(xyzzyaaae1),isnan,isinf,iszero
call timer('MATRIX_GENERATION',.true.)
xyzzyaaae8=0.d0
xyzzyaaas8=.false.
xyzzyaaan1=0.d0
xyzzyaaap1=0.d0
if(complex_wf)then
xyzzyaaau1=0.d0
xyzzyaaav1=0.d0
endif
do xyzzyaaaa8=1,xyzzyaaaa1
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Computing deriva&
&tives.',xyzzyaaaa8,xyzzyaaaa1)
call load_from_storage(xyzzyaaay1,xyzzyaaaa8)
call load_from_storage_energy(xyzzyaaay1,xyzzyaaaa8)
call prepare_ederivs(xyzzyaaay1,xyzzyaaai1,xyzzyaaad8,xyzzyaaar8,get_l&
&ogwfn=.true.,ei0=ei0,isnan=isnan,isinf=isinf,iszero=iszero)
if(isnan.or.isinf)call errstop('EMIN_MATRIX_GEN','Floating-point excep&
&tion reported when calculating energies. Don''t know what to do, so s&
&topping. Please file a bug report.')
if(iszero)call errstop('EMIN_MATRIX_GEN','Zero wave function encounter&
&ed. Don''t know what to do, so stopping. Please file a bug report.')
call save_to_storage(xyzzyaaay1,xyzzyaaaa8)
call save_to_storage_energy(xyzzyaaay1,xyzzyaaaa8)
xyzzyaaah8(1)=1.d0
xyzzyaaai8(1)=0.d0
xyzzyaaaf8(1)=0.d0
xyzzyaaag8(1)=0.d0
do xyzzyaaab8=1,xyzzyaaae1
xyzzyaaac8=xyzzyaaab8+1
call xyzzyaacd1(xyzzyaaaa8,xyzzyaaab8,xyzzyaaah8(xyzzyaaac8),xyzzyaaai&
&8(xyzzyaaac8),xyzzyaaaf8(xyzzyaaac8),xyzzyaaag8(xyzzyaaac8),xyzzyaaad&
&8,ei0,xyzzyaaar8,xyzzyaaae8,xyzzyaaas8,isnan,isinf,iszero)
if(isnan.or.isinf)call errstop('EMIN_MATRIX_GEN','Floating-point excep&
&tion reported when calculating energy derivatives. Don''t know what t&
&o do, so stopping. Please file a bug report.')
if(iszero)call errstop('EMIN_MATRIX_GEN','Zero wave function encounter&
&ed when calculating wave-function derivatives. Don''t know what to do&
&, so stopping. Please file a bug report.')
enddo
if(use_altsamp)then
xyzzyaaaq8=exp(dble(xyzzyaaar8)-logp_config(xyzzyaaaa8))
if(complex_wf)then
xyzzyaaaf8=xyzzyaaaq8*xyzzyaaaf8
xyzzyaaah8=xyzzyaaaq8*xyzzyaaah8
xyzzyaaag8=xyzzyaaaq8*xyzzyaaag8
xyzzyaaai8=xyzzyaaaq8*xyzzyaaai8
else
xyzzyaaaf8=xyzzyaaaq8*xyzzyaaaf8
xyzzyaaah8=xyzzyaaaq8*xyzzyaaah8
endif
endif
if(xyzzyaabq1)then
if(complex_wf)then
call dcopy(xyzzyaaaf1,xyzzyaaaf8(1),1,xyzzyaaaj8(1),1)
call daxpy(xyzzyaaaf1,xyzzyaaad8-energy0,xyzzyaaah8(1),1,xyzzyaaaj8(1)&
&,1)
call daxpy(xyzzyaaaf1,-ei0+energy0i,xyzzyaaai8(1),1,xyzzyaaaj8(1),1)
call dcopy(xyzzyaaaf1,xyzzyaaag8(1),1,xyzzyaaak8(1),1)
call daxpy(xyzzyaaaf1,xyzzyaaad8-energy0,xyzzyaaai8(1),1,xyzzyaaak8(1)&
&,1)
call daxpy(xyzzyaaaf1,ei0-energy0i,xyzzyaaah8(1),1,xyzzyaaak8(1),1)
do xyzzyaaac8=1,xyzzyaaaf1
xyzzyaaal8=xyzzyaaah8(xyzzyaaac8)
xyzzyaaam8=xyzzyaaai8(xyzzyaaac8)
xyzzyaaan8=xyzzyaaaj8(xyzzyaaac8)
xyzzyaaao8=xyzzyaaak8(xyzzyaaac8)
call daxpy(xyzzyaaaf1,xyzzyaaal8,xyzzyaaah8(1),1,xyzzyaaan1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,xyzzyaaam8,xyzzyaaai8(1),1,xyzzyaaan1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,xyzzyaaam8,xyzzyaaah8(1),1,xyzzyaaau1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,-xyzzyaaal8,xyzzyaaai8(1),1,xyzzyaaau1(1,xyzzyaa&
&ac8),1)
call daxpy(xyzzyaaaf1,xyzzyaaan8,xyzzyaaaj8(1),1,xyzzyaaap1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,xyzzyaaao8,xyzzyaaak8(1),1,xyzzyaaap1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,xyzzyaaao8,xyzzyaaaj8(1),1,xyzzyaaav1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,-xyzzyaaan8,xyzzyaaak8(1),1,xyzzyaaav1(1,xyzzyaa&
&ac8),1)
enddo
else
call dcopy(xyzzyaaaf1,xyzzyaaaf8(1),1,xyzzyaaaj8(1),1)
call daxpy(xyzzyaaaf1,xyzzyaaad8-energy0,xyzzyaaah8(1),1,xyzzyaaaj8(1)&
&,1)
do xyzzyaaac8=1,xyzzyaaaf1
call daxpy(xyzzyaaaf1,xyzzyaaah8(xyzzyaaac8),xyzzyaaah8(1),1,xyzzyaaan&
&1(1,xyzzyaaac8),1)
call daxpy(xyzzyaaaf1,xyzzyaaaj8(xyzzyaaac8),xyzzyaaaj8(1),1,xyzzyaaap&
&1(1,xyzzyaaac8),1)
enddo
endif
else
if(complex_wf)then
call dcopy(xyzzyaaaf1,xyzzyaaaf8(1),1,xyzzyaaaj8(1),1)
call daxpy(xyzzyaaaf1,xyzzyaaad8,xyzzyaaah8(1),1,xyzzyaaaj8(1),1)
call daxpy(xyzzyaaaf1,-ei0,xyzzyaaai8(1),1,xyzzyaaaj8(1),1)
call dcopy(xyzzyaaaf1,xyzzyaaag8(1),1,xyzzyaaak8(1),1)
call daxpy(xyzzyaaaf1,xyzzyaaad8,xyzzyaaai8(1),1,xyzzyaaak8(1),1)
call daxpy(xyzzyaaaf1,ei0,xyzzyaaah8(1),1,xyzzyaaak8(1),1)
do xyzzyaaac8=1,xyzzyaaaf1
xyzzyaaal8=xyzzyaaah8(xyzzyaaac8)
xyzzyaaam8=xyzzyaaai8(xyzzyaaac8)
xyzzyaaan8=xyzzyaaaj8(xyzzyaaac8)
xyzzyaaao8=xyzzyaaak8(xyzzyaaac8)
call daxpy(xyzzyaaaf1,xyzzyaaal8,xyzzyaaah8(1),1,xyzzyaaan1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,xyzzyaaam8,xyzzyaaai8(1),1,xyzzyaaan1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,xyzzyaaam8,xyzzyaaah8(1),1,xyzzyaaau1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,-xyzzyaaal8,xyzzyaaai8(1),1,xyzzyaaau1(1,xyzzyaa&
&ac8),1)
call daxpy(xyzzyaaaf1,xyzzyaaan8,xyzzyaaah8(1),1,xyzzyaaap1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,xyzzyaaao8,xyzzyaaai8(1),1,xyzzyaaap1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,xyzzyaaao8,xyzzyaaah8(1),1,xyzzyaaav1(1,xyzzyaaa&
&c8),1)
call daxpy(xyzzyaaaf1,-xyzzyaaan8,xyzzyaaai8(1),1,xyzzyaaav1(1,xyzzyaa&
&ac8),1)
enddo
else
call dcopy(xyzzyaaaf1,xyzzyaaaf8(1),1,xyzzyaaaj8(1),1)
call daxpy(xyzzyaaaf1,xyzzyaaad8,xyzzyaaah8(1),1,xyzzyaaaj8(1),1)
do xyzzyaaac8=1,xyzzyaaaf1
call daxpy(xyzzyaaaf1,xyzzyaaah8(xyzzyaaac8),xyzzyaaah8(1),1,xyzzyaaan&
&1(1,xyzzyaaac8),1)
call daxpy(xyzzyaaaf1,xyzzyaaaj8(xyzzyaaac8),xyzzyaaah8(1),1,xyzzyaaap&
&1(1,xyzzyaaac8),1)
enddo
endif
endif
if(xyzzyaaba1.and.xyzzyaaaa8==1.and.opt_info>2)then
call wout(' Step sizes for numerical derivatives:')
call labelled_list(xyzzyaaae8,xyzzyaaam1,xyzzyaaas8,'step sizes not op&
&timal')
endif
enddo
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Done.')
call xyzzyaace1(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaan1,xyzzyaaaq1)
call xyzzyaace1(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaap1,xyzzyaaaq1)
if(complex_wf)then
call xyzzyaace1(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaau1,xyzzyaaaq1)
call xyzzyaace1(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaav1,xyzzyaaaq1)
endif
xyzzyaaap8=1.d0/xyzzyaaan1(1,1)
call dscal(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaap8,xyzzyaaap1(1,1),1)
call dscal(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaap8,xyzzyaaan1(1,1),1)
if(complex_wf)then
call dscal(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaap8,xyzzyaaav1(1,1),1)
call dscal(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaap8,xyzzyaaau1(1,1),1)
endif
if(xyzzyaaba1.and.opt_info>3)then
call wout('S matrix:')
do xyzzyaaac8=1,xyzzyaaaf1
call wout(' Row '//trim(i2s(xyzzyaaac8))//': ',xyzzyaaan1(:,xyzzyaaac8&
&))
enddo
call wout
call wout('H matrix:')
do xyzzyaaac8=1,xyzzyaaaf1
call wout(' Row '//trim(i2s(xyzzyaaac8))//': ',xyzzyaaap1(:,xyzzyaaac8&
&))
enddo
call wout
endif
call timer('MATRIX_GENERATION',.false.)
end subroutine xyzzyaacc1
subroutine xyzzyaacd1(icfg,iparam,fderiv_psi,fderiv_psii,fderiv_el,fde&
&riv_eli,e0,ei0,logwfn0,param_stepsizes,flag_funny_param,isnan,isinf,i&
&szero)
implicit none
integer,intent(in) :: icfg,iparam
real(dp),intent(in) :: e0,ei0
real(dp),intent(out) :: fderiv_psi,fderiv_psii,fderiv_el,fderiv_eli
real(dp),intent(inout) :: param_stepsizes(xyzzyaaae1)
complex(dp),intent(in) :: logwfn0
logical,intent(out) :: isnan,isinf,iszero
logical,intent(inout) :: flag_funny_param(xyzzyaaae1)
integer xyzzyaaaa9
real(dp) xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,ei1,xyzzyaaae9,xyzzyaaaf9,xy&
&zzyaaag9,xyzzyaaah9,xyzzyaaai9,xyzzyaaaj9,xyzzyaaak9,xyzzyaaal9,xyzzy&
&aaam9,xyzzyaaan9,xyzzyaaao9,xyzzyaaap9,xyzzyaaaq9,xyzzyaaar9,xyzzyaaa&
&s9,xyzzyaaat9,xyzzyaaau9,xyzzyaaav9,xyzzyaaaw9,xyzzyaaax9,xyzzyaaay9,&
&xyzzyaaaz9,xyzzyaaba9,xyzzyaabb9
complex(dp) xyzzyaabc9,xyzzyaabd9
logical xyzzyaabe9,xyzzyaabf9,xyzzyaabg9,xyzzyaabh9,xyzzyaabi9,xyzzyaa&
&bj9
integer,parameter :: xyzzyaabk9=20
real(dp),parameter :: xyzzyaabl9=1.d-3,xyzzyaabm9=1.d0,xyzzyaabn9=0.1d&
&0
real(dp),parameter :: xyzzyaabo9=1.d-15
real(dp),parameter :: xyzzyaabp9=1.d-7,xyzzyaabq9=1.d-7
if(icfg==1.and..not.xyzzyaaal1(iparam))then
param_stepsizes(iparam)=0.d0
xyzzyaaac9=abs(xyzzyaaag1(iparam))
if(xyzzyaaac9<xyzzyaabo9)then
xyzzyaaab9=xyzzyaabm9
xyzzyaabb9=xyzzyaabq9
else
xyzzyaaab9=xyzzyaabl9*xyzzyaaac9
xyzzyaabb9=xyzzyaabp9*xyzzyaaac9
endif
xyzzyaabh9=.false.
xyzzyaaar9=0.d0
xyzzyaaba9=0.d0
xyzzyaaai9=0.d0
xyzzyaabg9=.false.
xyzzyaaap9=0.d0
xyzzyaaaq9=0.d0
xyzzyaaay9=0.d0
xyzzyaaaz9=0.d0
xyzzyaaah9=0.d0
xyzzyaabf9=.false.
xyzzyaaan9=0.d0
xyzzyaaao9=0.d0
xyzzyaaaw9=0.d0
xyzzyaaax9=0.d0
xyzzyaaag9=0.d0
xyzzyaabe9=.false.
xyzzyaaal9=0.d0
xyzzyaaam9=0.d0
xyzzyaaau9=0.d0
xyzzyaaav9=0.d0
xyzzyaaaf9=0.d0
xyzzyaaaj9=0.d0
xyzzyaaak9=0.d0
xyzzyaaas9=0.d0
xyzzyaaat9=0.d0
xyzzyaabi9=.false.
xyzzyaabj9=.false.
do xyzzyaaaa9=0,xyzzyaabk9
xyzzyaaae9=xyzzyaaab9+xyzzyaaag1(iparam)
call make_representable(xyzzyaaae9)
xyzzyaaab9=xyzzyaaae9-xyzzyaaag1(iparam)
if(xyzzyaaab9<epsilon(1.d0)*xyzzyaaac9)exit
call eval_energy_nderiv(xyzzyaaay1,xyzzyaaag1,iparam,xyzzyaaag1(iparam&
&),xyzzyaaae9,xyzzyaaai1,e0,logwfn0,xyzzyaaad9,xyzzyaabc9,is_first=.tr&
&ue.,get_logwfn=.true.,ei0=ei0,ei1=ei1,isnan=isnan,isinf=isinf,iszero=&
&iszero)
if(isnan.or.isinf.or.iszero)then
xyzzyaabh9=.false.
xyzzyaaar9=0.d0
xyzzyaaba9=0.d0
xyzzyaaai9=0.d0
xyzzyaabg9=.false.
xyzzyaaap9=0.d0
xyzzyaaaq9=0.d0
xyzzyaaay9=0.d0
xyzzyaaaz9=0.d0
xyzzyaaah9=0.d0
xyzzyaabf9=.false.
xyzzyaaan9=0.d0
xyzzyaaao9=0.d0
xyzzyaaaw9=0.d0
xyzzyaaax9=0.d0
xyzzyaaag9=0.d0
xyzzyaabe9=.false.
xyzzyaaal9=0.d0
xyzzyaaam9=0.d0
xyzzyaaau9=0.d0
xyzzyaaav9=0.d0
xyzzyaaaf9=0.d0
xyzzyaaaj9=0.d0
xyzzyaaak9=0.d0
xyzzyaaas9=0.d0
xyzzyaaat9=0.d0
xyzzyaabi9=.false.
xyzzyaabj9=.false.
xyzzyaaab9=xyzzyaaab9*xyzzyaabn9
cycle
endif
xyzzyaabd9=xyzzyaacm1(cmplx(e0,ei0,dp),cmplx(xyzzyaaad9,ei1,dp),xyzzya&
&aab9)
fderiv_el=dble(xyzzyaabd9)
fderiv_eli=aimag(xyzzyaabd9)
xyzzyaabd9=xyzzyaacn1(logwfn0,xyzzyaabc9,xyzzyaaab9)
fderiv_psi=dble(xyzzyaabd9)
fderiv_psii=aimag(xyzzyaabd9)
if(xyzzyaaaa9==0)then
xyzzyaaaf9=xyzzyaaab9
xyzzyaaaj9=fderiv_psi
xyzzyaaak9=fderiv_psii
xyzzyaaas9=fderiv_el
xyzzyaaat9=fderiv_eli
elseif(.not.xyzzyaabj9.and.abs(log(xyzzyaaaf9)-log(xyzzyaabb9))>abs(lo&
&g(xyzzyaaab9)-log(xyzzyaabb9)))then
xyzzyaaaf9=xyzzyaaab9
xyzzyaaaj9=fderiv_psi
xyzzyaaak9=fderiv_psii
xyzzyaaas9=fderiv_el
xyzzyaaat9=fderiv_eli
endif
xyzzyaabh9=xyzzyaabg9
xyzzyaaar9=xyzzyaaap9
xyzzyaaba9=xyzzyaaay9
xyzzyaaai9=xyzzyaaah9
xyzzyaabg9=xyzzyaabf9
xyzzyaaap9=xyzzyaaan9
xyzzyaaaq9=xyzzyaaao9
xyzzyaaay9=xyzzyaaaw9
xyzzyaaaz9=xyzzyaaax9
xyzzyaaah9=xyzzyaaag9
xyzzyaabf9=xyzzyaabe9
xyzzyaaan9=xyzzyaaal9
xyzzyaaao9=xyzzyaaam9
xyzzyaaaw9=xyzzyaaau9
xyzzyaaax9=xyzzyaaav9
xyzzyaaag9=xyzzyaaab9
xyzzyaabe9=.true.
xyzzyaaal9=fderiv_psi
xyzzyaaam9=fderiv_psii
xyzzyaaau9=fderiv_el
xyzzyaaav9=fderiv_eli
if(xyzzyaabh9)then
if(.not.xyzzyaaak1(iparam))then
if(abs(xyzzyaaar9-xyzzyaaap9)>abs(xyzzyaaap9-xyzzyaaan9).and.abs(xyzzy&
&aaap9-xyzzyaaan9)<abs(xyzzyaaan9-xyzzyaaal9))then
xyzzyaaab9=xyzzyaaai9
xyzzyaabi9=.true.
fderiv_psi=xyzzyaaap9
fderiv_psii=xyzzyaaaq9
fderiv_el=xyzzyaaay9
fderiv_eli=xyzzyaaaz9
exit
elseif(abs(xyzzyaaba9-xyzzyaaay9)>abs(xyzzyaaay9-xyzzyaaaw9).and.abs(x&
&yzzyaaay9-xyzzyaaaw9)<abs(xyzzyaaaw9-xyzzyaaau9))then
xyzzyaaaf9=xyzzyaaai9
xyzzyaabj9=.true.
xyzzyaaaj9=xyzzyaaap9
xyzzyaaak9=xyzzyaaaq9
xyzzyaaas9=xyzzyaaay9
xyzzyaaat9=xyzzyaaaz9
elseif(xyzzyaaar9==0.d0.and.xyzzyaaap9==0.d0.and.xyzzyaaan9==0.d0.and.&
&xyzzyaaal9==0.d0)then
if(.not.xyzzyaabj9.and.xyzzyaaba9==0.d0.and.xyzzyaaay9==0.d0.and.xyzzy&
&aaaw9==0.d0.and.xyzzyaaau9==0.d0)exit
endif
else
if(abs(xyzzyaaba9-xyzzyaaay9)>abs(xyzzyaaay9-xyzzyaaaw9).and.abs(xyzzy&
&aaay9-xyzzyaaaw9)<abs(xyzzyaaaw9-xyzzyaaau9))then
xyzzyaaaf9=xyzzyaaai9
xyzzyaabj9=.true.
xyzzyaaaj9=xyzzyaaap9
xyzzyaaak9=xyzzyaaaq9
xyzzyaaas9=xyzzyaaay9
xyzzyaaat9=xyzzyaaaz9
exit
elseif(xyzzyaaba9==0.d0.and.xyzzyaaay9==0.d0.and.xyzzyaaaw9==0.d0.and.&
&xyzzyaaau9==0.d0)then
exit
endif
endif
endif
xyzzyaaab9=xyzzyaaab9*xyzzyaabn9
enddo
if(.not.xyzzyaabi9)then
if(xyzzyaabj9)then
xyzzyaaab9=xyzzyaaaf9
fderiv_psi=xyzzyaaaj9
fderiv_psii=xyzzyaaak9
fderiv_el=xyzzyaaas9
fderiv_eli=xyzzyaaat9
else
flag_funny_param(iparam)=.true.
xyzzyaaab9=xyzzyaabb9
xyzzyaaae9=xyzzyaaab9+xyzzyaaag1(iparam)
call make_representable(xyzzyaaae9)
xyzzyaaab9=xyzzyaaae9-xyzzyaaag1(iparam)
call eval_energy_nderiv(xyzzyaaay1,xyzzyaaag1,iparam,xyzzyaaag1(iparam&
&),xyzzyaaae9, xyzzyaaai1,e0,logwfn0,xyzzyaaad9,xyzzyaabc9,is_first=.t&
&rue.,get_logwfn=.true.,ei1=ei1,isnan=isnan,isinf=isinf,iszero=iszero)
xyzzyaabd9=xyzzyaacm1(cmplx(e0,ei0,dp),cmplx(xyzzyaaad9,ei1,dp),xyzzya&
&aab9)
fderiv_el=dble(xyzzyaabd9)
fderiv_eli=aimag(xyzzyaabd9)
xyzzyaabd9=xyzzyaacn1(logwfn0,xyzzyaabc9,xyzzyaaab9)
fderiv_psi=dble(xyzzyaabd9)
fderiv_psii=aimag(xyzzyaabd9)
endif
endif
param_stepsizes(iparam)=xyzzyaaab9
elseif(.not.xyzzyaaal1(iparam))then
xyzzyaaab9=param_stepsizes(iparam)
xyzzyaaae9=xyzzyaaab9+xyzzyaaag1(iparam)
call make_representable(xyzzyaaae9)
xyzzyaaab9=xyzzyaaae9-xyzzyaaag1(iparam)
call eval_energy_nderiv(xyzzyaaay1,xyzzyaaag1,iparam,xyzzyaaag1(iparam&
&),xyzzyaaae9,xyzzyaaai1,e0,logwfn0,xyzzyaaad9,xyzzyaabc9,is_first=icf&
&g==2,get_logwfn=.true.,ei1=ei1,isnan=isnan,isinf=isinf,iszero=iszero)
xyzzyaabd9=xyzzyaacm1(cmplx(e0,ei0,dp),cmplx(xyzzyaaad9,ei1,dp),xyzzya&
&aab9)
fderiv_el=dble(xyzzyaabd9)
fderiv_eli=aimag(xyzzyaabd9)
xyzzyaabd9=xyzzyaacn1(logwfn0,xyzzyaabc9,xyzzyaaab9)
fderiv_psi=dble(xyzzyaabd9)
fderiv_psii=aimag(xyzzyaabd9)
else
if(icfg==1)param_stepsizes(iparam)=-1.d0
call eval_energy_aderiv(xyzzyaaay1,iparam,xyzzyaaai1,fderiv_el,fderiv_&
&eli,dlwfn=xyzzyaabd9)
fderiv_psi=dble(xyzzyaabd9)
fderiv_psii=aimag(xyzzyaabd9)
endif
end subroutine xyzzyaacd1
subroutine xyzzyaace1(n,a,a_temp)
implicit none
integer,intent(in) :: n
real(dp),intent(inout) :: a(n),a_temp(n)
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10,xyzzyaaad10
if(.not.xyzzyaaaz1)then
if(nnodes>1)then
call mpi_allreduce(a,a_temp,n,mpi_double_precision,mpi_sum,mpi_comm_wo&
&rld,ierror)
call checkmpi(ierror,'allreducing A in emin_allreduce_sum')
call dcopy(n,a_temp(1),1,a(1),1)
endif
else
if(nnodes>1)then
call mpi_reduce(a,a_temp,n,mpi_double_precision,mpi_sum,0,mpi_comm_wor&
&ld,ierror)
call checkmpi(ierror,'reducing A in emin_allreduce_sum')
if(am_master)then
call time_report('VP communication: reducing matrices',.true.)
if(xyzzyaaba1)then
call dcopy(n,a_temp(1),1,a(1),1)
do xyzzyaaab10=1,virtual_nnodes-1
write(6,*)xyzzyaaab10
read(5,*,iostat=xyzzyaaad10)(a_temp(xyzzyaaaa10),xyzzyaaaa10=1,n)
if(xyzzyaaad10/=0)call errstop('EMIN_ALLREDUCE_SUM','VP error (A_temp)&
&.')
a=a+a_temp
enddo
write(6,*)(a(xyzzyaaaa10),xyzzyaaaa10=1,n)
else
do xyzzyaaab10=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaad10)xyzzyaaac10
if(xyzzyaaad10/=0)call errstop('EMIN_ALLREDUCE_SUM','VP error (jnode).&
&')
if(xyzzyaaac10/=virtual_node)cycle
write(6,*)(a_temp(xyzzyaaaa10),xyzzyaaaa10=1,n)
enddo
read(5,*,iostat=xyzzyaaad10)(a(xyzzyaaaa10),xyzzyaaaa10=1,n)
if(xyzzyaaad10/=0)call errstop('EMIN_ALLREDUCE_SUM','VP error (A).')
endif
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(a,n,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting A in emin_allreduce_sum')
else
call time_report('VP communication: reducing matrices',.true.)
if(xyzzyaaba1)then
do xyzzyaaab10=1,virtual_nnodes-1
write(6,*)xyzzyaaab10
read(5,*,iostat=xyzzyaaad10)(a_temp(xyzzyaaaa10),xyzzyaaaa10=1,n)
if(xyzzyaaad10/=0)call errstop('EMIN_ALLREDUCE_SUM','VP error (A_temp)&
&.')
a=a+a_temp
enddo
write(6,*)(a(xyzzyaaaa10),xyzzyaaaa10=1,n)
else
do xyzzyaaab10=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaad10)xyzzyaaac10
if(xyzzyaaad10/=0)call errstop('EMIN_ALLREDUCE_SUM','VP error (jnode).&
&')
if(xyzzyaaac10/=virtual_node)cycle
write(6,*)(a(xyzzyaaaa10),xyzzyaaaa10=1,n)
enddo
read(5,*,iostat=xyzzyaaad10)(a(xyzzyaaaa10),xyzzyaaaa10=1,n)
if(xyzzyaaad10/=0)call errstop('EMIN_ALLREDUCE_SUM','VP error (A).')
endif
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
endif
end subroutine xyzzyaace1
subroutine xyzzyaacf1(eigvec,eigval,ierr)
implicit none
integer,intent(out) :: ierr
real(dp),intent(out) :: eigvec(xyzzyaaaf1),eigval
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11
real(dp) xyzzyaaad11,xyzzyaaae11,xyzzyaaaf11,xyzzyaaag11,xyzzyaaah11,x&
&yzzyaaai11,xyzzyaaaj11,xyzzyaaak11,xyzzyaaal11,xyzzyaaam11,xyzzyaaan1&
&1,xyzzyaaao11(xyzzyaaae1),xyzzyaaap11(xyzzyaaae1),xyzzyaaaq11(xyzzyaa&
&ae1),xyzzyaaar11,xyzzyaaas11,xyzzyaaat11(xyzzyaaae1),xyzzyaaau11(xyzz&
&yaaae1)
complex(dp) :: xyzzyaaav11,xyzzyaaaw11,xyzzyaaax11,xyzzyaaay11,xyzzyaa&
&az11,xyzzyaaba11,xyzzyaabb11
logical xyzzyaabc11
call timer('MATRIX_ALGEBRA',.true.)
ierr=0
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Performing matri&
&x algebra.',1,1)
if(.not.xyzzyaabq1)then
if(xyzzyaaba1.and.opt_info>2)call wout(' Performing first basis shift.&
&')
if(complex_wf)then
xyzzyaaat11(1:xyzzyaaae1)=-xyzzyaaan1(2:xyzzyaaaf1,1)
xyzzyaaau11(1:xyzzyaaae1)=-xyzzyaaau1(2:xyzzyaaaf1,1)
call xyzzyaach1(xyzzyaaan1,xyzzyaaau1,xyzzyaaat11,xyzzyaaau11)
call xyzzyaach1(xyzzyaaap1,xyzzyaaav1,xyzzyaaat11,xyzzyaaau11)
else
xyzzyaaat11(1:xyzzyaaae1)=-xyzzyaaan1(2:xyzzyaaaf1,1)
call xyzzyaacg1(xyzzyaaan1,xyzzyaaat11)
call xyzzyaacg1(xyzzyaaap1,xyzzyaaat11)
endif
endif
do xyzzyaaaa11=1,xyzzyaaae1
xyzzyaaan1(xyzzyaaaa11+1,xyzzyaaaa11+1)=max(0.d0,xyzzyaaan1(xyzzyaaaa1&
&1+1,xyzzyaaaa11+1))
enddo
if(xyzzyaaba1.and.opt_info>2)call wout(' Regularizing S and H.')
do xyzzyaaaa11=1,xyzzyaaae1
xyzzyaaaq11(xyzzyaaaa11)=sqrt(xyzzyaaan1(xyzzyaaaa11+1,xyzzyaaaa11+1))
enddo
xyzzyaaac11=0
if(all(xyzzyaaaq11==0.d0))then
xyzzyaaar11=1.d0
else
xyzzyaaar11=maxval(xyzzyaaaq11)*sqrt(epsilon(1.d0))
endif
do xyzzyaaaa11=1,xyzzyaaae1
if(xyzzyaaaq11(xyzzyaaaa11)<xyzzyaaar11)then
xyzzyaaac11=xyzzyaaac11+1
xyzzyaaah1(xyzzyaaaa11)=0.d0
else
xyzzyaaah1(xyzzyaaaa11)=1.d0/xyzzyaaaq11(xyzzyaaaa11)
endif
call dscal(xyzzyaaaf1,xyzzyaaah1(xyzzyaaaa11),xyzzyaaan1(xyzzyaaaa11+1&
&,1),xyzzyaaaf1)
call dscal(xyzzyaaaf1,xyzzyaaah1(xyzzyaaaa11),xyzzyaaan1(1,xyzzyaaaa11&
&+1),1)
call dscal(xyzzyaaaf1,xyzzyaaah1(xyzzyaaaa11),xyzzyaaap1(xyzzyaaaa11+1&
&,1),xyzzyaaaf1)
call dscal(xyzzyaaaf1,xyzzyaaah1(xyzzyaaaa11),xyzzyaaap1(1,xyzzyaaaa11&
&+1),1)
if(complex_wf)then
call dscal(xyzzyaaaf1,xyzzyaaah1(xyzzyaaaa11),xyzzyaaau1(xyzzyaaaa11+1&
&,1),xyzzyaaaf1)
call dscal(xyzzyaaaf1,xyzzyaaah1(xyzzyaaaa11),xyzzyaaau1(1,xyzzyaaaa11&
&+1),1)
call dscal(xyzzyaaaf1,xyzzyaaah1(xyzzyaaaa11),xyzzyaaav1(xyzzyaaaa11+1&
&,1),xyzzyaaaf1)
call dscal(xyzzyaaaf1,xyzzyaaah1(xyzzyaaaa11),xyzzyaaav1(1,xyzzyaaaa11&
&+1),1)
endif
enddo
if(xyzzyaaba1.and.opt_info>1)then
if(xyzzyaaac11==1)then
call wout(' Found 1 parameter with negligible derivative.')
elseif(xyzzyaaac11>1)then
call wout(' Found '//trim(i2s(xyzzyaaac11))//' parameters with negligi&
&ble derivatives.')
endif
endif
if(xyzzyaaba1.and.opt_info>2)call wout(' Building and solving eigensys&
&tem.')
xyzzyaaao1=inverse(xyzzyaaaf1,xyzzyaaaf1,xyzzyaaan1,nzeros=xyzzyaaac11&
&,lerror=xyzzyaabc11)
if(xyzzyaabc11)then
call errwarn('EMIN_MATRIX_ALGEBRA','Could not invert S (stage 1).')
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Done.')
call timer('MATRIX_ALGEBRA',.false.)
ierr=1
return
endif
if(xyzzyaaba1.and.opt_info>1)then
if(xyzzyaaac11==1)then
call wout(' Found 1 singularity inverting S (stage 1).')
elseif(xyzzyaaac11>1)then
call wout(' Found '//trim(i2s(xyzzyaaac11))//' singularities inverting&
& S (stage 1).')
endif
endif
call dmatmul(xyzzyaaao1,xyzzyaaap1,xyzzyaaaq1,xyzzyaaaf1,xyzzyaaaf1,xy&
&zzyaaaf1)
call xyzzyaack1(xyzzyaaaq1,eigvec,eigval,xyzzyaaab11)
if(xyzzyaaab11/=0)then
call errwarn('EMIN_MATRIX_ALGEBRA','Could not find new parameter vecto&
&r from E: '//trim(signal_text(xyzzyaaab11))//' (stage 1).')
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Done.')
call timer('MATRIX_ALGEBRA',.false.)
ierr=2
return
endif
if(.not.xyzzyaabq1)then
if(xyzzyaaba1.and.opt_info>2)call wout(' Performing second basis shift&
&.')
xyzzyaaao11=0.d0
xyzzyaaap11=0.d0
xyzzyaaag11=ddot(xyzzyaaae1,eigvec(2),1,xyzzyaaan1(2,1),1)
if(complex_wf)xyzzyaaah11=ddot(xyzzyaaae1,eigvec(2),1,xyzzyaaau1(2,1),&
&1)
xyzzyaaaf11=2.d0*sum(eigvec(2:xyzzyaaaf1)*xyzzyaaan1(2:xyzzyaaaf1,1),.&
&not.xyzzyaaak1(:))
do xyzzyaaaa11=1,xyzzyaaae1
if(.not.xyzzyaaak1(xyzzyaaaa11))xyzzyaaaf11=xyzzyaaaf11+eigvec(xyzzyaa&
&aa11+1)*sum(eigvec(2:xyzzyaaaf1)*xyzzyaaan1(2:xyzzyaaaf1,xyzzyaaaa11+&
&1),.not.xyzzyaaak1(:))
enddo
do xyzzyaaaa11=1,xyzzyaaae1
if(xyzzyaaak1(xyzzyaaaa11))then
xyzzyaaao11(xyzzyaaaa11)=-xyzzyaaat11(xyzzyaaaa11)*xyzzyaaah1(xyzzyaaa&
&a11)
if(complex_wf)xyzzyaaap11(xyzzyaaaa11)=-xyzzyaaau11(xyzzyaaaa11)*xyzzy&
&aaah1(xyzzyaaaa11)
else
xyzzyaaad11=sum(eigvec(2:xyzzyaaaf1)*xyzzyaaan1(xyzzyaaaa11+1,2:xyzzya&
&aaf1),.not.xyzzyaaak1(:))
if(complex_wf)xyzzyaaae11=sum(eigvec(2:xyzzyaaaf1)*xyzzyaaau1(xyzzyaaa&
&a11+1,2:xyzzyaaaf1),.not.xyzzyaaak1(:))
xyzzyaaas11=sqrt(1+xyzzyaaaf11)
if(complex_wf)then
xyzzyaaav11=xyzzyaaas11*cmplx(xyzzyaaan1(xyzzyaaaa11+1,1),xyzzyaaau1(x&
&yzzyaaaa11+1,1),dp)
xyzzyaaaw11=cmplx(xyzzyaaan1(xyzzyaaaa11+1,1)+xyzzyaaad11,xyzzyaaau1(x&
&yzzyaaaa11+1,1)+xyzzyaaae11,dp)
xyzzyaaax11=cmplx(xyzzyaaas11,0.d0,dp)
xyzzyaaay11=cmplx(1.d0+xyzzyaaag11,xyzzyaaah11,dp)
xyzzyaaaz11=emin_xi_value*xyzzyaaav11+(1.d0-emin_xi_value)*xyzzyaaaw11
xyzzyaaba11=emin_xi_value*xyzzyaaax11+(1.d0-emin_xi_value)*xyzzyaaay11
xyzzyaabb11=xyzzyaaaz11/xyzzyaaba11
xyzzyaaao11(xyzzyaaaa11)=dble(xyzzyaabb11)
xyzzyaaap11(xyzzyaaaa11)=aimag(xyzzyaabb11)
else
xyzzyaaai11=xyzzyaaas11*xyzzyaaan1(xyzzyaaaa11+1,1)
xyzzyaaaj11=xyzzyaaan1(xyzzyaaaa11+1,1)+xyzzyaaad11
xyzzyaaak11=xyzzyaaas11
xyzzyaaal11=1.d0+xyzzyaaag11
xyzzyaaam11=emin_xi_value*xyzzyaaai11+(1.d0-emin_xi_value)*xyzzyaaaj11
xyzzyaaan11=emin_xi_value*xyzzyaaak11+(1.d0-emin_xi_value)*xyzzyaaal11
xyzzyaaao11(xyzzyaaaa11)=xyzzyaaam11/xyzzyaaan11
endif
endif
enddo
if(complex_wf)then
call xyzzyaach1(xyzzyaaan1,xyzzyaaau1,xyzzyaaao11,xyzzyaaap11)
call xyzzyaach1(xyzzyaaap1,xyzzyaaav1,xyzzyaaao11,xyzzyaaap11)
else
call xyzzyaacg1(xyzzyaaan1,xyzzyaaao11)
call xyzzyaacg1(xyzzyaaap1,xyzzyaaao11)
endif
if(xyzzyaaba1.and.opt_info>2)call wout(' Building and solving eigensys&
&tem.')
xyzzyaaao1=inverse(xyzzyaaaf1,xyzzyaaaf1,xyzzyaaan1,nzeros=xyzzyaaac11&
&,lerror=xyzzyaabc11)
if(xyzzyaabc11)then
call errwarn('EMIN_MATRIX_ALGEBRA','Could not invert S (stage 2).')
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Done.')
call timer('MATRIX_ALGEBRA',.false.)
ierr=3
return
endif
if(xyzzyaaba1.and.opt_info>1)then
if(xyzzyaaac11==1)then
call wout(' Found 1 singularity inverting S (stage 2).')
elseif(xyzzyaaac11>1)then
call wout(' Found '//trim(i2s(xyzzyaaac11))//' singularities inverting&
& S (stage 2).')
endif
endif
call dmatmul(xyzzyaaao1,xyzzyaaap1,xyzzyaaaq1,xyzzyaaaf1,xyzzyaaaf1,xy&
&zzyaaaf1)
call xyzzyaack1(xyzzyaaaq1,eigvec,eigval,xyzzyaaab11)
if(xyzzyaaab11/=0)then
call errwarn('EMIN_MATRIX_ALGEBRA','Could not find new parameter vecto&
&r from E: '//trim(signal_text(xyzzyaaab11))//' (stage 2).')
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Done.')
call timer('MATRIX_ALGEBRA',.false.)
ierr=4
return
endif
endif
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Done.')
call timer('MATRIX_ALGEBRA',.false.)
end subroutine xyzzyaacf1
subroutine xyzzyaacg1(x,c)
implicit none
real(dp),intent(in) :: c(xyzzyaaae1)
real(dp),intent(inout) :: x(xyzzyaaaf1,xyzzyaaaf1)
integer xyzzyaaaa12,xyzzyaaab12
real(dp) xyzzyaaac12
xyzzyaaac12=x(1,1)
do xyzzyaaaa12=1,xyzzyaaae1
xyzzyaaab12=xyzzyaaaa12+1
call daxpy(xyzzyaaae1,x(1,xyzzyaaab12)+xyzzyaaac12*c(xyzzyaaaa12),c(1)&
&,1,x(2,xyzzyaaab12),1)
call daxpy(xyzzyaaae1,x(xyzzyaaab12,1),c(1),1,x(xyzzyaaab12,2),xyzzyaa&
&af1)
enddo
call daxpy(xyzzyaaae1,xyzzyaaac12,c(1),1,x(2,1),1)
call daxpy(xyzzyaaae1,xyzzyaaac12,c(1),1,x(1,2),xyzzyaaaf1)
end subroutine xyzzyaacg1
subroutine xyzzyaach1(x,xi,c,ci)
implicit none
real(dp),intent(in) :: c(xyzzyaaae1),ci(xyzzyaaae1)
real(dp),intent(inout) :: x(xyzzyaaaf1,xyzzyaaaf1),xi(xyzzyaaaf1,xyzzy&
&aaaf1)
integer xyzzyaaaa13,xyzzyaaab13
real(dp) xyzzyaaac13
xyzzyaaac13=x(1,1)
do xyzzyaaaa13=1,xyzzyaaae1
xyzzyaaab13=xyzzyaaaa13+1
call daxpy(xyzzyaaae1,x(1,xyzzyaaab13)+xyzzyaaac13*c(xyzzyaaaa13),c(1)&
&,1,x(2,xyzzyaaab13),1)
call daxpy(xyzzyaaae1,xi(1,xyzzyaaab13)+xyzzyaaac13*ci(xyzzyaaaa13),ci&
&(1),1,x(2,xyzzyaaab13),1)
call daxpy(xyzzyaaae1,xi(1,xyzzyaaab13)-xyzzyaaac13*ci(xyzzyaaaa13),c(&
&1),1,xi(2,xyzzyaaab13),1)
call daxpy(xyzzyaaae1,-x(1,xyzzyaaab13)-xyzzyaaac13*c(xyzzyaaaa13),ci(&
&1),1,xi(2,xyzzyaaab13),1)
call daxpy(xyzzyaaae1,x(xyzzyaaab13,1),c(1),1,x(xyzzyaaab13,2),xyzzyaa&
&af1)
call daxpy(xyzzyaaae1,-xi(xyzzyaaab13,1),ci(1),1,x(xyzzyaaab13,2),xyzz&
&yaaaf1)
call daxpy(xyzzyaaae1,xi(xyzzyaaab13,1),c(1),1,xi(xyzzyaaab13,2),xyzzy&
&aaaf1)
call daxpy(xyzzyaaae1,x(xyzzyaaab13,1),ci(1),1,xi(xyzzyaaab13,2),xyzzy&
&aaaf1)
enddo
call daxpy(xyzzyaaae1,xyzzyaaac13,c(1),1,x(2,1),1)
call daxpy(xyzzyaaae1,-xyzzyaaac13,ci(1),1,xi(2,1),1)
call daxpy(xyzzyaaae1,xyzzyaaac13,c(1),1,x(1,2),xyzzyaaaf1)
call daxpy(xyzzyaaae1,xyzzyaaac13,ci(1),1,xi(1,2),xyzzyaaaf1)
end subroutine xyzzyaach1
subroutine xyzzyaaci1(energy0,energy0i,errorbar0,errorbar0i,variance0,&
&variance0i,eigvec_best,energy_best,energyi_best,errorbar_best,errorba&
&ri_best,variance_best,variancei_best)
implicit none
real(dp),intent(in) :: energy0,energy0i,errorbar0,errorbar0i,variance0&
&,variance0i
real(dp),intent(out) :: eigvec_best(xyzzyaaaf1),energy_best,energyi_be&
&st,errorbar_best,errorbari_best,variance_best,variancei_best
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14
integer,parameter :: xyzzyaaad14=11
real(dp) xyzzyaaae14,xyzzyaaaf14,xyzzyaaag14,xyzzyaaah14,xyzzyaaai14,x&
&yzzyaaaj14,xyzzyaaak14,xyzzyaaal14,xyzzyaaam14,xyzzyaaan14,xyzzyaaao1&
&4,xyzzyaaap14,xyzzyaaaq14,xyzzyaaar14,xyzzyaaas14,xyzzyaaat14,xyzzyaa&
&au14,xyzzyaaav14,xyzzyaaaw14,xyzzyaaax14,xyzzyaaay14(xyzzyaaaf1),xyzz&
&yaaaz14,xyzzyaaba14(xyzzyaaad14),xyzzyaabb14(xyzzyaaad14),xyzzyaabc14&
&(xyzzyaaaf1,xyzzyaaad14),xyzzyaabd14(xyzzyaaad14),xyzzyaabe14(xyzzyaa&
&ad14),xyzzyaabf14(xyzzyaaad14),xyzzyaabg14(xyzzyaaad14),xyzzyaabh14(x&
&yzzyaaad14),xyzzyaabi14(xyzzyaaad14)
real(dp),parameter :: xyzzyaabj14=1.d-6,xyzzyaabk14=1.d6
logical xyzzyaabl14,xyzzyaabm14,xyzzyaabn14,xyzzyaabo14(xyzzyaaad14)
character(80) tmpr,tmpr2
integer,parameter :: xyzzyaabp14=2
integer xyzzyaabq14
call timer('MATRIX_MANIPULATION',.false.)
xyzzyaaas14=log(xyzzyaabk14)
xyzzyaaat14=log(xyzzyaabj14)
xyzzyaaau14=tan(pi*0.4d0)
xyzzyaaav14=(xyzzyaaas14+xyzzyaaat14)/2.d0
xyzzyaaaw14=(xyzzyaaas14-xyzzyaaat14)/(2.d0*xyzzyaaau14)
xyzzyaaax14=1.d0/dble(xyzzyaaad14-1)
xyzzyaaba14(1)=0.d0
xyzzyaaba14(xyzzyaaad14)=1.d0
do xyzzyaaaa14=2,xyzzyaaad14-1
xyzzyaaba14(xyzzyaaaa14)=xyzzyaaba14(xyzzyaaaa14-1)+xyzzyaaax14
enddo
xyzzyaabq14=0
do xyzzyaaaa14=1,xyzzyaaad14
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Optimizing manip&
&ulation constant.',xyzzyaaaa14,xyzzyaaad14+1)
xyzzyaaae14=xyzzyaaba14(xyzzyaaaa14)
if(xyzzyaaaa14==1)then
xyzzyaaam14=energy0
xyzzyaaan14=energy0i
xyzzyaaao14=errorbar0
xyzzyaaap14=errorbar0i
xyzzyaaaq14=variance0
xyzzyaaar14=variance0i
xyzzyaaab14=0
xyzzyaaay14=0.d0
else
call xyzzyaacj1(xyzzyaaam14,xyzzyaaan14,xyzzyaaao14,xyzzyaaap14,xyzzya&
&aaq14,xyzzyaaar14,xyzzyaaae14,xyzzyaaav14,xyzzyaaaw14,xyzzyaaab14,xyz&
&zyaaay14)
endif
call xyzzyaacp1(xyzzyaaam14,xyzzyaaan14,xyzzyaaao14,xyzzyaaap14,xyzzya&
&aaq14,xyzzyaaar14,xyzzyaaai14,target_f_print=xyzzyaaaz14)
xyzzyaabo14(xyzzyaaaa14)=xyzzyaaab14==0
xyzzyaabb14(xyzzyaaaa14)=xyzzyaaai14
xyzzyaabc14(1:xyzzyaaaf1,xyzzyaaaa14)=xyzzyaaay14(1:xyzzyaaaf1)
xyzzyaabd14(xyzzyaaaa14)=xyzzyaaam14
xyzzyaabe14(xyzzyaaaa14)=xyzzyaaan14
xyzzyaabf14(xyzzyaaaa14)=xyzzyaaao14
xyzzyaabg14(xyzzyaaaa14)=xyzzyaaap14
xyzzyaabh14(xyzzyaaaa14)=xyzzyaaaq14
xyzzyaabi14(xyzzyaaaa14)=xyzzyaaar14
if(xyzzyaaba1.and.opt_info>2)then
write(tmpr,'(f8.5)')xyzzyaaae14
write(tmpr2,*)xyzzyaaaz14
if(xyzzyaaab14==0)then
call wout(' At x = '//trim(adjustl(tmpr))//': f = '//trim(adjustl(tmpr&
&2)))
else
if(xyzzyaaab14<1.or.xyzzyaaab14>nsignal)then
tmpr2='unspecified problem'
else
tmpr2=signal_text(xyzzyaaab14)
endif
call wout(' At x = '//trim(adjustl(tmpr))//': '//trim(tmpr2))
endif
endif
if(xyzzyaabp14>0)then
if(xyzzyaaab14==0)then
xyzzyaabq14=0
else
xyzzyaabq14=xyzzyaabq14+1
if(xyzzyaabq14>=xyzzyaabp14.and.xyzzyaaaa14<xyzzyaaad14)then
xyzzyaabo14(xyzzyaaaa14+1:)=.false.
xyzzyaabb14(xyzzyaaaa14+1:)=0.d0
xyzzyaabc14(1:xyzzyaaaf1,xyzzyaaaa14+1:)=0.d0
xyzzyaabd14(xyzzyaaaa14+1:)=0.d0
xyzzyaabe14(xyzzyaaaa14+1:)=0.d0
xyzzyaabf14(xyzzyaaaa14+1:)=0.d0
xyzzyaabg14(xyzzyaaaa14+1:)=0.d0
xyzzyaabh14(xyzzyaaaa14+1:)=0.d0
xyzzyaabi14(xyzzyaaaa14+1:)=0.d0
if(xyzzyaaba1.and.opt_info>2)then
if(xyzzyaabp14==1)then
call wout(' Stopping line search after 1 bad point.')
else
call wout(' Stopping line search after '//trim(i2s(xyzzyaabp14))//' co&
&nsecutive bad points.')
endif
endif
exit
endif
endif
endif
enddo
xyzzyaabm14=count(xyzzyaabo14.and.xyzzyaabf14<2.d0*errorbar0)==0
if(.not.xyzzyaabm14)then
xyzzyaaac14=minloc(xyzzyaabb14,1,xyzzyaabo14.and.xyzzyaabf14<2.d0*erro&
&rbar0)
if(xyzzyaaac14==1)xyzzyaabm14=.true.
endif
if(.not.xyzzyaabm14)then
xyzzyaaag14=xyzzyaaba14(xyzzyaaac14)
xyzzyaaak14=xyzzyaabb14(xyzzyaaac14)
eigvec_best(1:xyzzyaaaf1)=xyzzyaabc14(1:xyzzyaaaf1,xyzzyaaac14)
energy_best=xyzzyaabd14(xyzzyaaac14)
energyi_best=xyzzyaabe14(xyzzyaaac14)
errorbar_best=xyzzyaabf14(xyzzyaaac14)
errorbari_best=xyzzyaabg14(xyzzyaaac14)
variance_best=xyzzyaabh14(xyzzyaaac14)
variancei_best=xyzzyaabi14(xyzzyaaac14)
xyzzyaabn14=.false.
if(xyzzyaaac14<=1.or.xyzzyaaac14>=xyzzyaaad14)then
xyzzyaabn14=.true.
else
if(.not.xyzzyaabo14(xyzzyaaac14-1).or..not.xyzzyaabo14(xyzzyaaac14+1))&
&xyzzyaabn14=.true.
endif
if(.not.xyzzyaabn14)then
xyzzyaaaf14=xyzzyaaba14(xyzzyaaac14-1)
xyzzyaaaj14=xyzzyaabb14(xyzzyaaac14-1)
xyzzyaaah14=xyzzyaaba14(xyzzyaaac14+1)
xyzzyaaal14=xyzzyaabb14(xyzzyaaac14+1)
if(xyzzyaaaj14<=xyzzyaaak14.or.xyzzyaaal14<=xyzzyaaak14)xyzzyaabn14=.t&
&rue.
endif
if(.not.xyzzyaabn14)then
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Optimizing manip&
&ulation constant.',xyzzyaaad14+1,xyzzyaaad14+1)
call parabolic_min(xyzzyaaaf14,xyzzyaaag14,xyzzyaaah14,xyzzyaaaj14,xyz&
&zyaaak14,xyzzyaaal14,xyzzyaaae14,xyzzyaaai14,xyzzyaabl14)
if(.not.(xyzzyaabl14.or.xyzzyaaae14>xyzzyaaah14.or.xyzzyaaae14<xyzzyaa&
&af14))then
call xyzzyaacj1(xyzzyaaam14,xyzzyaaan14,xyzzyaaao14,xyzzyaaap14,xyzzya&
&aaq14,xyzzyaaar14,xyzzyaaae14,xyzzyaaav14,xyzzyaaaw14,xyzzyaaab14,xyz&
&zyaaay14)
call xyzzyaacp1(xyzzyaaam14,xyzzyaaan14,xyzzyaaao14,xyzzyaaap14,xyzzya&
&aaq14,xyzzyaaar14,xyzzyaaai14,target_f_print=xyzzyaaaz14)
if(xyzzyaaab14==0.and.xyzzyaaai14<xyzzyaaak14.and.(xyzzyaaao14<2.d0*er&
&rorbar0))then
xyzzyaaag14=xyzzyaaae14
xyzzyaaak14=xyzzyaaai14
eigvec_best=xyzzyaaay14
energy_best=xyzzyaaam14
energyi_best=xyzzyaaan14
errorbar_best=xyzzyaaao14
errorbari_best=xyzzyaaap14
variance_best=xyzzyaaaq14
variancei_best=xyzzyaaar14
endif
if(xyzzyaaba1.and.opt_info>2)then
call wout(' Parabolic interpolation performed.')
write(tmpr,'(f8.5)')xyzzyaaae14
write(tmpr2,*)xyzzyaaaz14
if(xyzzyaaab14==0)then
call wout(' At x = '//trim(adjustl(tmpr))//': f = '//trim(adjustl(tmpr&
&2)))
else
if(xyzzyaaab14<1.or.xyzzyaaab14>nsignal)then
tmpr2='unspecified problem'
else
tmpr2=signal_text(xyzzyaaab14)
endif
call wout(' At x = '//trim(adjustl(tmpr))//': '//trim(tmpr2))
endif
endif
endif
endif
if(xyzzyaaba1)then
if(opt_info>2)then
write(tmpr,'(f8.5)')xyzzyaaag14
call wout(' Selecting x = '//trim(adjustl(tmpr)))
elseif(opt_info>1)then
call wout(' Succeeded.')
endif
endif
else
if(xyzzyaaba1)then
if(opt_info>2)then
write(tmpr,'(f8.5)')1.d0
call wout(' Reverting to x = '//trim(adjustl(tmpr)))
elseif(opt_info>1)then
call wout(' Failed.')
endif
endif
eigvec_best=0.d0
energy_best=energy0
energyi_best=energy0i
errorbar_best=errorbar0
errorbari_best=errorbar0i
variance_best=variance0
variancei_best=variance0i
endif
if(xyzzyaaba1.and.opt_info>1)call loop_time_estimate('Done.')
call timer('MATRIX_MANIPULATION',.false.)
end subroutine xyzzyaaci1
subroutine xyzzyaacj1(energy,energyi,errorbar,errorbari,variance,varia&
&ncei,x,k1,k2,isignal,eigvec)
implicit none
integer,intent(out) :: isignal
real(dp),intent(in) :: x,k1,k2
real(dp),intent(out) :: energy,energyi,errorbar,errorbari,variance,var&
&iancei,eigvec(xyzzyaaaf1)
real(dp) xyzzyaaaa15,xyzzyaaab15
isignal=0
energy=0.d0
energyi=0.d0
errorbar=-1.d0
errorbari=0.d0
variance=-1.d0
variancei=0.d0
if(x==0.d0)then
eigvec=0.d0
else
xyzzyaaaa15=xyzzyaaco1(x,k1,k2)
call dcopy(xyzzyaaaf1*xyzzyaaaf1,xyzzyaaaq1(1,1),1,xyzzyaaap1(1,1),1)
call daxpy(xyzzyaaae1*xyzzyaaaf1,xyzzyaaaa15,xyzzyaaao1(1,2),1,xyzzyaa&
&ap1(1,2),1)
call xyzzyaack1(xyzzyaaap1,eigvec,xyzzyaaab15,isignal)
if(isignal/=0)return
endif
call put_params(xyzzyaaag1(:)+eigvec(2:xyzzyaaaf1)*xyzzyaaah1(:),xyzzy&
&aaai1,.false.)
call invalidate_params(xyzzyaaai1)
call invalidate_params_energy
call xyzzyaacl1(energy,energyi,errorbar,errorbari,variance,variancei,i&
&signal)
if(isignal==0.and..not.xyzzyaabq1.and.energy<xyzzyaabb1)isignal=xyzzya&
&abj1
call put_params(xyzzyaaag1,xyzzyaaai1,.false.)
call invalidate_params(xyzzyaaai1)
call invalidate_params_energy
end subroutine xyzzyaacj1
subroutine xyzzyaack1(xyzzyaaaq1,eigvec,eigval,isignal)
implicit none
integer,intent(out) :: isignal
real(dp),intent(in) :: xyzzyaaaq1(xyzzyaaaf1,xyzzyaaaf1)
real(dp),intent(out) :: eigvec(xyzzyaaaf1),eigval
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16
real(dp) xyzzyaaad16
complex(dp) xyzzyaaae16,xyzzyaaaf16(xyzzyaaaf1)
logical xyzzyaaag16,xyzzyaaah16
isignal=0
eigvec=0.d0
eigval=0.d0
call eigenproblem(xyzzyaaaq1,xyzzyaaaf1,xyzzyaaaf16,xyzzyaaax1,xyzzyaa&
&ag16)
if(xyzzyaaag16)then
isignal=xyzzyaabg1
return
endif
xyzzyaaad16=xyzzyaabb1
if(xyzzyaabq1)xyzzyaaad16=0.d0
xyzzyaaaa16=0
do xyzzyaaac16=1,xyzzyaaaf1
xyzzyaaah16=.true.
if(xyzzyaaaa16>0)xyzzyaaah16=dble(xyzzyaaaf16(xyzzyaaac16))<dble(xyzzy&
&aaaf16(xyzzyaaaa16))
if(xyzzyaaah16)then
if(dble(xyzzyaaaf16(xyzzyaaac16))>=xyzzyaaad16)then
if(dble(xyzzyaaaf16(xyzzyaaac16))/=0.d0.or.dble(xyzzyaaax1(1,xyzzyaaac&
&16))/=0.d0)xyzzyaaaa16=xyzzyaaac16
endif
endif
enddo
if(xyzzyaaaa16==0)then
isignal=xyzzyaabh1
return
endif
eigval=dble(xyzzyaaaf16(xyzzyaaaa16))
xyzzyaaab16=maxloc(abs(xyzzyaaax1(:,xyzzyaaaa16)),1)
xyzzyaaae16=cmplx(abs(xyzzyaaax1(xyzzyaaab16,xyzzyaaaa16)),0.d0,dp)/xy&
&zzyaaax1(xyzzyaaab16,xyzzyaaaa16)
xyzzyaaax1(:,xyzzyaaaa16)=xyzzyaaax1(:,xyzzyaaaa16)*xyzzyaaae16
eigvec=real(xyzzyaaax1(:,xyzzyaaaa16),dp)
if(abs(eigvec(1))<sqrt(tiny(1.d0)))then
isignal=xyzzyaabi1
return
endif
call dscal(xyzzyaaaf1,1.d0/eigvec(1),eigvec(1),1)
end subroutine xyzzyaack1
subroutine xyzzyaacl1(energy,energyi,errorbar,errorbari,variance,varia&
&ncei,isignal)
implicit none
integer,intent(out) :: isignal
real(dp),intent(out) :: energy,energyi,errorbar,errorbari,variance,var&
&iancei
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17,xyzzyaaae17,xy&
&zzyaaaf17
real(dp) xyzzyaaag17
complex(dp) xyzzyaaah17
logical isnan,isinf,xyzzyaaai17,nan,inf,xyzzyaaaj17
energy=0.d0
energyi=0.d0
errorbar=0.d0
errorbari=0.d0
variance=0.d0
variancei=0.d0
isignal=0
xyzzyaaat1=0.d0
xyzzyaaar1=0.d0
isnan=.false.
isinf=.false.
xyzzyaaai17=.false.
do xyzzyaaaa17=1,xyzzyaaaa1
call load_from_storage(xyzzyaaay1,xyzzyaaaa17)
call load_from_storage_energy(xyzzyaaay1,xyzzyaaaa17)
if(xyzzyaabs1)call mc_twist_offset(twist_config(1:3,xyzzyaaaa17))
call eval_local_energy(xyzzyaaay1,etot=xyzzyaaar1(xyzzyaaaa17),etot_im&
&ag=xyzzyaaag17,fix_nl_grid=opt_fixnl,isnan=nan,isinf=inf)
isnan=isnan.or.nan
isinf=isinf.or.inf
if(xyzzyaabs1)xyzzyaaar1(xyzzyaaaa17)=xyzzyaaar1(xyzzyaaaa17)+xyzzyaab&
&t1*xyzzyaabv1(xyzzyaaaa17)+xyzzyaabu1*xyzzyaabw1(xyzzyaaaa17)
if(complex_wf)xyzzyaaas1(xyzzyaaaa17)=xyzzyaaag17
call wfn_logval(xyzzyaaay1,xyzzyaaah17,xyzzyaaaj17)
xyzzyaaai17=xyzzyaaai17.or.xyzzyaaaj17
call save_to_storage(xyzzyaaay1,xyzzyaaaa17)
call save_to_storage_energy(xyzzyaaay1,xyzzyaaaa17)
if(xyzzyaabr1)then
if(.not.use_altsamp)then
xyzzyaaat1(xyzzyaaaa17)=2.d0*dble(xyzzyaaah17-xyzzyaaaw1(xyzzyaaaa17))
else
xyzzyaaat1(xyzzyaaaa17)=2.d0*dble(xyzzyaaah17-cmplx(logp_config(xyzzya&
&aaa17),0.d0,dp))
endif
endif
enddo
call mpi_reduce(isnan,nan,1,mpi_logical,mpi_lor,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'reducing isNaN is emin_energy_recalc')
call mpi_reduce(isinf,inf,1,mpi_logical,mpi_lor,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'reducing isInf is emin_energy_recalc')
call mpi_reduce(xyzzyaaai17,xyzzyaaaj17,1,mpi_logical,mpi_lor,0,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'reducing isZero is emin_energy_recalc')
if(am_master)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaar1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_energies in emin_energy_recalc')
if(complex_wf)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaas1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_energies_imag in emin_energy_recalc&
&')
endif
if(xyzzyaabr1)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaat1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_weights in emin_energy_recalc')
endif
else
call mpi_gather(xyzzyaaar1,xyzzyaaaa1,mpi_double_precision,xyzzyaaar1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_energies in emin_energy_recalc')
if(complex_wf)then
call mpi_gather(xyzzyaaas1,xyzzyaaaa1,mpi_double_precision,xyzzyaaas1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_energies_imag in emin_energy_recalc&
&')
endif
if(xyzzyaabr1)then
call mpi_gather(xyzzyaaat1,xyzzyaaaa1,mpi_double_precision,xyzzyaaat1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering cs_weights in emin_energy_recalc')
endif
endif
if(xyzzyaaaz1.and.am_master)then
call time_report(' VP communication: gathering energies',.true.)
if(xyzzyaaba1)then
xyzzyaaab17=xyzzyaaab1
do xyzzyaaad17=1,virtual_nnodes-1
write(6,*)xyzzyaaad17
read(5,*,iostat=ierror)nan,inf,xyzzyaaaj17
if(ierror/=0)call errstop('EMIN_ENERGY_RECALC','VP error (NaN).')
if(nan)isnan=.true.
if(inf)isinf=.true.
if(xyzzyaaaj17)xyzzyaaai17=.true.
read(5,*,iostat=ierror)xyzzyaaac17
if(ierror/=0)call errstop('EMIN_ENERGY_RECALC','VP error (kcfg).')
read(5,*,iostat=ierror)(xyzzyaaar1(xyzzyaaaa17),xyzzyaaaa17=xyzzyaaab1&
&7+1,xyzzyaaab17+xyzzyaaac17)
if(ierror/=0)call errstop('EMIN_ENERGY_RECALC','VP error (cs_energies)&
&.')
if(complex_wf)then
read(5,*,iostat=ierror)(xyzzyaaas1(xyzzyaaaa17),xyzzyaaaa17=xyzzyaaab1&
&7+1,xyzzyaaab17+xyzzyaaac17)
if(ierror/=0)call errstop('EMIN_ENERGY_RECALC','VP error (cs_energies_&
&imag).')
endif
if(xyzzyaabr1)then
read(5,*,iostat=ierror)(xyzzyaaat1(xyzzyaaaa17),xyzzyaaaa17=xyzzyaaab1&
&7+1,xyzzyaaab17+xyzzyaaac17)
if(ierror/=0)call errstop('EMIN_ENERGY_RECALC','VP error (cs_weights).&
&')
endif
xyzzyaaab17=xyzzyaaab17+xyzzyaaac17
enddo
if(xyzzyaaab17/=virtual_nconfig)call errstop('EMIN_ENERGY_RECALC','VP &
&error (config count '//trim(i2s(xyzzyaaab17))//'/='//trim(i2s(virtual&
&_nconfig))//').')
else
do xyzzyaaad17=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaaf17)xyzzyaaae17
if(xyzzyaaaf17/=0)call errstop('EMIN_ENERGY_RECALC','VP error (jnode).&
&')
if(xyzzyaaae17/=virtual_node)cycle
write(6,*)isnan,isinf,xyzzyaaai17
write(6,*)xyzzyaaab1
write(6,*)(xyzzyaaar1(xyzzyaaaa17),xyzzyaaaa17=1,xyzzyaaab1)
if(complex_wf)write(6,*)(xyzzyaaas1(xyzzyaaaa17),xyzzyaaaa17=1,xyzzyaa&
&ab1)
if(xyzzyaabr1)write(6,*)(xyzzyaaat1(xyzzyaaaa17),xyzzyaaaa17=1,xyzzyaa&
&ab1)
enddo
endif
call time_report(' Done.',report_walltime=.true.,no_newline=.true.)
endif
if(xyzzyaaba1)then
if(nan.and.isignal==0)isignal=xyzzyaabn1
if(inf.and.isignal==0)isignal=xyzzyaabo1
if(xyzzyaaaj17.and.isignal==0)isignal=xyzzyaabp1
if(xyzzyaabr1.and.isignal==0)then
xyzzyaaag17=maxval(xyzzyaaat1)
xyzzyaaat1=xyzzyaaat1-xyzzyaaag17
do xyzzyaaaa17=1,xyzzyaaac1
if(xyzzyaaat1(xyzzyaaaa17)>maximum_exp_arg*0.8d0)then
xyzzyaaat1(xyzzyaaaa17)=0.d0
else
xyzzyaaat1(xyzzyaaaa17)=exp(xyzzyaaat1(xyzzyaaaa17))
endif
enddo
xyzzyaaag17=sum(xyzzyaaat1)
if(xyzzyaaag17==0.d0)then
isignal=xyzzyaabm1
else
call dscal(xyzzyaaac1,xyzzyaaad1/xyzzyaaag17,xyzzyaaat1(1),1)
endif
if(isignal==0)then
if(maxval(xyzzyaaat1(:))>xyzzyaabc1)isignal=xyzzyaabk1
endif
if(isignal==0)then
xyzzyaaag17=sum((xyzzyaaat1(:)-1.d0)**2)/(xyzzyaaad1*(xyzzyaaad1-1.d0)&
&)
if(xyzzyaaag17>xyzzyaabd1)isignal=xyzzyaabl1
endif
if(isignal==0)then
energy=sum(xyzzyaaat1(:)*xyzzyaaar1(:))/xyzzyaaad1
variance=sum((xyzzyaaat1(:)*(xyzzyaaar1(:)-energy))**2)/(xyzzyaaad1-1.&
&d0)
errorbar=sqrt(variance/xyzzyaaad1)
if(complex_wf)then
energyi=sum(xyzzyaaat1(:)*xyzzyaaas1(:))/xyzzyaaad1
variancei=sum((xyzzyaaat1(:)*(xyzzyaaas1(:)-energyi))**2)/(xyzzyaaad1-&
&1.d0)
errorbari=sqrt(variancei/xyzzyaaad1)
endif
endif
else
if(isignal==0)then
energy=sum(xyzzyaaar1(:))/xyzzyaaad1
variance=sum((xyzzyaaar1(:)-energy)**2)/(xyzzyaaad1-1.d0)
errorbar=sqrt(variance/xyzzyaaad1)
if(complex_wf)then
energyi=sum(xyzzyaaas1(:))/xyzzyaaad1
variance=sum((xyzzyaaas1(:)-energyi)**2)/(xyzzyaaad1-1.d0)
errorbari=sqrt(variancei/xyzzyaaad1)
endif
endif
endif
endif
if(xyzzyaaaz1.and.am_master)then
call time_report(' VP communication: broadcasting isignal',.true.)
if(xyzzyaaba1)then
write(6,*)isignal
else
read(5,*,iostat=xyzzyaaaf17)isignal
if(xyzzyaaaf17/=0)call errstop('EMIN_ENERGY_RECALC','VP error (isignal&
&).')
endif
call time_report(' Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(isignal,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting isignal in emin_energy_recalc.')
if(isignal==0)then
if(xyzzyaaaz1.and.am_master)then
call time_report(' VP communication: broadcasting energy',.true.)
if(xyzzyaaba1)then
write(6,*)energy,errorbar,variance
if(complex_wf)write(6,*)energyi,errorbari,variancei
else
read(5,*,iostat=xyzzyaaaf17)energy,errorbar,variance
if(xyzzyaaaf17/=0)call errstop('EMIN_ENERGY_RECALC','VP error (energy)&
&.')
if(complex_wf)then
read(5,*,iostat=xyzzyaaaf17)energyi,errorbari,variancei
if(xyzzyaaaf17/=0)call errstop('EMIN_ENERGY_RECALC','VP error (energyi&
&).')
endif
endif
call time_report(' Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(energy,1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting energy in emin_energy_recalc.')
call mpi_bcast(errorbar,1,mpi_double_precision,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'broadcasting errorbar in emin_energy_recalc.')
call mpi_bcast(variance,1,mpi_double_precision,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'broadcasting variance in emin_energy_recalc.')
if(complex_wf)then
call mpi_bcast(energyi,1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting energyi in emin_energy_recalc.')
call mpi_bcast(errorbari,1,mpi_double_precision,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'broadcasting errorbari in emin_energy_recalc.')
call mpi_bcast(variancei,1,mpi_double_precision,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'broadcasting variancei in emin_energy_recalc.')
endif
endif
end subroutine xyzzyaacl1
complex(dp) function xyzzyaacm1(f0,f1,step_size)
implicit none
real(dp),intent(in) :: step_size
complex(dp),intent(in) :: f0,f1
real(dp) xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,x&
&yzzyaaaf18
complex(dp) xyzzyaaag18
xyzzyaaag18=f1-f0
xyzzyaaae18=dble(xyzzyaaag18)
xyzzyaaaf18=aimag(xyzzyaaag18)
xyzzyaaaa18=dble(f0)
xyzzyaaac18=dble(f1)
if(abs(xyzzyaaae18)<=max(abs(xyzzyaaaa18),abs(xyzzyaaac18))*epsilon(1.&
&d0)*1.d2)xyzzyaaae18=0.d0
xyzzyaaab18=aimag(f0)
xyzzyaaad18=aimag(f1)
if(abs(xyzzyaaaf18)<=max(abs(xyzzyaaab18),abs(xyzzyaaad18))*epsilon(1.&
&d0)*1.d2)xyzzyaaaf18=0.d0
xyzzyaaag18=cmplx(xyzzyaaae18,xyzzyaaaf18,dp)
xyzzyaacm1=xyzzyaaag18/step_size
end function xyzzyaacm1
complex(dp) function xyzzyaacn1(logf0,logf1,step_size)
implicit none
real(dp),intent(in) :: step_size
complex(dp),intent(in) :: logf0,logf1
real(dp) xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,x&
&yzzyaaaf19
real(dp),parameter :: xyzzyaaag19=0.5d0
complex(dp) xyzzyaaah19,xyzzyaaai19
xyzzyaaah19=logf1-logf0
xyzzyaaae19=dble(xyzzyaaah19)
xyzzyaaaf19=aimag(xyzzyaaah19)
xyzzyaaaa19=dble(logf0)
xyzzyaaac19=dble(logf1)
if(abs(xyzzyaaae19)<=max(abs(xyzzyaaaa19),abs(xyzzyaaac19))*epsilon(1.&
&d0)*1.d2)xyzzyaaae19=0.d0
xyzzyaaab19=aimag(logf0)
xyzzyaaad19=aimag(logf1)
if(abs(xyzzyaaaf19)<=max(abs(xyzzyaaab19),abs(xyzzyaaad19))*epsilon(1.&
&d0)*1.d2)xyzzyaaaf19=0.d0
xyzzyaaah19=cmplx(xyzzyaaae19,xyzzyaaaf19,dp)
if(xyzzyaaae19<minimum_exp_arg)then
xyzzyaacn1=cmplx(-1.d0/step_size,0.d0,dp)
elseif(xyzzyaaae19==0.d0.and.xyzzyaaaf19==0.d0)then
xyzzyaacn1=cmplx(0.d0,0.d0,dp)
elseif(xyzzyaaae19<maximum_exp_arg*xyzzyaaag19)then
xyzzyaaai19=exp(xyzzyaaah19)
xyzzyaacn1=(xyzzyaaai19-1.d0)/step_size
else
xyzzyaacn1=cmplx(0.d0,0.d0,dp)
endif
end function xyzzyaacn1
real(dp) function xyzzyaaco1(x,k1,k2)
implicit none
real(dp),intent(in) :: x,k1,k2
if(x<0.d0.or.x>=1.d0)then
xyzzyaaco1=-1.d0
return
endif
if(x==1.d0)then
xyzzyaaco1=0.d0
else
xyzzyaaco1=exp(k1+k2*tan(pi*(.5d0-x)))
endif
end function xyzzyaaco1
subroutine xyzzyaacp1(energy,energyi,errorbar,errorbari,variance,varia&
&ncei,target_f,criterion,target_f_print)
implicit none
real(dp),intent(in) :: energy,energyi,errorbar,errorbari,variance,vari&
&ancei
real(dp),intent(out) :: target_f
real(dp),intent(out),optional :: criterion,target_f_print
real(dp) xyzzyaaaa21
if(xyzzyaabq1)then
target_f=variance+variancei
if(present(criterion))criterion=0.01d0*(variance+variancei)
if(present(target_f_print))target_f_print=target_f
else
xyzzyaaaa21=errorbar
target_f=energy+3.d0*xyzzyaaaa21
if(present(criterion))criterion=0.01d0*xyzzyaaaa21
if(present(target_f_print))target_f_print=target_f*xyzzyaabe1
endif
end subroutine xyzzyaacp1
end module slaarnaal

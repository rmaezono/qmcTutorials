module slaarnacm
use slaarnaaa
use slaarnaaf
use slaarnaag
use dsp
use slaarnaam
use slaarnabg
use parallel
use slaarnach
use store
use slaarnacs
use casl,          only : write_casl
use slaarnaah,   only : write_correlation_header
use slaarnaas,     only : mc_twist_offset,ppmcta_partial_sum,ppmcta_fi&
&t
use format_utils,  only : wout,i2s,r2s,byte2human,labelled_list
use slaarnabt,     only : ddot,dscal,median
use slaarnacc,only : put_random_state,get_random_state
use run_control,   only : errstop,errwarn,timer,exceeds_time_limit,loo&
&p_time_estimate,time_report,check_alloc,errstop_master
implicit none
private
public varmin_main,compute_lsf_master
public vm_reweight,vm_forgiving,vm_w_max,vm_w_min,vm_use_e_guess,vm_e_&
&guess,vm_filter,vm_filter_thres,vm_filter_width,vm_madmin
public e0_prop,ew_prop
real(dp) vm_w_max,vm_w_min,vm_e_guess,vm_filter_thres,vm_filter_width
logical vm_reweight,vm_forgiving,vm_use_e_guess,vm_madmin,vm_filter
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1
real(dp) xyzzyaaad1
integer xyzzyaaae1
real(dp),allocatable :: xyzzyaaaf1(:),xyzzyaaag1(:)
logical,allocatable :: xyzzyaaah1(:)
character(2),allocatable :: xyzzyaaai1(:)
integer,allocatable :: xyzzyaaaj1(:)
real(dp),allocatable :: xyzzyaaak1(:),xyzzyaaal1(:),xyzzyaaam1(:),xyzz&
&yaaan1(:,:),xyzzyaaao1(:,:),xyzzyaaap1(:,:),xyzzyaaaq1(:)
complex(dp),allocatable :: xyzzyaaar1(:)
integer xyzzyaaas1
integer no_fn_evals,xyzzyaaat1,xyzzyaaau1
real(dp) xyzzyaaav1,e0_prop,ew_prop,xyzzyaaaw1
real(dp) xyzzyaaax1
real(dp),parameter :: xyzzyaaay1=1.d-300
logical xyzzyaaaz1,xyzzyaaba1,xyzzyaabb1,xyzzyaabc1,xyzzyaabd1
character(80) tmpr,tmpr2,tmpr3,tmpr4
real(dp),allocatable :: xyzzyaabe1(:),xyzzyaabf1(:)
contains
subroutine varmin_main(abort)
implicit none
logical,intent(out) :: abort
call timer('VARMIN',.true.)
call xyzzyaabg1
call xyzzyaabi1(abort)
call xyzzyaabh1
call timer('VARMIN',.false.)
end subroutine varmin_main
subroutine xyzzyaabg1
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3,xyzzyaaal3&
&,xyzzyaaam3,xyzzyaaan3,xyzzyaaao3,xyzzyaaap3,xyzzyaaaq3,xyzzyaaar3,xy&
&zzyaaas3
real(dp) xyzzyaaat3
logical xyzzyaaau3(1),xyzzyaaav3(1)
logical,allocatable :: xyzzyaaaw3(:),xyzzyaaax3(:),xyzzyaaay3(:),xyzzy&
&aaaz3(:),xyzzyaaba3(:)
real(dp),allocatable :: xyzzyaabb3(:)
character(2),allocatable :: xyzzyaabc3(:)
character(20),allocatable :: xyzzyaabd3(:)
character(20) opt_item(1),req_extra(0),opt_extra(1)
call timer('SETUP',.true.,collapse=.true.)
if(am_master)then
if(vm_madmin)then
call wout('MAD minimization configuration')
call wout('==============================')
else
call wout('Variance minimization configuration')
call wout('===================================')
endif
call wout()
endif
xyzzyaaaa1=-1
xyzzyaaaa3=3
if(noncoll_spin)xyzzyaaaa3=xyzzyaaaa3+1
if(use_altsamp)xyzzyaaaa3=xyzzyaaaa3+1
allocate(xyzzyaabd3(xyzzyaaaa3),stat=xyzzyaaac3)
xyzzyaaab3=3
xyzzyaabd3(1:3)=(/'RELE ','ETOT ','NLTOT'/)
if(noncoll_spin)then
xyzzyaaab3=xyzzyaaab3+1
xyzzyaabd3(xyzzyaaab3)='SELE'
endif
if(use_altsamp)then
xyzzyaaab3=xyzzyaaab3+1
xyzzyaabd3(xyzzyaaab3)='LOGP'
endif
opt_item(1)='TWIST'
opt_extra(1)='RANDOM'
call load_configs(xyzzyaaaa1,'VMC',xyzzyaabd3,opt_item,xyzzyaaau3,req_&
&extra,opt_extra,xyzzyaaav3)
deallocate(xyzzyaabd3)
xyzzyaabd1=xyzzyaaau3(1)
if(xyzzyaaav3(1))call put_random_state(random_state_config)
xyzzyaaba1=virtual_nnodes>1
xyzzyaaab1=xyzzyaaaa1*nnodes
xyzzyaaac1=xyzzyaaab1
xyzzyaabb1=am_master.and.(.not.xyzzyaaba1.or.virtual_node==0)
if(xyzzyaaba1)then
if(xyzzyaabb1)xyzzyaaac1=virtual_nconfig
if(virtual_nconfig<xyzzyaaaa1)call errstop_master('SETUP_VARMIN','VIRT&
&UAL_NCONFIG < number of configs actually present.')
endif
if(xyzzyaaac1==0)call errstop_master('SETUP_VARMIN','No configurations&
&?')
xyzzyaaad1=1.d0/dble(xyzzyaaac1)
if(isperiodic.and.model_system)then
xyzzyaaax1=1.d0/dble(netot)
else
xyzzyaaax1=1.d0/dble(npcells)
endif
if(vm_use_e_guess)xyzzyaaav1=vm_e_guess/xyzzyaaax1
if(.not.use_altsamp)then
xyzzyaaaz1=opt_info>=4.or.vm_reweight
xyzzyaabc1=vm_reweight
else
xyzzyaaaz1=.true.
xyzzyaabc1=.true.
endif
if(am_master)then
if(.not.xyzzyaaba1)then
call wout('Number of nodes                                       : '//&
&trim(i2s(nnodes)))
if(nnodes>1)call wout('Number of configurations per node              &
&       : '//trim(i2s(xyzzyaaaa1)))
call wout('Total number of configurations                        : '//&
&trim(i2s(xyzzyaaac1)))
else
if(xyzzyaabb1)then
call wout('VP optimization: this is the VP master node.')
else
call wout('VP optimization: this is VP slave node #'//trim(i2s(virtual&
&_node))//'.')
endif
if(nnodes>1)then
call wout('Number of MPI nodes on this VP node                   : '//&
&trim(i2s(nnodes)))
call wout('Number of configurations per MPI node on this VP node : '//&
&trim(i2s(xyzzyaaaa1)))
endif
call wout('Number of configurations on this VP node              : '//&
&trim(i2s(xyzzyaaab1)))
call wout('Number of VP nodes                                    : '//&
&trim(i2s(virtual_nnodes)))
call wout('Total number of configurations                        : '//&
&trim(i2s(xyzzyaaac1)))
endif
call wout()
endif
call setup_wfn_params(xyzzyaaas3)
allocate(xyzzyaabb3(xyzzyaaas3),xyzzyaaah1(xyzzyaaas3),       xyzzyaaa&
&w3(xyzzyaaas3),xyzzyaaax3(xyzzyaaas3),   xyzzyaaay3(xyzzyaaas3),xyzzy&
&aaaz3(xyzzyaaas3),    xyzzyaaba3(xyzzyaaas3),xyzzyaabc3(xyzzyaaas3),s&
&tat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','temp_params,...')
call get_params(xyzzyaabb3,xyzzyaaaw3,xyzzyaaax3,xyzzyaaay3,xyzzyaaaz3&
&,xyzzyaaba3,xyzzyaabc3)
xyzzyaaah1(:)=.false.
if(fix_cutoffs)xyzzyaaah1(:)=xyzzyaaah1(:).or.xyzzyaaaw3(:)
xyzzyaaae1=count(.not.xyzzyaaah1)
if(xyzzyaaae1<1)call errstop_master('SETUP_VARMIN','No parameters to o&
&ptimize!')
allocate(xyzzyaaaf1(xyzzyaaae1),xyzzyaaai1(xyzzyaaae1),xyzzyaaag1(xyzz&
&yaaae1),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_params,...')
xyzzyaaaq3=0
do xyzzyaaar3=1,xyzzyaaas3
if(xyzzyaaah1(xyzzyaaar3))cycle
xyzzyaaaq3=xyzzyaaaq3+1
xyzzyaaaf1(xyzzyaaaq3)=xyzzyaabb3(xyzzyaaar3)
xyzzyaaai1(xyzzyaaaq3)=xyzzyaabc3(xyzzyaaar3)
enddo
deallocate(xyzzyaabb3,xyzzyaaaw3,xyzzyaaax3,xyzzyaaay3,xyzzyaaaz3,xyzz&
&yaaba3,xyzzyaabc3)
xyzzyaaas1=0
call scratch_protect(xyzzyaaas1)
call energy_scratch_request(xyzzyaaas1)
if(xyzzyaaaz1)call scratch_request(wfn_detail=xyzzyaaas1)
call ederiv_scratch_request(xyzzyaaas1,xyzzyaaah1,wfn_deriv=xyzzyaabc1&
&)
call setup_scratch
call which_scratch(xyzzyaaas1)
call setup_wfn_utils
if(use_altsamp)call setup_alt_utils
call setup_energy_utils
call setup_storage(xyzzyaaaa1,xyzzyaaah1)
call setup_storage_energy(xyzzyaaaa1)
if(xyzzyaabb1)then
allocate(xyzzyaaak1(xyzzyaaac1),xyzzyaaam1(xyzzyaaac1),xyzzyaaan1(xyzz&
&yaaac1,xyzzyaaae1),xyzzyaaap1(xyzzyaaac1,xyzzyaaae1),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_ecfg,...','Reduce no o&
&f parameters/configs.')
xyzzyaaak1=0.d0
xyzzyaaan1=0.d0
xyzzyaaam1=0.d0
xyzzyaaap1=0.d0
if(complex_wf)then
allocate(xyzzyaaal1(xyzzyaaac1),xyzzyaaao1(xyzzyaaac1,xyzzyaaae1),stat&
&=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_eicfg,...')
endif
if(xyzzyaabd1)then
allocate(xyzzyaabe1(xyzzyaaac1),xyzzyaabf1(xyzzyaaac1),stat=xyzzyaaac3&
&)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_hfke,...')
endif
elseif(am_master)then
allocate(xyzzyaaak1(xyzzyaaab1),xyzzyaaam1(xyzzyaaab1),xyzzyaaan1(xyzz&
&yaaab1,xyzzyaaae1),xyzzyaaap1(xyzzyaaab1,xyzzyaaae1),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_ecfg,...','Reduce no o&
&f parameters/configs.')
xyzzyaaak1=0.d0
xyzzyaaan1=0.d0
xyzzyaaam1=0.d0
xyzzyaaap1=0.d0
if(complex_wf)then
allocate(xyzzyaaal1(xyzzyaaab1),xyzzyaaao1(xyzzyaaab1,xyzzyaaae1),stat&
&=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_eicfg,...')
endif
if(xyzzyaabd1)then
allocate(xyzzyaabe1(xyzzyaaab1),xyzzyaabf1(xyzzyaaab1),stat=xyzzyaaac3&
&)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_hfke,...')
endif
else
allocate(xyzzyaaak1(xyzzyaaaa1),xyzzyaaam1(xyzzyaaaa1),xyzzyaaan1(xyzz&
&yaaaa1,xyzzyaaae1),xyzzyaaap1(xyzzyaaaa1,xyzzyaaae1),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_ecfg,...','Reduce no o&
&f parameters/configs.')
xyzzyaaak1=0.d0
xyzzyaaan1=0.d0
xyzzyaaam1=0.d0
xyzzyaaap1=0.d0
if(complex_wf)then
allocate(xyzzyaaal1(xyzzyaaaa1),xyzzyaaao1(xyzzyaaaa1,xyzzyaaae1),stat&
&=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_eicfg,...')
xyzzyaaal1=0.d0
endif
if(xyzzyaabd1)then
allocate(xyzzyaabe1(xyzzyaaaa1),xyzzyaabf1(xyzzyaaaa1),stat=xyzzyaaac3&
&)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','opt_hfke,...')
endif
endif
allocate(xyzzyaaar1(xyzzyaaaa1),stat=xyzzyaaac3)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','initial_logwfn')
xyzzyaaar1=czero
if(am_master)then
xyzzyaaat3=93.d0+dble(xyzzyaaac1)*dble(xyzzyaaae1+3)+dble(xyzzyaaae1)*&
&dble(3*xyzzyaaae1+33)*0.5d0
if(xyzzyaaat3>=dble(huge(1)))call errstop('SETUP_VARMIN','Allocation p&
&roblem: size of nl2sol arrays larger than largest representable integ&
&er. Reduce number of configurations or parameters (size of array goes&
& roughly as num_configs*num_params).')
xyzzyaaai3=93+xyzzyaaac1*(xyzzyaaae1+3)+xyzzyaaae1*(3*xyzzyaaae1+33)/2
xyzzyaaaj3=60+xyzzyaaae1
allocate(xyzzyaaaq1(xyzzyaaai3),xyzzyaaaj1(xyzzyaaaj3),stat=xyzzyaaac3&
&)
call check_alloc(xyzzyaaac3,'SETUP_VARMIN','NL2SOL arrays','Reduce no.&
& of configs or parameters.')
endif
if(xyzzyaabd1)then
xyzzyaabe1(1:xyzzyaaaa1)=twist_config(4,1:xyzzyaaaa1)*dble(netot)
xyzzyaabf1(1:xyzzyaaaa1)=twist_config(5,1:xyzzyaaaa1)*dble(netot)
if(nnodes>1)then
if(am_master)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaabe1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
else
call mpi_gather(xyzzyaabe1,xyzzyaaaa1,mpi_double_precision,xyzzyaabe1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
endif
call checkmpi(ierror,'gathering opt_hfke in setup_varmin')
if(am_master)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaabf1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
else
call mpi_gather(xyzzyaabf1,xyzzyaaaa1,mpi_double_precision,xyzzyaabf1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
endif
call checkmpi(ierror,'gathering opt_hfex in setup_varmin')
endif
if(xyzzyaaba1.and.am_master)then
if(xyzzyaabb1)then
xyzzyaaal3=xyzzyaaab1
do xyzzyaaan3=1,virtual_nnodes-1
write(6,*)xyzzyaaan3
read(5,*,iostat=ierror)xyzzyaaam3
if(ierror/=0)call errstop('SETUP_VARMIN','VP error (kcfg).')
read(5,*,iostat=ierror)(xyzzyaabe1(xyzzyaaak3),xyzzyaaak3=xyzzyaaal3+1&
&,xyzzyaaal3+xyzzyaaam3)
if(ierror/=0)call errstop('SETUP_VARMIN','VP error (opt_hfke).')
read(5,*,iostat=ierror)(xyzzyaabf1(xyzzyaaak3),xyzzyaaak3=xyzzyaaal3+1&
&,xyzzyaaal3+xyzzyaaam3)
if(ierror/=0)call errstop('SETUP_VARMIN','VP error (opt_hfex).')
xyzzyaaal3=xyzzyaaal3+xyzzyaaam3
enddo
if(xyzzyaaal3/=virtual_nconfig)call errstop('SETUP_VARMIN','VP error (&
&config count).')
else
do xyzzyaaan3=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaap3)xyzzyaaao3
if(xyzzyaaap3/=0)call errstop('SETUP_VARMIN','VP error (jnode).')
if(xyzzyaaao3/=virtual_node)cycle
write(6,*)xyzzyaaab1
write(6,*)(xyzzyaabe1(xyzzyaaak3),xyzzyaaak3=1,xyzzyaaab1)
write(6,*)(xyzzyaabf1(xyzzyaaak3),xyzzyaaak3=1,xyzzyaaab1)
enddo
endif
endif
endif
if(am_master)then
xyzzyaaag3=size(xyzzyaaak1)+size(xyzzyaaan1)
if(complex_wf)xyzzyaaag3=xyzzyaaag3+size(xyzzyaaal1)+size(xyzzyaaao1)
xyzzyaaaf3=size(etot_config)
if(xyzzyaabd1)xyzzyaaaf3=xyzzyaaaf3+size(xyzzyaabe1)+size(xyzzyaabf1)
xyzzyaaah3=0
if(xyzzyaaaz1)xyzzyaaah3=size(xyzzyaaam1)+size(xyzzyaaap1)+2*size(xyzz&
&yaaar1)
xyzzyaaad3=size(rele_config)
xyzzyaaae3=0
if(noncoll_spin)xyzzyaaae3=size(sele_config)
call wout('Optimization workspace:')
call wout(repeat('-',47))
call wout('Number of variable parameters  : ',trim(i2s(xyzzyaaae1)),fm&
&t='(1x,a,a15)')
call wout('Number of configurations       : ',trim(i2s(xyzzyaaaa1)),fm&
&t='(1x,a,a15)')
call wout(repeat('-',47))
call wout('Electron positions             : ',trim(byte2human(8.d0*dbl&
&e(xyzzyaaad3))),fmt='(1x,a,a15)')
if(noncoll_spin)call wout('Electron spins                 : ',trim(byt&
&e2human(4.d0*dble(xyzzyaaae3))),fmt='(1x,a,a15)')
call wout('Other energy buffers           : ',trim(byte2human(8.d0*dbl&
&e(xyzzyaaaf3))),fmt='(1x,a,a15)')
call wout('Local energies                 : ',trim(byte2human(8.d0*dbl&
&e(xyzzyaaag3))),fmt='(1x,a,a15)')
if(xyzzyaaaz1)call wout('Weights                        : ',trim(byte2&
&human(8.d0*dble(xyzzyaaah3))),fmt='(1x,a,a15)')
call wout('NL2SOL work array (real)       : ',trim(byte2human(8.d0*dbl&
&e(xyzzyaaai3))),fmt='(1x,a,a15)')
call wout('NL2SOL work array (int)        : ',trim(byte2human(4.d0*dbl&
&e(xyzzyaaaj3))),fmt='(1x,a,a15)')
call wout(repeat('-',47))
call wout('Total memory required          : ',trim(byte2human(8.d0*(db&
&le(xyzzyaaad3)+dble(xyzzyaaaf3)+dble(xyzzyaaag3)+dble(xyzzyaaah3)+dbl&
&e(xyzzyaaai3))+4.d0*(dble(xyzzyaaae3)+dble(xyzzyaaaj3)))),fmt='(1x,a,&
&a15)')
call wout(repeat('-',47))
call wout()
endif
call timer('SETUP',.false.)
end subroutine xyzzyaabg1
subroutine xyzzyaabh1
implicit none
logical :: xyzzyaaaa4=.true.
character(20) dum_config_item(0),extra_item(1)
if(.not.rng_restart_safe)xyzzyaaaa4=.false.
deallocate(xyzzyaaak1,xyzzyaaam1,xyzzyaaan1,xyzzyaaap1)
if(complex_wf)deallocate(xyzzyaaal1,xyzzyaaao1)
if(xyzzyaabd1)deallocate(xyzzyaabe1,xyzzyaabf1)
deallocate(xyzzyaaaf1,xyzzyaaag1,xyzzyaaah1,xyzzyaaai1,xyzzyaaar1)
if(am_master)deallocate(xyzzyaaaq1,xyzzyaaaj1)
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
call finish_energy_utils
call finish_wfn_utils
if(use_altsamp)call finish_alt_utils
call finish_scratch
call finish_wfn_params
end subroutine xyzzyaabh1
subroutine xyzzyaabi1(abort)
use machine_constants
use toms573
implicit none
logical,intent(out) :: abort
character(512) errmsg
integer xyzzyaaaa5,xyzzyaaab5
real(dp) xyzzyaaac5
logical xyzzyaaad5,xyzzyaaae5
interface
subroutine madr(n,p,x,nf,r,nl2sol_iteration,been_rst,stop_opt,x_h,bad_&
&point)
use machine_constants
integer,intent(in) :: n,p,nl2sol_iteration,been_rst
integer,intent(inout) :: nf
real(kind=kind(0.d0)),intent(in) :: x(:)
real(kind=kind(0.d0)),intent(in),optional :: x_h(:)
real(kind=kind(0.d0)),intent(out) :: r(:)
logical,intent(inout) :: stop_opt
logical,intent(inout),optional :: bad_point
end subroutine madr
end interface
abort=.false.
call mpi_bcast(xyzzyaaae1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Error broadcasting opt_nparam.')
if(xyzzyaabb1)then
call wout('Optimization start')
call wout('==================')
call wout()
if(xyzzyaaac1<xyzzyaaae1)call errstop('VARMIN_DRIVER','About to call N&
&L2SOL but no. of configs < no. of params.')
call dfault(xyzzyaaaj1,xyzzyaaaq1)
xyzzyaaaj1(17)=opt_maxeval
xyzzyaaaj1(18)=opt_maxiter
if(vm_madmin)then
xyzzyaaaq1(32)=min(0.1d0,0.5d0*sqrt(1.0d0/dble(2*xyzzyaaac1)))
else
xyzzyaaaq1(32)=min(0.1d0,0.5d0*sqrt(2.0d0/dble(xyzzyaaac1)))
endif
xyzzyaaaj1(14)=0
xyzzyaaaj1(15)=0
xyzzyaaaj1(19)=0
xyzzyaaaj1(20)=0
xyzzyaaaj1(21)=-1
xyzzyaaaj1(22)=0
xyzzyaaaj1(23)=0
xyzzyaaaj1(24)=0
no_fn_evals=0
xyzzyaaat1=0
xyzzyaaau1=-1
call nl2sno(xyzzyaaac1,xyzzyaaae1,xyzzyaaaf1,madr,xyzzyaaaj1,xyzzyaaaq&
&1)
xyzzyaaad5=.true.
if(xyzzyaaba1)then
call time_report('VP communication: broadcasting finish flag',.true.)
write(6,*)xyzzyaaad5
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(xyzzyaaad5,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Error broadcasting opt_finished.')
xyzzyaaaa5=xyzzyaaaj1(31)
xyzzyaaac5=xyzzyaaaq1(10)*2.d0
call wout('NL2SOL return code : '//trim(i2s(xyzzyaaaj1(1))))
select case(xyzzyaaaj1(1))
case(3)
call wout('Parameter convergence.')
xyzzyaaab5=0
case(4)
call wout('Relative function convergence.')
xyzzyaaab5=0
case(5)
call wout('Parameter and relative function convergence.')
xyzzyaaab5=0
case(6)
call wout('Absolute function convergence.')
xyzzyaaab5=0
case(7)
call wout('Singular convergence.')
xyzzyaaab5=1
case(8)
call wout('False convergence.')
xyzzyaaab5=1
case(9)
call wout('Exceeded OPT_MAXEVAL function evaluations.')
xyzzyaaab5=2
case(10)
call wout('Exceeded OPT_MAXITER iterations.')
xyzzyaaab5=2
case(11)
call wout('Internal stop request.')
xyzzyaaab5=3
abort=.true.
case(13)
call errstop('VARMIN_DRIVER','Couldn''t evaluate vector of residuals a&
&t initial parameters.')
case default
call errstop('VARMIN_DRIVER','Optimization failure -- bad choice of pa&
&rameters / bug?')
end select
call wout()
if(xyzzyaaab5==0)then
call wout('Successful optimization: some degree of confidence in minim&
&um.')
elseif(xyzzyaaab5==1)then
call wout('Optimization complete (numerical instabilities at minimum).&
&')
elseif(xyzzyaaab5==2)then
call wout('Optimization halted: minimum not reached within target accu&
&racy.')
elseif(xyzzyaaab5==3)then
call wout('Continue run as indicated above, if necessary.')
endif
call wout()
if(vm_madmin)then
call wout('Iterations                : '//trim(i2s(xyzzyaaaa5)))
call wout('Function evaluations      : '//trim(i2s(no_fn_evals)))
call wout('Mean abs. dev. reduced to : ',xyzzyaaac5,rfmt='(e16.8)',adj&
&ust=.true.)
else
call wout('Iterations           : '//trim(i2s(xyzzyaaaa5)))
call wout('Function evaluations : '//trim(i2s(no_fn_evals)))
call wout('Variance reduced to  : ',xyzzyaaac5,rfmt='(e16.8)',adjust=.&
&true.)
endif
call wout()
else
if(xyzzyaaba1.and.am_master)then
call wout('STARTING VIRTUAL SLAVE OPTIMIZATION')
call wout('===================================')
call wout()
endif
xyzzyaaae5=.true.
do
if(xyzzyaaba1.and.am_master)then
call time_report('VP communication: broadcasting finish flag',.true.)
read(5,*,iostat=ierror)xyzzyaaad5
if(ierror/=0)call errstop('VARMIN_DRIVER','VP error (opt_finished).')
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(xyzzyaaad5,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting opt_finished in varmin_driver')
if(xyzzyaaad5)then
if(xyzzyaaba1.and.am_master)then
call wout('Optimization finished.')
call wout()
endif
exit
endif
call compute_lsf_slave(xyzzyaaae5,abort)
xyzzyaaae5=.false.
if(abort)exit
enddo
endif
if(xyzzyaaba1.and.am_master)then
call time_report('VP communication: broadcasting params',.true.)
if(xyzzyaabb1)then
write(6,*)xyzzyaaaf1
else
read(5,*,iostat=ierror)xyzzyaaaf1
if(ierror/=0)call errstop('VARMIN_DRIVER','VP error (opt_params).')
endif
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(xyzzyaaaf1,xyzzyaaae1,mpi_double_precision,0,mpi_comm_w&
&orld,ierror)
call checkmpi(ierror,'broadcasting opt_params in varmin_driver')
call put_params(xyzzyaaaf1,xyzzyaaah1,.false.,params_print=xyzzyaaag1)
if(am_master)then
if(opt_info>=3)call param_monitor(xyzzyaaag1)
call write_correlation_header(.true.,.true.)
call update_wfn_casl
call write_casl(':parameters.casl','parameters.'//trim(i2s(max(1,opt_c&
&ycle)))//'.casl',errmsg)
if(len_trim(errmsg)>0)call errstop_master('VARMIN_DRIVER',trim(errmsg)&
&)
endif
if(use_altsamp.and.altsamp==1)then
if(am_master)then
call wout()
if(e0_prop>vmc_optimum_e0.or.ew_prop>vmc_optimum_ew.or.opt_cycle==1)th&
&en
call wout(' Unmodified PDF parameters.')
e0_prop=vmc_optimum_e0
ew_prop=vmc_optimum_ew
else
call wout(' Modified PDF parameters.')
endif
tmpr=r2s(vmc_optimum_e0,'(f10.4)')
tmpr2=r2s(vmc_optimum_ew,'(f10.4)')
tmpr3=r2s(e0_prop,'(f10.4)')
tmpr4=r2s(ew_prop,'(f10.4)')
call wout('(e0,ew)_old -> (e0,ew)_new = ('//trim(tmpr)//','//trim(tmpr&
&2)//') -> ('//trim(tmpr3)//','//trim(tmpr4)//')')
call wout()
vmc_optimum_e0=e0_prop
vmc_optimum_ew=ew_prop
endif
call mpi_bcast(vmc_optimum_e0,1,mpi_double_precision,0,mpi_comm_world,&
&ierror)
call mpi_bcast(vmc_optimum_ew,1,mpi_double_precision,0,mpi_comm_world,&
&ierror)
endif
end subroutine xyzzyaabi1
subroutine compute_lsf_master(params,residual,nl2sol_iteration,been_rs&
&t,abort,bad_point,params_h)
implicit none
integer,intent(in) :: nl2sol_iteration,been_rst
real(dp),intent(in) :: params(xyzzyaaae1)
real(dp),intent(in),optional :: params_h(:)
real(dp),intent(out) :: residual(:)
logical,intent(inout) :: abort,bad_point
character(512) errmsg
integer iparam,xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,indx,no_fil&
&tered
integer,save :: xyzzyaaae7=-1,xyzzyaaaf7=0,xyzzyaaag7=0
real(dp) spread,xyzzyaaah7,ecentre,mon_spread_array(xyzzyaaae1),mon_ec&
&entre_array(xyzzyaaae1),xyzzyaaai7,xyzzyaaaj7,xyzzyaaak7(9)
real(dp),save :: xyzzyaaal7=0.d0,xyzzyaaam7=0.d0,xyzzyaaan7=0.d0,xyzzy&
&aaao7=0.d0,xyzzyaaap7=0.d0,xyzzyaaaq7=0.d0,xyzzyaaar7=0.d0,xyzzyaaas7&
&=0.d0
logical xyzzyaaat7,xyzzyaaau7,xyzzyaaav7,xyzzyaaaw7,xyzzyaaax7,iszero,&
&xyzzyaaay7,xyzzyaaaz7,xyzzyaaba7,xyzzyaabb7
logical,save :: xyzzyaabc7=.false.,xyzzyaabd7=.false.
bad_point=.false.
xyzzyaaat7=.false.
ierror=0
if(xyzzyaaba1)then
call time_report('VP communication: broadcasting finish flag',.true.)
write(6,*)xyzzyaaat7
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(xyzzyaaat7,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'bcasting opt_finished in compute_lsf_master')
xyzzyaaau7=present(params_h)
xyzzyaaav7=.not.xyzzyaaau7.and.nl2sol_iteration==xyzzyaaae7
if(xyzzyaaba1)then
call time_report('VP communication: broadcasting params',.true.)
write(6,*)params,xyzzyaaau7,xyzzyaaav7
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(params,xyzzyaaae1,mpi_double_precision,0,mpi_comm_world&
&,ierror)
call checkmpi(ierror,'bcasting params in compute_lsf_master')
call mpi_bcast(xyzzyaaau7,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'bcasting is_deriv in compute_lsf_master')
call mpi_bcast(xyzzyaaav7,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'bcasting is_correction in compute_lsf_master')
no_fn_evals=no_fn_evals+1
if(no_fn_evals==1)then
xyzzyaabc7=.false.
xyzzyaabd7=.false.
xyzzyaaar7=0.d0
xyzzyaaas7=0.d0
endif
if(xyzzyaabc7)then
if(been_rst/=0)then
if(opt_info==2)then
call rejection_monitor(xyzzyaaal7,xyzzyaaam7,xyzzyaaao7,xyzzyaaap7)
elseif(opt_info>=3)then
call wout('Step rejected.')
call wout()
endif
else
if(opt_info==2)then
call wout('Accepted step:')
call lsf_monitor(xyzzyaaal7,xyzzyaaam7,xyzzyaaan7,xyzzyaaaf7,xyzzyaaao&
&7,xyzzyaaap7,xyzzyaaaq7,xyzzyaaag7)
call wout()
elseif(opt_info>=3)then
call wout('Step accepted.')
call wout()
endif
xyzzyaaao7=xyzzyaaal7
xyzzyaaap7=xyzzyaaam7
xyzzyaaaq7=xyzzyaaan7
xyzzyaaag7=xyzzyaaaf7
call write_correlation_header(.true.,.true.)
call update_wfn_casl
call write_casl(':parameters.casl','parameters.'//trim(i2s(max(1,opt_c&
&ycle)))//'.casl',errmsg)
if(len_trim(errmsg)>0)call errstop_master('COMPUTE_LSF_MASTER',trim(er&
&rmsg))
endif
endif
xyzzyaabc7=.false.
if(.not.xyzzyaaau7)then
if(.not.xyzzyaaav7)then
if(exceeds_time_limit(nl2sol_iteration==0))then
call wout()
call errwarn('COMPUTE_LSF_MASTER','Time limit exceeded or about to be &
&exceeded.')
call wout('CONTINUATION INFO:')
if(xyzzyaaae7>opt_maxiter/3.and.xyzzyaaae7>1)then
if(isopt_vmc.or.isvmc_opt)then
call wout(' Suggested action: continue run directly')
if(opt_cycle==opt_cycles)then
call wout(' Set RUNTYPE = vmc')
if(old_input)then
call wout(' Set NWRCON = 0')
else
call wout(' Set VMC_NCONFIG_WRITE = 0')
endif
call wout(' Set WRITEOUT_VMC_HIST = F')
call wout(' Move correlation.out.'//trim(i2s(opt_cycle))//' to correla&
&tion.data')
else
if(isopt_vmc)call wout(' Set RUNTYPE = vmc_opt')
call wout(' Move correlation.out.'//trim(i2s(opt_cycle))//' to correla&
&tion.data')
call wout(' Set OPT_CYCLES = '//trim(i2s(opt_cycles-opt_cycle)))
endif
call wout(' Move config.out to config.in')
else
call wout(' Suggested action: none.')
call wout(' No need to continue.')
endif
else
if(isopt_vmc.or.isvmc_opt)then
if(opt_cycle>1)then
call wout(' Suggested action: roll back and continue')
if(isvmc_opt)call wout(' Set RUNTYPE = opt_vmc')
call wout(' Move correlation.out.'//trim(i2s(opt_cycle-1))//' to corre&
&lation.data')
call wout(' Set OPT_CYCLES = '//trim(i2s(opt_cycles-opt_cycle+1)))
elseif(isvmc_opt)then
call wout(' Suggested action: roll back and continue')
call wout(' Set RUNTYPE = opt_vmc')
else
call wout(' Suggested action: restart with more time')
call wout(' PROBLEM: first optimization cycle won''t fit in time slot!&
&')
call wout(' You must re-run using a longer time slot.')
endif
else
call wout(' Suggested action: restart with more time')
call wout(' PROBLEM: optimization run won''t fit in time slot!')
call wout(' You must re-run using a longer time slot.')
endif
endif
call wout()
abort=.true.
return
endif
endif
if(opt_info>=2.and..not.xyzzyaaav7)call opt_monitor_header(nl2sol_iter&
&ation)
call put_params(params,xyzzyaaah1,.false.,params_print=xyzzyaaag1,bad_&
&params=bad_point)
if(.not.bad_point)then
if(opt_info>=2)call param_monitor(xyzzyaaag1)
call invalidate_params(xyzzyaaah1)
call invalidate_params_energy
xyzzyaaaw7=.false.
xyzzyaaax7=.false.
iszero=.false.
do xyzzyaaaa7=1,xyzzyaaaa1
if(.not.complex_wf)then
call compute_config_energy(xyzzyaaaa7,xyzzyaaak1(xyzzyaaaa7),xyzzyaaai&
&7,xyzzyaaam1(xyzzyaaaa7),no_fn_evals==1,xyzzyaaay7,xyzzyaaaz7,xyzzyaa&
&ba7)
else
call compute_config_energy(xyzzyaaaa7,xyzzyaaak1(xyzzyaaaa7),xyzzyaaal&
&1(xyzzyaaaa7),xyzzyaaam1(xyzzyaaaa7),no_fn_evals==1,xyzzyaaay7,xyzzya&
&aaz7,xyzzyaaba7)
endif
xyzzyaaaw7=xyzzyaaaw7.or.xyzzyaaay7
xyzzyaaax7=xyzzyaaax7.or.xyzzyaaaz7
iszero=iszero.or.xyzzyaaba7
enddo
xyzzyaaat1=xyzzyaaat1+1
if(xyzzyaaaw7.or.xyzzyaaax7.or.iszero)bad_point=.true.
xyzzyaabb7=bad_point
call mpi_reduce(xyzzyaabb7,bad_point,1,mpi_logical,mpi_lor,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'reducing bad_point is compute_lsf_master')
call mpi_bcast(bad_point,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bad_point in compute_lsf_master')
if(xyzzyaaba1)then
call time_report('VP communication: gathering bad_point',.true.)
do xyzzyaaab7=1,virtual_nnodes-1
write(6,*)xyzzyaaab7
read(5,*,iostat=ierror)xyzzyaabb7
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (bad_point1).&
&')
if(xyzzyaabb7)bad_point=.true.
enddo
write(6,*)bad_point
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
call mpi_bcast(bad_point,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bad_point in compute_lsf_master')
endif
if(.not.bad_point)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaak1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering opt_ecfg in compute_lsf_master')
if(complex_wf)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaal1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering opt_eicfg in compute_lsf_master')
endif
if(xyzzyaaaz1)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaam1,xy&
&zzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering opt_wcfg in compute_lsf_master')
endif
if(xyzzyaaba1)then
call time_report('VP communication: gathering energies',.true.)
xyzzyaaac7=xyzzyaaab1
do xyzzyaaab7=1,virtual_nnodes-1
write(6,*)xyzzyaaab7
read(5,*,iostat=ierror)xyzzyaaad7
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (kcfg).')
read(5,*,iostat=ierror)(xyzzyaaak1(xyzzyaaaa7),xyzzyaaaa7=xyzzyaaac7+1&
&,xyzzyaaac7+xyzzyaaad7)
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (opt_ecfg).')
if(complex_wf)then
read(5,*,iostat=ierror)(xyzzyaaal1(xyzzyaaaa7),xyzzyaaaa7=xyzzyaaac7+1&
&,xyzzyaaac7+xyzzyaaad7)
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (opt_eicfg).'&
&)
endif
if(xyzzyaaaz1)then
read(5,*,iostat=ierror)(xyzzyaaam1(xyzzyaaaa7),xyzzyaaaa7=xyzzyaaac7+1&
&,xyzzyaaac7+xyzzyaaad7)
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (opt_wcfg).')
endif
xyzzyaaac7=xyzzyaaac7+xyzzyaaad7
enddo
if(xyzzyaaac7/=virtual_nconfig)call errstop('COMPUTE_LSF_MASTER','VP e&
&rror (config count).')
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
endif
endif
if(bad_point)then
if(opt_info>=2)then
call wout('Bad parameter set.')
call wout()
endif
return
endif
if(xyzzyaabd1)then
if(.not.xyzzyaabd7)then
call ppmcta_partial_sum(xyzzyaaac1,xyzzyaaak1,xyzzyaabe1,xyzzyaabf1,xy&
&zzyaaak7)
call ppmcta_fit(xyzzyaaak7,xyzzyaaar7,xyzzyaaas7)
xyzzyaabd7=.true.
call wout('MCTA post-processing fit parameters for this run:')
call wout(' a = ',xyzzyaaar7)
call wout(' b = ',xyzzyaaas7)
call wout()
endif
xyzzyaaak1(1:xyzzyaaac1)=xyzzyaaak1(1:xyzzyaaac1)+xyzzyaaar7*xyzzyaabe&
&1(1:xyzzyaaac1)+xyzzyaaas7*xyzzyaabf1(1:xyzzyaaac1)
endif
xyzzyaaah7=0.d0
if(xyzzyaaaz1.and.opt_info>=2)then
xyzzyaaah7=0.d0
xyzzyaaaj7=dble(xyzzyaaac1)/sum(xyzzyaaam1(1:xyzzyaaac1))
do xyzzyaaaa7=1,xyzzyaaac1
xyzzyaaah7=xyzzyaaah7+(xyzzyaaam1(xyzzyaaaa7)*xyzzyaaaj7-1.d0)**2
enddo
xyzzyaaah7=xyzzyaaah7/dble(xyzzyaaac1-1)
endif
if(.not.complex_wf)then
call target_function(xyzzyaaak1,xyzzyaaam1,residual,spread,ecentre,no_&
&filtered,bad_point=bad_point,report=.true.)
else
call target_function(xyzzyaaak1,xyzzyaaam1,residual,spread,ecentre,no_&
&filtered,eivector=xyzzyaaal1,bad_point=bad_point,report=.true.)
endif
if(use_altsamp)then
xyzzyaabb7=bad_point
call mpi_reduce(xyzzyaabb7,bad_point,1,mpi_logical,mpi_lor,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'reducing bad_point is compute_lsf_master')
call mpi_bcast(bad_point,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bad_point in compute_lsf_master')
if(xyzzyaaba1)then
call time_report('VP communication: gathering bad_point',.true.)
do xyzzyaaab7=1,virtual_nnodes-1
write(6,*)xyzzyaaab7
read(5,*,iostat=ierror)xyzzyaabb7
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (bad_point1).&
&')
if(xyzzyaabb7)bad_point=.true.
enddo
write(6,*)bad_point
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
call mpi_bcast(bad_point,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bad_point in compute_lsf_master')
endif
endif
if(no_fn_evals==1)then
if(opt_info>=2)call lsf_monitor(spread,ecentre,xyzzyaaah7,no_filtered)
xyzzyaaao7=spread
xyzzyaaap7=ecentre
xyzzyaaaq7=xyzzyaaah7
else
if(.not.xyzzyaaav7)xyzzyaaae7=nl2sol_iteration
if(opt_info>=3)then
call lsf_monitor(spread,ecentre,xyzzyaaah7,no_filtered,xyzzyaaao7,xyzz&
&yaaap7,xyzzyaaaq7,xyzzyaaag7)
endif
xyzzyaaal7=spread
xyzzyaaam7=ecentre
xyzzyaaan7=xyzzyaaah7
xyzzyaaaf7=no_filtered
xyzzyaabc7=.true.
endif
else
if(xyzzyaaba1)then
call time_report('VP communication: broadcasting parameters',.true.)
write(6,*)params_h
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(params_h,xyzzyaaae1,mpi_double_precision,0,mpi_comm_wor&
&ld,ierror)
call checkmpi(ierror,'Error broadcasting params_h.')
xyzzyaaaw7=.false.
xyzzyaaax7=.false.
iszero=.false.
do xyzzyaaaa7=1,xyzzyaaaa1
call loop_time_estimate('Computing derivatives.',xyzzyaaaa7,xyzzyaaaa1&
&)
if(complex_wf)then
call compute_config_derivs(xyzzyaaaa7,params,params_h,xyzzyaaan1(xyzzy&
&aaaa7,:),xyzzyaaap1(xyzzyaaaa7,:),xyzzyaaay7,xyzzyaaaz7,xyzzyaaba7,ei&
&_h=xyzzyaaao1(xyzzyaaaa7,:))
else
call compute_config_derivs(xyzzyaaaa7,params,params_h,xyzzyaaan1(xyzzy&
&aaaa7,:),xyzzyaaap1(xyzzyaaaa7,:),xyzzyaaay7,xyzzyaaaz7,xyzzyaaba7)
endif
xyzzyaaaw7=xyzzyaaaw7.or.xyzzyaaay7
xyzzyaaax7=xyzzyaaax7.or.xyzzyaaaz7
iszero=iszero.or.xyzzyaaba7
enddo
call loop_time_estimate('Done.')
call qmc_barrier
call mpi_reduce(xyzzyaaaw7,xyzzyaaay7,1,mpi_logical,mpi_lor,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'reducing isNaN is compute_lsf_master')
call mpi_reduce(xyzzyaaax7,xyzzyaaaz7,1,mpi_logical,mpi_lor,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'reducing isInf is compute_lsf_master')
call mpi_reduce(iszero,xyzzyaaba7,1,mpi_logical,mpi_lor,0,mpi_comm_wor&
&ld,ierror)
call checkmpi(ierror,'reducing isZero is compute_lsf_master')
do iparam=1,xyzzyaaae1
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaan1(:,&
&iparam),xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Error gathering opt_decfg.')
if(complex_wf)then
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaao1(:,&
&iparam),xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Error gathering opt_deicfg.')
endif
enddo
if(xyzzyaaaz1)then
do iparam=1,xyzzyaaae1
call mpi_gather_in_place(xyzzyaaaa1,mpi_double_precision,xyzzyaaap1(:,&
&iparam),xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Error gathering opt_dwcfg.')
enddo
endif
if(xyzzyaaba1)then
call time_report('VP communication: gathering derivatives',.true.)
xyzzyaaac7=xyzzyaaab1
do xyzzyaaab7=1,virtual_nnodes-1
write(6,*)xyzzyaaab7
read(5,*,iostat=ierror)xyzzyaaay7
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (NaN).')
if(xyzzyaaay7)xyzzyaaaw7=.true.
read(5,*,iostat=ierror)xyzzyaaaz7
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (Inf).')
if(xyzzyaaaz7)xyzzyaaax7=.true.
read(5,*,iostat=ierror)xyzzyaaba7
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (Zero).')
if(xyzzyaaba7)iszero=.true.
read(5,*,iostat=ierror)xyzzyaaad7
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (kcfg).')
read(5,*,iostat=ierror)((xyzzyaaan1(xyzzyaaaa7,iparam),xyzzyaaaa7=xyzz&
&yaaac7+1,xyzzyaaac7+xyzzyaaad7),iparam=1,xyzzyaaae1)
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (opt_decfg).'&
&)
if(complex_wf)then
read(5,*,iostat=ierror)((xyzzyaaao1(xyzzyaaaa7,iparam),xyzzyaaaa7=xyzz&
&yaaac7+1,xyzzyaaac7+xyzzyaaad7),iparam=1,xyzzyaaae1)
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (opt_deicfg).&
&')
endif
if(xyzzyaaaz1)then
read(5,*,iostat=ierror)((xyzzyaaap1(xyzzyaaaa7,iparam),xyzzyaaaa7=xyzz&
&yaaac7+1,xyzzyaaac7+xyzzyaaad7),iparam=1,xyzzyaaae1)
if(ierror/=0)call errstop('COMPUTE_LSF_MASTER','VP error (opt_dwcfg).'&
&)
endif
xyzzyaaac7=xyzzyaaac7+xyzzyaaad7
enddo
if(xyzzyaaac7/=virtual_nconfig)call errstop('COMPUTE_LSF_MASTER','VP e&
&rror (config count).')
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
if(xyzzyaaaw7.or.xyzzyaaax7)call errstop('COMPUTE_LSF_MASTER','Floatin&
&g-point exception reported when calculating energy derivatives.  NL2S&
&OL does not support instability feedback for derivatives, so stopping&
&.')
if(iszero)call errstop('COMPUTE_LSF_MASTER','Zero wave function value &
&encountered when calculating energy derivatives.  NL2SOL does not sup&
&port instability feedback for derivatives, so stopping.')
if(xyzzyaabd1)then
do iparam=1,xyzzyaaae1
xyzzyaaan1(1:xyzzyaaac1,iparam)=xyzzyaaan1(1:xyzzyaaac1,iparam)+xyzzya&
&aar7*xyzzyaabe1(1:xyzzyaaac1)+xyzzyaaas7*xyzzyaabf1(1:xyzzyaaac1)
enddo
endif
indx=1
do iparam=1,xyzzyaaae1
if(.not.complex_wf)then
call target_function(xyzzyaaan1(:,iparam),xyzzyaaap1(:,iparam),residua&
&l(indx:indx+xyzzyaaac1-1),mon_spread_array(iparam),mon_ecentre_array(&
&iparam),no_filtered,report=.false.)
else
call target_function(xyzzyaaan1(:,iparam),xyzzyaaap1(:,iparam),residua&
&l(indx:indx+xyzzyaaac1-1),mon_spread_array(iparam),mon_ecentre_array(&
&iparam),no_filtered,eivector=xyzzyaaao1(:,iparam),report=.false.)
endif
indx=indx+xyzzyaaac1
enddo
if(opt_info>=3)call deriv_monitor(params,params_h,mon_spread_array,xyz&
&zyaaao7,mon_ecentre_array,xyzzyaaap7)
endif
end subroutine compute_lsf_master
subroutine target_function(evector,wvector,residual,spread,ecentre,no_&
&filtered,eivector,bad_point,report)
implicit none
integer,intent(out) :: no_filtered
real(dp),intent(in) :: evector(xyzzyaaac1)
real(dp),intent(in),optional :: eivector(xyzzyaaac1)
real(dp),intent(out) :: residual(xyzzyaaac1),spread,ecentre
real(dp),intent(inout) :: wvector(xyzzyaaac1)
logical,intent(in),optional :: report
logical,intent(inout),optional :: bad_point
integer xyzzyaaaa8
real(dp) xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaaaf8,xyzzya&
&aag8,xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaak8,xyzzyaaal8,xyzzyaaam&
&8,xyzzyaaan8
no_filtered=0
if(use_altsamp)then
xyzzyaaac8=sum(wvector(1:xyzzyaaac1))
if(xyzzyaaac8<xyzzyaaay1)call errstop('TARGET_FUNCTION','Problem summi&
&ng weights: total weight is zero.')
xyzzyaaan8=dble(xyzzyaaac1)/xyzzyaaac8
call dscal(xyzzyaaac1,xyzzyaaan8,wvector,1)
if(present(bad_point))then
xyzzyaaaw1=0.05d0*dble(xyzzyaaac1)
if(maxval(wvector(:))>xyzzyaaaw1)then
bad_point=.true.
call wout(' Reject move since Max[w] too large...')
endif
xyzzyaaal8=sum((wvector(:)-1.d0)**2)/dble(xyzzyaaac1*(xyzzyaaac1-1))
xyzzyaaam8=(0.3d0)**2
if(xyzzyaaal8>xyzzyaaam8)then
bad_point=.true.
call wout(' Reject move since Var[w] too large...')
endif
if(no_fn_evals==1)bad_point=.false.
endif
if(vm_w_max>0.d0)then
do xyzzyaaaa8=1,xyzzyaaac1
if(wvector(xyzzyaaaa8)>vm_w_max)wvector(xyzzyaaaa8)=vm_w_max
if(wvector(xyzzyaaaa8)<vm_w_min)wvector(xyzzyaaaa8)=vm_w_min
enddo
xyzzyaaac8=sum(wvector(1:xyzzyaaac1))
xyzzyaaab8=dble(xyzzyaaac1)/xyzzyaaac8
call dscal(xyzzyaaac1,xyzzyaaab8,wvector,1)
endif
xyzzyaaae8=1.d0/sqrt(dble(xyzzyaaac1-1))
ecentre=ddot(xyzzyaaac1,evector,1,wvector,1)*xyzzyaaad1
if(.not.vm_use_e_guess)then
do xyzzyaaaa8=1,xyzzyaaac1
if(.not.complex_wf)then
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-ecentre)*wvector(xyzzyaaaa8)
residual(xyzzyaaaa8)=residual(xyzzyaaaa8)*xyzzyaaae8
else
xyzzyaaak8=ddot(xyzzyaaac1,eivector,1,wvector,1)*xyzzyaaad1
xyzzyaaal8=evector(xyzzyaaaa8)-ecentre
xyzzyaaam8=eivector(xyzzyaaaa8)-xyzzyaaak8
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt(xyzzyaaal8*xyzzyaaal8+&
&xyzzyaaam8*xyzzyaaam8)*wvector(xyzzyaaaa8)
residual(xyzzyaaaa8)=residual(xyzzyaaaa8)*xyzzyaaae8
endif
enddo
else
do xyzzyaaaa8=1,xyzzyaaac1
if(.not.complex_wf)then
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-xyzzyaaav1)*wvector(xyzzyaaa&
&a8)
residual(xyzzyaaaa8)=residual(xyzzyaaaa8)*xyzzyaaae8
else
xyzzyaaak8=ddot(xyzzyaaac1,eivector,1,wvector,1)*xyzzyaaad1
xyzzyaaal8=evector(xyzzyaaaa8)-xyzzyaaav1
xyzzyaaam8=eivector(xyzzyaaaa8)
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt(xyzzyaaal8*xyzzyaaal8+&
&xyzzyaaam8*xyzzyaaam8)*wvector(xyzzyaaaa8)
residual(xyzzyaaaa8)=residual(xyzzyaaaa8)*xyzzyaaae8
endif
enddo
endif
if(altsamp==1)then
e0_prop=ecentre
ew_prop=sum(((evector-ecentre)*wvector)**2)
ew_prop=ew_prop/(dble(xyzzyaaac1*(xyzzyaaac1-1)))
ew_prop=sqrt(ew_prop)*2.d0
endif
elseif(vm_reweight)then
xyzzyaaac8=sum(wvector(1:xyzzyaaac1))
if(xyzzyaaac8<xyzzyaaay1)call errstop('TARGET_FUNCTION','Problem summi&
&ng weights: total weight is zero.')
xyzzyaaab8=dble(xyzzyaaac1)/xyzzyaaac8
call dscal(xyzzyaaac1,xyzzyaaab8,wvector,1)
if(vm_w_max>0.d0)then
do xyzzyaaaa8=1,xyzzyaaac1
if(wvector(xyzzyaaaa8)>vm_w_max)wvector(xyzzyaaaa8)=vm_w_max
if(wvector(xyzzyaaaa8)<vm_w_min)wvector(xyzzyaaaa8)=vm_w_min
enddo
xyzzyaaac8=sum(wvector(1:xyzzyaaac1))
xyzzyaaab8=dble(xyzzyaaac1)/xyzzyaaac8
call dscal(xyzzyaaac1,xyzzyaaab8,wvector,1)
endif
xyzzyaaad8=0.d0
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaad8=xyzzyaaad8+wvector(xyzzyaaaa8)**2
enddo
xyzzyaaad8=dble(xyzzyaaac1)-xyzzyaaad8*xyzzyaaad1
xyzzyaaae8=1.d0/sqrt(xyzzyaaad8)
ecentre=ddot(xyzzyaaac1,evector,1,wvector,1)*xyzzyaaad1
if(.not.vm_use_e_guess)then
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-ecentre)*sqrt(wvector(xyzzya&
&aaa8))*xyzzyaaae8
enddo
else
xyzzyaaak8=ddot(xyzzyaaac1,eivector,1,wvector,1)*xyzzyaaad1
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-ecentre
xyzzyaaam8=eivector(xyzzyaaaa8)-xyzzyaaak8
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt((xyzzyaaal8*xyzzyaaal8&
&+xyzzyaaam8*xyzzyaaam8)*wvector(xyzzyaaaa8))*xyzzyaaae8
enddo
endif
else
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-xyzzyaaav1)*sqrt(wvector(xyz&
&zyaaaa8))*xyzzyaaae8
enddo
else
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-xyzzyaaav1
xyzzyaaam8=eivector(xyzzyaaaa8)
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt((xyzzyaaal8*xyzzyaaal8&
&+xyzzyaaam8*xyzzyaaam8)*wvector(xyzzyaaaa8))*xyzzyaaae8
enddo
endif
endif
elseif(vm_filter)then
xyzzyaaae8=1.d0/sqrt(dble(xyzzyaaac1-1))
ecentre=sum(evector(1:xyzzyaaac1))*xyzzyaaad1
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-ecentre)*xyzzyaaae8
enddo
else
xyzzyaaak8=sum(eivector(1:xyzzyaaac1))*xyzzyaaad1
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-ecentre
xyzzyaaam8=eivector(xyzzyaaaa8)-xyzzyaaak8
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt(xyzzyaaal8*xyzzyaaal8+&
&xyzzyaaam8*xyzzyaaam8)*xyzzyaaae8
enddo
endif
spread=ddot(xyzzyaaac1,residual,1,residual,1)
wvector=1.d0
xyzzyaaai8=sqrt(spread)
xyzzyaaaf8=xyzzyaaai8*vm_filter_thres
xyzzyaaag8=-1.d0
if(vm_filter_width/=0.d0.and.xyzzyaaai8/=0.d0)xyzzyaaag8=1.d0 /(xyzzya&
&aai8*vm_filter_width)
do xyzzyaaaa8=1,xyzzyaaac1
if(.not.complex_wf)then
xyzzyaaaj8=abs(evector(xyzzyaaaa8)-ecentre)-xyzzyaaaf8
else
xyzzyaaal8=evector(xyzzyaaaa8)-ecentre
xyzzyaaam8=eivector(xyzzyaaaa8)-xyzzyaaak8
xyzzyaaaj8=sqrt(xyzzyaaal8*xyzzyaaal8+xyzzyaaam8*xyzzyaaam8)-xyzzyaaaf&
&8
endif
if(xyzzyaaaj8>0.d0)then
if(xyzzyaaag8>=0.d0)then
xyzzyaaah8=xyzzyaaaj8*xyzzyaaag8
wvector(xyzzyaaaa8)=exp(-xyzzyaaah8**2)
else
wvector(xyzzyaaaa8)=0.d0
endif
no_filtered=no_filtered+1
endif
enddo
xyzzyaaac8=sum(wvector(1:xyzzyaaac1))
if(xyzzyaaac8<xyzzyaaay1)call errstop('TARGET_FUNCTION','Problem summi&
&ng weights: total weight is zero <2>.')
xyzzyaaab8=real(xyzzyaaac1,dp)/xyzzyaaac8
call dscal(xyzzyaaac1,xyzzyaaab8,wvector,1)
xyzzyaaad8=real(xyzzyaaac1,dp)-ddot(xyzzyaaac1,wvector(1),1,wvector(1)&
&,1)*xyzzyaaad1
xyzzyaaae8=1.d0/sqrt(xyzzyaaad8)
ecentre=ddot(xyzzyaaac1,evector,1,wvector,1)*xyzzyaaad1
if(.not.vm_use_e_guess)then
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-ecentre)*sqrt(wvector(xyzzya&
&aaa8))*xyzzyaaae8
enddo
else
xyzzyaaak8=ddot(xyzzyaaac1,eivector,1,wvector,1)*xyzzyaaad1
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-ecentre
xyzzyaaam8=eivector(xyzzyaaaa8)-xyzzyaaak8
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt((xyzzyaaal8*xyzzyaaal8&
&+xyzzyaaam8*xyzzyaaam8)*wvector(xyzzyaaaa8))*xyzzyaaae8
enddo
endif
else
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-xyzzyaaav1)*sqrt(wvector(xyz&
&zyaaaa8))*xyzzyaaae8
enddo
else
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-xyzzyaaav1
xyzzyaaam8=eivector(xyzzyaaaa8)
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt((xyzzyaaal8*xyzzyaaal8&
&+xyzzyaaam8*xyzzyaaam8)*wvector(xyzzyaaaa8))*xyzzyaaae8
enddo
endif
endif
elseif(vm_madmin)then
if(.not.vm_use_e_guess)then
ecentre=median(xyzzyaaac1,evector)
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=sqrt(abs(evector(xyzzyaaaa8)-ecentre)*xyzzyaaad1)
enddo
else
xyzzyaaak8=median(xyzzyaaac1,eivector)
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-ecentre
xyzzyaaam8=eivector(xyzzyaaaa8)-xyzzyaaak8
residual(xyzzyaaaa8)=sqrt(sqrt(xyzzyaaal8*xyzzyaaal8+xyzzyaaam8*xyzzya&
&aam8)*xyzzyaaad1)
enddo
endif
else
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=sqrt(abs(evector(xyzzyaaaa8)-xyzzyaaav1)*xyzzyaaa&
&d1)
enddo
else
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-xyzzyaaav1
xyzzyaaam8=eivector(xyzzyaaaa8)
residual(xyzzyaaaa8)=sqrt(sqrt(xyzzyaaal8*xyzzyaaam8+xyzzyaaam8*xyzzya&
&aam8)*xyzzyaaad1)
enddo
endif
endif
else
xyzzyaaae8=1.d0/sqrt(dble(xyzzyaaac1-1))
ecentre=sum(evector(1:xyzzyaaac1))*xyzzyaaad1
if(.not.vm_use_e_guess)then
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-ecentre)*xyzzyaaae8
enddo
else
xyzzyaaak8=sum(eivector(1:xyzzyaaac1))*xyzzyaaad1
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-ecentre
xyzzyaaam8=eivector(xyzzyaaaa8)-xyzzyaaak8
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt(xyzzyaaal8*xyzzyaaal8+&
&xyzzyaaam8*xyzzyaaam8)*xyzzyaaae8
enddo
endif
else
if(.not.complex_wf)then
do xyzzyaaaa8=1,xyzzyaaac1
residual(xyzzyaaaa8)=(evector(xyzzyaaaa8)-xyzzyaaav1)*xyzzyaaae8
enddo
else
do xyzzyaaaa8=1,xyzzyaaac1
xyzzyaaal8=evector(xyzzyaaaa8)-xyzzyaaav1
xyzzyaaam8=eivector(xyzzyaaaa8)
residual(xyzzyaaaa8)=sign(1.d0,xyzzyaaal8)*sqrt(xyzzyaaal8*xyzzyaaal8+&
&xyzzyaaam8*xyzzyaaam8)*xyzzyaaae8
enddo
endif
endif
endif
spread=ddot(xyzzyaaac1,residual,1,residual,1)
end subroutine target_function
subroutine compute_lsf_slave(first_call,abort)
implicit none
logical,intent(in) :: first_call
logical,intent(out) :: abort
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,xyzzyaaae9
real(dp) params(xyzzyaaae1),xyzzyaaaf9(xyzzyaaae1),xyzzyaaag9
logical xyzzyaaah9,xyzzyaaai9,xyzzyaaaj9,xyzzyaaak9,xyzzyaaal9,xyzzyaa&
&am9,xyzzyaaan9,xyzzyaaao9,xyzzyaaap9,xyzzyaaaq9
xyzzyaaap9=.false.
abort=.false.
if(xyzzyaaba1.and.am_master)then
call time_report('VP communication: broadcasting parameters',.true.)
read(5,*,iostat=ierror)params,xyzzyaaah9,xyzzyaaai9
if(ierror/=0)call errstop('COMPUTE_LSF_SLAVE','VP error (params).')
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(params,xyzzyaaae1,mpi_double_precision,0,mpi_comm_world&
&,ierror)
call checkmpi(ierror,'bcasting params in compute_lsf_slave')
call mpi_bcast(xyzzyaaah9,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'bcasting is_deriv in compute_lsf_slave')
call mpi_bcast(xyzzyaaai9,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'bcasting is_correction in compute_lsf_slave')
if(.not.xyzzyaaah9)then
if(.not.xyzzyaaai9)then
if(exceeds_time_limit(.false.))then
abort=.true.
return
endif
endif
call put_params(params,xyzzyaaah1,.false.,bad_params=xyzzyaaap9)
if(.not.xyzzyaaap9)then
call invalidate_params(xyzzyaaah1)
call invalidate_params_energy
xyzzyaaaj9=.false.
xyzzyaaak9=.false.
xyzzyaaal9=.false.
do xyzzyaaab9=1,xyzzyaaaa1
if(.not.complex_wf)then
call compute_config_energy(xyzzyaaab9,xyzzyaaak1(xyzzyaaab9),xyzzyaaag&
&9,xyzzyaaam1(xyzzyaaab9),first_call,xyzzyaaam9,xyzzyaaan9,xyzzyaaao9)
else
call compute_config_energy(xyzzyaaab9,xyzzyaaak1(xyzzyaaab9),xyzzyaaal&
&1(xyzzyaaab9),xyzzyaaam1(xyzzyaaab9),first_call,xyzzyaaam9,xyzzyaaan9&
&,xyzzyaaao9)
endif
xyzzyaaaj9=xyzzyaaaj9.or.xyzzyaaam9
xyzzyaaak9=xyzzyaaak9.or.xyzzyaaan9
xyzzyaaal9=xyzzyaaal9.or.xyzzyaaao9
enddo
if(xyzzyaaaj9.or.xyzzyaaak9.or.xyzzyaaal9)xyzzyaaap9=.true.
xyzzyaaaq9=xyzzyaaap9
call mpi_reduce(xyzzyaaaq9,xyzzyaaap9,1,mpi_logical,mpi_lor,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'reducing bad_point is compute_lsf_slave')
call mpi_bcast(xyzzyaaap9,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bad_point in compute_lsf_slave')
if(xyzzyaaba1)then
if(am_master)then
call time_report('VP communication: gathering bad_point',.true.)
do xyzzyaaac9=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaae9)xyzzyaaad9
if(xyzzyaaae9/=0)call errstop('COMPUTE_LSF_SLAVE','VP error (jnode).')
if(xyzzyaaad9/=virtual_node)cycle
write(6,*)xyzzyaaap9
enddo
read(5,*,iostat=xyzzyaaae9)xyzzyaaap9
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(xyzzyaaap9,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bad_point in compute_lsf_slave')
endif
if(.not.xyzzyaaap9)then
call mpi_gather(xyzzyaaak1,xyzzyaaaa1,mpi_double_precision,xyzzyaaak1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering opt_ecfg in compute_lsf_slave')
if(complex_wf)then
call mpi_gather(xyzzyaaal1,xyzzyaaaa1,mpi_double_precision,xyzzyaaal1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering opt_eicfg in compute_lsf_slave')
endif
if(xyzzyaaaz1)then
call mpi_gather(xyzzyaaam1,xyzzyaaaa1,mpi_double_precision,xyzzyaaam1,&
&xyzzyaaaa1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'gathering opt_wcfg in compute_lsf_slave')
endif
if(xyzzyaaba1.and.am_master)then
call time_report('VP communication: gathering energies',.true.)
do xyzzyaaac9=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaae9)xyzzyaaad9
if(xyzzyaaae9/=0)call errstop('COMPUTE_LSF_SLAVE','VP error (jnode).')
if(xyzzyaaad9/=virtual_node)cycle
write(6,*)xyzzyaaaj9
write(6,*)xyzzyaaak9
write(6,*)xyzzyaaal9
write(6,*)xyzzyaaab1
write(6,*)(xyzzyaaak1(xyzzyaaab9),xyzzyaaab9=1,xyzzyaaab1)
if(complex_wf)write(6,*)(xyzzyaaal1(xyzzyaaab9),xyzzyaaab9=1,xyzzyaaab&
&1)
if(xyzzyaaaz1)write(6,*)(xyzzyaaam1(xyzzyaaab9),xyzzyaaab9=1,xyzzyaaab&
&1)
enddo
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
endif
endif
if(xyzzyaaap9)return
if(use_altsamp)then
xyzzyaaaq9=xyzzyaaap9
call mpi_reduce(xyzzyaaaq9,xyzzyaaap9,1,mpi_logical,mpi_lor,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'reducing bad_point is compute_lsf_slave')
call mpi_bcast(xyzzyaaap9,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bad_point in compute_lsf_slave')
if(xyzzyaaba1)then
if(am_master)then
call time_report('VP communication: gathering bad_point',.true.)
do xyzzyaaac9=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaae9)xyzzyaaad9
if(xyzzyaaae9/=0)call errstop('COMPUTE_LSF_SLAVE','VP error (jnode).')
if(xyzzyaaad9/=virtual_node)cycle
write(6,*)xyzzyaaap9
enddo
read(5,*,iostat=xyzzyaaae9)xyzzyaaap9
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(xyzzyaaap9,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcasting bad_point in compute_lsf_slave')
endif
endif
else
if(xyzzyaaba1.and.am_master)then
call time_report('VP communication: broadcasting parameters',.true.)
read(5,*,iostat=ierror)xyzzyaaaf9
if(ierror/=0)call errstop('COMPUTE_LSF_SLAVE','VP error (params_h).')
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
call mpi_bcast(xyzzyaaaf9,xyzzyaaae1,mpi_double_precision,0,mpi_comm_w&
&orld,ierror)
call checkmpi(ierror,'broadcasting params_h in compute_lsf_slave')
xyzzyaaaj9=.false.
xyzzyaaak9=.false.
xyzzyaaal9=.false.
do xyzzyaaab9=1,xyzzyaaaa1
if(.not.complex_wf)then
call compute_config_derivs(xyzzyaaab9,params,xyzzyaaaf9,xyzzyaaan1(xyz&
&zyaaab9,:),xyzzyaaap1(xyzzyaaab9,:),xyzzyaaam9,xyzzyaaan9,xyzzyaaao9)
else
call compute_config_derivs(xyzzyaaab9,params,xyzzyaaaf9,xyzzyaaan1(xyz&
&zyaaab9,:),xyzzyaaap1(xyzzyaaab9,:),xyzzyaaam9,xyzzyaaan9,xyzzyaaao9,&
&ei_h=xyzzyaaao1(xyzzyaaab9,:))
endif
xyzzyaaaj9=xyzzyaaaj9.or.xyzzyaaam9
xyzzyaaak9=xyzzyaaak9.or.xyzzyaaan9
xyzzyaaal9=xyzzyaaal9.or.xyzzyaaao9
enddo
call qmc_barrier
call mpi_reduce(xyzzyaaaj9,xyzzyaaam9,1,mpi_logical,mpi_lor,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'reducing isNaN is compute_lsf_slave')
call mpi_reduce(xyzzyaaak9,xyzzyaaan9,1,mpi_logical,mpi_lor,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'reducing isInf is compute_lsf_slave')
call mpi_reduce(xyzzyaaal9,xyzzyaaao9,1,mpi_logical,mpi_lor,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'reducing isZero is compute_lsf_slave')
do xyzzyaaaa9=1,xyzzyaaae1
call mpi_gather(xyzzyaaan1(:,xyzzyaaaa9),xyzzyaaaa1,mpi_double_precisi&
&on,xyzzyaaan1(:,xyzzyaaaa9),xyzzyaaaa1,mpi_double_precision,0,mpi_com&
&m_world,ierror)
call checkmpi(ierror,'Error gathering opt_decfg.')
if(complex_wf)then
call mpi_gather(xyzzyaaao1(:,xyzzyaaaa9),xyzzyaaaa1,mpi_double_precisi&
&on,xyzzyaaao1(:,xyzzyaaaa9),xyzzyaaaa1,mpi_double_precision,0,mpi_com&
&m_world,ierror)
call checkmpi(ierror,'Error gathering opt_deicfg.')
endif
enddo
if(xyzzyaaaz1)then
do xyzzyaaaa9=1,xyzzyaaae1
call mpi_gather(xyzzyaaap1(:,xyzzyaaaa9),xyzzyaaaa1,mpi_double_precisi&
&on,xyzzyaaap1(:,xyzzyaaaa9),xyzzyaaaa1,mpi_double_precision,0,mpi_com&
&m_world,ierror)
call checkmpi(ierror,'Error gathering opt_dwcfg.')
enddo
endif
if(xyzzyaaba1.and.am_master)then
call time_report('VP communication: gathering derivatives',.true.)
do xyzzyaaac9=1,virtual_nnodes-1
read(5,*,iostat=xyzzyaaae9)xyzzyaaad9
if(xyzzyaaae9/=0)call errstop('COMPUTE_LSF_SLAVE','VP error (jnode).')
if(xyzzyaaad9/=virtual_node)cycle
write(6,*)xyzzyaaaj9
write(6,*)xyzzyaaak9
write(6,*)xyzzyaaal9
write(6,*)xyzzyaaab1
write(6,*)((xyzzyaaan1(xyzzyaaab9,xyzzyaaaa9),xyzzyaaab9=1,xyzzyaaab1)&
&,xyzzyaaaa9=1,xyzzyaaae1)
if(complex_wf)then
write(6,*)((xyzzyaaao1(xyzzyaaab9,xyzzyaaaa9),xyzzyaaab9=1,xyzzyaaab1)&
&,xyzzyaaaa9=1,xyzzyaaae1)
endif
if(xyzzyaaaz1)then
write(6,*)((xyzzyaaap1(xyzzyaaab9,xyzzyaaaa9),xyzzyaaab9=1,xyzzyaaab1)&
&,xyzzyaaaa9=1,xyzzyaaae1)
endif
enddo
call time_report('Done.',report_walltime=.true.,no_newline=.true.)
endif
endif
end subroutine compute_lsf_slave
subroutine opt_monitor_header(nl2sol_iteration)
integer,intent(in) :: nl2sol_iteration
if(xyzzyaaau1/=nl2sol_iteration)then
call wout('Optimization monitor :')
call wout('----------------------')
call wout('Function evaluations : '//trim(i2s(xyzzyaaat1+1)))
call wout('NL2SOL iteration     : '//trim(i2s(nl2sol_iteration)))
call wout()
endif
xyzzyaaau1=nl2sol_iteration
end subroutine opt_monitor_header
subroutine param_monitor(params)
real(dp),intent(in) :: params(xyzzyaaae1)
call wout('Current parameters:')
call labelled_list(params,xyzzyaaai1)
call wout()
end subroutine param_monitor
subroutine lsf_monitor(spread,ecentre,weight_variance,no_filtered,spre&
&ad_previous,ecentre_previous,weight_variance_previous,no_filtered_pre&
&vious)
implicit none
integer,intent(in) :: no_filtered
integer,intent(in),optional :: no_filtered_previous
real(dp),intent(in) :: spread,ecentre,weight_variance
real(dp),intent(in),optional :: spread_previous,ecentre_previous,weigh&
&t_variance_previous
logical,save :: xyzzyaaaa12=.true.
character(80) ename,eunits
character(80),save :: xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12
if(xyzzyaaaa12)then
if(use_altsamp)then
xyzzyaaab12='Altsamp optimate (a.u.)'
ename='Altsamp energy centre'
elseif(vm_reweight)then
xyzzyaaab12='Reweighted variance of energy (a.u.)'
ename='Rew. mean energy'
elseif(vm_filter)then
xyzzyaaab12='Filt. unrew. var. of energy (a.u.)'
ename='Filt. unrew. mean energy'
elseif(vm_madmin)then
xyzzyaaab12='Mean abs. dev. from median energy (a.u.)'
ename='Median energy'
else
xyzzyaaab12='Unreweighted variance of energy (a.u.)'
ename='Unreweighted mean energy'
endif
if(model_system)then
eunits='(a.u./particle)'
elseif(isperiodic)then
eunits='(a.u./unit cell)'
else
eunits='(a.u.)'
endif
xyzzyaaac12=trim(ename)//' '//trim(eunits)
xyzzyaaab12=repeat(' ',max(0,41-len_trim(xyzzyaaab12)))//trim(xyzzyaaa&
&b12)//'  :'
xyzzyaaac12=repeat(' ',max(0,41-len_trim(xyzzyaaac12)))//trim(xyzzyaaa&
&c12)//'  :'
if(xyzzyaaaz1)then
xyzzyaaad12='Variance of weights'
xyzzyaaad12=repeat(' ',max(0,41-len_trim(xyzzyaaad12)))//trim(xyzzyaaa&
&d12)//'  :'
endif
if(vm_filter)then
xyzzyaaae12='Configs filtered/total'
xyzzyaaae12=repeat(' ',max(0,41-len_trim(xyzzyaaae12)))//trim(xyzzyaaa&
&e12)//'  :'
endif
xyzzyaaaa12=.false.
endif
if(present(spread_previous))then
call wout(trim(xyzzyaaab12)//' ',(/spread_previous,spread/),rfmt='(es1&
&4.6)',rsep=' -> ')
call wout(trim(xyzzyaaac12)//' ',(/ecentre_previous*xyzzyaaax1,ecentre&
&*xyzzyaaax1/),rfmt='(es14.6)',rsep=' -> ')
if(xyzzyaaaz1)call wout(trim(xyzzyaaad12)//' ',(/weight_variance_previ&
&ous,weight_variance/),rfmt='(es14.6)',rsep=' -> ')
if(vm_filter)call wout(trim(xyzzyaaae12)//'   '//trim(i2s(no_filtered_&
&previous))//' -> '//trim(i2s(no_filtered))//' / '//trim(i2s(xyzzyaaac&
&1)))
else
call wout(trim(xyzzyaaab12)//' ',spread,rfmt='(es14.6)')
call wout(trim(xyzzyaaac12)//' ',ecentre*xyzzyaaax1,rfmt='(es14.6)')
if(xyzzyaaaz1)call wout(trim(xyzzyaaad12)//' ',weight_variance,rfmt='(&
&es14.6)')
if(vm_filter)call wout(trim(xyzzyaaae12)//'   '//trim(i2s(no_filtered)&
&)//' / '//trim(i2s(xyzzyaaac1)))
endif
call wout()
end subroutine lsf_monitor
subroutine rejection_monitor(spread,ecentre,spread_previous,ecentre_pr&
&evious)
implicit none
real(dp),intent(in) :: spread,ecentre,spread_previous,ecentre_previous
call wout('Rejected step:')
if(vm_madmin)then
call wout('MAD       : ',(/spread_previous,spread/),rfmt='(es14.6)',rs&
&ep=' -> ')
else
call wout('Variance  : ',(/spread_previous,spread/),rfmt='(es14.6)',rs&
&ep=' -> ')
endif
call wout('Energy    : ',(/ecentre_previous*xyzzyaaax1,ecentre*xyzzyaa&
&ax1/),rfmt='(es14.6)',rsep=' -> ')
call wout()
end subroutine rejection_monitor
subroutine deriv_monitor(params,params_h,spread_array,spread,ecentre_a&
&rray,ecentre)
implicit none
real(dp),intent(in) :: params(xyzzyaaae1),params_h(xyzzyaaae1),spread_&
&array(xyzzyaaae1),ecentre_array(xyzzyaaae1),spread,ecentre
integer xyzzyaaaa14
real(dp) xyzzyaaab14(xyzzyaaae1),xyzzyaaac14(xyzzyaaae1),xyzzyaaad14
do xyzzyaaaa14=1,xyzzyaaae1
if(params_h(xyzzyaaaa14)==params(xyzzyaaaa14))then
xyzzyaaab14(xyzzyaaaa14)=0.d0
xyzzyaaac14(xyzzyaaaa14)=0.d0
cycle
endif
xyzzyaaad14=abs(params_h(xyzzyaaaa14)+params(xyzzyaaaa14))*0.5d0/(para&
&ms_h(xyzzyaaaa14)-params(xyzzyaaaa14))
xyzzyaaab14(xyzzyaaaa14)=(spread_array(xyzzyaaaa14)-spread)*xyzzyaaad1&
&4
xyzzyaaac14(xyzzyaaaa14)=(ecentre_array(xyzzyaaaa14)-ecentre)*xyzzyaaa&
&d14*xyzzyaaax1
enddo
call wout('Shifted parameter values (internal) for numerical gradient:&
&')
call labelled_list(params_h,xyzzyaaai1)
call wout()
call wout('Rescaled gradient of target function:')
call labelled_list(xyzzyaaab14,xyzzyaaai1)
call wout()
call wout('Rescaled gradient of central energy:')
call labelled_list(xyzzyaaac14,xyzzyaaai1)
call wout()
end subroutine deriv_monitor
subroutine compute_config_energy(icfg,e,ei,w,using_initial_params,isna&
&n,isinf,iszero)
implicit none
integer,intent(in) :: icfg
real(dp),intent(out) :: e,ei,w
logical,intent(in) :: using_initial_params
logical,intent(out) :: isnan,isinf,iszero
complex(dp) xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15
call load_from_storage(xyzzyaaas1,icfg)
call load_from_storage_energy(xyzzyaaas1,icfg)
if(xyzzyaabd1)call mc_twist_offset(twist_config(1:3,icfg))
call eval_local_energy(xyzzyaaas1,etot=e,etot_imag=ei,fix_nl_grid=opt_&
&fixnl,isnan=isnan,isinf=isinf)
if(xyzzyaaaz1)then
call wfn_logval(xyzzyaaas1,xyzzyaaaa15,iszero)
if(using_initial_params)then
if(iszero)call errstop('COMPUTE_CONFIG_ENERGY','Floating-point excepti&
&on reported when calculating wave function value. Don''t know what to&
& do, so stopping. Please file a bug report.')
if(use_altsamp)then
xyzzyaaar1(icfg)=cmplx(logp_config(icfg),0.d0,dp)
else
xyzzyaaar1(icfg)=xyzzyaaaa15
endif
endif
xyzzyaaac15=xyzzyaaaa15-xyzzyaaar1(icfg)
if(4*dble(xyzzyaaac15)>maximum_exp_arg.or.4*dble(xyzzyaaac15)>maximum_&
&exp_arg)then
w=sqrt(max_exp)
else
xyzzyaaab15=exp(xyzzyaaac15)
w=dble(xyzzyaaab15)**2+aimag(xyzzyaaab15)**2
endif
else
w=1.d0
iszero=.false.
endif
call save_to_storage(xyzzyaaas1,icfg)
call save_to_storage_energy(xyzzyaaas1,icfg)
end subroutine compute_config_energy
subroutine compute_config_derivs(icfg,params,params_h,e_h,w_h,isnan,is&
&inf,iszero,ei_h)
implicit none
integer,intent(in) :: icfg
real(dp),intent(in) :: params(xyzzyaaae1),params_h(xyzzyaaae1)
real(dp),intent(out) :: e_h(xyzzyaaae1),w_h(xyzzyaaae1)
real(dp),intent(out),optional :: ei_h(xyzzyaaae1)
logical,intent(out) :: isnan,isinf,iszero
integer xyzzyaaaa16
real(dp) xyzzyaaab16,ei0
complex(dp) xyzzyaaac16,logwfn,xyzzyaaad16,xyzzyaaae16
logical nan,inf,zero
isnan=.false.
isinf=.false.
iszero=.false.
call load_from_storage(xyzzyaaas1,icfg)
call load_from_storage_energy(xyzzyaaas1,icfg)
if(.not.complex_wf)then
call prepare_ederivs(xyzzyaaas1,xyzzyaaah1,xyzzyaaab16,xyzzyaaac16,get&
&_logwfn=xyzzyaabc1)
else
call prepare_ederivs(xyzzyaaas1,xyzzyaaah1,xyzzyaaab16,xyzzyaaac16,get&
&_logwfn=xyzzyaabc1,ei0=ei0)
endif
call save_to_storage(xyzzyaaas1,icfg)
call save_to_storage_energy(xyzzyaaas1,icfg)
do xyzzyaaaa16=1,xyzzyaaae1
if(.not.complex_wf)then
call eval_energy_nderiv(xyzzyaaas1,params,xyzzyaaaa16,params(xyzzyaaaa&
&16),params_h(xyzzyaaaa16),xyzzyaaah1,xyzzyaaab16,xyzzyaaac16,e_h(xyzz&
&yaaaa16),logwfn,icfg==1,get_logwfn=xyzzyaabc1,isnan=nan,isinf=inf,isz&
&ero=zero)
else
call eval_energy_nderiv(xyzzyaaas1,params,xyzzyaaaa16,params(xyzzyaaaa&
&16),params_h(xyzzyaaaa16),xyzzyaaah1,xyzzyaaab16,xyzzyaaac16,e_h(xyzz&
&yaaaa16),logwfn,icfg==1,get_logwfn=xyzzyaabc1,ei0=ei0,ei1=ei_h(xyzzya&
&aaa16),isnan=nan,isinf=inf,iszero=zero)
endif
isnan=isnan.or.nan
isinf=isinf.or.inf
iszero=iszero.or.zero
if(xyzzyaaaz1)then
xyzzyaaae16=logwfn-xyzzyaaar1(icfg)
if(4*dble(xyzzyaaae16)>maximum_exp_arg.or.4*dble(xyzzyaaae16)>maximum_&
&exp_arg)then
w_h(xyzzyaaaa16)=sqrt(max_exp)
else
xyzzyaaad16=exp(xyzzyaaae16)
w_h(xyzzyaaaa16)=dble(xyzzyaaad16)**2+aimag(xyzzyaaad16)**2
endif
endif
enddo
end subroutine compute_config_derivs
end module slaarnacm
subroutine madr(n,p,x,nf,r,nl2sol_iteration,been_rst,stop_opt,x_h,bad_&
&point)
use machine_constants
use slaarnacm, only : compute_lsf_master
implicit none
integer,intent(in) :: n,p,nl2sol_iteration,been_rst
integer,intent(inout) :: nf
real(kind=kind(0.d0)),intent(in) :: x(:)
real(kind=kind(0.d0)),intent(in),optional :: x_h(:)
real(kind=kind(0.d0)),intent(out) :: r(:)
logical,intent(inout) :: stop_opt
logical,intent(inout),optional :: bad_point
logical xyzzyaaaa17
if(.not.present(x_h))then
call compute_lsf_master(x,r,nl2sol_iteration,been_rst,stop_opt,xyzzyaa&
&aa17)
else
call compute_lsf_master(x,r,nl2sol_iteration,been_rst,stop_opt,xyzzyaa&
&aa17,params_h=x_h)
endif
if(present(bad_point))bad_point=xyzzyaaaa17
end subroutine madr

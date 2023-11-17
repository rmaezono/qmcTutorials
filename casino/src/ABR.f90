subroutine monte_carlo
use slaarnaag
use dsp
use slaarnaaq
use slaarnabg
use slaarnabk
use slaarnabp
!$ use openmp_base
use parallel
use slaarnaca
use store
use slaarnacm_linjas
use slaarnacq
use slaarnacr
use slaarnaab,    only : forces_info,initialize_forces
use slaarnaac,           only : readawf,awfdet_setup
use slaarnaad,         only : init_backflow,setup_backflow_plot,use_gb&
&ackflow
use slaarnaae,           only : single_precision_blips,readbwf,bwfdet_&
&setup,deshalloc_bwfdet_shm,write_binary_blips,conv_binary_blips
use casl,             only : read_casl,write_casl,query_casl_item,firs&
&t_unread_child,casl_keysize
use slaarnaaf,          only : copy_config_file,shift_config_files,con&
&_in,con_out,con_backup,con_loc,checkpoint_ncpu,make_chkpoint_groups,d&
&ismantle_configs
use slaarnaah,      only : read_correlation_header,write_correlation_h&
&eader,jastrow_in_corr
use slaarnaaj,              only : dmc_main,limdmc,tpdmc,max_rec_attem&
&pts,ebest_av_window,corper_dmc,ndmcave,cerefdmc,wdmcmin,wdmcmax,targ_&
&wt,trip_popn,dmc_init_eref,alimit,iaccum,ibran,lwdmc,dmc_equil_fixpop&
&,lwdmc_fixpop,growth_estimator,writeout_dmc_hist,dmc_twist_av,use_ini&
&t_eref,nc_dmc,poprenorm,nucleus_gf_mods,dmc_reweight_configs,dmc_spac&
&ewarping,nconfig_prelim,redist_grp_size
use slaarnaak,           only : dwfdet_setup
use esdf,             only : block_data,esdf_block,esdf_init,esdf_clos&
&e,esdf_integer,esdf_boolean,esdf_double,esdf_string,esdf_physical,esd&
&f_defined,esdf_warnout,esdf_dump,llength
use slaarnaal,             only : emin_main,emin_xi_value,emin_min_ene&
&rgy,emin_mine_present,emin_auto_varmin,emin_auto_varmin_present
use slaarnaan,    only : netot_nitot_products,init_geometry,print_geom&
&etry,points
use slaarnaao,            only : read_exmol
use slaarnaap,            only : read_expot,setup_expot
use slaarnabe,          only : read_geminal
use slaarnaar,          only : init_expot_wfdet
use file_utils,       only : open_units
use format_utils,     only : wout,i2s,r2s,r2ss,l2s,wordwrap,allow_slav&
&e_write,capitalize,switch_case
use slaarnaas,        only : init_free_orbs,corr_heg_required,heg_nbas&
&is,r_s,heg_cell,heg_wigner_basis,heg_crystal_type,heg_ferro,heg_nosit&
&es,heg_crystal_sites,heg_repeat,calc_hf_energies,mc_twist_offset,me_b&
&iex3,mh_biex3,xx_sep,harmwire_b,free_norb,eval_finite_hf_energies
use slaarnaat,        only : gautol,printgscreening,cusp_correction,cu&
&sp_control,cusp_info,molgscreening,cusp_threshold,s_plot
use slaarnaau,        only : readgw,gwfdet_setup,deshalloc_gauss_shm
use slaarnabf,only : run_mpc_generation
use slaarnabi,             only : use_gpcc
use slaarnabj,          only : check_hist_header
use slaarnabl,          only : read_jastrow,setup_jastrow_plot,use_gja&
&strow,finite_size_corr_ke_jastrow,check_varmin_linjas_jastrow,gen_gja&
&strow
use slaarnabm,         only : use_magnetic_field,read_magnetic_field
use slaarnabn,         only : read_mahan
use slaarnabs,        only : init_non_local,print_non_local_grid
use slaarnabt,        only : overflow_protection,inverse3,hypothenuse
use slaarnaby,             only : init_plot,setup_plot,plot_main
use slaarnabyter,          only : init_plotter,setup_plotter,qmc_plott&
&er
use slaarnacb,           only : readpwf,pwfdet_setup
use slaarnacc,   only : initialize_random,random_seed_kw
use slaarnace,       only : relativistic,relativity_setup
use slaarnacf,              only : rmc_main
use run_control,      only : errstop,errstop2,errstop_master,errwarn,t&
&imer,check_alloc,errwarn_silent,tcputime
use slaarnaci,              only : sdw_magvec
use shalloc,          only : shm_size,need_shm,shm_size_nproc
use slaarnacj,           only : dbarrc,small_transfer,small_buffers
use slaarnack,  only : setup_special_wfn
use slaarnacl,         only : stowfdet_read,stowfdet_setup
use slaarnacm,           only : varmin_main,vm_reweight,vm_forgiving,v&
&m_w_min,vm_w_max,vm_use_e_guess,vm_e_guess,vm_filter,vm_filter_thres,&
&vm_filter_width,vm_madmin
use slaarnaco,             only : compute_ecppii,compute_ii_fields
use slaarnacp,              only : vmc_main,vmc_ionjump,corper_default&
&_vmc,corper_default_opt
use slaarnacs,        only : psi_s,init_wfn,update_wfn_casl
implicit none
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae1,xyzzyaa&
&af1,xyzzyaaag1,xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1,xyzzyaaak1,xyzzyaaal1&
&,xyzzyaaam1,xyzzyaaan1,xyzzyaaao1,xyzzyaaap1,xyzzyaaaq1,xyzzyaaar1,xy&
&zzyaaas1,xyzzyaaat1,xyzzyaaau1,xyzzyaaav1,xyzzyaaaw1,xyzzyaaax1,xyzzy&
&aaay1,xyzzyaaaz1,xyzzyaaba1,xyzzyaabb1,xyzzyaabc1,xyzzyaabd1,xyzzyaab&
&e1,xyzzyaabf1,xyzzyaabg1,xyzzyaabh1,xyzzyaabi1,xyzzyaabj1,xyzzyaabk1,&
&xyzzyaabl1,xyzzyaabm1,xyzzyaabn1,xyzzyaabo1,xyzzyaabp1,xyzzyaabq1,xyz&
&zyaabr1,xyzzyaabs1,xyzzyaabt1,xyzzyaabu1,xyzzyaabv1,xyzzyaabw1,xyzzya&
&abx1,xyzzyaaby1,xyzzyaabz1,xyzzyaaca1,xyzzyaacb1,xyzzyaacc1
real(dp) xyzzyaacd1,xyzzyaace1,xyzzyaacf1,xyzzyaacg1,xyzzyaach1,xyzzya&
&aci1,xyzzyaacj1,xyzzyaack1,xyzzyaacl1
logical xyzzyaacm1,xyzzyaacn1,xyzzyaaco1,xyzzyaacp1,xyzzyaacq1,xyzzyaa&
&cr1,xyzzyaacs1,xyzzyaact1,xyzzyaacu1,writeout_vmc_hist,xyzzyaacv1,xyz&
&zyaacw1,xyzzyaacx1,xyzzyaacy1
character(20) opt_method,parallel_keywords,vmc_sampling
character(20),allocatable :: xyzzyaacz1(:)
logical, allocatable :: xyzzyaada1(:),xyzzyaadb1(:),xyzzyaadc1(:),xyzz&
&yaadd1(:),xyzzyaade1(:),xyzzyaadf1(:)
integer, allocatable :: xyzzyaadg1(:)
integer xyzzyaadh1,xyzzyaadi1,xyzzyaadj1
character(casl_keysize) casl_key
character(512) errmsg
integer xyzzyaadk1,xyzzyaadl1,xyzzyaadm1,xyzzyaadn1,xyzzyaado1(3),xyzz&
&yaadp1(5),xyzzyaadq1(1),xyzzyaadr1,xyzzyaads1,xyzzyaadt1,xyzzyaadu1,x&
&yzzyaadv1,xyzzyaadw1,nchildren
integer,allocatable :: xyzzyaadx1(:,:),xyzzyaady1(:,:),xyzzyaadz1(:,:,&
&:),xyzzyaaea1(:,:,:,:),xyzzyaaeb1(:,:,:,:),xyzzyaaec1(:)
real(sp) xyzzyaaed1,xyzzyaaee1
real(dp) xyzzyaaef1
real(dp),allocatable :: xyzzyaaeg1(:),xyzzyaaeh1(:),xyzzyaaei1(:),xyzz&
&yaaej1(:),xyzzyaaek1(:,:)
logical xyzzyaael1,xyzzyaaem1,xyzzyaaen1,xyzzyaaeo1,xyzzyaaep1,xyzzyaa&
&eq1,xyzzyaaer1,xyzzyaaes1,xyzzyaaet1
logical,allocatable :: xyzzyaaeu1(:)
call timer('SETUP',.true.)
xyzzyaaed1=tcputime()
call xyzzyaaev1
call overflow_protection
call esdf_init('input')
if(am_master)call esdf_warnout
call xyzzyaaew1
call xyzzyaafc1
call xyzzyaafd1
if(am_master)then
call xyzzyaaff1
call xyzzyaaey1
endif
call esdf_close
ranlog_unit=0
if(am_master.and.ranprint>0)call open_units(ranlog_unit,xyzzyaadl1)
call initialize_random(ranluxlevel,ranprint,ranlog_unit)
call make_chkpoint_groups
call read_correlation_header(.true.)
call read_casl('parameters.casl',errmsg)
if(len_trim(errmsg)>0)call errstop_master('MONTE_CARLO',trim(errmsg))
call query_casl_item('parameters.casl',nchildren=nchildren)
if(nchildren>0)then
if(am_master)then
call wout('Reading parameters.casl')
call wout('=======================')
call wout('Contents of parameters.casl:')
xyzzyaadl1=0
do while(xyzzyaadl1==0)
call first_unread_child('parameters.casl',casl_key,xyzzyaadl1,flag_as_&
&read=.true.)
if(xyzzyaadl1==0)then
call capitalize(casl_key)
call wout(' * '//trim(adjustl(casl_key)))
endif
enddo
call wout()
endif
endif
call read_mdet
call read_wave_function
if(am_master)call check_wave_function
call xyzzyaafg1
call print_determinant_info
call xyzzyaafs1
call xyzzyaafo1
call xyzzyaafn1
call setup_wfdet
if(xyzzyaacy1)then
if(am_master)call wfdet_write_mdet_casl
return
endif
call xyzzyaafr1
if(isgen_mpc)then
call timer('SETUP',.false.)
if(xyzzyaacm1)then
if(am_master)call wout('This is a TEST RUN -- skipping MPC generation.&
&')
else
call run_mpc_generation(xyzzyaacf1)
endif
return
endif
if(use_expot)then
call read_expot
call setup_expot
endif
if(use_magnetic_field)call read_magnetic_field
if(.not.xyzzyaaen1)call xyzzyaafq1
call xyzzyaafp1
call xyzzyaafu1
if(relativistic)call relativity_setup(xyzzyaacg1,isdmc.or.isvmc_dmc.or&
&.isdmc_dmc.or.isvmc_dmc_equil)
select case(trim(psi_s))
case('slater')
continue
case('exmol')
call read_exmol
case('pfaffian')
continue
case('geminal')
call read_geminal
case('mahan_ex')
call read_mahan
end select
call read_jastrow_function
call xyzzyaaft1
if(have_ppots)then
call init_non_local(nlrule1_ion)
if(am_master)call print_non_local_grid
endif
call xyzzyaafv1
if(forces)call initialize_forces
if(am_master.and.need_shm)call shm_size
if(am_master)then
call wout()
call wout('Setup complete.')
call wout()
xyzzyaaee1=tcputime()-xyzzyaaed1
call wout('Time taken in setup    : : : ',real(xyzzyaaee1,dp),rfmt='(f&
&13.4)')
call wout()
endif
if(xyzzyaacm1)then
if(am_master)then
call wout('TEST RUN only.')
call wout('Quitting.')
endif
elseif(xyzzyaaen1)then
if(xyzzyaadw1==1)then
call setup_plot
call plot_main(xyzzyaaaa1,xyzzyaaei1,xyzzyaaej1,xyzzyaaec1,xyzzyaaek1,&
&xyzzyaaeu1)
else
call setup_plotter
call qmc_plotter(xyzzyaaaa1,xyzzyaaei1,xyzzyaaej1,xyzzyaaec1,xyzzyaaek&
&1,xyzzyaaeu1)
endif
elseif(pair_corr)then
call xyzzyaafw1
elseif(isopt)then
opt_cycle=-1
opt_method=xyzzyaacz1(1)
opt_jastrow=xyzzyaada1(1)
opt_backflow=xyzzyaadb1(1)
opt_det_coeff=xyzzyaadc1(1)
opt_orbitals=xyzzyaadd1(1)
opt_geminal=xyzzyaade1(1)
opt_maxiter=xyzzyaadg1(1)
fix_cutoffs=xyzzyaadf1(1)
if(opt_jastrow)use_jastrow=.true.
if(opt_backflow)use_backflow=.true.
call init_wfn
call xyzzyaagd1('PERFORMING A SINGLE OPTIMIZATION RUN.')
call xyzzyaafh1(xyzzyaaeo1)
if(xyzzyaaeo1)return
elseif(isvmc)then
call xyzzyaagd1('PERFORMING A SINGLE VMC CALCULATION.')
call xyzzyaafj1(xyzzyaaeo1)
elseif(isdmc)then
if(.not.iaccum)then
call xyzzyaagd1('PERFORMING A SINGLE DMC EQUILIBRATION CALCULATION.')
else
call xyzzyaagd1('PERFORMING A SINGLE DMC STATISTICS-ACCUMULATION CALCU&
&LATION.')
endif
if(xyzzyaaav1==0.or..not.iaccum)then
call xyzzyaafk1(xyzzyaaeo1)
else
call xyzzyaafm1(.true.,xyzzyaaeo1)
endif
elseif(isvmc_dmc)then
call xyzzyaagd1('PERFORMING A VMC CONFIGURATION-GENERATION CALCULATION&
&.')
xyzzyaaem1=forces
forces=.false.
xyzzyaael1=expvals
expvals=.false.
call xyzzyaafj1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call xyzzyaagc1
call xyzzyaagd1('PERFORMING A DMC EQUILIBRATION CALCULATION.')
iaccum=.false.
newrun=.true.
call shift_config_files
forces=xyzzyaaem1
call xyzzyaafk1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call xyzzyaagc1
call xyzzyaagd1('PERFORMING A DMC STATISTICS-ACCUMULATION CALCULATION.&
&')
iaccum=.true.
newrun=.false.
expvals=xyzzyaael1
call shift_config_files
if(xyzzyaaav1==0)then
call xyzzyaafk1(xyzzyaaeo1)
else
call xyzzyaafm1(.false.,xyzzyaaeo1)
endif
elseif(isvmc_dmc_equil)then
call xyzzyaagd1('PERFORMING A VMC CONFIGURATION-GENERATION CALCULATION&
&.')
xyzzyaaem1=forces
forces=.false.
expvals=.false.
call xyzzyaafj1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call xyzzyaagc1
call xyzzyaagd1('PERFORMING A DMC EQUILIBRATION CALCULATION.')
iaccum=.false.
newrun=.true.
call shift_config_files
forces=xyzzyaaem1
call xyzzyaafk1(xyzzyaaeo1)
elseif(isdmc_dmc)then
call xyzzyaagd1('PERFORMING A DMC EQUILIBRATION CALCULATION.')
iaccum=.false.
xyzzyaael1=expvals
expvals=.false.
call xyzzyaafk1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call xyzzyaagc1
call xyzzyaagd1('PERFORMING A DMC STATISTICS-ACCUMULATION CALCULATION.&
&')
iaccum=.true.
newrun=.false.
expvals=xyzzyaael1
call shift_config_files
if(xyzzyaaav1==0)then
call xyzzyaafk1(xyzzyaaeo1)
else
call xyzzyaafm1(.false.,xyzzyaaeo1)
endif
elseif(isrmc)then
if(.not.iaccum)then
call xyzzyaagd1('PERFORMING A SINGLE RMC EQUILIBRATION CALCULATION.')
else
call xyzzyaagd1('PERFORMING A SINGLE RMC STATISTICS-ACCUMULATION CALCU&
&LATION.')
endif
call xyzzyaafl1(xyzzyaaeo1)
elseif(isrmc_rmc)then
call xyzzyaagd1('PERFORMING AN RMC EQUILIBRATION CALCULATION.')
iaccum=.false.
xyzzyaael1=expvals
expvals=.false.
call xyzzyaafl1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call xyzzyaagc1
call xyzzyaagd1('PERFORMING AN RMC STATISTICS-ACCUMULATION CALCULATION&
&.')
iaccum=.true.
newrun=.false.
expvals=xyzzyaael1
call shift_config_files
call xyzzyaafl1(xyzzyaaeo1)
elseif(isvmc_opt.or.isopt_vmc)then
expvals=.false.
forces=.false.
allocate(xyzzyaaeg1(opt_cycles),xyzzyaaeh1(opt_cycles),stat=xyzzyaadk1&
&)
call check_alloc(xyzzyaadk1,'MONTE_CARLO','vmc_E')
xyzzyaaeg1=0.d0
xyzzyaaeh1=0.d0
if(isopt_vmc)xyzzyaaeh1(1)=-1.d0
xyzzyaaes1=use_jastrow
xyzzyaaet1=use_backflow
if(any(xyzzyaada1))use_jastrow=.true.
if(any(xyzzyaadb1))use_backflow=.true.
opt_cycle=0
call write_correlation_header(.true.,.true.)
call update_wfn_casl
if(am_master)then
call write_casl(':parameters.casl','parameters.0.casl',errmsg)
if(len_trim(errmsg)>0)call errstop_master('MONTE_CARLO',trim(errmsg))
endif
use_jastrow=xyzzyaaes1
use_backflow=xyzzyaaet1
do opt_cycle=1,opt_cycles
opt_method=xyzzyaacz1(opt_cycle)
opt_jastrow=xyzzyaada1(opt_cycle)
opt_backflow=xyzzyaadb1(opt_cycle)
opt_det_coeff=xyzzyaadc1(opt_cycle)
opt_orbitals=xyzzyaadd1(opt_cycle)
opt_geminal=xyzzyaade1(opt_cycle)
opt_maxiter=xyzzyaadg1(opt_cycle)
fix_cutoffs=xyzzyaadf1(opt_cycle)
if(xyzzyaadh1==opt_cycle)then
use_jastrow=.true.
call init_wfn
endif
if(opt_cycle>1)then
call shift_config_files
call xyzzyaagc1
endif
if(isvmc_opt.or.opt_cycle>1)then
call xyzzyaagd1('PERFORMING VMC CONFIGURATION-GENERATION CALCULATION N&
&o. '//trim(i2s(opt_cycle)))
call xyzzyaafj1(xyzzyaaeo1,xyzzyaaeg1(opt_cycle),xyzzyaaeh1(opt_cycle)&
&)
if(xyzzyaaeo1)return
if(am_master)call xyzzyaafi1
call xyzzyaagc1
endif
call xyzzyaagd1('PERFORMING OPTIMIZATION CALCULATION No. '//trim(i2s(o&
&pt_cycle)))
if(xyzzyaadi1==opt_cycle)use_jastrow=.true.
if(xyzzyaadj1==opt_cycle)use_backflow=.true.
if(xyzzyaadi1==opt_cycle.or.xyzzyaadj1==opt_cycle)call init_wfn
if(am_master.and.any(xyzzyaada1).and.use_gjastrow.and..not.use_jastrow&
&)call write_casl(':parameters.casl','parameters.'//trim(i2s(opt_cycle&
&))//'.casl',errmsg)
if(isvmc_opt.or.opt_cycle>1)call shift_config_files
call xyzzyaafh1(xyzzyaaeo1)
if(xyzzyaaeo1)return
if(use_jastrow.and.finite_size_corr)call finite_size_corr_ke_jastrow
enddo
deallocate(xyzzyaaeg1,xyzzyaaeh1)
if(xyzzyaacr1)then
call xyzzyaagc1
call xyzzyaagd1('PERFORMING POST-FIT VMC CALCULATION.')
if(.not.xyzzyaacs1)xyzzyaaae1=0
opt_cycle=opt_cycles+1
call shift_config_files
call xyzzyaafj1(xyzzyaaeo1)
endif
else
call errstop_master('MONTE_CARLO','RUNTYPE problem - this is a bug.')
endif
call deshalloc_bwfdet_shm
call deshalloc_gauss_shm
call timer('SETUP',.false.)
contains
subroutine xyzzyaaev1
implicit none
calc_field=.false.
constant_energy=0.d0
corr_heg_required=.false.
dimensionality=3
fix_cutoffs=.false.
have_biex3pot=.false.
have_jastrow3=.false.
have_ppots=.false.
iaccum=.false.
ibran=.false.
init_by_ion=.false.
init_by_iontype=.false.
inversion_symmetry=.false.
isitcomplex=.true.
nbasis=0
nitype=0
noncoll_spin=.false.
npcells=0
scell_matrix=0
open_unit=.false.
end subroutine xyzzyaaev1
subroutine xyzzyaaew1
implicit none
integer xyzzyaaaa3,xyzzyaaab3,idum1,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,x&
&yzzyaaaf3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,nlines,xyzzyaaa&
&k3,xyzzyaaal3,xyzzyaaam3,xyzzyaaan3
real(dp) xyzzyaaao3,xyzzyaaap3(3),xyzzyaaaq3(3),xyzzyaaar3(3),xyzzyaaa&
&s3,xyzzyaaat3,xyzzyaaau3(3)
logical xyzzyaaav3
character(llength) token
integer xyzzyaaaw3,xyzzyaaax3,xyzzyaaay3
logical xyzzyaaaz3,xyzzyaaba3
character(10) special_wfn
xyzzyaaaf1=esdf_integer('neu',-1)
xyzzyaaag1=esdf_integer('ned',-1)
xyzzyaaah1=esdf_integer('nhu',0)
xyzzyaaai1=esdf_integer('nhd',0)
atom_basis_type=esdf_string('atom_basis_type','default')
isperiodic=esdf_boolean('periodic',.true.)
xyzzyaaar1=esdf_integer('blip_periodicity',-1)
complex_wf=esdf_boolean('complex_wf',.false.)
psi_s=esdf_string('psi_s','slater')
runtype=esdf_string('runtype','default')
newrun=esdf_boolean('newrun',.true.)
xyzzyaacm1=esdf_boolean('testrun',.false.)
block_time=real(esdf_physical('block_time',0.d0,'s'))
use_blocktime=block_time>0.
xyzzyaabd1=esdf_integer('vmc_equil_nstep',5000)
xyzzyaaay1=esdf_integer('vmc_nstep',-1)
xyzzyaabc1=esdf_integer('vmc_nblock',1)
xyzzyaaaz1=esdf_integer('vmc_nconfig_write',0)
if(trim(runtype)=='vmc')then
xyzzyaabb1=esdf_integer('vmc_decorr_period',corper_default_vmc)
else
xyzzyaabb1=esdf_integer('vmc_decorr_period',corper_default_opt)
endif
xyzzyaaba1=esdf_integer('vmc_ave_period',1)
xyzzyaacl1=esdf_double('dtvmc_shift',-1.d0)
xyzzyaaaq1=esdf_integer('opt_dtvmc',1)
writeout_vmc_hist=esdf_boolean('writeout_vmc_hist',.true.)
xyzzyaabn1=esdf_integer('vmc_ntwist',0)
xyzzyaabk1=esdf_integer('vmc_reequil_nstep',500)
vmc_ionjump=esdf_double('vmc_ionjump',0.d0)
xyzzyaaao1=esdf_integer('vmc_method',1)
if(xyzzyaaao1==3)then
xyzzyaacd1=esdf_double('dtvmc',0.01d0)
else
xyzzyaacd1=esdf_double('dtvmc',0.1d0)
endif
vmc_sampling=esdf_string('vmc_sampling','standard')
vmc_optimum_e0=esdf_physical('vmc_optimum_e0',0.d0,'hartree')
vmc_optimum_ew=esdf_physical('vmc_optimum_ew',100.d0,'hartree')
xyzzyaaaa1=esdf_integer('nequil',5000)
xyzzyaaab1=esdf_integer('nmove',-1)
xyzzyaaac1=esdf_integer('nblock',1)
if(trim(runtype)=='vmc')then
xyzzyaaad1=esdf_integer('corper',corper_default_vmc)
else
xyzzyaaad1=esdf_integer('corper',corper_default_opt)
corper_default_opt=xyzzyaaad1
endif
xyzzyaaaj1=esdf_integer('nvmcave',1)
xyzzyaaae1=esdf_integer('nwrcon',0)
xyzzyaacv1=esdf_boolean('vmc_twist_av',.false.)
xyzzyaaaw1=esdf_integer('nequil_ta',500)
xyzzyaabe1=esdf_integer('dmc_equil_nstep',-1)
xyzzyaabf1=esdf_integer('dmc_equil_nblock',1)
xyzzyaabg1=esdf_integer('dmc_stats_nstep',-1)
xyzzyaabh1=esdf_integer('dmc_stats_nblock',1)
dmc_md=esdf_boolean('dmc_md',.false.)
xyzzyaabz1=esdf_integer('dmcmd_equil_nstep',-1)
xyzzyaaca1=esdf_integer('dmcmd_stats_nstep',-1)
xyzzyaaci1=esdf_double('dmc_target_weight',0.d0)
xyzzyaace1=esdf_double('dtdmc',0.01d0)
xyzzyaabi1=esdf_integer('dmc_ave_period',1)
xyzzyaabj1=esdf_integer('dmc_decorr_period',1)
xyzzyaacj1=esdf_double('dmc_trip_weight',0.d0)
max_rec_attempts=esdf_integer('max_rec_attempts',5)
xyzzyaaby1=esdf_integer('dmc_nconf_prelim',-1)
lwdmc=esdf_boolean('lwdmc',.false.)
lwdmc_fixpop=esdf_boolean('lwdmc_fixpop',.false.)
wdmcmin=esdf_double('wdmcmin',0.5d0)
wdmcmax=esdf_double('wdmcmax',2.d0)
use_tmove=esdf_boolean('use_tmove',.false.)
xyzzyaabo1=esdf_integer('dmc_ntwist',0)
xyzzyaabl1=esdf_integer('dmc_reequil_nstep',200)
xyzzyaabm1=esdf_integer('dmc_reequil_nblock',1)
use_future=esdf_boolean('future_walking',.false.)
tpdmc=esdf_integer('tpdmc',0)
ibran=esdf_boolean('ibran',.true.)
xyzzyaaap1=esdf_integer('dmc_method',1)
dmc_equil_fixpop=esdf_double('dmc_equil_fixpop',0.d0)
nucleus_gf_mods=esdf_boolean('nucleus_gf_mods',.true.)
xyzzyaact1=esdf_defined('nucleus_gf_mods','L')
dmc_reweight_configs=esdf_boolean('dmc_reweight_conf',.false.)
dmc_spacewarping=esdf_boolean('dmc_spacewarping',.false.)
max_cpu_time=real(esdf_physical('max_cpu_time',0.d0,'s'))
max_real_time=real(esdf_physical('max_real_time',0.d0,'s'))
growth_estimator=esdf_boolean('growth_estimator',.false.)
limdmc=esdf_integer('limdmc',4)
alimit=esdf_double('alimit',0.5d0)
cerefdmc=esdf_double('cerefdmc',1.d0)
ebest_av_window=esdf_integer('ebest_av_window',25)
popstats=esdf_boolean('popstats',.false.)
writeout_dmc_hist=esdf_boolean('writeout_dmc_hist',.true.)
small_transfer=esdf_boolean('small_transfer',.false.)
small_buffers=esdf_boolean('opt_small_buffers',.false.)
dmc_init_eref=esdf_physical('dmc_init_eref',0.d0,'hartree')
use_init_eref=esdf_defined('dmc_init_eref','P')
nc_dmc=esdf_boolean('dmc_norm_conserve',.false.)
poprenorm=esdf_boolean('dmc_poprenorm',.false.)
redist_grp_size=esdf_integer('redist_grp_size',500)
xyzzyaaak1=esdf_integer('nmove_dmc_equil',-1)
xyzzyaaal1=esdf_integer('nmove_dmc_stats',-1)
xyzzyaaam1=esdf_integer('nblock_dmc_equil',-1)
xyzzyaaan1=esdf_integer('nblock_dmc_stats',-1)
xyzzyaacb1=esdf_integer('nmove_dmcmd_equil',-1)
xyzzyaacc1=esdf_integer('nmove_dmcmd_stats',-1)
targ_wt=esdf_double('nconfig',0.d0)
nconfig_prelim=esdf_integer('nconfig_prelim',-1)
trip_popn=esdf_double('trip_popn',0.d0)
corper_dmc=esdf_integer('corper_dmc',1)
ndmcave=esdf_integer('ndmcave',1)
xyzzyaaav1=esdf_integer('num_dmc_twists',0)
xyzzyaaat1=esdf_integer('nmove_dmct_equil',40)
xyzzyaaau1=esdf_integer('nblock_dmct_equil',4)
xyzzyaabs1=esdf_integer('rmc_equil_nstep',-1)
xyzzyaabu1=esdf_integer('rmc_equil_nblock',1)
xyzzyaabt1=esdf_integer('rmc_stats_nstep',-1)
xyzzyaabv1=esdf_integer('rmc_stats_nblock',1)
xyzzyaabq1=esdf_integer('rmc_decorr_period',1)
xyzzyaabr1=esdf_integer('rmc_ave_period',1)
xyzzyaack1=esdf_double('dtrmc',-1.d0)
xyzzyaabw1=esdf_integer('rmc_rep_length',-1)
xyzzyaabx1=esdf_integer('rmc_move_length',1)
xyzzyaacw1=esdf_boolean('rmc_bounce',.true.)
xyzzyaacx1=esdf_boolean('rmc_meas_pos',.false.)
opt_method=esdf_string('opt_method','varmin')
opt_jastrow=esdf_boolean('opt_jastrow',.true.)
opt_det_coeff=esdf_boolean('opt_det_coeff',.false.)
opt_backflow=esdf_boolean('opt_backflow',.false.)
opt_orbitals=esdf_boolean('opt_orbitals',.false.)
opt_geminal=esdf_boolean('opt_geminal',.false.)
opt_cycles=esdf_integer('opt_cycles',4)
opt_maxeval=esdf_integer('opt_maxeval',200)
opt_maxiter=esdf_integer('opt_maxiter',10)
opt_info=esdf_integer('opt_info',2)
opt_strict=esdf_boolean('opt_strict',.false.)
if(trim(opt_method)=='emin'.and.trim(vmc_sampling)=='standard')then
opt_fixnl=esdf_boolean('opt_fixnl',.false.)
else
opt_fixnl=esdf_boolean('opt_fixnl',.true.)
endif
xyzzyaaax1=esdf_integer('opt_noctf_cycles',0)
vm_forgiving=esdf_boolean('vm_forgiving',.true.)
vm_w_min=esdf_double('vm_w_min',0.d0)
vm_w_max=esdf_double('vm_w_max',0.d0)
vm_use_e_guess=esdf_boolean('vm_use_E_guess',.false.)
vm_e_guess=esdf_physical('vm_E_guess',-999.d0,'hartree')
vm_smooth_limits=esdf_boolean('vm_smooth_limits',.true.)
vm_linjas_method=esdf_string('vm_linjas_method','BFGS')
vm_linjas_its=esdf_integer('vm_linjas_its',-1)
vm_filter=esdf_boolean('vm_filter',.false.)
vm_filter_thres=esdf_double('vm_filter_thres',4.d0)
vm_filter_width=esdf_double('vm_filter_width',2.d0)
vm_reweight=esdf_boolean('vm_reweight',.false.)
emin_min_energy=esdf_physical('emin_min_energy',0.d0,'hartree')
emin_mine_present=esdf_defined('emin_min_energy','P')
emin_auto_varmin=esdf_boolean('emin_auto_varmin',.true.)
emin_auto_varmin_present=esdf_defined('emin_auto_varmin','L')
emin_xi_value=esdf_double('emin_xi_value',1.d0)
xyzzyaacr1=esdf_boolean('postfit_vmc',.true.)
xyzzyaacs1=esdf_boolean('postfit_keep_cfg',.false.)
interaction=esdf_string('interaction','default')
use_jastrow=esdf_boolean('use_jastrow',.true.)
xyzzyaaep1=esdf_defined('use_gjastrow','L')
use_gjastrow=esdf_boolean('use_gjastrow',.false.)
gen_gjastrow=esdf_boolean('gen_gjastrow',.false.)
use_backflow=esdf_boolean('backflow',.false.)
use_gbackflow=esdf_boolean('use_gbackflow',.false.)
use_gbackflow=.false.
bf_sparse=esdf_boolean('bf_sparse',.false.)
xyzzyaaeq1=esdf_defined('single_precision_blips','L')
xyzzyaaer1=esdf_defined('sp_blips','L')
if(xyzzyaaeq1.and.xyzzyaaer1)call errstop('SET_INPUT_PARAMETERS','Cann&
&ot have SINGLE_PRECISION_BLIPS and SP_BLIPS in the same input file. T&
&he former is deprecated.')
if(xyzzyaaeq1)single_precision_blips=esdf_boolean('single_precision_bl&
&ips',.false.)
if(xyzzyaaer1)single_precision_blips=esdf_boolean('sp_blips',.false.)
write_binary_blips=esdf_boolean('write_binary_blips',.true.)
conv_binary_blips=esdf_boolean('conv_binary_blips',.false.)
blip_mpc=esdf_boolean('blip_mpc',.false.)
use_orbmods=esdf_boolean('use_orbmods',.false.)
use_expot=esdf_boolean('expot',.false.)
xyzzyaacu1=esdf_boolean('timing_info',.false.)
esupercell=esdf_boolean('esupercell',.false.)
neighprint=esdf_integer('neighprint',0)
xyzzyaach1=esdf_physical('e_offset',0.d0,'hartree')
xyzzyaacf1=esdf_physical('mpc_cutoff',30.d0,'hartree')
gautol=esdf_double('gautol',7.d0)
printgscreening=esdf_boolean('printgscreening',.false.)
dbarrc=esdf_integer('dbarrc',100000)
con_loc=esdf_string('con_loc','.')
xyzzyaaas1=esdf_integer('non_local_grid',-1)
lcutofftol=esdf_double('lcutofftol',0.00001d0)
nlcutofftol=esdf_double('nlcutofftol',0.00001d0)
orbbuf=esdf_boolean('orbbuf',.true.)
jasbuf=esdf_boolean('jasbuf',.true.)
orb_norm=esdf_double('orb_norm',1.d0)
xyzzyaaco1=esdf_boolean('checkwfn',.false.)
xyzzyaacp1=esdf_boolean('kwarn',.false.)
ranprint=esdf_integer('ranprint',0)
random_seed_kw=esdf_string('random_seed','timer')
ranluxlevel=esdf_integer('ranluxlevel',3)
sparse=esdf_boolean('sparse',.false.)
sparse_threshold=esdf_double('sparse_threshold',1.d-12)
cusp_correction=esdf_boolean('cusp_correction',.true.)
cusp_info=esdf_boolean('cusp_info',.false.)
cusp_control=esdf_double('cusp_control',50.d0)
cusp_threshold=esdf_double('cusp_threshold',1.d-7)
use_gpcc=esdf_boolean('use_gpcc',.false.)
s_plot=esdf_boolean('splot',.false.)
ewald_control=esdf_double('ewald_control',0.d0)
xyzzyaacq1=esdf_boolean('ewald_check',.true.)
molgscreening=esdf_boolean('molgscreening',.false.)
bsmooth=esdf_boolean('bsmooth',.false.)
chkpoint_level=esdf_integer('checkpoint',1)
relativistic=esdf_boolean('relativistic',.false.)
xyzzyaacg1=esdf_double('isotope_mass',0.d0)
ke_verbose=esdf_boolean('ke_verbose',.false.)
ke_forgive=esdf_boolean('ke_forgive',.true.)
finite_size_corr=esdf_boolean('finite_size_corr',.false.)
xc_corr_method=esdf_integer('xc_corr_method',1)
forces=esdf_boolean('forces',.false.)
forces_info=esdf_integer('forces_info',2)
fix_holes=esdf_boolean('fix_holes',.false.)
allow_nochi_atoms=esdf_boolean('allow_nochi_atoms',.false.)
use_magnetic_field=esdf_boolean('magnetic_field',.false.)
checkpoint_ncpu=esdf_integer('checkpoint_ncpu',nnodes)
allow_slave_write=esdf_boolean('allow_slave_write',.true.)
rng_restart_safe=esdf_boolean('rng_restart_safe',.true.)
shm_size_nproc=esdf_integer('shm_size_nproc',0)
hartree_xc=esdf_boolean('hartree_xc',.false.)
makemovie=esdf_boolean('makemovie',.false.)
movieplot=esdf_integer('movieplot',1)
movienode=esdf_integer('movienode',0)
movie_supercells=esdf_boolean('moviecells',.false.)
density=esdf_boolean('density',.false.)
spin_density=esdf_boolean('spin_density',.false.)
pair_corr=esdf_boolean('pair_corr',.false.)
pair_corr_sph=esdf_boolean('pair_corr_sph',.false.)
loc_tensor=esdf_boolean('loc_tensor',.false.)
structure_factor=esdf_boolean('structure_factor',.false.)
structure_factor_sph=esdf_boolean('struc_factor_sph',.false.)
onep_density_mat=esdf_boolean('onep_density_mat',.false.)
twop_density_mat=esdf_boolean('twop_density_mat',.false.)
cond_fraction=esdf_boolean('cond_fraction',.false.)
eval_dipole_moment=esdf_boolean('dipole_moment',.false.)
mom_den=esdf_boolean('mom_den',.false.)
eval_contact_den=esdf_boolean('contact_den',.false.)
population=esdf_boolean('population',.false.)
twop_dm_mom=esdf_boolean('twop_dm_mom',.false.)
cond_frac_mom=esdf_boolean('cond_fraction_mom',.false.)
expval_cutoff=esdf_physical('expval_cutoff',30.d0,'hartree')
expval_error_bars=esdf_boolean('expval_error_bars',.false.)
permit_den_symm=esdf_boolean('permit_den_symm',.true.)
qmc_density_mpc=esdf_boolean('qmc_density_mpc',.false.)
int_sf=esdf_boolean('int_sf',.false.)
if(hartree_xc.and.structure_factor)int_sf=.true.
pcfs_nbins_in=esdf_integer('pcfs_nbins',-1)
pcfs_rcutoff_in=esdf_physical('pcfs_rcutoff',-1.d0,'bohr')
virtual_node=esdf_integer('virtual_node',0)
virtual_nconfig=esdf_integer('virtual_nconfig',0)
virtual_nnodes=esdf_integer('virtual_nnodes',1)
xyzzyaaax3=esdf_integer('btype',-1)
xyzzyaaaw3=esdf_integer('iterac',-1)
special_wfn=esdf_string('special_wfn','none')
xyzzyaaaz3=esdf_boolean('no_ee_int',.false.)
xyzzyaaay3=esdf_integer('nlrule1',-1)
iaccum=esdf_boolean('iaccumulate',.true.)
xyzzyaaba3=esdf_boolean('use_molorbmods',.false.)
xyzzyaabp1=esdf_integer('dmc_npops',1)
parallel_keywords=esdf_string('parallel_keywords','per_node')
if(esdf_defined('nlrule2','I'))call errwarn('SET_INPUT_PARAMETERS','ke&
&yword NLRULE2 is redundant and its value is ignored since v2.1. Reaso&
&n: dubious usefulness, functionality removed. Remove this keyword fro&
&m the input file to get rid of this message.')
if(esdf_defined('calc_variance','L'))call errwarn('SET_INPUT_PARAMETER&
&S','keyword CALC_VARIANCE is redundant and its value is ignored. Reas&
&on: functionality now always on. Remove this keyword from the input f&
&ile to get rid of this message.')
if(esdf_defined('vm_deriv_buffer','L'))call errwarn('SET_INPUT_PARAMET&
&ERS','keyword VM_DERIV_BUFFER is redundant and its value is ignored. &
&Reason: algorithm rewritten, functionality now always on. Remove this&
& keyword from the input file to get rid of this message.')
if(esdf_defined('vm_dist_buffer','L'))call errwarn('SET_INPUT_PARAMETE&
&RS','keyword VM_DIST_BUFFER is redundant and its value is ignored. Re&
&ason: algorithm rewritten, functionality not relevant. Remove this ke&
&yword from the input file to get rid of this message.')
if(esdf_defined('bf_save_memory','L'))call errwarn('SET_INPUT_PARAMETE&
&RS','keyword BF_SAVE_MEMORY is redundant and its value is ignored. Re&
&ason: algorithm rewritten without unneeded buffering, functionality n&
&ot relevant. Remove this keyword from the input file to get rid of th&
&is message.')
if(esdf_defined('emin_sampling','T'))call errwarn('SET_INPUT_PARAMETER&
&S','keyword EMIN_SAMPLING is redundant and its value is ignored. Reas&
&on: algorithm rewritten, functionality removed. Remove this keyword f&
&rom the input file to get rid of this message.')
if(esdf_defined('spin_density_mat','L'))call errwarn('SET_INPUT_PARAME&
&TERS','keyword SPIN_DENSITY_MAT is redundant and its value is ignored&
&. Reason: implied by DENSITY/SPIN_DENSITY=T in non-collinear system. &
&Remove this keyword from the input file to get rid of this message.')
if(esdf_defined('redist_period','I'))call errwarn('SET_INPUT_PARAMETER&
&S','keyword REDIST_PERIOD is redundant and its value is ignored. Reas&
&on: doesn''t work as intended - no benefit from any value other than &
&1. Remove this keyword from the input file to get rid of this message&
&.')
if(esdf_defined('num_cpus_in_group','I'))call errwarn('SET_INPUT_PARAM&
&ETERS','keyword NUM_CPUS_IN_GROUP is redundant and its value is ignor&
&ed. Reason: duplicate feature removed in favour of SHM.')
if(esdf_defined('use_mpiio','L'))call errwarn('SET_INPUT_PARAMETERS','&
&keyword USE_MPIIO is redundant and its value is ignored. Reason: feat&
&ure removed along with NUM_CPUS_IN_GROUP.')
select case(trim(atom_basis_type))
case('slater-type')
have_ae=.true.
case('non_int_he','h2','h3plus')
have_ae=.true.
case('nonint_he')
atom_basis_type='non_int_he'
have_ae=.true.
case default
have_ae=.false.
end select
xyzzyaacn1=esdf_defined('input_example','E')
if(xyzzyaacn1)then
if(am_master)then
call esdf_dump('input_example')
call wout()
call wordwrap("Written file 'input_example' containing all keywords th&
&at CASINO knows about in appropriate format.")
call wordwrap('NB: if block records not present in initial input file &
&then they will not appear in input example. Currently the possible bl&
&ock records are as follows (this list may be out of date): ')
call wout('custom_spair_dep')
call wout('custom_ssingle_dep')
call wout('custom_striplet_dep')
call wout('dtvmcs')
call wout('edist_by_ion')
call wout('edist_by_iontype')
call wout('expval_kgrid')
call wout('free_particles')
call wout('initial_config')
call wout('jastrow_plot')
call wout('manual_interaction')
call wout('npcell')
call wout('particles')
call wout('pcf_rfix')
call wout('on_top_pair')
call wout('plot_expval')
call wout('primitive_cell')
call wout('plot_backflow')
call wout('plot')
call wout('qmc_plot')
call wout()
call wout('See the manual or casinohelp for further details.')
call wout()
call wout(repeat('-',78))
call wout()
endif
call errstop_master('SET_INPUT_PARAMETERS','Quitting.')
endif
xyzzyaaan3=len_trim(con_loc)
if(xyzzyaaan3>=len(con_loc)-1)call errstop_master('SET_INPUT_PARAMETER&
&S','Directory name in CON_LOC is too long.')
if(con_loc(xyzzyaaan3:xyzzyaaan3)/="/")con_loc(xyzzyaaan3+1:xyzzyaaan3&
&+1)="/"
con_out=trim(con_loc)//trim(con_out)
con_in=trim(con_loc)//trim(con_in)
con_backup=trim(con_loc)//trim(con_backup)
old_input=esdf_defined('nmove','I').or.esdf_defined('nwrcon','I').or.e&
&sdf_defined('corper','I').or.esdf_defined('nvmcave','I').or.esdf_defi&
&ned('nblock','I').or.esdf_defined('nequil','I').or.esdf_defined('ncon&
&fig','D').or.esdf_defined('nmove_dmc_equil','I').or.esdf_defined('nbl&
&ock_dmc_equil','I').or.esdf_defined('nmove_dmc_stats','I').or.esdf_de&
&fined('nblock_dmc_stats','I').or.esdf_defined('ndmcave','I').or.esdf_&
&defined('corper_dmc','I').or.esdf_defined('trip_popn','D').or.esdf_de&
&fined('nequil_ta','I').or.esdf_defined('nblock_dmct_equil','I').or.es&
&df_defined('nmove_dmct_equil','I').or.esdf_defined('num_dmc_twists','&
&I').or.esdf_defined('vmc_twist_av','L').or.esdf_defined('nmove_dmcmd_&
&equil','I').or.esdf_defined('nmove_dmcmd_stats','I').or.esdf_defined(&
&'nconfig_prelim','I').or.esdf_defined('parallel_keywords','T')
new_input=esdf_defined('vmc_nstep','I').or.esdf_defined('vmc_nconfig_w&
&rite','I').or.esdf_defined('vmc_decorr_period','I').or.esdf_defined('&
&vmc_ave_period','I').or.esdf_defined('vmc_nblock','I').or.esdf_define&
&d('vmc_equil_nstep','I').or.esdf_defined('dmc_target_weight','D').or.&
&esdf_defined('dmc_equil_nstep','I').or.esdf_defined('dmc_stats_nstep'&
&,'I').or.esdf_defined('dmc_equil_nblock','I').or.esdf_defined('dmc_st&
&ats_nblock','I').or.esdf_defined('dmc_ave_period','I').or.esdf_define&
&d('dmc_decorr_period','I').or.esdf_defined('vmc_reequil_nstep','I').o&
&r.esdf_defined('dmc_reequil_nstep','I').or.esdf_defined('dmc_reequil_&
&nblock','I').or.esdf_defined('dmc_ntwist','I').or.esdf_defined('vmc_n&
&twist','I').or.esdf_defined('dmcmd_stats_nstep','I').or.esdf_defined(&
&'dmcmd_equil_nstep','I').or.esdf_defined('dmc_nconf_prelim','I').or.e&
&sdf_defined('dmc_trip_weight','D')
if(new_input.and.old_input)then
if(am_master)then
call wout('Keywords from the ''old'' set:')
call wordwrap('NMOVE, NWRCON, CORPER, NVMCAVE, NBLOCK, NEQUIL, NCONFIG&
&, NMOVE_DMC_EQUIL, NBLOCK_DMC_EQUIL, NMOVE_DMC_STATS, NBLOCK_DMC_STAT&
&S, NDMCAVE, CORPER_DMC, TRIP_POPN, NEQUIL_TA, NBLOCK_DMCT_EQUIL, NMOV&
&E_DMCT_EQUIL, NUM_DMC_TWISTS, VMC_TWIST_AV, NMOVE_DMCMD_EQUIL, NMOVE_&
&DMCMD_STATS, NCONFIG_PRELIM and PARALLEL_KEYWORDS')
call wout('cannot appear in the same input file as keywords from the '&
&'new'' set:')
call wordwrap('VMC_NSTEP, VMC_NCONFIG_WRITE, VMC_DECORR_PERIOD, VMC_AV&
&E_PERIOD, VMC_NBLOCK, VMC_EQUIL_NSTEP, DMC_TARGET_WEIGHT, DMC_EQUIL_N&
&STEP, DMC_STATS_NSTEP, DMC_EQUIL_NBLOCK, DMC_STATS_NBLOCK, DMC_AVE_PE&
&RIOD, DMC_DECORR_PERIOD, VMC_REEQUIL_NSTEP, DMC_REEQUIL_NSTEP, DMC_RE&
&EQUIL_NBLOCK, DMC_NTWIST, VMC_NTWIST, DMCMD_STATS_NSTEP, DMCMD_EQUIL_&
&NSTEP, DMC_NCONF_PRELIM and DMC_TRIP_WEIGHT')
endif
call errstop_master('SET_INPUT_PARAMETERS','Stopping.')
elseif(.not.(new_input.or.old_input))then
call errstop_master('SET_INPUT_PARAMETERS','Mandatory keywords missing!&
&')
endif
if(old_input)then
select case(trim(parallel_keywords))
case('per_node')
targ_wt=targ_wt*dble(nnodes)
trip_popn=trip_popn*dble(nnodes)
case('total')
xyzzyaaac3=xyzzyaaab1/nnodes
if(mod(xyzzyaaab1,nnodes)>0)xyzzyaaac3=xyzzyaaac3+1
xyzzyaaad3=xyzzyaaae1/nnodes
if(mod(xyzzyaaae1,nnodes)>0)xyzzyaaad3=xyzzyaaad3+1
if(xyzzyaaab1/=xyzzyaaac3*nnodes.or.xyzzyaaae1/=xyzzyaaad3*nnodes)call&
& errwarn('SET_INPUT_PARAMETERS','PARALLEL_KEYWORDS=total, but the val&
&ues of NMOVE and/or NWRCON are not divisible by the number of nodes. &
&These values have been increased by CASINO so that they are.')
xyzzyaaab1=xyzzyaaac3
xyzzyaaae1=xyzzyaaad3
nconfig_prelim=nconfig_prelim/nnodes
case default
call errstop_master('SET_INPUT_PARAMETERS','Unknown value for PARALLEL&
&_KEYWORDS.')
end select
if(xyzzyaacv1)then
xyzzyaabc1=1
block_time=0.
use_blocktime=.false.
xyzzyaaay1=xyzzyaaab1*xyzzyaaaj1*nnodes
xyzzyaabn1=xyzzyaaac1
else
xyzzyaabc1=xyzzyaaac1
xyzzyaaay1=xyzzyaaac1*xyzzyaaab1*xyzzyaaaj1*nnodes
xyzzyaabn1=1
endif
xyzzyaaba1=xyzzyaaaj1
xyzzyaabb1=xyzzyaaad1
xyzzyaaaz1=xyzzyaaae1*nnodes
xyzzyaabd1=xyzzyaaaa1
xyzzyaabk1=xyzzyaaaw1
if(dmc_md)then
xyzzyaaam1=1
xyzzyaaan1=1
xyzzyaaak1=xyzzyaacb1
xyzzyaaal1=xyzzyaacc1
block_time=0.
use_blocktime=.false.
endif
xyzzyaabf1=xyzzyaaam1
xyzzyaabm1=xyzzyaaau1
xyzzyaabh1=xyzzyaaan1
xyzzyaabe1=xyzzyaaam1*xyzzyaaak1*ndmcave/max(corper_dmc,1)
if(mod(xyzzyaaam1*xyzzyaaak1*ndmcave,max(corper_dmc,1))>0)xyzzyaabe1=x&
&yzzyaabe1+1
xyzzyaabl1=xyzzyaaau1*xyzzyaaat1*ndmcave/max(corper_dmc,1)
if(mod(xyzzyaaau1*xyzzyaaat1*ndmcave,max(corper_dmc,1))>0)xyzzyaabl1=x&
&yzzyaabl1+1
xyzzyaabg1=xyzzyaaan1*xyzzyaaal1*ndmcave/max(corper_dmc,1)
if(mod(xyzzyaaan1*xyzzyaaal1*ndmcave,max(corper_dmc,1))>0)xyzzyaabg1=x&
&yzzyaabg1+1
xyzzyaabi1=ndmcave/max(corper_dmc,1)
xyzzyaabj1=corper_dmc
xyzzyaaci1=targ_wt
xyzzyaacj1=trip_popn
xyzzyaabo1=xyzzyaaav1
if(use_blocktime)then
xyzzyaaak1=xyzzyaaak1*xyzzyaaam1
xyzzyaaal1=xyzzyaaal1*xyzzyaaan1
xyzzyaaam1=1
xyzzyaabf1=1
xyzzyaaan1=1
xyzzyaabh1=1
endif
else
xyzzyaaaj1=xyzzyaaba1
xyzzyaaad1=xyzzyaabb1
xyzzyaaaa1=xyzzyaabd1
xyzzyaaaw1=xyzzyaabk1
xyzzyaacv1=xyzzyaabn1>0
if(xyzzyaabn1>0)then
xyzzyaaac1=xyzzyaabn1
block_time=0.
use_blocktime=.false.
xyzzyaaab1=xyzzyaaay1/max(nnodes*xyzzyaaba1,1)
if(mod(xyzzyaaay1,max(nnodes*xyzzyaaba1,1))>0)xyzzyaaab1=xyzzyaaab1+1
xyzzyaaae1=xyzzyaaaz1/max(nnodes,1)
if(mod(xyzzyaaaz1,max(nnodes,1))>0)xyzzyaaae1=xyzzyaaae1+1
else
xyzzyaaac1=xyzzyaabc1
xyzzyaaab1=xyzzyaaay1/max(nnodes*xyzzyaabc1*xyzzyaaba1,1)
if(mod(xyzzyaaay1,max(nnodes*xyzzyaabc1*xyzzyaaba1,1))>0)xyzzyaaab1=xy&
&zzyaaab1+1
xyzzyaaae1=xyzzyaaaz1/max(nnodes,1)
if(mod(xyzzyaaaz1,max(nnodes,1))>0)xyzzyaaae1=xyzzyaaae1+1
endif
if(dmc_md)then
xyzzyaabf1=1
xyzzyaabh1=1
xyzzyaabe1=xyzzyaabz1
xyzzyaabg1=xyzzyaaca1
block_time=0.
use_blocktime=.false.
endif
xyzzyaaav1=xyzzyaabo1
xyzzyaaam1=xyzzyaabf1
xyzzyaaau1=xyzzyaabm1
xyzzyaaan1=xyzzyaabh1
corper_dmc=xyzzyaabj1
ndmcave=xyzzyaabi1*xyzzyaabj1
xyzzyaaak1=1+(xyzzyaabe1-1)/max(xyzzyaabf1*xyzzyaabi1,1)
xyzzyaaat1=1+(xyzzyaabl1-1)/max(xyzzyaabm1*xyzzyaabi1,1)
xyzzyaaal1=1+(xyzzyaabg1-1)/max(xyzzyaabh1*xyzzyaabi1,1)
targ_wt=xyzzyaaci1
trip_popn=xyzzyaacj1
xyzzyaabe1=xyzzyaaak1*max(xyzzyaabf1*xyzzyaabi1,1)
xyzzyaabl1=xyzzyaaat1*max(xyzzyaabm1*xyzzyaabi1,1)
xyzzyaabg1=xyzzyaaal1*max(xyzzyaabh1*xyzzyaabi1,1)
if(use_blocktime)then
xyzzyaaak1=xyzzyaaak1*xyzzyaaam1
xyzzyaaal1=xyzzyaaal1*xyzzyaaan1
xyzzyaaam1=1
xyzzyaabf1=1
xyzzyaaan1=1
xyzzyaabh1=1
endif
endif
if(nconfig_prelim>=0.and.xyzzyaaby1>=0)call errstop('SET_INPUT_PARAMET&
&ERS','Only one of NCONFIG_PRELIM and DMC_NCONF_PRELIM should be prese&
&nt in input.')
if(xyzzyaaby1>=0)nconfig_prelim=xyzzyaaby1/nnodes
if(nconfig_prelim<0)nconfig_prelim=0
xyzzyaadr1=xyzzyaabq1
xyzzyaadt1=xyzzyaabu1
xyzzyaadv1=xyzzyaabv1
xyzzyaads1=1+(xyzzyaabs1-1)/max(xyzzyaabu1*xyzzyaabr1,1)
xyzzyaadu1=1+(xyzzyaabt1-1)/max(xyzzyaabv1*xyzzyaabr1,1)
if(xyzzyaaax3/=-1)then
if(trim(atom_basis_type)/='default')call errstop_master('SET_INPUT_PAR&
&AMETERS','Cannot set BTYPE and ATOM_BASIS_TYPE at the same time. Remo&
&ve BTYPE from input and use ATOM_BASIS_TYPE instead.')
call errwarn('SET_INPUT_PARAMETERS','keyword BTYPE is redundant. Use A&
&TOM_BASIS_TYPE instead.')
select case(xyzzyaaax3)
case(0)
atom_basis_type='none'
case(1)
atom_basis_type='plane-wave'
case(2)
atom_basis_type='gaussian'
case(3)
atom_basis_type='numerical'
case(4)
atom_basis_type='blip'
case(5)
select case(trim(special_wfn))
case('noninthe')
atom_basis_type='non_int_he'
case default
call errstop_master('SET_INPUT_PARAMETERS','BTYPE=5 but SPECIAL_WFN va&
&lue not recognized. NB, BTYPE and SPECIAL_WFN are redundant. Use ATOM&
&_BASIS_TYPE instead.')
end select
case default
call errstop_master('SET_INPUT_PARAMETERS','Value of BTYPE not recogni&
&zed.')
end select
endif
select case(trim(atom_basis_type))
case('none','default')
atom_basis_type='none'
case('plane-wave','gaussian','slater-type','numerical','blip','dimer',&
&'non_int_he','h2','h3plus')
continue
case default
call errstop_master('SET_INPUT_PARAMETERS','Value of ATOM_BASIS_TYPE n&
&ot recognized.')
end select
if(xyzzyaaaw3/=-1)then
if(trim(interaction)/='default')call errstop_master('SET_INPUT_PARAMET&
&ERS','Cannot set ITERAC and INTERACTION at the same time. Remove ITER&
&AC from input and use INTERACTION instead.')
call errwarn('SET_INPUT_PARAMETERS','keyword ITERAC is redundant. Use &
&INTERACTION instead.')
select case(xyzzyaaaw3)
case(1)
if(isperiodic)then
interaction='ewald'
else
interaction='coulomb'
endif
case(2)
interaction='mpc'
case(3)
interaction='ewald_mpc'
case(4)
interaction='mpc_ewald'
case default
call errstop_master('SET_INPUT_PARAMETERS','Value of ITERAC not recogn&
&ized.')
end select
if(xyzzyaaaz3)then
interaction='none'
if(am_master)call errwarn('SET_INPUT_PARAMETERS','keyword NO_EE_INT is&
& redundant. Use INTERACTION instead.')
endif
endif
interaction1_present=.false.
interaction_mpc_present=.false.
interaction_mpc_use=.false.
int_name=""
select case(trim(interaction))
case('none')
continue
case('mpc')
if(.not.isperiodic)call errstop_master('SET_INPUT_PARAMETERS','Interac&
&tion type '//trim(interaction)//' cannot be used in aperiodic systems&
&.')
interaction_mpc_present=.true.
interaction_mpc_use=.true.
case('ewald_mpc','ewaldpp_mpc')
if(.not.isperiodic)call errstop_master('SET_INPUT_PARAMETERS','Interac&
&tion type '//trim(interaction)//' cannot be used in aperiodic systems&
&.')
if(trim(interaction)=='ewaldpp_mpc')then
interaction='ewald_mpc'
pseudo_ewald=.true.
interaction1_present=.true.
endif
interaction1_present=.true.
interaction_mpc_present=.true.
case('mpc_ewald','mpc_ewaldpp')
if(.not.isperiodic)call errstop_master('SET_INPUT_PARAMETERS','Interac&
&tion type '//trim(interaction)//' cannot be used in aperiodic systems&
&.')
if(trim(interaction)=='mpc_ewaldpp')then
interaction='mpc_ewald'
pseudo_ewald=.true.
interaction1_present=.true.
endif
interaction1_present=.true.
interaction_mpc_present=.true.
interaction_mpc_use=.true.
case('default','coulomb','ewald','ewaldpp')
if(trim(interaction)=='ewaldpp'.and.isperiodic)then
interaction='ewald'
pseudo_ewald=.true.
interaction1_present=.true.
endif
interaction='coulomb'
if(isperiodic)interaction='ewald'
interaction1_present=.true.
case('manual')
interaction1_present=.true.
interaction_mpc_present=.false.
case default
call errstop_master('SET_INPUT_PARAMETERS','Value of INTERACTION not r&
&ecognized.')
end select
if(xyzzyaaay3>=0)then
if(xyzzyaaas1>=0)call errstop_master('SET_INPUT_PARAMETERS','Cannot se&
&t NON_LOCAL_GRID and NLRULE1 at the same time. Use NON_LOCAL_GRID onl&
&y.')
call errwarn('SET_INPUT_PARAMETERS','keyword NLRULE1 is redundant. Use&
& NON_LOCAL_GRID instead.')
xyzzyaaas1=xyzzyaaay3
endif
if(xyzzyaaas1<0)xyzzyaaas1=4
if(esdf_defined('use_molorbmods','L'))then
if(esdf_defined('use_orbmods','L'))then
call errwarn('SET_INPUT_PARAMETERS','keyword USE_MOLORBMODS is redunda&
&nt. Use USE_ORBMODS instead. NB, value of USE_MOLORBMODS is being ign&
&ored for this run since USE_ORBMODS is also present in the input file&
&.')
else
call errwarn('SET_INPUT_PARAMETERS','keyword USE_MOLORBMODS is redunda&
&nt. Use USE_ORBMODS instead.')
use_orbmods=xyzzyaaba3
endif
endif
model_system=.false.
if(trim(atom_basis_type)=='none')then
if(xyzzyaaaf1>0.or.xyzzyaaag1>0)call errstop_master('SET_INPUT_PARAMET&
&ERS','NEU and NED refer to electrons in atomic orbitals. Hence ATOM_B&
&ASIS_TYPE=''none'' requires NEU=NED=0.')
if(xyzzyaaah1/=0.or.xyzzyaaai1/=0)call errstop_master('SET_INPUT_PARAM&
&ETERS','NHU and NHD refer to electrons in atomic orbitals. Hence ATOM&
&_BASIS_TYPE=''none'' requires NHU=NHD=0.')
xyzzyaaaf1=0
xyzzyaaag1=0
model_system=.true.
endif
select case(trim(psi_s))
case('none','slater','exmol','geminal','mahan_ex')
continue
case('pfaffian')
call errstop_master('SET_INPUT_PARAMETERS','Functional form '''//trim(&
&psi_s)//''' for Psi_S not implemented yet. Check back in a while.')
case default
call errstop_master('SET_INPUT_PARAMETERS','Functional form '''//trim(&
&psi_s)//''' for Psi_S unknown.')
end select
call init_wfn
if(complex_wf)then
real1_complex2=2
else
real1_complex2=1
endif
isvmc=.false.
isdmc=.false.
isvmc_dmc=.false.
isdmc_dmc=.false.
isvmc_dmc_equil=.false.
isopt=.false.
isvmc_opt=.false.
isopt_vmc=.false.
isgen_mpc=.false.
isrmc=.false.
isrmc_rmc=.false.
xyzzyaaen1=.false.
xyzzyaacy1=.false.
select case(trim(runtype))
case('gen_mdet_casl')
xyzzyaacy1=.true.
case('vmc')
isvmc=.true.
case('dmc')
if(esdf_defined('iaccumulate','L'))then
isdmc=.true.
isdmc_old=.true.
else
isdmc_dmc=.true.
endif
case('dmc_equil')
if(esdf_defined('iaccumulate','L'))call errstop_master('SET_INPUT_PARA&
&METERS','Cannot specify IACCUMULATE with RUNTYPE = dmc_equil. Remove &
&IACCUMULATE from the input file.')
isdmc=.true.
isdmc_equil=.true.
iaccum=.false.
case('dmc_stats')
if(esdf_defined('iaccumulate','L'))call errstop_master('SET_INPUT_PARA&
&METERS','Cannot specify IACCUMULATE with RUNTYPE = dmc_stats. Remove &
&IACCUMULATE from the input file.')
isdmc=.true.
isdmc_stats=.true.
iaccum=.true.
case('vmc_dmc')
isvmc_dmc=.true.
case('vmc_dmc_equil')
isvmc_dmc_equil=.true.
case('dmc_dmc')
isdmc_dmc=.true.
case('opt')
isopt=.true.
case('vmc_opt')
isvmc_opt=.true.
case('opt_vmc')
isopt_vmc=.true.
case('gen_mpc')
isgen_mpc=.true.
case('rmc')
isrmc=.true.
case('rmc_rmc')
isrmc_rmc=.true.
case('plot')
xyzzyaaen1=.true.
case default
call errstop_master('SET_INPUT_PARAMETERS','Run type not recognized.  &
&Check RUNTYPE value in the input file.')
end select
call xyzzyaaex1
call xyzzyaafb1
nelec=sum(nele(1:min(2,nspin)))
netot=sum(nele)
inv_netot=1.d0/dble(netot)
nemax=maxval(nele)
call xyzzyaafa1
if(esdf_block('dtvmcs',nlines))then
if(nlines/=no_difftypes)call errstop_master('SET_INPUT_PARAMETERS','Th&
&e number of lines in the DTVMCS block in the input file must equal th&
&e number of "diffusion types" in the system, which is '//trim(i2s(no_&
&difftypes))//' in the present case.')
do xyzzyaaah3=1,no_difftypes
read(block_data(xyzzyaaah3),*)xyzzyaaei1(xyzzyaaah3),xyzzyaaec1(xyzzya&
&aah3)
if(xyzzyaaec1(xyzzyaaah3)<0.or.xyzzyaaec1(xyzzyaaah3)>2)call errstop_m&
&aster('SET_INPUT_PARAMETERS','The OPT_DTVMC value must be 0, 1 or 2 i&
&n the DTVMCS block in input.')
enddo
xyzzyaacd1=xyzzyaaei1(1)
xyzzyaaaq1=xyzzyaaec1(1)
elseif(no_difftypes==1)then
xyzzyaaei1(:)=xyzzyaacd1
xyzzyaaec1(:)=xyzzyaaaq1
else
do xyzzyaaah3=1,no_difftypes
xyzzyaaei1(xyzzyaaah3)=xyzzyaacd1/difftype_mass(xyzzyaaah3)
enddo
xyzzyaaec1(:)=xyzzyaaaq1
endif
xyzzyaaej1(:)=xyzzyaacl1
allocate(xyzzyaaek1(3,netot),xyzzyaaeu1(netot),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'SET_INPUT_PARAMETERS','initial_rele')
xyzzyaaek1=0.d0
xyzzyaaeu1=.false.
if(esdf_block('initial_config',nlines))then
if(nlines<1.or.nlines>netot)call errstop_master('SET_INPUT_PARAMETERS'&
&,'The number of lines in the INITIAL_CONFIG block must be between 1 a&
&nd the total number of particles.')
do xyzzyaaah3=1,nlines
xyzzyaaau3=0.d0
read(block_data(xyzzyaaah3),*)xyzzyaaaf3,xyzzyaaai3,xyzzyaaau3(1:dimen&
&sionality)
if(xyzzyaaaf3<1.or.xyzzyaaaf3>nspin)call errstop_master('SET_INPUT_PAR&
&AMETERS','Spin index out of range at line '//trim(i2s(xyzzyaaah3))//'&
& of INITIAL_CONFIG block.')
if(xyzzyaaai3<1.or.xyzzyaaai3>nele(xyzzyaaaf3))call errstop_master('SE&
&T_INPUT_PARAMETERS','Particle index out of range at line '//trim(i2s(&
&xyzzyaaah3))//' of INITIAL_CONFIG block.')
xyzzyaaeu1(which_ii(xyzzyaaai3,xyzzyaaaf3))=.true.
xyzzyaaek1(1:3,which_ii(xyzzyaaai3,xyzzyaaaf3))=xyzzyaaau3(1:3)
enddo
endif
if(esdf_block('npcell',nlines))then
if(nlines/=1)call errstop_master('SET_INPUT_PARAMETERS','Too many line&
&s in NPCELL block in input file.')
read(block_data(1),*,iostat=xyzzyaadl1)scell_matrix(1,1),scell_matrix(&
&2,2),scell_matrix(3,3)
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding NPCELL block.')
endif
if(esdf_block('scell_matrix',nlines))then
if(nlines/=3)call errstop_master('SET_INPUT_PARAMETERS','There should &
&be three lines in the SCELL_MATRIX block in input file.')
do xyzzyaaah3=1,3
read(block_data(xyzzyaaah3),*,iostat=xyzzyaadl1)scell_matrix(xyzzyaaah&
&3,1:3)
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding line '//trim(i2s(xyzzyaaah3))//' of the SCELL_MATRIX block.')
enddo
if(esdf_block('npcell',nlines))call errstop_master('SET_INPUT_PARAMETE&
&RS','Should not have both NPCELL and SCELL_MATRIX in input.')
endif
if(trim(interaction)=='manual'.or.pseudo_ewald)then
if(esdf_block('manual_interaction',nlines))then
call xyzzyaafx1(nlines)
else
call errstop_master('SET_INPUT_PARAMETERS','INTERACTION is set to manu&
&al or pseudo-Ewald but there is no corresponding MANUAL_INTERACTION b&
&lock.')
endif
elseif(esdf_block('manual_interaction',nlines))then
call errstop_master('SET_INPUT_PARAMETERS','INTERACTION is not set to &
&manual or pseudo-Ewald but there is a MANUAL_INTERACTION block presen&
&t in input.')
endif
if(esdf_block('edist_by_ion',nlines))then
xyzzyaaav3=.false.
init_by_ion=.true.
xyzzyaadm1=nlines
allocate(neion(nlines,nspin),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'SET_INPUT_PARAMETERS','neion')
do xyzzyaaah3=1,nlines
read(block_data(xyzzyaaah3),*)xyzzyaaak3,neion(xyzzyaaak3,1:nspin)
enddo
if(sum(neion)/=(xyzzyaaaf1+xyzzyaaag1+xyzzyaaah1+xyzzyaaai1))then
if(am_master)call wout('Number of electrons in EDIST_BY_ION block not &
&equal to NEU+NED+NHU+NHD.')
xyzzyaaav3=.true.
endif
if(sum(neion(:,1))/=xyzzyaaaf1)then
if(am_master)call wout('Number of spin-up electrons in EDIST_BY_ION bl&
&ock not equal to NEU.')
xyzzyaaav3=.true.
endif
if(sum(neion(:,2))/=xyzzyaaag1)then
if(am_master)call wout('Number of spin-down electrons in EDIST_BY_ION &
&block not equal to NED.')
xyzzyaaav3=.true.
endif
if(xyzzyaaah1>0)then
if(sum(neion(:,3))/=xyzzyaaah1)then
if(am_master)call wout('Number of spin-up holes in EDIST_BY_ION block &
&not equal to NHU.')
xyzzyaaav3=.true.
endif
endif
if(xyzzyaaai1>0)then
if(sum(neion(:,4))/=xyzzyaaai1)then
if(am_master)call wout('Number of spin-down holes in EDIST_BY_ION bloc&
&k not equal to NHD.')
xyzzyaaav3=.true.
endif
endif
if(xyzzyaaav3)call errstop_master('SET_INPUT_PARAMETERS','Quitting.')
endif
if(esdf_block('edist_by_iontype',nlines))then
init_by_iontype=.true.
allocate(xyzzyaadx1(nlines,nspin),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'SET_INPUT_PARAMETERS','neiontype')
do xyzzyaaah3=1,nlines
read(block_data(xyzzyaaah3),*)xyzzyaaai3,xyzzyaadx1(xyzzyaaah3,1:nspin&
&)
enddo
endif
if(init_by_ion.and.init_by_iontype)call errstop_master('SET_INPUT_PARA&
&METERS','Choose one of EDIST_BY_ION and EDIST_BY_IONTYPE in input. Ca&
&nnot use both.')
if(xyzzyaaen1)then
if(esdf_block('plot',nlines))then
call init_plot(nlines,block_data)
xyzzyaadw1=1
elseif(esdf_block('qmc_plot',nlines))then
call init_plotter(nlines,block_data)
xyzzyaadw1=2
else
call errstop_master('SET_INPUT_PARAMETERS','RUNTYPE is ''plot'', but t&
&here is no PLOT or QMC_PLOT block in the input file.')
endif
endif
if(esdf_block('jastrow_plot',nlines))then
if(nlines/=6.and.nlines/=3)call errstop_master('SET_INPUT_PARAMETERS',&
&'There should be either six or three lines in the jastrow_plot block &
&in input.')
read(block_data(1),*,err=4,end=4)idum1
if(idum1<0.or.idum1>1)call errstop_master('SET_INPUT_PARAMETERS','Firs&
&t line of jastrow_plot should be 0 or 1')
read(block_data(2),*,err=4,end=4)xyzzyaaae3
read(block_data(3),*,err=4,end=4)xyzzyaaaf3
if(xyzzyaaae3<1.or.xyzzyaaae3>nspin.or.xyzzyaaaf3<1.or.xyzzyaaaf3>nspi&
&n)call errstop_master('SET_INPUT_PARAMETERS','Particle spins out of r&
&ange in jastrow_plot block in input.')
if(nlines==6)then
read(block_data(4),*,err=4,end=4)xyzzyaaap3(1:3)
read(block_data(5),*,err=4,end=4)xyzzyaaaq3(1:3)
read(block_data(6),*,err=4,end=4)xyzzyaaar3(1:3)
else
xyzzyaaap3=0.d0
xyzzyaaaq3=(/1.d0,0.d0,0.d0/)
xyzzyaaar3=0.d0
endif
call setup_jastrow_plot(idum1==1,xyzzyaaae3,xyzzyaaaf3,xyzzyaaap3,xyzz&
&yaaaq3,xyzzyaaar3)
else
call setup_jastrow_plot(.false.)
endif
if(esdf_block('plot_backflow',nlines))then
if(nlines<2.or.nlines>3)call errstop_master('SET_INPUT_PARAMETERS','Th&
&ere should be two/three lines in the plot_backflow block in input.')
read(block_data(1),*,err=5,end=5)xyzzyaaal3
if(xyzzyaaal3<0.or.xyzzyaaal3>1)call errstop_master('SET_INPUT_PARAMET&
&ERS','First line of plot_backflow should be 0 or 1')
if(xyzzyaaal3==1)then
read(block_data(2),*,err=5,end=5)xyzzyaaal3,xyzzyaaam3,xyzzyaaas3
if(xyzzyaaal3<1.or.xyzzyaaal3>nspin)call errstop_master('SET_INPUT_PAR&
&AMETERS','Quasiparticle spin out of range in plot_backflow block in i&
&nput.')
if(xyzzyaaam3<1.or.xyzzyaaam3>nele(xyzzyaaal3))call errstop_master('SE&
&T_INPUT_PARAMETERS','Quasiparticle number out of range in plot_backfl&
&ow block in input.')
if(nlines==3)then
read(block_data(3),*,err=5,end=5)xyzzyaaat3
if(xyzzyaaat3<0.d0)call errstop_master('SET_INPUT_PARAMETERS','Electro&
&n-ion distance out of range in plot_backflow block in input.')
endif
call setup_backflow_plot(.true.,nlines==3,xyzzyaaam3,xyzzyaaas3,xyzzya&
&aat3)
else
call setup_backflow_plot(.false.)
endif
else
call setup_backflow_plot(.false.)
endif
fixed_particle=0
input_type_fixed=-999
if(esdf_block('pcf_rfix',nlines))then
if(nlines/=2)call errstop_master('SET_INPUT_PARAMETERS','PCF_RFIX bloc&
&k in input should contain two lines.')
read(block_data(1),*,iostat=xyzzyaadl1)input_type_fixed
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding type of fixed particle in PCF_RFIX block in input.')
input_rfix=0.d0
read(block_data(2),*,iostat=xyzzyaadl1)input_rfix(1:dimensionality)
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding coordinates of fixed particle in PCF_RFIX block in input.')
endif
on_top_ii=0
on_top_jj=0
if(esdf_block('on_top_pair',nlines))then
if(nlines/=2)call errstop_master('SET_INPUT_PARAMETERS','ON_TOP_PAIR b&
&lock in input should contain two lines.')
read(block_data(1),*,iostat=xyzzyaadl1)xyzzyaaae3,xyzzyaaag3
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding first particle in ON_TOP_PAIR block in input.')
if(xyzzyaaae3<1.or.xyzzyaaae3>nspin)call errstop_master('SET_INPUT_PAR&
&AMETERS','Bad particle type index reading first line of ON_TOP_PAIR i&
&nput block.')
if(xyzzyaaag3<1.or.xyzzyaaag3>nele(xyzzyaaae3))call errstop_master('SE&
&T_INPUT_PARAMETERS','Bad particle index reading first line of ON_TOP_&
&PAIR input block.')
on_top_ii=which_ii(xyzzyaaag3,xyzzyaaae3)
read(block_data(2),*,iostat=xyzzyaadl1)xyzzyaaae3,xyzzyaaag3
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding second particle in ON_TOP_PAIR block in input.')
if(xyzzyaaae3<1.or.xyzzyaaae3>nspin)call errstop_master('SET_INPUT_PAR&
&AMETERS','Bad particle type index reading second line of ON_TOP_PAIR &
&input block.')
if(xyzzyaaag3<1.or.xyzzyaaag3>nele(xyzzyaaae3))call errstop_master('SE&
&T_INPUT_PARAMETERS','Bad particle index reading second line of ON_TOP&
&_PAIR input block.')
on_top_jj=which_ii(xyzzyaaag3,xyzzyaaae3)
if(on_top_ii==on_top_jj)call errstop_master('SET_INPUT_PARAMETERS','Th&
&e two lines in ON_TOP_PAIR block specify the same particle.')
endif
expval_nkdim(1:3)=-999
expval_nk(:,:)=0
if(esdf_block('expval_kgrid',nlines))then
if(nlines<4.or.nlines>1+nexpvals_need_kgrid*6)call errstop_master('SET&
&_INPUT_PARAMETERS','EXPVAL_KGRID block in input has wrong number of l&
&ines.')
read(block_data(1),*,iostat=xyzzyaadl1)expval_nkgrids
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding number of k grids in EXPVAL_KGRID input.')
if(expval_nkgrids<1)call errstop_master('SET_INPUT_PARAMETERS','Error &
&reading number of k grids in EXPVAL_KGRID input.')
if(expval_nkgrids>nexpvals_need_kgrid)call errstop2('SET_INPUT_PARAMET&
&ERS','Too many expval k grids defined. Maximum number required: ',nex&
&pvals_need_kgrid)
xyzzyaaaa3=2
do xyzzyaaah3=1,expval_nkgrids
read(block_data(xyzzyaaaa3),*,iostat=xyzzyaadl1)expval_kgrid_id(xyzzya&
&aah3)
xyzzyaaaa3=xyzzyaaaa3+1
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding grid ID in EXPVAL_KGRID input.')
xyzzyaaab3=expval_kgrid_id(xyzzyaaah3)
if(xyzzyaaab3<1.or.xyzzyaaab3>nexpvals_need_kgrid)call errstop_master(&
&'SET_INPUT_PARAMETERS','Invalid grid ID in EXPVAL_KGRID input.')
read(block_data(xyzzyaaaa3),*,iostat=xyzzyaadl1)expval_nkdim(xyzzyaaab&
&3)
xyzzyaaaa3=xyzzyaaaa3+1
if(expval_nkdim(xyzzyaaab3)<1.or.expval_nkdim(xyzzyaaab3)>3)call errst&
&op_master('SET_INPUT_PARAMETERS','Invalid grid dimensionality in EXPV&
&AL_KGRID input. Must be 1, 2 or 3.')
if(xyzzyaaab3==2)then
if(expval_nkdim(xyzzyaaab3)/=1)call errstop_master('SET_INPUT_PARAMETE&
&RS','Invalid grid dimensionality for radial k point grid in input EXP&
&VAL_KGRID. Must equal 1.')
read(block_data(xyzzyaaaa3),*,iostat=xyzzyaadl1)expval_ka(1,xyzzyaaab3&
&)
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding k grid coordinates for point A in EXPVAL_KGRID input.')
read(block_data(xyzzyaaaa3),*,err=6,end=6)expval_ka(1,xyzzyaaab3),expv&
&al_ka(2,xyzzyaaab3)
call errstop_master('SET_INPUT_PARAMETERS','Radial k grid required but&
& vector k points apparently given in EXPVAL_GRID input.')
6   xyzzyaaaa3=xyzzyaaaa3+1
else
read(block_data(xyzzyaaaa3),*,iostat=xyzzyaadl1)expval_ka(1,xyzzyaaab3&
&),expval_ka(2,xyzzyaaab3),expval_ka(3,xyzzyaaab3)
xyzzyaaaa3=xyzzyaaaa3+1
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding k grid coordinates for point A in EXPVAL_KGRID input.')
endif
if(xyzzyaaab3==2)then
read(block_data(xyzzyaaaa3),*,iostat=xyzzyaadl1)expval_kb(1,xyzzyaaab3&
&),expval_nk(1,xyzzyaaab3)
xyzzyaaaa3=xyzzyaaaa3+1
else
read(block_data(xyzzyaaaa3),*,iostat=xyzzyaadl1)expval_kb(1,xyzzyaaab3&
&),expval_kb(2,xyzzyaaab3),expval_kb(3,xyzzyaaab3),expval_nk(1,xyzzyaa&
&ab3)
xyzzyaaaa3=xyzzyaaaa3+1
endif
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding k grid coordinates for point B in EXPVAL_KGRID input.')
if(expval_nkdim(xyzzyaaah3)>1)then
read(block_data(xyzzyaaaa3),*,iostat=xyzzyaadl1)expval_kc(1,xyzzyaaab3&
&),expval_kc(2,xyzzyaaab3),expval_kc(3,xyzzyaaab3),expval_nk(2,xyzzyaa&
&ab3)
xyzzyaaaa3=xyzzyaaaa3+1
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding k grid coordinates for point C in EXPVAL_KGRID input.')
if(expval_nkdim(xyzzyaaah3)==2)then
xyzzyaaap3(1:3)=expval_kb(1:3,xyzzyaaab3)-expval_ka(1:3,xyzzyaaab3)
xyzzyaaaq3(1:3)=expval_kb(1:3,xyzzyaaab3)-expval_ka(1:3,xyzzyaaab3)
xyzzyaaao3=a1(1)*a2(2)-a1(2)*a2(1)
if(xyzzyaaao3==0.d0)call errstop_master('SET_INPUT_PARAMETERS','Vector&
&s AB and AC defining the expval k grid are collinear. Change EXPVAL_K&
&GRID in input.')
endif
endif
if(expval_nkdim(xyzzyaaah3)==3)then
read(block_data(xyzzyaaaa3),*,iostat=xyzzyaadl1)expval_kd(1,xyzzyaaab3&
&),expval_kd(2,xyzzyaaab3),expval_kd(3,xyzzyaaab3),expval_nk(3,xyzzyaa&
&ab3)
xyzzyaaaa3=xyzzyaaaa3+1
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Error rea&
&ding k grid coordinates for point D in EXPVAL_KGRID input.')
xyzzyaaap3(1:3)=expval_kb(1:3,xyzzyaaab3)-expval_ka(1:3,xyzzyaaab3)
xyzzyaaaq3(1:3)=expval_kc(1:3,xyzzyaaab3)-expval_ka(1:3,xyzzyaaab3)
xyzzyaaar3(1:3)=expval_kd(1:3,xyzzyaaab3)-expval_ka(1:3,xyzzyaaab3)
xyzzyaaao3=xyzzyaaap3(1)*xyzzyaaaq3(2)*xyzzyaaar3(3)+xyzzyaaap3(2)*xyz&
&zyaaaq3(3)*xyzzyaaar3(1)+xyzzyaaap3(3)*xyzzyaaaq3(1)*xyzzyaaar3(2)-xy&
&zzyaaap3(3)*xyzzyaaaq3(2)*xyzzyaaar3(1)-xyzzyaaap3(1)*xyzzyaaaq3(3)*x&
&yzzyaaar3(2)-xyzzyaaap3(2)*xyzzyaaaq3(1)*xyzzyaaar3(3)
if(xyzzyaaao3==0.d0)call errstop_master('SET_INPUT_PARAMETERS','Vector&
&s AB, AC and AD defining the expval k grid are linearly dependent. Ch&
&ange EXPVAL_KGRID in input.')
endif
enddo
if((xyzzyaaaa3-1)/=nlines)call errstop_master('SET_INPUT_PARAMETERS','&
&Wrong number of lines in EXPVAL_KGRID in input.')
endif
if(trim(runtype)=='opt')opt_cycles=1
if(esdf_block('opt_plan',nlines))opt_cycles=nlines
if(trim(runtype)=='opt'.and.opt_cycles/=1)call errstop_master('SET_INP&
&UT_PARAMETERS','OPT_PLAN input block must contain exactly one line if&
& RUNTYPE is "opt".')
allocate(xyzzyaacz1(opt_cycles),xyzzyaada1(opt_cycles),xyzzyaadb1(opt_&
&cycles),xyzzyaadc1(opt_cycles),xyzzyaadd1(opt_cycles),xyzzyaade1(opt_&
&cycles),xyzzyaadf1(opt_cycles),xyzzyaadg1(opt_cycles),stat=xyzzyaadk1&
&)
call check_alloc(xyzzyaadk1,'SET_INPUT_PARAMETERS','opt_plan')
xyzzyaacz1(:)=opt_method
xyzzyaada1(:)=opt_jastrow
xyzzyaadb1(:)=opt_backflow
xyzzyaadc1(:)=opt_det_coeff
xyzzyaadd1(:)=opt_orbitals
xyzzyaade1(:)=opt_geminal
xyzzyaadg1(:)=opt_maxiter
xyzzyaadf1(:)=.false.
xyzzyaadh1=0
xyzzyaadi1=0
xyzzyaadj1=0
if(xyzzyaaax1>0)xyzzyaadf1(max(1,opt_cycles-xyzzyaaax1+1):opt_cycles)=&
&.true.
if(esdf_block('opt_plan',nlines))then
do xyzzyaaah3=1,nlines
read(block_data(xyzzyaaah3),*,iostat=xyzzyaadl1)xyzzyaaaj3
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Problem r&
&eading cycle index at line '//trim(i2s(xyzzyaaah3))//' of OPT_PLAN in&
&put block.')
if(xyzzyaaaj3/=xyzzyaaah3)call errstop_master('SET_INPUT_PARAMETERS','&
&Cycle index mismatch at line '//trim(i2s(xyzzyaaah3))//' of OPT_PLAN &
&input block.')
xyzzyaaab3=0
do
read(block_data(xyzzyaaah3),*,iostat=xyzzyaadl1)xyzzyaaaj3,(token,xyzz&
&yaaai3=1,xyzzyaaab3+1)
if(xyzzyaadl1/=0)exit
if(len_trim(token)==0)exit
xyzzyaaab3=xyzzyaaab3+1
enddo
do xyzzyaaaa3=1,xyzzyaaab3
read(block_data(xyzzyaaah3),*)xyzzyaaaj3,(token,xyzzyaaai3=1,xyzzyaaaa&
&3)
xyzzyaaaj3=scan(token,'=')
if(xyzzyaaaj3<1)call errstop_master('SET_INPUT_PARAMETERS','"=" sign m&
&issing in token '//trim(i2s(xyzzyaaaa3))//' at line '//trim(i2s(xyzzy&
&aaah3))//' of OPT_PLAN input block.')
select case(token(1:xyzzyaaaj3-1))
case('method')
read(token(xyzzyaaaj3+1:),*,iostat=xyzzyaadl1)xyzzyaacz1(xyzzyaaah3)
select case(trim(xyzzyaacz1(xyzzyaaah3)))
case('varmin','varmin_linjas','madmin','emin')
continue
case default
call errstop_master('SET_INPUT_PARAMETERS','Unknown value for method i&
&n OPT_PLAN input block.')
end select
case('jastrow')
read(token(xyzzyaaaj3+1:),*,iostat=xyzzyaadl1)xyzzyaada1(xyzzyaaah3)
case('backflow')
read(token(xyzzyaaaj3+1:),*,iostat=xyzzyaadl1)xyzzyaadb1(xyzzyaaah3)
case('det_coeff')
read(token(xyzzyaaaj3+1:),*,iostat=xyzzyaadl1)xyzzyaadc1(xyzzyaaah3)
case('orbitals')
read(token(xyzzyaaaj3+1:),*,iostat=xyzzyaadl1)xyzzyaadd1(xyzzyaaah3)
case('geminal')
read(token(xyzzyaaaj3+1:),*,iostat=xyzzyaadl1)xyzzyaade1(xyzzyaaah3)
case('maxiter')
read(token(xyzzyaaaj3+1:),*,iostat=xyzzyaadl1)xyzzyaadg1(xyzzyaaah3)
case('fix_cutoffs')
read(token(xyzzyaaaj3+1:),*,iostat=xyzzyaadl1)xyzzyaadf1(xyzzyaaah3)
case default
call errstop_master('SET_INPUT_PARAMETERS','Unknown keyword in token '&
&//trim(i2s(xyzzyaaaa3))//' at line '//trim(i2s(xyzzyaaah3))//' of OPT&
&_PLAN input block.')
end select
if(xyzzyaadl1/=0)call errstop_master('SET_INPUT_PARAMETERS','Problem r&
&eading value of token '//trim(i2s(xyzzyaaaa3))//' at line '//trim(i2s&
&(xyzzyaaah3))//' of OPT_PLAN input block.')
enddo
enddo
endif
if(.not.model_system)then
npcells=abs(scell_matrix(1,1)*(scell_matrix(2,2)*scell_matrix(3,3)-sce&
&ll_matrix(2,3)*scell_matrix(3,2))+scell_matrix(2,1)*(scell_matrix(3,2&
&)*scell_matrix(1,3)-scell_matrix(1,2)*scell_matrix(3,3))+scell_matrix&
&(3,1)*(scell_matrix(1,2)*scell_matrix(2,3)-scell_matrix(1,3)*scell_ma&
&trix(2,2)))
if(npcells==0.and.isperiodic)call errstop_master('SET_INPUT_PARAMETERS&
&','You need to specify the npcell array in input.')
endif
return
4 call errstop_master('SET_INPUT_PARAMETERS','Problem reading the JAST&
&ROW_PLOT block in input.')
5 call errstop_master('SET_INPUT_PARAMETERS','Problem reading the PLOT&
&_BACKFLOW block in input.')
end subroutine xyzzyaaew1
subroutine xyzzyaaex1
implicit none
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4,xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4,xyzzyaaak4,xyzzyaaal4&
&,xyzzyaaam4,xyzzyaaan4
integer,allocatable :: xyzzyaaao4(:),xyzzyaaap4(:),xyzzyaaaq4(:,:),xyz&
&zyaaar4(:,:),xyzzyaaas4(:)
real(dp) xyzzyaaat4
real(dp),allocatable :: xyzzyaaau4(:,:,:),xyzzyaaav4(:),xyzzyaaaw4(:),&
&xyzzyaaax4(:),xyzzyaaay4(:),xyzzyaaaz4(:,:,:)
real(dp),parameter :: xyzzyaaba4=1.d-8
logical xyzzyaabb4,xyzzyaabc4,xyzzyaabd4,xyzzyaabe4,xyzzyaabf4
logical,allocatable :: xyzzyaabg4(:),xyzzyaabh4(:),xyzzyaabi4(:)
character(11) temp_slatt,temp_ferro,temp_off,temp_ctype,temp_sites,tem&
&p_orb,temp_det,temp_colon,temp_repeat,temp_mass
character(80) temp_string,temp_orbtype
character(64),allocatable :: xyzzyaabj4(:)
xyzzyaabc4=.false.
electron_system=.true.
xyzzyaaak4=2
heg_nlayers=1
allocate(heg_ylayer(heg_nlayers),heg_zlayer(heg_nlayers),stat=xyzzyaad&
&k1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_ylayer')
heg_ylayer=0.d0
heg_zlayer=0.d0
if(esdf_block('particles',xyzzyaaaj4))then
if(xyzzyaaaj4>0)then
xyzzyaaak4=xyzzyaaaj4+2
allocate(xyzzyaaav4(xyzzyaaak4),xyzzyaaaw4(xyzzyaaak4),xyzzyaabi4(xyzz&
&yaaak4),xyzzyaaay4(xyzzyaaak4),xyzzyaabj4(xyzzyaaak4),xyzzyaabg4(xyzz&
&yaaak4),xyzzyaaau4(3,3,xyzzyaaak4),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','PARTICLES block')
xyzzyaabg4=.true.
xyzzyaaau4=0.d0
xyzzyaabi4=.false.
xyzzyaaaa4=0
do
xyzzyaaaa4=xyzzyaaaa4+1
if(xyzzyaaaa4>xyzzyaaaj4)exit
temp_string=adjustl(block_data(xyzzyaaaa4))
xyzzyaaae4=scan(temp_string,' ')
read(temp_string(1:xyzzyaaae4-1),*,err=10,end=10)xyzzyaaah4
if(xyzzyaaah4<3)xyzzyaabc4=.true.
if(xyzzyaaah4>xyzzyaaak4)call errstop_master('READ_PARTICLES','Particl&
&e index must not skip numbers.  Detected at line '//trim(i2s(xyzzyaaa&
&a4))//' of PARTICLES block, when declaring particle '//trim(i2s(xyzzy&
&aaah4))//' of a maximum total of '//trim(i2s(xyzzyaaak4))//'.')
if(any(xyzzyaaau4(:,:,xyzzyaaah4)/=0.d0))call errstop_master('READ_PAR&
&TICLES','Particle '//trim(i2s(xyzzyaaah4))//' defined more than once &
&in PARTICLES block.')
temp_string=adjustl(temp_string(xyzzyaaae4+1:))
xyzzyaaae4=scan(temp_string,' ')
read(temp_string(1:xyzzyaaae4),*,err=10,end=10)xyzzyaaav4(xyzzyaaah4)
temp_string=adjustl(temp_string(xyzzyaaae4+1:))
xyzzyaaae4=scan(temp_string,' ')
read(temp_string(1:xyzzyaaae4),*,iostat=xyzzyaadl1)xyzzyaaaw4(xyzzyaaa&
&h4)
if(xyzzyaadl1/=0)then
read(temp_string(1:xyzzyaaae4),*,err=10,end=10)temp_mass
select case(trim(adjustl(temp_mass)))
case('heavy')
xyzzyaaaw4(xyzzyaaah4)=-2.d0
xyzzyaabi4(xyzzyaaah4)=.true.
case('anisotropic')
xyzzyaaaw4(xyzzyaaah4)=-1.d0
xyzzyaabg4(xyzzyaaah4)=.false.
case default
call errstop_master('READ_PARTICLES','Unrecognized value for the parti&
&cle mass in the PARTICLES block in input. Must be a mass (au) or one &
&of ''anisotropic'' or ''heavy''.')
end select
elseif(xyzzyaaaw4(xyzzyaaah4)<0.d0)then
xyzzyaabg4(xyzzyaaah4)=.false.
endif
temp_string=adjustl(temp_string(xyzzyaaae4+1:))
xyzzyaaae4=scan(temp_string,' ')
read(temp_string(1:xyzzyaaae4),*,err=10,end=10)xyzzyaaay4(xyzzyaaah4)
xyzzyaabj4(xyzzyaaah4)=adjustl(temp_string(xyzzyaaae4+1:))
if(.not.xyzzyaabg4(xyzzyaaah4))then
do xyzzyaaai4=1,3
xyzzyaaaa4=xyzzyaaaa4+1
if(xyzzyaaaa4>xyzzyaaaj4)call errstop_master('READ_PARTICLES','Reached&
& end of block while reading anisotropic mass tensor in PARTICLES bloc&
&k.')
read(block_data(xyzzyaaaa4),*,err=10,end=10)xyzzyaaau4(xyzzyaaai4,1:3,&
&xyzzyaaah4)
if(all(xyzzyaaau4(:,:,xyzzyaaah4)==0.d0))call errstop_master('READ_PAR&
&TICLES','Zero mass tensor found when reading PARTICLES block.')
enddo
elseif(xyzzyaaaw4(xyzzyaaah4)==0.d0)then
call errstop_master('READ_PARTICLES','Zero mass found when reading PAR&
&TICLES block.')
else
xyzzyaaau4(1,1,xyzzyaaah4)=xyzzyaaaw4(xyzzyaaah4)
xyzzyaaau4(2,2,xyzzyaaah4)=xyzzyaaaw4(xyzzyaaah4)
xyzzyaaau4(3,3,xyzzyaaah4)=xyzzyaaaw4(xyzzyaaah4)
endif
enddo
xyzzyaaaa4=2
if(xyzzyaabc4)xyzzyaaaa4=0
xyzzyaaai4=xyzzyaaaa4+1
do xyzzyaaah4=xyzzyaaai4,xyzzyaaak4
if(any(xyzzyaaau4(:,:,xyzzyaaah4)/=0.d0))xyzzyaaaa4=xyzzyaaaa4+1
enddo
do xyzzyaaah4=3,xyzzyaaaa4
if(all(xyzzyaaau4(:,:,xyzzyaaah4)==0.d0))call errstop_master('READ_PAR&
&TICLES','Particle index must not skip numbers (except 1 and 2).  Part&
&icle type '//trim(i2s(xyzzyaaah4))//' undefined in PARTICLES block in&
& input.')
enddo
xyzzyaaak4=xyzzyaaaa4
endif
endif
allocate(pcharge(xyzzyaaak4),pmass(xyzzyaaak4),inv_pmass(xyzzyaaak4),p&
&spin(xyzzyaaak4),pname(xyzzyaaak4),pisotropic(xyzzyaaak4),pmasstensor&
&(3,3,xyzzyaaak4),inv_pmasstensor(3,3,xyzzyaaak4),pfermion(xyzzyaaak4)&
&,popp_spin(xyzzyaaak4),pinfmass(xyzzyaaak4),heg_layer(xyzzyaaak4),sta&
&t=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','arrays for particle prop&
&erties')
pisotropic=.true.
pmasstensor=0.d0
inv_pmasstensor=0.d0
popp_spin=0
pinfmass=.false.
heg_layer=1
xyzzyaaal4=min(xyzzyaaak4,2)
pcharge(1:xyzzyaaal4)=-1.d0
pmass(1:xyzzyaaal4)=1.d0
inv_pmass(1:xyzzyaaal4)=1.d0
pspin(1)=0.5d0
pname(1)='Spin-up electron'
if(xyzzyaaak4>1)then
pspin(2)=-0.5d0
pname(2)='Spin-down electron'
popp_spin(1)=2
popp_spin(2)=1
else
popp_spin(1)=0
endif
pfermion(1:xyzzyaaal4)=.true.
pmasstensor(1,1,1:xyzzyaaal4)=1.d0
inv_pmasstensor(1,1,1:xyzzyaaal4)=1.d0
pmasstensor(2,2,1:xyzzyaaal4)=1.d0
inv_pmasstensor(2,2,1:xyzzyaaal4)=1.d0
pmasstensor(3,3,1:xyzzyaaal4)=1.d0
inv_pmasstensor(3,3,1:xyzzyaaal4)=1.d0
if(xyzzyaaak4>2.or.xyzzyaabc4)then
electron_system=.false.
xyzzyaaai4=3
if(xyzzyaabc4)xyzzyaaai4=1
do xyzzyaaah4=xyzzyaaai4,xyzzyaaak4
if(all(xyzzyaaau4(:,:,xyzzyaaah4)==0.d0))cycle
pcharge(xyzzyaaah4)=xyzzyaaav4(xyzzyaaah4)
pmass(xyzzyaaah4)=xyzzyaaaw4(xyzzyaaah4)
inv_pmass(xyzzyaaah4)=0.d0
if(pmass(xyzzyaaah4)>0.d0)inv_pmass(xyzzyaaah4)=1.d0/pmass(xyzzyaaah4)
pname(xyzzyaaah4)=xyzzyaabj4(xyzzyaaah4)
xyzzyaaaf4=nint(2*xyzzyaaay4(xyzzyaaah4))
pfermion(xyzzyaaah4)=.true.
if(mod(abs(xyzzyaaaf4),2)==0)pfermion(xyzzyaaah4)=.false.
pspin(xyzzyaaah4)=real(xyzzyaaaf4,dp)*0.5d0
pisotropic(xyzzyaaah4)=xyzzyaabg4(xyzzyaaah4)
pinfmass(xyzzyaaah4)=xyzzyaabi4(xyzzyaaah4)
pmasstensor(:,:,xyzzyaaah4)=xyzzyaaau4(:,:,xyzzyaaah4)
inv_pmasstensor(:,:,xyzzyaaah4)=inverse3(pmasstensor(:,:,xyzzyaaah4),x&
&yzzyaabe4)
if(xyzzyaabe4)call errstop_master('READ_PARTICLES','Mass tensor is sin&
&gular for particle type '//trim(i2s(xyzzyaaah4))//' in PARTICLES bloc&
&k in input.')
enddo
do xyzzyaaah4=1,xyzzyaaak4-1
do xyzzyaaai4=xyzzyaaah4+1,xyzzyaaak4
xyzzyaabb4=.true.
if(pcharge(xyzzyaaah4)/=pcharge(xyzzyaaai4))xyzzyaabb4=.false.
if(any(pmasstensor(:,:,xyzzyaaah4)/=pmasstensor(:,:,xyzzyaaai4)))xyzzy&
&aabb4=.false.
if(xyzzyaabb4.and.pspin(xyzzyaaah4)==-pspin(xyzzyaaai4))then
popp_spin(xyzzyaaah4)=xyzzyaaai4
popp_spin(xyzzyaaai4)=xyzzyaaah4
endif
if(pspin(xyzzyaaah4)/=pspin(xyzzyaaai4))xyzzyaabb4=.false.
if(xyzzyaabb4)call errstop_master('READ_PARTICLES','Two particle types&
& are the same. They should not be.')
enddo
enddo
if(.not.all(pisotropic))call errstop_master('READ_PARTICLES','Anisotro&
&pic mass tensors not yet implemented.')
endif
if(allocated(xyzzyaaaw4))deallocate(xyzzyaaav4,xyzzyaaaw4,xyzzyaaay4,x&
&yzzyaabj4,xyzzyaabg4,xyzzyaabi4,xyzzyaaau4)
if(esdf_block('free_particles',xyzzyaaaj4))then
if(xyzzyaaak4<nspin)call errstop_master('READ_PARTICLES','Fewer partic&
&le types than spins.')
allocate(heg_nele(xyzzyaaak4),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_nele')
xyzzyaadn1=0
do xyzzyaaah4=1,xyzzyaaaj4
read(block_data(xyzzyaaah4),*,iostat=xyzzyaadl1)temp_string,xyzzyaaaa4&
&,temp_det
if(xyzzyaadl1/=0)cycle
temp_string=adjustl(temp_string)
if(trim(temp_string)/='particle')cycle
temp_det=adjustl(temp_det)
if(trim(temp_det)=='det')then
if(xyzzyaadn1<0)call errstop_master('READ_PARTICLES','FREE_PARTICLES b&
&lock must contain either "particle i det j :" or "particle i :" const&
&ructs, not both.')
read(block_data(xyzzyaaah4),*,iostat=xyzzyaadl1)temp_string,xyzzyaaaa4&
&,temp_det,xyzzyaaac4
if(xyzzyaadl1/=0)call errstop_master('READ_PARTICLES','Problem reading&
& determinant index in FREE_PARTICLES block.')
xyzzyaadn1=max(xyzzyaaac4,xyzzyaadn1)
else
if(xyzzyaadn1>0)call errstop_master('READ_PARTICLES','FREE_PARTICLES b&
&lock must contain either "particle i det j :" or "particle i :" const&
&ructs, not both.')
xyzzyaadn1=-1
endif
enddo
if(xyzzyaadn1==0)call errstop_master('READ_PARTICLES','No orbital assi&
&gnments found in FREE_PARTICLES block.')
allocate(heg_orbtype(xyzzyaaak4,abs(xyzzyaadn1)),heg_slatt(xyzzyaaak4,&
&abs(xyzzyaadn1)),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_orbtype')
heg_orbtype=0
heg_nele=0
free_norb=0
r_s=-1.d0
dimensionality=0
xyzzyaabd4=.false.
heg_cell=0.d0
heg_cell(1,1)=5.d4
heg_cell(2,2)=5.d4
heg_cell(3,3)=5.d4
heg_crystal_type='unset'
heg_nbasis=0
heg_slatt=0
heg_repeat=0
pairing_wf=.false.
periodicity=0
sdw_magvec=0.d0
sdw_init_polarization=1.d0
xyzzyaabf4=.false.
xyzzyaaah4=0
do
xyzzyaaah4=xyzzyaaah4+1
if(xyzzyaaah4>xyzzyaaaj4)exit
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string
temp_string=adjustl(temp_string)
if(temp_string(1:3)=='r_s')then
if(r_s>0.d0)call errstop_master('READ_PARTICLES','Parameter "r_s" foun&
&d twice in FREE_PARTICLES block.')
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,r_s
if(r_s<=0.d0)call errstop_master('READ_PARTICLES','Bad r_s in FREE_PAR&
&TICLES block.')
elseif(temp_string(1:14)=='dimensionality')then
if(dimensionality/=0)call errstop_master('READ_PARTICLES','Parameter "&
&dimensionality" found twice in FREE_PARTICLES block.')
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,dimensionality
if((dimensionality<1.or.dimensionality>3))call errstop_master('SET_INP&
&UT_PARAMETER','Bad dimensionality in FREE_PARTICLES block. 1, 2 and 3&
& allowed.')
if(isperiodic)then
periodicity=dimensionality
else
xyzzyaabd4=.true.
endif
elseif(temp_string(1:11)=='periodicity')then
if(dimensionality==0)call errstop_master('READ_PARTICLES','Must specif&
&y dimensionality before periodicity in FREE_PARTICLES block.')
if(isperiodic.and.periodicity/=dimensionality)call errstop_master('REA&
&D_PARTICLES','Multiple conflicting definitions of periodicity in FREE&
&_PARTICLES block.')
if(r_s>0.d0)call errstop_master('READ_PARTICLES','Must specify r_s aft&
&er periodicity in FREE_PARTICLES block.')
if(xyzzyaabd4)call errstop_master('READ_PARTICLES','Cell geometry alre&
&ady defined when reading periodicity in FREE_PARTICLES block')
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,periodicity
if(isperiodic.and.periodicity==0)call errstop_master('READ_PARTICLES',&
&'Non-periodic systems must be defined by setting PERIODIC to F, rathe&
&r than by setting periodicity to 0 in the FREE_PARTICLES block.')
if(periodicity<0)call errstop_master('READ_PARTICLES','Bad periodicity&
& in FREE_PARTICLES block.')
if(periodicity>dimensionality)call errstop_master('READ_PARTICLES','Mu&
&st have 0 <= periodicity <= dimensionality in FREE_PARTICLES block.')
elseif(temp_string(1:13)=='cell_geometry')then
if(xyzzyaabd4)call errstop_master('READ_PARTICLES','Multiple or unnece&
&ssary definition for cell_geometry in FREE_PARTICLES block.')
if(dimensionality==0)call errstop_master('READ_PARTICLES','Need to def&
&ine dimensionality before cell_geometry in FREE_PARTICLES block.')
do xyzzyaaai4=1,periodicity
xyzzyaaah4=xyzzyaaah4+1
if(xyzzyaaah4>xyzzyaaaj4)call errstop_master('READ_PARTICLES','Reached&
& end of block while reading cell_geometry in FREE_PARTICLES.')
read(block_data(xyzzyaaah4),*,err=11,end=11)heg_cell(xyzzyaaai4,1:peri&
&odicity)
enddo
if(all(heg_cell(1:periodicity,1:periodicity)==0.d0))call errstop_maste&
&r('READ_PARTICLES','Invalid cell_geometry.')
xyzzyaabd4=.true.
elseif(temp_string(1:11)=='heg_nlayers')then
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,heg_nlayers
if(heg_nlayers<1)call errstop_master('READ_PARTICLES','Number of layer&
&s should be at least 1.')
deallocate(heg_ylayer,heg_zlayer)
allocate(heg_ylayer(heg_nlayers),heg_zlayer(heg_nlayers),stat=xyzzyaad&
&k1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_zlayer2')
heg_ylayer=0.d0
heg_zlayer=0.d0
elseif(temp_string(1:10)=='heg_ylayer')then
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,xyzzyaaai4,xyz&
&zyaaat4
if(xyzzyaaai4<1)call errstop_master('READ_PARTICLES','Layer index shou&
&ld be at least 1.')
if(xyzzyaaai4>heg_nlayers)call errstop_master('READ_PARTICLES','Layer &
&index should be at most the number of layers.')
heg_ylayer(xyzzyaaai4)=xyzzyaaat4
elseif(temp_string(1:10)=='heg_zlayer')then
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,xyzzyaaai4,xyz&
&zyaaat4
if(xyzzyaaai4<1)call errstop_master('READ_PARTICLES','Layer index shou&
&ld be at least 1.')
if(xyzzyaaai4>heg_nlayers)call errstop_master('READ_PARTICLES','Layer &
&index should be at most the number of layers.')
heg_zlayer(xyzzyaaai4)=xyzzyaaat4
elseif(temp_string(1:9)=='heg_layer')then
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,xyzzyaaai4,xyz&
&zyaaab4
if(xyzzyaaai4<1.or.xyzzyaaai4>xyzzyaaak4)call errstop_master('READ_PAR&
&TICLES','Particle type should be between 1 and the number of particle&
& types.')
if(xyzzyaaab4<1.or.xyzzyaaab4>heg_nlayers)call errstop_master('READ_PA&
&RTICLES','Layer should be between 1 and the number of layers.')
heg_layer(xyzzyaaai4)=xyzzyaaab4
elseif(temp_string(1:12)=='z-separation')then
if(heg_nlayers/=1)call errstop_master('READ_PARTICLES','Erroneous para&
&meter "z-separation" in FREE_PARTICLES block.')
if(dimensionality==0)call errstop_master('READ_PARTICLES','Need to def&
&ine dimensionality before z-separation in FREE_PARTICLES block.')
if(dimensionality==3)call errstop_master('READ_PARTICLES','z-separatio&
&n should only be given in 1D or 2D systems.')
heg_nlayers=2
deallocate(heg_ylayer,heg_zlayer)
allocate(heg_ylayer(heg_nlayers),heg_zlayer(heg_nlayers),stat=xyzzyaad&
&k1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_ylayer')
heg_ylayer=0.d0
heg_zlayer=0.d0
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,heg_zlayer(2)
if(heg_zlayer(2)<0.d0)call errstop_master('READ_PARTICLES','Negative z&
&-separation?')
elseif(temp_string(1:9)=='top-layer')then
if(heg_nlayers/=2)call errstop_master('READ_PARTICLES','Found "top-lay&
&er" in PARTICLES block, but number of layers is ' //trim(i2s(heg_nlay&
&ers))//'.')
temp_string=block_data(xyzzyaaah4)
temp_string=adjustl(temp_string(10:))
do
read(temp_string,*,iostat=xyzzyaadl1)xyzzyaaai4
if(xyzzyaadl1/=0)exit
if(xyzzyaaai4<1.or.xyzzyaaai4>xyzzyaaak4)call errstop_master('READ_PAR&
&TICLES','Bad particle number in "top-layer" assignment.')
heg_layer(xyzzyaaai4)=2
xyzzyaaae4=scan(temp_string,' ')
if(xyzzyaaae4<1)exit
temp_string=adjustl(temp_string(xyzzyaaae4:))
enddo
elseif(temp_string(1:16)=='crystal_type')then
if(trim(heg_crystal_type)/='unset')call errstop_master('READ_PARTICLES&
&','Parameter "crystal_type" found twice in FREE_PARTICLES block.')
if(.not.xyzzyaabd4)call errstop_master('READ_PARTICLES','Must define "&
&cell_geometry" before "crystal_type" in FREE_PARTICLES block.')
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,heg_crystal_ty&
&pe,heg_nbasis,temp_slatt
heg_crystal_type=adjustl(heg_crystal_type)
temp_slatt=adjustl(temp_slatt)
if(trim(temp_slatt(1:10))/='sublattice')call errstop_master('READ_PART&
&ICLES','Syntax: crystal_type <type> <n> sublattice[s] [repeat <r>].')
if(heg_nbasis<1)call errstop_master('READ_PARTICLES','Number of sublat&
&tices must be >0 .')
xyzzyaaae4=index(block_data(xyzzyaaah4),'sublattice')
temp_string=block_data(xyzzyaaah4)
temp_string=adjustl(temp_string(xyzzyaaae4:))
xyzzyaaae4=index(temp_string,' ')
if(xyzzyaaae4/=0)then
temp_string=temp_string(xyzzyaaae4:)
read(temp_string,*,iostat=xyzzyaadl1)temp_repeat,heg_repeat
if(xyzzyaadl1/=0)heg_repeat=1
endif
select case(trim(heg_crystal_type))
case('cubic','bcc','fcc','rectangular','hexagonal','triangular')
allocate(heg_wigner_basis(3,heg_nbasis),heg_ferro(heg_nbasis),stat=xyz&
&zyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_wigner_basis')
heg_wigner_basis=0.d0
heg_ferro=.true.
do xyzzyaaai4=1,heg_nbasis
xyzzyaaah4=xyzzyaaah4+1
if(xyzzyaaah4>xyzzyaaaj4)call errstop_master('READ_PARTICLES','Reached&
& end of block while reading crystal_type in FREE_PARTICLES.')
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_slatt,xyzzyaaab4,temp&
&_ferro
if(xyzzyaaab4/=xyzzyaaai4)call errstop_master('READ_PARTICLES','Need t&
&o supply sublattices in order in FREE_PARTICLES block.')
temp_ferro=adjustl(temp_ferro)
if(temp_ferro(1:9)=='antiferro')then
heg_ferro(xyzzyaaai4)=.false.
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_slatt,xyzzyaaab4,temp&
&_ferro,temp_off,heg_wigner_basis(1:3,xyzzyaaai4)
else
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_slatt,xyzzyaaab4,temp&
&_off,heg_wigner_basis(1:3,xyzzyaaai4)
endif
temp_slatt=adjustl(temp_slatt)
temp_off=adjustl(temp_off)
if(trim(temp_slatt)/='sublattice'.or.trim(temp_off)/='offset')call err&
&stop_master('READ_PARTICLES','Syntax: sublattice <i> [antiferro[magne&
&tic]] offset <x y z>')
enddo
case('manual')
if(heg_repeat/=1)call errstop_master('READ_PARTICLES','No "repeat" all&
&owed for manual crystal types in FREE_PARTICLES block.')
allocate(heg_nosites(heg_nbasis),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_nosites')
heg_nosites=0
xyzzyaaaa4=xyzzyaaah4+1
do xyzzyaaai4=1,heg_nbasis
if(xyzzyaaaa4>xyzzyaaaj4)call errstop_master('READ_PARTICLES','Reached&
& end of block while reading crystal_type in FREE_PARTICLES.')
temp_string=adjustl(block_data(xyzzyaaaa4))
xyzzyaaae4=index(temp_string,'manual')
if(xyzzyaaae4==0)call errstop_master('READ_PARTICLES','Syntax: sublatt&
&ice <i> manual <n> site[s]')
temp_string=adjustl(temp_string(xyzzyaaae4:))
xyzzyaaae4=scan(temp_string,' ')
read(temp_string,*,err=11,end=11)temp_ctype,xyzzyaaad4
if(xyzzyaaad4<1)call errstop_master('READ_PARTICLES','Empty sublattice&
&?')
heg_nosites(xyzzyaaai4)=xyzzyaaad4
xyzzyaaaa4=xyzzyaaaa4+xyzzyaaad4+1
enddo
allocate(heg_crystal_sites(3,maxval(heg_nosites),heg_nbasis),stat=xyzz&
&yaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_crystal_sites')
heg_crystal_sites=0.d0
do xyzzyaaai4=1,heg_nbasis
xyzzyaaah4=xyzzyaaah4+1
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_slatt,xyzzyaaab4,temp&
&_ctype,xyzzyaaad4,temp_sites
if(xyzzyaaab4/=xyzzyaaai4)call errstop_master('READ_PARTICLES','Need t&
&o supply lattice sites in order in FREE_PARTICLES block.')
temp_slatt=adjustl(temp_slatt)
temp_ctype=adjustl(temp_ctype)
temp_sites=adjustl(temp_sites)
if(trim(temp_slatt)/='sublattice'.or.trim(temp_ctype)/='manual'.or.tri&
&m(temp_sites(1:4))/='site')call errstop_master('READ_PARTICLES','Synt&
&ax: sublattice <i> manual <n> site[s]')
do xyzzyaaaa4=1,heg_nosites(xyzzyaaai4)
xyzzyaaah4=xyzzyaaah4+1
read(block_data(xyzzyaaah4),*,err=11,end=11)heg_crystal_sites(1:3,xyzz&
&yaaaa4,xyzzyaaai4)
enddo
enddo
case default
call errstop_master('READ_PARTICLES','Unrecognized crystalline type '/&
&/trim(heg_crystal_type)//'.')
end select
elseif(temp_string(1:10)=='sdw_magvec')then
if(xyzzyaabf4)call errstop_master('READ_PARTICLES','Parameter "sdw_mag&
&vec" found twice in FREE_PARTICLES block.')
xyzzyaabf4=.true.
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,sdw_magvec
if(.not.complex_wf)call errstop_master('READ_PARTICLES','Must have a c&
&omplex wave function for SDW calculations.')
elseif(temp_string(1:21)=='sdw_init_polarization')then
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,sdw_init_polar&
&ization
if(abs(sdw_init_polarization)>1.d0)call errstop_master('READ_PARTICLES&
&','SDW initial polarization must be between -1 and 1.')
elseif(temp_string(1:10)=='n_orbitals')then
if(free_norb>0)call errstop_master('READ_PARTICLES','Parameter "n_orbi&
&tals" found twice in FREE_PARTICLES block.')
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,free_norb
if(free_norb<=0)call errstop_master('READ_PARTICLES','Negative number &
&of orbitals read in FREE_PARTICLES block.')
elseif(temp_string(1:8)=='particle')then
if(xyzzyaadn1>0)then
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,xyzzyaaaa4,tem&
&p_det,xyzzyaaac4,temp_colon,xyzzyaaad4,temp_orb,temp_orbtype
temp_det=adjustl(temp_det)
temp_colon=adjustl(temp_colon)
temp_orb=adjustl(temp_orb)
temp_orbtype=adjustl(temp_orbtype)
if(trim(temp_det)/='det'.or.trim(temp_colon)/=':'.or.trim(temp_orb(1:7&
&))/='orbital')call errstop_master('READ_PARTICLES','Syntax: particle &
&<i> det <det> : <n> orbitals <orb_type> [orb_options]')
else
xyzzyaaac4=1
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,xyzzyaaaa4,tem&
&p_colon,xyzzyaaad4,temp_orb,temp_orbtype
temp_colon=adjustl(temp_colon)
temp_orb=adjustl(temp_orb)
temp_orbtype=adjustl(temp_orbtype)
if(trim(temp_colon)/=':'.or.trim(temp_orb(1:7))/='orbital')call errsto&
&p_master('READ_PARTICLES','Syntax: particle <i> : <n> orbitals <orb_t&
&ype> [orb_options]')
endif
if(xyzzyaaaa4<1.or.xyzzyaaaa4>xyzzyaaak4)call errstop_master('READ_PAR&
&TICLES','Particle type index not valid in FREE_PARTICLES block.')
if(xyzzyaaac4<1)call errstop_master('READ_PARTICLES','Bad multidetermi&
&nant index in FREE_PARTICLES block.')
if(xyzzyaaac4>abs(xyzzyaadn1))call errstop_master('READ_PARTICLES','Wo&
&uld have sworn there were '//trim(i2s(abs(xyzzyaadn1)))//' MD terms, &
&but a '//trim(i2s(xyzzyaaac4))//'th term was found. Terribly confused&
& (bug).')
if(heg_orbtype(xyzzyaaaa4,xyzzyaaac4)/=0)call errstop_master('READ_PAR&
&TICLES','Particle has been assigned orbitals twice in FREE_PARTICLES &
&block.')
if(xyzzyaaad4<0)call errstop_master('READ_PARTICLES','Negative number &
&of particles/orbitals?')
select case(trim(temp_orbtype))
case('free')
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=1
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants.')
heg_nele(xyzzyaaaa4)=xyzzyaaad4
case('crystal')
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=2
xyzzyaaae4=index(block_data(xyzzyaaah4),'crystal')
if(xyzzyaaae4<1)call errstop_master('READ_PARTICLES','Bug?')
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants.')
heg_nele(xyzzyaaaa4)=xyzzyaaad4
temp_string=block_data(xyzzyaaah4)
temp_string=adjustl(temp_string(xyzzyaaae4:))
read(temp_string,*,err=11,end=11)temp_orbtype,temp_slatt,xyzzyaaad4
if(xyzzyaaad4>heg_nbasis)call errstop_master('READ_PARTICLES','Sublatt&
&ice '//trim(i2s(xyzzyaaad4))//' must be defined before being assigned&
&.')
heg_slatt(xyzzyaaaa4,xyzzyaaac4)=xyzzyaaad4
case('pairing')
xyzzyaaae4=index(block_data(xyzzyaaah4),'pairing')
if(xyzzyaaae4<1)call errstop_master('READ_PARTICLES','Bug?')
temp_string=block_data(xyzzyaaah4)
temp_string=adjustl(temp_string(xyzzyaaae4:))
read(temp_string,*,err=11,end=11)temp_orbtype,xyzzyaaab4
if(xyzzyaaaa4==xyzzyaaab4)call errstop_master('READ_PARTICLES','Cannot&
& pair a particle type with itself')
if(heg_orbtype(xyzzyaaab4,xyzzyaaac4)/=0)call errstop_master('READ_PAR&
&TICLES','Particle to pair with has already been assigned an orbital i&
&n FREE_PARTICLES block.')
xyzzyaaai4=max(xyzzyaaaa4,xyzzyaaab4)
xyzzyaaaa4=min(xyzzyaaaa4,xyzzyaaab4)
xyzzyaaab4=xyzzyaaai4
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants <1>.')
if(heg_nele(xyzzyaaab4)/=0.and.heg_nele(xyzzyaaab4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants <2>.')
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=-xyzzyaaab4
heg_orbtype(xyzzyaaab4,xyzzyaaac4)=-xyzzyaaaa4
pairing_wf=.true.
heg_nele(xyzzyaaaa4)=xyzzyaaad4
heg_nele(xyzzyaaab4)=xyzzyaaad4
case('sdw')
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=3
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants.')
heg_nele(xyzzyaaaa4)=xyzzyaaad4
noncoll_spin=.true.
case('expot')
xyzzyaaae4=index(block_data(xyzzyaaah4),'expot')
if(xyzzyaaae4<1)call errstop_master('READ_PARTICLES','Bug?')
temp_string=block_data(xyzzyaaah4)
temp_string=adjustl(temp_string(xyzzyaaae4:))
read(temp_string,*,err=11,end=11)temp_orbtype,xyzzyaaab4
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=100+xyzzyaaab4
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants.')
heg_nele(xyzzyaaaa4)=xyzzyaaad4
case('biex1')
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=4
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants.')
heg_nele(xyzzyaaaa4)=xyzzyaaad4
case('biex2')
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=5
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants.')
heg_nele(xyzzyaaaa4)=xyzzyaaad4
case('biex3')
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=6
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants.')
heg_nele(xyzzyaaaa4)=xyzzyaaad4
case('exmol')
heg_orbtype(xyzzyaaaa4,xyzzyaaac4)=7
pairing_wf=.true.
if(heg_nele(xyzzyaaaa4)/=0.and.heg_nele(xyzzyaaaa4)/=xyzzyaaad4)call e&
&rrstop_master('READ_PARTICLES','Particle-number mismatch between dete&
&rminants.')
heg_nele(xyzzyaaaa4)=xyzzyaaad4
case default
call errstop_master('READ_PARTICLES','Bad orbital name "'//trim(temp_o&
&rbtype)//'".')
end select
elseif(temp_string(1:8)=='k_offset')then
if(periodicity==0)call errstop_master('READ_PARTICLES','Need to define&
& dimensionality before k_offset in FREE_PARTICLES block.')
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,k_offset(1:per&
&iodicity)
if(.not.complex_wf.and.any(k_offset/=0.d0))call errstop_master('READ_P&
&ARTICLES','K-vector offset must be zero if a real wave function is us&
&ed.')
elseif(temp_string(1:8)=='me_biex3')then
if(fix_holes)then
call errstop_master('READ_PARTICLES','Do not need to specify me_biex3 &
&for FIX_HOLES=T; the masses in the PARTICLES block are the masses of &
&the electrons.')
else
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,me_biex3
if(me_biex3<=0.d0)call errstop_master('READ_PARTICLES','Electron mass &
&me_biex3 should be positive.')
endif
elseif(temp_string(1:8)=='mh_biex3')then
if(fix_holes)then
call errstop_master('READ_PARTICLES','Hole mass should not be specifie&
&d when using the fixed hole constraint.')
else
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,mh_biex3
if(mh_biex3<=0.d0)call errstop_master('READ_PARTICLES','Hole mass mh_b&
&iex3 should be positive.')
endif
elseif(temp_string(1:8)=='xx_sep')then
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,xx_sep
if(xx_sep<0.d0)call errstop_master('READ_PARTICLES','Exciton-exciton s&
&eparation should be positive.')
elseif(temp_string(1:8)=='quasi_1d')then
if(dimensionality/=1)then
call errstop_master('READ_PARTICLES','Cannot specify QUASI_1D, the con&
&finement strength for a 1d wire, when the dimensionality is not 1.')
else
read(block_data(xyzzyaaah4),*,err=11,end=11)temp_string,harmwire_b
endif
else
call errstop_master('READ_PARTICLES','Unrecognized keyword '//trim(tem&
&p_string)//' in FREE_PARTICLES block.')
endif
enddo
if(trim(atom_basis_type)=='none')then
if(dimensionality==0.or.(isperiodic.and.(r_s<0.d0.or..not.xyzzyaabd4))&
&)call errstop_master('READ_PARTICLES','Geometry for atom-less calcula&
&tion not fully given.')
if(all(heg_orbtype/=6))then
do xyzzyaaam4=1,heg_nlayers
if(all(heg_layer(:)/=xyzzyaaam4))call errwarn('READ_PARTICLES','No par&
&ticles are assigned to layer '//trim(i2s(xyzzyaaam4))//'.')
enddo
endif
if(any(heg_layer<1.or.heg_layer>heg_nlayers))call errstop('READ_PARTIC&
&LES','Layer out of range.')
endif
if(dimensionality<=0)dimensionality=3
if(.not.allocated(heg_orbtype))call errstop_master('READ_PARTICLES','N&
&o free particles/orbitals defined <1>.')
if(all(heg_orbtype==0))call errstop_master('READ_PARTICLES','No free p&
&articles/orbitals defined <2> - bug?')
do xyzzyaaac4=1,abs(xyzzyaadn1)
do xyzzyaaah4=1,xyzzyaaak4
if(heg_orbtype(xyzzyaaah4,xyzzyaaac4)==0.and.heg_nele(xyzzyaaah4)>0)ca&
&ll errstop_master('READ_PARTICLES','Particle-number mismatch between &
&determinants <3>.')
enddo
enddo
if(any(heg_orbtype==3).and..not.xyzzyaabf4)call errstop_master('READ_P&
&ARTICLES','Magnetization wave vector not specified for spin density w&
&ave.')
nspin=xyzzyaaak4
if(any(heg_orbtype==2).or.any(heg_orbtype<0).or.any(heg_orbtype==3))co&
&rr_heg_required=.true.
else
if(xyzzyaaai1>0)then
nspin=4
elseif(xyzzyaaah1>0)then
nspin=3
else
nspin=2
endif
pairing_wf=.false.
allocate(heg_nele(nspin),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','heg_nele <2>')
heg_nele=0
endif
if(model_system)then
xyzzyaaah4=2
xyzzyaaag4=2
do
xyzzyaaah4=xyzzyaaah4+1
if(xyzzyaaah4>nspin)exit
xyzzyaaag4=xyzzyaaag4+1
if(heg_nele(xyzzyaaah4)/=0.and.any(heg_orbtype(xyzzyaaah4,:)/=0))cycle
call errwarn('READ_PARTICLES','user-defined particle type '//trim(i2s(&
&xyzzyaaag4))//' unused. Ignoring.')
do xyzzyaaai4=xyzzyaaah4,nspin-1
pcharge(xyzzyaaai4)=pcharge(xyzzyaaai4+1)
pmass(xyzzyaaai4)=pmass(xyzzyaaai4+1)
inv_pmass(xyzzyaaai4)=inv_pmass(xyzzyaaai4+1)
pname(xyzzyaaai4)=pname(xyzzyaaai4+1)
pfermion(xyzzyaaai4)=pfermion(xyzzyaaai4+1)
pspin(xyzzyaaai4)=pspin(xyzzyaaai4+1)
popp_spin(xyzzyaaai4)=popp_spin(xyzzyaaai4+1)
pisotropic(xyzzyaaai4)=pisotropic(xyzzyaaai4+1)
pinfmass(xyzzyaaai4)=pinfmass(xyzzyaaai4+1)
pmasstensor(:,:,xyzzyaaai4)=pmasstensor(:,:,xyzzyaaai4+1)
inv_pmasstensor(:,:,xyzzyaaai4)=inv_pmasstensor(:,:,xyzzyaaai4+1)
heg_nele(xyzzyaaai4)=heg_nele(xyzzyaaai4+1)
heg_orbtype(xyzzyaaai4,:)=heg_orbtype(xyzzyaaai4+1,:)
heg_layer(xyzzyaaai4)=heg_layer(xyzzyaaai4+1)
heg_slatt(xyzzyaaai4,:)=heg_slatt(xyzzyaaai4+1,:)
enddo
where(popp_spin(:)==xyzzyaaah4)popp_spin(:)=0
xyzzyaaah4=xyzzyaaah4-1
nspin=nspin-1
enddo
endif
if(nspin/=xyzzyaaak4)then
if(.not.model_system)call errstop('READ_PARTICLES','Please define your&
& holes in a PARTICLES block in input.')
xyzzyaaak4=nspin
allocate(xyzzyaaav4(nspin),xyzzyaaaw4(nspin),xyzzyaaax4(nspin),xyzzyaa&
&bj4(nspin),xyzzyaabh4(nspin),xyzzyaaay4(nspin),xyzzyaaao4(nspin),xyzz&
&yaabg4(nspin),xyzzyaaau4(3,3,nspin),xyzzyaaaz4(3,3,nspin),xyzzyaaap4(&
&nspin),xyzzyaaaq4(nspin,abs(xyzzyaadn1)),xyzzyaaas4(nspin),xyzzyaaar4&
&(nspin,abs(xyzzyaadn1)),xyzzyaabi4(nspin),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','resized particle propert&
&ies <1>')
xyzzyaaav4(1:nspin)=pcharge(1:nspin)
xyzzyaaaw4(1:nspin)=pmass(1:nspin)
xyzzyaaax4(1:nspin)=inv_pmass(1:nspin)
xyzzyaabj4(1:nspin)=pname(1:nspin)
xyzzyaabh4(1:nspin)=pfermion(1:nspin)
xyzzyaaay4(1:nspin)=pspin(1:nspin)
xyzzyaaao4(1:nspin)=popp_spin(1:nspin)
xyzzyaabg4(1:nspin)=pisotropic(1:nspin)
xyzzyaabi4(1:nspin)=pinfmass(1:nspin)
xyzzyaaau4(1:3,1:3,1:nspin)=pmasstensor(1:3,1:3,1:nspin)
xyzzyaaaz4(1:3,1:3,1:nspin)=inv_pmasstensor(1:3,1:3,1:nspin)
xyzzyaaap4(1:nspin)=heg_nele(1:nspin)
xyzzyaaaq4(1:nspin,1:abs(xyzzyaadn1))=heg_orbtype(1:nspin,1:abs(xyzzya&
&adn1))
xyzzyaaas4(1:nspin)=heg_layer(1:nspin)
xyzzyaaar4(1:nspin,1:abs(xyzzyaadn1))=heg_slatt(1:nspin,1:abs(xyzzyaad&
&n1))
deallocate(pcharge,pmass,inv_pmass,pname,pfermion,pspin,popp_spin,piso&
&tropic,pmasstensor,inv_pmasstensor,heg_nele,heg_orbtype,heg_layer,heg&
&_slatt,pinfmass)
allocate(pcharge(nspin),pmass(nspin),inv_pmass(nspin),pname(nspin),pfe&
&rmion(nspin),pspin(nspin),popp_spin(nspin),pisotropic(nspin),pmassten&
&sor(3,3,nspin),inv_pmasstensor(3,3,nspin),heg_nele(nspin),heg_orbtype&
&(nspin,abs(xyzzyaadn1)),heg_layer(nspin),heg_slatt(nspin,abs(xyzzyaad&
&n1)),pinfmass(nspin),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','resized particle propert&
&ies <2>')
pcharge(1:nspin)=xyzzyaaav4(1:nspin)
pmass(1:nspin)=xyzzyaaaw4(1:nspin)
inv_pmass(1:nspin)=xyzzyaaax4(1:nspin)
pname(1:nspin)=xyzzyaabj4(1:nspin)
pfermion(1:nspin)=xyzzyaabh4(1:nspin)
pspin(1:nspin)=xyzzyaaay4(1:nspin)
popp_spin(1:nspin)=xyzzyaaao4(1:nspin)
pisotropic(1:nspin)=xyzzyaabg4(1:nspin)
pinfmass(1:nspin)=xyzzyaabi4(1:nspin)
pmasstensor(1:3,1:3,1:nspin)=xyzzyaaau4(1:3,1:3,1:nspin)
inv_pmasstensor(1:3,1:3,1:nspin)=xyzzyaaaz4(1:3,1:3,1:nspin)
heg_nele(1:nspin)=xyzzyaaap4(1:nspin)
heg_orbtype(1:nspin,1:abs(xyzzyaadn1))=xyzzyaaaq4(1:nspin,1:abs(xyzzya&
&adn1))
heg_layer(1:nspin)=xyzzyaaas4(1:nspin)
heg_slatt(1:nspin,1:abs(xyzzyaadn1))=xyzzyaaar4(1:nspin,1:abs(xyzzyaad&
&n1))
deallocate(xyzzyaaav4,xyzzyaaaw4,xyzzyaaax4,xyzzyaabj4,xyzzyaabh4,xyzz&
&yaaay4,xyzzyaaao4,xyzzyaabg4,xyzzyaaau4,xyzzyaaaz4,xyzzyaaap4,xyzzyaa&
&aq4,xyzzyaaas4,xyzzyaaar4,xyzzyaabi4)
endif
allocate(nele(nspin),nuc_nele(nspin),three_nele(nspin),six_nele(nspin)&
&,stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PARTICLES','nele array')
nuc_nele=0
nuc_nele(1)=xyzzyaaaf1
if(nspin>1)then
nuc_nele(2)=xyzzyaaag1
if(nspin>2)then
nuc_nele(3)=xyzzyaaah1
if(nspin>3)nuc_nele(4)=xyzzyaaai1
endif
endif
nele=nuc_nele+heg_nele
three_nele=3*nele
six_nele=6*nele
ferromagnetic=.true.
do xyzzyaaah4=1,nspin-1
if(nele(xyzzyaaah4)>0.and.any(nele(xyzzyaaah4+1:nspin)>0))then
ferromagnetic=.false.
exit
endif
enddo
if(esdf_block('fixed_particles',xyzzyaaaj4))then
nbasis=xyzzyaaaj4
allocate(basis(3,nbasis),atno(nbasis),stat=xyzzyaadk1)
if(xyzzyaadk1/=0)call errstop('READ_PARTICLES','Allocation error: rion&
&.')
do xyzzyaaah4=1,nbasis
read(block_data(xyzzyaaah4),*,err=11,end=11)atno(xyzzyaaah4),basis(1:3&
&,xyzzyaaah4)
enddo
endif
if(am_master)then
do xyzzyaaam4=1,heg_nlayers-1
if(.not.any(heg_layer(:)==xyzzyaaam4.and.nele(:)>0))cycle
do xyzzyaaan4=xyzzyaaam4+1,heg_nlayers
if(.not.any(heg_layer(:)==xyzzyaaan4.and.nele(:)>0))cycle
if(abs(heg_ylayer(xyzzyaaam4)-heg_ylayer(xyzzyaaan4))<xyzzyaaba4 .and.&
&abs(heg_zlayer(xyzzyaaam4)-heg_zlayer(xyzzyaaan4))<xyzzyaaba4)call er&
&rwarn('READ_PARTICLES','Layers '//trim(i2s(xyzzyaaam4))//' and ' //tr&
&im(i2s(xyzzyaaan4))//' have zero separation.  The Kato cusp condition&
&s will not be enforced correctly.  I hope you know what you are doing&
&...')
enddo
enddo
endif
return
10 call errstop_master('READ_PARTICLES','Problem reading PARTICLES blo&
&ck in input.')
11 call errstop_master('READ_PARTICLES','Problem reading FREE_PARTICLE&
&S block in input.')
end subroutine xyzzyaaex1
subroutine xyzzyaaey1
implicit none
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5
logical xyzzyaaag5
character(80) temp1,temp2,temp3,temp4,temp5,char_80,append_char,comb_c&
&har,tmpr
call wout('Particles')
call wout('=========')
call wout('Particle name                 Charge        Mass         Sp&
&in   Type')
call wout(repeat('-',73))
do xyzzyaaaa5=1,nspin
temp1=trim(i2s(xyzzyaaaa5))//': '//trim(adjustl(pname(xyzzyaaaa5)))
temp2=trim(r2s(pcharge(xyzzyaaaa5),'(f10.5)'))
temp3=trim(r2s(pmass(xyzzyaaaa5),'(f14.7)'))
if(.not.pisotropic(xyzzyaaaa5))then
temp3='Anisotropic'
elseif(pinfmass(xyzzyaaaa5))then
temp3='Infinity'
endif
temp4=trim(r2s(pspin(xyzzyaaaa5),'(f6.1)'))
temp5='Fermion'
if(.not.pfermion(xyzzyaaaa5))temp5='Boson'
write(tmpr,'(a,t28,a10,t40,a14,t56,a6,t65,a)')trim(temp1(1:26)),trim(t&
&emp2),trim(temp3),trim(temp4),trim(temp5)
call wout(tmpr)
enddo
call wout(repeat('-',73))
call wout()
call wout('Number of diffusion types : '//trim(i2s(no_difftypes)))
call wout()
call wout('Single-particle groupings')
call wout('-------------------------')
do xyzzyaaaa5=-custom_ssingles,3
if(xyzzyaaaa5>0)then
xyzzyaaaf5=xyzzyaado1(xyzzyaaaa5)
if(xyzzyaaaf5==0)cycle
elseif(xyzzyaaaa5<0)then
xyzzyaaaf5=xyzzyaaaa5
else
xyzzyaaaf5=0
endif
if(xyzzyaaaf5<0)then
char_80=' Spin dep. '//trim(i2s(xyzzyaaaf5))//':'
else
char_80=' Spin dep. '//trim(i2s(xyzzyaaaf5))//' :'
endif
do xyzzyaaae5=1,no_ssingles(xyzzyaaaf5)
xyzzyaaag5=.true.
do xyzzyaaab5=1,nspin
if(which_ssingle(xyzzyaaab5,xyzzyaaaf5)==xyzzyaaae5)then
if(xyzzyaaag5)then
if(len_trim(char_80)+len_trim(i2s(xyzzyaaab5))+3>80)then
call wout(trim(char_80),fmt='(a)')
char_80='              ('//trim(i2s(xyzzyaaab5))
else
char_80=trim(char_80)//' ('//trim(i2s(xyzzyaaab5))
endif
xyzzyaaag5=.false.
else
if(len_trim(char_80)+len_trim(i2s(xyzzyaaab5))+2>80)then
call wout(trim(char_80)//',',fmt='(a)')
char_80='               '//trim(i2s(xyzzyaaab5))
else
char_80=trim(char_80)//','//trim(i2s(xyzzyaaab5))
endif
endif
endif
enddo
char_80=trim(char_80)//')'
enddo
append_char=''
if(xyzzyaaaf5==level_family)append_char=trim(append_char)//' [F]'
if(xyzzyaaaf5==level_eqvfamily)append_char=trim(append_char)//' [E]'
if(xyzzyaaaa5<0)append_char=trim(append_char)//' [custom]'
if(len_trim(char_80)+len_trim(append_char)>80)then
call wout(trim(char_80),fmt='(a)')
call wout(trim(append_char),fmt='(a)')
else
call wout(trim(char_80)//trim(append_char),fmt='(a)')
endif
enddo
call wout()
call wout('NB, partitions defining spin-[F]amilies and [E]quivalent pa&
&rticles flagged.')
call wout()
call wout('Particle-pair groupings')
call wout('-----------------------')
do xyzzyaaaa5=-custom_spairs,5
if(xyzzyaaaa5>0)then
xyzzyaaaf5=xyzzyaadp1(xyzzyaaaa5)
if(xyzzyaaaf5==0)cycle
elseif(xyzzyaaaa5<0)then
xyzzyaaaf5=xyzzyaaaa5
else
xyzzyaaaf5=0
endif
if(xyzzyaaaf5<0)then
char_80=' Spin-pair dep. '//trim(i2s(xyzzyaaaf5))//':'
else
char_80=' Spin-pair dep. '//trim(i2s(xyzzyaaaf5))//' :'
endif
do xyzzyaaae5=1,no_spairs(xyzzyaaaf5)
xyzzyaaag5=.true.
do xyzzyaaab5=1,nspin
do xyzzyaaac5=xyzzyaaab5,nspin
if(which_spair(xyzzyaaab5,xyzzyaaac5,xyzzyaaaf5)==xyzzyaaae5)then
comb_char=trim(i2s(xyzzyaaab5))//'-'//trim(i2s(xyzzyaaac5))
if(xyzzyaaag5)then
if(len_trim(char_80)+len_trim(comb_char)+3>80)then
call wout(trim(char_80),fmt='(a)')
char_80='                   ('//trim(comb_char)
else
char_80=trim(char_80)//' ('//trim(comb_char)
endif
xyzzyaaag5=.false.
else
if(len_trim(char_80)+len_trim(comb_char)+2>80)then
call wout(trim(char_80)//',',fmt='(a)')
char_80='                    '//trim(comb_char)
else
char_80=trim(char_80)//','//trim(comb_char)
endif
endif
endif
enddo
enddo
char_80=trim(char_80)//')'
enddo
append_char=''
if(xyzzyaaaa5<0)append_char=' [custom]'
if(len_trim(char_80)+len_trim(append_char)>80)then
call wout(trim(char_80),fmt='(a)')
call wout(trim(append_char),fmt='(a)')
else
call wout(trim(char_80)//trim(append_char),fmt='(a)')
endif
enddo
call wout()
if(custom_striplets>0)then
call wout('Particle-triplet groupings')
call wout('--------------------------')
do xyzzyaaaa5=-custom_striplets,1
if(xyzzyaaaa5>0)then
xyzzyaaaf5=xyzzyaadq1(xyzzyaaaa5)
if(xyzzyaaaf5==0)cycle
elseif(xyzzyaaaa5<0)then
xyzzyaaaf5=xyzzyaaaa5
else
xyzzyaaaf5=0
endif
if(xyzzyaaaf5<0)then
char_80=' Spin-triplet dep. '//trim(i2s(xyzzyaaaf5))//':'
else
char_80=' Spin-triplet dep. '//trim(i2s(xyzzyaaaf5))//' :'
endif
do xyzzyaaae5=1,no_striplets(xyzzyaaaf5)
xyzzyaaag5=.true.
do xyzzyaaab5=1,nspin
do xyzzyaaac5=xyzzyaaab5,nspin
do xyzzyaaad5=xyzzyaaac5,nspin
if(which_striplet(xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaaf5)==xyzzya&
&aae5)then
select case(eq_triplet(xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaaf5))
case(1)
comb_char=trim(i2s(xyzzyaaab5))//'='//trim(i2s(xyzzyaaac5))//'-'//trim&
&(i2s(xyzzyaaad5))
case(2)
comb_char=trim(i2s(xyzzyaaab5))//'-'//trim(i2s(xyzzyaaac5))//'='//trim&
&(i2s(xyzzyaaad5))
case(4)
comb_char=trim(i2s(xyzzyaaab5))//'='//trim(i2s(xyzzyaaad5))//'-'//trim&
&(i2s(xyzzyaaac5))
case(7)
comb_char=trim(i2s(xyzzyaaab5))//'='//trim(i2s(xyzzyaaac5))//'='//trim&
&(i2s(xyzzyaaad5))
case default
comb_char=trim(i2s(xyzzyaaab5))//'-'//trim(i2s(xyzzyaaac5))//'-'//trim&
&(i2s(xyzzyaaad5))
end select
if(xyzzyaaag5)then
if(len_trim(char_80)+len_trim(comb_char)+3>80)then
call wout(trim(char_80),fmt='(a)')
char_80='                       ('//trim(comb_char)
else
char_80=trim(char_80)//' ('//trim(comb_char)
endif
xyzzyaaag5=.false.
else
if(len_trim(char_80)+len_trim(comb_char)+2>80)then
call wout(trim(char_80)//',',fmt='(a)')
char_80='                        '//trim(comb_char)
else
char_80=trim(char_80)//','//trim(comb_char)
endif
endif
endif
enddo
enddo
enddo
char_80=trim(char_80)//')'
enddo
append_char=''
if(xyzzyaaaa5<0)append_char=' [custom]'
if(len_trim(char_80)+len_trim(append_char)>80)then
call wout(trim(char_80),fmt='(a)')
call wout(trim(append_char),fmt='(a)')
else
call wout(trim(char_80)//trim(append_char),fmt='(a)')
endif
enddo
call wout()
endif
end subroutine xyzzyaaey1
subroutine xyzzyaaez1
implicit none
integer xyzzyaaaa6,xyzzyaaab6
allocate(missing_det(nspin,ndet),update_by_column(nspin,ndet),upd_spin&
&(nspin,ndet),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'FLAG_MISSING_DETS','Missing determinant.'&
&)
missing_det=.false.
update_by_column=.true.
upd_spin=0
do xyzzyaaaa6=1,nspin
upd_spin(xyzzyaaaa6,:)=xyzzyaaaa6
do xyzzyaaab6=1,ndet
if(nele(xyzzyaaaa6)==0)missing_det(xyzzyaaaa6,xyzzyaaab6)=.true.
if(heg_orbtype(xyzzyaaaa6,xyzzyaaab6)<0.and.xyzzyaaaa6>-heg_orbtype(xy&
&zzyaaaa6,xyzzyaaab6))then
missing_det(xyzzyaaaa6,xyzzyaaab6)=.true.
update_by_column(xyzzyaaaa6,xyzzyaaab6)=.false.
upd_spin(xyzzyaaaa6,xyzzyaaab6)=-heg_orbtype(xyzzyaaaa6,xyzzyaaab6)
endif
enddo
enddo
end subroutine xyzzyaaez1
subroutine xyzzyaafa1
implicit none
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7,xyzzyaa&
&af7,xyzzyaaag7,xyzzyaaah7,xyzzyaaai7,xyzzyaaaj7,xyzzyaaak7,xyzzyaaal7&
&,xyzzyaaam7,xyzzyaaan7,xyzzyaaao7,xyzzyaaap7
integer,allocatable :: xyzzyaaaq7(:),xyzzyaaar7(:)
logical xyzzyaaas7,xyzzyaaat7,xyzzyaaau7,xyzzyaaav7,xyzzyaaaw7,xyzzyaa&
&ax7,xyzzyaaay7,xyzzyaaaz7,xyzzyaaba7,xyzzyaabb7,xyzzyaabc7,xyzzyaabd7
logical,allocatable :: xyzzyaabe7(:),xyzzyaabf7(:),xyzzyaabg7(:)
max_spin_singles=nspin
max_spin_pairs=(nspin*(nspin+1))/2
max_spin_triplets=(nspin*(nspin*(nspin+3)+2))/6
allocate(xyzzyaaaq7(nspin),xyzzyaaar7(nspin),which_difftype(nspin),sta&
&t=xyzzyaaah7)
call check_alloc(xyzzyaaah7,'ASSIGN_SPIN_DEPS','')
xyzzyaaaq7=0
xyzzyaaai7=0
xyzzyaaar7=0
xyzzyaaaj7=0
which_difftype=0
no_difftypes=0
do xyzzyaaaa7=1,nspin
do xyzzyaaab7=1,xyzzyaaaa7-1
xyzzyaaas7=pmass(xyzzyaaab7)==pmass(xyzzyaaaa7)
xyzzyaabd7=nele(xyzzyaaab7)==nele(xyzzyaaaa7).or.nele(xyzzyaaab7)==0.o&
&r.nele(xyzzyaaaa7)==0
xyzzyaaat7=pcharge(xyzzyaaab7)==pcharge(xyzzyaaaa7)
xyzzyaaau7=pcharge(xyzzyaaab7)==-pcharge(xyzzyaaaa7)
if(pfermion(xyzzyaaaa7).eqv.pfermion(xyzzyaaab7))then
if(xyzzyaaas7)then
if(xyzzyaaat7)then
xyzzyaaar7(xyzzyaaaa7)=xyzzyaaar7(xyzzyaaab7)
xyzzyaaaq7(xyzzyaaaa7)=xyzzyaaaq7(xyzzyaaab7)
if(xyzzyaabd7)which_difftype(xyzzyaaaa7)=which_difftype(xyzzyaaab7)
elseif(xyzzyaaau7)then
xyzzyaaar7(xyzzyaaaa7)=xyzzyaaar7(xyzzyaaab7)
if(xyzzyaabd7.and.nbasis==0.and.trim(atom_basis_type)=='none')which_di&
&fftype(xyzzyaaaa7)=which_difftype(xyzzyaaab7)
endif
endif
endif
enddo
if(xyzzyaaaq7(xyzzyaaaa7)==0)then
xyzzyaaai7=xyzzyaaai7+1
xyzzyaaaq7(xyzzyaaaa7)=xyzzyaaai7
endif
if(xyzzyaaar7(xyzzyaaaa7)==0)then
xyzzyaaaj7=xyzzyaaaj7+1
xyzzyaaar7(xyzzyaaaa7)=xyzzyaaaj7
endif
if(which_difftype(xyzzyaaaa7)==0)then
no_difftypes=no_difftypes+1
which_difftype(xyzzyaaaa7)=no_difftypes
endif
enddo
allocate(xyzzyaabe7(5),xyzzyaabf7(3),xyzzyaabg7(1),stat=xyzzyaaah7)
call check_alloc(xyzzyaaah7,'ASSIGN_SPIN_DEPS','')
xyzzyaabg7(1)=(nspin>1)
xyzzyaabe7(1)=(xyzzyaaaj7<nspin.and.xyzzyaaaj7>1)
xyzzyaabe7(2)=(xyzzyaaai7<nspin)
xyzzyaabe7(3)=(xyzzyaaai7<nspin.and.xyzzyaaai7>1)
xyzzyaabe7(4)=(xyzzyaaai7<nspin.and.xyzzyaaai7>xyzzyaaaj7)
xyzzyaabe7(5)=(nspin>1)
xyzzyaabf7(1)=(xyzzyaaaj7<nspin.and.xyzzyaaaj7>1.and.xyzzyaaaj7/=xyzzy&
&aaai7)
xyzzyaabf7(2)=(xyzzyaaai7<nspin.and.xyzzyaaai7>1)
xyzzyaabf7(3)=(nspin>1)
levels_striplets=count(xyzzyaabg7)
levels_spairs=count(xyzzyaabe7)
levels_ssingles=count(xyzzyaabf7)
allocate(which_ssingle(nspin,-custom_ssingles:levels_ssingles),which_s&
&pair(nspin,nspin,-custom_spairs:levels_spairs),which_striplet(nspin,n&
&spin,nspin,-custom_striplets:levels_striplets),eq_triplet(nspin,nspin&
&,nspin,-custom_striplets:levels_striplets),no_ssingles(-custom_ssingl&
&es:levels_ssingles),no_spairs(-custom_spairs:levels_spairs),no_stripl&
&ets(-custom_striplets:levels_striplets),stat=xyzzyaaah7)
call check_alloc(xyzzyaaah7,'ASSIGN_SPIN_DEPS','')
which_ssingle=0
which_spair=0
which_striplet=0
which_ssingle(:,0)=1
which_spair(:,:,0)=1
which_striplet(:,:,:,0)=1
no_ssingles=0
no_spairs=0
no_striplets=0
no_ssingles(0)=1
no_spairs(0)=1
no_striplets(0)=1
eq_triplet=0
eq_triplet(:,:,:,0)=7
xyzzyaaae7=0
do xyzzyaaaa7=1,1
if(.not.xyzzyaabg7(xyzzyaaaa7))cycle
xyzzyaaae7=xyzzyaaae7+1
xyzzyaadq1(xyzzyaaaa7)=xyzzyaaae7
enddo
xyzzyaaae7=0
xyzzyaadp1=0
do xyzzyaaaa7=1,5
if(.not.xyzzyaabe7(xyzzyaaaa7))cycle
xyzzyaaae7=xyzzyaaae7+1
xyzzyaadp1(xyzzyaaaa7)=xyzzyaaae7
enddo
xyzzyaaae7=0
xyzzyaado1=0
do xyzzyaaaa7=1,3
if(.not.xyzzyaabf7(xyzzyaaaa7))cycle
xyzzyaaae7=xyzzyaaae7+1
xyzzyaado1(xyzzyaaaa7)=xyzzyaaae7
enddo
if(levels_striplets>0)then
do xyzzyaaab7=1,nspin
do xyzzyaaac7=xyzzyaaab7,nspin
do xyzzyaaad7=xyzzyaaac7,nspin
no_striplets(1)=no_striplets(1)+1
xyzzyaaan7=0
xyzzyaaap7=0
xyzzyaaao7=0
if(xyzzyaaab7==xyzzyaaac7)xyzzyaaan7=1
if(xyzzyaaac7==xyzzyaaad7)xyzzyaaao7=1
if(xyzzyaaab7==xyzzyaaad7)xyzzyaaap7=1
which_striplet(xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,1)=no_striplets(1)
which_striplet(xyzzyaaac7,xyzzyaaad7,xyzzyaaab7,1)=no_striplets(1)
which_striplet(xyzzyaaad7,xyzzyaaab7,xyzzyaaac7,1)=no_striplets(1)
which_striplet(xyzzyaaad7,xyzzyaaac7,xyzzyaaab7,1)=no_striplets(1)
which_striplet(xyzzyaaac7,xyzzyaaab7,xyzzyaaad7,1)=no_striplets(1)
which_striplet(xyzzyaaab7,xyzzyaaad7,xyzzyaaac7,1)=no_striplets(1)
eq_triplet(xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,1)=xyzzyaaan7+2*xyzzyaaao7&
&+4*xyzzyaaap7
eq_triplet(xyzzyaaac7,xyzzyaaad7,xyzzyaaac7,1)=xyzzyaaao7+2*xyzzyaaap7&
&+4*xyzzyaaan7
eq_triplet(xyzzyaaad7,xyzzyaaab7,xyzzyaaac7,1)=xyzzyaaap7+2*xyzzyaaan7&
&+4*xyzzyaaao7
eq_triplet(xyzzyaaab7,xyzzyaaad7,xyzzyaaac7,1)=xyzzyaaap7+2*xyzzyaaao7&
&+4*xyzzyaaan7
eq_triplet(xyzzyaaad7,xyzzyaaac7,xyzzyaaab7,1)=xyzzyaaao7+2*xyzzyaaan7&
&+4*xyzzyaaap7
eq_triplet(xyzzyaaac7,xyzzyaaab7,xyzzyaaad7,1)=xyzzyaaan7+2*xyzzyaaap7&
&+4*xyzzyaaao7
enddo
enddo
enddo
endif
do xyzzyaaaa7=1,custom_striplets
which_striplet(1:nspin,1:nspin,1:nspin,-xyzzyaaaa7)=xyzzyaaea1(1:nspin&
&,1:nspin,1:nspin,xyzzyaaaa7)
eq_triplet(1:nspin,1:nspin,1:nspin,-xyzzyaaaa7)=xyzzyaaeb1(1:nspin,1:n&
&spin,1:nspin,xyzzyaaaa7)
no_striplets(-xyzzyaaaa7)=maxval(xyzzyaaea1(1:nspin,1:nspin,1:nspin,xy&
&zzyaaaa7))
enddo
if(allocated(xyzzyaaea1))deallocate(xyzzyaaea1,xyzzyaaeb1)
do xyzzyaaab7=1,nspin
do xyzzyaaac7=xyzzyaaab7,nspin
xyzzyaaav7=(xyzzyaaaq7(xyzzyaaab7)==xyzzyaaaq7(xyzzyaaac7))
xyzzyaaax7=(pspin(xyzzyaaab7)==pspin(xyzzyaaac7))
do xyzzyaaaf7=1,xyzzyaaab7
xyzzyaaak7=nspin
if(xyzzyaaaf7==xyzzyaaab7)xyzzyaaak7=xyzzyaaac7-1
do xyzzyaaag7=xyzzyaaaf7,xyzzyaaak7
xyzzyaaaw7=(xyzzyaaaq7(xyzzyaaaf7)==xyzzyaaaq7(xyzzyaaag7))
xyzzyaaay7=(pspin(xyzzyaaaf7)==pspin(xyzzyaaag7))
xyzzyaaaz7=(xyzzyaaav7.eqv.xyzzyaaaw7)
xyzzyaaba7=(xyzzyaaax7.eqv.xyzzyaaay7)
xyzzyaabb7=(xyzzyaaar7(xyzzyaaab7)==xyzzyaaar7(xyzzyaaaf7).and.xyzzyaa&
&ar7(xyzzyaaac7)==xyzzyaaar7(xyzzyaaag7)).or.(xyzzyaaar7(xyzzyaaab7)==&
&xyzzyaaar7(xyzzyaaag7).and.xyzzyaaar7(xyzzyaaac7)==xyzzyaaar7(xyzzyaa&
&af7))
xyzzyaabc7=(xyzzyaaaq7(xyzzyaaab7)==xyzzyaaaq7(xyzzyaaaf7).and.xyzzyaa&
&aq7(xyzzyaaac7)==xyzzyaaaq7(xyzzyaaag7)).or.(xyzzyaaaq7(xyzzyaaab7)==&
&xyzzyaaaq7(xyzzyaaag7).and.xyzzyaaaq7(xyzzyaaac7)==xyzzyaaaq7(xyzzyaa&
&af7))
if(.not.xyzzyaabb7.or..not.xyzzyaaaz7)cycle
xyzzyaaae7=xyzzyaadp1(1)
if(xyzzyaaae7/=0)then
if(which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)==0)then
which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)=which_spair(xyzzyaaaf7,x&
&yzzyaaag7,xyzzyaaae7)
which_spair(xyzzyaaac7,xyzzyaaab7,xyzzyaaae7)=which_spair(xyzzyaaaf7,x&
&yzzyaaag7,xyzzyaaae7)
endif
endif
if(xyzzyaaaz7.and.xyzzyaaav7.and..not.xyzzyaaba7)cycle
xyzzyaaae7=xyzzyaadp1(2)
if(xyzzyaaae7/=0)then
if(which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)==0)then
which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)=which_spair(xyzzyaaaf7,x&
&yzzyaaag7,xyzzyaaae7)
which_spair(xyzzyaaac7,xyzzyaaab7,xyzzyaaae7)=which_spair(xyzzyaaaf7,x&
&yzzyaaag7,xyzzyaaae7)
endif
endif
if(xyzzyaaax7.neqv.xyzzyaaay7)cycle
xyzzyaaae7=xyzzyaadp1(3)
if(xyzzyaaae7/=0)then
if(which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)==0)then
which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)=which_spair(xyzzyaaaf7,x&
&yzzyaaag7,xyzzyaaae7)
which_spair(xyzzyaaac7,xyzzyaaab7,xyzzyaaae7)=which_spair(xyzzyaaaf7,x&
&yzzyaaag7,xyzzyaaae7)
endif
endif
if(.not.xyzzyaabc7)cycle
xyzzyaaae7=xyzzyaadp1(4)
if(xyzzyaaae7/=0)then
if(which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)==0)then
which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)=which_spair(xyzzyaaaf7,x&
&yzzyaaag7,xyzzyaaae7)
which_spair(xyzzyaaac7,xyzzyaaab7,xyzzyaaae7)=which_spair(xyzzyaaaf7,x&
&yzzyaaag7,xyzzyaaae7)
endif
endif
enddo
if(all(which_spair(xyzzyaaab7,xyzzyaaac7,:)/=0))exit
enddo
do xyzzyaaaa7=1,5
xyzzyaaae7=xyzzyaadp1(xyzzyaaaa7)
if(xyzzyaaae7==0)cycle
if(which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)==0)then
no_spairs(xyzzyaaae7)=no_spairs(xyzzyaaae7)+1
which_spair(xyzzyaaab7,xyzzyaaac7,xyzzyaaae7)=no_spairs(xyzzyaaae7)
which_spair(xyzzyaaac7,xyzzyaaab7,xyzzyaaae7)=no_spairs(xyzzyaaae7)
endif
enddo
enddo
enddo
do xyzzyaaaa7=1,custom_spairs
which_spair(1:nspin,1:nspin,-xyzzyaaaa7)=xyzzyaadz1(1:nspin,1:nspin,xy&
&zzyaaaa7)
no_spairs(-xyzzyaaaa7)=maxval(xyzzyaadz1(1:nspin,1:nspin,xyzzyaaaa7))
enddo
if(allocated(xyzzyaadz1))deallocate(xyzzyaadz1)
do xyzzyaaab7=1,nspin
do xyzzyaaaf7=1,xyzzyaaab7-1
if(xyzzyaaar7(xyzzyaaab7)/=xyzzyaaar7(xyzzyaaaf7))cycle
xyzzyaaae7=xyzzyaado1(1)
if(xyzzyaaae7/=0)then
if(which_ssingle(xyzzyaaab7,xyzzyaaae7)==0)which_ssingle(xyzzyaaab7,xy&
&zzyaaae7)=which_ssingle(xyzzyaaaf7,xyzzyaaae7)
endif
if(xyzzyaaaq7(xyzzyaaab7)/=xyzzyaaaq7(xyzzyaaaf7))cycle
xyzzyaaae7=xyzzyaado1(2)
if(xyzzyaaae7/=0)then
if(which_ssingle(xyzzyaaab7,xyzzyaaae7)==0)which_ssingle(xyzzyaaab7,xy&
&zzyaaae7)=which_ssingle(xyzzyaaaf7,xyzzyaaae7)
endif
enddo
do xyzzyaaaa7=1,3
xyzzyaaae7=xyzzyaado1(xyzzyaaaa7)
if(xyzzyaaae7==0)cycle
if(which_ssingle(xyzzyaaab7,xyzzyaaae7)==0)then
no_ssingles(xyzzyaaae7)=no_ssingles(xyzzyaaae7)+1
which_ssingle(xyzzyaaab7,xyzzyaaae7)=no_ssingles(xyzzyaaae7)
endif
enddo
enddo
do xyzzyaaaa7=1,custom_ssingles
which_ssingle(1:nspin,-xyzzyaaaa7)=xyzzyaady1(1:nspin,xyzzyaaaa7)
no_ssingles(-xyzzyaaaa7)=maxval(xyzzyaady1(1:nspin,xyzzyaaaa7))
enddo
if(allocated(xyzzyaady1))deallocate(xyzzyaady1)
if(xyzzyaaai7==1)then
level_family=0
elseif(xyzzyaabf7(2))then
level_family=xyzzyaado1(2)
else
level_family=xyzzyaado1(3)
endif
no_families=no_ssingles(level_family)
allocate(fam_charge(no_families),which_fam(nspin),stat=xyzzyaaah7)
call check_alloc(xyzzyaaah7,'ASSIGN_SPIN_DEPS','fam_charge/which_fam')
do xyzzyaaaa7=1,nspin
which_fam(xyzzyaaaa7)=which_ssingle(xyzzyaaaa7,level_family)
fam_charge(which_fam(xyzzyaaaa7))=pcharge(xyzzyaaaa7)
enddo
if(xyzzyaaaj7==1)then
level_eqvfamily=0
elseif(xyzzyaabf7(1))then
level_eqvfamily=xyzzyaado1(1)
else
level_eqvfamily=level_family
endif
no_eqvfamilies=no_ssingles(level_eqvfamily)
allocate(which_eqvfam(nspin),stat=xyzzyaaah7)
call check_alloc(xyzzyaaah7,'ASSIGN_SPIN_DEPS','which_eqvfam')
do xyzzyaaaa7=1,nspin
which_eqvfam(xyzzyaaaa7)=which_ssingle(xyzzyaaaa7,level_eqvfamily)
enddo
allocate(difftype_mass(no_difftypes),xyzzyaaei1(no_difftypes),xyzzyaae&
&j1(no_difftypes),xyzzyaaec1(no_difftypes),stat=xyzzyaaah7)
call check_alloc(xyzzyaaah7,'ASSIGN_SPIN_DEPS','difftype_mass etc.')
do xyzzyaaaa7=1,nspin
difftype_mass(which_difftype(xyzzyaaaa7))=pmass(xyzzyaaaa7)
enddo
allocate(which_spin(netot),which_ie(netot),which_ii(nemax,nspin),which&
&_ee(netot,netot),stat=xyzzyaaah7)
call check_alloc(xyzzyaaah7,'ASSIGN_SPIN_DEPS','which_spin etc.')
which_ii=0
xyzzyaaae7=0
do xyzzyaaal7=1,nspin
do xyzzyaaam7=1,nele(xyzzyaaal7)
xyzzyaaae7=xyzzyaaae7+1
which_spin(xyzzyaaae7)=xyzzyaaal7
which_ie(xyzzyaaae7)=xyzzyaaam7
which_ii(xyzzyaaam7,xyzzyaaal7)=xyzzyaaae7
enddo
enddo
xyzzyaaac7=0
do xyzzyaaaa7=1,netot
which_ee(xyzzyaaaa7,xyzzyaaaa7)=0
do xyzzyaaab7=xyzzyaaaa7+1,netot
xyzzyaaac7=xyzzyaaac7+1
which_ee(xyzzyaaaa7,xyzzyaaab7)=xyzzyaaac7
which_ee(xyzzyaaab7,xyzzyaaaa7)=xyzzyaaac7
enddo
enddo
allocate(ee_cusp_in_orbital(nspin,nspin),stat=xyzzyaaah7)
call check_alloc(xyzzyaaah7,'ASSIGN_SPIN_DEPS','ee_cusp_in_orbital')
ee_cusp_in_orbital(:,:)=.false.
deallocate(xyzzyaaaq7,xyzzyaaar7,xyzzyaabe7,xyzzyaabf7)
end subroutine xyzzyaafa1
subroutine xyzzyaafb1
implicit none
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaa&
&af8,xyzzyaaag8,xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaak8,xyzzyaaal8&
&,xyzzyaaam8,xyzzyaaan8,xyzzyaaao8,xyzzyaaap8,xyzzyaaaq8,xyzzyaaar8,xy&
&zzyaaas8,xyzzyaaat8
logical xyzzyaaau8,xyzzyaaav8
character(80) char80
if(esdf_block('custom_striplet_dep',xyzzyaaaa8))then
xyzzyaaab8=1
read(block_data(1),*,iostat=xyzzyaaac8)char80,custom_striplets
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','First line&
& of custom_striplet_dep should be of the form "no_striplet_deps <inte&
&ger>".')
if(trim(char80)/='no_striplet_deps')call errstop_master('READ_CUSTOM_S&
&PINDEP',' Was expecting to find "no_striplet_deps" in first line of c&
&ustom_striplet_dep block.')
if(custom_striplets<0)call errstop_master('READ_CUSTOM_SPINDEP','Numbe&
&r of spin-triplet sets should be positive in first line of custom_str&
&iplet_dep.')
allocate(xyzzyaaea1(nspin,nspin,nspin,custom_striplets),xyzzyaaeb1(nsp&
&in,nspin,nspin,custom_striplets),stat=xyzzyaaat8)
call check_alloc(xyzzyaaat8,'READ_CUSTOM_SPINDEP','')
xyzzyaaea1=0
xyzzyaaeb1=0
do xyzzyaaad8=1,custom_striplets
xyzzyaaab8=xyzzyaaab8+1
read(block_data(xyzzyaaab8),*,iostat=xyzzyaaac8)char80,xyzzyaaaf8,xyzz&
&yaaag8
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_striplet_dep at line '//trim(i2s(xyzzyaaab8))//'.')
if(trim(char80)/='striplet_dep')call errstop_master('READ_CUSTOM_SPIND&
&EP','Was expecting to find "striplet_dep" in custom_striplet_dep line&
& '//trim(i2s(xyzzyaaab8))//'.')
if(xyzzyaaaf8/=-xyzzyaaad8)call errstop_master('READ_CUSTOM_SPINDEP','&
&Have found a spin-triplet label of '//trim(i2s(xyzzyaaaf8))//' instea&
&d of '//trim(i2s(-xyzzyaaad8))//' at line '//trim(i2s(xyzzyaaab8))//'&
& of custom_striplet_dep block.')
do xyzzyaaae8=1,xyzzyaaag8
xyzzyaaab8=xyzzyaaab8+1
read(block_data(xyzzyaaab8),'(a)',iostat=xyzzyaaac8)char80
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_striplet_dep at line '//trim(i2s(xyzzyaaab8))//'.')
xyzzyaaam8=index(char80,'#')
if(xyzzyaaam8>0)char80=char80(1:xyzzyaaam8-1)
xyzzyaaas8=-1
do
xyzzyaaao8=index(char80,'-')
xyzzyaaan8=index(char80,'=')
xyzzyaaau8=(xyzzyaaao8<1.or.xyzzyaaao8>xyzzyaaan8).and.xyzzyaaan8>0
if(xyzzyaaau8)xyzzyaaao8=xyzzyaaan8
if(xyzzyaaao8<=1)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_striplet_dep at line '//trim(i2s(xyzzyaaab8))//' <1>.')
read(char80(1:xyzzyaaao8-1),*,iostat=xyzzyaaac8)xyzzyaaah8
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_striplet_dep at line '//trim(i2s(xyzzyaaab8))//' <2>.')
char80=char80(xyzzyaaao8+1:)
xyzzyaaao8=index(char80,'-')
xyzzyaaan8=index(char80,'=')
xyzzyaaav8=(xyzzyaaao8<1.or.xyzzyaaao8>xyzzyaaan8).and.xyzzyaaan8>0
if(xyzzyaaav8)xyzzyaaao8=xyzzyaaan8
if(xyzzyaaao8<=1)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_striplet_dep at line '//trim(i2s(xyzzyaaab8))//' <3>.')
read(char80(1:xyzzyaaao8-1),*,iostat=xyzzyaaac8)xyzzyaaai8
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_striplet_dep at line '//trim(i2s(xyzzyaaab8))//' <4>.')
char80=char80(xyzzyaaao8+1:)
xyzzyaaal8=index(char80,',')
if(xyzzyaaal8==1)then
call errstop_master('READ_CUSTOM_SPINDEP','Problem reading custom_stri&
&plet_dep at line '//trim(i2s(xyzzyaaab8))//' <5>.')
elseif(xyzzyaaal8==0)then
read(char80,*,iostat=xyzzyaaac8)xyzzyaaaj8
else
read(char80(1:xyzzyaaal8-1),*,iostat=xyzzyaaac8)xyzzyaaaj8
endif
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_striplet_dep at line '//trim(i2s(xyzzyaaab8))//' <6>.')
if(xyzzyaaah8<1.or.xyzzyaaai8<1.or.xyzzyaaaj8<1.or.xyzzyaaah8>nspin.or&
&.xyzzyaaai8>nspin.or.xyzzyaaaj8>nspin)call errstop_master('READ_CUSTO&
&M_SPINDEP','Found a spin that was not in the range 1 -- '//trim(i2s(n&
&spin))//' in custom_striplet_block, line '//trim(i2s(xyzzyaaab8))//'.&
&')
if(xyzzyaaea1(xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaad8)/=0)call err&
&stop_master('READ_CUSTOM_SPINDEP','The spin dependence for '//trim(i2&
&s(xyzzyaaah8))//'-'//trim(i2s(xyzzyaaai8))//'-'//trim(i2s(xyzzyaaaj8)&
&)//' seems to be defined twice in the custom_striplet_dep block.')
xyzzyaaea1(xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaad8)=xyzzyaaae8
xyzzyaaea1(xyzzyaaai8,xyzzyaaaj8,xyzzyaaah8,xyzzyaaad8)=xyzzyaaae8
xyzzyaaea1(xyzzyaaaj8,xyzzyaaah8,xyzzyaaai8,xyzzyaaad8)=xyzzyaaae8
xyzzyaaea1(xyzzyaaah8,xyzzyaaaj8,xyzzyaaai8,xyzzyaaad8)=xyzzyaaae8
xyzzyaaea1(xyzzyaaaj8,xyzzyaaai8,xyzzyaaah8,xyzzyaaad8)=xyzzyaaae8
xyzzyaaea1(xyzzyaaai8,xyzzyaaah8,xyzzyaaaj8,xyzzyaaad8)=xyzzyaaae8
xyzzyaaap8=0
if(xyzzyaaau8.or.xyzzyaaah8==xyzzyaaai8)xyzzyaaap8=1
xyzzyaaaq8=0
if(xyzzyaaav8.or.xyzzyaaai8==xyzzyaaaj8)xyzzyaaaq8=1
xyzzyaaar8=0
if((xyzzyaaap8==1.and.xyzzyaaaq8==1).or.xyzzyaaah8==xyzzyaaaj8)xyzzyaa&
&ar8=1
xyzzyaaaf8=xyzzyaaap8+xyzzyaaaq8+xyzzyaaar8
if(xyzzyaaaf8==2)call errstop_master('READ_CUSTOM_SPINDEP','Bug in the&
& spin-triplet equivalence assignment algorithm.')
if(xyzzyaaas8<0)xyzzyaaas8=xyzzyaaaf8
if(xyzzyaaaf8/=xyzzyaaas8)call errstop_master('READ_CUSTOM_SPINDEP','A&
&ll triplets in a custom spin-triplet group should have the same numbe&
&r of equivalences.')
xyzzyaaeb1(xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaad8)=xyzzyaaap8+2*x&
&yzzyaaaq8+4*xyzzyaaar8
xyzzyaaeb1(xyzzyaaai8,xyzzyaaaj8,xyzzyaaah8,xyzzyaaad8)=xyzzyaaaq8+2*x&
&yzzyaaar8+4*xyzzyaaap8
xyzzyaaeb1(xyzzyaaaj8,xyzzyaaah8,xyzzyaaai8,xyzzyaaad8)=xyzzyaaar8+2*x&
&yzzyaaap8+4*xyzzyaaaq8
xyzzyaaeb1(xyzzyaaah8,xyzzyaaaj8,xyzzyaaai8,xyzzyaaad8)=xyzzyaaar8+2*x&
&yzzyaaaq8+4*xyzzyaaap8
xyzzyaaeb1(xyzzyaaaj8,xyzzyaaai8,xyzzyaaah8,xyzzyaaad8)=xyzzyaaaq8+2*x&
&yzzyaaap8+4*xyzzyaaar8
xyzzyaaeb1(xyzzyaaai8,xyzzyaaah8,xyzzyaaaj8,xyzzyaaad8)=xyzzyaaap8+2*x&
&yzzyaaar8+4*xyzzyaaaq8
if(xyzzyaaal8==0)exit
char80=char80(xyzzyaaal8+1:)
enddo
enddo
if(any(xyzzyaaea1(:,:,:,xyzzyaaad8)<=0))call errstop_master('READ_CUST&
&OM_SPINDEP','Some spin triplets are not placed in groups in the custo&
&m_striplet_dep block.')
enddo
if(xyzzyaaab8/=xyzzyaaaa8)call errstop_master('READ_CUSTOM_SPINDEP','P&
&roblem reading custom_striplet_dep block - too many lines.')
else
custom_striplets=0
endif
if(esdf_block('custom_spair_dep',xyzzyaaaa8))then
xyzzyaaab8=1
read(block_data(1),*,iostat=xyzzyaaac8)char80,custom_spairs
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','First line&
& of custom_spair_dep should be of the form "no_spair_deps <integer>".&
&')
if(trim(char80)/='no_spair_deps')call errstop_master('READ_CUSTOM_SPIN&
&DEP',' Was expecting to find "no_spair_deps" in first line of custom_&
&spair_dep block.')
if(custom_spairs<0)call errstop_master('READ_CUSTOM_SPINDEP','Number o&
&f spin-pair sets should be positive in first line of custom_spair_dep&
&.')
allocate(xyzzyaadz1(nspin,nspin,custom_spairs),stat=xyzzyaaat8)
call check_alloc(xyzzyaaat8,'READ_CUSTOM_SPINDEP','')
xyzzyaadz1=0
do xyzzyaaad8=1,custom_spairs
xyzzyaaab8=xyzzyaaab8+1
read(block_data(xyzzyaaab8),*,iostat=xyzzyaaac8)char80,xyzzyaaaf8,xyzz&
&yaaag8
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_spair_dep at line '//trim(i2s(xyzzyaaab8))//'.')
if(trim(char80)/='spair_dep')call errstop_master('READ_CUSTOM_SPINDEP'&
&,'Was expecting to find "spair_dep" in custom_spair_dep line '//trim(&
&i2s(xyzzyaaab8))//'.')
if(xyzzyaaaf8/=-xyzzyaaad8)call errstop_master('READ_CUSTOM_SPINDEP','&
&Have found a spin-pair label of '//trim(i2s(xyzzyaaaf8))//' instead o&
&f '//trim(i2s(-xyzzyaaad8))//' at line '//trim(i2s(xyzzyaaab8))//' of&
& custom_spair_dep block.')
do xyzzyaaae8=1,xyzzyaaag8
xyzzyaaab8=xyzzyaaab8+1
read(block_data(xyzzyaaab8),'(a)',iostat=xyzzyaaac8)char80
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_spair_dep at line '//trim(i2s(xyzzyaaab8))//'.')
xyzzyaaam8=index(char80,'#')
if(xyzzyaaam8>0)char80=char80(1:xyzzyaaam8-1)
do
xyzzyaaak8=index(char80,'-')
if(xyzzyaaak8<=1)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_spair_dep at line '//trim(i2s(xyzzyaaab8))//'.')
read(char80(1:xyzzyaaak8-1),*,iostat=xyzzyaaac8)xyzzyaaah8
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_spair_dep at line '//trim(i2s(xyzzyaaab8))//'.')
char80=char80(xyzzyaaak8+1:)
xyzzyaaal8=index(char80,',')
if(xyzzyaaal8==1)then
call errstop_master('READ_CUSTOM_SPINDEP','Problem reading custom_spai&
&r_dep at line '//trim(i2s(xyzzyaaab8))//'.')
elseif(xyzzyaaal8==0)then
read(char80,*,iostat=xyzzyaaac8)xyzzyaaai8
else
read(char80(1:xyzzyaaal8-1),*,iostat=xyzzyaaac8)xyzzyaaai8
endif
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_spair_dep at line '//trim(i2s(xyzzyaaab8))//'.')
if(xyzzyaaah8<1.or.xyzzyaaai8<1.or.xyzzyaaah8>nspin.or.xyzzyaaai8>nspi&
&n)call errstop_master('READ_CUSTOM_SPINDEP','Found a spin that was no&
&t in the range 1 -- '//trim(i2s(nspin))//' in custom_spair_block, lin&
&e '//trim(i2s(xyzzyaaab8))//'.')
if(xyzzyaadz1(xyzzyaaah8,xyzzyaaai8,xyzzyaaad8)/=0)call errstop_master&
&('READ_CUSTOM_SPINDEP','The spin dependence for '//trim(i2s(xyzzyaaah&
&8))//'-'//trim(i2s(xyzzyaaai8))//' seems to be defined twice in the c&
&ustom_spair_dep block.')
xyzzyaadz1(xyzzyaaah8,xyzzyaaai8,xyzzyaaad8)=xyzzyaaae8
xyzzyaadz1(xyzzyaaai8,xyzzyaaah8,xyzzyaaad8)=xyzzyaaae8
if(xyzzyaaal8==0)exit
char80=char80(xyzzyaaal8+1:)
enddo
enddo
if(any(xyzzyaadz1(:,:,xyzzyaaad8)<=0))call errstop_master('READ_CUSTOM&
&_SPINDEP','Some spin pairs are not placed in groups in the custom_spa&
&ir_dep block.')
enddo
if(xyzzyaaab8/=xyzzyaaaa8)call errstop_master('READ_CUSTOM_SPINDEP','P&
&roblem reading custom_spair_dep block - too many lines.')
else
custom_spairs=0
endif
if(esdf_block('custom_ssingle_dep',xyzzyaaaa8))then
xyzzyaaab8=1
read(block_data(1),*,iostat=xyzzyaaac8)char80,custom_ssingles
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','First line&
& of custom_ssingle_dep should be of the form "no_ssingle_deps <intege&
&r>".')
if(trim(char80)/='no_ssingle_deps')call errstop_master('READ_CUSTOM_SP&
&INDEP',' Was expecting to find "no_ssingle_deps" in first line of cus&
&tom_ssingle_dep block.')
if(custom_ssingles<0)call errstop_master('READ_CUSTOM_SPINDEP','Number&
& of spin singles should be positive.')
allocate(xyzzyaady1(nspin,custom_ssingles),stat=xyzzyaaat8)
call check_alloc(xyzzyaaat8,'READ_CUSTOM_SPINDEP','')
xyzzyaady1=0
do xyzzyaaad8=1,custom_ssingles
xyzzyaaab8=xyzzyaaab8+1
read(block_data(xyzzyaaab8),*,iostat=xyzzyaaac8)char80,xyzzyaaaf8,xyzz&
&yaaag8
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_ssingle_dep at line '//trim(i2s(xyzzyaaab8))//'.')
if(trim(char80)/='ssingle_dep')call errstop_master('READ_CUSTOM_SPINDE&
&P','Was expecting to find "ssingle_dep" at line '//trim(i2s(xyzzyaaab&
&8))//' of custom_ssingle_block.')
if(xyzzyaaaf8/=-xyzzyaaad8)call errstop_master('READ_CUSTOM_SPINDEP','&
&Have found a spin-single label of '//trim(i2s(xyzzyaaaf8))//' instead&
& of ' //trim(i2s(-xyzzyaaad8))//' at line '//trim(i2s(xyzzyaaab8))//'&
& of custom_ssingle_dep block.')
do xyzzyaaae8=1,xyzzyaaag8
xyzzyaaab8=xyzzyaaab8+1
read(block_data(xyzzyaaab8),'(a)',iostat=xyzzyaaac8)char80
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_ssingle_dep at line '//trim(i2s(xyzzyaaab8))//'.')
xyzzyaaam8=index(char80,'#')
if(xyzzyaaam8>0)char80=char80(1:xyzzyaaam8-1)
do
xyzzyaaal8=index(char80,',')
if(xyzzyaaal8==1)then
call errstop_master('READ_CUSTOM_SPINDEP','Problem reading custom_ssin&
&gle_dep at line '//trim(i2s(xyzzyaaab8))//'.')
elseif(xyzzyaaal8==0)then
read(char80,*,iostat=xyzzyaaac8)xyzzyaaah8
else
read(char80(1:xyzzyaaal8-1),*,iostat=xyzzyaaac8)xyzzyaaah8
endif
if(xyzzyaaac8/=0)call errstop_master('READ_CUSTOM_SPINDEP','Problem re&
&ading custom_ssingle_dep at line '//trim(i2s(xyzzyaaab8))//'.')
if(xyzzyaaah8<1.or.xyzzyaaah8>nspin)call errstop_master('READ_CUSTOM_S&
&PINDEP','Found a spin that was not in the range 1 -- '//trim(i2s(nspi&
&n))//' in custom_ssingle_block, line '//trim(i2s(xyzzyaaab8))//'.')
if(xyzzyaady1(xyzzyaaah8,xyzzyaaad8)/=0)call errstop_master('READ_CUST&
&OM_SPINDEP','The spin dependence for '//trim(i2s(xyzzyaaah8))//' seem&
&s to be defined twice in the custom_ssingle_block.')
xyzzyaady1(xyzzyaaah8,xyzzyaaad8)=xyzzyaaae8
if(xyzzyaaal8==0)exit
char80=char80(xyzzyaaal8+1:)
enddo
enddo
if(any(xyzzyaady1(:,xyzzyaaad8)<=0))call errstop_master('READ_CUSTOM_S&
&PINDEP','Some spins are not placed in groups in the custom_ssingle_de&
&p block.')
enddo
if(xyzzyaaab8/=xyzzyaaaa8)call errstop_master('READ_CUSTOM_SPINDEP','P&
&roblem reading custom_ssingle_dep block - too many lines.')
else
custom_ssingles=0
endif
end subroutine xyzzyaafb1
subroutine xyzzyaafc1
implicit none
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9
logical xyzzyaaad9
if(am_master)then
xyzzyaaad9=.false.
if(xyzzyaaad1<0)then
call wout('CORPER in input must be positive, or zero for automatic set&
&ting.')
call wout()
xyzzyaaad9=.true.
endif
if(dbarrc<1)then
call wout('DBARRC in input must be positive.')
call wout()
xyzzyaaad9=.true.
else
if(dbarrc<100000)then
call wordwrap('Invalid value of DBARRC in input. Since Feb 2014, CASIN&
&O will block any attempt to lower DBARRC below the default level of 1&
&00000, as this can lead to very significant increases in CPU time. If&
& you have a legitimate reason for wishing to do this, you must commen&
&t out this error trap in the source and in the runqmc script.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(xyzzyaaas1<1)then
call wout('NON_LOCAL_GRID in input must be positive.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaacd1<=0.d0)then
call wout('DTVMC in input must be positive.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaace1<=0.d0)then
call wout('DTDMC in input must be positive.')
call wout()
xyzzyaaad9=.true.
endif
if(movieplot<1)then
call wout('MOVIEPLOT in input must be positive.')
call wout()
xyzzyaaad9=.true.
endif
if(orb_norm<=0.d0)then
call wout('ORB_NORM in input must be a positive real number.')
call wout()
xyzzyaaad9=.true.
endif
if(movienode<0.or.movienode>nnodes-1)then
call wordwrap('MOVIENODE in input must be zero or positive integer les&
&s than the total number of nodes.')
call wout()
xyzzyaaad9=.true.
endif
if(trim(atom_basis_type)=='none'.and..not.allocated(heg_orbtype))then
call wout('ATOM_BASIS_TYPE=''none'' in input and no free-particle orbi&
&tals defined.')
xyzzyaaad9=.true.
endif
if(((trim(atom_basis_type)=='numerical'.or.trim(atom_basis_type)=='dim&
&er').and.isperiodic).or.(trim(atom_basis_type)=='plane-wave'.and..not&
&.isperiodic))then
call wout('Illegal combination of basis set and isperiodic flag.')
call wout('ATOM_BASIS_TYPE : '//trim(atom_basis_type))
call wout('PERIODIC        : '//l2s(isperiodic))
call wout()
xyzzyaaad9=.true.
endif
if(limdmc<0.or.limdmc>4)then
call wout('LIMDMC has an illegal value. Choose one of 0,1,2,3,4 (prefe&
&rably 4).')
call wout('Current value : '//trim(i2s(limdmc)))
call wout()
xyzzyaaad9=.true.
endif
if(alimit<=0.d0)then
call wout('The ALIMIT parameter in input should be strictly positive.'&
&)
call wout()
xyzzyaaad9=.true.
endif
if((xyzzyaact1.and.nucleus_gf_mods).and.limdmc/=2.and.limdmc/=4)then
call errwarn('CHECK_INPUT_PARAMETERS','LIMDMC in input should equal 2 &
&or 4 if the nucleus GF modifications are to be applied.')
endif
if(sparse_threshold<=0.d0)then
call wordwrap('The SPARSE_THRESHOLD input parameter must be greater th&
&an zero.  Change the code to call the MXMB routine rather than MXMB_W&
&ITH_THRESHOLD if you really want to do this.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaab1<=0)then
call wordwrap('The NMOVE input parameter is not allowed to be zero, si&
&nce various things need to be divided by it.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaak1<=0.and..not.iaccum)then
call wordwrap('The NMOVE_DMC_EQUIL input parameter is not allowed to b&
&e zero, since various things need to be divided by it.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaal1==0.and.iaccum)then
call wordwrap('The NMOVE_DMC_STATS in input parameter is not allowed t&
&o be zero, since various things need to be divided by it.')
call wout()
xyzzyaaad9=.true.
endif
if(opt_info<1.or.opt_info>5)then
call wout('OPT_INFO input parameter must be 1, 2, 3, 4 or 5.')
call wout()
xyzzyaaad9=.true.
endif
if(vm_w_max<1.d0.and.vm_w_max/=0.d0)then
call wout('VM_W_MAX input parameter must be greater than one or equal &
&to zero.')
call wout()
xyzzyaaad9=.true.
endif
if(vm_w_min>1.d0)then
call wout('VM_W_MIN input parameter must be less than one.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaaj1<1)then
call wout('Best to average over at least one successive move in VMC.')
call wout('Choose another value for NVMCAVE in input.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaao1/=1.and.xyzzyaaao1/=3)then
call wout('VMC_METHOD input parameter can only take values 1 or 3. The&
& default is 1.')
call wout()
xyzzyaaad9=.true.
endif
select case(trim(vmc_sampling))
case('standard','optimum','HF optimum','efficient','new','stdtest')
continue
case default
call wordwrap('VMC_SAMPLING can only take the values "standard", "opti&
&mum","HF optimum", "efficient", "new" and "stdtest".')
call wout()
xyzzyaaad9=.true.
end select
if(vmc_ionjump<0.d0.or.vmc_ionjump>1.d0)then
call wout('VMC_IONJUMP input parameter must be a probability (between &
&0.0 and 1.0).')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaap1<1.or.xyzzyaaap1>2)then
call wout('DMC_METHOD in input can only take values 1 or 2. The defaul&
&t is 1.')
call wout()
xyzzyaaad9=.true.
endif
if((isopt.or.isopt_vmc.or.isvmc_opt).and.any(xyzzyaadb1).and..not.use_&
&backflow)then
call wordwrap('The BACKFLOW keyword in input needs to be set to T when&
& optimizing backflow parameters (OPT_BACKFLOW=T).')
call wout()
xyzzyaaad9=.true.
endif
if(forces.and.trim(atom_basis_type)/='gaussian')then
call wordwrap('The calculation of atomic forces is only implemented fo&
&r the case of Gaussian basis sets.')
call wout()
xyzzyaaad9=.true.
endif
if(forces.and.isvmc.and..not.writeout_vmc_hist)then
call wout('When calculating forces in VMC, use writeout_vmc_hist=T in &
&input.')
call wout()
xyzzyaaad9=.true.
endif
if(forces.and..not.(isrmc.or.isrmc_rmc.or.isvmc.or.isdmc.or.isvmc_dmc.&
&or.isdmc_dmc.or.isvmc_dmc_equil))then
call wout('Set forces flag to F in input for this type of calculation.&
&')
call wout()
xyzzyaaad9=.true.
endif
if(forces.and.(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvmc_dmc_equil).and&
&..not.use_future)then
call wout('When calculating forces in DMC, choose use_future=T in inpu&
&t.')
call wout()
xyzzyaaad9=.true.
endif
if(noncoll_spin.and.(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvmc_dmc_equi&
&l))then
call wout('Non-collinear spin calculation is not possible with DMC.')
call wout()
xyzzyaaad9=.true.
endif
if(allocated(heg_orbtype).and.dimensionality/=3)then
if(any(heg_orbtype==3))then
call wout('Spin density wave orbitals can only be used in 3 dimensions&
&.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(ranluxlevel<0.or.ranluxlevel>4)then
call wout('RANLUXLEVEL input parameter out of range. Should be 0-4 (de&
&fault 3).')
call wout()
xyzzyaaad9=.true.
endif
if(ranprint<0)then
call wout('RANPRINT input parameter may not be negative.')
call wout()
xyzzyaaad9=.true.
endif
if(tpdmc<0)then
call wout('TPDMC input parameter may not be negative.')
call wout()
xyzzyaaad9=.true.
endif
if(neighprint<0)then
call wout('NEIGHPRINT input parameter may not be negative.')
call wout()
xyzzyaaad9=.true.
endif
if(max_cpu_time<0.)then
call wout('MAX_CPU_TIME input parameter may not be negative.')
call wout()
xyzzyaaad9=.true.
endif
if(max_real_time<0.)then
call wout('MAX_REAL_TIME input parameter may not be negative.')
call wout()
xyzzyaaad9=.true.
endif
if(block_time<0.)then
call wout('BLOCK_TIME input parameter may not be negative.')
call wout()
xyzzyaaad9=.true.
endif
if(use_blocktime.and.old_input)then
call wordwrap('The BLOCK_TIME algorithm is not implemented for the old&
& set of input keywords (NMOVE, etc.). Please convert your input files&
& to the new set - the input_kw_conv utility might help..')
call wout()
xyzzyaaad9=.true.
endif
if(isvmc.or.isvmc_dmc.or.isvmc_opt.or.isopt_vmc)then
if(xyzzyaaae1>xyzzyaaac1*xyzzyaaab1*xyzzyaaaj1)then
call wordwrap('Number of configs to be written is greater than the tot&
&al number of moves.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvmc_dmc_equil)then
if(targ_wt<1.d0)then
call wout('Target weight for DMC calculation is less than 1.')
xyzzyaaad9=.true.
endif
if(trip_popn<=targ_wt.and.trip_popn/=0.d0)then
call wordwrap('DMC_TRIP_WEIGHT is less than or equal to the target wei&
&ght DMC_TARGET_WEIGHT and hence the DMC run will stop immediately. Su&
&ggested value for DMC_TRIP_WEIGHT is three times DMC_TARGET_WEIGHT.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(xyzzyaaaq1<0.or.xyzzyaaaq1>2)then
call wout('OPT_DTVMC input parameter must be 0, 1 or 2.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaaq1/=0)then
if(xyzzyaaao1==1)then
xyzzyaaac9=max(100,500/netot)
if(xyzzyaaaa1<xyzzyaaac9)then
call wordwrap('NEQUIL < '//trim(i2s(xyzzyaaac9))//' and OPT_DTVMC > 0.&
& Must use at least '//trim(i2s(xyzzyaaac9))//' equilibration steps fo&
&r this system in order to optimize the VMC time step using VMC_METHOD&
&=1.')
call wout()
xyzzyaaad9=.true.
endif
elseif(xyzzyaaao1==3.and.xyzzyaaaa1<2000)then
call wordwrap('NEQUIL < 2000 and OPT_DTVMC > 0. Must use at least 2000&
& equilibration steps in order to optimize the VMC time step using VMC&
&_METHOD=3.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(isperiodic.and.use_orbmods)then
call wordwrap('Modifications of single particle orbitals not yet imple&
&mented for periodic systems.')
call wout()
xyzzyaaad9=.true.
endif
if(isopt.or.isvmc_opt.or.isopt_vmc)then
if(.not.all(xyzzyaadc1.or.xyzzyaada1.or.xyzzyaadd1.or.xyzzyaadb1.or.xy&
&zzyaade1))then
call wordwrap('During wave function optimization one of OPTJASTROW, OP&
&TDETCOEFF, OPTORBITALS or OPTBACKFLOW should be set to T in the input&
& file, otherwise you won''t have anything to optimize.  If using the &
&OPT_PLAN input block, all optimization cycles must have something to &
&optimize.')
call wout()
xyzzyaaad9=.true.
endif
if(vm_use_e_guess.and.vm_e_guess==-999.d0)then
call wout('If VM_USE_E_GUESS is set then VM_E_GUESS must be supplied.'&
&)
call wout()
xyzzyaaad9=.true.
endif
if(vm_filter.and.vm_reweight)then
call wout('If VM_FILTER is T then VM_REWEIGHT must be F.')
call wout()
xyzzyaaad9=.true.
endif
if(vm_filter.and.(vm_filter_thres<0.d0.or.vm_filter_width<0.d0))then
call wout('VM_FILTER_THRES and VM_FILTER_WIDTH must be >=0 in input.')
call wout()
xyzzyaaad9=.true.
endif
if(any(xyzzyaada1).and..not.use_jastrow)then
call wout('OPT_JASTROW is T but USE_JASTROW is F in input.')
call wout()
xyzzyaaad9=.true.
endif
endif
select case(trim(vm_linjas_method))
case('CG','MC','SD','LM','CG_MC','BFGS','BFGS_MC','GN','GN_MC')
case default
call wordwrap('VM_LINJAS_METHOD should be one of: CG, MC, SD, LM, CG_M&
&C, BFGS BFGS_MC, GN or GN_MC.')
call wout()
xyzzyaaad9=.true.
end select
if(vm_linjas_its==0)then
call wout('VM_LINJAS_ITS input parameter should not be zero.')
call wout()
xyzzyaaad9=.true.
endif
if(molgscreening.and.use_backflow)then
call wout('Backflow may not be used with MOLGSCREENING set to T in inp&
&ut.')
call wout()
xyzzyaaad9=.true.
endif
if(isopt.or.isvmc_opt.or.isopt_vmc)then
if(any(index(xyzzyaacz1,'varmin_linjas')==1).and..not.use_jastrow)then
call wout('OPT_METHOD=varmin_linjas in input can only be used if a Jas&
&trow factor is present.')
call wout()
xyzzyaaad9=.true.
endif
if(any(index(xyzzyaacz1,'varmin_linjas')==1.and.(xyzzyaadc1.or.xyzzyaa&
&dd1.or.xyzzyaadb1.or.xyzzyaade1)))then
call wout('OPT_METHOD=varmin_linjas can only be used to optimize Jastr&
&ow factors.')
call wout()
xyzzyaaad9=.true.
endif
if(any(index(xyzzyaacz1,'varmin_linjas')==1.and..not.xyzzyaada1))then
call wout('Must set OPT_JASTROW=T in input when OPT_METHOD=varmin_linj&
&as.')
call wout()
xyzzyaaad9=.true.
endif
if(any(index(xyzzyaacz1,'varmin_linjas')==1).and.complex_wf)then
call wordwrap('Varmin_linjas can only be used for real wave functions &
&at present.  Please either (i) choose a different OPT_METHOD, (ii) op&
&timize your wave function at a k_s such that the wave function is rea&
&l or (iii) ask Neil to implement Varmin_linjas for complex wave funct&
&ions.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(opt_cycles<1.and.(isopt_vmc.or.isvmc_opt))then
call wout('OPT_CYCLES input parameter should clearly be greater than 0&
&.')
call wout()
xyzzyaaad9=.true.
endif
if(makemovie.and..not.isvmc)then
call wordwrap('Can only make movies using VMC at present.  So MAKEMOVI&
&E should only be T if RUNTYPE is vmc in input.')
call wout()
xyzzyaaad9=.true.
endif
if(isvmc_opt.and.xyzzyaaae1<1)then
call wordwrap('Must write out some configurations to perform an optimi&
&zation.  You should increase NWRCON in input.')
call wout()
xyzzyaaad9=.true.
endif
if(isvmc_dmc.and.xyzzyaaae1<1)then
call wordwrap('Must write out some configurations to perform a DMC cal&
&culation. You should increase NWRCON in input.')
call wout()
xyzzyaaad9=.true.
endif
if(use_orbmods.and.trim(atom_basis_type)/='numerical'.and.trim(atom_ba&
&sis_type)/='gaussian'.and.trim(atom_basis_type)/='slater-type')then
call wordwrap('USE_ORBMODS currently only implemented for ATOM_BASIS_T&
&YPE=''gaussian'', ''slater-type'' or ''numerical''.')
call wout()
xyzzyaaad9=.true.
endif
if((isopt.or.isopt_vmc.or.isvmc_opt).and.any(xyzzyaadd1).and.(trim(ato&
&m_basis_type)=='numerical'.or.trim(atom_basis_type)=='gaussian'.or.tr&
&im(atom_basis_type)=='slater-type').and..not.use_orbmods)then
call wordwrap('The USE_ORBMODS keyword in input needs to be set to T w&
&hen optimizing Gaussian orbitals or numerical-orbital modification fu&
&nctions (OPT_ORBITALS=T with ATOM_BASIS_TYPE=''numerical'', ''gaussia&
&n'' or ''slater-type'').')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaar1<-1.or.xyzzyaaar1>3)then
call wout('BLIP_PERIODICITY should be -1, 0, 1, 2 or 3.')
call wout('-1 => use default.')
call wout(' 0 => finite system.')
call wout(' 1 => polymer, periodic in x direction.')
call wout(' 2 => slab, periodic in x and y directions.')
call wout(' 3 => crystalline solid, periodic in x, y and z directions.&
&')
call wout()
xyzzyaaad9=.true.
endif
if(trim(atom_basis_type)=='blip')then
if(isperiodic.and.xyzzyaaar1==0)then
call wout('BLIP_PERIODICITY must be 1, 2 or 3 if PERIODIC=T in input.'&
&)
call wout()
xyzzyaaad9=.true.
endif
if(.not.isperiodic.and.xyzzyaaar1>=1)then
call wout('BLIP_PERIODICITY must be 0 if PERIODIC=F in input.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(eval_dipole_moment.and.isperiodic)then
call wout('The dipole moment can only be accumulated for nonperiodic s&
&ystems.')
call wout()
xyzzyaaad9=.true.
endif
if(eval_contact_den)then
if(isperiodic)then
call wordwrap('The electron-positron contact density can only be accum&
&ulated for nonperiodic systems at present.')
call wout()
xyzzyaaad9=.true.
endif
if(nspin/=3)then
call wordwrap('The electron-positron contact density can only be accum&
&ulated for positronic molecule systems.')
call wout()
xyzzyaaad9=.true.
endif
if(any(pcharge(1:2)/=-1.d0).or.pcharge(3)/=1.d0.or.any(pmass(1:3)/=1.d&
&0))then
call wordwrap('The electron-positron contact density can only be accum&
&ulated for positronic molecule systems.')
call wout()
xyzzyaaad9=.true.
endif
endif
if((isopt_vmc.or.isvmc_opt.or.isopt).and.any(index(xyzzyaacz1,'emin')=&
&=1).and.(emin_xi_value>1.d0.or.emin_xi_value<0.d0))then
call wout('In energy minimization, the xi parameter must take a value &
&between 0 and 1. Reset EMIN_XI_VALUE to a valid value.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaacv1.and.(isvmc.or.isvmc_dmc.or.isvmc_opt.or.isopt_vmc))then
if(.not.complex_wf)then
call wordwrap('Monte Carlo twist averaging requires a complex wave fun&
&ction.  Set COMPLEX_WF to T.')
call wout()
xyzzyaaad9=.true.
endif
if(trim(atom_basis_type)/='none')then
call wout('Monte Carlo twist averaging can only be done for HEG fluids&
&.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaaw1<=0)then
call wordwrap('If Monte Carlo twist averaging is to be performed in VM&
&C then NEQUIL_TA must be given a positive value.')
call wout()
xyzzyaaad9=.true.
endif
if(pair_corr)then
call wordwrap('VMC twist averaging not implemented for PCF calculation&
&s with a fixed particle.')
call wout()
xyzzyaaad9=.true.
endif
if(mom_den)then
call wordwrap('VMC Monte Carlo twist averaging cannot be performed in &
&conjunction with momentum-density accumulation.')
call wout()
xyzzyaaad9=.true.
endif
if(twop_dm_mom)then
call wordwrap('VMC Monte Carlo twist averaging cannot be performed in &
&conjunction with two-particle density matrix accumulation in momentum&
& space.')
call wout()
xyzzyaaad9=.true.
endif
if(cond_frac_mom)then
call wordwrap('VMC Monte Carlo twist averaging cannot be performed in &
&conjunction with condensate fraction estimator accumulation in moment&
&um space.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(xyzzyaaav1>0.and.(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvmc_dmc_equi&
&l))then
if(.not.complex_wf)then
call wordwrap('Monte Carlo twist averaging requires a complex wave fun&
&ction.  Set COMPLEX_WF=T.')
call wout()
xyzzyaaad9=.true.
endif
if(trim(atom_basis_type)/='none')then
call wout('Monte Carlo twist averaging can only be done for HEG fluids&
&.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaaat1<=0.or.xyzzyaaau1<=0)then
call wordwrap('If Monte Carlo twist averaging is to be performed in DM&
&C then NMOVE_DMCT_EQUIL and NBLOCK_DMCT_EQUIL must be given positive &
&values.')
call wout()
xyzzyaaad9=.true.
endif
if(pair_corr)then
call wordwrap('DMC twist averaging not implemented for PCF calculation&
&s with a fixed particle.')
call wout()
xyzzyaaad9=.true.
endif
if(mom_den)then
call wordwrap('DMC Monte Carlo twist averaging cannot be performed in &
&conjunction with momentum-density accumulation.')
call wout()
xyzzyaaad9=.true.
endif
if(twop_dm_mom)then
call wordwrap('DMC Monte Carlo twist averaging cannot be performed in &
&conjunction with two-particle density matrix accumulation in momentum&
& space.')
call wout()
xyzzyaaad9=.true.
endif
if(cond_frac_mom)then
call wordwrap('DMC Monte Carlo twist averaging cannot be performed in &
&conjunction with condensate fraction estimator accumulation in moment&
&um space.')
call wout()
xyzzyaaad9=.true.
endif
elseif(xyzzyaaav1<0)then
call wout('Number of DMC twists (NUM_DMC_TWISTS) must be non-negative.&
&')
call wout()
xyzzyaaad9=.true.
endif
if(virtual_nconfig<0)then
call wout('VIRTUAL_NCONFIG input parameter must be 0 or positive.')
call wout()
xyzzyaaad9=.true.
endif
if(finite_size_corr.and..not.isperiodic)then
call wordwrap('Evaluation of finite size correction flagged, but this &
&is a finite system. Turn off FINITE_SIZE_CORR keyword in input.')
call wout()
xyzzyaaad9=.true.
endif
if(finite_size_corr.and.trim(interaction)=='manual')then
call wordwrap('Finite size correction not defined for INTERACTION=manu&
&al')
call wout()
xyzzyaaad9=.true.
endif
if(finite_size_corr.and.trim(interaction)=='none')then
call wordwrap('Finite size correction not defined for INTERACTION=none&
&')
call wout()
xyzzyaaad9=.true.
endif
if(periodicity>dimensionality.or.periodicity<0.or.dimensionality>3)the&
&n
call wout('PERIODICITY   : '//trim(i2s(periodicity)))
call wout('DIMENSIONALITY: '//trim(i2s(dimensionality)))
call wout('Should have 0<=periodicity<=dimensionality<=3.')
call wout()
xyzzyaaad9=.true.
endif
if(dimensionality==1.and.nspin>1.and.harmwire_b<0.d0)then
ispin_check : do xyzzyaaaa9=1,nspin-1
if(nele(xyzzyaaaa9)>0)then
do xyzzyaaab9=xyzzyaaaa9+1,nspin
if(nele(xyzzyaaab9)>0.and.heg_layer(xyzzyaaaa9)==heg_layer(xyzzyaaab9)&
&)then
call wordwrap('The wave function of interacting particles must be zero&
& at coalescence points.  Hence only ferromagnetic HEGs are allowed in&
& one dimension.  Multiple particle types are only allowed in differen&
&t wires (layers).')
call wout()
xyzzyaaad9=.true.
exit ispin_check
endif
enddo
endif
enddo ispin_check
endif
if(use_gpcc)then
if(cusp_correction)then
call wordwrap('Cannot use both GP cusp correction and Gaussian cusp co&
&rrector.  Set either USE_GPCC or CUSP_CORRECTION to F in input.')
call wout()
xyzzyaaad9=.true.
endif
if(trim(atom_basis_type)/="plane-wave".and.trim(atom_basis_type)/="bli&
&p" .and.trim(atom_basis_type)/="gaussian")then
call wordwrap('Should only use GP cusp correction for plane-wave, blip&
& or Gaussian bases.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(ebest_av_window<=0)then
call wout('EBEST_AV_WINDOW must be positive.')
call wout()
xyzzyaaad9=.true.
endif
if((isvmc_dmc.or.isdmc.or.isdmc_dmc).and.(corper_dmc<1.or.corper_dmc>x&
&yzzyaaal1*xyzzyaaan1*ndmcave))then
call wout('Must have 0 < CORPER_DMC <= NMOVE_DMC_STATS*NBLOCK_DMC_STAT&
&S*NDMCAVE.')
call wout()
xyzzyaaad9=.true.
endif
if(ndmcave<1)then
call wout('Must have NDMCAVE>0.')
call wout()
xyzzyaaad9=.true.
endif
if(isgen_mpc.and..not.isperiodic)then
call wout('No need to generate mpc.data for aperiodic calculations.')
call wout()
xyzzyaaad9=.true.
endif
if(dmc_reweight_configs.and..not.lwdmc)then
call wout('Need to set LWDMC to T if DMC_REWEIGHT_CONF is T in input.'&
&)
call wout()
xyzzyaaad9=.true.
endif
if(lwdmc_fixpop)then
lwdmc=.true.
endif
if(isvmc_dmc.or.isdmc.or.isdmc_dmc.or.isvmc_dmc_equil)then
if(lwdmc.and.nc_dmc)then
call wordwrap('Cannot perform norm-conserving weighted DMC.  Please se&
&t either DMC_NORM_CONSERVE or LWDMC to F in input.')
call wout()
xyzzyaaad9=.true.
endif
if(lwdmc.and.poprenorm)then
call wordwrap('Cannot perform population renormalization with weighted&
& DMC.  Please set either DMC_POPRENORM or LWDMC to F in input.')
call wout()
xyzzyaaad9=.true.
endif
if(nconfig_prelim>0)then
if(xyzzyaaav1>0)then
call wordwrap('Cannot perform a preliminary DMC calculation when twist&
& averaging.')
call wout()
xyzzyaaad9=.true.
endif
if(targ_wt*dble(xyzzyaaan1*xyzzyaaal1*ndmcave/corper_dmc) <dble(nconfi&
&g_prelim)*dble(nnodes))then
call wordwrap('Insufficient statistics-accumulation steps are requeste&
&d in the preliminary DMC run to generate the required population.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(redist_grp_size<1)then
call wordwrap('REDIST_GRP_SIZE must be a positive integer.')
call wout()
xyzzyaaad9=.true.
endif
if(redist_grp_size<500)then
call wordwrap('Setting REDIST_GRP_SIZE to less than 500 is neither adv&
&isable or necessary - delete this error trap in src/monte_carlo.f90 i&
&f you know what you are doing and you really want to.')
call wout()
xyzzyaaad9=.true.
endif
endif
if(isrmc.or.isrmc_rmc)then
if(xyzzyaack1<=0.d0)then
call wout('Must have DTRMC>0 in input.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaabw1<=4)then
call wout('RMC not supported with RMC_REP_LENGTH<=4 in input.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaabx1<=0)then
call wout('Must have RMC_MOVE_LENGTH>0 in input.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaadr1<=0)then
call wout('Must have CORPER_RMC>0 in input.')
call wout()
xyzzyaaad9=.true.
endif
if(iaccum.or.isrmc_rmc)then
if(xyzzyaadu1<=0)then
call wout('Must have NMOVE_RMC_STATS>0 in input for stats accumulation&
&.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaadv1<=0)then
call wout('Must have NBLOCK_RMC_STATS>0 in input for stats accumulatio&
&n.')
call wout()
xyzzyaaad9=.true.
endif
endif
if((.not.iaccum).or.isrmc_rmc)then
if(xyzzyaads1<=0)then
call wout('Must have NMOVE_RMC_EQUIL>0 in input for equilibration.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaadt1<=0)then
call wout('Must have NBLOCK_RMC_EQUIL>0 in input for equilibration.')
call wout()
xyzzyaaad9=.true.
endif
endif
endif
if(xyzzyaaah1<0.or.xyzzyaaai1<0)then
call wout('Number of up- and down-spin holes NHU and NHD in input shou&
&ld not be negative.')
call wout()
xyzzyaaad9=.true.
endif
if(xyzzyaabp1/=1)then
call wout('DMC_NPOPS keyword redundant - since there is no point inusi&
&ng this feature. Please remove this keyword from the input file.')
call wout()
xyzzyaaad9=.true.
endif
if(chkpoint_level<-1.or.chkpoint_level>2)then
call wout('Value of CHECKPOINT in input must be -1, 0, 1 or 2.')
call wout()
xyzzyaaad9=.true.
endif
if(shm_size_nproc<0)then
call wout('Invalid value of SHM_SIZE_NPROC in input; must be a positiv&
&einteger')
call wout()
xyzzyaaad9=.true.
endif
if(on_top_ii/=0)then
if(.not.isvmc)then
call wordwrap('ON_TOP_PAIR is only supported for VMC runs.')
xyzzyaaad9=.true.
endif
if(.not.mom_den)then
call errwarn('CHECK_INPUT_PARAMETERS','ON_TOP_PAIR is only useful for &
&VMC on-top momentum density calculations, but MOM_DEN has not been se&
&t.  Make sure you know what you are doing.  Continuing run.')
endif
call errwarn('PRINT_INPUT_PARAMETERS','ON_TOP_PAIR activated - ignore &
&energies from this run.')
endif
if(trim(atom_basis_type)/='none'.and..not.electron_system .and.xyzzyaa&
&ah1==0.and.xyzzyaaai1==0)call errwarn('CHECK_INPUT_PARAMETERS','Elect&
&rons have been redefined in PARTICLES block in input, and ATOM_BASIS_&
&TYPE is not ''none''. This is not expected to happen. Make sure you k&
&now what you are doing. Continuing run.')
if(limdmc/=4.and.(xyzzyaaae1/=0.or.isdmc.or.isvmc_dmc.or.isdmc_dmc.or.&
&isvmc_dmc_equil))call errwarn_silent('CHECK_INPUT_PARAMETERS','It is &
&recommended that the LIMDMC parameter be set to 4. Its current value &
&is '//trim(i2s(limdmc))//'.')
if(localized_orbitals.and.bsmooth)call errwarn('CHECK_INPUT_PARAMETERS&
&','It is recommended that BSMOOTH be set to FALSE.')
if(xyzzyaacu1)call errwarn('CHECK_INPUT_PARAMETERS','Setting TIMING_IN&
&FO=T may slow down the calculation. For best performance, set TIMING_&
&INFO=F.')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvmc_dmc_equil)then
if(targ_wt<5.d0.and.nnodes==1)call errwarn('CHECK_INPUT_PARAMETERS','I&
&nitial number of configs (DMC_TARGET_WEIGHT) is very small. This is i&
&nadvisable in single-processor DMC - the code will likely stop as soo&
&n as the natural population fluctuations lead to all the configs dyin&
&g out.')
if(nnodes>1.and.targ_wt<nnodes)call errwarn('CHECK_INPUT_PARAMETERS','&
&The initial number of configs (DMC_TARGET_WEIGHT) is less than  the n&
&umber of processors. This is highly inadvisable.')
if(trip_popn>0.d0.and.chkpoint_level<1)then
call errwarn('CHECK_INPUT_PARAMETERS','CHECKPOINT must be > 0 for DMC_&
&TARGET_WEIGHT > 0. Setting value of CHECKPOINT to 1.')
chkpoint_level=1
endif
if(nconfig_prelim>0.and.((dble(corper_dmc)*xyzzyaace1<1.d0.and.corper_&
&dmc<256) .or.corper_dmc<16))call errwarn_silent('CHECK_INPUT_PARAMETE&
&RS', 'CORPER_DMC is rather small for preliminary DMC.')
endif
if(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvmc.or.isvmc_dmc_equil)then
if(chkpoint_level==0.and.max_cpu_time==0.d0.and.max_real_time==0.d0)th&
&en
call errwarn('CHECK_INPUT_PARAMETERS','If CHECKPOINT=0 and your job ha&
&s a time limit, then it is recommended that you use the MAX_CPU_TIME &
&or MAX_REAL_TIME keywords to ensure that the config.out file is writt&
&en if the required time exceeds the limit.')
endif
if(chkpoint_level==-1)then
call errwarn('CHECK_INPUT_PARAMETERS','CHECKPOINT=-1 in input. Be awar&
&e this value should be chosen only if you *know* that the job will fi&
&t in any imposed time limit, and that such a run will be long enough &
&to give an acceptably small error bar, since it will be impossible to&
& subsequently continue the run.')
endif
endif
if(isopt.or.isvmc_opt.or.isopt_vmc)then
if(use_jastrow.and.any(xyzzyaadb1.and..not.xyzzyaada1))call errwarn_si&
&lent('CHECK_INPUT_PARAMETERS','If you optimize the backflow function &
&it is highly recommended that you also optimize the Jastrow factor.  &
&Set both OPT_JASTROW and OPT_BACKFLOW to T in the input file.')
endif
if(xyzzyaacw1.and.xyzzyaabx1>1)call errwarn('CHECK_INPUT_PARAMETERS','&
&If using bounce algorithm, you probably want to set RMC_MOVE_LENGTH t&
&o 1.')
else
if(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvmc_dmc_equil)then
if(trip_popn>0.d0.and.chkpoint_level<1)chkpoint_level=1
endif
endif
call mpi_bcast(xyzzyaaad9,1,mpi_logical,0,mpi_comm_world,ierror)
if(xyzzyaaad9)call errstop_master('CHECK_INPUT_PARAMETERS','Stopping.'&
&)
if(checkpoint_ncpu>nnodes.or.checkpoint_ncpu<1)checkpoint_ncpu=nnodes
end subroutine xyzzyaafc1
subroutine xyzzyaafd1
implicit none
if(model_system.and.nbasis==0.and.trim(heg_crystal_type)=='unset')then
homogeneous_system=all(heg_orbtype/=6)
else
homogeneous_system=.false.
endif
if(am_master.and..not.xyzzyaacu1)use_timer=.false.
if(tpdmc==9999)tpdmc=nint(10.d0/xyzzyaace1)
if(.not.isperiodic.and.npcells/=1)then
scell_matrix(1,1)=1
scell_matrix(2,2)=1
scell_matrix(3,3)=1
npcells=1
endif
if(xyzzyaach1/=0.d0)constant_energy=constant_energy-xyzzyaach1*netot
if(isvmc_dmc.or.isvmc_dmc_equil.or.isvmc_opt.or.isopt_vmc)writeout_vmc&
&_hist=.false.
if(.not.(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvmc_dmc_equil))use_tmove&
&=.false.
vmc_cfg_by_cfg=xyzzyaaao1==3
dmc_cfg_by_cfg=xyzzyaaap1==2
mc_twist_av=((xyzzyaacv1.and.(isvmc.or.isvmc_opt.or.isopt_vmc .or.isvm&
&c_dmc)).or.(xyzzyaaav1>0.and.(isdmc.or.isvmc_dmc.or.isdmc_dmc.or.isvm&
&c_dmc_equil)))
use_altsamp=.false.
altsamp=0
simplepdf=0
ndet_smp=0
select case(trim(vmc_sampling))
case('optimum')
use_altsamp=.true.
altsamp=1
simplepdf=0
case('HF optimum')
use_altsamp=.true.
altsamp=1
simplepdf=1
case('efficient')
use_altsamp=.true.
altsamp=2
simplepdf=1
case('new')
use_altsamp=.true.
altsamp=2
simplepdf=2
case('stdtest')
use_altsamp=.true.
altsamp=3
simplepdf=0
end select
end subroutine xyzzyaafd1
subroutine xyzzyaafe1(pname,pval)
implicit none
character(*),intent(in) :: pname,pval
character(80) tmpr
write(tmpr,'(a,t42,'':  '',a)')trim(pname),trim(pval)
call wout(tmpr)
end subroutine xyzzyaafe1
subroutine xyzzyaaff1
implicit none
character(80) tmpr,tmpr2
call wout('General input parameters')
call wout('========================')
call xyzzyaafe1('NEU (num up spin electrons)',i2s(xyzzyaaaf1))
call xyzzyaafe1('NED (num down spin electrons)',i2s(xyzzyaaag1))
if(xyzzyaaah1>0.or.xyzzyaaai1>0)then
call xyzzyaafe1('NHU (num up spin holes)',i2s(xyzzyaaah1))
call xyzzyaafe1('NHD (num down spin holes)',i2s(xyzzyaaai1))
endif
call xyzzyaafe1('RUNTYPE (type of run)',runtype)
call xyzzyaafe1('PSI_S  (form for [anti]symmetrizing wfn)',psi_s)
call xyzzyaafe1('ATOM_BASIS_TYPE (atom-centred orb basis)',atom_basis_&
&type)
call xyzzyaafe1('INTERACTION (interaction type)',interaction)
call xyzzyaafe1('TESTRUN (read input data,print and stop)',l2s(xyzzyaa&
&cm1))
call xyzzyaafe1('PERIODIC',l2s(isperiodic))
if(isperiodic)call xyzzyaafe1('NPCELLS (num primitive cells)',i2s(npce&
&lls))
call xyzzyaafe1('COMPLEX_WF (complex Slater wave fn.)',l2s(complex_wf)&
&)
call xyzzyaafe1('NEIGHPRINT (neighbour analysis)',i2s(neighprint))
call xyzzyaafe1('USE_JASTROW (use Jastrow factor)',l2s(use_jastrow))
call xyzzyaafe1('BACKFLOW (use backflow corrections)',l2s(use_backflow&
&))
if(trim(psi_s)=='slater')call xyzzyaafe1('DBARRC (DBAR recalculation p&
&eriod)',i2s(dbarrc))
if(trim(atom_basis_type)=='numerical')call xyzzyaafe1('USE_ORBMODS (us&
&e orbital modifications)',l2s(use_orbmods))
if(trim(atom_basis_type)=='gaussian')then
call xyzzyaafe1('USE_ORBMODS (use orbitals modifications)',l2s(use_orb&
&mods))
call xyzzyaafe1('CUSP_CORRECTION',l2s(cusp_correction))
call xyzzyaafe1('MOLGSCREENING',l2s(molgscreening))
endif
call xyzzyaafe1('NON_LOCAL_GRID (NL integration grid)',i2s(xyzzyaaas1)&
&)
tmpr=r2s(xyzzyaach1,'(f7.4)')
call xyzzyaafe1('E_OFFSET (energy offset)',tmpr)
call xyzzyaafe1('ESUPERCELL',l2s(esupercell))
if(interaction_mpc_present)then
call xyzzyaafe1('QMC_DENSITY_MPC',l2s(qmc_density_mpc))
call xyzzyaafe1('PERMIT_DEN_SYMM',l2s(permit_den_symm))
call xyzzyaafe1('BLIP_MPC',l2s(blip_mpc))
endif
if(trim(atom_basis_type)=='gaussian')then
tmpr=r2s(gautol,'(f4.1)')
call xyzzyaafe1('GAUTOL  (Gaussian evaluation tolerance)',tmpr)
elseif(trim(atom_basis_type)=='blip')then
if(xyzzyaaeq1)call xyzzyaafe1('SINGLE_PRECISION_BLIPS',l2s(single_prec&
&ision_blips))
if(xyzzyaaer1)call xyzzyaafe1('SP_BLIPS',l2s(single_precision_blips))
call xyzzyaafe1('WRITE_BINARY_BLIPS (create bwfn.data.bin)',l2s(write_&
&binary_blips))
call xyzzyaafe1('CONV_BINARY_BLIPS (blip b1 --> bin)',l2s(conv_binary_&
&blips))
elseif(allocated(heg_orbtype))then
if(any(heg_orbtype==2))then
tmpr=r2s(gautol,'(f4.1)')
call xyzzyaafe1('GAUTOL (Gaussian evaluation tolerance)',tmpr)
endif
if(any(heg_orbtype==6))call xyzzyaafe1('FIX_HOLES (constraint for BIEX&
&3)',l2s(fix_holes))
endif
call xyzzyaafe1('SPARSE',l2s(sparse))
if(sparse)then
tmpr=r2s(sparse_threshold,'(en10.1)')
call xyzzyaafe1('SPARSE_THRESHOLD',tmpr)
endif
if(localized_orbitals)call xyzzyaafe1('BSMOOTH (smoothing of loc. orbi&
&tal edge)',l2s(bsmooth))
if(.not.isperiodic)call xyzzyaafe1('DIPOLE_MOMENT',l2s(eval_dipole_mom&
&ent))
call xyzzyaafe1('CHECKPOINT (checkpointing level)',i2s(chkpoint_level)&
&)
call xyzzyaafe1('CHECKPOINT_NCPU (chkpnt group size)',i2s(checkpoint_n&
&cpu))
call xyzzyaafe1('CON_LOC (Dir to read/write config.*)',con_loc)
call xyzzyaafe1('RELATIVISTIC',l2s(relativistic))
if(on_top_ii/=0)then
call xyzzyaafe1('ON_TOP_PAIR - i',i2s(on_top_ii))
call xyzzyaafe1('ON_TOP_PAIR - j',i2s(on_top_jj))
endif
call wout()
if(isvmc)then
call wout('VMC input parameters')
call wout('====================')
elseif(isdmc.or.isdmc_dmc)then
call wout('DMC input parameters')
call wout('====================')
elseif(isvmc_dmc)then
call wout('VMC/DMC input parameters')
call wout('========================')
elseif(isrmc.or.isrmc_rmc)then
call wout('RMC input parameters')
call wout('====================')
elseif(isopt)then
call wout('Optimization input parameters')
call wout('=============================')
elseif(isvmc_opt.or.isopt_vmc)then
call wout('VMC/optimization input parameters')
call wout('=================================')
endif
if(.not.(isopt.or.isopt_vmc.or.isgen_mpc.or.xyzzyaaen1))call xyzzyaafe&
&1('NEWRUN (start new run)',l2s(newrun))
if(isvmc.or.isvmc_opt.or.isopt_vmc.or.isvmc_dmc.or.isvmc_dmc_equil)the&
&n
call xyzzyaafe1('VMC_METHOD (choice of VMC algorithm)',i2s(xyzzyaaao1)&
&)
tmpr=r2s(xyzzyaacd1,'(es12.4)')
call xyzzyaafe1('DTVMC (VMC time step)',tmpr)
if(xyzzyaacl1>0.d0)then
tmpr=r2s(xyzzyaacl1,'(es12.4)')
call xyzzyaafe1('DTVMC_SHIFT (rel. tr.prob.dist. shift)',tmpr)
endif
call xyzzyaafe1('OPT_DTVMC (VMC time-step optimization)',i2s(xyzzyaaaq&
&1))
if(use_blocktime)then
tmpr=r2ss(block_time,'(es16.6)')
call xyzzyaafe1('BLOCK_TIME (time per block in sec.)',tmpr)
endif
if(old_input)then
if(nnodes==1.or.trim(parallel_keywords)=='per_node')then
call xyzzyaafe1('NMOVE (num moves per block)',i2s(xyzzyaaab1))
call xyzzyaafe1('NWRCON (num configs to write)',i2s(xyzzyaaae1))
else
call xyzzyaafe1('NMOVE (num moves per block)',trim(i2s(xyzzyaaab1*nnod&
&es))//' ('//trim(i2s(xyzzyaaab1))//' per node)')
call xyzzyaafe1('NWRCON (num configs to write)',trim(i2s(xyzzyaaae1*nn&
&odes))//' ('//trim(i2s(xyzzyaaae1))//' per node)')
endif
call xyzzyaafe1('NBLOCK (num blocks)',i2s(xyzzyaaac1))
call xyzzyaafe1('NEQUIL (num equilibration steps)',i2s(xyzzyaaaa1))
if(xyzzyaacv1)then
call xyzzyaafe1('VMC_TWIST_AV',l2s(xyzzyaacv1))
call xyzzyaafe1('NEQUIL_TA (num equil moves for TA)',i2s(xyzzyaaaw1))
endif
if(xyzzyaaad1==0)then
call xyzzyaafe1('CORPER (correlation period)','0 (automatic)')
else
call xyzzyaafe1('CORPER (correlation period)',i2s(xyzzyaaad1))
endif
call xyzzyaafe1('NVMCAVE (num energies to ave betw write)',i2s(xyzzyaa&
&aj1))
call wout('Equivalents in new keyword set:')
call xyzzyaafe1('  VMC_NSTEP',i2s(xyzzyaaay1))
call xyzzyaafe1('  VMC_NCONFIG_WRITE',i2s(xyzzyaaaz1))
call xyzzyaafe1('  VMC_NBLOCK',i2s(xyzzyaabc1))
call xyzzyaafe1('  VMC_EQUIL_NSTEP',i2s(xyzzyaabd1))
if(xyzzyaacv1)then
call xyzzyaafe1('  VMC_NTWIST',i2s(xyzzyaabn1))
call xyzzyaafe1('  VMC_REEQUIL_NSTEP',i2s(xyzzyaabk1))
endif
if(xyzzyaabb1==0)then
call xyzzyaafe1('  VMC_DECORR_PERIOD','0 (automatic)')
else
call xyzzyaafe1('  VMC_DECORR_PERIOD',i2s(xyzzyaabb1))
endif
call xyzzyaafe1('  VMC_AVE_PERIOD',i2s(xyzzyaaba1))
else
call xyzzyaafe1('VMC_NSTEP (num VMC steps)',i2s(xyzzyaaay1))
call xyzzyaafe1('VMC_NCONFIG_WRITE (num configs to write)',i2s(xyzzyaa&
&az1))
if(.not.use_blocktime)call xyzzyaafe1('VMC_NBLOCK (num VMC blocks)',i2&
&s(xyzzyaabc1))
call xyzzyaafe1('VMC_EQUIL_NSTEP (num equil steps)',i2s(xyzzyaabd1))
if(xyzzyaacv1)then
call xyzzyaafe1('VMC_NTWIST (num twist angles)',i2s(xyzzyaabn1))
call xyzzyaafe1('VMC_REEQUIL_NSTEP (re-equil steps for TA)',i2s(xyzzya&
&abk1))
endif
if(xyzzyaabb1==0)then
call xyzzyaafe1('VMC_DECORR_PERIOD (length of inner loop)','0 (automat&
&ic)')
else
call xyzzyaafe1('VMC_DECORR_PERIOD (length of inner loop)',i2s(xyzzyaa&
&bb1))
endif
call xyzzyaafe1('VMC_AVE_PERIOD (hist reduction factor)',i2s(xyzzyaaba&
&1))
endif
call xyzzyaafe1('VMC_SAMPLING',vmc_sampling)
if(altsamp==1)then
tmpr=r2s(vmc_optimum_e0,'(f10.5)')
call xyzzyaafe1('VMC_OPTIMUM_E0',tmpr)
tmpr=r2s(vmc_optimum_ew,'(f10.5)')
call xyzzyaafe1('VMC_OPTIMUM_EW',tmpr)
endif
endif
if(isvmc_dmc.or.isdmc.or.isdmc_dmc.or.isvmc_dmc_equil)then
if(old_input)then
tmpr=r2s(targ_wt/dble(nnodes),'(f10.2)')
if(nnodes==1.or.trim(parallel_keywords)=='per_node')then
call xyzzyaafe1('NCONFIG (target weight)',tmpr)
else
tmpr2=r2s(targ_wt,'(f10.2)')
call xyzzyaafe1('NCONFIG (target weight)',trim(tmpr2)//' ('//trim(tmpr&
&)//' per node)')
endif
if(.not.(isdmc.and.iaccum))then
call xyzzyaafe1('NMOVE_DMC_EQUIL (num moves/block)',i2s(xyzzyaaak1))
if(.not.use_blocktime)call xyzzyaafe1('NBLOCK_DMC_EQUIL (num blocks)',&
&i2s(xyzzyaaam1))
endif
if(.not.(isdmc.and..not.iaccum).and..not.isvmc_dmc_equil)then
call xyzzyaafe1('NMOVE_DMC_STATS (num moves/block)',i2s(xyzzyaaal1))
if(.not.use_blocktime)call xyzzyaafe1('NBLOCK_DMC_STATS (num blocks)',&
&i2s(xyzzyaaan1))
if(xyzzyaaav1>0)then
call xyzzyaafe1('NMOVE_DMCT_EQUIL',i2s(xyzzyaaat1))
call xyzzyaafe1('NBLOCK_DMCT_EQUIL',i2s(xyzzyaaau1))
endif
endif
call xyzzyaafe1('CORPER_DMC',i2s(corper_dmc))
call xyzzyaafe1('NDMCAVE',i2s(ndmcave))
tmpr=r2s(trip_popn,'(f10.2)')
call xyzzyaafe1('TRIP_POPN (catastrophe thres)',tmpr)
call wout('Equivalents in new keyword set:')
tmpr=r2s(targ_wt,'(f10.2)')
call xyzzyaafe1('  DMC_TARGET_WEIGHT',tmpr)
if(.not.(isdmc.and.iaccum))then
call xyzzyaafe1('  DMC_EQUIL_NSTEP',i2s(xyzzyaabe1))
if(.not.use_blocktime)call xyzzyaafe1('  DMC_EQUIL_NBLOCK',i2s(xyzzyaa&
&bf1))
endif
if(.not.(isdmc.and..not.iaccum))then
call xyzzyaafe1('  DMC_STATS_NSTEP',i2s(xyzzyaabg1))
if(.not.use_blocktime)call xyzzyaafe1('  DMC_STATS_NBLOCK',i2s(xyzzyaa&
&bh1))
if(xyzzyaaav1>0)then
call xyzzyaafe1('  DMC_NTWIST',i2s(xyzzyaabo1))
call xyzzyaafe1('  DMC_REEQUIL_NSTEP',i2s(xyzzyaabl1))
call xyzzyaafe1('  DMC_REEQUIL_NBLOCK',i2s(xyzzyaabm1))
endif
endif
call xyzzyaafe1('  DMC_DECORR_PERIOD',i2s(xyzzyaabj1))
call xyzzyaafe1('  DMC_AVE_PERIOD',i2s(xyzzyaabi1))
tmpr=r2s(xyzzyaacj1,'(f10.2)')
call xyzzyaafe1('  DMC_TRIP_WEIGHT',tmpr)
else
tmpr=r2s(targ_wt,'(f10.2)')
call xyzzyaafe1('DMC_TARGET_WEIGHT',tmpr)
call xyzzyaafe1('DMC_MD',l2s(dmc_md))
if(use_blocktime)then
tmpr=r2ss(block_time,'(es16.6)')
call xyzzyaafe1('BLOCK_TIME (time per block in sec.)',tmpr)
endif
if(.not.(isdmc.and.iaccum))then
call xyzzyaafe1('DMC_EQUIL_NSTEP (num equil steps)',i2s(xyzzyaabe1))
if(.not.use_blocktime)call xyzzyaafe1('DMC_EQUIL_NBLOCK (num blocks)',&
&i2s(xyzzyaabf1))
if(dmc_md)call xyzzyaafe1('DMCMD_EQUIL_NSTEP (num equil steps)',i2s(xy&
&zzyaabe1))
endif
if(.not.(isdmc.and..not.iaccum).and..not.isvmc_dmc_equil)then
call xyzzyaafe1('DMC_STATS_NSTEP (num stats steps)',i2s(xyzzyaabg1))
if(.not.use_blocktime)call xyzzyaafe1('DMC_STATS_NBLOCK (num blocks)',&
&i2s(xyzzyaabh1))
if(dmc_md)call xyzzyaafe1('DMCMD_STATS_NSTEP (num stats steps)',i2s(xy&
&zzyaabg1))
if(xyzzyaaav1>0)then
call xyzzyaafe1('DMC_NTWIST (num twist angles)',i2s(xyzzyaabo1))
call xyzzyaafe1('DMC_REEQUIL_NSTEP (re-equil steps for TA)',i2s(xyzzya&
&abl1))
call xyzzyaafe1('DMC_REEQUIL_NBLOCK (re-equil blocks TA)',i2s(xyzzyaab&
&m1))
endif
endif
call xyzzyaafe1('DMC_DECORR_PERIOD (length of inner loop)',i2s(xyzzyaa&
&bj1))
call xyzzyaafe1('DMC_AVE_PERIOD (hist reduction factor)',i2s(xyzzyaabi&
&1))
tmpr=r2s(xyzzyaacj1,'(f10.2)')
call xyzzyaafe1('DMC_TRIP_WEIGHT (catastrophe thres)',tmpr)
endif
if(nconfig_prelim>0)call xyzzyaafe1('NCONFIG_PRELIM (DMC config gen)',&
&trim(i2s(nconfig_prelim)))
if(xyzzyaaby1>0)call xyzzyaafe1('DMC_NCONF_PRELIM (DMC config gen)',tr&
&im(i2s(xyzzyaaby1)))
call xyzzyaafe1('EBEST_AV_WINDOW (running av for energy)',i2s(ebest_av&
&_window))
call xyzzyaafe1('DMC_METHOD (choice of DMC algorithm)',i2s(xyzzyaaap1)&
&)
call xyzzyaafe1('DMC_REWEIGHT_CONF (Update weights)',l2s(dmc_reweight_&
&configs))
call xyzzyaafe1('DMC_SPACEWARPING (adjust e to new wfn)',l2s(dmc_space&
&warping))
call xyzzyaafe1('REDIST_GRP_SIZE (size of redist groups)',i2s(redist_g&
&rp_size))
tmpr=r2s(xyzzyaace1,'(f7.5)')
call xyzzyaafe1('DTDMC (DMC time step)',tmpr)
call xyzzyaafe1('TPDMC (DMC T_p parameter)',i2s(tpdmc))
tmpr=r2s(cerefdmc,'(f6.3)')
call xyzzyaafe1('CEREFDMC (constant for EREF [DMC])',tmpr)
call xyzzyaafe1('LIMDMC (limit type for drift vel/energy)',i2s(limdmc)&
&)
if(limdmc==2.or.limdmc==4)then
call xyzzyaafe1('NUCLEUS_GF_MODS (DMC GF mods for nuclei)',l2s(nucleus&
&_gf_mods))
endif
if(limdmc==3.or.(limdmc==2.and..not.nucleus_gf_mods).or.(limdmc==4.and&
&..not.nucleus_gf_mods))then
tmpr=r2s(alimit,'(f5.3)')
call xyzzyaafe1('ALIMIT',tmpr)
endif
if(trip_popn>0.d0)call xyzzyaafe1('MAX_REC_ATTEMPTS',i2s(max_rec_attem&
&pts))
call xyzzyaafe1('IACCUM (flag for statistics run [DMC])',l2s(iaccum))
call xyzzyaafe1('IBRAN (flag to enable branching [DMC])',l2s(ibran))
call xyzzyaafe1('LWDMC (flag for enabling weighted DMC)',l2s(lwdmc))
call xyzzyaafe1('LWDMC_FIXPOP (fixed population LWDMC)',l2s(lwdmc_fixp&
&op))
if(lwdmc)then
tmpr=r2s(wdmcmin,'(f5.3)')
call xyzzyaafe1('WDMCMIN (min walker weights)',tmpr)
tmpr=r2s(wdmcmax,'(f5.3)')
call xyzzyaafe1('WDMCMAX (max walker weights)',tmpr)
else
call xyzzyaafe1('DMC_NORM_CONSERVE',l2s(nc_dmc))
endif
call xyzzyaafe1('DMC_POPRENORM (renormalize config popn)',l2s(poprenor&
&m))
call xyzzyaafe1('GROWTH_ESTIMATOR (calc growth estimator)',l2s(growth_&
&estimator))
call xyzzyaafe1('USE_TMOVE',l2s(use_tmove))
call xyzzyaafe1('FUTURE WALKING',l2s(use_future))
call xyzzyaafe1('SMALL_TRANSFER (redist. transf. size)',l2s(small_tran&
&sfer))
call xyzzyaafe1('ORBBUF (orbital buffering)',l2s(orbbuf))
call xyzzyaafe1('JASBUF (Jastrow buffering)',l2s(jasbuf))
endif
if(isvmc.or.isvmc_opt.or.isopt_vmc.or.isvmc_dmc.or.isdmc.or.isdmc_dmc.&
&or.isvmc_dmc_equil)then
call xyzzyaafe1('MAKEMOVIE',l2s(makemovie))
if(makemovie)then
call xyzzyaafe1('MOVIEPLOT (on which move to plot)',i2s(movieplot))
call xyzzyaafe1('MOVIENODE (node to plot from)',i2s(movienode))
endif
call xyzzyaafe1('FORCES',l2s(forces))
endif
if(isopt.or.isvmc_opt.or.isopt_vmc)then
if(isvmc_opt.or.isopt_vmc)then
call xyzzyaafe1('OPT_CYCLES (num optimization cycles)',i2s(opt_cycles)&
&)
call xyzzyaafe1('POSTFIT_VMC (perform post-fit VMC calc)',l2s(xyzzyaac&
&r1))
call xyzzyaafe1('POSTFIT_KEEP_CFG (keep post-fit VMC cfgs)',l2s(xyzzya&
&acs1))
call xyzzyaafe1('OPT_NOCTF_CYCLES (fixed cutoff cycles)',i2s(xyzzyaaax&
&1))
endif
call xyzzyaafe1('OPT_INFO (information level)',i2s(opt_info))
call xyzzyaafe1('OPT_JASTROW (opt Jastrow factor)',l2s(opt_jastrow))
if(any(index(xyzzyaacz1,'varmin_linjas')/=1))then
call xyzzyaafe1('OPT_DET_COEFF (opt det coeffs)',l2s(opt_det_coeff))
call xyzzyaafe1('OPT_ORBITALS (opt orbitals)',l2s(opt_orbitals))
call xyzzyaafe1('OPT_BACKFLOW (opt backflow params)',l2s(opt_backflow)&
&)
call xyzzyaafe1('OPT_FIXNL (fix nonlocal energy)',l2s(opt_fixnl))
call xyzzyaafe1('OPT_MAXITER (max num iterations)',i2s(opt_maxiter))
if(any(index(xyzzyaacz1,'emin')/=1))call xyzzyaafe1('OPT_MAXEVAL (max &
&num evaluations)',i2s(opt_maxeval))
call xyzzyaafe1('VM_SMOOTH_LIMITS (smooth limiting)',l2s(vm_smooth_lim&
&its))
if(any(xyzzyaadc1.and..not.xyzzyaadb1))call xyzzyaafe1('OPT_SMALL_BUFF&
&ERS (detcoeff buffer size)',l2s(small_buffers))
endif
if(any(index(xyzzyaacz1,'varmin ')==1))then
call xyzzyaafe1('VM_REWEIGHT (reweighting)',l2s(vm_reweight))
if(vm_reweight.and.vm_w_max>0.d0)then
tmpr=r2s(vm_w_max,'(e14.6)')
call xyzzyaafe1('VM_W_MAX (max weight)',tmpr)
tmpr=r2s(vm_w_min,'(e14.6)')
call xyzzyaafe1('VM_W_MIN (min weight)',tmpr)
endif
call xyzzyaafe1('VM_FILTER (filter outlying configs)',l2s(vm_filter))
if(vm_filter)then
tmpr=r2s(vm_filter_thres,'(e14.6)')
call xyzzyaafe1('VM_FILTER_THRES (filter threshold)',tmpr)
tmpr=r2s(vm_filter_width,'(e14.6)')
call xyzzyaafe1('VM_FILTER_WIDTH (filter width)',tmpr)
endif
endif
if(any(index(xyzzyaacz1,'varmin ')==1).or.any(index(xyzzyaacz1,'madmin&
&')==1))then
call xyzzyaafe1('VM_USE_E_GUESS (use guess energy)',l2s(vm_use_e_guess&
&))
if(vm_use_e_guess)then
tmpr=r2s(vm_e_guess,'(e14.6)')
if(model_system)then
call xyzzyaafe1('VM_E_GUESS (guess E in au per particle)',tmpr)
else
if(isperiodic)then
call xyzzyaafe1('VM_E_GUESS (guess E in au per pr. cell)',tmpr)
else
call xyzzyaafe1('VM_E_GUESS (guess E in au)',tmpr)
endif
endif
endif
endif
if(any(index(xyzzyaacz1,'varmin_linjas')==1))then
call xyzzyaafe1('VM_LINJAS_METHOD',vm_linjas_method)
call xyzzyaafe1('VM_LINJAS_ITS',i2s(vm_linjas_its))
endif
if(any(index(xyzzyaacz1,'emin')==1))then
tmpr=r2s(emin_xi_value,'(f3.1)')
call xyzzyaafe1('EMIN_XI_VALUE (xi parameter)',tmpr)
if(emin_mine_present)then
tmpr=r2s(emin_min_energy,'(f10.5)')
call xyzzyaafe1('EMIN_MIN_ENERGY (energy threshold)',tmpr)
endif
endif
endif
if(isgen_mpc)then
select case(trim(atom_basis_type))
case('plane-wave')
call wout('Generating mpc.data file for geometry defined in pwfn.data.&
&')
case('gaussian')
call wout('Generating mpc.data file for geometry defined in gwfn.data.&
&')
case('slater-type')
call wout('Generating mpc.data file for geometry defined in stowfn.dat&
&a.')
case('blip')
call wout('Generating mpc.data file for geometry defined in bwfn.data.&
&')
case('dimer')
call wout('Generating mpc.data file for geometry defined in dwfn.data.&
&')
case default
call wout('Generating mpc.data file for geometry defined in input.')
end select
elseif(isrmc.or.isrmc_rmc)then
if(iaccum)then
call xyzzyaafe1('NMOVE_RMC_STATS (number of moves/block)',i2s(xyzzyaad&
&u1))
else
call xyzzyaafe1('NMOVE_RMC_EQUIL (number of moves/block)',i2s(xyzzyaad&
&s1))
endif
if(iaccum)then
call xyzzyaafe1('NBLOCK_RMC_STATS (Number of blocks)',i2s(xyzzyaadv1))
else
call xyzzyaafe1('NBLOCK_RMC_EQUIL (Number of blocks)',i2s(xyzzyaadt1))
endif
tmpr=r2s(xyzzyaack1,'(f7.5)')
call xyzzyaafe1('DTRMC (RMC time step)',tmpr)
if(.not.isrmc_rmc)call xyzzyaafe1('IACCUM (flag for statistics run)',l&
&2s(iaccum))
if(xyzzyaadr1>1)then
call xyzzyaafe1('CORPER_RMC',i2s(xyzzyaadr1))
endif
call xyzzyaafe1('RMC_REP_LENGTH (Length of reptile)',i2s(xyzzyaabw1))
call xyzzyaafe1('RMC_MOVE_LENGTH',i2s(xyzzyaabx1))
call xyzzyaafe1('RMC_BOUNCE',l2s(xyzzyaacw1))
call xyzzyaafe1('RMC_MEAS_POS',l2s(xyzzyaacx1))
endif
if(.not.xyzzyaaen1)call wout()
end subroutine xyzzyaaff1
subroutine read_wave_function
implicit none
integer xyzzyaaaa13
integer,allocatable :: xyzzyaaab13(:),xyzzyaaac13(:)
if(.not.model_system)then
select case(trim(atom_basis_type))
case('plane-wave')
if(am_master)then
call wout('Reading plane-wave function and associated data')
call wout('===============================================')
endif
call readpwf(xyzzyaaef1)
case('gaussian')
if(am_master)then
call wout('Reading Gaussian wave function and associated data')
call wout('==================================================')
endif
call readgw(xyzzyaaef1)
case('slater-type')
if(am_master)then
call wout('Reading STO wave function and associated data')
call wout('=============================================')
endif
call stowfdet_read(xyzzyaaef1)
case('numerical')
if(am_master)then
call wout('Reading atomic wave function and associated data')
call wout('================================================')
endif
call readawf
xyzzyaaef1=0.d0
periodicity=0
case('blip')
if(xyzzyaaar1==-1)then
if(isperiodic)then
periodicity=3
else
periodicity=0
endif
else
periodicity=xyzzyaaar1
endif
if(am_master)then
call wout('Reading nonlocalized blip wave function and associated data&
&')
call wout('===========================================================&
&')
endif
call readbwf(xyzzyaaef1)
case('dimer')
if(am_master)then
call wout('Reading dimer wave function and associated data')
call wout('===============================================')
endif
call dwfdet_setup(xyzzyaaef1)
periodicity=0
case('none')
continue
case default
if(am_master)then
call wout('Setting up special wave function')
call wout('================================')
endif
call setup_special_wfn(xyzzyaaef1)
end select
else
nitot=0
xyzzyaaef1=0.d0
endif
if(.not.allocated(heg_orbtype))then
allocate(heg_orbtype(nspin,ndet),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'INIT_HEG','heg_orbtype')
heg_orbtype=0
endif
if(any(heg_orbtype/=0))then
if(xyzzyaadn1<0)then
allocate(xyzzyaaab13(nspin),xyzzyaaac13(nspin))
xyzzyaaab13(1:nspin)=heg_orbtype(1:nspin,1)
xyzzyaaac13(1:nspin)=heg_slatt(1:nspin,1)
deallocate(heg_orbtype,heg_slatt)
allocate(heg_orbtype(nspin,ndet),heg_slatt(nspin,ndet))
do xyzzyaaaa13=1,nspin
heg_orbtype(xyzzyaaaa13,1:ndet)=xyzzyaaab13(xyzzyaaaa13)
heg_slatt(xyzzyaaaa13,1:ndet)=xyzzyaaac13(xyzzyaaaa13)
enddo
deallocate(xyzzyaaab13,xyzzyaaac13)
elseif(xyzzyaadn1/=ndet)then
call errstop_master('READ_WAVE_FUNCTION','Number of explicit MD terms &
&in FREE_PARTICLES block must match number given in MDET block in corr&
&elation.data.')
endif
endif
call xyzzyaaez1
end subroutine read_wave_function
subroutine check_wave_function
implicit none
if(.not.model_system)then
if(periodicity==1.and.(scell_matrix(2,2)/=1.or.scell_matrix(3,3)/=1.or&
&.any(scell_matrix(1,2:3)/=0).or.scell_matrix(2,3)/=0.or.any(scell_mat&
&rix(2:3,1)/=0).or.scell_matrix(3,2)/=0)) call errstop_master('CHECK_W&
&AVE_FUNCTION','Polymer geometry selected. NPCELL(2) and NPCELL(3) mus&
&t = 1.')
if(periodicity==2.and.(scell_matrix(3,3)/=1.or.any(scell_matrix(3,1:2)&
&/=0).or.any(scell_matrix(1:2,3)/=0)))call errstop_master('CHECK_WAVE_&
&FUNCTION','Slab geometry selected but NPCELL(3) not equal to 1.')
endif
if(periodicity==1.and.interaction_mpc_present)call errstop_master('CHE&
&CK_WAVE_FUNCTION','MPC interaction for 1D periodic system not yet imp&
&lemented.')
if(.not.electron_system.and.interaction_mpc_present)call errstop_maste&
&r('CHECK_WAVE_FUNCTION','MPC interaction for particles other than ele&
&ctrons not implemented.')
end subroutine check_wave_function
subroutine xyzzyaafg1
implicit none
integer xyzzyaaaa15
logical xyzzyaaab15
xyzzyaaab15=.false.
if(.not.model_system.or.trim(atom_basis_type)=='plane-wave'.or.trim(at&
&om_basis_type)=='gaussian'.or.trim(atom_basis_type)=='slater-type')th&
&en
call init_geometry
if(init_by_ion.and.nitot/=xyzzyaadm1)call errstop_master('CALCULATE_GE&
&OMETRY','Size of EDIST_BY_ION input block not equal to number of ions&
&.')
if(init_by_iontype)then
allocate(neion(nitot,nspin),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'CALCULATE_GEOMETRY','neion <2>')
do xyzzyaaaa15=1,nitot
neion(xyzzyaaaa15,1:nspin)=xyzzyaadx1(iontype(xyzzyaaaa15),1:nspin)
enddo
if(sum(neion)/=(xyzzyaaaf1+xyzzyaaag1+xyzzyaaah1+xyzzyaaai1))then
if(am_master)call wout('Number of electrons in EDIST_BY_IONTYPE block &
&not equal to NEU+NED+NHU+NHD.')
xyzzyaaab15=.true.
endif
if(sum(neion(:,1))/=xyzzyaaaf1)then
if(am_master)call wout('No. of spin-up electrons in EDIST_BY_IONTYPE b&
&lock not equal to NEU.')
xyzzyaaab15=.true.
endif
if(sum(neion(:,2))/=xyzzyaaag1)then
if(am_master)call wout('No. of spin-down electrons in EDIST_BY_IONTYPE&
& block not equal to NED.')
xyzzyaaab15=.true.
endif
if(xyzzyaaah1>0)then
if(sum(neion(:,3))/=xyzzyaaah1)then
if(am_master)call wout('No. of spin-up holes in EDIST_BY_IONTYPE block&
& not equal to NHU.')
xyzzyaaab15=.true.
endif
endif
if(xyzzyaaai1>0)then
if(sum(neion(:,4))/=xyzzyaaai1)then
if(am_master)call wout('No. of spin-down holes in EDIST_BY_IONTYPE blo&
&ck not equal to NHD.')
xyzzyaaab15=.true.
endif
endif
if(xyzzyaaab15)call errstop_master('CALCULATE_GEOMETRY','Quitting.')
deallocate(xyzzyaadx1)
endif
endif
call netot_nitot_products
end subroutine xyzzyaafg1
subroutine xyzzyaafh1(aborted)
implicit none
logical,intent(out) :: aborted
logical xyzzyaaaa16,xyzzyaaab16
aborted=.false.
select case(trim(opt_method))
case('varmin','madmin','varmin_linjas')
if(ranprint>0)then
if(am_master)then
call open_units(ranlog_unit,xyzzyaadl1)
else
ranlog_unit=0
endif
endif
end select
call timer('SETUP',.false.)
select case(trim(opt_method))
case('varmin','madmin')
if(neighprint==0.and.trim(atom_basis_type)=='gaussian'.and.cusp_correc&
&tion)neighprint=1
if(trim(opt_method)=='madmin')then
vm_madmin=.true.
xyzzyaaaa16=vm_reweight
xyzzyaaab16=vm_filter
vm_reweight=.false.
vm_filter=.false.
endif
call varmin_main(aborted)
if(vm_madmin)then
vm_madmin=.false.
vm_reweight=xyzzyaaaa16
vm_filter=xyzzyaaab16
endif
case('varmin_linjas')
call varmin_linjas_main
case('emin')
call emin_main
case default
call errstop_master('RUN_OPT','OPT_METHOD problem (bug).')
end select
call timer('SETUP',.true.)
end subroutine xyzzyaafh1
subroutine xyzzyaafi1
implicit none
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17
real(dp) xyzzyaaad17,xyzzyaaae17
logical xyzzyaaaf17,xyzzyaaag17
if(opt_cycle>2)then
xyzzyaaaa17=opt_cycle
xyzzyaaab17=xyzzyaaaa17-1
xyzzyaaac17=xyzzyaaaa17-2
if(xyzzyaaeh1(xyzzyaaac17)>0.d0)then
xyzzyaaad17=hypothenuse(xyzzyaaeh1(xyzzyaaac17),xyzzyaaeh1(xyzzyaaab17&
&))
xyzzyaaae17=xyzzyaaeg1(xyzzyaaac17)-xyzzyaaeg1(xyzzyaaab17)
xyzzyaaaf17=xyzzyaaae17+3.d0*xyzzyaaad17<0.d0
xyzzyaaad17=hypothenuse(xyzzyaaeh1(xyzzyaaab17),xyzzyaaeh1(xyzzyaaaa17&
&))
xyzzyaaae17=xyzzyaaeg1(xyzzyaaab17)-xyzzyaaeg1(xyzzyaaaa17)
xyzzyaaag17=xyzzyaaae17+3.d0*xyzzyaaad17<0.d0
if(xyzzyaaaf17.and.xyzzyaaag17)then
if(opt_strict)call errstop_master('CHECK_OPT_CONVERGENCE','Optimized e&
&nergies appear to increase for two consecutive cycles. [You can turn &
&this error into a warning by setting OPT_STRICT=F in input.]')
call errwarn('CHECK_OPT_CONVERGENCE','Optimized energies appear to inc&
&rease for at least two consecutive cycles. [You can turn this warning&
& into an error by setting OPT_STRICT=T in input.]')
endif
endif
endif
end subroutine xyzzyaafi1
subroutine xyzzyaafj1(aborted,e,de)
implicit none
real(dp),intent(out),optional :: e,de
logical,intent(out) :: aborted
real(dp) xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18
logical xyzzyaaad18
if(isvmc.and.makemovie.and.my_node==movienode)then
call open_units(io_mov,xyzzyaadl1)
if(xyzzyaadl1/=0)call errstop_master('MONTE CARLO','Unable to find fre&
&e i/o unit <movie.out>.')
open(io_mov,file='movie.out',status='unknown',iostat=xyzzyaadl1)
if(xyzzyaadl1/=0)call errstop_master('MONTE CARLO','Problem opening mo&
&vie.out.')
endif
if(isperiodic.and.model_system)then
xyzzyaaac18=inv_netot
else
xyzzyaaac18=1.d0
endif
call check_hist_header(interaction,netot,nbasis,npcells,atom_basis_typ&
&e,isperiodic,xyzzyaaef1*xyzzyaaac18,'VMC',have_ppots,have_veep,relati&
&vistic,iaccum,newrun,eval_dipole_moment,eval_contact_den,writeout_vmc&
&_hist,mc_twist_av)
xyzzyaaad18=use_tmove
use_tmove=.false.
call timer('SETUP',.false.)
call vmc_main(xyzzyaaei1,xyzzyaaej1,xyzzyaaec1,xyzzyaaaa1,xyzzyaaad1,x&
&yzzyaaab1,xyzzyaaac1,xyzzyaaae1,xyzzyaaaj1,writeout_vmc_hist,xyzzyaac&
&v1,xyzzyaaaw1,xyzzyaaaa18,xyzzyaaab18,aborted,xyzzyaaek1,xyzzyaaeu1)
call timer('SETUP',.true.)
if(isvmc.and.makemovie.and.my_node==movienode)then
close(io_mov)
open_unit(io_mov)=.false.
endif
use_tmove=xyzzyaaad18
if(present(e))e=xyzzyaaaa18
if(present(de))de=xyzzyaaab18
end subroutine xyzzyaafj1
subroutine xyzzyaafk1(aborted)
implicit none
logical,intent(out) :: aborted
integer xyzzyaaaa19,xyzzyaaab19
real(dp) xyzzyaaac19,xyzzyaaad19,xyzzyaaae19
if(iaccum)then
xyzzyaaaa19=xyzzyaaal1
xyzzyaaab19=xyzzyaaan1
else
xyzzyaaaa19=xyzzyaaak1
xyzzyaaab19=xyzzyaaam1
endif
dmc_twist_av=.false.
if(isperiodic.and.model_system)then
xyzzyaaae19=inv_netot
else
xyzzyaaae19=1.d0
endif
call check_hist_header(interaction,netot,nbasis,npcells,atom_basis_typ&
&e,isperiodic,xyzzyaaef1*xyzzyaaae19,'DMC',have_ppots,have_veep,relati&
&vistic,iaccum,newrun,eval_dipole_moment,eval_contact_den,writeout_dmc&
&_hist,mc_twist_av)
if(mc_twist_av)then
call eval_finite_hf_energies(xyzzyaaac19,xyzzyaaad19)
else
xyzzyaaac19=0.d0
xyzzyaaad19=0.d0
endif
call timer('SETUP',.false.)
call dmc_main(xyzzyaace1,xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad&
&19,aborted)
call timer('SETUP',.true.)
end subroutine xyzzyaafk1
subroutine xyzzyaafl1(aborted)
implicit none
logical,intent(out) :: aborted
integer xyzzyaaaa20,xyzzyaaab20
real(dp) xyzzyaaac20
if(iaccum)then
xyzzyaaaa20=xyzzyaadu1
xyzzyaaab20=xyzzyaadv1
else
xyzzyaaaa20=xyzzyaads1
xyzzyaaab20=xyzzyaadt1
endif
if(isperiodic.and.model_system)then
xyzzyaaac20=inv_netot
else
xyzzyaaac20=1.d0
endif
call timer('SETUP',.false.)
call rmc_main(iaccum,xyzzyaacw1,xyzzyaack1,xyzzyaaaa20,xyzzyaaab20,xyz&
&zyaabw1,xyzzyaabx1,xyzzyaadr1,xyzzyaacx1,aborted)
call timer('SETUP',.true.)
end subroutine xyzzyaafl1
subroutine xyzzyaafm1(twist_first,aborted)
implicit none
logical,intent(in) :: twist_first
logical,intent(out) :: aborted
integer xyzzyaaaa21,xyzzyaaab21
real(dp) xyzzyaaac21,xyzzyaaad21,xyzzyaaae21
logical xyzzyaaaf21
call timer('SETUP',.false.)
if(.not.iaccum)call errstop_master('RUN_DMC_TWIST_AV','This subroutine&
& should only be called to perform statistics accumulation.')
if(twist_first.and.trim(random_seed_kw)/='timer'.and.trim(random_seed_&
&kw)/='timer_reset')call errstop_master('RUN_DMC_TWIST_AV','Please set&
& RANDOM_SEED to "timer" to avoid duplication of twists.')
if(isperiodic.and.model_system)then
xyzzyaaae21=inv_netot
else
xyzzyaaae21=1.d0
endif
call check_hist_header(interaction,netot,nbasis,npcells,atom_basis_typ&
&e,isperiodic,xyzzyaaef1*xyzzyaaae21,'DMC',have_ppots,have_veep,relati&
&vistic,iaccum,newrun,eval_dipole_moment,eval_contact_den,writeout_dmc&
&_hist,mc_twist_av)
if(twist_first)then
xyzzyaaab21=1
else
dmc_twist_av=.false.
call eval_finite_hf_energies(xyzzyaaac21,xyzzyaaad21)
if(am_master)then
call wout('===========================================================&
&=========')
call wout('NEW K-VECTOR OFFSET: ')
call wout('   ',k_offset,rfmt='(es21.13)')
call wout('Twist number: 1')
call wout('===========================================================&
&=========')
endif
call dmc_main(xyzzyaace1,xyzzyaaal1,xyzzyaaan1,xyzzyaaac21,xyzzyaaad21&
&,aborted)
if(aborted)return
xyzzyaaab21=2
endif
dmc_twist_av=.true.
do xyzzyaaaa21=xyzzyaaab21,xyzzyaaav1
if(am_master)then
call mc_twist_offset
call eval_finite_hf_energies(xyzzyaaac21,xyzzyaaad21)
call wout('===========================================================&
&=========')
call wout('NEW K-VECTOR OFFSET: ')
call wout('   ',k_offset,rfmt='(es21.13)')
call wout('Twist number: '//trim(i2s(xyzzyaaaa21)))
call wout('===========================================================&
&=========')
call wout()
endif
call mpi_bcast(k_offset,3,mpi_double_precision,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'Broadcasting k_offset in run_dmc_twist_av.')
if(.not.am_master)call mc_twist_offset(k_offset)
if(am_master)then
call wout('Post-twist-change equilibration...')
call wout()
endif
if(xyzzyaaaa21>1)call shift_config_files
xyzzyaael1=expvals
expvals=.false.
xyzzyaaaf21=writeout_dmc_hist
writeout_dmc_hist=.false.
iaccum=.false.
call dmc_main(xyzzyaace1,xyzzyaaat1,xyzzyaaau1,xyzzyaaac21,xyzzyaaad21&
&,aborted)
if(aborted)return
if(am_master)then
call wout('More statistics accumulation...')
call wout()
endif
call shift_config_files
expvals=xyzzyaael1
writeout_dmc_hist=xyzzyaaaf21
iaccum=.true.
call dmc_main(xyzzyaace1,xyzzyaaal1,xyzzyaaan1,xyzzyaaac21,xyzzyaaad21&
&,aborted)
if(aborted)return
enddo
call timer('SETUP',.true.)
end subroutine xyzzyaafm1
subroutine xyzzyaafn1
implicit none
select case(trim(atom_basis_type))
case('plane-wave')
call pwfdet_setup(xyzzyaacp1)
case('gaussian')
call gwfdet_setup
case('slater-type')
call stowfdet_setup
case('numerical')
call awfdet_setup
case('blip')
call bwfdet_setup
end select
if(any(heg_orbtype/=0))call init_free_orbs
if(any(heg_orbtype>100))call init_expot_wfdet
if(sparse.and..not.localized_orbitals)call errstop_master('ORBITAL_SET&
&UP','The SPARSE parameter is T but there are no localized orbitals.')
end subroutine xyzzyaafn1
subroutine xyzzyaafo1
implicit none
if(am_master.and..not.model_system)then
select case(trim(atom_basis_type))
case('plane-wave')
call wout('Geometry derived from information in pwfn.data')
call wout('==============================================')
case('gaussian')
call wout('Geometry derived from information in gwfn.data')
call wout('==============================================')
case('slater-type')
call wout('Geometry derived from information in stowfn.data')
call wout('==============================================')
case('numerical')
call wout('Geometry derived from information in awfn.data')
call wout('==============================================')
case('blip')
call wout('Geometry derived from information in bwfn.data')
call wout('==============================================')
case('dimer')
call wout("Geometry derived from information in dwfn.data")
call wout('==============================================')
end select
call wout()
call print_geometry()
endif
if(am_slave)then
if(neighprint==0.and.trim(atom_basis_type)=='gaussian'.and.cusp_correc&
&tion)neighprint=1
endif
if(neighprint>0)then
if(am_slave)allocate(nearest_ion_d(nitype),nearest_ion_r(nbasis))
call mpi_bcast(nearest_ion_d,nitype,mpi_double_precision,0,mpi_comm_wo&
&rld,ierror)
call checkmpi(ierror,'Broadcasting nearest_ion_d in geometry_printout.&
&')
call mpi_bcast(nearest_ion_r,nbasis,mpi_double_precision,0,mpi_comm_wo&
&rld,ierror)
call checkmpi(ierror,'Broadcasting nearest_ion_r in geometry_printout.&
&')
endif
end subroutine xyzzyaafo1
subroutine xyzzyaafp1
implicit none
call init_interactions
allocate(ion_fields(3,nitot),electron_fields(3,nitot),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'SETUP_INTERACTIONS','electric field array&
&s')
ion_fields=0.d0
electron_fields=0.d0
if(homogeneous_system)call calc_hf_energies
if(finite_size_corr.and.am_master)call eval_fscorr_consts(finite_size_&
&const_c)
end subroutine xyzzyaafp1
subroutine xyzzyaafq1
implicit none
if(am_master)then
call wout('Expectation values')
call wout('==================')
endif
if(finite_size_corr.and..not.structure_factor)then
if(am_master)then
call wout('Finite size correction flagged but STRUCTURE_FACTOR = F in &
&input.')
call wout('Activating structure factor accumulation.')
call wout()
endif
structure_factor=.true.
endif
finite_density=.false.
mol_density=.false.
mol_spin_density=.false.
spin_density_mat=.false.
if(density.or.spin_density)then
if(.not.isperiodic)then
if(nitot==1)then
finite_density=.true.
else
mol_density=.true.
if(spin_density)mol_spin_density=.true.
endif
density=.false.
spin_density=.false.
elseif(noncoll_spin)then
spin_density_mat=.true.
density=.false.
spin_density=.false.
endif
endif
expvals=density.or.spin_density.or.pair_corr.or.pair_corr_sph.or.loc_t&
&ensor.or.structure_factor.or.structure_factor_sph.or.onep_density_mat&
&.or.twop_density_mat.or.cond_fraction.or.mom_den.or.spin_density_mat.&
&or.finite_density.or.population.or.mol_density.or.twop_dm_mom.or.cond&
&_frac_mom
xyzzyaael1=.false.
if(expvals.and..not.(isvmc.or.isdmc.or.isvmc_dmc.or.isdmc_dmc))then
expvals=.false.
if(am_master)then
call wout('Expectation value accumulation disabled since not VMC or DM&
&C run.')
call wout()
endif
endif
if(expvals.or.(periodicity/=0.and.interaction_mpc_present.and.qmc_dens&
&ity_mpc).or.(int_sf.and.(isvmc.or.isdmc.or.isvmc_dmc.or.isdmc_dmc)))t&
&hen
if(am_master)call read_expval(.false.)
call setup_expval(.false.)
else
if(am_master)then
if(eval_dipole_moment)then
call wout('Requested expectation values:')
call wout('- dipole moment')
call wout()
call wout('Data accumulated in hist file; analyse with reblock utility&
&.')
call wout()
else
call wout('None requested.')
call wout()
endif
endif
endif
particle_is_fixed=.false.
end subroutine xyzzyaafq1
subroutine xyzzyaafr1
implicit none
integer,allocatable :: xyzzyaaaa26(:)
real(dp) xyzzyaaab26(3)
real(dp),allocatable :: xyzzyaaac26(:,:),xyzzyaaad26(:,:,:),xyzzyaaae2&
&6(:,:,:,:),xyzzyaaaf26(:,:,:),xyzzyaaag26(:,:,:,:)
if(.not.xyzzyaaco1)return
allocate(xyzzyaaac26(3,netot),xyzzyaaaa26(netot),xyzzyaaad26(nemax,rea&
&l1_complex2,ndet),xyzzyaaae26(3,nemax,real1_complex2,ndet),xyzzyaaaf2&
&6(nemax,real1_complex2,ndet),xyzzyaaag26(6,nemax,real1_complex2,ndet)&
&,stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'CHECK_ORBITAL_DERIVATIVES','1')
if(.not.allocated(neion))then
allocate(neion(nitot,nspin),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'CHECK_ORBITAL_DERIVATIVES','2')
endif
call points(xyzzyaaac26,xyzzyaaaa26,xyzzyaaek1,xyzzyaaeu1)
xyzzyaaab26(1:3)=xyzzyaaac26(1:3,1)
if(.not.use_backflow)then
call orbital_check(xyzzyaaab26,xyzzyaaac26,xyzzyaaad26,xyzzyaaae26,xyz&
&zyaaaf26)
else
call orbital_check(xyzzyaaab26,xyzzyaaac26,xyzzyaaad26,xyzzyaaae26,xyz&
&zyaaaf26,xyzzyaaag26)
endif
deallocate(xyzzyaaac26,xyzzyaaaa26,xyzzyaaad26,xyzzyaaae26,xyzzyaaaf26&
&,xyzzyaaag26)
end subroutine xyzzyaafr1
subroutine xyzzyaafs1
implicit none
character(80) tmpr
allocate(nlrule1_ion(nitype),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'READ_PSEUDOPOTENTIALS','nlrule')
if(.not.model_system.or.trim(atom_basis_type)=='plane-wave'.or.trim(at&
&om_basis_type)=='gaussian')then
nlrule1_ion(:)=xyzzyaaas1
call read_ppots
if(have_veep)then
if(isperiodic)calc_field=.true.
if(periodicity==1)call errstop_master('READ_PSEUDOPOTENTIALS','Core po&
&larization potentials not implemented for 1D systems. Ask Mike')
endif
else
have_veep=.false.
endif
if(.not.have_ae)cusp_correction=.false.
if((.not.have_ae.or.model_system.or.dimensionality/=3).and.(isvmc_dmc.&
&or.isdmc.or.isdmc_dmc.or.isvmc_dmc_equil))then
if(nucleus_gf_mods.and.xyzzyaact1)then
tmpr=r2s(alimit,'(f5.3)')
call errwarn('READ_PSEUDOPOTENTIALS','NUCLEUS_GF_MODS set in input, bu&
&t there are no bare nuclei. Cannot therefore use DMC Green''s functio&
&n modifications. ALIMIT='//trim(tmpr)//'.')
endif
nucleus_gf_mods=.false.
endif
end subroutine xyzzyaafs1
subroutine read_jastrow_function
implicit none
integer xyzzyaaaa28
logical exists,xyzzyaaab28
if(use_jastrow)then
if(.not.xyzzyaaep1)then
call query_casl_item('parameters.casl:JASTROW',exists=exists)
use_gjastrow=.not.jastrow_in_corr.or.exists
endif
xyzzyaadi1=0
xyzzyaadh1=0
if(isvmc_opt.or.isopt_vmc.or.isopt)then
do xyzzyaaaa28=1,opt_cycles
if(xyzzyaada1(xyzzyaaaa28))exit
enddo
if(xyzzyaaaa28<=opt_cycles)then
if(trim(xyzzyaacz1(xyzzyaaaa28))/='emin')then
xyzzyaadi1=xyzzyaaaa28
else
xyzzyaadh1=xyzzyaaaa28
prefer_short_cusp=.true.
endif
endif
endif
call read_jastrow(xyzzyaaab28)
if(xyzzyaaab28)then
if(am_master)call wout('Initial Jastrow set is empty.')
if((xyzzyaadi1>0.or.xyzzyaadh1>1).and..not.use_altsamp)then
use_jastrow=.false.
if(am_master)call wout('Will not use a Jastrow factor until needed.')
elseif(isvmc.or.isvmc_dmc.or.isdmc.or.isdmc_dmc.or.isvmc_dmc_equil)the&
&n
use_jastrow=.false.
if(am_master)call wout('Will run entire calculation without a Jastrow.&
&')
else
if(am_master)call wout('Jastrow factor will be used throughout.')
endif
if(am_master)call wout()
endif
if(.not.use_gjastrow.and.(isopt.or.isvmc_opt.or.isopt_vmc).and.any(ind&
&ex(xyzzyaacz1,'varmin_linjas')==1))call check_varmin_linjas_jastrow
call init_wfn
endif
end subroutine read_jastrow_function
subroutine xyzzyaaft1
implicit none
integer xyzzyaaaa29
logical xyzzyaaab29
if(use_backflow)then
call init_backflow(xyzzyaaab29)
if(xyzzyaaab29)then
use_backflow=.false.
call init_wfn
xyzzyaadj1=0
if(isvmc_opt.or.isopt_vmc.or.isopt)then
do xyzzyaaaa29=1,opt_cycles
if(xyzzyaadb1(xyzzyaaaa29))exit
enddo
if(xyzzyaaaa29<=opt_cycles)xyzzyaadj1=xyzzyaaaa29
endif
if(am_master)then
call wout('Initial backflow set is empty.')
if(xyzzyaadj1>0)then
call wout('Will not use backflow until needed.')
else
call wout('Will run entire calculation without backflow.')
endif
call wout()
endif
endif
endif
end subroutine xyzzyaaft1
subroutine xyzzyaafu1
implicit none
integer xyzzyaaaa30,xyzzyaaab30
real(dp) xyzzyaaac30
real(dp),allocatable :: xyzzyaaad30(:)
logical xyzzyaaae30
character(80) tmpr
if(am_master)then
xyzzyaaae30=.false.
do xyzzyaaaa30=1,nitot-1
do xyzzyaaab30=xyzzyaaaa30+1,nitot
if(dot_product(rion(1:3,xyzzyaaaa30)-rion(1:3,xyzzyaaab30),rion(1:3,xy&
&zzyaaaa30)-rion(1:3,xyzzyaaab30))<1.d-10)then
call errwarn_silent('CHECK_NN','Ions '//trim(i2s(xyzzyaaaa30))//' and &
&'//trim(i2s(xyzzyaaab30))//' are on top of each other.  Using EIONION&
& value in xwfn.data.')
tmpr=r2s(xyzzyaaef1,'(f21.12)')
call wout('EIONION from file  : '//trim(tmpr))
call wout()
xyzzyaaae30=.true.
exit
endif
enddo
enddo
if(nitot>0.and..not.xyzzyaaae30.and..not.ignore_ionic_interactions)the&
&n
if(trim(atom_basis_type)/='numerical'.and..not.(trim(atom_basis_type)=&
&='gaussian'.and.nitot<2).and..not.(trim(atom_basis_type)=='slater-typ&
&e'.and.nitot<2).and.trim(atom_basis_type)/='non_int_he'.and.trim(atom&
&_basis_type)/='h2'.and.trim(atom_basis_type)/='h3plus')then
if(isperiodic)then
if(have_ppots)then
call wout('Ionic repulsion energy (au/primitive cell)')
call wout('==========================================')
else
call wout('Nuclear repulsion energy (au/primitive cell)')
call wout('============================================')
endif
else
if(have_ppots)then
call wout('Ionic repulsion energy (au)')
call wout('===========================')
else
call wout('Nuclear repulsion energy (au)')
call wout('=============================')
endif
endif
endif
allocate(xyzzyaaad30(nitot),stat=xyzzyaadk1)
call check_alloc(xyzzyaadk1,'CHECK_NN','ztemp')
xyzzyaaad30(1:nitot)=zion(iontype(1:nitot))
call nuclear_repulsion_energy(nitot,xyzzyaaad30,xyzzyaaac30)
deallocate(xyzzyaaad30)
if(trim(atom_basis_type)/='numerical'.and..not.(trim(atom_basis_type)=&
&='gaussian'.and.nitot<2).and.trim(atom_basis_type)/='non_int_he'.and.&
&xyzzyaaef1/=0.d0)then
tmpr=r2s(xyzzyaaef1,'(f21.12)')
call wout('EIONION from file  : '//trim(tmpr))
tmpr=r2s(xyzzyaaac30,'(f21.12)')
call wout('Calculated EIONION : '//trim(tmpr))
call wout()
if(abs(xyzzyaaef1-xyzzyaaac30)>1.d-5)then
tmpr=r2s(abs(xyzzyaaef1-xyzzyaaac30),'(f21.12)')
call wout('Difference         : '//trim(tmpr))
call wout()
select case(trim(atom_basis_type))
case('plane-wave')
call errwarn_silent('CHECK_NN','calculated nuclear repulsion energy di&
&ffers from value in pwfn.data. If the difference between the numbers &
&above is large, you should investigate this problem.')
case('gaussian')
call errwarn_silent('CHECK_NN','calculated nuclear repulsion energy di&
&ffers from value in gwfn.data. If the difference between the numbers &
&above is large, you should investigate this problem.')
case('slater-type')
call errwarn_silent('CHECK_NN','calculated nuclear repulsion energy di&
&ffers from value in stowfn.data. If the difference between the number&
&s above is large, you should investigate this problem.')
case('blip')
if(periodicity<3)then
call wordwrap('Calculated nuclear repulsion energy differs from value &
&in bwfn.data. This is probably because the orbital-generating code us&
&ed 3D periodicity to simulate the system, while CASINO is using ' //t&
&rim(i2s(periodicity))//'D periodicity.   However, if this is not the &
&case and the difference between the numbers above is large, you shoul&
&d investigate this problem.')
call wout()
else
call errwarn_silent('CHECK_NN','calculated nuclear repulsion energy di&
&ffers from value in bwfn.data. If the difference between the numbers &
&above is large, you should investigate this problem.')
endif
case('dimer')
call errwarn('CHECK_NN','calculated nuclear repulsion energy differs f&
&rom value in dwfn.data. If the difference between the numbers above i&
&s large, you should investigate this problem.')
case default
call errstop_master('CHECK_NN','Confused.')
end select
if(periodicity==1)then
call wordwrap('Note this is known to be an issue for 1D polymers with &
&long unit cells.')
call wout()
endif
if(have_ppots.and..not.(trim(atom_basis_type)=="blip" .and.periodicity&
&<3))then
call wordwrap('Please check that the xx_pp.data pseudopotential file b&
&oth exists and contains a valid pseudopotential (i.e. it''s not a lin&
&k to a now gzipped file, for example). If it does, then either CASINO&
& or the wave function generating program is calculating Ewald interac&
&tions inaccurately.')
call wout()
if(.not.xyzzyaacq1)then
call wordwrap("If you want CASINO to ignore this problem, set 'EWALD_C&
&HECK : T' in the input file.")
call errstop_master('CHECK_NN','Quitting.')
else
call wordwrap('CASINO is ignoring this problem because the EWALD_CHECK&
& flag has been set in the input file.')
call wout()
endif
endif
else
call wout('Calculated and input nuclear repulsion energies agree.')
call wout()
endif
endif
xyzzyaaef1=xyzzyaaac30
if(model_system)then
tmpr=r2s(xyzzyaaef1,'(f21.12)')
call wout('Calculated EIONION : '//trim(tmpr))
call wout()
endif
endif
endif
if(trim(atom_basis_type)/='numerical'.and.trim(atom_basis_type)/='non_&
&int_he'.and.nitot>0)then
call mpi_bcast(xyzzyaaef1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting eionion in monte_carlo.f90')
endif
constant_energy=constant_energy+xyzzyaaef1*npcells
end subroutine xyzzyaafu1
subroutine xyzzyaafv1
implicit none
real(dp) xyzzyaaaa31
character(80) tmpr
if(have_ppots)then
if(have_veep)call timer('CORE_POL_II',.true.)
if(isperiodic)then
call compute_ii_fields
xyzzyaaaa31=0.d0
else
call compute_ecppii(xyzzyaaaa31)
if(am_master.and.have_veep)then
call wout('Ion-ion core polarization contribution to total energy')
call wout('======================================================')
tmpr=r2s(xyzzyaaaa31,'(f21.12)')
call wout('ECPPII             : '//trim(tmpr))
call wout()
endif
endif
if(have_veep)call timer('CORE_POL_II',.false.)
else
xyzzyaaaa31=0.d0
endif
constant_energy=constant_energy+xyzzyaaaa31*npcells
end subroutine xyzzyaafv1
subroutine xyzzyaafw1
implicit none
integer xyzzyaaaa32
logical xyzzyaaab32,xyzzyaaac32
character(160) con_in_nofixed,con_in_fixed,con_out_nofixed,con_out_fix&
&ed
con_in_fixed=trim(con_in)//'_fixed'
con_in_nofixed=trim(con_in)//'_nofixed'
con_out_fixed=trim(con_out)//'_fixed'
con_out_nofixed=trim(con_out)//'_nofixed'
if(am_master)then
call wout()
call wout('ACCUMULATION OF PAIR CORRELATION FUNCTION REQUESTED')
call wout()
endif
if(pair_corr.and.pair_corr_sph)call errstop_master('RUN_QMC_PCF','Rout&
&ine should not be called with pair_corr and pair_corr_sph active. Bug&
&.')
turn_on_pair_corr=.false.
turn_on_pair_corr_sph=.false.
if(pair_corr)then
pair_corr=.false.
turn_on_pair_corr=.true.
endif
if(pair_corr_sph)then
pair_corr_sph=.false.
turn_on_pair_corr_sph=.true.
endif
if(isvmc)then
if(am_master)then
call wout('===========================================================&
&====')
call wout('PERFORMING FIRST OF TWO VMC CALCULATIONS (NO FIXED PARTICLE&
&) TO')
call wout('ACCUMULATE SPIN DENSITY BLOCK IN expval.data REQUIRED FOR P&
&CF.')
call wout('===========================================================&
&====')
call wout()
endif
if(.not.newrun)call copy_config_file(con_in_nofixed,con_in,move=.true.&
&)
call xyzzyaafj1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call copy_config_file(con_out,con_out_nofixed,move=.true.)
call dismantle_configs
if(am_master)then
call xyzzyaagc1
call wout('===========================================================&
&====')
call wout('PERFORMING SECOND OF TWO VMC CALCULATIONS (WITH FIXED PARTI&
&CLE)')
if(turn_on_pair_corr)then
call wout('TO ACCUMULATE RECIPROCAL-SPACE PCF BLOCK IN expval.data.')
else
call wout('TO ACCUMULATE SPHERICAL PCF BLOCK IN expval.data.')
endif
call wout('===========================================================&
&====')
call wout()
endif
if(turn_on_pair_corr)call switch_spin_density_to_pcf
if(turn_on_pair_corr_sph)call switch_spin_density_to_pcf_sph
if(.not.newrun)call copy_config_file(con_in_fixed,con_in,move=.true.)
call xyzzyaafj1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call copy_config_file(con_out,con_out_fixed,move=.true.)
call copy_config_file(con_out,con_in)
elseif(isdmc)then
inquire(file=trim(con_in_nofixed),exist=xyzzyaaab32)
inquire(file=trim(con_in_fixed),exist=xyzzyaaac32)
if(.not.(xyzzyaaab32.and.xyzzyaaac32))then
call errstop_master('RUN_QMC_PCF','Need both '//trim(con_in_nofixed)//&
&' and '//trim(con_in_fixed)//' files (normal/fixed particle) in order&
& to continue PCF accumulation run.')
endif
if(am_master)then
call wout('===========================================================&
&====')
call wout('PERFORMING FIRST OF TWO DMC CALCULATIONS (NO FIXED PARTICLE&
&) TO')
call wout('ACCUMULATE SPIN DENSITY BLOCK IN expval.data REQUIRED FOR P&
&CF.')
call wout('===========================================================&
&====')
call wout()
endif
expvals=.false.
call copy_config_file(con_in_nofixed,con_in,move=.true.)
call xyzzyaafk1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call copy_config_file(con_out,con_out_nofixed,move=.true.)
if(am_master)then
call xyzzyaagc1
call wout('===========================================================&
&====')
call wout('PERFORMING SECOND OF TWO DMC CALCULATIONS (WITH FIXED PARTI&
&CLE)')
if(turn_on_pair_corr)then
call wout('TO ACCUMULATE RECIPROCAL-SPACE PCF BLOCK IN expval.data.')
else
call wout('TO ACCUMULATE SPHERICAL PCF BLOCK IN expval.data.')
endif
call wout('(dmc.hist will not be written to during this step.)')
call wout('===========================================================&
&====')
call wout()
endif
expvals=.true.
if(turn_on_pair_corr)call switch_spin_density_to_pcf
if(turn_on_pair_corr_sph)call switch_spin_density_to_pcf_sph
call copy_config_file(con_in_fixed,con_in,move=.true.)
call xyzzyaafk1(xyzzyaaeo1)
if(xyzzyaaeo1)return
call copy_config_file(con_out,con_out_fixed,move=.true.)
elseif(isvmc_dmc)then
do xyzzyaaaa32=1,2
expvals=.false.
newrun=.true.
if(am_master)then
if(xyzzyaaaa32==1)then
call wout('===========================================================&
&=========')
call wout('PERFORMING FIRST OF TWO FULL DMC CALCULATIONS (NO FIXED PAR&
&TICLE) TO')
call wout('ACCUMULATE SPIN DENSITY BLOCK IN expval.data REQUIRED FOR P&
&CF.')
call wout('===========================================================&
&=========')
else
call xyzzyaagc1
call wout('===========================================================&
&=========')
call wout('PERFORMING SECOND OF TWO FULL DMC CALCULATIONS (WITH FIXED &
&PARTICLE)')
if(turn_on_pair_corr)then
call wout('TO ACCUMULATE RECIPROCAL-SPACE PCF BLOCK IN expval.data.')
else
call wout('TO ACCUMULATE SPHERICAL PCF BLOCK IN expval.data.')
endif
call wout('(dmc.hist will not be written to during this step.)')
call wout('===========================================================&
&=========')
endif
call wout()
call xyzzyaagd1('PERFORMING A VMC CONFIGURATION-GENERATION CALCULATION&
&.')
endif
call xyzzyaafj1(xyzzyaaeo1)
if(xyzzyaaeo1)return
iaccum=.false.
call xyzzyaagc1
call xyzzyaagd1('PERFORMING A DMC EQUILIBRATION CALCULATION.')
call shift_config_files
call xyzzyaafk1(xyzzyaaeo1)
if(xyzzyaaeo1)return
iaccum=.true.
expvals=.true.
newrun=.false.
call xyzzyaagc1
call xyzzyaagd1('PERFORMING A DMC STATISTICS-ACCUMULATION CALCULATION.&
&')
call shift_config_files
call xyzzyaafk1(xyzzyaaeo1)
if(xyzzyaaeo1)return
if(xyzzyaaaa32==1)then
if(turn_on_pair_corr)call switch_spin_density_to_pcf
if(turn_on_pair_corr_sph)call switch_spin_density_to_pcf_sph
if(xyzzyaael1)expvals=.false.
call copy_config_file(con_out,con_out_nofixed,move=.true.)
call copy_config_file(con_out,con_in)
else
call copy_config_file(con_out,con_out_fixed,move=.true.)
endif
enddo
elseif(isdmc_dmc)then
do xyzzyaaaa32=1,2
newrun=.true.
if(am_master)then
if(xyzzyaaaa32==1)then
call wout('===========================================================&
&=========')
call wout('PERFORMING FIRST OF TWO FULL DMC CALCULATIONS (NO FIXED PAR&
&TICLE) TO')
call wout('ACCUMULATE SPIN DENSITY BLOCK IN expval.data REQUIRED FOR P&
&CF.')
call wout('===========================================================&
&=========')
else
call xyzzyaagc1
call wout('===========================================================&
&=========')
call wout('PERFORMING SECOND OF TWO FULL DMC CALCULATIONS (WITH FIXED &
&PARTICLE)')
if(turn_on_pair_corr)then
call wout('TO ACCUMULATE RECIPROCAL-SPACE PCF BLOCK IN expval.data.')
else
call wout('TO ACCUMULATE SPHERICAL PCF BLOCK IN expval.data.')
endif
call wout('(dmc.hist will not be written to during this step.)')
call wout('===========================================================&
&=========')
endif
call wout()
call xyzzyaagd1('PERFORMING A DMC EQUILIBRATION CALCULATION.')
endif
expvals=.false.
iaccum=.false.
call xyzzyaafk1(xyzzyaaeo1)
if(xyzzyaaeo1)return
iaccum=.true.
newrun=.false.
expvals=.true.
call xyzzyaagc1
call xyzzyaagd1('PERFORMING A DMC STATISTICS-ACCUMULATION CALCULATION.&
&')
call shift_config_files
call xyzzyaafk1(xyzzyaaeo1)
if(xyzzyaaeo1)return
if(xyzzyaaaa32==1)then
if(turn_on_pair_corr)call switch_spin_density_to_pcf
if(turn_on_pair_corr_sph)call switch_spin_density_to_pcf_sph
if(xyzzyaael1)expvals=.false.
call copy_config_file(con_out,con_out_nofixed,move=.true.)
else
call copy_config_file(con_out,con_out_fixed,move=.true.)
endif
enddo
else
call errstop_master('RUN_QMC_PCF','Routine called with inappropriate r&
&un type.')
endif
end subroutine xyzzyaafw1
subroutine xyzzyaafx1(nlines)
implicit none
integer,intent(in) :: nlines
integer xyzzyaaaa33,xyzzyaaab33
logical xyzzyaaac33
real(dp) xyzzyaaad33
character(80) char_80,varname,dummy
man_int_params=0.d0
read(block_data(1),'(a)')char_80
int_name=trim(adjustl(char_80))
select case(trim(int_name))
case('square_well')
if(nlines/=3)call errstop_master('READ_MAN_INT','Wrong number of lines&
& in MAN_INT block in input file.')
do xyzzyaaaa33=2,nlines
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,dummy,xyzzya&
&aad33
if(xyzzyaaab33/=0)then
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,xyzzyaaad33
if(xyzzyaaab33/=0)call errstop_master('READ_MAN_INT','Syntax error rea&
&ding line '//trim(i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in&
& input.')
else
select case(trim(dummy))
case(':','=')
continue
case default
call errstop_master('READ_MAN_INT','Syntax error reading line '//trim(&
&i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in input (second tok&
&en not ":"/"=").')
end select
endif
select case(trim(varname))
case('width','Width','WIDTH')
man_int_params(1)=xyzzyaaad33
case('height','Height','HEIGHT')
man_int_params(2)=xyzzyaaad33
case default
call errstop_master('READ_MAN_INT','Unrecognized parameter in MANUAL_I&
&NTERACTION block. Allowed SQUARE_WELL parameters are ''width'' and ''&
&height''.')
end select
enddo
case('poschl_teller')
if(nlines/=3)call errstop_master('READ_MAN_INT','Wrong number of lines&
& in MANUAL_INTERACTION block in input file.')
do xyzzyaaaa33=2,nlines
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,dummy,xyzzya&
&aad33
if(xyzzyaaab33/=0)then
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,xyzzyaaad33
if(xyzzyaaab33/=0)call errstop_master('READ_MAN_INT','Syntax error rea&
&ding line '//trim(i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in&
& input (second token not ":"/"=").')
else
select case(trim(dummy))
case(':','=')
continue
case default
call errstop_master('READ_MAN_INT','Syntax error reading line '//trim(&
&i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in input (second tok&
&en not ":"/"=").')
end select
endif
select case(trim(varname))
case('mu','Mu')
man_int_params(1)=xyzzyaaad33
case('v_0','V_0')
man_int_params(2)=xyzzyaaad33
case default
call errstop_master('READ_MAN_INT','Unrecognized parameter in MANUAL_I&
&NTERACTION block. Allowed POSCHL_TELLER parameters are ''mu'' and ''v&
&_0''.')
end select
enddo
case('hard_sphere')
man_int_op_spins=.false.
if(nlines<2.or.nlines>4)call errstop_master('READ_MAN_INT','Wrong numb&
&er of lines in MANUAL_INTERACTION block in input file.')
do xyzzyaaaa33=2,nlines
xyzzyaaac33=.false.
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,dummy,xyzzya&
&aad33
if(xyzzyaaab33/=0)then
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,xyzzyaaad33
if(xyzzyaaab33/=0)then
xyzzyaaac33=.true.
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname
if(xyzzyaaab33/=0)call errstop_master('READ_MAN_INT','Syntax error rea&
&ding line '//trim(i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in&
& input.')
endif
else
select case(trim(dummy))
case(':','=')
continue
case default
call errstop_master('READ_MAN_INT','Syntax error reading line '//trim(&
&i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in input (second tok&
&en not ":"/"=").')
end select
endif
varname=switch_case(varname,to_lower=.true.)
select case(trim(varname))
case('op_spins')
man_int_op_spins=.true.
case default
if(xyzzyaaac33)call errstop_master('READ_MAN_INT','Invalid/missing val&
&ue for variable "'//trim(varname)//'" in line '//trim(i2s(xyzzyaaaa33&
&))//' of MANUAL_INTERACTION block in input.')
select case(trim(varname))
case('d','hard_diam')
man_int_params(1)=0.5d0*xyzzyaaad33
case('r','hard_radius')
man_int_params(1)=xyzzyaaad33
case('lambda')
man_int_params(2)=xyzzyaaad33
case default
call errstop_master('READ_MAN_INT','Unrecognized parameter in MANUAL_I&
&NTERACTION block. The first HARD_SPHERE parameter is the radius ''R''&
& (you can also give the diameter ''D''). You can optionally give the &
&weight of an r^(-3) tail, ''Lambda''.')
end select
end select
enddo
if(.not.use_jastrow.or.use_gjastrow)call errstop_master('READ_MAN_INT'&
&,'Hard-sphere interactions can only be used in conjunction with the s&
&tandard CASINO two-body Jastrow factor.')
case('polynomial')
call xyzzyaafy1(nlines)
case('logarithmic','2D_int')
if(nlines/=2)call errstop_master('READ_MAN_INT','Wrong number of lines&
& in MAN_INT block in input file.')
read(block_data(2),*,iostat=xyzzyaaab33)varname,dummy,xyzzyaaad33
if(xyzzyaaab33/=0)then
read(block_data(2),*,iostat=xyzzyaaab33)varname,xyzzyaaad33
if(xyzzyaaab33/=0)call errstop_master('READ_MAN_INT','Syntax error rea&
&ding line 2 of MANUAL_INTERACTION block in input.')
else
select case(trim(dummy))
case(':','=')
continue
case default
call errstop_master('READ_MAN_INT','Syntax error reading line 2 of MAN&
&_INT block in input (second token not ":"/"=").')
end select
endif
select case(trim(varname))
case('rstar')
man_int_params(1)=xyzzyaaad33
case default
call errstop_master('READ_MAN_INT','Unrecognized parameter in MANUAL_I&
&NTERACTION block. Allowed LOGARITHMIC / 2D_INT parameter: "rstar".')
end select
case('dipole')
if(nlines/=2)call errstop_master('READ_MAN_INT','Wrong number of lines&
& in MANUAL_INTERACTION block in input file.')
do xyzzyaaaa33=2,nlines
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,dummy,xyzzya&
&aad33
if(xyzzyaaab33/=0)then
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,xyzzyaaad33
if(xyzzyaaab33/=0)call errstop_master('READ_MAN_INT','Syntax error rea&
&ding line '//trim(i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in&
& input.')
else
select case(trim(dummy))
case(':','=')
continue
case default
call errstop_master('READ_MAN_INT','Syntax error reading line '//trim(&
&i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in input.')
end select
endif
select case(trim(varname))
case('d^2')
man_int_params(1)=xyzzyaaad33
case default
call errstop_master('READ_MAN_INT','Unrecognized parameter in MANUAL_I&
&NTERACTION block.  Allowed dipole parameter is d^2.')
end select
enddo
case('pseudodipole')
call xyzzyaafz1(nlines)
case('2D_tilt_pseudodipole')
call xyzzyaaga1(nlines)
case('2D_tilted_dipole')
if(nlines/=3)call errstop_master('READ_MAN_INT','Wrong number of lines&
& in MANUAL_INTERACTION block in input file.')
do xyzzyaaaa33=2,nlines
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,dummy,xyzzya&
&aad33
if(xyzzyaaab33/=0)then
read(block_data(xyzzyaaaa33),*,iostat=xyzzyaaab33)varname,xyzzyaaad33
if(xyzzyaaab33/=0)call errstop_master('READ_MAN_INT','Syntax error rea&
&ding line '//trim(i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in&
& input.')
else
select case(trim(dummy))
case(':','=')
continue
case default
call errstop_master('READ_MAN_INT','Syntax error reading line '//trim(&
&i2s(xyzzyaaaa33))//' of MANUAL_INTERACTION block in input.')
end select
endif
select case(trim(varname))
case('d^2')
man_int_params(1)=xyzzyaaad33
case('theta')
man_int_params(2)=xyzzyaaad33
case default
call errstop_master('READ_MAN_INT','Unrecognized parameter in MANUAL_I&
&NTERACTION block.  Allowed dipole parameters are d^2 and theta.')
end select
enddo
case default
call errstop_master('READ_MAN_INT','Unknown manual interaction thingy &
&'//trim(int_name)//'.')
end select
end subroutine xyzzyaafx1
subroutine xyzzyaafy1(nlines)
implicit none
integer, intent(in) :: nlines
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34
logical xyzzyaaad34
character(80) varname,data_buffer
if(nlines<3)call errstop_master('READ_MAN_INT_POLYNOMIAL', 'Wrong numb&
&er of lines in manual_interaction block in input file. Read '//trim(i&
&2s(nlines))//' lines but expected at least 2.')
xyzzyaaad34=.false.
do xyzzyaaaa34=2,nlines
if(.not.xyzzyaagb1(block_data(xyzzyaaaa34),varname,data_buffer))call e&
&rrstop_master('READ_MAN_INT_POLYNOMIAL','Syntax error reading line '/&
&/trim(i2s(xyzzyaaaa34))//' of manual interactions block in input. Syn&
&tax should be "kw : value".')
if(trim(switch_case(varname))=='ORDER')then
xyzzyaaad34=.true.
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_order
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_POLYNOMIAL','Syntax&
& error reading value for "order" in manual_interactions block. Could &
&not cast input to integer.')
allocate(man_int_poly_coeffs(0:man_int_poly_order),stat=xyzzyaaab34)
call check_alloc(xyzzyaaab34,'READ_MAN_INT_POLYNOMIAL','coefficients')
man_int_poly_coeffs=0.d0
exit
endif
enddo
if(.not.xyzzyaaad34)call errstop_master('READ_MAN_INT_POLYNOMIAL','Fai&
&led to find keyword "order" in manual_interactions block.')
do xyzzyaaaa34=2,nlines
if(.not.xyzzyaagb1(block_data(xyzzyaaaa34),varname,data_buffer))call e&
&rrstop_master('READ_MAN_INT_POLYNOMIAL','Syntax error reading line '/&
&/trim(i2s(xyzzyaaaa34))//' of manual_interactions block in input. Syn&
&tax should be "kw : value".')
select case(trim(switch_case(varname)))
case('ORDER')
continue
case('CUTOFF')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_cutoff
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_POLYNOMIAL','Syntax&
& error reading value for "order" in manual_interactions block. Could &
&not cast input to float.')
case default
if(varname(1:2)=='c_'.or.varname(1:2)=='C_')then
read(varname(3:),*,iostat=xyzzyaadl1)xyzzyaaac34
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_POLYNOMIAL','Syntax&
& error reading line '//trim(i2s(xyzzyaaaa34))//' of manual interactio&
&ns block in input. Failed to coerce coefficient to integer.')
if(xyzzyaaac34<0.or.xyzzyaaac34>man_int_poly_order)call errstop_master&
&('READ_MAN_INT_POLYNOMIAL','Invalid coefficient specified in manual_i&
&nteractions block : '//trim(varname)//'. The coefficient must be 0 <=&
& c <= order.')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_coeffs(xyzzyaaac34)
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_POLYNOMIAL','Syntax&
& error reading value for "'//trim(varname)//'" in manual_interactions&
& block. Could not cast input to float.')
else
call errstop_master('READ_MAN_INT_POLYNOMIAL','Unrecognized keyword in&
& manual_interactions block: '//trim(varname)//'. Allowed POLYNOMIAL p&
&arameters are "order", "cutoff" and "coefficients".')
endif
end select
enddo
end subroutine xyzzyaafy1
subroutine xyzzyaafz1(nlines)
implicit none
integer, intent(in) :: nlines
integer xyzzyaaaa35,xyzzyaaab35,xyzzyaaac35
logical xyzzyaaad35
character(80) varname,data_buffer
if(nlines<4)call errstop_master('READ_MAN_INT_PSEUDODIPOLE', 'Wrong nu&
&mber of lines in manual_interaction block in input file. Read '//trim&
&(i2s(nlines))//' lines but expected at least 3.')
xyzzyaaad35=.false.
do xyzzyaaaa35=2,nlines
if(.not.xyzzyaagb1(block_data(xyzzyaaaa35),varname,data_buffer))call e&
&rrstop_master('READ_MAN_INT_PSEUDODIPOLE','Syntax error reading line &
&'//trim(i2s(xyzzyaaaa35))//' of manual interactions block in input. S&
&yntax should be "keyword : value".')
if(trim(switch_case(varname))=='ORDER')then
xyzzyaaad35=.true.
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_order
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading value for "order" in manual_interactions block. Coul&
&d not cast input to integer.')
allocate(man_int_poly_coeffs(0:man_int_poly_order),stat=xyzzyaaab35)
call check_alloc(xyzzyaaab35,'READ_MAN_INT_PSEUDODIPOLE','coefficients&
&')
man_int_poly_coeffs=0.d0
exit
endif
enddo
if(.not.xyzzyaaad35)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','F&
&ailed to find keyword "order" in manual_interactions block.')
do xyzzyaaaa35=2,nlines
if(.not.xyzzyaagb1(block_data(xyzzyaaaa35),varname,data_buffer))call e&
&rrstop_master('READ_MAN_INT_PSEUDODIPOLE','Syntax error reading line &
&'//trim(i2s(xyzzyaaaa35))//' of manual interactions block in input. S&
&yntax should be "keyword : value".')
select case(trim(switch_case(varname)))
case('ORDER')
continue
case('CUTOFF')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_cutoff
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading value for "cutoff" in manual_interactions block. Cou&
&ld not cast input to float.')
case('D^2')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_params(1)
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading value for "d^2" in manual_interactions block. Could &
&not cast input to float.')
case default
if(varname(1:2)=='c_'.or.varname(1:2)=='C_')then
read(varname(3:),*,iostat=xyzzyaadl1)xyzzyaaac35
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading line '//trim(i2s(xyzzyaaaa35))//' of manual interact&
&ions block in input. Failed to coerce coefficient to integer.')
if(xyzzyaaac35<0.or.xyzzyaaac35>man_int_poly_order)call errstop_master&
&('READ_MAN_INT_PSEUDODIPOLE','Invalid coefficient specified in manual&
&_interactions_block : '//trim(varname)//'. The coefficient must be 0 &
&<= c <= order.')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_coeffs(xyzzyaaac35)
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading value for "'//trim(varname)//'" in manual_interactio&
&ns block. Could not cast input to float.')
else
call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Unrecognized keyword &
&in manual_interactions block: '//trim(varname)//'. Allowed PSEUDODIPO&
&LE parameters are "order", "cutoff", "d^2" and "coefficients".')
endif
end select
enddo
end subroutine xyzzyaafz1
subroutine xyzzyaaga1(nlines)
implicit none
integer, intent(in) :: nlines
integer xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36
logical xyzzyaaad36,xyzzyaaae36
character(80) varname,data_buffer
if(nlines<6)call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE', 'Wron&
&g number of lines in manual_interaction block in input file. Read '//&
&trim(i2s(nlines))//' lines but expected at least 5.')
xyzzyaaad36=.false.
do xyzzyaaaa36=2,nlines
if(.not.xyzzyaagb1(block_data(xyzzyaaaa36),varname,data_buffer))call e&
&rrstop_master('READ_MAN_INT_PSEUDODIPOLE','Syntax error reading line &
&'//trim(i2s(xyzzyaaaa36))//' of manual interactions block in input. S&
&yntax should be "keyword : value".')
if(trim(switch_case(varname))=='ORDER')then
xyzzyaaad36=.true.
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_order
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading value for "order" in manual_interactions block. Coul&
&d not cast input to integer.')
allocate(man_int_poly_coeffs(0:man_int_poly_order),stat=xyzzyaaab36)
call check_alloc(xyzzyaaab36,'READ_MAN_INT_PSEUDODIPOLE','coefficients&
&')
man_int_poly_coeffs=0.d0
exit
endif
enddo
xyzzyaaae36=.false.
do xyzzyaaaa36=2,nlines
if(.not.xyzzyaagb1(block_data(xyzzyaaaa36),varname,data_buffer))call e&
&rrstop_master('READ_MAN_INT_PSEUDODIPOLE','Syntax error reading line &
&'//trim(i2s(xyzzyaaaa36))//' of manual interactions block in input. S&
&yntax should be "keyword : value".')
if(trim(switch_case(varname))=='TILT_ORDER')then
xyzzyaaae36=.true.
read(data_buffer,*,iostat=xyzzyaadl1)man_int_tilt_poly_order
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE','&
&Syntax error reading value for "tilt_order" in manual_interactions bl&
&ock. Could not cast input to integer.')
allocate(man_int_tilt_poly_coeffs(0:man_int_tilt_poly_order),stat=xyzz&
&yaaab36)
call check_alloc(xyzzyaaab36,'READ_MAN_INT_PSEUDODIPOLE','coefficients&
&')
man_int_tilt_poly_coeffs=0.d0
exit
endif
enddo
if(.not.xyzzyaaad36)call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE&
&','Failed to find keyword "order" in manual_interactions block.')
if(.not.xyzzyaaae36)call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE&
&','Failed to find keyword "tilt_order" in manual_interactions block.'&
&)
do xyzzyaaaa36=2,nlines
if(.not.xyzzyaagb1(block_data(xyzzyaaaa36),varname,data_buffer))call e&
&rrstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE','Syntax error reading l&
&ine '//trim(i2s(xyzzyaaaa36))//' of manual interactions block in inpu&
&t. Syntax should be "keyword : value".')
select case(trim(switch_case(varname)))
case('ORDER')
continue
case('TILT_ORDER')
continue
case('CUTOFF')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_cutoff
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE','&
&Syntax error reading value for "cutoff" in manual_interactions block.&
& Could not cast input to float.')
case('D^2')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_params(1)
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE','&
&Syntax error reading value for "d^2" in manual_interactions block. Co&
&uld not cast input to float.')
case('THETA')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_params(2)
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE','&
&Syntax error reading value for "theta" in manual_interactions block. &
&Could not cast input to float.')
case default
if(varname(1:2)=='c_'.or.varname(1:2)=='C_')then
read(varname(3:),*,iostat=xyzzyaadl1)xyzzyaaac36
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading line '//trim(i2s(xyzzyaaaa36))//' of manual interact&
&ions block in input. Failed to coerce coefficient to integer.')
if(xyzzyaaac36<0.or.xyzzyaaac36>man_int_poly_order)call errstop_master&
&('READ_MAN_INT_PSEUDODIPOLE','Invalid coefficient specified in manual&
&_interactions_block : '//trim(varname)//'. The coefficient must be 0 &
&<= c <= order.')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_poly_coeffs(xyzzyaaac36)
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading value for "'//trim(varname)//'" in manual_interactio&
&ns block. Could not cast input to float.')
elseif(varname(1:2)=='b_'.or.varname(1:2)=='B_')then
read(varname(3:),*,iostat=xyzzyaadl1)xyzzyaaac36
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_PSEUDODIPOLE','Synt&
&ax error reading line '//trim(i2s(xyzzyaaaa36))//' of manual interact&
&ions block in input. Failed to coerce coefficient to integer.')
if(xyzzyaaac36<0.or.xyzzyaaac36>man_int_tilt_poly_order)call errstop_m&
&aster('READ_MAN_INT_2D_TILT_PSDIPOLE','Invalid coefficient specified &
&in manual_interactions_block : '//trim(varname)//'. The coefficient m&
&ust be 0 <= b <= tilt_order.')
read(data_buffer,*,iostat=xyzzyaadl1)man_int_tilt_poly_coeffs(xyzzyaaa&
&c36)
if(xyzzyaadl1/=0)call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE','&
&Syntax error reading value for "'//trim(varname)//'" in manual_intera&
&ctions block. Could not cast input to float.')
else
call errstop_master('READ_MAN_INT_2D_TILT_PSDIPOLE','Unrecognized keyw&
&ord in manual_interactions block: '//trim(varname)//'. Allowed 2D_til&
&t_pseudodipole parameters are "order", "cutoff", "d^2", "theta", "til&
&t_order" and "coefficients".')
endif
end select
enddo
end subroutine xyzzyaaga1
logical function xyzzyaagb1(line,varname,data_buffer) result(succ)
implicit none
character(80),intent(in) :: line
character(80),intent(out) :: varname, data_buffer
character(80) dummy
integer xyzzyaaaa37
succ=.true.
read(line,*,iostat=xyzzyaaaa37)varname,dummy,data_buffer
if(xyzzyaaaa37/=0)then
read(line,*,iostat=xyzzyaaaa37)varname,data_buffer
if(xyzzyaaaa37/=0)succ=.false.
else
select case(trim(dummy))
case(":","=")
continue
case default
succ=.false.
endselect
endif
varname=adjustl(varname)
data_buffer=adjustl(data_buffer)
end function xyzzyaagb1
subroutine xyzzyaagc1()
implicit none
if(am_master)then
call wout()
call wout('*     *     *     *     *     *     *     *     *     *    &
& *     *')
call wout()
endif
end subroutine xyzzyaagc1
subroutine xyzzyaagd1(char_var)
implicit none
character(*),intent(in) :: char_var
if(am_master)then
call wout(repeat('=',len_trim(char_var)))
call wout(trim(char_var))
call wout(repeat('=',len_trim(char_var)))
call wout()
endif
end subroutine xyzzyaagd1
end subroutine monte_carlo

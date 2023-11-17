module store
use dsp
implicit none
integer chkpoint_level
integer mc_twist_freq
integer movieplot,io_mov,movienode
integer nspin
integer netot
integer nelec
integer nemax
integer ndet
integer o
integer opt_cycle
integer opt_cycles
integer opt_info
integer opt_maxiter,opt_maxeval
logical opt_jastrow,opt_backflow,opt_det_coeff
logical opt_orbitals,opt_geminal
logical opt_fixnl
logical vmc_cfg_by_cfg
logical dmc_cfg_by_cfg
logical orbbuf,jasbuf
logical wout_inhibit_node
integer real1_complex2
integer,parameter :: re=1,im=2
real(dp) :: k_offset(3)=0.d0
integer virtual_node
integer virtual_nconfig
integer virtual_nnodes
integer ranluxlevel,ranprint,ranlog_unit
real(dp) inv_netot,sqrt_netot
integer nitot_netot,netot_netot,four_nitot_netot,four_netot_netot,four&
&_nitot,four_netot,three_nemax,three_netot,four_netot_3,size_dbar,size&
&_det,three_nemax_2,three_nitot,three_nitot_netot,three_nitot_2,three_&
&netot_2,three_nitot_netot_2,nemax2_nspin_ndet,nine_netot,nine_nitot,t&
&hree_netot_netot,nine_netot_netot,nine_netot_8,size_fi_prod_det,size_&
&fidet,size_prod_lapdet,nspin_ndet,size_rpsi,size_grad,size_sderivs,si&
&ze_dsmat,size_d2smat,size_onedbar,six_nemax,size_harray,size_farray
integer,allocatable :: nuc_nele(:)
integer,allocatable :: neion(:,:)
integer,allocatable :: nele(:),three_nele(:),six_nele(:)
real(sp) block_time
real(sp) max_cpu_time,max_real_time
real(dp) ewald_control,constant_energy,self_term
real(dp) ke_corr,xc_corr,chi_squared_sf,finite_size_const_c
real(dp) orb_norm
real(dp) man_int_params(2)
logical man_int_op_spins
character(20) interaction
character(20) int_name
character(20) output_file
character(20) runtype
character(20) stop_method
logical complex_wf
logical popstats
logical newrun
logical fix_cutoffs
logical use_future
logical forces
logical calc_pulay
logical electron_system
logical open_unit(99)
logical use_backflow,use_jastrow
logical use_orbmods
logical noncoll_spin
logical have_jastrow3
logical sparse
logical bf_sparse
logical use_tmove
logical finite_size_corr
logical mc_twist_av
logical prefer_short_cusp
logical opt_strict
logical have_biex3pot
logical fix_holes
logical allow_nochi_atoms
logical dmc_md
logical rng_restart_safe
logical hartree_xc
logical pairing_wf,calc_field,init_by_ion,init_by_iontype,esupercell, &
&vm_smooth_limits,use_expot,old_input,new_input,use_blocktime
logical isvmc,isdmc,isvmc_dmc,isopt,isvmc_opt,isgen_mpc,isopt_vmc,isdm&
&c_dmc,isrmc,isrmc_rmc,isdmc_stats,isdmc_equil,isdmc_old,isvmc_dmc_equ&
&il
logical makemovie,movie_supercells
logical,allocatable :: ee_cusp_in_orbital(:,:)
integer max_spin_singles,max_spin_pairs,max_spin_triplets,levels_ssing&
&les,levels_spairs,levels_striplets,level_eqvfamily,no_eqvfamilies,lev&
&el_family,no_families,custom_ssingles,custom_spairs,custom_striplets,&
&no_difftypes,heg_nlayers
real(dp),allocatable :: pmasstensor(:,:,:),inv_pmasstensor(:,:,:),pmas&
&s(:),inv_pmass(:),pcharge(:),pspin(:),difftype_mass(:),fam_charge(:)
logical,allocatable :: pisotropic(:),pfermion(:),missing_det(:,:),upda&
&te_by_column(:,:),pinfmass(:)
integer,allocatable :: heg_nele(:),heg_orbtype(:,:),popp_spin(:),upd_s&
&pin(:,:),no_ssingles(:),no_spairs(:),no_striplets(:),which_ssingle(:,&
&:),which_spair(:,:,:),which_striplet(:,:,:,:),which_eqvfam(:),which_f&
&am(:),which_spin(:),which_ie(:),which_ii(:,:),which_ee(:,:),eq_triple&
&t(:,:,:,:),which_difftype(:),heg_layer(:),heg_slatt(:,:)
character(20),allocatable :: pname(:)
integer spin_dep_wc,type_wc
logical ferromagnetic
real(dp) sdw_init_polarization
real(dp),allocatable :: site_pos_in_cell(:,:,:,:),wc_gauss_exp(:),heg_&
&ylayer(:),heg_zlayer(:)
real(dp) rstar
logical localized_orbitals
logical bsmooth
integer :: which_d2index(3,3)
data which_d2index /1,4,5,4,2,6,5,6,3/
logical use_altsamp
integer altsamp,simplepdf
real(dp) vmc_optimum_e0,vmc_optimum_ew
integer detstart,detstop,ndet_smp
logical isitcomplex
integer scr_nl_2
integer on_top_ii,on_top_jj
end module store

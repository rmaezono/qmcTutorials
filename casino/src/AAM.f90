module slaarnaam
use slaarnaag
use dsp
use slaarnaap
use slaarnabk
use slaarnabt
use slaarnaca
use slaarnace
use slaarnach
use store
use slaarnaco
use slaarnacs
use slaarnaaq,     only : setup_emt,eval_dipole_moment,dipole_calc,con&
&tact_den_calc,eval_contact_den
use slaarnabg,   only : nitot,isperiodic,iontype,ignore_ionic_interact&
&ions
use slaarnabm,   only : use_magnetic_field,eval_magnetic_energy
use slaarnabs,  only : v_non_local,calc_nl_projection
use run_control,only : timer,check_alloc,errstop,errstop_master
implicit none
private
public energy_scratch_request,eval_local_energy,setup_energy_utils,fin&
&ish_energy_utils,add_config_energy_items,ionic_potential,ederiv_scrat&
&ch_request,prepare_ederivs,eval_energy_nderiv,setup_storage_energy,fi&
&nish_storage_energy,load_from_storage_energy,save_to_storage_energy,i&
&nvalidate_params_energy,enumerate_plot_energy,query_plot_energy,get_p&
&lot_energy,finish_plot_energy,eval_energy_aderiv
public n_ecomp
integer n_ecomp
integer,public :: i_eloc_def,i_pei,i_kei,i_ti,i_fisq,i_pote,i_potil,i_&
&potinl,i_vcpp_ei,i_vcpp_e,i_vcpp_ee,i_ecpp_nl,i_ecpp_perfield,i_ecpp_&
&tot,i_mpc,i_short,i_long,i_eloc_mpc,i_emasspol,i_emassvel,i_edarwin_e&
&n,i_edarwin_ee,i_eretard,i_reltot,i_dipole1,i_dipole2,i_dipole3,i_dip&
&ole_sq,i_esqr,i_contact_den,i_wght
integer,parameter :: xyzzyaaaa1=50
integer :: xyzzyaaab1(xyzzyaaaa1)=0,xyzzyaaac1(xyzzyaaaa1)=0,xyzzyaaad&
&1(xyzzyaaaa1)=0,xyzzyaaae1(xyzzyaaaa1),xyzzyaaaf1=0
integer,allocatable :: xyzzyaaag1(:),xyzzyaaah1(:),xyzzyaaai1(:)
real(dp),allocatable :: xyzzyaaaj1(:),xyzzyaaak1(:),xyzzyaaal1(:,:),xy&
&zzyaaam1(:,:,:,:)
logical,allocatable :: xyzzyaaan1(:),xyzzyaaao1(:),xyzzyaaap1(:),xyzzy&
&aaaq1(:),xyzzyaaar1(:),xyzzyaaas1(:),xyzzyaaat1(:)
logical,allocatable :: xyzzyaaau1(:),xyzzyaaav1(:),xyzzyaaaw1(:),xyzzy&
&aaax1(:),xyzzyaaay1(:)
integer,parameter :: xyzzyaaaz1=5,xyzzyaaba1=1,xyzzyaabb1=2,xyzzyaabc1&
&=3,xyzzyaabd1=4,xyzzyaabe1=5
integer :: xyzzyaabf1(xyzzyaaaz1)=0
character(64),parameter :: xyzzyaabg1(xyzzyaaaz1)=(/'energy','expot ',&
&'kei   ','ti    ','fisq  '/)
character(64),parameter :: xyzzyaabh1(xyzzyaaaz1)=(/'local energy     &
&                    ',  'electron-ion/external local potential',  'lo&
&cal kinetic energy                 ',  'local kinetic energy componen&
&t TI    ',  'local kinetic energy component FISQ  '/)
contains
subroutine energy_scratch_request(is0,nl_with_drift)
implicit none
integer,intent(inout) :: is0
logical,intent(in),optional :: nl_with_drift
integer xyzzyaaaa2
xyzzyaaaf1=xyzzyaaaf1+1
call scratch_request(kinetic=is0)
xyzzyaaab1(xyzzyaaaf1)=is0
if(have_ppots)then
xyzzyaaaa2=0
call scratch_request(ratio1_from=is0,ratio1_to=xyzzyaaaa2)
if(present(nl_with_drift))then
if(nl_with_drift)call scratch_request(drift=xyzzyaaaa2)
endif
xyzzyaaac1(xyzzyaaaf1)=xyzzyaaaa2
endif
end subroutine energy_scratch_request
subroutine ederiv_scratch_request(is0,ignore,wfn_deriv)
implicit none
integer,intent(inout) :: is0
logical,intent(in) :: ignore(:)
logical,intent(in),optional :: wfn_deriv
integer xyzzyaaaa3,xyzzyaaab3
logical xyzzyaaac3
xyzzyaaac3=.false.
if(present(wfn_deriv))xyzzyaaac3=wfn_deriv
call energy_scratch_request(is0)
if(any_analytical_deriv(ignore))then
call scratch_request(kinetic_aderiv=is0)
if(xyzzyaaac3)call scratch_request(wfn_aderiv=is0)
endif
if(.true..or.any_numerical_deriv(ignore))then
call scratch_request(kinetic_detail=is0)
xyzzyaaaa3=0
call scratch_request(kinetic=xyzzyaaaa3)
call scratch_request(kinetic_detail=xyzzyaaaa3)
if(xyzzyaaac3)then
call scratch_request(wfn_detail=is0)
call scratch_request(wfn_detail=xyzzyaaaa3)
endif
xyzzyaaad1(xyzzyaaaf1)=xyzzyaaaa3
if(have_ppots.and..not.opt_fixnl)then
xyzzyaaab3=0
call scratch_request(ratio1_from=xyzzyaaaa3,ratio1_to=xyzzyaaab3)
call scratch_protect(is0,xyzzyaaab3)
xyzzyaaae1(xyzzyaaaf1)=xyzzyaaab3
endif
endif
end subroutine ederiv_scratch_request
subroutine setup_energy_utils
implicit none
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4
allocate(xyzzyaaag1(nscratch),xyzzyaaah1(nscratch),xyzzyaaai1(nscratch&
&),stat=xyzzyaaaf4)
xyzzyaaag1=0
xyzzyaaah1=0
xyzzyaaai1=0
do xyzzyaaaa4=1,xyzzyaaaf1
xyzzyaaab4=xyzzyaaab1(xyzzyaaaa4)
call which_scratch(xyzzyaaab4)
xyzzyaaac4=xyzzyaaac1(xyzzyaaaa4)
call which_scratch(xyzzyaaac4)
xyzzyaaad4=xyzzyaaad1(xyzzyaaaa4)
call which_scratch(xyzzyaaad4)
xyzzyaaae4=xyzzyaaae1(xyzzyaaaa4)
call which_scratch(xyzzyaaae4)
xyzzyaaag1(xyzzyaaab4)=xyzzyaaac4
xyzzyaaah1(xyzzyaaab4)=xyzzyaaad4
xyzzyaaai1(xyzzyaaab4)=xyzzyaaae4
enddo
i_eloc_def=1
i_pei=2
i_kei=3
i_ti=4
i_fisq=5
i_pote=6
i_potil=7
i_potinl=8
n_ecomp=8
i_esqr=0
i_mpc=0
i_short=0
i_long=0
i_eloc_mpc=0
i_ecpp_tot=0
i_vcpp_ei=0
i_vcpp_e=0
i_vcpp_ee=0
i_ecpp_nl=0
i_emasspol=0
i_emassvel=0
i_edarwin_en=0
i_edarwin_ee=0
i_eretard=0
i_reltot=0
i_dipole1=0
i_dipole2=0
i_dipole3=0
i_dipole_sq=0
i_contact_den=0
if(popstats)then
i_esqr=n_ecomp+1
n_ecomp=n_ecomp+1
endif
if(have_veep)then
i_vcpp_ei=n_ecomp+1
i_vcpp_e=n_ecomp+2
i_vcpp_ee=n_ecomp+3
i_ecpp_tot=n_ecomp+4
n_ecomp=n_ecomp+4
if(isperiodic)then
i_ecpp_nl=n_ecomp+1
i_ecpp_perfield=n_ecomp+2
n_ecomp=n_ecomp+2
endif
endif
if(interaction_mpc_present)then
i_mpc=n_ecomp+1
i_short=n_ecomp+2
i_long=n_ecomp+3
i_eloc_mpc=n_ecomp+4
n_ecomp=n_ecomp+4
endif
if(relativistic)then
if(eval_maspol)then
n_ecomp=n_ecomp+1
i_emasspol=n_ecomp
endif
if(eval_masvel)then
n_ecomp=n_ecomp+1
i_emassvel=n_ecomp
endif
if(eval_darwinen)then
n_ecomp=n_ecomp+1
i_edarwin_en=n_ecomp
endif
if(eval_darwinee)then
n_ecomp=n_ecomp+1
i_edarwin_ee=n_ecomp
endif
if(eval_retard)then
n_ecomp=n_ecomp+1
i_eretard=n_ecomp
endif
n_ecomp=n_ecomp+1
i_reltot=n_ecomp
endif
if(eval_dipole_moment)then
i_dipole1=n_ecomp+1
i_dipole2=n_ecomp+2
i_dipole3=n_ecomp+3
i_dipole_sq=n_ecomp+4
n_ecomp=n_ecomp+4
endif
if(eval_contact_den)then
i_contact_den=n_ecomp+1
n_ecomp=n_ecomp+1
endif
if(use_altsamp)then
i_wght=n_ecomp+1
n_ecomp=n_ecomp+1
endif
allocate(etot_scr(nscratch),keimag_scr(nscratch),etot_isnan_scr(nscrat&
&ch),etot_isinf_scr(nscratch),ecomps_scr(n_ecomp,nscratch),ke_isnan_sc&
&r(nscratch),ke_isinf_scr(nscratch),nl_isnan_scr(nscratch),nl_isinf_sc&
&r(nscratch),erest_isinf_scr(nscratch),stat=xyzzyaaaf4)
call check_alloc(xyzzyaaaf4,'SETUP_ENERGY_UTILS','energy')
etot_scr=0.d0
keimag_scr=0.d0
etot_isnan_scr=.false.
etot_isinf_scr=.false.
ecomps_scr=0.d0
ke_isnan_scr=.false.
ke_isinf_scr=.false.
nl_isnan_scr=.false.
nl_isinf_scr=.false.
erest_isinf_scr=.false.
allocate(etot_valid(nscratch),ke_valid(nscratch),nl_valid(nscratch),er&
&est_valid(nscratch),stat=xyzzyaaaf4)
call check_alloc(xyzzyaaaf4,'SETUP_ENERGY_UTILS','energy_valid')
etot_valid=.false.
ke_valid=.false.
nl_valid=.false.
erest_valid=.false.
if(eval_dipole_moment)call setup_emt
if(have_ppots)then
allocate(grid_angles_scr(3,nitot,netot,nscratch),stat=xyzzyaaaf4)
call check_alloc(xyzzyaaaf4,'SETUP_ENERGY_UTILS','grid_angles')
grid_angles_scr=0.d0
endif
allocate(grid_angles_valid(nscratch),stat=xyzzyaaaf4)
call check_alloc(xyzzyaaaf4,'SETUP_ENERGY_UTILS','grid_angles_valid')
grid_angles_valid=.false.
end subroutine setup_energy_utils
subroutine finish_energy_utils
implicit none
deallocate(etot_scr,keimag_scr,ecomps_scr,etot_isnan_scr,etot_isinf_sc&
&r,ke_isnan_scr,ke_isinf_scr,nl_isnan_scr,nl_isinf_scr,erest_isinf_scr&
&)
deallocate(etot_valid,ke_valid,nl_valid,erest_valid)
deallocate(xyzzyaaag1,xyzzyaaah1,xyzzyaaai1)
if(allocated(grid_angles_scr))deallocate(grid_angles_scr)
deallocate(grid_angles_valid)
xyzzyaaab1=0
xyzzyaaac1=0
xyzzyaaad1=0
xyzzyaaae1=0
xyzzyaaaf1=0
end subroutine finish_energy_utils
subroutine setup_storage_energy(nconfig)
use slaarnaaf
implicit none
integer,intent(in) :: nconfig
integer xyzzyaaaa6
allocate(xyzzyaaaj1(nconfig),xyzzyaaak1(nconfig),xyzzyaaan1(nconfig),x&
&yzzyaaao1(nconfig),xyzzyaaal1(n_ecomp,nconfig),xyzzyaaap1(nconfig),xy&
&zzyaaaq1(nconfig),xyzzyaaar1(nconfig),xyzzyaaas1(nconfig),xyzzyaaat1(&
&nconfig),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'SETUP_STORAGE_ENERGY','e_store')
xyzzyaaaj1=0.d0
xyzzyaaak1=0.d0
xyzzyaaan1=.false.
xyzzyaaao1=.false.
xyzzyaaal1=0.d0
xyzzyaaap1=.false.
xyzzyaaaq1=.false.
xyzzyaaar1=.false.
xyzzyaaas1=.false.
xyzzyaaat1=.false.
allocate(xyzzyaaau1(nconfig),xyzzyaaav1(nconfig),xyzzyaaaw1(nconfig),x&
&yzzyaaax1(nconfig),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'SETUP_STORAGE_ENERGY','e_valid')
xyzzyaaau1=.false.
xyzzyaaav1=.false.
xyzzyaaaw1=.false.
xyzzyaaax1=.false.
if(have_ppots.and..not.opt_fixnl)then
allocate(xyzzyaaam1(3,nitot,netot,nconfig),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'SETUP_STORAGE_ENERGY','grid_angles_store'&
&)
xyzzyaaam1=0.d0
endif
allocate(xyzzyaaay1(nconfig),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'SETUP_STORAGE_ENERGY','grid_angles_svalid&
&')
xyzzyaaay1=.false.
if(allocated(etot_config))then
call dcopy(nconfig,etot_config(1),1,xyzzyaaaj1(1),1)
xyzzyaaau1(1:nconfig)=.true.
endif
if(have_ppots.and.allocated(nltot_config))then
xyzzyaaal1(i_potinl,1:nconfig)=nltot_config(1:nconfig)
xyzzyaaaw1(1:nconfig)=.true.
endif
end subroutine setup_storage_energy
subroutine finish_storage_energy
implicit none
deallocate(xyzzyaaaj1,xyzzyaaak1,xyzzyaaan1,xyzzyaaao1,xyzzyaaal1,xyzz&
&yaaap1,xyzzyaaaq1,xyzzyaaar1,xyzzyaaas1,xyzzyaaat1)
deallocate(xyzzyaaau1,xyzzyaaav1,xyzzyaaaw1,xyzzyaaax1)
if(allocated(xyzzyaaam1))deallocate(xyzzyaaam1)
deallocate(xyzzyaaay1)
end subroutine finish_storage_energy
subroutine load_from_storage_energy(is,icfg)
implicit none
integer,intent(in) :: is,icfg
etot_scr(is)=xyzzyaaaj1(icfg)
keimag_scr(is)=xyzzyaaak1(icfg)
etot_isnan_scr(is)=xyzzyaaan1(icfg)
etot_isinf_scr(is)=xyzzyaaao1(icfg)
etot_valid(is)=xyzzyaaau1(icfg)
ecomps_scr(:,is)=xyzzyaaal1(:,icfg)
ke_isnan_scr(is)=xyzzyaaap1(icfg)
ke_isinf_scr(is)=xyzzyaaaq1(icfg)
ke_valid(is)=xyzzyaaav1(icfg)
nl_isnan_scr(is)=xyzzyaaar1(icfg)
nl_isinf_scr(is)=xyzzyaaas1(icfg)
nl_valid(is)=xyzzyaaaw1(icfg)
erest_isinf_scr(is)=xyzzyaaat1(icfg)
erest_valid(is)=xyzzyaaax1(icfg)
if(allocated(xyzzyaaam1))then
grid_angles_scr(:,:,:,is)=xyzzyaaam1(:,:,:,icfg)
grid_angles_valid(is)=xyzzyaaay1(icfg)
endif
end subroutine load_from_storage_energy
subroutine save_to_storage_energy(is,icfg)
implicit none
integer,intent(in) :: is,icfg
xyzzyaaaj1(icfg)=etot_scr(is)
xyzzyaaak1(icfg)=keimag_scr(is)
xyzzyaaan1(icfg)=etot_isnan_scr(is)
xyzzyaaao1(icfg)=etot_isinf_scr(is)
xyzzyaaau1(icfg)=etot_valid(is)
xyzzyaaal1(:,icfg)=ecomps_scr(:,is)
xyzzyaaap1(icfg)=ke_isnan_scr(is)
xyzzyaaaq1(icfg)=ke_isinf_scr(is)
xyzzyaaav1(icfg)=ke_valid(is)
xyzzyaaar1(icfg)=nl_isnan_scr(is)
xyzzyaaas1(icfg)=nl_isinf_scr(is)
xyzzyaaaw1(icfg)=nl_valid(is)
xyzzyaaat1(icfg)=erest_isinf_scr(is)
xyzzyaaax1(icfg)=erest_valid(is)
if(allocated(xyzzyaaam1))xyzzyaaam1(:,:,:,icfg)=grid_angles_scr(:,:,:,&
&is)
xyzzyaaay1(icfg)=grid_angles_valid(is)
end subroutine save_to_storage_energy
subroutine xyzzyaabi1(is,js)
implicit none
integer,intent(in) :: is,js
etot_scr(js)=etot_scr(is)
keimag_scr(js)=keimag_scr(is)
etot_isnan_scr(js)=etot_isnan_scr(is)
etot_isinf_scr(js)=etot_isinf_scr(is)
etot_valid(js)=etot_valid(is)
ecomps_scr(:,js)=ecomps_scr(:,is)
ke_isnan_scr(js)=ke_isnan_scr(is)
ke_isinf_scr(js)=ke_isinf_scr(is)
ke_valid(js)=ke_valid(is)
nl_isnan_scr(js)=nl_isnan_scr(is)
nl_isinf_scr(js)=nl_isinf_scr(is)
nl_valid(js)=nl_valid(is)
erest_isinf_scr(js)=erest_isinf_scr(is)
erest_valid(js)=erest_valid(is)
if(.not.opt_fixnl.and.grid_angles_valid(is).and..not.grid_angles_valid&
&(js))then
grid_angles_scr(:,:,:,js)=grid_angles_scr(:,:,:,is)
grid_angles_valid(js)=grid_angles_valid(is)
endif
end subroutine xyzzyaabi1
subroutine xyzzyaabj1(is)
implicit none
integer,intent(in) :: is
etot_valid(is)=.false.
ke_valid(is)=.false.
if(have_ppots.and..not.opt_fixnl)nl_valid(is)=.false.
end subroutine xyzzyaabj1
subroutine invalidate_params_energy
implicit none
xyzzyaaau1(:)=.false.
xyzzyaaav1(:)=.false.
if(have_ppots.and..not.opt_fixnl)xyzzyaaaw1(:)=.false.
end subroutine invalidate_params_energy
subroutine eval_local_energy(is,etot,etot_imag,isnan,isinf,ecomps,fix_&
&nl_grid,scr_nl_override,prefetch_aderiv)
implicit none
integer,intent(in) :: is
integer,intent(in),optional :: scr_nl_override
real(dp),intent(inout),optional :: etot,etot_imag,ecomps(n_ecomp)
logical,intent(in),optional :: fix_nl_grid,prefetch_aderiv
logical,intent(out),optional :: isnan,isinf
integer xyzzyaaaa13,xyzzyaaab13
real(dp) xyzzyaaac13(netot),xyzzyaaad13(n_ecomp),xyzzyaaae13(n_ecomp),&
&xyzzyaaaf13,xyzzyaaag13,xyzzyaaah13,xyzzyaaai13(3,netot,real1_complex&
&2),xyzzyaaaj13(3),xyzzyaaak13(3)
logical xyzzyaaal13,xyzzyaaam13,xyzzyaaan13,xyzzyaaao13,xyzzyaaap13,xy&
&zzyaaaq13,xyzzyaaar13,xyzzyaaas13,xyzzyaaat13,xyzzyaaau13,xyzzyaaav13&
&,xyzzyaaaw13,xyzzyaaax13,xyzzyaaay13
call timer('EVAL_LOCAL_ENERGY',.true.)
xyzzyaaaf13=0.d0
xyzzyaaah13=0.d0
xyzzyaaae13=0.d0
xyzzyaaam13=.false.
xyzzyaaan13=.false.
xyzzyaaai13=0.d0
xyzzyaaac13=0.d0
if(have_veep.and.isperiodic)electron_fields=0.d0
xyzzyaaab13=xyzzyaaag1(is)
if(present(scr_nl_override))xyzzyaaab13=scr_nl_override
xyzzyaaay13=.false.
if(present(prefetch_aderiv))xyzzyaaay13=prefetch_aderiv
if(have_ppots)then
xyzzyaaal13=.false.
if(present(fix_nl_grid))xyzzyaaal13=fix_nl_grid
if(.not.xyzzyaaal13)then
grid_angles_valid(is)=.false.
nl_valid(is)=.false.
elseif(have_veep.and..not.erest_valid(is))then
nl_valid(is)=.false.
endif
endif
call get_rsele(is)
call prefetch_wfn(is,.true.,.true.)
xyzzyaaao13=.false.
xyzzyaaap13=.false.
xyzzyaaaq13=.false.
xyzzyaaar13=.false.
xyzzyaaas13=.false.
do xyzzyaaaa13=1,netot
call xyzzyaabk1(xyzzyaaaa13,is,xyzzyaaab13,xyzzyaaay13,xyzzyaaad13,xyz&
&zyaaag13,electron_fields,xyzzyaaaj13,xyzzyaaak13,xyzzyaaat13,xyzzyaaa&
&u13,xyzzyaaav13,xyzzyaaaw13,xyzzyaaax13)
xyzzyaaao13=xyzzyaaao13.or.xyzzyaaat13
xyzzyaaap13=xyzzyaaap13.or.xyzzyaaau13
xyzzyaaaq13=xyzzyaaaq13.or.xyzzyaaav13
xyzzyaaar13=xyzzyaaar13.or.xyzzyaaaw13
xyzzyaaas13=xyzzyaaas13.or.xyzzyaaax13
xyzzyaaae13=xyzzyaaae13+xyzzyaaad13
xyzzyaaah13=xyzzyaaah13+xyzzyaaag13
if(relativistic.or.eval_contact_den)then
xyzzyaaai13(1:3,xyzzyaaaa13,1)=xyzzyaaaj13
if(complex_wf)xyzzyaaai13(1:3,xyzzyaaaa13,2)=xyzzyaaak13
xyzzyaaac13(xyzzyaaaa13)=xyzzyaaad13(i_kei)
endif
enddo
if(have_ppots.and..not.nl_valid(is))grid_angles_valid(is)=.true.
call xyzzyaabl1(is,xyzzyaaae13,xyzzyaaaf13,xyzzyaaah13,electron_fields&
&,xyzzyaaai13,xyzzyaaac13,xyzzyaaao13,xyzzyaaap13,xyzzyaaaq13,xyzzyaaa&
&r13,xyzzyaaas13,xyzzyaaam13,xyzzyaaan13)
if(present(etot))etot=xyzzyaaaf13
if(present(etot_imag))etot_imag=xyzzyaaah13
if(present(ecomps))ecomps=xyzzyaaae13
if(present(isnan))isnan=xyzzyaaam13
if(present(isinf))isinf=xyzzyaaan13
call timer('EVAL_LOCAL_ENERGY',.false.)
end subroutine eval_local_energy
subroutine xyzzyaabk1(ii,is,js,prefetch_aderiv,ecomps1,keimag,efields,&
&drift_r,drift_i,ke_isnan,ke_isinf,nl_isnan,nl_isinf,erest_isinf)
implicit none
integer,intent(in) :: ii,is,js
real(dp),intent(out) :: ecomps1(n_ecomp),drift_r(3),drift_i(3)
real(dp),intent(inout) :: efields(3,nitot),keimag
logical,intent(in) :: prefetch_aderiv
logical,intent(out) :: ke_isnan,ke_isinf,nl_isnan,nl_isinf,erest_isinf
integer xyzzyaaaa14,xyzzyaaab14
real(dp) xyzzyaaac14(3),xyzzyaaad14(3),xyzzyaaae14,xyzzyaaaf14,xyzzyaa&
&ag14,xyzzyaaah14(3,nitot),xyzzyaaai14,xyzzyaaaj14,xyzzyaaak14,xyzzyaa&
&al14
logical xyzzyaaam14,xyzzyaaan14
complex(dp) xyzzyaaao14(3),xyzzyaaap14
ecomps1=0.d0
keimag=0.d0
ke_isnan=.false.
ke_isinf=.false.
nl_isnan=.false.
nl_isinf=.false.
erest_isinf=.false.
xyzzyaaaa14=which_spin(ii)
xyzzyaaac14=rele_scr(1:3,ii,is)
xyzzyaaab14=is
if(buffer_move_ion_from(is)/=0)xyzzyaaab14=buffer_move_ion_from(is)
if(.not.ke_valid(is).or.((relativistic.or.eval_contact_den).and..not.e&
&rest_valid(is)))then
call wfn_loggrad(ii,is,0,xyzzyaaao14,prefetch_sd=.true.,prefetch_aderi&
&v=prefetch_aderiv,isnan=xyzzyaaam14,isinf=xyzzyaaan14)
ke_isnan=ke_isnan.or.xyzzyaaam14
ke_isinf=ke_isinf.or.xyzzyaaan14
call wfn_loglap(ii,is,0,xyzzyaaap14,prefetch_aderiv=prefetch_aderiv,is&
&nan=xyzzyaaam14,isinf=xyzzyaaan14)
ke_isnan=ke_isnan.or.xyzzyaaam14
ke_isinf=ke_isinf.or.xyzzyaaan14
drift_r=dble(xyzzyaaao14)
xyzzyaaaj14=ddot(3,drift_r(1),1,drift_r(1),1)
if(complex_wf)then
drift_i=aimag(xyzzyaaao14)
xyzzyaaaj14=xyzzyaaaj14+ddot(3,drift_i(1),1,drift_i(1),1)
endif
ecomps1(i_fisq)=0.5d0*inv_pmass(xyzzyaaaa14)*xyzzyaaaj14
ecomps1(i_ti)=-0.25d0*inv_pmass(xyzzyaaaa14)*dble(xyzzyaaap14)
ecomps1(i_kei)=2.d0*ecomps1(i_ti)-ecomps1(i_fisq)
keimag=0.5d0*inv_pmass(xyzzyaaaa14)*aimag(xyzzyaaap14)
if(use_magnetic_field)then
call eval_magnetic_energy(xyzzyaaac14,xyzzyaaaa14,drift_r,drift_i,xyzz&
&yaaak14,xyzzyaaal14)
ecomps1(i_fisq)=ecomps1(i_fisq)+xyzzyaaak14
ecomps1(i_ti)=ecomps1(i_ti)+xyzzyaaak14
ecomps1(i_kei)=ecomps1(i_kei)+xyzzyaaak14
keimag=keimag+xyzzyaaal14
endif
endif
if(have_ppots.and..not.nl_valid(is))then
if(.not.grid_angles_valid(is))grid_angles_scr(:,:,ii,is)=-1.d0
call calc_nl_projection(ii,is,js,grid_angles_scr(1,1,ii,is))
endif
if(nitot/=0.and.(.not.erest_valid(is).or.(have_ppots.and..not.nl_valid&
&(is))).and..not.ignore_ionic_interactions)then
call get_eivecs(xyzzyaaab14)
call ionic_potential(ii,eivecs_scr(1,1,1,xyzzyaaab14),xyzzyaaaf14,xyzz&
&yaaag14,xyzzyaaae14,xyzzyaaah14,nl_isnan,nl_isinf,xyzzyaaan14,skip_nl&
&=nl_valid(is))
if(have_ppots)ecomps1(i_potinl)=xyzzyaaag14
ecomps1(i_potil)=xyzzyaaaf14
erest_isinf=erest_isinf.or.xyzzyaaan14
endif
if(.not.erest_valid(is))then
if(interaction_mpc_present)call compute_fourier_basis_mpc(xyzzyaaac14,&
&1)
call get_eevecs(is)
call ee_potential(ii,1,eevecs_scr(1,1,ii,is),xyzzyaaad14,.true.,xyzzya&
&aan14)
ecomps1(i_pote)=xyzzyaaad14(1)
if(interaction_mpc_present)then
ecomps1(i_short)=xyzzyaaad14(2)
ecomps1(i_long)=xyzzyaaad14(3)
endif
erest_isinf=erest_isinf.or.xyzzyaaan14
endif
if(use_expot.and..not.erest_valid(is))then
call eval_expot(xyzzyaaac14,xyzzyaaaa14,xyzzyaaae14,xyzzyaaan14)
ecomps1(i_potil)=ecomps1(i_potil)+xyzzyaaae14
erest_isinf=erest_isinf.or.xyzzyaaan14
endif
if(have_veep.and.(.not.erest_valid(is).or.(have_ppots.and..not.nl_vali&
&d(is))))then
call get_eivecs(xyzzyaaab14)
call compute_vcpp(ii,rele_scr(1,1,is),eivecs_scr(1,1,1,xyzzyaaab14),ec&
&omps1(i_vcpp_ei),ecomps1(i_vcpp_e),ecomps1(i_vcpp_ee),xyzzyaaah14,xyz&
&zyaaai14)
if(isperiodic)then
efields=efields+xyzzyaaah14
ecomps1(i_ecpp_nl)=xyzzyaaai14
endif
endif
end subroutine xyzzyaabk1
subroutine xyzzyaabl1(is,ecomps,e_tot,ke_imag,efield,dvel,relkei,ke_is&
&nan,ke_isinf,nl_isnan,nl_isinf,erest_isinf,isnan,isinf)
implicit none
integer,intent(in) :: is
real(dp),intent(inout) :: ecomps(n_ecomp),e_tot,ke_imag,efield(3,nitot&
&),relkei(netot),dvel(3,netot,real1_complex2)
logical,intent(inout) :: ke_isnan,ke_isinf,nl_isnan,nl_isinf,erest_isi&
&nf,isnan,isinf
integer xyzzyaaaa15
real(dp) xyzzyaaab15,xyzzyaaac15(3),xyzzyaaad15(3),xyzzyaaae15,xyzzyaa&
&af15,xyzzyaaag15,xyzzyaaah15,xyzzyaaai15,xyzzyaaaj15,xyzzyaaak15
if(ke_valid(is))then
ecomps(i_fisq)=ecomps_scr(i_fisq,is)
ecomps(i_ti)=ecomps_scr(i_ti,is)
ecomps(i_kei)=ecomps_scr(i_kei,is)
ke_imag=keimag_scr(is)
ke_isnan=ke_isnan_scr(is)
ke_isinf=ke_isinf_scr(is)
endif
if(nl_valid(is))then
ecomps(i_potinl)=ecomps_scr(i_potinl,is)
nl_isnan=nl_isnan_scr(is)
nl_isinf=nl_isinf_scr(is)
endif
if(.not.erest_valid(is))then
if(have_biex3pot)then
call biex3pot(rele_scr(1,1,is),xyzzyaaab15)
ecomps(i_pote)=ecomps(i_pote)+xyzzyaaab15
if(interaction_mpc_present)ecomps(i_short)=ecomps(i_short)+xyzzyaaab15
endif
if(interaction_mpc_present.and..not.homogeneous_mpc)ecomps(i_long)=eco&
&mps(i_long)-mpc_correction
if(have_veep.and.isperiodic.and..not.ignore_ionic_interactions)then
do xyzzyaaaa15=1,nitot
xyzzyaaac15=efield(1:3,xyzzyaaaa15)+ion_fields(1:3,xyzzyaaaa15)
ecomps(i_ecpp_perfield)=ecomps(i_ecpp_perfield)+cpp_prefac(xyzzyaaaa15&
&)*ddot(3,xyzzyaaac15(1),1,xyzzyaaac15(1),1)
enddo
endif
else
ecomps(i_potil)=ecomps_scr(i_potil,is)
ecomps(i_pote)=ecomps_scr(i_pote,is)
if(interaction_mpc_present)then
ecomps(i_short)=ecomps_scr(i_short,is)
ecomps(i_long)=ecomps_scr(i_long,is)
endif
if(have_veep)then
ecomps(i_vcpp_ei)=ecomps_scr(i_vcpp_ei,is)
if(.not.have_ppots.or.nl_valid(is))ecomps(i_vcpp_e)=ecomps_scr(i_vcpp_&
&e,is)
ecomps(i_vcpp_ee)=ecomps_scr(i_vcpp_ee,is)
if(isperiodic)then
ecomps(i_ecpp_perfield)=ecomps_scr(i_ecpp_perfield,is)
if(.not.have_ppots.or.nl_valid(is))ecomps(i_ecpp_nl)=ecomps_scr(i_ecpp&
&_nl,is)
endif
endif
erest_isinf=erest_isinf_scr(is)
endif
if(interaction_mpc_present)ecomps(i_mpc)=ecomps(i_short)+ecomps(i_long&
&)
if(have_veep)then
ecomps(i_ecpp_tot)=ecomps(i_vcpp_ei)+ecomps(i_vcpp_e)+ecomps(i_vcpp_ee&
&)
if(isperiodic)ecomps(i_ecpp_tot)=ecomps(i_ecpp_tot)+ecomps(i_ecpp_perf&
&ield)+ecomps(i_ecpp_nl)
endif
xyzzyaaab15=ecomps(i_kei)+ecomps(i_potil)+ecomps(i_potinl)+constant_en&
&ergy
if(have_veep)xyzzyaaab15=xyzzyaaab15+ecomps(i_ecpp_tot)
ecomps(i_eloc_def)=xyzzyaaab15+ecomps(i_pote)
if(interaction_mpc_present)ecomps(i_eloc_mpc)=xyzzyaaab15+ecomps(i_mpc&
&)
isnan=ke_isnan.or.nl_isnan
isinf=ke_isinf.or.nl_isinf.or.erest_isinf
ecomps(i_pei)=ecomps(i_potil)+ecomps(i_potinl)
if(interaction_mpc_use)then
ecomps(i_pei)=ecomps(i_pei)+ecomps(i_mpc)
e_tot=ecomps(i_eloc_mpc)
else
ecomps(i_pei)=ecomps(i_pei)+ecomps(i_pote)
e_tot=ecomps(i_eloc_def)
endif
if(have_veep)ecomps(i_pei)=ecomps(i_pei)+ecomps(i_ecpp_tot)
if(popstats)ecomps(i_esqr)=e_tot*e_tot
if(eval_dipole_moment)then
if(.not.erest_valid(is))then
call dipole_calc(rele_scr(1,1,is),xyzzyaaad15)
ecomps(i_dipole1)=xyzzyaaad15(1)
ecomps(i_dipole2)=xyzzyaaad15(2)
ecomps(i_dipole3)=xyzzyaaad15(3)
ecomps(i_dipole_sq)=dot_product(xyzzyaaad15,xyzzyaaad15)
else
ecomps(i_dipole1)=ecomps_scr(i_dipole1,is)
ecomps(i_dipole2)=ecomps_scr(i_dipole2,is)
ecomps(i_dipole3)=ecomps_scr(i_dipole3,is)
ecomps(i_dipole_sq)=ecomps_scr(i_dipole_sq,is)
endif
endif
if((relativistic.or.eval_contact_den).and..not.erest_valid(is))then
call get_eevecs(is)
if(nitot>0.and..not.ignore_ionic_interactions)call get_eivecs(is)
endif
if(relativistic)then
if(.not.erest_valid(is))then
call timer('RELATIVISTIC',.true.)
call relativistic_terms(dvel,eivecs_scr(1,1,1,is),relkei,eevecs_scr(1,&
&1,1,is),xyzzyaaae15,xyzzyaaaf15,xyzzyaaag15,xyzzyaaah15,xyzzyaaai15,x&
&yzzyaaak15)
if(eval_maspol)ecomps(i_emasspol)=xyzzyaaae15
if(eval_masvel)ecomps(i_emassvel)=xyzzyaaaf15
if(eval_darwinen)ecomps(i_edarwin_en)=xyzzyaaag15
if(eval_darwinee)ecomps(i_edarwin_ee)=xyzzyaaah15
if(eval_retard)ecomps(i_eretard)=xyzzyaaai15
ecomps(i_reltot)=xyzzyaaak15
call timer('RELATIVISTIC',.false.)
else
if(eval_maspol)ecomps(i_emasspol)=ecomps_scr(i_emasspol,is)
if(eval_masvel)ecomps(i_emassvel)=ecomps_scr(i_emassvel,is)
if(eval_darwinen)ecomps(i_edarwin_en)=ecomps_scr(i_edarwin_en,is)
if(eval_darwinee)ecomps(i_edarwin_ee)=ecomps_scr(i_edarwin_ee,is)
if(eval_retard)ecomps(i_eretard)=ecomps_scr(i_eretard,is)
ecomps(i_reltot)=ecomps_scr(i_reltot,is)
endif
endif
if(eval_contact_den)then
if(.not.erest_valid(is))then
call contact_den_calc(dvel,relkei,eevecs_scr(1,1,1,is),xyzzyaaaj15)
ecomps(i_contact_den)=xyzzyaaaj15
else
ecomps(i_contact_den)=ecomps_scr(i_contact_den,is)
endif
endif
etot_scr(is)=e_tot
keimag_scr(is)=ke_imag
etot_isnan_scr(is)=isnan
etot_isinf_scr(is)=isinf
etot_valid(is)=.true.
ecomps_scr(:,is)=ecomps(:)
ke_isnan_scr(is)=ke_isnan
ke_isinf_scr(is)=ke_isinf
ke_valid(is)=.true.
nl_isnan_scr(is)=nl_isnan
nl_isinf_scr(is)=nl_isinf
nl_valid(is)=.true.
erest_isinf_scr(is)=erest_isinf
erest_valid(is)=.true.
end subroutine xyzzyaabl1
subroutine ionic_potential(ii,eivecs,poti_local,poti_non_local,eipot,f&
&ield_at_ion,nl_isnan,nl_isinf,pot_isinf,skip_nl)
implicit none
integer,intent(in) :: ii
real(dp),intent(in) ::  eivecs(4,nitot,netot)
real(dp),intent(out) ::  poti_local,poti_non_local,eipot,field_at_ion(&
&3,nitot)
logical,intent(in),optional :: skip_nl
logical,intent(out) :: nl_isnan,nl_isinf,pot_isinf
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16
real(dp) xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,xyzzyaaag16,xyzzyaaah16,x&
&yzzyaaai16
real(dp),allocatable,save :: xyzzyaaaj16(:)
logical xyzzyaaak16
logical,save :: xyzzyaaal16=.true.
if(xyzzyaaal16)then
if(isperiodic)then
allocate(xyzzyaaaj16(nitot))
do xyzzyaaaa16=1,nitot
xyzzyaaaj16(xyzzyaaaa16)=zion(iontype(xyzzyaaaa16))
enddo
endif
xyzzyaaal16=.false.
endif
xyzzyaaad16=0.d0
xyzzyaaae16=0.d0
nl_isnan=.false.
nl_isinf=.false.
pot_isinf=.false.
call timer('EI_POTENTIAL',.true.)
if(isperiodic)then
if(have_ppots)then
do xyzzyaaaa16=1,nitot
xyzzyaaaf16=eivecs(4,xyzzyaaaa16,ii)
xyzzyaaab16=iontype(xyzzyaaaa16)
if(xyzzyaaaf16<rcut_loc(xyzzyaaab16))then
call lookup(pp_radial_grid(1,xyzzyaaab16),ncoeff(xyzzyaaab16),xyzzyaaa&
&f16,xyzzyaaac16)
xyzzyaaac16=min(max(xyzzyaaac16-(npoly-1)/2,1),ncoeff(xyzzyaaab16)+1-n&
&poly)
call interp_nev(pp_radial_grid(xyzzyaaac16,xyzzyaaab16),fcoeff_loc(xyz&
&zyaaac16,xyzzyaaab16),npoly,xyzzyaaaf16,xyzzyaaah16,xyzzyaaag16)
xyzzyaaad16=xyzzyaaad16+(xyzzyaaah16+zion(xyzzyaaab16)/xyzzyaaaf16)
endif
enddo
else
xyzzyaaad16=0.d0
endif
call en_potential(nitot,xyzzyaaaj16,eivecs(1,1,ii),xyzzyaaae16,field_a&
&t_ion)
else
do xyzzyaaaa16=1,nitot
xyzzyaaaf16=eivecs(4,xyzzyaaaa16,ii)
xyzzyaaab16=iontype(xyzzyaaaa16)
if(ncoeff(xyzzyaaab16)/=0)then
if(xyzzyaaaf16<rcut_loc(xyzzyaaab16))then
call lookup(pp_radial_grid(1,xyzzyaaab16),ncoeff(xyzzyaaab16),xyzzyaaa&
&f16,xyzzyaaac16)
xyzzyaaac16=min(max(xyzzyaaac16-(npoly-1)/2,1),ncoeff(xyzzyaaab16)+1-n&
&poly)
call interp_nev(pp_radial_grid(xyzzyaaac16,xyzzyaaab16),fcoeff_loc(xyz&
&zyaaac16,xyzzyaaab16),npoly,xyzzyaaaf16,xyzzyaaah16,xyzzyaaag16)
xyzzyaaad16=xyzzyaaad16+xyzzyaaah16
else
xyzzyaaad16=xyzzyaaad16-zion(xyzzyaaab16)/xyzzyaaaf16
endif
else
select case(iinteraction)
case(icoulomb,inone)
if(xyzzyaaaf16==0.d0)then
pot_isinf=.true.
else
xyzzyaaad16=xyzzyaaad16-zion(xyzzyaaab16)/xyzzyaaaf16
endif
case(ilogarithmic)
if(eivecs(4,xyzzyaaaa16,ii)==0.d0)then
pot_isinf=.true.
else
xyzzyaaai16=rstar+abs(eivecs(3,xyzzyaaaa16,ii))
xyzzyaaad16=xyzzyaaad16-zion(xyzzyaaab16)*(log(2.d0*xyzzyaaai16/xyzzya&
&aaf16)-euler)/xyzzyaaai16
endif
case(i2d_int)
xyzzyaaad16=xyzzyaaad16-zion(xyzzyaaab16)*pot_2d_int(xyzzyaaaf16,eivec&
&s(3,xyzzyaaaa16,ii))
case default
call errstop('EI_POTENTIAL','Currently only the logarithmic and 2D_int&
& manual interactions are available for fixed charges.')
end select
endif
enddo
endif
poti_local=xyzzyaaad16+xyzzyaaae16
if(.not.electron_system)poti_local=-poti_local*pcharge(which_spin(ii))
call timer('EI_POTENTIAL',.false.)
xyzzyaaak16=.not.have_ppots
if(present(skip_nl))xyzzyaaak16=xyzzyaaak16.or.skip_nl
if(xyzzyaaak16)then
poti_non_local=0.d0
else
call timer('EI_POTENTIAL_NONLOCAL',.true.)
call v_non_local(ii,eivecs,poti_non_local,nl_isnan,nl_isinf)
call timer('EI_POTENTIAL_NONLOCAL',.false.)
endif
if(.not.electron_system)poti_non_local=-poti_non_local*pcharge(which_s&
&pin(ii))
eipot=poti_local+poti_non_local
end subroutine ionic_potential
subroutine prepare_ederivs(is,ignore,e0,logwfn0,get_logwfn,ei0,isnan,i&
&sinf,iszero)
implicit none
integer,intent(in) :: is
real(dp),intent(out) :: e0
real(dp),intent(out),optional :: ei0
complex(dp),intent(out) :: logwfn0
logical,intent(in) :: ignore(:)
logical,intent(in),optional :: get_logwfn
logical,intent(out),optional :: isnan,isinf,iszero
integer xyzzyaaaa17
real(dp) xyzzyaaab17
logical xyzzyaaac17,nan,inf,zero
call timer('PREPARE_EDERIVS',.true.)
xyzzyaaac17=.false.
if(present(get_logwfn))xyzzyaaac17=get_logwfn
e0=0.d0
xyzzyaaab17=0.d0
logwfn0=czero
call eval_local_energy(is,etot=e0,etot_imag=xyzzyaaab17,fix_nl_grid=.t&
&rue.,isnan=nan,isinf=inf,prefetch_aderiv=.true.)
if(present(ei0))ei0=xyzzyaaab17
if(xyzzyaaac17)call wfn_logval(is,logwfn0,iszero=zero)
if(any_numerical_deriv(ignore))then
xyzzyaaaa17=xyzzyaaah1(is)
call clone_scratch(is,xyzzyaaaa17,reset=.true.)
call xyzzyaabi1(is,xyzzyaaaa17)
endif
if(present(isnan))isnan=nan
if(present(isinf))isinf=inf
if(present(iszero))iszero=zero
call timer('PREPARE_EDERIVS',.false.)
end subroutine prepare_ederivs
subroutine eval_energy_nderiv(is,params,iparam,p0,ph,ignore,e0,logwfn0&
&,e1,logwfn1,is_first,get_logwfn,ei0,ei1,isnan,isinf,iszero)
implicit none
integer,intent(in) :: is,iparam
real(dp),intent(in) :: p0,ph,e0,params(:)
real(dp),intent(in),optional :: ei0
real(dp),intent(out) :: e1
real(dp),intent(out),optional :: ei1
complex(dp),intent(in) :: logwfn0
complex(dp),intent(out) :: logwfn1
logical,intent(in) :: ignore(:),is_first
logical,intent(in),optional :: get_logwfn
logical,intent(out),optional :: isnan,isinf,iszero
integer xyzzyaaaa18,xyzzyaaab18
real(dp) xyzzyaaac18,de,xyzzyaaad18,xyzzyaaae18
complex(dp) dlwfn
logical xyzzyaaaf18,nan,inf,zero
call timer('EVAL_ENERGY_NDERIV',.true.)
xyzzyaaaf18=.false.
if(present(get_logwfn))xyzzyaaaf18=get_logwfn
e1=0.d0
logwfn1=czero
nan=.false.
inf=.false.
zero=.false.
if(param_has_aderiv(iparam,ignore))then
de=0.d0
xyzzyaaad18=0.d0
dlwfn=czero
if(xyzzyaaaf18)then
call eval_energy_aderiv(is,iparam,ignore,de,xyzzyaaad18,dlwfn=dlwfn)
else
call eval_energy_aderiv(is,iparam,ignore,de,xyzzyaaad18)
endif
xyzzyaaac18=ph-p0
call make_representable(xyzzyaaac18)
e1=e0+xyzzyaaac18*de
if(xyzzyaaaf18)logwfn1=logwfn0+log(c_one+xyzzyaaac18*dlwfn)
if(present(ei1))then
if(present(ei0))then
ei1=ei0+xyzzyaaac18*xyzzyaaad18
else
ei1=0.d0
endif
endif
else
xyzzyaaaa18=xyzzyaaah1(is)
xyzzyaaab18=xyzzyaaai1(is)
call clone_scratch(is,xyzzyaaaa18,reset=.false.)
call xyzzyaabi1(is,xyzzyaaaa18)
call put_param1(params,iparam,ph,ignore,.not.is_first)
call invalidate_param1(xyzzyaaaa18,iparam,ignore)
call xyzzyaabj1(xyzzyaaaa18)
call eval_local_energy(xyzzyaaaa18,etot=e1,etot_imag=xyzzyaaae18,fix_n&
&l_grid=.true.,scr_nl_override=xyzzyaaab18,isnan=nan,isinf=inf)
if(present(ei1))ei1=xyzzyaaae18
if(xyzzyaaaf18)call wfn_logval(xyzzyaaaa18,logwfn1,iszero=zero)
call put_params(params,ignore,.true.,restore_hint=iparam,quirky=is_fir&
&st)
call invalidate_param1(xyzzyaaaa18,iparam,ignore)
call xyzzyaabj1(xyzzyaaaa18)
endif
if(present(isnan))isnan=nan
if(present(isinf))isinf=inf
if(present(iszero))iszero=zero
call timer('EVAL_ENERGY_NDERIV',.false.)
end subroutine eval_energy_nderiv
subroutine eval_energy_aderiv(is,iparam,ignore,de,dei,dlwfn)
implicit none
integer,intent(in) :: is,iparam
real(dp),intent(inout) :: de,dei
complex(dp),intent(inout),optional :: dlwfn
logical,intent(in) :: ignore(*)
integer xyzzyaaaa19
real(dp) xyzzyaaab19(3),xyzzyaaac19(3),xyzzyaaad19(3),xyzzyaaae19(3),x&
&yzzyaaaf19,xyzzyaaag19,xyzzyaaah19
complex(dp) xyzzyaaai19(3),xyzzyaaaj19(3,netot),xyzzyaaak19(netot)
call timer('EVAL_ENERGY_ADERIV',.true.)
if(present(dlwfn))then
call wfn_aderiv(is,iparam,ignore,xyzzyaaaj19,xyzzyaaak19,dlogval_psi=d&
&lwfn)
else
call wfn_aderiv(is,iparam,ignore,xyzzyaaaj19,xyzzyaaak19)
endif
de=0.d0
dei=0.d0
do xyzzyaaaa19=1,netot
call wfn_loggrad(xyzzyaaaa19,is,0,xyzzyaaai19)
xyzzyaaab19=dble(xyzzyaaai19)
xyzzyaaac19=dble(xyzzyaaaj19(:,xyzzyaaaa19))
xyzzyaaaf19=ddot(3,xyzzyaaab19(1),1,xyzzyaaac19(1),1)
if(complex_wf)then
xyzzyaaad19=aimag(xyzzyaaai19)
xyzzyaaae19=aimag(xyzzyaaaj19(:,xyzzyaaaa19))
xyzzyaaaf19=xyzzyaaaf19+ddot(3,xyzzyaaad19(1),1,xyzzyaaae19(1),1)
endif
xyzzyaaag19=-0.5d0*inv_pmass(which_spin(xyzzyaaaa19))*(dble(xyzzyaaak1&
&9(xyzzyaaaa19))+2.d0*xyzzyaaaf19)
if(complex_wf)xyzzyaaah19=-0.5d0*inv_pmass(which_spin(xyzzyaaaa19))*ai&
&mag(xyzzyaaak19(xyzzyaaaa19))
de=de+xyzzyaaag19
if(complex_wf)dei=dei+xyzzyaaah19
enddo
call timer('EVAL_ENERGY_ADERIV',.false.)
end subroutine eval_energy_aderiv
subroutine add_config_energy_items(is0)
use slaarnaaf, only : add_config
implicit none
integer,intent(in) :: is0
real(dp) etot,xyzzyaaaa20,xyzzyaaab20
if(.not.etot_valid(is0))call eval_local_energy(is0)
etot=etot_scr(is0)
xyzzyaaaa20=ecomps_scr(i_potinl,is0)
xyzzyaaab20=ecomps_scr(i_potil,is0)
if(have_veep)xyzzyaaab20=xyzzyaaab20+ecomps_scr(i_ecpp_tot,is0)
if(interaction_mpc_use)then
xyzzyaaab20=xyzzyaaab20+ecomps_scr(i_mpc,is0)
else
xyzzyaaab20=xyzzyaaab20+ecomps_scr(i_pote,is0)
endif
call add_config(modify=.true.,etot=etot,local_potential=xyzzyaaab20,nl&
&tot=xyzzyaaaa20)
end subroutine add_config_energy_items
subroutine enumerate_plot_energy(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
integer xyzzyaaaa21,xyzzyaaab21
logical xyzzyaaac21,xyzzyaaad21
xyzzyaaad21=.not.(present(keyword).and.present(description))
xyzzyaaab21=0
xyzzyaabf1=0
do xyzzyaaaa21=1,xyzzyaaaz1
xyzzyaaac21=.false.
select case(xyzzyaaaa21)
case(xyzzyaaba1)
xyzzyaaac21=.true.
case(xyzzyaabb1)
if(use_expot.or.(nitot>0.and..not.ignore_ionic_interactions))xyzzyaaac&
&21=.true.
case(xyzzyaabc1)
xyzzyaaac21=.true.
case(xyzzyaabd1)
xyzzyaaac21=.true.
case(xyzzyaabe1)
xyzzyaaac21=.true.
end select
if(xyzzyaaac21)then
xyzzyaaab21=xyzzyaaab21+1
xyzzyaabf1(xyzzyaaab21)=xyzzyaaaa21
endif
enddo
n=xyzzyaaab21
if(.not.xyzzyaaad21)then
do xyzzyaaaa21=1,xyzzyaaab21
keyword(xyzzyaaaa21)=xyzzyaabg1(xyzzyaabf1(xyzzyaaaa21))
description(xyzzyaaaa21)=xyzzyaabh1(xyzzyaabf1(xyzzyaaaa21))
enddo
endif
end subroutine enumerate_plot_energy
subroutine query_plot_energy(iplot,rank,is_complex,has_stderr,rot_tens&
&or,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
logical count_only
count_only=.not.present(function_name)
if(iplot<0.or.iplot>xyzzyaaaz1)call errstop_master('QUERY_PLOT_ENERGY'&
&,'IPLOT out of range. Bug in calling routine.')
if(xyzzyaabf1(iplot)==0)call errstop_master('QUERY_PLOT_ENERGY','IPLOT&
& out of range. Bug in calling routine.')
select case(xyzzyaabf1(iplot))
case(xyzzyaaba1)
rank=0
is_complex=complex_wf
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Total energy'
endif
case(xyzzyaabb1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='External potential'
endif
case(xyzzyaabc1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Kinetic energy'
endif
case(xyzzyaabd1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Kinetic energy component TI'
endif
case(xyzzyaabe1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Kinetic energy component FISQ'
endif
end select
end subroutine query_plot_energy
subroutine get_plot_energy(iplot,is0,f)
implicit none
integer,intent(in) :: iplot,is0
real(dp),intent(out) :: f(*)
real(dp) ecomps(n_ecomp)
if(iplot<0.or.iplot>xyzzyaaaz1)call errstop_master('GET_PLOT_ENERGY','&
&IPLOT out of range. Bug in calling routine.')
if(xyzzyaabf1(iplot)==0)call errstop_master('GET_PLOT_ENERGY','IPLOT o&
&ut of range. Bug in calling routine.')
select case(xyzzyaabf1(iplot))
case(xyzzyaaba1)
if(.not.complex_wf)then
call eval_local_energy(is0,ecomps=ecomps,fix_nl_grid=.true.)
f(1)=ecomps(i_eloc_def)
else
call eval_local_energy(is0,ecomps=ecomps,etot_imag=f(2),fix_nl_grid=.t&
&rue.)
f(1)=ecomps(i_eloc_def)
endif
case(xyzzyaabb1)
call eval_local_energy(is0,ecomps=ecomps,fix_nl_grid=.true.)
f(1)=ecomps(i_potil)
case(xyzzyaabc1)
call eval_local_energy(is0,ecomps=ecomps,fix_nl_grid=.true.)
f(1)=ecomps(i_kei)
case(xyzzyaabd1)
call eval_local_energy(is0,ecomps=ecomps,fix_nl_grid=.true.)
f(1)=ecomps(i_ti)
case(xyzzyaabe1)
call eval_local_energy(is0,ecomps=ecomps,fix_nl_grid=.true.)
f(1)=ecomps(i_fisq)
end select
end subroutine get_plot_energy
subroutine finish_plot_energy
implicit none
xyzzyaabf1=0
end subroutine finish_plot_energy
end module slaarnaam

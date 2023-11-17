module slaarnabl
use dsp
use slaarnabh
use slaarnabx
use run_control,only : check_alloc
implicit none
private
public read_jastrow,write_jastrow,update_jastrow_casl,get_linear_basis&
&,jastrow_assess_check_kinetic,finite_size_corr_ke_jastrow,setup_jastr&
&ow_plot,enumerate_plot_jastrow,query_plot_jastrow,get_plot_jastrow,fi&
&nish_plot_jastrow,check_varmin_linjas_jastrow
public query_jastrow_levels,query_jastrow_level_details,setup_jastrow,&
&finish_jastrow,wfn_ratio_jastrow,accept_move_jastrow,reset_config_jas&
&trow,wfn_logval_jastrow,wfn_loggrad_jastrow,wfn_loglap_jastrow,prefet&
&ch_wfn_jastrow,add_config_jastrow_items,clear_scratch_jastrow,setup_j&
&astrow_params,finish_jastrow_params,get_jastrow_params,put_jastrow_pa&
&rams,clone_scratch_jastrow,invalidate_params_jastrow,invalidate_param&
&1_jastrow,wfn_aderiv_jastrow,setup_storage_jastrow,finish_storage_jas&
&trow,load_from_storage_jastrow,save_to_storage_jastrow
public gen_config_jastrow,delete_config_jastrow,copy_config_jastrow,co&
&nfig_to_pt_jastrow,pt_to_config_jastrow,redist_allocations_jastrow,re&
&dist_load_jastrow,redist_send_jastrow,redist_recv_jastrow,redist_save&
&_jastrow,redist_deallocations_jastrow,load_from_pt_jastrow,save_to_pt&
&_jastrow
public config_wfn_jastrow
public use_gjastrow,gen_gjastrow
logical use_gjastrow
logical gen_gjastrow
type config_wfn_jastrow
private
type(config_wfn_gjastrow),pointer :: pt_gjastrow
type(config_wfn_pjastrow),pointer :: pt_pjastrow
end type config_wfn_jastrow
contains
subroutine query_jastrow_levels(nlevels)
implicit none
integer,intent(out) :: nlevels
if(use_gjastrow)then
call query_gjastrow_levels(nlevels)
else
call query_pjastrow_levels(nlevels)
endif
end subroutine query_jastrow_levels
subroutine query_jastrow_level_details(level_score,level_name)
implicit none
integer,intent(out) :: level_score(*)
character(80),intent(out) :: level_name(*)
if(use_gjastrow)then
call query_gjastrow_level_details(level_score,level_name)
else
call query_pjastrow_level_details(level_score,level_name)
endif
end subroutine query_jastrow_level_details
subroutine check_varmin_linjas_jastrow
implicit none
if(use_gjastrow)then
continue
else
call check_varmin_linjas_pjastrow
endif
end subroutine check_varmin_linjas_jastrow
subroutine setup_jastrow
implicit none
if(use_gjastrow)then
call setup_gjastrow
else
call setup_pjastrow
endif
end subroutine setup_jastrow
subroutine finish_jastrow
implicit none
if(use_gjastrow)then
call finish_gjastrow
else
call finish_pjastrow
endif
end subroutine finish_jastrow
subroutine wfn_ratio_jastrow(is,js,ilevel,ratio,fd,sd)
implicit none
integer,intent(in) :: is,js,ilevel
real(dp),intent(out) :: ratio
logical,intent(in) :: fd,sd
if(use_gjastrow)then
call wfn_ratio_gjastrow(is,js,ilevel,ratio,fd,sd)
else
call wfn_ratio_pjastrow(is,js,ilevel,ratio,fd,sd)
endif
end subroutine wfn_ratio_jastrow
subroutine accept_move_jastrow(is,js)
implicit none
integer,intent(in) :: is,js
if(use_gjastrow)then
call accept_move_gjastrow(is,js)
else
call accept_move_pjastrow(is,js)
endif
end subroutine accept_move_jastrow
subroutine reset_config_jastrow(is,js)
implicit none
integer,intent(in) :: is,js
if(use_gjastrow)then
call reset_config_gjastrow(is,js)
else
call reset_config_pjastrow(is,js)
endif
end subroutine reset_config_jastrow
subroutine wfn_logval_jastrow(is,logwfn)
implicit none
integer,intent(in) :: is
real(dp),intent(out) :: logwfn
if(use_gjastrow)then
call wfn_logval_gjastrow(is,logwfn)
else
call wfn_logval_pjastrow(is,logwfn)
endif
end subroutine wfn_logval_jastrow
subroutine wfn_loggrad_jastrow(ii,is,ilevel,val,sd,aderiv,gradjas)
implicit none
integer,intent(in) :: ii,is,ilevel
real(dp),intent(out) :: gradjas(3)
logical,intent(in) :: val,sd,aderiv
if(use_gjastrow)then
call wfn_loggrad_gjastrow(ii,is,ilevel,val,sd,aderiv,gradjas)
else
call wfn_loggrad_pjastrow(ii,is,ilevel,val,sd,gradjas)
endif
end subroutine wfn_loggrad_jastrow
subroutine wfn_loglap_jastrow(ii,is,ilevel,val,fd,aderiv,lapjas)
implicit none
integer,intent(in) :: ii,is,ilevel
real(dp),intent(out) :: lapjas
logical,intent(in) :: val,fd,aderiv
if(use_gjastrow)then
call wfn_loglap_gjastrow(ii,is,ilevel,val,fd,aderiv,lapjas)
else
call wfn_loglap_pjastrow(ii,is,ilevel,val,fd,lapjas)
endif
end subroutine wfn_loglap_jastrow
subroutine prefetch_wfn_jastrow(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
if(use_gjastrow)then
call prefetch_wfn_gjastrow(is,fd,sd)
else
call prefetch_wfn_pjastrow(is,fd,sd)
endif
end subroutine prefetch_wfn_jastrow
subroutine gen_config_jastrow(pt_config)
implicit none
type(config_wfn_jastrow),pointer :: pt_config
integer xyzzyaaaa14
allocate(pt_config,stat=xyzzyaaaa14)
call check_alloc(xyzzyaaaa14,'GEN_CONFIG_WFN_JASTROW','container')
if(use_gjastrow)then
call gen_config_gjastrow(pt_config%pt_gjastrow)
nullify(pt_config%pt_pjastrow)
else
call gen_config_pjastrow(pt_config%pt_pjastrow)
nullify(pt_config%pt_gjastrow)
endif
end subroutine gen_config_jastrow
subroutine delete_config_jastrow(pt_config)
implicit none
type(config_wfn_jastrow),pointer :: pt_config
if(use_gjastrow)then
call delete_config_gjastrow(pt_config%pt_gjastrow)
else
call delete_config_pjastrow(pt_config%pt_pjastrow)
endif
deallocate(pt_config)
end subroutine delete_config_jastrow
subroutine copy_config_jastrow(pt_from,pt_to)
implicit none
type(config_wfn_jastrow),pointer :: pt_from,pt_to
if(use_gjastrow)then
call copy_config_gjastrow(pt_from%pt_gjastrow,pt_to%pt_gjastrow)
else
call copy_config_pjastrow(pt_from%pt_pjastrow,pt_to%pt_pjastrow)
endif
end subroutine copy_config_jastrow
subroutine config_to_pt_jastrow(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_jastrow),pointer :: pt_config
if(use_gjastrow)then
call config_to_pt_gjastrow(pt_config%pt_gjastrow,k)
else
call config_to_pt_pjastrow(pt_config%pt_pjastrow,k)
endif
end subroutine config_to_pt_jastrow
subroutine pt_to_config_jastrow(pt_config)
implicit none
type(config_wfn_jastrow),pointer :: pt_config
if(use_gjastrow)then
call pt_to_config_gjastrow(pt_config%pt_gjastrow)
else
call pt_to_config_pjastrow(pt_config%pt_pjastrow)
endif
end subroutine pt_to_config_jastrow
subroutine redist_allocations_jastrow(kmax)
implicit none
integer,intent(in) :: kmax
if(use_gjastrow)then
call redist_allocations_gjastrow(kmax)
else
call redist_allocations_pjastrow(kmax)
endif
end subroutine redist_allocations_jastrow
subroutine redist_load_jastrow(pt_config,k,blocking)
implicit none
integer,intent(in) :: k
logical,intent(in) :: blocking
type(config_wfn_jastrow),pointer :: pt_config
if(use_gjastrow)then
call redist_load_gjastrow(pt_config%pt_gjastrow,k,blocking)
else
call redist_load_pjastrow(pt_config%pt_pjastrow,k,blocking)
endif
end subroutine redist_load_jastrow
subroutine redist_send_jastrow(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
if(use_gjastrow)then
call redist_send_gjastrow(jnode,kbase,k,nbt,reqbase,blocking)
else
call redist_send_pjastrow(jnode,kbase,k,nbt,reqbase,blocking)
endif
end subroutine redist_send_jastrow
subroutine redist_recv_jastrow(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
if(use_gjastrow)then
call redist_recv_gjastrow(jnode,kbase,k,nbt,reqbase,blocking)
else
call redist_recv_pjastrow(jnode,kbase,k,nbt,reqbase,blocking)
endif
end subroutine redist_recv_jastrow
subroutine redist_save_jastrow(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_jastrow),pointer :: pt_config
if(use_gjastrow)then
call redist_save_gjastrow(pt_config%pt_gjastrow,k)
else
call redist_save_pjastrow(pt_config%pt_pjastrow,k)
endif
end subroutine redist_save_jastrow
subroutine redist_deallocations_jastrow
implicit none
if(use_gjastrow)then
call redist_deallocations_gjastrow
else
call redist_deallocations_pjastrow
endif
end subroutine redist_deallocations_jastrow
subroutine load_from_pt_jastrow(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_jastrow),pointer :: pt_config
if(use_gjastrow)then
call load_from_pt_gjastrow(is,pt_config%pt_gjastrow)
else
call load_from_pt_pjastrow(is,pt_config%pt_pjastrow)
endif
end subroutine load_from_pt_jastrow
subroutine save_to_pt_jastrow(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_jastrow),pointer :: pt_config
if(use_gjastrow)then
call save_to_pt_gjastrow(is,pt_config%pt_gjastrow)
else
call save_to_pt_pjastrow(is,pt_config%pt_pjastrow)
endif
end subroutine save_to_pt_jastrow
subroutine add_config_jastrow_items(is)
implicit none
integer,intent(in) :: is
if(use_gjastrow)then
call add_config_gjastrow_items(is)
else
call add_config_pjastrow_items(is)
endif
end subroutine add_config_jastrow_items
subroutine setup_storage_jastrow(nconfig,ignore)
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(:)
if(use_gjastrow)then
call setup_storage_gjastrow(nconfig,ignore)
else
call setup_storage_pjastrow(nconfig,ignore)
endif
end subroutine setup_storage_jastrow
subroutine finish_storage_jastrow
implicit none
if(use_gjastrow)then
call finish_storage_gjastrow
else
call finish_storage_pjastrow
endif
end subroutine finish_storage_jastrow
subroutine load_from_storage_jastrow(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(use_gjastrow)then
call load_from_storage_gjastrow(is,icfg)
else
call load_from_storage_pjastrow(is,icfg)
endif
end subroutine load_from_storage_jastrow
subroutine save_to_storage_jastrow(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(use_gjastrow)then
call save_to_storage_gjastrow(is,icfg)
else
call save_to_storage_pjastrow(is,icfg)
endif
end subroutine save_to_storage_jastrow
subroutine clone_scratch_jastrow(is,js)
implicit none
integer,intent(in) :: is,js
if(use_gjastrow)then
call clone_scratch_gjastrow(is,js)
else
call clone_scratch_pjastrow(is,js)
endif
end subroutine clone_scratch_jastrow
subroutine setup_jastrow_params(nparam)
implicit none
integer,intent(inout) :: nparam
if(use_gjastrow)then
call setup_gjastrow_params(nparam)
else
call setup_pjastrow_params(nparam)
endif
end subroutine setup_jastrow_params
subroutine finish_jastrow_params
implicit none
if(use_gjastrow)then
call finish_gjastrow_params
else
call finish_pjastrow_params
endif
end subroutine finish_jastrow_params
subroutine get_jastrow_params(params,has_lolim,lolim,has_hilim,hilim,i&
&s_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,l&
&abel)
implicit none
real(dp),intent(inout) :: params(:),lolim(:),hilim(:)
logical,intent(inout) :: has_lolim(:),has_hilim(:),is_shallow(:),is_re&
&dundant(:),is_linear(:),is_loglinear(:),has_aderiv(:),affect_map(:,:)
character(2),intent(inout) :: label(:)
if(use_gjastrow)then
call get_gjastrow_params(params,has_lolim,lolim,has_hilim,hilim,is_sha&
&llow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,label)
else
call get_pjastrow_params(params,has_lolim,lolim,has_hilim,hilim,is_sha&
&llow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,label)
endif
end subroutine get_jastrow_params
subroutine put_jastrow_params(params,ignore,iparam_buffer,prestore,bad&
&_params,restore_hint,pjastrow_quirk)
implicit none
integer,intent(in) :: iparam_buffer,restore_hint
real(dp),intent(inout) :: params(:)
logical,intent(in) :: ignore(:),prestore,pjastrow_quirk
logical,intent(out) :: bad_params
if(use_gjastrow)then
call put_gjastrow_params(params,ignore,iparam_buffer,prestore,bad_para&
&ms,restore_hint)
else
call put_pjastrow_params(params,ignore,iparam_buffer,prestore.and..not&
&.pjastrow_quirk,bad_params,restore_hint)
endif
end subroutine put_jastrow_params
subroutine invalidate_param1_jastrow(is,iparam)
implicit none
integer,intent(in) :: is,iparam
if(use_gjastrow)then
call invalidate_param1_gjastrow(is,iparam)
else
call invalidate_param1_pjastrow(is,iparam)
endif
end subroutine invalidate_param1_jastrow
subroutine invalidate_params_jastrow(ignore)
implicit none
logical,intent(in) :: ignore(:)
if(use_gjastrow)then
call invalidate_params_gjastrow(ignore)
else
call invalidate_params_pjastrow(ignore)
endif
end subroutine invalidate_params_jastrow
subroutine clear_scratch_jastrow(is)
implicit none
integer,intent(in) :: is
if(use_gjastrow)then
call clear_scratch_gjastrow(is)
else
call clear_scratch_pjastrow(is)
endif
end subroutine clear_scratch_jastrow
subroutine read_jastrow(empty_jastrow)
implicit none
logical,intent(inout) :: empty_jastrow
if(use_gjastrow)then
call read_gjastrow(empty_jastrow)
else
call read_pjastrow(empty_jastrow,gen_gjastrow)
endif
end subroutine read_jastrow
subroutine write_jastrow(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
if(use_gjastrow)then
continue
else
call write_pjastrow(correlation_name)
endif
end subroutine write_jastrow
subroutine update_jastrow_casl
implicit none
if(use_gjastrow)then
call update_gjastrow_casl
else
continue
endif
end subroutine update_jastrow_casl
subroutine wfn_aderiv_jastrow(is,iparam,dloggrad,dloglap,dlogval)
implicit none
integer,intent(in) :: is,iparam
complex(dp),intent(inout) :: dloggrad(*),dloglap(*)
complex(dp),intent(inout),optional :: dlogval
if(use_gjastrow)then
call wfn_aderiv_gjastrow(is,iparam,dloggrad,dloglap,dlogval)
else
continue
endif
end subroutine wfn_aderiv_jastrow
subroutine jastrow_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,ve&
&rbose)
implicit none
integer,intent(out) :: ie,jspin
real(dp),intent(in) :: eevecs(*),eivecs(*)
logical,intent(in) :: verbose
logical,intent(out) :: fail
if(use_gjastrow)then
call gjastrow_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,verbose&
&)
else
call pjastrow_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,verbose&
&)
endif
end subroutine jastrow_assess_check_kinetic
subroutine setup_jastrow_plot(makeplot,ispin,jspin,r_j,dir_i,r_i)
implicit none
integer,intent(in),optional :: ispin,jspin
real(dp),intent(in),optional :: r_j(3),dir_i(3),r_i(3)
logical,intent(in) :: makeplot
if(use_gjastrow)then
continue
else
if(present(ispin))then
call setup_pjastrow_plot(makeplot,ispin,jspin,r_j,dir_i,r_i)
else
call setup_pjastrow_plot(makeplot)
endif
endif
end subroutine setup_jastrow_plot
subroutine get_linear_basis(is,nparam,nparam_all,ignore,lbasis_grad_f,&
&lbasis_lap_f,grad_j0,lap_j0)
implicit none
integer,intent(in) :: is,nparam,nparam_all
real(dp),intent(out) :: lbasis_grad_f(*),lbasis_lap_f(*),grad_j0(*),la&
&p_j0
logical,intent(in) :: ignore(*)
if(use_gjastrow)then
call get_linear_basis_gjastrow(is,nparam,nparam_all,ignore,lbasis_grad&
&_f,lbasis_lap_f,grad_j0,lap_j0)
else
call get_linear_basis_pjastrow(is,nparam,nparam_all,ignore,lbasis_grad&
&_f,lbasis_lap_f,grad_j0,lap_j0)
endif
end subroutine get_linear_basis
subroutine finite_size_corr_ke_jastrow
implicit none
if(use_gjastrow)then
call finite_size_corr_ke_gjastrow
else
call finite_size_corr_ke_pjastrow
endif
end subroutine finite_size_corr_ke_jastrow
subroutine enumerate_plot_jastrow(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
if(use_gjastrow)then
call enumerate_plot_gjastrow(n,keyword,description)
else
call enumerate_plot_pjastrow(n,keyword,description)
endif
end subroutine enumerate_plot_jastrow
subroutine query_plot_jastrow(iplot,ii,rank,is_complex,has_stderr,rot_&
&tensor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
if(use_gjastrow)then
call query_plot_gjastrow(iplot,ii,rank,is_complex,has_stderr,rot_tenso&
&r,transl_pos,nfunctions,function_name)
else
call query_plot_pjastrow(iplot,ii,rank,is_complex,has_stderr,rot_tenso&
&r,transl_pos,nfunctions,function_name)
endif
end subroutine query_plot_jastrow
subroutine get_plot_jastrow(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
if(use_gjastrow)then
call get_plot_gjastrow(iplot,ii,is0,is1,f)
else
call get_plot_pjastrow(iplot,ii,is0,is1,f)
endif
end subroutine get_plot_jastrow
subroutine finish_plot_jastrow
implicit none
if(use_gjastrow)then
call finish_plot_gjastrow
else
call finish_plot_pjastrow
endif
end subroutine finish_plot_jastrow
end module slaarnabl

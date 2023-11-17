module slaarnabw
use dsp
use slaarnach
use store
use slaarnaag  ,only : czero,c_one
use run_control,only : check_alloc
implicit none
private
public enumerate_plot_pfaff,query_plot_pfaff,get_plot_pfaff,finish_plo&
&t_pfaff
public query_pfaff_levels,query_pfaff_level_details,setup_pfaff,finish&
&_pfaff,wfn_ratio_pfaff,accept_move_pfaff,reset_config_pfaff,wfn_logva&
&l_pfaff,wfn_loggrad_pfaff,wfn_loglap_pfaff,prefetch_wfn_pfaff,clear_s&
&cratch_pfaff,add_config_pfaff_items,setup_pfaff_params,finish_pfaff_p&
&arams,get_pfaff_params,put_pfaff_params,clone_scratch_pfaff,invalidat&
&e_params_pfaff,invalidate_param1_pfaff,setup_storage_pfaff,finish_sto&
&rage_pfaff,load_from_storage_pfaff,save_to_storage_pfaff
public gen_config_pfaff,delete_config_pfaff,    copy_config_pfaff,conf&
&ig_to_pt_pfaff,         pt_to_config_pfaff,redist_allocations_pfaff, &
& redist_load_pfaff,redist_send_pfaff,          redist_recv_pfaff,redi&
&st_save_pfaff,          redist_deallocations_pfaff,load_from_pt_pfaff&
&,save_to_pt_pfaff
public config_wfn_pfaff
integer xyzzyaaaa1
type config_wfn_pfaff
private
logical dummy_variable
end type config_wfn_pfaff
contains
subroutine query_pfaff_levels(nlevels)
implicit none
integer,intent(out) :: nlevels
nlevels=1
end subroutine query_pfaff_levels
subroutine query_pfaff_level_details(level_score,level_name)
implicit none
integer,intent(out) :: level_score(*)
character(80),intent(out) :: level_name(*)
level_score(1)=200
level_name(1)='Pfaffians'
end subroutine query_pfaff_level_details
subroutine setup_pfaff
implicit none
integer xyzzyaaaa4
do xyzzyaaaa4=1,nscratch
call clear_scratch_pfaff(xyzzyaaaa4)
enddo
end subroutine setup_pfaff
subroutine finish_pfaff
implicit none
end subroutine finish_pfaff
subroutine wfn_ratio_pfaff(is,js,ilevel,ratio,fd,sd)
implicit none
integer,intent(in) :: is,js,ilevel
complex(dp),intent(out) :: ratio
logical,intent(in) :: fd,sd
ratio=c_one
end subroutine wfn_ratio_pfaff
subroutine accept_move_pfaff(is,js)
implicit none
integer,intent(in) :: is,js
end subroutine accept_move_pfaff
subroutine reset_config_pfaff(is,js)
implicit none
integer,intent(in) :: is,js
end subroutine reset_config_pfaff
subroutine wfn_logval_pfaff(is,logwfn)
implicit none
integer,intent(in) :: is
complex(dp),intent(out) :: logwfn
logwfn=czero
end subroutine wfn_logval_pfaff
subroutine wfn_loggrad_pfaff(ii,is,ilevel,val,sd,loggrad)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: loggrad(3)
logical,intent(in) :: val,sd
loggrad=czero
end subroutine wfn_loggrad_pfaff
subroutine wfn_loglap_pfaff(ii,is,ilevel,val,fd,loglap)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: loglap
logical,intent(in) :: val,fd
loglap=czero
end subroutine wfn_loglap_pfaff
subroutine prefetch_wfn_pfaff(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
end subroutine prefetch_wfn_pfaff
subroutine setup_pfaff_params(nparam)
implicit none
integer,intent(inout) :: nparam
nparam=0
xyzzyaaaa1=0
end subroutine setup_pfaff_params
subroutine finish_pfaff_params
implicit none
end subroutine finish_pfaff_params
subroutine get_pfaff_params(params,has_lolim,lolim,has_hilim,hilim,is_&
&shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,lab&
&el)
implicit none
real(dp),intent(inout) :: params(xyzzyaaaa1),lolim(xyzzyaaaa1),hilim(x&
&yzzyaaaa1)
logical,intent(inout) :: has_lolim(xyzzyaaaa1),has_hilim(xyzzyaaaa1),i&
&s_shallow(xyzzyaaaa1),is_redundant(xyzzyaaaa1),is_linear(xyzzyaaaa1),&
&is_loglinear(xyzzyaaaa1),has_aderiv(xyzzyaaaa1),affect_map(xyzzyaaaa1&
&,xyzzyaaaa1)
character(2),intent(inout) :: label(xyzzyaaaa1)
end subroutine get_pfaff_params
subroutine put_pfaff_params(params,ignore,iparam_buffer,prestore,bad_p&
&arams)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaaaa1)
logical,intent(in) :: ignore(xyzzyaaaa1),prestore
logical,intent(out) :: bad_params
bad_params=.false.
end subroutine put_pfaff_params
subroutine invalidate_param1_pfaff(is,iparam)
implicit none
integer,intent(in) :: is,iparam
end subroutine invalidate_param1_pfaff
subroutine invalidate_params_pfaff(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaaaa1)
end subroutine invalidate_params_pfaff
subroutine clear_scratch_pfaff(is)
implicit none
integer,intent(in) :: is
end subroutine clear_scratch_pfaff
subroutine gen_config_pfaff(pt_config)
implicit none
type(config_wfn_pfaff),pointer :: pt_config
integer xyzzyaaaa20
allocate(pt_config,stat=xyzzyaaaa20)
call check_alloc(xyzzyaaaa20,'GEN_CONFIG_GEMINAL','container')
end subroutine gen_config_pfaff
subroutine delete_config_pfaff(pt_config)
implicit none
type(config_wfn_pfaff),pointer :: pt_config
deallocate(pt_config)
end subroutine delete_config_pfaff
subroutine copy_config_pfaff(pt_from,pt_to)
implicit none
type(config_wfn_pfaff),pointer :: pt_from,pt_to
end subroutine copy_config_pfaff
subroutine config_to_pt_pfaff(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_pfaff),pointer :: pt_config
end subroutine config_to_pt_pfaff
subroutine pt_to_config_pfaff(pt_config)
implicit none
type(config_wfn_pfaff),pointer :: pt_config
end subroutine pt_to_config_pfaff
subroutine redist_allocations_pfaff(kmax)
implicit none
integer,intent(in) :: kmax
end subroutine redist_allocations_pfaff
subroutine redist_load_pfaff(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_pfaff),pointer :: pt_config
end subroutine redist_load_pfaff
subroutine redist_send_pfaff(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
end subroutine redist_send_pfaff
subroutine redist_recv_pfaff(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
end subroutine redist_recv_pfaff
subroutine redist_save_pfaff(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_pfaff),pointer :: pt_config
end subroutine redist_save_pfaff
subroutine redist_deallocations_pfaff
implicit none
end subroutine redist_deallocations_pfaff
subroutine load_from_pt_pfaff(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_pfaff),pointer :: pt_config
end subroutine load_from_pt_pfaff
subroutine save_to_pt_pfaff(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_pfaff),pointer :: pt_config
end subroutine save_to_pt_pfaff
subroutine clone_scratch_pfaff(is,js)
implicit none
integer,intent(in) :: is,js
end subroutine clone_scratch_pfaff
subroutine add_config_pfaff_items(is)
implicit none
integer,intent(in) :: is
end subroutine add_config_pfaff_items
subroutine setup_storage_pfaff(nconfig,ignore)
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(xyzzyaaaa1)
end subroutine setup_storage_pfaff
subroutine finish_storage_pfaff
implicit none
end subroutine finish_storage_pfaff
subroutine load_from_storage_pfaff(is,icfg)
implicit none
integer,intent(in) :: is,icfg
end subroutine load_from_storage_pfaff
subroutine save_to_storage_pfaff(is,icfg)
implicit none
integer,intent(in) :: is,icfg
end subroutine save_to_storage_pfaff
subroutine enumerate_plot_pfaff(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
n=0
end subroutine enumerate_plot_pfaff
subroutine query_plot_pfaff(iplot,ii,rank,is_complex,has_stderr,rot_te&
&nsor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
nfunctions=0
end subroutine query_plot_pfaff
subroutine get_plot_pfaff(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
f(1)=0.d0
end subroutine get_plot_pfaff
subroutine finish_plot_pfaff
implicit none
end subroutine finish_plot_pfaff
end module slaarnabw

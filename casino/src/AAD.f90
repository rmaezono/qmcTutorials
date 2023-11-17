module slaarnaad
use dsp
use slaarnabb
use slaarnabv
implicit none
private
public init_backflow,write_backflow,plot_backflow,bf_assess_check_kine&
&tic,backflow_stats,setup_backflow_plot,enumerate_plot_bf,query_plot_b&
&f,get_plot_bf,finish_plot_bf
public setup_bf,finish_bf,get_bf_x,get_eevecs_bf,accept_move_bf,reset_&
&config_bf,loggrad_bf,loglap_bf,setup_bf_params,finish_bf_params,get_b&
&f_params,put_bf_params,invalidate_param1_bf,clear_scratch_bf,clone_sc&
&ratch_bf
public use_gbackflow
logical use_gbackflow
contains
subroutine setup_bf
implicit none
if(use_gbackflow)then
call setup_gbf
else
call setup_pbf
endif
end subroutine setup_bf
subroutine finish_bf
implicit none
if(use_gbackflow)then
call finish_gbf
else
call finish_pbf
endif
end subroutine finish_bf
subroutine accept_move_bf(is,js)
implicit none
integer,intent(in) :: is,js
if(use_gbackflow)then
call accept_move_gbf(is,js)
else
call accept_move_pbf(is,js)
endif
end subroutine accept_move_bf
subroutine reset_config_bf(is,js)
implicit none
integer,intent(in) :: is,js
if(use_gbackflow)then
call reset_config_gbf(is,js)
else
call reset_config_pbf(is,js)
endif
end subroutine reset_config_bf
subroutine loggrad_bf(ii,is,sd,farray,loggrad_psi)
implicit none
integer,intent(in) :: ii,is
real(dp),intent(in) :: farray(*)
complex(dp),intent(out) :: loggrad_psi(3)
logical,intent(in) :: sd
if(use_gbackflow)then
call loggrad_gbf(ii,is,sd,farray,loggrad_psi)
else
call loggrad_pbf(ii,is,sd,farray,loggrad_psi)
endif
end subroutine loggrad_bf
subroutine loglap_bf(ii,is,farray,harray,loggrad_psi,loglap_psi)
implicit none
integer,intent(in) :: ii,is
real(dp),intent(in) :: farray(*),harray(*)
complex(dp),intent(in) :: loggrad_psi(3)
complex(dp),intent(out) :: loglap_psi
if(use_gbackflow)then
call loglap_gbf(ii,is,farray,harray,loggrad_psi,loglap_psi)
else
call loglap_pbf(ii,is,farray,harray,loggrad_psi,loglap_psi)
endif
end subroutine loglap_bf
subroutine get_bf_x(is,val,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: val,fd,sd
if(use_gbackflow)then
call get_gbf_x(is,val,fd,sd)
else
call get_pbf_x(is,val,fd,sd)
endif
end subroutine get_bf_x
subroutine get_eevecs_bf(is)
implicit none
integer,intent(in) :: is
if(use_gbackflow)then
call get_eevecs_gbf(is)
else
call get_eevecs_pbf(is)
endif
end subroutine get_eevecs_bf
subroutine invalidate_param1_bf(is,iparam)
implicit none
integer,intent(in) :: is,iparam
if(use_gbackflow)then
call invalidate_param1_gbf(is,iparam)
else
call invalidate_param1_pbf(is,iparam)
endif
end subroutine invalidate_param1_bf
subroutine clear_scratch_bf(is)
implicit none
integer,intent(in) :: is
if(use_gbackflow)then
call clear_scratch_gbf(is)
else
call clear_scratch_pbf(is)
endif
end subroutine clear_scratch_bf
subroutine clone_scratch_bf(is,js)
implicit none
integer,intent(in) :: is,js
if(use_gbackflow)then
call clone_scratch_gbf(is,js)
else
call clone_scratch_pbf(is,js)
endif
end subroutine clone_scratch_bf
subroutine init_backflow(empty_backflow)
implicit none
logical,intent(inout) :: empty_backflow
if(use_gbackflow)then
call init_gbackflow(empty_backflow)
else
call init_pbackflow(empty_backflow)
endif
end subroutine init_backflow
subroutine write_backflow(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
if(use_gbackflow)then
call write_gbackflow(correlation_name)
else
call write_pbackflow(correlation_name)
endif
end subroutine write_backflow
subroutine setup_bf_params(nparam)
implicit none
integer,intent(out) :: nparam
if(use_gbackflow)then
call setup_gbf_params(nparam)
else
call setup_pbf_params(nparam)
endif
end subroutine setup_bf_params
subroutine finish_bf_params
implicit none
if(use_gbackflow)then
call finish_gbf_params
else
call finish_pbf_params
endif
end subroutine finish_bf_params
subroutine get_bf_params(params,has_lolim,lolim,has_hilim,hilim,is_sha&
&llow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,label)
implicit none
real(dp),intent(inout) :: params(:),lolim(:),hilim(:)
logical,intent(inout) :: has_lolim(:),has_hilim(:),is_shallow(:),is_re&
&dundant(:),is_linear(:),is_loglinear(:),has_aderiv(:),affect_map(:,:)
character(2),intent(inout) :: label(:)
if(use_gbackflow)then
call get_gbf_params(params,has_lolim,lolim,has_hilim,hilim,is_shallow,&
&is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,label)
else
call get_pbf_params(params,has_lolim,lolim,has_hilim,hilim,is_shallow,&
&is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,label)
endif
end subroutine get_bf_params
subroutine put_bf_params(params,ignore,iparam_buffer,prestore,bad_para&
&ms)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(:)
logical,intent(in) :: ignore(:),prestore
logical,intent(out) :: bad_params
if(use_gbackflow)then
call put_gbf_params(params,ignore,iparam_buffer,prestore,bad_params)
else
call put_pbf_params(params,ignore,iparam_buffer,prestore,bad_params)
endif
end subroutine put_bf_params
subroutine bf_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,verbose&
&)
implicit none
integer,intent(out) :: ie,jspin
real(dp),intent(in) :: eevecs(*),eivecs(*)
logical,intent(in) :: verbose
logical,intent(out) :: fail
if(use_gbackflow)then
call gbf_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,verbose)
else
call pbf_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,verbose)
endif
end subroutine bf_assess_check_kinetic
subroutine setup_backflow_plot(makeplot,phiplot,ii,z,rii)
implicit none
integer,intent(in),optional :: ii
real(dp),intent(in),optional :: z,rii
logical,intent(in) :: makeplot
logical,intent(in),optional :: phiplot
if(use_gbackflow)then
if(present(phiplot))then
call setup_gbackflow_plot(makeplot,phiplot,ii,z,rii)
else
call setup_gbackflow_plot(makeplot)
endif
else
if(present(phiplot))then
call setup_pbackflow_plot(makeplot,phiplot,ii,z,rii)
else
call setup_pbackflow_plot(makeplot)
endif
endif
end subroutine setup_backflow_plot
subroutine plot_backflow(rele)
implicit none
real(dp),intent(inout) :: rele(*)
if(use_gbackflow)then
call plot_gbackflow(rele)
else
call plot_pbackflow(rele)
endif
end subroutine plot_backflow
subroutine backflow_stats(in_range,in_range1)
implicit none
real(dp),intent(out) :: in_range,in_range1
if(use_gbackflow)then
call gbackflow_stats(in_range,in_range1)
else
call pbackflow_stats(in_range,in_range1)
endif
end subroutine backflow_stats
subroutine enumerate_plot_bf(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
if(use_gbackflow)then
call enumerate_plot_gbf(n,keyword,description)
else
call enumerate_plot_pbf(n,keyword,description)
endif
end subroutine enumerate_plot_bf
subroutine query_plot_bf(iplot,ii,rank,is_complex,has_stderr,rot_tenso&
&r,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
if(use_gbackflow)then
call query_plot_gbf(iplot,ii,rank,is_complex,has_stderr,rot_tensor,tra&
&nsl_pos,nfunctions,function_name)
else
call query_plot_pbf(iplot,ii,rank,is_complex,has_stderr,rot_tensor,tra&
&nsl_pos,nfunctions,function_name)
endif
end subroutine query_plot_bf
subroutine get_plot_bf(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
if(use_gbackflow)then
call get_plot_gbf(iplot,ii,is0,is1,f)
else
call get_plot_pbf(iplot,ii,is0,is1,f)
endif
end subroutine get_plot_bf
subroutine finish_plot_bf
implicit none
if(use_gbackflow)then
call finish_plot_gbf
else
call finish_plot_pbf
endif
end subroutine finish_plot_bf
end module slaarnaad

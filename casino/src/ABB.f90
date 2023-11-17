module slaarnabb
use casl
use dsp
use slaarnaag,only : czero
use slaarnabg, only : nitot
use store,    only : netot,real1_complex2
implicit none
private
public init_gbackflow,write_gbackflow,plot_gbackflow,gbf_assess_check_&
&kinetic,gbackflow_stats,setup_gbackflow_plot,enumerate_plot_gbf,query&
&_plot_gbf,get_plot_gbf,finish_plot_gbf
public setup_gbf,finish_gbf,accept_move_gbf,reset_config_gbf,get_gbf_x&
&,get_eevecs_gbf,loggrad_gbf,loglap_gbf,setup_gbf_params,finish_gbf_pa&
&rams,get_gbf_params,put_gbf_params,invalidate_param1_gbf,clear_scratc&
&h_gbf,clone_scratch_gbf
contains
subroutine setup_gbf
implicit none
end subroutine setup_gbf
subroutine finish_gbf
implicit none
end subroutine finish_gbf
subroutine accept_move_gbf(is,js)
implicit none
integer,intent(in) :: is,js
end subroutine accept_move_gbf
subroutine reset_config_gbf(is,js)
implicit none
integer,intent(in) :: is,js
call clear_scratch_gbf(js)
end subroutine reset_config_gbf
subroutine loggrad_gbf(ii,is,sd,farray,loggrad_psi)
implicit none
integer,intent(in) :: ii,is
real(dp),intent(in) :: farray(3,real1_complex2,netot)
complex(dp),intent(out) :: loggrad_psi(3)
logical,intent(in) :: sd
loggrad_psi=czero
end subroutine loggrad_gbf
subroutine loglap_gbf(ii,is,farray,harray,loggrad_psi,loglap_psi)
implicit none
integer,intent(in) :: ii,is
real(dp),intent(in) :: farray(3,real1_complex2,netot),harray(3,3,real1&
&_complex2,netot,netot)
complex(dp),intent(in) :: loggrad_psi(3)
complex(dp),intent(out) :: loglap_psi
loglap_psi=czero
end subroutine loglap_gbf
subroutine get_gbf_x(is,val,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: val,fd,sd
end subroutine get_gbf_x
subroutine get_eevecs_gbf(is)
implicit none
integer,intent(in) :: is
end subroutine get_eevecs_gbf
subroutine invalidate_param1_gbf(is,iparam)
implicit none
integer,intent(in) :: is,iparam
end subroutine invalidate_param1_gbf
subroutine clear_scratch_gbf(is)
implicit none
integer,intent(in) :: is
end subroutine clear_scratch_gbf
subroutine clone_scratch_gbf(is,js)
implicit none
integer,intent(in) :: is,js
end subroutine clone_scratch_gbf
subroutine init_gbackflow(empty_backflow)
implicit none
logical,intent(inout) :: empty_backflow
end subroutine init_gbackflow
subroutine write_gbackflow(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
end subroutine write_gbackflow
subroutine setup_gbf_params(nparam)
implicit none
integer,intent(out) :: nparam
nparam=0
end subroutine setup_gbf_params
subroutine finish_gbf_params
implicit none
end subroutine finish_gbf_params
subroutine get_gbf_params(params,has_lolim,lolim,has_hilim,hilim,is_sh&
&allow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,label&
&)
implicit none
real(dp),intent(inout) :: params(:),lolim(:),hilim(:)
logical,intent(inout) :: has_lolim(:),has_hilim(:),is_shallow(:),is_re&
&dundant(:),is_linear(:),is_loglinear(:),has_aderiv(:),affect_map(:,:)
character(2),intent(inout) :: label(:)
end subroutine get_gbf_params
subroutine put_gbf_params(params,ignore,iparam_buffer,prestore,bad_par&
&ams)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(:)
logical,intent(in) :: ignore(:),prestore
logical,intent(out) :: bad_params
bad_params=.false.
end subroutine put_gbf_params
subroutine gbf_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,verbos&
&e)
implicit none
logical,intent(in) :: verbose
real(dp),intent(in) :: eevecs(4,netot,netot),eivecs(4,nitot,netot)
integer,intent(out) :: ie,jspin
logical,intent(out) :: fail
ie=1 
jspin=1 
fail=.false.
end subroutine gbf_assess_check_kinetic
subroutine setup_gbackflow_plot(makeplot,phiplot,ii,z,rii)
implicit none
integer,intent(in),optional :: ii
real(dp),intent(in),optional :: z,rii
logical,intent(in) :: makeplot
logical,intent(in),optional :: phiplot
end subroutine setup_gbackflow_plot
subroutine plot_gbackflow(rele)
implicit none
real(dp),intent(inout) :: rele(3,netot)
end subroutine plot_gbackflow
subroutine gbackflow_stats(in_range,in_range1)
implicit none
real(dp),intent(out) :: in_range,in_range1
in_range=1.d0/dble(netot) 
in_range1=in_range
end subroutine gbackflow_stats
subroutine enumerate_plot_gbf(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
n=0
end subroutine enumerate_plot_gbf
subroutine query_plot_gbf(iplot,ii,rank,is_complex,has_stderr,rot_tens&
&or,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
nfunctions=0
end subroutine query_plot_gbf
subroutine get_plot_gbf(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
f(1)=0.d0
end subroutine get_plot_gbf
subroutine finish_plot_gbf
implicit none
end subroutine finish_plot_gbf
end module slaarnabb

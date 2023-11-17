module slaarnabv
use dsp
use slaarnach
use store
use slaarnaag,    only : pi
use slaarnaan,only : ee_distances_all
use file_utils,   only : open_units
use format_utils, only : wout,i2s,d2s,r2s,write_list_int,byte2human
use slaarnabg,     only : dimensionality,wigner_seitz_radius,isperiodi&
&c,nitot,rion,nbasis,nitype
use slaarnabt,    only : dcopy,ddot,dcopy_ee,reduced_echelon
use parallel,     only : am_master
use slaarnaca,        only : is_ae
use run_control,  only : errstop,errwarn,timer,check_alloc,errstop_mas&
&ter
implicit none
private
public init_pbackflow,write_pbackflow,plot_pbackflow,pbf_assess_check_&
&kinetic,pbackflow_stats,setup_pbackflow_plot,enumerate_plot_pbf,query&
&_plot_pbf,get_plot_pbf,finish_plot_pbf
public setup_pbf,finish_pbf,get_pbf_x,get_eevecs_pbf,loggrad_pbf,logla&
&p_pbf,accept_move_pbf,reset_config_pbf,setup_pbf_params,finish_pbf_pa&
&rams,get_pbf_params,put_pbf_params,invalidate_param1_pbf,clear_scratc&
&h_pbf,clone_scratch_pbf
logical xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1
integer xyzzyaaad1,xyzzyaaae1,xyzzyaaaf1,xyzzyaaag1,xyzzyaaah1,xyzzyaa&
&ai1,xyzzyaaaj1,xyzzyaaak1,xyzzyaaal1,xyzzyaaam1
integer,allocatable :: xyzzyaaan1(:),xyzzyaaao1(:),xyzzyaaap1(:),xyzzy&
&aaaq1(:),xyzzyaaar1(:,:),xyzzyaaas1(:),xyzzyaaat1(:),xyzzyaaau1(:),xy&
&zzyaaav1(:),xyzzyaaaw1(:),xyzzyaaax1(:,:),xyzzyaaay1(:),xyzzyaaaz1(:)&
&,xyzzyaaba1(:)
integer :: xyzzyaabb1,xyzzyaabc1
integer,allocatable :: xyzzyaabd1(:,:),xyzzyaabe1(:,:)
integer,allocatable :: xyzzyaabf1(:),xyzzyaabg1(:)
real(dp),allocatable :: xyzzyaabh1(:),xyzzyaabi1(:,:),xyzzyaabj1(:),xy&
&zzyaabk1(:,:,:),xyzzyaabl1(:),xyzzyaabm1(:,:,:,:,:),xyzzyaabn1(:,:,:,&
&:,:)
integer,allocatable :: xyzzyaabo1(:),xyzzyaabp1(:,:),xyzzyaabq1(:),xyz&
&zyaabr1(:,:,:),xyzzyaabs1(:),xyzzyaabt1(:,:,:,:,:),xyzzyaabu1(:,:,:,:&
&,:)
integer xyzzyaabv1,xyzzyaabw1,xyzzyaabx1,xyzzyaaby1
real(dp),allocatable :: xyzzyaabz1(:),xyzzyaaca1(:),xyzzyaacb1(:),xyzz&
&yaacc1(:),xyzzyaacd1(:),xyzzyaace1(:),xyzzyaacf1(:),xyzzyaacg1(:),xyz&
&zyaach1(:),xyzzyaaci1(:)
real(dp),allocatable :: xyzzyaacj1(:,:),xyzzyaack1(:,:),xyzzyaacl1(:,:&
&),xyzzyaacm1(:,:,:),xyzzyaacn1(:,:,:),xyzzyaaco1(:,:,:)
real(dp) xyzzyaacp1,xyzzyaacq1
logical,allocatable :: xyzzyaacr1(:),xyzzyaacs1(:,:),xyzzyaact1(:),xyz&
&zyaacu1(:),xyzzyaacv1(:)
integer,allocatable :: xyzzyaacw1(:,:),xyzzyaacx1(:,:),xyzzyaacy1(:,:,&
&:)
real(dp),allocatable :: xyzzyaacz1(:),xyzzyaada1(:,:,:,:)
logical xyzzyaadb1
integer xyzzyaadc1
integer,allocatable :: xyzzyaadd1(:),xyzzyaade1(:)
real(dp),allocatable :: xyzzyaadf1(:),xyzzyaadg1(:),xyzzyaadh1(:),xyzz&
&yaadi1(:),xyzzyaadj1(:),xyzzyaadk1(:),xyzzyaadl1(:)
integer,allocatable :: xyzzyaadm1(:,:),xyzzyaadn1(:,:),xyzzyaado1(:,:,&
&:)
integer xyzzyaadp1
character(80) title
real(dp) xyzzyaadq1,xyzzyaadr1
logical xyzzyaads1,xyzzyaadt1
integer xyzzyaadu1
integer,parameter :: xyzzyaadv1=4
integer,allocatable :: xyzzyaadw1(:)
integer xyzzyaadx1,xyzzyaady1(xyzzyaadv1)
real(dp),allocatable :: xyzzyaadz1(:,:),xyzzyaaea1(:,:,:),xyzzyaaeb1(:&
&,:),xyzzyaaec1(:,:),xyzzyaaed1(:,:)
integer xyzzyaaee1
real(dp),allocatable :: xyzzyaaef1(:,:),xyzzyaaeg1(:,:,:,:),xyzzyaaeh1&
&(:,:),xyzzyaaei1(:,:),xyzzyaaej1(:,:)
integer xyzzyaaek1
real(dp),allocatable :: xyzzyaael1(:,:),xyzzyaaem1(:,:),xyzzyaaen1(:,:&
&),xyzzyaaeo1(:,:),xyzzyaaep1(:,:,:,:,:,:),xyzzyaaeq1(:,:,:,:,:,:),xyz&
&zyaaer1(:)
integer xyzzyaaes1
real(dp),allocatable :: xyzzyaaet1(:,:),xyzzyaaeu1(:,:),xyzzyaaev1(:,:&
&),xyzzyaaew1(:,:),xyzzyaaex1(:,:)
integer,parameter :: xyzzyaaey1=2,xyzzyaaez1=1,xyzzyaafa1=2
integer :: xyzzyaafb1(xyzzyaaey1)=0
character(64),parameter :: xyzzyaafc1(xyzzyaaey1)=(/'bfdisp  ','bfconf&
&ig'/)
character(64),parameter :: xyzzyaafd1(xyzzyaaey1)=(/'backflow displace&
&ment on selected particle',  'backflow configuration                 &
&   '/)
contains
subroutine setup_pbf
implicit none
integer xyzzyaaaa2(2),xyzzyaaab2(2),xyzzyaaac2(2)
integer xyzzyaaad2,xyzzyaaae2
xyzzyaabv1=0
xyzzyaabw1=0
xyzzyaabx1=0
xyzzyaaby1=0
xyzzyaaaa2=0
xyzzyaaab2=0
xyzzyaaac2=0
call include_range((/1,nscratch/),xyzzyaaaa2)
if(pairing_wf)call include_range((/1,nscratch/),xyzzyaaab2)
call include_range(ratio1_to_sz,xyzzyaaac2)
if(xyzzyaaaa2(1)/=0)then
allocate(bf_x_scr(3,netot,xyzzyaaaa2(1):xyzzyaaaa2(2)),bf_connect_scr(&
&netot,netot,xyzzyaaaa2(1):xyzzyaaaa2(2)),bf_m_scr(netot,xyzzyaaaa2(1)&
&:xyzzyaaaa2(2)),bf_rmap_scr(netot,netot,xyzzyaaaa2(1):xyzzyaaaa2(2)),&
&bf_m2_scr(netot,xyzzyaaaa2(1):xyzzyaaaa2(2)),bf_rmap2_scr(netot,netot&
&,xyzzyaaaa2(1):xyzzyaaaa2(2)),bf_dx_scr(3,3,netot,netot,xyzzyaaaa2(1)&
&:xyzzyaaaa2(2)),bf_d2x_scr(3,netot,netot,xyzzyaaaa2(1):xyzzyaaaa2(2))&
&,stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'SETUP_BF','bf_x')
bf_x_scr=0.d0
bf_connect_scr=.true.
bf_m_scr=0
bf_rmap_scr=0
bf_dx_scr=0.d0
bf_d2x_scr=0.d0
bf_m2_scr=0
bf_rmap2_scr=0
endif
allocate(bf_x_valid(nscratch),bf_dx_valid(nscratch),bf_d2x_valid(nscra&
&tch),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'SETUP_BF','bf_x_valid')
if(xyzzyaaab2(1)/=0)then
allocate(eevecs_bf_scr(4,netot,netot,xyzzyaaab2(1):xyzzyaaab2(2)),stat&
&=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'SETUP_BF','eevecs_bf')
eevecs_bf_scr=0.d0
endif
allocate(eevecs_bf_valid(nscratch),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'SETUP_BF','eevecs_bf_valid')
if(xyzzyaaac2(1)/=0)then
allocate(bf_m_chscr(nspin,xyzzyaaac2(1):xyzzyaaac2(2)),bf_rmap_chscr(n&
&emax,nspin,xyzzyaaac2(1):xyzzyaaac2(2)),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'SETUP_BF','bf_change')
bf_m_chscr=0
bf_rmap_chscr=0
endif
allocate(bf_rmap_chvalid(nscratch),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'SETUP_BF','bf_rmap_chvalid')
do xyzzyaaae2=1,nscratch
call clear_scratch_pbf(xyzzyaaae2)
enddo
end subroutine setup_pbf
subroutine finish_pbf
implicit none
if(allocated(bf_x_scr))deallocate(bf_x_scr,bf_connect_scr,bf_m_scr,bf_&
&rmap_scr,bf_dx_scr,bf_d2x_scr,bf_m2_scr,bf_rmap2_scr)
deallocate(bf_x_valid,bf_dx_valid,bf_d2x_valid)
if(allocated(eevecs_bf_scr))deallocate(eevecs_bf_scr)
deallocate(eevecs_bf_valid)
if(allocated(bf_m_chscr))deallocate(bf_m_chscr,bf_rmap_chscr)
deallocate(bf_rmap_chvalid)
end subroutine finish_pbf
subroutine accept_move_pbf(is,js)
implicit none
integer,intent(in) :: is,js
if(bf_x_valid(js))then
call dcopy(three_netot,bf_x_scr(1,1,js),1,bf_x_scr(1,1,is),1)
bf_connect_scr(:,:,is)=bf_connect_scr(:,:,js)
bf_m_scr(:,is)=bf_m_scr(:,js)
bf_rmap_scr(:,:,is)=bf_rmap_scr(:,:,js)
bf_x_valid(is)=.true.
else
bf_x_valid(is)=.false.
endif
if(bf_dx_valid(js))then
call dcopy(nine_netot_netot,bf_dx_scr(1,1,1,1,js),1,bf_dx_scr(1,1,1,1,&
&is),1)
bf_dx_valid(is)=.true.
else
bf_dx_valid(is)=.false.
endif
if(bf_d2x_valid(js))then
call dcopy(three_netot_netot,bf_d2x_scr(1,1,1,js),1,bf_d2x_scr(1,1,1,i&
&s),1)
bf_m2_scr(:,is)=bf_m2_scr(:,js)
bf_rmap2_scr(:,:,is)=bf_rmap2_scr(:,:,js)
bf_d2x_valid(is)=.true.
else
bf_d2x_valid(is)=.false.
endif
if(eevecs_bf_valid(js))then
call dcopy(four_netot_netot,eevecs_bf_scr(1,1,1,js),1,eevecs_bf_scr(1,&
&1,1,is),1)
eevecs_bf_valid(is)=.true.
else
eevecs_bf_valid(is)=.false.
endif
bf_rmap_chvalid(is)=.false.
end subroutine accept_move_pbf
subroutine reset_config_pbf(is,js)
implicit none
integer,intent(in) :: is,js
call clear_scratch_pbf(js)
end subroutine reset_config_pbf
subroutine loggrad_pbf(ii,is,sd,farray,loggrad_psi)
implicit none
integer,intent(in) :: ii,is
logical,intent(in) :: sd
real(dp),intent(in) :: farray(3,real1_complex2,netot)
complex(dp),intent(out) :: loggrad_psi(3)
call get_pbf_x(is,.true.,.true.,sd)
call xyzzyaagt1(farray,bf_dx_scr(1,1,1,ii,is),bf_m_scr(ii,is),bf_rmap_&
&scr(1,ii,is),loggrad_psi)
end subroutine loggrad_pbf
subroutine loglap_pbf(ii,is,farray,harray,loggrad_psi,loglap_psi)
implicit none
complex(dp),intent(in) :: loggrad_psi(3)
integer,intent(in) :: ii,is
real(dp),intent(in) :: farray(3,real1_complex2,netot),harray(3,3,real1&
&_complex2,netot,netot)
complex(dp),intent(out) :: loglap_psi
call get_pbf_x(is,.true.,.true.,.true.)
call xyzzyaagu1(farray,harray,bf_dx_scr(1,1,1,ii,is),bf_d2x_scr(1,1,ii&
&,is),bf_m_scr(ii,is),bf_rmap_scr(1,ii,is),loggrad_psi,loglap_psi)
end subroutine loglap_pbf
recursive subroutine get_pbf_x(xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzya&
&aad8)
implicit none
integer,intent(in) :: xyzzyaaaa8
logical,intent(in) :: xyzzyaaab8,xyzzyaaac8,xyzzyaaad8
integer xyzzyaaae8,xyzzyaaaf8
logical xyzzyaaag8,xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8
xyzzyaaag8=xyzzyaaab8.and..not.bf_x_valid(xyzzyaaaa8)
xyzzyaaah8=xyzzyaaac8.and..not.bf_dx_valid(xyzzyaaaa8)
xyzzyaaai8=xyzzyaaad8.and..not.bf_d2x_valid(xyzzyaaaa8)
xyzzyaaaj8=xyzzyaaag8.or.xyzzyaaah8.or.xyzzyaaai8
if(.not.xyzzyaaaj8)return
xyzzyaaae8=buffer_move1_from(xyzzyaaaa8)
if(xyzzyaaae8/=0.and..not.xyzzyaaai8)then
call get_pbf_x(xyzzyaaae8,xyzzyaaag8,xyzzyaaah8,.false.)
xyzzyaaaf8=buffer_move1_from_ii(xyzzyaaaa8)
bf_connect_scr(:,:,xyzzyaaaa8)=bf_connect_scr(:,:,xyzzyaaae8)
bf_m_scr(:,xyzzyaaaa8)=bf_m_scr(:,xyzzyaaae8)
bf_rmap_scr(:,:,xyzzyaaaa8)=bf_rmap_scr(:,:,xyzzyaaae8)
call dcopy(three_netot,bf_x_scr(1,1,xyzzyaaae8),1,bf_x_scr(1,1,xyzzyaa&
&aa8),1)
if(xyzzyaaah8)call dcopy(nine_netot_netot,bf_dx_scr(1,1,1,1,xyzzyaaae8&
&),1,bf_dx_scr(1,1,1,1,xyzzyaaaa8),1)
call get_rsele(xyzzyaaae8)
call get_eevecs(xyzzyaaae8)
call get_eivecs(xyzzyaaae8)
call get_eevecs1_ch(xyzzyaaaf8,xyzzyaaaa8)
call get_eivecs(xyzzyaaaa8)
if(xyzzyaaah8)then
call xyzzyaafr1(xyzzyaaaf8,rele1_chscr(1,xyzzyaaaa8),eevecs1_chscr(1,1&
&,xyzzyaaaa8),eivecs_scr(1,1,1,xyzzyaaaa8),rele_scr(1,xyzzyaaaf8,xyzzy&
&aaae8),eevecs_scr(1,1,xyzzyaaaf8,xyzzyaaae8),eivecs_scr(1,1,1,xyzzyaa&
&ae8),bf_connect_scr(1,1,xyzzyaaaa8),bf_m_scr(1,xyzzyaaaa8),bf_rmap_sc&
&r(1,1,xyzzyaaaa8),bf_m_chscr(1,xyzzyaaaa8),bf_rmap_chscr(1,1,xyzzyaaa&
&a8),bf_x_scr(1,1,xyzzyaaaa8),bf_dx=bf_dx_scr(1,1,1,1,xyzzyaaaa8))
bf_x_valid(xyzzyaaaa8)=.true.
bf_dx_valid(xyzzyaaaa8)=.true.
else
call xyzzyaafr1(xyzzyaaaf8,rele1_chscr(1,xyzzyaaaa8),eevecs1_chscr(1,1&
&,xyzzyaaaa8),eivecs_scr(1,1,1,xyzzyaaaa8),rele_scr(1,xyzzyaaaf8,xyzzy&
&aaae8),eevecs_scr(1,1,xyzzyaaaf8,xyzzyaaae8),eivecs_scr(1,1,1,xyzzyaa&
&ae8),bf_connect_scr(1,1,xyzzyaaaa8),bf_m_scr(1,xyzzyaaaa8),bf_rmap_sc&
&r(1,1,xyzzyaaaa8),bf_m_chscr(1,xyzzyaaaa8),bf_rmap_chscr(1,1,xyzzyaaa&
&a8),bf_x_scr(1,1,xyzzyaaaa8))
bf_x_valid(xyzzyaaaa8)=.true.
endif
bf_rmap_chvalid(xyzzyaaaa8)=.true.
else
call get_rsele(xyzzyaaaa8)
call get_eevecs(xyzzyaaaa8)
call get_eivecs(xyzzyaaaa8)
if(xyzzyaaah8.or.xyzzyaaai8)then
call xyzzyaafv1(rele_scr(1,1,xyzzyaaaa8),eevecs_scr(1,1,1,xyzzyaaaa8),&
&eivecs_scr(1,1,1,xyzzyaaaa8),bf_connect_scr(1,1,xyzzyaaaa8),bf_m_scr(&
&1,xyzzyaaaa8),bf_rmap_scr(1,1,xyzzyaaaa8),bf_x_scr(1,1,xyzzyaaaa8),bf&
&_dx=bf_dx_scr(1,1,1,1,xyzzyaaaa8),bf_d2x=bf_d2x_scr(1,1,1,xyzzyaaaa8)&
&,bf_m2=bf_m2_scr(1,xyzzyaaaa8),bf_rmap2=bf_rmap2_scr(1,1,xyzzyaaaa8))
bf_x_valid(xyzzyaaaa8)=.true.
bf_dx_valid(xyzzyaaaa8)=.true.
bf_d2x_valid(xyzzyaaaa8)=.true.
else
call xyzzyaafv1(rele_scr(1,1,xyzzyaaaa8),eevecs_scr(1,1,1,xyzzyaaaa8),&
&eivecs_scr(1,1,1,xyzzyaaaa8),bf_connect_scr(1,1,xyzzyaaaa8),bf_m_scr(&
&1,xyzzyaaaa8),bf_rmap_scr(1,1,xyzzyaaaa8),bf_x_scr(1,1,xyzzyaaaa8))
bf_x_valid(xyzzyaaaa8)=.true.
endif
endif
end subroutine get_pbf_x
subroutine get_eevecs_pbf(is)
implicit none
integer,intent(in) :: is
if(eevecs_bf_valid(is))return
call get_pbf_x(is,.true.,.false.,.false.)
call ee_distances_all(netot,bf_x_scr(1,1,is),eevecs_bf_scr(1,1,1,is))
eevecs_bf_valid(is)=.true.
end subroutine get_eevecs_pbf
subroutine invalidate_param1_pbf(is,iparam)
implicit none
integer,intent(in) :: is,iparam
bf_x_valid(is)=.false.
bf_dx_valid(is)=.false.
bf_d2x_valid(is)=.false.
bf_rmap_chvalid(is)=.false.
eevecs_bf_valid(is)=.false.
end subroutine invalidate_param1_pbf
subroutine clear_scratch_pbf(is)
implicit none
integer,intent(in) :: is
bf_x_valid(is)=.false.
bf_dx_valid(is)=.false.
bf_d2x_valid(is)=.false.
bf_rmap_chvalid(is)=.false.
eevecs_bf_valid(is)=.false.
end subroutine clear_scratch_pbf
subroutine clone_scratch_pbf(is,js)
implicit none
integer,intent(in) :: is,js
if(bf_x_valid(is).and..not.bf_x_valid(js))then
call dcopy(three_netot,bf_x_scr(1,1,is),1,bf_x_scr(1,1,js),1)
bf_connect_scr(:,:,js)=bf_connect_scr(:,:,is)
bf_m_scr(:,js)=bf_m_scr(:,is)
bf_rmap_scr(:,:,js)=bf_rmap_scr(:,:,is)
bf_x_valid(js)=.true.
endif
if(bf_dx_valid(is).and..not.bf_dx_valid(js))then
call dcopy(nine_netot_netot,bf_dx_scr(1,1,1,1,is),1,bf_dx_scr(1,1,1,1,&
&js),1)
bf_dx_valid(js)=.true.
endif
if(bf_d2x_valid(is).and..not.bf_d2x_valid(js))then
call dcopy(three_netot_netot,bf_d2x_scr(1,1,1,is),1,bf_d2x_scr(1,1,1,j&
&s),1)
bf_m2_scr(:,js)=bf_m2_scr(:,is)
bf_rmap2_scr(:,:,js)=bf_rmap2_scr(:,:,is)
bf_d2x_valid(js)=.true.
endif
if(eevecs_bf_valid(is).and..not.eevecs_bf_valid(js))then
call dcopy(four_netot_netot,eevecs_bf_scr(1,1,1,is),1,eevecs_bf_scr(1,&
&1,1,js),1)
eevecs_bf_valid(js)=.true.
endif
end subroutine clone_scratch_pbf
subroutine enumerate_plot_pbf(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
integer xyzzyaaaa13,xyzzyaaab13
logical xyzzyaaac13,xyzzyaaad13
xyzzyaaad13=.not.(present(keyword).and.present(description))
xyzzyaaab13=0
xyzzyaafb1=0
do xyzzyaaaa13=1,xyzzyaaey1
xyzzyaaac13=.false.
select case(xyzzyaaaa13)
case(xyzzyaaez1)
xyzzyaaac13=.true.
case(xyzzyaafa1)
xyzzyaaac13=.true.
end select
if(xyzzyaaac13)then
xyzzyaaab13=xyzzyaaab13+1
xyzzyaafb1(xyzzyaaab13)=xyzzyaaaa13
endif
enddo
n=xyzzyaaab13
if(.not.xyzzyaaad13)then
do xyzzyaaaa13=1,xyzzyaaab13
keyword(xyzzyaaaa13)=xyzzyaafc1(xyzzyaafb1(xyzzyaaaa13))
description(xyzzyaaaa13)=xyzzyaafd1(xyzzyaafb1(xyzzyaaaa13))
enddo
endif
end subroutine enumerate_plot_pbf
subroutine query_plot_pbf(iplot,ii,rank,is_complex,has_stderr,rot_tens&
&or,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
integer jj
logical count_only
count_only=.not.present(function_name)
if(iplot<0.or.iplot>xyzzyaaey1)call errstop_master('QUERY_PLOT_BF','IP&
&LOT out of range. Bug in calling routine.')
if(xyzzyaafb1(iplot)==0)call errstop_master('QUERY_PLOT_BF','IPLOT out&
& of range. Bug in calling routine.')
select case(xyzzyaafb1(iplot))
case(xyzzyaaez1)
rank=1
is_complex=.false.
has_stderr=.false.
rot_tensor=.true.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Backflow displacement on particle '//trim(i&
&2s(ii))
endif
case(xyzzyaafa1)
rank=1
is_complex=.false.
has_stderr=.false.
rot_tensor=.true.
transl_pos=.true.
nfunctions=0
do jj=1,netot
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Backflow coordinate of particle '//trim(i2s&
&(jj))
endif
enddo
end select
end subroutine query_plot_pbf
subroutine get_plot_pbf(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
integer xyzzyaaaa15,xyzzyaaab15
if(iplot<0.or.iplot>xyzzyaaey1)call errstop_master('GET_PLOT_BF','IPLO&
&T out of range. Bug in calling routine.')
if(xyzzyaafb1(iplot)==0)call errstop_master('GET_PLOT_BF','IPLOT out o&
&f range. Bug in calling routine.')
select case(xyzzyaafb1(iplot))
case(xyzzyaaez1)
call get_rsele(is1)
call get_pbf_x(is1,.true.,.false.,.false.)
f(1:dimensionality)=bf_x_scr(1:dimensionality,ii,is1)-rele_scr(1:dimen&
&sionality,ii,is1)
case(xyzzyaafa1)
call get_pbf_x(is1,.true.,.false.,.false.)
xyzzyaaaa15=0
do xyzzyaaab15=1,netot
f(xyzzyaaaa15+1:xyzzyaaaa15+dimensionality)=bf_x_scr(1:dimensionality,&
&xyzzyaaab15,is1)
xyzzyaaaa15=xyzzyaaaa15+dimensionality
enddo
end select
end subroutine get_plot_pbf
subroutine finish_plot_pbf
implicit none
xyzzyaafb1=0
end subroutine finish_plot_pbf
subroutine init_pbackflow(empty_backflow)
use slaarnaan,only : atoms_label_pcell,atoms_label_species
use slaarnabg, only : rion
implicit none
logical,intent(inout) :: empty_backflow
integer xyzzyaaaa17,xyzzyaaab17,s,xyzzyaaac17,xyzzyaaad17,xyzzyaaae17,&
&xyzzyaaaf17,xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,xyzzyaaaj17,xyzzyaaak&
&17,xyzzyaaal17,xyzzyaaam17(5),xyzzyaaan17(3),xyzzyaaao17,xyzzyaaap17,&
&xyzzyaaaq17,xyzzyaaar17,xyzzyaaas17,xyzzyaaat17
integer,allocatable :: xyzzyaaau17(:,:,:),xyzzyaaav17(:,:,:,:,:),xyzzy&
&aaaw17(:,:,:,:,:),xyzzyaaax17(:)
real(dp) xyzzyaaay17,xyzzyaaaz17,xyzzyaaba17,xyzzyaabb17
real(dp),allocatable :: xyzzyaabc17(:,:,:),xyzzyaabd17(:,:,:,:,:),xyzz&
&yaabe17(:,:,:,:,:)
logical xyzzyaabf17,xyzzyaabg17,xyzzyaabh17,xyzzyaabi17,xyzzyaabj17
logical,allocatable :: xyzzyaabk17(:,:),xyzzyaabl17(:,:)
character(80) lineread,tmpr,tmpr2,tmps,term_label
if(.not.use_backflow)return
xyzzyaabb17=0.d0
xyzzyaaba17=0.d0
xyzzyaaaa1=.false.
xyzzyaaac1=.false.
xyzzyaaab1=.false.
xyzzyaadb1=.false.
xyzzyaabi17=.false.
xyzzyaabh17=.false.
empty_backflow=.true.
if(am_master)then
call wout('Backflow setup')
call wout('==============')
endif
if(any(ee_cusp_in_orbital))call errwarn('INIT_BACKFLOW','Backflow shou&
&ldn''t be used when any particle-particle cusp conditions are applied&
& on the (pairing) orbitals. Using backflow will break the cusp condit&
&ions and the local energy will fluctuate wildly. Try applying the cus&
&p conditions on the Jastrow instead.')
if(nitot>0.and.allocated(is_ae))then
xyzzyaadb1=any(is_ae)
if(xyzzyaadb1)then
allocate(xyzzyaadk1(nitot),xyzzyaadf1(nitot),xyzzyaaax17(nitot),xyzzya&
&adg1(nitot),xyzzyaadh1(nitot),xyzzyaadi1(nitot),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','AE-cutoff arrays')
xyzzyaadk1(:)=0.d0
xyzzyaadf1(:)=0.d0
xyzzyaaax17(:)=1
do xyzzyaaae17=1,nitot
if(.not.is_ae(xyzzyaaae17))cycle
xyzzyaaaz17=0.d0
do xyzzyaaaf17=1,nitot
if(xyzzyaaae17==xyzzyaaaf17)cycle
if(is_ae(xyzzyaaaf17))then
xyzzyaaay17=.5d0*sqrt(sum((rion(1:3,xyzzyaaae17)-rion(1:3,xyzzyaaaf17)&
&)**2))
if(xyzzyaaaz17==0.d0)then
xyzzyaaaz17=xyzzyaaay17
else
xyzzyaaaz17=min(xyzzyaaaz17,xyzzyaaay17)
endif
endif
enddo
if(isperiodic)then
if(xyzzyaaaz17==0.d0)then
xyzzyaadk1(xyzzyaaae17)=0.999999d0*wigner_seitz_radius
else
xyzzyaadk1(xyzzyaaae17)=min(xyzzyaaaz17,0.999999d0*wigner_seitz_radius&
&)
endif
else
xyzzyaadk1(xyzzyaaae17)=xyzzyaaaz17
endif
enddo
endif
endif
inquire(file='correlation.data',exist=xyzzyaabf17)
if(am_master.and..not.xyzzyaabf17)call errstop('INIT_BACKFLOW','Cannot&
& find correlation.data file.')
call open_units(xyzzyaadp1,xyzzyaaab17)
open(unit=xyzzyaadp1,file='correlation.data',status='old',iostat=xyzzy&
&aaab17)
if(am_master.and.xyzzyaaab17/=0)call errstop('INIT_BACKFLOW','Problem &
&opening correlation.data .')
if(am_master)then
call wout('Reading correlation.data file.')
call wout()
endif
do
read(xyzzyaadp1,'(a)',iostat=xyzzyaaab17)lineread
if(trim(adjustl(lineread))=='START BACKFLOW')exit
if(am_master)then
if(xyzzyaaab17<0)call errstop('INIT_BACKFLOW','Cannot find "START BACK&
&FLOW" in correlation.data.')
if(xyzzyaaab17>0)call errstop('INIT_BACKFLOW','Problem reading correla&
&tion.data.')
endif
enddo
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,'(a)',err=100,end=101)title
if(am_master)then
call wout('Backflow function:')
call wout(' Title                  :  '//trim(adjustl(title)))
endif
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaad1
if(am_master.and.xyzzyaaad1<0)call errstop('INIT_BACKFLOW','Truncation&
& order must not be negative.')
xyzzyaaaf1=xyzzyaaad1-2
xyzzyaaae1=xyzzyaaad1-1
if(am_master)then
call wout(' Truncation order       :  '//trim(i2s(xyzzyaaad1)))
select case(xyzzyaaad1)
case(0)
tmps='Wave function, local energy'
case(1)
tmps='Derivatives of wave function, local energy'
case(2)
tmps='Local energy'
case default
tmps='None'
end select
call wout(' Discontinuities        :  '//trim(tmps))
endif
mainloop: do
read(xyzzyaadp1,'(a)',err=100,end=101)term_label
if(trim(adjustl(term_label))=="START ETA TERM")then
if(am_master)then
if(xyzzyaaaa1)call errstop('INIT_BACKFLOW','Only one eta term allowed.&
&')
call wout()
call wout(' Eta term:')
endif
xyzzyaaaa1=.true.
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaak17
if(am_master.and.xyzzyaaak17<1)call errstop('INIT_BACKFLOW','Expansion&
& order has to be greater than 0.')
xyzzyaaag1=xyzzyaaak17+1
if(am_master)then
call wout('  Expansion order       :  '//trim(i2s(xyzzyaaag1-1)))
endif
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaaj1
if(am_master)then
call wout('  Spin dependence       :  '//trim(i2s(xyzzyaaaj1)))
if(xyzzyaaaj1<-custom_spairs.or.xyzzyaaaj1>levels_spairs)call errstop(&
&'INIT_BACKFLOW','Spin-pair dependence has to be ' //trim(i2s(-custom_&
&spairs))//' -- '//trim(i2s(levels_spairs)) //' for eta.')
if(xyzzyaaaj1/=0.and.noncoll_spin)call errstop('INIT_BACKFLOW','Spin-d&
&ependence of eta should be zero for noncollinear-spin calculations.')
endif
xyzzyaaak1=no_spairs(xyzzyaaaj1)
allocate(xyzzyaabh1(xyzzyaaak1),xyzzyaabo1(xyzzyaaak1),xyzzyaacr1(xyzz&
&yaaak1),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','1')
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaak1)
xyzzyaaba17=xyzzyaaba17+dble(2*xyzzyaaak1)
call xyzzyaagy1(xyzzyaaaj1,xyzzyaacr1)
read(xyzzyaadp1,*,err=100,end=101)
read_cutoff_eta: do s=1,xyzzyaaak1
read(xyzzyaadp1,*,iostat=xyzzyaaab17)xyzzyaaay17,xyzzyaaak17
if(xyzzyaaab17/=0)then
if(am_master.and.s==1)call errstop('INIT_BACKFLOW','Error reading back&
&flow set in correlation.data')
backspace xyzzyaadp1
do xyzzyaaac17=s,xyzzyaaak1
if(xyzzyaabh17)then
xyzzyaabh1(xyzzyaaac17)=xyzzyaagz1(xyzzyaaac17)
if(am_master)then
tmpr=trim(r2s(xyzzyaabh1(xyzzyaaac17),'(f21.12)'))//' (default)'
if(xyzzyaabo1(1)==1)then
call wout('  Cutoff for s='//trim(i2s(xyzzyaaac17))//'  (opt) :  '//tr&
&im(tmpr))
else
call wout('  Cutoff for s='//trim(i2s(xyzzyaaac17))//' (fixd) :  '//tr&
&im(tmpr))
endif
endif
else
xyzzyaabh1(xyzzyaaac17)=xyzzyaabh1(1)
endif
xyzzyaabo1(xyzzyaaac17)=xyzzyaabo1(1)
enddo
if(am_master.and..not.xyzzyaabh17)then
call wout('  Unspecified cutoffs   :  Set to L_eta(1)')
endif
exit read_cutoff_eta
endif
if(am_master.and.s==1.and.(xyzzyaaak17<0.or.xyzzyaaak17>2))call errsto&
&p('INIT_BACKFLOW','Optimizable flag should be 0, 1 or 2 for L_eta(1).&
&')
if(am_master.and.s>1.and.(xyzzyaaak17<0.or.xyzzyaaak17>1))call errstop&
&('INIT_BACKFLOW','Optimizable flag should be 0 or 1 for L_eta(s), s>1&
&.')
if(abs(xyzzyaaay17)<1.d-8)then
xyzzyaaay17=xyzzyaagz1(s)
if(s==1)xyzzyaabh17=.true.
if(am_master)tmpr=trim(r2s(xyzzyaaay17,'(f21.12)'))//' (default)'
elseif(am_master)then
if(isperiodic)then
if(xyzzyaaay17>wigner_seitz_radius)call errstop('INIT_BACKFLOW','Cutof&
&f #'//trim(i2s(s))//' > radius of sphere inscribed in Wigner-Seitz ce&
&ll.')
endif
tmpr=r2s(xyzzyaaay17,'(f21.12)')
endif
if(xyzzyaaak17==2)then
if(am_master)then
call wout('  All cutoffs     (opt) :  '//trim(tmpr))
endif
xyzzyaabh1(:)=xyzzyaaay17
xyzzyaabo1(:)=xyzzyaaak17
exit read_cutoff_eta
else
if(am_master)then
if(xyzzyaaak17==1)then
call wout('  Cutoff for s='//trim(i2s(s))//'  (opt) :  '//trim(tmpr))
else
call wout('  Cutoff for s='//trim(i2s(s))//' (fixd) :  '//trim(tmpr))
endif
endif
xyzzyaabh1(s)=xyzzyaaay17
xyzzyaabo1(s)=xyzzyaaak17
endif
enddo read_cutoff_eta
if(am_master)then
xyzzyaaak17=xyzzyaaak1*xyzzyaaag1-count(xyzzyaacr1)
if(xyzzyaabo1(1)==2)then
call wout('  No. of free params    :  '//trim(i2s(xyzzyaaak17))//' + 1&
& cut-off length')
else
call wout('  No. of free params    :  '//trim(i2s(xyzzyaaak17))//' + '&
&//trim(i2s(xyzzyaaak1))//' cut-off lengths')
endif
endif
allocate(xyzzyaabi1(xyzzyaaag1,xyzzyaaak1),xyzzyaabp1(xyzzyaaag1,xyzzy&
&aaak1),stat=xyzzyaaaa17)
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaag1*xyzzyaaak1)
xyzzyaaba17=xyzzyaaba17+dble(xyzzyaaag1*xyzzyaaak1)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','2')
xyzzyaabi1(:,1:xyzzyaaak1)=0.d0
xyzzyaabp1(:,:)=1
read(xyzzyaadp1,*,err=100,end=101)
read_params_eta: do s=1,xyzzyaaak1
do xyzzyaaad17=1,xyzzyaaag1
if(xyzzyaaad17==2.and.xyzzyaacr1(s))cycle
read(xyzzyaadp1,*,iostat=xyzzyaaab17)xyzzyaaay17,xyzzyaaak17
if(xyzzyaaab17/=0)then
backspace xyzzyaadp1
if(am_master)call wout('  Unspecified params    :  Zeroed')
exit read_params_eta
endif
empty_backflow=.false.
if(am_master.and.xyzzyaaak17/=0.and.xyzzyaaak17/=1)call errstop('INIT_&
&BACKFLOW','Optimizable flag should be 0 or 1.')
xyzzyaabi1(xyzzyaaad17,s)=xyzzyaaay17
xyzzyaabp1(xyzzyaaad17,s)=xyzzyaaak17
if(am_master)then
tmpr=r2s(xyzzyaabi1(xyzzyaaad17,s),'(f21.12)')
if(xyzzyaabi1(xyzzyaaad17,s)<0.d0)then
if(xyzzyaabp1(xyzzyaaad17,s)==1)then
call wout('  c_'//trim(i2s(xyzzyaaad17-1))//','//trim(i2s(s))//'   (op&
&timizable) : '//trim(tmpr))
else
call wout('  c_'//trim(i2s(xyzzyaaad17-1))//','//trim(i2s(s))//'      &
&   (fixed) : '//trim(tmpr))
endif
else
if(xyzzyaabp1(xyzzyaaad17,s)==1)then
call wout('  c_'//trim(i2s(xyzzyaaad17-1))//','//trim(i2s(s))//'   (op&
&timizable) :  '//trim(tmpr))
else
call wout('  c_'//trim(i2s(xyzzyaaad17-1))//','//trim(i2s(s))//'      &
&   (fixed) :  '//trim(tmpr))
endif
endif
endif
enddo
enddo read_params_eta
read(xyzzyaadp1,'(a)',iostat=xyzzyaaab17)lineread
if(am_master.and.xyzzyaaab17/=0)call errstop('INIT_BACKFLOW','String "&
&END ETA TERM" not found.')
if(trim(adjustl(lineread))/='END ETA TERM')call errstop('INIT_BACKFLOW&
&','String "END ETA TERM" not found.')
elseif(trim(adjustl(term_label))=='START MU TERM')then
if(am_master)then
if(xyzzyaaab1)call errstop('INIT_BACKFLOW','Mu term found twice. Use m&
&ultiple sets instead of multiple terms.')
if(am_master.and.nitot==0)call errstop('INIT_BACKFLOW','No mu term nee&
&ded in a system without atoms.')
call wout()
call wout(' Mu term:')
endif
xyzzyaaab1=.true.
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,'(a)',err=100,end=101)lineread
read(lineread,*,iostat=xyzzyaaab17)xyzzyaaah1,xyzzyaabb1
if(xyzzyaaab17/=0)then
read(lineread,*,err=100,end=101)xyzzyaaah1
xyzzyaabb1=1
endif
if(am_master)then
call wout('  Number of sets        :  '//trim(i2s(xyzzyaaah1)))
if(xyzzyaabb1<1.or.xyzzyaabb1>3)call errstop('INIT_BACKFLOW','Label st&
&yle should be 1, 2 or 3.')
if(xyzzyaaah1<1.or.xyzzyaaah1>nitot.or.(xyzzyaabb1==2.and.xyzzyaaah1>n&
&basis).or.(xyzzyaabb1==3.and.xyzzyaaah1>nitype))call errstop('INIT_BA&
&CKFLOW','Problematic number of mu atom sets.')
endif
allocate(xyzzyaaan1(nitot),xyzzyaaap1(xyzzyaaah1),xyzzyaabj1(xyzzyaaah&
&1),xyzzyaabq1(xyzzyaaah1),xyzzyaaao1(xyzzyaaah1),xyzzyaaaq1(xyzzyaaah&
&1),xyzzyaaar1(nitot,xyzzyaaah1),xyzzyaaas1(xyzzyaaah1),xyzzyaact1(xyz&
&zyaaah1),xyzzyaabf1(xyzzyaaah1),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','3')
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaah1)
xyzzyaaba17=xyzzyaaba17+dble(nitot+xyzzyaaah1*(7+nitot))
xyzzyaaar1(1:nitot,1:xyzzyaaah1)=0
xyzzyaaan1(1:nitot)=0
if(xyzzyaabb1==1)then
allocate(xyzzyaabd1(nitot,xyzzyaaah1),stat=xyzzyaaaa17)
elseif(xyzzyaabb1==2)then
if(nbasis<=0)call errstop_master('INIT_BACKFLOW','Number of atoms in b&
&asis is zero.  If this is a Wigner crystal then use "atom" labels in &
&supercell.')
allocate(xyzzyaabd1(nbasis,xyzzyaaah1),stat=xyzzyaaaa17)
else
if(nitype<=0)call errstop_master('INIT_BACKFLOW','Number of atom speci&
&es is zero.  If this is a Wigner crystal then use "atom" labels in su&
&percell.')
allocate(xyzzyaabd1(nitype,xyzzyaaah1),stat=xyzzyaaaa17)
endif
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','3.1')
do xyzzyaaaj17=1,xyzzyaaah1
do
read(xyzzyaadp1,'(a)',err=100,end=101)lineread
if(trim(adjustl(lineread))=='START SET '//trim(i2s(xyzzyaaaj17)))exit
if(am_master.and.trim(adjustl(lineread))/='')call  errstop('INIT_BACKF&
&LOW','Was expecting to find "START SET '//trim(i2s(xyzzyaaaj17))//'".&
&')
enddo
if(am_master)call wout('  Set '//trim(i2s(xyzzyaaaj17)))
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaabf1(xyzzyaaaj17)
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaabd1(1:xyzzyaabf1(xyzzyaaaj17)&
&,xyzzyaaaj17)
if(xyzzyaabb1==1)then
xyzzyaaao1(xyzzyaaaj17)=xyzzyaabf1(xyzzyaaaj17)
xyzzyaaar1(1:xyzzyaabf1(xyzzyaaaj17),xyzzyaaaj17)=xyzzyaabd1(1:xyzzyaa&
&bf1(xyzzyaaaj17),xyzzyaaaj17)
elseif(xyzzyaabb1==2)then
if(xyzzyaabf1(xyzzyaaaj17)<1.or.xyzzyaabf1(xyzzyaaaj17)>nbasis+1-xyzzy&
&aaah1)call errstop_master('INIT_BACKFLOW','Problematic number of atom&
&s in set.')
if(any(xyzzyaabd1(1:xyzzyaabf1(xyzzyaaaj17),xyzzyaaaj17)<1) .or.any(xy&
&zzyaabd1(1:xyzzyaabf1(xyzzyaaaj17),xyzzyaaaj17)>nbasis))call errstop_&
&master('INIT_BACKFLOW','Problem with atom labels.')
call atoms_label_pcell(xyzzyaabf1(xyzzyaaaj17),xyzzyaabd1(1:xyzzyaabf1&
&(xyzzyaaaj17),xyzzyaaaj17),xyzzyaaao1(xyzzyaaaj17),xyzzyaaar1(:,xyzzy&
&aaaj17))
else
if(xyzzyaabf1(xyzzyaaaj17)<1.or.xyzzyaabf1(xyzzyaaaj17)>nitype+1-xyzzy&
&aaah1)call errstop_master('INIT_BACKFLOW','Problematic number of spec&
&ies in set.')
if(any(xyzzyaabd1(1:xyzzyaabf1(xyzzyaaaj17),xyzzyaaaj17)<1).or.any(xyz&
&zyaabd1(1:xyzzyaabf1(xyzzyaaaj17),xyzzyaaaj17)>nitype))call errstop_m&
&aster('INIT_BACKFLOW','Problem with species labels.')
call atoms_label_species(xyzzyaabf1(xyzzyaaaj17),xyzzyaabd1(1:xyzzyaab&
&f1(xyzzyaaaj17),xyzzyaaaj17),xyzzyaaao1(xyzzyaaaj17),xyzzyaaar1(:,xyz&
&zyaaaj17))
endif
if(am_master)then
call wout('   Atoms in set         :  '//trim(i2s(xyzzyaaao1(xyzzyaaaj&
&17))))
if(xyzzyaaao1(xyzzyaaaj17)<1.or.xyzzyaaao1(xyzzyaaaj17)>nitot+1-xyzzya&
&aah1)call errstop('INIT_BACKFLOW','Problem with number of atoms in se&
&t.')
xyzzyaaaf17=0
xyzzyaabg17=.true.
do while(xyzzyaaaf17<xyzzyaaao1(xyzzyaaaj17))
tmps=''
do xyzzyaaae17=xyzzyaaaf17+1,xyzzyaaao1(xyzzyaaaj17)
if(len_trim(adjustl(tmps))+len_trim(i2s(xyzzyaaar1(xyzzyaaae17,xyzzyaa&
&aj17)))+1>50)then
xyzzyaaaf17=xyzzyaaae17-1
exit
endif
tmps=trim(adjustl(tmps))//' '//trim(i2s(xyzzyaaar1(xyzzyaaae17,xyzzyaa&
&aj17)))
xyzzyaaaf17=xyzzyaaae17
enddo
if(xyzzyaabg17)then
call wout('   The atoms are        :  '//trim(adjustl(tmps)))
xyzzyaabg17=.false.
else
call wout('                           '//trim(adjustl(tmps)))
endif
enddo
if(any(xyzzyaaar1(1:xyzzyaaao1(xyzzyaaaj17),xyzzyaaaj17)<1) .or.any(xy&
&zzyaaar1(1:xyzzyaaao1(xyzzyaaaj17),xyzzyaaaj17)>nitot))call errstop('&
&INIT_BACKFLOW','Problem with atom labels.')
do xyzzyaaae17=1,xyzzyaaao1(xyzzyaaaj17)-1
do xyzzyaaaf17=xyzzyaaae17+1,xyzzyaaao1(xyzzyaaaj17)
if(xyzzyaaar1(xyzzyaaae17,xyzzyaaaj17)==xyzzyaaar1(xyzzyaaaf17,xyzzyaa&
&aj17))call errstop('INIT_BACKFLOW','Ion '//trim(i2s(xyzzyaaar1(xyzzya&
&aae17,xyzzyaaaj17)))// ' appears twice.')
enddo
enddo
endif
do xyzzyaaae17=1,xyzzyaaao1(xyzzyaaaj17)
if(xyzzyaaan1(xyzzyaaar1(xyzzyaaae17,xyzzyaaaj17))==0)then
xyzzyaaan1(xyzzyaaar1(xyzzyaaae17,xyzzyaaaj17))=xyzzyaaaj17
else
if(am_master)then
call wout('Ion '//trim(i2s(xyzzyaaar1(xyzzyaaae17,xyzzyaaaj17)))//' ap&
&pears in sets '// trim(i2s(xyzzyaaan1(xyzzyaaar1(xyzzyaaae17,xyzzyaaa&
&j17))))//' and '//trim(i2s(xyzzyaaaj17))//'.')
call errstop('INIT_BACKFLOW','Stopping.')
endif
endif
enddo
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaak17
if(am_master.and.(xyzzyaaak17<0.or.xyzzyaaak17>1))call errstop('INIT_B&
&ACKFLOW','Flag for e-N cusp type must be 0 or 1.')
if(xyzzyaaak17==0)then
xyzzyaact1(xyzzyaaaj17)=.false.
if(am_master)call wout('   Type of cusp conds.  :  PP / cuspless AE')
else
xyzzyaact1(xyzzyaaaj17)=.true.
if(am_master)call wout('   Type of cusp conds.  :  AE with cusp')
endif
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaap1(xyzzyaaaj17)
xyzzyaaap1(xyzzyaaaj17)=xyzzyaaap1(xyzzyaaaj17)+1
if(am_master)then
call wout('   Expansion order      :  '//trim(i2s(xyzzyaaap1(xyzzyaaaj&
&17)-1)))
if(xyzzyaaap1(xyzzyaaaj17)<2)call errstop('INIT_BACKFLOW','N_mu<1')
endif
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaaq1(xyzzyaaaj17)
if(am_master)then
call wout('   Spin dependence      :  '//trim(i2s(xyzzyaaaq1(xyzzyaaaj&
&17))))
if(xyzzyaaaq1(xyzzyaaaj17)<-custom_ssingles.or.xyzzyaaaq1(xyzzyaaaj17)&
& >levels_ssingles)call errstop ('INIT_BACKFLOW','Spin dependence shou&
&ld be '//trim(i2s(-custom_ssingles))//' -- '//trim(i2s(levels_ssingle&
&s))//' for mu.')
endif
xyzzyaaas1(xyzzyaaaj17)=no_ssingles(xyzzyaaaq1(xyzzyaaaj17))
if(xyzzyaaaj17==1)then
allocate(xyzzyaabk1(xyzzyaaap1(xyzzyaaaj17),xyzzyaaas1(xyzzyaaaj17),xy&
&zzyaaah1),xyzzyaabr1(xyzzyaaap1(xyzzyaaaj17),xyzzyaaas1(xyzzyaaaj17),&
&xyzzyaaah1),stat=xyzzyaaaa17)
xyzzyaaak17=xyzzyaaap1(xyzzyaaaj17)*xyzzyaaas1(xyzzyaaaj17)*xyzzyaaah1
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaak17)
xyzzyaaba17=xyzzyaaba17+dble(xyzzyaaak17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','4')
xyzzyaabk1(:,:,:)=0.d0
xyzzyaabr1(:,:,:)=1
else
xyzzyaaan17=shape(xyzzyaabk1)
if(xyzzyaaap1(xyzzyaaaj17)>xyzzyaaan17(1).or.xyzzyaaas1(xyzzyaaaj17)>x&
&yzzyaaan17(2))then
allocate(xyzzyaabc17(xyzzyaaan17(1),xyzzyaaan17(2),xyzzyaaaj17-1),xyzz&
&yaaau17(xyzzyaaan17(1),xyzzyaaan17(2),xyzzyaaaj17-1),stat=xyzzyaaaa17&
&)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','5.')
xyzzyaabc17(:,:,1:xyzzyaaaj17-1)=xyzzyaabk1(:,:,1:xyzzyaaaj17-1)
xyzzyaaau17(:,:,1:xyzzyaaaj17-1)=xyzzyaabr1(:,:,1:xyzzyaaaj17-1)
deallocate(xyzzyaabk1,xyzzyaabr1)
allocate(xyzzyaabk1(maxval(xyzzyaaap1(1:xyzzyaaaj17)),maxval(xyzzyaaas&
&1(1:xyzzyaaaj17)),xyzzyaaah1),xyzzyaabr1(maxval(xyzzyaaap1(1:xyzzyaaa&
&j17)),maxval(xyzzyaaas1(1:xyzzyaaaj17)),xyzzyaaah1),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','6')
xyzzyaabk1(:,:,:)=0.d0
xyzzyaabr1(:,:,:)=1
xyzzyaabk1(1:xyzzyaaan17(1),1:xyzzyaaan17(2),1:xyzzyaaaj17-1)=xyzzyaab&
&c17(1:xyzzyaaan17(1),1:xyzzyaaan17(2),1:xyzzyaaaj17-1)
xyzzyaabr1(1:xyzzyaaan17(1),1:xyzzyaaan17(2),1:xyzzyaaaj17-1)=xyzzyaaa&
&u17(1:xyzzyaaan17(1),1:xyzzyaaan17(2),1:xyzzyaaaj17-1)
deallocate(xyzzyaabc17,xyzzyaaau17)
xyzzyaaak17=(maxval(xyzzyaaap1(1:xyzzyaaaj17))*maxval(xyzzyaaas1(1:xyz&
&zyaaaj17))-xyzzyaaan17(1)*xyzzyaaan17(2))*xyzzyaaah1
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaak17)
xyzzyaaba17=xyzzyaaba17+dble(xyzzyaaak17)
endif
endif
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaabj1(xyzzyaaaj17),xyzzyaabq1(x&
&yzzyaaaj17)
if(am_master.and.xyzzyaabq1(xyzzyaaaj17)/=0.and.xyzzyaabq1(xyzzyaaaj17&
&)/=1)call errstop('INIT_BACKFLOW','Optimizable flag should be 0 or 1 &
&.')
if(abs(xyzzyaabj1(xyzzyaaaj17))<1.d-8)then
xyzzyaabj1(xyzzyaaaj17)=xyzzyaaha1()
if(am_master)tmpr=trim(r2s(xyzzyaabj1(xyzzyaaaj17),'(f21.12)'))//' (de&
&fault)'
else
if(isperiodic)then
if(am_master.and.xyzzyaabj1(xyzzyaaaj17)>wigner_seitz_radius)call errs&
&top('INIT_BACKFLOW','Cutoff > radius of sphere inscribed in Wigner-Se&
&itz cell.')
endif
if(am_master)tmpr=r2s(xyzzyaabj1(xyzzyaaaj17),'(f21.12)')
endif
if(am_master)then
if(xyzzyaabq1(xyzzyaaaj17)==1)then
call wout('   Cutoff (optimizable) :  '//trim(tmpr))
else
call wout('   Cutoff       (fixed) :  '//trim(tmpr))
endif
endif
if(am_master)then
xyzzyaaak17=(xyzzyaaap1(xyzzyaaaj17)-1)*xyzzyaaas1(xyzzyaaaj17)
if(xyzzyaact1(xyzzyaaaj17))xyzzyaaak17=xyzzyaaak17-xyzzyaaas1(xyzzyaaa&
&j17)
call wout('   No. of free params   :  '//trim(i2s(xyzzyaaak17))//' + c&
&ut-off length')
if(xyzzyaaak17<1)call errstop('INIT_BACKFLOW','No free parameters in s&
&et.')
endif
read(xyzzyaadp1,*,err=100,end=101)
read_params_mu: do s=1,xyzzyaaas1(xyzzyaaaj17)
do xyzzyaaad17=1,xyzzyaaap1(xyzzyaaaj17)
if(xyzzyaaad17==2.or.(xyzzyaaad17==1.and.xyzzyaact1(xyzzyaaaj17)))cycl&
&e
read(xyzzyaadp1,*,iostat=xyzzyaaab17)xyzzyaaay17,xyzzyaaak17
if(xyzzyaaab17/=0)then
backspace xyzzyaadp1
if(am_master)call wout('   Unspecified params   :  Zeroed')
exit read_params_mu
endif
empty_backflow=.false.
if(am_master.and.xyzzyaaak17/=0.and.xyzzyaaak17/=1)call errstop('INIT_&
&BACKFLOW','Optimizable flag should be 0 or 1 .')
xyzzyaabk1(xyzzyaaad17,s,xyzzyaaaj17)=xyzzyaaay17
xyzzyaabr1(xyzzyaaad17,s,xyzzyaaaj17)=xyzzyaaak17
if(am_master)then
tmpr=r2s(xyzzyaabk1(xyzzyaaad17,s,xyzzyaaaj17),'(f21.12)')
if(xyzzyaabk1(xyzzyaaad17,s,xyzzyaaaj17)<0.d0)then
if(xyzzyaabr1(xyzzyaaad17,s,xyzzyaaaj17)==1)then
call wout('   mu_'//trim(i2s(xyzzyaaad17-1))//','//trim(i2s(s))//' (op&
&timizable) : '//trim(tmpr))
else
call wout('   mu_'//trim(i2s(xyzzyaaad17-1))//','//trim(i2s(s))//'    &
&   (fixed) : '//trim(tmpr))
endif
else
if(xyzzyaabr1(xyzzyaaad17,s,xyzzyaaaj17)==1)then
call wout('   mu_'//trim(i2s(xyzzyaaad17-1))//','//trim(i2s(s))//' (op&
&timizable) :  '//trim(tmpr))
else
call wout('   mu_'//trim(i2s(xyzzyaaad17-1))//','//trim(i2s(s))//'    &
&   (fixed) :  '//trim(tmpr))
endif
endif
endif
enddo
enddo read_params_mu
read(xyzzyaadp1,'(a)',iostat=xyzzyaaab17)lineread
if(am_master.and.xyzzyaaab17/=0)call errstop('INIT_BACKFLOW','Can''t f&
&ind "END SET '//trim(i2s(xyzzyaaaj17))//'".')
if(trim(adjustl(lineread))/='END SET '//trim(i2s(xyzzyaaaj17)))call er&
&rstop('INIT_BACKFLOW','Can''t find "END SET '//trim(i2s(xyzzyaaaj17))&
&//'".')
enddo
if(am_master)then
xyzzyaaak17=count(xyzzyaaan1==0)
if(xyzzyaaak17==0)then
call wout('  Completeness of mu    :  All atoms included')
else
if(allow_nochi_atoms)then
if(xyzzyaaak17==1)then
call wout('  Completeness of mu    :  1 atom does not belong to any se&
&t')
else
call wout('  Completeness of mu    :  '//trim(i2s(xyzzyaaak17))//' ato&
&ms do not belong to any set')
endif
else
call errstop('INIT_BACKFLOW','Some atoms do not belong to any mu set.'&
&)
endif
endif
endif
do
read(xyzzyaadp1,'(a)',err=100,end=101)lineread
if(trim(adjustl(lineread))=='END MU TERM')exit
if(am_master.and.trim(adjustl(lineread))/='')call errstop('INIT_BACKFL&
&OW','String "END MU TERM" not found.')
enddo
elseif(trim(adjustl(term_label))=='START PHI TERM')then
if(am_master)then
if(xyzzyaaac1)call errstop('INIT_BACKFLOW','Phi term found twice. Use &
&multiple sets instead of multiple terms.')
if(nitot==0)call errstop('INIT_BACKFLOW','No phi term needed in a syst&
&em without atoms.')
call wout()
call wout(' Phi term:')
endif
xyzzyaaac1=.true.
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,'(a)',err=100,end=101)lineread
read(lineread,*,iostat=xyzzyaaab17)xyzzyaaai1,xyzzyaabc1
if(xyzzyaaab17/=0)then
read(lineread,*,err=100,end=101)xyzzyaaai1
xyzzyaabc1=1
endif
if(am_master)then
call wout('  Number of sets        :  '//trim(i2s(xyzzyaaai1)))
if(xyzzyaabc1<1.or.xyzzyaabc1>3)call errstop('INIT_BACKFLOW','Label st&
&yle should be 1, 2 or 3.')
if(xyzzyaaai1<1.or.xyzzyaaai1>nitot.or.(xyzzyaabc1==2 .and.xyzzyaaai1>&
&nbasis).or.(xyzzyaabc1==3.and.xyzzyaaai1>nitype))call errstop('INIT_B&
&ACKFLOW','Problematic number of phi atom sets.')
endif
allocate(xyzzyaaat1(nitot),xyzzyaaaw1(xyzzyaaai1),xyzzyaaav1(xyzzyaaai&
&1),xyzzyaabl1(xyzzyaaai1),xyzzyaabs1(xyzzyaaai1),xyzzyaaau1(xyzzyaaai&
&1),xyzzyaaay1(xyzzyaaai1),xyzzyaaax1(nitot,xyzzyaaai1),xyzzyaaaz1(xyz&
&zyaaai1),xyzzyaacu1(xyzzyaaai1),xyzzyaacv1(xyzzyaaai1),xyzzyaabg1(xyz&
&zyaaai1),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','7')
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaai1)
xyzzyaaba17=xyzzyaaba17+nitot+dble(xyzzyaaai1*(8+nitot))
xyzzyaaax1(1:nitot,1:xyzzyaaai1)=0
xyzzyaaat1(1:nitot)=0
if(xyzzyaabc1==1)then
allocate(xyzzyaabe1(nitot,xyzzyaaai1),stat=xyzzyaaaa17)
elseif(xyzzyaabc1==2)then
if(nbasis<=0)call errstop_master('INIT_BACKFLOW','Number of atoms in b&
&asis is zero.  If this is a Wigner crystal then use "atom" labels in &
&supercell.')
allocate(xyzzyaabe1(nbasis,xyzzyaaai1),stat=xyzzyaaaa17)
else
if(nitype<=0)call errstop_master('INIT_BACKFLOW','Number of atom speci&
&es is zero.  If this is a Wigner crystal then use "atom" labels in su&
&percell.')
allocate(xyzzyaabe1(nitype,xyzzyaaai1),stat=xyzzyaaaa17)
endif
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','7.1')
do xyzzyaaaj17=1,xyzzyaaai1
do
read(xyzzyaadp1,'(a)',err=100,end=101)lineread
if(trim(adjustl(lineread))=='START SET '//trim(i2s(xyzzyaaaj17)))exit
if(am_master.and.trim(adjustl(lineread))/='')call errstop('INIT_BACKFL&
&OW','Was expecting to find "START SET ' //trim(i2s(xyzzyaaaj17))//'".&
&')
enddo
if(am_master)call wout('  Set '//trim(i2s(xyzzyaaaj17)))
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaabg1(xyzzyaaaj17)
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaabe1(1:xyzzyaabg1(xyzzyaaaj17)&
&,xyzzyaaaj17)
if(xyzzyaabc1==1)then
xyzzyaaau1(xyzzyaaaj17)=xyzzyaabg1(xyzzyaaaj17)
xyzzyaaax1(1:xyzzyaabg1(xyzzyaaaj17),xyzzyaaaj17)=xyzzyaabe1(1:xyzzyaa&
&bg1(xyzzyaaaj17),xyzzyaaaj17)
elseif(xyzzyaabc1==2)then
if(xyzzyaabg1(xyzzyaaaj17)<1.or.xyzzyaabg1(xyzzyaaaj17)>nbasis+1-xyzzy&
&aaai1)call errstop_master('INIT_BACKFLOW','Problematic number of atom&
&s in set.')
if(any(xyzzyaabe1(1:xyzzyaabg1(xyzzyaaaj17),xyzzyaaaj17)<1) .or.any(xy&
&zzyaabe1(1:xyzzyaabg1(xyzzyaaaj17),xyzzyaaaj17)>nbasis))call errstop_&
&master('INIT_BACKFLOW','Problem with atom labels.')
call atoms_label_pcell(xyzzyaabg1(xyzzyaaaj17),xyzzyaabe1(1:xyzzyaabg1&
&(xyzzyaaaj17),xyzzyaaaj17),xyzzyaaau1(xyzzyaaaj17),xyzzyaaax1(:,xyzzy&
&aaaj17))
else
if(xyzzyaabg1(xyzzyaaaj17)<1.or.xyzzyaabg1(xyzzyaaaj17)>nitype+1-xyzzy&
&aaai1)call errstop_master('INIT_BACKFLOW','Problematic number of spec&
&ies in set.')
if(any(xyzzyaabe1(1:xyzzyaabg1(xyzzyaaaj17),xyzzyaaaj17)<1).or.any(xyz&
&zyaabe1(1:xyzzyaabg1(xyzzyaaaj17),xyzzyaaaj17)>nitype))call errstop_m&
&aster('INIT_BACKFLOW','Problematic label of species in set.')
call atoms_label_species(xyzzyaabg1(xyzzyaaaj17),xyzzyaabe1(1:xyzzyaab&
&g1(xyzzyaaaj17),xyzzyaaaj17),xyzzyaaau1(xyzzyaaaj17),xyzzyaaax1(:,xyz&
&zyaaaj17))
endif
if(am_master)then
call wout('   Atoms in set         :  '//trim(i2s(xyzzyaaau1(xyzzyaaaj&
&17))))
if(xyzzyaaau1(xyzzyaaaj17)<1.or.xyzzyaaau1(xyzzyaaaj17)>nitot+1-xyzzya&
&aai1)call errstop('INIT_BACKFLOW','Problem with number of atoms in se&
&t.')
xyzzyaaaf17=0
xyzzyaabg17=.true.
do while(xyzzyaaaf17<xyzzyaaau1(xyzzyaaaj17))
tmps=''
do xyzzyaaae17=xyzzyaaaf17+1,xyzzyaaau1(xyzzyaaaj17)
if(len_trim(adjustl(tmps))+len_trim(i2s(xyzzyaaax1(xyzzyaaae17,xyzzyaa&
&aj17)))+1>50)then
xyzzyaaaf17=xyzzyaaae17-1
exit
endif
tmps=trim(adjustl(tmps))//' '//trim(i2s(xyzzyaaax1(xyzzyaaae17,xyzzyaa&
&aj17)))
xyzzyaaaf17=xyzzyaaae17
enddo
if(xyzzyaabg17)then
call wout('   The atoms are        :  '//trim(adjustl(tmps)))
xyzzyaabg17=.false.
else
call wout('                           '//trim(adjustl(tmps)))
endif
enddo
if(any(xyzzyaaax1(1:xyzzyaaau1(xyzzyaaaj17),xyzzyaaaj17)<1) .or.any(xy&
&zzyaaax1(1:xyzzyaaau1(xyzzyaaaj17),xyzzyaaaj17)>nitot))call errstop('&
&INIT_BACKFLOW','Problem with atom labels.')
do xyzzyaaae17=1,xyzzyaaau1(xyzzyaaaj17)-1
do xyzzyaaaf17=xyzzyaaae17+1,xyzzyaaau1(xyzzyaaaj17)
if(xyzzyaaax1(xyzzyaaae17,xyzzyaaaj17)==xyzzyaaax1(xyzzyaaaf17,xyzzyaa&
&aj17))call errstop('INIT_BACKFLOW','Ion '//trim(i2s(xyzzyaaax1(xyzzya&
&aae17,xyzzyaaaj17)))// ' appears twice.')
enddo
enddo
endif
do xyzzyaaae17=1,xyzzyaaau1(xyzzyaaaj17)
if(xyzzyaaat1(xyzzyaaax1(xyzzyaaae17,xyzzyaaaj17))==0)then
xyzzyaaat1(xyzzyaaax1(xyzzyaaae17,xyzzyaaaj17))=xyzzyaaaj17
elseif(am_master)then
call wout('Ion '//trim(i2s(xyzzyaaax1(xyzzyaaae17,xyzzyaaaj17)))//' ap&
&pears in sets '//trim(i2s(xyzzyaaat1(xyzzyaaax1(xyzzyaaae17,xyzzyaaaj&
&17))))//' and '//trim(i2s(xyzzyaaaj17))//'.')
call errstop('INIT_BACKFLOW','Stopping.')
endif
enddo
read(xyzzyaadp1,'(a)',err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaak17
if(am_master.and.(xyzzyaaak17<0.or.xyzzyaaak17>1))call errstop('INIT_B&
&ACKFLOW','Flag for e-N cusp type must be 0 or 1.')
if(xyzzyaaak17==0)then
xyzzyaacu1(xyzzyaaaj17)=.false.
if(am_master)call wout('   Type of cusp conds.  :  PP / cuspless AE')
else
xyzzyaacu1(xyzzyaaaj17)=.true.
if(am_master)call wout('   Type of cusp conds.  :  AE with cusp')
endif
read(xyzzyaadp1,'(a)',err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaak17
if(am_master.and.(xyzzyaaak17<0.or.xyzzyaaak17>1))call errstop('INIT_B&
&ACKFLOW','Flag for irrotational backflow must be 0 or 1.')
if(xyzzyaaak17==0)then
xyzzyaacv1(xyzzyaaaj17)=.false.
if(am_master)call wout('   Irrotational constr. :  Not applied')
else
xyzzyaacv1(xyzzyaaaj17)=.true.
if(am_master)call wout('   Irrotational constr. :  Applied')
endif
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaaw1(xyzzyaaaj17)
xyzzyaaaw1(xyzzyaaaj17)=xyzzyaaaw1(xyzzyaaaj17)+1
if(am_master)then
call wout('   Expansion order e-N  :  '//trim(i2s(xyzzyaaaw1(xyzzyaaaj&
&17)-1)))
if(xyzzyaaaw1(xyzzyaaaj17)<2)call errstop('INIT_BACKFLOW','N_phi_eN<1'&
&)
endif
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaav1(xyzzyaaaj17)
xyzzyaaav1(xyzzyaaaj17)=xyzzyaaav1(xyzzyaaaj17)+1
if(am_master)then
call wout('   Expansion order e-e  :  '//trim(i2s(xyzzyaaav1(xyzzyaaaj&
&17)-1)))
if(xyzzyaaav1(xyzzyaaaj17)<2)call errstop('INIT_BACKFLOW','N_phi_ee<1'&
&)
endif
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaaay1(xyzzyaaaj17)
if(am_master)then
call wout('   Spin dependence      :  '//trim(i2s(xyzzyaaay1(xyzzyaaaj&
&17))))
if(xyzzyaaay1(xyzzyaaaj17)<-custom_spairs.or.xyzzyaaay1(xyzzyaaaj17)>l&
&evels_spairs) call errstop ('INIT_BACKFLOW','Spin-pair dependence sho&
&uld be ' //trim(i2s(-custom_spairs))//' -- '//trim(i2s(levels_spairs)&
&) //' for phi.')
endif
xyzzyaaaz1(xyzzyaaaj17)=no_spairs(xyzzyaaay1(xyzzyaaaj17))
xyzzyaaap17=xyzzyaaaw1(xyzzyaaaj17)
xyzzyaaaq17=xyzzyaaav1(xyzzyaaaj17)
xyzzyaaar17=xyzzyaaaz1(xyzzyaaaj17)
xyzzyaaas17=xyzzyaaai1
if(xyzzyaaaj17==1)then
allocate(xyzzyaabm1(xyzzyaaap17,xyzzyaaap17,xyzzyaaaq17,xyzzyaaar17,xy&
&zzyaaas17),xyzzyaabt1(xyzzyaaap17,xyzzyaaap17,xyzzyaaaq17,xyzzyaaar17&
&,xyzzyaaas17),xyzzyaabn1(xyzzyaaap17,xyzzyaaap17,xyzzyaaaq17,xyzzyaaa&
&r17,xyzzyaaas17),xyzzyaabu1(xyzzyaaap17,xyzzyaaap17,xyzzyaaaq17,xyzzy&
&aaar17,xyzzyaaas17),xyzzyaacs1(xyzzyaaar17,xyzzyaaas17),stat=xyzzyaaa&
&a17)
xyzzyaaak17=2*xyzzyaaap17*xyzzyaaap17*xyzzyaaaq17*xyzzyaaar17*xyzzyaaa&
&s17
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaak17)
xyzzyaaba17=xyzzyaaba17+dble(xyzzyaaak17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','8')
xyzzyaabm1(:,:,:,:,:)=0.d0
xyzzyaabn1(:,:,:,:,:)=0.d0
xyzzyaabt1(:,:,:,:,:)=1
xyzzyaabu1(:,:,:,:,:)=1
else
xyzzyaaam17=shape(xyzzyaabm1)
if(xyzzyaaap17>xyzzyaaam17(1).or.xyzzyaaaq17>xyzzyaaam17(3).or.xyzzyaa&
&ar17>xyzzyaaam17(4))then
allocate(xyzzyaabd17(xyzzyaaam17(1),xyzzyaaam17(2),xyzzyaaam17(3),xyzz&
&yaaam17(4),xyzzyaaaj17-1),xyzzyaaav17(xyzzyaaam17(1),xyzzyaaam17(2),x&
&yzzyaaam17(3),xyzzyaaam17(4),xyzzyaaaj17-1),xyzzyaabe17(xyzzyaaam17(1&
&),xyzzyaaam17(2),xyzzyaaam17(3),xyzzyaaam17(4),xyzzyaaaj17-1),xyzzyaa&
&aw17(xyzzyaaam17(1),xyzzyaaam17(2),xyzzyaaam17(3),xyzzyaaam17(4),xyzz&
&yaaaj17-1),xyzzyaabl17(xyzzyaaam17(4),xyzzyaaaj17-1),stat=xyzzyaaaa17&
&)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','9')
xyzzyaabd17(:,:,:,:,1:xyzzyaaaj17-1)=xyzzyaabm1(:,:,:,:,1:xyzzyaaaj17-&
&1)
xyzzyaaav17(:,:,:,:,1:xyzzyaaaj17-1)=xyzzyaabt1(:,:,:,:,1:xyzzyaaaj17-&
&1)
xyzzyaabe17(:,:,:,:,1:xyzzyaaaj17-1)=xyzzyaabn1(:,:,:,:,1:xyzzyaaaj17-&
&1)
xyzzyaaaw17(:,:,:,:,1:xyzzyaaaj17-1)=xyzzyaabu1(:,:,:,:,1:xyzzyaaaj17-&
&1)
xyzzyaabl17(:,1:xyzzyaaaj17-1)=xyzzyaacs1(:,1:xyzzyaaaj17-1)
deallocate(xyzzyaabm1,xyzzyaabt1,xyzzyaacs1,xyzzyaabn1,xyzzyaabu1)
xyzzyaaap17=maxval(xyzzyaaaw1(1:xyzzyaaaj17))
xyzzyaaaq17=maxval(xyzzyaaav1(1:xyzzyaaaj17))
xyzzyaaar17=maxval(xyzzyaaaz1(1:xyzzyaaaj17))
allocate(xyzzyaabm1(xyzzyaaap17,xyzzyaaap17,xyzzyaaaq17,xyzzyaaar17,xy&
&zzyaaas17),xyzzyaabt1(xyzzyaaap17,xyzzyaaap17,xyzzyaaaq17,xyzzyaaar17&
&,xyzzyaaas17),xyzzyaabn1(xyzzyaaap17,xyzzyaaap17,xyzzyaaaq17,xyzzyaaa&
&r17,xyzzyaaas17),xyzzyaabu1(xyzzyaaap17,xyzzyaaap17,xyzzyaaaq17,xyzzy&
&aaar17,xyzzyaaas17),xyzzyaacs1(xyzzyaaar17,xyzzyaaas17),stat=xyzzyaaa&
&a17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','10')
xyzzyaabm1(:,:,:,:,:)=0.d0
xyzzyaabn1(:,:,:,:,:)=0.d0
xyzzyaabt1(:,:,:,:,:)=1
xyzzyaabu1(:,:,:,:,:)=1
xyzzyaabm1(1:xyzzyaaam17(1),1:xyzzyaaam17(2),1:xyzzyaaam17(3),1:xyzzya&
&aam17(4),1:xyzzyaaaj17-1)=xyzzyaabd17(1:xyzzyaaam17(1),1:xyzzyaaam17(&
&2),1:xyzzyaaam17(3),1:xyzzyaaam17(4),1:xyzzyaaaj17-1)
xyzzyaabt1(1:xyzzyaaam17(1),1:xyzzyaaam17(2),1:xyzzyaaam17(3),1:xyzzya&
&aam17(4),1:xyzzyaaaj17-1)=xyzzyaaav17(1:xyzzyaaam17(1),1:xyzzyaaam17(&
&2),1:xyzzyaaam17(3),1:xyzzyaaam17(4),1:xyzzyaaaj17-1)
xyzzyaabn1(1:xyzzyaaam17(1),1:xyzzyaaam17(2),1:xyzzyaaam17(3),1:xyzzya&
&aam17(4),1:xyzzyaaaj17-1)=xyzzyaabe17(1:xyzzyaaam17(1),1:xyzzyaaam17(&
&2),1:xyzzyaaam17(3),1:xyzzyaaam17(4),1:xyzzyaaaj17-1)
xyzzyaabu1(1:xyzzyaaam17(1),1:xyzzyaaam17(2),1:xyzzyaaam17(3),1:xyzzya&
&aam17(4),1:xyzzyaaaj17-1)=xyzzyaaaw17(1:xyzzyaaam17(1),1:xyzzyaaam17(&
&2),1:xyzzyaaam17(3),1:xyzzyaaam17(4),1:xyzzyaaaj17-1)
xyzzyaacs1(1:xyzzyaaam17(4),1:xyzzyaaaj17-1)=xyzzyaabl17(1:xyzzyaaam17&
&(4),1:xyzzyaaaj17-1)
deallocate(xyzzyaabd17,xyzzyaaav17,xyzzyaabe17,xyzzyaaaw17,xyzzyaabl17&
&)
xyzzyaaak17=(xyzzyaaap17*xyzzyaaap17*xyzzyaaaq17*xyzzyaaar17-xyzzyaaam&
&17(1)*xyzzyaaam17(2)*xyzzyaaam17(3)*xyzzyaaam17(4))*xyzzyaaai1*2
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaak17)
xyzzyaaba17=xyzzyaaba17+dble(xyzzyaaak17)
endif
endif
call xyzzyaagy1(xyzzyaaay1(xyzzyaaaj17),xyzzyaacs1(:,xyzzyaaaj17))
read(xyzzyaadp1,*,err=100,end=101)
read(xyzzyaadp1,*,err=100,end=101)xyzzyaabl1(xyzzyaaaj17),xyzzyaabs1(x&
&yzzyaaaj17)
if(am_master.and.xyzzyaabs1(xyzzyaaaj17)/=0.and.xyzzyaabs1(xyzzyaaaj17&
&)/=1)call errstop('INIT_BACKFLOW','Optimizable flag should be 0 or 1 &
&.')
if(abs(xyzzyaabl1(xyzzyaaaj17))<1.d-8)then
xyzzyaabl1(xyzzyaaaj17)=xyzzyaahb1()
if(am_master)tmpr=trim(r2s(xyzzyaabl1(xyzzyaaaj17),'(f21.12)'))//' (de&
&fault)'
else
if(isperiodic)then
if(am_master.and.xyzzyaabl1(xyzzyaaaj17)>.5d0*wigner_seitz_radius)call&
& errstop('INIT_BACKFLOW','Cutoff > 0.5 * radius of sphere inscribed i&
&n Wigner-Seitz cell.')
endif
if(am_master)tmpr=r2s(xyzzyaabl1(xyzzyaaaj17),'(f21.12)')
endif
if(am_master)then
if(xyzzyaabs1(xyzzyaaaj17)==1)then
call wout('   Cutoff (optimizable) :  '//trim(tmpr))
else
call wout('   Cutoff       (fixed) :  '//trim(tmpr))
endif
endif
xyzzyaaat17=xyzzyaaav1(xyzzyaaaj17)*xyzzyaaaw1(xyzzyaaaj17)*xyzzyaaaw1&
&(xyzzyaaaj17)
if(allocated(xyzzyaabk17))deallocate(xyzzyaabk17)
allocate(xyzzyaabk17(xyzzyaaaz1(xyzzyaaaj17),2*xyzzyaaat17),stat=xyzzy&
&aaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','11')
xyzzyaaao17=0
do s=1,xyzzyaaaz1(xyzzyaaaj17)
call xyzzyaafl1(xyzzyaaaj17,2*xyzzyaaat17,xyzzyaacs1(s,xyzzyaaaj17),xy&
&zzyaabk17(s,:),xyzzyaaak17)
xyzzyaaao17=xyzzyaaao17+xyzzyaaak17
enddo
if(am_master)then
call wout('   No. of free params   :  '//trim(i2s(xyzzyaaao17))//' + c&
&ut-off length')
if(xyzzyaaao17<1)call errstop('INIT_BACKFLOW','No non-zero parameters &
&in set.')
endif
read(xyzzyaadp1,*,err=100,end=101)
read_params_phi: do s=1,xyzzyaaaz1(xyzzyaaaj17)
xyzzyaaad17=0
do xyzzyaaai17=1,xyzzyaaav1(xyzzyaaaj17)
do xyzzyaaah17=1,xyzzyaaaw1(xyzzyaaaj17)
do xyzzyaaag17=1,xyzzyaaaw1(xyzzyaaaj17)
xyzzyaaad17=xyzzyaaad17+1
if(xyzzyaabk17(s,xyzzyaaad17))cycle
read(xyzzyaadp1,*,iostat=xyzzyaaab17)xyzzyaaay17,xyzzyaaak17
if(xyzzyaaab17/=0)then
backspace xyzzyaadp1
if(am_master)call wout('   Unspecified params   :  Zeroed')
exit read_params_phi
endif
empty_backflow=.false.
if(am_master.and.xyzzyaaak17/=0.and.xyzzyaaak17/=1)call errstop('INIT_&
&BACKFLOW','Optimizable flag should be 0 or 1 .')
xyzzyaabm1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)=xyzzyaaa&
&y17
xyzzyaabt1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)=xyzzyaaa&
&k17
if(am_master)then
tmpr=r2s(xyzzyaabm1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)&
&,'(f21.12)')
if(xyzzyaabm1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)<0.d0)&
&then
if(xyzzyaabt1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)==1)th&
&en
call wout('   phi_'//trim(i2s(xyzzyaaag17-1))//','//trim(i2s(xyzzyaaah&
&17-1))//','//trim(i2s(xyzzyaaai17-1))//','//trim(i2s(s))//'    (opt) &
&: '//trim(tmpr))
else
call wout('   phi_'//trim(i2s(xyzzyaaag17-1))//','//trim(i2s(xyzzyaaah&
&17-1))//','//trim(i2s(xyzzyaaai17-1))//','//trim(i2s(s))//'  (fixed) &
&: '//trim(tmpr))
endif
else
if(xyzzyaabt1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)==1)th&
&en
call wout('   phi_'//trim(i2s(xyzzyaaag17-1))//','//trim(i2s(xyzzyaaah&
&17-1))//','//trim(i2s(xyzzyaaai17-1))//','//trim(i2s(s))//'    (opt) &
&:  '//trim(tmpr))
else
call wout('   phi_'//trim(i2s(xyzzyaaag17-1))//','//trim(i2s(xyzzyaaah&
&17-1))//','//trim(i2s(xyzzyaaai17-1))//','//trim(i2s(s))//'  (fixed) &
&:  '//trim(tmpr))
endif
endif
endif
enddo
enddo
enddo
xyzzyaaad17=0
do xyzzyaaai17=1,xyzzyaaav1(xyzzyaaaj17)
do xyzzyaaah17=1,xyzzyaaaw1(xyzzyaaaj17)
do xyzzyaaag17=1,xyzzyaaaw1(xyzzyaaaj17)
xyzzyaaad17=xyzzyaaad17+1
if(xyzzyaabk17(s,xyzzyaaad17+xyzzyaaat17))cycle
read(xyzzyaadp1,*,iostat=xyzzyaaab17)xyzzyaaay17,xyzzyaaak17
if(xyzzyaaab17/=0)then
backspace xyzzyaadp1
if(am_master)call wout('   Unspecified params   :  Zeroed')
exit read_params_phi
endif
empty_backflow=.false.
if(am_master.and.xyzzyaaak17/=0.and.xyzzyaaak17/=1)call errstop('INIT_&
&BACKFLOW','Optimizable flag should be 0 or 1 .')
xyzzyaabn1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)=xyzzyaaa&
&y17
xyzzyaabu1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)=xyzzyaaa&
&k17
if(am_master)then
tmpr=r2s(xyzzyaabn1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)&
&,'(f21.12)')
if(xyzzyaabn1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)<0.d0)&
&then
if(xyzzyaabu1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)==1)th&
&en
call wout('   theta_'//trim(i2s(xyzzyaaag17-1))//','//trim(i2s(xyzzyaa&
&ah17-1))//','//trim(i2s(xyzzyaaai17-1))//','//trim(i2s(s))//'  (opt) &
&: '//trim(tmpr))
else
call wout('   theta_'//trim(i2s(xyzzyaaag17-1))//','//trim(i2s(xyzzyaa&
&ah17-1))//','//trim(i2s(xyzzyaaai17-1))//','//trim(i2s(s))//' (fixd) &
&: '//trim(tmpr))
endif
else
if(xyzzyaabu1(xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,s,xyzzyaaaj17)==1)th&
&en
call wout('   theta_'//trim(i2s(xyzzyaaag17-1))//','//trim(i2s(xyzzyaa&
&ah17-1))//','//trim(i2s(xyzzyaaai17-1))//','//trim(i2s(s))//'  (opt) &
&:  '//trim(tmpr))
else
call wout('   theta_'//trim(i2s(xyzzyaaag17-1))//','//trim(i2s(xyzzyaa&
&ah17-1))//','//trim(i2s(xyzzyaaai17-1))//','//trim(i2s(s))//' (fixd) &
&:  '//trim(tmpr))
endif
endif
endif
enddo
enddo
enddo
enddo read_params_phi
read(xyzzyaadp1,'(a)',iostat=xyzzyaaab17)lineread
if(am_master.and.xyzzyaaab17/=0)call errstop('INIT_BACKFLOW','Can''t f&
&ind "END SET '//trim(i2s(xyzzyaaaj17))//'".')
if(trim(adjustl(lineread))/='END SET '//trim(i2s(xyzzyaaaj17)))call er&
&rstop('INIT_BACKFLOW','Can''t find "END SET '//trim(i2s(xyzzyaaaj17))&
&//'".')
enddo
if(am_master)then
xyzzyaaak17=count(xyzzyaaat1==0)
if(xyzzyaaak17==0)then
call wout('  Completeness of Phi   :  All atoms included')
else
if(allow_nochi_atoms)then
if(xyzzyaaak17==1)then
call wout('  Completeness of Phi   :  1 atom does not belong to any se&
&t')
else
call wout('  Completeness of Phi   :  '//trim(i2s(xyzzyaaak17))//' ato&
&ms do not belong to any set')
endif
else
call errstop('INIT_BACKFLOW','Some atoms do not belong to any Phi set.&
&')
endif
endif
endif
do
read(xyzzyaadp1,'(a)',err=100,end=101)lineread
if(trim(adjustl(lineread))=='END PHI TERM')exit
if(am_master.and.trim(adjustl(lineread))/='')call errstop('INIT_BACKFL&
&OW','String "END PHI TERM" not found.')
enddo
elseif(trim(adjustl(term_label))=='START AE CUTOFFS')then
if(xyzzyaabi17)call errstop_master('INIT_BACKFLOW',"Backflow block 'AE&
& CUTOFFS' found twice.")
if(.not.xyzzyaadb1)call errstop_master('INIT_BACKFLOW',"Found 'AE CUTO&
&FFS' block but there are no AE atoms in the system.")
if(am_master)then
call wout()
call wout(' AE cutoffs:')
endif
xyzzyaabi17=.true.
allocate(xyzzyaade1(nitot),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','in_set_ae_cut')
xyzzyaade1=0
read(xyzzyaadp1,*,err=100,end=101)
do xyzzyaaae17=1,nitot
if(.not.is_ae(xyzzyaaae17))cycle
read(xyzzyaadp1,*,iostat=xyzzyaaab17)xyzzyaaal17,xyzzyaaaj17,xyzzyaaay&
&17,xyzzyaaak17
if(xyzzyaaab17/=0)call errstop_master('INIT_BACKFLOW','Problem reading&
& AE CUTOFFS block in correlation.data. If you want to set the cutoffs&
& to their default values, remove the complete block, not just the dat&
&a.')
if(xyzzyaaal17/=xyzzyaaae17)call errstop_master('INIT_BACKFLOW','Seque&
&ntial order broken in AE CUTOFFS block.')
if(xyzzyaaaj17<1.or.xyzzyaaaj17>nitot)call errstop_master('INIT_BACKFL&
&OW','Set number out of range in AE CUTOFFS block.')
if(xyzzyaaak17<0.or.xyzzyaaak17>1)call errstop_master('INIT_BACKFLOW',&
&'Invalid "optimizable" flag in AE CUTOFFS block.')
xyzzyaade1(xyzzyaaae17)=xyzzyaaaj17
xyzzyaadf1(xyzzyaaae17)=xyzzyaaay17
xyzzyaaax17(xyzzyaaae17)=xyzzyaaak17
if(am_master)call wout('  Nucleus '//trim(i2s(xyzzyaaae17))//' in set &
&     :  '//trim(i2s(xyzzyaaaj17)))
enddo
xyzzyaadc1=maxval(xyzzyaade1)
if(am_master)then
call wout('  Total number of sets  :  '//trim(i2s(xyzzyaadc1)))
do xyzzyaaaj17=1,xyzzyaadc1
xyzzyaabj17=.false.
do xyzzyaaae17=1,nitot
if(.not.is_ae(xyzzyaaae17))cycle
if(xyzzyaade1(xyzzyaaae17)==xyzzyaaaj17)xyzzyaabj17=.true.
enddo
if(.not.xyzzyaabj17)call errstop('INIT_BACKFLOW','Set '//trim(i2s(xyzz&
&yaaaj17))//' empty in AE CUTOFFS block.')
enddo
endif
allocate(xyzzyaadj1(xyzzyaadc1),xyzzyaadd1(xyzzyaadc1),xyzzyaadl1(xyzz&
&yaadc1),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','ae_cutoff_set')
xyzzyaadj1=0.d0
xyzzyaadd1=1
xyzzyaadl1=0.d0
do xyzzyaaae17=1,nitot
if(.not.is_ae(xyzzyaaae17))cycle
xyzzyaaaj17=xyzzyaade1(xyzzyaaae17)
if(xyzzyaadk1(xyzzyaaae17)/=0.d0)then
if(xyzzyaadl1(xyzzyaaaj17)/=0.d0)then
xyzzyaadl1(xyzzyaaaj17)=min(xyzzyaadl1(xyzzyaaaj17),xyzzyaadk1(xyzzyaa&
&ae17))
else
xyzzyaadl1(xyzzyaaaj17)=xyzzyaadk1(xyzzyaaae17)
endif
endif
enddo
do xyzzyaaae17=1,nitot
if(.not.is_ae(xyzzyaaae17))cycle
xyzzyaaaj17=xyzzyaade1(xyzzyaaae17)
if(xyzzyaadf1(xyzzyaaae17)<=0.d0)then
if(xyzzyaadl1(xyzzyaaaj17)/=0.d0)then
xyzzyaadf1(xyzzyaaae17)=0.999d0*xyzzyaadl1(xyzzyaaaj17)
else
xyzzyaadf1(xyzzyaaae17)=1.d0
endif
endif
if(xyzzyaadl1(xyzzyaaaj17)/=0.d0.and.xyzzyaadf1(xyzzyaaae17)>xyzzyaadl&
&1(xyzzyaaaj17))then
tmpr=r2s(xyzzyaadl1(xyzzyaaaj17),'(f21.12)')
call errstop_master('INIT_BACKFLOW','Cutoff for nucleus no. '//trim(i2&
&s(xyzzyaaae17))//' exceeds set limit of '//trim(adjustl(tmpr))//' in &
&AE CUTOFFS block.')
endif
if(xyzzyaadj1(xyzzyaaaj17)/=0.d0)then
if(xyzzyaadj1(xyzzyaaaj17)/=xyzzyaadf1(xyzzyaaae17))call errstop_maste&
&r('INIT_BACKFLOW','The cutoffs for all nuclei belonging in the same s&
&et in AE CUTOFFS must be equal.')
if(xyzzyaadd1(xyzzyaaaj17)/=xyzzyaaax17(xyzzyaaae17))call errstop_mast&
&er('INIT_BACKFLOW','The "optimizable" flags for all nuclei belonging &
&in the same set in AE CUTOFFS must be equal.')
else
if(am_master)then
tmpr=r2s(xyzzyaadf1(xyzzyaaae17),'(f21.12)')
if(xyzzyaadl1(xyzzyaaaj17)/=0.d0)then
tmpr2=r2s(xyzzyaadl1(xyzzyaaaj17),'(f21.12)')
if(xyzzyaaax17(xyzzyaaae17)==1)then
call wout('  Cutoff_'//trim(i2s(xyzzyaaaj17))//' (optimizable):  '//tr&
&im(tmpr)//', limit: '//trim(tmpr2))
else
call wout('  Cutoff_'//trim(i2s(xyzzyaaaj17))//'       (fixed):  '//tr&
&im(tmpr)//', limit: '//trim(tmpr2))
endif
else
if(xyzzyaaax17(xyzzyaaae17)==1)then
call wout('  Cutoff_'//trim(i2s(xyzzyaaaj17))//' (optimizable):  '//tr&
&im(tmpr))
else
call wout('  Cutoff_'//trim(i2s(xyzzyaaaj17))//'       (fixed):  '//tr&
&im(tmpr))
endif
endif
endif
endif
xyzzyaadj1(xyzzyaaaj17)=xyzzyaadf1(xyzzyaaae17)
xyzzyaadd1(xyzzyaaaj17)=xyzzyaaax17(xyzzyaaae17)
enddo
do
read(xyzzyaadp1,'(a)',err=100,end=101)lineread
if(trim(adjustl(lineread))=='END AE CUTOFFS')exit
if(am_master.and.trim(adjustl(lineread))/='')call errstop('INIT_BACKFL&
&OW','String "END AE CUTOFFS" not found.')
enddo
elseif(trim(adjustl(term_label))=='END BACKFLOW')then
exit mainloop
elseif(trim(adjustl(term_label))/='')then
if(am_master)call errstop('INIT_BACKFLOW','Label not recognized.')
endif
enddo mainloop
do
read(xyzzyaadp1,'(a)',iostat=xyzzyaaab17)lineread
if(xyzzyaaab17<0)exit
if(xyzzyaaab17>0)call errstop('INIT_BACKFLOW','Problem reading correla&
&tion.data. Please check this file.')
if(trim(adjustl(lineread))=='START BACKFLOW')call errstop('INIT_BACKFL&
&OW','There seems to be more than one backflow block in correlation.da&
&ta.')
enddo
close(xyzzyaadp1)
if(am_master)then
call wout()
call wout('Finished reading backflow functions from correlation.data.'&
&)
call wout()
endif
xyzzyaabv1=0
xyzzyaabw1=0
xyzzyaabx1=0
xyzzyaaby1=0
xyzzyaaal1=0
xyzzyaaam1=0
if(xyzzyaaac1)then
allocate(xyzzyaado1(netot,netot,xyzzyaaai1),xyzzyaacc1(xyzzyaaai1),xyz&
&zyaacd1(xyzzyaaai1),xyzzyaace1(xyzzyaaai1),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','16')
xyzzyaabb17=xyzzyaabb17+dble(3*xyzzyaaai1+netot)
xyzzyaaba17=xyzzyaaba17+dble(netot*netot*xyzzyaaai1)
xyzzyaaam1=maxval(xyzzyaaaw1)
if(xyzzyaaab1)xyzzyaaam1=max(xyzzyaaam1,maxval(xyzzyaaap1))
allocate(xyzzyaacm1(xyzzyaaam1,nitot,netot),xyzzyaacn1(xyzzyaaam1,nito&
&t,netot),xyzzyaaco1(xyzzyaaam1,nitot,netot),xyzzyaaba1(nitot),xyzzyaa&
&ci1(nitot),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','17')
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaam1*nitot*(netot*3+2))
xyzzyaaba17=xyzzyaaba17+dble(nitot*2)
xyzzyaacm1=0.d0
xyzzyaacn1=0.d0
xyzzyaaco1=0.d0
do xyzzyaaae17=1,nitot
xyzzyaaba1(xyzzyaaae17)=0
if(xyzzyaaat1(xyzzyaaae17)/=0)xyzzyaaba1(xyzzyaaae17)=xyzzyaaaw1(xyzzy&
&aaat1(xyzzyaaae17))
if(xyzzyaaab1)then
if(xyzzyaaan1(xyzzyaaae17)/=0)xyzzyaaba1(xyzzyaaae17)=max(xyzzyaaba1(x&
&yzzyaaae17),xyzzyaaap1(xyzzyaaan1(xyzzyaaae17)))
endif
enddo
xyzzyaaak17=(netot*(netot-1))/2
xyzzyaaal1=maxval(xyzzyaaav1)
if(xyzzyaaaa1)xyzzyaaal1=max(xyzzyaaal1,xyzzyaaag1)
allocate(xyzzyaacj1(xyzzyaaal1,xyzzyaaak17),xyzzyaack1(xyzzyaaal1,xyzz&
&yaaak17),xyzzyaacl1(xyzzyaaal1,xyzzyaaak17),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','18')
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaak17*xyzzyaaal1*3)
xyzzyaacj1=0.d0
xyzzyaack1=0.d0
xyzzyaacl1=0.d0
endif
if(xyzzyaaab1)then
allocate(xyzzyaacf1(xyzzyaaah1),xyzzyaacg1(xyzzyaaah1),xyzzyaach1(xyzz&
&yaaah1),xyzzyaadn1(netot,xyzzyaaah1),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','19')
xyzzyaabb17=xyzzyaabb17+dble(3*xyzzyaaah1)
xyzzyaaba17=xyzzyaaba17+dble(xyzzyaaah1)
if(.not.xyzzyaaac1)then
xyzzyaaam1=maxval(xyzzyaaap1)
allocate(xyzzyaacm1(xyzzyaaam1,nitot,netot),xyzzyaacn1(xyzzyaaam1,nito&
&t,netot),xyzzyaaco1(xyzzyaaam1,nitot,netot),xyzzyaaba1(nitot),xyzzyaa&
&ci1(nitot),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','20')
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaam1*nitot*(netot*3+2))
xyzzyaacm1=0.d0
xyzzyaacn1=0.d0
xyzzyaaco1=0.d0
do xyzzyaaae17=1,nitot
xyzzyaaba1(xyzzyaaae17)=0
if(xyzzyaaan1(xyzzyaaae17)/=0)xyzzyaaba1(xyzzyaaae17)=xyzzyaaap1(xyzzy&
&aaan1(xyzzyaaae17))
enddo
endif
endif
if(xyzzyaaaa1)then
allocate(xyzzyaabz1(xyzzyaaak1),xyzzyaaca1(xyzzyaaak1),xyzzyaacb1(xyzz&
&yaaak1),xyzzyaadm1(netot,netot),stat=xyzzyaaaa17)
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaak1*7)
xyzzyaaba17=xyzzyaaba17+dble(netot*netot)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','21')
if(.not.xyzzyaaac1)then
xyzzyaaak17=(netot*(netot-1))/2
xyzzyaaal1=xyzzyaaag1
allocate(xyzzyaacj1(xyzzyaaal1,xyzzyaaak17),xyzzyaack1(xyzzyaaal1,xyzz&
&yaaak17),xyzzyaacl1(xyzzyaaal1,xyzzyaaak17),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'INIT_BACKFLOW','22')
xyzzyaabb17=xyzzyaabb17+dble(xyzzyaaak17*xyzzyaaal1*3)
endif
endif
call xyzzyaagv1
if(xyzzyaadb1)then
if(.not.xyzzyaabi17)then
xyzzyaadc1=count(is_ae)
allocate(xyzzyaadj1(xyzzyaadc1),xyzzyaadd1(xyzzyaadc1),xyzzyaadl1(xyzz&
&yaadc1),xyzzyaade1(nitot),stat=xyzzyaaaa17)
xyzzyaadj1=0.d0
xyzzyaadd1=1
xyzzyaadl1=0.d0
xyzzyaade1=0
xyzzyaaaf17=0
do xyzzyaaae17=1,nitot
if(.not.(is_ae(xyzzyaaae17)))cycle
xyzzyaaaf17=xyzzyaaaf17+1
xyzzyaade1(xyzzyaaae17)=xyzzyaaaf17
xyzzyaadj1(xyzzyaaaf17)=xyzzyaadf1(xyzzyaaae17)
xyzzyaadl1(xyzzyaaaf17)=xyzzyaadk1(xyzzyaaae17)
xyzzyaadd1(xyzzyaaaf17)=1
enddo
endif
do xyzzyaaae17=1,nitot
if(.not.(is_ae(xyzzyaaae17)))cycle
xyzzyaadg1(xyzzyaaae17)=1.d0/xyzzyaadf1(xyzzyaaae17)
xyzzyaadh1(xyzzyaaae17)=12*xyzzyaadg1(xyzzyaaae17)
xyzzyaadi1(xyzzyaaae17)=xyzzyaadh1(xyzzyaaae17)*xyzzyaadg1(xyzzyaaae17&
&)
enddo
endif
if(xyzzyaaaa1)then
call xyzzyaafi1
do s=1,xyzzyaaak1
xyzzyaacb1(s)=1.d0/xyzzyaabh1(s)
select case(xyzzyaaad1)
case(1)
xyzzyaabz1(s)=-xyzzyaacb1(s)
xyzzyaaca1(s)=0.d0
case(2)
xyzzyaabz1(s)=-2*xyzzyaacb1(s)
xyzzyaaca1(s)=2*xyzzyaacb1(s)**2
case(3)
xyzzyaabz1(s)=-3*xyzzyaacb1(s)
xyzzyaaca1(s)=6*xyzzyaacb1(s)**2
case default
xyzzyaabz1(s)=-xyzzyaaad1*xyzzyaacb1(s)
xyzzyaaca1(s)=xyzzyaaad1*xyzzyaaae1*xyzzyaacb1(s)**2
end select
enddo
endif
if(xyzzyaaab1)then
call xyzzyaafj1
do xyzzyaaaj17=1,xyzzyaaah1
xyzzyaach1(xyzzyaaaj17)=1.d0/xyzzyaabj1(xyzzyaaaj17)
select case(xyzzyaaad1)
case(1)
xyzzyaacf1(xyzzyaaaj17)=-xyzzyaach1(xyzzyaaaj17)
xyzzyaacg1(xyzzyaaaj17)=0.d0
case(2)
xyzzyaacf1(xyzzyaaaj17)=-2*xyzzyaach1(xyzzyaaaj17)
xyzzyaacg1(xyzzyaaaj17)=2*xyzzyaach1(xyzzyaaaj17)**2
case(3)
xyzzyaacf1(xyzzyaaaj17)=-3*xyzzyaach1(xyzzyaaaj17)
xyzzyaacg1(xyzzyaaaj17)=6*xyzzyaach1(xyzzyaaaj17)**2
case default
xyzzyaacf1(xyzzyaaaj17)=-xyzzyaaad1*xyzzyaach1(xyzzyaaaj17)
xyzzyaacg1(xyzzyaaaj17)=xyzzyaaad1*xyzzyaaae1*xyzzyaach1(xyzzyaaaj17)*&
&*2
end select
enddo
endif
if(xyzzyaaac1)then
call xyzzyaafk1
call xyzzyaafn1
do xyzzyaaaj17=1,xyzzyaaai1
xyzzyaace1(xyzzyaaaj17)=1.d0/xyzzyaabl1(xyzzyaaaj17)
select case(xyzzyaaad1)
case(1)
xyzzyaacc1(xyzzyaaaj17)=-xyzzyaace1(xyzzyaaaj17)
xyzzyaacd1(xyzzyaaaj17)=0.d0
case(2)
xyzzyaacc1(xyzzyaaaj17)=-2*xyzzyaace1(xyzzyaaaj17)
xyzzyaacd1(xyzzyaaaj17)=2*xyzzyaace1(xyzzyaaaj17)**2
case(3)
xyzzyaacc1(xyzzyaaaj17)=-3*xyzzyaace1(xyzzyaaaj17)
xyzzyaacd1(xyzzyaaaj17)=6*xyzzyaace1(xyzzyaaaj17)**2
case default
xyzzyaacc1(xyzzyaaaj17)=-xyzzyaaad1*xyzzyaace1(xyzzyaaaj17)
xyzzyaacd1(xyzzyaaaj17)=xyzzyaaad1*xyzzyaaae1*xyzzyaace1(xyzzyaaaj17)*&
&*2
end select
enddo
xyzzyaacp1=maxval(xyzzyaabl1(:))
endif
if(xyzzyaaaa1.or.xyzzyaaac1)then
xyzzyaacq1=0.d0
if(xyzzyaaaa1)xyzzyaacq1=maxval(xyzzyaabh1)
if(xyzzyaaac1)xyzzyaacq1=max(xyzzyaacq1,2*xyzzyaacp1)
endif
if(xyzzyaaab1.or.xyzzyaaac1)then
do xyzzyaaae17=1,nitot
xyzzyaaci1(xyzzyaaae17)=0.d0
if(xyzzyaaab1)then
if(xyzzyaaan1(xyzzyaaae17)/=0)xyzzyaaci1(xyzzyaaae17)=xyzzyaabj1(xyzzy&
&aaan1(xyzzyaaae17))
endif
if(xyzzyaaac1)then
if(xyzzyaaat1(xyzzyaaae17)/=0)xyzzyaaci1(xyzzyaaae17)=max(xyzzyaaci1(x&
&yzzyaaae17),xyzzyaabl1(xyzzyaaat1(xyzzyaaae17)))
endif
enddo
endif
if(xyzzyaaac1)call xyzzyaafp1
if(am_master)then
call wout('Allocated work arrays, using '//trim(byte2human(8.d0*xyzzya&
&abb17+4.d0*xyzzyaaba17))//'.')
if(xyzzyaaaa1.and..not.(xyzzyaaab1.or.xyzzyaaac1))then
call wout('Imposed e-e cusp conditions.')
elseif(xyzzyaaaa1.and.xyzzyaaab1.and..not.xyzzyaaac1)then
call wout('Imposed e-e and e-N cusp conditions.')
elseif(xyzzyaaac1)then
if(any(xyzzyaacv1(:)))then
call wout('Imposed e-e cusps, e-N cusps and irrotational-backflow cons&
&traints')
call wout('and checked them.')
else
call wout('Imposed e-e and e-N cusp conditions and checked them.')
endif
else
call wout('Imposed e-N cusp conditions.')
endif
if(xyzzyaads1)call wout('Plot of backflow function flagged.')
if(xyzzyaadb1)call wout('Will apply cut-offs around AE atoms.')
call wout()
call wout('Finished backflow setup.')
call wout()
endif
if(allocated(xyzzyaaax17))deallocate(xyzzyaaax17)
return
100 if(am_master)call errstop('INIT_BACKFLOW','Error reading correlati&
&on.data .')
101 if(am_master)call errstop('INIT_BACKFLOW','File correlation.data e&
&nded unexpectedly.')
end subroutine init_pbackflow
subroutine write_pbackflow(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xy&
&zzyaaaf18,xyzzyaaag18,xyzzyaaah18,xyzzyaaai18
logical xyzzyaaaj18
if(.not.am_master)return
inquire(file=trim(correlation_name),exist=xyzzyaaaj18)
if(xyzzyaaaj18)then
open(unit=xyzzyaadp1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaaa18)
else
open(unit=xyzzyaadp1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaaa18)
endif
if(xyzzyaaaa18/=0)call errstop('WRITE_BACKFLOW','Problem opening' //tr&
&im(correlation_name)//'.')
write(xyzzyaadp1,*)'START BACKFLOW'
write(xyzzyaadp1,*)'Title'
write(xyzzyaadp1,*)trim(adjustl(title))
write(xyzzyaadp1,*)'Truncation order'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaad1))
if(xyzzyaaaa1)then
write(xyzzyaadp1,*)'START ETA TERM'
write(xyzzyaadp1,*)'Expansion order'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaag1-1))
write(xyzzyaadp1,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud&
&)'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaaj1))
write(xyzzyaadp1,*)'Cut-off radii ;      Optimizable (0=NO; 1=YES; 2=Y&
&ES BUT NO SPIN-DEP)'
if(xyzzyaabo1(1)==2)then
write(xyzzyaadp1,*)xyzzyaabh1(1),xyzzyaabo1(1),'      ! L_1'
else
do xyzzyaaae18=1,xyzzyaaak1
write(xyzzyaadp1,*)xyzzyaabh1(xyzzyaaae18),xyzzyaabo1(xyzzyaaae18),'      !&
& L_',trim(i2s(xyzzyaaae18))
enddo
endif
write(xyzzyaadp1,*)'Parameter ;          Optimizable (0=NO; 1=YES)'
do xyzzyaaae18=1,xyzzyaaak1
do xyzzyaaaf18=1,xyzzyaaag1
if(xyzzyaabp1(xyzzyaaaf18,xyzzyaaae18)/=-1)then
write(xyzzyaadp1,*)xyzzyaabi1(xyzzyaaaf18,xyzzyaaae18),xyzzyaabp1(xyzz&
&yaaaf18,xyzzyaaae18),'      ! c_',trim(i2s(xyzzyaaaf18-1)),',',trim(i&
&2s(xyzzyaaae18))
endif
enddo
enddo
write(xyzzyaadp1,*)'END ETA TERM'
endif
if(xyzzyaaab1)then
write(xyzzyaadp1,*)'START MU TERM'
write(xyzzyaadp1,*)'Number of sets ; labelling (1->atom in s. cell; 2-&
&>atom in p. cell; 3->species)'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaah1))//' '//trim(i2s(xyzzyaa&
&bb1))
do xyzzyaaag18=1,xyzzyaaah1
write(xyzzyaadp1,*)'START SET ',trim(i2s(xyzzyaaag18))
if(xyzzyaabb1==3)then
write(xyzzyaadp1,*)'Number of species in set'
else
write(xyzzyaadp1,*)'Number of atoms in set'
endif
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaabf1(xyzzyaaag18)))
if(xyzzyaabf1(xyzzyaaag18)==1)then
if(xyzzyaabb1==3)then
write(xyzzyaadp1,*)'Label of the species in this set'
else
write(xyzzyaadp1,*)'Label of the atom in this set'
endif
else
if(xyzzyaabb1==3)then
write(xyzzyaadp1,*)'Labels of the species in this set'
else
write(xyzzyaadp1,*)'Labels of the atoms in this set'
endif
endif
call write_list_int(xyzzyaabf1(xyzzyaaag18),              xyzzyaabd1(1&
&:xyzzyaabf1(xyzzyaaag18),xyzzyaaag18),10,4,1,xyzzyaadp1)
xyzzyaaah18=0
if(xyzzyaact1(xyzzyaaag18))xyzzyaaah18=1
write(xyzzyaadp1,*)'Type of e-N cusp conditions (0->PP/cuspless AE; 1-&
&>AE with cusp)'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaah18))
write(xyzzyaadp1,*)'Expansion order'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaap1(xyzzyaaag18)-1))
write(xyzzyaadp1,*)'Spin dep (0->u=d; 1->u/=d)'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaaq1(xyzzyaaag18)))
write(xyzzyaadp1,*)'Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)'
write(xyzzyaadp1,*)xyzzyaabj1(xyzzyaaag18),xyzzyaabq1(xyzzyaaag18)
write(xyzzyaadp1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaae18=1,xyzzyaaas1(xyzzyaaag18)
do xyzzyaaaf18=1,xyzzyaaap1(xyzzyaaag18)
if(xyzzyaabr1(xyzzyaaaf18,xyzzyaaae18,xyzzyaaag18)/=-1)then
write(xyzzyaadp1,*)xyzzyaabk1(xyzzyaaaf18,xyzzyaaae18,xyzzyaaag18),xyz&
&zyaabr1(xyzzyaaaf18,xyzzyaaae18,xyzzyaaag18),'      ! mu_',trim(i2s(x&
&yzzyaaaf18-1)),',',trim(i2s(xyzzyaaae18))
endif
enddo
enddo
write(xyzzyaadp1,*)'END SET ',trim(i2s(xyzzyaaag18))
enddo
write(xyzzyaadp1,*)'END MU TERM'
endif
if(xyzzyaaac1)then
write(xyzzyaadp1,*)'START PHI TERM'
write(xyzzyaadp1,*)'Number of sets ; labelling (1->atom in s. cell; 2-&
&>atom in p. cell; 3->species)'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaai1))//' '//trim(i2s(xyzzyaa&
&bc1))
do xyzzyaaag18=1,xyzzyaaai1
write(xyzzyaadp1,*)'START SET ',trim(i2s(xyzzyaaag18))
if(xyzzyaabc1==3)then
write(xyzzyaadp1,*)'Number of species in set'
else
write(xyzzyaadp1,*)'Number of atoms in set'
endif
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaabg1(xyzzyaaag18)))
if(xyzzyaabg1(xyzzyaaag18)==1)then
if(xyzzyaabc1==3)then
write(xyzzyaadp1,*)'Label of the species in this set'
else
write(xyzzyaadp1,*)'Label of the atom in this set'
endif
else
if(xyzzyaabc1==3)then
write(xyzzyaadp1,*)'Labels of the species in this set'
else
write(xyzzyaadp1,*)'Labels of the atoms in this set'
endif
endif
call write_list_int(xyzzyaabg1(xyzzyaaag18),xyzzyaabe1(1:xyzzyaabg1(xy&
&zzyaaag18),xyzzyaaag18),10,4,1,xyzzyaadp1)
xyzzyaaah18=0
if(xyzzyaacu1(xyzzyaaag18))xyzzyaaah18=1
write(xyzzyaadp1,*)'Type of e-N cusp conditions (0=PP; 1=AE)'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaah18))
xyzzyaaah18=0
if(xyzzyaacv1(xyzzyaaag18))xyzzyaaah18=1
write(xyzzyaadp1,*)'Irrotational Phi term (0=NO; 1=YES)'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaah18))
write(xyzzyaadp1,*)'Electron-nucleus expansion order N_eN'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaaw1(xyzzyaaag18)-1))
write(xyzzyaadp1,*)'Electron-electron expansion order N_ee'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaav1(xyzzyaaag18)-1))
write(xyzzyaadp1,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud&
&)'
write(xyzzyaadp1,'(3x,a)')trim(i2s(xyzzyaaay1(xyzzyaaag18)))
write(xyzzyaadp1,*)'Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)'
write(xyzzyaadp1,*)xyzzyaabl1(xyzzyaaag18),xyzzyaabs1(xyzzyaaag18)
write(xyzzyaadp1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaae18=1,xyzzyaaaz1(xyzzyaaag18)
do xyzzyaaad18=1,xyzzyaaav1(xyzzyaaag18)
do xyzzyaaac18=1,xyzzyaaaw1(xyzzyaaag18)
do xyzzyaaab18=1,xyzzyaaaw1(xyzzyaaag18)
if(xyzzyaabt1(xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xyzzyaaa&
&g18)/=-1)then
write(xyzzyaadp1,*)xyzzyaabm1(xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzz&
&yaaae18,xyzzyaaag18),xyzzyaabt1(xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,x&
&yzzyaaae18,xyzzyaaag18),'      ! phi_',trim(i2s(xyzzyaaab18-1)),',',t&
&rim(i2s(xyzzyaaac18-1)),',',trim(i2s(xyzzyaaad18-1)),',',trim(i2s(xyz&
&zyaaae18))
endif
enddo
enddo
enddo
do xyzzyaaad18=1,xyzzyaaav1(xyzzyaaag18)
do xyzzyaaac18=1,xyzzyaaaw1(xyzzyaaag18)
do xyzzyaaab18=1,xyzzyaaaw1(xyzzyaaag18)
if(xyzzyaabu1(xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xyzzyaaa&
&g18)/=-1)then
write(xyzzyaadp1,*)xyzzyaabn1(xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzz&
&yaaae18,xyzzyaaag18),xyzzyaabu1(xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,x&
&yzzyaaae18,xyzzyaaag18),'      ! theta_',trim(i2s(xyzzyaaab18-1)),','&
&,trim(i2s(xyzzyaaac18-1)),',',trim(i2s(xyzzyaaad18-1)),',',trim(i2s(x&
&yzzyaaae18))
endif
enddo
enddo
enddo
enddo
write(xyzzyaadp1,*)'END SET ',trim(i2s(xyzzyaaag18))
enddo
write(xyzzyaadp1,*)'END PHI TERM'
endif
if(xyzzyaadb1)then
write(xyzzyaadp1,*)'START AE CUTOFFS'
write(xyzzyaadp1,*)'Nucleus ; Set ; Cutoff length     ;  Optimizable (&
&0=NO; 1=YES)'
do xyzzyaaai18=1,nitot
if(.not.is_ae(xyzzyaaai18))cycle
xyzzyaaag18=xyzzyaade1(xyzzyaaai18)
write(xyzzyaadp1,*)' ',trim(i2s(xyzzyaaai18)),'         ',trim(i2s(xyz&
&zyaaag18)),'   ',xyzzyaadj1(xyzzyaaag18),xyzzyaadd1(xyzzyaaag18)
enddo
write(xyzzyaadp1,*)'END AE CUTOFFS'
endif
write(xyzzyaadp1,*)'END BACKFLOW'
write(xyzzyaadp1,*)
close(xyzzyaadp1)
end subroutine write_pbackflow
subroutine setup_pbf_params(nparam)
implicit none
integer,intent(out) :: nparam
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,xy&
&zzyaaaf19,xyzzyaaag19,xyzzyaaah19,xyzzyaaai19,xyzzyaaaj19
nparam=0
xyzzyaady1=0
if(.not.opt_backflow)return
if(xyzzyaaaa1)then
xyzzyaaaa19=0
do xyzzyaaab19=1,xyzzyaaak1
if(xyzzyaabo1(xyzzyaaab19)==1.or.(xyzzyaabo1(xyzzyaaab19)==2.and.xyzzy&
&aaab19==1))xyzzyaaaa19=xyzzyaaaa19+1
do xyzzyaaac19=1,xyzzyaaag1
if(xyzzyaabp1(xyzzyaaac19,xyzzyaaab19)==1)xyzzyaaaa19=xyzzyaaaa19+1
enddo
enddo
if(xyzzyaaaa19==0)xyzzyaaaa19=-1
xyzzyaady1(1)=xyzzyaaaa19
endif
if(xyzzyaaab1)then
xyzzyaaaa19=0
do xyzzyaaag19=1,xyzzyaaah1
if(xyzzyaabq1(xyzzyaaag19)==1)xyzzyaaaa19=xyzzyaaaa19+1
do xyzzyaaab19=1,xyzzyaaas1(xyzzyaaag19)
do xyzzyaaac19=1,xyzzyaaap1(xyzzyaaag19)
if(xyzzyaabr1(xyzzyaaac19,xyzzyaaab19,xyzzyaaag19)==1)xyzzyaaaa19=xyzz&
&yaaaa19+1
enddo
enddo
enddo
if(xyzzyaaaa19==0)xyzzyaaaa19=-1
xyzzyaady1(2)=xyzzyaaaa19
endif
if(xyzzyaaac1)then
xyzzyaaaa19=0
do xyzzyaaag19=1,xyzzyaaai1
if(xyzzyaabs1(xyzzyaaag19)==1)xyzzyaaaa19=xyzzyaaaa19+1
do xyzzyaaab19=1,xyzzyaaaz1(xyzzyaaag19)
do xyzzyaaaf19=1,xyzzyaaav1(xyzzyaaag19)
do xyzzyaaae19=1,xyzzyaaaw1(xyzzyaaag19)
do xyzzyaaad19=1,xyzzyaaaw1(xyzzyaaag19)
if(xyzzyaabt1(xyzzyaaad19,xyzzyaaae19,xyzzyaaaf19,xyzzyaaab19,xyzzyaaa&
&g19)==1)xyzzyaaaa19=xyzzyaaaa19+1
enddo
enddo
enddo
do xyzzyaaaf19=1,xyzzyaaav1(xyzzyaaag19)
do xyzzyaaae19=1,xyzzyaaaw1(xyzzyaaag19)
do xyzzyaaad19=1,xyzzyaaaw1(xyzzyaaag19)
if(xyzzyaabu1(xyzzyaaad19,xyzzyaaae19,xyzzyaaaf19,xyzzyaaab19,xyzzyaaa&
&g19)==1)xyzzyaaaa19=xyzzyaaaa19+1
enddo
enddo
enddo
enddo
enddo
if(xyzzyaaaa19==0)xyzzyaaaa19=-1
xyzzyaady1(3)=xyzzyaaaa19
endif
if(xyzzyaadb1)then
xyzzyaaaa19=0
do xyzzyaaag19=1,xyzzyaadc1
if(xyzzyaadd1(xyzzyaaag19)==1)xyzzyaaaa19=xyzzyaaaa19+1
enddo
if(xyzzyaaaa19==0)xyzzyaaaa19=-1
xyzzyaady1(4)=xyzzyaaaa19
endif
xyzzyaadx1=sum(xyzzyaady1,xyzzyaady1>0)
nparam=xyzzyaadx1
allocate(xyzzyaadw1(xyzzyaadx1),stat=xyzzyaaah19)
call check_alloc(xyzzyaaah19,'SETUP_BF_PARAMS','bf_param_sec')
xyzzyaaaj19=0
do xyzzyaaaa19=1,xyzzyaadv1
if(xyzzyaady1(xyzzyaaaa19)<1)cycle
xyzzyaaai19=xyzzyaaaj19+1
xyzzyaaaj19=xyzzyaaaj19+xyzzyaady1(xyzzyaaaa19)
xyzzyaadw1(xyzzyaaai19:xyzzyaaaj19)=xyzzyaaaa19
enddo
call xyzzyaafe1
end subroutine setup_pbf_params
subroutine finish_pbf_params
implicit none
if(.not.opt_backflow)return
deallocate(xyzzyaadw1)
call xyzzyaaff1
end subroutine finish_pbf_params
subroutine xyzzyaafe1
implicit none
integer xyzzyaaaa21,xyzzyaaab21
if(xyzzyaady1(1)>0)then
xyzzyaaaa21=xyzzyaady1(1)
allocate(xyzzyaadz1(xyzzyaaak1,0:xyzzyaaaa21),xyzzyaaea1(xyzzyaaag1,xy&
&zzyaaak1,0:xyzzyaaaa21),xyzzyaaeb1(xyzzyaaak1,0:xyzzyaaaa21),xyzzyaae&
&c1(xyzzyaaak1,0:xyzzyaaaa21),xyzzyaaed1(xyzzyaaak1,0:xyzzyaaaa21),sta&
&t=xyzzyaaab21)
call check_alloc(xyzzyaaab21,'SETUP_BF_PBUFFER','eta')
xyzzyaadz1=0.d0
xyzzyaaea1=0.d0
xyzzyaaeb1=0.d0
xyzzyaaec1=0.d0
xyzzyaaed1=0.d0
xyzzyaaee1=xyzzyaaag1*xyzzyaaak1
endif
if(xyzzyaady1(2)>0)then
xyzzyaaaa21=xyzzyaady1(2)
allocate(xyzzyaaef1(xyzzyaaah1,0:xyzzyaaaa21),xyzzyaaeg1(maxval(xyzzya&
&aap1),maxval(xyzzyaaas1),xyzzyaaah1,0:xyzzyaaaa21),xyzzyaaeh1(xyzzyaa&
&ah1,0:xyzzyaaaa21),xyzzyaaei1(xyzzyaaah1,0:xyzzyaaaa21),xyzzyaaej1(xy&
&zzyaaah1,0:xyzzyaaaa21),stat=xyzzyaaab21)
call check_alloc(xyzzyaaab21,'SETUP_BF_PBUFFER','mu')
xyzzyaaef1=0.d0
xyzzyaaeg1=0.d0
xyzzyaaeh1=0.d0
xyzzyaaei1=0.d0
xyzzyaaej1=0.d0
xyzzyaaek1=maxval(xyzzyaaap1)*maxval(xyzzyaaas1)*xyzzyaaah1
endif
if(xyzzyaady1(3)>0)then
xyzzyaaaa21=xyzzyaady1(3)
allocate(xyzzyaael1(xyzzyaaai1,0:xyzzyaaaa21),xyzzyaaem1(xyzzyaaai1,0:&
&xyzzyaaaa21),xyzzyaaen1(xyzzyaaai1,0:xyzzyaaaa21),xyzzyaaeo1(xyzzyaaa&
&i1,0:xyzzyaaaa21),xyzzyaaep1(maxval(xyzzyaaaw1),maxval(xyzzyaaaw1),ma&
&xval(xyzzyaaav1),maxval(xyzzyaaaz1),xyzzyaaai1,0:xyzzyaaaa21),xyzzyaa&
&eq1(maxval(xyzzyaaaw1),maxval(xyzzyaaaw1),maxval(xyzzyaaav1),maxval(x&
&yzzyaaaz1),xyzzyaaai1,0:xyzzyaaaa21),xyzzyaaer1(0:xyzzyaaaa21),stat=x&
&yzzyaaab21)
call check_alloc(xyzzyaaab21,'SETUP_BF_PBUFFER','phi')
xyzzyaael1=0.d0
xyzzyaaem1=0.d0
xyzzyaaen1=0.d0
xyzzyaaeo1=0.d0
xyzzyaaep1=0.d0
xyzzyaaeq1=0.d0
xyzzyaaer1=0.d0
xyzzyaaes1=(maxval(xyzzyaaaw1)**2)*maxval(xyzzyaaav1)*maxval(xyzzyaaaz&
&1)*xyzzyaaai1
endif
if(xyzzyaady1(4)>0)then
xyzzyaaaa21=xyzzyaady1(4)
allocate(xyzzyaaet1(xyzzyaadc1,0:xyzzyaaaa21),xyzzyaaeu1(nitot,0:xyzzy&
&aaaa21),xyzzyaaev1(nitot,0:xyzzyaaaa21),xyzzyaaew1(nitot,0:xyzzyaaaa2&
&1),xyzzyaaex1(nitot,0:xyzzyaaaa21),stat=xyzzyaaab21)
call check_alloc(xyzzyaaab21,'SETUP_BF_PBUFFER','ae')
xyzzyaaet1=0.d0
xyzzyaaeu1=0.d0
xyzzyaaev1=0.d0
xyzzyaaew1=0.d0
xyzzyaaex1=0.d0
endif
end subroutine xyzzyaafe1
subroutine xyzzyaaff1
implicit none
if(xyzzyaady1(1)>0)deallocate(xyzzyaadz1,xyzzyaaea1,xyzzyaaeb1,xyzzyaa&
&ec1,xyzzyaaed1)
if(xyzzyaady1(2)>0)deallocate(xyzzyaaef1,xyzzyaaeg1,xyzzyaaeh1,xyzzyaa&
&ei1,xyzzyaaej1)
if(xyzzyaady1(3)>0)deallocate(xyzzyaael1,xyzzyaaem1,xyzzyaaen1,xyzzyaa&
&eo1,xyzzyaaep1,xyzzyaaeq1,xyzzyaaer1)
if(xyzzyaady1(4)>0)deallocate(xyzzyaaet1,xyzzyaaeu1,xyzzyaaev1,xyzzyaa&
&ew1,xyzzyaaex1)
end subroutine xyzzyaaff1
subroutine get_pbf_params(params,has_lolim,lolim,has_hilim,hilim,is_sh&
&allow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,label&
&)
implicit none
real(dp),intent(inout) :: params(xyzzyaadx1),lolim(xyzzyaadx1),hilim(x&
&yzzyaadx1)
logical,intent(inout) :: has_lolim(xyzzyaadx1),has_hilim(xyzzyaadx1),i&
&s_shallow(xyzzyaadx1),is_redundant(xyzzyaadx1),is_linear(xyzzyaadx1),&
&is_loglinear(xyzzyaadx1),has_aderiv(xyzzyaadx1),affect_map(xyzzyaadx1&
&,xyzzyaadx1)
character(2),intent(inout) :: label(xyzzyaadx1)
integer xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyzzyaaad23,xyzzyaaae23,xy&
&zzyaaaf23,xyzzyaaag23,xyzzyaaah23,xyzzyaaai23,xyzzyaaaj23,xyzzyaaak23
real(dp) xyzzyaaal23,xyzzyaaam23
logical xyzzyaaan23
has_lolim=.false.
lolim=0.d0
has_hilim=.false.
hilim=0.d0
is_shallow=.false.
is_redundant=.false.
is_linear=.false.
is_loglinear=.false.
has_aderiv=.false.
affect_map=.false.
do xyzzyaaaf23=1,xyzzyaadx1
affect_map(xyzzyaaaf23,xyzzyaaaf23)=.true.
enddo
xyzzyaaal23=1.1d-8
if(isperiodic)xyzzyaaam23=0.999999d0*wigner_seitz_radius
xyzzyaaaf23=0
xyzzyaaak23=0
if(xyzzyaady1(1)>0)then
xyzzyaaaj23=xyzzyaaak23+1
xyzzyaaak23=xyzzyaaak23+xyzzyaady1(1)
label(xyzzyaaaj23:xyzzyaaak23)='Be'
xyzzyaaah23=0
do xyzzyaaad23=1,xyzzyaaak1
if(xyzzyaabo1(xyzzyaaad23)==1.or.(xyzzyaabo1(xyzzyaaad23)==2.and.xyzzy&
&aaad23==1))then
xyzzyaaaf23=xyzzyaaaf23+1
params(xyzzyaaaf23)=xyzzyaabh1(xyzzyaaad23)
has_lolim(xyzzyaaaf23)=.true.
lolim(xyzzyaaaf23)=xyzzyaaal23
if(isperiodic)then
has_hilim(xyzzyaaaf23)=.true.
hilim(xyzzyaaaf23)=xyzzyaaam23
endif
is_shallow(xyzzyaaaf23)=.true.
if(xyzzyaabo1(xyzzyaaad23)==1)then
is_redundant(xyzzyaaaf23)=all(xyzzyaabi1(:,xyzzyaaad23)==0.d0)
else
is_redundant(xyzzyaaaf23)=all(xyzzyaabi1(:,:)==0.d0)
endif
xyzzyaaah23=xyzzyaaaf23
elseif(xyzzyaabo1(xyzzyaaad23)==0)then
xyzzyaaah23=0
endif
do xyzzyaaae23=1,xyzzyaaag1
if(xyzzyaabp1(xyzzyaaae23,xyzzyaaad23)==1)then
xyzzyaaaf23=xyzzyaaaf23+1
params(xyzzyaaaf23)=xyzzyaabi1(xyzzyaaae23,xyzzyaaad23)*xyzzyaabh1(xyz&
&zyaaad23)**(xyzzyaaae23-1)
if(xyzzyaaah23>0)affect_map(xyzzyaaaf23,xyzzyaaah23)=.true.
endif
enddo
enddo
if(xyzzyaaak23/=xyzzyaaaf23)call errstop('GET_BF_PARAMS','Bad param co&
&unt (eta).')
endif
if(xyzzyaady1(2)>0)then
xyzzyaaaj23=xyzzyaaak23+1
xyzzyaaak23=xyzzyaaak23+xyzzyaady1(2)
label(xyzzyaaaj23:xyzzyaaak23)='Bm'
do xyzzyaaag23=1,xyzzyaaah1
xyzzyaaah23=0
if(xyzzyaabq1(xyzzyaaag23)==1)then
xyzzyaaaf23=xyzzyaaaf23+1
params(xyzzyaaaf23)=xyzzyaabj1(xyzzyaaag23)
has_lolim(xyzzyaaaf23)=.true.
lolim(xyzzyaaaf23)=xyzzyaaal23
if(isperiodic)then
has_hilim(xyzzyaaaf23)=.true.
hilim(xyzzyaaaf23)=xyzzyaaam23
endif
is_shallow(xyzzyaaaf23)=.true.
is_redundant(xyzzyaaaf23)=all(xyzzyaabk1(:,:,xyzzyaaag23)==0.d0)
xyzzyaaah23=xyzzyaaaf23
endif
do xyzzyaaad23=1,xyzzyaaas1(xyzzyaaag23)
do xyzzyaaae23=1,xyzzyaaap1(xyzzyaaag23)
if(xyzzyaabr1(xyzzyaaae23,xyzzyaaad23,xyzzyaaag23)==1)then
xyzzyaaaf23=xyzzyaaaf23+1
params(xyzzyaaaf23)=xyzzyaabk1(xyzzyaaae23,xyzzyaaad23,xyzzyaaag23)*xy&
&zzyaabj1(xyzzyaaag23)**(xyzzyaaae23-1)
if(xyzzyaaah23>0)affect_map(xyzzyaaaf23,xyzzyaaah23)=.true.
endif
enddo
enddo
enddo
if(xyzzyaaak23/=xyzzyaaaf23)call errstop('GET_BF_PARAMS','Bad param co&
&unt (mu).')
endif
if(xyzzyaady1(3)>0)then
xyzzyaaaj23=xyzzyaaak23+1
xyzzyaaak23=xyzzyaaak23+xyzzyaady1(3)
label(xyzzyaaaj23:xyzzyaaak23)='BP'
do xyzzyaaag23=1,xyzzyaaai1
xyzzyaaah23=0
if(xyzzyaabs1(xyzzyaaag23)==1)then
xyzzyaaaf23=xyzzyaaaf23+1
params(xyzzyaaaf23)=xyzzyaabl1(xyzzyaaag23)
has_lolim(xyzzyaaaf23)=.true.
lolim(xyzzyaaaf23)=xyzzyaaal23
if(isperiodic)then
has_hilim(xyzzyaaaf23)=.true.
hilim(xyzzyaaaf23)=.5d0*xyzzyaaam23
endif
is_shallow(xyzzyaaaf23)=.true.
is_redundant(xyzzyaaaf23)=all(xyzzyaabm1(:,:,:,:,xyzzyaaag23)==0.d0).a&
&nd.all(xyzzyaabn1(:,:,:,:,xyzzyaaag23)==0.d0)
xyzzyaaah23=xyzzyaaaf23
endif
do xyzzyaaad23=1,xyzzyaaaz1(xyzzyaaag23)
do xyzzyaaac23=1,xyzzyaaav1(xyzzyaaag23)
do xyzzyaaab23=1,xyzzyaaaw1(xyzzyaaag23)
do xyzzyaaaa23=1,xyzzyaaaw1(xyzzyaaag23)
if(xyzzyaabt1(xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyzzyaaad23,xyzzyaaa&
&g23)==1)then
xyzzyaaaf23=xyzzyaaaf23+1
params(xyzzyaaaf23)=xyzzyaabm1(xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyz&
&zyaaad23,xyzzyaaag23)*xyzzyaabl1(xyzzyaaag23)**(xyzzyaaaa23+xyzzyaaab&
&23+xyzzyaaac23-3)
if(xyzzyaaah23>0)affect_map(xyzzyaaaf23,xyzzyaaah23)=.true.
endif
enddo
enddo
enddo
do xyzzyaaac23=1,xyzzyaaav1(xyzzyaaag23)
do xyzzyaaab23=1,xyzzyaaaw1(xyzzyaaag23)
do xyzzyaaaa23=1,xyzzyaaaw1(xyzzyaaag23)
if(xyzzyaabu1(xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyzzyaaad23,xyzzyaaa&
&g23)==1)then
xyzzyaaaf23=xyzzyaaaf23+1
params(xyzzyaaaf23)=xyzzyaabn1(xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyz&
&zyaaad23,xyzzyaaag23)*xyzzyaabl1(xyzzyaaag23)**(xyzzyaaaa23+xyzzyaaab&
&23+xyzzyaaac23-3)
if(xyzzyaaah23>0)affect_map(xyzzyaaaf23,xyzzyaaah23)=.true.
endif
enddo
enddo
enddo
enddo
enddo
if(xyzzyaaak23/=xyzzyaaaf23)call errstop('GET_BF_PARAMS','Bad param co&
&unt (Phi+Theta).')
endif
if(xyzzyaady1(4)>0)then
xyzzyaaaj23=xyzzyaaak23+1
xyzzyaaak23=xyzzyaaak23+xyzzyaady1(4)
label(xyzzyaaaj23:xyzzyaaak23)='BA'
do xyzzyaaag23=1,xyzzyaadc1
if(xyzzyaadd1(xyzzyaaag23)==1)then
xyzzyaaaf23=xyzzyaaaf23+1
params(xyzzyaaaf23)=xyzzyaadj1(xyzzyaaag23)
has_lolim(xyzzyaaaf23)=.true.
lolim(xyzzyaaaf23)=xyzzyaaal23
if(xyzzyaadl1(xyzzyaaag23)>0.d0)then
has_hilim(xyzzyaaaf23)=.true.
hilim(xyzzyaaaf23)=xyzzyaadl1(xyzzyaaag23)
endif
is_shallow(xyzzyaaaf23)=.true.
xyzzyaaan23=.true.
if(xyzzyaaaa1)xyzzyaaan23=all(xyzzyaabi1(:,:)==0.d0)
if(xyzzyaaan23)then
do xyzzyaaai23=1,nitot
if(xyzzyaade1(xyzzyaaai23)==xyzzyaaag23.and.count(xyzzyaade1==xyzzyaaa&
&g23)==1)cycle
if(xyzzyaaab1)then
if(xyzzyaaan1(xyzzyaaai23)/=0)xyzzyaaan23=xyzzyaaan23.and.all(xyzzyaab&
&k1(:,:,xyzzyaaan1(xyzzyaaai23))==0.d0)
endif
if(xyzzyaaac1)then
if(xyzzyaaat1(xyzzyaaai23)/=0)xyzzyaaan23=xyzzyaaan23.and.all(xyzzyaab&
&m1(:,:,:,:,xyzzyaaat1(xyzzyaaai23))==0.d0).and.all(xyzzyaabn1(:,:,:,:&
&,xyzzyaaat1(xyzzyaaai23))==0.d0)
endif
if(.not.xyzzyaaan23)exit
enddo
endif
is_redundant(xyzzyaaaf23)=xyzzyaaan23
endif
enddo
if(xyzzyaaak23/=xyzzyaaaf23)call errstop('GET_BF_PARAMS','Bad param co&
&unt (AE cutoff).')
endif
end subroutine get_pbf_params
subroutine put_pbf_params(params,ignore,iparam_buffer,prestore,bad_par&
&ams)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaadx1)
logical,intent(in) :: ignore(xyzzyaadx1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaae24,xy&
&zzyaaaf24,xyzzyaaag24,xyzzyaaah24,xyzzyaaai24,xyzzyaaaj24,xyzzyaaak24
real(dp),allocatable :: xyzzyaaal24(:,:)
bad_params=.false.
if(prestore)then
call xyzzyaafh1(iparam_buffer)
return
endif
xyzzyaaaj24=0
if(xyzzyaady1(1)>0)then
xyzzyaaaf24=xyzzyaaaj24
xyzzyaaai24=xyzzyaaaj24+1
xyzzyaaaj24=xyzzyaaaj24+xyzzyaady1(1)
if(any(.not.ignore(xyzzyaaai24:xyzzyaaaj24)))then
do xyzzyaaad24=1,xyzzyaaak1
if(xyzzyaabo1(xyzzyaaad24)==1.or.(xyzzyaabo1(xyzzyaaad24)==2.and.xyzzy&
&aaad24==1))then
xyzzyaaaf24=xyzzyaaaf24+1
if(.not.ignore(xyzzyaaaf24))then
xyzzyaabh1(xyzzyaaad24)=params(xyzzyaaaf24)
if(xyzzyaabo1(xyzzyaaad24)==2)xyzzyaabh1(:)=xyzzyaabh1(1)
endif
endif
do xyzzyaaae24=1,xyzzyaaag1
if(xyzzyaabp1(xyzzyaaae24,xyzzyaaad24)==1)then
xyzzyaaaf24=xyzzyaaaf24+1
if(.not.ignore(xyzzyaaaf24))then
params(xyzzyaaaf24)=params(xyzzyaaaf24)/xyzzyaabh1(xyzzyaaad24)**(xyzz&
&yaaae24-1)
xyzzyaabi1(xyzzyaaae24,xyzzyaaad24)=params(xyzzyaaaf24)
endif
endif
enddo
enddo
call xyzzyaafi1
do xyzzyaaad24=1,xyzzyaaak1
if(xyzzyaabo1(xyzzyaaad24)>0)then
xyzzyaacb1(xyzzyaaad24)=1.d0/xyzzyaabh1(xyzzyaaad24)
select case(xyzzyaaad1)
case(1)
xyzzyaabz1(xyzzyaaad24)=-xyzzyaacb1(xyzzyaaad24)
xyzzyaaca1(xyzzyaaad24)=0.d0
case(2)
xyzzyaabz1(xyzzyaaad24)=-2*xyzzyaacb1(xyzzyaaad24)
xyzzyaaca1(xyzzyaaad24)=2*xyzzyaacb1(xyzzyaaad24)**2
case(3)
xyzzyaabz1(xyzzyaaad24)=-3*xyzzyaacb1(xyzzyaaad24)
xyzzyaaca1(xyzzyaaad24)=6*xyzzyaacb1(xyzzyaaad24)**2
case default
xyzzyaabz1(xyzzyaaad24)=-xyzzyaaad1*xyzzyaacb1(xyzzyaaad24)
xyzzyaaca1(xyzzyaaad24)=xyzzyaaad1*xyzzyaaae1*xyzzyaacb1(xyzzyaaad24)*&
&*2
end select
endif
enddo
if(xyzzyaaaj24/=xyzzyaaaf24)call errstop('PUT_BF_PARAMS','Bad param co&
&unt (eta).')
endif
endif
if(xyzzyaady1(2)>0)then
xyzzyaaaf24=xyzzyaaaj24
xyzzyaaai24=xyzzyaaaj24+1
xyzzyaaaj24=xyzzyaaaj24+xyzzyaady1(2)
if(any(.not.ignore(xyzzyaaai24:xyzzyaaaj24)))then
do xyzzyaaag24=1,xyzzyaaah1
if(xyzzyaabq1(xyzzyaaag24)==1)then
xyzzyaaaf24=xyzzyaaaf24+1
if(.not.ignore(xyzzyaaaf24))xyzzyaabj1(xyzzyaaag24)=params(xyzzyaaaf24&
&)
endif
do xyzzyaaad24=1,xyzzyaaas1(xyzzyaaag24)
do xyzzyaaae24=1,xyzzyaaap1(xyzzyaaag24)
if(xyzzyaabr1(xyzzyaaae24,xyzzyaaad24,xyzzyaaag24)==1)then
xyzzyaaaf24=xyzzyaaaf24+1
if(.not.ignore(xyzzyaaaf24))then
params(xyzzyaaaf24)=params(xyzzyaaaf24)/xyzzyaabj1(xyzzyaaag24)**(xyzz&
&yaaae24-1)
xyzzyaabk1(xyzzyaaae24,xyzzyaaad24,xyzzyaaag24)=params(xyzzyaaaf24)
endif
endif
enddo
enddo
enddo
call xyzzyaafj1
do xyzzyaaag24=1,xyzzyaaah1
if(xyzzyaabq1(xyzzyaaag24)==1)then
xyzzyaach1(xyzzyaaag24)=1.d0/xyzzyaabj1(xyzzyaaag24)
select case(xyzzyaaad1)
case(1)
xyzzyaacf1(xyzzyaaag24)=-xyzzyaach1(xyzzyaaag24)
xyzzyaacg1(xyzzyaaag24)=0.d0
case(2)
xyzzyaacf1(xyzzyaaag24)=-2*xyzzyaach1(xyzzyaaag24)
xyzzyaacg1(xyzzyaaag24)=2*xyzzyaach1(xyzzyaaag24)**2
case(3)
xyzzyaacf1(xyzzyaaag24)=-3*xyzzyaach1(xyzzyaaag24)
xyzzyaacg1(xyzzyaaag24)=6*xyzzyaach1(xyzzyaaag24)**2
case default
xyzzyaacf1(xyzzyaaag24)=-xyzzyaaad1*xyzzyaach1(xyzzyaaag24)
xyzzyaacg1(xyzzyaaag24)=xyzzyaaad1*xyzzyaaae1*xyzzyaach1(xyzzyaaag24)*&
&*2
end select
endif
enddo
if(xyzzyaaaj24/=xyzzyaaaf24)call errstop('PUT_BF_PARAMS','Bad param co&
&unt (mu).')
endif
endif
if(xyzzyaady1(3)>0)then
xyzzyaaaf24=xyzzyaaaj24
xyzzyaaai24=xyzzyaaaj24+1
xyzzyaaaj24=xyzzyaaaj24+xyzzyaady1(3)
if(any(.not.ignore(xyzzyaaai24:xyzzyaaaj24)))then
do xyzzyaaag24=1,xyzzyaaai1
if(xyzzyaabs1(xyzzyaaag24)==1)then
xyzzyaaaf24=xyzzyaaaf24+1
if(.not.ignore(xyzzyaaaf24))then
xyzzyaabl1(xyzzyaaag24)=params(xyzzyaaaf24)
do xyzzyaaad24=1,xyzzyaaaz1(xyzzyaaag24)
allocate(xyzzyaaal24(xyzzyaacw1(xyzzyaaad24,xyzzyaaag24),xyzzyaacx1(xy&
&zzyaaad24,xyzzyaaag24)),stat=xyzzyaaah24)
call check_alloc(xyzzyaaah24,'PUT_BF_PARAMS','')
call xyzzyaafm1(xyzzyaaag24,xyzzyaacs1(xyzzyaaad24,xyzzyaaag24),xyzzya&
&acw1(xyzzyaaad24,xyzzyaaag24),xyzzyaacx1(xyzzyaaad24,xyzzyaaag24),xyz&
&zyaaal24)
xyzzyaada1(:,:,xyzzyaaad24,xyzzyaaag24)=0.d0
xyzzyaada1(1:xyzzyaacw1(xyzzyaaad24,xyzzyaaag24),1:xyzzyaacx1(xyzzyaaa&
&d24,xyzzyaaag24),xyzzyaaad24,xyzzyaaag24)=xyzzyaaal24
deallocate(xyzzyaaal24)
enddo
xyzzyaace1(xyzzyaaag24)=1.d0/xyzzyaabl1(xyzzyaaag24)
select case(xyzzyaaad1)
case(1)
xyzzyaacc1(xyzzyaaag24)=-xyzzyaace1(xyzzyaaag24)
xyzzyaacd1(xyzzyaaag24)=0.d0
case(2)
xyzzyaacc1(xyzzyaaag24)=-2*xyzzyaace1(xyzzyaaag24)
xyzzyaacd1(xyzzyaaag24)=2*xyzzyaace1(xyzzyaaag24)**2
case(3)
xyzzyaacc1(xyzzyaaag24)=-3*xyzzyaace1(xyzzyaaag24)
xyzzyaacd1(xyzzyaaag24)=6*xyzzyaace1(xyzzyaaag24)**2
case default
xyzzyaacc1(xyzzyaaag24)=-xyzzyaaad1*xyzzyaace1(xyzzyaaag24)
xyzzyaacd1(xyzzyaaag24)=xyzzyaaad1*xyzzyaaae1*xyzzyaace1(xyzzyaaag24)*&
&*2
end select
endif
endif
do xyzzyaaad24=1,xyzzyaaaz1(xyzzyaaag24)
do xyzzyaaac24=1,xyzzyaaav1(xyzzyaaag24)
do xyzzyaaab24=1,xyzzyaaaw1(xyzzyaaag24)
do xyzzyaaaa24=1,xyzzyaaaw1(xyzzyaaag24)
if(xyzzyaabt1(xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaa&
&g24)==1)then
xyzzyaaaf24=xyzzyaaaf24+1
if(.not.ignore(xyzzyaaaf24))then
params(xyzzyaaaf24)=params(xyzzyaaaf24)/xyzzyaabl1(xyzzyaaag24)**(xyzz&
&yaaaa24+xyzzyaaab24+xyzzyaaac24-3)
xyzzyaabm1(xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaag24&
&)=params(xyzzyaaaf24)
endif
endif
enddo
enddo
enddo
do xyzzyaaac24=1,xyzzyaaav1(xyzzyaaag24)
do xyzzyaaab24=1,xyzzyaaaw1(xyzzyaaag24)
do xyzzyaaaa24=1,xyzzyaaaw1(xyzzyaaag24)
if(xyzzyaabu1(xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaa&
&g24)==1)then
xyzzyaaaf24=xyzzyaaaf24+1
if(.not.ignore(xyzzyaaaf24))then
params(xyzzyaaaf24)=params(xyzzyaaaf24)/xyzzyaabl1(xyzzyaaag24)**(xyzz&
&yaaaa24+xyzzyaaab24+xyzzyaaac24-3)
xyzzyaabn1(xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaag24&
&)=params(xyzzyaaaf24)
endif
endif
enddo
enddo
enddo
enddo
enddo
call xyzzyaafn1
xyzzyaacp1=maxval(xyzzyaabl1(:))
if(xyzzyaaaj24/=xyzzyaaaf24)call errstop('PUT_BF_PARAMS','Bad param co&
&unt (Phi+Theta).')
endif
endif
if(xyzzyaaaa1.or.xyzzyaaac1)then
xyzzyaacq1=0.d0
if(xyzzyaaaa1)xyzzyaacq1=max(xyzzyaacq1,maxval(xyzzyaabh1))
if(xyzzyaaac1)xyzzyaacq1=max(xyzzyaacq1,2*xyzzyaacp1)
endif
if(xyzzyaaab1.or.xyzzyaaac1)then
do xyzzyaaak24=1,nitot
xyzzyaaci1(xyzzyaaak24)=0.d0
if(xyzzyaaab1)then
if(xyzzyaaan1(xyzzyaaak24)/=0)xyzzyaaci1(xyzzyaaak24)=max(xyzzyaaci1(x&
&yzzyaaak24),xyzzyaabj1(xyzzyaaan1(xyzzyaaak24)))
endif
if(xyzzyaaac1)then
if(xyzzyaaat1(xyzzyaaak24)/=0)xyzzyaaci1(xyzzyaaak24)=max(xyzzyaaci1(x&
&yzzyaaak24),xyzzyaabl1(xyzzyaaat1(xyzzyaaak24)))
endif
enddo
endif
if(xyzzyaady1(4)>0)then
xyzzyaaaf24=xyzzyaaaj24
xyzzyaaai24=xyzzyaaaj24+1
xyzzyaaaj24=xyzzyaaaj24+xyzzyaady1(4)
if(any(.not.ignore(xyzzyaaai24:xyzzyaaaj24)))then
do xyzzyaaag24=1,xyzzyaadc1
if(xyzzyaadd1(xyzzyaaag24)==1)then
xyzzyaaaf24=xyzzyaaaf24+1
if(.not.ignore(xyzzyaaaf24))xyzzyaadj1(xyzzyaaag24)=params(xyzzyaaaf24&
&)
endif
enddo
do xyzzyaaak24=1,nitot
if(.not.(is_ae(xyzzyaaak24)))cycle
xyzzyaaag24=xyzzyaade1(xyzzyaaak24)
xyzzyaadf1(xyzzyaaak24)=xyzzyaadj1(xyzzyaaag24)
xyzzyaadg1(xyzzyaaak24)=1.d0/xyzzyaadf1(xyzzyaaak24)
xyzzyaadh1(xyzzyaaak24)=12*xyzzyaadg1(xyzzyaaak24)
xyzzyaadi1(xyzzyaaak24)=xyzzyaadh1(xyzzyaaak24)*xyzzyaadg1(xyzzyaaak24&
&)
enddo
if(xyzzyaaaj24/=xyzzyaaaf24)call errstop('PUT_BF_PARAMS','Bad param co&
&unt (AE cutoff).')
endif
endif
call xyzzyaafg1(iparam_buffer)
end subroutine put_pbf_params
subroutine xyzzyaafg1(indx)
implicit none
integer,intent(in) :: indx
integer xyzzyaaaa25,xyzzyaaab25
xyzzyaaab25=0
xyzzyaaaa25=0
if(indx/=0)then
xyzzyaaab25=xyzzyaadw1(indx)
xyzzyaaaa25=indx-sum(xyzzyaady1(1:xyzzyaaab25-1))
endif
if(xyzzyaady1(1)>0.and.(xyzzyaaab25==1.or.xyzzyaaab25==0))then
call dcopy(xyzzyaaak1,xyzzyaabh1(1),1,xyzzyaadz1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaee1,xyzzyaabi1(1,1),1,xyzzyaaea1(1,1,xyzzyaaaa25),1)
call dcopy(xyzzyaaak1,xyzzyaacb1(1),1,xyzzyaaeb1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaak1,xyzzyaabz1(1),1,xyzzyaaec1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaak1,xyzzyaaca1(1),1,xyzzyaaed1(1,xyzzyaaaa25),1)
endif
if(xyzzyaady1(2)>0.and.(xyzzyaaab25==2.or.xyzzyaaab25==0))then
call dcopy(xyzzyaaah1,xyzzyaabj1(1),1,xyzzyaaef1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaek1,xyzzyaabk1(1,1,1),1,xyzzyaaeg1(1,1,1,xyzzyaaaa25&
&),1)
call dcopy(xyzzyaaah1,xyzzyaach1(1),1,xyzzyaaeh1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaah1,xyzzyaacf1(1),1,xyzzyaaei1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaah1,xyzzyaacg1(1),1,xyzzyaaej1(1,xyzzyaaaa25),1)
endif
if(xyzzyaady1(3)>0.and.(xyzzyaaab25==3.or.xyzzyaaab25==0))then
call dcopy(xyzzyaaai1,xyzzyaabl1(1),1,xyzzyaael1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaai1,xyzzyaace1(1),1,xyzzyaaem1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaai1,xyzzyaacc1(1),1,xyzzyaaen1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaai1,xyzzyaacd1(1),1,xyzzyaaeo1(1,xyzzyaaaa25),1)
call dcopy(xyzzyaaes1,xyzzyaabm1(1,1,1,1,1),1,xyzzyaaep1(1,1,1,1,1,xyz&
&zyaaaa25),1)
call dcopy(xyzzyaaes1,xyzzyaabn1(1,1,1,1,1),1,xyzzyaaeq1(1,1,1,1,1,xyz&
&zyaaaa25),1)
xyzzyaaer1(xyzzyaaaa25)=xyzzyaacp1
endif
if(xyzzyaady1(4)>0.and.(xyzzyaaab25==4.or.xyzzyaaab25==0))then
call dcopy(xyzzyaadc1,xyzzyaadj1(1),1,xyzzyaaet1(1,xyzzyaaaa25),1)
call dcopy(nitot,xyzzyaadf1(1),1,xyzzyaaeu1(1,xyzzyaaaa25),1)
call dcopy(nitot,xyzzyaadg1(1),1,xyzzyaaev1(1,xyzzyaaaa25),1)
call dcopy(nitot,xyzzyaadh1(1),1,xyzzyaaew1(1,xyzzyaaaa25),1)
call dcopy(nitot,xyzzyaadi1(1),1,xyzzyaaex1(1,xyzzyaaaa25),1)
endif
end subroutine xyzzyaafg1
subroutine xyzzyaafh1(indx)
implicit none
integer,intent(in) :: indx
integer xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26
xyzzyaaaa26=0
xyzzyaaab26=0
if(indx/=0)then
xyzzyaaaa26=xyzzyaadw1(indx)
xyzzyaaab26=indx-sum(xyzzyaady1(1:xyzzyaaaa26-1))
endif
if(xyzzyaady1(1)>0.and.(xyzzyaaaa26==1.or.xyzzyaaaa26==0))then
call dcopy(xyzzyaaak1,xyzzyaadz1(1,xyzzyaaab26),1,xyzzyaabh1(1),1)
call dcopy(xyzzyaaee1,xyzzyaaea1(1,1,xyzzyaaab26),1,xyzzyaabi1(1,1),1)
call dcopy(xyzzyaaak1,xyzzyaaeb1(1,xyzzyaaab26),1,xyzzyaacb1(1),1)
call dcopy(xyzzyaaak1,xyzzyaaec1(1,xyzzyaaab26),1,xyzzyaabz1(1),1)
call dcopy(xyzzyaaak1,xyzzyaaed1(1,xyzzyaaab26),1,xyzzyaaca1(1),1)
endif
if(xyzzyaady1(2)>0.and.(xyzzyaaaa26==2.or.xyzzyaaaa26==0))then
call dcopy(xyzzyaaah1,xyzzyaaef1(1,xyzzyaaab26),1,xyzzyaabj1(1),1)
call dcopy(xyzzyaaek1,xyzzyaaeg1(1,1,1,xyzzyaaab26),1,xyzzyaabk1(1,1,1&
&),1)
call dcopy(xyzzyaaah1,xyzzyaaeh1(1,xyzzyaaab26),1,xyzzyaach1(1),1)
call dcopy(xyzzyaaah1,xyzzyaaei1(1,xyzzyaaab26),1,xyzzyaacf1(1),1)
call dcopy(xyzzyaaah1,xyzzyaaej1(1,xyzzyaaab26),1,xyzzyaacg1(1),1)
endif
if(xyzzyaady1(3)>0.and.(xyzzyaaaa26==3.or.xyzzyaaaa26==0))then
call dcopy(xyzzyaaai1,xyzzyaael1(1,xyzzyaaab26),1,xyzzyaabl1(1),1)
call dcopy(xyzzyaaai1,xyzzyaaem1(1,xyzzyaaab26),1,xyzzyaace1(1),1)
call dcopy(xyzzyaaai1,xyzzyaaen1(1,xyzzyaaab26),1,xyzzyaacc1(1),1)
call dcopy(xyzzyaaai1,xyzzyaaeo1(1,xyzzyaaab26),1,xyzzyaacd1(1),1)
call dcopy(xyzzyaaes1,xyzzyaaep1(1,1,1,1,1,xyzzyaaab26),1,xyzzyaabm1(1&
&,1,1,1,1),1)
call dcopy(xyzzyaaes1,xyzzyaaeq1(1,1,1,1,1,xyzzyaaab26),1,xyzzyaabn1(1&
&,1,1,1,1),1)
xyzzyaacp1=xyzzyaaer1(xyzzyaaab26)
endif
if(xyzzyaady1(4)>0.and.(xyzzyaaaa26==4.or.xyzzyaaaa26==0))then
call dcopy(xyzzyaadc1,xyzzyaaet1(1,xyzzyaaab26),1,xyzzyaadj1(1),1)
call dcopy(nitot,xyzzyaaeu1(1,xyzzyaaab26),1,xyzzyaadf1(1),1)
call dcopy(nitot,xyzzyaaev1(1,xyzzyaaab26),1,xyzzyaadg1(1),1)
call dcopy(nitot,xyzzyaaew1(1,xyzzyaaab26),1,xyzzyaadh1(1),1)
call dcopy(nitot,xyzzyaaex1(1,xyzzyaaab26),1,xyzzyaadi1(1),1)
endif
if(xyzzyaaaa1.or.xyzzyaaac1)then
xyzzyaacq1=0.d0
if(xyzzyaaaa1)xyzzyaacq1=max(xyzzyaacq1,maxval(xyzzyaabh1))
if(xyzzyaaac1)xyzzyaacq1=max(xyzzyaacq1,2*xyzzyaacp1)
endif
if(xyzzyaaab1.or.xyzzyaaac1)then
do xyzzyaaac26=1,nitot
xyzzyaaci1(xyzzyaaac26)=0.d0
if(xyzzyaaab1)then
if(xyzzyaaan1(xyzzyaaac26)/=0)xyzzyaaci1(xyzzyaaac26)=max(xyzzyaaci1(x&
&yzzyaaac26),xyzzyaabj1(xyzzyaaan1(xyzzyaaac26)))
endif
if(xyzzyaaac1)then
if(xyzzyaaat1(xyzzyaaac26)/=0)xyzzyaaci1(xyzzyaaac26)=max(xyzzyaaci1(x&
&yzzyaaac26),xyzzyaabl1(xyzzyaaat1(xyzzyaaac26)))
endif
enddo
endif
end subroutine xyzzyaafh1
subroutine xyzzyaafi1
implicit none
integer xyzzyaaaa27
real(dp) xyzzyaaab27
do xyzzyaaaa27=1,xyzzyaaak1
if(xyzzyaacr1(xyzzyaaaa27))then
xyzzyaabp1(2,xyzzyaaaa27)=-1
xyzzyaaab27=real(xyzzyaaad1,dp)/xyzzyaabh1(xyzzyaaaa27)
xyzzyaabi1(2,xyzzyaaaa27)=xyzzyaaab27*xyzzyaabi1(1,xyzzyaaaa27)
endif
enddo
end subroutine xyzzyaafi1
subroutine xyzzyaafj1
implicit none
integer xyzzyaaaa28,xyzzyaaab28
real(dp) xyzzyaaac28,xyzzyaaad28
do xyzzyaaab28=1,xyzzyaaah1
xyzzyaaad28=real(xyzzyaaad1,dp)/xyzzyaabj1(xyzzyaaab28)
do xyzzyaaaa28=1,xyzzyaaas1(xyzzyaaab28)
if(xyzzyaact1(xyzzyaaab28))then
xyzzyaabk1(1,xyzzyaaaa28,xyzzyaaab28)=0.d0
xyzzyaabr1(1,xyzzyaaaa28,xyzzyaaab28)=-1
endif
xyzzyaaac28=xyzzyaabk1(1,xyzzyaaaa28,xyzzyaaab28)
xyzzyaabk1(2,xyzzyaaaa28,xyzzyaaab28)=xyzzyaaac28*xyzzyaaad28
xyzzyaabr1(2,xyzzyaaaa28,xyzzyaaab28)=-1
enddo
enddo
end subroutine xyzzyaafj1
subroutine xyzzyaafk1
implicit none
integer xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaae29,xy&
&zzyaaaf29,xyzzyaaag29,xyzzyaaah29,xyzzyaaai29,xyzzyaaaj29,xyzzyaaak29
real(dp),allocatable :: xyzzyaaal29(:,:)
logical,allocatable  :: xyzzyaaam29(:)
xyzzyaaak29=maxval(xyzzyaaaz1)
allocate(xyzzyaacw1(xyzzyaaak29,xyzzyaaai1),xyzzyaacx1(xyzzyaaak29,xyz&
&zyaaai1),stat=xyzzyaaaa29)
call check_alloc(xyzzyaaaa29,'CONSTRAINT_MATRIX','1')
xyzzyaacx1=0
xyzzyaacw1=0
do xyzzyaaab29=1,xyzzyaaai1
do xyzzyaaah29=1,xyzzyaaaz1(xyzzyaaab29)
xyzzyaacx1(xyzzyaaah29,xyzzyaaab29)=2*xyzzyaaav1(xyzzyaaab29)*xyzzyaaa&
&w1(xyzzyaaab29)*xyzzyaaaw1(xyzzyaaab29)
if(xyzzyaacs1(xyzzyaaah29,xyzzyaaab29))xyzzyaacw1(xyzzyaaah29,xyzzyaaa&
&b29)=2*xyzzyaaaw1(xyzzyaaab29)-1
xyzzyaacw1(xyzzyaaah29,xyzzyaaab29)=xyzzyaacw1(xyzzyaaah29,xyzzyaaab29&
&)+2*xyzzyaaaw1(xyzzyaaab29)-1
if(xyzzyaacu1(xyzzyaaab29))then
xyzzyaacw1(xyzzyaaah29,xyzzyaaab29)=xyzzyaacw1(xyzzyaaah29,xyzzyaaab29&
&)+11*(xyzzyaaav1(xyzzyaaab29)+xyzzyaaaw1(xyzzyaaab29))-15
else
xyzzyaacw1(xyzzyaaah29,xyzzyaaab29)=xyzzyaacw1(xyzzyaaah29,xyzzyaaab29&
&)+3*(xyzzyaaav1(xyzzyaaab29)+xyzzyaaaw1(xyzzyaaab29))-3
endif
if(xyzzyaacv1(xyzzyaaab29))then
if(xyzzyaaad1>0)then
xyzzyaacw1(xyzzyaaah29,xyzzyaaab29)=xyzzyaacw1(xyzzyaaah29,xyzzyaaab29&
&)+xyzzyaaav1(xyzzyaaab29)*xyzzyaaaw1(xyzzyaaab29)*xyzzyaaaw1(xyzzyaaa&
&b29)+2*xyzzyaaaw1(xyzzyaaab29)*(xyzzyaaav1(xyzzyaaab29)-1)+xyzzyaaaw1&
&(xyzzyaaab29)*xyzzyaaaw1(xyzzyaaab29)
else
xyzzyaacw1(xyzzyaaah29,xyzzyaaab29)=xyzzyaacw1(xyzzyaaah29,xyzzyaaab29&
&)+xyzzyaaav1(xyzzyaaab29)*xyzzyaaaw1(xyzzyaaab29)*xyzzyaaaw1(xyzzyaaa&
&b29)+xyzzyaaaw1(xyzzyaaab29)*(xyzzyaaaw1(xyzzyaaab29)+xyzzyaaav1(xyzz&
&yaaab29)-2)
endif
endif
enddo
enddo
xyzzyaaai29=maxval(xyzzyaacw1)
xyzzyaaaj29=maxval(xyzzyaacx1)
xyzzyaaak29=maxval(xyzzyaaaz1)
allocate(xyzzyaacz1(xyzzyaaaj29),xyzzyaada1(xyzzyaaai29,xyzzyaaaj29,xy&
&zzyaaak29,xyzzyaaai1),xyzzyaacy1(xyzzyaaai29,xyzzyaaak29,xyzzyaaai1),&
&stat=xyzzyaaaa29)
call check_alloc(xyzzyaaaa29,'CONSTRAINT_MATRIX','2')
do xyzzyaaab29=1,xyzzyaaai1
do xyzzyaaah29=1,xyzzyaaaz1(xyzzyaaab29)
allocate(xyzzyaaal29(xyzzyaacw1(xyzzyaaah29,xyzzyaaab29),xyzzyaacx1(xy&
&zzyaaah29,xyzzyaaab29)),stat=xyzzyaaaa29)
call check_alloc(xyzzyaaaa29,'CONSTRAINT_MATRIX','3')
call xyzzyaafm1(xyzzyaaab29,xyzzyaacs1(xyzzyaaah29,xyzzyaaab29),xyzzya&
&acw1(xyzzyaaah29,xyzzyaaab29),xyzzyaacx1(xyzzyaaah29,xyzzyaaab29),xyz&
&zyaaal29)
xyzzyaada1(:,:,xyzzyaaah29,xyzzyaaab29)=0.d0
xyzzyaada1(1:xyzzyaacw1(xyzzyaaah29,xyzzyaaab29),1:xyzzyaacx1(xyzzyaaa&
&h29,xyzzyaaab29),xyzzyaaah29,xyzzyaaab29)=xyzzyaaal29
deallocate(xyzzyaaal29)
allocate(xyzzyaaam29(xyzzyaacx1(xyzzyaaah29,xyzzyaaab29)),stat=xyzzyaa&
&aa29)
call check_alloc(xyzzyaaaa29,'CONSTRAINT_MATRIX','4')
xyzzyaaam29=.false.
do xyzzyaaac29=1,xyzzyaacw1(xyzzyaaah29,xyzzyaaab29)
xyzzyaacy1(xyzzyaaac29,xyzzyaaah29,xyzzyaaab29)=-1
do xyzzyaaad29=xyzzyaaac29,xyzzyaacx1(xyzzyaaah29,xyzzyaaab29)
if(xyzzyaada1(xyzzyaaac29,xyzzyaaad29,xyzzyaaah29,xyzzyaaab29)>0.5d0)t&
&hen
xyzzyaacy1(xyzzyaaac29,xyzzyaaah29,xyzzyaaab29)=xyzzyaaad29
exit
endif
enddo
if(xyzzyaacy1(xyzzyaaac29,xyzzyaaah29,xyzzyaaab29)>0)xyzzyaaam29(xyzzy&
&aacy1(xyzzyaaac29,xyzzyaaah29,xyzzyaaab29))=.true.
enddo
xyzzyaaac29=0
do xyzzyaaag29=1,xyzzyaaav1(xyzzyaaab29)
do xyzzyaaaf29=1,xyzzyaaaw1(xyzzyaaab29)
do xyzzyaaae29=1,xyzzyaaaw1(xyzzyaaab29)
xyzzyaaac29=xyzzyaaac29+1
if(xyzzyaaam29(xyzzyaaac29))xyzzyaabt1(xyzzyaaae29,xyzzyaaaf29,xyzzyaa&
&ag29,xyzzyaaah29,xyzzyaaab29)=-1
enddo
enddo
enddo
do xyzzyaaag29=1,xyzzyaaav1(xyzzyaaab29)
do xyzzyaaaf29=1,xyzzyaaaw1(xyzzyaaab29)
do xyzzyaaae29=1,xyzzyaaaw1(xyzzyaaab29)
xyzzyaaac29=xyzzyaaac29+1
if(xyzzyaaam29(xyzzyaaac29))xyzzyaabu1(xyzzyaaae29,xyzzyaaaf29,xyzzyaa&
&ag29,xyzzyaaah29,xyzzyaaab29)=-1
enddo
enddo
enddo
deallocate(xyzzyaaam29)
enddo
enddo
end subroutine xyzzyaafk1
subroutine xyzzyaafl1(set,no_vars,cee,x_determined,no_free)
implicit none
integer,intent(in) :: set,no_vars
integer,intent(out) :: no_free
logical,intent(in) :: cee
logical,intent(out) :: x_determined(:)
integer xyzzyaaaa30,xyzzyaaab30,xyzzyaaac30,xyzzyaaad30,xyzzyaaae30
real(dp),allocatable :: xyzzyaaaf30(:,:)
xyzzyaaae30=0
if(cee)xyzzyaaae30=2*xyzzyaaaw1(set)-1
xyzzyaaae30=xyzzyaaae30+2*xyzzyaaaw1(set)-1
if(xyzzyaacu1(set))then
xyzzyaaae30=xyzzyaaae30+11*(xyzzyaaav1(set)+xyzzyaaaw1(set))-15
else
xyzzyaaae30=xyzzyaaae30+3*(xyzzyaaav1(set)+xyzzyaaaw1(set))-3
endif
if(xyzzyaacv1(set))then
if(xyzzyaaad1>0)then
xyzzyaaae30=xyzzyaaae30+xyzzyaaav1(set)*xyzzyaaaw1(set)*xyzzyaaaw1(set&
&)+2*xyzzyaaaw1(set)*(xyzzyaaav1(set)-1)+xyzzyaaaw1(set)*xyzzyaaaw1(se&
&t)
else
xyzzyaaae30=xyzzyaaae30+xyzzyaaav1(set)*xyzzyaaaw1(set)*xyzzyaaaw1(set&
&)+xyzzyaaaw1(set)*(xyzzyaaaw1(set)+xyzzyaaav1(set)-2)
endif
endif
allocate(xyzzyaaaf30(xyzzyaaae30,no_vars),stat=xyzzyaaad30)
call check_alloc(xyzzyaaad30,'FIND_DETERMINED_PHI','')
call xyzzyaafm1(set,cee,xyzzyaaae30,no_vars,xyzzyaaaf30)
no_free=no_vars
x_determined=.false.
do xyzzyaaaa30=1,xyzzyaaae30
xyzzyaaac30=-1
do xyzzyaaab30=xyzzyaaaa30,no_vars
if(xyzzyaaaf30(xyzzyaaaa30,xyzzyaaab30)>0.5d0)then
xyzzyaaac30=xyzzyaaab30
exit
endif
enddo
if(xyzzyaaac30>0)then
x_determined(xyzzyaaac30)=.true.
no_free=no_free-1
endif
enddo
deallocate(xyzzyaaaf30)
end subroutine xyzzyaafl1
subroutine xyzzyaafm1(set,cee,no_eqns,no_vars,c)
implicit none
integer,intent(in) :: set,no_eqns,no_vars
real(dp),intent(out) :: c(no_eqns,no_vars)
logical,intent(in) :: cee
integer xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31,xyzzyaaad31,xyzzyaaae31,xy&
&zzyaaaf31,xyzzyaaag31,xyzzyaaah31,xyzzyaaai31,xyzzyaaaj31,xyzzyaaak31&
&,xyzzyaaal31,xyzzyaaam31,xyzzyaaan31,xyzzyaaao31,xyzzyaaap31,xyzzyaaa&
&q31,xyzzyaaar31,xyzzyaaas31,xyzzyaaat31,xyzzyaaau31
real(dp) xyzzyaaav31,xyzzyaaaw31
c=0.d0
xyzzyaaaw31=xyzzyaabl1(set)
xyzzyaaav31=real(xyzzyaaad1,dp)/xyzzyaaaw31
xyzzyaaag31=2*xyzzyaaaw1(set)-1
xyzzyaaah31=xyzzyaaav1(set)+xyzzyaaaw1(set)-1
xyzzyaaai31=xyzzyaaah31-1
xyzzyaaaa31=0
xyzzyaaap31=0
xyzzyaaaj31=xyzzyaaap31
if(cee)xyzzyaaaj31=xyzzyaaap31+xyzzyaaag31
xyzzyaaak31=xyzzyaaaj31+xyzzyaaah31
xyzzyaaal31=xyzzyaaak31+xyzzyaaah31
xyzzyaaam31=xyzzyaaal31+xyzzyaaah31
xyzzyaaan31=xyzzyaaam31+xyzzyaaah31
xyzzyaaao31=xyzzyaaan31+xyzzyaaai31
do xyzzyaaad31=1,xyzzyaaav1(set)
do xyzzyaaac31=1,xyzzyaaaw1(set)
do xyzzyaaab31=1,xyzzyaaaw1(set)
xyzzyaaaa31=xyzzyaaaa31+1
if(cee)then
if(xyzzyaaad31==2)then
xyzzyaaae31=xyzzyaaab31+xyzzyaaac31-1
c(xyzzyaaae31+xyzzyaaap31,xyzzyaaaa31)=1.d0
endif
endif
if(xyzzyaacu1(set))then
if(xyzzyaaac31==1)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaak31,xyzzyaaaa31)=1.d0
if(xyzzyaaad31>1)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-2
c(xyzzyaaaf31+xyzzyaaao31,xyzzyaaaa31)=real(xyzzyaaad31-1,dp)
endif
elseif(xyzzyaaac31==2)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaam31,xyzzyaaaa31)=1.d0
endif
if(xyzzyaaab31==1)then
xyzzyaaaf31=xyzzyaaac31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaaj31,xyzzyaaaa31)=1.d0
if(xyzzyaaad31>1)then
xyzzyaaaf31=xyzzyaaac31+xyzzyaaad31-2
c(xyzzyaaaf31+xyzzyaaan31,xyzzyaaaa31)=real(xyzzyaaad31-1,dp)
endif
elseif(xyzzyaaab31==2)then
xyzzyaaaf31=xyzzyaaac31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaal31,xyzzyaaaa31)=1.d0
endif
else
if(xyzzyaaac31==1)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaak31,xyzzyaaaa31)=-xyzzyaaav31
elseif(xyzzyaaac31==2)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaak31,xyzzyaaaa31)=1.d0
endif
if(xyzzyaaab31==1)then
xyzzyaaaf31=xyzzyaaac31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaaj31,xyzzyaaaa31)=-xyzzyaaav31
elseif(xyzzyaaab31==2)then
xyzzyaaaf31=xyzzyaaac31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaaj31,xyzzyaaaa31)=1.d0
endif
endif
enddo
enddo
enddo
if(xyzzyaacu1(set))then
xyzzyaaap31=xyzzyaaao31+xyzzyaaai31
else
xyzzyaaap31=xyzzyaaak31+xyzzyaaah31
endif
xyzzyaaaj31=xyzzyaaap31+xyzzyaaag31
xyzzyaaal31=xyzzyaaaj31+xyzzyaaah31
xyzzyaaam31=xyzzyaaal31+xyzzyaaah31
xyzzyaaan31=xyzzyaaam31+xyzzyaaah31
xyzzyaaao31=xyzzyaaan31+xyzzyaaai31
do xyzzyaaad31=1,xyzzyaaav1(set)
do xyzzyaaac31=1,xyzzyaaaw1(set)
do xyzzyaaab31=1,xyzzyaaaw1(set)
xyzzyaaaa31=xyzzyaaaa31+1
if(xyzzyaaad31==2)then
xyzzyaaae31=xyzzyaaab31+xyzzyaaac31-1
c(xyzzyaaae31+xyzzyaaap31,xyzzyaaaa31)=1.d0
endif
if(xyzzyaacu1(set))then
if(xyzzyaaac31==1)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaam31,xyzzyaaaa31)=-xyzzyaaav31
if(xyzzyaaad31>1)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-2
c(xyzzyaaaf31+xyzzyaaao31,xyzzyaaaa31)=real(xyzzyaaad31-1,dp)
endif
elseif(xyzzyaaac31==2)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaam31,xyzzyaaaa31)=1.d0
endif
if(xyzzyaaab31==1)then
xyzzyaaaf31=xyzzyaaac31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaaj31,xyzzyaaaa31)=1.d0
if(xyzzyaaad31>1)then
xyzzyaaaf31=xyzzyaaac31+xyzzyaaad31-2
c(xyzzyaaaf31+xyzzyaaan31,xyzzyaaaa31)=real(xyzzyaaad31-1,dp)
endif
elseif(xyzzyaaab31==2)then
xyzzyaaaf31=xyzzyaaac31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaal31,xyzzyaaaa31)=1.d0
endif
else
if(xyzzyaaac31==1)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaaj31,xyzzyaaaa31)=-xyzzyaaav31
elseif(xyzzyaaac31==2)then
xyzzyaaaf31=xyzzyaaab31+xyzzyaaad31-1
c(xyzzyaaaf31+xyzzyaaaj31,xyzzyaaaa31)=1.d0
endif
endif
enddo
enddo
enddo
if(xyzzyaacu1(set))then
xyzzyaaaq31=xyzzyaaao31+xyzzyaaai31
else
xyzzyaaaq31=xyzzyaaaj31+xyzzyaaah31
endif
if(xyzzyaacv1(set))then
xyzzyaaaa31=0
xyzzyaaae31=xyzzyaaaq31
xyzzyaaas31=1
xyzzyaaat31=xyzzyaaas31*xyzzyaaaw1(set)
xyzzyaaau31=xyzzyaaat31*xyzzyaaaw1(set)
xyzzyaaar31=xyzzyaaau31*xyzzyaaav1(set)
do xyzzyaaad31=1,xyzzyaaav1(set)
do xyzzyaaac31=1,xyzzyaaaw1(set)
do xyzzyaaab31=1,xyzzyaaaw1(set)
xyzzyaaaa31=xyzzyaaaa31+1
xyzzyaaae31=xyzzyaaae31+1
if(xyzzyaaad1>0)then
if(xyzzyaaad31>1)then
c(xyzzyaaae31,xyzzyaaaa31-1*xyzzyaaau31)=real(xyzzyaaad1+xyzzyaaab31-1&
&,dp)
if(xyzzyaaab31<xyzzyaaaw1(set))c(xyzzyaaae31,xyzzyaaaa31+1*xyzzyaaas31&
&-1*xyzzyaaau31)=-xyzzyaaaw31*real(xyzzyaaab31,dp)
endif
if(xyzzyaaad31<xyzzyaaav1(set))then
if(xyzzyaaab31>2)c(xyzzyaaae31,xyzzyaaaa31+xyzzyaaar31-2*xyzzyaaas31+1&
&*xyzzyaaau31)=-real(xyzzyaaad31,dp)
if(xyzzyaaab31>1)c(xyzzyaaae31,xyzzyaaaa31+xyzzyaaar31-1*xyzzyaaas31+1&
&*xyzzyaaau31)=xyzzyaaaw31*real(xyzzyaaad31,dp)
endif
else
if(xyzzyaaad31>1.and.xyzzyaaab31<xyzzyaaaw1(set))c(xyzzyaaae31,xyzzyaa&
&aa31+1*xyzzyaaas31-1*xyzzyaaau31)=real(xyzzyaaab31,dp)
if(xyzzyaaab31>1.and.xyzzyaaad31<xyzzyaaaw1(set))c(xyzzyaaae31,xyzzyaa&
&aa31+xyzzyaaar31-1*xyzzyaaas31+1*xyzzyaaau31)=-real(xyzzyaaad31,dp)
endif
enddo
enddo
enddo
if(xyzzyaaad1>0)then
xyzzyaaaa31=(xyzzyaaav1(set)-1)*xyzzyaaaw1(set)**2
do xyzzyaaac31=1,xyzzyaaaw1(set)
do xyzzyaaab31=1,xyzzyaaaw1(set)
xyzzyaaaa31=xyzzyaaaa31+1
xyzzyaaae31=xyzzyaaae31+1
c(xyzzyaaae31,xyzzyaaaa31)=real(xyzzyaaad1+xyzzyaaab31-1,dp)
if(xyzzyaaab31<xyzzyaaaw1(set))c(xyzzyaaae31,xyzzyaaaa31+1*xyzzyaaas31&
&)=-xyzzyaaaw31*real(xyzzyaaab31,dp)
enddo
enddo
xyzzyaaaa31=xyzzyaaaw1(set)-1-xyzzyaaat31
do xyzzyaaad31=1,xyzzyaaav1(set)-1
do xyzzyaaac31=1,xyzzyaaaw1(set)
xyzzyaaaa31=xyzzyaaaa31+xyzzyaaat31
xyzzyaaae31=xyzzyaaae31+1
c(xyzzyaaae31,xyzzyaaaa31+xyzzyaaar31+1*xyzzyaaau31)=-real(xyzzyaaad31&
&,dp)
c(xyzzyaaae31,xyzzyaaaa31+xyzzyaaar31+1*xyzzyaaas31+1*xyzzyaaau31)=xyz&
&zyaaaw31*real(xyzzyaaad31,dp)
enddo
enddo
xyzzyaaaa31=xyzzyaaaw1(set)-xyzzyaaat31
do xyzzyaaad31=1,xyzzyaaav1(set)-1
do xyzzyaaac31=1,xyzzyaaaw1(set)
xyzzyaaaa31=xyzzyaaaa31+xyzzyaaat31
xyzzyaaae31=xyzzyaaae31+1
c(xyzzyaaae31,xyzzyaaaa31+xyzzyaaar31+1*xyzzyaaau31)=-real(xyzzyaaad31&
&,dp)
enddo
enddo
else
xyzzyaaaa31=(xyzzyaaav1(set)-1)*xyzzyaaaw1(set)**2
do xyzzyaaac31=1,xyzzyaaaw1(set)
do xyzzyaaab31=1,xyzzyaaaw1(set)-1
xyzzyaaaa31=xyzzyaaaa31+1
xyzzyaaae31=xyzzyaaae31+1
c(xyzzyaaae31,xyzzyaaaa31+1*xyzzyaaas31)=1.d0
enddo
enddo
xyzzyaaaa31=xyzzyaaaw1(set)-1-xyzzyaaat31
do xyzzyaaad31=1,xyzzyaaav1(set)-1
do xyzzyaaac31=1,xyzzyaaaw1(set)
xyzzyaaaa31=xyzzyaaaa31+xyzzyaaat31
xyzzyaaae31=xyzzyaaae31+1
c(xyzzyaaae31,xyzzyaaaa31+xyzzyaaar31+1*xyzzyaaau31)=1.d0
enddo
enddo
endif
if(xyzzyaaae31/=no_eqns.and.am_master)call errstop('CONSTRUCT_C','Wron&
&g number of equations. Bug.')
endif
call reduced_echelon(no_eqns,no_vars,c)
end subroutine xyzzyaafm1
subroutine xyzzyaafn1
implicit none
integer xyzzyaaaa32,xyzzyaaab32
do xyzzyaaaa32=1,xyzzyaaai1
do xyzzyaaab32=1,xyzzyaaaz1(xyzzyaaaa32)
call xyzzyaafo1(xyzzyaaaa32,xyzzyaaab32)
enddo
enddo
end subroutine xyzzyaafn1
subroutine xyzzyaafo1(set,s)
implicit none
integer,intent(in) :: set,s
integer xyzzyaaaa33,xyzzyaaab33,xyzzyaaac33
real(dp) xyzzyaaad33
xyzzyaaac33=xyzzyaaav1(set)*xyzzyaaaw1(set)**2
call dcopy(xyzzyaaac33,xyzzyaabm1(1,1,1,s,set),1,xyzzyaacz1(1),1)
call dcopy(xyzzyaaac33,xyzzyaabn1(1,1,1,s,set),1,xyzzyaacz1(xyzzyaaac3&
&3+1),1)
do xyzzyaaaa33=1,xyzzyaacw1(s,set)
if(xyzzyaacy1(xyzzyaaaa33,s,set)>0)then
xyzzyaaad33=0.d0
do xyzzyaaab33=xyzzyaacy1(xyzzyaaaa33,s,set)+1,xyzzyaacx1(s,set)
xyzzyaaad33=xyzzyaaad33+xyzzyaada1(xyzzyaaaa33,xyzzyaaab33,s,set)*xyzz&
&yaacz1(xyzzyaaab33)
enddo
xyzzyaacz1(xyzzyaacy1(xyzzyaaaa33,s,set))=-xyzzyaaad33
endif
enddo
call dcopy(xyzzyaaac33,xyzzyaacz1(1),1,xyzzyaabm1(1,1,1,s,set),1)
call dcopy(xyzzyaaac33,xyzzyaacz1(xyzzyaaac33+1),1,xyzzyaabn1(1,1,1,s,&
&set),1)
end subroutine xyzzyaafo1
subroutine xyzzyaafp1
implicit none
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34,xyzzyaaad34,xyzzyaaae34,xy&
&zzyaaaf34,xyzzyaaag34,xyzzyaaah34,xyzzyaaai34
real(dp) xyzzyaaaj34
real(dp),parameter :: xyzzyaaak34=1.d-8
real(dp),allocatable :: xyzzyaaal34(:),xyzzyaaam34(:),xyzzyaaan34(:),x&
&yzzyaaao34(:),xyzzyaaap34(:),xyzzyaaaq34(:),xyzzyaaar34(:),xyzzyaaas3&
&4(:),xyzzyaaat34(:),xyzzyaaau34(:),xyzzyaaav34(:),xyzzyaaaw34(:),xyzz&
&yaaax34(:),xyzzyaaay34(:),xyzzyaaaz34(:)
logical :: xyzzyaaba34=.false.
if(.not.am_master)return
do xyzzyaaaa34=1,xyzzyaaai1
xyzzyaaaj34=real(xyzzyaaad1,dp)/xyzzyaabl1(xyzzyaaaa34)
xyzzyaaag34=2*xyzzyaaaw1(xyzzyaaaa34)-1
xyzzyaaah34=xyzzyaaav1(xyzzyaaaa34)+xyzzyaaaw1(xyzzyaaaa34)-1
xyzzyaaai34=xyzzyaaah34-1
do xyzzyaaab34=1,xyzzyaaaz1(xyzzyaaaa34)
allocate(xyzzyaaal34(xyzzyaaag34),xyzzyaaam34(xyzzyaaah34),xyzzyaaan34&
&(xyzzyaaah34),xyzzyaaao34(xyzzyaaah34),xyzzyaaap34(xyzzyaaah34),xyzzy&
&aaaq34(xyzzyaaah34),xyzzyaaar34(xyzzyaaah34),xyzzyaaas34(xyzzyaaai34)&
&,xyzzyaaat34(xyzzyaaai34),stat=xyzzyaaaf34)
call check_alloc(xyzzyaaaf34,'CHECK_CONSTRAINTS_PHI_THETA','1')
xyzzyaaba34=.false.
xyzzyaaal34=0.d0
xyzzyaaam34=0.d0
xyzzyaaan34=0.d0
xyzzyaaao34=0.d0
xyzzyaaap34=0.d0
xyzzyaaaq34=0.d0
xyzzyaaar34=0.d0
xyzzyaaas34=0.d0
xyzzyaaat34=0.d0
if(xyzzyaacs1(xyzzyaaab34,xyzzyaaaa34))then
do xyzzyaaad34=1,xyzzyaaaw1(xyzzyaaaa34)
do xyzzyaaac34=1,xyzzyaaaw1(xyzzyaaaa34)
xyzzyaaal34(xyzzyaaac34+xyzzyaaad34-1)=xyzzyaaal34(xyzzyaaac34+xyzzyaa&
&ad34-1)+xyzzyaabm1(xyzzyaaac34,xyzzyaaad34,2,xyzzyaaab34,xyzzyaaaa34)
enddo
enddo
endif
if(xyzzyaacu1(xyzzyaaaa34))then
do xyzzyaaae34=1,xyzzyaaav1(xyzzyaaaa34)
do xyzzyaaad34=1,xyzzyaaaw1(xyzzyaaaa34)
xyzzyaaao34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaao34(xyzzyaaae34+xyzzyaa&
&ad34-1)+xyzzyaabm1(1,xyzzyaaad34,xyzzyaaae34,xyzzyaaab34,xyzzyaaaa34)
xyzzyaaap34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaap34(xyzzyaaae34+xyzzyaa&
&ad34-1)+xyzzyaabm1(xyzzyaaad34,1,xyzzyaaae34,xyzzyaaab34,xyzzyaaaa34)
xyzzyaaaq34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaaq34(xyzzyaaae34+xyzzyaa&
&ad34-1)+xyzzyaabm1(2,xyzzyaaad34,xyzzyaaae34,xyzzyaaab34,xyzzyaaaa34)
xyzzyaaar34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaar34(xyzzyaaae34+xyzzyaa&
&ad34-1)+xyzzyaabm1(xyzzyaaad34,2,xyzzyaaae34,xyzzyaaab34,xyzzyaaaa34)
if(xyzzyaaae34>1)then
xyzzyaaas34(xyzzyaaae34+xyzzyaaad34-2)=xyzzyaaas34(xyzzyaaae34+xyzzyaa&
&ad34-2)+(xyzzyaaae34-1)*xyzzyaabm1(1,xyzzyaaad34,xyzzyaaae34,xyzzyaaa&
&b34,xyzzyaaaa34)
xyzzyaaat34(xyzzyaaae34+xyzzyaaad34-2)=xyzzyaaat34(xyzzyaaae34+xyzzyaa&
&ad34-2)+(xyzzyaaae34-1)*xyzzyaabm1(xyzzyaaad34,1,xyzzyaaae34,xyzzyaaa&
&b34,xyzzyaaaa34)
endif
enddo
enddo
else
do xyzzyaaae34=1,xyzzyaaav1(xyzzyaaaa34)
do xyzzyaaad34=1,xyzzyaaaw1(xyzzyaaaa34)
xyzzyaaam34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaam34(xyzzyaaae34+xyzzyaa&
&ad34-1)-xyzzyaaaj34*xyzzyaabm1(1,xyzzyaaad34,xyzzyaaae34,xyzzyaaab34,&
&xyzzyaaaa34)+xyzzyaabm1(2,xyzzyaaad34,xyzzyaaae34,xyzzyaaab34,xyzzyaa&
&aa34)
xyzzyaaan34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaan34(xyzzyaaae34+xyzzyaa&
&ad34-1)-xyzzyaaaj34*xyzzyaabm1(xyzzyaaad34,1,xyzzyaaae34,xyzzyaaab34,&
&xyzzyaaaa34)+xyzzyaabm1(xyzzyaaad34,2,xyzzyaaae34,xyzzyaaab34,xyzzyaa&
&aa34)
enddo
enddo
endif
if(any(abs(xyzzyaaal34)>xyzzyaaak34))then
call wout('Electron-electron cusp condition (phi) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaal34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaao34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (phi_1a) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaao34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaap34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (phi_1b) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaap34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaaq34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (phi_2a) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaaq34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaar34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (phi_2b) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaar34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaas34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (phi_3a) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaas34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaat34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (phi_3b) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaat34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaam34)>xyzzyaaak34))then
call wout('Electron-nucleus PP cusp condition (phi_1a) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaam34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaan34)>xyzzyaaak34))then
call wout('Electron-nucleus PP cusp condition (phi_1b) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaan34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
deallocate(xyzzyaaal34,xyzzyaaam34,xyzzyaaan34,xyzzyaaao34,xyzzyaaap34&
&,xyzzyaaaq34,xyzzyaaar34,xyzzyaaas34,xyzzyaaat34)
allocate(xyzzyaaau34(xyzzyaaag34),xyzzyaaav34(xyzzyaaah34),xyzzyaaaw34&
&(xyzzyaaah34),xyzzyaaax34(xyzzyaaah34),xyzzyaaay34(xyzzyaaai34),xyzzy&
&aaaz34(xyzzyaaai34),stat=xyzzyaaaf34)
call check_alloc(xyzzyaaaf34,'CHECK_CUSP_PHI_THETA','2')
xyzzyaaau34=0.d0
xyzzyaaav34=0.d0
xyzzyaaaw34=0.d0
xyzzyaaax34=0.d0
xyzzyaaay34=0.d0
xyzzyaaaz34=0.d0
do xyzzyaaad34=1,xyzzyaaaw1(xyzzyaaaa34)
do xyzzyaaac34=1,xyzzyaaaw1(xyzzyaaaa34)
xyzzyaaau34(xyzzyaaac34+xyzzyaaad34-1)=xyzzyaaau34(xyzzyaaac34+xyzzyaa&
&ad34-1)+xyzzyaabn1(xyzzyaaac34,xyzzyaaad34,2,xyzzyaaab34,xyzzyaaaa34)
enddo
enddo
if(xyzzyaacu1(xyzzyaaaa34))then
do xyzzyaaae34=1,xyzzyaaav1(xyzzyaaaa34)
do xyzzyaaad34=1,xyzzyaaaw1(xyzzyaaaa34)
xyzzyaaaw34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaaw34(xyzzyaaae34+xyzzyaa&
&ad34-1)+xyzzyaabn1(1,xyzzyaaad34,xyzzyaaae34,xyzzyaaab34,xyzzyaaaa34)
xyzzyaaax34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaax34(xyzzyaaae34+xyzzyaa&
&ad34-1)+xyzzyaabn1(2,xyzzyaaad34,xyzzyaaae34,xyzzyaaab34,xyzzyaaaa34)
xyzzyaaav34(xyzzyaaae34+xyzzyaaad34-1)=xyzzyaaav34(xyzzyaaae34+xyzzyaa&
&ad34-1)-xyzzyaaaj34*xyzzyaabn1(xyzzyaaad34,1,xyzzyaaae34,xyzzyaaab34,&
&xyzzyaaaa34)+xyzzyaabn1(xyzzyaaad34,2,xyzzyaaae34,xyzzyaaab34,xyzzyaa&
&aa34)
if(xyzzyaaae34>1)then
xyzzyaaay34(xyzzyaaae34+xyzzyaaad34-2)=xyzzyaaay34(xyzzyaaae34+xyzzyaa&
&ad34-2)+(xyzzyaaae34-1)*xyzzyaabn1(1,xyzzyaaad34,xyzzyaaae34,xyzzyaaa&
&b34,xyzzyaaaa34)
xyzzyaaaz34(xyzzyaaae34+xyzzyaaad34-2)=xyzzyaaaz34(xyzzyaaae34+xyzzyaa&
&ad34-2)+(xyzzyaaae34-1)*xyzzyaabn1(xyzzyaaad34,1,xyzzyaaae34,xyzzyaaa&
&b34,xyzzyaaaa34)
endif
enddo
enddo
else
do xyzzyaaae34=1,xyzzyaaav1(xyzzyaaaa34)
do xyzzyaaac34=1,xyzzyaaaw1(xyzzyaaaa34)
xyzzyaaav34(xyzzyaaae34+xyzzyaaac34-1)=xyzzyaaav34(xyzzyaaae34+xyzzyaa&
&ac34-1)-xyzzyaaaj34*xyzzyaabn1(xyzzyaaac34,1,xyzzyaaae34,xyzzyaaab34,&
&xyzzyaaaa34)+xyzzyaabn1(xyzzyaaac34,2,xyzzyaaae34,xyzzyaaab34,xyzzyaa&
&aa34)
enddo
enddo
endif
if(any(abs(xyzzyaaau34)>xyzzyaaak34))then
call wout('Electron-electron cusp condition (theta) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaau34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaaw34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (theta_1) not verified.'&
&)
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaaw34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaax34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (theta_2a) not verified.&
&')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaax34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaav34)>xyzzyaaak34))then
if(xyzzyaacu1(xyzzyaaaa34))then
call wout('Electron-nucleus cusp condition (theta_2b) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
else
call wout('Electron-nucleus cusp condition (theta_1a) not verified.')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
endif
call wout('',xyzzyaaav34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaay34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (theta_3a) not verified.&
&')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaay34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
if(any(abs(xyzzyaaaz34)>xyzzyaaak34))then
call wout('Electron-nucleus AE cusp condition (theta_3b) not verified.&
&')
call wout('Sums for set '//trim(i2s(xyzzyaaaa34))//', spin type '//tri&
&m(i2s(xyzzyaaab34))//':')
call wout('',xyzzyaaaz34,rfmt='(e21.12)')
xyzzyaaba34=.true.
endif
deallocate(xyzzyaaau34,xyzzyaaav34,xyzzyaaaw34,xyzzyaaax34,xyzzyaaay34&
&,xyzzyaaaz34)
if(xyzzyaacv1(xyzzyaaaa34))call xyzzyaafq1(xyzzyaaab34,xyzzyaaaa34,xyz&
&zyaaba34)
enddo
enddo
if(xyzzyaaba34)call errstop('CHECK_CUSP_PHI_THETA','Stopping.')
end subroutine xyzzyaafp1
subroutine xyzzyaafq1(s,set,die)
implicit none
integer,intent(in) :: s,set
logical,intent(inout) :: die
integer ii,jj,xyzzyaaaa35,xyzzyaaab35,xyzzyaaac35,ion,xyzzyaaad35
real(dp) rij,xyzzyaaae35,xyzzyaaaf35,xyzzyaaag35(4),xyzzyaaah35(4),xyz&
&zyaaai35(4),xyzzyaaaj35,xyzzyaaak35(3),xyzzyaaal35(3),xyzzyaaam35,xyz&
&zyaaan35(3),xyzzyaaao35(3),xyzzyaaap35(3)
real(dp),parameter :: xyzzyaaaq35=1.d-10
ii=0
jj=0
do xyzzyaaaa35=1,netot
do xyzzyaaab35=1,netot
if(xyzzyaaaa35==xyzzyaaab35)cycle
if(xyzzyaado1(xyzzyaaaa35,xyzzyaaab35,set)==s)then
ii=xyzzyaaaa35
jj=xyzzyaaab35
exit
endif
enddo
if(ii/=0)exit
enddo
if(ii==0)return
xyzzyaaac35=which_ee(ii,jj)
ion=0
do xyzzyaaad35=1,nitot
if(xyzzyaaat1(xyzzyaaad35)==set)then
ion=xyzzyaaad35
exit
endif
enddo
if(ion==0)return
xyzzyaaag35(1:3)=xyzzyaabl1(set)*(/.22d0,-.11d0,.01d0/)
xyzzyaaah35(1:3)=xyzzyaabl1(set)*(/-.05d0,.1d0,.3d0/)
xyzzyaaai35(1:3)=xyzzyaaag35(1:3)-xyzzyaaah35(1:3)
xyzzyaaae35=sqrt(sum(xyzzyaaag35(1:3)**2))
xyzzyaaag35(4)=xyzzyaaae35
xyzzyaaaf35=sqrt(sum(xyzzyaaah35(1:3)**2))
xyzzyaaah35(4)=xyzzyaaaf35
rij=sqrt(sum(xyzzyaaai35(1:3)**2))
xyzzyaaai35(4)=rij
xyzzyaacj1=0.d0
xyzzyaack1=0.d0
xyzzyaacm1=0.d0
xyzzyaacn1=0.d0
call xyzzyaagp1(rij,xyzzyaaav1(set),xyzzyaacj1(1,xyzzyaaac35),xyzzyaac&
&k1(1,xyzzyaaac35))
call xyzzyaagp1(xyzzyaaae35,xyzzyaaaw1(set),xyzzyaacm1(1,ion,ii),xyzzy&
&aacn1(1,ion,ii))
call xyzzyaagp1(xyzzyaaaf35,xyzzyaaaw1(set),xyzzyaacm1(1,ion,jj),xyzzy&
&aacn1(1,ion,jj))
call xyzzyaagf1(ii,jj,ion,xyzzyaaai35,xyzzyaaag35,xyzzyaaah35,xyzzyaaa&
&j35,xyzzyaaak35,xyzzyaaal35,xyzzyaaam35,xyzzyaaan35,xyzzyaaao35)
xyzzyaaap35(1)=xyzzyaaai35(3)*xyzzyaaak35(2)+xyzzyaaag35(3)*xyzzyaaan3&
&5(2)-xyzzyaaai35(2)*xyzzyaaak35(3)-xyzzyaaag35(2)*xyzzyaaan35(3)
xyzzyaaap35(2)=xyzzyaaai35(1)*xyzzyaaak35(3)+xyzzyaaag35(1)*xyzzyaaan3&
&5(3)-xyzzyaaai35(3)*xyzzyaaak35(1)-xyzzyaaag35(3)*xyzzyaaan35(1)
xyzzyaaap35(3)=xyzzyaaai35(2)*xyzzyaaak35(1)+xyzzyaaag35(2)*xyzzyaaan3&
&5(1)-xyzzyaaai35(1)*xyzzyaaak35(2)-xyzzyaaag35(1)*xyzzyaaan35(2)
if(any(abs(xyzzyaaap35)>xyzzyaaaq35))then
call wout('The irrotational-backflow set of constraints is not being o&
&beyed.')
call wout('Current set: s='//trim(i2s(s))//', set='//trim(i2s(set))//'&
&. Test particles: ii='//trim(i2s(ii))//', jj='//trim(i2s(jj))//   ', &
&ion='//trim(i2s(ion)))
call wout('Distances: riI='//trim(d2s(xyzzyaaae35,10))//', rjI='//trim&
&(d2s(xyzzyaaaf35,10))//', rij='//trim(d2s(rij,10)))
call wout('Tolerance: '//trim(d2s(xyzzyaaaq35,10)))
call wout('Curl of the backflow displacement at test point:')
call wout()
call wout('',xyzzyaaap35)
call wout()
die=.true.
endif
end subroutine xyzzyaafq1
subroutine xyzzyaafr1(ii,rnew,eevecs1,eivecs,rold,eevecs1_old,eivecs_o&
&ld,bf_connect,bf_m,bf_rmap,bf_m_change,bf_rmap_change,bf_x,bf_dx)
implicit none
integer,intent(in) :: ii
integer,intent(out) :: bf_m_change(nspin),bf_rmap_change(nemax,nspin)
integer,intent(inout) :: bf_m(netot),bf_rmap(netot,netot)
real(dp),intent(in) :: rnew(3),eevecs1(4,netot),eivecs(4,nitot,netot),&
&rold(3),eevecs1_old(4,netot),eivecs_old(4,nitot,netot)
real(dp),intent(inout) :: bf_x(3,netot)
real(dp),intent(inout),optional :: bf_dx(3,3,netot,netot)
logical,intent(inout) :: bf_connect(netot,netot)
logical xyzzyaaaa36
call timer('ONEELEC_BACKFLOW_R2X',.true.)
xyzzyaaaa36=present(bf_dx)
bf_x=-bf_x
if(xyzzyaaaa36)bf_dx=-bf_dx
if(xyzzyaaaa36)then
call xyzzyaafs1(ii,eevecs1_old,eivecs_old,bf_m(ii),bf_rmap(1,ii),bf_x,&
&bf_dx)
else
call xyzzyaafs1(ii,eevecs1_old,eivecs_old,bf_m(ii),bf_rmap(1,ii),bf_x)
endif
bf_x=-bf_x
if(xyzzyaaaa36)bf_dx=-bf_dx
call xyzzyaafu1(ii,eevecs1,eivecs,bf_m_change,bf_rmap_change,bf_connec&
&t,bf_m,bf_rmap)
if(xyzzyaaaa36)then
call xyzzyaafs1(ii,eevecs1,eivecs,bf_m(ii),bf_rmap(1,ii),bf_x,bf_dx)
else
call xyzzyaafs1(ii,eevecs1,eivecs,bf_m(ii),bf_rmap(1,ii),bf_x)
endif
bf_x(:,ii)=bf_x(:,ii)+rnew-rold
call timer('ONEELEC_BACKFLOW_R2X',.false.)
end subroutine xyzzyaafr1
subroutine xyzzyaafs1(ii,eevecs1,eivecs,bf_m_ii,bf_rmap_ii,bf_x,bf_dx)
implicit none
integer,intent(in) :: ii
integer,intent(in) :: bf_m_ii,bf_rmap_ii(bf_m_ii)
real(dp),intent(in) :: eevecs1(4,netot),eivecs(4,nitot,netot)
real(dp),intent(inout) :: bf_x(3,netot)
real(dp),intent(inout),optional :: bf_dx(3,3,netot,netot)
integer xyzzyaaaa37,xyzzyaaab37,xyzzyaaac37,xyzzyaaad37,xyzzyaaae37,xy&
&zzyaaaf37,xyzzyaaag37,xyzzyaaah37,xyzzyaaai37,xyzzyaaaj37(netot)
real(dp) xyzzyaaak37,xyzzyaaal37,xyzzyaaam37,xyzzyaaan37,xyzzyaaao37,x&
&yzzyaaap37(4),xyzzyaaaq37(4),xyzzyaaar37(4),xyzzyaaas37(4),xyzzyaaat3&
&7,xyzzyaaau37,xyzzyaaav37,xyzzyaaaw37(3),xyzzyaaax37(3),xyzzyaaay37(3&
&),xyzzyaaaz37(3),xyzzyaaba37(3),xyzzyaabb37(3),xyzzyaabc37,xyzzyaabd3&
&7,xyzzyaabe37(3),xyzzyaabf37(3),xyzzyaabg37(netot),xyzzyaabh37(3,neto&
&t)
logical xyzzyaabi37
xyzzyaabi37=present(bf_dx)
call xyzzyaafx1(ii,xyzzyaabi37,eevecs1,eivecs,bf_m_ii,bf_rmap_ii)
if(xyzzyaadb1)then
if(xyzzyaabi37)then
call xyzzyaagk1(eivecs,xyzzyaaaj37,xyzzyaabg37,xyzzyaabh37)
else
call xyzzyaagk1(eivecs,xyzzyaaaj37,xyzzyaabg37)
endif
endif
if(.not.xyzzyaabi37)then
if(.not.xyzzyaadb1)then
if(xyzzyaaaa1)then
do xyzzyaaaa37=1,bf_m_ii
xyzzyaaab37=bf_rmap_ii(xyzzyaaaa37)
if(xyzzyaaab37==ii)cycle
xyzzyaaac37=xyzzyaadm1(ii,xyzzyaaab37)
xyzzyaaao37=eevecs1(4,xyzzyaaab37)
if(xyzzyaaao37>=xyzzyaabh1(xyzzyaaac37))cycle
call xyzzyaaga1(ii,xyzzyaaab37,xyzzyaaao37,xyzzyaaak37)
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaak37*eevecs1(xyzzyaaae37,xyzzyaaab37)
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaat37
bf_x(xyzzyaaae37,xyzzyaaab37)=bf_x(xyzzyaaae37,xyzzyaaab37)-xyzzyaaat3&
&7
enddo
enddo
endif
if(xyzzyaaab1)then
do xyzzyaaag37=1,nitot
xyzzyaaad37=xyzzyaaan1(xyzzyaaag37)
if(xyzzyaaad37==0)cycle
if(eivecs(4,xyzzyaaag37,ii)>=xyzzyaabj1(xyzzyaaad37))cycle
call xyzzyaagd1(ii,xyzzyaaag37,xyzzyaaal37)
do xyzzyaaae37=1,dimensionality
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaal37*eivecs(xyzzyaaae&
&37,xyzzyaaag37,ii)
enddo
enddo
endif
if(xyzzyaaac1)then
do xyzzyaaaa37=1,bf_m_ii
xyzzyaaab37=bf_rmap_ii(xyzzyaaaa37)
if(xyzzyaaab37==ii)cycle
if(eevecs1(4,xyzzyaaab37)>=2*xyzzyaacp1)cycle
do xyzzyaaag37=1,nitot
xyzzyaaad37=xyzzyaaat1(xyzzyaaag37)
if(xyzzyaaad37==0)cycle
if(eivecs(4,xyzzyaaag37,ii)>=xyzzyaabl1(xyzzyaaad37).or.eivecs(4,xyzzy&
&aaag37,xyzzyaaab37)>=xyzzyaabl1(xyzzyaaad37))cycle
call xyzzyaagg1(ii,xyzzyaaab37,xyzzyaaag37,xyzzyaaam37,xyzzyaaan37)
do xyzzyaaae37=1,dimensionality
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaam37*eevecs1(xyzzyaaa&
&e37,xyzzyaaab37)+xyzzyaaan37*eivecs(xyzzyaaae37,xyzzyaaag37,ii)
enddo
call xyzzyaagg1(xyzzyaaab37,ii,xyzzyaaag37,xyzzyaaam37,xyzzyaaan37)
do xyzzyaaae37=1,dimensionality
bf_x(xyzzyaaae37,xyzzyaaab37)=bf_x(xyzzyaaae37,xyzzyaaab37)-xyzzyaaam3&
&7*eevecs1(xyzzyaaae37,xyzzyaaab37)+xyzzyaaan37*eivecs(xyzzyaaae37,xyz&
&zyaaag37,xyzzyaaab37)
enddo
enddo
enddo
endif
else
if(xyzzyaaaa1)then
xyzzyaabc37=xyzzyaabg37(ii)
do xyzzyaaaa37=1,bf_m_ii
xyzzyaaab37=bf_rmap_ii(xyzzyaaaa37)
if(xyzzyaaab37==ii)cycle
xyzzyaaac37=xyzzyaadm1(ii,xyzzyaaab37)
xyzzyaaao37=eevecs1(4,xyzzyaaab37)
if(xyzzyaaao37>=xyzzyaabh1(xyzzyaaac37))cycle
call xyzzyaaga1(ii,xyzzyaaab37,xyzzyaaao37,xyzzyaaak37)
xyzzyaabd37=xyzzyaabg37(xyzzyaaab37)
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaak37*eevecs1(xyzzyaaae37,xyzzyaaab37)
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaat37*xyzzyaabc37
bf_x(xyzzyaaae37,xyzzyaaab37)=bf_x(xyzzyaaae37,xyzzyaaab37)-xyzzyaaat3&
&7*xyzzyaabd37
enddo
enddo
endif
if(xyzzyaaab1)then
xyzzyaaah37=xyzzyaaaj37(ii)
do xyzzyaaag37=1,nitot
xyzzyaaad37=xyzzyaaan1(xyzzyaaag37)
if(xyzzyaaad37==0)cycle
if(eivecs(4,xyzzyaaag37,ii)>=xyzzyaabj1(xyzzyaaad37))cycle
xyzzyaabc37=1.d0
if(xyzzyaaah37/=xyzzyaaag37)xyzzyaabc37=xyzzyaabg37(ii)
call xyzzyaagd1(ii,xyzzyaaag37,xyzzyaaal37)
do xyzzyaaae37=1,dimensionality
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaal37*eivecs(xyzzyaaae&
&37,xyzzyaaag37,ii)*xyzzyaabc37
enddo
enddo
endif
if(xyzzyaaac1)then
xyzzyaaah37=xyzzyaaaj37(ii)
do xyzzyaaaa37=1,bf_m_ii
xyzzyaaab37=bf_rmap_ii(xyzzyaaaa37)
if(xyzzyaaab37==ii)cycle
if(eevecs1(4,xyzzyaaab37)>=2*xyzzyaacp1)cycle
xyzzyaaai37=xyzzyaaaj37(xyzzyaaab37)
do xyzzyaaag37=1,nitot
xyzzyaaad37=xyzzyaaat1(xyzzyaaag37)
if(xyzzyaaad37==0)cycle
if(eivecs(4,xyzzyaaag37,ii)>=xyzzyaabl1(xyzzyaaad37).or.eivecs(4,xyzzy&
&aaag37,xyzzyaaab37)>=xyzzyaabl1(xyzzyaaad37))cycle
xyzzyaabc37=1.d0
if(xyzzyaaah37/=xyzzyaaag37)xyzzyaabc37=xyzzyaabg37(ii)
xyzzyaabd37=1.d0
if(xyzzyaaai37/=xyzzyaaag37)xyzzyaabd37=xyzzyaabg37(xyzzyaaab37)
call xyzzyaagg1(ii,xyzzyaaab37,xyzzyaaag37,xyzzyaaam37,xyzzyaaan37)
do xyzzyaaae37=1,dimensionality
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+(xyzzyaaam37*eevecs1(xyzzyaa&
&ae37,xyzzyaaab37)+xyzzyaaan37*eivecs(xyzzyaaae37,xyzzyaaag37,ii))*xyz&
&zyaabc37
enddo
call xyzzyaagg1(xyzzyaaab37,ii,xyzzyaaag37,xyzzyaaam37,xyzzyaaan37)
do xyzzyaaae37=1,dimensionality
bf_x(xyzzyaaae37,xyzzyaaab37)=bf_x(xyzzyaaae37,xyzzyaaab37)+(-xyzzyaaa&
&m37*eevecs1(xyzzyaaae37,xyzzyaaab37)+xyzzyaaan37*eivecs(xyzzyaaae37,x&
&yzzyaaag37,xyzzyaaab37))*xyzzyaabd37
enddo
enddo
enddo
endif
endif
else
if(.not.xyzzyaadb1)then
if(xyzzyaaaa1)then
do xyzzyaaaa37=1,bf_m_ii
xyzzyaaab37=bf_rmap_ii(xyzzyaaaa37)
if(xyzzyaaab37==ii)cycle
xyzzyaaac37=xyzzyaadm1(ii,xyzzyaaab37)
xyzzyaaap37=eevecs1(:,xyzzyaaab37)
if(xyzzyaaap37(4)>=xyzzyaabh1(xyzzyaaac37))cycle
call xyzzyaafz1(ii,xyzzyaaab37,xyzzyaaap37,xyzzyaaak37,xyzzyaaaw37)
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaak37*xyzzyaaap37(xyzzyaaae37)
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaat37
bf_x(xyzzyaaae37,xyzzyaaab37)=bf_x(xyzzyaaae37,xyzzyaaab37)-xyzzyaaat3&
&7
do xyzzyaaaf37=1,dimensionality
xyzzyaaau37=xyzzyaaap37(xyzzyaaae37)*xyzzyaaaw37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,xyzzyaaab37)=bf_dx(xyzzyaaae37,xyzzya&
&aaf37,ii,xyzzyaaab37)-xyzzyaaau37
bf_dx(xyzzyaaae37,xyzzyaaaf37,xyzzyaaab37,ii)=bf_dx(xyzzyaaae37,xyzzya&
&aaf37,xyzzyaaab37,ii)-xyzzyaaau37
bf_dx(xyzzyaaae37,xyzzyaaaf37,xyzzyaaab37,xyzzyaaab37)=bf_dx(xyzzyaaae&
&37,xyzzyaaaf37,xyzzyaaab37,xyzzyaaab37)+xyzzyaaau37
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,&
&ii)+xyzzyaaau37
enddo
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,xyzzyaaab37)=bf_dx(xyzzyaaae37,xyzzya&
&aae37,ii,xyzzyaaab37)-xyzzyaaak37
bf_dx(xyzzyaaae37,xyzzyaaae37,xyzzyaaab37,ii)=bf_dx(xyzzyaaae37,xyzzya&
&aae37,xyzzyaaab37,ii)-xyzzyaaak37
bf_dx(xyzzyaaae37,xyzzyaaae37,xyzzyaaab37,xyzzyaaab37)=bf_dx(xyzzyaaae&
&37,xyzzyaaae37,xyzzyaaab37,xyzzyaaab37)+xyzzyaaak37
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaae37,ii,&
&ii)+xyzzyaaak37
enddo
enddo
endif
if(xyzzyaaab1)then
do xyzzyaaag37=1,nitot
xyzzyaaad37=xyzzyaaan1(xyzzyaaag37)
if(xyzzyaaad37==0)cycle
xyzzyaaar37=eivecs(:,xyzzyaaag37,ii)
if(xyzzyaaar37(4)>=xyzzyaabj1(xyzzyaaad37))cycle
call xyzzyaagc1(ii,xyzzyaaag37,xyzzyaaar37,xyzzyaaal37,xyzzyaaax37)
do xyzzyaaae37=1,dimensionality
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaal37*xyzzyaaar37(xyzz&
&yaaae37)
do xyzzyaaaf37=1,dimensionality
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,&
&ii)+xyzzyaaax37(xyzzyaaae37)*xyzzyaaar37(xyzzyaaaf37)
enddo
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaae37,ii,&
&ii)+xyzzyaaal37
enddo
enddo
endif
if(xyzzyaaac1)then
do xyzzyaaaa37=1,bf_m_ii
xyzzyaaab37=bf_rmap_ii(xyzzyaaaa37)
if(xyzzyaaab37==ii)cycle
xyzzyaaap37=eevecs1(:,xyzzyaaab37)
if(xyzzyaaap37(4)>=2*xyzzyaacp1)cycle
xyzzyaaaq37(1:3)=-xyzzyaaap37(1:3)
xyzzyaaaq37(4)=xyzzyaaap37(4)
do xyzzyaaag37=1,nitot
xyzzyaaad37=xyzzyaaat1(xyzzyaaag37)
if(xyzzyaaad37==0)cycle
xyzzyaaar37=eivecs(:,xyzzyaaag37,ii)
xyzzyaaas37=eivecs(:,xyzzyaaag37,xyzzyaaab37)
if(xyzzyaaar37(4)>=xyzzyaabl1(xyzzyaaad37).or.xyzzyaaas37(4)>=xyzzyaab&
&l1(xyzzyaaad37))cycle
call xyzzyaagf1(ii,xyzzyaaab37,xyzzyaaag37,xyzzyaaap37,xyzzyaaar37,xyz&
&zyaaas37,xyzzyaaam37,xyzzyaaay37,xyzzyaaaz37,xyzzyaaan37,xyzzyaaba37,&
&xyzzyaabb37)
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaam37*xyzzyaaap37(xyzzyaaae37)+xyzzyaaan37*xyzzyaaar&
&37(xyzzyaaae37)
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaat37
do xyzzyaaaf37=1,dimensionality
xyzzyaaau37=xyzzyaaap37(xyzzyaaae37)*xyzzyaaay37(xyzzyaaaf37)+xyzzyaaa&
&r37(xyzzyaaae37)*xyzzyaaba37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,&
&ii)+xyzzyaaau37
xyzzyaaav37=xyzzyaaap37(xyzzyaaae37)*xyzzyaaaz37(xyzzyaaaf37)+xyzzyaaa&
&r37(xyzzyaaae37)*xyzzyaabb37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,xyzzyaaab37)=bf_dx(xyzzyaaae37,xyzzya&
&aaf37,ii,xyzzyaaab37)+xyzzyaaav37
enddo
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaae37,ii,&
&ii)+xyzzyaaam37+xyzzyaaan37
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,xyzzyaaab37)=bf_dx(xyzzyaaae37,xyzzya&
&aae37,ii,xyzzyaaab37)-xyzzyaaam37
enddo
call xyzzyaagf1(xyzzyaaab37,ii,xyzzyaaag37,xyzzyaaaq37,xyzzyaaas37,xyz&
&zyaaar37,xyzzyaaam37,xyzzyaaaz37,xyzzyaaay37,xyzzyaaan37,xyzzyaabb37,&
&xyzzyaaba37)
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaam37*xyzzyaaaq37(xyzzyaaae37)+xyzzyaaan37*xyzzyaaas&
&37(xyzzyaaae37)
bf_x(xyzzyaaae37,xyzzyaaab37)=bf_x(xyzzyaaae37,xyzzyaaab37)+xyzzyaaat3&
&7
do xyzzyaaaf37=1,dimensionality
xyzzyaaau37=xyzzyaaaq37(xyzzyaaae37)*xyzzyaaaz37(xyzzyaaaf37)+xyzzyaaa&
&s37(xyzzyaaae37)*xyzzyaabb37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,xyzzyaaab37,xyzzyaaab37)=bf_dx(xyzzyaaae&
&37,xyzzyaaaf37,xyzzyaaab37,xyzzyaaab37)+xyzzyaaau37
xyzzyaaav37=xyzzyaaaq37(xyzzyaaae37)*xyzzyaaay37(xyzzyaaaf37)+xyzzyaaa&
&s37(xyzzyaaae37)*xyzzyaaba37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,xyzzyaaab37,ii)=bf_dx(xyzzyaaae37,xyzzya&
&aaf37,xyzzyaaab37,ii)+xyzzyaaav37
enddo
bf_dx(xyzzyaaae37,xyzzyaaae37,xyzzyaaab37,xyzzyaaab37)=bf_dx(xyzzyaaae&
&37,xyzzyaaae37,xyzzyaaab37,xyzzyaaab37)+xyzzyaaam37+xyzzyaaan37
bf_dx(xyzzyaaae37,xyzzyaaae37,xyzzyaaab37,ii)=bf_dx(xyzzyaaae37,xyzzya&
&aae37,xyzzyaaab37,ii)-xyzzyaaam37
enddo
enddo
enddo
endif
else
if(xyzzyaaaa1)then
xyzzyaabc37=xyzzyaabg37(ii)
xyzzyaabe37=xyzzyaabh37(:,ii)
do xyzzyaaaa37=1,bf_m_ii
xyzzyaaab37=bf_rmap_ii(xyzzyaaaa37)
if(xyzzyaaab37==ii)cycle
xyzzyaaac37=xyzzyaadm1(ii,xyzzyaaab37)
xyzzyaaap37=eevecs1(:,xyzzyaaab37)
if(xyzzyaaap37(4)>=xyzzyaabh1(xyzzyaaac37))cycle
call xyzzyaafz1(ii,xyzzyaaab37,xyzzyaaap37,xyzzyaaak37,xyzzyaaaw37)
xyzzyaabd37=xyzzyaabg37(xyzzyaaab37)
xyzzyaabf37=xyzzyaabh37(:,xyzzyaaab37)
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaak37*xyzzyaaap37(xyzzyaaae37)
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaat37*xyzzyaabc37
bf_x(xyzzyaaae37,xyzzyaaab37)=bf_x(xyzzyaaae37,xyzzyaaab37)-xyzzyaaat3&
&7*xyzzyaabd37
do xyzzyaaaf37=1,dimensionality
xyzzyaaau37=xyzzyaaap37(xyzzyaaae37)*xyzzyaaaw37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,xyzzyaaab37)=bf_dx(xyzzyaaae37,xyzzya&
&aaf37,ii,xyzzyaaab37)-xyzzyaaau37*xyzzyaabc37
bf_dx(xyzzyaaae37,xyzzyaaaf37,xyzzyaaab37,ii)=bf_dx(xyzzyaaae37,xyzzya&
&aaf37,xyzzyaaab37,ii)-xyzzyaaau37*xyzzyaabd37
bf_dx(xyzzyaaae37,xyzzyaaaf37,xyzzyaaab37,xyzzyaaab37)=bf_dx(xyzzyaaae&
&37,xyzzyaaaf37,xyzzyaaab37,xyzzyaaab37)+xyzzyaaau37*xyzzyaabd37-xyzzy&
&aaat37*xyzzyaabf37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,&
&ii)+xyzzyaaau37*xyzzyaabc37+xyzzyaaat37*xyzzyaabe37(xyzzyaaaf37)
enddo
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,xyzzyaaab37)=bf_dx(xyzzyaaae37,xyzzya&
&aae37,ii,xyzzyaaab37)-xyzzyaaak37*xyzzyaabc37
bf_dx(xyzzyaaae37,xyzzyaaae37,xyzzyaaab37,ii)=bf_dx(xyzzyaaae37,xyzzya&
&aae37,xyzzyaaab37,ii)-xyzzyaaak37*xyzzyaabd37
bf_dx(xyzzyaaae37,xyzzyaaae37,xyzzyaaab37,xyzzyaaab37)=bf_dx(xyzzyaaae&
&37,xyzzyaaae37,xyzzyaaab37,xyzzyaaab37)+xyzzyaaak37*xyzzyaabd37
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaae37,ii,&
&ii)+xyzzyaaak37*xyzzyaabc37
enddo
enddo
endif
if(xyzzyaaab1)then
xyzzyaaah37=xyzzyaaaj37(ii)
do xyzzyaaag37=1,nitot
xyzzyaaad37=xyzzyaaan1(xyzzyaaag37)
if(xyzzyaaad37==0)cycle
xyzzyaaar37=eivecs(:,xyzzyaaag37,ii)
if(xyzzyaaar37(4)>=xyzzyaabj1(xyzzyaaad37))cycle
call xyzzyaagc1(ii,xyzzyaaag37,xyzzyaaar37,xyzzyaaal37,xyzzyaaax37)
xyzzyaabc37=1.d0
xyzzyaabe37=0.d0
if(xyzzyaaah37/=xyzzyaaag37)then
xyzzyaabc37=xyzzyaabg37(ii)
xyzzyaabe37=xyzzyaabh37(:,ii)
endif
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaal37*xyzzyaaar37(xyzzyaaae37)
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaat37*xyzzyaabc37
do xyzzyaaaf37=1,dimensionality
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,&
&ii)+xyzzyaaax37(xyzzyaaae37)*xyzzyaaar37(xyzzyaaaf37)*xyzzyaabc37+xyz&
&zyaaat37*xyzzyaabe37(xyzzyaaaf37)
enddo
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaae37,ii,&
&ii)+xyzzyaaal37*xyzzyaabc37
enddo
enddo
endif
if(xyzzyaaac1)then
xyzzyaaah37=xyzzyaaaj37(ii)
do xyzzyaaaa37=1,bf_m_ii
xyzzyaaab37=bf_rmap_ii(xyzzyaaaa37)
if(xyzzyaaab37==ii)cycle
xyzzyaaap37=eevecs1(:,xyzzyaaab37)
if(xyzzyaaap37(4)>=2*xyzzyaacp1)cycle
xyzzyaaaq37(1:3)=-xyzzyaaap37(1:3)
xyzzyaaaq37(4)=xyzzyaaap37(4)
xyzzyaaai37=xyzzyaaaj37(xyzzyaaab37)
do xyzzyaaag37=1,nitot
xyzzyaaad37=xyzzyaaat1(xyzzyaaag37)
if(xyzzyaaad37==0)cycle
xyzzyaaar37=eivecs(:,xyzzyaaag37,ii)
xyzzyaaas37=eivecs(:,xyzzyaaag37,xyzzyaaab37)
if(xyzzyaaar37(4)>=xyzzyaabl1(xyzzyaaad37).or.xyzzyaaas37(4)>=xyzzyaab&
&l1(xyzzyaaad37))cycle
call xyzzyaagf1(ii,xyzzyaaab37,xyzzyaaag37,xyzzyaaap37,xyzzyaaar37,xyz&
&zyaaas37,xyzzyaaam37,xyzzyaaay37,xyzzyaaaz37,xyzzyaaan37,xyzzyaaba37,&
&xyzzyaabb37)
xyzzyaabc37=1.d0
xyzzyaabe37=0.d0
if(xyzzyaaah37/=xyzzyaaag37)then
xyzzyaabc37=xyzzyaabg37(ii)
xyzzyaabe37=xyzzyaabh37(:,ii)
endif
xyzzyaabd37=1.d0
xyzzyaabf37=0.d0
if(xyzzyaaai37/=xyzzyaaag37)then
xyzzyaabd37=xyzzyaabg37(xyzzyaaab37)
xyzzyaabf37=xyzzyaabh37(:,xyzzyaaab37)
endif
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaam37*xyzzyaaap37(xyzzyaaae37)+xyzzyaaan37*xyzzyaaar&
&37(xyzzyaaae37)
bf_x(xyzzyaaae37,ii)=bf_x(xyzzyaaae37,ii)+xyzzyaaat37*xyzzyaabc37
do xyzzyaaaf37=1,dimensionality
xyzzyaaau37=xyzzyaaap37(xyzzyaaae37)*xyzzyaaay37(xyzzyaaaf37)+xyzzyaaa&
&r37(xyzzyaaae37)*xyzzyaaba37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,&
&ii)+xyzzyaaau37*xyzzyaabc37+xyzzyaaat37*xyzzyaabe37(xyzzyaaaf37)
xyzzyaaav37=xyzzyaaap37(xyzzyaaae37)*xyzzyaaaz37(xyzzyaaaf37)+xyzzyaaa&
&r37(xyzzyaaae37)*xyzzyaabb37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,ii,xyzzyaaab37)=bf_dx(xyzzyaaae37,xyzzya&
&aaf37,ii,xyzzyaaab37)+xyzzyaaav37*xyzzyaabc37
enddo
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,ii)=bf_dx(xyzzyaaae37,xyzzyaaae37,ii,&
&ii)+(xyzzyaaam37+xyzzyaaan37)*xyzzyaabc37
bf_dx(xyzzyaaae37,xyzzyaaae37,ii,xyzzyaaab37)=bf_dx(xyzzyaaae37,xyzzya&
&aae37,ii,xyzzyaaab37)-xyzzyaaam37*xyzzyaabc37
enddo
call xyzzyaagf1(xyzzyaaab37,ii,xyzzyaaag37,xyzzyaaaq37,xyzzyaaas37,xyz&
&zyaaar37,xyzzyaaam37,xyzzyaaaz37,xyzzyaaay37,xyzzyaaan37,xyzzyaabb37,&
&xyzzyaaba37)
do xyzzyaaae37=1,dimensionality
xyzzyaaat37=xyzzyaaam37*xyzzyaaaq37(xyzzyaaae37)+xyzzyaaan37*xyzzyaaas&
&37(xyzzyaaae37)
bf_x(xyzzyaaae37,xyzzyaaab37)=bf_x(xyzzyaaae37,xyzzyaaab37)+xyzzyaaat3&
&7*xyzzyaabd37
do xyzzyaaaf37=1,dimensionality
xyzzyaaau37=xyzzyaaaq37(xyzzyaaae37)*xyzzyaaaz37(xyzzyaaaf37)+xyzzyaaa&
&s37(xyzzyaaae37)*xyzzyaabb37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,xyzzyaaab37,xyzzyaaab37)=bf_dx(xyzzyaaae&
&37,xyzzyaaaf37,xyzzyaaab37,xyzzyaaab37)+xyzzyaaau37*xyzzyaabd37+xyzzy&
&aaat37*xyzzyaabf37(xyzzyaaaf37)
xyzzyaaav37=xyzzyaaaq37(xyzzyaaae37)*xyzzyaaay37(xyzzyaaaf37)+xyzzyaaa&
&s37(xyzzyaaae37)*xyzzyaaba37(xyzzyaaaf37)
bf_dx(xyzzyaaae37,xyzzyaaaf37,xyzzyaaab37,ii)=bf_dx(xyzzyaaae37,xyzzya&
&aaf37,xyzzyaaab37,ii)+xyzzyaaav37*xyzzyaabd37
enddo
bf_dx(xyzzyaaae37,xyzzyaaae37,xyzzyaaab37,xyzzyaaab37)=bf_dx(xyzzyaaae&
&37,xyzzyaaae37,xyzzyaaab37,xyzzyaaab37)+(xyzzyaaam37+xyzzyaaan37)*xyz&
&zyaabd37
bf_dx(xyzzyaaae37,xyzzyaaae37,xyzzyaaab37,ii)=bf_dx(xyzzyaaae37,xyzzya&
&aae37,xyzzyaaab37,ii)-xyzzyaaam37*xyzzyaabd37
enddo
enddo
enddo
endif
endif
endif
end subroutine xyzzyaafs1
subroutine xyzzyaaft1(eevecs,eivecs,bf_connect,bf_m,bf_rmap,bf_m2,bf_r&
&map2)
implicit none
integer,intent(out) :: bf_m(netot),bf_rmap(netot,netot)
integer,intent(out),optional :: bf_m2(netot),bf_rmap2(netot,netot)
real(dp),intent(in) :: eevecs(4,netot,netot),eivecs(4,nitot,netot)
logical,intent(out) :: bf_connect(netot,netot)
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38,xyzzyaaad38,xyzzyaaae38,xy&
&zzyaaaf38,xyzzyaaag38,xyzzyaaah38
logical xyzzyaaai38,xyzzyaaaj38(netot,netot)
bf_m=0
bf_rmap=0
bf_connect=.false.
do xyzzyaaaa38=1,netot
bf_m(xyzzyaaaa38)=bf_m(xyzzyaaaa38)+1
bf_rmap(bf_m(xyzzyaaaa38),xyzzyaaaa38)=xyzzyaaaa38
bf_connect(xyzzyaaaa38,xyzzyaaaa38)=.true.
do xyzzyaaac38=xyzzyaaaa38+1,netot
xyzzyaaai38=.false.
if(xyzzyaaaa1)then
xyzzyaaaf38=xyzzyaadm1(xyzzyaaaa38,xyzzyaaac38)
if(eevecs(4,xyzzyaaac38,xyzzyaaaa38)<xyzzyaabh1(xyzzyaaaf38))xyzzyaaai&
&38=.true.
endif
if(.not.xyzzyaaai38.and.xyzzyaaac1)then
do xyzzyaaag38=1,nitot
xyzzyaaah38=xyzzyaaat1(xyzzyaaag38)
if(xyzzyaaah38==0)cycle
if(eivecs(4,xyzzyaaag38,xyzzyaaaa38)<xyzzyaabl1(xyzzyaaah38).and.eivec&
&s(4,xyzzyaaag38,xyzzyaaac38)<xyzzyaabl1(xyzzyaaah38))then
xyzzyaaai38=.true.
exit
endif
enddo
endif
if(xyzzyaaai38)then
bf_m(xyzzyaaaa38)=bf_m(xyzzyaaaa38)+1
bf_rmap(bf_m(xyzzyaaaa38),xyzzyaaaa38)=xyzzyaaac38
bf_connect(xyzzyaaac38,xyzzyaaaa38)=.true.
bf_m(xyzzyaaac38)=bf_m(xyzzyaaac38)+1
bf_rmap(bf_m(xyzzyaaac38),xyzzyaaac38)=xyzzyaaaa38
bf_connect(xyzzyaaaa38,xyzzyaaac38)=.true.
endif
enddo
enddo
if(present(bf_m2).and.present(bf_rmap2))then
xyzzyaaaj38=.false.
do xyzzyaaaa38=1,netot
xyzzyaaaj38(xyzzyaaaa38,xyzzyaaaa38)=.true.
do xyzzyaaab38=1,bf_m(xyzzyaaaa38)
xyzzyaaac38=bf_rmap(xyzzyaaab38,xyzzyaaaa38)
do xyzzyaaad38=1,bf_m(xyzzyaaac38)
xyzzyaaae38=bf_rmap(xyzzyaaad38,xyzzyaaac38)
if(xyzzyaaae38<=xyzzyaaaa38)cycle
xyzzyaaaj38(xyzzyaaae38,xyzzyaaaa38)=.true.
xyzzyaaaj38(xyzzyaaaa38,xyzzyaaae38)=.true.
enddo
enddo
enddo
bf_m2=0
bf_rmap2=0
do xyzzyaaaa38=1,netot
bf_m2(xyzzyaaaa38)=bf_m2(xyzzyaaaa38)+1
bf_rmap2(bf_m2(xyzzyaaaa38),xyzzyaaaa38)=xyzzyaaaa38
do xyzzyaaac38=xyzzyaaaa38+1,netot
if(xyzzyaaaj38(xyzzyaaac38,xyzzyaaaa38))then
bf_m2(xyzzyaaaa38)=bf_m2(xyzzyaaaa38)+1
bf_rmap2(bf_m2(xyzzyaaaa38),xyzzyaaaa38)=xyzzyaaac38
bf_m2(xyzzyaaac38)=bf_m2(xyzzyaaac38)+1
bf_rmap2(bf_m2(xyzzyaaac38),xyzzyaaac38)=xyzzyaaaa38
endif
enddo
enddo
endif
xyzzyaabv1=xyzzyaabv1+sum(bf_m)
xyzzyaabw1=xyzzyaabw1+netot
end subroutine xyzzyaaft1
subroutine xyzzyaafu1(ii,eevecs1,eivecs,bf_m_change,bf_rmap_change,bf_&
&connect,bf_m,bf_rmap)
implicit none
integer,intent(in) :: ii
integer,intent(out) :: bf_m_change(nspin),bf_rmap_change(nemax,nspin)
integer,intent(inout) :: bf_m(netot),bf_rmap(netot,netot)
real(dp),intent(in) :: eevecs1(4,netot),eivecs(4,nitot,netot)
logical,intent(inout) :: bf_connect(netot,netot)
integer xyzzyaaaa39,xyzzyaaab39,xyzzyaaac39,xyzzyaaad39,xyzzyaaae39,xy&
&zzyaaaf39,xyzzyaaag39,xyzzyaaah39
logical xyzzyaaai39,xyzzyaaaj39
bf_m_change=0
bf_rmap_change=0
do xyzzyaaaa39=1,netot
xyzzyaaag39=which_spin(xyzzyaaaa39)
xyzzyaaah39=which_ie(xyzzyaaaa39)
if(xyzzyaaaa39==ii)then
bf_m_change(xyzzyaaag39)=bf_m_change(xyzzyaaag39)+1
bf_rmap_change(bf_m_change(xyzzyaaag39),xyzzyaaag39)=xyzzyaaah39
cycle
endif
xyzzyaaaj39=bf_connect(xyzzyaaaa39,ii)
xyzzyaaai39=.false.
if(xyzzyaaaa1)then
xyzzyaaab39=xyzzyaadm1(ii,xyzzyaaaa39)
if(eevecs1(4,xyzzyaaaa39)<xyzzyaabh1(xyzzyaaab39))xyzzyaaai39=.true.
endif
if(.not.xyzzyaaai39.and.xyzzyaaac1)then
do xyzzyaaac39=1,nitot
xyzzyaaad39=xyzzyaaat1(xyzzyaaac39)
if(xyzzyaaad39==0)cycle
if(eivecs(4,xyzzyaaac39,ii)<xyzzyaabl1(xyzzyaaad39).and.eivecs(4,xyzzy&
&aaac39,xyzzyaaaa39)<xyzzyaabl1(xyzzyaaad39))then
xyzzyaaai39=.true.
exit
endif
enddo
endif
bf_connect(xyzzyaaaa39,ii)=xyzzyaaai39
bf_connect(ii,xyzzyaaaa39)=xyzzyaaai39
if(xyzzyaaai39.or.xyzzyaaaj39)then
bf_m_change(xyzzyaaag39)=bf_m_change(xyzzyaaag39)+1
bf_rmap_change(bf_m_change(xyzzyaaag39),xyzzyaaag39)=xyzzyaaah39
else
bf_rmap_change(nele(xyzzyaaag39)+1-xyzzyaaah39+bf_m_change(xyzzyaaag39&
&),xyzzyaaag39)=xyzzyaaah39
endif
if(xyzzyaaai39.eqv.xyzzyaaaj39)cycle
if(xyzzyaaai39)then
bf_m(ii)=bf_m(ii)+1
bf_rmap(bf_m(ii),ii)=xyzzyaaaa39
bf_m(xyzzyaaaa39)=bf_m(xyzzyaaaa39)+1
bf_rmap(bf_m(xyzzyaaaa39),xyzzyaaaa39)=ii
else
xyzzyaaaf39=maxloc(bf_rmap(1:bf_m(ii),ii),1,bf_rmap(1:bf_m(ii),ii)==xy&
&zzyaaaa39)
if(xyzzyaaaf39<bf_m(ii))then
bf_rmap(xyzzyaaaf39,ii)=bf_rmap(bf_m(ii),ii)
bf_rmap(bf_m(ii),ii)=0
endif
bf_m(ii)=bf_m(ii)-1
xyzzyaaae39=maxloc(bf_rmap(1:bf_m(xyzzyaaaa39),xyzzyaaaa39),1,bf_rmap(&
&1:bf_m(xyzzyaaaa39),xyzzyaaaa39)==ii)
if(xyzzyaaae39<bf_m(xyzzyaaaa39))then
bf_rmap(xyzzyaaae39,xyzzyaaaa39)=bf_rmap(bf_m(xyzzyaaaa39),xyzzyaaaa39&
&)
bf_rmap(bf_m(xyzzyaaaa39),xyzzyaaaa39)=0
endif
bf_m(xyzzyaaaa39)=bf_m(xyzzyaaaa39)-1
endif
enddo
xyzzyaabx1=xyzzyaabx1+sum(bf_m_change)
xyzzyaaby1=xyzzyaaby1+1
end subroutine xyzzyaafu1
subroutine pbackflow_stats(in_range,in_range1)
implicit none
real(dp),intent(out) :: in_range,in_range1
in_range=0.d0
if(xyzzyaabw1>0)in_range=100.d0*dble(xyzzyaabv1)/(dble(xyzzyaabw1)*dbl&
&e(netot))
in_range1=0.d0
if(xyzzyaaby1>0)in_range1=100.d0*dble(xyzzyaabx1)/(dble(xyzzyaaby1)*db&
&le(netot))
xyzzyaabv1=0
xyzzyaabw1=0
xyzzyaabx1=0
xyzzyaaby1=0
end subroutine pbackflow_stats
subroutine xyzzyaafv1(rele,eevecs,eivecs,bf_connect,bf_m,bf_rmap,bf_x,&
&bf_dx,bf_d2x,bf_m2,bf_rmap2)
implicit none
integer,intent(out) :: bf_m(netot),bf_rmap(netot,netot)
integer,intent(out),optional :: bf_m2(netot),bf_rmap2(netot,netot)
real(dp),intent(in) :: rele(3,netot),eevecs(4,netot,netot),eivecs(4,ni&
&tot,netot)
real(dp),intent(inout) :: bf_x(3,netot)
real(dp),intent(inout),optional :: bf_dx(3,3,netot,netot),bf_d2x(3,net&
&ot,netot)
logical,intent(out) :: bf_connect(netot,netot)
integer xyzzyaaaa41,xyzzyaaab41,xyzzyaaac41,xyzzyaaad41,xyzzyaaae41,xy&
&zzyaaaf41,xyzzyaaag41,xyzzyaaah41,xyzzyaaai41,xyzzyaaaj41,xyzzyaaak41&
&(netot)
real(dp) xyzzyaaal41,xyzzyaaam41(3),xyzzyaaan41,xyzzyaaao41,xyzzyaaap4&
&1(3),xyzzyaaaq41,xyzzyaaar41,xyzzyaaas41(3),xyzzyaaat41(3),xyzzyaaau4&
&1,xyzzyaaav41,xyzzyaaaw41,xyzzyaaax41(3),xyzzyaaay41(3),xyzzyaaaz41,x&
&yzzyaaba41,xyzzyaabb41(4),xyzzyaabc41(4),xyzzyaabd41(4),xyzzyaabe41,x&
&yzzyaabf41,xyzzyaabg41,xyzzyaabh41,xyzzyaabi41,xyzzyaabj41,xyzzyaabk4&
&1,xyzzyaabl41,xyzzyaabm41,xyzzyaabn41,xyzzyaabo41(3),xyzzyaabp41(3),x&
&yzzyaabq41,xyzzyaabr41,xyzzyaabs41(netot),xyzzyaabt41(3,netot),xyzzya&
&abu41(netot)
logical xyzzyaabv41
call timer('BACKFLOW_R2X',.true.)
xyzzyaabv41=present(bf_dx).and.present(bf_d2x)
if(present(bf_m2).and.present(bf_rmap2))then
call xyzzyaaft1(eevecs,eivecs,bf_connect,bf_m,bf_rmap,bf_m2,bf_rmap2)
else
call xyzzyaaft1(eevecs,eivecs,bf_connect,bf_m,bf_rmap)
endif
call xyzzyaafw1(xyzzyaabv41,eevecs,eivecs,bf_m,bf_rmap)
if(xyzzyaadb1)then
if(xyzzyaabv41)then
call xyzzyaagk1(eivecs,xyzzyaaak41,xyzzyaabs41,xyzzyaabt41,xyzzyaabu41&
&)
else
call xyzzyaagk1(eivecs,xyzzyaaak41,xyzzyaabs41)
endif
endif
if(.not.xyzzyaabv41)then
call dcopy(three_netot,rele(1,1),1,bf_x(1,1),1)
if(.not.xyzzyaadb1)then
if(xyzzyaaaa1)then
do xyzzyaaaa41=1,netot
do xyzzyaaab41=1,bf_m(xyzzyaaaa41)
xyzzyaaac41=bf_rmap(xyzzyaaab41,xyzzyaaaa41)
if(xyzzyaaac41<=xyzzyaaaa41)cycle
xyzzyaaag41=xyzzyaadm1(xyzzyaaaa41,xyzzyaaac41)
xyzzyaabl41=eevecs(4,xyzzyaaac41,xyzzyaaaa41)
if(xyzzyaabl41>=xyzzyaabh1(xyzzyaaag41))cycle
call xyzzyaaga1(xyzzyaaaa41,xyzzyaaac41,xyzzyaabl41,xyzzyaaal41)
do xyzzyaaae41=1,dimensionality
xyzzyaabe41=xyzzyaaal41*eevecs(xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaabe4&
&1
bf_x(xyzzyaaae41,xyzzyaaac41)=bf_x(xyzzyaaae41,xyzzyaaac41)-xyzzyaabe4&
&1
enddo
enddo
enddo
endif
if(xyzzyaaab1)then
do xyzzyaaaa41=1,netot
do xyzzyaaad41=1,nitot
xyzzyaaah41=xyzzyaaan1(xyzzyaaad41)
if(xyzzyaaah41==0)cycle
if(eivecs(4,xyzzyaaad41,xyzzyaaaa41)>=xyzzyaabj1(xyzzyaaah41))cycle
call xyzzyaagd1(xyzzyaaaa41,xyzzyaaad41,xyzzyaaao41)
do xyzzyaaae41=1,dimensionality
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaaao4&
&1*eivecs(xyzzyaaae41,xyzzyaaad41,xyzzyaaaa41)
enddo
enddo
enddo
endif
if(xyzzyaaac1)then
do xyzzyaaaa41=1,netot
do xyzzyaaab41=1,bf_m(xyzzyaaaa41)
xyzzyaaac41=bf_rmap(xyzzyaaab41,xyzzyaaaa41)
if(xyzzyaaaa41==xyzzyaaac41)cycle
if(eevecs(4,xyzzyaaac41,xyzzyaaaa41)>=2*xyzzyaacp1)cycle
do xyzzyaaad41=1,nitot
xyzzyaaah41=xyzzyaaat1(xyzzyaaad41)
if(xyzzyaaah41==0)cycle
if(eivecs(4,xyzzyaaad41,xyzzyaaaa41)>=xyzzyaabl1(xyzzyaaah41).or.eivec&
&s(4,xyzzyaaad41,xyzzyaaac41)>=xyzzyaabl1(xyzzyaaah41))cycle
call xyzzyaagg1(xyzzyaaaa41,xyzzyaaac41,xyzzyaaad41,xyzzyaaar41,xyzzya&
&aaw41)
do xyzzyaaae41=1,dimensionality
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaaar4&
&1*eevecs(xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)+xyzzyaaaw41*eivecs(xyzz&
&yaaae41,xyzzyaaad41,xyzzyaaaa41)
enddo
enddo
enddo
enddo
endif
else
if(xyzzyaaaa1)then
do xyzzyaaaa41=1,netot
xyzzyaabm41=xyzzyaabs41(xyzzyaaaa41)
do xyzzyaaab41=1,bf_m(xyzzyaaaa41)
xyzzyaaac41=bf_rmap(xyzzyaaab41,xyzzyaaaa41)
if(xyzzyaaac41<=xyzzyaaaa41)cycle
xyzzyaaag41=xyzzyaadm1(xyzzyaaaa41,xyzzyaaac41)
xyzzyaabl41=eevecs(4,xyzzyaaac41,xyzzyaaaa41)
if(xyzzyaabl41>=xyzzyaabh1(xyzzyaaag41))cycle
call xyzzyaaga1(xyzzyaaaa41,xyzzyaaac41,xyzzyaabl41,xyzzyaaal41)
xyzzyaabn41=xyzzyaabs41(xyzzyaaac41)
do xyzzyaaae41=1,dimensionality
xyzzyaabe41=xyzzyaaal41*eevecs(xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaabe4&
&1*xyzzyaabm41
bf_x(xyzzyaaae41,xyzzyaaac41)=bf_x(xyzzyaaae41,xyzzyaaac41)-xyzzyaabe4&
&1*xyzzyaabn41
enddo
enddo
enddo
endif
if(xyzzyaaab1)then
do xyzzyaaaa41=1,netot
xyzzyaaai41=xyzzyaaak41(xyzzyaaaa41)
do xyzzyaaad41=1,nitot
xyzzyaaah41=xyzzyaaan1(xyzzyaaad41)
if(xyzzyaaah41==0)cycle
if(eivecs(4,xyzzyaaad41,xyzzyaaaa41)>=xyzzyaabj1(xyzzyaaah41))cycle
xyzzyaabm41=1.d0
if(xyzzyaaai41/=xyzzyaaad41)xyzzyaabm41=xyzzyaabs41(xyzzyaaaa41)
call xyzzyaagd1(xyzzyaaaa41,xyzzyaaad41,xyzzyaaao41)
do xyzzyaaae41=1,dimensionality
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaaao4&
&1*eivecs(xyzzyaaae41,xyzzyaaad41,xyzzyaaaa41)*xyzzyaabm41
enddo
enddo
enddo
endif
if(xyzzyaaac1)then
do xyzzyaaaa41=1,netot
xyzzyaaai41=xyzzyaaak41(xyzzyaaaa41)
do xyzzyaaab41=1,bf_m(xyzzyaaaa41)
xyzzyaaac41=bf_rmap(xyzzyaaab41,xyzzyaaaa41)
if(xyzzyaaaa41==xyzzyaaac41)cycle
if(eevecs(4,xyzzyaaac41,xyzzyaaaa41)>=2*xyzzyaacp1)cycle
xyzzyaaaj41=xyzzyaaak41(xyzzyaaac41)
do xyzzyaaad41=1,nitot
xyzzyaaah41=xyzzyaaat1(xyzzyaaad41)
if(xyzzyaaah41==0)cycle
if(eivecs(4,xyzzyaaad41,xyzzyaaaa41)>=xyzzyaabl1(xyzzyaaah41).or.eivec&
&s(4,xyzzyaaad41,xyzzyaaac41)>=xyzzyaabl1(xyzzyaaah41))cycle
xyzzyaabm41=1.d0
if(xyzzyaaai41/=xyzzyaaad41)xyzzyaabm41=xyzzyaabs41(xyzzyaaaa41)
xyzzyaabn41=1.d0
if(xyzzyaaaj41/=xyzzyaaad41)xyzzyaabn41=xyzzyaabs41(xyzzyaaac41)
call xyzzyaagg1(xyzzyaaaa41,xyzzyaaac41,xyzzyaaad41,xyzzyaaar41,xyzzya&
&aaw41)
do xyzzyaaae41=1,dimensionality
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+(xyzzyaaar&
&41*eevecs(xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)+xyzzyaaaw41*eivecs(xyz&
&zyaaae41,xyzzyaaad41,xyzzyaaaa41))*xyzzyaabm41
enddo
enddo
enddo
enddo
endif
endif
else
call dcopy(three_netot,rele(1,1),1,bf_x(1,1),1)
bf_dx=0.d0
bf_d2x=0.d0
do xyzzyaaaa41=1,netot
do xyzzyaaae41=1,dimensionality
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=1.d0
enddo
enddo
if(.not.xyzzyaadb1)then
if(xyzzyaaaa1)then
do xyzzyaaaa41=1,netot
do xyzzyaaab41=1,bf_m(xyzzyaaaa41)
xyzzyaaac41=bf_rmap(xyzzyaaab41,xyzzyaaaa41)
if(xyzzyaaac41<=xyzzyaaaa41)cycle
xyzzyaaag41=xyzzyaadm1(xyzzyaaaa41,xyzzyaaac41)
xyzzyaabb41=eevecs(:,xyzzyaaac41,xyzzyaaaa41)
if(xyzzyaabb41(4)>=xyzzyaabh1(xyzzyaaag41))cycle
call xyzzyaafy1(xyzzyaaaa41,xyzzyaaac41,xyzzyaabb41,xyzzyaaal41,xyzzya&
&aam41,xyzzyaaan41)
do xyzzyaaae41=1,dimensionality
xyzzyaabe41=xyzzyaaal41*xyzzyaabb41(xyzzyaaae41)
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaabe4&
&1
bf_x(xyzzyaaae41,xyzzyaaac41)=bf_x(xyzzyaaae41,xyzzyaaac41)-xyzzyaabe4&
&1
do xyzzyaaaf41=1,dimensionality
xyzzyaabf41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaam41(xyzzyaaaf41)
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaac41)-xyzzyaabf41
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaac41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaac41,xyzzyaaaa41)-xyzzyaabf41
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaac41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaac41,xyzzyaaac41)+xyzzyaabf41
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaabf41
enddo
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)-xyzzyaaal41
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)-xyzzyaaal41
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaac41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaac41,xyzzyaaac41)+xyzzyaaal41
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaaal41
xyzzyaabg41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaan41+2*xyzzyaaam41(xyzzyaa&
&ae41)
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaac41)+xyzzyaabg41
bf_d2x(xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&ac41,xyzzyaaaa41)-xyzzyaabg41
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaaa41)+xyzzyaabg41
bf_d2x(xyzzyaaae41,xyzzyaaac41,xyzzyaaac41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&ac41,xyzzyaaac41)-xyzzyaabg41
enddo
enddo
enddo
endif
if(xyzzyaaab1)then
do xyzzyaaaa41=1,netot
do xyzzyaaad41=1,nitot
xyzzyaaah41=xyzzyaaan1(xyzzyaaad41)
if(xyzzyaaah41==0)cycle
xyzzyaabc41=eivecs(:,xyzzyaaad41,xyzzyaaaa41)
if(xyzzyaabc41(4)>=xyzzyaabj1(xyzzyaaah41))cycle
call xyzzyaagb1(xyzzyaaaa41,xyzzyaaad41,xyzzyaabc41,xyzzyaaao41,xyzzya&
&aap41,xyzzyaaaq41)
do xyzzyaaae41=1,dimensionality
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaaao4&
&1*xyzzyaabc41(xyzzyaaae41)
do xyzzyaaaf41=1,dimensionality
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaaap41(xyzzyaaae41)*xyzz&
&yaabc41(xyzzyaaaf41)
enddo
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaaao41
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaaa41)+xyzzyaaaq41*xyzzyaabc41(xyzzyaaae41)+2*xyzzyaaap41(&
&xyzzyaaae41)
enddo
enddo
enddo
endif
if(xyzzyaaac1)then
do xyzzyaaaa41=1,netot
do xyzzyaaab41=1,bf_m(xyzzyaaaa41)
xyzzyaaac41=bf_rmap(xyzzyaaab41,xyzzyaaaa41)
if(xyzzyaaaa41==xyzzyaaac41)cycle
xyzzyaabb41=eevecs(:,xyzzyaaac41,xyzzyaaaa41)
if(xyzzyaabb41(4)>=2*xyzzyaacp1)cycle
do xyzzyaaad41=1,nitot
xyzzyaaah41=xyzzyaaat1(xyzzyaaad41)
if(xyzzyaaah41==0)cycle
xyzzyaabc41=eivecs(:,xyzzyaaad41,xyzzyaaaa41)
xyzzyaabd41=eivecs(:,xyzzyaaad41,xyzzyaaac41)
if(xyzzyaabc41(4)>=xyzzyaabl1(xyzzyaaah41).or.xyzzyaabd41(4)>=xyzzyaab&
&l1(xyzzyaaah41))cycle
call xyzzyaage1(xyzzyaaaa41,xyzzyaaac41,xyzzyaaad41,xyzzyaabb41,xyzzya&
&abc41,xyzzyaabd41,xyzzyaaar41,xyzzyaaas41,xyzzyaaat41,   xyzzyaaau41,&
&xyzzyaaav41,xyzzyaaaw41,xyzzyaaax41,xyzzyaaay41,xyzzyaaaz41,xyzzyaaba&
&41)
do xyzzyaaae41=1,dimensionality
xyzzyaabe41=xyzzyaaar41*xyzzyaabb41(xyzzyaaae41)+xyzzyaaaw41*xyzzyaabc&
&41(xyzzyaaae41)
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaabe4&
&1
do xyzzyaaaf41=1,dimensionality
xyzzyaabf41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaas41(xyzzyaaaf41)+xyzzyaab&
&c41(xyzzyaaae41)*xyzzyaaax41(xyzzyaaaf41)
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaabf41
xyzzyaabg41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaat41(xyzzyaaaf41)+xyzzyaab&
&c41(xyzzyaaae41)*xyzzyaaay41(xyzzyaaaf41)
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaac41)+xyzzyaabg41
enddo
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaaar41+xyzzyaaaw41
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)-xyzzyaaar41
xyzzyaabf41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaau41+xyzzyaabc41(xyzzyaaae&
&41)*xyzzyaaaz41
xyzzyaabg41=xyzzyaaas41(xyzzyaaae41)+xyzzyaaax41(xyzzyaaae41)
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaaa41)+xyzzyaabf41+2*xyzzyaabg41
xyzzyaabh41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaav41+xyzzyaabc41(xyzzyaaae&
&41)*xyzzyaaba41
xyzzyaabi41=xyzzyaaat41(xyzzyaaae41)
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaac41)+xyzzyaabh41-2*xyzzyaabi41
enddo
enddo
enddo
enddo
endif
else
if(xyzzyaaaa1)then
do xyzzyaaaa41=1,netot
xyzzyaabm41=xyzzyaabs41(xyzzyaaaa41)
xyzzyaabo41=xyzzyaabt41(:,xyzzyaaaa41)
xyzzyaabq41=xyzzyaabu41(xyzzyaaaa41)
do xyzzyaaab41=1,bf_m(xyzzyaaaa41)
xyzzyaaac41=bf_rmap(xyzzyaaab41,xyzzyaaaa41)
if(xyzzyaaac41<=xyzzyaaaa41)cycle
xyzzyaaag41=xyzzyaadm1(xyzzyaaaa41,xyzzyaaac41)
xyzzyaabb41=eevecs(:,xyzzyaaac41,xyzzyaaaa41)
if(xyzzyaabb41(4)>=xyzzyaabh1(xyzzyaaag41))cycle
call xyzzyaafy1(xyzzyaaaa41,xyzzyaaac41,xyzzyaabb41,xyzzyaaal41,xyzzya&
&aam41,xyzzyaaan41)
xyzzyaabn41=xyzzyaabs41(xyzzyaaac41)
xyzzyaabp41=xyzzyaabt41(:,xyzzyaaac41)
xyzzyaabr41=xyzzyaabu41(xyzzyaaac41)
do xyzzyaaae41=1,dimensionality
xyzzyaabe41=xyzzyaaal41*xyzzyaabb41(xyzzyaaae41)
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaabe4&
&1*xyzzyaabm41
bf_x(xyzzyaaae41,xyzzyaaac41)=bf_x(xyzzyaaae41,xyzzyaaac41)-xyzzyaabe4&
&1*xyzzyaabn41
do xyzzyaaaf41=1,dimensionality
xyzzyaabf41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaam41(xyzzyaaaf41)
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaac41)-xyzzyaabf41*xyzzyaabm41
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaac41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaac41,xyzzyaaaa41)-xyzzyaabf41*xyzzyaabn41
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaac41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaac41,xyzzyaaac41)+xyzzyaabf41*xyzzyaabn41-xyzzy&
&aabe41*xyzzyaabp41(xyzzyaaaf41)
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaabf41*xyzzyaabm41+xyzzy&
&aabe41*xyzzyaabo41(xyzzyaaaf41)
enddo
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)-xyzzyaaal41*xyzzyaabm41
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)-xyzzyaaal41*xyzzyaabn41
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaac41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaac41,xyzzyaaac41)+xyzzyaaal41*xyzzyaabn41
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaaal41*xyzzyaabm41
xyzzyaabg41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaan41+2*xyzzyaaam41(xyzzyaa&
&ae41)
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaac41)+xyzzyaabg41*xyzzyaabm41
bf_d2x(xyzzyaaae41,xyzzyaaac41,xyzzyaaaa41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&ac41,xyzzyaaaa41)-xyzzyaabg41*xyzzyaabn41
xyzzyaabh41=(xyzzyaabq41*xyzzyaaal41+2*ddot(3,xyzzyaaam41(1),1,xyzzyaa&
&bo41(1),1))*xyzzyaabb41(xyzzyaaae41)+2*xyzzyaaal41*xyzzyaabo41(xyzzya&
&aae41)
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaaa41)+xyzzyaabg41*xyzzyaabm41+xyzzyaabh41
xyzzyaabi41=(-xyzzyaabr41*xyzzyaaal41+2*ddot(3,xyzzyaaam41(1),1,xyzzya&
&abp41(1),1))*xyzzyaabb41(xyzzyaaae41)+2*xyzzyaaal41*xyzzyaabp41(xyzzy&
&aaae41)
bf_d2x(xyzzyaaae41,xyzzyaaac41,xyzzyaaac41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&ac41,xyzzyaaac41)-xyzzyaabg41*xyzzyaabn41+xyzzyaabi41
enddo
enddo
enddo
endif
if(xyzzyaaab1)then
do xyzzyaaaa41=1,netot
xyzzyaaai41=xyzzyaaak41(xyzzyaaaa41)
do xyzzyaaad41=1,nitot
xyzzyaaah41=xyzzyaaan1(xyzzyaaad41)
if(xyzzyaaah41==0)cycle
xyzzyaabc41=eivecs(:,xyzzyaaad41,xyzzyaaaa41)
if(xyzzyaabc41(4)>=xyzzyaabj1(xyzzyaaah41))cycle
call xyzzyaagb1(xyzzyaaaa41,xyzzyaaad41,xyzzyaabc41,xyzzyaaao41,xyzzya&
&aap41,xyzzyaaaq41)
xyzzyaabm41=1.d0
xyzzyaabo41=0.d0
xyzzyaabq41=0.d0
if(xyzzyaaai41/=xyzzyaaad41)then
xyzzyaabm41=xyzzyaabs41(xyzzyaaaa41)
xyzzyaabo41=xyzzyaabt41(:,xyzzyaaaa41)
xyzzyaabq41=xyzzyaabu41(xyzzyaaaa41)
endif
do xyzzyaaae41=1,dimensionality
xyzzyaabe41=xyzzyaaao41*xyzzyaabc41(xyzzyaaae41)
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaabe4&
&1*xyzzyaabm41
do xyzzyaaaf41=1,dimensionality
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaaap41(xyzzyaaae41)*xyzz&
&yaabc41(xyzzyaaaf41)*xyzzyaabm41+xyzzyaabe41*xyzzyaabo41(xyzzyaaaf41)
enddo
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaaao41*xyzzyaabm41
xyzzyaabf41=(xyzzyaabq41*xyzzyaaao41+2*ddot(3,xyzzyaabo41(1),1,xyzzyaa&
&ap41(1),1))*xyzzyaabc41(xyzzyaaae41)+2*xyzzyaabo41(xyzzyaaae41)*xyzzy&
&aaao41
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaaa41)+(xyzzyaaaq41*xyzzyaabc41(xyzzyaaae41)+2*xyzzyaaap41&
&(xyzzyaaae41))*xyzzyaabm41+xyzzyaabf41
enddo
enddo
enddo
endif
if(xyzzyaaac1)then
do xyzzyaaaa41=1,netot
xyzzyaaai41=xyzzyaaak41(xyzzyaaaa41)
do xyzzyaaab41=1,bf_m(xyzzyaaaa41)
xyzzyaaac41=bf_rmap(xyzzyaaab41,xyzzyaaaa41)
if(xyzzyaaaa41==xyzzyaaac41)cycle
xyzzyaabb41=eevecs(:,xyzzyaaac41,xyzzyaaaa41)
if(xyzzyaabb41(4)>=2*xyzzyaacp1)cycle
xyzzyaaaj41=xyzzyaaak41(xyzzyaaac41)
do xyzzyaaad41=1,nitot
xyzzyaaah41=xyzzyaaat1(xyzzyaaad41)
if(xyzzyaaah41==0)cycle
xyzzyaabc41=eivecs(:,xyzzyaaad41,xyzzyaaaa41)
xyzzyaabd41=eivecs(:,xyzzyaaad41,xyzzyaaac41)
if(xyzzyaabc41(4)>=xyzzyaabl1(xyzzyaaah41).or.xyzzyaabd41(4)>=xyzzyaab&
&l1(xyzzyaaah41))cycle
call xyzzyaage1(xyzzyaaaa41,xyzzyaaac41,xyzzyaaad41,xyzzyaabb41,xyzzya&
&abc41,xyzzyaabd41,xyzzyaaar41,xyzzyaaas41,xyzzyaaat41,   xyzzyaaau41,&
&xyzzyaaav41,xyzzyaaaw41,xyzzyaaax41,xyzzyaaay41,xyzzyaaaz41,xyzzyaaba&
&41)
xyzzyaabm41=1.d0
xyzzyaabo41=0.d0
xyzzyaabq41=0.d0
if(xyzzyaaai41/=xyzzyaaad41)then
xyzzyaabm41=xyzzyaabs41(xyzzyaaaa41)
xyzzyaabo41=xyzzyaabt41(:,xyzzyaaaa41)
xyzzyaabq41=xyzzyaabu41(xyzzyaaaa41)
endif
xyzzyaabn41=1.d0
xyzzyaabp41=0.d0
xyzzyaabr41=0.d0
if(xyzzyaaaj41/=xyzzyaaad41)then
xyzzyaabn41=xyzzyaabs41(xyzzyaaac41)
xyzzyaabp41=xyzzyaabt41(:,xyzzyaaac41)
xyzzyaabr41=xyzzyaabu41(xyzzyaaac41)
endif
do xyzzyaaae41=1,dimensionality
xyzzyaabe41=xyzzyaaar41*xyzzyaabb41(xyzzyaaae41)+xyzzyaaaw41*xyzzyaabc&
&41(xyzzyaaae41)
bf_x(xyzzyaaae41,xyzzyaaaa41)=bf_x(xyzzyaaae41,xyzzyaaaa41)+xyzzyaabe4&
&1*xyzzyaabm41
do xyzzyaaaf41=1,dimensionality
xyzzyaabf41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaas41(xyzzyaaaf41)+xyzzyaab&
&c41(xyzzyaaae41)*xyzzyaaax41(xyzzyaaaf41)
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaaa41)+xyzzyaabf41*xyzzyaabm41+xyzzy&
&aabe41*xyzzyaabo41(xyzzyaaaf41)
xyzzyaabg41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaat41(xyzzyaaaf41)+xyzzyaab&
&c41(xyzzyaaae41)*xyzzyaaay41(xyzzyaaaf41)
bf_dx(xyzzyaaae41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaaf41,xyzzyaaaa41,xyzzyaaac41)+xyzzyaabg41*xyzzyaabm41
enddo
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)+(xyzzyaaar41+xyzzyaaaw41)*xyz&
&zyaabm41
bf_dx(xyzzyaaae41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)=bf_dx(xyzzyaaae&
&41,xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)-xyzzyaaar41*xyzzyaabm41
xyzzyaabf41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaau41+xyzzyaabc41(xyzzyaaae&
&41)*xyzzyaaaz41
xyzzyaabg41=xyzzyaaas41(xyzzyaaae41)+xyzzyaaax41(xyzzyaaae41)
xyzzyaabj41=(xyzzyaabq41*xyzzyaaar41+2*ddot(3,xyzzyaabo41(1),1,xyzzyaa&
&as41(1),1))*xyzzyaabb41(xyzzyaaae41)+2*xyzzyaabo41(xyzzyaaae41)*xyzzy&
&aaar41
xyzzyaabk41=(xyzzyaabq41*xyzzyaaaw41+2*ddot(3,xyzzyaabo41(1),1,xyzzyaa&
&ax41(1),1))*xyzzyaabc41(xyzzyaaae41)+2*xyzzyaabo41(xyzzyaaae41)*xyzzy&
&aaaw41
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaaa41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaaa41)+(xyzzyaabf41+2*xyzzyaabg41)*xyzzyaabm41+xyzzyaabj41&
&+xyzzyaabk41
xyzzyaabh41=xyzzyaabb41(xyzzyaaae41)*xyzzyaaav41+xyzzyaabc41(xyzzyaaae&
&41)*xyzzyaaba41
xyzzyaabi41=xyzzyaaat41(xyzzyaaae41)
bf_d2x(xyzzyaaae41,xyzzyaaaa41,xyzzyaaac41)=bf_d2x(xyzzyaaae41,xyzzyaa&
&aa41,xyzzyaaac41)+(xyzzyaabh41-2*xyzzyaabi41)*xyzzyaabm41
enddo
enddo
enddo
enddo
endif
endif
endif
call timer('BACKFLOW_R2X',.false.)
end subroutine xyzzyaafv1
subroutine xyzzyaafw1(fsd,eevecs,eivecs,bf_m,bf_rmap)
implicit none
integer,intent(in) :: bf_m(netot),bf_rmap(netot,netot)
real(dp),intent(in) :: eevecs(4,netot,netot),eivecs(4,nitot,netot)
logical,intent(in) :: fsd
integer xyzzyaaaa42,xyzzyaaab42,xyzzyaaac42,xyzzyaaad42,xyzzyaaae42,xy&
&zzyaaaf42
real(dp) xyzzyaaag42,xyzzyaaah42
if(xyzzyaaaa1.or.xyzzyaaac1)then
if(fsd)then
do xyzzyaaaa42=1,netot
do xyzzyaaab42=1,bf_m(xyzzyaaaa42)
xyzzyaaac42=bf_rmap(xyzzyaaab42,xyzzyaaaa42)
if(xyzzyaaaa42==xyzzyaaac42)cycle
xyzzyaaah42=eevecs(4,xyzzyaaac42,xyzzyaaaa42)
if(xyzzyaaah42>xyzzyaacq1)cycle
xyzzyaaaf42=which_ee(xyzzyaaaa42,xyzzyaaac42)
call xyzzyaago1(xyzzyaaah42,xyzzyaaal1,xyzzyaacj1(1,xyzzyaaaf42),xyzzy&
&aack1(1,xyzzyaaaf42),xyzzyaacl1(1,xyzzyaaaf42))
enddo
enddo
else
do xyzzyaaaa42=1,netot
do xyzzyaaab42=1,bf_m(xyzzyaaaa42)
xyzzyaaac42=bf_rmap(xyzzyaaab42,xyzzyaaaa42)
if(xyzzyaaaa42==xyzzyaaac42)cycle
xyzzyaaah42=eevecs(4,xyzzyaaac42,xyzzyaaaa42)
if(xyzzyaaah42>xyzzyaacq1)cycle
xyzzyaaaf42=which_ee(xyzzyaaaa42,xyzzyaaac42)
call xyzzyaagq1(xyzzyaaah42,xyzzyaaal1,xyzzyaacj1(1,xyzzyaaaf42))
enddo
enddo
endif
endif
if(xyzzyaaab1.or.xyzzyaaac1)then
if(fsd)then
do xyzzyaaaa42=1,netot
do xyzzyaaad42=1,nitot
xyzzyaaae42=xyzzyaaba1(xyzzyaaad42)
if(xyzzyaaae42==0)cycle
xyzzyaaag42=eivecs(4,xyzzyaaad42,xyzzyaaaa42)
if(xyzzyaaag42>xyzzyaaci1(xyzzyaaad42))cycle
call xyzzyaago1(xyzzyaaag42,xyzzyaaae42,xyzzyaacm1(1,xyzzyaaad42,xyzzy&
&aaaa42),xyzzyaacn1(1,xyzzyaaad42,xyzzyaaaa42),xyzzyaaco1(1,xyzzyaaad4&
&2,xyzzyaaaa42))
enddo
enddo
else
do xyzzyaaaa42=1,netot
do xyzzyaaad42=1,nitot
xyzzyaaae42=xyzzyaaba1(xyzzyaaad42)
if(xyzzyaaae42==0)cycle
xyzzyaaag42=eivecs(4,xyzzyaaad42,xyzzyaaaa42)
if(xyzzyaaag42>xyzzyaaci1(xyzzyaaad42))cycle
call xyzzyaagq1(xyzzyaaag42,xyzzyaaae42,xyzzyaacm1(1,xyzzyaaad42,xyzzy&
&aaaa42))
enddo
enddo
endif
endif
end subroutine xyzzyaafw1
subroutine xyzzyaafx1(ii,fd,eevecs1,eivecs,bf_m_ii,bf_rmap_ii)
implicit none
integer,intent(in) :: ii,bf_m_ii,bf_rmap_ii(bf_m_ii)
real(dp),intent(in) :: eevecs1(4,netot),eivecs(4,nitot,netot)
logical,intent(in) :: fd
real(dp) xyzzyaaaa43,xyzzyaaab43
integer xyzzyaaac43,xyzzyaaad43,xyzzyaaae43,xyzzyaaaf43,xyzzyaaag43
if(xyzzyaaaa1.or.xyzzyaaac1)then
do xyzzyaaac43=1,bf_m_ii
xyzzyaaad43=bf_rmap_ii(xyzzyaaac43)
if(ii==xyzzyaaad43)cycle
xyzzyaaab43=eevecs1(4,xyzzyaaad43)
if(xyzzyaaab43>xyzzyaacq1)cycle
xyzzyaaag43=which_ee(ii,xyzzyaaad43)
if(fd)then
call xyzzyaagp1(xyzzyaaab43,xyzzyaaal1,xyzzyaacj1(1,xyzzyaaag43),xyzzy&
&aack1(1,xyzzyaaag43))
else
call xyzzyaagq1(xyzzyaaab43,xyzzyaaal1,xyzzyaacj1(1,xyzzyaaag43))
endif
enddo
endif
if(xyzzyaaab1.or.xyzzyaaac1)then
do xyzzyaaac43=1,bf_m_ii
xyzzyaaad43=bf_rmap_ii(xyzzyaaac43)
do xyzzyaaae43=1,nitot
xyzzyaaaf43=xyzzyaaba1(xyzzyaaae43)
if(xyzzyaaaf43==0)cycle
xyzzyaaaa43=eivecs(4,xyzzyaaae43,xyzzyaaad43)
if(xyzzyaaaa43>xyzzyaaci1(xyzzyaaae43))cycle
if(fd)then
call xyzzyaagp1(xyzzyaaaa43,xyzzyaaaf43,xyzzyaacm1(1,xyzzyaaae43,xyzzy&
&aaad43),xyzzyaacn1(1,xyzzyaaae43,xyzzyaaad43))
else
call xyzzyaagq1(xyzzyaaaa43,xyzzyaaaf43,xyzzyaacm1(1,xyzzyaaae43,xyzzy&
&aaad43))
endif
enddo
enddo
endif
end subroutine xyzzyaafx1
subroutine xyzzyaafy1(ii,jj,vecij,eta,grad_eta,lap_eta)
implicit none
integer,intent(in) :: ii,jj
real(dp),intent(in) :: vecij(:)
real(dp),intent(out) :: eta,grad_eta(3),lap_eta
integer xyzzyaaaa44,xyzzyaaab44,xyzzyaaac44
real(dp) xyzzyaaad44,xyzzyaaae44(3),xyzzyaaaf44,xyzzyaaag44,xyzzyaaah4&
&4,xyzzyaaai44,xyzzyaaaj44,xyzzyaaak44,xyzzyaaal44
eta=0.d0
xyzzyaaak44=0.d0
xyzzyaaal44=0.d0
grad_eta=0.d0
lap_eta=0.d0
xyzzyaaab44=xyzzyaadm1(ii,jj)
xyzzyaaad44=vecij(4)
xyzzyaaac44=which_ee(ii,jj)
eta=xyzzyaabi1(1,xyzzyaaab44)
do xyzzyaaaa44=2,xyzzyaaag1
eta=eta+xyzzyaabi1(xyzzyaaaa44,xyzzyaaab44)*xyzzyaacj1(xyzzyaaaa44,xyz&
&zyaaac44)
xyzzyaaak44=xyzzyaaak44+xyzzyaabi1(xyzzyaaaa44,xyzzyaaab44)*xyzzyaack1&
&(xyzzyaaaa44,xyzzyaaac44)
xyzzyaaal44=xyzzyaaal44+xyzzyaabi1(xyzzyaaaa44,xyzzyaaab44)*xyzzyaacl1&
&(xyzzyaaaa44,xyzzyaaac44)
enddo
select case(xyzzyaaad1)
case(0)
case(1)
call xyzzyaagh1(xyzzyaaad44,xyzzyaacb1(xyzzyaaab44),xyzzyaabz1(xyzzyaa&
&ab44),xyzzyaaca1(xyzzyaaab44),xyzzyaaah44,xyzzyaaai44,xyzzyaaaj44)
xyzzyaaal44=xyzzyaaal44*xyzzyaaah44+2*xyzzyaaak44*xyzzyaaai44
xyzzyaaak44=xyzzyaaak44*xyzzyaaah44+eta*xyzzyaaai44
eta=eta*xyzzyaaah44
case default
call xyzzyaagh1(xyzzyaaad44,xyzzyaacb1(xyzzyaaab44),xyzzyaabz1(xyzzyaa&
&ab44),xyzzyaaca1(xyzzyaaab44),xyzzyaaah44,xyzzyaaai44,xyzzyaaaj44)
xyzzyaaal44=xyzzyaaal44*xyzzyaaah44+2*xyzzyaaak44*xyzzyaaai44+eta*xyzz&
&yaaaj44
xyzzyaaak44=xyzzyaaak44*xyzzyaaah44+eta*xyzzyaaai44
eta=eta*xyzzyaaah44
end select
call xyzzyaagr1(xyzzyaaad44,vecij(1:3),xyzzyaaae44,xyzzyaaaf44,xyzzyaa&
&ag44)
grad_eta(:)=xyzzyaaae44(:)*xyzzyaaak44
lap_eta=xyzzyaaaf44*xyzzyaaal44+xyzzyaaag44*xyzzyaaak44
end subroutine xyzzyaafy1
subroutine xyzzyaafz1(ii,jj,vecij,eta,grad_eta)
implicit none
integer,intent(in) :: ii,jj
real(dp),intent(in) :: vecij(:)
real(dp),intent(out) :: eta,grad_eta(3)
integer xyzzyaaaa45,xyzzyaaab45,xyzzyaaac45
real(dp) xyzzyaaad45,xyzzyaaae45(3),xyzzyaaaf45,xyzzyaaag45,xyzzyaaah4&
&5
eta=0.d0
xyzzyaaah45=0.d0
grad_eta=0.d0
xyzzyaaab45=xyzzyaadm1(ii,jj)
xyzzyaaad45=vecij(4)
xyzzyaaac45=which_ee(ii,jj)
eta=xyzzyaabi1(1,xyzzyaaab45)
do xyzzyaaaa45=2,xyzzyaaag1
eta=eta+xyzzyaabi1(xyzzyaaaa45,xyzzyaaab45)*xyzzyaacj1(xyzzyaaaa45,xyz&
&zyaaac45)
xyzzyaaah45=xyzzyaaah45+xyzzyaabi1(xyzzyaaaa45,xyzzyaaab45)*xyzzyaack1&
&(xyzzyaaaa45,xyzzyaaac45)
enddo
select case(xyzzyaaad1)
case(0)
case default
call xyzzyaagi1(xyzzyaaad45,xyzzyaacb1(xyzzyaaab45),xyzzyaabz1(xyzzyaa&
&ab45),xyzzyaaaf45,xyzzyaaag45)
xyzzyaaah45=xyzzyaaah45*xyzzyaaaf45+eta*xyzzyaaag45
eta=eta*xyzzyaaaf45
end select
call xyzzyaags1(xyzzyaaad45,vecij(1:3),xyzzyaaae45)
grad_eta(:)=xyzzyaaae45(:)*xyzzyaaah45
end subroutine xyzzyaafz1
subroutine xyzzyaaga1(ii,jj,rij,eta)
implicit none
integer,intent(in) :: ii,jj
real(dp),intent(in) :: rij
real(dp),intent(out) :: eta
integer xyzzyaaaa46,xyzzyaaab46,xyzzyaaac46
real(dp) xyzzyaaad46
eta=0.d0
xyzzyaaab46=xyzzyaadm1(ii,jj)
xyzzyaaac46=which_ee(ii,jj)
eta=xyzzyaabi1(1,xyzzyaaab46)
do xyzzyaaaa46=2,xyzzyaaag1
eta=eta+xyzzyaabi1(xyzzyaaaa46,xyzzyaaab46)*xyzzyaacj1(xyzzyaaaa46,xyz&
&zyaaac46)
enddo
if(xyzzyaaad1/=0)then
call xyzzyaagj1(rij,xyzzyaacb1(xyzzyaaab46),xyzzyaaad46)
eta=eta*xyzzyaaad46
endif
end subroutine xyzzyaaga1
subroutine xyzzyaagb1(ii,ion,veci,mu,grad_mu,lap_mu)
implicit none
integer,intent(in) :: ii,ion
real(dp),intent(in) :: veci(:)
real(dp),intent(out) :: mu,grad_mu(3),lap_mu
integer xyzzyaaaa47,xyzzyaaab47,xyzzyaaac47
real(dp) xyzzyaaad47,xyzzyaaae47(3),xyzzyaaaf47,xyzzyaaag47,xyzzyaaah4&
&7,xyzzyaaai47,xyzzyaaaj47,xyzzyaaak47,xyzzyaaal47
mu=0.d0
xyzzyaaak47=0.d0
xyzzyaaal47=0.d0
grad_mu=0.d0
lap_mu=0.d0
xyzzyaaac47=xyzzyaaan1(ion)
xyzzyaaab47=xyzzyaadn1(ii,xyzzyaaan1(ion))
xyzzyaaad47=veci(4)
mu=xyzzyaabk1(1,xyzzyaaab47,xyzzyaaac47)
do xyzzyaaaa47=2,xyzzyaaap1(xyzzyaaac47)
mu=mu+xyzzyaabk1(xyzzyaaaa47,xyzzyaaab47,xyzzyaaac47)*xyzzyaacm1(xyzzy&
&aaaa47,ion,ii)
xyzzyaaak47=xyzzyaaak47+xyzzyaabk1(xyzzyaaaa47,xyzzyaaab47,xyzzyaaac47&
&)*xyzzyaacn1(xyzzyaaaa47,ion,ii)
xyzzyaaal47=xyzzyaaal47+xyzzyaabk1(xyzzyaaaa47,xyzzyaaab47,xyzzyaaac47&
&)*xyzzyaaco1(xyzzyaaaa47,ion,ii)
enddo
select case(xyzzyaaad1)
case(0)
case(1)
call xyzzyaagh1(xyzzyaaad47,xyzzyaach1(xyzzyaaac47),xyzzyaacf1(xyzzyaa&
&ac47),xyzzyaacg1(xyzzyaaac47),xyzzyaaah47,xyzzyaaai47,xyzzyaaaj47)
xyzzyaaal47=xyzzyaaal47*xyzzyaaah47+2*xyzzyaaak47*xyzzyaaai47
xyzzyaaak47=xyzzyaaak47*xyzzyaaah47+mu*xyzzyaaai47
mu=mu*xyzzyaaah47
case default
call xyzzyaagh1(xyzzyaaad47,xyzzyaach1(xyzzyaaac47),xyzzyaacf1(xyzzyaa&
&ac47),xyzzyaacg1(xyzzyaaac47),xyzzyaaah47,xyzzyaaai47,xyzzyaaaj47)
xyzzyaaal47=xyzzyaaal47*xyzzyaaah47+2*xyzzyaaak47*xyzzyaaai47+mu*xyzzy&
&aaaj47
xyzzyaaak47=xyzzyaaak47*xyzzyaaah47+mu*xyzzyaaai47
mu=mu*xyzzyaaah47
end select
call xyzzyaagr1(xyzzyaaad47,veci(1:3),xyzzyaaae47,xyzzyaaaf47,xyzzyaaa&
&g47)
grad_mu(:)=xyzzyaaae47(:)*xyzzyaaak47
lap_mu=xyzzyaaaf47*xyzzyaaal47+xyzzyaaag47*xyzzyaaak47
end subroutine xyzzyaagb1
subroutine xyzzyaagc1(ii,ion,veci,mu,grad_mu)
implicit none
integer,intent(in) :: ii,ion
real(dp),intent(in) :: veci(:)
real(dp),intent(out) :: mu,grad_mu(3)
integer xyzzyaaaa48,xyzzyaaab48,xyzzyaaac48
real(dp) xyzzyaaad48,xyzzyaaae48(3),xyzzyaaaf48,xyzzyaaag48,xyzzyaaah4&
&8
mu=0.d0
xyzzyaaah48=0.d0
grad_mu=0.d0
xyzzyaaac48=xyzzyaaan1(ion)
xyzzyaaab48=xyzzyaadn1(ii,xyzzyaaan1(ion))
xyzzyaaad48=veci(4)
mu=xyzzyaabk1(1,xyzzyaaab48,xyzzyaaac48)
do xyzzyaaaa48=2,xyzzyaaap1(xyzzyaaac48)
mu=mu+xyzzyaabk1(xyzzyaaaa48,xyzzyaaab48,xyzzyaaac48)*xyzzyaacm1(xyzzy&
&aaaa48,ion,ii)
xyzzyaaah48=xyzzyaaah48+xyzzyaabk1(xyzzyaaaa48,xyzzyaaab48,xyzzyaaac48&
&)*xyzzyaacn1(xyzzyaaaa48,ion,ii)
enddo
select case(xyzzyaaad1)
case(0)
case default
call xyzzyaagi1(xyzzyaaad48,xyzzyaach1(xyzzyaaac48),xyzzyaacf1(xyzzyaa&
&ac48),xyzzyaaaf48,xyzzyaaag48)
xyzzyaaah48=xyzzyaaah48*xyzzyaaaf48+mu*xyzzyaaag48
mu=mu*xyzzyaaaf48
end select
call xyzzyaags1(xyzzyaaad48,veci(1:3),xyzzyaaae48)
grad_mu(:)=xyzzyaaae48(:)*xyzzyaaah48
end subroutine xyzzyaagc1
subroutine xyzzyaagd1(ii,ion,mu)
implicit none
integer,intent(in) :: ii,ion
real(dp),intent(out) :: mu
integer xyzzyaaaa49,xyzzyaaab49,xyzzyaaac49
real(dp) xyzzyaaad49
mu=0.d0
xyzzyaaac49=xyzzyaaan1(ion)
xyzzyaaab49=xyzzyaadn1(ii,xyzzyaaan1(ion))
mu=xyzzyaabk1(1,xyzzyaaab49,xyzzyaaac49)
do xyzzyaaaa49=2,xyzzyaaap1(xyzzyaaac49)
mu=mu+xyzzyaabk1(xyzzyaaaa49,xyzzyaaab49,xyzzyaaac49)*xyzzyaacm1(xyzzy&
&aaaa49,ion,ii)
enddo
if(xyzzyaaad1/=0)then
call xyzzyaagj1(xyzzyaacm1(2,ion,ii),xyzzyaach1(xyzzyaaac49),xyzzyaaad&
&49)
mu=mu*xyzzyaaad49
endif
end subroutine xyzzyaagd1
subroutine xyzzyaage1(ii,jj,ion,vecij,veci,vecj,phi,grad_phi_i,grad_ph&
&i_j,lap_phi_i,lap_phi_j,theta,grad_theta_i,grad_theta_j,lap_theta_i,l&
&ap_theta_j)
implicit none
integer,intent(in) :: ii,jj,ion
real(dp),intent(in) :: vecij(4),veci(4),vecj(4)
real(dp),intent(out) :: phi,theta,grad_phi_i(3),grad_phi_j(3),grad_the&
&ta_i(3),grad_theta_j(3),lap_phi_i,lap_phi_j,lap_theta_i,lap_theta_j
integer xyzzyaaaa50,xyzzyaaab50,xyzzyaaac50,xyzzyaaad50,xyzzyaaae50,xy&
&zzyaaaf50,xyzzyaaag50
real(dp) xyzzyaaah50,xyzzyaaai50,xyzzyaaaj50,xyzzyaaak50,xyzzyaaal50,x&
&yzzyaaam50,xyzzyaaan50,xyzzyaaao50,xyzzyaaap50,xyzzyaaaq50,xyzzyaaar5&
&0,xyzzyaaas50,xyzzyaaat50,xyzzyaaau50,xyzzyaaav50,xyzzyaaaw50,xyzzyaa&
&ax50,xyzzyaaay50,xyzzyaaaz50,xyzzyaaba50,xyzzyaabb50,xyzzyaabc50,xyzz&
&yaabd50,xyzzyaabe50,xyzzyaabf50,xyzzyaabg50,xyzzyaabh50,xyzzyaabi50,x&
&yzzyaabj50,xyzzyaabk50,xyzzyaabl50,xyzzyaabm50,xyzzyaabn50,xyzzyaabo5&
&0,xyzzyaabp50,xyzzyaabq50,xyzzyaabr50,xyzzyaabs50,xyzzyaabt50,xyzzyaa&
&bu50,xyzzyaabv50,xyzzyaabw50,xyzzyaabx50,xyzzyaaby50,xyzzyaabz50,xyzz&
&yaaca50,xyzzyaacb50,xyzzyaacc50,xyzzyaacd50,xyzzyaace50,xyzzyaacf50,x&
&yzzyaacg50,xyzzyaach50,xyzzyaaci50,xyzzyaacj50,xyzzyaack50,xyzzyaacl5&
&0,xyzzyaacm50,xyzzyaacn50,xyzzyaaco50,xyzzyaacp50,xyzzyaacq50,xyzzyaa&
&cr50,xyzzyaacs50,xyzzyaact50,xyzzyaacu50,xyzzyaacv50,xyzzyaacw50,xyzz&
&yaacx50,xyzzyaacy50,xyzzyaacz50,xyzzyaada50,xyzzyaadb50(3),xyzzyaadc5&
&0(3),xyzzyaadd50(3),xyzzyaade50,xyzzyaadf50,xyzzyaadg50,xyzzyaadh50,x&
&yzzyaadi50,xyzzyaadj50,xyzzyaadk50,xyzzyaadl50,xyzzyaadm50,xyzzyaadn5&
&0,xyzzyaado50,xyzzyaadp50,xyzzyaadq50
phi=0.d0
theta=0.d0
grad_phi_i(:)=0.d0
grad_phi_j(:)=0.d0
lap_phi_i=0.d0
lap_phi_j=0.d0
grad_theta_i(:)=0.d0
grad_theta_j(:)=0.d0
lap_theta_i=0.d0
lap_theta_j=0.d0
xyzzyaaaf50=xyzzyaaat1(ion)
xyzzyaaae50=xyzzyaado1(ii,jj,xyzzyaaaf50)
xyzzyaaah50=veci(4)
xyzzyaaai50=vecj(4)
xyzzyaaaj50=vecij(4)
call xyzzyaagr1(xyzzyaaah50,veci(1:3),xyzzyaadb50,xyzzyaade50,xyzzyaad&
&j50)
call xyzzyaagr1(xyzzyaaai50,vecj(1:3),xyzzyaadc50,xyzzyaadf50,xyzzyaad&
&k50)
call xyzzyaagr1(xyzzyaaaj50,vecij(1:3),xyzzyaadd50,xyzzyaadg50,xyzzyaa&
&dl50)
xyzzyaadh50=0.d0
xyzzyaadi50=0.d0
do xyzzyaaaa50=1,dimensionality
xyzzyaadh50=xyzzyaadh50+xyzzyaadb50(xyzzyaaaa50)*xyzzyaadd50(xyzzyaaaa&
&50)
xyzzyaadi50=xyzzyaadi50-xyzzyaadc50(xyzzyaaaa50)*xyzzyaadd50(xyzzyaaaa&
&50)
enddo
xyzzyaaag50=which_ee(ii,jj)
xyzzyaabe50=0.d0
xyzzyaabf50=0.d0
xyzzyaabg50=0.d0
xyzzyaabh50=0.d0
xyzzyaabi50=0.d0
xyzzyaabj50=0.d0
xyzzyaabk50=0.d0
xyzzyaabl50=0.d0
xyzzyaabm50=0.d0
xyzzyaabn50=0.d0
xyzzyaabo50=0.d0
xyzzyaabp50=0.d0
xyzzyaabq50=0.d0
xyzzyaabr50=0.d0
xyzzyaabs50=0.d0
xyzzyaabt50=0.d0
xyzzyaabu50=0.d0
xyzzyaabv50=0.d0
do xyzzyaaad50=1,xyzzyaaav1(xyzzyaaaf50)
xyzzyaaat50=xyzzyaacj1(xyzzyaaad50,xyzzyaaag50)
xyzzyaaau50=xyzzyaack1(xyzzyaaad50,xyzzyaaag50)
xyzzyaaav50=xyzzyaacl1(xyzzyaaad50,xyzzyaaag50)
do xyzzyaaac50=1,xyzzyaaaw1(xyzzyaaaf50)
xyzzyaaaq50=xyzzyaacm1(xyzzyaaac50,ion,jj)
xyzzyaaar50=xyzzyaacn1(xyzzyaaac50,ion,jj)
xyzzyaaas50=xyzzyaaco1(xyzzyaaac50,ion,jj)
xyzzyaaaw50=xyzzyaaaq50*xyzzyaaat50
xyzzyaaax50=xyzzyaaaq50*xyzzyaaau50
xyzzyaaay50=xyzzyaaar50*xyzzyaaat50
xyzzyaaaz50=xyzzyaaaq50*xyzzyaaav50
xyzzyaaba50=xyzzyaaar50*xyzzyaaau50
xyzzyaabb50=xyzzyaaas50*xyzzyaaat50
do xyzzyaaab50=1,xyzzyaaaw1(xyzzyaaaf50)
xyzzyaaan50=xyzzyaacm1(xyzzyaaab50,ion,ii)
xyzzyaaao50=xyzzyaacn1(xyzzyaaab50,ion,ii)
xyzzyaaap50=xyzzyaaco1(xyzzyaaab50,ion,ii)
xyzzyaabc50=xyzzyaabm1(xyzzyaaab50,xyzzyaaac50,xyzzyaaad50,xyzzyaaae50&
&,xyzzyaaaf50)
xyzzyaabd50=xyzzyaabn1(xyzzyaaab50,xyzzyaaac50,xyzzyaaad50,xyzzyaaae50&
&,xyzzyaaaf50)
xyzzyaabw50=xyzzyaaan50*xyzzyaaaw50
xyzzyaabx50=xyzzyaaao50*xyzzyaaaw50
xyzzyaaby50=xyzzyaaan50*xyzzyaaay50
xyzzyaabz50=xyzzyaaan50*xyzzyaaax50
xyzzyaaca50=xyzzyaaap50*xyzzyaaaw50
xyzzyaacb50=xyzzyaaan50*xyzzyaabb50
xyzzyaacc50=xyzzyaaan50*xyzzyaaaz50
xyzzyaacd50=xyzzyaaao50*xyzzyaaax50
xyzzyaace50=xyzzyaaan50*xyzzyaaba50
xyzzyaabe50=xyzzyaabe50+xyzzyaabc50*xyzzyaabw50
xyzzyaabf50=xyzzyaabf50+xyzzyaabc50*xyzzyaabx50
xyzzyaabg50=xyzzyaabg50+xyzzyaabc50*xyzzyaaby50
xyzzyaabh50=xyzzyaabh50+xyzzyaabc50*xyzzyaabz50
xyzzyaabi50=xyzzyaabi50+xyzzyaabc50*xyzzyaaca50
xyzzyaabj50=xyzzyaabj50+xyzzyaabc50*xyzzyaacb50
xyzzyaabk50=xyzzyaabk50+xyzzyaabc50*xyzzyaacc50
xyzzyaabl50=xyzzyaabl50+xyzzyaabc50*xyzzyaacd50
xyzzyaabm50=xyzzyaabm50+xyzzyaabc50*xyzzyaace50
xyzzyaabn50=xyzzyaabn50+xyzzyaabd50*xyzzyaabw50
xyzzyaabo50=xyzzyaabo50+xyzzyaabd50*xyzzyaabx50
xyzzyaabp50=xyzzyaabp50+xyzzyaabd50*xyzzyaaby50
xyzzyaabq50=xyzzyaabq50+xyzzyaabd50*xyzzyaabz50
xyzzyaabr50=xyzzyaabr50+xyzzyaabd50*xyzzyaaca50
xyzzyaabs50=xyzzyaabs50+xyzzyaabd50*xyzzyaacb50
xyzzyaabt50=xyzzyaabt50+xyzzyaabd50*xyzzyaacc50
xyzzyaabu50=xyzzyaabu50+xyzzyaabd50*xyzzyaacd50
xyzzyaabv50=xyzzyaabv50+xyzzyaabd50*xyzzyaace50
enddo
enddo
enddo
call xyzzyaagh1(xyzzyaaah50,xyzzyaace1(xyzzyaaaf50),xyzzyaacc1(xyzzyaa&
&af50),xyzzyaacd1(xyzzyaaaf50),xyzzyaaak50,xyzzyaacf50,xyzzyaach50)
call xyzzyaagh1(xyzzyaaai50,xyzzyaace1(xyzzyaaaf50),xyzzyaacc1(xyzzyaa&
&af50),xyzzyaacd1(xyzzyaaaf50),xyzzyaaal50,xyzzyaacg50,xyzzyaaci50)
xyzzyaaam50=xyzzyaaak50*xyzzyaaal50
phi=xyzzyaaam50*xyzzyaabe50
theta=xyzzyaaam50*xyzzyaabn50
xyzzyaadq50=xyzzyaaam50*xyzzyaadg50
xyzzyaacr50=xyzzyaadq50*xyzzyaabk50
xyzzyaada50=xyzzyaadq50*xyzzyaabt50
xyzzyaacj50=xyzzyaacf50*xyzzyaabe50+xyzzyaaak50*xyzzyaabf50
xyzzyaacs50=xyzzyaacf50*xyzzyaabn50+xyzzyaaak50*xyzzyaabo50
xyzzyaack50=xyzzyaacg50*xyzzyaabe50+xyzzyaaal50*xyzzyaabg50
xyzzyaact50=xyzzyaacg50*xyzzyaabn50+xyzzyaaal50*xyzzyaabp50
xyzzyaacl50=xyzzyaaak50*xyzzyaabh50
xyzzyaacu50=xyzzyaaak50*xyzzyaabq50
xyzzyaacm50=xyzzyaaal50*xyzzyaabh50
xyzzyaacv50=xyzzyaaal50*xyzzyaabq50
do xyzzyaaaa50=1,dimensionality
xyzzyaadm50=xyzzyaadb50(xyzzyaaaa50)*xyzzyaacj50+xyzzyaadd50(xyzzyaaaa&
&50)*xyzzyaacl50
xyzzyaadn50=xyzzyaadc50(xyzzyaaaa50)*xyzzyaack50-xyzzyaadd50(xyzzyaaaa&
&50)*xyzzyaacm50
grad_phi_i(xyzzyaaaa50)=xyzzyaadm50*xyzzyaaal50
grad_phi_j(xyzzyaaaa50)=xyzzyaadn50*xyzzyaaak50
xyzzyaadm50=xyzzyaadb50(xyzzyaaaa50)*xyzzyaacs50+xyzzyaadd50(xyzzyaaaa&
&50)*xyzzyaacu50
xyzzyaadn50=xyzzyaadc50(xyzzyaaaa50)*xyzzyaact50-xyzzyaadd50(xyzzyaaaa&
&50)*xyzzyaacv50
grad_theta_i(xyzzyaaaa50)=xyzzyaadm50*xyzzyaaal50
grad_theta_j(xyzzyaaaa50)=xyzzyaadn50*xyzzyaaak50
enddo
xyzzyaacn50=xyzzyaach50*xyzzyaabe50+2*xyzzyaacf50*xyzzyaabf50+xyzzyaaa&
&k50*xyzzyaabi50
xyzzyaacw50=xyzzyaach50*xyzzyaabn50+2*xyzzyaacf50*xyzzyaabo50+xyzzyaaa&
&k50*xyzzyaabr50
xyzzyaaco50=xyzzyaaci50*xyzzyaabe50+2*xyzzyaacg50*xyzzyaabg50+xyzzyaaa&
&l50*xyzzyaabj50
xyzzyaacx50=xyzzyaaci50*xyzzyaabn50+2*xyzzyaacg50*xyzzyaabp50+xyzzyaaa&
&l50*xyzzyaabs50
xyzzyaacp50=xyzzyaacf50*xyzzyaabh50+xyzzyaaak50*xyzzyaabl50
xyzzyaacy50=xyzzyaacf50*xyzzyaabq50+xyzzyaaak50*xyzzyaabu50
xyzzyaacq50=xyzzyaacg50*xyzzyaabh50+xyzzyaaal50*xyzzyaabm50
xyzzyaacz50=xyzzyaacg50*xyzzyaabq50+xyzzyaaal50*xyzzyaabv50
xyzzyaado50=xyzzyaade50*xyzzyaacn50+xyzzyaadj50*xyzzyaacj50+2*xyzzyaad&
&h50*xyzzyaacp50+xyzzyaadl50*xyzzyaacl50
xyzzyaadp50=xyzzyaadf50*xyzzyaaco50+xyzzyaadk50*xyzzyaack50+2*xyzzyaad&
&i50*xyzzyaacq50+xyzzyaadl50*xyzzyaacm50
lap_phi_i=xyzzyaado50*xyzzyaaal50+xyzzyaacr50
lap_phi_j=xyzzyaadp50*xyzzyaaak50+xyzzyaacr50
xyzzyaado50=xyzzyaade50*xyzzyaacw50+xyzzyaadj50*xyzzyaacs50+2*xyzzyaad&
&h50*xyzzyaacy50+xyzzyaadl50*xyzzyaacu50
xyzzyaadp50=xyzzyaadf50*xyzzyaacx50+xyzzyaadk50*xyzzyaact50+2*xyzzyaad&
&i50*xyzzyaacz50+xyzzyaadl50*xyzzyaacv50
lap_theta_i=xyzzyaado50*xyzzyaaal50+xyzzyaada50
lap_theta_j=xyzzyaadp50*xyzzyaaak50+xyzzyaada50
end subroutine xyzzyaage1
subroutine xyzzyaagf1(ii,jj,ion,vecij,veci,vecj,phi,grad_phi_i,grad_ph&
&i_j,theta,grad_theta_i,grad_theta_j)
implicit none
integer,intent(in) :: ii,jj,ion
real(dp),intent(in) :: vecij(4),veci(4),vecj(4)
real(dp),intent(out) :: phi,theta,grad_phi_i(3),grad_phi_j(3),grad_the&
&ta_i(3),grad_theta_j(3)
integer xyzzyaaaa51,xyzzyaaab51,xyzzyaaac51,xyzzyaaad51,xyzzyaaae51,xy&
&zzyaaaf51,xyzzyaaag51
real(dp) xyzzyaaah51,xyzzyaaai51,xyzzyaaaj51,xyzzyaaak51,xyzzyaaal51,x&
&yzzyaaam51,xyzzyaaan51,xyzzyaaao51,xyzzyaaap51,xyzzyaaaq51,xyzzyaaar5&
&1,xyzzyaaas51,xyzzyaaat51,xyzzyaaau51,xyzzyaaav51,xyzzyaaaw51,xyzzyaa&
&ax51,xyzzyaaay51,xyzzyaaaz51,xyzzyaaba51,xyzzyaabb51,xyzzyaabc51,xyzz&
&yaabd51,xyzzyaabe51,xyzzyaabf51,xyzzyaabg51,xyzzyaabh51,xyzzyaabi51,x&
&yzzyaabj51,xyzzyaabk51,xyzzyaabl51,xyzzyaabm51,xyzzyaabn51,xyzzyaabo5&
&1,xyzzyaabp51,xyzzyaabq51,xyzzyaabr51,xyzzyaabs51,xyzzyaabt51,xyzzyaa&
&bu51(3),xyzzyaabv51(3),xyzzyaabw51(3),xyzzyaabx51,xyzzyaaby51
phi=0.d0
theta=0.d0
grad_phi_i(:)=0.d0
grad_phi_j(:)=0.d0
grad_theta_i(:)=0.d0
grad_theta_j(:)=0.d0
xyzzyaaaf51=xyzzyaaat1(ion)
xyzzyaaae51=xyzzyaado1(ii,jj,xyzzyaaaf51)
xyzzyaaah51=veci(4)
xyzzyaaai51=vecj(4)
xyzzyaaaj51=vecij(4)
call xyzzyaags1(xyzzyaaah51,veci(1:3),xyzzyaabu51)
call xyzzyaags1(xyzzyaaai51,vecj(1:3),xyzzyaabv51)
call xyzzyaags1(xyzzyaaaj51,vecij(1:3),xyzzyaabw51)
xyzzyaaag51=which_ee(ii,jj)
xyzzyaaay51=0.d0
xyzzyaaaz51=0.d0
xyzzyaaba51=0.d0
xyzzyaabb51=0.d0
xyzzyaabc51=0.d0
xyzzyaabd51=0.d0
xyzzyaabe51=0.d0
xyzzyaabf51=0.d0
do xyzzyaaad51=1,xyzzyaaav1(xyzzyaaaf51)
xyzzyaaar51=xyzzyaacj1(xyzzyaaad51,xyzzyaaag51)
xyzzyaaas51=xyzzyaack1(xyzzyaaad51,xyzzyaaag51)
do xyzzyaaac51=1,xyzzyaaaw1(xyzzyaaaf51)
xyzzyaaap51=xyzzyaacm1(xyzzyaaac51,ion,jj)
xyzzyaaaq51=xyzzyaacn1(xyzzyaaac51,ion,jj)
xyzzyaaat51=xyzzyaaap51*xyzzyaaar51
xyzzyaaau51=xyzzyaaap51*xyzzyaaas51
xyzzyaaav51=xyzzyaaaq51*xyzzyaaar51
do xyzzyaaab51=1,xyzzyaaaw1(xyzzyaaaf51)
xyzzyaaan51=xyzzyaacm1(xyzzyaaab51,ion,ii)
xyzzyaaao51=xyzzyaacn1(xyzzyaaab51,ion,ii)
xyzzyaaaw51=xyzzyaabm1(xyzzyaaab51,xyzzyaaac51,xyzzyaaad51,xyzzyaaae51&
&,xyzzyaaaf51)
xyzzyaaax51=xyzzyaabn1(xyzzyaaab51,xyzzyaaac51,xyzzyaaad51,xyzzyaaae51&
&,xyzzyaaaf51)
xyzzyaabg51=xyzzyaaan51*xyzzyaaat51
xyzzyaabh51=xyzzyaaao51*xyzzyaaat51
xyzzyaabi51=xyzzyaaan51*xyzzyaaav51
xyzzyaabj51=xyzzyaaan51*xyzzyaaau51
xyzzyaaay51=xyzzyaaay51+xyzzyaaaw51*xyzzyaabg51
xyzzyaaaz51=xyzzyaaaz51+xyzzyaaaw51*xyzzyaabh51
xyzzyaaba51=xyzzyaaba51+xyzzyaaaw51*xyzzyaabi51
xyzzyaabb51=xyzzyaabb51+xyzzyaaaw51*xyzzyaabj51
xyzzyaabc51=xyzzyaabc51+xyzzyaaax51*xyzzyaabg51
xyzzyaabd51=xyzzyaabd51+xyzzyaaax51*xyzzyaabh51
xyzzyaabe51=xyzzyaabe51+xyzzyaaax51*xyzzyaabi51
xyzzyaabf51=xyzzyaabf51+xyzzyaaax51*xyzzyaabj51
enddo
enddo
enddo
call xyzzyaagi1(xyzzyaaah51,xyzzyaace1(xyzzyaaaf51),xyzzyaacc1(xyzzyaa&
&af51),xyzzyaaak51,xyzzyaabk51)
call xyzzyaagi1(xyzzyaaai51,xyzzyaace1(xyzzyaaaf51),xyzzyaacc1(xyzzyaa&
&af51),xyzzyaaal51,xyzzyaabl51)
xyzzyaaam51=xyzzyaaak51*xyzzyaaal51
phi=xyzzyaaam51*xyzzyaaay51
theta=xyzzyaaam51*xyzzyaabc51
xyzzyaabm51=xyzzyaabk51*xyzzyaaay51+xyzzyaaak51*xyzzyaaaz51
xyzzyaabq51=xyzzyaabk51*xyzzyaabc51+xyzzyaaak51*xyzzyaabd51
xyzzyaabn51=xyzzyaabl51*xyzzyaaay51+xyzzyaaal51*xyzzyaaba51
xyzzyaabr51=xyzzyaabl51*xyzzyaabc51+xyzzyaaal51*xyzzyaabe51
xyzzyaabo51=xyzzyaaak51*xyzzyaabb51
xyzzyaabs51=xyzzyaaak51*xyzzyaabf51
xyzzyaabp51=xyzzyaaal51*xyzzyaabb51
xyzzyaabt51=xyzzyaaal51*xyzzyaabf51
do xyzzyaaaa51=1,dimensionality
xyzzyaabx51=xyzzyaabu51(xyzzyaaaa51)*xyzzyaabm51+xyzzyaabw51(xyzzyaaaa&
&51)*xyzzyaabo51
xyzzyaaby51=xyzzyaabv51(xyzzyaaaa51)*xyzzyaabn51-xyzzyaabw51(xyzzyaaaa&
&51)*xyzzyaabp51
grad_phi_i(xyzzyaaaa51)=xyzzyaabx51*xyzzyaaal51
grad_phi_j(xyzzyaaaa51)=xyzzyaaby51*xyzzyaaak51
xyzzyaabx51=xyzzyaabu51(xyzzyaaaa51)*xyzzyaabq51+xyzzyaabw51(xyzzyaaaa&
&51)*xyzzyaabs51
xyzzyaaby51=xyzzyaabv51(xyzzyaaaa51)*xyzzyaabr51-xyzzyaabw51(xyzzyaaaa&
&51)*xyzzyaabt51
grad_theta_i(xyzzyaaaa51)=xyzzyaabx51*xyzzyaaal51
grad_theta_j(xyzzyaaaa51)=xyzzyaaby51*xyzzyaaak51
enddo
end subroutine xyzzyaagf1
subroutine xyzzyaagg1(ii,jj,ion,phi,theta)
implicit none
integer,intent(in) :: ii,jj,ion
real(dp),intent(out) :: phi,theta
integer xyzzyaaaa52,xyzzyaaab52,xyzzyaaac52,xyzzyaaad52,xyzzyaaae52,xy&
&zzyaaaf52
real(dp) xyzzyaaag52,xyzzyaaah52,xyzzyaaai52,xyzzyaaaj52,xyzzyaaak52,x&
&yzzyaaal52,xyzzyaaam52,xyzzyaaan52
phi=0.d0
theta=0.d0
xyzzyaaae52=xyzzyaaat1(ion)
xyzzyaaad52=xyzzyaado1(ii,jj,xyzzyaaae52)
xyzzyaaaf52=which_ee(ii,jj)
xyzzyaaaj52=0.d0
xyzzyaaak52=0.d0
do xyzzyaaac52=1,xyzzyaaav1(xyzzyaaae52)
xyzzyaaam52=xyzzyaacj1(xyzzyaaac52,xyzzyaaaf52)
do xyzzyaaab52=1,xyzzyaaaw1(xyzzyaaae52)
xyzzyaaan52=xyzzyaacm1(xyzzyaaab52,ion,jj)*xyzzyaaam52
do xyzzyaaaa52=1,xyzzyaaaw1(xyzzyaaae52)
xyzzyaaal52=xyzzyaacm1(xyzzyaaaa52,ion,ii)*xyzzyaaan52
xyzzyaaaj52=xyzzyaaaj52+xyzzyaabm1(xyzzyaaaa52,xyzzyaaab52,xyzzyaaac52&
&,xyzzyaaad52,xyzzyaaae52)*xyzzyaaal52
xyzzyaaak52=xyzzyaaak52+xyzzyaabn1(xyzzyaaaa52,xyzzyaaab52,xyzzyaaac52&
&,xyzzyaaad52,xyzzyaaae52)*xyzzyaaal52
enddo
enddo
enddo
call xyzzyaagj1(xyzzyaacm1(2,ion,ii),xyzzyaace1(xyzzyaaae52),xyzzyaaag&
&52)
call xyzzyaagj1(xyzzyaacm1(2,ion,jj),xyzzyaace1(xyzzyaaae52),xyzzyaaah&
&52)
xyzzyaaai52=xyzzyaaag52*xyzzyaaah52
phi=xyzzyaaai52*xyzzyaaaj52
theta=xyzzyaaai52*xyzzyaaak52
end subroutine xyzzyaagg1
subroutine xyzzyaagh1(r,linv,dfactor,d2factor,truncf,d_truncf,d2_trunc&
&f)
implicit none
real(dp),intent(in) :: r,linv,dfactor,d2factor
real(dp),intent(out) :: truncf,d_truncf,d2_truncf
real(dp) xyzzyaaaa53,xyzzyaaab53,xyzzyaaac53,xyzzyaaad53
if(xyzzyaaad1==0)then
truncf=1.d0
d_truncf=0.d0
d2_truncf=0.d0
return
endif
xyzzyaaaa53=1.d0-r*linv
select case(xyzzyaaad1)
case(1)
truncf=xyzzyaaaa53
d_truncf=dfactor
d2_truncf=0.d0
case(2)
truncf=xyzzyaaaa53*xyzzyaaaa53
d_truncf=dfactor*xyzzyaaaa53
d2_truncf=d2factor
case(3)
xyzzyaaab53=xyzzyaaaa53*xyzzyaaaa53
truncf=xyzzyaaab53*xyzzyaaaa53
d_truncf=dfactor*xyzzyaaab53
d2_truncf=d2factor*xyzzyaaaa53
case default
xyzzyaaac53=xyzzyaaaa53**xyzzyaaaf1
xyzzyaaad53=xyzzyaaac53*xyzzyaaaa53
truncf=xyzzyaaad53*xyzzyaaaa53
d_truncf=dfactor*xyzzyaaad53
d2_truncf=d2factor*xyzzyaaac53
end select
end subroutine xyzzyaagh1
subroutine xyzzyaagi1(r,linv,dfactor,truncf,d_truncf)
implicit none
real(dp),intent(in) :: r,linv,dfactor
real(dp),intent(out) :: truncf,d_truncf
real(dp) xyzzyaaaa54,xyzzyaaab54,xyzzyaaac54
if(xyzzyaaad1==0)then
truncf=1.d0
d_truncf=0.d0
return
endif
xyzzyaaaa54=1.d0-r*linv
select case(xyzzyaaad1)
case(1)
truncf=xyzzyaaaa54
d_truncf=dfactor
case(2)
truncf=xyzzyaaaa54*xyzzyaaaa54
d_truncf=dfactor*xyzzyaaaa54
case(3)
xyzzyaaab54=xyzzyaaaa54*xyzzyaaaa54
truncf=xyzzyaaab54*xyzzyaaaa54
d_truncf=dfactor*xyzzyaaab54
case default
xyzzyaaac54=xyzzyaaaa54**xyzzyaaae1
truncf=xyzzyaaac54*xyzzyaaaa54
d_truncf=dfactor*xyzzyaaac54
end select
end subroutine xyzzyaagi1
subroutine xyzzyaagj1(r,linv,truncf)
implicit none
real(dp),intent(in) :: r,linv
real(dp),intent(out) :: truncf
real(dp) xyzzyaaaa55,xyzzyaaab55
if(xyzzyaaad1==0)then
truncf=1.d0
return
endif
xyzzyaaaa55=1.d0-r*linv
select case(xyzzyaaad1)
case(1)
truncf=xyzzyaaaa55
case(2)
truncf=xyzzyaaaa55*xyzzyaaaa55
case(3)
xyzzyaaab55=xyzzyaaaa55*xyzzyaaaa55
truncf=xyzzyaaab55*xyzzyaaaa55
case default
truncf=xyzzyaaaa55**xyzzyaaad1
end select
end subroutine xyzzyaagj1
subroutine xyzzyaagk1(eivecs,ae_ion_close,f_cut_ae,grad_cut_ae,lap_cut&
&_ae)
implicit none
integer,intent(out) :: ae_ion_close(netot)
real(dp),intent(in) :: eivecs(4,nitot,netot)
real(dp),intent(out) :: f_cut_ae(netot)
real(dp),intent(out),optional :: grad_cut_ae(3,netot),lap_cut_ae(netot&
&)
integer xyzzyaaaa56,xyzzyaaab56
real(dp) xyzzyaaac56(4),xyzzyaaad56,xyzzyaaae56,xyzzyaaaf56(3)
logical xyzzyaaag56,xyzzyaaah56,xyzzyaaai56
xyzzyaaag56=present(grad_cut_ae)
xyzzyaaah56=present(lap_cut_ae)
xyzzyaaai56=xyzzyaaag56.or.xyzzyaaah56
f_cut_ae=1.d0
if(xyzzyaaag56)grad_cut_ae=0.d0
if(xyzzyaaah56)lap_cut_ae=0.d0
do xyzzyaaaa56=1,netot
xyzzyaaab56=xyzzyaagn1(eivecs(1,1,xyzzyaaaa56))
ae_ion_close(xyzzyaaaa56)=xyzzyaaab56
if(xyzzyaaab56==0)cycle
if(xyzzyaaai56)then
xyzzyaaac56(1:4)=eivecs(1:4,xyzzyaaab56,xyzzyaaaa56)
call xyzzyaagl1(xyzzyaaab56,xyzzyaaac56,f_cut_ae(xyzzyaaaa56),xyzzyaaa&
&f56,xyzzyaaae56)
if(xyzzyaaag56)grad_cut_ae(1:3,xyzzyaaaa56)=xyzzyaaaf56
if(xyzzyaaah56)lap_cut_ae(xyzzyaaaa56)=xyzzyaaae56
else
xyzzyaaad56=eivecs(4,xyzzyaaab56,xyzzyaaaa56)
call xyzzyaagm1(xyzzyaaab56,xyzzyaaad56,f_cut_ae(xyzzyaaaa56))
endif
enddo
end subroutine xyzzyaagk1
subroutine xyzzyaagl1(ion,veci,truncf,grad_truncf,lap_truncf)
implicit none
integer,intent(in) :: ion
real(dp),intent(in) :: veci(4)
real(dp),intent(out) :: truncf,grad_truncf(3),lap_truncf
real(dp) xyzzyaaaa57,xyzzyaaab57,xyzzyaaac57,xyzzyaaad57,xyzzyaaae57,x&
&yzzyaaaf57,xyzzyaaag57,xyzzyaaah57,xyzzyaaai57,xyzzyaaaj57(3),xyzzyaa&
&ak57,xyzzyaaal57
xyzzyaaaa57=veci(4)
if(xyzzyaaaa57>xyzzyaadf1(ion))then
truncf=1.d0
grad_truncf=0.d0
lap_truncf=0.d0
return
endif
call xyzzyaagr1(xyzzyaaaa57,veci(1:3),xyzzyaaaj57,xyzzyaaak57,xyzzyaaa&
&l57)
xyzzyaaab57=xyzzyaaaa57*xyzzyaadg1(ion)
xyzzyaaac57=xyzzyaaab57*xyzzyaaab57
xyzzyaaad57=xyzzyaaab57+xyzzyaaab57
xyzzyaaae57=xyzzyaaad57+xyzzyaaad57
xyzzyaaaf57=xyzzyaaae57+xyzzyaaae57
xyzzyaaag57=xyzzyaaac57+xyzzyaaac57+xyzzyaaac57
truncf=xyzzyaaac57*(6.d0-xyzzyaaaf57+xyzzyaaag57)
xyzzyaaah57=xyzzyaadh1(ion)*xyzzyaaab57*(1.d0-xyzzyaaad57+xyzzyaaac57)
xyzzyaaai57=xyzzyaadi1(ion)*(1.d0-xyzzyaaae57+xyzzyaaag57)
grad_truncf(1:3)=xyzzyaaah57*xyzzyaaaj57(1:3)
lap_truncf=xyzzyaaai57*xyzzyaaak57+xyzzyaaah57*xyzzyaaal57
end subroutine xyzzyaagl1
subroutine xyzzyaagm1(ion,r,truncf)
implicit none
integer,intent(in) :: ion
real(dp),intent(in) :: r
real(dp),intent(out) :: truncf
real(dp) xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58,xyzzyaaae58,x&
&yzzyaaaf58
if(r>xyzzyaadf1(ion))then
truncf=1.d0
return
endif
xyzzyaaaa58=r*xyzzyaadg1(ion)
xyzzyaaab58=xyzzyaaaa58*xyzzyaaaa58
xyzzyaaac58=xyzzyaaaa58+xyzzyaaaa58
xyzzyaaad58=xyzzyaaac58+xyzzyaaac58
xyzzyaaae58=xyzzyaaad58+xyzzyaaad58
xyzzyaaaf58=xyzzyaaab58+xyzzyaaab58+xyzzyaaab58
truncf=xyzzyaaab58*(6.d0-xyzzyaaae58+xyzzyaaaf58)
end subroutine xyzzyaagm1
integer function xyzzyaagn1(eivecs1)
implicit none
real(dp),intent(in) :: eivecs1(4,nitot)
integer xyzzyaaaa59
do xyzzyaaaa59=1,nitot
if(.not.is_ae(xyzzyaaaa59))cycle
if(eivecs1(4,xyzzyaaaa59)<xyzzyaadf1(xyzzyaaaa59))then
xyzzyaagn1=xyzzyaaaa59
return
endif
enddo
xyzzyaagn1=0
end function xyzzyaagn1
subroutine xyzzyaago1(r,n,pows,dpows,d2pows)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r
real(dp),intent(inout) :: pows(*),dpows(*),d2pows(*)
integer xyzzyaaaa60,xyzzyaaab60
select case(n)
case(0)
case(1)
pows(1)=1.d0
dpows(1)=0.d0
d2pows(1)=0.d0
case(2)
pows(1)=1.d0
dpows(1)=0.d0
d2pows(1)=0.d0
pows(2)=r
dpows(2)=1.d0
d2pows(2)=0.d0
case default
pows(1)=1.d0
dpows(1)=0.d0
d2pows(1)=0.d0
pows(2)=r
dpows(2)=1.d0
d2pows(2)=0.d0
xyzzyaaab60=2
do xyzzyaaaa60=3,n
d2pows(xyzzyaaaa60)=xyzzyaaab60*dpows(xyzzyaaab60)
dpows(xyzzyaaaa60)=xyzzyaaab60*pows(xyzzyaaab60)
pows(xyzzyaaaa60)=pows(xyzzyaaab60)*r
xyzzyaaab60=xyzzyaaaa60
enddo
end select
end subroutine xyzzyaago1
subroutine xyzzyaagp1(r,n,pows,dpows)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r
real(dp),intent(inout) :: pows(*),dpows(*)
integer xyzzyaaaa61,xyzzyaaab61
select case(n)
case(0)
case(1)
pows(1)=1.d0
dpows(1)=0.d0
case default
pows(1)=1.d0
dpows(1)=0.d0
pows(2)=r
dpows(2)=1.d0
xyzzyaaab61=2
do xyzzyaaaa61=3,n
dpows(xyzzyaaaa61)=xyzzyaaab61*pows(xyzzyaaab61)
pows(xyzzyaaaa61)=pows(xyzzyaaab61)*r
xyzzyaaab61=xyzzyaaaa61
enddo
end select
end subroutine xyzzyaagp1
subroutine xyzzyaagq1(r,n,pows)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r
real(dp),intent(inout) :: pows(*)
integer xyzzyaaaa62
real(dp) xyzzyaaab62
pows(1)=1.d0
xyzzyaaab62=1.d0
do xyzzyaaaa62=2,n
xyzzyaaab62=xyzzyaaab62*r
pows(xyzzyaaaa62)=xyzzyaaab62
enddo
end subroutine xyzzyaagq1
subroutine xyzzyaagr1(r,vecr,grad,grad2,lap)
implicit none
real(dp),intent(in) :: r,vecr(:)
real(dp),intent(out) :: grad(3),grad2,lap
real(dp) xyzzyaaaa63
grad2=1.d0
if(r==0.d0)then
grad(1:3)=0.d0
lap=0.d0
return
endif
select case(dimensionality)
case(3)
xyzzyaaaa63=1.d0/r
lap=2*xyzzyaaaa63
grad(1:3)=vecr(1:3)*xyzzyaaaa63
case(2)
xyzzyaaaa63=1.d0/r
lap=xyzzyaaaa63
grad(1:2)=vecr(1:2)*xyzzyaaaa63
grad(3)=0.d0
case(1)
grad(1)=sign(1.d0,vecr(1))
grad(2:3)=0.d0
lap=0.d0
end select
end subroutine xyzzyaagr1
subroutine xyzzyaags1(r,vecr,grad)
implicit none
real(dp),intent(in) :: r,vecr(:)
real(dp),intent(out) :: grad(3)
real(dp) xyzzyaaaa64
if(r==0.d0)then
grad(1:3)=0.d0
return
endif
select case(dimensionality)
case(3)
xyzzyaaaa64=1.d0/r
grad(1:3)=vecr(1:3)*xyzzyaaaa64
case(2)
xyzzyaaaa64=1.d0/r
grad(1:2)=vecr(1:2)*xyzzyaaaa64
grad(3)=0.d0
case(1)
grad(1)=sign(1.d0,vecr(1))
grad(2:3)=0.d0
end select
end subroutine xyzzyaags1
subroutine xyzzyaagt1(farray,bf_dx_ii,bf_m_ii,bf_rmap_ii,loggrad_psi)
implicit none
integer,intent(in) :: bf_m_ii,bf_rmap_ii(bf_m_ii)
real(dp),intent(in) :: farray(3,real1_complex2,netot),bf_dx_ii(3,3,net&
&ot)
complex(dp),intent(out) :: loggrad_psi(3)
integer xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65
real(dp) xyzzyaaad65(3),xyzzyaaae65(3)
if(.not.complex_wf)then
xyzzyaaad65=0.d0
do xyzzyaaaa65=1,bf_m_ii
xyzzyaaab65=bf_rmap_ii(xyzzyaaaa65)
do xyzzyaaac65=1,dimensionality
xyzzyaaad65(xyzzyaaac65)=xyzzyaaad65(xyzzyaaac65)+ddot(dimensionality,&
&farray(1,1,xyzzyaaab65),1,bf_dx_ii(1,xyzzyaaac65,xyzzyaaab65),1)
enddo
enddo
loggrad_psi=cmplx(xyzzyaaad65,0.d0,dp)
else
xyzzyaaad65=0.d0
xyzzyaaae65=0.d0
do xyzzyaaaa65=1,bf_m_ii
xyzzyaaab65=bf_rmap_ii(xyzzyaaaa65)
do xyzzyaaac65=1,dimensionality
xyzzyaaad65(xyzzyaaac65)=xyzzyaaad65(xyzzyaaac65)+ddot(dimensionality,&
&farray(1,1,xyzzyaaab65),1,bf_dx_ii(1,xyzzyaaac65,xyzzyaaab65),1)
xyzzyaaae65(xyzzyaaac65)=xyzzyaaae65(xyzzyaaac65)+ddot(dimensionality,&
&farray(1,2,xyzzyaaab65),1,bf_dx_ii(1,xyzzyaaac65,xyzzyaaab65),1)
enddo
enddo
loggrad_psi=cmplx(xyzzyaaad65,xyzzyaaae65,dp)
endif
end subroutine xyzzyaagt1
subroutine xyzzyaagu1(farray,harray,bf_dx_ii,bf_d2x_ii,bf_m_ii,bf_rmap&
&_ii,loggrad_psi,loglap_psi)
implicit none
integer,intent(in) :: bf_m_ii,bf_rmap_ii(bf_m_ii)
real(dp),intent(in) :: farray(3,real1_complex2,netot),harray(3,3,real1&
&_complex2,netot,netot),bf_dx_ii(3,3,netot),bf_d2x_ii(3,netot)
complex(dp),intent(in) :: loggrad_psi(3)
complex(dp),intent(out) :: loglap_psi
integer xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66,xyzzyaaad66,xyzzyaaae66,xy&
&zzyaaaf66
real(dp) xyzzyaaag66,xyzzyaaah66,xyzzyaaai66,xyzzyaaaj66
if(.not.complex_wf)then
xyzzyaaag66=0.d0
do xyzzyaaae66=1,dimensionality
xyzzyaaai66=dble(loggrad_psi(xyzzyaaae66))
xyzzyaaag66=xyzzyaaag66-xyzzyaaai66*xyzzyaaai66
enddo
do xyzzyaaaa66=1,bf_m_ii
xyzzyaaab66=bf_rmap_ii(xyzzyaaaa66)
xyzzyaaag66=xyzzyaaag66+ddot(dimensionality,farray(1,1,xyzzyaaab66),1,&
&bf_d2x_ii(1,xyzzyaaab66),1)
do xyzzyaaae66=1,dimensionality
xyzzyaaai66=ddot(dimensionality,bf_dx_ii(xyzzyaaae66,1,xyzzyaaab66),3,&
&bf_dx_ii(xyzzyaaae66,1,xyzzyaaab66),3)
xyzzyaaag66=xyzzyaaag66+harray(xyzzyaaae66,xyzzyaaae66,1,xyzzyaaab66,x&
&yzzyaaab66)*xyzzyaaai66
do xyzzyaaaf66=xyzzyaaae66+1,dimensionality
xyzzyaaai66=ddot(dimensionality,bf_dx_ii(xyzzyaaae66,1,xyzzyaaab66),3,&
&bf_dx_ii(xyzzyaaaf66,1,xyzzyaaab66),3)
xyzzyaaai66=xyzzyaaai66+xyzzyaaai66
xyzzyaaag66=xyzzyaaag66+harray(xyzzyaaae66,xyzzyaaaf66,1,xyzzyaaab66,x&
&yzzyaaab66)*xyzzyaaai66
enddo
enddo
do xyzzyaaac66=xyzzyaaaa66+1,bf_m_ii
xyzzyaaad66=bf_rmap_ii(xyzzyaaac66)
do xyzzyaaae66=1,dimensionality
do xyzzyaaaf66=1,dimensionality
xyzzyaaai66=ddot(dimensionality,bf_dx_ii(xyzzyaaae66,1,xyzzyaaab66),3,&
&bf_dx_ii(xyzzyaaaf66,1,xyzzyaaad66),3)
xyzzyaaag66=xyzzyaaag66+(harray(xyzzyaaae66,xyzzyaaaf66,1,xyzzyaaab66,&
&xyzzyaaad66)+harray(xyzzyaaaf66,xyzzyaaae66,1,xyzzyaaad66,xyzzyaaab66&
&))*xyzzyaaai66
enddo
enddo
enddo
enddo
loglap_psi=cmplx(xyzzyaaag66,0.d0,dp)
else
xyzzyaaag66=0.d0
xyzzyaaah66=0.d0
do xyzzyaaae66=1,dimensionality
xyzzyaaai66=dble(loggrad_psi(xyzzyaaae66))
xyzzyaaaj66=aimag(loggrad_psi(xyzzyaaae66))
xyzzyaaag66=xyzzyaaag66-xyzzyaaai66*xyzzyaaai66-xyzzyaaaj66*xyzzyaaaj6&
&6
enddo
do xyzzyaaaa66=1,bf_m_ii
xyzzyaaab66=bf_rmap_ii(xyzzyaaaa66)
xyzzyaaag66=xyzzyaaag66+ddot(dimensionality,farray(1,1,xyzzyaaab66),1,&
&bf_d2x_ii(1,xyzzyaaab66),1)
xyzzyaaah66=xyzzyaaah66+ddot(dimensionality,farray(1,2,xyzzyaaab66),1,&
&bf_d2x_ii(1,xyzzyaaab66),1)
do xyzzyaaae66=1,dimensionality
xyzzyaaai66=ddot(dimensionality,bf_dx_ii(xyzzyaaae66,1,xyzzyaaab66),3,&
&bf_dx_ii(xyzzyaaae66,1,xyzzyaaab66),3)
xyzzyaaag66=xyzzyaaag66+harray(xyzzyaaae66,xyzzyaaae66,1,xyzzyaaab66,x&
&yzzyaaab66)*xyzzyaaai66
xyzzyaaah66=xyzzyaaah66+harray(xyzzyaaae66,xyzzyaaae66,2,xyzzyaaab66,x&
&yzzyaaab66)*xyzzyaaai66
do xyzzyaaaf66=xyzzyaaae66+1,dimensionality
xyzzyaaai66=ddot(dimensionality,bf_dx_ii(xyzzyaaae66,1,xyzzyaaab66),3,&
&bf_dx_ii(xyzzyaaaf66,1,xyzzyaaab66),3)
xyzzyaaai66=xyzzyaaai66+xyzzyaaai66
xyzzyaaag66=xyzzyaaag66+harray(xyzzyaaae66,xyzzyaaaf66,1,xyzzyaaab66,x&
&yzzyaaab66)*xyzzyaaai66
xyzzyaaah66=xyzzyaaah66+harray(xyzzyaaae66,xyzzyaaaf66,2,xyzzyaaab66,x&
&yzzyaaab66)*xyzzyaaai66
enddo
enddo
do xyzzyaaac66=xyzzyaaaa66+1,bf_m_ii
xyzzyaaad66=bf_rmap_ii(xyzzyaaac66)
do xyzzyaaae66=1,dimensionality
do xyzzyaaaf66=1,dimensionality
xyzzyaaai66=ddot(dimensionality,bf_dx_ii(xyzzyaaae66,1,xyzzyaaab66),3,&
&bf_dx_ii(xyzzyaaaf66,1,xyzzyaaad66),3)
xyzzyaaag66=xyzzyaaag66+(harray(xyzzyaaae66,xyzzyaaaf66,1,xyzzyaaab66,&
&xyzzyaaad66)+harray(xyzzyaaaf66,xyzzyaaae66,1,xyzzyaaad66,xyzzyaaab66&
&))*xyzzyaaai66
xyzzyaaah66=xyzzyaaah66+(harray(xyzzyaaae66,xyzzyaaaf66,2,xyzzyaaab66,&
&xyzzyaaad66)+harray(xyzzyaaaf66,xyzzyaaae66,2,xyzzyaaad66,xyzzyaaab66&
&))*xyzzyaaai66
enddo
enddo
enddo
enddo
loglap_psi=cmplx(xyzzyaaag66,xyzzyaaah66,dp)
endif
end subroutine xyzzyaagu1
subroutine setup_pbackflow_plot(makeplot,phiplot,ii,z,rii)
implicit none
integer,intent(in),optional :: ii
real(dp),intent(in),optional :: z,rii
logical,intent(in) :: makeplot
logical,intent(in),optional :: phiplot
xyzzyaads1=makeplot
xyzzyaadt1=.false.
if(present(phiplot))xyzzyaadt1=phiplot
xyzzyaadu1=0
if(present(ii))xyzzyaadu1=ii
xyzzyaadq1=0.d0
if(present(z))xyzzyaadq1=z
xyzzyaadr1=0.d0
if(present(rii))xyzzyaadr1=rii
end subroutine setup_pbackflow_plot
subroutine plot_pbackflow(rele)
use slaarnaan, only : ee_distances_all,ee_distances,ei_distances_all,e&
&i_distances,map_to_simcell
implicit none
real(dp),intent(inout) :: rele(3,netot)
integer xyzzyaaaa68,xyzzyaaab68,xyzzyaaac68,xyzzyaaad68,xyzzyaaae68,xy&
&zzyaaaf68,xyzzyaaag68,xyzzyaaah68,xyzzyaaai68,xyzzyaaaj68,xyzzyaaak68&
&,xyzzyaaal68,xyzzyaaam68(2),xyzzyaaan68,xyzzyaaao68,xyzzyaaap68,xyzzy&
&aaaq68,xyzzyaaar68,xyzzyaaas68,xyzzyaaat68,xyzzyaaau68
real(dp) xyzzyaaav68(3),xyzzyaaaw68,xyzzyaaax68,xyzzyaaay68,xyzzyaaaz6&
&8,xyzzyaaba68,xyzzyaabb68,xyzzyaabc68,xyzzyaabd68,xyzzyaabe68,xyzzyaa&
&bf68,xyzzyaabg68,xyzzyaabh68,xyzzyaabi68,xyzzyaabj68,xyzzyaabk68,xyzz&
&yaabl68,xyzzyaabm68,xyzzyaabn68,xyzzyaabo68,xyzzyaabp68,xyzzyaabq68
integer,allocatable :: xyzzyaabr68(:),xyzzyaabs68(:,:)
logical,allocatable :: xyzzyaabt68(:,:)
real(dp),allocatable :: xyzzyaabu68(:,:,:),xyzzyaabv68(:,:,:),xyzzyaab&
&w68(:,:,:),xyzzyaabx68(:,:,:),xyzzyaaby68(:,:),xyzzyaabz68(:,:),xyzzy&
&aaca68(:,:)
if(.not.(am_master.and.xyzzyaads1))return
if(dimensionality==1)call errstop('PLOT_BACKFLOW','Backflow plots are &
&not implemented in 1D.')
allocate(xyzzyaabw68(4,netot,netot),xyzzyaabx68(4,nitot,netot),xyzzyaa&
&by68(3,netot),xyzzyaabz68(3,netot),xyzzyaabt68(netot,netot),xyzzyaabr&
&68(netot),xyzzyaabs68(netot,netot),stat=xyzzyaaau68)
call check_alloc(xyzzyaaau68,'PLOT_BACKFLOW','0')
xyzzyaaat68=xyzzyaadu1
xyzzyaabl68=xyzzyaadq1
if(xyzzyaaat68<1.or.xyzzyaaat68>netot)call errstop('PLOT_BACKFLOW','Ba&
&d particle number.')
xyzzyaaal68=0
xyzzyaaan68=200
xyzzyaaao68=40
xyzzyaaaq68=40
if(isperiodic)then
xyzzyaaaw68=wigner_seitz_radius
else
xyzzyaaaw68=0.d0
if(xyzzyaaab1.or.xyzzyaaac1)then
do xyzzyaaai68=1,nitot
xyzzyaabp68=0.d0
xyzzyaabq68=0.d0
xyzzyaabo68=sqrt(rion(1,xyzzyaaai68)**2+rion(2,xyzzyaaai68)**2+rion(3,&
&xyzzyaaai68)**2)
if(xyzzyaaab1)then
xyzzyaaah68=xyzzyaaan1(xyzzyaaai68)
if(xyzzyaaah68/=0)xyzzyaabp68=xyzzyaabo68+xyzzyaabj1(xyzzyaaah68)
endif
if(xyzzyaaac1)then
xyzzyaaah68=xyzzyaaat1(xyzzyaaai68)
if(xyzzyaaah68/=0)xyzzyaabq68=xyzzyaabo68+xyzzyaabl1(xyzzyaaah68)
endif
xyzzyaaaw68=max(xyzzyaabq68,xyzzyaabp68,xyzzyaaaw68)
enddo
endif
if(xyzzyaaaa1)then
do xyzzyaaac68=1,netot
if(xyzzyaaac68==xyzzyaaat68)cycle
xyzzyaabo68=sqrt(rele(1,xyzzyaaac68)**2+rele(2,xyzzyaaac68)**2+rele(3,&
&xyzzyaaac68)**2)
xyzzyaaaw68=max(xyzzyaabo68+maxval(xyzzyaabh1),xyzzyaaaw68)
enddo
endif
endif
if(xyzzyaaaw68<=0.d0)then
deallocate(xyzzyaabw68,xyzzyaabx68,xyzzyaaby68,xyzzyaabz68,xyzzyaabt68&
&,xyzzyaabr68,xyzzyaabs68)
return
endif
if(xyzzyaaac1.or.xyzzyaaab1)then
open(unit=xyzzyaadp1,file='bfions.dat',status='replace',iostat=xyzzyaa&
&ak68)
do xyzzyaaai68=1,nitot
if(xyzzyaaac1.and.xyzzyaaab1)then
if(xyzzyaaat1(xyzzyaaai68)==0.and.xyzzyaaan1(xyzzyaaai68)==0)cycle
elseif(xyzzyaaab1)then
if(xyzzyaaan1(xyzzyaaai68)==0)cycle
else
if(xyzzyaaat1(xyzzyaaai68)==0)cycle
endif
write(xyzzyaadp1,'(i3,1x,e20.12,1x,e20.12,1x,e20.12)')xyzzyaaai68,rion&
&(1:3,xyzzyaaai68)
enddo
close(xyzzyaadp1)
endif
open(unit=xyzzyaadp1,file='bfconfig.dat',status='replace',iostat=xyzzy&
&aaak68)
if(xyzzyaaak68/=0)call errstop('PLOT_BACKFLOW','Problem opening bfconf&
&ig.dat.')
xyzzyaaby68=rele
if(isperiodic)call map_to_simcell(3,netot,xyzzyaaby68)
do xyzzyaaac68=1,netot
if(xyzzyaaac68==xyzzyaaat68)then
write(xyzzyaadp1,'(i3,1x,e20.12,1x,e20.12,1x,e20.12)')0,xyzzyaaby68(:,&
&xyzzyaaac68)
else
write(xyzzyaadp1,'(i3,1x,e20.12,1x,e20.12,1x,e20.12)')xyzzyaaac68,xyzz&
&yaaby68(:,xyzzyaaac68)
endif
enddo
close(xyzzyaadp1)
open(unit=xyzzyaadp1,file='bfconfigx.dat',status='replace',iostat=xyzz&
&yaaak68)
if(xyzzyaaak68/=0)call errstop('PLOT_BACKFLOW','Problem opening bfconf&
&igx.dat .')
call ee_distances_all(netot,rele,xyzzyaabw68)
if(nitot>0)call ei_distances_all(netot,rele,nitot,rion,xyzzyaabx68)
call xyzzyaafv1(rele,xyzzyaabw68,xyzzyaabx68,xyzzyaabt68,xyzzyaabr68,x&
&yzzyaabs68,xyzzyaabz68)
xyzzyaaby68=xyzzyaabz68
if(isperiodic)call map_to_simcell(3,netot,xyzzyaaby68)
do xyzzyaaac68=1,netot
if(xyzzyaaac68==xyzzyaaat68)then
write(xyzzyaadp1,'(i3,1x,e20.12,1x,e20.12,1x,e20.12)')0,xyzzyaaby68(:,&
&xyzzyaaac68)
else
write(xyzzyaadp1,'(i3,1x,e20.12,1x,e20.12,1x,e20.12)')xyzzyaaac68,xyzz&
&yaaby68(:,xyzzyaaac68)
endif
enddo
close(xyzzyaadp1)
deallocate(xyzzyaaby68)
xyzzyaaav68(:)=rele(:,xyzzyaaat68)
if(dimensionality==3)rele(3,xyzzyaaat68)=xyzzyaabl68
xyzzyaabk68=xyzzyaaaw68/real(xyzzyaaaq68,dp)
open(unit=xyzzyaadp1,file='bffield.dat',status='replace',iostat=xyzzya&
&aak68)
allocate(xyzzyaabu68(-xyzzyaaaq68:xyzzyaaaq68,-xyzzyaaaq68:xyzzyaaaq68&
&,5))
allocate(xyzzyaaca68(4,netot))
xyzzyaabu68=0.d0
do xyzzyaaac68=-xyzzyaaaq68,xyzzyaaaq68
do xyzzyaaad68=-xyzzyaaaq68,xyzzyaaaq68
rele(1,xyzzyaaat68)=xyzzyaaac68*xyzzyaabk68
rele(2,xyzzyaaat68)=xyzzyaaad68*xyzzyaabk68
call ee_distances(netot,rele(1:3,xyzzyaaat68),rele,xyzzyaaca68)
call dcopy_ee(netot,xyzzyaaat68,xyzzyaaca68,xyzzyaabw68)
if(nitot>0)call ei_distances(rele(1,xyzzyaaat68),nitot,rion,xyzzyaabx6&
&8(1,1,xyzzyaaat68))
call xyzzyaafv1(rele,xyzzyaabw68,xyzzyaabx68,xyzzyaabt68,xyzzyaabr68,x&
&yzzyaabs68,xyzzyaabz68)
xyzzyaabu68(xyzzyaaac68,xyzzyaaad68,:)=(/rele(1,xyzzyaaat68),rele(2,xy&
&zzyaaat68),xyzzyaabz68(1,xyzzyaaat68)-rele(1,xyzzyaaat68),xyzzyaabz68&
&(2,xyzzyaaat68)-rele(2,xyzzyaaat68),xyzzyaabz68(3,xyzzyaaat68)-rele(3&
&,xyzzyaaat68)/)
enddo
enddo
deallocate(xyzzyaaca68)
xyzzyaabn68=maxval(xyzzyaabu68(:,:,3)**2+xyzzyaabu68(:,:,4)**2+xyzzyaa&
&bu68(:,:,5)**2)
if(xyzzyaabn68==0.d0)then
xyzzyaabn68=1.d0
else
xyzzyaabn68=2.d0*xyzzyaabk68/sqrt(xyzzyaabn68)
endif
do xyzzyaaac68=-xyzzyaaaq68,xyzzyaaaq68
do xyzzyaaad68=-xyzzyaaaq68,xyzzyaaaq68
if(any(xyzzyaabu68(xyzzyaaac68,xyzzyaaad68,3:5)/=0.d0))then
if(dimensionality==3)then
write(xyzzyaadp1,'(e15.7,4(1x,e15.7))')xyzzyaabu68(xyzzyaaac68,xyzzyaa&
&ad68,1:2),xyzzyaabu68(xyzzyaaac68,xyzzyaaad68,3:5)*xyzzyaabn68
elseif(dimensionality==2)then
write(xyzzyaadp1,'(e18.10,3(1x,e18.10))')xyzzyaabu68(xyzzyaaac68,xyzzy&
&aaad68,1:2),xyzzyaabu68(xyzzyaaac68,xyzzyaaad68,3:4)*xyzzyaabn68
else
call errstop('PLOT_BACKFLOW','1D backflow plots not coded.')
endif
endif
enddo
enddo
deallocate(xyzzyaabu68)
close(xyzzyaadp1)
rele(:,xyzzyaaat68)=xyzzyaaav68(:)
if(xyzzyaaaa1)then
do xyzzyaaag68=1,xyzzyaaak1
if(.not.any(xyzzyaadm1==xyzzyaaag68))cycle
xyzzyaaam68=minloc(xyzzyaadm1,xyzzyaadm1==xyzzyaaag68)
xyzzyaaae68=xyzzyaaam68(1)
xyzzyaaaf68=xyzzyaaam68(2)
if(xyzzyaaae68==xyzzyaaaf68)then
xyzzyaaaf68=xyzzyaaaf68+1
if(xyzzyaadm1(xyzzyaaae68,xyzzyaaaf68)/=xyzzyaaag68)cycle
endif
xyzzyaaal68=which_ee(xyzzyaaae68,xyzzyaaaf68)
open(unit=xyzzyaadp1,file='bfeta_'//trim(i2s(xyzzyaaag68))//'.dat',sta&
&tus='replace',iostat=xyzzyaaak68)
if(xyzzyaaak68/=0)call errstop('PLOT_BACKFLOW','Problem opening bfeta_&
&s.dat .')
do xyzzyaaac68=0,xyzzyaaan68
xyzzyaabb68=real(xyzzyaaac68,dp)*xyzzyaabh1(xyzzyaaag68)/real(xyzzyaaa&
&n68,dp)
call xyzzyaagq1(xyzzyaabb68,xyzzyaaag1,xyzzyaacj1(1,xyzzyaaal68))
call xyzzyaaga1(xyzzyaaae68,xyzzyaaaf68,xyzzyaabb68,xyzzyaaax68)
write(xyzzyaadp1,'(e20.12,1x,e20.12)')xyzzyaabb68,xyzzyaaax68
enddo
close(xyzzyaadp1)
enddo
endif
if(xyzzyaaab1)then
do xyzzyaaah68=1,xyzzyaaah1
if(xyzzyaaao1(xyzzyaaah68)==0)cycle
xyzzyaaai68=xyzzyaaar1(1,xyzzyaaah68)
do xyzzyaaag68=1,xyzzyaaas1(xyzzyaaah68)
if(.not.any(xyzzyaadn1(:,xyzzyaaah68)==xyzzyaaag68))cycle
xyzzyaaae68=minloc(xyzzyaadn1(:,xyzzyaaah68),1,xyzzyaadn1(:,xyzzyaaah6&
&8)==xyzzyaaag68)
open(unit=xyzzyaadp1,file='bfmu_'//trim(i2s(xyzzyaaag68))//'_'//trim(i&
&2s(xyzzyaaah68))//'.dat',status='replace',iostat=xyzzyaaak68)
if(xyzzyaaak68/=0)call errstop('PLOT_BACKFLOW','Problem opening bfmu_s&
&_set.dat')
do xyzzyaaac68=0,xyzzyaaan68
xyzzyaabc68=real(xyzzyaaac68,dp)*xyzzyaabj1(xyzzyaaah68)/real(xyzzyaaa&
&n68,dp)
call xyzzyaagq1(xyzzyaabc68,xyzzyaaap1(xyzzyaaah68),xyzzyaacm1(1,xyzzy&
&aaai68,xyzzyaaae68))
call xyzzyaagd1(xyzzyaaae68,xyzzyaaai68,xyzzyaaay68)
write(xyzzyaadp1,'(e20.12,1x,e20.12)')xyzzyaabc68,xyzzyaaay68
enddo
close(xyzzyaadp1)
enddo
enddo
endif
if(xyzzyaaac1.and.xyzzyaadt1)then
do xyzzyaaah68=1,xyzzyaaai1
if(xyzzyaaau1(xyzzyaaah68)==0)cycle
xyzzyaaai68=xyzzyaaax1(1,xyzzyaaah68)
do xyzzyaaag68=1,xyzzyaaaz1(xyzzyaaah68)
if(.not.any(xyzzyaado1==xyzzyaaag68))cycle
xyzzyaaam68=minloc(xyzzyaado1(:,:,xyzzyaaah68),xyzzyaado1(:,:,xyzzyaaa&
&h68)==xyzzyaaag68)
xyzzyaaae68=xyzzyaaam68(1)
xyzzyaaaf68=xyzzyaaam68(2)
if(xyzzyaaae68==xyzzyaaaf68)then
xyzzyaaaf68=xyzzyaaaf68+1
if(xyzzyaado1(xyzzyaaae68,xyzzyaaaf68,xyzzyaaah68)/=xyzzyaaag68)cycle
endif
xyzzyaaal68=which_ee(xyzzyaaae68,xyzzyaaaf68)
xyzzyaabd68=xyzzyaadr1
if(xyzzyaabd68>xyzzyaabl1(xyzzyaaah68))cycle
open(unit=xyzzyaadp1,file='bfphi_'//trim(i2s(xyzzyaaag68))//'_'//trim(&
&i2s(xyzzyaaah68))//'.dat',status='replace',iostat=xyzzyaaak68)
if(xyzzyaaak68/=0)call errstop('PLOT_BACKFLOW','Problem opening bfphi_&
&s_set.dat')
xyzzyaaaa68=xyzzyaaaw1(xyzzyaaah68)
xyzzyaaab68=xyzzyaaav1(xyzzyaaah68)
call xyzzyaagq1(xyzzyaabd68,xyzzyaaaa68,xyzzyaacm1(1,xyzzyaaai68,xyzzy&
&aaae68))
xyzzyaaap68=2**nint(log(xyzzyaaao68*pi)/log(2.d0))
allocate(xyzzyaabv68(0:xyzzyaaap68,0:xyzzyaaao68,5))
xyzzyaabv68=0.d0
do xyzzyaaac68=0,xyzzyaaap68
xyzzyaabh68=xyzzyaaac68*pi/real(xyzzyaaap68,dp)
xyzzyaaaj68=0
if(xyzzyaaac68/=0)then
xyzzyaaar68=xyzzyaaac68
xyzzyaaas68=nint(log(xyzzyaaao68*pi)/log(2.d0))
do
if(mod(xyzzyaaar68,2)==0)then
xyzzyaaar68=xyzzyaaar68/2
xyzzyaaas68=xyzzyaaas68-1
else
exit
endif
enddo
xyzzyaaaj68=nint((2**xyzzyaaas68)/pi)
endif
do xyzzyaaad68=xyzzyaaaj68,xyzzyaaao68
xyzzyaabe68=xyzzyaaad68*xyzzyaabl1(xyzzyaaah68)/real(xyzzyaaao68,dp)
xyzzyaabf68=xyzzyaabe68*cos(xyzzyaabh68)
xyzzyaabg68=xyzzyaabe68*sin(xyzzyaabh68)
xyzzyaabb68=sqrt(xyzzyaabd68**2+xyzzyaabe68**2-2*xyzzyaabd68*xyzzyaabf&
&68)
call xyzzyaagq1(xyzzyaabe68,xyzzyaaaa68,xyzzyaacm1(1,xyzzyaaai68,xyzzy&
&aaaf68))
call xyzzyaagq1(xyzzyaabb68,xyzzyaaab68,xyzzyaacj1(1,xyzzyaaal68))
call xyzzyaagg1(xyzzyaaaf68,xyzzyaaae68,xyzzyaaai68,xyzzyaaaz68,xyzzya&
&aba68)
xyzzyaabi68=xyzzyaaaz68*(xyzzyaabf68-xyzzyaabd68)+xyzzyaaba68*xyzzyaab&
&f68
xyzzyaabj68=(xyzzyaaaz68+xyzzyaaba68)*xyzzyaabg68
xyzzyaabv68(xyzzyaaac68,xyzzyaaad68,:)=(/xyzzyaabf68,xyzzyaabg68,xyzzy&
&aabi68,xyzzyaabj68,1.d0/)
enddo
enddo
xyzzyaabm68=maxval(xyzzyaabv68(:,:,3)**2+xyzzyaabv68(:,:,4)**2)
if(xyzzyaabm68==0.d0)then
xyzzyaabm68=1.d0
else
xyzzyaabm68=1.5d0*xyzzyaabl1(xyzzyaaah68)/(sqrt(xyzzyaabm68)*xyzzyaaao&
&68)
endif
do xyzzyaaac68=0,xyzzyaaap68
do xyzzyaaad68=0,xyzzyaaao68
if((xyzzyaaac68>0.and.xyzzyaaad68==0).or.xyzzyaabv68(xyzzyaaac68,xyzzy&
&aaad68,5)==0.d0)cycle
xyzzyaabf68=xyzzyaabv68(xyzzyaaac68,xyzzyaaad68,1)
xyzzyaabg68=xyzzyaabv68(xyzzyaaac68,xyzzyaaad68,2)
xyzzyaabi68=xyzzyaabv68(xyzzyaaac68,xyzzyaaad68,3)*xyzzyaabm68
xyzzyaabj68=xyzzyaabv68(xyzzyaaac68,xyzzyaaad68,4)*xyzzyaabm68
write(xyzzyaadp1,'(e18.10,1x,e18.10,1x,e18.10,1x,e18.10)')xyzzyaabf68,&
&xyzzyaabg68,xyzzyaabi68,xyzzyaabj68
if(xyzzyaaac68>0.and.xyzzyaaac68<xyzzyaaap68)write(xyzzyaadp1,'(e18.10&
&,1x,e18.10,1x,e18.10,1x,e18.10)')xyzzyaabf68,-xyzzyaabg68,xyzzyaabi68&
&,-xyzzyaabj68
enddo
enddo
deallocate(xyzzyaabv68)
close(xyzzyaadp1)
enddo
enddo
endif
deallocate(xyzzyaabw68,xyzzyaabx68,xyzzyaabz68,xyzzyaabt68,xyzzyaabr68&
&,xyzzyaabs68)
end subroutine plot_pbackflow
subroutine pbf_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,verbos&
&e)
implicit none
integer,intent(out) :: ie,jspin
real(dp),intent(in) :: eevecs(4,netot,netot),eivecs(4,nitot,netot)
logical,intent(in) :: verbose
logical,intent(out) :: fail
integer xyzzyaaaa69,xyzzyaaab69,xyzzyaaac69,xyzzyaaad69
real(dp),parameter :: xyzzyaaae69=1.d-2
logical xyzzyaaaf69,xyzzyaaag69,xyzzyaaah69,xyzzyaaai69
logical,allocatable :: xyzzyaaaj69(:)
allocate(xyzzyaaaj69(netot),stat=xyzzyaaad69)
call check_alloc(xyzzyaaad69,'BF_ASSESS_CHECK_KINETIC','0')
xyzzyaaaj69=.false.
fail=.false.
ie=0
jspin=1
do
if(jspin>nspin)then
do xyzzyaaaa69=1,netot
if(xyzzyaaaj69(xyzzyaaaa69))then
ie=which_ie(xyzzyaaaa69)
jspin=which_spin(xyzzyaaaa69)
deallocate(xyzzyaaaj69)
return
endif
enddo
fail=.true.
deallocate(xyzzyaaaj69)
return
endif
ie=ie+1
if(ie>nele(jspin))then
ie=0
jspin=jspin+1
cycle
endif
xyzzyaaaa69=which_ii(ie,jspin)
xyzzyaaaj69(xyzzyaaaa69)=.true.
xyzzyaaaf69=.false.
if(xyzzyaaaa1.and.xyzzyaaad1<3)then
do xyzzyaaab69=1,netot
if(abs(eevecs(4,xyzzyaaab69,xyzzyaaaa69)-xyzzyaabh1(xyzzyaadm1(xyzzyaa&
&aa69,xyzzyaaab69)))<=xyzzyaaae69)then
xyzzyaaaj69(xyzzyaaaa69)=.false.
if(verbose)then
call wout('For (i,jj)=('//trim(i2s(xyzzyaaaa69))//','//trim(i2s(xyzzya&
&aab69))//'), |rij - L_eta| =',abs(eevecs(4,xyzzyaaab69,xyzzyaaaa69)-x&
&yzzyaabh1(xyzzyaadm1(xyzzyaaaa69,xyzzyaaab69))))
endif
endif
enddo
endif
if(xyzzyaaab1.or.xyzzyaaac1)then
do xyzzyaaac69=1,nitot
if(xyzzyaaad1<3)then
if(xyzzyaaab1)then
if(xyzzyaaan1(xyzzyaaac69)/=0)then
if(abs(eivecs(4,xyzzyaaac69,xyzzyaaaa69)-xyzzyaabj1(xyzzyaaan1(xyzzyaa&
&ac69)))<=xyzzyaaae69)then
xyzzyaaaj69(xyzzyaaaa69)=.false.
if(verbose)then
call wout('For (i,I)=('//trim(i2s(xyzzyaaaa69))//','//trim(i2s(xyzzyaa&
&ac69))//'), |riI - L_mu| =',abs(eivecs(4,xyzzyaaac69,xyzzyaaaa69)-xyz&
&zyaabj1(xyzzyaaan1(xyzzyaaac69))))
endif
endif
endif
endif
if(xyzzyaaac1)then
if(xyzzyaaat1(xyzzyaaac69)/=0)then
if(abs(eivecs(4,xyzzyaaac69,xyzzyaaaa69)-xyzzyaabl1(xyzzyaaat1(xyzzyaa&
&ac69)))<=xyzzyaaae69)then
xyzzyaaaj69(xyzzyaaaa69)=.false.
if(verbose)then
call wout('For (i,I)=('//trim(i2s(xyzzyaaaa69))//','//trim(i2s(xyzzyaa&
&ac69))//'), |riI - L_phi| =',abs(eivecs(4,xyzzyaaac69,xyzzyaaaa69)-xy&
&zzyaabl1(xyzzyaaat1(xyzzyaaac69))))
endif
endif
endif
endif
endif
if(eivecs(4,xyzzyaaac69,xyzzyaaaa69)<=xyzzyaaae69)then
xyzzyaaaj69(xyzzyaaaa69)=.false.
if(verbose)then
call wout('For (i,I)=('//trim(i2s(xyzzyaaaa69))//','//trim(i2s(xyzzyaa&
&ac69))//'), r_iI =',eivecs(4,xyzzyaaac69,xyzzyaaaa69))
endif
endif
enddo
endif
if(.not.xyzzyaaaj69(xyzzyaaaa69))then
if(verbose)call wout('Particle '//trim(i2s(xyzzyaaaa69))//' skipped.')
cycle
endif
xyzzyaaag69=.true.
xyzzyaaah69=.true.
xyzzyaaai69=.true.
if(xyzzyaaaa1)then
xyzzyaaag69=.false.
do xyzzyaaab69=1,netot
if(xyzzyaaab69==xyzzyaaaa69)cycle
if(eevecs(4,xyzzyaaab69,xyzzyaaaa69)<xyzzyaabh1(xyzzyaadm1(xyzzyaaaa69&
&,xyzzyaaab69)))then
xyzzyaaag69=.true.
exit
endif
enddo
endif
if(xyzzyaaab1)then
xyzzyaaah69=.false.
do xyzzyaaac69=1,nitot
if(xyzzyaaan1(xyzzyaaac69)>0)then
if(eivecs(4,xyzzyaaac69,xyzzyaaaa69)<xyzzyaabj1(xyzzyaaan1(xyzzyaaac69&
&)))then
xyzzyaaah69=.true.
exit
endif
endif
enddo
endif
if(xyzzyaaac1)then
xyzzyaaai69=.false.
do xyzzyaaac69=1,nitot
if(xyzzyaaat1(xyzzyaaac69)>0)then
if(eivecs(4,xyzzyaaac69,xyzzyaaaa69)<xyzzyaabl1(xyzzyaaat1(xyzzyaaac69&
&)))then
do xyzzyaaab69=1,netot
if(xyzzyaaab69==xyzzyaaaa69)cycle
if(eivecs(4,xyzzyaaac69,xyzzyaaab69)<xyzzyaabl1(xyzzyaaat1(xyzzyaaac69&
&)))then
xyzzyaaai69=.true.
exit
endif
enddo
endif
if(xyzzyaaai69)exit
endif
enddo
endif
xyzzyaaaf69=xyzzyaaag69.and.xyzzyaaah69.and.xyzzyaaai69
if(xyzzyaaaf69)then
deallocate(xyzzyaaaj69)
return
endif
if(verbose)call wout('Particle '//trim(i2s(xyzzyaaaa69))//' temporaril&
&y skipped.')
enddo
deallocate(xyzzyaaaj69)
end subroutine pbf_assess_check_kinetic
subroutine xyzzyaagv1
implicit none
integer xyzzyaaaa70
if(xyzzyaaaa1)then
call xyzzyaagw1(xyzzyaadm1,xyzzyaaaj1)
endif
if(xyzzyaaab1)then
do xyzzyaaaa70=1,xyzzyaaah1
call xyzzyaagx1(xyzzyaadn1(:,xyzzyaaaa70),xyzzyaaaq1(xyzzyaaaa70))
enddo
endif
if(xyzzyaaac1)then
do xyzzyaaaa70=1,xyzzyaaai1
call xyzzyaagw1(xyzzyaado1(:,:,xyzzyaaaa70),xyzzyaaay1(xyzzyaaaa70))
enddo
endif
end subroutine xyzzyaagv1
subroutine xyzzyaagw1(which_spin_pair,spin_dep)
implicit none
integer,intent(in) :: spin_dep
integer,intent(out) :: which_spin_pair(netot,netot)
integer xyzzyaaaa71,xyzzyaaab71,xyzzyaaac71,xyzzyaaad71,xyzzyaaae71,xy&
&zzyaaaf71,xyzzyaaag71,xyzzyaaah71
if(spin_dep>levels_spairs.and.am_master)call errstop('ASSIGN_SPIN_PAIR&
&','Spin dependence too large for this system.')
xyzzyaaae71=0
do xyzzyaaac71=1,nspin
do xyzzyaaaa71=1,nele(xyzzyaaac71)
xyzzyaaae71=xyzzyaaae71+1
do xyzzyaaad71=xyzzyaaac71,nspin
xyzzyaaag71=which_spair(xyzzyaaac71,xyzzyaaad71,spin_dep)
xyzzyaaaf71=0
if(xyzzyaaad71>1)xyzzyaaaf71=sum(nele(:xyzzyaaad71-1))
xyzzyaaah71=nele(xyzzyaaad71)
if(xyzzyaaad71==xyzzyaaac71)xyzzyaaah71=xyzzyaaaa71
do xyzzyaaab71=1,xyzzyaaah71
xyzzyaaaf71=xyzzyaaaf71+1
which_spin_pair(xyzzyaaae71,xyzzyaaaf71)=xyzzyaaag71
which_spin_pair(xyzzyaaaf71,xyzzyaaae71)=xyzzyaaag71
enddo
enddo
enddo
enddo
end subroutine xyzzyaagw1
subroutine xyzzyaagx1(which_spin_single,spin_dep)
implicit none
integer,intent(in) :: spin_dep
integer,intent(out) :: which_spin_single(netot)
integer xyzzyaaaa72,xyzzyaaab72,xyzzyaaac72,xyzzyaaad72
if(spin_dep>levels_spairs.and.am_master)call errstop('ASSIGN_SPIN_PAIR&
&','Spin dependence too large for this system.')
xyzzyaaac72=0
do xyzzyaaab72=1,nspin
xyzzyaaad72=which_ssingle(xyzzyaaab72,spin_dep)
do xyzzyaaaa72=1,nele(xyzzyaaab72)
xyzzyaaac72=xyzzyaaac72+1
which_spin_single(xyzzyaaac72)=xyzzyaaad72
enddo
enddo
end subroutine xyzzyaagx1
subroutine xyzzyaagy1(spin_dep,flags)
implicit none
integer,intent(in) :: spin_dep
logical,intent(out) :: flags(:)
integer xyzzyaaaa73
flags=.false.
if(noncoll_spin.or.(spin_dep==0.and..not.ferromagnetic))return
do xyzzyaaaa73=1,nspin
flags(which_spair(xyzzyaaaa73,xyzzyaaaa73,spin_dep))=.true.
enddo
end subroutine xyzzyaagy1
real(dp) function xyzzyaagz1(s)
implicit none
integer,intent(in) :: s
if(isperiodic)then
if(xyzzyaacr1(s))then
xyzzyaagz1=wigner_seitz_radius*0.5d0
else
xyzzyaagz1=wigner_seitz_radius*0.9d0
endif
else
if(xyzzyaacr1(s))then
xyzzyaagz1=1.d0
else
xyzzyaagz1=4.d0
endif
endif
end function xyzzyaagz1
real(dp) function xyzzyaaha1(maxl)
implicit none
real(dp),intent(in),optional :: maxl
real(dp) xyzzyaaaa75
if(isperiodic)then
xyzzyaaaa75=wigner_seitz_radius*0.9d0
else
xyzzyaaaa75=4.5d0
endif
if(present(maxl))then
if(maxl/=0.d0.and.xyzzyaaaa75>maxl)xyzzyaaaa75=maxl
endif
xyzzyaaha1=xyzzyaaaa75
end function xyzzyaaha1
real(dp) function xyzzyaahb1(maxl)
implicit none
real(dp),intent(in),optional :: maxl
real(dp) xyzzyaaaa76
if(isperiodic)then
xyzzyaaaa76=wigner_seitz_radius*0.45d0
else
xyzzyaaaa76=4.5d0
endif
if(present(maxl))then
if(maxl/=0.d0.and.xyzzyaaaa76>maxl)xyzzyaaaa76=maxl
endif
xyzzyaahb1=xyzzyaaaa76
end function xyzzyaahb1
end module slaarnabv

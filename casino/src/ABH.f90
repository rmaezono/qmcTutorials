module slaarnabh
use casl
use dsp
use slaarnabc
use slaarnabc_noopt
use slaarnach
use slaarnaan,only : ee_kato_gamma,en_kato_gamma
use format_utils, only : wout,i2s
use slaarnabg,     only : nitot,iontype,dimensionality
use slaarnabt,    only : exp_protect,ddot,swap1,dcopy,resize_pointer,d&
&axpy,dsum
use parallel,     only : am_master
use slaarnaca,        only : zion
use run_control,  only : errstop_master,timer,check_alloc,errwarn
use store,        only : netot,three_netot,nele,which_spin,nspin,jasbu&
&f,opt_jastrow,ee_cusp_in_orbital,inv_pmass
implicit none
private
public read_gjastrow,update_gjastrow_casl,gjastrow_assess_check_kineti&
&c,finite_size_corr_ke_gjastrow,enumerate_plot_gjastrow,query_plot_gja&
&strow,get_plot_gjastrow,finish_plot_gjastrow
public query_gjastrow_levels,query_gjastrow_level_details,setup_gjastr&
&ow,finish_gjastrow,wfn_ratio_gjastrow,accept_move_gjastrow,reset_conf&
&ig_gjastrow,add_config_gjastrow_items,clear_scratch_gjastrow,wfn_logv&
&al_gjastrow,wfn_loggrad_gjastrow,wfn_loglap_gjastrow,prefetch_wfn_gja&
&strow,setup_gjastrow_params,finish_gjastrow_params,get_gjastrow_param&
&s,put_gjastrow_params,clone_scratch_gjastrow,invalidate_params_gjastr&
&ow,invalidate_param1_gjastrow,wfn_aderiv_gjastrow,setup_storage_gjast&
&row,finish_storage_gjastrow,load_from_storage_gjastrow,save_to_storag&
&e_gjastrow,get_linear_basis_gjastrow
public gen_config_gjastrow,delete_config_gjastrow,copy_config_gjastrow&
&,config_to_pt_gjastrow,pt_to_config_gjastrow,redist_allocations_gjast&
&row,redist_load_gjastrow,redist_send_gjastrow,redist_recv_gjastrow,re&
&dist_save_gjastrow,redist_deallocations_gjastrow,load_from_pt_gjastro&
&w,save_to_pt_gjastrow
public config_wfn_gjastrow
type config_wfn_gjastrow
private
real(dp) pt_valjas
real(dp),pointer :: pt_chi(:),pt_gradchi(:,:)
logical pt_valjas_valid
logical,pointer :: pt_chi_valid(:),pt_gradchi_valid(:)
end type config_wfn_gjastrow
integer :: xyzzyaaaa1=0,xyzzyaaab1=0,xyzzyaaac1=0,xyzzyaaad1=0,xyzzyaa&
&ae1=0,xyzzyaaaf1=0
integer,allocatable :: xyzzyaaag1(:),xyzzyaaah1(:),xyzzyaaai1(:),xyzzy&
&aaaj1(:)
integer,allocatable :: xyzzyaaak1(:),xyzzyaaal1(:),xyzzyaaam1(:),xyzzy&
&aaan1(:),xyzzyaaao1(:),xyzzyaaap1(:),xyzzyaaaq1(:),xyzzyaaar1(:),xyzz&
&yaaas1(:),xyzzyaaat1(:),xyzzyaaau1(:),xyzzyaaav1(:),xyzzyaaaw1(:),xyz&
&zyaaax1(:),xyzzyaaay1(:,:,:),xyzzyaaaz1(:,:,:),xyzzyaaba1(:),xyzzyaab&
&b1(:),xyzzyaabc1(:),xyzzyaabd1(:),xyzzyaabe1(:),xyzzyaabf1(:),xyzzyaa&
&bg1(:),xyzzyaabh1(:),xyzzyaabi1(:)
real(dp),allocatable :: xyzzyaabj1(:)
logical,allocatable :: xyzzyaabk1(:),xyzzyaabl1(:),xyzzyaabm1(:),xyzzy&
&aabn1(:),xyzzyaabo1(:),xyzzyaabp1(:),xyzzyaabq1(:),xyzzyaabr1(:),xyzz&
&yaabs1(:),xyzzyaabt1(:),xyzzyaabu1(:),xyzzyaabv1(:),xyzzyaabw1(:)
integer,allocatable :: xyzzyaabx1(:,:,:),xyzzyaaby1(:,:,:)
integer,allocatable :: xyzzyaabz1(:)
real(dp),allocatable :: xyzzyaaca1(:,:),xyzzyaacb1(:,:)
logical,allocatable :: xyzzyaacc1(:,:),xyzzyaacd1(:,:)
character(casl_keysize),pointer :: xyzzyaace1(:)=>null()
integer,allocatable :: xyzzyaacf1(:,:)
integer,allocatable :: xyzzyaacg1(:)
character(casl_keysize),pointer :: xyzzyaach1(:)=>null()
integer,allocatable :: xyzzyaaci1(:),xyzzyaacj1(:),xyzzyaack1(:),xyzzy&
&aacl1(:),xyzzyaacm1(:),xyzzyaacn1(:)
real(dp),allocatable :: param(:),xyzzyaaco1(:),xyzzyaacp1(:)
logical,allocatable :: xyzzyaacq1(:),xyzzyaacr1(:),xyzzyaacs1(:),xyzzy&
&aact1(:),xyzzyaacu1(:)
integer :: xyzzyaacv1=0,xyzzyaacw1=0
integer,allocatable :: xyzzyaacx1(:),xyzzyaacy1(:),xyzzyaacz1(:),xyzzy&
&aada1(:),xyzzyaadb1(:),xyzzyaadc1(:),xyzzyaadd1(:),xyzzyaade1(:),xyzz&
&yaadf1(:),xyzzyaadg1(:),xyzzyaadh1(:),xyzzyaadi1(:),xyzzyaadj1(:),xyz&
&zyaadk1(:)
integer,allocatable :: xyzzyaadl1(:),xyzzyaadm1(:),xyzzyaadn1(:),xyzzy&
&aado1(:)
integer,allocatable :: xyzzyaadp1(:),xyzzyaadq1(:)
real(dp),allocatable :: xyzzyaadr1(:),xyzzyaads1(:),eqn_rhs(:),xyzzyaa&
&dt1(:)
integer,parameter :: xyzzyaadu1=2,xyzzyaadv1=1,xyzzyaadw1=2
real(dp),allocatable :: xyzzyaadx1(:),xyzzyaady1(:,:),xyzzyaadz1(:,:,:&
&),xyzzyaaea1(:,:,:),xyzzyaaeb1(:,:),xyzzyaaec1(:,:),xyzzyaaed1(:,:),x&
&yzzyaaee1(:,:,:),xyzzyaaef1(:,:,:,:),xyzzyaaeg1(:,:,:),xyzzyaaeh1(:,:&
&,:),xyzzyaaei1(:,:,:,:),xyzzyaaej1(:,:,:)
logical,allocatable :: xyzzyaaek1(:,:),xyzzyaael1(:,:),xyzzyaaem1(:,:)&
&,xyzzyaaen1(:,:),xyzzyaaeo1(:,:),xyzzyaaep1(:),xyzzyaaeq1(:,:),xyzzya&
&aer1(:,:),xyzzyaaes1(:,:),xyzzyaaet1(:,:,:),xyzzyaaeu1(:,:,:)
real(dp),allocatable :: xyzzyaaev1(:)
logical :: xyzzyaaew1
real(dp),allocatable :: xyzzyaaex1(:),xyzzyaaey1(:,:,:),xyzzyaaez1(:,:&
&),xyzzyaafa1(:,:),xyzzyaafb1(:,:)
logical,allocatable :: xyzzyaafc1(:),xyzzyaafd1(:),xyzzyaafe1(:),xyzzy&
&aaff1(:)
integer,allocatable :: xyzzyaafg1(:),xyzzyaafh1(:),  xyzzyaafi1(:),xyz&
&zyaafj1(:),xyzzyaafk1(:),xyzzyaafl1(:),xyzzyaafm1(:),xyzzyaafn1(:),xy&
&zzyaafo1(:),xyzzyaafp1(:),xyzzyaafq1(:),xyzzyaafr1(:)
integer,allocatable :: xyzzyaafs1(:)
real(dp),allocatable :: xyzzyaaft1(:),xyzzyaafu1(:),xyzzyaafv1(:)
integer,parameter :: xyzzyaafw1        = 1,xyzzyaafx1        = 2,xyzzy&
&aafy1        = 3,xyzzyaafz1        = 4,xyzzyaaga1        = 5,xyzzyaag&
&b1        = 6,xyzzyaagc1    = 11,xyzzyaagd1    = 12,xyzzyaage1    = 1&
&3
integer,parameter :: xyzzyaagf1      = 1,xyzzyaagg1      = 2,xyzzyaagh&
&1      = 3
integer,parameter :: nthings_to_plot=19,xyzzyaagi1=1,xyzzyaagj1=2,xyzz&
&yaagk1=3,xyzzyaagl1=4,xyzzyaagm1=5,xyzzyaagn1=6,xyzzyaago1=7,xyzzyaag&
&p1=8,xyzzyaagq1=9,xyzzyaagr1=10,xyzzyaags1=11,xyzzyaagt1=12,xyzzyaagu&
&1=13,xyzzyaagv1=14,xyzzyaagw1=15,xyzzyaagx1=16,xyzzyaagy1=17,xyzzyaag&
&z1=18,xyzzyaaha1=19
integer :: xyzzyaahb1(nthings_to_plot)=0
character(64),parameter :: xyzzyaahc1(nthings_to_plot)=(/'gjastrow    &
&         ',  'gjastrow_oneelec     ',  'gjastrow_grad        ',  'gja&
&strow_lap         ',  'gjastrow_terms       ',  'gjastrow_terms_grad &
& ',  'gjastrow_terms_lap   ',  'gjastrow_eebasis     ',  'gjastrow_ee&
&basis_grad',  'gjastrow_eebasis_lap ',  'gjastrow_enbasis     ',  'gj&
&astrow_enbasis_grad',  'gjastrow_enbasis_lap ',  'gjastrow_eecut     &
&  ',  'gjastrow_eecut_grad  ',  'gjastrow_eecut_lap   ',  'gjastrow_e&
&ncut       ',  'gjastrow_encut_grad  ',  'gjastrow_encut_lap   '/)
character(64),parameter :: plot_description(nthings_to_plot)=(/'Jastro&
&w factor (gjastrow)                          ',  'one-electron Jastro&
&w factor (gjastrow)             ',  'gradient of Jastrow factor (gjas&
&trow)              ',  'Laplacian of Jastrow factor (gjastrow)       &
&      ',  'one-electron Jastrow factor terms (gjastrow)       ',  'gr&
&adient of Jastrow factor terms (gjastrow)        ',  'Laplacian of Ja&
&strow factor terms (gjastrow)       ',  'Jastrow e-e basis functions &
&(gjastrow)             ',  'gradient of Jastrow e-e basis functions (&
&gjastrow) ',  'Laplacian of Jastrow e-e basis functions (gjastrow)', &
& 'Jastrow e-n basis functions (gjastrow)             ',  'gradient of&
& Jastrow e-n basis functions (gjastrow) ',  'Laplacian of Jastrow e-n&
& basis functions (gjastrow)',  'Jastrow factor e-e cut-offs (gjastrow&
&)             ',  'gradient of Jastrow factor e-e cut-offs (gjastrow)&
& ',  'Laplacian of Jastrow factor e-e cut-offs (gjastrow)',  'Jastrow&
& factor e-n cut-offs (gjastrow)             ',  'gradient of Jastrow &
&factor e-n cut-offs (gjastrow) ',  'Laplacian of Jastrow factor e-n c&
&ut-offs (gjastrow)'/)
logical :: xyzzyaahd1=.false.
logical,parameter :: xyzzyaahe1=.false.
logical,parameter :: xyzzyaahf1=.false.
logical,parameter :: xyzzyaahg1=.false.
logical,parameter :: xyzzyaahh1=.false.
logical,parameter :: xyzzyaahi1=.false.
logical,parameter :: xyzzyaahj1=.false.
logical,parameter :: xyzzyaahk1=.false.
logical,parameter :: xyzzyaahl1=.false.
logical,parameter :: xyzzyaahm1=.false.
contains
subroutine query_gjastrow_levels(nlevels)
implicit none
integer,intent(out) :: nlevels
nlevels=1
end subroutine query_gjastrow_levels
subroutine query_gjastrow_level_details(level_score,level_name)
implicit none
integer,intent(out) :: level_score(*)
character(80),intent(out) :: level_name(*)
level_score(1)=70
level_name(1)='Jastrow factor (gjastrow)'
end subroutine query_gjastrow_level_details
subroutine setup_gjastrow
implicit none
integer xyzzyaaaa4(2),xyzzyaaab4(2),xyzzyaaac4(2),xyzzyaaad4(2),xyzzya&
&aae4(2),xyzzyaaaf4(2),xyzzyaaag4
xyzzyaaaa4=0
xyzzyaaab4=0
xyzzyaaac4=0
xyzzyaaad4=0
xyzzyaaae4=0
xyzzyaaaf4=0
call include_range(ratiocfg_from_sz,xyzzyaaaa4)
call include_range(ratiocfg_to_sz,xyzzyaaaa4)
call include_range(ratio1_from_sz,xyzzyaaaa4)
call include_range(ratio1_to_sz,xyzzyaaaa4)
call include_range(ratio2_from_sz,xyzzyaaaa4)
call include_range(ratio2_to_sz,xyzzyaaaa4)
call include_range(ratio_ion_from_sz,xyzzyaaaa4)
call include_range(ratio_ion_to_sz,xyzzyaaaa4)
call include_range(drift_sz,xyzzyaaab4)
call include_range(kinetic_sz,xyzzyaaab4)
call include_range(kinetic_sz,xyzzyaaac4)
call include_range(wfn_detail_sz,xyzzyaaaa4)
call include_range(wfn_detail_sz,xyzzyaaad4)
call include_range(wfn_detail_sz,xyzzyaaaf4)
call include_range(kinetic_detail_sz,xyzzyaaab4)
call include_range(kinetic_detail_sz,xyzzyaaac4)
call include_range(kinetic_detail_sz,xyzzyaaae4)
call include_range(kinetic_detail_sz,xyzzyaaaf4)
if(xyzzyaaaa4(1)/=0)then
allocate(xyzzyaadx1(xyzzyaaaa4(1):xyzzyaaaa4(2)),xyzzyaady1(netot,xyzz&
&yaaaa4(1):xyzzyaaaa4(2)),xyzzyaaec1(netot,xyzzyaaaa4(1):xyzzyaaaa4(2)&
&),xyzzyaaed1(netot,xyzzyaaaa4(1):xyzzyaaaa4(2)),stat=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','jas')
xyzzyaadx1=0.d0
xyzzyaady1=0.d0
xyzzyaaec1=0.d0
xyzzyaaed1=0.d0
endif
allocate(xyzzyaaep1(nscratch),xyzzyaaen1(netot,nscratch),xyzzyaaek1(ne&
&tot,nscratch),stat=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','jas_valid')
if(xyzzyaaab4(1)/=0)then
allocate(xyzzyaadz1(3,netot,xyzzyaaab4(1):xyzzyaaab4(2)),xyzzyaaea1(3,&
&netot,xyzzyaaab4(1):xyzzyaaab4(2)),stat=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','gradjas')
xyzzyaadz1=0.d0
xyzzyaaea1=0.d0
endif
allocate(xyzzyaael1(netot,nscratch),xyzzyaaeo1(netot,nscratch),stat=xy&
&zzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','gradjas_valid')
if(xyzzyaaac4(1)/=0)then
allocate(xyzzyaaeb1(netot,xyzzyaaac4(1):xyzzyaaac4(2)),stat=xyzzyaaag4&
&)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','lapjas')
xyzzyaaeb1=0.d0
endif
allocate(xyzzyaaem1(netot,nscratch),stat=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','lapjas_valid')
if(xyzzyaaad4(1)/=0)then
allocate(xyzzyaaeh1(netot,xyzzyaaaa1,xyzzyaaad4(1):xyzzyaaad4(2)),stat&
&=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','jas1term')
xyzzyaaeh1=0.d0
endif
allocate(xyzzyaaet1(netot,xyzzyaaaa1,nscratch),stat=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','jas1term_valid')
if(xyzzyaaae4(1)/=0)then
allocate(xyzzyaaei1(3,netot,xyzzyaaaa1,xyzzyaaae4(1):xyzzyaaae4(2)),xy&
&zzyaaej1(netot,xyzzyaaaa1,xyzzyaaae4(1):xyzzyaaae4(2)),stat=xyzzyaaag&
&4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','gljas1term')
xyzzyaaei1=0.d0
xyzzyaaej1=0.d0
endif
allocate(xyzzyaaeu1(netot,xyzzyaaaa1,nscratch),stat=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','gljas1term_valid')
if(xyzzyaaaf4(1)/=0)then
allocate(xyzzyaaee1(xyzzyaaad1,netot,xyzzyaaaf4(1):xyzzyaaaf4(2)),xyzz&
&yaaef1(3,xyzzyaaad1,netot,xyzzyaaaf4(1):xyzzyaaaf4(2)),xyzzyaaeg1(xyz&
&zyaaad1,netot,xyzzyaaaf4(1):xyzzyaaaf4(2)),stat=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','djas1')
xyzzyaaee1=0.d0
xyzzyaaef1=0.d0
xyzzyaaeg1=0.d0
endif
allocate(xyzzyaaeq1(netot,nscratch),xyzzyaaer1(netot,nscratch),xyzzyaa&
&es1(netot,nscratch),stat=xyzzyaaag4)
call check_alloc(xyzzyaaag4,'SETUP_JASTROW','djas1_valid')
end subroutine setup_gjastrow
subroutine finish_gjastrow
implicit none
if(allocated(xyzzyaadx1))deallocate(xyzzyaadx1,xyzzyaady1,xyzzyaaec1,x&
&yzzyaaed1)
deallocate(xyzzyaaep1,xyzzyaaen1,xyzzyaaek1)
if(allocated(xyzzyaadz1))deallocate(xyzzyaadz1,xyzzyaaea1)
deallocate(xyzzyaael1,xyzzyaaeo1)
if(allocated(xyzzyaaeb1))deallocate(xyzzyaaeb1)
deallocate(xyzzyaaem1)
if(allocated(xyzzyaaeh1))deallocate(xyzzyaaeh1)
deallocate(xyzzyaaet1)
if(allocated(xyzzyaaei1))deallocate(xyzzyaaei1,xyzzyaaej1)
deallocate(xyzzyaaeu1)
if(allocated(xyzzyaaee1))deallocate(xyzzyaaee1,xyzzyaaef1,xyzzyaaeg1)
deallocate(xyzzyaaeq1,xyzzyaaer1,xyzzyaaes1)
end subroutine finish_gjastrow
subroutine wfn_ratio_gjastrow(is,js,ilevel,ratio,fd,sd)
implicit none
integer,intent(in) :: is,js,ilevel
real(dp),intent(out) :: ratio
logical,intent(in) :: fd,sd
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaa&
&af6,xyzzyaaag6,xyzzyaaah6,xyzzyaaai6,xyzzyaaaj6,xyzzyaaak6
real(dp) xyzzyaaal6,xyzzyaaam6,xyzzyaaan6,xyzzyaaao6,jas
if(buffer_move1_from(js)==is)then
xyzzyaaaa6=buffer_move1_from_ii(js)
call xyzzyaahn1(xyzzyaaaa6,is,.true.,fd,sd,.false.)
call xyzzyaahn1(xyzzyaaaa6,js,.true.,fd,sd,.false.)
ratio=exp_protect(xyzzyaaec1(xyzzyaaaa6,js)-xyzzyaaec1(xyzzyaaaa6,is))
elseif(buffer_move2_from(js)==is)then
xyzzyaaaa6=buffer_move2_from_ii(js)
xyzzyaaab6=buffer_move2_from_jj(js)
call get_gbasis(is,.false.,.false.)
call get_gbasis(js,.false.,.false.)
call xyzzyaahn1(xyzzyaaaa6,is,.true.,.false.,.false.,.false.)
xyzzyaaal6=xyzzyaaec1(xyzzyaaaa6,is)
xyzzyaaao6=0.d0
xyzzyaaan6=0.d0
xyzzyaaam6=0.d0
do xyzzyaaac6=1,xyzzyaaaa1
call timer('TERM_'//trim(i2s(xyzzyaaac6)),.true.)
if(xyzzyaaao1(xyzzyaaac6)>0)then
xyzzyaaad6=xyzzyaaaq1(xyzzyaaac6)
xyzzyaaah6=ifn1_eebasis(xyzzyaaad6)
xyzzyaaae6=xyzzyaaar1(xyzzyaaac6)
xyzzyaaai6=ifn1_eebasis(xyzzyaaae6)
else
xyzzyaaad6=1
xyzzyaaah6=1
xyzzyaaae6=1
xyzzyaaai6=1
endif
if(xyzzyaaap1(xyzzyaaac6)>0)then
xyzzyaaaf6=xyzzyaaas1(xyzzyaaac6)
xyzzyaaaj6=ifn1_enbasis(xyzzyaaaf6)
xyzzyaaag6=xyzzyaaat1(xyzzyaaac6)
xyzzyaaak6=ifn1_enbasis(xyzzyaaag6)
else
xyzzyaaaf6=1
xyzzyaaaj6=1
xyzzyaaag6=1
xyzzyaaak6=1
endif
call xyzzyaaid1(xyzzyaaab6,xyzzyaaac6,eebasis_scr(xyzzyaaah6,1,xyzzyaa&
&ab6,js),eebasis_scr(xyzzyaaai6,1,xyzzyaaab6,js),nzeecut_scr(1,xyzzyaa&
&ab6,xyzzyaaae6,js),enbasis_scr(xyzzyaaaj6,1,xyzzyaaab6,js),enbasis_sc&
&r(xyzzyaaak6,1,xyzzyaaab6,js),nzencut_scr(1,xyzzyaaab6,xyzzyaaag6,js)&
&,eebasis_scr(xyzzyaaah6,1,1,js),eebasis_scr(xyzzyaaai6,1,1,js),nzeecu&
&t_scr(1,1,xyzzyaaae6,js),enbasis_scr(xyzzyaaaj6,1,1,js),enbasis_scr(x&
&yzzyaaak6,1,1,js),nzencut_scr(1,1,xyzzyaaag6,js),jas=jas)
xyzzyaaao6=xyzzyaaao6+jas
call xyzzyaaid1(xyzzyaaab6,xyzzyaaac6,eebasis_scr(xyzzyaaah6,1,xyzzyaa&
&ab6,is),eebasis_scr(xyzzyaaai6,1,xyzzyaaab6,is),nzeecut_scr(1,xyzzyaa&
&ab6,xyzzyaaae6,is),enbasis_scr(xyzzyaaaj6,1,xyzzyaaab6,is),enbasis_sc&
&r(xyzzyaaak6,1,xyzzyaaab6,is),nzencut_scr(1,xyzzyaaab6,xyzzyaaag6,is)&
&,eebasis_scr(xyzzyaaah6,1,1,js),eebasis_scr(xyzzyaaai6,1,1,js),nzeecu&
&t_scr(1,1,xyzzyaaae6,js),enbasis_scr(xyzzyaaaj6,1,1,js),enbasis_scr(x&
&yzzyaaak6,1,1,js),nzencut_scr(1,1,xyzzyaaag6,js),jas=jas)
xyzzyaaan6=xyzzyaaan6+jas
call xyzzyaaid1(xyzzyaaaa6,xyzzyaaac6,eebasis_scr(xyzzyaaah6,1,xyzzyaa&
&aa6,js),eebasis_scr(xyzzyaaai6,1,xyzzyaaaa6,js),nzeecut_scr(1,xyzzyaa&
&aa6,xyzzyaaae6,js),enbasis_scr(xyzzyaaaj6,1,xyzzyaaaa6,js),enbasis_sc&
&r(xyzzyaaak6,1,xyzzyaaaa6,js),nzencut_scr(1,xyzzyaaaa6,xyzzyaaag6,js)&
&,eebasis_scr(xyzzyaaah6,1,1,is),eebasis_scr(xyzzyaaai6,1,1,is),nzeecu&
&t_scr(1,1,xyzzyaaae6,is),enbasis_scr(xyzzyaaaj6,1,1,is),enbasis_scr(x&
&yzzyaaak6,1,1,is),nzencut_scr(1,1,xyzzyaaag6,is),jas=jas)
xyzzyaaam6=xyzzyaaam6+jas
call timer('TERM_'//trim(i2s(xyzzyaaac6)),.false.)
enddo
ratio=exp_protect(xyzzyaaam6-xyzzyaaal6+xyzzyaaao6-xyzzyaaan6)
else
call xyzzyaaho1(is)
call xyzzyaaho1(js)
ratio=exp_protect(xyzzyaadx1(js)-xyzzyaadx1(is))
endif
end subroutine wfn_ratio_gjastrow
subroutine accept_move_gjastrow(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa7,xyzzyaaab7
if(buffer_move1_from(js)==is)then
xyzzyaaaa7=buffer_move1_from_ii(js)
if(xyzzyaaen1(xyzzyaaaa7,js))xyzzyaady1(xyzzyaaaa7,is)=xyzzyaady1(xyzz&
&yaaaa7,js)
if(xyzzyaaeo1(xyzzyaaaa7,js))xyzzyaaea1(1:3,xyzzyaaaa7,is)=xyzzyaaea1(&
&1:3,xyzzyaaaa7,js)
xyzzyaaen1(xyzzyaaaa7,is)=xyzzyaaen1(xyzzyaaaa7,js)
xyzzyaaeo1(xyzzyaaaa7,is)=xyzzyaaeo1(xyzzyaaaa7,js)
xyzzyaaek1(:,is)=.false.
xyzzyaael1(:,is)=.false.
xyzzyaaem1(:,is)=.false.
if(scr_tasks(iwfn_detail,is).or.scr_tasks(ikinetic_detail,is))then
xyzzyaaet1(:,:,is)=.false.
xyzzyaaeu1(:,:,is)=.false.
endif
elseif(buffer_move2_from(js)==is)then
xyzzyaaaa7=buffer_move2_from_ii(js)
xyzzyaaab7=buffer_move2_from_jj(js)
if(xyzzyaaen1(xyzzyaaaa7,js))xyzzyaady1(xyzzyaaaa7,is)=xyzzyaady1(xyzz&
&yaaaa7,js)
if(xyzzyaaen1(xyzzyaaab7,js))xyzzyaady1(xyzzyaaab7,is)=xyzzyaady1(xyzz&
&yaaab7,js)
if(xyzzyaaeo1(xyzzyaaaa7,js))xyzzyaaea1(1:3,xyzzyaaaa7,is)=xyzzyaaea1(&
&1:3,xyzzyaaaa7,js)
if(xyzzyaaeo1(xyzzyaaab7,js))xyzzyaaea1(1:3,xyzzyaaab7,is)=xyzzyaaea1(&
&1:3,xyzzyaaab7,js)
xyzzyaaen1(xyzzyaaaa7,is)=xyzzyaaen1(xyzzyaaaa7,js)
xyzzyaaen1(xyzzyaaab7,is)=xyzzyaaen1(xyzzyaaab7,js)
xyzzyaaeo1(xyzzyaaaa7,is)=xyzzyaaeo1(xyzzyaaaa7,js)
xyzzyaaeo1(xyzzyaaab7,is)=xyzzyaaeo1(xyzzyaaab7,js)
xyzzyaaek1(:,is)=.false.
xyzzyaael1(:,is)=.false.
xyzzyaaem1(:,is)=.false.
if(scr_tasks(iwfn_detail,is).or.scr_tasks(ikinetic_detail,is))then
xyzzyaaet1(:,:,is)=.false.
xyzzyaaeu1(:,:,is)=.false.
endif
else
do xyzzyaaaa7=1,netot
if(xyzzyaaek1(xyzzyaaaa7,js))then
xyzzyaaec1(xyzzyaaaa7,is)=xyzzyaaec1(xyzzyaaaa7,js)
xyzzyaaed1(xyzzyaaaa7,is)=xyzzyaaed1(xyzzyaaaa7,js)
endif
if(xyzzyaael1(xyzzyaaaa7,js))xyzzyaadz1(1:3,xyzzyaaaa7,is)=xyzzyaadz1(&
&1:3,xyzzyaaaa7,js)
if(xyzzyaaem1(xyzzyaaaa7,js))xyzzyaaeb1(xyzzyaaaa7,is)=xyzzyaaeb1(xyzz&
&yaaaa7,js)
if(xyzzyaaen1(xyzzyaaaa7,js))xyzzyaady1(xyzzyaaaa7,is)=xyzzyaady1(xyzz&
&yaaaa7,js)
if(xyzzyaaeo1(xyzzyaaaa7,js))xyzzyaaea1(1:3,xyzzyaaaa7,is)=xyzzyaaea1(&
&1:3,xyzzyaaaa7,js)
xyzzyaaek1(xyzzyaaaa7,is)=xyzzyaaek1(xyzzyaaaa7,js)
xyzzyaael1(xyzzyaaaa7,is)=xyzzyaael1(xyzzyaaaa7,js)
xyzzyaaem1(xyzzyaaaa7,is)=xyzzyaaem1(xyzzyaaaa7,js)
xyzzyaaen1(xyzzyaaaa7,is)=xyzzyaaen1(xyzzyaaaa7,js)
xyzzyaaeo1(xyzzyaaaa7,is)=xyzzyaaeo1(xyzzyaaaa7,js)
if(scr_tasks(iwfn_detail,is).or.scr_tasks(ikinetic_detail,is))then
xyzzyaaet1(:,:,is)=.false.
xyzzyaaeu1(:,:,is)=.false.
endif
enddo
endif
if(xyzzyaaep1(js))then
xyzzyaadx1(is)=xyzzyaadx1(js)
xyzzyaaep1(is)=.true.
else
xyzzyaaep1(is)=.false.
endif
xyzzyaaeq1(:,is)=.false.
xyzzyaaer1(:,is)=.false.
xyzzyaaes1(:,is)=.false.
end subroutine accept_move_gjastrow
subroutine reset_config_gjastrow(is,js)
implicit none
integer,intent(in) :: is,js
call clear_scratch_gjastrow(js)
end subroutine reset_config_gjastrow
subroutine wfn_logval_gjastrow(is,logwfn)
implicit none
integer,intent(in) :: is
real(dp),intent(out) :: logwfn
call xyzzyaaho1(is)
logwfn=xyzzyaadx1(is)
end subroutine wfn_logval_gjastrow
subroutine wfn_loggrad_gjastrow(ii,is,ilevel,val,sd,aderiv,gradjas)
implicit none
integer,intent(in) :: ii,is,ilevel
real(dp),intent(out) :: gradjas(3)
logical,intent(in) :: val,sd,aderiv
call xyzzyaahn1(ii,is,val,.true.,sd,aderiv)
gradjas=xyzzyaadz1(1:3,ii,is)
end subroutine wfn_loggrad_gjastrow
subroutine wfn_loglap_gjastrow(ii,is,ilevel,val,fd,aderiv,lapjas)
implicit none
integer,intent(in) :: ii,is,ilevel
real(dp),intent(out) :: lapjas
logical,intent(in) :: val,fd,aderiv
call xyzzyaahn1(ii,is,val,fd,.true.,aderiv)
lapjas=xyzzyaaeb1(ii,is)
end subroutine wfn_loglap_gjastrow
subroutine prefetch_wfn_gjastrow(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
call get_gbasis(is,fd,sd)
end subroutine prefetch_wfn_gjastrow
subroutine xyzzyaahn1(ii,is,val,fd,sd,aderiv)
implicit none
integer,intent(in) :: ii,is
logical,intent(in) :: val,fd,sd,aderiv
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13,xyzzyaaae13,xy&
&zzyaaaf13,xyzzyaaag13,xyzzyaaah13,xyzzyaaai13,xyzzyaaaj13,xyzzyaaak13&
&,xyzzyaaal13,xyzzyaaam13
real(dp) xyzzyaaan13,xyzzyaaao13,xyzzyaaap13(3),xyzzyaaaq13,xyzzyaaar1&
&3,xyzzyaaas13(3),jas,gjas(3),ljas,djas(xyzzyaaad1),dgjas(3,xyzzyaaad1&
&),dljas(xyzzyaaad1)
logical need_val,xyzzyaaat13,xyzzyaaau13,xyzzyaaav13,xyzzyaaaw13,xyzzy&
&aaax13,xyzzyaaay13,detail,xyzzyaaaz13
need_val=val.and..not.xyzzyaaek1(ii,is)
xyzzyaaat13=fd.and..not.xyzzyaael1(ii,is)
xyzzyaaau13=sd.and..not.xyzzyaaem1(ii,is)
if(aderiv)then
need_val=val.and.(need_val.or..not.xyzzyaaeq1(ii,is))
xyzzyaaat13=fd.and.(xyzzyaaat13.or..not.xyzzyaaer1(ii,is))
xyzzyaaau13=sd.and.(xyzzyaaau13.or..not.xyzzyaaes1(ii,is))
endif
xyzzyaaav13=xyzzyaaat13.or.xyzzyaaau13
if(.not.(need_val.or.xyzzyaaav13))return
call timer('GET_JASTROW1',.true.)
xyzzyaaan13=0.d0
xyzzyaaao13=0.d0
xyzzyaaap13=0.d0
xyzzyaaaq13=0.d0
xyzzyaaar13=0.d0
xyzzyaaas13=0.d0
detail=xyzzyaaat13.and.xyzzyaaau13.and.scr_tasks(ikinetic_detail,is)
if(detail.and.need_val.and..not.scr_tasks(iwfn_detail,is))call errstop&
&_master('GET_JASTROW1','Attempted to get Jastrow value using detailed&
& version without having detail arrays for the wave function.')
xyzzyaaaz13=.not.xyzzyaaau13.and.buffer_move1_from_ii(is)==ii
if(xyzzyaaaz13)then
xyzzyaaab13=buffer_move1_from(is)
call get_gbasis(xyzzyaaab13,xyzzyaaat13,.false.)
call get_gbasis1_ch(ii,is,xyzzyaaat13)
else
call get_gbasis(is,xyzzyaaav13,xyzzyaaau13)
endif
do xyzzyaaaa13=1,xyzzyaaaa1
call timer('TERM_'//trim(i2s(xyzzyaaaa13)),.true.)
if(aderiv)then
xyzzyaaaw13=need_val
xyzzyaaax13=xyzzyaaav13
xyzzyaaay13=xyzzyaaav13
elseif(detail)then
xyzzyaaaw13=need_val.and..not.xyzzyaaet1(ii,xyzzyaaaa13,is)
xyzzyaaax13=xyzzyaaav13.and..not.xyzzyaaeu1(ii,xyzzyaaaa13,is)
xyzzyaaay13=xyzzyaaax13
elseif(xyzzyaaao1(xyzzyaaaa13)==0)then
jas=0.d0
gjas=0.d0
xyzzyaaaw13=need_val.and..not.xyzzyaaen1(ii,is)
xyzzyaaax13=xyzzyaaat13.and..not.xyzzyaaeo1(ii,is)
xyzzyaaay13=xyzzyaaau13
else
xyzzyaaaw13=need_val
xyzzyaaax13=xyzzyaaat13
xyzzyaaay13=xyzzyaaau13
endif
if(xyzzyaaao1(xyzzyaaaa13)>0)then
xyzzyaaac13=xyzzyaaaq1(xyzzyaaaa13)
xyzzyaaag13=ifn1_eebasis(xyzzyaaac13)
xyzzyaaad13=xyzzyaaar1(xyzzyaaaa13)
xyzzyaaah13=ifn1_eebasis(xyzzyaaad13)
else
xyzzyaaac13=1
xyzzyaaag13=1
xyzzyaaad13=1
xyzzyaaah13=1
endif
if(xyzzyaaap1(xyzzyaaaa13)>0)then
xyzzyaaae13=xyzzyaaas1(xyzzyaaaa13)
xyzzyaaai13=ifn1_enbasis(xyzzyaaae13)
xyzzyaaaf13=xyzzyaaat1(xyzzyaaaa13)
xyzzyaaaj13=ifn1_enbasis(xyzzyaaaf13)
else
xyzzyaaae13=1
xyzzyaaai13=1
xyzzyaaaf13=1
xyzzyaaaj13=1
endif
if(xyzzyaaaz13)then
if(xyzzyaaay13)then
call xyzzyaaid1(ii,xyzzyaaaa13,eebasis1_chscr(xyzzyaaag13,1,is),eebasi&
&s1_chscr(xyzzyaaah13,1,is),nzeecut1_chscr(1,xyzzyaaad13,is),enbasis1_&
&chscr(xyzzyaaai13,1,is),enbasis1_chscr(xyzzyaaaj13,1,is),nzencut1_chs&
&cr(1,xyzzyaaaf13,is),eebasis_scr(xyzzyaaag13,1,1,xyzzyaaab13),eebasis&
&_scr(xyzzyaaah13,1,1,xyzzyaaab13),nzeecut_scr(1,1,xyzzyaaad13,xyzzyaa&
&ab13),enbasis_scr(xyzzyaaai13,1,1,xyzzyaaab13),enbasis_scr(xyzzyaaaj1&
&3,1,1,xyzzyaaab13),nzencut_scr(1,1,xyzzyaaaf13,xyzzyaaab13),grad_eeba&
&sis1_chscr(1,xyzzyaaag13,1,is),grad_eebasis1_chscr(1,xyzzyaaah13,1,is&
&),grad_enbasis1_chscr(1,xyzzyaaai13,1,is),grad_enbasis1_chscr(1,xyzzy&
&aaaj13,1,is),deebasis1_chscr(xyzzyaaag13,1,is),deebasis1_chscr(xyzzya&
&aah13,1,is),denbasis1_chscr(xyzzyaaai13,1,is),denbasis1_chscr(xyzzyaa&
&aj13,1,is),gradr_eebasis1_chscr(1,1,is),gradr_enbasis1_chscr(1,1,is),&
&lap_eebasis_scr(xyzzyaaag13,1,ii,xyzzyaaab13),lap_eebasis_scr(xyzzyaa&
&ah13,1,ii,xyzzyaaab13),lap_enbasis_scr(xyzzyaaai13,1,ii,xyzzyaaab13),&
&lap_enbasis_scr(xyzzyaaaj13,1,ii,xyzzyaaab13),d2eebasis_scr(xyzzyaaag&
&13,1,ii,xyzzyaaab13),d2eebasis_scr(xyzzyaaah13,1,ii,xyzzyaaab13),d2en&
&basis_scr(xyzzyaaai13,1,ii,xyzzyaaab13),d2enbasis_scr(xyzzyaaaj13,1,i&
&i,xyzzyaaab13),lapr_eebasis_scr(1,ii,is),lapr_enbasis_scr(1,ii,is),ja&
&s=jas,gjas=gjas,ljas=ljas)
elseif(xyzzyaaax13)then
call xyzzyaaid1(ii,xyzzyaaaa13,eebasis1_chscr(xyzzyaaag13,1,is),eebasi&
&s1_chscr(xyzzyaaah13,1,is),nzeecut1_chscr(1,xyzzyaaad13,is),enbasis1_&
&chscr(xyzzyaaai13,1,is),enbasis1_chscr(xyzzyaaaj13,1,is),nzencut1_chs&
&cr(1,xyzzyaaaf13,is),eebasis_scr(xyzzyaaag13,1,1,xyzzyaaab13),eebasis&
&_scr(xyzzyaaah13,1,1,xyzzyaaab13),nzeecut_scr(1,1,xyzzyaaad13,xyzzyaa&
&ab13),enbasis_scr(xyzzyaaai13,1,1,xyzzyaaab13),enbasis_scr(xyzzyaaaj1&
&3,1,1,xyzzyaaab13),nzencut_scr(1,1,xyzzyaaaf13,xyzzyaaab13),grad_eeba&
&sis1_chscr(1,xyzzyaaag13,1,is),grad_eebasis1_chscr(1,xyzzyaaah13,1,is&
&),grad_enbasis1_chscr(1,xyzzyaaai13,1,is),grad_enbasis1_chscr(1,xyzzy&
&aaaj13,1,is),deebasis1_chscr(xyzzyaaag13,1,is),deebasis1_chscr(xyzzya&
&aah13,1,is),denbasis1_chscr(xyzzyaaai13,1,is),denbasis1_chscr(xyzzyaa&
&aj13,1,is),gradr_eebasis1_chscr(1,1,is),gradr_enbasis1_chscr(1,1,is),&
&jas=jas,gjas=gjas)
elseif(xyzzyaaaw13)then
call xyzzyaaid1(ii,xyzzyaaaa13,eebasis1_chscr(xyzzyaaag13,1,is),eebasi&
&s1_chscr(xyzzyaaah13,1,is),nzeecut1_chscr(1,xyzzyaaad13,is),enbasis1_&
&chscr(xyzzyaaai13,1,is),enbasis1_chscr(xyzzyaaaj13,1,is),nzencut1_chs&
&cr(1,xyzzyaaaf13,is),eebasis_scr(xyzzyaaag13,1,1,xyzzyaaab13),eebasis&
&_scr(xyzzyaaah13,1,1,xyzzyaaab13),nzeecut_scr(1,1,xyzzyaaad13,xyzzyaa&
&ab13),enbasis_scr(xyzzyaaai13,1,1,xyzzyaaab13),enbasis_scr(xyzzyaaaj1&
&3,1,1,xyzzyaaab13),nzencut_scr(1,1,xyzzyaaaf13,xyzzyaaab13),jas=jas)
endif
elseif(.not.aderiv)then
if(xyzzyaaay13)then
call xyzzyaaid1(ii,xyzzyaaaa13,eebasis_scr(xyzzyaaag13,1,ii,is),eebasi&
&s_scr(xyzzyaaah13,1,ii,is),nzeecut_scr(1,ii,xyzzyaaad13,is),enbasis_s&
&cr(xyzzyaaai13,1,ii,is),enbasis_scr(xyzzyaaaj13,1,ii,is),nzencut_scr(&
&1,ii,xyzzyaaaf13,is),eebasis_scr(xyzzyaaag13,1,1,is),eebasis_scr(xyzz&
&yaaah13,1,1,is),nzeecut_scr(1,1,xyzzyaaad13,is),enbasis_scr(xyzzyaaai&
&13,1,1,is),enbasis_scr(xyzzyaaaj13,1,1,is),nzencut_scr(1,1,xyzzyaaaf1&
&3,is),grad_eebasis_scr(1,xyzzyaaag13,1,ii,is),grad_eebasis_scr(1,xyzz&
&yaaah13,1,ii,is),grad_enbasis_scr(1,xyzzyaaai13,1,ii,is),grad_enbasis&
&_scr(1,xyzzyaaaj13,1,ii,is),deebasis_scr(xyzzyaaag13,1,ii,is),deebasi&
&s_scr(xyzzyaaah13,1,ii,is),denbasis_scr(xyzzyaaai13,1,ii,is),denbasis&
&_scr(xyzzyaaaj13,1,ii,is),gradr_eebasis_scr(1,1,ii,is),gradr_enbasis_&
&scr(1,1,ii,is),lap_eebasis_scr(xyzzyaaag13,1,ii,is),lap_eebasis_scr(x&
&yzzyaaah13,1,ii,is),lap_enbasis_scr(xyzzyaaai13,1,ii,is),lap_enbasis_&
&scr(xyzzyaaaj13,1,ii,is),d2eebasis_scr(xyzzyaaag13,1,ii,is),d2eebasis&
&_scr(xyzzyaaah13,1,ii,is),d2enbasis_scr(xyzzyaaai13,1,ii,is),d2enbasi&
&s_scr(xyzzyaaaj13,1,ii,is),lapr_eebasis_scr(1,ii,is),lapr_enbasis_scr&
&(1,ii,is),jas=jas,gjas=gjas,ljas=ljas)
elseif(xyzzyaaax13)then
call xyzzyaaid1(ii,xyzzyaaaa13,eebasis_scr(xyzzyaaag13,1,ii,is),eebasi&
&s_scr(xyzzyaaah13,1,ii,is),nzeecut_scr(1,ii,xyzzyaaad13,is),enbasis_s&
&cr(xyzzyaaai13,1,ii,is),enbasis_scr(xyzzyaaaj13,1,ii,is),nzencut_scr(&
&1,ii,xyzzyaaaf13,is),eebasis_scr(xyzzyaaag13,1,1,is),eebasis_scr(xyzz&
&yaaah13,1,1,is),nzeecut_scr(1,1,xyzzyaaad13,is),enbasis_scr(xyzzyaaai&
&13,1,1,is),enbasis_scr(xyzzyaaaj13,1,1,is),nzencut_scr(1,1,xyzzyaaaf1&
&3,is),grad_eebasis_scr(1,xyzzyaaag13,1,ii,is),grad_eebasis_scr(1,xyzz&
&yaaah13,1,ii,is),grad_enbasis_scr(1,xyzzyaaai13,1,ii,is),grad_enbasis&
&_scr(1,xyzzyaaaj13,1,ii,is),deebasis_scr(xyzzyaaag13,1,ii,is),deebasi&
&s_scr(xyzzyaaah13,1,ii,is),denbasis_scr(xyzzyaaai13,1,ii,is),denbasis&
&_scr(xyzzyaaaj13,1,ii,is),gradr_eebasis_scr(1,1,ii,is),gradr_enbasis_&
&scr(1,1,ii,is),jas=jas,gjas=gjas)
elseif(xyzzyaaaw13)then
call xyzzyaaid1(ii,xyzzyaaaa13,eebasis_scr(xyzzyaaag13,1,ii,is),eebasi&
&s_scr(xyzzyaaah13,1,ii,is),nzeecut_scr(1,ii,xyzzyaaad13,is),enbasis_s&
&cr(xyzzyaaai13,1,ii,is),enbasis_scr(xyzzyaaaj13,1,ii,is),nzencut_scr(&
&1,ii,xyzzyaaaf13,is),eebasis_scr(xyzzyaaag13,1,1,is),eebasis_scr(xyzz&
&yaaah13,1,1,is),nzeecut_scr(1,1,xyzzyaaad13,is),enbasis_scr(xyzzyaaai&
&13,1,1,is),enbasis_scr(xyzzyaaaj13,1,1,is),nzencut_scr(1,1,xyzzyaaaf1&
&3,is),jas=jas)
endif
elseif(aderiv)then
call xyzzyaaid1(ii,xyzzyaaaa13,eebasis_scr(xyzzyaaag13,1,ii,is),eebasi&
&s_scr(xyzzyaaah13,1,ii,is),nzeecut_scr(1,ii,xyzzyaaad13,is),enbasis_s&
&cr(xyzzyaaai13,1,ii,is),enbasis_scr(xyzzyaaaj13,1,ii,is),nzencut_scr(&
&1,ii,xyzzyaaaf13,is),eebasis_scr(xyzzyaaag13,1,1,is),eebasis_scr(xyzz&
&yaaah13,1,1,is),nzeecut_scr(1,1,xyzzyaaad13,is),enbasis_scr(xyzzyaaai&
&13,1,1,is),enbasis_scr(xyzzyaaaj13,1,1,is),nzencut_scr(1,1,xyzzyaaaf1&
&3,is),grad_eebasis_scr(1,xyzzyaaag13,1,ii,is),grad_eebasis_scr(1,xyzz&
&yaaah13,1,ii,is),grad_enbasis_scr(1,xyzzyaaai13,1,ii,is),grad_enbasis&
&_scr(1,xyzzyaaaj13,1,ii,is),deebasis_scr(xyzzyaaag13,1,ii,is),deebasi&
&s_scr(xyzzyaaah13,1,ii,is),denbasis_scr(xyzzyaaai13,1,ii,is),denbasis&
&_scr(xyzzyaaaj13,1,ii,is),gradr_eebasis_scr(1,1,ii,is),gradr_enbasis_&
&scr(1,1,ii,is),lap_eebasis_scr(xyzzyaaag13,1,ii,is),lap_eebasis_scr(x&
&yzzyaaah13,1,ii,is),lap_enbasis_scr(xyzzyaaai13,1,ii,is),lap_enbasis_&
&scr(xyzzyaaaj13,1,ii,is),d2eebasis_scr(xyzzyaaag13,1,ii,is),d2eebasis&
&_scr(xyzzyaaah13,1,ii,is),d2enbasis_scr(xyzzyaaai13,1,ii,is),d2enbasi&
&s_scr(xyzzyaaaj13,1,ii,is),lapr_eebasis_scr(1,ii,is),lapr_enbasis_scr&
&(1,ii,is),jas=jas,gjas=gjas,ljas=ljas,djas=djas,dgjas=dgjas,dljas=dlj&
&as)
endif
if(aderiv)then
xyzzyaaak13=xyzzyaaaj1(xyzzyaaag1(xyzzyaaaa13)+1)
xyzzyaaal13=xyzzyaaak13+1
xyzzyaaam13=xyzzyaaak13+xyzzyaaax1(xyzzyaaaa13)
if(xyzzyaaaw13.and..not.xyzzyaaeq1(ii,is))then
xyzzyaaee1(xyzzyaaal13:xyzzyaaam13,ii,is)=djas(xyzzyaaal13:xyzzyaaam13&
&)*xyzzyaabj1(xyzzyaaaa13)
endif
if(xyzzyaaax13.and..not.xyzzyaaer1(ii,is))then
xyzzyaaef1(1:3,xyzzyaaal13:xyzzyaaam13,ii,is)=dgjas(1:3,xyzzyaaal13:xy&
&zzyaaam13)
endif
if(xyzzyaaay13.and..not.xyzzyaaes1(ii,is))then
xyzzyaaeg1(xyzzyaaal13:xyzzyaaam13,ii,is)=dljas(xyzzyaaal13:xyzzyaaam1&
&3)
endif
endif
if(detail)then
if(xyzzyaaet1(ii,xyzzyaaaa13,is))then
jas=xyzzyaaeh1(ii,xyzzyaaaa13,is)
elseif(xyzzyaaaw13)then
xyzzyaaeh1(ii,xyzzyaaaa13,is)=jas
xyzzyaaet1(ii,xyzzyaaaa13,is)=.true.
endif
if(xyzzyaaeu1(ii,xyzzyaaaa13,is))then
gjas=xyzzyaaei1(1:3,ii,xyzzyaaaa13,is)
ljas=xyzzyaaej1(ii,xyzzyaaaa13,is)
elseif(xyzzyaaax13.and.xyzzyaaay13)then
xyzzyaaei1(1:3,ii,xyzzyaaaa13,is)=gjas
xyzzyaaej1(ii,xyzzyaaaa13,is)=ljas
xyzzyaaeu1(ii,xyzzyaaaa13,is)=.true.
endif
elseif(xyzzyaaao1(xyzzyaaaa13)==0)then
if(xyzzyaaaw13)then
xyzzyaaar13=xyzzyaaar13+jas
else
jas=0.d0
endif
if(xyzzyaaax13)then
xyzzyaaas13=xyzzyaaas13+gjas
else
gjas=0.d0
endif
endif
if(need_val)then
xyzzyaaan13=xyzzyaaan13+jas
xyzzyaaao13=xyzzyaaao13+jas*xyzzyaabj1(xyzzyaaaa13)
endif
if(xyzzyaaat13)xyzzyaaap13=xyzzyaaap13+gjas
if(xyzzyaaau13)xyzzyaaaq13=xyzzyaaaq13+ljas
call timer('TERM_'//trim(i2s(xyzzyaaaa13)),.false.)
enddo
if(aderiv)then
if(need_val.and..not.xyzzyaaeq1(ii,is))xyzzyaaeq1(ii,is)=.true.
if(xyzzyaaat13.and..not.xyzzyaaer1(ii,is))xyzzyaaer1(ii,is)=.true.
if(xyzzyaaau13.and..not.xyzzyaaes1(ii,is))xyzzyaaes1(ii,is)=.true.
endif
if(.not.detail)then
if(need_val)then
if(xyzzyaaen1(ii,is))then
xyzzyaaan13=xyzzyaaan13+xyzzyaady1(ii,is)
xyzzyaaao13=xyzzyaaao13+xyzzyaady1(ii,is)
else
xyzzyaady1(ii,is)=xyzzyaaar13
xyzzyaaen1(ii,is)=.true.
endif
endif
if(xyzzyaaat13)then
if(xyzzyaaeo1(ii,is))then
xyzzyaaap13=xyzzyaaap13+xyzzyaaea1(:,ii,is)
else
xyzzyaaea1(:,ii,is)=xyzzyaaas13
xyzzyaaeo1(ii,is)=.true.
endif
endif
endif
if(need_val)then
xyzzyaaec1(ii,is)=xyzzyaaan13
xyzzyaaed1(ii,is)=xyzzyaaao13
xyzzyaaek1(ii,is)=.true.
endif
if(xyzzyaaat13)then
xyzzyaadz1(:,ii,is)=xyzzyaaap13
xyzzyaael1(ii,is)=.true.
endif
if(xyzzyaaau13)then
xyzzyaaeb1(ii,is)=xyzzyaaaq13
xyzzyaaem1(ii,is)=.true.
endif
call timer('GET_JASTROW1',.false.)
end subroutine xyzzyaahn1
subroutine xyzzyaaho1(is)
implicit none
integer,intent(in) :: is
integer xyzzyaaaa14
real(dp) xyzzyaaab14,xyzzyaaac14
if(xyzzyaaep1(is))return
if(all(xyzzyaaek1(:,is)))then
xyzzyaadx1(is)=sum(xyzzyaaed1(:,is))
else
call timer('GET_JASTROW_VAL',.true.)
xyzzyaaab14=0.d0
call get_gbasis(is,.false.,.false.)
do xyzzyaaaa14=1,xyzzyaaaa1
call timer('TERM_'//trim(i2s(xyzzyaaaa14)),.true.)
call xyzzyaaie1(xyzzyaaaa14,eebasis_scr(1,1,1,is),nzeecut_scr(1,1,1,is&
&),enbasis_scr(1,1,1,is),nzencut_scr(1,1,1,is),xyzzyaaac14)
xyzzyaaab14=xyzzyaaab14+xyzzyaaac14
call timer('TERM_'//trim(i2s(xyzzyaaaa14)),.false.)
enddo
xyzzyaadx1(is)=xyzzyaaab14
call timer('GET_JASTROW_VAL',.false.)
endif
xyzzyaaep1(is)=.true.
end subroutine xyzzyaaho1
subroutine gen_config_gjastrow(pt_config)
implicit none
type(config_wfn_gjastrow),pointer :: pt_config
integer xyzzyaaaa15
allocate(pt_config,stat=xyzzyaaaa15)
call check_alloc(xyzzyaaaa15,'GEN_CONFIG_WFN_JASTROW','container')
if(jasbuf)then
allocate(pt_config%pt_chi(netot),pt_config%pt_gradchi(3,netot),pt_conf&
&ig%pt_chi_valid(netot),pt_config%pt_gradchi_valid(netot),stat=xyzzyaa&
&aa15)
call check_alloc(xyzzyaaaa15,'GEN_CONFIG_WFN_JASTROW','jasbuf')
pt_config%pt_valjas=0.d0
pt_config%pt_chi(:)=0.d0
pt_config%pt_gradchi(:,:)=0.d0
pt_config%pt_valjas_valid=.false.
pt_config%pt_chi_valid(:)=.false.
pt_config%pt_gradchi_valid(:)=.false.
endif
end subroutine gen_config_gjastrow
subroutine delete_config_gjastrow(pt_config)
implicit none
type(config_wfn_gjastrow),pointer :: pt_config
if(jasbuf)deallocate(pt_config%pt_chi,pt_config%pt_gradchi,pt_config%p&
&t_chi_valid,pt_config%pt_gradchi_valid)
deallocate(pt_config)
end subroutine delete_config_gjastrow
subroutine copy_config_gjastrow(pt_from,pt_to)
implicit none
type(config_wfn_gjastrow),pointer :: pt_from,pt_to
if(jasbuf)then
pt_to%pt_valjas=pt_from%pt_valjas
call dcopy(netot,pt_from%pt_chi,1,pt_to%pt_chi,1)
call dcopy(three_netot,pt_from%pt_gradchi,1,pt_to%pt_gradchi,1)
pt_to%pt_valjas_valid=pt_from%pt_valjas_valid
pt_to%pt_chi_valid=pt_from%pt_chi_valid
pt_to%pt_gradchi_valid=pt_from%pt_gradchi_valid
endif
end subroutine copy_config_gjastrow
subroutine config_to_pt_gjastrow(pt_config,k)
use slaarnaaf, only : valjas_config
implicit none
integer,intent(in) :: k
type(config_wfn_gjastrow),pointer :: pt_config
if(jasbuf)then
pt_config%pt_valjas_valid=allocated(valjas_config)
if(pt_config%pt_valjas_valid)pt_config%pt_valjas=valjas_config(k)
endif
end subroutine config_to_pt_gjastrow
subroutine pt_to_config_gjastrow(pt_config)
use slaarnaaf, only : add_config
implicit none
type(config_wfn_gjastrow),pointer :: pt_config
if(jasbuf)call add_config(modify=.true.,valjas=pt_config%pt_valjas)
end subroutine pt_to_config_gjastrow
subroutine redist_allocations_gjastrow(kmax)
implicit none
integer,intent(in) :: kmax
integer xyzzyaaaa20
xyzzyaaew1=jasbuf
if(jasbuf)then
allocate(xyzzyaaev1(kmax),stat=xyzzyaaaa20)
call check_alloc(xyzzyaaaa20,'REDIST_ALLOCATIONS_GJASTROW','valjas_trf&
&')
endif
end subroutine redist_allocations_gjastrow
subroutine redist_load_gjastrow(pt_config,k,blocking)
implicit none
integer,intent(in) :: k
logical,intent(in) :: blocking
type(config_wfn_gjastrow),pointer :: pt_config
xyzzyaaew1=xyzzyaaew1.and.pt_config%pt_valjas_valid
if(xyzzyaaew1.or..not.blocking)xyzzyaaev1(k)=pt_config%pt_valjas
end subroutine redist_load_gjastrow
subroutine redist_send_gjastrow(jnode,kbase,k,nbt,reqbase,blocking)
use parallel
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa22
if(jasbuf)then
if(blocking)then
call mpi_ssend(xyzzyaaew1,1,mpi_logical,jnode,move_msg,mpi_comm_world,&
&ierror)
call checkmpi(ierror,'ssend valjas_valid in redist_send')
if(xyzzyaaew1)then
call mpi_ssend(xyzzyaaev1(kbase+1),k,mpi_double_precision,jnode,move_m&
&sg,mpi_comm_world,ierror)
call checkmpi(ierror,'ssend valjas in redist_send')
endif
else
nbt=nbt+1
call mpi_isend(xyzzyaaev1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa22,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa22
call checkmpi(ierror,'ssend valjas in redist_send')
endif
endif
end subroutine redist_send_gjastrow
subroutine redist_recv_gjastrow(jnode,kbase,k,nbt,reqbase,blocking)
use parallel
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa23
if(jasbuf)then
if(blocking)then
call mpi_recv(xyzzyaaew1,1,mpi_logical,jnode,move_msg,mpi_comm_world,s&
&tatus,ierror)
call checkmpi(ierror,'recv valjas_valid in redist_recv')
if(xyzzyaaew1)then
call mpi_recv(xyzzyaaev1(kbase+1),k,mpi_double_precision,jnode,move_ms&
&g,mpi_comm_world,status,ierror)
call checkmpi(ierror,'recv valjas in redist_recv')
endif
else
nbt=nbt+1
call mpi_irecv(xyzzyaaev1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa23,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa23
call checkmpi(ierror,'recv valjas in redist_recv')
endif
endif
end subroutine redist_recv_gjastrow
subroutine redist_save_gjastrow(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_gjastrow),pointer :: pt_config
pt_config%pt_valjas_valid=xyzzyaaew1
if(xyzzyaaew1)pt_config%pt_valjas=xyzzyaaev1(k)
end subroutine redist_save_gjastrow
subroutine redist_deallocations_gjastrow
implicit none
if(jasbuf)deallocate(xyzzyaaev1)
end subroutine redist_deallocations_gjastrow
subroutine load_from_pt_gjastrow(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_gjastrow),pointer :: pt_config
if(jasbuf)then
xyzzyaadx1(is)=pt_config%pt_valjas
call dcopy(netot,pt_config%pt_chi,1,xyzzyaady1(1,is),1)
call dcopy(three_netot,pt_config%pt_gradchi,1,xyzzyaaea1(1,1,is),1)
xyzzyaaep1(is)=pt_config%pt_valjas_valid
xyzzyaaen1(1:netot,is)=pt_config%pt_chi_valid(1:netot)
xyzzyaaeo1(1:netot,is)=pt_config%pt_gradchi_valid(1:netot)
endif
end subroutine load_from_pt_gjastrow
subroutine save_to_pt_gjastrow(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_gjastrow),pointer :: pt_config
if(jasbuf)then
pt_config%pt_valjas=xyzzyaadx1(is)
call dcopy(netot,xyzzyaady1(1,is),1,pt_config%pt_chi,1)
call dcopy(three_netot,xyzzyaaea1(1,1,is),1,pt_config%pt_gradchi,1)
pt_config%pt_valjas_valid=xyzzyaaep1(is)
pt_config%pt_chi_valid(1:netot)=xyzzyaaen1(1:netot,is)
pt_config%pt_gradchi_valid(1:netot)=xyzzyaaeo1(1:netot,is)
endif
end subroutine save_to_pt_gjastrow
subroutine add_config_gjastrow_items(is)
use slaarnaaf,only : add_config
implicit none
integer,intent(in) :: is
call xyzzyaaho1(is)
call add_config(modify=.true.,valjas=xyzzyaadx1(is))
end subroutine add_config_gjastrow_items
subroutine setup_storage_gjastrow(nconfig,ignore)
use slaarnaaf,only : valjas_config
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(:)
integer xyzzyaaaa29
allocate(xyzzyaaex1(nconfig),xyzzyaaey1(3,netot,nconfig),xyzzyaaez1(ne&
&tot,nconfig),xyzzyaafa1(netot,nconfig),xyzzyaafb1(netot,nconfig),stat&
&=xyzzyaaaa29)
call check_alloc(xyzzyaaaa29,'SETUP_STORAGE_JASTROW','*jas_store')
xyzzyaaex1=0.d0
xyzzyaaey1=0.d0
xyzzyaaez1=0.d0
xyzzyaafa1=0.d0
xyzzyaafb1=0.d0
allocate(xyzzyaafc1(nconfig),xyzzyaafd1(nconfig),xyzzyaafe1(nconfig),x&
&yzzyaaff1(nconfig),stat=xyzzyaaaa29)
call check_alloc(xyzzyaaaa29,'SETUP_STORAGE_JASTROW','valjas_svalid')
xyzzyaafc1=.false.
xyzzyaafd1=.false.
xyzzyaafe1=.false.
xyzzyaaff1=.false.
if(allocated(valjas_config))then
call dcopy(nconfig,valjas_config(1),1,xyzzyaaex1(1),1)
xyzzyaafc1(1:nconfig)=.true.
endif
end subroutine setup_storage_gjastrow
subroutine finish_storage_gjastrow
implicit none
deallocate(xyzzyaaex1,xyzzyaaey1,xyzzyaaez1,xyzzyaafa1,xyzzyaafb1)
deallocate(xyzzyaafc1,xyzzyaafd1,xyzzyaafe1,xyzzyaaff1)
end subroutine finish_storage_gjastrow
subroutine load_from_storage_gjastrow(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(xyzzyaafc1(icfg))then
xyzzyaadx1(is)=xyzzyaaex1(icfg)
xyzzyaaep1(is)=.true.
endif
if(xyzzyaaff1(icfg))then
call dcopy(netot,xyzzyaafa1(1,icfg),1,xyzzyaaed1(1,is),1)
call dcopy(netot,xyzzyaafb1(1,icfg),1,xyzzyaaec1(1,is),1)
xyzzyaaek1(:,is)=.true.
endif
if(xyzzyaafd1(icfg))then
call dcopy(three_netot,xyzzyaaey1(1,1,icfg),1,xyzzyaadz1(1,1,is),1)
xyzzyaael1(:,is)=.true.
endif
if(xyzzyaafe1(icfg))then
call dcopy(netot,xyzzyaaez1(1,icfg),1,xyzzyaaeb1(1,is),1)
xyzzyaaem1(:,is)=.true.
endif
end subroutine load_from_storage_gjastrow
subroutine save_to_storage_gjastrow(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(xyzzyaaep1(is).and..not.xyzzyaafc1(icfg))then
xyzzyaaex1(icfg)=xyzzyaadx1(is)
xyzzyaafc1(icfg)=.true.
endif
if(all(xyzzyaaek1(:,is)))then
call dcopy(netot,xyzzyaaed1(1,is),1,xyzzyaafa1(1,icfg),1)
call dcopy(netot,xyzzyaaec1(1,is),1,xyzzyaafb1(1,icfg),1)
xyzzyaaff1(icfg)=.true.
endif
if(all(xyzzyaael1(:,is)))then
call dcopy(three_netot,xyzzyaadz1(1,1,is),1,xyzzyaaey1(1,1,icfg),1)
xyzzyaafd1(icfg)=.true.
endif
if(all(xyzzyaaem1(:,is)))then
call dcopy(netot,xyzzyaaeb1(1,is),1,xyzzyaaez1(1,icfg),1)
xyzzyaafe1(icfg)=.true.
endif
end subroutine save_to_storage_gjastrow
subroutine clone_scratch_gjastrow(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa33
if(scr_tasks(ikinetic,js))then
if(all(xyzzyaael1(:,is)).and..not.all(xyzzyaael1(:,js)))then
call dcopy(three_netot,xyzzyaadz1(1,1,is),1,xyzzyaadz1(1,1,js),1)
xyzzyaael1(:,js)=.true.
endif
if(all(xyzzyaaem1(:,is)).and..not.all(xyzzyaaem1(:,js)))then
call dcopy(netot,xyzzyaaeb1(1,is),1,xyzzyaaeb1(1,js),1)
xyzzyaaem1(:,js)=.true.
endif
endif
if(scr_tasks(ikinetic_detail,js))then
do xyzzyaaaa33=1,xyzzyaaaa1
if(all(xyzzyaaeu1(:,xyzzyaaaa33,is)).and..not.all(xyzzyaaeu1(:,xyzzyaa&
&aa33,js)))then
call dcopy(three_netot,xyzzyaaei1(1,1,xyzzyaaaa33,is),1,xyzzyaaei1(1,1&
&,xyzzyaaaa33,js),1)
call dcopy(netot,xyzzyaaej1(1,xyzzyaaaa33,is),1,xyzzyaaej1(1,xyzzyaaaa&
&33,js),1)
xyzzyaaeu1(:,xyzzyaaaa33,js)=.true.
endif
enddo
endif
if(scr_tasks(iwfn_detail,js))then
if(xyzzyaaep1(is).and..not.xyzzyaaep1(js))then
xyzzyaadx1(js)=xyzzyaadx1(is)
xyzzyaaep1(js)=.true.
endif
if(all(xyzzyaaek1(:,is)).and..not.all(xyzzyaaek1(:,js)))then
call dcopy(netot,xyzzyaaec1(1,is),1,xyzzyaaec1(1,js),1)
call dcopy(netot,xyzzyaaed1(1,is),1,xyzzyaaed1(1,js),1)
xyzzyaaek1(:,js)=.true.
endif
if(all(xyzzyaaen1(:,is)).and..not.all(xyzzyaaen1(:,js)))then
call dcopy(netot,xyzzyaady1(1,is),1,xyzzyaady1(1,js),1)
xyzzyaaen1(:,js)=.true.
endif
do xyzzyaaaa33=1,xyzzyaaaa1
if(all(xyzzyaaet1(:,xyzzyaaaa33,is)).and..not.all(xyzzyaaet1(:,xyzzyaa&
&aa33,js)))then
call dcopy(netot,xyzzyaaeh1(1,xyzzyaaaa33,is),1,xyzzyaaeh1(1,xyzzyaaaa&
&33,js),1)
xyzzyaaet1(:,xyzzyaaaa33,js)=.true.
endif
enddo
endif
end subroutine clone_scratch_gjastrow
subroutine setup_gjastrow_params(nparam)
implicit none
integer,intent(inout) :: nparam
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34,xyzzyaaad34,xyzzyaaae34,xy&
&zzyaaaf34,xyzzyaaag34,xyzzyaaah34,xyzzyaaai34,xyzzyaaaj34
nparam=0
if(.not.opt_jastrow)return
allocate(xyzzyaafh1(xyzzyaaaa1),xyzzyaafg1(xyzzyaaaa1),xyzzyaafj1(xyzz&
&yaaaa1),xyzzyaafi1(xyzzyaaaa1),xyzzyaafk1(xyzzyaaaa1),xyzzyaafl1(xyzz&
&yaaaa1),stat=xyzzyaaaa34)
call check_alloc(xyzzyaaaa34,'SETUP_GJASTROW_PARAMS','nopt')
xyzzyaafh1=0
xyzzyaafg1=0
xyzzyaafj1=0
xyzzyaafi1=0
xyzzyaafk1=0
xyzzyaafl1=0
do xyzzyaaab34=1,xyzzyaaaa1
call count_eebasis_params(xyzzyaaaq1(xyzzyaaab34),xyzzyaafh1(xyzzyaaab&
&34))
call count_eebasis_params(xyzzyaaar1(xyzzyaaab34),xyzzyaafg1(xyzzyaaab&
&34))
call count_enbasis_params(xyzzyaaas1(xyzzyaaab34),xyzzyaafj1(xyzzyaaab&
&34))
call count_enbasis_params(xyzzyaaat1(xyzzyaaab34),xyzzyaafi1(xyzzyaaab&
&34))
xyzzyaaac34=xyzzyaaag1(xyzzyaaab34)
do xyzzyaaae34=1,xyzzyaaav1(xyzzyaaab34)
xyzzyaaad34=xyzzyaaaj1(xyzzyaaac34+xyzzyaaae34)
do xyzzyaaaf34=1,xyzzyaaaw1(xyzzyaaab34)
if(xyzzyaaci1(xyzzyaaad34+xyzzyaaaf34)==pflag_opt)xyzzyaafk1(xyzzyaaab&
&34)=xyzzyaafk1(xyzzyaaab34)+1
enddo
enddo
xyzzyaafl1(xyzzyaaab34)=xyzzyaafh1(xyzzyaaab34)+        xyzzyaafg1(xyz&
&zyaaab34)+xyzzyaafj1(xyzzyaaab34)+xyzzyaafi1(xyzzyaaab34)+xyzzyaafk1(&
&xyzzyaaab34)
enddo
nparam=sum(xyzzyaafl1)
allocate(xyzzyaafm1(nparam),xyzzyaafn1(nparam),xyzzyaafo1(nparam),xyzz&
&yaafp1(nparam),xyzzyaafq1(nparam),xyzzyaafr1(nparam),stat=xyzzyaaaa34&
&)
call check_alloc(xyzzyaaaa34,'SETUP_GJASTROW_PARAMS','optparam_term')
xyzzyaaai34=0
do xyzzyaaab34=1,xyzzyaaaa1
xyzzyaaag34=xyzzyaaai34+1
if(xyzzyaafh1(xyzzyaaab34)>0)then
xyzzyaaah34=xyzzyaaai34+1
xyzzyaaai34=xyzzyaaai34+xyzzyaafh1(xyzzyaaab34)
xyzzyaafn1(xyzzyaaah34:xyzzyaaai34)=1
xyzzyaafp1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaag34
xyzzyaafo1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaah34
xyzzyaafq1(xyzzyaaah34:xyzzyaaai34)=0
xyzzyaafr1(xyzzyaaah34:xyzzyaaai34)=0
endif
if(xyzzyaafg1(xyzzyaaab34)>0)then
xyzzyaaah34=xyzzyaaai34+1
xyzzyaaai34=xyzzyaaai34+xyzzyaafg1(xyzzyaaab34)
xyzzyaafn1(xyzzyaaah34:xyzzyaaai34)=2
xyzzyaafp1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaag34
xyzzyaafo1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaah34
xyzzyaafq1(xyzzyaaah34:xyzzyaaai34)=0
xyzzyaafr1(xyzzyaaah34:xyzzyaaai34)=0
endif
if(xyzzyaafj1(xyzzyaaab34)>0)then
xyzzyaaah34=xyzzyaaai34+1
xyzzyaaai34=xyzzyaaai34+xyzzyaafj1(xyzzyaaab34)
xyzzyaafn1(xyzzyaaah34:xyzzyaaai34)=3
xyzzyaafp1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaag34
xyzzyaafo1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaah34
xyzzyaafq1(xyzzyaaah34:xyzzyaaai34)=0
xyzzyaafr1(xyzzyaaah34:xyzzyaaai34)=0
endif
if(xyzzyaafi1(xyzzyaaab34)>0)then
xyzzyaaah34=xyzzyaaai34+1
xyzzyaaai34=xyzzyaaai34+xyzzyaafi1(xyzzyaaab34)
xyzzyaafn1(xyzzyaaah34:xyzzyaaai34)=4
xyzzyaafp1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaag34
xyzzyaafo1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaah34
xyzzyaafq1(xyzzyaaah34:xyzzyaaai34)=0
xyzzyaafr1(xyzzyaaah34:xyzzyaaai34)=0
endif
if(xyzzyaafk1(xyzzyaaab34)>0)then
xyzzyaaah34=xyzzyaaai34+1
xyzzyaaai34=xyzzyaaai34+xyzzyaafk1(xyzzyaaab34)
xyzzyaafn1(xyzzyaaah34:xyzzyaaai34)=5
xyzzyaafp1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaag34
xyzzyaafo1(xyzzyaaah34:xyzzyaaai34)=xyzzyaaah34
xyzzyaaaj34=xyzzyaaah34-1
xyzzyaaac34=xyzzyaaag1(xyzzyaaab34)
do xyzzyaaae34=1,xyzzyaaav1(xyzzyaaab34)
xyzzyaaad34=xyzzyaaaj1(xyzzyaaac34+xyzzyaaae34)
do xyzzyaaaf34=1,xyzzyaaaw1(xyzzyaaab34)
if(xyzzyaaci1(xyzzyaaad34+xyzzyaaaf34)==pflag_opt)then
xyzzyaaaj34=xyzzyaaaj34+1
xyzzyaafq1(xyzzyaaaj34)=xyzzyaaae34
xyzzyaafr1(xyzzyaaaj34)=xyzzyaaaf34
endif
enddo
enddo
endif
xyzzyaafm1(xyzzyaaag34:xyzzyaaai34)=xyzzyaaab34
enddo
call xyzzyaahp1
end subroutine setup_gjastrow_params
subroutine finish_gjastrow_params
implicit none
if(.not.opt_jastrow)return
deallocate(xyzzyaafh1,xyzzyaafg1,xyzzyaafj1,xyzzyaafi1,xyzzyaafk1,xyzz&
&yaafl1)
deallocate(xyzzyaafm1,xyzzyaafn1,xyzzyaafo1,xyzzyaafp1,xyzzyaafq1,xyzz&
&yaafr1)
call xyzzyaahq1
end subroutine finish_gjastrow_params
subroutine get_gjastrow_params(params,has_lolim,lolim,has_hilim,hilim,&
&is_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,&
&label)
implicit none
real(dp),intent(inout) :: params(:),lolim(:),hilim(:)
logical,intent(inout) :: has_lolim(:),has_hilim(:),is_shallow(:),is_re&
&dundant(:),is_linear(:),is_loglinear(:),has_aderiv(:),affect_map(:,:)
character(2),intent(inout) :: label(:)
integer xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36,xyzzyaaad36,xyzzyaaae36,xy&
&zzyaaaf36,xyzzyaaag36,xyzzyaaah36,xyzzyaaai36,xyzzyaaaj36
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
do xyzzyaaaf36=1,size(affect_map,1)
affect_map(xyzzyaaaf36,xyzzyaaaf36)=.true.
enddo
xyzzyaaad36=0
do xyzzyaaae36=1,xyzzyaaaa1
xyzzyaaaa36=xyzzyaaad36+1
if(xyzzyaafh1(xyzzyaaae36)>0)then
xyzzyaaac36=xyzzyaaad36+1
xyzzyaaad36=xyzzyaaad36+xyzzyaafh1(xyzzyaaae36)
call get_eebasis_params(xyzzyaaaq1(xyzzyaaae36),params(xyzzyaaac36:xyz&
&zyaaad36),has_lolim(xyzzyaaac36:xyzzyaaad36),lolim(xyzzyaaac36:xyzzya&
&aad36),has_hilim(xyzzyaaac36:xyzzyaaad36),hilim(xyzzyaaac36:xyzzyaaad&
&36),is_shallow(xyzzyaaac36:xyzzyaaad36),is_redundant(xyzzyaaac36:xyzz&
&yaaad36),is_linear(xyzzyaaac36:xyzzyaaad36),is_loglinear(xyzzyaaac36:&
&xyzzyaaad36),has_aderiv(xyzzyaaac36:xyzzyaaad36),affect_map(xyzzyaaac&
&36:xyzzyaaad36,xyzzyaaac36:xyzzyaaad36),label(xyzzyaaac36:xyzzyaaad36&
&))
endif
if(xyzzyaafg1(xyzzyaaae36)>0)then
xyzzyaaac36=xyzzyaaad36+1
xyzzyaaad36=xyzzyaaad36+xyzzyaafg1(xyzzyaaae36)
call get_eebasis_params(xyzzyaaar1(xyzzyaaae36),params(xyzzyaaac36:xyz&
&zyaaad36),has_lolim(xyzzyaaac36:xyzzyaaad36),lolim(xyzzyaaac36:xyzzya&
&aad36),has_hilim(xyzzyaaac36:xyzzyaaad36),hilim(xyzzyaaac36:xyzzyaaad&
&36),is_shallow(xyzzyaaac36:xyzzyaaad36),is_redundant(xyzzyaaac36:xyzz&
&yaaad36),is_linear(xyzzyaaac36:xyzzyaaad36),is_loglinear(xyzzyaaac36:&
&xyzzyaaad36),has_aderiv(xyzzyaaac36:xyzzyaaad36),affect_map(xyzzyaaac&
&36:xyzzyaaad36,xyzzyaaac36:xyzzyaaad36),label(xyzzyaaac36:xyzzyaaad36&
&))
endif
if(xyzzyaafj1(xyzzyaaae36)>0)then
xyzzyaaac36=xyzzyaaad36+1
xyzzyaaad36=xyzzyaaad36+xyzzyaafj1(xyzzyaaae36)
call get_enbasis_params(xyzzyaaas1(xyzzyaaae36),params(xyzzyaaac36:xyz&
&zyaaad36),has_lolim(xyzzyaaac36:xyzzyaaad36),lolim(xyzzyaaac36:xyzzya&
&aad36),has_hilim(xyzzyaaac36:xyzzyaaad36),hilim(xyzzyaaac36:xyzzyaaad&
&36),is_shallow(xyzzyaaac36:xyzzyaaad36),is_redundant(xyzzyaaac36:xyzz&
&yaaad36),is_linear(xyzzyaaac36:xyzzyaaad36),is_loglinear(xyzzyaaac36:&
&xyzzyaaad36),has_aderiv(xyzzyaaac36:xyzzyaaad36),affect_map(xyzzyaaac&
&36:xyzzyaaad36,xyzzyaaac36:xyzzyaaad36),label(xyzzyaaac36:xyzzyaaad36&
&))
endif
if(xyzzyaafi1(xyzzyaaae36)>0)then
xyzzyaaac36=xyzzyaaad36+1
xyzzyaaad36=xyzzyaaad36+xyzzyaafi1(xyzzyaaae36)
call get_enbasis_params(xyzzyaaat1(xyzzyaaae36),params(xyzzyaaac36:xyz&
&zyaaad36),has_lolim(xyzzyaaac36:xyzzyaaad36),lolim(xyzzyaaac36:xyzzya&
&aad36),has_hilim(xyzzyaaac36:xyzzyaaad36),hilim(xyzzyaaac36:xyzzyaaad&
&36),is_shallow(xyzzyaaac36:xyzzyaaad36),is_redundant(xyzzyaaac36:xyzz&
&yaaad36),is_linear(xyzzyaaac36:xyzzyaaad36),is_loglinear(xyzzyaaac36:&
&xyzzyaaad36),has_aderiv(xyzzyaaac36:xyzzyaaad36),affect_map(xyzzyaaac&
&36:xyzzyaaad36,xyzzyaaac36:xyzzyaaad36),label(xyzzyaaac36:xyzzyaaad36&
&))
endif
xyzzyaaab36=xyzzyaaad36
if(xyzzyaafk1(xyzzyaaae36)>0)then
xyzzyaaaj36=xyzzyaaad36
xyzzyaaac36=xyzzyaaad36+1
xyzzyaaad36=xyzzyaaad36+xyzzyaafk1(xyzzyaaae36)
xyzzyaaag36=xyzzyaaag1(xyzzyaaae36)
do xyzzyaaai36=1,xyzzyaaav1(xyzzyaaae36)
xyzzyaaah36=xyzzyaaaj1(xyzzyaaag36+xyzzyaaai36)
do xyzzyaaaf36=1,xyzzyaaaw1(xyzzyaaae36)
if(xyzzyaaci1(xyzzyaaah36+xyzzyaaaf36)==pflag_opt)then
xyzzyaaaj36=xyzzyaaaj36+1
params(xyzzyaaaj36)=param(xyzzyaaah36+xyzzyaaaf36)
has_lolim(xyzzyaaaj36)=xyzzyaacq1(xyzzyaaah36+xyzzyaaaf36)
lolim(xyzzyaaaj36)=xyzzyaaco1(xyzzyaaah36+xyzzyaaaf36)
has_hilim(xyzzyaaaj36)=xyzzyaacr1(xyzzyaaah36+xyzzyaaaf36)
hilim(xyzzyaaaj36)=xyzzyaacp1(xyzzyaaah36+xyzzyaaaf36)
is_shallow(xyzzyaaaj36)=xyzzyaacs1(xyzzyaaah36+xyzzyaaaf36)
is_loglinear(xyzzyaaaj36)=.true.
if(.not.xyzzyaahl1)has_aderiv(xyzzyaaaj36)=.true.
endif
enddo
enddo
if(xyzzyaaaa36<=xyzzyaaab36.and.xyzzyaaac36<=xyzzyaaad36)then
affect_map(xyzzyaaac36:xyzzyaaad36,xyzzyaaaa36:xyzzyaaab36)=.true.
if(all(param(xyzzyaaaj1(xyzzyaaag36+1)+1:xyzzyaaaj1(xyzzyaaag36+1)+xyz&
&zyaaax1(xyzzyaaae36))==0.d0))is_redundant(xyzzyaaaa36:xyzzyaaab36)=.t&
&rue.
endif
endif
if(xyzzyaaaa36<=xyzzyaaad36)label(xyzzyaaaa36:xyzzyaaad36)='J'//trim(i&
&2s(xyzzyaaae36))
enddo
end subroutine get_gjastrow_params
subroutine put_gjastrow_params(params,ignore,iparam_buffer,prestore,ba&
&d_params,restore_hint)
implicit none
integer,intent(in) :: iparam_buffer,restore_hint
logical,intent(in) :: ignore(:),prestore
real(dp),intent(inout) :: params(:)
logical,intent(out) :: bad_params
integer xyzzyaaaa37,xyzzyaaab37,xyzzyaaac37,xyzzyaaad37,xyzzyaaae37,xy&
&zzyaaaf37,xyzzyaaag37,xyzzyaaah37,xyzzyaaai37,xyzzyaaaj37,xyzzyaaak37&
&,xyzzyaaal37,xyzzyaaam37,xyzzyaaan37,xyzzyaaao37,xyzzyaaap37,xyzzyaaa&
&q37,xyzzyaaar37,xyzzyaaas37,xyzzyaaat37
logical xyzzyaaau37
bad_params=.false.
if(prestore)then
call xyzzyaahs1(iparam_buffer,restore_hint)
return
endif
if(iparam_buffer==0)then
xyzzyaaab37=0
do xyzzyaaac37=1,xyzzyaaaa1
xyzzyaaad37=xyzzyaaag1(xyzzyaaac37)
xyzzyaaat37=xyzzyaaai1(xyzzyaaac37)
if(xyzzyaafh1(xyzzyaaac37)>0)then
xyzzyaaaa37=xyzzyaaab37+1
xyzzyaaab37=xyzzyaaab37+xyzzyaafh1(xyzzyaaac37)
call put_eebasis_params(xyzzyaaaq1(xyzzyaaac37),params(xyzzyaaaa37:xyz&
&zyaaab37),ignore(xyzzyaaaa37:xyzzyaaab37))
endif
if(xyzzyaafg1(xyzzyaaac37)>0)then
xyzzyaaaa37=xyzzyaaab37+1
xyzzyaaab37=xyzzyaaab37+xyzzyaafg1(xyzzyaaac37)
call put_eebasis_params(xyzzyaaar1(xyzzyaaac37),params(xyzzyaaaa37:xyz&
&zyaaab37),ignore(xyzzyaaaa37:xyzzyaaab37))
endif
if(xyzzyaafj1(xyzzyaaac37)>0)then
xyzzyaaaa37=xyzzyaaab37+1
xyzzyaaab37=xyzzyaaab37+xyzzyaafj1(xyzzyaaac37)
call put_enbasis_params(xyzzyaaas1(xyzzyaaac37),params(xyzzyaaaa37:xyz&
&zyaaab37),ignore(xyzzyaaaa37:xyzzyaaab37))
endif
if(xyzzyaafi1(xyzzyaaac37)>0)then
xyzzyaaaa37=xyzzyaaab37+1
xyzzyaaab37=xyzzyaaab37+xyzzyaafi1(xyzzyaaac37)
call put_enbasis_params(xyzzyaaat1(xyzzyaaac37),params(xyzzyaaaa37:xyz&
&zyaaab37),ignore(xyzzyaaaa37:xyzzyaaab37))
endif
if(xyzzyaafk1(xyzzyaaac37)>0)then
xyzzyaaai37=xyzzyaaab37
xyzzyaaaa37=xyzzyaaab37+1
xyzzyaaab37=xyzzyaaab37+xyzzyaafk1(xyzzyaaac37)
do xyzzyaaae37=1,xyzzyaaav1(xyzzyaaac37)
xyzzyaaaf37=xyzzyaaaj1(xyzzyaaad37+xyzzyaaae37)
do xyzzyaaah37=1,xyzzyaaaw1(xyzzyaaac37)
if(xyzzyaaci1(xyzzyaaaf37+xyzzyaaah37)==pflag_opt)then
xyzzyaaai37=xyzzyaaai37+1
if(.not.ignore(xyzzyaaai37))param(xyzzyaaaf37+xyzzyaaah37)=params(xyzz&
&yaaai37)
endif
enddo
enddo
endif
if((xyzzyaafh1(xyzzyaaac37)>0.or.xyzzyaafg1(xyzzyaaac37)>0.or.xyzzyaaf&
&j1(xyzzyaaac37)>0.or.xyzzyaafi1(xyzzyaaac37)>0))then
do xyzzyaaae37=1,xyzzyaaav1(xyzzyaaac37)
if(xyzzyaacy1(xyzzyaaad37+xyzzyaaae37)==0)cycle
xyzzyaaap37=xyzzyaacz1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaaq37=xyzzyaadg1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaar37=xyzzyaadk1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaaj37=xyzzyaacx1(xyzzyaaad37+xyzzyaaae37)
if(xyzzyaaaj37>0)then
xyzzyaaak37=xyzzyaadb1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaal37=xyzzyaade1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaam37=xyzzyaadi1(xyzzyaaad37+xyzzyaaae37)
call dcopy(xyzzyaaak37,xyzzyaads1(xyzzyaaam37),1,xyzzyaadr1(xyzzyaaar3&
&7),1)
call dcopy(xyzzyaaaj37,xyzzyaadt1(xyzzyaaal37),1,eqn_rhs(xyzzyaaaq37),&
&1)
endif
xyzzyaaan37=xyzzyaadd1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaao37=xyzzyaadh1(xyzzyaaad37+xyzzyaaae37)
xyzzyaadr1(xyzzyaaao37:xyzzyaaao37+xyzzyaacy1(xyzzyaaad37+xyzzyaaae37)&
&*xyzzyaaaw1(xyzzyaaac37)-1)=0.d0
eqn_rhs(xyzzyaaan37:xyzzyaaan37+xyzzyaacy1(xyzzyaaad37+xyzzyaaae37)-1)&
&=0.d0
xyzzyaaag37=xyzzyaaaj1(xyzzyaaad37+xyzzyaaae37)+1
xyzzyaaas37=xyzzyaado1(xyzzyaaad37+xyzzyaaae37)
call xyzzyaaia1(xyzzyaaao1(xyzzyaaac37),xyzzyaaap1(xyzzyaaac37),xyzzya&
&aau1(xyzzyaaac37),xyzzyaaam1(xyzzyaaac37),xyzzyaaan1(xyzzyaaac37),xyz&
&zyaaaw1(xyzzyaaac37),xyzzyaacg1(xyzzyaaat37),xyzzyaacf1(1,xyzzyaaad37&
&+xyzzyaaae37),xyzzyaaaq1(xyzzyaaac37),xyzzyaaar1(xyzzyaaac37),xyzzyaa&
&as1(xyzzyaaac37),xyzzyaaat1(xyzzyaaac37),xyzzyaabf1(xyzzyaaac37),xyzz&
&yaabh1(xyzzyaaac37),xyzzyaabg1(xyzzyaaac37),xyzzyaabi1(xyzzyaaac37),x&
&yzzyaabr1(xyzzyaaac37),xyzzyaabs1(xyzzyaaac37),xyzzyaadn1(xyzzyaaas37&
&),xyzzyaact1(xyzzyaaag37),.false.,xyzzyaaca1(1,xyzzyaaad37+xyzzyaaae3&
&7),xyzzyaacb1(1,xyzzyaaad37+xyzzyaaae37),xyzzyaacc1(1,xyzzyaaad37+xyz&
&zyaaae37),xyzzyaacd1(1,xyzzyaaad37+xyzzyaaae37),xyzzyaadr1(xyzzyaaao3&
&7),eqn_rhs(xyzzyaaan37))
if(xyzzyaahg1)then
if(am_master)then
call wout(repeat('*',78))
call wout('DEBUG_OPT_CONSTRAINTS - before Gaussian elimination:')
call wout(repeat('*',78))
call wout('Term '//trim(i2s(xyzzyaaac37))//', channel '//trim(i2s(xyzz&
&yaaae37))//':')
call print_eqns(xyzzyaaaw1(xyzzyaaac37),xyzzyaaap37,xyzzyaadr1(xyzzyaa&
&ar37),eqn_rhs(xyzzyaaaq37),xyzzyaach1(xyzzyaaah1(xyzzyaaac37)+1:xyzzy&
&aaah1(xyzzyaaac37)+xyzzyaaaw1(xyzzyaaac37)))
call wout(repeat('*',78))
endif
endif
call redo_gaussian_elimination(xyzzyaaaw1(xyzzyaaac37),xyzzyaaap37,xyz&
&zyaadr1(xyzzyaaar37),eqn_rhs(xyzzyaaaq37),xyzzyaada1(xyzzyaaad37+xyzz&
&yaaae37),xyzzyaadp1(xyzzyaaaq37),xyzzyaadq1(xyzzyaaaq37),bad_params=x&
&yzzyaaau37)
bad_params=bad_params.or.xyzzyaaau37
if(xyzzyaahg1)then
if(am_master)then
call wout(repeat('*',78))
call wout('DEBUG_OPT_CONSTRAINTS - after Gaussian elimination:')
call wout(repeat('*',78))
call wout('Term '//trim(i2s(xyzzyaaac37))//', channel '//trim(i2s(xyzz&
&yaaae37))//':')
call print_eqns(xyzzyaaaw1(xyzzyaaac37),xyzzyaaap37,xyzzyaadr1(xyzzyaa&
&ar37),eqn_rhs(xyzzyaaaq37),xyzzyaach1(xyzzyaaah1(xyzzyaaac37)+1:xyzzy&
&aaah1(xyzzyaaac37)+xyzzyaaaw1(xyzzyaaac37)),xyzzyaada1(xyzzyaaad37+xy&
&zzyaaae37),xyzzyaadq1(xyzzyaaaq37),xyzzyaadp1(xyzzyaaaq37))
call wout(repeat('*',78))
endif
endif
enddo
endif
do xyzzyaaae37=1,xyzzyaaav1(xyzzyaaac37)
xyzzyaaaf37=xyzzyaaaj1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaag37=xyzzyaaaf37+1
xyzzyaaap37=xyzzyaacz1(xyzzyaaad37+xyzzyaaae37)
if(xyzzyaaap37>0)then
xyzzyaaaq37=xyzzyaadg1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaar37=xyzzyaadk1(xyzzyaaad37+xyzzyaaae37)
call backwards_substitution(xyzzyaaaw1(xyzzyaaac37),xyzzyaaap37,xyzzya&
&ada1(xyzzyaaad37+xyzzyaaae37),xyzzyaadr1(xyzzyaaar37),eqn_rhs(xyzzyaa&
&aq37),xyzzyaadp1(xyzzyaaaq37),xyzzyaadq1(xyzzyaaaq37),xyzzyaacu1(xyzz&
&yaaag37),param(xyzzyaaag37))
if(xyzzyaahh1)then
if(am_master)then
call wout(repeat('*',78))
call wout('DEBUG_OPT_PARAMETERS - visible parameters:')
call wout(repeat('*',78))
call wout('Term '//trim(i2s(xyzzyaaac37))//', channel '//trim(i2s(xyzz&
&yaaae37))//':')
do xyzzyaaah37=1,xyzzyaaaw1(xyzzyaaac37)
if(xyzzyaaci1(xyzzyaaaf37+xyzzyaaah37)/=pflag_det)then
call wout(trim(xyzzyaach1(xyzzyaaah1(xyzzyaaac37)+xyzzyaaah37))//' = '&
&,param(xyzzyaaaf37+xyzzyaaah37))
endif
enddo
call wout(repeat('*',78))
call wout(repeat('*',78))
call wout('DEBUG_OPT_PARAMETERS - determined parameters:')
call wout(repeat('*',78))
call wout('Term '//trim(i2s(xyzzyaaac37))//', channel '//trim(i2s(xyzz&
&yaaae37))//':')
do xyzzyaaah37=1,xyzzyaaaw1(xyzzyaaac37)
if(xyzzyaaci1(xyzzyaaaf37+xyzzyaaah37)==pflag_det)then
call wout(trim(xyzzyaach1(xyzzyaaah1(xyzzyaaac37)+xyzzyaaah37))//' = '&
&,param(xyzzyaaaf37+xyzzyaaah37))
endif
enddo
call wout(repeat('*',78))
endif
endif
endif
enddo
enddo
else
xyzzyaaac37=xyzzyaafm1(iparam_buffer)
xyzzyaaad37=xyzzyaaag1(xyzzyaaac37)
xyzzyaaaa37=xyzzyaafo1(iparam_buffer)
xyzzyaaat37=xyzzyaaai1(xyzzyaaac37)
if(xyzzyaafn1(iparam_buffer)/=5)then
select case(xyzzyaafn1(iparam_buffer))
case(1)
xyzzyaaab37=xyzzyaaaa37+xyzzyaafh1(xyzzyaaac37)-1
call put_eebasis_params(xyzzyaaaq1(xyzzyaaac37),params(xyzzyaaaa37:xyz&
&zyaaab37),ignore(xyzzyaaaa37:xyzzyaaab37))
case(2)
xyzzyaaab37=xyzzyaaaa37+xyzzyaafg1(xyzzyaaac37)-1
call put_eebasis_params(xyzzyaaar1(xyzzyaaac37),params(xyzzyaaaa37:xyz&
&zyaaab37),ignore(xyzzyaaaa37:xyzzyaaab37))
case(3)
xyzzyaaab37=xyzzyaaaa37+xyzzyaafj1(xyzzyaaac37)-1
call put_enbasis_params(xyzzyaaas1(xyzzyaaac37),params(xyzzyaaaa37:xyz&
&zyaaab37),ignore(xyzzyaaaa37:xyzzyaaab37))
case(4)
xyzzyaaab37=xyzzyaaaa37+xyzzyaafi1(xyzzyaaac37)-1
call put_enbasis_params(xyzzyaaat1(xyzzyaaac37),params(xyzzyaaaa37:xyz&
&zyaaab37),ignore(xyzzyaaaa37:xyzzyaaab37))
end select
do xyzzyaaae37=1,xyzzyaaav1(xyzzyaaac37)
if(xyzzyaacy1(xyzzyaaad37+xyzzyaaae37)==0)cycle
xyzzyaaap37=xyzzyaacz1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaaq37=xyzzyaadg1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaar37=xyzzyaadk1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaaj37=xyzzyaacx1(xyzzyaaad37+xyzzyaaae37)
if(xyzzyaaaj37>0)then
xyzzyaaak37=xyzzyaadb1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaal37=xyzzyaade1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaam37=xyzzyaadi1(xyzzyaaad37+xyzzyaaae37)
call dcopy(xyzzyaaak37,xyzzyaads1(xyzzyaaam37),1,xyzzyaadr1(xyzzyaaar3&
&7),1)
call dcopy(xyzzyaaaj37,xyzzyaadt1(xyzzyaaal37),1,eqn_rhs(xyzzyaaaq37),&
&1)
endif
xyzzyaaan37=xyzzyaadd1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaao37=xyzzyaadh1(xyzzyaaad37+xyzzyaaae37)
xyzzyaadr1(xyzzyaaao37:xyzzyaaao37+xyzzyaacy1(xyzzyaaad37+xyzzyaaae37)&
&*xyzzyaaaw1(xyzzyaaac37)-1)=0.d0
eqn_rhs(xyzzyaaan37:xyzzyaaan37+xyzzyaacy1(xyzzyaaad37+xyzzyaaae37)-1)&
&=0.d0
xyzzyaaag37=xyzzyaaaj1(xyzzyaaad37+xyzzyaaae37)+1
xyzzyaaas37=xyzzyaado1(xyzzyaaad37+xyzzyaaae37)
call xyzzyaaia1(xyzzyaaao1(xyzzyaaac37),xyzzyaaap1(xyzzyaaac37),xyzzya&
&aau1(xyzzyaaac37),xyzzyaaam1(xyzzyaaac37),xyzzyaaan1(xyzzyaaac37),xyz&
&zyaaaw1(xyzzyaaac37),xyzzyaacg1(xyzzyaaat37),xyzzyaacf1(1,xyzzyaaad37&
&+xyzzyaaae37),xyzzyaaaq1(xyzzyaaac37),xyzzyaaar1(xyzzyaaac37),xyzzyaa&
&as1(xyzzyaaac37),xyzzyaaat1(xyzzyaaac37),xyzzyaabf1(xyzzyaaac37),xyzz&
&yaabh1(xyzzyaaac37),xyzzyaabg1(xyzzyaaac37),xyzzyaabi1(xyzzyaaac37),x&
&yzzyaabr1(xyzzyaaac37),xyzzyaabs1(xyzzyaaac37),xyzzyaadn1(xyzzyaaas37&
&),xyzzyaact1(xyzzyaaag37),.false.,xyzzyaaca1(1,xyzzyaaad37+xyzzyaaae3&
&7),xyzzyaacb1(1,xyzzyaaad37+xyzzyaaae37),xyzzyaacc1(1,xyzzyaaad37+xyz&
&zyaaae37),xyzzyaacd1(1,xyzzyaaad37+xyzzyaaae37),xyzzyaadr1(xyzzyaaao3&
&7),eqn_rhs(xyzzyaaan37))
call redo_gaussian_elimination(xyzzyaaaw1(xyzzyaaac37),xyzzyaaap37,xyz&
&zyaadr1(xyzzyaaar37),eqn_rhs(xyzzyaaaq37),xyzzyaada1(xyzzyaaad37+xyzz&
&yaaae37),xyzzyaadp1(xyzzyaaaq37),xyzzyaadq1(xyzzyaaaq37),bad_params=x&
&yzzyaaau37)
bad_params=bad_params.or.xyzzyaaau37
enddo
do xyzzyaaae37=1,xyzzyaaav1(xyzzyaaac37)
xyzzyaaaf37=xyzzyaaaj1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaag37=xyzzyaaaf37+1
xyzzyaaap37=xyzzyaacz1(xyzzyaaad37+xyzzyaaae37)
if(xyzzyaaap37>0)then
xyzzyaaaq37=xyzzyaadg1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaar37=xyzzyaadk1(xyzzyaaad37+xyzzyaaae37)
call backwards_substitution(xyzzyaaaw1(xyzzyaaac37),xyzzyaaap37,xyzzya&
&ada1(xyzzyaaad37+xyzzyaaae37),xyzzyaadr1(xyzzyaaar37),eqn_rhs(xyzzyaa&
&aq37),xyzzyaadp1(xyzzyaaaq37),xyzzyaadq1(xyzzyaaaq37),xyzzyaacu1(xyzz&
&yaaag37),param(xyzzyaaag37))
endif
enddo
else
xyzzyaaae37=xyzzyaafq1(iparam_buffer)
xyzzyaaah37=xyzzyaafr1(iparam_buffer)
xyzzyaaaf37=xyzzyaaaj1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaag37=xyzzyaaaf37+1
param(xyzzyaaaf37+xyzzyaaah37)=params(iparam_buffer)
xyzzyaaap37=xyzzyaacz1(xyzzyaaad37+xyzzyaaae37)
if(xyzzyaaap37>0)then
xyzzyaaaq37=xyzzyaadg1(xyzzyaaad37+xyzzyaaae37)
xyzzyaaar37=xyzzyaadk1(xyzzyaaad37+xyzzyaaae37)
call backwards_substitution(xyzzyaaaw1(xyzzyaaac37),xyzzyaaap37,xyzzya&
&ada1(xyzzyaaad37+xyzzyaaae37),xyzzyaadr1(xyzzyaaar37),eqn_rhs(xyzzyaa&
&aq37),xyzzyaadp1(xyzzyaaaq37),xyzzyaadq1(xyzzyaaaq37),xyzzyaacu1(xyzz&
&yaaag37),param(xyzzyaaag37))
endif
endif
endif
call xyzzyaahr1(iparam_buffer)
end subroutine put_gjastrow_params
subroutine xyzzyaahp1
implicit none
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38
allocate(xyzzyaafs1(xyzzyaaaa1),stat=xyzzyaaaa38)
call check_alloc(xyzzyaaaa38,'SETUP_JASTROW_PBUFFER','term_ipbuffer0')
xyzzyaafs1=0
xyzzyaaac38=0
do xyzzyaaab38=1,xyzzyaaaa1
xyzzyaafs1(xyzzyaaab38)=xyzzyaaac38
xyzzyaaac38=xyzzyaaac38+xyzzyaaax1(xyzzyaaab38)*(1+xyzzyaafl1(xyzzyaaa&
&b38))
enddo
allocate(xyzzyaaft1(xyzzyaaac38),xyzzyaafu1(xyzzyaacw1),xyzzyaafv1(xyz&
&zyaacv1),stat=xyzzyaaaa38)
call check_alloc(xyzzyaaaa38,'SETUP_JASTROW_PBUFFER','param_pbuffer')
xyzzyaaft1=0.d0
xyzzyaafu1=0.d0
xyzzyaafv1=0.d0
end subroutine xyzzyaahp1
subroutine xyzzyaahq1
implicit none
deallocate(xyzzyaafs1)
deallocate(xyzzyaaft1,xyzzyaafu1,xyzzyaafv1)
end subroutine xyzzyaahq1
subroutine xyzzyaahr1(iparam)
implicit none
integer,intent(in) :: iparam
integer xyzzyaaaa40,xyzzyaaab40,xyzzyaaac40,xyzzyaaad40,xyzzyaaae40,xy&
&zzyaaaf40,xyzzyaaag40,xyzzyaaah40,xyzzyaaai40
if(iparam==0)then
do xyzzyaaaa40=1,xyzzyaaaa1
if(xyzzyaafh1(xyzzyaaaa40)>0)call save_eebasis_pbuffer(xyzzyaaaq1(xyzz&
&yaaaa40),0)
if(xyzzyaafg1(xyzzyaaaa40)>0)call save_eebasis_pbuffer(xyzzyaaar1(xyzz&
&yaaaa40),0)
if(xyzzyaafj1(xyzzyaaaa40)>0)call save_enbasis_pbuffer(xyzzyaaas1(xyzz&
&yaaaa40),0)
if(xyzzyaafi1(xyzzyaaaa40)>0)call save_enbasis_pbuffer(xyzzyaaat1(xyzz&
&yaaaa40),0)
xyzzyaaae40=xyzzyaaag1(xyzzyaaaa40)
xyzzyaaah40=xyzzyaaaj1(xyzzyaaae40+1)
xyzzyaaai40=xyzzyaaah40+1
xyzzyaaaf40=xyzzyaafs1(xyzzyaaaa40)
xyzzyaaag40=xyzzyaaaf40+1
call dcopy(xyzzyaaax1(xyzzyaaaa40),param(xyzzyaaai40),1,xyzzyaaft1(xyz&
&zyaaag40),1)
enddo
if(xyzzyaacv1>0)then
call dcopy(xyzzyaacw1,xyzzyaadr1(1),1,xyzzyaafu1(1),1)
call dcopy(xyzzyaacv1,eqn_rhs(1),1,xyzzyaafv1(1),1)
endif
else
xyzzyaaaa40=xyzzyaafm1(iparam)
xyzzyaaab40=xyzzyaafp1(iparam)
xyzzyaaac40=xyzzyaafo1(iparam)
if(xyzzyaafn1(iparam)/=5)then
select case(xyzzyaafn1(iparam))
case(1)
call save_eebasis_pbuffer(xyzzyaaaq1(xyzzyaaaa40),iparam-xyzzyaaac40+1&
&)
case(2)
call save_eebasis_pbuffer(xyzzyaaar1(xyzzyaaaa40),iparam-xyzzyaaac40+1&
&)
case(3)
call save_enbasis_pbuffer(xyzzyaaas1(xyzzyaaaa40),iparam-xyzzyaaac40+1&
&)
case(4)
call save_enbasis_pbuffer(xyzzyaaat1(xyzzyaaaa40),iparam-xyzzyaaac40+1&
&)
end select
xyzzyaaae40=xyzzyaaag1(xyzzyaaaa40)
xyzzyaaah40=xyzzyaaaj1(xyzzyaaae40+1)
xyzzyaaai40=xyzzyaaah40+1
xyzzyaaaf40=xyzzyaafs1(xyzzyaaaa40)
xyzzyaaag40=xyzzyaaaf40+xyzzyaaax1(xyzzyaaaa40)*(iparam-xyzzyaaab40+1)&
&+1
call dcopy(xyzzyaaax1(xyzzyaaaa40),param(xyzzyaaai40),1,xyzzyaaft1(xyz&
&zyaaag40),1)
else
xyzzyaaad40=xyzzyaafq1(iparam)
xyzzyaaae40=xyzzyaaag1(xyzzyaaaa40)
xyzzyaaah40=xyzzyaaaj1(xyzzyaaae40+xyzzyaaad40)
xyzzyaaai40=xyzzyaaah40+1
xyzzyaaaf40=xyzzyaafs1(xyzzyaaaa40)
xyzzyaaag40=xyzzyaaaf40+xyzzyaaax1(xyzzyaaaa40)*(iparam-xyzzyaaab40+1)&
&+xyzzyaaaw1(xyzzyaaaa40)*(xyzzyaaad40-1)+1
call dcopy(xyzzyaaaw1(xyzzyaaaa40),param(xyzzyaaai40),1,xyzzyaaft1(xyz&
&zyaaag40),1)
endif
endif
end subroutine xyzzyaahr1
subroutine xyzzyaahs1(iparam,restore_hint)
implicit none
integer,intent(in) :: iparam,restore_hint
integer xyzzyaaaa41,xyzzyaaab41,xyzzyaaac41,xyzzyaaad41,xyzzyaaae41,xy&
&zzyaaaf41,xyzzyaaag41,xyzzyaaah41,xyzzyaaai41
if(iparam==0)then
if(restore_hint==0)then
do xyzzyaaaa41=1,xyzzyaaaa1
if(xyzzyaafh1(xyzzyaaaa41)>0)call restore_eebasis_pbuffer(xyzzyaaaq1(x&
&yzzyaaaa41),0)
if(xyzzyaafg1(xyzzyaaaa41)>0)call restore_eebasis_pbuffer(xyzzyaaar1(x&
&yzzyaaaa41),0)
if(xyzzyaafj1(xyzzyaaaa41)>0)call restore_enbasis_pbuffer(xyzzyaaas1(x&
&yzzyaaaa41),0)
if(xyzzyaafi1(xyzzyaaaa41)>0)call restore_enbasis_pbuffer(xyzzyaaat1(x&
&yzzyaaaa41),0)
xyzzyaaae41=xyzzyaaag1(xyzzyaaaa41)
xyzzyaaaf41=xyzzyaafs1(xyzzyaaaa41)
xyzzyaaag41=xyzzyaaaf41+1
xyzzyaaah41=xyzzyaaaj1(xyzzyaaae41+1)
xyzzyaaai41=xyzzyaaah41+1
call dcopy(xyzzyaaax1(xyzzyaaaa41),xyzzyaaft1(xyzzyaaag41),1,param(xyz&
&zyaaai41),1)
enddo
if(xyzzyaacv1>0)then
call dcopy(xyzzyaacw1,xyzzyaafu1(1),1,xyzzyaadr1(1),1)
call dcopy(xyzzyaacv1,xyzzyaafv1(1),1,eqn_rhs(1),1)
endif
else
xyzzyaaaa41=xyzzyaafm1(restore_hint)
if(xyzzyaafn1(restore_hint)/=5)then
select case(xyzzyaafn1(restore_hint))
case(1)
call restore_eebasis_pbuffer(xyzzyaaaq1(xyzzyaaaa41),0)
case(2)
call restore_eebasis_pbuffer(xyzzyaaar1(xyzzyaaaa41),0)
case(3)
call restore_enbasis_pbuffer(xyzzyaaas1(xyzzyaaaa41),0)
case(4)
call restore_enbasis_pbuffer(xyzzyaaat1(xyzzyaaaa41),0)
end select
xyzzyaaae41=xyzzyaaag1(xyzzyaaaa41)
xyzzyaaaf41=xyzzyaafs1(xyzzyaaaa41)
xyzzyaaag41=xyzzyaaaf41+1
xyzzyaaah41=xyzzyaaaj1(xyzzyaaae41+1)
xyzzyaaai41=xyzzyaaah41+1
call dcopy(xyzzyaaax1(xyzzyaaaa41),xyzzyaaft1(xyzzyaaag41),1,param(xyz&
&zyaaai41),1)
do xyzzyaaad41=xyzzyaaae41+1,xyzzyaaae41+xyzzyaaav1(xyzzyaaaa41)
if(xyzzyaacz1(xyzzyaaad41)>0)then
call dcopy(xyzzyaacz1(xyzzyaaad41)*xyzzyaaaw1(xyzzyaaaa41),xyzzyaafu1(&
&xyzzyaadk1(xyzzyaaad41)),1,xyzzyaadr1(xyzzyaadk1(xyzzyaaad41)),1)
call dcopy(xyzzyaacz1(xyzzyaaad41),xyzzyaafv1(xyzzyaadg1(xyzzyaaad41))&
&,1,eqn_rhs(xyzzyaadg1(xyzzyaaad41)),1)
endif
enddo
else
xyzzyaaad41=xyzzyaafq1(restore_hint)
xyzzyaaae41=xyzzyaaag1(xyzzyaaaa41)
xyzzyaaah41=xyzzyaaaj1(xyzzyaaae41+xyzzyaaad41)
xyzzyaaai41=xyzzyaaah41+1
xyzzyaaaf41=xyzzyaafs1(xyzzyaaaa41)
xyzzyaaag41=xyzzyaaaf41+xyzzyaaaw1(xyzzyaaaa41)*(xyzzyaaad41-1)+1
call dcopy(xyzzyaaaw1(xyzzyaaaa41),xyzzyaaft1(xyzzyaaag41),1,param(xyz&
&zyaaai41),1)
endif
endif
else
xyzzyaaaa41=xyzzyaafm1(iparam)
xyzzyaaab41=xyzzyaafp1(iparam)
xyzzyaaac41=xyzzyaafo1(iparam)
if(xyzzyaafn1(iparam)/=5)then
select case(xyzzyaafn1(iparam))
case(1)
call restore_eebasis_pbuffer(xyzzyaaaq1(xyzzyaaaa41),iparam-xyzzyaaac4&
&1+1)
case(2)
call restore_eebasis_pbuffer(xyzzyaaar1(xyzzyaaaa41),iparam-xyzzyaaac4&
&1+1)
case(3)
call restore_enbasis_pbuffer(xyzzyaaas1(xyzzyaaaa41),iparam-xyzzyaaac4&
&1+1)
case(4)
call restore_enbasis_pbuffer(xyzzyaaat1(xyzzyaaaa41),iparam-xyzzyaaac4&
&1+1)
end select
xyzzyaaae41=xyzzyaaag1(xyzzyaaaa41)
xyzzyaaah41=xyzzyaaaj1(xyzzyaaae41+1)
xyzzyaaai41=xyzzyaaah41+1
xyzzyaaaf41=xyzzyaafs1(xyzzyaaaa41)
xyzzyaaag41=xyzzyaaaf41+xyzzyaaax1(xyzzyaaaa41)*(iparam-xyzzyaaab41+1)&
&+1
call dcopy(xyzzyaaax1(xyzzyaaaa41),xyzzyaaft1(xyzzyaaag41),1,param(xyz&
&zyaaai41),1)
else
xyzzyaaad41=xyzzyaafq1(iparam)
xyzzyaaae41=xyzzyaaag1(xyzzyaaaa41)
xyzzyaaah41=xyzzyaaaj1(xyzzyaaae41+xyzzyaaad41)
xyzzyaaai41=xyzzyaaah41+1
xyzzyaaaf41=xyzzyaafs1(xyzzyaaaa41)
xyzzyaaag41=xyzzyaaaf41+xyzzyaaax1(xyzzyaaaa41)*(iparam-xyzzyaaab41+1)&
&+xyzzyaaaw1(xyzzyaaaa41)*(xyzzyaaad41-1)+1
call dcopy(xyzzyaaaw1(xyzzyaaaa41),xyzzyaaft1(xyzzyaaag41),1,param(xyz&
&zyaaai41),1)
endif
endif
end subroutine xyzzyaahs1
subroutine invalidate_param1_gjastrow(is,iparam)
implicit none
integer,intent(in) :: is,iparam
integer xyzzyaaaa42,xyzzyaaab42
xyzzyaaek1(:,is)=.false.
xyzzyaael1(:,is)=.false.
xyzzyaaem1(:,is)=.false.
xyzzyaaen1(:,is)=.false.
xyzzyaaeo1(:,is)=.false.
xyzzyaaep1(is)=.false.
xyzzyaaeq1(:,is)=.false.
xyzzyaaer1(:,is)=.false.
xyzzyaaes1(:,is)=.false.
xyzzyaaab42=xyzzyaafm1(iparam)
xyzzyaaaa42=xyzzyaafo1(iparam)
select case(xyzzyaafn1(iparam))
case(1)
call invalidate_param1_eebasis(is,xyzzyaaaq1(xyzzyaaab42),iparam-xyzzy&
&aaaa42+1)
xyzzyaaet1(:,xyzzyaaab42,is)=.false.
xyzzyaaeu1(:,xyzzyaaab42,is)=.false.
case(2)
call invalidate_param1_eebasis(is,xyzzyaaar1(xyzzyaaab42),iparam-xyzzy&
&aaaa42+1)
xyzzyaaet1(:,xyzzyaaab42,is)=.false.
xyzzyaaeu1(:,xyzzyaaab42,is)=.false.
case(3)
call invalidate_param1_enbasis(is,xyzzyaaas1(xyzzyaaab42),iparam-xyzzy&
&aaaa42+1)
xyzzyaaet1(:,xyzzyaaab42,is)=.false.
xyzzyaaeu1(:,xyzzyaaab42,is)=.false.
case(4)
call invalidate_param1_enbasis(is,xyzzyaaat1(xyzzyaaab42),iparam-xyzzy&
&aaaa42+1)
xyzzyaaet1(:,xyzzyaaab42,is)=.false.
xyzzyaaeu1(:,xyzzyaaab42,is)=.false.
case(5)
xyzzyaaet1(:,xyzzyaaab42,is)=.false.
xyzzyaaeu1(:,xyzzyaaab42,is)=.false.
end select
end subroutine invalidate_param1_gjastrow
subroutine invalidate_params_gjastrow(ignore)
implicit none
logical,intent(in) :: ignore(:)
if(any(.not.ignore))then
xyzzyaafc1(:)=.false.
xyzzyaafd1(:)=.false.
xyzzyaafe1(:)=.false.
xyzzyaaff1(:)=.false.
endif
end subroutine invalidate_params_gjastrow
subroutine get_linear_basis_gjastrow(is,nparam,nparam_all,ignore,lbasi&
&s_grad_f,lbasis_lap_f,grad_j0,lap_j0)
implicit none
integer,intent(in) :: is,nparam,nparam_all
real(dp),intent(out) :: lbasis_grad_f(3,netot,nparam),lbasis_lap_f(npa&
&ram),grad_j0(3,netot),lap_j0
logical,intent(in) :: ignore(nparam_all)
integer xyzzyaaaa44,xyzzyaaab44,xyzzyaaac44,xyzzyaaad44,xyzzyaaae44,xy&
&zzyaaaf44
real(dp) xyzzyaaag44,xyzzyaaah44,xyzzyaaai44(3),xyzzyaaaj44
complex(dp) xyzzyaaak44(3,netot),xyzzyaaal44(netot)
call timer('JASTROW',.true.)
lbasis_grad_f=0.d0
lbasis_lap_f=0.d0
grad_j0=0.d0
lap_j0=0.d0
xyzzyaaaa44=0
do xyzzyaaac44=1,nspin
xyzzyaaag44=inv_pmass(xyzzyaaac44)
xyzzyaaah44=sqrt(xyzzyaaag44)
do xyzzyaaab44=1,nele(xyzzyaaac44)
xyzzyaaaa44=xyzzyaaaa44+1
call xyzzyaahn1(xyzzyaaaa44,is,.false.,.true.,.true.,.true.)
lap_j0=lap_j0+sum(xyzzyaaej1(xyzzyaaaa44,1:xyzzyaaaa1,is))*xyzzyaaag44
grad_j0(1:3,xyzzyaaaa44)=xyzzyaadz1(1:3,xyzzyaaaa44,is)*xyzzyaaah44
enddo
enddo
xyzzyaaad44=0
do xyzzyaaae44=1,nparam_all
if(ignore(xyzzyaaae44))cycle
xyzzyaaad44=xyzzyaaad44+1
xyzzyaaaf44=xyzzyaafr1(xyzzyaaae44)+xyzzyaaaj1(xyzzyaaag1(xyzzyaafm1(x&
&yzzyaaae44))+xyzzyaafq1(xyzzyaaae44))
call wfn_aderiv_gjastrow(is,xyzzyaaae44,xyzzyaaak44,xyzzyaaal44)
xyzzyaaaa44=0
do xyzzyaaac44=1,nspin
xyzzyaaag44=inv_pmass(xyzzyaaac44)
xyzzyaaah44=sqrt(xyzzyaaag44)
do xyzzyaaab44=1,nele(xyzzyaaac44)
xyzzyaaaa44=xyzzyaaaa44+1
xyzzyaaai44(1:3)=dble(xyzzyaaak44(1:3,xyzzyaaaa44))*xyzzyaaah44
xyzzyaaaj44=dble(xyzzyaaal44(xyzzyaaaa44))*xyzzyaaag44
grad_j0(1:3,xyzzyaaaa44)=grad_j0(:,xyzzyaaaa44)-param(xyzzyaaaf44)*xyz&
&zyaaai44(1:3)
lap_j0=lap_j0-param(xyzzyaaaf44)*xyzzyaaaj44
lbasis_grad_f(1:3,xyzzyaaaa44,xyzzyaaad44)=xyzzyaaai44(1:3)
lbasis_lap_f(xyzzyaaad44)=lbasis_lap_f(xyzzyaaad44)+xyzzyaaaj44
enddo
enddo
enddo
call timer('JASTROW',.false.)
end subroutine get_linear_basis_gjastrow
subroutine wfn_aderiv_gjastrow(is,iparam_in,dloggrad,dloglap,dlogval)
implicit none
integer,intent(in) :: is,iparam_in
complex(dp),intent(inout) :: dloggrad(3,netot),dloglap(netot)
complex(dp),intent(inout),optional :: dlogval
integer xyzzyaaaa45,xyzzyaaab45,xyzzyaaac45,xyzzyaaad45,xyzzyaaae45,xy&
&zzyaaaf45,xyzzyaaag45,xyzzyaaah45,xyzzyaaai45,xyzzyaaaj45
integer,allocatable,save :: xyzzyaaak45(:)
real(dp) xyzzyaaal45,xyzzyaaam45(3,netot),xyzzyaaan45(netot)
real(dp),allocatable,save :: xyzzyaaao45(:)
logical xyzzyaaap45
xyzzyaaap45=present(dlogval)
xyzzyaaai45=maxval(xyzzyaaaw1)
if(allocated(xyzzyaaak45))then
if(size(xyzzyaaak45,1)<xyzzyaaai45)then
deallocate(xyzzyaaak45,xyzzyaaao45)
else
xyzzyaaai45=-1
endif
endif
if(xyzzyaaai45>0)then
allocate(xyzzyaaak45(xyzzyaaai45),xyzzyaaao45(xyzzyaaai45),stat=xyzzya&
&aaj45)
call check_alloc(xyzzyaaaj45,'WFN_ADERIV_GJASTROW','plist, pcoeff')
endif
xyzzyaaaa45=xyzzyaafm1(iparam_in)
xyzzyaaab45=xyzzyaaag1(xyzzyaaaa45)+xyzzyaafq1(iparam_in)
xyzzyaaac45=xyzzyaafr1(iparam_in)
xyzzyaaad45=xyzzyaaaj1(xyzzyaaab45)
xyzzyaaae45=xyzzyaadg1(xyzzyaaab45)
xyzzyaaaf45=xyzzyaadk1(xyzzyaaab45)
call timer('GET_JASTROW_ADERIV_'//trim(i2s(xyzzyaaaa45)),.true.)
if(xyzzyaacz1(xyzzyaaab45)>0)then
call param_dep_list(xyzzyaaac45,xyzzyaaaw1(xyzzyaaaa45),xyzzyaacz1(xyz&
&zyaaab45),xyzzyaada1(xyzzyaaab45),xyzzyaadr1(xyzzyaaaf45),xyzzyaadp1(&
&xyzzyaaae45),xyzzyaadq1(xyzzyaaae45),xyzzyaacu1(xyzzyaaad45+1),xyzzya&
&aai45,xyzzyaaak45,xyzzyaaao45)
else
xyzzyaaai45=1
xyzzyaaak45(1)=xyzzyaaac45
xyzzyaaao45(1)=1.d0
endif
do xyzzyaaag45=1,netot
call xyzzyaahn1(xyzzyaaag45,is,xyzzyaaap45,.true.,.true.,.true.)
enddo
if(xyzzyaaap45)xyzzyaaal45=0.d0
xyzzyaaam45=0.d0
xyzzyaaan45=0.d0
do xyzzyaaah45=1,xyzzyaaai45
if(xyzzyaaap45)xyzzyaaal45=xyzzyaaal45+xyzzyaaao45(xyzzyaaah45)*dsum(n&
&etot,xyzzyaaee1(xyzzyaaad45+xyzzyaaak45(xyzzyaaah45),1,is),xyzzyaaad1&
&)
call daxpy(netot,xyzzyaaao45(xyzzyaaah45),xyzzyaaef1(1,xyzzyaaad45+xyz&
&zyaaak45(xyzzyaaah45),1,is),3*xyzzyaaad1,xyzzyaaam45(1,1),3)
call daxpy(netot,xyzzyaaao45(xyzzyaaah45),xyzzyaaef1(2,xyzzyaaad45+xyz&
&zyaaak45(xyzzyaaah45),1,is),3*xyzzyaaad1,xyzzyaaam45(2,1),3)
call daxpy(netot,xyzzyaaao45(xyzzyaaah45),xyzzyaaef1(3,xyzzyaaad45+xyz&
&zyaaak45(xyzzyaaah45),1,is),3*xyzzyaaad1,xyzzyaaam45(3,1),3)
call daxpy(netot,xyzzyaaao45(xyzzyaaah45),xyzzyaaeg1(xyzzyaaad45+xyzzy&
&aaak45(xyzzyaaah45),1,is),xyzzyaaad1,xyzzyaaan45,1)
enddo
if(xyzzyaaap45)dlogval=cmplx(xyzzyaaal45,0.d0,dp)
dloggrad=cmplx(xyzzyaaam45,0.d0,dp)
dloglap=cmplx(xyzzyaaan45,0.d0,dp)
call timer('GET_JASTROW_ADERIV_'//trim(i2s(xyzzyaaaa45)),.false.)
end subroutine wfn_aderiv_gjastrow
subroutine clear_scratch_gjastrow(is)
implicit none
integer,intent(in) :: is
xyzzyaaek1(:,is)=.false.
xyzzyaael1(:,is)=.false.
xyzzyaaem1(:,is)=.false.
xyzzyaaen1(:,is)=.false.
xyzzyaaeo1(:,is)=.false.
xyzzyaaep1(is)=.false.
xyzzyaaet1(:,:,is)=.false.
xyzzyaaeu1(:,:,is)=.false.
xyzzyaaeq1(:,is)=.false.
xyzzyaaer1(:,is)=.false.
xyzzyaaes1(:,is)=.false.
end subroutine clear_scratch_gjastrow
subroutine read_gjastrow(empty_jastrow)
implicit none
logical,intent(inout) :: empty_jastrow
integer iterm,xyzzyaaaa47,xyzzyaaab47,xyzzyaaac47,xyzzyaaad47,xyzzyaaa&
&e47,xyzzyaaaf47,xyzzyaaag47,xyzzyaaah47,xyzzyaaai47,xyzzyaaaj47,xyzzy&
&aaak47,ispin,jspin,xyzzyaaal47,xyzzyaaam47,xyzzyaaan47,xyzzyaaao47,xy&
&zzyaaap47,xyzzyaaaq47,xyzzyaaar47,xyzzyaaas47,xyzzyaaat47,xyzzyaaau47&
&,xyzzyaaav47,xyzzyaaaw47,xyzzyaaax47,xyzzyaaay47,xyzzyaaaz47,xyzzyaab&
&a47,xyzzyaabb47,xyzzyaabc47,xyzzyaabd47,xyzzyaabe47,xyzzyaabf47,xyzzy&
&aabg47,xyzzyaabh47,xyzzyaabi47,xyzzyaabj47,xyzzyaabk47,xyzzyaabl47,xy&
&zzyaabm47,xyzzyaabn47,xyzzyaabo47,xyzzyaabp47,xyzzyaabq47,xyzzyaabr47&
&,xyzzyaabs47,xyzzyaabt47,xyzzyaabu47,xyzzyaabv47,xyzzyaabw47,xyzzyaab&
&x47,xyzzyaaby47,xyzzyaabz47,xyzzyaaca47,xyzzyaacb47,xyzzyaacc47,xyzzy&
&aacd47,xyzzyaace47(0),xyzzyaacf47(1),xyzzyaacg47(1),xyzzyaach47(2)
character(8),allocatable :: xyzzyaaci47(:),xyzzyaacj47(:),xyzzyaack47(&
&:),xyzzyaacl47(:)
integer,allocatable :: xyzzyaacm47(:),xyzzyaacn47(:),xyzzyaaco47(:),xy&
&zzyaacp47(:),xyzzyaacq47(:),xyzzyaacr47(:),xyzzyaacs47(:),xyzzyaact47&
&(:),xyzzyaacu47(:),xyzzyaacv47(:),xyzzyaacw47(:),xyzzyaacx47(:),xyzzy&
&aacy47(:),xyzzyaacz47(:),xyzzyaada47(:),xyzzyaadb47(:),xyzzyaadc47(:)
integer,pointer :: xyzzyaadd47(:,:)=>null(),xyzzyaade47(:,:)=>null(),x&
&yzzyaadf47(:,:)=>null(),xyzzyaadg47(:,:)=>null(),integer_model(:)=>nu&
&ll(),xyzzyaadh47(:)=>null()
real(dp) xyzzyaadi47
real(dp),allocatable :: xyzzyaadj47(:)
integer,allocatable :: xyzzyaadk47(:)
real(dp),parameter :: xyzzyaadl47=1.d-5
real(dp),parameter :: xyzzyaadm47=1.d-20
logical exists,is_block,xyzzyaadn47,xyzzyaado47,xyzzyaadp47,xyzzyaadq4&
&7,xyzzyaadr47,xyzzyaads47,xyzzyaadt47,xyzzyaadu47,xyzzyaadv47,xyzzyaa&
&dw47,xyzzyaadx47,xyzzyaady47,xyzzyaadz47,cut_ee_is_none,cut_en_is_non&
&e,ee_needs_cutoff,en_needs_cutoff,xyzzyaaea47
logical,allocatable :: xyzzyaaeb47(:,:),xyzzyaaec47(:,:),xyzzyaaed47(:&
&,:),xyzzyaaee47(:,:),xyzzyaaef47(:),xyzzyaaeg47(:),xyzzyaaeh47(:,:),x&
&yzzyaaei47(:),xyzzyaaej47(:),xyzzyaaek47(:,:),xyzzyaael47(:,:)
character(8) f_zero
character(20) word1,word2,explain
character(80) tmpr
character(512) errmsg
character(casl_keysize) chname
character(casl_valsize) title
character(casl_valsize),pointer :: xyzzyaaem47(:)=>null(),xyzzyaaen47(&
&:)=>null()
call timer('READ_JASTROW',.true.)
if(am_master)then
call wout('General Jastrow setup')
call wout('=====================')
call wout()
endif
call query_casl_item(':parameters.casl:JASTROW',exists=exists)
if(.not.exists)then
if(am_master)then
call wout('Jastrow factor not found.')
call wout('Generating a default one.')
call wout()
endif
call xyzzyaaic1
endif
call push_casl_context(':parameters.casl:JASTROW')
call get_casl_item('Title',title,xyzzyaaab47)
if(xyzzyaaab47/=0)title='No title'
if(am_master)call wout('Title: '//trim(title))
call get_casl_item('Print determined',xyzzyaahd1,xyzzyaaab47)
if(xyzzyaaab47/=0)xyzzyaahd1=.false.
xyzzyaaaa1=0
do
call query_casl_item('TERM '//trim(i2s(xyzzyaaaa1+1)),exists=exists,is&
&_block=is_block)
if(.not.exists)exit
if(.not.is_block)call errstop_master('READ_GJASTROW','Non-block "term &
&'//trim(i2s(xyzzyaaaa1+1))//'" not implemented.')
xyzzyaaaa1=xyzzyaaaa1+1
enddo
if(xyzzyaaaa1==0)call errstop_master('READ_GJASTROW','No terms found i&
&n the JASTROW block - perhaps they are numbered wrongly?')
if(am_master)then
if(xyzzyaaaa1==1)then
call wout('Reading '//trim(i2s(xyzzyaaaa1))//' Jastrow factor term.')
else
call wout('Reading '//trim(i2s(xyzzyaaaa1))//' Jastrow factor terms.')
endif
call wout()
endif
allocate(xyzzyaaak1(xyzzyaaaa1),xyzzyaaal1(xyzzyaaaa1),xyzzyaaao1(xyzz&
&yaaaa1),xyzzyaaap1(xyzzyaaaa1),xyzzyaaau1(xyzzyaaaa1),xyzzyaaam1(xyzz&
&yaaaa1),xyzzyaaan1(xyzzyaaaa1),xyzzyaabe1(xyzzyaaaa1),xyzzyaabo1(xyzz&
&yaaaa1),xyzzyaabk1(xyzzyaaaa1),xyzzyaabl1(xyzzyaaaa1),xyzzyaabm1(xyzz&
&yaaaa1),xyzzyaabn1(xyzzyaaaa1),xyzzyaaeh47(xyzzyaadu1,xyzzyaaaa1),xyz&
&zyaabp1(xyzzyaaaa1),xyzzyaabq1(xyzzyaaaa1),xyzzyaabr1(xyzzyaaaa1),xyz&
&zyaabs1(xyzzyaaaa1),xyzzyaaaw1(xyzzyaaaa1),xyzzyaaax1(xyzzyaaaa1),xyz&
&zyaabj1(xyzzyaaaa1),xyzzyaaav1(xyzzyaaaa1),xyzzyaaag1(xyzzyaaaa1),xyz&
&zyaaay1(nspin,nspin,xyzzyaaaa1),xyzzyaaaz1(max(nitot,1),nspin,xyzzyaa&
&aa1),xyzzyaaba1(xyzzyaaaa1),xyzzyaabb1(xyzzyaaaa1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','ranks,...')
xyzzyaaak1=0
xyzzyaaal1=0
xyzzyaaao1=0
xyzzyaaap1=0
xyzzyaaau1=0
xyzzyaaam1=0
xyzzyaaan1=0
xyzzyaabe1=0
xyzzyaabo1=.false.
xyzzyaabk1=.false.
xyzzyaabl1=.false.
xyzzyaabm1=.false.
xyzzyaabn1=.false.
xyzzyaaeh47=.false.
xyzzyaabp1=.false.
xyzzyaabq1=.false.
xyzzyaabr1=.false.
xyzzyaabs1=.false.
xyzzyaaaw1=0
xyzzyaaax1=0
xyzzyaabj1=0.d0
xyzzyaaav1=0
xyzzyaaag1=0
xyzzyaaba1=0
xyzzyaabb1=0
call init_groups_ee(xyzzyaaay1(1,1,1))
call init_groups_en(xyzzyaaaz1(1,1,1))
do iterm=1,xyzzyaaaa1
xyzzyaaay1(:,:,iterm)=xyzzyaaay1(:,:,1)
xyzzyaaaz1(:,:,iterm)=xyzzyaaaz1(:,:,1)
enddo
do iterm=1,xyzzyaaaa1
call push_casl_context('TERM '//trim(i2s(iterm)))
call query_casl_item('Rank',exists=exists,is_block=is_block)
if(.not.exists)call errstop_master('READ_GJASTROW','No declared rank f&
&or term #'//trim(i2s(iterm))//'.')
if(.not.is_block)call errstop_master('READ_GJASTROW','Non-block "rank"&
& not implemented.')
call get_casl_item('Rank:%u1',xyzzyaaak1(iterm),xyzzyaaab47)
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Problem getting&
& Rank_e for term #'//trim(i2s(iterm))//'.')
call get_casl_item('Rank:%u2',xyzzyaaal1(iterm),xyzzyaaab47)
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Problem getting&
& Rank_n for term #'//trim(i2s(iterm))//'.')
if(xyzzyaaak1(iterm)<0.or.xyzzyaaal1(iterm)<0)call errstop_master('REA&
&D_GJASTROW','Negative ranks given for term #'//trim(i2s(iterm))//'.')
if(xyzzyaaak1(iterm)<1)call errstop_master('READ_GJASTROW','Zero Rank_&
&e given for term #'//trim(i2s(iterm))//'.')
if(xyzzyaaak1(iterm)+xyzzyaaal1(iterm)<2)call errstop_master('READ_GJA&
&STROW','Sum of ranks of Jastrow term #'//trim(i2s(iterm))//' is less &
&than two. Such term cannot be generated.')
xyzzyaabj1(iterm)=1.d0/dble(xyzzyaaak1(iterm))
xyzzyaaao1(iterm)=(xyzzyaaak1(iterm)*(xyzzyaaak1(iterm)-1))/2
xyzzyaaap1(iterm)=xyzzyaaak1(iterm)*xyzzyaaal1(iterm)
xyzzyaaau1(iterm)=xyzzyaaao1(iterm)+xyzzyaaap1(iterm)
call get_casl_item('e-e basis:Order',xyzzyaaac47,xyzzyaaab47)
if(xyzzyaaab47/=0)then
if(xyzzyaaao1(iterm)>0)xyzzyaaam1(iterm)=1
else
if(xyzzyaaao1(iterm)==0)call errstop_master('READ_GJASTROW','Rank_e = &
&1 for term #'//trim(i2s(iterm))//', but an e-e basis expansion order &
&was found.')
xyzzyaaam1(iterm)=xyzzyaaac47
endif
call get_casl_item('e-n basis:Order',xyzzyaaac47,xyzzyaaab47)
if(xyzzyaaab47/=0)then
if(xyzzyaaap1(iterm)>0)xyzzyaaan1(iterm)=1
else
if(xyzzyaaap1(iterm)==0)call errstop_master('READ_GJASTROW','Rank_n = &
&0 for term #'//trim(i2s(iterm))//', but an e-n basis expansion order &
&was found.')
xyzzyaaan1(iterm)=xyzzyaaac47
endif
xyzzyaaea47=.true.
call query_casl_item('Indexing',exists=exists,is_block=is_block)
if(exists.and.is_block)then
call get_casl_item('Indexing:Maximum sum',xyzzyaaac47,xyzzyaaab47)
if(xyzzyaaab47==0)then
if(xyzzyaaao1(iterm)+xyzzyaaap1(iterm)<2)call errstop_master('READ_GJA&
&STROW','Cannot use "Maximum sum" on a term whose linear parameters on&
&ly have one index.')
if(xyzzyaaac47<1)call errstop_master('READ_GJASTROW','Maximum sum cann&
&ot be less than 1.')
if(xyzzyaaac47<xyzzyaaao1(iterm)*xyzzyaaam1(iterm)+xyzzyaaap1(iterm)*x&
&yzzyaaan1(iterm))then
xyzzyaabe1(iterm)=xyzzyaaac47
xyzzyaaea47=.false.
endif
endif
call get_casl_item('Indexing:All masks',xyzzyaadz47,xyzzyaaab47)
if(xyzzyaaab47==0)then
if(xyzzyaaao1(iterm)+xyzzyaaap1(iterm)<2)call errstop_master('READ_GJA&
&STROW','Cannot use "All masks" on a term whose linear parameters only&
& have one index.')
xyzzyaabo1(iterm)=xyzzyaadz47
if(xyzzyaadz47)xyzzyaaea47=.false.
endif
call get_casl_item('e-e dot product',xyzzyaadz47,xyzzyaaab47)
if(xyzzyaaab47==0)then
if(xyzzyaadz47)then
word1='T'
else
word1='F'
endif
else
call get_casl_item('e-e dot product',word1,xyzzyaaab47)
endif
if(xyzzyaaab47==0)then
select case(trim(word1))
case('T','no repeat')
if(xyzzyaaao1(iterm)<2)call errstop_master('READ_GJASTROW','e-e dot pr&
&oduct indexing requested for term #'//trim(i2s(iterm))//', but this t&
&erm does not have two or more e-e indices as this constraint requires&
&.')
xyzzyaabk1(iterm)=.true.
xyzzyaaea47=.false.
xyzzyaabp1(iterm)=.true.
xyzzyaaeh47(xyzzyaadv1,iterm)=.true.
xyzzyaabr1(iterm)=.true.
xyzzyaabm1(iterm)=trim(word1)=='no repeat'
case('F')
continue
case default
call errstop_master('READ_GJASTROW','Invalid value given for e-e dot p&
&roduct in term #'//trim(i2s(iterm))//'.  Valid values are "normal" an&
&d "no repeat".')
end select
endif
call get_casl_item('e-n dot product',word1,xyzzyaaab47)
if(xyzzyaaab47==0)then
if(xyzzyaadz47)then
word1='T'
else
word1='F'
endif
else
call get_casl_item('e-n dot product',word1,xyzzyaaab47)
endif
if(xyzzyaaab47==0)then
select case(trim(word1))
case('T','no repeat')
if(xyzzyaaap1(iterm)<2)call errstop_master('READ_GJASTROW','e-n dot pr&
&oduct indexing requested for term #'//trim(i2s(iterm))//', but this t&
&erm does not have two or more e-n indices as this constraint requires&
&.')
xyzzyaabl1(iterm)=.true.
xyzzyaaea47=.false.
xyzzyaabq1(iterm)=.true.
xyzzyaaeh47(xyzzyaadw1,iterm)=.true.
xyzzyaabs1(iterm)=.true.
xyzzyaabn1(iterm)=trim(word1)=='no repeat'
case('F')
continue
case default
call errstop_master('READ_GJASTROW','Invalid value given for e-n dot p&
&roduct in term #'//trim(i2s(iterm))//'.  Valid values are "normal" an&
&d "no repeat".')
end select
endif
endif
call query_casl_item('e-e indexing',exists=exists,is_block=is_block)
if(exists.and.is_block)then
call get_casl_item('e-e indexing:Dot product',xyzzyaadz47,xyzzyaaab47)
if(xyzzyaaab47==0)then
if(xyzzyaadz47.and.xyzzyaaao1(iterm)<2)call errstop_master('READ_GJAST&
&ROW','e-e dot-product indexing requested for term #'//trim(i2s(iterm)&
&)//', but this term does not have two or more e-e indices as this con&
&straint requires.')
xyzzyaabk1(iterm)=xyzzyaadz47
if(xyzzyaadz47)then
xyzzyaaea47=.false.
xyzzyaabp1(iterm)=.true.
xyzzyaaeh47(xyzzyaadv1,iterm)=.true.
xyzzyaabr1(iterm)=.true.
endif
endif
call get_casl_item('e-e indexing:No repeat',xyzzyaadz47,xyzzyaaab47)
if(xyzzyaaab47==0)then
if(xyzzyaadz47.and..not.xyzzyaabk1(iterm))call errstop_master('READ_GJ&
&ASTROW','e-e no-repeat indexing constraint requested for term #'//tri&
&m(i2s(iterm))//', but this depends on the e-e dot-product constraint,&
& which has not been requested for this term.')
xyzzyaabm1(iterm)=xyzzyaadz47
endif
endif
call query_casl_item('e-n indexing',exists=exists,is_block=is_block)
if(exists.and.is_block)then
call get_casl_item('e-n indexing:Dot product',xyzzyaadz47,xyzzyaaab47)
if(xyzzyaaab47==0)then
if(xyzzyaadz47.and.xyzzyaaap1(iterm)<2)call errstop_master('READ_GJAST&
&ROW','e-n dot-product indexing requested for term #'//trim(i2s(iterm)&
&)//', but this term does not have two or more e-n indices as this con&
&straint requires.')
xyzzyaabl1(iterm)=xyzzyaadz47
if(xyzzyaadz47)then
xyzzyaaea47=.false.
xyzzyaabq1(iterm)=.true.
xyzzyaaeh47(xyzzyaadw1,iterm)=.true.
xyzzyaabs1(iterm)=.true.
endif
endif
call get_casl_item('e-n indexing:No repeat',xyzzyaadz47,xyzzyaaab47)
if(xyzzyaaab47==0)then
if(xyzzyaadz47.and..not.xyzzyaabl1(iterm))call errstop_master('READ_GJ&
&ASTROW','e-n no-repeat indexing constraint requested for term #'//tri&
&m(i2s(iterm))//', but this depends on the e-n dot-product constraint,&
& which has not been requested for this term.')
xyzzyaabn1(iterm)=xyzzyaadz47
endif
endif
tmpr='unidentified'
if(xyzzyaaap1(iterm)==0)then
if(xyzzyaaea47.and.xyzzyaaao1(iterm)==1.and..not.xyzzyaahi1)then
xyzzyaaba1(iterm)=xyzzyaafw1
tmpr='1_0'
else
xyzzyaaba1(iterm)=xyzzyaafz1
tmpr='x_0'
endif
elseif(xyzzyaaao1(iterm)==0)then
if(xyzzyaaea47.and.xyzzyaaap1(iterm)==1.and..not.xyzzyaahj1)then
xyzzyaaba1(iterm)=xyzzyaafx1
tmpr='0_1'
else
xyzzyaaba1(iterm)=xyzzyaaga1
tmpr='0_x'
endif
else
if(xyzzyaaea47.and.xyzzyaaao1(iterm)==1.and.xyzzyaaap1(iterm)==2.and..&
&not.xyzzyaahk1)then
xyzzyaaba1(iterm)=xyzzyaafy1
tmpr='1_2'
else
xyzzyaaba1(iterm)=xyzzyaagb1
tmpr='x_x'
endif
endif
if(xyzzyaaba1(iterm)==0)call errstop_master('READ_GJASTROW','One-elect&
&ron evaluation routine for term #'//trim(i2s(iterm))//' is "'//trim(t&
&mpr)//'", but it is not coded. Please report this bug.')
tmpr='unidentified'
if(xyzzyaaap1(iterm)==0)then
xyzzyaabb1(iterm)=xyzzyaagf1
tmpr='x_0'
elseif(xyzzyaaao1(iterm)==0)then
xyzzyaabb1(iterm)=xyzzyaagg1
tmpr='0_x'
else
xyzzyaabb1(iterm)=xyzzyaagh1
tmpr='x_x'
endif
if(xyzzyaabb1(iterm)==0)call errstop_master('READ_GJASTROW','All-elect&
&ron evaluation routine for term #'//trim(i2s(iterm))//' is "'//trim(t&
&mpr)//'", but it is not coded. Please report this bug.')
call pop_casl_context()
enddo
allocate(xyzzyaaah1(xyzzyaaaa1),xyzzyaaai1(xyzzyaaaa1),stat=xyzzyaaaa4&
&7)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','term_iparam0,...')
xyzzyaaah1=0
xyzzyaaai1=0
xyzzyaaab1=0
xyzzyaaad1=0
xyzzyaaac1=0
xyzzyaaax47=0
xyzzyaaae1=maxval(xyzzyaaau1)
xyzzyaaaf1=maxval(xyzzyaaak1+xyzzyaaal1)
allocate(xyzzyaabx1(2,max(1,maxval(xyzzyaaao1)),xyzzyaaaa1),xyzzyaaby1&
&(2,max(1,maxval(xyzzyaaap1)),xyzzyaaaa1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','which_e*_pair')
xyzzyaabx1=0
xyzzyaaby1=0
do iterm=1,xyzzyaaaa1
xyzzyaaan47=0
do xyzzyaaao47=1,xyzzyaaak1(iterm)
do xyzzyaaap47=xyzzyaaao47+1,xyzzyaaak1(iterm)
xyzzyaaan47=xyzzyaaan47+1
xyzzyaabx1(1,xyzzyaaan47,iterm)=xyzzyaaao47
xyzzyaabx1(2,xyzzyaaan47,iterm)=xyzzyaaap47
enddo
enddo
xyzzyaaan47=0
do xyzzyaaao47=1,xyzzyaaak1(iterm)
do xyzzyaaap47=1,xyzzyaaal1(iterm)
xyzzyaaan47=xyzzyaaan47+1
xyzzyaaby1(1,xyzzyaaan47,iterm)=xyzzyaaao47
xyzzyaaby1(2,xyzzyaaan47,iterm)=xyzzyaaap47
enddo
enddo
call generate_all_indices(xyzzyaaak1(iterm),xyzzyaaal1(iterm),xyzzyaaa&
&o1(iterm),xyzzyaaap1(iterm),xyzzyaaau1(iterm),xyzzyaaam1(iterm),xyzzy&
&aaan1(iterm),xyzzyaabk1(iterm),xyzzyaabl1(iterm),xyzzyaabm1(iterm),xy&
&zzyaabn1(iterm),xyzzyaabe1(iterm),xyzzyaabo1(iterm),xyzzyaabx1(1,1,it&
&erm),xyzzyaaby1(1,1,iterm),nparam=xyzzyaaaw1(iterm))
xyzzyaaah1(iterm)=xyzzyaaac1
xyzzyaaac1=xyzzyaaac1+xyzzyaaaw1(iterm)
xyzzyaaai1(iterm)=xyzzyaaax47+1
xyzzyaaax47=xyzzyaaax47+xyzzyaaaw1(iterm)*xyzzyaaau1(iterm)
enddo
allocate(xyzzyaacg1(xyzzyaaax47),xyzzyaach1(xyzzyaaac1),stat=xyzzyaaaa&
&47)
xyzzyaacg1=0
do iterm=1,xyzzyaaaa1
xyzzyaaay47=xyzzyaaai1(iterm)
xyzzyaaag47=xyzzyaaah1(iterm)+1
xyzzyaaah47=xyzzyaaah1(iterm)+xyzzyaaaw1(iterm)
call generate_all_indices(xyzzyaaak1(iterm),xyzzyaaal1(iterm),xyzzyaaa&
&o1(iterm),xyzzyaaap1(iterm),xyzzyaaau1(iterm),xyzzyaaam1(iterm),xyzzy&
&aaan1(iterm),xyzzyaabk1(iterm),xyzzyaabl1(iterm),xyzzyaabm1(iterm),xy&
&zzyaabn1(iterm),xyzzyaabe1(iterm),xyzzyaabo1(iterm),xyzzyaabx1(1,1,it&
&erm),xyzzyaaby1(1,1,iterm),index_list=xyzzyaacg1(xyzzyaaay47),param_n&
&ame=xyzzyaach1(xyzzyaaag47:xyzzyaaah47))
enddo
call gen_sensible_ruleset(xyzzyaaen47)
do iterm=1,xyzzyaaaa1
call push_casl_context('TERM '//trim(i2s(iterm)))
call query_casl_item('Rules',exists=exists,is_block=is_block,nchildren&
&=xyzzyaaad47)
nullify(xyzzyaaem47)
xyzzyaadn47=.not.exists
xyzzyaado47=.false.
if(xyzzyaadn47)then
xyzzyaaem47=>xyzzyaaen47
xyzzyaaad47=0
if(associated(xyzzyaaem47))xyzzyaaad47=size(xyzzyaaem47)
else
if(.not.is_block)call errstop_master('READ_GJASTROW','Non-block "Rules&
&" not implemented.')
allocate(xyzzyaaem47(xyzzyaaad47),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','ruleset')
xyzzyaaat47=0
do xyzzyaaae47=1,xyzzyaaad47
xyzzyaaat47=xyzzyaaat47+1
call get_casl_item('Rules:%u'//trim(i2s(xyzzyaaae47)),xyzzyaaem47(xyzz&
&yaaat47),xyzzyaaab47)
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Problem getting&
& rule #'//trim(i2s(xyzzyaaae47))//' in term #'//trim(i2s(iterm))//'.'&
&)
if(trim(xyzzyaaem47(xyzzyaaat47))=='default')then
xyzzyaado47=.true.
xyzzyaaat47=xyzzyaaat47-1
endif
enddo
if(xyzzyaado47.and.associated(xyzzyaaen47))then
xyzzyaaad47=xyzzyaaat47+size(xyzzyaaen47)
call resize_pointer(casl_valsize,(/xyzzyaaad47/),xyzzyaaem47)
do xyzzyaaae47=xyzzyaaat47+1,xyzzyaaad47
xyzzyaaem47(xyzzyaaae47)=xyzzyaaen47(xyzzyaaae47-xyzzyaaat47)
enddo
endif
endif
do xyzzyaaae47=1,xyzzyaaad47
call digest_rule(xyzzyaaem47(xyzzyaaae47),xyzzyaaak1(iterm)>1,xyzzyaaa&
&l1(iterm)>0,xyzzyaaay1(1,1,iterm),xyzzyaaaz1(1,1,iterm))
enddo
if(.not.xyzzyaadn47)deallocate(xyzzyaaem47)
nullify(xyzzyaaem47)
if(xyzzyaadn47.or.xyzzyaado47)then
call reconstruct_ruleset(xyzzyaaak1(iterm)>1,xyzzyaaal1(iterm)>0,xyzzy&
&aaay1(1,1,iterm),xyzzyaaaz1(1,1,iterm),xyzzyaaem47)
call delete_casl_item('Rules')
call set_casl_block('Rules',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('READ_GJASTROW',trim(errmsg)&
&)
if(associated(xyzzyaaem47))then
do xyzzyaaae47=1,size(xyzzyaaem47)
call set_casl_item('Rules:%u',trim(xyzzyaaem47(xyzzyaaae47)),errmsg)
if(len_trim(errmsg)>0)call errstop_master('READ_GJASTROW',trim(errmsg)&
&)
enddo
deallocate(xyzzyaaem47)
endif
nullify(xyzzyaaem47)
endif
call build_channels(xyzzyaaay1(1,1,iterm),xyzzyaaaz1(1,1,iterm),xyzzya&
&aak1(iterm),xyzzyaaal1(iterm),xyzzyaaav1(iterm),xyzzyaadd47,xyzzyaade&
&47)
xyzzyaaag1(iterm)=xyzzyaaab1
xyzzyaaab1=xyzzyaaab1+xyzzyaaav1(iterm)
call resize_pointer((/xyzzyaaae1,xyzzyaaab1/),xyzzyaadf47)
xyzzyaadf47(1:xyzzyaaau1(iterm),xyzzyaaag1(iterm)+1:xyzzyaaag1(iterm)+&
&xyzzyaaav1(iterm))=xyzzyaadd47
call resize_pointer((/xyzzyaaaf1,xyzzyaaab1/),xyzzyaadg47)
xyzzyaadg47(1:xyzzyaaak1(iterm)+xyzzyaaal1(iterm),xyzzyaaag1(iterm)+1:&
&xyzzyaaag1(iterm)+xyzzyaaav1(iterm))=xyzzyaade47
deallocate(xyzzyaadd47,xyzzyaade47)
xyzzyaaax1(iterm)=xyzzyaaaw1(iterm)*xyzzyaaav1(iterm)
xyzzyaaad1=xyzzyaaad1+xyzzyaaax1(iterm)
call pop_casl_context()
enddo
if(associated(xyzzyaaen47))deallocate(xyzzyaaen47)
allocate(xyzzyaacf1(xyzzyaaae1,xyzzyaaab1),xyzzyaade47(xyzzyaaaf1,xyzz&
&yaaab1),xyzzyaabz1(xyzzyaaab1),xyzzyaacj1(xyzzyaaad1),xyzzyaack1(xyzz&
&yaaad1),xyzzyaacl1(xyzzyaaad1),xyzzyaacm1(xyzzyaaad1),xyzzyaacn1(xyzz&
&yaaad1),xyzzyaacu1(xyzzyaaad1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','channel_sig')
xyzzyaacf1=xyzzyaadf47
xyzzyaade47=xyzzyaadg47
deallocate(xyzzyaadf47,xyzzyaadg47)
xyzzyaabz1=0
xyzzyaacj1=0
xyzzyaack1=0
xyzzyaacl1=0
xyzzyaacm1=0
xyzzyaacn1=0
xyzzyaacu1=.false.
allocate(xyzzyaaaj1(xyzzyaaab1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','iparam_unfold0')
xyzzyaaai47=0
do iterm=1,xyzzyaaaa1
xyzzyaaam47=xyzzyaaag1(iterm)
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
xyzzyaaaj1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaaai47
xyzzyaaai47=xyzzyaaai47+xyzzyaaaw1(iterm)
enddo
enddo
allocate(xyzzyaace1(xyzzyaaab1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','channel_name')
allocate(integer_model(xyzzyaaaf1),xyzzyaadh47(xyzzyaaae1),stat=xyzzya&
&aaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','tmps')
xyzzyaace1=''
do iterm=1,xyzzyaaaa1
call push_casl_context('TERM '//trim(i2s(iterm)))
do
call first_unread_child('Linear parameters',chname,xyzzyaaab47,flag_as&
&_read=.true.)
if(xyzzyaaab47/=0)exit
if(chname(1:7)/='channel')call errstop_master('READ_GJASTROW','Misname&
&d item under "Linear parameters" for term #'//trim(i2s(iterm))//': ex&
&pected "Channel <model>".')
chname=trim(chname(8:))
call get_sig_from_model(chname,xyzzyaaay1(1,1,iterm),xyzzyaaaz1(1,1,it&
&erm),xyzzyaaak1(iterm),xyzzyaaal1(iterm),xyzzyaadh47,xyzzyaaab47,inte&
&ger_model=integer_model)
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Could not get s&
&ignature from channel name in "Linear parameters:Channel '//trim(chna&
&me)//'" for term #'//trim(i2s(iterm))//'.')
if(any(xyzzyaadh47(1:xyzzyaaau1(iterm))==0))call errstop_master('READ_&
&GJASTROW','Channel name "'//trim(chname)//'" contains particles or pa&
&rticle pairs which do not occur in this system (e.g., you have specif&
&ied "1-1" but there is only one up-spin electron).')
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
if(all(xyzzyaadh47==xyzzyaacf1(:,xyzzyaaag1(iterm)+xyzzyaaaf47)))then
if(trim(xyzzyaace1(xyzzyaaag1(iterm)+xyzzyaaaf47))/='')call errstop_ma&
&ster('READ_GJASTROW','In "Linear parameters" for term #'//trim(i2s(it&
&erm))//': channel "'//trim(xyzzyaace1(xyzzyaaag1(iterm)+xyzzyaaaf47))&
&//'" and channel "'//trim(chname)//'" are redundant, since they corre&
&spond to the same particle groups according to the term''s rules.')
xyzzyaace1(xyzzyaaag1(iterm)+xyzzyaaaf47)=trim(chname)
xyzzyaade47(:,xyzzyaaag1(iterm)+xyzzyaaaf47)=integer_model
endif
enddo
enddo
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
if(trim(xyzzyaace1(xyzzyaaag1(iterm)+xyzzyaaaf47))/='')cycle
xyzzyaace1(xyzzyaaag1(iterm)+xyzzyaaaf47)=trim(model_string(xyzzyaaak1&
&(iterm),xyzzyaaal1(iterm),xyzzyaade47(:,xyzzyaaag1(iterm)+xyzzyaaaf47&
&)))
enddo
call pop_casl_context()
enddo
deallocate(integer_model,xyzzyaadh47)
allocate(xyzzyaaaq1(xyzzyaaaa1),xyzzyaaar1(xyzzyaaaa1),xyzzyaaas1(xyzz&
&yaaaa1),xyzzyaaat1(xyzzyaaaa1),xyzzyaabf1(xyzzyaaaa1),xyzzyaabg1(xyzz&
&yaaaa1),xyzzyaabh1(xyzzyaaaa1),xyzzyaabi1(xyzzyaaaa1),xyzzyaabt1(xyzz&
&yaaaa1),xyzzyaabv1(xyzzyaaaa1),xyzzyaabu1(xyzzyaaaa1),xyzzyaabw1(xyzz&
&yaaaa1),xyzzyaabc1(xyzzyaaaa1),xyzzyaabd1(xyzzyaaaa1),xyzzyaaci47(xyz&
&zyaaaa1),xyzzyaacj47(xyzzyaaaa1),xyzzyaack47(xyzzyaaaa1),xyzzyaacl47(&
&xyzzyaaaa1),xyzzyaaef47(xyzzyaaaa1),xyzzyaaeg47(xyzzyaaaa1),param(xyz&
&zyaaad1),xyzzyaaci1(xyzzyaaad1),xyzzyaacq1(xyzzyaaad1),xyzzyaaco1(xyz&
&zyaaad1),xyzzyaacr1(xyzzyaaad1),xyzzyaacp1(xyzzyaaad1),xyzzyaacs1(xyz&
&zyaaad1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','param,...')
xyzzyaaaq1=0
xyzzyaaar1=0
xyzzyaaas1=0
xyzzyaaat1=0
xyzzyaabf1=0
xyzzyaabg1=0
xyzzyaabh1=0
xyzzyaabi1=0
xyzzyaabt1=.true.
xyzzyaabv1=.true.
xyzzyaabu1=.true.
xyzzyaabw1=.true.
param=0.d0
xyzzyaaci1=pflag_unset
xyzzyaacq1=.false.
xyzzyaaco1=0.d0
xyzzyaacr1=.false.
xyzzyaacp1=0.d0
xyzzyaacs1=.false.
empty_jastrow=.true.
do iterm=1,xyzzyaaaa1
call push_casl_context('TERM '//trim(i2s(iterm)))
ee_needs_cutoff=.false.
cut_ee_is_none=.true.
xyzzyaadt47=.true.
xyzzyaadx47=.true.
xyzzyaabc1(iterm)=1
if(xyzzyaaao1(iterm)>0)then
call read_eebasis('e-e basis','basis',xyzzyaaay1(1,1,iterm),xyzzyaaam1&
&(iterm),xyzzyaaaq1(iterm),xyzzyaabt1(iterm),xyzzyaabf1(iterm),ee_need&
&s_cutoff,xyzzyaadv47,xyzzyaaci47(iterm),xyzzyaads47,xyzzyaabf47)
call read_eebasis('e-e cutoff','cutoff',xyzzyaaay1(1,1,iterm),1,xyzzya&
&aar1(iterm),xyzzyaabu1(iterm),xyzzyaabh1(iterm),xyzzyaadz47,xyzzyaadw&
&47,xyzzyaack47(iterm),cut_ee_is_none,xyzzyaabg47)
xyzzyaadt47=xyzzyaads47.and.cut_ee_is_none
xyzzyaadx47=xyzzyaadv47.and.xyzzyaadw47
xyzzyaabc1(iterm)=xyzzyaabf47*xyzzyaabg47
endif
en_needs_cutoff=.false.
cut_en_is_none=.true.
xyzzyaadu47=.true.
xyzzyaady47=.true.
xyzzyaabd1(iterm)=1
if(xyzzyaaap1(iterm)>0)then
call read_enbasis('e-n basis','basis',xyzzyaaaz1(1,1,iterm),xyzzyaaan1&
&(iterm),xyzzyaaas1(iterm),xyzzyaabv1(iterm),xyzzyaabg1(iterm),en_need&
&s_cutoff,xyzzyaadv47,xyzzyaaci47(iterm),xyzzyaads47,xyzzyaabf47)
call read_enbasis('e-n cutoff','cutoff',xyzzyaaaz1(1,1,iterm),1,xyzzya&
&aat1(iterm),xyzzyaabw1(iterm),xyzzyaabi1(iterm),xyzzyaadz47,xyzzyaadw&
&47,xyzzyaack47(iterm),cut_en_is_none,xyzzyaabg47)
xyzzyaadu47=xyzzyaads47.and.cut_en_is_none
xyzzyaady47=xyzzyaadv47.and.xyzzyaadw47
xyzzyaabd1(iterm)=xyzzyaabf47*xyzzyaabg47
endif
if(ee_needs_cutoff.and.cut_ee_is_none.and.cut_en_is_none)call errstop_&
&master('READ_GJASTROW','In term '//trim(i2s(iterm))//' selected e-e b&
&asis requires use of an e-e (or e-n) cut-off for this system, but no &
&cut-off function was used.  Possibly periodic system with a non-perio&
&dic, uncut basis function?')
if(en_needs_cutoff.and.cut_ee_is_none.and.cut_en_is_none)call errstop_&
&master('READ_GJASTROW','In term '//trim(i2s(iterm))//' selected e-n b&
&asis requires use of an e-n (or e-e) cut-off for this system, but no &
&cut-off function was used.  Possibly periodic system with a non-perio&
&dic, uncut basis function?')
if(xyzzyaadt47.and.xyzzyaadu47)call errstop_master('READ_GJASTROW','Te&
&rm '//trim(i2s(iterm))//' is empty - the term should have at least on&
&e e-e basis or e-e cut-off or e-n basis or e-n cut-off.')
if(xyzzyaabt1(iterm).and.xyzzyaabu1(iterm).and.xyzzyaabv1(iterm).and.x&
&yzzyaabw1(iterm))then
select case(xyzzyaaba1(iterm))
case(xyzzyaafw1)
xyzzyaaba1(iterm)=xyzzyaagc1
case(xyzzyaafx1)
xyzzyaaba1(iterm)=xyzzyaagd1
case(xyzzyaafy1)
xyzzyaaba1(iterm)=xyzzyaage1
end select
endif
xyzzyaaef47(iterm)=.false.
xyzzyaaeg47(iterm)=.false.
xyzzyaaay47=xyzzyaaai1(iterm)
do xyzzyaabd47=1,xyzzyaaaw1(iterm)
xyzzyaaao47=0
if(xyzzyaaao1(iterm)>0)then
xyzzyaaaz47=xyzzyaaay47+(xyzzyaabd47-1)*xyzzyaaau1(iterm)
xyzzyaaba47=xyzzyaaay47+(xyzzyaabd47-1)*xyzzyaaau1(iterm)+xyzzyaaao1(i&
&term)-1
if(cut_ee_is_none)then
xyzzyaaao47=xyzzyaaao47+count(xyzzyaacg1(xyzzyaaaz47:xyzzyaaba47)/=0.a&
&nd.xyzzyaacg1(xyzzyaaaz47:xyzzyaaba47)/=xyzzyaabf1(iterm))
else
xyzzyaaao47=xyzzyaaao47+count(xyzzyaacg1(xyzzyaaaz47:xyzzyaaba47)/=0)
endif
endif
xyzzyaaap47=0
if(xyzzyaaap1(iterm)>0)then
xyzzyaabb47=xyzzyaaay47+(xyzzyaabd47-1)*xyzzyaaau1(iterm)+xyzzyaaao1(i&
&term)
xyzzyaabc47=xyzzyaaay47+(xyzzyaabd47-1)*xyzzyaaau1(iterm)+xyzzyaaau1(i&
&term)-1
if(cut_en_is_none)then
xyzzyaaap47=xyzzyaaap47+count(xyzzyaacg1(xyzzyaabb47:xyzzyaabc47)/=0.a&
&nd.xyzzyaacg1(xyzzyaabb47:xyzzyaabc47)/=xyzzyaabg1(iterm))
else
xyzzyaaap47=xyzzyaaap47+count(xyzzyaacg1(xyzzyaabb47:xyzzyaabc47)/=0)
endif
endif
if(xyzzyaadx47)then
if(xyzzyaaao47==1.and.xyzzyaaap47==0)then
xyzzyaaef47(iterm)=.true.
endif
endif
if(xyzzyaady47)then
if(xyzzyaaao47==0.and.xyzzyaaap47==1)then
xyzzyaaeg47(iterm)=.true.
endif
endif
if(xyzzyaaef47(iterm).and.xyzzyaaeg47(iterm))exit
enddo
xyzzyaabe47=0
xyzzyaaay47=xyzzyaaai1(iterm)
do xyzzyaabd47=1,xyzzyaaaw1(iterm)
xyzzyaabe47=xyzzyaabd47
if(xyzzyaaao1(iterm)>0)then
xyzzyaaaz47=xyzzyaaay47+(xyzzyaabd47-1)*xyzzyaaau1(iterm)
xyzzyaaba47=xyzzyaaay47+(xyzzyaabd47-1)*xyzzyaaau1(iterm)+xyzzyaaao1(i&
&term)-1
if(.not.((xyzzyaabh1(iterm)>0.and.all(xyzzyaacg1(xyzzyaaaz47:xyzzyaaba&
&47)==xyzzyaabf1(iterm))).or.all(xyzzyaacg1(xyzzyaaaz47:xyzzyaaba47)==&
&0)))xyzzyaabe47=0
endif
if(xyzzyaaap1(iterm)>0)then
xyzzyaabb47=xyzzyaaay47+(xyzzyaabd47-1)*xyzzyaaau1(iterm)+xyzzyaaao1(i&
&term)
xyzzyaabc47=xyzzyaaay47+(xyzzyaabd47-1)*xyzzyaaau1(iterm)+xyzzyaaau1(i&
&term)-1
if(.not.((xyzzyaabi1(iterm)>0.and.all(xyzzyaacg1(xyzzyaabb47:xyzzyaabc&
&47)==xyzzyaabg1(iterm))).or.all(xyzzyaacg1(xyzzyaabb47:xyzzyaabc47)==&
&0)))xyzzyaabe47=0
endif
if(xyzzyaabe47/=0)exit
enddo
xyzzyaaam47=xyzzyaaag1(iterm)
xyzzyaaag47=xyzzyaaah1(iterm)+1
xyzzyaaah47=xyzzyaaah1(iterm)+xyzzyaaaw1(iterm)
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
xyzzyaaaj47=xyzzyaaaj1(xyzzyaaam47+xyzzyaaaf47)+1
xyzzyaaak47=xyzzyaaaj1(xyzzyaaam47+xyzzyaaaf47)+xyzzyaaaw1(iterm)
call read_gparam_channel('Linear parameters:Channel '//trim(xyzzyaace1&
&(xyzzyaaam47+xyzzyaaaf47)),xyzzyaaaw1(iterm),xyzzyaach1(xyzzyaaag47:x&
&yzzyaaah47),param(xyzzyaaaj47:xyzzyaaak47),xyzzyaaci1(xyzzyaaaj47:xyz&
&zyaaak47),xyzzyaacq1(xyzzyaaaj47:xyzzyaaak47),xyzzyaaco1(xyzzyaaaj47:&
&xyzzyaaak47),xyzzyaacr1(xyzzyaaaj47:xyzzyaaak47),xyzzyaacp1(xyzzyaaaj&
&47:xyzzyaaak47),xyzzyaacs1(xyzzyaaaj47:xyzzyaaak47))
if(any(xyzzyaaci1(xyzzyaaaj47:xyzzyaaak47)/=pflag_unset))empty_jastrow&
&=.false.
if(xyzzyaabe47>0)xyzzyaaci1(xyzzyaaaj47+xyzzyaabe47-1)=pflag_fix
enddo
call pop_casl_context()
enddo
call gbasis_reset_aniso()
do iterm=1,xyzzyaaaa1
select case(xyzzyaaba1(iterm))
case(xyzzyaagc1,xyzzyaagd1,xyzzyaage1)
continue
case default
call eebasis_set_aniso(xyzzyaaaq1(iterm))
call eebasis_set_aniso(xyzzyaaar1(iterm))
call enbasis_set_aniso(xyzzyaaas1(iterm))
call enbasis_set_aniso(xyzzyaaat1(iterm))
end select
enddo
allocate(xyzzyaaei47(xyzzyaaab1),xyzzyaaej47(xyzzyaaab1),xyzzyaaca1(ma&
&x(1,maxval(xyzzyaacf1)),xyzzyaaab1),xyzzyaacb1(max(1,maxval(xyzzyaacf&
&1)),xyzzyaaab1),xyzzyaacc1(max(1,maxval(xyzzyaacf1)),xyzzyaaab1),xyzz&
&yaacd1(max(1,maxval(xyzzyaacf1)),xyzzyaaab1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','cusp')
xyzzyaaei47=.false.
xyzzyaaej47=.false.
xyzzyaaca1=0.d0
xyzzyaacb1=0.d0
xyzzyaacc1=.false.
xyzzyaacd1=.false.
allocate(xyzzyaaeb47(nspin,nspin),xyzzyaaec47(nspin,nitot),xyzzyaaed47&
&(nspin,nspin),xyzzyaaee47(nspin,nitot),xyzzyaaek47(nspin,nspin),xyzzy&
&aael47(nspin,nitot),stat=xyzzyaaaa47)
xyzzyaaeb47=.false.
xyzzyaaec47=.false.
xyzzyaaek47=.false.
xyzzyaael47=.false.
xyzzyaaed47=.false.
xyzzyaaee47=.false.
xyzzyaadp47=.false.
do iterm=1,xyzzyaaaa1
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
xyzzyaaca47=0
xyzzyaacb47=0
do xyzzyaaao47=1,xyzzyaaak1(iterm)
ispin=xyzzyaade47(xyzzyaaao47,xyzzyaaag1(iterm)+xyzzyaaaf47)
do xyzzyaaap47=xyzzyaaao47+1,xyzzyaaak1(iterm)
jspin=xyzzyaade47(xyzzyaaap47,xyzzyaaag1(iterm)+xyzzyaaaf47)
xyzzyaaca47=xyzzyaaca47+1
xyzzyaacd47=xyzzyaacf1(xyzzyaaca47,xyzzyaaag1(iterm)+xyzzyaaaf47)
if(.not.ee_cusp_in_orbital(ispin,jspin))then
xyzzyaaca1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)=ee_kato_gamma(is&
&pin,jspin,f_zero=f_zero)
if(.not.(xyzzyaaci47(iterm)==f_zero.or.xyzzyaack47(iterm)==f_zero))xyz&
&zyaaca1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)=0.d0
endif
enddo
do xyzzyaaap47=1,xyzzyaaal1(iterm)
xyzzyaaal47=xyzzyaade47(xyzzyaaak1(iterm)+xyzzyaaap47,xyzzyaaag1(iterm&
&)+xyzzyaaaf47)
xyzzyaacb47=xyzzyaacb47+1
xyzzyaacd47=xyzzyaacf1(xyzzyaaao1(iterm)+xyzzyaacb47,xyzzyaaag1(iterm)&
&+xyzzyaaaf47)
if(allocated(zion))then
xyzzyaacb1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)=en_kato_gamma(is&
&pin,zion(iontype(xyzzyaaal47)))
else
xyzzyaacb1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)=en_kato_gamma(is&
&pin,0.d0)
endif
enddo
enddo
enddo
enddo
do iterm=1,xyzzyaaaa1
call push_casl_context('TERM '//trim(i2s(iterm)))
call query_casl_item('e-e cusp',exists=exists)
if(exists)then
call get_casl_item('e-e cusp',xyzzyaadq47,xyzzyaaab47)
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Problem reading&
& "e-e cusp" item.')
if(xyzzyaadq47.and..not.xyzzyaaef47(iterm))call errstop_master('READ_G&
&JASTROW','e-e cusps requested for term #'//trim(i2s(iterm))//', but t&
&his term does not allow e-e cusps to be applied.')
else
xyzzyaadq47=xyzzyaaef47(iterm).and..not.xyzzyaadp47
endif
if(xyzzyaadq47)xyzzyaadp47=.true.
call query_casl_item('e-n cusp',exists=exists)
if(exists)then
call get_casl_item('e-n cusp',xyzzyaadr47,xyzzyaaab47)
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Problem reading&
& "e-n cusp" item.')
if(xyzzyaadr47.and..not.xyzzyaaeg47(iterm))call errstop_master('READ_G&
&JASTROW','e-n cusps requested for term #'//trim(i2s(iterm))//', but t&
&his term does not allow e-n cusps to be applied.')
else
xyzzyaadr47=.false.
endif
xyzzyaaek47=.false.
call query_casl_item('Waive e-e cusp',exists=exists,is_block=is_block)
if(exists.and.is_block)then
xyzzyaaao47=0
do
xyzzyaaao47=xyzzyaaao47+1
call get_casl_item('Waive e-e cusp:%u'//trim(i2s(xyzzyaaao47)),tmpr,xy&
&zzyaaab47)
if(xyzzyaaab47/=0)exit
call parse_model(trim(tmpr),2,0,xyzzyaach47,xyzzyaace47,xyzzyaaab47)
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Syntax error in&
& element '//trim(i2s(xyzzyaaao47))//' of "Waive e-e cusp".')
if(any(xyzzyaach47<0).or.any(xyzzyaach47>nspin))call errstop_master('R&
&EAD_GJASTROW','Particle type indices out of range in element '//trim(&
&i2s(xyzzyaaao47))//' of "Waive e-e cusp".')
where(xyzzyaaay1(:,:,iterm)==xyzzyaaay1(xyzzyaach47(1),xyzzyaach47(2),&
&iterm))xyzzyaaek47(:,:)=.true.
enddo
endif
xyzzyaael47=.false.
call query_casl_item('Waive e-n cusp',exists=exists,is_block=is_block)
if(exists.and.is_block)then
xyzzyaaao47=0
do
xyzzyaaao47=xyzzyaaao47+1
call get_casl_item('Waive e-n cusp:%u'//trim(i2s(xyzzyaaao47)),tmpr,xy&
&zzyaaab47)
if(xyzzyaaab47/=0)exit
call parse_model(trim(tmpr),1,1,xyzzyaacf47,xyzzyaacg47,xyzzyaaab47)
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Syntax error in&
& element '//trim(i2s(xyzzyaaao47))//' of "Waive e-n cusp".')
if(xyzzyaacf47(1)<0.or.xyzzyaacf47(1)>nspin)call errstop_master('READ_&
&GJASTROW','Particle type indices out of range in element '//trim(i2s(&
&xyzzyaaao47))//' of "Waive e-n cusp".')
if(xyzzyaacg47(1)<0.or.xyzzyaacg47(1)>nitot)call errstop_master('READ_&
&GJASTROW','Nucleus indices out of range in element '//trim(i2s(xyzzya&
&aao47))//' of "Waive e-n cusp".')
where(xyzzyaaaz1(:,:,iterm)==xyzzyaaaz1(xyzzyaacf47(1),xyzzyaacg47(1),&
&iterm))xyzzyaael47(:,:)=.true.
enddo
endif
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
call get_casl_item('Linear parameters:Channel '//trim(xyzzyaace1(xyzzy&
&aaag1(iterm)+xyzzyaaaf47))//':e-e cusp',xyzzyaaei47(xyzzyaaag1(iterm)&
&+xyzzyaaaf47),xyzzyaaab47)
if(xyzzyaaab47/=0)xyzzyaaei47(xyzzyaaag1(iterm)+xyzzyaaaf47)=xyzzyaadq&
&47
call get_casl_item('Linear parameters:Channel '//trim(xyzzyaace1(xyzzy&
&aaag1(iterm)+xyzzyaaaf47))//':e-n cusp',xyzzyaaej47(xyzzyaaag1(iterm)&
&+xyzzyaaaf47),xyzzyaaab47)
if(xyzzyaaab47/=0)xyzzyaaej47(xyzzyaaag1(iterm)+xyzzyaaaf47)=xyzzyaadr&
&47
enddo
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
xyzzyaaed47=.false.
xyzzyaaee47=.false.
if(.not.xyzzyaaei47(xyzzyaaag1(iterm)+xyzzyaaaf47))then
xyzzyaaca1(:,xyzzyaaag1(iterm)+xyzzyaaaf47)=0.d0
do ispin=1,nspin
if(nele(ispin)==0)cycle
do jspin=ispin,nspin
if(nele(jspin)==0.or.(jspin==ispin.and.nele(ispin)<2))cycle
do xyzzyaacc47=1,xyzzyaaao1(iterm)
if(xyzzyaaay1(ispin,jspin,iterm)/=xyzzyaacd47)cycle
if(xyzzyaaek47(ispin,jspin))then
xyzzyaacc1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)=.true.
endif
enddo
enddo
enddo
else
do ispin=1,nspin
if(nele(ispin)==0)cycle
do jspin=ispin,nspin
if(nele(jspin)==0.or.(jspin==ispin.and.nele(ispin)<2))cycle
do xyzzyaacc47=1,xyzzyaaao1(iterm)
xyzzyaacd47=xyzzyaacf1(xyzzyaacc47,xyzzyaaag1(iterm)+xyzzyaaaf47)
if(xyzzyaaay1(ispin,jspin,iterm)/=xyzzyaacd47)cycle
if(xyzzyaaed47(ispin,jspin))cycle
if(xyzzyaaek47(ispin,jspin).or.xyzzyaaeb47(ispin,jspin))then
xyzzyaaca1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)=0.d0
if(xyzzyaaek47(ispin,jspin))xyzzyaacc1(xyzzyaacd47,xyzzyaaag1(iterm)+x&
&yzzyaaaf47)=.true.
cycle
endif
xyzzyaadi47=0.d0
if(.not.ee_cusp_in_orbital(ispin,jspin))then
xyzzyaadi47=ee_kato_gamma(ispin,jspin,f_zero=f_zero)
if(.not.(xyzzyaaci47(iterm)==f_zero.or.xyzzyaack47(iterm)==f_zero))xyz&
&zyaadi47=0.d0
endif
if(xyzzyaadi47/=xyzzyaaca1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47))&
&call errwarn('READ_GJASTROW','In term #'//trim(i2s(iterm))//': partic&
&le pair '//trim(i2s(ispin))//'-'//trim(i2s(jspin))//' will have the w&
&rong cusp since it falls in the channel whose model is '//trim(xyzzya&
&ace1(xyzzyaaag1(iterm)+xyzzyaaaf47))//'. Will assume this is intended&
&.')
xyzzyaadi47=xyzzyaaca1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)
xyzzyaaed47(ispin,jspin)=abs(xyzzyaadi47)>1.d-11
xyzzyaaed47(jspin,ispin)=abs(xyzzyaadi47)>1.d-11
enddo
enddo
enddo
xyzzyaaeb47=xyzzyaaeb47.or.xyzzyaaed47
endif
if(.not.xyzzyaaej47(xyzzyaaag1(iterm)+xyzzyaaaf47))then
xyzzyaacb1(:,xyzzyaaag1(iterm)+xyzzyaaaf47)=0.d0
do xyzzyaaal47=1,nitot
do ispin=1,nspin
if(nele(ispin)==0)cycle
do xyzzyaacc47=1,xyzzyaaap1(iterm)
xyzzyaacd47=xyzzyaacf1(xyzzyaaao1(iterm)+xyzzyaacc47,xyzzyaaag1(iterm)&
&+xyzzyaaaf47)
if(xyzzyaaaz1(xyzzyaaal47,ispin,iterm)/=xyzzyaacd47)cycle
if(xyzzyaael47(ispin,xyzzyaaal47))then
xyzzyaacd1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)=.true.
endif
enddo
enddo
enddo
else
do xyzzyaaal47=1,nitot
do ispin=1,nspin
if(nele(ispin)==0)cycle
do xyzzyaacc47=1,xyzzyaaap1(iterm)
xyzzyaacd47=xyzzyaacf1(xyzzyaaao1(iterm)+xyzzyaacc47,xyzzyaaag1(iterm)&
&+xyzzyaaaf47)
if(xyzzyaaaz1(xyzzyaaal47,ispin,iterm)/=xyzzyaacd47)cycle
if(xyzzyaaee47(ispin,xyzzyaaal47))cycle
if(xyzzyaael47(ispin,xyzzyaaal47).or.xyzzyaaec47(ispin,xyzzyaaal47))th&
&en
xyzzyaacb1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47)=0.d0
if(xyzzyaael47(ispin,xyzzyaaal47))xyzzyaacd1(xyzzyaacd47,xyzzyaaag1(it&
&erm)+xyzzyaaaf47)=.true.
cycle
endif
if(allocated(zion))then
xyzzyaadi47=en_kato_gamma(ispin,zion(iontype(xyzzyaaal47)))
else
xyzzyaadi47=en_kato_gamma(ispin,0.d0)
endif
if(xyzzyaadi47/=xyzzyaacb1(xyzzyaacd47,xyzzyaaag1(iterm)+xyzzyaaaf47))&
&call errwarn('READ_GJASTROW','In term #'//trim(i2s(iterm))//': partic&
&le-ion pair '//trim(i2s(ispin))//'-n'//trim(i2s(xyzzyaaal47))//' will&
& have the wrong cusp since it falls in the channel whose model is '//&
&trim(xyzzyaace1(xyzzyaaag1(iterm)+xyzzyaaaf47))//'. Will assume this &
&is intended.')
xyzzyaaee47(ispin,xyzzyaaal47)=.true.
enddo
enddo
enddo
xyzzyaaec47=xyzzyaaec47.or.xyzzyaaee47
endif
enddo
call pop_casl_context()
enddo
deallocate(xyzzyaaei47,xyzzyaaej47,xyzzyaaed47,xyzzyaaee47,xyzzyaaek47&
&,xyzzyaael47,xyzzyaaci47,xyzzyaacj47,xyzzyaack47,xyzzyaacl47)
call check_unread_casl('JASTROW',errmsg)
if(len_trim(errmsg)>0)call errstop_master('READ_GJASTROW',errmsg)
allocate(xyzzyaact1(xyzzyaaad1),xyzzyaadb47(xyzzyaaad1),xyzzyaadc47(xy&
&zzyaaad1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','channel_param_equal')
xyzzyaact1=.false.
xyzzyaadb47=0
xyzzyaadc47=1
allocate(xyzzyaadl1(xyzzyaaab1),xyzzyaadm1(xyzzyaaab1),       xyzzyaad&
&o1(xyzzyaaab1),  stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','channel_icons_*')
xyzzyaadl1=1
xyzzyaadm1=1
xyzzyaado1=1
xyzzyaaau47=0
xyzzyaaav47=0
do iterm=1,xyzzyaaaa1
xyzzyaaam47=xyzzyaaag1(iterm)
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
if(xyzzyaaao1(iterm)==0)then
xyzzyaadl1(xyzzyaaam47+xyzzyaaaf47)=1
else
xyzzyaadl1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaaau47+1
xyzzyaaau47=xyzzyaaau47+xyzzyaaao1(iterm)
endif
if(xyzzyaaap1(iterm)==0)then
xyzzyaadm1(xyzzyaaam47+xyzzyaaaf47)=1
else
xyzzyaadm1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaaav47+1
xyzzyaaav47=xyzzyaaav47+xyzzyaaap1(iterm)
endif
enddo
enddo
xyzzyaaao47=0
do iterm=1,xyzzyaaaa1
xyzzyaaam47=xyzzyaaag1(iterm)
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
xyzzyaado1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaaao47+1
xyzzyaaao47=xyzzyaaao47+2*xyzzyaaau1(iterm)*xyzzyaaaw1(iterm)
enddo
enddo
allocate(xyzzyaadn1(max(1,xyzzyaaao47)),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','channel_param_in_eqn')
xyzzyaadn1=0
allocate(xyzzyaacm47(xyzzyaaab1),xyzzyaacn47(xyzzyaaab1),  xyzzyaaco47&
&(xyzzyaaab1),xyzzyaacp47(xyzzyaaab1),xyzzyaacq47(xyzzyaaab1),xyzzyaac&
&r47(xyzzyaaab1),  xyzzyaacs47(xyzzyaaab1),xyzzyaact47(xyzzyaaab1),xyz&
&zyaacu47(xyzzyaaab1),xyzzyaacv47(xyzzyaaab1),  xyzzyaacw47(xyzzyaaab1&
&),xyzzyaacx47(xyzzyaaab1),xyzzyaacy47(xyzzyaaab1),xyzzyaacz47(xyzzyaa&
&ab1),  xyzzyaadd1(xyzzyaaab1),xyzzyaadh1(xyzzyaaab1),xyzzyaacx1(xyzzy&
&aaab1),xyzzyaadb1(xyzzyaaab1),xyzzyaade1(xyzzyaaab1),xyzzyaadi1(xyzzy&
&aaab1),xyzzyaacy1(xyzzyaaab1),xyzzyaadc1(xyzzyaaab1),xyzzyaadf1(xyzzy&
&aaab1),xyzzyaadj1(xyzzyaaab1),xyzzyaacz1(xyzzyaaab1),xyzzyaada47(xyzz&
&yaaab1),xyzzyaadg1(xyzzyaaab1),xyzzyaadk1(xyzzyaaab1),xyzzyaada1(xyzz&
&yaaab1),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','eqn')
xyzzyaacm47=0
xyzzyaaco47=0
xyzzyaacn47=0
xyzzyaacp47=0
xyzzyaacq47=0
xyzzyaacs47=0
xyzzyaacr47=0
xyzzyaact47=0
xyzzyaacu47=0
xyzzyaacw47=0
xyzzyaacv47=0
xyzzyaacx47=0
xyzzyaacy47=0
xyzzyaadd1=0
xyzzyaacz47=0
xyzzyaadh1=0
xyzzyaacx1=0
xyzzyaade1=0
xyzzyaadb1=0
xyzzyaadi1=0
xyzzyaacy1=0
xyzzyaadf1=0
xyzzyaadc1=0
xyzzyaadj1=0
xyzzyaacz1=0
xyzzyaadg1=0
xyzzyaada47=0
xyzzyaadk1=0
xyzzyaada1=0
xyzzyaabx47=0
xyzzyaaby47=0
xyzzyaacv1=0
xyzzyaacw1=0
do iterm=1,xyzzyaaaa1
xyzzyaaam47=xyzzyaaag1(iterm)
xyzzyaaay47=xyzzyaaai1(iterm)
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
xyzzyaadg1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacv1+1
xyzzyaadk1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacw1+1
xyzzyaaai47=xyzzyaaaj1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaaaj47=xyzzyaaai47+1
xyzzyaaak47=xyzzyaaai47+xyzzyaaaw1(iterm)
xyzzyaaaw47=xyzzyaado1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaade1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaabx47+1
xyzzyaadi1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaaby47+1
call xyzzyaahu1(xyzzyaaak1(iterm),xyzzyaaal1(iterm),xyzzyaaau1(iterm),&
&xyzzyaaaw1(iterm),xyzzyaacg1(xyzzyaaay47),xyzzyaacf1(1,xyzzyaaam47+xy&
&zzyaaaf47),xyzzyaabc1(iterm),xyzzyaact1(xyzzyaaaj47),xyzzyaadb47(xyzz&
&yaaaj47),xyzzyaadc47(xyzzyaaaj47),xyzzyaacq47(xyzzyaaam47+xyzzyaaaf47&
&))
xyzzyaaco47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaabx47+1
xyzzyaacp47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaaby47+1
xyzzyaacm47(xyzzyaaam47+xyzzyaaaf47)=count(xyzzyaact1(xyzzyaaaj47:xyzz&
&yaaak47))
if(dble(xyzzyaacm47(xyzzyaaam47+xyzzyaaaf47))*dble(xyzzyaaaw1(iterm))>&
&dble(huge(1)))call errstop_master('READ_GJASTROW','Size of parameter-&
&removal constraints matrix for channel '//trim(i2s(xyzzyaaaf47))//' i&
&n term '//trim(i2s(iterm))//' exceeds the largest representable integ&
&er.  Try using smaller expansion orders in this term.')
xyzzyaacn47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacm47(xyzzyaaam47+xyzzyaaaf&
&47)*xyzzyaaaw1(iterm)
xyzzyaabx47=xyzzyaabx47+xyzzyaacm47(xyzzyaaam47+xyzzyaaaf47)
if(dble(xyzzyaaby47)+dble(xyzzyaacn47(xyzzyaaam47+xyzzyaaaf47))>dble(h&
&uge(1)))call errstop_master('READ_GJASTROW','Size of matrix containin&
&g invariant constraints up to the parameter-removal constraints in ch&
&annel '//trim(i2s(xyzzyaaaf47))//' in term '//trim(i2s(iterm))//' exc&
&eeds the largest representable integer.  Try using smaller expansion &
&orders or fewer Jastrow terms overall.')
xyzzyaaby47=xyzzyaaby47+xyzzyaacn47(xyzzyaaam47+xyzzyaaaf47)
xyzzyaacs47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaabx47+1
xyzzyaact47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaaby47+1
if(dble(xyzzyaacq47(xyzzyaaam47+xyzzyaaaf47))*dble(xyzzyaaaw1(iterm))>&
&dble(huge(1)))call errstop_master('READ_GJASTROW','Size of symmetry c&
&onstraints matrix for channel '//trim(i2s(xyzzyaaaf47))//' in term '/&
&/trim(i2s(iterm))//' exceeds the largest representable integer.  Try &
&using smaller expansion orders in this term.')
xyzzyaacr47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacq47(xyzzyaaam47+xyzzyaaaf&
&47)*xyzzyaaaw1(iterm)
xyzzyaabx47=xyzzyaabx47+xyzzyaacq47(xyzzyaaam47+xyzzyaaaf47)
if(dble(xyzzyaaby47)+dble(xyzzyaacr47(xyzzyaaam47+xyzzyaaaf47))>dble(h&
&uge(1)))call errstop_master('READ_GJASTROW','Size of matrix containin&
&g invariant constraints up to the symmetry constraints in channel '//&
&trim(i2s(xyzzyaaaf47))//' in term '//trim(i2s(iterm))//' exceeds the &
&largest representable integer.  Try using smaller expansion orders or&
& fewer Jastrow terms overall.')
xyzzyaaby47=xyzzyaaby47+xyzzyaacr47(xyzzyaaam47+xyzzyaaaf47)
xyzzyaacw47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaabx47+1
xyzzyaacx47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaaby47+1
call xyzzyaahw1(xyzzyaaao1(iterm),xyzzyaaap1(iterm),xyzzyaaau1(iterm),&
&xyzzyaaaw1(iterm),xyzzyaacg1(xyzzyaaay47),                     xyzzya&
&acf1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaact1(xyzzyaaaj47),      xyzzyaa&
&eh47(:,iterm),xyzzyaacu47(xyzzyaaam47+xyzzyaaaf47))
if(dble(xyzzyaacu47(xyzzyaaam47+xyzzyaaaf47))*dble(xyzzyaaaw1(iterm))>&
&dble(huge(1)))call errstop_master('READ_GJASTROW','Size of user const&
&raints matrix for channel '//trim(i2s(xyzzyaaaf47))//' in term '//tri&
&m(i2s(iterm))//' exceeds the largest representable integer.  Try usin&
&g smaller expansion orders in this term.')
xyzzyaacv47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacu47(xyzzyaaam47+xyzzyaaaf&
&47)*xyzzyaaaw1(iterm)
xyzzyaabx47=xyzzyaabx47+xyzzyaacu47(xyzzyaaam47+xyzzyaaaf47)
if(dble(xyzzyaaby47)+dble(xyzzyaacv47(xyzzyaaam47+xyzzyaaaf47))>dble(h&
&uge(1)))call errstop_master('READ_GJASTROW','Size of matrix containin&
&g invariant constraints up to the user constraints in channel '//trim&
&(i2s(xyzzyaaaf47))//' in term '//trim(i2s(iterm))//' exceeds the larg&
&est representable integer.  Try using smaller expansion orders or few&
&er Jastrow terms overall.')
xyzzyaaby47=xyzzyaaby47+xyzzyaacv47(xyzzyaaam47+xyzzyaaaf47)
xyzzyaacx1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacm47(xyzzyaaam47+xyzzyaaaf4&
&7)+xyzzyaacq47(xyzzyaaam47+xyzzyaaaf47)+xyzzyaacu47(xyzzyaaam47+xyzzy&
&aaaf47)
xyzzyaadb1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacn47(xyzzyaaam47+xyzzyaaaf4&
&7)+xyzzyaacr47(xyzzyaaam47+xyzzyaaaf47)+xyzzyaacv47(xyzzyaaam47+xyzzy&
&aaaf47)
xyzzyaacz1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacz1(xyzzyaaam47+xyzzyaaaf47&
&)+xyzzyaacx1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaada47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaada47(xyzzyaaam47+xyzzyaaaf&
&47)+xyzzyaadb1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaacv1=xyzzyaacv1+xyzzyaacx1(xyzzyaaam47+xyzzyaaaf47)
if(dble(xyzzyaacw1)+dble(xyzzyaadb1(xyzzyaaam47+xyzzyaaaf47))>dble(hug&
&e(1)))call errstop_master('READ_GJASTROW','Size of matrix containing &
&all constraints up to the invariant constraints in channel '//trim(i2&
&s(xyzzyaaaf47))//' in term '//trim(i2s(iterm))//' exceeds the largest&
& representable integer.  Try using smaller expansion orders or fewer &
&Jastrow terms overall.')
xyzzyaacw1=xyzzyaacw1+xyzzyaadb1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaadf1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacv1+1
xyzzyaadj1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacw1+1
xyzzyaadd1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacv1+1
xyzzyaadh1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacw1+1
call xyzzyaahz1(xyzzyaaak1(iterm),xyzzyaaal1(iterm),xyzzyaaao1(iterm),&
&xyzzyaaap1(iterm),xyzzyaaau1(iterm),xyzzyaaam1(iterm),xyzzyaaan1(iter&
&m),xyzzyaaaw1(iterm),xyzzyaacg1(xyzzyaaay47),xyzzyaacf1(1,xyzzyaaam47&
&+xyzzyaaaf47),xyzzyaaaq1(iterm),xyzzyaaar1(iterm),xyzzyaaas1(iterm),x&
&yzzyaaat1(iterm),xyzzyaabf1(iterm),xyzzyaabh1(iterm),xyzzyaabg1(iterm&
&),xyzzyaabi1(iterm),xyzzyaabr1(iterm),xyzzyaabs1(iterm),xyzzyaabx1(1,&
&1,iterm),xyzzyaaby1(1,1,iterm),xyzzyaact1(xyzzyaaaj47),xyzzyaaca1(1,x&
&yzzyaaam47+xyzzyaaaf47),xyzzyaacb1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaa&
&cc1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaacd1(1,xyzzyaaam47+xyzzyaaaf47),&
&xyzzyaadn1(xyzzyaaaw47),xyzzyaacy47(xyzzyaaam47+xyzzyaaaf47))
if(dble(xyzzyaacy47(xyzzyaaam47+xyzzyaaaf47))*dble(xyzzyaaaw1(iterm))>&
&dble(huge(1)))call errstop_master('READ_GJASTROW','Size of cusp const&
&raints matrix for channel '//trim(i2s(xyzzyaaaf47))//' in term '//tri&
&m(i2s(iterm))//' exceeds the largest representable integer.  Try usin&
&g smaller expansion orders in this term.')
xyzzyaacz47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacy47(xyzzyaaam47+xyzzyaaaf&
&47)*xyzzyaaaw1(iterm)
xyzzyaacy1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacy47(xyzzyaaam47+xyzzyaaaf4&
&7)
xyzzyaadc1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacz47(xyzzyaaam47+xyzzyaaaf4&
&7)
xyzzyaacz1(xyzzyaaam47+xyzzyaaaf47)=xyzzyaacz1(xyzzyaaam47+xyzzyaaaf47&
&)+xyzzyaacy1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaada47(xyzzyaaam47+xyzzyaaaf47)=xyzzyaada47(xyzzyaaam47+xyzzyaaaf&
&47)+xyzzyaadc1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaacv1=xyzzyaacv1+xyzzyaacy1(xyzzyaaam47+xyzzyaaaf47)
if(dble(xyzzyaacw1)+dble(xyzzyaadc1(xyzzyaaam47+xyzzyaaaf47))>dble(hug&
&e(1)))call errstop_master('READ_GJASTROW','Size of matrix containing &
&all constraints up to the cusp constraints in channel '//trim(i2s(xyz&
&zyaaaf47))//' in term '//trim(i2s(iterm))//' exceeds the largest repr&
&esentable integer.  Try using smaller expansion orders or fewer Jastr&
&ow terms overall.')
xyzzyaacw1=xyzzyaacw1+xyzzyaadc1(xyzzyaaam47+xyzzyaaaf47)
enddo
enddo
allocate(xyzzyaadr1(xyzzyaacw1),xyzzyaads1(xyzzyaaby47),eqn_rhs(xyzzya&
&acv1),xyzzyaadt1(xyzzyaabx47),stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','constraint matrices')
xyzzyaadr1=0.d0
xyzzyaads1=0.d0
eqn_rhs=0.d0
xyzzyaadt1=0.d0
allocate(xyzzyaadp1(xyzzyaacv1),xyzzyaadq1(xyzzyaacv1),stat=xyzzyaaaa4&
&7)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','eqn_solves_for')
xyzzyaadp1=0
xyzzyaadq1=0
do iterm=1,xyzzyaaaa1
call push_casl_context('TERM '//trim(i2s(iterm)))
if(am_master)then
call wout('Term '//trim(i2s(iterm)))
call wout('------')
call wout('Rank(e)             : '//trim(i2s(xyzzyaaak1(iterm))))
call wout('Rank(n)             : '//trim(i2s(xyzzyaaal1(iterm))))
endif
if(xyzzyaaao1(iterm)>0)then
if(am_master)call wout('e-e basis:')
call rewrite_eebasis('e-e basis',xyzzyaaaq1(iterm),xyzzyaaam1(iterm),x&
&yzzyaahd1)
if(am_master)call wout('e-e cutoff:')
call rewrite_eebasis('e-e cutoff',xyzzyaaar1(iterm),1,xyzzyaahd1)
endif
if(xyzzyaaap1(iterm)>0)then
if(am_master)call wout('e-n basis:')
call rewrite_enbasis('e-n basis',xyzzyaaas1(iterm),xyzzyaaan1(iterm),x&
&yzzyaahd1)
if(am_master)call wout('e-n cutoff:')
call rewrite_enbasis('e-n cutoff',xyzzyaaat1(iterm),1,xyzzyaahd1)
endif
if(am_master)then
call wout('Linear parameters:')
call wout(' Number of linear parameters per channel: '//trim(i2s(xyzzy&
&aaaw1(iterm))))
endif
xyzzyaaay47=xyzzyaaai1(iterm)
xyzzyaaam47=xyzzyaaag1(iterm)
do xyzzyaaaf47=1,xyzzyaaav1(iterm)
if(am_master)call wout(' Channel '//trim(xyzzyaace1(xyzzyaaam47+xyzzya&
&aaf47))//':')
xyzzyaaai47=xyzzyaaaj1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaaaj47=xyzzyaaai47+1
xyzzyaaak47=xyzzyaaai47+xyzzyaaaw1(iterm)
xyzzyaaag47=xyzzyaaah1(iterm)+1
xyzzyaaah47=xyzzyaaah1(iterm)+xyzzyaaaw1(iterm)
xyzzyaaaq47=xyzzyaacz1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaaar47=xyzzyaadg1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaaas47=xyzzyaadk1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaabh47=xyzzyaacm47(xyzzyaaam47+xyzzyaaaf47)
if(xyzzyaabh47>0)then
xyzzyaabi47=xyzzyaaco47(xyzzyaaam47+xyzzyaaaf47)
xyzzyaabj47=xyzzyaacp47(xyzzyaaam47+xyzzyaaaf47)
call xyzzyaaht1(xyzzyaaaw1(iterm),xyzzyaact1(xyzzyaaaj47),xyzzyaads1(x&
&yzzyaabj47),xyzzyaadt1(xyzzyaabi47))
endif
xyzzyaabk47=xyzzyaacq47(xyzzyaaam47+xyzzyaaaf47)
if(xyzzyaabk47>0)then
xyzzyaabl47=xyzzyaacs47(xyzzyaaam47+xyzzyaaaf47)
xyzzyaabm47=xyzzyaact47(xyzzyaaam47+xyzzyaaaf47)
call xyzzyaahv1(xyzzyaaaw1(iterm),xyzzyaadb47(xyzzyaaaj47),xyzzyaadc47&
&(xyzzyaaaj47),xyzzyaads1(xyzzyaabm47),xyzzyaadt1(xyzzyaabl47))
endif
xyzzyaabn47=xyzzyaacu47(xyzzyaaam47+xyzzyaaaf47)
if(xyzzyaabn47>0)then
xyzzyaabo47=xyzzyaacw47(xyzzyaaam47+xyzzyaaaf47)
xyzzyaabp47=xyzzyaacx47(xyzzyaaam47+xyzzyaaaf47)
call xyzzyaahx1(xyzzyaaao1(iterm),xyzzyaaap1(iterm),xyzzyaaau1(iterm),&
&xyzzyaaaw1(iterm),xyzzyaacg1(xyzzyaaay47),xyzzyaacf1(1,xyzzyaaam47+xy&
&zzyaaaf47),xyzzyaact1(xyzzyaaaj47),xyzzyaaeh47(:,iterm),xyzzyaads1(xy&
&zzyaabp47),xyzzyaadt1(xyzzyaabo47))
endif
xyzzyaabt47=xyzzyaacx1(xyzzyaaam47+xyzzyaaaf47)
if(xyzzyaabt47>0)then
xyzzyaabu47=xyzzyaadb1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaabv47=xyzzyaade1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaabw47=xyzzyaadi1(xyzzyaaam47+xyzzyaaaf47)
call dcopy(xyzzyaabu47,xyzzyaads1(xyzzyaabw47),1,xyzzyaadr1(xyzzyaaas4&
&7),1)
call dcopy(xyzzyaabt47,xyzzyaadt1(xyzzyaabv47),1,eqn_rhs(xyzzyaaar47),&
&1)
endif
xyzzyaabq47=xyzzyaacy47(xyzzyaaam47+xyzzyaaaf47)
if(xyzzyaabq47>0)then
xyzzyaabr47=xyzzyaadd1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaabs47=xyzzyaadh1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaadr1(xyzzyaabs47:xyzzyaabs47+xyzzyaacy1(xyzzyaaam47+xyzzyaaaf47)&
&*xyzzyaaaw1(iterm)-1)=0.d0
eqn_rhs(xyzzyaabr47:xyzzyaabr47+xyzzyaacy1(xyzzyaaam47+xyzzyaaaf47)-1)&
&=0.d0
xyzzyaaaw47=xyzzyaado1(xyzzyaaam47+xyzzyaaaf47)
call xyzzyaaia1(xyzzyaaao1(iterm),xyzzyaaap1(iterm),xyzzyaaau1(iterm),&
&xyzzyaaam1(iterm),xyzzyaaan1(iterm),xyzzyaaaw1(iterm),xyzzyaacg1(xyzz&
&yaaay47),xyzzyaacf1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaaaq1(iterm),xyzz&
&yaaar1(iterm),xyzzyaaas1(iterm),xyzzyaaat1(iterm),xyzzyaabf1(iterm),x&
&yzzyaabh1(iterm),xyzzyaabg1(iterm),xyzzyaabi1(iterm),xyzzyaabr1(iterm&
&),xyzzyaabs1(iterm),xyzzyaadn1(xyzzyaaaw47),xyzzyaact1(xyzzyaaaj47),.&
&true.,xyzzyaaca1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaacb1(1,xyzzyaaam47+&
&xyzzyaaaf47),xyzzyaacc1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaacd1(1,xyzzy&
&aaam47+xyzzyaaaf47),xyzzyaadr1(xyzzyaabs47),eqn_rhs(xyzzyaabr47))
endif
if(xyzzyaaaq47>0)then
if(xyzzyaahe1)then
if(am_master)then
call wout(repeat('*',78))
if(xyzzyaaaq47>xyzzyaabt47)then
call wout('DEBUG_INITIAL_CONSTRAINTS - before Gaussian elim. (FAKE):')
else
call wout('DEBUG_INITIAL_CONSTRAINTS - before Gaussian elim.:')
endif
call wout(repeat('*',78))
call wout('Term '//trim(i2s(iterm))//', channel '//trim(i2s(xyzzyaaaf4&
&7))//':')
call print_eqns(xyzzyaaaw1(iterm),xyzzyaaaq47,xyzzyaadr1(xyzzyaaas47),&
&eqn_rhs(xyzzyaaar47),xyzzyaach1(xyzzyaaag47:xyzzyaaah47))
call wout(repeat('*',78))
endif
endif
call gaussian_elimination(xyzzyaaaw1(iterm),xyzzyaaaq47,xyzzyaadr1(xyz&
&zyaaas47),eqn_rhs(xyzzyaaar47),xyzzyaada1(xyzzyaaam47+xyzzyaaaf47),xy&
&zzyaadp1(xyzzyaaar47),xyzzyaadq1(xyzzyaaar47),xyzzyaaab47)
if(xyzzyaahe1)then
if(am_master)then
call wout(repeat('*',78))
if(xyzzyaaaq47>xyzzyaabt47)then
call wout('DEBUG_INITIAL_CONSTRAINTS - after Gaussian elim. (FAKE):')
else
call wout('DEBUG_INITIAL_CONSTRAINTS - after Gaussian elim.:')
endif
call wout(repeat('*',78))
call wout('Term '//trim(i2s(iterm))//', channel '//trim(i2s(xyzzyaaaf4&
&7))//':')
call print_eqns(xyzzyaaaw1(iterm),xyzzyaaaq47,xyzzyaadr1(xyzzyaaas47),&
&eqn_rhs(xyzzyaaar47),xyzzyaach1(xyzzyaaag47:xyzzyaaah47),xyzzyaada1(x&
&yzzyaaam47+xyzzyaaaf47),xyzzyaadq1(xyzzyaaar47),xyzzyaadp1(xyzzyaaar4&
&7))
call wout(repeat('*',78))
endif
endif
if(xyzzyaaab47/=0)call errstop_master('READ_GJASTROW','Constraints for&
& term #'//trim(i2s(iterm))//', channel #'//trim(i2s(xyzzyaaaf47))//' &
&result in a linear system of equations with no solution.')
xyzzyaabz47=xyzzyaacy1(xyzzyaaam47+xyzzyaaaf47)
if(xyzzyaabz47>0)then
if(xyzzyaabt47>0)then
call dcopy(xyzzyaabu47,xyzzyaads1(xyzzyaabw47),1,xyzzyaadr1(xyzzyaaas4&
&7),1)
call dcopy(xyzzyaabt47,xyzzyaadt1(xyzzyaabv47),1,eqn_rhs(xyzzyaaar47),&
&1)
endif
xyzzyaadr1(xyzzyaabs47:xyzzyaabs47+xyzzyaacy1(xyzzyaaam47+xyzzyaaaf47)&
&*xyzzyaaaw1(iterm)-1)=0.d0
eqn_rhs(xyzzyaabr47:xyzzyaabr47+xyzzyaacy1(xyzzyaaam47+xyzzyaaaf47)-1)&
&=0.d0
call xyzzyaaia1(xyzzyaaao1(iterm),xyzzyaaap1(iterm),xyzzyaaau1(iterm),&
&xyzzyaaam1(iterm),xyzzyaaan1(iterm),xyzzyaaaw1(iterm),xyzzyaacg1(xyzz&
&yaaay47),xyzzyaacf1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaaaq1(iterm),xyzz&
&yaaar1(iterm),xyzzyaaas1(iterm),xyzzyaaat1(iterm),xyzzyaabf1(iterm),x&
&yzzyaabh1(iterm),xyzzyaabg1(iterm),xyzzyaabi1(iterm),xyzzyaabr1(iterm&
&),xyzzyaabs1(iterm),xyzzyaadn1(xyzzyaaaw47),xyzzyaact1(xyzzyaaaj47),.&
&false.,xyzzyaaca1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaacb1(1,xyzzyaaam47&
&+xyzzyaaaf47),xyzzyaacc1(1,xyzzyaaam47+xyzzyaaaf47),xyzzyaacd1(1,xyzz&
&yaaam47+xyzzyaaaf47),xyzzyaadr1(xyzzyaabs47),eqn_rhs(xyzzyaabr47))
if(xyzzyaahe1)then
if(am_master)then
call wout(repeat('*',78))
call wout('DEBUG_INITIAL_CONSTRAINTS - before Gaussian elim.:')
call wout(repeat('*',78))
call wout('Term '//trim(i2s(iterm))//', channel '//trim(i2s(xyzzyaaaf4&
&7))//':')
call print_eqns(xyzzyaaaw1(iterm),xyzzyaaaq47,xyzzyaadr1(xyzzyaaas47),&
&eqn_rhs(xyzzyaaar47),xyzzyaach1(xyzzyaaag47:xyzzyaaah47))
call wout(repeat('*',78))
endif
endif
call redo_gaussian_elimination(xyzzyaaaw1(iterm),xyzzyaaaq47,xyzzyaadr&
&1(xyzzyaaas47),eqn_rhs(xyzzyaaar47),xyzzyaada1(xyzzyaaam47+xyzzyaaaf4&
&7),xyzzyaadp1(xyzzyaaar47),xyzzyaadq1(xyzzyaaar47))
if(xyzzyaahe1)then
if(am_master)then
call wout(repeat('*',78))
call wout('DEBUG_INITIAL_CONSTRAINTS - after Gaussian elim.:')
call wout(repeat('*',78))
call wout('Term '//trim(i2s(iterm))//', channel '//trim(i2s(xyzzyaaaf4&
&7))//':')
call print_eqns(xyzzyaaaw1(iterm),xyzzyaaaq47,xyzzyaadr1(xyzzyaaas47),&
&eqn_rhs(xyzzyaaar47),xyzzyaach1(xyzzyaaag47:xyzzyaaah47),xyzzyaada1(x&
&yzzyaaam47+xyzzyaaaf47),xyzzyaadq1(xyzzyaaar47),xyzzyaadp1(xyzzyaaar4&
&7))
call wout(repeat('*',78))
endif
endif
endif
if(am_master)then
word1='equations'
word2='parameters'
explain=trim(i2s(xyzzyaabh47))//'r+'//trim(i2s(xyzzyaabk47))//'s+'//tr&
&im(i2s(xyzzyaabn47))//'u+'//trim(i2s(xyzzyaabq47))//'c'
if(xyzzyaaaq47==1)word1='equation'
if(xyzzyaada1(xyzzyaaam47+xyzzyaaaf47)==1)word2='parameter'
call wout('  Constraints: '//trim(i2s(xyzzyaaaq47))//' '//trim(word1)/&
&/' ('//trim(explain)//'), '//trim(i2s(xyzzyaada1(xyzzyaaam47+xyzzyaaa&
&f47)))//' '//trim(word2)//' determined')
endif
allocate(xyzzyaadj47(xyzzyaaaw1(iterm)),xyzzyaadk47(xyzzyaaaw1(iterm))&
&,stat=xyzzyaaaa47)
call check_alloc(xyzzyaaaa47,'READ_GJASTROW','tmp_param')
xyzzyaadj47=param(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))
xyzzyaadk47=xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))
where(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))==pflag_d&
&et)xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))=pflag_opt
do xyzzyaaao47=1,xyzzyaada1(xyzzyaaam47+xyzzyaaaf47)
xyzzyaaci1(xyzzyaaai47+xyzzyaadp1(xyzzyaaar47+xyzzyaaao47-1))=pflag_de&
&t
enddo
call construct_forced_zero_mask(xyzzyaaaw1(iterm),xyzzyaacz1(xyzzyaaam&
&47+xyzzyaaaf47),xyzzyaada1(xyzzyaaam47+xyzzyaaaf47),xyzzyaadr1(xyzzya&
&aas47),eqn_rhs(xyzzyaaar47),xyzzyaadp1(xyzzyaaar47),xyzzyaadq1(xyzzya&
&aar47),xyzzyaacu1(xyzzyaaai47+1))
call backwards_substitution(xyzzyaaaw1(iterm),xyzzyaacz1(xyzzyaaam47+x&
&yzzyaaaf47),xyzzyaada1(xyzzyaaam47+xyzzyaaaf47),xyzzyaadr1(xyzzyaaas4&
&7),eqn_rhs(xyzzyaaar47),xyzzyaadp1(xyzzyaaar47),xyzzyaadq1(xyzzyaaar4&
&7),xyzzyaacu1(xyzzyaaai47+1),param(xyzzyaaai47+1))
xyzzyaaao47=count(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iter&
&m))==pflag_unset)
if(xyzzyaaao47>0)then
if(am_master)then
word1='parameters'
if(xyzzyaaao47==1)word1='parameter'
call wout('  Info: '//trim(i2s(xyzzyaaao47))//' unspecified '//trim(wo&
&rd1)//' added as optimizable')
endif
where(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))==pflag_u&
&nset)xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))=pflag_op&
&t
endif
xyzzyaaao47=count(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iter&
&m))==pflag_opt.and.xyzzyaadk47(1:xyzzyaaaw1(iterm))==pflag_det)
if(xyzzyaaao47>0)then
if(am_master)then
word1='parameters'
if(xyzzyaaao47==1)word1='parameter'
call wout('  Info: '//trim(i2s(xyzzyaaao47))//' determined '//trim(wor&
&d1)//' reclassified as optimizable')
endif
endif
call get_all_index_properties(xyzzyaaao1(iterm),xyzzyaaap1(iterm),xyzz&
&yaaau1(iterm),xyzzyaaaw1(iterm),xyzzyaacg1(xyzzyaaay47),xyzzyaacu1(xy&
&zzyaaaj47),xyzzyaaai47,xyzzyaabz1(xyzzyaaam47+xyzzyaaaf47),xyzzyaacj1&
&(xyzzyaaaj47),xyzzyaack1(xyzzyaaaj47),xyzzyaacl1(xyzzyaaaj47),xyzzyaa&
&cm1(xyzzyaaaj47),xyzzyaacn1(xyzzyaaaj47))
xyzzyaaao47=count(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iter&
&m))==pflag_det.and.xyzzyaadk47(1:xyzzyaaaw1(iterm))==pflag_opt)
if(xyzzyaaao47>0)then
if(am_master)then
word1='parameters'
if(xyzzyaaao47==1)word1='parameter'
call wout('  Info: '//trim(i2s(xyzzyaaao47))//' optimizable '//trim(wo&
&rd1)//' reclassified as determined')
endif
xyzzyaaap47=count(abs(param(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm&
&))-xyzzyaadj47(1:xyzzyaaaw1(iterm)))>max(xyzzyaadm47,xyzzyaadl47*abs(&
&param(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))),xyzzyaadl47*abs(x&
&yzzyaadj47(1:xyzzyaaaw1(iterm)))).and.xyzzyaaci1(xyzzyaaai47+1:xyzzya&
&aai47+xyzzyaaaw1(iterm))==pflag_det.and.xyzzyaadk47(1:xyzzyaaaw1(iter&
&m))==pflag_opt)
if(xyzzyaaap47>0)then
if(am_master)then
do xyzzyaaao47=1,xyzzyaaaw1(iterm)
if(xyzzyaaci1(xyzzyaaai47+xyzzyaaao47)==pflag_det.and.xyzzyaadk47(xyzz&
&yaaao47)==pflag_opt.and.abs(param(xyzzyaaai47+xyzzyaaao47)-xyzzyaadj4&
&7(xyzzyaaao47))>max(xyzzyaadm47,xyzzyaadl47*abs(param(xyzzyaaai47+xyz&
&zyaaao47)),xyzzyaadl47*abs(xyzzyaadj47(xyzzyaaao47))))then
write(tmpr,'(2x,a,es16.8,a,es16.8)')trim(xyzzyaach1(xyzzyaaag47-1+xyzz&
&yaaao47))//' changed from ',xyzzyaadj47(xyzzyaaao47),' to ',param(xyz&
&zyaaai47+xyzzyaaao47)
call wout(tmpr)
endif
enddo
endif
word2='differ'
if(xyzzyaaap47==1)word2='differs'
if(xyzzyaahm1)then
do xyzzyaaao47=1,xyzzyaaaw1(iterm)
if(xyzzyaaci1(xyzzyaaai47+xyzzyaaao47)==pflag_det.and.xyzzyaadk47(xyzz&
&yaaao47)==pflag_opt.and.abs(param(xyzzyaaai47+xyzzyaaao47)-xyzzyaadj4&
&7(xyzzyaaao47))>max(xyzzyaadm47,xyzzyaadl47*abs(param(xyzzyaaai47+xyz&
&zyaaao47)),xyzzyaadl47*abs(xyzzyaadj47(xyzzyaaao47))))param(xyzzyaaai&
&47+xyzzyaaao47)=xyzzyaadj47(xyzzyaaao47)
enddo
call errwarn('READ_GJASTROW',trim(i2s(xyzzyaaap47))//' of the reclassi&
&fied parameters '//trim(word2)//' from input value by more than the r&
&elative and/or absolute tolerance.  Assuming input is correct and con&
&tinuing.')
else
call errstop_master('READ_GJASTROW',trim(i2s(xyzzyaaap47))//' of the r&
&eclassified parameters '//trim(word2)//' from input value by more tha&
&n the relative and/or absolute tolerance.')
endif
endif
endif
xyzzyaaao47=count(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iter&
&m))==pflag_det.and.xyzzyaadk47(1:xyzzyaaaw1(iterm))==pflag_fix)
if(xyzzyaaao47>0)then
if(am_master)then
word1='parameters'
if(xyzzyaaao47==1)word1='parameter'
call wout('  Info: '//trim(i2s(xyzzyaaao47))//' fixed '//trim(word1)//&
&' reclassified as determined')
endif
xyzzyaaap47=count(abs(param(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm&
&))-xyzzyaadj47(1:xyzzyaaaw1(iterm)))>max(xyzzyaadm47,xyzzyaadl47*abs(&
&param(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))),xyzzyaadl47*abs(x&
&yzzyaadj47(1:xyzzyaaaw1(iterm)))).and.xyzzyaaci1(xyzzyaaai47+1:xyzzya&
&aai47+xyzzyaaaw1(iterm))==pflag_det.and.xyzzyaadk47(1:xyzzyaaaw1(iter&
&m))==pflag_fix)
if(xyzzyaaap47>0)then
if(am_master)then
do xyzzyaaao47=1,xyzzyaaaw1(iterm)
if(xyzzyaaci1(xyzzyaaai47+xyzzyaaao47)==pflag_det.and.xyzzyaadk47(xyzz&
&yaaao47)==pflag_fix.and.abs(param(xyzzyaaai47+xyzzyaaao47)-xyzzyaadj4&
&7(xyzzyaaao47))>max(xyzzyaadm47,xyzzyaadl47*abs(param(xyzzyaaai47+xyz&
&zyaaao47)),xyzzyaadl47*abs(xyzzyaadj47(xyzzyaaao47))))then
write(tmpr,'(2x,a,es16.8,a,es16.8)')trim(xyzzyaach1(xyzzyaaag47-1+xyzz&
&yaaao47))//' changed from ',xyzzyaadj47(xyzzyaaao47),' to ',param(xyz&
&zyaaai47+xyzzyaaao47)
call wout(tmpr)
endif
enddo
endif
word2='differ'
if(xyzzyaaap47==1)word2='differs'
if(xyzzyaahm1)then
do xyzzyaaao47=1,xyzzyaaaw1(iterm)
if(xyzzyaaci1(xyzzyaaai47+xyzzyaaao47)==pflag_det.and.xyzzyaadk47(xyzz&
&yaaao47)==pflag_fix.and.abs(param(xyzzyaaai47+xyzzyaaao47)-xyzzyaadj4&
&7(xyzzyaaao47))>max(xyzzyaadm47,xyzzyaadl47*abs(param(xyzzyaaai47+xyz&
&zyaaao47)),xyzzyaadl47*abs(xyzzyaadj47(xyzzyaaao47))))param(xyzzyaaai&
&47+xyzzyaaao47)=xyzzyaadj47(xyzzyaaao47)
enddo
call errwarn('READ_GJASTROW',trim(i2s(xyzzyaaap47))//' of the reclassi&
&fied parameters '//trim(word2)//' from input value by more than the r&
&elative and/or absolute tolerance.  Assuming input is correct and con&
&tinuing.')
else
call errstop_master('READ_GJASTROW',trim(i2s(xyzzyaaap47))//' of the r&
&eclassified parameters '//trim(word2)//' from input value by more tha&
&n the relative and/or absolute tolerance.')
endif
endif
endif
xyzzyaaao47=count(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iter&
&m))==pflag_det.and.xyzzyaadk47(1:xyzzyaaaw1(iterm))==pflag_det)
if(xyzzyaaao47>0)then
xyzzyaaap47=count(abs(param(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm&
&))-xyzzyaadj47(1:xyzzyaaaw1(iterm)))>max(xyzzyaadm47,xyzzyaadl47*abs(&
&param(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))),xyzzyaadl47*abs(x&
&yzzyaadj47(1:xyzzyaaaw1(iterm)))).and.xyzzyaaci1(xyzzyaaai47+1:xyzzya&
&aai47+xyzzyaaaw1(iterm))==pflag_det.and.xyzzyaadk47(1:xyzzyaaaw1(iter&
&m))==pflag_det)
if(xyzzyaaap47>0)then
if(am_master)then
do xyzzyaaao47=1,xyzzyaaaw1(iterm)
if(xyzzyaaci1(xyzzyaaai47+xyzzyaaao47)==pflag_det.and.xyzzyaadk47(xyzz&
&yaaao47)==pflag_det.and.abs(param(xyzzyaaai47+xyzzyaaao47)-xyzzyaadj4&
&7(xyzzyaaao47))>max(xyzzyaadm47,xyzzyaadl47*abs(param(xyzzyaaai47+xyz&
&zyaaao47)),xyzzyaadl47*abs(xyzzyaadj47(xyzzyaaao47))))then
write(tmpr,'(2x,a,es16.8,a,es16.8)')trim(xyzzyaach1(xyzzyaaag47-1+xyzz&
&yaaao47))//' changed from ',xyzzyaadj47(xyzzyaaao47),' to ',param(xyz&
&zyaaai47+xyzzyaaao47)
call wout(tmpr)
endif
enddo
endif
word2='differ'
if(xyzzyaaap47==1)word2='differs'
if(xyzzyaahm1)then
do xyzzyaaao47=1,xyzzyaaaw1(iterm)
if(xyzzyaaci1(xyzzyaaai47+xyzzyaaao47)==pflag_det.and.xyzzyaadk47(xyzz&
&yaaao47)==pflag_det.and.abs(param(xyzzyaaai47+xyzzyaaao47)-xyzzyaadj4&
&7(xyzzyaaao47))>max(xyzzyaadm47,xyzzyaadl47*abs(param(xyzzyaaai47+xyz&
&zyaaao47)),xyzzyaadl47*abs(xyzzyaadj47(xyzzyaaao47))))param(xyzzyaaai&
&47+xyzzyaaao47)=xyzzyaadj47(xyzzyaaao47)
enddo
call errwarn('READ_GJASTROW',trim(i2s(xyzzyaaap47))//' of the determin&
&ed parameters '//trim(word2)//' from input value by more than the rel&
&ative and/or absolute tolerance.  Assuming input is correct and conti&
&nuing.')
else
call errstop_master('READ_GJASTROW',trim(i2s(xyzzyaaap47))//' of the d&
&etermined parameters '//trim(word2)//' from input value by more than &
&the relative and/or absolute tolerance.')
endif
endif
endif
deallocate(xyzzyaadj47,xyzzyaadk47)
endif
where(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))==pflag_u&
&nset)xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))=pflag_op&
&t
if(all(xyzzyaaci1(xyzzyaaai47+1:xyzzyaaai47+xyzzyaaaw1(iterm))==pflag_&
&det))empty_jastrow=.false.
call rewrite_gparam_channel('Linear parameters:Channel '//         tri&
&m(xyzzyaace1(xyzzyaaam47+xyzzyaaaf47)),xyzzyaaaw1(iterm),xyzzyaach1(x&
&yzzyaaag47:xyzzyaaah47),param(xyzzyaaaj47:xyzzyaaak47),xyzzyaaci1(xyz&
&zyaaaj47:xyzzyaaak47),xyzzyaacq1(xyzzyaaaj47:xyzzyaaak47),xyzzyaaco1(&
&xyzzyaaaj47:xyzzyaaak47),xyzzyaacr1(xyzzyaaaj47:xyzzyaaak47),xyzzyaac&
&p1(xyzzyaaaj47:xyzzyaaak47),xyzzyaacs1(xyzzyaaaj47:xyzzyaaak47),xyzzy&
&aahd1)
if(xyzzyaahf1)then
if(am_master)then
call wout(repeat('*',78))
call wout('DEBUG_INITIAL_PARAMETERS - determined parameters:')
call wout(repeat('*',78))
call wout('Term '//trim(i2s(iterm))//', channel '//trim(i2s(xyzzyaaaf4&
&7))//':')
do xyzzyaaao47=1,xyzzyaaaw1(iterm)
if(xyzzyaaci1(xyzzyaaai47+xyzzyaaao47)==pflag_det)then
call wout(trim(xyzzyaach1(xyzzyaaag47-1+xyzzyaaao47))//' = ',param(xyz&
&zyaaai47+xyzzyaaao47))
endif
enddo
call wout(repeat('*',78))
endif
endif
enddo
if(am_master)call wout()
call pop_casl_context()
enddo
deallocate(xyzzyaaef47,xyzzyaaeg47,xyzzyaacm47,xyzzyaacn47,xyzzyaaco47&
&,xyzzyaacp47,xyzzyaacq47,xyzzyaacr47,xyzzyaacs47,xyzzyaact47,xyzzyaac&
&u47,xyzzyaacv47,xyzzyaacw47,xyzzyaacx47,xyzzyaacy47,xyzzyaacz47,xyzzy&
&aada47)
if(am_master)then
do ispin=1,nspin
if(nele(ispin)==0)cycle
do jspin=ispin,nspin
if(nele(jspin)==0.or.(jspin==ispin.and.nele(ispin)<2))cycle
if(.not.xyzzyaaeb47(ispin,jspin).and..not.ee_cusp_in_orbital(ispin,jsp&
&in))then
if(abs(ee_kato_gamma(ispin,jspin,f_zero=f_zero))>1.d-11)then
if(trim(f_zero)/='quasi')then
call errwarn('READ_GJASTROW','Particle pair '//trim(i2s(ispin))//'-'//&
&trim(i2s(jspin))//' does not have cusp conditions applied anywhere in&
& the Jastrow factor.  This particle pair requires a function that goe&
&s as "a0 + a1 '//trim(f_zero)//'" near r=0 to keep the local energy f&
&inite at coalescence points.')
endif
endif
endif
enddo
enddo
call wout()
call wout('Finished General Jastrow setup.')
call wout()
endif
deallocate(xyzzyaaeb47,xyzzyaaec47)
if(xyzzyaaaa1==0)call errstop_master('READ_GJASTROW','Jastrow is empty!&
&')
call pop_casl_context()
call timer('READ_JASTROW',.false.)
end subroutine read_gjastrow
subroutine xyzzyaaht1(nparam,p_removed,eqn_matrix,eqn_rhs)
implicit none
integer,intent(in) :: nparam
real(dp),intent(inout) :: eqn_matrix(nparam,*),eqn_rhs(*)
logical,intent(in) :: p_removed(nparam)
integer xyzzyaaaa48,xyzzyaaab48
xyzzyaaaa48=0
do xyzzyaaab48=1,nparam
if(p_removed(xyzzyaaab48))then
xyzzyaaaa48=xyzzyaaaa48+1
eqn_matrix(xyzzyaaab48,xyzzyaaaa48)=1.d0
eqn_rhs(xyzzyaaaa48)=0.d0
endif
enddo
end subroutine xyzzyaaht1
subroutine xyzzyaahu1(rank_e,rank_n,size_sig,nparam,index_list,sig,sym&
&m_ee,p_removed,param_equal,equality_sign,neqn)
implicit none
integer,intent(in) :: rank_e,rank_n,size_sig,nparam,index_list(size_si&
&g,nparam),sig(size_sig),symm_ee
integer,intent(inout) :: neqn,param_equal(nparam),equality_sign(nparam&
&)
logical,intent(inout) :: p_removed(nparam)
integer,parameter :: xyzzyaaaa49=1024
integer xyzzyaaab49,xyzzyaaac49,xyzzyaaad49(size_sig),xyzzyaaae49(size&
&_sig),xyzzyaaaf49,xyzzyaaag49,xyzzyaaah49,xyzzyaaai49(size_sig,xyzzya&
&aaa49),xyzzyaaaj49,xyzzyaaak49,xyzzyaaal49,xyzzyaaam49(size_sig,xyzzy&
&aaaa49)
logical xyzzyaaan49(nparam),xyzzyaaao49,xyzzyaaap49,xyzzyaaaq49
neqn=0
param_equal(:)=0
equality_sign(:)=0
xyzzyaaan49(:)=.true.
do xyzzyaaab49=1,nparam
if(p_removed(xyzzyaaab49))cycle
if(param_equal(xyzzyaaab49)/=0)cycle
xyzzyaaad49(1:size_sig)=index_list(1:size_sig,xyzzyaaab49)
call sort_indices_matrix(rank_e,rank_n,size_sig,symm_ee,sig,xyzzyaaad4&
&9,xyzzyaaaf49,xyzzyaaai49,xyzzyaaah49,xyzzyaaao49)
if(xyzzyaaao49.and.xyzzyaaah49==0)cycle
xyzzyaaan49(xyzzyaaab49)=xyzzyaaan49(xyzzyaaab49).and.xyzzyaaao49
do xyzzyaaac49=xyzzyaaab49+1,nparam
if(p_removed(xyzzyaaac49))cycle
if(param_equal(xyzzyaaac49)/=0)cycle
xyzzyaaae49(1:size_sig)=index_list(1:size_sig,xyzzyaaac49)
call sort_indices_matrix(rank_e,rank_n,size_sig,symm_ee,sig,xyzzyaaae4&
&9,xyzzyaaaj49,xyzzyaaam49,xyzzyaaal49,xyzzyaaap49)
if(xyzzyaaap49.and.xyzzyaaal49==0)cycle
xyzzyaaan49(xyzzyaaac49)=xyzzyaaan49(xyzzyaaac49).and.xyzzyaaap49
if(all(xyzzyaaad49==xyzzyaaae49))then
xyzzyaaaq49=.false.
do xyzzyaaag49=1,xyzzyaaaf49
do xyzzyaaak49=1,xyzzyaaaj49
if(all(xyzzyaaai49(:,xyzzyaaag49)==xyzzyaaam49(:,xyzzyaaak49)))then
xyzzyaaaq49=.true.
exit
endif
enddo
if(xyzzyaaaq49)exit
enddo
if(xyzzyaaaq49)then
param_equal(xyzzyaaab49)=xyzzyaaab49
param_equal(xyzzyaaac49)=xyzzyaaab49
equality_sign(xyzzyaaac49)=xyzzyaaah49*xyzzyaaal49
xyzzyaaan49(xyzzyaaab49)=xyzzyaaan49(xyzzyaaab49).and.xyzzyaaao49.and.&
&xyzzyaaap49
neqn=neqn+1
endif
endif
enddo
enddo
do xyzzyaaab49=1,nparam
if(p_removed(xyzzyaaab49))cycle
if(xyzzyaaan49(xyzzyaaab49))cycle
p_removed(xyzzyaaab49)=.true.
do xyzzyaaac49=xyzzyaaab49+1,nparam
if(p_removed(xyzzyaaac49))cycle
p_removed(xyzzyaaac49)=.false.
enddo
enddo
end subroutine xyzzyaahu1
subroutine xyzzyaahv1(nparam,param_equal,equality_sign,eqn_matrix,eqn_&
&rhs)
implicit none
integer,intent(in) :: nparam,param_equal(nparam),equality_sign(nparam)
real(dp),intent(inout) :: eqn_matrix(nparam,*),eqn_rhs(*)
integer xyzzyaaaa50,xyzzyaaab50,xyzzyaaac50
xyzzyaaaa50=0
do xyzzyaaac50=1,nparam
xyzzyaaab50=param_equal(xyzzyaaac50)
if(xyzzyaaab50==0.or.xyzzyaaab50==xyzzyaaac50)cycle
xyzzyaaaa50=xyzzyaaaa50+1
eqn_matrix(xyzzyaaab50,xyzzyaaaa50)=dble(equality_sign(xyzzyaaac50))
eqn_matrix(xyzzyaaac50,xyzzyaaaa50)=-1.d0
eqn_rhs(xyzzyaaaa50)=0.d0
enddo
end subroutine xyzzyaahv1
subroutine xyzzyaahw1(size_ee,size_en,size_sig,nparam,index_list,sig,p&
&_removed,constraint_active,neqn)
implicit none
integer,intent(in) :: size_ee,size_en,size_sig,nparam,index_list(size_&
&sig,nparam),sig(size_sig)
integer,intent(inout) :: neqn
logical,intent(in) :: p_removed(nparam),constraint_active(xyzzyaadu1)
call xyzzyaahy1(size_ee,size_en,size_sig,nparam,index_list,sig,p_remov&
&ed,constraint_active,neqn=neqn)
end subroutine xyzzyaahw1
subroutine xyzzyaahx1(size_ee,size_en,size_sig,nparam,index_list,sig,p&
&_removed,constraint_active,eqn_matrix,eqn_rhs)
implicit none
integer,intent(in) :: size_ee,size_en,size_sig,nparam,index_list(size_&
&sig,nparam),sig(size_sig)
real(dp),intent(inout) :: eqn_matrix(nparam,*),eqn_rhs(*)
logical,intent(in) :: p_removed(nparam),constraint_active(xyzzyaadu1)
call xyzzyaahy1(size_ee,size_en,size_sig,nparam,index_list,sig,p_remov&
&ed,constraint_active,eqn_matrix=eqn_matrix,eqn_rhs=eqn_rhs)
end subroutine xyzzyaahx1
subroutine xyzzyaahy1(size_ee,size_en,size_sig,nparam,index_list,sig,p&
&_removed,constraint_active,neqn,eqn_matrix,eqn_rhs)
implicit none
integer,intent(in) :: size_ee,size_en,size_sig,nparam,index_list(size_&
&sig,nparam),sig(size_sig)
integer,intent(inout),optional :: neqn
real(dp),intent(inout),optional :: eqn_matrix(nparam,*),eqn_rhs(*)
logical,intent(in) :: p_removed(nparam),constraint_active(xyzzyaadu1)
integer xyzzyaaaa53,xyzzyaaab53,xyzzyaaac53,xyzzyaaad53(2),xyzzyaaae53&
&(2),xyzzyaaaf53(2),xyzzyaaag53(2),xyzzyaaah53
logical xyzzyaaai53
xyzzyaaai53=present(neqn)
xyzzyaaah53=0
do xyzzyaaaa53=1,xyzzyaadu1
if(.not.constraint_active(xyzzyaaaa53))cycle
select case(xyzzyaaaa53)
case(xyzzyaadv1)
do xyzzyaaab53=1,nparam
xyzzyaaad53=pack(mod(index_list(1:size_ee,xyzzyaaab53),dimensionality)&
&,index_list(1:size_ee,xyzzyaaab53)>0)
if(.not.all(xyzzyaaad53==1))cycle
xyzzyaaaf53=pack((index_list(1:size_ee,xyzzyaaab53)-1)/dimensionality,&
&index_list(1:size_ee,xyzzyaaab53)>0)
do xyzzyaaac53=xyzzyaaab53+1,nparam
if(any(index_list(1:size_ee,xyzzyaaab53)>0.neqv.index_list(1:size_ee,x&
&yzzyaaac53)>0))cycle
xyzzyaaag53=pack((index_list(1:size_ee,xyzzyaaac53)-1)/dimensionality,&
&index_list(1:size_ee,xyzzyaaac53)>0)
if(.not.all(xyzzyaaag53==xyzzyaaaf53))cycle
xyzzyaaae53=pack(mod(index_list(1:size_ee,xyzzyaaac53),dimensionality)&
&,index_list(1:size_ee,xyzzyaaac53)>0)
if(xyzzyaaae53(1)/=xyzzyaaae53(2))cycle
if(all(xyzzyaaae53==xyzzyaaad53))cycle
if(size_en>0)then
if(.not.(all(index_list(size_ee+1:size_ee+size_en,xyzzyaaac53)==index_&
&list(size_ee+1:size_ee+size_en,xyzzyaaab53))))cycle
endif
if(p_removed(xyzzyaaab53).and.p_removed(xyzzyaaac53))cycle
xyzzyaaah53=xyzzyaaah53+1
if(.not.xyzzyaaai53)then
if(.not.p_removed(xyzzyaaab53))eqn_matrix(xyzzyaaab53,xyzzyaaah53)=1.d&
&0
if(.not.p_removed(xyzzyaaac53))eqn_matrix(xyzzyaaac53,xyzzyaaah53)=-1.&
&d0
eqn_rhs(xyzzyaaah53)=0.d0
endif
enddo
enddo
case(xyzzyaadw1)
do xyzzyaaab53=1,nparam
xyzzyaaad53=pack(mod(index_list(size_ee+1:size_ee+size_en,xyzzyaaab53)&
&,dimensionality),index_list(size_ee+1:size_ee+size_en,xyzzyaaab53)>0)
if(.not.all(xyzzyaaad53==1))cycle
xyzzyaaaf53=pack((index_list(size_ee+1:size_ee+size_en,xyzzyaaab53)-1)&
&/dimensionality,index_list(size_ee+1:size_ee+size_en,xyzzyaaab53)>0)
do xyzzyaaac53=xyzzyaaab53+1,nparam
if(any(index_list(size_ee+1:size_ee+size_en,xyzzyaaab53)>0.neqv.index_&
&list(size_ee+1:size_ee+size_en,xyzzyaaac53)>0))cycle
xyzzyaaag53=pack((index_list(size_ee+1:size_ee+size_en,xyzzyaaac53)-1)&
&/dimensionality,index_list(size_ee+1:size_ee+size_en,xyzzyaaac53)>0)
if(.not.all(xyzzyaaag53==xyzzyaaaf53))cycle
xyzzyaaae53=pack(mod(index_list(size_ee+1:size_ee+size_en,xyzzyaaac53)&
&,dimensionality),index_list(size_ee+1:size_ee+size_en,xyzzyaaac53)>0)
if(xyzzyaaae53(1)/=xyzzyaaae53(2))cycle
if(all(xyzzyaaae53==xyzzyaaad53))cycle
if(size_ee>0)then
if(.not.(all(index_list(1:size_ee,xyzzyaaac53)==index_list(1:size_ee,x&
&yzzyaaab53))))cycle
endif
if(p_removed(xyzzyaaab53).and.p_removed(xyzzyaaac53))cycle
xyzzyaaah53=xyzzyaaah53+1
if(.not.xyzzyaaai53)then
if(.not.p_removed(xyzzyaaab53))eqn_matrix(xyzzyaaab53,xyzzyaaah53)=1.d&
&0
if(.not.p_removed(xyzzyaaac53))eqn_matrix(xyzzyaaac53,xyzzyaaah53)=-1.&
&d0
eqn_rhs(xyzzyaaah53)=0.d0
endif
enddo
enddo
end select
enddo
if(xyzzyaaai53)neqn=xyzzyaaah53
end subroutine xyzzyaahy1
subroutine xyzzyaahz1(rank_e,rank_n,size_ee,size_en,size_sig,order_ee,&
&order_en,nparam,index_list,sig,iset_eebasis,iset_eecut,       iset_en&
&basis,iset_encut,which_unity_eebasis,which_unity_eecut,         which&
&_unity_enbasis,which_unity_encut,fold_anisotropy_ee,              fol&
&d_anisotropy_en,which_ee_pair,which_en_pair,p_removed,cusp_value_ee,c&
&usp_value_en,nocoalesce_ee,nocoalesce_en,param_in_eqn,neqn)
implicit none
integer,intent(in) :: rank_e,rank_n,size_ee,size_en,size_sig,order_ee,&
&order_en,nparam,index_list(size_sig,nparam),sig(size_sig),iset_eebasi&
&s,iset_eecut,iset_enbasis,iset_encut,which_unity_eebasis,which_unity_&
&eecut,which_unity_enbasis,which_unity_encut,which_ee_pair(2,*),which_&
&en_pair(2,*)
integer,intent(inout) :: neqn,param_in_eqn(nparam,2,size_sig)
real(dp),intent(in) :: cusp_value_ee(*),cusp_value_en(*)
logical,intent(in) :: fold_anisotropy_ee,fold_anisotropy_en,p_removed(&
&nparam),nocoalesce_ee(*),nocoalesce_en(*)
integer xyzzyaaaa54,xyzzyaaab54,xyzzyaaac54,xyzzyaaad54,xyzzyaaae54,xy&
&zzyaaaf54,xyzzyaaag54,xyzzyaaah54,xyzzyaaai54,xyzzyaaaj54,xyzzyaaak54&
&,xyzzyaaal54,xyzzyaaam54,xyzzyaaan54,xyzzyaaao54,xyzzyaaap54(rank_e,r&
&ank_e),xyzzyaaaq54(rank_e,rank_n),xyzzyaaar54(rank_e,rank_e),xyzzyaaa&
&s54(rank_e,rank_n),xyzzyaaat54(rank_e),xyzzyaaau54(rank_n),xyzzyaaav5&
&4,xyzzyaaaw54
integer,allocatable :: xyzzyaaax54(:,:)
logical xyzzyaaay54(2)
logical,allocatable :: xyzzyaaaz54(:,:),xyzzyaaba54(:,:)
integer xyzzyaabb54(0:order_ee,0:order_ee),xyzzyaabc54(0:order_ee,0:or&
&der_ee),xyzzyaabd54(0:order_en,0:order_en),xyzzyaabe54(0:order_en,0:o&
&rder_en),xyzzyaabf54(0:order_ee,0:order_en)
integer xyzzyaabg54(order_ee),xyzzyaabh54(1),xyzzyaabi54(order_ee),xyz&
&zyaabj54(order_en),xyzzyaabk54(1),xyzzyaabl54(order_en)
real(dp) xyzzyaabm54(order_en),xyzzyaabn54(order_en),xyzzyaabo54(1),xy&
&zzyaabp54(1),xyzzyaabq54(order_ee),xyzzyaabr54(order_ee),xyzzyaabs54(&
&1),xyzzyaabt54(1)
real(dp),allocatable :: xyzzyaabu54(:,:,:),xyzzyaabv54(:,:,:),xyzzyaab&
&w54(:,:,:),xyzzyaabx54(:,:,:)
param_in_eqn(:,:,:)=0
neqn=0
allocate(xyzzyaabu54(order_ee,2,size_ee),xyzzyaabv54(order_ee,2,size_e&
&e),xyzzyaabw54(order_en,2,size_en),xyzzyaabx54(order_en,2,size_en),xy&
&zzyaaax54(size_sig,nparam),xyzzyaaaz54(2,nparam),xyzzyaaba54(2,nparam&
&),stat=xyzzyaaaa54)
call check_alloc(xyzzyaaaa54,'INIT_VALUE_CONSTRAINTS','target_at_zero_&
&*, ...')
xyzzyaabu54=0.d0
xyzzyaabv54=0.d0
xyzzyaabw54=0.d0
xyzzyaabx54=0.d0
xyzzyaaax54=0
xyzzyaaaz54=.false.
xyzzyaaba54=.false.
do xyzzyaaae54=1,size_ee
xyzzyaaan54=sig(xyzzyaaae54)
call query_cusp_eebasis(iset_eebasis,xyzzyaaan54,order_ee,.true.,xyzzy&
&aabq54(1),xyzzyaabr54(1),xyzzyaabg54(1))
call query_cusp_eebasis(iset_eecut,xyzzyaaan54,1,.true.,xyzzyaabs54(1)&
&,xyzzyaabt54(1),xyzzyaabh54(1))
call xyzzyaaib1(fold_anisotropy_ee,xyzzyaabh54(1),xyzzyaabg54,xyzzyaab&
&i54)
if(nocoalesce_ee(xyzzyaaan54))cycle
xyzzyaabu54(1:order_ee,1,xyzzyaaae54)=xyzzyaabs54(1)*xyzzyaabq54(1:ord&
&er_ee)
xyzzyaabu54(1:order_ee,2,xyzzyaaae54)=xyzzyaabt54(1)*xyzzyaabq54(1:ord&
&er_ee)+xyzzyaabs54(1)*xyzzyaabr54(1:order_ee)
where(xyzzyaabi54(1:order_ee)==0)xyzzyaabv54(1:order_ee,2,xyzzyaaae54)&
&=cusp_value_ee(xyzzyaaan54)
enddo
do xyzzyaaae54=1,size_en
xyzzyaaan54=sig(xyzzyaaae54+size_ee)
call query_cusp_enbasis(iset_enbasis,xyzzyaaan54,order_en,.true.,xyzzy&
&aabm54(1),xyzzyaabn54(1),xyzzyaabj54(1))
call query_cusp_enbasis(iset_encut,xyzzyaaan54,1,.true.,xyzzyaabo54(1)&
&,xyzzyaabp54(1),xyzzyaabk54(1))
call xyzzyaaib1(fold_anisotropy_en,xyzzyaabk54(1),xyzzyaabj54,xyzzyaab&
&l54)
if(nocoalesce_en(xyzzyaaan54))cycle
xyzzyaabw54(1:order_en,1,xyzzyaaae54)=xyzzyaabo54(1)*xyzzyaabm54(1:ord&
&er_en)
xyzzyaabw54(1:order_en,2,xyzzyaaae54)=xyzzyaabp54(1)*xyzzyaabm54(1:ord&
&er_en)+xyzzyaabo54(1)*xyzzyaabn54(1:order_en)
where(xyzzyaabl54(1:order_en)==0)xyzzyaabx54(1:order_en,2,xyzzyaaae54)&
&=cusp_value_en(xyzzyaaan54)
enddo
if(size_ee>1)then
call query_eqprod_ee_ee(iset_eebasis,iset_eecut,order_ee,.true.,xyzzya&
&abb54)
call query_eqprod_ee_ee(iset_eebasis,iset_eecut,order_ee,.false.,xyzzy&
&aabc54)
endif
if(size_en>1)then
call query_eqprod_en_en(iset_enbasis,iset_encut,order_en,.true.,xyzzya&
&abd54)
call query_eqprod_en_en(iset_enbasis,iset_encut,order_en,.false.,xyzzy&
&aabe54)
endif
if(size_en>0.and.size_ee>0)call query_eqprod_ee_en(iset_eebasis,iset_e&
&nbasis,iset_eecut,iset_encut,order_ee,order_en,xyzzyaabf54)
do xyzzyaaaj54=1,size_ee
if(all(xyzzyaabu54(1:order_ee,1:2,xyzzyaaaj54)==0.d0).and.all(xyzzyaab&
&v54(1:order_ee,1:2,xyzzyaaaj54)==0.d0))cycle
xyzzyaaba54=.false.
xyzzyaaaz54=.false.
xyzzyaaax54=0
xyzzyaaab54=which_ee_pair(1,xyzzyaaaj54)
xyzzyaaac54=which_ee_pair(2,xyzzyaaaj54)
call expindx2matrices(rank_e,rank_n,sig,xyzzyaaap54,xyzzyaaaq54)
do xyzzyaaae54=1,rank_e
xyzzyaaat54(xyzzyaaae54)=xyzzyaaae54
enddo
do xyzzyaaae54=1,rank_n
xyzzyaaau54(xyzzyaaae54)=xyzzyaaae54
enddo
call swap1(xyzzyaaat54(1),xyzzyaaat54(xyzzyaaab54))
call swap1(xyzzyaaat54(2),xyzzyaaat54(xyzzyaaac54))
do xyzzyaaag54=1,nparam
if(p_removed(xyzzyaaag54))then
xyzzyaaax54(:,xyzzyaaag54)=0
xyzzyaaba54(1:2,xyzzyaaag54)=(/.false.,.false./)
xyzzyaaaz54(1:2,xyzzyaaag54)=(/.false.,.false./)
cycle
endif
call expindx2matrices(rank_e,rank_n,index_list(1,xyzzyaaag54),xyzzyaaa&
&r54,xyzzyaaas54)
do xyzzyaaae54=3,rank_e
if(xyzzyaaap54(xyzzyaaat54(1),xyzzyaaat54(xyzzyaaae54))==xyzzyaaap54(x&
&yzzyaaat54(2),xyzzyaaat54(xyzzyaaae54)))then
xyzzyaaar54(xyzzyaaat54(2),xyzzyaaat54(xyzzyaaae54))=xyzzyaabb54(xyzzy&
&aaar54(xyzzyaaat54(1),xyzzyaaat54(xyzzyaaae54)),xyzzyaaar54(xyzzyaaat&
&54(2),xyzzyaaat54(xyzzyaaae54)))
else
xyzzyaaar54(xyzzyaaat54(2),xyzzyaaat54(xyzzyaaae54))=xyzzyaabc54(xyzzy&
&aaar54(xyzzyaaat54(1),xyzzyaaat54(xyzzyaaae54)),xyzzyaaar54(xyzzyaaat&
&54(2),xyzzyaaat54(xyzzyaaae54)))
endif
xyzzyaaar54(xyzzyaaat54(xyzzyaaae54),xyzzyaaat54(2))=xyzzyaaar54(xyzzy&
&aaat54(2),xyzzyaaat54(xyzzyaaae54))
enddo
do xyzzyaaae54=1,rank_n
if(xyzzyaaaq54(xyzzyaaat54(1),xyzzyaaau54(xyzzyaaae54))==xyzzyaaaq54(x&
&yzzyaaat54(2),xyzzyaaau54(xyzzyaaae54)))then
xyzzyaaas54(xyzzyaaat54(2),xyzzyaaau54(xyzzyaaae54))=xyzzyaabd54(xyzzy&
&aaas54(xyzzyaaat54(1),xyzzyaaau54(xyzzyaaae54)),xyzzyaaas54(xyzzyaaat&
&54(2),xyzzyaaau54(xyzzyaaae54)))
else
xyzzyaaas54(xyzzyaaat54(2),xyzzyaaau54(xyzzyaaae54))=xyzzyaabe54(xyzzy&
&aaas54(xyzzyaaat54(1),xyzzyaaau54(xyzzyaaae54)),xyzzyaaas54(xyzzyaaat&
&54(2),xyzzyaaau54(xyzzyaaae54)))
endif
enddo
xyzzyaaax54(1:size_sig,xyzzyaaag54)=0
xyzzyaaao54=1
xyzzyaaal54=xyzzyaaar54(xyzzyaaat54(1),xyzzyaaat54(2))
if(xyzzyaaal54>0)then
xyzzyaaax54(xyzzyaaao54,xyzzyaaag54)=xyzzyaabi54(xyzzyaaal54)
else
xyzzyaaax54(xyzzyaaao54,xyzzyaaag54)=0
endif
do xyzzyaaae54=2,rank_e
do xyzzyaaaf54=xyzzyaaae54+1,rank_e
xyzzyaaao54=xyzzyaaao54+1
xyzzyaaax54(xyzzyaaao54,xyzzyaaag54)=xyzzyaaar54(xyzzyaaat54(xyzzyaaaf&
&54),xyzzyaaat54(xyzzyaaae54))
enddo
enddo
do xyzzyaaae54=2,rank_e
do xyzzyaaaf54=1,rank_n
xyzzyaaao54=xyzzyaaao54+1
xyzzyaaax54(xyzzyaaao54,xyzzyaaag54)=xyzzyaaas54(xyzzyaaat54(xyzzyaaae&
&54),xyzzyaaau54(xyzzyaaaf54))
enddo
enddo
if(xyzzyaaal54==0)then
xyzzyaaba54(1:2,xyzzyaaag54)=(/.false.,.false./)
xyzzyaaaz54(1:2,xyzzyaaag54)=(/.false.,.false./)
else
if(xyzzyaabi54(xyzzyaaal54)==0)then
xyzzyaaba54(1:2,xyzzyaaag54)=(/.false.,.true./)
elseif(xyzzyaabi54(xyzzyaaal54)>0)then
xyzzyaaba54(1:2,xyzzyaaag54)=(/.true.,.false./)
else
xyzzyaaba54(1:2,xyzzyaaag54)=(/.true.,.true./)
endif
xyzzyaaaz54(1:2,xyzzyaaag54)=xyzzyaabu54(xyzzyaaal54,1:2,xyzzyaaaj54)/&
&=0.d0
xyzzyaaay54(1:2)=xyzzyaabv54(xyzzyaaal54,1:2,xyzzyaaaj54)/=0.d0
if(any(xyzzyaaay54))then
do xyzzyaaav54=1,size_ee
if(xyzzyaaav54==xyzzyaaaj54)cycle
xyzzyaaaw54=index_list(xyzzyaaav54,xyzzyaaag54)
if(xyzzyaaaw54==0.or.(which_unity_eebasis==xyzzyaaaw54.and.which_unity&
&_eecut==1))cycle
xyzzyaaay54(1:2)=.false.
exit
enddo
endif
if(any(xyzzyaaay54))then
do xyzzyaaav54=size_ee+1,size_ee+size_en
xyzzyaaaw54=index_list(xyzzyaaav54,xyzzyaaag54)
if(xyzzyaaaw54==0.or.(which_unity_enbasis==xyzzyaaaw54.and.which_unity&
&_encut==1))cycle
xyzzyaaay54(1:2)=.false.
exit
enddo
endif
xyzzyaaba54(1:2,xyzzyaaag54)=xyzzyaaba54(1:2,xyzzyaaag54).and.(xyzzyaa&
&az54(1:2,xyzzyaaag54).or.xyzzyaaay54(1:2))
endif
enddo
do xyzzyaaai54=1,2
do xyzzyaaag54=1,nparam
if(.not.xyzzyaaba54(xyzzyaaai54,xyzzyaaag54))cycle
neqn=neqn+1
xyzzyaaba54(xyzzyaaai54,xyzzyaaag54)=.false.
if(xyzzyaaaz54(xyzzyaaai54,xyzzyaaag54))param_in_eqn(xyzzyaaag54,xyzzy&
&aaai54,xyzzyaaaj54)=neqn
do xyzzyaaah54=xyzzyaaag54+1,nparam
if(.not.xyzzyaaba54(xyzzyaaai54,xyzzyaaah54))cycle
if(any(xyzzyaaax54(:,xyzzyaaag54)/=xyzzyaaax54(:,xyzzyaaah54)))cycle
xyzzyaaba54(xyzzyaaai54,xyzzyaaah54)=.false.
if(xyzzyaaaz54(xyzzyaaai54,xyzzyaaah54))param_in_eqn(xyzzyaaah54,xyzzy&
&aaai54,xyzzyaaaj54)=neqn
enddo
enddo
enddo
enddo
do xyzzyaaak54=1,size_en
if(all(xyzzyaabw54(1:order_en,1:2,xyzzyaaak54)==0.d0).and.all(xyzzyaab&
&x54(1:order_en,1:2,xyzzyaaak54)==0.d0))cycle
xyzzyaaba54=.false.
xyzzyaaaz54=.false.
xyzzyaaax54=0
xyzzyaaab54=which_en_pair(1,xyzzyaaak54)
xyzzyaaad54=which_en_pair(2,xyzzyaaak54)
call expindx2matrices(rank_e,rank_n,sig,xyzzyaaap54,xyzzyaaaq54)
do xyzzyaaae54=1,rank_e
xyzzyaaat54(xyzzyaaae54)=xyzzyaaae54
enddo
do xyzzyaaae54=1,rank_n
xyzzyaaau54(xyzzyaaae54)=xyzzyaaae54
enddo
call swap1(xyzzyaaat54(1),xyzzyaaat54(xyzzyaaab54))
call swap1(xyzzyaaau54(1),xyzzyaaau54(xyzzyaaad54))
do xyzzyaaag54=1,nparam
if(p_removed(xyzzyaaag54))then
xyzzyaaax54(:,xyzzyaaag54)=0
xyzzyaaba54(1:2,xyzzyaaag54)=(/.false.,.false./)
xyzzyaaaz54(1:2,xyzzyaaag54)=(/.false.,.false./)
cycle
endif
call expindx2matrices(rank_e,rank_n,index_list(1,xyzzyaaag54),xyzzyaaa&
&r54,xyzzyaaas54)
do xyzzyaaae54=2,rank_e
xyzzyaaar54(xyzzyaaat54(1),xyzzyaaat54(xyzzyaaae54))=xyzzyaabf54(xyzzy&
&aaar54(xyzzyaaat54(1),xyzzyaaat54(xyzzyaaae54)),xyzzyaaas54(xyzzyaaat&
&54(xyzzyaaae54),xyzzyaaau54(1)))
xyzzyaaar54(xyzzyaaat54(xyzzyaaae54),xyzzyaaat54(1))=xyzzyaaar54(xyzzy&
&aaat54(1),xyzzyaaat54(xyzzyaaae54))
enddo
xyzzyaaax54(1:size_sig,xyzzyaaag54)=0
xyzzyaaao54=1
xyzzyaaam54=xyzzyaaas54(xyzzyaaat54(1),xyzzyaaau54(1))
if(xyzzyaaam54>0)then
xyzzyaaax54(xyzzyaaao54,xyzzyaaag54)=xyzzyaabl54(xyzzyaaam54)
else
xyzzyaaax54(xyzzyaaao54,xyzzyaaag54)=0
endif
do xyzzyaaae54=1,rank_e
do xyzzyaaaf54=2,rank_n
xyzzyaaao54=xyzzyaaao54+1
xyzzyaaax54(xyzzyaaao54,xyzzyaaag54)=xyzzyaaas54(xyzzyaaat54(xyzzyaaae&
&54),xyzzyaaau54(xyzzyaaaf54))
enddo
enddo
do xyzzyaaae54=1,rank_e
do xyzzyaaaf54=xyzzyaaae54+1,rank_e
xyzzyaaao54=xyzzyaaao54+1
xyzzyaaax54(xyzzyaaao54,xyzzyaaag54)=xyzzyaaar54(xyzzyaaat54(xyzzyaaaf&
&54),xyzzyaaat54(xyzzyaaae54))
enddo
enddo
if(xyzzyaaam54==0)then
xyzzyaaba54(1:2,xyzzyaaag54)=(/.false.,.false./)
xyzzyaaaz54(1:2,xyzzyaaag54)=(/.false.,.false./)
else
if(xyzzyaabl54(xyzzyaaam54)==0)then
xyzzyaaba54(1:2,xyzzyaaag54)=(/.false.,.true./)
elseif(xyzzyaabl54(xyzzyaaam54)>0)then
xyzzyaaba54(1:2,xyzzyaaag54)=(/.true.,.false./)
else
xyzzyaaba54(1:2,xyzzyaaag54)=(/.true.,.true./)
endif
xyzzyaaaz54(1:2,xyzzyaaag54)=xyzzyaabw54(xyzzyaaam54,1:2,xyzzyaaak54)/&
&=0.d0
xyzzyaaay54(1:2)=xyzzyaabx54(xyzzyaaam54,1:2,xyzzyaaak54)/=0.d0
if(any(xyzzyaaay54))then
do xyzzyaaav54=1,size_ee
xyzzyaaaw54=index_list(xyzzyaaav54,xyzzyaaag54)
if(xyzzyaaaw54==0.or.(which_unity_eebasis==xyzzyaaaw54.and.which_unity&
&_eecut==1))cycle
xyzzyaaay54(1:2)=.false.
exit
enddo
endif
if(any(xyzzyaaay54))then
do xyzzyaaav54=size_ee+1,size_ee+size_en
if(xyzzyaaav54==size_ee+xyzzyaaak54)cycle
xyzzyaaaw54=index_list(xyzzyaaav54,xyzzyaaag54)
if(xyzzyaaaw54==0.or.(which_unity_enbasis==xyzzyaaaw54.and.which_unity&
&_encut==1))cycle
xyzzyaaay54(1:2)=.false.
exit
enddo
endif
xyzzyaaba54(1:2,xyzzyaaag54)=xyzzyaaba54(1:2,xyzzyaaag54).and.(xyzzyaa&
&az54(1:2,xyzzyaaag54).or.xyzzyaaay54(1:2))
endif
enddo
do xyzzyaaai54=1,2
do xyzzyaaag54=1,nparam
if(.not.xyzzyaaba54(xyzzyaaai54,xyzzyaaag54))cycle
neqn=neqn+1
xyzzyaaba54(xyzzyaaai54,xyzzyaaag54)=.false.
if(xyzzyaaaz54(xyzzyaaai54,xyzzyaaag54))param_in_eqn(xyzzyaaag54,xyzzy&
&aaai54,xyzzyaaak54+size_ee)=neqn
do xyzzyaaah54=xyzzyaaag54+1,nparam
if(.not.xyzzyaaba54(xyzzyaaai54,xyzzyaaah54))cycle
if(any(xyzzyaaax54(:,xyzzyaaag54)/=xyzzyaaax54(:,xyzzyaaah54)))cycle
xyzzyaaba54(xyzzyaaai54,xyzzyaaah54)=.false.
if(xyzzyaaaz54(xyzzyaaai54,xyzzyaaah54))param_in_eqn(xyzzyaaah54,xyzzy&
&aaai54,xyzzyaaak54+size_ee)=neqn
enddo
enddo
enddo
enddo
end subroutine xyzzyaahz1
subroutine xyzzyaaia1(size_ee,size_en,size_sig,order_ee,order_en,npara&
&m,index_list,sig,iset_eebasis,iset_eecut,iset_enbasis,iset_encut,whic&
&h_unity_eebasis,which_unity_eecut,which_unity_enbasis,which_unity_enc&
&ut,fold_anisotropy_ee,fold_anisotropy_en,param_in_eqn,p_removed,fake,&
&cusp_value_ee,cusp_value_en,nocoalesce_ee,nocoalesce_en,eqn_matrix,eq&
&n_rhs)
implicit none
integer,intent(in) :: size_ee,size_en,size_sig,order_ee,order_en,npara&
&m,index_list(size_sig,nparam),sig(size_sig),iset_eebasis,iset_eecut,i&
&set_enbasis,iset_encut,which_unity_eebasis,which_unity_eecut,which_un&
&ity_enbasis,which_unity_encut,param_in_eqn(nparam,2,size_sig)
real(dp),intent(in) :: cusp_value_ee(*),cusp_value_en(*)
real(dp),intent(inout) :: eqn_matrix(nparam,*),eqn_rhs(*)
logical,intent(in) :: fold_anisotropy_ee,fold_anisotropy_en,p_removed(&
&nparam),fake,nocoalesce_ee(*),nocoalesce_en(*)
integer xyzzyaaaa55,xyzzyaaab55,xyzzyaaac55,xyzzyaaad55,xyzzyaaae55,xy&
&zzyaaaf55,xyzzyaaag55,xyzzyaaah55,xyzzyaaai55
integer xyzzyaaaj55(order_ee),xyzzyaaak55(1),xyzzyaaal55(order_en),xyz&
&zyaaam55(1),xyzzyaaan55(order_ee),xyzzyaaao55(order_en),xyzzyaaap55(s&
&ize_sig)
real(dp) xyzzyaaaq55,xyzzyaaar55(order_ee),xyzzyaaas55(order_ee),xyzzy&
&aaat55(1),xyzzyaaau55(1),xyzzyaaav55(order_en),xyzzyaaaw55(order_en),&
&xyzzyaaax55(1),xyzzyaaay55(1),xyzzyaaaz55(max(order_ee,order_en),2,si&
&ze_sig),xyzzyaaba55(max(order_ee,order_en),2,size_sig)
xyzzyaaaz55=0.d0
xyzzyaaba55=0.d0
xyzzyaaap55=0
do xyzzyaaae55=1,size_ee
xyzzyaaag55=sig(xyzzyaaae55)
call query_cusp_eebasis(iset_eebasis,xyzzyaaag55,order_ee,fake,xyzzyaa&
&ar55,xyzzyaaas55,xyzzyaaaj55)
call query_cusp_eebasis(iset_eecut,xyzzyaaag55,1,fake,xyzzyaaat55,xyzz&
&yaaau55,xyzzyaaak55)
call xyzzyaaib1(fold_anisotropy_ee,xyzzyaaak55(1),xyzzyaaaj55,xyzzyaaa&
&n55)
if(nocoalesce_ee(xyzzyaaag55))cycle
xyzzyaaaz55(1:order_ee,1,xyzzyaaae55)=xyzzyaaat55(1)*xyzzyaaar55(1:ord&
&er_ee)
xyzzyaaaz55(1:order_ee,2,xyzzyaaae55)=xyzzyaaau55(1)*xyzzyaaar55(1:ord&
&er_ee)+xyzzyaaat55(1)*xyzzyaaas55(1:order_ee)
where(xyzzyaaan55(1:order_ee)==0)xyzzyaaba55(1:order_ee,2,xyzzyaaae55)&
&=cusp_value_ee(xyzzyaaag55)
if(which_unity_eecut==1)xyzzyaaap55(xyzzyaaae55)=which_unity_eebasis
enddo
do xyzzyaaaf55=1,size_en
xyzzyaaag55=sig(size_ee+xyzzyaaaf55)
call query_cusp_enbasis(iset_enbasis,xyzzyaaag55,order_en,fake,xyzzyaa&
&av55,xyzzyaaaw55,xyzzyaaal55)
call query_cusp_enbasis(iset_encut,xyzzyaaag55,1,fake,xyzzyaaax55,xyzz&
&yaaay55,xyzzyaaam55)
call xyzzyaaib1(fold_anisotropy_en,xyzzyaaam55(1),xyzzyaaal55,xyzzyaaa&
&o55)
if(nocoalesce_en(xyzzyaaag55))cycle
xyzzyaaaz55(1:order_en,1,size_ee+xyzzyaaaf55)=xyzzyaaax55(1)*xyzzyaaav&
&55(1:order_en)
xyzzyaaaz55(1:order_en,2,size_ee+xyzzyaaaf55)=xyzzyaaay55(1)*xyzzyaaav&
&55(1:order_en)+xyzzyaaax55(1)*xyzzyaaaw55(1:order_en)
where(xyzzyaaao55(1:order_en)==0)xyzzyaaba55(1:order_en,2,size_ee+xyzz&
&yaaaf55)=cusp_value_en(xyzzyaaag55)
if(which_unity_encut==1)xyzzyaaap55(size_ee+xyzzyaaaf55)=which_unity_e&
&nbasis
enddo
do xyzzyaaac55=1,size_sig
do xyzzyaaad55=1,2
do xyzzyaaaa55=1,nparam
if(p_removed(xyzzyaaaa55))cycle
xyzzyaaab55=param_in_eqn(xyzzyaaaa55,xyzzyaaad55,xyzzyaaac55)
if(xyzzyaaab55<1)cycle
eqn_matrix(xyzzyaaaa55,xyzzyaaab55)=xyzzyaaaz55(index_list(xyzzyaaac55&
&,xyzzyaaaa55),xyzzyaaad55,xyzzyaaac55)
xyzzyaaaq55=xyzzyaaba55(index_list(xyzzyaaac55,xyzzyaaaa55),xyzzyaaad5&
&5,xyzzyaaac55)
if(xyzzyaaaq55/=0.d0)then
do xyzzyaaai55=1,size_sig
if(xyzzyaaai55==xyzzyaaac55)cycle
xyzzyaaah55=index_list(xyzzyaaai55,xyzzyaaaa55)
if(xyzzyaaah55==0.or.xyzzyaaap55(xyzzyaaac55)==xyzzyaaah55)cycle
xyzzyaaaq55=0.d0
exit
enddo
endif
eqn_rhs(xyzzyaaab55)=xyzzyaaaq55
enddo
enddo
enddo
end subroutine xyzzyaaia1
subroutine xyzzyaaib1(fold_anisotropy,aniso_index_cut,aniso_index_basi&
&s,aniso_index)
implicit none
integer,intent(in) :: aniso_index_cut,aniso_index_basis(:)
integer,intent(inout) :: aniso_index(:)
logical,intent(in) :: fold_anisotropy
integer xyzzyaaaa56
if(aniso_index_cut==0)then
aniso_index=aniso_index_basis
elseif(aniso_index_cut/=0)then
xyzzyaaaa56=min(0,minval(aniso_index_basis))
where(aniso_index_basis<0)aniso_index=aniso_index_basis
where(aniso_index_basis>0)aniso_index=xyzzyaaaa56-aniso_index_basis
if(aniso_index_cut>0)then
where(aniso_index_basis==0)aniso_index=aniso_index_cut
else
xyzzyaaaa56=min(0,minval(aniso_index_basis))
where(aniso_index_basis==0)aniso_index=xyzzyaaaa56+aniso_index_cut
endif
endif
if(fold_anisotropy)then
where(aniso_index>0)aniso_index=1
where(aniso_index<0)aniso_index=-1
endif
end subroutine xyzzyaaib1
subroutine xyzzyaaic1
implicit none
character(512) errmsg
call set_casl_block(':parameters.casl:JASTROW',errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call push_casl_context(':parameters.casl:JASTROW')
call set_casl_block('TERM 1:Rank',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 1:Rank:%u',2,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 1:Rank:%u',0,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_block('TERM 1:e-e basis',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 1:e-e basis:Type','natural power',errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 1:e-e basis:Order',9,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_block('TERM 1:e-e cutoff',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 1:e-e cutoff:Type','alt polynomial',errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
if(nitot>0)then
call set_casl_block('TERM 2:Rank',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 2:Rank:%u',1,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 2:Rank:%u',1,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_block('TERM 2:e-n basis',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 2:e-n basis:Type','natural power',errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 2:e-n basis:Order',9,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_block('TERM 2:e-n cutoff',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 2:e-n cutoff:Type','alt polynomial',errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_block('TERM 3:Rank',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 3:Rank:%u',2,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 3:Rank:%u',1,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_block('TERM 3:e-e basis',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 3:e-e basis:Type','natural power',errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 3:e-e basis:Order',4,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_block('TERM 3:e-n basis',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 3:e-n basis:Type','natural power',errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 3:e-n basis:Order',4,errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_block('TERM 3:e-n cutoff',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
call set_casl_item('TERM 3:e-n cutoff:Type','alt polynomial',errmsg)
if(len_trim(errmsg)>0)call errstop_master('GENERATE_INITIAL_JASTROW',t&
&rim(errmsg))
endif
call pop_casl_context()
end subroutine xyzzyaaic1
subroutine update_gjastrow_casl
implicit none
integer xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58,xyzzyaaae58,xy&
&zzyaaaf58,xyzzyaaag58,xyzzyaaah58
if(xyzzyaaaa1==0)return
call push_casl_context(':parameters.casl:JASTROW')
do xyzzyaaaa58=1,xyzzyaaaa1
call push_casl_context('TERM '//trim(i2s(xyzzyaaaa58)))
call update_eebasis_casl('e-e basis',xyzzyaaaq1(xyzzyaaaa58))
call update_eebasis_casl('e-e cutoff',xyzzyaaar1(xyzzyaaaa58))
call update_enbasis_casl('e-n basis',xyzzyaaas1(xyzzyaaaa58))
call update_enbasis_casl('e-n cutoff',xyzzyaaat1(xyzzyaaaa58))
xyzzyaaad58=xyzzyaaaw1(xyzzyaaaa58)
xyzzyaaab58=xyzzyaaag1(xyzzyaaaa58)
xyzzyaaag58=xyzzyaaah1(xyzzyaaaa58)+1
xyzzyaaah58=xyzzyaaag58+xyzzyaaad58-1
do xyzzyaaac58=1,xyzzyaaav1(xyzzyaaaa58)
xyzzyaaae58=xyzzyaaaj1(xyzzyaaab58+xyzzyaaac58)+1
xyzzyaaaf58=xyzzyaaae58+xyzzyaaad58-1
call update_gparam_channel_casl('Linear parameters:Channel '//trim(xyz&
&zyaace1(xyzzyaaab58+xyzzyaaac58)),xyzzyaaad58,xyzzyaach1(xyzzyaaag58:&
&xyzzyaaah58),param(xyzzyaaae58:xyzzyaaaf58),xyzzyaaci1(xyzzyaaae58:xy&
&zzyaaaf58),xyzzyaahd1)
enddo
call pop_casl_context()
enddo
call pop_casl_context()
end subroutine update_gjastrow_casl
subroutine xyzzyaaid1(ii,iterm,eebasis1,eecut1,nzeecut1,enbasis1,encut&
&1,nzencut1,eebasis,eecut,nzeecut,enbasis,encut,nzencut,grad_eebasis1,&
&grad_eecut1,grad_enbasis1,grad_encut1,deebasis1,deecut1,denbasis1,den&
&cut1,gradr_eebasis1,gradr_enbasis1,lap_eebasis1,lap_eecut1,lap_enbasi&
&s1,lap_encut1,d2eebasis1,d2eecut1,d2enbasis1,d2encut1,lapr_eebasis1,l&
&apr_enbasis1,jas,gjas,ljas,djas,dgjas,dljas)
implicit none
integer,intent(in) :: ii,iterm
real(dp),intent(in) :: eebasis(nfn_eebasis,netot,*),eecut(nfn_eebasis,&
&netot,*),enbasis(nfn_enbasis,nitot,*),encut(nfn_enbasis,nitot,*),eeba&
&sis1(nfn_eebasis,*),eecut1(nfn_eebasis,*),enbasis1(nfn_enbasis,*),enc&
&ut1(nfn_enbasis,*)
real(dp),intent(in),optional :: grad_eebasis1(3,nfn_eebasis,*),grad_ee&
&cut1(3,nfn_eebasis,*),grad_enbasis1(3,nfn_enbasis,*),grad_encut1(3,nf&
&n_enbasis,*),deebasis1(nfn_eebasis,*),deecut1(nfn_eebasis,*),denbasis&
&1(nfn_enbasis,*),dencut1(nfn_enbasis,*),gradr_eebasis1(3,*),gradr_enb&
&asis1(3,*),lap_eebasis1(nfn_eebasis,*),lap_eecut1(nfn_eebasis,*),lap_&
&enbasis1(nfn_enbasis,*),lap_encut1(nfn_enbasis,*),d2eebasis1(nfn_eeba&
&sis,*),d2eecut1(nfn_eebasis,*),d2enbasis1(nfn_enbasis,*),d2encut1(nfn&
&_enbasis,*),lapr_eebasis1(*),lapr_enbasis1(*)
real(dp),intent(out),optional :: jas,gjas(3),ljas,djas(*),dgjas(3,*),d&
&ljas(*)
logical,intent(in) :: nzeecut(netot,netot),nzencut(nitot,netot),nzeecu&
&t1(netot),nzencut1(nitot)
integer xyzzyaaaa59,xyzzyaaab59
if(present(jas))jas=0.d0
if(present(gjas))gjas=0.d0
if(present(ljas))ljas=0.d0
xyzzyaaaa59=xyzzyaaag1(iterm)+1
xyzzyaaab59=xyzzyaaai1(iterm)
select case(xyzzyaaba1(iterm))
case(xyzzyaafw1)
call xyzzyaaih1(ii,xyzzyaaam1(iterm),xyzzyaaay1(1,which_spin(ii),iterm&
&),param(1),xyzzyaaaj1(xyzzyaaaa59),eebasis1,eecut1,nzeecut1,grad_eeba&
&sis1,grad_eecut1,lap_eebasis1,lap_eecut1,jas,gjas,ljas,djas,dgjas,dlj&
&as)
case(xyzzyaagc1)
call xyzzyaaii1(ii,xyzzyaaam1(iterm),xyzzyaaay1(1,which_spin(ii),iterm&
&),param(1),xyzzyaaaj1(xyzzyaaaa59),eebasis1,eecut1,nzeecut1,deebasis1&
&,deecut1,gradr_eebasis1,d2eebasis1,d2eecut1,lapr_eebasis1,jas,gjas,lj&
&as,djas,dgjas,dljas)
case(xyzzyaafx1)
call xyzzyaaik1(xyzzyaaan1(iterm),xyzzyaaaz1(1,which_spin(ii),iterm),p&
&aram(1),xyzzyaaaj1(xyzzyaaaa59),enbasis1,encut1,nzencut1,grad_enbasis&
&1,grad_encut1,lap_enbasis1,lap_encut1,jas,gjas,ljas,djas,dgjas,dljas)
case(xyzzyaagd1)
call xyzzyaail1(xyzzyaaan1(iterm),xyzzyaaaz1(1,which_spin(ii),iterm),p&
&aram(1),xyzzyaaaj1(xyzzyaaaa59),enbasis1,encut1,nzencut1,denbasis1,de&
&ncut1,gradr_enbasis1,d2enbasis1,d2encut1,lapr_enbasis1,jas,gjas,ljas,&
&djas,dgjas,dljas)
case(xyzzyaafy1)
if(.not.present(dljas))then
call xyzzyaaim1(ii,xyzzyaaam1(iterm),xyzzyaaan1(iterm),xyzzyaaay1(1,wh&
&ich_spin(ii),iterm),xyzzyaaaz1(1,1,iterm),param(xyzzyaaaj1(xyzzyaaaa5&
&9)+1),xyzzyaaav1(iterm),xyzzyaacf1(1,xyzzyaaaa59),eebasis1,eecut1,nze&
&ecut1,enbasis1,encut1,nzencut1,enbasis,encut,nzencut,grad_eebasis1,gr&
&ad_eecut1,grad_enbasis1,grad_encut1,lap_eebasis1,lap_eecut1,lap_enbas&
&is1,lap_encut1,jas,gjas,ljas)
else
call xyzzyaaim1(ii,xyzzyaaam1(iterm),xyzzyaaan1(iterm),xyzzyaaay1(1,wh&
&ich_spin(ii),iterm),xyzzyaaaz1(1,1,iterm),param(xyzzyaaaj1(xyzzyaaaa5&
&9)+1),xyzzyaaav1(iterm),xyzzyaacf1(1,xyzzyaaaa59),eebasis1,eecut1,nze&
&ecut1,enbasis1,encut1,nzencut1,enbasis,encut,nzencut,grad_eebasis1,gr&
&ad_eecut1,grad_enbasis1,grad_encut1,lap_eebasis1,lap_eecut1,lap_enbas&
&is1,lap_encut1,jas,gjas,ljas,djas(xyzzyaaaj1(xyzzyaaaa59)+1),dgjas(1,&
&xyzzyaaaj1(xyzzyaaaa59)+1),dljas(xyzzyaaaj1(xyzzyaaaa59)+1))
endif
case(xyzzyaage1)
if(.not.present(dljas))then
call xyzzyaain1(ii,xyzzyaaam1(iterm),xyzzyaaan1(iterm),xyzzyaaay1(1,wh&
&ich_spin(ii),iterm),xyzzyaaaz1(1,1,iterm),xyzzyaabf1(iterm)==1.and.xy&
&zzyaabg1(iterm)==1,param(xyzzyaaaj1(xyzzyaaaa59)+1),xyzzyaaav1(iterm)&
&,xyzzyaacf1(1,xyzzyaaaa59),eebasis1,eecut1,nzeecut1,enbasis1,encut1,n&
&zencut1,enbasis,encut,nzencut,deebasis1,deecut1,denbasis1,dencut1,gra&
&dr_eebasis1,gradr_enbasis1,d2eebasis1,d2eecut1,d2enbasis1,d2encut1,la&
&pr_eebasis1,lapr_enbasis1,jas,gjas,ljas)
else
call xyzzyaain1(ii,xyzzyaaam1(iterm),xyzzyaaan1(iterm),xyzzyaaay1(1,wh&
&ich_spin(ii),iterm),xyzzyaaaz1(1,1,iterm),xyzzyaabf1(iterm)==1.and.xy&
&zzyaabg1(iterm)==1,param(xyzzyaaaj1(xyzzyaaaa59)+1),xyzzyaaav1(iterm)&
&,xyzzyaacf1(1,xyzzyaaaa59),eebasis1,eecut1,nzeecut1,enbasis1,encut1,n&
&zencut1,enbasis,encut,nzencut,deebasis1,deecut1,denbasis1,dencut1,gra&
&dr_eebasis1,gradr_enbasis1,d2eebasis1,d2eecut1,d2enbasis1,d2encut1,la&
&pr_eebasis1,lapr_enbasis1,jas,gjas,ljas,djas(xyzzyaaaj1(xyzzyaaaa59)+&
&1),dgjas(1,xyzzyaaaj1(xyzzyaaaa59)+1),dljas(xyzzyaaaj1(xyzzyaaaa59)+1&
&))
endif
case(xyzzyaafz1)
call xyzzyaaig1(ii,xyzzyaaak1(iterm),xyzzyaaao1(iterm),xyzzyaacg1(xyzz&
&yaaab59),xyzzyaabz1(xyzzyaaaa59),xyzzyaacj1(1),xyzzyaack1(1),xyzzyaac&
&m1(1),xyzzyaabp1(iterm),xyzzyaaay1(1,1,iterm),param(1),xyzzyaaav1(ite&
&rm),xyzzyaacf1(1,xyzzyaaaa59),xyzzyaaaj1(xyzzyaaaa59),xyzzyaabx1(1,1,&
&iterm),xyzzyaabc1(iterm),eebasis1,eecut1,nzeecut1,eebasis,eecut,nzeec&
&ut,grad_eebasis1,grad_eecut1,lap_eebasis1,lap_eecut1,jas,gjas,ljas,dj&
&as,dgjas,dljas)
case(xyzzyaaga1)
call xyzzyaaij1(ii,xyzzyaaal1(iterm),xyzzyaaap1(iterm),xyzzyaacg1(xyzz&
&yaaab59),xyzzyaabz1(xyzzyaaaa59),xyzzyaacj1(1),xyzzyaacl1(1),xyzzyaac&
&n1(1),xyzzyaabq1(iterm),xyzzyaaaz1(1,1,iterm),param(1),xyzzyaaav1(ite&
&rm),xyzzyaacf1(1,xyzzyaaaa59),xyzzyaaaj1(xyzzyaaaa59),xyzzyaaby1(1,1,&
&iterm),enbasis1,encut1,nzencut1,grad_enbasis1,grad_encut1,lap_enbasis&
&1,lap_encut1,jas,gjas,ljas,djas,dgjas,dljas)
case(xyzzyaagb1)
call xyzzyaaif1(ii,xyzzyaaak1(iterm),xyzzyaaal1(iterm),xyzzyaaao1(iter&
&m),xyzzyaaap1(iterm),xyzzyaaau1(iterm),xyzzyaacg1(xyzzyaaab59),xyzzya&
&abz1(xyzzyaaaa59),xyzzyaacj1(1),xyzzyaack1(1),xyzzyaacm1(1),xyzzyaacl&
&1(1),xyzzyaacn1(1),xyzzyaabp1(iterm),xyzzyaabq1(iterm),xyzzyaaay1(1,1&
&,iterm),xyzzyaaaz1(1,1,iterm),param(1),xyzzyaaav1(iterm),xyzzyaacf1(1&
&,xyzzyaaaa59),xyzzyaaaj1(xyzzyaaaa59),xyzzyaabx1(1,1,iterm),xyzzyaaby&
&1(1,1,iterm),xyzzyaabc1(iterm),eebasis1,eecut1,nzeecut1,enbasis1,encu&
&t1,nzencut1,eebasis,eecut,nzeecut,enbasis,encut,nzencut,grad_eebasis1&
&,grad_eecut1,grad_enbasis1,grad_encut1,lap_eebasis1,lap_eecut1,lap_en&
&basis1,lap_encut1,jas,gjas,ljas,djas,dgjas,dljas)
end select
end subroutine xyzzyaaid1
subroutine xyzzyaaie1(iterm,eebasis,nzeecut,enbasis,nzencut,jas)
implicit none
integer,intent(in) :: iterm
real(dp),intent(in) :: eebasis(nfn_eebasis,netot,netot),enbasis(nfn_en&
&basis,nitot,netot)
real(dp),intent(out) :: jas
logical,intent(in) :: nzeecut(netot,netot,*),nzencut(nitot,netot,*)
integer xyzzyaaaa60,xyzzyaaab60,xyzzyaaac60,xyzzyaaad60,xyzzyaaae60,xy&
&zzyaaaf60,xyzzyaaag60,xyzzyaaah60,xyzzyaaai60,xyzzyaaaj60,xyzzyaaak60&
&,xyzzyaaal60,xyzzyaaam60,xyzzyaaan60,xyzzyaaao60
jas=0.d0
xyzzyaaac60=xyzzyaaag1(iterm)+1
xyzzyaaag60=xyzzyaaai1(iterm)
select case(xyzzyaabb1(iterm))
case(xyzzyaagf1)
xyzzyaaaa60=xyzzyaaao1(iterm)
xyzzyaaah60=xyzzyaaaq1(iterm)
xyzzyaaai60=ifn1_eebasis(xyzzyaaah60)
xyzzyaaaj60=xyzzyaaar1(iterm)
xyzzyaaak60=ifn1_eebasis(xyzzyaaaj60)
xyzzyaaad60=xyzzyaaak1(iterm)
call xyzzyaair1(xyzzyaaad60,xyzzyaaaa60,xyzzyaacg1(xyzzyaaag60),xyzzya&
&abz1(xyzzyaaac60),xyzzyaacj1(1),xyzzyaack1(1),xyzzyaacm1(1),xyzzyaabp&
&1(iterm),xyzzyaaay1(1,1,iterm),param(1),xyzzyaaav1(iterm),xyzzyaacf1(&
&1,xyzzyaaac60),xyzzyaaaj1(xyzzyaaac60),xyzzyaabx1(1,1,iterm),eebasis(&
&xyzzyaaai60,1,1),eebasis(xyzzyaaak60,1,1),nzeecut(1,1,xyzzyaaaj60),ja&
&s)
case(xyzzyaagg1)
xyzzyaaab60=xyzzyaaap1(iterm)
xyzzyaaal60=xyzzyaaas1(iterm)
xyzzyaaam60=ifn1_enbasis(xyzzyaaal60)
xyzzyaaan60=xyzzyaaat1(iterm)
xyzzyaaao60=ifn1_enbasis(xyzzyaaan60)
xyzzyaaae60=xyzzyaaal1(iterm)
call xyzzyaais1(xyzzyaaae60,xyzzyaaab60,xyzzyaacg1(xyzzyaaag60),xyzzya&
&abz1(xyzzyaaac60),xyzzyaacj1(1),xyzzyaacl1(1),xyzzyaacn1(1),xyzzyaabq&
&1(iterm),xyzzyaaaz1(1,1,iterm),param(1),xyzzyaaav1(iterm),xyzzyaacf1(&
&1,xyzzyaaac60),xyzzyaaaj1(xyzzyaaac60),xyzzyaaby1(1,1,iterm),enbasis(&
&xyzzyaaam60,1,1),enbasis(xyzzyaaao60,1,1),nzencut(1,1,xyzzyaaan60),ja&
&s)
case(xyzzyaagh1)
xyzzyaaaa60=xyzzyaaao1(iterm)
xyzzyaaab60=xyzzyaaap1(iterm)
xyzzyaaah60=xyzzyaaaq1(iterm)
xyzzyaaai60=ifn1_eebasis(xyzzyaaah60)
xyzzyaaaj60=xyzzyaaar1(iterm)
xyzzyaaak60=ifn1_eebasis(xyzzyaaaj60)
xyzzyaaal60=xyzzyaaas1(iterm)
xyzzyaaam60=ifn1_enbasis(xyzzyaaal60)
xyzzyaaan60=xyzzyaaat1(iterm)
xyzzyaaao60=ifn1_enbasis(xyzzyaaan60)
xyzzyaaad60=xyzzyaaak1(iterm)
xyzzyaaae60=xyzzyaaal1(iterm)
xyzzyaaaf60=xyzzyaaau1(iterm)
call xyzzyaaiq1(xyzzyaaad60,xyzzyaaae60,xyzzyaaaa60,xyzzyaaab60,xyzzya&
&aaf60,xyzzyaacg1(xyzzyaaag60),xyzzyaabz1(xyzzyaaac60),xyzzyaacj1(1),x&
&yzzyaack1(1),xyzzyaacm1(1),xyzzyaacl1(1),xyzzyaacn1(1),xyzzyaabp1(ite&
&rm),xyzzyaabq1(iterm),xyzzyaaay1(1,1,iterm),xyzzyaaaz1(1,1,iterm),par&
&am(1),xyzzyaaav1(iterm),xyzzyaacf1(1,xyzzyaaac60),xyzzyaaaj1(xyzzyaaa&
&c60),xyzzyaabx1(1,1,iterm),xyzzyaaby1(1,1,iterm),eebasis(xyzzyaaai60,&
&1,1),eebasis(xyzzyaaak60,1,1),nzeecut(1,1,xyzzyaaaj60),enbasis(xyzzya&
&aam60,1,1),enbasis(xyzzyaaao60,1,1),nzencut(1,1,xyzzyaaan60),jas)
end select
end subroutine xyzzyaaie1
subroutine xyzzyaaif1(ii,rank_e,rank_n,size_ee,size_en,size_sig,index_&
&list,nparam_remap,index_remap,leftmost_change_ee,mask_leftmost_change&
&_ee,leftmost_change_en,mask_leftmost_change_en,uncluster_ee,uncluster&
&_en,groups_ee,groups_en,param,nchannel,chsig,iparam0,which_ee_pair,wh&
&ich_en_pair,symm_ee,eebasis1,eecut1,nzeecut1,enbasis1,encut1,nzencut1&
&,eebasis,eecut,nzeecut,enbasis,encut,nzencut,grad_eebasis1,grad_eecut&
&1,grad_enbasis1,grad_encut1,lap_eebasis1,lap_eecut1,lap_enbasis1,lap_&
&encut1,jas,gjas,ljas,djas,dgjas,dljas)
implicit none
integer,intent(in) :: ii,rank_e,rank_n,size_ee,size_en,size_sig,nchann&
&el,groups_ee(nspin,nspin),groups_en(nitot,nspin)
integer,intent(in) :: index_list(size_sig,*),nparam_remap(nchannel),in&
&dex_remap(*),leftmost_change_ee(*),mask_leftmost_change_ee(*),leftmos&
&t_change_en(*),mask_leftmost_change_en(*)
integer,intent(in) :: chsig(xyzzyaaae1,*),iparam0(*),which_ee_pair(2,*&
&),which_en_pair(2,*),symm_ee
real(dp),intent(in) :: eebasis1(nfn_eebasis,*),eecut1(nfn_eebasis,*),e&
&nbasis1(nfn_enbasis,*),encut1(nfn_enbasis,*),eebasis(nfn_eebasis,neto&
&t,*),eecut(nfn_eebasis,netot,*),enbasis(nfn_enbasis,nitot,*),encut(nf&
&n_enbasis,nitot,*),param(*)
real(dp),intent(in),optional :: grad_eebasis1(3,nfn_eebasis,*),grad_ee&
&cut1(3,nfn_eebasis,*),grad_enbasis1(3,nfn_enbasis,*),grad_encut1(3,nf&
&n_enbasis,*),lap_eebasis1(nfn_eebasis,*),lap_eecut1(nfn_eebasis,*),  &
&    lap_enbasis1(nfn_enbasis,*),lap_encut1(nfn_enbasis,*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(xyzzyaaad1),d&
&gjas(3,xyzzyaaad1),dljas(xyzzyaaad1)
logical,intent(in) :: nzeecut1(netot),nzencut1(nitot),nzeecut(netot,ne&
&tot),nzencut(nitot,netot),uncluster_ee,uncluster_en
integer xyzzyaaaa61,xyzzyaaab61(rank_e),xyzzyaaac61(rank_e),xyzzyaaad6&
&1(rank_n),xyzzyaaae61(size_ee),xyzzyaaaf61(size_en),xyzzyaaag61(size_&
&sig),xyzzyaaah61,xyzzyaaai61,xyzzyaaaj61,xyzzyaaak61,xyzzyaaal61,xyzz&
&yaaam61,xyzzyaaan61(rank_e),xyzzyaaao61(rank_n),xyzzyaaap61(netot,ran&
&k_n),xyzzyaaaq61(nitot,rank_n),xyzzyaaar61(netot,rank_e),xyzzyaaas61(&
&rank_n),xyzzyaaat61(rank_n),xyzzyaaau61(rank_e),xyzzyaaav61(rank_n),x&
&yzzyaaaw61(rank_e)
real(dp) xyzzyaaax61,xyzzyaaay61,xyzzyaaaz61,xyzzyaaba61,xyzzyaabb61(3&
&),xyzzyaabc61(3),xyzzyaabd61(3),xyzzyaabe61(3),xyzzyaabf61,xyzzyaabg6&
&1,xyzzyaabh61,xyzzyaabi61,xyzzyaabj61(size_ee),xyzzyaabk61(size_ee),x&
&yzzyaabl61(size_en),xyzzyaabm61(size_en),xyzzyaabn61(3,size_ee),xyzzy&
&aabo61(3,size_ee),xyzzyaabp61(3,size_en),xyzzyaabq61(3,size_en),xyzzy&
&aabr61(size_ee),xyzzyaabs61(size_ee),xyzzyaabt61(size_en),xyzzyaabu61&
&(size_en),xyzzyaabv61,xyzzyaabw61(3),xyzzyaabx61,xyzzyaaby61
logical xyzzyaabz61,xyzzyaaca61,xyzzyaacb61,xyzzyaacc61,xyzzyaacd61,xy&
&zzyaace61
xyzzyaacc61=present(dljas)
xyzzyaabz61=.true.
if(xyzzyaacc61.or.present(ljas))then
xyzzyaaca61=.true.
xyzzyaacb61=.true.
elseif(present(gjas))then
xyzzyaaca61=.true.
xyzzyaacb61=.false.
else
xyzzyaaca61=.false.
xyzzyaacb61=.false.
endif
if(xyzzyaabz61)jas=0.d0
if(xyzzyaaca61)then
gjas=0.d0
if(xyzzyaacb61)ljas=0.d0
endif
if(xyzzyaacc61)then
if(xyzzyaabz61)djas=0.d0
if(xyzzyaaca61)then
dgjas=0.d0
if(xyzzyaacb61)dljas=0.d0
endif
endif
xyzzyaaao61=0
xyzzyaaad61(1)=0
do while(iterate_nuclei_indices_fix(rank_e,rank_n,netot,nitot,ii,nzenc&
&ut1,nzencut,uncluster_en,xyzzyaaap61,xyzzyaaaq61,xyzzyaaas61,xyzzyaaa&
&t61,xyzzyaaav61,xyzzyaaad61))
xyzzyaacd61=all(xyzzyaaad61==xyzzyaaao61)
xyzzyaaao61=xyzzyaaad61
xyzzyaaan61=0
xyzzyaaac61(1)=0
do while(iterate_electron_indices_fix(rank_e,netot,ii,nzeecut1,nzeecut&
&,uncluster_ee,xyzzyaaar61,xyzzyaaau61,xyzzyaaaw61,xyzzyaaac61,preinit&
&=xyzzyaaap61(1,rank_n),preinit_max=xyzzyaaas61(rank_n)))
do xyzzyaaai61=1,rank_e
xyzzyaaab61(xyzzyaaai61)=which_spin(xyzzyaaac61(xyzzyaaai61))
enddo
xyzzyaace61=all(xyzzyaaab61==xyzzyaaan61)
xyzzyaaan61=xyzzyaaab61
if(.not.(xyzzyaace61.and.xyzzyaacd61))then
call get_sig(rank_e,rank_n,groups_ee,groups_en,xyzzyaaab61,xyzzyaaad61&
&,xyzzyaaag61,xyzzyaaae61,xyzzyaaaf61)
call match_signature(size_sig,xyzzyaaae1,nchannel,xyzzyaaag61,chsig,xy&
&zzyaaaj61)
xyzzyaacd61=.true.
endif
if(xyzzyaaaj61<1)cycle
if(symm_ee==-1)call insertion_point(rank_e,xyzzyaaac61,xyzzyaaam61)
xyzzyaaax61=1.d0
xyzzyaaaz61=1.d0
xyzzyaaay61=1.d0
xyzzyaaba61=1.d0
xyzzyaabb61=0.d0
xyzzyaabd61=0.d0
xyzzyaabc61=0.d0
xyzzyaabe61=0.d0
xyzzyaabf61=0.d0
xyzzyaabh61=0.d0
xyzzyaabg61=0.d0
xyzzyaabi61=0.d0
xyzzyaaal61=iparam0(xyzzyaaaj61)
do xyzzyaaah61=1,nparam_remap(xyzzyaaaj61)
xyzzyaaaa61=index_remap(xyzzyaaal61+xyzzyaaah61)
xyzzyaaak61=xyzzyaaaa61-xyzzyaaal61
xyzzyaaby61=param(xyzzyaaaa61)
if(symm_ee==-1)call apply_insertion_sign(xyzzyaaam61,xyzzyaaae61,index&
&_list(1,xyzzyaaak61),xyzzyaaby61)
if(xyzzyaacb61)then
if(mask_leftmost_change_ee(xyzzyaaaa61)>0)call xyzzyaait1(nfn_eebasis,&
&netot,rank_e,rank_e,size_ee,mask_leftmost_change_ee(xyzzyaaaa61),xyzz&
&yaaae61,index_list(1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaac61,.true.,.&
&true.,which_ee_pair,eecut,eecut1,grad_eecut1,lap_eecut1,xyzzyaabk61,x&
&yzzyaabo61,xyzzyaabs61,xyzzyaaay61,xyzzyaabc61,xyzzyaabg61)
if(mask_leftmost_change_en(xyzzyaaaa61)>0)call xyzzyaait1(nfn_enbasis,&
&nitot,rank_e,rank_n,size_en,mask_leftmost_change_en(xyzzyaaaa61),xyzz&
&yaaaf61,index_list(size_ee+1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaad61,&
&.false.,.true.,which_en_pair,encut,encut1,grad_encut1,lap_encut1,xyzz&
&yaabm61,xyzzyaabq61,xyzzyaabu61,xyzzyaaba61,xyzzyaabe61,xyzzyaabi61)
if(leftmost_change_ee(xyzzyaaaa61)>0)call xyzzyaait1(nfn_eebasis,netot&
&,rank_e,rank_e,size_ee,leftmost_change_ee(xyzzyaaaa61),xyzzyaaae61,in&
&dex_list(1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaac61,.true.,.false.,whi&
&ch_ee_pair,eebasis,eebasis1,grad_eebasis1,lap_eebasis1,xyzzyaabj61,xy&
&zzyaabn61,xyzzyaabr61,xyzzyaaax61,xyzzyaabb61,xyzzyaabf61)
if(leftmost_change_en(xyzzyaaaa61)>0)call xyzzyaait1(nfn_enbasis,nitot&
&,rank_e,rank_n,size_en,leftmost_change_en(xyzzyaaaa61),xyzzyaaaf61,in&
&dex_list(size_ee+1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaad61,.false.,.f&
&alse.,which_en_pair,enbasis,enbasis1,grad_enbasis1,lap_enbasis1,xyzzy&
&aabl61,xyzzyaabp61,xyzzyaabt61,xyzzyaaaz61,xyzzyaabd61,xyzzyaabh61)
xyzzyaabv61=xyzzyaaay61*xyzzyaaba61*xyzzyaaax61*xyzzyaaaz61
xyzzyaabw61=xyzzyaabc61+xyzzyaabe61+xyzzyaabb61+xyzzyaabd61
xyzzyaabx61=xyzzyaabg61+xyzzyaabi61+xyzzyaabf61+xyzzyaabh61+ddot(3,xyz&
&zyaabw61,1,xyzzyaabw61,1)
if(xyzzyaacc61)then
djas(xyzzyaaaa61)=djas(xyzzyaaaa61)+xyzzyaabv61
dgjas(:,xyzzyaaaa61)=dgjas(:,xyzzyaaaa61)+xyzzyaabv61*xyzzyaabw61
dljas(xyzzyaaaa61)=dljas(xyzzyaaaa61)+xyzzyaabv61*xyzzyaabx61
endif
xyzzyaabv61=xyzzyaabv61*xyzzyaaby61
jas=jas+xyzzyaabv61
gjas=gjas+xyzzyaabv61*xyzzyaabw61
ljas=ljas+xyzzyaabv61*xyzzyaabx61
elseif(xyzzyaaca61)then
if(mask_leftmost_change_ee(xyzzyaaaa61)>0)call xyzzyaaiu1(nfn_eebasis,&
&netot,rank_e,rank_e,size_ee,mask_leftmost_change_ee(xyzzyaaaa61),xyzz&
&yaaae61,index_list(1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaac61,.true.,.&
&true.,which_ee_pair,eecut,eecut1,grad_eecut1,xyzzyaabk61,xyzzyaabo61,&
&xyzzyaaay61,xyzzyaabc61)
if(mask_leftmost_change_en(xyzzyaaaa61)>0)call xyzzyaaiu1(nfn_enbasis,&
&nitot,rank_e,rank_n,size_en,mask_leftmost_change_en(xyzzyaaaa61),xyzz&
&yaaaf61,index_list(size_ee+1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaad61,&
&.false.,.true.,which_en_pair,encut,encut1,grad_encut1,xyzzyaabm61,xyz&
&zyaabq61,xyzzyaaba61,xyzzyaabe61)
if(leftmost_change_ee(xyzzyaaaa61)>0)call xyzzyaaiu1(nfn_eebasis,netot&
&,rank_e,rank_e,size_ee,leftmost_change_ee(xyzzyaaaa61),xyzzyaaae61,in&
&dex_list(1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaac61,.true.,.false.,whi&
&ch_ee_pair,eebasis,eebasis1,grad_eebasis1,xyzzyaabj61,xyzzyaabn61,xyz&
&zyaaax61,xyzzyaabb61)
if(leftmost_change_en(xyzzyaaaa61)>0)call xyzzyaaiu1(nfn_enbasis,nitot&
&,rank_e,rank_n,size_en,leftmost_change_en(xyzzyaaaa61),xyzzyaaaf61,in&
&dex_list(size_ee+1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaad61,.false.,.f&
&alse.,which_en_pair,enbasis,enbasis1,grad_enbasis1,xyzzyaabl61,xyzzya&
&abp61,xyzzyaaaz61,xyzzyaabd61)
xyzzyaabv61=xyzzyaaby61*xyzzyaaay61*xyzzyaaba61*xyzzyaaax61*xyzzyaaaz6&
&1
jas=jas+xyzzyaabv61
gjas=gjas+xyzzyaabv61*(xyzzyaabc61+xyzzyaabe61+xyzzyaabb61+xyzzyaabd61&
&)
else
if(mask_leftmost_change_ee(xyzzyaaaa61)>0)call xyzzyaaiv1(nfn_eebasis,&
&netot,rank_e,rank_e,size_ee,mask_leftmost_change_ee(xyzzyaaaa61),xyzz&
&yaaae61,index_list(1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaac61,.true.,.&
&true.,which_ee_pair,eecut,eecut1,xyzzyaabk61,xyzzyaaay61)
if(mask_leftmost_change_en(xyzzyaaaa61)>0)call xyzzyaaiv1(nfn_enbasis,&
&nitot,rank_e,rank_n,size_en,mask_leftmost_change_en(xyzzyaaaa61),xyzz&
&yaaaf61,index_list(size_ee+1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaad61,&
&.false.,.true.,which_en_pair,encut,encut1,xyzzyaabm61,xyzzyaaba61)
if(leftmost_change_ee(xyzzyaaaa61)>0)call xyzzyaaiv1(nfn_eebasis,netot&
&,rank_e,rank_e,size_ee,leftmost_change_ee(xyzzyaaaa61),xyzzyaaae61,in&
&dex_list(1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaac61,.true.,.false.,whi&
&ch_ee_pair,eebasis,eebasis1,xyzzyaabj61,xyzzyaaax61)
if(leftmost_change_en(xyzzyaaaa61)>0)call xyzzyaaiv1(nfn_enbasis,nitot&
&,rank_e,rank_n,size_en,leftmost_change_en(xyzzyaaaa61),xyzzyaaaf61,in&
&dex_list(size_ee+1,xyzzyaaak61),ii,xyzzyaaac61,xyzzyaaad61,.false.,.f&
&alse.,which_en_pair,enbasis,enbasis1,xyzzyaabl61,xyzzyaaaz61)
jas=jas+xyzzyaaby61*xyzzyaaay61*xyzzyaaba61*xyzzyaaax61*xyzzyaaaz61
endif
enddo
enddo
enddo
end subroutine xyzzyaaif1
subroutine xyzzyaaig1(ii,rank_e,size_ee,index_list,nparam_remap,index_&
&remap,leftmost_change_ee,mask_leftmost_change_ee,uncluster_ee,groups_&
&ee,param,nchannel,chsig,iparam0,which_ee_pair,symm_ee,eebasis1,eecut1&
&,nzeecut1,eebasis,eecut,nzeecut,grad_eebasis1,grad_eecut1,lap_eebasis&
&1,lap_eecut1,jas,gjas,ljas,djas,dgjas,dljas)
implicit none
integer,intent(in) :: ii,rank_e,size_ee,nchannel,groups_ee(nspin,nspin&
&)
integer,intent(in) :: index_list(size_ee,*),nparam_remap(nchannel),ind&
&ex_remap(*),leftmost_change_ee(*),mask_leftmost_change_ee(*)
logical,intent(in) :: nzeecut1(netot),nzeecut(netot,netot),uncluster_e&
&e
integer,intent(in) :: chsig(xyzzyaaae1,*),iparam0(*),which_ee_pair(2,*&
&),symm_ee
real(dp),intent(in) :: eebasis1(nfn_eebasis,*),eecut1(nfn_eebasis,*), &
&         eebasis(nfn_eebasis,netot,*),eecut(nfn_eebasis,netot,*),para&
&m(*)
real(dp),intent(in),optional :: grad_eebasis1(3,nfn_eebasis,*),grad_ee&
&cut1(3,nfn_eebasis,*),lap_eebasis1(nfn_eebasis,*),lap_eecut1(nfn_eeba&
&sis,*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(xyzzyaaad1),d&
&gjas(3,xyzzyaaad1),dljas(xyzzyaaad1)
integer xyzzyaaaa62,xyzzyaaab62(rank_e),xyzzyaaac62(rank_e),xyzzyaaad6&
&2(size_ee),xyzzyaaae62(size_ee),xyzzyaaaf62,xyzzyaaag62,xyzzyaaah62,x&
&yzzyaaai62,xyzzyaaaj62,xyzzyaaak62,xyzzyaaal62(rank_e)
integer xyzzyaaam62(netot,rank_e),xyzzyaaan62(rank_e),xyzzyaaao62(rank&
&_e)
real(dp) xyzzyaaap62,xyzzyaaaq62,xyzzyaaar62(3),xyzzyaaas62(3),xyzzyaa&
&at62,xyzzyaaau62,xyzzyaaav62(size_ee),xyzzyaaaw62(size_ee),xyzzyaaax6&
&2(3,size_ee),xyzzyaaay62(3,size_ee),xyzzyaaaz62(size_ee),xyzzyaaba62(&
&size_ee),xyzzyaabb62,xyzzyaabc62(3),xyzzyaabd62,xyzzyaabe62
logical xyzzyaabf62,xyzzyaabg62,xyzzyaabh62,xyzzyaabi62,xyzzyaabj62
xyzzyaabi62=present(dljas)
xyzzyaabf62=.true.
if(xyzzyaabi62.or.present(ljas))then
xyzzyaabg62=.true.
xyzzyaabh62=.true.
elseif(present(gjas))then
xyzzyaabg62=.true.
xyzzyaabh62=.false.
else
xyzzyaabg62=.false.
xyzzyaabh62=.false.
endif
if(xyzzyaabf62)jas=0.d0
if(xyzzyaabg62)then
gjas=0.d0
if(xyzzyaabh62)ljas=0.d0
endif
if(xyzzyaabi62)then
if(xyzzyaabf62)djas=0.d0
if(xyzzyaabg62)then
dgjas=0.d0
if(xyzzyaabh62)dljas=0.d0
endif
endif
xyzzyaaal62=0
xyzzyaaac62(1)=0
do while(iterate_electron_indices_fix(rank_e,netot,ii,nzeecut1,nzeecut&
&,uncluster_ee,xyzzyaaam62,xyzzyaaan62,xyzzyaaao62,xyzzyaaac62))
do xyzzyaaag62=1,rank_e
xyzzyaaab62(xyzzyaaag62)=which_spin(xyzzyaaac62(xyzzyaaag62))
enddo
xyzzyaabj62=all(xyzzyaaab62==xyzzyaaal62)
xyzzyaaal62=xyzzyaaab62
if(.not.xyzzyaabj62)then
call get_sig_ee_only(rank_e,groups_ee,xyzzyaaab62,xyzzyaaad62,xyzzyaaa&
&e62)
call match_signature(size_ee,xyzzyaaae1,nchannel,xyzzyaaad62,chsig,xyz&
&zyaaah62)
endif
if(xyzzyaaah62<1)cycle
if(symm_ee==-1)call insertion_point(rank_e,xyzzyaaac62,xyzzyaaak62)
xyzzyaaap62=1.d0
xyzzyaaaq62=1.d0
xyzzyaaar62=0.d0
xyzzyaaas62=0.d0
xyzzyaaat62=0.d0
xyzzyaaau62=0.d0
xyzzyaaaj62=iparam0(xyzzyaaah62)
do xyzzyaaaf62=1,nparam_remap(xyzzyaaah62)
xyzzyaaaa62=index_remap(xyzzyaaaj62+xyzzyaaaf62)
xyzzyaaai62=xyzzyaaaa62-xyzzyaaaj62
xyzzyaabe62=param(xyzzyaaaa62)
if(symm_ee==-1)call apply_insertion_sign(xyzzyaaak62,xyzzyaaae62,index&
&_list(1,xyzzyaaai62),xyzzyaabe62)
if(xyzzyaabh62)then
if(mask_leftmost_change_ee(xyzzyaaaa62)>0)call xyzzyaait1(nfn_eebasis,&
&netot,rank_e,rank_e,size_ee,mask_leftmost_change_ee(xyzzyaaaa62),xyzz&
&yaaae62,index_list(1,xyzzyaaai62),ii,xyzzyaaac62,xyzzyaaac62,.true.,.&
&true.,which_ee_pair,eecut,eecut1,grad_eecut1,lap_eecut1,xyzzyaaaw62,x&
&yzzyaaay62,xyzzyaaba62,xyzzyaaaq62,xyzzyaaas62,xyzzyaaau62)
if(leftmost_change_ee(xyzzyaaaa62)>0)call xyzzyaait1(nfn_eebasis,netot&
&,rank_e,rank_e,size_ee,leftmost_change_ee(xyzzyaaaa62),xyzzyaaae62,in&
&dex_list(1,xyzzyaaai62),ii,xyzzyaaac62,xyzzyaaac62,.true.,.false.,whi&
&ch_ee_pair,eebasis,eebasis1,grad_eebasis1,lap_eebasis1,xyzzyaaav62,xy&
&zzyaaax62,xyzzyaaaz62,xyzzyaaap62,xyzzyaaar62,xyzzyaaat62)
xyzzyaabb62=xyzzyaaaq62*xyzzyaaap62
xyzzyaabc62=xyzzyaaas62+xyzzyaaar62
xyzzyaabd62=xyzzyaaau62+xyzzyaaat62+ddot(3,xyzzyaabc62,1,xyzzyaabc62,1&
&)
if(xyzzyaabi62)then
djas(xyzzyaaaa62)=djas(xyzzyaaaa62)+xyzzyaabb62
dgjas(:,xyzzyaaaa62)=dgjas(:,xyzzyaaaa62)+xyzzyaabb62*xyzzyaabc62
dljas(xyzzyaaaa62)=dljas(xyzzyaaaa62)+xyzzyaabb62*xyzzyaabd62
endif
xyzzyaabb62=xyzzyaabb62*xyzzyaabe62
jas=jas+xyzzyaabb62
gjas=gjas+xyzzyaabb62*xyzzyaabc62
ljas=ljas+xyzzyaabb62*xyzzyaabd62
elseif(xyzzyaabg62)then
if(mask_leftmost_change_ee(xyzzyaaaa62)>0)call xyzzyaaiu1(nfn_eebasis,&
&netot,rank_e,rank_e,size_ee,mask_leftmost_change_ee(xyzzyaaaa62),xyzz&
&yaaae62,index_list(1,xyzzyaaai62),ii,xyzzyaaac62,xyzzyaaac62,.true.,.&
&true.,which_ee_pair,eecut,eecut1,grad_eecut1,xyzzyaaaw62,xyzzyaaay62,&
&xyzzyaaaq62,xyzzyaaas62)
if(leftmost_change_ee(xyzzyaaaa62)>0)call xyzzyaaiu1(nfn_eebasis,netot&
&,rank_e,rank_e,size_ee,leftmost_change_ee(xyzzyaaaa62),xyzzyaaae62,in&
&dex_list(1,xyzzyaaai62),ii,xyzzyaaac62,xyzzyaaac62,.true.,.false.,whi&
&ch_ee_pair,eebasis,eebasis1,grad_eebasis1,xyzzyaaav62,xyzzyaaax62,xyz&
&zyaaap62,xyzzyaaar62)
xyzzyaabb62=xyzzyaabe62*xyzzyaaaq62*xyzzyaaap62
jas=jas+xyzzyaabb62
gjas=gjas+xyzzyaabb62*(xyzzyaaas62+xyzzyaaar62)
else
if(mask_leftmost_change_ee(xyzzyaaaa62)>0)call xyzzyaaiv1(nfn_eebasis,&
&netot,rank_e,rank_e,size_ee,mask_leftmost_change_ee(xyzzyaaaa62),xyzz&
&yaaae62,index_list(1,xyzzyaaai62),ii,xyzzyaaac62,xyzzyaaac62,.true.,.&
&true.,which_ee_pair,eecut,eecut1,xyzzyaaaw62,xyzzyaaaq62)
if(leftmost_change_ee(xyzzyaaaa62)>0)call xyzzyaaiv1(nfn_eebasis,netot&
&,rank_e,rank_e,size_ee,leftmost_change_ee(xyzzyaaaa62),xyzzyaaae62,in&
&dex_list(1,xyzzyaaai62),ii,xyzzyaaac62,xyzzyaaac62,.true.,.false.,whi&
&ch_ee_pair,eebasis,eebasis1,xyzzyaaav62,xyzzyaaap62)
jas=jas+xyzzyaabe62*xyzzyaaaq62*xyzzyaaap62
endif
enddo
enddo
end subroutine xyzzyaaig1
subroutine xyzzyaaih1(ii,order_ee,groups_ee1,param,iparam0,eebasis1,ee&
&cut1,nzeecut1,grad_eebasis1,grad_eecut1,lap_eebasis1,lap_eecut1,jas,g&
&jas,ljas,djas,dgjas,dljas)
implicit none
integer,intent(in) :: ii,order_ee,groups_ee1(nspin)
integer,intent(in) :: iparam0(*)
real(dp),intent(in) :: eebasis1(nfn_eebasis,*),eecut1(nfn_eebasis,*),p&
&aram(*)
real(dp),intent(in),optional :: grad_eebasis1(3,nfn_eebasis,*),grad_ee&
&cut1(3,nfn_eebasis,*),lap_eebasis1(nfn_eebasis,*),lap_eecut1(nfn_eeba&
&sis,*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(xyzzyaaad1),d&
&gjas(3,xyzzyaaad1),dljas(xyzzyaaad1)
logical,intent(in) :: nzeecut1(netot)
integer xyzzyaaaa63,xyzzyaaab63,xyzzyaaac63,xyzzyaaad63
real(dp) xyzzyaaae63,xyzzyaaaf63(3),xyzzyaaag63
logical xyzzyaaah63,xyzzyaaai63,xyzzyaaaj63,xyzzyaaak63
xyzzyaaak63=present(dljas)
xyzzyaaah63=.true.
if(xyzzyaaak63.or.present(ljas))then
xyzzyaaai63=.true.
xyzzyaaaj63=.true.
elseif(present(gjas))then
xyzzyaaai63=.true.
xyzzyaaaj63=.false.
else
xyzzyaaai63=.false.
xyzzyaaaj63=.false.
endif
if(xyzzyaaah63)jas=0.d0
if(xyzzyaaai63)then
gjas=0.d0
if(xyzzyaaaj63)ljas=0.d0
endif
if(xyzzyaaak63)then
if(xyzzyaaah63)djas=0.d0
if(xyzzyaaai63)then
dgjas=0.d0
if(xyzzyaaaj63)dljas=0.d0
endif
endif
do xyzzyaaaa63=1,netot
if(ii==xyzzyaaaa63.or..not.nzeecut1(xyzzyaaaa63))cycle
xyzzyaaac63=groups_ee1(which_spin(xyzzyaaaa63))
if(xyzzyaaac63<1)cycle
xyzzyaaab63=iparam0(xyzzyaaac63)+1
xyzzyaaae63=ddot(order_ee,param(xyzzyaaab63),1,eebasis1(1,xyzzyaaaa63)&
&,1)
jas=jas+xyzzyaaae63*eecut1(1,xyzzyaaaa63)
if(xyzzyaaai63)then
xyzzyaaaf63(:)=param(xyzzyaaab63)*grad_eebasis1(:,1,xyzzyaaaa63)
do xyzzyaaad63=2,order_ee
xyzzyaaaf63(:)=xyzzyaaaf63(:)+param(xyzzyaaab63+xyzzyaaad63-1)*grad_ee&
&basis1(:,xyzzyaaad63,xyzzyaaaa63)
enddo
gjas=gjas+xyzzyaaae63*grad_eecut1(1:3,1,xyzzyaaaa63)+xyzzyaaaf63*eecut&
&1(1,xyzzyaaaa63)
if(xyzzyaaaj63)then
xyzzyaaag63=ddot(order_ee,param(xyzzyaaab63),1,lap_eebasis1(1,xyzzyaaa&
&a63),1)
ljas=ljas+xyzzyaaae63*lap_eecut1(1,xyzzyaaaa63)+2*ddot(3,xyzzyaaaf63(1&
&),1,grad_eecut1(1,1,xyzzyaaaa63),1)+xyzzyaaag63*eecut1(1,xyzzyaaaa63)
if(xyzzyaaak63)then
djas(xyzzyaaab63:xyzzyaaab63+order_ee-1)=djas(xyzzyaaab63:xyzzyaaab63+&
&order_ee-1)+eebasis1(1:order_ee,xyzzyaaaa63)*eecut1(1,xyzzyaaaa63)
dgjas(1,xyzzyaaab63:xyzzyaaab63+order_ee-1)=dgjas(1,xyzzyaaab63:xyzzya&
&aab63+order_ee-1)+eebasis1(1:order_ee,xyzzyaaaa63)*grad_eecut1(1,1,xy&
&zzyaaaa63)+grad_eebasis1(1,1:order_ee,xyzzyaaaa63)*eecut1(1,xyzzyaaaa&
&63)
dgjas(2,xyzzyaaab63:xyzzyaaab63+order_ee-1)=dgjas(2,xyzzyaaab63:xyzzya&
&aab63+order_ee-1)+eebasis1(1:order_ee,xyzzyaaaa63)*grad_eecut1(2,1,xy&
&zzyaaaa63)+grad_eebasis1(2,1:order_ee,xyzzyaaaa63)*eecut1(1,xyzzyaaaa&
&63)
dgjas(3,xyzzyaaab63:xyzzyaaab63+order_ee-1)=dgjas(3,xyzzyaaab63:xyzzya&
&aab63+order_ee-1)+eebasis1(1:order_ee,xyzzyaaaa63)*grad_eecut1(3,1,xy&
&zzyaaaa63)+grad_eebasis1(3,1:order_ee,xyzzyaaaa63)*eecut1(1,xyzzyaaaa&
&63)
dljas(xyzzyaaab63:xyzzyaaab63+order_ee-1)=dljas(xyzzyaaab63:xyzzyaaab6&
&3+order_ee-1)+eebasis1(1:order_ee,xyzzyaaaa63)*lap_eecut1(1,xyzzyaaaa&
&63)+2*grad_eebasis1(1,1:order_ee,xyzzyaaaa63)*grad_eecut1(1,1,xyzzyaa&
&aa63)+2*grad_eebasis1(2,1:order_ee,xyzzyaaaa63)*grad_eecut1(2,1,xyzzy&
&aaaa63)+2*grad_eebasis1(3,1:order_ee,xyzzyaaaa63)*grad_eecut1(3,1,xyz&
&zyaaaa63)+lap_eebasis1(1:order_ee,xyzzyaaaa63)*eecut1(1,xyzzyaaaa63)
endif
endif
endif
enddo
end subroutine xyzzyaaih1
subroutine xyzzyaaii1(ii,order_ee,groups_ee1,param,iparam0,eebasis1,ee&
&cut1,nzeecut1,deebasis1,deecut1,gradr_eebasis1,d2eebasis1,d2eecut1,la&
&pr_eebasis1,jas,gjas,ljas,djas,dgjas,dljas)
implicit none
integer,intent(in) :: ii,order_ee,groups_ee1(nspin)
integer,intent(in) :: iparam0(*)
real(dp),intent(in) :: eebasis1(nfn_eebasis,*),eecut1(nfn_eebasis,*),p&
&aram(*)
real(dp),intent(in),optional :: deebasis1(nfn_eebasis,*),deecut1(nfn_e&
&ebasis,*),gradr_eebasis1(3,*),d2eebasis1(nfn_eebasis,*),d2eecut1(nfn_&
&eebasis,*),lapr_eebasis1(*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(xyzzyaaad1),d&
&gjas(3,xyzzyaaad1),dljas(xyzzyaaad1)
logical,intent(in) :: nzeecut1(netot)
integer xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64
real(dp) xyzzyaaad64,xyzzyaaae64,xyzzyaaaf64,xyzzyaaag64,xyzzyaaah64,x&
&yzzyaaai64,xyzzyaaaj64,xyzzyaaak64(order_ee)
logical xyzzyaaal64,xyzzyaaam64,xyzzyaaan64,xyzzyaaao64
xyzzyaaao64=present(dljas)
xyzzyaaal64=.true.
if(xyzzyaaao64.or.present(ljas))then
xyzzyaaam64=.true.
xyzzyaaan64=.true.
elseif(present(gjas))then
xyzzyaaam64=.true.
xyzzyaaan64=.false.
else
xyzzyaaam64=.false.
xyzzyaaan64=.false.
endif
if(xyzzyaaal64)jas=0.d0
if(xyzzyaaam64)then
gjas=0.d0
if(xyzzyaaan64)ljas=0.d0
endif
if(xyzzyaaao64)then
if(xyzzyaaal64)djas=0.d0
if(xyzzyaaam64)then
dgjas=0.d0
if(xyzzyaaan64)dljas=0.d0
endif
endif
do xyzzyaaaa64=1,netot
if(ii==xyzzyaaaa64.or..not.nzeecut1(xyzzyaaaa64))cycle
xyzzyaaac64=groups_ee1(which_spin(xyzzyaaaa64))
if(xyzzyaaac64<1)cycle
xyzzyaaab64=iparam0(xyzzyaaac64)+1
xyzzyaaad64=ddot(order_ee,param(xyzzyaaab64),1,eebasis1(1,xyzzyaaaa64)&
&,1)
xyzzyaaag64=eecut1(1,xyzzyaaaa64)
jas=jas+xyzzyaaad64*xyzzyaaag64
if(xyzzyaaam64)then
xyzzyaaae64=ddot(order_ee,param(xyzzyaaab64),1,deebasis1(1,xyzzyaaaa64&
&),1)
xyzzyaaah64=deecut1(1,xyzzyaaaa64)
xyzzyaaaj64=xyzzyaaae64*xyzzyaaag64+xyzzyaaad64*xyzzyaaah64
gjas=gjas+xyzzyaaaj64*gradr_eebasis1(:,xyzzyaaaa64)
if(xyzzyaaan64)then
xyzzyaaaf64=ddot(order_ee,param(xyzzyaaab64),1,d2eebasis1(1,xyzzyaaaa6&
&4),1)
xyzzyaaai64=d2eecut1(1,xyzzyaaaa64)
ljas=ljas+xyzzyaaaf64*xyzzyaaag64+2*xyzzyaaae64*xyzzyaaah64+xyzzyaaad6&
&4*xyzzyaaai64+xyzzyaaaj64*lapr_eebasis1(xyzzyaaaa64)
if(xyzzyaaao64)then
djas(xyzzyaaab64:xyzzyaaab64+order_ee-1)=djas(xyzzyaaab64:xyzzyaaab64+&
&order_ee-1)+eebasis1(1:order_ee,xyzzyaaaa64)*xyzzyaaag64
xyzzyaaak64(1:order_ee)=eebasis1(1:order_ee,xyzzyaaaa64)*xyzzyaaah64+d&
&eebasis1(1:order_ee,xyzzyaaaa64)*xyzzyaaag64
dgjas(1,xyzzyaaab64:xyzzyaaab64+order_ee-1)=dgjas(1,xyzzyaaab64:xyzzya&
&aab64+order_ee-1)+xyzzyaaak64(1:order_ee)*gradr_eebasis1(1,xyzzyaaaa6&
&4)
dgjas(2,xyzzyaaab64:xyzzyaaab64+order_ee-1)=dgjas(2,xyzzyaaab64:xyzzya&
&aab64+order_ee-1)+xyzzyaaak64(1:order_ee)*gradr_eebasis1(2,xyzzyaaaa6&
&4)
dgjas(3,xyzzyaaab64:xyzzyaaab64+order_ee-1)=dgjas(3,xyzzyaaab64:xyzzya&
&aab64+order_ee-1)+xyzzyaaak64(1:order_ee)*gradr_eebasis1(3,xyzzyaaaa6&
&4)
dljas(xyzzyaaab64:xyzzyaaab64+order_ee-1)=dljas(xyzzyaaab64:xyzzyaaab6&
&4+order_ee-1)+eebasis1(1:order_ee,xyzzyaaaa64)*xyzzyaaai64+2*deebasis&
&1(1:order_ee,xyzzyaaaa64)*xyzzyaaah64+d2eebasis1(1:order_ee,xyzzyaaaa&
&64)*xyzzyaaag64+xyzzyaaak64(1:order_ee)*lapr_eebasis1(xyzzyaaaa64)
endif
endif
endif
enddo
end subroutine xyzzyaaii1
subroutine xyzzyaaij1(ii,rank_n,size_en,index_list,nparam_remap,index_&
&remap,leftmost_change_en,mask_leftmost_change_en,uncluster_en,groups_&
&en,param,nchannel,chsig,iparam0,which_en_pair,enbasis1,encut1,nzencut&
&1,grad_enbasis1,grad_encut1,lap_enbasis1,lap_encut1,jas,gjas,ljas,dja&
&s,dgjas,dljas)
implicit none
integer,intent(in) :: ii,rank_n,size_en,nchannel,groups_en(nitot,nspin&
&)
integer,intent(in) :: index_list(size_en,*),nparam_remap(nchannel),ind&
&ex_remap(*),leftmost_change_en(*),mask_leftmost_change_en(*)
integer,intent(in) :: chsig(xyzzyaaae1,*),iparam0(*),which_en_pair(2,*&
&)
real(dp),intent(in) :: enbasis1(nfn_enbasis,*),encut1(nfn_enbasis,*),p&
&aram(*)
real(dp),intent(in),optional :: grad_enbasis1(3,nfn_enbasis,*),grad_en&
&cut1(3,nfn_enbasis,*),lap_enbasis1(nfn_enbasis,*),lap_encut1(nfn_enba&
&sis,*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(xyzzyaaad1),d&
&gjas(3,xyzzyaaad1),dljas(xyzzyaaad1)
logical,intent(in) :: nzencut1(nitot),uncluster_en
integer xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65(rank_n),xyzzyaaad65(size_e&
&n),xyzzyaaae65(size_en),xyzzyaaaf65,xyzzyaaag65,xyzzyaaah65,xyzzyaaai&
&65,xyzzyaaaj65(netot,rank_n),xyzzyaaak65(nitot,rank_n),xyzzyaaal65(ra&
&nk_n),xyzzyaaam65(rank_n),xyzzyaaan65(rank_n)
real(dp) xyzzyaaao65,xyzzyaaap65,xyzzyaaaq65(3),xyzzyaaar65(3),xyzzyaa&
&as65,xyzzyaaat65,xyzzyaaau65(size_en),xyzzyaaav65(size_en),xyzzyaaaw6&
&5(3,size_en),xyzzyaaax65(3,size_en),xyzzyaaay65(size_en),xyzzyaaaz65(&
&size_en),xyzzyaaba65,xyzzyaabb65(3),xyzzyaabc65,xyzzyaabd65
logical xyzzyaabe65,xyzzyaabf65,xyzzyaabg65,xyzzyaabh65
xyzzyaabh65=present(dljas)
xyzzyaabe65=.true.
if(xyzzyaabh65.or.present(ljas))then
xyzzyaabf65=.true.
xyzzyaabg65=.true.
elseif(present(gjas))then
xyzzyaabf65=.true.
xyzzyaabg65=.false.
else
xyzzyaabf65=.false.
xyzzyaabg65=.false.
endif
if(xyzzyaabe65)jas=0.d0
if(xyzzyaabf65)then
gjas=0.d0
if(xyzzyaabg65)ljas=0.d0
endif
if(xyzzyaabh65)then
if(xyzzyaabe65)djas=0.d0
if(xyzzyaabf65)then
dgjas=0.d0
if(xyzzyaabg65)dljas=0.d0
endif
endif
xyzzyaaaa65=which_spin(ii)
xyzzyaaac65(1)=0
do while(iterate_nuclei_indices_fix(1,rank_n,netot,nitot,ii,nzencut1,n&
&zencut1,uncluster_en,xyzzyaaaj65,xyzzyaaak65,xyzzyaaal65,xyzzyaaam65,&
&xyzzyaaan65,xyzzyaaac65))
call get_sig_en_only(rank_n,groups_en,xyzzyaaaa65,xyzzyaaac65,xyzzyaaa&
&d65,xyzzyaaae65)
call match_signature(size_en,xyzzyaaae1,nchannel,xyzzyaaad65,chsig,xyz&
&zyaaaf65)
if(xyzzyaaaf65<1)cycle
xyzzyaaao65=1.d0
xyzzyaaap65=1.d0
xyzzyaaaq65=0.d0
xyzzyaaar65=0.d0
xyzzyaaas65=0.d0
xyzzyaaat65=0.d0
xyzzyaaai65=iparam0(xyzzyaaaf65)
do xyzzyaaag65=1,nparam_remap(xyzzyaaaf65)
xyzzyaaab65=index_remap(xyzzyaaai65+xyzzyaaag65)
xyzzyaaah65=xyzzyaaab65-xyzzyaaai65
xyzzyaabd65=param(xyzzyaaab65)
if(xyzzyaabg65)then
if(mask_leftmost_change_en(xyzzyaaab65)>0)call xyzzyaait1(nfn_enbasis,&
&nitot,1,rank_n,size_en,mask_leftmost_change_en(xyzzyaaab65),xyzzyaaae&
&65,index_list(1,xyzzyaaah65),ii,(/ii/),xyzzyaaac65,.false.,.true.,whi&
&ch_en_pair,encut1,encut1,grad_encut1,lap_encut1,xyzzyaaav65,xyzzyaaax&
&65,xyzzyaaaz65,xyzzyaaap65,xyzzyaaar65,xyzzyaaat65)
if(leftmost_change_en(xyzzyaaab65)>0)call xyzzyaait1(nfn_enbasis,nitot&
&,1,rank_n,size_en,leftmost_change_en(xyzzyaaab65),xyzzyaaae65,index_l&
&ist(1,xyzzyaaah65),ii,(/ii/),xyzzyaaac65,.false.,.false.,which_en_pai&
&r,enbasis1,enbasis1,grad_enbasis1,lap_enbasis1,xyzzyaaau65,xyzzyaaaw6&
&5,xyzzyaaay65,xyzzyaaao65,xyzzyaaaq65,xyzzyaaas65)
xyzzyaaba65=xyzzyaaap65*xyzzyaaao65
xyzzyaabb65=xyzzyaaar65+xyzzyaaaq65
xyzzyaabc65=xyzzyaaat65+xyzzyaaas65+ddot(3,xyzzyaabb65,1,xyzzyaabb65,1&
&)
if(xyzzyaabh65)then
djas(xyzzyaaab65)=djas(xyzzyaaab65)+xyzzyaaba65
dgjas(:,xyzzyaaab65)=dgjas(:,xyzzyaaab65)+xyzzyaaba65*xyzzyaabb65
dljas(xyzzyaaab65)=dljas(xyzzyaaab65)+xyzzyaaba65*xyzzyaabc65
endif
xyzzyaaba65=xyzzyaaba65*xyzzyaabd65
jas=jas+xyzzyaaba65
gjas=gjas+xyzzyaaba65*xyzzyaabb65
ljas=ljas+xyzzyaaba65*xyzzyaabc65
elseif(xyzzyaabf65)then
if(mask_leftmost_change_en(xyzzyaaab65)>0)call xyzzyaaiu1(nfn_enbasis,&
&nitot,1,rank_n,size_en,mask_leftmost_change_en(xyzzyaaab65),xyzzyaaae&
&65,index_list(1,xyzzyaaah65),ii,(/ii/),xyzzyaaac65,.false.,.true.,whi&
&ch_en_pair,encut1,encut1,grad_encut1,xyzzyaaav65,xyzzyaaax65,xyzzyaaa&
&p65,xyzzyaaar65)
if(leftmost_change_en(xyzzyaaab65)>0)call xyzzyaaiu1(nfn_enbasis,nitot&
&,1,rank_n,size_en,leftmost_change_en(xyzzyaaab65),xyzzyaaae65,index_l&
&ist(1,xyzzyaaah65),ii,(/ii/),xyzzyaaac65,.false.,.false.,which_en_pai&
&r,enbasis1,enbasis1,grad_enbasis1,xyzzyaaau65,xyzzyaaaw65,xyzzyaaao65&
&,xyzzyaaaq65)
xyzzyaaba65=xyzzyaabd65*xyzzyaaap65*xyzzyaaao65
jas=jas+xyzzyaaba65
gjas=gjas+xyzzyaaba65*(xyzzyaaar65+xyzzyaaaq65)
else
if(mask_leftmost_change_en(xyzzyaaab65)>0)call xyzzyaaiv1(nfn_enbasis,&
&nitot,1,rank_n,size_en,mask_leftmost_change_en(xyzzyaaab65),xyzzyaaae&
&65,index_list(1,xyzzyaaah65),ii,(/ii/),xyzzyaaac65,.false.,.true.,whi&
&ch_en_pair,encut1,encut1,xyzzyaaav65,xyzzyaaap65)
if(leftmost_change_en(xyzzyaaab65)>0)call xyzzyaaiv1(nfn_enbasis,nitot&
&,1,rank_n,size_en,leftmost_change_en(xyzzyaaab65),xyzzyaaae65,index_l&
&ist(1,xyzzyaaah65),ii,(/ii/),xyzzyaaac65,.false.,.false.,which_en_pai&
&r,enbasis1,enbasis1,xyzzyaaau65,xyzzyaaao65)
jas=jas+xyzzyaabd65*xyzzyaaap65*xyzzyaaao65
endif
enddo
enddo
end subroutine xyzzyaaij1
subroutine xyzzyaaik1(order_en,groups_en1,param,iparam0,enbasis1,encut&
&1,nzencut1,grad_enbasis1,grad_encut1,lap_enbasis1,lap_encut1,jas,gjas&
&,ljas,djas,dgjas,dljas)
implicit none
integer,intent(in) :: order_en,groups_en1(nitot)
integer,intent(in) :: iparam0(*)
real(dp),intent(in) :: enbasis1(nfn_enbasis,*),encut1(nfn_enbasis,*),p&
&aram(*)
real(dp),intent(in),optional :: grad_enbasis1(3,nfn_enbasis,*),grad_en&
&cut1(3,nfn_enbasis,*),lap_enbasis1(nfn_enbasis,*),lap_encut1(nfn_enba&
&sis,*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(xyzzyaaad1),d&
&gjas(3,xyzzyaaad1),dljas(xyzzyaaad1)
logical,intent(in) :: nzencut1(nitot)
integer xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66,xyzzyaaad66
real(dp) xyzzyaaae66,xyzzyaaaf66(3),xyzzyaaag66
logical xyzzyaaah66,xyzzyaaai66,xyzzyaaaj66,xyzzyaaak66
xyzzyaaak66=present(dljas)
xyzzyaaah66=.true.
if(xyzzyaaak66.or.present(ljas))then
xyzzyaaai66=.true.
xyzzyaaaj66=.true.
elseif(present(gjas))then
xyzzyaaai66=.true.
xyzzyaaaj66=.false.
else
xyzzyaaai66=.false.
xyzzyaaaj66=.false.
endif
if(xyzzyaaah66)jas=0.d0
if(xyzzyaaai66)then
gjas=0.d0
if(xyzzyaaaj66)ljas=0.d0
endif
if(xyzzyaaak66)then
if(xyzzyaaah66)djas=0.d0
if(xyzzyaaai66)then
dgjas=0.d0
if(xyzzyaaaj66)dljas=0.d0
endif
endif
do xyzzyaaab66=1,nitot
if(.not.nzencut1(xyzzyaaab66))cycle
xyzzyaaac66=groups_en1(xyzzyaaab66)
if(xyzzyaaac66<1)cycle
xyzzyaaaa66=iparam0(xyzzyaaac66)+1
xyzzyaaae66=ddot(order_en,param(xyzzyaaaa66),1,enbasis1(1,xyzzyaaab66)&
&,1)
jas=jas+xyzzyaaae66*encut1(1,xyzzyaaab66)
if(xyzzyaaai66)then
xyzzyaaaf66(:)=param(xyzzyaaaa66)*grad_enbasis1(:,1,xyzzyaaab66)
do xyzzyaaad66=2,order_en
xyzzyaaaf66(:)=xyzzyaaaf66(:)+param(xyzzyaaaa66+xyzzyaaad66-1)*grad_en&
&basis1(:,xyzzyaaad66,xyzzyaaab66)
enddo
gjas=gjas+xyzzyaaae66*grad_encut1(1:3,1,xyzzyaaab66)+xyzzyaaaf66*encut&
&1(1,xyzzyaaab66)
if(xyzzyaaaj66)then
xyzzyaaag66=ddot(order_en,param(xyzzyaaaa66),1,lap_enbasis1(1,xyzzyaaa&
&b66),1)
ljas=ljas+xyzzyaaae66*lap_encut1(1,xyzzyaaab66)+2*ddot(3,xyzzyaaaf66(1&
&),1,grad_encut1(1,1,xyzzyaaab66),1)+xyzzyaaag66*encut1(1,xyzzyaaab66)
if(xyzzyaaak66)then
djas(xyzzyaaaa66:xyzzyaaaa66+order_en-1)=djas(xyzzyaaaa66:xyzzyaaaa66+&
&order_en-1)+enbasis1(1:order_en,xyzzyaaab66)*encut1(1,xyzzyaaab66)
dgjas(1,xyzzyaaaa66:xyzzyaaaa66+order_en-1)=dgjas(1,xyzzyaaaa66:xyzzya&
&aaa66+order_en-1)+enbasis1(1:order_en,xyzzyaaab66)*grad_encut1(1,1,xy&
&zzyaaab66)+grad_enbasis1(1,1:order_en,xyzzyaaab66)*encut1(1,xyzzyaaab&
&66)
dgjas(2,xyzzyaaaa66:xyzzyaaaa66+order_en-1)=dgjas(2,xyzzyaaaa66:xyzzya&
&aaa66+order_en-1)+enbasis1(1:order_en,xyzzyaaab66)*grad_encut1(2,1,xy&
&zzyaaab66)+grad_enbasis1(2,1:order_en,xyzzyaaab66)*encut1(1,xyzzyaaab&
&66)
dgjas(3,xyzzyaaaa66:xyzzyaaaa66+order_en-1)=dgjas(3,xyzzyaaaa66:xyzzya&
&aaa66+order_en-1)+enbasis1(1:order_en,xyzzyaaab66)*grad_encut1(3,1,xy&
&zzyaaab66)+grad_enbasis1(3,1:order_en,xyzzyaaab66)*encut1(1,xyzzyaaab&
&66)
dljas(xyzzyaaaa66:xyzzyaaaa66+order_en-1)=dljas(xyzzyaaaa66:xyzzyaaaa6&
&6+order_en-1)+enbasis1(1:order_en,xyzzyaaab66)*lap_encut1(1,xyzzyaaab&
&66)+2*grad_enbasis1(1,1:order_en,xyzzyaaab66)*grad_encut1(1,1,xyzzyaa&
&ab66)+2*grad_enbasis1(2,1:order_en,xyzzyaaab66)*grad_encut1(2,1,xyzzy&
&aaab66)+2*grad_enbasis1(3,1:order_en,xyzzyaaab66)*grad_encut1(3,1,xyz&
&zyaaab66)+lap_enbasis1(1:order_en,xyzzyaaab66)*encut1(1,xyzzyaaab66)
endif
endif
endif
enddo
end subroutine xyzzyaaik1
subroutine xyzzyaail1(order_en,groups_en1,param,iparam0,enbasis1,encut&
&1,nzencut1,denbasis1,dencut1,gradr_enbasis1,d2enbasis1,d2encut1,lapr_&
&enbasis1,jas,gjas,ljas,djas,dgjas,dljas)
implicit none
integer,intent(in) :: order_en,groups_en1(nitot)
integer,intent(in) :: iparam0(*)
real(dp),intent(in) :: enbasis1(nfn_enbasis,*),encut1(nfn_enbasis,*),p&
&aram(*)
real(dp),intent(in),optional :: denbasis1(nfn_enbasis,*),dencut1(nfn_e&
&nbasis,*),gradr_enbasis1(3,*),d2enbasis1(nfn_enbasis,*),d2encut1(nfn_&
&enbasis,*),lapr_enbasis1(*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(xyzzyaaad1),d&
&gjas(3,xyzzyaaad1),dljas(xyzzyaaad1)
logical,intent(in) :: nzencut1(nitot)
integer xyzzyaaaa67,xyzzyaaab67,xyzzyaaac67
real(dp) xyzzyaaad67,xyzzyaaae67,xyzzyaaaf67,xyzzyaaag67,xyzzyaaah67,x&
&yzzyaaai67,xyzzyaaaj67,xyzzyaaak67(order_en)
logical xyzzyaaal67,xyzzyaaam67,xyzzyaaan67,xyzzyaaao67
xyzzyaaao67=present(dljas)
xyzzyaaal67=.true.
if(xyzzyaaao67.or.present(ljas))then
xyzzyaaam67=.true.
xyzzyaaan67=.true.
elseif(present(gjas))then
xyzzyaaam67=.true.
xyzzyaaan67=.false.
else
xyzzyaaam67=.false.
xyzzyaaan67=.false.
endif
if(xyzzyaaal67)jas=0.d0
if(xyzzyaaam67)then
gjas=0.d0
if(xyzzyaaan67)ljas=0.d0
endif
if(xyzzyaaao67)then
if(xyzzyaaal67)djas=0.d0
if(xyzzyaaam67)then
dgjas=0.d0
if(xyzzyaaan67)dljas=0.d0
endif
endif
do xyzzyaaab67=1,nitot
if(.not.nzencut1(xyzzyaaab67))cycle
xyzzyaaac67=groups_en1(xyzzyaaab67)
if(xyzzyaaac67<1)cycle
xyzzyaaaa67=iparam0(xyzzyaaac67)+1
xyzzyaaad67=ddot(order_en,param(xyzzyaaaa67),1,enbasis1(1,xyzzyaaab67)&
&,1)
xyzzyaaag67=encut1(1,xyzzyaaab67)
jas=jas+xyzzyaaad67*xyzzyaaag67
if(xyzzyaaam67)then
xyzzyaaae67=ddot(order_en,param(xyzzyaaaa67),1,denbasis1(1,xyzzyaaab67&
&),1)
xyzzyaaah67=dencut1(1,xyzzyaaab67)
xyzzyaaaj67=xyzzyaaae67*xyzzyaaag67+xyzzyaaad67*xyzzyaaah67
gjas=gjas+xyzzyaaaj67*gradr_enbasis1(:,xyzzyaaab67)
if(xyzzyaaan67)then
xyzzyaaaf67=ddot(order_en,param(xyzzyaaaa67),1,d2enbasis1(1,xyzzyaaab6&
&7),1)
xyzzyaaai67=d2encut1(1,xyzzyaaab67)
ljas=ljas+xyzzyaaaf67*xyzzyaaag67+2*xyzzyaaae67*xyzzyaaah67+xyzzyaaad6&
&7*xyzzyaaai67+xyzzyaaaj67*lapr_enbasis1(xyzzyaaab67)
if(xyzzyaaao67)then
djas(xyzzyaaaa67:xyzzyaaaa67+order_en-1)=djas(xyzzyaaaa67:xyzzyaaaa67+&
&order_en-1)+enbasis1(1:order_en,xyzzyaaab67)*xyzzyaaag67
xyzzyaaak67(1:order_en)=enbasis1(1:order_en,xyzzyaaab67)*xyzzyaaah67+d&
&enbasis1(1:order_en,xyzzyaaab67)*xyzzyaaag67
dgjas(1,xyzzyaaaa67:xyzzyaaaa67+order_en-1)=dgjas(1,xyzzyaaaa67:xyzzya&
&aaa67+order_en-1)+xyzzyaaak67(1:order_en)*gradr_enbasis1(1,xyzzyaaab6&
&7)
dgjas(2,xyzzyaaaa67:xyzzyaaaa67+order_en-1)=dgjas(2,xyzzyaaaa67:xyzzya&
&aaa67+order_en-1)+xyzzyaaak67(1:order_en)*gradr_enbasis1(2,xyzzyaaab6&
&7)
dgjas(3,xyzzyaaaa67:xyzzyaaaa67+order_en-1)=dgjas(3,xyzzyaaaa67:xyzzya&
&aaa67+order_en-1)+xyzzyaaak67(1:order_en)*gradr_enbasis1(3,xyzzyaaab6&
&7)
dljas(xyzzyaaaa67:xyzzyaaaa67+order_en-1)=dljas(xyzzyaaaa67:xyzzyaaaa6&
&7+order_en-1)+enbasis1(1:order_en,xyzzyaaab67)*xyzzyaaai67+2*denbasis&
&1(1:order_en,xyzzyaaab67)*xyzzyaaah67+d2enbasis1(1:order_en,xyzzyaaab&
&67)*xyzzyaaag67+xyzzyaaak67(1:order_en)*lapr_enbasis1(xyzzyaaab67)
endif
endif
endif
enddo
end subroutine xyzzyaail1
subroutine xyzzyaaim1(ii,order_ee,order_en,groups_ee1,groups_en,param,&
&nchannel,chsig,eebasis1,eecut1,nzeecut1,enbasis1,encut1,nzencut1,enba&
&sis,encut,nzencut,grad_eebasis1,grad_eecut1,grad_enbasis1,grad_encut1&
&,lap_eebasis1,lap_eecut1,lap_enbasis1,lap_encut1,jas,gjas,ljas,djas,d&
&gjas,dljas)
implicit none
integer,intent(in) :: ii,order_ee,order_en,nchannel,groups_ee1(nspin),&
&groups_en(nitot,nspin)
integer,intent(in) :: chsig(xyzzyaaae1,*)
real(dp),intent(in) :: eebasis1(nfn_eebasis,*),eecut1(nfn_eebasis,*),e&
&nbasis1(nfn_enbasis,*),encut1(nfn_enbasis,*),enbasis(nfn_enbasis,nito&
&t,*),encut(nfn_enbasis,nitot,*),param(order_en,order_en,order_ee,*)
real(dp),intent(in),optional :: grad_eebasis1(3,nfn_eebasis,*),grad_ee&
&cut1(3,nfn_eebasis,*),grad_enbasis1(3,nfn_enbasis,*),grad_encut1(3,nf&
&n_enbasis,*),lap_eebasis1(nfn_eebasis,*),lap_eecut1(nfn_eebasis,*),la&
&p_enbasis1(nfn_enbasis,*),lap_encut1(nfn_enbasis,*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(order_en,orde&
&r_en,order_ee,nchannel),dgjas(3,order_en,order_en,order_ee,nchannel),&
&dljas(order_en,order_en,order_ee,nchannel)
logical,intent(in) :: nzeecut1(netot),nzencut1(nitot),nzencut(nitot,ne&
&tot)
integer xyzzyaaaa68,xyzzyaaab68,xyzzyaaac68,xyzzyaaad68,xyzzyaaae68,xy&
&zzyaaaf68(3),xyzzyaaag68,xyzzyaaah68,xyzzyaaai68,xyzzyaaaj68,xyzzyaaa&
&k68
real(dp) xyzzyaaal68,xyzzyaaam68,xyzzyaaan68(3),xyzzyaaao68(3),xyzzyaa&
&ap68,xyzzyaaaq68,xyzzyaaar68,xyzzyaaas68,xyzzyaaat68(3),xyzzyaaau68,x&
&yzzyaaav68,xyzzyaaaw68(3),xyzzyaaax68,xyzzyaaay68,xyzzyaaaz68(3),xyzz&
&yaaba68,xyzzyaabb68,xyzzyaabc68(3),xyzzyaabd68,xyzzyaabe68,xyzzyaabf6&
&8(3),xyzzyaabg68,xyzzyaabh68
logical xyzzyaabi68,xyzzyaabj68,xyzzyaabk68,xyzzyaabl68,xyzzyaabm68
xyzzyaabm68=present(dljas)
xyzzyaabi68=.true.
if(xyzzyaabm68.or.present(ljas))then
xyzzyaabj68=.true.
xyzzyaabk68=.true.
elseif(present(gjas))then
xyzzyaabj68=.true.
xyzzyaabk68=.false.
else
xyzzyaabj68=.false.
xyzzyaabk68=.false.
endif
if(xyzzyaabi68)jas=0.d0
if(xyzzyaabj68)then
gjas=0.d0
if(xyzzyaabk68)ljas=0.d0
endif
if(xyzzyaabm68)then
if(xyzzyaabi68)djas=0.d0
if(xyzzyaabj68)then
dgjas=0.d0
if(xyzzyaabk68)dljas=0.d0
endif
endif
xyzzyaaaa68=which_spin(ii)
do xyzzyaaad68=1,nitot
if(.not.nzencut1(xyzzyaaad68))cycle
xyzzyaaaj68=0
do xyzzyaaab68=1,nspin
xyzzyaaaf68=(/groups_ee1(xyzzyaaab68),groups_en(xyzzyaaad68,xyzzyaaaa6&
&8),groups_en(xyzzyaaad68,xyzzyaaab68)/)
xyzzyaabl68=xyzzyaaaf68(2)>xyzzyaaaf68(3)
if(xyzzyaabl68)call swap1(xyzzyaaaf68(2),xyzzyaaaf68(3))
call match_signature(3,xyzzyaaae1,nchannel,xyzzyaaaf68,chsig,xyzzyaaae&
&68)
if(xyzzyaaae68<1)then
xyzzyaaaj68=xyzzyaaaj68+nele(xyzzyaaab68)
cycle
endif
xyzzyaaac68=xyzzyaaaj68
if(.not.xyzzyaabj68)then
if(.not.xyzzyaabl68)then
do xyzzyaaak68=1,nele(xyzzyaaab68)
xyzzyaaac68=xyzzyaaac68+1
if(xyzzyaaac68==ii)cycle
if(.not.nzeecut1(xyzzyaaac68).or..not.nzencut(xyzzyaaad68,xyzzyaaac68)&
&)cycle
xyzzyaaam68=eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaaad68)*encut(1,xyzzya&
&aad68,xyzzyaaac68)
xyzzyaaal68=0.d0
do xyzzyaaag68=1,order_ee
xyzzyaaas68=0.d0
do xyzzyaaah68=1,order_en
xyzzyaaas68=xyzzyaaas68+ddot(order_en,param(1,xyzzyaaah68,xyzzyaaag68,&
&xyzzyaaae68),1,enbasis(1,xyzzyaaad68,xyzzyaaac68),1)*enbasis1(xyzzyaa&
&ah68,xyzzyaaad68)
enddo
xyzzyaaal68=xyzzyaaal68+xyzzyaaas68*eebasis1(xyzzyaaag68,xyzzyaaac68)
enddo
jas=jas+xyzzyaaal68*xyzzyaaam68
enddo
else
do xyzzyaaak68=1,nele(xyzzyaaab68)
xyzzyaaac68=xyzzyaaac68+1
if(xyzzyaaac68==ii)cycle
if(.not.nzeecut1(xyzzyaaac68).or..not.nzencut(xyzzyaaad68,xyzzyaaac68)&
&)cycle
xyzzyaaam68=eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaaad68)*encut(1,xyzzya&
&aad68,xyzzyaaac68)
xyzzyaaal68=0.d0
do xyzzyaaag68=1,order_ee
xyzzyaaas68=0.d0
do xyzzyaaai68=1,order_en
xyzzyaaas68=xyzzyaaas68+ddot(order_en,param(1,xyzzyaaai68,xyzzyaaag68,&
&xyzzyaaae68),1,enbasis1(1,xyzzyaaad68),1)*enbasis(xyzzyaaai68,xyzzyaa&
&ad68,xyzzyaaac68)
enddo
xyzzyaaal68=xyzzyaaal68+xyzzyaaas68*eebasis1(xyzzyaaag68,xyzzyaaac68)
enddo
jas=jas+xyzzyaaal68*xyzzyaaam68
enddo
endif
elseif(.not.xyzzyaabm68)then
if(.not.xyzzyaabl68)then
do xyzzyaaak68=1,nele(xyzzyaaab68)
xyzzyaaac68=xyzzyaaac68+1
if(xyzzyaaac68==ii)cycle
if(.not.nzeecut1(xyzzyaaac68).or..not.nzencut(xyzzyaaad68,xyzzyaaac68)&
&)cycle
xyzzyaaam68=eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaaad68)*encut(1,xyzzya&
&aad68,xyzzyaaac68)
xyzzyaaao68=(grad_eecut1(:,1,xyzzyaaac68)*encut1(1,xyzzyaaad68)+eecut1&
&(1,xyzzyaaac68)*grad_encut1(:,1,xyzzyaaad68))*encut(1,xyzzyaaad68,xyz&
&zyaaac68)
if(xyzzyaabk68)xyzzyaaaq68=(lap_eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaa&
&ad68)+2*ddot(3,grad_eecut1(1,1,xyzzyaaac68),1,grad_encut1(1,1,xyzzyaa&
&ad68),1)+eecut1(1,xyzzyaaac68)*lap_encut1(1,xyzzyaaad68))*encut(1,xyz&
&zyaaad68,xyzzyaaac68)
xyzzyaaal68=0.d0
xyzzyaaan68=0.d0
if(xyzzyaabk68)xyzzyaaap68=0.d0
do xyzzyaaag68=1,order_ee
xyzzyaaas68=0.d0
xyzzyaaat68=0.d0
if(xyzzyaabk68)xyzzyaaau68=0.d0
do xyzzyaaah68=1,order_en
xyzzyaaar68=ddot(order_en,param(1,xyzzyaaah68,xyzzyaaag68,xyzzyaaae68)&
&,1,enbasis(1,xyzzyaaad68,xyzzyaaac68),1)
xyzzyaaas68=xyzzyaaas68+xyzzyaaar68*enbasis1(xyzzyaaah68,xyzzyaaad68)
xyzzyaaat68=xyzzyaaat68+xyzzyaaar68*grad_enbasis1(:,xyzzyaaah68,xyzzya&
&aad68)
if(xyzzyaabk68)xyzzyaaau68=xyzzyaaau68+xyzzyaaar68*lap_enbasis1(xyzzya&
&aah68,xyzzyaaad68)
enddo
xyzzyaaal68=xyzzyaaal68+xyzzyaaas68*eebasis1(xyzzyaaag68,xyzzyaaac68)
xyzzyaaan68=xyzzyaaan68+xyzzyaaat68*eebasis1(xyzzyaaag68,xyzzyaaac68)+&
&xyzzyaaas68*grad_eebasis1(:,xyzzyaaag68,xyzzyaaac68)
if(xyzzyaabk68)xyzzyaaap68=xyzzyaaap68+xyzzyaaau68*eebasis1(xyzzyaaag6&
&8,xyzzyaaac68)+2*ddot(3,xyzzyaaat68(1),1,grad_eebasis1(1,xyzzyaaag68,&
&xyzzyaaac68),1)+xyzzyaaas68*lap_eebasis1(xyzzyaaag68,xyzzyaaac68)
enddo
jas=jas+xyzzyaaal68*xyzzyaaam68
gjas=gjas+xyzzyaaan68*xyzzyaaam68+xyzzyaaal68*xyzzyaaao68
if(xyzzyaabk68)ljas=ljas+xyzzyaaap68*xyzzyaaam68+2*ddot(3,xyzzyaaan68(&
&1),1,xyzzyaaao68(1),1)+xyzzyaaal68*xyzzyaaaq68
enddo
else
do xyzzyaaak68=1,nele(xyzzyaaab68)
xyzzyaaac68=xyzzyaaac68+1
if(xyzzyaaac68==ii)cycle
if(.not.nzeecut1(xyzzyaaac68).or..not.nzencut(xyzzyaaad68,xyzzyaaac68)&
&)cycle
xyzzyaaam68=eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaaad68)*encut(1,xyzzya&
&aad68,xyzzyaaac68)
xyzzyaaao68=(grad_eecut1(:,1,xyzzyaaac68)*encut1(1,xyzzyaaad68)+eecut1&
&(1,xyzzyaaac68)*grad_encut1(:,1,xyzzyaaad68))*encut(1,xyzzyaaad68,xyz&
&zyaaac68)
if(xyzzyaabk68)xyzzyaaaq68=(lap_eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaa&
&ad68)+2*ddot(3,grad_eecut1(1,1,xyzzyaaac68),1,grad_encut1(1,1,xyzzyaa&
&ad68),1)+eecut1(1,xyzzyaaac68)*lap_encut1(1,xyzzyaaad68))*encut(1,xyz&
&zyaaad68,xyzzyaaac68)
xyzzyaaal68=0.d0
xyzzyaaan68=0.d0
if(xyzzyaabk68)xyzzyaaap68=0.d0
do xyzzyaaag68=1,order_ee
xyzzyaaas68=0.d0
xyzzyaaat68=0.d0
if(xyzzyaabk68)xyzzyaaau68=0.d0
do xyzzyaaah68=1,order_en
xyzzyaaar68=ddot(order_en,param(xyzzyaaah68,1,xyzzyaaag68,xyzzyaaae68)&
&,order_en,enbasis(1,xyzzyaaad68,xyzzyaaac68),1)
xyzzyaaas68=xyzzyaaas68+xyzzyaaar68*enbasis1(xyzzyaaah68,xyzzyaaad68)
xyzzyaaat68=xyzzyaaat68+xyzzyaaar68*grad_enbasis1(:,xyzzyaaah68,xyzzya&
&aad68)
if(xyzzyaabk68)xyzzyaaau68=xyzzyaaau68+xyzzyaaar68*lap_enbasis1(xyzzya&
&aah68,xyzzyaaad68)
enddo
xyzzyaaal68=xyzzyaaal68+xyzzyaaas68*eebasis1(xyzzyaaag68,xyzzyaaac68)
xyzzyaaan68=xyzzyaaan68+xyzzyaaat68*eebasis1(xyzzyaaag68,xyzzyaaac68)+&
&xyzzyaaas68*grad_eebasis1(:,xyzzyaaag68,xyzzyaaac68)
if(xyzzyaabk68)xyzzyaaap68=xyzzyaaap68+xyzzyaaau68*eebasis1(xyzzyaaag6&
&8,xyzzyaaac68)+2*ddot(3,xyzzyaaat68(1),1,grad_eebasis1(1,xyzzyaaag68,&
&xyzzyaaac68),1)+xyzzyaaas68*lap_eebasis1(xyzzyaaag68,xyzzyaaac68)
enddo
jas=jas+xyzzyaaal68*xyzzyaaam68
gjas=gjas+xyzzyaaan68*xyzzyaaam68+xyzzyaaal68*xyzzyaaao68
if(xyzzyaabk68)ljas=ljas+xyzzyaaap68*xyzzyaaam68+2*ddot(3,xyzzyaaan68(&
&1),1,xyzzyaaao68(1),1)+xyzzyaaal68*xyzzyaaaq68
enddo
endif
else
if(.not.xyzzyaabl68)then
do xyzzyaaak68=1,nele(xyzzyaaab68)
xyzzyaaac68=xyzzyaaac68+1
if(xyzzyaaac68==ii)cycle
if(.not.nzeecut1(xyzzyaaac68).or..not.nzencut(xyzzyaaad68,xyzzyaaac68)&
&)cycle
xyzzyaaam68=eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaaad68)*encut(1,xyzzya&
&aad68,xyzzyaaac68)
xyzzyaaao68=(grad_eecut1(:,1,xyzzyaaac68)*encut1(1,xyzzyaaad68)+eecut1&
&(1,xyzzyaaac68)*grad_encut1(:,1,xyzzyaaad68))*encut(1,xyzzyaaad68,xyz&
&zyaaac68)
xyzzyaaaq68=(lap_eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaaad68)+2*ddot(3,&
&grad_eecut1(1,1,xyzzyaaac68),1,grad_encut1(1,1,xyzzyaaad68),1)+eecut1&
&(1,xyzzyaaac68)*lap_encut1(1,xyzzyaaad68))*encut(1,xyzzyaaad68,xyzzya&
&aac68)
do xyzzyaaag68=1,order_ee
xyzzyaaav68=eebasis1(xyzzyaaag68,xyzzyaaac68)
xyzzyaaaw68=grad_eebasis1(1:3,xyzzyaaag68,xyzzyaaac68)
xyzzyaaax68=lap_eebasis1(xyzzyaaag68,xyzzyaaac68)
xyzzyaaay68=xyzzyaaam68*xyzzyaaav68
xyzzyaaaz68=xyzzyaaam68*xyzzyaaaw68+xyzzyaaao68*xyzzyaaav68
xyzzyaaba68=xyzzyaaam68*xyzzyaaax68+2*ddot(3,xyzzyaaao68,1,xyzzyaaaw68&
&,1)+xyzzyaaaq68*xyzzyaaav68
do xyzzyaaah68=1,order_en
xyzzyaaav68=enbasis1(xyzzyaaah68,xyzzyaaad68)
xyzzyaaaw68=grad_enbasis1(1:3,xyzzyaaah68,xyzzyaaad68)
xyzzyaaax68=lap_enbasis1(xyzzyaaah68,xyzzyaaad68)
xyzzyaabb68=xyzzyaaay68*xyzzyaaav68
xyzzyaabc68=xyzzyaaay68*xyzzyaaaw68+xyzzyaaaz68*xyzzyaaav68
xyzzyaabd68=xyzzyaaay68*xyzzyaaax68+2*ddot(3,xyzzyaaaz68,1,xyzzyaaaw68&
&,1)+xyzzyaaba68*xyzzyaaav68
do xyzzyaaai68=1,order_en
xyzzyaaav68=enbasis(xyzzyaaai68,xyzzyaaad68,xyzzyaaac68)
xyzzyaabe68=xyzzyaabb68*xyzzyaaav68
xyzzyaabf68=xyzzyaabc68*xyzzyaaav68
xyzzyaabg68=xyzzyaabd68*xyzzyaaav68
xyzzyaabh68=param(xyzzyaaai68,xyzzyaaah68,xyzzyaaag68,xyzzyaaae68)
jas=jas+xyzzyaabh68*xyzzyaabe68
gjas=gjas+xyzzyaabh68*xyzzyaabf68
ljas=ljas+xyzzyaabh68*xyzzyaabg68
djas(xyzzyaaai68,xyzzyaaah68,xyzzyaaag68,xyzzyaaae68)=djas(xyzzyaaai68&
&,xyzzyaaah68,xyzzyaaag68,xyzzyaaae68)+xyzzyaabe68
dgjas(1:3,xyzzyaaai68,xyzzyaaah68,xyzzyaaag68,xyzzyaaae68)=dgjas(1:3,x&
&yzzyaaai68,xyzzyaaah68,xyzzyaaag68,xyzzyaaae68)+xyzzyaabf68
dljas(xyzzyaaai68,xyzzyaaah68,xyzzyaaag68,xyzzyaaae68)=dljas(xyzzyaaai&
&68,xyzzyaaah68,xyzzyaaag68,xyzzyaaae68)+xyzzyaabg68
enddo
enddo
enddo
enddo
else
do xyzzyaaak68=1,nele(xyzzyaaab68)
xyzzyaaac68=xyzzyaaac68+1
if(xyzzyaaac68==ii)cycle
if(.not.nzeecut1(xyzzyaaac68).or..not.nzencut(xyzzyaaad68,xyzzyaaac68)&
&)cycle
xyzzyaaam68=eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaaad68)*encut(1,xyzzya&
&aad68,xyzzyaaac68)
xyzzyaaao68=(grad_eecut1(:,1,xyzzyaaac68)*encut1(1,xyzzyaaad68)+eecut1&
&(1,xyzzyaaac68)*grad_encut1(:,1,xyzzyaaad68))*encut(1,xyzzyaaad68,xyz&
&zyaaac68)
xyzzyaaaq68=(lap_eecut1(1,xyzzyaaac68)*encut1(1,xyzzyaaad68)+2*ddot(3,&
&grad_eecut1(1,1,xyzzyaaac68),1,grad_encut1(1,1,xyzzyaaad68),1)+eecut1&
&(1,xyzzyaaac68)*lap_encut1(1,xyzzyaaad68))*encut(1,xyzzyaaad68,xyzzya&
&aac68)
do xyzzyaaag68=1,order_ee
xyzzyaaav68=eebasis1(xyzzyaaag68,xyzzyaaac68)
xyzzyaaaw68=grad_eebasis1(1:3,xyzzyaaag68,xyzzyaaac68)
xyzzyaaax68=lap_eebasis1(xyzzyaaag68,xyzzyaaac68)
xyzzyaaay68=xyzzyaaam68*xyzzyaaav68
xyzzyaaaz68=xyzzyaaam68*xyzzyaaaw68+xyzzyaaao68*xyzzyaaav68
xyzzyaaba68=xyzzyaaam68*xyzzyaaax68+2*ddot(3,xyzzyaaao68,1,xyzzyaaaw68&
&,1)+xyzzyaaaq68*xyzzyaaav68
do xyzzyaaai68=1,order_en
xyzzyaaav68=enbasis(xyzzyaaai68,xyzzyaaad68,xyzzyaaac68)
xyzzyaabb68=xyzzyaaay68*xyzzyaaav68
xyzzyaabc68=xyzzyaaaz68*xyzzyaaav68
xyzzyaabd68=xyzzyaaba68*xyzzyaaav68
do xyzzyaaah68=1,order_en
xyzzyaaav68=enbasis1(xyzzyaaah68,xyzzyaaad68)
xyzzyaaaw68=grad_enbasis1(1:3,xyzzyaaah68,xyzzyaaad68)
xyzzyaaax68=lap_enbasis1(xyzzyaaah68,xyzzyaaad68)
xyzzyaabe68=xyzzyaabb68*xyzzyaaav68
xyzzyaabf68=xyzzyaabb68*xyzzyaaaw68+xyzzyaabc68*xyzzyaaav68
xyzzyaabg68=xyzzyaabb68*xyzzyaaax68+2*ddot(3,xyzzyaabc68(1),1,xyzzyaaa&
&w68,1)+xyzzyaabd68*xyzzyaaav68
xyzzyaabh68=param(xyzzyaaah68,xyzzyaaai68,xyzzyaaag68,xyzzyaaae68)
jas=jas+xyzzyaabh68*xyzzyaabe68
gjas=gjas+xyzzyaabh68*xyzzyaabf68
ljas=ljas+xyzzyaabh68*xyzzyaabg68
djas(xyzzyaaah68,xyzzyaaai68,xyzzyaaag68,xyzzyaaae68)=djas(xyzzyaaah68&
&,xyzzyaaai68,xyzzyaaag68,xyzzyaaae68)+xyzzyaabe68
dgjas(1:3,xyzzyaaah68,xyzzyaaai68,xyzzyaaag68,xyzzyaaae68)=dgjas(1:3,x&
&yzzyaaah68,xyzzyaaai68,xyzzyaaag68,xyzzyaaae68)+xyzzyaabf68
dljas(xyzzyaaah68,xyzzyaaai68,xyzzyaaag68,xyzzyaaae68)=dljas(xyzzyaaah&
&68,xyzzyaaai68,xyzzyaaag68,xyzzyaaae68)+xyzzyaabg68
enddo
enddo
enddo
enddo
endif
endif
xyzzyaaaj68=xyzzyaaaj68+nele(xyzzyaaab68)
enddo
enddo
end subroutine xyzzyaaim1
subroutine xyzzyaain1(ii,order_ee,order_en,groups_ee1,groups_en,is_u1,&
&param,nchannel,chsig,eebasis1,eecut1,nzeecut1,enbasis1,encut1,nzencut&
&1,enbasis,encut,nzencut,deebasis1,deecut1,denbasis1,dencut1,gradr_eeb&
&asis1,gradr_enbasis1,d2eebasis1,d2eecut1,d2enbasis1,d2encut1,lapr_eeb&
&asis1,lapr_enbasis1,jas,gjas,ljas,djas,dgjas,dljas)
implicit none
integer,intent(in) :: ii,order_ee,order_en,nchannel,groups_ee1(nspin),&
&groups_en(nitot,nspin)
logical,intent(in) :: is_u1
integer,intent(in) :: chsig(xyzzyaaae1,*)
real(dp),intent(in) :: eebasis1(nfn_eebasis,*),eecut1(nfn_eebasis,*),e&
&nbasis1(nfn_enbasis,*),encut1(nfn_enbasis,*),enbasis(nfn_enbasis,nito&
&t,*),encut(nfn_enbasis,nitot,*),param(order_en,order_en,order_ee,*)
real(dp),intent(in),optional :: deebasis1(nfn_eebasis,*),deecut1(nfn_e&
&ebasis,*),denbasis1(nfn_enbasis,*),dencut1(nfn_enbasis,*),gradr_eebas&
&is1(3,*),gradr_enbasis1(3,*),d2eebasis1(nfn_eebasis,*),d2eecut1(nfn_e&
&ebasis,*),d2enbasis1(nfn_enbasis,*),d2encut1(nfn_enbasis,*),lapr_eeba&
&sis1(*),lapr_enbasis1(*)
real(dp),intent(inout),optional :: jas,gjas(3),ljas,djas(order_en,orde&
&r_en,order_ee,nchannel),dgjas(3,order_en,order_en,order_ee,nchannel),&
&dljas(order_en,order_en,order_ee,nchannel)
logical,intent(in) :: nzeecut1(netot),nzencut1(nitot),nzencut(nitot,ne&
&tot)
integer xyzzyaaaa69,xyzzyaaab69,xyzzyaaac69,xyzzyaaad69,xyzzyaaae69,xy&
&zzyaaaf69(3),xyzzyaaag69,xyzzyaaah69,xyzzyaaai69,xyzzyaaaj69,xyzzyaaa&
&k69
real(dp) xyzzyaaal69,xyzzyaaam69,xyzzyaaan69,xyzzyaaao69,xyzzyaaap69,x&
&yzzyaaaq69,xyzzyaaar69,xyzzyaaas69,xyzzyaaat69,xyzzyaaau69,xyzzyaaav6&
&9,xyzzyaaaw69,xyzzyaaax69,xyzzyaaay69,xyzzyaaaz69,xyzzyaaba69,xyzzyaa&
&bb69(3),xyzzyaabc69,xyzzyaabd69,xyzzyaabe69(3),xyzzyaabf69,xyzzyaabg6&
&9,xyzzyaabh69(3),xyzzyaabi69,xyzzyaabj69,xyzzyaabk69(3),xyzzyaabl69,x&
&yzzyaabm69,xyzzyaabn69(3),xyzzyaabo69,xyzzyaabp69
logical xyzzyaabq69,xyzzyaabr69,xyzzyaabs69,xyzzyaabt69,xyzzyaabu69
xyzzyaabu69=present(dljas)
xyzzyaabq69=.true.
if(xyzzyaabu69.or.present(ljas))then
xyzzyaabr69=.true.
xyzzyaabs69=.true.
elseif(present(gjas))then
xyzzyaabr69=.true.
xyzzyaabs69=.false.
else
xyzzyaabr69=.false.
xyzzyaabs69=.false.
endif
if(xyzzyaabq69)jas=0.d0
if(xyzzyaabr69)then
gjas=0.d0
if(xyzzyaabs69)ljas=0.d0
endif
if(xyzzyaabu69)then
if(xyzzyaabq69)djas=0.d0
if(xyzzyaabr69)then
dgjas=0.d0
if(xyzzyaabs69)dljas=0.d0
endif
endif
xyzzyaaaa69=which_spin(ii)
do xyzzyaaad69=1,nitot
if(.not.nzencut1(xyzzyaaad69))cycle
xyzzyaaag69=0
do xyzzyaaab69=1,nspin
xyzzyaaaf69=(/groups_ee1(xyzzyaaab69),groups_en(xyzzyaaad69,xyzzyaaaa6&
&9),groups_en(xyzzyaaad69,xyzzyaaab69)/)
xyzzyaabt69=xyzzyaaaf69(2)>xyzzyaaaf69(3)
if(xyzzyaabt69)call swap1(xyzzyaaaf69(2),xyzzyaaaf69(3))
call match_signature(3,xyzzyaaae1,nchannel,xyzzyaaaf69,chsig,xyzzyaaae&
&69)
if(xyzzyaaae69<1)then
xyzzyaaag69=xyzzyaaag69+nele(xyzzyaaab69)
cycle
endif
xyzzyaaac69=xyzzyaaag69
if(.not.xyzzyaabr69)then
if(.not.xyzzyaabt69)then
do xyzzyaaah69=1,nele(xyzzyaaab69)
xyzzyaaac69=xyzzyaaac69+1
if(xyzzyaaac69==ii.or..not.nzeecut1(xyzzyaaac69).or..not.nzencut(xyzzy&
&aaad69,xyzzyaaac69))cycle
xyzzyaaam69=eecut1(1,xyzzyaaac69)*encut1(1,xyzzyaaad69)*encut(1,xyzzya&
&aad69,xyzzyaaac69)
call xyzzyaaio1(is_u1,order_ee,order_en,eebasis1(1,xyzzyaaac69),enbasi&
&s1(1,xyzzyaaad69),enbasis(1,xyzzyaaad69,xyzzyaaac69),param(1,1,1,xyzz&
&yaaae69),xyzzyaaal69)
jas=jas+xyzzyaaal69*xyzzyaaam69
enddo
else
do xyzzyaaah69=1,nele(xyzzyaaab69)
xyzzyaaac69=xyzzyaaac69+1
if(xyzzyaaac69==ii.or..not.nzeecut1(xyzzyaaac69).or..not.nzencut(xyzzy&
&aaad69,xyzzyaaac69))cycle
xyzzyaaam69=eecut1(1,xyzzyaaac69)*encut1(1,xyzzyaaad69)*encut(1,xyzzya&
&aad69,xyzzyaaac69)
call xyzzyaaio1(is_u1,order_ee,order_en,eebasis1(1,xyzzyaaac69),enbasi&
&s(1,xyzzyaaad69,xyzzyaaac69),enbasis1(1,xyzzyaaad69),param(1,1,1,xyzz&
&yaaae69),xyzzyaaal69)
jas=jas+xyzzyaaal69*xyzzyaaam69
enddo
endif
elseif(.not.xyzzyaabu69)then
do xyzzyaaah69=1,nele(xyzzyaaab69)
xyzzyaaac69=xyzzyaaac69+1
if(xyzzyaaac69==ii.or..not.nzeecut1(xyzzyaaac69).or..not.nzencut(xyzzy&
&aaad69,xyzzyaaac69))cycle
xyzzyaaaz69=encut1(1,xyzzyaaad69)*encut(1,xyzzyaaad69,xyzzyaaac69)
xyzzyaaba69=eecut1(1,xyzzyaaac69)*encut(1,xyzzyaaad69,xyzzyaaac69)
xyzzyaaam69=eecut1(1,xyzzyaaac69)*xyzzyaaaz69
xyzzyaaap69=deecut1(1,xyzzyaaac69)*xyzzyaaaz69
xyzzyaaaq69=xyzzyaaba69*dencut1(1,xyzzyaaad69)
if(xyzzyaabs69)then
xyzzyaaaw69=d2eecut1(1,xyzzyaaac69)*xyzzyaaaz69
xyzzyaaax69=xyzzyaaba69*d2encut1(1,xyzzyaaad69)
xyzzyaaay69=deecut1(1,xyzzyaaac69)*dencut1(1,xyzzyaaad69)*encut(1,xyzz&
&yaaad69,xyzzyaaac69)
endif
call xyzzyaaip1(is_u1,xyzzyaabs69,xyzzyaabt69,order_ee,order_en,eebasi&
&s1(1,xyzzyaaac69),enbasis1(1,xyzzyaaad69),enbasis(1,xyzzyaaad69,xyzzy&
&aaac69),deebasis1(1,xyzzyaaac69),denbasis1(1,xyzzyaaad69),d2eebasis1(&
&1,xyzzyaaac69),d2enbasis1(1,xyzzyaaad69),param(1,1,1,xyzzyaaae69),xyz&
&zyaaal69,xyzzyaaan69,xyzzyaaao69,xyzzyaaat69,xyzzyaaau69,xyzzyaaav69)
jas=jas+xyzzyaaal69*xyzzyaaam69
xyzzyaaar69=xyzzyaaan69*xyzzyaaam69+xyzzyaaal69*xyzzyaaap69
xyzzyaaas69=xyzzyaaao69*xyzzyaaam69+xyzzyaaal69*xyzzyaaaq69
gjas=gjas+xyzzyaaar69*gradr_eebasis1(:,xyzzyaaac69)+xyzzyaaas69*gradr_&
&enbasis1(:,xyzzyaaad69)
if(xyzzyaabs69)ljas=ljas+xyzzyaaat69*xyzzyaaam69+2*xyzzyaaan69*xyzzyaa&
&ap69+xyzzyaaal69*xyzzyaaaw69+xyzzyaaau69*xyzzyaaam69+2*xyzzyaaao69*xy&
&zzyaaaq69+xyzzyaaal69*xyzzyaaax69+2*(xyzzyaaav69*xyzzyaaam69+xyzzyaaa&
&n69*xyzzyaaaq69+xyzzyaaao69*xyzzyaaap69+xyzzyaaal69*xyzzyaaay69)*ddot&
&(3,gradr_eebasis1(1,xyzzyaaac69),1,gradr_enbasis1(1,xyzzyaaad69),1)+x&
&yzzyaaar69*lapr_eebasis1(xyzzyaaac69)+xyzzyaaas69*lapr_enbasis1(xyzzy&
&aaad69)
enddo
else
if(.not.xyzzyaabt69)then
do xyzzyaaah69=1,nele(xyzzyaaab69)
xyzzyaaac69=xyzzyaaac69+1
if(xyzzyaaac69==ii)cycle
if(.not.nzeecut1(xyzzyaaac69).or..not.nzencut(xyzzyaaad69,xyzzyaaac69)&
&)cycle
xyzzyaaaz69=encut1(1,xyzzyaaad69)*encut(1,xyzzyaaad69,xyzzyaaac69)
xyzzyaaba69=eecut1(1,xyzzyaaac69)*encut(1,xyzzyaaad69,xyzzyaaac69)
xyzzyaaam69=eecut1(1,xyzzyaaac69)*xyzzyaaaz69
xyzzyaaap69=deecut1(1,xyzzyaaac69)*xyzzyaaaz69
xyzzyaaaq69=xyzzyaaba69*dencut1(1,xyzzyaaad69)
xyzzyaabb69=xyzzyaaap69*gradr_eebasis1(:,xyzzyaaac69)+xyzzyaaaq69*grad&
&r_enbasis1(:,xyzzyaaad69)
xyzzyaabc69=d2eecut1(1,xyzzyaaac69)*xyzzyaaaz69+xyzzyaaba69*d2encut1(1&
&,xyzzyaaad69)+2*deecut1(1,xyzzyaaac69)*dencut1(1,xyzzyaaad69)*encut(1&
&,xyzzyaaad69,xyzzyaaac69)*ddot(3,gradr_eebasis1(1,xyzzyaaac69),1,grad&
&r_enbasis1(1,xyzzyaaad69),1)+xyzzyaaap69*lapr_eebasis1(xyzzyaaac69)+x&
&yzzyaaaq69*lapr_enbasis1(xyzzyaaad69)
do xyzzyaaai69=1,order_ee
xyzzyaabd69=eebasis1(xyzzyaaai69,xyzzyaaac69)
xyzzyaabe69=deebasis1(xyzzyaaai69,xyzzyaaac69)*gradr_eebasis1(:,xyzzya&
&aac69)
xyzzyaabf69=d2eebasis1(xyzzyaaai69,xyzzyaaac69)+deebasis1(xyzzyaaai69,&
&xyzzyaaac69)*lapr_eebasis1(xyzzyaaac69)
xyzzyaabg69=xyzzyaaam69*xyzzyaabd69
xyzzyaabh69=xyzzyaaam69*xyzzyaabe69+xyzzyaabb69*xyzzyaabd69
xyzzyaabi69=xyzzyaaam69*xyzzyaabf69+2*ddot(3,xyzzyaabb69,1,xyzzyaabe69&
&,1)+xyzzyaabc69*xyzzyaabd69
do xyzzyaaaj69=1,order_en
xyzzyaabd69=enbasis1(xyzzyaaaj69,xyzzyaaad69)
xyzzyaabe69=denbasis1(xyzzyaaaj69,xyzzyaaad69)*gradr_enbasis1(:,xyzzya&
&aad69)
xyzzyaabf69=d2enbasis1(xyzzyaaaj69,xyzzyaaad69)+denbasis1(xyzzyaaaj69,&
&xyzzyaaad69)*lapr_enbasis1(xyzzyaaad69)
xyzzyaabj69=xyzzyaabg69*xyzzyaabd69
xyzzyaabk69=xyzzyaabg69*xyzzyaabe69+xyzzyaabh69*xyzzyaabd69
xyzzyaabl69=xyzzyaabg69*xyzzyaabf69+2*ddot(3,xyzzyaabh69,1,xyzzyaabe69&
&,1)+xyzzyaabi69*xyzzyaabd69
do xyzzyaaak69=1,order_en
xyzzyaabd69=enbasis(xyzzyaaak69,xyzzyaaad69,xyzzyaaac69)
xyzzyaabm69=xyzzyaabj69*xyzzyaabd69
xyzzyaabn69=xyzzyaabk69*xyzzyaabd69
xyzzyaabo69=xyzzyaabl69*xyzzyaabd69
xyzzyaabp69=param(xyzzyaaak69,xyzzyaaaj69,xyzzyaaai69,xyzzyaaae69)
jas=jas+xyzzyaabp69*xyzzyaabm69
gjas=gjas+xyzzyaabp69*xyzzyaabn69
ljas=ljas+xyzzyaabp69*xyzzyaabo69
djas(xyzzyaaak69,xyzzyaaaj69,xyzzyaaai69,xyzzyaaae69)=djas(xyzzyaaak69&
&,xyzzyaaaj69,xyzzyaaai69,xyzzyaaae69)+xyzzyaabm69
dgjas(1:3,xyzzyaaak69,xyzzyaaaj69,xyzzyaaai69,xyzzyaaae69)=dgjas(1:3,x&
&yzzyaaak69,xyzzyaaaj69,xyzzyaaai69,xyzzyaaae69)+xyzzyaabn69
dljas(xyzzyaaak69,xyzzyaaaj69,xyzzyaaai69,xyzzyaaae69)=dljas(xyzzyaaak&
&69,xyzzyaaaj69,xyzzyaaai69,xyzzyaaae69)+xyzzyaabo69
enddo
enddo
enddo
enddo
else
do xyzzyaaah69=1,nele(xyzzyaaab69)
xyzzyaaac69=xyzzyaaac69+1
if(xyzzyaaac69==ii)cycle
if(.not.nzeecut1(xyzzyaaac69).or..not.nzencut(xyzzyaaad69,xyzzyaaac69)&
&)cycle
xyzzyaaaz69=encut1(1,xyzzyaaad69)*encut(1,xyzzyaaad69,xyzzyaaac69)
xyzzyaaba69=eecut1(1,xyzzyaaac69)*encut(1,xyzzyaaad69,xyzzyaaac69)
xyzzyaaam69=eecut1(1,xyzzyaaac69)*xyzzyaaaz69
xyzzyaaap69=deecut1(1,xyzzyaaac69)*xyzzyaaaz69
xyzzyaaaq69=xyzzyaaba69*dencut1(1,xyzzyaaad69)
xyzzyaabb69=xyzzyaaap69*gradr_eebasis1(:,xyzzyaaac69)+xyzzyaaaq69*grad&
&r_enbasis1(:,xyzzyaaad69)
xyzzyaabc69=d2eecut1(1,xyzzyaaac69)*xyzzyaaaz69+xyzzyaaba69*d2encut1(1&
&,xyzzyaaad69)+2*deecut1(1,xyzzyaaac69)*dencut1(1,xyzzyaaad69)*encut(1&
&,xyzzyaaad69,xyzzyaaac69)*ddot(3,gradr_eebasis1(1,xyzzyaaac69),1,grad&
&r_enbasis1(1,xyzzyaaad69),1)+xyzzyaaap69*lapr_eebasis1(xyzzyaaac69)+x&
&yzzyaaaq69*lapr_enbasis1(xyzzyaaad69)
do xyzzyaaai69=1,order_ee
xyzzyaabd69=eebasis1(xyzzyaaai69,xyzzyaaac69)
xyzzyaabe69=deebasis1(xyzzyaaai69,xyzzyaaac69)*gradr_eebasis1(:,xyzzya&
&aac69)
xyzzyaabf69=d2eebasis1(xyzzyaaai69,xyzzyaaac69)+deebasis1(xyzzyaaai69,&
&xyzzyaaac69)*lapr_eebasis1(xyzzyaaac69)
xyzzyaabg69=xyzzyaaam69*xyzzyaabd69
xyzzyaabh69=xyzzyaaam69*xyzzyaabe69+xyzzyaabb69*xyzzyaabd69
xyzzyaabi69=xyzzyaaam69*xyzzyaabf69+2*ddot(3,xyzzyaabb69,1,xyzzyaabe69&
&,1)+xyzzyaabc69*xyzzyaabd69
do xyzzyaaak69=1,order_en
xyzzyaabd69=enbasis(xyzzyaaak69,xyzzyaaad69,xyzzyaaac69)
xyzzyaabj69=xyzzyaabg69*xyzzyaabd69
xyzzyaabk69=xyzzyaabh69*xyzzyaabd69
xyzzyaabl69=xyzzyaabi69*xyzzyaabd69
do xyzzyaaaj69=1,order_en
xyzzyaabd69=enbasis1(xyzzyaaaj69,xyzzyaaad69)
xyzzyaabe69=denbasis1(xyzzyaaaj69,xyzzyaaad69)*gradr_enbasis1(:,xyzzya&
&aad69)
xyzzyaabf69=d2enbasis1(xyzzyaaaj69,xyzzyaaad69)+denbasis1(xyzzyaaaj69,&
&xyzzyaaad69)*lapr_enbasis1(xyzzyaaad69)
xyzzyaabm69=xyzzyaabj69*xyzzyaabd69
xyzzyaabn69=xyzzyaabj69*xyzzyaabe69+xyzzyaabk69*xyzzyaabd69
xyzzyaabo69=xyzzyaabj69*xyzzyaabf69+2*ddot(3,xyzzyaabk69(1),1,xyzzyaab&
&e69,1)+xyzzyaabl69*xyzzyaabd69
xyzzyaabp69=param(xyzzyaaaj69,xyzzyaaak69,xyzzyaaai69,xyzzyaaae69)
jas=jas+xyzzyaabp69*xyzzyaabm69
gjas=gjas+xyzzyaabp69*xyzzyaabn69
ljas=ljas+xyzzyaabp69*xyzzyaabo69
djas(xyzzyaaaj69,xyzzyaaak69,xyzzyaaai69,xyzzyaaae69)=djas(xyzzyaaaj69&
&,xyzzyaaak69,xyzzyaaai69,xyzzyaaae69)+xyzzyaabm69
dgjas(1:3,xyzzyaaaj69,xyzzyaaak69,xyzzyaaai69,xyzzyaaae69)=dgjas(1:3,x&
&yzzyaaaj69,xyzzyaaak69,xyzzyaaai69,xyzzyaaae69)+xyzzyaabn69
dljas(xyzzyaaaj69,xyzzyaaak69,xyzzyaaai69,xyzzyaaae69)=dljas(xyzzyaaaj&
&69,xyzzyaaak69,xyzzyaaai69,xyzzyaaae69)+xyzzyaabo69
enddo
enddo
enddo
enddo
endif
endif
xyzzyaaag69=xyzzyaaag69+nele(xyzzyaaab69)
enddo
enddo
end subroutine xyzzyaain1
subroutine xyzzyaaio1(is_u1,order_ee,order_en,pow_rij,pow_rii,pow_rji,&
&param_ichannel,f_core)
implicit none
logical,intent(in) :: is_u1
integer,intent(in) :: order_ee,order_en
real(dp),intent(in) :: pow_rij(*),pow_rii(*),pow_rji(*),param_ichannel&
&(order_en,order_en,order_ee)
real(dp),intent(out) :: f_core
integer xyzzyaaaa70,xyzzyaaab70,xyzzyaaac70
real(dp) xyzzyaaad70,xyzzyaaae70
if(is_u1.and.order_ee==3.and.order_en==3)then
f_core=param_ichannel(1,1,1)+param_ichannel(2,1,1)*pow_rji(2)+param_ic&
&hannel(3,1,1)*pow_rji(3)+(param_ichannel(1,2,1)+ param_ichannel(2,2,1&
&)*pow_rji(2)+ param_ichannel(3,2,1)*pow_rji(3))*pow_rii(2)+(param_ich&
&annel(1,3,1)+ param_ichannel(2,3,1)*pow_rji(2)+ param_ichannel(3,3,1)&
&*pow_rji(3))*pow_rii(3)+(param_ichannel(1,1,2)+ param_ichannel(2,1,2)&
&*pow_rji(2)+ param_ichannel(3,1,2)*pow_rji(3)+ (param_ichannel(1,2,2)&
&+  param_ichannel(2,2,2)*pow_rji(2)+  param_ichannel(3,2,2)*pow_rji(3&
&))*pow_rii(2)+ (param_ichannel(1,3,2)+  param_ichannel(2,3,2)*pow_rji&
&(2)+  param_ichannel(3,3,2)*pow_rji(3))*pow_rii(3))*pow_rij(2)+(param&
&_ichannel(1,1,3)+ param_ichannel(2,1,3)*pow_rji(2)+ param_ichannel(3,&
&1,3)*pow_rji(3)+ (param_ichannel(1,2,3)+  param_ichannel(2,2,3)*pow_r&
&ji(2)+  param_ichannel(3,2,3)*pow_rji(3))*pow_rii(2)+ (param_ichannel&
&(1,3,3)+  param_ichannel(2,3,3)*pow_rji(2)+  param_ichannel(3,3,3)*po&
&w_rji(3))*pow_rii(3))*pow_rij(3)
elseif(is_u1)then
f_core=param_ichannel(1,1,1)
do xyzzyaaaa70=2,order_en
f_core=f_core+param_ichannel(xyzzyaaaa70,1,1)*pow_rji(xyzzyaaaa70)
enddo
do xyzzyaaab70=2,order_en
xyzzyaaad70=param_ichannel(1,xyzzyaaab70,1)
do xyzzyaaaa70=2,order_en
xyzzyaaad70=xyzzyaaad70+param_ichannel(xyzzyaaaa70,xyzzyaaab70,1)*pow_&
&rji(xyzzyaaaa70)
enddo
f_core=f_core+xyzzyaaad70*pow_rii(xyzzyaaab70)
enddo
do xyzzyaaac70=2,order_ee
xyzzyaaad70=param_ichannel(1,1,xyzzyaaac70)
do xyzzyaaaa70=2,order_en
xyzzyaaad70=xyzzyaaad70+param_ichannel(xyzzyaaaa70,1,xyzzyaaac70)*pow_&
&rji(xyzzyaaaa70)
enddo
do xyzzyaaab70=2,order_en
xyzzyaaae70=param_ichannel(1,xyzzyaaab70,xyzzyaaac70)
do xyzzyaaaa70=2,order_en
xyzzyaaae70=xyzzyaaae70+param_ichannel(xyzzyaaaa70,xyzzyaaab70,xyzzyaa&
&ac70)*pow_rji(xyzzyaaaa70)
enddo
xyzzyaaad70=xyzzyaaad70+xyzzyaaae70*pow_rii(xyzzyaaab70)
enddo
f_core=f_core+xyzzyaaad70*pow_rij(xyzzyaaac70)
enddo
else
f_core=0.d0
do xyzzyaaac70=1,order_ee
xyzzyaaad70=0.d0
do xyzzyaaab70=1,order_en
xyzzyaaae70=0.d0
do xyzzyaaaa70=1,order_en
xyzzyaaae70=xyzzyaaae70+param_ichannel(xyzzyaaaa70,xyzzyaaab70,xyzzyaa&
&ac70)*pow_rji(xyzzyaaaa70)
enddo
xyzzyaaad70=xyzzyaaad70+xyzzyaaae70*pow_rii(xyzzyaaab70)
enddo
f_core=f_core+xyzzyaaad70*pow_rij(xyzzyaaac70)
enddo
endif
end subroutine xyzzyaaio1
subroutine xyzzyaaip1(is_u1,sd,perm,order_ee,order_en,pow_rij,pow_rii,&
&pow_rji,dpow_rij,dpow_rii,d2pow_rij,d2pow_rii,param_ichannel,f_core,d&
&f_core_rij,df_core_rii,d2f_core_rij,d2f_core_rii,d2f_core_rij_rii)
implicit none
logical,intent(in) :: is_u1,sd,perm
integer,intent(in) :: order_ee,order_en
real(dp),intent(in) :: pow_rij(*),pow_rii(*),pow_rji(*),dpow_rij(*),dp&
&ow_rii(*),d2pow_rij(*),d2pow_rii(*),param_ichannel(order_en,order_en,&
&order_ee)
real(dp),intent(inout) :: f_core,df_core_rij,df_core_rii,d2f_core_rij,&
&d2f_core_rii,d2f_core_rij_rii
integer xyzzyaaaa71,xyzzyaaab71,xyzzyaaac71
real(dp) xyzzyaaad71,xyzzyaaae71,xyzzyaaaf71,xyzzyaaag71,xyzzyaaah71,x&
&yzzyaaai71,xyzzyaaaj71,xyzzyaaak71,xyzzyaaal71,xyzzyaaam71,xyzzyaaan7&
&1,xyzzyaaao71,xyzzyaaap71,xyzzyaaaq71,xyzzyaaar71,xyzzyaaas71
if(is_u1.and.order_ee==3.and.order_en==3)then
if(.not.perm)then
xyzzyaaad71=param_ichannel(1,1,1)+param_ichannel(2,1,1)*pow_rji(2)+par&
&am_ichannel(3,1,1)*pow_rji(3)
xyzzyaaae71=param_ichannel(1,2,1)+param_ichannel(2,2,1)*pow_rji(2)+par&
&am_ichannel(3,2,1)*pow_rji(3)
xyzzyaaaf71=param_ichannel(1,3,1)+param_ichannel(2,3,1)*pow_rji(2)+par&
&am_ichannel(3,3,1)*pow_rji(3)
xyzzyaaag71=param_ichannel(1,1,2)+param_ichannel(2,1,2)*pow_rji(2)+par&
&am_ichannel(3,1,2)*pow_rji(3)
xyzzyaaah71=param_ichannel(1,2,2)+param_ichannel(2,2,2)*pow_rji(2)+par&
&am_ichannel(3,2,2)*pow_rji(3)
xyzzyaaai71=param_ichannel(1,3,2)+param_ichannel(2,3,2)*pow_rji(2)+par&
&am_ichannel(3,3,2)*pow_rji(3)
xyzzyaaaj71=param_ichannel(1,1,3)+param_ichannel(2,1,3)*pow_rji(2)+par&
&am_ichannel(3,1,3)*pow_rji(3)
xyzzyaaak71=param_ichannel(1,2,3)+param_ichannel(2,2,3)*pow_rji(2)+par&
&am_ichannel(3,2,3)*pow_rji(3)
xyzzyaaal71=param_ichannel(1,3,3)+param_ichannel(2,3,3)*pow_rji(2)+par&
&am_ichannel(3,3,3)*pow_rji(3)
else
xyzzyaaad71=param_ichannel(1,1,1)+param_ichannel(1,2,1)*pow_rji(2)+par&
&am_ichannel(1,3,1)*pow_rji(3)
xyzzyaaae71=param_ichannel(2,1,1)+param_ichannel(2,2,1)*pow_rji(2)+par&
&am_ichannel(2,3,1)*pow_rji(3)
xyzzyaaaf71=param_ichannel(3,1,1)+param_ichannel(3,2,1)*pow_rji(2)+par&
&am_ichannel(3,3,1)*pow_rji(3)
xyzzyaaag71=param_ichannel(1,1,2)+param_ichannel(1,2,2)*pow_rji(2)+par&
&am_ichannel(1,3,2)*pow_rji(3)
xyzzyaaah71=param_ichannel(2,1,2)+param_ichannel(2,2,2)*pow_rji(2)+par&
&am_ichannel(2,3,2)*pow_rji(3)
xyzzyaaai71=param_ichannel(3,1,2)+param_ichannel(3,2,2)*pow_rji(2)+par&
&am_ichannel(3,3,2)*pow_rji(3)
xyzzyaaaj71=param_ichannel(1,1,3)+param_ichannel(1,2,3)*pow_rji(2)+par&
&am_ichannel(1,3,3)*pow_rji(3)
xyzzyaaak71=param_ichannel(2,1,3)+param_ichannel(2,2,3)*pow_rji(2)+par&
&am_ichannel(2,3,3)*pow_rji(3)
xyzzyaaal71=param_ichannel(3,1,3)+param_ichannel(3,2,3)*pow_rji(2)+par&
&am_ichannel(3,3,3)*pow_rji(3)
endif
xyzzyaaam71=xyzzyaaag71+xyzzyaaah71*pow_rii(2)+xyzzyaaai71*pow_rii(3)
xyzzyaaan71=xyzzyaaaj71+xyzzyaaak71*pow_rii(2)+xyzzyaaal71*pow_rii(3)
xyzzyaaao71=xyzzyaaah71*dpow_rii(2)+xyzzyaaai71*dpow_rii(3)
xyzzyaaap71=xyzzyaaak71*dpow_rii(2)+xyzzyaaal71*dpow_rii(3)
f_core=xyzzyaaad71+xyzzyaaae71*pow_rii(2)+xyzzyaaaf71*pow_rii(3)+xyzzy&
&aaam71*pow_rij(2)+xyzzyaaan71*pow_rij(3)
df_core_rij=xyzzyaaam71*dpow_rij(2)+xyzzyaaan71*dpow_rij(3)
df_core_rii=xyzzyaaae71*dpow_rii(2)+xyzzyaaaf71*dpow_rii(3)+xyzzyaaao7&
&1*pow_rij(2)+xyzzyaaap71*pow_rij(3)
if(sd)then
d2f_core_rij=xyzzyaaam71*d2pow_rij(2)+xyzzyaaan71*d2pow_rij(3)
d2f_core_rii=xyzzyaaae71*d2pow_rii(2)+xyzzyaaaf71*d2pow_rii(3)+(xyzzya&
&aah71*d2pow_rii(2)+xyzzyaaai71*d2pow_rii(3))*pow_rij(2)+(xyzzyaaak71*&
&d2pow_rii(2)+xyzzyaaal71*d2pow_rii(3))*pow_rij(3)
d2f_core_rij_rii=xyzzyaaao71*dpow_rij(2)+xyzzyaaap71*dpow_rij(3)
endif
elseif(is_u1)then
f_core=param_ichannel(1,1,1)
df_core_rij=0.d0
df_core_rii=0.d0
if(sd)then
d2f_core_rij=0.d0
d2f_core_rii=0.d0
d2f_core_rij_rii=0.d0
endif
if(.not.perm)then
do xyzzyaaaa71=2,order_en
f_core=f_core+param_ichannel(xyzzyaaaa71,1,1)*pow_rji(xyzzyaaaa71)
enddo
do xyzzyaaab71=2,order_en
xyzzyaaad71=param_ichannel(1,xyzzyaaab71,1)
do xyzzyaaaa71=2,order_en
xyzzyaaad71=xyzzyaaad71+param_ichannel(xyzzyaaaa71,xyzzyaaab71,1)*pow_&
&rji(xyzzyaaaa71)
enddo
f_core=f_core+xyzzyaaad71*pow_rii(xyzzyaaab71)
df_core_rii=df_core_rii+xyzzyaaad71*dpow_rii(xyzzyaaab71)
if(sd)d2f_core_rii=d2f_core_rii+xyzzyaaad71*d2pow_rii(xyzzyaaab71)
enddo
do xyzzyaaac71=2,order_ee
xyzzyaaaq71=param_ichannel(1,1,xyzzyaaac71)
xyzzyaaar71=0.d0
if(sd)xyzzyaaas71=0.d0
do xyzzyaaaa71=2,order_en
xyzzyaaaq71=xyzzyaaaq71+param_ichannel(xyzzyaaaa71,1,xyzzyaaac71)*pow_&
&rji(xyzzyaaaa71)
enddo
do xyzzyaaab71=2,order_en
xyzzyaaad71=param_ichannel(1,xyzzyaaab71,xyzzyaaac71)
do xyzzyaaaa71=2,order_en
xyzzyaaad71=xyzzyaaad71+param_ichannel(xyzzyaaaa71,xyzzyaaab71,xyzzyaa&
&ac71)*pow_rji(xyzzyaaaa71)
enddo
xyzzyaaaq71=xyzzyaaaq71+xyzzyaaad71*pow_rii(xyzzyaaab71)
xyzzyaaar71=xyzzyaaar71+xyzzyaaad71*dpow_rii(xyzzyaaab71)
if(sd)xyzzyaaas71=xyzzyaaas71+xyzzyaaad71*d2pow_rii(xyzzyaaab71)
enddo
f_core=f_core+xyzzyaaaq71*pow_rij(xyzzyaaac71)
df_core_rij=df_core_rij+xyzzyaaaq71*dpow_rij(xyzzyaaac71)
df_core_rii=df_core_rii+xyzzyaaar71*pow_rij(xyzzyaaac71)
if(sd)then
d2f_core_rij=d2f_core_rij+xyzzyaaaq71*d2pow_rij(xyzzyaaac71)
d2f_core_rii=d2f_core_rii+xyzzyaaas71*pow_rij(xyzzyaaac71)
d2f_core_rij_rii=d2f_core_rij_rii+xyzzyaaar71*dpow_rij(xyzzyaaac71)
endif
enddo
else
do xyzzyaaaa71=2,order_en
f_core=f_core+param_ichannel(1,xyzzyaaaa71,1)*pow_rji(xyzzyaaaa71)
enddo
do xyzzyaaab71=2,order_en
xyzzyaaad71=param_ichannel(xyzzyaaab71,1,1)
do xyzzyaaaa71=2,order_en
xyzzyaaad71=xyzzyaaad71+param_ichannel(xyzzyaaab71,xyzzyaaaa71,1)*pow_&
&rji(xyzzyaaaa71)
enddo
f_core=f_core+xyzzyaaad71*pow_rii(xyzzyaaab71)
df_core_rii=df_core_rii+xyzzyaaad71*dpow_rii(xyzzyaaab71)
if(sd)d2f_core_rii=d2f_core_rii+xyzzyaaad71*d2pow_rii(xyzzyaaab71)
enddo
do xyzzyaaac71=2,order_ee
xyzzyaaaq71=param_ichannel(1,1,xyzzyaaac71)
xyzzyaaar71=0.d0
if(sd)xyzzyaaas71=0.d0
do xyzzyaaaa71=2,order_en
xyzzyaaaq71=xyzzyaaaq71+param_ichannel(1,xyzzyaaaa71,xyzzyaaac71)*pow_&
&rji(xyzzyaaaa71)
enddo
do xyzzyaaab71=2,order_en
xyzzyaaad71=param_ichannel(xyzzyaaab71,1,xyzzyaaac71)
do xyzzyaaaa71=2,order_en
xyzzyaaad71=xyzzyaaad71+param_ichannel(xyzzyaaab71,xyzzyaaaa71,xyzzyaa&
&ac71)*pow_rji(xyzzyaaaa71)
enddo
xyzzyaaaq71=xyzzyaaaq71+xyzzyaaad71*pow_rii(xyzzyaaab71)
xyzzyaaar71=xyzzyaaar71+xyzzyaaad71*dpow_rii(xyzzyaaab71)
if(sd)xyzzyaaas71=xyzzyaaas71+xyzzyaaad71*d2pow_rii(xyzzyaaab71)
enddo
f_core=f_core+xyzzyaaaq71*pow_rij(xyzzyaaac71)
df_core_rij=df_core_rij+xyzzyaaaq71*dpow_rij(xyzzyaaac71)
df_core_rii=df_core_rii+xyzzyaaar71*pow_rij(xyzzyaaac71)
if(sd)then
d2f_core_rij=d2f_core_rij+xyzzyaaaq71*d2pow_rij(xyzzyaaac71)
d2f_core_rii=d2f_core_rii+xyzzyaaas71*pow_rij(xyzzyaaac71)
d2f_core_rij_rii=d2f_core_rij_rii+xyzzyaaar71*dpow_rij(xyzzyaaac71)
endif
enddo
endif
else
if(.not.perm)then
do xyzzyaaac71=1,order_ee
xyzzyaaaq71=0.d0
xyzzyaaar71=0.d0
if(sd)xyzzyaaas71=0.d0
do xyzzyaaab71=1,order_en
xyzzyaaad71=0.d0
do xyzzyaaaa71=1,order_en
xyzzyaaad71=xyzzyaaad71+param_ichannel(xyzzyaaaa71,xyzzyaaab71,xyzzyaa&
&ac71)*pow_rji(xyzzyaaaa71)
enddo
xyzzyaaaq71=xyzzyaaaq71+xyzzyaaad71*pow_rii(xyzzyaaab71)
xyzzyaaar71=xyzzyaaar71+xyzzyaaad71*dpow_rii(xyzzyaaab71)
if(sd)xyzzyaaas71=xyzzyaaas71+xyzzyaaad71*d2pow_rii(xyzzyaaab71)
enddo
f_core=f_core+xyzzyaaaq71*pow_rij(xyzzyaaac71)
df_core_rij=df_core_rij+xyzzyaaaq71*dpow_rij(xyzzyaaac71)
df_core_rii=df_core_rii+xyzzyaaar71*pow_rij(xyzzyaaac71)
if(sd)then
d2f_core_rij=d2f_core_rij+xyzzyaaaq71*d2pow_rij(xyzzyaaac71)
d2f_core_rii=d2f_core_rii+xyzzyaaas71*pow_rij(xyzzyaaac71)
d2f_core_rij_rii=d2f_core_rij_rii+xyzzyaaar71*dpow_rij(xyzzyaaac71)
endif
enddo
else
do xyzzyaaac71=1,order_ee
xyzzyaaaq71=0.d0
xyzzyaaar71=0.d0
if(sd)xyzzyaaas71=0.d0
do xyzzyaaab71=1,order_en
xyzzyaaad71=0.d0
do xyzzyaaaa71=1,order_en
xyzzyaaad71=xyzzyaaad71+param_ichannel(xyzzyaaab71,xyzzyaaaa71,xyzzyaa&
&ac71)*pow_rji(xyzzyaaaa71)
enddo
xyzzyaaaq71=xyzzyaaaq71+xyzzyaaad71*pow_rii(xyzzyaaab71)
xyzzyaaar71=xyzzyaaar71+xyzzyaaad71*dpow_rii(xyzzyaaab71)
if(sd)xyzzyaaas71=xyzzyaaas71+xyzzyaaad71*d2pow_rii(xyzzyaaab71)
enddo
f_core=f_core+xyzzyaaaq71*pow_rij(xyzzyaaac71)
df_core_rij=df_core_rij+xyzzyaaaq71*dpow_rij(xyzzyaaac71)
df_core_rii=df_core_rii+xyzzyaaar71*pow_rij(xyzzyaaac71)
if(sd)then
d2f_core_rij=d2f_core_rij+xyzzyaaaq71*d2pow_rij(xyzzyaaac71)
d2f_core_rii=d2f_core_rii+xyzzyaaas71*pow_rij(xyzzyaaac71)
d2f_core_rij_rii=d2f_core_rij_rii+xyzzyaaar71*dpow_rij(xyzzyaaac71)
endif
enddo
endif
endif
end subroutine xyzzyaaip1
subroutine xyzzyaaiq1(rank_e,rank_n,size_ee,size_en,size_sig,index_lis&
&t,nparam_remap,index_remap,leftmost_change_ee,mask_leftmost_change_ee&
&,leftmost_change_en,mask_leftmost_change_en,uncluster_ee,uncluster_en&
&,groups_ee,groups_en,param,nchannel,chsig,iparam0,which_ee_pair,which&
&_en_pair,eebasis,eecut,nzeecut,enbasis,encut,nzencut,jas)
implicit none
integer,intent(in) :: rank_e,rank_n,size_ee,size_en,size_sig,groups_ee&
&(nspin,nspin),groups_en(nitot,nspin),nchannel
integer,intent(in) :: index_list(size_sig,*),nparam_remap(nchannel),in&
&dex_remap(*),leftmost_change_ee(*),mask_leftmost_change_ee(*),leftmos&
&t_change_en(*),mask_leftmost_change_en(*)
integer,intent(in) :: chsig(xyzzyaaae1,*),iparam0(*),which_ee_pair(2,*&
&),which_en_pair(2,*)
real(dp),intent(in) :: eebasis(nfn_eebasis,netot,*),eecut(nfn_eebasis,&
&netot,*),enbasis(nfn_enbasis,nitot,*),encut(nfn_enbasis,nitot,*),para&
&m(*)
real(dp),intent(inout) :: jas
logical,intent(in) :: nzeecut(netot,netot),nzencut(nitot,netot),unclus&
&ter_ee,uncluster_en
integer xyzzyaaaa72,xyzzyaaab72(rank_e),xyzzyaaac72(rank_e),xyzzyaaad7&
&2(rank_n),xyzzyaaae72(size_ee),xyzzyaaaf72(size_en),xyzzyaaag72,xyzzy&
&aaah72(size_sig),xyzzyaaai72,xyzzyaaaj72,xyzzyaaak72,xyzzyaaal72,xyzz&
&yaaam72(netot,rank_n),xyzzyaaan72(nitot,rank_n),xyzzyaaao72(netot,ran&
&k_e),xyzzyaaap72(rank_n),xyzzyaaaq72(rank_n),xyzzyaaar72(rank_e),xyzz&
&yaaas72(rank_n),xyzzyaaat72(rank_e)
real(dp) xyzzyaaau72,xyzzyaaav72,xyzzyaaaw72,xyzzyaaax72,xyzzyaaay72(s&
&ize_ee),xyzzyaaaz72(size_ee),xyzzyaaba72(size_en),xyzzyaabb72(size_en&
&),xyzzyaabc72
jas=0.d0
xyzzyaaad72(1)=0
do while(iterate_nuclei_indices(rank_e,rank_n,netot,nitot,nzencut,uncl&
&uster_en,xyzzyaaam72,xyzzyaaan72,xyzzyaaap72,xyzzyaaaq72,xyzzyaaas72,&
&xyzzyaaad72))
xyzzyaaac72(1)=0
do while(iterate_electron_indices(rank_e,netot,nzeecut,uncluster_ee,xy&
&zzyaaao72,xyzzyaaar72,xyzzyaaat72,xyzzyaaac72,preinit=xyzzyaaam72(1,r&
&ank_n),preinit_max=xyzzyaaap72(rank_n)))
do xyzzyaaag72=1,rank_e
xyzzyaaab72(xyzzyaaag72)=which_spin(xyzzyaaac72(xyzzyaaag72))
enddo
call get_sig(rank_e,rank_n,groups_ee,groups_en,xyzzyaaab72,xyzzyaaad72&
&,xyzzyaaah72,xyzzyaaae72,xyzzyaaaf72)
call match_signature(size_sig,xyzzyaaae1,nchannel,xyzzyaaah72,chsig,xy&
&zzyaaai72)
if(xyzzyaaai72<1)cycle
xyzzyaaau72=1.d0
xyzzyaaav72=1.d0
xyzzyaaaw72=1.d0
xyzzyaaax72=1.d0
xyzzyaaal72=iparam0(xyzzyaaai72)
do xyzzyaaaj72=1,nparam_remap(xyzzyaaai72)
xyzzyaaaa72=index_remap(xyzzyaaal72+xyzzyaaaj72)
xyzzyaaak72=xyzzyaaaa72-xyzzyaaal72
xyzzyaabc72=param(xyzzyaaaa72)
if(mask_leftmost_change_ee(xyzzyaaaa72)>0)call xyzzyaaiw1(nfn_eebasis,&
&netot,rank_e,rank_e,size_ee,mask_leftmost_change_ee(xyzzyaaaa72),xyzz&
&yaaae72,index_list(1,xyzzyaaak72),xyzzyaaac72,xyzzyaaac72,.true.,whic&
&h_ee_pair,eecut,xyzzyaaaz72,xyzzyaaav72)
if(mask_leftmost_change_en(xyzzyaaaa72)>0)call xyzzyaaiw1(nfn_enbasis,&
&nitot,rank_e,rank_n,size_en,mask_leftmost_change_en(xyzzyaaaa72),xyzz&
&yaaaf72,index_list(size_ee+1,xyzzyaaak72),xyzzyaaac72,xyzzyaaad72,.tr&
&ue.,which_en_pair,encut,xyzzyaabb72,xyzzyaaax72)
if(leftmost_change_ee(xyzzyaaaa72)>0)call xyzzyaaiw1(nfn_eebasis,netot&
&,rank_e,rank_e,size_ee,leftmost_change_ee(xyzzyaaaa72),xyzzyaaae72,in&
&dex_list(1,xyzzyaaak72),xyzzyaaac72,xyzzyaaac72,.false.,which_ee_pair&
&,eebasis,xyzzyaaay72,xyzzyaaau72)
if(leftmost_change_en(xyzzyaaaa72)>0)call xyzzyaaiw1(nfn_enbasis,nitot&
&,rank_e,rank_n,size_en,leftmost_change_en(xyzzyaaaa72),xyzzyaaaf72,in&
&dex_list(size_ee+1,xyzzyaaak72),xyzzyaaac72,xyzzyaaad72,.false.,which&
&_en_pair,enbasis,xyzzyaaba72,xyzzyaaaw72)
jas=jas+xyzzyaabc72*xyzzyaaav72*xyzzyaaax72*xyzzyaaau72*xyzzyaaaw72
enddo
enddo
enddo
end subroutine xyzzyaaiq1
subroutine xyzzyaair1(rank_e,size_ee,index_list,nparam_remap,index_rem&
&ap,leftmost_change_ee,mask_leftmost_change_ee,uncluster_ee,groups_ee,&
&param,nchannel,chsig,iparam0,which_ee_pair,eebasis,eecut,nzeecut,jas)
implicit none
integer,intent(in) :: rank_e,size_ee,groups_ee(nspin,nspin),nchannel
integer,intent(in) :: index_list(size_ee,*),nparam_remap(nchannel),ind&
&ex_remap(*),leftmost_change_ee(*),mask_leftmost_change_ee(*)
integer,intent(in) :: chsig(xyzzyaaae1,*),iparam0(*),which_ee_pair(2,*&
&)
real(dp),intent(in) :: eebasis(nfn_eebasis,netot,*),eecut(nfn_eebasis,&
&netot,*),param(*)
real(dp),intent(inout) :: jas
logical,intent(in) :: nzeecut(netot,netot),uncluster_ee
integer xyzzyaaaa73,xyzzyaaab73(rank_e),xyzzyaaac73(rank_e),xyzzyaaad7&
&3(size_ee),xyzzyaaae73(size_ee),xyzzyaaaf73,xyzzyaaag73,xyzzyaaah73,x&
&yzzyaaai73,xyzzyaaaj73,xyzzyaaak73(netot,rank_e),xyzzyaaal73(rank_e),&
&xyzzyaaam73(rank_e)
real(dp) xyzzyaaan73,xyzzyaaao73,xyzzyaaap73(size_ee),xyzzyaaaq73(size&
&_ee),xyzzyaaar73
jas=0.d0
xyzzyaaac73(1)=0
do while(iterate_electron_indices(rank_e,netot,nzeecut,uncluster_ee,xy&
&zzyaaak73,xyzzyaaal73,xyzzyaaam73,xyzzyaaac73))
do xyzzyaaaf73=1,rank_e
xyzzyaaab73(xyzzyaaaf73)=which_spin(xyzzyaaac73(xyzzyaaaf73))
enddo
call get_sig_ee_only(rank_e,groups_ee,xyzzyaaab73,xyzzyaaad73,xyzzyaaa&
&e73)
call match_signature(size_ee,xyzzyaaae1,nchannel,xyzzyaaad73,chsig,xyz&
&zyaaag73)
if(xyzzyaaag73<1)cycle
xyzzyaaan73=1.d0
xyzzyaaao73=1.d0
xyzzyaaaj73=iparam0(xyzzyaaag73)
do xyzzyaaah73=1,nparam_remap(xyzzyaaag73)
xyzzyaaaa73=index_remap(xyzzyaaaj73+xyzzyaaah73)
xyzzyaaai73=xyzzyaaaa73-xyzzyaaaj73
xyzzyaaar73=param(xyzzyaaaa73)
if(mask_leftmost_change_ee(xyzzyaaaa73)>0)call xyzzyaaiw1(nfn_eebasis,&
&netot,rank_e,rank_e,size_ee,mask_leftmost_change_ee(xyzzyaaaa73),xyzz&
&yaaae73,index_list(1,xyzzyaaai73),xyzzyaaac73,xyzzyaaac73,.true.,whic&
&h_ee_pair,eecut,xyzzyaaaq73,xyzzyaaao73)
if(leftmost_change_ee(xyzzyaaaa73)>0)call xyzzyaaiw1(nfn_eebasis,netot&
&,rank_e,rank_e,size_ee,leftmost_change_ee(xyzzyaaaa73),xyzzyaaae73,in&
&dex_list(1,xyzzyaaai73),xyzzyaaac73,xyzzyaaac73,.false.,which_ee_pair&
&,eebasis,xyzzyaaap73,xyzzyaaan73)
jas=jas+xyzzyaaar73*xyzzyaaao73*xyzzyaaan73
enddo
enddo
end subroutine xyzzyaair1
subroutine xyzzyaais1(rank_n,size_en,index_list,nparam_remap,index_rem&
&ap,leftmost_change_en,mask_leftmost_change_en,uncluster_en,groups_en,&
&param,nchannel,chsig,iparam0,which_en_pair,enbasis,encut,nzencut,jas)
implicit none
integer,intent(in) :: rank_n,size_en,groups_en(nitot,nspin),nchannel
integer,intent(in) :: index_list(size_en,*),nparam_remap(nchannel),ind&
&ex_remap(*),leftmost_change_en(*),mask_leftmost_change_en(*)
integer,intent(in) :: chsig(xyzzyaaae1,*),iparam0(*),which_en_pair(2,*&
&)
real(dp),intent(in) :: enbasis(nfn_enbasis,nitot,*),encut(nfn_enbasis,&
&nitot,*),param(*)
real(dp),intent(inout) :: jas
logical,intent(in) :: nzencut(nitot,netot),uncluster_en
integer xyzzyaaaa74,xyzzyaaab74,xyzzyaaac74,xyzzyaaad74,xyzzyaaae74(ra&
&nk_n),xyzzyaaaf74(size_en),xyzzyaaag74(size_en),xyzzyaaah74,xyzzyaaai&
&74,xyzzyaaaj74,xyzzyaaak74
integer xyzzyaaal74(netot,rank_n),xyzzyaaam74(nitot,rank_n),xyzzyaaan7&
&4(rank_n),xyzzyaaao74(rank_n),xyzzyaaap74(rank_n)
real(dp) xyzzyaaaq74,xyzzyaaar74,xyzzyaaas74(size_en),xyzzyaaat74(size&
&_en),xyzzyaaau74
jas=0.d0
xyzzyaaae74(1)=0
do while(iterate_nuclei_indices(1,rank_n,netot,nitot,nzencut,uncluster&
&_en,xyzzyaaal74,xyzzyaaam74,xyzzyaaan74,xyzzyaaao74,xyzzyaaap74,xyzzy&
&aaae74))
do xyzzyaaac74=1,xyzzyaaan74(rank_n)
xyzzyaaad74=xyzzyaaal74(xyzzyaaac74,rank_n)
xyzzyaaab74=which_spin(xyzzyaaad74)
call get_sig_en_only(rank_n,groups_en,xyzzyaaab74,xyzzyaaae74,xyzzyaaa&
&f74,xyzzyaaag74)
call match_signature(size_en,xyzzyaaae1,nchannel,xyzzyaaaf74,chsig,xyz&
&zyaaah74)
if(xyzzyaaah74<1)cycle
xyzzyaaaq74=1.d0
xyzzyaaar74=1.d0
xyzzyaaak74=iparam0(xyzzyaaah74)
do xyzzyaaai74=1,nparam_remap(xyzzyaaah74)
xyzzyaaaa74=index_remap(xyzzyaaak74+xyzzyaaai74)
xyzzyaaaj74=xyzzyaaaa74-xyzzyaaak74
xyzzyaaau74=param(xyzzyaaaa74)
if(mask_leftmost_change_en(xyzzyaaaa74)>0)call xyzzyaaiw1(nfn_enbasis,&
&nitot,1,rank_n,size_en,mask_leftmost_change_en(xyzzyaaaa74),xyzzyaaag&
&74,index_list(1,xyzzyaaaj74),(/xyzzyaaad74/),xyzzyaaae74,.true.,which&
&_en_pair,encut,xyzzyaaat74,xyzzyaaar74)
if(leftmost_change_en(xyzzyaaaa74)>0)call xyzzyaaiw1(nfn_enbasis,nitot&
&,1,rank_n,size_en,leftmost_change_en(xyzzyaaaa74),xyzzyaaag74,index_l&
&ist(1,xyzzyaaaj74),(/xyzzyaaad74/),xyzzyaaae74,.false.,which_en_pair,&
&enbasis,xyzzyaaas74,xyzzyaaaq74)
jas=jas+xyzzyaaau74*xyzzyaaar74*xyzzyaaaq74
enddo
enddo
enddo
end subroutine xyzzyaais1
subroutine xyzzyaait1(nfn,ne2,rank1,rank2,size12,leftmost_change,perm,&
&indx,ii,i1_vector,i2_vector,same12,is01,which_pair,basis,basis1,grad_&
&basis1,lap_basis1,incr_prod,incr_grad,incr_lap,prod,grad,lap)
implicit none
integer,intent(in) :: nfn,ne2,rank1,rank2,size12,leftmost_change,perm(&
&size12),indx(size12),ii,i1_vector(rank1),i2_vector(rank2),which_pair(&
&2,size12)
real(dp),intent(in) :: basis(nfn,ne2,*),basis1(nfn,*),grad_basis1(3,nf&
&n,*),lap_basis1(nfn,*)
real(dp),intent(inout) :: prod,incr_prod(size12),grad(3),incr_grad(3,s&
&ize12),lap,incr_lap(size12)
logical,intent(in) :: same12,is01
integer xyzzyaaaa75,xyzzyaaab75,xyzzyaaac75,xyzzyaaad75,xyzzyaaae75
real(dp) xyzzyaaaf75,xyzzyaaag75,xyzzyaaah75(3),xyzzyaaai75
if(leftmost_change<1)return
if(leftmost_change>1)then
prod=incr_prod(leftmost_change-1)
grad=incr_grad(1:3,leftmost_change-1)
lap=incr_lap(leftmost_change-1)
else
prod=1.d0
grad=0.d0
lap=0.d0
endif
do xyzzyaaaa75=leftmost_change,size12
xyzzyaaac75=indx(xyzzyaaaa75)
if(xyzzyaaac75>0)then
if(is01)xyzzyaaac75=1
xyzzyaaab75=perm(xyzzyaaaa75)
xyzzyaaad75=i1_vector(which_pair(1,xyzzyaaab75))
xyzzyaaae75=i2_vector(which_pair(2,xyzzyaaab75))
if(xyzzyaaad75==ii)then
xyzzyaaaf75=basis1(xyzzyaaac75,xyzzyaaae75)
if(xyzzyaaaf75/=0.d0)then
xyzzyaaai75=1.d0/xyzzyaaaf75
xyzzyaaah75=grad_basis1(1:3,xyzzyaaac75,xyzzyaaae75)*xyzzyaaai75
grad=grad+xyzzyaaah75
xyzzyaaag75=lap_basis1(xyzzyaaac75,xyzzyaaae75)
lap=lap+xyzzyaaag75*xyzzyaaai75-ddot(3,xyzzyaaah75,1,xyzzyaaah75,1)
endif
elseif(same12.and.xyzzyaaae75==ii)then
xyzzyaaaf75=basis1(xyzzyaaac75,xyzzyaaad75)
if(xyzzyaaaf75/=0.d0)then
xyzzyaaai75=1.d0/xyzzyaaaf75
xyzzyaaah75=-grad_basis1(1:3,xyzzyaaac75,xyzzyaaad75)*xyzzyaaai75
grad=grad+xyzzyaaah75
xyzzyaaag75=lap_basis1(xyzzyaaac75,xyzzyaaad75)
lap=lap+xyzzyaaag75*xyzzyaaai75-ddot(3,xyzzyaaah75,1,xyzzyaaah75,1)
endif
else
xyzzyaaaf75=basis(xyzzyaaac75,xyzzyaaae75,xyzzyaaad75)
endif
prod=prod*xyzzyaaaf75
endif
incr_prod(xyzzyaaaa75)=prod
incr_grad(1:3,xyzzyaaaa75)=grad(1:3)
incr_lap(xyzzyaaaa75)=lap
enddo
end subroutine xyzzyaait1
subroutine xyzzyaaiu1(nfn,ne2,rank1,rank2,size12,leftmost_change,perm,&
&indx,ii,i1_vector,i2_vector,same12,is01,which_pair,basis,basis1,grad_&
&basis1,incr_prod,incr_grad,prod,grad)
implicit none
integer,intent(in) :: nfn,ne2,rank1,rank2,size12,leftmost_change,perm(&
&size12),indx(size12),ii,i1_vector(rank1),i2_vector(rank2),which_pair(&
&2,size12)
real(dp),intent(in) :: basis(nfn,ne2,*),basis1(nfn,*),grad_basis1(3,nf&
&n,*)
real(dp),intent(inout) :: prod,incr_prod(size12),grad(3),incr_grad(3,s&
&ize12)
logical,intent(in) :: same12,is01
integer xyzzyaaaa76,xyzzyaaab76,xyzzyaaac76,xyzzyaaad76,xyzzyaaae76
real(dp) xyzzyaaaf76,xyzzyaaag76
if(leftmost_change<1)return
if(leftmost_change>1)then
prod=incr_prod(leftmost_change-1)
grad=incr_grad(1:3,leftmost_change-1)
else
prod=1.d0
grad=0.d0
endif
do xyzzyaaaa76=leftmost_change,size12
xyzzyaaac76=indx(xyzzyaaaa76)
if(xyzzyaaac76>0)then
if(is01)xyzzyaaac76=1
xyzzyaaab76=perm(xyzzyaaaa76)
xyzzyaaad76=i1_vector(which_pair(1,xyzzyaaab76))
xyzzyaaae76=i2_vector(which_pair(2,xyzzyaaab76))
if(xyzzyaaad76==ii)then
xyzzyaaaf76=basis1(xyzzyaaac76,xyzzyaaae76)
if(xyzzyaaaf76/=0.d0)then
xyzzyaaag76=1.d0/xyzzyaaaf76
grad=grad+grad_basis1(1:3,xyzzyaaac76,xyzzyaaae76)*xyzzyaaag76
endif
elseif(same12.and.xyzzyaaae76==ii)then
xyzzyaaaf76=basis1(xyzzyaaac76,xyzzyaaad76)
if(xyzzyaaaf76/=0.d0)then
xyzzyaaag76=1.d0/xyzzyaaaf76
grad=grad-grad_basis1(1:3,xyzzyaaac76,xyzzyaaad76)*xyzzyaaag76
endif
else
xyzzyaaaf76=basis(xyzzyaaac76,xyzzyaaae76,xyzzyaaad76)
endif
prod=prod*xyzzyaaaf76
endif
incr_prod(xyzzyaaaa76)=prod
incr_grad(1:3,xyzzyaaaa76)=grad(1:3)
enddo
end subroutine xyzzyaaiu1
subroutine xyzzyaaiv1(nfn,ne2,rank1,rank2,size12,leftmost_change,perm,&
&indx,ii,i1_vector,i2_vector,same12,is01,which_pair,basis,basis1,incr_&
&prod,prod)
implicit none
integer,intent(in) :: nfn,ne2,rank1,rank2,size12,leftmost_change,perm(&
&size12),indx(size12),ii,i1_vector(rank1),i2_vector(rank2),which_pair(&
&2,size12)
real(dp),intent(in) :: basis(nfn,ne2,*),basis1(nfn,*)
real(dp),intent(inout) :: prod,incr_prod(size12)
logical,intent(in) :: same12,is01
integer xyzzyaaaa77,xyzzyaaab77,xyzzyaaac77,xyzzyaaad77,xyzzyaaae77
real(dp) xyzzyaaaf77
if(leftmost_change<1)return
if(leftmost_change>1)then
prod=incr_prod(leftmost_change-1)
else
prod=1.d0
endif
do xyzzyaaaa77=leftmost_change,size12
xyzzyaaac77=indx(xyzzyaaaa77)
if(xyzzyaaac77>0)then
if(is01)xyzzyaaac77=1
xyzzyaaab77=perm(xyzzyaaaa77)
xyzzyaaad77=i1_vector(which_pair(1,xyzzyaaab77))
xyzzyaaae77=i2_vector(which_pair(2,xyzzyaaab77))
if(xyzzyaaad77==ii)then
xyzzyaaaf77=basis1(xyzzyaaac77,xyzzyaaae77)
elseif(same12.and.xyzzyaaae77==ii)then
xyzzyaaaf77=basis1(xyzzyaaac77,xyzzyaaad77)
else
xyzzyaaaf77=basis(xyzzyaaac77,xyzzyaaae77,xyzzyaaad77)
endif
prod=prod*xyzzyaaaf77
endif
incr_prod(xyzzyaaaa77)=prod
enddo
end subroutine xyzzyaaiv1
subroutine xyzzyaaiw1(nfn,ne2,rank1,rank2,size12,leftmost_change,perm,&
&indx,i1_vector,i2_vector,is01,which_pair,basis,incr_prod,prod)
implicit none
integer,intent(in) :: nfn,ne2,rank1,rank2,size12,leftmost_change,perm(&
&size12),indx(size12),i1_vector(rank1),i2_vector(rank2),which_pair(2,s&
&ize12)
real(dp),intent(in) :: basis(nfn,ne2,*)
real(dp),intent(inout) :: prod,incr_prod(size12)
logical,intent(in) :: is01
integer xyzzyaaaa78,xyzzyaaab78,xyzzyaaac78,xyzzyaaad78,xyzzyaaae78
if(leftmost_change<1)return
if(leftmost_change>1)then
prod=incr_prod(leftmost_change-1)
else
prod=1.d0
endif
do xyzzyaaaa78=leftmost_change,size12
xyzzyaaac78=indx(xyzzyaaaa78)
if(xyzzyaaac78>0)then
if(is01)xyzzyaaac78=1
xyzzyaaab78=perm(xyzzyaaaa78)
xyzzyaaad78=i1_vector(which_pair(1,xyzzyaaab78))
xyzzyaaae78=i2_vector(which_pair(2,xyzzyaaab78))
prod=prod*basis(xyzzyaaac78,xyzzyaaae78,xyzzyaaad78)
endif
incr_prod(xyzzyaaaa78)=prod
enddo
end subroutine xyzzyaaiw1
subroutine gjastrow_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,v&
&erbose)
implicit none
real(dp),intent(in) :: eevecs(4,netot,netot),eivecs(4,nitot,netot)
integer,intent(out) :: ie,jspin
logical,intent(in) :: verbose
logical,intent(out) :: fail
ie=1
jspin=1
fail=.true.
end subroutine gjastrow_assess_check_kinetic
subroutine finite_size_corr_ke_gjastrow
implicit none
end subroutine finite_size_corr_ke_gjastrow
subroutine enumerate_plot_gjastrow(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
integer xyzzyaaaa81,xyzzyaaab81
logical xyzzyaaac81,xyzzyaaad81
xyzzyaaad81=.not.(present(keyword).and.present(description))
xyzzyaaab81=0
xyzzyaahb1=0
do xyzzyaaaa81=1,nthings_to_plot
xyzzyaaac81=.false.
select case(xyzzyaaaa81)
case(xyzzyaagi1)
xyzzyaaac81=.true.
case(xyzzyaagj1)
xyzzyaaac81=.true.
case(xyzzyaagk1)
xyzzyaaac81=.true.
case(xyzzyaagl1)
xyzzyaaac81=.true.
case(xyzzyaagm1)
xyzzyaaac81=.true.
case(xyzzyaagn1)
xyzzyaaac81=.true.
case(xyzzyaago1)
xyzzyaaac81=.true.
case(xyzzyaagp1)
if(any(xyzzyaaao1>0))xyzzyaaac81=.true.
case(xyzzyaagq1)
if(any(xyzzyaaao1>0))xyzzyaaac81=.true.
case(xyzzyaagr1)
if(any(xyzzyaaao1>0))xyzzyaaac81=.true.
case(xyzzyaags1)
if(any(xyzzyaaap1>0))xyzzyaaac81=.true.
case(xyzzyaagt1)
if(any(xyzzyaaap1>0))xyzzyaaac81=.true.
case(xyzzyaagu1)
if(any(xyzzyaaap1>0))xyzzyaaac81=.true.
case(xyzzyaagv1)
if(any(xyzzyaaao1>0))xyzzyaaac81=.true.
case(xyzzyaagw1)
if(any(xyzzyaaao1>0))xyzzyaaac81=.true.
case(xyzzyaagx1)
if(any(xyzzyaaao1>0))xyzzyaaac81=.true.
case(xyzzyaagy1)
if(any(xyzzyaaap1>0))xyzzyaaac81=.true.
case(xyzzyaagz1)
if(any(xyzzyaaap1>0))xyzzyaaac81=.true.
case(xyzzyaaha1)
if(any(xyzzyaaap1>0))xyzzyaaac81=.true.
end select
if(xyzzyaaac81)then
xyzzyaaab81=xyzzyaaab81+1
xyzzyaahb1(xyzzyaaab81)=xyzzyaaaa81
endif
enddo
n=xyzzyaaab81
if(.not.xyzzyaaad81)then
do xyzzyaaaa81=1,xyzzyaaab81
keyword(xyzzyaaaa81)=xyzzyaahc1(xyzzyaahb1(xyzzyaaaa81))
description(xyzzyaaaa81)=plot_description(xyzzyaahb1(xyzzyaaaa81))
enddo
endif
end subroutine enumerate_plot_gjastrow
subroutine query_plot_gjastrow(iplot,ii,rank,is_complex,has_stderr,rot&
&_tensor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
integer iterm,indx,jj,ion
logical count_only
count_only=.not.present(function_name)
if(iplot<0.or.iplot>nthings_to_plot)call errstop_master('QUERY_PLOT_GJ&
&ASTROW','IPLOT out of range. Bug in calling routine.')
if(xyzzyaahb1(iplot)==0)call errstop_master('QUERY_PLOT_GJASTROW','IPL&
&OT out of range. Bug in calling routine.')
select case(xyzzyaahb1(iplot))
case(xyzzyaagi1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Jastrow factor value'
endif
case(xyzzyaagj1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='One-electron Jastrow factor value'
endif
case(xyzzyaagk1)
rank=1
is_complex=.false.
has_stderr=.false.
rot_tensor=.true.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Jastrow factor gradient'
endif
case(xyzzyaagl1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Jastrow factor Laplacian'
endif
case(xyzzyaagm1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Term '//trim(i2s(iterm))//' value'
endif
enddo
case(xyzzyaagn1)
rank=1
is_complex=.false.
has_stderr=.false.
rot_tensor=.true.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Term '//trim(i2s(iterm))//' gradient'
endif
enddo
case(xyzzyaago1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Term '//trim(i2s(iterm))//' Laplacian'
endif
enddo
case(xyzzyaagp1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaao1(iterm)==0)cycle
do indx=1,xyzzyaaam1(iterm)
do jj=1,netot
if(jj==ii)cycle
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='value of e-e basis function #'//trim(i2s(in&
&dx))//' of term '//trim(i2s(iterm))//' for particle '//trim(i2s(jj))
endif
enddo
enddo
enddo
case(xyzzyaagq1)
rank=1
is_complex=.false.
has_stderr=.false.
rot_tensor=.true.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaao1(iterm)==0)cycle
do indx=1,xyzzyaaam1(iterm)
do jj=1,netot
if(jj==ii)cycle
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='gradient of e-e basis function #'//trim(i2s&
&(indx))//' of term '//trim(i2s(iterm))//' for particle '//trim(i2s(jj&
&))
endif
enddo
enddo
enddo
case(xyzzyaagr1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaao1(iterm)==0)cycle
do indx=1,xyzzyaaam1(iterm)
do jj=1,netot
if(jj==ii)cycle
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Laplacian of e-e basis function #'//trim(i2&
&s(indx))//' of term '//trim(i2s(iterm))//' for particle '//trim(i2s(j&
&j))
endif
enddo
enddo
enddo
case(xyzzyaags1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaap1(iterm)==0)cycle
do indx=1,xyzzyaaan1(iterm)
do ion=1,nitot
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='value of e-n basis function #'//trim(i2s(in&
&dx))//' of term '//trim(i2s(iterm))//' for nucleus '//trim(i2s(ion))
endif
enddo
enddo
enddo
case(xyzzyaagt1)
rank=1
is_complex=.false.
has_stderr=.false.
rot_tensor=.true.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaap1(iterm)==0)cycle
do indx=1,xyzzyaaan1(iterm)
do ion=1,nitot
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='gradient of e-n basis function #'//trim(i2s&
&(indx))//' of term '//trim(i2s(iterm))//' for nucleus '//trim(i2s(ion&
&))
endif
enddo
enddo
enddo
case(xyzzyaagu1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaap1(iterm)==0)cycle
do indx=1,xyzzyaaan1(iterm)
do ion=1,nitot
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='laplacian of e-n basis function #'//trim(i2&
&s(indx))//' of term '//trim(i2s(iterm))//' for nucleus '//trim(i2s(io&
&n))
endif
enddo
enddo
enddo
case(xyzzyaagv1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaao1(iterm)==0)cycle
do jj=1,netot
if(jj==ii)cycle
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='value of e-e cut-off of term '//trim(i2s(it&
&erm))//' for particle '//trim(i2s(jj))
endif
enddo
enddo
case(xyzzyaagw1)
rank=1
is_complex=.false.
has_stderr=.false.
rot_tensor=.true.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaao1(iterm)==0)cycle
do jj=1,netot
if(jj==ii)cycle
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='gradient of e-e cut-off of term '//trim(i2s&
&(iterm))//' for particle '//trim(i2s(jj))
endif
enddo
enddo
case(xyzzyaagx1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaao1(iterm)==0)cycle
do jj=1,netot
if(jj==ii)cycle
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Laplacian of e-e cut-off of term '//trim(i2&
&s(iterm))//' for particle '//trim(i2s(jj))
endif
enddo
enddo
case(xyzzyaagy1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaap1(iterm)==0)cycle
do ion=1,nitot
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='value of e-n cut-off of term '//trim(i2s(it&
&erm))//' for nucleus '//trim(i2s(ion))
endif
enddo
enddo
case(xyzzyaagz1)
rank=1
is_complex=.false.
has_stderr=.false.
rot_tensor=.true.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaap1(iterm)==0)cycle
do ion=1,nitot
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='gradient of e-n cut-off of term '//trim(i2s&
&(iterm))//' for nucleus '//trim(i2s(ion))
endif
enddo
enddo
case(xyzzyaaha1)
rank=0
is_complex=.false.
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
do iterm=1,xyzzyaaaa1
if(xyzzyaaap1(iterm)==0)cycle
do ion=1,nitot
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Laplacian of e-n cut-off of term '//trim(i2&
&s(iterm))//' for nucleus '//trim(i2s(ion))
endif
enddo
enddo
end select
end subroutine query_plot_gjastrow
subroutine get_plot_gjastrow(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
integer xyzzyaaaa83,xyzzyaaab83,xyzzyaaac83,xyzzyaaad83,xyzzyaaae83,xy&
&zzyaaaf83,xyzzyaaag83,xyzzyaaah83
if(iplot<0.or.iplot>nthings_to_plot)call errstop_master('GET_PLOT_GJAS&
&TROW','IPLOT out of range. Bug in calling routine.')
if(xyzzyaahb1(iplot)==0)call errstop_master('GET_PLOT_GJASTROW','IPLOT&
& out of range. Bug in calling routine.')
select case(xyzzyaahb1(iplot))
case(xyzzyaagi1)
call xyzzyaaho1(is1)
f(1)=xyzzyaadx1(is1)
case(xyzzyaagj1)
call xyzzyaahn1(ii,is0,.true.,.false.,.false.,.false.)
call xyzzyaahn1(ii,is1,.true.,.false.,.false.,.false.)
f(1)=xyzzyaaec1(ii,is1)
case(xyzzyaagk1)
call xyzzyaahn1(ii,is0,.true.,.true.,.false.,.false.)
call xyzzyaahn1(ii,is1,.true.,.true.,.false.,.false.)
f(1:dimensionality)=xyzzyaadz1(1:dimensionality,ii,is1)
case(xyzzyaagl1)
call xyzzyaahn1(ii,is0,.true.,.true.,.true.,.false.)
call xyzzyaahn1(ii,is1,.true.,.true.,.true.,.false.)
f(1)=xyzzyaaeb1(ii,is1)
case(xyzzyaagm1)
call xyzzyaahn1(ii,is0,.true.,.true.,.true.,.false.)
call xyzzyaahn1(ii,is1,.true.,.true.,.true.,.false.)
do xyzzyaaaa83=1,xyzzyaaaa1
f(xyzzyaaaa83)=xyzzyaaeh1(ii,xyzzyaaaa83,is1)
enddo
case(xyzzyaagn1)
call xyzzyaahn1(ii,is0,.true.,.true.,.true.,.false.)
call xyzzyaahn1(ii,is1,.true.,.true.,.true.,.false.)
do xyzzyaaaa83=1,xyzzyaaaa1
f((xyzzyaaaa83-1)*dimensionality+1:xyzzyaaaa83*dimensionality)=xyzzyaa&
&ei1(1:dimensionality,ii,xyzzyaaaa83,is1)
enddo
case(xyzzyaago1)
call xyzzyaahn1(ii,is0,.true.,.true.,.true.,.false.)
call xyzzyaahn1(ii,is1,.true.,.true.,.true.,.false.)
do xyzzyaaaa83=1,xyzzyaaaa1
f(xyzzyaaaa83)=xyzzyaaej1(ii,xyzzyaaaa83,is1)
enddo
case(xyzzyaagp1)
call get_gbasis(is0,.false.,.false.)
call get_gbasis(is1,.false.,.false.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaao1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaaq1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaag83=ifn1_eebasis(xyzzyaaab83)-1
do xyzzyaaac83=1,xyzzyaaam1(xyzzyaaaa83)
do xyzzyaaad83=1,netot
if(xyzzyaaad83/=ii)cycle
xyzzyaaaf83=xyzzyaaaf83+1
f(xyzzyaaaf83)=eebasis_scr(xyzzyaaag83+xyzzyaaac83,xyzzyaaad83,ii,is1)
enddo
enddo
enddo
case(xyzzyaagq1)
call get_gbasis(is0,.true.,.false.)
call get_gbasis(is1,.true.,.false.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaao1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaaq1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaag83=ifn1_eebasis(xyzzyaaab83)-1
do xyzzyaaac83=1,xyzzyaaam1(xyzzyaaaa83)
do xyzzyaaad83=1,netot
if(xyzzyaaad83/=ii)cycle
f(xyzzyaaaf83+1:xyzzyaaaf83+dimensionality)=grad_eebasis_scr(1:dimensi&
&onality,xyzzyaaag83+xyzzyaaac83,xyzzyaaad83,ii,is1)
xyzzyaaaf83=xyzzyaaaf83+dimensionality
enddo
enddo
enddo
case(xyzzyaagr1)
call get_gbasis(is0,.true.,.true.)
call get_gbasis(is1,.true.,.true.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaao1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaaq1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaag83=ifn1_eebasis(xyzzyaaab83)-1
do xyzzyaaac83=1,xyzzyaaam1(xyzzyaaaa83)
do xyzzyaaad83=1,netot
if(xyzzyaaad83/=ii)cycle
xyzzyaaaf83=xyzzyaaaf83+1
f(xyzzyaaaf83)=lap_eebasis_scr(xyzzyaaag83+xyzzyaaac83,xyzzyaaad83,ii,&
&is1)
enddo
enddo
enddo
case(xyzzyaags1)
call get_gbasis(is0,.false.,.false.)
call get_gbasis(is1,.false.,.false.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaap1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaas1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaag83=ifn1_enbasis(xyzzyaaab83)-1
do xyzzyaaac83=1,xyzzyaaan1(xyzzyaaaa83)
do xyzzyaaae83=1,nitot
xyzzyaaaf83=xyzzyaaaf83+1
f(xyzzyaaaf83)=enbasis_scr(xyzzyaaag83+xyzzyaaac83,xyzzyaaae83,ii,is1)
enddo
enddo
enddo
case(xyzzyaagt1)
call get_gbasis(is0,.true.,.false.)
call get_gbasis(is1,.true.,.false.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaap1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaas1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaag83=ifn1_enbasis(xyzzyaaab83)-1
do xyzzyaaac83=1,xyzzyaaan1(xyzzyaaaa83)
do xyzzyaaae83=1,nitot
f(xyzzyaaaf83+1:xyzzyaaaf83+dimensionality)=grad_enbasis_scr(1:dimensi&
&onality,xyzzyaaag83+xyzzyaaac83,xyzzyaaae83,ii,is1)
xyzzyaaaf83=xyzzyaaaf83+dimensionality
enddo
enddo
enddo
case(xyzzyaagu1)
call get_gbasis(is0,.true.,.true.)
call get_gbasis(is1,.true.,.true.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaap1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaas1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaag83=ifn1_enbasis(xyzzyaaab83)-1
do xyzzyaaac83=1,xyzzyaaan1(xyzzyaaaa83)
do xyzzyaaae83=1,nitot
xyzzyaaaf83=xyzzyaaaf83+1
f(xyzzyaaaf83)=lap_enbasis_scr(xyzzyaaag83+xyzzyaaac83,xyzzyaaae83,ii,&
&is1)
enddo
enddo
enddo
case(xyzzyaagv1)
call get_gbasis(is0,.false.,.false.)
call get_gbasis(is1,.false.,.false.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaao1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaar1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaah83=ifn1_eebasis(xyzzyaaab83)
do xyzzyaaad83=1,netot
if(xyzzyaaad83/=ii)cycle
xyzzyaaaf83=xyzzyaaaf83+1
f(xyzzyaaaf83)=eebasis_scr(xyzzyaaah83,xyzzyaaad83,ii,is1)
enddo
enddo
case(xyzzyaagw1)
call get_gbasis(is0,.true.,.false.)
call get_gbasis(is1,.true.,.false.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaao1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaar1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaah83=ifn1_eebasis(xyzzyaaab83)
do xyzzyaaad83=1,netot
if(xyzzyaaad83/=ii)cycle
f(xyzzyaaaf83+1:xyzzyaaaf83+dimensionality)=grad_eebasis_scr(1:dimensi&
&onality,xyzzyaaah83,xyzzyaaad83,ii,is1)
xyzzyaaaf83=xyzzyaaaf83+dimensionality
enddo
enddo
case(xyzzyaagx1)
call get_gbasis(is0,.true.,.true.)
call get_gbasis(is1,.true.,.true.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaao1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaar1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaah83=ifn1_eebasis(xyzzyaaab83)
do xyzzyaaad83=1,netot
if(xyzzyaaad83/=ii)cycle
xyzzyaaaf83=xyzzyaaaf83+1
f(xyzzyaaaf83)=lap_eebasis_scr(xyzzyaaah83,xyzzyaaad83,ii,is1)
enddo
enddo
case(xyzzyaagy1)
call get_gbasis(is0,.false.,.false.)
call get_gbasis(is1,.false.,.false.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaap1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaat1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaah83=ifn1_enbasis(xyzzyaaab83)
do xyzzyaaae83=1,nitot
xyzzyaaaf83=xyzzyaaaf83+1
f(xyzzyaaaf83)=enbasis_scr(xyzzyaaah83,xyzzyaaae83,ii,is1)
enddo
enddo
case(xyzzyaagz1)
call get_gbasis(is0,.true.,.false.)
call get_gbasis(is1,.true.,.false.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaap1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaat1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaah83=ifn1_enbasis(xyzzyaaab83)
do xyzzyaaae83=1,nitot
f(xyzzyaaaf83+1:xyzzyaaaf83+dimensionality)=grad_enbasis_scr(1:dimensi&
&onality,xyzzyaaah83,xyzzyaaae83,ii,is1)
xyzzyaaaf83=xyzzyaaaf83+dimensionality
enddo
enddo
case(xyzzyaaha1)
call get_gbasis(is0,.true.,.true.)
call get_gbasis(is1,.true.,.true.)
xyzzyaaaf83=0
do xyzzyaaaa83=1,xyzzyaaaa1
if(xyzzyaaap1(xyzzyaaaa83)==0)cycle
xyzzyaaab83=xyzzyaaat1(xyzzyaaaa83)
if(xyzzyaaab83==0)cycle
xyzzyaaah83=ifn1_enbasis(xyzzyaaab83)
do xyzzyaaae83=1,nitot
xyzzyaaaf83=xyzzyaaaf83+1
f(xyzzyaaaf83)=lap_enbasis_scr(xyzzyaaah83,xyzzyaaae83,ii,is1)
enddo
enddo
end select
end subroutine get_plot_gjastrow
subroutine finish_plot_gjastrow
implicit none
xyzzyaahb1=0
end subroutine finish_plot_gjastrow
end module slaarnabh

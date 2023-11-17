module slaarnabx
use slaarnaag
use dsp
use parallel
use slaarnach
use store
use file_utils,   only : open_units
use format_utils, only : wout,i2s,l2s,r2s,write_list_int,wordwrap,disp&
&lay_param
use slaarnaas,    only : me_biex3,mh_biex3,mu_biex3,xx_sep
use slaarnabg,     only : nitot,wigner_seitz_radius,isperiodic,dimensi&
&onality,pb1,pb2,pb3,b1,b2,b3,iontype,periodicity,npcells,model_system&
&,homogeneous_system,rionion,hard_sphere,hard_op_spins,hard_diam,nityp&
&e,nbasis
use slaarnabt,    only : exp_protect,dcopy,swap1
use slaarnaca,        only : zion
use run_control,  only : errstop,errwarn,errstop_master,errwarn_silent&
&,timer,check_alloc
implicit none
private
public read_pjastrow,write_pjastrow,get_linear_basis_pjastrow,pjastrow&
&_assess_check_kinetic,finite_size_corr_ke_pjastrow,setup_pjastrow_plo&
&t,enumerate_plot_pjastrow,query_plot_pjastrow,get_plot_pjastrow,finis&
&h_plot_pjastrow,check_varmin_linjas_pjastrow
public query_pjastrow_levels,query_pjastrow_level_details,setup_pjastr&
&ow,finish_pjastrow,wfn_ratio_pjastrow,accept_move_pjastrow,reset_conf&
&ig_pjastrow,add_config_pjastrow_items,clear_scratch_pjastrow,wfn_logv&
&al_pjastrow,wfn_loggrad_pjastrow,wfn_loglap_pjastrow,prefetch_wfn_pja&
&strow,setup_pjastrow_params,finish_pjastrow_params,get_pjastrow_param&
&s,put_pjastrow_params,clone_scratch_pjastrow,invalidate_params_pjastr&
&ow,invalidate_param1_pjastrow,setup_storage_pjastrow,finish_storage_p&
&jastrow,load_from_storage_pjastrow,save_to_storage_pjastrow
public gen_config_pjastrow,delete_config_pjastrow,copy_config_pjastrow&
&,config_to_pt_pjastrow,pt_to_config_pjastrow,redist_allocations_pjast&
&row,redist_load_pjastrow,redist_send_pjastrow,redist_recv_pjastrow,re&
&dist_save_pjastrow,redist_deallocations_pjastrow,load_from_pt_pjastro&
&w,save_to_pt_pjastrow
public config_wfn_pjastrow
logical,parameter :: xyzzyaaaa1=.false.
real(dp),allocatable :: xyzzyaaab1(:,:),xyzzyaaac1(:,:,:),xyzzyaaad1(:&
&,:,:),xyzzyaaae1(:,:,:,:,:),xyzzyaaaf1(:,:),xyzzyaaag1(:,:),xyzzyaaah&
&1(:,:),xyzzyaaai1(:,:),xyzzyaaaj1(:,:),xyzzyaaak1(:,:,:),xyzzyaaal1(:&
&,:,:),xyzzyaaam1(:,:,:),xyzzyaaan1(:),xyzzyaaao1(:),xyzzyaaap1(:,:),x&
&yzzyaaaq1(:,:),xyzzyaaar1(:,:,:,:),xyzzyaaas1(:,:,:,:)
integer,allocatable :: xyzzyaaat1(:,:),xyzzyaaau1(:,:,:),xyzzyaaav1(:,&
&:,:),xyzzyaaaw1(:,:,:,:,:),xyzzyaaax1(:,:),xyzzyaaay1(:,:),xyzzyaaaz1&
&(:),xyzzyaaba1(:),xyzzyaabb1(:,:),xyzzyaabc1(:,:,:,:),xyzzyaabd1(:,:,&
&:,:),xyzzyaabe1(:,:,:,:)
integer :: xyzzyaabf1,xyzzyaabg1
integer xyzzyaabh1,xyzzyaabi1,xyzzyaabj1,xyzzyaabk1,xyzzyaabl1
integer xyzzyaabm1,xyzzyaabn1
integer,allocatable :: xyzzyaabo1(:),xyzzyaabp1(:)
integer,allocatable :: xyzzyaabq1(:,:),xyzzyaabr1(:,:),xyzzyaabs1(:,:)&
&,xyzzyaabt1(:,:)
integer,allocatable :: xyzzyaabu1(:,:),xyzzyaabv1(:)
integer,allocatable :: xyzzyaabw1(:)
integer xyzzyaabx1,xyzzyaaby1,xyzzyaabz1,xyzzyaaca1,xyzzyaacb1,xyzzyaa&
&cc1,xyzzyaacd1
integer,allocatable :: xyzzyaace1(:),xyzzyaacf1(:),xyzzyaacg1(:)
real(dp) xyzzyaach1,xyzzyaaci1,xyzzyaacj1,xyzzyaack1,xyzzyaacl1,xyzzya&
&acm1,xyzzyaacn1,xyzzyaaco1,xyzzyaacp1,xyzzyaacq1,xyzzyaacr1(3),xyzzya&
&acs1,xyzzyaact1
real(dp),allocatable :: xyzzyaacu1(:),xyzzyaacv1(:),xyzzyaacw1(:,:,:)
integer xyzzyaacx1,xyzzyaacy1,xyzzyaacz1,xyzzyaada1,xyzzyaadb1,xyzzyaa&
&dc1,xyzzyaadd1,xyzzyaade1
integer,allocatable :: xyzzyaadf1(:),xyzzyaadg1(:),xyzzyaadh1(:,:,:)
integer,allocatable :: xyzzyaadi1(:),xyzzyaadj1(:),xyzzyaadk1(:),xyzzy&
&aadl1(:),xyzzyaadm1(:)
integer,allocatable :: xyzzyaadn1(:)
real(dp),allocatable :: xyzzyaado1(:),xyzzyaadp1(:),xyzzyaadq1(:,:),xy&
&zzyaadr1(:),xyzzyaads1(:,:,:),xyzzyaadt1(:),xyzzyaadu1(:),xyzzyaadv1(&
&:),xyzzyaadw1(:),xyzzyaadx1(:,:,:)
integer xyzzyaady1
!$omp threadprivate(xyzzyaado1,xyzzyaadp1,xyzzyaadq1,xyzzyaadr1,xyzzya&
!$omp &adt1,xyzzyaadu1,xyzzyaadv1,xyzzyaadw1,xyzzyaadx1)
real(dp),allocatable :: xyzzyaadz1(:,:),xyzzyaaea1(:,:),xyzzyaaeb1(:,:&
&,:),xyzzyaaec1(:,:),xyzzyaaed1(:),xyzzyaaee1(:),xyzzyaaef1(:,:,:),xyz&
&zyaaeg1(:),xyzzyaaeh1(:),xyzzyaaei1(:,:),xyzzyaaej1(:,:,:),xyzzyaaek1&
&(:,:),xyzzyaael1(:),xyzzyaaem1(:),xyzzyaaen1(:,:),xyzzyaaeo1(:)
logical xyzzyaaep1,xyzzyaaeq1,xyzzyaaer1,xyzzyaaes1,xyzzyaaet1,xyzzyaa&
&eu1,xyzzyaaev1,xyzzyaaew1,xyzzyaaex1,xyzzyaaey1,xyzzyaaez1,xyzzyaafa1&
&,xyzzyaafb1,xyzzyaafc1,xyzzyaafd1,xyzzyaafe1,xyzzyaaff1,xyzzyaafg1
integer xyzzyaafh1,xyzzyaafi1,xyzzyaafj1,xyzzyaafk1,xyzzyaafl1,xyzzyaa&
&fm1
integer,allocatable :: xyzzyaafn1(:),xyzzyaafo1(:)
character(80) title
integer xyzzyaafp1,xyzzyaafq1,xyzzyaafr1,xyzzyaafs1,xyzzyaaft1,xyzzyaa&
&fu1
integer,allocatable :: xyzzyaafv1(:),xyzzyaafw1(:)
real(dp),allocatable :: xyzzyaafx1(:)
integer xyzzyaafy1,xyzzyaafz1
integer,allocatable :: xyzzyaaga1(:,:),xyzzyaagb1(:),xyzzyaagc1(:,:),x&
&yzzyaagd1(:)
integer xyzzyaage1(3),xyzzyaagf1(3)
real(dp),allocatable :: xyzzyaagg1(:,:),xyzzyaagh1(:),xyzzyaagi1(:,:),&
&xyzzyaagj1(:)
complex(dp),allocatable :: xyzzyaagk1(:),xyzzyaagl1(:),xyzzyaagm1(:),x&
&yzzyaagn1(:),xyzzyaago1(:),xyzzyaagp1(:)
integer xyzzyaagq1,xyzzyaagr1,xyzzyaags1
integer,allocatable :: xyzzyaagt1(:),xyzzyaagu1(:)
integer,allocatable :: xyzzyaagv1(:)
integer xyzzyaagw1
integer xyzzyaagx1
integer,allocatable :: xyzzyaagy1(:),xyzzyaagz1(:)
integer xyzzyaaha1,xyzzyaahb1,xyzzyaahc1,xyzzyaahd1,xyzzyaahe1,xyzzyaa&
&hf1
integer,allocatable :: xyzzyaahg1(:),xyzzyaahh1(:)
real(dp),parameter :: xyzzyaahi1=1.d-8
integer xyzzyaahj1,xyzzyaahk1
real(dp) xyzzyaahl1(3),xyzzyaahm1(3),xyzzyaahn1(3)
logical xyzzyaaho1
integer xyzzyaahp1,xyzzyaahq1
integer,allocatable :: xyzzyaahr1(:,:),xyzzyaahs1(:),xyzzyaaht1(:)
real(dp),allocatable :: xyzzyaahu1(:,:,:),xyzzyaahv1(:)
integer xyzzyaahw1,xyzzyaahx1
integer,allocatable :: xyzzyaahy1(:,:),xyzzyaahz1(:),xyzzyaaia1(:)
real(dp),allocatable :: xyzzyaaib1(:,:,:),xyzzyaaic1(:)
logical xyzzyaaid1
logical,allocatable :: xyzzyaaie1(:),xyzzyaaif1(:)
real(dp),allocatable :: xyzzyaaig1(:,:,:)
real(dp),allocatable :: xyzzyaaih1(:),xyzzyaaii1(:,:)
integer,allocatable :: xyzzyaaij1(:),xyzzyaaik1(:,:)
real(dp) xyzzyaail1(2),xyzzyaaim1(2),xyzzyaain1(2)
integer xyzzyaaio1(2),xyzzyaaip1(2),xyzzyaaiq1,xyzzyaair1,xyzzyaais1,x&
&yzzyaait1
logical xyzzyaaiu1
real(dp),allocatable :: xyzzyaaiv1(:,:),xyzzyaaiw1(:)
integer xyzzyaaix1
integer,parameter :: xyzzyaaiy1=9999
logical,parameter :: xyzzyaaiz1=.false.
integer xyzzyaaja1,xyzzyaajb1,xyzzyaajc1
real(dp),allocatable :: xyzzyaajd1(:),xyzzyaaje1(:,:,:),xyzzyaajf1(:,:&
&),xyzzyaajg1(:,:),xyzzyaajh1(:,:,:),xyzzyaaji1(:,:),xyzzyaajj1(:,:),x&
&yzzyaajk1(:,:),xyzzyaajl1(:,:,:),xyzzyaajm1(:,:),xyzzyaajn1(:,:),xyzz&
&yaajo1(:,:,:),xyzzyaajp1(:,:),xyzzyaajq1(:,:),xyzzyaajr1(:,:,:),xyzzy&
&aajs1(:,:),xyzzyaajt1(:,:),xyzzyaaju1(:,:,:),xyzzyaajv1(:,:),xyzzyaaj&
&w1(:,:),xyzzyaajx1(:,:,:),xyzzyaajy1(:,:),xyzzyaajz1(:,:),xyzzyaaka1(&
&:,:,:),xyzzyaakb1(:,:),xyzzyaakc1(:,:),xyzzyaakd1(:,:,:),xyzzyaake1(:&
&,:),xyzzyaakf1(:,:),xyzzyaakg1(:,:,:),xyzzyaakh1(:,:),xyzzyaaki1(:),x&
&yzzyaakj1(:,:,:),xyzzyaakk1(:,:),xyzzyaakl1(:,:),xyzzyaakm1(:,:,:),xy&
&zzyaakn1(:,:)
logical,allocatable :: xyzzyaako1(:,:),xyzzyaakp1(:),xyzzyaakq1(:,:),x&
&yzzyaakr1(:,:),xyzzyaaks1(:,:),xyzzyaakt1(:,:),xyzzyaaku1(:,:),xyzzya&
&akv1(:,:),xyzzyaakw1(:,:),xyzzyaakx1(:,:),xyzzyaaky1(:,:),xyzzyaakz1(&
&:,:),xyzzyaala1(:,:),xyzzyaalb1(:,:),xyzzyaalc1(:,:),xyzzyaald1(:,:),&
&xyzzyaale1(:,:),xyzzyaalf1(:,:),xyzzyaalg1(:,:),xyzzyaalh1(:,:),xyzzy&
&aali1(:,:),xyzzyaalj1(:,:),xyzzyaalk1(:),xyzzyaall1(:,:),xyzzyaalm1(:&
&,:),xyzzyaaln1(:,:)
integer,parameter :: xyzzyaalo1=10
integer xyzzyaalp1,xyzzyaalq1(xyzzyaalo1)
integer,parameter :: xyzzyaalr1=1,xyzzyaals1=2,xyzzyaalt1=3,xyzzyaalu1&
&=4,xyzzyaalv1=5,xyzzyaalw1=6,xyzzyaalx1=7,xyzzyaaly1=8,xyzzyaalz1=9,x&
&yzzyaama1=10
integer,allocatable :: xyzzyaamb1(:)
character(4),parameter :: xyzzyaamc1(1:2)=(/"para","perp"/)
type config_wfn_pjastrow
private
real(dp) pt_valjas
real(dp),pointer :: pt_chi(:),pt_gradchi(:,:)
logical pt_valjas_valid
logical,pointer :: pt_chi_valid(:),pt_gradchi_valid(:)
end type config_wfn_pjastrow
real(dp),allocatable :: xyzzyaamd1(:)
logical xyzzyaame1
real(dp),allocatable :: xyzzyaamf1(:),xyzzyaamg1(:,:,:),xyzzyaamh1(:,:&
&),xyzzyaami1(:,:),xyzzyaamj1(:,:)
logical,allocatable :: xyzzyaamk1(:),xyzzyaaml1(:),xyzzyaamm1(:),xyzzy&
&aamn1(:)
real(dp),allocatable :: xyzzyaamo1(:),xyzzyaamp1(:,:,:),xyzzyaamq1(:,:&
&),xyzzyaamr1(:,:)
integer xyzzyaams1
real(dp),allocatable :: xyzzyaamt1(:)
real(dp),allocatable :: xyzzyaamu1(:),xyzzyaamv1(:),xyzzyaamw1(:),xyzz&
&yaamx1(:),xyzzyaamy1(:,:,:,:)
integer xyzzyaamz1
real(dp),allocatable :: xyzzyaana1(:),xyzzyaanb1(:,:,:),xyzzyaanc1(:),&
&xyzzyaand1(:,:,:,:,:),xyzzyaane1(:),xyzzyaanf1(:)
integer xyzzyaang1,xyzzyaanh1
real(dp),allocatable :: xyzzyaani1(:,:),xyzzyaanj1(:,:,:,:)
integer xyzzyaank1
real(dp),allocatable :: xyzzyaanl1(:,:),xyzzyaanm1(:,:,:,:)
integer xyzzyaann1
real(dp),allocatable :: xyzzyaano1(:,:,:)
integer xyzzyaanp1
real(dp),allocatable :: xyzzyaanq1(:,:,:)
integer xyzzyaanr1
real(dp),allocatable :: xyzzyaans1(:,:),xyzzyaant1(:,:),xyzzyaanu1(:,:&
&),xyzzyaanv1(:,:),xyzzyaanw1(:,:,:)
real(dp),allocatable :: xyzzyaanx1(:,:,:,:),xyzzyaany1(:,:,:,:,:)
real(dp) xyzzyaanz1,xyzzyaaoa1,xyzzyaaob1,xyzzyaaoc1
contains
subroutine query_pjastrow_levels(nlevels)
implicit none
integer,intent(out) :: nlevels
nlevels=1
end subroutine query_pjastrow_levels
subroutine query_pjastrow_level_details(level_score,level_name)
implicit none
integer,intent(out) :: level_score(*)
character(80),intent(out) :: level_name(*)
level_score(1)=1
level_name(1)='Jastrow factor'
end subroutine query_pjastrow_level_details
subroutine check_varmin_linjas_pjastrow
implicit none
if(xyzzyaaeu1.and.xyzzyaaer1)call errstop_master('CHECK_VARMIN_LINJAS_&
&PJASTROW','The varmin_linjas optimization method may be not used for &
&the non-linear Jastrow terms U and W.')
if(xyzzyaaeu1)call errstop_master('CHECK_VARMIN_LINJAS_PJASTROW','The &
&varmin_linjas optimization method cannot used for the three-body W te&
&rm, since the linear parameters in the core functions of W occur quad&
&ratically in W.')
if(xyzzyaaer1)call errstop_master('CHECK_VARMIN_LINJAS_PJASTROW','The &
&varmin_linjas optimization method cannot be used for the non-linear u&
&_RPA term in the Jastrow factor. Consider using an ordinary u term, s&
&ince u_RPA is useless anyway.')
end subroutine check_varmin_linjas_pjastrow
subroutine setup_pjastrow
implicit none
integer xyzzyaaaa5(2),xyzzyaaab5(2),xyzzyaaac5(2),xyzzyaaad5(2),xyzzya&
&aae5(2)
integer xyzzyaaaf5
xyzzyaaaa5=0
xyzzyaaab5=0
xyzzyaaac5=0
xyzzyaaad5=0
xyzzyaaae5=0
call include_range(ratiocfg_from_sz,xyzzyaaaa5)
call include_range(ratiocfg_to_sz,xyzzyaaaa5)
call include_range(ratio1_from_sz,xyzzyaaaa5)
call include_range(ratio1_to_sz,xyzzyaaaa5)
call include_range(ratio2_from_sz,xyzzyaaaa5)
call include_range(ratio2_to_sz,xyzzyaaaa5)
call include_range(ratio_ion_from_sz,xyzzyaaaa5)
call include_range(ratio_ion_to_sz,xyzzyaaaa5)
call include_range(drift_sz,xyzzyaaab5)
call include_range(kinetic_sz,xyzzyaaab5)
call include_range(kinetic_sz,xyzzyaaac5)
call include_range(wfn_detail_sz,xyzzyaaaa5)
call include_range(wfn_detail_sz,xyzzyaaad5)
call include_range(kinetic_detail_sz,xyzzyaaaa5)
call include_range(kinetic_detail_sz,xyzzyaaab5)
call include_range(kinetic_detail_sz,xyzzyaaac5)
call include_range(kinetic_detail_sz,xyzzyaaae5)
if(use_altsamp)then
call include_range(ratio1_to_sz,xyzzyaaab5)
endif
if(xyzzyaaaa5(1)/=0)then
allocate(xyzzyaajd1(xyzzyaaaa5(1):xyzzyaaaa5(2)),xyzzyaajg1(netot,xyzz&
&yaaaa5(1):xyzzyaaaa5(2)),xyzzyaaji1(netot,xyzzyaaaa5(1):xyzzyaaaa5(2)&
&),xyzzyaajj1(netot,xyzzyaaaa5(1):xyzzyaaaa5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas')
xyzzyaajd1=0.d0
xyzzyaajg1=0.d0
xyzzyaaji1=0.d0
xyzzyaajj1=0.d0
endif
allocate(xyzzyaakp1(nscratch),xyzzyaaks1(netot,nscratch),xyzzyaako1(ne&
&tot,nscratch),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas_valid')
if(xyzzyaaab5(1)/=0)then
allocate(xyzzyaaje1(3,netot,xyzzyaaab5(1):xyzzyaaab5(2)),xyzzyaajh1(3,&
&netot,xyzzyaaab5(1):xyzzyaaab5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','djas')
xyzzyaaje1=0.d0
xyzzyaajh1=0.d0
endif
allocate(xyzzyaakq1(netot,nscratch),xyzzyaakt1(netot,nscratch),stat=xy&
&zzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','djas_valid')
if(xyzzyaaac5(1)/=0)then
allocate(xyzzyaajf1(netot,xyzzyaaac5(1):xyzzyaaac5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','d2jas')
xyzzyaajf1=0.d0
endif
allocate(xyzzyaakr1(netot,nscratch),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','d2jas_valid')
if(xyzzyaaad5(1)/=0)then
if(xyzzyaaep1)then
allocate(xyzzyaajk1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1u')
xyzzyaajk1=0.d0
endif
if(xyzzyaaes1)then
allocate(xyzzyaajn1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1ucyl')
xyzzyaajn1=0.d0
endif
if(xyzzyaaet1)then
allocate(xyzzyaajq1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1qcusp')
xyzzyaajq1=0.d0
endif
if(xyzzyaafe1)then
allocate(xyzzyaajt1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1w')
xyzzyaajt1=0.d0
endif
if(xyzzyaaew1)then
allocate(xyzzyaajw1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1chi')
xyzzyaajw1=0.d0
endif
if(xyzzyaaex1)then
allocate(xyzzyaajz1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1f')
xyzzyaajz1=0.d0
endif
if(xyzzyaaey1)then
allocate(xyzzyaakc1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1p')
xyzzyaakc1=0.d0
endif
if(xyzzyaaez1)then
allocate(xyzzyaakf1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1q')
xyzzyaakf1=0.d0
endif
if(xyzzyaafd1.or.xyzzyaafg1)then
allocate(xyzzyaaki1(xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1b')
xyzzyaaki1=0.d0
endif
if(xyzzyaaff1)then
allocate(xyzzyaakl1(netot,xyzzyaaad5(1):xyzzyaaad5(2)),stat=xyzzyaaaf5&
&)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1d')
xyzzyaakl1=0.d0
endif
endif
allocate(xyzzyaaku1(netot,nscratch),xyzzyaakw1(netot,nscratch),xyzzyaa&
&ky1(netot,nscratch),xyzzyaala1(netot,nscratch),xyzzyaalc1(netot,nscra&
&tch),xyzzyaale1(netot,nscratch),xyzzyaalg1(netot,nscratch),xyzzyaali1&
&(netot,nscratch),xyzzyaalk1(nscratch),xyzzyaalm1(netot,nscratch),stat&
&=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','jas1*_valid')
if(xyzzyaaae5(1)/=0)then
if(xyzzyaaep1)then
allocate(xyzzyaajl1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaajm1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1u')
endif
if(xyzzyaaes1)then
allocate(xyzzyaajo1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaajp1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1ucyl')
endif
if(xyzzyaaet1)then
allocate(xyzzyaajr1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaajs1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1qcusp')
endif
if(xyzzyaafe1)then
allocate(xyzzyaaju1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaajv1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1w')
endif
if(xyzzyaaew1)then
allocate(xyzzyaajx1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaajy1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1chi')
endif
if(xyzzyaaex1)then
allocate(xyzzyaaka1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaakb1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1f')
endif
if(xyzzyaaey1)then
allocate(xyzzyaakd1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaake1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1p')
endif
if(xyzzyaaez1)then
allocate(xyzzyaakg1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaakh1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1q')
endif
if(xyzzyaafd1.or.xyzzyaafg1)then
allocate(xyzzyaakj1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaakk1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1b')
endif
if(xyzzyaaff1)then
allocate(xyzzyaakm1(3,netot,xyzzyaaae5(1):xyzzyaaae5(2)),xyzzyaakn1(ne&
&tot,xyzzyaaae5(1):xyzzyaaae5(2)),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1d')
endif
endif
allocate(xyzzyaakv1(netot,nscratch),xyzzyaakx1(netot,nscratch),xyzzyaa&
&kz1(netot,nscratch),xyzzyaalb1(netot,nscratch),xyzzyaald1(netot,nscra&
&tch),xyzzyaalf1(netot,nscratch),xyzzyaalh1(netot,nscratch),xyzzyaalj1&
&(netot,nscratch),xyzzyaall1(netot,nscratch),xyzzyaaln1(netot,nscratch&
&),stat=xyzzyaaaf5)
call check_alloc(xyzzyaaaf5,'SETUP_JASTROW','gljas1*_valid')
xyzzyaaja1=0
xyzzyaajb1=0
xyzzyaajc1=0
end subroutine setup_pjastrow
subroutine finish_pjastrow
implicit none
if(allocated(xyzzyaajd1))deallocate(xyzzyaajd1,xyzzyaajg1,xyzzyaaji1,x&
&yzzyaajj1)
deallocate(xyzzyaakp1,xyzzyaaks1,xyzzyaako1)
if(allocated(xyzzyaaje1))deallocate(xyzzyaaje1,xyzzyaajh1)
deallocate(xyzzyaakq1,xyzzyaakt1)
if(allocated(xyzzyaajf1))deallocate(xyzzyaajf1)
deallocate(xyzzyaakr1)
if(allocated(xyzzyaajk1))deallocate(xyzzyaajk1)
if(allocated(xyzzyaajn1))deallocate(xyzzyaajn1)
if(allocated(xyzzyaajq1))deallocate(xyzzyaajq1)
if(allocated(xyzzyaajt1))deallocate(xyzzyaajt1)
if(allocated(xyzzyaajw1))deallocate(xyzzyaajw1)
if(allocated(xyzzyaajz1))deallocate(xyzzyaajz1)
if(allocated(xyzzyaakc1))deallocate(xyzzyaakc1)
if(allocated(xyzzyaakf1))deallocate(xyzzyaakf1)
if(allocated(xyzzyaaki1))deallocate(xyzzyaaki1)
if(allocated(xyzzyaakl1))deallocate(xyzzyaakl1)
deallocate(xyzzyaaku1,xyzzyaakw1,xyzzyaaky1,xyzzyaala1,xyzzyaalc1,xyzz&
&yaale1,xyzzyaalg1,xyzzyaali1,xyzzyaalk1,xyzzyaalm1)
if(allocated(xyzzyaajl1))deallocate(xyzzyaajl1,xyzzyaajm1)
if(allocated(xyzzyaajo1))deallocate(xyzzyaajo1,xyzzyaajp1)
if(allocated(xyzzyaajr1))deallocate(xyzzyaajr1,xyzzyaajs1)
if(allocated(xyzzyaaju1))deallocate(xyzzyaaju1,xyzzyaajv1)
if(allocated(xyzzyaajx1))deallocate(xyzzyaajx1,xyzzyaajy1)
if(allocated(xyzzyaaka1))deallocate(xyzzyaaka1,xyzzyaakb1)
if(allocated(xyzzyaakd1))deallocate(xyzzyaakd1,xyzzyaake1)
if(allocated(xyzzyaakg1))deallocate(xyzzyaakg1,xyzzyaakh1)
if(allocated(xyzzyaakj1))deallocate(xyzzyaakj1,xyzzyaakk1)
if(allocated(xyzzyaakm1))deallocate(xyzzyaakm1,xyzzyaakn1)
deallocate(xyzzyaakv1,xyzzyaakx1,xyzzyaakz1,xyzzyaalb1,xyzzyaald1,xyzz&
&yaalf1,xyzzyaalh1,xyzzyaalj1,xyzzyaall1,xyzzyaaln1)
end subroutine finish_pjastrow
subroutine wfn_ratio_pjastrow(is,js,ilevel,ratio,fd,sd)
implicit none
integer,intent(in) :: is,js,ilevel
real(dp),intent(out) :: ratio
logical,intent(in) :: fd,sd
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7
real(dp) xyzzyaaae7,xyzzyaaaf7,xyzzyaaag7,xyzzyaaah7,xyzzyaaai7,xyzzya&
&aaj7,xyzzyaaak7,xyzzyaaal7,xyzzyaaam7(3),xyzzyaaan7(3)
logical xyzzyaaao7
if(buffer_move1_from(js)==is)then
xyzzyaaaa7=buffer_move1_from_ii(js)
call xyzzyaaod1(xyzzyaaaa7,is,.true.,fd,sd)
call xyzzyaaod1(xyzzyaaaa7,js,.true.,fd,sd,xyzzyaaao7)
if(xyzzyaaao7)then
ratio=0.d0
else
ratio=exp_protect(xyzzyaaji1(xyzzyaaaa7,js)-xyzzyaaji1(xyzzyaaaa7,is))
endif
elseif(buffer_move2_from(js)==is)then
if(hard_sphere)call errstop('WFN_RATIO_PJASTROW','Need to work out wha&
&t this is all about.')
xyzzyaaaa7=buffer_move2_from_ii(js)
xyzzyaaac7=buffer_move2_from_jj(js)
call get_rsele(is)
call get_eevecs(is)
call get_eevecs(js)
if(.not.homogeneous_system)then
call get_eivecs(is)
call get_eivecs(js)
endif
call xyzzyaaod1(xyzzyaaaa7,is,.true.,.false.,.false.)
xyzzyaaae7=xyzzyaaji1(xyzzyaaaa7,is)
xyzzyaaab7=which_spin(xyzzyaaaa7)
xyzzyaaad7=which_spin(xyzzyaaac7)
call xyzzyaapp1(xyzzyaaac7,xyzzyaaad7,rele2jj_chscr(1,js),eevecs_scr(1&
&,1,1,js),eevecs_scr(1,1,xyzzyaaac7,js),eivecs_scr(1,1,1,js),eivecs_sc&
&r(1,1,xyzzyaaac7,js),.true.,.false.,.false.,.true.,.false.,xyzzyaaah7&
&,xyzzyaaam7,xyzzyaaai7,xyzzyaaaj7,xyzzyaaan7,xyzzyaaak7,.false.,xyzzy&
&aaal7)
call xyzzyaapp1(xyzzyaaac7,xyzzyaaad7,rele_scr(1,xyzzyaaac7,is),eevecs&
&_scr(1,1,1,js),eevecs_scr(1,1,xyzzyaaac7,is),eivecs_scr(1,1,1,js),eiv&
&ecs_scr(1,1,xyzzyaaac7,is),.true.,.false.,.false.,.true.,.false.,xyzz&
&yaaag7,xyzzyaaam7,xyzzyaaai7,xyzzyaaaj7,xyzzyaaan7,xyzzyaaak7,.false.&
&,xyzzyaaal7)
call xyzzyaapp1(xyzzyaaaa7,xyzzyaaab7,rele2ii_chscr(1,js),eevecs_scr(1&
&,1,1,is),eevecs_scr(1,1,xyzzyaaaa7,js),eivecs_scr(1,1,1,is),eivecs_sc&
&r(1,1,xyzzyaaaa7,js),.true.,.false.,.false.,.true.,.false.,xyzzyaaaf7&
&,xyzzyaaam7,xyzzyaaai7,xyzzyaaaj7,xyzzyaaan7,xyzzyaaak7,.false.,xyzzy&
&aaal7)
ratio=exp_protect(xyzzyaaaf7-xyzzyaaae7+xyzzyaaah7-xyzzyaaag7)
else
call xyzzyaaoe1(is)
call xyzzyaaoe1(js,xyzzyaaao7)
if(xyzzyaaao7)then
ratio=0.d0
else
ratio=exp_protect(xyzzyaajd1(js)-xyzzyaajd1(is))
endif
endif
end subroutine wfn_ratio_pjastrow
subroutine accept_move_pjastrow(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa8,xyzzyaaab8
if(buffer_move1_from(js)==is)then
xyzzyaaaa8=buffer_move1_from_ii(js)
if(xyzzyaaks1(xyzzyaaaa8,js))xyzzyaajg1(xyzzyaaaa8,is)=xyzzyaajg1(xyzz&
&yaaaa8,js)
if(xyzzyaakt1(xyzzyaaaa8,js))xyzzyaajh1(1:3,xyzzyaaaa8,is)=xyzzyaajh1(&
&1:3,xyzzyaaaa8,js)
xyzzyaaks1(xyzzyaaaa8,is)=xyzzyaaks1(xyzzyaaaa8,js)
xyzzyaakt1(xyzzyaaaa8,is)=xyzzyaakt1(xyzzyaaaa8,js)
xyzzyaako1(:,is)=.false.
xyzzyaakq1(:,is)=.false.
xyzzyaakr1(:,is)=.false.
xyzzyaakp1(is)=.false.
if(scr_tasks(iwfn_detail,is).or.scr_tasks(ikinetic_detail,is))then
xyzzyaaku1(:,is)=.false.
xyzzyaakv1(:,is)=.false.
xyzzyaakw1(:,is)=.false.
xyzzyaakx1(:,is)=.false.
xyzzyaaky1(:,is)=.false.
xyzzyaakz1(:,is)=.false.
xyzzyaala1(:,is)=.false.
xyzzyaalb1(:,is)=.false.
xyzzyaalc1(:,is)=.false.
xyzzyaald1(:,is)=.false.
xyzzyaale1(:,is)=.false.
xyzzyaalf1(:,is)=.false.
xyzzyaalg1(:,is)=.false.
xyzzyaalh1(:,is)=.false.
xyzzyaali1(:,is)=.false.
xyzzyaalj1(:,is)=.false.
xyzzyaalk1(is)=.false.
xyzzyaall1(:,is)=.false.
xyzzyaalm1(:,is)=.false.
xyzzyaaln1(:,is)=.false.
endif
elseif(buffer_move2_from(js)==is)then
xyzzyaaaa8=buffer_move2_from_ii(js)
xyzzyaaab8=buffer_move2_from_jj(js)
if(xyzzyaaks1(xyzzyaaaa8,js))xyzzyaajg1(xyzzyaaaa8,is)=xyzzyaajg1(xyzz&
&yaaaa8,js)
if(xyzzyaaks1(xyzzyaaab8,js))xyzzyaajg1(xyzzyaaab8,is)=xyzzyaajg1(xyzz&
&yaaab8,js)
if(xyzzyaakt1(xyzzyaaaa8,js))xyzzyaajh1(1:3,xyzzyaaaa8,is)=xyzzyaajh1(&
&1:3,xyzzyaaaa8,js)
if(xyzzyaakt1(xyzzyaaab8,js))xyzzyaajh1(1:3,xyzzyaaab8,is)=xyzzyaajh1(&
&1:3,xyzzyaaab8,js)
xyzzyaaks1(xyzzyaaaa8,is)=xyzzyaaks1(xyzzyaaaa8,js)
xyzzyaaks1(xyzzyaaab8,is)=xyzzyaaks1(xyzzyaaab8,js)
xyzzyaakt1(xyzzyaaaa8,is)=xyzzyaakt1(xyzzyaaaa8,js)
xyzzyaakt1(xyzzyaaab8,is)=xyzzyaakt1(xyzzyaaab8,js)
xyzzyaako1(:,is)=.false.
xyzzyaakq1(:,is)=.false.
xyzzyaakr1(:,is)=.false.
xyzzyaakp1(is)=.false.
if(scr_tasks(iwfn_detail,is).or.scr_tasks(ikinetic_detail,is))then
xyzzyaaku1(:,is)=.false.
xyzzyaakv1(:,is)=.false.
xyzzyaakw1(:,is)=.false.
xyzzyaakx1(:,is)=.false.
xyzzyaaky1(:,is)=.false.
xyzzyaakz1(:,is)=.false.
xyzzyaala1(:,is)=.false.
xyzzyaalb1(:,is)=.false.
xyzzyaalc1(:,is)=.false.
xyzzyaald1(:,is)=.false.
xyzzyaale1(:,is)=.false.
xyzzyaalf1(:,is)=.false.
xyzzyaalg1(:,is)=.false.
xyzzyaalh1(:,is)=.false.
xyzzyaali1(:,is)=.false.
xyzzyaalj1(:,is)=.false.
xyzzyaalk1(is)=.false.
xyzzyaall1(:,is)=.false.
xyzzyaalm1(:,is)=.false.
xyzzyaaln1(:,is)=.false.
endif
else
do xyzzyaaaa8=1,netot
if(xyzzyaako1(xyzzyaaaa8,js))then
xyzzyaaji1(xyzzyaaaa8,is)=xyzzyaaji1(xyzzyaaaa8,js)
xyzzyaajj1(xyzzyaaaa8,is)=xyzzyaajj1(xyzzyaaaa8,js)
endif
if(xyzzyaakq1(xyzzyaaaa8,js))xyzzyaaje1(1:3,xyzzyaaaa8,is)=xyzzyaaje1(&
&1:3,xyzzyaaaa8,js)
if(xyzzyaakr1(xyzzyaaaa8,js))xyzzyaajf1(xyzzyaaaa8,is)=xyzzyaajf1(xyzz&
&yaaaa8,js)
if(xyzzyaaks1(xyzzyaaaa8,js))xyzzyaajg1(xyzzyaaaa8,is)=xyzzyaajg1(xyzz&
&yaaaa8,js)
if(xyzzyaakt1(xyzzyaaaa8,js))xyzzyaajh1(1:3,xyzzyaaaa8,is)=xyzzyaajh1(&
&1:3,xyzzyaaaa8,js)
xyzzyaako1(xyzzyaaaa8,is)=xyzzyaako1(xyzzyaaaa8,js)
xyzzyaakq1(xyzzyaaaa8,is)=xyzzyaakq1(xyzzyaaaa8,js)
xyzzyaakr1(xyzzyaaaa8,is)=xyzzyaakr1(xyzzyaaaa8,js)
xyzzyaaks1(xyzzyaaaa8,is)=xyzzyaaks1(xyzzyaaaa8,js)
xyzzyaakt1(xyzzyaaaa8,is)=xyzzyaakt1(xyzzyaaaa8,js)
enddo
if(scr_tasks(iwfn_detail,is).or.scr_tasks(ikinetic_detail,is))then
xyzzyaaku1(:,is)=.false.
xyzzyaakv1(:,is)=.false.
xyzzyaakw1(:,is)=.false.
xyzzyaakx1(:,is)=.false.
xyzzyaaky1(:,is)=.false.
xyzzyaakz1(:,is)=.false.
xyzzyaala1(:,is)=.false.
xyzzyaalb1(:,is)=.false.
xyzzyaalc1(:,is)=.false.
xyzzyaald1(:,is)=.false.
xyzzyaale1(:,is)=.false.
xyzzyaalf1(:,is)=.false.
xyzzyaalg1(:,is)=.false.
xyzzyaalh1(:,is)=.false.
xyzzyaali1(:,is)=.false.
xyzzyaalj1(:,is)=.false.
xyzzyaalk1(is)=.false.
xyzzyaall1(:,is)=.false.
xyzzyaalm1(:,is)=.false.
xyzzyaaln1(:,is)=.false.
endif
endif
if(xyzzyaakp1(js))then
xyzzyaajd1(is)=xyzzyaajd1(js)
xyzzyaakp1(is)=.true.
else
xyzzyaakp1(is)=.false.
endif
if(xyzzyaaja1==is)then
xyzzyaaja1=0
elseif(xyzzyaaja1==js)then
xyzzyaaja1=is
endif
if(xyzzyaajb1==is)then
xyzzyaajb1=0
elseif(xyzzyaajb1==js)then
xyzzyaajb1=is
endif
if(xyzzyaajc1==is)then
xyzzyaajc1=0
elseif(xyzzyaajc1==js)then
xyzzyaajc1=is
endif
end subroutine accept_move_pjastrow
subroutine reset_config_pjastrow(is,js)
implicit none
integer,intent(in) :: is,js
call clear_scratch_pjastrow(js)
end subroutine reset_config_pjastrow
subroutine wfn_logval_pjastrow(is,logwfn)
implicit none
integer,intent(in) :: is
real(dp),intent(out) :: logwfn
call xyzzyaaoe1(is)
logwfn=xyzzyaajd1(is)
end subroutine wfn_logval_pjastrow
subroutine wfn_loggrad_pjastrow(ii,is,ilevel,val,sd,gradjas)
implicit none
integer,intent(in) :: ii,is,ilevel
real(dp),intent(out) :: gradjas(3)
logical,intent(in) :: val,sd
call xyzzyaaod1(ii,is,val,.true.,sd)
gradjas=xyzzyaaje1(1:3,ii,is)
end subroutine wfn_loggrad_pjastrow
subroutine wfn_loglap_pjastrow(ii,is,ilevel,val,fd,lapjas)
implicit none
integer,intent(in) :: ii,is,ilevel
real(dp),intent(out) :: lapjas
logical,intent(in) :: val,fd
call xyzzyaaod1(ii,is,val,fd,.true.)
lapjas=xyzzyaajf1(ii,is)
end subroutine wfn_loglap_pjastrow
subroutine prefetch_wfn_pjastrow(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
call get_eevecs(is)
call get_eivecs(is)
if(have_jastrow3)then
call xyzzyaapo1(eevecs_scr(1,1,1,is),fd.or.sd,sd)
xyzzyaaja1=is
if(fd)xyzzyaajb1=is
if(sd)xyzzyaajc1=is
endif
end subroutine prefetch_wfn_pjastrow
subroutine xyzzyaaod1(ii,is,val,fd,sd,minus_infty)
implicit none
integer,intent(in) :: ii,is
logical,intent(in) :: val,fd,sd
logical,intent(out),optional :: minus_infty
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14
real(dp) xyzzyaaad14,xyzzyaaae14,xyzzyaaaf14,xyzzyaaag14(3),xyzzyaaah1&
&4,xyzzyaaai14(3),xyzzyaaaj14,xyzzyaaak14,xyzzyaaal14(3),xyzzyaaam14,x&
&yzzyaaan14,xyzzyaaao14(3),xyzzyaaap14,xyzzyaaaq14,xyzzyaaar14(3),xyzz&
&yaaas14,xyzzyaaat14,xyzzyaaau14(3),xyzzyaaav14,xyzzyaaaw14,xyzzyaaax1&
&4(3),xyzzyaaay14,xyzzyaaaz14,xyzzyaaba14(3),xyzzyaabb14,xyzzyaabc14,x&
&yzzyaabd14(3),xyzzyaabe14,xyzzyaabf14,xyzzyaabg14(3),xyzzyaabh14,xyzz&
&yaabi14,xyzzyaabj14(3),xyzzyaabk14,xyzzyaabl14,xyzzyaabm14(3),xyzzyaa&
&bn14
real(dp),pointer :: xyzzyaabo14(:,:,:),xyzzyaabp14(:,:,:),xyzzyaabq14(&
&:,:),xyzzyaabr14(:,:),xyzzyaabs14(:)
logical need_val,xyzzyaabt14,xyzzyaabu14,xyzzyaabv14,xyzzyaabw14,xyzzy&
&aabx14,xyzzyaaby14,xyzzyaabz14,xyzzyaaca14,xyzzyaacb14,xyzzyaacc14,xy&
&zzyaacd14,xyzzyaace14,xyzzyaacf14,xyzzyaacg14,xyzzyaach14
need_val=val.and..not.xyzzyaako1(ii,is)
xyzzyaabt14=fd.and..not.xyzzyaakq1(ii,is)
xyzzyaabu14=sd.and..not.xyzzyaakr1(ii,is)
xyzzyaabv14=need_val.or.xyzzyaabt14.or.xyzzyaabu14
if(.not.xyzzyaabv14)return
if(ii==buffer_move1_from_ii(is))then
xyzzyaabs14=>rele1_chscr(:,is)
else
call get_rsele(is)
xyzzyaabs14=>rele_scr(:,ii,is)
endif
xyzzyaabw14=.false.
if(have_jastrow3)then
xyzzyaabw14=xyzzyaaja1==is
if(xyzzyaabt14.and.xyzzyaabw14)xyzzyaabw14=xyzzyaajb1==is
if(xyzzyaabu14.and.xyzzyaabw14)xyzzyaabw14=xyzzyaajc1==is
call get_eevecs(is)
endif
if(eevecs_valid(is))then
xyzzyaabq14=>eevecs_scr(:,:,ii,is)
xyzzyaabo14=>eevecs_scr(:,:,:,is)
else
xyzzyaabo14=>eevecs_scr(:,:,:,is)
xyzzyaaab14=buffer_move1_from(is)
xyzzyaach14=xyzzyaaab14/=0.and.ii==buffer_move1_from_ii(is)
if(xyzzyaach14)then
call get_eevecs1_ch(ii,is)
xyzzyaabq14=>eevecs1_chscr(:,:,is)
else
call get_eevecs(is)
xyzzyaabq14=>eevecs_scr(:,:,ii,is)
endif
endif
if(present(minus_infty))then
minus_infty=.false.
if(hard_sphere)then
xyzzyaaaa14=which_spin(ii)
do xyzzyaaac14=1,netot
if(xyzzyaaac14==ii)cycle
if(.not.hard_op_spins.or.which_spin(xyzzyaaac14)/=xyzzyaaaa14)then
if(xyzzyaabq14(4,xyzzyaaac14)<=hard_diam)then
minus_infty=.true.
return
endif
endif
enddo
endif
endif
if(.not.homogeneous_system)then
xyzzyaaab14=buffer_move1_from(is)
xyzzyaach14=xyzzyaaab14/=0.and.ii==buffer_move1_from_ii(is).and..not.e&
&ivecs_valid(is)
if(xyzzyaach14)then
call get_eivecs(xyzzyaaab14)
xyzzyaabp14=>eivecs_scr(:,:,:,xyzzyaaab14)
call get_eivecs1_ch(ii,is)
xyzzyaabr14=>eivecs1_chscr(:,:,is)
else
call get_eivecs(is)
xyzzyaabp14=>eivecs_scr(:,:,:,is)
xyzzyaabr14=>eivecs_scr(:,:,ii,is)
endif
else
xyzzyaabp14=>eivecs_scr(:,:,:,is)
xyzzyaabr14=>eivecs1
endif
xyzzyaaaa14=which_spin(ii)
if(.not.(xyzzyaabt14.and.xyzzyaabu14.and.scr_tasks(ikinetic_detail,is)&
&))then
call xyzzyaapp1(ii,xyzzyaaaa14,xyzzyaabs14,xyzzyaabo14,xyzzyaabq14,xyz&
&zyaabp14,xyzzyaabr14,need_val,xyzzyaabt14,xyzzyaabu14,.not.xyzzyaaks1&
&(ii,is),.not.xyzzyaakt1(ii,is),xyzzyaaad14,xyzzyaaag14,xyzzyaaah14,xy&
&zzyaaae14,xyzzyaaai14,xyzzyaaaf14,xyzzyaabw14,xyzzyaaaj14)
if(need_val)then
if(xyzzyaaks1(ii,is))then
xyzzyaaae14=xyzzyaajg1(ii,is)
xyzzyaaad14=xyzzyaaad14+xyzzyaaae14
else
xyzzyaajg1(ii,is)=xyzzyaaae14
xyzzyaaks1(ii,is)=.true.
endif
xyzzyaaji1(ii,is)=xyzzyaaad14
xyzzyaajj1(ii,is)=0.5d0*(xyzzyaaad14+xyzzyaaae14-third*xyzzyaaaf14-xyz&
&zyaaaj14)+xyzzyaaaj14*inv_netot
xyzzyaako1(ii,is)=.true.
endif
if(xyzzyaabt14)then
if(xyzzyaakt1(ii,is))then
xyzzyaaag14=xyzzyaaag14+xyzzyaajh1(:,ii,is)
else
xyzzyaajh1(:,ii,is)=xyzzyaaai14
xyzzyaakt1(ii,is)=.true.
endif
xyzzyaaje1(:,ii,is)=xyzzyaaag14
xyzzyaakq1(ii,is)=.true.
endif
if(xyzzyaabu14)then
xyzzyaajf1(ii,is)=xyzzyaaah14
xyzzyaakr1(ii,is)=.true.
endif
else
if(need_val.and..not.scr_tasks(iwfn_detail,is))call errstop('GET_JASTR&
&OW1','Attempted to get Jastrow value using detailed version without h&
&aving detail arrays for the wave function. This is easily fixed in at&
& least two ways. Just ask. P.')
if(need_val)then
xyzzyaabx14=xyzzyaaep1.and..not.(xyzzyaaku1(ii,is).and.xyzzyaakv1(ii,i&
&s))
xyzzyaaby14=xyzzyaaes1.and..not.(xyzzyaakw1(ii,is).and.xyzzyaakx1(ii,i&
&s))
xyzzyaabz14=xyzzyaaet1.and..not.(xyzzyaaky1(ii,is).and.xyzzyaakz1(ii,i&
&s))
xyzzyaaca14=xyzzyaafe1.and..not.(xyzzyaala1(ii,is).and.xyzzyaalb1(ii,i&
&s))
xyzzyaacb14=xyzzyaaew1.and..not.(xyzzyaalc1(ii,is).and.xyzzyaald1(ii,i&
&s))
xyzzyaacc14=xyzzyaaex1.and..not.(xyzzyaale1(ii,is).and.xyzzyaalf1(ii,i&
&s))
xyzzyaacd14=xyzzyaaey1.and..not.(xyzzyaalg1(ii,is).and.xyzzyaalh1(ii,i&
&s))
xyzzyaace14=xyzzyaaez1.and..not.(xyzzyaali1(ii,is).and.xyzzyaalj1(ii,i&
&s))
xyzzyaacf14=(xyzzyaafd1.or.xyzzyaafg1).and..not.(xyzzyaalk1(is) .and.x&
&yzzyaall1(ii,is))
xyzzyaacg14=xyzzyaaff1.and..not.(xyzzyaalm1(ii,is).and.xyzzyaaln1(ii,i&
&s))
else
xyzzyaabx14=xyzzyaaep1.and..not.xyzzyaakv1(ii,is)
xyzzyaaby14=xyzzyaaes1.and..not.xyzzyaakx1(ii,is)
xyzzyaabz14=xyzzyaaet1.and..not.xyzzyaakz1(ii,is)
xyzzyaaca14=xyzzyaafe1.and..not.xyzzyaalb1(ii,is)
xyzzyaacb14=xyzzyaaew1.and..not.xyzzyaald1(ii,is)
xyzzyaacc14=xyzzyaaex1.and..not.xyzzyaalf1(ii,is)
xyzzyaacd14=xyzzyaaey1.and..not.xyzzyaalh1(ii,is)
xyzzyaace14=xyzzyaaez1.and..not.xyzzyaalj1(ii,is)
xyzzyaacf14=(xyzzyaafd1.or.xyzzyaafg1).and..not.xyzzyaall1(ii,is)
xyzzyaacg14=xyzzyaaff1.and..not.xyzzyaaln1(ii,is)
endif
call xyzzyaapq1(ii,xyzzyaaaa14,xyzzyaabs14,xyzzyaabo14,xyzzyaabq14,xyz&
&zyaabp14,need_val,xyzzyaabt14,xyzzyaabu14,xyzzyaabw14,xyzzyaaak14,xyz&
&zyaaal14,xyzzyaaam14,xyzzyaaan14,xyzzyaaao14,xyzzyaaap14,xyzzyaaaw14,&
&xyzzyaaax14,xyzzyaaay14,xyzzyaaaq14,xyzzyaaar14,xyzzyaaas14,xyzzyaaat&
&14,xyzzyaaau14,xyzzyaaav14,xyzzyaaaz14,xyzzyaaba14,xyzzyaabb14,xyzzya&
&abc14,xyzzyaabd14,xyzzyaabe14,xyzzyaabf14,xyzzyaabg14,xyzzyaabh14,xyz&
&zyaabi14,xyzzyaabj14,xyzzyaabk14,xyzzyaabl14,xyzzyaabm14,xyzzyaabn14,&
&xyzzyaabx14,xyzzyaaby14,xyzzyaabz14,xyzzyaaca14,xyzzyaacb14,xyzzyaacc&
&14,xyzzyaacd14,xyzzyaace14,xyzzyaacf14,xyzzyaacg14)
if(need_val)then
if(xyzzyaaep1)then
if(xyzzyaabx14)then
xyzzyaajk1(ii,is)=xyzzyaaak14
xyzzyaajl1(1:3,ii,is)=xyzzyaaal14(1:3)
xyzzyaajm1(ii,is)=xyzzyaaam14
xyzzyaaku1(ii,is)=.true.
xyzzyaakv1(ii,is)=.true.
else
xyzzyaaak14=xyzzyaajk1(ii,is)
xyzzyaaal14(1:3)=xyzzyaajl1(1:3,ii,is)
xyzzyaaam14=xyzzyaajm1(ii,is)
endif
endif
if(xyzzyaaes1)then
if(xyzzyaaby14)then
xyzzyaajn1(ii,is)=xyzzyaaan14
xyzzyaajo1(1:3,ii,is)=xyzzyaaao14(1:3)
xyzzyaajp1(ii,is)=xyzzyaaap14
xyzzyaakw1(ii,is)=.true.
xyzzyaakx1(ii,is)=.true.
else
xyzzyaaan14=xyzzyaajn1(ii,is)
xyzzyaaao14(1:3)=xyzzyaajo1(1:3,ii,is)
xyzzyaaap14=xyzzyaajp1(ii,is)
endif
endif
if(xyzzyaaet1)then
if(xyzzyaabz14)then
xyzzyaajq1(ii,is)=xyzzyaaaw14
xyzzyaajr1(1:3,ii,is)=xyzzyaaax14(1:3)
xyzzyaajs1(ii,is)=xyzzyaaay14
xyzzyaaky1(ii,is)=.true.
xyzzyaakz1(ii,is)=.true.
else
xyzzyaaaw14=xyzzyaajq1(ii,is)
xyzzyaaax14(1:3)=xyzzyaajr1(1:3,ii,is)
xyzzyaaay14=xyzzyaajs1(ii,is)
endif
endif
if(xyzzyaafe1)then
if(xyzzyaaca14)then
xyzzyaajt1(ii,is)=xyzzyaaaq14
xyzzyaaju1(1:3,ii,is)=xyzzyaaar14(1:3)
xyzzyaajv1(ii,is)=xyzzyaaas14
xyzzyaala1(ii,is)=.true.
xyzzyaalb1(ii,is)=.true.
else
xyzzyaaaq14=xyzzyaajt1(ii,is)
xyzzyaaar14(1:3)=xyzzyaaju1(1:3,ii,is)
xyzzyaaas14=xyzzyaajv1(ii,is)
endif
endif
if(xyzzyaaew1)then
if(xyzzyaacb14)then
xyzzyaajw1(ii,is)=xyzzyaaat14
xyzzyaajx1(1:3,ii,is)=xyzzyaaau14(1:3)
xyzzyaajy1(ii,is)=xyzzyaaav14
xyzzyaalc1(ii,is)=.true.
xyzzyaald1(ii,is)=.true.
else
xyzzyaaat14=xyzzyaajw1(ii,is)
xyzzyaaau14(1:3)=xyzzyaajx1(1:3,ii,is)
xyzzyaaav14=xyzzyaajy1(ii,is)
endif
endif
if(xyzzyaaex1)then
if(xyzzyaacc14)then
xyzzyaajz1(ii,is)=xyzzyaaaz14
xyzzyaaka1(1:3,ii,is)=xyzzyaaba14(1:3)
xyzzyaakb1(ii,is)=xyzzyaabb14
xyzzyaale1(ii,is)=.true.
xyzzyaalf1(ii,is)=.true.
else
xyzzyaaaz14=xyzzyaajz1(ii,is)
xyzzyaaba14(1:3)=xyzzyaaka1(1:3,ii,is)
xyzzyaabb14=xyzzyaakb1(ii,is)
endif
endif
if(xyzzyaaey1)then
if(xyzzyaacd14)then
xyzzyaakc1(ii,is)=xyzzyaabc14
xyzzyaakd1(1:3,ii,is)=xyzzyaabd14(1:3)
xyzzyaake1(ii,is)=xyzzyaabe14
xyzzyaalg1(ii,is)=.true.
xyzzyaalh1(ii,is)=.true.
else
xyzzyaabc14=xyzzyaakc1(ii,is)
xyzzyaabd14(1:3)=xyzzyaakd1(1:3,ii,is)
xyzzyaabe14=xyzzyaake1(ii,is)
endif
endif
if(xyzzyaaez1)then
if(xyzzyaace14)then
xyzzyaakf1(ii,is)=xyzzyaabf14
xyzzyaakg1(1:3,ii,is)=xyzzyaabg14(1:3)
xyzzyaakh1(ii,is)=xyzzyaabh14
xyzzyaali1(ii,is)=.true.
xyzzyaalj1(ii,is)=.true.
else
xyzzyaabf14=xyzzyaakf1(ii,is)
xyzzyaabg14(1:3)=xyzzyaakg1(1:3,ii,is)
xyzzyaabh14=xyzzyaakh1(ii,is)
endif
endif
if(xyzzyaafd1.or.xyzzyaafg1)then
if(xyzzyaacf14)then
xyzzyaaki1(is)=xyzzyaabi14
xyzzyaakj1(1:3,ii,is)=xyzzyaabj14(1:3)
xyzzyaakk1(ii,is)=xyzzyaabk14
xyzzyaalk1(is)=.true.
xyzzyaall1(ii,is)=.true.
else
xyzzyaabi14=xyzzyaaki1(is)
xyzzyaabj14(1:3)=xyzzyaakj1(1:3,ii,is)
xyzzyaabk14=xyzzyaakk1(ii,is)
endif
endif
if(xyzzyaaff1)then
if(xyzzyaacg14)then
xyzzyaakl1(ii,is)=xyzzyaabl14
xyzzyaakm1(1:3,ii,is)=xyzzyaabm14(1:3)
xyzzyaakn1(ii,is)=xyzzyaabn14
xyzzyaalm1(ii,is)=.true.
xyzzyaaln1(ii,is)=.true.
else
xyzzyaabl14=xyzzyaakl1(ii,is)
xyzzyaabm14(1:3)=xyzzyaakm1(1:3,ii,is)
xyzzyaabn14=xyzzyaakn1(ii,is)
endif
endif
xyzzyaaji1(ii,is)=xyzzyaaak14+xyzzyaaaw14+xyzzyaaaq14+xyzzyaaat14+xyzz&
&yaaaz14+xyzzyaabc14+xyzzyaabf14+xyzzyaabi14+xyzzyaabl14
xyzzyaajj1(ii,is)=xyzzyaaat14+xyzzyaabf14+0.5d0*(xyzzyaaak14+xyzzyaaaw&
&14+xyzzyaaaz14+xyzzyaabc14+xyzzyaabl14)+third*xyzzyaaaq14+inv_netot*x&
&yzzyaabi14
xyzzyaako1(ii,is)=.true.
else
if(xyzzyaaep1)then
if(xyzzyaabx14)then
xyzzyaajl1(1:3,ii,is)=xyzzyaaal14(1:3)
xyzzyaajm1(ii,is)=xyzzyaaam14
xyzzyaakv1(ii,is)=.true.
else
xyzzyaaal14(1:3)=xyzzyaajl1(1:3,ii,is)
xyzzyaaam14=xyzzyaajm1(ii,is)
endif
endif
if(xyzzyaaes1)then
if(xyzzyaaby14)then
xyzzyaajo1(1:3,ii,is)=xyzzyaaao14(1:3)
xyzzyaajp1(ii,is)=xyzzyaaap14
xyzzyaakx1(ii,is)=.true.
else
xyzzyaaao14(1:3)=xyzzyaajo1(1:3,ii,is)
xyzzyaaap14=xyzzyaajp1(ii,is)
endif
endif
if(xyzzyaaet1)then
if(xyzzyaabz14)then
xyzzyaajr1(1:3,ii,is)=xyzzyaaax14(1:3)
xyzzyaajs1(ii,is)=xyzzyaaay14
xyzzyaakz1(ii,is)=.true.
else
xyzzyaaax14(1:3)=xyzzyaajr1(1:3,ii,is)
xyzzyaaay14=xyzzyaajs1(ii,is)
endif
endif
if(xyzzyaafe1)then
if(xyzzyaaca14)then
xyzzyaaju1(1:3,ii,is)=xyzzyaaar14(1:3)
xyzzyaajv1(ii,is)=xyzzyaaas14
xyzzyaalb1(ii,is)=.true.
else
xyzzyaaar14(1:3)=xyzzyaaju1(1:3,ii,is)
xyzzyaaas14=xyzzyaajv1(ii,is)
endif
endif
if(xyzzyaaew1)then
if(xyzzyaacb14)then
xyzzyaajx1(1:3,ii,is)=xyzzyaaau14(1:3)
xyzzyaajy1(ii,is)=xyzzyaaav14
xyzzyaald1(ii,is)=.true.
else
xyzzyaaau14(1:3)=xyzzyaajx1(1:3,ii,is)
xyzzyaaav14=xyzzyaajy1(ii,is)
endif
endif
if(xyzzyaaex1)then
if(xyzzyaacc14)then
xyzzyaaka1(1:3,ii,is)=xyzzyaaba14(1:3)
xyzzyaakb1(ii,is)=xyzzyaabb14
xyzzyaalf1(ii,is)=.true.
else
xyzzyaaba14(1:3)=xyzzyaaka1(1:3,ii,is)
xyzzyaabb14=xyzzyaakb1(ii,is)
endif
endif
if(xyzzyaaey1)then
if(xyzzyaacd14)then
xyzzyaakd1(1:3,ii,is)=xyzzyaabd14(1:3)
xyzzyaake1(ii,is)=xyzzyaabe14
xyzzyaalh1(ii,is)=.true.
else
xyzzyaabd14(1:3)=xyzzyaakd1(1:3,ii,is)
xyzzyaabe14=xyzzyaake1(ii,is)
endif
endif
if(xyzzyaaez1)then
if(xyzzyaace14)then
xyzzyaakg1(1:3,ii,is)=xyzzyaabg14(1:3)
xyzzyaakh1(ii,is)=xyzzyaabh14
xyzzyaalj1(ii,is)=.true.
else
xyzzyaabg14(1:3)=xyzzyaakg1(1:3,ii,is)
xyzzyaabh14=xyzzyaakh1(ii,is)
endif
endif
if(xyzzyaafd1.or.xyzzyaafg1)then
if(xyzzyaacf14)then
xyzzyaakj1(1:3,ii,is)=xyzzyaabj14(1:3)
xyzzyaakk1(ii,is)=xyzzyaabk14
xyzzyaall1(ii,is)=.true.
else
xyzzyaabj14(1:3)=xyzzyaakj1(1:3,ii,is)
xyzzyaabk14=xyzzyaakk1(ii,is)
endif
endif
if(xyzzyaaff1)then
if(xyzzyaacg14)then
xyzzyaakm1(1:3,ii,is)=xyzzyaabm14(1:3)
xyzzyaakn1(ii,is)=xyzzyaabn14
xyzzyaaln1(ii,is)=.true.
else
xyzzyaabm14(1:3)=xyzzyaakm1(1:3,ii,is)
xyzzyaabn14=xyzzyaakn1(ii,is)
endif
endif
endif
xyzzyaaje1(1:3,ii,is)=xyzzyaaal14(1:3)+xyzzyaaao14(1:3)+xyzzyaaax14(1:&
&3)+xyzzyaaar14(1:3)+xyzzyaaau14(1:3)+xyzzyaaba14(1:3)+xyzzyaabd14(1:3&
&)+xyzzyaabg14(1:3)+xyzzyaabj14(1:3)+xyzzyaabm14(1:3)
xyzzyaakq1(ii,is)=.true.
xyzzyaajf1(ii,is)=xyzzyaaam14+xyzzyaaap14+xyzzyaaay14+xyzzyaaas14+xyzz&
&yaaav14+xyzzyaabb14+xyzzyaabe14+xyzzyaabh14+xyzzyaabk14+xyzzyaabn14
xyzzyaakr1(ii,is)=.true.
endif
end subroutine xyzzyaaod1
subroutine xyzzyaaoe1(is,minus_infty)
implicit none
integer,intent(in) :: is
logical,intent(out),optional :: minus_infty
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15
if(xyzzyaakp1(is))return
if(all(xyzzyaako1(:,is)))then
xyzzyaajd1(is)=sum(xyzzyaajj1(:,is))
else
call get_rsele(is)
call get_eevecs(is)
if(present(minus_infty))then
minus_infty=.false.
if(hard_sphere)then
do xyzzyaaaa15=2,netot
xyzzyaaac15=which_spin(xyzzyaaaa15)
do xyzzyaaab15=1,xyzzyaaaa15-1
if(.not.hard_op_spins.or.xyzzyaaac15/=which_spin(xyzzyaaab15))then
if(eevecs_scr(4,xyzzyaaaa15,xyzzyaaab15,is)<=hard_diam)then
minus_infty=.true.
return
endif
endif
enddo
enddo
endif
endif
call get_eivecs(is)
call xyzzyaapn1(rele_scr(1,1,is),eevecs_scr(1,1,1,is),eivecs_scr(1,1,1&
&,is),xyzzyaajd1(is))
endif
xyzzyaakp1(is)=.true.
end subroutine xyzzyaaoe1
subroutine gen_config_pjastrow(pt_config)
implicit none
type(config_wfn_pjastrow),pointer :: pt_config
integer xyzzyaaaa16
allocate(pt_config,stat=xyzzyaaaa16)
call check_alloc(xyzzyaaaa16,'GEN_CONFIG_WFN_JASTROW','container')
if(jasbuf)then
allocate(pt_config%pt_chi(netot),pt_config%pt_gradchi(3,netot),pt_conf&
&ig%pt_chi_valid(netot),pt_config%pt_gradchi_valid(netot),stat=xyzzyaa&
&aa16)
call check_alloc(xyzzyaaaa16,'GEN_CONFIG_WFN_JASTROW','jasbuf')
pt_config%pt_valjas=0.d0
pt_config%pt_chi(:)=0.d0
pt_config%pt_gradchi(:,:)=0.d0
pt_config%pt_valjas_valid=.false.
pt_config%pt_chi_valid(:)=.false.
pt_config%pt_gradchi_valid(:)=.false.
endif
end subroutine gen_config_pjastrow
subroutine delete_config_pjastrow(pt_config)
implicit none
type(config_wfn_pjastrow),pointer :: pt_config
if(jasbuf)deallocate(pt_config%pt_chi,pt_config%pt_gradchi,pt_config%p&
&t_chi_valid,pt_config%pt_gradchi_valid)
deallocate(pt_config)
end subroutine delete_config_pjastrow
subroutine copy_config_pjastrow(pt_from,pt_to)
implicit none
type(config_wfn_pjastrow),pointer :: pt_from,pt_to
if(jasbuf)then
pt_to%pt_valjas=pt_from%pt_valjas
call dcopy(netot,pt_from%pt_chi,1,pt_to%pt_chi,1)
call dcopy(three_netot,pt_from%pt_gradchi,1,pt_to%pt_gradchi,1)
pt_to%pt_valjas_valid=pt_from%pt_valjas_valid
pt_to%pt_chi_valid=pt_from%pt_chi_valid
pt_to%pt_gradchi_valid=pt_from%pt_gradchi_valid
endif
end subroutine copy_config_pjastrow
subroutine config_to_pt_pjastrow(pt_config,k)
use slaarnaaf, only : valjas_config
implicit none
integer,intent(in) :: k
type(config_wfn_pjastrow),pointer :: pt_config
if(jasbuf)then
pt_config%pt_valjas_valid=allocated(valjas_config)
if(pt_config%pt_valjas_valid)pt_config%pt_valjas=valjas_config(k)
endif
end subroutine config_to_pt_pjastrow
subroutine pt_to_config_pjastrow(pt_config)
use slaarnaaf, only : add_config
implicit none
type(config_wfn_pjastrow),pointer :: pt_config
if(jasbuf)call add_config(modify=.true.,valjas=pt_config%pt_valjas)
end subroutine pt_to_config_pjastrow
subroutine redist_allocations_pjastrow(kmax)
implicit none
integer,intent(in) :: kmax
integer xyzzyaaaa21
xyzzyaame1=jasbuf
if(jasbuf)then
allocate(xyzzyaamd1(kmax),stat=xyzzyaaaa21)
call check_alloc(xyzzyaaaa21,'REDIST_ALLOCATIONS_PJASTROW','valjas_trf&
&')
endif
end subroutine redist_allocations_pjastrow
subroutine redist_load_pjastrow(pt_config,k,blocking)
implicit none
integer,intent(in) :: k
logical,intent(in) :: blocking
type(config_wfn_pjastrow),pointer :: pt_config
xyzzyaame1=xyzzyaame1.and.pt_config%pt_valjas_valid
if(xyzzyaame1.or..not.blocking)xyzzyaamd1(k)=pt_config%pt_valjas
end subroutine redist_load_pjastrow
subroutine redist_send_pjastrow(jnode,kbase,k,nbt,reqbase,blocking)
use parallel
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa23
if(jasbuf)then
if(blocking)then
call mpi_ssend(xyzzyaame1,1,mpi_logical,jnode,move_msg,mpi_comm_world,&
&ierror)
call checkmpi(ierror,'ssend valjas_valid in redist_send')
if(xyzzyaame1)then
call mpi_ssend(xyzzyaamd1(kbase+1),k,mpi_double_precision,jnode,move_m&
&sg,mpi_comm_world,ierror)
call checkmpi(ierror,'ssend valjas in redist_send')
endif
else
nbt=nbt+1
call mpi_isend(xyzzyaamd1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa23,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa23
call checkmpi(ierror,'ssend valjas in redist_send')
endif
endif
end subroutine redist_send_pjastrow
subroutine redist_recv_pjastrow(jnode,kbase,k,nbt,reqbase,blocking)
use parallel
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa24
if(jasbuf)then
if(blocking)then
call mpi_recv(xyzzyaame1,1,mpi_logical,jnode,move_msg,mpi_comm_world,s&
&tatus,ierror)
call checkmpi(ierror,'recv valjas_valid in redist_recv')
if(xyzzyaame1)then
call mpi_recv(xyzzyaamd1(kbase+1),k,mpi_double_precision,jnode,move_ms&
&g,mpi_comm_world,status,ierror)
call checkmpi(ierror,'recv valjas in redist_recv')
endif
else
nbt=nbt+1
call mpi_irecv(xyzzyaamd1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa24,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa24
call checkmpi(ierror,'recv valjas in redist_recv')
endif
endif
end subroutine redist_recv_pjastrow
subroutine redist_save_pjastrow(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_pjastrow),pointer :: pt_config
pt_config%pt_valjas_valid=xyzzyaame1
if(xyzzyaame1)pt_config%pt_valjas=xyzzyaamd1(k)
end subroutine redist_save_pjastrow
subroutine redist_deallocations_pjastrow
implicit none
if(jasbuf)deallocate(xyzzyaamd1)
end subroutine redist_deallocations_pjastrow
subroutine load_from_pt_pjastrow(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_pjastrow),pointer :: pt_config
if(jasbuf)then
xyzzyaajd1(is)=pt_config%pt_valjas
call dcopy(netot,pt_config%pt_chi,1,xyzzyaajg1(1,is),1)
call dcopy(three_netot,pt_config%pt_gradchi,1,xyzzyaajh1(1,1,is),1)
xyzzyaakp1(is)=pt_config%pt_valjas_valid
xyzzyaaks1(1:netot,is)=pt_config%pt_chi_valid(1:netot)
xyzzyaakt1(1:netot,is)=pt_config%pt_gradchi_valid(1:netot)
endif
end subroutine load_from_pt_pjastrow
subroutine save_to_pt_pjastrow(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_pjastrow),pointer :: pt_config
if(jasbuf)then
pt_config%pt_valjas=xyzzyaajd1(is)
call dcopy(netot,xyzzyaajg1(1,is),1,pt_config%pt_chi,1)
call dcopy(three_netot,xyzzyaajh1(1,1,is),1,pt_config%pt_gradchi,1)
pt_config%pt_valjas_valid=xyzzyaakp1(is)
pt_config%pt_chi_valid(1:netot)=xyzzyaaks1(1:netot,is)
pt_config%pt_gradchi_valid(1:netot)=xyzzyaakt1(1:netot,is)
endif
end subroutine save_to_pt_pjastrow
subroutine add_config_pjastrow_items(is)
use slaarnaaf, only : add_config
implicit none
integer,intent(in) :: is
call xyzzyaaoe1(is)
call add_config(modify=.true.,valjas=xyzzyaajd1(is))
end subroutine add_config_pjastrow_items
subroutine setup_storage_pjastrow(nconfig,ignore)
use slaarnaaf
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(xyzzyaalp1)
integer xyzzyaaaa30
allocate(xyzzyaamf1(nconfig),xyzzyaamg1(3,netot,nconfig),xyzzyaamh1(ne&
&tot,nconfig),xyzzyaami1(netot,nconfig),xyzzyaamj1(netot,nconfig),stat&
&=xyzzyaaaa30)
call check_alloc(xyzzyaaaa30,'SETUP_STORAGE_JASTROW','*jas_store')
xyzzyaamf1=0.d0
xyzzyaamg1=0.d0
xyzzyaamh1=0.d0
xyzzyaami1=0.d0
xyzzyaamj1=0.d0
allocate(xyzzyaamk1(nconfig),xyzzyaaml1(nconfig),xyzzyaamm1(nconfig),x&
&yzzyaamn1(nconfig),stat=xyzzyaaaa30)
call check_alloc(xyzzyaaaa30,'SETUP_STORAGE_JASTROW','valjas_svalid')
xyzzyaamk1=.false.
xyzzyaaml1=.false.
xyzzyaamm1=.false.
xyzzyaamn1=.false.
if(allocated(valjas_config))then
call dcopy(nconfig,valjas_config(1),1,xyzzyaamf1(1),1)
xyzzyaamk1(1:nconfig)=.true.
endif
end subroutine setup_storage_pjastrow
subroutine finish_storage_pjastrow
implicit none
deallocate(xyzzyaamf1,xyzzyaamg1,xyzzyaamh1,xyzzyaami1,xyzzyaamj1)
deallocate(xyzzyaamk1,xyzzyaaml1,xyzzyaamm1,xyzzyaamn1)
end subroutine finish_storage_pjastrow
subroutine load_from_storage_pjastrow(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(xyzzyaamk1(icfg))then
xyzzyaajd1(is)=xyzzyaamf1(icfg)
xyzzyaakp1(is)=.true.
endif
if(xyzzyaamn1(icfg))then
call dcopy(netot,xyzzyaami1(1,icfg),1,xyzzyaajj1(1,is),1)
call dcopy(netot,xyzzyaamj1(1,icfg),1,xyzzyaaji1(1,is),1)
xyzzyaako1(:,is)=.true.
endif
if(xyzzyaaml1(icfg))then
call dcopy(three_netot,xyzzyaamg1(1,1,icfg),1,xyzzyaaje1(1,1,is),1)
xyzzyaakq1(:,is)=.true.
endif
if(xyzzyaamm1(icfg))then
call dcopy(netot,xyzzyaamh1(1,icfg),1,xyzzyaajf1(1,is),1)
xyzzyaakr1(:,is)=.true.
endif
end subroutine load_from_storage_pjastrow
subroutine save_to_storage_pjastrow(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(xyzzyaakp1(is).and..not.xyzzyaamk1(icfg))then
xyzzyaamf1(icfg)=xyzzyaajd1(is)
xyzzyaamk1(icfg)=.true.
endif
if(all(xyzzyaako1(:,is)))then
call dcopy(netot,xyzzyaajj1(1,is),1,xyzzyaami1(1,icfg),1)
call dcopy(netot,xyzzyaaji1(1,is),1,xyzzyaamj1(1,icfg),1)
xyzzyaamn1(icfg)=.true.
endif
if(all(xyzzyaakq1(:,is)))then
call dcopy(three_netot,xyzzyaaje1(1,1,is),1,xyzzyaamg1(1,1,icfg),1)
xyzzyaaml1(icfg)=.true.
endif
if(all(xyzzyaakr1(:,is)))then
call dcopy(netot,xyzzyaajf1(1,is),1,xyzzyaamh1(1,icfg),1)
xyzzyaamm1(icfg)=.true.
endif
end subroutine save_to_storage_pjastrow
subroutine clone_scratch_pjastrow(is,js)
implicit none
integer,intent(in) :: is,js
if(scr_tasks(ikinetic,js))then
if(all(xyzzyaakq1(:,is)).and..not.all(xyzzyaakq1(:,js)))then
call dcopy(three_netot,xyzzyaaje1(1,1,is),1,xyzzyaaje1(1,1,js),1)
xyzzyaakq1(:,js)=.true.
endif
if(all(xyzzyaakr1(:,is)).and..not.all(xyzzyaakr1(:,js)))then
call dcopy(netot,xyzzyaajf1(1,is),1,xyzzyaajf1(1,js),1)
xyzzyaakr1(:,js)=.true.
endif
endif
if(scr_tasks(ikinetic_detail,js))then
if(xyzzyaaep1)then
if(all(xyzzyaakv1(:,is)).and..not.all(xyzzyaakv1(:,js)))then
call dcopy(three_netot,xyzzyaajl1(1,1,is),1,xyzzyaajl1(1,1,js),1)
call dcopy(netot,xyzzyaajm1(1,is),1,xyzzyaajm1(1,js),1)
xyzzyaakv1(:,js)=.true.
endif
endif
if(xyzzyaaes1)then
if(all(xyzzyaakx1(:,is)).and..not.all(xyzzyaakx1(:,js)))then
call dcopy(three_netot,xyzzyaajo1(1,1,is),1,xyzzyaajo1(1,1,js),1)
call dcopy(netot,xyzzyaajp1(1,is),1,xyzzyaajp1(1,js),1)
xyzzyaakx1(:,js)=.true.
endif
endif
if(xyzzyaaet1)then
if(all(xyzzyaakz1(:,is)).and..not.all(xyzzyaakz1(:,js)))then
call dcopy(three_netot,xyzzyaajr1(1,1,is),1,xyzzyaajr1(1,1,js),1)
call dcopy(netot,xyzzyaajs1(1,is),1,xyzzyaajs1(1,js),1)
xyzzyaakz1(:,js)=.true.
endif
endif
if(xyzzyaafe1)then
if(all(xyzzyaalb1(:,is)).and..not.all(xyzzyaalb1(:,js)))then
call dcopy(three_netot,xyzzyaaju1(1,1,is),1,xyzzyaaju1(1,1,js),1)
call dcopy(netot,xyzzyaajv1(1,is),1,xyzzyaajv1(1,js),1)
xyzzyaalb1(:,js)=.true.
endif
endif
if(xyzzyaaew1)then
if(all(xyzzyaald1(:,is)).and..not.all(xyzzyaald1(:,js)))then
call dcopy(three_netot,xyzzyaajx1(1,1,is),1,xyzzyaajx1(1,1,js),1)
call dcopy(netot,xyzzyaajy1(1,is),1,xyzzyaajy1(1,js),1)
xyzzyaald1(:,js)=.true.
endif
endif
if(xyzzyaaex1)then
if(all(xyzzyaalf1(:,is)).and..not.all(xyzzyaalf1(:,js)))then
call dcopy(three_netot,xyzzyaaka1(1,1,is),1,xyzzyaaka1(1,1,js),1)
call dcopy(netot,xyzzyaakb1(1,is),1,xyzzyaakb1(1,js),1)
xyzzyaalf1(:,js)=.true.
endif
endif
if(xyzzyaaey1)then
if(all(xyzzyaalh1(:,is)).and..not.all(xyzzyaalh1(:,js)))then
call dcopy(three_netot,xyzzyaakd1(1,1,is),1,xyzzyaakd1(1,1,js),1)
call dcopy(netot,xyzzyaake1(1,is),1,xyzzyaake1(1,js),1)
xyzzyaalh1(:,js)=.true.
endif
endif
if(xyzzyaaez1)then
if(all(xyzzyaalj1(:,is)).and..not.all(xyzzyaalj1(:,js)))then
call dcopy(three_netot,xyzzyaakg1(1,1,is),1,xyzzyaakg1(1,1,js),1)
call dcopy(netot,xyzzyaakh1(1,is),1,xyzzyaakh1(1,js),1)
xyzzyaalj1(:,js)=.true.
endif
endif
if(xyzzyaafd1.or.xyzzyaafg1)then
if(all(xyzzyaall1(:,is)).and..not.all(xyzzyaall1(:,js)))then
call dcopy(three_netot,xyzzyaakj1(1,1,is),1,xyzzyaakj1(1,1,js),1)
call dcopy(netot,xyzzyaakk1(1,is),1,xyzzyaakk1(1,js),1)
xyzzyaall1(:,js)=.true.
endif
endif
if(xyzzyaaff1)then
if(all(xyzzyaaln1(:,is)).and..not.all(xyzzyaaln1(:,js)))then
call dcopy(three_netot,xyzzyaakm1(1,1,is),1,xyzzyaakm1(1,1,js),1)
call dcopy(netot,xyzzyaakn1(1,is),1,xyzzyaakn1(1,js),1)
xyzzyaaln1(:,js)=.true.
endif
endif
endif
if(scr_tasks(iwfn_detail,js))then
if(xyzzyaakp1(is).and..not.xyzzyaakp1(js))then
xyzzyaajd1(js)=xyzzyaajd1(is)
xyzzyaakp1(js)=.true.
endif
if(all(xyzzyaako1(:,is)).and..not.all(xyzzyaako1(:,js)))then
call dcopy(netot,xyzzyaaji1(1,is),1,xyzzyaaji1(1,js),1)
call dcopy(netot,xyzzyaajj1(1,is),1,xyzzyaajj1(1,js),1)
xyzzyaako1(:,js)=.true.
endif
if(all(xyzzyaaks1(:,is)).and..not.all(xyzzyaaks1(:,js)))then
call dcopy(netot,xyzzyaajg1(1,is),1,xyzzyaajg1(1,js),1)
xyzzyaaks1(:,js)=.true.
endif
if(xyzzyaaep1)then
if(all(xyzzyaaku1(:,is)).and..not.all(xyzzyaaku1(:,js)))then
call dcopy(netot,xyzzyaajk1(1,is),1,xyzzyaajk1(1,js),1)
xyzzyaaku1(:,js)=.true.
endif
endif
if(xyzzyaaes1)then
if(all(xyzzyaakw1(:,is)).and..not.all(xyzzyaakw1(:,js)))then
call dcopy(netot,xyzzyaajn1(1,is),1,xyzzyaajn1(1,js),1)
xyzzyaakw1(:,js)=.true.
endif
endif
if(xyzzyaaet1)then
if(all(xyzzyaaky1(:,is)).and..not.all(xyzzyaaky1(:,js)))then
call dcopy(netot,xyzzyaajq1(1,is),1,xyzzyaajq1(1,js),1)
xyzzyaaky1(:,js)=.true.
endif
endif
if(xyzzyaafe1)then
if(all(xyzzyaala1(:,is)).and..not.all(xyzzyaala1(:,js)))then
call dcopy(netot,xyzzyaajt1(1,is),1,xyzzyaajt1(1,js),1)
xyzzyaala1(:,js)=.true.
endif
endif
if(xyzzyaaew1)then
if(all(xyzzyaalc1(:,is)).and..not.all(xyzzyaalc1(:,js)))then
call dcopy(netot,xyzzyaajw1(1,is),1,xyzzyaajw1(1,js),1)
xyzzyaalc1(:,js)=.true.
endif
endif
if(xyzzyaaex1)then
if(all(xyzzyaale1(:,is)).and..not.all(xyzzyaale1(:,js)))then
call dcopy(netot,xyzzyaajz1(1,is),1,xyzzyaajz1(1,js),1)
xyzzyaale1(:,js)=.true.
endif
endif
if(xyzzyaaey1)then
if(all(xyzzyaalg1(:,is)).and..not.all(xyzzyaalg1(:,js)))then
call dcopy(netot,xyzzyaakc1(1,is),1,xyzzyaakc1(1,js),1)
xyzzyaalg1(:,js)=.true.
endif
endif
if(xyzzyaaez1)then
if(all(xyzzyaali1(:,is)).and..not.all(xyzzyaali1(:,js)))then
call dcopy(netot,xyzzyaakf1(1,is),1,xyzzyaakf1(1,js),1)
xyzzyaali1(:,js)=.true.
endif
endif
if(xyzzyaafd1.or.xyzzyaafg1)then
if(xyzzyaalk1(is).and..not.xyzzyaalk1(js))then
xyzzyaaki1(js)=xyzzyaaki1(is)
xyzzyaalk1(js)=.true.
endif
endif
if(xyzzyaaff1)then
if(all(xyzzyaalm1(:,is)).and..not.all(xyzzyaalm1(:,js)))then
call dcopy(netot,xyzzyaakl1(1,is),1,xyzzyaakl1(1,js),1)
xyzzyaalm1(:,js)=.true.
endif
endif
endif
end subroutine clone_scratch_pjastrow
subroutine setup_pjastrow_params(nparam)
implicit none
integer,intent(inout) :: nparam
integer xyzzyaaaa35,xyzzyaaab35,xyzzyaaac35,xyzzyaaad35,xyzzyaaae35,xy&
&zzyaaaf35,xyzzyaaag35,xyzzyaaah35,xyzzyaaai35,xyzzyaaaj35,xyzzyaaak35
nparam=0
xyzzyaalq1=0
if(.not.opt_jastrow)return
if(xyzzyaaep1)then
xyzzyaaaa35=0
if(xyzzyaacx1==1)xyzzyaaaa35=xyzzyaaaa35+1
if(xyzzyaaeq1)then
do xyzzyaaae35=1,xyzzyaafp1
if(xyzzyaaat1(0,xyzzyaaae35)==1)xyzzyaaaa35=xyzzyaaaa35+1
do xyzzyaaaf35=2,xyzzyaabx1
if(xyzzyaaat1(xyzzyaaaf35,xyzzyaaae35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
enddo
else
do xyzzyaaae35=1,xyzzyaafp1
if(xyzzyaaba1(xyzzyaaae35)==1)xyzzyaaaa35=xyzzyaaaa35+1
if(xyzzyaaaz1(xyzzyaaae35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
endif
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalr1)=xyzzyaaaa35
endif
if(xyzzyaaes1)then
xyzzyaaaa35=0
if(xyzzyaadd1==1)xyzzyaaaa35=xyzzyaaaa35+1
if(xyzzyaade1==1)xyzzyaaaa35=xyzzyaaaa35+1
if(xyzzyaacy1==1)xyzzyaaaa35=xyzzyaaaa35+1
if(xyzzyaacz1==1)xyzzyaaaa35=xyzzyaaaa35+1
do xyzzyaaae35=1,xyzzyaafq1
do xyzzyaaah35=0,xyzzyaabz1
do xyzzyaaaf35=0,xyzzyaaby1
if(xyzzyaaau1(xyzzyaaaf35,xyzzyaaah35,xyzzyaaae35)==1)xyzzyaaaa35=xyzz&
&yaaaa35+1
enddo
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaals1)=xyzzyaaaa35
endif
if(xyzzyaaet1)then
xyzzyaaaa35=0
if(xyzzyaadc1==1)xyzzyaaaa35=xyzzyaaaa35+1
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalt1)=xyzzyaaaa35
endif
if(xyzzyaaeu1)then
xyzzyaaaa35=0
if(xyzzyaada1==1)xyzzyaaaa35=xyzzyaaaa35+1
do xyzzyaaae35=1,xyzzyaaft1
do xyzzyaaaf35=0,xyzzyaacb1
if(xyzzyaabb1(xyzzyaaaf35,xyzzyaaae35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalu1)=xyzzyaaaa35
elseif(xyzzyaaev1)then
xyzzyaaaa35=0
if(xyzzyaadb1==1)xyzzyaaaa35=xyzzyaaaa35+1
do xyzzyaaae35=1,xyzzyaafu1
do xyzzyaaaf35=0,xyzzyaacc1
do xyzzyaaah35=0,xyzzyaacc1
do xyzzyaaai35=0,xyzzyaacc1
if(xyzzyaabc1(xyzzyaaaf35,xyzzyaaah35,xyzzyaaai35,xyzzyaaae35)==1)xyzz&
&yaaaa35=xyzzyaaaa35+1
enddo
enddo
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalu1)=xyzzyaaaa35
endif
if(xyzzyaaew1)then
xyzzyaaaa35=0
do xyzzyaaag35=1,xyzzyaabi1
if(xyzzyaadf1(xyzzyaaag35)==1)xyzzyaaaa35=xyzzyaaaa35+1
do xyzzyaaae35=1,xyzzyaafv1(xyzzyaaag35)
if(xyzzyaaav1(0,xyzzyaaae35,xyzzyaaag35)==1)xyzzyaaaa35=xyzzyaaaa35+1
do xyzzyaaah35=2,xyzzyaace1(xyzzyaaag35)
if(xyzzyaaav1(xyzzyaaah35,xyzzyaaae35,xyzzyaaag35)==1)xyzzyaaaa35=xyzz&
&yaaaa35+1
enddo
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalv1)=xyzzyaaaa35
endif
if(xyzzyaaex1)then
xyzzyaaaa35=0
do xyzzyaaag35=1,xyzzyaabj1
if(xyzzyaadg1(xyzzyaaag35)==1)xyzzyaaaa35=xyzzyaaaa35+1
do xyzzyaaae35=1,xyzzyaafw1(xyzzyaaag35)
do xyzzyaaai35=0,xyzzyaacf1(xyzzyaaag35)
do xyzzyaaah35=0,xyzzyaacg1(xyzzyaaag35)
do xyzzyaaaf35=xyzzyaaah35,xyzzyaacg1(xyzzyaaag35)
if(xyzzyaaaw1(xyzzyaaaf35,xyzzyaaah35,xyzzyaaai35,xyzzyaaae35,xyzzyaaa&
&g35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
enddo
enddo
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalw1)=xyzzyaaaa35
endif
if(xyzzyaaey1)then
xyzzyaaaa35=0
do xyzzyaaae35=1,xyzzyaafr1
do xyzzyaaaj35=1,xyzzyaagq1
if(xyzzyaaax1(xyzzyaaaj35,xyzzyaaae35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalx1)=xyzzyaaaa35
endif
if(xyzzyaaez1)then
xyzzyaaaa35=0
do xyzzyaaae35=1,xyzzyaafs1
do xyzzyaaaj35=1,xyzzyaagr1
if(xyzzyaaay1(xyzzyaaaj35,xyzzyaaae35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaaly1)=xyzzyaaaa35
endif
if(xyzzyaafa1.or.xyzzyaafc1)then
xyzzyaaaa35=0
do xyzzyaaaj35=1,9
if(xyzzyaaij1(xyzzyaaaj35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalz1)=xyzzyaaaa35
elseif(xyzzyaafb1)then
xyzzyaaaa35=0
do xyzzyaaaj35=1,2
if(xyzzyaaio1(xyzzyaaaj35)==1)xyzzyaaaa35=xyzzyaaaa35+1
if(xyzzyaaip1(xyzzyaaaj35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalz1)=xyzzyaaaa35
elseif(xyzzyaafg1)then
xyzzyaaaa35=0
do xyzzyaaae35=1,xyzzyaair1
do xyzzyaaaj35=1,xyzzyaaiq1
if(xyzzyaaik1(xyzzyaaaj35,xyzzyaaae35)==1)xyzzyaaaa35=xyzzyaaaa35+1
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaalz1)=xyzzyaaaa35
endif
if(xyzzyaaff1)then
xyzzyaaaa35=0
do xyzzyaaag35=1,xyzzyaabn1
do xyzzyaaak35=1,2
if(xyzzyaadh1(xyzzyaaag35,xyzzyaaak35,1)==1)xyzzyaaaa35=xyzzyaaaa35+1
if(xyzzyaabw1(xyzzyaaag35)==0.and.xyzzyaadh1(xyzzyaaag35,xyzzyaaak35,2&
&)==1)xyzzyaaaa35=xyzzyaaaa35+1
do xyzzyaaah35=0,xyzzyaacd1
do xyzzyaaai35=xyzzyaaah35*xyzzyaabw1(xyzzyaaag35),xyzzyaacd1
if(xyzzyaabe1(xyzzyaaag35,xyzzyaaak35,xyzzyaaah35,xyzzyaaai35)==1)xyzz&
&yaaaa35=xyzzyaaaa35+1
enddo
enddo
enddo
enddo
if(xyzzyaaaa35==0)xyzzyaaaa35=-1
xyzzyaalq1(xyzzyaama1)=xyzzyaaaa35
endif
xyzzyaalp1=sum(xyzzyaalq1,xyzzyaalq1>0)
nparam=xyzzyaalp1
allocate(xyzzyaamb1(xyzzyaalp1),stat=xyzzyaaad35)
call check_alloc(xyzzyaaad35,'SETUP_JASTROW_PARAMS','jastrow_param_sec&
&')
xyzzyaaac35=0
do xyzzyaaaa35=1,xyzzyaalo1
if(xyzzyaalq1(xyzzyaaaa35)<1)cycle
xyzzyaaab35=xyzzyaaac35+1
xyzzyaaac35=xyzzyaaac35+xyzzyaalq1(xyzzyaaaa35)
xyzzyaamb1(xyzzyaaab35:xyzzyaaac35)=xyzzyaaaa35
enddo
call xyzzyaaof1
end subroutine setup_pjastrow_params
subroutine finish_pjastrow_params
implicit none
if(.not.opt_jastrow)return
deallocate(xyzzyaamb1)
call xyzzyaaog1
end subroutine finish_pjastrow_params
subroutine xyzzyaaof1
implicit none
integer xyzzyaaaa37,xyzzyaaab37
if(xyzzyaalq1(xyzzyaalr1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaalr1)
allocate(xyzzyaamo1(0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','Lu')
xyzzyaamo1=0.d0
if(xyzzyaaeq1)then
allocate(xyzzyaamp1(0:xyzzyaabx1,max_spin_pairs,0:xyzzyaaaa37),stat=xy&
&zzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','u')
xyzzyaamp1=0.d0
xyzzyaams1=(xyzzyaabx1+1)*max_spin_pairs
else
allocate(xyzzyaamq1(xyzzyaafp1,0:xyzzyaaaa37),xyzzyaamr1(xyzzyaafp1,0:&
&xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','u')
xyzzyaamq1=0.d0
xyzzyaamr1=0.d0
endif
endif
if(xyzzyaalq1(xyzzyaals1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaals1)
allocate(xyzzyaamw1(0:xyzzyaaaa37),xyzzyaamx1(0:xyzzyaaaa37),xyzzyaamu&
&1(0:xyzzyaaaa37),xyzzyaamv1(0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','Lucyl')
xyzzyaamw1=0.d0
xyzzyaamx1=0.d0
xyzzyaamu1=0.d0
xyzzyaamv1=0.d0
allocate(xyzzyaamy1(0:xyzzyaaby1,0:xyzzyaabz1,max_spin_pairs,0:xyzzyaa&
&aa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','ucyl')
xyzzyaamy1=0.d0
xyzzyaamz1=(xyzzyaaby1+1)*(xyzzyaabz1+1)*max_spin_pairs
endif
if(xyzzyaalq1(xyzzyaalt1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaalt1)
allocate(xyzzyaamt1(0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','Lqcusp')
xyzzyaamt1=0.d0
endif
if(xyzzyaalq1(xyzzyaalu1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaalu1)
if(xyzzyaaeu1)then
allocate(xyzzyaana1(0:xyzzyaaaa37),xyzzyaanb1(0:xyzzyaacb1,max_spin_pa&
&irs,0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','w')
xyzzyaana1=0.d0
xyzzyaanb1=0.d0
xyzzyaang1=(xyzzyaacb1+1)*max_spin_pairs
elseif(xyzzyaaev1)then
allocate(xyzzyaanc1(0:xyzzyaaaa37),xyzzyaand1(0:xyzzyaacc1,0:xyzzyaacc&
&1,0:xyzzyaacc1,xyzzyaafu1,0:xyzzyaaaa37),xyzzyaane1(0:xyzzyaaaa37),xy&
&zzyaanf1(0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','w')
xyzzyaanc1=0.d0
xyzzyaand1=0.d0
xyzzyaane1=0.d0
xyzzyaanf1=0.d0
xyzzyaanh1=(xyzzyaacc1+1)*(xyzzyaacc1+1)*(xyzzyaacc1+1)*xyzzyaafu1
endif
endif
if(xyzzyaalq1(xyzzyaalv1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaalv1)
allocate(xyzzyaani1(xyzzyaabi1,0:xyzzyaaaa37),xyzzyaanj1(0:maxval(xyzz&
&yaace1),nspin,xyzzyaabi1,0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','X')
xyzzyaani1=0.d0
xyzzyaanj1=0.d0
xyzzyaank1=(maxval(xyzzyaace1)+1)*nspin*xyzzyaabi1
endif
if(xyzzyaalq1(xyzzyaalw1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaalw1)
allocate(xyzzyaanl1(xyzzyaabj1,0:xyzzyaaaa37),xyzzyaanm1((maxval(xyzzy&
&aacg1)+1)*(maxval(xyzzyaacg1)+1)*(maxval(xyzzyaacf1)+1),max_spin_pair&
&s,xyzzyaabj1,0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','f')
xyzzyaanl1=0.d0
xyzzyaanm1=0.d0
xyzzyaann1=(maxval(xyzzyaacg1)+1)*(maxval(xyzzyaacg1)+1)*(maxval(xyzzy&
&aacf1)+1)*max_spin_pairs*xyzzyaabj1
endif
if(xyzzyaalq1(xyzzyaalx1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaalx1)
allocate(xyzzyaano1(xyzzyaagq1,max_spin_pairs,0:xyzzyaaaa37),stat=xyzz&
&yaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','p')
xyzzyaano1=0.d0
xyzzyaanp1=xyzzyaagq1*max_spin_pairs
endif
if(xyzzyaalq1(xyzzyaaly1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaaly1)
allocate(xyzzyaanq1(xyzzyaags1,nspin,0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','q')
xyzzyaanq1=0.d0
xyzzyaanr1=xyzzyaags1*nspin
endif
if(xyzzyaalq1(xyzzyaalz1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaalz1)
if(xyzzyaafa1.or.xyzzyaafc1)then
allocate(xyzzyaans1(9,0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','b')
xyzzyaans1=0.d0
elseif(xyzzyaafb1)then
allocate(xyzzyaant1(2,0:xyzzyaaaa37),xyzzyaanu1(2,0:xyzzyaaaa37),xyzzy&
&aanv1(2,0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','b')
xyzzyaant1=0.d0
xyzzyaanu1=0.d0
xyzzyaanv1=0.d0
elseif(xyzzyaafg1)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaalz1)
allocate(xyzzyaanw1(xyzzyaaiq1,xyzzyaair1,0:xyzzyaaaa37),stat=xyzzyaaa&
&b37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','ex2d')
xyzzyaanw1=0.d0
endif
endif
if(xyzzyaalq1(xyzzyaama1)>0)then
xyzzyaaaa37=xyzzyaalq1(xyzzyaama1)
allocate(xyzzyaany1(xyzzyaabn1,2,0:xyzzyaacd1,0:xyzzyaacd1,0:xyzzyaaaa&
&37),xyzzyaanx1(xyzzyaabn1,2,2,0:xyzzyaaaa37),stat=xyzzyaaab37)
call check_alloc(xyzzyaaab37,'SETUP_JASTROW_PBUFFER','d')
xyzzyaany1=0.d0
xyzzyaanx1=0.d0
endif
end subroutine xyzzyaaof1
subroutine xyzzyaaog1
implicit none
if(xyzzyaalq1(xyzzyaalr1)>0)then
if(xyzzyaaeq1)then
deallocate(xyzzyaamo1,xyzzyaamp1)
else
deallocate(xyzzyaamq1,xyzzyaamr1,xyzzyaamo1)
endif
endif
if(xyzzyaalq1(xyzzyaals1)>0)deallocate(xyzzyaamu1,xyzzyaamv1,xyzzyaamw&
&1,xyzzyaamx1,xyzzyaamy1)
if(xyzzyaalq1(xyzzyaalt1)>0)deallocate(xyzzyaamt1)
if(xyzzyaalq1(xyzzyaalu1)>0)then
if(xyzzyaaeu1)then
deallocate(xyzzyaana1,xyzzyaanb1)
elseif(xyzzyaaev1)then
deallocate(xyzzyaanc1,xyzzyaand1,xyzzyaane1,xyzzyaanf1)
endif
endif
if(xyzzyaalq1(xyzzyaalv1)>0)deallocate(xyzzyaani1,xyzzyaanj1)
if(xyzzyaalq1(xyzzyaalw1)>0)deallocate(xyzzyaanl1,xyzzyaanm1)
if(xyzzyaalq1(xyzzyaalx1)>0)deallocate(xyzzyaano1)
if(xyzzyaalq1(xyzzyaaly1)>0)deallocate(xyzzyaanq1)
if(xyzzyaalq1(xyzzyaalz1)>0)then
if(xyzzyaafa1.or.xyzzyaafc1)then
deallocate(xyzzyaans1)
elseif(xyzzyaafb1)then
deallocate(xyzzyaant1,xyzzyaanu1,xyzzyaanv1)
elseif(xyzzyaafg1)then
deallocate(xyzzyaanw1)
endif
endif
if(xyzzyaalq1(xyzzyaama1)>0)deallocate(xyzzyaanx1,xyzzyaany1)
end subroutine xyzzyaaog1
subroutine get_pjastrow_params(params,has_lolim,lolim,has_hilim,hilim,&
&is_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,&
&label)
implicit none
real(dp),intent(inout) :: params(xyzzyaalp1),lolim(xyzzyaalp1),hilim(x&
&yzzyaalp1)
logical,intent(inout) :: has_lolim(xyzzyaalp1),has_hilim(xyzzyaalp1),i&
&s_shallow(xyzzyaalp1),is_redundant(xyzzyaalp1),is_linear(xyzzyaalp1),&
&is_loglinear(xyzzyaalp1),has_aderiv(xyzzyaalp1),affect_map(xyzzyaalp1&
&,xyzzyaalp1)
character(2),intent(inout) :: label(xyzzyaalp1)
integer xyzzyaaaa39,xyzzyaaab39,xyzzyaaac39,xyzzyaaad39,xyzzyaaae39,xy&
&zzyaaaf39,xyzzyaaag39,xyzzyaaah39,xyzzyaaai39,xyzzyaaaj39,xyzzyaaak39&
&,xyzzyaaal39,xyzzyaaam39
real(dp) xyzzyaaan39,xyzzyaaao39,xyzzyaaap39,xyzzyaaaq39
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
do xyzzyaaag39=1,xyzzyaalp1
affect_map(xyzzyaaag39,xyzzyaaag39)=.true.
enddo
xyzzyaaan39=1.1d-8
if(isperiodic)xyzzyaaao39=0.999999d0*wigner_seitz_radius
xyzzyaaag39=0
xyzzyaaaj39=0
if(xyzzyaalq1(xyzzyaalr1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaalr1)
label(xyzzyaaai39:xyzzyaaaj39)='Ju'
if(xyzzyaacx1==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaach1
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
if(isperiodic)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=xyzzyaaao39
endif
is_shallow(xyzzyaaag39)=.true.
if(allocated(xyzzyaaaj1))then
is_redundant(xyzzyaaag39)=all(xyzzyaaaj1(:,:)==0.d0)
else
is_redundant(xyzzyaaag39)=.true.
endif
endif
if(xyzzyaaeq1)then
do xyzzyaaab39=1,xyzzyaafp1
if(xyzzyaaat1(0,xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaab1(0,xyzzyaaab39)
is_loglinear(xyzzyaaag39)=.true.
endif
do xyzzyaaac39=2,xyzzyaabx1
if(xyzzyaaat1(xyzzyaaac39,xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaab1(xyzzyaaac39,xyzzyaaab39)
is_loglinear(xyzzyaaag39)=.true.
endif
enddo
enddo
else
do xyzzyaaab39=1,xyzzyaafp1
if(xyzzyaaba1(xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaao1(xyzzyaaab39)
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
endif
if(xyzzyaaaz1(xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaan1(xyzzyaaab39)
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
endif
enddo
endif
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (U).')
endif
if(xyzzyaalq1(xyzzyaals1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaals1)
label(xyzzyaaai39:xyzzyaaaj39)='Jc'
if(xyzzyaadd1==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaacs1
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=0.d0
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=pi
is_shallow(xyzzyaaag39)=.true.
if(allocated(xyzzyaaak1))then
is_redundant(xyzzyaaag39)=all(xyzzyaaak1(:,:,:)==0.d0)
else
is_redundant(xyzzyaaag39)=.true.
endif
endif
if(xyzzyaade1==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaact1
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=0.d0
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=twopi
is_shallow(xyzzyaaag39)=.true.
if(allocated(xyzzyaaak1))then
is_redundant(xyzzyaaag39)=all(xyzzyaaak1(:,:,:)==0.d0)
else
is_redundant(xyzzyaaag39)=.true.
endif
endif
if(xyzzyaacy1==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaci1
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
if(periodicity==2.or.periodicity==3)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=xyzzyaaao39
endif
is_shallow(xyzzyaaag39)=.true.
if(allocated(xyzzyaaak1))then
is_redundant(xyzzyaaag39)=all(xyzzyaaak1(:,:,:)==0.d0)
else
is_redundant(xyzzyaaag39)=.true.
endif
endif
if(xyzzyaacz1==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaack1
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
if(periodicity==1.or.periodicity==3)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=xyzzyaaao39
endif
is_shallow(xyzzyaaag39)=.true.
if(allocated(xyzzyaaak1))then
is_redundant(xyzzyaaag39)=all(xyzzyaaak1(:,:,:)==0.d0)
else
is_redundant(xyzzyaaag39)=.true.
endif
endif
do xyzzyaaab39=1,xyzzyaafq1
do xyzzyaaad39=0,xyzzyaabz1
do xyzzyaaac39=0,xyzzyaaby1
if(xyzzyaaau1(xyzzyaaac39,xyzzyaaad39,xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaac1(xyzzyaaac39,xyzzyaaad39,xyzzyaaab39)
is_loglinear(xyzzyaaag39)=.true.
endif
enddo
enddo
enddo
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (Ucyl).')
endif
if(xyzzyaalq1(xyzzyaalt1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaalt1)
label(xyzzyaaai39:xyzzyaaaj39)='J^'
if(xyzzyaadc1==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaacq1
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
if(isperiodic)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=xyzzyaaao39
endif
is_shallow(xyzzyaaag39)=.true.
endif
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (QCUSP).')
endif
if(xyzzyaalq1(xyzzyaalu1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaalu1)
label(xyzzyaaai39:xyzzyaaaj39)='J3'
if(xyzzyaaeu1)then
if(xyzzyaada1==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaacl1
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
if(isperiodic)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=xyzzyaaao39
endif
is_shallow(xyzzyaaag39)=.true.
is_redundant(xyzzyaaag39)=all(xyzzyaaap1(:,:)==0.d0)
endif
do xyzzyaaab39=1,xyzzyaaft1
do xyzzyaaac39=0,xyzzyaacb1
if(xyzzyaabb1(xyzzyaaac39,xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaap1(xyzzyaaac39,xyzzyaaab39)
endif
enddo
enddo
elseif(xyzzyaaev1)then
xyzzyaaak39=0
if(xyzzyaadb1==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaacm1
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
if(isperiodic)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=xyzzyaaao39
endif
is_shallow(xyzzyaaag39)=.true.
is_redundant(xyzzyaaag39)=all(xyzzyaaar1(:,:,:,:)==0.d0)
xyzzyaaak39=xyzzyaaag39
endif
do xyzzyaaab39=1,xyzzyaafu1
do xyzzyaaae39=0,xyzzyaacc1
do xyzzyaaad39=0,xyzzyaacc1
do xyzzyaaac39=0,xyzzyaacc1
if(xyzzyaabc1(xyzzyaaac39,xyzzyaaad39,xyzzyaaae39,xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaar1(xyzzyaaac39,xyzzyaaad39,xyzzyaaae39,xyz&
&zyaaab39)*xyzzyaacm1**(xyzzyaaae39+xyzzyaaad39+xyzzyaaac39)
is_loglinear(xyzzyaaag39)=.true.
if(xyzzyaaak39>0)affect_map(xyzzyaaag39,xyzzyaaak39)=.true.
endif
enddo
enddo
enddo
enddo
endif
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (W/H).')
endif
if(xyzzyaalq1(xyzzyaalv1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaalv1)
label(xyzzyaaai39:xyzzyaaaj39)='JX'
do xyzzyaaaa39=1,xyzzyaabi1
if(xyzzyaadf1(xyzzyaaaa39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaacu1(xyzzyaaaa39)
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
if(isperiodic)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=xyzzyaaao39
endif
is_shallow(xyzzyaaag39)=.true.
is_redundant(xyzzyaaag39)=all(xyzzyaaad1(:,:,xyzzyaaaa39)==0.d0)
endif
do xyzzyaaab39=1,xyzzyaafv1(xyzzyaaaa39)
if(xyzzyaaav1(0,xyzzyaaab39,xyzzyaaaa39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaad1(0,xyzzyaaab39,xyzzyaaaa39)
is_loglinear(xyzzyaaag39)=.true.
endif
do xyzzyaaad39=2,xyzzyaace1(xyzzyaaaa39)
if(xyzzyaaav1(xyzzyaaad39,xyzzyaaab39,xyzzyaaaa39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaad1(xyzzyaaad39,xyzzyaaab39,xyzzyaaaa39)
is_loglinear(xyzzyaaag39)=.true.
endif
enddo
enddo
enddo
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (Chi).')
endif
if(xyzzyaalq1(xyzzyaalw1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaalw1)
label(xyzzyaaai39:xyzzyaaaj39)='Jf'
do xyzzyaaaa39=1,xyzzyaabj1
if(xyzzyaadg1(xyzzyaaaa39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaacv1(xyzzyaaaa39)
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
if(isperiodic)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=0.5*xyzzyaaao39
endif
is_shallow(xyzzyaaag39)=.true.
is_redundant(xyzzyaaag39)=all(xyzzyaaae1(:,:,:,:,xyzzyaaaa39)==0.d0)
endif
do xyzzyaaab39=1,xyzzyaafw1(xyzzyaaaa39)
do xyzzyaaae39=0,xyzzyaacf1(xyzzyaaaa39)
do xyzzyaaad39=0,xyzzyaacg1(xyzzyaaaa39)
do xyzzyaaac39=xyzzyaaad39,xyzzyaacg1(xyzzyaaaa39)
if(xyzzyaaaw1(xyzzyaaac39,xyzzyaaad39,xyzzyaaae39,xyzzyaaab39,xyzzyaaa&
&a39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaae1(xyzzyaaac39,xyzzyaaad39,xyzzyaaae39,xyz&
&zyaaab39,xyzzyaaaa39)
is_loglinear(xyzzyaaag39)=.true.
endif
enddo
enddo
enddo
enddo
enddo
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (F).')
endif
if(xyzzyaalq1(xyzzyaalx1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaalx1)
label(xyzzyaaai39:xyzzyaaaj39)='Jp'
do xyzzyaaab39=1,xyzzyaafr1
do xyzzyaaaf39=1,xyzzyaagq1
if(xyzzyaaax1(xyzzyaaaf39,xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaaf1(xyzzyaaaf39,xyzzyaaab39)
is_loglinear(xyzzyaaag39)=.true.
endif
enddo
enddo
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (P).')
endif
if(xyzzyaalq1(xyzzyaaly1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaaly1)
label(xyzzyaaai39:xyzzyaaaj39)='Jq'
do xyzzyaaab39=1,xyzzyaafs1
do xyzzyaaaf39=1,xyzzyaagr1
if(xyzzyaaay1(xyzzyaaaf39,xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaag1(xyzzyaaaf39,xyzzyaaab39)
is_loglinear(xyzzyaaag39)=.true.
endif
enddo
enddo
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (Q).')
endif
if(xyzzyaalq1(xyzzyaalz1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaalz1)
label(xyzzyaaai39:xyzzyaaaj39)='Jb'
if(xyzzyaafa1.or.xyzzyaafc1)then
do xyzzyaaah39=1,9
if(xyzzyaaij1(xyzzyaaah39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaih1(xyzzyaaah39)
select case(xyzzyaaah39)
case(1,2,3,4,7,9)
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=0.d0
case(5,6,8)
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=0.d0
end select
endif
enddo
elseif(xyzzyaafb1)then
do xyzzyaaah39=1,2
if(xyzzyaaio1(xyzzyaaah39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaail1(xyzzyaaah39)
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
endif
if(xyzzyaaip1(xyzzyaaah39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaim1(xyzzyaaah39)
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
endif
enddo
elseif(xyzzyaafg1)then
do xyzzyaaab39=1,xyzzyaair1
do xyzzyaaaf39=1,xyzzyaaiq1
if(xyzzyaaik1(xyzzyaaaf39,xyzzyaaab39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaii1(xyzzyaaaf39,xyzzyaaab39)
if(xyzzyaaiu1)then
if(xyzzyaaaf39==2)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=0.d0
elseif(xyzzyaaaf39==3)then
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=0.d0
endif
else
if(xyzzyaaaf39==1)then
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=0.d0
elseif(xyzzyaaaf39==2)then
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=0.d0
endif
endif
endif
enddo
enddo
endif
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (BIEX/EX2D).')
endif
if(xyzzyaalq1(xyzzyaama1)>0)then
xyzzyaaai39=xyzzyaaaj39+1
xyzzyaaaj39=xyzzyaaaj39+xyzzyaalq1(xyzzyaama1)
label(xyzzyaaai39:xyzzyaaaj39)='Jd'
do xyzzyaaaa39=1,xyzzyaabn1
xyzzyaaaq39=0.d0
xyzzyaaap39=0.d0
do xyzzyaaal39=1,xyzzyaabm1
if(xyzzyaabv1(xyzzyaaal39)==xyzzyaaaa39)then
xyzzyaaap39=rionion(4,xyzzyaabu1(xyzzyaaal39,1),xyzzyaabu1(xyzzyaaal39&
&,2))
if(xyzzyaaaq39==0.d0.or.xyzzyaaaq39>xyzzyaaap39)xyzzyaaaq39=xyzzyaaap3&
&9
endif
enddo
do xyzzyaaam39=1,2
do xyzzyaaah39=1,2-xyzzyaabw1(xyzzyaaaa39)
if(xyzzyaadh1(xyzzyaaaa39,xyzzyaaam39,xyzzyaaah39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaacw1(xyzzyaaaa39,xyzzyaaam39,xyzzyaaah39)
has_lolim(xyzzyaaag39)=.true.
lolim(xyzzyaaag39)=xyzzyaaan39
has_hilim(xyzzyaaag39)=.true.
hilim(xyzzyaaag39)=0.5d0*xyzzyaaaq39
is_shallow(xyzzyaaag39)=.true.
is_redundant(xyzzyaaag39)=all(xyzzyaaas1(xyzzyaaaa39,xyzzyaaam39,:,:)=&
&=0.d0)
endif
enddo
do xyzzyaaad39=0,xyzzyaacd1
do xyzzyaaae39=xyzzyaaad39*xyzzyaabw1(xyzzyaaaa39),xyzzyaacd1
if(xyzzyaabe1(xyzzyaaaa39,xyzzyaaam39,xyzzyaaad39,xyzzyaaae39)==1)then
xyzzyaaag39=xyzzyaaag39+1
params(xyzzyaaag39)=xyzzyaaas1(xyzzyaaaa39,xyzzyaaam39,xyzzyaaad39,xyz&
&zyaaae39)
is_loglinear(xyzzyaaag39)=.true.
endif
enddo
enddo
enddo
enddo
if(xyzzyaaaj39/=xyzzyaaag39)call errstop('GET_JASTROW_PARAMS','Bad par&
&am count (D).')
endif
end subroutine get_pjastrow_params
subroutine put_pjastrow_params(params,ignore,iparam_buffer,prestore,ba&
&d_params,restore_hint)
implicit none
integer,intent(in) :: iparam_buffer,restore_hint
real(dp),intent(inout) :: params(xyzzyaalp1)
logical,intent(in) :: ignore(xyzzyaalp1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa40,xyzzyaaab40,xyzzyaaac40,xyzzyaaad40,xyzzyaaae40,xy&
&zzyaaaf40,xyzzyaaag40,xyzzyaaah40,xyzzyaaai40,xyzzyaaaj40,xyzzyaaak40
logical xyzzyaaal40,xyzzyaaam40
bad_params=.false.
if(prestore)then
call xyzzyaaoi1(iparam_buffer,restore_hint)
return
endif
xyzzyaaaj40=0
if(xyzzyaalq1(xyzzyaalr1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaalr1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
if(xyzzyaacx1==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaach1=params(xyzzyaaah40)
endif
if(xyzzyaaeq1)then
do xyzzyaaab40=1,xyzzyaafp1
if(xyzzyaaat1(0,xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaab1(0,xyzzyaaab40)=params(xyzzyaaah&
&40)
endif
do xyzzyaaac40=2,xyzzyaabx1
if(xyzzyaaat1(xyzzyaaac40,xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaab1(xyzzyaaac40,xyzzyaaab40)=params&
&(xyzzyaaah40)
endif
enddo
enddo
call xyzzyaapd1
else
do xyzzyaaab40=1,xyzzyaafp1
if(xyzzyaaba1(xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaao1(xyzzyaaab40)=params(xyzzyaaah40&
&)
endif
if(xyzzyaaaz1(xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaan1(xyzzyaaab40)=params(xyzzyaaah40&
&)
endif
enddo
call xyzzyaapd1
endif
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (U).')
endif
endif
if(xyzzyaalq1(xyzzyaals1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaals1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
if(xyzzyaadd1==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaacs1=params(xyzzyaaah40)
endif
if(xyzzyaade1==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaact1=params(xyzzyaaah40)
endif
if(xyzzyaacy1==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaci1=params(xyzzyaaah40)
endif
if(xyzzyaacz1==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaack1=params(xyzzyaaah40)
endif
do xyzzyaaab40=1,xyzzyaafq1
do xyzzyaaad40=0,xyzzyaabz1
do xyzzyaaac40=0,xyzzyaaby1
if(xyzzyaaau1(xyzzyaaac40,xyzzyaaad40,xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaac1(xyzzyaaac40,xyzzyaaad40,xyzzyaa&
&ab40)=params(xyzzyaaah40)
endif
enddo
enddo
enddo
call xyzzyaape1
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (Ucyl).')
endif
endif
if(xyzzyaalq1(xyzzyaalt1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaalt1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
if(xyzzyaadc1==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaacq1=params(xyzzyaaah40)
endif
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (QCUSP).')
endif
endif
if(xyzzyaalq1(xyzzyaalu1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaalu1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
if(xyzzyaaeu1)then
if(xyzzyaada1==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaacl1=params(xyzzyaaah40)
endif
do xyzzyaaab40=1,xyzzyaaft1
do xyzzyaaac40=0,xyzzyaacb1
if(xyzzyaabb1(xyzzyaaac40,xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaap1(xyzzyaaac40,xyzzyaaab40)=params&
&(xyzzyaaah40)
endif
enddo
enddo
call xyzzyaapf1
elseif(xyzzyaaev1)then
if(xyzzyaadb1==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))then
xyzzyaacm1=params(xyzzyaaah40)
xyzzyaacn1=1.d0/xyzzyaacm1
xyzzyaaco1=xyzzyaacn1*xyzzyaacn1
do xyzzyaaab40=1,xyzzyaafu1
call xyzzyaaov1(xyzzyaaab40)
enddo
endif
endif
do xyzzyaaab40=1,xyzzyaafu1
do xyzzyaaae40=0,xyzzyaacc1
do xyzzyaaad40=0,xyzzyaacc1
do xyzzyaaac40=0,xyzzyaacc1
if(xyzzyaabc1(xyzzyaaac40,xyzzyaaad40,xyzzyaaae40,xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))then
params(xyzzyaaah40)=params(xyzzyaaah40)*xyzzyaacn1**(xyzzyaaae40+xyzzy&
&aaad40+xyzzyaaac40)
xyzzyaaar1(xyzzyaaac40,xyzzyaaad40,xyzzyaaae40,xyzzyaaab40)=params(xyz&
&zyaaah40)
endif
endif
enddo
enddo
enddo
enddo
call xyzzyaaox1
endif
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (W/H).')
endif
endif
if(xyzzyaalq1(xyzzyaalv1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaalv1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
do xyzzyaaaa40=1,xyzzyaabi1
if(xyzzyaadf1(xyzzyaaaa40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaacu1(xyzzyaaaa40)=params(xyzzyaaah40&
&)
endif
do xyzzyaaab40=1,xyzzyaafv1(xyzzyaaaa40)
if(xyzzyaaav1(0,xyzzyaaab40,xyzzyaaaa40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaad1(0,xyzzyaaab40,xyzzyaaaa40)=para&
&ms(xyzzyaaah40)
endif
do xyzzyaaad40=2,xyzzyaace1(xyzzyaaaa40)
if(xyzzyaaav1(xyzzyaaad40,xyzzyaaab40,xyzzyaaaa40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaad1(xyzzyaaad40,xyzzyaaab40,xyzzyaa&
&aa40)=params(xyzzyaaah40)
endif
enddo
enddo
enddo
call xyzzyaapg1
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (Chi).')
endif
endif
if(xyzzyaalq1(xyzzyaalw1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaalw1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
do xyzzyaaaa40=1,xyzzyaabj1
xyzzyaaam40=.false.
if(xyzzyaadg1(xyzzyaaaa40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))then
xyzzyaacv1(xyzzyaaaa40)=params(xyzzyaaah40)
call xyzzyaaol1(xyzzyaaaa40)
xyzzyaaam40=.true.
endif
endif
do xyzzyaaab40=1,xyzzyaafw1(xyzzyaaaa40)
xyzzyaaal40=.false.
do xyzzyaaae40=0,xyzzyaacf1(xyzzyaaaa40)
do xyzzyaaad40=0,xyzzyaacg1(xyzzyaaaa40)
do xyzzyaaac40=xyzzyaaad40,xyzzyaacg1(xyzzyaaaa40)
if(xyzzyaaaw1(xyzzyaaac40,xyzzyaaad40,xyzzyaaae40,xyzzyaaab40,xyzzyaaa&
&a40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))then
xyzzyaaae1(xyzzyaaac40,xyzzyaaad40,xyzzyaaae40,xyzzyaaab40,xyzzyaaaa40&
&)=params(xyzzyaaah40)
xyzzyaaal40=.true.
endif
endif
enddo
enddo
enddo
if(xyzzyaaal40.or.xyzzyaaam40)then
call xyzzyaaoo1(xyzzyaaaa40,xyzzyaaab40)
call xyzzyaapi1(xyzzyaaaa40,s_only=xyzzyaaab40)
endif
enddo
enddo
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (F).')
endif
endif
if(xyzzyaalq1(xyzzyaalx1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaalx1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
do xyzzyaaab40=1,xyzzyaafr1
do xyzzyaaaf40=1,xyzzyaagq1
if(xyzzyaaax1(xyzzyaaaf40,xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaaf1(xyzzyaaaf40,xyzzyaaab40)=params&
&(xyzzyaaah40)
endif
enddo
enddo
call xyzzyaapl1
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (P).')
endif
endif
if(xyzzyaalq1(xyzzyaaly1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaaly1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
do xyzzyaaab40=1,xyzzyaafs1
do xyzzyaaaf40=1,xyzzyaagr1
if(xyzzyaaay1(xyzzyaaaf40,xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaag1(xyzzyaaaf40,xyzzyaaab40)=params&
&(xyzzyaaah40)
endif
enddo
enddo
call xyzzyaapm1
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (Q).')
endif
endif
if(xyzzyaalq1(xyzzyaalz1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaalz1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
if(xyzzyaafa1.or.xyzzyaafc1)then
do xyzzyaaag40=1,9
if(xyzzyaaij1(xyzzyaaag40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaih1(xyzzyaaag40)=params(xyzzyaaah40&
&)
endif
enddo
if(.not.fix_holes)then
if(xyzzyaaij1(4)==-1)xyzzyaaih1(4)=xyzzyaaih1(2)
endif
elseif(xyzzyaafb1)then
do xyzzyaaag40=1,2
if(xyzzyaaio1(xyzzyaaag40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))then
xyzzyaail1(xyzzyaaag40)=params(xyzzyaaah40)
xyzzyaain1(xyzzyaaag40)=1.d0/xyzzyaail1(xyzzyaaag40)
endif
endif
if(xyzzyaaip1(xyzzyaaag40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaim1(xyzzyaaag40)=params(xyzzyaaah40&
&)
endif
enddo
elseif(xyzzyaafg1)then
do xyzzyaaab40=1,xyzzyaair1
do xyzzyaaaf40=1,xyzzyaaiq1
if(xyzzyaaik1(xyzzyaaaf40,xyzzyaaab40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))xyzzyaaii1(xyzzyaaaf40,xyzzyaaab40)=params&
&(xyzzyaaah40)
endif
enddo
enddo
endif
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (BIEX/EX2D).')
endif
endif
if(xyzzyaalq1(xyzzyaama1)>0)then
xyzzyaaah40=xyzzyaaaj40
xyzzyaaai40=xyzzyaaaj40+1
xyzzyaaaj40=xyzzyaaaj40+xyzzyaalq1(xyzzyaama1)
if(any(.not.ignore(xyzzyaaai40:xyzzyaaaj40)))then
do xyzzyaaaa40=1,xyzzyaabn1
do xyzzyaaak40=1,2
do xyzzyaaag40=1,2-xyzzyaabw1(xyzzyaaaa40)
if(xyzzyaadh1(xyzzyaaaa40,xyzzyaaak40,xyzzyaaag40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))then
xyzzyaacw1(xyzzyaaaa40,xyzzyaaak40,xyzzyaaag40)=params(xyzzyaaah40)
endif
endif
enddo
do xyzzyaaad40=0,xyzzyaacd1
do xyzzyaaae40=xyzzyaaad40*xyzzyaabw1(xyzzyaaaa40),xyzzyaacd1
if(xyzzyaabe1(xyzzyaaaa40,xyzzyaaak40,xyzzyaaad40,xyzzyaaae40)==1)then
xyzzyaaah40=xyzzyaaah40+1
if(.not.ignore(xyzzyaaah40))then
xyzzyaaas1(xyzzyaaaa40,xyzzyaaak40,xyzzyaaad40,xyzzyaaae40)=params(xyz&
&zyaaah40)
endif
endif
enddo
enddo
if(xyzzyaabw1(xyzzyaaaa40)==1)then
xyzzyaacw1(xyzzyaaaa40,xyzzyaaak40,2)=xyzzyaacw1(xyzzyaaaa40,xyzzyaaak&
&40,1)
do xyzzyaaad40=0,xyzzyaacd1
do xyzzyaaae40=xyzzyaaad40+1,xyzzyaacd1
xyzzyaaas1(xyzzyaaaa40,xyzzyaaak40,xyzzyaaae40,xyzzyaaad40)=xyzzyaaas1&
&(xyzzyaaaa40,xyzzyaaak40,xyzzyaaad40,xyzzyaaae40)
enddo
enddo
endif
enddo
enddo
if(xyzzyaaaj40/=xyzzyaaah40)call errstop('PUT_JASTROW_PARAMS','Bad par&
&am count (D).')
endif
endif
call xyzzyaaoh1(iparam_buffer)
end subroutine put_pjastrow_params
subroutine xyzzyaaoh1(indx)
implicit none
integer,intent(in) :: indx
integer xyzzyaaaa41,xyzzyaaab41
xyzzyaaab41=0
xyzzyaaaa41=0
if(indx/=0)then
xyzzyaaab41=xyzzyaamb1(indx)
xyzzyaaaa41=indx-sum(xyzzyaalq1(1:xyzzyaaab41-1),xyzzyaalq1(1:xyzzyaaa&
&b41-1)>0)
endif
if(xyzzyaalq1(xyzzyaalr1)>0.and.(xyzzyaaab41==xyzzyaalr1.or.xyzzyaaab4&
&1==0))then
xyzzyaamo1(xyzzyaaaa41)=xyzzyaach1
if(xyzzyaaeq1)then
call dcopy(xyzzyaams1,xyzzyaaaj1(0,1),1,xyzzyaamp1(0,1,xyzzyaaaa41),1)
else
call dcopy(xyzzyaafp1,xyzzyaaao1(1),1,xyzzyaamq1(1,xyzzyaaaa41),1)
call dcopy(xyzzyaafp1,xyzzyaaan1(1),1,xyzzyaamr1(1,xyzzyaaaa41),1)
endif
endif
if(xyzzyaalq1(xyzzyaals1)>0.and.(xyzzyaaab41==xyzzyaals1.or.xyzzyaaab4&
&1==0))then
xyzzyaamw1(xyzzyaaaa41)=xyzzyaacs1
xyzzyaamx1(xyzzyaaaa41)=xyzzyaact1
xyzzyaamu1(xyzzyaaaa41)=xyzzyaaci1
xyzzyaamv1(xyzzyaaaa41)=xyzzyaack1
call dcopy(xyzzyaamz1,xyzzyaaak1(0,0,1),1,xyzzyaamy1(0,0,1,xyzzyaaaa41&
&),1)
endif
if(xyzzyaalq1(xyzzyaalt1)>0.and.(xyzzyaaab41==xyzzyaalt1.or.xyzzyaaab4&
&1==0))xyzzyaamt1(xyzzyaaaa41)=xyzzyaacq1
if(xyzzyaalq1(xyzzyaalu1)>0.and.(xyzzyaaab41==xyzzyaalu1.or.xyzzyaaab4&
&1==0))then
if(xyzzyaaeu1)then
xyzzyaana1(xyzzyaaaa41)=xyzzyaacl1
call dcopy(xyzzyaang1,xyzzyaaaq1(0,1),1,xyzzyaanb1(0,1,xyzzyaaaa41),1)
elseif(xyzzyaaev1)then
xyzzyaanc1(xyzzyaaaa41)=xyzzyaacm1
call dcopy(xyzzyaanh1,xyzzyaaar1(0,0,0,1),1,xyzzyaand1(0,0,0,1,xyzzyaa&
&aa41),1)
xyzzyaane1(xyzzyaaaa41)=xyzzyaacn1
xyzzyaanf1(xyzzyaaaa41)=xyzzyaaco1
endif
endif
if(xyzzyaalq1(xyzzyaalv1)>0.and.(xyzzyaaab41==xyzzyaalv1.or.xyzzyaaab4&
&1==0))then
call dcopy(xyzzyaabi1,xyzzyaacu1(1),1,xyzzyaani1(1,xyzzyaaaa41),1)
call dcopy(xyzzyaank1,xyzzyaaal1(0,1,1),1,xyzzyaanj1(0,1,1,xyzzyaaaa41&
&),1)
endif
if(xyzzyaalq1(xyzzyaalw1)>0.and.(xyzzyaaab41==xyzzyaalw1.or.xyzzyaaab4&
&1==0))then
call dcopy(xyzzyaabj1,xyzzyaacv1(1),1,xyzzyaanl1(1,xyzzyaaaa41),1)
call dcopy(xyzzyaann1,xyzzyaaam1(1,1,1),1,xyzzyaanm1(1,1,1,xyzzyaaaa41&
&),1)
endif
if(xyzzyaalq1(xyzzyaalx1)>0.and.(xyzzyaaab41==xyzzyaalx1.or.xyzzyaaab4&
&1==0))then
call dcopy(xyzzyaanp1,xyzzyaaah1(1,1),1,xyzzyaano1(1,1,xyzzyaaaa41),1)
endif
if(xyzzyaalq1(xyzzyaaly1)>0.and.(xyzzyaaab41==xyzzyaaly1.or.xyzzyaaab4&
&1==0))then
call dcopy(xyzzyaanr1,xyzzyaaai1(1,1),1,xyzzyaanq1(1,1,xyzzyaaaa41),1)
endif
if(xyzzyaalq1(xyzzyaalz1)>0.and.(xyzzyaaab41==xyzzyaalz1.or.xyzzyaaab4&
&1==0))then
if(xyzzyaafa1.or.xyzzyaafc1)then
call dcopy(9,xyzzyaaih1(1),1,xyzzyaans1(1,xyzzyaaaa41),1)
elseif(xyzzyaafb1)then
call dcopy(2,xyzzyaail1(1),1,xyzzyaant1(1,xyzzyaaaa41),1)
call dcopy(2,xyzzyaain1(1),1,xyzzyaanv1(1,xyzzyaaaa41),1)
call dcopy(2,xyzzyaaim1(1),1,xyzzyaanu1(1,xyzzyaaaa41),1)
elseif(xyzzyaafg1)then
call dcopy(xyzzyaaiq1*xyzzyaair1,xyzzyaaii1(1,1),1,xyzzyaanw1(1,1,xyzz&
&yaaaa41),1)
endif
endif
if(xyzzyaalq1(xyzzyaama1)>0.and.(xyzzyaaab41==xyzzyaama1.or.xyzzyaaab4&
&1==0))then
xyzzyaanx1(:,:,:,xyzzyaaaa41)=xyzzyaacw1(:,:,:)
xyzzyaany1(:,:,:,:,xyzzyaaaa41)=xyzzyaaas1(:,:,:,:)
endif
end subroutine xyzzyaaoh1
subroutine xyzzyaaoi1(indx,restore_hint)
implicit none
integer,intent(in) :: indx,restore_hint
integer xyzzyaaaa42,xyzzyaaab42
xyzzyaaab42=0
xyzzyaaaa42=0
if(indx/=0)then
xyzzyaaab42=xyzzyaamb1(indx)
xyzzyaaaa42=indx-sum(xyzzyaalq1(1:xyzzyaaab42-1),xyzzyaalq1(1:xyzzyaaa&
&b42-1)>0)
elseif(restore_hint/=0)then
xyzzyaaab42=xyzzyaamb1(restore_hint)
xyzzyaaaa42=0
endif
if(xyzzyaalq1(xyzzyaalr1)>0.and.(xyzzyaaab42==xyzzyaalr1.or.xyzzyaaab4&
&2==0))then
xyzzyaach1=xyzzyaamo1(xyzzyaaaa42)
if(xyzzyaaeq1)then
call dcopy(xyzzyaams1,xyzzyaamp1(0,1,xyzzyaaaa42),1,xyzzyaaaj1(0,1),1)
else
call dcopy(xyzzyaafp1,xyzzyaamq1(1,xyzzyaaaa42),1,xyzzyaaao1(1),1)
call dcopy(xyzzyaafp1,xyzzyaamr1(1,xyzzyaaaa42),1,xyzzyaaan1(1),1)
endif
endif
if(xyzzyaalq1(xyzzyaals1)>0.and.(xyzzyaaab42==xyzzyaals1.or.xyzzyaaab4&
&2==0))then
xyzzyaacs1=xyzzyaamw1(xyzzyaaaa42)
xyzzyaact1=xyzzyaamx1(xyzzyaaaa42)
xyzzyaacr1=(/sin(xyzzyaacs1)*cos(xyzzyaact1),sin(xyzzyaacs1)*sin(xyzzy&
&aact1),cos(xyzzyaacs1)/)
xyzzyaaci1=xyzzyaamu1(xyzzyaaaa42)
xyzzyaacj1=xyzzyaaci1**2
xyzzyaack1=xyzzyaamv1(xyzzyaaaa42)
call dcopy(xyzzyaamz1,xyzzyaamy1(0,0,1,xyzzyaaaa42),1,xyzzyaaak1(0,0,1&
&),1)
endif
if(xyzzyaalq1(xyzzyaalt1)>0.and.(xyzzyaaab42==xyzzyaalt1.or.xyzzyaaab4&
&2==0))xyzzyaacq1=xyzzyaamt1(xyzzyaaaa42)
if(xyzzyaalq1(xyzzyaalu1)>0.and.(xyzzyaaab42==xyzzyaalu1.or.xyzzyaaab4&
&2==0))then
if(xyzzyaaeu1)then
xyzzyaacl1=xyzzyaana1(xyzzyaaaa42)
call dcopy(xyzzyaang1,xyzzyaanb1(0,1,xyzzyaaaa42),1,xyzzyaaaq1(0,1),1)
elseif(xyzzyaaev1)then
xyzzyaacm1=xyzzyaanc1(xyzzyaaaa42)
call dcopy(xyzzyaanh1,xyzzyaand1(0,0,0,1,xyzzyaaaa42),1,xyzzyaaar1(0,0&
&,0,1),1)
xyzzyaacn1=xyzzyaane1(xyzzyaaaa42)
xyzzyaaco1=xyzzyaanf1(xyzzyaaaa42)
endif
endif
if(xyzzyaalq1(xyzzyaalv1)>0.and.(xyzzyaaab42==xyzzyaalv1.or.xyzzyaaab4&
&2==0))then
call dcopy(xyzzyaabi1,xyzzyaani1(1,xyzzyaaaa42),1,xyzzyaacu1(1),1)
call dcopy(xyzzyaank1,xyzzyaanj1(0,1,1,xyzzyaaaa42),1,xyzzyaaal1(0,1,1&
&),1)
endif
if(xyzzyaalq1(xyzzyaalw1)>0.and.(xyzzyaaab42==xyzzyaalw1.or.xyzzyaaab4&
&2==0))then
call dcopy(xyzzyaabj1,xyzzyaanl1(1,xyzzyaaaa42),1,xyzzyaacv1(1),1)
call dcopy(xyzzyaann1,xyzzyaanm1(1,1,1,xyzzyaaaa42),1,xyzzyaaam1(1,1,1&
&),1)
endif
if(xyzzyaalq1(xyzzyaalx1)>0.and.(xyzzyaaab42==xyzzyaalx1.or.xyzzyaaab4&
&2==0))then
call dcopy(xyzzyaanp1,xyzzyaano1(1,1,xyzzyaaaa42),1,xyzzyaaah1(1,1),1)
endif
if(xyzzyaalq1(xyzzyaaly1)>0.and.(xyzzyaaab42==xyzzyaaly1.or.xyzzyaaab4&
&2==0))then
call dcopy(xyzzyaanr1,xyzzyaanq1(1,1,xyzzyaaaa42),1,xyzzyaaai1(1,1),1)
endif
if(xyzzyaalq1(xyzzyaalz1)>0.and.(xyzzyaaab42==xyzzyaalz1.or.xyzzyaaab4&
&2==0))then
if(xyzzyaafa1.or.xyzzyaafc1)then
call dcopy(9,xyzzyaans1(1,xyzzyaaaa42),1,xyzzyaaih1(1),1)
elseif(xyzzyaafb1)then
call dcopy(2,xyzzyaant1(1,xyzzyaaaa42),1,xyzzyaail1(1),1)
call dcopy(2,xyzzyaanv1(1,xyzzyaaaa42),1,xyzzyaain1(1),1)
call dcopy(2,xyzzyaanu1(1,xyzzyaaaa42),1,xyzzyaaim1(1),1)
elseif(xyzzyaafg1)then
call dcopy(xyzzyaaiq1*xyzzyaair1,xyzzyaanw1(1,1,xyzzyaaaa42),1,xyzzyaa&
&ii1(1,1),1)
endif
endif
if(xyzzyaalq1(xyzzyaama1)>0.and.(xyzzyaaab42==xyzzyaama1.or.xyzzyaaab4&
&2==0))then
xyzzyaacw1(:,:,:)=xyzzyaanx1(:,:,:,xyzzyaaaa42)
xyzzyaaas1(:,:,:,:)=xyzzyaany1(:,:,:,:,xyzzyaaaa42)
endif
end subroutine xyzzyaaoi1
subroutine invalidate_param1_pjastrow(is,iparam)
implicit none
integer,intent(in) :: is,iparam
integer xyzzyaaaa43
xyzzyaako1(:,is)=.false.
xyzzyaakq1(:,is)=.false.
xyzzyaakr1(:,is)=.false.
xyzzyaaks1(:,is)=.false.
xyzzyaakt1(:,is)=.false.
xyzzyaakp1(is)=.false.
xyzzyaaaa43=xyzzyaamb1(iparam)
select case(xyzzyaaaa43)
case(xyzzyaalr1)
xyzzyaaku1(:,is)=.false.
xyzzyaakv1(:,is)=.false.
case(xyzzyaals1)
xyzzyaakw1(:,is)=.false.
xyzzyaakx1(:,is)=.false.
case(xyzzyaalt1)
xyzzyaaky1(:,is)=.false.
xyzzyaakz1(:,is)=.false.
case(xyzzyaalu1)
xyzzyaala1(:,is)=.false.
xyzzyaalb1(:,is)=.false.
case(xyzzyaalv1)
xyzzyaalc1(:,is)=.false.
xyzzyaald1(:,is)=.false.
case(xyzzyaalw1)
xyzzyaale1(:,is)=.false.
xyzzyaalf1(:,is)=.false.
case(xyzzyaalx1)
xyzzyaalg1(:,is)=.false.
xyzzyaalh1(:,is)=.false.
case(xyzzyaaly1)
xyzzyaali1(:,is)=.false.
xyzzyaalj1(:,is)=.false.
case(xyzzyaalz1)
xyzzyaalk1(is)=.false.
xyzzyaall1(:,is)=.false.
case(xyzzyaama1)
xyzzyaalm1(:,is)=.false.
xyzzyaaln1(:,is)=.false.
end select
end subroutine invalidate_param1_pjastrow
subroutine invalidate_params_pjastrow(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaalp1)
integer xyzzyaaaa44,xyzzyaaab44,xyzzyaaac44
logical xyzzyaaad44(xyzzyaalo1)
xyzzyaaad44=.false.
xyzzyaaab44=0
do xyzzyaaac44=1,xyzzyaalo1
if(xyzzyaalq1(xyzzyaaac44)>0)then
xyzzyaaaa44=xyzzyaaab44+1
xyzzyaaab44=xyzzyaaab44+xyzzyaalq1(xyzzyaaac44)
xyzzyaaad44(xyzzyaaac44)=any(.not.ignore(xyzzyaaaa44:xyzzyaaab44))
endif
enddo
if(any(xyzzyaaad44))then
xyzzyaamk1(:)=.false.
xyzzyaaml1(:)=.false.
xyzzyaamm1(:)=.false.
xyzzyaamn1(:)=.false.
endif
end subroutine invalidate_params_pjastrow
subroutine clear_scratch_pjastrow(is)
implicit none
integer,intent(in) :: is
xyzzyaako1(:,is)=.false.
xyzzyaakq1(:,is)=.false.
xyzzyaakr1(:,is)=.false.
xyzzyaaks1(:,is)=.false.
xyzzyaakt1(:,is)=.false.
xyzzyaakp1(is)=.false.
xyzzyaaku1(:,is)=.false.
xyzzyaakv1(:,is)=.false.
xyzzyaakw1(:,is)=.false.
xyzzyaakx1(:,is)=.false.
xyzzyaaky1(:,is)=.false.
xyzzyaakz1(:,is)=.false.
xyzzyaala1(:,is)=.false.
xyzzyaalb1(:,is)=.false.
xyzzyaalc1(:,is)=.false.
xyzzyaald1(:,is)=.false.
xyzzyaale1(:,is)=.false.
xyzzyaalf1(:,is)=.false.
xyzzyaalg1(:,is)=.false.
xyzzyaalh1(:,is)=.false.
xyzzyaali1(:,is)=.false.
xyzzyaalj1(:,is)=.false.
xyzzyaalk1(is)=.false.
xyzzyaall1(:,is)=.false.
xyzzyaalm1(:,is)=.false.
xyzzyaaln1(:,is)=.false.
if(xyzzyaaja1==is)xyzzyaaja1=0
if(xyzzyaajb1==is)xyzzyaajb1=0
if(xyzzyaajc1==is)xyzzyaajc1=0
end subroutine clear_scratch_pjastrow
subroutine enumerate_plot_pjastrow(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
n=0
end subroutine enumerate_plot_pjastrow
subroutine query_plot_pjastrow(iplot,ii,rank,is_complex,has_stderr,rot&
&_tensor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
nfunctions=0
end subroutine query_plot_pjastrow
subroutine get_plot_pjastrow(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
f(1)=0.d0
end subroutine get_plot_pjastrow
subroutine finish_plot_pjastrow
implicit none
end subroutine finish_plot_pjastrow
subroutine read_pjastrow(empty_jastrow,gen_gjastrow)
use slaarnaan,only : atoms_label_species,atoms_label_pcell
implicit none
logical,intent(inout) :: empty_jastrow
logical,intent(in) :: gen_gjastrow
integer xyzzyaaaa50,xyzzyaaab50,xyzzyaaac50,xyzzyaaad50,xyzzyaaae50,l,&
&m,n,xyzzyaaaf50,xyzzyaaag50,xyzzyaaah50(3),xyzzyaaai50(5),xyzzyaaaj50&
&,xyzzyaaak50,xyzzyaaal50,xyzzyaaam50,xyzzyaaan50,xyzzyaaao50,xyzzyaaa&
&p50,xyzzyaaaq50,xyzzyaaar50,xyzzyaaas50,xyzzyaaat50,xyzzyaaau50
integer,allocatable :: xyzzyaaav50(:,:,:),xyzzyaaaw50(:,:,:,:,:)
real(dp) xyzzyaaax50
real(dp),allocatable :: xyzzyaaay50(:,:,:),xyzzyaaaz50(:,:,:,:,:)
real(dp),parameter :: xyzzyaaba50=1.d-10
logical xyzzyaabb50,xyzzyaabc50
logical,allocatable :: xyzzyaabd50(:),xyzzyaabe50(:,:)
character(18) char_18
character(80) tmpr,char_80
xyzzyaaep1=.false.
xyzzyaaeq1=.false.
xyzzyaaer1=.false.
xyzzyaaew1=.false.
xyzzyaaex1=.false.
xyzzyaaey1=.false.
xyzzyaaez1=.false.
xyzzyaaeu1=.false.
xyzzyaaev1=.false.
xyzzyaafe1=.false.
xyzzyaafd1=.false.
xyzzyaaff1=.false.
xyzzyaaet1=.false.
xyzzyaaes1=.false.
xyzzyaafg1=.false.
empty_jastrow=.true.
call open_units(xyzzyaagw1,xyzzyaaaa50)
if(xyzzyaaaa50/=0)call errstop('READ_PJASTROW','Unable to find free i/&
&o unit.')
open(unit=xyzzyaagw1,file='correlation.data',status='old',iostat=xyzzy&
&aaaa50)
if(xyzzyaaaa50/=0)call errstop('READ_PJASTROW','Problem opening correl&
&ation.data')
if(am_master)then
call wout('Jastrow factor')
call wout('==============')
call wout('Reading Jastrow factor from correlation.data file.')
call wout()
endif
do
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50>0)call errstop_master('READ_PJASTROW','Problem reading &
&correlation.data.  Please check this file.')
if(xyzzyaaaa50<0)call errstop_master('READ_PJASTROW','Could not find "&
&START JASTROW" in correlation.data.')
if(trim(adjustl(char_80))=='START JASTROW')exit
enddo
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,'(a)',err=666,end=666)title
if(am_master)then
call wout('Title: '//trim(adjustl(title)))
call wout()
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaagx1
if(am_master)then
call wout('Truncation order                   :  '//trim(i2s(xyzzyaagx&
&1)))
if(xyzzyaagx1<3.and.xyzzyaagx1>0)call wout('Local energy is discontinu&
&ous at cutoffs.')
if(xyzzyaagx1==1)call errwarn('READ_PJASTROW','derivative discontinuit&
&y at Jastrow cutoffs!')
if(xyzzyaagx1==0)call errwarn('READ_PJASTROW','no cutoff will be appli&
&ed on the Jastrow functions!')
if(xyzzyaagx1<0.and..not.isperiodic)call errstop('READ_PJASTROW','The &
&truncation order should be non-negative.')
if(xyzzyaagx1<2.and.isperiodic)call errstop('READ_PJASTROW','The trunc&
&ation order should be at least 2 in a periodic system.')
call wout()
endif
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='START U TERM')then
if(xyzzyaaeq1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about u terms?')
if(xyzzyaaer1)call errstop_master('READ_PJASTROW','Cannot have both U &
&and U_RPA.')
xyzzyaaep1=.true.
xyzzyaaeq1=.true.
xyzzyaaer1=.false.
if(am_master)call wout('U term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabh1
if(xyzzyaabh1/=1)call errstop_master('READ_PJASTROW','Should only have&
& one set of u terms.')
if(am_master)call wout(' SET 1')
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='START SET 1')exit
if(trim(adjustl(char_80))/="".and.am_master)call errstop('READ_PJASTRO&
&W','Expecting to find "START SET 1".')
enddo
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaaha1,xyzzyaahb1
if(am_master)then
call wout('  Spherical harmonic l             :  '//trim(i2s(xyzzyaaha&
&1)))
call wout('  Spherical harmonic m             :  '//trim(i2s(xyzzyaahb&
&1)))
if(xyzzyaaha1/=0.or.xyzzyaahb1/=0)call errstop('READ_PJASTROW','Need t&
&o have l=m=0 at present.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabx1
if(am_master)then
call wout('  Expansion order (N_u)            :  '//trim(i2s(xyzzyaabx&
&1)))
if(xyzzyaabx1<1)call errstop('READ_PJASTROW','N_u<1.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafh1
if(am_master)then
call wout('  Spin dependence                  :  '//trim(i2s(xyzzyaafh&
&1)))
if(xyzzyaafh1<-custom_spairs.or.xyzzyaafh1>levels_spairs)call errstop(&
&'READ_PJASTROW','Spin dep should be '//trim(i2s(-custom_spairs))//' -&
&- '//trim(i2s(levels_spairs))//' for U.')
if(xyzzyaafh1/=0.and.noncoll_spin)call errstop('READ_PJASTROW','Spin-d&
&ependence of U should be zero for noncollinear-spin calculations.')
endif
xyzzyaafp1=no_spairs(xyzzyaafh1)
allocate(xyzzyaaab1(0:xyzzyaabx1,xyzzyaafp1),xyzzyaaat1(0:xyzzyaabx1,x&
&yzzyaafp1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','1')
xyzzyaaab1(:,:)=0.d0
xyzzyaaat1(:,:)=1
xyzzyaaat1(1,:)=-1
if(xyzzyaagx1==0)xyzzyaaat1(0,:)=0
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaach1,xyzzyaacx1
if(abs(xyzzyaach1)<xyzzyaahi1)then
xyzzyaach1=xyzzyaaqr1()
if(am_master)call wout('  Using default cutoff length L_u.')
endif
if(am_master)then
call display_param(xyzzyaach1,xyzzyaacx1,'Cutoff')
if(xyzzyaach1<0.d0)call errstop('READ_PJASTROW','L_u<0.')
if(isperiodic.and.xyzzyaach1>wigner_seitz_radius)call errstop('READ_PJ&
&ASTROW','L_u > radius of sphere inscribed in Wigner-Seitz cell.')
if(xyzzyaacx1/=0.and.xyzzyaacx1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
if(xyzzyaacx1==1.and.xyzzyaagx1<2)call errstop('READ_PJASTROW','Cannot&
& optimize L_u if truncation order < 2.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read_params_alpha: do xyzzyaaae50=1,xyzzyaafp1
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaab1(0,xyzzyaaae50),xyzzyaa&
&at1(0,xyzzyaaae50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaab1(0,xyzzyaaae50)=0.d0
xyzzyaaat1(0,xyzzyaaae50)=1
exit read_params_alpha
endif
empty_jastrow=.false.
if(xyzzyaagx1==0.and.xyzzyaaat1(0,xyzzyaaae50)==1)then
call errwarn('READ_PJASTROW','alpha_0 is not an optimizable parameter &
&if trunc. order is 0.')
xyzzyaaat1(0,xyzzyaaae50)=0
endif
if(am_master)then
call display_param(xyzzyaaab1(0,xyzzyaaae50),xyzzyaaat1(0,xyzzyaaae50)&
&,'alpha_0,'//trim(i2s(xyzzyaaae50)))
if(xyzzyaaat1(0,xyzzyaaae50)/=0.and.xyzzyaaat1(0,xyzzyaaae50)/=1)call &
&errstop('READ_PJASTROW','The optimizable flag should be 0 or 1.')
endif
do l=2,xyzzyaabx1
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaab1(l,xyzzyaaae50),xyzzyaa&
&at1(l,xyzzyaaae50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaab1(l,xyzzyaaae50)=0.d0
xyzzyaaat1(l,xyzzyaaae50)=1
exit read_params_alpha
endif
if(am_master)then
if(l==1)then
call display_param(xyzzyaaab1(l,xyzzyaaae50),xyzzyaaat1(l,xyzzyaaae50)&
&,'alpha_'//trim(i2s(l))//','//trim(i2s(xyzzyaaae50)),comment_in='*')
else
call display_param(xyzzyaaab1(l,xyzzyaaae50),xyzzyaaat1(l,xyzzyaaae50)&
&,'alpha_'//trim(i2s(l))//','//trim(i2s(xyzzyaaae50)))
endif
if(xyzzyaaat1(l,xyzzyaaae50)/=0.and.xyzzyaaat1(l,xyzzyaaae50)/=1)call &
&errstop('READ_PJASTROW','The optimizable flag should be 0 or 1.')
endif
enddo
enddo read_params_alpha
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END SET 1".')
if(trim(adjustl(char_80))/='END SET 1')call errstop_master('READ_PJAST&
&ROW','Was expecting to find "END SET 1".')
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END U TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "END U TERM".')
enddo
if(am_master)then
if(xyzzyaabx1>=2)then
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaabx&
&1*xyzzyaafp1)))
else
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaafp&
&1)))
endif
if(xyzzyaagx1>0)call wout('  (In addition to the cutoff length.)')
call wout()
endif
allocate(xyzzyaaaj1(0:xyzzyaabx1,max_spin_pairs),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','1.5')
xyzzyaaaj1=0.d0
if(hard_sphere)empty_jastrow=.false.
elseif(trim(adjustl(char_80))=='START UCYL TERM')then
if(xyzzyaaes1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about ucyl terms?')
if(dimensionality/=3)call errstop('READ_PJASTROW','Can only have Ucyl &
&term in a 3D system.')
xyzzyaaes1=.true.
if(am_master)call wout('Ucyl term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaaby1
if(am_master)then
call wout('  Expansion order (N_ucylrho)       :  '//trim(i2s(xyzzyaab&
&y1)))
if(xyzzyaaby1<1)call errstop('READ_PJASTROW','N_ucylrho<1.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabz1
if(am_master)then
call wout('  Expansion order (N_ucylz)         :  '//trim(i2s(xyzzyaab&
&z1)))
if(xyzzyaabz1<1)call errstop('READ_PJASTROW','N_ucylz<1.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafi1
if(am_master)then
call wout('  Spin dependence                  :  '//trim(i2s(xyzzyaafi&
&1)))
if(xyzzyaafi1<-custom_spairs.or.xyzzyaafi1>levels_spairs)call errstop(&
&'READ_PJASTROW','Spin dep should be '//trim(i2s(-custom_spairs))//' -&
&- '//trim(i2s(levels_spairs))//' for Ucyl.')
if(xyzzyaafi1/=0.and.noncoll_spin)call errstop('READ_PJASTROW','Spin-d&
&ependence of Ucyl should be zero for noncollinear-spin calculations.'&
&)
endif
xyzzyaafq1=no_spairs(xyzzyaafi1)
allocate(xyzzyaaac1(0:xyzzyaaby1,0:xyzzyaabz1,xyzzyaafq1),xyzzyaaau1(0&
&:xyzzyaaby1,0:xyzzyaabz1,xyzzyaafq1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','1ucyl')
xyzzyaaac1=0.d0
xyzzyaaau1=1
xyzzyaaau1(1,:,:)=-1
xyzzyaaau1(:,1,:)=-1
if(xyzzyaagx1==0)xyzzyaaau1(0,0,:)=0
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacs1,xyzzyaadd1
if(am_master)then
call display_param(xyzzyaacs1,xyzzyaadd1,'Polar angle')
if(xyzzyaacs1<0.d0.or.xyzzyaacs1>pi)call errstop('READ_PJASTROW','Pola&
&r angle theta for ucyl term should be between 0 and pi.')
if(xyzzyaadd1/=0.and.xyzzyaadd1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
if(periodicity==2.and.xyzzyaacs1/=0.d0)call errstop('READ_PJASTROW','P&
&olar angle for ucyl term should be zero in a 2D-periodic calculation.&
&')
if(periodicity==1.and.abs(xyzzyaacs1-pi_over_two)>xyzzyaaba50)call err&
&stop('READ_PJASTROW','Polar angle for ucyl term should be pi/2 in a 1&
&D-periodic calculation.')
if((periodicity==1.or.periodicity==2).and.xyzzyaadd1==1)call errstop('&
&READ_PJASTROW','The polar angle should not be optimisable in a '//tri&
&m(i2s(periodicity))//'D-periodic calculation.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaact1,xyzzyaade1
if(am_master)then
call display_param(xyzzyaact1,xyzzyaade1,'Azimuthal angle')
if(xyzzyaact1<0.d0.or.xyzzyaact1>=twopi)call errstop('READ_PJASTROW','&
&Azimuthal angle phi for ucyl term should be between 0 and 2.pi.')
if(xyzzyaade1/=0.and.xyzzyaade1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
if(periodicity==1.and.xyzzyaact1/=0.d0)call errstop('READ_PJASTROW','A&
&zimuthal angle for ucyl term should be 0 in a 1D-periodic calculation&
&.')
if((periodicity==1.or.periodicity==2).and.xyzzyaade1==1)call errstop('&
&READ_PJASTROW','The azimuthal angle should not be optimisable in a '/&
&/trim(i2s(periodicity))//'D-periodic calculation.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaaci1,xyzzyaacy1
if(abs(xyzzyaaci1)<xyzzyaahi1)then
xyzzyaaci1=xyzzyaaqs1()
if(am_master)call wout('  Using default cutoff length L_ucylrho.')
endif
if(am_master)then
call display_param(xyzzyaaci1,xyzzyaacy1,'Radial cutoff')
if(xyzzyaaci1<0.d0)call errstop('READ_PJASTROW','L_ucylrho<0.')
if((periodicity==2.or.periodicity==3).and.xyzzyaaci1>wigner_seitz_radi&
&us)call errstop('READ_PJASTROW','L_ucylrho > radius of circle inscrib&
&ed in Wigner-Seitz cell.')
if(xyzzyaacy1/=0.and.xyzzyaacy1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
if(xyzzyaacy1==1.and.xyzzyaagx1<2)call errstop('READ_PJASTROW','Cannot&
& optimize L_ucylrho if truncation order < 2.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaack1,xyzzyaacz1
if(abs(xyzzyaack1)<xyzzyaahi1)then
xyzzyaack1=xyzzyaaqt1()
if(am_master)call wout('  Using default cutoff length L_ucylz.')
endif
if(am_master)then
call display_param(xyzzyaack1,xyzzyaacz1,'Axial cutoff')
if(xyzzyaack1<0.d0)call errstop('READ_PJASTROW','L_ucylz<0.')
if((periodicity==1.or.periodicity==3).and.xyzzyaaci1>wigner_seitz_radi&
&us)call errstop('READ_PJASTROW','L_ucylrho > half-length of Wigner-Se&
&itz cell.')
if(xyzzyaacz1/=0.and.xyzzyaacz1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
if(xyzzyaacz1==1.and.xyzzyaagx1<2)call errstop('READ_PJASTROW','Cannot&
& optimize L_ucylz if truncation order < 2.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read_params_epsil: do xyzzyaaae50=1,xyzzyaafq1
do m=0,xyzzyaabz1
if(m==1)cycle
do l=0,xyzzyaaby1
if(l==1)cycle
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaac1(l,m,xyzzyaaae50),xyzzy&
&aaau1(l,m,xyzzyaaae50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaac1(l,m,xyzzyaaae50)=0.d0
xyzzyaaau1(l,m,xyzzyaaae50)=1
exit read_params_epsil
endif
empty_jastrow=.false.
if(am_master)then
if(xyzzyaagx1==0.and.l==0.and.m==0.and.xyzzyaaau1(0,0,xyzzyaaae50)==1)&
&then
call errwarn('READ_PJASTROW','epsil_0,0 is not an optimizable paramete&
&r if trunc. order is 0.')
xyzzyaaau1(0,0,xyzzyaaae50)=0
endif
call display_param(xyzzyaaac1(l,m,xyzzyaaae50),xyzzyaaau1(l,m,xyzzyaaa&
&e50),'epsil_'//trim(i2s(l))//','//trim(i2s(m))//','//trim(i2s(xyzzyaa&
&ae50)))
if(xyzzyaaau1(l,m,xyzzyaaae50)/=0.and.xyzzyaaau1(l,m,xyzzyaaae50)/=1)c&
&all errstop('READ_PJASTROW','The optimizable flag should be 0 or 1.')
endif
enddo
enddo
enddo read_params_epsil
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END UCYL TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "END UCYL TERM".')
enddo
if(am_master)then
if(xyzzyaaby1>1)then
if(xyzzyaabz1>1)then
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaaby&
&1*xyzzyaabz1*xyzzyaafq1)))
else
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaaby&
&1*xyzzyaafq1)))
endif
else
if(xyzzyaabz1>1)then
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaabz&
&1*xyzzyaafq1)))
else
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaafq&
&1)))
endif
endif
if(xyzzyaagx1>0)call wout('  (In addition to the 2 angles and the 2 cu&
&toff lengths.)')
call wout()
endif
allocate(xyzzyaaak1(0:xyzzyaaby1,0:xyzzyaabz1,max_spin_pairs),stat=xyz&
&zyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','1.5ucyl')
xyzzyaaak1=0.d0
elseif(trim(adjustl(char_80))=='START QCUSP TERM')then
if(heg_nlayers==1)call errstop_master('READ_PJASTROW','Cannot use a QC&
&USP term except on multilayers and multiwires.')
if(xyzzyaaet1)call errstop_master('READ_PJASTROW','Can only have one Q&
&CUSP term.')
xyzzyaaet1=.true.
if(am_master)call wout('QCUSP term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacq1,xyzzyaadc1
empty_jastrow=.false.
if(abs(xyzzyaacq1)<xyzzyaahi1)then
xyzzyaacq1=1.d0
if(isperiodic)xyzzyaacq1=min(xyzzyaacq1,0.999999*wigner_seitz_radius)
if(am_master)call wout('  Using default cutoff length L_qcusp.')
endif
if(am_master)then
call display_param(xyzzyaacq1,xyzzyaadc1,'Cutoff')
if(xyzzyaacq1<0.d0)call errstop('READ_PJASTROW','L_qcusp<0.')
if(isperiodic.and.xyzzyaacq1>wigner_seitz_radius)call errstop('READ_PJ&
&ASTROW','L_qcusp > radius of sphere inscribed in Wigner-Seitz cell.')
if(xyzzyaadc1/=0.and.xyzzyaadc1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
endif
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END QCUSP TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "END QCUSP TERM".')
enddo
if(am_master)then
call wout('  No. of parameters in set         :  '//trim(i2s(1)))
call wout()
endif
elseif(trim(adjustl(char_80))=='START U_RPA TERM')then
if(xyzzyaaer1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about U_RPA terms?')
if(xyzzyaaeq1)call errstop_master('READ_PJASTROW','Cannot have both U &
&and U_RPA.')
xyzzyaaep1=.true.
xyzzyaaeq1=.false.
xyzzyaaer1=.true.
if(am_master)call wout('U (RPA) term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabh1
if(xyzzyaabh1/=1)call errstop_master('READ_PJASTROW','Should only have&
& one set of u_RPA terms.')
if(am_master)call wout(' SET 1')
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='START SET 1')exit
if(trim(adjustl(char_80))/="".and.am_master)call errstop('READ_PJASTRO&
&W','Expecting to find "START SET 1".')
enddo
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaaha1,xyzzyaahb1
if(am_master)then
call wout('  Spherical harmonic l             :  '//trim(i2s(xyzzyaaha&
&1)))
call wout('  Spherical harmonic m             :  '//trim(i2s(xyzzyaahb&
&1)))
if(xyzzyaaha1/=0.or.xyzzyaahb1/=0)call errstop('READ_PJASTROW','Need t&
&o have l=m=0 at present.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafh1
if(am_master)then
call wout('  Spin dependence                  :  '//trim(i2s(xyzzyaafh&
&1)))
if(xyzzyaafh1<-custom_spairs.or.xyzzyaafh1>levels_spairs)call errstop(&
&'READ_PJASTROW','Spin dep should be '//trim(i2s(-custom_spairs))//' -&
&- '//trim(i2s(levels_spairs))//' for U_RPA.')
if(xyzzyaafh1/=0.and.noncoll_spin)call errstop('READ_PJASTROW','Spin-d&
&ependence of U_RPA should be zero for noncollinear-spin calculations.&
&')
endif
xyzzyaafp1=no_spairs(xyzzyaafh1)
allocate(xyzzyaaao1(xyzzyaafp1),xyzzyaaan1(xyzzyaafp1),xyzzyaaba1(xyzz&
&yaafp1),xyzzyaaaz1(xyzzyaafp1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','1rpa')
xyzzyaaao1(:)=1.d0
xyzzyaaba1(:)=1
xyzzyaaan1(:)=1.d0
xyzzyaaaz1(:)=-1
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaach1,xyzzyaacx1
if(abs(xyzzyaach1)<xyzzyaahi1)then
xyzzyaach1=xyzzyaaqr1()
if(am_master)call wout('  Using default cutoff length L_u.')
endif
if(am_master)then
call display_param(xyzzyaach1,xyzzyaacx1,'Cutoff')
if(xyzzyaach1<0.d0)call errstop('READ_PJASTROW','L_u<0.')
if(isperiodic.and.xyzzyaach1>wigner_seitz_radius)call errstop('READ_PJ&
&ASTROW','L_u > radius of sphere inscribed in Wigner-Seitz cell.')
if(xyzzyaacx1/=0.and.xyzzyaacx1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
if(xyzzyaacx1==1.and.xyzzyaagx1<2)call errstop('READ_PJASTROW','Cannot&
& optimize L_u if truncation order < 2.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read_params_rpa: do xyzzyaaae50=1,xyzzyaafp1
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaao1(xyzzyaaae50),xyzzyaaba&
&1(xyzzyaaae50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaao1(xyzzyaaae50)=1.d0
xyzzyaaao1(xyzzyaaae50)=1
exit read_params_rpa
endif
empty_jastrow=.false.
if(am_master)then
call display_param(xyzzyaaao1(xyzzyaaae50),xyzzyaaba1(xyzzyaaae50),'F_&
&rpa_'//trim(i2s(xyzzyaaae50)))
if(xyzzyaaba1(xyzzyaaae50)/=0.and.xyzzyaaba1(xyzzyaaae50)/=1)call errs&
&top('READ_PJASTROW','The optimizable flag should be 0 or 1.')
endif
enddo read_params_rpa
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END SET 1".')
if(trim(adjustl(char_80))/='END SET 1')call errstop_master('READ_PJAST&
&ROW','Was expecting to find "END SET 1".')
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END U_RPA TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "END U_RPA TERM".')
enddo
if(am_master)then
xyzzyaaan50=xyzzyaafp1
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaaan&
&50)))
if(xyzzyaagx1>0)call wout('  (In addition to the cutoff length.)')
call wout()
endif
elseif(trim(adjustl(char_80))=='START W TERM')then
if(xyzzyaaeu1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about W terms?')
if(xyzzyaaev1)call errstop_master('READ_PJASTROW','Cannot use W and H &
&terms at the same time. Try using H, which is more general.')
xyzzyaaeu1=.true.
if(am_master)call wout('W term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabk1
if(xyzzyaabk1/=1)call errstop_master('READ_PJASTROW','Should only have&
& one set of W terms.')
if(am_master)call wout(' SET 1')
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='START SET 1')exit
if(trim(adjustl(char_80))/="".and.am_master)call errstop('READ_PJASTRO&
&W','Expecting to find "START SET 1".')
enddo
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaahc1,xyzzyaahd1
if(am_master)then
call wout('  Spherical harmonic l             :  '//trim(i2s(xyzzyaahc&
&1)))
call wout('  Spherical harmonic m             :  '//trim(i2s(xyzzyaahd&
&1)))
if(xyzzyaahc1/=0.or.xyzzyaahd1/=0)call errstop('READ_PJASTROW','Need t&
&o have l=m=0 at present.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacb1
if(am_master)then
call wout('  Expansion order (N_w)            :  '//trim(i2s(xyzzyaacb&
&1)))
if(xyzzyaacb1<1)call errstop('READ_PJASTROW','N_w<1.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafl1
if(am_master)then
call wout('  Spin dependence                  :  '//trim(i2s(xyzzyaafl&
&1)))
if(xyzzyaafl1<-custom_spairs.or.xyzzyaafl1>levels_spairs)call errstop(&
&'READ_PJASTROW','Spin dep should be '//trim(i2s(-custom_spairs))//' -&
&- '//trim(i2s(levels_spairs))//' for W.')
if(xyzzyaafl1/=0.and.noncoll_spin)call errstop('READ_PJASTROW','Spin-d&
&ependence of W should be zero for noncollinear-spin calculations.')
endif
xyzzyaaft1=no_spairs(xyzzyaafl1)
allocate(xyzzyaaap1(0:xyzzyaacb1,xyzzyaaft1),xyzzyaabb1(0:xyzzyaacb1,x&
&yzzyaaft1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','1w')
xyzzyaaap1(:,:)=0.d0
xyzzyaabb1(:,:)=1
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacl1,xyzzyaada1
if(abs(xyzzyaacl1)<xyzzyaahi1)then
xyzzyaacl1=xyzzyaaqu1()
if(am_master)call wout('  Using default cutoff length L_w.')
endif
if(am_master)then
call display_param(xyzzyaacl1,xyzzyaada1,'Cutoff')
if(xyzzyaacl1<0.d0)call errstop('READ_PJASTROW','L_w<0.')
if(isperiodic.and.xyzzyaacl1>wigner_seitz_radius)call errstop('READ_PJ&
&ASTROW','L_w > radius of sphere inscribed in Wigner-Seitz cell.')
if(xyzzyaada1/=0.and.xyzzyaada1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
if(xyzzyaada1==1.and.xyzzyaagx1<2)call errstop('READ_PJASTROW','Cannot&
& optimize L_w if truncation order < 2.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read_params_w: do xyzzyaaae50=1,xyzzyaaft1
do l=0,xyzzyaacb1
if(xyzzyaabb1(l,xyzzyaaae50)<0)cycle
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaap1(l,xyzzyaaae50),xyzzyaa&
&bb1(l,xyzzyaaae50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaap1(l,xyzzyaaae50)=0.d0
xyzzyaabb1(l,xyzzyaaae50)=1
exit read_params_w
endif
empty_jastrow=.false.
if(am_master)then
call display_param(xyzzyaaap1(l,xyzzyaaae50),xyzzyaabb1(l,xyzzyaaae50)&
&,'w_'//trim(i2s(l))//','//trim(i2s(xyzzyaaae50)))
if(xyzzyaabb1(l,xyzzyaaae50)/=0.and.xyzzyaabb1(l,xyzzyaaae50)/=1)call &
&errstop('READ_PJASTROW','The optimizable flag should be 0 or 1.')
endif
enddo
enddo read_params_w
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END SET 1".')
if(trim(adjustl(char_80))/='END SET 1')call errstop_master('READ_PJAST&
&ROW','Was expecting to find "END SET 1".')
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END W TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "END W TERM".')
enddo
if(am_master)then
xyzzyaaan50=count(xyzzyaabb1>=0)
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaaan&
&50)))
if(xyzzyaagx1>0)call wout('  (In addition to the cutoff length.)')
call wout()
endif
xyzzyaaac50=netot*(netot-1)/2
allocate(xyzzyaaaq1(0:xyzzyaacb1,max_spin_pairs),xyzzyaadz1(3,netot),x&
&yzzyaaea1(3,xyzzyaaac50),xyzzyaaeb1(3,3,xyzzyaaac50),xyzzyaaec1(3,xyz&
&zyaaac50),xyzzyaael1(netot),xyzzyaaem1(xyzzyaaac50),xyzzyaaen1(3,xyzz&
&yaaac50),xyzzyaaeo1(xyzzyaaac50),xyzzyaaie1(netot),xyzzyaaed1(netot),&
&xyzzyaaee1(netot),xyzzyaaef1(3,3,netot),xyzzyaaeg1(xyzzyaaac50),xyzzy&
&aaeh1(netot),xyzzyaaei1(netot,netot),xyzzyaaej1(3,netot,netot),xyzzya&
&aek1(netot,netot),xyzzyaaif1(xyzzyaaac50),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','1.5w')
xyzzyaaif1=.false.
xyzzyaaea1=0.d0
xyzzyaaeb1=0.d0
xyzzyaaec1=0.d0
xyzzyaaen1=0.d0
elseif(trim(adjustl(char_80))=='START H TERM')then
if(xyzzyaaev1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about H terms?')
if(xyzzyaaeu1)call errstop_master('READ_PJASTROW','Cannot use W and H &
&terms at the same time. Try using H, which is more general.')
xyzzyaaev1=.true.
if(am_master)call wout('H term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabl1
if(xyzzyaabl1/=1)call errstop_master('READ_PJASTROW','Should only have&
& one set of w terms.')
if(am_master)call wout(' SET 1')
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='START SET 1')exit
if(trim(adjustl(char_80))/="".and.am_master)call errstop('READ_PJASTRO&
&W','Expecting to find "START SET 1".')
enddo
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaahe1,xyzzyaahf1
if(am_master)then
call wout('  Spherical harmonic l             :  '//trim(i2s(xyzzyaahe&
&1)))
call wout('  Spherical harmonic m             :  '//trim(i2s(xyzzyaahf&
&1)))
if(xyzzyaahe1/=0.or.xyzzyaahf1/=0)call errstop('READ_PJASTROW','Need t&
&o have l=m=0 at present.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacc1
xyzzyaahw1=(xyzzyaacc1+1)**3
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','0.9h')
if(am_master)then
call wout('  Expansion order (N_h)            :  '//trim(i2s(xyzzyaacc&
&1)))
if(xyzzyaacc1<1)call errstop('READ_PJASTROW','N_h<1.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafm1
if(am_master)then
call wout('  Spin dependence                  :  '//trim(i2s(xyzzyaafm&
&1)))
if(xyzzyaafm1<-custom_striplets.or.xyzzyaafm1>levels_striplets)call er&
&rstop('READ_PJASTROW','Spin dep should be '//trim(i2s(-custom_striple&
&ts))//' -- '//trim(i2s(levels_striplets))//' for H.')
if(xyzzyaafm1/=0.and.noncoll_spin)call errstop('READ_PJASTROW','Spin-d&
&ependence of H should be zero for noncollinear-spin calculations.')
endif
xyzzyaafu1=no_striplets(xyzzyaafm1)
allocate(xyzzyaaar1(0:xyzzyaacc1,0:xyzzyaacc1,0:xyzzyaacc1,xyzzyaafu1)&
&,xyzzyaabc1(0:xyzzyaacc1,0:xyzzyaacc1,0:xyzzyaacc1,xyzzyaafu1),xyzzya&
&aia1(xyzzyaafu1),xyzzyaahz1(xyzzyaafu1),xyzzyaabe50(xyzzyaahw1,xyzzya&
&afu1),xyzzyaabd1(3,nspin,nspin,nspin),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','1h')
xyzzyaaar1=0.d0
xyzzyaabc1=1
xyzzyaabe50=.false.
xyzzyaabd1=0
xyzzyaaia1=0
do xyzzyaaao50=1,nspin
do xyzzyaaap50=xyzzyaaao50,nspin
do xyzzyaaaq50=xyzzyaaap50,nspin
xyzzyaaae50=which_striplet(xyzzyaaao50,xyzzyaaap50,xyzzyaaaq50,xyzzyaa&
&fm1)
select case(eq_triplet(xyzzyaaao50,xyzzyaaap50,xyzzyaaaq50,xyzzyaafm1)&
&)
case(0,1,7)
xyzzyaabd1(1:3,xyzzyaaao50,xyzzyaaap50,xyzzyaaaq50)=(/1,2,3/)
xyzzyaabd1(1:3,xyzzyaaao50,xyzzyaaaq50,xyzzyaaap50)=(/1,3,2/)
xyzzyaabd1(1:3,xyzzyaaap50,xyzzyaaao50,xyzzyaaaq50)=(/2,1,3/)
xyzzyaabd1(1:3,xyzzyaaap50,xyzzyaaaq50,xyzzyaaao50)=(/3,1,2/)
xyzzyaabd1(1:3,xyzzyaaaq50,xyzzyaaao50,xyzzyaaap50)=(/2,3,1/)
xyzzyaabd1(1:3,xyzzyaaaq50,xyzzyaaap50,xyzzyaaao50)=(/3,2,1/)
case(2)
xyzzyaabd1(1:3,xyzzyaaao50,xyzzyaaap50,xyzzyaaaq50)=(/2,3,1/)
xyzzyaabd1(1:3,xyzzyaaao50,xyzzyaaaq50,xyzzyaaap50)=(/3,2,1/)
xyzzyaabd1(1:3,xyzzyaaap50,xyzzyaaao50,xyzzyaaaq50)=(/1,3,2/)
xyzzyaabd1(1:3,xyzzyaaap50,xyzzyaaaq50,xyzzyaaao50)=(/1,2,3/)
xyzzyaabd1(1:3,xyzzyaaaq50,xyzzyaaao50,xyzzyaaap50)=(/3,1,2/)
xyzzyaabd1(1:3,xyzzyaaaq50,xyzzyaaap50,xyzzyaaao50)=(/2,1,3/)
case(4)
xyzzyaabd1(1:3,xyzzyaaao50,xyzzyaaap50,xyzzyaaaq50)=(/1,3,2/)
xyzzyaabd1(1:3,xyzzyaaao50,xyzzyaaaq50,xyzzyaaap50)=(/1,2,3/)
xyzzyaabd1(1:3,xyzzyaaap50,xyzzyaaao50,xyzzyaaaq50)=(/2,3,1/)
xyzzyaabd1(1:3,xyzzyaaap50,xyzzyaaaq50,xyzzyaaao50)=(/3,2,1/)
xyzzyaabd1(1:3,xyzzyaaaq50,xyzzyaaao50,xyzzyaaap50)=(/2,1,3/)
xyzzyaabd1(1:3,xyzzyaaaq50,xyzzyaaap50,xyzzyaaao50)=(/3,1,2/)
end select
select case(eq_triplet(xyzzyaaao50,xyzzyaaap50,xyzzyaaaq50,xyzzyaafm1)&
&)
case(0)
xyzzyaaia1(xyzzyaaae50)=0
case(1,2,4)
xyzzyaaia1(xyzzyaaae50)=1
case(7)
xyzzyaaia1(xyzzyaaae50)=3
end select
enddo
enddo
enddo
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacm1,xyzzyaadb1
if(abs(xyzzyaacm1)<xyzzyaahi1)then
xyzzyaacm1=xyzzyaaqv1()
if(am_master)call wout('  Using default cutoff length L_h.')
endif
xyzzyaacn1=1.d0/xyzzyaacm1
xyzzyaaco1=xyzzyaacn1*xyzzyaacn1
if(am_master)then
call display_param(xyzzyaacm1,xyzzyaadb1,'Cutoff')
if(xyzzyaacm1<0.d0)call errstop('READ_PJASTROW','L_h<0.')
if(isperiodic.and.xyzzyaacm1>wigner_seitz_radius)call errstop('READ_PJ&
&ASTROW','L_h > radius of sphere inscribed in Wigner-Seitz cell.')
if(xyzzyaadb1/=0.and.xyzzyaadb1/=1)call errstop('READ_PJASTROW','The o&
&ptimizable flag should be 0 or 1.')
if(xyzzyaadb1==1.and.xyzzyaagx1<2)call errstop('READ_PJASTROW','Cannot&
& optimize L_h if truncation order < 2.')
endif
read(xyzzyaagw1,*,err=666,end=666)
call xyzzyaaoy1(xyzzyaabe50)
read_params_h: do xyzzyaaae50=1,xyzzyaafu1
xyzzyaaaf50=0
do n=0,xyzzyaacc1
do m=0,xyzzyaacc1
do l=0,xyzzyaacc1
xyzzyaaaf50=xyzzyaaaf50+1
if(xyzzyaabe50(xyzzyaaaf50,xyzzyaaae50))cycle
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaar1(l,m,n,xyzzyaaae50),xyz&
&zyaabc1(l,m,n,xyzzyaaae50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaar1(l,m,n,xyzzyaaae50)=0.d0
xyzzyaabc1(l,m,n,xyzzyaaae50)=1
exit read_params_h
endif
empty_jastrow=.false.
if(am_master)then
call display_param(xyzzyaaar1(l,m,n,xyzzyaaae50),xyzzyaabc1(l,m,n,xyzz&
&yaaae50),'h_'//trim(i2s(l))//','//trim(i2s(m)) //','//trim(i2s(n))//'&
&,'//trim(i2s(xyzzyaaae50)))
if(xyzzyaabc1(l,m,n,xyzzyaaae50)/=0.and.xyzzyaabc1(l,m,n,xyzzyaaae50)/&
&=1)call errstop('READ_PJASTROW','The optimizable flag should be 0 or &
&1.')
endif
enddo
enddo
enddo
enddo read_params_h
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END SET 1".')
if(trim(adjustl(char_80))/='END SET 1')call errstop_master('READ_PJAST&
&ROW','Was expecting to find "END SET 1".')
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END H TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "END H TERM".')
enddo
if(am_master)then
call wout('  No. of parameters in set         :  '//trim(i2s(count(.no&
&t.xyzzyaabe50))))
if(xyzzyaagx1>0)call wout('  (In addition to the cutoff length.)')
call wout()
endif
deallocate(xyzzyaabe50)
elseif(trim(adjustl(char_80))=='START CHI TERM')then
if(xyzzyaaew1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about chi terms?')
if(nitot<1)call errstop_master('READ_PJASTROW','No chi term needed in &
&a system without atoms.')
xyzzyaaew1=.true.
if(am_master)call wout('Chi term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
read(char_80,*,iostat=xyzzyaaaa50)xyzzyaabi1,xyzzyaabf1
if(xyzzyaaaa50/=0)then
read(char_80,*,err=666,end=666)xyzzyaabi1
xyzzyaabf1=1
endif
if(am_master)then
call wout(' Number of sets                    :  '//trim(i2s(xyzzyaabi&
&1)))
if(xyzzyaabf1<1.or.xyzzyaabf1>3)call errstop('READ_PJASTROW','Label st&
&yle should be 1, 2 or 3.')
if(xyzzyaabi1<1.or.xyzzyaabi1>nitot.or.(xyzzyaabf1==2 .and.xyzzyaabi1>&
&nbasis).or.(xyzzyaabf1==3.and.xyzzyaabi1>nitype)) call errstop('READ_&
&PJASTROW','Problematic number of chi sets of ions.')
endif
allocate(xyzzyaabo1(nitot),xyzzyaace1(xyzzyaabi1),xyzzyaacu1(xyzzyaabi&
&1),xyzzyaadf1(xyzzyaabi1),xyzzyaadi1(xyzzyaabi1),xyzzyaadn1(xyzzyaabi&
&1),xyzzyaafn1(xyzzyaabi1),xyzzyaabq1(nitot,xyzzyaabi1),xyzzyaafv1(xyz&
&zyaabi1),xyzzyaahg1(xyzzyaabi1),xyzzyaahh1(xyzzyaabi1),xyzzyaadl1(xyz&
&zyaabi1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','2')
xyzzyaabq1(1:nitot,1:xyzzyaabi1)=0
xyzzyaabo1(1:nitot)=0
if(xyzzyaabf1==1)then
allocate(xyzzyaabr1(nitot,xyzzyaabi1),stat=xyzzyaaaj50)
elseif(xyzzyaabf1==2)then
if(nbasis<=0)call errstop_master('READ_PJASTROW','Number of atoms in b&
&asis is zero.  If this is a Wigner crystal then use "atom" labels in &
&supercell.')
allocate(xyzzyaabr1(nbasis,xyzzyaabi1),stat=xyzzyaaaj50)
else
if(nitype<=0)call errstop_master('READ_PJASTROW','Number of atom speci&
&es is zero.  If this is a Wigner crystal then use "atom" labels in su&
&percell.')
allocate(xyzzyaabr1(nitype,xyzzyaabi1),stat=xyzzyaaaj50)
endif
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','2.1')
do xyzzyaaab50=1,xyzzyaabi1
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='START SET '//trim(i2s(xyzzyaaab50)))exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "START SET '//trim(i2s(xyzzyaaab50))//'".')
enddo
if(am_master)call wout(' SET '//trim(i2s(xyzzyaaab50))//':')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaahg1(xyzzyaaab50),xyzzyaahh1(x&
&yzzyaaab50)
if(am_master)then
call wout('  Spherical harmonic l             :  '//trim(i2s(xyzzyaahg&
&1(xyzzyaaab50))))
call wout('  Spherical harmonic m             :  '//trim(i2s(xyzzyaahh&
&1(xyzzyaaab50))))
if(xyzzyaahg1(xyzzyaaab50)/=0.or.xyzzyaahh1(xyzzyaaab50)/=0)call errst&
&op('READ_PJASTROW','Need to have l=m=0 at present.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaadl1(xyzzyaaab50)
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabr1(1:xyzzyaadl1(xyzzyaaab50)&
&,xyzzyaaab50)
if(xyzzyaabf1==1)then
xyzzyaadi1(xyzzyaaab50)=xyzzyaadl1(xyzzyaaab50)
xyzzyaabq1(1:xyzzyaadl1(xyzzyaaab50),xyzzyaaab50) =xyzzyaabr1(1:xyzzya&
&adl1(xyzzyaaab50),xyzzyaaab50)
elseif(xyzzyaabf1==2)then
if(xyzzyaadl1(xyzzyaaab50)<1.or.xyzzyaadl1(xyzzyaaab50)>nbasis+1-xyzzy&
&aabi1)call errstop_master('READ_PJASTROW','Problematic number of atom&
&s in set.')
if(any(xyzzyaabr1(1:xyzzyaadl1(xyzzyaaab50),xyzzyaaab50)<1).or.any(xyz&
&zyaabr1(1:xyzzyaadl1(xyzzyaaab50),xyzzyaaab50)>nbasis))call errstop_m&
&aster('READ_PJASTROW','Problem with atom labels.')
call atoms_label_pcell(xyzzyaadl1(xyzzyaaab50),xyzzyaabr1(1:xyzzyaadl1&
&(xyzzyaaab50),xyzzyaaab50),xyzzyaadi1(xyzzyaaab50),xyzzyaabq1(:,xyzzy&
&aaab50))
else
if(xyzzyaadl1(xyzzyaaab50)<1.or.xyzzyaadl1(xyzzyaaab50)>nitype+1-xyzzy&
&aabi1)call errstop_master('READ_PJASTROW','Problematic number of spec&
&ies in set.')
if(any(xyzzyaabr1(1:xyzzyaadl1(xyzzyaaab50),xyzzyaaab50)<1).or.any(xyz&
&zyaabr1(1:xyzzyaadl1(xyzzyaaab50),xyzzyaaab50)>nitype))call errstop_m&
&aster('READ_PJASTROW','Problem with species labels.')
call atoms_label_species(xyzzyaadl1(xyzzyaaab50),xyzzyaabr1(1:xyzzyaad&
&l1(xyzzyaaab50),xyzzyaaab50),xyzzyaadi1(xyzzyaaab50),xyzzyaabq1(:,xyz&
&zyaaab50))
endif
if(am_master)then
call wout('  Number of atoms in set           :  '//trim(i2s(xyzzyaadi&
&1(xyzzyaaab50))))
if(xyzzyaadi1(xyzzyaaab50)<1.or.xyzzyaadi1(xyzzyaaab50)>nitot+1-xyzzya&
&abi1)call errstop('READ_PJASTROW','Problematic number of atoms in set&
&.')
call wout('  The atoms are:')
call write_list_int(xyzzyaadi1(xyzzyaaab50),xyzzyaabq1(1:xyzzyaadi1(xy&
&zzyaaab50),xyzzyaaab50),10,4,1)
if(any(xyzzyaabq1(1:xyzzyaadi1(xyzzyaaab50),xyzzyaaab50)<1).or.any(xyz&
&zyaabq1(1:xyzzyaadi1(xyzzyaaab50),xyzzyaaab50)>nitot))call errstop('R&
&EAD_PJASTROW','Problem with atom labels.')
do xyzzyaaac50=1,xyzzyaadi1(xyzzyaaab50)-1
do xyzzyaaad50=xyzzyaaac50+1,xyzzyaadi1(xyzzyaaab50)
if(xyzzyaabq1(xyzzyaaac50,xyzzyaaab50)==xyzzyaabq1(xyzzyaaad50,xyzzyaa&
&ab50))call errstop('READ_PJASTROW','Ion '//trim(i2s(xyzzyaabq1(xyzzya&
&aac50,xyzzyaaab50)))//' appears twice.')
enddo
enddo
endif
do xyzzyaaac50=1,xyzzyaadi1(xyzzyaaab50)
if(xyzzyaabo1(xyzzyaabq1(xyzzyaaac50,xyzzyaaab50))==0)then
xyzzyaabo1(xyzzyaabq1(xyzzyaaac50,xyzzyaaab50))=xyzzyaaab50
else
if(am_master)then
call wout('Ion '//trim(i2s(xyzzyaabq1(xyzzyaaac50,xyzzyaaab50)))//' ap&
&pears in sets '//trim(i2s(xyzzyaabo1(xyzzyaabq1(xyzzyaaac50,xyzzyaaab&
&50))))//' and '//trim(i2s(xyzzyaaab50))//'.')
call errstop('READ_PJASTROW','Stopping.')
endif
endif
enddo
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaadn1(xyzzyaaab50)
if(am_master)then
if(xyzzyaadn1(xyzzyaaab50)==1)then
call wout('  Electron-nucleus cusp imposed in Jastrow')
else
call wout('  Electron-nucleus cusp not imposed in Jastrow')
endif
if(xyzzyaadn1(xyzzyaaab50)/=0.and.xyzzyaadn1(xyzzyaaab50)/=1)call errs&
&top('READ_PJASTROW','Impose e-N cusp flag should be 0 or 1.')
endif
if(xyzzyaadn1(xyzzyaaab50)==1)then
if(.not.allocated(xyzzyaafx1))then
allocate(xyzzyaafx1(xyzzyaabi1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','2.1.5')
endif
xyzzyaafx1(xyzzyaaab50)=xyzzyaaqy1(xyzzyaaab50)
tmpr=r2s(xyzzyaafx1(xyzzyaaab50),'(f21.12)')
if(am_master)call wout('  Z for the set                    : '//trim(t&
&mpr))
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaace1(xyzzyaaab50)
if(am_master)then
call wout('  Expansion order (N_chi)          :  '//trim(i2s(xyzzyaace&
&1(xyzzyaaab50))))
if(xyzzyaace1(xyzzyaaab50)<1)call errstop('READ_PJASTROW','N_chi<1.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafn1(xyzzyaaab50)
if(am_master)then
call wout('  Spin dependence                  :  '//trim(i2s(xyzzyaafn&
&1(xyzzyaaab50))))
if(xyzzyaafn1(xyzzyaaab50)<-custom_ssingles.or.xyzzyaafn1(xyzzyaaab50)&
&>levels_ssingles)call errstop('READ_PJASTROW','Spin-dep chi should be&
& '//trim(i2s(-custom_ssingles))//' -- '//trim(i2s(levels_ssingles))//&
&'.')
endif
xyzzyaafv1(xyzzyaaab50)=no_ssingles(xyzzyaafn1(xyzzyaaab50))
if(xyzzyaaab50==1)then
allocate(xyzzyaaad1(0:xyzzyaace1(xyzzyaaab50),xyzzyaafv1(xyzzyaaab50),&
&xyzzyaabi1),xyzzyaaav1(0:xyzzyaace1(xyzzyaaab50),xyzzyaafv1(xyzzyaaab&
&50),xyzzyaabi1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','3')
xyzzyaaad1(:,:,:)=0.d0
xyzzyaaav1(:,:,:)=1
xyzzyaaav1(1,:,:)=-1
if(xyzzyaagx1==0)xyzzyaaav1(0,:,:)=0
else
xyzzyaaah50=shape(xyzzyaaad1)
if(xyzzyaace1(xyzzyaaab50)>xyzzyaaah50(1)-1.or.xyzzyaafv1(xyzzyaaab50)&
&>xyzzyaaah50(2))then
allocate(xyzzyaaay50(0:xyzzyaaah50(1)-1,xyzzyaaah50(2),xyzzyaaab50-1),&
&xyzzyaaav50(0:xyzzyaaah50(1)-1,xyzzyaaah50(2),xyzzyaaab50-1),stat=xyz&
&zyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','4')
xyzzyaaay50(0:xyzzyaaah50(1)-1,1:xyzzyaaah50(2),1:xyzzyaaab50-1)=xyzzy&
&aaad1(0:xyzzyaaah50(1)-1,1:xyzzyaaah50(2),1:xyzzyaaab50-1)
xyzzyaaav50(0:xyzzyaaah50(1)-1,1:xyzzyaaah50(2),1:xyzzyaaab50-1)=xyzzy&
&aaav1(0:xyzzyaaah50(1)-1,1:xyzzyaaah50(2),1:xyzzyaaab50-1)
deallocate(xyzzyaaad1,xyzzyaaav1)
allocate(xyzzyaaad1(0:max(xyzzyaace1(xyzzyaaab50),xyzzyaaah50(1)-1),ma&
&x(xyzzyaafv1(xyzzyaaab50),xyzzyaaah50(2)),xyzzyaaah50(3)),xyzzyaaav1(&
&0:max(xyzzyaace1(xyzzyaaab50),xyzzyaaah50(1)-1),max(xyzzyaafv1(xyzzya&
&aab50),xyzzyaaah50(2)),xyzzyaaah50(3)),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','5')
xyzzyaaad1(:,:,:)=0.d0
xyzzyaaav1(:,:,:)=1
xyzzyaaav1(1,:,:)=-1
if(xyzzyaagx1==0)xyzzyaaav1(0,:,:)=0
xyzzyaaad1(0:xyzzyaaah50(1)-1,1:xyzzyaaah50(2),1:xyzzyaaab50-1)=xyzzya&
&aay50(0:xyzzyaaah50(1)-1,1:xyzzyaaah50(2),1:xyzzyaaab50-1)
xyzzyaaav1(0:xyzzyaaah50(1)-1,1:xyzzyaaah50(2),1:xyzzyaaab50-1)=xyzzya&
&aav50(0:xyzzyaaah50(1)-1,1:xyzzyaaah50(2),1:xyzzyaaab50-1)
deallocate(xyzzyaaay50,xyzzyaaav50)
endif
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacu1(xyzzyaaab50),xyzzyaadf1(x&
&yzzyaaab50)
if(abs(xyzzyaacu1(xyzzyaaab50))<xyzzyaahi1)then
xyzzyaacu1(xyzzyaaab50)=xyzzyaaqw1()
if(am_master)call wout('  Using default cutoff length L_chi.')
endif
if(am_master)then
call display_param(xyzzyaacu1(xyzzyaaab50),xyzzyaadf1(xyzzyaaab50),'Cu&
&toff')
if(xyzzyaacu1(xyzzyaaab50)<0.d0)call errstop('READ_PJASTROW','L_chi<0'&
&)
if(isperiodic.and.xyzzyaacu1(xyzzyaaab50)>wigner_seitz_radius)call err&
&stop('READ_PJASTROW','L_chi > radius of sphere inscribed in Wigner-Se&
&itz cell.')
if(xyzzyaadf1(xyzzyaaab50)<0.or.xyzzyaadf1(xyzzyaaab50)>1)call errstop&
&('READ_PJASTROW','Optimizable flag should be 0 or 1.')
if(xyzzyaadf1(xyzzyaaab50)==1.and.xyzzyaagx1<2)call errstop('READ_PJAS&
&TROW','Cannot optimize L_chi if truncation order < 2.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read_params_beta: do xyzzyaaae50=1,xyzzyaafv1(xyzzyaaab50)
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaad1(0,xyzzyaaae50,xyzzyaaa&
&b50),xyzzyaaav1(0,xyzzyaaae50,xyzzyaaab50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaad1(0,xyzzyaaae50,xyzzyaaab50)=0.d0
xyzzyaaav1(0,xyzzyaaae50,xyzzyaaab50)=1
exit read_params_beta
endif
empty_jastrow=.false.
if(xyzzyaagx1==0.and.xyzzyaaav1(0,xyzzyaaae50,xyzzyaaab50)==1)then
call errwarn('READ_PJASTROW','beta_0 is not an optimizable parameter i&
&f trunc. order is 0.')
xyzzyaaav1(0,xyzzyaaae50,xyzzyaaab50)=0
endif
if(am_master)then
call display_param(xyzzyaaad1(0,xyzzyaaae50,xyzzyaaab50),xyzzyaaav1(0,&
&xyzzyaaae50,xyzzyaaab50),'beta_0,'//trim(i2s(xyzzyaaae50))//','//trim&
&(i2s(xyzzyaaab50)))
if(xyzzyaaav1(0,xyzzyaaae50,xyzzyaaab50)<0.or.xyzzyaaav1(0,xyzzyaaae50&
&,xyzzyaaab50)>1)call errstop('READ_PJASTROW','Optimizable flag should&
& be 0 or 1.')
endif
do m=2,xyzzyaace1(xyzzyaaab50)
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaad1(m,xyzzyaaae50,xyzzyaaa&
&b50),xyzzyaaav1(m,xyzzyaaae50,xyzzyaaab50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaad1(m,xyzzyaaae50,xyzzyaaab50)=0.d0
xyzzyaaav1(m,xyzzyaaae50,xyzzyaaab50)=1
exit read_params_beta
endif
if(am_master)then
call display_param(xyzzyaaad1(m,xyzzyaaae50,xyzzyaaab50),xyzzyaaav1(m,&
&xyzzyaaae50,xyzzyaaab50),'beta_'//trim(i2s(m))//','//trim(i2s(xyzzyaa&
&ae50))//','//trim(i2s(xyzzyaaab50)))
if(xyzzyaaav1(m,xyzzyaaae50,xyzzyaaab50)<0.or.xyzzyaaav1(m,xyzzyaaae50&
&,xyzzyaaab50)>1)call errstop('READ_PJASTROW','Optimizable flag should&
& be 0 or 1.')
endif
enddo
enddo read_params_beta
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END SET '//trim(i2s(xyzzyaaab50))//'".')
if(trim(adjustl(char_80))/='END SET '//trim(i2s(xyzzyaaab50)))call err&
&stop_master('READ_PJASTROW','Was expecting to find "END SET ' //trim(&
&i2s(xyzzyaaab50))//'".')
if(am_master)then
if(xyzzyaace1(xyzzyaaab50)>=2)then
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaace&
&1(xyzzyaaab50)*xyzzyaafv1(xyzzyaaab50))))
else
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaafv&
&1(xyzzyaaab50))))
endif
call wout('  (In addition to the cutoff length.)')
endif
enddo
if(any(xyzzyaabo1==0).and.am_master)then
if(allow_nochi_atoms)then
call errwarn_silent('READ_PJASTROW','Some atoms do not belong to any c&
&hi set.')
else
call errstop('READ_PJASTROW','Some atoms do not belong to any chi set.&
&')
endif
endif
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END CHI TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find END CHI TERM.')
enddo
if(am_master)call wout()
allocate(xyzzyaaal1(0:maxval(xyzzyaace1),nspin,xyzzyaabi1),stat=xyzzya&
&aaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','5.5')
xyzzyaaal1=0.d0
elseif(trim(adjustl(char_80))=='START F TERM')then
if(xyzzyaaex1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about f terms?')
if(nitot<1)call errstop_master('READ_PJASTROW','No f term needed in a &
&system without atoms.')
xyzzyaaex1=.true.
if(am_master)call wout('F term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
read(char_80,*,iostat=xyzzyaaaa50)xyzzyaabj1,xyzzyaabg1
if(xyzzyaaaa50/=0)then
read(char_80,*,err=666,end=666)xyzzyaabj1
xyzzyaabg1=1
endif
if(am_master)then
call wout(' Number of sets                    :  '//trim(i2s(xyzzyaabj&
&1)))
if(xyzzyaabg1<1.or.xyzzyaabg1>3)call errstop('READ_PJASTROW','Label st&
&yle should be 1, 2 or 3.')
if(xyzzyaabj1<1.or.xyzzyaabj1>nitot.or.(xyzzyaabg1==2 .and.xyzzyaabj1>&
&nbasis).or.(xyzzyaabg1==3.and.xyzzyaabj1>nitype))call errstop('READ_P&
&JASTROW','Problematic number of f sets of ions.')
endif
allocate(xyzzyaabp1(nitot),xyzzyaacg1(xyzzyaabj1),xyzzyaacf1(xyzzyaabj&
&1),xyzzyaacv1(xyzzyaabj1),xyzzyaadg1(xyzzyaabj1),xyzzyaadj1(xyzzyaabj&
&1),xyzzyaafo1(xyzzyaabj1),xyzzyaabs1(nitot,xyzzyaabj1),xyzzyaafw1(xyz&
&zyaabj1),xyzzyaagy1(xyzzyaabj1),xyzzyaagz1(xyzzyaabj1),xyzzyaadm1(xyz&
&zyaabj1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','6')
xyzzyaabs1(1:nitot,1:xyzzyaabj1)=0
xyzzyaabp1(1:nitot)=0
if(xyzzyaabg1==1)then
allocate(xyzzyaabt1(nitot,xyzzyaabj1),stat=xyzzyaaaj50)
elseif(xyzzyaabg1==2)then
if(nbasis<=0)call errstop_master('READ_PJASTROW','Number of atoms in b&
&asis is zero.  If this is a Wigner crystal then use "atom" labels in &
&supercell.')
allocate(xyzzyaabt1(nbasis,xyzzyaabj1),stat=xyzzyaaaj50)
else
if(nitype<=0)call errstop_master('READ_PJASTROW','Number of atom speci&
&es is zero.  If this is a Wigner crystal then use "atom" labels in su&
&percell.')
allocate(xyzzyaabt1(nitype,xyzzyaabj1),stat=xyzzyaaaj50)
endif
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','6.1')
do xyzzyaaab50=1,xyzzyaabj1
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='START SET '//trim(i2s(xyzzyaaab50)))exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "START SET '//trim(i2s(xyzzyaaab50))//'".')
enddo
if(am_master)call wout(' SET '//trim(i2s(xyzzyaaab50))//':')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaadm1(xyzzyaaab50)
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabt1(1:xyzzyaadm1(xyzzyaaab50)&
&,xyzzyaaab50)
if(xyzzyaabg1==1)then
xyzzyaadj1(xyzzyaaab50)=xyzzyaadm1(xyzzyaaab50)
xyzzyaabs1(1:xyzzyaadm1(xyzzyaaab50),xyzzyaaab50)=xyzzyaabt1(1:xyzzyaa&
&dm1(xyzzyaaab50),xyzzyaaab50)
elseif(xyzzyaabg1==2)then
if(xyzzyaadm1(xyzzyaaab50)<1.or.xyzzyaadm1(xyzzyaaab50)>nbasis+1-xyzzy&
&aabj1)call errstop_master('READ_PJASTROW','Problematic number of atom&
&s in set.')
if(any(xyzzyaabt1(1:xyzzyaadm1(xyzzyaaab50),xyzzyaaab50)<1).or.any(xyz&
&zyaabt1(1:xyzzyaadm1(xyzzyaaab50),xyzzyaaab50)>nbasis))call errstop_m&
&aster('READ_PJASTROW','Problem with atom labels.')
call atoms_label_pcell(xyzzyaadm1(xyzzyaaab50),xyzzyaabt1(1:xyzzyaadm1&
&(xyzzyaaab50),xyzzyaaab50),xyzzyaadj1(xyzzyaaab50),xyzzyaabs1(:,xyzzy&
&aaab50))
else
if(xyzzyaadm1(xyzzyaaab50)<1.or.xyzzyaadm1(xyzzyaaab50)>nitype+1-xyzzy&
&aabj1)call errstop_master('READ_PJASTROW','Problematic number of spec&
&ies in set.')
if(any(xyzzyaabt1(1:xyzzyaadm1(xyzzyaaab50),xyzzyaaab50)<1).or.any(xyz&
&zyaabt1(1:xyzzyaadm1(xyzzyaaab50),xyzzyaaab50)>nitype))call errstop_m&
&aster('READ_PJASTROW','Problem with species labels.')
call atoms_label_species(xyzzyaadm1(xyzzyaaab50),xyzzyaabt1(1:xyzzyaad&
&m1(xyzzyaaab50),xyzzyaaab50),xyzzyaadj1(xyzzyaaab50),xyzzyaabs1(:,xyz&
&zyaaab50))
endif
if(am_master)then
call wout('  Number of atoms in set           :  '//trim(i2s(xyzzyaadj&
&1(xyzzyaaab50))))
if(xyzzyaadj1(xyzzyaaab50)<1.or.xyzzyaadj1(xyzzyaaab50)>nitot+1-xyzzya&
&abj1)call errstop('READ_PJASTROW','Problematic number of atoms in set&
&.')
call wout('  The atoms are:')
call write_list_int(xyzzyaadj1(xyzzyaaab50),xyzzyaabs1(1:xyzzyaadj1(xy&
&zzyaaab50),xyzzyaaab50),10,4,1)
if(any(xyzzyaabs1(1:xyzzyaadj1(xyzzyaaab50),xyzzyaaab50)<1).or.any(xyz&
&zyaabs1(1:xyzzyaadj1(xyzzyaaab50),xyzzyaaab50)>nitot))call errstop('R&
&EAD_PJASTROW','Problem with atom labels.')
do xyzzyaaac50=1,xyzzyaadj1(xyzzyaaab50)-1
do xyzzyaaad50=xyzzyaaac50+1,xyzzyaadj1(xyzzyaaab50)
if(xyzzyaabs1(xyzzyaaac50,xyzzyaaab50)==xyzzyaabs1(xyzzyaaad50,xyzzyaa&
&ab50))call errstop('READ_PJASTROW','Ion '//trim(i2s(xyzzyaabs1(xyzzya&
&aac50,xyzzyaaab50)))//' appears twice.')
enddo
enddo
endif
do xyzzyaaac50=1,xyzzyaadj1(xyzzyaaab50)
if(xyzzyaabp1(xyzzyaabs1(xyzzyaaac50,xyzzyaaab50))==0)then
xyzzyaabp1(xyzzyaabs1(xyzzyaaac50,xyzzyaaab50))=xyzzyaaab50
else
if(am_master)then
call wout('Ion '//trim(i2s(xyzzyaabs1(xyzzyaaac50,xyzzyaaab50)))//' ap&
&pears in sets '//trim(i2s(xyzzyaabp1(xyzzyaabs1(xyzzyaaac50,xyzzyaaab&
&50))))//' and '//trim(i2s(xyzzyaaab50))//'.')
call errstop('READ_PJASTROW','Stopping.')
endif
endif
enddo
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaagy1(xyzzyaaab50)
if(am_master)then
if(xyzzyaagy1(xyzzyaaab50)==1)then
call wout('  Additional constraints applied to avoid duplicating u.')
elseif(xyzzyaagy1(xyzzyaaab50)==0)then
call wout('  No constraints applied to avoid duplicating u.')
else
call errstop('READ_PJASTROW','Flag for preventing duplication of u sho&
&uld be either 0 or 1.')
endif
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaagz1(xyzzyaaab50)
if(am_master)then
if(xyzzyaagz1(xyzzyaaab50)==1)then
call wout('  Additional constraints applied to avoid duplicating chi.'&
&)
elseif(xyzzyaagz1(xyzzyaaab50)==0)then
call wout('  No constraints applied to avoid duplicating chi.')
else
call errstop('READ_PJASTROW','Flag for preventing duplication of chi s&
&hould be either 0 or 1.')
endif
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacg1(xyzzyaaab50)
if(am_master)then
call wout('  Expansion order (N_f_eN)         :  '//trim(i2s(xyzzyaacg&
&1(xyzzyaaab50))))
if(xyzzyaacg1(xyzzyaaab50)<0)call errstop('READ_PJASTROW','N_f_eN<0')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacf1(xyzzyaaab50)
if(am_master)then
call wout('  Expansion order (N_f_ee)         :  '//trim(i2s(xyzzyaacf&
&1(xyzzyaaab50))))
if(xyzzyaacf1(xyzzyaaab50)<0)call errstop('READ_PJASTROW','N_f_ee<0')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafo1(xyzzyaaab50)
if(am_master)then
call wout('  Spin dependence                       :  '//trim(i2s(xyzz&
&yaafo1(xyzzyaaab50))))
if(xyzzyaafo1(xyzzyaaab50)<-custom_spairs.or.xyzzyaafo1(xyzzyaaab50)>l&
&evels_spairs)call errstop('READ_PJASTROW','Spin dep should be '//trim&
&(i2s(-custom_spairs))//' -- '//trim(i2s(levels_spairs))//' for f.')
endif
xyzzyaafw1(xyzzyaaab50)=no_spairs(xyzzyaafo1(xyzzyaaab50))
if(xyzzyaaab50==1)then
allocate(xyzzyaaae1(0:xyzzyaacg1(xyzzyaaab50),0:xyzzyaacg1(xyzzyaaab50&
&),0:xyzzyaacf1(xyzzyaaab50),xyzzyaafw1(xyzzyaaab50),xyzzyaabj1),xyzzy&
&aaaw1(0:xyzzyaacg1(xyzzyaaab50),0:xyzzyaacg1(xyzzyaaab50),0:xyzzyaacf&
&1(xyzzyaaab50),xyzzyaafw1(xyzzyaaab50),xyzzyaabj1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','7')
xyzzyaaae1(:,:,:,:,:)=0.d0
xyzzyaaaw1(:,:,:,:,:)=1
else
xyzzyaaai50=shape(xyzzyaaae1)
if(xyzzyaacg1(xyzzyaaab50)>xyzzyaaai50(1)-1.or.xyzzyaacf1(xyzzyaaab50)&
&>xyzzyaaai50(3)-1.or.xyzzyaafw1(xyzzyaaab50)>xyzzyaaai50(4))then
allocate(xyzzyaaaz50(0:xyzzyaaai50(1)-1,0:xyzzyaaai50(2)-1,0:xyzzyaaai&
&50(3)-1,xyzzyaaai50(4),xyzzyaaab50-1),xyzzyaaaw50(0:xyzzyaaai50(1)-1,&
&0:xyzzyaaai50(2)-1,0:xyzzyaaai50(3)-1,xyzzyaaai50(4),xyzzyaaab50-1),s&
&tat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','8')
xyzzyaaaz50(0:xyzzyaaai50(1)-1,0:xyzzyaaai50(2)-1,0:xyzzyaaai50(3)-1,1&
&:xyzzyaaai50(4),1:xyzzyaaab50-1)=xyzzyaaae1(0:xyzzyaaai50(1)-1,0:xyzz&
&yaaai50(2)-1,0:xyzzyaaai50(3)-1,1:xyzzyaaai50(4),1:xyzzyaaab50-1)
xyzzyaaaw50(0:xyzzyaaai50(1)-1,0:xyzzyaaai50(2)-1,0:xyzzyaaai50(3)-1,1&
&:xyzzyaaai50(4),1:xyzzyaaab50-1) =xyzzyaaaw1(0:xyzzyaaai50(1)-1,0:xyz&
&zyaaai50(2)-1,0:xyzzyaaai50(3)-1,1:xyzzyaaai50(4),1:xyzzyaaab50-1)
deallocate(xyzzyaaae1,xyzzyaaaw1)
allocate(xyzzyaaae1(0:max(xyzzyaacg1(xyzzyaaab50),xyzzyaaai50(1)-1),0:&
&max(xyzzyaacg1(xyzzyaaab50),xyzzyaaai50(2)-1),0:max(xyzzyaacf1(xyzzya&
&aab50),xyzzyaaai50(3)-1),max(xyzzyaafw1(xyzzyaaab50),xyzzyaaai50(4)),&
&xyzzyaaai50(5)),xyzzyaaaw1(0:max(xyzzyaacg1(xyzzyaaab50),xyzzyaaai50(&
&1)-1),0:max(xyzzyaacg1(xyzzyaaab50),xyzzyaaai50(2)-1),0:max(xyzzyaacf&
&1(xyzzyaaab50),xyzzyaaai50(3)-1),max(xyzzyaafw1(xyzzyaaab50),xyzzyaaa&
&i50(4)),xyzzyaaai50(5)),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','9')
xyzzyaaae1(:,:,:,:,:)=0.d0
xyzzyaaaw1(:,:,:,:,:)=1
xyzzyaaae1(0:xyzzyaaai50(1)-1,0:xyzzyaaai50(2)-1,0:xyzzyaaai50(3)-1,1:&
&xyzzyaaai50(4),1:xyzzyaaab50-1)=xyzzyaaaz50(0:xyzzyaaai50(1)-1,0:xyzz&
&yaaai50(2)-1,0:xyzzyaaai50(3)-1,1:xyzzyaaai50(4),1:xyzzyaaab50-1)
xyzzyaaaw1(0:xyzzyaaai50(1)-1,0:xyzzyaaai50(2)-1,0:xyzzyaaai50(3)-1,1:&
&xyzzyaaai50(4),1:xyzzyaaab50-1)=xyzzyaaaw50(0:xyzzyaaai50(1)-1,0:xyzz&
&yaaai50(2)-1,0:xyzzyaaai50(3)-1,1:xyzzyaaai50(4),1:xyzzyaaab50-1)
deallocate(xyzzyaaaz50,xyzzyaaaw50)
endif
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacv1(xyzzyaaab50),xyzzyaadg1(x&
&yzzyaaab50)
if(abs(xyzzyaacv1(xyzzyaaab50))<xyzzyaahi1)then
xyzzyaacv1(xyzzyaaab50)=xyzzyaaqx1()
if(am_master)call wout('  Using default cutoff length L_f.')
endif
if(am_master)then
call display_param(xyzzyaacv1(xyzzyaaab50),xyzzyaadg1(xyzzyaaab50),'Cu&
&toff')
if(xyzzyaacv1(xyzzyaaab50)<0.d0)call errstop('READ_PJASTROW','L_f<0')
if(isperiodic.and.xyzzyaacv1(xyzzyaaab50)>0.5d0*wigner_seitz_radius)ca&
&ll errstop('READ_PJASTROW','L_f > 0.5 * radius of sphere inscribed in&
& Wigner-Seitz cell.')
if(xyzzyaadg1(xyzzyaaab50)<0.or.xyzzyaadg1(xyzzyaaab50)>1)call errstop&
&('READ_PJASTROW','Optimizable flag should be 0 or 1.')
if(xyzzyaadg1(xyzzyaaab50)==1.and.xyzzyaagx1<2)call errstop('READ_PJAS&
&TROW','Cannot optimize L_f if truncation order < 2.')
endif
read(xyzzyaagw1,*,err=666,end=666)
xyzzyaaak50=xyzzyaaor1(xyzzyaaab50)
allocate(xyzzyaabd50(xyzzyaaak50),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','10')
call xyzzyaaop1(xyzzyaaab50,xyzzyaaak50,xyzzyaabd50,xyzzyaaag50)
if(xyzzyaaag50<1)call errstop_master('READ_PJASTROW','No free gamma pa&
&rameters in this set. Stopping.')
read_params_gamma: do xyzzyaaae50=1,xyzzyaafw1(xyzzyaaab50)
xyzzyaaaf50=0
do n=0,xyzzyaacf1(xyzzyaaab50)
do m=0,xyzzyaacg1(xyzzyaaab50)
do l=m,xyzzyaacg1(xyzzyaaab50)
xyzzyaaaf50=xyzzyaaaf50+1
if(.not.xyzzyaabd50(xyzzyaaaf50))then
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaae1(l,m,n,xyzzyaaae50,xyzz&
&yaaab50),xyzzyaaaw1(l,m,n,xyzzyaaae50,xyzzyaaab50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
xyzzyaaae1(l,m,n,xyzzyaaae50,xyzzyaaab50)=0.d0
xyzzyaaaw1(l,m,n,xyzzyaaae50,xyzzyaaab50)=1
exit read_params_gamma
endif
empty_jastrow=.false.
if(am_master)then
call display_param(xyzzyaaae1(l,m,n,xyzzyaaae50,xyzzyaaab50),xyzzyaaaw&
&1(l,m,n,xyzzyaaae50,xyzzyaaab50),'gamma_'//trim(i2s(l))//','//trim(i2&
&s(m))//','//trim(i2s(n))//','//trim(i2s(xyzzyaaae50))//','//trim(i2s(&
&xyzzyaaab50)))
if(xyzzyaaaw1(l,m,n,xyzzyaaae50,xyzzyaaab50)<0.or.xyzzyaaaw1(l,m,n,xyz&
&zyaaae50,xyzzyaaab50)>1)call errstop('READ_PJASTROW','Optimizable fla&
&g should be 0 or 1.')
endif
endif
enddo
enddo
enddo
enddo read_params_gamma
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END SET '//trim(i2s(xyzzyaaab50))//'".')
if(trim(adjustl(char_80))/='END SET '//trim(i2s(xyzzyaaab50)))call err&
&stop_master('READ_PJASTROW','Was expecting to find "END SET '//trim(i&
&2s(xyzzyaaab50))//'".')
if(am_master)then
call wout('  No. of parameters in set         :  '//trim(i2s(xyzzyaaag&
&50*xyzzyaafw1(xyzzyaaab50))))
call wout('  (In addition to the cutoff length.)')
endif
deallocate(xyzzyaabd50)
enddo
if(any(xyzzyaabp1==0).and.am_master)then
if(allow_nochi_atoms)then
call errwarn('READ_PJASTROW','Some atoms do not belong to any f set.')
else
call errstop('READ_PJASTROW','Some atoms do not belong to any f set.')
endif
endif
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END F TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find "END F TERM".')
enddo
if(am_master)call wout()
allocate(xyzzyaaam1((maxval(xyzzyaacg1)+1)*(maxval(xyzzyaacg1)+1)*(max&
&val(xyzzyaacf1)+1),max_spin_pairs,xyzzyaabj1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','10.3')
xyzzyaaca1=maxval(xyzzyaacg1)*netot
allocate(xyzzyaads1(maxval(xyzzyaacg1),netot,nitot),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','10.5')
xyzzyaads1=0.d0
!$omp parallel default(none) private(xyzzyaaaj50) shared(xyzzyaacg1,xy&
!$omp &zzyaacf1,netot,max_spin_pairs)
allocate(xyzzyaado1(0:maxval(xyzzyaacg1)),xyzzyaadp1(0:maxval(xyzzyaac&
&g1)),xyzzyaadq1(maxval(xyzzyaacg1),netot),xyzzyaadr1(0:maxval(xyzzyaa&
&cf1)),xyzzyaadt1(maxval(xyzzyaacg1)),xyzzyaadu1(2:maxval(xyzzyaacg1))&
&,xyzzyaadv1(maxval(xyzzyaacg1)),xyzzyaadw1(2:maxval(xyzzyaacg1)),xyzz&
&yaadx1(0:maxval(xyzzyaacg1),0:maxval(xyzzyaacf1),max_spin_pairs),stat&
&=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','10.6')
xyzzyaado1(0)=1.d0
xyzzyaadp1(0)=1.d0
xyzzyaadr1(0)=1.d0
xyzzyaadq1=0.d0
!$omp end parallel
xyzzyaady1=0
if(all(xyzzyaacg1(:)==1).and.all(xyzzyaacf1(:)==1))xyzzyaady1=1
if(all(xyzzyaacg1(:)==2).and.all(xyzzyaacf1(:)==2))xyzzyaady1=2
if(all(xyzzyaacg1(:)==3).and.all(xyzzyaacf1(:)==3))xyzzyaady1=3
elseif(trim(adjustl(char_80))=='START P TERM')then
if(xyzzyaaey1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about p terms?')
if(.not.isperiodic)call errstop_master('READ_PJASTROW','Found a p-term&
& in a non-periodic system.')
xyzzyaaey1=.true.
if(am_master)call wout('P term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafj1
if(am_master)then
call wout(' Spin dependence                   :  '//trim(i2s(xyzzyaafj&
&1)))
if(xyzzyaafj1<-custom_spairs.or.xyzzyaafj1>levels_spairs)call errstop(&
&'READ_PJASTROW','Spin dep should be '//trim(i2s(-custom_spairs))//' -&
&- '//trim(i2s(levels_spairs))//' for P.')
if(xyzzyaafj1/=0.and.noncoll_spin)call errstop('READ_PJASTROW','Spin-d&
&ependence of P should be zero for noncollinear-spin calculations.')
endif
xyzzyaafr1=no_spairs(xyzzyaafj1)
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafy1
if(am_master)then
call wout(' Number of G-vectors               :  '//trim(i2s(xyzzyaafy&
&1)))
if(xyzzyaafy1<=0)call errstop('READ_PJASTROW','Need more G-vectors.')
endif
read(xyzzyaagw1,*,err=666,end=666)
allocate(xyzzyaaga1(3,xyzzyaafy1),xyzzyaagb1(xyzzyaafy1),stat=xyzzyaaa&
&j50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','11')
if(am_master)call wout(' G-vector           Label')
do xyzzyaaal50=1,xyzzyaafy1
read(xyzzyaagw1,*,err=666,end=666)xyzzyaaga1(1:3,xyzzyaaal50),xyzzyaag&
&b1(xyzzyaaal50)
if(am_master)then
char_18='('//trim(i2s(xyzzyaaga1(1,xyzzyaaal50)))//','//trim(i2s(xyzzy&
&aaga1(2,xyzzyaaal50)))//','//trim(i2s(xyzzyaaga1(3,xyzzyaaal50)))//')&
&'
call wout('  '//char_18//'  '//trim(i2s(xyzzyaagb1(xyzzyaaal50))))
if(xyzzyaagb1(xyzzyaaal50)<1.or.xyzzyaagb1(xyzzyaaal50)>xyzzyaafy1)cal&
&l errstop('READ_PJASTROW','Label for G-vector should be between 1 and&
& the number of G-vectors.')
if(dimensionality<3.and.xyzzyaaga1(3,xyzzyaaal50)/=0)call errstop('REA&
&D_PJASTROW','Component of G vector in direction of third lattice vect&
&or should be 0 in a '//trim(i2s(dimensionality))//'D system.')
if(dimensionality==1.and.xyzzyaaga1(2,xyzzyaaal50)/=0)call errstop('RE&
&AD_PJASTROW','Component of G vector in direction of second lattice ve&
&ctor should be 0 in a 1D system.')
endif
enddo
if(am_master)then
do xyzzyaaal50=1,xyzzyaafy1
do xyzzyaaam50=1,xyzzyaaal50
if((xyzzyaaga1(1,xyzzyaaal50)+xyzzyaaga1(1,xyzzyaaam50)==0).and.(xyzzy&
&aaga1(2,xyzzyaaal50)+xyzzyaaga1(2,xyzzyaaam50)==0).and.(xyzzyaaga1(3,&
&xyzzyaaal50)+xyzzyaaga1(3,xyzzyaaam50)==0))then
call wout('Both ('//trim(i2s(xyzzyaaga1(1,xyzzyaaal50)))//','//trim(i2&
&s(xyzzyaaga1(2,xyzzyaaal50)))//','//trim(i2s(xyzzyaaga1(3,xyzzyaaal50&
&)))//') and ('//trim(i2s(xyzzyaaga1(1,xyzzyaaam50)))//','//trim(i2s(x&
&yzzyaaga1(2,xyzzyaaam50)))//','//trim(i2s(xyzzyaaga1(3,xyzzyaaam50)))&
&//') are in the list of G-vectors.')
call wout('If G is present then -G should not be present.')
call wout('NB, G=(0,0,0) is not allowed.')
call errstop('READ_PJASTROW','Stopping.')
endif
enddo
enddo
endif
xyzzyaage1=0
do xyzzyaaal50=1,xyzzyaafy1
if(xyzzyaaga1(1,xyzzyaaal50)<0)xyzzyaaga1(1:3,xyzzyaaal50)=-xyzzyaaga1&
&(1:3,xyzzyaaal50)
if(xyzzyaaga1(1,xyzzyaaal50)>xyzzyaage1(1))xyzzyaage1(1)=xyzzyaaga1(1,&
&xyzzyaaal50)
if(abs(xyzzyaaga1(2,xyzzyaaal50))>xyzzyaage1(2))xyzzyaage1(2)=abs(xyzz&
&yaaga1(2,xyzzyaaal50))
if(abs(xyzzyaaga1(3,xyzzyaaal50))>xyzzyaage1(3))xyzzyaage1(3)=abs(xyzz&
&yaaga1(3,xyzzyaaal50))
enddo
allocate(xyzzyaagk1(0:xyzzyaage1(1)),xyzzyaagl1(-xyzzyaage1(2):xyzzyaa&
&ge1(2)),xyzzyaagm1(-xyzzyaage1(3):xyzzyaage1(3)),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','11.5')
xyzzyaagq1=1
do xyzzyaaal50=2,xyzzyaafy1
xyzzyaabb50=.true.
do xyzzyaaam50=1,xyzzyaaal50-1
if(xyzzyaagb1(xyzzyaaam50)==xyzzyaagb1(xyzzyaaal50))xyzzyaabb50=.false&
&.
enddo
if(xyzzyaabb50)xyzzyaagq1=xyzzyaagq1+1
enddo
if(am_master)then
if(xyzzyaagq1>1)then
call wout(' There are '//trim(i2s(xyzzyaagq1))//' independent expansio&
&n coefficients a_A.')
else
call wout(' There is 1 independent expansion coefficient a_A.')
endif
do xyzzyaaal50=1,xyzzyaafy1
if(xyzzyaagb1(xyzzyaaal50)>xyzzyaagq1)call errstop('READ_PJASTROW','G-&
&vector labels should be 1,2,3,... ')
enddo
endif
allocate(xyzzyaaaf1(xyzzyaagq1,xyzzyaafr1),xyzzyaaax1(xyzzyaagq1,xyzzy&
&aafr1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','12')
xyzzyaaaf1(:,:)=0.d0
xyzzyaaax1(:,:)=1
read(xyzzyaagw1,*,err=666,end=666)
read_params_a: do xyzzyaaae50=1,xyzzyaafr1
do xyzzyaaal50=1,xyzzyaagq1
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaaf1(xyzzyaaal50,xyzzyaaae5&
&0),xyzzyaaax1(xyzzyaaal50,xyzzyaaae50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout(' Not all coefficients supplied: rest assumed t&
&o be zero.')
xyzzyaaaf1(xyzzyaaal50,xyzzyaaae50)=0.d0
xyzzyaaax1(xyzzyaaal50,xyzzyaaae50)=1
exit read_params_a
endif
empty_jastrow=.false.
if(am_master)call display_param(xyzzyaaaf1(xyzzyaaal50,xyzzyaaae50),xy&
&zzyaaax1(xyzzyaaal50,xyzzyaaae50),'a_'//trim(i2s(xyzzyaaal50))//','//&
&trim(i2s(xyzzyaaae50)))
enddo
enddo read_params_a
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END P TERM".')
if(trim(adjustl(char_80))/='END P TERM')call errstop_master('READ_PJAS&
&TROW','Was expecting to find "END P TERM".')
if(am_master)then
call wout(' No. of parameters                 :  '//trim(i2s(xyzzyaagq&
&1*xyzzyaafr1)))
call wout()
endif
allocate(xyzzyaagg1(3,xyzzyaafy1),xyzzyaagh1(xyzzyaafy1),xyzzyaagt1(xy&
&zzyaagq1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','13')
allocate(xyzzyaaah1(xyzzyaagq1,max_spin_pairs),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'CONSTRUCT_EXP_P','')
xyzzyaaah1=0.d0
elseif(trim(adjustl(char_80))=='START Q TERM')then
if(xyzzyaaez1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about q terms?')
if(.not.isperiodic)call errstop_master('READ_PJASTROW','Found a q-term&
& in a non-periodic system.')
xyzzyaaez1=.true.
if(am_master)call wout('Q term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafk1
if(am_master)then
call wout(' Spin dependence                   :  '//trim(i2s(xyzzyaafk&
&1)))
if(xyzzyaafk1<-custom_ssingles.or.xyzzyaafk1>levels_ssingles)call errs&
&top('READ_PJASTROW','Spin dep should be '//trim(i2s(-custom_ssingles)&
&)//' -- '//trim(i2s(levels_ssingles))//' for q.')
endif
xyzzyaafs1=no_ssingles(xyzzyaafk1)
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaafz1
if(am_master)then
call wout(' Number of G-vectors               :  '//trim(i2s(xyzzyaafz&
&1)))
if(xyzzyaafz1<=0)call errstop('READ_PJASTROW','Need more G-vectors.')
endif
read(xyzzyaagw1,*,err=666,end=666)
allocate(xyzzyaagc1(3,xyzzyaafz1),xyzzyaagd1(xyzzyaafz1),stat=xyzzyaaa&
&j50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','15')
if(am_master)call wout(' G-vector           Label')
do xyzzyaaal50=1,xyzzyaafz1
read(xyzzyaagw1,*,err=666,end=666)xyzzyaagc1(1:3,xyzzyaaal50),xyzzyaag&
&d1(xyzzyaaal50)
if(am_master)then
char_18='('//trim(i2s(xyzzyaagc1(1,xyzzyaaal50)))//','//trim(i2s(xyzzy&
&aagc1(2,xyzzyaaal50)))//','//trim(i2s(xyzzyaagc1(3,xyzzyaaal50)))//')&
&'
call wout('  '//char_18//'  '//trim(i2s(xyzzyaagd1(xyzzyaaal50))))
if(xyzzyaagd1(xyzzyaaal50)==0.or.abs(xyzzyaagd1(xyzzyaaal50))>xyzzyaaf&
&z1)call errstop('READ_PJASTROW','Label for G-vector should be between&
& 1 and the number of G-vectors (or the negative of another G-vector l&
&abel).')
if(dimensionality<3.and.xyzzyaagc1(3,xyzzyaaal50)/=0)call errstop('REA&
&D_PJASTROW','Component of G vector in direction of third lattice vect&
&or should be 0 in a '//trim(i2s(dimensionality))//'D system.')
if(dimensionality==1.and.xyzzyaagc1(2,xyzzyaaal50)/=0)call errstop('RE&
&AD_PJASTROW','Component of G vector in direction of second lattice ve&
&ctor should be 0 in a 1D system.')
endif
enddo
if(am_master)then
do xyzzyaaal50=1,xyzzyaafz1
do xyzzyaaam50=1,xyzzyaaal50
if((xyzzyaagc1(1,xyzzyaaal50)+xyzzyaagc1(1,xyzzyaaam50)==0).and.(xyzzy&
&aagc1(2,xyzzyaaal50)+xyzzyaagc1(2,xyzzyaaam50)==0).and.(xyzzyaagc1(3,&
&xyzzyaaal50)+xyzzyaagc1(3,xyzzyaaam50)==0))then
call wout('Both ('//trim(i2s(xyzzyaagc1(1,xyzzyaaal50)))//','//trim(i2&
&s(xyzzyaagc1(2,xyzzyaaal50)))//','//trim(i2s(xyzzyaagc1(3,xyzzyaaal50&
&)))//') and ('//trim(i2s(xyzzyaagc1(1,xyzzyaaam50)))//','//trim(i2s(x&
&yzzyaagc1(2,xyzzyaaam50)))//','//trim(i2s(xyzzyaagc1(3,xyzzyaaam50)))&
&//') are in the list of G-vectors.')
call wout('If G is present then -G should not be present.')
call wout('NB, G=(0,0,0) is not allowed.')
call errstop('READ_PJASTROW','Stopping.')
endif
enddo
enddo
endif
xyzzyaagf1=0
do xyzzyaaal50=1,xyzzyaafz1
if(xyzzyaagc1(1,xyzzyaaal50)<0)xyzzyaagc1(1:3,xyzzyaaal50)=-xyzzyaagc1&
&(1:3,xyzzyaaal50)
if(xyzzyaagc1(1,xyzzyaaal50)>xyzzyaagf1(1))xyzzyaagf1(1)=xyzzyaagc1(1,&
&xyzzyaaal50)
if(abs(xyzzyaagc1(2,xyzzyaaal50))>xyzzyaagf1(2))xyzzyaagf1(2)=abs(xyzz&
&yaagc1(2,xyzzyaaal50))
if(abs(xyzzyaagc1(3,xyzzyaaal50))>xyzzyaagf1(3))xyzzyaagf1(3)=abs(xyzz&
&yaagc1(3,xyzzyaaal50))
enddo
allocate(xyzzyaagn1(0:xyzzyaagf1(1)),xyzzyaago1(-xyzzyaagf1(2):xyzzyaa&
&gf1(2)),xyzzyaagp1(-xyzzyaagf1(3):xyzzyaagf1(3)),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','17.5')
xyzzyaagr1=1
do xyzzyaaal50=2,xyzzyaafz1
xyzzyaabb50=.true.
do xyzzyaaam50=1,xyzzyaaal50-1
if(abs(xyzzyaagd1(xyzzyaaam50))==abs(xyzzyaagd1(xyzzyaaal50)))xyzzyaab&
&b50=.false.
enddo
if(xyzzyaabb50)xyzzyaagr1=xyzzyaagr1+1
enddo
if(am_master)then
if(xyzzyaagr1>1)then
call wout(' There are '//trim(i2s(xyzzyaagr1))//' independent expansio&
&n coefficients b_B.')
else
call wout(' There is 1 independent expansion coefficient b_B.')
endif
do xyzzyaaal50=1,xyzzyaafz1
if(abs(xyzzyaagd1(xyzzyaaal50))>xyzzyaagr1)call errstop('READ_PJASTROW&
&','G-vector labels should be 1,2,3,... ')
enddo
do xyzzyaaal50=1,xyzzyaafz1
if(xyzzyaagd1(xyzzyaaal50)<0)then
xyzzyaabc50=.false.
do xyzzyaaam50=1,xyzzyaafz1
if(xyzzyaagd1(xyzzyaaal50)+xyzzyaagd1(xyzzyaaam50)==0)then
xyzzyaabc50=.true.
exit
endif
enddo
if(.not.xyzzyaabc50)call errstop('READ_PJASTROW','Have found the label&
& '//trim(i2s(xyzzyaagd1(xyzzyaaal50)))//' but cannot find '//trim(i2s&
&(-xyzzyaagd1(xyzzyaaal50)))//'.')
endif
enddo
endif
xyzzyaags1=1
do xyzzyaaal50=2,xyzzyaafz1
xyzzyaabb50=.true.
do xyzzyaaam50=1,xyzzyaaal50-1
if(xyzzyaagd1(xyzzyaaam50)==xyzzyaagd1(xyzzyaaal50))xyzzyaabb50=.false&
&.
enddo
if(xyzzyaabb50)xyzzyaags1=xyzzyaags1+1
enddo
allocate(xyzzyaaag1(xyzzyaagr1,xyzzyaafs1),xyzzyaaay1(xyzzyaagr1,xyzzy&
&aafs1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','18')
xyzzyaaag1(:,:)=0.d0
xyzzyaaay1(:,:)=1
read(xyzzyaagw1,*,err=666,end=666)
read_params_b: do xyzzyaaae50=1,xyzzyaafs1
do xyzzyaaal50=1,xyzzyaagr1
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaag1(xyzzyaaal50,xyzzyaaae5&
&0),xyzzyaaay1(xyzzyaaal50,xyzzyaaae50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout(' Not all coefficients supplied: rest assumed t&
&o be zero.')
xyzzyaaag1(xyzzyaaal50,xyzzyaaae50)=0.d0
xyzzyaaay1(xyzzyaaal50,xyzzyaaae50)=1
exit read_params_b
endif
empty_jastrow=.false.
if(am_master)call display_param(xyzzyaaag1(xyzzyaaal50,xyzzyaaae50),xy&
&zzyaaay1(xyzzyaaal50,xyzzyaaae50),'b_'//trim(i2s(xyzzyaaal50))//','//&
&trim(i2s(xyzzyaaae50)))
enddo
enddo read_params_b
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END Q TERM".')
if(trim(adjustl(char_80))/='END Q TERM')call errstop_master('READ_PJAS&
&TROW','Was expecting to find "END Q TERM".')
if(am_master)call wout(' No. of q parameters               :  '//trim(&
&i2s(xyzzyaagr1*xyzzyaafs1)))
if(am_master)call wout()
allocate(xyzzyaagi1(3,xyzzyaafz1),xyzzyaagj1(xyzzyaafz1),xyzzyaagu1(xy&
&zzyaags1),xyzzyaagv1(-xyzzyaagr1:xyzzyaagr1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','21')
xyzzyaagv1(:)=0
allocate(xyzzyaaai1(xyzzyaags1,nspin),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'CONSTRUCT_EXP_Q','')
elseif(trim(adjustl(char_80))=='START BIEX1 TERM'.or.trim(adjustl(char&
&_80))=='START BIEX3 TERM')then
if((xyzzyaafa1.or.xyzzyaafc1).and.am_master)call errstop('READ_PJASTRO&
&W','Two sets of information about biex1 / biex3 terms?')
if(xyzzyaafb1)call errstop_master('READ_PJASTROW','Cannot use BIEX1 an&
&d BIEX2 at the same time.')
if(netot==4)then
xyzzyaafa1=.true.
elseif(netot==2)then
xyzzyaafc1=.true.
else
call errstop_master('READ_PJASTROW','Bug - biexciton terms.')
endif
if(am_master)then
if(xyzzyaafa1)then
call wout('BIEX1 term:')
else
call wout('BIEX3 term:')
endif
endif
if(xyzzyaafg1)call errstop_master('READ_PJASTROW','You cannot have bot&
&h "biex" and "ex2d" terms.')
if(isperiodic)call errstop_master('READ_PJASTROW','Found a biex1/3 ter&
&m in a periodic system.')
if(xyzzyaafa1)then
if(nspin/=4.and.any(nele/=1).and.am_master)call errstop('READ_PJASTROW&
&','Should have one spin-up electron, one spin-down electron, one spin&
&-up hole and one spin-down hole in biexciton system.')
if(abs(pmass(2)-pmass(1))>xyzzyaahi1*abs(pmass(1)).and.am_master)call &
&errstop('READ_PJASTROW','Masses of the two electrons should be the sa&
&me.')
if(abs(pmass(4)-pmass(3))>xyzzyaahi1*abs(pmass(3)).and.am_master)call &
&errstop('READ_PJASTROW','Masses of the two holes should be the same.'&
&)
else
if(nspin/=2.and.any(nele/=1).and.am_master)call errstop('READ_PJASTROW&
&','Should have one spin-up particle and one spin-down particle in eff&
&ective biexciton system.')
if(abs(pmass(2)-pmass(1))>xyzzyaahi1*abs(pmass(1)).and.am_master)call &
&errstop('READ_PJASTROW','Masses of the two particles should be the sa&
&me.')
endif
allocate(xyzzyaaih1(9),xyzzyaaij1(9),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','biex1/3 arrays')
xyzzyaaih1=1.d0
xyzzyaaih1(6)=-1.d0
xyzzyaaih1(8)=-1.d0
xyzzyaaij1=1
empty_jastrow=.false.
if(xyzzyaafa1)then
xyzzyaaih1(1)=pmass(1)*pcharge(1)**2
xyzzyaaij1(1)=-1
xyzzyaaih1(3)=pmass(3)*pcharge(3)**2
xyzzyaaij1(3)=-1
if(heg_nlayers==1)then
xyzzyaaih1(5)=2.d0*pcharge(1)*pcharge(3)*pmass(1)*pmass(3)/(pmass(1)+p&
&mass(3))
xyzzyaaij1(5)=-1
else
xyzzyaaih1(5)=0.d0
xyzzyaaij1(5)=-1
endif
if(abs(pmass(1)-pmass(3))<xyzzyaahi1*abs(pmass(1)))then
xyzzyaaih1(3)=xyzzyaaih1(1)
xyzzyaaij1(3)=-1
xyzzyaaij1(4)=-1
endif
else
if(.not.fix_holes)then
xyzzyaaih1(1)=me_biex3
xyzzyaaij1(1)=-1
xyzzyaaih1(3)=mh_biex3
xyzzyaaij1(3)=-1
else
xyzzyaaih1(1)=pmass(1)
xyzzyaaij1(1)=-1
xyzzyaaih1(3)=0.d0
xyzzyaaij1(3)=-1
xyzzyaaih1(4)=0.d0
xyzzyaaij1(4)=-1
endif
if(heg_nlayers==1)then
if(.not.fix_holes)then
xyzzyaaih1(5)=-2.d0*mu_biex3
xyzzyaaij1(5)=-1
else
xyzzyaaih1(5)=-2.d0*pmass(1)
xyzzyaaij1(5)=-1
endif
else
xyzzyaaih1(5)=0.d0
xyzzyaaij1(5)=-1
endif
if(.not.fix_holes)then
if(abs(me_biex3-mh_biex3)<xyzzyaahi1*abs(me_biex3))then
xyzzyaaih1(3)=xyzzyaaih1(1)
xyzzyaaij1(3)=-1
xyzzyaaij1(4)=-1
endif
endif
endif
read(xyzzyaagw1,*,err=666,end=666)
do xyzzyaaac50=1,9
if(xyzzyaaij1(xyzzyaaac50)>=0)then
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaih1(xyzzyaaac50),xyzzyaaij&
&1(xyzzyaaac50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
xyzzyaaih1(xyzzyaaac50)=0.d0
xyzzyaaij1(xyzzyaaac50)=1
exit
endif
endif
enddo
if(.not.fix_holes)then
if(xyzzyaaij1(4)==-1)xyzzyaaih1(4)=xyzzyaaih1(2)
endif
if(am_master)then
call wout('Parameters in biexciton wave function: ')
do xyzzyaaac50=1,9
call display_param(xyzzyaaih1(xyzzyaaac50),xyzzyaaij1(xyzzyaaac50),'c_&
&'//trim(i2s(xyzzyaaac50)),fmtstring_in='(es20.12)')
if(xyzzyaaij1(xyzzyaaac50)<-1.or.xyzzyaaij1(xyzzyaaac50)>1)call errsto&
&p('READ_PJASTROW','Problem with optimizable flag in biex1.')
enddo
call wout()
if(xyzzyaaih1(1)<0.d0)call errstop('READ_PJASTROW','Should have c_1>=0&
&.')
if(xyzzyaaih1(2)<0.d0)call errstop('READ_PJASTROW','Should have c_2>=0&
&.')
if(xyzzyaaih1(3)<0.d0)call errstop('READ_PJASTROW','Should have c_3>=0&
&.')
if(xyzzyaaih1(4)<0.d0)call errstop('READ_PJASTROW','Should have c_4>=0&
&.')
if(xyzzyaaih1(5)>0.d0)call errstop('READ_PJASTROW','Should have c_5<=0&
&.')
if(xyzzyaaih1(6)>0.d0)call errstop('READ_PJASTROW','Should have c_6<=0&
&.')
if(xyzzyaaih1(7)<0.d0)call errstop('READ_PJASTROW','Should have c_7>=0&
&.')
if(xyzzyaaih1(8)>0.d0)call errstop('READ_PJASTROW','Should have c_8<=0&
&.')
if(xyzzyaaih1(9)<0.d0)call errstop('READ_PJASTROW','Should have c_9>=0&
&.')
endif
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END BIEX1 TERM".')
if(trim(adjustl(char_80))/='END BIEX1 TERM'.and.trim(adjustl(char_80))&
&/='END BIEX3 TERM')call errstop_master('READ_PJASTROW','Was expecting&
& to find "END BIEX1 TERM".')
elseif(trim(adjustl(char_80))=='START BIEX2 TERM')then
if(xyzzyaafb1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about biex2 terms?')
if(xyzzyaafa1)call errstop_master('READ_PJASTROW','Cannot use BIEX1 an&
&d BIEX2 at the same time.')
xyzzyaafb1=.true.
if(am_master)call wout('BIEX2 term:')
if(xyzzyaafg1)call errstop_master('READ_PJASTROW','You cannot have bot&
&h "biex" and "ex2d" terms.')
if(isperiodic)call errstop_master('READ_PJASTROW','Found a biex2 term &
&in a periodic system.')
if(nspin/=4.and.any(nele/=1))call errstop_master('READ_PJASTROW','Shou&
&ld have exactly four particles of different types in order to use a B&
&IEX2 Jastrow term.')
xyzzyaail1=1.d0
xyzzyaaim1=1.d0
xyzzyaaio1=1
xyzzyaaip1=1
empty_jastrow=.false.
read(xyzzyaagw1,*,err=666,end=666)
do xyzzyaaac50=1,2
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaail1(xyzzyaaac50),xyzzyaaio&
&1(xyzzyaaac50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
xyzzyaail1(xyzzyaaac50)=1.d0
xyzzyaaio1(xyzzyaaac50)=1
exit
endif
if(xyzzyaaio1(xyzzyaaac50)/=0.and.xyzzyaaio1(xyzzyaaac50)/=1)call errs&
&top_master('READ_PJASTROW','Problem with optimizable flag in BIEX2 Ja&
&strow term.')
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaim1(xyzzyaaac50),xyzzyaaip&
&1(xyzzyaaac50)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
xyzzyaaim1(xyzzyaaac50)=1.d0
xyzzyaaip1(xyzzyaaac50)=1
exit
endif
if(xyzzyaaip1(xyzzyaaac50)/=0.and.xyzzyaaip1(xyzzyaaac50)/=1)call errs&
&top_master('READ_PJASTROW','Problem with optimizable flag in BIEX2 Ja&
&strow term.')
enddo
if(any(xyzzyaail1<0.d0).or.any(xyzzyaaim1<0.d0))call errstop_master('R&
&EAD_PJASTROW','All parameters in the BIEX2 term should be positive.')
xyzzyaain1(:)=1.d0/xyzzyaail1(:)
if(am_master)then
call wout('Parameters in BIEX2 Jastrow term:')
do xyzzyaaac50=1,2
call display_param(xyzzyaaim1(xyzzyaaac50),xyzzyaaip1(xyzzyaaac50),'b_&
&'//trim(i2s(xyzzyaaac50)),fmtstring_in='(es20.12)')
enddo
call wout()
endif
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END BIEX2 TERM".')
if(trim(adjustl(char_80))/='END BIEX2 TERM')call errstop_master('READ_&
&PJASTROW','Was expecting to find "END BIEX2 TERM".')
elseif(trim(adjustl(char_80))=='START EX2D TERM')then
if(xyzzyaafg1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about ex2d terms?')
xyzzyaafg1=.true.
if(am_master)call wout('EX2D term:')
if(xyzzyaafa1.or.xyzzyaafb1.or.xyzzyaafc1)call errstop_master('READ_PJ&
&ASTROW','You cannot have both "biex" and "ex2d" terms.')
if(isperiodic)call errstop_master('READ_PJASTROW','Found an EX2D term &
&in a periodic system.')
if(any(nele/=1))call errstop_master('READ_PJASTROW','Should have one p&
&article of each species with an EX2D term.')
if(dimensionality/=2)call errstop_master('READ_PJASTROW','The dimensio&
&nality ought to be 2 if you are using an EX2D term.')
xyzzyaaiu1=(trim(int_name)=="logarithmic".or.trim(int_name)=="2D_int")
if(xyzzyaaiu1)then
xyzzyaaiq1=3
if(am_master)call wout(' Logarithmic Kato cusp conditions will be impo&
&sed in the EX2D term.')
else
xyzzyaaiq1=2
if(am_master)call wout(' Coulomb Kato cusp conditions will be imposed &
&in the EX2D term.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaais1
if(am_master)then
call wout(' Pair spin-dependence  : '//trim(i2s(xyzzyaais1)))
if(xyzzyaais1<-custom_spairs.or.xyzzyaais1>levels_spairs)call errstop(&
&'READ_PJASTROW','Pair spin dependence should be between '//trim(i2s(-&
&custom_spairs))//' and '//trim(i2s(levels_spairs))//'.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaait1
if(am_master)then
call wout(' Single spin-dependence: '//trim(i2s(xyzzyaait1)))
if(xyzzyaait1<-custom_ssingles.or.xyzzyaait1>levels_ssingles)call errs&
&top('READ_PJASTROW','Single spin dependence should be between '//trim&
&(i2s(-custom_ssingles))//' and '//trim(i2s(levels_ssingles))//'.')
endif
xyzzyaair1=no_spairs(xyzzyaais1)+nitot*no_ssingles(xyzzyaait1)
allocate(xyzzyaaii1(xyzzyaaiq1,xyzzyaair1),xyzzyaaik1(xyzzyaaiq1,xyzzy&
&aair1),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','ex2d arrays')
xyzzyaaii1=0.d0
if(xyzzyaaiu1)then
xyzzyaaii1(1,:)=0.13d0
xyzzyaaii1(2,:)=-0.37d0
xyzzyaaii1(3,:)=0.62d0
else
xyzzyaaii1(1,:)=-0.41d0
xyzzyaaii1(2,:)=0.78d0
endif
xyzzyaaik1(:,1:no_spairs(xyzzyaais1))=-1
do xyzzyaaac50=1,netot-1
do xyzzyaaad50=xyzzyaaac50+1,netot
xyzzyaaik1(:,which_spair(xyzzyaaac50,xyzzyaaad50,xyzzyaais1))=1
enddo
enddo
xyzzyaaik1(:,no_spairs(xyzzyaais1)+1:)=1
empty_jastrow=.false.
read(xyzzyaagw1,*,err=666,end=666)
do xyzzyaaae50=1,xyzzyaair1
do xyzzyaaac50=1,xyzzyaaiq1
if(xyzzyaaik1(xyzzyaaac50,xyzzyaaae50)>=0)then
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaax50,xyzzyaaau50
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
exit
endif
xyzzyaaii1(xyzzyaaac50,xyzzyaaae50)=xyzzyaaax50
xyzzyaaik1(xyzzyaaac50,xyzzyaaae50)=xyzzyaaau50
endif
enddo
enddo
if(am_master)then
call wout(' Parameters in EX2D wave function: ')
do xyzzyaaae50=1,xyzzyaair1
do xyzzyaaac50=1,xyzzyaaiq1
if(xyzzyaaik1(xyzzyaaac50,xyzzyaaae50)>=0)then
call display_param(xyzzyaaii1(xyzzyaaac50,xyzzyaaae50),xyzzyaaik1(xyzz&
&yaaac50,xyzzyaaae50),'c_'//trim(i2s(xyzzyaaac50))//','//trim(i2s(xyzz&
&yaaae50)),fmtstring_in='(es20.12)')
if(xyzzyaaik1(xyzzyaaac50,xyzzyaaae50)<0.or.xyzzyaaik1(xyzzyaaac50,xyz&
&zyaaae50)>1)call errstop('READ_PJASTROW','Problem with optimizable fl&
&ag in ex2d.')
endif
enddo
if(xyzzyaaiu1)then
if(xyzzyaaii1(2,xyzzyaaae50)>0.d0)call errstop('READ_PJASTROW','Should&
& have c_2<=0.')
if(xyzzyaaii1(3,xyzzyaaae50)<0.d0)call errstop('READ_PJASTROW','Should&
& have c_3>=0.')
else
if(xyzzyaaii1(1,xyzzyaaae50)>0.d0)call errstop('READ_PJASTROW','Should&
& have c_1<=0.')
if(xyzzyaaii1(2,xyzzyaaae50)<0.d0)call errstop('READ_PJASTROW','Should&
& have c_2>=0.')
endif
enddo
call wout()
endif
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END EX2D TERM".')
if(trim(adjustl(char_80))/='END EX2D TERM'.and.trim(adjustl(char_80))/&
&='END EX2D TERM')call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END EX2D TERM".')
elseif(trim(adjustl(char_80))=='START D TERM')then
if(xyzzyaaff1)call errstop_master('READ_PJASTROW','Two sets of informa&
&tion about D terms?')
if(nitot<2)call errstop_master('READ_PJASTROW','No D term needed in a &
&system with less than two atoms.')
xyzzyaaff1=.true.
if(am_master)call wout('D term:')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabm1
if(am_master)call wout(' Number of pairs           :  '//trim(i2s(xyzz&
&yaabm1)))
if(xyzzyaabm1<1.or.xyzzyaabm1>(nitot*(nitot-1))/2)call errstop_master(&
&'READ_PJASTROW','Problematic number of ion pairs in D term.')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabn1
if(am_master)call wout(' Number of sets            :  '//trim(i2s(xyzz&
&yaabn1)))
if(xyzzyaabn1<1.or.xyzzyaabn1>xyzzyaabm1)call errstop_master('READ_PJA&
&STROW','Problematic number of sets of ion pairs in D term.')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacd1
if(am_master)call wout('  Expansion order (N_d)    :  '//trim(i2s(xyzz&
&yaacd1)))
if(xyzzyaacd1<0)call errstop_master('READ_PJASTROW','N_d<0.')
allocate(xyzzyaabu1(xyzzyaabm1,2),xyzzyaabv1(xyzzyaabm1),xyzzyaadk1(xy&
&zzyaabn1),xyzzyaabw1(xyzzyaabn1),xyzzyaaas1(xyzzyaabn1,2,0:xyzzyaacd1&
&,0:xyzzyaacd1),xyzzyaacw1(xyzzyaabn1,2,2),xyzzyaabe1(xyzzyaabn1,2,0:x&
&yzzyaacd1,0:xyzzyaacd1),xyzzyaadh1(xyzzyaabn1,2,2),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','2')
xyzzyaaar50=0
do xyzzyaaab50=1,xyzzyaabn1
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='START SET '//trim(i2s(xyzzyaaab50)))exit
if(trim(adjustl(char_80))/='')call errstop_master('READ_PJASTROW','Was&
& expecting to find "START SET '//trim(i2s(xyzzyaaab50))//'".')
enddo
if(am_master)call wout(' SET '//trim(i2s(xyzzyaaab50))//':')
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabw1(xyzzyaaab50)
if(xyzzyaabw1(xyzzyaaab50)==0)then
if(am_master)call wout('  Pairs in set are asymmetric')
elseif(xyzzyaabw1(xyzzyaaab50)==1)then
if(am_master)call wout('  Pairs in set are symmetric')
else
call errstop_master('READ_PJASTROW','Pair symmetry flag should be 0 or&
& 1.')
endif
read(xyzzyaagw1,*,err=666,end=666)
read(xyzzyaagw1,*,err=666,end=666)xyzzyaadk1(xyzzyaaab50)
if(am_master)call wout('  Number of atom pairs in set   :  '//trim(i2s&
&(xyzzyaadk1(xyzzyaaab50))))
if(xyzzyaadk1(xyzzyaaab50)<1.or.xyzzyaadk1(xyzzyaaab50)>nitot+1-xyzzya&
&abn1)call errstop_master('READ_PJASTROW','Problematic number of atom &
&pairs in the set.')
if(xyzzyaaar50+xyzzyaadk1(xyzzyaaab50)>nitot)call errstop_master('READ&
&_PJASTROW','Number of atom pairs in the set too large.')
read(xyzzyaagw1,*,err=666,end=666)
do xyzzyaaaf50=1,xyzzyaadk1(xyzzyaaab50)
xyzzyaaar50=xyzzyaaar50+1
read(xyzzyaagw1,*,err=666,end=666)xyzzyaabu1(xyzzyaaar50,1:2)
xyzzyaabv1(xyzzyaaar50)=xyzzyaaab50
if(xyzzyaabu1(xyzzyaaar50,1)<1.or.xyzzyaabu1(xyzzyaaar50,1)>nitot.or.x&
&yzzyaabu1(xyzzyaaar50,2)<1.or.xyzzyaabu1(xyzzyaaar50,2)>nitot)call er&
&rstop_master('READ_PJASTROW','Illegal atom index')
if(xyzzyaabu1(xyzzyaaar50,1)==xyzzyaabu1(xyzzyaaar50,2))call errstop_m&
&aster('READ_PJASTROW','Pair of identical atoms ('//trim(i2s(xyzzyaabu&
&1(xyzzyaaat50,1)))//','//trim(i2s(xyzzyaabu1(xyzzyaaat50,2)))//')')
do xyzzyaaat50=1,xyzzyaaar50-1
if(xyzzyaabu1(xyzzyaaat50,1)==xyzzyaabu1(xyzzyaaar50,1).and.xyzzyaabu1&
&(xyzzyaaat50,2)==xyzzyaabu1(xyzzyaaar50,2))call errstop_master('READ_&
&PJASTROW','Pair of atoms ('//trim(i2s(xyzzyaabu1(xyzzyaaat50,1)))//',&
&'//trim(i2s(xyzzyaabu1(xyzzyaaat50,2)))//'appears twice.')
if(xyzzyaabu1(xyzzyaaat50,1)==xyzzyaabu1(xyzzyaaar50,2).and.xyzzyaabu1&
&(xyzzyaaat50,2)==xyzzyaabu1(xyzzyaaar50,1))call errstop_master('READ_&
&PJASTROW','Pair of atoms ('//trim(i2s(xyzzyaabu1(xyzzyaaat50,1)))//',&
&'//trim(i2s(xyzzyaabu1(xyzzyaaat50,2)))//') appears twice in reversed&
& order.')
enddo
enddo
if(xyzzyaaar50<xyzzyaabm1)call errstop_master('READ_PJASTROW','Too few&
& pairs of atoms specified')
if(am_master)then
call wout('  The atom pairs are:')
do xyzzyaaaf50=xyzzyaaar50-xyzzyaadk1(xyzzyaaab50)+1,xyzzyaaar50
call wout('    '//trim(i2s(xyzzyaabu1(xyzzyaaaf50,1)))//' '//trim(i2s(&
&xyzzyaabu1(xyzzyaaaf50,2))))
enddo
endif
xyzzyaacp1=0.d0
read(xyzzyaagw1,*,err=666,end=666)
do xyzzyaaas50=1,2
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,1&
&),xyzzyaadh1(xyzzyaaab50,xyzzyaaas50,1)
if(xyzzyaabw1(xyzzyaaab50)==1)then
xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,2)=xyzzyaacw1(xyzzyaaab50,xyzzyaaas&
&50,1)
xyzzyaadh1(xyzzyaaab50,xyzzyaaas50,2)=xyzzyaadh1(xyzzyaaab50,xyzzyaaas&
&50,1)
else
read(xyzzyaagw1,*,err=666,end=666)xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,2&
&),xyzzyaadh1(xyzzyaaab50,xyzzyaaas50,2)
endif
if(xyzzyaadh1(xyzzyaaab50,xyzzyaaas50,1)<0.or.xyzzyaadh1(xyzzyaaab50,x&
&yzzyaaas50,1)>1)call errstop_master('READ_PJASTROW','Optimizable flag&
& should be 0 or 1.')
if(xyzzyaadh1(xyzzyaaab50,xyzzyaaas50,2)<0.or.xyzzyaadh1(xyzzyaaab50,x&
&yzzyaaas50,2)>1)call errstop_master('READ_PJASTROW','Optimizable flag&
& should be 0 or 1.')
if(abs(xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,1))<xyzzyaahi1.or.abs(xyzzya&
&acw1(xyzzyaaab50,xyzzyaaas50,2))<xyzzyaahi1)then
xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,1)=0.d0
xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,2)=0.d0
xyzzyaadh1(xyzzyaaab50,xyzzyaaas50,1)=0
xyzzyaadh1(xyzzyaaab50,xyzzyaaas50,2)=0
if(am_master)call wout('D_'//xyzzyaamc1(xyzzyaaas50)//' term deactivat&
&ed by zero cutoff')
elseif(am_master)then
call display_param(xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,1),xyzzyaadh1(xy&
&zzyaaab50,xyzzyaaas50,1),'Cutoff_'//xyzzyaamc1(xyzzyaaas50)//'_(I,J)'&
&)
call display_param(xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,2),xyzzyaadh1(xy&
&zzyaaab50,xyzzyaaas50,2),'Cutoff_'//xyzzyaamc1(xyzzyaaas50)//'_(J,I)'&
&)
endif
if(xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,1)<0.d0)call errstop_master('REA&
&D_PJASTROW','L_D_I,J<0')
if(xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,2)<0.d0)call errstop_master('REA&
&D_PJASTROW','L_D_J,I<0')
if(xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,1)>xyzzyaacp1)xyzzyaacp1=xyzzyaa&
&cw1(xyzzyaaab50,xyzzyaaas50,1)
if(xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,2)>xyzzyaacp1)xyzzyaacp1=xyzzyaa&
&cw1(xyzzyaaab50,xyzzyaaas50,2)
do xyzzyaaaf50=xyzzyaaar50-xyzzyaadk1(xyzzyaaab50)+1,xyzzyaaar50
if(xyzzyaacw1(xyzzyaaab50,xyzzyaaas50,1)+xyzzyaacw1(xyzzyaaab50,xyzzya&
&aas50,2)>rionion(4,xyzzyaabu1(xyzzyaaaf50,1),xyzzyaabu1(xyzzyaaaf50,2&
&)))call errstop_master('READ_PJASTROW','sum of cutoff lengths larger &
&than ion distance.')
enddo
enddo
xyzzyaaas1(xyzzyaaab50,:,:,:)=0.d0
xyzzyaabe1(xyzzyaaab50,:,:,:)=1
read(xyzzyaagw1,*,err=666,end=666)
read_d_param: do xyzzyaaas50=1,2
do m=0,xyzzyaacd1
do n=xyzzyaabw1(xyzzyaaab50)*m,xyzzyaacd1
read(xyzzyaagw1,*,iostat=xyzzyaaaa50)xyzzyaaas1(xyzzyaaab50,xyzzyaaas5&
&0,m,n),xyzzyaabe1(xyzzyaaab50,xyzzyaaas50,m,n)
if(xyzzyaaaa50/=0)then
backspace xyzzyaagw1
if(am_master)call wout('  Not all coefficients supplied: rest assumed &
&to be zero.')
exit read_d_param
endif
empty_jastrow=.false.
if(xyzzyaabw1(xyzzyaaab50)==1.and.n>m)then
xyzzyaaas1(xyzzyaaab50,xyzzyaaas50,n,m)=xyzzyaaas1(xyzzyaaab50,xyzzyaa&
&as50,m,n)
xyzzyaabe1(xyzzyaaab50,xyzzyaaas50,n,m)=xyzzyaabe1(xyzzyaaab50,xyzzyaa&
&as50,m,n)
endif
if(xyzzyaabe1(xyzzyaaab50,xyzzyaaas50,m,n)<0.or.xyzzyaabe1(xyzzyaaab50&
&,xyzzyaaas50,m,n)>1)call errstop_master('READ_PJASTROW','Optimizable &
&flag should be 0 or 1.')
if(am_master)call display_param(xyzzyaaas1(xyzzyaaab50,xyzzyaaas50,m,n&
&),xyzzyaabe1(xyzzyaaab50,xyzzyaaas50,m,n),'D_'//xyzzyaamc1(xyzzyaaas5&
&0)//'_('//trim(i2s(xyzzyaaab50))//')'//trim(i2s(m))//','//trim(i2s(n)&
&))
enddo
enddo
enddo read_d_param
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50/=0)call errstop_master('READ_PJASTROW','Was expecting t&
&o find "END SET '//trim(i2s(xyzzyaaab50))//'".')
if(trim(adjustl(char_80))/='END SET '//trim(i2s(xyzzyaaab50)))call err&
&stop_master('READ_PJASTROW','Was expecting to find "END SET '//trim(i&
&2s(xyzzyaaab50))//'".')
if(am_master)then
call wout('  No. of parameters in set :  '//trim(i2s((xyzzyaacd1+1)*xy&
&zzyaacd1/(xyzzyaabw1(xyzzyaaab50)+1)+(xyzzyaacd1+1))))
call wout('  (In addition to the cutoff length.)')
endif
enddo
do
read(xyzzyaagw1,'(a)',err=666,end=666)char_80
if(trim(adjustl(char_80))=='END D TERM')exit
if(trim(adjustl(char_80))/=''.and.am_master)call errstop('READ_PJASTRO&
&W','Was expecting to find END D TERM.')
enddo
if(am_master)call wout()
elseif(trim(adjustl(char_80))=='END JASTROW')then
exit
elseif(trim(adjustl(char_80))/="")then
if(am_master)call errstop('READ_PJASTROW','Expecting either "END JASTR&
&OW" or e.g. "START U TERM" in correlation.data.')
endif
enddo
do
read(xyzzyaagw1,'(a)',iostat=xyzzyaaaa50)char_80
if(xyzzyaaaa50<0)exit
if(xyzzyaaaa50>0)call errstop('READ_PJASTROW','Problem reading correla&
&tion.data. Please check this file.')
if(trim(adjustl(char_80))=='START JASTROW')call errstop('READ_PJASTROW&
&','There seems to be more than one Jastrow factor in correlation.data&
&.')
enddo
close(xyzzyaagw1)
if(am_master)then
call wout('Finished reading Jastrow factor from correlation.data.')
call wout()
endif
if(empty_jastrow.and.xyzzyaaeq1.and.all(pcharge/=0.d0).and.(prefer_sho&
&rt_cusp.or.homogeneous_system.or.periodicity==3))then
call xyzzyaapc1
if(any(xyzzyaaab1(0,:)/=0.d0))then
if(am_master)then
call wout('Have inserted default values for alpha_0 in the u term:')
do xyzzyaaae50=1,xyzzyaafp1
call display_param(xyzzyaaab1(0,xyzzyaaae50),xyzzyaaat1(0,xyzzyaaae50)&
&,'alpha_0,'//trim(i2s(xyzzyaaae50)))
enddo
call wout()
endif
endif
endif
if(hard_sphere)then
if(.not.xyzzyaaep1.and.am_master)call errstop('READ_PJASTROW','Need to&
& have u term if hard-sphere interaction is used.')
xyzzyaanz1=wigner_seitz_radius
if(hard_diam<=0.d0)call errstop('READ_PJASTROW','Must have hard diamet&
&er D>0.')
xyzzyaaoa1=1.d0/hard_diam
xyzzyaaob1=1.d0/xyzzyaanz1
xyzzyaaoc1=xyzzyaaoa1-xyzzyaaob1
endif
xyzzyaafd1=xyzzyaafa1.or.xyzzyaafb1.or.xyzzyaafc1
xyzzyaafe1=xyzzyaaeu1.or.xyzzyaaev1
have_jastrow3=xyzzyaafe1.or.xyzzyaafd1.or.xyzzyaafg1
if(xyzzyaaev1.or.xyzzyaafd1.or.xyzzyaafg1)then
allocate(xyzzyaaig1(4,netot,netot),stat=xyzzyaaaj50)
call check_alloc(xyzzyaaaj50,'READ_PJASTROW','eevecs_precomp')
endif
if(xyzzyaaex1)then
call xyzzyaaok1
call xyzzyaaon1
if(am_master)then
call wout('Imposed symmetry and no-cusp constraints on gamma array.')
if(.not.xyzzyaaot1(.true.))call errstop('READ_PJASTROW','Parameters of&
& F term do not satisfy gamma constraints.')
call wout('Checked that gamma array satisfies its constraints.')
endif
endif
if(xyzzyaaev1)then
call xyzzyaaou1
call xyzzyaaox1
if(am_master)then
call wout('Imposed symmetry and no-cusp constraints on H parameters.')
if(.not.xyzzyaapb1(.true.))call errstop('READ_PJASTROW','Parameters of&
& H term do not satisfy their constraints.')
call wout('Checked that H parameters satisfy their constraints.')
endif
endif
if(xyzzyaaey1)call xyzzyaapj1
if(xyzzyaaez1)call xyzzyaapk1
if(am_master.and.(xyzzyaaey1.or.xyzzyaaez1))call wout('G-vector arrays&
& set up.')
if(xyzzyaaep1)call xyzzyaapd1
if(xyzzyaaes1)call xyzzyaape1
if(xyzzyaaeu1)call xyzzyaapf1
if(xyzzyaaew1)call xyzzyaapg1
if(xyzzyaaex1)call xyzzyaaph1
if(am_master.and.(xyzzyaaeq1.or.xyzzyaaes1.or.xyzzyaaeu1.or.xyzzyaaew1&
&.or.xyzzyaaex1))call wout('Polynomials constructed.')
if(xyzzyaaey1)call xyzzyaapl1
if(xyzzyaaez1)call xyzzyaapm1
if(am_master.and.(xyzzyaaey1.or.xyzzyaaez1))call wout('Plane-wave expa&
&nsion coefficients constructed.')
call xyzzyaaqz1
if(am_master)then
if(xyzzyaaep1.or.xyzzyaaew1.or.xyzzyaaex1)then
call wout('Checked that cusp and cutoff conditions are satisfied.')
endif
call wout()
endif
if(gen_gjastrow)call xyzzyaaoj1
if(xyzzyaaho1)call xyzzyaara1(xyzzyaaep1,xyzzyaaes1,xyzzyaaeu1,xyzzyaa&
&ew1,xyzzyaaex1,xyzzyaaey1,xyzzyaaez1,xyzzyaahj1,xyzzyaahk1)
if(am_master)then
call wout('Finished Jastrow setup.')
call wout()
endif
if(.not.empty_jastrow.and.finite_size_corr)call finite_size_corr_ke_pj&
&astrow
return
666 if(am_master)call errstop('READ_PJASTROW','Problem reading correla&
&tion.data. Please check this file.')
end subroutine read_pjastrow
subroutine xyzzyaaoj1
use slaarnaan, only : ee_kato_gamma
implicit none
integer xyzzyaaaa51,xyzzyaaab51,xyzzyaaac51,xyzzyaaad51,xyzzyaaae51,xy&
&zzyaaaf51,xyzzyaaag51,xyzzyaaah51,xyzzyaaai51,xyzzyaaaj51,xyzzyaaak51&
&(2),xyzzyaaal51(3),xyzzyaaam51,xyzzyaaan51,xyzzyaaao51,xyzzyaaap51,xy&
&zzyaaaq51,xyzzyaaar51,xyzzyaaas51,xyzzyaaat51,xyzzyaaau51,xyzzyaaav51&
&,xyzzyaaaw51,xyzzyaaax51,xyzzyaaay51(3),xyzzyaaaz51,xyzzyaaba51,xyzzy&
&aabb51,xyzzyaabc51,xyzzyaabd51,xyzzyaabe51,xyzzyaabf51,xyzzyaabg51
integer,allocatable :: xyzzyaabh51(:,:),xyzzyaabi51(:,:),xyzzyaabj51(:&
&)
logical xyzzyaabk51,xyzzyaabl51
character(8) f_zero
character(80) char_80,char_80_1,char_80_2,char_80_3
character(20),parameter :: xyzzyaabm51(-2:1)=(/'determined ','determin&
&ed ','fixed      ','optimizable'/)
if(.not.am_master)return
call wout('Writing GJASTROW conversion to "parameters_converted.casl".&
&')
call open_units(xyzzyaaaa51,xyzzyaaab51)
if(xyzzyaaab51/=0)call errstop('CONVERT_TO_GJASTROW','Unable to find f&
&ree i/o unit.')
open(unit=xyzzyaaaa51,file='parameters_converted.casl',status='replace&
&',iostat=xyzzyaaab51)
if(xyzzyaaab51/=0)call errstop('CONVERT_TO_GJASTROW','Problem opening &
&parameters_converted.casl.')
write(xyzzyaaaa51,'(a)')'JASTROW:'
write(xyzzyaaaa51,*)' Title: automatic conversion by CASINO'
xyzzyaaac51=0
if(xyzzyaaep1.and.xyzzyaaeq1)then
xyzzyaaac51=xyzzyaaac51+1
write(xyzzyaaaa51,*)' TERM '//trim(i2s(xyzzyaaac51))//':'
write(xyzzyaaaa51,*)'   Rank: [ 2, 0 ]'
write(xyzzyaaaa51,*)'   e-e cusp: T'
if(maxval(which_spair(:,:,xyzzyaafh1))==(nspin*(nspin+1))/2)then
write(xyzzyaaaa51,*)'   Rules: [ ]'
else
write(xyzzyaaaa51,*)'   Rules:'
do xyzzyaaad51=1,no_spairs(xyzzyaafh1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafh1),which_spair(:,:,xyzzyaa&
&fh1)==xyzzyaaad51)
xyzzyaaam51=minval(xyzzyaaak51)
xyzzyaaan51=maxval(xyzzyaaak51)
do xyzzyaaao51=xyzzyaaam51,nspin
xyzzyaaaq51=xyzzyaaao51
if(xyzzyaaao51==xyzzyaaam51)xyzzyaaaq51=xyzzyaaan51+1
do xyzzyaaap51=xyzzyaaaq51,nspin
if(nele(xyzzyaaam51)<1.or.nele(xyzzyaaan51)<1.or.xyzzyaaam51==xyzzyaaa&
&n51.and.nele(xyzzyaaam51)<2)then
xyzzyaaam51=xyzzyaaao51
xyzzyaaan51=xyzzyaaap51
cycle
endif
if(nele(xyzzyaaao51)<1.or.nele(xyzzyaaap51)<1.or.xyzzyaaao51==xyzzyaaa&
&p51.and.nele(xyzzyaaao51)<2)cycle
if(which_spair(xyzzyaaao51,xyzzyaaap51,xyzzyaafh1)==which_spair(xyzzya&
&aam51,xyzzyaaan51,xyzzyaafh1))write(xyzzyaaaa51,*)'     ',trim(i2s(xy&
&zzyaaam51)),'-',trim(i2s(xyzzyaaan51)),'=',trim(i2s(xyzzyaaao51)),'-'&
&,trim(i2s(xyzzyaaap51))
enddo
enddo
enddo
endif
write(xyzzyaaaa51,*)'   e-e basis: [ Type: natural power, Order: '//tr&
&im(i2s(xyzzyaabx1+1))//' ]'
write(xyzzyaaaa51,*)'   e-e cutoff:'
write(xyzzyaaaa51,*)'     Type: alt polynomial'
write(xyzzyaaaa51,*)'     Constants: [ C: '//trim(i2s(xyzzyaagx1))//' &
&]'
write(xyzzyaaaa51,*)'     Parameters:'
write(char_80,*)xyzzyaach1
do xyzzyaaad51=1,no_spairs(xyzzyaafh1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafh1),which_spair(:,:,xyzzyaa&
&fh1)==xyzzyaaad51)
if(nele(xyzzyaaak51(1))<1.or.nele(xyzzyaaak51(2))<1.or.(xyzzyaaak51(1)&
&==xyzzyaaak51(2).and.nele(xyzzyaaak51(1))<2))cycle
write(xyzzyaaaa51,*)'       Channel '//trim(i2s(minval(xyzzyaaak51)))/&
&/'-'//trim(i2s(maxval(xyzzyaaak51)))//': [ L: [ ',trim(adjustl(char_8&
&0)),', '//trim(xyzzyaabm51(xyzzyaacx1))//' ] ]'
enddo
write(xyzzyaaaa51,*)'   Linear parameters:'
do xyzzyaaad51=1,no_spairs(xyzzyaafh1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafh1),which_spair(:,:,xyzzyaa&
&fh1)==xyzzyaaad51)
if(nele(xyzzyaaak51(1))<1.or.nele(xyzzyaaak51(2))<1.or.(xyzzyaaak51(1)&
&==xyzzyaaak51(2).and.nele(xyzzyaaak51(1))<2))cycle
write(xyzzyaaaa51,*)'     Channel '//trim(i2s(minval(xyzzyaaak51)))//'&
&-'//trim(i2s(maxval(xyzzyaaak51)))//':'
xyzzyaaae51=which_spair(xyzzyaaak51(1),xyzzyaaak51(2),levels_spairs)
if(xyzzyaafh1==0.and.nspin==2)xyzzyaaae51=which_spair(1,2,levels_spair&
&s)
do xyzzyaaaf51=0,xyzzyaabx1
write(char_80,*)xyzzyaaaj1(xyzzyaaaf51,xyzzyaaae51)
write(xyzzyaaaa51,*)'       c_'//trim(i2s(xyzzyaaaf51+1))//': [ ',trim&
&(adjustl(char_80)),', '//trim(xyzzyaabm51(xyzzyaaat1(xyzzyaaaf51,xyzz&
&yaaad51)))//' ]'
enddo
enddo
call wout('Successfully converted U term.')
elseif(xyzzyaaep1)then
xyzzyaaac51=xyzzyaaac51+1
write(xyzzyaaaa51,*)' TERM '//trim(i2s(xyzzyaaac51))//':'
write(xyzzyaaaa51,*)'   Rank: [ 2, 0 ]'
write(xyzzyaaaa51,*)'   e-e cusp: T'
if(maxval(which_spair(:,:,xyzzyaafh1))==(nspin*(nspin+1))/2)then
write(xyzzyaaaa51,*)'   Rules: [ ]'
else
write(xyzzyaaaa51,*)'   Rules:'
do xyzzyaaad51=1,no_spairs(xyzzyaafh1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafh1),which_spair(:,:,xyzzyaa&
&fh1)==xyzzyaaad51)
xyzzyaaam51=minval(xyzzyaaak51)
xyzzyaaan51=maxval(xyzzyaaak51)
do xyzzyaaao51=xyzzyaaam51,nspin
xyzzyaaaq51=xyzzyaaao51
if(xyzzyaaao51==xyzzyaaam51)xyzzyaaaq51=xyzzyaaan51+1
do xyzzyaaap51=xyzzyaaaq51,nspin
if(nele(xyzzyaaam51)<1.or.nele(xyzzyaaan51)<1.or.xyzzyaaam51==xyzzyaaa&
&n51.and.nele(xyzzyaaam51)<2)then
xyzzyaaam51=xyzzyaaao51
xyzzyaaan51=xyzzyaaap51
cycle
endif
if(nele(xyzzyaaao51)<1.or.nele(xyzzyaaap51)<1.or.xyzzyaaao51==xyzzyaaa&
&p51.and.nele(xyzzyaaao51)<2)cycle
if(which_spair(xyzzyaaao51,xyzzyaaap51,xyzzyaafh1)==which_spair(xyzzya&
&aam51,xyzzyaaan51,xyzzyaafh1))write(xyzzyaaaa51,*)'     ',trim(i2s(xy&
&zzyaaam51)),'-',trim(i2s(xyzzyaaan51)),'=',trim(i2s(xyzzyaaao51)),'-'&
&,trim(i2s(xyzzyaaap51))
enddo
enddo
enddo
endif
write(xyzzyaaaa51,*)'   e-e basis:'
write(xyzzyaaaa51,*)'     Type: RPA'
write(xyzzyaaaa51,*)'     Order: 1'
write(xyzzyaaaa51,*)'     Parameters:'
write(char_80,*)xyzzyaach1
do xyzzyaaad51=1,no_spairs(xyzzyaafh1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafh1),which_spair(:,:,xyzzyaa&
&fh1)==xyzzyaaad51)
if(nele(xyzzyaaak51(1))<1.or.nele(xyzzyaaak51(2))<1.or.(xyzzyaaak51(1)&
&==xyzzyaaak51(2).and.nele(xyzzyaaak51(1))<2))cycle
write(char_80,*)xyzzyaaao1(xyzzyaaad51)
write(xyzzyaaaa51,*)'       Channel '//trim(i2s(minval(xyzzyaaak51)))/&
&/'-'//trim(i2s(maxval(xyzzyaaak51)))//': [ F: [ ',trim(adjustl(char_8&
&0)),', '//trim(xyzzyaabm51(xyzzyaaba1(xyzzyaaad51)))//' ] ]'
enddo
write(xyzzyaaaa51,*)'   e-e cutoff:'
write(xyzzyaaaa51,*)'     Type: polynomial'
write(xyzzyaaaa51,*)'     Constants: [ C: '//trim(i2s(xyzzyaagx1))//' &
&]'
write(xyzzyaaaa51,*)'     Parameters:'
write(char_80,*)xyzzyaach1
do xyzzyaaad51=1,no_spairs(xyzzyaafh1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafh1),which_spair(:,:,xyzzyaa&
&fh1)==xyzzyaaad51)
write(xyzzyaaaa51,*)'       Channel '//trim(i2s(minval(xyzzyaaak51)))/&
&/'-'//trim(i2s(maxval(xyzzyaaak51)))//': [ L: [ ',trim(adjustl(char_8&
&0)),', '//trim(xyzzyaabm51(xyzzyaacx1))//' ] ]'
enddo
write(xyzzyaaaa51,*)'   Linear parameters:'
do xyzzyaaad51=1,no_spairs(xyzzyaafh1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafh1),which_spair(:,:,xyzzyaa&
&fh1)==xyzzyaaad51)
if(nele(xyzzyaaak51(1))<1.or.nele(xyzzyaaak51(2))<1.or.(xyzzyaaak51(1)&
&==xyzzyaaak51(2).and.nele(xyzzyaaak51(1))<2))cycle
write(xyzzyaaaa51,*)'     Channel '//trim(i2s(minval(xyzzyaaak51)))//'&
&-'//trim(i2s(maxval(xyzzyaaak51)))//':'
xyzzyaaae51=which_spair(xyzzyaaak51(1),xyzzyaaak51(2),levels_spairs)
write(char_80,*)xyzzyaaan1(xyzzyaaad51)
write(xyzzyaaaa51,*)'       c_1: [ ',trim(adjustl(char_80)),', determi&
&ned ]'
enddo
call wout('Successfully converted RPA U term.')
endif
if(xyzzyaaes1)call errstop('CONVERT_TO_GJASTROW','Conversion of Ucyl n&
&ot implemented.  Sorry.')
if(xyzzyaaet1)then
xyzzyaaac51=xyzzyaaac51+1
write(xyzzyaaaa51,*)' TERM '//trim(i2s(xyzzyaaac51))//':'
write(xyzzyaaaa51,*)'   Rank: [ 2, 0 ]'
write(xyzzyaaaa51,*)'   e-e cusp: T'
allocate(xyzzyaabj51(heg_nlayers))
xyzzyaabj51=0
xyzzyaabf51=0
do xyzzyaabd51=1,heg_nlayers
xyzzyaabg51=0
char_80=''
do xyzzyaaam51=1,nspin
if(nele(xyzzyaaam51)<1)cycle
if(heg_layer(xyzzyaaam51)==xyzzyaabd51)then
xyzzyaabg51=xyzzyaabg51+1
if(xyzzyaabj51(xyzzyaabd51)==0)then
xyzzyaabj51(xyzzyaabd51)=xyzzyaaam51
char_80=trim(i2s(xyzzyaaam51))
else
char_80=trim(char_80)//'='//trim(i2s(xyzzyaaam51))
endif
endif
enddo
if(xyzzyaabg51>1)then
if(xyzzyaabf51==0)write(xyzzyaaaa51,*)'   Rules:'
write(xyzzyaaaa51,*)'      '//trim(char_80)
xyzzyaabf51=xyzzyaabf51+1
endif
if(xyzzyaabj51(xyzzyaabd51)>0)write(xyzzyaaaa51,*)'      !'//trim(i2s(&
&xyzzyaabj51(xyzzyaabd51)))//'-'//trim(i2s(xyzzyaabj51(xyzzyaabd51)))
enddo
if(xyzzyaabf51==0)write(xyzzyaaaa51,*)'   Rules: [ ]'
write(xyzzyaaaa51,*)'   e-e basis: [ Type: none ]'
write(xyzzyaaaa51,*)'   e-e cutoff:'
write(xyzzyaaaa51,*)'     Type: quasicusp'
write(xyzzyaaaa51,*)'     Parameters:'
write(char_80,*)xyzzyaacq1
do xyzzyaabd51=1,heg_nlayers
if(xyzzyaabj51(xyzzyaabd51)==0)cycle
do xyzzyaabe51=xyzzyaabd51+1,heg_nlayers
if(xyzzyaabj51(xyzzyaabe51)==0)cycle
write(xyzzyaaaa51,*)'       Channel '//trim(i2s(xyzzyaabj51(xyzzyaabd5&
&1)))//'-'//trim(i2s(xyzzyaabj51(xyzzyaabe51)))//': [ L: [ ',trim(adju&
&stl(char_80)),', '//trim(xyzzyaabm51(xyzzyaadc1))//' ] ]'
enddo
enddo
write(xyzzyaaaa51,*)'   Linear parameters:'
do xyzzyaabd51=1,heg_nlayers
if(xyzzyaabj51(xyzzyaabd51)==0)cycle
do xyzzyaabe51=xyzzyaabd51+1,heg_nlayers
if(xyzzyaabj51(xyzzyaabe51)==0)cycle
write(xyzzyaaaa51,*)'     Channel '//trim(i2s(xyzzyaabj51(xyzzyaabd51)&
&))//'-'//trim(i2s(xyzzyaabj51(xyzzyaabe51)))//':'
write(char_80,*)ee_kato_gamma(xyzzyaabj51(xyzzyaabd51),xyzzyaabj51(xyz&
&zyaabe51),f_zero=f_zero)
write(xyzzyaaaa51,*)'       c_1: [ '//trim(adjustl(char_80))//', fixed&
& ]'
enddo
enddo
deallocate(xyzzyaabj51)
call wout('Successfully converted QCUSP term.')
endif
if(xyzzyaaeu1)then
xyzzyaaac51=xyzzyaaac51+1
write(xyzzyaaaa51,*)' TERM '//trim(i2s(xyzzyaaac51))//':'
write(xyzzyaaaa51,*)'   Rank: [ 3, 0 ]'
if(maxval(which_spair(:,:,xyzzyaafl1))==(nspin*(nspin+1))/2)then
write(xyzzyaaaa51,*)'   Rules: [ ]'
else
write(xyzzyaaaa51,*)'   Rules:'
do xyzzyaaad51=1,no_spairs(xyzzyaafl1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafl1),which_spair(:,:,xyzzyaa&
&fl1)==xyzzyaaad51)
xyzzyaaam51=minval(xyzzyaaak51)
xyzzyaaan51=maxval(xyzzyaaak51)
do xyzzyaaao51=xyzzyaaam51,nspin
xyzzyaaaq51=xyzzyaaao51
if(xyzzyaaao51==xyzzyaaam51)xyzzyaaaq51=xyzzyaaan51+1
do xyzzyaaap51=xyzzyaaaq51,nspin
if(nele(xyzzyaaam51)<1.or.nele(xyzzyaaan51)<1.or.xyzzyaaam51==xyzzyaaa&
&n51.and.nele(xyzzyaaam51)<2)then
xyzzyaaam51=xyzzyaaao51
xyzzyaaan51=xyzzyaaap51
cycle
endif
if(nele(xyzzyaaao51)<1.or.nele(xyzzyaaap51)<1.or.xyzzyaaao51==xyzzyaaa&
&p51.and.nele(xyzzyaaao51)<2)cycle
if(which_spair(xyzzyaaao51,xyzzyaaap51,xyzzyaafl1)==which_spair(xyzzya&
&aam51,xyzzyaaan51,xyzzyaafl1))write(xyzzyaaaa51,*)'     ',trim(i2s(xy&
&zzyaaam51)),'-',trim(i2s(xyzzyaaan51)),'=',trim(i2s(xyzzyaaao51)),'-'&
&,trim(i2s(xyzzyaaap51))
enddo
enddo
enddo
endif
write(xyzzyaaaa51,*)'   Indexing:'
write(xyzzyaaaa51,*)'     e-e dot product: T'
write(xyzzyaaaa51,*)'   e-e basis:'
write(xyzzyaaaa51,*)'     Type: natural polynomial vectorial'
write(xyzzyaaaa51,*)'     Order: '//trim(i2s(dimensionality))
write(xyzzyaaaa51,*)'     Constants:'
write(xyzzyaaaa51,*)'       k0: 1'
write(xyzzyaaaa51,*)'       Split: [ '//trim(i2s(xyzzyaacb1+1))//' ]'
write(xyzzyaaaa51,*)'     Parameters:'
do xyzzyaaad51=1,no_spairs(xyzzyaafl1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafl1),which_spair(:,:,xyzzyaa&
&fl1)==xyzzyaaad51)
if(nele(xyzzyaaak51(1))<1.or.nele(xyzzyaaak51(2))<1.or.(xyzzyaaak51(1)&
&==xyzzyaaak51(2).and.nele(xyzzyaaak51(1))<2))cycle
write(xyzzyaaaa51,*)'       Channel '//trim(i2s(minval(xyzzyaaak51)))/&
&/'-'//trim(i2s(maxval(xyzzyaaak51)))//':'
do xyzzyaaai51=1,xyzzyaacb1
if(xyzzyaaap1(0,xyzzyaaad51)/=0.d0)then
write(char_80,*)xyzzyaaap1(xyzzyaaai51,xyzzyaaad51)/xyzzyaaap1(0,xyzzy&
&aaad51)
else
write(char_80,*)0.d0
endif
write(xyzzyaaaa51,*)'         c_'//trim(i2s(xyzzyaaai51+1))//': [ '//t&
&rim(adjustl(char_80))//', '//trim(xyzzyaabm51(xyzzyaabb1(xyzzyaaai51,&
&xyzzyaaad51)))//' ]'
enddo
enddo
write(xyzzyaaaa51,*)'   e-e cutoff:'
write(xyzzyaaaa51,*)'     Type: alt polynomial'
write(xyzzyaaaa51,*)'     Constants: [ C: '//trim(i2s(xyzzyaagx1))//' &
&]'
write(xyzzyaaaa51,*)'     Parameters:'
write(char_80,*)xyzzyaacl1
do xyzzyaaad51=1,no_spairs(xyzzyaafl1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafl1),which_spair(:,:,xyzzyaa&
&fl1)==xyzzyaaad51)
if(nele(xyzzyaaak51(1))<1.or.nele(xyzzyaaak51(2))<1.or.(xyzzyaaak51(1)&
&==xyzzyaaak51(2).and.nele(xyzzyaaak51(1))<2))cycle
write(xyzzyaaaa51,*)'       Channel '//trim(i2s(minval(xyzzyaaak51)))/&
&/'-'//trim(i2s(maxval(xyzzyaaak51)))//':'
write(xyzzyaaaa51,*)'         L: [ ',trim(adjustl(char_80)),', '//trim&
&(xyzzyaabm51(xyzzyaada1))//' ]'
enddo
write(xyzzyaaaa51,*)'   Linear parameters:'
xyzzyaaau51=nspin**3-2*nspin*(nspin-1)
allocate(xyzzyaabh51(3,xyzzyaaau51),xyzzyaabi51(3,xyzzyaaau51),stat=xy&
&zzyaaav51)
call check_alloc(xyzzyaaav51,'CONVERT_T_GJASTROW','channel_signature')
xyzzyaabh51=0
xyzzyaabi51=0
xyzzyaaax51=0
do xyzzyaaas51=1,nspin
if(nele(xyzzyaaas51)<1)cycle
do xyzzyaaao51=xyzzyaaas51,nspin
if(nele(xyzzyaaao51)<1.or.(xyzzyaaas51==xyzzyaaao51.and.nele(xyzzyaaas&
&51)<2))cycle
do xyzzyaaat51=xyzzyaaao51,nspin
if(nele(xyzzyaaat51)<1.or.(xyzzyaaas51==xyzzyaaat51.and.nele(xyzzyaaas&
&51)<2).or.(xyzzyaaao51==xyzzyaaat51.and.nele(xyzzyaaao51)<2).or.(xyzz&
&yaaas51==xyzzyaaao51.and.xyzzyaaas51==xyzzyaaat51.and.nele(xyzzyaaas5&
&1)<3))cycle
xyzzyaaay51=(/which_spair(xyzzyaaas51,xyzzyaaao51,xyzzyaafl1),which_sp&
&air(xyzzyaaas51,xyzzyaaat51,xyzzyaafl1),which_spair(xyzzyaaao51,xyzzy&
&aaat51,xyzzyaafl1)/)
if(xyzzyaaay51(1)>xyzzyaaay51(2))call swap1(xyzzyaaay51(1),xyzzyaaay51&
&(2))
if(xyzzyaaay51(1)>xyzzyaaay51(3))call swap1(xyzzyaaay51(1),xyzzyaaay51&
&(3))
if(xyzzyaaay51(2)>xyzzyaaay51(3))call swap1(xyzzyaaay51(2),xyzzyaaay51&
&(3))
xyzzyaabl51=.false.
do xyzzyaaaw51=1,xyzzyaaax51
if(all(xyzzyaaay51==xyzzyaabh51(:,xyzzyaaaw51)))then
xyzzyaabl51=.true.
exit
endif
enddo
if(.not.xyzzyaabl51)then
xyzzyaaax51=xyzzyaaax51+1
xyzzyaabh51(:,xyzzyaaax51)=xyzzyaaay51
xyzzyaabi51(:,xyzzyaaax51)=(/xyzzyaaas51,xyzzyaaao51,xyzzyaaat51/)
endif
enddo
enddo
enddo
do xyzzyaaaw51=1,xyzzyaaax51
write(xyzzyaaaa51,*)'     Channel '//trim(i2s(xyzzyaabi51(1,xyzzyaaaw5&
&1)))//'-'//trim(i2s(xyzzyaabi51(2,xyzzyaaaw51)))//'-'//trim(i2s(xyzzy&
&aabi51(3,xyzzyaaaw51)))//':'
xyzzyaaaz51=xyzzyaabh51(1,xyzzyaaaw51)
xyzzyaaba51=xyzzyaabh51(2,xyzzyaaaw51)
xyzzyaabb51=xyzzyaabh51(3,xyzzyaaaw51)
write(char_80_1,*)2*xyzzyaaap1(0,xyzzyaaaz51)*xyzzyaaap1(0,xyzzyaaba51&
&)
write(char_80_2,*)-2*xyzzyaaap1(0,xyzzyaaaz51)*xyzzyaaap1(0,xyzzyaabb5&
&1)
write(char_80_3,*)2*xyzzyaaap1(0,xyzzyaaba51)*xyzzyaaap1(0,xyzzyaabb51&
&)
do xyzzyaabc51=1,dimensionality
write(xyzzyaaaa51,*)'       c_'//trim(i2s(xyzzyaabc51))//','//trim(i2s&
&(xyzzyaabc51))//',0: [ '//trim(adjustl(char_80_1))//', optimizable ]'
write(xyzzyaaaa51,*)'       c_'//trim(i2s(xyzzyaabc51))//',0,'//trim(i&
&2s(xyzzyaabc51))//': [ '//trim(adjustl(char_80_2))//', optimizable ]'
write(xyzzyaaaa51,*)'       c_0,'//trim(i2s(xyzzyaabc51))//','//trim(i&
&2s(xyzzyaabc51))//': [ '//trim(adjustl(char_80_3))//', optimizable ]'
enddo
enddo
deallocate(xyzzyaabh51,xyzzyaabi51)
call wout('Successfully converted W term.')
endif
if(xyzzyaaev1)then
xyzzyaaac51=xyzzyaaac51+1
write(xyzzyaaaa51,*)' TERM '//trim(i2s(xyzzyaaac51))//':'
write(xyzzyaaaa51,*)'   Rank: [ 3, 0 ]'
write(xyzzyaaaa51,*)'   # CASINO can''t convert triplets into GJASTROW&
& rules!'
write(xyzzyaaaa51,*)'   # Write the appropriate rules here:'
write(xyzzyaaaa51,*)'   Rules: [ xxx ]'
write(xyzzyaaaa51,*)'   e-e basis: [ Type: natural power, Order: '//tr&
&im(i2s(xyzzyaacc1+1))//' ]'
write(xyzzyaaaa51,*)'   e-e cutoff:'
write(xyzzyaaaa51,*)'     Type: polynomial'
write(xyzzyaaaa51,*)'     Constants: [ C: '//trim(i2s(xyzzyaagx1))//' &
&]'
write(xyzzyaaaa51,*)'     Parameters:'
write(xyzzyaaaa51,*)'       # CASINO can''t convert triplets into GJAS&
&TROW e-e cutoff channels!'
write(xyzzyaaaa51,*)'       # Copy this line for each e-e channel gene&
&rated by the above rules:'
write(char_80,*)xyzzyaacm1
write(xyzzyaaaa51,*)'       Channel xxx: [ L: [ ',trim(adjustl(char_80&
&)),', '//trim(xyzzyaabm51(xyzzyaadb1))//' ] ]'
write(xyzzyaaaa51,*)'   Linear parameters:'
do xyzzyaaad51=1,no_striplets(xyzzyaafm1)
xyzzyaaal51=minloc(which_striplet(:,:,:,xyzzyaafm1),which_striplet(:,:&
&,:,xyzzyaafm1)==xyzzyaaad51)
if(xyzzyaaal51(1)>xyzzyaaal51(2))xyzzyaaal51=(/xyzzyaaal51(2),xyzzyaaa&
&l51(1),xyzzyaaal51(3)/)
if(xyzzyaaal51(1)>xyzzyaaal51(3))xyzzyaaal51=(/xyzzyaaal51(3),xyzzyaaa&
&l51(1),xyzzyaaal51(2)/)
if(xyzzyaaal51(2)>xyzzyaaal51(3))xyzzyaaal51=(/xyzzyaaal51(1),xyzzyaaa&
&l51(3),xyzzyaaal51(2)/)
if(nele(xyzzyaaal51(1))<1.or.nele(xyzzyaaal51(2))<1.or.nele(xyzzyaaal5&
&1(3))<1.or.(xyzzyaaal51(1)==xyzzyaaal51(2).and.nele(xyzzyaaal51(1))<2&
&).or.(xyzzyaaal51(1)==xyzzyaaal51(3).and.nele(xyzzyaaal51(1))<2).or.(&
&xyzzyaaal51(2)==xyzzyaaal51(3).and.nele(xyzzyaaal51(2))<2).or.(xyzzya&
&aal51(1)==xyzzyaaal51(2).and.xyzzyaaal51(1)==xyzzyaaal51(3).and.nele(&
&xyzzyaaal51(1))<3))cycle
write(xyzzyaaaa51,*)'     Channel '//trim(i2s(xyzzyaaal51(1)))//'-'//t&
&rim(i2s(xyzzyaaal51(2)))//'-'//trim(i2s(xyzzyaaal51(3)))//':'
do xyzzyaaaj51=0,xyzzyaacc1
do xyzzyaaai51=0,xyzzyaacc1
do xyzzyaaah51=0,xyzzyaacc1
write(char_80,*)6*xyzzyaaar1(xyzzyaaah51,xyzzyaaai51,xyzzyaaaj51,xyzzy&
&aaad51)
write(xyzzyaaaa51,*)'       c_'//trim(i2s(xyzzyaaaj51+1))//','//trim(i&
&2s(xyzzyaaai51+1)),','//trim(i2s(xyzzyaaah51+1))//': [ ',trim(adjustl&
&(char_80)),', '//trim(xyzzyaabm51(xyzzyaabc1(xyzzyaaah51,xyzzyaaai51,&
&xyzzyaaaj51,xyzzyaaad51)))//' ]'
enddo
enddo
enddo
enddo
call wout('Successfully converted H term (manual edits needed).')
endif
if(xyzzyaaew1)then
xyzzyaaac51=xyzzyaaac51+1
write(xyzzyaaaa51,*)' TERM '//trim(i2s(xyzzyaaac51))//':'
write(xyzzyaaaa51,*)'   Rank: [ 1, 1 ]'
xyzzyaabk51=.true.
do xyzzyaaag51=1,xyzzyaabi1
if(no_ssingles(xyzzyaafn1(xyzzyaaag51))>1)xyzzyaabk51=.false.
enddo
if(all(xyzzyaadi1==1))then
if(.not.xyzzyaabk51.or.nele(1)<1.or.nele(2)<1)then
write(xyzzyaaaa51,*)'   Rules: [ ]'
else
write(xyzzyaaaa51,*)'   Rules: [ 1=2 ]'
endif
else
write(xyzzyaaaa51,*)'   Rules:'
do xyzzyaaag51=1,xyzzyaabi1
if(xyzzyaadi1(xyzzyaaag51)>1)then
do xyzzyaaah51=2,xyzzyaadi1(xyzzyaaag51)
write(xyzzyaaaa51,*)'     n'//trim(i2s(xyzzyaabq1(1,xyzzyaaag51)))//'=&
&n'//trim(i2s(xyzzyaabq1(xyzzyaaah51,xyzzyaaag51)))
enddo
endif
enddo
if(xyzzyaabk51.and.nele(1)>0.and.nele(2)>0)write(xyzzyaaaa51,*)'     1&
&=2'
endif
write(xyzzyaaaa51,*)'   e-n basis: [ Type: natural power, Order: '//tr&
&im(i2s(xyzzyaace1(1)+1))//' ]'
write(xyzzyaaaa51,*)'   e-n cutoff:'
write(xyzzyaaaa51,*)'     Type: alt polynomial'
write(xyzzyaaaa51,*)'     Constants: [ C: '//trim(i2s(xyzzyaagx1))//' &
&]'
write(xyzzyaaaa51,*)'     Parameters:'
do xyzzyaaag51=1,xyzzyaabi1
write(char_80,*)xyzzyaacu1(xyzzyaaag51)
if(nele(1)>0)write(xyzzyaaaa51,*)'       Channel 1-n',trim(i2s(xyzzyaa&
&bq1(1,xyzzyaaag51))),': [ L: [ ',trim(adjustl(char_80)),', '//trim(xy&
&zzyaabm51(xyzzyaadf1(xyzzyaaag51)))//' ] ]'
if((.not.xyzzyaabk51.and.nele(2)>0).or.nele(1)<1)write(xyzzyaaaa51,*)'&
&       Channel 2-n',trim(i2s(xyzzyaabq1(1,xyzzyaaag51))),': [ L: [ ',&
&trim(adjustl(char_80)),', '//trim(xyzzyaabm51(xyzzyaadf1(xyzzyaaag51)&
&))//' ] ]'
enddo
write(xyzzyaaaa51,*)'   Linear parameters:'
do xyzzyaaag51=1,xyzzyaabi1
do xyzzyaaad51=1,no_ssingles(xyzzyaafn1(xyzzyaaag51))
if(nele(xyzzyaaad51)<1)cycle
write(xyzzyaaaa51,*)'     Channel '//trim(i2s(xyzzyaaad51))//'-n'//tri&
&m(i2s(xyzzyaabq1(1,xyzzyaaag51)))//':'
if(xyzzyaadn1(xyzzyaaag51)==1)write(xyzzyaaaa51,*)'       e-n cusp: T'
do xyzzyaaaf51=0,xyzzyaace1(xyzzyaaag51)
write(char_80,*)xyzzyaaal1(xyzzyaaaf51,xyzzyaaad51,xyzzyaaag51)
write(xyzzyaaaa51,*)'       c_'//trim(i2s(xyzzyaaaf51+1))//': [ ',trim&
&(adjustl(char_80)),', '//trim(xyzzyaabm51(xyzzyaaav1(xyzzyaaaf51,xyzz&
&yaaad51,xyzzyaaag51)))//' ]'
enddo
enddo
enddo
call wout('Successfully converted Chi term.')
endif
if(xyzzyaaex1)then
xyzzyaaac51=xyzzyaaac51+1
write(xyzzyaaaa51,*)' TERM '//trim(i2s(xyzzyaaac51))//':'
write(xyzzyaaaa51,*)'   Rank: [ 2, 1 ]'
write(xyzzyaaaa51,*)'   Rules:'
do xyzzyaaag51=1,xyzzyaabj1
if(nele(1)>0.and.nele(2)>0)write(xyzzyaaaa51,*)'     1-n'//trim(i2s(xy&
&zzyaabs1(1,xyzzyaaag51)))//'='//'2-n'//trim(i2s(xyzzyaabs1(1,xyzzyaaa&
&g51)))
if(xyzzyaadj1(xyzzyaaag51)>1)then
do xyzzyaaah51=2,xyzzyaadj1(xyzzyaaag51)
write(xyzzyaaaa51,*)'     n'//trim(i2s(xyzzyaabs1(1,xyzzyaaag51)))//'=&
&n'//trim(i2s(xyzzyaabs1(xyzzyaaah51,xyzzyaaag51)))
enddo
endif
enddo
xyzzyaaar51=maxval(xyzzyaafo1)
do xyzzyaaad51=1,no_spairs(xyzzyaaar51)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaaar51),which_spair(:,:,xyzzya&
&aar51)==xyzzyaaad51)
xyzzyaaam51=minval(xyzzyaaak51)
xyzzyaaan51=maxval(xyzzyaaak51)
do xyzzyaaao51=xyzzyaaam51,nspin
xyzzyaaaq51=xyzzyaaao51
if(xyzzyaaao51==xyzzyaaam51)xyzzyaaaq51=xyzzyaaan51+1
do xyzzyaaap51=xyzzyaaaq51,nspin
if(nele(xyzzyaaam51)<1.or.nele(xyzzyaaan51)<1.or.xyzzyaaam51==xyzzyaaa&
&n51.and.nele(xyzzyaaam51)<2)then
xyzzyaaam51=xyzzyaaao51
xyzzyaaan51=xyzzyaaap51
cycle
endif
if(nele(xyzzyaaao51)<1.or.nele(xyzzyaaap51)<1.or.xyzzyaaao51==xyzzyaaa&
&p51.and.nele(xyzzyaaao51)<2)cycle
if(which_spair(xyzzyaaao51,xyzzyaaap51,xyzzyaaar51)==which_spair(xyzzy&
&aaam51,xyzzyaaan51,xyzzyaaar51))write(xyzzyaaaa51,*)'     ',trim(i2s(&
&xyzzyaaam51)),'-',trim(i2s(xyzzyaaan51)),'=',trim(i2s(xyzzyaaao51)),'&
&-',trim(i2s(xyzzyaaap51))
enddo
enddo
enddo
write(xyzzyaaaa51,*)'   e-e basis: [ Type: natural power, Order: '//tr&
&im(i2s(xyzzyaacf1(1)+1))//' ]'
write(xyzzyaaaa51,*)'   e-n basis: [ Type: natural power, Order: '//tr&
&im(i2s(xyzzyaacg1(1)+1))//' ]'
write(xyzzyaaaa51,*)'   e-n cutoff:'
write(xyzzyaaaa51,*)'     Type: alt polynomial'
write(xyzzyaaaa51,*)'     Constants: [ C: '//trim(i2s(xyzzyaagx1))//' &
&]'
write(xyzzyaaaa51,*)'     Parameters:'
do xyzzyaaag51=1,xyzzyaabj1
write(char_80,*)xyzzyaacv1(xyzzyaaag51)
write(xyzzyaaaa51,*)'       Channel 1-n',trim(i2s(xyzzyaabs1(1,xyzzyaa&
&ag51))),': [ L: [ ',trim(adjustl(char_80)),', '//trim(xyzzyaabm51(xyz&
&zyaadg1(xyzzyaaag51)))//' ] ]'
enddo
write(xyzzyaaaa51,*)'   Linear parameters:'
do xyzzyaaag51=1,xyzzyaabj1
do xyzzyaaad51=1,no_spairs(xyzzyaafo1(xyzzyaaag51))
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafo1(xyzzyaaag51)),which_spai&
&r(:,:,xyzzyaafo1(xyzzyaaag51))==xyzzyaaad51)
if(nele(xyzzyaaak51(1))<1.or.nele(xyzzyaaak51(2))<1.or.(xyzzyaaak51(1)&
&==xyzzyaaak51(2).and.nele(xyzzyaaak51(1))<2))cycle
write(xyzzyaaaa51,*)'     Channel '//trim(i2s(minval(xyzzyaaak51)))//'&
&-'//trim(i2s(maxval(xyzzyaaak51)))//'-n'//trim(i2s(xyzzyaabs1(1,xyzzy&
&aaag51)))//':'
xyzzyaaae51=which_spair(xyzzyaaak51(1),xyzzyaaak51(2),levels_spairs)
xyzzyaaaf51=0
do xyzzyaaaj51=0,xyzzyaacf1(xyzzyaaag51)
do xyzzyaaai51=0,xyzzyaacg1(xyzzyaaag51)
do xyzzyaaah51=0,xyzzyaacg1(xyzzyaaag51)
xyzzyaaaf51=xyzzyaaaf51+1
write(char_80,*)xyzzyaaam1(xyzzyaaaf51,xyzzyaaad51,xyzzyaaag51)
write(xyzzyaaaa51,*)'       c_'//trim(i2s(xyzzyaaaj51+1))//','//trim(i&
&2s(xyzzyaaai51+1)),','//trim(i2s(xyzzyaaah51+1))//': [ ',trim(adjustl&
&(char_80)),', '//trim(xyzzyaabm51(xyzzyaaaw1(xyzzyaaah51,xyzzyaaai51,&
&xyzzyaaaj51,xyzzyaaad51,xyzzyaaag51)))//' ]'
enddo
enddo
enddo
enddo
enddo
call wout('Successfully converted f term.')
endif
if(xyzzyaaey1)then
xyzzyaaac51=xyzzyaaac51+1
write(xyzzyaaaa51,*)' TERM '//trim(i2s(xyzzyaaac51))//':'
write(xyzzyaaaa51,*)'   Rank: [ 2, 0 ]'
if(maxval(which_spair(:,:,xyzzyaafj1))==(nspin*(nspin+1))/2)then
write(xyzzyaaaa51,*)'   Rules: [ ]'
else
write(xyzzyaaaa51,*)'   Rules:'
do xyzzyaaad51=1,no_spairs(xyzzyaafj1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafj1),which_spair(:,:,xyzzyaa&
&fj1)==xyzzyaaad51)
xyzzyaaam51=minval(xyzzyaaak51)
xyzzyaaan51=maxval(xyzzyaaak51)
do xyzzyaaao51=xyzzyaaam51,nspin
xyzzyaaaq51=xyzzyaaao51
if(xyzzyaaao51==xyzzyaaam51)xyzzyaaaq51=xyzzyaaan51+1
do xyzzyaaap51=xyzzyaaaq51,nspin
if(nele(xyzzyaaam51)<1.or.nele(xyzzyaaan51)<1.or.xyzzyaaam51==xyzzyaaa&
&n51.and.nele(xyzzyaaam51)<2)then
xyzzyaaam51=xyzzyaaao51
xyzzyaaan51=xyzzyaaap51
cycle
endif
if(nele(xyzzyaaao51)<1.or.nele(xyzzyaaap51)<1.or.xyzzyaaao51==xyzzyaaa&
&p51.and.nele(xyzzyaaao51)<2)cycle
if(which_spair(xyzzyaaao51,xyzzyaaap51,xyzzyaafj1)==which_spair(xyzzya&
&aam51,xyzzyaaan51,xyzzyaafj1))write(xyzzyaaaa51,*)'     ',trim(i2s(xy&
&zzyaaam51)),'-',trim(i2s(xyzzyaaan51)),'=',trim(i2s(xyzzyaaao51)),'-'&
&,trim(i2s(xyzzyaaap51))
enddo
enddo
enddo
endif
write(xyzzyaaaa51,*)'   e-e basis:'
write(xyzzyaaaa51,*)'     Type: cosine'
write(xyzzyaaaa51,*)'     Order: ',trim(i2s(xyzzyaagq1))
write(xyzzyaaaa51,*)'     Constants:'
do xyzzyaaah51=1,xyzzyaagq1
write(xyzzyaaaa51,*)'       Star '//trim(i2s(xyzzyaaah51))//':'
xyzzyaaaj51=0
do xyzzyaaai51=1,xyzzyaafy1
if(xyzzyaagb1(xyzzyaaai51)/=xyzzyaaah51)cycle
xyzzyaaaj51=xyzzyaaaj51+1
char_80=trim(i2s(xyzzyaaga1(1,xyzzyaaai51)))
if(periodicity>1)char_80=trim(char_80)//', '//trim(i2s(xyzzyaaga1(2,xy&
&zzyaaai51)))
if(periodicity>2)char_80=trim(char_80)//', '//trim(i2s(xyzzyaaga1(3,xy&
&zzyaaai51)))
write(xyzzyaaaa51,*)'         G_',trim(i2s(xyzzyaaaj51)),': [ '//trim(&
&char_80)//' ]'
enddo
enddo
write(xyzzyaaaa51,*)'   Linear parameters:'
do xyzzyaaad51=1,no_spairs(xyzzyaafj1)
xyzzyaaak51=minloc(which_spair(:,:,xyzzyaafj1),which_spair(:,:,xyzzyaa&
&fj1)==xyzzyaaad51)
if(nele(xyzzyaaak51(1))<1.or.nele(xyzzyaaak51(2))<1.or.(xyzzyaaak51(1)&
&==xyzzyaaak51(2).and.nele(xyzzyaaak51(1))<2))cycle
write(xyzzyaaaa51,*)'     Channel '//trim(i2s(minval(xyzzyaaak51)))//'&
&-'//trim(i2s(maxval(xyzzyaaak51)))//':'
xyzzyaaae51=which_spair(xyzzyaaak51(1),xyzzyaaak51(2),levels_spairs)
do xyzzyaaaf51=1,xyzzyaagq1
write(char_80,*)xyzzyaaaf1(xyzzyaaaf51,xyzzyaaad51)
write(xyzzyaaaa51,*)'       c_'//trim(i2s(xyzzyaaaf51))//': [ ',trim(a&
&djustl(char_80)),', '//trim(xyzzyaabm51(xyzzyaaax1(xyzzyaaaf51,xyzzya&
&aad51)))//' ]'
enddo
enddo
call wout('Successfully converted P term.')
endif
if(xyzzyaaez1)call wout('Cannot convert Q term.')
if(xyzzyaaff1)call wout('Cannot convert D term.')
if(xyzzyaaes1)call wout('Cannot convert UCYL term.')
if(xyzzyaafa1.or.xyzzyaafb1.or.xyzzyaafc1.or.xyzzyaafg1)call wout('Can&
&not convert BIEX/EX2D terms. Use the EXMOL facility instead.')
close(xyzzyaaaa51)
open_unit(xyzzyaaaa51)=.false.
call wout('Conversion done.')
call wout()
end subroutine xyzzyaaoj1
subroutine write_pjastrow(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
integer xyzzyaaaa52,xyzzyaaab52,xyzzyaaac52,xyzzyaaad52,xyzzyaaae52,xy&
&zzyaaaf52,xyzzyaaag52,xyzzyaaah52,xyzzyaaai52,xyzzyaaaj52
logical xyzzyaaak52
character(64) blurb
if(am_master)then
inquire(file=trim(correlation_name),exist=xyzzyaaak52)
if(xyzzyaaak52)then
open(unit=xyzzyaagw1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaaa52)
else
open(unit=xyzzyaagw1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaaa52)
endif
if(xyzzyaaaa52/=0)call errstop('WRITE_JASTROW','Problem opening '//tri&
&m(correlation_name)//'.')
write(xyzzyaagw1,*)'START JASTROW'
write(xyzzyaagw1,*)'Title'
write(xyzzyaagw1,*)trim(adjustl(title))
write(xyzzyaagw1,*)'Truncation order C'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaagx1))
if(xyzzyaaep1)then
if(xyzzyaaeq1)then
write(xyzzyaagw1,*)'START U TERM'
else
write(xyzzyaagw1,*)'START U_RPA TERM'
endif
write(xyzzyaagw1,*)'Number of sets'
write(xyzzyaagw1,*)'  1'
write(xyzzyaagw1,*)'START SET 1'
write(xyzzyaagw1,*)'Spherical harmonic l,m'
write(xyzzyaagw1,*)'  0 0'
if(xyzzyaaeq1)then
write(xyzzyaagw1,*)'Expansion order N_u'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaabx1))
endif
write(xyzzyaagw1,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud&
&)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafh1))
write(xyzzyaagw1,*)'Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)'
write(xyzzyaagw1,*)xyzzyaach1,xyzzyaacx1
write(xyzzyaagw1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
if(xyzzyaaeq1)then
do xyzzyaaac52=1,xyzzyaafp1
do xyzzyaaad52=0,xyzzyaabx1
if(xyzzyaaad52/=1)then
blurb='      ! alpha_'//trim(i2s(xyzzyaaad52))//','//trim(i2s(xyzzyaaa&
&c52))
write(xyzzyaagw1,*)xyzzyaaab1(xyzzyaaad52,xyzzyaaac52),xyzzyaaat1(xyzz&
&yaaad52,xyzzyaaac52),trim(blurb)
endif
enddo
enddo
else
do xyzzyaaac52=1,xyzzyaafp1
blurb='      ! F_rpa'//','//trim(i2s(xyzzyaaac52))
write(xyzzyaagw1,*)xyzzyaaao1(xyzzyaaac52),xyzzyaaba1(xyzzyaaac52),tri&
&m(blurb)
enddo
endif
write(xyzzyaagw1,*)'END SET 1'
if(xyzzyaaeq1)then
write(xyzzyaagw1,*)'END U TERM'
else
write(xyzzyaagw1,*)'END U_RPA TERM'
endif
endif
if(xyzzyaaes1)then
write(xyzzyaagw1,*)'START UCYL TERM'
write(xyzzyaagw1,*)'Expansion order N_ucylrho'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaaby1))
write(xyzzyaagw1,*)'Expansion order N_ucylz'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaabz1))
write(xyzzyaagw1,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud&
&)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafi1))
write(xyzzyaagw1,*)'Axis polar angle theta   ;  Optimizable (0=NO; 1=Y&
&ES)'
write(xyzzyaagw1,*)xyzzyaacs1,xyzzyaadd1
write(xyzzyaagw1,*)'Axis azimuthal angle phi ;  Optimizable (0=NO; 1=Y&
&ES)'
write(xyzzyaagw1,*)xyzzyaact1,xyzzyaade1
write(xyzzyaagw1,*)'Radial cutoff (a.u.)     ;  Optimizable (0=NO; 1=Y&
&ES)'
write(xyzzyaagw1,*)xyzzyaaci1,xyzzyaacy1
write(xyzzyaagw1,*)'Axial cutoff (a.u.)      ;  Optimizable (0=NO; 1=Y&
&ES)'
write(xyzzyaagw1,*)xyzzyaack1,xyzzyaacz1
write(xyzzyaagw1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaac52=1,xyzzyaafq1
do xyzzyaaae52=0,xyzzyaabz1
if(xyzzyaaae52==1)cycle
do xyzzyaaad52=0,xyzzyaaby1
if(xyzzyaaad52==1)cycle
blurb='      ! epsil_'//trim(i2s(xyzzyaaad52))//','//trim(i2s(xyzzyaaa&
&e52))//','//trim(i2s(xyzzyaaac52))
write(xyzzyaagw1,*)xyzzyaaac1(xyzzyaaad52,xyzzyaaae52,xyzzyaaac52),xyz&
&zyaaau1(xyzzyaaad52,xyzzyaaae52,xyzzyaaac52),trim(blurb)
enddo
enddo
enddo
write(xyzzyaagw1,*)'END UCYL TERM'
endif
if(xyzzyaaet1)then
write(xyzzyaagw1,*)'START QCUSP TERM'
write(xyzzyaagw1,*)'Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)'
write(xyzzyaagw1,*)xyzzyaacq1,xyzzyaadc1
write(xyzzyaagw1,*)'END QCUSP TERM'
endif
if(xyzzyaaeu1)then
write(xyzzyaagw1,*)'START W TERM'
write(xyzzyaagw1,*)'Number of sets'
write(xyzzyaagw1,*)'  1'
write(xyzzyaagw1,*)'START SET 1'
write(xyzzyaagw1,*)'Spherical harmonic l,m'
write(xyzzyaagw1,*)'  0 0'
write(xyzzyaagw1,*)'Expansion order N_w'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaacb1))
write(xyzzyaagw1,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud&
&)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafl1))
write(xyzzyaagw1,*)'Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)'
write(xyzzyaagw1,*)xyzzyaacl1,xyzzyaada1
write(xyzzyaagw1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaac52=1,xyzzyaaft1
do xyzzyaaad52=0,xyzzyaacb1
if(xyzzyaabb1(xyzzyaaad52,xyzzyaaac52)<0)cycle
blurb='      ! w_'//trim(i2s(xyzzyaaad52))//','//trim(i2s(xyzzyaaac52)&
&)
write(xyzzyaagw1,*)xyzzyaaap1(xyzzyaaad52,xyzzyaaac52),xyzzyaabb1(xyzz&
&yaaad52,xyzzyaaac52),trim(blurb)
enddo
enddo
write(xyzzyaagw1,*)'END SET 1'
write(xyzzyaagw1,*)'END W TERM'
endif
if(xyzzyaaev1)then
write(xyzzyaagw1,*)'START H TERM'
write(xyzzyaagw1,*)'Number of sets'
write(xyzzyaagw1,*)'  1'
write(xyzzyaagw1,*)'START SET 1'
write(xyzzyaagw1,*)'Spherical harmonic l,m'
write(xyzzyaagw1,*)'  0 0'
write(xyzzyaagw1,*)'Expansion order N_h'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaacc1))
write(xyzzyaagw1,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud&
&)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafm1))
write(xyzzyaagw1,*)'Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)'
write(xyzzyaagw1,*)xyzzyaacm1,xyzzyaadb1
write(xyzzyaagw1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaac52=1,xyzzyaafu1
do xyzzyaaaf52=0,xyzzyaacc1
do xyzzyaaae52=0,xyzzyaacc1
do xyzzyaaad52=0,xyzzyaacc1
if(xyzzyaabc1(xyzzyaaad52,xyzzyaaae52,xyzzyaaaf52,xyzzyaaac52)<0)cycle
blurb='      ! h_'//trim(i2s(xyzzyaaad52))//','//trim(i2s(xyzzyaaae52)&
&)//','//trim(i2s(xyzzyaaaf52))//','//trim(i2s(xyzzyaaac52))
write(xyzzyaagw1,*)xyzzyaaar1(xyzzyaaad52,xyzzyaaae52,xyzzyaaaf52,xyzz&
&yaaac52),xyzzyaabc1(xyzzyaaad52,xyzzyaaae52,xyzzyaaaf52,xyzzyaaac52),&
&trim(blurb)
enddo
enddo
enddo
enddo
write(xyzzyaagw1,*)'END SET 1'
write(xyzzyaagw1,*)'END H TERM'
endif
if(xyzzyaaew1)then
write(xyzzyaagw1,*)'START CHI TERM'
write(xyzzyaagw1,*)'Number of sets ; labelling (1->atom in s. cell; 2-&
&>atom in p. cell; 3->species)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaabi1))//' ' //trim(i2s(xyzzya&
&abf1))
do xyzzyaaab52=1,xyzzyaabi1
write(xyzzyaagw1,*)'START SET ',trim(i2s(xyzzyaaab52))
write(xyzzyaagw1,*)'Spherical harmonic l,m'
write(xyzzyaagw1,*)'  0 0'
if(xyzzyaabf1==3)then
write(xyzzyaagw1,*)'Number of species in set'
else
write(xyzzyaagw1,*)'Number of atoms in set'
endif
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaadl1(xyzzyaaab52)))
if(xyzzyaadl1(xyzzyaaab52)==1)then
if(xyzzyaabf1==3)then
write(xyzzyaagw1,*)'Label of the species in this set'
else
write(xyzzyaagw1,*)'Label of the atom in this set'
endif
else
if(xyzzyaabf1==3)then
write(xyzzyaagw1,*)'Labels of the species in this set'
else
write(xyzzyaagw1,*)'Labels of the atoms in this set'
endif
endif
call write_list_int(xyzzyaadl1(xyzzyaaab52),xyzzyaabr1(1:xyzzyaadl1(xy&
&zzyaaab52),xyzzyaaab52),10,4,1,xyzzyaagw1)
write(xyzzyaagw1,*)'Impose electron-nucleus cusp (0=NO; 1=YES)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaadn1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Expansion order N_chi'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaace1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Spin dep (0->u=d; 1->u/=d)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafn1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)'
write(xyzzyaagw1,*)xyzzyaacu1(xyzzyaaab52),xyzzyaadf1(xyzzyaaab52)
write(xyzzyaagw1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaac52=1,xyzzyaafv1(xyzzyaaab52)
do xyzzyaaae52=0,xyzzyaace1(xyzzyaaab52)
if(xyzzyaaae52/=1)then
blurb='      ! beta_'//trim(i2s(xyzzyaaae52))//','//trim(i2s(xyzzyaaac&
&52))//','//trim(i2s(xyzzyaaab52))
write(xyzzyaagw1,*)xyzzyaaad1(xyzzyaaae52,xyzzyaaac52,xyzzyaaab52),xyz&
&zyaaav1(xyzzyaaae52,xyzzyaaac52,xyzzyaaab52),trim(blurb)
endif
enddo
enddo
write(xyzzyaagw1,*)'END SET ',trim(i2s(xyzzyaaab52))
enddo
write(xyzzyaagw1,*)'END CHI TERM'
endif
if(xyzzyaaex1)then
write(xyzzyaagw1,*)'START F TERM'
write(xyzzyaagw1,*)'Number of sets ; labelling (1->atom in s. cell; 2-&
&>atom in p. cell; 3->species)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaabj1))//' '//trim(i2s(xyzzyaa&
&bg1))
do xyzzyaaab52=1,xyzzyaabj1
write(xyzzyaagw1,*)'START SET ',trim(i2s(xyzzyaaab52))
if(xyzzyaabg1==3)then
write(xyzzyaagw1,*)'Number of species in set'
else
write(xyzzyaagw1,*)'Number of atoms in set'
endif
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaadm1(xyzzyaaab52)))
if(xyzzyaadm1(xyzzyaaab52)==1)then
if(xyzzyaabg1==3)then
write(xyzzyaagw1,*)'Label of the species in this set'
else
write(xyzzyaagw1,*)'Label of the atom in this set'
endif
else
if(xyzzyaabg1==3)then
write(xyzzyaagw1,*)'Labels of the species in this set'
else
write(xyzzyaagw1,*)'Labels of the atoms in this set'
endif
endif
call write_list_int(xyzzyaadm1(xyzzyaaab52),xyzzyaabt1(1:xyzzyaadm1(xy&
&zzyaaab52),xyzzyaaab52),10,4,1,xyzzyaagw1)
write(xyzzyaagw1,*)'Prevent duplication of u term (0=NO; 1=YES)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaagy1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Prevent duplication of chi term (0=NO; 1=YES)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaagz1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Electron-nucleus expansion order N_f_eN'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaacg1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Electron-electron expansion order N_f_ee'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaacf1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud&
&)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafo1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)'
write(xyzzyaagw1,*)xyzzyaacv1(xyzzyaaab52),xyzzyaadg1(xyzzyaaab52)
write(xyzzyaagw1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaac52=1,xyzzyaafw1(xyzzyaaab52)
do xyzzyaaaf52=0,xyzzyaacf1(xyzzyaaab52)
do xyzzyaaae52=0,xyzzyaacg1(xyzzyaaab52)
do xyzzyaaad52=xyzzyaaae52,xyzzyaacg1(xyzzyaaab52)
if(xyzzyaaaw1(xyzzyaaad52,xyzzyaaae52,xyzzyaaaf52,xyzzyaaac52,xyzzyaaa&
&b52)>=0)then
blurb='      ! gamma_'//trim(i2s(xyzzyaaad52))//','//trim(i2s(xyzzyaaa&
&e52))//','//trim(i2s(xyzzyaaaf52))//','//trim(i2s(xyzzyaaac52))//','/&
&/trim(i2s(xyzzyaaab52))
write(xyzzyaagw1,*)xyzzyaaae1(xyzzyaaad52,xyzzyaaae52,xyzzyaaaf52,xyzz&
&yaaac52,xyzzyaaab52),xyzzyaaaw1(xyzzyaaad52,xyzzyaaae52,xyzzyaaaf52,x&
&yzzyaaac52,xyzzyaaab52),trim(blurb)
endif
enddo
enddo
enddo
enddo
write(xyzzyaagw1,*)'END SET ',trim(i2s(xyzzyaaab52))
enddo
write(xyzzyaagw1,*)'END F TERM'
endif
if(xyzzyaaey1)then
write(xyzzyaagw1,*)'START P TERM'
write(xyzzyaagw1,*)'Spin dep (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud&
&)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafj1))
write(xyzzyaagw1,*)'Number of simulation-cell G-vectors (NB, cannot ha&
&ve both G & -G)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafy1))
write(xyzzyaagw1,*)'G-vector (in terms of rec latt vects) ; label'
do xyzzyaaag52=1,xyzzyaafy1
write(xyzzyaagw1,'(3x,3(i5,1x),4x,i5)')xyzzyaaga1(1:3,xyzzyaaag52),xyz&
&zyaagb1(xyzzyaaag52)
enddo
write(xyzzyaagw1,*)'Parameter value ; Optimizable (0=NO; 1=YES)'
do xyzzyaaac52=1,xyzzyaafr1
do xyzzyaaag52=1,xyzzyaagq1
blurb='      ! a_'//trim(i2s(xyzzyaaag52))//','//trim(i2s(xyzzyaaac52)&
&)
write(xyzzyaagw1,*)xyzzyaaaf1(xyzzyaaag52,xyzzyaaac52),xyzzyaaax1(xyzz&
&yaaag52,xyzzyaaac52),trim(blurb)
enddo
enddo
write(xyzzyaagw1,*)'END P TERM'
endif
if(xyzzyaaez1)then
write(xyzzyaagw1,*)'START Q TERM'
write(xyzzyaagw1,*)'Spin dep (0->u=d; 1->u/=d)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafk1))
write(xyzzyaagw1,*)'Number of primitive-cell G-vectors (NB, cannot hav&
&e both G & -G)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaafz1))
write(xyzzyaagw1,*)'G-vector (in terms of rec latt vects) ; label'
do xyzzyaaag52=1,xyzzyaafz1
write(xyzzyaagw1,'(3x,3(i5,1x),4x,i5)')xyzzyaagc1(1:3,xyzzyaaag52),xyz&
&zyaagd1(xyzzyaaag52)
enddo
write(xyzzyaagw1,*)'Parameter value ; Optimizable (0=NO; 1=YES)'
do xyzzyaaac52=1,xyzzyaafs1
do xyzzyaaag52=1,xyzzyaagr1
blurb='      ! b_'//trim(i2s(xyzzyaaag52))//','//trim(i2s(xyzzyaaac52)&
&)
write(xyzzyaagw1,*)xyzzyaaag1(xyzzyaaag52,xyzzyaaac52),xyzzyaaay1(xyzz&
&yaaag52,xyzzyaaac52),trim(blurb)
enddo
enddo
write(xyzzyaagw1,*)'END Q TERM'
endif
if(xyzzyaafa1)then
write(xyzzyaagw1,*)'START BIEX1 TERM'
write(xyzzyaagw1,*)'Parameter value ; Optimizable (0=NO; 1=YES)'
do xyzzyaaah52=1,9
if(xyzzyaaij1(xyzzyaaah52)>=0)then
blurb='      ! c_'//trim(i2s(xyzzyaaah52))
write(xyzzyaagw1,*,iostat=xyzzyaaaa52)xyzzyaaih1(xyzzyaaah52),xyzzyaai&
&j1(xyzzyaaah52),trim(blurb)
endif
enddo
write(xyzzyaagw1,*)'END BIEX1 TERM'
endif
if(xyzzyaafb1)then
write(xyzzyaagw1,*)'START BIEX2 TERM'
write(xyzzyaagw1,*)'Parameter value ; Optimizable (0=NO; 1=YES)'
do xyzzyaaah52=1,2
write(xyzzyaagw1,*)xyzzyaail1(xyzzyaaah52),xyzzyaaio1(xyzzyaaah52),'      !&
& a_',trim(i2s(xyzzyaaah52))
write(xyzzyaagw1,*)xyzzyaaim1(xyzzyaaah52),xyzzyaaip1(xyzzyaaah52),'      !&
& b_',trim(i2s(xyzzyaaah52))
enddo
write(xyzzyaagw1,*)'END BIEX2 TERM'
endif
if(xyzzyaafc1)then
write(xyzzyaagw1,*)'START BIEX3 TERM'
write(xyzzyaagw1,*)'Parameter value ; Optimizable (0=NO; 1=YES)'
do xyzzyaaah52=1,9
if(xyzzyaaij1(xyzzyaaah52)>=0)then
blurb='      ! c_'//trim(i2s(xyzzyaaah52))
write(xyzzyaagw1,*,iostat=xyzzyaaaa52)xyzzyaaih1(xyzzyaaah52),xyzzyaai&
&j1(xyzzyaaah52),trim(blurb)
endif
enddo
write(xyzzyaagw1,*)'END BIEX3 TERM'
endif
if(xyzzyaafg1)then
write(xyzzyaagw1,*)'START EX2D TERM'
write(xyzzyaagw1,*)'Pair spin-dependence'
write(xyzzyaagw1,*)xyzzyaais1
write(xyzzyaagw1,*)'Single spin-dependence'
write(xyzzyaagw1,*)xyzzyaait1
write(xyzzyaagw1,*)'Parameter value ; Optimizable (0=NO; 1=YES)'
do xyzzyaaac52=1,xyzzyaair1
do xyzzyaaah52=1,xyzzyaaiq1
if(xyzzyaaik1(xyzzyaaah52,xyzzyaaac52)>=0)then
blurb='      ! c_'//trim(i2s(xyzzyaaah52))//','//trim(i2s(xyzzyaaac52)&
&)
write(xyzzyaagw1,*,iostat=xyzzyaaaa52)xyzzyaaii1(xyzzyaaah52,xyzzyaaac&
&52),xyzzyaaik1(xyzzyaaah52,xyzzyaaac52),trim(blurb)
endif
enddo
enddo
write(xyzzyaagw1,*)'END EX2D TERM'
endif
if(xyzzyaaff1)then
write(xyzzyaagw1,*)'START D TERM'
write(xyzzyaagw1,*)'Number of atom pairs'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaabm1))
write(xyzzyaagw1,*)'Number of sets'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaabn1))
write(xyzzyaagw1,*)'Expansion order N_d'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaacd1))
do xyzzyaaab52=1,xyzzyaabn1
write(xyzzyaagw1,*)'START SET ',trim(i2s(xyzzyaaab52))
write(xyzzyaagw1,*)'Pairs symmetric (0=NO; 1=YES)'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaabw1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Number of atom pairs in set'
write(xyzzyaagw1,'(3x,a)')trim(i2s(xyzzyaadk1(xyzzyaaab52)))
write(xyzzyaagw1,*)'Labels of the atom pairs in this set'
do xyzzyaaaj52=1,xyzzyaabm1
if(xyzzyaabv1(xyzzyaaaj52)==xyzzyaaab52)write(xyzzyaagw1,*)'    ',xyzz&
&yaabu1(xyzzyaaaj52,1:2)
enddo
write(xyzzyaagw1,*)'Cutoffs (a.u.)     ;  Optimizable (0=NO; 1=YES)'
if(xyzzyaabw1(xyzzyaaab52)==1)then
write(xyzzyaagw1,*)xyzzyaacw1(xyzzyaaab52,1,1),xyzzyaadh1(xyzzyaaab52,&
&1,1),' ! L_D_para_I,J (= J,I)'
write(xyzzyaagw1,*)xyzzyaacw1(xyzzyaaab52,2,1),xyzzyaadh1(xyzzyaaab52,&
&2,1),' ! L_D_perp_I,J (= J,I)'
else
write(xyzzyaagw1,*)xyzzyaacw1(xyzzyaaab52,1,1),xyzzyaadh1(xyzzyaaab52,&
&1,1),' ! L_D_para_I,J'
write(xyzzyaagw1,*)xyzzyaacw1(xyzzyaaab52,1,2),xyzzyaadh1(xyzzyaaab52,&
&1,2),' ! L_D_para_J,I'
write(xyzzyaagw1,*)xyzzyaacw1(xyzzyaaab52,2,1),xyzzyaadh1(xyzzyaaab52,&
&2,1),' ! L_D_perp_I,J'
write(xyzzyaagw1,*)xyzzyaacw1(xyzzyaaab52,2,2),xyzzyaadh1(xyzzyaaab52,&
&2,2),' ! L_D_perp_J,I'
endif
write(xyzzyaagw1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaai52=1,2
do xyzzyaaae52=0,xyzzyaacd1
do xyzzyaaaf52=xyzzyaabw1(xyzzyaaab52)*xyzzyaaae52,xyzzyaacd1
blurb='      ! D_('//trim(i2s(xyzzyaaab52))//')_'//xyzzyaamc1(xyzzyaaa&
&i52)//'_'//trim(i2s(xyzzyaaae52))//','//trim(i2s(xyzzyaaaf52))
if(xyzzyaaae52/=xyzzyaaaf52.and.xyzzyaabw1(xyzzyaaab52)==1)blurb=blurb&
&//' (='//trim(i2s(xyzzyaaaf52))//','//trim(i2s(xyzzyaaae52))//')'
write(xyzzyaagw1,*)xyzzyaaas1(xyzzyaaab52,xyzzyaaai52,xyzzyaaae52,xyzz&
&yaaaf52),xyzzyaabe1(xyzzyaaab52,xyzzyaaai52,xyzzyaaae52,xyzzyaaaf52),&
&trim(blurb)
enddo
enddo
enddo
write(xyzzyaagw1,*)'END SET ',trim(i2s(xyzzyaaab52))
enddo
write(xyzzyaagw1,*)'END D TERM'
endif
write(xyzzyaagw1,*)'END JASTROW'
write(xyzzyaagw1,*)
close(xyzzyaagw1)
endif
end subroutine write_pjastrow
subroutine xyzzyaaok1
implicit none
integer xyzzyaaaa53,xyzzyaaab53
if(.not.allocated(xyzzyaahu1))then
allocate(xyzzyaahs1(xyzzyaabj1),xyzzyaaht1(xyzzyaabj1),stat=xyzzyaaaa5&
&3)
call check_alloc(xyzzyaaaa53,'PREPARE_A_MATRICES','1')
do xyzzyaaab53=1,xyzzyaabj1
xyzzyaahs1(xyzzyaaab53)=xyzzyaaoq1(xyzzyaaab53)
xyzzyaaht1(xyzzyaaab53)=xyzzyaaor1(xyzzyaaab53)
enddo
xyzzyaahp1=maxval(xyzzyaahs1)
xyzzyaahq1=maxval(xyzzyaaht1)
allocate(xyzzyaahu1(xyzzyaahp1,xyzzyaahq1,xyzzyaabj1),xyzzyaahv1(xyzzy&
&aahq1),xyzzyaahr1(xyzzyaahp1,xyzzyaabj1),stat=xyzzyaaaa53)
call check_alloc(xyzzyaaaa53,'PREPARE_A_MATRICES','2')
endif
do xyzzyaaab53=1,xyzzyaabj1
call xyzzyaaol1(xyzzyaaab53)
call xyzzyaaom1(xyzzyaaab53)
enddo
end subroutine xyzzyaaok1
subroutine xyzzyaaol1(set)
implicit none
integer,intent(in) :: set
integer xyzzyaaaa54
real(dp),allocatable :: xyzzyaaab54(:,:)
allocate(xyzzyaaab54(xyzzyaahs1(set),xyzzyaaht1(set)),stat=xyzzyaaaa54&
&)
call check_alloc(xyzzyaaaa54,'PREPARE_A_MATRICES','3')
call xyzzyaaos1(set,xyzzyaahs1(set),xyzzyaaht1(set),xyzzyaaab54)
xyzzyaahu1(:,:,set)=0.d0
xyzzyaahu1(1:xyzzyaahs1(set),1:xyzzyaaht1(set),set)=xyzzyaaab54
deallocate(xyzzyaaab54)
end subroutine xyzzyaaol1
subroutine xyzzyaaom1(set)
implicit none
integer,intent(in) :: set
integer xyzzyaaaa55,xyzzyaaab55,xyzzyaaac55,xyzzyaaad55,xyzzyaaae55,xy&
&zzyaaaf55,xyzzyaaag55
logical,allocatable :: xyzzyaaah55(:)
allocate(xyzzyaaah55(xyzzyaaht1(set)),stat=xyzzyaaaa55)
call check_alloc(xyzzyaaaa55,'PREPARE_A_MATRICES','4')
xyzzyaaah55(:)=.false.
do xyzzyaaab55=1,xyzzyaahs1(set)
xyzzyaahr1(xyzzyaaab55,set)=-1
do xyzzyaaac55=xyzzyaaab55,xyzzyaaht1(set)
if(xyzzyaahu1(xyzzyaaab55,xyzzyaaac55,set)>0.5d0)then
xyzzyaahr1(xyzzyaaab55,set)=xyzzyaaac55
exit
endif
enddo
if(xyzzyaahr1(xyzzyaaab55,set)>0)xyzzyaaah55(xyzzyaahr1(xyzzyaaab55,se&
&t))=.true.
enddo
do xyzzyaaag55=1,xyzzyaafw1(set)
xyzzyaaab55=0
do xyzzyaaaf55=0,xyzzyaacf1(set)
do xyzzyaaae55=0,xyzzyaacg1(set)
do xyzzyaaad55=xyzzyaaae55,xyzzyaacg1(set)
xyzzyaaab55=xyzzyaaab55+1
if(xyzzyaaah55(xyzzyaaab55))then
if(xyzzyaaaw1(xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,set)==1)&
&xyzzyaaaw1(xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,set)=-2
if(xyzzyaaaw1(xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,set)==0)&
&xyzzyaaaw1(xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,set)=-1
else
if(xyzzyaaaw1(xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,set)==-2&
&)xyzzyaaaw1(xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,set)=1
if(xyzzyaaaw1(xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,set)==-1&
&)xyzzyaaaw1(xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,set)=0
endif
if(xyzzyaaad55>xyzzyaaae55)xyzzyaaaw1(xyzzyaaae55,xyzzyaaad55,xyzzyaaa&
&f55,xyzzyaaag55,set)=-1
enddo
enddo
enddo
enddo
deallocate(xyzzyaaah55)
end subroutine xyzzyaaom1
subroutine xyzzyaaon1
implicit none
integer xyzzyaaaa56,xyzzyaaab56
do xyzzyaaaa56=1,xyzzyaabj1
do xyzzyaaab56=1,xyzzyaafw1(xyzzyaaaa56)
call xyzzyaaoo1(xyzzyaaaa56,xyzzyaaab56)
enddo
enddo
end subroutine xyzzyaaon1
subroutine xyzzyaaoo1(set,s)
use slaarnabt, only : ddot
implicit none
integer,intent(in) :: set,s
integer xyzzyaaaa57,xyzzyaaab57,xyzzyaaac57,xyzzyaaad57,xyzzyaaae57
xyzzyaaaa57=0
do xyzzyaaad57=0,xyzzyaacf1(set)
do xyzzyaaac57=0,xyzzyaacg1(set)
do xyzzyaaab57=xyzzyaaac57,xyzzyaacg1(set)
xyzzyaaaa57=xyzzyaaaa57+1
xyzzyaahv1(xyzzyaaaa57)=xyzzyaaae1(xyzzyaaab57,xyzzyaaac57,xyzzyaaad57&
&,s,set)
enddo
enddo
enddo
do xyzzyaaaa57=1,xyzzyaahs1(set)
if(xyzzyaahr1(xyzzyaaaa57,set)>0)then
xyzzyaaae57=xyzzyaahr1(xyzzyaaaa57,set)+1
xyzzyaahv1(xyzzyaahr1(xyzzyaaaa57,set))=-ddot(xyzzyaaht1(set)-xyzzyaah&
&r1(xyzzyaaaa57,set),xyzzyaahu1(xyzzyaaaa57,xyzzyaaae57,set),xyzzyaahp&
&1,xyzzyaahv1(xyzzyaaae57),1)
endif
enddo
xyzzyaaaa57=0
do xyzzyaaad57=0,xyzzyaacf1(set)
do xyzzyaaac57=0,xyzzyaacg1(set)
xyzzyaaaa57=xyzzyaaaa57+1
xyzzyaaae1(xyzzyaaac57,xyzzyaaac57,xyzzyaaad57,s,set)=xyzzyaahv1(xyzzy&
&aaaa57)
do xyzzyaaab57=xyzzyaaac57+1,xyzzyaacg1(set)
xyzzyaaaa57=xyzzyaaaa57+1
xyzzyaaae1(xyzzyaaab57,xyzzyaaac57,xyzzyaaad57,s,set)=xyzzyaahv1(xyzzy&
&aaaa57)
xyzzyaaae1(xyzzyaaac57,xyzzyaaab57,xyzzyaaad57,s,set)=xyzzyaahv1(xyzzy&
&aaaa57)
enddo
enddo
enddo
end subroutine xyzzyaaoo1
subroutine xyzzyaaop1(set,no_gamma_vars,determined,no_free)
implicit none
integer,intent(in) :: set,no_gamma_vars
integer,intent(out) :: no_free
logical,intent(out) :: determined(no_gamma_vars)
integer xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58,xyzzyaaae58
real(dp),allocatable :: xyzzyaaaf58(:,:)
xyzzyaaae58=xyzzyaaoq1(set)
allocate(xyzzyaaaf58(xyzzyaaae58,no_gamma_vars),stat=xyzzyaaad58)
call check_alloc(xyzzyaaad58,'FIND_DETERMINED_GAMMA','')
call xyzzyaaos1(set,xyzzyaaae58,no_gamma_vars,xyzzyaaaf58)
no_free=no_gamma_vars
determined(:)=.false.
do xyzzyaaaa58=1,xyzzyaaae58
xyzzyaaac58=-1
do xyzzyaaab58=xyzzyaaaa58,no_gamma_vars
if(xyzzyaaaf58(xyzzyaaaa58,xyzzyaaab58)>0.5d0)then
xyzzyaaac58=xyzzyaaab58
exit
endif
enddo
if(xyzzyaaac58>0)then
determined(xyzzyaaac58)=.true.
no_free=no_free-1
endif
enddo
deallocate(xyzzyaaaf58)
end subroutine xyzzyaaop1
integer function xyzzyaaoq1(set)
implicit none
integer,intent(in) :: set
xyzzyaaoq1=2+3*xyzzyaacg1(set)+xyzzyaacf1(set)
if(xyzzyaagy1(set)==1)xyzzyaaoq1=xyzzyaaoq1+xyzzyaacf1(set)+1
if(xyzzyaagz1(set)==1)xyzzyaaoq1=xyzzyaaoq1+xyzzyaacg1(set)+1
end function xyzzyaaoq1
integer function xyzzyaaor1(set)
implicit none
integer,intent(in) :: set
xyzzyaaor1=(xyzzyaacf1(set)+1)*(xyzzyaacg1(set)+1+(xyzzyaacg1(set)*(xy&
&zzyaacg1(set)+1))/2)
end function xyzzyaaor1
subroutine xyzzyaaos1(set,xyzzyaahs1,no_gamma_vars,a)
use slaarnabt, only : reduced_echelon
implicit none
integer,intent(in) :: set,xyzzyaahs1,no_gamma_vars
real(dp),intent(out) :: a(xyzzyaahs1,no_gamma_vars)
integer xyzzyaaaa61,xyzzyaaab61,xyzzyaaac61,xyzzyaaad61,xyzzyaaae61,xy&
&zzyaaaf61,xyzzyaaag61
real(dp) xyzzyaaah61
xyzzyaaae61=2*xyzzyaacg1(set)+2
xyzzyaaah61=dble(xyzzyaagx1)
xyzzyaaaf61=3*xyzzyaacg1(set)+xyzzyaacf1(set)+3
xyzzyaaag61=3*xyzzyaacg1(set)+2*xyzzyaacf1(set)+4
a(1:xyzzyaahs1,1:no_gamma_vars)=0.d0
xyzzyaaaa61=0
do xyzzyaaad61=0,xyzzyaacf1(set)
do xyzzyaaac61=0,xyzzyaacg1(set)
do xyzzyaaab61=xyzzyaaac61,xyzzyaacg1(set)
xyzzyaaaa61=xyzzyaaaa61+1
if(xyzzyaaad61==1)then
if(xyzzyaaab61/=xyzzyaaac61)then
a(xyzzyaaab61+xyzzyaaac61+1,xyzzyaaaa61)=2.d0
else
a(xyzzyaaab61+xyzzyaaac61+1,xyzzyaaaa61)=1.d0
endif
endif
if(xyzzyaaac61==1)then
a(xyzzyaaab61+xyzzyaaad61+xyzzyaaae61,xyzzyaaaa61)=-xyzzyaacv1(set)
elseif(xyzzyaaac61==0)then
a(xyzzyaaab61+xyzzyaaad61+xyzzyaaae61,xyzzyaaaa61)=xyzzyaaah61
if(xyzzyaaab61==1)then
a(xyzzyaaad61+xyzzyaaae61,xyzzyaaaa61)=-xyzzyaacv1(set)
elseif(xyzzyaaab61==0)then
a(xyzzyaaad61+xyzzyaaae61,xyzzyaaaa61)=xyzzyaaah61
endif
if(xyzzyaagy1(set)==1)then
if(xyzzyaaab61==0)a(xyzzyaaaf61+xyzzyaaad61,xyzzyaaaa61)=1.d0
if(xyzzyaagz1(set)==1.and.xyzzyaaad61==0)a(xyzzyaaag61+xyzzyaaab61,xyz&
&zyaaaa61)=1.d0
else
if(xyzzyaagz1(set)==1.and.xyzzyaaad61==0)a(xyzzyaaaf61+xyzzyaaab61,xyz&
&zyaaaa61)=1.d0
endif
endif
enddo
enddo
enddo
call reduced_echelon(xyzzyaahs1,no_gamma_vars,a)
end subroutine xyzzyaaos1
logical function xyzzyaaot1(printout)
implicit none
logical,intent(in) :: printout
integer xyzzyaaaa62,xyzzyaaab62,l,m,n,xyzzyaaac62,xyzzyaaad62,xyzzyaaa&
&e62,xyzzyaaaf62
real(dp),parameter :: xyzzyaaag62=1.d-8
real(dp),allocatable :: xyzzyaaah62(:),xyzzyaaai62(:),xyzzyaaaj62(:)
logical :: xyzzyaaak62=.false.
character(80) tmpr
xyzzyaaot1=.true.
if(.not.am_master)return
xyzzyaaak62=.false.
do xyzzyaaaa62=1,xyzzyaabj1
allocate(xyzzyaaah62(0:2*xyzzyaacg1(xyzzyaaaa62)),xyzzyaaai62(0:xyzzya&
&acg1(xyzzyaaaa62)+xyzzyaacf1(xyzzyaaaa62)),xyzzyaaaj62(0:xyzzyaacg1(x&
&yzzyaaaa62)+xyzzyaacf1(xyzzyaaaa62)),stat=xyzzyaaaf62)
call check_alloc(xyzzyaaaf62,'TEST_GAMMA','')
do xyzzyaaab62=1,xyzzyaafw1(xyzzyaaaa62)
xyzzyaaah62(:)=0.d0
xyzzyaaai62(:)=0.d0
xyzzyaaaj62(:)=0.d0
do n=0,xyzzyaacf1(xyzzyaaaa62)
do m=0,xyzzyaacg1(xyzzyaaaa62)
do l=0,xyzzyaacg1(xyzzyaaaa62)
if(abs(xyzzyaaae1(m,l,n,xyzzyaaab62,xyzzyaaaa62)-xyzzyaaae1(l,m,n,xyzz&
&yaaab62,xyzzyaaaa62))>xyzzyaaag62)then
if(printout)then
call wout('Electron exchange symmetry not obeyed:')
call wout('  gamma_'//trim(i2s(l))//','//trim(i2s(m))//','//trim(i2s(n&
&))//','//trim(i2s(xyzzyaaab62))//','//trim(i2s(xyzzyaaaa62))//'  =  '&
&,xyzzyaaae1(l,m,n,xyzzyaaab62,xyzzyaaaa62),rfmt='(f21.12)')
call wout('  gamma_'//trim(i2s(m))//','//trim(i2s(l))//','//trim(i2s(n&
&))//','//trim(i2s(xyzzyaaab62))//','//trim(i2s(xyzzyaaaa62))//'  =  '&
&,xyzzyaaae1(m,l,n,xyzzyaaab62,xyzzyaaaa62),rfmt='(f21.12)')
endif
xyzzyaaak62=.true.
endif
xyzzyaaac62=l+m
xyzzyaaad62=m+n
xyzzyaaae62=l+n
if(n==1)xyzzyaaah62(xyzzyaaac62)=xyzzyaaah62(xyzzyaaac62)+xyzzyaaae1(l&
&,m,n,xyzzyaaab62,xyzzyaaaa62)
if(l==1)then
xyzzyaaai62(xyzzyaaad62)=xyzzyaaai62(xyzzyaaad62)+xyzzyaacv1(xyzzyaaaa&
&62)*xyzzyaaae1(l,m,n,xyzzyaaab62,xyzzyaaaa62)-dble(xyzzyaagx1)*xyzzya&
&aae1(0,m,n,xyzzyaaab62,xyzzyaaaa62)
endif
if(m==1)then
xyzzyaaaj62(xyzzyaaae62)=xyzzyaaaj62(xyzzyaaae62)+xyzzyaacv1(xyzzyaaaa&
&62)*xyzzyaaae1(l,m,n,xyzzyaaab62,xyzzyaaaa62)-dble(xyzzyaagx1)*xyzzya&
&aae1(l,0,n,xyzzyaaab62,xyzzyaaaa62)
endif
enddo
enddo
enddo
do xyzzyaaac62=0,2*xyzzyaacg1(xyzzyaaaa62)
if(abs(xyzzyaaah62(xyzzyaaac62))>xyzzyaaag62)then
if(printout)then
call wout('Condition on n=1 not obeyed.')
call wout('Set                       = '//trim(i2s(xyzzyaaaa62)))
call wout('Spin type                 = '//trim(i2s(xyzzyaaab62)))
call wout('k                         = '//trim(i2s(xyzzyaaac62)))
tmpr=r2s(xyzzyaaah62(xyzzyaaac62),'(f21.12)')
call wout('Sum_{l,m:l+m=k} gamma_lm1 = '//trim(tmpr))
endif
xyzzyaaak62=.true.
endif
enddo
do xyzzyaaad62=0,xyzzyaacf1(xyzzyaaaa62)+xyzzyaacg1(xyzzyaaaa62)
if(abs(xyzzyaaai62(xyzzyaaad62))>xyzzyaaag62)then
if(printout)then
call wout('Condition on l=1 not obeyed.')
call wout('Set                                       = '//trim(i2s(xyz&
&zyaaaa62)))
call wout('Spin type                                 = '//trim(i2s(xyz&
&zyaaab62)))
call wout('k                                         = '//trim(i2s(xyz&
&zyaaad62)))
tmpr=r2s(xyzzyaaai62(xyzzyaaad62),'(f21.12)')
call wout('Sum_{m,n:m+n=k} (L*gamma_1mn-2*gamma_0mn) = '//trim(tmpr))
endif
xyzzyaaak62=.true.
endif
enddo
do xyzzyaaae62=0,xyzzyaacf1(xyzzyaaaa62)+xyzzyaacg1(xyzzyaaaa62)
if(abs(xyzzyaaaj62(xyzzyaaae62))>xyzzyaaag62)then
if(printout)then
call wout('Condition on m=1 not obeyed.')
call wout('Set                                       = '//trim(i2s(xyz&
&zyaaaa62)))
call wout('Spin type                                 = '//trim(i2s(xyz&
&zyaaab62)))
call wout('k                                         = '//trim(i2s(xyz&
&zyaaae62)))
tmpr=r2s(xyzzyaaaj62(xyzzyaaae62),'(f21.12)')
call wout('Sum_{l,n:l+n=k} (L*gamma_l1n-2*gamma_l0n) = '//trim(tmpr))
endif
xyzzyaaak62=.true.
endif
enddo
enddo
deallocate(xyzzyaaah62,xyzzyaaai62,xyzzyaaaj62)
if(xyzzyaagy1(xyzzyaaaa62)==1)then
do xyzzyaaab62=1,xyzzyaafw1(xyzzyaaaa62)
do n=0,xyzzyaacf1(xyzzyaaaa62)
if(abs(xyzzyaaae1(0,0,n,xyzzyaaab62,xyzzyaaaa62))>xyzzyaaag62)then
if(printout)then
call wout('Condition that u must not be duplicated failed:')
call wout('Set       = '//trim(i2s(xyzzyaaaa62)))
call wout('Spin type = '//trim(i2s(xyzzyaaab62)))
call wout('n         = '//trim(i2s(n)))
tmpr=r2s(xyzzyaaae1(0,0,n,xyzzyaaab62,xyzzyaaaa62),'(f21.12)')
call wout('gamma_00n = '//trim(tmpr))
endif
xyzzyaaak62=.true.
endif
enddo
enddo
endif
if(xyzzyaagz1(xyzzyaaaa62)==1)then
do xyzzyaaab62=1,xyzzyaafw1(xyzzyaaaa62)
do l=0,xyzzyaacg1(xyzzyaaaa62)
if(abs(xyzzyaaae1(l,0,0,xyzzyaaab62,xyzzyaaaa62))>xyzzyaaag62)then
if(printout)then
call wout('Condition that chi must not be duplicated failed:')
call wout('Set       = '//trim(i2s(xyzzyaaaa62)))
call wout('Spin type = '//trim(i2s(xyzzyaaab62)))
call wout('l         = '//trim(i2s(l)))
tmpr=r2s(xyzzyaaae1(l,0,0,xyzzyaaab62,xyzzyaaaa62),'(f21.12)')
call wout('gamma_l00 = '//trim(tmpr))
endif
xyzzyaaak62=.true.
endif
enddo
enddo
endif
enddo
xyzzyaaot1=.not.xyzzyaaak62
end function xyzzyaaot1
subroutine xyzzyaaou1
implicit none
integer xyzzyaaaa63,xyzzyaaab63
if(.not.allocated(xyzzyaaib1))then
do xyzzyaaab63=1,xyzzyaafu1
xyzzyaahz1(xyzzyaaab63)=xyzzyaaoz1(xyzzyaaab63)
enddo
xyzzyaahx1=maxval(xyzzyaahz1)
allocate(xyzzyaaib1(xyzzyaahx1,xyzzyaahw1,xyzzyaafu1),xyzzyaaic1(xyzzy&
&aahw1),xyzzyaahy1(xyzzyaahx1,xyzzyaafu1),stat=xyzzyaaaa63)
call check_alloc(xyzzyaaaa63,'PREPARE_H_MATRICES','1')
xyzzyaaib1=0.d0
xyzzyaahy1=-1
endif
do xyzzyaaab63=1,xyzzyaafu1
call xyzzyaaov1(xyzzyaaab63)
call xyzzyaaow1(xyzzyaaab63)
enddo
end subroutine xyzzyaaou1
subroutine xyzzyaaov1(s)
implicit none
integer,intent(in) :: s
integer xyzzyaaaa64
real(dp),allocatable :: xyzzyaaab64(:,:)
allocate(xyzzyaaab64(xyzzyaahz1(s),xyzzyaahw1),stat=xyzzyaaaa64)
call check_alloc(xyzzyaaaa64,'PREPARE_H_MATRIX','2')
call xyzzyaapa1(s,xyzzyaahz1(s),xyzzyaaab64)
xyzzyaaib1(:,:,s)=0.d0
xyzzyaaib1(1:xyzzyaahz1(s),1:xyzzyaahw1,s)=xyzzyaaab64
deallocate(xyzzyaaab64)
end subroutine xyzzyaaov1
subroutine xyzzyaaow1(s)
implicit none
integer,intent(in) :: s
integer xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65,xyzzyaaad65,xyzzyaaae65,xy&
&zzyaaaf65
logical,allocatable :: xyzzyaaag65(:)
allocate(xyzzyaaag65(xyzzyaahw1),stat=xyzzyaaaa65)
call check_alloc(xyzzyaaaa65,'GET_H_PIVOT','3')
xyzzyaaag65(:)=.false.
do xyzzyaaab65=1,xyzzyaahz1(s)
do xyzzyaaac65=xyzzyaaab65,xyzzyaahw1
if(xyzzyaaib1(xyzzyaaab65,xyzzyaaac65,s)>0.5d0)then
xyzzyaahy1(xyzzyaaab65,s)=xyzzyaaac65
exit
endif
enddo
if(xyzzyaahy1(xyzzyaaab65,s)>0)xyzzyaaag65(xyzzyaahy1(xyzzyaaab65,s))=&
&.true.
enddo
xyzzyaaab65=0
do xyzzyaaaf65=0,xyzzyaacc1
do xyzzyaaae65=0,xyzzyaacc1
do xyzzyaaad65=0,xyzzyaacc1
xyzzyaaab65=xyzzyaaab65+1
if(xyzzyaaag65(xyzzyaaab65))then
if(xyzzyaabc1(xyzzyaaad65,xyzzyaaae65,xyzzyaaaf65,s)==1)xyzzyaabc1(xyz&
&zyaaad65,xyzzyaaae65,xyzzyaaaf65,s)=-2
if(xyzzyaabc1(xyzzyaaad65,xyzzyaaae65,xyzzyaaaf65,s)==0)xyzzyaabc1(xyz&
&zyaaad65,xyzzyaaae65,xyzzyaaaf65,s)=-1
else
if(xyzzyaabc1(xyzzyaaad65,xyzzyaaae65,xyzzyaaaf65,s)==-2)xyzzyaabc1(xy&
&zzyaaad65,xyzzyaaae65,xyzzyaaaf65,s)=1
if(xyzzyaabc1(xyzzyaaad65,xyzzyaaae65,xyzzyaaaf65,s)==-1)xyzzyaabc1(xy&
&zzyaaad65,xyzzyaaae65,xyzzyaaaf65,s)=0
endif
enddo
enddo
enddo
deallocate(xyzzyaaag65)
end subroutine xyzzyaaow1
subroutine xyzzyaaox1
use slaarnabt, only : ddot
implicit none
integer xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66,xyzzyaaad66,xyzzyaaae66,xy&
&zzyaaaf66
do xyzzyaaaa66=1,xyzzyaafu1
xyzzyaaab66=0
do xyzzyaaae66=0,xyzzyaacc1
do xyzzyaaad66=0,xyzzyaacc1
do xyzzyaaac66=0,xyzzyaacc1
xyzzyaaab66=xyzzyaaab66+1
xyzzyaaic1(xyzzyaaab66)=xyzzyaaar1(xyzzyaaac66,xyzzyaaad66,xyzzyaaae66&
&,xyzzyaaaa66)
enddo
enddo
enddo
do xyzzyaaab66=1,xyzzyaahz1(xyzzyaaaa66)
if(xyzzyaahy1(xyzzyaaab66,xyzzyaaaa66)>0)then
xyzzyaaaf66=xyzzyaahy1(xyzzyaaab66,xyzzyaaaa66)+1
xyzzyaaic1(xyzzyaahy1(xyzzyaaab66,xyzzyaaaa66))=-ddot(xyzzyaahw1-xyzzy&
&aahy1(xyzzyaaab66,xyzzyaaaa66),xyzzyaaib1(xyzzyaaab66,xyzzyaaaf66,xyz&
&zyaaaa66),xyzzyaahx1,xyzzyaaic1(xyzzyaaaf66),1)
endif
enddo
xyzzyaaab66=0
do xyzzyaaae66=0,xyzzyaacc1
do xyzzyaaad66=0,xyzzyaacc1
do xyzzyaaac66=0,xyzzyaacc1
xyzzyaaab66=xyzzyaaab66+1
xyzzyaaar1(xyzzyaaac66,xyzzyaaad66,xyzzyaaae66,xyzzyaaaa66)=xyzzyaaic1&
&(xyzzyaaab66)
enddo
enddo
enddo
enddo
end subroutine xyzzyaaox1
subroutine xyzzyaaoy1(determined)
implicit none
logical,intent(out) :: determined(xyzzyaahw1,xyzzyaafu1)
integer xyzzyaaaa67,xyzzyaaab67,xyzzyaaac67,xyzzyaaad67,xyzzyaaae67,xy&
&zzyaaaf67
real(dp),allocatable :: xyzzyaaag67(:,:)
determined=.false.
do xyzzyaaaa67=1,xyzzyaafu1
xyzzyaaaf67=xyzzyaaoz1(xyzzyaaaa67)
allocate(xyzzyaaag67(xyzzyaaaf67,xyzzyaahw1),stat=xyzzyaaae67)
call check_alloc(xyzzyaaae67,'FIND_DETERMINED_H','')
call xyzzyaapa1(xyzzyaaaa67,xyzzyaaaf67,xyzzyaaag67)
do xyzzyaaab67=1,xyzzyaaaf67
xyzzyaaad67=-1
do xyzzyaaac67=xyzzyaaab67,xyzzyaahw1
if(xyzzyaaag67(xyzzyaaab67,xyzzyaaac67)>0.5d0)then
xyzzyaaad67=xyzzyaaac67
exit
endif
enddo
if(xyzzyaaad67>0)determined(xyzzyaaad67,xyzzyaaaa67)=.true.
enddo
deallocate(xyzzyaaag67)
enddo
end subroutine xyzzyaaoy1
integer function xyzzyaaoz1(s)
implicit none
integer,intent(in) :: s
integer xyzzyaaaa68,xyzzyaaab68
xyzzyaaaa68=3*(2*xyzzyaacc1+1)
xyzzyaaab68=(xyzzyaacc1+1)**3
select case(xyzzyaaia1(s))
case(1)
continue
case(3)
xyzzyaaab68=xyzzyaaab68*2
case default
xyzzyaaab68=0
end select
xyzzyaaoz1=xyzzyaaaa68+xyzzyaaab68
end function xyzzyaaoz1
subroutine xyzzyaapa1(s,ncons_h,h)
use slaarnabt, only : reduced_echelon
implicit none
integer,intent(in) :: s,ncons_h
real(dp),intent(out) :: h(ncons_h,xyzzyaahw1)
integer xyzzyaaaa69,xyzzyaaab69,xyzzyaaac69,xyzzyaaad69,xyzzyaaae69,xy&
&zzyaaaf69,xyzzyaaag69,xyzzyaaah69,xyzzyaaai69,xyzzyaaaj69,xyzzyaaak69&
&,xyzzyaaal69
real(dp) xyzzyaaam69
xyzzyaaak69=xyzzyaacc1+1
xyzzyaaal69=xyzzyaaak69*xyzzyaaak69
xyzzyaaam69=dble(xyzzyaagx1)
h(1:ncons_h,1:xyzzyaahw1)=0.d0
do xyzzyaaaa69=0,2*xyzzyaacc1
xyzzyaaad69=1+xyzzyaaaa69
xyzzyaaae69=xyzzyaaad69+2*xyzzyaacc1+1
xyzzyaaaf69=xyzzyaaae69+2*xyzzyaacc1+1
do xyzzyaaaj69=max(0,xyzzyaaaa69-xyzzyaacc1),min(xyzzyaacc1,xyzzyaaaa6&
&9)
xyzzyaaai69=xyzzyaaaa69-xyzzyaaaj69
xyzzyaaab69=1+0+xyzzyaaai69*xyzzyaaak69+xyzzyaaaj69*xyzzyaaal69
h(xyzzyaaad69,xyzzyaaab69)=h(xyzzyaaad69,xyzzyaaab69)+xyzzyaaam69
xyzzyaaab69=1+1+xyzzyaaai69*xyzzyaaak69+xyzzyaaaj69*xyzzyaaal69
h(xyzzyaaad69,xyzzyaaab69)=h(xyzzyaaad69,xyzzyaaab69)-xyzzyaacm1
xyzzyaaab69=1+xyzzyaaai69+0*xyzzyaaak69+xyzzyaaaj69*xyzzyaaal69
h(xyzzyaaae69,xyzzyaaab69)=h(xyzzyaaae69,xyzzyaaab69)+xyzzyaaam69
xyzzyaaab69=1+xyzzyaaai69+1*xyzzyaaak69+xyzzyaaaj69*xyzzyaaal69
h(xyzzyaaae69,xyzzyaaab69)=h(xyzzyaaae69,xyzzyaaab69)-xyzzyaacm1
xyzzyaaab69=1+xyzzyaaai69+xyzzyaaaj69*xyzzyaaak69+0*xyzzyaaal69
h(xyzzyaaaf69,xyzzyaaab69)=h(xyzzyaaaf69,xyzzyaaab69)+xyzzyaaam69
xyzzyaaab69=1+xyzzyaaai69+xyzzyaaaj69*xyzzyaaak69+1*xyzzyaaal69
h(xyzzyaaaf69,xyzzyaaab69)=h(xyzzyaaaf69,xyzzyaaab69)-xyzzyaacm1
enddo
enddo
xyzzyaaag69=3*(2*xyzzyaacc1+1)
if(xyzzyaaia1(s)>0)then
do xyzzyaaaj69=0,xyzzyaacc1
do xyzzyaaai69=0,xyzzyaacc1
do xyzzyaaah69=0,xyzzyaacc1
xyzzyaaac69=1+xyzzyaaah69+xyzzyaaai69*xyzzyaaak69+xyzzyaaaj69*xyzzyaaa&
&l69
select case(xyzzyaaia1(s))
case(1)
xyzzyaaag69=xyzzyaaag69+1
xyzzyaaab69=1+xyzzyaaai69+xyzzyaaah69*xyzzyaaak69+xyzzyaaaj69*xyzzyaaa&
&l69
h(xyzzyaaag69,xyzzyaaac69)=h(xyzzyaaag69,xyzzyaaac69)+1.d0
h(xyzzyaaag69,xyzzyaaab69)=h(xyzzyaaag69,xyzzyaaab69)-1.d0
case(3)
xyzzyaaag69=xyzzyaaag69+1
xyzzyaaab69=1+xyzzyaaai69+xyzzyaaah69*xyzzyaaak69+xyzzyaaaj69*xyzzyaaa&
&l69
h(xyzzyaaag69,xyzzyaaac69)=h(xyzzyaaag69,xyzzyaaac69)+1.d0
h(xyzzyaaag69,xyzzyaaab69)=h(xyzzyaaag69,xyzzyaaab69)-1.d0
xyzzyaaag69=xyzzyaaag69+1
xyzzyaaab69=1+xyzzyaaah69+xyzzyaaaj69*xyzzyaaak69+xyzzyaaai69*xyzzyaaa&
&l69
h(xyzzyaaag69,xyzzyaaac69)=h(xyzzyaaag69,xyzzyaaac69)+1.d0
h(xyzzyaaag69,xyzzyaaab69)=h(xyzzyaaag69,xyzzyaaab69)-1.d0
end select
enddo
enddo
enddo
endif
if(xyzzyaaag69/=ncons_h)call errstop('CONSTRUCT_H','Number of equation&
&s doesn''t match size of matrix ('//trim(i2s(xyzzyaaag69))//'/='//tri&
&m(i2s(ncons_h))//').')
call reduced_echelon(ncons_h,xyzzyaahw1,h)
end subroutine xyzzyaapa1
logical function xyzzyaapb1(printout)
implicit none
logical,intent(in) :: printout
integer xyzzyaaaa70,xyzzyaaab70,xyzzyaaac70,xyzzyaaad70,xyzzyaaae70
real(dp) xyzzyaaaf70,xyzzyaaag70,xyzzyaaah70,xyzzyaaai70,xyzzyaaaj70,x&
&yzzyaaak70,xyzzyaaal70,xyzzyaaam70,xyzzyaaan70,xyzzyaaao70
real(dp),parameter :: xyzzyaaap70=1.d-8
logical :: xyzzyaaaq70=.false.
xyzzyaapb1=.true.
if(.not.am_master)return
xyzzyaaaf70=dble(xyzzyaagx1)
xyzzyaaaq70=.false.
do xyzzyaaaa70=1,xyzzyaafu1
do xyzzyaaab70=0,2*xyzzyaacc1
xyzzyaaaj70=0.d0
xyzzyaaak70=0.d0
xyzzyaaal70=0.d0
xyzzyaaam70=0.d0
xyzzyaaan70=0.d0
xyzzyaaao70=0.d0
do xyzzyaaae70=max(0,xyzzyaaab70-xyzzyaacc1),min(xyzzyaacc1,xyzzyaaab7&
&0)
xyzzyaaad70=xyzzyaaab70-xyzzyaaae70
xyzzyaaaj70=xyzzyaaaj70+xyzzyaaar1(0,xyzzyaaad70,xyzzyaaae70,xyzzyaaaa&
&70)
xyzzyaaak70=xyzzyaaak70+xyzzyaaar1(1,xyzzyaaad70,xyzzyaaae70,xyzzyaaaa&
&70)
xyzzyaaal70=xyzzyaaal70+xyzzyaaar1(xyzzyaaad70,0,xyzzyaaae70,xyzzyaaaa&
&70)
xyzzyaaam70=xyzzyaaam70+xyzzyaaar1(xyzzyaaad70,1,xyzzyaaae70,xyzzyaaaa&
&70)
xyzzyaaan70=xyzzyaaan70+xyzzyaaar1(xyzzyaaad70,xyzzyaaae70,0,xyzzyaaaa&
&70)
xyzzyaaao70=xyzzyaaao70+xyzzyaaar1(xyzzyaaad70,xyzzyaaae70,1,xyzzyaaaa&
&70)
enddo
xyzzyaaag70=xyzzyaaaf70*xyzzyaaaj70-xyzzyaacm1*xyzzyaaak70
xyzzyaaah70=xyzzyaaaf70*xyzzyaaal70-xyzzyaacm1*xyzzyaaam70
xyzzyaaai70=xyzzyaaaf70*xyzzyaaan70-xyzzyaacm1*xyzzyaaao70
if(abs(xyzzyaaag70)>xyzzyaaap70.or.abs(xyzzyaaah70)>xyzzyaaap70.or.abs&
&(xyzzyaaai70)>xyzzyaaap70)then
if(printout)then
call wout('No-cusp condition for H term not obeyed:')
call wout('Spin type                 = '//trim(i2s(xyzzyaaaa70)))
call wout('a                         = '//trim(i2s(xyzzyaaab70)))
call wout('Sum (xmn)                 = ',xyzzyaaag70)
call wout('Sum (mxn)                 = ',xyzzyaaah70)
call wout('Sum (mnx)                 = ',xyzzyaaai70)
endif
xyzzyaaaq70=.true.
endif
enddo
if(xyzzyaaia1(xyzzyaaaa70)==3)then
do xyzzyaaac70=0,xyzzyaacc1
do xyzzyaaad70=0,xyzzyaacc1
do xyzzyaaae70=0,xyzzyaacc1
if(abs(xyzzyaaar1(xyzzyaaac70,xyzzyaaad70,xyzzyaaae70,xyzzyaaaa70)-xyz&
&zyaaar1(xyzzyaaad70,xyzzyaaac70,xyzzyaaae70,xyzzyaaaa70))>xyzzyaaap70&
&.or.abs(xyzzyaaar1(xyzzyaaac70,xyzzyaaad70,xyzzyaaae70,xyzzyaaaa70)-x&
&yzzyaaar1(xyzzyaaae70,xyzzyaaad70,xyzzyaaac70,xyzzyaaaa70))>xyzzyaaap&
&70.or.abs(xyzzyaaar1(xyzzyaaac70,xyzzyaaad70,xyzzyaaae70,xyzzyaaaa70)&
&-xyzzyaaar1(xyzzyaaae70,xyzzyaaac70,xyzzyaaad70,xyzzyaaaa70))>xyzzyaa&
&ap70.or.abs(xyzzyaaar1(xyzzyaaac70,xyzzyaaad70,xyzzyaaae70,xyzzyaaaa7&
&0)-xyzzyaaar1(xyzzyaaac70,xyzzyaaae70,xyzzyaaad70,xyzzyaaaa70))>xyzzy&
&aaap70.or.abs(xyzzyaaar1(xyzzyaaac70,xyzzyaaad70,xyzzyaaae70,xyzzyaaa&
&a70)-xyzzyaaar1(xyzzyaaad70,xyzzyaaae70,xyzzyaaac70,xyzzyaaaa70))>xyz&
&zyaaap70)then
if(printout)then
call wout('Symmetry conditions for H term not obeyed.')
call wout('Spin type                 = '//trim(i2s(xyzzyaaaa70)))
call wout('l                         = '//trim(i2s(xyzzyaaac70)))
call wout('m                         = '//trim(i2s(xyzzyaaad70)))
call wout('n                         = '//trim(i2s(xyzzyaaae70)))
endif
xyzzyaaaq70=.true.
endif
enddo
enddo
enddo
elseif(xyzzyaaia1(xyzzyaaaa70)==1)then
do xyzzyaaac70=0,xyzzyaacc1
do xyzzyaaad70=0,xyzzyaacc1
do xyzzyaaae70=0,xyzzyaacc1
xyzzyaaag70=xyzzyaaar1(xyzzyaaac70,xyzzyaaad70,xyzzyaaae70,xyzzyaaaa70&
&)-xyzzyaaar1(xyzzyaaad70,xyzzyaaac70,xyzzyaaae70,xyzzyaaaa70)
if(abs(xyzzyaaag70)>xyzzyaaap70)then
if(printout)then
call wout('Symmetry conditions for H term not obeyed.')
call wout('Spin type                 = '//trim(i2s(xyzzyaaaa70)))
call wout('l                         = '//trim(i2s(xyzzyaaac70)))
call wout('m                         = '//trim(i2s(xyzzyaaad70)))
call wout('n                         = '//trim(i2s(xyzzyaaae70)))
call wout('h_lmn - h_mln             = ',abs(xyzzyaaag70))
endif
xyzzyaaaq70=.true.
endif
enddo
enddo
enddo
endif
enddo
xyzzyaapb1=.not.xyzzyaaaq70
end function xyzzyaapb1
subroutine xyzzyaapc1
use slaarnaan, only : ee_kato_gamma
use slaarnaas, only : harmwire_b
implicit none
integer xyzzyaaaa71,xyzzyaaab71,xyzzyaaac71
real(dp) xyzzyaaad71
if(xyzzyaagx1<=0)return
xyzzyaaad71=1.d0/(dble(xyzzyaagx1)*(-xyzzyaach1)**(xyzzyaagx1-1))
do xyzzyaaab71=1,nspin
do xyzzyaaac71=xyzzyaaab71,nspin
xyzzyaaaa71=which_spair(xyzzyaaab71,xyzzyaaac71,xyzzyaafh1)
if(.not.ee_cusp_in_orbital(xyzzyaaab71,xyzzyaaac71).and..not.xyzzyaafa&
&1.and..not.xyzzyaafc1.and..not.xyzzyaafg1.and.harmwire_b<0.d0)then
xyzzyaaab1(0,xyzzyaaaa71)=ee_kato_gamma(xyzzyaaab71,xyzzyaaac71,.not.n&
&oncoll_spin.and.(xyzzyaafh1/=0.or.ferromagnetic))*xyzzyaaad71
else
xyzzyaaab1(0,xyzzyaaaa71)=0.d0
endif
enddo
enddo
end subroutine xyzzyaapc1
subroutine xyzzyaapd1
use slaarnaan, only : ee_kato_gamma
use slaarnabt, only : dcopy
use slaarnaas, only : harmwire_b
implicit none
integer xyzzyaaaa72,xyzzyaaab72,xyzzyaaac72,xyzzyaaad72,xyzzyaaae72
real(dp) xyzzyaaaf72,xyzzyaaag72
xyzzyaaae72=xyzzyaabx1+1
if(xyzzyaaeq1)then
do xyzzyaaaa72=1,xyzzyaafp1
xyzzyaaab1(1,xyzzyaaaa72)=xyzzyaaab1(0,xyzzyaaaa72)*dble(xyzzyaagx1)/x&
&yzzyaach1
enddo
xyzzyaaag72=1.d0/(-xyzzyaach1)**xyzzyaagx1
endif
do xyzzyaaac72=1,nspin
do xyzzyaaad72=xyzzyaaac72,nspin
xyzzyaaaa72=which_spair(xyzzyaaac72,xyzzyaaad72,xyzzyaafh1)
if(xyzzyaaeq1)then
xyzzyaaab72=which_spair(xyzzyaaac72,xyzzyaaad72,levels_spairs)
call dcopy(xyzzyaaae72,xyzzyaaab1(0,xyzzyaaaa72),1,xyzzyaaaj1(0,xyzzya&
&aab72),1)
endif
xyzzyaaaf72=0.d0
if(.not.ee_cusp_in_orbital(xyzzyaaac72,xyzzyaaad72).and..not.xyzzyaafa&
&1.and..not.xyzzyaafc1.and..not.xyzzyaafg1.and.harmwire_b<0.d0)xyzzyaa&
&af72=ee_kato_gamma(xyzzyaaac72,xyzzyaaad72,.not.noncoll_spin.and.(xyz&
&zyaafh1/=0.or.ferromagnetic))
if(xyzzyaaeq1)then
xyzzyaaaj1(1,xyzzyaaab72)=xyzzyaaaj1(1,xyzzyaaab72)+xyzzyaaaf72*xyzzya&
&aag72
else
if(dimensionality==3)then
xyzzyaaan1(xyzzyaaaa72)=2*xyzzyaaaf72*xyzzyaach1*xyzzyaaao1(xyzzyaaaa7&
&2)**2/(xyzzyaach1+2*xyzzyaagx1*xyzzyaaao1(xyzzyaaaa72))
elseif(dimensionality==2)then
xyzzyaaan1(xyzzyaaaa72)=3*xyzzyaaaf72*xyzzyaach1*xyzzyaaao1(xyzzyaaaa7&
&2)**1.5d0/(xyzzyaach1+3*xyzzyaagx1*xyzzyaaao1(xyzzyaaaa72))
else
call errstop('CONSTRUCT_U','RPA Jastrow not implemented in 1D.')
endif
endif
enddo
enddo
end subroutine xyzzyaapd1
subroutine xyzzyaape1
use slaarnabt, only : dcopy,dscal
implicit none
integer xyzzyaaaa73,xyzzyaaab73,xyzzyaaac73,xyzzyaaad73,xyzzyaaae73
real(dp) xyzzyaaaf73,xyzzyaaag73
xyzzyaacr1=(/sin(xyzzyaacs1)*cos(xyzzyaact1),sin(xyzzyaacs1)*sin(xyzzy&
&aact1),cos(xyzzyaacs1)/)
xyzzyaacj1=xyzzyaaci1**2
if(xyzzyaaci1<=0.d0)call errstop('CONSTRUCT_2D','L_ucylrho<=0.')
xyzzyaaaf73=dble(xyzzyaagx1)/xyzzyaaci1
if(xyzzyaack1<=0.d0)call errstop('CONSTRUCT_2D','L_ucylz<=0.')
xyzzyaaag73=dble(xyzzyaagx1)/xyzzyaack1
do xyzzyaaaa73=1,xyzzyaafq1
call dcopy(xyzzyaabz1+1,xyzzyaaac1(0,0,xyzzyaaaa73),xyzzyaaby1+1,xyzzy&
&aaac1(1,0,xyzzyaaaa73),xyzzyaaby1+1)
call dscal(xyzzyaabz1+1,xyzzyaaaf73,xyzzyaaac1(1,0,xyzzyaaaa73),xyzzya&
&aby1+1)
call dcopy(xyzzyaaby1+1,xyzzyaaac1(0,0,xyzzyaaaa73),1,xyzzyaaac1(0,1,x&
&yzzyaaaa73),1)
call dscal(xyzzyaaby1+1,xyzzyaaag73,xyzzyaaac1(0,1,xyzzyaaaa73),1)
enddo
xyzzyaaae73=(xyzzyaaby1+1)*(xyzzyaabz1+1)
do xyzzyaaac73=1,nspin
do xyzzyaaad73=xyzzyaaac73,nspin
xyzzyaaaa73=which_spair(xyzzyaaad73,xyzzyaaac73,xyzzyaafi1)
xyzzyaaab73=which_spair(xyzzyaaad73,xyzzyaaac73,levels_spairs)
call dcopy(xyzzyaaae73,xyzzyaaac1(0,0,xyzzyaaaa73),1,xyzzyaaak1(0,0,xy&
&zzyaaab73),1)
enddo
enddo
end subroutine xyzzyaape1
subroutine xyzzyaapf1
implicit none
integer xyzzyaaaa74,xyzzyaaab74,xyzzyaaac74,xyzzyaaad74
do xyzzyaaac74=1,nspin
do xyzzyaaad74=xyzzyaaac74,nspin
xyzzyaaaa74=which_spair(xyzzyaaac74,xyzzyaaad74,xyzzyaafl1)
xyzzyaaab74=which_spair(xyzzyaaac74,xyzzyaaad74,levels_spairs)
xyzzyaaaq1(:,xyzzyaaab74)=xyzzyaaap1(:,xyzzyaaaa74)
enddo
enddo
end subroutine xyzzyaapf1
subroutine xyzzyaapg1
use slaarnaan, only : en_kato_gamma
use slaarnabt, only : dcopy
implicit none
integer xyzzyaaaa75,xyzzyaaab75,xyzzyaaac75,xyzzyaaad75
real(dp) xyzzyaaae75,xyzzyaaaf75
do xyzzyaaaa75=1,xyzzyaabi1
xyzzyaaaf75=dble(xyzzyaagx1)/xyzzyaacu1(xyzzyaaaa75)
do xyzzyaaab75=1,xyzzyaafv1(xyzzyaaaa75)
xyzzyaaad1(1,xyzzyaaab75,xyzzyaaaa75)=xyzzyaaad1(0,xyzzyaaab75,xyzzyaa&
&aa75)*xyzzyaaaf75
enddo
xyzzyaaad75=xyzzyaace1(xyzzyaaaa75)+1
do xyzzyaaac75=1,nspin
xyzzyaaab75=which_ssingle(xyzzyaaac75,xyzzyaafn1(xyzzyaaaa75))
call dcopy(xyzzyaaad75,xyzzyaaad1(0,xyzzyaaab75,xyzzyaaaa75),1,xyzzyaa&
&al1(0,xyzzyaaac75,xyzzyaaaa75),1)
enddo
if(xyzzyaadn1(xyzzyaaaa75)==1)then
do xyzzyaaac75=1,nspin
xyzzyaaae75=en_kato_gamma(xyzzyaaac75,xyzzyaafx1(xyzzyaaaa75))
xyzzyaaal1(1,xyzzyaaac75,xyzzyaaaa75)=xyzzyaaal1(1,xyzzyaaac75,xyzzyaa&
&aa75)+xyzzyaaae75/(-xyzzyaacu1(xyzzyaaaa75))**xyzzyaagx1
enddo
endif
enddo
end subroutine xyzzyaapg1
subroutine xyzzyaaph1
implicit none
integer xyzzyaaaa76
do xyzzyaaaa76=1,xyzzyaabj1
call xyzzyaapi1(xyzzyaaaa76)
enddo
end subroutine xyzzyaaph1
subroutine xyzzyaapi1(set,s_only)
implicit none
integer,intent(in) :: set
integer,intent(in),optional :: s_only
integer xyzzyaaaa77,xyzzyaaab77,xyzzyaaac77,xyzzyaaad77,xyzzyaaae77,xy&
&zzyaaaf77,xyzzyaaag77,xyzzyaaah77
do xyzzyaaag77=1,nspin
do xyzzyaaah77=xyzzyaaag77,nspin
xyzzyaaaa77=which_spair(xyzzyaaah77,xyzzyaaag77,xyzzyaafo1(set))
if(present(s_only))then
if(xyzzyaaaa77/=s_only)cycle
endif
xyzzyaaab77=which_spair(xyzzyaaah77,xyzzyaaag77,levels_spairs)
xyzzyaaaf77=0
do xyzzyaaae77=0,xyzzyaacf1(set)
do xyzzyaaad77=0,xyzzyaacg1(set)
do xyzzyaaac77=0,xyzzyaacg1(set)
xyzzyaaaf77=xyzzyaaaf77+1
xyzzyaaam1(xyzzyaaaf77,xyzzyaaab77,set)=xyzzyaaae1(xyzzyaaac77,xyzzyaa&
&ad77,xyzzyaaae77,xyzzyaaaa77,set)
enddo
enddo
enddo
enddo
enddo
end subroutine xyzzyaapi1
subroutine xyzzyaapj1
implicit none
integer xyzzyaaaa78,xyzzyaaab78,xyzzyaaac78,xyzzyaaad78(3)
do xyzzyaaaa78=1,xyzzyaafy1-1
do xyzzyaaab78=xyzzyaaaa78+1,xyzzyaafy1
if(xyzzyaagb1(xyzzyaaaa78)>xyzzyaagb1(xyzzyaaab78))then
xyzzyaaac78=xyzzyaagb1(xyzzyaaaa78)
xyzzyaagb1(xyzzyaaaa78)=xyzzyaagb1(xyzzyaaab78)
xyzzyaagb1(xyzzyaaab78)=xyzzyaaac78
xyzzyaaad78(1:3)=xyzzyaaga1(1:3,xyzzyaaaa78)
xyzzyaaga1(1:3,xyzzyaaaa78)=xyzzyaaga1(1:3,xyzzyaaab78)
xyzzyaaga1(1:3,xyzzyaaab78)=xyzzyaaad78(1:3)
endif
enddo
enddo
xyzzyaagt1(1:xyzzyaagq1)=0
do xyzzyaaaa78=1,xyzzyaafy1
xyzzyaagg1(1:3,xyzzyaaaa78)=xyzzyaaga1(1,xyzzyaaaa78)*b1(1:3)+xyzzyaag&
&a1(2,xyzzyaaaa78)*b2(1:3)+xyzzyaaga1(3,xyzzyaaaa78)*b3(1:3)
xyzzyaagh1(xyzzyaaaa78)=xyzzyaagg1(1,xyzzyaaaa78)**2+xyzzyaagg1(2,xyzz&
&yaaaa78)**2+xyzzyaagg1(3,xyzzyaaaa78)**2
xyzzyaagt1(xyzzyaagb1(xyzzyaaaa78))=xyzzyaagt1(xyzzyaagb1(xyzzyaaaa78)&
&)+1
enddo
end subroutine xyzzyaapj1
subroutine xyzzyaapk1
implicit none
integer xyzzyaaaa79,xyzzyaaab79,xyzzyaaac79,xyzzyaaad79,xyzzyaaae79(3)
do xyzzyaaaa79=1,xyzzyaafz1-1
do xyzzyaaab79=xyzzyaaaa79+1,xyzzyaafz1
if(xyzzyaagd1(xyzzyaaaa79)>xyzzyaagd1(xyzzyaaab79))then
xyzzyaaad79=xyzzyaagd1(xyzzyaaaa79)
xyzzyaagd1(xyzzyaaaa79)=xyzzyaagd1(xyzzyaaab79)
xyzzyaagd1(xyzzyaaab79)=xyzzyaaad79
xyzzyaaae79(1:3)=xyzzyaagc1(1:3,xyzzyaaaa79)
xyzzyaagc1(1:3,xyzzyaaaa79)=xyzzyaagc1(1:3,xyzzyaaab79)
xyzzyaagc1(1:3,xyzzyaaab79)=xyzzyaaae79(1:3)
endif
enddo
enddo
xyzzyaaab79=1
xyzzyaaac79=-xyzzyaagr1-1
do xyzzyaaaa79=1,xyzzyaafz1
if(xyzzyaagd1(xyzzyaaaa79)/=xyzzyaaac79)then
xyzzyaaac79=xyzzyaagd1(xyzzyaaaa79)
xyzzyaagv1(xyzzyaaac79)=xyzzyaaab79
xyzzyaaab79=xyzzyaaab79+1
endif
enddo
if(xyzzyaaab79-1/=xyzzyaags1)call errstop('CONSTRUCT_GVECS_Q','Bug.')
xyzzyaagu1(1:xyzzyaags1)=0
do xyzzyaaaa79=1,xyzzyaafz1
xyzzyaagi1(1:3,xyzzyaaaa79)=xyzzyaagc1(1,xyzzyaaaa79)*pb1(1:3)+xyzzyaa&
&gc1(2,xyzzyaaaa79)*pb2(1:3)+xyzzyaagc1(3,xyzzyaaaa79)*pb3(1:3)
xyzzyaagj1(xyzzyaaaa79)=xyzzyaagi1(1,xyzzyaaaa79)**2+xyzzyaagi1(2,xyzz&
&yaaaa79)**2+xyzzyaagi1(3,xyzzyaaaa79)**2
xyzzyaagu1(xyzzyaagv1(xyzzyaagd1(xyzzyaaaa79)))=xyzzyaagu1(xyzzyaagv1(&
&xyzzyaagd1(xyzzyaaaa79)))+1
enddo
end subroutine xyzzyaapk1
subroutine xyzzyaapl1
use slaarnabt, only : dcopy
implicit none
integer xyzzyaaaa80,xyzzyaaab80,xyzzyaaac80,xyzzyaaad80
do xyzzyaaac80=1,nspin
do xyzzyaaad80=xyzzyaaac80,nspin
xyzzyaaaa80=which_spair(xyzzyaaac80,xyzzyaaad80,xyzzyaafj1)
xyzzyaaab80=which_spair(xyzzyaaac80,xyzzyaaad80,levels_spairs)
call dcopy(xyzzyaagq1,xyzzyaaaf1(1,xyzzyaaaa80),1,xyzzyaaah1(1,xyzzyaa&
&ab80),1)
enddo
enddo
end subroutine xyzzyaapl1
subroutine xyzzyaapm1
implicit none
integer xyzzyaaaa81,xyzzyaaab81,xyzzyaaac81
do xyzzyaaac81=1,nspin
xyzzyaaaa81=which_ssingle(xyzzyaaac81,xyzzyaafk1)
do xyzzyaaab81=-xyzzyaagr1,-1
if(xyzzyaagv1(xyzzyaaab81)/=0)xyzzyaaai1(xyzzyaagv1(xyzzyaaab81),xyzzy&
&aaac81)=-xyzzyaaag1(-xyzzyaaab81,xyzzyaaaa81)
enddo
do xyzzyaaab81=1,xyzzyaagr1
xyzzyaaai1(xyzzyaagv1(xyzzyaaab81),xyzzyaaac81)=xyzzyaaag1(xyzzyaaab81&
&,xyzzyaaaa81)
enddo
enddo
end subroutine xyzzyaapm1
subroutine xyzzyaapn1(rele,eevecs,eivecs,value_jas)
use slaarnabt, only : dcopy
implicit none
real(dp),intent(in) :: rele(3,netot),eevecs(4,netot,netot),eivecs(4,ni&
&tot,netot)
real(dp),intent(out) :: value_jas
integer xyzzyaaaa82,xyzzyaaab82,xyzzyaaac82,xyzzyaaad82,xyzzyaaae82,xy&
&zzyaaaf82,xyzzyaaag82,xyzzyaaah82,xyzzyaaai82,xyzzyaaaj82,xyzzyaaak82&
&,xyzzyaaal82,xyzzyaaam82
integer,save :: xyzzyaaan82
real(dp) xyzzyaaao82,xyzzyaaap82,xyzzyaaaq82,xyzzyaaar82,xyzzyaaas82,x&
&yzzyaaat82,xyzzyaaau82,xyzzyaaav82,xyzzyaaaw82,xyzzyaaax82,xyzzyaaay8&
&2,xyzzyaaaz82,xyzzyaaba82,xyzzyaabb82,xyzzyaabc82,xyzzyaabd82,xyzzyaa&
&be82,xyzzyaabf82,xyzzyaabg82,xyzzyaabh82,xyzzyaabi82,xyzzyaabj82(3),x&
&yzzyaabk82(3),xyzzyaabl82,xyzzyaabm82,xyzzyaabn82(3),xyzzyaabo82(4),x&
&yzzyaabp82(4),xyzzyaabq82,xyzzyaabr82,xyzzyaabs82,xyzzyaabt82,xyzzyaa&
&bu82,xyzzyaabv82,xyzzyaabw82,xyzzyaabx82,xyzzyaaby82,xyzzyaabz82
real(dp),allocatable,save :: xyzzyaaca82(:,:),xyzzyaacb82(:)
logical,save :: xyzzyaacc82=.true.,xyzzyaacd82=.false.
if(xyzzyaacc82)then
if(xyzzyaaep1.or.xyzzyaaes1.or.xyzzyaaex1)then
xyzzyaacd82=.true.
xyzzyaaan82=netot**2
allocate(xyzzyaaca82(netot,netot),xyzzyaacb82(netot),stat=xyzzyaaak82)
call check_alloc(xyzzyaaak82,'FULL_JASTROW_VAL','')
endif
xyzzyaacc82=.false.
endif
call timer('JASTROW',.true.)
xyzzyaabq82=0.d0
xyzzyaabt82=0.d0
xyzzyaabu82=0.d0
xyzzyaabv82=0.d0
xyzzyaabw82=0.d0
xyzzyaabx82=0.d0
xyzzyaaby82=0.d0
xyzzyaabz82=0.d0
xyzzyaabs82=0.d0
xyzzyaabr82=0.d0
if(xyzzyaaeu1)call xyzzyaapw1(0,eevecs,.false.,.false.)
if(xyzzyaaex1)then
do xyzzyaaae82=1,nitot
xyzzyaaah82=xyzzyaabp1(xyzzyaaae82)
if(xyzzyaaah82==0)cycle
xyzzyaaag82=xyzzyaacg1(xyzzyaaah82)
xyzzyaabl82=xyzzyaacv1(xyzzyaaah82)
do xyzzyaaaa82=1,netot
xyzzyaaax82=eivecs(4,xyzzyaaae82,xyzzyaaaa82)
if(xyzzyaaax82>=xyzzyaabl82)cycle
xyzzyaabm82=xyzzyaaax82
xyzzyaads1(1,xyzzyaaaa82,xyzzyaaae82)=xyzzyaabm82
do xyzzyaaaf82=2,xyzzyaaag82
xyzzyaabm82=xyzzyaabm82*xyzzyaaax82
xyzzyaads1(xyzzyaaaf82,xyzzyaaaa82,xyzzyaaae82)=xyzzyaabm82
enddo
enddo
enddo
endif
if(xyzzyaacd82)then
call dcopy(xyzzyaaan82,eevecs(4,1,1),4,xyzzyaaca82(1,1),1)
xyzzyaaaj82=1
xyzzyaaai82=netot
endif
do xyzzyaaaa82=1,netot
xyzzyaaab82=which_spin(xyzzyaaaa82)
if(xyzzyaacd82)then
xyzzyaaaj82=xyzzyaaaj82+1
xyzzyaaai82=xyzzyaaai82-1
if(xyzzyaaai82>0)call dcopy(xyzzyaaai82,xyzzyaaca82(xyzzyaaaj82,xyzzya&
&aaa82),1,xyzzyaacb82(xyzzyaaaj82),1)
endif
if(xyzzyaaep1)then
do xyzzyaaac82=xyzzyaaaa82+1,netot
xyzzyaaad82=which_spin(xyzzyaaac82)
xyzzyaaaw82=xyzzyaacb82(xyzzyaaac82)
call xyzzyaapr1(xyzzyaaaw82,xyzzyaaab82,xyzzyaaad82,.true.,.false.,.fa&
&lse.,xyzzyaaao82,xyzzyaabb82,xyzzyaabc82)
xyzzyaabq82=xyzzyaabq82+xyzzyaaao82
enddo
endif
if(xyzzyaaes1)then
do xyzzyaaac82=xyzzyaaaa82+1,netot
xyzzyaaad82=which_spin(xyzzyaaac82)
xyzzyaaba82=dot_product(eevecs(1:3,xyzzyaaac82,xyzzyaaaa82),xyzzyaacr1&
&)
xyzzyaaaz82=xyzzyaacb82(xyzzyaaac82)*xyzzyaacb82(xyzzyaaac82)-xyzzyaab&
&a82*xyzzyaaba82
if(xyzzyaaaz82<xyzzyaacj1.or.xyzzyaagx1==0)then
xyzzyaaaz82=sqrt(max(xyzzyaaaz82,0.d0))
call xyzzyaapu1(xyzzyaaaz82,abs(xyzzyaaba82),xyzzyaaab82,xyzzyaaad82,.&
&true.,.false.,.false.,xyzzyaaap82,xyzzyaabb82,xyzzyaabc82,xyzzyaabd82&
&,xyzzyaabe82)
xyzzyaabr82=xyzzyaabr82+xyzzyaaap82
endif
enddo
endif
if(xyzzyaaet1)then
do xyzzyaaac82=xyzzyaaaa82+1,netot
xyzzyaaad82=which_spin(xyzzyaaac82)
xyzzyaaaw82=xyzzyaacb82(xyzzyaaac82)
call xyzzyaapv1(xyzzyaaaw82,xyzzyaaab82,xyzzyaaad82,.true.,.false.,.fa&
&lse.,xyzzyaaaq82,xyzzyaabb82,xyzzyaabc82)
xyzzyaabs82=xyzzyaabs82+xyzzyaaaq82
enddo
endif
if(xyzzyaaeu1)then
if(xyzzyaaie1(xyzzyaaaa82))then
xyzzyaabb82=xyzzyaael1(xyzzyaaaa82)
xyzzyaabt82=xyzzyaabt82+xyzzyaabb82*xyzzyaabb82-2*xyzzyaaed1(xyzzyaaaa&
&82)-3*xyzzyaaee1(xyzzyaaaa82)
endif
endif
if(xyzzyaaev1)then
call xyzzyaapy1(xyzzyaaaa82,eevecs,.true.,.false.,.false.,xyzzyaabb82,&
&xyzzyaabj82,xyzzyaabc82)
xyzzyaabt82=xyzzyaabt82+xyzzyaabb82
endif
if(xyzzyaaew1)then
do xyzzyaaae82=1,nitot
xyzzyaaah82=xyzzyaabo1(xyzzyaaae82)
if(xyzzyaaah82==0)cycle
xyzzyaaax82=eivecs(4,xyzzyaaae82,xyzzyaaaa82)
call xyzzyaaqa1(xyzzyaaax82,xyzzyaaab82,xyzzyaaah82,.true.,.false.,.fa&
&lse.,xyzzyaaar82,xyzzyaabb82,xyzzyaabc82)
xyzzyaabu82=xyzzyaabu82+xyzzyaaar82
enddo
endif
if(xyzzyaaex1.and.xyzzyaaaa82<netot)then
do xyzzyaaae82=1,nitot
xyzzyaaah82=xyzzyaabp1(xyzzyaaae82)
if(xyzzyaaah82==0)cycle
xyzzyaaax82=eivecs(4,xyzzyaaae82,xyzzyaaaa82)
xyzzyaabl82=xyzzyaacv1(xyzzyaaah82)
if(xyzzyaaax82>=xyzzyaabl82)cycle
xyzzyaaag82=xyzzyaacg1(xyzzyaaah82)
call dcopy(xyzzyaaca1,xyzzyaads1(1,1,xyzzyaaae82),1,xyzzyaadq1(1,1),1)
call dcopy(xyzzyaaag82,xyzzyaadq1(1,xyzzyaaaa82),1,xyzzyaado1(1),1)
do xyzzyaaac82=xyzzyaaaa82+1,netot
xyzzyaaad82=which_spin(xyzzyaaac82)
xyzzyaaay82=eivecs(4,xyzzyaaae82,xyzzyaaac82)
xyzzyaaaw82=xyzzyaacb82(xyzzyaaac82)
call xyzzyaaqb1(xyzzyaaac82,xyzzyaaax82,xyzzyaaay82,xyzzyaaaw82,xyzzya&
&aab82,xyzzyaaad82,xyzzyaaah82,.true.,.false.,.false.,.false.,xyzzyaaa&
&s82,xyzzyaabb82,xyzzyaabc82,xyzzyaabd82,xyzzyaabe82,xyzzyaabf82,xyzzy&
&aabg82,xyzzyaabh82,xyzzyaabi82)
xyzzyaabv82=xyzzyaabv82+xyzzyaaas82
enddo
enddo
endif
if(xyzzyaaey1)then
do xyzzyaaac82=xyzzyaaaa82+1,netot
xyzzyaaad82=which_spin(xyzzyaaac82)
xyzzyaabk82=eevecs(1:3,xyzzyaaac82,xyzzyaaaa82)
call xyzzyaaqd1(xyzzyaabk82,xyzzyaaab82,xyzzyaaad82,.true.,.false.,.fa&
&lse.,xyzzyaaat82,xyzzyaabj82,xyzzyaabb82)
xyzzyaabw82=xyzzyaabw82+xyzzyaaat82
enddo
endif
if(xyzzyaaez1)then
xyzzyaabn82=rele(1:3,xyzzyaaaa82)
call xyzzyaaqe1(xyzzyaabn82,xyzzyaaab82,.true.,.false.,.false.,xyzzyaa&
&au82,xyzzyaabj82,xyzzyaabb82)
xyzzyaabx82=xyzzyaabx82+xyzzyaaau82
endif
if(xyzzyaaff1)then
do xyzzyaaam82=1,xyzzyaabm1
xyzzyaaae82=xyzzyaabu1(xyzzyaaam82,1)
xyzzyaabo82=eivecs(:,xyzzyaaae82,xyzzyaaaa82)
if(xyzzyaabo82(4)>xyzzyaacp1)cycle
xyzzyaaal82=xyzzyaabu1(xyzzyaaam82,2)
do xyzzyaaac82=1,netot
if(xyzzyaaaa82==xyzzyaaac82)cycle
xyzzyaabp82=eivecs(:,xyzzyaaal82,xyzzyaaac82)
if(xyzzyaabp82(4)<xyzzyaacp1) then
call xyzzyaaqq1(xyzzyaabo82,xyzzyaabp82,xyzzyaaam82,1,.true.,.false.,.&
&false.,xyzzyaaav82,xyzzyaabj82,xyzzyaabb82)
xyzzyaabz82=xyzzyaabz82+xyzzyaaav82
endif
enddo
enddo
endif
enddo
if(xyzzyaafe1)xyzzyaabt82=xyzzyaabt82*third
if(xyzzyaafa1)then
call xyzzyaaqf1(.true.,.false.,.false.,1,eevecs(:,:,1),xyzzyaaby82,xyz&
&zyaabj82,xyzzyaabb82)
elseif(xyzzyaafb1)then
call xyzzyaaqg1(.true.,.false.,.false.,1,eevecs,xyzzyaaby82,xyzzyaabj8&
&2,xyzzyaabb82)
elseif(xyzzyaafc1)then
if(fix_holes)then
call xyzzyaaql1(.true.,.false.,.false.,1,rele,xyzzyaaby82,xyzzyaabj82,&
&xyzzyaabb82)
else
call xyzzyaaqi1(.true.,.false.,.false.,1,rele,xyzzyaaby82,xyzzyaabj82,&
&xyzzyaabb82)
endif
elseif(xyzzyaafg1)then
call xyzzyaaqp1(.true.,.false.,.false.,1,eevecs(:,:,1),eevecs,eivecs(:&
&,:,1),eivecs,xyzzyaaby82,xyzzyaabj82,xyzzyaabb82)
endif
value_jas=xyzzyaabq82+xyzzyaabs82+xyzzyaabt82+xyzzyaabu82+xyzzyaabv82+&
&xyzzyaabw82+xyzzyaabx82+xyzzyaaby82+xyzzyaabz82+xyzzyaabr82
call timer('JASTROW',.false.)
if(xyzzyaaaa1)call wout('full_jastrow_val: value_jas=',value_jas)
end subroutine xyzzyaapn1
subroutine xyzzyaapo1(eevecs,fd,sd)
use slaarnabt, only : dcopy
implicit none
real(dp),intent(in) :: eevecs(4,netot,netot)
logical,intent(in) :: fd,sd
call timer('JASTROW',.true.)
if(xyzzyaaeu1)call xyzzyaapw1(0,eevecs,fd.or.sd,sd)
if(xyzzyaaev1.or.xyzzyaafd1.or.xyzzyaafg1)call dcopy(four_netot_netot,&
&eevecs(1,1,1),1,xyzzyaaig1(1,1,1),1)
call timer('JASTROW',.false.)
end subroutine xyzzyaapo1
subroutine xyzzyaapp1(i,ispin,rvec,eevecs,eevecs1,eivecs,eivecs1,val,f&
&d,sd,val_x_q,fd_x_q,value_jas,grad_jas,lap_jas,value_x_q,grad_x_q,val&
&ue_3body,precomp_3body,total_b)
use slaarnabt, only : ddot,dcopy,dsum,dsum3
implicit none
integer,intent(in) :: i,ispin
real(dp),intent(in) :: rvec(3),eevecs(4,netot,netot),eevecs1(4,netot),&
&eivecs(4,nitot,netot),eivecs1(4,nitot)
real(dp),intent(out) :: value_jas,lap_jas,value_3body,total_b
real(dp),intent(inout) :: value_x_q,grad_jas(3),grad_x_q(3)
logical,intent(in) :: val,fd,sd,val_x_q,fd_x_q,precomp_3body
integer xyzzyaaaa84,xyzzyaaab84,xyzzyaaac84,xyzzyaaad84,xyzzyaaae84,xy&
&zzyaaaf84,xyzzyaaag84,xyzzyaaah84,xyzzyaaai84,xyzzyaaaj84,xyzzyaaak84&
&,xyzzyaaal84,xyzzyaaam84,xyzzyaaan84,s,xyzzyaaao84(nspin),xyzzyaaap84
real(dp) xyzzyaaaq84,xyzzyaaar84(3),xyzzyaaas84,xyzzyaaat84,xyzzyaaau8&
&4(3),xyzzyaaav84,xyzzyaaaw84,xyzzyaaax84(3),xyzzyaaay84,xyzzyaaaz84,x&
&yzzyaaba84(3),xyzzyaabb84,xyzzyaabc84,xyzzyaabd84(3),xyzzyaabe84,xyzz&
&yaabf84,xyzzyaabg84(3),xyzzyaabh84,xyzzyaabi84,xyzzyaabj84(3),xyzzyaa&
&bk84,xyzzyaabl84,xyzzyaabm84(3),xyzzyaabn84,xyzzyaabo84(3),xyzzyaabp8&
&4,xyzzyaabq84,xyzzyaabr84(3),xyzzyaabs84,xyzzyaabt84,xyzzyaabu84,xyzz&
&yaabv84,xyzzyaabw84,xyzzyaabx84,xyzzyaaby84,xyzzyaabz84,xyzzyaaca84,x&
&yzzyaacb84,xyzzyaacc84,xyzzyaacd84,xyzzyaace84,xyzzyaacf84,xyzzyaacg8&
&4,xyzzyaach84,xyzzyaaci84,xyzzyaacj84(3),xyzzyaack84,xyzzyaacl84,xyzz&
&yaacm84(3),xyzzyaacn84,xyzzyaaco84,xyzzyaacp84,xyzzyaacq84,xyzzyaacr8&
&4,xyzzyaacs84,xyzzyaact84,xyzzyaacu84,xyzzyaacv84,xyzzyaacw84,xyzzyaa&
&cx84,xyzzyaacy84,xyzzyaacz84,xyzzyaada84(3),xyzzyaadb84(3),xyzzyaadc8&
&4(3),xyzzyaadd84,xyzzyaade84,xyzzyaadf84,xyzzyaadg84(4),xyzzyaadh84(4&
&),xyzzyaadi84(3),xyzzyaadj84,xyzzyaadk84,xyzzyaadl84,xyzzyaadm84,xyzz&
&yaadn84,xyzzyaado84,xyzzyaadp84,xyzzyaadq84,xyzzyaadr84,xyzzyaads84,x&
&yzzyaadt84,xyzzyaadu84,xyzzyaadv84,xyzzyaadw84,xyzzyaadx84,xyzzyaady8&
&4,xyzzyaadz84,xyzzyaaea84,xyzzyaaeb84,xyzzyaaec84,xyzzyaaed84
real(dp),allocatable,save :: xyzzyaaee84(:,:)
logical xyzzyaaef84,xyzzyaaeg84,xyzzyaaeh84,xyzzyaaei84,xyzzyaaej84(ma&
&x_spin_pairs)
logical,save :: xyzzyaaek84=.true.
call timer('JASTROW',.true.)
if(xyzzyaaek84)then
if(xyzzyaafc1)then
allocate(xyzzyaaee84(3,netot),stat=xyzzyaaai84)
call check_alloc(xyzzyaaai84,'ONEELEC_JASTROW','')
endif
xyzzyaaek84=.false.
endif
xyzzyaaaq84=0.d0
xyzzyaaat84=0.d0
xyzzyaaaz84=0.d0
xyzzyaabc84=0.d0
xyzzyaaaw84=0.d0
xyzzyaabf84=0.d0
xyzzyaabi84=0.d0
xyzzyaabl84=0.d0
total_b=0.d0
xyzzyaaar84=0.d0
xyzzyaaau84=0.d0
xyzzyaaba84=0.d0
xyzzyaabd84=0.d0
xyzzyaaax84=0.d0
xyzzyaabg84=0.d0
xyzzyaabj84=0.d0
xyzzyaabm84=0.d0
xyzzyaabo84=0.d0
xyzzyaaas84=0.d0
xyzzyaaav84=0.d0
xyzzyaabb84=0.d0
xyzzyaabe84=0.d0
xyzzyaaay84=0.d0
xyzzyaabh84=0.d0
xyzzyaabk84=0.d0
xyzzyaabn84=0.d0
xyzzyaabp84=0.d0
xyzzyaabq84=0.d0
xyzzyaabr84=0.d0
xyzzyaabs84=0.d0
xyzzyaaef84=fd.or.sd
xyzzyaaeh84=fd_x_q.or.sd
xyzzyaaeg84=val_x_q.or.xyzzyaaeh84
if(dimensionality==3)then
!$omp parallel default(none) shared(i,netot,xyzzyaaep1,xyzzyaaeu1,xyzz&
!$omp &yaaex1,which_spin,val,xyzzyaaef84,fd,sd,eevecs1,ispin,xyzzyaaaq&
!$omp &84,xyzzyaaar84,xyzzyaaas84,xyzzyaabo1,eivecs1,xyzzyaaew1,xyzzya&
!$omp &abc84,xyzzyaabd84,xyzzyaabe84,nitot,xyzzyaaeg84,val_x_q,xyzzyaa&
!$omp &eh84,fd_x_q,xyzzyaabp1,xyzzyaacv1,xyzzyaacg1,eivecs,xyzzyaabf84&
!$omp &,xyzzyaabg84,xyzzyaabh84,nspin,xyzzyaady1,max_spin_pairs,xyzzya&
!$omp &acf1,hard_diam,xyzzyaaoa1,xyzzyaaob1,xyzzyaanz1,hard_sphere,har&
!$omp &d_op_spins,xyzzyaabx1,xyzzyaagx1,xyzzyaaaj1,which_spair,xyzzyaa&
!$omp &eq1,levels_spairs,xyzzyaafh1,xyzzyaaam1,xyzzyaach1,xyzzyaaes1,x&
!$omp &yzzyaacj1,xyzzyaacr1,xyzzyaaat84,xyzzyaaau84,xyzzyaaav84) priva&
!$omp &te(xyzzyaaaa84,xyzzyaaab84,xyzzyaacx84,xyzzyaabt84,xyzzyaabu84,&
!$omp &xyzzyaabv84,xyzzyaadl84,xyzzyaadj84,xyzzyaaah84,xyzzyaadi84,xyz&
!$omp &zyaadk84,xyzzyaaac84,xyzzyaaag84,xyzzyaacy84,xyzzyaace84,xyzzya&
!$omp &acf84,xyzzyaacg84,xyzzyaaad84,xyzzyaadf84,xyzzyaadc84,xyzzyaacz&
!$omp &84,xyzzyaadd84,xyzzyaaae84,xyzzyaaaf84,xyzzyaach84,xyzzyaaco84,&
!$omp &xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
!$omp &zyaacs84,xyzzyaacv84,xyzzyaade84,xyzzyaacw84,xyzzyaadb84,xyzzya&
!$omp &aej84,xyzzyaaao84,xyzzyaadn84,xyzzyaaal84,xyzzyaaam84,xyzzyaads&
!$omp &84,xyzzyaadt84,xyzzyaadv84,xyzzyaadw84,xyzzyaady84,xyzzyaadx84,&
!$omp &xyzzyaaei84,xyzzyaaan84,xyzzyaadq84,xyzzyaadp84,xyzzyaadr84,xyz&
!$omp &zyaado84,xyzzyaadu84,s,xyzzyaaeb84,xyzzyaaea84,xyzzyaaap84,xyzz&
!$omp &yaadz84,xyzzyaaec84,xyzzyaaed84,xyzzyaabw84,xyzzyaabx84,xyzzyaa&
!$omp &by84,xyzzyaabz84,xyzzyaaca84,xyzzyaadm84)
xyzzyaaej84(:)=.true.
do xyzzyaaab84=1,nspin
xyzzyaaao84(xyzzyaaab84)=which_spair(ispin,xyzzyaaab84,levels_spairs)
xyzzyaaej84(xyzzyaaao84(xyzzyaaab84))=.false.
enddo
if(xyzzyaaep1)then
!$omp  master
call timer('U TERM',.true.)
!$omp  end master
if(xyzzyaaef84)then
!$omp  do reduction(+:xyzzyaaaq84,xyzzyaaar84,xyzzyaaas84)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
call xyzzyaapr1(xyzzyaacx84,ispin,xyzzyaaab84,val,xyzzyaaef84,sd,xyzzy&
&aabt84,xyzzyaabu84,xyzzyaabv84)
if(val)xyzzyaaaq84=xyzzyaaaq84+xyzzyaabt84
xyzzyaadl84=xyzzyaabu84/xyzzyaacx84
if(fd)xyzzyaaar84(1:3)=xyzzyaaar84(1:3)+xyzzyaadl84*eevecs1(1:3,xyzzya&
&aaa84)
if(sd)xyzzyaaas84=xyzzyaaas84+(xyzzyaabv84+2*xyzzyaadl84)
enddo
!$omp  enddo nowait
elseif(fd)then
!$omp  do reduction(+:xyzzyaaaq84,xyzzyaaar84)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
call xyzzyaapr1(xyzzyaacx84,ispin,xyzzyaaab84,val,xyzzyaaef84,sd,xyzzy&
&aabt84,xyzzyaabu84,xyzzyaabv84)
if(val)xyzzyaaaq84=xyzzyaaaq84+xyzzyaabt84
xyzzyaadl84=xyzzyaabu84/xyzzyaacx84
xyzzyaaar84(1:3)=xyzzyaaar84(1:3)+xyzzyaadl84*eevecs1(1:3,xyzzyaaaa84)
enddo
!$omp  enddo nowait
else
if(xyzzyaaeq1)then
!$omp  do reduction(+:xyzzyaaaq84)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
if(xyzzyaacx84>xyzzyaach1.and.xyzzyaagx1>0)then
xyzzyaabt84=0.d0
else
s=xyzzyaaao84(xyzzyaaab84)
xyzzyaadv84=xyzzyaacx84-xyzzyaach1
if(xyzzyaagx1==2)then
xyzzyaady84=xyzzyaaaj1(0,s)+xyzzyaacx84*xyzzyaaaj1(1,s)
xyzzyaadx84=xyzzyaacx84
do xyzzyaaad84=2,xyzzyaabx1
xyzzyaadx84=xyzzyaadx84*xyzzyaacx84
xyzzyaady84=xyzzyaady84+xyzzyaaaj1(xyzzyaaad84,s)*xyzzyaadx84
enddo
xyzzyaabt84=xyzzyaady84*xyzzyaadv84*xyzzyaadv84
elseif(xyzzyaagx1==3)then
xyzzyaady84=xyzzyaaaj1(0,s)+xyzzyaacx84*xyzzyaaaj1(1,s)
xyzzyaadx84=xyzzyaacx84
do xyzzyaaad84=2,xyzzyaabx1
xyzzyaadx84=xyzzyaadx84*xyzzyaacx84
xyzzyaady84=xyzzyaady84+xyzzyaaaj1(xyzzyaaad84,s)*xyzzyaadx84
enddo
xyzzyaadw84=xyzzyaadv84*xyzzyaadv84
xyzzyaabt84=xyzzyaady84*xyzzyaadw84*xyzzyaadv84
else
xyzzyaady84=xyzzyaaaj1(0,s)+xyzzyaacx84*xyzzyaaaj1(1,s)
xyzzyaadx84=xyzzyaacx84
do xyzzyaaad84=2,xyzzyaabx1
xyzzyaadx84=xyzzyaadx84*xyzzyaacx84
xyzzyaady84=xyzzyaady84+xyzzyaaaj1(xyzzyaaad84,s)*xyzzyaadx84
enddo
xyzzyaabt84=xyzzyaady84*xyzzyaadv84**xyzzyaagx1
endif
endif
if(hard_sphere)then
if(.not.hard_op_spins.or.ispin/=xyzzyaaab84)then
if(xyzzyaacx84<xyzzyaanz1)then
if(xyzzyaacx84<=hard_diam)call errstop('COMPUTE_U','Hard spheres are o&
&verlapping.  This is not supposed to happen.')
xyzzyaads84=xyzzyaacx84*xyzzyaaoa1
xyzzyaadt84=xyzzyaacx84*xyzzyaaob1
xyzzyaadu84=tanh((xyzzyaads84-1.d0)/(1.d0-xyzzyaadt84))
xyzzyaabt84=xyzzyaabt84+log(xyzzyaadu84)
endif
endif
endif
xyzzyaaaq84=xyzzyaaaq84+xyzzyaabt84
enddo
!$omp  enddo nowait
else
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
if(xyzzyaacx84>xyzzyaach1.and.xyzzyaagx1>0)then
xyzzyaabt84=0.d0
else
s=which_spair(ispin,xyzzyaaab84,xyzzyaafh1)
call xyzzyaapt1(xyzzyaacx84,s,val,fd,sd,xyzzyaabt84,xyzzyaabu84,xyzzya&
&abv84)
endif
if(hard_sphere)then
if(.not.hard_op_spins.or.ispin/=xyzzyaaab84)then
if(xyzzyaacx84<xyzzyaanz1)then
if(xyzzyaacx84<=hard_diam)call errstop('COMPUTE_U','Hard spheres are o&
&verlapping.  This is not supposed to happen.')
xyzzyaads84=xyzzyaacx84*xyzzyaaoa1
xyzzyaadt84=xyzzyaacx84*xyzzyaaob1
xyzzyaadu84=tanh((xyzzyaads84-1.d0)/(1.d0-xyzzyaadt84))
xyzzyaabt84=xyzzyaabt84+log(xyzzyaadu84)
endif
endif
endif
xyzzyaaaq84=xyzzyaaaq84+xyzzyaabt84
enddo
endif
endif
!$omp  master
call timer('U TERM',.false.)
!$omp  end master
endif
if(xyzzyaaes1)then
call timer('Ucyl TERM',.true.)
if(xyzzyaaef84)then
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
if(eevecs1(4,xyzzyaaaa84)<=0.d0)cycle
xyzzyaaed84=dot_product(eevecs1(1:3,xyzzyaaaa84),xyzzyaacr1)
xyzzyaaec84=eevecs1(4,xyzzyaaaa84)*eevecs1(4,xyzzyaaaa84)-xyzzyaaed84*&
&xyzzyaaed84
if(xyzzyaaec84<xyzzyaacj1.or.xyzzyaagx1==0)then
xyzzyaaec84=sqrt(max(xyzzyaaec84,0.d0))
call xyzzyaapu1(xyzzyaaec84,abs(xyzzyaaed84),ispin,xyzzyaaab84,val,xyz&
&zyaaef84,sd,xyzzyaabw84,xyzzyaabx84,xyzzyaaby84,xyzzyaabz84,xyzzyaaca&
&84)
if(val)xyzzyaaat84=xyzzyaaat84+xyzzyaabw84
if(xyzzyaaec84>0.d0)then
xyzzyaadm84=xyzzyaabx84/xyzzyaaec84
else
xyzzyaadm84=0.d0
endif
if(fd)xyzzyaaau84=xyzzyaaau84+xyzzyaadm84*(eevecs1(1:3,xyzzyaaaa84)-xy&
&zzyaaed84*xyzzyaacr1)+xyzzyaaby84*sign(1.d0,xyzzyaaed84)*xyzzyaacr1
xyzzyaaav84=xyzzyaaav84+xyzzyaabz84+xyzzyaadm84+xyzzyaaca84
endif
enddo
elseif(fd)then
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
if(eevecs1(4,xyzzyaaaa84)<=0.d0)cycle
xyzzyaaed84=dot_product(eevecs1(1:3,xyzzyaaaa84),xyzzyaacr1)
xyzzyaaec84=eevecs1(4,xyzzyaaaa84)*eevecs1(4,xyzzyaaaa84)-xyzzyaaed84*&
&xyzzyaaed84
if(xyzzyaaec84<xyzzyaacj1.or.xyzzyaagx1==0)then
xyzzyaaec84=sqrt(xyzzyaaec84)
call xyzzyaapu1(xyzzyaaec84,abs(xyzzyaaed84),ispin,xyzzyaaab84,val,xyz&
&zyaaef84,sd,xyzzyaabw84,xyzzyaabx84,xyzzyaaby84,xyzzyaabz84,xyzzyaaca&
&84)
if(val)xyzzyaaat84=xyzzyaaat84+xyzzyaabw84
if(xyzzyaaec84>0.d0)then
xyzzyaadm84=xyzzyaabx84/xyzzyaaec84
else
xyzzyaadm84=0.d0
endif
xyzzyaaau84=xyzzyaaau84+xyzzyaadm84*(eevecs1(1:3,xyzzyaaaa84)-xyzzyaae&
&d84*xyzzyaacr1)+xyzzyaaby84*sign(1.d0,xyzzyaaed84)*xyzzyaacr1
endif
enddo
else
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
if(eevecs1(4,xyzzyaaaa84)<=0.d0)cycle
xyzzyaaed84=dot_product(eevecs1(1:3,xyzzyaaaa84),xyzzyaacr1)
xyzzyaaec84=eevecs1(4,xyzzyaaaa84)*eevecs1(4,xyzzyaaaa84)-xyzzyaaed84*&
&xyzzyaaed84
if(xyzzyaaec84<xyzzyaacj1.or.xyzzyaagx1==0)then
xyzzyaaec84=sqrt(xyzzyaaec84)
call xyzzyaapu1(xyzzyaaec84,abs(xyzzyaaed84),ispin,xyzzyaaab84,val,xyz&
&zyaaef84,sd,xyzzyaabw84,xyzzyaabx84,xyzzyaaby84,xyzzyaabz84,xyzzyaaca&
&84)
xyzzyaaat84=xyzzyaaat84+xyzzyaabw84
endif
enddo
endif
call timer('Ucyl TERM',.false.)
endif
if(xyzzyaaew1.and.xyzzyaaeg84)then
!$omp  master
call timer('CHI TERM',.true.)
!$omp  end master
!$omp  do reduction(+:xyzzyaabc84,xyzzyaabd84,xyzzyaabe84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabo1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacy84<=0.d0)cycle
call xyzzyaaqa1(xyzzyaacy84,ispin,xyzzyaaag84,val_x_q,xyzzyaaeh84,sd,x&
&yzzyaace84,xyzzyaacf84,xyzzyaacg84)
if(val_x_q)xyzzyaabc84=xyzzyaabc84+xyzzyaace84
if(xyzzyaaeh84)then
xyzzyaadl84=xyzzyaacf84/xyzzyaacy84
if(fd_x_q)xyzzyaabd84(1:3)=xyzzyaabd84(1:3)+xyzzyaadl84*eivecs1(1:3,xy&
&zzyaaac84)
if(sd)xyzzyaabe84=xyzzyaabe84+(xyzzyaacg84+2*xyzzyaadl84)
endif
enddo
!$omp  enddo nowait
!$omp  master
call timer('CHI TERM',.false.)
!$omp  end master
endif
if(xyzzyaaex1)then
!$omp  master
call timer('F TERM',.true.)
!$omp  end master
select case(xyzzyaady1)
case(0)
if(xyzzyaaef84)then
!$omp  do reduction(+:xyzzyaabf84,xyzzyaabg84,xyzzyaabh84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
xyzzyaaaf84=xyzzyaacg1(xyzzyaaag84)
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaaaa84==i)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84.or.xyzzyaacz84<=0.d0)cycle
xyzzyaadd84=xyzzyaacz84
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaadd84
do xyzzyaaae84=2,xyzzyaaaf84
xyzzyaadd84=xyzzyaadd84*xyzzyaacz84
xyzzyaadq1(xyzzyaaae84,xyzzyaaaa84)=xyzzyaadd84
enddo
enddo
xyzzyaado1(1:xyzzyaaaf84)=xyzzyaadq1(1:xyzzyaaaf84,i)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
call xyzzyaaqb1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
if(fd)xyzzyaabg84(1:3)=xyzzyaabg84(1:3)+(xyzzyaacp84*xyzzyaadb84(1:3)+&
&xyzzyaaco84*xyzzyaadc84(1:3))
if(sd)then
xyzzyaacw84=xyzzyaadc84(1)*xyzzyaadb84(1)+xyzzyaadc84(2)*xyzzyaadb84(2&
&)+xyzzyaadc84(3)*xyzzyaadb84(3)
xyzzyaabh84=xyzzyaabh84+(2.d0*xyzzyaaco84*xyzzyaadf84+xyzzyaacq84+2.d0&
&*xyzzyaacp84*xyzzyaade84+xyzzyaacr84+2.d0*xyzzyaacw84*xyzzyaacs84)
endif
enddo
enddo
!$omp  enddo nowait
elseif(fd)then
!$omp  do reduction(+:xyzzyaabf84,xyzzyaabg84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
xyzzyaaaf84=xyzzyaacg1(xyzzyaaag84)
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaaaa84==i)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84.or.xyzzyaacz84<=0.d0)cycle
xyzzyaadd84=xyzzyaacz84
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaadd84
do xyzzyaaae84=2,xyzzyaaaf84
xyzzyaadd84=xyzzyaadd84*xyzzyaacz84
xyzzyaadq1(xyzzyaaae84,xyzzyaaaa84)=xyzzyaadd84
enddo
enddo
xyzzyaado1(1:xyzzyaaaf84)=xyzzyaadq1(1:xyzzyaaaf84,i)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
call xyzzyaaqb1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
xyzzyaabg84(1:3)=xyzzyaabg84(1:3)+(xyzzyaacp84*xyzzyaadb84(1:3)+xyzzya&
&aco84*xyzzyaadc84(1:3))
enddo
enddo
!$omp  enddo nowait
else
!$omp  do reduction(+:xyzzyaabf84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
xyzzyaaei84=xyzzyaacy84>xyzzyaadn84
if(xyzzyaaei84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaado84=xyzzyaacy84-xyzzyaadn84
xyzzyaaal84=xyzzyaacf1(xyzzyaaag84)
xyzzyaaam84=xyzzyaacg1(xyzzyaaag84)
do xyzzyaaae84=1,xyzzyaaam84
xyzzyaado1(xyzzyaaae84)=xyzzyaado1(xyzzyaaae84-1)*xyzzyaacy84
enddo
do s=1,max_spin_pairs
if(xyzzyaaej84(s))cycle
xyzzyaaan84=0
do xyzzyaaaf84=0,xyzzyaaal84
do xyzzyaaae84=0,xyzzyaaam84
xyzzyaadq84=0.d0
do xyzzyaaad84=0,xyzzyaaam84
xyzzyaaan84=xyzzyaaan84+1
xyzzyaadq84=xyzzyaadq84+xyzzyaaam1(xyzzyaaan84,s,xyzzyaaag84)*xyzzyaad&
&o1(xyzzyaaad84)
enddo
xyzzyaadx1(xyzzyaaae84,xyzzyaaaf84,s)=xyzzyaadq84
enddo
enddo
enddo
select case(xyzzyaagx1)
case(3)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaacz84<xyzzyaadn84)then
xyzzyaadr84=(xyzzyaado84*(xyzzyaacz84-xyzzyaadn84))**3
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
do xyzzyaaae84=1,xyzzyaaal84
xyzzyaadr1(xyzzyaaae84)=xyzzyaadr1(xyzzyaaae84-1)*xyzzyaacx84
enddo
do xyzzyaaae84=1,xyzzyaaam84
xyzzyaadp1(xyzzyaaae84)=xyzzyaadp1(xyzzyaaae84-1)*xyzzyaacz84
enddo
xyzzyaach84=0.d0
do xyzzyaaaf84=0,xyzzyaaal84
xyzzyaadp84=sum(xyzzyaadx1(0:xyzzyaaam84,xyzzyaaaf84,s)*xyzzyaadp1(0:x&
&yzzyaaam84))
xyzzyaach84=xyzzyaach84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
xyzzyaabf84=xyzzyaabf84+xyzzyaach84*xyzzyaadr84
endif
enddo
case(0)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
do xyzzyaaae84=1,xyzzyaaal84
xyzzyaadr1(xyzzyaaae84)=xyzzyaadr1(xyzzyaaae84-1)*xyzzyaacx84
enddo
do xyzzyaaae84=1,xyzzyaaam84
xyzzyaadp1(xyzzyaaae84)=xyzzyaadp1(xyzzyaaae84-1)*xyzzyaacz84
enddo
do xyzzyaaaf84=0,xyzzyaaal84
xyzzyaadp84=sum(xyzzyaadx1(0:xyzzyaaam84,xyzzyaaaf84,s)*xyzzyaadp1(0:x&
&yzzyaaam84))
xyzzyaabf84=xyzzyaabf84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
enddo
case default
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaacz84<xyzzyaadn84)then
xyzzyaadr84=(xyzzyaado84*(xyzzyaacz84-xyzzyaadn84))**xyzzyaagx1
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
do xyzzyaaae84=1,xyzzyaaal84
xyzzyaadr1(xyzzyaaae84)=xyzzyaadr1(xyzzyaaae84-1)*xyzzyaacx84
enddo
do xyzzyaaae84=1,xyzzyaaam84
xyzzyaadp1(xyzzyaaae84)=xyzzyaadp1(xyzzyaaae84-1)*xyzzyaacz84
enddo
xyzzyaach84=0.d0
do xyzzyaaaf84=0,xyzzyaaal84
xyzzyaadp84=sum(xyzzyaadx1(0:xyzzyaaam84,xyzzyaaaf84,s)*xyzzyaadp1(0:x&
&yzzyaaam84))
xyzzyaach84=xyzzyaach84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
xyzzyaabf84=xyzzyaabf84+xyzzyaach84*xyzzyaadr84
endif
enddo
end select
enddo
!$omp  enddo nowait
endif
case(1)
if(xyzzyaaef84)then
!$omp  do reduction(+:xyzzyaabf84,xyzzyaabg84,xyzzyaabh84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaaaa84==i)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84.or.xyzzyaacz84<=0.d0)cycle
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaacz84
enddo
xyzzyaado1(1)=xyzzyaadq1(1,i)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
call xyzzyaaqb1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
if(fd)xyzzyaabg84(1:3)=xyzzyaabg84(1:3)+(xyzzyaacp84*xyzzyaadb84(1:3)+&
&xyzzyaaco84*xyzzyaadc84(1:3))
if(sd)then
xyzzyaacw84=xyzzyaadc84(1)*xyzzyaadb84(1)+xyzzyaadc84(2)*xyzzyaadb84(2&
&)+xyzzyaadc84(3)*xyzzyaadb84(3)
xyzzyaabh84=xyzzyaabh84+(2.d0*xyzzyaaco84*xyzzyaadf84+xyzzyaacq84+2.d0&
&*xyzzyaacp84*xyzzyaade84+xyzzyaacr84+2.d0*xyzzyaacw84*xyzzyaacs84)
endif
enddo
enddo
!$omp  enddo nowait
elseif(fd)then
!$omp  do reduction(+:xyzzyaabf84,xyzzyaabg84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaaaa84==i)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84.or.xyzzyaacz84<=0.d0)cycle
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaacz84
enddo
xyzzyaado1(1)=xyzzyaadq1(1,i)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
call xyzzyaaqb1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
xyzzyaabg84(1:3)=xyzzyaabg84(1:3)+(xyzzyaacp84*xyzzyaadb84(1:3)+xyzzya&
&aco84*xyzzyaadc84(1:3))
enddo
enddo
!$omp  enddo nowait
else
!$omp  do reduction(+:xyzzyaabf84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
xyzzyaaei84=xyzzyaacy84>xyzzyaadn84
if(xyzzyaaei84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaado84=xyzzyaacy84-xyzzyaadn84
xyzzyaado1(1)=xyzzyaado1(0)*xyzzyaacy84
do s=1,max_spin_pairs
if(xyzzyaaej84(s))cycle
xyzzyaaan84=0
do xyzzyaaaf84=0,1
do xyzzyaaae84=0,1
xyzzyaadq84=0.d0
do xyzzyaaad84=0,1
xyzzyaaan84=xyzzyaaan84+1
xyzzyaadq84=xyzzyaadq84+xyzzyaaam1(xyzzyaaan84,s,xyzzyaaag84)*xyzzyaad&
&o1(xyzzyaaad84)
enddo
xyzzyaadx1(xyzzyaaae84,xyzzyaaaf84,s)=xyzzyaadq84
enddo
enddo
enddo
select case(xyzzyaagx1)
case(3)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaacz84<xyzzyaadn84)then
xyzzyaadr84=(xyzzyaado84*(xyzzyaacz84-xyzzyaadn84))**3
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaadr1(1)=xyzzyaadr1(0)*xyzzyaacx84
xyzzyaadp1(1)=xyzzyaadp1(0)*xyzzyaacz84
xyzzyaach84=0.d0
do xyzzyaaaf84=0,1
xyzzyaadp84=xyzzyaadx1(0,xyzzyaaaf84,s)*xyzzyaadp1(0)+xyzzyaadx1(1,xyz&
&zyaaaf84,s)*xyzzyaadp1(1)
xyzzyaach84=xyzzyaach84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
xyzzyaabf84=xyzzyaabf84+xyzzyaach84*xyzzyaadr84
endif
enddo
case(0)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaadr1(1)=xyzzyaadr1(0)*xyzzyaacx84
xyzzyaadp1(1)=xyzzyaadp1(0)*xyzzyaacz84
do xyzzyaaaf84=0,1
xyzzyaadp84=xyzzyaadx1(0,xyzzyaaaf84,s)*xyzzyaadp1(0)+xyzzyaadx1(1,xyz&
&zyaaaf84,s)*xyzzyaadp1(1)
xyzzyaabf84=xyzzyaabf84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
enddo
case default
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaacz84<xyzzyaadn84)then
xyzzyaadr84=(xyzzyaado84*(xyzzyaacz84-xyzzyaadn84))**xyzzyaagx1
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaadr1(1)=xyzzyaadr1(0)*xyzzyaacx84
xyzzyaadp1(1)=xyzzyaadp1(0)*xyzzyaacz84
xyzzyaach84=0.d0
do xyzzyaaaf84=0,1
xyzzyaadp84=xyzzyaadx1(0,xyzzyaaaf84,s)*xyzzyaadp1(0)+xyzzyaadx1(1,xyz&
&zyaaaf84,s)*xyzzyaadp1(1)
xyzzyaach84=xyzzyaach84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
xyzzyaabf84=xyzzyaabf84+xyzzyaach84*xyzzyaadr84
endif
enddo
end select
enddo
!$omp  enddo nowait
endif
case(2)
if(xyzzyaaef84)then
!$omp  do reduction(+:xyzzyaabf84,xyzzyaabg84,xyzzyaabh84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaaaa84==i)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>xyzzyaadn84.or.xyzzyaacz84<=0.d0)cycle
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaacz84
xyzzyaadq1(2,xyzzyaaaa84)=xyzzyaacz84*xyzzyaacz84
enddo
xyzzyaado1(1)=xyzzyaadq1(1,i)
xyzzyaado1(2)=xyzzyaadq1(2,i)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
call xyzzyaaqc1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
if(fd)xyzzyaabg84(1:3)=xyzzyaabg84(1:3)+(xyzzyaacp84*xyzzyaadb84(1:3)+&
&xyzzyaaco84*xyzzyaadc84(1:3))
if(sd)then
xyzzyaacw84=xyzzyaadc84(1)*xyzzyaadb84(1)+xyzzyaadc84(2)*xyzzyaadb84(2&
&)+xyzzyaadc84(3)*xyzzyaadb84(3)
xyzzyaabh84=xyzzyaabh84+(2.d0*xyzzyaaco84*xyzzyaadf84+xyzzyaacq84+2.d0&
&*xyzzyaacp84*xyzzyaade84+xyzzyaacr84+2.d0*xyzzyaacw84*xyzzyaacs84)
endif
enddo
enddo
!$omp  enddo nowait
elseif(fd)then
!$omp  do reduction(+:xyzzyaabf84,xyzzyaabg84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaaaa84==i)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84.or.xyzzyaacz84<=0.d0)cycle
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaacz84
xyzzyaadq1(2,xyzzyaaaa84)=xyzzyaacz84*xyzzyaacz84
enddo
xyzzyaado1(1:2)=xyzzyaadq1(1:2,i)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
call xyzzyaaqc1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
xyzzyaabg84(1:3)=xyzzyaabg84(1:3)+(xyzzyaacp84*xyzzyaadb84(1:3)+xyzzya&
&aco84*xyzzyaadc84(1:3))
enddo
enddo
!$omp  enddo nowait
else
!$omp  do reduction(+:xyzzyaabf84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
xyzzyaaei84=xyzzyaacy84>xyzzyaadn84
if(xyzzyaaei84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaado84=xyzzyaacy84-xyzzyaadn84
xyzzyaadz84=xyzzyaacy84*xyzzyaacy84
do s=1,max_spin_pairs
if(xyzzyaaej84(s))cycle
xyzzyaaan84=1
do xyzzyaaaf84=0,2
do xyzzyaaae84=0,2
xyzzyaadx1(xyzzyaaae84,xyzzyaaaf84,s)=xyzzyaaam1(xyzzyaaan84,s,xyzzyaa&
&ag84)+xyzzyaaam1(xyzzyaaan84+1,s,xyzzyaaag84)*xyzzyaacy84+xyzzyaaam1(&
&xyzzyaaan84+2,s,xyzzyaaag84)*xyzzyaadz84
xyzzyaaan84=xyzzyaaan84+3
enddo
enddo
enddo
select case(xyzzyaagx1)
case(3)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaacz84<xyzzyaadn84)then
xyzzyaadr84=(xyzzyaado84*(xyzzyaacz84-xyzzyaadn84))**3
xyzzyaaap84=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaaeb84=xyzzyaacx84*xyzzyaacx84
xyzzyaaea84=xyzzyaacz84*xyzzyaacz84
xyzzyaach84=xyzzyaadx1(0,0,xyzzyaaap84)+xyzzyaadx1(1,0,xyzzyaaap84)*xy&
&zzyaacz84+xyzzyaadx1(2,0,xyzzyaaap84)*xyzzyaaea84
xyzzyaach84=xyzzyaach84+(xyzzyaadx1(0,1,xyzzyaaap84)+xyzzyaadx1(1,1,xy&
&zzyaaap84)*xyzzyaacz84+xyzzyaadx1(2,1,xyzzyaaap84)*xyzzyaaea84)*xyzzy&
&aacx84
xyzzyaach84=xyzzyaach84+(xyzzyaadx1(0,2,xyzzyaaap84)+xyzzyaadx1(1,2,xy&
&zzyaaap84)*xyzzyaacz84+xyzzyaadx1(2,2,xyzzyaaap84)*xyzzyaaea84)*xyzzy&
&aaeb84
xyzzyaabf84=xyzzyaabf84+xyzzyaach84*xyzzyaadr84
endif
enddo
case(0)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
xyzzyaaap84=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaaeb84=xyzzyaacx84*xyzzyaacx84
xyzzyaaea84=xyzzyaacz84*xyzzyaacz84
xyzzyaach84=xyzzyaadx1(0,0,xyzzyaaap84)+xyzzyaadx1(1,0,xyzzyaaap84)*xy&
&zzyaacz84+xyzzyaadx1(2,0,xyzzyaaap84)*xyzzyaaea84
xyzzyaach84=xyzzyaach84+(xyzzyaadx1(0,1,xyzzyaaap84)+xyzzyaadx1(1,1,xy&
&zzyaaap84)*xyzzyaacz84+xyzzyaadx1(2,1,xyzzyaaap84)*xyzzyaaea84)*xyzzy&
&aacx84
xyzzyaach84=xyzzyaach84+(xyzzyaadx1(0,2,xyzzyaaap84)+xyzzyaadx1(1,2,xy&
&zzyaaap84)*xyzzyaacz84+xyzzyaadx1(2,2,xyzzyaaap84)*xyzzyaaea84)*xyzzy&
&aaeb84
xyzzyaabf84=xyzzyaabf84+xyzzyaach84
enddo
case default
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaacz84<xyzzyaadn84)then
xyzzyaadr84=(xyzzyaado84*(xyzzyaacz84-xyzzyaadn84))**xyzzyaagx1
xyzzyaaap84=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaaeb84=xyzzyaacx84*xyzzyaacx84
xyzzyaaea84=xyzzyaacz84*xyzzyaacz84
xyzzyaach84=xyzzyaadx1(0,0,xyzzyaaap84)+xyzzyaadx1(1,0,xyzzyaaap84)*xy&
&zzyaacz84+xyzzyaadx1(2,0,xyzzyaaap84)*xyzzyaaea84
xyzzyaach84=xyzzyaach84+(xyzzyaadx1(0,1,xyzzyaaap84)+xyzzyaadx1(1,1,xy&
&zzyaaap84)*xyzzyaacz84+xyzzyaadx1(2,1,xyzzyaaap84)*xyzzyaaea84)*xyzzy&
&aacx84
xyzzyaach84=xyzzyaach84+(xyzzyaadx1(0,2,xyzzyaaap84)+xyzzyaadx1(1,2,xy&
&zzyaaap84)*xyzzyaacz84+xyzzyaadx1(2,2,xyzzyaaap84)*xyzzyaaea84)*xyzzy&
&aaeb84
xyzzyaabf84=xyzzyaabf84+xyzzyaach84*xyzzyaadr84
endif
enddo
end select
enddo
!$omp  enddo nowait
endif
case(3)
if(xyzzyaaef84)then
!$omp  do reduction(+:xyzzyaabf84,xyzzyaabg84,xyzzyaabh84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaaaa84==i)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84.or.xyzzyaacz84<=0.d0)cycle
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaacz84
xyzzyaadd84=xyzzyaacz84*xyzzyaacz84
xyzzyaadq1(2,xyzzyaaaa84)=xyzzyaadd84
xyzzyaadq1(3,xyzzyaaaa84)=xyzzyaadd84*xyzzyaacz84
enddo
xyzzyaado1(1:3)=xyzzyaadq1(1:3,i)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
call xyzzyaaqb1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
if(fd)xyzzyaabg84(1:3)=xyzzyaabg84(1:3)+(xyzzyaacp84*xyzzyaadb84(1:3)+&
&xyzzyaaco84*xyzzyaadc84(1:3))
if(sd)then
xyzzyaacw84=xyzzyaadc84(1)*xyzzyaadb84(1)+xyzzyaadc84(2)*xyzzyaadb84(2&
&)+xyzzyaadc84(3)*xyzzyaadb84(3)
xyzzyaabh84=xyzzyaabh84+(2.d0*xyzzyaaco84*xyzzyaadf84+xyzzyaacq84+2.d0&
&*xyzzyaacp84*xyzzyaade84+xyzzyaacr84+2.d0*xyzzyaacw84*xyzzyaacs84)
endif
enddo
enddo
!$omp  enddo nowait
elseif(fd)then
!$omp  do reduction(+:xyzzyaabf84,xyzzyaabg84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaaaa84==i)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84.or.xyzzyaacz84<=0.d0)cycle
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaacz84
xyzzyaadd84=xyzzyaacz84*xyzzyaacz84
xyzzyaadq1(2,xyzzyaaaa84)=xyzzyaadd84
xyzzyaadq1(3,xyzzyaaaa84)=xyzzyaadd84*xyzzyaacz84
enddo
xyzzyaado1(1:3)=xyzzyaadq1(1:3,i)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
call xyzzyaaqb1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
xyzzyaabg84(1:3)=xyzzyaabg84(1:3)+(xyzzyaacp84*xyzzyaadb84(1:3)+xyzzya&
&aco84*xyzzyaadc84(1:3))
enddo
enddo
!$omp  enddo nowait
else
!$omp  do reduction(+:xyzzyaabf84)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
xyzzyaaei84=xyzzyaacy84>xyzzyaadn84
if(xyzzyaaei84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaado84=xyzzyaacy84-xyzzyaadn84
xyzzyaado1(1)=xyzzyaado1(0)*xyzzyaacy84
xyzzyaado1(2)=xyzzyaado1(1)*xyzzyaacy84
xyzzyaado1(3)=xyzzyaado1(2)*xyzzyaacy84
do s=1,max_spin_pairs
if(xyzzyaaej84(s))cycle
xyzzyaaan84=0
do xyzzyaaaf84=0,3
do xyzzyaaae84=0,3
xyzzyaadq84=0.d0
do xyzzyaaad84=0,3
xyzzyaaan84=xyzzyaaan84+1
xyzzyaadq84=xyzzyaadq84+xyzzyaaam1(xyzzyaaan84,s,xyzzyaaag84)*xyzzyaad&
&o1(xyzzyaaad84)
enddo
xyzzyaadx1(xyzzyaaae84,xyzzyaaaf84,s)=xyzzyaadq84
enddo
enddo
enddo
select case(xyzzyaagx1)
case(3)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaacz84<xyzzyaadn84)then
xyzzyaadr84=(xyzzyaado84*(xyzzyaacz84-xyzzyaadn84))**3
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaadr1(1)=xyzzyaadr1(0)*xyzzyaacx84
xyzzyaadr1(2)=xyzzyaadr1(1)*xyzzyaacx84
xyzzyaadr1(3)=xyzzyaadr1(2)*xyzzyaacx84
xyzzyaadp1(1)=xyzzyaadp1(0)*xyzzyaacz84
xyzzyaadp1(2)=xyzzyaadp1(1)*xyzzyaacz84
xyzzyaadp1(3)=xyzzyaadp1(2)*xyzzyaacz84
xyzzyaach84=0.d0
do xyzzyaaaf84=0,3
xyzzyaadp84=sum(xyzzyaadx1(0:3,xyzzyaaaf84,s)*xyzzyaadp1(0:3))
xyzzyaach84=xyzzyaach84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
xyzzyaabf84=xyzzyaabf84+xyzzyaach84*xyzzyaadr84
endif
enddo
case(0)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaadr1(1)=xyzzyaadr1(0)*xyzzyaacx84
xyzzyaadr1(2)=xyzzyaadr1(1)*xyzzyaacx84
xyzzyaadr1(3)=xyzzyaadr1(2)*xyzzyaacx84
xyzzyaadp1(1)=xyzzyaadp1(0)*xyzzyaacz84
xyzzyaadp1(2)=xyzzyaadp1(1)*xyzzyaacz84
xyzzyaadp1(3)=xyzzyaadp1(2)*xyzzyaacz84
do xyzzyaaaf84=0,3
xyzzyaadp84=sum(xyzzyaadx1(0:3,xyzzyaaaf84,s)*xyzzyaadp1(0:3))
xyzzyaabf84=xyzzyaabf84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
enddo
case default
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(xyzzyaacz84<xyzzyaadn84)then
xyzzyaadr84=(xyzzyaado84*(xyzzyaacz84-xyzzyaadn84))**xyzzyaagx1
s=xyzzyaaao84(which_spin(xyzzyaaaa84))
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
xyzzyaadr1(1)=xyzzyaadr1(0)*xyzzyaacx84
xyzzyaadr1(2)=xyzzyaadr1(1)*xyzzyaacx84
xyzzyaadr1(3)=xyzzyaadr1(2)*xyzzyaacx84
xyzzyaadp1(1)=xyzzyaadp1(0)*xyzzyaacz84
xyzzyaadp1(2)=xyzzyaadp1(1)*xyzzyaacz84
xyzzyaadp1(3)=xyzzyaadp1(2)*xyzzyaacz84
xyzzyaach84=0.d0
do xyzzyaaaf84=0,3
xyzzyaadp84=sum(xyzzyaadx1(0:3,xyzzyaaaf84,s)*xyzzyaadp1(0:3))
xyzzyaach84=xyzzyaach84+xyzzyaadp84*xyzzyaadr1(xyzzyaaaf84)
enddo
xyzzyaabf84=xyzzyaabf84+xyzzyaach84*xyzzyaadr84
endif
enddo
end select
enddo
!$omp  enddo nowait
endif
end select
!$omp  master
call timer('F TERM',.false.)
!$omp  end master
endif
!$omp end parallel
elseif(dimensionality==2)then
if(xyzzyaaep1)then
call timer('U TERM',.true.)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
call xyzzyaapr1(xyzzyaacx84,ispin,xyzzyaaab84,val,xyzzyaaef84,sd,xyzzy&
&aabt84,xyzzyaabu84,xyzzyaabv84)
if(val)xyzzyaaaq84=xyzzyaaaq84+xyzzyaabt84
if(xyzzyaaef84)then
xyzzyaadl84=xyzzyaabu84/xyzzyaacx84
if(fd)xyzzyaaar84(1:2)=xyzzyaaar84(1:2)+xyzzyaadl84*eevecs1(1:2,xyzzya&
&aaa84)
if(sd)xyzzyaaas84=xyzzyaaas84+xyzzyaabv84+xyzzyaadl84
endif
enddo
call timer('U TERM',.false.)
endif
if(xyzzyaaet1)then
call timer('QCUSP TERM',.true.)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
call xyzzyaapv1(xyzzyaacx84,ispin,xyzzyaaab84,val,xyzzyaaef84,sd,xyzzy&
&aacb84,xyzzyaacc84,xyzzyaacd84)
if(val)xyzzyaaaw84=xyzzyaaaw84+xyzzyaacb84
if(xyzzyaaef84)then
xyzzyaadl84=xyzzyaacc84/xyzzyaacx84
if(fd)xyzzyaaax84(1:2)=xyzzyaaax84(1:2)+xyzzyaadl84*eevecs1(1:2,xyzzya&
&aaa84)
if(sd)xyzzyaaay84=xyzzyaaay84+xyzzyaacd84+xyzzyaadl84
endif
enddo
call timer('QCUSP TERM',.false.)
endif
if(xyzzyaaew1.and.xyzzyaaeg84)then
call timer('CHI TERM',.true.)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabo1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacy84<=0.d0)cycle
call xyzzyaaqa1(xyzzyaacy84,ispin,xyzzyaaag84,val_x_q,xyzzyaaeh84,sd,x&
&yzzyaace84,xyzzyaacf84,xyzzyaacg84)
if(val_x_q)xyzzyaabc84=xyzzyaabc84+xyzzyaace84
if(xyzzyaaeh84)then
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
if(fd_x_q)xyzzyaabd84(1:2)=xyzzyaabd84(1:2)+xyzzyaacf84*xyzzyaadc84(1:&
&2)
if(sd)xyzzyaabe84=xyzzyaabe84+xyzzyaacg84+xyzzyaacf84*xyzzyaadf84
endif
enddo
call timer('CHI TERM',.false.)
endif
if(xyzzyaaex1)then
call timer('F TERM',.true.)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84.or.xyzzyaacy84<=0.d0)cycle
xyzzyaaaf84=xyzzyaacg1(xyzzyaaag84)
if(xyzzyaaef84)then
xyzzyaadf84=1.d0/xyzzyaacy84
xyzzyaadc84=eivecs1(1:3,xyzzyaaac84)*xyzzyaadf84
endif
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(i==xyzzyaaaa84)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84)cycle
xyzzyaadd84=xyzzyaacz84
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaadd84
do xyzzyaaae84=2,xyzzyaaaf84
xyzzyaadd84=xyzzyaadd84*xyzzyaacz84
xyzzyaadq1(xyzzyaaae84,xyzzyaaaa84)=xyzzyaadd84
enddo
enddo
call dcopy(xyzzyaaaf84,xyzzyaadq1(1,i),1,xyzzyaado1(1),1)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
call xyzzyaaqb1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,xyzzyaaef84,sd,.false.,xyzzyaach84,xyzzya&
&aco84,xyzzyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyz&
&zyaacs84,xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
if(xyzzyaaef84)then
xyzzyaade84=1.d0/xyzzyaacx84
xyzzyaadb84=eevecs1(1:3,xyzzyaaaa84)*xyzzyaade84
if(fd)xyzzyaabg84(1:2)=xyzzyaabg84(1:2)+xyzzyaacp84*xyzzyaadb84(1:2)+x&
&yzzyaaco84*xyzzyaadc84(1:2)
if(sd)xyzzyaabh84=xyzzyaabh84+xyzzyaaco84*xyzzyaadf84+xyzzyaacq84+xyzz&
&yaacp84*xyzzyaade84+xyzzyaacr84+2*xyzzyaacs84*(xyzzyaadc84(1)*xyzzyaa&
&db84(1)+xyzzyaadc84(2)*xyzzyaadb84(2))
endif
enddo
enddo
call timer('F TERM',.false.)
endif
else
if(xyzzyaaep1)then
call timer('U TERM',.true.)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
call xyzzyaapr1(xyzzyaacx84,ispin,xyzzyaaab84,val,fd,sd,xyzzyaabt84,xy&
&zzyaabu84,xyzzyaabv84)
if(val)xyzzyaaaq84=xyzzyaaaq84+xyzzyaabt84
if(fd)xyzzyaaar84(1)=xyzzyaaar84(1)+xyzzyaabu84*sign(1.d0,eevecs1(1,xy&
&zzyaaaa84))
if(sd)xyzzyaaas84=xyzzyaaas84+xyzzyaabv84
enddo
call timer('U TERM',.false.)
endif
if(xyzzyaaet1)then
call timer('QCUSP TERM',.true.)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacx84=eevecs1(4,xyzzyaaaa84)
if(xyzzyaacx84<=0.d0)cycle
call xyzzyaapv1(xyzzyaacx84,ispin,xyzzyaaab84,val,fd,sd,xyzzyaacb84,xy&
&zzyaacc84,xyzzyaacd84)
if(val)xyzzyaaaw84=xyzzyaaaw84+xyzzyaacb84
if(xyzzyaaef84)then
if(fd)xyzzyaaax84(1)=xyzzyaaax84(1)+xyzzyaacc84*sign(1.d0,eevecs1(1,xy&
&zzyaaaa84))
if(sd)xyzzyaaay84=xyzzyaaay84+xyzzyaacd84
endif
enddo
call timer('QCUSP TERM',.false.)
endif
if(xyzzyaaew1.and.xyzzyaaeg84)then
call timer('CHI TERM',.true.)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabo1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
call xyzzyaaqa1(xyzzyaacy84,ispin,xyzzyaaag84,val_x_q,fd_x_q,sd,xyzzya&
&ace84,xyzzyaacf84,xyzzyaacg84)
if(val_x_q)xyzzyaabc84=xyzzyaabc84+xyzzyaace84
if(fd_x_q)xyzzyaabd84(1)=xyzzyaabd84(1)+xyzzyaacf84*sign(1.d0,eivecs1(&
&1,xyzzyaaac84))
if(sd)xyzzyaabe84=xyzzyaabe84+xyzzyaacg84
enddo
call timer('CHI TERM',.false.)
endif
if(xyzzyaaex1)then
call timer('F TERM',.true.)
do xyzzyaaac84=1,nitot
xyzzyaaag84=xyzzyaabp1(xyzzyaaac84)
if(xyzzyaaag84==0)cycle
xyzzyaacy84=eivecs1(4,xyzzyaaac84)
xyzzyaadn84=xyzzyaacv1(xyzzyaaag84)
if(xyzzyaacy84>=xyzzyaadn84)cycle
xyzzyaaaf84=xyzzyaacg1(xyzzyaaag84)
do xyzzyaaaa84=1,netot
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
if(i==xyzzyaaaa84)xyzzyaacz84=eivecs1(4,xyzzyaaac84)
if(xyzzyaacz84>=xyzzyaadn84)cycle
xyzzyaadd84=xyzzyaacz84
xyzzyaadq1(1,xyzzyaaaa84)=xyzzyaadd84
do xyzzyaaae84=2,xyzzyaaaf84
xyzzyaadd84=xyzzyaadd84*xyzzyaacz84
xyzzyaadq1(xyzzyaaae84,xyzzyaaaa84)=xyzzyaadd84
enddo
enddo
call dcopy(xyzzyaaaf84,xyzzyaadq1(1,i),1,xyzzyaado1(1),1)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaacz84=eivecs(4,xyzzyaaac84,xyzzyaaaa84)
xyzzyaacx84=abs(eevecs1(1,xyzzyaaaa84))
call xyzzyaaqb1(xyzzyaaaa84,xyzzyaacy84,xyzzyaacz84,xyzzyaacx84,ispin,&
&xyzzyaaab84,xyzzyaaag84,val,fd,sd,.false.,xyzzyaach84,xyzzyaaco84,xyz&
&zyaact84,xyzzyaacp84,xyzzyaacq84,xyzzyaacu84,xyzzyaacr84,xyzzyaacs84,&
&xyzzyaacv84)
if(val)xyzzyaabf84=xyzzyaabf84+xyzzyaach84
if(fd)xyzzyaabg84(1)=xyzzyaabg84(1)+xyzzyaacp84*sign(1.d0,eevecs1(1,xy&
&zzyaaaa84))+xyzzyaaco84*sign(1.d0,eivecs1(1,xyzzyaaac84))
if(sd)xyzzyaabh84=xyzzyaabh84+xyzzyaacq84+xyzzyaacr84+xyzzyaacs84*sign&
&(2.d0,eevecs1(1,xyzzyaaaa84))*sign(1.d0,eivecs1(1,xyzzyaaac84))
enddo
enddo
call timer('F TERM',.false.)
endif
endif
if(xyzzyaaeu1)then
call timer('W TERM',.true.)
if(.not.precomp_3body)call xyzzyaapw1(i,eevecs,xyzzyaaef84,sd)
if(xyzzyaaie1(i))then
if(val)then
xyzzyaadj84=xyzzyaael1(i)
xyzzyaaaz84=xyzzyaadj84*xyzzyaadj84-2*xyzzyaaed1(i)-3*xyzzyaaee1(i)
endif
if(fd)then
xyzzyaaba84=0.5d0*dsum3(netot,xyzzyaaej1(1,1,i),1)
do xyzzyaaaa84=1,i-1
xyzzyaaah84=which_ee(i,xyzzyaaaa84)
xyzzyaadi84=xyzzyaaem1(xyzzyaaah84)*xyzzyaaen1(:,xyzzyaaah84)
xyzzyaaba84=xyzzyaaba84+xyzzyaadi84
enddo
do xyzzyaaaa84=i+1,netot
xyzzyaaah84=which_ee(i,xyzzyaaaa84)
xyzzyaadi84=-xyzzyaaem1(xyzzyaaah84)*xyzzyaaen1(:,xyzzyaaah84)
xyzzyaaba84=xyzzyaaba84+xyzzyaadi84
enddo
xyzzyaaba84=4*xyzzyaaba84
endif
if(sd)then
xyzzyaabb84=0.5d0*dsum(netot,xyzzyaaek1(1,i),1)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaah84=which_ee(i,xyzzyaaaa84)
xyzzyaadk84=xyzzyaaem1(xyzzyaaah84)*xyzzyaaeo1(xyzzyaaah84)+ddot(dimen&
&sionality,xyzzyaaen1(1,xyzzyaaah84),1,xyzzyaaen1(1,xyzzyaaah84),1)
xyzzyaabb84=xyzzyaabb84-xyzzyaadk84
enddo
xyzzyaabb84=4*xyzzyaabb84
endif
endif
call timer('W TERM',.false.)
endif
if(xyzzyaaev1)then
call timer('H TERM',.true.)
if(precomp_3body)then
call xyzzyaapy1(i,xyzzyaaig1,val,fd,sd,xyzzyaaaz84,xyzzyaaba84,xyzzyaa&
&bb84)
else
call xyzzyaapy1(i,eevecs,val,fd,sd,xyzzyaaaz84,xyzzyaaba84,xyzzyaabb84&
&)
endif
call timer('H TERM',.false.)
endif
if(xyzzyaaey1)then
call timer('P TERM',.true.)
do xyzzyaaaa84=1,netot
if(xyzzyaaaa84==i)cycle
xyzzyaaab84=which_spin(xyzzyaaaa84)
xyzzyaada84=eevecs1(1:3,xyzzyaaaa84)
call xyzzyaaqd1(xyzzyaada84,ispin,xyzzyaaab84,val,fd.or.sd,sd,xyzzyaac&
&i84,xyzzyaacj84,xyzzyaack84)
if(val)xyzzyaabi84=xyzzyaabi84+xyzzyaaci84
if(fd)xyzzyaabj84(1:3)=xyzzyaabj84(1:3)+xyzzyaacj84(1:3)
if(sd)xyzzyaabk84=xyzzyaabk84+xyzzyaack84
enddo
call timer('P TERM',.false.)
endif
if(xyzzyaaez1.and.xyzzyaaeg84)then
call timer('Q TERM',.true.)
call xyzzyaaqe1(rvec,ispin,val_x_q,fd_x_q,sd,xyzzyaabl84,xyzzyaabm84,x&
&yzzyaabn84)
call timer('Q TERM',.false.)
endif
if(xyzzyaafa1)then
call timer('BIEX1',.true.)
call xyzzyaaqf1(val,fd,sd,ispin,eevecs1,total_b,xyzzyaabo84,xyzzyaabp8&
&4)
call timer('BIEX1',.false.)
elseif(xyzzyaafb1)then
call timer('BIEX2',.true.)
if(precomp_3body)then
call xyzzyaaqg1(val,fd,sd,ispin,xyzzyaaig1,total_b,xyzzyaabo84,xyzzyaa&
&bp84)
else
call xyzzyaaqg1(val,fd,sd,ispin,eevecs,total_b,xyzzyaabo84,xyzzyaabp84&
&)
endif
call timer('BIEX2',.false.)
elseif(xyzzyaafc1)then
call timer('BIEX3',.true.)
do xyzzyaaaa84=1,i-1
xyzzyaaee84(1:3,xyzzyaaaa84)=rvec-eevecs1(1:3,xyzzyaaaa84)
enddo
xyzzyaaee84(1:3,i)=rvec
do xyzzyaaaa84=i+1,netot
xyzzyaaee84(1:3,xyzzyaaaa84)=rvec-eevecs1(1:3,xyzzyaaaa84)
enddo
if(fix_holes)then
call xyzzyaaql1(val,fd,sd,ispin,xyzzyaaee84,total_b,xyzzyaabo84,xyzzya&
&abp84)
else
call xyzzyaaqi1(val,fd,sd,ispin,xyzzyaaee84,total_b,xyzzyaabo84,xyzzya&
&abp84)
endif
call timer('BIEX3',.false.)
elseif(xyzzyaafg1)then
call timer('EX2D',.true.)
if(precomp_3body)then
call xyzzyaaqp1(val,fd,sd,ispin,eevecs1,xyzzyaaig1,eivecs1,eivecs,tota&
&l_b,xyzzyaabo84,xyzzyaabp84)
else
call xyzzyaaqp1(val,fd,sd,ispin,eevecs1,eevecs,eivecs1,eivecs,total_b,&
&xyzzyaabo84,xyzzyaabp84)
endif
call timer('EX2D',.false.)
endif
if(xyzzyaaff1)then
call timer('D TERM',.true.)
do xyzzyaaac84=1,nitot
xyzzyaadg84=eivecs1(:,xyzzyaaac84)
if(xyzzyaadg84(4)>xyzzyaacp1)cycle
do xyzzyaaak84=1,xyzzyaabm1
if(xyzzyaabu1(xyzzyaaak84,1)==xyzzyaaac84)then
xyzzyaaaj84=xyzzyaabu1(xyzzyaaak84,2)
do xyzzyaaaa84=1,netot
if(i==xyzzyaaaa84)cycle
if(eivecs(4,xyzzyaaaj84,xyzzyaaaa84)>=xyzzyaacp1)cycle
xyzzyaadh84=eivecs(:,xyzzyaaaj84,xyzzyaaaa84)
call xyzzyaaqq1(xyzzyaadg84,xyzzyaadh84,xyzzyaaak84,1,val,fd,sd,xyzzya&
&acl84,xyzzyaacm84,xyzzyaacn84)
if(val)xyzzyaabq84=xyzzyaabq84+xyzzyaacl84
if(fd)xyzzyaabr84=xyzzyaabr84+xyzzyaacm84
if(sd)xyzzyaabs84=xyzzyaabs84+xyzzyaacn84
enddo
elseif(xyzzyaabu1(xyzzyaaak84,2)==xyzzyaaac84)then
xyzzyaaaj84=xyzzyaabu1(xyzzyaaak84,1)
do xyzzyaaaa84=1,netot
if(i==xyzzyaaaa84)cycle
if(eivecs(4,xyzzyaaaj84,xyzzyaaaa84)>=xyzzyaacp1)cycle
xyzzyaadh84=eivecs(:,xyzzyaaaj84,xyzzyaaaa84)
call xyzzyaaqq1(xyzzyaadh84,xyzzyaadg84,xyzzyaaak84,2,val,fd,sd,xyzzya&
&acl84,xyzzyaacm84,xyzzyaacn84)
if(val)xyzzyaabq84=xyzzyaabq84+xyzzyaacl84
if(fd)xyzzyaabr84=xyzzyaabr84+xyzzyaacm84
if(sd)xyzzyaabs84=xyzzyaabs84+xyzzyaacn84
enddo
endif
enddo
enddo
call timer('D TERM',.false.)
endif
if(val)value_jas=xyzzyaaaq84+xyzzyaaat84+xyzzyaaaw84+xyzzyaaaz84+xyzzy&
&aabc84+xyzzyaabf84+xyzzyaabi84+xyzzyaabl84+total_b+xyzzyaabq84
if(fd)grad_jas=xyzzyaaar84+xyzzyaaau84+xyzzyaaax84+xyzzyaaba84+xyzzyaa&
&bd84+xyzzyaabg84+xyzzyaabj84+xyzzyaabm84+xyzzyaabo84+xyzzyaabr84
if(sd)lap_jas=xyzzyaaas84+xyzzyaaav84+xyzzyaaay84+xyzzyaabb84+xyzzyaab&
&e84+xyzzyaabh84+xyzzyaabk84+xyzzyaabn84+xyzzyaabp84+xyzzyaabs84
if(val_x_q)value_x_q=xyzzyaabc84+xyzzyaabl84
if(fd_x_q)grad_x_q(1:3)=xyzzyaabd84(1:3)+xyzzyaabm84(1:3)
value_3body=xyzzyaaaz84
if(xyzzyaaaa1.and.val)call wout('oneelec_jastrow: value_jas=',value_ja&
&s)
if(xyzzyaaaa1.and.fd)call wout('oneelec_jastrow: grad_jas=',grad_jas)
if(xyzzyaaaa1.and.sd)call wout('oneelec_jastrow: lap_jas=',lap_jas)
call timer('JASTROW',.false.)
end subroutine xyzzyaapp1
subroutine xyzzyaapq1(i,ispin,rvec,eevecs,eevecs1,eivecs,val,fd,sd,pre&
&comp_3body,total_u,grad_u,lap_u,total_ucyl,grad_ucyl,lap_ucyl,total_q&
&cusp,grad_qcusp,lap_qcusp,total_w,grad_w,lap_w,total_x,grad_x,lap_x,t&
&otal_f,grad_f,lap_f,total_p,grad_p,lap_p,total_q,grad_q,lap_q,total_b&
&,grad_b,lap_b,total_d,grad_d,lap_d,want_u,want_ucyl,want_qcusp,want_3&
&body,want_chi,want_f,want_p,want_q,want_b,want_d)
use slaarnabt, only : ddot,dcopy,dsum,dsum3
implicit none
integer,intent(in) :: i,ispin
real(dp),intent(in) :: rvec(3),eevecs(4,netot,netot),eevecs1(4,netot),&
&eivecs(4,nitot,netot)
real(dp),intent(out) :: total_u,grad_u(3),lap_u,total_ucyl,grad_ucyl(3&
&),lap_ucyl,total_qcusp,grad_qcusp(3),lap_qcusp,total_w,grad_w(3),lap_&
&w,total_x,grad_x(3),lap_x,total_f,grad_f(3),lap_f,total_p,grad_p(3),l&
&ap_p,total_q,grad_q(3),lap_q,total_b,grad_b(3),lap_b,total_d,grad_d(3&
&),lap_d
logical,intent(in) :: val,fd,sd,want_u,want_ucyl,want_qcusp,want_3body&
&,want_chi,want_f,want_p,want_q,want_b,want_d,precomp_3body
integer xyzzyaaaa85,xyzzyaaab85,xyzzyaaac85,xyzzyaaad85,xyzzyaaae85,xy&
&zzyaaaf85,xyzzyaaag85,xyzzyaaah85,xyzzyaaai85,xyzzyaaaj85
real(dp) xyzzyaaak85,xyzzyaaal85,xyzzyaaam85,xyzzyaaan85,xyzzyaaao85,x&
&yzzyaaap85,xyzzyaaaq85,xyzzyaaar85,xyzzyaaas85,xyzzyaaat85,xyzzyaaau8&
&5,xyzzyaaav85,xyzzyaaaw85,xyzzyaaax85,xyzzyaaay85,xyzzyaaaz85,xyzzyaa&
&ba85(3),xyzzyaabb85,xyzzyaabc85,xyzzyaabd85,xyzzyaabe85,xyzzyaabf85,x&
&yzzyaabg85,xyzzyaabh85,xyzzyaabi85,xyzzyaabj85,xyzzyaabk85,xyzzyaabl8&
&5,xyzzyaabm85,xyzzyaabn85,xyzzyaabo85,xyzzyaabp85(3),xyzzyaabq85(3),x&
&yzzyaabr85(3),xyzzyaabs85,xyzzyaabt85,xyzzyaabu85,xyzzyaabv85(3),xyzz&
&yaabw85,xyzzyaabx85,xyzzyaaby85,xyzzyaabz85,xyzzyaaca85,xyzzyaacb85,x&
&yzzyaacc85(4),xyzzyaacd85(4),xyzzyaace85,xyzzyaacf85(3),xyzzyaacg85
real(dp),allocatable,save :: xyzzyaach85(:,:)
logical xyzzyaaci85
logical,save :: xyzzyaacj85=.true.
call timer('JASTROW',.true.)
if(xyzzyaacj85)then
if(xyzzyaafc1)then
allocate(xyzzyaach85(3,netot),stat=xyzzyaaai85)
call check_alloc(xyzzyaaai85,'ONEELEC_JASTROW_SPLIT','')
endif
xyzzyaacj85=.false.
endif
total_u=0.d0
total_ucyl=0.d0
total_w=0.d0
total_x=0.d0
total_qcusp=0.d0
total_f=0.d0
total_p=0.d0
total_q=0.d0
total_b=0.d0
grad_u=0.d0
grad_ucyl=0.d0
grad_w=0.d0
grad_x=0.d0
grad_qcusp=0.d0
grad_f=0.d0
grad_p=0.d0
grad_q=0.d0
grad_b=0.d0
lap_u=0.d0
lap_ucyl=0.d0
lap_w=0.d0
lap_x=0.d0
lap_qcusp=0.d0
lap_f=0.d0
lap_p=0.d0
lap_q=0.d0
lap_b=0.d0
total_d=0.d0
grad_d=0.d0
lap_d=0.d0
xyzzyaaci85=fd.or.sd
if(dimensionality==3)then
if(xyzzyaaep1.and.want_u)then
call timer('U TERM',.true.)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabl85=eevecs1(4,xyzzyaaaa85)
if(xyzzyaabl85<=0.d0)cycle
call xyzzyaapr1(xyzzyaabl85,ispin,xyzzyaaab85,val,xyzzyaaci85,sd,xyzzy&
&aaak85,xyzzyaaal85,xyzzyaaam85)
if(val)total_u=total_u+xyzzyaaak85
if(xyzzyaaci85)then
xyzzyaaby85=xyzzyaaal85/xyzzyaabl85
if(fd)grad_u(1:3)=grad_u(1:3)+xyzzyaaby85*eevecs1(1:3,xyzzyaaaa85)
if(sd)lap_u=lap_u+xyzzyaaam85+2*xyzzyaaby85
endif
enddo
call timer('U TERM',.false.)
endif
if(xyzzyaaes1.and.want_ucyl)then
call timer('Ucyl TERM',.true.)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
if(eevecs1(4,xyzzyaaaa85)<=0.d0)cycle
xyzzyaaca85=dot_product(eevecs1(1:3,xyzzyaaaa85),xyzzyaacr1)
xyzzyaabz85=eevecs1(4,xyzzyaaaa85)*eevecs1(4,xyzzyaaaa85)-xyzzyaaca85*&
&xyzzyaaca85
if(xyzzyaabz85<xyzzyaacj1.or.xyzzyaagx1==0)then
xyzzyaabz85=sqrt(xyzzyaabz85)
call xyzzyaapu1(xyzzyaabz85,abs(xyzzyaaca85),ispin,xyzzyaaab85,val,xyz&
&zyaaci85,sd,xyzzyaaan85,xyzzyaaao85,xyzzyaaap85,xyzzyaaaq85,xyzzyaaar&
&85)
if(val)total_ucyl=total_ucyl+xyzzyaaan85
if(xyzzyaaci85)then
if(xyzzyaabz85>0.d0)then
xyzzyaacb85=xyzzyaaao85/xyzzyaabz85
else
xyzzyaacb85=0.d0
endif
if(fd)grad_ucyl=grad_ucyl+xyzzyaacb85*(eevecs1(1:3,xyzzyaaaa85)-xyzzya&
&aca85*xyzzyaacr1)+xyzzyaaap85*sign(1.d0,xyzzyaaca85)*xyzzyaacr1
if(sd)lap_ucyl=lap_ucyl+xyzzyaaaq85+xyzzyaacb85+xyzzyaaar85
endif
endif
enddo
call timer('Ucyl TERM',.false.)
endif
if(xyzzyaaew1.and.want_chi)then
call timer('CHI TERM',.true.)
do xyzzyaaac85=1,nitot
xyzzyaaag85=xyzzyaabo1(xyzzyaaac85)
if(xyzzyaaag85==0)cycle
xyzzyaabm85=eivecs(4,xyzzyaaac85,i)
if(xyzzyaabm85<=0.d0)cycle
call xyzzyaaqa1(xyzzyaabm85,ispin,xyzzyaaag85,val,xyzzyaaci85,sd,xyzzy&
&aaav85,xyzzyaaaw85,xyzzyaaax85)
if(val)total_x=total_x+xyzzyaaav85
if(xyzzyaaci85)then
xyzzyaaby85=xyzzyaaaw85/xyzzyaabm85
if(fd)grad_x(1:3)=grad_x(1:3)+xyzzyaaby85*eivecs(1:3,xyzzyaaac85,i)
if(sd)lap_x=lap_x+xyzzyaaax85+2*xyzzyaaby85
endif
enddo
call timer('CHI TERM',.false.)
endif
if(xyzzyaaex1.and.want_f)then
call timer('F TERM',.true.)
do xyzzyaaac85=1,nitot
xyzzyaaag85=xyzzyaabp1(xyzzyaaac85)
if(xyzzyaaag85==0)cycle
xyzzyaabm85=eivecs(4,xyzzyaaac85,i)
xyzzyaabo85=xyzzyaacv1(xyzzyaaag85)
if(xyzzyaabm85>=xyzzyaabo85.or.xyzzyaabm85<=0.d0)cycle
if(xyzzyaaci85)then
xyzzyaabu85=1.d0/xyzzyaabm85
xyzzyaabr85=eivecs(1:3,xyzzyaaac85,i)*xyzzyaabu85
endif
xyzzyaaae85=xyzzyaacg1(xyzzyaaag85)
do xyzzyaaaa85=1,netot
xyzzyaabn85=eivecs(4,xyzzyaaac85,xyzzyaaaa85)
if(xyzzyaabn85>=xyzzyaabo85)cycle
xyzzyaabs85=xyzzyaabn85
xyzzyaadq1(1,xyzzyaaaa85)=xyzzyaabs85
do xyzzyaaad85=2,xyzzyaaae85
xyzzyaabs85=xyzzyaabs85*xyzzyaabn85
xyzzyaadq1(xyzzyaaad85,xyzzyaaaa85)=xyzzyaabs85
enddo
enddo
call dcopy(xyzzyaaae85,xyzzyaadq1(1,i),1,xyzzyaado1(1),1)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabn85=eivecs(4,xyzzyaaac85,xyzzyaaaa85)
xyzzyaabl85=eevecs1(4,xyzzyaaaa85)
if(xyzzyaabl85<=0.d0)cycle
call xyzzyaaqb1(xyzzyaaaa85,xyzzyaabm85,xyzzyaabn85,xyzzyaabl85,ispin,&
&xyzzyaaab85,xyzzyaaag85,val,xyzzyaaci85,sd,.false.,xyzzyaaay85,xyzzya&
&abc85,xyzzyaabh85,xyzzyaabd85,xyzzyaabe85,xyzzyaabi85,xyzzyaabf85,xyz&
&zyaabg85,xyzzyaabj85)
if(val)total_f=total_f+xyzzyaaay85
if(xyzzyaaci85)then
xyzzyaabt85=1.d0/xyzzyaabl85
xyzzyaabq85=eevecs1(1:3,xyzzyaaaa85)*xyzzyaabt85
if(fd)grad_f(1:3)=grad_f(1:3)+xyzzyaabd85*xyzzyaabq85(1:3)+xyzzyaabc85&
&*xyzzyaabr85(1:3)
if(sd)then
xyzzyaabk85=xyzzyaabr85(1)*xyzzyaabq85(1)+xyzzyaabr85(2)*xyzzyaabq85(2&
&)+xyzzyaabr85(3)*xyzzyaabq85(3)
lap_f=lap_f+2*xyzzyaabc85*xyzzyaabu85+xyzzyaabe85+2*xyzzyaabd85*xyzzya&
&abt85+xyzzyaabf85+2*xyzzyaabk85*xyzzyaabg85
endif
endif
enddo
enddo
call timer('F TERM',.false.)
endif
elseif(dimensionality==2)then
if(xyzzyaaep1.and.want_u)then
call timer('U TERM',.true.)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabl85=eevecs1(4,xyzzyaaaa85)
if(xyzzyaabl85<=0.d0)cycle
call xyzzyaapr1(xyzzyaabl85,ispin,xyzzyaaab85,val,xyzzyaaci85,sd,xyzzy&
&aaak85,xyzzyaaal85,xyzzyaaam85)
if(val)total_u=total_u+xyzzyaaak85
if(xyzzyaaci85)then
xyzzyaaby85=xyzzyaaal85/xyzzyaabl85
if(fd)grad_u(1:2)=grad_u(1:2)+xyzzyaaby85*eevecs1(1:2,xyzzyaaaa85)
if(sd)lap_u=lap_u+xyzzyaaam85+xyzzyaaby85
endif
enddo
call timer('U TERM',.false.)
endif
if(xyzzyaaet1.and.want_qcusp)then
call timer('QCUSP TERM',.true.)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabl85=eevecs1(4,xyzzyaaaa85)
if(xyzzyaabl85<=0.d0)cycle
call xyzzyaapv1(xyzzyaabl85,ispin,xyzzyaaab85,val,xyzzyaaci85,sd,xyzzy&
&aaas85,xyzzyaaat85,xyzzyaaau85)
if(val)total_qcusp=total_qcusp+xyzzyaaas85
if(xyzzyaaci85)then
xyzzyaaby85=xyzzyaaat85/xyzzyaabl85
if(fd)grad_qcusp(1:2)=grad_qcusp(1:2)+xyzzyaaby85*eevecs1(1:2,xyzzyaaa&
&a85)
if(sd)lap_qcusp=lap_qcusp+xyzzyaaau85+xyzzyaaby85
endif
enddo
call timer('QCUSP TERM',.false.)
endif
if(xyzzyaaew1.and.want_chi)then
call timer('CHI TERM',.true.)
do xyzzyaaac85=1,nitot
xyzzyaaag85=xyzzyaabo1(xyzzyaaac85)
if(xyzzyaaag85==0)cycle
xyzzyaabm85=eivecs(4,xyzzyaaac85,i)
if(xyzzyaabm85<=0.d0)cycle
call xyzzyaaqa1(xyzzyaabm85,ispin,xyzzyaaag85,val,xyzzyaaci85,sd,xyzzy&
&aaav85,xyzzyaaaw85,xyzzyaaax85)
if(val)total_x=total_x+xyzzyaaav85
if(xyzzyaaci85)then
xyzzyaabu85=1.d0/xyzzyaabm85
xyzzyaabr85=eivecs(1:3,xyzzyaaac85,i)*xyzzyaabu85
if(fd)grad_x(1:2)=grad_x(1:2)+xyzzyaaaw85*xyzzyaabr85(1:2)
if(sd)lap_x=lap_x+xyzzyaaax85+xyzzyaaaw85*xyzzyaabu85
endif
enddo
call timer('CHI TERM',.false.)
endif
if(xyzzyaaex1.and.want_f)then
call timer('F TERM',.true.)
do xyzzyaaac85=1,nitot
xyzzyaaag85=xyzzyaabp1(xyzzyaaac85)
if(xyzzyaaag85==0)cycle
xyzzyaabm85=eivecs(4,xyzzyaaac85,i)
xyzzyaabo85=xyzzyaacv1(xyzzyaaag85)
if(xyzzyaabm85>=xyzzyaabo85.or.xyzzyaabm85<=0.d0)cycle
xyzzyaaae85=xyzzyaacg1(xyzzyaaag85)
if(xyzzyaaci85)then
xyzzyaabu85=1.d0/xyzzyaabm85
xyzzyaabr85=eivecs(1:3,xyzzyaaac85,i)*xyzzyaabu85
endif
do xyzzyaaaa85=1,netot
xyzzyaabn85=eivecs(4,xyzzyaaac85,xyzzyaaaa85)
if(xyzzyaabn85>=xyzzyaabo85)cycle
xyzzyaabs85=xyzzyaabn85
xyzzyaadq1(1,xyzzyaaaa85)=xyzzyaabs85
do xyzzyaaad85=2,xyzzyaaae85
xyzzyaabs85=xyzzyaabs85*xyzzyaabn85
xyzzyaadq1(xyzzyaaad85,xyzzyaaaa85)=xyzzyaabs85
enddo
enddo
call dcopy(xyzzyaaae85,xyzzyaadq1(1,i),1,xyzzyaado1(1),1)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabn85=eivecs(4,xyzzyaaac85,xyzzyaaaa85)
xyzzyaabl85=eevecs1(4,xyzzyaaaa85)
if(xyzzyaabl85<=0.d0)cycle
call xyzzyaaqb1(xyzzyaaaa85,xyzzyaabm85,xyzzyaabn85,xyzzyaabl85,ispin,&
&xyzzyaaab85,xyzzyaaag85,val,xyzzyaaci85,sd,.false.,xyzzyaaay85,xyzzya&
&abc85,xyzzyaabh85,xyzzyaabd85,xyzzyaabe85,xyzzyaabi85,xyzzyaabf85,xyz&
&zyaabg85,xyzzyaabj85)
if(val)total_f=total_f+xyzzyaaay85
if(xyzzyaaci85)then
xyzzyaabt85=1.d0/xyzzyaabl85
xyzzyaabq85=eevecs1(1:3,xyzzyaaaa85)*xyzzyaabt85
if(fd)grad_f(1:2)=grad_f(1:2)+xyzzyaabd85*xyzzyaabq85(1:2)+xyzzyaabc85&
&*xyzzyaabr85(1:2)
if(sd)lap_f=lap_f+xyzzyaabc85*xyzzyaabu85+xyzzyaabe85+xyzzyaabd85*xyzz&
&yaabt85+xyzzyaabf85+2*xyzzyaabg85*(xyzzyaabr85(1)*xyzzyaabq85(1)+xyzz&
&yaabr85(2)*xyzzyaabq85(2))
endif
enddo
enddo
call timer('F TERM',.false.)
endif
else
if(xyzzyaaep1.and.want_u)then
call timer('U TERM',.true.)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabl85=abs(eevecs1(1,xyzzyaaaa85))
call xyzzyaapr1(xyzzyaabl85,ispin,xyzzyaaab85,val,fd,sd,xyzzyaaak85,xy&
&zzyaaal85,xyzzyaaam85)
if(val)total_u=total_u+xyzzyaaak85
if(fd)grad_u(1)=grad_u(1)+xyzzyaaal85*sign(1.d0,eevecs1(1,xyzzyaaaa85)&
&)
if(sd)lap_u=lap_u+xyzzyaaam85
enddo
call timer('U TERM',.false.)
endif
if(xyzzyaaet1.and.want_qcusp)then
call timer('QCUSP TERM',.true.)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabl85=eevecs1(4,xyzzyaaaa85)
if(xyzzyaabl85<=0.d0)cycle
call xyzzyaapv1(xyzzyaabl85,ispin,xyzzyaaab85,val,fd,sd,xyzzyaaas85,xy&
&zzyaaat85,xyzzyaaau85)
if(val)total_qcusp=total_qcusp+xyzzyaaas85
if(xyzzyaaci85)then
if(fd)grad_qcusp(1)=grad_qcusp(1)+xyzzyaaat85*sign(1.d0,eevecs1(1,xyzz&
&yaaaa85))
if(sd)lap_qcusp=lap_qcusp+xyzzyaaau85
endif
enddo
call timer('QCUSP TERM',.false.)
endif
if(xyzzyaaew1.and.want_chi)then
call timer('CHI TERM',.true.)
do xyzzyaaac85=1,nitot
xyzzyaaag85=xyzzyaabo1(xyzzyaaac85)
if(xyzzyaaag85==0)cycle
xyzzyaabm85=eivecs(4,xyzzyaaac85,i)
call xyzzyaaqa1(xyzzyaabm85,ispin,xyzzyaaag85,val,fd,sd,xyzzyaaav85,xy&
&zzyaaaw85,xyzzyaaax85)
if(val)total_x=total_x+xyzzyaaav85
if(fd)grad_x(1)=grad_x(1)+xyzzyaaaw85*sign(1.d0,eivecs(1,xyzzyaaac85,i&
&))
if(sd)lap_x=lap_x+xyzzyaaax85
enddo
call timer('CHI TERM',.false.)
endif
if(xyzzyaaex1.and.want_f)then
call timer('F TERM',.true.)
do xyzzyaaac85=1,nitot
xyzzyaaag85=xyzzyaabp1(xyzzyaaac85)
if(xyzzyaaag85==0)cycle
xyzzyaabm85=eivecs(4,xyzzyaaac85,i)
xyzzyaabo85=xyzzyaacv1(xyzzyaaag85)
if(xyzzyaabm85>=xyzzyaabo85)cycle
xyzzyaaae85=xyzzyaacg1(xyzzyaaag85)
do xyzzyaaaa85=1,netot
xyzzyaabn85=eivecs(4,xyzzyaaac85,xyzzyaaaa85)
if(xyzzyaabn85>=xyzzyaabo85)cycle
xyzzyaabs85=xyzzyaabn85
xyzzyaadq1(1,xyzzyaaaa85)=xyzzyaabs85
do xyzzyaaad85=2,xyzzyaaae85
xyzzyaabs85=xyzzyaabs85*xyzzyaabn85
xyzzyaadq1(xyzzyaaad85,xyzzyaaaa85)=xyzzyaabs85
enddo
enddo
call dcopy(xyzzyaaae85,xyzzyaadq1(1,i),1,xyzzyaado1(1),1)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabn85=eivecs(4,xyzzyaaac85,xyzzyaaaa85)
xyzzyaabl85=abs(eevecs1(1,xyzzyaaaa85))
call xyzzyaaqb1(xyzzyaaaa85,xyzzyaabm85,xyzzyaabn85,xyzzyaabl85,ispin,&
&xyzzyaaab85,xyzzyaaag85,val,fd,sd,.false.,xyzzyaaay85,xyzzyaabc85,xyz&
&zyaabh85,xyzzyaabd85,xyzzyaabe85,xyzzyaabi85,xyzzyaabf85,xyzzyaabg85,&
&xyzzyaabj85)
if(val)total_f=total_f+xyzzyaaay85
if(fd)grad_f(1)=grad_f(1)+xyzzyaabd85*sign(1.d0,eevecs1(1,xyzzyaaaa85)&
&)+xyzzyaabc85*sign(1.d0,eivecs(1,xyzzyaaac85,i))
if(sd)lap_f=lap_f+xyzzyaabe85+xyzzyaabf85+xyzzyaabg85*sign(2.d0,eevecs&
&1(1,xyzzyaaaa85))*sign(1.d0,eivecs(1,xyzzyaaac85,i))
enddo
enddo
call timer('F TERM',.false.)
endif
endif
if(xyzzyaaeu1.and.want_3body)then
call timer('W TERM',.true.)
if(.not.precomp_3body)call xyzzyaapw1(i,eevecs,xyzzyaaci85,sd)
if(xyzzyaaie1(i))then
if(val)then
xyzzyaabw85=xyzzyaael1(i)
total_w=xyzzyaabw85*xyzzyaabw85-2*xyzzyaaed1(i)-3*xyzzyaaee1(i)
endif
if(fd)then
grad_w=0.5d0*dsum3(netot,xyzzyaaej1(1,1,i),1)
do xyzzyaaaa85=1,i-1
xyzzyaaah85=which_ee(i,xyzzyaaaa85)
xyzzyaabv85=xyzzyaaem1(xyzzyaaah85)*xyzzyaaen1(:,xyzzyaaah85)
grad_w=grad_w+xyzzyaabv85
enddo
do xyzzyaaaa85=i+1,netot
xyzzyaaah85=which_ee(i,xyzzyaaaa85)
xyzzyaabv85=-xyzzyaaem1(xyzzyaaah85)*xyzzyaaen1(:,xyzzyaaah85)
grad_w=grad_w+xyzzyaabv85
enddo
grad_w=4*grad_w
endif
if(sd)then
lap_w=0.5d0*dsum(netot,xyzzyaaek1(1,i),1)
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaah85=which_ee(i,xyzzyaaaa85)
xyzzyaabx85=xyzzyaaem1(xyzzyaaah85)*xyzzyaaeo1(xyzzyaaah85)+ddot(dimen&
&sionality,xyzzyaaen1(1,xyzzyaaah85),1,xyzzyaaen1(1,xyzzyaaah85),1)
lap_w=lap_w-xyzzyaabx85
enddo
lap_w=4*lap_w
endif
endif
call timer('W TERM',.false.)
endif
if(xyzzyaaev1.and.want_3body)then
call timer('H TERM',.true.)
if(precomp_3body)then
call xyzzyaapy1(i,xyzzyaaig1,val,fd,sd,total_w,grad_w,lap_w)
else
call xyzzyaapy1(i,eevecs,val,fd,sd,total_w,grad_w,lap_w)
endif
call timer('H TERM',.false.)
endif
if(xyzzyaaey1.and.want_p)then
do xyzzyaaaa85=1,netot
if(xyzzyaaaa85==i)cycle
xyzzyaaab85=which_spin(xyzzyaaaa85)
xyzzyaabp85=eevecs1(1:3,xyzzyaaaa85)
call xyzzyaaqd1(xyzzyaabp85,ispin,xyzzyaaab85,val,fd.or.sd,sd,xyzzyaaa&
&z85,xyzzyaaba85,xyzzyaabb85)
if(val)total_p=total_p+xyzzyaaaz85
if(fd)grad_p(1:3)=grad_p(1:3)+xyzzyaaba85(1:3)
if(sd)lap_p=lap_p+xyzzyaabb85
enddo
endif
if(xyzzyaaez1.and.want_q)then
call timer('Q TERM',.true.)
call xyzzyaaqe1(rvec,ispin,val,fd,sd,total_q,grad_q,lap_q)
call timer('Q TERM',.false.)
endif
if(want_b)then
if(xyzzyaafa1)then
call timer('BIEX1',.true.)
call xyzzyaaqf1(val,fd,sd,ispin,eevecs1,total_b,grad_b,lap_b)
call timer('BIEX1',.false.)
elseif(xyzzyaafb1)then
call timer('BIEX2',.true.)
call xyzzyaaqg1(val,fd,sd,ispin,xyzzyaaig1,total_b,grad_b,lap_b)
call timer('BIEX2',.false.)
elseif(xyzzyaafc1)then
call timer('BIEX3',.true.)
do xyzzyaaaa85=1,i-1
xyzzyaach85(1:3,xyzzyaaaa85)=rvec-eevecs1(1:3,xyzzyaaaa85)
enddo
xyzzyaach85(1:3,i)=rvec
do xyzzyaaaa85=i+1,netot
xyzzyaach85(1:3,xyzzyaaaa85)=rvec-eevecs1(1:3,xyzzyaaaa85)
enddo
if(fix_holes)then
call xyzzyaaql1(val,fd,sd,ispin,xyzzyaach85,total_b,grad_b,lap_b)
else
call xyzzyaaqi1(val,fd,sd,ispin,xyzzyaach85,total_b,grad_b,lap_b)
endif
call timer('BIEX3',.false.)
elseif(xyzzyaafg1)then
call timer('EX2D',.true.)
call xyzzyaaqp1(val,fd,sd,ispin,eevecs1,xyzzyaaig1,eivecs(:,:,i),eivec&
&s,total_b,grad_b,lap_b)
call timer('EX2D',.false.)
endif
endif
if(xyzzyaaff1.and.want_d)then
call timer('D TERM',.true.)
do xyzzyaaac85=1,nitot
xyzzyaacc85=eivecs(:,xyzzyaaac85,i)
if(xyzzyaacc85(4)>xyzzyaacp1)cycle
do xyzzyaaaf85=1,xyzzyaabm1
if(xyzzyaabu1(xyzzyaaaf85,1)==xyzzyaaac85)then
xyzzyaaaj85=xyzzyaabu1(xyzzyaaaf85,2)
do xyzzyaaaa85=1,netot
if(i==xyzzyaaaa85)cycle
if(eivecs(4,xyzzyaaaj85,xyzzyaaaa85)>=xyzzyaacp1)cycle
xyzzyaacd85=eivecs(:,xyzzyaaaj85,xyzzyaaaa85)
call xyzzyaaqq1(xyzzyaacc85,xyzzyaacd85,xyzzyaaaf85,1,val,fd,sd,xyzzya&
&ace85,xyzzyaacf85,xyzzyaacg85)
if(val)total_d=total_d+xyzzyaace85
if(fd)grad_d=grad_d+xyzzyaacf85
if(sd)lap_d=lap_d+xyzzyaacg85
enddo
elseif(xyzzyaabu1(xyzzyaaaf85,2)==xyzzyaaac85)then
xyzzyaaaj85=xyzzyaabu1(xyzzyaaaf85,1)
do xyzzyaaaa85=1,netot
if(i==xyzzyaaaa85)cycle
if(eivecs(4,xyzzyaaaj85,xyzzyaaaa85)>=xyzzyaacp1)cycle
xyzzyaacd85=eivecs(:,xyzzyaaaj85,xyzzyaaaa85)
call xyzzyaaqq1(xyzzyaacd85,xyzzyaacc85,xyzzyaaaf85,2,val,fd,sd,xyzzya&
&ace85,xyzzyaacf85,xyzzyaacg85)
if(val)total_d=total_d+xyzzyaace85
if(fd)grad_d=grad_d+xyzzyaacf85
if(sd)lap_d=lap_d+xyzzyaacg85
enddo
endif
enddo
enddo
call timer('D TERM',.false.)
endif
call timer('JASTROW',.false.)
end subroutine xyzzyaapq1
subroutine xyzzyaapr1(rij,ispin,jspin,val,fd,sd,value_u,deriv_u,sderiv&
&_u)
implicit none
integer,intent(in) :: ispin,jspin
real(dp),intent(in) :: rij
real(dp),intent(out) :: value_u,deriv_u,sderiv_u
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa86
real(dp) xyzzyaaab86,xyzzyaaac86,xyzzyaaad86,xyzzyaaae86,xyzzyaaaf86,x&
&yzzyaaag86,xyzzyaaah86
if(rij>=xyzzyaach1.and.xyzzyaagx1>0)then
value_u=0.d0
deriv_u=0.d0
sderiv_u=0.d0
else
if(xyzzyaaeq1)then
xyzzyaaaa86=which_spair(ispin,jspin,levels_spairs)
call xyzzyaaps1(rij,xyzzyaaaa86,val,fd,sd,value_u,deriv_u,sderiv_u)
else
xyzzyaaaa86=which_spair(ispin,jspin,xyzzyaafh1)
call xyzzyaapt1(rij,xyzzyaaaa86,val,fd,sd,value_u,deriv_u,sderiv_u)
endif
endif
if(hard_sphere)then
if(.not.hard_op_spins.or.ispin/=jspin)then
if(rij<xyzzyaanz1)then
if(rij<=hard_diam)call errstop('COMPUTE_U','Hard spheres are overlappi&
&ng.  This is not supposed to happen.')
xyzzyaaae86=rij*xyzzyaaoa1
xyzzyaaaf86=rij*xyzzyaaob1
xyzzyaaab86=tanh((xyzzyaaae86-1.d0)/(1.d0-xyzzyaaaf86))
if(val)value_u=value_u+log(xyzzyaaab86)
if(fd.or.sd)then
xyzzyaaac86=1.d0/xyzzyaaab86
xyzzyaaad86=xyzzyaaac86-xyzzyaaab86
xyzzyaaag86=1.d0/(1.d0-xyzzyaaaf86)
xyzzyaaah86=xyzzyaaag86*xyzzyaaag86
if(fd)deriv_u=deriv_u+xyzzyaaad86*xyzzyaaoc1*xyzzyaaah86
if(sd)sderiv_u=sderiv_u+xyzzyaaoc1*xyzzyaaah86*xyzzyaaag86*(-xyzzyaaad&
&86*(xyzzyaaac86+xyzzyaaab86)*xyzzyaaoc1*xyzzyaaag86+2.d0*xyzzyaaob1*x&
&yzzyaaad86)
endif
endif
endif
endif
end subroutine xyzzyaapr1
subroutine xyzzyaaps1(rij,s,val,fd,sd,value_u,deriv_u,sderiv_u)
implicit none
integer,intent(in) :: s
real(dp),intent(in) :: rij
real(dp),intent(out) :: value_u,deriv_u,sderiv_u
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa87
real(dp) xyzzyaaab87,xyzzyaaac87,xyzzyaaad87,xyzzyaaae87,xyzzyaaaf87,x&
&yzzyaaag87,xyzzyaaah87,xyzzyaaai87,xyzzyaaaj87,xyzzyaaak87
xyzzyaaah87=rij-xyzzyaach1
if(xyzzyaagx1==2)then
if(sd)then
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaaf87=xyzzyaaaj1(1,s)
xyzzyaaag87=0.d0
xyzzyaaac87=1.d0
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaad87=xyzzyaaac87*xyzzyaaaa87
xyzzyaaac87=xyzzyaaab87*xyzzyaaaa87
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
xyzzyaaaf87=xyzzyaaaf87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaac87
xyzzyaaag87=xyzzyaaag87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaad87
enddo
sderiv_u=2*xyzzyaaae87+xyzzyaaah87*(4*xyzzyaaaf87+xyzzyaaah87*xyzzyaaa&
&g87)
if(fd)deriv_u=xyzzyaaah87*(2*xyzzyaaae87+xyzzyaaah87*xyzzyaaaf87)
if(val)value_u=xyzzyaaae87*xyzzyaaah87*xyzzyaaah87
elseif(fd)then
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaaf87=xyzzyaaaj1(1,s)
xyzzyaaac87=1.d0
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaac87=xyzzyaaab87*xyzzyaaaa87
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
xyzzyaaaf87=xyzzyaaaf87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaac87
enddo
deriv_u=xyzzyaaah87*(2*xyzzyaaae87+xyzzyaaah87*xyzzyaaaf87)
if(val)value_u=xyzzyaaae87*xyzzyaaah87*xyzzyaaah87
else
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
enddo
value_u=xyzzyaaae87*xyzzyaaah87*xyzzyaaah87
endif
elseif(xyzzyaagx1==3)then
xyzzyaaai87=xyzzyaaah87*xyzzyaaah87
if(sd)then
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaaf87=xyzzyaaaj1(1,s)
xyzzyaaag87=0.d0
xyzzyaaac87=1.d0
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaad87=xyzzyaaac87*xyzzyaaaa87
xyzzyaaac87=xyzzyaaab87*xyzzyaaaa87
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
xyzzyaaaf87=xyzzyaaaf87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaac87
xyzzyaaag87=xyzzyaaag87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaad87
enddo
sderiv_u=xyzzyaaah87*(6*xyzzyaaae87+xyzzyaaah87*(6*xyzzyaaaf87+xyzzyaa&
&ah87*xyzzyaaag87))
if(fd)deriv_u=xyzzyaaai87*(3*xyzzyaaae87+(rij-xyzzyaach1)*xyzzyaaaf87)
if(val)value_u=xyzzyaaae87*xyzzyaaai87*xyzzyaaah87
elseif(fd)then
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaaf87=xyzzyaaaj1(1,s)
xyzzyaaac87=1.d0
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaac87=xyzzyaaab87*xyzzyaaaa87
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
xyzzyaaaf87=xyzzyaaaf87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaac87
enddo
deriv_u=xyzzyaaai87*(3*xyzzyaaae87+xyzzyaaah87*xyzzyaaaf87)
if(val)value_u=xyzzyaaae87*xyzzyaaai87*xyzzyaaah87
else
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
enddo
value_u=xyzzyaaae87*xyzzyaaai87*xyzzyaaah87
endif
else
if(sd)then
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaaf87=xyzzyaaaj1(1,s)
xyzzyaaag87=0.d0
xyzzyaaac87=1.d0
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaad87=xyzzyaaac87*xyzzyaaaa87
xyzzyaaac87=xyzzyaaab87*xyzzyaaaa87
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
xyzzyaaaf87=xyzzyaaaf87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaac87
xyzzyaaag87=xyzzyaaag87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaad87
enddo
if(xyzzyaagx1>1)then
xyzzyaaak87=xyzzyaaah87**(xyzzyaagx1-2)
xyzzyaaaj87=xyzzyaaak87*xyzzyaaah87
sderiv_u=xyzzyaaak87*(xyzzyaagx1*(xyzzyaagx1-1)*xyzzyaaae87+xyzzyaaah8&
&7*(2*xyzzyaagx1*xyzzyaaaf87+xyzzyaaah87*xyzzyaaag87))
if(fd)deriv_u=xyzzyaaaj87*(xyzzyaagx1*xyzzyaaae87+xyzzyaaah87*xyzzyaaa&
&f87)
if(val)value_u=xyzzyaaae87*xyzzyaaaj87*xyzzyaaah87
elseif(xyzzyaagx1==1)then
sderiv_u=2*xyzzyaaaf87+xyzzyaaah87*xyzzyaaag87
if(fd)deriv_u=xyzzyaaae87+xyzzyaaah87*xyzzyaaaf87
if(val)value_u=xyzzyaaae87*xyzzyaaah87
else
sderiv_u=xyzzyaaag87
deriv_u=xyzzyaaaf87
value_u=xyzzyaaae87
endif
elseif(fd)then
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaaf87=xyzzyaaaj1(1,s)
xyzzyaaac87=1.d0
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaac87=xyzzyaaab87*xyzzyaaaa87
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
xyzzyaaaf87=xyzzyaaaf87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaac87
enddo
if(xyzzyaagx1>0)then
xyzzyaaaj87=xyzzyaaah87**(xyzzyaagx1-1)
deriv_u=xyzzyaaaj87*(xyzzyaagx1*xyzzyaaae87+xyzzyaaah87*xyzzyaaaf87)
if(val)value_u=xyzzyaaae87*xyzzyaaaj87*xyzzyaaah87
else
deriv_u=xyzzyaaaf87
value_u=xyzzyaaae87
endif
else
xyzzyaaae87=xyzzyaaaj1(0,s)+rij*xyzzyaaaj1(1,s)
xyzzyaaab87=rij
do xyzzyaaaa87=2,xyzzyaabx1
xyzzyaaab87=xyzzyaaab87*rij
xyzzyaaae87=xyzzyaaae87+xyzzyaaaj1(xyzzyaaaa87,s)*xyzzyaaab87
enddo
value_u=xyzzyaaae87*xyzzyaaah87**xyzzyaagx1
endif
endif
end subroutine xyzzyaaps1
subroutine xyzzyaapt1(rij,s,val,fd,sd,value_u,deriv_u,sderiv_u)
implicit none
integer,intent(in) :: s
real(dp),intent(in) :: rij
real(dp),intent(out) :: value_u,deriv_u,sderiv_u
logical,intent(in) :: val,fd,sd
real(dp) xyzzyaaaa88,xyzzyaaab88,xyzzyaaac88,xyzzyaaad88,xyzzyaaae88,x&
&yzzyaaaf88,xyzzyaaag88,xyzzyaaah88,xyzzyaaai88,xyzzyaaaj88,xyzzyaaak8&
&8,xyzzyaaal88,xyzzyaaam88,xyzzyaaan88,xyzzyaaao88,xyzzyaaap88,xyzzyaa&
&aq88,xyzzyaaar88,xyzzyaaas88
logical xyzzyaaat88
xyzzyaaat88=fd.or.sd
if(rij/=0.d0)xyzzyaaak88=1.d0/rij
xyzzyaaal88=1.d0/xyzzyaach1
xyzzyaaam88=1.d0/xyzzyaaao1(s)
xyzzyaaaq88=rij*xyzzyaaam88
if(dimensionality==3)then
if(rij/=0.d0)then
xyzzyaaaa88=xyzzyaaan1(s)*xyzzyaaak88
xyzzyaaad88=exp(-xyzzyaaaq88)
xyzzyaaae88=-xyzzyaaaa88*(1.d0-xyzzyaaad88)
if(xyzzyaaat88)then
xyzzyaaas88=1.d0+xyzzyaaaq88
xyzzyaaac88=xyzzyaaaa88*xyzzyaaak88
xyzzyaaaf88=xyzzyaaac88*(1.d0-xyzzyaaas88*xyzzyaaad88)
endif
if(sd)xyzzyaaag88=-xyzzyaaac88*xyzzyaaak88*(2.d0-(1.d0+xyzzyaaas88*xyz&
&zyaaas88)*xyzzyaaad88)
else
xyzzyaaae88=-xyzzyaaan1(s)*xyzzyaaam88
if(xyzzyaaat88)xyzzyaaaf88=0.5d0*xyzzyaaan1(s)*xyzzyaaam88**2
if(sd)xyzzyaaag88=-xyzzyaaan1(s)*xyzzyaaam88**3*third
endif
elseif(dimensionality==2)then
if(rij/=0.d0)then
xyzzyaaan88=sqrt(xyzzyaaak88)
xyzzyaaar88=sqrt(xyzzyaaaq88)
xyzzyaaab88=xyzzyaaan1(s)*xyzzyaaan88
xyzzyaaad88=exp(-0.5d0*xyzzyaaaq88-xyzzyaaar88)
xyzzyaaae88=-xyzzyaaab88*(1.d0-xyzzyaaad88)
if(xyzzyaaat88)xyzzyaaaf88=0.5d0*xyzzyaaab88*xyzzyaaak88*(1.d0-xyzzyaa&
&ad88*(1.d0+xyzzyaaaq88+xyzzyaaar88))
if(sd)xyzzyaaag88=-0.25d0*xyzzyaaab88*xyzzyaaak88*xyzzyaaak88*(3.d0-(3&
&.d0+3*(xyzzyaaar88+xyzzyaaaq88)+(2*xyzzyaaar88+xyzzyaaaq88)*xyzzyaaaq&
&88)*xyzzyaaad88)
else
xyzzyaaae88=-xyzzyaaan1(s)*sqrt(xyzzyaaam88)
if(xyzzyaaat88)xyzzyaaaf88=xyzzyaaan1(s)*sqrt(xyzzyaaam88**3)*third
if(sd)xyzzyaaag88=0.d0
endif
else
call errstop('COMPUTE_U_RPA','RPA Jastrow not coded for 1D systems.')
endif
select case(xyzzyaagx1)
case(0)
xyzzyaaah88=1.d0
xyzzyaaai88=0.d0
xyzzyaaaj88=0.d0
case(1)
xyzzyaaao88=xyzzyaaal88*(xyzzyaach1-rij)
xyzzyaaah88=xyzzyaaao88
if(xyzzyaaat88)xyzzyaaai88=-xyzzyaaal88
if(sd)xyzzyaaaj88=0.d0
case(2)
xyzzyaaao88=xyzzyaaal88*(xyzzyaach1-rij)
xyzzyaaah88=xyzzyaaao88**2
if(xyzzyaaat88)xyzzyaaai88=-2*xyzzyaaal88*xyzzyaaao88
if(sd)xyzzyaaaj88=2*xyzzyaaal88**2
case(3)
xyzzyaaao88=xyzzyaaal88*(xyzzyaach1-rij)
xyzzyaaap88=xyzzyaaao88*xyzzyaaao88
xyzzyaaah88=xyzzyaaao88*xyzzyaaap88
if(xyzzyaaat88)xyzzyaaai88=-3*xyzzyaaal88*xyzzyaaap88
if(sd)xyzzyaaaj88=6*xyzzyaaao88*xyzzyaaal88**2
case default
xyzzyaaao88=xyzzyaaal88*(xyzzyaach1-rij)
xyzzyaaap88=xyzzyaaao88**(xyzzyaagx1-2)
xyzzyaaah88=xyzzyaaap88*xyzzyaaao88*xyzzyaaao88
if(xyzzyaaat88)xyzzyaaai88=-xyzzyaagx1*xyzzyaaal88*xyzzyaaap88*xyzzyaa&
&ao88
if(sd)xyzzyaaaj88=xyzzyaagx1*(xyzzyaagx1-1)*xyzzyaaap88*xyzzyaaal88*xy&
&zzyaaal88
end select
if(val)value_u=xyzzyaaae88*xyzzyaaah88
if(fd)deriv_u=xyzzyaaaf88*xyzzyaaah88+xyzzyaaae88*xyzzyaaai88
if(sd)sderiv_u=xyzzyaaag88*xyzzyaaah88+2*xyzzyaaaf88*xyzzyaaai88+xyzzy&
&aaae88*xyzzyaaaj88
end subroutine xyzzyaapt1
subroutine xyzzyaapu1(rhoij,modzij,ispin,jspin,val,fd,sd,value_ucyl,de&
&riv_ucyl_rho,deriv_ucyl_z,sderiv_ucyl_rho,sderiv_ucyl_z)
use slaarnabt,only : ddot
implicit none
integer,intent(in) :: ispin,jspin
real(dp),intent(in) :: rhoij,modzij
real(dp),intent(out) :: value_ucyl,deriv_ucyl_rho,deriv_ucyl_z,sderiv_&
&ucyl_rho,sderiv_ucyl_z
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa89,xyzzyaaab89,xyzzyaaac89
real(dp) xyzzyaaad89,xyzzyaaae89,xyzzyaaaf89,xyzzyaaag89,xyzzyaaah89,x&
&yzzyaaai89,xyzzyaaaj89,xyzzyaaak89,xyzzyaaal89,xyzzyaaam89(xyzzyaaby1&
&),xyzzyaaan89(xyzzyaaby1),xyzzyaaao89(2:xyzzyaaby1),xyzzyaaap89,xyzzy&
&aaaq89,xyzzyaaar89
value_ucyl=0.d0
deriv_ucyl_rho=0.d0
deriv_ucyl_z=0.d0
sderiv_ucyl_rho=0.d0
sderiv_ucyl_z=0.d0
if((rhoij<xyzzyaaci1.and.modzij<xyzzyaack1).or.xyzzyaagx1==0)then
xyzzyaaaa89=which_spair(ispin,jspin,levels_spairs)
if(sd)then
xyzzyaaam89(1)=rhoij
xyzzyaaan89(1)=1.d0
do xyzzyaaab89=2,xyzzyaaby1
xyzzyaaao89(xyzzyaaab89)=xyzzyaaab89*xyzzyaaan89(xyzzyaaab89-1)
xyzzyaaan89(xyzzyaaab89)=xyzzyaaab89*xyzzyaaam89(xyzzyaaab89-1)
xyzzyaaam89(xyzzyaaab89)=xyzzyaaam89(xyzzyaaab89-1)*rhoij
enddo
value_ucyl=value_ucyl+(xyzzyaaak1(0,0,xyzzyaaaa89)+ddot(xyzzyaaby1,xyz&
&zyaaak1(1,0,xyzzyaaaa89),1,xyzzyaaam89(1),1))+(xyzzyaaak1(0,1,xyzzyaa&
&aa89)+ddot(xyzzyaaby1,xyzzyaaak1(1,1,xyzzyaaaa89),1,xyzzyaaam89(1),1)&
&)*modzij
deriv_ucyl_rho=deriv_ucyl_rho+ddot(xyzzyaaby1,xyzzyaaak1(1,0,xyzzyaaaa&
&89),1,xyzzyaaan89(1),1)+ddot(xyzzyaaby1,xyzzyaaak1(1,1,xyzzyaaaa89),1&
&,xyzzyaaan89(1),1)*modzij
deriv_ucyl_z=deriv_ucyl_z+(xyzzyaaak1(0,1,xyzzyaaaa89)+ddot(xyzzyaaby1&
&,xyzzyaaak1(1,1,xyzzyaaaa89),1,xyzzyaaam89(1),1))
sderiv_ucyl_rho=sderiv_ucyl_rho+ddot(xyzzyaaby1-1,xyzzyaaak1(2,0,xyzzy&
&aaaa89),1,xyzzyaaao89(2),1)+ddot(xyzzyaaby1-1,xyzzyaaak1(2,1,xyzzyaaa&
&a89),1,xyzzyaaao89(2),1)*modzij
xyzzyaaap89=modzij
xyzzyaaaq89=1.d0
do xyzzyaaac89=2,xyzzyaabz1
xyzzyaaar89=xyzzyaaac89*xyzzyaaaq89
xyzzyaaaq89=xyzzyaaac89*xyzzyaaap89
xyzzyaaap89=xyzzyaaap89*modzij
value_ucyl=value_ucyl+(xyzzyaaak1(0,xyzzyaaac89,xyzzyaaaa89)+ddot(xyzz&
&yaaby1,xyzzyaaak1(1,xyzzyaaac89,xyzzyaaaa89),1,xyzzyaaam89(1),1))*xyz&
&zyaaap89
deriv_ucyl_rho=deriv_ucyl_rho+ddot(xyzzyaaby1,xyzzyaaak1(1,xyzzyaaac89&
&,xyzzyaaaa89),1,xyzzyaaan89(1),1)*xyzzyaaap89
deriv_ucyl_z=deriv_ucyl_z+(xyzzyaaak1(0,xyzzyaaac89,xyzzyaaaa89)+ddot(&
&xyzzyaaby1,xyzzyaaak1(1,xyzzyaaac89,xyzzyaaaa89),1,xyzzyaaam89(1),1))&
&*xyzzyaaaq89
sderiv_ucyl_rho=sderiv_ucyl_rho+ddot(xyzzyaaby1-1,xyzzyaaak1(2,xyzzyaa&
&ac89,xyzzyaaaa89),1,xyzzyaaao89(2),1)*xyzzyaaap89
sderiv_ucyl_z=sderiv_ucyl_z+(xyzzyaaak1(0,xyzzyaaac89,xyzzyaaaa89)+ddo&
&t(xyzzyaaby1,xyzzyaaak1(1,xyzzyaaac89,xyzzyaaaa89),1,xyzzyaaam89(1),1&
&))*xyzzyaaar89
enddo
elseif(fd)then
xyzzyaaam89(1)=rhoij
xyzzyaaan89(1)=1.d0
do xyzzyaaab89=2,xyzzyaaby1
xyzzyaaan89(xyzzyaaab89)=xyzzyaaab89*xyzzyaaam89(xyzzyaaab89-1)
xyzzyaaam89(xyzzyaaab89)=xyzzyaaam89(xyzzyaaab89-1)*rhoij
enddo
value_ucyl=value_ucyl+(xyzzyaaak1(0,0,xyzzyaaaa89)+ddot(xyzzyaaby1,xyz&
&zyaaak1(1,0,xyzzyaaaa89),1,xyzzyaaam89(1),1))+(xyzzyaaak1(0,1,xyzzyaa&
&aa89)+ddot(xyzzyaaby1,xyzzyaaak1(1,1,xyzzyaaaa89),1,xyzzyaaam89(1),1)&
&)*modzij
deriv_ucyl_rho=deriv_ucyl_rho+ddot(xyzzyaaby1,xyzzyaaak1(1,0,xyzzyaaaa&
&89),1,xyzzyaaan89(1),1)+ddot(xyzzyaaby1,xyzzyaaak1(1,1,xyzzyaaaa89),1&
&,xyzzyaaan89(1),1)*modzij
deriv_ucyl_z=deriv_ucyl_z+(xyzzyaaak1(0,1,xyzzyaaaa89)+ddot(xyzzyaaby1&
&,xyzzyaaak1(1,1,xyzzyaaaa89),1,xyzzyaaam89(1),1))
xyzzyaaap89=modzij
xyzzyaaaq89=1.d0
do xyzzyaaac89=2,xyzzyaabz1
xyzzyaaaq89=xyzzyaaac89*xyzzyaaap89
xyzzyaaap89=xyzzyaaap89*modzij
value_ucyl=value_ucyl+(xyzzyaaak1(0,xyzzyaaac89,xyzzyaaaa89)+ddot(xyzz&
&yaaby1,xyzzyaaak1(1,xyzzyaaac89,xyzzyaaaa89),1,xyzzyaaam89(1),1))*xyz&
&zyaaap89
deriv_ucyl_rho=deriv_ucyl_rho+ddot(xyzzyaaby1,xyzzyaaak1(1,xyzzyaaac89&
&,xyzzyaaaa89),1,xyzzyaaan89(1),1)*xyzzyaaap89
deriv_ucyl_z=deriv_ucyl_z+(xyzzyaaak1(0,xyzzyaaac89,xyzzyaaaa89)+ddot(&
&xyzzyaaby1,xyzzyaaak1(1,xyzzyaaac89,xyzzyaaaa89),1,xyzzyaaam89(1),1))&
&*xyzzyaaaq89
enddo
else
xyzzyaaam89(1)=rhoij
do xyzzyaaab89=2,xyzzyaaby1
xyzzyaaam89(xyzzyaaab89)=xyzzyaaam89(xyzzyaaab89-1)*rhoij
enddo
value_ucyl=value_ucyl+(xyzzyaaak1(0,0,xyzzyaaaa89)+ddot(xyzzyaaby1,xyz&
&zyaaak1(1,0,xyzzyaaaa89),1,xyzzyaaam89(1),1))+(xyzzyaaak1(0,1,xyzzyaa&
&aa89)+ddot(xyzzyaaby1,xyzzyaaak1(1,1,xyzzyaaaa89),1,xyzzyaaam89(1),1)&
&)*modzij
xyzzyaaap89=modzij
do xyzzyaaac89=2,xyzzyaabz1
xyzzyaaap89=xyzzyaaap89*modzij
value_ucyl=value_ucyl+(xyzzyaaak1(0,xyzzyaaac89,xyzzyaaaa89)+ddot(xyzz&
&yaaby1,xyzzyaaak1(1,xyzzyaaac89,xyzzyaaaa89),1,xyzzyaaam89(1),1))*xyz&
&zyaaap89
enddo
endif
xyzzyaaad89=rhoij-xyzzyaaci1
xyzzyaaae89=modzij-xyzzyaack1
if(xyzzyaagx1==2)then
xyzzyaaah89=xyzzyaaad89*xyzzyaaad89
xyzzyaaak89=xyzzyaaae89*xyzzyaaae89
xyzzyaaal89=xyzzyaaah89*xyzzyaaak89
if(sd)then
sderiv_ucyl_rho=sderiv_ucyl_rho*xyzzyaaal89+4*deriv_ucyl_rho*xyzzyaaad&
&89*xyzzyaaak89+2*value_ucyl*xyzzyaaak89
sderiv_ucyl_z=sderiv_ucyl_z*xyzzyaaal89+4*deriv_ucyl_z*xyzzyaaah89*xyz&
&zyaaae89+2*value_ucyl*xyzzyaaah89
endif
if(fd)then
deriv_ucyl_rho=deriv_ucyl_rho*xyzzyaaal89+2*value_ucyl*xyzzyaaad89*xyz&
&zyaaak89
deriv_ucyl_z=deriv_ucyl_z*xyzzyaaal89+2*value_ucyl*xyzzyaaah89*xyzzyaa&
&ae89
endif
if(val)value_ucyl=value_ucyl*xyzzyaaal89
elseif(xyzzyaagx1==3)then
xyzzyaaag89=xyzzyaaad89*xyzzyaaad89
xyzzyaaah89=xyzzyaaag89*xyzzyaaad89
xyzzyaaaj89=xyzzyaaae89*xyzzyaaae89
xyzzyaaak89=xyzzyaaaj89*xyzzyaaae89
xyzzyaaal89=xyzzyaaah89*xyzzyaaak89
if(sd)then
sderiv_ucyl_rho=sderiv_ucyl_rho*xyzzyaaal89+6*(deriv_ucyl_rho*xyzzyaaa&
&g89*xyzzyaaak89+value_ucyl*xyzzyaaad89*xyzzyaaak89)
sderiv_ucyl_z=sderiv_ucyl_z*xyzzyaaal89+6*(deriv_ucyl_z*xyzzyaaah89*xy&
&zzyaaaj89+value_ucyl*xyzzyaaah89*xyzzyaaae89)
endif
if(fd)then
deriv_ucyl_rho=deriv_ucyl_rho*xyzzyaaal89+3*value_ucyl*xyzzyaaag89*xyz&
&zyaaak89
deriv_ucyl_z=deriv_ucyl_z*xyzzyaaal89+3*value_ucyl*xyzzyaaah89*xyzzyaa&
&aj89
endif
if(val)value_ucyl=value_ucyl*xyzzyaaal89
else
xyzzyaaaf89=xyzzyaaad89**(xyzzyaagx1-2)
xyzzyaaag89=xyzzyaaaf89*xyzzyaaad89
xyzzyaaah89=xyzzyaaag89*xyzzyaaad89
xyzzyaaai89=xyzzyaaae89**(xyzzyaagx1-2)
xyzzyaaaj89=xyzzyaaai89*xyzzyaaae89
xyzzyaaak89=xyzzyaaaj89*xyzzyaaae89
xyzzyaaal89=xyzzyaaah89*xyzzyaaak89
if(sd)then
sderiv_ucyl_rho=sderiv_ucyl_rho*xyzzyaaal89+deriv_ucyl_rho*(2*xyzzyaag&
&x1)*xyzzyaaag89*xyzzyaaak89+value_ucyl*(xyzzyaagx1*(xyzzyaagx1-1))*xy&
&zzyaaaf89*xyzzyaaak89
sderiv_ucyl_z=sderiv_ucyl_z*xyzzyaaal89+deriv_ucyl_z*(2*xyzzyaagx1)*xy&
&zzyaaah89*xyzzyaaaj89+value_ucyl*(xyzzyaagx1*(xyzzyaagx1-1))*xyzzyaaa&
&h89*xyzzyaaai89
endif
if(fd)then
deriv_ucyl_rho=deriv_ucyl_rho*xyzzyaaal89+value_ucyl*xyzzyaagx1*xyzzya&
&aag89*xyzzyaaak89
deriv_ucyl_z=deriv_ucyl_z*xyzzyaaal89+value_ucyl*xyzzyaagx1*xyzzyaaah8&
&9*xyzzyaaaj89
endif
if(val)value_ucyl=value_ucyl*xyzzyaaal89
endif
endif
end subroutine xyzzyaapu1
subroutine xyzzyaapv1(rij,ispin,jspin,val,fd,sd,value_qcusp,deriv_qcus&
&p,sderiv_qcusp)
use slaarnaan,only : ee_kato_gamma
implicit none
integer,intent(in) :: ispin,jspin
real(dp),intent(in) :: rij
real(dp),intent(out) :: value_qcusp,deriv_qcusp,sderiv_qcusp
logical,intent(in) :: val,fd,sd
real(dp) xyzzyaaaa90,xyzzyaaab90,xyzzyaaac90,xyzzyaaad90,xyzzyaaae90,x&
&yzzyaaaf90,xyzzyaaag90,xyzzyaaah90,xyzzyaaai90,xyzzyaaaj90,xyzzyaaak9&
&0,xyzzyaaal90,xyzzyaaam90,xyzzyaaan90
character(8) f_zero
if(rij>xyzzyaacq1.or.heg_layer(ispin)==heg_layer(jspin))then
value_qcusp=0.d0
deriv_qcusp=0.d0
sderiv_qcusp=0.d0
else
xyzzyaaaa90=ee_kato_gamma(ispin,jspin,f_zero=f_zero)
xyzzyaaab90=(heg_ylayer(heg_layer(ispin))-heg_ylayer(heg_layer(jspin))&
&)**2+(heg_zlayer(heg_layer(ispin))-heg_zlayer(heg_layer(jspin)))**2
xyzzyaaac90=1.d0/xyzzyaacq1
xyzzyaaad90=rij*xyzzyaaac90
xyzzyaaae90=sqrt(xyzzyaacq1*xyzzyaacq1+xyzzyaaab90)
xyzzyaaai90=1.d0+xyzzyaaad90*xyzzyaaad90*(-6.d0+xyzzyaaad90*(8.d0-xyzz&
&yaaad90-xyzzyaaad90-xyzzyaaad90))
xyzzyaaaf90=sqrt(rij*rij+xyzzyaaab90)
xyzzyaaal90=xyzzyaaaa90*(xyzzyaaaf90-xyzzyaaae90)
if(val)value_qcusp=xyzzyaaal90*xyzzyaaai90
if(fd.or.sd)then
xyzzyaaag90=rij/xyzzyaaaf90
xyzzyaaam90=xyzzyaaaa90*xyzzyaaag90
xyzzyaaaj90=12.d0*xyzzyaaac90*xyzzyaaad90*(-1.d0+xyzzyaaad90*(2.d0-xyz&
&zyaaad90))
deriv_qcusp=xyzzyaaam90*xyzzyaaai90+xyzzyaaal90*xyzzyaaaj90
if(sd)then
xyzzyaaah90=(xyzzyaaaf90-rij*xyzzyaaag90)/(xyzzyaaaf90*xyzzyaaaf90)
xyzzyaaan90=xyzzyaaaa90*xyzzyaaah90
xyzzyaaak90=12.d0*xyzzyaaac90*xyzzyaaac90*(-1.d0+xyzzyaaad90*(4.d0-xyz&
&zyaaad90-xyzzyaaad90-xyzzyaaad90))
sderiv_qcusp=xyzzyaaan90*xyzzyaaai90+2.d0*xyzzyaaam90*xyzzyaaaj90+xyzz&
&yaaal90*xyzzyaaak90
endif
endif
endif
end subroutine xyzzyaapv1
subroutine xyzzyaapw1(id,eevecs,fd,sd)
use slaarnabt, only : ddot
implicit none
integer,intent(in) :: id
real(dp),intent(in) :: eevecs(4,netot,netot)
logical,intent(in) :: fd,sd
integer xyzzyaaaa91,xyzzyaaab91,xyzzyaaac91,xyzzyaaad91,xyzzyaaae91,xy&
&zzyaaaf91,xyzzyaaag91,xyzzyaaah91
real(dp) xyzzyaaai91,xyzzyaaaj91,xyzzyaaak91,xyzzyaaal91,xyzzyaaam91,x&
&yzzyaaan91(3),xyzzyaaao91(3),xyzzyaaap91,xyzzyaaaq91,xyzzyaaar91
logical xyzzyaaas91
logical,save :: xyzzyaaat91=.true.
integer,allocatable,save :: xyzzyaaau91(:),xyzzyaaav91(:)
if(xyzzyaaat91)then
allocate(xyzzyaaau91(netot),xyzzyaaav91(netot),stat=xyzzyaaah91)
call check_alloc(xyzzyaaah91,'COMPUTE_SG_BASIS','0')
xyzzyaaat91=.false.
endif
xyzzyaaas91=fd.or.sd
if(.not.xyzzyaaas91)then
xyzzyaaee1=0.d0
do xyzzyaaaa91=1,netot
do xyzzyaaab91=xyzzyaaaa91+1,netot
xyzzyaaae91=which_ee(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaac91=which_spin(xyzzyaaaa91)
xyzzyaaad91=which_spin(xyzzyaaab91)
xyzzyaaai91=eevecs(4,xyzzyaaab91,xyzzyaaaa91)
call xyzzyaapx1(xyzzyaaai91,xyzzyaaac91,xyzzyaaad91,.true.,.false.,.fa&
&lse.,xyzzyaaap91,xyzzyaaaq91,xyzzyaaar91,xyzzyaaif1(xyzzyaaae91))
xyzzyaaea1(1:dimensionality,xyzzyaaae91)=xyzzyaaap91*eevecs(1:dimensio&
&nality,xyzzyaaab91,xyzzyaaaa91)
xyzzyaaem1(xyzzyaaae91)=xyzzyaaap91*xyzzyaaai91
xyzzyaaak91=xyzzyaaem1(xyzzyaaae91)*xyzzyaaem1(xyzzyaaae91)
xyzzyaaee1(xyzzyaaaa91)=xyzzyaaee1(xyzzyaaaa91)+xyzzyaaak91
xyzzyaaee1(xyzzyaaab91)=xyzzyaaee1(xyzzyaaab91)+xyzzyaaak91
enddo
enddo
elseif(.not.sd)then
xyzzyaaee1=0.d0
xyzzyaaef1=0.d0
do xyzzyaaaa91=1,netot
do xyzzyaaab91=xyzzyaaaa91+1,netot
xyzzyaaae91=which_ee(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaac91=which_spin(xyzzyaaaa91)
xyzzyaaad91=which_spin(xyzzyaaab91)
xyzzyaaai91=eevecs(4,xyzzyaaab91,xyzzyaaaa91)
xyzzyaaaj91=1.d0/xyzzyaaai91
xyzzyaaao91(1:dimensionality)=eevecs(1:dimensionality,xyzzyaaab91,xyzz&
&yaaaa91)
call xyzzyaapx1(xyzzyaaai91,xyzzyaaac91,xyzzyaaad91,.true.,.true.,.fal&
&se.,xyzzyaaap91,xyzzyaaaq91,xyzzyaaar91,xyzzyaaif1(xyzzyaaae91))
xyzzyaaea1(1:dimensionality,xyzzyaaae91)=xyzzyaaap91*eevecs(1:dimensio&
&nality,xyzzyaaab91,xyzzyaaaa91)
xyzzyaaem1(xyzzyaaae91)=xyzzyaaap91*xyzzyaaai91
xyzzyaaak91=xyzzyaaem1(xyzzyaaae91)*xyzzyaaem1(xyzzyaaae91)
xyzzyaaee1(xyzzyaaaa91)=xyzzyaaee1(xyzzyaaaa91)+xyzzyaaak91
xyzzyaaee1(xyzzyaaab91)=xyzzyaaee1(xyzzyaaab91)+xyzzyaaak91
do xyzzyaaaf91=1,dimensionality
xyzzyaaak91=xyzzyaaaq91*xyzzyaaao91(xyzzyaaaf91)
xyzzyaaal91=xyzzyaaaj91*xyzzyaaao91(xyzzyaaaf91)
xyzzyaaam91=xyzzyaaak91*xyzzyaaal91+xyzzyaaap91
xyzzyaaeb1(xyzzyaaaf91,xyzzyaaaf91,xyzzyaaae91)=xyzzyaaam91
xyzzyaaen1(xyzzyaaaf91,xyzzyaaae91)=xyzzyaaak91+xyzzyaaap91*xyzzyaaal9&
&1
xyzzyaaef1(xyzzyaaaf91,xyzzyaaaf91,xyzzyaaaa91)=xyzzyaaef1(xyzzyaaaf91&
&,xyzzyaaaf91,xyzzyaaaa91)+xyzzyaaam91
xyzzyaaef1(xyzzyaaaf91,xyzzyaaaf91,xyzzyaaab91)=xyzzyaaef1(xyzzyaaaf91&
&,xyzzyaaaf91,xyzzyaaab91)+xyzzyaaam91
do xyzzyaaag91=xyzzyaaaf91+1,dimensionality
xyzzyaaam91=xyzzyaaak91*xyzzyaaao91(xyzzyaaag91)*xyzzyaaaj91
xyzzyaaeb1(xyzzyaaag91,xyzzyaaaf91,xyzzyaaae91)=xyzzyaaam91
xyzzyaaeb1(xyzzyaaaf91,xyzzyaaag91,xyzzyaaae91)=xyzzyaaam91
xyzzyaaef1(xyzzyaaaf91,xyzzyaaag91,xyzzyaaaa91)=xyzzyaaef1(xyzzyaaaf91&
&,xyzzyaaag91,xyzzyaaaa91)+xyzzyaaam91
xyzzyaaef1(xyzzyaaag91,xyzzyaaaf91,xyzzyaaaa91)=xyzzyaaef1(xyzzyaaag91&
&,xyzzyaaaf91,xyzzyaaaa91)+xyzzyaaam91
xyzzyaaef1(xyzzyaaaf91,xyzzyaaag91,xyzzyaaab91)=xyzzyaaef1(xyzzyaaaf91&
&,xyzzyaaag91,xyzzyaaab91)+xyzzyaaam91
xyzzyaaef1(xyzzyaaag91,xyzzyaaaf91,xyzzyaaab91)=xyzzyaaef1(xyzzyaaag91&
&,xyzzyaaaf91,xyzzyaaab91)+xyzzyaaam91
enddo
enddo
enddo
enddo
else
xyzzyaaee1=0.d0
xyzzyaaef1=0.d0
xyzzyaaeg1=0.d0
xyzzyaaeh1=0.d0
xyzzyaaei1=0.d0
do xyzzyaaaa91=1,netot
do xyzzyaaab91=xyzzyaaaa91+1,netot
xyzzyaaae91=which_ee(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaac91=which_spin(xyzzyaaaa91)
xyzzyaaad91=which_spin(xyzzyaaab91)
xyzzyaaai91=eevecs(4,xyzzyaaab91,xyzzyaaaa91)
xyzzyaaaj91=1.d0/xyzzyaaai91
xyzzyaaao91(1:dimensionality)=eevecs(1:dimensionality,xyzzyaaab91,xyzz&
&yaaaa91)
call xyzzyaapx1(xyzzyaaai91,xyzzyaaac91,xyzzyaaad91,.true.,.true.,.tru&
&e.,xyzzyaaap91,xyzzyaaaq91,xyzzyaaar91,xyzzyaaif1(xyzzyaaae91))
xyzzyaaea1(1:dimensionality,xyzzyaaae91)=xyzzyaaap91*eevecs(1:dimensio&
&nality,xyzzyaaab91,xyzzyaaaa91)
xyzzyaaem1(xyzzyaaae91)=xyzzyaaap91*xyzzyaaai91
xyzzyaaak91=xyzzyaaem1(xyzzyaaae91)*xyzzyaaem1(xyzzyaaae91)
xyzzyaaee1(xyzzyaaaa91)=xyzzyaaee1(xyzzyaaaa91)+xyzzyaaak91
xyzzyaaee1(xyzzyaaab91)=xyzzyaaee1(xyzzyaaab91)+xyzzyaaak91
do xyzzyaaaf91=1,dimensionality
xyzzyaaak91=xyzzyaaaq91*xyzzyaaao91(xyzzyaaaf91)
xyzzyaaal91=xyzzyaaaj91*xyzzyaaao91(xyzzyaaaf91)
xyzzyaaam91=xyzzyaaak91*xyzzyaaal91+xyzzyaaap91
xyzzyaaeb1(xyzzyaaaf91,xyzzyaaaf91,xyzzyaaae91)=xyzzyaaam91
xyzzyaaen1(xyzzyaaaf91,xyzzyaaae91)=xyzzyaaak91+xyzzyaaap91*xyzzyaaal9&
&1
xyzzyaaef1(xyzzyaaaf91,xyzzyaaaf91,xyzzyaaaa91)=xyzzyaaef1(xyzzyaaaf91&
&,xyzzyaaaf91,xyzzyaaaa91)+xyzzyaaam91
xyzzyaaef1(xyzzyaaaf91,xyzzyaaaf91,xyzzyaaab91)=xyzzyaaef1(xyzzyaaaf91&
&,xyzzyaaaf91,xyzzyaaab91)+xyzzyaaam91
do xyzzyaaag91=xyzzyaaaf91+1,dimensionality
xyzzyaaam91=xyzzyaaak91*xyzzyaaao91(xyzzyaaag91)*xyzzyaaaj91
xyzzyaaeb1(xyzzyaaag91,xyzzyaaaf91,xyzzyaaae91)=xyzzyaaam91
xyzzyaaeb1(xyzzyaaaf91,xyzzyaaag91,xyzzyaaae91)=xyzzyaaam91
xyzzyaaef1(xyzzyaaaf91,xyzzyaaag91,xyzzyaaaa91)=xyzzyaaef1(xyzzyaaaf91&
&,xyzzyaaag91,xyzzyaaaa91)+xyzzyaaam91
xyzzyaaef1(xyzzyaaag91,xyzzyaaaf91,xyzzyaaaa91)=xyzzyaaef1(xyzzyaaag91&
&,xyzzyaaaf91,xyzzyaaaa91)+xyzzyaaam91
xyzzyaaef1(xyzzyaaaf91,xyzzyaaag91,xyzzyaaab91)=xyzzyaaef1(xyzzyaaaf91&
&,xyzzyaaag91,xyzzyaaab91)+xyzzyaaam91
xyzzyaaef1(xyzzyaaag91,xyzzyaaaf91,xyzzyaaab91)=xyzzyaaef1(xyzzyaaag91&
&,xyzzyaaaf91,xyzzyaaab91)+xyzzyaaam91
enddo
enddo
xyzzyaaeg1(xyzzyaaae91)=sum(xyzzyaaeb1(1:dimensionality,1:dimensionali&
&ty,xyzzyaaae91)**2)
xyzzyaaak91=xyzzyaaar91+dble(dimensionality+1)*xyzzyaaaq91*xyzzyaaaj91
xyzzyaaan91(1:dimensionality)=xyzzyaaak91*xyzzyaaao91(1:dimensionality&
&)
xyzzyaaec1(1:dimensionality,xyzzyaaae91)=xyzzyaaan91(1:dimensionality)
xyzzyaaeo1(xyzzyaaae91)=xyzzyaaak91*xyzzyaaai91+dble(dimensionality-1)&
&*xyzzyaaap91*xyzzyaaaj91
xyzzyaaei1(1:dimensionality,xyzzyaaaa91)=xyzzyaaei1(1:dimensionality,x&
&yzzyaaaa91)+xyzzyaaan91(1:dimensionality)
xyzzyaaei1(1:dimensionality,xyzzyaaab91)=xyzzyaaei1(1:dimensionality,x&
&yzzyaaab91)-xyzzyaaan91(1:dimensionality)
enddo
xyzzyaaeh1(xyzzyaaaa91)=sum(xyzzyaaef1(1:dimensionality,1:dimensionali&
&ty,xyzzyaaaa91)**2)
enddo
endif
xyzzyaaau91=0
xyzzyaaav91=0
xyzzyaadz1=0.d0
xyzzyaael1=0.d0
do xyzzyaaaa91=1,netot
do xyzzyaaab91=xyzzyaaaa91+1,netot
xyzzyaaae91=which_ee(xyzzyaaaa91,xyzzyaaab91)
if(xyzzyaaif1(xyzzyaaae91))cycle
xyzzyaaau91(xyzzyaaaa91)=xyzzyaaau91(xyzzyaaaa91)+1
xyzzyaaav91(xyzzyaaaa91)=xyzzyaaab91
xyzzyaaau91(xyzzyaaab91)=xyzzyaaau91(xyzzyaaab91)+1
xyzzyaaav91(xyzzyaaab91)=xyzzyaaaa91
xyzzyaaan91=xyzzyaaea1(:,xyzzyaaae91)
xyzzyaadz1(:,xyzzyaaaa91)=xyzzyaadz1(:,xyzzyaaaa91)+xyzzyaaan91
xyzzyaadz1(:,xyzzyaaab91)=xyzzyaadz1(:,xyzzyaaab91)-xyzzyaaan91
enddo
xyzzyaael1(xyzzyaaaa91)=sqrt(ddot(3,xyzzyaadz1(1,xyzzyaaaa91),1,xyzzya&
&adz1(1,xyzzyaaaa91),1))
enddo
if(id==0)then
xyzzyaaed1=0.d0
do xyzzyaaaa91=1,netot
do xyzzyaaab91=xyzzyaaaa91+1,netot
xyzzyaaae91=which_ee(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaak91=ddot(3,xyzzyaaea1(1,xyzzyaaae91),1,xyzzyaadz1(1,xyzzyaaab9&
&1),1)
xyzzyaaed1(xyzzyaaaa91)=xyzzyaaed1(xyzzyaaaa91)+xyzzyaaak91
xyzzyaaak91=-ddot(3,xyzzyaaea1(1,xyzzyaaae91),1,xyzzyaadz1(1,xyzzyaaaa&
&91),1)
xyzzyaaed1(xyzzyaaab91)=xyzzyaaed1(xyzzyaaab91)+xyzzyaaak91
enddo
enddo
else
xyzzyaaed1(id)=0.d0
do xyzzyaaab91=1,id-1
xyzzyaaae91=which_ee(id,xyzzyaaab91)
xyzzyaaak91=-ddot(3,xyzzyaaea1(1,xyzzyaaae91),1,xyzzyaadz1(1,xyzzyaaab&
&91),1)
xyzzyaaed1(id)=xyzzyaaed1(id)+xyzzyaaak91
enddo
do xyzzyaaab91=id+1,netot
xyzzyaaae91=which_ee(id,xyzzyaaab91)
xyzzyaaak91=ddot(3,xyzzyaaea1(1,xyzzyaaae91),1,xyzzyaadz1(1,xyzzyaaab9&
&1),1)
xyzzyaaed1(id)=xyzzyaaed1(id)+xyzzyaaak91
enddo
endif
xyzzyaaie1=.false.
xyzzyaaid1=.false.
if(any(xyzzyaaau91>1))then
xyzzyaaid1=.true.
do xyzzyaaaa91=1,netot
if(xyzzyaaau91(xyzzyaaaa91)>1)then
xyzzyaaie1(xyzzyaaaa91)=.true.
elseif(xyzzyaaau91(xyzzyaaaa91)==1)then
xyzzyaaab91=xyzzyaaav91(xyzzyaaaa91)
if(xyzzyaaau91(xyzzyaaab91)>1)then
xyzzyaaie1(xyzzyaaaa91)=.true.
exit
endif
endif
enddo
endif
if(id/=0)then
if(.not.xyzzyaaie1(id))then
if(xyzzyaaas91)xyzzyaaej1(:,:,id)=0.d0
if(sd)xyzzyaaek1(:,id)=0.d0
return
endif
else
if(.not.xyzzyaaid1)then
if(xyzzyaaas91)xyzzyaaej1=0.d0
if(sd)xyzzyaaek1=0.d0
return
endif
endif
if(xyzzyaaas91)then
if(id<1)then
xyzzyaaej1=0.d0
do xyzzyaaaa91=1,netot
xyzzyaaej1(1,xyzzyaaaa91,xyzzyaaaa91)=ddot(3,xyzzyaadz1(1,xyzzyaaaa91)&
&,1,xyzzyaaef1(1,1,xyzzyaaaa91),1)
xyzzyaaej1(2,xyzzyaaaa91,xyzzyaaaa91)=ddot(3,xyzzyaadz1(1,xyzzyaaaa91)&
&,1,xyzzyaaef1(1,2,xyzzyaaaa91),1)
xyzzyaaej1(3,xyzzyaaaa91,xyzzyaaaa91)=ddot(3,xyzzyaadz1(1,xyzzyaaaa91)&
&,1,xyzzyaaef1(1,3,xyzzyaaaa91),1)
do xyzzyaaab91=xyzzyaaaa91+1,netot
xyzzyaaae91=which_ee(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaej1(1,xyzzyaaab91,xyzzyaaaa91)=-ddot(3,xyzzyaadz1(1,xyzzyaaab91&
&),1,xyzzyaaeb1(1,1,xyzzyaaae91),1)
xyzzyaaej1(2,xyzzyaaab91,xyzzyaaaa91)=-ddot(3,xyzzyaadz1(1,xyzzyaaab91&
&),1,xyzzyaaeb1(1,2,xyzzyaaae91),1)
xyzzyaaej1(3,xyzzyaaab91,xyzzyaaaa91)=-ddot(3,xyzzyaadz1(1,xyzzyaaab91&
&),1,xyzzyaaeb1(1,3,xyzzyaaae91),1)
xyzzyaaej1(1,xyzzyaaaa91,xyzzyaaab91)=-ddot(3,xyzzyaadz1(1,xyzzyaaaa91&
&),1,xyzzyaaeb1(1,1,xyzzyaaae91),1)
xyzzyaaej1(2,xyzzyaaaa91,xyzzyaaab91)=-ddot(3,xyzzyaadz1(1,xyzzyaaaa91&
&),1,xyzzyaaeb1(1,2,xyzzyaaae91),1)
xyzzyaaej1(3,xyzzyaaaa91,xyzzyaaab91)=-ddot(3,xyzzyaadz1(1,xyzzyaaaa91&
&),1,xyzzyaaeb1(1,3,xyzzyaaae91),1)
enddo
enddo
else
xyzzyaaej1(:,:,id)=0.d0
xyzzyaaej1(1,id,id)=ddot(3,xyzzyaadz1(1,id),1,xyzzyaaef1(1,1,id),1)
xyzzyaaej1(2,id,id)=ddot(3,xyzzyaadz1(1,id),1,xyzzyaaef1(1,2,id),1)
xyzzyaaej1(3,id,id)=ddot(3,xyzzyaadz1(1,id),1,xyzzyaaef1(1,3,id),1)
do xyzzyaaab91=1,netot
if(xyzzyaaab91==id)cycle
xyzzyaaae91=which_ee(id,xyzzyaaab91)
xyzzyaaej1(1,xyzzyaaab91,id)=-ddot(3,xyzzyaadz1(1,xyzzyaaab91),1,xyzzy&
&aaeb1(1,1,xyzzyaaae91),1)
xyzzyaaej1(2,xyzzyaaab91,id)=-ddot(3,xyzzyaadz1(1,xyzzyaaab91),1,xyzzy&
&aaeb1(1,2,xyzzyaaae91),1)
xyzzyaaej1(3,xyzzyaaab91,id)=-ddot(3,xyzzyaadz1(1,xyzzyaaab91),1,xyzzy&
&aaeb1(1,3,xyzzyaaae91),1)
enddo
endif
endif
if(sd)then
if(id<1)then
xyzzyaaek1=0.d0
do xyzzyaaaa91=1,netot
xyzzyaaek1(xyzzyaaaa91,xyzzyaaaa91)=ddot(3,xyzzyaadz1(1,xyzzyaaaa91),1&
&,xyzzyaaei1(1,xyzzyaaaa91),1)+xyzzyaaeh1(xyzzyaaaa91)
do xyzzyaaab91=xyzzyaaaa91+1,netot
xyzzyaaae91=which_ee(xyzzyaaaa91,xyzzyaaab91)
xyzzyaaek1(xyzzyaaab91,xyzzyaaaa91)=-ddot(3,xyzzyaadz1(1,xyzzyaaab91),&
&1,xyzzyaaec1(1,xyzzyaaae91),1)+xyzzyaaeg1(xyzzyaaae91)
xyzzyaaek1(xyzzyaaaa91,xyzzyaaab91)=ddot(3,xyzzyaadz1(1,xyzzyaaaa91),1&
&,xyzzyaaec1(1,xyzzyaaae91),1)+xyzzyaaeg1(xyzzyaaae91)
enddo
enddo
else
xyzzyaaek1(:,id)=0.d0
xyzzyaaek1(id,id)=ddot(3,xyzzyaadz1(1,id),1,xyzzyaaei1(1,id),1)+xyzzya&
&aeh1(id)
do xyzzyaaab91=1,id-1
xyzzyaaae91=which_ee(id,xyzzyaaab91)
xyzzyaaek1(xyzzyaaab91,id)=ddot(3,xyzzyaadz1(1,xyzzyaaab91),1,xyzzyaae&
&c1(1,xyzzyaaae91),1)+xyzzyaaeg1(xyzzyaaae91)
enddo
do xyzzyaaab91=id+1,netot
xyzzyaaae91=which_ee(id,xyzzyaaab91)
xyzzyaaek1(xyzzyaaab91,id)=-ddot(3,xyzzyaadz1(1,xyzzyaaab91),1,xyzzyaa&
&ec1(1,xyzzyaaae91),1)+xyzzyaaeg1(xyzzyaaae91)
enddo
endif
endif
end subroutine xyzzyaapw1
subroutine xyzzyaapx1(rij,ispin,jspin,val,fd,sd,value_w,deriv_w,sderiv&
&_w,beyond_cutoff)
implicit none
integer,intent(in) :: ispin,jspin
real(dp),intent(in) :: rij
real(dp),intent(out) :: value_w,deriv_w,sderiv_w
logical,intent(in) :: val,fd,sd
logical,intent(out) :: beyond_cutoff
integer xyzzyaaaa92,xyzzyaaab92
real(dp) xyzzyaaac92,xyzzyaaad92,xyzzyaaae92,xyzzyaaaf92,xyzzyaaag92,x&
&yzzyaaah92
if(rij>xyzzyaacl1.and.xyzzyaagx1>0)then
value_w=0.d0
deriv_w=0.d0
sderiv_w=0.d0
beyond_cutoff=.true.
return
endif
beyond_cutoff=.false.
xyzzyaaaa92=which_spair(ispin,jspin,levels_spairs)
if(xyzzyaagx1==2)then
if(.not.fd.and..not.sd)then
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
enddo
value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**2
elseif(.not.sd)then
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaag92=xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaad92=1.d0
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaad92=xyzzyaaac92*xyzzyaaab92
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
xyzzyaaag92=xyzzyaaag92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaad9&
&2
enddo
deriv_w=(rij-xyzzyaacl1)*(2*xyzzyaaaf92+(rij-xyzzyaacl1)*xyzzyaaag92)
if(val)value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**2
else
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaag92=xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaah92=0.d0
xyzzyaaad92=1.d0
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaae92=xyzzyaaad92*xyzzyaaab92
xyzzyaaad92=xyzzyaaac92*xyzzyaaab92
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
xyzzyaaag92=xyzzyaaag92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaad9&
&2
xyzzyaaah92=xyzzyaaah92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaae9&
&2
enddo
sderiv_w=2*xyzzyaaaf92+(rij-xyzzyaacl1)*(4*xyzzyaaag92+(rij-xyzzyaacl1&
&)*xyzzyaaah92)
if(fd)deriv_w=(rij-xyzzyaacl1)*(2*xyzzyaaaf92+(rij-xyzzyaacl1)*xyzzyaa&
&ag92)
if(val)value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**2
endif
elseif(xyzzyaagx1==3)then
if(.not.fd.and..not.sd)then
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
enddo
value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**3
elseif(.not.sd)then
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaag92=xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaad92=1.d0
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaad92=xyzzyaaac92*xyzzyaaab92
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
xyzzyaaag92=xyzzyaaag92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaad9&
&2
enddo
deriv_w=(rij-xyzzyaacl1)**2*(3*xyzzyaaaf92+(rij-xyzzyaacl1)*xyzzyaaag9&
&2)
if(val)value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**3
else
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaag92=xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaah92=0.d0
xyzzyaaad92=1.d0
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaae92=xyzzyaaad92*xyzzyaaab92
xyzzyaaad92=xyzzyaaac92*xyzzyaaab92
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
xyzzyaaag92=xyzzyaaag92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaad9&
&2
xyzzyaaah92=xyzzyaaah92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaae9&
&2
enddo
sderiv_w=(rij-xyzzyaacl1)*(6*xyzzyaaaf92+(rij-xyzzyaacl1)*(6*xyzzyaaag&
&92+(rij-xyzzyaacl1)*xyzzyaaah92))
if(fd)deriv_w=(rij-xyzzyaacl1)**2*(3*xyzzyaaaf92+(rij-xyzzyaacl1)*xyzz&
&yaaag92)
if(val)value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**3
endif
else
if(.not.fd.and..not.sd)then
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
enddo
value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**xyzzyaagx1
elseif(.not.sd)then
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaag92=xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaad92=1.d0
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaad92=xyzzyaaac92*xyzzyaaab92
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
xyzzyaaag92=xyzzyaaag92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaad9&
&2
enddo
if(xyzzyaagx1>0)then
deriv_w=xyzzyaagx1*(rij-xyzzyaacl1)**(xyzzyaagx1-1)*xyzzyaaaf92+(rij-x&
&yzzyaacl1)**xyzzyaagx1*xyzzyaaag92
else
deriv_w=xyzzyaaag92
endif
if(val)value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**xyzzyaagx1
else
xyzzyaaaf92=xyzzyaaaq1(0,xyzzyaaaa92)+rij*xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaag92=xyzzyaaaq1(1,xyzzyaaaa92)
xyzzyaaah92=0.d0
xyzzyaaad92=1.d0
xyzzyaaac92=rij
do xyzzyaaab92=2,xyzzyaacb1
xyzzyaaae92=xyzzyaaad92*xyzzyaaab92
xyzzyaaad92=xyzzyaaac92*xyzzyaaab92
xyzzyaaac92=xyzzyaaac92*rij
xyzzyaaaf92=xyzzyaaaf92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaac9&
&2
xyzzyaaag92=xyzzyaaag92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaad9&
&2
xyzzyaaah92=xyzzyaaah92+xyzzyaaaq1(xyzzyaaab92,xyzzyaaaa92)*xyzzyaaae9&
&2
enddo
if(xyzzyaagx1>1)then
sderiv_w=(rij-xyzzyaacl1)**(xyzzyaagx1-2)*(xyzzyaagx1*(xyzzyaagx1-1)*x&
&yzzyaaaf92+(rij-xyzzyaacl1)*(2*xyzzyaagx1*xyzzyaaag92+(rij-xyzzyaacl1&
&)*xyzzyaaah92))
elseif(xyzzyaagx1==1)then
sderiv_w=2*xyzzyaaag92+(rij-xyzzyaacl1)*xyzzyaaah92
else
sderiv_w=xyzzyaaah92
endif
if(fd)then
if(xyzzyaagx1>0)then
deriv_w=(rij-xyzzyaacl1)**(xyzzyaagx1-1)*(xyzzyaagx1*xyzzyaaaf92+(rij-&
&xyzzyaacl1)*xyzzyaaag92)
else
deriv_w=xyzzyaaag92
endif
endif
if(val)value_w=xyzzyaaaf92*(rij-xyzzyaacl1)**xyzzyaagx1
endif
endif
end subroutine xyzzyaapx1
subroutine xyzzyaapy1(ii,eevecs,val,fd,sd,value_h,grad_h,lap_h)
use slaarnabt, only : ddot
implicit none
integer,intent(in) :: ii
real(dp),intent(in) :: eevecs(4,netot,netot)
real(dp),intent(inout) :: value_h,grad_h(3),lap_h
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa93,xyzzyaaab93,xyzzyaaac93,xyzzyaaad93,xyzzyaaae93,xy&
&zzyaaaf93,xyzzyaaag93(3),xyzzyaaah93(3),xyzzyaaai93,xyzzyaaaj93,xyzzy&
&aaak93
real(dp) xyzzyaaal93,xyzzyaaam93,xyzzyaaan93,xyzzyaaao93(3),xyzzyaaap9&
&3(3),xyzzyaaaq93,xyzzyaaar93,xyzzyaaas93,xyzzyaaat93,xyzzyaaau93,xyzz&
&yaaav93,xyzzyaaaw93,xyzzyaaax93,xyzzyaaay93,xyzzyaaaz93,xyzzyaaba93,x&
&yzzyaabb93(3),xyzzyaabc93,xyzzyaabd93,xyzzyaabe93(3),xyzzyaabf93,xyzz&
&yaabg93,xyzzyaabh93,xyzzyaabi93,xyzzyaabj93,xyzzyaabk93,xyzzyaabl93,x&
&yzzyaabm93,xyzzyaabn93,xyzzyaabo93,xyzzyaabp93,xyzzyaabq93,xyzzyaabr9&
&3,xyzzyaabs93
logical xyzzyaabt93
value_h=0.d0
grad_h=0.d0
lap_h=0.d0
xyzzyaaab93=which_spin(ii)
xyzzyaabt93=fd.or.sd
if(sd)xyzzyaabm93=dble(dimensionality-1)
do xyzzyaaae93=1,netot
if(xyzzyaaae93==ii)cycle
xyzzyaaac93=which_spin(xyzzyaaae93)
do xyzzyaaaf93=xyzzyaaae93+1,netot
if(xyzzyaaaf93==ii)cycle
xyzzyaaal93=eevecs(4,ii,xyzzyaaae93)
xyzzyaaam93=eevecs(4,xyzzyaaae93,xyzzyaaaf93)
xyzzyaaan93=eevecs(4,ii,xyzzyaaaf93)
if(xyzzyaaal93>xyzzyaacm1.or.xyzzyaaan93>xyzzyaacm1.or.xyzzyaaam93>xyz&
&zyaacm1)cycle
xyzzyaaad93=which_spin(xyzzyaaaf93)
xyzzyaaag93(1:3)=xyzzyaabd1(1:3,xyzzyaaab93,xyzzyaaac93,xyzzyaaad93)
xyzzyaaaa93=which_striplet(xyzzyaaab93,xyzzyaaac93,xyzzyaaad93,xyzzyaa&
&fm1)
call xyzzyaapz1(xyzzyaaam93,xyzzyaaaw93,xyzzyaabk93,xyzzyaabl93,.false&
&.)
call xyzzyaapz1(xyzzyaaan93,xyzzyaaax93,xyzzyaaay93,xyzzyaaaz93,xyzzya&
&abt93)
call xyzzyaapz1(xyzzyaaal93,xyzzyaaat93,xyzzyaaau93,xyzzyaaav93,xyzzya&
&abt93)
xyzzyaabd93=xyzzyaaaw93*xyzzyaaax93*xyzzyaaat93
if(xyzzyaabt93)then
xyzzyaabn93=0.d0
if(xyzzyaaan93/=0.d0)xyzzyaabn93=1.d0/xyzzyaaan93
xyzzyaabo93=0.d0
if(xyzzyaaal93/=0.d0)xyzzyaabo93=1.d0/xyzzyaaal93
xyzzyaaap93=0.d0
xyzzyaaao93=0.d0
xyzzyaaap93(1:dimensionality)=-eevecs(1:dimensionality,ii,xyzzyaaaf93)&
&*xyzzyaabn93
xyzzyaaao93(1:dimensionality)=-eevecs(1:dimensionality,ii,xyzzyaaae93)&
&*xyzzyaabo93
xyzzyaabe93(:)=xyzzyaaaw93*(xyzzyaaau93*xyzzyaaax93*xyzzyaaao93(:)+xyz&
&zyaaat93*xyzzyaaay93*xyzzyaaap93(:))
if(sd)then
xyzzyaaar93=xyzzyaabm93*xyzzyaabn93
xyzzyaaaq93=xyzzyaabm93*xyzzyaabo93
xyzzyaaas93=ddot(dimensionality,xyzzyaaao93(1),1,xyzzyaaap93(1),1)
xyzzyaabf93=xyzzyaaaw93*(xyzzyaaax93*(xyzzyaaav93+xyzzyaaau93*xyzzyaaa&
&q93)+xyzzyaaat93*(xyzzyaaaz93+xyzzyaaay93*xyzzyaaar93)+2.d0*xyzzyaaau&
&93*xyzzyaaay93*xyzzyaaas93)
endif
endif
xyzzyaaba93=0.d0
xyzzyaabb93=0.d0
xyzzyaabc93=0.d0
xyzzyaabi93=1.d0
do xyzzyaaai93=0,xyzzyaacc1
xyzzyaabq93=0.d0
xyzzyaabp93=0.d0
xyzzyaabh93=1.d0
xyzzyaaah93(1)=xyzzyaaai93
do xyzzyaaaj93=0,xyzzyaacc1
xyzzyaabs93=0.d0
xyzzyaabr93=0.d0
xyzzyaabg93=1.d0
xyzzyaaah93(2)=xyzzyaaaj93
do xyzzyaaak93=0,xyzzyaacc1
xyzzyaaah93(3)=xyzzyaaak93
xyzzyaabj93=6.d0*xyzzyaaar1(xyzzyaaah93(xyzzyaaag93(1)),xyzzyaaah93(xy&
&zzyaaag93(2)),xyzzyaaah93(xyzzyaaag93(3)),xyzzyaaaa93)*xyzzyaabi93
xyzzyaaba93=xyzzyaaba93+xyzzyaabj93*xyzzyaabh93*xyzzyaabg93
if(xyzzyaabt93)then
xyzzyaabb93(:)=xyzzyaabb93(:)+xyzzyaabj93*(xyzzyaabp93*xyzzyaabg93*xyz&
&zyaaap93(:)+xyzzyaabh93*xyzzyaabr93*xyzzyaaao93(:))
if(sd)xyzzyaabc93=xyzzyaabc93+xyzzyaabj93*(xyzzyaabg93*(xyzzyaabq93+xy&
&zzyaabp93*xyzzyaaar93)+xyzzyaabh93*(xyzzyaabs93+xyzzyaabr93*xyzzyaaaq&
&93)+2.d0*xyzzyaabp93*xyzzyaabr93*xyzzyaaas93)
endif
if(xyzzyaaak93<xyzzyaacc1)then
if(xyzzyaabt93)then
if(sd)xyzzyaabs93=dble(xyzzyaaak93+1)*xyzzyaabr93
xyzzyaabr93=dble(xyzzyaaak93+1)*xyzzyaabg93
endif
xyzzyaabg93=xyzzyaabg93*xyzzyaaal93
endif
enddo
if(xyzzyaaaj93<xyzzyaacc1)then
if(xyzzyaabt93)then
if(sd)xyzzyaabq93=dble(xyzzyaaaj93+1)*xyzzyaabp93
xyzzyaabp93=dble(xyzzyaaaj93+1)*xyzzyaabh93
endif
xyzzyaabh93=xyzzyaabh93*xyzzyaaan93
endif
enddo
if(xyzzyaaai93<xyzzyaacc1)xyzzyaabi93=xyzzyaabi93*xyzzyaaam93
enddo
if(val)value_h=value_h+xyzzyaabd93*xyzzyaaba93
if(fd)grad_h=grad_h+xyzzyaabd93*xyzzyaabb93+xyzzyaabe93*xyzzyaaba93
if(sd)lap_h=lap_h+xyzzyaabd93*xyzzyaabc93+2.d0*ddot(dimensionality,xyz&
&zyaabe93(1),1,xyzzyaabb93(1),1)+xyzzyaabf93*xyzzyaaba93
enddo
enddo
end subroutine xyzzyaapy1
subroutine xyzzyaapz1(r,f,df,d2f,fsd)
implicit none
real(dp),intent(in) :: r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fsd
real(dp) xyzzyaaaa94,xyzzyaaab94,xyzzyaaac94,xyzzyaaad94
f=0.d0
df=0.d0
d2f=0.d0
xyzzyaaaa94=1.d0-r*xyzzyaacn1
if(fsd)then
xyzzyaaab94=xyzzyaaaa94*xyzzyaaaa94
select case(xyzzyaagx1)
case(0)
f=1.d0
case(1)
f=xyzzyaaaa94
df=-xyzzyaacn1
case(2)
f=xyzzyaaab94
df=-2.d0*xyzzyaacn1*xyzzyaaaa94
d2f=2.d0*xyzzyaaco1
case(3)
f=xyzzyaaab94*xyzzyaaaa94
df=-3.d0*xyzzyaacn1*xyzzyaaab94
d2f=6.d0*xyzzyaaco1*xyzzyaaaa94
case default
xyzzyaaac94=xyzzyaaaa94**(xyzzyaagx1-2)
xyzzyaaad94=xyzzyaaac94*xyzzyaaaa94
f=xyzzyaaad94*xyzzyaaaa94
df=-xyzzyaagx1*xyzzyaacn1*xyzzyaaad94
d2f=xyzzyaagx1*(xyzzyaagx1-1)*xyzzyaaco1*xyzzyaaac94
end select
else
select case(xyzzyaagx1)
case(0)
f=1.d0
case(1)
f=xyzzyaaaa94
case(2)
f=xyzzyaaaa94*xyzzyaaaa94
case default
f=xyzzyaaaa94**xyzzyaagx1
end select
endif
end subroutine xyzzyaapz1
subroutine xyzzyaaqa1(ri,ispin,set,val,fd,sd,value_chi,deriv_chi,sderi&
&v_chi)
implicit none
integer,intent(in) :: ispin,set
real(dp),intent(in) :: ri
real(dp),intent(out) :: value_chi,deriv_chi,sderiv_chi
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa95
real(dp) xyzzyaaab95,xyzzyaaac95,xyzzyaaad95,xyzzyaaae95,xyzzyaaaf95,x&
&yzzyaaag95,xyzzyaaah95,xyzzyaaai95,xyzzyaaaj95,xyzzyaaak95
if(ri>xyzzyaacu1(set).and.xyzzyaagx1>0)then
value_chi=0.d0
deriv_chi=0.d0
sderiv_chi=0.d0
else
xyzzyaaah95=ri-xyzzyaacu1(set)
if(xyzzyaagx1==2)then
if(sd)then
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaaf95=xyzzyaaal1(1,ispin,set)
xyzzyaaag95=0.d0
xyzzyaaac95=1.d0
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaad95=xyzzyaaac95*xyzzyaaaa95
xyzzyaaac95=xyzzyaaab95*xyzzyaaaa95
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
xyzzyaaaf95=xyzzyaaaf95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaac95
xyzzyaaag95=xyzzyaaag95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaad95
enddo
sderiv_chi=2*xyzzyaaae95+xyzzyaaah95*(4*xyzzyaaaf95+xyzzyaaah95*xyzzya&
&aag95)
if(fd)deriv_chi=xyzzyaaah95*(2*xyzzyaaae95+xyzzyaaah95*xyzzyaaaf95)
if(val)value_chi=xyzzyaaah95*xyzzyaaah95*xyzzyaaae95
elseif(fd)then
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaaf95=xyzzyaaal1(1,ispin,set)
xyzzyaaac95=1.d0
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaac95=xyzzyaaab95*xyzzyaaaa95
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
xyzzyaaaf95=xyzzyaaaf95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaac95
enddo
deriv_chi=xyzzyaaah95*(2*xyzzyaaae95+xyzzyaaah95*xyzzyaaaf95)
if(val)value_chi=xyzzyaaah95*xyzzyaaah95*xyzzyaaae95
else
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
enddo
value_chi=xyzzyaaah95*xyzzyaaah95*xyzzyaaae95
endif
elseif(xyzzyaagx1==3)then
xyzzyaaai95=xyzzyaaah95*xyzzyaaah95
if(sd)then
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaaf95=xyzzyaaal1(1,ispin,set)
xyzzyaaag95=0.d0
xyzzyaaac95=1.d0
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaad95=xyzzyaaac95*xyzzyaaaa95
xyzzyaaac95=xyzzyaaab95*xyzzyaaaa95
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
xyzzyaaaf95=xyzzyaaaf95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaac95
xyzzyaaag95=xyzzyaaag95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaad95
enddo
sderiv_chi=xyzzyaaah95*(6*xyzzyaaae95+xyzzyaaah95*(6*xyzzyaaaf95+xyzzy&
&aaah95*xyzzyaaag95))
if(fd)deriv_chi=xyzzyaaai95*(3*xyzzyaaae95+xyzzyaaah95*xyzzyaaaf95)
if(val)value_chi=xyzzyaaai95*xyzzyaaah95*xyzzyaaae95
elseif(fd)then
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaaf95=xyzzyaaal1(1,ispin,set)
xyzzyaaac95=1.d0
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaac95=xyzzyaaab95*xyzzyaaaa95
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
xyzzyaaaf95=xyzzyaaaf95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaac95
enddo
deriv_chi=xyzzyaaai95*(3*xyzzyaaae95+xyzzyaaah95*xyzzyaaaf95)
if(val)value_chi=xyzzyaaai95*xyzzyaaah95*xyzzyaaae95
else
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
enddo
value_chi=xyzzyaaai95*xyzzyaaah95*xyzzyaaae95
endif
else
if(sd)then
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaaf95=xyzzyaaal1(1,ispin,set)
xyzzyaaag95=0.d0
xyzzyaaac95=1.d0
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaad95=xyzzyaaac95*xyzzyaaaa95
xyzzyaaac95=xyzzyaaab95*xyzzyaaaa95
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
xyzzyaaaf95=xyzzyaaaf95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaac95
xyzzyaaag95=xyzzyaaag95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaad95
enddo
if(xyzzyaagx1>1)then
xyzzyaaak95=xyzzyaaah95**(xyzzyaagx1-2)
xyzzyaaaj95=xyzzyaaak95*xyzzyaaah95
sderiv_chi=xyzzyaaak95*(xyzzyaagx1*(xyzzyaagx1-1)*xyzzyaaae95+xyzzyaaa&
&h95*(2*xyzzyaagx1*xyzzyaaaf95+xyzzyaaah95*xyzzyaaag95))
if(fd)deriv_chi=xyzzyaaaj95*(xyzzyaagx1*xyzzyaaae95+xyzzyaaah95*xyzzya&
&aaf95)
if(val)value_chi=xyzzyaaaj95*xyzzyaaah95*xyzzyaaae95
elseif(xyzzyaagx1==1)then
sderiv_chi=2*xyzzyaaaf95+xyzzyaaah95*xyzzyaaag95
if(fd)deriv_chi=xyzzyaaae95+xyzzyaaah95*xyzzyaaaf95
if(val)value_chi=xyzzyaaah95*xyzzyaaae95
else
sderiv_chi=xyzzyaaag95
deriv_chi=xyzzyaaaf95
value_chi=xyzzyaaae95
endif
elseif(fd)then
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaaf95=xyzzyaaal1(1,ispin,set)
xyzzyaaac95=1.d0
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaac95=xyzzyaaab95*xyzzyaaaa95
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
xyzzyaaaf95=xyzzyaaaf95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaac95
enddo
if(xyzzyaagx1>0)then
xyzzyaaaj95=xyzzyaaah95**(xyzzyaagx1-1)
deriv_chi=xyzzyaaaj95*(xyzzyaagx1*xyzzyaaae95+xyzzyaaah95*xyzzyaaaf95)
if(val)value_chi=xyzzyaaaj95*xyzzyaaah95*xyzzyaaae95
else
deriv_chi=xyzzyaaaf95
value_chi=xyzzyaaae95
endif
else
xyzzyaaae95=xyzzyaaal1(0,ispin,set)+ri*xyzzyaaal1(1,ispin,set)
xyzzyaaab95=ri
do xyzzyaaaa95=2,xyzzyaace1(set)
xyzzyaaab95=xyzzyaaab95*ri
xyzzyaaae95=xyzzyaaae95+xyzzyaaal1(xyzzyaaaa95,ispin,set)*xyzzyaaab95
enddo
value_chi=xyzzyaaah95**xyzzyaagx1*xyzzyaaae95
endif
endif
endif
end subroutine xyzzyaaqa1
subroutine xyzzyaaqb1(j,ri,rj,rij,ispin,jspin,set,val,fd,sd,flag_j,val&
&ue_f,df_drii,df_drji,df_drij,d2f_drii2,d2f_drji2,d2f_drij2,d2f_driidr&
&ij,d2f_drjidrij)
implicit none
integer,intent(in) :: j,ispin,jspin,set
real(dp),intent(in) :: ri,rj,rij
real(dp),intent(out) :: value_f,df_drii,df_drij,d2f_drii2,d2f_drij2,d2&
&f_driidrij,df_drji,d2f_drji2,d2f_drjidrij
logical,intent(in) :: val,fd,sd,flag_j
integer xyzzyaaaa96,xyzzyaaab96,xyzzyaaac96,xyzzyaaad96,xyzzyaaae96,xy&
&zzyaaaf96,xyzzyaaag96,xyzzyaaah96
real(dp) xyzzyaaai96,xyzzyaaaj96,xyzzyaaak96,xyzzyaaal96,xyzzyaaam96,x&
&yzzyaaan96,xyzzyaaao96,xyzzyaaap96,xyzzyaaaq96,xyzzyaaar96,xyzzyaaas9&
&6,xyzzyaaat96,xyzzyaaau96,xyzzyaaav96,xyzzyaaaw96,xyzzyaaax96,xyzzyaa&
&ay96,xyzzyaaaz96,xyzzyaaba96
xyzzyaaal96=xyzzyaacv1(set)
if(ri>xyzzyaaal96.or.rj>xyzzyaaal96.and.xyzzyaagx1>0)then
if(val)value_f=0.d0
if(fd)then
df_drii=0.d0
df_drji=0.d0
df_drij=0.d0
endif
if(sd)then
d2f_drij2=0.d0
d2f_drii2=0.d0
d2f_drji2=0.d0
d2f_driidrij=0.d0
d2f_drjidrij=0.d0
endif
else
xyzzyaaap96=ri-xyzzyaaal96
xyzzyaaaq96=rj-xyzzyaaal96
xyzzyaaar96=xyzzyaaap96*xyzzyaaaq96
xyzzyaaaa96=which_spair(ispin,jspin,levels_spairs)
xyzzyaaaf96=xyzzyaacf1(set)
xyzzyaaag96=xyzzyaacg1(set)
do xyzzyaaad96=1,xyzzyaaaf96
xyzzyaadr1(xyzzyaaad96)=xyzzyaadr1(xyzzyaaad96-1)*rij
enddo
if(j==0)then
do xyzzyaaac96=1,xyzzyaaag96
xyzzyaado1(xyzzyaaac96)=xyzzyaado1(xyzzyaaac96-1)*ri
xyzzyaadp1(xyzzyaaac96)=xyzzyaadp1(xyzzyaaac96-1)*rj
enddo
else
do xyzzyaaac96=1,xyzzyaaag96
xyzzyaadp1(xyzzyaaac96)=xyzzyaadq1(xyzzyaaac96,j)
enddo
endif
if(sd)then
if(flag_j)then
xyzzyaadt1(2)=2*xyzzyaado1(1)
xyzzyaadv1(2)=2*xyzzyaadp1(1)
xyzzyaadu1(2)=2.d0
xyzzyaadw1(2)=2.d0
do xyzzyaaab96=3,xyzzyaaag96
xyzzyaadt1(xyzzyaaab96)=xyzzyaaab96*xyzzyaado1(xyzzyaaab96-1)
xyzzyaadv1(xyzzyaaab96)=xyzzyaaab96*xyzzyaadp1(xyzzyaaab96-1)
xyzzyaadu1(xyzzyaaab96)=xyzzyaaab96*xyzzyaadt1(xyzzyaaab96-1)
xyzzyaadw1(xyzzyaaab96)=xyzzyaaab96*xyzzyaadv1(xyzzyaaab96-1)
enddo
value_f=xyzzyaaam1(1,xyzzyaaaa96,set)
xyzzyaaae96=2
value_f=value_f+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzyaado1(1)
df_drii=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
d2f_drii2=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
value_f=value_f+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzyaado1(xyz&
&zyaaab96)
df_drii=df_drii+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzyaadt1(xyz&
&zyaaab96)
d2f_drii2=d2f_drii2+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzyaadu1&
&(xyzzyaaab96)
enddo
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
value_f=value_f+xyzzyaaat96*xyzzyaadp1(1)
df_drii=df_drii+xyzzyaaav96*xyzzyaadp1(1)
df_drji=xyzzyaaat96
d2f_drii2=d2f_drii2+xyzzyaaax96*xyzzyaadp1(1)
d2f_drji2=0.d0
do xyzzyaaac96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
value_f=value_f+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
df_drii=df_drii+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
df_drji=df_drji+xyzzyaaat96*xyzzyaadv1(xyzzyaaac96)
d2f_drii2=d2f_drii2+xyzzyaaax96*xyzzyaadp1(xyzzyaaac96)
d2f_drji2=d2f_drji2+xyzzyaaat96*xyzzyaadw1(xyzzyaaac96)
enddo
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaas96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaas96=xyzzyaaas96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaau96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaaw96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaas96=xyzzyaaas96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaaw96=xyzzyaaaw96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(1)
xyzzyaaau96=xyzzyaaau96+xyzzyaaav96*xyzzyaadp1(1)
xyzzyaaay96=xyzzyaaat96
xyzzyaaaw96=xyzzyaaaw96+xyzzyaaax96*xyzzyaadp1(1)
xyzzyaaaz96=0.d0
do xyzzyaaac96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaay96=xyzzyaaay96+xyzzyaaat96*xyzzyaadv1(xyzzyaaac96)
xyzzyaaaw96=xyzzyaaaw96+xyzzyaaax96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaaz96=xyzzyaaaz96+xyzzyaaat96*xyzzyaadw1(xyzzyaaac96)
enddo
value_f=value_f+xyzzyaaas96*xyzzyaadr1(1)
df_drii=df_drii+xyzzyaaau96*xyzzyaadr1(1)
df_drij=xyzzyaaas96
df_drji=df_drji+xyzzyaaay96*xyzzyaadr1(1)
d2f_drii2=d2f_drii2+xyzzyaaaw96*xyzzyaadr1(1)
d2f_driidrij=xyzzyaaau96
d2f_drji2=d2f_drji2+xyzzyaaaz96*xyzzyaadr1(1)
d2f_drjidrij=xyzzyaaay96
d2f_drij2=0.d0
do xyzzyaaad96=2,xyzzyaaaf96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaas96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaas96=xyzzyaaas96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaau96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaaw96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaas96=xyzzyaaas96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaaw96=xyzzyaaaw96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(1)
xyzzyaaau96=xyzzyaaau96+xyzzyaaav96*xyzzyaadp1(1)
xyzzyaaay96=xyzzyaaat96
xyzzyaaaw96=xyzzyaaaw96+xyzzyaaax96*xyzzyaadp1(1)
xyzzyaaaz96=0.d0
do xyzzyaaac96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaay96=xyzzyaaay96+xyzzyaaat96*xyzzyaadv1(xyzzyaaac96)
xyzzyaaaw96=xyzzyaaaw96+xyzzyaaax96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaaz96=xyzzyaaaz96+xyzzyaaat96*xyzzyaadw1(xyzzyaaac96)
enddo
value_f=value_f+xyzzyaaas96*xyzzyaadr1(xyzzyaaad96)
df_drii=df_drii+xyzzyaaau96*xyzzyaadr1(xyzzyaaad96)
df_drji=df_drji+xyzzyaaay96*xyzzyaadr1(xyzzyaaad96)
d2f_drii2=d2f_drii2+xyzzyaaaw96*xyzzyaadr1(xyzzyaaad96)
xyzzyaaba96=xyzzyaaad96*xyzzyaadr1(xyzzyaaad96-1)
df_drij=df_drij+xyzzyaaas96*xyzzyaaba96
d2f_driidrij=d2f_driidrij+xyzzyaaau96*xyzzyaaba96
d2f_drjidrij=d2f_drjidrij+xyzzyaaay96*xyzzyaaba96
d2f_drij2=d2f_drij2+xyzzyaaas96*xyzzyaaad96*(xyzzyaaad96-1)*xyzzyaadr1&
&(xyzzyaaad96-2)
d2f_drji2=d2f_drji2+xyzzyaaaz96*xyzzyaadr1(xyzzyaaad96)
enddo
else
xyzzyaadt1(2)=2*xyzzyaado1(1)
xyzzyaadu1(2)=2.d0
do xyzzyaaab96=3,xyzzyaaag96
xyzzyaadt1(xyzzyaaab96)=xyzzyaaab96*xyzzyaado1(xyzzyaaab96-1)
xyzzyaadu1(xyzzyaaab96)=xyzzyaaab96*xyzzyaadt1(xyzzyaaab96-1)
enddo
xyzzyaaae96=0
value_f=0.d0
df_drii=0.d0
d2f_drii2=0.d0
do xyzzyaaac96=0,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
value_f=value_f+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
df_drii=df_drii+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
d2f_drii2=d2f_drii2+xyzzyaaax96*xyzzyaadp1(xyzzyaaac96)
enddo
xyzzyaaas96=0.d0
xyzzyaaau96=0.d0
xyzzyaaaw96=0.d0
do xyzzyaaac96=0,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaaw96=xyzzyaaaw96+xyzzyaaax96*xyzzyaadp1(xyzzyaaac96)
enddo
value_f=value_f+xyzzyaaas96*xyzzyaadr1(1)
df_drii=df_drii+xyzzyaaau96*xyzzyaadr1(1)
df_drij=xyzzyaaas96
d2f_drii2=d2f_drii2+xyzzyaaaw96*xyzzyaadr1(1)
d2f_driidrij=xyzzyaaau96
d2f_drij2=0.d0
do xyzzyaaad96=2,xyzzyaaaf96
xyzzyaaas96=0.d0
xyzzyaaau96=0.d0
xyzzyaaaw96=0.d0
do xyzzyaaac96=0,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(1)
xyzzyaaav96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaax96=0.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
xyzzyaaax96=xyzzyaaax96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adu1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaaw96=xyzzyaaaw96+xyzzyaaax96*xyzzyaadp1(xyzzyaaac96)
enddo
value_f=value_f+xyzzyaaas96*xyzzyaadr1(xyzzyaaad96)
df_drii=df_drii+xyzzyaaau96*xyzzyaadr1(xyzzyaaad96)
d2f_drii2=d2f_drii2+xyzzyaaaw96*xyzzyaadr1(xyzzyaaad96)
xyzzyaaba96=xyzzyaaad96*xyzzyaadr1(xyzzyaaad96-1)
df_drij=df_drij+xyzzyaaas96*xyzzyaaba96
d2f_driidrij=d2f_driidrij+xyzzyaaau96*xyzzyaaba96
d2f_drij2=d2f_drij2+xyzzyaaas96*xyzzyaaad96*(xyzzyaaad96-1)*xyzzyaadr1&
&(xyzzyaaad96-2)
enddo
endif
if(xyzzyaagx1==1)then
d2f_drij2=d2f_drij2*xyzzyaaar96
d2f_drii2=2*xyzzyaaaq96*df_drii+xyzzyaaar96*d2f_drii2
d2f_driidrij=xyzzyaaaq96*df_drij+xyzzyaaar96*d2f_driidrij
if(flag_j)then
d2f_drji2=2*xyzzyaaap96*df_drji+xyzzyaaar96*d2f_drji2
d2f_drjidrij=xyzzyaaap96*df_drij+xyzzyaaar96*d2f_drjidrij
endif
if(fd)then
df_drij=df_drij*xyzzyaaar96
df_drii=xyzzyaaaq96*value_f+xyzzyaaar96*df_drii
if(flag_j)df_drji=xyzzyaaar96*value_f+xyzzyaaar96*df_drji
endif
if(val)value_f=value_f*xyzzyaaar96
elseif(xyzzyaagx1>1)then
select case(xyzzyaagx1)
case(2)
xyzzyaaah96=2
xyzzyaaao96=1.d0
case(3)
xyzzyaaah96=6
xyzzyaaao96=xyzzyaaar96
case default
xyzzyaaah96=xyzzyaagx1*(xyzzyaagx1-1)
xyzzyaaao96=xyzzyaaar96**(xyzzyaagx1-2)
end select
xyzzyaaan96=xyzzyaaao96*xyzzyaaar96
xyzzyaaam96=xyzzyaaan96*xyzzyaaar96
xyzzyaaai96=xyzzyaaah96*value_f
xyzzyaaaj96=2*xyzzyaagx1*xyzzyaaar96
xyzzyaaak96=xyzzyaagx1*df_drij
d2f_drij2=d2f_drij2*xyzzyaaam96
d2f_drii2=xyzzyaaao96*((xyzzyaaai96*xyzzyaaaq96+xyzzyaaaj96*df_drii)*x&
&yzzyaaaq96+xyzzyaaar96**2*d2f_drii2)
d2f_driidrij=xyzzyaaan96*(xyzzyaaak96*xyzzyaaaq96+xyzzyaaar96*d2f_drii&
&drij)
if(flag_j)then
d2f_drji2=xyzzyaaao96*((xyzzyaaai96*xyzzyaaap96+xyzzyaaaj96*df_drji)*x&
&yzzyaaap96+xyzzyaaar96**2*d2f_drji2)
d2f_drjidrij=xyzzyaaan96*(xyzzyaaak96*xyzzyaaap96+xyzzyaaar96*d2f_drji&
&drij)
endif
if(fd)then
xyzzyaaai96=xyzzyaagx1*value_f
df_drij=df_drij*xyzzyaaam96
df_drii=xyzzyaaan96*(xyzzyaaai96*xyzzyaaaq96+xyzzyaaar96*df_drii)
if(flag_j)df_drji=xyzzyaaan96*(xyzzyaaai96*xyzzyaaap96+xyzzyaaar96*df_&
&drji)
endif
if(val)value_f=value_f*xyzzyaaam96
endif
elseif(fd)then
if(flag_j)then
xyzzyaadt1(1)=1.d0
xyzzyaadv1(1)=1.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaadt1(xyzzyaaab96)=xyzzyaaab96*xyzzyaado1(xyzzyaaab96-1)
xyzzyaadv1(xyzzyaaab96)=xyzzyaaab96*xyzzyaadp1(xyzzyaaab96-1)
enddo
xyzzyaaae96=1
value_f=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
df_drii=0.d0
do xyzzyaaab96=1,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
value_f=value_f+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzyaado1(xyz&
&zyaaab96)
df_drii=df_drii+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzyaadt1(xyz&
&zyaaab96)
enddo
df_drji=0.d0
do xyzzyaaac96=1,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaav96=0.d0
do xyzzyaaab96=1,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
enddo
value_f=value_f+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
df_drii=df_drii+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
df_drji=df_drji+xyzzyaaat96*xyzzyaadv1(xyzzyaaac96)
enddo
df_drij=0.d0
do xyzzyaaad96=1,xyzzyaaaf96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaas96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaau96=0.d0
do xyzzyaaab96=1,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaas96=xyzzyaaas96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
enddo
xyzzyaaay96=0.d0
do xyzzyaaac96=1,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaav96=0.d0
do xyzzyaaab96=1,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaay96=xyzzyaaay96+xyzzyaaat96*xyzzyaadv1(xyzzyaaac96)
enddo
value_f=value_f+xyzzyaaas96*xyzzyaadr1(xyzzyaaad96)
df_drii=df_drii+xyzzyaaau96*xyzzyaadr1(xyzzyaaad96)
df_drij=df_drij+xyzzyaaad96*xyzzyaaas96*xyzzyaadr1(xyzzyaaad96-1)
df_drji=df_drji+xyzzyaaay96*xyzzyaadr1(xyzzyaaad96)
enddo
else
xyzzyaadt1(1)=1.d0
do xyzzyaaab96=2,xyzzyaaag96
xyzzyaadt1(xyzzyaaab96)=xyzzyaaab96*xyzzyaado1(xyzzyaaab96-1)
enddo
xyzzyaaae96=0
value_f=0.d0
df_drii=0.d0
do xyzzyaaac96=0,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaav96=0.d0
do xyzzyaaab96=1,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
enddo
value_f=value_f+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
df_drii=df_drii+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
enddo
df_drij=0.d0
do xyzzyaaad96=1,xyzzyaaaf96
xyzzyaaas96=0.d0
xyzzyaaau96=0.d0
do xyzzyaaac96=0,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)
xyzzyaaav96=0.d0
do xyzzyaaab96=1,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
xyzzyaaav96=xyzzyaaav96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&adt1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
xyzzyaaau96=xyzzyaaau96+xyzzyaaav96*xyzzyaadp1(xyzzyaaac96)
enddo
value_f=value_f+xyzzyaaas96*xyzzyaadr1(xyzzyaaad96)
df_drii=df_drii+xyzzyaaau96*xyzzyaadr1(xyzzyaaad96)
df_drij=df_drij+xyzzyaaad96*xyzzyaaas96*xyzzyaadr1(xyzzyaaad96-1)
enddo
endif
if(xyzzyaagx1>0)then
select case(xyzzyaagx1)
case(1)
xyzzyaaan96=1.d0
case(2)
xyzzyaaan96=xyzzyaaar96
case(3)
xyzzyaaan96=xyzzyaaar96*xyzzyaaar96
case default
xyzzyaaan96=xyzzyaaar96**(xyzzyaagx1-1)
end select
xyzzyaaam96=xyzzyaaan96*xyzzyaaar96
df_drij=df_drij*xyzzyaaam96
df_drii=xyzzyaaan96*xyzzyaaaq96*xyzzyaagx1*value_f+xyzzyaaam96*df_drii
if(flag_j)df_drji=xyzzyaaan96*xyzzyaaap96*(xyzzyaagx1*value_f+xyzzyaaa&
&q96*df_drji)
if(val)value_f=value_f*xyzzyaaam96
endif
else
value_f=0.d0
xyzzyaaae96=0
do xyzzyaaad96=0,xyzzyaaaf96
xyzzyaaas96=0.d0
do xyzzyaaac96=0,xyzzyaaag96
xyzzyaaat96=0.d0
do xyzzyaaab96=0,xyzzyaaag96
xyzzyaaae96=xyzzyaaae96+1
xyzzyaaat96=xyzzyaaat96+xyzzyaaam1(xyzzyaaae96,xyzzyaaaa96,set)*xyzzya&
&ado1(xyzzyaaab96)
enddo
xyzzyaaas96=xyzzyaaas96+xyzzyaaat96*xyzzyaadp1(xyzzyaaac96)
enddo
value_f=value_f+xyzzyaaas96*xyzzyaadr1(xyzzyaaad96)
enddo
value_f=value_f*xyzzyaaar96**xyzzyaagx1
endif
endif
end subroutine xyzzyaaqb1
subroutine xyzzyaaqc1(j,ri,rj,rij,ispin,jspin,set,val,fd,sd,flag_j,val&
&ue_f,df_drii,df_drji,df_drij,d2f_drii2,d2f_drji2,d2f_drij2,d2f_driidr&
&ij,d2f_drjidrij)
implicit none
integer,intent(in) :: j,ispin,jspin,set
real(dp),intent(in) :: ri,rj,rij
real(dp),intent(out) :: value_f,df_drii,df_drij,d2f_drii2,d2f_drij2,d2&
&f_driidrij,df_drji,d2f_drji2,d2f_drjidrij
logical,intent(in) :: val,fd,sd,flag_j
integer xyzzyaaaa97,xyzzyaaab97,xyzzyaaac97,xyzzyaaad97,xyzzyaaae97,xy&
&zzyaaaf97
real(dp) xyzzyaaag97,xyzzyaaah97,xyzzyaaai97,xyzzyaaaj97,xyzzyaaak97,x&
&yzzyaaal97,xyzzyaaam97,xyzzyaaan97,xyzzyaaao97,xyzzyaaap97,xyzzyaaaq9&
&7,xyzzyaaar97,xyzzyaaas97,xyzzyaaat97,xyzzyaaau97,xyzzyaaav97,xyzzyaa&
&aw97,xyzzyaaax97,xyzzyaaay97
xyzzyaaaj97=xyzzyaacv1(set)
if(ri>xyzzyaaaj97.or.rj>xyzzyaaaj97.and.xyzzyaagx1>0)then
if(val)value_f=0.d0
if(fd)then
df_drii=0.d0
df_drji=0.d0
df_drij=0.d0
endif
if(sd)then
d2f_drij2=0.d0
d2f_drii2=0.d0
d2f_drji2=0.d0
d2f_driidrij=0.d0
d2f_drjidrij=0.d0
endif
else
xyzzyaaan97=ri-xyzzyaaaj97
xyzzyaaao97=rj-xyzzyaaaj97
xyzzyaaap97=xyzzyaaan97*xyzzyaaao97
xyzzyaaaa97=which_spair(ispin,jspin,levels_spairs)
xyzzyaadr1(1)=xyzzyaadr1(0)*rij
xyzzyaadr1(2)=xyzzyaadr1(1)*rij
if(j==0)then
xyzzyaado1(1)=xyzzyaado1(0)*ri
xyzzyaado1(2)=xyzzyaado1(1)*ri
xyzzyaadp1(1)=xyzzyaadp1(0)*rj
xyzzyaadp1(2)=xyzzyaadp1(1)*rj
else
xyzzyaadp1(1)=xyzzyaadq1(1,j)
xyzzyaadp1(2)=xyzzyaadq1(2,j)
endif
if(sd)then
if(flag_j)then
xyzzyaadt1(2)=2*xyzzyaado1(1)
xyzzyaadv1(2)=2*xyzzyaadp1(1)
xyzzyaadu1(2)=2.d0
xyzzyaadw1(2)=2.d0
value_f=xyzzyaaam1(1,xyzzyaaaa97,set)
value_f=value_f+xyzzyaaam1(2,xyzzyaaaa97,set)*xyzzyaado1(1)
df_drii=xyzzyaaam1(2,xyzzyaaaa97,set)
d2f_drii2=0.d0
value_f=value_f+xyzzyaaam1(3,xyzzyaaaa97,set)*xyzzyaado1(2)
df_drii=df_drii+xyzzyaaam1(3,xyzzyaaaa97,set)*xyzzyaadt1(2)
d2f_drii2=d2f_drii2+xyzzyaaam1(3,xyzzyaaaa97,set)*xyzzyaadu1(2)
xyzzyaaar97=xyzzyaaam1(4,xyzzyaaaa97,set)
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(5,xyzzyaaaa97,set)*xyzzyaado1(1)
xyzzyaaat97=xyzzyaaam1(5,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(6,xyzzyaaaa97,set)*xyzzyaado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(6,xyzzyaaaa97,set)*xyzzyaadt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(6,xyzzyaaaa97,set)*xyzzyaadu1(2)
value_f=value_f+xyzzyaaar97*xyzzyaadp1(1)
df_drii=df_drii+xyzzyaaat97*xyzzyaadp1(1)
df_drji=xyzzyaaar97
d2f_drii2=d2f_drii2+xyzzyaaav97*xyzzyaadp1(1)
d2f_drji2=0.d0
xyzzyaaar97=xyzzyaaam1(7,xyzzyaaaa97,set)
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(8,xyzzyaaaa97,set)*xyzzyaado1(1)
xyzzyaaat97=xyzzyaaam1(8,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(9,xyzzyaaaa97,set)*xyzzyaado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(9,xyzzyaaaa97,set)*xyzzyaadt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(9,xyzzyaaaa97,set)*xyzzyaadu1(2)
value_f=value_f+xyzzyaaar97*xyzzyaadp1(2)
df_drii=df_drii+xyzzyaaat97*xyzzyaadp1(2)
df_drji=df_drji+xyzzyaaar97*xyzzyaadv1(2)
d2f_drii2=d2f_drii2+xyzzyaaav97*xyzzyaadp1(2)
d2f_drji2=d2f_drji2+xyzzyaaar97*xyzzyaadw1(2)
xyzzyaaaq97=xyzzyaaam1(10,xyzzyaaaa97,set)
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaam1(11,xyzzyaaaa97,set)*xyzzyaado1(1)
xyzzyaaas97=xyzzyaaam1(11,xyzzyaaaa97,set)
xyzzyaaau97=0.d0
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaam1(12,xyzzyaaaa97,set)*xyzzyaado1(2)
xyzzyaaas97=xyzzyaaas97+xyzzyaaam1(12,xyzzyaaaa97,set)*xyzzyaadt1(2)
xyzzyaaau97=xyzzyaaau97+xyzzyaaam1(12,xyzzyaaaa97,set)*xyzzyaadu1(2)
xyzzyaaar97=xyzzyaaam1(13,xyzzyaaaa97,set)
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(14,xyzzyaaaa97,set)*xyzzyaado1(1)
xyzzyaaat97=xyzzyaaam1(14,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(15,xyzzyaaaa97,set)*xyzzyaado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(15,xyzzyaaaa97,set)*xyzzyaadt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(15,xyzzyaaaa97,set)*xyzzyaadu1(2)
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(1)
xyzzyaaas97=xyzzyaaas97+xyzzyaaat97*xyzzyaadp1(1)
xyzzyaaaw97=xyzzyaaar97
xyzzyaaau97=xyzzyaaau97+xyzzyaaav97*xyzzyaadp1(1)
xyzzyaaax97=0.d0
xyzzyaaar97=xyzzyaaam1(16,xyzzyaaaa97,set)
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(17,xyzzyaaaa97,set)*xyzzyaado1(1)
xyzzyaaat97=xyzzyaaam1(17,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(18,xyzzyaaaa97,set)*xyzzyaado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(18,xyzzyaaaa97,set)*xyzzyaadt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(18,xyzzyaaaa97,set)*xyzzyaadu1(2)
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(2)
xyzzyaaas97=xyzzyaaas97+xyzzyaaat97*xyzzyaadp1(2)
xyzzyaaaw97=xyzzyaaaw97+xyzzyaaar97*xyzzyaadv1(2)
xyzzyaaau97=xyzzyaaau97+xyzzyaaav97*xyzzyaadp1(2)
xyzzyaaax97=xyzzyaaax97+xyzzyaaar97*xyzzyaadw1(2)
value_f=value_f+xyzzyaaaq97*xyzzyaadr1(1)
df_drii=df_drii+xyzzyaaas97*xyzzyaadr1(1)
df_drij=xyzzyaaaq97
df_drji=df_drji+xyzzyaaaw97*xyzzyaadr1(1)
d2f_drii2=d2f_drii2+xyzzyaaau97*xyzzyaadr1(1)
d2f_driidrij=xyzzyaaas97
d2f_drji2=d2f_drji2+xyzzyaaax97*xyzzyaadr1(1)
d2f_drjidrij=xyzzyaaaw97
d2f_drij2=0.d0
xyzzyaaaq97=xyzzyaaam1(19,xyzzyaaaa97,set)
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaam1(20,xyzzyaaaa97,set)*xyzzyaado1(1)
xyzzyaaas97=xyzzyaaam1(20,xyzzyaaaa97,set)
xyzzyaaau97=0.d0
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaam1(21,xyzzyaaaa97,set)*xyzzyaado1(2)
xyzzyaaas97=xyzzyaaas97+xyzzyaaam1(21,xyzzyaaaa97,set)*xyzzyaadt1(2)
xyzzyaaau97=xyzzyaaau97+xyzzyaaam1(21,xyzzyaaaa97,set)*xyzzyaadu1(2)
xyzzyaaar97=xyzzyaaam1(22,xyzzyaaaa97,set)
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(23,xyzzyaaaa97,set)*xyzzyaado1(1)
xyzzyaaat97=xyzzyaaam1(23,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(24,xyzzyaaaa97,set)*xyzzyaado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(24,xyzzyaaaa97,set)*xyzzyaadt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(24,xyzzyaaaa97,set)*xyzzyaadu1(2)
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(1)
xyzzyaaas97=xyzzyaaas97+xyzzyaaat97*xyzzyaadp1(1)
xyzzyaaaw97=xyzzyaaar97
xyzzyaaau97=xyzzyaaau97+xyzzyaaav97*xyzzyaadp1(1)
xyzzyaaax97=0.d0
xyzzyaaar97=xyzzyaaam1(25,xyzzyaaaa97,set)
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(26,xyzzyaaaa97,set)*xyzzyaado1(1)
xyzzyaaat97=xyzzyaaam1(26,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(27,xyzzyaaaa97,set)*xyzzyaado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(27,xyzzyaaaa97,set)*xyzzyaadt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(27,xyzzyaaaa97,set)*xyzzyaadu1(2)
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(2)
xyzzyaaas97=xyzzyaaas97+xyzzyaaat97*xyzzyaadp1(2)
xyzzyaaaw97=xyzzyaaaw97+xyzzyaaar97*xyzzyaadv1(2)
xyzzyaaau97=xyzzyaaau97+xyzzyaaav97*xyzzyaadp1(2)
xyzzyaaax97=xyzzyaaax97+xyzzyaaar97*xyzzyaadw1(2)
value_f=value_f+xyzzyaaaq97*xyzzyaadr1(2)
df_drii=df_drii+xyzzyaaas97*xyzzyaadr1(2)
df_drji=df_drji+xyzzyaaaw97*xyzzyaadr1(2)
d2f_drii2=d2f_drii2+xyzzyaaau97*xyzzyaadr1(2)
xyzzyaaay97=2.d0*xyzzyaadr1(1)
df_drij=df_drij+xyzzyaaaq97*xyzzyaaay97
d2f_driidrij=d2f_driidrij+xyzzyaaas97*xyzzyaaay97
d2f_drjidrij=d2f_drjidrij+xyzzyaaaw97*xyzzyaaay97
d2f_drij2=d2f_drij2+xyzzyaaaq97*2.d0*xyzzyaadr1(0)
d2f_drji2=d2f_drji2+xyzzyaaax97*xyzzyaadr1(2)
else
xyzzyaadt1(2)=2*xyzzyaado1(1)
xyzzyaadu1(2)=2.d0
xyzzyaaae97=0
value_f=0.d0
df_drii=0.d0
d2f_drii2=0.d0
do xyzzyaaac97=0,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(1)
xyzzyaaat97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adu1(2)
value_f=value_f+xyzzyaaar97*xyzzyaadp1(xyzzyaaac97)
df_drii=df_drii+xyzzyaaat97*xyzzyaadp1(xyzzyaaac97)
d2f_drii2=d2f_drii2+xyzzyaaav97*xyzzyaadp1(xyzzyaaac97)
enddo
xyzzyaaaq97=0.d0
xyzzyaaas97=0.d0
xyzzyaaau97=0.d0
do xyzzyaaac97=0,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(1)
xyzzyaaat97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adu1(2)
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(xyzzyaaac97)
xyzzyaaas97=xyzzyaaas97+xyzzyaaat97*xyzzyaadp1(xyzzyaaac97)
xyzzyaaau97=xyzzyaaau97+xyzzyaaav97*xyzzyaadp1(xyzzyaaac97)
enddo
value_f=value_f+xyzzyaaaq97*xyzzyaadr1(1)
df_drii=df_drii+xyzzyaaas97*xyzzyaadr1(1)
df_drij=xyzzyaaaq97
d2f_drii2=d2f_drii2+xyzzyaaau97*xyzzyaadr1(1)
d2f_driidrij=xyzzyaaas97
d2f_drij2=0.d0
xyzzyaaaq97=0.d0
xyzzyaaas97=0.d0
xyzzyaaau97=0.d0
do xyzzyaaac97=0,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(1)
xyzzyaaat97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaav97=0.d0
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(2)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adt1(2)
xyzzyaaav97=xyzzyaaav97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adu1(2)
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(xyzzyaaac97)
xyzzyaaas97=xyzzyaaas97+xyzzyaaat97*xyzzyaadp1(xyzzyaaac97)
xyzzyaaau97=xyzzyaaau97+xyzzyaaav97*xyzzyaadp1(xyzzyaaac97)
enddo
value_f=value_f+xyzzyaaaq97*xyzzyaadr1(2)
df_drii=df_drii+xyzzyaaas97*xyzzyaadr1(2)
d2f_drii2=d2f_drii2+xyzzyaaau97*xyzzyaadr1(2)
xyzzyaaay97=2.d0*xyzzyaadr1(1)
df_drij=df_drij+xyzzyaaaq97*xyzzyaaay97
d2f_driidrij=d2f_driidrij+xyzzyaaas97*xyzzyaaay97
d2f_drij2=d2f_drij2+xyzzyaaaq97*2.d0*xyzzyaadr1(0)
endif
if(xyzzyaagx1==1)then
d2f_drij2=d2f_drij2*xyzzyaaap97
d2f_drii2=2*xyzzyaaao97*df_drii+xyzzyaaap97*d2f_drii2
d2f_driidrij=xyzzyaaao97*df_drij+xyzzyaaap97*d2f_driidrij
if(flag_j)then
d2f_drji2=2*xyzzyaaan97*df_drji+xyzzyaaap97*d2f_drji2
d2f_drjidrij=xyzzyaaan97*df_drij+xyzzyaaap97*d2f_drjidrij
endif
if(fd)then
df_drij=df_drij*xyzzyaaap97
df_drii=xyzzyaaao97*value_f+xyzzyaaap97*df_drii
if(flag_j)df_drji=xyzzyaaap97*value_f+xyzzyaaap97*df_drji
endif
if(val)value_f=value_f*xyzzyaaap97
elseif(xyzzyaagx1>1)then
select case(xyzzyaagx1)
case(2)
xyzzyaaaf97=2
xyzzyaaam97=1.d0
case(3)
xyzzyaaaf97=6
xyzzyaaam97=xyzzyaaap97
case default
xyzzyaaaf97=xyzzyaagx1*(xyzzyaagx1-1)
xyzzyaaam97=xyzzyaaap97**(xyzzyaagx1-2)
end select
xyzzyaaal97=xyzzyaaam97*xyzzyaaap97
xyzzyaaak97=xyzzyaaal97*xyzzyaaap97
xyzzyaaag97=xyzzyaaaf97*value_f
xyzzyaaah97=2.d0*xyzzyaagx1*xyzzyaaap97
xyzzyaaai97=xyzzyaagx1*df_drij
d2f_drij2=d2f_drij2*xyzzyaaak97
d2f_drii2=xyzzyaaam97*((xyzzyaaag97*xyzzyaaao97+xyzzyaaah97*df_drii)*x&
&yzzyaaao97+xyzzyaaap97**2*d2f_drii2)
d2f_driidrij=xyzzyaaal97*(xyzzyaaai97*xyzzyaaao97+xyzzyaaap97*d2f_drii&
&drij)
if(flag_j)then
d2f_drji2=xyzzyaaam97*((xyzzyaaag97*xyzzyaaan97+xyzzyaaah97*df_drji)*x&
&yzzyaaan97+xyzzyaaap97**2*d2f_drji2)
d2f_drjidrij=xyzzyaaal97*(xyzzyaaai97*xyzzyaaan97+xyzzyaaap97*d2f_drji&
&drij)
endif
if(fd)then
xyzzyaaag97=xyzzyaagx1*value_f
df_drij=df_drij*xyzzyaaak97
df_drii=xyzzyaaal97*(xyzzyaaag97*xyzzyaaao97+xyzzyaaap97*df_drii)
if(flag_j)df_drji=xyzzyaaal97*(xyzzyaaag97*xyzzyaaan97+xyzzyaaap97*df_&
&drji)
endif
if(val)value_f=value_f*xyzzyaaak97
endif
elseif(fd)then
if(flag_j)then
xyzzyaadt1(1)=1.d0
xyzzyaadv1(1)=1.d0
xyzzyaadt1(2)=2.d0*xyzzyaado1(1)
xyzzyaadv1(2)=2.d0*xyzzyaadp1(1)
xyzzyaaae97=1
value_f=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
df_drii=0.d0
do xyzzyaaab97=1,2
xyzzyaaae97=xyzzyaaae97+1
value_f=value_f+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzyaado1(xyz&
&zyaaab97)
df_drii=df_drii+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzyaadt1(xyz&
&zyaaab97)
enddo
df_drji=0.d0
do xyzzyaaac97=1,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaat97=0.d0
do xyzzyaaab97=1,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(xyzzyaaab97)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adt1(xyzzyaaab97)
enddo
value_f=value_f+xyzzyaaar97*xyzzyaadp1(xyzzyaaac97)
df_drii=df_drii+xyzzyaaat97*xyzzyaadp1(xyzzyaaac97)
df_drji=df_drji+xyzzyaaar97*xyzzyaadv1(xyzzyaaac97)
enddo
df_drij=0.d0
do xyzzyaaad97=1,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaaq97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaas97=0.d0
do xyzzyaaab97=1,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(xyzzyaaab97)
xyzzyaaas97=xyzzyaaas97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adt1(xyzzyaaab97)
enddo
xyzzyaaaw97=0.d0
do xyzzyaaac97=1,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaat97=0.d0
do xyzzyaaab97=1,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(xyzzyaaab97)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adt1(xyzzyaaab97)
enddo
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(xyzzyaaac97)
xyzzyaaas97=xyzzyaaas97+xyzzyaaat97*xyzzyaadp1(xyzzyaaac97)
xyzzyaaaw97=xyzzyaaaw97+xyzzyaaar97*xyzzyaadv1(xyzzyaaac97)
enddo
value_f=value_f+xyzzyaaaq97*xyzzyaadr1(xyzzyaaad97)
df_drii=df_drii+xyzzyaaas97*xyzzyaadr1(xyzzyaaad97)
df_drij=df_drij+xyzzyaaad97*xyzzyaaaq97*xyzzyaadr1(xyzzyaaad97-1)
df_drji=df_drji+xyzzyaaaw97*xyzzyaadr1(xyzzyaaad97)
enddo
else
xyzzyaadt1(1)=1.d0
xyzzyaadt1(2)=2.d0*xyzzyaado1(1)
xyzzyaaae97=0
value_f=0.d0
df_drii=0.d0
do xyzzyaaac97=0,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaat97=0.d0
do xyzzyaaab97=1,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(xyzzyaaab97)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adt1(xyzzyaaab97)
enddo
value_f=value_f+xyzzyaaar97*xyzzyaadp1(xyzzyaaac97)
df_drii=df_drii+xyzzyaaat97*xyzzyaadp1(xyzzyaaac97)
enddo
df_drij=0.d0
do xyzzyaaad97=1,2
xyzzyaaaq97=0.d0
xyzzyaaas97=0.d0
do xyzzyaaac97=0,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)
xyzzyaaat97=0.d0
do xyzzyaaab97=1,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(xyzzyaaab97)
xyzzyaaat97=xyzzyaaat97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&adt1(xyzzyaaab97)
enddo
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(xyzzyaaac97)
xyzzyaaas97=xyzzyaaas97+xyzzyaaat97*xyzzyaadp1(xyzzyaaac97)
enddo
value_f=value_f+xyzzyaaaq97*xyzzyaadr1(xyzzyaaad97)
df_drii=df_drii+xyzzyaaas97*xyzzyaadr1(xyzzyaaad97)
df_drij=df_drij+xyzzyaaad97*xyzzyaaaq97*xyzzyaadr1(xyzzyaaad97-1)
enddo
endif
if(xyzzyaagx1>0)then
select case(xyzzyaagx1)
case(1)
xyzzyaaal97=1.d0
case(2)
xyzzyaaal97=xyzzyaaap97
case(3)
xyzzyaaal97=xyzzyaaap97*xyzzyaaap97
case default
xyzzyaaal97=xyzzyaaap97**(xyzzyaagx1-1)
end select
xyzzyaaak97=xyzzyaaal97*xyzzyaaap97
df_drij=df_drij*xyzzyaaak97
df_drii=xyzzyaaal97*xyzzyaaao97*xyzzyaagx1*value_f+xyzzyaaak97*df_drii
if(flag_j)df_drji=xyzzyaaal97*xyzzyaaan97*(xyzzyaagx1*value_f+xyzzyaaa&
&o97*df_drji)
if(val)value_f=value_f*xyzzyaaak97
endif
else
value_f=0.d0
xyzzyaaae97=0
do xyzzyaaad97=0,2
xyzzyaaaq97=0.d0
do xyzzyaaac97=0,2
xyzzyaaar97=0.d0
do xyzzyaaab97=0,2
xyzzyaaae97=xyzzyaaae97+1
xyzzyaaar97=xyzzyaaar97+xyzzyaaam1(xyzzyaaae97,xyzzyaaaa97,set)*xyzzya&
&ado1(xyzzyaaab97)
enddo
xyzzyaaaq97=xyzzyaaaq97+xyzzyaaar97*xyzzyaadp1(xyzzyaaac97)
enddo
value_f=value_f+xyzzyaaaq97*xyzzyaadr1(xyzzyaaad97)
enddo
value_f=value_f*xyzzyaaap97**xyzzyaagx1
endif
endif
end subroutine xyzzyaaqc1
subroutine xyzzyaaqd1(r_ij,ispin,jspin,val,fd,sd,value_p,grad_p,lap_p)
implicit none
integer,intent(in) :: ispin,jspin
real(dp),intent(in) :: r_ij(3)
real(dp),intent(out) :: value_p,grad_p(3),lap_p
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa98,xyzzyaaab98,xyzzyaaac98,xyzzyaaad98
real(dp) xyzzyaaae98,xyzzyaaaf98(3),xyzzyaaag98,xyzzyaaah98,xyzzyaaai9&
&8,xyzzyaaaj98
complex(dp) xyzzyaaak98,xyzzyaaal98,xyzzyaaam98,xyzzyaaan98,xyzzyaaao9&
&8
xyzzyaaaa98=which_spair(ispin,jspin,levels_spairs)
xyzzyaaah98=b1(1)*r_ij(1)+b1(2)*r_ij(2)+b1(3)*r_ij(3)
xyzzyaaal98=cmplx(cos(xyzzyaaah98),sin(xyzzyaaah98),dp)
xyzzyaaah98=b2(1)*r_ij(1)+b2(2)*r_ij(2)+b2(3)*r_ij(3)
xyzzyaaam98=cmplx(cos(xyzzyaaah98),sin(xyzzyaaah98),dp)
xyzzyaaah98=b3(1)*r_ij(1)+b3(2)*r_ij(2)+b3(3)*r_ij(3)
xyzzyaaan98=cmplx(cos(xyzzyaaah98),sin(xyzzyaaah98),dp)
xyzzyaagk1(0)=c_one
do xyzzyaaab98=1,xyzzyaage1(1)
xyzzyaagk1(xyzzyaaab98)=xyzzyaagk1(xyzzyaaab98-1)*xyzzyaaal98
enddo
xyzzyaagl1(0)=c_one
do xyzzyaaab98=1,xyzzyaage1(2)
xyzzyaagl1(xyzzyaaab98)=xyzzyaagl1(xyzzyaaab98-1)*xyzzyaaam98
xyzzyaagl1(-xyzzyaaab98)=conjg(xyzzyaagl1(xyzzyaaab98))
enddo
xyzzyaagm1(0)=c_one
do xyzzyaaab98=1,xyzzyaage1(3)
xyzzyaagm1(xyzzyaaab98)=xyzzyaagm1(xyzzyaaab98-1)*xyzzyaaan98
xyzzyaagm1(-xyzzyaaab98)=conjg(xyzzyaagm1(xyzzyaaab98))
enddo
if(.not.fd.and..not.sd)then
value_p=0.d0
xyzzyaaad98=0
do xyzzyaaab98=1,xyzzyaagq1
xyzzyaaae98=0.d0
do xyzzyaaac98=1,xyzzyaagt1(xyzzyaaab98)
xyzzyaaad98=xyzzyaaad98+1
xyzzyaaao98=xyzzyaagl1(xyzzyaaga1(2,xyzzyaaad98))*xyzzyaagm1(xyzzyaaga&
&1(3,xyzzyaaad98))
xyzzyaaae98=xyzzyaaae98+dble(xyzzyaagk1(xyzzyaaga1(1,xyzzyaaad98)))*db&
&le(xyzzyaaao98)-aimag(xyzzyaagk1(xyzzyaaga1(1,xyzzyaaad98)))*aimag(xy&
&zzyaaao98)
enddo
value_p=value_p+xyzzyaaae98*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
enddo
elseif(.not.val.and..not.sd)then
grad_p(1)=0.d0
grad_p(2)=0.d0
grad_p(3)=0.d0
xyzzyaaad98=0
do xyzzyaaab98=1,xyzzyaagq1
xyzzyaaaf98(1)=0.d0
xyzzyaaaf98(2)=0.d0
xyzzyaaaf98(3)=0.d0
do xyzzyaaac98=1,xyzzyaagt1(xyzzyaaab98)
xyzzyaaad98=xyzzyaaad98+1
xyzzyaaao98=xyzzyaagl1(xyzzyaaga1(2,xyzzyaaad98))*xyzzyaagm1(xyzzyaaga&
&1(3,xyzzyaaad98))
xyzzyaaaj98=dble(xyzzyaagk1(xyzzyaaga1(1,xyzzyaaad98)))*aimag(xyzzyaaa&
&o98)+aimag(xyzzyaagk1(xyzzyaaga1(1,xyzzyaaad98)))*dble(xyzzyaaao98)
xyzzyaaaf98(1)=xyzzyaaaf98(1)-xyzzyaagg1(1,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaaf98(2)=xyzzyaaaf98(2)-xyzzyaagg1(2,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaaf98(3)=xyzzyaaaf98(3)-xyzzyaagg1(3,xyzzyaaad98)*xyzzyaaaj98
enddo
grad_p(1)=grad_p(1)+xyzzyaaaf98(1)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(2)=grad_p(2)+xyzzyaaaf98(2)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(3)=grad_p(3)+xyzzyaaaf98(3)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
enddo
elseif(val.and.fd.and..not.sd)then
value_p=0.d0
grad_p(1)=0.d0
grad_p(2)=0.d0
grad_p(3)=0.d0
xyzzyaaad98=0
do xyzzyaaab98=1,xyzzyaagq1
xyzzyaaae98=0.d0
xyzzyaaaf98(1)=0.d0
xyzzyaaaf98(2)=0.d0
xyzzyaaaf98(3)=0.d0
do xyzzyaaac98=1,xyzzyaagt1(xyzzyaaab98)
xyzzyaaad98=xyzzyaaad98+1
xyzzyaaak98=xyzzyaagk1(xyzzyaaga1(1,xyzzyaaad98))*xyzzyaagl1(xyzzyaaga&
&1(2,xyzzyaaad98))*xyzzyaagm1(xyzzyaaga1(3,xyzzyaaad98))
xyzzyaaaj98=aimag(xyzzyaaak98)
xyzzyaaae98=xyzzyaaae98+dble(xyzzyaaak98)
xyzzyaaaf98(1)=xyzzyaaaf98(1)-xyzzyaagg1(1,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaaf98(2)=xyzzyaaaf98(2)-xyzzyaagg1(2,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaaf98(3)=xyzzyaaaf98(3)-xyzzyaagg1(3,xyzzyaaad98)*xyzzyaaaj98
enddo
value_p=value_p+xyzzyaaae98*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(1)=grad_p(1)+xyzzyaaaf98(1)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(2)=grad_p(2)+xyzzyaaaf98(2)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(3)=grad_p(3)+xyzzyaaaf98(3)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
enddo
elseif(.not.val.and.sd)then
grad_p(1)=0.d0
grad_p(2)=0.d0
grad_p(3)=0.d0
lap_p=0.d0
xyzzyaaad98=0
do xyzzyaaab98=1,xyzzyaagq1
xyzzyaaaf98(1)=0.d0
xyzzyaaaf98(2)=0.d0
xyzzyaaaf98(3)=0.d0
xyzzyaaag98=0.d0
do xyzzyaaac98=1,xyzzyaagt1(xyzzyaaab98)
xyzzyaaad98=xyzzyaaad98+1
xyzzyaaak98=xyzzyaagk1(xyzzyaaga1(1,xyzzyaaad98))*xyzzyaagl1(xyzzyaaga&
&1(2,xyzzyaaad98))*xyzzyaagm1(xyzzyaaga1(3,xyzzyaaad98))
xyzzyaaaj98=aimag(xyzzyaaak98)
xyzzyaaaf98(1)=xyzzyaaaf98(1)-xyzzyaagg1(1,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaaf98(2)=xyzzyaaaf98(2)-xyzzyaagg1(2,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaaf98(3)=xyzzyaaaf98(3)-xyzzyaagg1(3,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaag98=xyzzyaaag98-xyzzyaagh1(xyzzyaaad98)*dble(xyzzyaaak98)
enddo
grad_p(1)=grad_p(1)+xyzzyaaaf98(1)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(2)=grad_p(2)+xyzzyaaaf98(2)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(3)=grad_p(3)+xyzzyaaaf98(3)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
lap_p=lap_p+xyzzyaaag98*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
enddo
elseif(val.and.fd.and.sd)then
value_p=0.d0
grad_p(1)=0.d0
grad_p(2)=0.d0
grad_p(3)=0.d0
lap_p=0.d0
xyzzyaaad98=0
do xyzzyaaab98=1,xyzzyaagq1
xyzzyaaae98=0.d0
xyzzyaaaf98(1)=0.d0
xyzzyaaaf98(2)=0.d0
xyzzyaaaf98(3)=0.d0
xyzzyaaag98=0.d0
do xyzzyaaac98=1,xyzzyaagt1(xyzzyaaab98)
xyzzyaaad98=xyzzyaaad98+1
xyzzyaaak98=xyzzyaagk1(xyzzyaaga1(1,xyzzyaaad98))*xyzzyaagl1(xyzzyaaga&
&1(2,xyzzyaaad98))*xyzzyaagm1(xyzzyaaga1(3,xyzzyaaad98))
xyzzyaaai98=dble(xyzzyaaak98)
xyzzyaaaj98=aimag(xyzzyaaak98)
xyzzyaaae98=xyzzyaaae98+xyzzyaaai98
xyzzyaaaf98(1)=xyzzyaaaf98(1)-xyzzyaagg1(1,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaaf98(2)=xyzzyaaaf98(2)-xyzzyaagg1(2,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaaf98(3)=xyzzyaaaf98(3)-xyzzyaagg1(3,xyzzyaaad98)*xyzzyaaaj98
xyzzyaaag98=xyzzyaaag98-xyzzyaagh1(xyzzyaaad98)*xyzzyaaai98
enddo
value_p=value_p+xyzzyaaae98*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(1)=grad_p(1)+xyzzyaaaf98(1)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(2)=grad_p(2)+xyzzyaaaf98(2)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
grad_p(3)=grad_p(3)+xyzzyaaaf98(3)*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
lap_p=lap_p+xyzzyaaag98*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
enddo
elseif(val.and..not.fd.and.sd)then
value_p=0.d0
lap_p=0.d0
xyzzyaaad98=0
do xyzzyaaab98=1,xyzzyaagq1
xyzzyaaae98=0.d0
xyzzyaaag98=0.d0
do xyzzyaaac98=1,xyzzyaagt1(xyzzyaaab98)
xyzzyaaad98=xyzzyaaad98+1
xyzzyaaai98=dble(xyzzyaagk1(xyzzyaaga1(1,xyzzyaaad98))*xyzzyaagl1(xyzz&
&yaaga1(2,xyzzyaaad98))*xyzzyaagm1(xyzzyaaga1(3,xyzzyaaad98)))
xyzzyaaae98=xyzzyaaae98+xyzzyaaai98
xyzzyaaag98=xyzzyaaag98-xyzzyaagh1(xyzzyaaad98)*xyzzyaaai98
enddo
value_p=value_p+xyzzyaaae98*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
lap_p=lap_p+xyzzyaaag98*xyzzyaaah1(xyzzyaaab98,xyzzyaaaa98)
enddo
else
call errstop('COMPUTE_P','Combination of val ('//l2s(val)//'), fd ('//&
&l2s(fd)//') and sd ('//l2s(sd)//') not yet catered for.')
endif
end subroutine xyzzyaaqd1
subroutine xyzzyaaqe1(r_i,ispin,val,fd,sd,value_q,grad_q,lap_q)
implicit none
integer,intent(in) :: ispin
real(dp),intent(in) :: r_i(3)
real(dp),intent(out) :: value_q,grad_q(3),lap_q
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa99,xyzzyaaab99,xyzzyaaac99
real(dp) xyzzyaaad99,xyzzyaaae99(3),xyzzyaaaf99,xyzzyaaag99,xyzzyaaah9&
&9,xyzzyaaai99
complex(dp) xyzzyaaaj99,xyzzyaaak99,xyzzyaaal99,xyzzyaaam99,xyzzyaaan9&
&9
xyzzyaaag99=pb1(1)*r_i(1)+pb1(2)*r_i(2)+pb1(3)*r_i(3)
xyzzyaaak99=cmplx(cos(xyzzyaaag99),sin(xyzzyaaag99),dp)
xyzzyaaag99=pb2(1)*r_i(1)+pb2(2)*r_i(2)+pb2(3)*r_i(3)
xyzzyaaal99=cmplx(cos(xyzzyaaag99),sin(xyzzyaaag99),dp)
xyzzyaaag99=pb3(1)*r_i(1)+pb3(2)*r_i(2)+pb3(3)*r_i(3)
xyzzyaaam99=cmplx(cos(xyzzyaaag99),sin(xyzzyaaag99),dp)
xyzzyaagn1(0)=c_one
do xyzzyaaaa99=1,xyzzyaagf1(1)
xyzzyaagn1(xyzzyaaaa99)=xyzzyaagn1(xyzzyaaaa99-1)*xyzzyaaak99
enddo
xyzzyaago1(0)=c_one
do xyzzyaaaa99=1,xyzzyaagf1(2)
xyzzyaago1(xyzzyaaaa99)=xyzzyaago1(xyzzyaaaa99-1)*xyzzyaaal99
xyzzyaago1(-xyzzyaaaa99)=conjg(xyzzyaago1(xyzzyaaaa99))
enddo
xyzzyaagp1(0)=c_one
do xyzzyaaaa99=1,xyzzyaagf1(3)
xyzzyaagp1(xyzzyaaaa99)=xyzzyaagp1(xyzzyaaaa99-1)*xyzzyaaam99
xyzzyaagp1(-xyzzyaaaa99)=conjg(xyzzyaagp1(xyzzyaaaa99))
enddo
value_q=0.d0
lap_q=0.d0
grad_q(1)=0.d0
grad_q(2)=0.d0
grad_q(3)=0.d0
if(.not.fd.and..not.sd)then
xyzzyaaac99=0
do xyzzyaaaa99=1,xyzzyaags1
xyzzyaaad99=0.d0
do xyzzyaaab99=1,xyzzyaagu1(xyzzyaaaa99)
xyzzyaaac99=xyzzyaaac99+1
xyzzyaaan99=xyzzyaago1(xyzzyaagc1(2,xyzzyaaac99))*xyzzyaagp1(xyzzyaagc&
&1(3,xyzzyaaac99))
xyzzyaaad99=xyzzyaaad99+dble(xyzzyaagn1(xyzzyaagc1(1,xyzzyaaac99)))*db&
&le(xyzzyaaan99)-aimag(xyzzyaagn1(xyzzyaagc1(1,xyzzyaaac99)))*aimag(xy&
&zzyaaan99)
enddo
value_q=value_q+xyzzyaaad99*xyzzyaaai1(xyzzyaaaa99,ispin)
enddo
elseif(.not.val.and..not.sd)then
xyzzyaaac99=0
do xyzzyaaaa99=1,xyzzyaags1
xyzzyaaae99(1)=0.d0
xyzzyaaae99(2)=0.d0
xyzzyaaae99(3)=0.d0
do xyzzyaaab99=1,xyzzyaagu1(xyzzyaaaa99)
xyzzyaaac99=xyzzyaaac99+1
xyzzyaaan99=xyzzyaago1(xyzzyaagc1(2,xyzzyaaac99))*xyzzyaagp1(xyzzyaagc&
&1(3,xyzzyaaac99))
xyzzyaaai99=dble(xyzzyaagn1(xyzzyaagc1(1,xyzzyaaac99)))*aimag(xyzzyaaa&
&n99)+aimag(xyzzyaagn1(xyzzyaagc1(1,xyzzyaaac99)))*dble(xyzzyaaan99)
xyzzyaaae99(1)=xyzzyaaae99(1)-xyzzyaagi1(1,xyzzyaaac99)*xyzzyaaai99
xyzzyaaae99(2)=xyzzyaaae99(2)-xyzzyaagi1(2,xyzzyaaac99)*xyzzyaaai99
xyzzyaaae99(3)=xyzzyaaae99(3)-xyzzyaagi1(3,xyzzyaaac99)*xyzzyaaai99
enddo
grad_q(1)=grad_q(1)+xyzzyaaae99(1)*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(2)=grad_q(2)+xyzzyaaae99(2)*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(3)=grad_q(3)+xyzzyaaae99(3)*xyzzyaaai1(xyzzyaaaa99,ispin)
enddo
elseif(.not.val.and..not.fd)then
xyzzyaaac99=0
do xyzzyaaaa99=1,xyzzyaags1
xyzzyaaaf99=0.d0
do xyzzyaaab99=1,xyzzyaagu1(xyzzyaaaa99)
xyzzyaaac99=xyzzyaaac99+1
xyzzyaaan99=xyzzyaago1(xyzzyaagc1(2,xyzzyaaac99))*xyzzyaagp1(xyzzyaagc&
&1(3,xyzzyaaac99))
xyzzyaaaf99=xyzzyaaaf99-xyzzyaagj1(xyzzyaaac99)*dble(xyzzyaagn1(xyzzya&
&agc1(1,xyzzyaaac99)))*dble(xyzzyaaan99)-aimag(xyzzyaagn1(xyzzyaagc1(1&
&,xyzzyaaac99)))*aimag(xyzzyaaan99)
enddo
lap_q=lap_q+xyzzyaaaf99*xyzzyaaai1(xyzzyaaaa99,ispin)
enddo
elseif(val.and.fd.and..not.sd)then
xyzzyaaac99=0
do xyzzyaaaa99=1,xyzzyaags1
xyzzyaaad99=0.d0
xyzzyaaae99(1)=0.d0
xyzzyaaae99(2)=0.d0
xyzzyaaae99(3)=0.d0
do xyzzyaaab99=1,xyzzyaagu1(xyzzyaaaa99)
xyzzyaaac99=xyzzyaaac99+1
xyzzyaaaj99=xyzzyaagn1(xyzzyaagc1(1,xyzzyaaac99))*xyzzyaago1(xyzzyaagc&
&1(2,xyzzyaaac99))*xyzzyaagp1(xyzzyaagc1(3,xyzzyaaac99))
xyzzyaaai99=aimag(xyzzyaaaj99)
xyzzyaaad99=xyzzyaaad99+dble(xyzzyaaaj99)
xyzzyaaae99(1)=xyzzyaaae99(1)-xyzzyaagi1(1,xyzzyaaac99)*xyzzyaaai99
xyzzyaaae99(2)=xyzzyaaae99(2)-xyzzyaagi1(2,xyzzyaaac99)*xyzzyaaai99
xyzzyaaae99(3)=xyzzyaaae99(3)-xyzzyaagi1(3,xyzzyaaac99)*xyzzyaaai99
enddo
value_q=value_q+xyzzyaaad99*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(1)=grad_q(1)+xyzzyaaae99(1)*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(2)=grad_q(2)+xyzzyaaae99(2)*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(3)=grad_q(3)+xyzzyaaae99(3)*xyzzyaaai1(xyzzyaaaa99,ispin)
enddo
elseif(sd.and.fd.and..not.val)then
xyzzyaaac99=0
do xyzzyaaaa99=1,xyzzyaags1
xyzzyaaae99(1)=0.d0
xyzzyaaae99(2)=0.d0
xyzzyaaae99(3)=0.d0
xyzzyaaaf99=0.d0
do xyzzyaaab99=1,xyzzyaagu1(xyzzyaaaa99)
xyzzyaaac99=xyzzyaaac99+1
xyzzyaaaj99=xyzzyaagn1(xyzzyaagc1(1,xyzzyaaac99))*xyzzyaago1(xyzzyaagc&
&1(2,xyzzyaaac99))*xyzzyaagp1(xyzzyaagc1(3,xyzzyaaac99))
xyzzyaaai99=aimag(xyzzyaaaj99)
xyzzyaaaf99=xyzzyaaaf99-xyzzyaagj1(xyzzyaaac99)*dble(xyzzyaaaj99)
xyzzyaaae99(1)=xyzzyaaae99(1)-xyzzyaagi1(1,xyzzyaaac99)*xyzzyaaai99
xyzzyaaae99(2)=xyzzyaaae99(2)-xyzzyaagi1(2,xyzzyaaac99)*xyzzyaaai99
xyzzyaaae99(3)=xyzzyaaae99(3)-xyzzyaagi1(3,xyzzyaaac99)*xyzzyaaai99
enddo
grad_q(1)=grad_q(1)+xyzzyaaae99(1)*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(2)=grad_q(2)+xyzzyaaae99(2)*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(3)=grad_q(3)+xyzzyaaae99(3)*xyzzyaaai1(xyzzyaaaa99,ispin)
lap_q=lap_q+xyzzyaaaf99*xyzzyaaai1(xyzzyaaaa99,ispin)
enddo
elseif(val.and.fd.and.sd)then
xyzzyaaac99=0
do xyzzyaaaa99=1,xyzzyaags1
xyzzyaaad99=0.d0
xyzzyaaae99(1)=0.d0
xyzzyaaae99(2)=0.d0
xyzzyaaae99(3)=0.d0
xyzzyaaaf99=0.d0
do xyzzyaaab99=1,xyzzyaagu1(xyzzyaaaa99)
xyzzyaaac99=xyzzyaaac99+1
xyzzyaaaj99=xyzzyaagn1(xyzzyaagc1(1,xyzzyaaac99))*xyzzyaago1(xyzzyaagc&
&1(2,xyzzyaaac99))*xyzzyaagp1(xyzzyaagc1(3,xyzzyaaac99))
xyzzyaaah99=dble(xyzzyaaaj99)
xyzzyaaai99=aimag(xyzzyaaaj99)
xyzzyaaad99=xyzzyaaad99+xyzzyaaah99
xyzzyaaae99(1)=xyzzyaaae99(1)-xyzzyaagi1(1,xyzzyaaac99)*xyzzyaaai99
xyzzyaaae99(2)=xyzzyaaae99(2)-xyzzyaagi1(2,xyzzyaaac99)*xyzzyaaai99
xyzzyaaae99(3)=xyzzyaaae99(3)-xyzzyaagi1(3,xyzzyaaac99)*xyzzyaaai99
xyzzyaaaf99=xyzzyaaaf99-xyzzyaagj1(xyzzyaaac99)*xyzzyaaah99
enddo
value_q=value_q+xyzzyaaad99*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(1)=grad_q(1)+xyzzyaaae99(1)*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(2)=grad_q(2)+xyzzyaaae99(2)*xyzzyaaai1(xyzzyaaaa99,ispin)
grad_q(3)=grad_q(3)+xyzzyaaae99(3)*xyzzyaaai1(xyzzyaaaa99,ispin)
lap_q=lap_q+xyzzyaaaf99*xyzzyaaai1(xyzzyaaaa99,ispin)
enddo
elseif(val.and..not.fd.and.sd)then
xyzzyaaac99=0
do xyzzyaaaa99=1,xyzzyaags1
xyzzyaaad99=0.d0
xyzzyaaaf99=0.d0
do xyzzyaaab99=1,xyzzyaagu1(xyzzyaaaa99)
xyzzyaaac99=xyzzyaaac99+1
xyzzyaaah99=dble(xyzzyaagn1(xyzzyaagc1(1,xyzzyaaac99))*xyzzyaago1(xyzz&
&yaagc1(2,xyzzyaaac99))*xyzzyaagp1(xyzzyaagc1(3,xyzzyaaac99)))
xyzzyaaad99=xyzzyaaad99+xyzzyaaah99
xyzzyaaaf99=xyzzyaaaf99-xyzzyaagj1(xyzzyaaac99)*xyzzyaaah99
enddo
value_q=value_q+xyzzyaaad99*xyzzyaaai1(xyzzyaaaa99,ispin)
lap_q=lap_q+xyzzyaaaf99*xyzzyaaai1(xyzzyaaaa99,ispin)
enddo
else
call errstop('COMPUTE_Q','Combination of val ('//l2s(val)//'), fd ('//&
&l2s(fd)//') and sd ('//l2s(sd)//') not yet catered for.')
endif
end subroutine xyzzyaaqe1
subroutine xyzzyaaqf1(val,fd,sd,jspin,eevecs,val_jas,grad_jas,lap_jas)
use slaarnabt,only : exp_limit
implicit none
integer,intent(in) :: jspin
real(dp),intent(out) :: val_jas,grad_jas(3),lap_jas
real(dp),intent(in) :: eevecs(4,netot)
logical,intent(in) :: val,fd,sd
real(dp) xyzzyaaaa100,xyzzyaaab100,xyzzyaaac100,xyzzyaaad100,xyzzyaaae&
&100,xyzzyaaaf100,xyzzyaaag100,xyzzyaaah100,xyzzyaaai100
real(dp) xyzzyaaaj100,xyzzyaaak100,xyzzyaaal100,xyzzyaaam100,xyzzyaaan&
&100,xyzzyaaao100
real(dp) xyzzyaaap100,xyzzyaaaq100,xyzzyaaar100,xyzzyaaas100,xyzzyaaat&
&100,xyzzyaaau100
real(dp) xyzzyaaav100,xyzzyaaaw100,xyzzyaaax100,xyzzyaaay100,xyzzyaaaz&
&100,xyzzyaaba100
real(dp) xyzzyaabb100,xyzzyaabc100,xyzzyaabd100,xyzzyaabe100,xyzzyaabf&
&100,xyzzyaabg100,xyzzyaabh100,xyzzyaabi100,xyzzyaabj100,xyzzyaabk100
real(dp) xyzzyaabl100,xyzzyaabm100,xyzzyaabn100,xyzzyaabo100,xyzzyaabp&
&100,xyzzyaabq100,xyzzyaabr100,xyzzyaabs100,xyzzyaabt100,xyzzyaabu100,&
&xyzzyaabv100,xyzzyaabw100,xyzzyaabx100,xyzzyaaby100,xyzzyaabz100,xyzz&
&yaaca100,xyzzyaacb100,xyzzyaacc100,xyzzyaacd100,xyzzyaace100,xyzzyaac&
&f100,xyzzyaacg100,xyzzyaach100(2),xyzzyaaci100(2),xyzzyaacj100,xyzzya&
&ack100,xyzzyaacl100,xyzzyaacm100,xyzzyaacn100(2),xyzzyaaco100,xyzzyaa&
&cp100
xyzzyaaaa100=xyzzyaaih1(1)
xyzzyaaab100=xyzzyaaih1(2)
xyzzyaaac100=xyzzyaaih1(3)
xyzzyaaad100=xyzzyaaih1(4)
xyzzyaaae100=xyzzyaaih1(5)
xyzzyaaaf100=xyzzyaaih1(6)
xyzzyaaag100=xyzzyaaih1(7)
xyzzyaaah100=xyzzyaaih1(8)
xyzzyaaai100=xyzzyaaih1(9)
if(jspin==1)then
xyzzyaaaj100=eevecs(1,2)
xyzzyaaap100=eevecs(2,2)
xyzzyaaav100=eevecs(4,2)
xyzzyaaal100=eevecs(1,3)
xyzzyaaar100=eevecs(2,3)
xyzzyaaax100=sqrt(xyzzyaaal100*xyzzyaaal100+xyzzyaaar100*xyzzyaaar100)
xyzzyaaam100=eevecs(1,4)
xyzzyaaas100=eevecs(2,4)
xyzzyaaaz100=sqrt(xyzzyaaam100*xyzzyaaam100+xyzzyaaas100*xyzzyaaas100)
xyzzyaaan100=xyzzyaaal100-xyzzyaaaj100
xyzzyaaat100=xyzzyaaar100-xyzzyaaap100
xyzzyaaay100=sqrt(xyzzyaaan100*xyzzyaaan100+xyzzyaaat100*xyzzyaaat100)
xyzzyaaao100=xyzzyaaam100-xyzzyaaaj100
xyzzyaaau100=xyzzyaaas100-xyzzyaaap100
xyzzyaaba100=sqrt(xyzzyaaao100*xyzzyaaao100+xyzzyaaau100*xyzzyaaau100)
xyzzyaabb100=1.d0/(1.d0+xyzzyaaab100*xyzzyaaav100)
if(val)then
xyzzyaaak100=xyzzyaaam100-xyzzyaaal100
xyzzyaaaq100=xyzzyaaas100-xyzzyaaar100
xyzzyaaaw100=sqrt(xyzzyaaak100*xyzzyaaak100+xyzzyaaaq100*xyzzyaaaq100)
xyzzyaabc100=1.d0/(1.d0+xyzzyaaad100*xyzzyaaaw100)
endif
elseif(jspin==2)then
xyzzyaaaj100=-eevecs(1,1)
xyzzyaaap100=-eevecs(2,1)
xyzzyaaav100=eevecs(4,1)
xyzzyaaan100=eevecs(1,3)
xyzzyaaat100=eevecs(2,3)
xyzzyaaay100=sqrt(xyzzyaaan100*xyzzyaaan100+xyzzyaaat100*xyzzyaaat100)
xyzzyaaao100=eevecs(1,4)
xyzzyaaau100=eevecs(2,4)
xyzzyaaba100=sqrt(xyzzyaaao100*xyzzyaaao100+xyzzyaaau100*xyzzyaaau100)
xyzzyaaal100=xyzzyaaan100+xyzzyaaaj100
xyzzyaaar100=xyzzyaaat100+xyzzyaaap100
xyzzyaaax100=sqrt(xyzzyaaal100*xyzzyaaal100+xyzzyaaar100*xyzzyaaar100)
xyzzyaaam100=xyzzyaaao100+xyzzyaaaj100
xyzzyaaas100=xyzzyaaau100+xyzzyaaap100
xyzzyaaaz100=sqrt(xyzzyaaam100*xyzzyaaam100+xyzzyaaas100*xyzzyaaas100)
xyzzyaabb100=1.d0/(1.d0+xyzzyaaab100*xyzzyaaav100)
if(val)then
xyzzyaaak100=xyzzyaaao100-xyzzyaaan100
xyzzyaaaq100=xyzzyaaau100-xyzzyaaat100
xyzzyaaaw100=sqrt(xyzzyaaak100*xyzzyaaak100+xyzzyaaaq100*xyzzyaaaq100)
xyzzyaabc100=1.d0/(1.d0+xyzzyaaad100*xyzzyaaaw100)
endif
elseif(jspin==3)then
xyzzyaaak100=eevecs(1,4)
xyzzyaaaq100=eevecs(2,4)
xyzzyaaaw100=eevecs(4,4)
xyzzyaaal100=-eevecs(1,1)
xyzzyaaar100=-eevecs(2,1)
xyzzyaaax100=sqrt(xyzzyaaal100*xyzzyaaal100+xyzzyaaar100*xyzzyaaar100)
xyzzyaaan100=-eevecs(1,2)
xyzzyaaat100=-eevecs(2,2)
xyzzyaaay100=sqrt(xyzzyaaan100*xyzzyaaan100+xyzzyaaat100*xyzzyaaat100)
xyzzyaaam100=xyzzyaaak100+xyzzyaaal100
xyzzyaaas100=xyzzyaaaq100+xyzzyaaar100
xyzzyaaaz100=sqrt(xyzzyaaam100*xyzzyaaam100+xyzzyaaas100*xyzzyaaas100)
xyzzyaaao100=xyzzyaaak100+xyzzyaaan100
xyzzyaaau100=xyzzyaaaq100+xyzzyaaat100
xyzzyaaba100=sqrt(xyzzyaaao100*xyzzyaaao100+xyzzyaaau100*xyzzyaaau100)
xyzzyaabc100=1.d0/(1.d0+xyzzyaaad100*xyzzyaaaw100)
if(val)then
xyzzyaaaj100=xyzzyaaal100-xyzzyaaan100
xyzzyaaap100=xyzzyaaar100-xyzzyaaat100
xyzzyaaav100=sqrt(xyzzyaaaj100*xyzzyaaaj100+xyzzyaaap100*xyzzyaaap100)
xyzzyaabb100=1.d0/(1.d0+xyzzyaaab100*xyzzyaaav100)
endif
else
xyzzyaaak100=-eevecs(1,3)
xyzzyaaaq100=-eevecs(2,3)
xyzzyaaaw100=eevecs(4,3)
xyzzyaaam100=-eevecs(1,1)
xyzzyaaas100=-eevecs(2,1)
xyzzyaaaz100=sqrt(xyzzyaaam100*xyzzyaaam100+xyzzyaaas100*xyzzyaaas100)
xyzzyaaao100=-eevecs(1,2)
xyzzyaaau100=-eevecs(2,2)
xyzzyaaba100=sqrt(xyzzyaaao100*xyzzyaaao100+xyzzyaaau100*xyzzyaaau100)
xyzzyaaal100=xyzzyaaam100-xyzzyaaak100
xyzzyaaar100=xyzzyaaas100-xyzzyaaaq100
xyzzyaaax100=sqrt(xyzzyaaal100*xyzzyaaal100+xyzzyaaar100*xyzzyaaar100)
xyzzyaaan100=xyzzyaaao100-xyzzyaaak100
xyzzyaaat100=xyzzyaaau100-xyzzyaaaq100
xyzzyaaay100=sqrt(xyzzyaaan100*xyzzyaaan100+xyzzyaaat100*xyzzyaaat100)
xyzzyaabc100=1.d0/(1.d0+xyzzyaaad100*xyzzyaaaw100)
if(val)then
xyzzyaaaj100=xyzzyaaam100-xyzzyaaao100
xyzzyaaap100=xyzzyaaas100-xyzzyaaau100
xyzzyaaav100=sqrt(xyzzyaaaj100*xyzzyaaaj100+xyzzyaaap100*xyzzyaaap100)
xyzzyaabb100=1.d0/(1.d0+xyzzyaaab100*xyzzyaaav100)
endif
endif
xyzzyaabd100=1.d0/(1.d0+xyzzyaaag100*xyzzyaaax100)
xyzzyaabh100=1.d0/(1.d0+xyzzyaaai100*xyzzyaaax100)
xyzzyaabe100=1.d0/(1.d0+xyzzyaaag100*xyzzyaaaz100)
xyzzyaabi100=1.d0/(1.d0+xyzzyaaai100*xyzzyaaaz100)
xyzzyaabf100=1.d0/(1.d0+xyzzyaaag100*xyzzyaaay100)
xyzzyaabj100=1.d0/(1.d0+xyzzyaaai100*xyzzyaaay100)
xyzzyaabg100=1.d0/(1.d0+xyzzyaaag100*xyzzyaaba100)
xyzzyaabk100=1.d0/(1.d0+xyzzyaaai100*xyzzyaaba100)
xyzzyaaco100=(xyzzyaaae100+xyzzyaaaf100*xyzzyaaax100)*xyzzyaaax100*xyz&
&zyaabd100+(xyzzyaaae100+xyzzyaaah100*xyzzyaaaz100)*xyzzyaaaz100*xyzzy&
&aabi100+(xyzzyaaae100+xyzzyaaah100*xyzzyaaay100)*xyzzyaaay100*xyzzyaa&
&bj100+(xyzzyaaae100+xyzzyaaaf100*xyzzyaaba100)*xyzzyaaba100*xyzzyaabg&
&100
xyzzyaabl100=exp_limit(xyzzyaaco100)
xyzzyaaco100=(xyzzyaaae100+xyzzyaaah100*xyzzyaaax100)*xyzzyaaax100*xyz&
&zyaabh100+(xyzzyaaae100+xyzzyaaaf100*xyzzyaaaz100)*xyzzyaaaz100*xyzzy&
&aabe100+(xyzzyaaae100+xyzzyaaaf100*xyzzyaaay100)*xyzzyaaay100*xyzzyaa&
&bf100+(xyzzyaaae100+xyzzyaaah100*xyzzyaaba100)*xyzzyaaba100*xyzzyaabk&
&100
xyzzyaabm100=exp_limit(xyzzyaaco100)
xyzzyaacp100=xyzzyaabl100+xyzzyaabm100
if(val)val_jas=xyzzyaaaa100*xyzzyaaav100*xyzzyaabb100+xyzzyaaac100*xyz&
&zyaaaw100*xyzzyaabc100+log(xyzzyaacp100)
if(fd.or.sd)then
if(jspin==1)then
xyzzyaabn100=xyzzyaaaa100*xyzzyaabb100*(1.d0/xyzzyaaav100-xyzzyaaab100&
&*xyzzyaabb100)
xyzzyaabr100=xyzzyaabd100*(xyzzyaaae100/xyzzyaaax100+2*xyzzyaaaf100-xy&
&zzyaaag100*xyzzyaabd100*(xyzzyaaae100+xyzzyaaaf100*xyzzyaaax100))
xyzzyaabs100=xyzzyaabi100*(xyzzyaaae100/xyzzyaaaz100+2*xyzzyaaah100-xy&
&zzyaaai100*xyzzyaabi100*(xyzzyaaae100+xyzzyaaah100*xyzzyaaaz100))
xyzzyaabv100=xyzzyaabh100*(xyzzyaaae100/xyzzyaaax100+2*xyzzyaaah100-xy&
&zzyaaai100*xyzzyaabh100*(xyzzyaaae100+xyzzyaaah100*xyzzyaaax100))
xyzzyaabw100=xyzzyaabe100*(xyzzyaaae100/xyzzyaaaz100+2*xyzzyaaaf100-xy&
&zzyaaag100*xyzzyaabe100*(xyzzyaaae100+xyzzyaaaf100*xyzzyaaaz100))
if(fd)then
grad_jas(1)=xyzzyaabn100*xyzzyaaaj100
grad_jas(2)=xyzzyaabn100*xyzzyaaap100
endif
xyzzyaach100(1)=xyzzyaabr100*xyzzyaaal100+xyzzyaabs100*xyzzyaaam100
xyzzyaach100(2)=xyzzyaabr100*xyzzyaaar100+xyzzyaabs100*xyzzyaaas100
xyzzyaaci100(1)=xyzzyaabv100*xyzzyaaal100+xyzzyaabw100*xyzzyaaam100
xyzzyaaci100(2)=xyzzyaabv100*xyzzyaaar100+xyzzyaabw100*xyzzyaaas100
if(sd)then
xyzzyaabp100=2*xyzzyaaaa100*xyzzyaaab100*xyzzyaabb100**2*(xyzzyaaab100&
&*xyzzyaabb100*xyzzyaaav100-1.d0)
xyzzyaabz100=2*xyzzyaabd100*(xyzzyaaaf100-xyzzyaaag100*xyzzyaabd100*(x&
&yzzyaaae100+xyzzyaaax100*(2*xyzzyaaaf100-xyzzyaaag100*xyzzyaabd100*(x&
&yzzyaaae100+xyzzyaaaf100*xyzzyaaax100))))
xyzzyaaca100=2*xyzzyaabi100*(xyzzyaaah100-xyzzyaaai100*xyzzyaabi100*(x&
&yzzyaaae100+xyzzyaaaz100*(2*xyzzyaaah100-xyzzyaaai100*xyzzyaabi100*(x&
&yzzyaaae100+xyzzyaaah100*xyzzyaaaz100))))
xyzzyaacd100=2*xyzzyaabh100*(xyzzyaaah100-xyzzyaaai100*xyzzyaabh100*(x&
&yzzyaaae100+xyzzyaaax100*(2*xyzzyaaah100-xyzzyaaai100*xyzzyaabh100*(x&
&yzzyaaae100+xyzzyaaah100*xyzzyaaax100))))
xyzzyaace100=2*xyzzyaabe100*(xyzzyaaaf100-xyzzyaaag100*xyzzyaabe100*(x&
&yzzyaaae100+xyzzyaaaz100*(2*xyzzyaaaf100-xyzzyaaag100*xyzzyaabe100*(x&
&yzzyaaae100+xyzzyaaaf100*xyzzyaaaz100))))
lap_jas=xyzzyaabp100+xyzzyaabn100
xyzzyaacj100=xyzzyaabz100+xyzzyaabr100+xyzzyaaca100+xyzzyaabs100
xyzzyaack100=xyzzyaacd100+xyzzyaabv100+xyzzyaace100+xyzzyaabw100
endif
elseif(jspin==2)then
xyzzyaabn100=xyzzyaaaa100*xyzzyaabb100*(1.d0/xyzzyaaav100-xyzzyaaab100&
&*xyzzyaabb100)
xyzzyaabt100=xyzzyaabj100*(xyzzyaaae100/xyzzyaaay100+2*xyzzyaaah100-xy&
&zzyaaai100*xyzzyaabj100*(xyzzyaaae100+xyzzyaaah100*xyzzyaaay100))
xyzzyaabu100=xyzzyaabg100*(xyzzyaaae100/xyzzyaaba100+2*xyzzyaaaf100-xy&
&zzyaaag100*xyzzyaabg100*(xyzzyaaae100+xyzzyaaaf100*xyzzyaaba100))
xyzzyaabx100=xyzzyaabf100*(xyzzyaaae100/xyzzyaaay100+2*xyzzyaaaf100-xy&
&zzyaaag100*xyzzyaabf100*(xyzzyaaae100+xyzzyaaaf100*xyzzyaaay100))
xyzzyaaby100=xyzzyaabk100*(xyzzyaaae100/xyzzyaaba100+2*xyzzyaaah100-xy&
&zzyaaai100*xyzzyaabk100*(xyzzyaaae100+xyzzyaaah100*xyzzyaaba100))
if(fd)then
grad_jas(1)=-xyzzyaabn100*xyzzyaaaj100
grad_jas(2)=-xyzzyaabn100*xyzzyaaap100
endif
xyzzyaach100(1)=xyzzyaabt100*xyzzyaaan100+xyzzyaabu100*xyzzyaaao100
xyzzyaach100(2)=xyzzyaabt100*xyzzyaaat100+xyzzyaabu100*xyzzyaaau100
xyzzyaaci100(1)=xyzzyaabx100*xyzzyaaan100+xyzzyaaby100*xyzzyaaao100
xyzzyaaci100(2)=xyzzyaabx100*xyzzyaaat100+xyzzyaaby100*xyzzyaaau100
if(sd)then
xyzzyaabp100=2*xyzzyaaaa100*xyzzyaaab100*xyzzyaabb100**2*(xyzzyaaab100&
&*xyzzyaabb100*xyzzyaaav100-1.d0)
xyzzyaacb100=2*xyzzyaabj100*(xyzzyaaah100-xyzzyaaai100*xyzzyaabj100*(x&
&yzzyaaae100+xyzzyaaay100*(2*xyzzyaaah100-xyzzyaaai100*xyzzyaabj100*(x&
&yzzyaaae100+xyzzyaaah100*xyzzyaaay100))))
xyzzyaacc100=2*xyzzyaabg100*(xyzzyaaaf100-xyzzyaaag100*xyzzyaabg100*(x&
&yzzyaaae100+xyzzyaaba100*(2*xyzzyaaaf100-xyzzyaaag100*xyzzyaabg100*(x&
&yzzyaaae100+xyzzyaaaf100*xyzzyaaba100))))
xyzzyaacf100=2*xyzzyaabf100*(xyzzyaaaf100-xyzzyaaag100*xyzzyaabf100*(x&
&yzzyaaae100+xyzzyaaay100*(2*xyzzyaaaf100-xyzzyaaag100*xyzzyaabf100*(x&
&yzzyaaae100+xyzzyaaaf100*xyzzyaaay100))))
xyzzyaacg100=2*xyzzyaabk100*(xyzzyaaah100-xyzzyaaai100*xyzzyaabk100*(x&
&yzzyaaae100+xyzzyaaba100*(2*xyzzyaaah100-xyzzyaaai100*xyzzyaabk100*(x&
&yzzyaaae100+xyzzyaaah100*xyzzyaaba100))))
lap_jas=xyzzyaabp100+xyzzyaabn100
xyzzyaacj100=xyzzyaacb100+xyzzyaabt100+xyzzyaacc100+xyzzyaabu100
xyzzyaack100=xyzzyaacf100+xyzzyaabx100+xyzzyaacg100+xyzzyaaby100
endif
elseif(jspin==3)then
xyzzyaabo100=xyzzyaaac100*xyzzyaabc100*(1.d0/xyzzyaaaw100-xyzzyaaad100&
&*xyzzyaabc100)
xyzzyaabr100=xyzzyaabd100*(xyzzyaaae100/xyzzyaaax100+2*xyzzyaaaf100-xy&
&zzyaaag100*xyzzyaabd100*(xyzzyaaae100+xyzzyaaaf100*xyzzyaaax100))
xyzzyaabt100=xyzzyaabj100*(xyzzyaaae100/xyzzyaaay100+2*xyzzyaaah100-xy&
&zzyaaai100*xyzzyaabj100*(xyzzyaaae100+xyzzyaaah100*xyzzyaaay100))
xyzzyaabv100=xyzzyaabh100*(xyzzyaaae100/xyzzyaaax100+2*xyzzyaaah100-xy&
&zzyaaai100*xyzzyaabh100*(xyzzyaaae100+xyzzyaaah100*xyzzyaaax100))
xyzzyaabx100=xyzzyaabf100*(xyzzyaaae100/xyzzyaaay100+2*xyzzyaaaf100-xy&
&zzyaaag100*xyzzyaabf100*(xyzzyaaae100+xyzzyaaaf100*xyzzyaaay100))
if(fd)then
grad_jas(1)=xyzzyaabo100*xyzzyaaak100
grad_jas(2)=xyzzyaabo100*xyzzyaaaq100
endif
xyzzyaach100(1)=-xyzzyaabr100*xyzzyaaal100-xyzzyaabt100*xyzzyaaan100
xyzzyaach100(2)=-xyzzyaabr100*xyzzyaaar100-xyzzyaabt100*xyzzyaaat100
xyzzyaaci100(1)=-xyzzyaabv100*xyzzyaaal100-xyzzyaabx100*xyzzyaaan100
xyzzyaaci100(2)=-xyzzyaabv100*xyzzyaaar100-xyzzyaabx100*xyzzyaaat100
if(sd)then
xyzzyaabq100=2*xyzzyaaac100*xyzzyaaad100*xyzzyaabc100**2*(xyzzyaaad100&
&*xyzzyaabc100*xyzzyaaaw100-1.d0)
xyzzyaabz100=2*xyzzyaabd100*(xyzzyaaaf100-xyzzyaaag100*xyzzyaabd100*(x&
&yzzyaaae100+xyzzyaaax100*(2*xyzzyaaaf100-xyzzyaaag100*xyzzyaabd100*(x&
&yzzyaaae100+xyzzyaaaf100*xyzzyaaax100))))
xyzzyaacb100=2*xyzzyaabj100*(xyzzyaaah100-xyzzyaaai100*xyzzyaabj100*(x&
&yzzyaaae100+xyzzyaaay100*(2*xyzzyaaah100-xyzzyaaai100*xyzzyaabj100*(x&
&yzzyaaae100+xyzzyaaah100*xyzzyaaay100))))
xyzzyaacd100=2*xyzzyaabh100*(xyzzyaaah100-xyzzyaaai100*xyzzyaabh100*(x&
&yzzyaaae100+xyzzyaaax100*(2*xyzzyaaah100-xyzzyaaai100*xyzzyaabh100*(x&
&yzzyaaae100+xyzzyaaah100*xyzzyaaax100))))
xyzzyaacf100=2*xyzzyaabf100*(xyzzyaaaf100-xyzzyaaag100*xyzzyaabf100*(x&
&yzzyaaae100+xyzzyaaay100*(2*xyzzyaaaf100-xyzzyaaag100*xyzzyaabf100*(x&
&yzzyaaae100+xyzzyaaaf100*xyzzyaaay100))))
lap_jas=xyzzyaabq100+xyzzyaabo100
xyzzyaacj100=xyzzyaabz100+xyzzyaabr100+xyzzyaacb100+xyzzyaabt100
xyzzyaack100=xyzzyaacd100+xyzzyaabv100+xyzzyaacf100+xyzzyaabx100
endif
else
xyzzyaabo100=xyzzyaaac100*xyzzyaabc100*(1.d0/xyzzyaaaw100-xyzzyaaad100&
&*xyzzyaabc100)
xyzzyaabs100=xyzzyaabi100*(xyzzyaaae100/xyzzyaaaz100+2*xyzzyaaah100-xy&
&zzyaaai100*xyzzyaabi100*(xyzzyaaae100+xyzzyaaah100*xyzzyaaaz100))
xyzzyaabu100=xyzzyaabg100*(xyzzyaaae100/xyzzyaaba100+2*xyzzyaaaf100-xy&
&zzyaaag100*xyzzyaabg100*(xyzzyaaae100+xyzzyaaaf100*xyzzyaaba100))
xyzzyaabw100=xyzzyaabe100*(xyzzyaaae100/xyzzyaaaz100+2*xyzzyaaaf100-xy&
&zzyaaag100*xyzzyaabe100*(xyzzyaaae100+xyzzyaaaf100*xyzzyaaaz100))
xyzzyaaby100=xyzzyaabk100*(xyzzyaaae100/xyzzyaaba100+2*xyzzyaaah100-xy&
&zzyaaai100*xyzzyaabk100*(xyzzyaaae100+xyzzyaaah100*xyzzyaaba100))
if(fd)then
grad_jas(1)=-xyzzyaabo100*xyzzyaaak100
grad_jas(2)=-xyzzyaabo100*xyzzyaaaq100
endif
xyzzyaach100(1)=-xyzzyaabs100*xyzzyaaam100-xyzzyaabu100*xyzzyaaao100
xyzzyaach100(2)=-xyzzyaabs100*xyzzyaaas100-xyzzyaabu100*xyzzyaaau100
xyzzyaaci100(1)=-xyzzyaabw100*xyzzyaaam100-xyzzyaaby100*xyzzyaaao100
xyzzyaaci100(2)=-xyzzyaabw100*xyzzyaaas100-xyzzyaaby100*xyzzyaaau100
if(sd)then
xyzzyaabq100=2*xyzzyaaac100*xyzzyaaad100*xyzzyaabc100**2*(xyzzyaaad100&
&*xyzzyaabc100*xyzzyaaaw100-1.d0)
xyzzyaaca100=2*xyzzyaabi100*(xyzzyaaah100-xyzzyaaai100*xyzzyaabi100*(x&
&yzzyaaae100+xyzzyaaaz100*(2*xyzzyaaah100-xyzzyaaai100*xyzzyaabi100*(x&
&yzzyaaae100+xyzzyaaah100*xyzzyaaaz100))))
xyzzyaacc100=2*xyzzyaabg100*(xyzzyaaaf100-xyzzyaaag100*xyzzyaabg100*(x&
&yzzyaaae100+xyzzyaaba100*(2*xyzzyaaaf100-xyzzyaaag100*xyzzyaabg100*(x&
&yzzyaaae100+xyzzyaaaf100*xyzzyaaba100))))
xyzzyaace100=2*xyzzyaabe100*(xyzzyaaaf100-xyzzyaaag100*xyzzyaabe100*(x&
&yzzyaaae100+xyzzyaaaz100*(2*xyzzyaaaf100-xyzzyaaag100*xyzzyaabe100*(x&
&yzzyaaae100+xyzzyaaaf100*xyzzyaaaz100))))
xyzzyaacg100=2*xyzzyaabk100*(xyzzyaaah100-xyzzyaaai100*xyzzyaabk100*(x&
&yzzyaaae100+xyzzyaaba100*(2*xyzzyaaah100-xyzzyaaai100*xyzzyaabk100*(x&
&yzzyaaae100+xyzzyaaah100*xyzzyaaba100))))
lap_jas=xyzzyaabq100+xyzzyaabo100
xyzzyaacj100=xyzzyaaca100+xyzzyaabs100+xyzzyaacc100+xyzzyaabu100
xyzzyaack100=xyzzyaace100+xyzzyaabw100+xyzzyaacg100+xyzzyaaby100
endif
endif
xyzzyaacl100=xyzzyaabl100/xyzzyaacp100
xyzzyaacm100=xyzzyaabm100/xyzzyaacp100
xyzzyaacn100=xyzzyaacl100*xyzzyaach100+xyzzyaacm100*xyzzyaaci100
if(fd)then
grad_jas(1:2)=grad_jas(1:2)+xyzzyaacn100
grad_jas(3)=0.d0
endif
if(sd)lap_jas=lap_jas+xyzzyaacl100*(xyzzyaacj100+xyzzyaach100(1)**2+xy&
&zzyaach100(2)**2)+xyzzyaacm100*(xyzzyaack100+xyzzyaaci100(1)**2+xyzzy&
&aaci100(2)**2)-xyzzyaacn100(1)**2-xyzzyaacn100(2)**2
endif
end subroutine xyzzyaaqf1
subroutine xyzzyaaqg1(val,fd,sd,jspin,eevecs,val_jas,grad_jas,lap_jas)
use slaarnabt, only : ddot,exp_limit
implicit none
integer,intent(in) :: jspin
real(dp),intent(out) :: val_jas,grad_jas(3),lap_jas
real(dp),intent(in) :: eevecs(4,netot,netot)
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa101
real(dp) xyzzyaaab101(2),xyzzyaaac101(2),xyzzyaaad101(2),xyzzyaaae101(&
&2),xyzzyaaaf101(2),xyzzyaaag101(2),xyzzyaaah101(2),xyzzyaaai101(2),xy&
&zzyaaaj101(2),xyzzyaaak101(2),xyzzyaaal101(2),xyzzyaaam101(2),xyzzyaa&
&an101,xyzzyaaao101,xyzzyaaap101,xyzzyaaaq101,xyzzyaaar101,xyzzyaaas10&
&1,xyzzyaaat101,xyzzyaaau101,xyzzyaaav101,xyzzyaaaw101,xyzzyaaax101,xy&
&zzyaaay101,xyzzyaaaz101,xyzzyaaba101,xyzzyaabb101(3),xyzzyaabc101(3),&
&xyzzyaabd101,xyzzyaabe101,xyzzyaabf101,xyzzyaabg101,xyzzyaabh101,xyzz&
&yaabi101,xyzzyaabj101,xyzzyaabk101,xyzzyaabl101
logical xyzzyaabm101,xyzzyaabn101,xyzzyaabo101,xyzzyaabp101,xyzzyaabq1&
&01
xyzzyaabm101=fd.or.sd
if(xyzzyaabm101)then
xyzzyaabn101=jspin==1.or.jspin==3
xyzzyaabo101=.not.xyzzyaabn101
xyzzyaabp101=jspin==1.or.jspin==4
xyzzyaabq101=.not.xyzzyaabp101
else
xyzzyaabn101=.false.
xyzzyaabo101=.false.
xyzzyaabp101=.false.
xyzzyaabq101=.false.
endif
xyzzyaaan101=eevecs(4,3,1)
xyzzyaaao101=eevecs(4,4,2)
xyzzyaaap101=eevecs(4,4,1)
xyzzyaaaq101=eevecs(4,3,2)
do xyzzyaaaa101=1,2
call xyzzyaaqh1(xyzzyaaaa101,xyzzyaaan101,xyzzyaabn101,sd,xyzzyaaab101&
&(xyzzyaaaa101),xyzzyaaaf101(xyzzyaaaa101),xyzzyaaaj101(xyzzyaaaa101))
call xyzzyaaqh1(xyzzyaaaa101,xyzzyaaao101,xyzzyaabo101,sd,xyzzyaaac101&
&(xyzzyaaaa101),xyzzyaaag101(xyzzyaaaa101),xyzzyaaak101(xyzzyaaaa101))
call xyzzyaaqh1(xyzzyaaaa101,xyzzyaaap101,xyzzyaabp101,sd,xyzzyaaad101&
&(xyzzyaaaa101),xyzzyaaah101(xyzzyaaaa101),xyzzyaaal101(xyzzyaaaa101))
call xyzzyaaqh1(xyzzyaaaa101,xyzzyaaaq101,xyzzyaabq101,sd,xyzzyaaae101&
&(xyzzyaaaa101),xyzzyaaai101(xyzzyaaaa101),xyzzyaaam101(xyzzyaaaa101))
enddo
xyzzyaabi101=xyzzyaaab101(1)+xyzzyaaac101(1)+xyzzyaaad101(2)+xyzzyaaae&
&101(2)
xyzzyaaar101=exp_limit(xyzzyaabi101)
xyzzyaabj101=xyzzyaaab101(2)+xyzzyaaac101(2)+xyzzyaaad101(1)+xyzzyaaae&
&101(1)
xyzzyaaas101=exp_limit(xyzzyaabj101)
xyzzyaaat101=xyzzyaaar101+xyzzyaaas101
if(xyzzyaabm101)then
xyzzyaaay101=1.d0/xyzzyaaat101
if(xyzzyaabn101)then
xyzzyaaau101=xyzzyaaar101*xyzzyaaaf101(1)+xyzzyaaas101*xyzzyaaaf101(2)
xyzzyaaaz101=1.d0/xyzzyaaan101
xyzzyaabb101=eevecs(1:3,3,1)*xyzzyaaaz101
if(sd)xyzzyaaaw101=xyzzyaaar101*(xyzzyaaaf101(1)*xyzzyaaaf101(1)+xyzzy&
&aaaj101(1))+xyzzyaaas101*(xyzzyaaaf101(2)*xyzzyaaaf101(2)+xyzzyaaaj10&
&1(2))
else
xyzzyaaau101=xyzzyaaar101*xyzzyaaag101(1)+xyzzyaaas101*xyzzyaaag101(2)
xyzzyaaaz101=1.d0/xyzzyaaao101
xyzzyaabb101=eevecs(1:3,4,2)*xyzzyaaaz101
if(sd)xyzzyaaaw101=xyzzyaaar101*(xyzzyaaag101(1)*xyzzyaaag101(1)+xyzzy&
&aaak101(1))+xyzzyaaas101*(xyzzyaaag101(2)*xyzzyaaag101(2)+xyzzyaaak10&
&1(2))
endif
if(xyzzyaabp101)then
xyzzyaaav101=xyzzyaaar101*xyzzyaaah101(2)+xyzzyaaas101*xyzzyaaah101(1)
xyzzyaaba101=1.d0/xyzzyaaap101
xyzzyaabc101=eevecs(1:3,4,1)*xyzzyaaba101
if(sd)xyzzyaaax101=xyzzyaaar101*(xyzzyaaah101(2)*xyzzyaaah101(2)+xyzzy&
&aaal101(2))+xyzzyaaas101*(xyzzyaaah101(1)*xyzzyaaah101(1)+xyzzyaaal10&
&1(1))
else
xyzzyaaav101=xyzzyaaar101*xyzzyaaai101(2)+xyzzyaaas101*xyzzyaaai101(1)
xyzzyaaba101=1.d0/xyzzyaaaq101
xyzzyaabc101=eevecs(1:3,3,2)*xyzzyaaba101
if(sd)xyzzyaaax101=xyzzyaaar101*(xyzzyaaai101(2)*xyzzyaaai101(2)+xyzzy&
&aaam101(2))+xyzzyaaas101*(xyzzyaaai101(1)*xyzzyaaai101(1)+xyzzyaaam10&
&1(1))
endif
xyzzyaabe101=xyzzyaaau101*xyzzyaaay101
xyzzyaabf101=xyzzyaaav101*xyzzyaaay101
if(jspin>2)then
xyzzyaabb101=-xyzzyaabb101
xyzzyaabc101=-xyzzyaabc101
endif
if(sd)then
select case(jspin)
case(1)
xyzzyaabk101=xyzzyaaar101*xyzzyaaaf101(1)*xyzzyaaah101(2)+xyzzyaaas101&
&*xyzzyaaaf101(2)*xyzzyaaah101(1)
case(2)
xyzzyaabk101=xyzzyaaar101*xyzzyaaag101(1)*xyzzyaaai101(2)+xyzzyaaas101&
&*xyzzyaaag101(2)*xyzzyaaai101(1)
case(3)
xyzzyaabk101=xyzzyaaar101*xyzzyaaaf101(1)*xyzzyaaai101(2)+xyzzyaaas101&
&*xyzzyaaaf101(2)*xyzzyaaai101(1)
case(4)
xyzzyaabk101=xyzzyaaar101*xyzzyaaag101(1)*xyzzyaaah101(2)+xyzzyaaas101&
&*xyzzyaaag101(2)*xyzzyaaah101(1)
end select
xyzzyaabg101=xyzzyaaay101*(xyzzyaaaw101-xyzzyaaay101*xyzzyaaau101*xyzz&
&yaaau101)
xyzzyaabh101=xyzzyaaay101*(xyzzyaaax101-xyzzyaaay101*xyzzyaaav101*xyzz&
&yaaav101)
xyzzyaabl101=xyzzyaaay101*(xyzzyaabk101-xyzzyaaay101*xyzzyaaau101*xyzz&
&yaaav101)
xyzzyaabd101=ddot(dimensionality,xyzzyaabb101(1),1,xyzzyaabc101(1),1)
endif
endif
if(val)val_jas=log(xyzzyaaat101)
if(fd)grad_jas(1:dimensionality)=xyzzyaabe101*xyzzyaabb101(1:dimension&
&ality)+xyzzyaabf101*xyzzyaabc101(1:dimensionality)
if(sd)lap_jas=xyzzyaabg101+2.d0*xyzzyaabl101*xyzzyaabd101+xyzzyaabh101&
&+dble(dimensionality-1)*(xyzzyaabe101*xyzzyaaaz101+xyzzyaabf101*xyzzy&
&aaba101)
end subroutine xyzzyaaqg1
subroutine xyzzyaaqh1(k,r,fsd,sd,f,df,d2f)
implicit none
integer,intent(in) :: k
real(dp),intent(in) :: r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fsd,sd
real(dp) xyzzyaaaa102,xyzzyaaab102,xyzzyaaac102,xyzzyaaad102,xyzzyaaae&
&102
xyzzyaaab102=xyzzyaain1(k)
xyzzyaaaa102=xyzzyaaim1(k)
xyzzyaaac102=1.d0/(xyzzyaaaa102+r)
xyzzyaaad102=xyzzyaaab102*xyzzyaaac102
f=-xyzzyaaad102*r*r
if(fsd)then
xyzzyaaae102=xyzzyaaad102*xyzzyaaac102
df=-xyzzyaaae102*(xyzzyaaaa102+xyzzyaaaa102+r)*r
if(sd)d2f=-2.d0*xyzzyaaae102*xyzzyaaac102*xyzzyaaaa102*xyzzyaaaa102
endif
end subroutine xyzzyaaqh1
subroutine xyzzyaaqi1(val,fd,sd,jspin,rele,val_jas,grad_jas,lap_jas)
use slaarnabt,only : exp_limit
integer,intent(in) :: jspin
real(dp),intent(out) :: val_jas,grad_jas(3),lap_jas
real(dp),intent(in) :: rele(3,netot)
logical,intent(in) :: val,fd,sd
real(dp) xyzzyaaaa103,xyzzyaaab103,xyzzyaaac103,xyzzyaaad103,xyzzyaaae&
&103,xyzzyaaaf103,xyzzyaaag103,xyzzyaaah103,xyzzyaaai103
real(dp) xyzzyaaaj103,xyzzyaaak103,xyzzyaaal103,xyzzyaaam103,xyzzyaaan&
&103,xyzzyaaao103
real(dp) xyzzyaaap103,xyzzyaaaq103,xyzzyaaar103,xyzzyaaas103,xyzzyaaat&
&103,xyzzyaaau103,xyzzyaaav103,xyzzyaaaw103
real(dp) xyzzyaaax103,xyzzyaaay103,xyzzyaaaz103,xyzzyaaba103,xyzzyaabb&
&103,xyzzyaabc103,xyzzyaabd103,xyzzyaabe103
real(dp) xyzzyaabf103,xyzzyaabg103,xyzzyaabh103,xyzzyaabi103,xyzzyaabj&
&103,xyzzyaabk103,xyzzyaabl103,xyzzyaabm103,xyzzyaabn103,xyzzyaabo103
real(dp) xyzzyaabp103,xyzzyaabq103,xyzzyaabr103,xyzzyaabs103,xyzzyaabt&
&103,xyzzyaabu103,xyzzyaabv103,xyzzyaabw103
real(dp) xyzzyaabx103,xyzzyaaby103,xyzzyaabz103,xyzzyaaca103,xyzzyaacb&
&103,xyzzyaacc103,xyzzyaacd103,xyzzyaace103,xyzzyaacf103,xyzzyaacg103,&
&xyzzyaach103,xyzzyaaci103,xyzzyaacj103,xyzzyaack103,xyzzyaacl103,xyzz&
&yaacm103,xyzzyaacn103,xyzzyaaco103,xyzzyaacp103(2),xyzzyaacq103(2),xy&
&zzyaacr103,xyzzyaacs103,xyzzyaact103,xyzzyaacu103,xyzzyaacv103,xyzzya&
&acw103,xyzzyaacx103,xyzzyaacy103,xyzzyaacz103,xyzzyaada103,xyzzyaadb1&
&03,xyzzyaadc103,xyzzyaadd103,xyzzyaade103
real(dp),save :: xyzzyaadf103,xyzzyaadg103,xyzzyaadh103,xyzzyaadi103,x&
&yzzyaadj103
logical,save :: xyzzyaadk103=.true.
xyzzyaaaa103=xyzzyaaih1(1)
xyzzyaaab103=xyzzyaaih1(2)
xyzzyaaac103=xyzzyaaih1(3)
xyzzyaaad103=xyzzyaaih1(4)
xyzzyaaae103=xyzzyaaih1(5)
xyzzyaaaf103=xyzzyaaih1(6)
xyzzyaaag103=xyzzyaaih1(7)
xyzzyaaah103=xyzzyaaih1(8)
xyzzyaaai103=xyzzyaaih1(9)
if(xyzzyaadk103)then
xyzzyaadf103=mu_biex3/me_biex3
xyzzyaadg103=mu_biex3/mh_biex3
xyzzyaadh103=-xyzzyaaaa103*xx_sep/(1.d0+xyzzyaaab103*xx_sep)
xyzzyaadi103=-xyzzyaaac103*xx_sep/(1.d0+xyzzyaaad103*xx_sep)
xyzzyaadj103=-(xyzzyaaae103*xx_sep+xyzzyaaah103*(xx_sep**2))/(1.d0+xyz&
&zyaaai103*xx_sep)
xyzzyaadj103=xyzzyaadj103-(xyzzyaaae103*xx_sep+xyzzyaaaf103*(xx_sep**2&
&))/(1.d0+xyzzyaaag103*xx_sep)
xyzzyaadk103=.false.
endif
xyzzyaaaj103=rele(1,1)
xyzzyaaam103=rele(2,1)
xyzzyaaar103=sqrt(xyzzyaaaj103*xyzzyaaaj103+xyzzyaaam103*xyzzyaaam103)
xyzzyaaak103=rele(1,2)
xyzzyaaan103=rele(2,2)
xyzzyaaau103=sqrt(xyzzyaaak103*xyzzyaaak103+xyzzyaaan103*xyzzyaaan103)
xyzzyaaal103=xyzzyaaaj103-xyzzyaaak103
xyzzyaaao103=xyzzyaaam103-xyzzyaaan103
if(xyzzyaaar103>0.d0)then
xyzzyaaav103=1.d0/xyzzyaaar103
else
xyzzyaaav103=0.d0
endif
if(xyzzyaaau103>0.d0)then
xyzzyaaaw103=1.d0/xyzzyaaau103
else
xyzzyaaaw103=0.d0
endif
xyzzyaaax103=xx_sep+xyzzyaadf103*xyzzyaaal103
xyzzyaaay103=xyzzyaadf103*xyzzyaaao103
xyzzyaaaz103=xx_sep+xyzzyaadg103*(-xyzzyaaal103)
xyzzyaaba103=xyzzyaadg103*(-xyzzyaaao103)
xyzzyaabb103=xx_sep+xyzzyaadf103*xyzzyaaaj103+xyzzyaadg103*xyzzyaaak10&
&3
xyzzyaabc103=xyzzyaadf103*xyzzyaaam103+xyzzyaadg103*xyzzyaaan103
xyzzyaabd103=xx_sep-xyzzyaadg103*xyzzyaaaj103-xyzzyaadf103*xyzzyaaak10&
&3
xyzzyaabe103=-xyzzyaadg103*xyzzyaaam103-xyzzyaadf103*xyzzyaaan103
xyzzyaaap103=sqrt(xyzzyaaax103*xyzzyaaax103+xyzzyaaay103*xyzzyaaay103)
xyzzyaaaq103=sqrt(xyzzyaaaz103*xyzzyaaaz103+xyzzyaaba103*xyzzyaaba103)
xyzzyaaas103=sqrt(xyzzyaabb103*xyzzyaabb103+xyzzyaabc103*xyzzyaabc103)
xyzzyaaat103=sqrt(xyzzyaabd103*xyzzyaabd103+xyzzyaabe103*xyzzyaabe103)
xyzzyaabf103=1.d0/(1.d0+xyzzyaaab103*xyzzyaaap103)
xyzzyaabg103=1.d0/(1.d0+xyzzyaaad103*xyzzyaaaq103)
xyzzyaabh103=1.d0/(1.d0+xyzzyaaag103*xyzzyaaar103)
xyzzyaabp103=(xyzzyaaae103+xyzzyaaaf103*xyzzyaaar103)*xyzzyaaar103
xyzzyaabl103=1.d0/(1.d0+xyzzyaaai103*xyzzyaaar103)
xyzzyaabt103=(xyzzyaaae103+xyzzyaaah103*xyzzyaaar103)*xyzzyaaar103
xyzzyaabi103=1.d0/(1.d0+xyzzyaaag103*xyzzyaaas103)
xyzzyaabq103=(xyzzyaaae103+xyzzyaaaf103*xyzzyaaas103)*xyzzyaaas103
xyzzyaabm103=1.d0/(1.d0+xyzzyaaai103*xyzzyaaas103)
xyzzyaabu103=(xyzzyaaae103+xyzzyaaah103*xyzzyaaas103)*xyzzyaaas103
xyzzyaabj103=1.d0/(1.d0+xyzzyaaag103*xyzzyaaat103)
xyzzyaabr103=(xyzzyaaae103+xyzzyaaaf103*xyzzyaaat103)*xyzzyaaat103
xyzzyaabn103=1.d0/(1.d0+xyzzyaaai103*xyzzyaaat103)
xyzzyaabv103=(xyzzyaaae103+xyzzyaaah103*xyzzyaaat103)*xyzzyaaat103
xyzzyaabk103=1.d0/(1.d0+xyzzyaaag103*xyzzyaaau103)
xyzzyaabs103=(xyzzyaaae103+xyzzyaaaf103*xyzzyaaau103)*xyzzyaaau103
xyzzyaabo103=1.d0/(1.d0+xyzzyaaai103*xyzzyaaau103)
xyzzyaabw103=(xyzzyaaae103+xyzzyaaah103*xyzzyaaau103)*xyzzyaaau103
xyzzyaacb103=xyzzyaabp103*xyzzyaabh103+xyzzyaabu103*xyzzyaabm103+xyzzy&
&aabv103*xyzzyaabn103+xyzzyaabs103*xyzzyaabk103+xyzzyaadj103
xyzzyaabx103=exp_limit(xyzzyaacb103)
xyzzyaacb103=xyzzyaabt103*xyzzyaabl103+xyzzyaabq103*xyzzyaabi103+xyzzy&
&aabr103*xyzzyaabj103+xyzzyaabw103*xyzzyaabo103+xyzzyaadj103
xyzzyaaby103=exp_limit(xyzzyaacb103)
xyzzyaacc103=xyzzyaabx103+xyzzyaaby103
if(val)val_jas=xyzzyaaaa103*xyzzyaaap103*xyzzyaabf103+xyzzyaadh103+xyz&
&zyaaac103*xyzzyaaaq103*xyzzyaabg103+xyzzyaadi103+log(xyzzyaacc103)
if(fd.or.sd)then
if(jspin==1)then
xyzzyaabz103=xyzzyaaaa103*xyzzyaabf103*xyzzyaadf103*(1.d0/xyzzyaaap103&
&-xyzzyaaab103*xyzzyaabf103)
xyzzyaaca103=xyzzyaaac103*xyzzyaabg103*xyzzyaadg103*(xyzzyaabg103*xyzz&
&yaaad103-1.d0/xyzzyaaaq103)
xyzzyaacd103=(-xyzzyaaai103*xyzzyaadf103*xyzzyaabu103/xyzzyaaas103)*xy&
&zzyaabm103**2
xyzzyaace103=(2.d0*xyzzyaaah103*xyzzyaadf103+xyzzyaaae103*xyzzyaadf103&
&/xyzzyaaas103)*xyzzyaabm103
xyzzyaacf103=(xyzzyaaai103*xyzzyaadg103*xyzzyaabv103/xyzzyaaat103)*xyz&
&zyaabn103**2
xyzzyaacg103=(-2.d0*xyzzyaaah103*xyzzyaadg103-xyzzyaaae103*xyzzyaadg10&
&3/xyzzyaaat103)*xyzzyaabn103
xyzzyaach103=(-xyzzyaaag103*xyzzyaabp103*xyzzyaaav103)*xyzzyaabh103**2
xyzzyaaci103=(2.d0*xyzzyaaaf103+xyzzyaaae103*xyzzyaaav103)*xyzzyaabh10&
&3
xyzzyaacj103=(-xyzzyaaag103*xyzzyaadf103*xyzzyaabq103/xyzzyaaas103)*xy&
&zzyaabi103**2
xyzzyaack103=(2.d0*xyzzyaaaf103*xyzzyaadf103+xyzzyaaae103*xyzzyaadf103&
&/xyzzyaaas103)*xyzzyaabi103
xyzzyaacl103=(xyzzyaaag103*xyzzyaadg103*xyzzyaabr103/xyzzyaaat103)*xyz&
&zyaabj103**2
xyzzyaacm103=(-2.d0*xyzzyaaaf103*xyzzyaadg103-xyzzyaaae103*xyzzyaadg10&
&3/xyzzyaaat103)*xyzzyaabj103
xyzzyaacn103=(-xyzzyaaai103*xyzzyaabt103*xyzzyaaav103)*xyzzyaabl103**2
xyzzyaaco103=(2.d0*xyzzyaaah103+xyzzyaaae103*xyzzyaaav103)*xyzzyaabl10&
&3
xyzzyaacp103(1)=(xyzzyaacd103+xyzzyaace103)*xyzzyaabb103+(xyzzyaacf103&
&+xyzzyaacg103)*xyzzyaabd103+(xyzzyaach103+xyzzyaaci103)*xyzzyaaaj103
xyzzyaacq103(1)=(xyzzyaacj103+xyzzyaack103)*xyzzyaabb103+(xyzzyaacl103&
&+xyzzyaacm103)*xyzzyaabd103+(xyzzyaacn103+xyzzyaaco103)*xyzzyaaaj103
xyzzyaacp103(2)=(xyzzyaacd103+xyzzyaace103)*xyzzyaabc103+(xyzzyaacf103&
&+xyzzyaacg103)*xyzzyaabe103+(xyzzyaach103+xyzzyaaci103)*xyzzyaaam103
xyzzyaacq103(2)=(xyzzyaacj103+xyzzyaack103)*xyzzyaabc103+(xyzzyaacl103&
&+xyzzyaacm103)*xyzzyaabe103+(xyzzyaacn103+xyzzyaaco103)*xyzzyaaam103
xyzzyaacp103(:)=xyzzyaacp103(:)*xyzzyaabx103
xyzzyaacq103(:)=xyzzyaacq103(:)*xyzzyaaby103
if(fd)then
grad_jas(1)=xyzzyaabz103*xyzzyaaax103+xyzzyaaca103*xyzzyaaaz103
grad_jas(1)=(xyzzyaacp103(1)+xyzzyaacq103(1))/(xyzzyaacc103)+grad_jas(&
&1)
grad_jas(2)=xyzzyaabz103*xyzzyaaay103+xyzzyaaca103*xyzzyaaba103
grad_jas(2)=(xyzzyaacp103(2)+xyzzyaacq103(2))/(xyzzyaacc103)+grad_jas(&
&2)
endif
if(sd)then
call xyzzyaaqj1(xyzzyaabz103,xyzzyaadf103,xyzzyaaaa103,xyzzyaaab103,xy&
&zzyaabf103,xyzzyaaap103,xyzzyaaax103,xyzzyaacv103)
call xyzzyaaqj1(xyzzyaabz103,xyzzyaadf103,xyzzyaaaa103,xyzzyaaab103,xy&
&zzyaabf103,xyzzyaaap103,xyzzyaaay103,xyzzyaacw103)
call xyzzyaaqj1(xyzzyaaca103,-xyzzyaadg103,xyzzyaaac103,xyzzyaaad103,x&
&yzzyaabg103,xyzzyaaaq103,xyzzyaaaz103,xyzzyaacx103)
call xyzzyaaqj1(xyzzyaaca103,-xyzzyaadg103,xyzzyaaac103,xyzzyaaad103,x&
&yzzyaabg103,xyzzyaaaq103,xyzzyaaba103,xyzzyaacy103)
call xyzzyaaqk1(xyzzyaaar103,xyzzyaaaj103,xyzzyaaas103,xyzzyaabb103,xy&
&zzyaaat103,xyzzyaabd103,xyzzyaabp103,xyzzyaabu103,xyzzyaabv103,xyzzya&
&abh103,xyzzyaabm103,xyzzyaabn103,xyzzyaadf103,xyzzyaadg103,xyzzyaaae1&
&03,xyzzyaaaf103,xyzzyaaag103,xyzzyaaah103,xyzzyaaai103,xyzzyaacz103)
xyzzyaacz103=xyzzyaabx103*xyzzyaacz103
call xyzzyaaqk1(xyzzyaaar103,xyzzyaaaj103,xyzzyaaas103,xyzzyaabb103,xy&
&zzyaaat103,xyzzyaabd103,xyzzyaabt103,xyzzyaabq103,xyzzyaabr103,xyzzya&
&abl103,xyzzyaabi103,xyzzyaabj103,xyzzyaadf103,xyzzyaadg103,xyzzyaaae1&
&03,xyzzyaaah103,xyzzyaaai103,xyzzyaaaf103,xyzzyaaag103,xyzzyaada103)
xyzzyaada103=xyzzyaaby103*xyzzyaada103
xyzzyaadd103=xyzzyaacz103+xyzzyaada103
call xyzzyaaqk1(xyzzyaaar103,xyzzyaaam103,xyzzyaaas103,xyzzyaabc103,xy&
&zzyaaat103,xyzzyaabe103,xyzzyaabp103,xyzzyaabu103,xyzzyaabv103,xyzzya&
&abh103,xyzzyaabm103,xyzzyaabn103,xyzzyaadf103,xyzzyaadg103,xyzzyaaae1&
&03,xyzzyaaaf103,xyzzyaaag103,xyzzyaaah103,xyzzyaaai103,xyzzyaadb103)
xyzzyaadb103=xyzzyaabx103*xyzzyaadb103
call xyzzyaaqk1(xyzzyaaar103,xyzzyaaam103,xyzzyaaas103,xyzzyaabc103,xy&
&zzyaaat103,xyzzyaabe103,xyzzyaabt103,xyzzyaabq103,xyzzyaabr103,xyzzya&
&abl103,xyzzyaabi103,xyzzyaabj103,xyzzyaadf103,xyzzyaadg103,xyzzyaaae1&
&03,xyzzyaaah103,xyzzyaaai103,xyzzyaaaf103,xyzzyaaag103,xyzzyaadc103)
xyzzyaadc103=xyzzyaaby103*xyzzyaadc103
xyzzyaade103=xyzzyaadb103+xyzzyaadc103
lap_jas=(xyzzyaadd103+xyzzyaade103+(-(xyzzyaacp103(1)+xyzzyaacq103(1))&
&**2-(xyzzyaacp103(2)+xyzzyaacq103(2))**2)/xyzzyaacc103)/xyzzyaacc103
lap_jas=lap_jas+xyzzyaacv103+xyzzyaacw103+xyzzyaacx103+xyzzyaacy103
endif
else
xyzzyaabz103=-xyzzyaaaa103*xyzzyaabf103*xyzzyaadf103*(1.d0/xyzzyaaap10&
&3-xyzzyaaab103*xyzzyaabf103)
xyzzyaaca103=-xyzzyaaac103*xyzzyaabg103*xyzzyaadg103*(xyzzyaabg103*xyz&
&zyaaad103-1.d0/xyzzyaaaq103)
xyzzyaacd103=(-xyzzyaaai103*xyzzyaadg103*xyzzyaabu103/xyzzyaaas103)*xy&
&zzyaabm103**2
xyzzyaace103=(2.d0*xyzzyaaah103*xyzzyaadg103+xyzzyaaae103*xyzzyaadg103&
&/xyzzyaaas103)*xyzzyaabm103
xyzzyaacf103=(xyzzyaaai103*xyzzyaadf103*xyzzyaabv103/xyzzyaaat103)*xyz&
&zyaabn103**2
xyzzyaacg103=(-2.d0*xyzzyaaah103*xyzzyaadf103-xyzzyaaae103*xyzzyaadf10&
&3/xyzzyaaat103)*xyzzyaabn103
xyzzyaact103=(-xyzzyaaag103*xyzzyaabs103*xyzzyaaaw103)*xyzzyaabk103**2
xyzzyaacu103=(2.d0*xyzzyaaaf103+xyzzyaaae103*xyzzyaaaw103)*xyzzyaabk10&
&3
xyzzyaacj103=(-xyzzyaaag103*xyzzyaadg103*xyzzyaabq103/xyzzyaaas103)*xy&
&zzyaabi103**2
xyzzyaack103=(2.d0*xyzzyaaaf103*xyzzyaadg103+xyzzyaaae103*xyzzyaadg103&
&/xyzzyaaas103)*xyzzyaabi103
xyzzyaacl103=(xyzzyaaag103*xyzzyaadf103*xyzzyaabr103/xyzzyaaat103)*xyz&
&zyaabj103**2
xyzzyaacm103=(-2.d0*xyzzyaaaf103*xyzzyaadf103-xyzzyaaae103*xyzzyaadf10&
&3/xyzzyaaat103)*xyzzyaabj103
xyzzyaacr103=(-xyzzyaaai103*xyzzyaabw103*xyzzyaaaw103)*xyzzyaabo103**2
xyzzyaacs103=(2.d0*xyzzyaaah103+xyzzyaaae103*xyzzyaaaw103)*xyzzyaabo10&
&3
xyzzyaacp103(1)=(xyzzyaacd103+xyzzyaace103)*xyzzyaabb103+(xyzzyaacf103&
&+xyzzyaacg103)*xyzzyaabd103+(xyzzyaact103+xyzzyaacu103)*xyzzyaaak103
xyzzyaacq103(1)=(xyzzyaacj103+xyzzyaack103)*xyzzyaabb103+(xyzzyaacl103&
&+xyzzyaacm103)*xyzzyaabd103+(xyzzyaacr103+xyzzyaacs103)*xyzzyaaak103
xyzzyaacp103(2)=(xyzzyaacd103+xyzzyaace103)*xyzzyaabc103+(xyzzyaacf103&
&+xyzzyaacg103)*xyzzyaabe103+(xyzzyaact103+xyzzyaacu103)*xyzzyaaan103
xyzzyaacq103(2)=(xyzzyaacj103+xyzzyaack103)*xyzzyaabc103+(xyzzyaacl103&
&+xyzzyaacm103)*xyzzyaabe103+(xyzzyaacr103+xyzzyaacs103)*xyzzyaaan103
xyzzyaacp103(:)=xyzzyaacp103(:)*xyzzyaabx103
xyzzyaacq103(:)=xyzzyaacq103(:)*xyzzyaaby103
if(fd)then
grad_jas(1)=xyzzyaabz103*xyzzyaaax103+xyzzyaaca103*xyzzyaaaz103
grad_jas(1)=(xyzzyaacp103(1)+xyzzyaacq103(1))/(xyzzyaacc103)+grad_jas(&
&1)
grad_jas(2)=xyzzyaabz103*xyzzyaaay103+xyzzyaaca103*xyzzyaaba103
grad_jas(2)=(xyzzyaacp103(2)+xyzzyaacq103(2))/(xyzzyaacc103)+grad_jas(&
&2)
endif
if(sd)then
call xyzzyaaqj1(xyzzyaabz103,-xyzzyaadf103,xyzzyaaaa103,xyzzyaaab103,x&
&yzzyaabf103,xyzzyaaap103,xyzzyaaax103,xyzzyaacv103)
call xyzzyaaqj1(xyzzyaabz103,-xyzzyaadf103,xyzzyaaaa103,xyzzyaaab103,x&
&yzzyaabf103,xyzzyaaap103,xyzzyaaay103,xyzzyaacw103)
call xyzzyaaqj1(xyzzyaaca103,xyzzyaadg103,xyzzyaaac103,xyzzyaaad103,xy&
&zzyaabg103,xyzzyaaaq103,xyzzyaaaz103,xyzzyaacx103)
call xyzzyaaqj1(xyzzyaaca103,xyzzyaadg103,xyzzyaaac103,xyzzyaaad103,xy&
&zzyaabg103,xyzzyaaaq103,xyzzyaaba103,xyzzyaacy103)
call xyzzyaaqk1(xyzzyaaau103,xyzzyaaak103,xyzzyaaas103,xyzzyaabb103,xy&
&zzyaaat103,xyzzyaabd103,xyzzyaabs103,xyzzyaabu103,xyzzyaabv103,xyzzya&
&abk103,xyzzyaabm103,xyzzyaabn103,xyzzyaadg103,xyzzyaadf103,xyzzyaaae1&
&03,xyzzyaaaf103,xyzzyaaag103,xyzzyaaah103,xyzzyaaai103,xyzzyaacz103)
xyzzyaacz103=xyzzyaabx103*xyzzyaacz103
call xyzzyaaqk1(xyzzyaaau103,xyzzyaaak103,xyzzyaaas103,xyzzyaabb103,xy&
&zzyaaat103,xyzzyaabd103,xyzzyaabw103,xyzzyaabq103,xyzzyaabr103,xyzzya&
&abo103,xyzzyaabi103,xyzzyaabj103,xyzzyaadg103,xyzzyaadf103,xyzzyaaae1&
&03,xyzzyaaah103,xyzzyaaai103,xyzzyaaaf103,xyzzyaaag103,xyzzyaada103)
xyzzyaada103=xyzzyaaby103*xyzzyaada103
xyzzyaadd103=xyzzyaacz103+xyzzyaada103
call xyzzyaaqk1(xyzzyaaau103,xyzzyaaan103,xyzzyaaas103,xyzzyaabc103,xy&
&zzyaaat103,xyzzyaabe103,xyzzyaabs103,xyzzyaabu103,xyzzyaabv103,xyzzya&
&abk103,xyzzyaabm103,xyzzyaabn103,xyzzyaadg103,xyzzyaadf103,xyzzyaaae1&
&03,xyzzyaaaf103,xyzzyaaag103,xyzzyaaah103,xyzzyaaai103,xyzzyaadb103)
xyzzyaadb103=xyzzyaabx103*xyzzyaadb103
call xyzzyaaqk1(xyzzyaaau103,xyzzyaaan103,xyzzyaaas103,xyzzyaabc103,xy&
&zzyaaat103,xyzzyaabe103,xyzzyaabw103,xyzzyaabq103,xyzzyaabr103,xyzzya&
&abo103,xyzzyaabi103,xyzzyaabj103,xyzzyaadg103,xyzzyaadf103,xyzzyaaae1&
&03,xyzzyaaah103,xyzzyaaai103,xyzzyaaaf103,xyzzyaaag103,xyzzyaadc103)
xyzzyaadc103=xyzzyaaby103*xyzzyaadc103
xyzzyaade103=xyzzyaadb103+xyzzyaadc103
lap_jas=(xyzzyaadd103+xyzzyaade103+(-(xyzzyaacp103(1)+xyzzyaacq103(1))&
&**2-(xyzzyaacp103(2)+xyzzyaacq103(2))**2)/xyzzyaacc103)/xyzzyaacc103
lap_jas=lap_jas+xyzzyaacv103+xyzzyaacw103+xyzzyaacx103+xyzzyaacy103
endif
endif
endif
end subroutine xyzzyaaqi1
subroutine xyzzyaaqj1(dj,m,c1,c2,xyzzyaaab1,r,r_d,d2j)
implicit none
real(dp),intent(in) :: dj,m,c1,c2,xyzzyaaab1,r,r_d
real(dp),intent(out) :: d2j
real(dp) xyzzyaaaa104,xyzzyaaab104,xyzzyaaac104,xyzzyaaad104,xyzzyaaae&
&104
if(r/=0.d0)then
xyzzyaaae104=1.d0/r
else
xyzzyaaae104=0.d0
endif
xyzzyaaaa104=dj*m
xyzzyaaab104=2*c1*c2*c2*(m**2)*(xyzzyaaab1**3)*r_d*r_d*xyzzyaaae104
xyzzyaaac104=-c1*c2*(m**2)*r_d*r_d*(xyzzyaaab1**2)*xyzzyaaae104**2
xyzzyaaad104=-c1*(m**2)*r_d*r_d*xyzzyaaab1*xyzzyaaae104**3
d2j=xyzzyaaaa104+xyzzyaaab104+xyzzyaaac104+xyzzyaaad104
end subroutine xyzzyaaqj1
subroutine xyzzyaaqk1(r,r_d,u,u_d,v,v_d,beta_r,beta_u,beta_v,alpha_r,a&
&lpha_u,alpha_v,m_u,m_v,c5,c6,c7,c8,c9,term_out)
implicit none
real(dp),intent(in) :: alpha_r,alpha_u,alpha_v,beta_r,beta_u,beta_v,r,&
&r_d,u,u_d,v,v_d,m_u,m_v,c5,c6,c7,c8,c9
real(dp),intent(out) :: term_out
real(dp) xyzzyaaaa105,xyzzyaaab105,xyzzyaaac105(6),xyzzyaaad105(16),xy&
&zzyaaae105,xyzzyaaaf105,xyzzyaaag105
if(r/=0.d0)then
xyzzyaaae105=1.d0/r
else
xyzzyaaae105=0.d0
endif
if(u/=0.d0)then
xyzzyaaaf105=1.d0/u
else
xyzzyaaaf105=0.d0
endif
if(v/=0.d0)then
xyzzyaaag105=1.d0/v
else
xyzzyaaag105=0.d0
endif
xyzzyaaaa105=2*c8*m_u*u_d+c5*m_u*u_d*xyzzyaaaf105
xyzzyaaab105=-2*c8*m_v*v_d-c5*m_v*v_d*xyzzyaaag105
xyzzyaaac105(1)=xyzzyaaaa105*alpha_u
xyzzyaaac105(2)=xyzzyaaab105*alpha_v
xyzzyaaac105(3)=-beta_u*c9*m_u*u_d*(alpha_u**2)*xyzzyaaaf105
xyzzyaaac105(4)=beta_v*c9*m_v*v_d*(alpha_v**2)*xyzzyaaag105
xyzzyaaac105(5)=-beta_r*c7*r_d*(alpha_r**2)*xyzzyaaae105
xyzzyaaac105(6)=alpha_r*(2*c6*r_d+c5*r_d*xyzzyaaae105)
xyzzyaaad105(1)=(sum(xyzzyaaac105(:)))**2
xyzzyaaad105(2)=-beta_r*c7*(alpha_r**2)*xyzzyaaae105
xyzzyaaad105(3)=-beta_u*c9*((m_u*alpha_u)**2)*xyzzyaaaf105
xyzzyaaad105(4)=-2*c9*xyzzyaaaa105*m_u*u_d*(alpha_u**2)*xyzzyaaaf105
xyzzyaaad105(5)=2*beta_u*c9*c9*(m_u**2)*u_d*u_d*(alpha_u**3)*xyzzyaaaf&
&105**2
xyzzyaaad105(6)=beta_u*c9*(m_u**2)*u_d*u_d*(alpha_u**2)*xyzzyaaaf105**&
&3
xyzzyaaad105(7)=alpha_u*(m_u**2)*(2*c8+c5*xyzzyaaaf105-c5*u_d*u_d*xyzz&
&yaaaf105**3)
xyzzyaaad105(8)=-beta_v*c9*(m_v**2)*(alpha_v**2)*xyzzyaaag105
xyzzyaaad105(9)=2*c9*xyzzyaaab105*m_v*v_d*(alpha_v**2)*xyzzyaaag105
xyzzyaaad105(10)=2*beta_v*c9*c9*(m_v**2)*v_d*v_d*(alpha_v**3)*xyzzyaaa&
&g105**2
xyzzyaaad105(11)=beta_v*c9*(m_v**2)*v_d*v_d*(alpha_v**2)*xyzzyaaag105*&
&*3
xyzzyaaad105(12)=alpha_v*(m_v**2)*(2*c8+c5*xyzzyaaag105-c5*v_d*v_d*xyz&
&zyaaag105**3)
xyzzyaaad105(13)=2*beta_r*c7*c7*r_d*r_d*(alpha_r**3)*xyzzyaaae105**2
xyzzyaaad105(14)=beta_r*c7*r_d*r_d*(alpha_r**2)*xyzzyaaae105**3
xyzzyaaad105(15)=-2*c7*r_d*r_d*(alpha_r**2)*(2*c6+c5*xyzzyaaae105)*xyz&
&zyaaae105
xyzzyaaad105(16)=alpha_r*(2*c6+c5*xyzzyaaae105-c5*r_d*r_d*xyzzyaaae105&
&**3)
term_out=sum(xyzzyaaad105(:))
end subroutine xyzzyaaqk1
subroutine xyzzyaaql1(val,fd,sd,jspin,rele,val_jas,grad_jas,lap_jas)
use slaarnabt,only : exp_limit
integer,intent(in) :: jspin
real(dp),intent(out) :: val_jas,grad_jas(3),lap_jas
real(dp),intent(in) :: rele(3,netot)
logical,intent(in) :: val,fd,sd
real(dp) xyzzyaaaa106,xyzzyaaab106,xyzzyaaac106,xyzzyaaad106,xyzzyaaae&
&106,xyzzyaaaf106,xyzzyaaag106,xyzzyaaah106,xyzzyaaai106
real(dp) xyzzyaaaj106,xyzzyaaak106,xyzzyaaal106,xyzzyaaam106,xyzzyaaan&
&106,xyzzyaaao106
real(dp) xyzzyaaap106,xyzzyaaaq106,xyzzyaaar106,xyzzyaaas106,xyzzyaaat&
&106,xyzzyaaau106,xyzzyaaav106,xyzzyaaaw106
real(dp) xyzzyaaax106,xyzzyaaay106,xyzzyaaaz106,xyzzyaaba106,xyzzyaabb&
&106,xyzzyaabc106,xyzzyaabd106,xyzzyaabe106
real(dp) xyzzyaabf106,xyzzyaabg106,xyzzyaabh106,xyzzyaabi106,xyzzyaabj&
&106,xyzzyaabk106,xyzzyaabl106,xyzzyaabm106,xyzzyaabn106,xyzzyaabo106
real(dp) xyzzyaabp106,xyzzyaabq106,xyzzyaabr106,xyzzyaabs106,xyzzyaabt&
&106,xyzzyaabu106,xyzzyaabv106,xyzzyaabw106
real(dp) xyzzyaabx106,xyzzyaaby106,xyzzyaabz106,xyzzyaaca106,xyzzyaacb&
&106,xyzzyaacc106,xyzzyaacd106,xyzzyaace106,xyzzyaacf106,xyzzyaacg106,&
&xyzzyaach106,xyzzyaaci106,xyzzyaacj106,xyzzyaack106,xyzzyaacl106,xyzz&
&yaacm106,xyzzyaacn106,xyzzyaaco106,xyzzyaacp106,xyzzyaacq106,xyzzyaac&
&r106,xyzzyaacs106,xyzzyaact106,xyzzyaacu106
real(dp),save :: xyzzyaacv106,xyzzyaacw106,xyzzyaacx106
logical,save :: xyzzyaacy106=.true.
xyzzyaaaa106=xyzzyaaih1(1)
xyzzyaaab106=xyzzyaaih1(2)
xyzzyaaac106=xyzzyaaih1(3)
xyzzyaaad106=xyzzyaaih1(4)
xyzzyaaae106=xyzzyaaih1(5)
xyzzyaaaf106=xyzzyaaih1(6)
xyzzyaaag106=xyzzyaaih1(7)
xyzzyaaah106=xyzzyaaih1(8)
xyzzyaaai106=xyzzyaaih1(9)
if(xyzzyaacy106)then
xyzzyaacv106=-xyzzyaaaa106*xx_sep/(1.d0+xyzzyaaab106*xx_sep)
xyzzyaacw106=-xyzzyaaac106*xx_sep/(1.d0+xyzzyaaad106*xx_sep)
xyzzyaacx106=-(xyzzyaaae106*xx_sep+xyzzyaaah106*(xx_sep**2))/(1.d0+xyz&
&zyaaai106*xx_sep)
xyzzyaacx106=xyzzyaacx106-(xyzzyaaae106*xx_sep+xyzzyaaaf106*(xx_sep**2&
&))/(1.d0+xyzzyaaag106*xx_sep)
xyzzyaacy106=.false.
endif
xyzzyaaaj106=rele(1,1)
xyzzyaaam106=rele(2,1)
xyzzyaaar106=sqrt(xyzzyaaaj106*xyzzyaaaj106+xyzzyaaam106*xyzzyaaam106)
xyzzyaaak106=rele(1,2)
xyzzyaaan106=rele(2,2)
xyzzyaaau106=sqrt(xyzzyaaak106*xyzzyaaak106+xyzzyaaan106*xyzzyaaan106)
xyzzyaaal106=xyzzyaaaj106-xyzzyaaak106
xyzzyaaao106=xyzzyaaam106-xyzzyaaan106
if(xyzzyaaar106>0.d0)then
xyzzyaaav106=1.d0/xyzzyaaar106
else
xyzzyaaav106=0.d0
endif
if(xyzzyaaau106>0.d0)then
xyzzyaaaw106=1.d0/xyzzyaaau106
else
xyzzyaaaw106=0.d0
endif
xyzzyaaax106=xx_sep-xyzzyaaal106
xyzzyaaay106=-xyzzyaaao106
xyzzyaaaz106=xx_sep
xyzzyaaba106=0.d0
xyzzyaabb106=xx_sep-xyzzyaaaj106
xyzzyaabc106=xyzzyaaam106
xyzzyaabd106=xx_sep+xyzzyaaak106
xyzzyaabe106=xyzzyaaan106
xyzzyaaap106=sqrt(xyzzyaaax106*xyzzyaaax106+xyzzyaaay106*xyzzyaaay106)
xyzzyaaaq106=sqrt(xyzzyaaaz106*xyzzyaaaz106+xyzzyaaba106*xyzzyaaba106)
xyzzyaaas106=sqrt(xyzzyaabb106*xyzzyaabb106+xyzzyaabc106*xyzzyaabc106)
xyzzyaaat106=sqrt(xyzzyaabd106*xyzzyaabd106+xyzzyaabe106*xyzzyaabe106)
xyzzyaabf106=1.d0/(1.d0+xyzzyaaab106*xyzzyaaap106)
xyzzyaabg106=1.d0/(1.d0+xyzzyaaad106*xyzzyaaaq106)
xyzzyaabh106=1.d0/(1.d0+xyzzyaaag106*xyzzyaaar106)
xyzzyaabp106=(xyzzyaaae106+xyzzyaaaf106*xyzzyaaar106)*xyzzyaaar106
xyzzyaabl106=1.d0/(1.d0+xyzzyaaai106*xyzzyaaar106)
xyzzyaabt106=(xyzzyaaae106+xyzzyaaah106*xyzzyaaar106)*xyzzyaaar106
xyzzyaabi106=1.d0/(1.d0+xyzzyaaag106*xyzzyaaas106)
xyzzyaabq106=(xyzzyaaae106+xyzzyaaaf106*xyzzyaaas106)*xyzzyaaas106
xyzzyaabm106=1.d0/(1.d0+xyzzyaaai106*xyzzyaaas106)
xyzzyaabu106=(xyzzyaaae106+xyzzyaaah106*xyzzyaaas106)*xyzzyaaas106
xyzzyaabj106=1.d0/(1.d0+xyzzyaaag106*xyzzyaaat106)
xyzzyaabr106=(xyzzyaaae106+xyzzyaaaf106*xyzzyaaat106)*xyzzyaaat106
xyzzyaabn106=1.d0/(1.d0+xyzzyaaai106*xyzzyaaat106)
xyzzyaabv106=(xyzzyaaae106+xyzzyaaah106*xyzzyaaat106)*xyzzyaaat106
xyzzyaabk106=1.d0/(1.d0+xyzzyaaag106*xyzzyaaau106)
xyzzyaabs106=(xyzzyaaae106+xyzzyaaaf106*xyzzyaaau106)*xyzzyaaau106
xyzzyaabo106=1.d0/(1.d0+xyzzyaaai106*xyzzyaaau106)
xyzzyaabw106=(xyzzyaaae106+xyzzyaaah106*xyzzyaaau106)*xyzzyaaau106
xyzzyaabz106=xyzzyaabp106*xyzzyaabh106+xyzzyaabu106*xyzzyaabm106+xyzzy&
&aabv106*xyzzyaabn106+xyzzyaabs106*xyzzyaabk106+xyzzyaacx106
xyzzyaabx106=exp_limit(xyzzyaabz106)
xyzzyaabz106=xyzzyaabt106*xyzzyaabl106+xyzzyaabq106*xyzzyaabi106+xyzzy&
&aabr106*xyzzyaabj106+xyzzyaabw106*xyzzyaabo106+xyzzyaacx106
xyzzyaaby106=exp_limit(xyzzyaabz106)
xyzzyaaca106=xyzzyaabx106+xyzzyaaby106
if(val)val_jas=xyzzyaaaa106*xyzzyaaap106*xyzzyaabf106+xyzzyaacv106+xyz&
&zyaaac106*xyzzyaaaq106*xyzzyaabg106+xyzzyaacw106+log(xyzzyaaca106)
if(fd.or.sd)then
if(jspin==1)then
xyzzyaacu106=(xyzzyaaaa106*xyzzyaabf106)*(-1.d0/xyzzyaaap106+xyzzyaaab&
&106*xyzzyaabf106)
call xyzzyaaqm1(xyzzyaaar106,xyzzyaaaj106*xyzzyaaav106,xyzzyaaae106,xy&
&zzyaaaf106,xyzzyaaag106,xyzzyaacb106)
call xyzzyaaqm1(xyzzyaaas106,-xyzzyaabb106/xyzzyaaas106,xyzzyaaae106,x&
&yzzyaaah106,xyzzyaaai106,xyzzyaacd106)
call xyzzyaaqm1(xyzzyaaar106,xyzzyaaaj106*xyzzyaaav106,xyzzyaaae106,xy&
&zzyaaah106,xyzzyaaai106,xyzzyaacc106)
call xyzzyaaqm1(xyzzyaaas106,-xyzzyaabb106/xyzzyaaas106,xyzzyaaae106,x&
&yzzyaaaf106,xyzzyaaag106,xyzzyaace106)
grad_jas(1)=((xyzzyaacb106+xyzzyaacd106)*xyzzyaabx106+(xyzzyaacc106+xy&
&zzyaace106)*xyzzyaaby106)/xyzzyaaca106+xyzzyaacu106*xyzzyaaax106
call xyzzyaaqm1(xyzzyaaar106,xyzzyaaam106*xyzzyaaav106,xyzzyaaae106,xy&
&zzyaaaf106,xyzzyaaag106,xyzzyaacf106)
call xyzzyaaqm1(xyzzyaaas106,xyzzyaabc106/xyzzyaaas106,xyzzyaaae106,xy&
&zzyaaah106,xyzzyaaai106,xyzzyaach106)
call xyzzyaaqm1(xyzzyaaar106,xyzzyaaam106*xyzzyaaav106,xyzzyaaae106,xy&
&zzyaaah106,xyzzyaaai106,xyzzyaacg106)
call xyzzyaaqm1(xyzzyaaas106,xyzzyaabc106/xyzzyaaas106,xyzzyaaae106,xy&
&zzyaaaf106,xyzzyaaag106,xyzzyaaci106)
grad_jas(2)=((xyzzyaacf106+xyzzyaach106)*xyzzyaabx106+(xyzzyaacg106+xy&
&zzyaaci106)*xyzzyaaby106)/xyzzyaaca106+xyzzyaacu106*xyzzyaaay106
if(sd)then
call xyzzyaaqo1(xyzzyaaap106,-xyzzyaaax106/xyzzyaaap106,-(xyzzyaaax106&
&**2)/(xyzzyaaap106**3)+1/xyzzyaaap106,xyzzyaaaa106,xyzzyaaab106,xyzzy&
&aacj106)
call xyzzyaaqo1(xyzzyaaap106,-xyzzyaaay106/xyzzyaaap106,-(xyzzyaaay106&
&**2)/(xyzzyaaap106**3)+1/xyzzyaaap106,xyzzyaaaa106,xyzzyaaab106,xyzzy&
&aack106)
xyzzyaacs106=-(((xyzzyaacb106+xyzzyaacd106)*xyzzyaabx106+(xyzzyaacc106&
&+xyzzyaace106)*xyzzyaaby106)/xyzzyaaca106)**2
xyzzyaact106=-(xyzzyaaaj106**2)/(xyzzyaaar106**3)+xyzzyaaav106
call xyzzyaaqn1(xyzzyaaar106,xyzzyaaaj106*xyzzyaaav106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaaf106,xyzzyaaag106,xyzzyaacl106)
xyzzyaact106=-(xyzzyaabb106**2)/(xyzzyaaas106**3)+1/xyzzyaaas106
call xyzzyaaqn1(xyzzyaaas106,-xyzzyaabb106/xyzzyaaas106,xyzzyaact106,x&
&yzzyaaae106,xyzzyaaah106,xyzzyaaai106,xyzzyaacm106)
xyzzyaacp106=xyzzyaacl106+xyzzyaacm106
xyzzyaact106=-(xyzzyaaaj106**2)/(xyzzyaaar106**3)+xyzzyaaav106
call xyzzyaaqn1(xyzzyaaar106,xyzzyaaaj106*xyzzyaaav106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaah106,xyzzyaaai106,xyzzyaacl106)
xyzzyaact106=-(xyzzyaabb106**2)/(xyzzyaaas106**3)+1/xyzzyaaas106
call xyzzyaaqn1(xyzzyaaas106,-xyzzyaabb106/xyzzyaaas106,xyzzyaact106,x&
&yzzyaaae106,xyzzyaaaf106,xyzzyaaag106,xyzzyaacm106)
xyzzyaacq106=xyzzyaacl106+xyzzyaacm106
xyzzyaacs106=xyzzyaacs106+(((xyzzyaacb106+xyzzyaacd106)**2)*xyzzyaabx1&
&06+xyzzyaacp106*xyzzyaabx106+((xyzzyaacc106+xyzzyaace106)**2)*xyzzyaa&
&by106+xyzzyaacq106*xyzzyaaby106)/xyzzyaaca106
xyzzyaacr106=-(((xyzzyaacf106+xyzzyaach106)*xyzzyaabx106+(xyzzyaacg106&
&+xyzzyaaci106)*xyzzyaaby106)/xyzzyaaca106)**2
xyzzyaact106=-(xyzzyaaam106**2)/(xyzzyaaar106**3)+xyzzyaaav106
call xyzzyaaqn1(xyzzyaaar106,xyzzyaaam106*xyzzyaaav106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaaf106,xyzzyaaag106,xyzzyaacl106)
xyzzyaact106=-(xyzzyaabc106**2)/(xyzzyaaas106**3)+1/xyzzyaaas106
call xyzzyaaqn1(xyzzyaaas106,xyzzyaabc106/xyzzyaaas106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaah106,xyzzyaaai106,xyzzyaacm106)
xyzzyaacp106=xyzzyaacl106+xyzzyaacm106
xyzzyaact106=-(xyzzyaaam106**2)/(xyzzyaaar106**3)+xyzzyaaav106
call xyzzyaaqn1(xyzzyaaar106,xyzzyaaam106*xyzzyaaav106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaah106,xyzzyaaai106,xyzzyaacl106)
xyzzyaact106=-(xyzzyaabc106**2)/(xyzzyaaas106**3)+1/xyzzyaaas106
call xyzzyaaqn1(xyzzyaaas106,xyzzyaabc106/xyzzyaaas106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaaf106,xyzzyaaag106,xyzzyaacm106)
xyzzyaacq106=xyzzyaacl106+xyzzyaacm106
xyzzyaacr106=xyzzyaacr106+(((xyzzyaacf106+xyzzyaach106)**2)*xyzzyaabx1&
&06+xyzzyaacp106*xyzzyaabx106+((xyzzyaacg106+xyzzyaaci106)**2)*xyzzyaa&
&by106+xyzzyaacq106*xyzzyaaby106)/xyzzyaaca106
lap_jas=xyzzyaacj106+xyzzyaack106+xyzzyaacs106+xyzzyaacr106
endif
else
xyzzyaacu106=-(xyzzyaaaa106*xyzzyaabf106)*(-1.d0/xyzzyaaap106+xyzzyaaa&
&b106*xyzzyaabf106)
call xyzzyaaqm1(xyzzyaaau106,xyzzyaaak106*xyzzyaaaw106,xyzzyaaae106,xy&
&zzyaaaf106,xyzzyaaag106,xyzzyaacb106)
call xyzzyaaqm1(xyzzyaaat106,xyzzyaabd106/xyzzyaaat106,xyzzyaaae106,xy&
&zzyaaah106,xyzzyaaai106,xyzzyaacd106)
call xyzzyaaqm1(xyzzyaaau106,xyzzyaaak106*xyzzyaaaw106,xyzzyaaae106,xy&
&zzyaaah106,xyzzyaaai106,xyzzyaacc106)
call xyzzyaaqm1(xyzzyaaat106,xyzzyaabd106/xyzzyaaat106,xyzzyaaae106,xy&
&zzyaaaf106,xyzzyaaag106,xyzzyaace106)
grad_jas(1)=((xyzzyaacb106+xyzzyaacd106)*xyzzyaabx106+(xyzzyaacc106+xy&
&zzyaace106)*xyzzyaaby106)/xyzzyaaca106+xyzzyaacu106*xyzzyaaax106
call xyzzyaaqm1(xyzzyaaau106,xyzzyaaan106*xyzzyaaaw106,xyzzyaaae106,xy&
&zzyaaaf106,xyzzyaaag106,xyzzyaacf106)
call xyzzyaaqm1(xyzzyaaat106,xyzzyaabe106/xyzzyaaat106,xyzzyaaae106,xy&
&zzyaaah106,xyzzyaaai106,xyzzyaach106)
call xyzzyaaqm1(xyzzyaaau106,xyzzyaaan106*xyzzyaaaw106,xyzzyaaae106,xy&
&zzyaaah106,xyzzyaaai106,xyzzyaacg106)
call xyzzyaaqm1(xyzzyaaat106,xyzzyaabe106/xyzzyaaat106,xyzzyaaae106,xy&
&zzyaaaf106,xyzzyaaag106,xyzzyaaci106)
grad_jas(2)=((xyzzyaacf106+xyzzyaach106)*xyzzyaabx106+(xyzzyaacg106+xy&
&zzyaaci106)*xyzzyaaby106)/xyzzyaaca106+xyzzyaacu106*xyzzyaaay106
if(sd)then
call xyzzyaaqo1(xyzzyaaap106,xyzzyaaax106/xyzzyaaap106,-(xyzzyaaax106*&
&*2)/(xyzzyaaap106**3)+1/xyzzyaaap106,xyzzyaaaa106,xyzzyaaab106,xyzzya&
&acj106)
call xyzzyaaqo1(xyzzyaaap106,xyzzyaaay106/xyzzyaaap106,-(xyzzyaaay106*&
&*2)/(xyzzyaaap106**3)+1/xyzzyaaap106,xyzzyaaaa106,xyzzyaaab106,xyzzya&
&ack106)
xyzzyaacs106=-(((xyzzyaacb106+xyzzyaacd106)*xyzzyaabx106+(xyzzyaacc106&
&+xyzzyaace106)*xyzzyaaby106)/xyzzyaaca106)**2
xyzzyaact106=-(xyzzyaaak106**2)/(xyzzyaaau106**3)+xyzzyaaaw106
call xyzzyaaqn1(xyzzyaaau106,xyzzyaaak106*xyzzyaaaw106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaaf106,xyzzyaaag106,xyzzyaacn106)
xyzzyaact106=-(xyzzyaabd106**2)/(xyzzyaaat106**3)+1/xyzzyaaat106
call xyzzyaaqn1(xyzzyaaat106,xyzzyaabd106/xyzzyaaat106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaah106,xyzzyaaai106,xyzzyaaco106)
xyzzyaacp106=xyzzyaacn106+xyzzyaaco106
xyzzyaact106=-(xyzzyaaak106**2)/(xyzzyaaau106**3)+xyzzyaaaw106
call xyzzyaaqn1(xyzzyaaau106,xyzzyaaak106*xyzzyaaaw106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaah106,xyzzyaaai106,xyzzyaacn106)
xyzzyaact106=-(xyzzyaabd106**2)/(xyzzyaaat106**3)+1/xyzzyaaat106
call xyzzyaaqn1(xyzzyaaat106,xyzzyaabd106/xyzzyaaat106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaaf106,xyzzyaaag106,xyzzyaaco106)
xyzzyaacq106=xyzzyaacn106+xyzzyaaco106
xyzzyaacs106=xyzzyaacs106+(((xyzzyaacb106+xyzzyaacd106)**2)*xyzzyaabx1&
&06+xyzzyaacp106*xyzzyaabx106+((xyzzyaacc106+xyzzyaace106)**2)*xyzzyaa&
&by106+xyzzyaacq106*xyzzyaaby106)/xyzzyaaca106
xyzzyaacr106=-(((xyzzyaacf106+xyzzyaach106)*xyzzyaabx106+(xyzzyaacg106&
&+xyzzyaaci106)*xyzzyaaby106)/xyzzyaaca106)**2
xyzzyaact106=-(xyzzyaaan106**2)/(xyzzyaaau106**3)+xyzzyaaaw106
call xyzzyaaqn1(xyzzyaaau106,xyzzyaaan106*xyzzyaaaw106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaaf106,xyzzyaaag106,xyzzyaacn106)
xyzzyaact106=-(xyzzyaabe106**2)/(xyzzyaaat106**3)+1/xyzzyaaat106
call xyzzyaaqn1(xyzzyaaat106,xyzzyaabe106/xyzzyaaat106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaah106,xyzzyaaai106,xyzzyaaco106)
xyzzyaacp106=xyzzyaacn106+xyzzyaaco106
xyzzyaact106=-(xyzzyaaan106**2)/(xyzzyaaau106**3)+xyzzyaaaw106
call xyzzyaaqn1(xyzzyaaau106,xyzzyaaan106*xyzzyaaaw106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaah106,xyzzyaaai106,xyzzyaacn106)
xyzzyaact106=-(xyzzyaabe106**2)/(xyzzyaaat106**3)+1/xyzzyaaat106
call xyzzyaaqn1(xyzzyaaat106,xyzzyaabe106/xyzzyaaat106,xyzzyaact106,xy&
&zzyaaae106,xyzzyaaaf106,xyzzyaaag106,xyzzyaaco106)
xyzzyaacq106=xyzzyaacn106+xyzzyaaco106
xyzzyaacr106=xyzzyaacr106+(((xyzzyaacf106+xyzzyaach106)**2)*xyzzyaabx1&
&06+xyzzyaacp106*xyzzyaabx106+((xyzzyaacg106+xyzzyaaci106)**2)*xyzzyaa&
&by106+xyzzyaacq106*xyzzyaaby106)/xyzzyaaca106
lap_jas=xyzzyaacj106+xyzzyaack106+xyzzyaacs106+xyzzyaacr106
endif
endif
endif
end subroutine xyzzyaaql1
subroutine xyzzyaaqm1(r,r_prime,c5,c_a,c_b,dj_eh)
implicit none
real(dp),intent(in) :: r,r_prime,c5,c_a,c_b
real(dp),intent(out) :: dj_eh
dj_eh=-(c5*r+c_a*c_b*r_prime*r**2)/(1+c_b*r)**2
dj_eh=dj_eh+r_prime*(c5+2*c_a*r)/(1+c_b*r)
end subroutine xyzzyaaqm1
subroutine xyzzyaaqn1(r,r_prime,r_2prime,c5,c_a,c_b,d2j_eh)
implicit none
real(dp),intent(in) :: r,r_prime,r_2prime,c5,c_a,c_b
real(dp),intent(out) :: d2j_eh
real(dp) xyzzyaaaa108,xyzzyaaab108,xyzzyaaac108,xyzzyaaad108,xyzzyaaae&
&108
xyzzyaaae108=1.d0/(1+c_b*r)
xyzzyaaaa108=((-c5*r-c_a*(r**2))*c_b*r_2prime+c_b*r_prime*(-c5*r_prime&
&-2*c_a*r*r_prime))*(xyzzyaaae108**2)
xyzzyaaab108=(2*(c5*r+c_a*(r**2))*(c_b*r_prime)**2)*(xyzzyaaae108**3)
xyzzyaaac108=-(c5*r_prime+2*c_a*r*r_prime)*(c_b*r_prime)*(xyzzyaaae108&
&**2)
xyzzyaaad108=(c5*r_2prime+2*c_a*(r*r_2prime+(r_prime)**2))*xyzzyaaae10&
&8
d2j_eh=xyzzyaaaa108+xyzzyaaab108+xyzzyaaac108+xyzzyaaad108
end subroutine xyzzyaaqn1
subroutine xyzzyaaqo1(r,r_prime,r_2prime,c1,c2,d2j_ee)
implicit none
real(dp),intent(in) :: r,r_prime,r_2prime,c1,c2
real(dp),intent(out) :: d2j_ee
real(dp) xyzzyaaaa109
xyzzyaaaa109=1.d0/(1+c2*r)
d2j_ee=-c1*c2*r_prime*r_prime*xyzzyaaaa109*xyzzyaaaa109
d2j_ee=d2j_ee+c1*r_2prime*xyzzyaaaa109
d2j_ee=d2j_ee-(-2*c1*(c2**2)*r*(r_prime**2)*(xyzzyaaaa109**3)+c1*c2*r*&
&r_2prime*(xyzzyaaaa109**2)+c1*c2*(r_prime**2)*(xyzzyaaaa109**2))
end subroutine xyzzyaaqo1
subroutine xyzzyaaqp1(val,fd,sd,kk,eevecs1,eevecs,eivecs1,eivecs,val_j&
&as,grad_jas,lap_jas)
use slaarnabg,only : rion
implicit none
integer,intent(in) :: kk
real(dp),intent(out) :: val_jas,grad_jas(3),lap_jas
real(dp),intent(in) :: eevecs1(4,netot),eevecs(4,netot,netot),eivecs1(&
&4,nitot),eivecs(4,nitot,netot)
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa110,xyzzyaaab110,xyzzyaaac110,xyzzyaaad110,xyzzyaaae1&
&10,xyzzyaaaf110,xyzzyaaag110,xyzzyaaah110,xyzzyaaai110,xyzzyaaaj110,x&
&yzzyaaak110
real(dp) xyzzyaaal110,xyzzyaaam110
real(dp),allocatable,save :: xyzzyaaan110(:,:),xyzzyaaao110(:,:)
real(dp),parameter :: xyzzyaaap110=1.d-8
logical,save :: xyzzyaaaq110=.true.
if(xyzzyaaaq110)then
allocate(xyzzyaaan110(nspin,nspin),xyzzyaaao110(nitot,nspin),stat=xyzz&
&yaaae110)
if(xyzzyaaae110/=0)call errstop('EX2D_JAS','Allocation error.')
do xyzzyaaaa110=1,netot
xyzzyaaaf110=which_spin(xyzzyaaaa110)
xyzzyaaag110=heg_layer(xyzzyaaaf110)
do xyzzyaaab110=xyzzyaaaa110+1,netot
xyzzyaaah110=which_spin(xyzzyaaab110)
xyzzyaaai110=heg_layer(xyzzyaaah110)
if(xyzzyaaiu1)then
xyzzyaaan110(xyzzyaaah110,xyzzyaaaf110)=-pcharge(xyzzyaaah110)*pcharge&
&(xyzzyaaaf110)*pmass(xyzzyaaah110)*pmass(xyzzyaaaf110)/(2.d0*(rstar+a&
&bs(heg_zlayer(xyzzyaaag110)-heg_zlayer(xyzzyaaai110)))*(pmass(xyzzyaa&
&af110)+pmass(xyzzyaaah110)))
else
if(xyzzyaaag110==xyzzyaaai110)then
xyzzyaaan110(xyzzyaaah110,xyzzyaaaf110)=2.d0*pcharge(xyzzyaaah110)*pch&
&arge(xyzzyaaaf110)*pmass(xyzzyaaah110)*pmass(xyzzyaaaf110)/(pmass(xyz&
&zyaaaf110)+pmass(xyzzyaaah110))
else
xyzzyaaan110(xyzzyaaah110,xyzzyaaaf110)=0.d0
endif
endif
xyzzyaaan110(xyzzyaaaf110,xyzzyaaah110)=xyzzyaaan110(xyzzyaaah110,xyzz&
&yaaaf110)
enddo
do xyzzyaaac110=1,nitot
if(xyzzyaaiu1)then
xyzzyaaao110(xyzzyaaac110,xyzzyaaaf110)=-zion(iontype(xyzzyaaac110))*p&
&charge(xyzzyaaaf110)*pmass(xyzzyaaaf110)/(2.d0*(rstar+abs(heg_zlayer(&
&xyzzyaaag110)-rion(3,xyzzyaaac110))))
else
if(abs(heg_zlayer(xyzzyaaag110)-rion(3,xyzzyaaac110))<xyzzyaaap110)the&
&n
xyzzyaaao110(xyzzyaaac110,xyzzyaaaf110)=2.d0*zion(iontype(xyzzyaaac110&
&))*pcharge(xyzzyaaaf110)*pmass(xyzzyaaaf110)
else
xyzzyaaao110(xyzzyaaac110,xyzzyaaaf110)=0.d0
endif
endif
enddo
enddo
xyzzyaaaq110=.false.
endif
val_jas=0.d0
grad_jas=0.d0
lap_jas=0.d0
if(xyzzyaaiu1)then
if(val)then
do xyzzyaaaa110=1,netot
xyzzyaaaf110=which_spin(xyzzyaaaa110)
xyzzyaaag110=heg_layer(xyzzyaaaf110)
do xyzzyaaab110=xyzzyaaaa110+1,netot
xyzzyaaah110=which_spin(xyzzyaaab110)
xyzzyaaal110=eevecs(4,xyzzyaaab110,xyzzyaaaa110)
xyzzyaaad110=which_spair(xyzzyaaah110,xyzzyaaaf110,xyzzyaais1)
val_jas=val_jas+(xyzzyaaii1(1,xyzzyaaad110)+xyzzyaaan110(xyzzyaaah110,&
&xyzzyaaaf110)*log(xyzzyaaal110)+xyzzyaaii1(2,xyzzyaaad110)*xyzzyaaal1&
&10)*xyzzyaaal110**2/(1.d0+xyzzyaaii1(3,xyzzyaaad110)*xyzzyaaal110**2)
enddo
xyzzyaaad110=no_spairs(xyzzyaais1)+(which_ssingle(xyzzyaaaf110,xyzzyaa&
&it1)-1)*nitot
if(xyzzyaaaa110==kk)then
do xyzzyaaac110=1,nitot
xyzzyaaal110=eivecs1(4,xyzzyaaac110)
xyzzyaaad110=xyzzyaaad110+1
val_jas=val_jas+(xyzzyaaii1(1,xyzzyaaad110)+xyzzyaaao110(xyzzyaaac110,&
&xyzzyaaaf110)*log(xyzzyaaal110)+xyzzyaaii1(2,xyzzyaaad110)*xyzzyaaal1&
&10)*xyzzyaaal110**2/(1.d0+xyzzyaaii1(3,xyzzyaaad110)*xyzzyaaal110**2)
enddo
else
do xyzzyaaac110=1,nitot
xyzzyaaal110=eivecs(4,xyzzyaaac110,xyzzyaaaa110)
xyzzyaaad110=xyzzyaaad110+1
val_jas=val_jas+(xyzzyaaii1(1,xyzzyaaad110)+xyzzyaaao110(xyzzyaaac110,&
&xyzzyaaaf110)*log(xyzzyaaal110)+xyzzyaaii1(2,xyzzyaaad110)*xyzzyaaal1&
&10)*xyzzyaaal110**2/(1.d0+xyzzyaaii1(3,xyzzyaaad110)*xyzzyaaal110**2)
enddo
endif
enddo
endif
if(fd.or.sd)then
xyzzyaaaj110=which_spin(kk)
xyzzyaaak110=heg_layer(xyzzyaaaj110)
do xyzzyaaaa110=1,netot
if(xyzzyaaaa110==kk)cycle
xyzzyaaaf110=which_spin(xyzzyaaaa110)
xyzzyaaal110=eevecs1(4,xyzzyaaaa110)
xyzzyaaam110=log(xyzzyaaal110)
xyzzyaaad110=which_spair(xyzzyaaaf110,xyzzyaaaj110,xyzzyaais1)
if(fd)grad_jas(1:2)=grad_jas(1:2)+(((xyzzyaaii1(3,xyzzyaaad110)*(xyzzy&
&aaii1(2,xyzzyaaad110)*xyzzyaaal110+xyzzyaaan110(xyzzyaaaf110,xyzzyaaa&
&j110))*xyzzyaaal110+3.d0*xyzzyaaii1(2,xyzzyaaad110))*xyzzyaaal110+2.d&
&0*(xyzzyaaan110(xyzzyaaaf110,xyzzyaaaj110)*xyzzyaaam110+xyzzyaaii1(1,&
&xyzzyaaad110))+xyzzyaaan110(xyzzyaaaf110,xyzzyaaaj110))/(1.d0+xyzzyaa&
&ii1(3,xyzzyaaad110)*xyzzyaaal110*xyzzyaaal110)**2)*eevecs1(1:2,xyzzya&
&aaa110)
if(sd)lap_jas=lap_jas+(((xyzzyaaii1(2,xyzzyaaad110)*xyzzyaaii1(3,xyzzy&
&aaad110)*(xyzzyaaii1(3,xyzzyaaad110)*xyzzyaaal110*xyzzyaaal110+2.d0)*&
&xyzzyaaal110-4.d0*(xyzzyaaii1(3,xyzzyaaad110)*xyzzyaaan110(xyzzyaaaf1&
&10,xyzzyaaaj110)*xyzzyaaam110+xyzzyaaii1(1,xyzzyaaad110)*xyzzyaaii1(3&
&,xyzzyaaad110)))*xyzzyaaal110+4.d0*xyzzyaaii1(3,xyzzyaaad110)*xyzzyaa&
&an110(xyzzyaaaf110,xyzzyaaaj110)*xyzzyaaal110+9.d0*xyzzyaaii1(2,xyzzy&
&aaad110))*xyzzyaaal110+4.d0*(xyzzyaaan110(xyzzyaaaf110,xyzzyaaaj110)*&
&(xyzzyaaam110+1.d0)+xyzzyaaii1(1,xyzzyaaad110)))/(1.d0+xyzzyaaii1(3,x&
&yzzyaaad110)*xyzzyaaal110*xyzzyaaal110)**3
enddo
xyzzyaaad110=no_spairs(xyzzyaais1)+(which_ssingle(xyzzyaaaj110,xyzzyaa&
&it1)-1)*nitot
do xyzzyaaac110=1,nitot
xyzzyaaal110=eivecs1(4,xyzzyaaac110)
xyzzyaaam110=log(xyzzyaaal110)
xyzzyaaad110=xyzzyaaad110+1
if(fd)grad_jas(1:2)=grad_jas(1:2)+(((xyzzyaaii1(3,xyzzyaaad110)*(xyzzy&
&aaii1(2,xyzzyaaad110)*xyzzyaaal110+xyzzyaaao110(xyzzyaaac110,xyzzyaaa&
&j110))*xyzzyaaal110+3.d0*xyzzyaaii1(2,xyzzyaaad110))*xyzzyaaal110+2.d&
&0*(xyzzyaaao110(xyzzyaaac110,xyzzyaaaj110)*xyzzyaaam110+xyzzyaaii1(1,&
&xyzzyaaad110))+xyzzyaaao110(xyzzyaaac110,xyzzyaaaj110))/(1.d0+xyzzyaa&
&ii1(3,xyzzyaaad110)*xyzzyaaal110*xyzzyaaal110)**2)*eivecs1(1:2,xyzzya&
&aac110)
if(sd)lap_jas=lap_jas+(((xyzzyaaii1(2,xyzzyaaad110)*xyzzyaaii1(3,xyzzy&
&aaad110)*(xyzzyaaii1(3,xyzzyaaad110)*xyzzyaaal110*xyzzyaaal110+2.d0)*&
&xyzzyaaal110-4.d0*(xyzzyaaii1(3,xyzzyaaad110)*xyzzyaaao110(xyzzyaaac1&
&10,xyzzyaaaj110)*xyzzyaaam110+xyzzyaaii1(1,xyzzyaaad110)*xyzzyaaii1(3&
&,xyzzyaaad110)))*xyzzyaaal110+4.d0*xyzzyaaii1(3,xyzzyaaad110)*xyzzyaa&
&ao110(xyzzyaaac110,xyzzyaaaj110)*xyzzyaaal110+9.d0*xyzzyaaii1(2,xyzzy&
&aaad110))*xyzzyaaal110+4.d0*(xyzzyaaao110(xyzzyaaac110,xyzzyaaaj110)*&
&(xyzzyaaam110+1.d0)+xyzzyaaii1(1,xyzzyaaad110)))/(1.d0+xyzzyaaii1(3,x&
&yzzyaaad110)*xyzzyaaal110*xyzzyaaal110)**3
enddo
endif
else
if(val)then
do xyzzyaaaa110=1,netot
xyzzyaaaf110=which_spin(xyzzyaaaa110)
xyzzyaaag110=heg_layer(xyzzyaaaf110)
do xyzzyaaab110=xyzzyaaaa110+1,netot
xyzzyaaah110=which_spin(xyzzyaaab110)
xyzzyaaal110=eevecs(4,xyzzyaaab110,xyzzyaaaa110)
xyzzyaaad110=which_spair(xyzzyaaah110,xyzzyaaaf110,xyzzyaais1)
val_jas=val_jas+(xyzzyaaan110(xyzzyaaah110,xyzzyaaaf110)+xyzzyaaii1(1,&
&xyzzyaaad110)*xyzzyaaal110)*xyzzyaaal110/(1.d0+xyzzyaaii1(2,xyzzyaaad&
&110)*xyzzyaaal110)
enddo
xyzzyaaad110=no_spairs(xyzzyaais1)+(which_ssingle(xyzzyaaaf110,xyzzyaa&
&it1)-1)*nitot
if(xyzzyaaaa110==kk)then
do xyzzyaaac110=1,nitot
xyzzyaaal110=eivecs1(4,xyzzyaaac110)
xyzzyaaad110=xyzzyaaad110+1
val_jas=val_jas+(xyzzyaaao110(xyzzyaaac110,xyzzyaaaf110)+xyzzyaaii1(1,&
&xyzzyaaad110)*xyzzyaaal110)*xyzzyaaal110 /(1.d0+xyzzyaaii1(2,xyzzyaaa&
&d110)*xyzzyaaal110)
enddo
else
do xyzzyaaac110=1,nitot
xyzzyaaal110=eivecs(4,xyzzyaaac110,xyzzyaaaa110)
xyzzyaaad110=xyzzyaaad110+1
val_jas=val_jas+(xyzzyaaao110(xyzzyaaac110,xyzzyaaaf110)+xyzzyaaii1(1,&
&xyzzyaaad110)*xyzzyaaal110)*xyzzyaaal110/(1.d0+xyzzyaaii1(2,xyzzyaaad&
&110)*xyzzyaaal110)
enddo
endif
enddo
endif
if(fd.or.sd)then
xyzzyaaaj110=which_spin(kk)
xyzzyaaak110=heg_layer(xyzzyaaaj110)
do xyzzyaaaa110=1,netot
if(xyzzyaaaa110==kk)cycle
xyzzyaaaf110=which_spin(xyzzyaaaa110)
xyzzyaaal110=eevecs1(4,xyzzyaaaa110)
xyzzyaaad110=which_spair(xyzzyaaaj110,xyzzyaaaf110,xyzzyaais1)
if(fd)grad_jas(1:2)=grad_jas(1:2)+(xyzzyaaan110(xyzzyaaaf110,xyzzyaaaj&
&110)+xyzzyaaii1(1,xyzzyaaad110)*xyzzyaaal110*(2.d0+xyzzyaaii1(2,xyzzy&
&aaad110)*xyzzyaaal110))/(xyzzyaaal110*(1.d0+xyzzyaaii1(2,xyzzyaaad110&
&)*xyzzyaaal110)**2)*eevecs1(1:2,xyzzyaaaa110)
if(sd)lap_jas=lap_jas+(xyzzyaaan110(xyzzyaaaf110,xyzzyaaaj110)-xyzzyaa&
&ii1(2,xyzzyaaad110)*xyzzyaaan110(xyzzyaaaf110,xyzzyaaaj110)*xyzzyaaal&
&110+xyzzyaaii1(1,xyzzyaaad110)*xyzzyaaal110*(4.d0+3.d0*xyzzyaaii1(2,x&
&yzzyaaad110)*xyzzyaaal110+xyzzyaaii1(2,xyzzyaaad110)**2*xyzzyaaal110*&
&xyzzyaaal110))/(xyzzyaaal110*(1.d0+xyzzyaaii1(2,xyzzyaaad110)*xyzzyaa&
&al110)**3)
enddo
xyzzyaaad110=no_spairs(xyzzyaais1)+(which_ssingle(xyzzyaaaj110,xyzzyaa&
&it1)-1) *nitot
do xyzzyaaac110=1,nitot
xyzzyaaal110=eivecs1(4,xyzzyaaac110)
xyzzyaaad110=xyzzyaaad110+1
if(fd)grad_jas(1:2)=grad_jas(1:2)+(xyzzyaaao110(xyzzyaaac110,xyzzyaaaj&
&110)+xyzzyaaii1(1,xyzzyaaad110)*xyzzyaaal110*(2.d0+xyzzyaaii1(2,xyzzy&
&aaad110)*xyzzyaaal110))/(xyzzyaaal110*(1.d0+xyzzyaaii1(2,xyzzyaaad110&
&)*xyzzyaaal110)**2)*eivecs1(1:2,xyzzyaaac110)
if(sd)lap_jas=lap_jas+(xyzzyaaao110(xyzzyaaac110,xyzzyaaaj110)-xyzzyaa&
&ii1(2,xyzzyaaad110)*xyzzyaaao110(xyzzyaaac110,xyzzyaaaj110)*xyzzyaaal&
&110+xyzzyaaii1(1,xyzzyaaad110)*xyzzyaaal110*(4.d0+3.d0*xyzzyaaii1(2,x&
&yzzyaaad110)*xyzzyaaal110+xyzzyaaii1(2,xyzzyaaad110)**2*xyzzyaaal110*&
&xyzzyaaal110))/(xyzzyaaal110*(1.d0+xyzzyaaii1(2,xyzzyaaad110)*xyzzyaa&
&al110)**3)
enddo
endif
endif
end subroutine xyzzyaaqp1
subroutine xyzzyaaqq1(rii,rjj,pair,derivearg,val,fd,sd,value_d,grad_d,&
&lap_d)
implicit none
real(dp),intent(in) :: rii(4),rjj(4)
integer,intent(in) :: pair
real(dp),intent(out) :: value_d,grad_d(3),lap_d
logical,intent(in) :: val,fd,sd
integer,intent(in) :: derivearg
real(dp) xyzzyaaaa111(4),xyzzyaaab111
real(dp) xyzzyaaac111(0:xyzzyaacd1,0:xyzzyaacd1)
real(dp) xyzzyaaad111
real(dp) xyzzyaaae111
real(dp) xyzzyaaaf111
real(dp) xyzzyaaag111
real(dp) xyzzyaaah111
real(dp) xyzzyaaai111
real(dp) xyzzyaaaj111
real(dp) xyzzyaaak111
real(dp) xyzzyaaal111
real(dp) xyzzyaaam111
real(dp) xyzzyaaan111
real(dp) xyzzyaaao111,xyzzyaaap111,xyzzyaaaq111
real(dp) xyzzyaaar111,xyzzyaaas111,xyzzyaaat111
logical xyzzyaaau111(2)
integer xyzzyaaav111,xyzzyaaaw111,xyzzyaaax111,xyzzyaaay111,xyzzyaaaz1&
&11
value_d=0.d0
grad_d=0.d0
lap_d=0.d0
do xyzzyaaax111=1,2
xyzzyaaau111(xyzzyaaax111)=rii(4)>=xyzzyaacw1(pair,xyzzyaaax111,1).or.&
&rjj(4)>=xyzzyaacw1(pair,xyzzyaaax111,2)
enddo
if(all(xyzzyaaau111))return
xyzzyaaav111=xyzzyaagx1
xyzzyaaaw111=xyzzyaabv1(pair)
xyzzyaaaa111=rionion(:,xyzzyaabu1(pair,1),xyzzyaabu1(pair,2))
xyzzyaaad111=(rii(1)*xyzzyaaaa111(1)+rii(2)*xyzzyaaaa111(2)+rii(3)*xyz&
&zyaaaa111(3))
xyzzyaaae111=(rjj(1)*xyzzyaaaa111(1)+rjj(2)*xyzzyaaaa111(2)+rjj(3)*xyz&
&zyaaaa111(3))
xyzzyaaac111(0,0)=1.d0
do xyzzyaaaz111=1,xyzzyaacd1
xyzzyaaac111(0,xyzzyaaaz111)=xyzzyaaac111(0,xyzzyaaaz111-1)*rjj(4)
enddo
do xyzzyaaay111=1,xyzzyaacd1
xyzzyaaac111(xyzzyaaay111,0)=xyzzyaaac111(xyzzyaaay111-1,0)*rii(4)
do xyzzyaaaz111=1,xyzzyaacd1
xyzzyaaac111(xyzzyaaay111,xyzzyaaaz111)=xyzzyaaac111(xyzzyaaay111,xyzz&
&yaaaz111-1)*rjj(4)
enddo
enddo
do xyzzyaaax111=1,2
if(xyzzyaaau111(xyzzyaaax111))cycle
xyzzyaaag111=rii(4)-xyzzyaacw1(pair,xyzzyaaax111,1)
xyzzyaaah111=rjj(4)-xyzzyaacw1(pair,xyzzyaaax111,2)
if(derivearg==1)then
xyzzyaaam111=xyzzyaaag111**(xyzzyaaav111-2)
xyzzyaaak111=xyzzyaaam111*xyzzyaaag111
xyzzyaaai111=xyzzyaaak111*xyzzyaaag111
xyzzyaaaj111=xyzzyaaah111**xyzzyaaav111
else
xyzzyaaai111=xyzzyaaag111**xyzzyaaav111
xyzzyaaan111=xyzzyaaah111**(xyzzyaaav111-2)
xyzzyaaal111=xyzzyaaan111*xyzzyaaah111
xyzzyaaaj111=xyzzyaaal111*xyzzyaaah111
endif
xyzzyaaar111=0.d0
do xyzzyaaay111=0,xyzzyaacd1
do xyzzyaaaz111=0,xyzzyaacd1
xyzzyaaar111=xyzzyaaar111+xyzzyaaas1(xyzzyaaaw111,xyzzyaaax111,xyzzyaa&
&ay111,xyzzyaaaz111)*xyzzyaaac111(xyzzyaaay111,xyzzyaaaz111)
enddo
enddo
xyzzyaaao111=xyzzyaaai111*xyzzyaaaj111*xyzzyaaar111
if(fd.or.sd)then
if(derivearg==1)then
xyzzyaaas111=0.d0
do xyzzyaaay111=1,xyzzyaacd1
do xyzzyaaaz111=0,xyzzyaacd1
xyzzyaaas111=xyzzyaaas111+xyzzyaaay111*xyzzyaaas1(xyzzyaaaw111,xyzzyaa&
&ax111,xyzzyaaay111,xyzzyaaaz111)*xyzzyaaac111(xyzzyaaay111-1,xyzzyaaa&
&z111)
enddo
enddo
xyzzyaaap111=xyzzyaaav111*xyzzyaaak111*xyzzyaaaj111*xyzzyaaar111+xyzzy&
&aaai111*xyzzyaaaj111*xyzzyaaas111
if(sd)then
xyzzyaaat111=0.d0
do xyzzyaaay111=2,xyzzyaacd1
do xyzzyaaaz111=0,xyzzyaacd1
xyzzyaaat111=xyzzyaaat111+xyzzyaaay111*(xyzzyaaay111-1)*xyzzyaaas1(xyz&
&zyaaaw111,xyzzyaaax111,xyzzyaaay111,xyzzyaaaz111)*xyzzyaaac111(xyzzya&
&aay111-2,xyzzyaaaz111)
enddo
enddo
xyzzyaaaq111=xyzzyaaav111*(xyzzyaaav111-1)*xyzzyaaam111*xyzzyaaaj111*x&
&yzzyaaar111+2*xyzzyaaav111*xyzzyaaak111*xyzzyaaaj111*xyzzyaaas111+xyz&
&zyaaai111*xyzzyaaaj111*xyzzyaaat111
endif
else
xyzzyaaas111=0.d0
do xyzzyaaay111=0,xyzzyaacd1
do xyzzyaaaz111=1,xyzzyaacd1
xyzzyaaas111=xyzzyaaas111+xyzzyaaaz111*xyzzyaaas1(xyzzyaaaw111,xyzzyaa&
&ax111,xyzzyaaay111,xyzzyaaaz111)*xyzzyaaac111(xyzzyaaay111,xyzzyaaaz1&
&11-1)
enddo
enddo
xyzzyaaap111=xyzzyaaai111*xyzzyaaav111*xyzzyaaal111*xyzzyaaar111+xyzzy&
&aaai111*xyzzyaaaj111*xyzzyaaas111
if(sd)then
xyzzyaaat111=0.d0
do xyzzyaaay111=0,xyzzyaacd1
do xyzzyaaaz111=2,xyzzyaacd1
xyzzyaaat111=xyzzyaaat111+xyzzyaaaz111*(xyzzyaaaz111-1)*xyzzyaaas1(xyz&
&zyaaaw111,xyzzyaaax111,xyzzyaaay111,xyzzyaaaz111)*xyzzyaaac111(xyzzya&
&aay111,xyzzyaaaz111-2)
enddo
enddo
xyzzyaaaq111=xyzzyaaai111*xyzzyaaav111*(xyzzyaaav111-1)*xyzzyaaan111*x&
&yzzyaaar111+2*xyzzyaaai111*xyzzyaaav111*xyzzyaaal111*xyzzyaaas111+xyz&
&zyaaai111*xyzzyaaaj111*xyzzyaaat111
endif
endif
endif
if(xyzzyaaax111==1)then
if(val)then
value_d=value_d+xyzzyaaad111*xyzzyaaae111*xyzzyaaao111
endif
if(fd)then
if(derivearg==1)then
grad_d=grad_d+xyzzyaaaa111(1:3)*xyzzyaaae111*xyzzyaaao111+xyzzyaaad111&
&*xyzzyaaae111*rii(1:3)/rii(4)*xyzzyaaap111
else
grad_d=grad_d+xyzzyaaaa111(1:3)*xyzzyaaad111*xyzzyaaao111+xyzzyaaad111&
&*xyzzyaaae111*rjj(1:3)/rjj(4)*xyzzyaaap111
endif
endif
if(sd)then
if(derivearg==1)then
lap_d=lap_d+xyzzyaaad111*xyzzyaaae111*((dimensionality+1)/rii(4)*xyzzy&
&aaap111+xyzzyaaaq111)
else
lap_d=lap_d+xyzzyaaad111*xyzzyaaae111*((dimensionality+1)/rjj(4)*xyzzy&
&aaap111+xyzzyaaaq111)
endif
endif
else
xyzzyaaaf111=(rii(1)*rjj(1)+rii(2)*rjj(2)+rii(3)*rjj(3))
xyzzyaaab111=xyzzyaaaa111(4)*xyzzyaaaa111(4)
if(val)value_d=value_d+(xyzzyaaaf111*xyzzyaaab111-xyzzyaaad111*xyzzyaa&
&ae111)*xyzzyaaao111
if(derivearg==1)then
if(fd)grad_d=grad_d+(rjj(1:3)*xyzzyaaab111-xyzzyaaaa111(1:3)*xyzzyaaae&
&111)*xyzzyaaao111+(xyzzyaaaf111*xyzzyaaab111-xyzzyaaad111*xyzzyaaae11&
&1)*rii(1:3)/rii(4)*xyzzyaaap111
if(sd)lap_d=lap_d+(xyzzyaaaf111*xyzzyaaab111-xyzzyaaad111*xyzzyaaae111&
&)*((dimensionality+1)/rii(4)*xyzzyaaap111+xyzzyaaaq111)
else
if(fd)grad_d=grad_d+(rii(1:3)*xyzzyaaab111-xyzzyaaaa111(1:3)*xyzzyaaad&
&111)*xyzzyaaao111+(xyzzyaaaf111*xyzzyaaab111-xyzzyaaad111*xyzzyaaae11&
&1)*rjj(1:3)/rjj(4)*xyzzyaaap111
if(sd)lap_d=lap_d+(xyzzyaaaf111*xyzzyaaab111-xyzzyaaad111*xyzzyaaae111&
&)*((dimensionality+1)/rjj(4)*xyzzyaaap111+xyzzyaaaq111)
endif
endif
enddo
if(xyzzyaaaa1)then
call wout('compute_D: rIi,rIj=',(/rii,rjj/))
if(val)call wout('compute_D: value=',value_d)
if(fd)call wout('compute_D: grad=',grad_d)
if(sd)call wout('compute_D: lap=',lap_d)
endif
end subroutine xyzzyaaqq1
real(dp) function xyzzyaaqr1()
implicit none
real(dp),parameter :: xyzzyaaaa112=0.999999d0
if(isperiodic)then
xyzzyaaqr1=xyzzyaaaa112*wigner_seitz_radius
else
if(xyzzyaagx1>=2)then
if(nitot==1)then
xyzzyaaqr1=2.d0
else
xyzzyaaqr1=5.d0
endif
elseif(xyzzyaagx1==1)then
xyzzyaaqr1=1.d2
else
xyzzyaaqr1=1.d0
endif
endif
end function xyzzyaaqr1
real(dp) function xyzzyaaqs1()
implicit none
real(dp),parameter :: xyzzyaaaa113=0.999999d0
if(periodicity==2)then
xyzzyaaqs1=xyzzyaaaa113*wigner_seitz_radius
elseif(periodicity==0.or.periodicity==1)then
xyzzyaaqs1=5.d0
else
xyzzyaaqs1=min(2.d0,xyzzyaaaa113*wigner_seitz_radius)
endif
end function xyzzyaaqs1
real(dp) function xyzzyaaqt1()
implicit none
real(dp),parameter :: xyzzyaaaa114=0.999999d0
if(periodicity==2.or.periodicity==0)then
xyzzyaaqt1=5.d0
elseif(periodicity==1)then
xyzzyaaqt1=xyzzyaaaa114*wigner_seitz_radius
else
xyzzyaaqt1=min(2.d0,xyzzyaaaa114*wigner_seitz_radius)
endif
end function xyzzyaaqt1
real(dp) function xyzzyaaqu1()
implicit none
real(dp),parameter :: xyzzyaaaa115=0.999999d0
if(isperiodic)then
xyzzyaaqu1=xyzzyaaaa115*wigner_seitz_radius
else
if(xyzzyaagx1>=2)then
if(nitot==1)then
xyzzyaaqu1=2.d0
else
xyzzyaaqu1=5.d0
endif
elseif(xyzzyaagx1==1)then
xyzzyaaqu1=1.d2
else
xyzzyaaqu1=1.d0
endif
endif
end function xyzzyaaqu1
real(dp) function xyzzyaaqv1()
implicit none
real(dp),parameter :: xyzzyaaaa116=0.999999d0
if(isperiodic)then
xyzzyaaqv1=xyzzyaaaa116*wigner_seitz_radius
else
if(xyzzyaagx1>=2)then
if(nitot==1)then
xyzzyaaqv1=2.d0
else
xyzzyaaqv1=5.d0
endif
elseif(xyzzyaagx1==1)then
xyzzyaaqv1=1.d2
else
xyzzyaaqv1=1.d0
endif
endif
end function xyzzyaaqv1
real(dp) function xyzzyaaqw1()
implicit none
real(dp),parameter :: xyzzyaaaa117=0.999999d0
if(isperiodic)then
xyzzyaaqw1=min(xyzzyaaaa117*wigner_seitz_radius,4.d0)
else
if(xyzzyaagx1>=2)then
xyzzyaaqw1=4.d0
elseif(xyzzyaagx1==1)then
xyzzyaaqw1=1.d2
else
xyzzyaaqw1=1.d0
endif
endif
end function xyzzyaaqw1
real(dp) function xyzzyaaqx1()
implicit none
real(dp),parameter :: xyzzyaaaa118=0.999999d0
if(isperiodic)then
xyzzyaaqx1=min(0.5d0*xyzzyaaaa118*wigner_seitz_radius,3.d0)
else
if(xyzzyaagx1>=2)then
xyzzyaaqx1=3.d0
elseif(xyzzyaagx1==1)then
xyzzyaaqx1=1.d2
else
xyzzyaaqx1=1.d0
endif
endif
end function xyzzyaaqx1
real(dp) function xyzzyaaqy1(set)
use slaarnaat,only : cusp_correction
use slaarnabg,only : atom_basis_type
use slaarnabi,only : use_gpcc
use slaarnaca,only : is_ae
implicit none
integer,intent(in) :: set
integer xyzzyaaaa119
real(dp),parameter :: xyzzyaaab119=1.d-7
xyzzyaaqy1=zion(iontype(xyzzyaabq1(1,set)))
do xyzzyaaaa119=2,xyzzyaadi1(set)
if(abs(zion(iontype(xyzzyaabq1(xyzzyaaaa119,set)))-xyzzyaaqy1)>xyzzyaa&
&ab119)call errstop_master('COMPUTE_Z_NUC','All atoms in a given set s&
&hould have the same charge.')
enddo
if(.not.is_ae(xyzzyaabq1(1,set)))then
call errstop_master('COMPUTE_Z_NUC','You are attempting to enforce the&
& electron-nucleus cusp conditions for set '//trim(i2s(set))//' of the&
& chi term, but the ions in this set are described by pseudopotentials&
&.')
elseif((trim(atom_basis_type)=='gaussian'.and.cusp_correction) .or.(tr&
&im(atom_basis_type)/='none'.and.use_gpcc) .or.trim(atom_basis_type)==&
&'slater' .or.trim(atom_basis_type)=='numerical')then
call errstop_master('COMPUTE_Z_NUC','You are attempting to enforce the&
& electron-nucleus cusp conditions for set '//trim(i2s(set)) //' of th&
&e chi term, but you already have orbitals that satisfy the cusp condi&
&tions.')
endif
end function xyzzyaaqy1
subroutine xyzzyaaqz1
use slaarnaan, only : ee_kato_gamma,en_kato_gamma
use slaarnaas, only : harmwire_b
implicit none
integer xyzzyaaaa120,xyzzyaaab120,xyzzyaaac120
real(dp) xyzzyaaad120,ri,rj,xyzzyaaae120,xyzzyaaaf120,xyzzyaaag120,xyz&
&zyaaah120,xyzzyaaai120,xyzzyaaaj120,xyzzyaaak120,xyzzyaaal120,xyzzyaa&
&am120,xyzzyaaan120,xyzzyaaao120,xyzzyaaap120,xyzzyaaaq120,xyzzyaaar12&
&0,xyzzyaaas120,xyzzyaaat120
real(dp),parameter :: xyzzyaaau120=1.d-8
character(80) tmpr
if(am_master.and..not.hard_sphere)then
if(xyzzyaaep1)then
xyzzyaaad120=0.d0
do xyzzyaaaa120=1,nspin
do xyzzyaaab120=xyzzyaaaa120,nspin
call xyzzyaapr1(xyzzyaaad120,xyzzyaaaa120,xyzzyaaab120,.false.,.true.,&
&.false.,xyzzyaaae120,xyzzyaaaf120,xyzzyaaag120)
xyzzyaaah120=0.d0
if(.not.ee_cusp_in_orbital(xyzzyaaaa120,xyzzyaaab120).and..not.xyzzyaa&
&fa1.and..not.xyzzyaafc1.and..not.xyzzyaafg1.and.harmwire_b<0.d0)xyzzy&
&aaah120=ee_kato_gamma(xyzzyaaaa120,xyzzyaaab120,.not.noncoll_spin.and&
&.(xyzzyaafh1/=0.or.ferromagnetic))
if(abs(xyzzyaaaf120-xyzzyaaah120)>xyzzyaaau120)then
call wout('ispin    = '//trim(i2s(xyzzyaaaa120)))
call wout('jspin    = '//trim(i2s(xyzzyaaab120)))
tmpr=r2s(pmass(xyzzyaaaa120),'(f21.12)')
call wout('m_i      = '//trim(tmpr))
tmpr=r2s(pmass(xyzzyaaab120),'(f21.12)')
call wout('m_j      = '//trim(tmpr))
call wout('r_ij     = 0')
tmpr=r2s(xyzzyaaaf120,'(f21.12)')
call wout('du/dr_ij = '//trim(tmpr))
tmpr=r2s(xyzzyaaah120,'(f21.12)')
call wout('Gamma    = '//trim(tmpr))
call wout('So cusp condition on u not satisfied.  The derivative shoul&
&d equal Gamma.')
call errstop('TEST_CUSPS_AND_CUTOFFS','Stopping.')
endif
enddo
enddo
if(xyzzyaagx1>=1)then
xyzzyaaad120=xyzzyaach1
do xyzzyaaaa120=1,nspin
do xyzzyaaab120=1,nspin
call xyzzyaapr1(xyzzyaaad120,xyzzyaaaa120,xyzzyaaab120,.true.,.true.,.&
&true.,xyzzyaaae120,xyzzyaaaf120,xyzzyaaag120)
if(abs(xyzzyaaae120)>xyzzyaaau120.or.(abs(xyzzyaaaf120)>xyzzyaaau120.a&
&nd.xyzzyaagx1>=2).or.(abs(xyzzyaaag120)>xyzzyaaau120.and.xyzzyaagx1>=&
&3))then
call wout('ispin                    = '//trim(i2s(xyzzyaaaa120)))
call wout('jspin                    = '//trim(i2s(xyzzyaaab120)))
tmpr=r2s(xyzzyaaae120,'(f21.12)')
call wout('When r_ij=L_u, u         = '//trim(tmpr))
tmpr=r2s(xyzzyaaaf120,'(f21.12)')
call wout('              du/dr_ij   = '//trim(tmpr))
tmpr=r2s(xyzzyaaag120,'(f21.12)')
call wout('             d2u/dr_ij^2 = '//trim(tmpr))
call wout('So smooth cutoff conditions on u are not satisfied.')
call errstop('TEST_CUSPS_AND_CUTOFFS','Stopping.')
endif
enddo
enddo
endif
endif
if(xyzzyaaew1)then
do xyzzyaaac120=1,xyzzyaabi1
ri=0.d0
do xyzzyaaaa120=1,nspin
call xyzzyaaqa1(ri,xyzzyaaaa120,xyzzyaaac120,.false.,.true.,.false.,xy&
&zzyaaai120,xyzzyaaaj120,xyzzyaaak120)
if(xyzzyaadn1(xyzzyaaac120)==1)then
xyzzyaaah120=en_kato_gamma(xyzzyaaaa120,xyzzyaafx1(xyzzyaaac120))
else
xyzzyaaah120=0.d0
endif
if(abs(xyzzyaaaj120-xyzzyaaah120)>xyzzyaaau120)then
call wout('Chi set   = '//trim(i2s(xyzzyaaac120)))
call wout('ispin     = '//trim(i2s(xyzzyaaaa120)))
call wout('r_i       = 0')
tmpr=r2s(xyzzyaaaj120,'(f21.12)')
call wout('dchi/dr_i = '//trim(tmpr))
tmpr=r2s(xyzzyaaah120,'(f21.12)')
call wout('Gamma     = '//trim(tmpr))
call wout('So cusp condition on chi not satisfied.  The derivative sho&
&uld equal Gamma.')
call errstop('TEST_CUSPS_AND_CUTOFFS','Stopping.')
endif
enddo
if(xyzzyaagx1>=1)then
ri=xyzzyaacu1(xyzzyaaac120)
do xyzzyaaaa120=1,nspin
call xyzzyaaqa1(ri,xyzzyaaaa120,xyzzyaaac120,.true.,.true.,.true.,xyzz&
&yaaai120,xyzzyaaaj120,xyzzyaaak120)
if(abs(xyzzyaaai120)>xyzzyaaau120.or.(abs(xyzzyaaaj120)>xyzzyaaau120.a&
&nd.xyzzyaagx1>=2).or.(abs(xyzzyaaak120)>xyzzyaaau120.and.xyzzyaagx1>=&
&3))then
call wout('Chi set                   = '//trim(i2s(xyzzyaaac120)))
call wout('ispin                     = '//trim(i2s(xyzzyaaaa120)))
tmpr=r2s(xyzzyaaai120,'(f21.12)')
call wout('When r_i=L_chi, chi       = '//trim(tmpr))
tmpr=r2s(xyzzyaaaj120,'(f21.12)')
call wout('              dchi/dr_i   = '//trim(tmpr))
tmpr=r2s(xyzzyaaak120,'(f21.12)')
call wout('             d2chi/dr_i^2 = '//trim(tmpr))
call wout('So smooth cutoff conditions on chi are not satisfied.')
call errstop('TEST_CUSPS_AND_CUTOFFS','Stopping.')
endif
enddo
endif
enddo
endif
if(xyzzyaaex1)then
do xyzzyaaac120=1,xyzzyaabj1
ri=xyzzyaacv1(xyzzyaaac120)*0.5d0
rj=ri
xyzzyaaad120=0.d0
do xyzzyaaaa120=1,nspin
do xyzzyaaab120=xyzzyaaaa120,nspin
call xyzzyaaqb1(0,ri,rj,xyzzyaaad120,xyzzyaaaa120,xyzzyaaab120,xyzzyaa&
&ac120,.false.,.true.,.false.,.false.,xyzzyaaal120,xyzzyaaam120,xyzzya&
&aar120,xyzzyaaan120,xyzzyaaao120,xyzzyaaas120,xyzzyaaap120,xyzzyaaaq1&
&20,xyzzyaaat120)
if(abs(xyzzyaaan120)>xyzzyaaau120)then
call wout('F set    = '//trim(i2s(xyzzyaaac120)))
call wout('ispin    = '//trim(i2s(xyzzyaaaa120)))
call wout('jspin    = '//trim(i2s(xyzzyaaab120)))
tmpr=r2s(ri,'(f21.12)')
call wout('ri=rj    = '//trim(tmpr))
call wout('rij      = 0')
tmpr=r2s(xyzzyaaan120,'(f21.12)')
call wout('df/dr_ij = '//trim(tmpr))
call wout('So no-e-e-cusp condition on f is not satisfied.')
call errstop('TEST_CUSPS_AND_CUTOFFS','Stopping.')
endif
enddo
enddo
ri=0.d0
rj=xyzzyaacv1(xyzzyaaac120)*0.5d0
xyzzyaaad120=rj
do xyzzyaaaa120=1,nspin
do xyzzyaaab120=1,nspin
call xyzzyaaqb1(0,ri,rj,xyzzyaaad120,xyzzyaaaa120,xyzzyaaab120,xyzzyaa&
&ac120,.false.,.true.,.false.,.false.,xyzzyaaal120,xyzzyaaam120,xyzzya&
&aar120,xyzzyaaan120,xyzzyaaao120,xyzzyaaas120,xyzzyaaap120,xyzzyaaaq1&
&20,xyzzyaaat120)
if(abs(xyzzyaaam120)>xyzzyaaau120)then
call wout('F set   = '//trim(i2s(xyzzyaaac120)))
call wout('ispin   = '//trim(i2s(xyzzyaaaa120)))
call wout('jspin   = '//trim(i2s(xyzzyaaab120)))
tmpr=r2s(rj,'(f21.12)')
call wout('rj=rij  = '//trim(tmpr))
call wout('ri      = 0')
tmpr=r2s(xyzzyaaam120,'(f21.12)')
call wout('df/dr_i = '//trim(tmpr))
call wout('So no-e-N-cusp condition on f is not satisfied.')
call errstop('TEST_CUSPS_AND_CUTOFFS','Stopping.')
endif
enddo
enddo
if(xyzzyaagx1>=1)then
ri=xyzzyaacv1(xyzzyaaac120)
rj=xyzzyaacv1(xyzzyaaac120)*0.5d0
xyzzyaaad120=xyzzyaacv1(xyzzyaaac120)*1.25d0
do xyzzyaaaa120=1,nspin
do xyzzyaaab120=1,nspin
call xyzzyaaqb1(0,ri,rj,xyzzyaaad120,xyzzyaaaa120,xyzzyaaab120,xyzzyaa&
&ac120,.true.,.true.,.true.,.false.,xyzzyaaal120,xyzzyaaam120,xyzzyaaa&
&r120,xyzzyaaan120,xyzzyaaao120,xyzzyaaas120,xyzzyaaap120,xyzzyaaaq120&
&,xyzzyaaat120)
if(abs(xyzzyaaal120)>xyzzyaaau120.or.((abs(xyzzyaaam120)>xyzzyaaau120.&
&or.abs(xyzzyaaaq120)>xyzzyaaau120).and.xyzzyaagx1>=2).or.(abs(xyzzyaa&
&ao120)>xyzzyaaau120.and.xyzzyaagx1>=3))then
call wout('F set                     = '//trim(i2s(xyzzyaaac120)))
call wout('ispin                     = '//trim(i2s(xyzzyaaaa120)))
call wout('jspin                     = '//trim(i2s(xyzzyaaab120)))
tmpr=r2s(xyzzyaaal120,'(f21.12)')
call wout('When r_i=L_f, f           = '//trim(tmpr))
tmpr=r2s(xyzzyaaam120,'(f21.12)')
call wout('             df/dr_i      = '//trim(tmpr))
tmpr=r2s(xyzzyaaaq120,'(f21.12)')
call wout('            d2f/dr_idr_ij = '//trim(tmpr))
tmpr=r2s(xyzzyaaao120,'(f21.12)')
call wout('            d2f/dr_i^2    = '//trim(tmpr))
call wout('So smooth cutoff conditions on f are not satisfied.')
call errstop('TEST_CUSPS_AND_CUTOFFS','Stopping.')
endif
call xyzzyaaqb1(0,rj,ri,xyzzyaaad120,xyzzyaaaa120,xyzzyaaab120,xyzzyaa&
&ac120,.true.,.true.,.true.,.false.,xyzzyaaal120,xyzzyaaam120,xyzzyaaa&
&r120,xyzzyaaan120,xyzzyaaao120,xyzzyaaas120,xyzzyaaap120,xyzzyaaaq120&
&,xyzzyaaat120)
if(abs(xyzzyaaal120)>xyzzyaaau120.or.(abs(xyzzyaaaq120)>xyzzyaaau120.a&
&nd.xyzzyaagx1>=2))then
call wout('F set                     = '//trim(i2s(xyzzyaaac120)))
call wout('ispin                     = '//trim(i2s(xyzzyaaaa120)))
call wout('jspin                     = '//trim(i2s(xyzzyaaab120)))
tmpr=r2s(xyzzyaaal120,'(f21.12)')
call wout('When r_j=L_f, f           = '//trim(tmpr))
tmpr=r2s(xyzzyaaaq120,'(f21.12)')
call wout('            d2f/dr_idr_ij = '//trim(tmpr))
call wout('So smooth cutoff conditions on f are not satisfied.')
call errstop('TEST_CUSPS_AND_CUTOFFS','Stopping.')
endif
enddo
enddo
endif
enddo
endif
endif
end subroutine xyzzyaaqz1
subroutine pjastrow_assess_check_kinetic(eevecs,eivecs,ie,jspin,fail,v&
&erbose)
use slaarnacc,only : ranx
implicit none
real(dp),intent(in) :: eevecs(4,netot,netot),eivecs(4,nitot,netot)
integer,intent(out) :: ie,jspin
logical,intent(in) :: verbose
logical,intent(out) :: fail
integer xyzzyaaaa121,xyzzyaaab121,xyzzyaaac121,xyzzyaaad121
real(dp) xyzzyaaae121
real(dp),parameter :: xyzzyaaaf121=1.d-2
logical xyzzyaaag121,xyzzyaaah121,xyzzyaaai121,xyzzyaaaj121,xyzzyaaak1&
&21
logical,allocatable :: xyzzyaaal121(:)
if(xyzzyaafd1.or.xyzzyaafg1)then
ie=1
jspin=ceiling(nspin*ranx())
fail=.false.
return
endif
allocate(xyzzyaaal121(netot),stat=xyzzyaaad121)
call check_alloc(xyzzyaaad121,'JASTROW_ASSESS_CHECK_KINETIC','0.')
xyzzyaaal121=.false.
fail=.false.
ie=0
jspin=1
do
if(jspin>nspin)then
do xyzzyaaaa121=1,netot
if(xyzzyaaal121(xyzzyaaaa121))then
ie=which_ie(xyzzyaaaa121)
jspin=which_spin(xyzzyaaaa121)
deallocate(xyzzyaaal121)
return
endif
enddo
fail=.true.
deallocate(xyzzyaaal121)
return
endif
ie=ie+1
if(ie>nele(jspin))then
ie=0
jspin=jspin+1
cycle
endif
xyzzyaaaa121=which_ii(ie,jspin)
xyzzyaaal121(xyzzyaaaa121)=.true.
xyzzyaaag121=.false.
if(xyzzyaaep1.and.xyzzyaagx1<3)then
do xyzzyaaab121=1,netot
if(abs(eevecs(4,xyzzyaaab121,xyzzyaaaa121)-xyzzyaach1)<=xyzzyaaaf121)t&
&hen
xyzzyaaal121(xyzzyaaaa121)=.false.
if(verbose)call wout('For (i,j)=('//trim(i2s(xyzzyaaaa121))//','//trim&
&(i2s(xyzzyaaab121))//'), |rij - L_u| =',abs(eevecs(4,xyzzyaaab121,xyz&
&zyaaaa121)-xyzzyaach1))
endif
enddo
endif
if(xyzzyaaes1.and.xyzzyaagx1<3)then
do xyzzyaaab121=1,netot
xyzzyaaae121=dot_product(eevecs(1:3,xyzzyaaab121,xyzzyaaaa121),xyzzyaa&
&cr1)
if(abs(abs(xyzzyaaae121)-xyzzyaack1)<=xyzzyaaaf121.or.abs(sqrt(max(eev&
&ecs(4,xyzzyaaab121,xyzzyaaaa121)**2-xyzzyaaae121**2,0.d0))-xyzzyaack1&
&)<=xyzzyaaaf121)then
xyzzyaaal121(xyzzyaaaa121)=.false.
if(verbose)then
call wout('For (i,j)=('//trim(i2s(xyzzyaaaa121))//','//trim(i2s(xyzzya&
&aab121))//'), |rhoij - L_ucylrho| =',abs(sqrt(max(eevecs(4,xyzzyaaab1&
&21,xyzzyaaaa121)**2-xyzzyaaae121**2,0.d0))-xyzzyaack1))
call wout('For (i,j)=('//trim(i2s(xyzzyaaaa121))//','//trim(i2s(xyzzya&
&aab121))//'), ||zij| - L_ucylz| =',abs(abs(xyzzyaaae121)-xyzzyaack1))
endif
endif
enddo
endif
if(xyzzyaaew1.or.xyzzyaaex1)then
do xyzzyaaac121=1,nitot
if(xyzzyaagx1<3)then
if(xyzzyaaew1)then
if(xyzzyaabo1(xyzzyaaac121)/=0)then
if(abs(eivecs(4,xyzzyaaac121,xyzzyaaaa121)-xyzzyaacu1(xyzzyaabo1(xyzzy&
&aaac121)))<=xyzzyaaaf121)then
xyzzyaaal121(xyzzyaaaa121)=.false.
if(verbose)call wout('For (i,I)=('//trim(i2s(xyzzyaaaa121))//','//trim&
&(i2s(xyzzyaaac121))//'), |riI - L_chi| =',abs(eivecs(4,xyzzyaaac121,x&
&yzzyaaaa121)-xyzzyaacu1(xyzzyaabo1(xyzzyaaac121))))
endif
endif
endif
if(xyzzyaaex1)then
if(xyzzyaabp1(xyzzyaaac121)/=0)then
if(abs(eivecs(4,xyzzyaaac121,xyzzyaaaa121)-xyzzyaacv1(xyzzyaabp1(xyzzy&
&aaac121)))<=xyzzyaaaf121)then
xyzzyaaal121(xyzzyaaaa121)=.false.
if(verbose)call wout('For (i,I)=('//trim(i2s(xyzzyaaaa121))//','//trim&
&(i2s(xyzzyaaac121))//'), |riI - L_f| =',abs(eivecs(4,xyzzyaaac121,xyz&
&zyaaaa121)-xyzzyaacv1(xyzzyaabp1(xyzzyaaac121))))
endif
endif
endif
endif
if(eivecs(4,xyzzyaaac121,xyzzyaaaa121)<=xyzzyaaaf121)then
xyzzyaaal121(xyzzyaaaa121)=.false.
if(verbose)call wout('For (i,I)=('//trim(i2s(xyzzyaaaa121))//','//trim&
&(i2s(xyzzyaaac121))//'), r_iI =',eivecs(4,xyzzyaaac121,xyzzyaaaa121))
endif
enddo
endif
if(.not.xyzzyaaal121(xyzzyaaaa121))then
if(verbose)call wout('Particle '//trim(i2s(xyzzyaaaa121))//' skipped.'&
&)
cycle
endif
xyzzyaaah121=.true.
xyzzyaaai121=.true.
xyzzyaaaj121=.true.
xyzzyaaak121=.true.
if(xyzzyaaep1)then
xyzzyaaah121=.false.
do xyzzyaaab121=1,netot
if(xyzzyaaab121==xyzzyaaaa121)cycle
if(eevecs(4,xyzzyaaab121,xyzzyaaaa121)<xyzzyaach1)then
xyzzyaaah121=.true.
exit
endif
enddo
endif
if(xyzzyaaes1)then
xyzzyaaai121=.false.
do xyzzyaaab121=1,netot
if(xyzzyaaab121==xyzzyaaaa121)cycle
xyzzyaaae121=dot_product(eevecs(1:3,xyzzyaaab121,xyzzyaaaa121),xyzzyaa&
&cr1)
if(abs(xyzzyaaae121)<xyzzyaack1.and.eevecs(4,xyzzyaaab121,xyzzyaaaa121&
&)**2-xyzzyaaae121**2<xyzzyaacj1)then
xyzzyaaai121=.true.
exit
endif
enddo
endif
if(xyzzyaaew1)then
xyzzyaaaj121=.false.
do xyzzyaaac121=1,nitot
if(xyzzyaabo1(xyzzyaaac121)==0)cycle
if(eivecs(4,xyzzyaaac121,xyzzyaaaa121)<xyzzyaacu1(xyzzyaabo1(xyzzyaaac&
&121)))then
xyzzyaaaj121=.true.
exit
endif
enddo
endif
if(xyzzyaaex1)then
xyzzyaaak121=.false.
do xyzzyaaac121=1,nitot
if(xyzzyaabp1(xyzzyaaac121)==0)cycle
if(eivecs(4,xyzzyaaac121,xyzzyaaaa121)<xyzzyaacv1(xyzzyaabp1(xyzzyaaac&
&121)))then
do xyzzyaaab121=1,netot
if(xyzzyaaab121==xyzzyaaaa121)cycle
if(eivecs(4,xyzzyaaac121,xyzzyaaab121)<xyzzyaacv1(xyzzyaabp1(xyzzyaaac&
&121)))then
xyzzyaaak121=.true.
exit
endif
enddo
endif
if(xyzzyaaak121)exit
enddo
endif
xyzzyaaag121=xyzzyaaah121.and.xyzzyaaai121.and.xyzzyaaaj121.and.xyzzya&
&aak121
if(xyzzyaaag121)then
deallocate(xyzzyaaal121)
return
endif
if(verbose)call wout('Particle '//trim(i2s(xyzzyaaaa121))//' temporari&
&ly skipped.')
enddo
deallocate(xyzzyaaal121)
end subroutine pjastrow_assess_check_kinetic
subroutine setup_pjastrow_plot(makeplot,ispin,jspin,r_j,dir_i,r_i)
implicit none
integer,intent(in),optional :: ispin,jspin
real(dp),intent(in),optional :: r_j(3),dir_i(3),r_i(3)
logical,intent(in) :: makeplot
xyzzyaaho1=makeplot
xyzzyaahj1=0
if(present(ispin))xyzzyaahj1=ispin
xyzzyaahk1=0
if(present(jspin))xyzzyaahk1=jspin
xyzzyaahl1=0.d0
if(present(r_j))xyzzyaahl1=r_j
xyzzyaahm1=0.d0
if(present(dir_i))xyzzyaahm1=dir_i
xyzzyaahn1=0.d0
if(present(r_i))xyzzyaahn1=r_i
end subroutine setup_pjastrow_plot
subroutine xyzzyaara1(plot_u,plot_ucyl,plot_w,plot_chi,plot_f,plot_p,p&
&lot_q,ispin,jspin)
implicit none
integer,intent(in) :: ispin,jspin
logical,intent(in) :: plot_u,plot_ucyl,plot_w,plot_chi,plot_f,plot_p,p&
&lot_q
integer xyzzyaaaa123,xyzzyaaab123,xyzzyaaac123,xyzzyaaad123,xyzzyaaae1&
&23,xyzzyaaaf123
real(dp) xyzzyaaag123,xyzzyaaah123,xyzzyaaai123,xyzzyaaaj123,xyzzyaaak&
&123,xyzzyaaal123,xyzzyaaam123,xyzzyaaan123,xyzzyaaao123,xyzzyaaap123,&
&xyzzyaaaq123,xyzzyaaar123,xyzzyaaas123,xyzzyaaat123,xyzzyaaau123,xyzz&
&yaaav123,xyzzyaaaw123,xyzzyaaax123,xyzzyaaay123,xyzzyaaaz123,xyzzyaab&
&a123,xyzzyaabb123,xyzzyaabc123,xyzzyaabd123(3),xyzzyaabe123(3),xyzzya&
&abf123,xyzzyaabg123(3),xyzzyaabh123,xyzzyaabi123,xyzzyaabj123(3),xyzz&
&yaabk123,xyzzyaabl123,xyzzyaabm123,xyzzyaabn123,xyzzyaabo123,xyzzyaab&
&p123,xyzzyaabq123,xyzzyaabr123,xyzzyaabs123,xyzzyaabt123,xyzzyaabu123&
&,xyzzyaabv123,xyzzyaabw123
logical xyzzyaabx123
character(80) tmpr,tmpr2,tmpr3
integer,parameter :: xyzzyaaby123=2000
if(am_master)then
if(plot_u.and..not.xyzzyaaep1)call errstop('PLOT_JASTROW','Trying to p&
&lot the u term when it''s not available.')
if(plot_ucyl.and..not.xyzzyaaes1)call errstop('PLOT_JASTROW','Trying t&
&o plot the ucyl term when it''s not available.')
if(plot_w.and..not.xyzzyaaeu1)call errstop('PLOT_JASTROW','Trying to p&
&lot the w term when it''s not available.')
if(plot_chi.and..not.xyzzyaaew1)call errstop('PLOT_JASTROW','Trying to&
& plot chi terms when they''re not available.')
if(plot_f.and..not.xyzzyaaex1)call errstop('PLOT_JASTROW','Trying to p&
&lot f terms when they''re not available.')
if(plot_p.and..not.xyzzyaaey1)call errstop('PLOT_JASTROW','Trying to p&
&lot p terms when they''re not available.')
if(plot_q.and..not.xyzzyaaez1)call errstop('PLOT_JASTROW','Trying to p&
&lot q terms when they''re not available.')
if(plot_ucyl.or.plot_f.or.plot_p.or.plot_q)then
if(all(xyzzyaahm1==0.d0))call errstop('PLOT_JASTROW','The direction ve&
&ctor for electron i is zero.  Please correct your jastrow_plot block &
&in input.')
xyzzyaahm1(1:3)=xyzzyaahm1(1:3)/sqrt(dot_product(xyzzyaahm1,xyzzyaahm1&
&))
endif
call open_units(xyzzyaaad123,xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Unable to find free i/&
&o unit [1].')
call open_units(xyzzyaaae123,xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Unable to find free i/&
&o unit [2].')
call open_units(xyzzyaaaf123,xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Unable to find free i/&
&o unit [3].')
if(plot_u)then
open(unit=xyzzyaaad123,file='jastrow_value_u.dat',status='replace',ios&
&tat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_value_u.dat.')
open(unit=xyzzyaaae123,file='jastrow_deriv_u.dat',status='replace',ios&
&tat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_deriv_u.dat.')
open(unit=xyzzyaaaf123,file='jastrow_sderiv_u.dat',status='replace',io&
&stat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_sderiv_u.dat.')
if(hard_sphere.and.(.not.hard_op_spins.or.ispin/=jspin))then
xyzzyaabo123=(max(xyzzyaach1,xyzzyaanz1)-hard_diam)/dble(xyzzyaaby123)
else
xyzzyaabo123=xyzzyaach1/dble(xyzzyaaby123)
endif
do xyzzyaaaa123=0,xyzzyaaby123
if(hard_sphere.and.(.not.hard_op_spins.or.ispin/=jspin))then
if(xyzzyaaaa123==0)cycle
xyzzyaaag123=hard_diam+dble(xyzzyaaaa123)*xyzzyaabo123
else
xyzzyaaag123=dble(xyzzyaaaa123)*xyzzyaabo123
endif
call xyzzyaapr1(xyzzyaaag123,ispin,jspin,.true.,.true.,.true.,xyzzyaaa&
&h123,xyzzyaaai123,xyzzyaaaj123)
write(xyzzyaaad123,*)xyzzyaaag123,xyzzyaaah123
write(xyzzyaaae123,*)xyzzyaaag123,xyzzyaaai123
write(xyzzyaaaf123,*)xyzzyaaag123,xyzzyaaaj123
enddo
close(xyzzyaaad123)
open_unit(xyzzyaaad123)=.false.
close(xyzzyaaae123)
open_unit(xyzzyaaae123)=.false.
close(xyzzyaaaf123)
open_unit(xyzzyaaaf123)=.false.
call wout('Have made a plot of u(rij) and its first and second derivat&
&ives in the files:')
call wout('jastrow_value_u.dat, jastrow_deriv_u.dat and jastrow_sderiv&
&_u.dat.')
call wout('Spin of electron i : '//trim(i2s(ispin)))
call wout('Spin of electron j : '//trim(i2s(jspin)))
call wout()
endif
if(plot_ucyl)then
open(unit=xyzzyaaad123,file='jastrow_value_ucyl.dat',status='replace',&
&iostat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_value_ucyl.dat.')
open(unit=xyzzyaaaf123,file='jastrow_lap_ucyl.dat',status='replace',io&
&stat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_lap_ucyl.dat.')
xyzzyaaba123=-sqrt(xyzzyaacj1+xyzzyaack1**2)
xyzzyaabb123=-xyzzyaaba123
do xyzzyaaaa123=0,xyzzyaaby123
xyzzyaabc123=xyzzyaaba123+dble(xyzzyaaaa123)*(xyzzyaabb123-xyzzyaaba12&
&3)/dble(xyzzyaaby123)
xyzzyaabd123(1:3)=xyzzyaahn1(1:3)+xyzzyaabc123*xyzzyaahm1(1:3)
xyzzyaabe123=xyzzyaabd123-xyzzyaahl1
xyzzyaabv123=dot_product(xyzzyaabe123,xyzzyaacr1)
xyzzyaabu123=dot_product(xyzzyaabe123,xyzzyaabe123)-xyzzyaabv123*xyzzy&
&aabv123
if(xyzzyaabu123<xyzzyaacj1.or.xyzzyaagx1==0)then
xyzzyaabu123=sqrt(max(xyzzyaabu123,0.d0))
call xyzzyaapu1(xyzzyaabu123,abs(xyzzyaabv123),ispin,jspin,.true.,.tru&
&e.,.true.,xyzzyaabp123,xyzzyaabq123,xyzzyaabr123,xyzzyaabs123,xyzzyaa&
&bt123)
if(xyzzyaabu123>0.d0)then
xyzzyaabw123=xyzzyaabq123/xyzzyaabu123
else
xyzzyaabw123=0.d0
endif
else
xyzzyaabp123=0.d0
xyzzyaabq123=0.d0
xyzzyaabr123=0.d0
xyzzyaabs123=0.d0
xyzzyaabt123=0.d0
xyzzyaabw123=0.d0
endif
write(xyzzyaaad123,*)xyzzyaabc123,xyzzyaabp123
write(xyzzyaaaf123,*)xyzzyaabc123,xyzzyaabs123+xyzzyaabw123+xyzzyaabt1&
&23
enddo
close(xyzzyaaad123)
close(xyzzyaaaf123)
call wordwrap('Have made a plot of ucyl(r_ij) in the file jastrow_valu&
&e_ucyl.dat and lap[ucyl(r_ij)] in the file jastrow_lap_ucyl.dat.')
call wout('Spin of electron i : '//trim(i2s(ispin)))
call wout('Spin of electron j : '//trim(i2s(jspin)))
tmpr=r2s(xyzzyaahl1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahl1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahl1(3),'(f21.12)')
call wout('Position of j      : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
tmpr=r2s(xyzzyaahm1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahm1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahm1(3),'(f21.12)')
call wout('Direction of i     : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
tmpr=r2s(xyzzyaahn1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahn1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahn1(3),'(f21.12)')
call wout('Point on path of i : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
call wout('(NB, distance along line measured from this point.)')
call wout()
endif
if(plot_w)then
open(unit=xyzzyaaad123,file='jastrow_value_w.dat',status='replace',ios&
&tat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_value_u.dat.')
open(unit=xyzzyaaae123,file='jastrow_deriv_w.dat',status='replace',ios&
&tat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_deriv_u.dat.')
open(unit=xyzzyaaaf123,file='jastrow_sderiv_w.dat',status='replace',io&
&stat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_sderiv_u.dat.')
do xyzzyaaaa123=0,xyzzyaaby123
xyzzyaaag123=dble(xyzzyaaaa123)*xyzzyaacl1/dble(xyzzyaaby123)
call xyzzyaapx1(xyzzyaaag123,ispin,jspin,.true.,.true.,.true.,xyzzyaab&
&l123,xyzzyaabm123,xyzzyaabn123,xyzzyaabx123)
write(xyzzyaaad123,*)xyzzyaaag123,xyzzyaabl123
write(xyzzyaaae123,*)xyzzyaaag123,xyzzyaabm123
write(xyzzyaaaf123,*)xyzzyaaag123,xyzzyaabn123
enddo
close(xyzzyaaad123)
open_unit(xyzzyaaad123)=.false.
close(xyzzyaaae123)
open_unit(xyzzyaaae123)=.false.
close(xyzzyaaaf123)
open_unit(xyzzyaaaf123)=.false.
call wout('Have made a plot of w(rij) and its first and second derivat&
&ives in the files:')
call wout('jastrow_value_w.dat, jastrow_deriv_w.dat and jastrow_sderiv&
&_w.dat.')
call wout('Spin of electron i : '//trim(i2s(ispin)))
call wout('Spin of electron j : '//trim(i2s(jspin)))
call wout()
endif
if(plot_chi)then
do xyzzyaaac123=1,xyzzyaabi1
open(unit=xyzzyaaad123,file='jastrow_value_chi_'//trim(i2s(xyzzyaaac12&
&3))//'.dat',status='replace',iostat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_value_chi_'//trim(i2s(xyzzyaaac123))//'.dat.')
open(unit=xyzzyaaae123,file='jastrow_deriv_chi_'//trim(i2s(xyzzyaaac12&
&3))//'.dat',status='replace',iostat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_deriv_chi_'//trim(i2s(xyzzyaaac123))//'.dat.')
open(unit=xyzzyaaaf123,file='jastrow_sderiv_chi_'//trim(i2s(xyzzyaaac1&
&23))//'.dat',status='replace',iostat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_sderiv_chi_'//trim(i2s(xyzzyaaac123))//'.dat.')
do xyzzyaaaa123=0,xyzzyaaby123
xyzzyaaak123=dble(xyzzyaaaa123)*xyzzyaacu1(xyzzyaaac123)/dble(xyzzyaab&
&y123)
call xyzzyaaqa1(xyzzyaaak123,ispin,xyzzyaaac123,.true.,.true.,.true.,x&
&yzzyaaal123,xyzzyaaam123,xyzzyaaan123)
write(xyzzyaaad123,*)xyzzyaaak123,xyzzyaaal123
write(xyzzyaaae123,*)xyzzyaaak123,xyzzyaaam123
write(xyzzyaaaf123,*)xyzzyaaak123,xyzzyaaan123
enddo
close(xyzzyaaad123)
close(xyzzyaaae123)
close(xyzzyaaaf123)
call wout('Have made a plot of chi(ri) and its first and second deriva&
&tives in the files:')
call wout('jastrow_value_chi_'//trim(i2s(xyzzyaaac123))//'.dat, jastro&
&w_deriv_chi_'//trim(i2s(xyzzyaaac123))//'.dat and jastrow_sderiv_chi_&
&'//trim(i2s(xyzzyaaac123))//'.dat.')
call wout('Spin of electron i : '//trim(i2s(ispin)))
call wout('Set of chi terms   : '//trim(i2s(xyzzyaaac123)))
call wout()
enddo
endif
if(plot_f)then
xyzzyaaap123=sqrt(dot_product(xyzzyaahl1,xyzzyaahl1))
do xyzzyaaac123=1,xyzzyaabj1
open(unit=xyzzyaaad123,file='jastrow_value_f_'//trim(i2s(xyzzyaaac123)&
&)//'.dat',status='replace',iostat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_value_f_'//trim(i2s(xyzzyaaac123))//'.dat.')
xyzzyaaba123=-xyzzyaacv1(xyzzyaaac123)-abs(dot_product(xyzzyaahn1,xyzz&
&yaahm1))
xyzzyaabb123=-xyzzyaaba123
do xyzzyaaaa123=0,xyzzyaaby123
xyzzyaabc123=xyzzyaaba123+dble(xyzzyaaaa123)*(xyzzyaabb123-xyzzyaaba12&
&3)/dble(xyzzyaaby123)
xyzzyaabd123(1:3)=xyzzyaahn1(1:3)+xyzzyaabc123*xyzzyaahm1(1:3)
xyzzyaaaq123=sqrt(dot_product(xyzzyaabd123-xyzzyaahl1,xyzzyaabd123-xyz&
&zyaahl1))
xyzzyaaao123=sqrt(dot_product(xyzzyaabd123,xyzzyaabd123))
call xyzzyaaqb1(0,xyzzyaaao123,xyzzyaaap123,xyzzyaaaq123,ispin,jspin,x&
&yzzyaaac123,.true.,.false.,.false.,.false.,xyzzyaaar123,xyzzyaaas123,&
&xyzzyaaax123,xyzzyaaat123,xyzzyaaau123,xyzzyaaay123,xyzzyaaav123,xyzz&
&yaaaw123,xyzzyaaaz123)
write(xyzzyaaad123,*)xyzzyaabc123,xyzzyaaar123
enddo
close(xyzzyaaad123)
call wout('Have made a plot of f(ri,rj,rij) in the file jastrow_value_&
&f_'//trim(i2s(xyzzyaaac123))//'.dat.')
call wout('Spin of electron i : '//trim(i2s(ispin)))
call wout('Spin of electron j : '//trim(i2s(jspin)))
call wout('Set of f terms     : '//trim(i2s(xyzzyaaac123)))
tmpr=r2s(xyzzyaahl1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahl1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahl1(3),'(f21.12)')
call wout('Position of j      : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
tmpr=r2s(xyzzyaahm1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahm1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahm1(3),'(f21.12)')
call wout('Direction of i     : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
tmpr=r2s(xyzzyaahn1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahn1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahn1(3),'(f21.12)')
call wout('Point on path of i : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
call wout('(NB, distance along line measured from this point.)')
call wout()
enddo
endif
if(plot_p)then
open(unit=xyzzyaaad123,file='jastrow_value_p.dat',status='replace',ios&
&tat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_value_p.dat.')
open(unit=xyzzyaaaf123,file='jastrow_lap_p.dat',status='replace',iosta&
&t=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_lap_p.dat.')
xyzzyaaba123=-2.d0*wigner_seitz_radius
xyzzyaabb123=-xyzzyaaba123
do xyzzyaaaa123=0,xyzzyaaby123
xyzzyaabc123=xyzzyaaba123+dble(xyzzyaaaa123)*(xyzzyaabb123-xyzzyaaba12&
&3)/dble(xyzzyaaby123)
xyzzyaabd123(1:3)=xyzzyaahn1(1:3)+xyzzyaabc123*xyzzyaahm1(1:3)
xyzzyaabe123=xyzzyaabd123-xyzzyaahl1
call xyzzyaaqd1(xyzzyaabe123,ispin,jspin,.true.,.false.,.true.,xyzzyaa&
&bf123,xyzzyaabg123,xyzzyaabh123)
write(xyzzyaaad123,*)xyzzyaabc123,xyzzyaabf123
write(xyzzyaaaf123,*)xyzzyaabc123,xyzzyaabh123
enddo
close(xyzzyaaad123)
close(xyzzyaaaf123)
call wordwrap('Have made a plot of p(r_ij) in the file jastrow_value_p&
&.dat and lap[p(r_ij)] in the file jastrow_lap_p.dat.')
call wout('Spin of electron i : '//trim(i2s(ispin)))
call wout('Spin of electron j : '//trim(i2s(jspin)))
tmpr=r2s(xyzzyaahl1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahl1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahl1(3),'(f21.12)')
call wout('Position of j      : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
tmpr=r2s(xyzzyaahm1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahm1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahm1(3),'(f21.12)')
call wout('Direction of i     : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
tmpr=r2s(xyzzyaahn1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahn1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahn1(3),'(f21.12)')
call wout('Point on path of i : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
call wout('(NB, distance along line measured from this point.)')
call wout()
endif
if(plot_q)then
open(unit=xyzzyaaad123,file='jastrow_value_q.dat',status='replace',ios&
&tat=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_value_q.dat.')
open(unit=xyzzyaaaf123,file='jastrow_lap_q.dat',status='replace',iosta&
&t=xyzzyaaab123)
if(xyzzyaaab123/=0)call errstop('PLOT_JASTROW','Problem opening jastro&
&w_lap_q.dat.')
xyzzyaaba123=-2.d0*wigner_seitz_radius
xyzzyaabb123=-xyzzyaaba123
do xyzzyaaaa123=0,xyzzyaaby123
xyzzyaabc123=xyzzyaaba123+dble(xyzzyaaaa123)*(xyzzyaabb123-xyzzyaaba12&
&3)/dble(xyzzyaaby123)
xyzzyaabd123(1:3)=xyzzyaahn1(1:3)+xyzzyaabc123*xyzzyaahm1(1:3)
call xyzzyaaqe1(xyzzyaabd123,ispin,.true.,.false.,.true.,xyzzyaabi123,&
&xyzzyaabj123,xyzzyaabk123)
write(xyzzyaaad123,*)xyzzyaabc123,xyzzyaabi123
write(xyzzyaaaf123,*)xyzzyaabc123,xyzzyaabk123
enddo
close(xyzzyaaad123)
close(xyzzyaaaf123)
call wordwrap('Have made a plot of q(r) in the file jastrow_value_q.da&
&t and lap[q(r)] in the file jastrow_lap_q.dat.')
call wout('Spin of electron i : '//trim(i2s(ispin)))
tmpr=r2s(xyzzyaahm1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahm1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahm1(3),'(f21.12)')
call wout('Direction of i     : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
tmpr=r2s(xyzzyaahn1(1),'(f21.12)')
tmpr2=r2s(xyzzyaahn1(2),'(f21.12)')
tmpr3=r2s(xyzzyaahn1(3),'(f21.12)')
call wout('Point on path of i : ('//trim(tmpr)//','//trim(tmpr2)//','/&
&/trim(tmpr3)//')')
call wout('(NB, distance along line measured from this point.)')
call wout()
endif
call wout('Have finished plotting the Jastrow factor.')
call wout()
endif
end subroutine xyzzyaara1
subroutine get_linear_basis_pjastrow(is,no_lin_param,nparam_all,ignore&
&,lbasis_grad_f,lbasis_lap_f,grad_j0,lap_j0)
use slaarnabt, only : dcopy
implicit none
integer,intent(in) :: is,no_lin_param,nparam_all
real(dp),intent(out) :: lbasis_grad_f(three_netot,no_lin_param),lbasis&
&_lap_f(no_lin_param),grad_j0(three_netot),lap_j0
logical,intent(in) :: ignore(nparam_all)
integer xyzzyaaaa124,xyzzyaaab124,xyzzyaaac124,xyzzyaaad124,xyzzyaaae1&
&24,xyzzyaaaf124,xyzzyaaag124,xyzzyaaah124,xyzzyaaai124,xyzzyaaaj124,x&
&yzzyaaak124,xyzzyaaal124,xyzzyaaam124
real(dp) xyzzyaaan124,xyzzyaaao124(3),xyzzyaaap124,xyzzyaaaq124,xyzzya&
&aar124,xyzzyaaas124(3),xyzzyaaat124,xyzzyaaau124,xyzzyaaav124,xyzzyaa&
&aw124(3),xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124(3),xyzzyaaba124,xyzzy&
&aabb124,xyzzyaabc124(3),xyzzyaabd124,xyzzyaabe124,xyzzyaabf124,xyzzya&
&abg124(3),xyzzyaabh124,xyzzyaabi124,xyzzyaabj124,xyzzyaabk124(3),xyzz&
&yaabl124,xyzzyaabm124,xyzzyaabn124,xyzzyaabo124(3),xyzzyaabp124,xyzzy&
&aabq124,xyzzyaabr124,xyzzyaabs124(3),xyzzyaabt124,xyzzyaabu124,xyzzya&
&abv124(3),xyzzyaabw124,xyzzyaabx124,xyzzyaaby124,xyzzyaabz124
integer xyzzyaaca124,xyzzyaacb124(3),xyzzyaacc124(5)
real(dp),allocatable :: xyzzyaacd124(:,:),xyzzyaace124(:,:,:),xyzzyaac&
&f124(:,:,:),xyzzyaacg124(:,:,:,:,:),xyzzyaach124(:,:),xyzzyaaci124(:,&
&:),xyzzyaacj124(:,:,:,:),xyzzyaack124(:),xyzzyaacl124(:),xyzzyaacm124&
&(:),xyzzyaacn124(:),xyzzyaaco124(:),xyzzyaacp124(:),xyzzyaacq124(:),x&
&yzzyaacr124(:)
integer xyzzyaacs124,xyzzyaact124,xyzzyaacu124,xyzzyaacv124,xyzzyaacw1&
&24,xyzzyaacx124,xyzzyaacy124
call get_rsele(is)
call get_eevecs(is)
call get_eivecs(is)
call timer('JASTROW',.true.)
if(xyzzyaaeq1)then
allocate(xyzzyaacd124(0:xyzzyaabx1,xyzzyaafp1),stat=xyzzyaaca124)
call check_alloc(xyzzyaaca124,'INIT_LINEAR_BASIS','1')
xyzzyaacs124=(xyzzyaabx1+1)*xyzzyaafp1
call dcopy(xyzzyaacs124,xyzzyaaab1(0,1),1,xyzzyaacd124(0,1),1)
endif
if(xyzzyaaes1)then
allocate(xyzzyaace124(0:xyzzyaaby1,0:xyzzyaabz1,xyzzyaafq1),stat=xyzzy&
&aaca124)
call check_alloc(xyzzyaaca124,'INIT_LINEAR_BASIS','1.5')
xyzzyaact124=(xyzzyaaby1+1)*(xyzzyaabz1+1)*xyzzyaafq1
call dcopy(xyzzyaact124,xyzzyaaac1(0,0,1),1,xyzzyaace124(0,0,1),1)
endif
if(xyzzyaaew1)then
xyzzyaacb124=shape(xyzzyaaad1)
allocate(xyzzyaacf124(0:xyzzyaacb124(1)-1,xyzzyaacb124(2),xyzzyaacb124&
&(3)),stat=xyzzyaaca124)
call check_alloc(xyzzyaaca124,'INIT_LINEAR_BASIS','2')
xyzzyaacu124=xyzzyaacb124(1)*xyzzyaacb124(2)*xyzzyaacb124(3)
call dcopy(xyzzyaacu124,xyzzyaaad1(0,1,1),1,xyzzyaacf124(0,1,1),1)
endif
if(xyzzyaaex1)then
xyzzyaacc124=shape(xyzzyaaae1)
allocate(xyzzyaacg124(0:xyzzyaacc124(1)-1,0:xyzzyaacc124(2)-1,0:xyzzya&
&acc124(3)-1,xyzzyaacc124(4),xyzzyaacc124(5)),stat=xyzzyaaca124)
call check_alloc(xyzzyaaca124,'INIT_LINEAR_BASIS','3')
xyzzyaacv124=xyzzyaacc124(1)*xyzzyaacc124(2)*xyzzyaacc124(3)*xyzzyaacc&
&124(4)*xyzzyaacc124(5)
call dcopy(xyzzyaacv124,xyzzyaaae1(0,0,0,1,1),1,xyzzyaacg124(0,0,0,1,1&
&),1)
endif
if(xyzzyaaey1)then
allocate(xyzzyaach124(xyzzyaagq1,xyzzyaafr1),stat=xyzzyaaca124)
call check_alloc(xyzzyaaca124,'INIT_LINEAR_BASIS','4')
xyzzyaacw124=xyzzyaagq1*xyzzyaafr1
call dcopy(xyzzyaacw124,xyzzyaaaf1(1,1),1,xyzzyaach124(1,1),1)
endif
if(xyzzyaaez1)then
allocate(xyzzyaaci124(xyzzyaagr1,xyzzyaafs1),stat=xyzzyaaca124)
call check_alloc(xyzzyaaca124,'INIT_LINEAR_BASIS','5')
xyzzyaacx124=xyzzyaagr1*xyzzyaafs1
call dcopy(xyzzyaacx124,xyzzyaaag1(1,1),1,xyzzyaaci124(1,1),1)
endif
if(xyzzyaaff1)then
allocate(xyzzyaacj124(xyzzyaabn1,2,xyzzyaacd1,xyzzyaacd1),stat=xyzzyaa&
&ca124)
call check_alloc(xyzzyaaca124,'INIT_LINEAR_BASIS','6')
xyzzyaacy124=xyzzyaabn1*2*xyzzyaacd1*xyzzyaacd1
call dcopy(xyzzyaacy124,xyzzyaaas1(1,1,1,1),1,xyzzyaacj124(1,1,1,1),1)
endif
allocate(xyzzyaacl124(three_netot),xyzzyaacm124(three_netot),xyzzyaacn&
&124(three_netot),xyzzyaaco124(three_netot),xyzzyaacp124(three_netot),&
&xyzzyaacq124(three_netot),xyzzyaacr124(three_netot),stat=xyzzyaaca124&
&)
call check_alloc(xyzzyaaca124,'INIT_LINEAR_BASIS','7')
allocate(xyzzyaack124(nspin),stat=xyzzyaaca124)
call check_alloc(xyzzyaaca124,'GET_LINEAR_BASIS','8')
xyzzyaack124=sqrt(inv_pmass)
if(xyzzyaaeq1)then
do xyzzyaaad124=1,xyzzyaafp1
do xyzzyaaai124=0,xyzzyaabx1
if(xyzzyaaat1(xyzzyaaai124,xyzzyaaad124)==1)xyzzyaaab1(xyzzyaaai124,xy&
&zzyaaad124)=0.d0
enddo
enddo
call xyzzyaapd1
endif
if(xyzzyaaes1)then
do xyzzyaaad124=1,xyzzyaafq1
do xyzzyaaaj124=0,xyzzyaabz1
do xyzzyaaai124=0,xyzzyaaby1
if(xyzzyaaau1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaad124)==1)xyzzyaaac1(xy&
&zzyaaai124,xyzzyaaaj124,xyzzyaaad124)=0.d0
enddo
enddo
enddo
call xyzzyaape1
endif
if(xyzzyaaew1)then
do xyzzyaaae124=1,xyzzyaabi1
do xyzzyaaad124=1,xyzzyaafv1(xyzzyaaae124)
do xyzzyaaaj124=0,xyzzyaace1(xyzzyaaae124)
if(xyzzyaaav1(xyzzyaaaj124,xyzzyaaad124,xyzzyaaae124)==1)xyzzyaaad1(xy&
&zzyaaaj124,xyzzyaaad124,xyzzyaaae124)=0.d0
enddo
enddo
enddo
call xyzzyaapg1
endif
if(xyzzyaaex1)then
do xyzzyaaae124=1,xyzzyaabj1
do xyzzyaaad124=1,xyzzyaafw1(xyzzyaaae124)
do xyzzyaaak124=0,xyzzyaacf1(xyzzyaaae124)
do xyzzyaaaj124=0,xyzzyaacg1(xyzzyaaae124)
do xyzzyaaai124=0,xyzzyaacg1(xyzzyaaae124)
if(xyzzyaaaw1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaak124,xyzzyaaad124,xyzz&
&yaaae124)==1)xyzzyaaae1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaak124,xyzzya&
&aad124,xyzzyaaae124)=0.d0
enddo
enddo
enddo
enddo
enddo
call xyzzyaaon1
call xyzzyaaph1
endif
if(xyzzyaaey1)then
do xyzzyaaad124=1,xyzzyaafr1
do xyzzyaaab124=1,xyzzyaagq1
if(xyzzyaaax1(xyzzyaaab124,xyzzyaaad124)==1)xyzzyaaaf1(xyzzyaaab124,xy&
&zzyaaad124)=0.d0
enddo
enddo
call xyzzyaapl1
endif
if(xyzzyaaez1)then
do xyzzyaaad124=1,xyzzyaafs1
do xyzzyaaab124=1,xyzzyaagr1
if(xyzzyaaay1(xyzzyaaab124,xyzzyaaad124)==1)xyzzyaaag1(xyzzyaaab124,xy&
&zzyaaad124)=0.d0
enddo
enddo
call xyzzyaapm1
endif
if(xyzzyaaff1)then
where(xyzzyaabe1(:,:,:,:)==1)xyzzyaaas1(:,:,:,:)=0.d0
endif
grad_j0=0.d0
lap_j0=0.d0
xyzzyaacl124=0.d0
xyzzyaacm124=0.d0
xyzzyaacn124=0.d0
xyzzyaaco124=0.d0
xyzzyaacp124=0.d0
xyzzyaacq124=0.d0
xyzzyaacr124=0.d0
xyzzyaaaq124=0.d0
xyzzyaaau124=0.d0
xyzzyaabe124=0.d0
xyzzyaabi124=0.d0
xyzzyaabm124=0.d0
xyzzyaabq124=0.d0
xyzzyaabx124=0.d0
if(have_jastrow3)call xyzzyaapo1(eevecs_scr(1,1,1,is),.true.,.true.)
xyzzyaaaa124=0
xyzzyaaac124=0
do xyzzyaaaf124=1,nspin
xyzzyaabz124=inv_pmass(xyzzyaaaf124)
xyzzyaaby124=xyzzyaack124(xyzzyaaaf124)
do xyzzyaaag124=1,nele(xyzzyaaaf124)
xyzzyaaaa124=xyzzyaaaa124+1
call xyzzyaapq1(xyzzyaaaa124,xyzzyaaaf124,rele_scr(1,xyzzyaaaa124,is),&
&eevecs_scr(1,1,1,is),eevecs_scr(1,1,xyzzyaaaa124,is),eivecs_scr(1,1,1&
&,is),.false.,.true.,.true.,.true.,xyzzyaaan124,xyzzyaaao124,xyzzyaaap&
&124,xyzzyaaar124,xyzzyaaas124,xyzzyaaat124,xyzzyaaav124,xyzzyaaaw124,&
&xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124,xyzzyaaba124,xyzzyaabb124,xyzz&
&yaabc124,xyzzyaabd124,xyzzyaabf124,xyzzyaabg124,xyzzyaabh124,xyzzyaab&
&j124,xyzzyaabk124,xyzzyaabl124,xyzzyaabn124,xyzzyaabo124,xyzzyaabp124&
&,xyzzyaabr124,xyzzyaabs124,xyzzyaabt124,xyzzyaabu124,xyzzyaabv124,xyz&
&zyaabw124,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.tr&
&ue.,.true.)
do xyzzyaaab124=1,3
xyzzyaaac124=xyzzyaaac124+1
grad_j0(xyzzyaaac124)=(xyzzyaaao124(xyzzyaaab124)+xyzzyaaas124(xyzzyaa&
&ab124)+xyzzyaabc124(xyzzyaaab124)+xyzzyaabg124(xyzzyaaab124)+xyzzyaab&
&k124(xyzzyaaab124)+xyzzyaabo124(xyzzyaaab124)+xyzzyaabs124(xyzzyaaab1&
&24)+xyzzyaabv124(xyzzyaaab124))*xyzzyaaby124
xyzzyaacl124(xyzzyaaac124)=xyzzyaaao124(xyzzyaaab124)
xyzzyaacm124(xyzzyaaac124)=xyzzyaaas124(xyzzyaaab124)
xyzzyaacn124(xyzzyaaac124)=xyzzyaabc124(xyzzyaaab124)
xyzzyaaco124(xyzzyaaac124)=xyzzyaabg124(xyzzyaaab124)
xyzzyaacp124(xyzzyaaac124)=xyzzyaabk124(xyzzyaaab124)
xyzzyaacq124(xyzzyaaac124)=xyzzyaabo124(xyzzyaaab124)
xyzzyaacr124(xyzzyaaac124)=xyzzyaabv124(xyzzyaaab124)
enddo
lap_j0=lap_j0+(xyzzyaaap124+xyzzyaaat124+xyzzyaabd124+xyzzyaabh124+xyz&
&zyaabl124+xyzzyaabp124+xyzzyaabt124+xyzzyaabw124)*xyzzyaabz124
xyzzyaaaq124=xyzzyaaaq124+xyzzyaaap124*xyzzyaabz124
xyzzyaaau124=xyzzyaaau124+xyzzyaaat124*xyzzyaabz124
xyzzyaabe124=xyzzyaabe124+xyzzyaabd124*xyzzyaabz124
xyzzyaabi124=xyzzyaabi124+xyzzyaabh124*xyzzyaabz124
xyzzyaabm124=xyzzyaabm124+xyzzyaabl124*xyzzyaabz124
xyzzyaabq124=xyzzyaabq124+xyzzyaabp124*xyzzyaabz124
xyzzyaabx124=xyzzyaabx124+xyzzyaabw124*xyzzyaabz124
enddo
enddo
xyzzyaaah124=0
lbasis_grad_f=0.d0
lbasis_lap_f=0.d0
if(xyzzyaaeq1)then
do xyzzyaaad124=1,xyzzyaafp1
do xyzzyaaai124=0,xyzzyaabx1
if(xyzzyaaat1(xyzzyaaai124,xyzzyaaad124)==1)then
xyzzyaaah124=xyzzyaaah124+1
xyzzyaaab1(xyzzyaaai124,xyzzyaaad124)=1.d0
call xyzzyaapd1
xyzzyaaaa124=0
xyzzyaaac124=0
do xyzzyaaaf124=1,nspin
xyzzyaabz124=inv_pmass(xyzzyaaaf124)
xyzzyaaby124=xyzzyaack124(xyzzyaaaf124)
do xyzzyaaag124=1,nele(xyzzyaaaf124)
xyzzyaaaa124=xyzzyaaaa124+1
call xyzzyaapq1(xyzzyaaaa124,xyzzyaaaf124,rele_scr(1,xyzzyaaaa124,is),&
&eevecs_scr(1,1,1,is),eevecs_scr(1,1,xyzzyaaaa124,is),eivecs_scr(1,1,1&
&,is),.false.,.true.,.true.,.true.,xyzzyaaan124,xyzzyaaao124,xyzzyaaap&
&124,xyzzyaaar124,xyzzyaaas124,xyzzyaaat124,xyzzyaaav124,xyzzyaaaw124,&
&xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124,xyzzyaaba124,xyzzyaabb124,xyzz&
&yaabc124,xyzzyaabd124,xyzzyaabf124,xyzzyaabg124,xyzzyaabh124,xyzzyaab&
&j124,xyzzyaabk124,xyzzyaabl124,xyzzyaabn124,xyzzyaabo124,xyzzyaabp124&
&,xyzzyaabr124,xyzzyaabs124,xyzzyaabt124,xyzzyaabu124,xyzzyaabv124,xyz&
&zyaabw124,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.fal&
&se.,.false.,.false.)
do xyzzyaaab124=1,3
xyzzyaaac124=xyzzyaaac124+1
lbasis_grad_f(xyzzyaaac124,xyzzyaaah124)=(xyzzyaaao124(xyzzyaaab124)-x&
&yzzyaacl124(xyzzyaaac124))*xyzzyaaby124
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)+xyzzyaaap124*xyz&
&zyaabz124
enddo
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)-xyzzyaaaq124
xyzzyaaab1(xyzzyaaai124,xyzzyaaad124)=0.d0
endif
enddo
enddo
endif
if(xyzzyaaes1)then
do xyzzyaaad124=1,xyzzyaafq1
do xyzzyaaaj124=0,xyzzyaabz1
do xyzzyaaai124=0,xyzzyaaby1
if(xyzzyaaau1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaad124)==1)then
xyzzyaaah124=xyzzyaaah124+1
xyzzyaaac1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaad124)=1.d0
call xyzzyaape1
xyzzyaaaa124=0
xyzzyaaac124=0
do xyzzyaaaf124=1,nspin
xyzzyaabz124=inv_pmass(xyzzyaaaf124)
xyzzyaaby124=xyzzyaack124(xyzzyaaaf124)
do xyzzyaaag124=1,nele(xyzzyaaaf124)
xyzzyaaaa124=xyzzyaaaa124+1
call xyzzyaapq1(xyzzyaaaa124,xyzzyaaaf124,rele_scr(1,xyzzyaaaa124,is),&
&eevecs_scr(1,1,1,is),eevecs_scr(1,1,xyzzyaaaa124,is),eivecs_scr(1,1,1&
&,is),.false.,.true.,.true.,.true.,xyzzyaaan124,xyzzyaaao124,xyzzyaaap&
&124,xyzzyaaar124,xyzzyaaas124,xyzzyaaat124,xyzzyaaav124,xyzzyaaaw124,&
&xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124,xyzzyaaba124,xyzzyaabb124,xyzz&
&yaabc124,xyzzyaabd124,xyzzyaabf124,xyzzyaabg124,xyzzyaabh124,xyzzyaab&
&j124,xyzzyaabk124,xyzzyaabl124,xyzzyaabn124,xyzzyaabo124,xyzzyaabp124&
&,xyzzyaabr124,xyzzyaabs124,xyzzyaabt124,xyzzyaabu124,xyzzyaabv124,xyz&
&zyaabw124,.false.,.true.,.false.,.false.,.false.,.false.,.false.,.fal&
&se.,.false.,.false.)
do xyzzyaaab124=1,3
xyzzyaaac124=xyzzyaaac124+1
lbasis_grad_f(xyzzyaaac124,xyzzyaaah124)=(xyzzyaaas124(xyzzyaaab124)-x&
&yzzyaacm124(xyzzyaaac124))*xyzzyaaby124
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)+xyzzyaaat124*xyz&
&zyaabz124
enddo
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)-xyzzyaaau124
xyzzyaaac1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaad124)=0.d0
endif
enddo
enddo
enddo
endif
if(xyzzyaaew1)then
do xyzzyaaae124=1,xyzzyaabi1
do xyzzyaaad124=1,xyzzyaafv1(xyzzyaaae124)
do xyzzyaaaj124=0,xyzzyaace1(xyzzyaaae124)
if(xyzzyaaav1(xyzzyaaaj124,xyzzyaaad124,xyzzyaaae124)==1)then
xyzzyaaah124=xyzzyaaah124+1
xyzzyaaad1(xyzzyaaaj124,xyzzyaaad124,xyzzyaaae124)=1.d0
call xyzzyaapg1
xyzzyaaaa124=0
xyzzyaaac124=0
do xyzzyaaaf124=1,nspin
xyzzyaabz124=inv_pmass(xyzzyaaaf124)
xyzzyaaby124=xyzzyaack124(xyzzyaaaf124)
do xyzzyaaag124=1,nele(xyzzyaaaf124)
xyzzyaaaa124=xyzzyaaaa124+1
call xyzzyaapq1(xyzzyaaaa124,xyzzyaaaf124,rele_scr(1,xyzzyaaaa124,is),&
&eevecs_scr(1,1,1,is),eevecs_scr(1,1,xyzzyaaaa124,is),eivecs_scr(1,1,1&
&,is),.false.,.true.,.true.,.true.,xyzzyaaan124,xyzzyaaao124,xyzzyaaap&
&124,xyzzyaaar124,xyzzyaaas124,xyzzyaaat124,xyzzyaaav124,xyzzyaaaw124,&
&xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124,xyzzyaaba124,xyzzyaabb124,xyzz&
&yaabc124,xyzzyaabd124,xyzzyaabf124,xyzzyaabg124,xyzzyaabh124,xyzzyaab&
&j124,xyzzyaabk124,xyzzyaabl124,xyzzyaabn124,xyzzyaabo124,xyzzyaabp124&
&,xyzzyaabr124,xyzzyaabs124,xyzzyaabt124,xyzzyaabu124,xyzzyaabv124,xyz&
&zyaabw124,.false.,.false.,.false.,.false.,.true.,.false.,.false.,.fal&
&se.,.false.,.false.)
do xyzzyaaab124=1,3
xyzzyaaac124=xyzzyaaac124+1
lbasis_grad_f(xyzzyaaac124,xyzzyaaah124)=(xyzzyaabc124(xyzzyaaab124)-x&
&yzzyaacn124(xyzzyaaac124))*xyzzyaaby124
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)+xyzzyaabd124*xyz&
&zyaabz124
enddo
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)-xyzzyaabe124
xyzzyaaad1(xyzzyaaaj124,xyzzyaaad124,xyzzyaaae124)=0.d0
endif
enddo
enddo
enddo
endif
if(xyzzyaaex1)then
do xyzzyaaae124=1,xyzzyaabj1
do xyzzyaaad124=1,xyzzyaafw1(xyzzyaaae124)
do xyzzyaaak124=0,xyzzyaacf1(xyzzyaaae124)
do xyzzyaaaj124=0,xyzzyaacg1(xyzzyaaae124)
do xyzzyaaai124=xyzzyaaaj124,xyzzyaacg1(xyzzyaaae124)
if(xyzzyaaaw1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaak124,xyzzyaaad124,xyzz&
&yaaae124)==1)then
xyzzyaaah124=xyzzyaaah124+1
xyzzyaaae1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaak124,xyzzyaaad124,xyzzyaa&
&ae124)=1.d0
call xyzzyaaoo1(xyzzyaaae124,xyzzyaaad124)
call xyzzyaaph1
xyzzyaaaa124=0
xyzzyaaac124=0
do xyzzyaaaf124=1,nspin
xyzzyaabz124=inv_pmass(xyzzyaaaf124)
xyzzyaaby124=xyzzyaack124(xyzzyaaaf124)
do xyzzyaaag124=1,nele(xyzzyaaaf124)
xyzzyaaaa124=xyzzyaaaa124+1
call xyzzyaapq1(xyzzyaaaa124,xyzzyaaaf124,rele_scr(1,xyzzyaaaa124,is),&
&eevecs_scr(1,1,1,is),eevecs_scr(1,1,xyzzyaaaa124,is),eivecs_scr(1,1,1&
&,is),.false.,.true.,.true.,.true.,xyzzyaaan124,xyzzyaaao124,xyzzyaaap&
&124,xyzzyaaar124,xyzzyaaas124,xyzzyaaat124,xyzzyaaav124,xyzzyaaaw124,&
&xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124,xyzzyaaba124,xyzzyaabb124,xyzz&
&yaabc124,xyzzyaabd124,xyzzyaabf124,xyzzyaabg124,xyzzyaabh124,xyzzyaab&
&j124,xyzzyaabk124,xyzzyaabl124,xyzzyaabn124,xyzzyaabo124,xyzzyaabp124&
&,xyzzyaabr124,xyzzyaabs124,xyzzyaabt124,xyzzyaabu124,xyzzyaabv124,xyz&
&zyaabw124,.false.,.false.,.false.,.false.,.false.,.true.,.false.,.fal&
&se.,.false.,.false.)
do xyzzyaaab124=1,3
xyzzyaaac124=xyzzyaaac124+1
lbasis_grad_f(xyzzyaaac124,xyzzyaaah124)=(xyzzyaabg124(xyzzyaaab124)-x&
&yzzyaaco124(xyzzyaaac124))*xyzzyaaby124
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)+xyzzyaabh124*xyz&
&zyaabz124
enddo
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)-xyzzyaabi124
xyzzyaaae1(xyzzyaaai124,xyzzyaaaj124,xyzzyaaak124,xyzzyaaad124,xyzzyaa&
&ae124)=0.d0
endif
enddo
enddo
enddo
enddo
enddo
endif
if(xyzzyaaey1)then
do xyzzyaaad124=1,xyzzyaafr1
do xyzzyaaab124=1,xyzzyaagq1
if(xyzzyaaax1(xyzzyaaab124,xyzzyaaad124)==1)then
xyzzyaaah124=xyzzyaaah124+1
xyzzyaaaf1(xyzzyaaab124,xyzzyaaad124)=1.d0
call xyzzyaapl1
xyzzyaaaa124=0
xyzzyaaac124=0
do xyzzyaaaf124=1,nspin
xyzzyaabz124=inv_pmass(xyzzyaaaf124)
xyzzyaaby124=xyzzyaack124(xyzzyaaaf124)
do xyzzyaaag124=1,nele(xyzzyaaaf124)
xyzzyaaaa124=xyzzyaaaa124+1
call xyzzyaapq1(xyzzyaaaa124,xyzzyaaaf124,rele_scr(1,xyzzyaaaa124,is),&
&eevecs_scr(1,1,1,is),eevecs_scr(1,1,xyzzyaaaa124,is),eivecs_scr(1,1,1&
&,is),.false.,.true.,.true.,.true.,xyzzyaaan124,xyzzyaaao124,xyzzyaaap&
&124,xyzzyaaar124,xyzzyaaas124,xyzzyaaat124,xyzzyaaav124,xyzzyaaaw124,&
&xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124,xyzzyaaba124,xyzzyaabb124,xyzz&
&yaabc124,xyzzyaabd124,xyzzyaabf124,xyzzyaabg124,xyzzyaabh124,xyzzyaab&
&j124,xyzzyaabk124,xyzzyaabl124,xyzzyaabn124,xyzzyaabo124,xyzzyaabp124&
&,xyzzyaabr124,xyzzyaabs124,xyzzyaabt124,xyzzyaabu124,xyzzyaabv124,xyz&
&zyaabw124,.false.,.false.,.false.,.false.,.false.,.false.,.true.,.fal&
&se.,.false.,.false.)
do xyzzyaaal124=1,3
xyzzyaaac124=xyzzyaaac124+1
lbasis_grad_f(xyzzyaaac124,xyzzyaaah124)=(xyzzyaabk124(xyzzyaaal124)-x&
&yzzyaacp124(xyzzyaaac124))*xyzzyaaby124
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)+xyzzyaabl124*xyz&
&zyaabz124
enddo
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)-xyzzyaabm124
xyzzyaaaf1(xyzzyaaab124,xyzzyaaad124)=0.d0
endif
enddo
enddo
endif
if(xyzzyaaez1)then
do xyzzyaaad124=1,xyzzyaafs1
do xyzzyaaab124=1,xyzzyaagr1
if(xyzzyaaay1(xyzzyaaab124,xyzzyaaad124)==1)then
xyzzyaaah124=xyzzyaaah124+1
xyzzyaaag1(xyzzyaaab124,xyzzyaaad124)=1.d0
call xyzzyaapm1
xyzzyaaaa124=0
xyzzyaaac124=0
do xyzzyaaaf124=1,nspin
xyzzyaabz124=inv_pmass(xyzzyaaaf124)
xyzzyaaby124=xyzzyaack124(xyzzyaaaf124)
do xyzzyaaag124=1,nele(xyzzyaaaf124)
xyzzyaaaa124=xyzzyaaaa124+1
call xyzzyaapq1(xyzzyaaaa124,xyzzyaaaf124,rele_scr(1,xyzzyaaaa124,is),&
&eevecs_scr(1,1,1,is),eevecs_scr(1,1,xyzzyaaaa124,is),eivecs_scr(1,1,1&
&,is),.false.,.true.,.true.,.true.,xyzzyaaan124,xyzzyaaao124,xyzzyaaap&
&124,xyzzyaaar124,xyzzyaaas124,xyzzyaaat124,xyzzyaaav124,xyzzyaaaw124,&
&xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124,xyzzyaaba124,xyzzyaabb124,xyzz&
&yaabc124,xyzzyaabd124,xyzzyaabf124,xyzzyaabg124,xyzzyaabh124,xyzzyaab&
&j124,xyzzyaabk124,xyzzyaabl124,xyzzyaabn124,xyzzyaabo124,xyzzyaabp124&
&,xyzzyaabr124,xyzzyaabs124,xyzzyaabt124,xyzzyaabu124,xyzzyaabv124,xyz&
&zyaabw124,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.tr&
&ue.,.false.,.false.)
do xyzzyaaal124=1,3
xyzzyaaac124=xyzzyaaac124+1
lbasis_grad_f(xyzzyaaac124,xyzzyaaah124)=(xyzzyaabo124(xyzzyaaal124)-x&
&yzzyaacq124(xyzzyaaac124))*xyzzyaaby124
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)+xyzzyaabp124*xyz&
&zyaabz124
enddo
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)-xyzzyaabq124
xyzzyaaag1(xyzzyaaab124,xyzzyaaad124)=0.d0
endif
enddo
enddo
endif
if(xyzzyaaff1)then
do xyzzyaaae124=1,xyzzyaabn1
do xyzzyaaam124=1,2
do xyzzyaaaj124=0,xyzzyaacd1
do xyzzyaaak124=xyzzyaaaj124*xyzzyaabw1(xyzzyaaae124),xyzzyaacd1
if(xyzzyaabe1(xyzzyaaae124,xyzzyaaam124,xyzzyaaaj124,xyzzyaaak124)==1)&
&then
xyzzyaaah124=xyzzyaaah124+1
xyzzyaaas1(xyzzyaaae124,xyzzyaaam124,xyzzyaaaj124,xyzzyaaak124)=1.d0
if(xyzzyaabw1(xyzzyaaae124)==1)xyzzyaaas1(xyzzyaaae124,xyzzyaaam124,xy&
&zzyaaak124,xyzzyaaaj124)=1.d0
xyzzyaaaa124=0
xyzzyaaac124=0
do xyzzyaaaf124=1,nspin
xyzzyaabz124=inv_pmass(xyzzyaaaf124)
xyzzyaaby124=xyzzyaack124(xyzzyaaaf124)
do xyzzyaaag124=1,nele(xyzzyaaaf124)
xyzzyaaaa124=xyzzyaaaa124+1
call xyzzyaapq1(xyzzyaaaa124,xyzzyaaaf124,rele_scr(1,xyzzyaaaa124,is),&
&eevecs_scr(1,1,1,is),eevecs_scr(1,1,xyzzyaaaa124,is),eivecs_scr(1,1,1&
&,is),.false.,.true.,.true.,.true.,xyzzyaaan124,xyzzyaaao124,xyzzyaaap&
&124,xyzzyaaar124,xyzzyaaas124,xyzzyaaat124,xyzzyaaav124,xyzzyaaaw124,&
&xyzzyaaax124,xyzzyaaay124,xyzzyaaaz124,xyzzyaaba124,xyzzyaabb124,xyzz&
&yaabc124,xyzzyaabd124,xyzzyaabf124,xyzzyaabg124,xyzzyaabh124,xyzzyaab&
&j124,xyzzyaabk124,xyzzyaabl124,xyzzyaabn124,xyzzyaabo124,xyzzyaabp124&
&,xyzzyaabr124,xyzzyaabs124,xyzzyaabt124,xyzzyaabu124,xyzzyaabv124,xyz&
&zyaabw124,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.fa&
&lse.,.false.,.true.)
do xyzzyaaal124=1,3
xyzzyaaac124=xyzzyaaac124+1
lbasis_grad_f(xyzzyaaac124,xyzzyaaah124)=(xyzzyaabv124(xyzzyaaal124)-x&
&yzzyaacr124(xyzzyaaac124))*xyzzyaaby124
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)+xyzzyaabw124*xyz&
&zyaabz124
enddo
enddo
lbasis_lap_f(xyzzyaaah124)=lbasis_lap_f(xyzzyaaah124)-xyzzyaabx124
xyzzyaaas1(xyzzyaaae124,xyzzyaaam124,xyzzyaaaj124,xyzzyaaak124)=0.d0
if(xyzzyaabw1(xyzzyaaae124)==1)xyzzyaaas1(xyzzyaaae124,xyzzyaaam124,xy&
&zzyaaak124,xyzzyaaaj124)=0.d0
endif
enddo
enddo
enddo
enddo
endif
if(xyzzyaaeq1)then
call dcopy(xyzzyaacs124,xyzzyaacd124(0,1),1,xyzzyaaab1(0,1),1)
call xyzzyaapd1
endif
if(xyzzyaaes1)then
call dcopy(xyzzyaact124,xyzzyaace124(0,0,1),1,xyzzyaaac1(0,0,1),1)
call xyzzyaape1
endif
if(xyzzyaaew1)then
call dcopy(xyzzyaacu124,xyzzyaacf124(0,1,1),1,xyzzyaaad1(0,1,1),1)
call xyzzyaapg1
endif
if(xyzzyaaex1)then
call dcopy(xyzzyaacv124,xyzzyaacg124(0,0,0,1,1),1,xyzzyaaae1(0,0,0,1,1&
&),1)
call xyzzyaaon1
call xyzzyaaph1
endif
if(xyzzyaaey1)then
call dcopy(xyzzyaacw124,xyzzyaach124(1,1),1,xyzzyaaaf1(1,1),1)
call xyzzyaapl1
endif
if(xyzzyaaez1)then
call dcopy(xyzzyaacx124,xyzzyaaci124(1,1),1,xyzzyaaag1(1,1),1)
call xyzzyaapm1
endif
if(xyzzyaaff1)then
call dcopy(xyzzyaacy124,xyzzyaacj124(1,1,1,1),1,xyzzyaaas1(1,1,1,1),1)
call xyzzyaapm1
endif
if(xyzzyaaeq1)deallocate(xyzzyaacd124)
if(xyzzyaaes1)deallocate(xyzzyaace124)
if(xyzzyaaew1)deallocate(xyzzyaacf124)
if(xyzzyaaex1)deallocate(xyzzyaacg124)
if(xyzzyaaey1)deallocate(xyzzyaach124)
if(xyzzyaaez1)deallocate(xyzzyaaci124)
if(xyzzyaaff1)deallocate(xyzzyaacj124)
deallocate(xyzzyaacl124,xyzzyaacm124,xyzzyaacn124,xyzzyaaco124,xyzzyaa&
&cp124,xyzzyaacq124,xyzzyaacr124,xyzzyaack124)
call timer('JASTROW',.false.)
end subroutine get_linear_basis_pjastrow
subroutine xyzzyaarb1
use slaarnabt, only : choose,ddot
implicit none
integer xyzzyaaaa125,xyzzyaaab125,xyzzyaaac125,xyzzyaaad125,xyzzyaaae1&
&25,xyzzyaaaf125,xyzzyaaag125,xyzzyaaah125,xyzzyaaai125
real(dp) xyzzyaaaj125,xyzzyaaak125,xyzzyaaal125,xyzzyaaam125,xyzzyaaan&
&125,xyzzyaaao125,xyzzyaaap125,xyzzyaaaq125
real(dp),allocatable :: xyzzyaaar125(:),xyzzyaaas125(:),xyzzyaaat125(:&
&)
if(periodicity/=3)call errstop('COMPUTE_FT_U_3D','This subroutine shou&
&ld only be called for 3D-periodic systems.')
if(xyzzyaaer1)call errstop('COMPUTE_FT_U_3D','Fourier transform of RPA&
& u not implemented.  Use the standard version of u.')
if(.not.xyzzyaaeq1)call errstop('COMPUTE_FT_U_3D','Requested Fourier t&
&ransform of u, so u term must be present.')
call xyzzyaare1
if(allocated(xyzzyaaiv1))deallocate(xyzzyaaiv1)
allocate(xyzzyaaiv1(xyzzyaaix1,max_spin_pairs),stat=xyzzyaaab125)
call check_alloc(xyzzyaaab125,'COMPUTE_FT_U_3D','1')
xyzzyaaiv1=0.d0
allocate(xyzzyaaas125(0:xyzzyaagx1),stat=xyzzyaaab125)
call check_alloc(xyzzyaaab125,'COMPUTE_FT_U_3D','2')
do xyzzyaaae125=0,xyzzyaagx1
xyzzyaaas125(xyzzyaaae125)=dble(choose(xyzzyaagx1,xyzzyaaae125))*(-xyz&
&zyaach1)**(xyzzyaagx1-xyzzyaaae125)*fourpi
enddo
xyzzyaaaf125=xyzzyaabx1+1
xyzzyaaaa125=xyzzyaabx1+xyzzyaagx1
xyzzyaaai125=xyzzyaagx1+1
allocate(xyzzyaaar125(-1:xyzzyaaaa125),xyzzyaaat125(0:xyzzyaabx1),stat&
&=xyzzyaaab125)
call check_alloc(xyzzyaaab125,'COMPUTE_FT_U_3D','3')
xyzzyaaap125=xyzzyaach1*xyzzyaach1
do xyzzyaaad125=0,xyzzyaaaa125
xyzzyaaap125=xyzzyaaap125*xyzzyaach1
xyzzyaaar125(xyzzyaaad125)=xyzzyaaap125/dble(xyzzyaaad125+3)
enddo
do xyzzyaaah125=0,xyzzyaabx1
xyzzyaaat125(xyzzyaaah125)=ddot(xyzzyaaai125,xyzzyaaas125(0),1,xyzzyaa&
&ar125(xyzzyaaah125),1)
enddo
do xyzzyaaag125=1,max_spin_pairs
xyzzyaaiv1(1,xyzzyaaag125)=xyzzyaaiv1(1,xyzzyaaag125)+ddot(xyzzyaaaf12&
&5,xyzzyaaaj1(0,xyzzyaaag125),1,xyzzyaaat125(0),1)
enddo
do xyzzyaaac125=2,xyzzyaaix1
xyzzyaaaq125=sqrt(xyzzyaaiw1(xyzzyaaac125))
xyzzyaaal125=xyzzyaaaq125*xyzzyaach1
xyzzyaaao125=1.d0/xyzzyaaaq125
xyzzyaaam125=sin(xyzzyaaal125)
xyzzyaaan125=cos(xyzzyaaal125)
xyzzyaaar125(-1)=(1.d0-xyzzyaaan125)*xyzzyaaao125
xyzzyaaar125(0)=(xyzzyaaao125*xyzzyaaam125-xyzzyaach1*xyzzyaaan125)*xy&
&zzyaaao125
xyzzyaaak125=xyzzyaach1
do xyzzyaaad125=1,xyzzyaaaa125
xyzzyaaaj125=xyzzyaaak125
xyzzyaaak125=xyzzyaaak125*xyzzyaach1
xyzzyaaar125(xyzzyaaad125)=((xyzzyaaaj125*xyzzyaaam125-dble(xyzzyaaad1&
&25)*xyzzyaaar125(xyzzyaaad125-2))*dble(xyzzyaaad125+1)*xyzzyaaao125-x&
&yzzyaaak125*xyzzyaaan125)*xyzzyaaao125
enddo
do xyzzyaaah125=0,xyzzyaabx1
xyzzyaaat125(xyzzyaaah125)=ddot(xyzzyaaai125,xyzzyaaas125(0),1,xyzzyaa&
&ar125(xyzzyaaah125),1)*xyzzyaaao125
enddo
do xyzzyaaag125=1,max_spin_pairs
xyzzyaaiv1(xyzzyaaac125,xyzzyaaag125)=xyzzyaaiv1(xyzzyaaac125,xyzzyaaa&
&g125)+ddot(xyzzyaaaf125,xyzzyaaaj1(0,xyzzyaaag125),1,xyzzyaaat125(0),&
&1)
enddo
enddo
deallocate(xyzzyaaar125,xyzzyaaat125,xyzzyaaas125)
end subroutine xyzzyaarb1
subroutine xyzzyaarc1
use slaarnabg,  only : a1,a2,area
use slaarnabq, only : minimum_image
use singleton, only : fftn
implicit none
integer xyzzyaaaa126(2),xyzzyaaab126,xyzzyaaac126,xyzzyaaad126,xyzzyaa&
&ae126,xyzzyaaaf126,xyzzyaaag126,xyzzyaaah126,xyzzyaaai126,xyzzyaaaj12&
&6,xyzzyaaak126(2)
integer,allocatable :: xyzzyaaal126(:)
real(dp) xyzzyaaam126(2),xyzzyaaan126,xyzzyaaao126(2),xyzzyaaap126(2),&
&xyzzyaaaq126,xyzzyaaar126,xyzzyaaas126,xyzzyaaat126,xyzzyaaau126(2),x&
&yzzyaaav126,xyzzyaaaw126,xyzzyaaax126,xyzzyaaay126(2),xyzzyaaaz126(2)&
&,xyzzyaaba126,xyzzyaabb126,xyzzyaabc126,xyzzyaabd126
real(dp),parameter :: xyzzyaabe126=1.d-4
complex(dp),allocatable :: xyzzyaabf126(:,:,:)
if(periodicity/=2)call errstop('COMPUTE_FT_U_2D','This subroutine shou&
&ld only be called for 2D-periodic systems.')
if(dimensionality/=2)call errstop('COMPUTE_FT_U_2D','The Fourier trans&
&form of u(r), which is needed for the finite-size correction to the k&
&inetic energy, is not implemented for quasi-2D systems (i.e., 3D syst&
&ems with 2D periodicity).  If you really think this would be worthwhi&
&le then please discuss with Neil.')
if(xyzzyaaer1)call errstop('COMPUTE_FT_U_2D','Fourier transform of RPA&
& u not implemented.  Use the standard version of u.')
xyzzyaaak126=(/1024,1024/)
xyzzyaaaa126=xyzzyaaak126/2
xyzzyaaay126=a1(1:2)*(1.d0/dble(xyzzyaaak126(1)))
xyzzyaaaz126=a2(1:2)*(1.d0/dble(xyzzyaaak126(2)))
allocate(xyzzyaabf126(0:xyzzyaaak126(1)-1,0:xyzzyaaak126(2)-1,max_spin&
&_pairs),stat=xyzzyaaaf126)
call check_alloc(xyzzyaaaf126,'COMPUTE_FT_U_2D','1')
if(xyzzyaaeq1.and.xyzzyaaes1)then
xyzzyaaaw126=max(xyzzyaach1**2,xyzzyaacj1)
elseif(xyzzyaaeq1)then
xyzzyaaaw126=xyzzyaach1**2
elseif(xyzzyaaes1)then
xyzzyaaaw126=xyzzyaacj1
else
xyzzyaaaw126=0.d0
endif
xyzzyaaap126=0.d0
do xyzzyaaac126=0,xyzzyaaak126(2)-1
xyzzyaaao126=xyzzyaaap126
xyzzyaaap126=xyzzyaaap126+xyzzyaaaz126
do xyzzyaaab126=0,xyzzyaaak126(1)-1
xyzzyaaau126=xyzzyaaao126
xyzzyaaao126=xyzzyaaao126+xyzzyaaay126
call minimum_image(2,1,xyzzyaaau126)
xyzzyaaav126=sum(xyzzyaaau126**2)
if(xyzzyaaav126<xyzzyaaaw126)then
xyzzyaaav126=sqrt(xyzzyaaav126)
do xyzzyaaah126=1,nspin
do xyzzyaaai126=xyzzyaaah126,nspin
xyzzyaaaj126=which_spair(xyzzyaaai126,xyzzyaaah126,levels_spairs)
xyzzyaabd126=0.d0
if(xyzzyaaeq1)then
call xyzzyaapr1(xyzzyaaav126,xyzzyaaah126,xyzzyaaai126,.true.,.false.,&
&.false.,xyzzyaabc126,xyzzyaaaq126,xyzzyaaar126)
xyzzyaabd126=xyzzyaabd126+xyzzyaabc126
endif
if(xyzzyaaes1)then
call xyzzyaapu1(xyzzyaaav126,0.d0,xyzzyaaah126,xyzzyaaai126,.true.,.fa&
&lse.,.false.,xyzzyaabc126,xyzzyaaaq126,xyzzyaaar126,xyzzyaaas126,xyzz&
&yaaat126)
xyzzyaabd126=xyzzyaabd126+xyzzyaabc126
endif
xyzzyaabf126(xyzzyaaab126,xyzzyaaac126,xyzzyaaaj126)=cmplx(xyzzyaabd12&
&6,0.d0,dp)
enddo
enddo
else
xyzzyaabf126(xyzzyaaab126,xyzzyaaac126,:)=czero
endif
enddo
enddo
do xyzzyaaaj126=1,max_spin_pairs
call fftn(xyzzyaabf126(:,:,xyzzyaaaj126),xyzzyaaak126)
enddo
call xyzzyaare1
xyzzyaaba126=xyzzyaaiw1(xyzzyaaix1)
xyzzyaaax126=sqrt(max(dot_product(b1,b1),dot_product(b2,b2)))*xyzzyaab&
&e126
if(allocated(xyzzyaaiv1))deallocate(xyzzyaaiv1)
allocate(xyzzyaaiv1(xyzzyaaix1,max_spin_pairs),xyzzyaaal126(xyzzyaaix1&
&),stat=xyzzyaaaf126)
call check_alloc(xyzzyaaaf126,'COMPUTE_FT_U_2D','2')
xyzzyaaal126=0
xyzzyaaiv1=0.d0
do xyzzyaaac126=0,xyzzyaaak126(2)-1
if(xyzzyaaac126>xyzzyaaaa126(2))then
xyzzyaaap126=dble(xyzzyaaac126-xyzzyaaak126(2))*b2(1:2)
else
xyzzyaaap126=dble(xyzzyaaac126)*b2(1:2)
endif
do xyzzyaaab126=0,xyzzyaaak126(1)-1
if(xyzzyaaab126>xyzzyaaaa126(1))then
xyzzyaaam126=xyzzyaaap126+dble(xyzzyaaab126-xyzzyaaak126(1))*b1(1:2)
else
xyzzyaaam126=xyzzyaaap126+dble(xyzzyaaab126)*b1(1:2)
endif
xyzzyaaan126=xyzzyaaam126(1)*xyzzyaaam126(1)+xyzzyaaam126(2)*xyzzyaaam&
&126(2)
if(xyzzyaaan126<=xyzzyaaba126)then
xyzzyaaad126=1
xyzzyaaae126=xyzzyaaix1
do
xyzzyaaag126=(xyzzyaaad126+xyzzyaaae126)/2
if(xyzzyaaan126<xyzzyaaiw1(xyzzyaaag126))then
xyzzyaaae126=xyzzyaaag126
else
xyzzyaaad126=xyzzyaaag126
endif
if(xyzzyaaae126-xyzzyaaad126<=1)exit
enddo
if(abs(xyzzyaaan126-xyzzyaaiw1(xyzzyaaad126))<abs(xyzzyaaan126-xyzzyaa&
&iw1(xyzzyaaae126)))then
xyzzyaaag126=xyzzyaaad126
else
xyzzyaaag126=xyzzyaaae126
endif
if(abs(xyzzyaaan126-xyzzyaaiw1(xyzzyaaag126))>xyzzyaaax126)call errsto&
&p('COMPUTE_FT_U_2D','Error in bisection search.')
xyzzyaaiv1(xyzzyaaag126,1:max_spin_pairs)=xyzzyaaiv1(xyzzyaaag126,1:ma&
&x_spin_pairs)+dble(xyzzyaabf126(xyzzyaaab126,xyzzyaaac126,1:max_spin_&
&pairs))
xyzzyaaal126(xyzzyaaag126)=xyzzyaaal126(xyzzyaaag126)+1
endif
enddo
enddo
xyzzyaabb126=area/sqrt(dble(product(xyzzyaaak126)))
do xyzzyaaab126=1,xyzzyaaix1
if(xyzzyaaal126(xyzzyaaab126)>0)xyzzyaaiv1(xyzzyaaab126,:)=xyzzyaaiv1(&
&xyzzyaaab126,:)*xyzzyaabb126/dble(xyzzyaaal126(xyzzyaaab126))
enddo
deallocate(xyzzyaabf126,xyzzyaaal126)
end subroutine xyzzyaarc1
subroutine xyzzyaard1
use slaarnabt, only : choose,ddot
implicit none
integer xyzzyaaaa127,xyzzyaaab127,xyzzyaaac127,xyzzyaaad127,xyzzyaaae1&
&27,xyzzyaaaf127,xyzzyaaag127,xyzzyaaah127,xyzzyaaai127
real(dp) xyzzyaaaj127,xyzzyaaak127,xyzzyaaal127,xyzzyaaam127,xyzzyaaan&
&127,xyzzyaaao127,xyzzyaaap127,xyzzyaaaq127
real(dp),allocatable :: xyzzyaaar127(:),xyzzyaaas127(:),xyzzyaaat127(:&
&)
if(periodicity/=1)call errstop('COMPUTE_FT_U_1D','This subroutine shou&
&ld only be called for 1D-periodic systems.')
if(dimensionality/=1)call errstop('COMPUTE_FT_U_2D','The Fourier trans&
&form of u(r), which is needed for the finite-size correction to the k&
&inetic energy, is not implemented for quasi-1D systems (i.e., 3D syst&
&ems with 1D periodicity).  If you really think this would be worthwhi&
&le then please discuss with Neil.')
if(xyzzyaaer1)call errstop('COMPUTE_FT_U_1D','Fourier transform of RPA&
& u not implemented.  Use the standard version of u.')
if(.not.xyzzyaaeq1)call errstop('COMPUTE_FT_U_1D','Requested Fourier t&
&ransform of u, so u term must be present.')
call xyzzyaare1
if(allocated(xyzzyaaiv1))deallocate(xyzzyaaiv1)
allocate(xyzzyaaiv1(xyzzyaaix1,max_spin_pairs),stat=xyzzyaaab127)
call check_alloc(xyzzyaaab127,'COMPUTE_FT_U_1D','1')
xyzzyaaiv1=0.d0
allocate(xyzzyaaas127(0:xyzzyaagx1),stat=xyzzyaaab127)
call check_alloc(xyzzyaaab127,'COMPUTE_FT_U_1D','2')
do xyzzyaaae127=0,xyzzyaagx1
xyzzyaaas127(xyzzyaaae127)=dble(2*choose(xyzzyaagx1,xyzzyaaae127))*(-x&
&yzzyaach1)**(xyzzyaagx1-xyzzyaaae127)
enddo
xyzzyaaaf127=xyzzyaabx1+1
xyzzyaaaa127=xyzzyaabx1+xyzzyaagx1
xyzzyaaai127=xyzzyaagx1+1
allocate(xyzzyaaar127(0:xyzzyaaaa127),xyzzyaaat127(0:xyzzyaabx1),stat=&
&xyzzyaaab127)
call check_alloc(xyzzyaaab127,'COMPUTE_FT_U_1D','3')
xyzzyaaar127(0)=xyzzyaach1
xyzzyaaal127=xyzzyaach1
do xyzzyaaad127=1,xyzzyaaaa127
xyzzyaaal127=xyzzyaaal127*xyzzyaach1
xyzzyaaar127(xyzzyaaad127)=xyzzyaaal127/dble(xyzzyaaad127+1)
enddo
do xyzzyaaah127=0,xyzzyaabx1
xyzzyaaat127(xyzzyaaah127)=ddot(xyzzyaaai127,xyzzyaaas127(0),1,xyzzyaa&
&ar127(xyzzyaaah127),1)
enddo
do xyzzyaaag127=1,max_spin_pairs
xyzzyaaiv1(1,xyzzyaaag127)=xyzzyaaiv1(1,xyzzyaaag127)+ddot(xyzzyaaaf12&
&7,xyzzyaaaj1(0,xyzzyaaag127),1,xyzzyaaat127(0),1)
enddo
do xyzzyaaac127=2,xyzzyaaix1
xyzzyaaaq127=sqrt(xyzzyaaiw1(xyzzyaaac127))
xyzzyaaam127=xyzzyaaaq127*xyzzyaach1
xyzzyaaap127=1.d0/xyzzyaaaq127
xyzzyaaan127=sin(xyzzyaaam127)
xyzzyaaao127=cos(xyzzyaaam127)
xyzzyaaar127(0)=xyzzyaaap127*xyzzyaaan127
xyzzyaaar127(1)=xyzzyaaap127*(xyzzyaaap127*(xyzzyaaao127-1.d0)+xyzzyaa&
&ch1*xyzzyaaan127)
xyzzyaaaj127=xyzzyaach1
do xyzzyaaad127=2,xyzzyaaaa127
xyzzyaaak127=xyzzyaaaj127
xyzzyaaaj127=xyzzyaaaj127*xyzzyaach1
xyzzyaaar127(xyzzyaaad127)=xyzzyaaap127*(dble(xyzzyaaad127)*xyzzyaaap1&
&27*(xyzzyaaak127*xyzzyaaao127-dble(xyzzyaaad127-1)*xyzzyaaar127(xyzzy&
&aaad127-2))+xyzzyaaaj127*xyzzyaaan127)
enddo
do xyzzyaaah127=0,xyzzyaabx1
xyzzyaaat127(xyzzyaaah127)=ddot(xyzzyaaai127,xyzzyaaas127(0),1,xyzzyaa&
&ar127(xyzzyaaah127),1)
enddo
do xyzzyaaag127=1,max_spin_pairs
xyzzyaaiv1(xyzzyaaac127,xyzzyaaag127)=xyzzyaaiv1(xyzzyaaac127,xyzzyaaa&
&g127)+ddot(xyzzyaaaf127,xyzzyaaaj1(0,xyzzyaaag127),1,xyzzyaaat127(0),&
&1)
enddo
enddo
deallocate(xyzzyaaar127,xyzzyaaat127,xyzzyaaas127)
end subroutine xyzzyaard1
subroutine xyzzyaare1
use slaarnaan, only : lattice_generator
use slaarnabt, only : dcopy
implicit none
integer xyzzyaaaa128
integer,allocatable :: xyzzyaaab128(:,:),xyzzyaaac128(:)
real(dp) xyzzyaaad128(3,3)
real(dp),allocatable :: xyzzyaaae128(:,:),xyzzyaaaf128(:)
allocate(xyzzyaaae128(3,xyzzyaaiy1),xyzzyaaaf128((xyzzyaaiy1+3)/2),xyz&
&zyaaab128(3,xyzzyaaiy1),xyzzyaaac128((xyzzyaaiy1+3)/2),stat=xyzzyaaaa&
&128)
call check_alloc(xyzzyaaaa128,'PREPARE_G_STARS','1')
call lattice_generator(xyzzyaaiy1,periodicity,b1,b2,b3,xyzzyaaae128,xy&
&zzyaaaf128,xyzzyaaac128,xyzzyaaad128,xyzzyaaab128,xyzzyaaix1)
deallocate(xyzzyaaae128,xyzzyaaab128,xyzzyaaac128)
if(allocated(xyzzyaaiw1))deallocate(xyzzyaaiw1)
allocate(xyzzyaaiw1(xyzzyaaix1),stat=xyzzyaaaa128)
call check_alloc(xyzzyaaaa128,'PREPARE_G_STARS','2')
call dcopy(xyzzyaaix1,xyzzyaaaf128(1),1,xyzzyaaiw1(1),1)
deallocate(xyzzyaaaf128)
end subroutine xyzzyaare1
subroutine xyzzyaarf1
implicit none
if(allocated(xyzzyaaiw1))deallocate(xyzzyaaiw1)
xyzzyaaix1=-1
if(allocated(xyzzyaaiv1))deallocate(xyzzyaaiv1)
end subroutine xyzzyaarf1
subroutine finite_size_corr_ke_pjastrow
use slaarnabg, only : area,volume,a1
use slaarnabt, only : dcopy
implicit none
integer xyzzyaaaa130,xyzzyaaab130,xyzzyaaac130,xyzzyaaad130,xyzzyaaae1&
&30,xyzzyaaaf130,xyzzyaaag130,xyzzyaaah130,xyzzyaaai130,xyzzyaaaj130
integer,allocatable :: xyzzyaaak130(:)
integer,parameter :: xyzzyaaal130=8
real(dp) xyzzyaaam130,xyzzyaaan130,xyzzyaaao130,xyzzyaaap130,xyzzyaaaq&
&130,xyzzyaaar130,xyzzyaaas130,xyzzyaaat130,xyzzyaaau130,xyzzyaaav130,&
&xyzzyaaaw130,xyzzyaaax130,xyzzyaaay130,xyzzyaaaz130
real(dp),parameter :: xyzzyaaba130=1.d-4
real(dp),allocatable :: xyzzyaabb130(:,:)
logical,parameter :: xyzzyaabc130=.true.
if(.not.(am_master.and.isperiodic))return
if(.not.(xyzzyaaeq1.or.xyzzyaaey1))return
call wout('Finite-size correction to the kinetic energy')
call wout('============================================')
call wout()
if(xyzzyaaeq1)then
if(.not.allocated(xyzzyaaiv1))then
if(periodicity==3)then
call xyzzyaarb1
elseif(periodicity==2)then
call xyzzyaarc1
else
call xyzzyaard1
endif
endif
else
call xyzzyaare1
endif
allocate(xyzzyaabb130(xyzzyaaix1,max_spin_pairs),stat=xyzzyaaaa130)
call check_alloc(xyzzyaaaa130,'FINITE_SIZE_CORR_KE','1')
if(xyzzyaaey1)then
allocate(xyzzyaaak130(xyzzyaaix1),stat=xyzzyaaaa130)
call check_alloc(xyzzyaaaa130,'FINITE_SIZE_CORR_KE','2')
xyzzyaaak130=0
xyzzyaabb130=0.d0
if(periodicity==3)then
xyzzyaaao130=sqrt(max(dot_product(b1,b1),dot_product(b2,b2),dot_produc&
&t(b3,b3)))*xyzzyaaba130
elseif(periodicity==2)then
xyzzyaaao130=sqrt(max(dot_product(b1,b1),dot_product(b2,b2)))*xyzzyaab&
&a130
else
xyzzyaaao130=abs(b1(1))*xyzzyaaba130
endif
xyzzyaaae130=0
do xyzzyaaac130=1,xyzzyaagq1
do xyzzyaaad130=1,xyzzyaagt1(xyzzyaaac130)
xyzzyaaae130=xyzzyaaae130+1
xyzzyaaam130=xyzzyaagh1(xyzzyaaae130)
xyzzyaaai130=1
xyzzyaaah130=xyzzyaaix1
do
xyzzyaaaf130=(xyzzyaaai130+xyzzyaaah130)/2
if(xyzzyaaam130<xyzzyaaiw1(xyzzyaaaf130))then
xyzzyaaah130=xyzzyaaaf130
else
xyzzyaaai130=xyzzyaaaf130
endif
if(xyzzyaaah130-xyzzyaaai130<=1)exit
enddo
if(abs(xyzzyaaam130-xyzzyaaiw1(xyzzyaaai130))<abs(xyzzyaaam130-xyzzyaa&
&iw1(xyzzyaaah130)))then
xyzzyaaaf130=xyzzyaaai130
else
xyzzyaaaf130=xyzzyaaah130
endif
if(abs(xyzzyaaam130-xyzzyaaiw1(xyzzyaaaf130))>xyzzyaaao130)then
call errstop('FINITE_SIZE_CORR_KE','Error finding G vector in FT of p.&
&')
endif
xyzzyaaak130(xyzzyaaaf130)=xyzzyaaak130(xyzzyaaaf130)+2
do xyzzyaaag130=1,max_spin_pairs
xyzzyaabb130(xyzzyaaaf130,xyzzyaaag130)=xyzzyaabb130(xyzzyaaaf130,xyzz&
&yaaag130)+xyzzyaaah1(xyzzyaaac130,xyzzyaaag130)
enddo
enddo
enddo
if(periodicity==3)then
xyzzyaaav130=volume
elseif(periodicity==2)then
xyzzyaaav130=area
else
xyzzyaaav130=abs(a1(1))
endif
do xyzzyaaac130=1,xyzzyaaix1
if(xyzzyaaak130(xyzzyaaac130)>0)then
xyzzyaaan130=xyzzyaaav130/dble(xyzzyaaak130(xyzzyaaac130))
do xyzzyaaag130=1,max_spin_pairs
xyzzyaabb130(xyzzyaaac130,xyzzyaaag130)=xyzzyaabb130(xyzzyaaac130,xyzz&
&yaaag130)*xyzzyaaan130
enddo
endif
enddo
deallocate(xyzzyaaak130)
if(xyzzyaaeq1)xyzzyaabb130=xyzzyaabb130+xyzzyaaiv1
else
if(xyzzyaaeq1)then
call dcopy(xyzzyaaix1*max_spin_pairs,xyzzyaaiv1(1,1),1,xyzzyaabb130(1,&
&1),1)
else
xyzzyaabb130=0.d0
endif
endif
if(xyzzyaaiz1)then
open(unit=xyzzyaagw1,file='ft_of_jastrow.dat',status='replace',iostat=&
&xyzzyaaab130)
if(xyzzyaaab130/=0)call errstop('FINITE_SIZE_CORR_KE','Error opening f&
&t_of_jastrow.dat.')
do xyzzyaaag130=1,max_spin_pairs
if(xyzzyaaag130/=1)write(xyzzyaagw1,*)'&'
do xyzzyaaac130=2,xyzzyaaix1
write(xyzzyaagw1,*)sqrt(xyzzyaaiw1(xyzzyaaac130)),xyzzyaabb130(xyzzyaa&
&ac130,xyzzyaaag130)
enddo
enddo
close(xyzzyaagw1)
call wordwrap('Have written out the spherically averaged Fourier trans&
&form of the two-body Jastrow factor to ft_of_jastrow.dat.')
call wout()
endif
ke_corr=0.d0
xyzzyaaas130=0.d0
do xyzzyaaaj130=1,nspin
xyzzyaaag130=which_spair(xyzzyaaaj130,xyzzyaaaj130,levels_spairs)
call wout('Calculating kinetic energy correction for particles of type&
& '//trim(i2s(xyzzyaaaj130))//'.')
if(periodicity==3)then
if(xyzzyaaal130+1>xyzzyaaix1)call errstop('FINITE_SIZE_CORR_KE','Not e&
&nough stars of G vectors for KE correction.')
xyzzyaaaw130=0.d0
xyzzyaaax130=0.d0
xyzzyaaaz130=0.d0
xyzzyaaay130=0.d0
do xyzzyaaac130=2,xyzzyaaal130+1
xyzzyaaaw130=xyzzyaaaw130+xyzzyaaiw1(xyzzyaaac130)*xyzzyaabb130(xyzzya&
&aac130,xyzzyaaag130)
xyzzyaaax130=xyzzyaaax130+sqrt(xyzzyaaiw1(xyzzyaaac130))
xyzzyaaay130=xyzzyaaay130+xyzzyaaiw1(xyzzyaaac130)
xyzzyaaaz130=xyzzyaaaz130+sqrt(xyzzyaaiw1(xyzzyaaac130))*xyzzyaaiw1(xy&
&zzyaaac130)*xyzzyaabb130(xyzzyaaac130,xyzzyaaag130)
enddo
xyzzyaaaq130=(xyzzyaaay130*xyzzyaaaw130-xyzzyaaax130*xyzzyaaaz130)/(fo&
&urpi*(xyzzyaaax130**2-dble(xyzzyaaal130)*xyzzyaaay130))
xyzzyaaar130=(dble(xyzzyaaal130)*xyzzyaaaz130-xyzzyaaax130*xyzzyaaaw13&
&0)/(fourpi*(xyzzyaaax130**2-dble(xyzzyaaal130)*xyzzyaaay130))
xyzzyaaat130=dble(nele(xyzzyaaaj130))*inv_pmass(xyzzyaaaj130)*(pi/volu&
&me*xyzzyaaaq130+finite_size_const_c*xyzzyaaar130/volume**fourthirds)
ke_corr=ke_corr+xyzzyaaat130
call wout('  Parameter A in model for u(k) : ',xyzzyaaaq130)
call wout('  Parameter B in model for u(k) : ',xyzzyaaar130)
call wout('  Correction to KE (a.u.)       : ',xyzzyaaat130)
xyzzyaaap130=(6.d0*pi**2/volume)**third
xyzzyaaau130=2.d0*inv_pmass(xyzzyaaaj130)/pcharge(xyzzyaaaj130)**2
if(fourpi+xyzzyaaau130*xyzzyaabb130(2,xyzzyaaag130)*xyzzyaaiw1(2)**2/=&
&0.d0)then
xyzzyaaaq130=-xyzzyaabb130(2,xyzzyaaag130)*xyzzyaaiw1(2)/(fourpi+xyzzy&
&aaau130*xyzzyaabb130(2,xyzzyaaag130)*xyzzyaaiw1(2)**2)
else
xyzzyaaaq130=-1.d10
endif
if(xyzzyaaaq130>0.d0.and.xyzzyaaau130/=0.d0)then
xyzzyaaat130=one_over_twopi/xyzzyaaau130*(xyzzyaaap130-atan(sqrt(xyzzy&
&aaau130*xyzzyaaaq130)*xyzzyaaap130)/sqrt(xyzzyaaau130*xyzzyaaaq130))*&
&dble(nele(xyzzyaaaj130))*inv_pmass(xyzzyaaaj130)
xyzzyaaas130=xyzzyaaas130+xyzzyaaat130
endif
elseif(periodicity==2)then
if(xyzzyaaal130+1>xyzzyaaix1)call errstop('FINITE_SIZE_CORR_KE','Not e&
&nough stars of G vectors for KE correction.')
xyzzyaaaw130=0.d0
xyzzyaaax130=0.d0
xyzzyaaaz130=0.d0
xyzzyaaay130=0.d0
do xyzzyaaac130=2,xyzzyaaal130+1
xyzzyaaaw130=xyzzyaaaw130+xyzzyaaiw1(xyzzyaaac130)**0.75d0*xyzzyaabb13&
&0(xyzzyaaac130,xyzzyaaag130)
xyzzyaaax130=xyzzyaaax130+xyzzyaaiw1(xyzzyaaac130)**0.25d0
xyzzyaaay130=xyzzyaaay130+sqrt(xyzzyaaiw1(xyzzyaaac130))
xyzzyaaaz130=xyzzyaaaz130+xyzzyaaiw1(xyzzyaaac130)*xyzzyaabb130(xyzzya&
&aac130,xyzzyaaag130)
enddo
xyzzyaaaq130=(xyzzyaaay130*xyzzyaaaw130-xyzzyaaax130*xyzzyaaaz130)/(xy&
&zzyaaax130**2-dble(xyzzyaaal130)*xyzzyaaay130)
xyzzyaaar130=(dble(xyzzyaaal130)*xyzzyaaaz130-xyzzyaaax130*xyzzyaaaw13&
&0)/(xyzzyaaax130**2-dble(xyzzyaaal130)*xyzzyaaay130)
xyzzyaaat130=dble(nele(xyzzyaaaj130))*inv_pmass(xyzzyaaaj130)*finite_s&
&ize_const_c*xyzzyaaaq130/(fourpi*area**1.25d0)
ke_corr=ke_corr+xyzzyaaat130
call wout('  Parameter a in model for u(k) : ',xyzzyaaaq130)
call wout('  Parameter b in model for u(k) : ',xyzzyaaar130)
call wout('  Correction to KE (a.u.)       : ',xyzzyaaat130)
else
call errstop('FINITE_SIZE_CORR_KE','Extrapolation not implemented for &
&1D-periodic systems.')
endif
enddo
call wout()
if(model_system)then
call wout('Total KE correction (a.u. per particle)    : ',ke_corr/dble&
&(netot))
if(xyzzyaaas130/=0.d0)call wout('which should be approximately equal t&
&o     : ',xyzzyaaas130/dble(netot))
else
call wout('Total KE correction (a.u. per prim. cell)  : ',ke_corr/dble&
&(npcells))
if(xyzzyaaas130/=0.d0)call wout('which should be approximately equal t&
&o     : ',xyzzyaaas130/dble(npcells))
endif
call wout()
deallocate(xyzzyaabb130)
call xyzzyaarf1
if(xyzzyaaiz1.and.xyzzyaabc130)call xyzzyaarg1
end subroutine finite_size_corr_ke_pjastrow
subroutine xyzzyaarg1
use slaarnabg,  only : a1,a2,a3,area,volume
use slaarnabq, only : minimum_image
use singleton, only : fftn
implicit none
integer xyzzyaaaa131(3),xyzzyaaab131,xyzzyaaac131,xyzzyaaad131,xyzzyaa&
&ae131,xyzzyaaaf131,xyzzyaaag131,xyzzyaaah131,xyzzyaaai131,xyzzyaaaj13&
&1,xyzzyaaak131,xyzzyaaal131,xyzzyaaam131(3)
integer,allocatable :: xyzzyaaan131(:)
real(dp) xyzzyaaao131(3),xyzzyaaap131,xyzzyaaaq131(3),xyzzyaaar131(3),&
&xyzzyaaas131(3),xyzzyaaat131,xyzzyaaau131,xyzzyaaav131,xyzzyaaaw131,x&
&yzzyaaax131(3),xyzzyaaay131,xyzzyaaaz131,xyzzyaaba131(3),xyzzyaabb131&
&(3),xyzzyaabc131(3),xyzzyaabd131(3),xyzzyaabe131
real(dp),parameter :: xyzzyaabf131=1.d-4
real(dp),allocatable :: xyzzyaabg131(:,:)
complex(dp),allocatable :: xyzzyaabh131(:,:,:,:)
if(.not.(am_master.and.isperiodic))return
if(.not.(xyzzyaaeq1.or.xyzzyaaey1))return
call wout('Performing FFT of long-ranged two-body Jastrow factor.')
call wout('======================================================')
if(periodicity==3)then
xyzzyaaam131=(/128,128,128/)
call wout('FFT grid dimensions: '//trim(i2s(xyzzyaaam131(1)))//' by '/&
&/trim(i2s(xyzzyaaam131(2)))//' by '//trim(i2s(xyzzyaaam131(3)))//'.')
elseif(periodicity==2)then
xyzzyaaam131=(/1024,1024,1/)
call wout('FFT grid dimensions: '//trim(i2s(xyzzyaaam131(1)))//' by '/&
&/trim(i2s(xyzzyaaam131(2)))//'.')
else
xyzzyaaam131=(/2097152,1,1/)
call wout('FFT grid dimension : '//trim(i2s(xyzzyaaam131(1)))//'.')
endif
call xyzzyaare1
allocate(xyzzyaabh131(0:xyzzyaaam131(1)-1,0:xyzzyaaam131(2)-1,0:xyzzya&
&aam131(3)-1,max_spin_pairs),xyzzyaaan131(xyzzyaaix1),xyzzyaabg131(xyz&
&zyaaix1,max_spin_pairs),stat=xyzzyaaag131)
call check_alloc(xyzzyaaag131,'COMPUTE_FFT_U_AND_P','')
xyzzyaabg131=0.d0
if(periodicity==3)then
xyzzyaaaz131=sqrt(max(dot_product(b1,b1),dot_product(b2,b2),dot_produc&
&t(b3,b3)))*xyzzyaabf131
elseif(periodicity==2)then
xyzzyaaaz131=sqrt(max(dot_product(b1,b1),dot_product(b2,b2)))*xyzzyaab&
&f131
else
xyzzyaaaz131=abs(b1(1))*xyzzyaabf131
endif
xyzzyaaaa131=xyzzyaaam131/2
xyzzyaaba131=a1*(1.d0/dble(xyzzyaaam131(1)))
xyzzyaabb131=a2*(1.d0/dble(xyzzyaaam131(2)))
xyzzyaabc131=a3*(1.d0/dble(xyzzyaaam131(3)))
xyzzyaaav131=0.d0
xyzzyaaaw131=0.d0
xyzzyaaas131=0.d0
do xyzzyaaad131=0,xyzzyaaam131(3)-1
xyzzyaaar131=xyzzyaaas131
xyzzyaaas131=xyzzyaaas131+xyzzyaabc131
do xyzzyaaac131=0,xyzzyaaam131(2)-1
xyzzyaaaq131=xyzzyaaar131
xyzzyaaar131=xyzzyaaar131+xyzzyaabb131
do xyzzyaaab131=0,xyzzyaaam131(1)-1
xyzzyaaax131=xyzzyaaaq131
xyzzyaaaq131=xyzzyaaaq131+xyzzyaaba131
call minimum_image(3,1,xyzzyaaax131)
if(xyzzyaaeq1)xyzzyaaay131=sqrt(sum(xyzzyaaax131**2))
do xyzzyaaaj131=1,nspin
do xyzzyaaak131=xyzzyaaaj131,nspin
xyzzyaaal131=which_spair(xyzzyaaak131,xyzzyaaaj131,levels_spairs)
if(xyzzyaaeq1)call xyzzyaapr1(xyzzyaaay131,xyzzyaaaj131,xyzzyaaak131,.&
&true.,.false.,.false.,xyzzyaaav131,xyzzyaaat131,xyzzyaaau131)
if(xyzzyaaey1)call xyzzyaaqd1(xyzzyaaax131,xyzzyaaaj131,xyzzyaaak131,.&
&true.,.false.,.false.,xyzzyaaaw131,xyzzyaabd131,xyzzyaaat131)
xyzzyaabh131(xyzzyaaab131,xyzzyaaac131,xyzzyaaad131,xyzzyaaal131)=cmpl&
&x(xyzzyaaav131+xyzzyaaaw131,0.d0,dp)
enddo
enddo
enddo
enddo
enddo
do xyzzyaaal131=1,max_spin_pairs
call fftn(xyzzyaabh131(:,:,:,xyzzyaaal131),xyzzyaaam131)
enddo
xyzzyaaan131=0
do xyzzyaaad131=0,xyzzyaaam131(3)-1
if(xyzzyaaad131>xyzzyaaaa131(3))then
xyzzyaaas131=dble(xyzzyaaad131-xyzzyaaam131(3))*b3
else
xyzzyaaas131=dble(xyzzyaaad131)*b3
endif
do xyzzyaaac131=0,xyzzyaaam131(2)-1
if(xyzzyaaac131>xyzzyaaaa131(2))then
xyzzyaaar131=xyzzyaaas131+dble(xyzzyaaac131-xyzzyaaam131(2))*b2
else
xyzzyaaar131=xyzzyaaas131+dble(xyzzyaaac131)*b2
endif
do xyzzyaaab131=0,xyzzyaaam131(1)-1
if(xyzzyaaab131>xyzzyaaaa131(1))then
xyzzyaaao131=xyzzyaaar131+dble(xyzzyaaab131-xyzzyaaam131(1))*b1
else
xyzzyaaao131=xyzzyaaar131+dble(xyzzyaaab131)*b1
endif
xyzzyaaap131=dot_product(xyzzyaaao131,xyzzyaaao131)
if(xyzzyaaap131<=xyzzyaaiw1(xyzzyaaix1))then
xyzzyaaae131=1
xyzzyaaaf131=xyzzyaaix1
do
xyzzyaaah131=(xyzzyaaae131+xyzzyaaaf131)/2
if(xyzzyaaap131<xyzzyaaiw1(xyzzyaaah131))then
xyzzyaaaf131=xyzzyaaah131
else
xyzzyaaae131=xyzzyaaah131
endif
if(xyzzyaaaf131-xyzzyaaae131<=1)exit
enddo
if(abs(xyzzyaaap131-xyzzyaaiw1(xyzzyaaae131))<abs(xyzzyaaap131-xyzzyaa&
&iw1(xyzzyaaaf131)))then
xyzzyaaah131=xyzzyaaae131
else
xyzzyaaah131=xyzzyaaaf131
endif
if(abs(xyzzyaaap131-xyzzyaaiw1(xyzzyaaah131))>xyzzyaaaz131)call errsto&
&p('COMPUTE_FFT_U_AND_P','Error in bisection search.')
xyzzyaabg131(xyzzyaaah131,1:max_spin_pairs)=xyzzyaabg131(xyzzyaaah131,&
&1:max_spin_pairs)+dble(xyzzyaabh131(xyzzyaaab131,xyzzyaaac131,xyzzyaa&
&ad131,1:max_spin_pairs))
xyzzyaaan131(xyzzyaaah131)=xyzzyaaan131(xyzzyaaah131)+1
endif
enddo
enddo
enddo
if(periodicity==3)then
xyzzyaabe131=volume/sqrt(dble(product(xyzzyaaam131)))
elseif(periodicity==2)then
xyzzyaabe131=area/sqrt(dble(product(xyzzyaaam131)))
else
xyzzyaabe131=abs(a1(1))/sqrt(dble(product(xyzzyaaam131)))
endif
do xyzzyaaab131=1,xyzzyaaix1
if(xyzzyaaan131(xyzzyaaab131)>0)xyzzyaabg131(xyzzyaaab131,:)=xyzzyaabg&
&131(xyzzyaaab131,:)*xyzzyaabe131/dble(xyzzyaaan131(xyzzyaaab131))
enddo
if(xyzzyaaiz1)then
open(unit=xyzzyaagw1,file='fft_of_jastrow.dat',status='replace',iostat&
&=xyzzyaaai131)
if(xyzzyaaai131/=0)call errstop('COMPUTE_FFT_U_AND_P','Error opening f&
&ft_of_jastrow.dat.')
do xyzzyaaal131=1,max_spin_pairs
if(xyzzyaaal131>1)write(xyzzyaagw1,*)'&'
do xyzzyaaab131=2,xyzzyaaix1
write(xyzzyaagw1,*)sqrt(xyzzyaaiw1(xyzzyaaab131)),xyzzyaabg131(xyzzyaa&
&ab131,xyzzyaaal131)
enddo
enddo
close(xyzzyaagw1)
call wordwrap('Have written out the FFT of the two-body Jastrow factor&
& to fft_of_jastrow.dat.  This should be compared with the analytical &
&result in ft_of_jastrow.dat.')
call wout()
endif
deallocate(xyzzyaabh131,xyzzyaabg131,xyzzyaaan131)
call xyzzyaarf1
end subroutine xyzzyaarg1
end module slaarnabx

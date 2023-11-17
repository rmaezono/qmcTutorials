module slaarnabc
use casl
use dsp
use slaarnach
use format_utils,only : wout,i2s,r2s,l2s,r2s2,wordwrap
use slaarnabg,    only : dimensionality,periodicity,nitot,wigner_seitz&
&_radius,isperiodic,bmat,volume,area,amat,iontype,nbasis,npcells
use slaarnabt,   only : swap1,ddot,quicksort,dcopy,dscal,resize_pointe&
&r,sort_matrix_symm,sort_matrix_rect,amb2cand_sort_matrix_rect
use parallel,    only : am_master
use run_control, only : errstop_master,errwarn,timer,check_alloc
use store,       only : nspin,no_families,which_fam,which_eqvfam,nele,&
&netot,netot_netot,nitot_netot,three_netot_netot,three_nitot_netot,whi&
&ch_spin,which_ie,heg_nlayers
use slaarnaag,   only : pi,twopi,one_over_pi
implicit none
private
public read_eebasis,read_enbasis,read_gparam_channel,rewrite_eebasis,r&
&ewrite_enbasis,rewrite_gparam_channel,update_eebasis_casl,update_enba&
&sis_casl,update_gparam_channel_casl,get_sig,match_signature,get_sig_f&
&rom_model,parse_model,model_string,init_groups_ee,init_groups_en,gen_&
&sensible_ruleset,digest_rule,build_channels,print_eqns,query_cusp_eeb&
&asis,query_cusp_enbasis,query_eqprod_ee_ee,query_eqprod_en_en,query_e&
&qprod_ee_en,reconstruct_ruleset,sort_indices_matrix,expindx2matrices,&
&insertion_point,apply_insertion_sign,get_sig_ee_only,get_sig_en_only,&
&gbasis_reset_aniso,eebasis_set_aniso,enbasis_set_aniso
public iterate_nuclei_indices,iterate_electron_indices,iterate_nuclei_&
&indices_fix,iterate_electron_indices_fix,sort_indices_grouped_ge,grou&
&p_by_sig,generate_all_indices,get_all_index_properties
public setup_gbasis,finish_gbasis,get_gbasis1_ch,get_gbasis,accept_mov&
&e_gbasis,reset_config_gbasis,setup_gbasis_params,finish_gbasis_params&
&,count_eebasis_params,count_enbasis_params,get_eebasis_params,get_enb&
&asis_params,put_eebasis_params,put_enbasis_params,invalidate_param1_e&
&ebasis,invalidate_param1_enbasis,clear_scratch_gbasis,save_eebasis_pb&
&uffer,save_enbasis_pbuffer,restore_eebasis_pbuffer,restore_enbasis_pb&
&uffer
public ifn1_eebasis,ifn1_enbasis,nfn_eebasis,nfn_enbasis
public eebasis_scr,enbasis_scr,grad_eebasis_scr,grad_enbasis_scr,lap_e&
&ebasis_scr,lap_enbasis_scr,eebasis1_chscr,enbasis1_chscr,grad_eebasis&
&1_chscr,grad_enbasis1_chscr,nzeecut_scr,nzencut_scr,nzeecut1_chscr,nz&
&encut1_chscr,deebasis_scr,denbasis_scr,d2eebasis_scr,d2enbasis_scr,gr&
&adr_eebasis_scr,gradr_enbasis_scr,lapr_eebasis_scr,lapr_enbasis_scr,d&
&eebasis1_chscr,denbasis1_chscr,gradr_eebasis1_chscr,gradr_enbasis1_ch&
&scr
public pflag_det,pflag_fix,pflag_opt,pflag_unset
integer,parameter :: xyzzyaaaa1=1,xyzzyaaab1=2,xyzzyaaac1=3,xyzzyaaad1&
&=4,xyzzyaaae1=5,xyzzyaaaf1=6,xyzzyaaag1=7,xyzzyaaah1=8,xyzzyaaai1=9,x&
&yzzyaaaj1=10,xyzzyaaak1=11,xyzzyaaal1=12,xyzzyaaam1=13,xyzzyaaan1=14,&
&xyzzyaaao1=15,xyzzyaaap1=16,xyzzyaaaq1=17,xyzzyaaar1=18,xyzzyaaas1=19&
&,xyzzyaaat1=20,xyzzyaaau1=21
character(casl_fullkeysize),parameter :: xyzzyaaav1(21)=(/'basis:natur&
&al power               ','basis:cosine                      ','basis:&
&cosine with k-cutoff        ','basis:r/(r+a) power               ','b&
&asis:r/(r^b+a) power             ','basis:1/(r+a) power              &
& ','basis:natural power vectorial     ','basis:natural polynomial    &
&      ','basis:natural polynomial vectorial','basis:RPA              &
&           ','basis:logarithmic cusp            ','basis:dipole cusp &
&                ','basis:half-integer power          ','basis:tilted &
&dipole cusp          ','basis:nu                          ','cutoff:p&
&olynomial                 ','cutoff:alt polynomial             ','cut&
&off:gaussian                   ','cutoff:anisotropic polynomial     '&
&,'cutoff:quasicusp                  ','cutoff:spline                 &
&    '/)
integer :: xyzzyaaaw1=0,xyzzyaaax1=0,nfn_eebasis=0,nfn_enbasis=0
integer,allocatable :: xyzzyaaay1(:),xyzzyaaaz1(:),ifn1_eebasis(:),ifn&
&1_enbasis(:)
integer,allocatable :: xyzzyaaba1(:),xyzzyaabb1(:),xyzzyaabc1(:),xyzzy&
&aabd1(:),xyzzyaabe1(:),xyzzyaabf1(:)
integer,allocatable :: xyzzyaabg1(:),xyzzyaabh1(:),xyzzyaabi1(:),xyzzy&
&aabj1(:),xyzzyaabk1(:),xyzzyaabl1(:)
logical xyzzyaabm1,xyzzyaabn1,xyzzyaabo1,xyzzyaabp1
integer,allocatable :: xyzzyaabq1(:),xyzzyaabr1(:),xyzzyaabs1(:),xyzzy&
&aabt1(:),xyzzyaabu1(:,:,:),xyzzyaabv1(:,:,:),xyzzyaabw1(:),xyzzyaabx1&
&(:),xyzzyaaby1(:),xyzzyaabz1(:)
logical,allocatable :: xyzzyaaca1(:),xyzzyaacb1(:),xyzzyaacc1(:),xyzzy&
&aacd1(:)
integer,allocatable :: xyzzyaace1(:),xyzzyaacf1(:),xyzzyaacg1(:),xyzzy&
&aach1(:),xyzzyaaci1(:),xyzzyaacj1(:)
integer,allocatable :: xyzzyaack1(:),xyzzyaacl1(:),xyzzyaacm1(:),xyzzy&
&aacn1(:),xyzzyaaco1(:),xyzzyaacp1(:)
character(casl_keysize),pointer :: xyzzyaacq1(:),xyzzyaacr1(:)
character(casl_keysize),pointer :: xyzzyaacs1(:),xyzzyaact1(:)
integer,allocatable :: xyzzyaacu1(:),xyzzyaacv1(:)
real(dp),allocatable :: xyzzyaacw1(:),xyzzyaacx1(:),xyzzyaacy1(:),xyzz&
&yaacz1(:),xyzzyaada1(:),xyzzyaadb1(:)
logical,allocatable :: xyzzyaadc1(:),xyzzyaadd1(:),xyzzyaade1(:),xyzzy&
&aadf1(:),xyzzyaadg1(:),xyzzyaadh1(:)
real(dp),allocatable :: xyzzyaadi1(:),xyzzyaadj1(:)
integer,allocatable :: xyzzyaadk1(:),xyzzyaadl1(:)
real(dp),allocatable :: xyzzyaadm1(:),xyzzyaadn1(:)
character(casl_fullkeysize),pointer :: xyzzyaado1(:),xyzzyaadp1(:),xyz&
&zyaadq1(:),xyzzyaadr1(:)
integer,allocatable :: xyzzyaads1(:),xyzzyaadt1(:)
real(dp),allocatable :: xyzzyaadu1(:),xyzzyaadv1(:)
integer,allocatable :: xyzzyaadw1(:),xyzzyaadx1(:)
real(dp),allocatable :: eebasis_scr(:,:,:,:),enbasis_scr(:,:,:,:),grad&
&_eebasis_scr(:,:,:,:,:),grad_enbasis_scr(:,:,:,:,:),lap_eebasis_scr(:&
&,:,:,:),lap_enbasis_scr(:,:,:,:),deebasis_scr(:,:,:,:),denbasis_scr(:&
&,:,:,:),d2eebasis_scr(:,:,:,:),d2enbasis_scr(:,:,:,:),gradr_eebasis_s&
&cr(:,:,:,:),gradr_enbasis_scr(:,:,:,:),lapr_eebasis_scr(:,:,:),lapr_e&
&nbasis_scr(:,:,:)
real(dp),allocatable :: eebasis1_chscr(:,:,:),enbasis1_chscr(:,:,:),gr&
&ad_eebasis1_chscr(:,:,:,:),grad_enbasis1_chscr(:,:,:,:),deebasis1_chs&
&cr(:,:,:),denbasis1_chscr(:,:,:),gradr_eebasis1_chscr(:,:,:),gradr_en&
&basis1_chscr(:,:,:)
logical,allocatable :: nzeecut_scr(:,:,:,:),nzencut_scr(:,:,:,:)
logical,allocatable :: xyzzyaady1(:),xyzzyaadz1(:),xyzzyaaea1(:)
logical,allocatable :: nzeecut1_chscr(:,:,:),nzencut1_chscr(:,:,:)
logical,allocatable :: xyzzyaaeb1(:),xyzzyaaec1(:)
integer,allocatable :: xyzzyaaed1(:),xyzzyaaee1(:),xyzzyaaef1(:),xyzzy&
&aaeg1(:)
real(dp),allocatable :: xyzzyaaeh1(:),xyzzyaaei1(:),xyzzyaaej1(:),xyzz&
&yaaek1(:)
integer,parameter :: pflag_det=-1,pflag_fix=0,pflag_opt=1,pflag_unset=&
&2
logical,parameter :: xyzzyaael1=.false.
contains
subroutine setup_gbasis
implicit none
integer xyzzyaaaa2(2),xyzzyaaab2(2),xyzzyaaac2(2),xyzzyaaad2(2),xyzzya&
&aae2
xyzzyaaaa2=0
xyzzyaaab2=0
xyzzyaaac2=0
xyzzyaaad2=0
call include_range(ratiocfg_from_sz,xyzzyaaaa2)
call include_range(ratiocfg_to_sz,xyzzyaaaa2)
call include_range(ratio1_from_sz,xyzzyaaaa2)
call include_range(ratio1_to_sz,xyzzyaaad2)
call include_range(ratio2_from_sz,xyzzyaaaa2)
call include_range(ratio2_to_sz,xyzzyaaaa2)
call include_range(ratio_ion_from_sz,xyzzyaaaa2)
call include_range(ratio_ion_to_sz,xyzzyaaaa2)
call include_range(drift_sz,xyzzyaaaa2)
call include_range(drift_sz,xyzzyaaab2)
call include_range(kinetic_sz,xyzzyaaaa2)
call include_range(kinetic_sz,xyzzyaaab2)
call include_range(kinetic_sz,xyzzyaaac2)
if(xyzzyaaaa2(1)/=0)then
allocate(eebasis_scr(max(1,nfn_eebasis),netot,netot,xyzzyaaaa2(1):xyzz&
&yaaaa2(2)),enbasis_scr(max(1,nfn_enbasis),max(1,nitot),netot,xyzzyaaa&
&a2(1):xyzzyaaaa2(2)),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','e*basis')
eebasis_scr=1.d0
enbasis_scr=1.d0
allocate(nzeecut_scr(netot,netot,max(1,xyzzyaaaw1),xyzzyaaaa2(1):xyzzy&
&aaaa2(2)),nzencut_scr(max(1,nitot),netot,max(1,xyzzyaaax1),          &
&  xyzzyaaaa2(1):xyzzyaaaa2(2)),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','nze*cut')
nzeecut_scr=.true.
nzencut_scr=.true.
endif
allocate(xyzzyaady1(nscratch),xyzzyaadw1(nscratch),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','gbasis_valid')
if(xyzzyaaad2(1)/=0)then
allocate(eebasis1_chscr(max(1,nfn_eebasis),netot,xyzzyaaad2(1):xyzzyaa&
&ad2(2)),enbasis1_chscr(max(1,nfn_enbasis),max(1,nitot),              &
& xyzzyaaad2(1):xyzzyaaad2(2)),grad_eebasis1_chscr(3,max(1,nfn_eebasis&
&),netot,                    xyzzyaaad2(1):xyzzyaaad2(2)),grad_enbasis&
&1_chscr(3,max(1,nfn_enbasis),max(1,nitot),                    xyzzyaa&
&ad2(1):xyzzyaaad2(2)),deebasis1_chscr(max(1,nfn_eebasis),netot,xyzzya&
&aad2(1):xyzzyaaad2(2)),denbasis1_chscr(max(1,nfn_enbasis),max(1,nitot&
&),                xyzzyaaad2(1):xyzzyaaad2(2)),gradr_eebasis1_chscr(3&
&,netot,xyzzyaaad2(1):xyzzyaaad2(2)),gradr_enbasis1_chscr(3,max(1,nito&
&t),xyzzyaaad2(1):xyzzyaaad2(2)),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','e*basis1')
eebasis1_chscr=0.d0
enbasis1_chscr=0.d0
grad_eebasis1_chscr=0.d0
grad_enbasis1_chscr=0.d0
deebasis1_chscr=0.d0
denbasis1_chscr=0.d0
gradr_eebasis1_chscr=0.d0
gradr_enbasis1_chscr=0.d0
allocate(nzeecut1_chscr(netot,max(1,xyzzyaaaw1),xyzzyaaad2(1):xyzzyaaa&
&d2(2)),nzencut1_chscr(max(1,nitot),max(1,xyzzyaaax1),xyzzyaaad2(1):xy&
&zzyaaad2(2)),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','nze*cut1')
nzeecut1_chscr=.true.
nzencut1_chscr=.true.
endif
allocate(xyzzyaaeb1(nscratch),xyzzyaaec1(nscratch),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','gbasis1_chvalid')
if(xyzzyaaab2(1)/=0)then
allocate(grad_eebasis_scr(3,max(1,nfn_eebasis),netot,netot,           &
&      xyzzyaaab2(1):xyzzyaaab2(2)),grad_enbasis_scr(3,max(1,nfn_enbas&
&is),max(1,nitot),netot,                 xyzzyaaab2(1):xyzzyaaab2(2)),&
&deebasis_scr(max(1,nfn_eebasis),netot,netot,             xyzzyaaab2(1&
&):xyzzyaaab2(2)),denbasis_scr(max(1,nfn_enbasis),max(1,nitot),netot, &
&            xyzzyaaab2(1):xyzzyaaab2(2)),gradr_eebasis_scr(3,netot,ne&
&tot,xyzzyaaab2(1):xyzzyaaab2(2)),gradr_enbasis_scr(3,max(1,nitot),net&
&ot,xyzzyaaab2(1):xyzzyaaab2(2)),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','grad_e*basis')
grad_eebasis_scr=0.d0
grad_enbasis_scr=0.d0
deebasis_scr=0.d0
denbasis_scr=0.d0
gradr_eebasis_scr=0.d0
gradr_enbasis_scr=0.d0
endif
allocate(xyzzyaadz1(nscratch),xyzzyaadx1(nscratch),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','grad_gbasis_valid')
if(xyzzyaaac2(1)/=0)then
allocate(lap_eebasis_scr(max(1,nfn_eebasis),netot,netot,xyzzyaaac2(1):&
&xyzzyaaac2(2)),lap_enbasis_scr(max(1,nfn_enbasis),max(1,nitot),netot,&
&                xyzzyaaac2(1):xyzzyaaac2(2)),d2eebasis_scr(max(1,nfn_&
&eebasis),netot,netot,xyzzyaaac2(1):xyzzyaaac2(2)),d2enbasis_scr(max(1&
&,nfn_enbasis),max(1,nitot),netot,              xyzzyaaac2(1):xyzzyaaa&
&c2(2)),lapr_eebasis_scr(netot,netot,xyzzyaaac2(1):xyzzyaaac2(2)),lapr&
&_enbasis_scr(max(1,nitot),netot,xyzzyaaac2(1):xyzzyaaac2(2)),stat=xyz&
&zyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','lap_e*basis')
lap_eebasis_scr=0.d0
lap_enbasis_scr=0.d0
d2eebasis_scr=0.d0
d2enbasis_scr=0.d0
lapr_eebasis_scr=0.d0
lapr_enbasis_scr=0.d0
endif
allocate(xyzzyaaea1(nscratch),stat=xyzzyaaae2)
call check_alloc(xyzzyaaae2,'SETUP_GBASIS','lap_gbasis_valid')
end subroutine setup_gbasis
subroutine finish_gbasis
implicit none
if(allocated(eebasis_scr))deallocate(eebasis_scr,enbasis_scr)
if(allocated(nzeecut_scr))deallocate(nzeecut_scr,nzencut_scr)
deallocate(xyzzyaady1,xyzzyaadw1)
if(allocated(eebasis1_chscr))deallocate(eebasis1_chscr,enbasis1_chscr)
if(allocated(nzeecut1_chscr))deallocate(nzeecut1_chscr,nzencut1_chscr)
deallocate(xyzzyaaeb1)
if(allocated(grad_eebasis1_chscr))deallocate(grad_eebasis1_chscr,grad_&
&enbasis1_chscr,deebasis1_chscr,denbasis1_chscr,gradr_eebasis1_chscr,g&
&radr_enbasis1_chscr)
deallocate(xyzzyaaec1)
if(allocated(grad_eebasis_scr))deallocate(grad_eebasis_scr,grad_enbasi&
&s_scr,deebasis_scr,denbasis_scr,gradr_eebasis_scr,gradr_enbasis_scr)
deallocate(xyzzyaadz1,xyzzyaadx1)
if(allocated(lap_eebasis_scr))deallocate(lap_eebasis_scr,lap_enbasis_s&
&cr,d2eebasis_scr,d2enbasis_scr,lapr_eebasis_scr,lapr_enbasis_scr)
deallocate(xyzzyaaea1)
end subroutine finish_gbasis
subroutine get_gbasis1_ch(ii,is,fd)
implicit none
integer,intent(in) :: ii,is
logical,intent(in) :: fd
logical xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4
xyzzyaaaa4=.not.xyzzyaaeb1(is)
xyzzyaaab4=fd.and..not.xyzzyaaec1(is)
if(.not.(xyzzyaaaa4.or.xyzzyaaab4))return
call timer('GET_GBASIS1_CH',.true.)
xyzzyaaac4=.true.
if(xyzzyaaaa4)xyzzyaaac4=xyzzyaaac4.and.xyzzyaady1(is)
if(xyzzyaaab4)xyzzyaaac4=xyzzyaaac4.and.xyzzyaadz1(is)
if(xyzzyaaac4)then
if(xyzzyaaaa4)then
if(nfn_eebasis>0)then
call dcopy(nfn_eebasis*netot,eebasis_scr(1,1,ii,is),1,eebasis1_chscr(1&
&,1,is),1)
nzeecut1_chscr(:,:,is)=nzeecut_scr(:,ii,:,is)
endif
if(nfn_enbasis>0)then
call dcopy(nfn_enbasis*nitot,enbasis_scr(1,1,ii,is),1,enbasis1_chscr(1&
&,1,is),1)
nzencut1_chscr(:,:,is)=nzencut_scr(:,ii,:,is)
endif
xyzzyaaeb1(is)=.true.
endif
if(xyzzyaaab4)then
if(nfn_eebasis>0)then
if(xyzzyaabm1)then
call dcopy(3*nfn_eebasis*netot,grad_eebasis_scr(1,1,1,ii,is),1,grad_ee&
&basis1_chscr(1,1,1,is),1)
endif
if(xyzzyaabn1)then
call dcopy(3*netot,gradr_eebasis_scr(1,1,ii,is),1,gradr_eebasis1_chscr&
&(1,1,is),1)
call dcopy(nfn_eebasis*netot,deebasis_scr(1,1,ii,is),1,deebasis1_chscr&
&(1,1,is),1)
endif
endif
if(nfn_enbasis>0)then
if(xyzzyaabo1)then
call dcopy(3*nfn_enbasis*nitot,grad_enbasis_scr(1,1,1,ii,is),1,grad_en&
&basis1_chscr(1,1,1,is),1)
endif
if(xyzzyaabp1)then
call dcopy(3*nitot,gradr_enbasis_scr(1,1,ii,is),1,gradr_enbasis1_chscr&
&(1,1,is),1)
call dcopy(nfn_enbasis*nitot,denbasis_scr(1,1,ii,is),1,denbasis1_chscr&
&(1,1,is),1)
endif
endif
xyzzyaaec1(is)=.true.
endif
else
if(.not.xyzzyaaab4)then
if(xyzzyaaaw1>0)then
call get_eevecs1_ch(ii,is)
call xyzzyaaga1(ii,1,eevecs1_chscr(1,1,is),eebasis1_chscr(1,1,is),nzee&
&cut1_chscr(1,1,is))
endif
if(xyzzyaaax1>0)then
call get_eivecs1_ch(ii,is)
call xyzzyaagb1(ii,1,eivecs1_chscr(1,1,is),enbasis1_chscr(1,1,is),nzen&
&cut1_chscr(1,1,is))
endif
xyzzyaaeb1(is)=.true.
else
if(xyzzyaaaw1>0)then
call get_eevecs1_ch(ii,is)
call xyzzyaagc1(ii,1,eevecs1_chscr(1,1,is),   eebasis1_chscr(1,1,is),g&
&rad_eebasis1_chscr(1,1,1,is),deebasis1_chscr(1,1,is),gradr_eebasis1_c&
&hscr(1,1,is),nzeecut1_chscr(1,1,is))
endif
if(xyzzyaaax1>0)then
call get_eivecs1_ch(ii,is)
call xyzzyaagd1(ii,1,eivecs1_chscr(1,1,is),   enbasis1_chscr(1,1,is),g&
&rad_enbasis1_chscr(1,1,1,is),denbasis1_chscr(1,1,is),gradr_enbasis1_c&
&hscr(1,1,is),nzencut1_chscr(1,1,is))
endif
xyzzyaaeb1(is)=.true.
xyzzyaaec1(is)=.true.
endif
endif
call timer('GET_GBASIS1_CH',.false.)
end subroutine get_gbasis1_ch
recursive subroutine get_gbasis(xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5)
implicit none
integer,intent(in) :: xyzzyaaaa5
logical,intent(in) :: xyzzyaaab5,xyzzyaaac5
integer xyzzyaaad5,xyzzyaaae5,xyzzyaaaf5,xyzzyaaag5,xyzzyaaah5
logical xyzzyaaai5,xyzzyaaaj5,xyzzyaaak5,xyzzyaaal5
xyzzyaaai5=.not.xyzzyaady1(xyzzyaaaa5)
xyzzyaaaj5=xyzzyaaab5.and..not.xyzzyaadz1(xyzzyaaaa5)
xyzzyaaak5=xyzzyaaac5.and..not.xyzzyaaea1(xyzzyaaaa5)
xyzzyaaal5=xyzzyaaaj5.or.xyzzyaaak5
if(.not.(xyzzyaaai5.or.xyzzyaaal5))return
call timer('GET_GBASIS',.true.)
xyzzyaaag5=buffer_move1_from(xyzzyaaaa5)
xyzzyaaah5=buffer_move2_from(xyzzyaaaa5)
if(.not.xyzzyaaak5.and.xyzzyaaag5/=0)then
if((xyzzyaaai5.and.xyzzyaadw1(xyzzyaaaa5)/=xyzzyaaag5).or.(xyzzyaaaj5.&
&and.xyzzyaadx1(xyzzyaaaa5)/=xyzzyaaag5))then
call get_gbasis(xyzzyaaag5,xyzzyaaaj5.and.xyzzyaadx1(xyzzyaaaa5)/=xyzz&
&yaaag5,.false.)
endif
if(xyzzyaaai5.and.xyzzyaadw1(xyzzyaaaa5)/=xyzzyaaag5)then
if(nfn_eebasis>0)then
call dcopy(nfn_eebasis*netot_netot,eebasis_scr(1,1,1,xyzzyaaag5),1,eeb&
&asis_scr(1,1,1,xyzzyaaaa5),1)
nzeecut_scr(:,:,:,xyzzyaaaa5)=nzeecut_scr(:,:,:,xyzzyaaag5)
endif
if(nfn_enbasis>0)then
call dcopy(nfn_enbasis*nitot_netot,enbasis_scr(1,1,1,xyzzyaaag5),1,enb&
&asis_scr(1,1,1,xyzzyaaaa5),1)
nzencut_scr(:,:,:,xyzzyaaaa5)=nzencut_scr(:,:,:,xyzzyaaag5)
endif
endif
if(xyzzyaaaj5.and.xyzzyaadx1(xyzzyaaaa5)/=xyzzyaaag5)then
if(nfn_eebasis>0)then
if(xyzzyaabm1)then
call dcopy(nfn_eebasis*three_netot_netot,grad_eebasis_scr(1,1,1,1,xyzz&
&yaaag5),1,grad_eebasis_scr(1,1,1,1,xyzzyaaaa5),1)
endif
if(xyzzyaabn1)then
call dcopy(three_netot_netot,gradr_eebasis_scr(1,1,1,xyzzyaaag5),1,gra&
&dr_eebasis_scr(1,1,1,xyzzyaaaa5),1)
call dcopy(nfn_eebasis*netot_netot,deebasis_scr(1,1,1,xyzzyaaag5),1,de&
&ebasis_scr(1,1,1,xyzzyaaaa5),1)
endif
endif
if(nfn_enbasis>0)then
if(xyzzyaabo1)then
call dcopy(nfn_enbasis*three_nitot_netot,grad_enbasis_scr(1,1,1,1,xyzz&
&yaaag5),1,grad_enbasis_scr(1,1,1,1,xyzzyaaaa5),1)
endif
if(xyzzyaabp1)then
call dcopy(three_nitot_netot,gradr_enbasis_scr(1,1,1,xyzzyaaag5),1,gra&
&dr_enbasis_scr(1,1,1,xyzzyaaaa5),1)
call dcopy(nfn_enbasis*nitot_netot,denbasis_scr(1,1,1,xyzzyaaag5),1,de&
&nbasis_scr(1,1,1,xyzzyaaaa5),1)
endif
endif
endif
xyzzyaaad5=buffer_move1_from_ii(xyzzyaaaa5)
call get_gbasis1_ch(xyzzyaaad5,xyzzyaaaa5,xyzzyaaaj5)
if(xyzzyaaai5)then
if(nfn_eebasis>0)then
call dcopy(nfn_eebasis*netot,eebasis1_chscr(1,1,xyzzyaaaa5),1,eebasis_&
&scr(1,1,xyzzyaaad5,xyzzyaaaa5),1)
eebasis_scr(:,xyzzyaaad5,:,xyzzyaaaa5)=eebasis1_chscr(:,:,xyzzyaaaa5)
nzeecut_scr(:,xyzzyaaad5,:,xyzzyaaaa5)=nzeecut1_chscr(:,:,xyzzyaaaa5)
nzeecut_scr(xyzzyaaad5,:,:,xyzzyaaaa5)=nzeecut1_chscr(:,:,xyzzyaaaa5)
endif
if(nfn_enbasis>0)then
call dcopy(nfn_enbasis*nitot,enbasis1_chscr(1,1,xyzzyaaaa5),1,enbasis_&
&scr(1,1,xyzzyaaad5,xyzzyaaaa5),1)
nzencut_scr(:,xyzzyaaad5,:,xyzzyaaaa5)=nzencut1_chscr(:,:,xyzzyaaaa5)
endif
xyzzyaady1(xyzzyaaaa5)=.true.
xyzzyaadw1(xyzzyaaaa5)=0
endif
if(xyzzyaaaj5)then
if(nfn_eebasis>0)then
if(xyzzyaabm1)then
call dcopy(3*nfn_eebasis*netot,grad_eebasis1_chscr(1,1,1,xyzzyaaaa5),1&
&,grad_eebasis_scr(1,1,1,xyzzyaaad5,xyzzyaaaa5),1)
do xyzzyaaae5=1,netot
call dcopy(3*nfn_eebasis,grad_eebasis1_chscr(1,1,xyzzyaaae5,xyzzyaaaa5&
&),1,grad_eebasis_scr(1,1,xyzzyaaad5,xyzzyaaae5,xyzzyaaaa5),1)
call dscal(3*nfn_eebasis,-1.d0,grad_eebasis_scr(1,1,xyzzyaaad5,xyzzyaa&
&ae5,xyzzyaaaa5),1)
enddo
endif
if(xyzzyaabn1)then
call dcopy(3*netot,gradr_eebasis1_chscr(1,1,xyzzyaaaa5),1,gradr_eebasi&
&s_scr(1,1,xyzzyaaad5,xyzzyaaaa5),1)
gradr_eebasis_scr(:,xyzzyaaad5,:,xyzzyaaaa5)=-gradr_eebasis1_chscr(:,:&
&,xyzzyaaaa5)
call dcopy(nfn_eebasis*netot,deebasis1_chscr(1,1,xyzzyaaaa5),1,deebasi&
&s_scr(1,1,xyzzyaaad5,xyzzyaaaa5),1)
deebasis_scr(:,xyzzyaaad5,:,xyzzyaaaa5)=deebasis1_chscr(:,:,xyzzyaaaa5&
&)
endif
endif
if(nfn_enbasis>0)then
if(xyzzyaabo1)then
call dcopy(3*nfn_enbasis*nitot,grad_enbasis1_chscr(1,1,1,xyzzyaaaa5),1&
&,grad_enbasis_scr(1,1,1,xyzzyaaad5,xyzzyaaaa5),1)
endif
if(xyzzyaabp1)then
call dcopy(3*nitot,gradr_enbasis1_chscr(1,1,xyzzyaaaa5),1,gradr_enbasi&
&s_scr(1,1,xyzzyaaad5,xyzzyaaaa5),1)
call dcopy(nfn_enbasis*nitot,denbasis1_chscr(1,1,xyzzyaaaa5),1,denbasi&
&s_scr(1,1,xyzzyaaad5,xyzzyaaaa5),1)
endif
endif
xyzzyaadz1(xyzzyaaaa5)=.true.
xyzzyaadx1(xyzzyaaaa5)=0
endif
elseif(.not.(xyzzyaaaj5.or.xyzzyaaak5).and.xyzzyaaah5/=0)then
if(xyzzyaadw1(xyzzyaaaa5)/=xyzzyaaah5)then
call get_gbasis(xyzzyaaah5,.false.,.false.)
if(nfn_eebasis>0)then
call dcopy(nfn_eebasis*netot_netot,eebasis_scr(1,1,1,xyzzyaaah5),1,eeb&
&asis_scr(1,1,1,xyzzyaaaa5),1)
nzeecut_scr(:,:,:,xyzzyaaaa5)=nzeecut_scr(:,:,:,xyzzyaaah5)
endif
if(nfn_enbasis>0)then
call dcopy(nfn_enbasis*nitot_netot,enbasis_scr(1,1,1,xyzzyaaah5),1,enb&
&asis_scr(1,1,1,xyzzyaaaa5),1)
nzencut_scr(:,:,:,xyzzyaaaa5)=nzencut_scr(:,:,:,xyzzyaaah5)
endif
endif
xyzzyaaad5=buffer_move2_from_ii(xyzzyaaaa5)
xyzzyaaae5=buffer_move2_from_jj(xyzzyaaaa5)
if(xyzzyaaaw1>0)then
call get_eevecs(xyzzyaaaa5)
call xyzzyaaga1(xyzzyaaad5,netot,eevecs_scr(1,1,xyzzyaaad5,xyzzyaaaa5)&
&,eebasis_scr(1,1,xyzzyaaad5,xyzzyaaaa5),nzeecut_scr(1,xyzzyaaad5,1,xy&
&zzyaaaa5))
call xyzzyaaga1(xyzzyaaae5,netot,eevecs_scr(1,1,xyzzyaaae5,xyzzyaaaa5)&
&,eebasis_scr(1,1,xyzzyaaae5,xyzzyaaaa5),nzeecut_scr(1,xyzzyaaae5,1,xy&
&zzyaaaa5))
do xyzzyaaaf5=1,netot
call dcopy(nfn_eebasis,eebasis_scr(1,xyzzyaaaf5,xyzzyaaad5,xyzzyaaaa5)&
&,1,eebasis_scr(1,xyzzyaaad5,xyzzyaaaf5,xyzzyaaaa5),1)
enddo
nzeecut_scr(xyzzyaaad5,:,:,xyzzyaaaa5)=nzeecut_scr(:,xyzzyaaad5,:,xyzz&
&yaaaa5)
do xyzzyaaaf5=1,netot
call dcopy(nfn_eebasis,eebasis_scr(1,xyzzyaaaf5,xyzzyaaae5,xyzzyaaaa5)&
&,1,eebasis_scr(1,xyzzyaaae5,xyzzyaaaf5,xyzzyaaaa5),1)
enddo
nzeecut_scr(xyzzyaaae5,:,:,xyzzyaaaa5)=nzeecut_scr(:,xyzzyaaae5,:,xyzz&
&yaaaa5)
endif
if(xyzzyaaax1>0)then
call get_eivecs(xyzzyaaaa5)
call xyzzyaagb1(xyzzyaaad5,netot,eivecs_scr(1,1,xyzzyaaad5,xyzzyaaaa5)&
&,enbasis_scr(1,1,xyzzyaaad5,xyzzyaaaa5),nzencut_scr(1,xyzzyaaad5,1,xy&
&zzyaaaa5))
call xyzzyaagb1(xyzzyaaae5,netot,eivecs_scr(1,1,xyzzyaaae5,xyzzyaaaa5)&
&,enbasis_scr(1,1,xyzzyaaae5,xyzzyaaaa5),nzencut_scr(1,xyzzyaaae5,1,xy&
&zzyaaaa5))
endif
xyzzyaady1(xyzzyaaaa5)=.true.
xyzzyaadw1(xyzzyaaaa5)=0
else
if(.not.xyzzyaaal5)then
if(xyzzyaaaw1>0)then
call get_eevecs(xyzzyaaaa5)
call xyzzyaage1(eevecs_scr(1,1,1,xyzzyaaaa5),eebasis_scr(1,1,1,xyzzyaa&
&aa5),nzeecut_scr(1,1,1,xyzzyaaaa5))
endif
if(xyzzyaaax1>0)then
call get_eivecs(xyzzyaaaa5)
call xyzzyaagf1(eivecs_scr(1,1,1,xyzzyaaaa5),enbasis_scr(1,1,1,xyzzyaa&
&aa5),nzencut_scr(1,1,1,xyzzyaaaa5))
endif
xyzzyaady1(xyzzyaaaa5)=.true.
xyzzyaadw1(xyzzyaaaa5)=0
elseif(.not.xyzzyaaak5)then
if(xyzzyaaaw1>0)then
call get_eevecs(xyzzyaaaa5)
call xyzzyaagg1(eevecs_scr(1,1,1,xyzzyaaaa5),eebasis_scr(1,1,1,xyzzyaa&
&aa5),grad_eebasis_scr(1,1,1,1,xyzzyaaaa5),deebasis_scr(1,1,1,xyzzyaaa&
&a5),gradr_eebasis_scr(1,1,1,xyzzyaaaa5),nzeecut_scr(1,1,1,xyzzyaaaa5)&
&)
endif
if(xyzzyaaax1>0)then
call get_eivecs(xyzzyaaaa5)
call xyzzyaagh1(eivecs_scr(1,1,1,xyzzyaaaa5),enbasis_scr(1,1,1,xyzzyaa&
&aa5),grad_enbasis_scr(1,1,1,1,xyzzyaaaa5),denbasis_scr(1,1,1,xyzzyaaa&
&a5),gradr_enbasis_scr(1,1,1,xyzzyaaaa5),nzencut_scr(1,1,1,xyzzyaaaa5)&
&)
endif
xyzzyaady1(xyzzyaaaa5)=.true.
xyzzyaadw1(xyzzyaaaa5)=0
xyzzyaadz1(xyzzyaaaa5)=.true.
xyzzyaadx1(xyzzyaaaa5)=0
else
if(xyzzyaaaw1>0)then
call get_eevecs(xyzzyaaaa5)
call xyzzyaagi1(eevecs_scr(1,1,1,xyzzyaaaa5),eebasis_scr(1,1,1,xyzzyaa&
&aa5),grad_eebasis_scr(1,1,1,1,xyzzyaaaa5),deebasis_scr(1,1,1,xyzzyaaa&
&a5),gradr_eebasis_scr(1,1,1,xyzzyaaaa5),lap_eebasis_scr(1,1,1,xyzzyaa&
&aa5),d2eebasis_scr(1,1,1,xyzzyaaaa5),lapr_eebasis_scr(1,1,xyzzyaaaa5)&
&,nzeecut_scr(1,1,1,xyzzyaaaa5))
endif
if(xyzzyaaax1>0)then
call get_eivecs(xyzzyaaaa5)
call xyzzyaagj1(eivecs_scr(1,1,1,xyzzyaaaa5),enbasis_scr(1,1,1,xyzzyaa&
&aa5),grad_enbasis_scr(1,1,1,1,xyzzyaaaa5),denbasis_scr(1,1,1,xyzzyaaa&
&a5),gradr_enbasis_scr(1,1,1,xyzzyaaaa5),lap_enbasis_scr(1,1,1,xyzzyaa&
&aa5),d2enbasis_scr(1,1,1,xyzzyaaaa5),lapr_enbasis_scr(1,1,xyzzyaaaa5)&
&,nzencut_scr(1,1,1,xyzzyaaaa5))
endif
xyzzyaady1(xyzzyaaaa5)=.true.
xyzzyaadw1(xyzzyaaaa5)=0
xyzzyaadz1(xyzzyaaaa5)=.true.
xyzzyaadx1(xyzzyaaaa5)=0
xyzzyaaea1(xyzzyaaaa5)=.true.
endif
endif
call timer('GET_GBASIS',.false.)
end subroutine get_gbasis
subroutine accept_move_gbasis(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa6,xyzzyaaab6
if(buffer_move1_from(js)==is)then
xyzzyaaaa6=buffer_move1_from_ii(js)
if(xyzzyaady1(is).and.xyzzyaaeb1(js))then
if(nfn_eebasis>0)then
call dcopy(nfn_eebasis*netot,eebasis1_chscr(1,1,js),1,eebasis_scr(1,1,&
&xyzzyaaaa6,is),1)
eebasis_scr(:,xyzzyaaaa6,:,is)=eebasis1_chscr(:,:,js)
nzeecut_scr(:,xyzzyaaaa6,:,is)=nzeecut1_chscr(:,:,js)
nzeecut_scr(xyzzyaaaa6,:,:,is)=nzeecut1_chscr(:,:,js)
endif
if(nfn_enbasis>0)then
call dcopy(nfn_enbasis*nitot,enbasis1_chscr(1,1,js),1,enbasis_scr(1,1,&
&xyzzyaaaa6,is),1)
nzencut_scr(:,xyzzyaaaa6,:,is)=nzencut1_chscr(:,:,js)
endif
else
xyzzyaady1(is)=.false.
xyzzyaadw1(is)=0
endif
if(xyzzyaadz1(is).and.xyzzyaaec1(js))then
if(nfn_eebasis>0)then
if(xyzzyaabm1)then
call dcopy(3*nfn_eebasis*netot,grad_eebasis1_chscr(1,1,1,js),1,grad_ee&
&basis_scr(1,1,1,xyzzyaaaa6,is),1)
do xyzzyaaab6=1,netot
call dcopy(3*nfn_eebasis,grad_eebasis1_chscr(1,1,xyzzyaaab6,js),1,grad&
&_eebasis_scr(1,1,xyzzyaaaa6,xyzzyaaab6,is),1)
call dscal(3*nfn_eebasis,-1.d0,grad_eebasis_scr(1,1,xyzzyaaaa6,xyzzyaa&
&ab6,is),1)
enddo
endif
if(xyzzyaabn1)then
call dcopy(3*netot,gradr_eebasis1_chscr(1,1,js),1,gradr_eebasis_scr(1,&
&1,xyzzyaaaa6,is),1)
gradr_eebasis_scr(:,xyzzyaaaa6,:,is)=-gradr_eebasis1_chscr(:,:,js)
call dcopy(nfn_eebasis*netot,deebasis1_chscr(1,1,js),1,deebasis_scr(1,&
&1,xyzzyaaaa6,is),1)
deebasis_scr(:,xyzzyaaaa6,:,is)=deebasis1_chscr(:,:,js)
endif
endif
if(nfn_enbasis>0)then
if(xyzzyaabo1)then
call dcopy(3*nfn_enbasis*nitot,grad_enbasis1_chscr(1,1,1,js),1,grad_en&
&basis_scr(1,1,1,xyzzyaaaa6,is),1)
endif
if(xyzzyaabp1)then
call dcopy(3*nitot,gradr_enbasis1_chscr(1,1,js),1,gradr_enbasis_scr(1,&
&1,xyzzyaaaa6,is),1)
call dcopy(nfn_enbasis*nitot,denbasis1_chscr(1,1,js),1,denbasis_scr(1,&
&1,xyzzyaaaa6,is),1)
endif
endif
else
xyzzyaadz1(is)=.false.
xyzzyaadx1(is)=0
endif
xyzzyaaea1(is)=.false.
else
if(xyzzyaady1(js))then
if(xyzzyaaaw1>0)then
call dcopy(nfn_eebasis*netot_netot,eebasis_scr(1,1,1,js),1,eebasis_scr&
&(1,1,1,is),1)
nzeecut_scr(:,:,:,is)=nzeecut_scr(:,:,:,js)
endif
if(xyzzyaaax1>0)then
call dcopy(nfn_enbasis*nitot_netot,enbasis_scr(1,1,1,js),1,enbasis_scr&
&(1,1,1,is),1)
nzencut_scr(:,:,:,is)=nzencut_scr(:,:,:,js)
endif
xyzzyaady1(is)=.true.
xyzzyaadw1(is)=0
else
xyzzyaady1(is)=.false.
xyzzyaadw1(is)=0
endif
if(xyzzyaadz1(js))then
if(xyzzyaaaw1>0)then
if(xyzzyaabm1)then
call dcopy(nfn_eebasis*three_netot_netot,grad_eebasis_scr(1,1,1,1,js),&
&1,grad_eebasis_scr(1,1,1,1,is),1)
endif
if(xyzzyaabn1)then
call dcopy(three_netot_netot,gradr_eebasis_scr(1,1,1,js),1,gradr_eebas&
&is_scr(1,1,1,is),1)
call dcopy(nfn_eebasis*netot_netot,deebasis_scr(1,1,1,js),1,deebasis_s&
&cr(1,1,1,is),1)
endif
endif
if(xyzzyaaax1>0)then
if(xyzzyaabo1)then
call dcopy(nfn_enbasis*three_nitot_netot,grad_enbasis_scr(1,1,1,1,js),&
&1,grad_enbasis_scr(1,1,1,1,is),1)
endif
if(xyzzyaabp1)then
call dcopy(three_nitot_netot,gradr_enbasis_scr(1,1,1,js),1,gradr_enbas&
&is_scr(1,1,1,is),1)
call dcopy(nfn_enbasis*nitot_netot,denbasis_scr(1,1,1,js),1,denbasis_s&
&cr(1,1,1,is),1)
endif
endif
xyzzyaadz1(is)=.true.
xyzzyaadx1(is)=0
else
xyzzyaadz1(is)=.false.
xyzzyaadx1(is)=0
endif
if(xyzzyaaea1(js))then
if(xyzzyaaaw1>0)then
if(xyzzyaabm1)then
call dcopy(nfn_eebasis*netot_netot,lap_eebasis_scr(1,1,1,js),1,lap_eeb&
&asis_scr(1,1,1,is),1)
endif
if(xyzzyaabn1)then
call dcopy(netot_netot,lapr_eebasis_scr(1,1,js),1,lapr_eebasis_scr(1,1&
&,is),1)
call dcopy(nfn_eebasis*netot_netot,d2eebasis_scr(1,1,1,js),1,d2eebasis&
&_scr(1,1,1,is),1)
endif
endif
if(xyzzyaaax1>0)then
if(xyzzyaabo1)then
call dcopy(nfn_enbasis*nitot_netot,lap_enbasis_scr(1,1,1,js),1,lap_enb&
&asis_scr(1,1,1,is),1)
endif
if(xyzzyaabp1)then
call dcopy(nitot_netot,lapr_enbasis_scr(1,1,js),1,lapr_enbasis_scr(1,1&
&,is),1)
call dcopy(nfn_enbasis*nitot_netot,d2enbasis_scr(1,1,1,js),1,d2enbasis&
&_scr(1,1,1,is),1)
endif
endif
xyzzyaaea1(is)=.true.
else
xyzzyaaea1(is)=.false.
endif
endif
xyzzyaaeb1(is)=.false.
xyzzyaaec1(is)=.false.
end subroutine accept_move_gbasis
subroutine reset_config_gbasis(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7
if(xyzzyaadw1(js)==is)then
xyzzyaady1(js)=.false.
xyzzyaadw1(js)=js
elseif(xyzzyaady1(is).and.xyzzyaady1(js))then
if(buffer_move1_from(js)==is)then
xyzzyaaaa7=buffer_move1_from_ii(js)
if(nfn_eebasis>0)then
call dcopy(nfn_eebasis*netot,eebasis_scr(1,1,xyzzyaaaa7,is),1,eebasis_&
&scr(1,1,xyzzyaaaa7,js),1)
eebasis_scr(:,xyzzyaaaa7,:,js)=eebasis_scr(:,xyzzyaaaa7,:,is)
nzeecut_scr(:,xyzzyaaaa7,:,js)=nzeecut_scr(:,xyzzyaaaa7,:,is)
nzeecut_scr(xyzzyaaaa7,:,:,js)=nzeecut_scr(xyzzyaaaa7,:,:,is)
endif
if(nfn_enbasis>0)then
call dcopy(nfn_enbasis*nitot,enbasis_scr(1,1,xyzzyaaaa7,is),1,enbasis_&
&scr(1,1,xyzzyaaaa7,js),1)
nzencut_scr(:,xyzzyaaaa7,:,is)=nzencut_scr(:,xyzzyaaaa7,:,js)
endif
xyzzyaady1(js)=.false.
xyzzyaadw1(js)=js
elseif(buffer_move2_from(js)==is)then
xyzzyaaaa7=buffer_move2_from_ii(js)
xyzzyaaab7=buffer_move2_from_jj(js)
if(nfn_eebasis>0)then
call dcopy(nfn_eebasis*netot,eebasis_scr(1,1,xyzzyaaaa7,is),1,eebasis_&
&scr(1,1,xyzzyaaaa7,js),1)
call dcopy(nfn_eebasis*netot,eebasis_scr(1,1,xyzzyaaab7,is),1,eebasis_&
&scr(1,1,xyzzyaaab7,js),1)
do xyzzyaaac7=1,netot
call dcopy(nfn_eebasis,eebasis_scr(1,xyzzyaaaa7,xyzzyaaac7,is),1,eebas&
&is_scr(1,xyzzyaaaa7,xyzzyaaac7,js),1)
call dcopy(nfn_eebasis,eebasis_scr(1,xyzzyaaab7,xyzzyaaac7,is),1,eebas&
&is_scr(1,xyzzyaaab7,xyzzyaaac7,js),1)
enddo
nzeecut_scr(:,xyzzyaaaa7,:,js)=nzeecut_scr(:,xyzzyaaaa7,:,is)
nzeecut_scr(xyzzyaaaa7,:,:,js)=nzeecut_scr(xyzzyaaaa7,:,:,is)
nzeecut_scr(:,xyzzyaaab7,:,js)=nzeecut_scr(:,xyzzyaaab7,:,is)
nzeecut_scr(xyzzyaaab7,:,:,js)=nzeecut_scr(xyzzyaaab7,:,:,is)
endif
if(nfn_enbasis>0)then
call dcopy(nfn_enbasis*nitot,enbasis_scr(1,1,xyzzyaaaa7,is),1,enbasis_&
&scr(1,1,xyzzyaaaa7,js),1)
call dcopy(nfn_enbasis*nitot,enbasis_scr(1,1,xyzzyaaab7,is),1,enbasis_&
&scr(1,1,xyzzyaaab7,js),1)
nzencut_scr(:,xyzzyaaaa7,:,is)=nzencut_scr(:,xyzzyaaaa7,:,js)
nzencut_scr(:,xyzzyaaab7,:,is)=nzencut_scr(:,xyzzyaaab7,:,js)
endif
xyzzyaady1(js)=.false.
xyzzyaadw1(js)=js
else
xyzzyaady1(js)=.false.
xyzzyaadw1(js)=0
endif
else
xyzzyaady1(js)=.false.
xyzzyaadw1(js)=0
endif
if(xyzzyaadx1(js)==is)then
xyzzyaadz1(js)=.false.
xyzzyaadx1(js)=js
elseif(xyzzyaadz1(is).and.xyzzyaadz1(js))then
if(buffer_move1_from(js)==is)then
xyzzyaaaa7=buffer_move1_from_ii(js)
if(nfn_eebasis>0)then
if(xyzzyaabm1)then
call dcopy(3*nfn_eebasis*netot,grad_eebasis_scr(1,1,1,xyzzyaaaa7,is),1&
&,grad_eebasis_scr(1,1,1,xyzzyaaaa7,js),1)
do xyzzyaaac7=1,netot
call dcopy(3*nfn_eebasis,grad_eebasis_scr(1,1,xyzzyaaaa7,xyzzyaaac7,is&
&),1,grad_eebasis_scr(1,1,xyzzyaaaa7,xyzzyaaac7,js),1)
enddo
endif
if(xyzzyaabn1)then
call dcopy(3*netot,gradr_eebasis_scr(1,1,xyzzyaaaa7,is),1,gradr_eebasi&
&s_scr(1,1,xyzzyaaaa7,js),1)
call dcopy(nfn_eebasis*netot,deebasis_scr(1,1,xyzzyaaaa7,is),1,deebasi&
&s_scr(1,1,xyzzyaaaa7,js),1)
do xyzzyaaac7=1,netot
call dcopy(3,gradr_eebasis_scr(1,xyzzyaaaa7,xyzzyaaac7,is),1,gradr_eeb&
&asis_scr(1,xyzzyaaaa7,xyzzyaaac7,js),1)
call dcopy(nfn_eebasis,deebasis_scr(1,xyzzyaaaa7,xyzzyaaac7,is),1,deeb&
&asis_scr(1,xyzzyaaaa7,xyzzyaaac7,js),1)
enddo
endif
endif
if(nfn_enbasis>0)then
if(xyzzyaabo1)then
call dcopy(3*nfn_enbasis*nitot,grad_enbasis_scr(1,1,1,xyzzyaaaa7,is),1&
&,grad_enbasis_scr(1,1,1,xyzzyaaaa7,js),1)
endif
if(xyzzyaabp1)then
call dcopy(3*nitot,gradr_enbasis_scr(1,1,xyzzyaaaa7,is),1,gradr_enbasi&
&s_scr(1,1,xyzzyaaaa7,js),1)
call dcopy(nfn_enbasis*nitot,denbasis_scr(1,1,xyzzyaaaa7,is),1,denbasi&
&s_scr(1,1,xyzzyaaaa7,js),1)
endif
endif
xyzzyaadz1(js)=.false.
xyzzyaadx1(js)=js
elseif(buffer_move2_from(js)==is)then
xyzzyaaaa7=buffer_move2_from_ii(js)
xyzzyaaab7=buffer_move2_from_jj(js)
if(nfn_eebasis>0)then
if(xyzzyaabm1)then
call dcopy(3*nfn_eebasis*netot,grad_eebasis_scr(1,1,1,xyzzyaaaa7,is),1&
&,grad_eebasis_scr(1,1,1,xyzzyaaaa7,js),1)
call dcopy(3*nfn_eebasis*netot,grad_eebasis_scr(1,1,1,xyzzyaaab7,is),1&
&,grad_eebasis_scr(1,1,1,xyzzyaaab7,js),1)
do xyzzyaaac7=1,netot
call dcopy(3*nfn_eebasis,grad_eebasis_scr(1,1,xyzzyaaaa7,xyzzyaaac7,is&
&),1,grad_eebasis_scr(1,1,xyzzyaaaa7,xyzzyaaac7,js),1)
call dcopy(3*nfn_eebasis,grad_eebasis_scr(1,1,xyzzyaaab7,xyzzyaaac7,is&
&),1,grad_eebasis_scr(1,1,xyzzyaaab7,xyzzyaaac7,js),1)
enddo
endif
if(xyzzyaabn1)then
call dcopy(3*netot,gradr_eebasis_scr(1,1,xyzzyaaaa7,is),1,gradr_eebasi&
&s_scr(1,1,xyzzyaaaa7,js),1)
call dcopy(nfn_eebasis*netot,deebasis_scr(1,1,xyzzyaaaa7,is),1,deebasi&
&s_scr(1,1,xyzzyaaaa7,js),1)
call dcopy(3*netot,gradr_eebasis_scr(1,1,xyzzyaaab7,is),1,gradr_eebasi&
&s_scr(1,1,xyzzyaaab7,js),1)
call dcopy(nfn_eebasis*netot,deebasis_scr(1,1,xyzzyaaab7,is),1,deebasi&
&s_scr(1,1,xyzzyaaab7,js),1)
do xyzzyaaac7=1,netot
call dcopy(3,gradr_eebasis_scr(1,xyzzyaaaa7,xyzzyaaac7,is),1,gradr_eeb&
&asis_scr(1,xyzzyaaaa7,xyzzyaaac7,js),1)
call dcopy(nfn_eebasis,deebasis_scr(1,xyzzyaaaa7,xyzzyaaac7,is),1,deeb&
&asis_scr(1,xyzzyaaaa7,xyzzyaaac7,js),1)
call dcopy(3,gradr_eebasis_scr(1,xyzzyaaab7,xyzzyaaac7,is),1,gradr_eeb&
&asis_scr(1,xyzzyaaab7,xyzzyaaac7,js),1)
call dcopy(nfn_eebasis,deebasis_scr(1,xyzzyaaab7,xyzzyaaac7,is),1,deeb&
&asis_scr(1,xyzzyaaab7,xyzzyaaac7,js),1)
enddo
endif
endif
if(nfn_enbasis>0)then
if(xyzzyaabo1)then
call dcopy(3*nfn_enbasis*nitot,grad_enbasis_scr(1,1,1,xyzzyaaaa7,is),1&
&,grad_enbasis_scr(1,1,1,xyzzyaaaa7,js),1)
call dcopy(3*nfn_enbasis*nitot,grad_enbasis_scr(1,1,1,xyzzyaaab7,is),1&
&,grad_enbasis_scr(1,1,1,xyzzyaaab7,js),1)
endif
if(xyzzyaabp1)then
call dcopy(3*nitot,gradr_enbasis_scr(1,1,xyzzyaaaa7,is),1,gradr_enbasi&
&s_scr(1,1,xyzzyaaaa7,js),1)
call dcopy(nfn_enbasis*nitot,denbasis_scr(1,1,xyzzyaaaa7,is),1,denbasi&
&s_scr(1,1,xyzzyaaaa7,js),1)
call dcopy(3*nitot,gradr_enbasis_scr(1,1,xyzzyaaab7,is),1,gradr_enbasi&
&s_scr(1,1,xyzzyaaab7,js),1)
call dcopy(nfn_enbasis*nitot,denbasis_scr(1,1,xyzzyaaab7,is),1,denbasi&
&s_scr(1,1,xyzzyaaab7,js),1)
endif
endif
xyzzyaadz1(js)=.false.
xyzzyaadx1(js)=js
else
xyzzyaadz1(js)=.false.
xyzzyaadx1(js)=0
endif
else
xyzzyaadz1(js)=.false.
xyzzyaadx1(js)=0
endif
xyzzyaaea1(js)=.false.
xyzzyaaeb1(js)=.false.
xyzzyaaec1(js)=.false.
end subroutine reset_config_gbasis
subroutine count_eebasis_params(iset,nparam)
implicit none
integer,intent(in) :: iset
integer,intent(inout) :: nparam
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8
nparam=0
if(iset<1)return
if(xyzzyaace1(iset)==0)return
xyzzyaaaa8=xyzzyaaay1(iset)
do xyzzyaaac8=1,xyzzyaabw1(iset)
xyzzyaaab8=xyzzyaaba1(xyzzyaaaa8+xyzzyaaac8)
do xyzzyaaad8=1,xyzzyaace1(iset)
if(xyzzyaacu1(xyzzyaaab8+xyzzyaaad8)==pflag_opt)nparam=nparam+1
enddo
enddo
end subroutine count_eebasis_params
subroutine count_enbasis_params(iset,nparam)
implicit none
integer,intent(in) :: iset
integer,intent(inout) :: nparam
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9
nparam=0
if(iset<1)return
if(xyzzyaacf1(iset)==0)return
xyzzyaaaa9=xyzzyaaaz1(iset)
do xyzzyaaac9=1,xyzzyaabx1(iset)
xyzzyaaab9=xyzzyaabb1(xyzzyaaaa9+xyzzyaaac9)
do xyzzyaaad9=1,xyzzyaacf1(iset)
if(xyzzyaacv1(xyzzyaaab9+xyzzyaaad9)==pflag_opt)nparam=nparam+1
enddo
enddo
end subroutine count_enbasis_params
subroutine setup_gbasis_params
implicit none
call xyzzyaaem1
end subroutine setup_gbasis_params
subroutine finish_gbasis_params
implicit none
call xyzzyaaen1
end subroutine finish_gbasis_params
subroutine get_eebasis_params(iset,params,has_lolim,lolim,has_hilim,hi&
&lim,is_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_&
&map,label)
implicit none
integer,intent(in) :: iset
real(dp),intent(inout) :: params(:),lolim(:),hilim(:)
logical,intent(inout) :: has_lolim(:),has_hilim(:),is_shallow(:),is_re&
&dundant(:),is_linear(:),is_loglinear(:),has_aderiv(:),affect_map(:,:)
character(2),intent(inout) :: label(:)
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12
if(iset<1)return
if(xyzzyaace1(iset)==0)return
label=''
is_linear=.false.
is_loglinear=.false.
has_aderiv=.false.
is_redundant=.false.
affect_map=.false.
do xyzzyaaaa12=1,size(affect_map,1)
affect_map(xyzzyaaaa12,xyzzyaaaa12)=.true.
enddo
xyzzyaaab12=0
xyzzyaaac12=xyzzyaaay1(iset)
do xyzzyaaae12=1,xyzzyaabw1(iset)
xyzzyaaad12=xyzzyaaba1(xyzzyaaac12+xyzzyaaae12)
do xyzzyaaaa12=1,xyzzyaace1(iset)
if(xyzzyaacu1(xyzzyaaad12+xyzzyaaaa12)==pflag_opt)then
xyzzyaaab12=xyzzyaaab12+1
params(xyzzyaaab12)=xyzzyaacw1(xyzzyaaad12+xyzzyaaaa12)
has_lolim(xyzzyaaab12)=xyzzyaadc1(xyzzyaaad12+xyzzyaaaa12)
lolim(xyzzyaaab12)=xyzzyaacy1(xyzzyaaad12+xyzzyaaaa12)
has_hilim(xyzzyaaab12)=xyzzyaade1(xyzzyaaad12+xyzzyaaaa12)
hilim(xyzzyaaab12)=xyzzyaada1(xyzzyaaad12+xyzzyaaaa12)
is_shallow(xyzzyaaab12)=xyzzyaadg1(xyzzyaaad12+xyzzyaaaa12)
endif
enddo
enddo
end subroutine get_eebasis_params
subroutine get_enbasis_params(iset,params,has_lolim,lolim,has_hilim,hi&
&lim,is_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_&
&map,label)
implicit none
integer,intent(in) :: iset
real(dp),intent(inout) :: params(:),lolim(:),hilim(:)
logical,intent(inout) :: has_lolim(:),has_hilim(:),is_shallow(:),is_re&
&dundant(:),is_linear(:),is_loglinear(:),has_aderiv(:),affect_map(:,:)
character(2),intent(inout) :: label(:)
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13,xyzzyaaae13
if(iset<1)return
if(xyzzyaacf1(iset)==0)return
label=''
is_linear=.false.
is_loglinear=.false.
has_aderiv=.false.
is_redundant=.false.
affect_map=.false.
do xyzzyaaaa13=1,size(affect_map,1)
affect_map(xyzzyaaaa13,xyzzyaaaa13)=.true.
enddo
xyzzyaaab13=0
xyzzyaaac13=xyzzyaaaz1(iset)
do xyzzyaaae13=1,xyzzyaabx1(iset)
xyzzyaaad13=xyzzyaabb1(xyzzyaaac13+xyzzyaaae13)
do xyzzyaaaa13=1,xyzzyaacf1(iset)
if(xyzzyaacv1(xyzzyaaad13+xyzzyaaaa13)==pflag_opt)then
xyzzyaaab13=xyzzyaaab13+1
params(xyzzyaaab13)=xyzzyaacx1(xyzzyaaad13+xyzzyaaaa13)
has_lolim(xyzzyaaab13)=xyzzyaadd1(xyzzyaaad13+xyzzyaaaa13)
lolim(xyzzyaaab13)=xyzzyaacz1(xyzzyaaad13+xyzzyaaaa13)
has_hilim(xyzzyaaab13)=xyzzyaadf1(xyzzyaaad13+xyzzyaaaa13)
hilim(xyzzyaaab13)=xyzzyaadb1(xyzzyaaad13+xyzzyaaaa13)
is_shallow(xyzzyaaab13)=xyzzyaadh1(xyzzyaaad13+xyzzyaaaa13)
endif
enddo
enddo
end subroutine get_enbasis_params
subroutine put_eebasis_params(iset,params,ignore)
implicit none
integer,intent(in) :: iset
real(dp),intent(inout) :: params(:)
logical,intent(in) :: ignore(:)
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14,xyzzyaaad14,xyzzyaaae14,xy&
&zzyaaaf14,xyzzyaaag14,xyzzyaaah14,xyzzyaaai14
if(iset<1)return
if(xyzzyaace1(iset)==0)return
xyzzyaaaa14=0
xyzzyaaab14=xyzzyaaay1(iset)
xyzzyaaaf14=xyzzyaabi1(iset)+1
xyzzyaaag14=xyzzyaabk1(iset)+1
do xyzzyaaah14=1,xyzzyaabw1(iset)
xyzzyaaac14=xyzzyaaba1(xyzzyaaab14+xyzzyaaah14)
do xyzzyaaai14=1,xyzzyaace1(iset)
if(xyzzyaacu1(xyzzyaaac14+xyzzyaaai14)==pflag_opt)then
xyzzyaaaa14=xyzzyaaaa14+1
if(.not.ignore(xyzzyaaaa14))xyzzyaacw1(xyzzyaaac14+xyzzyaaai14)=params&
&(xyzzyaaaa14)
endif
enddo
xyzzyaaad14=xyzzyaaba1(xyzzyaaab14+xyzzyaaah14)+1
xyzzyaaae14=xyzzyaabg1(xyzzyaaab14+xyzzyaaah14)+1
call xyzzyaafc1(xyzzyaabq1(iset),xyzzyaabs1(iset),xyzzyaace1(iset),xyz&
&zyaack1(iset),xyzzyaacw1(xyzzyaaad14),xyzzyaaco1(iset),xyzzyaads1(xyz&
&zyaaaf14),xyzzyaacm1(iset),xyzzyaadu1(xyzzyaaag14),xyzzyaadi1(xyzzyaa&
&ae14))
enddo
end subroutine put_eebasis_params
subroutine put_enbasis_params(iset,params,ignore)
implicit none
integer,intent(in) :: iset
logical,intent(in) :: ignore(:)
real(dp),intent(inout) :: params(:)
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15,xyzzyaaad15,xyzzyaaae15,xy&
&zzyaaaf15,xyzzyaaag15,xyzzyaaah15,xyzzyaaai15
if(iset<1)return
if(xyzzyaacf1(iset)==0)return
xyzzyaaaa15=0
xyzzyaaab15=xyzzyaaaz1(iset)
xyzzyaaaf15=xyzzyaabj1(iset)+1
xyzzyaaag15=xyzzyaabl1(iset)+1
do xyzzyaaah15=1,xyzzyaabx1(iset)
xyzzyaaac15=xyzzyaabb1(xyzzyaaab15+xyzzyaaah15)
do xyzzyaaai15=1,xyzzyaacf1(iset)
if(xyzzyaacv1(xyzzyaaac15+xyzzyaaai15)==pflag_opt)then
xyzzyaaaa15=xyzzyaaaa15+1
if(.not.ignore(xyzzyaaaa15))xyzzyaacx1(xyzzyaaac15+xyzzyaaai15)=params&
&(xyzzyaaaa15)
endif
enddo
xyzzyaaad15=xyzzyaabb1(xyzzyaaab15+xyzzyaaah15)+1
xyzzyaaae15=xyzzyaabh1(xyzzyaaab15+xyzzyaaah15)+1
call xyzzyaafc1(xyzzyaabr1(iset),xyzzyaabt1(iset),xyzzyaacf1(iset),xyz&
&zyaacl1(iset),xyzzyaacx1(xyzzyaaad15),xyzzyaacp1(iset),xyzzyaadt1(xyz&
&zyaaaf15),xyzzyaacn1(iset),xyzzyaadv1(xyzzyaaag15),xyzzyaadj1(xyzzyaa&
&ae15))
enddo
end subroutine put_enbasis_params
subroutine xyzzyaaem1
implicit none
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16,xyzzyaaae16
if(xyzzyaaaw1>0)then
allocate(xyzzyaaed1(xyzzyaaaw1),xyzzyaaef1(xyzzyaaaw1),stat=xyzzyaaab1&
&6)
call check_alloc(xyzzyaaab16,'SETUP_GBASIS_PBUFFER','iparam0_pbuffer_e&
&ebasis')
xyzzyaaac16=0
xyzzyaaad16=0
do xyzzyaaaa16=1,xyzzyaaaw1
xyzzyaaed1(xyzzyaaaa16)=xyzzyaaac16
xyzzyaaef1(xyzzyaaaa16)=xyzzyaaad16
call count_eebasis_params(xyzzyaaaa16,xyzzyaaae16)
xyzzyaaac16=xyzzyaaac16+(1+xyzzyaaae16)*xyzzyaack1(xyzzyaaaa16)*xyzzya&
&abw1(xyzzyaaaa16)
xyzzyaaad16=xyzzyaaad16+(1+xyzzyaaae16)*xyzzyaace1(xyzzyaaaa16)*xyzzya&
&abw1(xyzzyaaaa16)
enddo
allocate(xyzzyaaeh1(xyzzyaaac16),xyzzyaaej1(xyzzyaaad16),stat=xyzzyaaa&
&b16)
call check_alloc(xyzzyaaab16,'SETUP_GBASIS_PBUFFER','ee')
xyzzyaaeh1=0.d0
xyzzyaaej1=0.d0
endif
if(xyzzyaaax1>0)then
allocate(xyzzyaaee1(xyzzyaaax1),xyzzyaaeg1(xyzzyaaax1),stat=xyzzyaaab1&
&6)
call check_alloc(xyzzyaaab16,'SETUP_GBASIS_PBUFFER','iparam0_calc_pbuf&
&fer_enbasis')
xyzzyaaac16=0
xyzzyaaad16=0
do xyzzyaaaa16=1,xyzzyaaax1
xyzzyaaee1(xyzzyaaaa16)=xyzzyaaac16
xyzzyaaeg1(xyzzyaaaa16)=xyzzyaaad16
call count_enbasis_params(xyzzyaaaa16,xyzzyaaae16)
xyzzyaaac16=xyzzyaaac16+(1+xyzzyaaae16)*xyzzyaacl1(xyzzyaaaa16)*xyzzya&
&abx1(xyzzyaaaa16)
xyzzyaaad16=xyzzyaaad16+(1+xyzzyaaae16)*xyzzyaacf1(xyzzyaaaa16)*xyzzya&
&abx1(xyzzyaaaa16)
enddo
allocate(xyzzyaaei1(xyzzyaaac16),xyzzyaaek1(xyzzyaaad16),stat=xyzzyaaa&
&b16)
call check_alloc(xyzzyaaab16,'SETUP_GBASIS_PBUFFER','en')
xyzzyaaei1=0.d0
xyzzyaaek1=0.d0
endif
end subroutine xyzzyaaem1
subroutine xyzzyaaen1
implicit none
if(xyzzyaaaw1>0)deallocate(xyzzyaaed1,xyzzyaaef1,xyzzyaaeh1,xyzzyaaej1&
&)
if(xyzzyaaax1>0)deallocate(xyzzyaaee1,xyzzyaaeg1,xyzzyaaei1,xyzzyaaek1&
&)
end subroutine xyzzyaaen1
subroutine save_eebasis_pbuffer(iset,iparam)
implicit none
integer,intent(in) :: iset,iparam
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18
if(iset<1)return
if(xyzzyaace1(iset)==0)return
xyzzyaaaa18=xyzzyaaay1(iset)
xyzzyaaab18=xyzzyaabg1(xyzzyaaaa18+1)+1
xyzzyaaac18=xyzzyaaba1(xyzzyaaaa18+1)+1
xyzzyaaad18=xyzzyaaed1(iset)+iparam*xyzzyaack1(iset)*xyzzyaabw1(iset)+&
&1
xyzzyaaae18=xyzzyaaef1(iset)+iparam*xyzzyaace1(iset)*xyzzyaabw1(iset)+&
&1
call dcopy(xyzzyaack1(iset)*xyzzyaabw1(iset),xyzzyaadi1(xyzzyaaab18),1&
&,xyzzyaaeh1(xyzzyaaad18),1)
call dcopy(xyzzyaace1(iset)*xyzzyaabw1(iset),xyzzyaacw1(xyzzyaaac18),1&
&,xyzzyaaej1(xyzzyaaae18),1)
end subroutine save_eebasis_pbuffer
subroutine save_enbasis_pbuffer(iset,iparam)
implicit none
integer,intent(in) :: iset,iparam
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19
if(iset<1)return
if(xyzzyaacf1(iset)==0)return
xyzzyaaaa19=xyzzyaaaz1(iset)
xyzzyaaab19=xyzzyaabh1(xyzzyaaaa19+1)+1
xyzzyaaac19=xyzzyaabb1(xyzzyaaaa19+1)+1
xyzzyaaad19=xyzzyaaee1(iset)+iparam*xyzzyaacl1(iset)*xyzzyaabx1(iset)+&
&1
xyzzyaaae19=xyzzyaaeg1(iset)+iparam*xyzzyaacf1(iset)*xyzzyaabx1(iset)+&
&1
call dcopy(xyzzyaacl1(iset)*xyzzyaabx1(iset),xyzzyaadj1(xyzzyaaab19),1&
&,xyzzyaaei1(xyzzyaaad19),1)
call dcopy(xyzzyaacf1(iset)*xyzzyaabx1(iset),xyzzyaacx1(xyzzyaaac19),1&
&,xyzzyaaek1(xyzzyaaae19),1)
end subroutine save_enbasis_pbuffer
subroutine restore_eebasis_pbuffer(iset,iparam)
implicit none
integer,intent(in) :: iset,iparam
integer xyzzyaaaa20,xyzzyaaab20,xyzzyaaac20,xyzzyaaad20,xyzzyaaae20
if(iset<1)return
if(xyzzyaace1(iset)==0)return
xyzzyaaaa20=xyzzyaaay1(iset)
xyzzyaaab20=xyzzyaabg1(xyzzyaaaa20+1)+1
xyzzyaaac20=xyzzyaaba1(xyzzyaaaa20+1)+1
xyzzyaaad20=xyzzyaaed1(iset)+iparam*xyzzyaack1(iset)*xyzzyaabw1(iset)+&
&1
xyzzyaaae20=xyzzyaaef1(iset)+iparam*xyzzyaace1(iset)*xyzzyaabw1(iset)+&
&1
call dcopy(xyzzyaack1(iset)*xyzzyaabw1(iset),xyzzyaaeh1(xyzzyaaad20),1&
&,xyzzyaadi1(xyzzyaaab20),1)
call dcopy(xyzzyaace1(iset)*xyzzyaabw1(iset),xyzzyaaej1(xyzzyaaae20),1&
&,xyzzyaacw1(xyzzyaaac20),1)
end subroutine restore_eebasis_pbuffer
subroutine restore_enbasis_pbuffer(iset,iparam)
implicit none
integer,intent(in) :: iset,iparam
integer xyzzyaaaa21,xyzzyaaab21,xyzzyaaac21,xyzzyaaad21,xyzzyaaae21
if(iset<1)return
if(xyzzyaacf1(iset)==0)return
xyzzyaaaa21=xyzzyaaaz1(iset)
xyzzyaaab21=xyzzyaabh1(xyzzyaaaa21+1)+1
xyzzyaaac21=xyzzyaabb1(xyzzyaaaa21+1)+1
xyzzyaaad21=xyzzyaaee1(iset)+iparam*xyzzyaacl1(iset)*xyzzyaabx1(iset)+&
&1
xyzzyaaae21=xyzzyaaeg1(iset)+iparam*xyzzyaacf1(iset)*xyzzyaabx1(iset)+&
&1
call dcopy(xyzzyaacl1(iset)*xyzzyaabx1(iset),xyzzyaaei1(xyzzyaaad21),1&
&,xyzzyaadj1(xyzzyaaab21),1)
call dcopy(xyzzyaacf1(iset)*xyzzyaabx1(iset),xyzzyaaek1(xyzzyaaae21),1&
&,xyzzyaacx1(xyzzyaaac21),1)
end subroutine restore_enbasis_pbuffer
subroutine invalidate_param1_eebasis(is,iset,iparam)
implicit none
integer,intent(in) :: is,iset,iparam
if(iset<1)return
xyzzyaady1(is)=.false.
xyzzyaadw1(is)=0
xyzzyaadz1(is)=.false.
xyzzyaadx1(is)=0
xyzzyaaea1(is)=.false.
xyzzyaaeb1(is)=.false.
xyzzyaaec1(is)=.false.
end subroutine invalidate_param1_eebasis
subroutine invalidate_param1_enbasis(is,iset,iparam)
implicit none
integer,intent(in) :: is,iset,iparam
if(iset<1)return
xyzzyaady1(is)=.false.
xyzzyaadw1(is)=0
xyzzyaadz1(is)=.false.
xyzzyaadx1(is)=0
xyzzyaaea1(is)=.false.
xyzzyaaeb1(is)=.false.
xyzzyaaec1(is)=.false.
end subroutine invalidate_param1_enbasis
subroutine clear_scratch_gbasis(is)
implicit none
integer,intent(in) :: is
xyzzyaady1(is)=.false.
xyzzyaadw1(is)=0
xyzzyaadz1(is)=.false.
xyzzyaadx1(is)=0
xyzzyaaea1(is)=.false.
xyzzyaaeb1(is)=.false.
xyzzyaaec1(is)=.false.
end subroutine clear_scratch_gbasis
subroutine read_eebasis(label,namespace,groups,iorder,iref,is_isotropi&
&c,which_unity,require_cutoff,is_cusp_friendly,f_zero,is_none,basis_sy&
&mmetry)
implicit none
integer,intent(in) :: groups(nspin,nspin),iorder
integer,intent(out) :: iref,which_unity,basis_symmetry
character(8),intent(out) :: f_zero
logical,intent(out) :: is_isotropic,require_cutoff,is_cusp_friendly,is&
&_none
character(*),intent(in) :: label,namespace
integer xyzzyaaaa25,xyzzyaaab25,xyzzyaaac25,xyzzyaaad25,xyzzyaaae25,xy&
&zzyaaaf25,xyzzyaaag25,xyzzyaaah25,xyzzyaaai25,xyzzyaaaj25,xyzzyaaak25&
&,xyzzyaaal25,xyzzyaaam25,xyzzyaaan25,xyzzyaaao25(2),xyzzyaaap25(0),xy&
&zzyaaaq25
logical xyzzyaaar25,xyzzyaaas25
character(casl_fullkeysize) btype
character(casl_keysize) chname
integer xyzzyaaat25,xyzzyaaau25,xyzzyaaav25,xyzzyaaaw25,xyzzyaaax25,xy&
&zzyaaay25,xyzzyaaaz25,xyzzyaaba25
integer,pointer :: xyzzyaabb25(:),xyzzyaabc25(:)
real(dp),pointer :: xyzzyaabd25(:),xyzzyaabe25(:),xyzzyaabf25(:),xyzzy&
&aabg25(:),xyzzyaabh25(:)
logical,pointer :: xyzzyaabi25(:),xyzzyaabj25(:),xyzzyaabk25(:)
character(casl_fullkeysize),pointer :: xyzzyaabl25(:),xyzzyaabm25(:)
character(casl_keysize),pointer :: xyzzyaabn25(:)
integer xyzzyaabo25,xyzzyaabp25,xyzzyaabq25,xyzzyaabr25,xyzzyaabs25,xy&
&zzyaabt25,xyzzyaabu25,xyzzyaabv25,xyzzyaabw25,xyzzyaabx25,xyzzyaaby25&
&,xyzzyaabz25,xyzzyaaca25,xyzzyaacb25,xyzzyaacc25,xyzzyaacd25,xyzzyaac&
&e25,xyzzyaacf25,xyzzyaacg25,xyzzyaach25
integer,allocatable :: xyzzyaaci25(:),xyzzyaacj25(:),xyzzyaack25(:,:,:&
&),xyzzyaacl25(:),xyzzyaacm25(:),xyzzyaacn25(:),xyzzyaaco25(:),xyzzyaa&
&cp25(:),xyzzyaacq25(:),xyzzyaacr25(:),xyzzyaacs25(:),xyzzyaact25(:),x&
&yzzyaacu25(:),xyzzyaacv25(:),xyzzyaacw25(:),xyzzyaacx25(:),xyzzyaacy2&
&5(:),xyzzyaacz25(:),xyzzyaada25(:),xyzzyaadb25(:),xyzzyaadc25(:),xyzz&
&yaadd25(:)
logical,allocatable :: xyzzyaade25(:),xyzzyaadf25(:),xyzzyaadg25(:),xy&
&zzyaadh25(:)
real(dp),allocatable :: xyzzyaadi25(:),xyzzyaadj25(:),xyzzyaadk25(:),x&
&yzzyaadl25(:),xyzzyaadm25(:),xyzzyaadn25(:)
character(casl_fullkeysize),pointer :: xyzzyaado25(:),xyzzyaadp25(:)
character(casl_keysize),pointer :: xyzzyaadq25(:),xyzzyaadr25(:)
call get_casl_item(trim(label)//':Type',btype,xyzzyaaak25)
if(xyzzyaaak25/=0)then
if(iorder>1)call errstop_master('READ_EEBASIS','Problem getting "'//tr&
&im(label)//':Type".')
btype='none'
endif
call xyzzyaaez1(trim(namespace)//':'//trim(btype),xyzzyaaat25,is_isotr&
&opic,xyzzyaaar25,xyzzyaaas25,is_cusp_friendly,f_zero,is_none)
if(xyzzyaaat25<0)call errstop_master('READ_EEBASIS','Unrecognized basi&
&s set "'//trim(btype)//'".')
if(xyzzyaaar25.and..not.isperiodic)call errstop_master('READ_EEBASIS',&
&'Basis set "'//trim(btype)//'" cannot be used in non-periodic systems&
&.')
require_cutoff=xyzzyaaas25.and.isperiodic
call xyzzyaafa1(trim(label)//':Constants',xyzzyaaat25,iorder,xyzzyaaav&
&25,xyzzyaabl25,xyzzyaabb25,xyzzyaaaw25,xyzzyaabm25,xyzzyaabd25,xyzzya&
&aay25,xyzzyaabc25,xyzzyaaaz25,xyzzyaabe25,which_unity,basis_symmetry)
call xyzzyaafb1(xyzzyaaat25,iorder,xyzzyaaay25,xyzzyaabc25,xyzzyaaaz25&
&,xyzzyaabe25,xyzzyaaau25,xyzzyaaax25,xyzzyaabn25,xyzzyaabf25,xyzzyaab&
&i25,xyzzyaabg25,xyzzyaabj25,xyzzyaabh25,xyzzyaabk25)
if(xyzzyaaau25==0)then
do xyzzyaaal25=1,xyzzyaaaw1
if(xyzzyaabq1(xyzzyaaal25)/=xyzzyaaat25)cycle
if(xyzzyaaav25/=0)then
xyzzyaaaa25=min(xyzzyaaav25,xyzzyaaci1(xyzzyaaal25))
if(any(xyzzyaabb25(1:xyzzyaaaa25)/=xyzzyaadk1(xyzzyaabc1(xyzzyaaal25)+&
&1:xyzzyaabc1(xyzzyaaal25)+xyzzyaaaa25)))cycle
endif
if(xyzzyaaaw25/=0)then
xyzzyaaaa25=min(xyzzyaaaw25,xyzzyaacg1(xyzzyaaal25))
if(any(xyzzyaabd25(1:xyzzyaaaa25)/=xyzzyaadm1(xyzzyaabe1(xyzzyaaal25)+&
&1:xyzzyaabe1(xyzzyaaal25)+xyzzyaaaa25)))cycle
endif
iref=xyzzyaaal25
if(iorder>xyzzyaabs1(xyzzyaaal25))then
nfn_eebasis=nfn_eebasis+iorder-xyzzyaabs1(xyzzyaaal25)
do xyzzyaaaq25=xyzzyaaal25+1,xyzzyaaaw1
ifn1_eebasis(xyzzyaaaq25)=ifn1_eebasis(xyzzyaaaq25)+iorder-xyzzyaabs1(&
&xyzzyaaal25)
enddo
xyzzyaabs1(xyzzyaaal25)=iorder
endif
if(xyzzyaaav25>0)then
if(xyzzyaaci1(xyzzyaaal25)<xyzzyaaav25)then
xyzzyaabx25=sum(xyzzyaaci1)
xyzzyaabw25=xyzzyaabx25+xyzzyaaav25-xyzzyaaci1(xyzzyaaal25)
allocate(xyzzyaadc25(xyzzyaabx25),xyzzyaado25(xyzzyaabx25),stat=xyzzya&
&aam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','(tmp_const_int')
xyzzyaadc25=xyzzyaadk1
xyzzyaado25=xyzzyaado1
deallocate(xyzzyaadk1,xyzzyaado1)
allocate(xyzzyaadk1(xyzzyaabw25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','const_int')
xyzzyaadk1(1:xyzzyaabc1(xyzzyaaal25))=xyzzyaadc25(1:xyzzyaabc1(xyzzyaa&
&al25))
xyzzyaadk1(xyzzyaabc1(xyzzyaaal25)+1:xyzzyaabc1(xyzzyaaal25)+xyzzyaaav&
&25)=xyzzyaabb25(:)
xyzzyaadk1(xyzzyaabc1(xyzzyaaal25)+xyzzyaaav25:xyzzyaabw25)=xyzzyaadc2&
&5(xyzzyaabc1(xyzzyaaal25)+xyzzyaaci1(xyzzyaaal25):xyzzyaabx25)
xyzzyaado1(1:xyzzyaabc1(xyzzyaaal25))=xyzzyaado25(1:xyzzyaabc1(xyzzyaa&
&al25))
xyzzyaado1(xyzzyaabc1(xyzzyaaal25)+1:xyzzyaabc1(xyzzyaaal25)+xyzzyaaav&
&25)=xyzzyaabl25(:)
xyzzyaado1(xyzzyaabc1(xyzzyaaal25)+xyzzyaaav25:xyzzyaabw25)=xyzzyaado2&
&5(xyzzyaabc1(xyzzyaaal25)+xyzzyaaci1(xyzzyaaal25):xyzzyaabx25)
deallocate(xyzzyaadc25,xyzzyaado25)
xyzzyaabc1(xyzzyaaal25+1:)=xyzzyaabc1(xyzzyaaal25+1:)+xyzzyaaav25-xyzz&
&yaaci1(xyzzyaaal25)
xyzzyaaci1(xyzzyaaal25)=xyzzyaaav25
endif
deallocate(xyzzyaabb25,xyzzyaabl25)
if(xyzzyaaco1(xyzzyaaal25)<xyzzyaaay25)then
xyzzyaacd25=sum(xyzzyaaco1)
xyzzyaacc25=xyzzyaacd25+xyzzyaaay25-xyzzyaaco1(xyzzyaaal25)
allocate(xyzzyaadd25(xyzzyaacd25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp_const_int_calc')
xyzzyaadd25=xyzzyaads1
deallocate(xyzzyaads1)
allocate(xyzzyaads1(xyzzyaacc25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','const_int_calc')
xyzzyaads1(1:xyzzyaabi1(xyzzyaaal25))=xyzzyaadd25(1:xyzzyaabi1(xyzzyaa&
&al25))
xyzzyaads1(xyzzyaabi1(xyzzyaaal25)+1:xyzzyaabi1(xyzzyaaal25)+xyzzyaaay&
&25)=xyzzyaabc25(:)
xyzzyaads1(xyzzyaabi1(xyzzyaaal25)+xyzzyaaay25:xyzzyaacc25)=xyzzyaadd2&
&5(xyzzyaabi1(xyzzyaaal25)+xyzzyaaco1(xyzzyaaal25):xyzzyaacd25)
deallocate(xyzzyaadd25)
xyzzyaabi1(xyzzyaaal25+1:)=xyzzyaabi1(xyzzyaaal25+1:)+xyzzyaaay25-xyzz&
&yaaco1(xyzzyaaal25)
xyzzyaaco1(xyzzyaaal25)=xyzzyaaay25
endif
deallocate(xyzzyaabc25)
endif
if(xyzzyaaaw25>0)then
if(xyzzyaacg1(xyzzyaaal25)<xyzzyaaaw25)then
xyzzyaabz25=sum(xyzzyaacg1)
xyzzyaaby25=xyzzyaabz25+xyzzyaaaw25-xyzzyaacg1(xyzzyaaal25)
allocate(xyzzyaadm25(xyzzyaabz25),xyzzyaadp25(xyzzyaabz25),stat=xyzzya&
&aam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp_const_dble')
xyzzyaadm25=xyzzyaadm1
xyzzyaadp25=xyzzyaadq1
deallocate(xyzzyaadm1,xyzzyaadq1)
allocate(xyzzyaadm1(xyzzyaaby25),xyzzyaadq1(xyzzyaaby25),stat=xyzzyaaa&
&m25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','const_dble')
xyzzyaadm1(1:xyzzyaabe1(xyzzyaaal25))=xyzzyaadm25(1:xyzzyaabe1(xyzzyaa&
&al25))
xyzzyaadm1(xyzzyaabe1(xyzzyaaal25)+1:xyzzyaabe1(xyzzyaaal25)+xyzzyaaaw&
&25)=xyzzyaabd25(:)
xyzzyaadm1(xyzzyaabe1(xyzzyaaal25)+xyzzyaaaw25:xyzzyaaby25)=xyzzyaadm2&
&5(xyzzyaabe1(xyzzyaaal25)+xyzzyaacg1(xyzzyaaal25):xyzzyaabz25)
xyzzyaadq1(1:xyzzyaabe1(xyzzyaaal25))=xyzzyaadp25(1:xyzzyaabe1(xyzzyaa&
&al25))
xyzzyaadq1(xyzzyaabe1(xyzzyaaal25)+1:xyzzyaabe1(xyzzyaaal25)+xyzzyaaaw&
&25)=xyzzyaabm25(:)
xyzzyaadq1(xyzzyaabe1(xyzzyaaal25)+xyzzyaaaw25:xyzzyaaby25)=xyzzyaadp2&
&5(xyzzyaabe1(xyzzyaaal25)+xyzzyaacg1(xyzzyaaal25):xyzzyaabz25)
deallocate(xyzzyaadm25,xyzzyaadp25)
xyzzyaabe1(xyzzyaaal25+1:)=xyzzyaabe1(xyzzyaaal25+1:)+xyzzyaaaw25-xyzz&
&yaacg1(xyzzyaaal25)
xyzzyaacg1(xyzzyaaal25)=xyzzyaaaw25
endif
deallocate(xyzzyaabd25,xyzzyaabm25)
if(xyzzyaacm1(xyzzyaaal25)<xyzzyaaaz25)then
xyzzyaacf25=sum(xyzzyaacm1)
xyzzyaace25=xyzzyaacf25+xyzzyaaaz25-xyzzyaacm1(xyzzyaaal25)
allocate(xyzzyaadn25(xyzzyaacf25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp_const_dble_calc')
xyzzyaadn25=xyzzyaadu1
deallocate(xyzzyaadu1)
allocate(xyzzyaadu1(xyzzyaace25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','const_dble_calc')
xyzzyaadu1(1:xyzzyaabk1(xyzzyaaal25))=xyzzyaadn25(1:xyzzyaabk1(xyzzyaa&
&al25))
xyzzyaadu1(xyzzyaabk1(xyzzyaaal25)+1:xyzzyaabk1(xyzzyaaal25)+xyzzyaaaz&
&25)=xyzzyaabe25(:)
xyzzyaadu1(xyzzyaabk1(xyzzyaaal25)+xyzzyaaaz25:xyzzyaace25)=xyzzyaadn2&
&5(xyzzyaabk1(xyzzyaaal25)+xyzzyaacm1(xyzzyaaal25):xyzzyaacf25)
deallocate(xyzzyaadn25)
xyzzyaabk1(xyzzyaaal25+1:)=xyzzyaabk1(xyzzyaaal25+1:)+xyzzyaaaz25-xyzz&
&yaacm1(xyzzyaaal25)
xyzzyaacm1(xyzzyaaal25)=xyzzyaaaz25
endif
deallocate(xyzzyaabd25)
endif
return
enddo
endif
xyzzyaaba25=0
if(xyzzyaaau25>0)xyzzyaaba25=maxval(groups)
xyzzyaabp25=xyzzyaaaw1
xyzzyaach25=nfn_eebasis
if(xyzzyaabp25==0)then
xyzzyaabr25=0
xyzzyaabt25=0
xyzzyaabv25=0
xyzzyaacb25=0
xyzzyaabx25=0
xyzzyaabz25=0
xyzzyaacd25=0
xyzzyaacf25=0
else
xyzzyaabr25=sum(xyzzyaabw1)
xyzzyaabt25=sum(xyzzyaace1)
xyzzyaabv25=sum(xyzzyaabw1(:)*xyzzyaace1(:))
xyzzyaacb25=sum(xyzzyaabw1(:)*xyzzyaack1(:))
xyzzyaabx25=sum(xyzzyaaci1)
xyzzyaabz25=sum(xyzzyaacg1)
xyzzyaacd25=sum(xyzzyaaco1)
xyzzyaacf25=sum(xyzzyaacm1)
endif
xyzzyaabo25=xyzzyaabp25+1
xyzzyaacg25=xyzzyaach25+iorder
xyzzyaabq25=xyzzyaabr25+xyzzyaaba25
xyzzyaabs25=xyzzyaabt25+xyzzyaaau25
xyzzyaabu25=xyzzyaabv25+xyzzyaaau25*xyzzyaaba25
xyzzyaaca25=xyzzyaacb25+xyzzyaaax25*xyzzyaaba25
xyzzyaabw25=xyzzyaabx25+xyzzyaaav25
xyzzyaaby25=xyzzyaabz25+xyzzyaaaw25
xyzzyaacc25=xyzzyaacd25+xyzzyaaay25
xyzzyaace25=xyzzyaacf25+xyzzyaaaz25
iref=xyzzyaabo25
if(xyzzyaabp25>0)then
allocate(xyzzyaact25(xyzzyaabp25),xyzzyaacw25(xyzzyaabp25),xyzzyaacy25&
&(xyzzyaabp25),xyzzyaacx25(xyzzyaabp25),xyzzyaacz25(xyzzyaabp25),xyzzy&
&aaci25(xyzzyaabp25),xyzzyaacj25(xyzzyaabp25),xyzzyaack25(nspin,nspin,&
&xyzzyaabp25),xyzzyaacl25(xyzzyaabp25),xyzzyaacm25(xyzzyaabp25),xyzzya&
&acn25(xyzzyaabp25),xyzzyaaco25(xyzzyaabp25),xyzzyaacq25(xyzzyaabp25),&
&xyzzyaacp25(xyzzyaabp25),xyzzyaacr25(xyzzyaabp25),xyzzyaada25(xyzzyaa&
&bp25),xyzzyaadb25(xyzzyaabp25),xyzzyaadh25(xyzzyaabp25),stat=xyzzyaaa&
&m25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp of size nset')
xyzzyaact25=xyzzyaaay1
xyzzyaacw25=xyzzyaabc1
xyzzyaacy25=xyzzyaabe1
xyzzyaacx25=xyzzyaabi1
xyzzyaacz25=xyzzyaabk1
xyzzyaaci25=xyzzyaabq1
xyzzyaacj25=xyzzyaabs1
xyzzyaack25=xyzzyaabu1
xyzzyaacl25=xyzzyaabw1
xyzzyaacm25=xyzzyaace1
xyzzyaacn25=xyzzyaack1
xyzzyaaco25=xyzzyaaci1
xyzzyaacq25=xyzzyaacg1
xyzzyaacp25=xyzzyaaco1
xyzzyaacr25=xyzzyaacm1
xyzzyaada25=ifn1_eebasis
xyzzyaadb25=xyzzyaaby1
xyzzyaadh25=xyzzyaaca1
deallocate(xyzzyaaay1,xyzzyaabc1,xyzzyaabe1,xyzzyaabi1,xyzzyaabk1,xyzz&
&yaabq1,  xyzzyaabs1,xyzzyaabu1,xyzzyaabw1,xyzzyaace1,xyzzyaack1,xyzzy&
&aaci1,xyzzyaacg1,xyzzyaaco1,xyzzyaacm1,ifn1_eebasis,xyzzyaaby1,xyzzya&
&aca1)
endif
allocate(xyzzyaaay1(xyzzyaabo25),xyzzyaabc1(xyzzyaabo25),xyzzyaabe1(xy&
&zzyaabo25),xyzzyaabi1(xyzzyaabo25),xyzzyaabk1(xyzzyaabo25),xyzzyaabq1&
&(xyzzyaabo25),xyzzyaabs1(xyzzyaabo25),xyzzyaabu1(nspin,nspin,xyzzyaab&
&o25),xyzzyaabw1(xyzzyaabo25),xyzzyaace1(xyzzyaabo25),xyzzyaack1(xyzzy&
&aabo25),xyzzyaaci1(xyzzyaabo25),xyzzyaacg1(xyzzyaabo25),xyzzyaaco1(xy&
&zzyaabo25),xyzzyaacm1(xyzzyaabo25),ifn1_eebasis(xyzzyaabo25),xyzzyaab&
&y1(xyzzyaabo25),xyzzyaaca1(xyzzyaabo25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','data of size nset')
if(xyzzyaabp25>0)then
xyzzyaaay1(1:xyzzyaabp25)=xyzzyaact25
xyzzyaabc1(1:xyzzyaabp25)=xyzzyaacw25
xyzzyaabe1(1:xyzzyaabp25)=xyzzyaacy25
xyzzyaabi1(1:xyzzyaabp25)=xyzzyaacx25
xyzzyaabk1(1:xyzzyaabp25)=xyzzyaacz25
xyzzyaabq1(1:xyzzyaabp25)=xyzzyaaci25
xyzzyaabs1(1:xyzzyaabp25)=xyzzyaacj25
xyzzyaabu1(:,:,1:xyzzyaabp25)=xyzzyaack25
xyzzyaabw1(1:xyzzyaabp25)=xyzzyaacl25
xyzzyaace1(1:xyzzyaabp25)=xyzzyaacm25
xyzzyaack1(1:xyzzyaabp25)=xyzzyaacn25
xyzzyaaci1(1:xyzzyaabp25)=xyzzyaaco25
xyzzyaacg1(1:xyzzyaabp25)=xyzzyaacq25
xyzzyaaco1(1:xyzzyaabp25)=xyzzyaacp25
xyzzyaacm1(1:xyzzyaabp25)=xyzzyaacr25
ifn1_eebasis(1:xyzzyaabp25)=xyzzyaada25
xyzzyaaby1(1:xyzzyaabp25)=xyzzyaadb25
xyzzyaaca1(1:xyzzyaabp25)=xyzzyaadh25
deallocate(xyzzyaact25,xyzzyaacw25,xyzzyaacy25,xyzzyaacx25,xyzzyaacz25&
&,xyzzyaaci25,xyzzyaacj25,xyzzyaack25,xyzzyaacm25,xyzzyaacn25,xyzzyaac&
&o25,xyzzyaacq25,xyzzyaacp25,xyzzyaacr25,xyzzyaada25,xyzzyaadb25,xyzzy&
&aadh25)
endif
if(xyzzyaaau25>0)then
if(xyzzyaabt25>0)then
allocate(xyzzyaacu25(xyzzyaabr25),xyzzyaadr25(xyzzyaabr25),xyzzyaadq25&
&(xyzzyaabt25),xyzzyaacs25(xyzzyaabv25),xyzzyaade25(xyzzyaabv25),xyzzy&
&aadf25(xyzzyaabv25),xyzzyaadg25(xyzzyaabv25),xyzzyaadi25(xyzzyaabv25)&
&,xyzzyaadk25(xyzzyaabv25),xyzzyaadl25(xyzzyaabv25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp of size nparam')
xyzzyaacu25=xyzzyaaba1
xyzzyaadr25=xyzzyaacq1
xyzzyaadq25=xyzzyaacs1
xyzzyaacs25=xyzzyaacu1
xyzzyaade25=xyzzyaadc1
xyzzyaadf25=xyzzyaade1
xyzzyaadg25=xyzzyaadg1
xyzzyaadi25=xyzzyaacw1
xyzzyaadk25=xyzzyaacy1
xyzzyaadl25=xyzzyaada1
deallocate(xyzzyaaba1,xyzzyaacq1,xyzzyaacs1,xyzzyaacu1,xyzzyaadc1,xyzz&
&yaade1,xyzzyaadg1,xyzzyaacw1,xyzzyaacy1,xyzzyaada1)
endif
allocate(xyzzyaaba1(xyzzyaabq25),xyzzyaacq1(xyzzyaabq25),xyzzyaacs1(xy&
&zzyaabs25),xyzzyaacu1(xyzzyaabu25),xyzzyaadc1(xyzzyaabu25),xyzzyaade1&
&(xyzzyaabu25),xyzzyaadg1(xyzzyaabu25),xyzzyaacw1(xyzzyaabu25),xyzzyaa&
&cy1(xyzzyaabu25),xyzzyaada1(xyzzyaabu25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','data of size nparam')
if(xyzzyaabt25>0)then
xyzzyaaba1(1:xyzzyaabr25)=xyzzyaacu25
xyzzyaacq1(1:xyzzyaabr25)=xyzzyaadr25
xyzzyaacs1(1:xyzzyaabt25)=xyzzyaadq25
xyzzyaacu1(1:xyzzyaabv25)=xyzzyaacs25
xyzzyaadc1(1:xyzzyaabv25)=xyzzyaade25
xyzzyaade1(1:xyzzyaabv25)=xyzzyaadf25
xyzzyaadg1(1:xyzzyaabv25)=xyzzyaadg25
xyzzyaacw1(1:xyzzyaabv25)=xyzzyaadi25
xyzzyaacy1(1:xyzzyaabv25)=xyzzyaadk25
xyzzyaada1(1:xyzzyaabv25)=xyzzyaadl25
deallocate(xyzzyaacu25,xyzzyaadr25,xyzzyaadq25,xyzzyaacs25,xyzzyaade25&
&,xyzzyaadf25,xyzzyaadg25,xyzzyaadi25,xyzzyaadk25,xyzzyaadl25)
endif
endif
if(xyzzyaaax25>0)then
if(xyzzyaacb25>0)then
allocate(xyzzyaacv25(xyzzyaabr25),xyzzyaadj25(xyzzyaacb25),stat=xyzzya&
&aam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp of size nparam_calc')
xyzzyaacv25=xyzzyaabg1
xyzzyaadj25=xyzzyaadi1
endif
if(allocated(xyzzyaadi1))deallocate(xyzzyaabg1,xyzzyaadi1)
allocate(xyzzyaabg1(xyzzyaabq25),xyzzyaadi1(xyzzyaaca25),stat=xyzzyaaa&
&m25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','data of size nparam_calc'&
&)
if(xyzzyaacb25>0)then
xyzzyaabg1(1:xyzzyaabr25)=xyzzyaacv25
xyzzyaadi1(1:xyzzyaacb25)=xyzzyaadj25
deallocate(xyzzyaacv25,xyzzyaadj25)
endif
endif
if(.not.allocated(xyzzyaadi1))then
allocate(xyzzyaabg1(1),xyzzyaadi1(1),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','param_calc')
xyzzyaabg1=0
xyzzyaadi1=0.d0
endif
if(xyzzyaaav25>0)then
if(xyzzyaabx25>0)then
allocate(xyzzyaadc25(xyzzyaabx25),xyzzyaado25(xyzzyaabx25),stat=xyzzya&
&aam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp of size nconst_int')
xyzzyaadc25=xyzzyaadk1
xyzzyaado25=xyzzyaado1
deallocate(xyzzyaadk1,xyzzyaado1)
endif
allocate(xyzzyaadk1(xyzzyaabw25),xyzzyaado1(xyzzyaabw25),stat=xyzzyaaa&
&m25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','data of size nconst_int')
if(xyzzyaabx25>0)then
xyzzyaadk1(1:xyzzyaabx25)=xyzzyaadc25
xyzzyaado1(1:xyzzyaabx25)=xyzzyaado25
deallocate(xyzzyaadc25,xyzzyaado25)
endif
endif
if(xyzzyaaaw25>0)then
if(xyzzyaabz25>0)then
allocate(xyzzyaadm25(xyzzyaabz25),xyzzyaadp25(xyzzyaabz25),stat=xyzzya&
&aam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp of size nconst_dble')
xyzzyaadm25=xyzzyaadm1
xyzzyaadp25=xyzzyaadq1
deallocate(xyzzyaadm1,xyzzyaadq1)
endif
allocate(xyzzyaadm1(xyzzyaaby25),xyzzyaadq1(xyzzyaaby25),stat=xyzzyaaa&
&m25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','data of size nconst_dble'&
&)
if(xyzzyaabz25>0)then
xyzzyaadm1(1:xyzzyaabz25)=xyzzyaadm25
xyzzyaadq1(1:xyzzyaabz25)=xyzzyaadp25
deallocate(xyzzyaadm25,xyzzyaadp25)
endif
endif
if(xyzzyaaay25>0)then
if(xyzzyaacd25>0)then
allocate(xyzzyaadd25(xyzzyaacd25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp of size nconst_int_ca&
&lc')
xyzzyaadd25=xyzzyaads1
endif
if(allocated(xyzzyaads1))deallocate(xyzzyaads1)
allocate(xyzzyaads1(xyzzyaacc25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','data of size nconst_int_c&
&alc')
if(xyzzyaacd25>0)then
xyzzyaads1(1:xyzzyaacd25)=xyzzyaadd25
deallocate(xyzzyaadd25)
endif
endif
if(.not.allocated(xyzzyaads1))then
allocate(xyzzyaads1(1),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','const_int_calc')
xyzzyaads1=0
endif
if(xyzzyaaaz25>0)then
if(xyzzyaacf25>0)then
allocate(xyzzyaadn25(xyzzyaacf25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','tmp of size nconst_dble_c&
&alc')
xyzzyaadn25=xyzzyaadu1
endif
if(allocated(xyzzyaadu1))deallocate(xyzzyaadu1)
allocate(xyzzyaadu1(xyzzyaace25),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','data of size nconst_dble_&
&calc')
if(xyzzyaacf25>0)then
xyzzyaadu1(1:xyzzyaacf25)=xyzzyaadn25
deallocate(xyzzyaadn25)
endif
endif
if(.not.allocated(xyzzyaadu1))then
allocate(xyzzyaadu1(1),stat=xyzzyaaam25)
call check_alloc(xyzzyaaam25,'READ_EEBASIS','const_dble_calc')
xyzzyaadu1=0.d0
endif
xyzzyaaaw1=xyzzyaabo25
nfn_eebasis=xyzzyaacg25
xyzzyaaay1(xyzzyaabo25)=xyzzyaabr25
xyzzyaabq1(xyzzyaabo25)=xyzzyaaat25
xyzzyaabs1(xyzzyaabo25)=iorder
xyzzyaabu1(:,:,xyzzyaabo25)=groups
xyzzyaabw1(xyzzyaabo25)=xyzzyaaba25
xyzzyaace1(xyzzyaabo25)=xyzzyaaau25
xyzzyaack1(xyzzyaabo25)=xyzzyaaax25
xyzzyaaci1(xyzzyaabo25)=xyzzyaaav25
xyzzyaacg1(xyzzyaabo25)=xyzzyaaaw25
xyzzyaaco1(xyzzyaabo25)=xyzzyaaay25
xyzzyaacm1(xyzzyaabo25)=xyzzyaaaz25
ifn1_eebasis(xyzzyaabo25)=xyzzyaach25+1
xyzzyaaby1(xyzzyaabo25)=which_unity
xyzzyaaca1(xyzzyaabo25)=is_isotropic
if(xyzzyaaav25>0)then
xyzzyaadk1(xyzzyaabx25+1:xyzzyaabw25)=xyzzyaabb25
xyzzyaado1(xyzzyaabx25+1:xyzzyaabw25)=xyzzyaabl25
deallocate(xyzzyaabb25,xyzzyaabl25)
xyzzyaabc1(xyzzyaabo25)=xyzzyaabx25
else
xyzzyaabc1(xyzzyaabo25)=0
endif
if(xyzzyaaaw25>0)then
xyzzyaadm1(xyzzyaabz25+1:xyzzyaaby25)=xyzzyaabd25
xyzzyaadq1(xyzzyaabz25+1:xyzzyaaby25)=xyzzyaabm25
deallocate(xyzzyaabd25,xyzzyaabm25)
xyzzyaabe1(xyzzyaabo25)=xyzzyaabz25
else
xyzzyaabe1(xyzzyaabo25)=0
endif
if(xyzzyaaay25>0)then
xyzzyaads1(xyzzyaacd25+1:xyzzyaacc25)=xyzzyaabc25
deallocate(xyzzyaabc25)
xyzzyaabi1(xyzzyaabo25)=xyzzyaacd25
else
xyzzyaabi1(xyzzyaabo25)=0
endif
if(xyzzyaaaz25>0)then
xyzzyaadu1(xyzzyaacf25+1:xyzzyaace25)=xyzzyaabe25
deallocate(xyzzyaabe25)
xyzzyaabk1(xyzzyaabo25)=xyzzyaacf25
else
xyzzyaabk1(xyzzyaabo25)=0
endif
if(xyzzyaaau25>0)then
xyzzyaacs1(xyzzyaabt25+1:xyzzyaabt25+xyzzyaaau25)=xyzzyaabn25(1:xyzzya&
&aau25)
xyzzyaaad25=0
xyzzyaaae25=xyzzyaabv25
xyzzyaaaf25=0
xyzzyaaag25=xyzzyaacb25
do xyzzyaaan25=1,xyzzyaaba25
xyzzyaaba1(xyzzyaabr25+xyzzyaaan25)=xyzzyaaae25
xyzzyaabg1(xyzzyaabr25+xyzzyaaan25)=xyzzyaaag25
xyzzyaaad25=xyzzyaaae25+1
xyzzyaaae25=xyzzyaaae25+xyzzyaaau25
xyzzyaaaf25=xyzzyaaag25+1
xyzzyaaag25=xyzzyaaag25+xyzzyaaax25
xyzzyaacu1(xyzzyaaad25:xyzzyaaae25)=pflag_unset
xyzzyaadc1(xyzzyaaad25:xyzzyaaae25)=xyzzyaabi25
xyzzyaade1(xyzzyaaad25:xyzzyaaae25)=xyzzyaabj25
xyzzyaadg1(xyzzyaaad25:xyzzyaaae25)=xyzzyaabk25
xyzzyaacw1(xyzzyaaad25:xyzzyaaae25)=xyzzyaabf25
xyzzyaacy1(xyzzyaaad25:xyzzyaaae25)=xyzzyaabg25
xyzzyaada1(xyzzyaaad25:xyzzyaaae25)=xyzzyaabh25
enddo
deallocate(xyzzyaabn25,xyzzyaabf25,xyzzyaabi25,xyzzyaabg25,xyzzyaabj25&
&,xyzzyaabh25,xyzzyaabk25)
endif
if(xyzzyaaau25>0)then
xyzzyaaah25=xyzzyaaay1(xyzzyaabo25)
xyzzyaacq1(xyzzyaaah25+1:xyzzyaaah25+xyzzyaaba25)=''
do
call first_unread_child(trim(label)//':Parameters',chname,xyzzyaaak25,&
&flag_as_read=.true.)
if(xyzzyaaak25/=0)exit
if(chname(1:7)/='channel')call errstop_master('READ_EEBASIS','Misnamed&
& item under "'//trim(label)//':Parameters": expected "Channel <model>&
&".')
chname=chname(8:)
call parse_model(chname,2,0,xyzzyaaao25,xyzzyaaap25,xyzzyaaak25)
if(xyzzyaaak25/=0)call errstop_master('READ_EEBASIS','Could not parse &
&model in "'//trim(label)//':Parameters:Channel '//trim(chname)//'".')
xyzzyaaan25=groups(xyzzyaaao25(1),xyzzyaaao25(2))
if(xyzzyaaan25==0)call errstop_master('READ_EEBASIS','Channel "'//trim&
&(chname)//'" does not exist in this system.')
if(xyzzyaaan25<0)call errstop_master('READ_EEBASIS','Channel "'//trim(&
&chname)//'" has been removed by user-defined rules.')
if(trim(xyzzyaacq1(xyzzyaaah25+xyzzyaaan25))/='')call errstop_master('&
&READ_EEBASIS','In "'//trim(label)//':Parameters": channel "'//trim(xy&
&zzyaacq1(xyzzyaaah25+xyzzyaaan25))//'" and channel "'//trim(chname)//&
&'" are redundant, since they correspond to the same particle groups a&
&ccording to the provided rules.')
xyzzyaacq1(xyzzyaaah25+xyzzyaaan25)=chname
enddo
do xyzzyaaan25=1,xyzzyaaba25
if(trim(xyzzyaacq1(xyzzyaaah25+xyzzyaaan25))/='')cycle
xyzzyaaao25=minloc(groups,groups==xyzzyaaan25)
xyzzyaacq1(xyzzyaaah25+xyzzyaaan25)=model_string(2,0,(/minval(xyzzyaaa&
&o25),maxval(xyzzyaaao25)/))
enddo
xyzzyaaab25=sum(xyzzyaace1(1:xyzzyaabo25-1))+1
xyzzyaaac25=xyzzyaaab25+xyzzyaaau25-1
do xyzzyaaan25=1,xyzzyaaba25
xyzzyaaad25=xyzzyaaba1(xyzzyaaah25+xyzzyaaan25)+1
xyzzyaaae25=xyzzyaaad25+xyzzyaaau25-1
call read_gparam_channel(trim(label)//':Parameters:Channel '//trim(xyz&
&zyaacq1(xyzzyaaah25+xyzzyaaan25)),xyzzyaace1(xyzzyaabo25),xyzzyaacs1(&
&xyzzyaaab25:xyzzyaaac25),xyzzyaacw1(xyzzyaaad25:xyzzyaaae25),xyzzyaac&
&u1(xyzzyaaad25:xyzzyaaae25),xyzzyaadc1(xyzzyaaad25:xyzzyaaae25),xyzzy&
&aacy1(xyzzyaaad25:xyzzyaaae25),xyzzyaade1(xyzzyaaad25:xyzzyaaae25),xy&
&zzyaada1(xyzzyaaad25:xyzzyaaae25),xyzzyaadg1(xyzzyaaad25:xyzzyaaae25)&
&)
xyzzyaaaf25=xyzzyaabg1(xyzzyaaah25+xyzzyaaan25)+1
xyzzyaaai25=xyzzyaabi1(xyzzyaabo25)+1
xyzzyaaaj25=xyzzyaabk1(xyzzyaabo25)+1
call xyzzyaafc1(xyzzyaabq1(xyzzyaabo25),xyzzyaabs1(xyzzyaabo25),xyzzya&
&ace1(xyzzyaabo25),xyzzyaack1(xyzzyaabo25),xyzzyaacw1(xyzzyaaad25),xyz&
&zyaaco1(xyzzyaabo25),xyzzyaads1(xyzzyaaai25),xyzzyaacm1(xyzzyaabo25),&
&xyzzyaadu1(xyzzyaaaj25),xyzzyaadi1(xyzzyaaaf25))
enddo
endif
end subroutine read_eebasis
subroutine read_enbasis(label,namespace,groups,iorder,iref,is_isotropi&
&c,which_unity,require_cutoff,is_cusp_friendly,f_zero,is_none,basis_sy&
&mmetry)
implicit none
integer,intent(in) :: groups(nitot,nspin),iorder
integer,intent(out) :: iref,which_unity,basis_symmetry
character(8),intent(out) :: f_zero
logical,intent(out) :: is_isotropic,require_cutoff,is_cusp_friendly,is&
&_none
character(*),intent(in) :: label,namespace
integer xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26,xyzzyaaad26,xyzzyaaae26,xy&
&zzyaaaf26,xyzzyaaag26,xyzzyaaah26,xyzzyaaai26,xyzzyaaaj26,xyzzyaaak26&
&,xyzzyaaal26,xyzzyaaam26,xyzzyaaan26,xyzzyaaao26(1),xyzzyaaap26(1),xy&
&zzyaaaq26(2),xyzzyaaar26
logical xyzzyaaas26,xyzzyaaat26
character(casl_fullkeysize) btype
character(casl_keysize) chname
integer xyzzyaaau26,xyzzyaaav26,xyzzyaaaw26,xyzzyaaax26,xyzzyaaay26,xy&
&zzyaaaz26,xyzzyaaba26,xyzzyaabb26
integer,pointer :: xyzzyaabc26(:),xyzzyaabd26(:)
real(dp),pointer :: xyzzyaabe26(:),xyzzyaabf26(:),xyzzyaabg26(:),xyzzy&
&aabh26(:),xyzzyaabi26(:)
logical,pointer :: xyzzyaabj26(:),xyzzyaabk26(:),xyzzyaabl26(:)
character(casl_fullkeysize),pointer :: xyzzyaabm26(:),xyzzyaabn26(:)
character(casl_keysize),pointer :: xyzzyaabo26(:)
integer xyzzyaabp26,xyzzyaabq26,xyzzyaabr26,xyzzyaabs26,xyzzyaabt26,xy&
&zzyaabu26,xyzzyaabv26,xyzzyaabw26,xyzzyaabx26,xyzzyaaby26,xyzzyaabz26&
&,xyzzyaaca26,xyzzyaacb26,xyzzyaacc26,xyzzyaacd26,xyzzyaace26,xyzzyaac&
&f26,xyzzyaacg26,xyzzyaach26,xyzzyaaci26
integer,allocatable :: xyzzyaacj26(:),xyzzyaack26(:),xyzzyaacl26(:,:,:&
&),xyzzyaacm26(:),xyzzyaacn26(:),xyzzyaaco26(:),xyzzyaacp26(:),xyzzyaa&
&cq26(:),xyzzyaacr26(:),xyzzyaacs26(:),xyzzyaact26(:),xyzzyaacu26(:),x&
&yzzyaacv26(:),xyzzyaacw26(:),xyzzyaacx26(:),xyzzyaacy26(:),xyzzyaacz2&
&6(:),xyzzyaada26(:),xyzzyaadb26(:),xyzzyaadc26(:),xyzzyaadd26(:),xyzz&
&yaade26(:)
logical,allocatable :: xyzzyaadf26(:),xyzzyaadg26(:),xyzzyaadh26(:),xy&
&zzyaadi26(:)
real(dp),allocatable :: xyzzyaadj26(:),xyzzyaadk26(:),xyzzyaadl26(:),x&
&yzzyaadm26(:),xyzzyaadn26(:),xyzzyaado26(:)
character(casl_fullkeysize),pointer :: xyzzyaadp26(:),xyzzyaadq26(:)
character(casl_keysize),pointer :: xyzzyaadr26(:),xyzzyaads26(:)
call get_casl_item(trim(label)//':Type',btype,xyzzyaaak26)
if(xyzzyaaak26/=0)then
if(iorder>1)call errstop_master('READ_ENBASIS','Problem getting "'//tr&
&im(label)//':Type".')
btype='none'
endif
call xyzzyaaez1(trim(namespace)//':'//trim(btype),xyzzyaaau26,is_isotr&
&opic,xyzzyaaas26,xyzzyaaat26,is_cusp_friendly,f_zero,is_none)
if(xyzzyaaau26<0)call errstop_master('READ_ENBASIS','Unrecognized basi&
&s set "'//trim(btype)//'".')
if(xyzzyaaas26.and..not.isperiodic)call errstop_master('READ_ENBASIS',&
&'Basis set "'//trim(btype)//'" cannot be used in non-periodic systems&
&.')
require_cutoff=xyzzyaaat26.and.isperiodic
call xyzzyaafa1(trim(label)//':Constants',xyzzyaaau26,iorder,xyzzyaaaw&
&26,xyzzyaabm26,xyzzyaabc26,xyzzyaaax26,xyzzyaabn26,xyzzyaabe26,xyzzya&
&aaz26,xyzzyaabd26,xyzzyaaba26,xyzzyaabf26,which_unity,basis_symmetry)
call xyzzyaafb1(xyzzyaaau26,iorder,xyzzyaaaz26,xyzzyaabd26,xyzzyaaba26&
&,xyzzyaabf26,xyzzyaaav26,xyzzyaaay26,xyzzyaabo26,xyzzyaabg26,xyzzyaab&
&j26,xyzzyaabh26,xyzzyaabk26,xyzzyaabi26,xyzzyaabl26)
if(xyzzyaaav26==0)then
do xyzzyaaal26=1,xyzzyaaax1
if(xyzzyaabr1(xyzzyaaal26)/=xyzzyaaau26)cycle
if(xyzzyaaaw26/=0)then
xyzzyaaaa26=min(xyzzyaaaw26,xyzzyaacj1(xyzzyaaal26))
if(any(xyzzyaabc26(1:xyzzyaaaa26)/=xyzzyaadl1(xyzzyaabd1(xyzzyaaal26)+&
&1:xyzzyaabd1(xyzzyaaal26)+xyzzyaaaa26)))cycle
endif
if(xyzzyaaax26/=0)then
xyzzyaaaa26=min(xyzzyaaax26,xyzzyaach1(xyzzyaaal26))
if(any(xyzzyaabe26(1:xyzzyaaaa26)/=xyzzyaadn1(xyzzyaabf1(xyzzyaaal26)+&
&1:xyzzyaabf1(xyzzyaaal26)+xyzzyaaaa26)))cycle
endif
iref=xyzzyaaal26
if(iorder>xyzzyaabt1(xyzzyaaal26))then
nfn_enbasis=nfn_enbasis+iorder-xyzzyaabt1(xyzzyaaal26)
do xyzzyaaar26=xyzzyaaal26+1,xyzzyaaax1
ifn1_enbasis(xyzzyaaar26)=ifn1_enbasis(xyzzyaaar26)+iorder-xyzzyaabt1(&
&xyzzyaaal26)
enddo
xyzzyaabt1(xyzzyaaal26)=iorder
endif
if(xyzzyaaaw26>0)then
if(xyzzyaacj1(xyzzyaaal26)<xyzzyaaaw26)then
xyzzyaaby26=sum(xyzzyaacj1)
xyzzyaabx26=xyzzyaaby26+xyzzyaaaw26-xyzzyaacj1(xyzzyaaal26)
allocate(xyzzyaadd26(xyzzyaaby26),xyzzyaadp26(xyzzyaaby26),stat=xyzzya&
&aam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp_const_int')
xyzzyaadd26=xyzzyaadl1
xyzzyaadp26=xyzzyaadp1
deallocate(xyzzyaadl1,xyzzyaadp1)
allocate(xyzzyaadl1(xyzzyaabx26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','const_int')
xyzzyaadl1(1:xyzzyaabd1(xyzzyaaal26))=xyzzyaadd26(1:xyzzyaabd1(xyzzyaa&
&al26))
xyzzyaadl1(xyzzyaabd1(xyzzyaaal26)+1:xyzzyaabd1(xyzzyaaal26)+xyzzyaaaw&
&26)=xyzzyaabc26(:)
xyzzyaadl1(xyzzyaabd1(xyzzyaaal26)+xyzzyaaaw26:xyzzyaabx26)=xyzzyaadd2&
&6(xyzzyaabd1(xyzzyaaal26)+xyzzyaacj1(xyzzyaaal26):xyzzyaaby26)
xyzzyaadp1(1:xyzzyaabd1(xyzzyaaal26))=xyzzyaadp26(1:xyzzyaabd1(xyzzyaa&
&al26))
xyzzyaadp1(xyzzyaabd1(xyzzyaaal26)+1:xyzzyaabd1(xyzzyaaal26)+xyzzyaaaw&
&26)=xyzzyaabm26(:)
xyzzyaadp1(xyzzyaabd1(xyzzyaaal26)+xyzzyaaaw26:xyzzyaabx26)=xyzzyaadp2&
&6(xyzzyaabd1(xyzzyaaal26)+xyzzyaacj1(xyzzyaaal26):xyzzyaaby26)
deallocate(xyzzyaadd26,xyzzyaadp26)
xyzzyaabd1(xyzzyaaal26+1:)=xyzzyaabd1(xyzzyaaal26+1:)+xyzzyaaaw26-xyzz&
&yaacj1(xyzzyaaal26)
xyzzyaacj1(xyzzyaaal26)=xyzzyaaaw26
endif
deallocate(xyzzyaabc26,xyzzyaabm26)
if(xyzzyaacp1(xyzzyaaal26)<xyzzyaaaz26)then
xyzzyaace26=sum(xyzzyaacp1)
xyzzyaacd26=xyzzyaace26+xyzzyaaaz26-xyzzyaacp1(xyzzyaaal26)
allocate(xyzzyaade26(xyzzyaace26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp_const_int_calc')
xyzzyaade26=xyzzyaadt1
deallocate(xyzzyaadt1)
allocate(xyzzyaadt1(xyzzyaacd26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','const_int_calc')
xyzzyaadt1(1:xyzzyaabj1(xyzzyaaal26))=xyzzyaade26(1:xyzzyaabj1(xyzzyaa&
&al26))
xyzzyaadt1(xyzzyaabj1(xyzzyaaal26)+1:xyzzyaabj1(xyzzyaaal26)+xyzzyaaaz&
&26)=xyzzyaabd26(:)
xyzzyaadt1(xyzzyaabj1(xyzzyaaal26)+xyzzyaaaz26:xyzzyaacd26)=xyzzyaade2&
&6(xyzzyaabj1(xyzzyaaal26)+xyzzyaacp1(xyzzyaaal26):xyzzyaace26)
deallocate(xyzzyaade26)
xyzzyaabj1(xyzzyaaal26+1:)=xyzzyaabj1(xyzzyaaal26+1:)+xyzzyaaaz26-xyzz&
&yaacp1(xyzzyaaal26)
xyzzyaacp1(xyzzyaaal26)=xyzzyaaaz26
endif
deallocate(xyzzyaabd26)
endif
if(xyzzyaaax26>0)then
if(xyzzyaach1(xyzzyaaal26)<xyzzyaaax26)then
xyzzyaaca26=sum(xyzzyaach1)
xyzzyaabz26=xyzzyaaca26+xyzzyaaax26-xyzzyaach1(xyzzyaaal26)
allocate(xyzzyaadn26(xyzzyaaca26),xyzzyaadq26(xyzzyaaca26),stat=xyzzya&
&aam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp_const_dble')
xyzzyaadn26=xyzzyaadn1
xyzzyaadq26=xyzzyaadr1
deallocate(xyzzyaadn1,xyzzyaadr1)
allocate(xyzzyaadn1(xyzzyaabz26),xyzzyaadr1(xyzzyaabz26),stat=xyzzyaaa&
&m26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','const_dble')
xyzzyaadn1(1:xyzzyaabf1(xyzzyaaal26))=xyzzyaadn26(1:xyzzyaabf1(xyzzyaa&
&al26))
xyzzyaadn1(xyzzyaabf1(xyzzyaaal26)+1:xyzzyaabf1(xyzzyaaal26)+xyzzyaaax&
&26)=xyzzyaabe26(:)
xyzzyaadn1(xyzzyaabf1(xyzzyaaal26)+xyzzyaaax26:xyzzyaabz26)=xyzzyaadn2&
&6(xyzzyaabf1(xyzzyaaal26)+xyzzyaach1(xyzzyaaal26):xyzzyaaca26)
xyzzyaadr1(1:xyzzyaabf1(xyzzyaaal26))=xyzzyaadq26(1:xyzzyaabf1(xyzzyaa&
&al26))
xyzzyaadr1(xyzzyaabf1(xyzzyaaal26)+1:xyzzyaabf1(xyzzyaaal26)+xyzzyaaax&
&26)=xyzzyaabn26(:)
xyzzyaadr1(xyzzyaabf1(xyzzyaaal26)+xyzzyaaax26:xyzzyaabz26)=xyzzyaadq2&
&6(xyzzyaabf1(xyzzyaaal26)+xyzzyaach1(xyzzyaaal26):xyzzyaaca26)
deallocate(xyzzyaadn26,xyzzyaadq26)
xyzzyaabf1(xyzzyaaal26+1:)=xyzzyaabf1(xyzzyaaal26+1:)+xyzzyaaax26-xyzz&
&yaach1(xyzzyaaal26)
xyzzyaach1(xyzzyaaal26)=xyzzyaaax26
endif
deallocate(xyzzyaabe26,xyzzyaabn26)
if(xyzzyaacn1(xyzzyaaal26)<xyzzyaaba26)then
xyzzyaacg26=sum(xyzzyaacn1)
xyzzyaacf26=xyzzyaacg26+xyzzyaaba26-xyzzyaacn1(xyzzyaaal26)
allocate(xyzzyaado26(xyzzyaacg26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp_const_dble_calc')
xyzzyaado26=xyzzyaadv1
deallocate(xyzzyaadv1)
allocate(xyzzyaadv1(xyzzyaacf26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','const_dble_calc')
xyzzyaadv1(1:xyzzyaabl1(xyzzyaaal26))=xyzzyaado26(1:xyzzyaabl1(xyzzyaa&
&al26))
xyzzyaadv1(xyzzyaabl1(xyzzyaaal26)+1:xyzzyaabl1(xyzzyaaal26)+xyzzyaaba&
&26)=xyzzyaabf26(:)
xyzzyaadv1(xyzzyaabl1(xyzzyaaal26)+xyzzyaaba26:xyzzyaacf26)=xyzzyaado2&
&6(xyzzyaabl1(xyzzyaaal26)+xyzzyaacn1(xyzzyaaal26):xyzzyaacg26)
deallocate(xyzzyaado26)
xyzzyaabl1(xyzzyaaal26+1:)=xyzzyaabl1(xyzzyaaal26+1:)+xyzzyaaba26-xyzz&
&yaacn1(xyzzyaaal26)
xyzzyaacn1(xyzzyaaal26)=xyzzyaaba26
endif
deallocate(xyzzyaabe26)
endif
return
enddo
endif
xyzzyaabb26=0
if(xyzzyaaav26>0)xyzzyaabb26=maxval(groups)
xyzzyaabq26=xyzzyaaax1
xyzzyaaci26=nfn_enbasis
if(xyzzyaabq26==0)then
xyzzyaabs26=0
xyzzyaabu26=0
xyzzyaabw26=0
xyzzyaacc26=0
xyzzyaaby26=0
xyzzyaaca26=0
xyzzyaace26=0
xyzzyaacg26=0
else
xyzzyaabs26=sum(xyzzyaabx1)
xyzzyaabu26=sum(xyzzyaacf1)
xyzzyaabw26=sum(xyzzyaabx1(:)*xyzzyaacf1(:))
xyzzyaacc26=sum(xyzzyaabx1(:)*xyzzyaacl1(:))
xyzzyaaby26=sum(xyzzyaacj1)
xyzzyaaca26=sum(xyzzyaach1)
xyzzyaace26=sum(xyzzyaacp1)
xyzzyaacg26=sum(xyzzyaacn1)
endif
xyzzyaabp26=xyzzyaabq26+1
xyzzyaach26=xyzzyaaci26+iorder
xyzzyaabr26=xyzzyaabs26+xyzzyaabb26
xyzzyaabt26=xyzzyaabu26+xyzzyaaav26
xyzzyaabv26=xyzzyaabw26+xyzzyaaav26*xyzzyaabb26
xyzzyaacb26=xyzzyaacc26+xyzzyaaay26*xyzzyaabb26
xyzzyaabx26=xyzzyaaby26+xyzzyaaaw26
xyzzyaabz26=xyzzyaaca26+xyzzyaaax26
xyzzyaacd26=xyzzyaace26+xyzzyaaaz26
xyzzyaacf26=xyzzyaacg26+xyzzyaaba26
iref=xyzzyaabp26
if(xyzzyaabq26>0)then
allocate(xyzzyaacu26(xyzzyaabq26),xyzzyaacx26(xyzzyaabq26),xyzzyaacz26&
&(xyzzyaabq26),xyzzyaacy26(xyzzyaabq26),xyzzyaada26(xyzzyaabq26),xyzzy&
&aacj26(xyzzyaabq26),xyzzyaack26(xyzzyaabq26),xyzzyaacl26(nitot,nspin,&
&xyzzyaabq26),xyzzyaacm26(xyzzyaabq26),xyzzyaacn26(xyzzyaabq26),xyzzya&
&aco26(xyzzyaabq26),xyzzyaacp26(xyzzyaabq26),xyzzyaacr26(xyzzyaabq26),&
&xyzzyaacq26(xyzzyaabq26),xyzzyaacs26(xyzzyaabq26),xyzzyaadb26(xyzzyaa&
&bq26),xyzzyaadc26(xyzzyaabq26),xyzzyaadi26(xyzzyaabq26),stat=xyzzyaaa&
&m26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp of size nset')
xyzzyaacu26=xyzzyaaaz1
xyzzyaacx26=xyzzyaabd1
xyzzyaacz26=xyzzyaabf1
xyzzyaacy26=xyzzyaabj1
xyzzyaada26=xyzzyaabl1
xyzzyaacj26=xyzzyaabr1
xyzzyaack26=xyzzyaabt1
xyzzyaacl26=xyzzyaabv1
xyzzyaacm26=xyzzyaabx1
xyzzyaacn26=xyzzyaacf1
xyzzyaaco26=xyzzyaacl1
xyzzyaacp26=xyzzyaacj1
xyzzyaacr26=xyzzyaach1
xyzzyaacq26=xyzzyaacp1
xyzzyaacs26=xyzzyaacn1
xyzzyaadb26=ifn1_enbasis
xyzzyaadc26=xyzzyaabz1
xyzzyaadi26=xyzzyaacb1
deallocate(xyzzyaaaz1,xyzzyaabd1,xyzzyaabf1,xyzzyaabj1,xyzzyaabl1,xyzz&
&yaabr1,xyzzyaabt1,xyzzyaabv1,xyzzyaabx1,xyzzyaacf1,xyzzyaacl1,xyzzyaa&
&cj1,xyzzyaach1,xyzzyaacp1,xyzzyaacn1,ifn1_enbasis,xyzzyaabz1,xyzzyaac&
&b1)
endif
allocate(xyzzyaaaz1(xyzzyaabp26),xyzzyaabd1(xyzzyaabp26),xyzzyaabf1(xy&
&zzyaabp26),xyzzyaabj1(xyzzyaabp26),xyzzyaabl1(xyzzyaabp26),xyzzyaabr1&
&(xyzzyaabp26),xyzzyaabt1(xyzzyaabp26),xyzzyaabv1(nitot,nspin,xyzzyaab&
&p26),xyzzyaabx1(xyzzyaabp26),xyzzyaacf1(xyzzyaabp26),xyzzyaacl1(xyzzy&
&aabp26),xyzzyaacj1(xyzzyaabp26),xyzzyaach1(xyzzyaabp26),xyzzyaacp1(xy&
&zzyaabp26),xyzzyaacn1(xyzzyaabp26),ifn1_enbasis(xyzzyaabp26),xyzzyaab&
&z1(xyzzyaabp26),xyzzyaacb1(xyzzyaabp26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','data of size nset')
if(xyzzyaabq26>0)then
xyzzyaaaz1(1:xyzzyaabq26)=xyzzyaacu26
xyzzyaabd1(1:xyzzyaabq26)=xyzzyaacx26
xyzzyaabf1(1:xyzzyaabq26)=xyzzyaacz26
xyzzyaabj1(1:xyzzyaabq26)=xyzzyaacy26
xyzzyaabl1(1:xyzzyaabq26)=xyzzyaada26
xyzzyaabr1(1:xyzzyaabq26)=xyzzyaacj26
xyzzyaabt1(1:xyzzyaabq26)=xyzzyaack26
xyzzyaabv1(:,:,1:xyzzyaabq26)=xyzzyaacl26
xyzzyaabx1(1:xyzzyaabq26)=xyzzyaacm26
xyzzyaacf1(1:xyzzyaabq26)=xyzzyaacn26
xyzzyaacl1(1:xyzzyaabq26)=xyzzyaaco26
xyzzyaacj1(1:xyzzyaabq26)=xyzzyaacp26
xyzzyaach1(1:xyzzyaabq26)=xyzzyaacr26
xyzzyaacp1(1:xyzzyaabq26)=xyzzyaacq26
xyzzyaacn1(1:xyzzyaabq26)=xyzzyaacs26
ifn1_enbasis(1:xyzzyaabq26)=xyzzyaadb26
xyzzyaabz1(1:xyzzyaabq26)=xyzzyaadc26
xyzzyaacb1(1:xyzzyaabq26)=xyzzyaadi26
deallocate(xyzzyaacu26,xyzzyaacx26,xyzzyaacz26,xyzzyaacy26,xyzzyaada26&
&,xyzzyaacj26,xyzzyaack26,xyzzyaacl26,xyzzyaacn26,xyzzyaaco26,xyzzyaac&
&p26,xyzzyaacr26,xyzzyaacq26,xyzzyaacs26,xyzzyaadb26,xyzzyaadc26,xyzzy&
&aadi26)
endif
if(xyzzyaaav26>0)then
if(xyzzyaabu26>0)then
allocate(xyzzyaacv26(xyzzyaabs26),xyzzyaads26(xyzzyaabs26),xyzzyaadr26&
&(xyzzyaabu26),xyzzyaact26(xyzzyaabw26),xyzzyaadf26(xyzzyaabw26),xyzzy&
&aadg26(xyzzyaabw26),xyzzyaadh26(xyzzyaabw26),xyzzyaadj26(xyzzyaabw26)&
&,xyzzyaadl26(xyzzyaabw26),xyzzyaadm26(xyzzyaabw26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp of size nparam')
xyzzyaacv26=xyzzyaabb1
xyzzyaads26=xyzzyaacr1
xyzzyaadr26=xyzzyaact1
xyzzyaact26=xyzzyaacv1
xyzzyaadf26=xyzzyaadd1
xyzzyaadg26=xyzzyaadf1
xyzzyaadh26=xyzzyaadh1
xyzzyaadj26=xyzzyaacx1
xyzzyaadl26=xyzzyaacz1
xyzzyaadm26=xyzzyaadb1
deallocate(xyzzyaabb1,xyzzyaacr1,xyzzyaact1,xyzzyaacv1,xyzzyaadd1,xyzz&
&yaadf1,xyzzyaadh1,xyzzyaacx1,xyzzyaacz1,xyzzyaadb1)
endif
allocate(xyzzyaabb1(xyzzyaabr26),xyzzyaacr1(xyzzyaabr26),xyzzyaact1(xy&
&zzyaabt26),xyzzyaacv1(xyzzyaabv26),xyzzyaadd1(xyzzyaabv26),xyzzyaadf1&
&(xyzzyaabv26),xyzzyaadh1(xyzzyaabv26),xyzzyaacx1(xyzzyaabv26),xyzzyaa&
&cz1(xyzzyaabv26),xyzzyaadb1(xyzzyaabv26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','data of size nparam')
if(xyzzyaabu26>0)then
xyzzyaabb1(1:xyzzyaabs26)=xyzzyaacv26
xyzzyaacr1(1:xyzzyaabs26)=xyzzyaads26
xyzzyaact1(1:xyzzyaabu26)=xyzzyaadr26
xyzzyaacv1(1:xyzzyaabw26)=xyzzyaact26
xyzzyaadd1(1:xyzzyaabw26)=xyzzyaadf26
xyzzyaadf1(1:xyzzyaabw26)=xyzzyaadg26
xyzzyaadh1(1:xyzzyaabw26)=xyzzyaadh26
xyzzyaacx1(1:xyzzyaabw26)=xyzzyaadj26
xyzzyaacz1(1:xyzzyaabw26)=xyzzyaadl26
xyzzyaadb1(1:xyzzyaabw26)=xyzzyaadm26
deallocate(xyzzyaacv26,xyzzyaads26,xyzzyaadr26,xyzzyaact26,xyzzyaadf26&
&,xyzzyaadg26,xyzzyaadh26,xyzzyaadj26,xyzzyaadl26,xyzzyaadm26)
endif
endif
if(xyzzyaaay26>0)then
if(xyzzyaacc26>0)then
allocate(xyzzyaacw26(xyzzyaabs26),xyzzyaadk26(xyzzyaacc26),stat=xyzzya&
&aam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp of size nparam_calc')
xyzzyaacw26=xyzzyaabh1
xyzzyaadk26=xyzzyaadj1
endif
if(allocated(xyzzyaadj1))deallocate(xyzzyaabh1,xyzzyaadj1)
allocate(xyzzyaabh1(xyzzyaabr26),xyzzyaadj1(xyzzyaacb26),stat=xyzzyaaa&
&m26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','data of size nparam_calc'&
&)
if(xyzzyaacc26>0)then
xyzzyaabh1(1:xyzzyaabs26)=xyzzyaacw26
xyzzyaadj1(1:xyzzyaacc26)=xyzzyaadk26
deallocate(xyzzyaacw26,xyzzyaadk26)
endif
endif
if(.not.allocated(xyzzyaadj1))then
allocate(xyzzyaabh1(1),xyzzyaadj1(1),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','param_calc')
xyzzyaabh1=0
xyzzyaadj1=0.d0
endif
if(xyzzyaaaw26>0)then
if(xyzzyaaby26>0)then
allocate(xyzzyaadd26(xyzzyaaby26),xyzzyaadp26(xyzzyaaby26),stat=xyzzya&
&aam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp of size nconst_int')
xyzzyaadd26=xyzzyaadl1
xyzzyaadp26=xyzzyaadp1
deallocate(xyzzyaadl1,xyzzyaadp1)
endif
allocate(xyzzyaadl1(xyzzyaabx26),xyzzyaadp1(xyzzyaabx26),stat=xyzzyaaa&
&m26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','data of size nconst_int')
if(xyzzyaaby26>0)then
xyzzyaadl1(1:xyzzyaaby26)=xyzzyaadd26
xyzzyaadp1(1:xyzzyaaby26)=xyzzyaadp26
deallocate(xyzzyaadd26,xyzzyaadp26)
endif
endif
if(xyzzyaaax26>0)then
if(xyzzyaaca26>0)then
allocate(xyzzyaadn26(xyzzyaaca26),xyzzyaadq26(xyzzyaaca26),stat=xyzzya&
&aam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp of size nconst_dble')
xyzzyaadn26=xyzzyaadn1
xyzzyaadq26=xyzzyaadr1
deallocate(xyzzyaadn1,xyzzyaadr1)
endif
allocate(xyzzyaadn1(xyzzyaabz26),xyzzyaadr1(xyzzyaabz26),stat=xyzzyaaa&
&m26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','data of size nconst_dble'&
&)
if(xyzzyaaca26>0)then
xyzzyaadn1(1:xyzzyaaca26)=xyzzyaadn26
xyzzyaadr1(1:xyzzyaaca26)=xyzzyaadq26
deallocate(xyzzyaadn26,xyzzyaadq26)
endif
endif
if(xyzzyaaaz26>0)then
if(xyzzyaace26>0)then
allocate(xyzzyaade26(xyzzyaace26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp of size nconst_int_ca&
&lc')
xyzzyaade26=xyzzyaadt1
endif
if(allocated(xyzzyaadt1))deallocate(xyzzyaadt1)
allocate(xyzzyaadt1(xyzzyaacd26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','data of size nconst_int_c&
&alc')
if(xyzzyaace26>0)then
xyzzyaadt1(1:xyzzyaace26)=xyzzyaade26
deallocate(xyzzyaade26)
endif
endif
if(.not.allocated(xyzzyaadt1))then
allocate(xyzzyaadt1(1),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','const_int_calc')
xyzzyaadt1=0
endif
if(xyzzyaaba26>0)then
if(xyzzyaacg26>0)then
allocate(xyzzyaado26(xyzzyaacg26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','tmp of size nconst_dble_c&
&alc')
xyzzyaado26=xyzzyaadv1
deallocate(xyzzyaadv1)
endif
if(allocated(xyzzyaadv1))deallocate(xyzzyaadv1)
allocate(xyzzyaadv1(xyzzyaacf26),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','data of size nconst_dble_&
&calc')
if(xyzzyaacg26>0)then
xyzzyaadv1(1:xyzzyaacg26)=xyzzyaado26
deallocate(xyzzyaado26)
endif
endif
if(.not.allocated(xyzzyaadv1))then
allocate(xyzzyaadv1(1),stat=xyzzyaaam26)
call check_alloc(xyzzyaaam26,'READ_ENBASIS','const_dble_calc')
xyzzyaadv1=0.d0
endif
xyzzyaaax1=xyzzyaabp26
nfn_enbasis=xyzzyaach26
xyzzyaaaz1(xyzzyaabp26)=xyzzyaabs26
xyzzyaabr1(xyzzyaabp26)=xyzzyaaau26
xyzzyaabt1(xyzzyaabp26)=iorder
xyzzyaabv1(:,:,xyzzyaabp26)=groups
xyzzyaabx1(xyzzyaabp26)=xyzzyaabb26
xyzzyaacf1(xyzzyaabp26)=xyzzyaaav26
xyzzyaacl1(xyzzyaabp26)=xyzzyaaay26
xyzzyaacj1(xyzzyaabp26)=xyzzyaaaw26
xyzzyaach1(xyzzyaabp26)=xyzzyaaax26
xyzzyaacp1(xyzzyaabp26)=xyzzyaaaz26
xyzzyaacn1(xyzzyaabp26)=xyzzyaaba26
ifn1_enbasis(xyzzyaabp26)=xyzzyaaci26+1
xyzzyaabz1(xyzzyaabp26)=which_unity
xyzzyaacb1(xyzzyaabp26)=is_isotropic
if(xyzzyaaaw26>0)then
xyzzyaadl1(xyzzyaaby26+1:xyzzyaabx26)=xyzzyaabc26
xyzzyaadp1(xyzzyaaby26+1:xyzzyaabx26)=xyzzyaabm26
deallocate(xyzzyaabc26,xyzzyaabm26)
xyzzyaabd1(xyzzyaabp26)=xyzzyaaby26
else
xyzzyaabd1(xyzzyaabp26)=0
endif
if(xyzzyaaax26>0)then
xyzzyaadn1(xyzzyaaca26+1:xyzzyaabz26)=xyzzyaabe26
xyzzyaadr1(xyzzyaaca26+1:xyzzyaabz26)=xyzzyaabn26
deallocate(xyzzyaabe26,xyzzyaabn26)
xyzzyaabf1(xyzzyaabp26)=xyzzyaaca26
else
xyzzyaabf1(xyzzyaabp26)=0
endif
if(xyzzyaaaz26>0)then
xyzzyaadt1(xyzzyaace26+1:xyzzyaacd26)=xyzzyaabd26
deallocate(xyzzyaabd26)
xyzzyaabj1(xyzzyaabp26)=xyzzyaace26
else
xyzzyaabj1(xyzzyaabp26)=0
endif
if(xyzzyaaba26>0)then
xyzzyaadv1(xyzzyaacg26+1:xyzzyaacf26)=xyzzyaabf26
deallocate(xyzzyaabf26)
xyzzyaabl1(xyzzyaabp26)=xyzzyaacg26
else
xyzzyaabl1(xyzzyaabp26)=0
endif
if(xyzzyaaav26>0)then
xyzzyaact1(xyzzyaabu26+1:xyzzyaabu26+xyzzyaaav26)=xyzzyaabo26(1:xyzzya&
&aav26)
xyzzyaaad26=0
xyzzyaaae26=xyzzyaabw26
xyzzyaaaf26=0
xyzzyaaag26=xyzzyaacc26
do xyzzyaaan26=1,xyzzyaabb26
xyzzyaabb1(xyzzyaabs26+xyzzyaaan26)=xyzzyaaae26
xyzzyaabh1(xyzzyaabs26+xyzzyaaan26)=xyzzyaaag26
xyzzyaaad26=xyzzyaaae26+1
xyzzyaaae26=xyzzyaaae26+xyzzyaaav26
xyzzyaaaf26=xyzzyaaag26+1
xyzzyaaag26=xyzzyaaag26+xyzzyaaay26
xyzzyaacv1(xyzzyaaad26:xyzzyaaae26)=pflag_unset
xyzzyaadd1(xyzzyaaad26:xyzzyaaae26)=xyzzyaabj26
xyzzyaadf1(xyzzyaaad26:xyzzyaaae26)=xyzzyaabk26
xyzzyaadh1(xyzzyaaad26:xyzzyaaae26)=xyzzyaabl26
xyzzyaacx1(xyzzyaaad26:xyzzyaaae26)=xyzzyaabg26
xyzzyaacz1(xyzzyaaad26:xyzzyaaae26)=xyzzyaabh26
xyzzyaadb1(xyzzyaaad26:xyzzyaaae26)=xyzzyaabi26
enddo
deallocate(xyzzyaabo26,xyzzyaabg26,xyzzyaabj26,xyzzyaabh26,xyzzyaabk26&
&,xyzzyaabi26,xyzzyaabl26)
endif
if(xyzzyaaav26>0)then
xyzzyaaah26=xyzzyaaaz1(xyzzyaabp26)
xyzzyaacr1(xyzzyaaah26+1:xyzzyaaah26+xyzzyaabb26)=''
do
call first_unread_child(trim(label)//':Parameters',chname,xyzzyaaak26,&
&flag_as_read=.true.)
if(xyzzyaaak26/=0)exit
if(chname(1:7)/='channel')call errstop_master('READ_ENBASIS','Misnamed&
& item under "'//trim(label)//':Parameters": expected "Channel <model>&
&".')
chname=chname(8:)
call parse_model(chname,1,1,xyzzyaaao26,xyzzyaaap26,xyzzyaaak26)
if(xyzzyaaak26/=0)call errstop_master('READ_ENBASIS','Could not parse &
&model in "'//trim(label)//':Parameters:Channel '//trim(chname)//'".')
xyzzyaaan26=groups(xyzzyaaap26(1),xyzzyaaao26(1))
if(xyzzyaaan26==0)call errstop_master('READ_ENBASIS','Channel "'//trim&
&(chname)//'" does not exist in this system.')
if(xyzzyaaan26<0)call errstop_master('READ_ENBASIS','Channel "'//trim(&
&chname)//'" has been removed by user-defined rules.')
if(trim(xyzzyaacr1(xyzzyaaah26+xyzzyaaan26))/='')call errstop_master('&
&READ_ENBASIS','In "'//trim(label)//':Parameters": channel "'//trim(xy&
&zzyaacr1(xyzzyaaah26+xyzzyaaan26))//'" and channel "'//trim(chname)//&
&'" are redundant, since they correspond to the same particle-nucleus &
&groups according to the provided rules.')
xyzzyaacr1(xyzzyaaah26+xyzzyaaan26)=chname
enddo
do xyzzyaaan26=1,xyzzyaabb26
if(trim(xyzzyaacr1(xyzzyaaah26+xyzzyaaan26))/='')cycle
xyzzyaaaq26=minloc(groups,groups==xyzzyaaan26)
xyzzyaacr1(xyzzyaaah26+xyzzyaaan26)=model_string(1,1,(/xyzzyaaaq26(2),&
&xyzzyaaaq26(1)/))
enddo
xyzzyaaab26=sum(xyzzyaacf1(1:xyzzyaabp26-1))+1
xyzzyaaac26=xyzzyaaab26+xyzzyaaav26-1
do xyzzyaaan26=1,xyzzyaabb26
xyzzyaaad26=xyzzyaabb1(xyzzyaaah26+xyzzyaaan26)+1
xyzzyaaae26=xyzzyaaad26+xyzzyaaav26-1
call read_gparam_channel(trim(label)//':Parameters:Channel '//trim(xyz&
&zyaacr1(xyzzyaaah26+xyzzyaaan26)),xyzzyaacf1(xyzzyaabp26),xyzzyaact1(&
&xyzzyaaab26:xyzzyaaac26),xyzzyaacx1(xyzzyaaad26:xyzzyaaae26),xyzzyaac&
&v1(xyzzyaaad26:xyzzyaaae26),xyzzyaadd1(xyzzyaaad26:xyzzyaaae26),xyzzy&
&aacz1(xyzzyaaad26:xyzzyaaae26),xyzzyaadf1(xyzzyaaad26:xyzzyaaae26),xy&
&zzyaadb1(xyzzyaaad26:xyzzyaaae26),xyzzyaadh1(xyzzyaaad26:xyzzyaaae26)&
&)
xyzzyaaaf26=xyzzyaabh1(xyzzyaaah26+xyzzyaaan26)+1
xyzzyaaai26=xyzzyaabj1(xyzzyaabp26)+1
xyzzyaaaj26=xyzzyaabl1(xyzzyaabp26)+1
call xyzzyaafc1(xyzzyaabr1(xyzzyaabp26),xyzzyaabt1(xyzzyaabp26),xyzzya&
&acf1(xyzzyaabp26),xyzzyaacl1(xyzzyaabp26),xyzzyaacx1(xyzzyaaad26),xyz&
&zyaacp1(xyzzyaabp26),xyzzyaadt1(xyzzyaaai26),xyzzyaacn1(xyzzyaabp26),&
&xyzzyaadv1(xyzzyaaaj26),xyzzyaadj1(xyzzyaaaf26))
enddo
endif
end subroutine read_enbasis
subroutine rewrite_eebasis(label,iset,iorder,print_determined)
implicit none
integer,intent(in) :: iset,iorder
logical,intent(in) :: print_determined
character(*),intent(in) :: label
integer xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xyzzyaaad27,xyzzyaaae27,xy&
&zzyaaaf27,xyzzyaaag27,xyzzyaaah27,xyzzyaaai27,xyzzyaaaj27,xyzzyaaak27&
&,xyzzyaaal27
character(casl_fullkeysize) namespace,btype
character(512) errmsg
if(iset<1)return
xyzzyaaaa27=xyzzyaabq1(iset)
if(xyzzyaaaa27/=0)then
btype=trim(xyzzyaaav1(xyzzyaaaa27))
xyzzyaaab27=index(btype,':')
namespace=''
if(xyzzyaaab27>0)then
namespace=btype(1:xyzzyaaab27-1)
btype=btype(xyzzyaaab27+1:)
endif
call set_casl_item(trim(label)//':Type',trim(btype),errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_EEBASIS',trim(errms&
&g))
else
namespace=''
btype='none'
endif
if(am_master)then
call wout(' Functional basis   : '//trim(btype))
if(trim(namespace)=='basis')then
if(iorder/=xyzzyaabs1(iset))then
call wout(' Expansion order    : '//trim(i2s(iorder))//' (of '//trim(i&
&2s(xyzzyaabs1(iset)))//' in shared basis)')
else
call wout(' Expansion order    : '//trim(i2s(iorder)))
endif
endif
endif
if(xyzzyaaci1(iset)>0)then
if(am_master)then
if(iorder/=xyzzyaabs1(iset))then
call wout(' Integer constants (of full shared basis):')
else
call wout(' Integer constants:')
endif
endif
xyzzyaaac27=xyzzyaabc1(iset)+1
xyzzyaaad27=xyzzyaaac27+xyzzyaaci1(iset)-1
call xyzzyaaeo1(trim(label)//':Constants',xyzzyaaci1(iset),xyzzyaado1(&
&xyzzyaaac27:xyzzyaaad27),xyzzyaadk1(xyzzyaaac27:xyzzyaaad27))
endif
if(xyzzyaacg1(iset)>0)then
if(am_master)then
if(iorder/=xyzzyaabs1(iset))then
call wout(' Real constants (of full shared basis):')
else
call wout(' Real constants:')
endif
endif
xyzzyaaae27=xyzzyaabe1(iset)+1
xyzzyaaaf27=xyzzyaaae27+xyzzyaacg1(iset)-1
call xyzzyaaep1(trim(label)//':Constants',xyzzyaacg1(iset),xyzzyaadq1(&
&xyzzyaaae27:xyzzyaaaf27),xyzzyaadm1(xyzzyaaae27:xyzzyaaaf27))
endif
if(xyzzyaace1(iset)>0)then
xyzzyaaag27=xyzzyaaay1(iset)
xyzzyaaah27=sum(xyzzyaace1(1:iset-1))+1
xyzzyaaai27=sum(xyzzyaace1(1:iset))
do xyzzyaaal27=1,xyzzyaabw1(iset)
xyzzyaaaj27=xyzzyaaba1(xyzzyaaag27+xyzzyaaal27)+1
xyzzyaaak27=xyzzyaaaj27+xyzzyaace1(iset)-1
if(am_master)call wout(' Channel '//trim(xyzzyaacq1(xyzzyaaag27+xyzzya&
&aal27))//':')
where(xyzzyaacu1(xyzzyaaaj27:xyzzyaaak27)==pflag_unset)xyzzyaacu1(xyzz&
&yaaaj27:xyzzyaaak27)=pflag_opt
call rewrite_gparam_channel(trim(label)//':Parameters:Channel '//trim(&
&xyzzyaacq1(xyzzyaaag27+xyzzyaaal27)),xyzzyaace1(iset),xyzzyaacs1(xyzz&
&yaaah27:xyzzyaaai27),xyzzyaacw1(xyzzyaaaj27:xyzzyaaak27),xyzzyaacu1(x&
&yzzyaaaj27:xyzzyaaak27),xyzzyaadc1(xyzzyaaaj27:xyzzyaaak27),xyzzyaacy&
&1(xyzzyaaaj27:xyzzyaaak27),xyzzyaade1(xyzzyaaaj27:xyzzyaaak27),xyzzya&
&ada1(xyzzyaaaj27:xyzzyaaak27),xyzzyaadg1(xyzzyaaaj27:xyzzyaaak27),pri&
&nt_determined)
enddo
endif
end subroutine rewrite_eebasis
subroutine rewrite_enbasis(label,iset,iorder,print_determined)
implicit none
integer,intent(in) :: iset,iorder
logical,intent(in) :: print_determined
character(*),intent(in) :: label
integer xyzzyaaaa28,xyzzyaaab28,xyzzyaaac28,xyzzyaaad28,xyzzyaaae28,xy&
&zzyaaaf28,xyzzyaaag28,xyzzyaaah28,xyzzyaaai28,xyzzyaaaj28,xyzzyaaak28&
&,xyzzyaaal28
character(casl_fullkeysize) namespace,btype
character(512) errmsg
if(iset<1)return
xyzzyaaaa28=xyzzyaabr1(iset)
if(xyzzyaaaa28/=0)then
btype=trim(xyzzyaaav1(xyzzyaaaa28))
xyzzyaaab28=index(btype,':')
namespace=''
if(xyzzyaaab28>0)then
namespace=btype(1:xyzzyaaab28-1)
btype=btype(xyzzyaaab28+1:)
endif
call set_casl_item(trim(label)//':Type',trim(btype),errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_ENBASIS',trim(errms&
&g))
else
namespace=''
btype='none'
endif
if(am_master)then
call wout(' Functional basis   : '//trim(btype))
if(trim(namespace)=='basis')then
if(iorder/=xyzzyaabt1(iset))then
call wout(' Expansion order    : '//trim(i2s(iorder))//' (of '//trim(i&
&2s(xyzzyaabt1(iset)))//' in shared basis)')
else
call wout(' Expansion order    : '//trim(i2s(iorder)))
endif
endif
endif
if(xyzzyaacj1(iset)>0)then
if(am_master)then
if(iorder/=xyzzyaabt1(iset))then
call wout(' Integer constants (of full shared basis):')
else
call wout(' Integer constants:')
endif
endif
xyzzyaaac28=xyzzyaabd1(iset)+1
xyzzyaaad28=xyzzyaaac28+xyzzyaacj1(iset)-1
call xyzzyaaeo1(trim(label)//':Constants',xyzzyaacj1(iset),xyzzyaadp1(&
&xyzzyaaac28:xyzzyaaad28),xyzzyaadl1(xyzzyaaac28:xyzzyaaad28))
endif
if(xyzzyaach1(iset)>0)then
if(am_master)then
if(iorder/=xyzzyaabt1(iset))then
call wout(' Real constants (of full shared basis):')
else
call wout(' Real constants:')
endif
endif
xyzzyaaae28=xyzzyaabf1(iset)+1
xyzzyaaaf28=xyzzyaaae28+xyzzyaach1(iset)-1
call xyzzyaaep1(trim(label)//':Constants',xyzzyaach1(iset),xyzzyaadr1(&
&xyzzyaaae28:xyzzyaaaf28),xyzzyaadn1(xyzzyaaae28:xyzzyaaaf28))
endif
if(xyzzyaacf1(iset)>0)then
xyzzyaaag28=xyzzyaaaz1(iset)
xyzzyaaah28=sum(xyzzyaacf1(1:iset-1))+1
xyzzyaaai28=sum(xyzzyaacf1(1:iset))
do xyzzyaaal28=1,xyzzyaabx1(iset)
xyzzyaaaj28=xyzzyaabb1(xyzzyaaag28+xyzzyaaal28)+1
xyzzyaaak28=xyzzyaaaj28+xyzzyaacf1(iset)-1
if(am_master)call wout(' Channel '//trim(xyzzyaacr1(xyzzyaaag28+xyzzya&
&aal28))//':')
where(xyzzyaacv1(xyzzyaaaj28:xyzzyaaak28)==pflag_unset)xyzzyaacv1(xyzz&
&yaaaj28:xyzzyaaak28)=pflag_opt
call rewrite_gparam_channel(trim(label)//':Parameters:Channel '//trim(&
&xyzzyaacr1(xyzzyaaag28+xyzzyaaal28)),xyzzyaacf1(iset),xyzzyaact1(xyzz&
&yaaah28:xyzzyaaai28),xyzzyaacx1(xyzzyaaaj28:xyzzyaaak28),xyzzyaacv1(x&
&yzzyaaaj28:xyzzyaaak28),xyzzyaadd1(xyzzyaaaj28:xyzzyaaak28),xyzzyaacz&
&1(xyzzyaaaj28:xyzzyaaak28),xyzzyaadf1(xyzzyaaaj28:xyzzyaaak28),xyzzya&
&adb1(xyzzyaaaj28:xyzzyaaak28),xyzzyaadh1(xyzzyaaaj28:xyzzyaaak28),pri&
&nt_determined)
enddo
endif
end subroutine rewrite_enbasis
subroutine xyzzyaaeo1(label,nconst_int,const_int_name,const_int)
implicit none
integer,intent(in) :: nconst_int
integer,intent(in) :: const_int(nconst_int)
character(*),intent(in) :: label,const_int_name(nconst_int)
integer xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaae29
logical xyzzyaaaf29,xyzzyaaag29
character(80) tmpr
character(casl_fullkeysize) cname,pname,prev_pname,tstr,tstr2,tprint,t&
&print2
character(512) errmsg
do xyzzyaaaa29=1,nconst_int
cname=const_int_name(xyzzyaaaa29)
xyzzyaaac29=scan(cname,':',back=.true.)
if(xyzzyaaac29<1)xyzzyaaac29=0
if(cname(xyzzyaaac29+1:xyzzyaaac29+2)=='%u')then
call set_casl_block(trim(label)//':'//cname(1:xyzzyaaac29-1),errmsg,pr&
&efer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GCONST_INT',trim(er&
&rmsg))
endif
call set_casl_item(trim(label)//':'//trim(cname),const_int(xyzzyaaaa29&
&),errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GCONST_INT',trim(er&
&rmsg))
enddo
if(am_master)then
xyzzyaaab29=0
do xyzzyaaaa29=1,nconst_int
cname=const_int_name(xyzzyaaaa29)
xyzzyaaac29=scan(cname,':',back=.true.)
if(xyzzyaaac29<1)xyzzyaaac29=0
if(cname(xyzzyaaac29+1:xyzzyaaac29+2)=='%u')then
cname=cname(1:xyzzyaaac29-1)
xyzzyaaac29=scan(cname,':',back=.true.)
if(xyzzyaaac29<1)xyzzyaaac29=0
endif
tstr=cname(1:xyzzyaaac29-1)
cname=cname(xyzzyaaac29+1:)
xyzzyaaae29=0
if(xyzzyaaac29>0)then
xyzzyaaae29=1
do
xyzzyaaac29=scan(tstr,':')
if(xyzzyaaac29<1)exit
xyzzyaaae29=xyzzyaaae29+1
tstr=tstr(xyzzyaaac29+1:)
enddo
endif
xyzzyaaab29=max(xyzzyaaab29,len_trim(cname)+xyzzyaaae29)
enddo
xyzzyaaaa29=0
prev_pname=''
do while(xyzzyaaaa29<nconst_int)
xyzzyaaaa29=xyzzyaaaa29+1
cname=const_int_name(xyzzyaaaa29)
xyzzyaaac29=scan(cname,':',back=.true.)
if(xyzzyaaac29<1)then
xyzzyaaaf29=.false.
pname=''
prev_pname=''
else
xyzzyaaaf29=cname(xyzzyaaac29+1:xyzzyaaac29+2)=='%u'
if(xyzzyaaaf29)then
cname=cname(1:xyzzyaaac29-1)
xyzzyaaac29=scan(cname,':',back=.true.)
if(xyzzyaaac29<1)then
pname=''
else
pname=cname(1:xyzzyaaac29-1)
cname=cname(xyzzyaaac29+1:)
endif
else
pname=cname(1:xyzzyaaac29-1)
cname=cname(xyzzyaaac29+1:)
endif
endif
xyzzyaaae29=0
if(len_trim(pname)>0)then
tstr=pname
xyzzyaaag29=.false.
if(len_trim(prev_pname)>0)then
tstr2=prev_pname
xyzzyaaag29=.true.
endif
do
xyzzyaaac29=scan(tstr,':')
if(xyzzyaaac29<1)xyzzyaaac29=len_trim(tstr)+1
tprint=tstr(1:xyzzyaaac29-1)
if(xyzzyaaag29)then
xyzzyaaad29=scan(tstr2,':')
if(xyzzyaaad29<1)xyzzyaaad29=len_trim(tstr2)+1
tprint2=tstr2(1:xyzzyaaad29-1)
if(trim(tprint2)/=trim(tprint))xyzzyaaag29=.false.
endif
if(.not.xyzzyaaag29)call wout(repeat(' ',2+xyzzyaaae29)//trim(tprint)/&
&/':')
xyzzyaaae29=xyzzyaaae29+1
if(xyzzyaaac29>=len_trim(tstr))exit
tstr=tstr(xyzzyaaac29+1:)
enddo
prev_pname=pname
endif
if(xyzzyaaaf29)then
tprint='[ '//trim(i2s(const_int(xyzzyaaaa29)))
do
if(xyzzyaaaa29==nconst_int)exit
tstr=const_int_name(xyzzyaaaa29+1)
tstr2=''
if(len_trim(prev_pname)>0)tstr2=trim(prev_pname)//':'
tstr2=trim(tstr2)//trim(cname)//':%u'
tstr=tstr(1:len_trim(tstr2))
if(trim(tstr)/=trim(tstr2))exit
xyzzyaaaa29=xyzzyaaaa29+1
tprint=trim(tprint)//', '//trim(i2s(const_int(xyzzyaaaa29)))
enddo
tprint=trim(tprint)//' ]'
write(tmpr,'('//trim(i2s(2+xyzzyaaae29))//'x,a,t'//trim(i2s(3+xyzzyaaa&
&b29))//'," = ",a)')trim(cname),trim(tprint)
call wout(trim(tmpr))
else
write(tmpr,'('//trim(i2s(2+xyzzyaaae29))//'x,a,t'//trim(i2s(3+xyzzyaaa&
&b29))//'," = ",a)')trim(cname),trim(i2s(const_int(xyzzyaaaa29)))
call wout(trim(tmpr))
endif
enddo
endif
end subroutine xyzzyaaeo1
subroutine xyzzyaaep1(label,nconst_dble,const_dble_name,const_dble)
implicit none
integer,intent(in) :: nconst_dble
real(dp),intent(in) :: const_dble(nconst_dble)
character(*),intent(in) :: label,const_dble_name(nconst_dble)
integer xyzzyaaaa30,xyzzyaaab30,xyzzyaaac30,xyzzyaaad30,xyzzyaaae30
logical xyzzyaaaf30,xyzzyaaag30
character(80) tmpr,tmpr2
character(casl_fullkeysize) cname,pname,prev_pname,tstr,tstr2,tprint,t&
&print2
character(512) errmsg
do xyzzyaaaa30=1,nconst_dble
cname=const_dble_name(xyzzyaaaa30)
xyzzyaaac30=scan(cname,':',back=.true.)
if(xyzzyaaac30<1)xyzzyaaac30=0
if(cname(xyzzyaaac30+1:xyzzyaaac30+2)=='%u')then
call set_casl_block(trim(label)//':'//cname(1:xyzzyaaac30-1),errmsg,pr&
&efer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GCONST_DBLE',trim(e&
&rrmsg))
endif
call set_casl_item(trim(label)//':'//trim(cname),const_dble(xyzzyaaaa3&
&0),errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GCONST_DBLE',trim(e&
&rrmsg))
enddo
if(am_master)then
xyzzyaaab30=0
do xyzzyaaaa30=1,nconst_dble
cname=const_dble_name(xyzzyaaaa30)
xyzzyaaac30=scan(cname,':',back=.true.)
if(xyzzyaaac30<1)xyzzyaaac30=0
if(cname(xyzzyaaac30+1:xyzzyaaac30+2)=='%u')then
cname=cname(1:xyzzyaaac30-1)
xyzzyaaac30=scan(cname,':',back=.true.)
if(xyzzyaaac30<1)xyzzyaaac30=0
endif
tstr=cname(1:xyzzyaaac30-1)
cname=cname(xyzzyaaac30+1:)
xyzzyaaae30=0
if(xyzzyaaac30>0)then
xyzzyaaae30=1
do
xyzzyaaac30=scan(tstr,':')
if(xyzzyaaac30<1)exit
xyzzyaaae30=xyzzyaaae30+1
tstr=tstr(xyzzyaaac30+1:)
enddo
endif
xyzzyaaab30=max(xyzzyaaab30,len_trim(cname)+xyzzyaaae30)
enddo
xyzzyaaaa30=0
prev_pname=''
do while(xyzzyaaaa30<nconst_dble)
xyzzyaaaa30=xyzzyaaaa30+1
cname=const_dble_name(xyzzyaaaa30)
xyzzyaaac30=scan(cname,':',back=.true.)
if(xyzzyaaac30<1)then
xyzzyaaaf30=.false.
pname=''
prev_pname=''
else
xyzzyaaaf30=cname(xyzzyaaac30+1:xyzzyaaac30+2)=='%u'
if(xyzzyaaaf30)then
cname=cname(1:xyzzyaaac30-1)
xyzzyaaac30=scan(cname,':',back=.true.)
if(xyzzyaaac30<1)then
pname=''
else
pname=cname(1:xyzzyaaac30-1)
cname=cname(xyzzyaaac30+1:)
endif
else
pname=cname(1:xyzzyaaac30-1)
cname=cname(xyzzyaaac30+1:)
endif
endif
xyzzyaaae30=0
if(len_trim(pname)>0)then
tstr=pname
xyzzyaaag30=.false.
if(len_trim(prev_pname)>0)then
tstr2=prev_pname
xyzzyaaag30=.true.
endif
do
xyzzyaaac30=scan(tstr,':')
if(xyzzyaaac30<1)xyzzyaaac30=len_trim(tstr)+1
tprint=tstr(1:xyzzyaaac30-1)
if(xyzzyaaag30)then
xyzzyaaad30=scan(tstr2,':')
if(xyzzyaaad30<1)xyzzyaaad30=len_trim(tstr2)+1
tprint2=tstr2(1:xyzzyaaad30-1)
if(trim(tprint2)/=trim(tprint))xyzzyaaag30=.false.
endif
if(.not.xyzzyaaag30)call wout(repeat(' ',2+xyzzyaaae30)//trim(tprint)/&
&/':')
xyzzyaaae30=xyzzyaaae30+1
if(xyzzyaaac30>=len_trim(tstr))exit
tstr=tstr(xyzzyaaac30+1:)
enddo
prev_pname=pname
endif
if(xyzzyaaaf30)then
tmpr=r2s2(const_dble(xyzzyaaaa30),'(es16.8)')
tprint='[ '//trim(tmpr)
do
if(xyzzyaaaa30==nconst_dble)exit
tstr=const_dble_name(xyzzyaaaa30+1)
tstr2=''
if(len_trim(prev_pname)>0)tstr2=trim(prev_pname)//':'
tstr2=trim(tstr2)//trim(cname)//':%u'
tstr=tstr(1:len_trim(tstr2))
if(trim(tstr)/=trim(tstr2))exit
xyzzyaaaa30=xyzzyaaaa30+1
tmpr=r2s2(const_dble(xyzzyaaaa30),'(es16.8)')
tprint=trim(tprint)//', '//trim(tmpr)
enddo
tprint=trim(tprint)//' ]'
write(tmpr2,'('//trim(i2s(2+xyzzyaaae30))//'x,a,t'//trim(i2s(3+xyzzyaa&
&ab30))//'," = ",a)')trim(cname),trim(tprint)
call wout(trim(tmpr2))
else
tmpr=r2s2(const_dble(xyzzyaaaa30),'(es16.8)')
write(tmpr2,'('//trim(i2s(2+xyzzyaaae30))//'x,a,t'//trim(i2s(3+xyzzyaa&
&ab30))//'," = ",a)')trim(cname),trim(tmpr)
call wout(trim(tmpr2))
endif
enddo
endif
end subroutine xyzzyaaep1
subroutine read_gparam_channel(label,nparam,param_name,param,param_opt&
&able,has_lolim,lolim,has_hilim,hilim,is_shallow)
implicit none
integer,intent(in) :: nparam
integer,intent(inout) :: param_optable(nparam)
real(dp),intent(inout) :: param(nparam),lolim(nparam),hilim(nparam)
logical,intent(inout) :: has_lolim(nparam),has_hilim(nparam),is_shallo&
&w(nparam)
character(*),intent(in) :: label,param_name(nparam)
integer xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31
real(dp) xyzzyaaad31,xyzzyaaae31,xyzzyaaaf31
logical exists,is_block,xyzzyaaag31,xyzzyaaah31,xyzzyaaai31,xyzzyaaaj3&
&1
character(casl_keysize) param_test
character(casl_fullkeysize) optflag,limstring
do xyzzyaaaa31=1,nparam
call get_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa31))//':&
&%u1',param_test,xyzzyaaab31)
xyzzyaaah31=.false.
xyzzyaaag31=.false.
if(xyzzyaaab31==0)then
xyzzyaaag31=.true.
if(trim(param_test)=='default')then
continue
else
xyzzyaaac31=scan(param_test,'%')
if(xyzzyaaac31>0)then
if(len_trim(param_test(xyzzyaaac31+1:))/=0)call errstop_master('READ_G&
&PARAM_CHANNEL','Trailing characters after "%" parsing value of parame&
&ter "'//trim(param_name(xyzzyaaaa31))//'".')
param_test=param_test(:xyzzyaaac31-1)
read(param_test,*,iostat=xyzzyaaab31)xyzzyaaad31
if(xyzzyaaab31/=0)call errstop_master('READ_GPARAM_CHANNEL','Could not&
& parse percent initialization string for parameter "'//trim(param_nam&
&e(xyzzyaaaa31))//'".')
if(xyzzyaaad31<0.d0.or.xyzzyaaad31>100.d0)call errstop_master('READ_GP&
&ARAM_CHANNEL','Percent initialization for parameter "'//trim(param_na&
&me(xyzzyaaaa31))//'" out of range.')
xyzzyaaah31=.true.
else
call get_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa31))//':&
&%u1',xyzzyaaad31,xyzzyaaab31)
if(xyzzyaaab31/=0)call errstop_master('READ_GPARAM_CHANNEL','Could not&
& parse value of parameter "'//trim(param_name(xyzzyaaaa31))//'".')
param(xyzzyaaaa31)=xyzzyaaad31
endif
endif
endif
call get_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa31))//':&
&%u2',optflag,xyzzyaaab31)
if(xyzzyaaab31==0)then
select case(trim(optflag))
case('determined','-1')
param_optable(xyzzyaaaa31)=pflag_det
case('fixed','0')
param_optable(xyzzyaaaa31)=pflag_fix
case('opt','optimizable','1')
param_optable(xyzzyaaaa31)=pflag_opt
case default
call errstop_master('READ_GPARAM_CHANNEL','Bad optflag "'//trim(optfla&
&g)//'" for parameter "'//trim(param_name(xyzzyaaaa31))//'".')
end select
elseif(xyzzyaaag31)then
param_optable(xyzzyaaaa31)=pflag_opt
endif
call query_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa31))//&
&':limits',exists=exists,is_block=is_block)
if(exists)then
if(.not.is_block)call errstop_master('READ_GPARAM_CHANNEL','Non-block &
&"limits" not implemented.')
xyzzyaaae31=hilim(xyzzyaaaa31)
xyzzyaaaf31=lolim(xyzzyaaaa31)
xyzzyaaai31=has_hilim(xyzzyaaaa31)
xyzzyaaaj31=has_lolim(xyzzyaaaa31)
call get_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa31))//':&
&limits:%u1',limstring,xyzzyaaab31)
if(xyzzyaaab31/=0)call errstop_master('READ_GPARAM_CHANNEL','Could not&
& find lower limit for parameter "'//trim(param_name(xyzzyaaaa31))//'"&
&.')
select case(trim(limstring))
case('default')
continue
case('none','-Inf','Inf','-inf','inf')
if(has_lolim(xyzzyaaaa31))call errstop_master('READ_GPARAM_CHANNEL','A&
&ttempt to override requirement of lower limit for parameter "'//trim(&
&param_name(xyzzyaaaa31))//'"')
case default
xyzzyaaac31=scan(limstring,'%')
if(xyzzyaaac31>0)then
if(.not.xyzzyaaaj31.or..not.xyzzyaaai31)call errstop_master('READ_GPAR&
&AM_CHANNEL','Percent initialization used for lower limit of parameter&
& '//trim(param_name(xyzzyaaaa31))//', but this requires parameter to &
&have default upper and lower limits.')
if(len_trim(limstring(xyzzyaaac31+1:))/=0)call errstop_master('READ_GP&
&ARAM_CHANNEL','Trailing characters after "%" parsing lower limit of p&
&arameter "'//trim(param_name(xyzzyaaaa31))//'".')
limstring=limstring(:xyzzyaaac31-1)
read(limstring,*,iostat=xyzzyaaab31)lolim(xyzzyaaaa31)
if(xyzzyaaab31/=0)call errstop_master('READ_GPARAM_CHANNEL','Could not&
& parse percent initialization string for lower limit of parameter "'/&
&/trim(param_name(xyzzyaaaa31))//'".')
if(lolim(xyzzyaaaa31)<0.d0.or.lolim(xyzzyaaaa31)>100.d0)call errstop_m&
&aster('READ_GPARAM_CHANNEL','Percent initialization for lower limit o&
&f parameter "'//trim(param_name(xyzzyaaaa31))//'" out of range.')
if(lolim(xyzzyaaaa31)<=0.d0)then
lolim(xyzzyaaaa31)=xyzzyaaaf31
elseif(lolim(xyzzyaaaa31)>=100.d0)then
lolim(xyzzyaaaa31)=xyzzyaaae31
else
lolim(xyzzyaaaa31)=xyzzyaaaf31+(lolim(xyzzyaaaa31)/100.d0)*(xyzzyaaae3&
&1-xyzzyaaaf31)
endif
else
call get_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa31))//':&
&limits:%u1',lolim(xyzzyaaaa31),xyzzyaaab31)
if(xyzzyaaab31/=0)call errstop_master('READ_GPARAM_CHANNEL','Specified&
& lower limit for parameter "'//trim(param_name(xyzzyaaaa31))//'" is n&
&ot a real number or identifiable keyword.')
if(xyzzyaaaj31)lolim(xyzzyaaaa31)=max(xyzzyaaaf31,lolim(xyzzyaaaa31))
has_lolim(xyzzyaaaa31)=.true.
endif
end select
call get_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa31))//':&
&limits:%u2',limstring,xyzzyaaab31)
if(xyzzyaaab31/=0)call errstop_master('READ_GPARAM_CHANNEL','Could not&
& find upper limit for parameter "'//trim(param_name(xyzzyaaaa31))//'"&
&.')
select case(trim(limstring))
case('default')
continue
case('none','+Inf','Inf','+inf','inf')
if(has_hilim(xyzzyaaaa31))call errstop_master('READ_GPARAM_CHANNEL','A&
&ttempt to override requirement of upper limit for parameter "'//trim(&
&param_name(xyzzyaaaa31))//'"')
case default
xyzzyaaac31=scan(limstring,'%')
if(xyzzyaaac31>0)then
if(.not.xyzzyaaaj31.or..not.xyzzyaaai31)call errstop_master('READ_GPAR&
&AM_CHANNEL','Percent initialization used for upper limit of parameter&
& '//trim(param_name(xyzzyaaaa31))//', but this requires the parameter&
& to have default upper and lower limits.')
if(len_trim(limstring(xyzzyaaac31+1:))/=0)call errstop_master('READ_GP&
&ARAM_CHANNEL','Trailing characters after "%" parsing upper limit of p&
&arameter "'//trim(param_name(xyzzyaaaa31))//'".')
limstring=limstring(:xyzzyaaac31-1)
read(limstring,*,iostat=xyzzyaaab31)hilim(xyzzyaaaa31)
if(xyzzyaaab31/=0)call errstop_master('READ_GPARAM_CHANNEL','Could not&
& parse percent initialization string for upper limit of parameter "'/&
&/trim(param_name(xyzzyaaaa31))//'".')
if(hilim(xyzzyaaaa31)<0.d0.or.hilim(xyzzyaaaa31)>100.d0)call errstop_m&
&aster('READ_GPARAM_CHANNEL','Percent initialization for upper limit o&
&f parameter "'//trim(param_name(xyzzyaaaa31))//'" out of range.')
if(hilim(xyzzyaaaa31)<=0.d0)then
hilim(xyzzyaaaa31)=xyzzyaaaf31
elseif(hilim(xyzzyaaaa31)>=100.d0)then
hilim(xyzzyaaaa31)=xyzzyaaae31
else
hilim(xyzzyaaaa31)=xyzzyaaaf31+(hilim(xyzzyaaaa31)/100.d0)*(xyzzyaaae3&
&1-xyzzyaaaf31)
endif
else
call get_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa31))//':&
&limits:%u2',hilim(xyzzyaaaa31),xyzzyaaab31)
if(xyzzyaaab31/=0)call errstop_master('READ_GPARAM_CHANNEL','Specified&
& upper limit for parameter "'//trim(param_name(xyzzyaaaa31))//'" is n&
&ot a real number or identifiable keyword.')
if(xyzzyaaai31)hilim(xyzzyaaaa31)=min(xyzzyaaae31,hilim(xyzzyaaaa31))
has_hilim(xyzzyaaaa31)=.true.
endif
end select
endif
if(has_lolim(xyzzyaaaa31).and.has_hilim(xyzzyaaaa31).and.lolim(xyzzyaa&
&aa31)>=hilim(xyzzyaaaa31))call errstop_master('READ_GPARAM_CHANNEL','&
&Lower limit for parameter "'//trim(param_name(xyzzyaaaa31))//'" >= up&
&per limit.')
if(xyzzyaaah31)then
if(.not.has_lolim(xyzzyaaaa31).or..not.has_hilim(xyzzyaaaa31))call err&
&stop_master('READ_GPARAM_CHANNEL','Percent initialization used for pa&
&rameter '//trim(param_name(xyzzyaaaa31))//', but this requires the pa&
&rameter to have upper and lower limits defined.')
if(xyzzyaaad31<=0.d0)then
param(xyzzyaaaa31)=lolim(xyzzyaaaa31)
elseif(xyzzyaaad31>=100.d0)then
param(xyzzyaaaa31)=hilim(xyzzyaaaa31)
else
param(xyzzyaaaa31)=lolim(xyzzyaaaa31)+(xyzzyaaad31/100.d0)*(hilim(xyzz&
&yaaaa31)-lolim(xyzzyaaaa31))
endif
elseif(has_lolim(xyzzyaaaa31).and.param(xyzzyaaaa31)<lolim(xyzzyaaaa31&
&))then
param(xyzzyaaaa31)=lolim(xyzzyaaaa31)
elseif(has_hilim(xyzzyaaaa31).and.param(xyzzyaaaa31)>hilim(xyzzyaaaa31&
&))then
param(xyzzyaaaa31)=hilim(xyzzyaaaa31)
endif
enddo
end subroutine read_gparam_channel
subroutine rewrite_gparam_channel(label,nparam,param_name,param,param_&
&optable,has_lolim,lolim,has_hilim,hilim,is_shallow,print_determined)
implicit none
integer,intent(in) :: nparam
integer,intent(inout) :: param_optable(nparam)
real(dp),intent(inout) :: param(nparam),lolim(nparam),hilim(nparam)
logical,intent(in) :: print_determined
logical,intent(inout) :: has_lolim(nparam),has_hilim(nparam),is_shallo&
&w(nparam)
character(*),intent(in) :: label,param_name(nparam)
integer xyzzyaaaa32,xyzzyaaab32
character(80) tmpr,tmpr2
character(casl_fullkeysize) optflag,limstring
character(512) errmsg
call set_casl_block(trim(label),errmsg,force_noninline=.true.)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GPARAM_CHANNEL',tri&
&m(errmsg))
if(am_master)then
xyzzyaaab32=0
do xyzzyaaaa32=1,nparam
xyzzyaaab32=max(xyzzyaaab32,len_trim(param_name(xyzzyaaaa32)))
enddo
endif
do xyzzyaaaa32=1,nparam
select case(param_optable(xyzzyaaaa32))
case(pflag_det)
optflag='determined'
if(.not.print_determined)cycle
case(pflag_fix)
optflag='fixed'
case(pflag_opt)
optflag='optimizable'
case default
cycle
end select
call set_casl_block(trim(label)//':'//trim(param_name(xyzzyaaaa32)),er&
&rmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GPARAM_CHANNEL',tri&
&m(errmsg))
call set_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa32))//':&
&%u1',param(xyzzyaaaa32),errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GPARAM_CHANNEL',tri&
&m(errmsg))
call set_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa32))//':&
&%u2',trim(optflag),errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GPARAM_CHANNEL',tri&
&m(errmsg))
if(am_master)limstring=''
if(has_lolim(xyzzyaaaa32).or.has_hilim(xyzzyaaaa32))then
if(am_master)then
tmpr=r2s(lolim(xyzzyaaaa32),'(es12.4)')
tmpr2=r2s(hilim(xyzzyaaaa32),'(es12.4)')
if(has_lolim(xyzzyaaaa32).and.has_hilim(xyzzyaaaa32))then
limstring=' in ['//trim(tmpr)//', '//trim(tmpr2)//']'
elseif(has_lolim(xyzzyaaaa32))then
limstring=' in ['//trim(tmpr)//', +Inf)'
else
limstring=' in (-Inf, '//trim(tmpr2)//']'
endif
endif
if(has_lolim(xyzzyaaaa32))then
call set_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa32))//':&
&limits:%u1',lolim(xyzzyaaaa32),errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GPARAM_CHANNEL',tri&
&m(errmsg))
else
call set_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa32))//':&
&limits:%u1','-Inf',errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GPARAM_CHANNEL',tri&
&m(errmsg))
endif
if(has_hilim(xyzzyaaaa32))then
call set_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa32))//':&
&limits:%u2',hilim(xyzzyaaaa32),errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GPARAM_CHANNEL',tri&
&m(errmsg))
else
call set_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa32))//':&
&limits:%u2','+Inf',errmsg)
if(len_trim(errmsg)>0)call errstop_master('REWRITE_GPARAM_CHANNEL',tri&
&m(errmsg))
endif
endif
if(am_master)then
tmpr=r2s2(param(xyzzyaaaa32),'(es16.8)')
write(tmpr2,'(2x,a,t'//trim(i2s(4+xyzzyaaab32))//',a13," = ",a,a)')tri&
&m(param_name(xyzzyaaaa32)),'('//trim(optflag)//')',trim(tmpr),     tr&
&im(limstring)
call wout(tmpr2)
endif
enddo
end subroutine rewrite_gparam_channel
subroutine update_eebasis_casl(label,iset)
implicit none
integer,intent(in) :: iset
character(*),intent(in) :: label
integer xyzzyaaaa33,xyzzyaaab33,xyzzyaaac33,xyzzyaaad33,xyzzyaaae33,xy&
&zzyaaaf33,xyzzyaaag33,xyzzyaaah33
if(iset<1)return
xyzzyaaaa33=xyzzyaace1(iset)
if(xyzzyaaaa33==0)return
xyzzyaaab33=xyzzyaabw1(iset)
xyzzyaaad33=xyzzyaaay1(iset)
xyzzyaaae33=sum(xyzzyaace1(1:iset-1))+1
xyzzyaaaf33=xyzzyaaae33+xyzzyaaaa33-1
do xyzzyaaac33=1,xyzzyaaab33
xyzzyaaag33=xyzzyaaba1(xyzzyaaad33+xyzzyaaac33)+1
xyzzyaaah33=xyzzyaaag33+xyzzyaaaa33-1
call update_gparam_channel_casl(trim(label)//':Parameters:Channel '//t&
&rim(xyzzyaacq1(xyzzyaaad33+xyzzyaaac33)),xyzzyaaaa33,xyzzyaacs1(xyzzy&
&aaae33:xyzzyaaaf33),xyzzyaacw1(xyzzyaaag33:xyzzyaaah33),xyzzyaacu1(xy&
&zzyaaag33:xyzzyaaah33),.false.)
enddo
end subroutine update_eebasis_casl
subroutine update_enbasis_casl(label,iset)
implicit none
integer,intent(in) :: iset
character(*),intent(in) :: label
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34,xyzzyaaad34,xyzzyaaae34,xy&
&zzyaaaf34,xyzzyaaag34,xyzzyaaah34
if(iset<1)return
xyzzyaaaa34=xyzzyaacf1(iset)
if(xyzzyaaaa34==0)return
xyzzyaaab34=xyzzyaabx1(iset)
xyzzyaaad34=xyzzyaaaz1(iset)
xyzzyaaae34=sum(xyzzyaacf1(1:iset-1))+1
xyzzyaaaf34=xyzzyaaae34+xyzzyaaaa34-1
do xyzzyaaac34=1,xyzzyaaab34
xyzzyaaag34=xyzzyaabb1(xyzzyaaad34+xyzzyaaac34)+1
xyzzyaaah34=xyzzyaaag34+xyzzyaaaa34-1
call update_gparam_channel_casl(trim(label)//':Parameters:Channel '//t&
&rim(xyzzyaacr1(xyzzyaaad34+xyzzyaaac34)),xyzzyaaaa34,xyzzyaact1(xyzzy&
&aaae34:xyzzyaaaf34),xyzzyaacx1(xyzzyaaag34:xyzzyaaah34),xyzzyaacv1(xy&
&zzyaaag34:xyzzyaaah34),.false.)
enddo
end subroutine update_enbasis_casl
subroutine update_gparam_channel_casl(label,nparam,param_name,param,pa&
&ram_optable,print_determined)
implicit none
integer,intent(in) :: nparam,param_optable(nparam)
real(dp),intent(in) :: param(nparam)
logical,intent(in) :: print_determined
character(*),intent(in) :: label,param_name(nparam)
character(512) errmsg
integer xyzzyaaaa35
character(20),parameter :: xyzzyaaab35(-1:1)=(/'determined ','fixed   &
&   ','optimizable'/)
do xyzzyaaaa35=1,nparam
if(param_optable(xyzzyaaaa35)==pflag_det.and..not.print_determined)the&
&n
call delete_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa35)))
cycle
endif
call set_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa35))//':&
&%u1',param(xyzzyaaaa35),errmsg)
if(len_trim(errmsg)>0)call errstop_master('UPDATE_GPARAM_CHANNEL_CASL'&
&,trim(errmsg))
call set_casl_item(trim(label)//':'//trim(param_name(xyzzyaaaa35))//':&
&%u2',trim(xyzzyaaab35(param_optable(xyzzyaaaa35))),errmsg)
if(len_trim(errmsg)>0)call errstop_master('UPDATE_GPARAM_CHANNEL_CASL'&
&,trim(errmsg))
enddo
end subroutine update_gparam_channel_casl
subroutine init_groups_ee(groups_ee)
implicit none
integer,intent(out) :: groups_ee(nspin,nspin)
integer xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36
groups_ee=0
xyzzyaaac36=0
do xyzzyaaaa36=1,nspin
if(nele(xyzzyaaaa36)==0)cycle
if(nele(xyzzyaaaa36)>1)then
xyzzyaaac36=xyzzyaaac36+1
groups_ee(xyzzyaaaa36,xyzzyaaaa36)=xyzzyaaac36
endif
do xyzzyaaab36=xyzzyaaaa36+1,nspin
if(nele(xyzzyaaab36)==0)cycle
xyzzyaaac36=xyzzyaaac36+1
groups_ee(xyzzyaaaa36,xyzzyaaab36)=xyzzyaaac36
groups_ee(xyzzyaaab36,xyzzyaaaa36)=xyzzyaaac36
enddo
enddo
end subroutine init_groups_ee
subroutine init_groups_en(groups_en)
implicit none
integer,intent(out) :: groups_en(nitot,nspin)
integer xyzzyaaaa37,xyzzyaaab37,xyzzyaaac37
groups_en=0
xyzzyaaac37=0
do xyzzyaaaa37=1,nspin
if(nele(xyzzyaaaa37)==0)cycle
do xyzzyaaab37=1,nitot
xyzzyaaac37=xyzzyaaac37+1
groups_en(xyzzyaaab37,xyzzyaaaa37)=xyzzyaaac37
enddo
enddo
end subroutine init_groups_en
subroutine digest_rule(rule,has_ee,has_en,groups_ee,groups_en)
implicit none
integer,intent(inout) :: groups_ee(nspin,nspin),groups_en(nitot,nspin)
logical,intent(in) :: has_ee,has_en
character(*),intent(in) :: rule
character(len(rule)) trule
character(1) ion1_match,ion2_match
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38,xyzzyaaad38,xyzzyaaae38,xy&
&zzyaaaf38,xyzzyaaag38,xyzzyaaah38,xyzzyaaai38,xyzzyaaaj38,xyzzyaaak38&
&,xyzzyaaal38,xyzzyaaam38(nitot),xyzzyaaan38
logical xyzzyaaao38,xyzzyaaap38,xyzzyaaaq38,xyzzyaaar38,xyzzyaaas38
trule=rule
xyzzyaaao38=.false.
xyzzyaaap38=.false.
xyzzyaaaq38=.false.
xyzzyaaar38=.false.
xyzzyaaas38=.false.
if(trim(adjustl(trule))=='N')then
xyzzyaaao38=.true.
elseif(trim(adjustl(trule))=='Z')then
xyzzyaaap38=.true.
else
xyzzyaaaa38=index(trule,'=')
if(xyzzyaaaa38>0)then
xyzzyaaar38=.true.
else
xyzzyaaaa38=index(trule,'!')
if(xyzzyaaaa38>0)then
if(len_trim(trule(1:xyzzyaaaa38-1))>0)call errstop_master('DIGEST_RULE&
&S','Error parsing removal rule: leading characters before "!".')
xyzzyaaaq38=.true.
else
xyzzyaaas38=.true.
endif
endif
endif
if(xyzzyaaao38)then
if(.not.isperiodic)return
if(.not.has_en)return
do xyzzyaaad38=1,nbasis
do xyzzyaaan38=2,npcells
xyzzyaaak38=xyzzyaaad38+(xyzzyaaan38-1)*nbasis
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
enddo
enddo
elseif(xyzzyaaap38)then
if(.not.has_en)return
xyzzyaaal38=0
xyzzyaaam38=0
do xyzzyaaad38=1,nitot
if(any(xyzzyaaam38(1:xyzzyaaal38)==iontype(xyzzyaaad38)))cycle
do xyzzyaaak38=xyzzyaaad38+1,nitot
if(iontype(xyzzyaaad38)==iontype(xyzzyaaak38))then
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
endif
enddo
xyzzyaaal38=xyzzyaaal38+1
xyzzyaaam38(xyzzyaaal38)=iontype(xyzzyaaad38)
enddo
elseif(xyzzyaaas38)then
if(.not.has_en)return
call xyzzyaaeq1(trule,xyzzyaaab38,xyzzyaaac38,xyzzyaaad38,ion1_match,x&
&yzzyaaah38)
if(xyzzyaaah38/=0)call errstop_master('DIGEST_RULES','Error parsing ru&
&le (assumed grouping rule).')
if(xyzzyaaab38/=0.or.xyzzyaaac38/=0)call errstop_master('DIGEST_RULES'&
&,'Found rule involving electrons when expecting grouping rule.')
if(xyzzyaaad38==0)call errstop_master('DIGEST_RULES','Found rule not i&
&nvolving nuclei when expecting grouping rule.')
if(ion1_match=='n')call errstop_master('DIGEST_RULES','Grouping rules &
&cannot use "n<ion>" syntax; use "N<ion>" or "Z<Z>" instead.')
if(ion1_match=='N')then
do xyzzyaaan38=2,npcells
xyzzyaaak38=xyzzyaaad38+(xyzzyaaan38-1)*nbasis
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
enddo
elseif(ion1_match=='Z')then
do xyzzyaaak38=xyzzyaaad38+1,nitot
if(iontype(xyzzyaaad38)==iontype(xyzzyaaak38))then
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
endif
enddo
endif
elseif(xyzzyaaaq38)then
trule=trim(trule(xyzzyaaaa38+1:))
call xyzzyaaeq1(trule,xyzzyaaab38,xyzzyaaac38,xyzzyaaad38,ion1_match,x&
&yzzyaaah38)
if(xyzzyaaah38/=0)call errstop_master('DIGEST_RULES','Error parsing re&
&moval rule.')
if(xyzzyaaab38==0.and.xyzzyaaac38==0.and.xyzzyaaad38==0)call errstop_m&
&aster('DIGEST_RULES','Empty rule? Bug.')
if(xyzzyaaab38/=0.and.xyzzyaaac38/=0.and.xyzzyaaad38/=0)call errstop_m&
&aster('DIGEST_RULES','Three-way rule? Bug.')
if(xyzzyaaab38>0.and.xyzzyaaac38>0)then
if(.not.has_ee)return
elseif(xyzzyaaab38>0.and.xyzzyaaad38>0)then
if(.not.has_en)return
elseif(xyzzyaaad38>0)then
if(.not.has_en)return
endif
if(xyzzyaaad38/=0)then
if(ion1_match=='N')then
do xyzzyaaan38=2,npcells
xyzzyaaak38=xyzzyaaad38+(xyzzyaaan38-1)*nbasis
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
enddo
elseif(ion1_match=='Z')then
do xyzzyaaak38=xyzzyaaad38+1,nitot
if(iontype(xyzzyaaad38)==iontype(xyzzyaaak38))then
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
endif
enddo
endif
endif
if(xyzzyaaab38>0.and.xyzzyaaac38>0)then
call xyzzyaaer1(groups_ee(xyzzyaaab38,xyzzyaaac38),-1,groups_ee)
elseif(xyzzyaaab38>0.and.xyzzyaaad38>0)then
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaab38),-1,groups_en)
elseif(xyzzyaaab38>0)then
if(has_ee)then
do xyzzyaaaj38=1,nspin
call xyzzyaaer1(groups_ee(xyzzyaaab38,xyzzyaaaj38),-1,groups_ee)
enddo
endif
if(has_en)then
do xyzzyaaak38=1,nitot
call xyzzyaaes1(groups_en(xyzzyaaak38,xyzzyaaab38),-1,groups_en)
enddo
endif
elseif(xyzzyaaad38>0)then
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),-1,groups_en)
enddo
else
call errstop_master('DIGEST_RULES','Rule type not caught. Bug.')
endif
elseif(xyzzyaaar38)then
call xyzzyaaeq1(trule(1:xyzzyaaaa38-1),xyzzyaaab38,xyzzyaaac38,xyzzyaa&
&ad38,ion1_match,xyzzyaaah38)
if(xyzzyaaah38/=0)call errstop_master('DIGEST_RULES','Error parsing fi&
&rst clause of rule.')
if(xyzzyaaab38==0.and.xyzzyaaac38==0.and.xyzzyaaad38==0)call errstop_m&
&aster('DIGEST_RULES','Empty rule? Bug.')
if(xyzzyaaab38/=0.and.xyzzyaaac38/=0.and.xyzzyaaad38/=0)call errstop_m&
&aster('DIGEST_RULES','Three-way rule? Bug.')
if(xyzzyaaab38>0.and.xyzzyaaac38>0)then
if(.not.has_ee)return
elseif(xyzzyaaab38>0.and.xyzzyaaad38>0)then
if(.not.has_en)return
elseif(xyzzyaaad38>0)then
if(.not.has_en)return
endif
trule=trule(xyzzyaaaa38+1:)
if(xyzzyaaad38/=0)then
if(ion1_match=='N')then
do xyzzyaaan38=2,npcells
xyzzyaaak38=xyzzyaaad38+(xyzzyaaan38-1)*nbasis
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
enddo
elseif(ion1_match=='Z')then
do xyzzyaaak38=xyzzyaaad38+1,nitot
if(iontype(xyzzyaaad38)==iontype(xyzzyaaak38))then
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
endif
enddo
endif
endif
xyzzyaaai38=1
do while(len_trim(trule)>0)
xyzzyaaai38=xyzzyaaai38+1
xyzzyaaaa38=index(trule,'=')
if(xyzzyaaaa38<1)xyzzyaaaa38=len_trim(trule)+1
call xyzzyaaeq1(trule(1:xyzzyaaaa38-1),xyzzyaaae38,xyzzyaaaf38,xyzzyaa&
&ag38,ion2_match,xyzzyaaah38)
if(xyzzyaaah38/=0)call errstop_master('DIGEST_RULES','Error parsing cl&
&ause #'//trim(i2s(xyzzyaaai38))//' of rule.')
if((xyzzyaaab38==0.neqv.xyzzyaaae38==0).or.(xyzzyaaac38==0.neqv.xyzzya&
&aaf38==0).or.(xyzzyaaad38==0.neqv.xyzzyaaag38==0))call errstop_master&
&('DIGEST_RULES','Clause #'//trim(i2s(xyzzyaaai38))//' of rule does no&
&t conform with first clause.')
if(xyzzyaaaa38<len_trim(trule))then
trule=trule(xyzzyaaaa38+1:)
else
trule=''
endif
if(xyzzyaaag38/=0)then
if(ion2_match=='N')then
do xyzzyaaan38=2,npcells
xyzzyaaak38=xyzzyaaag38+(xyzzyaaan38-1)*nbasis
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaag38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
enddo
elseif(ion2_match=='Z')then
do xyzzyaaak38=xyzzyaaag38+1,nitot
if(iontype(xyzzyaaag38)==iontype(xyzzyaaak38))then
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaag38,xyzzyaaaj38),groups_en(xyzzyaaak&
&38,xyzzyaaaj38),groups_en)
enddo
endif
enddo
endif
endif
if(xyzzyaaab38>0.and.xyzzyaaac38>0)then
call xyzzyaaer1(groups_ee(xyzzyaaab38,xyzzyaaac38),groups_ee(xyzzyaaae&
&38,xyzzyaaaf38),groups_ee)
elseif(xyzzyaaab38>0.and.xyzzyaaad38>0)then
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaab38),groups_en(xyzzyaaag&
&38,xyzzyaaae38),groups_en)
elseif(xyzzyaaab38>0)then
if(xyzzyaaab38==xyzzyaaae38)cycle
if(has_ee)then
do xyzzyaaaj38=1,nspin
call xyzzyaaer1(groups_ee(xyzzyaaab38,xyzzyaaaj38),groups_ee(xyzzyaaae&
&38,xyzzyaaaj38),groups_ee)
enddo
endif
if(has_en)then
do xyzzyaaak38=1,nitot
call xyzzyaaes1(groups_en(xyzzyaaak38,xyzzyaaab38),groups_en(xyzzyaaak&
&38,xyzzyaaae38),groups_en)
enddo
endif
elseif(xyzzyaaad38>0)then
if(xyzzyaaad38==xyzzyaaag38)cycle
do xyzzyaaaj38=1,nspin
call xyzzyaaes1(groups_en(xyzzyaaad38,xyzzyaaaj38),groups_en(xyzzyaaag&
&38,xyzzyaaaj38),groups_en)
enddo
else
call errstop_master('DIGEST_RULES','Rule type not caught. Bug.')
endif
enddo
endif
end subroutine digest_rule
subroutine xyzzyaaeq1(clause,ispin,jspin,ion,ion_match,ierr)
implicit none
character(*),intent(in) :: clause
integer,intent(out) :: ispin,jspin,ion,ierr
character(1),intent(out) :: ion_match
integer xyzzyaaaa39,xyzzyaaab39
character(20) :: xyzzyaaac39
ispin=0
jspin=0
ion=0
ion_match='n'
ierr=0
xyzzyaaaa39=index(clause,'-')
if(xyzzyaaaa39>0)then
xyzzyaaac39=adjustl(clause(1:xyzzyaaaa39-1))
else
xyzzyaaac39=adjustl(clause)
endif
select case(xyzzyaaac39(1:1))
case('n','N','Z')
ion_match=xyzzyaaac39(1:1)
if(ion_match=='N'.and..not.isperiodic)ion_match='n'
read(xyzzyaaac39(2:),*,iostat=ierr)ion
if(ierr/=0)return
if(trim(xyzzyaaac39(2:))/=trim(i2s(ion)))then
ierr=1
return
endif
select case(ion_match)
case('n')
if(ion<1.or.ion>nitot)then
ierr=1
return
endif
case('N')
if(ion<1.or.ion>nbasis)then
ierr=1
return
endif
case('Z')
xyzzyaaab39=ion
do ion=1,nitot
if(iontype(ion)==xyzzyaaab39)exit
enddo
if(ion>nitot)then
ierr=1
return
endif
end select
case default
read(xyzzyaaac39,*,iostat=ierr)ispin
if(ierr/=0)return
if(trim(xyzzyaaac39)/=trim(i2s(ispin)))then
ierr=1
return
endif
if(ispin<1.or.ispin>nspin)then
ierr=1
return
endif
end select
if(xyzzyaaaa39<1)return
xyzzyaaac39=adjustl(clause(xyzzyaaaa39+1:))
select case(xyzzyaaac39(1:1))
case('n','N','Z')
if(ion/=0)then
ierr=1
return
endif
ion_match=xyzzyaaac39(1:1)
if(ion_match=='N'.and..not.isperiodic)ion_match='n'
read(xyzzyaaac39(2:),*,iostat=ierr)ion
if(ierr/=0)return
if(trim(xyzzyaaac39(2:))/=trim(i2s(ion)))then
ierr=1
return
endif
select case(ion_match)
case('n')
if(ion<1.or.ion>nitot)then
ierr=1
return
endif
case('N')
if(ion<1.or.ion>nbasis)then
ierr=1
return
endif
case('Z')
xyzzyaaab39=ion
do ion=1,nitot
if(iontype(ion)==xyzzyaaab39)exit
enddo
if(ion>nitot)then
ierr=1
return
endif
end select
case default
if(ispin/=0)then
read(xyzzyaaac39,*,iostat=ierr)jspin
if(trim(xyzzyaaac39)/=trim(i2s(jspin)))then
ierr=1
return
endif
if(ierr/=0)return
if(jspin<1.or.jspin>nspin)then
ierr=1
return
endif
if(jspin<ispin)call swap1(ispin,jspin)
else
read(xyzzyaaac39,*,iostat=ierr)ispin
if(ierr/=0)return
if(trim(xyzzyaaac39)/=trim(i2s(ispin)))then
ierr=1
return
endif
if(ispin<1.or.ispin>nspin)then
ierr=1
return
endif
endif
end select
end subroutine xyzzyaaeq1
subroutine xyzzyaaer1(g1,g2,groups_ee)
implicit none
integer,intent(in) :: g1,g2
integer,intent(inout) :: groups_ee(nspin,nspin)
integer xyzzyaaaa40,xyzzyaaab40,xyzzyaaac40,xyzzyaaad40,xyzzyaaae40
if(g1==0.or.g2==0)return
if(g1==g2)return
xyzzyaaaa40=g1
xyzzyaaab40=g2
if(xyzzyaaaa40<xyzzyaaab40)call swap1(xyzzyaaaa40,xyzzyaaab40)
do xyzzyaaac40=1,nspin
xyzzyaaae40=groups_ee(xyzzyaaac40,xyzzyaaac40)
if(xyzzyaaae40==xyzzyaaaa40)then
groups_ee(xyzzyaaac40,xyzzyaaac40)=xyzzyaaab40
elseif(xyzzyaaae40>xyzzyaaaa40)then
groups_ee(xyzzyaaac40,xyzzyaaac40)=xyzzyaaae40-1
endif
do xyzzyaaad40=xyzzyaaac40+1,nspin
xyzzyaaae40=groups_ee(xyzzyaaad40,xyzzyaaac40)
if(xyzzyaaae40==xyzzyaaaa40)then
groups_ee(xyzzyaaad40,xyzzyaaac40)=xyzzyaaab40
groups_ee(xyzzyaaac40,xyzzyaaad40)=xyzzyaaab40
elseif(xyzzyaaae40>xyzzyaaaa40)then
groups_ee(xyzzyaaad40,xyzzyaaac40)=xyzzyaaae40-1
groups_ee(xyzzyaaac40,xyzzyaaad40)=xyzzyaaae40-1
endif
enddo
enddo
end subroutine xyzzyaaer1
subroutine xyzzyaaes1(g1,g2,groups_en)
implicit none
integer,intent(in) :: g1,g2
integer,intent(inout) :: groups_en(nitot,nspin)
integer xyzzyaaaa41,xyzzyaaab41,xyzzyaaac41,xyzzyaaad41,xyzzyaaae41
if(g1==0.or.g2==0)return
if(g1==g2)return
xyzzyaaaa41=g1
xyzzyaaab41=g2
if(xyzzyaaaa41<xyzzyaaab41)call swap1(xyzzyaaaa41,xyzzyaaab41)
do xyzzyaaac41=1,nspin
do xyzzyaaad41=1,nitot
xyzzyaaae41=groups_en(xyzzyaaad41,xyzzyaaac41)
if(xyzzyaaae41==xyzzyaaaa41)then
groups_en(xyzzyaaad41,xyzzyaaac41)=xyzzyaaab41
elseif(xyzzyaaae41>xyzzyaaaa41)then
groups_en(xyzzyaaad41,xyzzyaaac41)=xyzzyaaae41-1
endif
enddo
enddo
end subroutine xyzzyaaes1
subroutine reconstruct_ruleset(has_ee,has_en,groups_ee,groups_en,full_&
&ruleset)
implicit none
integer,intent(in) :: groups_ee(nspin,nspin),groups_en(max(1,nitot),ns&
&pin)
logical,intent(in) :: has_ee,has_en
integer xyzzyaaaa42,xyzzyaaab42,xyzzyaaac42,xyzzyaaad42,xyzzyaaae42,xy&
&zzyaaaf42,xyzzyaaag42,xyzzyaaah42,xyzzyaaai42,xyzzyaaaj42,xyzzyaaak42&
&,xyzzyaaal42(nitot),xyzzyaaam42(nitot),xyzzyaaan42,xyzzyaaao42
integer,allocatable :: xyzzyaaap42(:,:),xyzzyaaaq42(:,:),xyzzyaaar42(:&
&,:),xyzzyaaas42(:,:)
logical xyzzyaaat42,xyzzyaaau42(nitot),xyzzyaaav42(nitot),xyzzyaaaw42(&
&nitot),xyzzyaaax42(nitot),xyzzyaaay42(nitot),xyzzyaaaz42,xyzzyaaba42
character(20) lhs,rhs
character(casl_valsize) trule
character(casl_valsize),pointer :: full_ruleset(:)
allocate(xyzzyaaar42(nspin,nspin),xyzzyaaap42(nspin,nspin),xyzzyaaas42&
&(max(1,nitot),nspin),xyzzyaaaq42(max(1,nitot),nspin),stat=xyzzyaaaa42&
&)
call check_alloc(xyzzyaaaa42,'RECONSTRUCT_RULESET','gcurr_*')
call init_groups_ee(xyzzyaaar42)
call init_groups_en(xyzzyaaas42)
nullify(full_ruleset)
xyzzyaaai42=0
jump_loop: do
if(has_ee.and.has_en)then
if(all(xyzzyaaar42==groups_ee).and.all(xyzzyaaas42==groups_en))exit ju&
&mp_loop
elseif(has_ee)then
if(all(xyzzyaaar42==groups_ee))exit jump_loop
elseif(has_en)then
if(all(xyzzyaaas42==groups_en))exit jump_loop
endif
do xyzzyaaab42=1,nspin
lhs=trim(i2s(xyzzyaaab42))
xyzzyaaat42=.false.
do xyzzyaaac42=xyzzyaaab42+1,nspin
if(has_ee)then
if(any(groups_ee(:,xyzzyaaab42)/=groups_ee(:,xyzzyaaac42).and.groups_e&
&e(:,xyzzyaaab42)/=0.and.groups_ee(:,xyzzyaaac42)/=0))cycle
endif
if(has_en)then
if(any(groups_en(:,xyzzyaaab42)/=groups_en(:,xyzzyaaac42).and.groups_e&
&n(:,xyzzyaaab42)/=0.and.groups_en(:,xyzzyaaac42)/=0))cycle
endif
rhs=trim(i2s(xyzzyaaac42))
trule=trim(lhs)//'='//trim(rhs)
if(has_ee)xyzzyaaap42=xyzzyaaar42
if(has_en)xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaap42,xyzzyaaaq42)
if(has_ee.and.has_en)then
if(all(xyzzyaaap42==xyzzyaaar42).and.all(xyzzyaaaq42==xyzzyaaas42))cyc&
&le
elseif(has_ee)then
if(all(xyzzyaaap42==xyzzyaaar42))cycle
elseif(has_en)then
if(all(xyzzyaaaq42==xyzzyaaas42))cycle
endif
if(has_ee)xyzzyaaar42=xyzzyaaap42
if(has_en)xyzzyaaas42=xyzzyaaaq42
if(.not.xyzzyaaat42)then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaat42=.true.
else
full_ruleset(xyzzyaaai42)=trim(full_ruleset(xyzzyaaai42))//'='//trim(r&
&hs)
endif
if(has_ee.and.any(xyzzyaaar42/=groups_ee))cycle
if(has_en.and.any(xyzzyaaas42/=groups_en))cycle
exit jump_loop
enddo
enddo
if(has_ee)then
if(any(xyzzyaaar42/=groups_ee))then
particle_pair: do xyzzyaaab42=1,nspin
do xyzzyaaac42=xyzzyaaab42,nspin
lhs=trim(i2s(xyzzyaaab42))//'-'//trim(i2s(xyzzyaaac42))
xyzzyaaat42=.false.
do xyzzyaaad42=xyzzyaaab42,nspin
xyzzyaaah42=xyzzyaaad42
if(xyzzyaaab42==xyzzyaaad42)xyzzyaaah42=xyzzyaaac42+1
do xyzzyaaae42=xyzzyaaah42,nspin
if(groups_ee(xyzzyaaab42,xyzzyaaac42)/=groups_ee(xyzzyaaad42,xyzzyaaae&
&42).and.groups_ee(xyzzyaaab42,xyzzyaaac42)/=0.and.groups_ee(xyzzyaaad&
&42,xyzzyaaae42)/=0)cycle
rhs=trim(i2s(xyzzyaaad42))//'-'//trim(i2s(xyzzyaaae42))
trule=trim(lhs)//'='//trim(rhs)
xyzzyaaap42=xyzzyaaar42
call digest_rule(trule,has_ee,has_en,xyzzyaaap42,xyzzyaaas42)
if(all(xyzzyaaap42==xyzzyaaar42))cycle
xyzzyaaar42=xyzzyaaap42
if(.not.xyzzyaaat42)then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaat42=.true.
else
full_ruleset(xyzzyaaai42)=trim(full_ruleset(xyzzyaaai42))//'='//trim(r&
&hs)
endif
if(any(xyzzyaaar42/=groups_ee))cycle
if(has_en.and.any(xyzzyaaas42/=groups_en))exit particle_pair
exit jump_loop
enddo
enddo
enddo
enddo particle_pair
endif
endif
if(has_en)then
xyzzyaaak42=0
xyzzyaaal42=0
do xyzzyaaaf42=1,nitot
do xyzzyaaaj42=1,xyzzyaaak42
if(iontype(xyzzyaaaf42)==xyzzyaaam42(xyzzyaaaj42))exit
enddo
if(xyzzyaaaj42>xyzzyaaak42)then
xyzzyaaak42=xyzzyaaak42+1
xyzzyaaam42(xyzzyaaak42)=iontype(xyzzyaaaf42)
xyzzyaaal42(xyzzyaaak42)=xyzzyaaal42(xyzzyaaak42)+1
else
xyzzyaaal42(xyzzyaaak42)=xyzzyaaal42(xyzzyaaak42)+1
endif
enddo
xyzzyaaay42=.false.
xyzzyaaau42=.false.
do xyzzyaaaf42=1,nitot
do xyzzyaaaj42=1,xyzzyaaak42
if(iontype(xyzzyaaaf42)==xyzzyaaam42(xyzzyaaaj42))exit
enddo
if(xyzzyaaay42(xyzzyaaaj42))cycle
xyzzyaaao42=1
do xyzzyaaag42=xyzzyaaaf42+1,nitot
if(iontype(xyzzyaaag42)==xyzzyaaam42(xyzzyaaaj42))then
if(all(groups_en(xyzzyaaaf42,:)==groups_en(xyzzyaaag42,:)))xyzzyaaao42&
&=xyzzyaaao42+1
endif
enddo
xyzzyaaau42(xyzzyaaaj42)=xyzzyaaao42==xyzzyaaal42(xyzzyaaaj42)
xyzzyaaay42(xyzzyaaaj42)=.true.
enddo
do xyzzyaaaf42=1,nitot
xyzzyaaav42(xyzzyaaaf42)=xyzzyaaau42(iontype(xyzzyaaaf42))
enddo
xyzzyaaaz42=all(xyzzyaaau42(1:xyzzyaaak42))
if(xyzzyaaaz42)then
trule='Z'
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaar42,xyzzyaaaq42)
if(.not.all(xyzzyaaaq42==xyzzyaaas42))then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaas42=xyzzyaaaq42
if(all(xyzzyaaas42==groups_en))then
if(.not.has_ee.or.all(xyzzyaaar42==groups_ee))exit jump_loop
endif
endif
else
zrule_loop: do xyzzyaaaj42=1,xyzzyaaak42
if(.not.xyzzyaaau42(xyzzyaaaj42))cycle
trule='Z'//trim(i2s(xyzzyaaam42(xyzzyaaaj42)))
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaar42,xyzzyaaaq42)
if(.not.all(xyzzyaaaq42==xyzzyaaas42))then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaas42=xyzzyaaaq42
if(any(xyzzyaaas42/=groups_en))cycle
if(has_ee.and.any(xyzzyaaar42/=groups_ee))exit zrule_loop
exit jump_loop
endif
enddo zrule_loop
endif
if(isperiodic.and..not.xyzzyaaaz42)then
xyzzyaaaw42=.false.
do xyzzyaaaf42=1,nbasis
xyzzyaaao42=1
do xyzzyaaan42=2,npcells
xyzzyaaag42=xyzzyaaaf42+(xyzzyaaan42-1)*nbasis
if(all(groups_en(xyzzyaaaf42,:)==groups_en(xyzzyaaag42,:)))xyzzyaaao42&
&=xyzzyaaao42+1
enddo
xyzzyaaaw42(xyzzyaaaf42)=xyzzyaaao42==npcells
enddo
do xyzzyaaaf42=1,nitot
xyzzyaaax42(xyzzyaaaf42)=xyzzyaaaw42(mod(xyzzyaaaf42-1,nbasis)+1)
enddo
xyzzyaaba42=all(xyzzyaaaw42(1:nbasis))
if(xyzzyaaba42)then
trule='N'
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaar42,xyzzyaaaq42)
if(.not.all(xyzzyaaaq42==xyzzyaaas42))then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaas42=xyzzyaaaq42
if(all(xyzzyaaas42==groups_en))then
if(.not.has_ee.or.all(xyzzyaaar42==groups_ee))exit jump_loop
endif
endif
else
nrule_loop: do xyzzyaaaf42=1,nbasis
if(.not.xyzzyaaaw42(xyzzyaaaf42))cycle
trule='N'//trim(i2s(xyzzyaaaf42))
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaar42,xyzzyaaaq42)
if(.not.all(xyzzyaaaq42==xyzzyaaas42))then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaas42=xyzzyaaaq42
if(any(xyzzyaaas42/=groups_en))cycle
if(has_ee.and.any(xyzzyaaar42/=groups_ee))exit nrule_loop
exit jump_loop
endif
enddo nrule_loop
endif
endif
if(any(xyzzyaaas42/=groups_en))then
single_ion: do xyzzyaaaf42=1,nitot
if(xyzzyaaav42(xyzzyaaaf42))then
lhs='Z'//trim(i2s(iontype(xyzzyaaaf42)))
elseif(xyzzyaaax42(xyzzyaaaf42))then
lhs='N'//trim(i2s(mod(xyzzyaaaf42-1,nbasis)+1))
else
lhs='n'//trim(i2s(xyzzyaaaf42))
endif
xyzzyaaat42=.false.
do xyzzyaaag42=xyzzyaaaf42+1,nitot
if(any(groups_en(xyzzyaaaf42,:)/=groups_en(xyzzyaaag42,:).and.groups_e&
&n(xyzzyaaaf42,:)/=0.and.groups_en(xyzzyaaag42,:)/=0))cycle
if(xyzzyaaav42(xyzzyaaag42))then
rhs='Z'//trim(i2s(iontype(xyzzyaaag42)))
elseif(xyzzyaaax42(xyzzyaaag42))then
rhs='N'//trim(i2s(mod(xyzzyaaag42-1,nbasis)+1))
else
rhs='n'//trim(i2s(xyzzyaaag42))
endif
if(lhs==rhs)cycle
trule=trim(lhs)//'='//trim(rhs)
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaar42,xyzzyaaaq42)
if(all(xyzzyaaaq42==xyzzyaaas42))cycle
xyzzyaaas42=xyzzyaaaq42
if(.not.xyzzyaaat42)then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaat42=.true.
else
full_ruleset(xyzzyaaai42)=trim(full_ruleset(xyzzyaaai42))//'='//trim(r&
&hs)
endif
if(any(xyzzyaaas42/=groups_en))cycle
if(has_ee.and.any(xyzzyaaar42/=groups_ee))exit single_ion
exit jump_loop
enddo
enddo single_ion
endif
if(any(xyzzyaaas42/=groups_en))then
particle_ion: do xyzzyaaab42=1,nspin
do xyzzyaaaf42=1,nitot
if(xyzzyaaav42(xyzzyaaaf42))then
lhs=trim(i2s(xyzzyaaab42))//'-Z'//trim(i2s(iontype(xyzzyaaaf42)))
elseif(xyzzyaaax42(xyzzyaaaf42))then
lhs=trim(i2s(xyzzyaaab42))//'-N'//trim(i2s(mod(xyzzyaaaf42-1,nbasis)+1&
&))
else
lhs=trim(i2s(xyzzyaaab42))//'-n'//trim(i2s(xyzzyaaaf42))
endif
xyzzyaaat42=.false.
do xyzzyaaad42=xyzzyaaab42,nspin
xyzzyaaah42=1
if(xyzzyaaab42==xyzzyaaad42)xyzzyaaah42=xyzzyaaaf42+1
do xyzzyaaag42=xyzzyaaah42,nitot
if(groups_en(xyzzyaaaf42,xyzzyaaab42)/=groups_en(xyzzyaaag42,xyzzyaaad&
&42).and.groups_en(xyzzyaaaf42,xyzzyaaab42)/=0.and.groups_en(xyzzyaaag&
&42,xyzzyaaad42)/=0)cycle
if(xyzzyaaav42(xyzzyaaaf42))then
rhs=trim(i2s(xyzzyaaad42))//'-Z'//trim(i2s(iontype(xyzzyaaag42)))
elseif(xyzzyaaax42(xyzzyaaaf42))then
rhs=trim(i2s(xyzzyaaad42))//'-N'//trim(i2s(mod(xyzzyaaag42-1,nbasis)+1&
&))
else
rhs=trim(i2s(xyzzyaaad42))//'-n'//trim(i2s(xyzzyaaag42))
endif
trule=trim(lhs)//'='//trim(rhs)
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaar42,xyzzyaaaq42)
if(all(xyzzyaaaq42==xyzzyaaas42))cycle
xyzzyaaas42=xyzzyaaaq42
if(.not.xyzzyaaat42)then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaat42=.true.
else
full_ruleset(xyzzyaaai42)=trim(full_ruleset(xyzzyaaai42))//'='//trim(r&
&hs)
endif
if(any(xyzzyaaas42/=groups_en))cycle
if(has_ee.and.any(xyzzyaaar42/=groups_ee))exit particle_ion
exit jump_loop
enddo
enddo
enddo
enddo particle_ion
endif
endif
e_removal: do xyzzyaaab42=1,nspin
if(all(groups_ee(:,xyzzyaaab42)<=0).and.any(xyzzyaaar42(:,xyzzyaaab42)&
&>0))then
trule='!'//trim(i2s(xyzzyaaab42))
xyzzyaaap42=xyzzyaaar42
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaap42,xyzzyaaaq42)
if(.not.(all(xyzzyaaap42==xyzzyaaar42).and.all(xyzzyaaaq42==xyzzyaaas4&
&2)))then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaar42=xyzzyaaap42
xyzzyaaas42=xyzzyaaaq42
if(has_ee.and.any(xyzzyaaar42/=groups_ee))cycle
if(any(xyzzyaaas42/=groups_en))exit e_removal
exit jump_loop
endif
endif
enddo e_removal
if(has_ee)then
ee_removal: do xyzzyaaab42=1,nspin
do xyzzyaaad42=xyzzyaaab42,nspin
if(groups_ee(xyzzyaaab42,xyzzyaaad42)<=0.and.xyzzyaaar42(xyzzyaaad42,x&
&yzzyaaab42)>0)then
trule='!'//trim(i2s(xyzzyaaab42))//'-'//trim(i2s(xyzzyaaad42))
xyzzyaaap42=xyzzyaaar42
call digest_rule(trule,has_ee,has_en,xyzzyaaap42,xyzzyaaas42)
if(.not.all(xyzzyaaap42==xyzzyaaar42))then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaar42=xyzzyaaap42
if(any(xyzzyaaar42/=groups_ee))cycle
if(has_en.and.any(xyzzyaaas42/=groups_en))exit ee_removal
exit jump_loop
endif
endif
enddo
enddo ee_removal
endif
if(has_en)then
n_removal: do xyzzyaaaf42=1,nitot
if(all(groups_en(xyzzyaaaf42,:)<=0).and.any(xyzzyaaas42(xyzzyaaaf42,:)&
&>0))then
if(xyzzyaaav42(xyzzyaaaf42))then
trule='!Z'//trim(i2s(iontype(xyzzyaaaf42)))
elseif(xyzzyaaax42(xyzzyaaaf42))then
trule='!N'//trim(i2s(mod(xyzzyaaaf42-1,nbasis)+1))
else
trule='!n'//trim(i2s(xyzzyaaaf42))
endif
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaar42,xyzzyaaaq42)
if(.not.all(xyzzyaaaq42==xyzzyaaas42))then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaas42=xyzzyaaaq42
if(any(xyzzyaaas42/=groups_en))cycle
if(has_ee.and.any(xyzzyaaar42/=groups_ee))exit n_removal
exit jump_loop
endif
endif
enddo n_removal
en_removal: do xyzzyaaaf42=1,nitot
do xyzzyaaab42=1,nspin
if(groups_en(xyzzyaaaf42,xyzzyaaab42)<=0.and.xyzzyaaas42(xyzzyaaaf42,x&
&yzzyaaab42)>0)then
if(xyzzyaaav42(xyzzyaaaf42))then
trule='!'//trim(i2s(xyzzyaaab42))//'-Z'//trim(i2s(iontype(xyzzyaaaf42)&
&))
elseif(xyzzyaaax42(xyzzyaaaf42))then
trule='!'//trim(i2s(xyzzyaaab42))//'-N'//trim(i2s(mod(xyzzyaaaf42-1,nb&
&asis)+1))
else
trule='!'//trim(i2s(xyzzyaaab42))//'-n'//trim(i2s(xyzzyaaaf42))
endif
xyzzyaaaq42=xyzzyaaas42
call digest_rule(trule,has_ee,has_en,xyzzyaaar42,xyzzyaaaq42)
if(.not.all(xyzzyaaaq42==xyzzyaaas42))then
xyzzyaaai42=xyzzyaaai42+1
call resize_pointer(casl_valsize,(/xyzzyaaai42/),full_ruleset)
full_ruleset(xyzzyaaai42)=trule
xyzzyaaas42=xyzzyaaaq42
if(any(xyzzyaaas42/=groups_en))cycle
if(has_ee.and.any(xyzzyaaar42/=groups_ee))exit en_removal
exit jump_loop
endif
endif
enddo
enddo en_removal
endif
call wout('Current ruleset:')
if(xyzzyaaai42==0)then
call wout('  (empty)')
else
do xyzzyaaao42=1,xyzzyaaai42
call wout('  '//trim(full_ruleset(xyzzyaaao42)))
enddo
endif
call wout()
call wout('Target groups_ee:')
do xyzzyaaab42=1,nspin
call wout(' ',groups_ee(:,xyzzyaaab42))
enddo
call wout('Current groups_ee:')
do xyzzyaaab42=1,nspin
call wout(' ',xyzzyaaar42(:,xyzzyaaab42))
enddo
call wout()
call wout('Target groups_en:')
do xyzzyaaab42=1,nspin
call wout(' ',groups_en(:,xyzzyaaab42))
enddo
call wout('Current groups_en:')
do xyzzyaaab42=1,nspin
call wout(' ',xyzzyaaas42(:,xyzzyaaab42))
enddo
call wout()
call errstop_master('RECONSTRUCT_RULESET','Unable to reconstruct the t&
&arget groups. This is a bug.')
enddo jump_loop
deallocate(xyzzyaaar42,xyzzyaaap42,xyzzyaaas42,xyzzyaaaq42)
end subroutine reconstruct_ruleset
subroutine gen_sensible_ruleset(ruleset)
implicit none
character(casl_valsize),pointer :: ruleset(:)
integer xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43,xyzzyaaad43,xyzzyaaae43,xy&
&zzyaaaf43,xyzzyaaag43,xyzzyaaah43,xyzzyaaai43,xyzzyaaaj43
logical xyzzyaaak43
nullify(ruleset)
xyzzyaaah43=0
if(isperiodic.and.nitot>0)then
xyzzyaaah43=xyzzyaaah43+1
call resize_pointer(casl_valsize,(/xyzzyaaah43/),ruleset)
ruleset(xyzzyaaah43)='N'
endif
do xyzzyaaaa43=1,no_families
do xyzzyaaab43=1,nspin
if(which_fam(xyzzyaaab43)/=xyzzyaaaa43)cycle
do xyzzyaaac43=xyzzyaaab43,nspin
if(which_fam(xyzzyaaac43)/=xyzzyaaaa43.or.nele(xyzzyaaab43)/=nele(xyzz&
&yaaac43).or.nele(xyzzyaaab43)==0)cycle
do xyzzyaaad43=xyzzyaaaa43,no_families
xyzzyaaai43=1
if(xyzzyaaaa43==xyzzyaaad43)xyzzyaaai43=xyzzyaaab43
do xyzzyaaae43=xyzzyaaai43,nspin
if(which_fam(xyzzyaaae43)/=xyzzyaaad43)cycle
xyzzyaaaj43=xyzzyaaae43
if(xyzzyaaae43==xyzzyaaab43)xyzzyaaaj43=xyzzyaaac43+1
do xyzzyaaaf43=xyzzyaaaj43,nspin
if(which_fam(xyzzyaaaf43)/=xyzzyaaad43.or.nele(xyzzyaaae43)/=nele(xyzz&
&yaaaf43).or.nele(xyzzyaaae43)==0)cycle
xyzzyaaak43=.true.
xyzzyaaak43=xyzzyaaak43.and.(xyzzyaaaa43/=xyzzyaaad43.or.(xyzzyaaab43=&
&=xyzzyaaae43.eqv.xyzzyaaac43==xyzzyaaaf43))
xyzzyaaak43=xyzzyaaak43.and..not.(xyzzyaaab43==xyzzyaaae43.and.nele(xy&
&zzyaaab43)<2).and..not.(xyzzyaaac43==xyzzyaaaf43.and.nele(xyzzyaaac43&
&)<2)
xyzzyaaak43=xyzzyaaak43.and..not.(xyzzyaaab43==xyzzyaaac43.and.xyzzyaa&
&ae43==xyzzyaaaf43)
if(xyzzyaaak43)then
xyzzyaaah43=xyzzyaaah43+1
call resize_pointer(casl_valsize,(/xyzzyaaah43/),ruleset)
ruleset(xyzzyaaah43)=trim(i2s(xyzzyaaab43))//'-'//trim(i2s(xyzzyaaae43&
&))//'='//trim(i2s(xyzzyaaac43))//'-'//trim(i2s(xyzzyaaaf43))
endif
xyzzyaaak43=.true.
xyzzyaaak43=xyzzyaaak43.and.(xyzzyaaab43==xyzzyaaac43.eqv.xyzzyaaae43=&
&=xyzzyaaaf43)
xyzzyaaak43=xyzzyaaak43.and.((xyzzyaaaa43/=xyzzyaaad43.and.which_eqvfa&
&m(xyzzyaaab43)==which_eqvfam(xyzzyaaae43)).or.(xyzzyaaaa43==xyzzyaaad&
&43.and..not.xyzzyaaac43==xyzzyaaae43))
xyzzyaaak43=xyzzyaaak43.and.nele(xyzzyaaab43)==nele(xyzzyaaae43)
xyzzyaaak43=xyzzyaaak43.and..not.(xyzzyaaab43==xyzzyaaac43.and.nele(xy&
&zzyaaab43)<2).and..not.(xyzzyaaae43==xyzzyaaaf43.and.nele(xyzzyaaae43&
&)<2)
xyzzyaaak43=xyzzyaaak43.and..not.(xyzzyaaab43==xyzzyaaae43.and.xyzzyaa&
&ac43==xyzzyaaaf43)
if(xyzzyaaak43)then
xyzzyaaah43=xyzzyaaah43+1
call resize_pointer(casl_valsize,(/xyzzyaaah43/),ruleset)
ruleset(xyzzyaaah43)=trim(i2s(xyzzyaaab43))//'-'//trim(i2s(xyzzyaaac43&
&))//'='//trim(i2s(xyzzyaaae43))//'-'//trim(i2s(xyzzyaaaf43))
endif
enddo
enddo
enddo
if(xyzzyaaab43/=xyzzyaaac43)then
do xyzzyaaag43=1,nitot
xyzzyaaah43=xyzzyaaah43+1
call resize_pointer(casl_valsize,(/xyzzyaaah43/),ruleset)
ruleset(xyzzyaaah43)=trim(i2s(xyzzyaaab43))//'-n'//trim(i2s(xyzzyaaag4&
&3))//'='//trim(i2s(xyzzyaaac43))//'-n'//trim(i2s(xyzzyaaag43))
enddo
endif
enddo
enddo
enddo
end subroutine gen_sensible_ruleset
subroutine get_sig_from_model(model,groups_ee,groups_en,rank_e,rank_n,&
&signature,ierr,integer_model)
implicit none
integer,intent(in) :: rank_e,rank_n,groups_ee(nspin,nspin),groups_en(n&
&itot,nspin)
integer,intent(out) :: signature(:),ierr
integer,intent(out),optional :: integer_model(rank_e+rank_n)
character(*),intent(in) :: model
integer xyzzyaaaa44(rank_e),xyzzyaaab44(rank_n)
signature=0
call parse_model(model,rank_e,rank_n,xyzzyaaaa44,xyzzyaaab44,ierr)
if(ierr/=0)return
call get_sig(rank_e,rank_n,groups_ee,groups_en,xyzzyaaaa44,xyzzyaaab44&
&,signature)
if(present(integer_model))then
integer_model=0
integer_model(1:rank_e)=xyzzyaaaa44
if(rank_n>0)integer_model(rank_e+1:rank_e+rank_n)=xyzzyaaab44
endif
end subroutine get_sig_from_model
subroutine parse_model(model,rank_e,rank_n,ispin_vector,ion_vector,ier&
&r)
implicit none
integer,intent(in) :: rank_e,rank_n
integer,intent(out) :: ispin_vector(rank_e),ion_vector(rank_n),ierr
character(*),intent(in) :: model
integer xyzzyaaaa45,xyzzyaaab45,xyzzyaaac45,xyzzyaaad45,xyzzyaaae45,xy&
&zzyaaaf45
character(20) :: xyzzyaaag45
character(casl_keysize) tmodel
ierr=0
xyzzyaaae45=0
xyzzyaaaf45=0
ispin_vector=0
if(rank_n>0)ion_vector=0
tmodel=model
do xyzzyaaaa45=1,rank_e+rank_n
xyzzyaaab45=index(tmodel,'-')
if(xyzzyaaab45>0)then
xyzzyaaag45=adjustl(tmodel(1:xyzzyaaab45-1))
tmodel=trim(adjustl(tmodel(xyzzyaaab45+1:)))
else
xyzzyaaag45=adjustl(tmodel)
tmodel=''
endif
xyzzyaaad45=0
read(xyzzyaaag45,*,iostat=ierr)xyzzyaaac45
if(ierr/=0)then
if(xyzzyaaag45(1:1)/='n')return
xyzzyaaac45=0
read(xyzzyaaag45(2:),*,iostat=ierr)xyzzyaaad45
if(ierr/=0)return
endif
ierr=-2
if(xyzzyaaac45>0)then
xyzzyaaae45=xyzzyaaae45+1
if(xyzzyaaae45>rank_e.or.xyzzyaaac45>nspin)return
ispin_vector(xyzzyaaae45)=xyzzyaaac45
elseif(xyzzyaaad45>0)then
xyzzyaaaf45=xyzzyaaaf45+1
if(xyzzyaaaf45>rank_n.or.xyzzyaaad45>nitot)return
ion_vector(xyzzyaaaf45)=xyzzyaaad45
else
ierr=-1
return
endif
enddo
ierr=0
end subroutine parse_model
character(casl_keysize) function model_string(rank_e,rank_n,model)
implicit none
integer,intent(in) :: rank_e,rank_n,model(:)
integer xyzzyaaaa46
model_string=trim(i2s(model(1)))
do xyzzyaaaa46=2,rank_e
model_string=trim(model_string)//'-'//trim(i2s(model(xyzzyaaaa46)))
enddo
do xyzzyaaaa46=1,rank_n
model_string=trim(model_string)//'-n'//trim(i2s(model(rank_e+xyzzyaaaa&
&46)))
enddo
end function model_string
subroutine build_channels(groups_ee,groups_en,rank_e,rank_n,nchannel,c&
&hannel_sig,channel_model)
implicit none
integer,intent(in) :: rank_e,rank_n,groups_ee(nspin,nspin),groups_en(n&
&itot,nspin)
integer,intent(out) :: nchannel
integer xyzzyaaaa47(rank_e),xyzzyaaab47(rank_n),xyzzyaaac47,xyzzyaaad4&
&7,xyzzyaaae47,xyzzyaaaf47,xyzzyaaag47
integer,pointer :: channel_sig(:,:),channel_model(:,:)
integer,allocatable :: xyzzyaaah47(:)
logical xyzzyaaai47,xyzzyaaaj47
type signature
integer,pointer :: sig(:)=>null(),model(:)=>null()
type(signature),pointer :: next=>null()
end type signature
type(signature),pointer :: xyzzyaaak47,xyzzyaaal47,sig
xyzzyaaad47=(rank_e*(rank_e-1))/2+rank_e*rank_n
xyzzyaaae47=rank_e+rank_n
allocate(xyzzyaaah47(xyzzyaaad47),stat=xyzzyaaac47)
call check_alloc(xyzzyaaac47,'BUILD_CHANNELS','sig_full')
xyzzyaaah47=0
nullify(xyzzyaaak47)
nchannel=0
xyzzyaaaa47=0
do while(xyzzyaaew1(rank_e,nspin,xyzzyaaaa47))
xyzzyaaaj47=.false.
do xyzzyaaag47=1,nspin
if(count(xyzzyaaaa47==xyzzyaaag47)>nele(xyzzyaaag47))then
xyzzyaaaj47=.true.
exit
endif
enddo
if(xyzzyaaaj47)cycle
if(rank_n>0)xyzzyaaab47=0
do
if(rank_n>0)then
if(.not.xyzzyaaev1(rank_n,nitot,xyzzyaaab47))exit
endif
call get_sig(rank_e,rank_n,groups_ee,groups_en,xyzzyaaaa47,xyzzyaaab47&
&,xyzzyaaah47)
if(all(xyzzyaaah47>0))then
if(associated(xyzzyaaak47))then
xyzzyaaal47=>xyzzyaaak47
do
xyzzyaaai47=all(xyzzyaaah47==xyzzyaaal47%sig)
if(xyzzyaaai47)exit
if(.not.associated(xyzzyaaal47%next))then
nchannel=nchannel+1
allocate(sig,stat=xyzzyaaac47)
call check_alloc(xyzzyaaac47,'BUILD_CHANNELS','sig')
allocate(sig%sig(xyzzyaaad47),sig%model(xyzzyaaae47),stat=xyzzyaaac47)
call check_alloc(xyzzyaaac47,'BUILD_CHANNELS','sig%*')
sig%sig=xyzzyaaah47
sig%model=(/xyzzyaaaa47,xyzzyaaab47/)
xyzzyaaal47%next=>sig
exit
endif
xyzzyaaal47=>xyzzyaaal47%next
enddo
else
nchannel=nchannel+1
allocate(xyzzyaaak47,stat=xyzzyaaac47)
call check_alloc(xyzzyaaac47,'BUILD_CHANNELS','sig1')
allocate(xyzzyaaak47%sig(xyzzyaaad47),xyzzyaaak47%model(xyzzyaaae47),s&
&tat=xyzzyaaac47)
call check_alloc(xyzzyaaac47,'BUILD_CHANNELS','sig1%*')
xyzzyaaak47%sig=xyzzyaaah47
xyzzyaaak47%model=(/xyzzyaaaa47,xyzzyaaab47/)
endif
endif
if(rank_n==0)exit
enddo
enddo
deallocate(xyzzyaaah47)
allocate(channel_sig(xyzzyaaad47,nchannel),channel_model(xyzzyaaae47,n&
&channel),stat=xyzzyaaac47)
call check_alloc(xyzzyaaac47,'BUILD_CHANNEL','channel_sig')
do xyzzyaaaf47=1,nchannel
channel_sig(1:xyzzyaaad47,xyzzyaaaf47)=xyzzyaaak47%sig(1:xyzzyaaad47)
channel_model(:,xyzzyaaaf47)=xyzzyaaak47%model(:)
deallocate(xyzzyaaak47%sig,xyzzyaaak47%model)
if(xyzzyaaaf47<nchannel)sig=>xyzzyaaak47%next
deallocate(xyzzyaaak47)
if(xyzzyaaaf47<nchannel)xyzzyaaak47=>sig
enddo
end subroutine build_channels
subroutine get_sig(rank_e,rank_n,groups_ee,groups_en,ispin_vector,ion_&
&vector,sig,perm_ee,perm_en)
implicit none
integer,intent(in) :: rank_e,rank_n,groups_ee(nspin,nspin),groups_en(n&
&itot,nspin),ispin_vector(rank_e),ion_vector(rank_n)
integer,intent(inout) :: sig(*)
integer,intent(inout),optional :: perm_ee(*),perm_en(*)
integer,parameter :: xyzzyaaaa48=1024
integer xyzzyaaab48,xyzzyaaac48,xyzzyaaad48,xyzzyaaae48,xyzzyaaaf48,xy&
&zzyaaag48(2,rank_e),xyzzyaaah48(rank_e,rank_e),xyzzyaaai48(rank_e,ran&
&k_n),xyzzyaaaj48(rank_e,xyzzyaaaa48),xyzzyaaak48(rank_n,xyzzyaaaa48)
do xyzzyaaab48=1,rank_e
xyzzyaaah48(xyzzyaaab48,xyzzyaaab48)=0
do xyzzyaaac48=xyzzyaaab48+1,rank_e
xyzzyaaah48(xyzzyaaac48,xyzzyaaab48)=groups_ee(ispin_vector(xyzzyaaac4&
&8),ispin_vector(xyzzyaaab48))
xyzzyaaah48(xyzzyaaab48,xyzzyaaac48)=xyzzyaaah48(xyzzyaaac48,xyzzyaaab&
&48)
enddo
enddo
do xyzzyaaab48=1,rank_e
do xyzzyaaac48=1,rank_n
xyzzyaaai48(xyzzyaaab48,xyzzyaaac48)=groups_en(ion_vector(xyzzyaaac48)&
&,ispin_vector(xyzzyaaab48))
enddo
enddo
xyzzyaaae48=1
xyzzyaaaf48=1
xyzzyaaag48(1:2,1)=(/1,rank_e/)
do xyzzyaaab48=1,rank_e
xyzzyaaaj48(xyzzyaaab48,1)=xyzzyaaab48
enddo
call sort_matrix_symm(rank_e,xyzzyaaae48,xyzzyaaah48,xyzzyaaaj48,xyzzy&
&aaaf48,xyzzyaaag48)
do xyzzyaaad48=1,xyzzyaaae48
do xyzzyaaab48=1,rank_n
xyzzyaaak48(xyzzyaaab48,xyzzyaaad48)=xyzzyaaab48
enddo
enddo
call sort_matrix_rect(rank_e,rank_n,xyzzyaaae48,xyzzyaaai48,xyzzyaaaj4&
&8,xyzzyaaak48,xyzzyaaaf48,xyzzyaaag48)
call xyzzyaaeu1(rank_e,rank_n,sig,xyzzyaaaj48,xyzzyaaak48,xyzzyaaah48,&
&xyzzyaaai48,perm_ee,perm_en)
end subroutine get_sig
subroutine get_sig_ee_only(rank_e,groups_ee,ispin_vector,sig,perm_ee)
implicit none
integer,intent(in) :: rank_e,groups_ee(nspin,nspin),ispin_vector(rank_&
&e)
integer,intent(inout) :: sig(*)
integer,intent(inout),optional :: perm_ee(*)
integer,parameter :: xyzzyaaaa49=1024
integer xyzzyaaab49,xyzzyaaac49,xyzzyaaad49,xyzzyaaae49,xyzzyaaaf49(2,&
&rank_e),ee_matrix(rank_e,rank_e),xyzzyaaag49(rank_e,xyzzyaaaa49)
do xyzzyaaab49=1,rank_e
ee_matrix(xyzzyaaab49,xyzzyaaab49)=0
do xyzzyaaac49=xyzzyaaab49+1,rank_e
ee_matrix(xyzzyaaac49,xyzzyaaab49)=groups_ee(ispin_vector(xyzzyaaac49)&
&,ispin_vector(xyzzyaaab49))
ee_matrix(xyzzyaaab49,xyzzyaaac49)=ee_matrix(xyzzyaaac49,xyzzyaaab49)
enddo
enddo
xyzzyaaad49=1
xyzzyaaae49=1
xyzzyaaaf49(1:2,1)=(/1,rank_e/)
do xyzzyaaab49=1,rank_e
xyzzyaaag49(xyzzyaaab49,1)=xyzzyaaab49
enddo
call sort_matrix_symm(rank_e,xyzzyaaad49,ee_matrix,xyzzyaaag49,xyzzyaa&
&ae49,xyzzyaaaf49)
call xyzzyaaeu1(rank_e,0,sig,xyzzyaaag49,ee_matrix=ee_matrix,perm_ee=p&
&erm_ee)
end subroutine get_sig_ee_only
subroutine get_sig_en_only(rank_n,groups_en,ispin,ion_vector,sig,perm_&
&en)
implicit none
integer,intent(in) :: rank_n,groups_en(nitot,nspin),ispin,ion_vector(r&
&ank_n)
integer,intent(inout) :: sig(*)
integer,intent(inout),optional :: perm_en(*)
integer,parameter :: xyzzyaaaa50=1024
integer xyzzyaaab50,xyzzyaaac50,xyzzyaaad50,xyzzyaaae50(2,1),en_matrix&
&(1,rank_n),xyzzyaaaf50(1,xyzzyaaaa50),xyzzyaaag50(rank_n,xyzzyaaaa50)
do xyzzyaaab50=1,rank_n
en_matrix(1,xyzzyaaab50)=groups_en(ion_vector(xyzzyaaab50),ispin)
enddo
xyzzyaaac50=1
xyzzyaaad50=1
xyzzyaaae50(1:2,1)=(/1,1/)
xyzzyaaaf50(1,1)=1
do xyzzyaaab50=1,rank_n
xyzzyaaag50(xyzzyaaab50,1)=xyzzyaaab50
enddo
call sort_matrix_rect(1,rank_n,xyzzyaaac50,en_matrix,xyzzyaaaf50,xyzzy&
&aaag50,xyzzyaaad50,xyzzyaaae50)
call xyzzyaaeu1(1,rank_n,sig,xyzzyaaaf50,xyzzyaaag50,en_matrix=en_matr&
&ix,perm_en=perm_en)
end subroutine get_sig_en_only
subroutine sort_indices_matrix(rank_e,rank_n,size_sig,symm_eebasis,sig&
&,expindx,nperm,sig_perm,perm_common_sign,perm_sign_consistent)
implicit none
integer,intent(in) :: rank_e,rank_n,size_sig,symm_eebasis,sig(size_sig&
&)
integer,intent(inout) :: expindx(size_sig),nperm,sig_perm(size_sig,*),&
&perm_common_sign
logical,intent(inout) :: perm_sign_consistent
integer,parameter :: xyzzyaaaa51=1024
integer xyzzyaaab51,xyzzyaaac51,xyzzyaaad51,xyzzyaaae51,xyzzyaaaf51(ra&
&nk_e,rank_e),xyzzyaaag51(rank_e,rank_n),xyzzyaaah51,xyzzyaaai51(2,ran&
&k_e),xyzzyaaaj51(max(1,rank_e),xyzzyaaaa51),xyzzyaaak51(max(1,rank_n)&
&,xyzzyaaaa51)
logical xyzzyaaal51
call expindx2matrices(rank_e,rank_n,expindx,xyzzyaaaf51,xyzzyaaag51)
xyzzyaaad51=1
xyzzyaaah51=1
xyzzyaaai51(1:2,1)=(/1,rank_e/)
do xyzzyaaab51=1,rank_e
xyzzyaaaj51(xyzzyaaab51,1)=xyzzyaaab51
enddo
call sort_matrix_symm(rank_e,xyzzyaaad51,xyzzyaaaf51,xyzzyaaaj51,xyzzy&
&aaah51,xyzzyaaai51)
do xyzzyaaac51=1,xyzzyaaad51
do xyzzyaaab51=1,rank_n
xyzzyaaak51(xyzzyaaab51,xyzzyaaac51)=xyzzyaaab51
enddo
enddo
call sort_matrix_rect(rank_e,rank_n,xyzzyaaad51,xyzzyaaag51,xyzzyaaaj5&
&1,xyzzyaaak51,xyzzyaaah51,xyzzyaaai51)
call xyzzyaaeu1(rank_e,rank_n,expindx,xyzzyaaaj51,xyzzyaaak51,xyzzyaaa&
&f51,xyzzyaaag51)
call amb2cand_sort_matrix_rect(rank_e,rank_n,xyzzyaaad51,xyzzyaaaj51,x&
&yzzyaaak51,xyzzyaaah51,xyzzyaaai51)
if(rank_e>1)then
perm_common_sign=0
xyzzyaaal51=.false.
perm_sign_consistent=.true.
do xyzzyaaac51=1,xyzzyaaad51
call xyzzyaaet1(rank_e,xyzzyaaaj51(1,xyzzyaaac51),xyzzyaaaf51,symm_eeb&
&asis,xyzzyaaae51)
if(xyzzyaaae51==0)cycle
if(xyzzyaaal51)then
if(xyzzyaaae51==-perm_common_sign)then
perm_common_sign=0
perm_sign_consistent=.false.
exit
endif
else
perm_common_sign=xyzzyaaae51
xyzzyaaal51=.true.
endif
enddo
else
perm_common_sign=1
perm_sign_consistent=.true.
endif
nperm=xyzzyaaad51
call expindx2matrices(rank_e,rank_n,sig,xyzzyaaaf51,xyzzyaaag51)
do xyzzyaaac51=1,xyzzyaaad51
call xyzzyaaeu1(rank_e,rank_n,sig_perm(1,xyzzyaaac51),xyzzyaaaj51(1,xy&
&zzyaaac51),xyzzyaaak51(1,xyzzyaaac51),xyzzyaaaf51,xyzzyaaag51)
enddo
end subroutine sort_indices_matrix
subroutine xyzzyaaet1(rank_e,e_indx,ee_matrix,symm_ee,perm_sign)
implicit none
integer,intent(in) :: rank_e,e_indx(rank_e),ee_matrix(rank_e,rank_e),s&
&ymm_ee
integer,intent(out) :: perm_sign
integer xyzzyaaaa52,xyzzyaaab52,xyzzyaaac52,xyzzyaaad52
if(rank_e==1)then
perm_sign=1
return
endif
select case(symm_ee)
case(-1)
perm_sign=1
do xyzzyaaaa52=1,rank_e
xyzzyaaac52=e_indx(xyzzyaaaa52)
do xyzzyaaab52=xyzzyaaaa52+1,rank_e
xyzzyaaad52=e_indx(xyzzyaaab52)
if(xyzzyaaad52<xyzzyaaac52)then
if(ee_matrix(xyzzyaaad52,xyzzyaaac52)/=0)perm_sign=-perm_sign
endif
enddo
enddo
case(0)
perm_sign=0
case(1)
perm_sign=1
end select
end subroutine xyzzyaaet1
subroutine insertion_point(rank_e,ii_vector,ii_loc)
implicit none
integer,intent(in) :: rank_e,ii_vector(rank_e)
integer,intent(out) :: ii_loc
integer xyzzyaaaa53
xyzzyaaaa53=ii_vector(1)
do ii_loc=1,rank_e-1
if(xyzzyaaaa53<ii_vector(ii_loc+1))exit
enddo
end subroutine insertion_point
subroutine apply_insertion_sign(ii_loc,perm_ee,expindx,p)
implicit none
integer,intent(in) :: ii_loc,perm_ee(*),expindx(*)
real(dp),intent(inout) :: p
integer xyzzyaaaa54
do xyzzyaaaa54=1,ii_loc-1
if(expindx(perm_ee(xyzzyaaaa54))/=0)p=-p
enddo
end subroutine apply_insertion_sign
subroutine expindx2matrices(rank_e,rank_n,expindx,ee_matrix,en_matrix)
implicit none
integer,intent(in) :: rank_e,rank_n,expindx(*)
integer,intent(inout) :: ee_matrix(rank_e,rank_e),en_matrix(rank_e,ran&
&k_n)
integer xyzzyaaaa55,xyzzyaaab55,xyzzyaaac55,xyzzyaaad55
xyzzyaaaa55=0
do xyzzyaaac55=1,rank_e
ee_matrix(xyzzyaaac55,xyzzyaaac55)=0
do xyzzyaaad55=xyzzyaaac55+1,rank_e
xyzzyaaaa55=xyzzyaaaa55+1
ee_matrix(xyzzyaaad55,xyzzyaaac55)=expindx(xyzzyaaaa55)
ee_matrix(xyzzyaaac55,xyzzyaaad55)=expindx(xyzzyaaaa55)
enddo
enddo
xyzzyaaab55=0
do xyzzyaaac55=1,rank_e
do xyzzyaaad55=1,rank_n
xyzzyaaab55=xyzzyaaab55+1
en_matrix(xyzzyaaac55,xyzzyaaad55)=expindx(xyzzyaaaa55+xyzzyaaab55)
enddo
enddo
end subroutine expindx2matrices
subroutine xyzzyaaeu1(rank_e,rank_n,expindx,e_indx,n_indx,ee_matrix,en&
&_matrix,perm_ee,perm_en)
implicit none
integer,intent(in) :: rank_e,rank_n,e_indx(rank_e)
integer,intent(in),optional :: n_indx(rank_n),ee_matrix(rank_e,rank_e)&
&,en_matrix(rank_e,rank_n)
integer,intent(inout) :: expindx(*)
integer,intent(inout),optional :: perm_ee(*),perm_en(*)
integer xyzzyaaaa56,xyzzyaaab56,xyzzyaaac56,xyzzyaaad56
integer xyzzyaaae56(rank_e,rank_e),xyzzyaaaf56(rank_e,rank_n)
xyzzyaaaa56=0
do xyzzyaaac56=1,rank_e
do xyzzyaaad56=xyzzyaaac56+1,rank_e
xyzzyaaaa56=xyzzyaaaa56+1
expindx(xyzzyaaaa56)=ee_matrix(e_indx(xyzzyaaad56),e_indx(xyzzyaaac56)&
&)
enddo
enddo
xyzzyaaab56=0
do xyzzyaaac56=1,rank_e
do xyzzyaaad56=1,rank_n
xyzzyaaab56=xyzzyaaab56+1
expindx(xyzzyaaaa56+xyzzyaaab56)=en_matrix(e_indx(xyzzyaaac56),n_indx(&
&xyzzyaaad56))
enddo
enddo
if(present(perm_ee))then
xyzzyaaaa56=0
do xyzzyaaac56=1,rank_e
do xyzzyaaad56=xyzzyaaac56+1,rank_e
xyzzyaaaa56=xyzzyaaaa56+1
xyzzyaaae56(xyzzyaaad56,xyzzyaaac56)=xyzzyaaaa56
xyzzyaaae56(xyzzyaaac56,xyzzyaaad56)=xyzzyaaaa56
enddo
enddo
xyzzyaaaa56=0
do xyzzyaaac56=1,rank_e
do xyzzyaaad56=xyzzyaaac56+1,rank_e
xyzzyaaaa56=xyzzyaaaa56+1
perm_ee(xyzzyaaaa56)=xyzzyaaae56(e_indx(xyzzyaaad56),e_indx(xyzzyaaac5&
&6))
enddo
enddo
endif
if(present(perm_en))then
xyzzyaaab56=0
do xyzzyaaac56=1,rank_e
do xyzzyaaad56=1,rank_n
xyzzyaaab56=xyzzyaaab56+1
xyzzyaaaf56(xyzzyaaac56,xyzzyaaad56)=xyzzyaaab56
enddo
enddo
xyzzyaaab56=0
do xyzzyaaac56=1,rank_e
do xyzzyaaad56=1,rank_n
xyzzyaaab56=xyzzyaaab56+1
perm_en(xyzzyaaab56)=xyzzyaaaf56(e_indx(xyzzyaaac56),n_indx(xyzzyaaad5&
&6))
enddo
enddo
endif
end subroutine xyzzyaaeu1
subroutine match_signature(size_sig,max_size_sig,nchannel,signature,ch&
&annel_sig,ichannel)
implicit none
integer,intent(in) :: size_sig,max_size_sig,nchannel,signature(size_si&
&g),channel_sig(max_size_sig,nchannel)
integer,intent(out) :: ichannel
if(size_sig==1)then
ichannel=signature(1)
elseif(nchannel==1)then
ichannel=1
else
do ichannel=1,nchannel
if(all(signature(1:size_sig)==channel_sig(1:size_sig,ichannel)))return
enddo
ichannel=0
endif
end subroutine match_signature
logical function iterate_nuclei_indices_fix(rank_e,rank_n,ne,nn,ii_fix&
&,nzencut1,nzencut,ignore_nzencut,nzmap_e,nzmap_n,nzmax_e,nzmax_n,nzin&
&dex_n,ion_vector)
implicit none
integer,intent(in) :: rank_e,rank_n,ne,nn,ii_fix
integer,intent(inout) :: nzmap_e(ne,rank_n),nzmap_n(nn,rank_n),nzmax_e&
&(rank_n),nzmax_n(rank_n),nzindex_n(rank_n),ion_vector(rank_n)
logical,intent(in) :: nzencut1(nn),nzencut(nn,*),ignore_nzencut
integer xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58,xyzzyaaae58
iterate_nuclei_indices_fix=.true.
if(ion_vector(1)==0)then
xyzzyaaac58=0
do xyzzyaaab58=1,nn
if(.not.nzencut1(xyzzyaaab58).and..not.ignore_nzencut)cycle
xyzzyaaac58=xyzzyaaac58+1
nzmap_n(xyzzyaaac58,1)=xyzzyaaab58
enddo
nzmax_n(1)=xyzzyaaac58
if(xyzzyaaac58==0)then
iterate_nuclei_indices_fix=.false.
return
endif
xyzzyaaaa58=1
nzindex_n(1)=0
else
xyzzyaaaa58=rank_n
endif
do
xyzzyaaac58=nzindex_n(xyzzyaaaa58)+1
if(xyzzyaaac58-xyzzyaaaa58+rank_n>nzmax_n(xyzzyaaaa58))then
if(xyzzyaaaa58==1)then
iterate_nuclei_indices_fix=.false.
return
endif
xyzzyaaaa58=xyzzyaaaa58-1
cycle
endif
nzindex_n(xyzzyaaaa58)=xyzzyaaac58
xyzzyaaad58=nzmap_n(xyzzyaaac58,xyzzyaaaa58)
ion_vector(xyzzyaaaa58)=xyzzyaaad58
if(rank_e>1)then
xyzzyaaac58=0
if(xyzzyaaaa58==1)then
xyzzyaaac58=xyzzyaaac58+1
nzmap_e(xyzzyaaac58,xyzzyaaaa58)=ii_fix
do xyzzyaaab58=1,ne
if(xyzzyaaab58==ii_fix)cycle
if(.not.nzencut(xyzzyaaad58,xyzzyaaab58).and..not.ignore_nzencut)cycle
xyzzyaaac58=xyzzyaaac58+1
nzmap_e(xyzzyaaac58,xyzzyaaaa58)=xyzzyaaab58
enddo
else
xyzzyaaac58=xyzzyaaac58+1
nzmap_e(xyzzyaaac58,xyzzyaaaa58)=ii_fix
do xyzzyaaab58=2,nzmax_e(xyzzyaaaa58-1)
xyzzyaaae58=nzmap_e(xyzzyaaab58,xyzzyaaaa58-1)
if(.not.nzencut(xyzzyaaad58,xyzzyaaae58).and..not.ignore_nzencut)cycle
xyzzyaaac58=xyzzyaaac58+1
nzmap_e(xyzzyaaac58,xyzzyaaaa58)=xyzzyaaae58
enddo
endif
if(xyzzyaaac58<rank_e)cycle
nzmax_e(xyzzyaaaa58)=xyzzyaaac58
endif
if(xyzzyaaaa58==rank_n)return
xyzzyaaac58=0
do xyzzyaaab58=nzindex_n(xyzzyaaaa58)+1,nzmax_n(xyzzyaaaa58)
xyzzyaaac58=xyzzyaaac58+1
nzmap_n(xyzzyaaac58,xyzzyaaaa58+1)=nzmap_n(xyzzyaaab58,xyzzyaaaa58)
enddo
nzmax_n(xyzzyaaaa58+1)=xyzzyaaac58
xyzzyaaaa58=xyzzyaaaa58+1
nzindex_n(xyzzyaaaa58)=0
enddo
end function iterate_nuclei_indices_fix
logical function iterate_nuclei_indices(rank_e,rank_n,ne,nn,nzencut,ig&
&nore_nzencut,nzmap_e,nzmap_n,nzmax_e,nzmax_n,nzindex_n,ion_vector)
implicit none
integer,intent(in) :: rank_e,rank_n,ne,nn
integer,intent(inout) :: nzmap_e(ne,rank_n),nzmap_n(nn,rank_n),nzmax_e&
&(rank_n),nzmax_n(rank_n),nzindex_n(rank_n),ion_vector(rank_n)
logical,intent(in) :: nzencut(nn,*),ignore_nzencut
integer xyzzyaaaa59,xyzzyaaab59,xyzzyaaac59,xyzzyaaad59,xyzzyaaae59
iterate_nuclei_indices=.true.
if(ion_vector(1)==0)then
xyzzyaaac59=0
do xyzzyaaab59=1,nn
xyzzyaaac59=xyzzyaaac59+1
nzmap_n(xyzzyaaac59,1)=xyzzyaaab59
enddo
nzmax_n(1)=xyzzyaaac59
if(xyzzyaaac59==0)then
iterate_nuclei_indices=.false.
return
endif
xyzzyaaaa59=1
nzindex_n(1)=0
else
xyzzyaaaa59=rank_n
endif
do
xyzzyaaac59=nzindex_n(xyzzyaaaa59)+1
if(xyzzyaaac59-xyzzyaaaa59+rank_n>nzmax_n(xyzzyaaaa59))then
if(xyzzyaaaa59==1)then
iterate_nuclei_indices=.false.
return
endif
xyzzyaaaa59=xyzzyaaaa59-1
cycle
endif
nzindex_n(xyzzyaaaa59)=xyzzyaaac59
xyzzyaaad59=nzmap_n(xyzzyaaac59,xyzzyaaaa59)
ion_vector(xyzzyaaaa59)=xyzzyaaad59
xyzzyaaac59=0
if(xyzzyaaaa59==1)then
do xyzzyaaab59=1,ne
if(.not.nzencut(xyzzyaaad59,xyzzyaaab59).and..not.ignore_nzencut)cycle
xyzzyaaac59=xyzzyaaac59+1
nzmap_e(xyzzyaaac59,xyzzyaaaa59)=xyzzyaaab59
enddo
else
do xyzzyaaab59=1,nzmax_e(xyzzyaaaa59-1)
xyzzyaaae59=nzmap_e(xyzzyaaab59,xyzzyaaaa59-1)
if(.not.nzencut(xyzzyaaad59,xyzzyaaae59).and..not.ignore_nzencut)cycle
xyzzyaaac59=xyzzyaaac59+1
nzmap_e(xyzzyaaac59,xyzzyaaaa59)=xyzzyaaae59
enddo
endif
if(xyzzyaaac59<rank_e)cycle
nzmax_e(xyzzyaaaa59)=xyzzyaaac59
if(xyzzyaaaa59==rank_n)return
xyzzyaaac59=0
do xyzzyaaab59=nzindex_n(xyzzyaaaa59)+1,nzmax_n(xyzzyaaaa59)
xyzzyaaac59=xyzzyaaac59+1
nzmap_n(xyzzyaaac59,xyzzyaaaa59+1)=nzmap_n(xyzzyaaab59,xyzzyaaaa59)
enddo
nzmax_n(xyzzyaaaa59+1)=xyzzyaaac59
xyzzyaaaa59=xyzzyaaaa59+1
nzindex_n(xyzzyaaaa59)=0
enddo
end function iterate_nuclei_indices
logical function iterate_electron_indices_fix(rank_e,ne,ii_fix,nzeecut&
&1,nzeecut,ignore_nzeecut,nzmap_e,nzmax_e,nzindex_e,ii_vector,preinit,&
&preinit_max)
implicit none
integer,intent(in) :: rank_e,ne,ii_fix
integer,intent(in),optional :: preinit(ne),preinit_max
integer,intent(inout) :: nzmap_e(ne,rank_e),nzmax_e(rank_e),nzindex_e(&
&rank_e),ii_vector(rank_e)
logical,intent(in) :: nzeecut1(ne),nzeecut(ne,*),ignore_nzeecut
integer xyzzyaaaa60,xyzzyaaab60,xyzzyaaac60,xyzzyaaad60,xyzzyaaae60
iterate_electron_indices_fix=.true.
if(ii_vector(1)==0)then
nzmap_e(1,1)=ii_fix
nzindex_e(1)=1
nzmax_e(1)=1
ii_vector(1)=ii_fix
if(rank_e==1)return
if(present(preinit))then
xyzzyaaac60=0
do xyzzyaaab60=1,preinit_max
xyzzyaaae60=preinit(xyzzyaaab60)
if(xyzzyaaae60==ii_fix)cycle
if(.not.nzeecut1(xyzzyaaae60).and..not.ignore_nzeecut)cycle
xyzzyaaac60=xyzzyaaac60+1
nzmap_e(xyzzyaaac60,2)=xyzzyaaae60
enddo
nzmax_e(2)=xyzzyaaac60
else
xyzzyaaac60=0
do xyzzyaaae60=1,ne
if(xyzzyaaae60==ii_fix)cycle
if(.not.nzeecut1(xyzzyaaae60).and..not.ignore_nzeecut)cycle
xyzzyaaac60=xyzzyaaac60+1
nzmap_e(xyzzyaaac60,2)=xyzzyaaae60
enddo
nzmax_e(2)=xyzzyaaac60
endif
if(nzmax_e(2)==0)then
iterate_electron_indices_fix=.false.
return
endif
xyzzyaaaa60=2
nzindex_e(2)=0
else
xyzzyaaaa60=rank_e
endif
do
xyzzyaaac60=nzindex_e(xyzzyaaaa60)+1
if(xyzzyaaac60-xyzzyaaaa60+rank_e>nzmax_e(xyzzyaaaa60))then
if(xyzzyaaaa60==1)then
iterate_electron_indices_fix=.false.
return
endif
xyzzyaaaa60=xyzzyaaaa60-1
cycle
endif
nzindex_e(xyzzyaaaa60)=xyzzyaaac60
xyzzyaaad60=nzmap_e(xyzzyaaac60,xyzzyaaaa60)
ii_vector(xyzzyaaaa60)=xyzzyaaad60
if(xyzzyaaaa60==rank_e)return
xyzzyaaac60=0
do xyzzyaaab60=nzindex_e(xyzzyaaaa60)+1,nzmax_e(xyzzyaaaa60)
xyzzyaaae60=nzmap_e(xyzzyaaab60,xyzzyaaaa60)
if(.not.nzeecut(xyzzyaaae60,xyzzyaaad60).and..not.ignore_nzeecut)cycle
xyzzyaaac60=xyzzyaaac60+1
nzmap_e(xyzzyaaac60,xyzzyaaaa60+1)=xyzzyaaae60
enddo
nzmax_e(xyzzyaaaa60+1)=xyzzyaaac60
xyzzyaaaa60=xyzzyaaaa60+1
nzindex_e(xyzzyaaaa60)=0
enddo
end function iterate_electron_indices_fix
logical function iterate_electron_indices(rank_e,ne,nzeecut,ignore_nze&
&ecut,nzmap_e,nzmax_e,nzindex_e,ii_vector,preinit,preinit_max)
implicit none
integer,intent(in) :: rank_e,ne
integer,intent(in),optional :: preinit(ne),preinit_max
integer,intent(inout) :: nzmap_e(ne,rank_e),nzmax_e(rank_e),nzindex_e(&
&rank_e),ii_vector(rank_e)
logical,intent(in) :: nzeecut(ne,*),ignore_nzeecut
integer xyzzyaaaa61,xyzzyaaab61,xyzzyaaac61,xyzzyaaad61,xyzzyaaae61
iterate_electron_indices=.true.
if(ii_vector(1)==0)then
if(present(preinit))then
nzmap_e(1:preinit_max,1)=preinit(1:preinit_max)
nzmax_e(1)=preinit_max
else
xyzzyaaac61=0
do xyzzyaaab61=1,ne
xyzzyaaac61=xyzzyaaac61+1
nzmap_e(xyzzyaaac61,1)=xyzzyaaab61
enddo
nzmax_e(1)=xyzzyaaac61
endif
if(nzmax_e(1)==0)then
iterate_electron_indices=.false.
return
endif
xyzzyaaaa61=1
nzindex_e(1)=0
else
xyzzyaaaa61=rank_e
endif
do
xyzzyaaac61=nzindex_e(xyzzyaaaa61)+1
if(xyzzyaaac61-xyzzyaaaa61+rank_e>nzmax_e(xyzzyaaaa61))then
if(xyzzyaaaa61==1)then
iterate_electron_indices=.false.
return
endif
xyzzyaaaa61=xyzzyaaaa61-1
cycle
endif
nzindex_e(xyzzyaaaa61)=xyzzyaaac61
xyzzyaaad61=nzmap_e(xyzzyaaac61,xyzzyaaaa61)
ii_vector(xyzzyaaaa61)=xyzzyaaad61
if(xyzzyaaaa61==rank_e)return
xyzzyaaac61=0
do xyzzyaaab61=nzindex_e(xyzzyaaaa61)+1,nzmax_e(xyzzyaaaa61)
xyzzyaaae61=nzmap_e(xyzzyaaab61,xyzzyaaaa61)
if(.not.nzeecut(xyzzyaaae61,xyzzyaaad61).and..not.ignore_nzeecut)cycle
xyzzyaaac61=xyzzyaaac61+1
nzmap_e(xyzzyaaac61,xyzzyaaaa61+1)=xyzzyaaae61
enddo
nzmax_e(xyzzyaaaa61+1)=xyzzyaaac61
xyzzyaaaa61=xyzzyaaaa61+1
nzindex_e(xyzzyaaaa61)=0
enddo
end function iterate_electron_indices
logical function xyzzyaaev1(n,order,indices)
integer,intent(in) :: n,order
integer,intent(inout) :: indices(n)
integer xyzzyaaaa62,xyzzyaaab62,xyzzyaaac62
xyzzyaaev1=.true.
if(indices(1)==0)then
if(n>order)then
xyzzyaaev1=.false.
return
endif
do xyzzyaaaa62=1,n
indices(xyzzyaaaa62)=xyzzyaaaa62
enddo
return
endif
do xyzzyaaaa62=n,1,-1
xyzzyaaac62=indices(xyzzyaaaa62)+1
if(xyzzyaaac62+n-xyzzyaaaa62>order)cycle
indices(xyzzyaaaa62)=xyzzyaaac62
do xyzzyaaab62=xyzzyaaaa62+1,n
xyzzyaaac62=xyzzyaaac62+1
indices(xyzzyaaab62)=xyzzyaaac62
enddo
return
enddo
xyzzyaaev1=.false.
end function xyzzyaaev1
logical function xyzzyaaew1(n,order,indices)
integer,intent(in) :: n,order
integer,intent(inout) :: indices(n)
integer xyzzyaaaa63,xyzzyaaab63
xyzzyaaew1=.true.
if(indices(1)==0)then
if(order<1)then
xyzzyaaew1=.false.
return
endif
indices=1
return
endif
do xyzzyaaaa63=n,1,-1
xyzzyaaab63=indices(xyzzyaaaa63)+1
if(xyzzyaaab63>order)cycle
indices(xyzzyaaaa63:n)=xyzzyaaab63
return
enddo
xyzzyaaew1=.false.
end function xyzzyaaew1
subroutine group_by_sig(size_ee,size_en,size_sig,sig,ngroup,length_gro&
&up)
implicit none
integer,intent(in) :: size_ee,size_en,size_sig,sig(size_sig)
integer,intent(out) :: ngroup,length_group(size_sig)
integer xyzzyaaaa64,xyzzyaaab64
ngroup=0
length_group=0
xyzzyaaab64=0
do xyzzyaaaa64=1,size_ee
if(sig(xyzzyaaaa64)>xyzzyaaab64)then
ngroup=ngroup+1
length_group(ngroup)=1
xyzzyaaab64=sig(xyzzyaaaa64)
else
length_group(ngroup)=length_group(ngroup)+1
endif
enddo
xyzzyaaab64=0
do xyzzyaaaa64=size_ee+1,size_ee+size_en
if(sig(xyzzyaaaa64)>xyzzyaaab64)then
ngroup=ngroup+1
length_group(ngroup)=1
xyzzyaaab64=sig(xyzzyaaaa64)
else
length_group(ngroup)=length_group(ngroup)+1
endif
enddo
end subroutine group_by_sig
subroutine sort_indices_grouped_ge(ngroup,length_group,indices)
implicit none
integer,intent(in) :: ngroup,length_group(ngroup)
integer,intent(inout) :: indices(*)
integer xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65,xyzzyaaad65,xyzzyaaae65
xyzzyaaac65=0
do xyzzyaaaa65=1,ngroup
xyzzyaaab65=xyzzyaaac65+1
xyzzyaaac65=xyzzyaaac65+length_group(xyzzyaaaa65)
do xyzzyaaad65=xyzzyaaab65,xyzzyaaac65
do xyzzyaaae65=xyzzyaaad65+1,xyzzyaaac65
if(indices(xyzzyaaae65)<indices(xyzzyaaad65))call swap1(indices(xyzzya&
&aad65),indices(xyzzyaaae65))
enddo
enddo
enddo
end subroutine sort_indices_grouped_ge
logical function xyzzyaaex1(n,nn,mask,is_first)
implicit none
integer,intent(in) :: n,nn
logical,intent(inout) :: mask(n),is_first
integer xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66,xyzzyaaad66
logical xyzzyaaae66
xyzzyaaex1=.true.
xyzzyaaae66=n<1
if(is_first)then
is_first=.false.
if(.not.xyzzyaaae66)then
xyzzyaaaa66=n-nn
if(xyzzyaaaa66>0)mask(1:xyzzyaaaa66)=.false.
if(xyzzyaaaa66<n)mask(xyzzyaaaa66+1:n)=.true.
endif
return
endif
if(.not.xyzzyaaae66.and.nn>0)then
xyzzyaaab66=nn
xyzzyaaac66=0
do xyzzyaaad66=n,2,-1
if(mask(xyzzyaaad66))then
if(.not.mask(xyzzyaaad66-1))then
xyzzyaaac66=xyzzyaaad66-1
exit
endif
xyzzyaaab66=xyzzyaaab66-1
endif
enddo
if(xyzzyaaac66>0)then
mask(xyzzyaaac66)=.true.
xyzzyaaaa66=n-nn+xyzzyaaab66
if(xyzzyaaac66<xyzzyaaaa66)mask(xyzzyaaac66+1:xyzzyaaaa66)=.false.
if(xyzzyaaaa66<n)mask(xyzzyaaaa66+1:n)=.true.
return
endif
endif
xyzzyaaex1=.false.
end function xyzzyaaex1
logical function xyzzyaaey1(n,order,indices,is_first)
implicit none
integer,intent(in) :: n,order
integer,intent(inout) :: indices(n)
logical,intent(inout) :: is_first
integer xyzzyaaaa67,xyzzyaaab67
xyzzyaaey1=.true.
if(is_first)then
is_first=.false.
indices(1:n)=1
return
endif
do xyzzyaaaa67=n,1,-1
xyzzyaaab67=indices(xyzzyaaaa67)+1
if(xyzzyaaab67>order)cycle
indices(xyzzyaaaa67)=xyzzyaaab67
if(xyzzyaaaa67<n)indices(xyzzyaaaa67+1:n)=1
return
enddo
xyzzyaaey1=.false.
end function xyzzyaaey1
subroutine generate_all_indices(rank_e,rank_n,size_ee,size_en,size_sig&
&,order_ee,order_en,is_dot_product_ee,is_dot_product_en,no_repeat_ee,n&
&o_repeat_en,maxsum,all_masks,which_ee_pair,which_en_pair,nparam,index&
&_list,param_name)
implicit none
integer,intent(in) :: rank_e,rank_n,size_ee,size_en,size_sig,order_ee,&
&order_en,maxsum,which_ee_pair(2,size_ee),which_en_pair(2,size_en)
integer,intent(inout),optional :: nparam,index_list(size_sig,*)
logical,intent(in) :: is_dot_product_ee,is_dot_product_en,no_repeat_ee&
&,no_repeat_en,all_masks
character(*),optional :: param_name(:)
integer xyzzyaaaa68,xyzzyaaab68,xyzzyaaac68,xyzzyaaad68,xyzzyaaae68,xy&
&zzyaaaf68,xyzzyaaag68,xyzzyaaah68(size_ee),xyzzyaaai68(size_en),xyzzy&
&aaaj68(size_ee),xyzzyaaak68(size_en),xyzzyaaal68,xyzzyaaam68,xyzzyaaa&
&n68,xyzzyaaao68,xyzzyaaap68,xyzzyaaaq68
logical xyzzyaaar68,xyzzyaaas68(size_ee),xyzzyaaat68(size_en),xyzzyaaa&
&u68(rank_e),xyzzyaaav68(rank_n),xyzzyaaaw68,xyzzyaaax68,xyzzyaaay68,x&
&yzzyaaaz68
xyzzyaaal68=size_ee
xyzzyaaaf68=1
if(is_dot_product_ee)then
xyzzyaaal68=2
xyzzyaaaf68=dimensionality
endif
xyzzyaaao68=size_en
xyzzyaaag68=1
if(is_dot_product_en)then
xyzzyaaao68=2
xyzzyaaag68=dimensionality
endif
if(all_masks)then
if(xyzzyaaal68/=size_ee.or.xyzzyaaao68/=size_en)call errstop_master('G&
&ENERATE_ALL_INDICES','Unsupported combination of "All masks" and "Dot&
& product".')
xyzzyaaam68=min(1,size_ee)
xyzzyaaan68=size_ee
xyzzyaaap68=min(1,size_en)
xyzzyaaaq68=size_en
else
xyzzyaaam68=xyzzyaaal68
xyzzyaaan68=xyzzyaaal68
xyzzyaaap68=xyzzyaaao68
xyzzyaaaq68=xyzzyaaao68
endif
if(maxsum>0)then
xyzzyaaae68=maxsum
else
xyzzyaaae68=0
endif
xyzzyaaad68=0
do xyzzyaaal68=xyzzyaaam68,xyzzyaaan68
xyzzyaaaw68=.true.
do while(xyzzyaaex1(size_ee,xyzzyaaal68,xyzzyaaas68,xyzzyaaaw68))
if(no_repeat_ee)then
xyzzyaaar68=.true.
xyzzyaaau68=.false.
do xyzzyaaaa68=1,size_ee
if(.not.xyzzyaaas68(xyzzyaaaa68))cycle
xyzzyaaab68=which_ee_pair(1,xyzzyaaaa68)
xyzzyaaac68=which_ee_pair(2,xyzzyaaaa68)
if(xyzzyaaau68(xyzzyaaab68).or.xyzzyaaau68(xyzzyaaac68))then
xyzzyaaar68=.false.
exit
endif
xyzzyaaau68(xyzzyaaab68)=.true.
xyzzyaaau68(xyzzyaaac68)=.true.
enddo
if(.not.xyzzyaaar68)cycle
endif
xyzzyaaay68=.true.
do while(xyzzyaaey1(xyzzyaaal68,order_ee,xyzzyaaaj68,xyzzyaaay68))
if(xyzzyaaaf68>1)then
if(any(modulo(xyzzyaaaj68(2:xyzzyaaal68),xyzzyaaaf68)/=modulo(xyzzyaaa&
&j68(1),xyzzyaaaf68)))cycle
endif
xyzzyaaah68=unpack(xyzzyaaaj68(1:xyzzyaaal68),xyzzyaaas68,0)
do xyzzyaaao68=xyzzyaaap68,xyzzyaaaq68
xyzzyaaax68=.true.
do while(xyzzyaaex1(size_en,xyzzyaaao68,xyzzyaaat68,xyzzyaaax68))
if(no_repeat_en)then
xyzzyaaar68=.true.
xyzzyaaau68=.false.
xyzzyaaav68=.false.
do xyzzyaaaa68=1,size_en
if(.not.xyzzyaaat68(xyzzyaaaa68))cycle
xyzzyaaab68=which_en_pair(1,xyzzyaaaa68)
xyzzyaaac68=which_en_pair(2,xyzzyaaaa68)
if(xyzzyaaau68(xyzzyaaab68).or.xyzzyaaav68(xyzzyaaac68))then
xyzzyaaar68=.false.
exit
endif
xyzzyaaau68(xyzzyaaab68)=.true.
xyzzyaaav68(xyzzyaaac68)=.true.
enddo
if(.not.xyzzyaaar68)cycle
endif
xyzzyaaaz68=.true.
do while(xyzzyaaey1(xyzzyaaao68,order_en,xyzzyaaak68,xyzzyaaaz68))
if(xyzzyaaag68>1)then
if(any(modulo(xyzzyaaak68(2:xyzzyaaao68),xyzzyaaag68)/=modulo(xyzzyaaa&
&k68(1),xyzzyaaag68)))cycle
endif
xyzzyaaai68=unpack(xyzzyaaak68(1:xyzzyaaao68),xyzzyaaat68,0)
if(xyzzyaaae68>0.and.sum(xyzzyaaah68)+sum(xyzzyaaai68)>xyzzyaaae68)cyc&
&le
xyzzyaaad68=xyzzyaaad68+1
if(present(index_list))then
index_list(1:size_sig,xyzzyaaad68)=(/xyzzyaaah68,xyzzyaaai68/)
if(present(param_name))then
param_name(xyzzyaaad68)='c_'//trim(i2s(index_list(1,xyzzyaaad68)))
do xyzzyaaaa68=2,size_sig
param_name(xyzzyaaad68)=trim(param_name(xyzzyaaad68))//','//trim(i2s(i&
&ndex_list(xyzzyaaaa68,xyzzyaaad68)))
enddo
endif
endif
enddo
enddo
enddo
enddo
enddo
enddo
if(present(nparam))nparam=xyzzyaaad68
end subroutine generate_all_indices
subroutine get_all_index_properties(size_ee,size_en,size_sig,nparam,in&
&dex_list,param_iszero,nparam_offset,nparam_remap,index_remap,leftmost&
&_change_ee,leftmost_change_en,mask_leftmost_change_ee,mask_leftmost_c&
&hange_en)
implicit none
integer,intent(in) :: size_ee,size_en,size_sig,nparam,index_list(size_&
&sig,nparam),nparam_offset
integer,intent(inout) :: nparam_remap,index_remap(nparam),leftmost_cha&
&nge_ee(nparam),leftmost_change_en(nparam),mask_leftmost_change_ee(npa&
&ram),mask_leftmost_change_en(nparam)
logical,intent(in) :: param_iszero(nparam)
integer xyzzyaaaa69,xyzzyaaab69(size_sig),xyzzyaaac69(size_sig),xyzzya&
&aad69,xyzzyaaae69,xyzzyaaaf69,xyzzyaaag69,xyzzyaaah69
do xyzzyaaad69=1,max(size_ee,size_en)
xyzzyaaac69(xyzzyaaad69)=xyzzyaaad69
enddo
xyzzyaaae69=1
xyzzyaaaf69=size_ee
xyzzyaaag69=size_ee+1
xyzzyaaah69=size_sig
nparam_remap=0
do xyzzyaaaa69=1,nparam
if(param_iszero(xyzzyaaaa69))cycle
nparam_remap=nparam_remap+1
index_remap(nparam_remap)=xyzzyaaaa69+nparam_offset
if(size_ee==0)then
leftmost_change_ee(xyzzyaaaa69)=0
mask_leftmost_change_ee(xyzzyaaaa69)=0
elseif(nparam_remap==1)then
leftmost_change_ee(xyzzyaaaa69)=1
mask_leftmost_change_ee(xyzzyaaaa69)=1
else
leftmost_change_ee(xyzzyaaaa69)=sum(minloc(xyzzyaaac69(1:size_ee),inde&
&x_list(xyzzyaaae69:xyzzyaaaf69,xyzzyaaaa69)/=xyzzyaaab69(xyzzyaaae69:&
&xyzzyaaaf69)))
mask_leftmost_change_ee(xyzzyaaaa69)=sum(minloc(xyzzyaaac69(1:size_ee)&
&,index_list(xyzzyaaae69:xyzzyaaaf69,xyzzyaaaa69)<1.neqv.xyzzyaaab69(x&
&yzzyaaae69:xyzzyaaaf69)<1))
endif
if(size_en==0)then
leftmost_change_en(xyzzyaaaa69)=0
mask_leftmost_change_en(xyzzyaaaa69)=0
elseif(nparam_remap==1)then
leftmost_change_en(xyzzyaaaa69)=1
mask_leftmost_change_en(xyzzyaaaa69)=1
else
leftmost_change_en(xyzzyaaaa69)=sum(minloc(xyzzyaaac69(1:size_en),inde&
&x_list(xyzzyaaag69:xyzzyaaah69,xyzzyaaaa69)/=xyzzyaaab69(xyzzyaaag69:&
&xyzzyaaah69)))
mask_leftmost_change_en(xyzzyaaaa69)=sum(minloc(xyzzyaaac69(1:size_en)&
&,index_list(xyzzyaaag69:xyzzyaaah69,xyzzyaaaa69)<1.neqv.xyzzyaaab69(x&
&yzzyaaag69:xyzzyaaah69)<1))
endif
xyzzyaaab69(1:size_sig)=index_list(1:size_sig,xyzzyaaaa69)
enddo
end subroutine get_all_index_properties
subroutine print_eqns(nparam,neqn,matrix,rhs,pname,neqn_actual,eqn_piv&
&ot,eqn_solves_for)
implicit none
integer,intent(in) :: nparam,neqn
integer,intent(in),optional :: neqn_actual,eqn_pivot(neqn),eqn_solves_&
&for(neqn)
real(dp),intent(in) :: matrix(nparam,neqn),rhs(neqn)
character(*),intent(in),optional :: pname(nparam)
integer xyzzyaaaa70,xyzzyaaab70,xyzzyaaac70,xyzzyaaad70,xyzzyaaae70
character(20) char_20
character(4096) text
if(present(neqn_actual))then
xyzzyaaac70=neqn_actual
else
xyzzyaaac70=neqn
endif
call wout('Printing '//trim(i2s(xyzzyaaac70))//' equations for '//trim&
&(i2s(nparam))//' parameters.')
do xyzzyaaaa70=1,xyzzyaaac70
if(present(neqn_actual))then
xyzzyaaad70=eqn_pivot(xyzzyaaaa70)
else
xyzzyaaad70=xyzzyaaaa70
endif
if(present(eqn_solves_for))then
xyzzyaaae70=eqn_solves_for(xyzzyaaaa70)
else
xyzzyaaae70=1
endif
text=''
do xyzzyaaab70=xyzzyaaae70,nparam
if(matrix(xyzzyaaab70,xyzzyaaad70)==0.d0)cycle
if(matrix(xyzzyaaab70,xyzzyaaad70)>0.d0)then
write(char_20,'(es16.8)')matrix(xyzzyaaab70,xyzzyaaad70)
if(len_trim(text)==0)then
if(present(pname))then
text=trim(adjustl(char_20))//'*'//trim(pname(xyzzyaaab70))
else
text=trim(adjustl(char_20))//'*[#'//trim(i2s(xyzzyaaab70))//']'
endif
else
if(present(pname))then
text=trim(text)//' + '//trim(adjustl(char_20))//'*'//trim(pname(xyzzya&
&aab70))
else
text=trim(text)//' + '//trim(adjustl(char_20))//'*[#'//trim(i2s(xyzzya&
&aab70))//']'
endif
endif
else
write(char_20,'(es16.8)')-matrix(xyzzyaaab70,xyzzyaaad70)
if(len_trim(text)==0)then
if(present(pname))then
text='-'//trim(adjustl(char_20))//'*'//trim(pname(xyzzyaaab70))
else
text='-'//trim(adjustl(char_20))//'*[#'//trim(i2s(xyzzyaaab70))//']'
endif
else
if(present(pname))then
text=trim(text)//' - '//trim(adjustl(char_20))//'*'//trim(pname(xyzzya&
&aab70))
else
text=trim(text)//' - '//trim(adjustl(char_20))//'*[#'//trim(i2s(xyzzya&
&aab70))//']'
endif
endif
endif
enddo
if(len_trim(text)==0)then
if(rhs(xyzzyaaad70)==0.d9)cycle
text='0'
endif
write(char_20,'(es16.8)')rhs(xyzzyaaad70)
text=trim(text)//' = '//trim(adjustl(char_20))
call wordwrap('#'//trim(i2s(xyzzyaaaa70))//': '//trim(text))
enddo
end subroutine print_eqns
subroutine xyzzyaaez1(string,indx,is_isotropic,periodic_only,require_c&
&utoff_periodic,is_cusp_friendly,f_zero,is_none)
implicit none
integer,intent(out) :: indx
logical,intent(out) :: is_isotropic,periodic_only,require_cutoff_perio&
&dic,is_cusp_friendly,is_none
character(8),intent(out) :: f_zero
character(*),intent(in) :: string
integer xyzzyaaaa71
indx=-1
select case(string)
case('basis:none','cutoff:none')
indx=0
case default
do xyzzyaaaa71=1,size(xyzzyaaav1)
if(unique_casl_string(trim(string))==unique_casl_string(trim(xyzzyaaav&
&1(xyzzyaaaa71))))then
indx=xyzzyaaaa71
exit
endif
enddo
end select
is_isotropic=.true.
periodic_only=.false.
require_cutoff_periodic=.false.
is_cusp_friendly=.true.
is_none=.false.
f_zero='r'
select case(indx)
case(-1)
continue
case(0)
is_none=.true.
case(xyzzyaaaa1)
require_cutoff_periodic=.true.
case(xyzzyaaab1)
is_isotropic=.false.
periodic_only=.true.
case(xyzzyaaac1)
is_isotropic=.false.
periodic_only=.true.
case(xyzzyaaad1)
require_cutoff_periodic=.true.
case(xyzzyaaae1)
require_cutoff_periodic=.true.
case(xyzzyaaaf1)
require_cutoff_periodic=.true.
case(xyzzyaaag1)
require_cutoff_periodic=.true.
is_cusp_friendly=.false.
case(xyzzyaaah1)
require_cutoff_periodic=.true.
case(xyzzyaaai1)
require_cutoff_periodic=.true.
is_cusp_friendly=.false.
case(xyzzyaaaj1)
require_cutoff_periodic=.true.
case(xyzzyaaak1)
require_cutoff_periodic=.true.
f_zero='r2log'
case(xyzzyaaal1)
require_cutoff_periodic=.true.
f_zero='1/sqrtr'
case(xyzzyaaam1)
require_cutoff_periodic=.true.
case(xyzzyaaan1)
require_cutoff_periodic=.true.
f_zero='1/sqrtr'
case(xyzzyaaao1)
is_isotropic=.false.
periodic_only=.true.
case(xyzzyaaap1)
continue
case(xyzzyaaaq1)
continue
case(xyzzyaaar1)
continue
case(xyzzyaaas1)
is_isotropic=.false.
is_cusp_friendly=.false.
case(xyzzyaaat1)
f_zero='quasi'
case(xyzzyaaau1)
continue
end select
end subroutine xyzzyaaez1
subroutine xyzzyaafa1(label,ibasis,iorder,nconst_int,const_int_name,co&
&nst_int,nconst_dble,const_dble_name,const_dble,nconst_int_calc,const_&
&int_calc,nconst_dble_calc,const_dble_calc,which_unity,basis_symmetry)
implicit none
integer,intent(in) :: ibasis,iorder
integer,intent(out) :: nconst_int,nconst_dble,nconst_int_calc,nconst_d&
&ble_calc,which_unity,basis_symmetry
character(*),intent(in) :: label
integer,pointer :: const_int(:),const_int_calc(:)
real(dp),pointer :: const_dble(:),const_dble_calc(:)
character(casl_fullkeysize),pointer :: const_int_name(:),const_dble_na&
&me(:)
integer iconst,xyzzyaaaa72,xyzzyaaab72,xyzzyaaac72,xyzzyaaad72,nchildr&
&en,xyzzyaaae72,ivector,xyzzyaaaf72,xyzzyaaag72,xyzzyaaah72,xyzzyaaai7&
&2,xyzzyaaaj72,xyzzyaaak72,xyzzyaaal72,xyzzyaaam72,xyzzyaaan72,xyzzyaa&
&ao72,xyzzyaaap72,xyzzyaaaq72,xyzzyaaar72,xyzzyaaas72,xyzzyaaat72,xyzz&
&yaaau72,xyzzyaaav72,xyzzyaaaw72,xyzzyaaax72,k
integer,allocatable :: xyzzyaaay72(:)
integer,pointer :: xyzzyaaaz72(:)=>null(),xyzzyaaba72(:)=>null(),xyzzy&
&aabb72(:)=>null()
real(dp) xyzzyaabc72,xyzzyaabd72,xyzzyaabe72
real(dp),pointer :: xyzzyaabf72(:,:)=>null(),xyzzyaabg72(:,:)=>null()
real(dp),parameter :: xyzzyaabh72=1.d-13,xyzzyaabi72=1.d-10
real(dp),allocatable :: xyzzyaabj72(:,:,:)
logical exists,is_block
nullify(const_int_name,const_int,const_int_calc,const_dble_name,const_&
&dble,const_dble_calc)
nconst_int=0
nconst_int_calc=0
nconst_dble=0
nconst_dble_calc=0
which_unity=0
basis_symmetry=1
select case(ibasis)
case(0)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for basis "none" must be 1.')
which_unity=1
case(xyzzyaaaa1)
which_unity=1
case(xyzzyaaab1)
call query_casl_item(trim(label),exists=exists,is_block=is_block,nchil&
&dren=nchildren)
if(.not.exists)then
xyzzyaaam72=iorder
xyzzyaabc72=-1.d0
call xyzzyaaja1(xyzzyaaam72,xyzzyaabc72,const_int_calc,const_int,const&
&_int_name)
nconst_int_calc=iorder
nconst_int=periodicity*sum(const_int_calc(1:iorder))
else
if(nchildren<iorder)call errstop_master('SETUP_BASIS_CONSTANTS','Numbe&
&r of stars of G-vectors listed for the cosine basis must be enough to&
& cover the requested expansion order; no default constants are genera&
&ted for this basis when any portion of them is provided.')
nconst_int_calc=iorder
allocate(const_int_calc(iorder),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:cosine')
const_int_calc=0
do xyzzyaaae72=1,iorder
call query_casl_item(trim(label)//':Star '//trim(i2s(xyzzyaaae72)),exi&
&sts=exists,is_block=is_block,nchildren=nchildren)
if(.not.exists)call errstop_master('SETUP_BASIS_CONSTANTS','Expected t&
&o find CASL item "'//trim(label)//':Star '//trim(i2s(xyzzyaaae72))//'&
&", but it''s not there.')
if(.not.is_block.or.nchildren==0)call errstop_master('SETUP_BASIS_CONS&
&TANTS','CASL item "'//trim(label)//':Star '//trim(i2s(xyzzyaaae72))//&
&'" expected to be a non-empty block, but no.')
const_int_calc(xyzzyaaae72)=nchildren
enddo
nconst_int=periodicity*sum(const_int_calc(1:iorder))
allocate(const_int(nconst_int),const_int_name(nconst_int),stat=xyzzyaa&
&ab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:cosine')
const_int=0
iconst=0
do xyzzyaaae72=1,iorder
do ivector=1,const_int_calc(xyzzyaaae72)
do xyzzyaaaf72=1,periodicity
iconst=iconst+1
const_int_name(iconst)='Star '//trim(i2s(xyzzyaaae72))//':G_'//trim(i2&
&s(ivector))//':%u'//trim(i2s(xyzzyaaaf72))
enddo
enddo
enddo
do iconst=1,nconst_int
call get_casl_item(trim(label)//':'//trim(const_int_name(iconst)),xyzz&
&yaaad72,xyzzyaaac72)
if(xyzzyaaac72/=0)call errstop_master('SETUP_BASIS_CONSTANTS','Problem&
& reading CASL item "'//trim(label)//':'//trim(const_int_name(iconst))&
&//'".')
const_int(iconst)=xyzzyaaad72
enddo
endif
nconst_dble_calc=(periodicity+1)*sum(const_int_calc(1:iorder))
allocate(const_dble_calc(nconst_dble_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:cosine')
const_dble_calc=0.d0
iconst=0
xyzzyaaaa72=0
do xyzzyaaae72=1,iorder
do ivector=1,const_int_calc(xyzzyaaae72)
xyzzyaabc72=0.d0
do xyzzyaaaf72=1,periodicity
iconst=iconst+1
do xyzzyaaag72=1,periodicity
const_dble_calc(iconst)=const_dble_calc(iconst)+const_int(xyzzyaaaa72+&
&xyzzyaaag72)*bmat(xyzzyaaaf72,xyzzyaaag72)
enddo
xyzzyaabc72=xyzzyaabc72+const_dble_calc(iconst)*const_dble_calc(iconst&
&)
enddo
xyzzyaaaa72=xyzzyaaaa72+periodicity
iconst=iconst+1
const_dble_calc(iconst)=xyzzyaabc72
enddo
enddo
call xyzzyaajb1(iorder,const_int_calc,sum(const_int_calc(1:iorder)),co&
&nst_dble_calc,which_unity)
case(xyzzyaaac1)
nconst_dble=3
nconst_int=1
allocate(const_dble_name(nconst_dble),const_dble(nconst_dble),const_in&
&t_name(nconst_int),const_int(nconst_int),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:cosine wit&
&h k-cutoff')
const_dble_name(1)='k_cut'
const_dble(1)=1.d0
const_dble_name(2)='p_0'
const_dble_name(3)='delta_p'
select case(periodicity)
case(3)
const_dble(2)=-2.d0
const_dble(3)=1.d0
case(2)
const_dble(2)=-1.5d0
const_dble(3)=0.5d0
case default
const_dble(2)=0.d0
const_dble(3)=1.d0
end select
const_int_name(1)='C'
const_int(1)=1
do iconst=1,nconst_int
call get_casl_item(trim(label)//':'//trim(const_int_name(iconst)),xyzz&
&yaaad72,xyzzyaaac72)
if(xyzzyaaac72==0)const_int(iconst)=xyzzyaaad72
enddo
do iconst=1,nconst_dble
call get_casl_item(trim(label)//':'//trim(const_dble_name(iconst)),xyz&
&zyaabc72,xyzzyaaac72)
if(xyzzyaaac72==0)const_dble(iconst)=xyzzyaabc72
enddo
xyzzyaaam72=-1
call xyzzyaaja1(xyzzyaaam72,const_dble(1),xyzzyaaaz72,xyzzyaaba72,excl&
&ude_zero=.true.)
nconst_int_calc=1+xyzzyaaam72
nconst_dble_calc=iorder*xyzzyaaam72+(periodicity+1)*sum(xyzzyaaaz72)
allocate(const_int_calc(nconst_int_calc),const_dble_calc(nconst_dble_c&
&alc),stat=xyzzyaaab72)
const_int_calc=0
const_dble_calc=0.d0
const_int_calc(1)=xyzzyaaam72
iconst=1
do xyzzyaaae72=1,xyzzyaaam72
iconst=iconst+1
const_int_calc(iconst)=xyzzyaaaz72(xyzzyaaae72)
enddo
select case(dimensionality)
case(3)
xyzzyaabd72=1.d0/volume
case(2)
xyzzyaabd72=1.d0/area
case(1)
xyzzyaabd72=1.d0/amat(1,1)
end select
iconst=iorder*xyzzyaaam72
xyzzyaaaa72=0
do xyzzyaaae72=1,xyzzyaaam72
do ivector=1,xyzzyaaaz72(xyzzyaaae72)
xyzzyaabc72=0.d0
do xyzzyaaaf72=1,periodicity
iconst=iconst+1
do xyzzyaaag72=1,periodicity
const_dble_calc(iconst)=const_dble_calc(iconst)+xyzzyaaba72(xyzzyaaaa7&
&2+xyzzyaaag72)*bmat(xyzzyaaaf72,xyzzyaaag72)
enddo
xyzzyaabc72=xyzzyaabc72+const_dble_calc(iconst)*const_dble_calc(iconst&
&)
enddo
xyzzyaaaa72=xyzzyaaaa72+periodicity
iconst=iconst+1
const_dble_calc(iconst)=xyzzyaabc72
enddo
enddo
if(iconst>nconst_dble_calc)call errstop_master('SETUP_BASIS_CONSTANTS'&
&,'Exceeded size of CONST_DBLE_CALC vector. Bug.')
iconst=0
xyzzyaaaa72=iorder*xyzzyaaam72
do xyzzyaaae72=1,xyzzyaaam72
xyzzyaabc72=sqrt(const_dble_calc(xyzzyaaaa72+periodicity+1))
do xyzzyaaan72=1,iorder
iconst=iconst+1
if(xyzzyaabc72<const_dble(1))const_dble_calc(iconst)=xyzzyaabd72*xyzzy&
&aabc72**(const_dble(2)+const_dble(3)*dble(xyzzyaaan72-1))*(xyzzyaabc7&
&2-const_dble(1))**const_int(1)
enddo
xyzzyaaaa72=xyzzyaaaa72+periodicity+1
enddo
if(iconst>iorder*xyzzyaaam72)call errstop_master('SETUP_BASIS_CONSTANT&
&S','Exceeded size of first part of CONST_DBLE_CALC vector. Bug.')
deallocate(xyzzyaaaz72,xyzzyaaba72)
case(xyzzyaaad1)
which_unity=1
nconst_dble=1
allocate(const_dble_name(nconst_dble),const_dble(nconst_dble),stat=xyz&
&zyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:r/(r+a) po&
&wer')
const_dble_name(1)='L'
const_dble(1)=1.d100
do iconst=1,nconst_dble
call get_casl_item(trim(label)//':'//trim(const_dble_name(iconst)),xyz&
&zyaabc72,xyzzyaaac72)
if(xyzzyaaac72==0)const_dble(iconst)=xyzzyaabc72
enddo
nconst_dble_calc=1
allocate(const_dble_calc(nconst_dble_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:r/(r+a) po&
&wer')
const_dble_calc(1)=const_dble(1)
case(xyzzyaaae1)
which_unity=1
case(xyzzyaaaf1)
which_unity=1
case(xyzzyaaag1)
if(mod(iorder,dimensionality)/=0)call errstop_master('SETUP_BASIS_CONS&
&TANTS','Expansion order for basis:natural power vectorial should be a&
& multiple of the dimensionality.')
basis_symmetry=-1
case(xyzzyaaah1)
which_unity=0
nconst_int=1+iorder
allocate(const_int_name(nconst_int),const_int(nconst_int),stat=xyzzyaa&
&ab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:natural po&
&lynomial')
const_int_name(1)='k0'
const_int(1)=0
do xyzzyaaaw72=1,iorder
const_int_name(xyzzyaaaw72+1)='Split:%u'//trim(i2s(xyzzyaaaw72))
const_int(xyzzyaaaw72+1)=1
enddo
do iconst=1,nconst_int
call get_casl_item(trim(label)//':'//trim(const_int_name(iconst)),xyzz&
&yaaad72,xyzzyaaac72)
if(xyzzyaaac72==0)const_int(iconst)=xyzzyaaad72
enddo
if(const_int(1)<0)call errstop_master('SETUP_BASIS_CONSTANTS','"k0" co&
&nstant for natural polynomial must be >= 0.')
if(any(const_int(2:)<1))call errstop_master('SETUP_BASIS_CONSTANTS','C&
&omponents of "Split" constant for natural polynomial must be >= 1.')
nconst_int_calc=nconst_int
allocate(const_int_calc(nconst_int_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:natural po&
&lynomial')
const_int_calc(1:nconst_int)=const_int(1:nconst_int)
case(xyzzyaaai1)
if(mod(iorder,dimensionality)/=0)call errstop_master('SETUP_BASIS_CONS&
&TANTS','Expansion order for basis:natural polynomial vectorial should&
& be a multiple of the dimensionality.')
basis_symmetry=-1
which_unity=0
xyzzyaaav72=iorder/dimensionality
nconst_int=1+xyzzyaaav72
allocate(const_int_name(nconst_int),const_int(nconst_int),stat=xyzzyaa&
&ab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:natural po&
&lynomial vectorial')
const_int_name(1)='k0'
const_int(1)=0
do xyzzyaaaw72=1,xyzzyaaav72
const_int_name(xyzzyaaaw72+1)='Split:%u'//trim(i2s(xyzzyaaaw72))
const_int(xyzzyaaaw72+1)=1
enddo
do iconst=1,nconst_int
call get_casl_item(trim(label)//':'//trim(const_int_name(iconst)),xyzz&
&yaaad72,xyzzyaaac72)
if(xyzzyaaac72==0)const_int(iconst)=xyzzyaaad72
enddo
if(const_int(1)<0)call errstop_master('SETUP_BASIS_CONSTANTS','"k0" co&
&nstant for natural polynomial vectorial must be >= 0.')
if(any(const_int(2:)<1))call errstop_master('SETUP_BASIS_CONSTANTS','C&
&omponents of "Split" constant for natural polynomial vectorial must b&
&e >= 1.')
nconst_int_calc=nconst_int
allocate(const_int_calc(nconst_int_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:natural po&
&lynomial vectorial')
const_int_calc(1:nconst_int)=const_int(1:nconst_int)
case(xyzzyaaaj1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for basis:rpa must be 1.')
if(dimensionality==1)call errstop_master('SETUP_BASIS_CONSTANTS','basi&
&s:rpa not supported in 1D systems.')
case(xyzzyaaak1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for basis:logarithmic cusp must be 1.')
case(xyzzyaaal1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for basis:dipole cusp must be 1.')
case(xyzzyaaam1)
nconst_dble=iorder
allocate(const_dble_name(nconst_dble),const_dble(nconst_dble),stat=xyz&
&zyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:hipow')
do iconst=1,nconst_dble
const_dble_name(iconst)='k'//trim(i2s(iconst))
const_dble(iconst)=0.5d0*dble(-2+iconst)
enddo
do iconst=1,nconst_dble
call get_casl_item(trim(label)//':'//trim(const_dble_name(iconst)),xyz&
&zyaabc72,xyzzyaaac72)
if(xyzzyaaac72==0)const_dble(iconst)=xyzzyaabc72
enddo
nconst_int_calc=iorder
allocate(const_int_calc(nconst_int_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:hipow')
const_int_calc(:)=nint(2.d0*const_dble(:))
case(xyzzyaaan1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for basis:tilted dipole cusp must be 1.')
case(xyzzyaaao1)
call query_casl_item(trim(label),exists=exists,is_block=is_block,nchil&
&dren=nchildren)
if(.not.exists)then
nullify(xyzzyaabf72,xyzzyaabg72)
call xyzzyaajc1(xyzzyaaah72,xyzzyaabf72,xyzzyaabg72)
nconst_dble=6*xyzzyaaah72
allocate(const_dble_name(nconst_dble),const_dble(nconst_dble),stat=xyz&
&zyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:nu')
iconst=0
do ivector=1,xyzzyaaah72
do xyzzyaaaw72=1,3
iconst=iconst+1
const_dble_name(iconst)='a_'//trim(i2s(ivector))//':%u'//trim(i2s(xyzz&
&yaaaw72))
const_dble(iconst)=xyzzyaabf72(xyzzyaaaw72,ivector)
enddo
enddo
do ivector=1,xyzzyaaah72
do xyzzyaaaw72=1,3
iconst=iconst+1
const_dble_name(iconst)='b_'//trim(i2s(ivector))//':%u'//trim(i2s(xyzz&
&yaaaw72))
const_dble(iconst)=xyzzyaabg72(xyzzyaaaw72,ivector)
enddo
enddo
deallocate(xyzzyaabf72,xyzzyaabg72)
else
if(mod(nchildren,2)/=0)call errstop_master('SETUP_BASIS_CONSTANTS','Nu&
&mber of vectors in basis:nu constants block must be an even number.')
xyzzyaaah72=nchildren/2
nconst_dble=6*xyzzyaaah72
allocate(const_dble_name(nconst_dble),const_dble(nconst_dble),stat=xyz&
&zyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:nu')
iconst=0
do ivector=1,xyzzyaaah72
do xyzzyaaaw72=1,3
iconst=iconst+1
const_dble_name(iconst)='a_'//trim(i2s(ivector))//':%u'//trim(i2s(xyzz&
&yaaaw72))
enddo
enddo
do ivector=1,xyzzyaaah72
do xyzzyaaaw72=1,3
iconst=iconst+1
const_dble_name(iconst)='b_'//trim(i2s(ivector))//':%u'//trim(i2s(xyzz&
&yaaaw72))
enddo
enddo
do iconst=1,nconst_dble
call get_casl_item(trim(label)//':'//trim(const_dble_name(iconst)),con&
&st_dble(iconst),xyzzyaaac72)
if(xyzzyaaac72/=0)call errstop_master('SETUP_BASIS_CONSTANTS','Constan&
&ts in basis:nu must be either entirely omitted or specified in full a&
&s three-dimensional vectors regardless of system dimensionality.')
enddo
endif
nconst_dble_calc=1+6*xyzzyaaah72+2*xyzzyaaah72*xyzzyaaah72
nconst_int_calc=5
allocate(const_dble_calc(nconst_dble_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:nu')
allocate(const_int_calc(nconst_int_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','basis:nu')
const_int_calc(1)=xyzzyaaah72
if(abs(ddot(3,amat(1,1),3,amat(2,1),3))<=xyzzyaabi72.and.abs(ddot(3,am&
&at(1,1),3,amat(3,1),3))<=xyzzyaabi72.and.abs(ddot(3,amat(2,1),3,amat(&
&3,1),3))<=xyzzyaabi72)then
const_int_calc(2)=1
else
const_int_calc(2)=0
endif
const_int_calc(3)=3*xyzzyaaah72+2
const_int_calc(4)=6*xyzzyaaah72+2
const_int_calc(5)=6*xyzzyaaah72+xyzzyaaah72*xyzzyaaah72+2
const_dble_calc(1)=1.d0/wigner_seitz_radius
const_dble_calc(2:nconst_dble+1)=const_dble(1:nconst_dble)
k=1+nconst_dble
do xyzzyaaaw72=1,xyzzyaaah72
do xyzzyaaax72=1,xyzzyaaah72
k=k+1
const_dble_calc(k)=dot_product(const_dble(1+3*(xyzzyaaaw72-1):3*xyzzya&
&aaw72),const_dble(1+3*(xyzzyaaax72-1):3*xyzzyaaax72))
enddo
enddo
do xyzzyaaaw72=1,xyzzyaaah72
do xyzzyaaax72=1,xyzzyaaah72
k=k+1
const_dble_calc(k)=dot_product(const_dble(1+3*(xyzzyaaah72+xyzzyaaaw72&
&-1):3*(xyzzyaaah72+xyzzyaaaw72)),const_dble(1+3*(xyzzyaaah72+xyzzyaaa&
&x72-1):3*(xyzzyaaah72+xyzzyaaax72)))
enddo
enddo
case(xyzzyaaap1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for cut-off functions must be 1.')
nconst_int=1
allocate(const_int_name(nconst_int),const_int(nconst_int),stat=xyzzyaa&
&ab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:polynomia&
&l')
const_int_name(1)='C'
const_int(1)=3
do iconst=1,nconst_int
call get_casl_item(trim(label)//':'//trim(const_int_name(iconst)),xyzz&
&yaaad72,xyzzyaaac72)
if(xyzzyaaac72==0)const_int(iconst)=xyzzyaaad72
enddo
nconst_int_calc=1
allocate(const_int_calc(nconst_int_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:polynomia&
&l')
const_int_calc(1)=const_int(1)
case(xyzzyaaaq1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for cut-off functions must be 1.')
nconst_int=1
allocate(const_int_name(nconst_int),const_int(nconst_int),stat=xyzzyaa&
&ab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:alt polyn&
&omial')
const_int_name(1)='C'
const_int(1)=3
do iconst=1,nconst_int
call get_casl_item(trim(label)//':'//trim(const_int_name(iconst)),xyzz&
&yaaad72,xyzzyaaac72)
if(xyzzyaaac72==0)const_int(iconst)=xyzzyaaad72
enddo
nconst_int_calc=1
allocate(const_int_calc(nconst_int_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:alt polyn&
&omial')
const_int_calc(1)=const_int(1)
case(xyzzyaaar1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for cut-off functions must be 1.')
nconst_dble=1
allocate(const_dble_name(nconst_dble),const_dble(nconst_dble),stat=xyz&
&zyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:gaussian'&
&)
const_dble_name(1)='L_hard'
const_dble(1)=6.d0
do iconst=1,nconst_dble
call get_casl_item(trim(label)//':'//trim(const_dble_name(iconst)),xyz&
&zyaabc72,xyzzyaaac72)
if(xyzzyaaac72==0)const_dble(iconst)=xyzzyaabc72
enddo
nconst_dble_calc=1
allocate(const_dble_calc(nconst_dble_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:gaussian'&
&)
const_dble_calc(1)=const_dble(1)
case(xyzzyaaas1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for cut-off functions must be 1.')
call query_casl_item(label,exists=exists,is_block=is_block,nchildren=n&
&children)
if(.not.exists.or..not.is_block.or.nchildren<1)call errstop_master('SE&
&TUP_BASIS_CONSTANTS','Need to specify xyz_i terms explicitly (cutoff:&
&anisotropic polynomial')
xyzzyaaah72=nchildren
call query_casl_item(trim(label)//':C',exists=exists)
if(exists)xyzzyaaah72=xyzzyaaah72-1
xyzzyaaar72=0
do
call query_casl_item(trim(label)//':Frame '//trim(i2s(xyzzyaaar72+1)),&
&exists=exists)
if(.not.exists)exit
xyzzyaaar72=xyzzyaaar72+1
xyzzyaaah72=xyzzyaaah72-1
call resize_pointer((/xyzzyaaar72/),xyzzyaabb72)
call query_casl_item(trim(label)//':Frame '//trim(i2s(xyzzyaaar72))//'&
&:Atoms',exists=exists,nchildren=nchildren)
xyzzyaabb72(xyzzyaaar72)=nchildren
enddo
if(xyzzyaaar72==1)then
if(xyzzyaabb72(1)==0)xyzzyaaar72=-1
endif
if(xyzzyaaah72<1)call errstop_master('SETUP_BASIS_CONSTANTS','You need&
& to provide '//trim(i2s(xyzzyaaah72))//' '//trim(i2s(dimensionality))&
&//'-dimensional exponent vectors. None were found (cutoff:anisotropic&
& polynomial).')
nconst_int=1+dimensionality*xyzzyaaah72
nconst_dble=xyzzyaaah72
if(xyzzyaaar72/=0)then
nconst_int=nconst_int+sum(xyzzyaabb72)
nconst_dble=nconst_dble+abs(xyzzyaaar72)*dimensionality*dimensionality
endif
allocate(const_int_name(nconst_int),const_int(nconst_int),const_dble_n&
&ame(nconst_dble),const_dble(nconst_dble),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:anisotrop&
&ic polynomial')
const_int=-1
const_int_name(1)='C'
const_int(1)=3
iconst=1
do ivector=1,xyzzyaaah72
do xyzzyaaai72=1,dimensionality
iconst=iconst+1
const_int_name(iconst)='xyz_'//trim(i2s(ivector))//':p:%u'//trim(i2s(x&
&yzzyaaai72))
enddo
enddo
do xyzzyaaas72=1,abs(xyzzyaaar72)
do xyzzyaaat72=1,xyzzyaabb72(xyzzyaaas72)
iconst=iconst+1
const_int_name(iconst)='Frame '//trim(i2s(xyzzyaaas72))//':Atoms:%u'//&
&trim(i2s(xyzzyaaat72))
enddo
enddo
const_dble=0.d0
iconst=0
do ivector=1,xyzzyaaah72
iconst=iconst+1
const_dble_name(iconst)='xyz_'//trim(i2s(ivector))//':c'
enddo
do xyzzyaaas72=1,abs(xyzzyaaar72)
do xyzzyaaai72=1,dimensionality
do xyzzyaaaj72=1,dimensionality
iconst=iconst+1
const_dble_name(iconst)='Frame '//trim(i2s(xyzzyaaas72))//':u_'//trim(&
&i2s(xyzzyaaai72))//':%u'//trim(i2s(xyzzyaaaj72))
enddo
enddo
enddo
do iconst=1,nconst_int
call get_casl_item(trim(label)//':'//trim(const_int_name(iconst)),xyzz&
&yaaad72,xyzzyaaac72)
if(xyzzyaaac72==0)const_int(iconst)=xyzzyaaad72
enddo
do iconst=1,nconst_dble
call get_casl_item(trim(label)//':'//trim(const_dble_name(iconst)),xyz&
&zyaabc72,xyzzyaaac72)
if(xyzzyaaac72==0)const_dble(iconst)=xyzzyaabc72
enddo
iconst=1
xyzzyaaaq72=0
xyzzyaaap72=0
do ivector=1,xyzzyaaah72
if(any(const_int(iconst+1:iconst+dimensionality)<0))call errstop_maste&
&r('SETUP_BASIS_CONSTANTS','Exponent vector in xyz_'//trim(i2s(ivector&
&))//' not (fully) specified or has negative components. You need to p&
&rovide '//trim(i2s(xyzzyaaah72))//' '//trim(i2s(dimensionality))//'-d&
&imensional exponent vectors (cutoff:anisotropic polynomial).')
if(all(const_int(iconst+1:iconst+dimensionality)==0))call errstop_mast&
&er('SETUP_BASIS_CONSTANTS','Exponent vector in xyz_'//trim(i2s(ivecto&
&r))//' is zero. Use the isotropic cutoff function for that (cutoff:an&
&isotropic polynomial).')
xyzzyaaao72=sum(const_int(iconst+1:iconst+dimensionality))
if(xyzzyaaap72==0)xyzzyaaap72=xyzzyaaao72
if(xyzzyaaap72/=xyzzyaaao72)call errstop_master('SETUP_BASIS_CONSTANTS&
&','Sum of exponents in exponent vector of xyz_'//trim(i2s(ivector))//&
&' differs from that in previous exponent vectors.  The sum of exponen&
&ts must be the same for all exponent vectors in a term.')
xyzzyaaaa72=iconst+dimensionality
do xyzzyaaal72=ivector+1,xyzzyaaah72
if(all(const_int(iconst+1:iconst+dimensionality)==const_int(xyzzyaaaa7&
&2+1:xyzzyaaaa72+dimensionality)))call errstop_master('SETUP_BASIS_CON&
&STANTS','Exponent vectors in xyz_'//trim(i2s(ivector))//' and n xyz_'&
&//trim(i2s(xyzzyaaal72))//' are equal.  All vectors should be differe&
&nt (cutoff:anisotropic polynomial).')
xyzzyaaaa72=xyzzyaaaa72+dimensionality
enddo
iconst=iconst+dimensionality
xyzzyaaaq72=xyzzyaaaq72+1
if(const_dble(xyzzyaaaq72)==0.d0)call errstop_master('SETUP_BASIS_CONS&
&TANTS','Coefficient in xyz_'//trim(i2s(ivector))//' is zero.  This is&
& not allowed.')
enddo
if(mod(xyzzyaaap72,2)==0)then
basis_symmetry=1
else
basis_symmetry=-1
endif
if(xyzzyaaar72/=0)then
if(xyzzyaaar72>0)then
if(nitot==0)call errstop_master('SETUP_BASIS_CONSTANTS','Cannot use mo&
&dified reference frames in systems without atoms (cutoff:anisotropic &
&polynomial).')
allocate(xyzzyaaay72(nitot),stat=xyzzyaaab72)
endif
allocate(xyzzyaabj72(dimensionality,dimensionality,abs(xyzzyaaar72)),s&
&tat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','tmats')
xyzzyaabj72=0.d0
xyzzyaaay72=0
iconst=xyzzyaaah72
xyzzyaaaa72=1+xyzzyaaah72*dimensionality
do xyzzyaaas72=1,abs(xyzzyaaar72)
do xyzzyaaat72=1,xyzzyaabb72(xyzzyaaas72)
xyzzyaaaa72=xyzzyaaaa72+1
xyzzyaaau72=const_int(xyzzyaaaa72)
if(xyzzyaaay72(xyzzyaaau72)/=0)call errstop_master('SETUP_BASIS_CONSTA&
&NTS','Atom '//trim(i2s(xyzzyaaau72))//' has been assigned reference f&
&rames twice (cutoff:anisotropic polynomial).')
xyzzyaaay72(xyzzyaaau72)=xyzzyaaas72
enddo
do xyzzyaaai72=1,dimensionality
do xyzzyaaaj72=1,dimensionality
iconst=iconst+1
xyzzyaabj72(xyzzyaaaj72,xyzzyaaai72,xyzzyaaas72)=const_dble(iconst)
enddo
enddo
do xyzzyaaai72=1,dimensionality
xyzzyaabc72=sqrt(sum(xyzzyaabj72(:,xyzzyaaai72,xyzzyaaas72)**2))
if(xyzzyaabc72==0.d0)call errstop_master('SETUP_BASIS_CONSTANTS','Refe&
&rence frame #'//trim(i2s(xyzzyaaas72))//' ill-defined because vector &
&u'//trim(i2s(xyzzyaaai72))//' has zero norm (cutoff:anisotropic polyn&
&omial).')
xyzzyaabj72(:,xyzzyaaai72,xyzzyaaas72)=xyzzyaabj72(:,xyzzyaaai72,xyzzy&
&aaas72)/xyzzyaabc72
enddo
select case(dimensionality)
case(1)
xyzzyaabe72=1.d0
case(2)
xyzzyaabe72=xyzzyaabj72(1,1,xyzzyaaas72)*xyzzyaabj72(2,2,xyzzyaaas72)-&
&xyzzyaabj72(1,2,xyzzyaaas72)*xyzzyaabj72(2,1,xyzzyaaas72)
case(3)
xyzzyaabe72=xyzzyaabj72(1,1,xyzzyaaas72)*(xyzzyaabj72(2,2,xyzzyaaas72)&
&*xyzzyaabj72(3,3,xyzzyaaas72)-xyzzyaabj72(2,3,xyzzyaaas72)*xyzzyaabj7&
&2(3,2,xyzzyaaas72))-xyzzyaabj72(1,2,xyzzyaaas72)*(xyzzyaabj72(2,1,xyz&
&zyaaas72)*xyzzyaabj72(3,3,xyzzyaaas72)-xyzzyaabj72(2,3,xyzzyaaas72)*x&
&yzzyaabj72(3,1,xyzzyaaas72))+xyzzyaabj72(1,3,xyzzyaaas72)*(xyzzyaabj7&
&2(2,1,xyzzyaaas72)*xyzzyaabj72(3,2,xyzzyaaas72)-xyzzyaabj72(2,2,xyzzy&
&aaas72)*xyzzyaabj72(3,1,xyzzyaaas72))
end select
if(abs(xyzzyaabe72)<xyzzyaabh72)call errstop_master('SETUP_BASIS_CONST&
&ANTS','Reference frame #'//trim(i2s(xyzzyaaas72))//' ill-defined beca&
&use vectors are collinear (cutoff:anisotropic polynomial).')
if(abs(abs(xyzzyaabe72)-1.d0)>xyzzyaabh72)call errstop_master('SETUP_B&
&ASIS_CONSTANTS','Reference frame #'//trim(i2s(xyzzyaaas72))//' is not&
& orthogonal.  Reference frames are required to be orthogonal (cutoff:&
&anisotropic polynomial)')
if(xyzzyaabe72<0.d0)call errstop_master('SETUP_BASIS_CONSTANTS','Refer&
&ence frame #'//trim(i2s(xyzzyaaas72))//' is left-handed.  Reference f&
&rames are required to be right-handed (cutoff:anisotropic polynomial)&
&.')
enddo
do xyzzyaaat72=1,nitot
if(xyzzyaaay72(xyzzyaaat72)==0)call errstop_master('SETUP_BASIS_CONSTA&
&NTS','Atom '//trim(i2s(xyzzyaaat72))//' has been asigned no reference&
& frame (cutoff:anisotropic polynomial).')
enddo
endif
nconst_int_calc=4+xyzzyaaah72*dimensionality
nconst_dble_calc=dimensionality*dimensionality+xyzzyaaah72
if(xyzzyaaar72>0)then
nconst_int_calc=nconst_int_calc+nitot
nconst_dble_calc=nconst_dble_calc+2*dimensionality*dimensionality*xyzz&
&yaaar72
elseif(xyzzyaaar72<0)then
nconst_dble_calc=nconst_dble_calc+2*dimensionality*dimensionality
endif
allocate(const_int_calc(nconst_int_calc),const_dble_calc(nconst_dble_c&
&alc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:anisotrop&
&ic polynomial')
const_int_calc(1)=const_int(1)
const_int_calc(2)=xyzzyaaap72
const_int_calc(3)=xyzzyaaah72
const_int_calc(4)=xyzzyaaar72
iconst=4
xyzzyaaaa72=1
do ivector=1,xyzzyaaah72
do xyzzyaaai72=1,dimensionality
iconst=iconst+1
xyzzyaaaa72=xyzzyaaaa72+1
const_int_calc(iconst)=const_int(xyzzyaaaa72)
enddo
enddo
if(xyzzyaaar72>0)then
do xyzzyaaat72=1,nitot
iconst=iconst+1
const_int_calc(iconst)=xyzzyaaay72(xyzzyaaat72)
enddo
endif
iconst=0
do ivector=1,xyzzyaaah72
iconst=iconst+1
const_dble_calc(iconst)=const_dble(iconst)
enddo
if(xyzzyaaar72/=0)then
do xyzzyaaas72=1,abs(xyzzyaaar72)
do xyzzyaaai72=1,dimensionality
do xyzzyaaaj72=1,dimensionality
iconst=iconst+1
const_dble_calc(iconst)=xyzzyaabj72(xyzzyaaaj72,xyzzyaaai72,xyzzyaaas7&
&2)
enddo
enddo
do xyzzyaaai72=1,dimensionality
do xyzzyaaaj72=1,dimensionality
xyzzyaabc72=0.d0
do xyzzyaaak72=1,dimensionality
xyzzyaabc72=xyzzyaabc72+xyzzyaabj72(xyzzyaaak72,xyzzyaaai72,xyzzyaaas7&
&2)*xyzzyaabj72(xyzzyaaak72,xyzzyaaaj72,xyzzyaaas72)
enddo
iconst=iconst+1
const_dble_calc(iconst)=xyzzyaabc72
enddo
enddo
enddo
endif
if(associated(xyzzyaabb72))deallocate(xyzzyaabb72)
if(allocated(xyzzyaabj72))deallocate(xyzzyaabj72)
if(allocated(xyzzyaaay72))deallocate(xyzzyaaay72)
case(xyzzyaaat1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for cut-off functions must be 1.')
if(dimensionality==3.or.heg_nlayers==1)call errstop_master('SETUP_BASI&
&S_CONSTANTS','Cut-off function "quasicusp" only available on 1D bi-/m&
&ulti-wires and 2D bi-/multi-layers.')
case(xyzzyaaau1)
if(iorder/=1)call errstop_master('SETUP_BASIS_CONSTANTS','Expansion or&
&der for cut-off functions must be 1.')
nconst_int=1
allocate(const_int_name(nconst_int),const_int(nconst_int),stat=xyzzyaa&
&ab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:spline')
const_int_name(1)='C'
const_int(1)=2
do iconst=1,nconst_int
call get_casl_item(trim(label)//':'//trim(const_int_name(iconst)),xyzz&
&yaaad72,xyzzyaaac72)
if(xyzzyaaac72==0)const_int(iconst)=xyzzyaaad72
enddo
if(const_int(1)<0.or.const_int(1)>3)call errstop_master('SETUP_BASIS_C&
&ONSTANTS','Truncation order C for cutoff:spline must be 0--3.')
nconst_int_calc=1
allocate(const_int_calc(nconst_int_calc),stat=xyzzyaaab72)
call check_alloc(xyzzyaaab72,'SETUP_BASIS_CONSTANTS','cutoff:spline')
const_int_calc(1)=const_int(1)
end select
end subroutine xyzzyaafa1
subroutine xyzzyaafb1(ibasis,iorder,nconst_int_calc,const_int_calc,nco&
&nst_dble_calc,const_dble_calc,nparam,nparam_calc,param_name,param,has&
&_lolim,lolim,has_hilim,hilim,is_shallow)
implicit none
integer,intent(in) :: ibasis,iorder,nconst_int_calc,nconst_dble_calc
integer,intent(out) :: nparam,nparam_calc
integer,pointer :: const_int_calc(:)
real(dp),pointer :: const_dble_calc(:),param(:),lolim(:),hilim(:)
logical,pointer :: has_lolim(:),has_hilim(:),is_shallow(:)
character(casl_keysize),pointer :: param_name(:)
integer xyzzyaaaa73,xyzzyaaab73,xyzzyaaac73,xyzzyaaad73,xyzzyaaae73,xy&
&zzyaaaf73
real(dp) xyzzyaaag73
nullify(param_name,param)
select case(ibasis)
case(0)
nparam=0
nparam_calc=0
case(xyzzyaaaa1)
nparam=0
nparam_calc=0
case(xyzzyaaab1)
nparam=0
nparam_calc=0
case(xyzzyaaac1)
nparam=0
nparam_calc=0
case(xyzzyaaad1)
nparam=1
nparam_calc=1
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_PARAMETERS','basis:r/(r+a) p&
&ower')
param(1)=3.d0
has_lolim(1)=.true.
lolim(1)=1.1d-8
has_hilim(1)=.false.
hilim(1)=0.d0
is_shallow(1)=.true.
param_name(1)='a'
case(xyzzyaaae1)
nparam=2
nparam_calc=2
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_PARAMETERS','basis:r/(r^b+a)&
& power')
param(1:2)=(/3.d0,1.3d0/)
has_lolim(1:2)=.true.
lolim(1:2)=(/1.1d-8,1.d0/)
has_hilim(1:2)=.false.
hilim(1:2)=0.d0
is_shallow(1:2)=.true.
param_name(1:2)=(/'a','b'/)
case(xyzzyaaaf1)
nparam=1
nparam_calc=1
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_PARAMETERS','basis:1/(r+a) p&
&ower')
param(1)=3.d0
has_lolim(1)=.true.
lolim(1)=1.1d-8
has_hilim(1)=.false.
hilim(1)=0.d0
is_shallow(1)=.true.
param_name(1)='a'
case(xyzzyaaag1)
nparam=0
nparam_calc=0
case(xyzzyaaah1)
nparam=sum(const_int_calc(2:iorder+1))-iorder
nparam_calc=nparam
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_PARAMETERS','basis:natural p&
&ower expansion')
xyzzyaaae73=0
xyzzyaaaf73=const_int_calc(1)-1
do xyzzyaaac73=1,iorder
xyzzyaaaf73=xyzzyaaaf73+1
do xyzzyaaad73=2,const_int_calc(1+xyzzyaaac73)
xyzzyaaae73=xyzzyaaae73+1
xyzzyaaaf73=xyzzyaaaf73+1
param(xyzzyaaae73)=0.d0
has_lolim(xyzzyaaae73)=.false.
lolim(xyzzyaaae73)=0.d0
has_hilim(xyzzyaaae73)=.false.
hilim(xyzzyaaae73)=0.d0
is_shallow(xyzzyaaae73)=.false.
param_name(xyzzyaaae73)='c_'//trim(i2s(xyzzyaaaf73))
enddo
enddo
case(xyzzyaaai1)
xyzzyaaab73=iorder/dimensionality
nparam=sum(const_int_calc(2:xyzzyaaab73+1))-xyzzyaaab73
nparam_calc=nparam
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_PARAMETERS','basis:natural p&
&ower expansion vectorial')
xyzzyaaae73=0
xyzzyaaaf73=const_int_calc(1)-1
do xyzzyaaac73=1,xyzzyaaab73
xyzzyaaaf73=xyzzyaaaf73+1
do xyzzyaaad73=2,const_int_calc(1+xyzzyaaac73)
xyzzyaaae73=xyzzyaaae73+1
xyzzyaaaf73=xyzzyaaaf73+1
param(xyzzyaaae73)=0.d0
has_lolim(xyzzyaaae73)=.false.
lolim(xyzzyaaae73)=0.d0
has_hilim(xyzzyaaae73)=.false.
hilim(xyzzyaaae73)=0.d0
is_shallow(xyzzyaaae73)=.false.
param_name(xyzzyaaae73)='c_'//trim(i2s(xyzzyaaaf73))
enddo
enddo
case(xyzzyaaaj1)
nparam=1
nparam_calc=3
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_CONSTANTS','basis:rpa')
param_name(1)='F'
has_lolim(1)=.true.
lolim(1)=1.d-7
has_hilim(1)=.false.
hilim(1)=0.d0
param(1)=1.d0
case(xyzzyaaak1)
nparam=0
nparam_calc=0
case(xyzzyaaal1)
nparam=0
nparam_calc=0
case(xyzzyaaam1)
nparam=0
nparam_calc=0
case(xyzzyaaan1)
nparam=1
nparam_calc=1
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_CONSTANTS','basis: tilted di&
&pole cusp')
param_name(1)='theta'
has_lolim(1)=.false.
lolim(1)=0.d0
has_hilim(1)=.true.
hilim(1)=0.4253407848436730058
param(1)=0.d0
is_shallow(1)=.true.
case(xyzzyaaao1)
nparam=0
nparam_calc=0
case(xyzzyaaap1)
nparam=1
nparam_calc=4
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_CONSTANTS','cutoff:polynomia&
&l')
param_name(1)='L'
has_lolim(1)=.true.
lolim(1)=0.5d0
has_hilim(1)=.false.
hilim(1)=0.d0
param(1)=3.d0
if(isperiodic)then
lolim(1)=min(0.5d0,0.1d0*wigner_seitz_radius)
has_hilim(1)=.true.
hilim(1)=0.999999d0*wigner_seitz_radius
param(1)=0.999d0*hilim(1)
endif
is_shallow(1)=.true.
case(xyzzyaaaq1)
nparam=1
nparam_calc=1
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_CONSTANTS','cutoff:alt polyn&
&omial')
param_name(1)='L'
has_lolim(1)=.true.
lolim(1)=0.5d0
has_hilim(1)=.false.
hilim(1)=0.d0
param(1)=3.d0
if(isperiodic)then
lolim(1)=min(0.5d0,0.1d0*wigner_seitz_radius)
has_hilim(1)=.true.
hilim(1)=0.999999d0*wigner_seitz_radius
param(1)=0.999d0*hilim(1)
endif
is_shallow(1)=.true.
case(xyzzyaaar1)
nparam=1
nparam_calc=3
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_CONSTANTS','cutoff:gaussian'&
&)
param_name(1)='L'
has_lolim(1)=.true.
lolim(1)=0.5d0
has_hilim(1)=.false.
hilim(1)=0.d0
xyzzyaaag73=const_dble_calc(1)
param(1)=3.d0
if(isperiodic)then
lolim(1)=min(0.5d0,0.1d0*wigner_seitz_radius/xyzzyaaag73)
has_hilim(1)=.true.
hilim(1)=0.999999d0*wigner_seitz_radius/xyzzyaaag73
param(1)=0.999d0*hilim(1)
endif
is_shallow(1)=.true.
case(xyzzyaaas1)
nparam=1
nparam_calc=4
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_CONSTANTS','cutoff:anisotrop&
&ic polynomial')
param_name(1)='L'
has_lolim(1)=.true.
lolim(1)=0.5d0
has_hilim(1)=.false.
hilim(1)=0.d0
param(1)=3.d0
if(isperiodic)then
lolim(1)=min(0.5d0,0.1d0*wigner_seitz_radius)
has_hilim(1)=.true.
hilim(1)=0.999999d0*wigner_seitz_radius
param(1)=0.999d0*hilim(1)
endif
is_shallow(1)=.true.
case(xyzzyaaat1)
nparam=1
nparam_calc=3
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_CONSTANTS','cutoff:quasicusp&
&')
param_name(1)='L'
has_lolim(1)=.true.
lolim(1)=0.5d0
has_hilim(1)=.false.
hilim(1)=0.d0
param(1)=1.d0
if(isperiodic)then
lolim(1)=min(0.5d0,0.1d0*wigner_seitz_radius)
has_hilim(1)=.true.
hilim(1)=0.999999d0*wigner_seitz_radius
param(1)=min(1.d0,0.999d0*hilim(1))
endif
is_shallow(1)=.true.
case(xyzzyaaau1)
nparam=2
nparam_calc=4
allocate(param_name(nparam),param(nparam),has_lolim(nparam),lolim(npar&
&am),has_hilim(nparam),hilim(nparam),is_shallow(nparam),stat=xyzzyaaaa&
&73)
call check_alloc(xyzzyaaaa73,'SETUP_BASIS_CONSTANTS','cutoff:quasicusp&
&')
param_name(1:2)=(/'L','x'/)
has_lolim(1:2)=(/.true.,.true./)
lolim(1:2)=(/0.1d0,0.05d0/)
has_hilim(1:2)=(/.false.,.true./)
hilim(1:2)=(/0.d0,0.95d0/)
param(1:2)=(/1.d0,0.5d0/)
if(isperiodic)then
lolim(1)=min(0.1d0,0.1d0*wigner_seitz_radius)
has_hilim(1)=.true.
hilim(1)=0.999999d0*wigner_seitz_radius
param(1)=min(1.d0,0.999d0*hilim(1))
endif
is_shallow(1:2)=(/.true.,.true./)
end select
end subroutine xyzzyaafb1
subroutine xyzzyaafc1(ibasis,iorder,nparam,nparam_calc,param,nconst_in&
&t_calc,const_int_calc,nconst_dble_calc,const_dble_calc,param_calc)
implicit none
integer,intent(in) :: ibasis,iorder,nparam,nparam_calc,nconst_int_calc&
&,nconst_dble_calc,const_int_calc(nconst_int_calc)
real(dp),intent(in) :: param(nparam),const_dble_calc(nconst_dble_calc)
real(dp),intent(inout) :: param_calc(nparam_calc)
select case(ibasis)
case(0)
continue
case(xyzzyaaaa1)
continue
case(xyzzyaaab1)
continue
case(xyzzyaaac1)
continue
case(xyzzyaaad1)
param_calc(1)=param(1)
case(xyzzyaaae1)
param_calc(1:2)=param(1:2)
case(xyzzyaaaf1)
param_calc(1)=param(1)
case(xyzzyaaag1)
continue
case(xyzzyaaah1)
if(nparam>0)param_calc(1:nparam)=param(1:nparam)
case(xyzzyaaai1)
if(nparam>0)param_calc(1:nparam)=param(1:nparam)
case(xyzzyaaaj1)
param_calc(1)=param(1)
param_calc(2)=1.d0/param(1)
param_calc(3)=sqrt(1.d0/param(1))
case(xyzzyaaak1)
continue
case(xyzzyaaal1)
continue
case(xyzzyaaam1)
continue
case(xyzzyaaan1)
param_calc(1)=param(1)
case(xyzzyaaao1)
continue
case(xyzzyaaap1)
param_calc(1)=param(1)
param_calc(2)=1.d0/param(1)
param_calc(3)=-const_int_calc(1)/param(1)
param_calc(4)=(const_int_calc(1)-1)*const_int_calc(1)/(param(1)*param(&
&1))
case(xyzzyaaaq1)
param_calc(1)=param(1)
case(xyzzyaaar1)
param_calc(1)=1.d0/param(1)
param_calc(2)=1.d0/(param(1)*param(1))
param_calc(3)=param(1)*const_dble_calc(1)
case(xyzzyaaas1)
param_calc(1)=param(1)
param_calc(2)=1.d0/param(1)
param_calc(3)=-const_int_calc(1)/param(1)
param_calc(4)=(const_int_calc(1)-1)*const_int_calc(1)/(param(1)*param(&
&1))
case(xyzzyaaat1)
param_calc(1)=param(1)
param_calc(2)=1.d0/param(1)
param_calc(3)=param(1)*param(1)
case(xyzzyaaau1)
param_calc(1)=param(1)
param_calc(2)=param(2)*param(1)
param_calc(3)=1.d0/(param(1)*(1.d0-param(2)))
param_calc(4)=1.d0/(param(1)*(1.d0-param(2)))**2
end select
end subroutine xyzzyaafc1
subroutine gbasis_reset_aniso()
implicit none
integer xyzzyaaaa75
if(allocated(xyzzyaacc1))deallocate(xyzzyaacc1)
xyzzyaabm1=.false.
xyzzyaabn1=.false.
if(xyzzyaaaw1>0)then
allocate(xyzzyaacc1(xyzzyaaaw1),stat=xyzzyaaaa75)
call check_alloc(xyzzyaaaa75,'GBASIS_RESET_ANISO','deriv_aniso_eebasis&
&')
xyzzyaacc1=.not.xyzzyaaca1
xyzzyaabm1=any(xyzzyaacc1)
xyzzyaabn1=any(.not.xyzzyaacc1)
endif
if(allocated(xyzzyaacd1))deallocate(xyzzyaacd1)
xyzzyaabo1=.false.
xyzzyaabp1=.false.
if(xyzzyaaax1>0)then
allocate(xyzzyaacd1(xyzzyaaax1),stat=xyzzyaaaa75)
call check_alloc(xyzzyaaaa75,'GBASIS_RESET_ANISO','deriv_aniso_enbasis&
&')
xyzzyaacd1=.not.xyzzyaacb1
xyzzyaabo1=any(xyzzyaacd1)
xyzzyaabp1=any(.not.xyzzyaacd1)
endif
end subroutine gbasis_reset_aniso
subroutine eebasis_set_aniso(iset)
implicit none
integer,intent(in) :: iset
if(iset>0)then
xyzzyaacc1(iset)=.true.
xyzzyaabm1=any(xyzzyaacc1)
xyzzyaabn1=any(.not.xyzzyaacc1)
endif
end subroutine eebasis_set_aniso
subroutine enbasis_set_aniso(iset)
implicit none
integer,intent(in) :: iset
if(iset>0)then
xyzzyaacd1(iset)=.true.
xyzzyaabo1=any(xyzzyaacd1)
xyzzyaabp1=any(.not.xyzzyaacd1)
endif
end subroutine enbasis_set_aniso
subroutine query_cusp_eebasis(iset,ichannel,n,fake,f0,dfdr0,anisotropy&
&_index)
implicit none
integer,intent(in) :: iset,n,ichannel
real(dp),intent(inout) :: f0(n),dfdr0(n)
integer,intent(inout) :: anisotropy_index(n)
logical,intent(in) :: fake
integer xyzzyaaaa78,xyzzyaaab78,xyzzyaaac78,xyzzyaaad78,xyzzyaaae78,xy&
&zzyaaaf78,xyzzyaaag78,xyzzyaaah78,xyzzyaaai78
real(dp) xyzzyaaaj78
real(dp) :: xyzzyaaak78(1)=0.d0
real(dp),parameter :: xyzzyaaal78=0.9d0,xyzzyaaam78=1.1d0
xyzzyaaaa78=xyzzyaabq1(iset)
xyzzyaaab78=xyzzyaack1(iset)
xyzzyaaac78=xyzzyaabi1(iset)+1
xyzzyaaad78=xyzzyaabk1(iset)+1
if(xyzzyaaab78>0)then
xyzzyaaae78=xyzzyaaay1(iset)
xyzzyaaaf78=xyzzyaabg1(xyzzyaaae78+ichannel)+1
if(fake)then
xyzzyaaai78=0
if(allocated(xyzzyaabx1))xyzzyaaai78=sum(xyzzyaabx1)
xyzzyaaag78=xyzzyaaai78+sum(xyzzyaabw1)
xyzzyaaah78=xyzzyaaae78+ichannel
xyzzyaaaj78=sqrt(xyzzyaaal78**2+dble(xyzzyaaah78)*(xyzzyaaam78**2-xyzz&
&yaaal78**2)/dble(xyzzyaaag78+1))
else
xyzzyaaaj78=1.d0
endif
call xyzzyaafd1(xyzzyaaaa78,xyzzyaads1(xyzzyaaac78),             xyzzy&
&aadu1(xyzzyaaad78),xyzzyaadi1(xyzzyaaaf78),n,fake,xyzzyaaaj78,f0(1),d&
&fdr0(1),anisotropy_index)
else
call xyzzyaafd1(xyzzyaaaa78,xyzzyaads1(xyzzyaaac78),xyzzyaadu1(xyzzyaa&
&ad78),xyzzyaaak78(1),n,fake,1.d0,f0(1),dfdr0(1),anisotropy_index)
endif
end subroutine query_cusp_eebasis
subroutine query_cusp_enbasis(iset,ichannel,n,fake,f0,dfdr0,anisotropy&
&_index)
implicit none
integer,intent(in) :: iset,n,ichannel
integer,intent(inout) :: anisotropy_index(n)
real(dp),intent(inout) :: f0(n),dfdr0(n)
logical,intent(in) :: fake
integer xyzzyaaaa79,xyzzyaaab79,xyzzyaaac79,xyzzyaaad79,xyzzyaaae79,xy&
&zzyaaaf79,xyzzyaaag79,xyzzyaaah79,xyzzyaaai79
real(dp) xyzzyaaaj79
real(dp) :: xyzzyaaak79(1)=0.d0
real(dp),parameter :: xyzzyaaal79=0.9d0,xyzzyaaam79=1.1d0
xyzzyaaaa79=xyzzyaabr1(iset)
xyzzyaaab79=xyzzyaacl1(iset)
xyzzyaaac79=xyzzyaabj1(iset)+1
xyzzyaaad79=xyzzyaabl1(iset)+1
if(xyzzyaaab79>0)then
xyzzyaaae79=xyzzyaaaz1(iset)
xyzzyaaaf79=xyzzyaabh1(xyzzyaaae79+ichannel)+1
if(fake)then
xyzzyaaai79=0
if(allocated(xyzzyaabw1))xyzzyaaai79=sum(xyzzyaabw1)
xyzzyaaag79=xyzzyaaai79+sum(xyzzyaabx1)
xyzzyaaah79=xyzzyaaai79+xyzzyaaae79+ichannel
xyzzyaaaj79=sqrt(xyzzyaaal79**2+dble(xyzzyaaah79)*(xyzzyaaam79**2-xyzz&
&yaaal79**2)/dble(xyzzyaaag79+1))
else
xyzzyaaaj79=1.d0
endif
call xyzzyaafd1(xyzzyaaaa79,xyzzyaadt1(xyzzyaaac79),xyzzyaadv1(xyzzyaa&
&ad79),xyzzyaadj1(xyzzyaaaf79),n,fake,xyzzyaaaj79,f0(1),dfdr0(1),aniso&
&tropy_index)
else
call xyzzyaafd1(xyzzyaaaa79,xyzzyaadt1(xyzzyaaac79),xyzzyaadv1(xyzzyaa&
&ad79),xyzzyaaak79(1),n,fake,1.d0,f0(1),dfdr0(1),anisotropy_index)
endif
end subroutine query_cusp_enbasis
subroutine xyzzyaafd1(ibasis,ci,cd,p,n,fake,fakectt,f0,dfdr0,anisotrop&
&y_index)
implicit none
integer,intent(in) :: ibasis,n,ci(*)
integer,intent(inout) :: anisotropy_index(n)
real(dp),intent(in) :: cd(*),p(*),fakectt
real(dp),intent(inout) :: f0(n),dfdr0(n)
logical,intent(in) :: fake
f0(1:n)=1.d0
dfdr0(1:n)=0.d0
anisotropy_index=0
select case(ibasis)
case(0)
continue
case(xyzzyaaaa1)
call xyzzyaafe1(n,f0(1),dfdr0(1))
case(xyzzyaaab1)
continue
case(xyzzyaaac1)
continue
case(xyzzyaaad1)
call xyzzyaaff1(n,cd(1),p(1),fake,fakectt,f0(1),dfdr0(1))
case(xyzzyaaae1)
call xyzzyaafg1(n,p(1),fake,fakectt,f0(1),dfdr0(1))
case(xyzzyaaaf1)
call xyzzyaafh1(n,p(1),fake,fakectt,f0(1),dfdr0(1))
case(xyzzyaaag1)
call xyzzyaafi1(n,f0(1),dfdr0(1),anisotropy_index)
case(xyzzyaaah1)
call xyzzyaafj1(n,ci(1),ci(2),p(1),fake,fakectt,f0(1),dfdr0(1))
case(xyzzyaaai1)
call xyzzyaafk1(n,ci(1),ci(2),p(1),fake,fakectt,f0(1),dfdr0(1),anisotr&
&opy_index)
case(xyzzyaaaj1)
call xyzzyaafl1(p(1),fake,fakectt,f0(1),dfdr0(1))
case(xyzzyaaak1)
call xyzzyaafm1(dfdr0(1))
case(xyzzyaaal1)
call xyzzyaafn1(dfdr0(1))
case(xyzzyaaam1)
continue
case(xyzzyaaan1)
call xyzzyaafo1(dfdr0(1))
case(xyzzyaaao1)
call xyzzyaafs1(n,cd(1),f0(1),dfdr0(1))
case(xyzzyaaap1)
call xyzzyaafp1(ci(1),p(1),fake,fakectt,f0(1),dfdr0(1))
case(xyzzyaaaq1)
call xyzzyaafq1(ci(1),p(1),fake,fakectt,f0(1),dfdr0(1))
case(xyzzyaaar1)
continue
case(xyzzyaaas1)
call xyzzyaafp1(ci(1),p(1),fake,fakectt,f0(1),dfdr0(1))
anisotropy_index(1)=xyzzyaaft1(ci(2))
case(xyzzyaaat1)
call xyzzyaafr1(f0(1),dfdr0(1))
case(xyzzyaaau1)
continue
end select
end subroutine xyzzyaafd1
subroutine xyzzyaafe1(n,f0,dfdr0)
implicit none
integer,intent(in) :: n
real(dp),intent(inout) :: f0(n),dfdr0(n)
f0(1)=1.d0
f0(2:n)=0.d0
dfdr0(1:n)=0.d0
if(n>1)dfdr0(2)=1.d0
end subroutine xyzzyaafe1
subroutine xyzzyaaff1(n,l,a,fake,fakectt,f0,dfdr0)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: l,a,fakectt
real(dp),intent(inout) :: f0(n),dfdr0(n)
logical,intent(in) :: fake
real(dp) xyzzyaaaa82
f0(1)=1.d0
f0(2:n)=0.d0
dfdr0(1:n)=0.d0
if(n>1)then
if(.not.fake)then
dfdr0(2)=1.d0/a
else
xyzzyaaaa82=a*fakectt
dfdr0(2)=1.d0/xyzzyaaaa82
endif
endif
end subroutine xyzzyaaff1
subroutine xyzzyaafg1(n,p,fake,fakectt,f0,dfdr0)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: p(2),fakectt
real(dp),intent(inout) :: f0(n),dfdr0(n)
logical,intent(in) :: fake
real(dp) xyzzyaaaa83
f0(1)=1.d0
f0(2:n)=0.d0
dfdr0(1:n)=0.d0
if(n>1)then
if(.not.fake)then
dfdr0(2)=1.d0/p(1)
else
xyzzyaaaa83=p(1)*fakectt
dfdr0(2)=1.d0/xyzzyaaaa83
endif
endif
end subroutine xyzzyaafg1
subroutine xyzzyaafh1(n,a,fake,fakectt,f0,dfdr0)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: a,fakectt
real(dp),intent(inout) :: f0(n),dfdr0(n)
logical,intent(in) :: fake
real(dp) xyzzyaaaa84
integer xyzzyaaab84
f0(1)=1.d0
dfdr0(1)=0.d0
if(n>1)then
do xyzzyaaab84=2,n
if(.not.fake)then
f0(xyzzyaaab84)=1/(a**(xyzzyaaab84-1))
dfdr0(xyzzyaaab84)=-(xyzzyaaab84-1)/(a**xyzzyaaab84)
else
xyzzyaaaa84=a*fakectt
f0(xyzzyaaab84)=1/(xyzzyaaaa84**(xyzzyaaab84-1))
dfdr0(xyzzyaaab84)=-(xyzzyaaab84-1)/(xyzzyaaaa84**xyzzyaaab84)
endif
enddo
endif
end subroutine xyzzyaafh1
subroutine xyzzyaafi1(n,f0,dfdr0,anisotropy_index)
implicit none
integer,intent(in) :: n
real(dp),intent(inout) :: f0(n),dfdr0(n)
integer,intent(inout) :: anisotropy_index(n)
integer xyzzyaaaa85
f0(1:n)=0.d0
dfdr0(1:n)=0.d0
anisotropy_index(1:n)=0
if(n>0)then
do xyzzyaaaa85=1,dimensionality
f0(xyzzyaaaa85)=1.d0
dfdr0(xyzzyaaaa85)=0.d0
anisotropy_index(xyzzyaaaa85)=xyzzyaaaa85
if(n>dimensionality)then
f0(dimensionality+xyzzyaaaa85)=0.d0
dfdr0(dimensionality+xyzzyaaaa85)=1.d0
anisotropy_index(dimensionality+xyzzyaaaa85)=xyzzyaaaa85
endif
enddo
endif
end subroutine xyzzyaafi1
subroutine xyzzyaafj1(n,k0,nk,p,fake,fakectt,f0,dfdr0)
implicit none
integer,intent(in) :: n,k0,nk(*)
real(dp),intent(in) :: p(*),fakectt
real(dp),intent(inout) :: f0(n),dfdr0(n)
logical,intent(in) :: fake
f0(1:n)=0.d0
dfdr0(1:n)=0.d0
if(k0==0)then
f0(1)=1.d0
if(nk(1)>1)then
if(fake)then
dfdr0(1)=p(1)*fakectt
else
dfdr0(1)=p(1)
endif
elseif(n>1)then
dfdr0(2)=1.d0
endif
elseif(k0==1)then
dfdr0(1)=1.d0
endif
end subroutine xyzzyaafj1
subroutine xyzzyaafk1(n,k0,nk,p,fake,fakectt,f0,dfdr0,anisotropy_index&
&)
implicit none
integer,intent(in) :: n,k0,nk(*)
integer,intent(inout) :: anisotropy_index(n)
real(dp),intent(in) :: p(*),fakectt
real(dp),intent(inout) :: f0(n),dfdr0(n)
logical,intent(in) :: fake
integer xyzzyaaaa87
f0(1:n)=0.d0
dfdr0(1:n)=0.d0
anisotropy_index(1:n)=0
if(k0==0)then
do xyzzyaaaa87=1,dimensionality
f0(xyzzyaaaa87)=1.d0
anisotropy_index(xyzzyaaaa87)=xyzzyaaaa87
if(nk(1)>1)then
if(fake)then
dfdr0(xyzzyaaaa87)=p(1)*fakectt
else
dfdr0(xyzzyaaaa87)=p(1)
endif
elseif(n>dimensionality)then
dfdr0(dimensionality+xyzzyaaaa87)=1.d0
anisotropy_index(dimensionality+xyzzyaaaa87)=xyzzyaaaa87
endif
enddo
elseif(k0==1)then
do xyzzyaaaa87=1,dimensionality
dfdr0(xyzzyaaaa87)=1.d0
anisotropy_index(xyzzyaaaa87)=xyzzyaaaa87
enddo
endif
end subroutine xyzzyaafk1
subroutine xyzzyaafl1(p,fake,fakectt,f0,dfdr0)
implicit none
real(dp),intent(in) :: p(4),fakectt
real(dp),intent(inout) :: f0,dfdr0
logical,intent(in) :: fake
real(dp) xyzzyaaaa88
xyzzyaaaa88=p(1)
if(fake)xyzzyaaaa88=xyzzyaaaa88*fakectt
if(dimensionality==3)then
f0=-1.d0/xyzzyaaaa88
dfdr0=1.d0/(2*xyzzyaaaa88*xyzzyaaaa88)
else
f0=-1.d0/sqrt(xyzzyaaaa88)
dfdr0=1.d0/(3*sqrt(xyzzyaaaa88*xyzzyaaaa88*xyzzyaaaa88))
endif
end subroutine xyzzyaafl1
subroutine xyzzyaafm1(prefactor_at_zero)
implicit none
real(dp),intent(inout) :: prefactor_at_zero
prefactor_at_zero=1.d0
end subroutine xyzzyaafm1
subroutine xyzzyaafn1(prefactor_at_zero)
implicit none
real(dp),intent(inout) :: prefactor_at_zero
prefactor_at_zero=1.d0
end subroutine xyzzyaafn1
subroutine xyzzyaafo1(prefactor_at_zero)
implicit none
real(dp),intent(inout) :: prefactor_at_zero
prefactor_at_zero=1.d0
end subroutine xyzzyaafo1
subroutine xyzzyaafp1(ctrunc,p,fake,fakectt,f0,dfdr0)
implicit none
integer,intent(in) :: ctrunc
real(dp),intent(in) :: p(4),fakectt
real(dp),intent(inout) :: f0,dfdr0
logical,intent(in) :: fake
real(dp) xyzzyaaaa92
f0=1.d0
if(.not.fake)then
dfdr0=-ctrunc/p(1)
else
xyzzyaaaa92=p(1)*fakectt
dfdr0=-ctrunc/xyzzyaaaa92
endif
end subroutine xyzzyaafp1
subroutine xyzzyaafq1(ctrunc,lcut,fake,fakectt,f0,dfdr0)
implicit none
integer,intent(in) :: ctrunc
real(dp),intent(in) :: lcut,fakectt
real(dp),intent(inout) :: f0,dfdr0
logical,intent(in) :: fake
real(dp) xyzzyaaaa93
if(.not.fake)then
f0=(-lcut)**ctrunc
dfdr0=ctrunc*(-lcut)**(ctrunc-1)
else
xyzzyaaaa93=lcut*fakectt
f0=(-xyzzyaaaa93)**ctrunc
dfdr0=ctrunc*(-xyzzyaaaa93)**(ctrunc-1)
endif
end subroutine xyzzyaafq1
subroutine xyzzyaafr1(f0,dfdr0)
implicit none
real(dp),intent(inout) :: f0,dfdr0
f0=-1.d0
dfdr0=1.d0
end subroutine xyzzyaafr1
subroutine xyzzyaafs1(n,l,f0,dfdr0)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: l
real(dp),intent(inout) :: f0(n),dfdr0(n)
f0(1:n)=0.d0
dfdr0(1:n)=0.d0
dfdr0(1)=1.d0
end subroutine xyzzyaafs1
integer function xyzzyaaft1(l)
implicit none
integer,intent(in) :: l
select case(l)
case(0)
xyzzyaaft1=0
case(1)
xyzzyaaft1=1
case default
xyzzyaaft1=-1
end select
end function xyzzyaaft1
subroutine query_eqprod_ee_ee(iset,iset_cut,n,same_group,eqprod_array)
implicit none
integer,intent(in) :: iset,iset_cut,n
integer,intent(inout) :: eqprod_array(0:n,0:n)
logical,intent(in) :: same_group
integer xyzzyaaaa97,xyzzyaaab97,xyzzyaaac97,xyzzyaaad97
xyzzyaaaa97=xyzzyaabq1(iset)
if(iset_cut==0)then
xyzzyaaab97=xyzzyaaby1(iset)
else
xyzzyaaab97=0
endif
xyzzyaaac97=xyzzyaabi1(iset)+1
xyzzyaaad97=xyzzyaabk1(iset)+1
call xyzzyaafu1(xyzzyaaaa97,n,n,xyzzyaaab97,xyzzyaaab97,same_group,xyz&
&zyaads1(xyzzyaaac97),xyzzyaadu1(xyzzyaaad97),xyzzyaads1(xyzzyaaac97),&
&xyzzyaadu1(xyzzyaaad97),eqprod_array)
end subroutine query_eqprod_ee_ee
subroutine query_eqprod_en_en(iset,iset_cut,n,same_group,eqprod_array)
implicit none
integer,intent(in) :: iset,iset_cut,n
integer,intent(inout) :: eqprod_array(0:n,0:n)
logical,intent(in) :: same_group
integer xyzzyaaaa98,xyzzyaaab98,xyzzyaaac98,xyzzyaaad98
xyzzyaaaa98=xyzzyaabr1(iset)
if(iset_cut==0)then
xyzzyaaab98=xyzzyaabz1(iset)
else
xyzzyaaab98=0
endif
xyzzyaaac98=xyzzyaabj1(iset)+1
xyzzyaaad98=xyzzyaabl1(iset)+1
call xyzzyaafu1(xyzzyaaaa98,n,n,xyzzyaaab98,xyzzyaaab98,same_group,xyz&
&zyaads1(xyzzyaaac98),xyzzyaadu1(xyzzyaaad98),xyzzyaads1(xyzzyaaac98),&
&xyzzyaadu1(xyzzyaaad98),eqprod_array)
end subroutine query_eqprod_en_en
subroutine query_eqprod_ee_en(iset_ee,iset_en,iset_eecut,iset_encut,n_&
&ee,n_en,eqprod_array)
implicit none
integer,intent(in) :: iset_ee,iset_en,iset_eecut,iset_encut,n_ee,n_en
integer,intent(inout) :: eqprod_array(0:n_ee,0:n_en)
integer xyzzyaaaa99,xyzzyaaab99,xyzzyaaac99,xyzzyaaad99,xyzzyaaae99,xy&
&zzyaaaf99,xyzzyaaag99,xyzzyaaah99
xyzzyaaaa99=xyzzyaabq1(iset_ee)
xyzzyaaab99=xyzzyaabr1(iset_en)
if(iset_eecut==0)then
xyzzyaaac99=xyzzyaaby1(iset_ee)
else
xyzzyaaac99=0
endif
if(iset_encut==0)then
xyzzyaaad99=xyzzyaabz1(iset_en)
else
xyzzyaaad99=0
endif
xyzzyaaae99=xyzzyaabi1(iset_ee)+1
xyzzyaaaf99=xyzzyaabk1(iset_ee)+1
xyzzyaaag99=xyzzyaabj1(iset_en)+1
xyzzyaaah99=xyzzyaabl1(iset_en)+1
if(xyzzyaaaa99==xyzzyaaab99)then
call xyzzyaafu1(xyzzyaaaa99,n_ee,n_en,xyzzyaaac99,xyzzyaaad99,.false.,&
&xyzzyaads1(xyzzyaaae99),xyzzyaadu1(xyzzyaaaf99),xyzzyaadt1(xyzzyaaag9&
&9),xyzzyaadv1(xyzzyaaah99),eqprod_array)
else
call xyzzyaafw1(n_ee,n_en,xyzzyaaac99,xyzzyaaad99,eqprod_array)
endif
end subroutine query_eqprod_ee_en
subroutine xyzzyaafu1(ibasis,n,m,n1,m1,same_group,cin,cdn,cim,cdm,eqpr&
&od_array)
implicit none
integer,intent(in) :: ibasis,n,m,n1,m1,cin(*),cim(*)
real(dp),intent(in) :: cdn(*),cdm(*)
integer,intent(inout) :: eqprod_array(0:n,0:m)
logical,intent(in) :: same_group
eqprod_array(0:n,0:m)=0
select case(ibasis)
case(0)
eqprod_array(0:n,0:m)=1
case(xyzzyaaaa1)
call xyzzyaafv1(n,m,eqprod_array)
case(xyzzyaaab1)
if(same_group)then
call xyzzyaafx1(n,m,n1,m1,eqprod_array)
else
call xyzzyaafw1(n,m,n1,m1,eqprod_array)
endif
case(xyzzyaaac1)
if(same_group)then
call xyzzyaafx1(n,m,n1,m1,eqprod_array)
else
call xyzzyaafw1(n,m,n1,m1,eqprod_array)
endif
case(xyzzyaaad1)
if(same_group)then
call xyzzyaafv1(n,m,eqprod_array)
else
call xyzzyaafw1(n,m,n1,m1,eqprod_array)
endif
case(xyzzyaaae1)
if(same_group)then
call xyzzyaafv1(n,m,eqprod_array)
else
call xyzzyaafw1(n,m,n1,m1,eqprod_array)
endif
case(xyzzyaaaf1)
if(same_group)then
call xyzzyaafv1(n,m,eqprod_array)
else
call xyzzyaafw1(n,m,n1,m1,eqprod_array)
endif
case(xyzzyaaag1)
if(same_group)then
call xyzzyaafy1(n,m,eqprod_array)
else
call xyzzyaafw1(n,m,n1,m1,eqprod_array)
endif
case(xyzzyaaah1)
if(same_group)then
call xyzzyaafx1(n,m,n1,m1,eqprod_array)
else
call xyzzyaafw1(n,m,n1,m1,eqprod_array)
endif
case(xyzzyaaai1)
if(same_group)then
call xyzzyaafx1(n,m,n1,m1,eqprod_array)
else
call xyzzyaafw1(n,m,n1,m1,eqprod_array)
endif
case(xyzzyaaam1)
call xyzzyaafz1(n,m,cin,cim,eqprod_array)
case(xyzzyaaao1)
call xyzzyaafv1(n,m,eqprod_array)
end select
if(any(eqprod_array==0))call errstop_master('EVAL_EQPROD_GBASIS','Inte&
&rnal error: there are zeroes in EQPROD_ARRAY for basis function "'//t&
&rim(xyzzyaaav1(ibasis))//'" with SAME_GROUP='//trim(l2s(same_group))/&
&/'. Please email Pablo (pl275 at cam dot ac dot uk) about this.')
end subroutine xyzzyaafu1
subroutine xyzzyaafv1(n,m,eqprod_array)
implicit none
integer,intent(in) :: n,m
integer,intent(inout) :: eqprod_array(0:n,0:m)
integer xyzzyaaaa101,xyzzyaaab101
do xyzzyaaab101=1,m
do xyzzyaaaa101=1,n
eqprod_array(xyzzyaaaa101,xyzzyaaab101)=xyzzyaaaa101+xyzzyaaab101-1
enddo
enddo
eqprod_array(0,0)=1
do xyzzyaaaa101=1,n
eqprod_array(xyzzyaaaa101,0)=xyzzyaaaa101
enddo
do xyzzyaaab101=1,m
eqprod_array(0,xyzzyaaab101)=xyzzyaaab101
enddo
end subroutine xyzzyaafv1
subroutine xyzzyaafw1(n,m,n1,m1,eqprod_array)
implicit none
integer,intent(in) :: n,m,n1,m1
integer,intent(inout) :: eqprod_array(0:n,0:m)
integer xyzzyaaaa102,xyzzyaaab102,xyzzyaaac102
xyzzyaaac102=0
do xyzzyaaab102=1,m
do xyzzyaaaa102=1,n
xyzzyaaac102=xyzzyaaac102+1
eqprod_array(xyzzyaaaa102,xyzzyaaab102)=xyzzyaaac102
enddo
enddo
if(n1>0.and.n1<=n.and.m1>0.and.m1<=m)then
eqprod_array(0,0)=eqprod_array(n1,m1)
else
xyzzyaaac102=xyzzyaaac102+1
eqprod_array(0,0)=xyzzyaaac102
endif
if(n1>0.and.n1<=n)then
do xyzzyaaab102=1,m
eqprod_array(0,xyzzyaaab102)=eqprod_array(n1,xyzzyaaab102)
enddo
else
do xyzzyaaab102=1,m
xyzzyaaac102=xyzzyaaac102+1
eqprod_array(0,xyzzyaaab102)=xyzzyaaac102
enddo
endif
if(m1>0.and.m1<=m)then
do xyzzyaaaa102=1,n
eqprod_array(xyzzyaaaa102,0)=eqprod_array(xyzzyaaaa102,m1)
enddo
else
do xyzzyaaaa102=1,n
xyzzyaaac102=xyzzyaaac102+1
eqprod_array(xyzzyaaaa102,0)=xyzzyaaac102
enddo
endif
end subroutine xyzzyaafw1
subroutine xyzzyaafx1(n,m,n1,m1,eqprod_array)
implicit none
integer,intent(in) :: n,m,n1,m1
integer,intent(inout) :: eqprod_array(0:n,0:m)
integer xyzzyaaaa103,xyzzyaaab103,xyzzyaaac103
xyzzyaaac103=0
do xyzzyaaab103=1,min(n,m)
xyzzyaaac103=xyzzyaaac103+1
eqprod_array(xyzzyaaab103,xyzzyaaab103)=xyzzyaaac103
do xyzzyaaaa103=xyzzyaaab103+1,max(n,m)
xyzzyaaac103=xyzzyaaac103+1
if(xyzzyaaaa103<=n)eqprod_array(xyzzyaaaa103,xyzzyaaab103)=xyzzyaaac10&
&3
if(xyzzyaaaa103<=m)eqprod_array(xyzzyaaab103,xyzzyaaaa103)=xyzzyaaac10&
&3
enddo
enddo
if(n1>0.and.n1<=n.and.m1>0.and.m1<=m)then
eqprod_array(0,0)=eqprod_array(n1,m1)
else
xyzzyaaac103=xyzzyaaac103+1
eqprod_array(0,0)=xyzzyaaac103
endif
if(n1>0.and.n1<=n)then
do xyzzyaaab103=1,m
eqprod_array(0,xyzzyaaab103)=eqprod_array(n1,xyzzyaaab103)
enddo
else
do xyzzyaaab103=1,m
xyzzyaaac103=xyzzyaaac103+1
eqprod_array(0,xyzzyaaab103)=xyzzyaaac103
enddo
endif
if(m1>0.and.m1<=m)then
do xyzzyaaaa103=1,n
eqprod_array(xyzzyaaaa103,0)=eqprod_array(xyzzyaaaa103,m1)
enddo
else
do xyzzyaaaa103=1,n
xyzzyaaac103=xyzzyaaac103+1
eqprod_array(xyzzyaaaa103,0)=xyzzyaaac103
enddo
endif
end subroutine xyzzyaafx1
subroutine xyzzyaafy1(n,m,eqprod_array)
implicit none
integer,intent(in) :: n,m
integer,intent(inout) :: eqprod_array(0:n,0:m)
integer xyzzyaaaa104,xyzzyaaab104,xyzzyaaac104,xyzzyaaad104,xyzzyaaae1&
&04,xyzzyaaaf104,xyzzyaaag104,xyzzyaaah104,xyzzyaaai104,xyzzyaaaj104,x&
&yzzyaaak104
xyzzyaaaj104=m/dimensionality
xyzzyaaak104=n/dimensionality
eqprod_array(0,0)=1
xyzzyaaac104=1
do xyzzyaaaa104=1,max(n,m)
xyzzyaaac104=xyzzyaaac104+1
if(xyzzyaaaa104<=n)eqprod_array(xyzzyaaaa104,0)=xyzzyaaac104
if(xyzzyaaaa104<=m)eqprod_array(0,xyzzyaaaa104)=xyzzyaaac104
enddo
xyzzyaaaf104=xyzzyaaac104
xyzzyaaab104=0
do xyzzyaaae104=1,xyzzyaaaj104
xyzzyaaaa104=0
do xyzzyaaad104=1,xyzzyaaak104
xyzzyaaac104=xyzzyaaaf104+(xyzzyaaad104+xyzzyaaae104-2)*(((dimensional&
&ity-1)*dimensionality)/2)
xyzzyaaai104=0
do xyzzyaaag104=1,dimensionality
do xyzzyaaah104=xyzzyaaag104,dimensionality
xyzzyaaai104=xyzzyaaai104+1
eqprod_array(xyzzyaaaa104+xyzzyaaag104,xyzzyaaab104+xyzzyaaah104)=xyzz&
&yaaac104+xyzzyaaai104
eqprod_array(xyzzyaaaa104+xyzzyaaah104,xyzzyaaab104+xyzzyaaag104)=xyzz&
&yaaac104+xyzzyaaai104
enddo
enddo
xyzzyaaaa104=xyzzyaaaa104+dimensionality
enddo
xyzzyaaab104=xyzzyaaab104+dimensionality
enddo
end subroutine xyzzyaafy1
subroutine xyzzyaafz1(n,m,cin,cim,eqprod_array)
implicit none
integer,intent(in) :: n,m,cin(n),cim(m)
integer,intent(inout) :: eqprod_array(0:n,0:m)
integer xyzzyaaaa105,xyzzyaaab105,xyzzyaaac105,xyzzyaaad105,xyzzyaaae1&
&05,xyzzyaaaf105,xyzzyaaag105(0:n,0:m),xyzzyaaah105(0:n,0:m)
xyzzyaaag105=0
do xyzzyaaab105=1,m
do xyzzyaaaa105=1,n
xyzzyaaah105(xyzzyaaaa105,xyzzyaaab105)=cin(xyzzyaaaa105)+cin(m)
if(cin(xyzzyaaaa105)==0)xyzzyaaag105(xyzzyaaaa105,xyzzyaaab105)=xyzzya&
&aag105(xyzzyaaaa105,xyzzyaaab105)+1
if(cim(xyzzyaaab105)==0)xyzzyaaag105(xyzzyaaaa105,xyzzyaaab105)=xyzzya&
&aag105(xyzzyaaaa105,xyzzyaaab105)+1
enddo
enddo
xyzzyaaah105(0,0)=1
do xyzzyaaaa105=1,n
xyzzyaaah105(xyzzyaaaa105,0)=cin(xyzzyaaaa105)
if(cin(xyzzyaaaa105)==0)xyzzyaaag105(xyzzyaaaa105,xyzzyaaab105)=xyzzya&
&aag105(xyzzyaaaa105,xyzzyaaab105)+1
enddo
do xyzzyaaab105=1,m
xyzzyaaah105(0,xyzzyaaab105)=cim(xyzzyaaab105)
if(cim(xyzzyaaab105)==0)xyzzyaaag105(xyzzyaaaa105,xyzzyaaab105)=xyzzya&
&aag105(xyzzyaaaa105,xyzzyaaab105)+1
enddo
eqprod_array(0:n,0:m)=0
xyzzyaaac105=0
do xyzzyaaab105=0,m
do xyzzyaaaa105=0,n
do xyzzyaaae105=0,xyzzyaaab105
xyzzyaaaf105=n
if(xyzzyaaae105==xyzzyaaab105)xyzzyaaaf105=xyzzyaaaa105-1
if(xyzzyaaaf105<0)cycle
do xyzzyaaad105=0,xyzzyaaaf105
if(xyzzyaaah105(xyzzyaaaa105,xyzzyaaab105)==xyzzyaaah105(xyzzyaaad105,&
&xyzzyaaae105).and.xyzzyaaag105(xyzzyaaaa105,xyzzyaaab105)==xyzzyaaag1&
&05(xyzzyaaad105,xyzzyaaae105))exit
enddo
if(xyzzyaaad105<=xyzzyaaaf105)exit
enddo
if(xyzzyaaae105<=xyzzyaaab105)then
eqprod_array(xyzzyaaaa105,xyzzyaaab105)=eqprod_array(xyzzyaaad105,xyzz&
&yaaae105)
else
xyzzyaaac105=xyzzyaaac105+1
eqprod_array(xyzzyaaaa105,xyzzyaaab105)=xyzzyaaac105
endif
enddo
enddo
end subroutine xyzzyaafz1
subroutine xyzzyaaga1(ii,nzstride,eevecs1,eebasis1,nzeecut1)
implicit none
integer,intent(in) :: ii,nzstride
real(dp),intent(in) :: eevecs1(4,netot)
real(dp),intent(inout) :: eebasis1(nfn_eebasis,netot)
logical,intent(inout) :: nzeecut1(netot,nzstride,*)
integer xyzzyaaaa106,xyzzyaaab106,xyzzyaaac106,xyzzyaaad106,xyzzyaaae1&
&06,xyzzyaaaf106,xyzzyaaag106,xyzzyaaah106,xyzzyaaai106,xyzzyaaaj106,x&
&yzzyaaak106,xyzzyaaal106,xyzzyaaam106,xyzzyaaan106,xyzzyaaao106
real(dp) xyzzyaaap106(1)
xyzzyaaab106=which_spin(ii)
xyzzyaaaa106=which_ie(ii)
do xyzzyaaac106=1,xyzzyaaaw1
xyzzyaaad106=xyzzyaack1(xyzzyaaac106)
xyzzyaaae106=xyzzyaabq1(xyzzyaaac106)
xyzzyaaaf106=xyzzyaabs1(xyzzyaaac106)
xyzzyaaag106=ifn1_eebasis(xyzzyaaac106)
xyzzyaaai106=xyzzyaabk1(xyzzyaaac106)+1
xyzzyaaaj106=xyzzyaabi1(xyzzyaaac106)+1
if(xyzzyaaad106>0)then
xyzzyaaah106=xyzzyaaay1(xyzzyaaac106)
xyzzyaaak106=1
do xyzzyaaal106=1,nspin
if(.not.(nele(xyzzyaaal106)==0.or.((xyzzyaaal106==xyzzyaaab106).and.ne&
&le(xyzzyaaab106)==1)))then
if(xyzzyaaal106==xyzzyaaab106)then
xyzzyaaam106=xyzzyaaaa106
elseif(xyzzyaaal106<xyzzyaaab106)then
xyzzyaaam106=nele(xyzzyaaal106)+1
else
xyzzyaaam106=0
endif
xyzzyaaao106=xyzzyaabu1(xyzzyaaab106,xyzzyaaal106,xyzzyaaac106)
if(xyzzyaaao106>0)then
xyzzyaaan106=xyzzyaabg1(xyzzyaaah106+xyzzyaaao106)+1
call xyzzyaagk1(xyzzyaaaf106,xyzzyaads1(xyzzyaaaj106),           xyzzy&
&aadu1(xyzzyaaai106),xyzzyaadi1(xyzzyaaan106),nele(xyzzyaaal106),xyzzy&
&aaam106,0,eevecs1(1,xyzzyaaak106),xyzzyaaae106,nfn_eebasis,eebasis1(x&
&yzzyaaag106,xyzzyaaak106),      nzeecut1(xyzzyaaak106,1,xyzzyaaac106)&
&)
endif
endif
xyzzyaaak106=xyzzyaaak106+nele(xyzzyaaal106)
enddo
else
call xyzzyaagk1(xyzzyaaaf106,xyzzyaads1(xyzzyaaaj106),       xyzzyaadu&
&1(xyzzyaaai106),xyzzyaaap106(1),netot,ii,0,eevecs1(1,1),xyzzyaaae106,&
&nfn_eebasis,eebasis1(xyzzyaaag106,1),nzeecut1(1,1,xyzzyaaac106))
endif
enddo
end subroutine xyzzyaaga1
subroutine xyzzyaagb1(ii,nzstride,eivecs1,enbasis1,nzencut1)
implicit none
integer,intent(in) :: ii,nzstride
real(dp),intent(in) :: eivecs1(4,nitot)
real(dp),intent(inout) :: enbasis1(nfn_enbasis,nitot)
logical,intent(inout) :: nzencut1(nitot,nzstride,*)
integer xyzzyaaaa107,xyzzyaaab107,xyzzyaaac107,xyzzyaaad107,xyzzyaaae1&
&07,xyzzyaaaf107,xyzzyaaag107,xyzzyaaah107,xyzzyaaai107,xyzzyaaaj107,x&
&yzzyaaak107,xyzzyaaal107
real(dp) xyzzyaaam107(1)
xyzzyaaaa107=which_spin(ii)
do xyzzyaaab107=1,xyzzyaaax1
xyzzyaaac107=xyzzyaacl1(xyzzyaaab107)
xyzzyaaad107=xyzzyaabr1(xyzzyaaab107)
xyzzyaaae107=xyzzyaabt1(xyzzyaaab107)
xyzzyaaaf107=ifn1_enbasis(xyzzyaaab107)
xyzzyaaah107=xyzzyaabl1(xyzzyaaab107)+1
xyzzyaaai107=xyzzyaabj1(xyzzyaaab107)+1
if(xyzzyaaac107>0)then
xyzzyaaag107=xyzzyaaaz1(xyzzyaaab107)
do xyzzyaaak107=1,nitot
xyzzyaaal107=xyzzyaabv1(xyzzyaaak107,xyzzyaaaa107,xyzzyaaab107)
if(xyzzyaaal107<1)cycle
xyzzyaaaj107=xyzzyaabh1(xyzzyaaag107+xyzzyaaal107)+1
call xyzzyaagk1(xyzzyaaae107,xyzzyaadt1(xyzzyaaai107),xyzzyaadv1(xyzzy&
&aaah107),xyzzyaadj1(xyzzyaaaj107),1,0,xyzzyaaak107,eivecs1(1,xyzzyaaa&
&k107),xyzzyaaad107,nfn_enbasis,enbasis1(xyzzyaaaf107,xyzzyaaak107),nz&
&encut1(xyzzyaaak107,1,xyzzyaaab107))
enddo
else
call xyzzyaagk1(xyzzyaaae107,xyzzyaadt1(xyzzyaaai107),xyzzyaadv1(xyzzy&
&aaah107),xyzzyaaam107(1),nitot,0,1,eivecs1(1,1),xyzzyaaad107,nfn_enba&
&sis,enbasis1(xyzzyaaaf107,1),nzencut1(1,1,xyzzyaaab107))
endif
enddo
end subroutine xyzzyaagb1
subroutine xyzzyaagc1(ii,nzstride,eevecs1,eebasis1,grad_eebasis1,deeba&
&sis1,gradr_eebasis1,nzeecut1)
implicit none
integer,intent(in) :: ii,nzstride
real(dp),intent(in) :: eevecs1(4,netot)
real(dp),intent(inout) :: eebasis1(nfn_eebasis,netot),grad_eebasis1(3,&
&nfn_eebasis,netot),deebasis1(nfn_eebasis,netot),gradr_eebasis1(3,neto&
&t)
logical,intent(inout) :: nzeecut1(netot,nzstride,*)
integer xyzzyaaaa108,xyzzyaaab108,xyzzyaaac108,xyzzyaaad108,xyzzyaaae1&
&08,xyzzyaaaf108,xyzzyaaag108,xyzzyaaah108,xyzzyaaai108,xyzzyaaaj108,x&
&yzzyaaak108,xyzzyaaal108,xyzzyaaam108,xyzzyaaan108,xyzzyaaao108
logical xyzzyaaap108
real(dp) xyzzyaaaq108(1)
xyzzyaaab108=which_spin(ii)
xyzzyaaaa108=which_ie(ii)
do xyzzyaaak108=1,netot
if(ii==xyzzyaaak108)cycle
call xyzzyaaiz1(eevecs1(4,xyzzyaaak108),eevecs1(1,xyzzyaaak108),gradr_&
&eebasis1(1,xyzzyaaak108))
enddo
do xyzzyaaac108=1,xyzzyaaaw1
xyzzyaaap108=xyzzyaaca1(xyzzyaaac108).and.xyzzyaacc1(xyzzyaaac108)
xyzzyaaad108=xyzzyaack1(xyzzyaaac108)
xyzzyaaae108=xyzzyaabq1(xyzzyaaac108)
xyzzyaaaf108=xyzzyaabs1(xyzzyaaac108)
xyzzyaaag108=ifn1_eebasis(xyzzyaaac108)
xyzzyaaai108=xyzzyaabk1(xyzzyaaac108)+1
xyzzyaaaj108=xyzzyaabi1(xyzzyaaac108)+1
if(xyzzyaaad108>0)then
xyzzyaaah108=xyzzyaaay1(xyzzyaaac108)
xyzzyaaak108=1
do xyzzyaaal108=1,nspin
if(.not.(nele(xyzzyaaal108)==0.or.((xyzzyaaal108==xyzzyaaab108).and.ne&
&le(xyzzyaaab108)==1)))then
if(xyzzyaaal108==xyzzyaaab108)then
xyzzyaaam108=xyzzyaaaa108
elseif(xyzzyaaal108<xyzzyaaab108)then
xyzzyaaam108=nele(xyzzyaaal108)+1
else
xyzzyaaam108=0
endif
xyzzyaaao108=xyzzyaabu1(xyzzyaaab108,xyzzyaaal108,xyzzyaaac108)
if(xyzzyaaao108>0)then
xyzzyaaan108=xyzzyaabg1(xyzzyaaah108+xyzzyaaao108)+1
call xyzzyaagl1(xyzzyaaaf108,xyzzyaads1(xyzzyaaaj108),xyzzyaadu1(xyzzy&
&aaai108),xyzzyaadi1(xyzzyaaan108),nele(xyzzyaaal108),xyzzyaaam108,0,e&
&evecs1(1,xyzzyaaak108),xyzzyaaae108,nfn_eebasis,eebasis1(xyzzyaaag108&
&,xyzzyaaak108),grad_eebasis1(1,xyzzyaaag108,xyzzyaaak108),deebasis1(x&
&yzzyaaag108,xyzzyaaak108),gradr_eebasis1(1,xyzzyaaak108),nzeecut1(xyz&
&zyaaak108,1,xyzzyaaac108),xyzzyaaap108)
endif
endif
xyzzyaaak108=xyzzyaaak108+nele(xyzzyaaal108)
enddo
else
call xyzzyaagl1(xyzzyaaaf108,xyzzyaads1(xyzzyaaaj108),xyzzyaadu1(xyzzy&
&aaai108),xyzzyaaaq108(1),netot,ii,0,eevecs1(1,1),xyzzyaaae108,nfn_eeb&
&asis,eebasis1(xyzzyaaag108,1),grad_eebasis1(1,xyzzyaaag108,1),  deeba&
&sis1(xyzzyaaag108,1),gradr_eebasis1(1,1),nzeecut1(1,1,xyzzyaaac108),x&
&yzzyaaap108)
endif
enddo
end subroutine xyzzyaagc1
subroutine xyzzyaagd1(ii,nzstride,eivecs1,enbasis1,grad_enbasis1,denba&
&sis1,gradr_enbasis1,nzencut1)
implicit none
integer,intent(in) :: ii,nzstride
real(dp),intent(in) :: eivecs1(4,nitot)
real(dp),intent(inout) :: enbasis1(nfn_enbasis,nitot),grad_enbasis1(3,&
&nfn_enbasis,nitot),denbasis1(nfn_enbasis,nitot),gradr_enbasis1(3,nito&
&t)
logical,intent(inout) :: nzencut1(nitot,nzstride,*)
integer xyzzyaaaa109,xyzzyaaab109,xyzzyaaac109,xyzzyaaad109,xyzzyaaae1&
&09,xyzzyaaaf109,xyzzyaaag109,xyzzyaaah109,xyzzyaaai109,xyzzyaaaj109,x&
&yzzyaaak109,xyzzyaaal109
logical xyzzyaaam109
real(dp) xyzzyaaan109(1)
xyzzyaaaa109=which_spin(ii)
do xyzzyaaak109=1,nitot
call xyzzyaaiz1(eivecs1(4,xyzzyaaak109),eivecs1(1,xyzzyaaak109),gradr_&
&enbasis1(1,xyzzyaaak109))
enddo
do xyzzyaaab109=1,xyzzyaaax1
xyzzyaaam109=xyzzyaacb1(xyzzyaaab109).and.xyzzyaacd1(xyzzyaaab109)
xyzzyaaac109=xyzzyaacl1(xyzzyaaab109)
xyzzyaaad109=xyzzyaabr1(xyzzyaaab109)
xyzzyaaae109=xyzzyaabt1(xyzzyaaab109)
xyzzyaaaf109=ifn1_enbasis(xyzzyaaab109)
xyzzyaaah109=xyzzyaabl1(xyzzyaaab109)+1
xyzzyaaai109=xyzzyaabj1(xyzzyaaab109)+1
if(xyzzyaaac109>0)then
xyzzyaaag109=xyzzyaaaz1(xyzzyaaab109)
do xyzzyaaak109=1,nitot
xyzzyaaal109=xyzzyaabv1(xyzzyaaak109,xyzzyaaaa109,xyzzyaaab109)
if(xyzzyaaal109<1)cycle
xyzzyaaaj109=xyzzyaabh1(xyzzyaaag109+xyzzyaaal109)+1
call xyzzyaagl1(xyzzyaaae109,xyzzyaadt1(xyzzyaaai109),xyzzyaadv1(xyzzy&
&aaah109),xyzzyaadj1(xyzzyaaaj109),1,0,xyzzyaaak109,eivecs1(1,xyzzyaaa&
&k109),xyzzyaaad109,nfn_enbasis,enbasis1(xyzzyaaaf109,xyzzyaaak109),gr&
&ad_enbasis1(1,xyzzyaaaf109,xyzzyaaak109),denbasis1(xyzzyaaaf109,xyzzy&
&aaak109),gradr_enbasis1(1,xyzzyaaak109),nzencut1(xyzzyaaak109,1,xyzzy&
&aaab109),xyzzyaaam109)
enddo
else
call xyzzyaagl1(xyzzyaaae109,xyzzyaadt1(xyzzyaaai109),xyzzyaadv1(xyzzy&
&aaah109),xyzzyaaan109(1),nitot,0,1,eivecs1(1,1),xyzzyaaad109,nfn_enba&
&sis,enbasis1(xyzzyaaaf109,1),grad_enbasis1(1,xyzzyaaaf109,1),denbasis&
&1(xyzzyaaaf109,1),gradr_enbasis1(1,1),nzencut1(1,1,xyzzyaaab109),xyzz&
&yaaam109)
endif
enddo
end subroutine xyzzyaagd1
subroutine xyzzyaage1(eevecs,eebasis,nzeecut)
implicit none
real(dp),intent(in) :: eevecs(4,netot,netot)
real(dp),intent(inout) :: eebasis(nfn_eebasis,netot,netot)
logical,intent(inout) :: nzeecut(netot,netot,xyzzyaaaw1)
integer xyzzyaaaa110,xyzzyaaab110,xyzzyaaac110,xyzzyaaad110,xyzzyaaae1&
&10,xyzzyaaaf110,xyzzyaaag110,xyzzyaaah110,xyzzyaaai110,xyzzyaaaj110,x&
&yzzyaaak110,xyzzyaaal110,xyzzyaaam110,xyzzyaaan110,xyzzyaaao110,xyzzy&
&aaap110,xyzzyaaaq110
real(dp) xyzzyaaar110(1)
do xyzzyaaad110=1,xyzzyaaaw1
xyzzyaaah110=ifn1_eebasis(xyzzyaaad110)
xyzzyaaae110=xyzzyaack1(xyzzyaaad110)
xyzzyaaaf110=xyzzyaabq1(xyzzyaaad110)
xyzzyaaag110=xyzzyaabs1(xyzzyaaad110)
xyzzyaaaj110=xyzzyaabk1(xyzzyaaad110)+1
xyzzyaaak110=xyzzyaabi1(xyzzyaaad110)+1
if(xyzzyaaae110>0)then
xyzzyaaai110=xyzzyaaay1(xyzzyaaad110)
xyzzyaaap110=0
do xyzzyaaac110=1,nspin
if(nele(xyzzyaaac110)==0)cycle
if(nele(xyzzyaaac110)>1)then
xyzzyaaao110=xyzzyaabu1(xyzzyaaac110,xyzzyaaac110,xyzzyaaad110)
if(xyzzyaaao110>0)then
xyzzyaaan110=xyzzyaabg1(xyzzyaaai110+xyzzyaaao110)+1
xyzzyaaaa110=xyzzyaaap110
xyzzyaaal110=xyzzyaaap110+1
xyzzyaaaq110=nele(xyzzyaaac110)
do xyzzyaaab110=1,nele(xyzzyaaac110)-1
xyzzyaaaa110=xyzzyaaaa110+1
xyzzyaaal110=xyzzyaaal110+1
xyzzyaaaq110=xyzzyaaaq110-1
call xyzzyaagk1(xyzzyaaag110,xyzzyaads1(xyzzyaaak110),xyzzyaadu1(xyzzy&
&aaaj110),xyzzyaadi1(xyzzyaaan110),xyzzyaaaq110,0,0,eevecs(1,xyzzyaaal&
&110,xyzzyaaaa110),xyzzyaaaf110,nfn_eebasis,eebasis(xyzzyaaah110,xyzzy&
&aaal110,xyzzyaaaa110),nzeecut(xyzzyaaal110,xyzzyaaaa110,xyzzyaaad110)&
&)
enddo
endif
endif
xyzzyaaal110=xyzzyaaap110+nele(xyzzyaaac110)+1
do xyzzyaaam110=xyzzyaaac110+1,nspin
if(nele(xyzzyaaam110)==0)cycle
xyzzyaaao110=xyzzyaabu1(xyzzyaaac110,xyzzyaaam110,xyzzyaaad110)
if(xyzzyaaao110>0)then
xyzzyaaan110=xyzzyaabg1(xyzzyaaai110+xyzzyaaao110)+1
xyzzyaaaa110=xyzzyaaap110
xyzzyaaaq110=nele(xyzzyaaam110)
do xyzzyaaab110=1,nele(xyzzyaaac110)
xyzzyaaaa110=xyzzyaaaa110+1
call xyzzyaagk1(xyzzyaaag110,xyzzyaads1(xyzzyaaak110),xyzzyaadu1(xyzzy&
&aaaj110),xyzzyaadi1(xyzzyaaan110),xyzzyaaaq110,0,0,eevecs(1,xyzzyaaal&
&110,xyzzyaaaa110),xyzzyaaaf110,nfn_eebasis,eebasis(xyzzyaaah110,xyzzy&
&aaal110,xyzzyaaaa110),nzeecut(xyzzyaaal110,xyzzyaaaa110,xyzzyaaad110)&
&)
enddo
endif
xyzzyaaal110=xyzzyaaal110+nele(xyzzyaaam110)
enddo
xyzzyaaap110=xyzzyaaap110+nele(xyzzyaaac110)
enddo
else
xyzzyaaal110=1
xyzzyaaaq110=netot
do xyzzyaaaa110=1,netot-1
xyzzyaaal110=xyzzyaaal110+1
xyzzyaaaq110=xyzzyaaaq110-1
call xyzzyaagk1(xyzzyaaag110,xyzzyaads1(xyzzyaaak110),xyzzyaadu1(xyzzy&
&aaaj110),xyzzyaaar110(1),xyzzyaaaq110,0,0,eevecs(1,xyzzyaaal110,xyzzy&
&aaaa110),xyzzyaaaf110,nfn_eebasis,eebasis(xyzzyaaah110,xyzzyaaal110,x&
&yzzyaaaa110),nzeecut(xyzzyaaal110,xyzzyaaaa110,xyzzyaaad110))
enddo
endif
enddo
do xyzzyaaaa110=1,netot-1
eebasis(:,xyzzyaaaa110,xyzzyaaaa110+1:netot)=eebasis(:,xyzzyaaaa110+1:&
&netot,xyzzyaaaa110)
nzeecut(xyzzyaaaa110,xyzzyaaaa110+1:netot,:)=nzeecut(xyzzyaaaa110+1:ne&
&tot,xyzzyaaaa110,:)
enddo
end subroutine xyzzyaage1
subroutine xyzzyaagf1(eivecs,enbasis,nzencut)
implicit none
real(dp),intent(in) :: eivecs(4,nitot,netot)
logical,intent(inout) :: nzencut(nitot,netot,xyzzyaaax1)
real(dp),intent(inout) :: enbasis(nfn_enbasis,nitot,netot)
integer xyzzyaaaa111,xyzzyaaab111,xyzzyaaac111,xyzzyaaad111,xyzzyaaae1&
&11,xyzzyaaaf111,xyzzyaaag111,xyzzyaaah111,xyzzyaaai111,xyzzyaaaj111,x&
&yzzyaaak111,xyzzyaaal111,xyzzyaaam111,xyzzyaaan111
real(dp) xyzzyaaao111(1)
do xyzzyaaad111=1,xyzzyaaax1
xyzzyaaah111=ifn1_enbasis(xyzzyaaad111)
xyzzyaaae111=xyzzyaacl1(xyzzyaaad111)
xyzzyaaaf111=xyzzyaabr1(xyzzyaaad111)
xyzzyaaag111=xyzzyaabt1(xyzzyaaad111)
xyzzyaaaj111=xyzzyaabl1(xyzzyaaad111)+1
xyzzyaaak111=xyzzyaabj1(xyzzyaaad111)+1
if(xyzzyaaae111>0)then
xyzzyaaai111=xyzzyaaaz1(xyzzyaaad111)
do xyzzyaaan111=1,nitot
xyzzyaaaa111=0
do xyzzyaaac111=1,nspin
if(nele(xyzzyaaac111)==0)cycle
xyzzyaaam111=xyzzyaabv1(xyzzyaaan111,xyzzyaaac111,xyzzyaaad111)
if(xyzzyaaam111<1)then
xyzzyaaaa111=xyzzyaaaa111+nele(xyzzyaaac111)
cycle
endif
xyzzyaaal111=xyzzyaabh1(xyzzyaaai111+xyzzyaaam111)+1
do xyzzyaaab111=1,nele(xyzzyaaac111)
xyzzyaaaa111=xyzzyaaaa111+1
call xyzzyaagk1(xyzzyaaag111,xyzzyaadt1(xyzzyaaak111),xyzzyaadv1(xyzzy&
&aaaj111),xyzzyaadj1(xyzzyaaal111),1,0,xyzzyaaan111,eivecs(1,xyzzyaaan&
&111,xyzzyaaaa111),xyzzyaaaf111,nfn_enbasis,enbasis(xyzzyaaah111,xyzzy&
&aaan111,xyzzyaaaa111),nzencut(xyzzyaaan111,xyzzyaaaa111,xyzzyaaad111)&
&)
enddo
enddo
enddo
else
do xyzzyaaaa111=1,netot
call xyzzyaagk1(xyzzyaaag111,xyzzyaadt1(xyzzyaaak111),xyzzyaadv1(xyzzy&
&aaaj111),xyzzyaaao111(1),nitot,0,1,eivecs(1,1,xyzzyaaaa111),xyzzyaaaf&
&111,nfn_enbasis,enbasis(xyzzyaaah111,1,xyzzyaaaa111),nzencut(1,xyzzya&
&aaa111,xyzzyaaad111))
enddo
endif
enddo
end subroutine xyzzyaagf1
subroutine xyzzyaagg1(eevecs,eebasis,grad_eebasis,deebasis,gradr_eebas&
&is,nzeecut)
implicit none
real(dp),intent(in) :: eevecs(4,netot,netot)
logical,intent(inout) :: nzeecut(netot,netot,xyzzyaaaw1)
real(dp),intent(inout) :: eebasis(nfn_eebasis,netot,netot),grad_eebasi&
&s(3,nfn_eebasis,netot,netot),deebasis(nfn_eebasis,netot,netot),gradr_&
&eebasis(3,netot,netot)
integer xyzzyaaaa112,xyzzyaaab112,xyzzyaaac112,xyzzyaaad112,xyzzyaaae1&
&12,xyzzyaaaf112,xyzzyaaag112,xyzzyaaah112,xyzzyaaai112,xyzzyaaaj112,x&
&yzzyaaak112,xyzzyaaal112,xyzzyaaam112,xyzzyaaan112,xyzzyaaao112,xyzzy&
&aaap112,xyzzyaaaq112
logical xyzzyaaar112
real(dp) xyzzyaaas112(1)
do xyzzyaaaa112=1,netot
do xyzzyaaal112=xyzzyaaaa112+1,netot
call xyzzyaaiz1(eevecs(4,xyzzyaaal112,xyzzyaaaa112),eevecs(1,xyzzyaaal&
&112,xyzzyaaaa112),gradr_eebasis(1,xyzzyaaal112,xyzzyaaaa112))
gradr_eebasis(1:3,xyzzyaaaa112,xyzzyaaal112)=-gradr_eebasis(1:3,xyzzya&
&aal112,xyzzyaaaa112)
enddo
enddo
do xyzzyaaad112=1,xyzzyaaaw1
xyzzyaaar112=xyzzyaaca1(xyzzyaaad112).and.xyzzyaacc1(xyzzyaaad112)
xyzzyaaah112=ifn1_eebasis(xyzzyaaad112)
xyzzyaaae112=xyzzyaack1(xyzzyaaad112)
xyzzyaaaf112=xyzzyaabq1(xyzzyaaad112)
xyzzyaaag112=xyzzyaabs1(xyzzyaaad112)
xyzzyaaaj112=xyzzyaabk1(xyzzyaaad112)+1
xyzzyaaak112=xyzzyaabi1(xyzzyaaad112)+1
if(xyzzyaaae112>0)then
xyzzyaaai112=xyzzyaaay1(xyzzyaaad112)
xyzzyaaap112=0
do xyzzyaaac112=1,nspin
if(nele(xyzzyaaac112)==0)cycle
if(nele(xyzzyaaac112)>1)then
xyzzyaaao112=xyzzyaabu1(xyzzyaaac112,xyzzyaaac112,xyzzyaaad112)
if(xyzzyaaao112>0)then
xyzzyaaan112=xyzzyaabg1(xyzzyaaai112+xyzzyaaao112)+1
xyzzyaaaa112=xyzzyaaap112
xyzzyaaal112=xyzzyaaap112+1
xyzzyaaaq112=nele(xyzzyaaac112)
do xyzzyaaab112=1,nele(xyzzyaaac112)-1
xyzzyaaaa112=xyzzyaaaa112+1
xyzzyaaal112=xyzzyaaal112+1
xyzzyaaaq112=xyzzyaaaq112-1
call xyzzyaagl1(xyzzyaaag112,xyzzyaads1(xyzzyaaak112),xyzzyaadu1(xyzzy&
&aaaj112),xyzzyaadi1(xyzzyaaan112),xyzzyaaaq112,0,0,eevecs(1,xyzzyaaal&
&112,xyzzyaaaa112),xyzzyaaaf112,nfn_eebasis,eebasis(xyzzyaaah112,xyzzy&
&aaal112,xyzzyaaaa112),grad_eebasis(1,xyzzyaaah112,xyzzyaaal112,xyzzya&
&aaa112),deebasis(xyzzyaaah112,xyzzyaaal112,xyzzyaaaa112),gradr_eebasi&
&s(1,xyzzyaaal112,xyzzyaaaa112),nzeecut(xyzzyaaal112,xyzzyaaaa112,xyzz&
&yaaad112),xyzzyaaar112)
enddo
endif
endif
xyzzyaaal112=xyzzyaaap112+nele(xyzzyaaac112)+1
do xyzzyaaam112=xyzzyaaac112+1,nspin
if(nele(xyzzyaaam112)==0)cycle
xyzzyaaao112=xyzzyaabu1(xyzzyaaac112,xyzzyaaam112,xyzzyaaad112)
if(xyzzyaaao112>0)then
xyzzyaaan112=xyzzyaabg1(xyzzyaaai112+xyzzyaaao112)+1
xyzzyaaaa112=xyzzyaaap112
xyzzyaaaq112=nele(xyzzyaaam112)
do xyzzyaaab112=1,nele(xyzzyaaac112)
xyzzyaaaa112=xyzzyaaaa112+1
call xyzzyaagl1(xyzzyaaag112,xyzzyaads1(xyzzyaaak112),xyzzyaadu1(xyzzy&
&aaaj112),xyzzyaadi1(xyzzyaaan112),xyzzyaaaq112,0,0,eevecs(1,xyzzyaaal&
&112,xyzzyaaaa112),xyzzyaaaf112,nfn_eebasis,eebasis(xyzzyaaah112,xyzzy&
&aaal112,xyzzyaaaa112),grad_eebasis(1,xyzzyaaah112,xyzzyaaal112,xyzzya&
&aaa112),deebasis(xyzzyaaah112,xyzzyaaal112,xyzzyaaaa112),gradr_eebasi&
&s(1,xyzzyaaal112,xyzzyaaaa112),nzeecut(xyzzyaaal112,xyzzyaaaa112,xyzz&
&yaaad112),xyzzyaaar112)
enddo
endif
xyzzyaaal112=xyzzyaaal112+nele(xyzzyaaam112)
enddo
xyzzyaaap112=xyzzyaaap112+nele(xyzzyaaac112)
enddo
else
xyzzyaaal112=1
xyzzyaaaq112=netot
do xyzzyaaaa112=1,netot-1
xyzzyaaal112=xyzzyaaal112+1
xyzzyaaaq112=xyzzyaaaq112-1
call xyzzyaagl1(xyzzyaaag112,xyzzyaads1(xyzzyaaak112),xyzzyaadu1(xyzzy&
&aaaj112),xyzzyaaas112(1),xyzzyaaaq112,0,0,eevecs(1,xyzzyaaal112,xyzzy&
&aaaa112),xyzzyaaaf112,nfn_eebasis,eebasis(xyzzyaaah112,xyzzyaaal112,x&
&yzzyaaaa112),grad_eebasis(1,xyzzyaaah112,xyzzyaaal112,xyzzyaaaa112),d&
&eebasis(xyzzyaaah112,xyzzyaaal112,xyzzyaaaa112),gradr_eebasis(1,xyzzy&
&aaal112,xyzzyaaaa112),nzeecut(xyzzyaaal112,xyzzyaaaa112,xyzzyaaad112)&
&,xyzzyaaar112)
enddo
endif
enddo
do xyzzyaaaa112=1,netot-1
eebasis(:,xyzzyaaaa112,xyzzyaaaa112+1:netot)=eebasis(:,xyzzyaaaa112+1:&
&netot,xyzzyaaaa112)
if(xyzzyaabm1)grad_eebasis(:,:,xyzzyaaaa112,xyzzyaaaa112+1:netot)=-gra&
&d_eebasis(:,:,xyzzyaaaa112+1:netot,xyzzyaaaa112)
if(xyzzyaabn1)deebasis(:,xyzzyaaaa112,xyzzyaaaa112+1:netot)=deebasis(:&
&,xyzzyaaaa112+1:netot,xyzzyaaaa112)
nzeecut(xyzzyaaaa112,xyzzyaaaa112+1:netot,:)=nzeecut(xyzzyaaaa112+1:ne&
&tot,xyzzyaaaa112,:)
enddo
end subroutine xyzzyaagg1
subroutine xyzzyaagh1(eivecs,enbasis,grad_enbasis,denbasis,gradr_enbas&
&is,nzencut)
implicit none
real(dp),intent(in) :: eivecs(4,nitot,netot)
real(dp),intent(inout) :: enbasis(nfn_enbasis,nitot,netot),grad_enbasi&
&s(3,nfn_enbasis,nitot,netot),denbasis(nfn_enbasis,nitot,netot),gradr_&
&enbasis(3,nitot,netot)
logical,intent(inout) :: nzencut(nitot,netot,xyzzyaaax1)
integer xyzzyaaaa113,xyzzyaaab113,xyzzyaaac113,xyzzyaaad113,xyzzyaaae1&
&13,xyzzyaaaf113,xyzzyaaag113,xyzzyaaah113,xyzzyaaai113,xyzzyaaaj113,x&
&yzzyaaak113,xyzzyaaal113,xyzzyaaam113,xyzzyaaan113
logical xyzzyaaao113
real(dp) xyzzyaaap113(1)
do xyzzyaaaa113=1,netot
do xyzzyaaan113=1,nitot
call xyzzyaaiz1(eivecs(4,xyzzyaaan113,xyzzyaaaa113),eivecs(1,xyzzyaaan&
&113,xyzzyaaaa113),gradr_enbasis(1,xyzzyaaan113,xyzzyaaaa113))
enddo
enddo
do xyzzyaaad113=1,xyzzyaaax1
xyzzyaaao113=xyzzyaacb1(xyzzyaaad113).and.xyzzyaacd1(xyzzyaaad113)
xyzzyaaah113=ifn1_enbasis(xyzzyaaad113)
xyzzyaaae113=xyzzyaacl1(xyzzyaaad113)
xyzzyaaaf113=xyzzyaabr1(xyzzyaaad113)
xyzzyaaag113=xyzzyaabt1(xyzzyaaad113)
xyzzyaaaj113=xyzzyaabl1(xyzzyaaad113)+1
xyzzyaaak113=xyzzyaabj1(xyzzyaaad113)+1
if(xyzzyaaae113>0)then
xyzzyaaai113=xyzzyaaaz1(xyzzyaaad113)
do xyzzyaaan113=1,nitot
xyzzyaaaa113=0
do xyzzyaaac113=1,nspin
if(nele(xyzzyaaac113)==0)cycle
xyzzyaaam113=xyzzyaabv1(xyzzyaaan113,xyzzyaaac113,xyzzyaaad113)
if(xyzzyaaam113<1)then
xyzzyaaaa113=xyzzyaaaa113+nele(xyzzyaaac113)
cycle
endif
xyzzyaaal113=xyzzyaabh1(xyzzyaaai113+xyzzyaaam113)+1
do xyzzyaaab113=1,nele(xyzzyaaac113)
xyzzyaaaa113=xyzzyaaaa113+1
call xyzzyaagl1(xyzzyaaag113,xyzzyaadt1(xyzzyaaak113),xyzzyaadv1(xyzzy&
&aaaj113),xyzzyaadj1(xyzzyaaal113),1,0,xyzzyaaan113,eivecs(1,xyzzyaaan&
&113,xyzzyaaaa113),xyzzyaaaf113,nfn_enbasis,enbasis(xyzzyaaah113,xyzzy&
&aaan113,xyzzyaaaa113),grad_enbasis(1,xyzzyaaah113,xyzzyaaan113,xyzzya&
&aaa113),denbasis(xyzzyaaah113,xyzzyaaan113,xyzzyaaaa113),gradr_enbasi&
&s(1,xyzzyaaan113,xyzzyaaaa113),nzencut(xyzzyaaan113,xyzzyaaaa113,xyzz&
&yaaad113),xyzzyaaao113)
enddo
enddo
enddo
else
do xyzzyaaaa113=1,netot
call xyzzyaagl1(xyzzyaaag113,xyzzyaadt1(xyzzyaaak113),xyzzyaadv1(xyzzy&
&aaaj113),xyzzyaaap113(1),nitot,0,1,eivecs(1,1,xyzzyaaaa113),xyzzyaaaf&
&113,nfn_enbasis,enbasis(xyzzyaaah113,1,xyzzyaaaa113),grad_enbasis(1,x&
&yzzyaaah113,1,xyzzyaaaa113),denbasis(xyzzyaaah113,1,xyzzyaaaa113),gra&
&dr_enbasis(1,1,xyzzyaaaa113),nzencut(1,xyzzyaaaa113,xyzzyaaad113),xyz&
&zyaaao113)
enddo
endif
enddo
end subroutine xyzzyaagh1
subroutine xyzzyaagi1(eevecs,eebasis,grad_eebasis,deebasis,gradr_eebas&
&is,lap_eebasis,d2eebasis,lapr_eebasis,nzeecut)
implicit none
real(dp),intent(in) :: eevecs(4,netot,netot)
real(dp),intent(inout) :: eebasis(nfn_eebasis,netot,netot),grad_eebasi&
&s(3,nfn_eebasis,netot,netot),deebasis(nfn_eebasis,netot,netot),gradr_&
&eebasis(3,netot,netot),lap_eebasis(nfn_eebasis,netot,netot),lapr_eeba&
&sis(netot,netot),d2eebasis(nfn_eebasis,netot,netot)
logical,intent(inout) :: nzeecut(netot,netot,xyzzyaaaw1)
integer xyzzyaaaa114,xyzzyaaab114,xyzzyaaac114,xyzzyaaad114,xyzzyaaae1&
&14,xyzzyaaaf114,xyzzyaaag114,xyzzyaaah114,xyzzyaaai114,xyzzyaaaj114,x&
&yzzyaaak114,xyzzyaaal114,xyzzyaaam114,xyzzyaaan114,xyzzyaaao114,xyzzy&
&aaap114,xyzzyaaaq114
logical xyzzyaaar114
real(dp) xyzzyaaas114(1)
do xyzzyaaaa114=1,netot
do xyzzyaaal114=xyzzyaaaa114+1,netot
call xyzzyaaiy1(eevecs(4,xyzzyaaal114,xyzzyaaaa114),eevecs(1,xyzzyaaal&
&114,xyzzyaaaa114),gradr_eebasis(1,xyzzyaaal114,xyzzyaaaa114),lapr_eeb&
&asis(xyzzyaaal114,xyzzyaaaa114))
gradr_eebasis(1:3,xyzzyaaaa114,xyzzyaaal114)=-gradr_eebasis(1:3,xyzzya&
&aal114,xyzzyaaaa114)
lapr_eebasis(xyzzyaaaa114,xyzzyaaal114)=lapr_eebasis(xyzzyaaal114,xyzz&
&yaaaa114)
enddo
enddo
do xyzzyaaad114=1,xyzzyaaaw1
xyzzyaaar114=xyzzyaaca1(xyzzyaaad114).and.xyzzyaacc1(xyzzyaaad114)
xyzzyaaah114=ifn1_eebasis(xyzzyaaad114)
xyzzyaaae114=xyzzyaack1(xyzzyaaad114)
xyzzyaaaf114=xyzzyaabq1(xyzzyaaad114)
xyzzyaaag114=xyzzyaabs1(xyzzyaaad114)
xyzzyaaaj114=xyzzyaabk1(xyzzyaaad114)+1
xyzzyaaak114=xyzzyaabi1(xyzzyaaad114)+1
if(xyzzyaaae114>0)then
xyzzyaaai114=xyzzyaaay1(xyzzyaaad114)
xyzzyaaap114=0
do xyzzyaaac114=1,nspin
if(nele(xyzzyaaac114)==0)cycle
if(nele(xyzzyaaac114)>1)then
xyzzyaaao114=xyzzyaabu1(xyzzyaaac114,xyzzyaaac114,xyzzyaaad114)
if(xyzzyaaao114>0)then
xyzzyaaan114=xyzzyaabg1(xyzzyaaai114+xyzzyaaao114)+1
xyzzyaaaa114=xyzzyaaap114
xyzzyaaal114=xyzzyaaap114+1
xyzzyaaaq114=nele(xyzzyaaac114)
do xyzzyaaab114=1,nele(xyzzyaaac114)-1
xyzzyaaaa114=xyzzyaaaa114+1
xyzzyaaal114=xyzzyaaal114+1
xyzzyaaaq114=xyzzyaaaq114-1
call xyzzyaagm1(xyzzyaaag114,xyzzyaads1(xyzzyaaak114),xyzzyaadu1(xyzzy&
&aaaj114),xyzzyaadi1(xyzzyaaan114),xyzzyaaaq114,0,0,eevecs(1,xyzzyaaal&
&114,xyzzyaaaa114),xyzzyaaaf114,nfn_eebasis,eebasis(xyzzyaaah114,xyzzy&
&aaal114,xyzzyaaaa114),grad_eebasis(1,xyzzyaaah114,xyzzyaaal114,xyzzya&
&aaa114),deebasis(xyzzyaaah114,xyzzyaaal114,xyzzyaaaa114),gradr_eebasi&
&s(1,xyzzyaaal114,xyzzyaaaa114),lap_eebasis(xyzzyaaah114,xyzzyaaal114,&
&xyzzyaaaa114),d2eebasis(xyzzyaaah114,xyzzyaaal114,xyzzyaaaa114),lapr_&
&eebasis(xyzzyaaal114,xyzzyaaaa114),nzeecut(xyzzyaaal114,xyzzyaaaa114,&
&xyzzyaaad114),xyzzyaaar114)
enddo
endif
endif
xyzzyaaal114=xyzzyaaap114+nele(xyzzyaaac114)+1
do xyzzyaaam114=xyzzyaaac114+1,nspin
if(nele(xyzzyaaam114)==0)cycle
xyzzyaaao114=xyzzyaabu1(xyzzyaaac114,xyzzyaaam114,xyzzyaaad114)
if(xyzzyaaao114>0)then
xyzzyaaan114=xyzzyaabg1(xyzzyaaai114+xyzzyaaao114)+1
xyzzyaaaa114=xyzzyaaap114
xyzzyaaaq114=nele(xyzzyaaam114)
do xyzzyaaab114=1,nele(xyzzyaaac114)
xyzzyaaaa114=xyzzyaaaa114+1
call xyzzyaagm1(xyzzyaaag114,xyzzyaads1(xyzzyaaak114),xyzzyaadu1(xyzzy&
&aaaj114),xyzzyaadi1(xyzzyaaan114),xyzzyaaaq114,0,0,eevecs(1,xyzzyaaal&
&114,xyzzyaaaa114),xyzzyaaaf114,nfn_eebasis,eebasis(xyzzyaaah114,xyzzy&
&aaal114,xyzzyaaaa114),grad_eebasis(1,xyzzyaaah114,xyzzyaaal114,xyzzya&
&aaa114),deebasis(xyzzyaaah114,xyzzyaaal114,xyzzyaaaa114),gradr_eebasi&
&s(1,xyzzyaaal114,xyzzyaaaa114),lap_eebasis(xyzzyaaah114,xyzzyaaal114,&
&xyzzyaaaa114),d2eebasis(xyzzyaaah114,xyzzyaaal114,xyzzyaaaa114),lapr_&
&eebasis(xyzzyaaal114,xyzzyaaaa114),nzeecut(xyzzyaaal114,xyzzyaaaa114,&
&xyzzyaaad114),xyzzyaaar114)
enddo
endif
xyzzyaaal114=xyzzyaaal114+nele(xyzzyaaam114)
enddo
xyzzyaaap114=xyzzyaaap114+nele(xyzzyaaac114)
enddo
else
xyzzyaaal114=1
xyzzyaaaq114=netot
do xyzzyaaaa114=1,netot-1
xyzzyaaal114=xyzzyaaal114+1
xyzzyaaaq114=xyzzyaaaq114-1
call xyzzyaagm1(xyzzyaaag114,xyzzyaads1(xyzzyaaak114),xyzzyaadu1(xyzzy&
&aaaj114),xyzzyaaas114(1),xyzzyaaaq114,0,0,eevecs(1,xyzzyaaal114,xyzzy&
&aaaa114),xyzzyaaaf114,nfn_eebasis,eebasis(xyzzyaaah114,xyzzyaaal114,x&
&yzzyaaaa114),grad_eebasis(1,xyzzyaaah114,xyzzyaaal114,xyzzyaaaa114),d&
&eebasis(xyzzyaaah114,xyzzyaaal114,xyzzyaaaa114),gradr_eebasis(1,xyzzy&
&aaal114,xyzzyaaaa114),lap_eebasis(xyzzyaaah114,xyzzyaaal114,xyzzyaaaa&
&114),d2eebasis(xyzzyaaah114,xyzzyaaal114,xyzzyaaaa114),lapr_eebasis(x&
&yzzyaaal114,xyzzyaaaa114),nzeecut(xyzzyaaal114,xyzzyaaaa114,xyzzyaaad&
&114),xyzzyaaar114)
enddo
endif
enddo
do xyzzyaaaa114=1,netot-1
eebasis(:,xyzzyaaaa114,xyzzyaaaa114+1:netot)=eebasis(:,xyzzyaaaa114+1:&
&netot,xyzzyaaaa114)
if(xyzzyaabm1)then
grad_eebasis(:,:,xyzzyaaaa114,xyzzyaaaa114+1:netot)=-grad_eebasis(:,:,&
&xyzzyaaaa114+1:netot,xyzzyaaaa114)
lap_eebasis(:,xyzzyaaaa114,xyzzyaaaa114+1:netot)=lap_eebasis(:,xyzzyaa&
&aa114+1:netot,xyzzyaaaa114)
endif
if(xyzzyaabn1)then
deebasis(:,xyzzyaaaa114,xyzzyaaaa114+1:netot)=deebasis(:,xyzzyaaaa114+&
&1:netot,xyzzyaaaa114)
d2eebasis(:,xyzzyaaaa114,xyzzyaaaa114+1:netot)=d2eebasis(:,xyzzyaaaa11&
&4+1:netot,xyzzyaaaa114)
endif
nzeecut(xyzzyaaaa114,xyzzyaaaa114+1:netot,:)=nzeecut(xyzzyaaaa114+1:ne&
&tot,xyzzyaaaa114,:)
enddo
end subroutine xyzzyaagi1
subroutine xyzzyaagj1(eivecs,enbasis,grad_enbasis,denbasis,gradr_enbas&
&is,lap_enbasis,d2enbasis,lapr_enbasis,nzencut)
implicit none
real(dp),intent(in) :: eivecs(4,nitot,netot)
real(dp),intent(inout) :: enbasis(nfn_enbasis,nitot,netot),grad_enbasi&
&s(3,nfn_enbasis,nitot,netot),denbasis(nfn_enbasis,nitot,netot),gradr_&
&enbasis(3,nitot,netot),lap_enbasis(nfn_enbasis,nitot,netot),d2enbasis&
&(nfn_enbasis,nitot,netot),lapr_enbasis(nitot,netot)
logical,intent(inout) :: nzencut(nitot,netot,xyzzyaaax1)
integer xyzzyaaaa115,xyzzyaaab115,xyzzyaaac115,xyzzyaaad115,xyzzyaaae1&
&15,xyzzyaaaf115,xyzzyaaag115,xyzzyaaah115,xyzzyaaai115,xyzzyaaaj115,x&
&yzzyaaak115,xyzzyaaal115,xyzzyaaam115,xyzzyaaan115
logical xyzzyaaao115
real(dp) xyzzyaaap115(1)
do xyzzyaaaa115=1,netot
do xyzzyaaan115=1,nitot
call xyzzyaaiy1(eivecs(4,xyzzyaaan115,xyzzyaaaa115),eivecs(1,xyzzyaaan&
&115,xyzzyaaaa115),gradr_enbasis(1,xyzzyaaan115,xyzzyaaaa115),lapr_enb&
&asis(xyzzyaaan115,xyzzyaaaa115))
enddo
enddo
do xyzzyaaad115=1,xyzzyaaax1
xyzzyaaao115=xyzzyaacb1(xyzzyaaad115).and.xyzzyaacd1(xyzzyaaad115)
xyzzyaaah115=ifn1_enbasis(xyzzyaaad115)
xyzzyaaae115=xyzzyaacl1(xyzzyaaad115)
xyzzyaaaf115=xyzzyaabr1(xyzzyaaad115)
xyzzyaaag115=xyzzyaabt1(xyzzyaaad115)
xyzzyaaaj115=xyzzyaabl1(xyzzyaaad115)+1
xyzzyaaak115=xyzzyaabj1(xyzzyaaad115)+1
if(xyzzyaaae115>0)then
xyzzyaaai115=xyzzyaaaz1(xyzzyaaad115)
do xyzzyaaan115=1,nitot
xyzzyaaaa115=0
do xyzzyaaac115=1,nspin
if(nele(xyzzyaaac115)==0)cycle
xyzzyaaam115=xyzzyaabv1(xyzzyaaan115,xyzzyaaac115,xyzzyaaad115)
if(xyzzyaaam115<1)then
xyzzyaaaa115=xyzzyaaaa115+nele(xyzzyaaac115)
cycle
endif
xyzzyaaal115=xyzzyaabh1(xyzzyaaai115+xyzzyaaam115)+1
do xyzzyaaab115=1,nele(xyzzyaaac115)
xyzzyaaaa115=xyzzyaaaa115+1
call xyzzyaagm1(xyzzyaaag115,xyzzyaadt1(xyzzyaaak115),xyzzyaadv1(xyzzy&
&aaaj115),xyzzyaadj1(xyzzyaaal115),1,0,xyzzyaaan115,eivecs(1,xyzzyaaan&
&115,xyzzyaaaa115),xyzzyaaaf115,nfn_enbasis,enbasis(xyzzyaaah115,xyzzy&
&aaan115,xyzzyaaaa115),grad_enbasis(1,xyzzyaaah115,xyzzyaaan115,xyzzya&
&aaa115),denbasis(xyzzyaaah115,xyzzyaaan115,xyzzyaaaa115),gradr_enbasi&
&s(1,xyzzyaaan115,xyzzyaaaa115),lap_enbasis(xyzzyaaah115,xyzzyaaan115,&
&xyzzyaaaa115),d2enbasis(xyzzyaaah115,xyzzyaaan115,xyzzyaaaa115),lapr_&
&enbasis(xyzzyaaan115,xyzzyaaaa115),nzencut(xyzzyaaan115,xyzzyaaaa115,&
&xyzzyaaad115),xyzzyaaao115)
enddo
enddo
enddo
else
do xyzzyaaaa115=1,netot
call xyzzyaagm1(xyzzyaaag115,xyzzyaadt1(xyzzyaaak115),xyzzyaadv1(xyzzy&
&aaaj115),xyzzyaaap115(1),nitot,0,1,eivecs(1,1,xyzzyaaaa115),xyzzyaaaf&
&115,nfn_enbasis,enbasis(xyzzyaaah115,1,xyzzyaaaa115),grad_enbasis(1,x&
&yzzyaaah115,1,xyzzyaaaa115),denbasis(xyzzyaaah115,1,xyzzyaaaa115),gra&
&dr_enbasis(1,1,xyzzyaaaa115),lap_enbasis(xyzzyaaah115,1,xyzzyaaaa115)&
&,d2enbasis(xyzzyaaah115,1,xyzzyaaaa115),lapr_enbasis(1,xyzzyaaaa115),&
&nzencut(1,xyzzyaaaa115,xyzzyaaad115),xyzzyaaao115)
enddo
endif
enddo
end subroutine xyzzyaagj1
subroutine xyzzyaagk1(n,ci,cd,p,ne,ieskip,ion1,exvecs,ibasis,fstride,f&
&,nzcut)
implicit none
integer,intent(in) :: n,ne,ieskip,ion1,ibasis,fstride,ci(*)
real(dp),intent(in) :: exvecs(4,ne),cd(*),p(*)
real(dp),intent(inout) :: f(fstride,*)
logical,intent(inout) :: nzcut(ne)
select case(ibasis)
case(0)
f(1:n,1:ne)=1.d0
nzcut(1:ne)=.true.
case(xyzzyaaaa1)
call xyzzyaagn1(n,ne,ieskip,exvecs,fstride,f)
case(xyzzyaaab1)
call xyzzyaagq1(n,ci(1),cd(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaac1)
call xyzzyaagt1(n,ci(1),cd(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaad1)
call xyzzyaagw1(n,cd(1),p(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaae1)
call xyzzyaagz1(n,p(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaaf1)
call xyzzyaahc1(n,p(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaag1)
call xyzzyaahf1(n,ne,ieskip,exvecs,fstride,f)
case(xyzzyaaah1)
call xyzzyaahi1(n,ci(1),ci(2),p(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaai1)
call xyzzyaahl1(n,ci(1),ci(2),p(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaaj1)
call xyzzyaaho1(p(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaak1)
call xyzzyaahr1(ne,ieskip,exvecs,fstride,f)
case(xyzzyaaal1)
call xyzzyaahu1(ne,ieskip,exvecs,fstride,f)
case(xyzzyaaam1)
call xyzzyaahx1(n,ci(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaan1)
call xyzzyaaia1(p(1),ne,ieskip,exvecs,fstride,f)
case(xyzzyaaao1)
call xyzzyaaid1(n,ci(1),ci(2),cd(1),cd(2),cd(ci(3)),cd(ci(4)),ne,ieski&
&p,exvecs,fstride,f)
case(xyzzyaaap1)
call xyzzyaaig1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,nzcut)
case(xyzzyaaaq1)
call xyzzyaaij1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,nzcut)
case(xyzzyaaar1)
call xyzzyaaim1(p(1),ne,ieskip,exvecs,fstride,f,nzcut)
case(xyzzyaaas1)
call xyzzyaaip1(ci(1),cd(1),p(1),ne,ieskip,ion1,exvecs,fstride,f,nzcut&
&)
case(xyzzyaaat1)
call xyzzyaais1(p(1),ne,ieskip,exvecs,fstride,f,nzcut)
case(xyzzyaaau1)
call xyzzyaaiv1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,nzcut)
end select
end subroutine xyzzyaagk1
subroutine xyzzyaagl1(n,ci,cd,p,ne,ieskip,ion1,exvecs,ibasis,fstride,f&
&,gradf,df,gradr,nzcut,iso_want_aniso)
implicit none
integer,intent(in) :: n,ne,ieskip,ion1,ibasis,fstride,ci(*)
logical,intent(in) :: iso_want_aniso
real(dp),intent(in) :: exvecs(4,ne),cd(*),p(*),gradr(3,*)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),df(fstride,*&
&)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa117,xyzzyaaab117
select case(ibasis)
case(0)
f(1:n,1:ne)=1.d0
df(1:n,1:ne)=0.d0
nzcut(1:ne)=.true.
case(xyzzyaaaa1)
call xyzzyaago1(n,ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaab1)
call xyzzyaagr1(n,ci(1),cd(1),ne,ieskip,exvecs,fstride,f,gradf)
case(xyzzyaaac1)
call xyzzyaagu1(n,ci(1),cd(1),ne,ieskip,exvecs,fstride,f,gradf)
case(xyzzyaaad1)
call xyzzyaagx1(n,cd(1),p(1),ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaae1)
call xyzzyaaha1(n,p(1),ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaaf1)
call xyzzyaahd1(n,p(1),ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaag1)
call xyzzyaahg1(n,ne,ieskip,exvecs,fstride,f,gradf)
case(xyzzyaaah1)
call xyzzyaahj1(n,ci(1),ci(2),p(1),ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaai1)
call xyzzyaahm1(n,ci(1),ci(2),p(1),ne,ieskip,exvecs,fstride,f,gradf)
case(xyzzyaaaj1)
call xyzzyaahp1(p(1),ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaak1)
call xyzzyaahs1(ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaal1)
call xyzzyaahv1(ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaam1)
call xyzzyaahy1(n,ci(1),ne,ieskip,exvecs,fstride,f,df)
case(xyzzyaaan1)
call xyzzyaaib1(p(1),ne,ieskip,exvecs,fstride,f,gradf)
case(xyzzyaaao1)
call xyzzyaaie1(n,ci(1),ci(2),cd(1),cd(2),cd(ci(3)),cd(ci(4)),ne,ieski&
&p,exvecs,fstride,f,gradf)
case(xyzzyaaap1)
call xyzzyaaih1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,df,nzcut)
case(xyzzyaaaq1)
call xyzzyaaik1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,df,nzcut)
case(xyzzyaaar1)
call xyzzyaain1(p(1),ne,ieskip,exvecs,fstride,f,df,nzcut)
case(xyzzyaaas1)
call xyzzyaaiq1(ci(1),cd(1),p(1),ne,ieskip,ion1,exvecs,fstride,f,gradf&
&,nzcut)
case(xyzzyaaat1)
call xyzzyaait1(p(1),ne,ieskip,exvecs,fstride,f,df,nzcut)
case(xyzzyaaau1)
call xyzzyaaiw1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,df,nzcut)
end select
if(iso_want_aniso)then
if(ibasis==0)then
gradf(1:3,1:n,1:ne)=0.d0
else
do xyzzyaaab117=1,ne
do xyzzyaaaa117=1,n
gradf(1:3,xyzzyaaaa117,xyzzyaaab117)=df(xyzzyaaaa117,xyzzyaaab117)*gra&
&dr(1:3,xyzzyaaab117)
enddo
enddo
endif
endif
end subroutine xyzzyaagl1
subroutine xyzzyaagm1(n,ci,cd,p,ne,ieskip,ion1,exvecs,ibasis,fstride,f&
&,gradf,df,gradr,lapf,d2f,lapr,nzcut,iso_want_aniso)
implicit none
integer,intent(in) :: n,ne,ieskip,ion1,ibasis,fstride,ci(*)
logical,intent(in) :: iso_want_aniso
real(dp),intent(in) :: exvecs(4,ne),cd(*),p(*),gradr(3,*),lapr(*)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),df(fstride,*&
&),lapf(fstride,*),d2f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa118,xyzzyaaab118
select case(ibasis)
case(0)
f(1:n,1:ne)=1.d0
df(1:n,1:ne)=0.d0
d2f(1:n,1:ne)=0.d0
nzcut(1:ne)=.true.
case(xyzzyaaaa1)
call xyzzyaagp1(n,ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaab1)
call xyzzyaags1(n,ci(1),cd(1),ne,ieskip,exvecs,fstride,f,gradf,lapf)
case(xyzzyaaac1)
call xyzzyaagv1(n,ci(1),cd(1),ne,ieskip,exvecs,fstride,f,gradf,lapf)
case(xyzzyaaad1)
call xyzzyaagy1(n,cd(1),p(1),ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaae1)
call xyzzyaahb1(n,p(1),ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaaf1)
call xyzzyaahe1(n,p(1),ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaag1)
call xyzzyaahh1(n,ne,ieskip,exvecs,fstride,f,gradf,lapf)
case(xyzzyaaah1)
call xyzzyaahk1(n,ci(1),ci(2),p(1),ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaai1)
call xyzzyaahn1(n,ci(1),ci(2),p(1),ne,ieskip,exvecs,fstride,f,gradf,la&
&pf)
case(xyzzyaaaj1)
call xyzzyaahq1(p(1),ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaak1)
call xyzzyaaht1(ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaal1)
call xyzzyaahw1(ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaam1)
call xyzzyaahz1(n,ci(1),ne,ieskip,exvecs,fstride,f,df,d2f)
case(xyzzyaaan1)
call xyzzyaaic1(p(1),ne,ieskip,exvecs,fstride,f,gradf,lapf)
case(xyzzyaaao1)
call xyzzyaaif1(n,ci(1),ci(2),cd(1),cd(2),cd(ci(3)),cd(ci(4)),cd(ci(5)&
&),ne,ieskip,exvecs,fstride,f,gradf,lapf)
case(xyzzyaaap1)
call xyzzyaaii1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,df,d2f,nzcut)
case(xyzzyaaaq1)
call xyzzyaail1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,df,d2f,nzcut)
case(xyzzyaaar1)
call xyzzyaaio1(p(1),ne,ieskip,exvecs,fstride,f,df,d2f,nzcut)
case(xyzzyaaas1)
call xyzzyaair1(ci(1),cd(1),p(1),ne,ieskip,ion1,exvecs,fstride,f,gradf&
&,lapf,nzcut)
case(xyzzyaaat1)
call xyzzyaaiu1(p(1),ne,ieskip,exvecs,fstride,f,df,d2f,nzcut)
case(xyzzyaaau1)
call xyzzyaaix1(ci(1),p(1),ne,ieskip,exvecs,fstride,f,df,d2f,nzcut)
end select
if(iso_want_aniso)then
if(ibasis==0)then
gradf(1:3,1:n,1:ne)=0.d0
lapf(1:n,1:ne)=0.d0
else
do xyzzyaaab118=1,ne
do xyzzyaaaa118=1,n
gradf(1:3,xyzzyaaaa118,xyzzyaaab118)=df(xyzzyaaaa118,xyzzyaaab118)*gra&
&dr(1:3,xyzzyaaab118)
lapf(xyzzyaaaa118,xyzzyaaab118)=d2f(xyzzyaaaa118,xyzzyaaab118)+df(xyzz&
&yaaaa118,xyzzyaaab118)*lapr(xyzzyaaab118)
enddo
enddo
endif
endif
end subroutine xyzzyaagm1
subroutine xyzzyaagn1(n,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa119,xyzzyaaab119
real(dp) xyzzyaaac119,xyzzyaaad119
do xyzzyaaaa119=1,ne
if(xyzzyaaaa119==ieskip)cycle
xyzzyaaac119=exvecs(4,xyzzyaaaa119)
f(1,xyzzyaaaa119)=1.d0
f(2,xyzzyaaaa119)=xyzzyaaac119
xyzzyaaad119=xyzzyaaac119
do xyzzyaaab119=3,n
xyzzyaaad119=xyzzyaaad119*xyzzyaaac119
f(xyzzyaaab119,xyzzyaaaa119)=xyzzyaaad119
enddo
enddo
end subroutine xyzzyaagn1
subroutine xyzzyaago1(n,ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa120,xyzzyaaab120
real(dp) xyzzyaaac120,xyzzyaaad120,xyzzyaaae120
do xyzzyaaaa120=1,ne
if(xyzzyaaaa120==ieskip)cycle
xyzzyaaac120=exvecs(4,xyzzyaaaa120)
f(1,xyzzyaaaa120)=1.d0
df(1,xyzzyaaaa120)=0.d0
f(2,xyzzyaaaa120)=xyzzyaaac120
df(2,xyzzyaaaa120)=1.d0
xyzzyaaad120=xyzzyaaac120
xyzzyaaae120=1.d0
do xyzzyaaab120=3,n
xyzzyaaae120=xyzzyaaad120*(xyzzyaaab120-1)
xyzzyaaad120=xyzzyaaad120*xyzzyaaac120
f(xyzzyaaab120,xyzzyaaaa120)=xyzzyaaad120
df(xyzzyaaab120,xyzzyaaaa120)=xyzzyaaae120
enddo
enddo
end subroutine xyzzyaago1
subroutine xyzzyaagp1(n,ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa121,xyzzyaaab121
real(dp) xyzzyaaac121,xyzzyaaad121,xyzzyaaae121,xyzzyaaaf121
do xyzzyaaaa121=1,ne
if(xyzzyaaaa121==ieskip)cycle
xyzzyaaac121=exvecs(4,xyzzyaaaa121)
f(1,xyzzyaaaa121)=1.d0
df(1,xyzzyaaaa121)=0.d0
d2f(1,xyzzyaaaa121)=0.d0
f(2,xyzzyaaaa121)=xyzzyaaac121
df(2,xyzzyaaaa121)=1.d0
d2f(2,xyzzyaaaa121)=0.d0
xyzzyaaad121=xyzzyaaac121
xyzzyaaae121=1.d0
xyzzyaaaf121=0.d0
do xyzzyaaab121=3,n
xyzzyaaaf121=xyzzyaaae121*(xyzzyaaab121-1)
xyzzyaaae121=xyzzyaaad121*(xyzzyaaab121-1)
xyzzyaaad121=xyzzyaaad121*xyzzyaaac121
f(xyzzyaaab121,xyzzyaaaa121)=xyzzyaaad121
df(xyzzyaaab121,xyzzyaaaa121)=xyzzyaaae121
d2f(xyzzyaaab121,xyzzyaaaa121)=xyzzyaaaf121
enddo
enddo
end subroutine xyzzyaagp1
subroutine xyzzyaagq1(n,nk,k,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,nk(n)
real(dp),intent(in) :: k(periodicity+1,*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa122,xyzzyaaab122,xyzzyaaac122,xyzzyaaad122
real(dp) xyzzyaaae122(3),xyzzyaaaf122,xyzzyaaag122
do xyzzyaaaa122=1,ne
if(xyzzyaaaa122==ieskip)cycle
xyzzyaaae122=exvecs(1:3,xyzzyaaaa122)
xyzzyaaad122=0
do xyzzyaaab122=1,n
xyzzyaaaf122=0.d0
do xyzzyaaac122=1,nk(xyzzyaaab122)
xyzzyaaad122=xyzzyaaad122+1
xyzzyaaag122=ddot(periodicity,k(1,xyzzyaaad122),1,xyzzyaaae122(1),1)
xyzzyaaaf122=xyzzyaaaf122+cos(xyzzyaaag122)
enddo
f(xyzzyaaab122,xyzzyaaaa122)=xyzzyaaaf122
enddo
enddo
end subroutine xyzzyaagq1
subroutine xyzzyaagr1(n,nk,k,ne,ieskip,exvecs,fstride,f,gradf)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,nk(n)
real(dp),intent(in) :: k(periodicity+1,*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*)
integer xyzzyaaaa123,xyzzyaaab123,xyzzyaaac123,xyzzyaaad123
real(dp) xyzzyaaae123(3),xyzzyaaaf123,xyzzyaaag123(3),xyzzyaaah123
do xyzzyaaaa123=1,ne
if(xyzzyaaaa123==ieskip)cycle
xyzzyaaae123=exvecs(1:3,xyzzyaaaa123)
xyzzyaaad123=0
do xyzzyaaab123=1,n
xyzzyaaaf123=0.d0
xyzzyaaag123=0.d0
do xyzzyaaac123=1,nk(xyzzyaaab123)
xyzzyaaad123=xyzzyaaad123+1
xyzzyaaah123=ddot(periodicity,k(1,xyzzyaaad123),1,xyzzyaaae123(1),1)
xyzzyaaaf123=xyzzyaaaf123+cos(xyzzyaaah123)
xyzzyaaag123(1:periodicity)=xyzzyaaag123(1:periodicity)-k(1:periodicit&
&y,xyzzyaaad123)*sin(xyzzyaaah123)
enddo
f(xyzzyaaab123,xyzzyaaaa123)=xyzzyaaaf123
gradf(1:3,xyzzyaaab123,xyzzyaaaa123)=xyzzyaaag123(1:3)
enddo
enddo
end subroutine xyzzyaagr1
subroutine xyzzyaags1(n,nk,k,ne,ieskip,exvecs,fstride,f,gradf,lapf)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,nk(n)
real(dp),intent(in) :: k(periodicity+1,*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),lapf(fstride&
&,*)
integer xyzzyaaaa124,xyzzyaaab124,xyzzyaaac124,xyzzyaaad124
real(dp) xyzzyaaae124(3),xyzzyaaaf124,xyzzyaaag124(3),xyzzyaaah124,xyz&
&zyaaai124,xyzzyaaaj124
do xyzzyaaaa124=1,ne
if(xyzzyaaaa124==ieskip)cycle
xyzzyaaae124=exvecs(1:3,xyzzyaaaa124)
xyzzyaaad124=0
do xyzzyaaab124=1,n
xyzzyaaaf124=0.d0
xyzzyaaag124=0.d0
xyzzyaaah124=0.d0
do xyzzyaaac124=1,nk(xyzzyaaab124)
xyzzyaaad124=xyzzyaaad124+1
xyzzyaaai124=ddot(periodicity,k(1,xyzzyaaad124),1,xyzzyaaae124(1),1)
xyzzyaaaj124=cos(xyzzyaaai124)
xyzzyaaaf124=xyzzyaaaf124+xyzzyaaaj124
xyzzyaaag124(1:periodicity)=xyzzyaaag124(1:periodicity)-k(1:periodicit&
&y,xyzzyaaad124)*sin(xyzzyaaai124)
xyzzyaaah124=xyzzyaaah124-k(periodicity+1,xyzzyaaad124)*xyzzyaaaj124
enddo
f(xyzzyaaab124,xyzzyaaaa124)=xyzzyaaaf124
gradf(1:3,xyzzyaaab124,xyzzyaaaa124)=xyzzyaaag124(1:3)
lapf(xyzzyaaab124,xyzzyaaaa124)=xyzzyaaah124
enddo
enddo
end subroutine xyzzyaags1
subroutine xyzzyaagt1(n,const_int,const_dble,ne,ieskip,exvecs,fstride,&
&f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,const_int(*)
real(dp),intent(in) :: const_dble(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa125,xyzzyaaab125,xyzzyaaac125,xyzzyaaad125,xyzzyaaae1&
&25,xyzzyaaaf125,xyzzyaaag125,xyzzyaaah125,xyzzyaaai125,xyzzyaaaj125
real(dp) xyzzyaaak125(3),xyzzyaaal125,xyzzyaaam125
xyzzyaaae125=const_int(1)
xyzzyaaaj125=xyzzyaaae125*n
do xyzzyaaaa125=1,ne
if(xyzzyaaaa125==ieskip)cycle
xyzzyaaak125(1:3)=exvecs(1:3,xyzzyaaaa125)
f(1:n,xyzzyaaaa125)=0.d0
xyzzyaaag125=1
xyzzyaaah125=xyzzyaaaj125
xyzzyaaai125=0
do xyzzyaaad125=1,xyzzyaaae125
xyzzyaaag125=xyzzyaaag125+1
xyzzyaaaf125=const_int(xyzzyaaag125)
xyzzyaaal125=0.d0
do xyzzyaaac125=1,xyzzyaaaf125
xyzzyaaam125=ddot(periodicity,const_dble(xyzzyaaah125+1),1,xyzzyaaak12&
&5(1),1)
xyzzyaaal125=xyzzyaaal125+cos(xyzzyaaam125)
xyzzyaaah125=xyzzyaaah125+periodicity+1
enddo
do xyzzyaaab125=1,n
xyzzyaaai125=xyzzyaaai125+1
f(xyzzyaaab125,xyzzyaaaa125)=f(xyzzyaaab125,xyzzyaaaa125)+const_dble(x&
&yzzyaaai125)*xyzzyaaal125
enddo
enddo
enddo
end subroutine xyzzyaagt1
subroutine xyzzyaagu1(n,const_int,const_dble,ne,ieskip,exvecs,fstride,&
&f,gradf)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,const_int(*)
real(dp),intent(in) :: const_dble(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*)
integer xyzzyaaaa126,xyzzyaaab126,xyzzyaaac126,xyzzyaaad126,xyzzyaaae1&
&26,xyzzyaaaf126,xyzzyaaag126,xyzzyaaah126,xyzzyaaai126,xyzzyaaaj126
real(dp) xyzzyaaak126(3),xyzzyaaal126,xyzzyaaam126(3),xyzzyaaan126
xyzzyaaae126=const_int(1)
xyzzyaaaj126=xyzzyaaae126*n
do xyzzyaaaa126=1,ne
if(xyzzyaaaa126==ieskip)cycle
xyzzyaaak126(1:3)=exvecs(1:3,xyzzyaaaa126)
f(1:n,xyzzyaaaa126)=0.d0
gradf(1:3,1:n,xyzzyaaaa126)=0.d0
xyzzyaaag126=1
xyzzyaaah126=xyzzyaaaj126
xyzzyaaai126=0
do xyzzyaaad126=1,xyzzyaaae126
xyzzyaaag126=xyzzyaaag126+1
xyzzyaaaf126=const_int(xyzzyaaag126)
xyzzyaaal126=0.d0
xyzzyaaam126=0.d0
do xyzzyaaac126=1,xyzzyaaaf126
xyzzyaaan126=ddot(periodicity,const_dble(xyzzyaaah126+1),1,xyzzyaaak12&
&6(1),1)
xyzzyaaal126=xyzzyaaal126+cos(xyzzyaaan126)
xyzzyaaam126(1:periodicity)=xyzzyaaam126(1:periodicity)-const_dble(xyz&
&zyaaah126+1:xyzzyaaah126+periodicity)*sin(xyzzyaaan126)
xyzzyaaah126=xyzzyaaah126+periodicity+1
enddo
do xyzzyaaab126=1,n
xyzzyaaai126=xyzzyaaai126+1
f(xyzzyaaab126,xyzzyaaaa126)=f(xyzzyaaab126,xyzzyaaaa126)+const_dble(x&
&yzzyaaai126)*xyzzyaaal126
gradf(1:periodicity,xyzzyaaab126,xyzzyaaaa126)=gradf(1:periodicity,xyz&
&zyaaab126,xyzzyaaaa126)+const_dble(xyzzyaaai126)*xyzzyaaam126(1:perio&
&dicity)
enddo
enddo
enddo
end subroutine xyzzyaagu1
subroutine xyzzyaagv1(n,const_int,const_dble,ne,ieskip,exvecs,fstride,&
&f,gradf,lapf)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,const_int(*)
real(dp),intent(in) :: const_dble(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),lapf(fstride&
&,*)
integer xyzzyaaaa127,xyzzyaaab127,xyzzyaaac127,xyzzyaaad127,xyzzyaaae1&
&27,xyzzyaaaf127,xyzzyaaag127,xyzzyaaah127,xyzzyaaai127,xyzzyaaaj127
real(dp) xyzzyaaak127(3),xyzzyaaal127,xyzzyaaam127(3),xyzzyaaan127,xyz&
&zyaaao127,xyzzyaaap127
xyzzyaaae127=const_int(1)
xyzzyaaaj127=xyzzyaaae127*n
do xyzzyaaaa127=1,ne
if(xyzzyaaaa127==ieskip)cycle
xyzzyaaak127(1:3)=exvecs(1:3,xyzzyaaaa127)
f(1:n,xyzzyaaaa127)=0.d0
gradf(1:3,1:n,xyzzyaaaa127)=0.d0
lapf(1:n,xyzzyaaaa127)=0.d0
xyzzyaaag127=1
xyzzyaaah127=xyzzyaaaj127
xyzzyaaai127=0
do xyzzyaaad127=1,xyzzyaaae127
xyzzyaaag127=xyzzyaaag127+1
xyzzyaaaf127=const_int(xyzzyaaag127)
xyzzyaaal127=0.d0
xyzzyaaam127=0.d0
do xyzzyaaac127=1,xyzzyaaaf127
xyzzyaaan127=ddot(periodicity,const_dble(xyzzyaaah127+1),1,xyzzyaaak12&
&7(1),1)
xyzzyaaal127=xyzzyaaal127+cos(xyzzyaaan127)
xyzzyaaam127(1:periodicity)=xyzzyaaam127(1:periodicity)-const_dble(xyz&
&zyaaah127+1:xyzzyaaah127+periodicity)*sin(xyzzyaaan127)
xyzzyaaah127=xyzzyaaah127+periodicity+1
enddo
xyzzyaaao127=const_dble(xyzzyaaah127)
do xyzzyaaab127=1,n
xyzzyaaai127=xyzzyaaai127+1
xyzzyaaap127=const_dble(xyzzyaaai127)*xyzzyaaal127
f(xyzzyaaab127,xyzzyaaaa127)=f(xyzzyaaab127,xyzzyaaaa127)+xyzzyaaap127
gradf(1:periodicity,xyzzyaaab127,xyzzyaaaa127)=gradf(1:periodicity,xyz&
&zyaaab127,xyzzyaaaa127)+const_dble(xyzzyaaai127)*xyzzyaaam127(1:perio&
&dicity)
lapf(xyzzyaaab127,xyzzyaaaa127)=lapf(xyzzyaaab127,xyzzyaaaa127)-xyzzya&
&aao127*xyzzyaaap127
enddo
enddo
enddo
end subroutine xyzzyaagv1
subroutine xyzzyaagw1(n,l,a,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: l,a,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa128,xyzzyaaab128
real(dp) xyzzyaaac128,xyzzyaaad128,xyzzyaaae128,xyzzyaaaf128
do xyzzyaaaa128=1,ne
if(xyzzyaaaa128==ieskip)cycle
xyzzyaaac128=exvecs(4,xyzzyaaaa128)
if(xyzzyaaac128>=l)xyzzyaaac128=l
xyzzyaaad128=1.d0/(xyzzyaaac128+a)
xyzzyaaaf128=xyzzyaaac128*xyzzyaaad128
f(1,xyzzyaaaa128)=1.d0
f(2,xyzzyaaaa128)=xyzzyaaaf128
xyzzyaaae128=xyzzyaaaf128
do xyzzyaaab128=3,n
xyzzyaaae128=xyzzyaaae128*xyzzyaaaf128
f(xyzzyaaab128,xyzzyaaaa128)=xyzzyaaae128
enddo
enddo
end subroutine xyzzyaagw1
subroutine xyzzyaagx1(n,l,a,ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: l,a,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa129,xyzzyaaab129
real(dp) xyzzyaaac129,xyzzyaaad129,xyzzyaaae129,xyzzyaaaf129,xyzzyaaag&
&129,xyzzyaaah129,xyzzyaaai129
do xyzzyaaaa129=1,ne
if(xyzzyaaaa129==ieskip)cycle
xyzzyaaac129=exvecs(4,xyzzyaaaa129)
if(xyzzyaaac129>=l)then
xyzzyaaad129=1.d0/(l+a)
xyzzyaaag129=l*xyzzyaaad129
f(1,xyzzyaaaa129)=1.d0
df(1,xyzzyaaaa129)=0.d0
f(2,xyzzyaaaa129)=xyzzyaaag129
df(2,xyzzyaaaa129)=0.d0
xyzzyaaae129=xyzzyaaag129
do xyzzyaaab129=3,n
xyzzyaaae129=xyzzyaaae129*xyzzyaaag129
f(xyzzyaaab129,xyzzyaaaa129)=xyzzyaaae129
df(xyzzyaaab129,xyzzyaaaa129)=0.d0
enddo
else
xyzzyaaad129=1.d0/(xyzzyaaac129+a)
xyzzyaaag129=xyzzyaaac129*xyzzyaaad129
xyzzyaaah129=a*xyzzyaaad129*xyzzyaaad129
f(1,xyzzyaaaa129)=1.d0
df(1,xyzzyaaaa129)=0.d0
f(2,xyzzyaaaa129)=xyzzyaaag129
df(2,xyzzyaaaa129)=xyzzyaaah129
xyzzyaaae129=xyzzyaaag129
xyzzyaaaf129=1.d0
xyzzyaaai129=xyzzyaaah129
do xyzzyaaab129=3,n
xyzzyaaaf129=xyzzyaaae129
xyzzyaaae129=xyzzyaaae129*xyzzyaaag129
f(xyzzyaaab129,xyzzyaaaa129)=xyzzyaaae129
xyzzyaaai129=xyzzyaaai129+xyzzyaaah129
df(xyzzyaaab129,xyzzyaaaa129)=xyzzyaaai129*xyzzyaaaf129
enddo
endif
enddo
end subroutine xyzzyaagx1
subroutine xyzzyaagy1(n,l,a,ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: l,a,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa130,xyzzyaaab130
real(dp) xyzzyaaac130,xyzzyaaad130,xyzzyaaae130,xyzzyaaaf130,xyzzyaaag&
&130,xyzzyaaah130,xyzzyaaai130,xyzzyaaaj130,xyzzyaaak130,xyzzyaaal130
do xyzzyaaaa130=1,ne
if(xyzzyaaaa130==ieskip)cycle
xyzzyaaac130=exvecs(4,xyzzyaaaa130)
if(xyzzyaaac130>=l)then
xyzzyaaad130=1.d0/(l+a)
xyzzyaaai130=l*xyzzyaaad130
f(1,xyzzyaaaa130)=1.d0
df(1,xyzzyaaaa130)=0.d0
d2f(1,xyzzyaaaa130)=0.d0
f(2,xyzzyaaaa130)=xyzzyaaai130
df(2,xyzzyaaaa130)=0.d0
d2f(2,xyzzyaaaa130)=0.d0
xyzzyaaaf130=xyzzyaaai130
do xyzzyaaab130=3,n
xyzzyaaaf130=xyzzyaaaf130*xyzzyaaai130
f(xyzzyaaab130,xyzzyaaaa130)=xyzzyaaaf130
df(xyzzyaaab130,xyzzyaaaa130)=0.d0
d2f(xyzzyaaab130,xyzzyaaaa130)=0.d0
enddo
else
xyzzyaaad130=1.d0/(xyzzyaaac130+a)
xyzzyaaai130=xyzzyaaac130*xyzzyaaad130
xyzzyaaae130=xyzzyaaad130*xyzzyaaad130
xyzzyaaaj130=a*xyzzyaaae130
f(1,xyzzyaaaa130)=1.d0
df(1,xyzzyaaaa130)=0.d0
d2f(1,xyzzyaaaa130)=0.d0
f(2,xyzzyaaaa130)=xyzzyaaai130
df(2,xyzzyaaaa130)=xyzzyaaaj130
d2f(2,xyzzyaaaa130)=-(xyzzyaaaj130+xyzzyaaaj130)*xyzzyaaad130
xyzzyaaaf130=xyzzyaaai130
xyzzyaaag130=1.d0
xyzzyaaah130=0.d0
xyzzyaaak130=xyzzyaaaj130
xyzzyaaal130=-(xyzzyaaac130+xyzzyaaac130)*xyzzyaaae130
do xyzzyaaab130=3,n
xyzzyaaah130=xyzzyaaag130
xyzzyaaag130=xyzzyaaaf130
xyzzyaaaf130=xyzzyaaaf130*xyzzyaaai130
f(xyzzyaaab130,xyzzyaaaa130)=xyzzyaaaf130
xyzzyaaak130=xyzzyaaak130+xyzzyaaaj130
df(xyzzyaaab130,xyzzyaaaa130)=xyzzyaaak130*xyzzyaaag130
xyzzyaaal130=xyzzyaaal130+xyzzyaaaj130
d2f(xyzzyaaab130,xyzzyaaaa130)=xyzzyaaak130*xyzzyaaal130*xyzzyaaah130
enddo
endif
enddo
end subroutine xyzzyaagy1
subroutine xyzzyaagz1(n,p,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: p(2),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa131,xyzzyaaab131
real(dp) xyzzyaaac131,xyzzyaaad131,xyzzyaaae131,xyzzyaaaf131,xyzzyaaag&
&131,xyzzyaaah131,xyzzyaaai131
xyzzyaaac131=p(1)
xyzzyaaad131=p(2)
do xyzzyaaaa131=1,ne
if(xyzzyaaaa131==ieskip)cycle
xyzzyaaae131=exvecs(4,xyzzyaaaa131)
xyzzyaaaf131=xyzzyaaae131**xyzzyaaad131
xyzzyaaag131=1.d0/(xyzzyaaaf131+xyzzyaaac131)
xyzzyaaah131=xyzzyaaae131*xyzzyaaag131
f(1,xyzzyaaaa131)=1.d0
f(2,xyzzyaaaa131)=xyzzyaaah131
xyzzyaaai131=xyzzyaaah131
do xyzzyaaab131=3,n
xyzzyaaai131=xyzzyaaai131*xyzzyaaah131
f(xyzzyaaab131,xyzzyaaaa131)=xyzzyaaai131
enddo
enddo
end subroutine xyzzyaagz1
subroutine xyzzyaaha1(n,p,ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: p(2),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa132,xyzzyaaab132
real(dp) xyzzyaaac132,xyzzyaaad132,xyzzyaaae132,xyzzyaaaf132,xyzzyaaag&
&132,xyzzyaaah132,xyzzyaaai132,xyzzyaaaj132,xyzzyaaak132,xyzzyaaal132
xyzzyaaac132=p(1)
xyzzyaaad132=p(2)
do xyzzyaaaa132=1,ne
if(xyzzyaaaa132==ieskip)cycle
xyzzyaaae132=exvecs(4,xyzzyaaaa132)
xyzzyaaaf132=xyzzyaaae132**xyzzyaaad132
xyzzyaaag132=1.d0/(xyzzyaaaf132+xyzzyaaac132)
xyzzyaaah132=xyzzyaaag132*xyzzyaaag132
xyzzyaaai132=xyzzyaaae132*xyzzyaaag132
xyzzyaaaj132=((1.d0-xyzzyaaad132)*xyzzyaaaf132+xyzzyaaac132)*xyzzyaaah&
&132
f(1,xyzzyaaaa132)=1.d0
df(1,xyzzyaaaa132)=0.d0
f(2,xyzzyaaaa132)=xyzzyaaai132
df(2,xyzzyaaaa132)=xyzzyaaaj132
xyzzyaaak132=xyzzyaaai132
xyzzyaaal132=1.d0
do xyzzyaaab132=3,n
xyzzyaaal132=xyzzyaaak132*(xyzzyaaab132-1)
xyzzyaaak132=xyzzyaaak132*xyzzyaaai132
f(xyzzyaaab132,xyzzyaaaa132)=xyzzyaaak132
df(xyzzyaaab132,xyzzyaaaa132)=xyzzyaaal132*xyzzyaaaj132
enddo
enddo
end subroutine xyzzyaaha1
subroutine xyzzyaahb1(n,p,ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: p(2),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa133,xyzzyaaab133
real(dp) xyzzyaaac133,xyzzyaaad133,xyzzyaaae133,xyzzyaaaf133,xyzzyaaag&
&133,xyzzyaaah133,xyzzyaaai133,xyzzyaaaj133,xyzzyaaak133,xyzzyaaal133,&
&xyzzyaaam133,xyzzyaaan133,xyzzyaaao133,xyzzyaaap133,xyzzyaaaq133
xyzzyaaac133=p(1)
xyzzyaaad133=p(2)
do xyzzyaaaa133=1,ne
if(xyzzyaaaa133==ieskip)cycle
xyzzyaaae133=exvecs(4,xyzzyaaaa133)
xyzzyaaaf133=xyzzyaaae133**(xyzzyaaad133-1.d0)
xyzzyaaag133=xyzzyaaaf133*xyzzyaaae133
xyzzyaaah133=1.d0/(xyzzyaaag133+xyzzyaaac133)
xyzzyaaai133=xyzzyaaah133*xyzzyaaah133
xyzzyaaaj133=xyzzyaaai133*xyzzyaaah133
xyzzyaaak133=xyzzyaaae133*xyzzyaaah133
xyzzyaaal133=((1.d0-xyzzyaaad133)*xyzzyaaag133+xyzzyaaac133)*xyzzyaaai&
&133
xyzzyaaam133=xyzzyaaal133*xyzzyaaal133
xyzzyaaan133=xyzzyaaad133*xyzzyaaaf133*xyzzyaaaj133*(xyzzyaaad133*(xyz&
&zyaaag133-xyzzyaaac133)-xyzzyaaag133-xyzzyaaac133)
f(1,xyzzyaaaa133)=1.d0
df(1,xyzzyaaaa133)=0.d0
d2f(1,xyzzyaaaa133)=0.d0
f(2,xyzzyaaaa133)=xyzzyaaak133
df(2,xyzzyaaaa133)=xyzzyaaal133
d2f(2,xyzzyaaaa133)=xyzzyaaan133
xyzzyaaao133=xyzzyaaak133
xyzzyaaap133=1.d0
xyzzyaaaq133=0.d0
do xyzzyaaab133=3,n
xyzzyaaaq133=xyzzyaaap133*(xyzzyaaab133-1)
xyzzyaaap133=xyzzyaaao133*(xyzzyaaab133-1)
xyzzyaaao133=xyzzyaaao133*xyzzyaaak133
f(xyzzyaaab133,xyzzyaaaa133)=xyzzyaaao133
df(xyzzyaaab133,xyzzyaaaa133)=xyzzyaaap133*xyzzyaaal133
d2f(xyzzyaaab133,xyzzyaaaa133)=xyzzyaaaq133*xyzzyaaam133+xyzzyaaap133*&
&xyzzyaaan133
enddo
enddo
end subroutine xyzzyaahb1
subroutine xyzzyaahc1(n,a,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: a,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa134,xyzzyaaab134
real(dp) xyzzyaaac134,xyzzyaaad134,xyzzyaaae134
do xyzzyaaaa134=1,ne
if(xyzzyaaaa134==ieskip)cycle
xyzzyaaac134=exvecs(4,xyzzyaaaa134)
xyzzyaaae134=1.d0/(xyzzyaaac134+a)
f(1,xyzzyaaaa134)=1.d0
f(2,xyzzyaaaa134)=xyzzyaaae134
xyzzyaaad134=xyzzyaaae134
do xyzzyaaab134=3,n
xyzzyaaad134=xyzzyaaad134*xyzzyaaae134
f(xyzzyaaab134,xyzzyaaaa134)=xyzzyaaad134
enddo
enddo
end subroutine xyzzyaahc1
subroutine xyzzyaahd1(n,a,ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: a,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa135,xyzzyaaab135
real(dp) xyzzyaaac135,xyzzyaaad135,xyzzyaaae135,xyzzyaaaf135,xyzzyaaag&
&135
do xyzzyaaaa135=1,ne
if(xyzzyaaaa135==ieskip)cycle
xyzzyaaac135=exvecs(4,xyzzyaaaa135)
xyzzyaaaf135=1.d0/(xyzzyaaac135+a)
xyzzyaaag135=-xyzzyaaaf135*xyzzyaaaf135
f(1,xyzzyaaaa135)=1.d0
df(1,xyzzyaaaa135)=0.d0
f(2,xyzzyaaaa135)=xyzzyaaaf135
df(2,xyzzyaaaa135)=xyzzyaaag135
xyzzyaaad135=xyzzyaaaf135
xyzzyaaae135=1.d0
do xyzzyaaab135=3,n
xyzzyaaae135=xyzzyaaad135*(xyzzyaaab135-1)
xyzzyaaad135=xyzzyaaad135*xyzzyaaaf135
f(xyzzyaaab135,xyzzyaaaa135)=xyzzyaaad135
df(xyzzyaaab135,xyzzyaaaa135)=xyzzyaaae135*xyzzyaaag135
enddo
enddo
end subroutine xyzzyaahd1
subroutine xyzzyaahe1(n,a,ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: a,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa136,xyzzyaaab136
real(dp) xyzzyaaac136,xyzzyaaad136,xyzzyaaae136,xyzzyaaaf136,xyzzyaaag&
&136,xyzzyaaah136,xyzzyaaai136,xyzzyaaaj136
do xyzzyaaaa136=1,ne
if(xyzzyaaaa136==ieskip)cycle
xyzzyaaac136=exvecs(4,xyzzyaaaa136)
xyzzyaaag136=1.d0/(xyzzyaaac136+a)
xyzzyaaah136=-xyzzyaaag136*xyzzyaaag136
xyzzyaaai136=xyzzyaaah136*xyzzyaaah136
xyzzyaaaj136=-2.d0*xyzzyaaag136*xyzzyaaah136
f(1,xyzzyaaaa136)=1.d0
df(1,xyzzyaaaa136)=0.d0
d2f(1,xyzzyaaaa136)=0.d0
f(2,xyzzyaaaa136)=xyzzyaaag136
df(2,xyzzyaaaa136)=xyzzyaaah136
d2f(2,xyzzyaaaa136)=xyzzyaaaj136
xyzzyaaad136=xyzzyaaag136
xyzzyaaae136=1.d0
xyzzyaaaf136=0.d0
do xyzzyaaab136=3,n
xyzzyaaaf136=xyzzyaaae136*(xyzzyaaab136-1)
xyzzyaaae136=xyzzyaaad136*(xyzzyaaab136-1)
xyzzyaaad136=xyzzyaaad136*xyzzyaaag136
f(xyzzyaaab136,xyzzyaaaa136)=xyzzyaaad136
df(xyzzyaaab136,xyzzyaaaa136)=xyzzyaaae136*xyzzyaaah136
d2f(xyzzyaaab136,xyzzyaaaa136)=xyzzyaaaf136*xyzzyaaai136+xyzzyaaae136*&
&xyzzyaaaj136
enddo
enddo
end subroutine xyzzyaahe1
subroutine xyzzyaahf1(n,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa137,xyzzyaaab137,xyzzyaaac137,xyzzyaaad137,xyzzyaaae1&
&37
real(dp) xyzzyaaaf137,xyzzyaaag137(3),xyzzyaaah137(3)
xyzzyaaac137=n/dimensionality
call xyzzyaagn1(xyzzyaaac137,ne,ieskip,exvecs,fstride,f)
do xyzzyaaaa137=1,ne
if(xyzzyaaaa137==ieskip)cycle
if(xyzzyaaaa137>ieskip)then
xyzzyaaag137(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa137)
else
xyzzyaaag137(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa137)
endif
xyzzyaaaf137=exvecs(4,xyzzyaaaa137)
if(xyzzyaaaf137==0.d0)then
xyzzyaaah137=0.d0
else
xyzzyaaah137(1:dimensionality)=xyzzyaaag137(1:dimensionality)/xyzzyaaa&
&f137
endif
xyzzyaaad137=n
do xyzzyaaab137=xyzzyaaac137,1,-1
do xyzzyaaae137=dimensionality,1,-1
f(xyzzyaaad137,xyzzyaaaa137)=f(xyzzyaaab137,xyzzyaaaa137)*xyzzyaaah137&
&(xyzzyaaae137)
xyzzyaaad137=xyzzyaaad137-1
enddo
enddo
enddo
end subroutine xyzzyaahf1
subroutine xyzzyaahg1(n,ne,ieskip,exvecs,fstride,f,gradf)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*)
integer xyzzyaaaa138,xyzzyaaab138,xyzzyaaac138,xyzzyaaad138,xyzzyaaae1&
&38
real(dp) xyzzyaaaf138,xyzzyaaag138(3),xyzzyaaah138(3),xyzzyaaai138(3,3&
&),xyzzyaaaj138,xyzzyaaak138(n/dimensionality,ne),xyzzyaaal138(3)
xyzzyaaac138=n/dimensionality
call xyzzyaago1(xyzzyaaac138,ne,ieskip,exvecs,fstride,f,xyzzyaaak138)
do xyzzyaaaa138=1,ne
if(xyzzyaaaa138==ieskip)cycle
if(xyzzyaaaa138>ieskip)then
xyzzyaaag138(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa138)
else
xyzzyaaag138(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa138)
endif
xyzzyaaag138(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa138)
xyzzyaaaf138=exvecs(4,xyzzyaaaa138)
call xyzzyaaiz1(xyzzyaaaf138,exvecs(1,xyzzyaaaa138),xyzzyaaal138)
if(xyzzyaaaf138==0.d0)then
xyzzyaaah138=0.d0
xyzzyaaai138(1:dimensionality,1:dimensionality)=0.d0
else
xyzzyaaaj138=1.d0/xyzzyaaaf138
xyzzyaaah138(1:dimensionality)=xyzzyaaag138(1:dimensionality)*xyzzyaaa&
&j138
do xyzzyaaae138=1,dimensionality
xyzzyaaai138(1:dimensionality,xyzzyaaae138)=-xyzzyaaah138(1:dimensiona&
&lity)*xyzzyaaah138(xyzzyaaae138)*xyzzyaaaj138
xyzzyaaai138(xyzzyaaae138,xyzzyaaae138)=xyzzyaaai138(xyzzyaaae138,xyzz&
&yaaae138)+xyzzyaaaj138
enddo
endif
xyzzyaaad138=n
do xyzzyaaab138=xyzzyaaac138,1,-1
do xyzzyaaae138=dimensionality,1,-1
gradf(1:dimensionality,xyzzyaaad138,xyzzyaaaa138)=xyzzyaaal138(1:dimen&
&sionality)*xyzzyaaak138(xyzzyaaab138,xyzzyaaaa138)*xyzzyaaah138(xyzzy&
&aaae138)+f(xyzzyaaab138,xyzzyaaaa138)*xyzzyaaai138(1:dimensionality,x&
&yzzyaaae138)
f(xyzzyaaad138,xyzzyaaaa138)=f(xyzzyaaab138,xyzzyaaaa138)*xyzzyaaah138&
&(xyzzyaaae138)
xyzzyaaad138=xyzzyaaad138-1
enddo
enddo
enddo
end subroutine xyzzyaahg1
subroutine xyzzyaahh1(n,ne,ieskip,exvecs,fstride,f,gradf,lapf)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),lapf(fstride&
&,*)
integer xyzzyaaaa139,xyzzyaaab139,xyzzyaaac139,xyzzyaaad139,xyzzyaaae1&
&39
real(dp) xyzzyaaaf139,xyzzyaaag139(3),xyzzyaaah139(3),xyzzyaaai139(3,3&
&),xyzzyaaaj139(3),xyzzyaaak139,xyzzyaaal139(n/dimensionality,ne),xyzz&
&yaaam139(n/dimensionality,ne),xyzzyaaan139(3),xyzzyaaao139
xyzzyaaac139=n/dimensionality
call xyzzyaagp1(xyzzyaaac139,ne,ieskip,exvecs,fstride,f,xyzzyaaal139,x&
&yzzyaaam139)
do xyzzyaaaa139=1,ne
if(xyzzyaaaa139==ieskip)cycle
if(xyzzyaaaa139>ieskip)then
xyzzyaaag139(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa139)
else
xyzzyaaag139(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa139)
endif
xyzzyaaaf139=exvecs(4,xyzzyaaaa139)
call xyzzyaaiy1(xyzzyaaaf139,exvecs(1,xyzzyaaaa139),xyzzyaaan139,xyzzy&
&aaao139)
if(xyzzyaaaf139==0.d0)then
xyzzyaaah139=0.d0
xyzzyaaai139(1:dimensionality,1:dimensionality)=0.d0
xyzzyaaaj139=0.d0
else
xyzzyaaak139=1.d0/xyzzyaaaf139
xyzzyaaah139(1:dimensionality)=xyzzyaaag139(1:dimensionality)*xyzzyaaa&
&k139
do xyzzyaaae139=1,dimensionality
xyzzyaaai139(1:dimensionality,xyzzyaaae139)=-xyzzyaaah139(1:dimensiona&
&lity)*xyzzyaaah139(xyzzyaaae139)*xyzzyaaak139
xyzzyaaai139(xyzzyaaae139,xyzzyaaae139)=xyzzyaaai139(xyzzyaaae139,xyzz&
&yaaae139)+xyzzyaaak139
enddo
xyzzyaaaj139(1:dimensionality)=dble(1-dimensionality)*xyzzyaaak139*xyz&
&zyaaak139*xyzzyaaah139(1:dimensionality)
endif
xyzzyaaad139=n
do xyzzyaaab139=xyzzyaaac139,1,-1
do xyzzyaaae139=dimensionality,1,-1
lapf(xyzzyaaad139,xyzzyaaaa139)=(xyzzyaaam139(xyzzyaaab139,xyzzyaaaa13&
&9)+xyzzyaaal139(xyzzyaaab139,xyzzyaaaa139)*xyzzyaaao139)*xyzzyaaah139&
&(xyzzyaaae139)+f(xyzzyaaab139,xyzzyaaaa139)*xyzzyaaaj139(xyzzyaaae139&
&)
gradf(1:dimensionality,xyzzyaaad139,xyzzyaaaa139)=xyzzyaaan139(1:dimen&
&sionality)*xyzzyaaal139(xyzzyaaab139,xyzzyaaaa139)*xyzzyaaah139(xyzzy&
&aaae139)+f(xyzzyaaab139,xyzzyaaaa139)*xyzzyaaai139(1:dimensionality,x&
&yzzyaaae139)
f(xyzzyaaad139,xyzzyaaaa139)=f(xyzzyaaab139,xyzzyaaaa139)*xyzzyaaah139&
&(xyzzyaaae139)
xyzzyaaad139=xyzzyaaad139-1
enddo
enddo
enddo
end subroutine xyzzyaahh1
subroutine xyzzyaahi1(n,k0,nk,p,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,k0,nk(*),ne,ieskip,fstride
real(dp),intent(in) :: p(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa140,xyzzyaaab140,xyzzyaaac140,xyzzyaaad140
real(dp) xyzzyaaae140,xyzzyaaaf140
do xyzzyaaaa140=1,ne
if(xyzzyaaaa140==ieskip)cycle
xyzzyaaae140=exvecs(4,xyzzyaaaa140)
select case(k0)
case(0)
xyzzyaaaf140=1.d0
case(1)
xyzzyaaaf140=xyzzyaaae140
case default
xyzzyaaaf140=xyzzyaaae140**k0
end select
xyzzyaaac140=0
do xyzzyaaab140=1,n
f(xyzzyaaab140,xyzzyaaaa140)=xyzzyaaaf140
xyzzyaaaf140=xyzzyaaaf140*xyzzyaaae140
do xyzzyaaad140=2,nk(xyzzyaaab140)
xyzzyaaac140=xyzzyaaac140+1
f(xyzzyaaab140,xyzzyaaaa140)=f(xyzzyaaab140,xyzzyaaaa140)+p(xyzzyaaac1&
&40)*xyzzyaaaf140
xyzzyaaaf140=xyzzyaaaf140*xyzzyaaae140
enddo
enddo
enddo
end subroutine xyzzyaahi1
subroutine xyzzyaahj1(n,k0,nk,p,ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: n,k0,nk(*),ne,ieskip,fstride
real(dp),intent(in) :: p(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa141,xyzzyaaab141,xyzzyaaac141,xyzzyaaad141,xyzzyaaae1&
&41
real(dp) xyzzyaaaf141,xyzzyaaag141,xyzzyaaah141,xyzzyaaai141
do xyzzyaaaa141=1,ne
if(xyzzyaaaa141==ieskip)cycle
xyzzyaaaf141=exvecs(4,xyzzyaaaa141)
select case(k0)
case(0)
xyzzyaaah141=0.d0
xyzzyaaag141=1.d0
case(1)
xyzzyaaah141=1.d0
xyzzyaaag141=xyzzyaaaf141
case default
xyzzyaaah141=xyzzyaaaf141**(k0-1)
xyzzyaaag141=xyzzyaaah141*xyzzyaaaf141
xyzzyaaah141=xyzzyaaah141*k0
end select
xyzzyaaac141=0
xyzzyaaae141=k0
do xyzzyaaab141=1,n
f(xyzzyaaab141,xyzzyaaaa141)=xyzzyaaag141
df(xyzzyaaab141,xyzzyaaaa141)=xyzzyaaah141
xyzzyaaae141=xyzzyaaae141+1
xyzzyaaah141=xyzzyaaag141*xyzzyaaae141
xyzzyaaag141=xyzzyaaag141*xyzzyaaaf141
do xyzzyaaad141=2,nk(xyzzyaaab141)
xyzzyaaac141=xyzzyaaac141+1
xyzzyaaai141=p(xyzzyaaac141)
f(xyzzyaaab141,xyzzyaaaa141)=f(xyzzyaaab141,xyzzyaaaa141)+xyzzyaaai141&
&*xyzzyaaag141
df(xyzzyaaab141,xyzzyaaaa141)=df(xyzzyaaab141,xyzzyaaaa141)+xyzzyaaai1&
&41*xyzzyaaah141
xyzzyaaae141=xyzzyaaae141+1
xyzzyaaah141=xyzzyaaag141*xyzzyaaae141
xyzzyaaag141=xyzzyaaag141*xyzzyaaaf141
enddo
enddo
enddo
end subroutine xyzzyaahj1
subroutine xyzzyaahk1(n,k0,nk,p,ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: n,k0,nk(*),ne,ieskip,fstride
real(dp),intent(in) :: p(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa142,xyzzyaaab142,xyzzyaaac142,xyzzyaaad142,xyzzyaaae1&
&42
real(dp) xyzzyaaaf142,xyzzyaaag142,xyzzyaaah142,xyzzyaaai142,xyzzyaaaj&
&142
do xyzzyaaaa142=1,ne
if(xyzzyaaaa142==ieskip)cycle
xyzzyaaaf142=exvecs(4,xyzzyaaaa142)
select case(k0)
case(0)
xyzzyaaai142=0.d0
xyzzyaaah142=0.d0
xyzzyaaag142=1.d0
case(1)
xyzzyaaai142=0.d0
xyzzyaaah142=1.d0
xyzzyaaag142=xyzzyaaaf142
case default
xyzzyaaai142=xyzzyaaaf142**(k0-2)
xyzzyaaah142=xyzzyaaai142*xyzzyaaaf142
xyzzyaaai142=xyzzyaaai142*k0*(k0-1)
xyzzyaaag142=xyzzyaaah142*xyzzyaaaf142
xyzzyaaah142=xyzzyaaah142*k0
end select
xyzzyaaac142=0
xyzzyaaae142=k0
do xyzzyaaab142=1,n
f(xyzzyaaab142,xyzzyaaaa142)=xyzzyaaag142
df(xyzzyaaab142,xyzzyaaaa142)=xyzzyaaah142
d2f(xyzzyaaab142,xyzzyaaaa142)=xyzzyaaai142
xyzzyaaae142=xyzzyaaae142+1
xyzzyaaai142=xyzzyaaah142*xyzzyaaae142
xyzzyaaah142=xyzzyaaag142*xyzzyaaae142
xyzzyaaag142=xyzzyaaag142*xyzzyaaaf142
do xyzzyaaad142=2,nk(xyzzyaaab142)
xyzzyaaac142=xyzzyaaac142+1
xyzzyaaaj142=p(xyzzyaaac142)
f(xyzzyaaab142,xyzzyaaaa142)=f(xyzzyaaab142,xyzzyaaaa142)+xyzzyaaaj142&
&*xyzzyaaag142
df(xyzzyaaab142,xyzzyaaaa142)=df(xyzzyaaab142,xyzzyaaaa142)+xyzzyaaaj1&
&42*xyzzyaaah142
d2f(xyzzyaaab142,xyzzyaaaa142)=d2f(xyzzyaaab142,xyzzyaaaa142)+xyzzyaaa&
&j142*xyzzyaaai142
xyzzyaaae142=xyzzyaaae142+1
xyzzyaaai142=xyzzyaaah142*xyzzyaaae142
xyzzyaaah142=xyzzyaaag142*xyzzyaaae142
xyzzyaaag142=xyzzyaaag142*xyzzyaaaf142
enddo
enddo
enddo
end subroutine xyzzyaahk1
subroutine xyzzyaahl1(n,k0,nk,p,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,k0,nk(*),ne,ieskip,fstride
real(dp),intent(in) :: p(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa143,xyzzyaaab143,xyzzyaaac143,xyzzyaaad143,xyzzyaaae1&
&43
real(dp) xyzzyaaaf143,xyzzyaaag143(3),xyzzyaaah143(3)
xyzzyaaac143=n/dimensionality
call xyzzyaahi1(xyzzyaaac143,k0,nk,p,ne,ieskip,exvecs,fstride,f)
do xyzzyaaaa143=1,ne
if(xyzzyaaaa143==ieskip)cycle
if(xyzzyaaaa143>ieskip)then
xyzzyaaag143(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa143)
else
xyzzyaaag143(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa143)
endif
xyzzyaaaf143=exvecs(4,xyzzyaaaa143)
if(xyzzyaaaf143==0.d0)then
xyzzyaaah143=0.d0
else
xyzzyaaah143(1:dimensionality)=xyzzyaaag143(1:dimensionality)/xyzzyaaa&
&f143
endif
xyzzyaaad143=n
do xyzzyaaab143=xyzzyaaac143,1,-1
do xyzzyaaae143=dimensionality,1,-1
f(xyzzyaaad143,xyzzyaaaa143)=f(xyzzyaaab143,xyzzyaaaa143)*xyzzyaaah143&
&(xyzzyaaae143)
xyzzyaaad143=xyzzyaaad143-1
enddo
enddo
enddo
end subroutine xyzzyaahl1
subroutine xyzzyaahm1(n,k0,nk,p,ne,ieskip,exvecs,fstride,f,gradf)
implicit none
integer,intent(in) :: n,k0,nk(*),ne,ieskip,fstride
real(dp),intent(in) :: p(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*)
integer xyzzyaaaa144,xyzzyaaab144,xyzzyaaac144,xyzzyaaad144,xyzzyaaae1&
&44
real(dp) xyzzyaaaf144,xyzzyaaag144(3),xyzzyaaah144(3),xyzzyaaai144(3,3&
&),xyzzyaaaj144,xyzzyaaak144(n/dimensionality,ne),xyzzyaaal144(3)
xyzzyaaac144=n/dimensionality
call xyzzyaahj1(xyzzyaaac144,k0,nk,p,ne,ieskip,exvecs,fstride,f,xyzzya&
&aak144)
do xyzzyaaaa144=1,ne
if(xyzzyaaaa144==ieskip)cycle
if(xyzzyaaaa144>ieskip)then
xyzzyaaag144(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa144)
else
xyzzyaaag144(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa144)
endif
xyzzyaaag144(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa144)
xyzzyaaaf144=exvecs(4,xyzzyaaaa144)
call xyzzyaaiz1(xyzzyaaaf144,exvecs(1,xyzzyaaaa144),xyzzyaaal144)
if(xyzzyaaaf144==0.d0)then
xyzzyaaah144=0.d0
xyzzyaaai144(1:dimensionality,1:dimensionality)=0.d0
else
xyzzyaaaj144=1.d0/xyzzyaaaf144
xyzzyaaah144(1:dimensionality)=xyzzyaaag144(1:dimensionality)*xyzzyaaa&
&j144
do xyzzyaaae144=1,dimensionality
xyzzyaaai144(1:dimensionality,xyzzyaaae144)=-xyzzyaaah144(1:dimensiona&
&lity)*xyzzyaaah144(xyzzyaaae144)*xyzzyaaaj144
xyzzyaaai144(xyzzyaaae144,xyzzyaaae144)=xyzzyaaai144(xyzzyaaae144,xyzz&
&yaaae144)+xyzzyaaaj144
enddo
endif
xyzzyaaad144=n
do xyzzyaaab144=xyzzyaaac144,1,-1
do xyzzyaaae144=dimensionality,1,-1
gradf(1:dimensionality,xyzzyaaad144,xyzzyaaaa144)=xyzzyaaal144(1:dimen&
&sionality)*xyzzyaaak144(xyzzyaaab144,xyzzyaaaa144)*xyzzyaaah144(xyzzy&
&aaae144)+f(xyzzyaaab144,xyzzyaaaa144)*xyzzyaaai144(1:dimensionality,x&
&yzzyaaae144)
f(xyzzyaaad144,xyzzyaaaa144)=f(xyzzyaaab144,xyzzyaaaa144)*xyzzyaaah144&
&(xyzzyaaae144)
xyzzyaaad144=xyzzyaaad144-1
enddo
enddo
enddo
end subroutine xyzzyaahm1
subroutine xyzzyaahn1(n,k0,nk,p,ne,ieskip,exvecs,fstride,f,gradf,lapf)
implicit none
integer,intent(in) :: n,k0,nk(*),ne,ieskip,fstride
real(dp),intent(in) :: p(*),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),lapf(fstride&
&,*)
integer xyzzyaaaa145,xyzzyaaab145,xyzzyaaac145,xyzzyaaad145,xyzzyaaae1&
&45
real(dp) xyzzyaaaf145,xyzzyaaag145(3),xyzzyaaah145(3),xyzzyaaai145(3,3&
&),xyzzyaaaj145(3),xyzzyaaak145,xyzzyaaal145(n/dimensionality,ne),xyzz&
&yaaam145(n/dimensionality,ne),xyzzyaaan145(3),xyzzyaaao145
xyzzyaaac145=n/dimensionality
call xyzzyaahk1(xyzzyaaac145,k0,nk,p,ne,ieskip,exvecs,fstride,f,xyzzya&
&aal145,xyzzyaaam145)
do xyzzyaaaa145=1,ne
if(xyzzyaaaa145==ieskip)cycle
if(xyzzyaaaa145>ieskip)then
xyzzyaaag145(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa145)
else
xyzzyaaag145(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa145)
endif
xyzzyaaaf145=exvecs(4,xyzzyaaaa145)
call xyzzyaaiy1(xyzzyaaaf145,exvecs(1,xyzzyaaaa145),xyzzyaaan145,xyzzy&
&aaao145)
if(xyzzyaaaf145==0.d0)then
xyzzyaaah145=0.d0
xyzzyaaai145(1:dimensionality,1:dimensionality)=0.d0
xyzzyaaaj145=0.d0
else
xyzzyaaak145=1.d0/xyzzyaaaf145
xyzzyaaah145(1:dimensionality)=xyzzyaaag145(1:dimensionality)*xyzzyaaa&
&k145
do xyzzyaaae145=1,dimensionality
xyzzyaaai145(1:dimensionality,xyzzyaaae145)=-xyzzyaaah145(1:dimensiona&
&lity)*xyzzyaaah145(xyzzyaaae145)*xyzzyaaak145
xyzzyaaai145(xyzzyaaae145,xyzzyaaae145)=xyzzyaaai145(xyzzyaaae145,xyzz&
&yaaae145)+xyzzyaaak145
enddo
xyzzyaaaj145(1:dimensionality)=dble(1-dimensionality)*xyzzyaaak145*xyz&
&zyaaak145*xyzzyaaah145(1:dimensionality)
endif
xyzzyaaad145=n
do xyzzyaaab145=xyzzyaaac145,1,-1
do xyzzyaaae145=dimensionality,1,-1
lapf(xyzzyaaad145,xyzzyaaaa145)=(xyzzyaaam145(xyzzyaaab145,xyzzyaaaa14&
&5)+xyzzyaaal145(xyzzyaaab145,xyzzyaaaa145)*xyzzyaaao145)*xyzzyaaah145&
&(xyzzyaaae145)+f(xyzzyaaab145,xyzzyaaaa145)*xyzzyaaaj145(xyzzyaaae145&
&)
gradf(1:dimensionality,xyzzyaaad145,xyzzyaaaa145)=xyzzyaaan145(1:dimen&
&sionality)*xyzzyaaal145(xyzzyaaab145,xyzzyaaaa145)*xyzzyaaah145(xyzzy&
&aaae145)+f(xyzzyaaab145,xyzzyaaaa145)*xyzzyaaai145(1:dimensionality,x&
&yzzyaaae145)
f(xyzzyaaad145,xyzzyaaaa145)=f(xyzzyaaab145,xyzzyaaaa145)*xyzzyaaah145&
&(xyzzyaaae145)
xyzzyaaad145=xyzzyaaad145-1
enddo
enddo
enddo
end subroutine xyzzyaahn1
subroutine xyzzyaaho1(p,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa146
real(dp) xyzzyaaab146,xyzzyaaac146,xyzzyaaad146,xyzzyaaae146,xyzzyaaaf&
&146
xyzzyaaad146=p(2)
xyzzyaaae146=p(3)
do xyzzyaaaa146=1,ne
if(xyzzyaaaa146==ieskip)cycle
xyzzyaaab146=exvecs(4,xyzzyaaaa146)
if(xyzzyaaab146==0.d0)then
select case(dimensionality)
case(3)
xyzzyaaaf146=-xyzzyaaad146
case(2)
xyzzyaaaf146=-xyzzyaaae146
end select
else
select case(dimensionality)
case(3)
xyzzyaaaf146=-(1.d0-exp(-xyzzyaaab146*xyzzyaaad146))/xyzzyaaab146
case(2)
xyzzyaaac146=sqrt(xyzzyaaab146)
xyzzyaaaf146=-(1.d0-exp(-0.5d0*xyzzyaaab146*xyzzyaaad146-xyzzyaaac146*&
&xyzzyaaae146))/xyzzyaaac146
end select
endif
f(1,xyzzyaaaa146)=xyzzyaaaf146
enddo
end subroutine xyzzyaaho1
subroutine xyzzyaahp1(p,ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa147
real(dp) xyzzyaaab147,xyzzyaaac147,xyzzyaaad147,xyzzyaaae147,xyzzyaaaf&
&147,xyzzyaaag147,xyzzyaaah147,xyzzyaaai147,xyzzyaaaj147,xyzzyaaak147,&
&xyzzyaaal147
xyzzyaaaf147=p(2)
xyzzyaaag147=p(3)
do xyzzyaaaa147=1,ne
if(xyzzyaaaa147==ieskip)cycle
xyzzyaaab147=exvecs(4,xyzzyaaaa147)
if(xyzzyaaab147==0.d0)then
select case(dimensionality)
case(3)
xyzzyaaak147=-xyzzyaaaf147
xyzzyaaal147=0.5d0*xyzzyaaaf147*xyzzyaaaf147
case(2)
xyzzyaaak147=-xyzzyaaag147
xyzzyaaal147=1.d0/(3.d0*xyzzyaaag147*xyzzyaaag147*xyzzyaaag147)
end select
else
select case(dimensionality)
case(3)
xyzzyaaac147=1.d0/xyzzyaaab147
xyzzyaaah147=xyzzyaaab147*xyzzyaaaf147
xyzzyaaaj147=exp(-xyzzyaaah147)
xyzzyaaak147=-xyzzyaaac147*(1.d0-xyzzyaaaj147)
xyzzyaaal147=xyzzyaaac147*xyzzyaaac147*(1.d0-(1.d0+xyzzyaaah147)*xyzzy&
&aaaj147)
case(2)
xyzzyaaac147=1.d0/xyzzyaaab147
xyzzyaaad147=sqrt(xyzzyaaab147)
xyzzyaaae147=1.d0/xyzzyaaad147
xyzzyaaah147=xyzzyaaab147*xyzzyaaaf147
xyzzyaaai147=xyzzyaaad147*xyzzyaaag147
xyzzyaaaj147=exp(-0.5d0*xyzzyaaah147-xyzzyaaai147)
xyzzyaaak147=-xyzzyaaae147*(1.d0-xyzzyaaaj147)
xyzzyaaal147=0.5d0*xyzzyaaae147*xyzzyaaac147*(1.d0-xyzzyaaaj147*(1.d0+&
&xyzzyaaah147+xyzzyaaai147))
end select
endif
f(1,xyzzyaaaa147)=xyzzyaaak147
df(1,xyzzyaaaa147)=xyzzyaaal147
enddo
end subroutine xyzzyaahp1
subroutine xyzzyaahq1(p,ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa148
real(dp) xyzzyaaab148,xyzzyaaac148,xyzzyaaad148,xyzzyaaae148,xyzzyaaaf&
&148,xyzzyaaag148,xyzzyaaah148,xyzzyaaai148,xyzzyaaaj148,xyzzyaaak148,&
&xyzzyaaal148,xyzzyaaam148,xyzzyaaan148,xyzzyaaao148
xyzzyaaag148=p(2)
xyzzyaaah148=p(3)
do xyzzyaaaa148=1,ne
if(xyzzyaaaa148==ieskip)cycle
xyzzyaaab148=exvecs(4,xyzzyaaaa148)
if(xyzzyaaab148==0.d0)then
select case(dimensionality)
case(3)
xyzzyaaam148=-xyzzyaaag148
xyzzyaaan148=0.5d0*xyzzyaaag148*xyzzyaaag148
xyzzyaaao148=-1.d0/(6*xyzzyaaag148*xyzzyaaag148*xyzzyaaag148)
case(2)
xyzzyaaam148=-xyzzyaaah148
xyzzyaaan148=1.d0/(3.d0*xyzzyaaah148*xyzzyaaah148*xyzzyaaah148)
xyzzyaaao148=0.d0
end select
else
select case(dimensionality)
case(3)
xyzzyaaac148=1.d0/xyzzyaaab148
xyzzyaaad148=xyzzyaaac148*xyzzyaaac148
xyzzyaaai148=xyzzyaaab148*xyzzyaaag148
xyzzyaaaj148=1.d0+xyzzyaaai148
xyzzyaaal148=exp(-xyzzyaaai148)
xyzzyaaam148=-xyzzyaaac148*(1.d0-xyzzyaaal148)
xyzzyaaan148=xyzzyaaad148*(1.d0-xyzzyaaaj148*xyzzyaaal148)
xyzzyaaao148=-xyzzyaaad148*xyzzyaaac148*(2.d0-(1.d0+xyzzyaaaj148*xyzzy&
&aaaj148)*xyzzyaaal148)
case(2)
xyzzyaaac148=1.d0/xyzzyaaab148
xyzzyaaae148=sqrt(xyzzyaaab148)
xyzzyaaaf148=1.d0/xyzzyaaae148
xyzzyaaai148=xyzzyaaab148*xyzzyaaag148
xyzzyaaak148=xyzzyaaae148*xyzzyaaah148
xyzzyaaal148=exp(-0.5d0*xyzzyaaai148-xyzzyaaak148)
xyzzyaaam148=-xyzzyaaaf148*(1.d0-xyzzyaaal148)
xyzzyaaan148=0.5d0*xyzzyaaaf148*xyzzyaaac148*(1.d0-xyzzyaaal148*(1.d0+&
&xyzzyaaai148+xyzzyaaak148))
xyzzyaaao148=-0.25d0*xyzzyaaaf148*xyzzyaaac148*xyzzyaaac148*(3.d0-(3.d&
&0+3*(xyzzyaaak148+xyzzyaaai148)+(2*xyzzyaaak148+xyzzyaaai148)*xyzzyaa&
&ai148)*xyzzyaaal148)
end select
endif
f(1,xyzzyaaaa148)=xyzzyaaam148
df(1,xyzzyaaaa148)=xyzzyaaan148
d2f(1,xyzzyaaaa148)=xyzzyaaao148
enddo
end subroutine xyzzyaahq1
subroutine xyzzyaahr1(ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa149
real(dp) xyzzyaaab149
do xyzzyaaaa149=1,ne
if(xyzzyaaaa149==ieskip)cycle
xyzzyaaab149=exvecs(4,xyzzyaaaa149)
if(xyzzyaaab149==0.d0)then
f(1,xyzzyaaaa149)=0.d0
else
f(1,xyzzyaaaa149)=xyzzyaaab149*xyzzyaaab149*log(xyzzyaaab149)
endif
enddo
end subroutine xyzzyaahr1
subroutine xyzzyaahs1(ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa150
real(dp) xyzzyaaab150,xyzzyaaac150
do xyzzyaaaa150=1,ne
if(xyzzyaaaa150==ieskip)cycle
xyzzyaaab150=exvecs(4,xyzzyaaaa150)
if(xyzzyaaab150==0.d0)then
f(1,xyzzyaaaa150)=0.d0
df(1,xyzzyaaaa150)=0.d0
else
xyzzyaaac150=xyzzyaaab150*log(xyzzyaaab150)
f(1,xyzzyaaaa150)=xyzzyaaab150*xyzzyaaac150
df(1,xyzzyaaaa150)=2.d0*xyzzyaaac150+xyzzyaaab150
endif
enddo
end subroutine xyzzyaahs1
subroutine xyzzyaaht1(ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa151
real(dp) xyzzyaaab151,xyzzyaaac151,xyzzyaaad151
do xyzzyaaaa151=1,ne
if(xyzzyaaaa151==ieskip)cycle
xyzzyaaab151=exvecs(4,xyzzyaaaa151)
if(xyzzyaaab151==0.d0)then
f(1,xyzzyaaaa151)=0.d0
df(1,xyzzyaaaa151)=0.d0
d2f(1,xyzzyaaaa151)=0.d0
else
xyzzyaaac151=log(xyzzyaaab151)
xyzzyaaad151=xyzzyaaab151*xyzzyaaac151
f(1,xyzzyaaaa151)=xyzzyaaab151*xyzzyaaad151
df(1,xyzzyaaaa151)=2.d0*xyzzyaaad151+xyzzyaaab151
d2f(1,xyzzyaaaa151)=2.d0*xyzzyaaac151+3.d0
endif
enddo
end subroutine xyzzyaaht1
subroutine xyzzyaahu1(ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa152
real(dp) xyzzyaaab152
do xyzzyaaaa152=1,ne
if(xyzzyaaaa152==ieskip)cycle
xyzzyaaab152=exvecs(4,xyzzyaaaa152)
if(xyzzyaaab152==0.d0)then
f(1,xyzzyaaaa152)=0.d0
else
f(1,xyzzyaaaa152)=1.d0/sqrt(xyzzyaaab152)
endif
enddo
end subroutine xyzzyaahu1
subroutine xyzzyaahv1(ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa153
real(dp) xyzzyaaab153,xyzzyaaac153
do xyzzyaaaa153=1,ne
if(xyzzyaaaa153==ieskip)cycle
xyzzyaaab153=exvecs(4,xyzzyaaaa153)
if(xyzzyaaab153==0.d0)then
f(1,xyzzyaaaa153)=0.d0
df(1,xyzzyaaaa153)=0.d0
else
xyzzyaaac153=1.d0/sqrt(xyzzyaaab153*xyzzyaaab153*xyzzyaaab153)
f(1,xyzzyaaaa153)=xyzzyaaac153*xyzzyaaab153
df(1,xyzzyaaaa153)=-0.5d0*xyzzyaaac153
endif
enddo
end subroutine xyzzyaahv1
subroutine xyzzyaahw1(ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa154
real(dp) xyzzyaaab154,xyzzyaaac154,xyzzyaaad154
do xyzzyaaaa154=1,ne
if(xyzzyaaaa154==ieskip)cycle
xyzzyaaab154=exvecs(4,xyzzyaaaa154)
if(xyzzyaaab154==0.d0)then
f(1,xyzzyaaaa154)=0.d0
df(1,xyzzyaaaa154)=0.d0
d2f(1,xyzzyaaaa154)=0.d0
else
xyzzyaaad154=1.d0/sqrt(xyzzyaaab154*xyzzyaaab154*xyzzyaaab154*xyzzyaaa&
&b154*xyzzyaaab154)
xyzzyaaac154=xyzzyaaad154*xyzzyaaab154
f(1,xyzzyaaaa154)=xyzzyaaac154*xyzzyaaab154
df(1,xyzzyaaaa154)=-0.5d0*xyzzyaaac154
d2f(1,xyzzyaaaa154)=0.75d0*xyzzyaaad154
endif
enddo
end subroutine xyzzyaahw1
subroutine xyzzyaahx1(n,iexp,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,iexp(n),ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa155,xyzzyaaab155,xyzzyaaac155
real(dp) xyzzyaaad155,xyzzyaaae155
do xyzzyaaaa155=1,ne
if(xyzzyaaaa155==ieskip)cycle
xyzzyaaad155=exvecs(4,xyzzyaaaa155)
if(xyzzyaaad155==0.d0)then
f(1:n,xyzzyaaaa155)=0.d0
else
xyzzyaaae155=sqrt(xyzzyaaad155)
do xyzzyaaab155=1,n
xyzzyaaac155=iexp(xyzzyaaab155)
select case(xyzzyaaac155)
case(-1)
f(xyzzyaaab155,xyzzyaaaa155)=1.d0/xyzzyaaae155
case(0)
f(xyzzyaaab155,xyzzyaaaa155)=log(xyzzyaaad155)
case(1)
f(xyzzyaaab155,xyzzyaaaa155)=xyzzyaaae155
case(2)
f(xyzzyaaab155,xyzzyaaaa155)=xyzzyaaad155
case(3)
f(xyzzyaaab155,xyzzyaaaa155)=xyzzyaaad155*xyzzyaaae155
case(4)
f(xyzzyaaab155,xyzzyaaaa155)=xyzzyaaad155*xyzzyaaad155
case default
if(mod(xyzzyaaac155,2)==0)then
f(xyzzyaaab155,xyzzyaaaa155)=xyzzyaaad155**(xyzzyaaac155/2)
else
f(xyzzyaaab155,xyzzyaaaa155)=xyzzyaaae155**xyzzyaaac155
endif
end select
enddo
endif
enddo
end subroutine xyzzyaahx1
subroutine xyzzyaahy1(n,iexp,ne,ieskip,exvecs,fstride,f,df)
implicit none
integer,intent(in) :: n,iexp(n),ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
integer xyzzyaaaa156,xyzzyaaab156,xyzzyaaac156
real(dp) xyzzyaaad156,xyzzyaaae156,xyzzyaaaf156,xyzzyaaag156,xyzzyaaah&
&156
do xyzzyaaaa156=1,ne
if(xyzzyaaaa156==ieskip)cycle
xyzzyaaad156=exvecs(4,xyzzyaaaa156)
if(xyzzyaaad156==0.d0)then
f(1:n,xyzzyaaaa156)=0.d0
df(1:n,xyzzyaaaa156)=0.d0
else
xyzzyaaae156=sqrt(xyzzyaaad156)
xyzzyaaaf156=1.d0/xyzzyaaae156
xyzzyaaag156=1.d0/xyzzyaaad156
do xyzzyaaab156=1,n
xyzzyaaac156=iexp(xyzzyaaab156)
select case(xyzzyaaac156)
case(-1)
f(xyzzyaaab156,xyzzyaaaa156)=xyzzyaaaf156
df(xyzzyaaab156,xyzzyaaaa156)=-0.5d0*xyzzyaaag156*xyzzyaaaf156
case(0)
f(xyzzyaaab156,xyzzyaaaa156)=log(xyzzyaaad156)
df(xyzzyaaab156,xyzzyaaaa156)=xyzzyaaag156
case(1)
f(xyzzyaaab156,xyzzyaaaa156)=xyzzyaaae156
df(xyzzyaaab156,xyzzyaaaa156)=0.5d0*xyzzyaaaf156
case(2)
f(xyzzyaaab156,xyzzyaaaa156)=xyzzyaaad156
df(xyzzyaaab156,xyzzyaaaa156)=1.d0
case(3)
f(xyzzyaaab156,xyzzyaaaa156)=xyzzyaaad156*xyzzyaaae156
df(xyzzyaaab156,xyzzyaaaa156)=1.5d0*xyzzyaaae156
case(4)
f(xyzzyaaab156,xyzzyaaaa156)=xyzzyaaad156*xyzzyaaad156
df(xyzzyaaab156,xyzzyaaaa156)=xyzzyaaad156+xyzzyaaad156
case default
if(mod(xyzzyaaac156,2)==0)then
xyzzyaaah156=xyzzyaaad156**(xyzzyaaac156/2-1)
else
xyzzyaaah156=xyzzyaaae156**(xyzzyaaac156-2)
endif
f(xyzzyaaab156,xyzzyaaaa156)=xyzzyaaah156*xyzzyaaad156
df(xyzzyaaab156,xyzzyaaaa156)=0.5d0*dble(xyzzyaaac156)*xyzzyaaah156
end select
enddo
endif
enddo
end subroutine xyzzyaahy1
subroutine xyzzyaahz1(n,iexp,ne,ieskip,exvecs,fstride,f,df,d2f)
implicit none
integer,intent(in) :: n,iexp(n),ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
integer xyzzyaaaa157,xyzzyaaab157,xyzzyaaac157
real(dp) xyzzyaaad157,xyzzyaaae157,xyzzyaaaf157,xyzzyaaag157,xyzzyaaah&
&157,xyzzyaaai157,xyzzyaaaj157
do xyzzyaaaa157=1,ne
if(xyzzyaaaa157==ieskip)cycle
xyzzyaaad157=exvecs(4,xyzzyaaaa157)
if(xyzzyaaad157==0.d0)then
f(1:n,xyzzyaaaa157)=0.d0
df(1:n,xyzzyaaaa157)=0.d0
else
xyzzyaaae157=sqrt(xyzzyaaad157)
xyzzyaaaf157=1.d0/xyzzyaaae157
xyzzyaaag157=1.d0/xyzzyaaad157
xyzzyaaah157=xyzzyaaaf157*xyzzyaaag157
do xyzzyaaab157=1,n
xyzzyaaac157=iexp(xyzzyaaab157)
select case(xyzzyaaac157)
case(-1)
f(xyzzyaaab157,xyzzyaaaa157)=xyzzyaaaf157
df(xyzzyaaab157,xyzzyaaaa157)=-0.5d0*xyzzyaaah157
d2f(xyzzyaaab157,xyzzyaaaa157)=0.75d0*xyzzyaaah157*xyzzyaaag157
case(0)
f(xyzzyaaab157,xyzzyaaaa157)=log(xyzzyaaad157)
df(xyzzyaaab157,xyzzyaaaa157)=xyzzyaaag157
d2f(xyzzyaaab157,xyzzyaaaa157)=-xyzzyaaag157*xyzzyaaag157
case(1)
f(xyzzyaaab157,xyzzyaaaa157)=xyzzyaaae157
df(xyzzyaaab157,xyzzyaaaa157)=0.5d0*xyzzyaaaf157
d2f(xyzzyaaab157,xyzzyaaaa157)=-0.25d0*xyzzyaaah157
case(2)
f(xyzzyaaab157,xyzzyaaaa157)=xyzzyaaad157
df(xyzzyaaab157,xyzzyaaaa157)=1.d0
d2f(xyzzyaaab157,xyzzyaaaa157)=0.d0
case(3)
f(xyzzyaaab157,xyzzyaaaa157)=xyzzyaaad157*xyzzyaaae157
df(xyzzyaaab157,xyzzyaaaa157)=1.5d0*xyzzyaaae157
d2f(xyzzyaaab157,xyzzyaaaa157)=0.75d0*xyzzyaaaf157
case(4)
f(xyzzyaaab157,xyzzyaaaa157)=xyzzyaaad157*xyzzyaaad157
df(xyzzyaaab157,xyzzyaaaa157)=xyzzyaaad157+xyzzyaaad157
d2f(xyzzyaaab157,xyzzyaaaa157)=2.d0
case default
if(mod(xyzzyaaac157,2)==0)then
xyzzyaaai157=xyzzyaaad157**(xyzzyaaac157/2-2)
else
xyzzyaaai157=xyzzyaaae157**(xyzzyaaac157-4)
endif
xyzzyaaaj157=xyzzyaaai157*xyzzyaaad157
f(xyzzyaaab157,xyzzyaaaa157)=xyzzyaaaj157*xyzzyaaad157
df(xyzzyaaab157,xyzzyaaaa157)=0.5d0*dble(xyzzyaaac157)*xyzzyaaaj157
d2f(xyzzyaaab157,xyzzyaaaa157)=0.25d0*dble((xyzzyaaac157-2)*xyzzyaaac1&
&57)*xyzzyaaai157
end select
enddo
endif
enddo
end subroutine xyzzyaahz1
subroutine xyzzyaaia1(theta,ne,ieskip,exvecs,fstride,f)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne),theta
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa158
real(dp) r,x,rootfunction,a0,a2
do xyzzyaaaa158=1,ne
if(xyzzyaaaa158==ieskip)cycle
r=exvecs(4,xyzzyaaaa158)
x=exvecs(1,xyzzyaaaa158)
rootfunction=sqrt(636.d0*cos(2.d0*theta)-437.d0-135.d0*cos(4.d0*theta)&
&)
a0=0.25d0*sqrt(2.d0+6.d0*cos(2.d0*theta)+rootfunction)
a2=a0*(rootfunction-2.d0-6.d0*cos(2.d0*theta))/(102.d0*sin(theta)*sin(&
&theta))
if(r==0.d0)then
f(1,xyzzyaaaa158)=0.d0
else
f(1,xyzzyaaaa158)=(a0+a2*(2.d0*x*x/(r*r)-1.d0))/sqrt(r)
endif
enddo
end subroutine xyzzyaaia1
subroutine xyzzyaaib1(theta,ne,ieskip,exvecs,fstride,f,gradf)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne),theta
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*)
integer xyzzyaaaa159
real(dp) r,one_r_3_2,x,rootfunction,a0,a2
do xyzzyaaaa159=1,ne
if(xyzzyaaaa159==ieskip)cycle
r=exvecs(4,xyzzyaaaa159)
x=exvecs(1,xyzzyaaaa159)
rootfunction=sqrt(636.d0*cos(2.d0*theta)-437.d0-135.d0*cos(4.d0*theta)&
&)
a0=0.25d0*sqrt(2.d0+6.d0*cos(2.d0*theta)+rootfunction)
a2=a0*(rootfunction-2.d0-6.d0*cos(2.d0*theta))/(102.d0*sin(theta)*sin(&
&theta))
if(r==0.d0)then
f(1,xyzzyaaaa159)=0.d0
gradf(1:3,1,xyzzyaaaa159)=0.d0
else
one_r_3_2=1.d0/sqrt(r*r*r)
f(1,xyzzyaaaa159)=(a0+a2*(2.d0*x*x/(r*r)-1.d0))*one_r_3_2*r
gradf(1,1,xyzzyaaaa159)=-0.5d0*(x/r)*(a0-4.d0*a2+5.d0*a2*(2.d0*x*x/(r*&
&r)-1.d0))*one_r_3_2
gradf(2,1,xyzzyaaaa159)=-0.5d0*(sqrt(1.d0-x*x/(r*r)))*(a0+4.d0*a2+5.d0&
&*a2*(2.d0*x*x/(r*r)-1.d0))*one_r_3_2
gradf(3,1,xyzzyaaaa159)=0.d0
endif
enddo
end subroutine xyzzyaaib1
subroutine xyzzyaaic1(theta,ne,ieskip,exvecs,fstride,f,gradf,lapf)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: exvecs(4,ne),theta
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),lapf(fstride&
&,*)
integer xyzzyaaaa160
real(dp) r,one_r_3_2,one_r_5_2,x,rootfunction,a0,a2
do xyzzyaaaa160=1,ne
if(xyzzyaaaa160==ieskip)cycle
r=exvecs(4,xyzzyaaaa160)
x=exvecs(1,xyzzyaaaa160)
rootfunction=sqrt(636.d0*cos(2.d0*theta)-437.d0-135.d0*cos(4.d0*theta)&
&)
a0=0.25d0*sqrt(2.d0+6.d0*cos(2.d0*theta)+rootfunction)
a2=a0*(rootfunction-2.d0-6.d0*cos(2.d0*theta))/(102.d0*sin(theta)*sin(&
&theta))
if(r==0.d0)then
f(1,xyzzyaaaa160)=0.d0
gradf(1:3,1,xyzzyaaaa160)=0.d0
lapf(1,xyzzyaaaa160)=0.d0
else
one_r_5_2=1.d0/sqrt(r*r*r*r*r)
one_r_3_2=one_r_5_2*r
f(1,xyzzyaaaa160)=(a0+a2*(2.d0*x*x/(r*r)-1.d0))*one_r_3_2*r
gradf(1,1,xyzzyaaaa160)=-0.5d0*(x/r)*(a0-4.d0*a2+5.d0*a2*(2.d0*x*x/(r*&
&r)-1.d0))*one_r_3_2
gradf(2,1,xyzzyaaaa160)=-0.5d0*(sqrt(1.d0-x*x/(r*r)))*(a0+4.d0*a2+5.d0&
&*a2*(2.d0*x*x/(r*r)-1.d0))*one_r_3_2
gradf(3,1,xyzzyaaaa160)=0.d0
lapf(1,xyzzyaaaa160)=0.25d0*(a0-15.d0*a2*(2.d0*x*x/(r*r)-1.d0))*one_r_&
&5_2
endif
enddo
end subroutine xyzzyaaic1
subroutine xyzzyaaid1(n,nvector,orthorhombic,inv_l,an,bn,an_an,ne,iesk&
&ip,exvecs,fstride,f)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,nvector,orthorhombic
real(dp),intent(in) :: exvecs(4,ne),inv_l,an(3,nvector),bn(3,nvector),&
&an_an(nvector,nvector)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa161,xyzzyaaab161
real(dp) xyzzyaaac161,xyzzyaaad161(nvector),xyzzyaaae161(nvector),xyzz&
&yaaaf161,xyzzyaaag161,xyzzyaaah161,xyzzyaaai161(3)
do xyzzyaaaa161=1,ne
if(xyzzyaaaa161==ieskip)cycle
if(exvecs(4,xyzzyaaaa161)==0.d0)then
f(1:n,xyzzyaaaa161)=0.d0
else
xyzzyaaai161(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa161)
do xyzzyaaab161=1,nvector
xyzzyaaac161=ddot(dimensionality,bn(1,xyzzyaaab161),1,xyzzyaaai161,1)
call xyzzyaajf1(xyzzyaaac161,xyzzyaaad161(xyzzyaaab161))
if(orthorhombic==0)call xyzzyaajg1(xyzzyaaac161,xyzzyaaae161(xyzzyaaab&
&161))
enddo
xyzzyaaag161=0.d0
do xyzzyaaab161=1,nvector
xyzzyaaag161=xyzzyaaag161+an_an(xyzzyaaab161,xyzzyaaab161)*xyzzyaaad16&
&1(xyzzyaaab161)*xyzzyaaad161(xyzzyaaab161)
if(orthorhombic==0.and.xyzzyaaab161<nvector)xyzzyaaag161=xyzzyaaag161+&
&2.d0*xyzzyaaae161(xyzzyaaab161)*ddot(nvector-xyzzyaaab161,an_an(xyzzy&
&aaab161,xyzzyaaab161+1),nvector,xyzzyaaae161(xyzzyaaab161+1),1)
enddo
xyzzyaaaf161=sqrt(xyzzyaaag161)
xyzzyaaah161=xyzzyaaaf161
f(1,xyzzyaaaa161)=xyzzyaaah161
do xyzzyaaab161=2,n
xyzzyaaah161=xyzzyaaah161*xyzzyaaaf161*inv_l
f(xyzzyaaab161,xyzzyaaaa161)=xyzzyaaah161
enddo
endif
enddo
end subroutine xyzzyaaid1
subroutine xyzzyaaie1(n,nvector,orthorhombic,inv_l,an,bn,an_an,ne,iesk&
&ip,exvecs,fstride,f,gradf)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,nvector,orthorhombic
real(dp),intent(in) :: exvecs(4,ne),inv_l,an(3,nvector),bn(3,nvector),&
&an_an(nvector,nvector)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*)
integer xyzzyaaaa162,xyzzyaaab162,xyzzyaaac162
real(dp) xyzzyaaad162,xyzzyaaae162(nvector),xyzzyaaaf162(nvector),xyzz&
&yaaag162,xyzzyaaah162,xyzzyaaai162,xyzzyaaaj162(nvector),xyzzyaaak162&
&(nvector),xyzzyaaal162(3),xyzzyaaam162(3)
do xyzzyaaaa162=1,ne
if(xyzzyaaaa162==ieskip)cycle
if(exvecs(4,xyzzyaaaa162)==0.d0)then
f(1:n,xyzzyaaaa162)=0.d0
gradf(1:3,1:n,xyzzyaaaa162)=0.d0
else
xyzzyaaam162(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa162)
do xyzzyaaab162=1,nvector
xyzzyaaad162=ddot(dimensionality,bn(1,xyzzyaaab162),1,xyzzyaaam162,1)
call xyzzyaajf1(xyzzyaaad162,xyzzyaaae162(xyzzyaaab162),xyzzyaaaj162(x&
&yzzyaaab162))
if(orthorhombic==0)call xyzzyaajg1(xyzzyaaad162,xyzzyaaaf162(xyzzyaaab&
&162),xyzzyaaak162(xyzzyaaab162))
enddo
xyzzyaaah162=0.d0
xyzzyaaal162=0.d0
do xyzzyaaab162=1,nvector
xyzzyaaah162=xyzzyaaah162+an_an(xyzzyaaab162,xyzzyaaab162)*xyzzyaaae16&
&2(xyzzyaaab162)*xyzzyaaae162(xyzzyaaab162)
xyzzyaaal162(1:dimensionality)=xyzzyaaal162(1:dimensionality)+an_an(xy&
&zzyaaab162,xyzzyaaab162)*xyzzyaaae162(xyzzyaaab162)*xyzzyaaaj162(xyzz&
&yaaab162)*bn(1:dimensionality,xyzzyaaab162)
if(orthorhombic==0)then
if(xyzzyaaab162<nvector)xyzzyaaah162=xyzzyaaah162+2.d0*xyzzyaaaf162(xy&
&zzyaaab162)*ddot(nvector-xyzzyaaab162,an_an(xyzzyaaab162,xyzzyaaab162&
&+1),nvector,xyzzyaaaf162(xyzzyaaab162+1),1)
do xyzzyaaac162=xyzzyaaab162+1,nvector
xyzzyaaal162(1:dimensionality)=xyzzyaaal162(1:dimensionality)+an_an(xy&
&zzyaaab162,xyzzyaaac162)*(xyzzyaaaf162(xyzzyaaab162)*xyzzyaaak162(xyz&
&zyaaac162)*bn(1:dimensionality,xyzzyaaac162)+xyzzyaaaf162(xyzzyaaac16&
&2)*xyzzyaaak162(xyzzyaaab162)*bn(1:dimensionality,xyzzyaaab162))
enddo
endif
enddo
xyzzyaaag162=sqrt(xyzzyaaah162)
xyzzyaaal162=xyzzyaaal162/xyzzyaaah162
xyzzyaaai162=xyzzyaaag162
f(1,xyzzyaaaa162)=xyzzyaaai162
gradf(1:dimensionality,1,xyzzyaaaa162)=xyzzyaaai162*xyzzyaaal162(1:dim&
&ensionality)
do xyzzyaaab162=2,n
xyzzyaaai162=xyzzyaaai162*xyzzyaaag162*inv_l
f(xyzzyaaab162,xyzzyaaaa162)=xyzzyaaai162
gradf(1:dimensionality,xyzzyaaab162,xyzzyaaaa162)=xyzzyaaai162*dble(xy&
&zzyaaab162)*xyzzyaaal162(1:dimensionality)
enddo
endif
enddo
end subroutine xyzzyaaie1
subroutine xyzzyaaif1(n,nvector,orthorhombic,inv_l,an,bn,an_an,bn_bn,n&
&e,ieskip,exvecs,fstride,f,gradf,lapf)
implicit none
integer,intent(in) :: n,ne,ieskip,fstride,nvector,orthorhombic
real(dp),intent(in) :: exvecs(4,ne),inv_l,an(3,nvector),bn(3,nvector),&
&an_an(nvector,nvector),bn_bn(nvector,nvector)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),lapf(fstride&
&,*)
integer xyzzyaaaa163,xyzzyaaab163,xyzzyaaac163
real(dp) xyzzyaaad163,xyzzyaaae163(nvector),xyzzyaaaf163(nvector),xyzz&
&yaaag163,xyzzyaaah163,xyzzyaaai163,xyzzyaaaj163,xyzzyaaak163(nvector)&
&,xyzzyaaal163(nvector),xyzzyaaam163(nvector),xyzzyaaan163(nvector),xy&
&zzyaaao163,xyzzyaaap163,xyzzyaaaq163(3),xyzzyaaar163(3)
do xyzzyaaaa163=1,ne
if(xyzzyaaaa163==ieskip)cycle
if(exvecs(4,xyzzyaaaa163)==0.d0)then
f(1:n,xyzzyaaaa163)=0.d0
gradf(1:3,1:n,xyzzyaaaa163)=0.d0
lapf(1:n,xyzzyaaaa163)=0.d0
else
xyzzyaaar163(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa163)
do xyzzyaaab163=1,nvector
xyzzyaaad163=ddot(dimensionality,bn(1,xyzzyaaab163),1,xyzzyaaar163,1)
call xyzzyaajf1(xyzzyaaad163,xyzzyaaae163(xyzzyaaab163),xyzzyaaak163(x&
&yzzyaaab163),xyzzyaaam163(xyzzyaaab163))
if(orthorhombic==0)call xyzzyaajg1(xyzzyaaad163,xyzzyaaaf163(xyzzyaaab&
&163),xyzzyaaal163(xyzzyaaab163),xyzzyaaan163(xyzzyaaab163))
enddo
xyzzyaaah163=0.d0
xyzzyaaaq163=0.d0
xyzzyaaap163=0.d0
do xyzzyaaab163=1,nvector
xyzzyaaah163=xyzzyaaah163+an_an(xyzzyaaab163,xyzzyaaab163)*xyzzyaaae16&
&3(xyzzyaaab163)*xyzzyaaae163(xyzzyaaab163)
xyzzyaaaq163(1:dimensionality)=xyzzyaaaq163(1:dimensionality)+an_an(xy&
&zzyaaab163,xyzzyaaab163)*xyzzyaaae163(xyzzyaaab163)*xyzzyaaak163(xyzz&
&yaaab163)*bn(1:dimensionality,xyzzyaaab163)
xyzzyaaap163=xyzzyaaap163+an_an(xyzzyaaab163,xyzzyaaab163)*(xyzzyaaak1&
&63(xyzzyaaab163)*xyzzyaaak163(xyzzyaaab163)+xyzzyaaae163(xyzzyaaab163&
&)*xyzzyaaam163(xyzzyaaab163))*bn_bn(xyzzyaaab163,xyzzyaaab163)
if(orthorhombic==0)then
if(xyzzyaaab163<nvector)xyzzyaaah163=xyzzyaaah163+2.d0*xyzzyaaaf163(xy&
&zzyaaab163)*ddot(nvector-xyzzyaaab163,an_an(xyzzyaaab163,xyzzyaaab163&
&+1),nvector,xyzzyaaaf163(xyzzyaaab163+1),1)
do xyzzyaaac163=xyzzyaaab163+1,nvector
xyzzyaaaq163(1:dimensionality)=xyzzyaaaq163(1:dimensionality)+an_an(xy&
&zzyaaab163,xyzzyaaac163)*(xyzzyaaaf163(xyzzyaaab163)*xyzzyaaal163(xyz&
&zyaaac163)*bn(1:dimensionality,xyzzyaaac163)+xyzzyaaaf163(xyzzyaaac16&
&3)*xyzzyaaal163(xyzzyaaab163)*bn(1:dimensionality,xyzzyaaab163))
xyzzyaaap163=xyzzyaaap163+an_an(xyzzyaaab163,xyzzyaaac163)*(2.d0*xyzzy&
&aaal163(xyzzyaaab163)*xyzzyaaal163(xyzzyaaac163)*bn_bn(xyzzyaaab163,x&
&yzzyaaac163)+xyzzyaaaf163(xyzzyaaab163)*xyzzyaaan163(xyzzyaaac163)*bn&
&_bn(xyzzyaaac163,xyzzyaaac163)+xyzzyaaaf163(xyzzyaaac163)*xyzzyaaan16&
&3(xyzzyaaab163)*bn_bn(xyzzyaaab163,xyzzyaaab163))
enddo
endif
enddo
xyzzyaaag163=sqrt(xyzzyaaah163)
xyzzyaaai163=1.d0/xyzzyaaah163
xyzzyaaaq163=xyzzyaaaq163*xyzzyaaai163
xyzzyaaap163=xyzzyaaap163*xyzzyaaai163
xyzzyaaao163=ddot(dimensionality,xyzzyaaaq163,1,xyzzyaaaq163,1)
xyzzyaaaj163=xyzzyaaag163
f(1,xyzzyaaaa163)=xyzzyaaaj163
gradf(1:dimensionality,1,xyzzyaaaa163)=xyzzyaaaj163*xyzzyaaaq163(1:dim&
&ensionality)
lapf(1,xyzzyaaaa163)=xyzzyaaaj163*(-xyzzyaaao163+xyzzyaaap163)
do xyzzyaaab163=2,n
xyzzyaaaj163=xyzzyaaaj163*xyzzyaaag163*inv_l
f(xyzzyaaab163,xyzzyaaaa163)=xyzzyaaaj163
gradf(1:dimensionality,xyzzyaaab163,xyzzyaaaa163)=xyzzyaaaj163*dble(xy&
&zzyaaab163)*xyzzyaaaq163(1:dimensionality)
lapf(xyzzyaaab163,xyzzyaaaa163)=xyzzyaaaj163*dble(xyzzyaaab163)*(dble(&
&xyzzyaaab163-2)*xyzzyaaao163+xyzzyaaap163)
enddo
endif
enddo
end subroutine xyzzyaaif1
subroutine xyzzyaaig1(ctrunc,p,ne,ieskip,exvecs,fstride,f,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride,ctrunc
real(dp),intent(in) :: p(4),exvecs(4,ne)
logical,intent(inout) :: nzcut(ne)
real(dp),intent(inout) :: f(fstride,*)
integer xyzzyaaaa164
real(dp) xyzzyaaab164,xyzzyaaac164,xyzzyaaad164,xyzzyaaae164
xyzzyaaad164=p(1)
xyzzyaaae164=p(2)
select case(ctrunc)
case(0)
do xyzzyaaaa164=1,ne
if(xyzzyaaaa164==ieskip)cycle
if(exvecs(4,xyzzyaaaa164)>=xyzzyaaad164)then
f(1,xyzzyaaaa164)=0.d0
nzcut(xyzzyaaaa164)=.false.
else
f(1,xyzzyaaaa164)=1.d0
nzcut(xyzzyaaaa164)=.true.
endif
enddo
case(1)
do xyzzyaaaa164=1,ne
if(xyzzyaaaa164==ieskip)cycle
xyzzyaaab164=exvecs(4,xyzzyaaaa164)
if(xyzzyaaab164>=xyzzyaaad164)then
f(1,xyzzyaaaa164)=0.d0
nzcut(xyzzyaaaa164)=.false.
else
f(1,xyzzyaaaa164)=1.d0-xyzzyaaab164*xyzzyaaae164
nzcut(xyzzyaaaa164)=.true.
endif
enddo
case(2)
do xyzzyaaaa164=1,ne
if(xyzzyaaaa164==ieskip)cycle
xyzzyaaab164=exvecs(4,xyzzyaaaa164)
if(xyzzyaaab164>=xyzzyaaad164)then
f(1,xyzzyaaaa164)=0.d0
nzcut(xyzzyaaaa164)=.false.
else
xyzzyaaac164=1.d0-xyzzyaaab164*xyzzyaaae164
f(1,xyzzyaaaa164)=xyzzyaaac164*xyzzyaaac164
nzcut(xyzzyaaaa164)=.true.
endif
enddo
case(3)
do xyzzyaaaa164=1,ne
if(xyzzyaaaa164==ieskip)cycle
xyzzyaaab164=exvecs(4,xyzzyaaaa164)
if(xyzzyaaab164>=xyzzyaaad164)then
f(1,xyzzyaaaa164)=0.d0
nzcut(xyzzyaaaa164)=.false.
else
xyzzyaaac164=1.d0-xyzzyaaab164*xyzzyaaae164
f(1,xyzzyaaaa164)=xyzzyaaac164*xyzzyaaac164*xyzzyaaac164
nzcut(xyzzyaaaa164)=.true.
endif
enddo
case default
do xyzzyaaaa164=1,ne
if(xyzzyaaaa164==ieskip)cycle
xyzzyaaab164=exvecs(4,xyzzyaaaa164)
if(xyzzyaaab164>=xyzzyaaad164)then
f(1,xyzzyaaaa164)=0.d0
nzcut(xyzzyaaaa164)=.false.
else
xyzzyaaac164=1.d0-xyzzyaaab164*xyzzyaaae164
f(1,xyzzyaaaa164)=xyzzyaaac164**ctrunc
nzcut(xyzzyaaaa164)=.true.
endif
enddo
end select
end subroutine xyzzyaaig1
subroutine xyzzyaaih1(ctrunc,p,ne,ieskip,exvecs,fstride,f,df,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride,ctrunc
real(dp),intent(in) :: p(4),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa165
real(dp) xyzzyaaab165,xyzzyaaac165,xyzzyaaad165,xyzzyaaae165,xyzzyaaaf&
&165,xyzzyaaag165
xyzzyaaac165=p(1)
xyzzyaaad165=p(2)
xyzzyaaag165=p(3)
select case(ctrunc)
case(0)
do xyzzyaaaa165=1,ne
if(xyzzyaaaa165==ieskip)cycle
if(exvecs(4,xyzzyaaaa165)>=xyzzyaaac165)then
f(1,xyzzyaaaa165)=0.d0
nzcut(xyzzyaaaa165)=.false.
else
f(1,xyzzyaaaa165)=1.d0
nzcut(xyzzyaaaa165)=.true.
endif
df(1,xyzzyaaaa165)=0.d0
enddo
case(1)
do xyzzyaaaa165=1,ne
if(xyzzyaaaa165==ieskip)cycle
xyzzyaaab165=exvecs(4,xyzzyaaaa165)
if(xyzzyaaab165>=xyzzyaaac165)then
f(1,xyzzyaaaa165)=0.d0
df(1,xyzzyaaaa165)=0.d0
nzcut(xyzzyaaaa165)=.false.
else
f(1,xyzzyaaaa165)=1.d0-xyzzyaaab165*xyzzyaaad165
df(1,xyzzyaaaa165)=xyzzyaaag165
nzcut(xyzzyaaaa165)=.true.
endif
enddo
case(2)
do xyzzyaaaa165=1,ne
if(xyzzyaaaa165==ieskip)cycle
xyzzyaaab165=exvecs(4,xyzzyaaaa165)
if(xyzzyaaab165>=xyzzyaaac165)then
f(1,xyzzyaaaa165)=0.d0
df(1,xyzzyaaaa165)=0.d0
nzcut(xyzzyaaaa165)=.false.
else
xyzzyaaae165=1.d0-xyzzyaaab165*xyzzyaaad165
f(1,xyzzyaaaa165)=xyzzyaaae165*xyzzyaaae165
df(1,xyzzyaaaa165)=xyzzyaaag165*xyzzyaaae165
nzcut(xyzzyaaaa165)=.true.
endif
enddo
case(3)
do xyzzyaaaa165=1,ne
if(xyzzyaaaa165==ieskip)cycle
xyzzyaaab165=exvecs(4,xyzzyaaaa165)
if(xyzzyaaab165>=xyzzyaaac165)then
f(1,xyzzyaaaa165)=0.d0
df(1,xyzzyaaaa165)=0.d0
nzcut(xyzzyaaaa165)=.false.
else
xyzzyaaae165=1.d0-xyzzyaaab165*xyzzyaaad165
xyzzyaaaf165=xyzzyaaae165*xyzzyaaae165
f(1,xyzzyaaaa165)=xyzzyaaaf165*xyzzyaaae165
df(1,xyzzyaaaa165)=xyzzyaaag165*xyzzyaaaf165
nzcut(xyzzyaaaa165)=.true.
endif
enddo
case default
do xyzzyaaaa165=1,ne
if(xyzzyaaaa165==ieskip)cycle
xyzzyaaab165=exvecs(4,xyzzyaaaa165)
if(xyzzyaaab165>=xyzzyaaac165)then
f(1,xyzzyaaaa165)=0.d0
df(1,xyzzyaaaa165)=0.d0
nzcut(xyzzyaaaa165)=.false.
else
xyzzyaaae165=1.d0-xyzzyaaab165*xyzzyaaad165
xyzzyaaaf165=xyzzyaaae165**(ctrunc-1)
f(1,xyzzyaaaa165)=xyzzyaaaf165*xyzzyaaae165
df(1,xyzzyaaaa165)=xyzzyaaag165*xyzzyaaaf165
nzcut(xyzzyaaaa165)=.true.
endif
enddo
end select
end subroutine xyzzyaaih1
subroutine xyzzyaaii1(ctrunc,p,ne,ieskip,exvecs,fstride,f,df,d2f,nzcut&
&)
implicit none
integer,intent(in) :: ne,ieskip,fstride,ctrunc
real(dp),intent(in) :: p(4),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa166
real(dp) xyzzyaaab166,xyzzyaaac166,xyzzyaaad166,xyzzyaaae166,xyzzyaaaf&
&166,xyzzyaaag166,xyzzyaaah166,xyzzyaaai166
xyzzyaaac166=p(1)
xyzzyaaad166=p(2)
xyzzyaaah166=p(3)
xyzzyaaai166=p(4)
select case(ctrunc)
case(0)
do xyzzyaaaa166=1,ne
if(xyzzyaaaa166==ieskip)cycle
if(exvecs(4,xyzzyaaaa166)>=xyzzyaaac166)then
f(1,xyzzyaaaa166)=0.d0
nzcut(xyzzyaaaa166)=.false.
else
f(1,xyzzyaaaa166)=1.d0
nzcut(xyzzyaaaa166)=.true.
endif
df(1,xyzzyaaaa166)=0.d0
d2f(1,xyzzyaaaa166)=0.d0
enddo
case(1)
do xyzzyaaaa166=1,ne
if(xyzzyaaaa166==ieskip)cycle
xyzzyaaab166=exvecs(4,xyzzyaaaa166)
if(xyzzyaaab166>=xyzzyaaac166)then
f(1,xyzzyaaaa166)=0.d0
df(1,xyzzyaaaa166)=0.d0
nzcut(xyzzyaaaa166)=.false.
else
f(1,xyzzyaaaa166)=1.d0-xyzzyaaab166*xyzzyaaad166
df(1,xyzzyaaaa166)=xyzzyaaah166
nzcut(xyzzyaaaa166)=.true.
endif
d2f(1,xyzzyaaaa166)=0.d0
enddo
case(2)
do xyzzyaaaa166=1,ne
if(xyzzyaaaa166==ieskip)cycle
xyzzyaaab166=exvecs(4,xyzzyaaaa166)
if(xyzzyaaab166>=xyzzyaaac166)then
f(1,xyzzyaaaa166)=0.d0
df(1,xyzzyaaaa166)=0.d0
d2f(1,xyzzyaaaa166)=0.d0
nzcut(xyzzyaaaa166)=.false.
else
xyzzyaaae166=1.d0-xyzzyaaab166*xyzzyaaad166
f(1,xyzzyaaaa166)=xyzzyaaae166*xyzzyaaae166
df(1,xyzzyaaaa166)=xyzzyaaah166*xyzzyaaae166
d2f(1,xyzzyaaaa166)=xyzzyaaai166
nzcut(xyzzyaaaa166)=.true.
endif
enddo
case(3)
do xyzzyaaaa166=1,ne
if(xyzzyaaaa166==ieskip)cycle
xyzzyaaab166=exvecs(4,xyzzyaaaa166)
if(xyzzyaaab166>=xyzzyaaac166)then
f(1,xyzzyaaaa166)=0.d0
df(1,xyzzyaaaa166)=0.d0
d2f(1,xyzzyaaaa166)=0.d0
nzcut(xyzzyaaaa166)=.false.
else
xyzzyaaae166=1.d0-xyzzyaaab166*xyzzyaaad166
xyzzyaaaf166=xyzzyaaae166*xyzzyaaae166
f(1,xyzzyaaaa166)=xyzzyaaaf166*xyzzyaaae166
df(1,xyzzyaaaa166)=xyzzyaaah166*xyzzyaaaf166
d2f(1,xyzzyaaaa166)=xyzzyaaai166*xyzzyaaae166
nzcut(xyzzyaaaa166)=.true.
endif
enddo
case default
do xyzzyaaaa166=1,ne
if(xyzzyaaaa166==ieskip)cycle
xyzzyaaab166=exvecs(4,xyzzyaaaa166)
if(xyzzyaaab166>=xyzzyaaac166)then
f(1,xyzzyaaaa166)=0.d0
df(1,xyzzyaaaa166)=0.d0
d2f(1,xyzzyaaaa166)=0.d0
nzcut(xyzzyaaaa166)=.false.
else
xyzzyaaae166=1.d0-xyzzyaaab166*xyzzyaaad166
xyzzyaaag166=xyzzyaaae166**(ctrunc-2)
xyzzyaaaf166=xyzzyaaag166*xyzzyaaae166
f(1,xyzzyaaaa166)=xyzzyaaaf166*xyzzyaaae166
df(1,xyzzyaaaa166)=xyzzyaaah166*xyzzyaaaf166
d2f(1,xyzzyaaaa166)=xyzzyaaai166*xyzzyaaag166
nzcut(xyzzyaaaa166)=.true.
endif
enddo
end select
end subroutine xyzzyaaii1
subroutine xyzzyaaij1(ctrunc,lcut,ne,ieskip,exvecs,fstride,f,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride,ctrunc
real(dp),intent(in) :: lcut,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa167
real(dp) xyzzyaaab167,xyzzyaaac167
select case(ctrunc)
case(0)
do xyzzyaaaa167=1,ne
if(xyzzyaaaa167==ieskip)cycle
if(exvecs(4,xyzzyaaaa167)>=lcut)then
f(1,xyzzyaaaa167)=0.d0
nzcut(xyzzyaaaa167)=.false.
else
f(1,xyzzyaaaa167)=1.d0
nzcut(xyzzyaaaa167)=.true.
endif
enddo
case(1)
do xyzzyaaaa167=1,ne
if(xyzzyaaaa167==ieskip)cycle
xyzzyaaab167=exvecs(4,xyzzyaaaa167)
if(xyzzyaaab167>=lcut)then
f(1,xyzzyaaaa167)=0.d0
nzcut(xyzzyaaaa167)=.false.
else
f(1,xyzzyaaaa167)=xyzzyaaab167-lcut
nzcut(xyzzyaaaa167)=.true.
endif
enddo
case(2)
do xyzzyaaaa167=1,ne
if(xyzzyaaaa167==ieskip)cycle
xyzzyaaab167=exvecs(4,xyzzyaaaa167)
if(xyzzyaaab167>=lcut)then
f(1,xyzzyaaaa167)=0.d0
nzcut(xyzzyaaaa167)=.false.
else
xyzzyaaac167=xyzzyaaab167-lcut
f(1,xyzzyaaaa167)=xyzzyaaac167*xyzzyaaac167
nzcut(xyzzyaaaa167)=.true.
endif
enddo
case(3)
do xyzzyaaaa167=1,ne
if(xyzzyaaaa167==ieskip)cycle
xyzzyaaab167=exvecs(4,xyzzyaaaa167)
if(xyzzyaaab167>=lcut)then
f(1,xyzzyaaaa167)=0.d0
nzcut(xyzzyaaaa167)=.false.
else
xyzzyaaac167=xyzzyaaab167-lcut
f(1,xyzzyaaaa167)=xyzzyaaac167*xyzzyaaac167*xyzzyaaac167
nzcut(xyzzyaaaa167)=.true.
endif
enddo
case default
do xyzzyaaaa167=1,ne
if(xyzzyaaaa167==ieskip)cycle
xyzzyaaab167=exvecs(4,xyzzyaaaa167)
if(xyzzyaaab167>=lcut)then
f(1,xyzzyaaaa167)=0.d0
nzcut(xyzzyaaaa167)=.false.
else
xyzzyaaac167=xyzzyaaab167-lcut
f(1,xyzzyaaaa167)=xyzzyaaac167**ctrunc
nzcut(xyzzyaaaa167)=.true.
endif
enddo
end select
end subroutine xyzzyaaij1
subroutine xyzzyaaik1(ctrunc,lcut,ne,ieskip,exvecs,fstride,f,df,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride,ctrunc
real(dp),intent(in) :: lcut,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa168
real(dp) xyzzyaaab168,xyzzyaaac168,xyzzyaaad168
select case(ctrunc)
case(0)
do xyzzyaaaa168=1,ne
if(xyzzyaaaa168==ieskip)cycle
if(exvecs(4,xyzzyaaaa168)>=lcut)then
f(1,xyzzyaaaa168)=0.d0
nzcut(xyzzyaaaa168)=.false.
else
f(1,xyzzyaaaa168)=1.d0
nzcut(xyzzyaaaa168)=.true.
endif
df(1,xyzzyaaaa168)=0.d0
enddo
case(1)
do xyzzyaaaa168=1,ne
if(xyzzyaaaa168==ieskip)cycle
xyzzyaaab168=exvecs(4,xyzzyaaaa168)
if(xyzzyaaab168>=lcut)then
f(1,xyzzyaaaa168)=0.d0
df(1,xyzzyaaaa168)=0.d0
nzcut(xyzzyaaaa168)=.false.
else
f(1,xyzzyaaaa168)=xyzzyaaab168-lcut
df(1,xyzzyaaaa168)=1.d0
nzcut(xyzzyaaaa168)=.true.
endif
enddo
case(2)
do xyzzyaaaa168=1,ne
if(xyzzyaaaa168==ieskip)cycle
xyzzyaaab168=exvecs(4,xyzzyaaaa168)
if(xyzzyaaab168>=lcut)then
f(1,xyzzyaaaa168)=0.d0
df(1,xyzzyaaaa168)=0.d0
nzcut(xyzzyaaaa168)=.false.
else
xyzzyaaac168=xyzzyaaab168-lcut
f(1,xyzzyaaaa168)=xyzzyaaac168*xyzzyaaac168
df(1,xyzzyaaaa168)=2*xyzzyaaac168
nzcut(xyzzyaaaa168)=.true.
endif
enddo
case(3)
do xyzzyaaaa168=1,ne
if(xyzzyaaaa168==ieskip)cycle
xyzzyaaab168=exvecs(4,xyzzyaaaa168)
if(xyzzyaaab168>=lcut)then
f(1,xyzzyaaaa168)=0.d0
df(1,xyzzyaaaa168)=0.d0
nzcut(xyzzyaaaa168)=.false.
else
xyzzyaaac168=xyzzyaaab168-lcut
xyzzyaaad168=xyzzyaaac168*xyzzyaaac168
f(1,xyzzyaaaa168)=xyzzyaaad168*xyzzyaaac168
df(1,xyzzyaaaa168)=3*xyzzyaaad168
nzcut(xyzzyaaaa168)=.true.
endif
enddo
case default
do xyzzyaaaa168=1,ne
if(xyzzyaaaa168==ieskip)cycle
xyzzyaaab168=exvecs(4,xyzzyaaaa168)
if(xyzzyaaab168>=lcut)then
f(1,xyzzyaaaa168)=0.d0
df(1,xyzzyaaaa168)=0.d0
nzcut(xyzzyaaaa168)=.false.
else
xyzzyaaac168=xyzzyaaab168-lcut
xyzzyaaad168=xyzzyaaac168**(ctrunc-1)
f(1,xyzzyaaaa168)=xyzzyaaad168*xyzzyaaac168
df(1,xyzzyaaaa168)=ctrunc*xyzzyaaad168
nzcut(xyzzyaaaa168)=.true.
endif
enddo
end select
end subroutine xyzzyaaik1
subroutine xyzzyaail1(ctrunc,lcut,ne,ieskip,exvecs,fstride,f,df,d2f,nz&
&cut)
implicit none
integer,intent(in) :: ne,ieskip,fstride,ctrunc
real(dp),intent(in) :: lcut,exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa169
real(dp) xyzzyaaab169,xyzzyaaac169,xyzzyaaad169,xyzzyaaae169
select case(ctrunc)
case(0)
do xyzzyaaaa169=1,ne
if(xyzzyaaaa169==ieskip)cycle
if(exvecs(4,xyzzyaaaa169)>=lcut)then
f(1,xyzzyaaaa169)=0.d0
nzcut(xyzzyaaaa169)=.false.
else
f(1,xyzzyaaaa169)=1.d0
nzcut(xyzzyaaaa169)=.true.
endif
df(1,xyzzyaaaa169)=0.d0
d2f(1,xyzzyaaaa169)=0.d0
enddo
case(1)
do xyzzyaaaa169=1,ne
if(xyzzyaaaa169==ieskip)cycle
xyzzyaaab169=exvecs(4,xyzzyaaaa169)
if(xyzzyaaab169>=lcut)then
f(1,xyzzyaaaa169)=0.d0
df(1,xyzzyaaaa169)=0.d0
nzcut(xyzzyaaaa169)=.false.
else
f(1,xyzzyaaaa169)=xyzzyaaab169-lcut
df(1,xyzzyaaaa169)=1.d0
nzcut(xyzzyaaaa169)=.true.
endif
d2f(1,xyzzyaaaa169)=0.d0
enddo
case(2)
do xyzzyaaaa169=1,ne
if(xyzzyaaaa169==ieskip)cycle
xyzzyaaab169=exvecs(4,xyzzyaaaa169)
if(xyzzyaaab169>=lcut)then
f(1,xyzzyaaaa169)=0.d0
df(1,xyzzyaaaa169)=0.d0
d2f(1,xyzzyaaaa169)=0.d0
nzcut(xyzzyaaaa169)=.false.
else
xyzzyaaac169=xyzzyaaab169-lcut
f(1,xyzzyaaaa169)=xyzzyaaac169*xyzzyaaac169
df(1,xyzzyaaaa169)=2*xyzzyaaac169
d2f(1,xyzzyaaaa169)=2.d0
nzcut(xyzzyaaaa169)=.true.
endif
enddo
case(3)
do xyzzyaaaa169=1,ne
if(xyzzyaaaa169==ieskip)cycle
xyzzyaaab169=exvecs(4,xyzzyaaaa169)
if(xyzzyaaab169>=lcut)then
f(1,xyzzyaaaa169)=0.d0
df(1,xyzzyaaaa169)=0.d0
d2f(1,xyzzyaaaa169)=0.d0
nzcut(xyzzyaaaa169)=.false.
else
xyzzyaaac169=xyzzyaaab169-lcut
xyzzyaaad169=xyzzyaaac169*xyzzyaaac169
f(1,xyzzyaaaa169)=xyzzyaaad169*xyzzyaaac169
df(1,xyzzyaaaa169)=3*xyzzyaaad169
d2f(1,xyzzyaaaa169)=6*xyzzyaaac169
nzcut(xyzzyaaaa169)=.true.
endif
enddo
case default
do xyzzyaaaa169=1,ne
if(xyzzyaaaa169==ieskip)cycle
xyzzyaaab169=exvecs(4,xyzzyaaaa169)
if(xyzzyaaab169>=lcut)then
f(1,xyzzyaaaa169)=0.d0
df(1,xyzzyaaaa169)=0.d0
d2f(1,xyzzyaaaa169)=0.d0
nzcut(xyzzyaaaa169)=.false.
else
xyzzyaaac169=xyzzyaaab169-lcut
xyzzyaaae169=xyzzyaaac169**(ctrunc-2)
xyzzyaaad169=xyzzyaaae169*xyzzyaaac169
f(1,xyzzyaaaa169)=xyzzyaaad169*xyzzyaaac169
df(1,xyzzyaaaa169)=ctrunc*xyzzyaaad169
d2f(1,xyzzyaaaa169)=ctrunc*(ctrunc-1)*xyzzyaaae169
nzcut(xyzzyaaaa169)=.true.
endif
enddo
end select
end subroutine xyzzyaail1
subroutine xyzzyaaim1(p,ne,ieskip,exvecs,fstride,f,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa170
real(dp) xyzzyaaab170,xyzzyaaac170,xyzzyaaad170,xyzzyaaae170
xyzzyaaad170=p(1)
xyzzyaaae170=p(3)
do xyzzyaaaa170=1,ne
if(xyzzyaaaa170==ieskip)cycle
xyzzyaaab170=exvecs(4,ne)
if(xyzzyaaab170>=xyzzyaaae170)then
f(1,xyzzyaaaa170)=0.d0
nzcut(xyzzyaaaa170)=.false.
else
xyzzyaaac170=xyzzyaaab170*xyzzyaaad170
f(1,xyzzyaaaa170)=exp(-xyzzyaaac170*xyzzyaaac170)
nzcut(xyzzyaaaa170)=.true.
endif
enddo
end subroutine xyzzyaaim1
subroutine xyzzyaain1(p,ne,ieskip,exvecs,fstride,f,df,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa171
real(dp) xyzzyaaab171,xyzzyaaac171,xyzzyaaad171,xyzzyaaae171,xyzzyaaaf&
&171,xyzzyaaag171
xyzzyaaad171=p(1)
xyzzyaaae171=p(3)
do xyzzyaaaa171=1,ne
if(xyzzyaaaa171==ieskip)cycle
xyzzyaaab171=exvecs(4,ne)
if(xyzzyaaab171>=xyzzyaaae171)then
f(1,xyzzyaaaa171)=0.d0
df(1,xyzzyaaaa171)=0.d0
nzcut(xyzzyaaaa171)=.false.
else
xyzzyaaac171=xyzzyaaab171*xyzzyaaad171
xyzzyaaaf171=exp(-xyzzyaaac171*xyzzyaaac171)
xyzzyaaag171=-xyzzyaaad171*xyzzyaaac171*(xyzzyaaaf171+xyzzyaaaf171)
f(1,xyzzyaaaa171)=xyzzyaaaf171
df(1,xyzzyaaaa171)=xyzzyaaag171
nzcut(xyzzyaaaa171)=.true.
endif
enddo
end subroutine xyzzyaain1
subroutine xyzzyaaio1(p,ne,ieskip,exvecs,fstride,f,df,d2f,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa172
real(dp) xyzzyaaab172,xyzzyaaac172,xyzzyaaad172,xyzzyaaae172,xyzzyaaaf&
&172,xyzzyaaag172,xyzzyaaah172,xyzzyaaai172,xyzzyaaaj172
xyzzyaaae172=p(1)
xyzzyaaaf172=p(2)
xyzzyaaag172=p(3)
do xyzzyaaaa172=1,ne
if(xyzzyaaaa172==ieskip)cycle
xyzzyaaab172=exvecs(4,ne)
if(xyzzyaaab172>=xyzzyaaag172)then
f(1,xyzzyaaaa172)=0.d0
df(1,xyzzyaaaa172)=0.d0
d2f(1,xyzzyaaaa172)=0.d0
nzcut(xyzzyaaaa172)=.false.
else
xyzzyaaac172=xyzzyaaab172*xyzzyaaae172
xyzzyaaad172=xyzzyaaac172*xyzzyaaac172
xyzzyaaah172=exp(-xyzzyaaad172)
xyzzyaaai172=-xyzzyaaae172*xyzzyaaac172*(xyzzyaaah172+xyzzyaaah172)
xyzzyaaaj172=xyzzyaaaf172*(xyzzyaaad172+xyzzyaaad172-1.d0)*(xyzzyaaah1&
&72+xyzzyaaah172)
f(1,xyzzyaaaa172)=xyzzyaaah172
df(1,xyzzyaaaa172)=xyzzyaaai172
d2f(1,xyzzyaaaa172)=xyzzyaaaj172
nzcut(xyzzyaaaa172)=.true.
endif
enddo
end subroutine xyzzyaaio1
subroutine xyzzyaaip1(const_int,const_dble,p,ne,ieskip,ion1,exvecs,fst&
&ride,f,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,ion1,fstride,const_int(*)
real(dp),intent(in) :: const_dble(*),p(4),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa173,xyzzyaaab173,xyzzyaaac173,xyzzyaaad173,xyzzyaaae1&
&73,xyzzyaaaf173,xyzzyaaag173,xyzzyaaah173,xyzzyaaai173,xyzzyaaaj173,x&
&yzzyaaak173,xyzzyaaal173,xyzzyaaam173,xyzzyaaan173
real(dp) xyzzyaaao173(3),xyzzyaaap173(3),xyzzyaaaq173(3),xyzzyaaar173,&
&xyzzyaaas173,xyzzyaaat173,xyzzyaaau173,xyzzyaaav173,xyzzyaaaw173
call xyzzyaaig1(const_int(1),p(1),ne,ieskip,exvecs(1,1),fstride,f(1,1)&
&,nzcut(1))
xyzzyaaah173=dimensionality*dimensionality
xyzzyaaai173=const_int(2)
xyzzyaaab173=const_int(3)
xyzzyaaal173=const_int(4)
xyzzyaaak173=4+xyzzyaaab173*dimensionality
do xyzzyaaaa173=1,ne
if(xyzzyaaaa173==ieskip)cycle
if(.not.nzcut(xyzzyaaaa173))cycle
if(xyzzyaaaa173>ieskip)then
xyzzyaaao173(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa173)
else
xyzzyaaao173(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa173)
endif
xyzzyaaau173=exvecs(4,xyzzyaaaa173)
if(xyzzyaaau173==0.d0)then
xyzzyaaav173=0.d0
else
xyzzyaaav173=xyzzyaaau173**(-xyzzyaaai173)
endif
if(xyzzyaaal173/=0)then
if(xyzzyaaal173<0.or.ion1==0)then
xyzzyaaaj173=1
xyzzyaaam173=xyzzyaaab173
else
xyzzyaaaj173=const_int(xyzzyaaak173+ion1+xyzzyaaaa173-1)
xyzzyaaam173=xyzzyaaab173+(xyzzyaaaj173-1)*(xyzzyaaah173+xyzzyaaah173)
endif
xyzzyaaae173=xyzzyaaam173+1
do xyzzyaaad173=1,dimensionality
xyzzyaaap173(xyzzyaaad173)=ddot(dimensionality,xyzzyaaao173(1),1,const&
&_dble(xyzzyaaae173),1)
xyzzyaaae173=xyzzyaaae173+dimensionality
enddo
else
xyzzyaaap173(1:dimensionality)=xyzzyaaao173(1:dimensionality)
endif
xyzzyaaas173=0.d0
xyzzyaaaf173=4
xyzzyaaag173=0
do xyzzyaaac173=1,xyzzyaaab173
xyzzyaaag173=xyzzyaaag173+1
xyzzyaaat173=const_dble(xyzzyaaag173)
do xyzzyaaad173=1,dimensionality
xyzzyaaaf173=xyzzyaaaf173+1
xyzzyaaan173=const_int(xyzzyaaaf173)
xyzzyaaaw173=xyzzyaaap173(xyzzyaaad173)
select case(xyzzyaaan173)
case(0)
xyzzyaaaq173(xyzzyaaad173)=1.d0
case(1)
xyzzyaaaq173(xyzzyaaad173)=xyzzyaaaw173
case(2)
xyzzyaaaq173(xyzzyaaad173)=xyzzyaaaw173*xyzzyaaaw173
case default
xyzzyaaaq173(xyzzyaaad173)=xyzzyaaaw173**xyzzyaaan173
end select
enddo
xyzzyaaas173=xyzzyaaas173+xyzzyaaat173*product(xyzzyaaaq173(1:dimensio&
&nality))
enddo
xyzzyaaar173=xyzzyaaav173*xyzzyaaas173
f(1,xyzzyaaaa173)=f(1,xyzzyaaaa173)*xyzzyaaar173
enddo
end subroutine xyzzyaaip1
subroutine xyzzyaaiq1(const_int,const_dble,p,ne,ieskip,ion1,exvecs,fst&
&ride,f,gradf,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,ion1,fstride,const_int(*)
real(dp),intent(in) :: const_dble(*),p(4),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa174,xyzzyaaab174,xyzzyaaac174,xyzzyaaad174,xyzzyaaae1&
&74,xyzzyaaaf174,xyzzyaaag174,xyzzyaaah174,xyzzyaaai174,xyzzyaaaj174,x&
&yzzyaaak174,xyzzyaaal174,xyzzyaaam174,xyzzyaaan174,xyzzyaaao174,xyzzy&
&aaap174
real(dp) xyzzyaaaq174(3),xyzzyaaar174(3),xyzzyaaas174,xyzzyaaat174(3),&
&xyzzyaaau174,xyzzyaaav174(3),xyzzyaaaw174(3),xyzzyaaax174,xyzzyaaay17&
&4,xyzzyaaaz174(3),xyzzyaaba174(3),xyzzyaabb174(3),xyzzyaabc174(3),xyz&
&zyaabd174,xyzzyaabe174,xyzzyaabf174,xyzzyaabg174,xyzzyaabh174,xyzzyaa&
&bi174(fstride,ne),xyzzyaabj174(3)
call xyzzyaaih1(const_int(1),p(1),ne,ieskip,exvecs(1,1),fstride,f(1,1)&
&,xyzzyaabi174(1,1),nzcut(1))
xyzzyaaaj174=dimensionality*dimensionality
xyzzyaaak174=const_int(2)
xyzzyaaab174=const_int(3)
xyzzyaaan174=const_int(4)
xyzzyaaam174=4+xyzzyaaab174*dimensionality
xyzzyaabf174=dble(xyzzyaaak174)
do xyzzyaaaa174=1,ne
if(xyzzyaaaa174==ieskip)cycle
if(.not.nzcut(xyzzyaaaa174))cycle
if(xyzzyaaaa174>ieskip)then
xyzzyaaaq174(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa174)
else
xyzzyaaaq174(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa174)
endif
xyzzyaabe174=exvecs(4,xyzzyaaaa174)
call xyzzyaaiz1(xyzzyaabe174,exvecs(1,xyzzyaaaa174),xyzzyaabj174)
if(xyzzyaabe174==0.d0)then
xyzzyaaay174=0.d0
xyzzyaaax174=0.d0
xyzzyaaaz174=0.d0
else
xyzzyaaay174=xyzzyaabe174**(-xyzzyaaak174-2)
xyzzyaaax174=xyzzyaaay174*xyzzyaabe174*xyzzyaabe174
xyzzyaaaz174(1:dimensionality)=-xyzzyaabf174*xyzzyaaay174*xyzzyaaaq174&
&(1:dimensionality)
endif
if(xyzzyaaan174/=0)then
if(xyzzyaaan174<0.or.ion1==0)then
xyzzyaaal174=1
xyzzyaaao174=xyzzyaaab174
else
xyzzyaaal174=const_int(xyzzyaaam174+ion1+xyzzyaaaa174-1)
xyzzyaaao174=xyzzyaaab174+(xyzzyaaal174-1)*(xyzzyaaaj174+xyzzyaaaj174)
endif
xyzzyaaaf174=xyzzyaaao174+1
do xyzzyaaad174=1,dimensionality
xyzzyaaar174(xyzzyaaad174)=ddot(dimensionality,xyzzyaaaq174(1),1,const&
&_dble(xyzzyaaaf174),1)
xyzzyaaaf174=xyzzyaaaf174+dimensionality
enddo
else
xyzzyaaar174(1:dimensionality)=xyzzyaaaq174(1:dimensionality)
endif
xyzzyaaau174=0.d0
xyzzyaaaw174=0.d0
xyzzyaaah174=4
xyzzyaaai174=0
do xyzzyaaac174=1,xyzzyaaab174
xyzzyaaai174=xyzzyaaai174+1
xyzzyaabd174=const_dble(xyzzyaaai174)
do xyzzyaaad174=1,dimensionality
xyzzyaaah174=xyzzyaaah174+1
xyzzyaaap174=const_int(xyzzyaaah174)
xyzzyaabg174=xyzzyaaar174(xyzzyaaad174)
select case(xyzzyaaap174)
case(0)
xyzzyaabc174(xyzzyaaad174)=0.d0
xyzzyaabb174(xyzzyaaad174)=1.d0
case(1)
xyzzyaabc174(xyzzyaaad174)=1.d0
xyzzyaabb174(xyzzyaaad174)=xyzzyaabg174
case(2)
xyzzyaabc174(xyzzyaaad174)=2.d0*xyzzyaabg174
xyzzyaabb174(xyzzyaaad174)=xyzzyaabg174*xyzzyaabg174
case default
xyzzyaabh174=xyzzyaabg174**(xyzzyaaap174-1)
xyzzyaabc174(xyzzyaaad174)=dble(xyzzyaaap174)*xyzzyaabh174
xyzzyaabb174(xyzzyaaad174)=xyzzyaabh174*xyzzyaabg174
end select
enddo
xyzzyaaau174=xyzzyaaau174+xyzzyaabd174*product(xyzzyaabb174(1:dimensio&
&nality))
do xyzzyaaad174=1,dimensionality
xyzzyaaba174(xyzzyaaad174)=xyzzyaabc174(xyzzyaaad174)
do xyzzyaaae174=1,dimensionality
if(xyzzyaaae174/=xyzzyaaad174)xyzzyaaba174(xyzzyaaad174)=xyzzyaaba174(&
&xyzzyaaad174)*xyzzyaabb174(xyzzyaaae174)
enddo
enddo
xyzzyaaaw174(1:dimensionality)=xyzzyaaaw174(1:dimensionality)+xyzzyaab&
&d174*xyzzyaaba174(1:dimensionality)
enddo
if(xyzzyaaan174/=0)then
xyzzyaaag174=xyzzyaaao174
do xyzzyaaad174=1,dimensionality
xyzzyaaag174=xyzzyaaag174+1
xyzzyaaav174(xyzzyaaad174)=ddot(dimensionality,xyzzyaaaw174(1),1,const&
&_dble(xyzzyaaag174),dimensionality)
enddo
else
xyzzyaaav174(1:dimensionality)=xyzzyaaaw174(1:dimensionality)
endif
xyzzyaaat174(1:dimensionality)=xyzzyaaaz174(1:dimensionality)*xyzzyaaa&
&u174+xyzzyaaax174*xyzzyaaav174(1:dimensionality)
xyzzyaaas174=xyzzyaaax174*xyzzyaaau174
gradf(1:dimensionality,1,xyzzyaaaa174)=xyzzyaabj174(1:dimensionality)*&
&xyzzyaabi174(1,xyzzyaaaa174)*xyzzyaaas174+f(1,xyzzyaaaa174)*xyzzyaaat&
&174(1:dimensionality)
f(1,xyzzyaaaa174)=f(1,xyzzyaaaa174)*xyzzyaaas174
enddo
end subroutine xyzzyaaiq1
subroutine xyzzyaair1(const_int,const_dble,p,ne,ieskip,ion1,exvecs,fst&
&ride,f,gradf,lapf,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,ion1,fstride,const_int(*)
real(dp),intent(in) :: const_dble(*),p(4),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),gradf(3,fstride,*),lapf(fstride&
&,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa175,xyzzyaaab175,xyzzyaaac175,xyzzyaaad175,xyzzyaaae1&
&75,xyzzyaaaf175,xyzzyaaag175,xyzzyaaah175,xyzzyaaai175,xyzzyaaaj175,x&
&yzzyaaak175,xyzzyaaal175,xyzzyaaam175,xyzzyaaan175,xyzzyaaao175,xyzzy&
&aaap175,xyzzyaaaq175,xyzzyaaar175,xyzzyaaas175
real(dp) xyzzyaaat175(3),xyzzyaaau175(3),xyzzyaaav175,xyzzyaaaw175(3),&
&xyzzyaaax175,xyzzyaaay175,xyzzyaaaz175(3),xyzzyaaba175(3),xyzzyaabb17&
&5,xyzzyaabc175(3,3),xyzzyaabd175,xyzzyaabe175(3),xyzzyaabf175,xyzzyaa&
&bg175(3),xyzzyaabh175(3,3),xyzzyaabi175(3),xyzzyaabj175(3),xyzzyaabk1&
&75(3),xyzzyaabl175,xyzzyaabm175,xyzzyaabn175,xyzzyaabo175,xyzzyaabp17&
&5,xyzzyaabq175,xyzzyaabr175,xyzzyaabs175,xyzzyaabt175(fstride,ne),xyz&
&zyaabu175(fstride,ne),xyzzyaabv175(3),xyzzyaabw175
call xyzzyaaii1(const_int(1),p(1),ne,ieskip,exvecs(1,1),fstride,f(1,1)&
&,xyzzyaabt175(1,1),xyzzyaabu175(1,1),nzcut(1))
xyzzyaaal175=dimensionality*dimensionality
xyzzyaaam175=const_int(2)
xyzzyaaab175=const_int(3)
xyzzyaaap175=const_int(4)
xyzzyaaao175=4+xyzzyaaab175*dimensionality
xyzzyaabo175=dble(xyzzyaaam175)
xyzzyaabp175=dble(xyzzyaaam175*(xyzzyaaam175+2-dimensionality))
do xyzzyaaaa175=1,ne
if(xyzzyaaaa175==ieskip)cycle
if(.not.nzcut(xyzzyaaaa175))cycle
if(xyzzyaaaa175>ieskip)then
xyzzyaaat175(1:dimensionality)=exvecs(1:dimensionality,xyzzyaaaa175)
else
xyzzyaaat175(1:dimensionality)=-exvecs(1:dimensionality,xyzzyaaaa175)
endif
xyzzyaabm175=exvecs(4,xyzzyaaaa175)
call xyzzyaaiy1(xyzzyaabm175,exvecs(1,xyzzyaaaa175),xyzzyaabv175,xyzzy&
&aabw175)
if(xyzzyaabm175==0.d0)then
xyzzyaabd175=0.d0
xyzzyaabe175=0.d0
xyzzyaabf175=0.d0
else
xyzzyaabn175=xyzzyaabm175**(-xyzzyaaam175-2)
xyzzyaabd175=xyzzyaabn175*xyzzyaabm175*xyzzyaabm175
xyzzyaabe175(1:dimensionality)=-xyzzyaabo175*xyzzyaabn175*xyzzyaaat175&
&(1:dimensionality)
xyzzyaabf175=xyzzyaabp175*xyzzyaabn175
endif
if(xyzzyaaap175/=0)then
if(xyzzyaaap175<0.or.ion1==0)then
xyzzyaaan175=1
xyzzyaaaq175=xyzzyaaab175
else
xyzzyaaan175=const_int(xyzzyaaao175+ion1+xyzzyaaaa175-1)
xyzzyaaaq175=xyzzyaaab175+(xyzzyaaan175-1)*(xyzzyaaal175+xyzzyaaal175)
endif
xyzzyaaar175=xyzzyaaaq175+xyzzyaaal175
xyzzyaaag175=xyzzyaaaq175+1
do xyzzyaaad175=1,dimensionality
xyzzyaaau175(xyzzyaaad175)=ddot(dimensionality,xyzzyaaat175(1),1,const&
&_dble(xyzzyaaag175),1)
xyzzyaaag175=xyzzyaaag175+dimensionality
enddo
else
xyzzyaaau175(1:dimensionality)=xyzzyaaat175(1:dimensionality)
endif
xyzzyaaay175=0.d0
xyzzyaaba175=0.d0
xyzzyaabc175=0.d0
xyzzyaaaj175=4
xyzzyaaak175=0
do xyzzyaaac175=1,xyzzyaaab175
xyzzyaaak175=xyzzyaaak175+1
xyzzyaabl175=const_dble(xyzzyaaak175)
do xyzzyaaad175=1,dimensionality
xyzzyaaaj175=xyzzyaaaj175+1
xyzzyaaas175=const_int(xyzzyaaaj175)
xyzzyaabq175=xyzzyaaau175(xyzzyaaad175)
select case(xyzzyaaas175)
case(0)
xyzzyaabk175(xyzzyaaad175)=0.d0
xyzzyaabj175(xyzzyaaad175)=0.d0
xyzzyaabi175(xyzzyaaad175)=1.d0
case(1)
xyzzyaabk175(xyzzyaaad175)=0.d0
xyzzyaabj175(xyzzyaaad175)=1.d0
xyzzyaabi175(xyzzyaaad175)=xyzzyaabq175
case(2)
xyzzyaabk175(xyzzyaaad175)=2.d0
xyzzyaabj175(xyzzyaaad175)=2.d0*xyzzyaabq175
xyzzyaabi175(xyzzyaaad175)=xyzzyaabq175*xyzzyaabq175
case default
xyzzyaabr175=xyzzyaabq175**(xyzzyaaas175-2)
xyzzyaabs175=xyzzyaabq175*xyzzyaabr175
xyzzyaabk175(xyzzyaaad175)=dble(xyzzyaaas175*(xyzzyaaas175-1))*xyzzyaa&
&br175
xyzzyaabj175(xyzzyaaad175)=dble(xyzzyaaas175)*xyzzyaabs175
xyzzyaabi175(xyzzyaaad175)=xyzzyaabq175*xyzzyaabs175
end select
enddo
xyzzyaaay175=xyzzyaaay175+xyzzyaabl175*product(xyzzyaabi175(1:dimensio&
&nality))
do xyzzyaaad175=1,dimensionality
xyzzyaabg175(xyzzyaaad175)=xyzzyaabj175(xyzzyaaad175)
do xyzzyaaae175=1,dimensionality
if(xyzzyaaae175/=xyzzyaaad175)xyzzyaabg175(xyzzyaaad175)=xyzzyaabg175(&
&xyzzyaaad175)*xyzzyaabi175(xyzzyaaae175)
enddo
enddo
xyzzyaaba175(1:dimensionality)=xyzzyaaba175(1:dimensionality)+xyzzyaab&
&l175*xyzzyaabg175(1:dimensionality)
if(xyzzyaaam175>1)then
do xyzzyaaad175=1,dimensionality
xyzzyaabh175(xyzzyaaad175,xyzzyaaad175)=xyzzyaabk175(xyzzyaaad175)
do xyzzyaaae175=1,xyzzyaaad175-1
xyzzyaabh175(xyzzyaaad175,xyzzyaaad175)=xyzzyaabh175(xyzzyaaad175,xyzz&
&yaaad175)*xyzzyaabi175(xyzzyaaae175)
enddo
do xyzzyaaae175=xyzzyaaad175+1,dimensionality
xyzzyaabh175(xyzzyaaad175,xyzzyaaad175)=xyzzyaabh175(xyzzyaaad175,xyzz&
&yaaad175)*xyzzyaabi175(xyzzyaaae175)
xyzzyaabh175(xyzzyaaae175,xyzzyaaad175)=xyzzyaabj175(xyzzyaaad175)*xyz&
&zyaabj175(xyzzyaaae175)
do xyzzyaaaf175=1,dimensionality
if(xyzzyaaaf175/=xyzzyaaad175.and.xyzzyaaaf175/=xyzzyaaae175)xyzzyaabh&
&175(xyzzyaaae175,xyzzyaaad175)=xyzzyaabh175(xyzzyaaae175,xyzzyaaad175&
&)*xyzzyaabi175(xyzzyaaaf175)
enddo
xyzzyaabh175(xyzzyaaad175,xyzzyaaae175)=xyzzyaabh175(xyzzyaaae175,xyzz&
&yaaad175)
enddo
enddo
xyzzyaabc175(1:dimensionality,1:dimensionality)=xyzzyaabc175(1:dimensi&
&onality,1:dimensionality)+xyzzyaabl175*xyzzyaabh175(1:dimensionality,&
&1:dimensionality)
endif
enddo
if(xyzzyaaap175/=0)then
xyzzyaaah175=xyzzyaaaq175
do xyzzyaaad175=1,dimensionality
xyzzyaaah175=xyzzyaaah175+1
xyzzyaaaz175(xyzzyaaad175)=ddot(dimensionality,xyzzyaaba175(1),1,const&
&_dble(xyzzyaaah175),dimensionality)
enddo
xyzzyaabb175=0.d0
if(xyzzyaaam175>1)then
xyzzyaaai175=xyzzyaaar175
do xyzzyaaad175=1,dimensionality
do xyzzyaaae175=1,dimensionality
xyzzyaaai175=xyzzyaaai175+1
xyzzyaabb175=xyzzyaabb175+xyzzyaabc175(xyzzyaaae175,xyzzyaaad175)*cons&
&t_dble(xyzzyaaai175)
enddo
enddo
endif
else
xyzzyaaaz175(1:dimensionality)=xyzzyaaba175(1:dimensionality)
xyzzyaabb175=0.d0
if(xyzzyaaam175>1)then
do xyzzyaaad175=1,dimensionality
xyzzyaabb175=xyzzyaabb175+xyzzyaabc175(xyzzyaaad175,xyzzyaaad175)
enddo
endif
endif
xyzzyaaax175=xyzzyaabf175*xyzzyaaay175+2.d0*ddot(dimensionality,xyzzya&
&abe175(1),1,xyzzyaaaz175(1),1)+xyzzyaabd175*xyzzyaabb175
xyzzyaaaw175(1:dimensionality)=xyzzyaabe175(1:dimensionality)*xyzzyaaa&
&y175+xyzzyaabd175*xyzzyaaaz175(1:dimensionality)
xyzzyaaav175=xyzzyaabd175*xyzzyaaay175
lapf(1,xyzzyaaaa175)=(xyzzyaabu175(1,xyzzyaaaa175)+xyzzyaabt175(1,xyzz&
&yaaaa175)*xyzzyaabw175)*xyzzyaaav175+2.d0*xyzzyaabt175(1,xyzzyaaaa175&
&)*ddot(dimensionality,xyzzyaabv175(1),1,xyzzyaaaw175(1),1)+f(1,xyzzya&
&aaa175)*xyzzyaaax175
gradf(1:dimensionality,1,xyzzyaaaa175)=xyzzyaabv175(1:dimensionality)*&
&xyzzyaabt175(1,xyzzyaaaa175)*xyzzyaaav175+f(1,xyzzyaaaa175)*xyzzyaaaw&
&175(1:dimensionality)
f(1,xyzzyaaaa175)=f(1,xyzzyaaaa175)*xyzzyaaav175
enddo
end subroutine xyzzyaair1
subroutine xyzzyaais1(p,ne,ieskip,exvecs,fstride,f,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa176
real(dp) xyzzyaaab176,xyzzyaaac176,xyzzyaaad176,xyzzyaaae176,xyzzyaaaf&
&176,xyzzyaaag176,xyzzyaaah176,xyzzyaaai176,xyzzyaaaj176,xyzzyaaak176,&
&xyzzyaaal176,xyzzyaaam176,xyzzyaaan176
xyzzyaaal176=p(1)
xyzzyaaam176=p(2)
xyzzyaaan176=p(3)
do xyzzyaaaa176=1,ne
if(xyzzyaaaa176==ieskip)cycle
xyzzyaaab176=exvecs(4,xyzzyaaaa176)
xyzzyaaac176=0.d0
if(dimensionality==1)xyzzyaaac176=exvecs(2,xyzzyaaaa176)
xyzzyaaad176=exvecs(3,xyzzyaaaa176)
if((xyzzyaaac176==0.d0.and.xyzzyaaad176==0.d0).or.xyzzyaaab176>=xyzzya&
&aal176)then
f(1,xyzzyaaaa176)=0.d0
nzcut(xyzzyaaaa176)=.false.
else
xyzzyaaah176=xyzzyaaab176*xyzzyaaam176
xyzzyaaai176=xyzzyaaah176*xyzzyaaah176
xyzzyaaaj176=xyzzyaaai176*xyzzyaaah176
xyzzyaaae176=xyzzyaaad176*xyzzyaaad176
if(dimensionality==1)xyzzyaaae176=xyzzyaaae176+xyzzyaaac176*xyzzyaaac1&
&76
xyzzyaaaf176=sqrt(xyzzyaaab176*xyzzyaaab176+xyzzyaaae176)
xyzzyaaag176=sqrt(xyzzyaaan176+xyzzyaaae176)
xyzzyaaak176=1.d0-6*xyzzyaaai176+8*xyzzyaaaj176-3*xyzzyaaai176*xyzzyaa&
&ai176
f(1,xyzzyaaaa176)=(xyzzyaaaf176-xyzzyaaag176)*xyzzyaaak176
nzcut(xyzzyaaaa176)=.true.
endif
enddo
end subroutine xyzzyaais1
subroutine xyzzyaait1(p,ne,ieskip,exvecs,fstride,f,df,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa177
real(dp) xyzzyaaab177,xyzzyaaac177,xyzzyaaad177,xyzzyaaae177,xyzzyaaaf&
&177,xyzzyaaag177,xyzzyaaah177,xyzzyaaai177,xyzzyaaaj177,xyzzyaaak177,&
&xyzzyaaal177,xyzzyaaam177,xyzzyaaan177,xyzzyaaao177,xyzzyaaap177,xyzz&
&yaaaq177,xyzzyaaar177,xyzzyaaas177
xyzzyaaaq177=p(1)
xyzzyaaar177=p(2)
xyzzyaaas177=p(3)
do xyzzyaaaa177=1,ne
if(xyzzyaaaa177==ieskip)cycle
xyzzyaaab177=exvecs(4,xyzzyaaaa177)
xyzzyaaac177=0.d0
if(dimensionality==1)xyzzyaaac177=exvecs(2,xyzzyaaaa177)
xyzzyaaad177=exvecs(3,xyzzyaaaa177)
if((xyzzyaaac177==0.d0.and.xyzzyaaad177==0.d0).or.xyzzyaaab177>=xyzzya&
&aaq177)then
f(1,xyzzyaaaa177)=0.d0
df(1,xyzzyaaaa177)=0.d0
nzcut(xyzzyaaaa177)=.false.
else
xyzzyaaaj177=xyzzyaaab177*xyzzyaaar177
xyzzyaaak177=xyzzyaaaj177*xyzzyaaaj177
xyzzyaaal177=xyzzyaaak177*xyzzyaaaj177
xyzzyaaae177=xyzzyaaad177*xyzzyaaad177
if(dimensionality==1)xyzzyaaae177=xyzzyaaae177+xyzzyaaac177*xyzzyaaac1&
&77
xyzzyaaaf177=sqrt(xyzzyaaab177*xyzzyaaab177+xyzzyaaae177)
xyzzyaaai177=sqrt(xyzzyaaas177+xyzzyaaae177)
xyzzyaaah177=1.d0/xyzzyaaaf177
xyzzyaaag177=xyzzyaaab177*xyzzyaaah177
xyzzyaaam177=1.d0-6*xyzzyaaak177+8*xyzzyaaal177-3*xyzzyaaak177*xyzzyaa&
&ak177
xyzzyaaan177=xyzzyaaar177*12*(-xyzzyaaaj177+2*xyzzyaaak177-xyzzyaaal17&
&7)
xyzzyaaao177=(xyzzyaaaf177-xyzzyaaai177)*xyzzyaaam177
xyzzyaaap177=xyzzyaaag177*xyzzyaaam177+(xyzzyaaaf177-xyzzyaaai177)*xyz&
&zyaaan177
f(1,xyzzyaaaa177)=xyzzyaaao177
df(1,xyzzyaaaa177)=xyzzyaaap177
nzcut(xyzzyaaaa177)=.true.
endif
enddo
end subroutine xyzzyaait1
subroutine xyzzyaaiu1(p,ne,ieskip,exvecs,fstride,f,df,d2f,nzcut)
implicit none
integer,intent(in) :: ne,ieskip,fstride
real(dp),intent(in) :: p(3),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa178
real(dp) xyzzyaaab178,xyzzyaaac178,xyzzyaaad178,xyzzyaaae178,xyzzyaaaf&
&178,xyzzyaaag178,xyzzyaaah178,xyzzyaaai178,xyzzyaaaj178,xyzzyaaak178,&
&xyzzyaaal178,xyzzyaaam178,xyzzyaaan178,xyzzyaaao178,xyzzyaaap178,xyzz&
&yaaaq178,xyzzyaaar178,xyzzyaaas178,xyzzyaaat178,xyzzyaaau178,xyzzyaaa&
&v178
xyzzyaaat178=p(1)
xyzzyaaau178=p(2)
xyzzyaaav178=p(3)
do xyzzyaaaa178=1,ne
if(xyzzyaaaa178==ieskip)cycle
xyzzyaaab178=exvecs(4,xyzzyaaaa178)
xyzzyaaac178=0.d0
if(dimensionality==1)xyzzyaaac178=exvecs(2,xyzzyaaaa178)
xyzzyaaad178=exvecs(3,xyzzyaaaa178)
if((xyzzyaaac178==0.d0.and.xyzzyaaad178==0.d0).or.xyzzyaaab178>=xyzzya&
&aat178)then
f(1,xyzzyaaaa178)=0.d0
df(1,xyzzyaaaa178)=0.d0
d2f(1,xyzzyaaaa178)=0.d0
nzcut(xyzzyaaaa178)=.false.
else
xyzzyaaak178=xyzzyaaab178*xyzzyaaau178
xyzzyaaal178=xyzzyaaak178*xyzzyaaak178
xyzzyaaam178=xyzzyaaal178*xyzzyaaak178
xyzzyaaae178=xyzzyaaad178*xyzzyaaad178
if(dimensionality==1)xyzzyaaae178=xyzzyaaae178+xyzzyaaac178*xyzzyaaac1&
&78
xyzzyaaaf178=sqrt(xyzzyaaab178*xyzzyaaab178+xyzzyaaae178)
xyzzyaaaj178=sqrt(xyzzyaaav178+xyzzyaaae178)
xyzzyaaai178=1.d0/xyzzyaaaf178
xyzzyaaag178=xyzzyaaab178*xyzzyaaai178
xyzzyaaah178=(1.d0-xyzzyaaab178*xyzzyaaag178*xyzzyaaai178)*xyzzyaaai17&
&8
xyzzyaaan178=1.d0-6*xyzzyaaal178+8*xyzzyaaam178-3*xyzzyaaal178*xyzzyaa&
&al178
xyzzyaaao178=xyzzyaaau178*12*(-xyzzyaaak178+2*xyzzyaaal178-xyzzyaaam17&
&8)
xyzzyaaap178=xyzzyaaau178*xyzzyaaau178*12*(-1.d0+4*xyzzyaaak178-3*xyzz&
&yaaal178)
xyzzyaaaq178=(xyzzyaaaf178-xyzzyaaaj178)*xyzzyaaan178
xyzzyaaar178=xyzzyaaag178*xyzzyaaan178+(xyzzyaaaf178-xyzzyaaaj178)*xyz&
&zyaaao178
xyzzyaaas178=xyzzyaaah178*xyzzyaaan178+2*xyzzyaaag178*xyzzyaaao178+(xy&
&zzyaaaf178-xyzzyaaaj178)*xyzzyaaap178
f(1,xyzzyaaaa178)=xyzzyaaaq178
df(1,xyzzyaaaa178)=xyzzyaaar178
d2f(1,xyzzyaaaa178)=xyzzyaaas178
nzcut(xyzzyaaaa178)=.true.
endif
enddo
end subroutine xyzzyaaiu1
subroutine xyzzyaaiv1(ctrunc,p,ne,ieskip,exvecs,fstride,f,nzcut)
implicit none
integer,intent(in) :: ctrunc,ne,ieskip,fstride
real(dp),intent(in) :: p(4),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa179
real(dp) xyzzyaaab179,xyzzyaaac179,xyzzyaaad179,xyzzyaaae179,xyzzyaaaf&
&179,xyzzyaaag179,xyzzyaaah179
xyzzyaaae179=p(1)
xyzzyaaaf179=p(2)
xyzzyaaag179=p(3)
do xyzzyaaaa179=1,ne
if(xyzzyaaaa179==ieskip)cycle
xyzzyaaab179=exvecs(4,xyzzyaaaa179)
if(xyzzyaaab179>=xyzzyaaae179)then
f(1,xyzzyaaaa179)=0.d0
nzcut(xyzzyaaaa179)=.false.
elseif(xyzzyaaab179<xyzzyaaaf179)then
f(1,xyzzyaaaa179)=1.d0
nzcut(xyzzyaaaa179)=.true.
else
xyzzyaaac179=(xyzzyaaab179-xyzzyaaaf179)*xyzzyaaag179
select case(ctrunc)
case(0)
xyzzyaaah179=1.d0-xyzzyaaac179
case(1)
xyzzyaaah179=1.d0-xyzzyaaac179*xyzzyaaac179*(3.d0-2.d0*xyzzyaaac179)
case(2)
xyzzyaaah179=1.d0+xyzzyaaac179*xyzzyaaac179*xyzzyaaac179*(-10.d0+xyzzy&
&aaac179*(15.d0-6.d0*xyzzyaaac179))
case(3)
xyzzyaaad179=xyzzyaaac179*xyzzyaaac179
xyzzyaaah179=1.d0+xyzzyaaad179*xyzzyaaad179*(-35.d0+xyzzyaaac179*(84.d&
&0+xyzzyaaac179*(-70.d0+xyzzyaaac179*20.d0)))
end select
f(1,xyzzyaaaa179)=xyzzyaaah179
nzcut(xyzzyaaaa179)=.true.
endif
enddo
end subroutine xyzzyaaiv1
subroutine xyzzyaaiw1(ctrunc,p,ne,ieskip,exvecs,fstride,f,df,nzcut)
implicit none
integer,intent(in) :: ctrunc,ne,ieskip,fstride
real(dp),intent(in) :: p(4),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa180
real(dp) xyzzyaaab180,xyzzyaaac180,xyzzyaaad180,xyzzyaaae180,xyzzyaaaf&
&180,xyzzyaaag180,xyzzyaaah180,xyzzyaaai180
xyzzyaaae180=p(1)
xyzzyaaaf180=p(2)
xyzzyaaag180=p(3)
do xyzzyaaaa180=1,ne
if(xyzzyaaaa180==ieskip)cycle
xyzzyaaab180=exvecs(4,xyzzyaaaa180)
if(xyzzyaaab180>=xyzzyaaae180)then
f(1,xyzzyaaaa180)=0.d0
df(1,xyzzyaaaa180)=0.d0
nzcut(xyzzyaaaa180)=.false.
elseif(xyzzyaaab180<xyzzyaaaf180)then
f(1,xyzzyaaaa180)=1.d0
df(1,xyzzyaaaa180)=0.d0
nzcut(xyzzyaaaa180)=.true.
else
xyzzyaaac180=(xyzzyaaab180-xyzzyaaaf180)*xyzzyaaag180
select case(ctrunc)
case(0)
xyzzyaaah180=1.d0-xyzzyaaac180
xyzzyaaai180=-1.d0
case(1)
xyzzyaaah180=1.d0-xyzzyaaac180*xyzzyaaac180*(3.d0-2.d0*xyzzyaaac180)
xyzzyaaai180=-6.d0*xyzzyaaac180*(1.d0-xyzzyaaac180)
case(2)
xyzzyaaad180=xyzzyaaac180*xyzzyaaac180
xyzzyaaah180=1.d0+xyzzyaaad180*xyzzyaaac180*(-10.d0+xyzzyaaac180*(15.d&
&0-6.d0*xyzzyaaac180))
xyzzyaaai180=-30.d0*xyzzyaaad180*(1.d0-xyzzyaaac180)*(1.d0-xyzzyaaac18&
&0)
case(3)
xyzzyaaad180=xyzzyaaac180*xyzzyaaac180
xyzzyaaah180=1.d0+xyzzyaaad180*xyzzyaaad180*(-35.d0+xyzzyaaac180*(84.d&
&0+xyzzyaaac180*(-70.d0+xyzzyaaac180*20.d0)))
xyzzyaaai180=-140.d0*xyzzyaaad180*xyzzyaaac180*(1.d0-xyzzyaaac180)*(1.&
&d0-xyzzyaaac180)*(1.d0-xyzzyaaac180)
end select
f(1,xyzzyaaaa180)=xyzzyaaah180
df(1,xyzzyaaaa180)=xyzzyaaai180*xyzzyaaag180
nzcut(xyzzyaaaa180)=.true.
endif
enddo
end subroutine xyzzyaaiw1
subroutine xyzzyaaix1(ctrunc,p,ne,ieskip,exvecs,fstride,f,df,d2f,nzcut&
&)
implicit none
integer,intent(in) :: ctrunc,ne,ieskip,fstride
real(dp),intent(in) :: p(4),exvecs(4,ne)
real(dp),intent(inout) :: f(fstride,*),df(fstride,*),d2f(fstride,*)
logical,intent(inout) :: nzcut(ne)
integer xyzzyaaaa181
real(dp) xyzzyaaab181,xyzzyaaac181,xyzzyaaad181,xyzzyaaae181,xyzzyaaaf&
&181,xyzzyaaag181,xyzzyaaah181,xyzzyaaai181,xyzzyaaaj181,xyzzyaaak181
xyzzyaaae181=p(1)
xyzzyaaaf181=p(2)
xyzzyaaag181=p(3)
xyzzyaaah181=p(4)
do xyzzyaaaa181=1,ne
if(xyzzyaaaa181==ieskip)cycle
xyzzyaaab181=exvecs(4,xyzzyaaaa181)
if(xyzzyaaab181>=xyzzyaaae181)then
f(1,xyzzyaaaa181)=0.d0
df(1,xyzzyaaaa181)=0.d0
d2f(1,xyzzyaaaa181)=0.d0
nzcut(xyzzyaaaa181)=.false.
elseif(xyzzyaaab181<xyzzyaaaf181)then
f(1,xyzzyaaaa181)=1.d0
df(1,xyzzyaaaa181)=0.d0
d2f(1,xyzzyaaaa181)=0.d0
nzcut(xyzzyaaaa181)=.true.
else
xyzzyaaac181=(xyzzyaaab181-xyzzyaaaf181)*xyzzyaaag181
select case(ctrunc)
case(0)
xyzzyaaai181=1.d0-xyzzyaaac181
xyzzyaaaj181=-1.d0
xyzzyaaak181=0.d0
case(1)
xyzzyaaai181=1.d0-xyzzyaaac181*xyzzyaaac181*(3.d0-2.d0*xyzzyaaac181)
xyzzyaaaj181=-6.d0*xyzzyaaac181*(1.d0-xyzzyaaac181)
xyzzyaaak181=-6.d0*(1.d0-xyzzyaaac181-xyzzyaaac181)
case(2)
xyzzyaaad181=xyzzyaaac181*xyzzyaaac181
xyzzyaaai181=1.d0+xyzzyaaad181*xyzzyaaac181*(-10.d0+xyzzyaaac181*(15.d&
&0-6.d0*xyzzyaaac181))
xyzzyaaaj181=-30.d0*xyzzyaaad181*(1.d0-xyzzyaaac181)*(1.d0-xyzzyaaac18&
&1)
xyzzyaaak181=-60.d0*xyzzyaaac181*(1.d0-xyzzyaaac181)*(1.d0-xyzzyaaac18&
&1-xyzzyaaac181)
case(3)
xyzzyaaad181=xyzzyaaac181*xyzzyaaac181
xyzzyaaai181=1.d0+xyzzyaaad181*xyzzyaaad181*(-35.d0+xyzzyaaac181*(84.d&
&0+xyzzyaaac181*(-70.d0+xyzzyaaac181*20.d0)))
xyzzyaaaj181=-140.d0*xyzzyaaad181*xyzzyaaac181*(1.d0-xyzzyaaac181)*(1.&
&d0-xyzzyaaac181)*(1.d0-xyzzyaaac181)
xyzzyaaak181=420.d0*xyzzyaaad181*(-1.d0+xyzzyaaac181*(4.d0+xyzzyaaac18&
&1*(-5.d0+xyzzyaaac181*2.d0)))
end select
f(1,xyzzyaaaa181)=xyzzyaaai181
df(1,xyzzyaaaa181)=xyzzyaaaj181*xyzzyaaag181
d2f(1,xyzzyaaaa181)=xyzzyaaak181*xyzzyaaah181
nzcut(xyzzyaaaa181)=.true.
endif
enddo
end subroutine xyzzyaaix1
subroutine xyzzyaaiy1(r,vecr,grad,lap)
implicit none
real(dp),intent(in) :: r,vecr(3)
real(dp),intent(out) :: grad(3),lap
real(dp) xyzzyaaaa182
if(r==0.d0)then
grad(1:3)=0.d0
lap=0.d0
return
endif
select case(dimensionality)
case(3)
xyzzyaaaa182=1.d0/r
lap=2*xyzzyaaaa182
grad(1:3)=vecr(1:3)*xyzzyaaaa182
case(2)
xyzzyaaaa182=1.d0/r
lap=xyzzyaaaa182
grad(1:2)=vecr(1:2)*xyzzyaaaa182
grad(3)=0.d0
case(1)
grad(1)=sign(1.d0,vecr(1))
grad(2:3)=0.d0
lap=0.d0
end select
end subroutine xyzzyaaiy1
subroutine xyzzyaaiz1(r,vecr,grad)
implicit none
real(dp),intent(in) :: r,vecr(3)
real(dp),intent(out) :: grad(3)
real(dp) xyzzyaaaa183
if(r==0.d0)then
grad(1:3)=0.d0
return
endif
select case(dimensionality)
case(3)
xyzzyaaaa183=1.d0/r
grad(1:3)=vecr(1:3)*xyzzyaaaa183
case(2)
xyzzyaaaa183=1.d0/r
grad(1:2)=vecr(1:2)*xyzzyaaaa183
grad(3)=0.d0
case(1)
grad(1)=sign(1.d0,vecr(1))
grad(2:3)=0.d0
end select
end subroutine xyzzyaaiz1
subroutine xyzzyaaja1(no_stars,max_gmod,no_g_in_star,g_int,g_name,excl&
&ude_zero)
implicit none
integer,intent(inout) :: no_stars
real(dp),intent(inout) :: max_gmod
logical,intent(in),optional :: exclude_zero
integer,pointer :: no_g_in_star(:)
integer,pointer,optional :: g_int(:)
character(casl_fullkeysize),pointer,optional :: g_name(:)
integer xyzzyaaaa184,xyzzyaaab184,xyzzyaaac184,xyzzyaaad184,xyzzyaaae1&
&84,xyzzyaaaf184,xyzzyaaag184,xyzzyaaah184(3),xyzzyaaai184,xyzzyaaaj18&
&4,xyzzyaaak184,xyzzyaaal184,xyzzyaaam184,xyzzyaaan184,xyzzyaaao184,xy&
&zzyaaap184,xyzzyaaaq184
integer,allocatable :: xyzzyaaar184(:,:),xyzzyaaas184(:)
real(dp) xyzzyaaat184(3),xyzzyaaau184,xyzzyaaav184,xyzzyaaaw184(3,3),x&
&yzzyaaax184
real(dp),parameter :: xyzzyaaay184=1.d-8
real(dp),allocatable :: xyzzyaaaz184(:)
logical xyzzyaaba184
character(80) char_80
xyzzyaaba184=.false.
if(present(exclude_zero))xyzzyaaba184=exclude_zero
if(no_stars<0.eqv.max_gmod<0.d0)call errstop_master('MAKE_G_STARS','Ne&
&ither/both a maximum number of stars or/and maximum |G| have been sup&
&plied. Need only one of them. Bug in calling routine.')
if(.not.no_stars<0)then
if(xyzzyaael1)call wout('DEBUG: Generating fixed number of stars '//tr&
&im(i2s(no_stars))//'.')
xyzzyaaah184=0
xyzzyaaah184(1:periodicity)=no_stars-1
else
if(xyzzyaael1)call wout('DEBUG: Generating G-stars up to cut-off ',max&
&_gmod,'.')
xyzzyaaax184=max_gmod*max_gmod
xyzzyaaah184=0
xyzzyaaaw184=0.d0
select case(periodicity)
case(3)
do xyzzyaaap184=1,3
do xyzzyaaaq184=1,3
xyzzyaaaw184(xyzzyaaaq184,xyzzyaaap184)=bmat(mod(1+xyzzyaaaq184-1,3)+1&
&,mod(1+xyzzyaaap184-1,3)+1)*bmat(mod(2+xyzzyaaaq184-1,3)+1,mod(2+xyzz&
&yaaap184-1,3)+1)-bmat(mod(2+xyzzyaaaq184-1,3)+1,mod(1+xyzzyaaap184-1,&
&3)+1)*bmat(mod(1+xyzzyaaaq184-1,3)+1,mod(2+xyzzyaaap184-1,3)+1)
enddo
enddo
case(2)
do xyzzyaaap184=1,2
xyzzyaaaw184(1:2,xyzzyaaap184)=(/-bmat(2,3-xyzzyaaap184),bmat(1,3-xyzz&
&yaaap184)/)
enddo
case(1)
xyzzyaaaw184(1,1)=1.d0
end select
do xyzzyaaap184=1,periodicity
xyzzyaaaw184(1:periodicity,xyzzyaaap184)=xyzzyaaaw184(1:periodicity,xy&
&zzyaaap184)/sqrt(ddot(periodicity,xyzzyaaaw184(1,xyzzyaaap184),1,xyzz&
&yaaaw184(1,xyzzyaaap184),1))
if(xyzzyaael1)call wout('DEBUG: '//trim(i2s(periodicity))//'D uvec_'//&
&trim(i2s(xyzzyaaap184))//' = ',xyzzyaaaw184(1:periodicity,xyzzyaaap18&
&4))
enddo
do xyzzyaaap184=1,periodicity
xyzzyaaah184(xyzzyaaap184)=ceiling(max_gmod/abs(ddot(periodicity,xyzzy&
&aaaw184(1,xyzzyaaap184),1,bmat(1,xyzzyaaap184),1)))
if(xyzzyaael1)call wout('DEBUG: no_shells_'//trim(i2s(xyzzyaaap184))//&
&' = '//trim(i2s(xyzzyaaah184(xyzzyaaap184))))
enddo
endif
xyzzyaaai184=((2*xyzzyaaah184(1)+1)*(2*xyzzyaaah184(2)+1)*(2*xyzzyaaah&
&184(3)+1))/2
if(.not.xyzzyaaba184)xyzzyaaai184=xyzzyaaai184+1
if(xyzzyaael1)call wout('DEBUG: no_G_tot = '//trim(i2s(xyzzyaaai184)))
allocate(xyzzyaaar184(3,xyzzyaaai184),xyzzyaaaz184(xyzzyaaai184),xyzzy&
&aaas184(xyzzyaaai184),stat=xyzzyaaaa184)
call check_alloc(xyzzyaaaa184,'MAKE_G_STARS','work')
xyzzyaaab184=0
do xyzzyaaad184=0,xyzzyaaah184(1)
xyzzyaaaj184=-xyzzyaaah184(2)
if(xyzzyaaad184==0)xyzzyaaaj184=0
do xyzzyaaae184=xyzzyaaaj184,xyzzyaaah184(2)
xyzzyaaak184=-xyzzyaaah184(3)
if(xyzzyaaad184==0.and.xyzzyaaae184==0)xyzzyaaak184=0
do xyzzyaaaf184=xyzzyaaak184,xyzzyaaah184(3)
if(xyzzyaaba184.and.xyzzyaaad184==0.and.xyzzyaaae184==0.and.xyzzyaaaf1&
&84==0)cycle
xyzzyaaab184=xyzzyaaab184+1
xyzzyaaat184(1:3)=xyzzyaaad184*bmat(1:3,1)+xyzzyaaae184*bmat(1:3,2)+xy&
&zzyaaaf184*bmat(1:3,3)
xyzzyaaaz184(xyzzyaaab184)=ddot(3,xyzzyaaat184,1,xyzzyaaat184,1)
xyzzyaaar184(1:3,xyzzyaaab184)=(/xyzzyaaad184,xyzzyaaae184,xyzzyaaaf18&
&4/)
if(xyzzyaael1)then
char_80=trim(i2s(xyzzyaaar184(1,xyzzyaaab184)))
do xyzzyaaan184=2,periodicity
char_80=trim(char_80)//', '//trim(i2s(xyzzyaaar184(xyzzyaaan184,xyzzya&
&aab184)))
enddo
call wout('DEBUG: unsorted integer vector #'//trim(i2s(xyzzyaaab184))/&
&/' = [ '//trim(char_80)//' ] ; sq = ',xyzzyaaaz184(xyzzyaaab184))
endif
enddo
enddo
enddo
call quicksort(xyzzyaaai184,xyzzyaaaz184,xyzzyaaas184)
if(xyzzyaael1)then
do xyzzyaaab184=1,xyzzyaaai184
char_80=trim(i2s(xyzzyaaar184(1,xyzzyaaas184(xyzzyaaab184))))
do xyzzyaaan184=2,periodicity
char_80=trim(char_80)//', '//trim(i2s(xyzzyaaar184(xyzzyaaan184,xyzzya&
&aas184(xyzzyaaab184))))
enddo
call wout('DEBUG: sorted integer vector #'//trim(i2s(xyzzyaaab184))//'&
& = [ '//trim(char_80)//' ] ; sq = ',xyzzyaaaz184(xyzzyaaas184(xyzzyaa&
&ab184)))
enddo
endif
if(periodicity==3)then
xyzzyaaav184=xyzzyaaay184*max(ddot(3,bmat(1,1),1,bmat(1,1),1),ddot(3,b&
&mat(1,2),1,bmat(1,2),1),ddot(3,bmat(1,3),1,bmat(1,3),1))
elseif(periodicity==2)then
xyzzyaaav184=xyzzyaaay184*max(ddot(2,bmat(1,1),1,bmat(1,1),1),ddot(2,b&
&mat(1,2),1,bmat(1,2),1))
else
xyzzyaaav184=xyzzyaaay184*abs(bmat(1,1))
endif
if(xyzzyaael1)call wout('DEBUG: determined G^2 comparison tolerance = &
&',xyzzyaaav184)
allocate(no_g_in_star(xyzzyaaai184),stat=xyzzyaaaa184)
call check_alloc(xyzzyaaaa184,'MAKE_G_STARS','no_G_in_star')
xyzzyaaag184=0
xyzzyaaau184=-10.d0*xyzzyaaav184
no_g_in_star=0
do xyzzyaaab184=1,xyzzyaaai184
xyzzyaaac184=xyzzyaaas184(xyzzyaaab184)
if(xyzzyaaaz184(xyzzyaaac184)-xyzzyaaau184>xyzzyaaav184)then
if(.not.no_stars<0)then
if(xyzzyaaag184==no_stars)then
max_gmod=sqrt(xyzzyaaau184)
if(xyzzyaael1)call wout('DEBUG: determined max |G| = ',max_gmod)
exit
endif
else
if(xyzzyaaaz184(xyzzyaaac184)>=xyzzyaaax184)then
no_stars=xyzzyaaag184
if(xyzzyaael1)call wout('DEBUG: determined no_stars = '//trim(i2s(no_s&
&tars)))
exit
endif
endif
xyzzyaaag184=xyzzyaaag184+1
xyzzyaaau184=xyzzyaaaz184(xyzzyaaac184)
endif
no_g_in_star(xyzzyaaag184)=no_g_in_star(xyzzyaaag184)+1
enddo
call resize_pointer((/no_stars/),no_g_in_star)
xyzzyaaal184=sum(no_g_in_star)
if(xyzzyaael1)call wout('DEBUG: total number of G vectors = '//trim(i2&
&s(xyzzyaaal184)))
if(present(g_int))then
allocate(g_int(xyzzyaaal184*periodicity),stat=xyzzyaaaa184)
call check_alloc(xyzzyaaaa184,'MAKE_G_STARS','G_int')
g_int=0
endif
if(present(g_name))then
allocate(g_name(xyzzyaaal184*periodicity),stat=xyzzyaaaa184)
call check_alloc(xyzzyaaaa184,'MAKE_G_STARS','G_name')
endif
xyzzyaaam184=0
xyzzyaaab184=0
do xyzzyaaag184=1,no_stars
do xyzzyaaao184=1,no_g_in_star(xyzzyaaag184)
xyzzyaaab184=xyzzyaaab184+1
xyzzyaaac184=xyzzyaaas184(xyzzyaaab184)
do xyzzyaaan184=1,periodicity
xyzzyaaam184=xyzzyaaam184+1
if(present(g_name))g_name(xyzzyaaam184)='Star '//trim(i2s(xyzzyaaag184&
&))//':G_'//trim(i2s(xyzzyaaao184))//':%u'//trim(i2s(xyzzyaaan184))
if(present(g_int))g_int(xyzzyaaam184)=xyzzyaaar184(xyzzyaaan184,xyzzya&
&aac184)
enddo
if(xyzzyaael1)then
char_80=trim(i2s(xyzzyaaar184(1,xyzzyaaas184(xyzzyaaab184))))
do xyzzyaaan184=2,periodicity
char_80=trim(char_80)//', '//trim(i2s(xyzzyaaar184(xyzzyaaan184,xyzzya&
&aas184(xyzzyaaab184))))
enddo
call wout('DEBUG: star #'//trim(i2s(xyzzyaaag184))//', vector #'//trim&
&(i2s(xyzzyaaao184))//' = '//trim(char_80))
endif
enddo
enddo
deallocate(xyzzyaaar184,xyzzyaaaz184,xyzzyaaas184)
end subroutine xyzzyaaja1
subroutine xyzzyaajb1(no_stars,no_g_in_star,ngvec,g,which_unity)
implicit none
integer,intent(in) :: no_stars,no_g_in_star(no_stars),ngvec
integer,intent(inout) :: which_unity
real(dp),intent(in) :: g(periodicity+1,ngvec)
integer xyzzyaaaa185,xyzzyaaab185,xyzzyaaac185,xyzzyaaad185,xyzzyaaae1&
&85
real(dp),parameter :: xyzzyaaaf185=1.d-8
logical xyzzyaaag185
which_unity=0
xyzzyaaab185=0
do xyzzyaaaa185=1,no_stars
xyzzyaaag185=.true.
do xyzzyaaad185=1,no_g_in_star(xyzzyaaaa185)
xyzzyaaab185=xyzzyaaab185+1
xyzzyaaac185=xyzzyaaab185
do xyzzyaaae185=xyzzyaaad185+1,no_g_in_star(xyzzyaaaa185)
xyzzyaaac185=xyzzyaaac185+1
if(abs(g(periodicity+1,xyzzyaaab185)-g(periodicity+1,xyzzyaaac185))>xy&
&zzyaaaf185)then
call errwarn('CHECK_G_STARS','In star #'//trim(i2s(xyzzyaaaa185))//': &
&G vectors '//trim(i2s(xyzzyaaad185))//' and '//trim(i2s(xyzzyaaae185)&
&)//' have different sizes.')
elseif(sqrt(sum((g(1:periodicity,xyzzyaaab185)+g(1:periodicity,xyzzyaa&
&ac185))**2))<xyzzyaaaf185)then
call errwarn('CHECK_G_STARS','In star #'//trim(i2s(xyzzyaaaa185))//': &
&G vectors '//trim(i2s(xyzzyaaad185))//' and '//trim(i2s(xyzzyaaae185)&
&)//' are negatives of each other.')
elseif(abs(g(periodicity+1,xyzzyaaab185))>xyzzyaaaf185)then
xyzzyaaag185=.false.
endif
enddo
enddo
if(xyzzyaaag185.and.which_unity<1)which_unity=xyzzyaaaa185
enddo
end subroutine xyzzyaajb1
subroutine xyzzyaajc1(nvector,an,bn)
implicit none
integer, intent(inout) :: nvector
double precision, pointer :: an(:,:),bn(:,:)
real(dp) xyzzyaaaa186,xyzzyaaab186,xyzzyaaac186,xyzzyaaad186(3),xyzzya&
&aae186,xyzzyaaaf186,xyzzyaaag186,xyzzyaaah186,xyzzyaaai186,xyzzyaaaj1&
&86,xyzzyaaak186,xyzzyaaal186,xyzzyaaam186,xyzzyaaan186,xyzzyaaao186,x&
&yzzyaaap186,xyzzyaaaq186
integer xyzzyaaar186
real(dp),parameter :: xyzzyaaas186=1.d-13,xyzzyaaat186=1.d-6
xyzzyaaaf186=ddot(3,amat(1,1),3,amat(1,1),3)
xyzzyaaag186=ddot(3,amat(1,1),3,amat(2,1),3)
xyzzyaaah186=ddot(3,amat(1,1),3,amat(3,1),3)
xyzzyaaai186=ddot(3,amat(2,1),3,amat(2,1),3)
xyzzyaaaj186=ddot(3,amat(2,1),3,amat(3,1),3)
xyzzyaaak186=ddot(3,amat(3,1),3,amat(3,1),3)
if(dimensionality==3.and.abs(xyzzyaaaf186-xyzzyaaai186)<=xyzzyaaat186.&
&and.abs(xyzzyaaaf186-xyzzyaaak186)<=xyzzyaaat186.and.abs(abs(xyzzyaaa&
&g186)-0.5d0*xyzzyaaaf186)<=xyzzyaaat186.and.abs(abs(xyzzyaaah186)-0.5&
&d0*xyzzyaaaf186)<=xyzzyaaat186.and.abs(abs(xyzzyaaaj186)-0.5d0*xyzzya&
&aaf186)<=xyzzyaaat186)then
nvector=4
allocate(an(3,nvector),bn(3,nvector),stat=xyzzyaaar186)
call check_alloc(xyzzyaaar186,'MAKE_AB_VECTORS','an,bn')
bn(:,1)=bmat(1,:)
bn(:,2)=bmat(2,:)
bn(:,3)=bmat(3,:)
bn(:,4)=bmat(1,:)+bmat(2,:)+bmat(3,:)
xyzzyaaae186=0.25d0/xyzzyaaje1(bn(:,1),bn(:,2),bn(:,3))
an(:,1)=(xyzzyaajd1(bn(:,1),bn(:,3))-xyzzyaajd1(bn(:,1),bn(:,2))+3.d0*&
&xyzzyaajd1(bn(:,2),bn(:,3)))*xyzzyaaae186
an(:,2)=(xyzzyaajd1(bn(:,2),bn(:,1))-xyzzyaajd1(bn(:,2),bn(:,3))+3.d0*&
&xyzzyaajd1(bn(:,3),bn(:,1)))*xyzzyaaae186
an(:,3)=(xyzzyaajd1(bn(:,3),bn(:,2))-xyzzyaajd1(bn(:,3),bn(:,1))+3.d0*&
&xyzzyaajd1(bn(:,1),bn(:,2)))*xyzzyaaae186
an(:,4)=(xyzzyaajd1(bn(:,1),bn(:,2))-xyzzyaajd1(bn(:,1),bn(:,3))+xyzzy&
&aajd1(bn(:,2),bn(:,3)))*xyzzyaaae186
elseif(dimensionality==3.and.abs(xyzzyaaaf186-xyzzyaaai186)<=xyzzyaaat&
&186.and.abs(xyzzyaaaf186-xyzzyaaak186)<=xyzzyaaat186.and.abs(abs(xyzz&
&yaaag186)-xyzzyaaaf186/3.d0)<=xyzzyaaat186.and.abs(abs(xyzzyaaah186)-&
&xyzzyaaaf186/3.d0)<=xyzzyaaat186.and.abs(abs(xyzzyaaaj186)-xyzzyaaaf1&
&86/3.d0)<=xyzzyaaat186)then
nvector=6
allocate(an(3,nvector),bn(3,nvector),stat=xyzzyaaar186)
call check_alloc(xyzzyaaar186,'MAKE_AB_VECTORS','an,bn')
xyzzyaaaa186=sign(1.d0,sum(bmat(1,:)))
xyzzyaaab186=sign(1.d0,sum(bmat(2,:)))
xyzzyaaac186=sign(1.d0,sum(bmat(3,:)))
bn(:,1)=xyzzyaaaa186*bmat(1,:)
bn(:,2)=xyzzyaaab186*bmat(2,:)
bn(:,3)=xyzzyaaac186*bmat(3,:)
bn(:,4)=xyzzyaaaa186*bmat(1,:)-xyzzyaaab186*bmat(2,:)
bn(:,5)=xyzzyaaaa186*bmat(1,:)-xyzzyaaac186*bmat(3,:)
bn(:,6)=xyzzyaaab186*bmat(2,:)-xyzzyaaac186*bmat(3,:)
xyzzyaaae186=0.250/xyzzyaaje1(bn(:,1),bn(:,2),bn(:,3))
an(:,1)=(2.d0*xyzzyaajd1(bn(:,2),bn(:,3))-xyzzyaajd1(bn(:,1),bn(:,3))+&
&xyzzyaajd1(bn(:,1),bn(:,2)))*xyzzyaaae186
an(:,2)=(2.d0*xyzzyaajd1(bn(:,3),bn(:,1))-xyzzyaajd1(bn(:,2),bn(:,1))+&
&xyzzyaajd1(bn(:,2),bn(:,3)))*xyzzyaaae186
an(:,3)=(2.d0*xyzzyaajd1(bn(:,1),bn(:,2))-xyzzyaajd1(bn(:,3),bn(:,2))+&
&xyzzyaajd1(bn(:,3),bn(:,1)))*xyzzyaaae186
an(:,4)=(xyzzyaajd1(bn(:,1),bn(:,3))+xyzzyaajd1(bn(:,2),bn(:,3)))*xyzz&
&yaaae186
an(:,5)=(xyzzyaajd1(bn(:,2),bn(:,3))+xyzzyaajd1(bn(:,2),bn(:,1)))*xyzz&
&yaaae186
an(:,6)=(xyzzyaajd1(bn(:,2),bn(:,1))+xyzzyaajd1(bn(:,3),bn(:,1)))*xyzz&
&yaaae186
elseif(dimensionality==2.and.abs(abs(xyzzyaaag186)-0.5d0*xyzzyaaaf186)&
&<xyzzyaaat186.and.abs(abs(xyzzyaaag186)-0.5d0*xyzzyaaai186)<xyzzyaaat&
&186)then
nvector=3
allocate(an(3,nvector),bn(3,nvector),stat=xyzzyaaar186)
call check_alloc(xyzzyaaar186,'MAKE_AB_VECTORS','an,bn')
bn(:,1)=bmat(1,:)
bn(:,2)=bmat(2,:)
bn(:,3)=bmat(1,:)+bmat(2,:)
xyzzyaaad186=(/0.d0,0.d0,1.d0/)
xyzzyaaae186=1.d0/(3.d0*xyzzyaaje1(xyzzyaaad186,bn(:,1),bn(:,2)))
an(:,1)=xyzzyaajd1(bn(:,1)+2.d0*bn(:,2),xyzzyaaad186)*xyzzyaaae186
an(:,2)=xyzzyaajd1(xyzzyaaad186,bn(:,2)+2.d0*bn(:,1))*xyzzyaaae186
an(:,3)=xyzzyaajd1(bn(:,2)-bn(:,1),xyzzyaaad186)*xyzzyaaae186
elseif(dimensionality==3.and.(abs(xyzzyaaaf186-xyzzyaaai186)<xyzzyaaat&
&186.or.abs(xyzzyaaaf186-xyzzyaaak186)<xyzzyaaat186.or.abs(xyzzyaaai18&
&6-xyzzyaaak186)<xyzzyaaat186).and.(abs(abs(xyzzyaaag186)-0.5d0*xyzzya&
&aaf186)<xyzzyaaat186.or.abs(abs(xyzzyaaah186)-0.5d0*xyzzyaaaf186)<xyz&
&zyaaat186.or.abs(abs(xyzzyaaaj186)-0.5d0*xyzzyaaai186)<xyzzyaaat186))&
&then
nvector=4
allocate(an(3,nvector),bn(3,nvector),stat=xyzzyaaar186)
call check_alloc(xyzzyaaar186,'MAKE_AB_VECTORS','an,bn')
xyzzyaaal186=ddot(3,bmat(1,1),3,bmat(1,1),3)
xyzzyaaam186=ddot(3,bmat(1,1),3,bmat(2,1),3)
xyzzyaaan186=ddot(3,bmat(1,1),3,bmat(3,1),3)
xyzzyaaao186=ddot(3,bmat(2,1),3,bmat(2,1),3)
xyzzyaaap186=ddot(3,bmat(2,1),3,bmat(3,1),3)
xyzzyaaaq186=ddot(3,bmat(3,1),3,bmat(3,1),3)
if(abs(xyzzyaaan186)<xyzzyaaat186.and.abs(xyzzyaaal186-xyzzyaaao186)<x&
&yzzyaaat186)then
bn(:,1)=bmat(1,:)
bn(:,2)=bmat(2,:)
bn(:,3)=bmat(3,:)
bn(:,4)=bn(:,1)+bn(:,2)
elseif(abs(xyzzyaaam186)<xyzzyaaat186.and.abs(xyzzyaaal186-xyzzyaaaq18&
&6)<xyzzyaaat186)then
bn(:,1)=bmat(1,:)
bn(:,2)=bmat(3,:)
bn(:,3)=bmat(2,:)
bn(:,4)=bn(:,1)+bn(:,2)
elseif(abs(xyzzyaaam186)<xyzzyaaat186.and.abs(xyzzyaaao186-xyzzyaaaq18&
&6)<xyzzyaaat186)then
bn(:,1)=bmat(2,:)
bn(:,2)=bmat(3,:)
bn(:,3)=bmat(1,:)
bn(:,4)=bn(:,1)+bn(:,2)
else
call errstop_master('MAKE_AB_VECTORS','Could not identify precise 3D h&
&exagonal geometry.')
endif
xyzzyaaae186=1.d0/(3.d0*xyzzyaaje1(bn(:,1),bn(:,2),bn(:,3)))
an(:,1)=(xyzzyaajd1(bn(:,1),bn(:,3))+2.d0*xyzzyaajd1(bn(:,2),bn(:,3)))&
&*xyzzyaaae186
an(:,2)=(xyzzyaajd1(bn(:,3),bn(:,2))+2.d0*xyzzyaajd1(bn(:,3),bn(:,1)))&
&*xyzzyaaae186
an(:,3)=3.d0*xyzzyaajd1(bn(:,1),bn(:,2))*xyzzyaaae186
an(:,4)=(xyzzyaajd1(bn(:,2)-bn(:,1),bn(:,3)))*xyzzyaaae186
else
nvector=dimensionality
allocate(an(3,nvector),bn(3,nvector),stat=xyzzyaaar186)
call check_alloc(xyzzyaaar186,'MAKE_AB_VECTORS','an,bn')
bn(:,1)=bmat(1,:)
if(dimensionality>1)bn(:,2)=bmat(2,:)
if(dimensionality>2)bn(:,3)=bmat(3,:)
an(:,1)=amat(1,:)/twopi
if(dimensionality>1)an(:,2)=amat(2,:)/twopi
if(dimensionality>2)an(:,3)=amat(3,:)/twopi
endif
end subroutine xyzzyaajc1
function xyzzyaajd1(a,b) result(c)
implicit none
real(dp),intent(in) :: a(3),b(3)
real(dp) c(3)
c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)
end function xyzzyaajd1
real(dp) function xyzzyaaje1(a,b,c)
implicit none
real(dp),intent(in) :: a(3),b(3),c(3)
real(dp) xyzzyaaaa188(3)
xyzzyaaaa188=xyzzyaajd1(b,c)
xyzzyaaje1=ddot(3,a(1),1,xyzzyaaaa188(1),1)
end function xyzzyaaje1
subroutine xyzzyaajf1(phi,f,df,d2f)
implicit none
real(dp),intent(in) :: phi
real(dp),intent(inout) :: f
real(dp),intent(inout),optional :: df,d2f
integer,parameter :: xyzzyaaaa189=3
real(dp),parameter :: xyzzyaaab189=1.d0/dble(xyzzyaaaa189+1)
real(dp) xyzzyaaac189,xyzzyaaad189,xyzzyaaae189,xyzzyaaaf189,xyzzyaaag&
&189
xyzzyaaac189=modulo(phi+pi,twopi)-pi
xyzzyaaad189=abs(xyzzyaaac189)
xyzzyaaae189=xyzzyaaad189*one_over_pi
if(.not.present(df))then
f=xyzzyaaad189*(1.d0-xyzzyaaab189*xyzzyaaae189**xyzzyaaaa189)
elseif(.not.present(d2f))then
xyzzyaaaf189=xyzzyaaae189**xyzzyaaaa189
f=xyzzyaaad189*(1.d0-xyzzyaaab189*xyzzyaaaf189)
df=sign(1.d0,xyzzyaaac189)*(1.d0-xyzzyaaaf189)
else
xyzzyaaag189=xyzzyaaae189**(xyzzyaaaa189-1)
xyzzyaaaf189=xyzzyaaag189*xyzzyaaae189
f=xyzzyaaad189*(1.d0-xyzzyaaab189*xyzzyaaaf189)
df=sign(1.d0,xyzzyaaac189)*(1.d0-xyzzyaaaf189)
d2f=-dble(xyzzyaaaa189)*one_over_pi*xyzzyaaag189
endif
end subroutine xyzzyaajf1
subroutine xyzzyaajg1(phi,g,dg,d2g)
implicit none
real(dp),intent(in) :: phi
real(dp),intent(inout) :: g
real(dp),intent(inout),optional :: dg,d2g
integer,parameter :: xyzzyaaaa190=1,xyzzyaaab190=2
real(dp),parameter :: xyzzyaaac190=dble(xyzzyaaab190*(xyzzyaaab190+1))&
&/dble(xyzzyaaaa190*(xyzzyaaaa190+1)-xyzzyaaab190*(xyzzyaaab190+1)),xy&
&zzyaaad190=dble(xyzzyaaaa190*(xyzzyaaaa190+1))/dble(xyzzyaaab190*(xyz&
&zyaaab190+1)-xyzzyaaaa190*(xyzzyaaaa190+1))
real(dp) xyzzyaaae190,xyzzyaaaf190,xyzzyaaag190,xyzzyaaah190,xyzzyaaai&
&190,xyzzyaaaj190
xyzzyaaae190=modulo(phi+pi,twopi)-pi
xyzzyaaaf190=abs(xyzzyaaae190)*one_over_pi
if(.not.present(dg))then
g=xyzzyaaae190*(1.d0+xyzzyaaac190*xyzzyaaaf190**xyzzyaaaa190+xyzzyaaad&
&190*xyzzyaaaf190**xyzzyaaab190)
elseif(.not.present(d2g))then
xyzzyaaag190=xyzzyaaac190*xyzzyaaaf190**xyzzyaaaa190
xyzzyaaai190=xyzzyaaad190*xyzzyaaaf190**xyzzyaaab190
g=xyzzyaaae190*(1.d0+xyzzyaaag190+xyzzyaaai190)
dg=1.d0+dble(xyzzyaaaa190+1)*xyzzyaaag190+dble(xyzzyaaab190+1)*xyzzyaa&
&ai190
else
xyzzyaaah190=xyzzyaaac190*xyzzyaaaf190**(xyzzyaaaa190-1)
xyzzyaaag190=xyzzyaaah190*xyzzyaaaf190
xyzzyaaaj190=xyzzyaaad190*xyzzyaaaf190**(xyzzyaaab190-1)
xyzzyaaai190=xyzzyaaaj190*xyzzyaaaf190
g=xyzzyaaae190*(1.d0+xyzzyaaag190+xyzzyaaai190)
dg=1.d0+dble(xyzzyaaaa190+1)*xyzzyaaag190+dble(xyzzyaaab190+1)*xyzzyaa&
&ai190
d2g=sign(one_over_pi,xyzzyaaae190)*(dble(xyzzyaaaa190*(xyzzyaaaa190+1)&
&)*xyzzyaaah190+dble(xyzzyaaab190*(xyzzyaaab190+1))*xyzzyaaaj190)
endif
end subroutine xyzzyaajg1
end module slaarnabc

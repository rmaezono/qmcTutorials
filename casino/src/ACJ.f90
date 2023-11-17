module slaarnacj
use slaarnaad
use dsp
use slaarnach
use store
use slaarnaag,   only : c_one,czero
use format_utils,only : i2s
use slaarnabg,    only : dimensionality
use slaarnabp,        only : detcoef,setup_mdet_params,finish_mdet_par&
&ams,get_mdet_params,put_mdet_params,detcoef_label
use slaarnabt,   only : dcopy,ddot,daxpy,dscal,ddot_s,dcopy3,dcopy6,lu&
&_decom,lu_decom_cmplx_dz,lu_solve_n,lu_solve_once_cmplx_dz,lu_logdet,&
&lu_logdet_cmplx_dz,zcopy,dswap
use run_control, only : timer,check_alloc,errstop_master
use slaarnacq, only : wfdet,setup_wfdet_params,finish_wfdet_params,get&
&_wfdet_params,put_wfdet_params,get_wfdet_rmax,enumerate_plot_wfdet,qu&
&ery_plot_wfdet,get_plot_wfdet,finish_plot_wfdet,copy_orb_to_det,copy_&
&rmap_orb_to_det,wfdet_orbmap,wfdet_orbmask,wfdet_norb,wfdet_mdet_to_c&
&mdet,wfdet_detcoef_affect_orbs
implicit none
private
public enumerate_plot_slater,query_plot_slater,get_plot_slater,finish_&
&plot_slater
public query_slater_levels,query_slater_level_details,setup_slater,fin&
&ish_slater,wfn_ratio_slater,accept_move_slater,reset_config_slater,wf&
&n_logval_slater,wfn_loggrad_slater,wfn_loglap_slater,prefetch_wfn_sla&
&ter,clear_scratch_slater,add_config_slater_items,setup_slater_params,&
&finish_slater_params,get_slater_params,put_slater_params,clone_scratc&
&h_slater,invalidate_params_slater,invalidate_param1_slater,setup_stor&
&age_slater,finish_storage_slater,load_from_storage_slater,save_to_sto&
&rage_slater
public gen_config_slater,delete_config_slater,copy_config_slater,confi&
&g_to_pt_slater,pt_to_config_slater,redist_allocations_slater,redist_l&
&oad_slater,redist_send_slater,redist_recv_slater,redist_save_slater,r&
&edist_deallocations_slater,load_from_pt_slater,save_to_pt_slater
public get_slater_rmax
public config_wfn_slater
public dbarrc,small_transfer,small_buffers
integer dbarrc,xyzzyaaaa1
logical small_transfer,small_buffers
real(dp),parameter :: xyzzyaaab1=-690.d0
complex(dp),parameter :: xyzzyaaac1=(-1.d6,0.d0),xyzzyaaad1=(-2.d6,0.d&
&0)
integer,parameter :: xyzzyaaae1=0,xyzzyaaaf1=1,xyzzyaaag1=2
integer,allocatable,target :: xyzzyaaah1(:,:,:),xyzzyaaai1(:,:)
complex(dp),allocatable,target :: xyzzyaaaj1(:,:,:),xyzzyaaak1(:,:)
complex(dp),allocatable :: xyzzyaaal1(:)
logical,allocatable :: xyzzyaaam1(:)
integer,allocatable,target :: xyzzyaaan1(:)
complex(dp),allocatable :: xyzzyaaao1(:,:),xyzzyaaap1(:)
logical,allocatable :: xyzzyaaaq1(:)
complex(dp),allocatable,target :: xyzzyaaar1(:,:)
logical,allocatable :: xyzzyaaas1(:)
integer,allocatable :: xyzzyaaat1(:,:,:,:)
real(dp),allocatable :: xyzzyaaau1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaaav1(:)
integer,allocatable :: xyzzyaaaw1(:)
real(dp),allocatable,target :: xyzzyaaax1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaaay1(:)
real(dp),allocatable,target :: xyzzyaaaz1(:,:,:,:,:,:)
integer,allocatable :: dbar_age(:)
logical,allocatable :: xyzzyaaba1(:)
integer,allocatable :: xyzzyaabb1(:)
real(dp),allocatable,target :: xyzzyaabc1(:,:,:,:,:,:,:)
logical,allocatable :: xyzzyaabd1(:)
real(dp),allocatable,target :: xyzzyaabe1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaabf1(:)
real(dp),allocatable,target :: xyzzyaabg1(:,:,:,:,:,:,:)
logical,allocatable :: xyzzyaabh1(:)
real(dp),allocatable :: xyzzyaabi1(:,:,:,:)
logical,allocatable :: xyzzyaabj1(:,:)
real(dp),allocatable :: xyzzyaabk1(:,:,:)
logical,allocatable :: xyzzyaabl1(:,:)
real(dp),allocatable :: xyzzyaabm1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaabn1(:,:)
real(dp),allocatable :: xyzzyaabo1(:,:,:,:,:)
logical,allocatable :: xyzzyaabp1(:,:)
integer,allocatable :: xyzzyaabq1(:,:),xyzzyaabr1(:,:,:)
real(dp),allocatable :: xyzzyaabs1(:,:,:,:)
logical,allocatable :: xyzzyaabt1(:)
integer xyzzyaabu1,xyzzyaabv1
real(dp),allocatable :: xyzzyaabw1(:,:,:)
integer xyzzyaabx1,xyzzyaaby1
real(dp),allocatable :: xyzzyaabz1(:,:,:,:)
real(dp),allocatable,target :: xyzzyaaca1(:,:,:,:),xyzzyaacb1(:,:,:)
real(dp),allocatable :: xyzzyaacc1(:,:,:,:),xyzzyaacd1(:,:,:,:,:)
logical,allocatable :: xyzzyaace1(:)
real(dp),allocatable :: xyzzyaacf1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaacg1(:)
real(dp),allocatable,target :: xyzzyaach1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaaci1(:)
real(dp),allocatable,target :: xyzzyaacj1(:,:,:,:,:,:,:)
logical,allocatable :: xyzzyaack1(:)
integer,parameter :: xyzzyaacl1=3
integer xyzzyaacm1,xyzzyaacn1(xyzzyaacl1)
integer,allocatable :: xyzzyaaco1(:)
type config_wfn_slater
private
complex(dp),pointer :: pt_logdet(:,:)=>null()
real(dp),pointer :: pt_dbar(:,:,:,:,:)=>null(),pt_dsmat(:,:,:,:,:,:)=>&
&null(),pt_lapsmat(:,:,:,:,:)=>null(),pt_d2smat(:,:,:,:,:,:)=>null()
integer dbar_age
logical orb_dvalid,orb_d2valid
end type config_wfn_slater
integer,allocatable :: xyzzyaacp1(:)
real(dp),allocatable :: xyzzyaacq1(:,:,:,:,:,:)
complex(dp),allocatable :: xyzzyaacr1(:,:,:)
integer xyzzyaacs1
integer,allocatable :: xyzzyaact1(:),xyzzyaacu1(:)
real(dp),allocatable :: xyzzyaacv1(:,:),xyzzyaacw1(:,:,:)
real(dp),allocatable,target :: xyzzyaacx1(:,:,:),xyzzyaacy1(:,:,:,:),x&
&yzzyaacz1(:,:,:),xyzzyaada1(:,:,:,:),xyzzyaadb1(:,:,:,:)
real(dp),allocatable,target :: xyzzyaadc1(:,:),xyzzyaadd1(:,:,:),xyzzy&
&aade1(:,:),xyzzyaadf1(:,:,:)
complex(dp),allocatable,target :: xyzzyaadg1(:)
complex(dp),allocatable :: xyzzyaadh1(:,:),xyzzyaadi1(:,:,:)
real(dp),allocatable :: xyzzyaadj1(:,:,:,:),xyzzyaadk1(:,:,:),xyzzyaad&
&l1(:,:,:,:,:,:),xyzzyaadm1(:,:,:,:,:)
complex(dp),allocatable :: xyzzyaadn1(:,:,:)
logical,allocatable :: xyzzyaado1(:),xyzzyaadp1(:),xyzzyaadq1(:),xyzzy&
&aadr1(:),xyzzyaads1(:)
integer xyzzyaadt1(2),xyzzyaadu1(2),xyzzyaadv1(2),xyzzyaadw1(2),xyzzya&
&adx1(2),xyzzyaady1(2),xyzzyaadz1(2),xyzzyaaea1(2),xyzzyaaeb1(2),xyzzy&
&aaec1(2),xyzzyaaed1(2),xyzzyaaee1(2),xyzzyaaef1(2),xyzzyaaeg1(2),xyzz&
&yaaeh1(2),xyzzyaaei1(2),xyzzyaaej1(2)
integer,parameter :: xyzzyaaek1=3,xyzzyaael1=1,xyzzyaaem1=2,xyzzyaaen1&
&=3
integer :: xyzzyaaeo1(xyzzyaaek1)=0,xyzzyaaep1(0:2)=0
integer,allocatable :: xyzzyaaeq1(:)
character(64),parameter :: xyzzyaaer1(xyzzyaaek1)=(/'orb     ','orb_gr&
&ad','orb_lap '/)
character(64),parameter :: xyzzyaaes1(xyzzyaaek1)=(/'orbital values   &
&      ','gradient of orbitals   ',  'Laplacian of orbitals  '/)
contains
subroutine query_slater_levels(nlevels)
implicit none
integer,intent(out) :: nlevels
nlevels=1
end subroutine query_slater_levels
subroutine query_slater_level_details(level_score,level_name)
implicit none
integer,intent(out) :: level_score(*)
character(80),intent(out) :: level_name(*)
level_score(1)=200
level_name(1)='Slater determinants'
end subroutine query_slater_level_details
subroutine setup_slater
implicit none
integer xyzzyaaaa4,xyzzyaaab4
xyzzyaaaa1=dbarrc
xyzzyaadt1=0
xyzzyaadu1=0
xyzzyaaea1=0
xyzzyaaei1=0
xyzzyaaeb1=0
xyzzyaaej1=0
xyzzyaadv1=0
xyzzyaadw1=0
xyzzyaadx1=0
xyzzyaady1=0
xyzzyaadz1=0
xyzzyaaec1=0
xyzzyaaef1=0
xyzzyaaed1=0
xyzzyaaee1=0
xyzzyaaeg1=0
xyzzyaaeh1=0
call include_range((/1,nscratch/),xyzzyaadt1)
call include_range(ratiocfg_from_sz,xyzzyaadv1)
call include_range(ratiocfg_to_sz,xyzzyaadv1)
call include_range(ratio2_from_sz,xyzzyaadv1)
call include_range(ratio2_to_sz,xyzzyaadv1)
call include_range(ratio_ion_from_sz,xyzzyaadv1)
call include_range(ratio_ion_to_sz,xyzzyaadv1)
call include_range(drift_sz,xyzzyaaea1)
call include_range(drift_sz,xyzzyaadw1)
call include_range(drift_sz,xyzzyaadu1)
call include_range(kinetic_sz,xyzzyaaea1)
call include_range(kinetic_sz,xyzzyaaeb1)
call include_range(kinetic_sz,xyzzyaadu1)
call include_range(kinetic_detail_sz,xyzzyaaea1)
call include_range(kinetic_detail_sz,xyzzyaaeb1)
call include_range(kinetic_detail_sz,xyzzyaadu1)
call include_range(wfn_detail_sz,xyzzyaadv1)
if(.not.use_backflow)then
call include_range(ratio1_from_sz,xyzzyaadu1)
call include_range(ratio1_to_sz,xyzzyaady1)
call include_range(ratio1_to_sz,xyzzyaadz1)
call include_range(drift_sz,xyzzyaadx1)
call include_range(drift_sz,xyzzyaaei1)
call include_range(kinetic_sz,xyzzyaaei1)
call include_range(kinetic_sz,xyzzyaaej1)
call include_range(kinetic_detail_sz,xyzzyaaei1)
call include_range(kinetic_detail_sz,xyzzyaaej1)
else
call include_range(ratio1_from_sz,xyzzyaadv1)
if(bf_sparse)call include_range(ratio1_from_sz,xyzzyaadu1)
call include_range(ratio1_to_sz,xyzzyaadv1)
call include_range(ratio1_to_sz,xyzzyaaed1)
if(bf_sparse)call include_range(ratio1_to_sz,xyzzyaaef1)
if(bf_sparse.and.pairing_wf)call include_range(ratio1_to_sz,xyzzyaadu1&
&)
call include_range(drift_sz,xyzzyaadv1)
call include_range(drift_sz,xyzzyaaec1)
call include_range(drift_sz,xyzzyaaee1)
call include_range(drift_sz,xyzzyaaeg1)
call include_range(kinetic_sz,xyzzyaadv1)
call include_range(kinetic_sz,xyzzyaadw1)
call include_range(kinetic_sz,xyzzyaaec1)
call include_range(kinetic_sz,xyzzyaaeg1)
call include_range(kinetic_sz,xyzzyaaeh1)
call include_range(kinetic_detail_sz,xyzzyaadw1)
call include_range(kinetic_detail_sz,xyzzyaaec1)
call include_range(kinetic_detail_sz,xyzzyaaeg1)
call include_range(kinetic_detail_sz,xyzzyaaeh1)
endif
if(use_altsamp)then
call include_range(ratio1_to_sz,xyzzyaaei1)
call include_range(ratio1_to_sz,xyzzyaaea1)
endif
if(use_altsamp.and.use_backflow.and.simplepdf==1)then
call include_range(ratio1_from_sz,xyzzyaadu1)
call include_range(ratio1_to_sz,xyzzyaady1)
call include_range(ratio1_to_sz,xyzzyaadz1)
call include_range(drift_sz,xyzzyaadx1)
call include_range(drift_sz,xyzzyaaei1)
call include_range(kinetic_sz,xyzzyaaei1)
call include_range(kinetic_sz,xyzzyaaej1)
call include_range(kinetic_detail_sz,xyzzyaaei1)
call include_range(kinetic_detail_sz,xyzzyaaej1)
endif
if(xyzzyaadu1(1)/=0)then
allocate(dbar_age(xyzzyaadu1(1):xyzzyaadu1(2)),xyzzyaaaz1(nemax,nemax,&
&real1_complex2,nspin,ndet,xyzzyaadu1(1):xyzzyaadu1(2)),stat=xyzzyaaaa&
&4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','dbar')
xyzzyaaaz1=0.d0
dbar_age=0
endif
allocate(xyzzyaaba1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','dbar_valid')
if(xyzzyaadv1(1)/=0)then
allocate(xyzzyaaax1(nemax,real1_complex2,ndet,nemax,nspin,xyzzyaadv1(1&
&):xyzzyaadv1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','smat')
xyzzyaaax1=0.d0
endif
allocate(xyzzyaaay1(nscratch),xyzzyaaaw1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','smat_valid')
xyzzyaaaw1=0
if(xyzzyaadw1(1)/=0)then
allocate(xyzzyaabc1(3,nemax,real1_complex2,ndet,nemax,nspin,xyzzyaadw1&
&(1):xyzzyaadw1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','dsmat')
xyzzyaabc1=0.d0
endif
allocate(xyzzyaabd1(nscratch),xyzzyaabb1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','dsmat_valid')
xyzzyaabb1=0
if(xyzzyaadx1(1)/=0)then
allocate(xyzzyaabe1(nemax,real1_complex2,ndet,nemax,nspin,            &
&xyzzyaadx1(1):xyzzyaadx1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','lapsmat')
xyzzyaabe1=0.d0
endif
allocate(xyzzyaabf1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','lapsmat_valid')
if(xyzzyaaec1(1)/=0)then
allocate(xyzzyaabg1(6,nemax,real1_complex2,ndet,nemax,nspin,          &
& xyzzyaaec1(1):xyzzyaaec1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','d2smat')
xyzzyaabg1=0.d0
endif
allocate(xyzzyaabh1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','d2smat_valid')
if(xyzzyaadt1(1)/=0)then
allocate(xyzzyaaaj1(nspin,ndet,xyzzyaadt1(1):xyzzyaadt1(2)),xyzzyaaah1&
&(nspin,ndet,xyzzyaadt1(1):xyzzyaadt1(2)),xyzzyaaak1(ndet,xyzzyaadt1(1&
&):xyzzyaadt1(2)),xyzzyaaal1(xyzzyaadt1(1):xyzzyaadt1(2)),xyzzyaaai1(n&
&det,xyzzyaadt1(1):xyzzyaadt1(2)),xyzzyaaao1(ndet,xyzzyaadt1(1):xyzzya&
&adt1(2)),xyzzyaaap1(xyzzyaadt1(1):xyzzyaadt1(2)),xyzzyaaan1(xyzzyaadt&
&1(1):xyzzyaadt1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','logdet')
xyzzyaaaj1=czero
xyzzyaaah1=xyzzyaaae1
xyzzyaaak1=czero
xyzzyaaal1=czero
xyzzyaaai1=xyzzyaaae1
xyzzyaaao1=czero
xyzzyaaap1=czero
xyzzyaaan1=xyzzyaaae1
endif
allocate(xyzzyaaam1(nscratch),xyzzyaaaq1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','logdet_valid')
if(xyzzyaady1(1)/=0)then
allocate(xyzzyaaar1(ndet,xyzzyaady1(1):xyzzyaady1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','q')
xyzzyaaar1=czero
endif
allocate(xyzzyaaas1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','q_chvalid')
if(xyzzyaadz1(1)/=0)then
allocate(xyzzyaabs1(nemax,real1_complex2,ndet,xyzzyaadz1(1):xyzzyaadz1&
&(2)),xyzzyaabq1(ndet,xyzzyaadz1(1):xyzzyaadz1(2)),xyzzyaabr1(nemax,nd&
&et,xyzzyaadz1(1):xyzzyaadz1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','smat1')
xyzzyaabs1=0.d0
xyzzyaabq1=0
xyzzyaabr1=0
endif
allocate(xyzzyaabt1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','smat1_chvalid')
if(xyzzyaaea1(1)/=0)then
allocate(xyzzyaabi1(3,netot,real1_complex2,xyzzyaaea1(1):xyzzyaaea1(2)&
&),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','fidet')
xyzzyaabi1=0.d0
endif
allocate(xyzzyaabj1(netot,nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','fidet_valid')
if(xyzzyaaei1(1)/=0)then
allocate(xyzzyaabm1(3,ndet,nemax,real1_complex2,nspin,                &
&xyzzyaaei1(1):xyzzyaaei1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','fi_pdet')
xyzzyaabm1=0.d0
endif
allocate(xyzzyaabn1(netot,nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','fi_pdet_valid')
if(xyzzyaaeb1(1)/=0)then
allocate(xyzzyaabk1(netot,real1_complex2,xyzzyaaeb1(1):xyzzyaaeb1(2)),&
&stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','lapdet')
xyzzyaabk1=0.d0
endif
allocate(xyzzyaabl1(netot,nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','lapdet_valid')
if(xyzzyaaej1(1)/=0)then
allocate(xyzzyaabo1(ndet,nemax,real1_complex2,nspin,                xy&
&zzyaaej1(1):xyzzyaaej1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','lap_pdet')
xyzzyaabo1=0.d0
endif
allocate(xyzzyaabp1(netot,nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','lap_pdet_valid')
if(xyzzyaaef1(1)/=0)then
allocate(xyzzyaaau1(nemax,nemax,real1_complex2,nspin,ndet,           x&
&yzzyaaef1(1):xyzzyaaef1(2)),xyzzyaaat1(nemax,nspin,ndet,xyzzyaaef1(1)&
&:xyzzyaaef1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','sbar')
xyzzyaaau1=0.d0
xyzzyaaat1=0
endif
allocate(xyzzyaaav1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','sbar_chvalid')
if(xyzzyaaed1(1)/=0)then
allocate(xyzzyaach1(nemax,real1_complex2,ndet,nemax,nspin,            &
& xyzzyaaed1(1):xyzzyaaed1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','schunk')
xyzzyaach1=0.d0
endif
allocate(xyzzyaaci1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','schunk_chvalid')
if(xyzzyaaee1(1)/=0)then
allocate(xyzzyaacj1(3,nemax,real1_complex2,ndet,nemax,nspin,          &
&    xyzzyaaee1(1):xyzzyaaee1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','dschunk')
xyzzyaacj1=0.d0
endif
allocate(xyzzyaack1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','dschunk_chvalid')
if(xyzzyaaeg1(1)/=0)then
allocate(xyzzyaacc1(3,real1_complex2,netot,xyzzyaaeg1(1):xyzzyaaeg1(2)&
&),xyzzyaacd1(3,real1_complex2,ndet,netot,xyzzyaaeg1(1):xyzzyaaeg1(2))&
&,stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','Farray')
xyzzyaacc1=0.d0
xyzzyaacd1=0.d0
endif
allocate(xyzzyaace1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','Farray_valid')
if(xyzzyaaeh1(1)/=0)then
allocate(xyzzyaacf1(3,3,real1_complex2,netot,netot,xyzzyaaeh1(1):xyzzy&
&aaeh1(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','Harray')
xyzzyaacf1=0.d0
endif
allocate(xyzzyaacg1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','Harray_valid')
allocate(xyzzyaabw1(nemax,real1_complex2,ndet),xyzzyaaca1(3,nemax,real&
&1_complex2,ndet),xyzzyaacb1(nemax,real1_complex2,ndet),xyzzyaabz1(6,n&
&emax,real1_complex2,ndet),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SLATER','rpsi')
xyzzyaabw1=0.d0
xyzzyaabu1=0
xyzzyaabv1=0
xyzzyaaca1=0.d0
xyzzyaacb1=0.d0
xyzzyaabz1=0.d0
xyzzyaabx1=0
xyzzyaaby1=0
if(use_backflow)call setup_bf
do xyzzyaaab4=1,nscratch
call clear_scratch_slater(xyzzyaaab4)
enddo
call xyzzyaaet1
end subroutine setup_slater
subroutine xyzzyaaet1
implicit none
integer xyzzyaaaa5
allocate(xyzzyaacv1(nemax,real1_complex2),xyzzyaacx1(nemax,real1_compl&
&ex2,ndet),xyzzyaacy1(3,nemax,real1_complex2,ndet),xyzzyaacz1(nemax,re&
&al1_complex2,ndet),xyzzyaada1(6,nemax,real1_complex2,ndet),xyzzyaadb1&
&(nemax,nemax,real1_complex2,ndet),xyzzyaact1(nemax),xyzzyaadg1(ndet),&
&stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_SLATER_TEMPS','1a')
xyzzyaacv1=0.d0
xyzzyaacx1=0.d0
xyzzyaacy1=0.d0
xyzzyaacz1=0.d0
xyzzyaada1=0.d0
xyzzyaadb1=0.d0
xyzzyaact1=0
xyzzyaadg1=c_one
allocate(xyzzyaadc1(wfdet_norb,real1_complex2),xyzzyaadd1(3,wfdet_norb&
&,real1_complex2),xyzzyaade1(wfdet_norb,real1_complex2),xyzzyaadf1(6,w&
&fdet_norb,real1_complex2),xyzzyaacu1(wfdet_norb),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_SLATER_TEMPS','1b')
xyzzyaadc1=0.d0
xyzzyaadd1=0.d0
xyzzyaade1=0.d0
xyzzyaadf1=0.d0
xyzzyaacu1=0
if(complex_wf)then
allocate(xyzzyaadh1(nemax,ndet),xyzzyaadi1(nemax,nemax,ndet),stat=xyzz&
&yaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_SLATER_TEMPS','2')
xyzzyaadh1=czero
xyzzyaadi1=czero
endif
if(use_backflow)then
allocate(xyzzyaacw1(nemax,nemax,real1_complex2),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_SLATER_TEMPS','3')
xyzzyaacw1=0.d0
endif
end subroutine xyzzyaaet1
subroutine finish_slater
implicit none
if(allocated(xyzzyaaaz1))deallocate(dbar_age,xyzzyaaaz1)
deallocate(xyzzyaaba1)
if(allocated(xyzzyaaax1))deallocate(xyzzyaaax1)
deallocate(xyzzyaaay1,xyzzyaaaw1)
if(allocated(xyzzyaabc1))deallocate(xyzzyaabc1)
deallocate(xyzzyaabd1,xyzzyaabb1)
if(allocated(xyzzyaabe1))deallocate(xyzzyaabe1)
deallocate(xyzzyaabf1)
if(allocated(xyzzyaabg1))deallocate(xyzzyaabg1)
deallocate(xyzzyaabh1)
if(allocated(xyzzyaaaj1))deallocate(xyzzyaaaj1,xyzzyaaah1,xyzzyaaak1,x&
&yzzyaaai1,xyzzyaaal1,xyzzyaaao1,xyzzyaaap1,xyzzyaaan1)
deallocate(xyzzyaaam1,xyzzyaaaq1)
if(allocated(xyzzyaaar1))deallocate(xyzzyaaar1)
deallocate(xyzzyaaas1)
if(allocated(xyzzyaabs1))deallocate(xyzzyaabs1,xyzzyaabq1,xyzzyaabr1)
deallocate(xyzzyaabt1)
if(allocated(xyzzyaabi1))deallocate(xyzzyaabi1)
deallocate(xyzzyaabj1)
if(allocated(xyzzyaabm1))deallocate(xyzzyaabm1)
deallocate(xyzzyaabn1)
if(allocated(xyzzyaabk1))deallocate(xyzzyaabk1)
deallocate(xyzzyaabl1)
if(allocated(xyzzyaabo1))deallocate(xyzzyaabo1)
deallocate(xyzzyaabp1)
if(allocated(xyzzyaaau1))deallocate(xyzzyaaau1,xyzzyaaat1)
deallocate(xyzzyaaav1)
if(allocated(xyzzyaach1))deallocate(xyzzyaach1)
deallocate(xyzzyaaci1)
if(allocated(xyzzyaacj1))deallocate(xyzzyaacj1)
deallocate(xyzzyaack1)
if(allocated(xyzzyaacc1))deallocate(xyzzyaacc1,xyzzyaacd1)
deallocate(xyzzyaace1)
if(allocated(xyzzyaacf1))deallocate(xyzzyaacf1)
deallocate(xyzzyaacg1)
deallocate(xyzzyaabw1,xyzzyaaca1,xyzzyaacb1,xyzzyaabz1)
if(use_backflow)call finish_bf
call xyzzyaaeu1
end subroutine finish_slater
subroutine xyzzyaaeu1
implicit none
deallocate(xyzzyaacv1,xyzzyaacx1,xyzzyaacy1,xyzzyaacz1,xyzzyaada1,xyzz&
&yaadb1,xyzzyaact1,xyzzyaadg1)
deallocate(xyzzyaadc1,xyzzyaadd1,xyzzyaade1,xyzzyaadf1,xyzzyaacu1)
if(complex_wf)deallocate(xyzzyaadh1,xyzzyaadi1)
if(use_backflow)deallocate(xyzzyaacw1)
end subroutine xyzzyaaeu1
subroutine wfn_ratio_slater(is,js,ilevel,ratio,fd,sd,isnan,isinf)
implicit none
integer,intent(in) :: is,js,ilevel
complex(dp),intent(out) :: ratio
logical,intent(in) :: fd,sd
logical,intent(out) :: isnan,isinf
complex(dp) xyzzyaaaa8
call timer('WFN_RATIO_SLATER',.true.)
isnan=.false.
isinf=.false.
call xyzzyaaev1(is,fd,sd)
call xyzzyaaev1(js,fd,sd)
if(xyzzyaaan1(is)==xyzzyaaae1)then
if(xyzzyaaan1(js)==xyzzyaaae1)then
xyzzyaaaa8=(xyzzyaaap1(js)+xyzzyaaal1(js))-(xyzzyaaap1(is)+xyzzyaaal1(&
&is))
if(dble(xyzzyaaaa8)<xyzzyaaab1)then
ratio=czero
else
ratio=exp(xyzzyaaaa8)
endif
else
ratio=czero
endif
else
ratio=czero
if(xyzzyaaan1(js)==xyzzyaaae1)then
isinf=.true.
else
isnan=.true.
endif
endif
call timer('WFN_RATIO_SLATER',.false.)
end subroutine wfn_ratio_slater
subroutine accept_move_slater(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,xyzzyaaae9,xyzzyaa&
&af9,xyzzyaaag9,xyzzyaaah9,xyzzyaaai9
call timer('ACCEPT_MOVE_SLATER',.true.)
if(xyzzyaaaq1(js))then
xyzzyaaao1(detstart:detstop,is)=xyzzyaaao1(detstart:detstop,js)
xyzzyaaap1(is)=xyzzyaaap1(js)
xyzzyaaaq1(is)=.true.
else
xyzzyaaaq1(is)=.false.
endif
if(xyzzyaaam1(js))then
xyzzyaaaj1(:,detstart:detstop,is)=xyzzyaaaj1(:,detstart:detstop,js)
xyzzyaaak1(detstart:detstop,is)=xyzzyaaak1(detstart:detstop,js)
xyzzyaaal1(is)=xyzzyaaal1(js)
xyzzyaaam1(is)=.true.
else
xyzzyaaam1(is)=.false.
endif
if(xyzzyaaay1(is).and.xyzzyaaci1(js))then
do xyzzyaaai9=1,nspin
do xyzzyaaaf9=1,bf_m_chscr(xyzzyaaai9,js)
xyzzyaaag9=bf_rmap_chscr(xyzzyaaaf9,xyzzyaaai9,js)
xyzzyaaah9=which_ii(xyzzyaaag9,xyzzyaaai9)
call xyzzyaafu1(xyzzyaaah9,xyzzyaach1(1,1,1,xyzzyaaag9,xyzzyaaai9,js),&
&xyzzyaaax1(1,1,1,1,1,is))
enddo
enddo
xyzzyaaay1(is)=.true.
xyzzyaaaw1(is)=0
elseif(xyzzyaaay1(js))then
call dcopy(size_dbar,xyzzyaaax1(1,1,1,1,1,js),1,xyzzyaaax1(1,1,1,1,1,i&
&s),1)
xyzzyaaay1(is)=.true.
xyzzyaaaw1(is)=0
else
xyzzyaaay1(is)=.false.
xyzzyaaaw1(is)=0
endif
if(xyzzyaabd1(is).and.xyzzyaack1(js))then
do xyzzyaaai9=1,nspin
do xyzzyaaaf9=1,bf_m_chscr(xyzzyaaai9,js)
xyzzyaaag9=bf_rmap_chscr(xyzzyaaaf9,xyzzyaaai9,js)
xyzzyaaah9=which_ii(xyzzyaaag9,xyzzyaaai9)
call xyzzyaafv1(xyzzyaaah9,xyzzyaacj1(1,1,1,1,xyzzyaaag9,xyzzyaaai9,js&
&),xyzzyaabc1(1,1,1,1,1,1,is))
enddo
enddo
xyzzyaabd1(is)=.true.
xyzzyaabb1(is)=0
elseif(xyzzyaabd1(is).and.buffer_move1_from(js)==is.and.xyzzyaabx1==js&
&.and.xyzzyaaby1==buffer_move1_from_ii(js))then
xyzzyaaaa9=buffer_move1_from_ii(js)
call xyzzyaafv1(xyzzyaaaa9,xyzzyaaca1,xyzzyaabc1(1,1,1,1,1,1,is))
xyzzyaabd1(is)=.true.
xyzzyaabb1(is)=0
elseif(xyzzyaabd1(js))then
call dcopy(size_dsmat,xyzzyaabc1(1,1,1,1,1,1,js),1,xyzzyaabc1(1,1,1,1,&
&1,1,is),1)
xyzzyaabd1(is)=.true.
xyzzyaabb1(is)=0
else
xyzzyaabd1(is)=.false.
xyzzyaabb1(is)=0
endif
if(xyzzyaabh1(js))then
call dcopy(size_d2smat,xyzzyaabg1(1,1,1,1,1,1,js),1,xyzzyaabg1(1,1,1,1&
&,1,1,is),1)
xyzzyaabh1(is)=.true.
else
xyzzyaabh1(is)=.false.
endif
if(xyzzyaabf1(is).and.buffer_move1_from(js)==is.and.xyzzyaabx1==js.and&
&.xyzzyaaby1==buffer_move1_from_ii(js))then
xyzzyaaaa9=buffer_move1_from_ii(js)
call xyzzyaafw1(xyzzyaaaa9,xyzzyaacb1,xyzzyaabe1(1,1,1,1,1,is))
xyzzyaabf1(is)=.true.
elseif(xyzzyaabf1(js))then
call dcopy(size_dbar,xyzzyaabe1(1,1,1,1,1,js),1,xyzzyaabe1(1,1,1,1,1,i&
&s),1)
xyzzyaabf1(is)=.true.
else
xyzzyaabf1(is)=.false.
endif
if(xyzzyaace1(js))then
call dcopy(size_farray,xyzzyaacc1(1,1,1,js),1,xyzzyaacc1(1,1,1,is),1)
call dcopy(size_farray*ndet,xyzzyaacd1(1,1,1,1,js),1,xyzzyaacd1(1,1,1,&
&1,is),1)
xyzzyaace1(is)=.true.
else
xyzzyaace1(is)=.false.
endif
if(xyzzyaacg1(js))then
call dcopy(size_harray,xyzzyaacf1(1,1,1,1,1,js),1,xyzzyaacf1(1,1,1,1,1&
&,is),1)
xyzzyaacg1(is)=.true.
else
xyzzyaacg1(is)=.false.
endif
if(all(xyzzyaabj1(:,js)))then
call dcopy(size_fidet,xyzzyaabi1(1,1,1,js),1,xyzzyaabi1(1,1,1,is),1)
xyzzyaabj1(:,is)=.true.
else
xyzzyaabj1(:,is)=.false.
endif
if(all(xyzzyaabl1(:,js)))then
call dcopy(netot*real1_complex2,xyzzyaabk1(1,1,js),1,xyzzyaabk1(1,1,is&
&),1)
xyzzyaabl1(:,is)=.true.
else
xyzzyaabl1(:,is)=.false.
endif
if(all(xyzzyaabn1(:,js)))then
call dcopy(size_fi_prod_det,xyzzyaabm1(1,1,1,1,1,js),1,xyzzyaabm1(1,1,&
&1,1,1,is),1)
xyzzyaabn1(:,is)=.true.
else
xyzzyaabn1(:,is)=.false.
endif
if(all(xyzzyaabp1(:,js)))then
call dcopy(size_prod_lapdet,xyzzyaabo1(1,1,1,1,js),1,xyzzyaabo1(1,1,1,&
&1,is),1)
xyzzyaabp1(:,is)=.true.
else
xyzzyaabp1(:,is)=.false.
endif
if(xyzzyaaba1(js))then
call dcopy(size_dbar,xyzzyaaaz1(1,1,1,1,1,js),1,xyzzyaaaz1(1,1,1,1,1,i&
&s),1)
xyzzyaaba1(is)=.true.
dbar_age(is)=dbar_age(js)
elseif(xyzzyaaba1(is).and.xyzzyaaas1(js).and.dbar_age(is)<xyzzyaaaa1)t&
&hen
if(any(xyzzyaaah1(:,detstart:detstop,is)/=xyzzyaaae1.and.xyzzyaaah1(:,&
&detstart:detstop,js)==xyzzyaaae1))then
xyzzyaaba1(is)=.false.
else
xyzzyaaaa9=buffer_move1_from_ii(js)
xyzzyaaab9=which_ie(xyzzyaaaa9)
xyzzyaaac9=which_spin(xyzzyaaaa9)
do xyzzyaaae9=detstart,detstop
xyzzyaaad9=upd_spin(xyzzyaaac9,xyzzyaaae9)
if(xyzzyaaah1(xyzzyaaad9,xyzzyaaae9,is)/=xyzzyaaae1)cycle
if(xyzzyaaah1(xyzzyaaad9,xyzzyaaae9,js)/=xyzzyaaae1)cycle
call xyzzyaafq1(xyzzyaaab9,xyzzyaaac9,xyzzyaaae9,xyzzyaabs1(1,1,xyzzya&
&aae9,js),xyzzyaaaz1(1,1,1,xyzzyaaad9,xyzzyaaae9,is),xyzzyaaar1(xyzzya&
&aae9,js),xyzzyaabq1(xyzzyaaae9,js),xyzzyaabr1(1,xyzzyaaae9,js))
enddo
xyzzyaaba1(is)=.true.
dbar_age(is)=dbar_age(is)+1
endif
elseif(xyzzyaaba1(is).and.xyzzyaaav1(js).and.dbar_age(is)<xyzzyaaaa1)t&
&hen
if(any(xyzzyaaah1(:,detstart:detstop,is)/=xyzzyaaae1.and.xyzzyaaah1(:,&
&detstart:detstop,js)==xyzzyaaae1))then
xyzzyaaba1(is)=.false.
else
if(pairing_wf)then
do xyzzyaaae9=detstart,detstop
do xyzzyaaai9=1,nspin
if(update_by_column(xyzzyaaai9,xyzzyaaae9))cycle
xyzzyaaad9=upd_spin(xyzzyaaai9,xyzzyaaae9)
if(xyzzyaaah1(xyzzyaaad9,xyzzyaaae9,is)/=xyzzyaaae1)cycle
if(xyzzyaaah1(xyzzyaaad9,xyzzyaaae9,js)/=xyzzyaaae1)cycle
call dcopy(size_onedbar,xyzzyaaaz1(1,1,1,xyzzyaaad9,xyzzyaaae9,js),1,x&
&yzzyaaaz1(1,1,1,xyzzyaaad9,xyzzyaaae9,is),1)
call xyzzyaafr1(bf_m_chscr(xyzzyaaad9,js),bf_rmap_chscr(1,xyzzyaaad9,j&
&s),xyzzyaaad9,xyzzyaaae9,xyzzyaaaz1(1,1,1,xyzzyaaad9,xyzzyaaae9,is),x&
&yzzyaach1(1,1,1,1,1,js),xyzzyaaau1(1,1,1,xyzzyaaad9,xyzzyaaae9,js),xy&
&zzyaaat1(1,xyzzyaaad9,xyzzyaaae9,js))
enddo
enddo
do xyzzyaaae9=detstart,detstop
do xyzzyaaai9=1,nspin
if(count(upd_spin(:,xyzzyaaae9)==xyzzyaaai9)/=1)cycle
if(xyzzyaaah1(xyzzyaaai9,xyzzyaaae9,is)/=xyzzyaaae1)cycle
if(xyzzyaaah1(xyzzyaaai9,xyzzyaaae9,js)/=xyzzyaaae1)cycle
call xyzzyaafr1(bf_m_chscr(xyzzyaaai9,js),bf_rmap_chscr(1,xyzzyaaai9,j&
&s),xyzzyaaai9,xyzzyaaae9,xyzzyaaaz1(1,1,1,xyzzyaaai9,xyzzyaaae9,is),x&
&yzzyaach1(1,1,1,1,1,js),xyzzyaaau1(1,1,1,xyzzyaaai9,xyzzyaaae9,js),xy&
&zzyaaat1(1,xyzzyaaai9,xyzzyaaae9,js))
enddo
enddo
else
do xyzzyaaae9=detstart,detstop
do xyzzyaaai9=1,nspin
if(xyzzyaaah1(xyzzyaaai9,xyzzyaaae9,is)/=xyzzyaaae1)cycle
if(xyzzyaaah1(xyzzyaaai9,xyzzyaaae9,js)/=xyzzyaaae1)cycle
call xyzzyaafr1(bf_m_chscr(xyzzyaaai9,js),bf_rmap_chscr(1,xyzzyaaai9,j&
&s),xyzzyaaai9,xyzzyaaae9,xyzzyaaaz1(1,1,1,xyzzyaaai9,xyzzyaaae9,is),x&
&yzzyaach1(1,1,1,1,1,js),xyzzyaaau1(1,1,1,xyzzyaaai9,xyzzyaaae9,js),xy&
&zzyaaat1(1,xyzzyaaai9,xyzzyaaae9,js))
enddo
enddo
endif
xyzzyaaba1(is)=.true.
dbar_age(is)=dbar_age(is)+1
endif
else
xyzzyaaba1(is)=.false.
endif
xyzzyaaah1(:,detstart:detstop,is)=xyzzyaaah1(:,detstart:detstop,js)
xyzzyaaai1(detstart:detstop,is)=xyzzyaaai1(detstart:detstop,js)
xyzzyaaan1(is)=xyzzyaaan1(js)
if(use_backflow)call accept_move_bf(is,js)
if(xyzzyaabu1==is)then
xyzzyaabu1=0
xyzzyaabv1=0
elseif(xyzzyaabu1==js)then
xyzzyaabu1=is
endif
if(xyzzyaabx1==is)then
xyzzyaabx1=0
xyzzyaaby1=0
elseif(xyzzyaabx1==js)then
xyzzyaabx1=is
endif
xyzzyaaas1(is)=.false.
xyzzyaaav1(is)=.false.
xyzzyaabt1(is)=.false.
xyzzyaaci1(is)=.false.
xyzzyaack1(is)=.false.
call timer('ACCEPT_MOVE_SLATER',.false.)
end subroutine accept_move_slater
subroutine clear_scratch_slater(is)
implicit none
integer,intent(in) :: is
call timer('CLEAR_SCRATCH_SLATER',.true.)
xyzzyaaay1(is)=.false.
xyzzyaaaw1(is)=0
xyzzyaabd1(is)=.false.
xyzzyaabb1(is)=0
xyzzyaabf1(is)=.false.
xyzzyaabh1(is)=.false.
xyzzyaaam1(is)=.false.
xyzzyaaaq1(is)=.false.
xyzzyaabj1(:,is)=.false.
xyzzyaabl1(:,is)=.false.
xyzzyaabn1(:,is)=.false.
xyzzyaabp1(:,is)=.false.
xyzzyaaas1(is)=.false.
xyzzyaabt1(is)=.false.
xyzzyaaba1(is)=.false.
xyzzyaaci1(is)=.false.
xyzzyaack1(is)=.false.
xyzzyaaav1(is)=.false.
xyzzyaace1(is)=.false.
xyzzyaacg1(is)=.false.
if(xyzzyaabu1==is)then
xyzzyaabu1=0
xyzzyaabv1=0
endif
if(xyzzyaabx1==is)then
xyzzyaabx1=0
xyzzyaaby1=0
endif
if(use_backflow)call clear_scratch_bf(is)
call timer('CLEAR_SCRATCH_SLATER',.false.)
end subroutine clear_scratch_slater
subroutine reset_config_slater(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11,xyzzyaaae11,xy&
&zzyaaaf11
call timer('RESET_CONFIG_SLATER',.true.)
xyzzyaaah1(:,:,js)=xyzzyaaah1(:,:,is)
xyzzyaaan1(js)=xyzzyaaan1(is)
if(xyzzyaaaw1(js)==is)then
xyzzyaaay1(js)=.false.
xyzzyaaaw1(js)=is
elseif(.not.use_backflow.and.xyzzyaaay1(is).and.xyzzyaaay1(js))then
if(buffer_move1_from(js)==is)then
xyzzyaaaa11=buffer_move1_from_ii(js)
if(pairing_wf)then
call xyzzyaafu1(xyzzyaaaa11,xyzzyaacx1,xyzzyaaax1(1,1,1,1,1,is),back=.&
&true.)
call xyzzyaafu1(xyzzyaaaa11,xyzzyaacx1,xyzzyaaax1(1,1,1,1,1,js))
else
xyzzyaaab11=which_ie(xyzzyaaaa11)
xyzzyaaac11=which_spin(xyzzyaaaa11)
call xyzzyaafu1(xyzzyaaaa11,xyzzyaaax1(1,1,1,xyzzyaaab11,xyzzyaaac11,i&
&s),xyzzyaaax1(1,1,1,1,1,js))
endif
xyzzyaaay1(js)=.false.
xyzzyaaaw1(js)=is
elseif(buffer_move2_from(js)==is)then
xyzzyaaaa11=buffer_move2_from_ii(js)
xyzzyaaad11=buffer_move2_from_jj(js)
if(pairing_wf)then
call xyzzyaafu1(xyzzyaaaa11,xyzzyaacx1,xyzzyaaax1(1,1,1,1,1,is),back=.&
&true.)
call xyzzyaafu1(xyzzyaaaa11,xyzzyaacx1,xyzzyaaax1(1,1,1,1,1,js))
call xyzzyaafu1(xyzzyaaad11,xyzzyaacx1,xyzzyaaax1(1,1,1,1,1,is),back=.&
&true.)
call xyzzyaafu1(xyzzyaaad11,xyzzyaacx1,xyzzyaaax1(1,1,1,1,1,js))
else
xyzzyaaab11=which_ie(xyzzyaaaa11)
xyzzyaaac11=which_spin(xyzzyaaaa11)
call xyzzyaafu1(xyzzyaaaa11,xyzzyaaax1(1,1,1,xyzzyaaab11,xyzzyaaac11,i&
&s),xyzzyaaax1(1,1,1,1,1,js))
xyzzyaaae11=which_ie(xyzzyaaad11)
xyzzyaaaf11=which_spin(xyzzyaaad11)
call xyzzyaafu1(xyzzyaaad11,xyzzyaaax1(1,1,1,xyzzyaaae11,xyzzyaaaf11,i&
&s),xyzzyaaax1(1,1,1,1,1,js))
endif
xyzzyaaay1(js)=.false.
xyzzyaaaw1(js)=is
else
xyzzyaaay1(js)=.false.
xyzzyaaaw1(js)=0
endif
elseif(xyzzyaaay1(is).and.in_range(js,xyzzyaadv1))then
call dcopy(size_dbar,xyzzyaaax1(1,1,1,1,1,is),1,xyzzyaaax1(1,1,1,1,1,j&
&s),1)
xyzzyaaay1(js)=.false.
xyzzyaaaw1(js)=is
else
xyzzyaaay1(js)=.false.
xyzzyaaaw1(js)=0
endif
if(xyzzyaabb1(js)==is)then
xyzzyaabd1(js)=.false.
xyzzyaabb1(js)=is
elseif(.not.use_backflow.and.xyzzyaabd1(is).and.xyzzyaabd1(js))then
if(buffer_move1_from(js)==is)then
xyzzyaaaa11=buffer_move1_from_ii(js)
if(pairing_wf)then
call xyzzyaafv1(xyzzyaaaa11,xyzzyaacy1,xyzzyaabc1(1,1,1,1,1,1,is),back&
&=.true.)
call xyzzyaafv1(xyzzyaaaa11,xyzzyaacy1,xyzzyaabc1(1,1,1,1,1,1,js))
else
xyzzyaaab11=which_ie(xyzzyaaaa11)
xyzzyaaac11=which_spin(xyzzyaaaa11)
call xyzzyaafv1(xyzzyaaaa11,xyzzyaabc1(1,1,1,1,xyzzyaaab11,xyzzyaaac11&
&,is),xyzzyaabc1(1,1,1,1,1,1,js))
endif
xyzzyaabd1(js)=.false.
xyzzyaabb1(js)=is
elseif(buffer_move2_from(js)==is)then
xyzzyaaaa11=buffer_move2_from_ii(js)
xyzzyaaad11=buffer_move2_from_jj(js)
if(pairing_wf)then
call xyzzyaafv1(xyzzyaaaa11,xyzzyaacy1,xyzzyaabc1(1,1,1,1,1,1,is),back&
&=.true.)
call xyzzyaafv1(xyzzyaaaa11,xyzzyaacy1,xyzzyaabc1(1,1,1,1,1,1,js))
else
xyzzyaaab11=which_ie(xyzzyaaaa11)
xyzzyaaac11=which_spin(xyzzyaaaa11)
call xyzzyaafv1(xyzzyaaaa11,xyzzyaabc1(1,1,1,1,xyzzyaaab11,xyzzyaaac11&
&,is),xyzzyaabc1(1,1,1,1,1,1,js))
xyzzyaaae11=which_ie(xyzzyaaad11)
xyzzyaaaf11=which_spin(xyzzyaaad11)
call xyzzyaafv1(xyzzyaaad11,xyzzyaabc1(1,1,1,1,xyzzyaaae11,xyzzyaaaf11&
&,is),xyzzyaabc1(1,1,1,1,1,1,js))
endif
xyzzyaabd1(js)=.false.
xyzzyaabb1(js)=is
else
xyzzyaabd1(js)=.false.
xyzzyaabb1(js)=0
endif
elseif(xyzzyaabd1(is).and.in_range(js,xyzzyaadw1))then
call dcopy(size_dsmat,xyzzyaabc1(1,1,1,1,1,1,is),1,xyzzyaabc1(1,1,1,1,&
&1,1,js),1)
xyzzyaabd1(js)=.false.
xyzzyaabb1(js)=is
else
xyzzyaabd1(js)=.false.
xyzzyaabb1(js)=0
endif
xyzzyaabf1(js)=.false.
xyzzyaabh1(js)=.false.
xyzzyaaam1(js)=.false.
xyzzyaaaq1(js)=.false.
xyzzyaabj1(:,js)=.false.
xyzzyaabl1(:,js)=.false.
xyzzyaabn1(:,js)=.false.
xyzzyaabp1(:,js)=.false.
xyzzyaaas1(js)=.false.
xyzzyaabt1(js)=.false.
xyzzyaaba1(js)=.false.
xyzzyaaci1(js)=.false.
xyzzyaack1(js)=.false.
xyzzyaaav1(js)=.false.
xyzzyaace1(js)=.false.
xyzzyaacg1(js)=.false.
if(xyzzyaabu1==js)then
xyzzyaabu1=0
xyzzyaabv1=0
endif
if(xyzzyaabx1==js)then
xyzzyaabx1=0
xyzzyaaby1=0
endif
if(use_backflow)call reset_config_bf(is,js)
call timer('RESET_CONFIG_SLATER',.false.)
end subroutine reset_config_slater
subroutine clone_scratch_slater(is,js)
implicit none
integer,intent(in) :: is,js
call timer('CLONE_SCRATCH_SLATER',.true.)
if(xyzzyaaaq1(is).and..not.xyzzyaaaq1(js))then
xyzzyaaao1(:,js)=xyzzyaaao1(:,is)
xyzzyaaap1(js)=xyzzyaaap1(is)
xyzzyaaan1(js)=xyzzyaaan1(is)
xyzzyaaaq1(js)=.true.
else
xyzzyaaaq1(js)=.false.
endif
if(xyzzyaaam1(is).and..not.xyzzyaaam1(js))then
xyzzyaaaj1(:,:,js)=xyzzyaaaj1(:,:,is)
xyzzyaaak1(:,js)=xyzzyaaak1(:,is)
xyzzyaaal1(js)=xyzzyaaal1(is)
xyzzyaaah1(:,:,js)=xyzzyaaah1(:,:,is)
xyzzyaaam1(js)=.true.
else
xyzzyaaam1(js)=.false.
endif
if(scr_tasks(ikinetic,js))then
if(all(xyzzyaabj1(:,is)).and..not.all(xyzzyaabj1(:,js)))then
call dcopy(size_fidet,xyzzyaabi1(1,1,1,is),1,xyzzyaabi1(1,1,1,js),1)
xyzzyaabj1(:,js)=.true.
else
xyzzyaabj1(:,js)=.false.
endif
if(all(xyzzyaabn1(:,is)).and..not.all(xyzzyaabn1(:,js)))then
call dcopy(size_fi_prod_det,xyzzyaabm1(1,1,1,1,1,is),1,xyzzyaabm1(1,1,&
&1,1,1,js),1)
xyzzyaabn1(:,js)=.true.
else
xyzzyaabn1(:,js)=.false.
endif
if(all(xyzzyaabl1(:,is)).and..not.all(xyzzyaabl1(:,js)))then
call dcopy(netot*real1_complex2,xyzzyaabk1(1,1,is),1,xyzzyaabk1(1,1,js&
&),1)
xyzzyaabl1(:,js)=.true.
else
xyzzyaabl1(:,js)=.false.
endif
if(all(xyzzyaabp1(:,is)).and..not.all(xyzzyaabp1(:,js)))then
call dcopy(size_prod_lapdet,xyzzyaabo1(1,1,1,1,is),1,xyzzyaabo1(1,1,1,&
&1,js),1)
xyzzyaabp1(:,js)=.true.
else
xyzzyaabp1(:,js)=.false.
endif
if(xyzzyaaba1(is).and..not.xyzzyaaba1(js))then
call dcopy(size_dbar,xyzzyaaaz1(1,1,1,1,1,is),1,xyzzyaaaz1(1,1,1,1,1,j&
&s),1)
dbar_age(js)=dbar_age(is)
xyzzyaaba1(js)=.true.
endif
if(xyzzyaabd1(is).and..not.xyzzyaabd1(js))then
call dcopy(size_dsmat,xyzzyaabc1(1,1,1,1,1,1,is),1,xyzzyaabc1(1,1,1,1,&
&1,1,js),1)
xyzzyaabd1(js)=.true.
xyzzyaabb1(js)=0
else
xyzzyaabd1(js)=.false.
xyzzyaabb1(js)=0
endif
if(xyzzyaabh1(is).and..not.xyzzyaabh1(js))then
call dcopy(size_d2smat,xyzzyaabg1(1,1,1,1,1,1,is),1,xyzzyaabg1(1,1,1,1&
&,1,1,js),1)
xyzzyaabh1(js)=.true.
else
xyzzyaabh1(js)=.false.
endif
if(xyzzyaace1(is).and..not.xyzzyaace1(js))then
call dcopy(size_farray,xyzzyaacc1(1,1,1,is),1,xyzzyaacc1(1,1,1,js),1)
call dcopy(size_farray*ndet,xyzzyaacd1(1,1,1,1,is),1,xyzzyaacd1(1,1,1,&
&1,js),1)
xyzzyaace1(js)=.true.
else
xyzzyaace1(js)=.false.
endif
if(xyzzyaacg1(is).and..not.xyzzyaacg1(js))then
call dcopy(size_harray,xyzzyaacf1(1,1,1,1,1,is),1,xyzzyaacf1(1,1,1,1,1&
&,js),1)
xyzzyaacg1(js)=.true.
else
xyzzyaacg1(js)=.false.
endif
endif
if(scr_tasks(iwfn_detail,js))then
if(xyzzyaaay1(is).and..not.xyzzyaaay1(js))then
call dcopy(size_dbar,xyzzyaaax1(1,1,1,1,1,is),1,xyzzyaaax1(1,1,1,1,1,j&
&s),1)
xyzzyaaay1(js)=.true.
xyzzyaaaw1(js)=0
else
xyzzyaaay1(js)=.false.
xyzzyaaaw1(js)=0
endif
endif
call timer('CLONE_SCRATCH_SLATER',.false.)
end subroutine clone_scratch_slater
subroutine wfn_logval_slater(is,logwfn,iszero)
implicit none
integer,intent(in) :: is
complex(dp),intent(out) :: logwfn
logical,intent(out) :: iszero
call timer('WFN_LOGVAL_SLATER',.true.)
call xyzzyaaev1(is,.false.,.false.)
logwfn=xyzzyaaap1(is)+xyzzyaaal1(is)
iszero=xyzzyaaan1(is)/=xyzzyaaae1
call timer('WFN_LOGVAL_SLATER',.false.)
end subroutine wfn_logval_slater
subroutine wfn_loggrad_slater(ii,is,ilevel,val,sd,fidet,isnan,isinf)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: fidet(3)
logical,intent(in) :: val,sd
logical,intent(out) :: isnan,isinf
call timer('WFN_LOGGRAD_SLATER',.true.)
isnan=.false.
isinf=.false.
fidet=czero
if(xyzzyaabj1(ii,is))then
if(complex_wf)then
fidet=cmplx(xyzzyaabi1(1:3,ii,1,is),xyzzyaabi1(1:3,ii,2,is),dp)
else
fidet=cmplx(xyzzyaabi1(1:3,ii,1,is),0.d0,dp)
endif
else
if(.not.use_backflow)then
call xyzzyaafg1(ii,is,val)
call xyzzyaaev1(is,.false.,.false.)
call xyzzyaage1(ii,xyzzyaaap1(is),xyzzyaaan1(is),xyzzyaabm1(1,1,1,1,1,&
&is),xyzzyaaai1(1,is),fidet)
else
if(sd)then
call xyzzyaafk1(is)
else
call xyzzyaafj1(ii,is)
endif
if(xyzzyaaan1(is)==xyzzyaaae1)call loggrad_bf(ii,is,sd,xyzzyaacc1(1,1,&
&1,is),fidet)
endif
xyzzyaabi1(1:3,ii,1,is)=dble(fidet(1:3))
if(complex_wf)xyzzyaabi1(1:3,ii,2,is)=aimag(fidet(1:3))
xyzzyaabj1(ii,is)=.true.
endif
if(any(xyzzyaaai1(:,is)==xyzzyaaaf1))then
isnan=.true.
elseif(xyzzyaaan1(is)/=xyzzyaaae1)then
if(all(xyzzyaaai1(:,is)/=xyzzyaaae1))then
isnan=.true.
else
isinf=.true.
endif
endif
call timer('WFN_LOGGRAD_SLATER',.false.)
end subroutine wfn_loggrad_slater
subroutine wfn_loglap_slater(ii,is,ilevel,val,lapdet,isnan,isinf)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: lapdet
logical,intent(in) :: val
logical,intent(out) :: isnan,isinf
complex(dp) xyzzyaaaa15(3)
logical xyzzyaaab15,xyzzyaaac15
call timer('WFN_LOGLAP_SLATER',.true.)
lapdet=czero
if(xyzzyaabl1(ii,is))then
if(complex_wf)then
lapdet=cmplx(xyzzyaabk1(ii,1,is),xyzzyaabk1(ii,2,is),dp)
else
lapdet=cmplx(xyzzyaabk1(ii,1,is),0.d0,dp)
endif
else
if(.not.xyzzyaabj1(ii,is))then
call wfn_loggrad_slater(ii,is,ilevel,val,.true.,xyzzyaaaa15,xyzzyaaab1&
&5,xyzzyaaac15)
else
if(complex_wf)then
xyzzyaaaa15=cmplx(xyzzyaabi1(1:3,ii,1,is),xyzzyaabi1(1:3,ii,2,is),dp)
else
xyzzyaaaa15=cmplx(xyzzyaabi1(1:3,ii,1,is),0.d0,dp)
endif
endif
if(.not.use_backflow)then
call xyzzyaafh1(ii,is,val)
call xyzzyaaev1(is,.false.,.false.)
call xyzzyaagg1(ii,xyzzyaaap1(is),xyzzyaaan1(is),xyzzyaabo1(1,1,1,1,is&
&),xyzzyaaai1(1,is),xyzzyaaaa15,lapdet)
else
call xyzzyaafk1(is)
if(xyzzyaaan1(is)==xyzzyaaae1)call loglap_bf(ii,is,xyzzyaacc1(1,1,1,is&
&),xyzzyaacf1(1,1,1,1,1,is),xyzzyaaaa15,lapdet)
endif
xyzzyaabk1(ii,1,is)=dble(lapdet)
if(complex_wf)xyzzyaabk1(ii,2,is)=aimag(lapdet)
xyzzyaabl1(ii,is)=.true.
endif
if(any(xyzzyaaai1(:,is)==xyzzyaaaf1))then
isnan=.true.
elseif(xyzzyaaan1(is)/=xyzzyaaae1)then
if(all(xyzzyaaai1(:,is)/=xyzzyaaae1))then
isnan=.true.
else
isinf=.true.
endif
endif
call timer('WFN_LOGLAP_SLATER',.false.)
end subroutine wfn_loglap_slater
subroutine prefetch_wfn_slater(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
call timer('PREFETCH_WFN_SLATER',.true.)
if(.not.use_backflow)then
continue
else
if(sd)then
call xyzzyaafk1(is)
elseif(fd)then
call xyzzyaafi1(is)
endif
endif
call timer('PREFETCH_WFN_SLATER',.false.)
end subroutine prefetch_wfn_slater
subroutine xyzzyaaev1(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
if(xyzzyaaaq1(is))return
call timer('GET_LOGWFN',.true.)
call xyzzyaaew1(is,fd,sd)
if(use_altsamp.and.altsamp==2.and..not.isitcomplex)then
call xyzzyaafn1(xyzzyaaak1(1,is),xyzzyaaai1(1,is),xyzzyaaao1(1,is),xyz&
&zyaaap1(is),xyzzyaaan1(is))
else
call xyzzyaafm1(xyzzyaaak1(1,is),xyzzyaaai1(1,is),xyzzyaaao1(1,is),xyz&
&zyaaap1(is),xyzzyaaan1(is))
endif
xyzzyaaaq1(is)=.true.
call timer('GET_LOGWFN',.false.)
end subroutine xyzzyaaev1
recursive subroutine xyzzyaaew1(xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18)
implicit none
integer,intent(in) :: xyzzyaaaa18
logical,intent(in) :: xyzzyaaab18,xyzzyaaac18
integer xyzzyaaad18
if(xyzzyaaam1(xyzzyaaaa18))return
call timer('GET_LOGDET',.true.)
xyzzyaaad18=buffer_move1_from(xyzzyaaaa18)
if(xyzzyaaad18/=0)then
if(.not.use_backflow)then
call xyzzyaaex1(xyzzyaaad18,xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18)
else
if(bf_sparse)then
call xyzzyaafa1(xyzzyaaad18,xyzzyaaaa18,xyzzyaaab18)
else
call xyzzyaaff1(xyzzyaaaa18,.true.,xyzzyaaab18,xyzzyaaac18)
endif
endif
else
if(.not.use_backflow)then
if(in_range(xyzzyaaaa18,xyzzyaadv1))then
call xyzzyaafe1(xyzzyaaaa18,.true.,xyzzyaaab18.or.xyzzyaaac18)
else
call xyzzyaafc1(xyzzyaaaa18)
endif
else
call xyzzyaaff1(xyzzyaaaa18,.true.,xyzzyaaab18,xyzzyaaac18)
endif
endif
call timer('GET_LOGDET',.false.)
end subroutine xyzzyaaew1
recursive subroutine xyzzyaaex1(xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xy&
&zzyaaad19)
implicit none
integer,intent(in) :: xyzzyaaaa19,xyzzyaaab19
logical,intent(in) :: xyzzyaaac19,xyzzyaaad19
integer xyzzyaaae19,xyzzyaaaf19,xyzzyaaag19,xyzzyaaah19,xyzzyaaai19
if(xyzzyaaas1(xyzzyaaab19))return
call timer('GET_Q',.true.)
call xyzzyaaew1(xyzzyaaaa19,.false.,.false.)
call xyzzyaafc1(xyzzyaaaa19)
call xyzzyaaey1(xyzzyaaab19,xyzzyaaac19.or.xyzzyaaad19)
xyzzyaaae19=buffer_move1_from_ii(xyzzyaaab19)
xyzzyaaaf19=which_ie(xyzzyaaae19)
xyzzyaaag19=which_spin(xyzzyaaae19)
xyzzyaaaj1(:,detstart:detstop,xyzzyaaab19)=xyzzyaaaj1(:,detstart:detst&
&op,xyzzyaaaa19)
do xyzzyaaai19=detstart,detstop
xyzzyaaah19=upd_spin(xyzzyaaag19,xyzzyaaai19)
if(xyzzyaaah1(xyzzyaaah19,xyzzyaaai19,xyzzyaaaa19)/=xyzzyaaae1)then
call get_rsele(xyzzyaaab19)
if(pairing_wf)call get_eevecs(xyzzyaaab19)
call xyzzyaagb1(xyzzyaaah19,xyzzyaaai19,rele_scr(1,1,xyzzyaaab19),sele&
&_scr(1,xyzzyaaab19),xyzzyaaaj1(xyzzyaaah19,xyzzyaaai19,xyzzyaaab19),x&
&yzzyaaah1(xyzzyaaah19,xyzzyaaai19,xyzzyaaab19),eevecs_scr(1,1,1,xyzzy&
&aaab19))
cycle
endif
call xyzzyaafo1(xyzzyaaaf19,xyzzyaaag19,xyzzyaaai19,xyzzyaaaz1(1,1,1,x&
&yzzyaaah19,xyzzyaaai19,xyzzyaaaa19),xyzzyaabs1(1,1,xyzzyaaai19,xyzzya&
&aab19),xyzzyaaar1(xyzzyaaai19,xyzzyaaab19),xyzzyaaah1(xyzzyaaah19,xyz&
&zyaaai19,xyzzyaaab19),xyzzyaabq1(xyzzyaaai19,xyzzyaaab19),xyzzyaabr1(&
&1,xyzzyaaai19,xyzzyaaab19))
if(xyzzyaaah1(xyzzyaaah19,xyzzyaaai19,xyzzyaaab19)/=xyzzyaaae1)then
xyzzyaaaj1(xyzzyaaah19,xyzzyaaai19,xyzzyaaab19)=czero
cycle
endif
xyzzyaaaj1(xyzzyaaah19,xyzzyaaai19,xyzzyaaab19)=log(xyzzyaaar1(xyzzyaa&
&ai19,xyzzyaaab19))+xyzzyaaaj1(xyzzyaaah19,xyzzyaaai19,xyzzyaaab19)
enddo
xyzzyaaas1(xyzzyaaab19)=.true.
call xyzzyaafl1(xyzzyaaaj1(1,1,xyzzyaaab19),xyzzyaaah1(1,1,xyzzyaaab19&
&),xyzzyaaak1(1,xyzzyaaab19),xyzzyaaal1(xyzzyaaab19),xyzzyaaai1(1,xyzz&
&yaaab19))
xyzzyaaam1(xyzzyaaab19)=.true.
call timer('GET_Q',.false.)
end subroutine xyzzyaaex1
subroutine xyzzyaaey1(is,fsd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fsd
integer xyzzyaaaa20,xyzzyaaab20
logical xyzzyaaac20,xyzzyaaad20,xyzzyaaae20
xyzzyaaaa20=buffer_move1_from_ii(is)
xyzzyaaab20=which_spin(xyzzyaaaa20)
xyzzyaaac20=.not.xyzzyaabt1(is)
xyzzyaaad20=fsd.and.(xyzzyaabx1/=is.or.xyzzyaaby1/=xyzzyaaaa20)
xyzzyaaae20=xyzzyaaac20.or.xyzzyaaad20
if(.not.xyzzyaaae20)return
call timer('GET_SMAT1',.true.)
if(pairing_wf)call get_eevecs1_ch(xyzzyaaaa20,is)
call push_rion_hack(is)
call wfdet(rele1_chscr(1,is),sele1_chscr(is),xyzzyaaab20,wfdet_norb,wf&
&det_orbmask(1,xyzzyaaab20),xyzzyaaac20,xyzzyaaad20,xyzzyaadc1,xyzzyaa&
&dd1,xyzzyaade1,eevecs1_chscr(1,1,is),orb_m=xyzzyaacs1,orb_rmap=xyzzya&
&acu1)
if(sparse)call copy_rmap_orb_to_det(xyzzyaaab20,wfdet_orbmap,xyzzyaacs&
&1,xyzzyaacu1,xyzzyaabq1(1,is),xyzzyaabr1(1,1,is))
if(xyzzyaaac20)call copy_orb_to_det(1,xyzzyaaab20,wfdet_orbmap,xyzzyaa&
&dc1,xyzzyaabs1(1,1,1,is))
if(xyzzyaaad20)then
call copy_orb_to_det(3,xyzzyaaab20,wfdet_orbmap,xyzzyaadd1,xyzzyaaca1)
call copy_orb_to_det(1,xyzzyaaab20,wfdet_orbmap,xyzzyaade1,xyzzyaacb1)
endif
if(xyzzyaaac20)xyzzyaabt1(is)=.true.
if(xyzzyaaad20)then
xyzzyaabx1=is
xyzzyaaby1=xyzzyaaaa20
endif
call timer('GET_SMAT1',.false.)
end subroutine xyzzyaaey1
subroutine xyzzyaaez1(ii,ispin,is,val,fsd)
implicit none
integer,intent(in) :: ii,ispin,is
logical,intent(in) :: val,fsd
logical xyzzyaaaa21,xyzzyaaab21
call timer('GET_SMAT1',.true.)
xyzzyaaaa21=val.and.(xyzzyaabu1/=is.or.xyzzyaabv1/=ii)
if(xyzzyaaaa21)then
if(xyzzyaabt1(is).and.ii==buffer_move1_from_ii(is))then
call dcopy(size_rpsi,xyzzyaabs1(1,1,1,is),1,xyzzyaabw1,1)
xyzzyaabu1=is
xyzzyaabv1=ii
xyzzyaaaa21=.false.
elseif(xyzzyaaay1(is))then
call xyzzyaafu1(ii,xyzzyaabw1,xyzzyaaax1(1,1,1,1,1,is),back=.true.)
xyzzyaabu1=is
xyzzyaabv1=ii
xyzzyaaaa21=.false.
endif
endif
xyzzyaaab21=fsd.and.(xyzzyaabx1/=is.or.xyzzyaaby1/=ii)
if(xyzzyaaab21.and.xyzzyaabd1(is).and.xyzzyaabf1(is))then
call xyzzyaafv1(ii,xyzzyaaca1,xyzzyaabc1(1,1,1,1,1,1,is),back=.true.)
call xyzzyaafw1(ii,xyzzyaacb1,xyzzyaabe1(1,1,1,1,1,is),back=.true.)
xyzzyaabx1=is
xyzzyaaby1=ii
xyzzyaaab21=.false.
endif
if(xyzzyaaaa21.or.xyzzyaaab21)then
call get_rsele(is)
if(pairing_wf)call get_eevecs(is)
call push_rion_hack(is)
call wfdet(rele_scr(1,ii,is),sele_scr(ii,is),ispin,wfdet_norb,wfdet_or&
&bmask(1,ispin),xyzzyaaaa21,xyzzyaaab21,xyzzyaadc1,xyzzyaadd1,xyzzyaad&
&e1,eevecs_scr(1,1,ii,is))
if(xyzzyaaaa21)call copy_orb_to_det(1,ispin,wfdet_orbmap,xyzzyaadc1,xy&
&zzyaabw1)
if(xyzzyaaab21)then
call copy_orb_to_det(3,ispin,wfdet_orbmap,xyzzyaadd1,xyzzyaaca1)
call copy_orb_to_det(1,ispin,wfdet_orbmap,xyzzyaade1,xyzzyaacb1)
endif
if(xyzzyaaaa21)then
xyzzyaabu1=is
xyzzyaabv1=ii
endif
if(xyzzyaaab21)then
xyzzyaabx1=is
xyzzyaaby1=ii
endif
endif
call timer('GET_SMAT1',.false.)
end subroutine xyzzyaaez1
recursive subroutine xyzzyaafa1(xyzzyaaaa22,xyzzyaaab22,xyzzyaaac22)
implicit none
integer,intent(in) :: xyzzyaaaa22,xyzzyaaab22
logical,intent(in) :: xyzzyaaac22
integer xyzzyaaad22,xyzzyaaae22,xyzzyaaaf22
complex(dp) xyzzyaaag22,xyzzyaaah22
if(xyzzyaaav1(xyzzyaaab22))return
call timer('GET_SBAR',.true.)
call xyzzyaaew1(xyzzyaaaa22,.false.,.false.)
call xyzzyaafd1(xyzzyaaaa22)
call xyzzyaafb1(xyzzyaaab22,.true.,xyzzyaaac22)
xyzzyaaaj1(:,:,xyzzyaaab22)=xyzzyaaaj1(:,:,xyzzyaaaa22)
xyzzyaaah1(:,:,xyzzyaaab22)=xyzzyaaah1(:,:,xyzzyaaaa22)
if(pairing_wf)then
do xyzzyaaae22=detstart,detstop
do xyzzyaaad22=1,nspin
if(update_by_column(xyzzyaaad22,xyzzyaaae22))cycle
xyzzyaaaf22=upd_spin(xyzzyaaad22,xyzzyaaae22)
if(xyzzyaaah1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaaa22)/=xyzzyaaae1)then
call get_bf_x(xyzzyaaab22,.true.,.false.,.false.)
call get_eevecs_bf(xyzzyaaab22)
call xyzzyaagb1(xyzzyaaaf22,xyzzyaaae22,bf_x_scr(1,1,xyzzyaaab22),sele&
&_scr(1,xyzzyaaab22),xyzzyaaaj1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22),x&
&yzzyaaah1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22),eevecs_bf_scr(1,1,1,xy&
&zzyaaab22))
cycle
endif
call dcopy(size_onedbar,xyzzyaaaz1(1,1,1,xyzzyaaaf22,xyzzyaaae22,xyzzy&
&aaaa22),1,xyzzyaaaz1(1,1,1,xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22),1)
call xyzzyaafp1(bf_m_chscr(xyzzyaaad22,xyzzyaaab22),bf_rmap_chscr(1,xy&
&zzyaaad22,xyzzyaaab22),xyzzyaaad22,xyzzyaaae22,xyzzyaaaz1(1,1,1,xyzzy&
&aaaf22,xyzzyaaae22,xyzzyaaab22),xyzzyaach1(1,1,1,1,1,xyzzyaaab22),xyz&
&zyaaau1(1,1,1,xyzzyaaad22,xyzzyaaae22,xyzzyaaab22),xyzzyaaat1(1,xyzzy&
&aaad22,xyzzyaaae22,xyzzyaaab22),xyzzyaaag22,xyzzyaaah1(xyzzyaaaf22,xy&
&zzyaaae22,xyzzyaaab22))
if(xyzzyaaah1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22)/=xyzzyaaae1)then
call get_bf_x(xyzzyaaab22,.true.,.false.,.false.)
call get_eevecs_bf(xyzzyaaab22)
call xyzzyaagb1(xyzzyaaaf22,xyzzyaaae22,bf_x_scr(1,1,xyzzyaaab22),sele&
&_scr(1,xyzzyaaab22),xyzzyaaaj1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22),x&
&yzzyaaah1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22),eevecs_bf_scr(1,1,1,xy&
&zzyaaab22))
endif
call xyzzyaafr1(bf_m_chscr(xyzzyaaad22,xyzzyaaab22),bf_rmap_chscr(1,xy&
&zzyaaad22,xyzzyaaab22),xyzzyaaad22,xyzzyaaae22,xyzzyaaaz1(1,1,1,xyzzy&
&aaaf22,xyzzyaaae22,xyzzyaaab22),xyzzyaach1(1,1,1,1,1,xyzzyaaab22),xyz&
&zyaaau1(1,1,1,xyzzyaaad22,xyzzyaaae22,xyzzyaaab22),xyzzyaaat1(1,xyzzy&
&aaad22,xyzzyaaae22,xyzzyaaab22))
call xyzzyaafp1(bf_m_chscr(xyzzyaaaf22,xyzzyaaab22),bf_rmap_chscr(1,xy&
&zzyaaaf22,xyzzyaaab22),xyzzyaaaf22,xyzzyaaae22,xyzzyaaaz1(1,1,1,xyzzy&
&aaaf22,xyzzyaaae22,xyzzyaaab22),xyzzyaach1(1,1,1,1,1,xyzzyaaab22),xyz&
&zyaaau1(1,1,1,xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22),xyzzyaaat1(1,xyzzy&
&aaaf22,xyzzyaaae22,xyzzyaaab22),xyzzyaaah22,xyzzyaaah1(xyzzyaaaf22,xy&
&zzyaaae22,xyzzyaaab22))
xyzzyaaah22=xyzzyaaah22*xyzzyaaag22
if(xyzzyaaah1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22)/=xyzzyaaae1)then
xyzzyaaaj1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22)=czero
cycle
endif
xyzzyaaaj1(xyzzyaaaf22,xyzzyaaae22,xyzzyaaab22)=xyzzyaaaj1(xyzzyaaaf22&
&,xyzzyaaae22,xyzzyaaab22)+log(xyzzyaaah22)
enddo
enddo
do xyzzyaaae22=detstart,detstop
do xyzzyaaad22=1,nspin
if(count(upd_spin(:,xyzzyaaae22)==xyzzyaaad22)/=1)cycle
if(xyzzyaaah1(xyzzyaaad22,xyzzyaaae22,xyzzyaaaa22)/=xyzzyaaae1)then
call get_bf_x(xyzzyaaab22,.true.,.false.,.false.)
call get_eevecs_bf(xyzzyaaab22)
call xyzzyaagb1(xyzzyaaad22,xyzzyaaae22,bf_x_scr(1,1,xyzzyaaab22),sele&
&_scr(1,xyzzyaaab22),xyzzyaaaj1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22),x&
&yzzyaaah1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22),eevecs_bf_scr(1,1,1,xy&
&zzyaaab22))
cycle
endif
call xyzzyaafp1(bf_m_chscr(xyzzyaaad22,xyzzyaaab22),bf_rmap_chscr(1,xy&
&zzyaaad22,xyzzyaaab22),xyzzyaaad22,xyzzyaaae22,xyzzyaaaz1(1,1,1,xyzzy&
&aaad22,xyzzyaaae22,xyzzyaaaa22),xyzzyaach1(1,1,1,1,1,xyzzyaaab22),xyz&
&zyaaau1(1,1,1,xyzzyaaad22,xyzzyaaae22,xyzzyaaab22),xyzzyaaat1(1,xyzzy&
&aaad22,xyzzyaaae22,xyzzyaaab22),xyzzyaaah22,xyzzyaaah1(xyzzyaaad22,xy&
&zzyaaae22,xyzzyaaab22))
if(xyzzyaaah1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22)/=xyzzyaaae1)then
xyzzyaaaj1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22)=czero
cycle
endif
xyzzyaaaj1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22)=xyzzyaaaj1(xyzzyaaad22&
&,xyzzyaaae22,xyzzyaaab22)+log(xyzzyaaah22)
enddo
enddo
else
do xyzzyaaae22=detstart,detstop
do xyzzyaaad22=1,nspin
if(xyzzyaaah1(xyzzyaaad22,xyzzyaaae22,xyzzyaaaa22)/=xyzzyaaae1)then
call get_bf_x(xyzzyaaab22,.true.,.false.,.false.)
call xyzzyaagb1(xyzzyaaad22,xyzzyaaae22,bf_x_scr(1,1,xyzzyaaab22),sele&
&_scr(1,xyzzyaaab22),xyzzyaaaj1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22),x&
&yzzyaaah1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22),eevecs_scr(1,1,1,xyzzy&
&aaab22))
cycle
endif
call xyzzyaafp1(bf_m_chscr(xyzzyaaad22,xyzzyaaab22),bf_rmap_chscr(1,xy&
&zzyaaad22,xyzzyaaab22),xyzzyaaad22,xyzzyaaae22,xyzzyaaaz1(1,1,1,xyzzy&
&aaad22,xyzzyaaae22,xyzzyaaaa22),xyzzyaach1(1,1,1,1,1,xyzzyaaab22),xyz&
&zyaaau1(1,1,1,xyzzyaaad22,xyzzyaaae22,xyzzyaaab22),xyzzyaaat1(1,xyzzy&
&aaad22,xyzzyaaae22,xyzzyaaab22),xyzzyaaah22,xyzzyaaah1(xyzzyaaad22,xy&
&zzyaaae22,xyzzyaaab22))
if(xyzzyaaah1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22)/=xyzzyaaae1)then
xyzzyaaaj1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22)=czero
cycle
endif
xyzzyaaaj1(xyzzyaaad22,xyzzyaaae22,xyzzyaaab22)=xyzzyaaaj1(xyzzyaaad22&
&,xyzzyaaae22,xyzzyaaab22)+log(xyzzyaaah22)
enddo
enddo
endif
xyzzyaaav1(xyzzyaaab22)=.true.
call xyzzyaafl1(xyzzyaaaj1(1,1,xyzzyaaab22),xyzzyaaah1(1,1,xyzzyaaab22&
&),xyzzyaaak1(1,xyzzyaaab22),xyzzyaaal1(xyzzyaaab22),xyzzyaaai1(1,xyzz&
&yaaab22))
xyzzyaaam1(xyzzyaaab22)=.true.
call timer('GET_SBAR',.false.)
end subroutine xyzzyaafa1
subroutine xyzzyaafb1(is,val,fd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: val,fd
real(dp),pointer :: xyzzyaaaa23(:,:,:)
logical xyzzyaaab23,xyzzyaaac23,xyzzyaaad23
xyzzyaaab23=val.and..not.xyzzyaaci1(is)
xyzzyaaac23=fd.and..not.xyzzyaack1(is)
xyzzyaaad23=xyzzyaaab23.or.xyzzyaaac23
if(.not.xyzzyaaad23)return
call timer('GET_SCHUNK',.true.)
call get_bf_x(is,.true.,fd,.false.)
if(noncoll_spin)then
sele_scr(:,is)=sele_scr(:,buffer_move1_from(is))
sele_scr(buffer_move1_from_ii(is),is)=sele1_chscr(is)
endif
if(pairing_wf)then
call get_eevecs_bf(is)
xyzzyaaaa23=>eevecs_bf_scr(:,:,:,is)
else
xyzzyaaaa23=>eevecs_scr(:,:,:,is)
endif
call push_rion_hack(is)
if(xyzzyaaab23.and.xyzzyaaac23)then
call xyzzyaaft1(bf_x_scr(1,1,is),sele_scr(1,is),bf_m_chscr(1,is),bf_rm&
&ap_chscr(1,1,is),val,fd,xyzzyaaaa23,smat=xyzzyaach1(1,1,1,1,1,is),dsm&
&at=xyzzyaacj1(1,1,1,1,1,1,is),untwist=.true.)
elseif(xyzzyaaac23)then
call xyzzyaaft1(bf_x_scr(1,1,is),sele_scr(1,is),bf_m_chscr(1,is),bf_rm&
&ap_chscr(1,1,is),val,fd,xyzzyaaaa23,dsmat=xyzzyaacj1(1,1,1,1,1,1,is),&
&untwist=.true.)
else
call xyzzyaaft1(bf_x_scr(1,1,is),sele_scr(1,is),bf_m_chscr(1,is),bf_rm&
&ap_chscr(1,1,is),val,fd,xyzzyaaaa23,smat=xyzzyaach1(1,1,1,1,1,is),unt&
&wist=.true.)
endif
if(xyzzyaaab23)xyzzyaaci1(is)=.true.
if(xyzzyaaac23)xyzzyaack1(is)=.true.
call timer('GET_SCHUNK',.false.)
end subroutine xyzzyaafb1
recursive subroutine xyzzyaafc1(xyzzyaaaa24)
implicit none
integer,intent(in) :: xyzzyaaaa24
integer xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaae24,xyzzyaaaf24,xy&
&zzyaaag24
logical xyzzyaaah24
if(xyzzyaaba1(xyzzyaaaa24))return
call timer('GET_DBAR',.true.)
xyzzyaaah24=.false.
xyzzyaaac24=buffer_move1_from(xyzzyaaaa24)
if(xyzzyaaac24/=0)xyzzyaaah24=xyzzyaaba1(xyzzyaaac24).and.dbar_age(xyz&
&zyaaac24)<xyzzyaaaa1
if(xyzzyaaah24)then
call xyzzyaaex1(xyzzyaaac24,xyzzyaaaa24,.true.,.false.)
xyzzyaaah24=.not.any(xyzzyaaah1(:,:,xyzzyaaac24)/=xyzzyaaae1.and.xyzzy&
&aaah1(:,:,xyzzyaaaa24)==xyzzyaaae1)
endif
if(xyzzyaaah24)then
call xyzzyaaex1(xyzzyaaac24,xyzzyaaaa24,.true.,.false.)
call dcopy(size_dbar,xyzzyaaaz1(1,1,1,1,1,xyzzyaaac24),1,xyzzyaaaz1(1,&
&1,1,1,1,xyzzyaaaa24),1)
xyzzyaaab24=buffer_move1_from_ii(xyzzyaaaa24)
xyzzyaaad24=which_ie(xyzzyaaab24)
xyzzyaaae24=which_spin(xyzzyaaab24)
do xyzzyaaaf24=detstart,detstop
xyzzyaaag24=upd_spin(xyzzyaaae24,xyzzyaaaf24)
if(xyzzyaaah1(xyzzyaaag24,xyzzyaaaf24,xyzzyaaac24)/=xyzzyaaae1)cycle
if(xyzzyaaah1(xyzzyaaag24,xyzzyaaaf24,xyzzyaaaa24)/=xyzzyaaae1)cycle
call xyzzyaafq1(xyzzyaaad24,xyzzyaaae24,xyzzyaaaf24,xyzzyaabs1(1,1,xyz&
&zyaaaf24,xyzzyaaaa24),xyzzyaaaz1(1,1,1,xyzzyaaag24,xyzzyaaaf24,xyzzya&
&aaa24),xyzzyaaar1(xyzzyaaaf24,xyzzyaaaa24),xyzzyaabq1(xyzzyaaaf24,xyz&
&zyaaaa24),xyzzyaabr1(1,xyzzyaaaf24,xyzzyaaaa24))
enddo
xyzzyaaba1(xyzzyaaaa24)=.true.
dbar_age(xyzzyaaaa24)=dbar_age(xyzzyaaac24)+1
else
if(xyzzyaaay1(xyzzyaaaa24))then
call xyzzyaafz1(xyzzyaaax1(1,1,1,1,1,xyzzyaaaa24),xyzzyaaaz1(1,1,1,1,1&
&,xyzzyaaaa24),xyzzyaaaj1(1,1,xyzzyaaaa24),xyzzyaaah1(1,1,xyzzyaaaa24)&
&)
else
call get_rsele(xyzzyaaaa24)
if(pairing_wf)call get_eevecs(xyzzyaaaa24)
call push_rion_hack(xyzzyaaaa24)
call xyzzyaafy1(rele_scr(1,1,xyzzyaaaa24),sele_scr(1,xyzzyaaaa24),xyzz&
&yaaaz1(1,1,1,1,1,xyzzyaaaa24),xyzzyaaaj1(1,1,xyzzyaaaa24),xyzzyaaah1(&
&1,1,xyzzyaaaa24),eevecs_scr(1,1,1,xyzzyaaaa24))
endif
xyzzyaaba1(xyzzyaaaa24)=.true.
dbar_age(xyzzyaaaa24)=0
call xyzzyaafl1(xyzzyaaaj1(1,1,xyzzyaaaa24),xyzzyaaah1(1,1,xyzzyaaaa24&
&),xyzzyaaak1(1,xyzzyaaaa24),xyzzyaaal1(xyzzyaaaa24),xyzzyaaai1(1,xyzz&
&yaaaa24))
xyzzyaaam1(xyzzyaaaa24)=.true.
endif
call timer('GET_DBAR',.false.)
end subroutine xyzzyaafc1
subroutine xyzzyaafd1(is)
implicit none
integer,intent(in) :: is
integer xyzzyaaaa25,xyzzyaaab25,xyzzyaaac25,xyzzyaaad25
logical xyzzyaaae25
if(xyzzyaaba1(is))return
call timer('GET_DBAR_BF',.true.)
xyzzyaaae25=.false.
xyzzyaaaa25=buffer_move1_from(is)
if(bf_sparse.and.xyzzyaaaa25/=0)xyzzyaaae25=xyzzyaaba1(xyzzyaaaa25).an&
&d.dbar_age(xyzzyaaaa25)<xyzzyaaaa1
if(xyzzyaaae25)then
call xyzzyaafa1(xyzzyaaaa25,is,.false.)
xyzzyaaae25=.not.any(xyzzyaaah1(:,:,xyzzyaaaa25)/=xyzzyaaae1.and.xyzzy&
&aaah1(:,:,is)==xyzzyaaae1)
endif
if(xyzzyaaae25)then
if(pairing_wf)then
do xyzzyaaab25=detstart,detstop
do xyzzyaaac25=1,nspin
if(update_by_column(xyzzyaaac25,xyzzyaaab25))cycle
xyzzyaaad25=upd_spin(xyzzyaaac25,xyzzyaaab25)
if(xyzzyaaah1(xyzzyaaad25,xyzzyaaab25,xyzzyaaaa25)/=xyzzyaaae1)cycle
if(xyzzyaaah1(xyzzyaaad25,xyzzyaaab25,is)/=xyzzyaaae1)cycle
call xyzzyaafr1(bf_m_chscr(xyzzyaaad25,is),bf_rmap_chscr(1,xyzzyaaad25&
&,is),xyzzyaaad25,xyzzyaaab25,xyzzyaaaz1(1,1,1,xyzzyaaad25,xyzzyaaab25&
&,is),xyzzyaach1(1,1,1,1,1,is),xyzzyaaau1(1,1,1,xyzzyaaad25,xyzzyaaab2&
&5,is),xyzzyaaat1(1,xyzzyaaad25,xyzzyaaab25,is))
enddo
enddo
do xyzzyaaab25=detstart,detstop
do xyzzyaaac25=1,nspin
if(count(upd_spin(:,xyzzyaaab25)==xyzzyaaac25)/=1)cycle
if(xyzzyaaah1(xyzzyaaac25,xyzzyaaab25,xyzzyaaaa25)/=xyzzyaaae1)cycle
if(xyzzyaaah1(xyzzyaaac25,xyzzyaaab25,is)/=xyzzyaaae1)cycle
call dcopy(size_onedbar,xyzzyaaaz1(1,1,1,xyzzyaaac25,xyzzyaaab25,xyzzy&
&aaaa25),1,xyzzyaaaz1(1,1,1,xyzzyaaac25,xyzzyaaab25,is),1)
call xyzzyaafr1(bf_m_chscr(xyzzyaaac25,is),bf_rmap_chscr(1,xyzzyaaac25&
&,is),xyzzyaaac25,xyzzyaaab25,xyzzyaaaz1(1,1,1,xyzzyaaac25,xyzzyaaab25&
&,is),xyzzyaach1(1,1,1,1,1,is),xyzzyaaau1(1,1,1,xyzzyaaac25,xyzzyaaab2&
&5,is),xyzzyaaat1(1,xyzzyaaac25,xyzzyaaab25,is))
enddo
enddo
else
call dcopy(size_dbar,xyzzyaaaz1(1,1,1,1,1,xyzzyaaaa25),1,xyzzyaaaz1(1,&
&1,1,1,1,is),1)
do xyzzyaaab25=detstart,detstop
do xyzzyaaac25=1,nspin
if(xyzzyaaah1(xyzzyaaac25,xyzzyaaab25,xyzzyaaaa25)/=xyzzyaaae1)cycle
if(xyzzyaaah1(xyzzyaaac25,xyzzyaaab25,is)/=xyzzyaaae1)cycle
call xyzzyaafr1(bf_m_chscr(xyzzyaaac25,is),bf_rmap_chscr(1,xyzzyaaac25&
&,is),xyzzyaaac25,xyzzyaaab25,xyzzyaaaz1(1,1,1,xyzzyaaac25,xyzzyaaab25&
&,is),xyzzyaach1(1,1,1,1,1,is),xyzzyaaau1(1,1,1,xyzzyaaac25,xyzzyaaab2&
&5,is),xyzzyaaat1(1,xyzzyaaac25,xyzzyaaab25,is))
enddo
enddo
endif
xyzzyaaba1(is)=.true.
dbar_age(is)=dbar_age(xyzzyaaaa25)+1
else
call xyzzyaaff1(is,.true.,.false.,.false.)
call xyzzyaafz1(xyzzyaaax1(1,1,1,1,1,is),xyzzyaaaz1(1,1,1,1,1,is),xyzz&
&yaaaj1(1,1,is),xyzzyaaah1(1,1,is))
xyzzyaaba1(is)=.true.
dbar_age(is)=0
call xyzzyaafl1(xyzzyaaaj1(1,1,is),xyzzyaaah1(1,1,is),xyzzyaaak1(1,is)&
&,xyzzyaaal1(is),xyzzyaaai1(1,is))
xyzzyaaam1(is)=.true.
endif
call timer('GET_DBAR_BF',.false.)
end subroutine xyzzyaafd1
recursive subroutine xyzzyaafe1(xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26)
implicit none
integer,intent(in) :: xyzzyaaaa26
logical,intent(in) :: xyzzyaaab26,xyzzyaaac26
integer xyzzyaaad26,xyzzyaaae26,xyzzyaaaf26,xyzzyaaag26,xyzzyaaah26,xy&
&zzyaaai26,xyzzyaaaj26,xyzzyaaak26
logical xyzzyaaal26,xyzzyaaam26,xyzzyaaan26
call timer('GET_SMAT',.true.)
xyzzyaaal26=xyzzyaaab26.and..not.xyzzyaaay1(xyzzyaaaa26)
xyzzyaaam26=xyzzyaaac26.and..not.(xyzzyaabd1(xyzzyaaaa26).and.xyzzyaab&
&f1(xyzzyaaaa26))
xyzzyaaan26=xyzzyaaal26.or.xyzzyaaam26
if(xyzzyaaan26)then
xyzzyaaaj26=buffer_move1_from(xyzzyaaaa26)
xyzzyaaak26=buffer_move2_from(xyzzyaaaa26)
if(xyzzyaaaj26/=0.and..not.xyzzyaaac26)then
if(xyzzyaaaw1(xyzzyaaaa26)/=xyzzyaaaj26)then
call xyzzyaafe1(xyzzyaaaj26,.true.,.false.)
call dcopy(size_dbar,xyzzyaaax1(1,1,1,1,1,xyzzyaaaj26),1,xyzzyaaax1(1,&
&1,1,1,1,xyzzyaaaa26),1)
endif
call xyzzyaaey1(xyzzyaaaa26,.false.)
xyzzyaaaf26=buffer_move1_from_ii(xyzzyaaaa26)
call xyzzyaafu1(xyzzyaaaf26,xyzzyaabs1(1,1,1,xyzzyaaaa26),xyzzyaaax1(1&
&,1,1,1,1,xyzzyaaaa26))
xyzzyaaay1(xyzzyaaaa26)=.true.
xyzzyaaaw1(xyzzyaaaa26)=0
elseif(xyzzyaaak26/=0.and..not.xyzzyaaac26)then
xyzzyaaaf26=buffer_move2_from_ii(xyzzyaaaa26)
xyzzyaaah26=which_spin(xyzzyaaaf26)
xyzzyaaag26=buffer_move2_from_jj(xyzzyaaaa26)
xyzzyaaai26=which_spin(xyzzyaaag26)
if(xyzzyaaaw1(xyzzyaaaa26)/=xyzzyaaak26)then
call xyzzyaafe1(xyzzyaaak26,.true.,.false.)
call dcopy(size_dbar,xyzzyaaax1(1,1,1,1,1,xyzzyaaak26),1,xyzzyaaax1(1,&
&1,1,1,1,xyzzyaaaa26),1)
endif
call xyzzyaaez1(xyzzyaaaf26,xyzzyaaah26,xyzzyaaaa26,.true.,.false.)
call xyzzyaafu1(xyzzyaaaf26,xyzzyaabw1,xyzzyaaax1(1,1,1,1,1,xyzzyaaaa2&
&6))
call xyzzyaaez1(xyzzyaaag26,xyzzyaaai26,xyzzyaaaa26,.true.,.false.)
call xyzzyaafu1(xyzzyaaag26,xyzzyaabw1,xyzzyaaax1(1,1,1,1,1,xyzzyaaaa2&
&6))
xyzzyaaay1(xyzzyaaaa26)=.true.
xyzzyaaaw1(xyzzyaaaa26)=0
else
call get_rsele(xyzzyaaaa26)
if(pairing_wf)call get_eevecs(xyzzyaaaa26)
call push_rion_hack(xyzzyaaaa26)
if(xyzzyaaal26.and.xyzzyaaam26)then
call xyzzyaafs1(rele_scr(1,1,xyzzyaaaa26),sele_scr(1,xyzzyaaaa26),xyzz&
&yaaal26,xyzzyaaam26,xyzzyaaam26,eevecs_scr(1,1,1,xyzzyaaaa26),smat=xy&
&zzyaaax1(1,1,1,1,1,xyzzyaaaa26),dsmat=xyzzyaabc1(1,1,1,1,1,1,xyzzyaaa&
&a26),lapsmat=xyzzyaabe1(1,1,1,1,1,xyzzyaaaa26),logdet=xyzzyaaaj1(1,1,&
&xyzzyaaaa26),fpeinfo_det=xyzzyaaah1(1,1,xyzzyaaaa26))
elseif(xyzzyaaam26)then
call xyzzyaafs1(rele_scr(1,1,xyzzyaaaa26),sele_scr(1,xyzzyaaaa26),xyzz&
&yaaal26,xyzzyaaam26,xyzzyaaam26,eevecs_scr(1,1,1,xyzzyaaaa26),dsmat=x&
&yzzyaabc1(1,1,1,1,1,1,xyzzyaaaa26),lapsmat=xyzzyaabe1(1,1,1,1,1,xyzzy&
&aaaa26))
else
call xyzzyaafs1(rele_scr(1,1,xyzzyaaaa26),sele_scr(1,xyzzyaaaa26),xyzz&
&yaaal26,xyzzyaaam26,xyzzyaaam26,eevecs_scr(1,1,1,xyzzyaaaa26),smat=xy&
&zzyaaax1(1,1,1,1,1,xyzzyaaaa26),logdet=xyzzyaaaj1(1,1,xyzzyaaaa26),fp&
&einfo_det=xyzzyaaah1(1,1,xyzzyaaaa26))
endif
if(xyzzyaaal26)then
xyzzyaaay1(xyzzyaaaa26)=.true.
xyzzyaaaw1(xyzzyaaaa26)=0
call xyzzyaafl1(xyzzyaaaj1(1,1,xyzzyaaaa26),xyzzyaaah1(1,1,xyzzyaaaa26&
&),xyzzyaaak1(1,xyzzyaaaa26),xyzzyaaal1(xyzzyaaaa26),xyzzyaaai1(1,xyzz&
&yaaaa26))
xyzzyaaam1(xyzzyaaaa26)=.true.
endif
if(xyzzyaaam26)then
xyzzyaabd1(xyzzyaaaa26)=.true.
xyzzyaabb1(xyzzyaaaa26)=0
xyzzyaabf1(xyzzyaaaa26)=.true.
endif
endif
endif
if(xyzzyaaab26.and..not.xyzzyaaam1(xyzzyaaaa26))then
xyzzyaaaj1(:,:,xyzzyaaaa26)=czero
xyzzyaaah1(:,:,xyzzyaaaa26)=xyzzyaaae1
do xyzzyaaae26=detstart,detstop
do xyzzyaaad26=1,nspin
if(missing_det(xyzzyaaad26,xyzzyaaae26))cycle
call xyzzyaaga1(xyzzyaaax1(1,1,1,1,1,xyzzyaaaa26),xyzzyaaad26,xyzzyaaa&
&e26,xyzzyaaaj1(xyzzyaaad26,xyzzyaaae26,xyzzyaaaa26),xyzzyaaah1(xyzzya&
&aad26,xyzzyaaae26,xyzzyaaaa26))
enddo
enddo
call xyzzyaafl1(xyzzyaaaj1(1,1,xyzzyaaaa26),xyzzyaaah1(1,1,xyzzyaaaa26&
&),xyzzyaaak1(1,xyzzyaaaa26),xyzzyaaal1(xyzzyaaaa26),xyzzyaaai1(1,xyzz&
&yaaaa26))
xyzzyaaam1(xyzzyaaaa26)=.true.
endif
call timer('GET_SMAT',.false.)
end subroutine xyzzyaafe1
recursive subroutine xyzzyaaff1(xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xy&
&zzyaaad27)
implicit none
integer,intent(in) :: xyzzyaaaa27
logical,intent(in) :: xyzzyaaab27,xyzzyaaac27,xyzzyaaad27
integer xyzzyaaae27,xyzzyaaaf27,xyzzyaaag27,xyzzyaaah27,xyzzyaaai27,xy&
&zzyaaaj27
logical xyzzyaaak27,xyzzyaaal27,xyzzyaaam27,xyzzyaaan27,xyzzyaaao27
real(dp),pointer :: xyzzyaaap27(:,:,:)
call timer('GET_SMAT_BF',.true.)
xyzzyaaak27=xyzzyaaab27.and..not.xyzzyaaay1(xyzzyaaaa27)
xyzzyaaal27=xyzzyaaac27.and..not.xyzzyaabd1(xyzzyaaaa27)
xyzzyaaam27=xyzzyaaad27.and..not.xyzzyaabh1(xyzzyaaaa27)
xyzzyaaao27=xyzzyaaak27.or.xyzzyaaal27.or.xyzzyaaam27
if(xyzzyaaao27)then
xyzzyaaaf27=buffer_move1_from(xyzzyaaaa27)
if(xyzzyaaaf27/=0.and..not.xyzzyaaam27)then
if(xyzzyaaaw1(xyzzyaaaa27)/=xyzzyaaaf27)then
call xyzzyaaff1(xyzzyaaaf27,xyzzyaaak27,xyzzyaaal27,.false.)
if(xyzzyaaak27)call dcopy(size_dbar,xyzzyaaax1(1,1,1,1,1,xyzzyaaaf27),&
&1,xyzzyaaax1(1,1,1,1,1,xyzzyaaaa27),1)
if(xyzzyaaal27)call dcopy(size_dsmat,xyzzyaabc1(1,1,1,1,1,1,xyzzyaaaf2&
&7),1,xyzzyaabc1(1,1,1,1,1,1,xyzzyaaaa27),1)
endif
call xyzzyaafb1(xyzzyaaaa27,xyzzyaaak27,xyzzyaaal27)
do xyzzyaaai27=1,nspin
do xyzzyaaae27=1,bf_m_chscr(xyzzyaaai27,xyzzyaaaa27)
xyzzyaaah27=bf_rmap_chscr(xyzzyaaae27,xyzzyaaai27,xyzzyaaaa27)
xyzzyaaag27=which_ii(xyzzyaaah27,xyzzyaaai27)
if(xyzzyaaak27)call xyzzyaafu1(xyzzyaaag27,xyzzyaach1(1,1,1,xyzzyaaah2&
&7,xyzzyaaai27,xyzzyaaaa27),xyzzyaaax1(1,1,1,1,1,xyzzyaaaa27))
if(xyzzyaaal27)call xyzzyaafv1(xyzzyaaag27,xyzzyaacj1(1,1,1,1,xyzzyaaa&
&h27,xyzzyaaai27,xyzzyaaaa27),xyzzyaabc1(1,1,1,1,1,1,xyzzyaaaa27))
enddo
enddo
if(xyzzyaaak27)then
xyzzyaaay1(xyzzyaaaa27)=.true.
xyzzyaaaw1(xyzzyaaaa27)=0
endif
if(xyzzyaaal27)then
xyzzyaabd1(xyzzyaaaa27)=.true.
xyzzyaabb1(xyzzyaaaa27)=0
endif
else
call get_bf_x(xyzzyaaaa27,.true.,.false.,.false.)
if(pairing_wf)then
call get_eevecs_bf(xyzzyaaaa27)
xyzzyaaap27=>eevecs_bf_scr(:,:,:,xyzzyaaaa27)
else
xyzzyaaap27=>eevecs_scr(:,:,:,xyzzyaaaa27)
endif
xyzzyaaan27=xyzzyaaal27.or.xyzzyaaam27
call push_rion_hack(xyzzyaaaa27)
if(xyzzyaaak27.and.xyzzyaaan27)then
call xyzzyaafs1(bf_x_scr(1,1,xyzzyaaaa27),sele_scr(1,xyzzyaaaa27),xyzz&
&yaaak27,xyzzyaaal27,xyzzyaaam27,xyzzyaaap27,smat=xyzzyaaax1(1,1,1,1,1&
&,xyzzyaaaa27),dsmat=xyzzyaabc1(1,1,1,1,1,1,xyzzyaaaa27),d2smat=xyzzya&
&abg1(1,1,1,1,1,1,xyzzyaaaa27),logdet=xyzzyaaaj1(1,1,xyzzyaaaa27),fpei&
&nfo_det=xyzzyaaah1(1,1,xyzzyaaaa27))
elseif(xyzzyaaan27)then
call xyzzyaafs1(bf_x_scr(1,1,xyzzyaaaa27),sele_scr(1,xyzzyaaaa27),xyzz&
&yaaak27,xyzzyaaal27,xyzzyaaam27,xyzzyaaap27,dsmat=xyzzyaabc1(1,1,1,1,&
&1,1,xyzzyaaaa27),d2smat=xyzzyaabg1(1,1,1,1,1,1,xyzzyaaaa27))
else
call xyzzyaafs1(bf_x_scr(1,1,xyzzyaaaa27),sele_scr(1,xyzzyaaaa27),xyzz&
&yaaak27,xyzzyaaal27,xyzzyaaam27,xyzzyaaap27,smat=xyzzyaaax1(1,1,1,1,1&
&,xyzzyaaaa27),logdet=xyzzyaaaj1(1,1,xyzzyaaaa27),fpeinfo_det=xyzzyaaa&
&h1(1,1,xyzzyaaaa27))
endif
if(xyzzyaaak27)then
xyzzyaaay1(xyzzyaaaa27)=.true.
xyzzyaaaw1(xyzzyaaaa27)=0
call xyzzyaafl1(xyzzyaaaj1(1,1,xyzzyaaaa27),xyzzyaaah1(1,1,xyzzyaaaa27&
&),xyzzyaaak1(1,xyzzyaaaa27),xyzzyaaal1(xyzzyaaaa27),xyzzyaaai1(1,xyzz&
&yaaaa27))
xyzzyaaam1(xyzzyaaaa27)=.true.
endif
if(xyzzyaaal27)then
xyzzyaabd1(xyzzyaaaa27)=.true.
xyzzyaabb1(xyzzyaaaa27)=0
endif
if(xyzzyaaam27)xyzzyaabh1(xyzzyaaaa27)=.true.
endif
endif
if(xyzzyaaab27.and..not.xyzzyaaam1(xyzzyaaaa27))then
xyzzyaaaj1(:,:,xyzzyaaaa27)=czero
xyzzyaaah1(:,:,xyzzyaaaa27)=xyzzyaaae1
do xyzzyaaaj27=detstart,detstop
do xyzzyaaai27=1,nspin
if(missing_det(xyzzyaaai27,xyzzyaaaj27))cycle
call xyzzyaaga1(xyzzyaaax1(1,1,1,1,1,xyzzyaaaa27),xyzzyaaai27,xyzzyaaa&
&j27,xyzzyaaaj1(xyzzyaaai27,xyzzyaaaj27,xyzzyaaaa27),xyzzyaaah1(xyzzya&
&aai27,xyzzyaaaj27,xyzzyaaaa27))
enddo
enddo
call xyzzyaafl1(xyzzyaaaj1(1,1,xyzzyaaaa27),xyzzyaaah1(1,1,xyzzyaaaa27&
&),xyzzyaaak1(1,xyzzyaaaa27),xyzzyaaal1(xyzzyaaaa27),xyzzyaaai1(1,xyzz&
&yaaaa27))
xyzzyaaam1(xyzzyaaaa27)=.true.
endif
call timer('GET_SMAT_BF',.false.)
end subroutine xyzzyaaff1
subroutine xyzzyaafg1(ii,is,val)
implicit none
integer,intent(in) :: ii,is
logical,intent(in) :: val
integer xyzzyaaaa28,xyzzyaaab28,xyzzyaaac28
logical xyzzyaaad28
real(dp),pointer :: xyzzyaaae28(:,:,:,:)=>null()
if(xyzzyaabn1(ii,is))return
call timer('GET_FI_PROD_DET',.true.)
xyzzyaaab28=which_ie(ii)
xyzzyaaac28=which_spin(ii)
if(xyzzyaabd1(is).and..not.pairing_wf)then
xyzzyaaae28=>xyzzyaabc1(:,:,:,:,xyzzyaaab28,xyzzyaaac28,is)
else
call xyzzyaaez1(ii,xyzzyaaac28,is,val,.true.)
xyzzyaaae28=>xyzzyaaca1(:,:,:,:)
endif
xyzzyaaad28=.false.
if(.not.xyzzyaaba1(is))then
xyzzyaaaa28=buffer_move1_from(is)
if(xyzzyaaaa28/=0)xyzzyaaad28=buffer_move1_from_ii(is)==ii.and.xyzzyaa&
&ba1(xyzzyaaaa28).and.xyzzyaaas1(is).and.all(xyzzyaaah1(:,:,is)==xyzzy&
&aaae1)
if(.not.xyzzyaaad28)call xyzzyaafc1(is)
endif
if(.not.xyzzyaaad28)then
call xyzzyaagd1(xyzzyaaab28,xyzzyaaac28,xyzzyaaae28,xyzzyaaaz1(1,1,1,1&
&,1,is),xyzzyaadg1,xyzzyaaak1(1,is),xyzzyaabm1(1,1,1,1,1,is),xyzzyaaai&
&1(1,is))
else
call xyzzyaagd1(xyzzyaaab28,xyzzyaaac28,xyzzyaaae28,xyzzyaaaz1(1,1,1,1&
&,1,xyzzyaaaa28),xyzzyaaar1(1,is),xyzzyaaak1(1,is),xyzzyaabm1(1,1,1,1,&
&1,is),xyzzyaaai1(1,is))
endif
xyzzyaabn1(ii,is)=.true.
call timer('GET_FI_PROD_DET',.false.)
end subroutine xyzzyaafg1
subroutine xyzzyaafh1(ii,is,val)
implicit none
integer,intent(in) :: ii,is
logical,intent(in) :: val
integer xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29
real(dp),pointer :: xyzzyaaad29(:,:,:)=>null()
logical xyzzyaaae29
if(xyzzyaabp1(ii,is))return
call timer('GET_PROD_LAPDET',.true.)
xyzzyaaab29=which_ie(ii)
xyzzyaaac29=which_spin(ii)
if(xyzzyaabf1(is).and..not.pairing_wf)then
xyzzyaaad29=>xyzzyaabe1(:,:,:,xyzzyaaab29,xyzzyaaac29,is)
else
call xyzzyaaez1(ii,xyzzyaaac29,is,val,.true.)
xyzzyaaad29=>xyzzyaacb1(:,:,:)
endif
xyzzyaaae29=.false.
if(.not.xyzzyaaba1(is))then
xyzzyaaaa29=buffer_move1_from(is)
if(xyzzyaaaa29/=0)xyzzyaaae29=buffer_move1_from_ii(is)==ii.and.xyzzyaa&
&ba1(xyzzyaaaa29).and.xyzzyaaas1(is).and.all(xyzzyaaah1(:,:,is)==xyzzy&
&aaae1)
if(.not.xyzzyaaae29)call xyzzyaafc1(is)
endif
if(.not.xyzzyaaae29)then
call xyzzyaagf1(xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaaz1(1,1,1,1&
&,1,is),xyzzyaadg1,xyzzyaaak1(1,is),xyzzyaabo1(1,1,1,1,is),xyzzyaaai1(&
&1,is))
else
call xyzzyaagf1(xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaaz1(1,1,1,1&
&,1,xyzzyaaaa29),xyzzyaaar1(1,is),xyzzyaaak1(1,is),xyzzyaabo1(1,1,1,1,&
&is),xyzzyaaai1(1,is))
endif
xyzzyaabp1(ii,is)=.true.
call timer('GET_PROD_LAPDET',.false.)
end subroutine xyzzyaafh1
subroutine xyzzyaafi1(is)
implicit none
integer,intent(in) :: is
if(xyzzyaace1(is))return
call timer('GET_FARRAY',.true.)
call xyzzyaaff1(is,.true.,.true.,.false.)
call xyzzyaafd1(is)
call xyzzyaaev1(is,.false.,.false.)
call xyzzyaagh1(xyzzyaaaz1(1,1,1,1,1,is),xyzzyaaao1(1,is),xyzzyaaap1(i&
&s),xyzzyaabc1(1,1,1,1,1,1,is),xyzzyaacd1(1,1,1,1,is),xyzzyaacc1(1,1,1&
&,is),xyzzyaaai1(1,is))
xyzzyaace1(is)=.true.
call timer('GET_FARRAY',.false.)
end subroutine xyzzyaafi1
subroutine xyzzyaafj1(ii,is)
implicit none
integer,intent(in) :: ii,is
if(xyzzyaace1(is))return
call timer('GET_FARRAY_II',.true.)
call xyzzyaaff1(is,.true.,.true.,.false.)
call xyzzyaafd1(is)
call xyzzyaaev1(is,.false.,.false.)
call xyzzyaagi1(xyzzyaaaz1(1,1,1,1,1,is),xyzzyaaao1(1,is),xyzzyaaap1(i&
&s),xyzzyaabc1(1,1,1,1,1,1,is),bf_m_scr(ii,is),bf_rmap_scr(1,ii,is),xy&
&zzyaacd1(1,1,1,1,is),xyzzyaacc1(1,1,1,is),xyzzyaaai1(1,is))
call timer('GET_FARRAY_II',.false.)
end subroutine xyzzyaafj1
subroutine xyzzyaafk1(is)
implicit none
integer,intent(in) :: is
if(xyzzyaace1(is).and.xyzzyaacg1(is))return
call timer('GET_HARRAY',.true.)
call get_bf_x(is,.true.,.true.,.true.)
call xyzzyaaff1(is,.true.,.true.,.true.)
call xyzzyaafd1(is)
call xyzzyaaev1(is,.false.,.false.)
if(.not.xyzzyaace1(is))then
call xyzzyaagh1(xyzzyaaaz1(1,1,1,1,1,is),xyzzyaaao1(1,is),xyzzyaaap1(i&
&s),xyzzyaabc1(1,1,1,1,1,1,is),xyzzyaacd1(1,1,1,1,is),xyzzyaacc1(1,1,1&
&,is),xyzzyaaai1(1,is))
xyzzyaace1(is)=.true.
endif
call xyzzyaagj1(xyzzyaaaz1(1,1,1,1,1,is),xyzzyaaao1(1,is),xyzzyaaap1(i&
&s),xyzzyaabc1(1,1,1,1,1,1,is),xyzzyaabg1(1,1,1,1,1,1,is),bf_m2_scr(1,&
&is),bf_rmap2_scr(1,1,is),xyzzyaacd1(1,1,1,1,is),xyzzyaacf1(1,1,1,1,1,&
&is),xyzzyaaai1(1,is))
xyzzyaacg1(is)=.true.
call timer('GET_HARRAY',.false.)
end subroutine xyzzyaafk1
subroutine setup_slater_params(nparam)
implicit none
integer,intent(inout) :: nparam
integer xyzzyaaaa33,xyzzyaaab33,xyzzyaaac33,xyzzyaaad33
xyzzyaacn1=0
call setup_mdet_params(xyzzyaacn1(1))
call setup_wfdet_params(xyzzyaacn1(2))
call setup_bf_params(xyzzyaacn1(3))
xyzzyaacm1=sum(xyzzyaacn1,xyzzyaacn1>0)
nparam=xyzzyaacm1
allocate(xyzzyaaco1(xyzzyaacm1),stat=xyzzyaaad33)
call check_alloc(xyzzyaaad33,'SETUP_SLATER_PARAMS','slater_param_sec')
xyzzyaaac33=0
do xyzzyaaaa33=1,xyzzyaacl1
if(xyzzyaacn1(xyzzyaaaa33)<1)cycle
xyzzyaaab33=xyzzyaaac33+1
xyzzyaaac33=xyzzyaaac33+xyzzyaacn1(xyzzyaaaa33)
xyzzyaaco1(xyzzyaaab33:xyzzyaaac33)=xyzzyaaaa33
enddo
end subroutine setup_slater_params
subroutine finish_slater_params
implicit none
call finish_mdet_params
call finish_wfdet_params
call finish_bf_params
deallocate(xyzzyaaco1)
end subroutine finish_slater_params
subroutine get_slater_params(params,has_lolim,lolim,has_hilim,hilim,is&
&_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,la&
&bel)
implicit none
character(2),intent(inout) :: label(xyzzyaacm1)
real(dp),intent(inout) :: params(xyzzyaacm1),lolim(xyzzyaacm1),hilim(x&
&yzzyaacm1)
logical,intent(inout) :: has_lolim(xyzzyaacm1),has_hilim(xyzzyaacm1),i&
&s_shallow(xyzzyaacm1),is_redundant(xyzzyaacm1),is_linear(xyzzyaacm1),&
&is_loglinear(xyzzyaacm1),has_aderiv(xyzzyaacm1),affect_map(xyzzyaacm1&
&,xyzzyaacm1)
integer xyzzyaaaa35,xyzzyaaab35
label='S '
xyzzyaaab35=0
if(xyzzyaacn1(1)>0)then
xyzzyaaaa35=xyzzyaaab35+1
xyzzyaaab35=xyzzyaaab35+xyzzyaacn1(1)
call get_mdet_params(params(xyzzyaaaa35:xyzzyaaab35),has_lolim(xyzzyaa&
&aa35:xyzzyaaab35),lolim(xyzzyaaaa35:xyzzyaaab35),has_hilim(xyzzyaaaa3&
&5:xyzzyaaab35),hilim(xyzzyaaaa35:xyzzyaaab35),is_shallow(xyzzyaaaa35:&
&xyzzyaaab35),is_redundant(xyzzyaaaa35:xyzzyaaab35),is_linear(xyzzyaaa&
&a35:xyzzyaaab35),is_loglinear(xyzzyaaaa35:xyzzyaaab35),has_aderiv(xyz&
&zyaaaa35:xyzzyaaab35),affect_map(xyzzyaaaa35:xyzzyaaab35,xyzzyaaaa35:&
&xyzzyaaab35),label(xyzzyaaaa35:xyzzyaaab35))
endif
if(xyzzyaacn1(2)>0)then
xyzzyaaaa35=xyzzyaaab35+1
xyzzyaaab35=xyzzyaaab35+xyzzyaacn1(2)
call get_wfdet_params(params(xyzzyaaaa35:xyzzyaaab35),has_lolim(xyzzya&
&aaa35:xyzzyaaab35),lolim(xyzzyaaaa35:xyzzyaaab35),has_hilim(xyzzyaaaa&
&35:xyzzyaaab35),hilim(xyzzyaaaa35:xyzzyaaab35),is_shallow(xyzzyaaaa35&
&:xyzzyaaab35),is_redundant(xyzzyaaaa35:xyzzyaaab35),is_linear(xyzzyaa&
&aa35:xyzzyaaab35),is_loglinear(xyzzyaaaa35:xyzzyaaab35),has_aderiv(xy&
&zzyaaaa35:xyzzyaaab35),affect_map(xyzzyaaaa35:xyzzyaaab35,xyzzyaaaa35&
&:xyzzyaaab35),label(xyzzyaaaa35:xyzzyaaab35))
endif
if(xyzzyaacn1(3)>0)then
xyzzyaaaa35=xyzzyaaab35+1
xyzzyaaab35=xyzzyaaab35+xyzzyaacn1(3)
call get_bf_params(params(xyzzyaaaa35:xyzzyaaab35),has_lolim(xyzzyaaaa&
&35:xyzzyaaab35),lolim(xyzzyaaaa35:xyzzyaaab35),has_hilim(xyzzyaaaa35:&
&xyzzyaaab35),hilim(xyzzyaaaa35:xyzzyaaab35),is_shallow(xyzzyaaaa35:xy&
&zzyaaab35),is_redundant(xyzzyaaaa35:xyzzyaaab35),is_linear(xyzzyaaaa3&
&5:xyzzyaaab35),is_loglinear(xyzzyaaaa35:xyzzyaaab35),has_aderiv(xyzzy&
&aaaa35:xyzzyaaab35),affect_map(xyzzyaaaa35:xyzzyaaab35,xyzzyaaaa35:xy&
&zzyaaab35),label(xyzzyaaaa35:xyzzyaaab35))
endif
end subroutine get_slater_params
subroutine put_slater_params(params,ignore,iparam_buffer,prestore,bad_&
&params)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaacm1)
logical,intent(in) :: ignore(xyzzyaacm1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36,xyzzyaaad36
logical xyzzyaaae36
bad_params=.false.
xyzzyaaac36=0
xyzzyaaad36=0
if(iparam_buffer>0)then
xyzzyaaac36=xyzzyaaco1(iparam_buffer)
xyzzyaaad36=iparam_buffer-sum(xyzzyaacn1(1:xyzzyaaac36-1))
endif
xyzzyaaab36=0
if(xyzzyaacn1(1)>0)then
xyzzyaaaa36=xyzzyaaab36+1
xyzzyaaab36=xyzzyaaab36+xyzzyaacn1(1)
if(xyzzyaaac36==1.or.xyzzyaaac36==0)then
call put_mdet_params(params(xyzzyaaaa36:xyzzyaaab36),ignore(xyzzyaaaa3&
&6:xyzzyaaab36),xyzzyaaad36,prestore,xyzzyaaae36)
bad_params=bad_params.or.xyzzyaaae36
call wfdet_mdet_to_cmdet
endif
endif
if(xyzzyaacn1(2)>0)then
xyzzyaaaa36=xyzzyaaab36+1
xyzzyaaab36=xyzzyaaab36+xyzzyaacn1(2)
if(xyzzyaaac36==2.or.xyzzyaaac36==0)then
call put_wfdet_params(params(xyzzyaaaa36:xyzzyaaab36),ignore(xyzzyaaaa&
&36:xyzzyaaab36),xyzzyaaad36,prestore,xyzzyaaae36)
bad_params=bad_params.or.xyzzyaaae36
endif
endif
if(xyzzyaacn1(3)>0)then
xyzzyaaaa36=xyzzyaaab36+1
xyzzyaaab36=xyzzyaaab36+xyzzyaacn1(3)
if(xyzzyaaac36==3.or.xyzzyaaac36==0)then
call put_bf_params(params(xyzzyaaaa36:xyzzyaaab36),ignore(xyzzyaaaa36:&
&xyzzyaaab36),xyzzyaaad36,prestore,xyzzyaaae36)
bad_params=bad_params.or.xyzzyaaae36
endif
endif
end subroutine put_slater_params
subroutine invalidate_param1_slater(is,iparam)
implicit none
integer,intent(in) :: is,iparam
integer xyzzyaaaa37,xyzzyaaab37
xyzzyaaaa37=xyzzyaaco1(iparam)
xyzzyaaab37=iparam-sum(xyzzyaacn1(1:xyzzyaaaa37-1))
select case(xyzzyaaaa37)
case(1)
xyzzyaaaq1(is)=.false.
xyzzyaabj1(:,is)=.false.
xyzzyaabl1(:,is)=.false.
xyzzyaace1(is)=.false.
xyzzyaacg1(is)=.false.
if(.not.wfdet_detcoef_affect_orbs())then
continue
else
xyzzyaaay1(is)=.false.
xyzzyaaaw1(is)=0
xyzzyaabd1(is)=.false.
xyzzyaabb1(is)=0
xyzzyaabf1(is)=.false.
xyzzyaabh1(is)=.false.
xyzzyaaam1(is)=.false.
xyzzyaabn1(:,is)=.false.
xyzzyaabp1(:,is)=.false.
xyzzyaaas1(is)=.false.
xyzzyaabt1(is)=.false.
xyzzyaaba1(is)=.false.
xyzzyaaci1(is)=.false.
xyzzyaack1(is)=.false.
xyzzyaaav1(is)=.false.
endif
case(2,3)
xyzzyaaay1(is)=.false.
xyzzyaaaw1(is)=0
xyzzyaabd1(is)=.false.
xyzzyaabb1(is)=0
xyzzyaabf1(is)=.false.
xyzzyaabh1(is)=.false.
xyzzyaaam1(is)=.false.
xyzzyaaaq1(is)=.false.
xyzzyaabj1(:,is)=.false.
xyzzyaabl1(:,is)=.false.
xyzzyaabn1(:,is)=.false.
xyzzyaabp1(:,is)=.false.
xyzzyaaas1(is)=.false.
xyzzyaabt1(is)=.false.
xyzzyaaba1(is)=.false.
xyzzyaaci1(is)=.false.
xyzzyaack1(is)=.false.
xyzzyaaav1(is)=.false.
xyzzyaace1(is)=.false.
xyzzyaacg1(is)=.false.
if(xyzzyaaaa37==3)call invalidate_param1_bf(is,xyzzyaaab37)
end select
end subroutine invalidate_param1_slater
subroutine invalidate_params_slater(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaacm1)
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38
logical xyzzyaaad38(xyzzyaacl1)
xyzzyaaad38=.false.
xyzzyaaab38=0
do xyzzyaaac38=1,xyzzyaacl1
if(xyzzyaacn1(xyzzyaaac38)>0)then
xyzzyaaaa38=xyzzyaaab38+1
xyzzyaaab38=xyzzyaaab38+xyzzyaacn1(xyzzyaaac38)
xyzzyaaad38(xyzzyaaac38)=any(.not.ignore(xyzzyaaaa38:xyzzyaaab38))
endif
enddo
xyzzyaado1(:)=.false.
xyzzyaadp1(:)=.false.
if(xyzzyaaad38(3).or.xyzzyaaad38(2).or.wfdet_detcoef_affect_orbs())the&
&n
xyzzyaadq1(:)=.false.
xyzzyaadr1(:)=.false.
endif
end subroutine invalidate_params_slater
subroutine gen_config_slater(pt_config)
implicit none
type(config_wfn_slater),pointer :: pt_config
integer xyzzyaaaa39
allocate(pt_config,stat=xyzzyaaaa39)
call check_alloc(xyzzyaaaa39,'GEN_CONFIG_WFN_SLATER','container')
allocate(pt_config%pt_dbar(nemax,nemax,real1_complex2,nspin,ndet),pt_c&
&onfig%pt_logdet(nspin,ndet),stat=xyzzyaaaa39)
call check_alloc(xyzzyaaaa39,'GEN_CONFIG_WFN_SLATER','dbar')
pt_config%dbar_age=-1
if(orbbuf)then
allocate(pt_config%pt_dsmat(3,nemax,real1_complex2,ndet,nemax,nspin),s&
&tat=xyzzyaaaa39)
call check_alloc(xyzzyaaaa39,'GEN_CONFIG_WFN_SLATER','dsmat')
if(.not.use_backflow)then
allocate(pt_config%pt_lapsmat(nemax,real1_complex2,ndet,nemax,nspin),s&
&tat=xyzzyaaaa39)
call check_alloc(xyzzyaaaa39,'GEN_CONFIG_WFN_SLATER','lapsmat')
else
allocate(pt_config%pt_d2smat(6,nemax,real1_complex2,ndet,nemax,nspin),&
&stat=xyzzyaaaa39)
call check_alloc(xyzzyaaaa39,'GEN_CONFIG_WFN_SLATER','d2smat')
endif
else
nullify(pt_config%pt_dsmat,pt_config%pt_lapsmat,pt_config%pt_d2smat)
endif
pt_config%orb_dvalid=.false.
pt_config%orb_d2valid=.false.
end subroutine gen_config_slater
subroutine delete_config_slater(pt_config)
implicit none
type(config_wfn_slater),pointer :: pt_config
deallocate(pt_config%pt_dbar,pt_config%pt_logdet)
if(orbbuf)then
deallocate(pt_config%pt_dsmat)
if(.not.use_backflow)then
deallocate(pt_config%pt_lapsmat)
else
deallocate(pt_config%pt_d2smat)
endif
endif
deallocate(pt_config)
end subroutine delete_config_slater
subroutine copy_config_slater(pt_from,pt_to)
implicit none
type(config_wfn_slater),pointer :: pt_from,pt_to
call zcopy(size_det,pt_from%pt_logdet,1,pt_to%pt_logdet,1)
if(pt_from%dbar_age>=0)call dcopy(size_dbar,pt_from%pt_dbar,1,pt_to%pt&
&_dbar,1)
pt_to%dbar_age=pt_from%dbar_age
if(orbbuf)then
if(pt_from%orb_dvalid)call dcopy(size_dsmat,pt_from%pt_dsmat,1,pt_to%p&
&t_dsmat,1)
pt_to%orb_dvalid=pt_from%orb_dvalid
if(pt_from%orb_d2valid)then
if(.not.use_backflow)then
call dcopy(size_dbar,pt_from%pt_lapsmat,1,pt_to%pt_lapsmat,1)
else
call dcopy(size_d2smat,pt_from%pt_d2smat,1,pt_to%pt_d2smat,1)
endif
endif
pt_to%orb_d2valid=pt_from%orb_d2valid
endif
end subroutine copy_config_slater
subroutine config_to_pt_slater(pt_config,k)
use slaarnaaf, only : logdet_config
implicit none
integer,intent(in) :: k
type(config_wfn_slater),pointer :: pt_config
call zcopy(size_det,logdet_config(1,1,k),1,pt_config%pt_logdet,1)
pt_config%pt_dbar=0.d0
pt_config%dbar_age=-1
end subroutine config_to_pt_slater
subroutine pt_to_config_slater(pt_config)
use slaarnaaf, only : add_config
implicit none
type(config_wfn_slater),pointer :: pt_config
call add_config(modify=.true.,logdet=pt_config%pt_logdet)
end subroutine pt_to_config_slater
subroutine redist_allocations_slater(kmax)
implicit none
integer,intent(in) :: kmax
integer xyzzyaaaa44
allocate(xyzzyaacr1(nspin,ndet,kmax),stat=xyzzyaaaa44)
call check_alloc(xyzzyaaaa44,'REDIST_ALLOCATIONS_WFN_SLATER','1')
if(.not.small_transfer)then
allocate(xyzzyaacq1(nemax,nemax,real1_complex2,nspin,ndet,kmax),xyzzya&
&acp1(kmax),stat=xyzzyaaaa44)
call check_alloc(xyzzyaaaa44,'REDIST_ALLOCATIONS_WFN_SLATER','2')
endif
end subroutine redist_allocations_slater
subroutine redist_load_slater(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_slater),pointer :: pt_config
call zcopy(size_det,pt_config%pt_logdet,1,xyzzyaacr1(1,1,k),1)
if(.not.small_transfer)then
call dcopy(size_dbar,pt_config%pt_dbar,1,xyzzyaacq1(1,1,1,1,1,k),1)
xyzzyaacp1(k)=pt_config%dbar_age
endif
end subroutine redist_load_slater
subroutine redist_send_slater(jnode,kbase,k,nbt,reqbase,blocking)
use parallel
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa46
if(blocking)then
call mpi_ssend(xyzzyaacr1(1,1,kbase+1),k*size_det,mpi_double_complex,j&
&node,move_msg,mpi_comm_world,ierror)
call checkmpi(ierror,'ssend logdet in redist_send')
if(.not.small_transfer)then
call mpi_ssend(xyzzyaacq1(1,1,1,1,1,kbase+1),k*size_dbar,mpi_double_pr&
&ecision,jnode,move_msg,mpi_comm_world,ierror)
call checkmpi(ierror,'ssend dbar in redist_send')
call mpi_ssend(xyzzyaacp1(kbase+1),k,mpi_integer,jnode,move_msg,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'ssend dbar_age in redist_send')
endif
else
nbt=nbt+1
call mpi_isend(xyzzyaacr1(1,1,kbase+1),k*size_det,mpi_double_complex,j&
&node,nbt,mpi_comm_world,xyzzyaaaa46,ierror)
call checkmpi(ierror,'ssend logdet in redist_send')
nbreq(reqbase+nbt)=xyzzyaaaa46
if(.not.small_transfer)then
nbt=nbt+1
call mpi_isend(xyzzyaacq1(1,1,1,1,1,kbase+1),k*size_dbar,mpi_double_pr&
&ecision,jnode,nbt,mpi_comm_world,xyzzyaaaa46,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa46
call checkmpi(ierror,'ssend dbar in redist_send')
nbt=nbt+1
call mpi_isend(xyzzyaacp1(kbase+1),k,mpi_integer,jnode,nbt,mpi_comm_wo&
&rld,xyzzyaaaa46,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa46
call checkmpi(ierror,'ssend dbar_age in redist_send')
endif
endif
end subroutine redist_send_slater
subroutine redist_recv_slater(jnode,kbase,k,nbt,reqbase,blocking)
use parallel
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa47
if(blocking)then
call mpi_recv(xyzzyaacr1(1,1,kbase+1),k*size_det,mpi_double_complex,jn&
&ode,move_msg,mpi_comm_world,status,ierror)
call checkmpi(ierror,'recv logdet in redist_recv')
if(.not.small_transfer)then
call mpi_recv(xyzzyaacq1(1,1,1,1,1,kbase+1),k*size_dbar,mpi_double_pre&
&cision,jnode,move_msg,mpi_comm_world,status,ierror)
call checkmpi(ierror,'recv dbar in redist_recv')
call mpi_recv(xyzzyaacp1(kbase+1),k,mpi_integer,jnode,move_msg,mpi_com&
&m_world,status,ierror)
call checkmpi(ierror,'recv dbar_age in redist_recv')
endif
else
nbt=nbt+1
call mpi_irecv(xyzzyaacr1(1,1,kbase+1),k*size_det,mpi_double_complex,j&
&node,nbt,mpi_comm_world,xyzzyaaaa47,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa47
call checkmpi(ierror,'recv logdet in redist_recv')
if(.not.small_transfer)then
nbt=nbt+1
call mpi_irecv(xyzzyaacq1(1,1,1,1,1,kbase+1),k*size_dbar,mpi_double_pr&
&ecision,jnode,nbt,mpi_comm_world,xyzzyaaaa47,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa47
call checkmpi(ierror,'recv dbar in redist_recv')
nbt=nbt+1
call mpi_irecv(xyzzyaacp1(kbase+1),k,mpi_integer,jnode,nbt,mpi_comm_wo&
&rld,xyzzyaaaa47,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa47
call checkmpi(ierror,'recv dbar_age in redist_recv')
endif
endif
end subroutine redist_recv_slater
subroutine redist_save_slater(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_slater),pointer :: pt_config
call zcopy(size_det,xyzzyaacr1(1,1,k),1,pt_config%pt_logdet,1)
if(.not.small_transfer)then
call dcopy(size_dbar,xyzzyaacq1(1,1,1,1,1,k),1,pt_config%pt_dbar,1)
pt_config%dbar_age=xyzzyaacp1(k)
else
pt_config%dbar_age=-1
endif
end subroutine redist_save_slater
subroutine redist_deallocations_slater
implicit none
deallocate(xyzzyaacr1)
if(.not.small_transfer)deallocate(xyzzyaacq1,xyzzyaacp1)
end subroutine redist_deallocations_slater
subroutine load_from_pt_slater(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_slater),pointer :: pt_config
call zcopy(size_det,pt_config%pt_logdet,1,xyzzyaaaj1(1,1,is),1)
where(dble(xyzzyaaaj1(:,:,is))<real(xyzzyaaac1,dp))xyzzyaaah1(:,:,is)=&
&xyzzyaaaf1
where(dble(xyzzyaaaj1(:,:,is))<real(xyzzyaaad1,dp))xyzzyaaah1(:,:,is)=&
&xyzzyaaag1
where(xyzzyaaah1(:,:,is)/=xyzzyaaae1)xyzzyaaaj1(:,:,is)=czero
call xyzzyaafl1(xyzzyaaaj1(1,1,is),xyzzyaaah1(1,1,is),xyzzyaaak1(1,is)&
&,xyzzyaaal1(is),xyzzyaaai1(1,is))
xyzzyaaam1(is)=.true.
dbar_age(is)=pt_config%dbar_age
if(dbar_age(is)>=0)then
call dcopy(size_dbar,pt_config%pt_dbar,1,xyzzyaaaz1(1,1,1,1,1,is),1)
xyzzyaaba1(is)=.true.
else
dbar_age(is)=0
xyzzyaaba1(is)=.false.
endif
if(orbbuf)then
xyzzyaabd1(is)=pt_config%orb_dvalid
xyzzyaabb1(is)=0
if(xyzzyaabd1(is))call dcopy(size_dsmat,pt_config%pt_dsmat,1,xyzzyaabc&
&1(1,1,1,1,1,1,is),1)
if(.not.use_backflow)then
xyzzyaabf1(is)=pt_config%orb_d2valid
if(xyzzyaabf1(is))call dcopy(size_dbar,pt_config%pt_lapsmat,1,xyzzyaab&
&e1(1,1,1,1,1,is),1)
call xyzzyaafe1(is,.false.,.true.)
else
xyzzyaabh1(is)=pt_config%orb_d2valid
if(xyzzyaabh1(is))call dcopy(size_d2smat,pt_config%pt_d2smat,1,xyzzyaa&
&bg1(1,1,1,1,1,1,is),1)
call xyzzyaaff1(is,.false.,.true.,.true.)
endif
endif
end subroutine load_from_pt_slater
subroutine save_to_pt_slater(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_slater),pointer :: pt_config
call xyzzyaaew1(is,.false.,.false.)
call zcopy(size_det,xyzzyaaaj1(1,1,is),1,pt_config%pt_logdet,1)
if(xyzzyaaba1(is))then
pt_config%dbar_age=dbar_age(is)
call dcopy(size_dbar,xyzzyaaaz1(1,1,1,1,1,is),1,pt_config%pt_dbar,1)
else
pt_config%dbar_age=-1
endif
if(orbbuf)then
pt_config%orb_dvalid=xyzzyaabd1(is)
if(xyzzyaabd1(is))call dcopy(size_dsmat,xyzzyaabc1(1,1,1,1,1,1,is),1,p&
&t_config%pt_dsmat,1)
if(.not.use_backflow)then
pt_config%orb_d2valid=xyzzyaabf1(is)
if(xyzzyaabf1(is))call dcopy(size_dbar,xyzzyaabe1(1,1,1,1,1,is),1,pt_c&
&onfig%pt_lapsmat,1)
else
pt_config%orb_d2valid=xyzzyaabh1(is)
if(xyzzyaabh1(is))call dcopy(size_d2smat,xyzzyaabg1(1,1,1,1,1,1,is),1,&
&pt_config%pt_d2smat,1)
endif
endif
end subroutine save_to_pt_slater
subroutine add_config_slater_items(is)
use slaarnaaf, only : add_config
implicit none
integer,intent(in) :: is
integer xyzzyaaaa52,xyzzyaaab52
real(dp) lapdet
complex(dp) xyzzyaaac52(3),xyzzyaaad52,xyzzyaaae52(nspin,ndet)
logical xyzzyaaaf52,xyzzyaaag52
if(.not.(all(xyzzyaabj1(:,is)).and.all(xyzzyaabl1(:,is))))then
do xyzzyaaaa52=1,netot
call wfn_loggrad_slater(xyzzyaaaa52,is,0,.false.,.true.,xyzzyaaac52,xy&
&zzyaaaf52,xyzzyaaag52)
call wfn_loglap_slater(xyzzyaaaa52,is,0,.false.,xyzzyaaad52,xyzzyaaaf5&
&2,xyzzyaaag52)
enddo
endif
lapdet=0.d0
xyzzyaaaa52=0
do xyzzyaaab52=1,nspin
lapdet=lapdet+inv_pmass(xyzzyaaab52)*sum(xyzzyaabk1(xyzzyaaaa52+1:xyzz&
&yaaaa52+nele(xyzzyaaab52),1,is))
xyzzyaaaa52=xyzzyaaaa52+nele(xyzzyaaab52)
enddo
xyzzyaaae52(:,:)=xyzzyaaaj1(:,:,is)
where(xyzzyaaah1(:,:,is)==xyzzyaaaf1)xyzzyaaae52=xyzzyaaac1
where(xyzzyaaah1(:,:,is)==xyzzyaaag1)xyzzyaaae52=xyzzyaaad1
call add_config(modify=.true.,logdet=xyzzyaaae52,fidet=xyzzyaabi1(1,1,&
&1,is),lapdet=lapdet)
end subroutine add_config_slater_items
subroutine setup_storage_slater(nconfig,ignore)
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(xyzzyaacm1)
integer xyzzyaaaa53,xyzzyaaab53,xyzzyaaac53,xyzzyaaad53
logical xyzzyaaae53(xyzzyaacl1)
allocate(xyzzyaadn1(nspin,ndet,nconfig),xyzzyaadj1(3,netot,real1_compl&
&ex2,nconfig),xyzzyaadk1(netot,real1_complex2,nconfig),stat=xyzzyaaad5&
&3)
call check_alloc(xyzzyaaad53,'SETUP_STORAGE_SLATER','*det_store')
xyzzyaadn1=czero
xyzzyaadj1=0.d0
xyzzyaadk1=0.d0
allocate(xyzzyaads1(nconfig),xyzzyaado1(nconfig),xyzzyaadp1(nconfig),s&
&tat=xyzzyaaad53)
call check_alloc(xyzzyaaad53,'SETUP_STORAGE_SLATER','*det_svalid')
xyzzyaads1=.false.
xyzzyaado1=.false.
xyzzyaadp1=.false.
xyzzyaaae53=.false.
xyzzyaaac53=0
do xyzzyaaaa53=1,xyzzyaacl1
if(xyzzyaacn1(xyzzyaaaa53)>0)then
xyzzyaaab53=xyzzyaaac53+1
xyzzyaaac53=xyzzyaaac53+xyzzyaacn1(xyzzyaaaa53)
xyzzyaaae53(xyzzyaaaa53)=any(.not.ignore(xyzzyaaab53:xyzzyaaac53))
endif
enddo
if(xyzzyaaae53(1).and..not.use_backflow.and..not.small_buffers)then
allocate(xyzzyaadl1(3,ndet,nemax,real1_complex2,nspin,nconfig),xyzzyaa&
&dm1(ndet,nemax,real1_complex2,nspin,nconfig),stat=xyzzyaaad53)
call check_alloc(xyzzyaaad53,'SETUP_STORAGE_SLATER','fi_prod_det_store&
&')
xyzzyaadl1=0.d0
xyzzyaadm1=0.d0
endif
allocate(xyzzyaadq1(nconfig),xyzzyaadr1(nconfig),stat=xyzzyaaad53)
call check_alloc(xyzzyaaad53,'SETUP_STORAGE_SLATER','fi_prod_det_svali&
&d')
xyzzyaadq1=.false.
xyzzyaadr1=.false.
end subroutine setup_storage_slater
subroutine finish_storage_slater
implicit none
deallocate(xyzzyaadn1,xyzzyaadj1,xyzzyaadk1)
deallocate(xyzzyaads1,xyzzyaado1,xyzzyaadp1)
if(allocated(xyzzyaadl1))deallocate(xyzzyaadl1)
if(allocated(xyzzyaadm1))deallocate(xyzzyaadm1)
deallocate(xyzzyaadq1,xyzzyaadr1)
end subroutine finish_storage_slater
subroutine load_from_storage_slater(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(xyzzyaads1(icfg))then
xyzzyaaaj1(:,:,is)=xyzzyaadn1(:,:,icfg)
xyzzyaaah1(:,:,is)=xyzzyaaae1
where(dble(xyzzyaaaj1(:,:,is))<real(xyzzyaaac1,dp))xyzzyaaah1(:,:,is)=&
&xyzzyaaaf1
where(dble(xyzzyaaaj1(:,:,is))<real(xyzzyaaad1,dp))xyzzyaaah1(:,:,is)=&
&xyzzyaaag1
where(xyzzyaaah1(:,:,is)/=xyzzyaaae1)xyzzyaaaj1(:,:,is)=czero
call xyzzyaafl1(xyzzyaaaj1(1,1,is),xyzzyaaah1(1,1,is),xyzzyaaak1(1,is)&
&,xyzzyaaal1(is),xyzzyaaai1(1,is))
xyzzyaaam1(is)=.true.
endif
if(xyzzyaado1(icfg))then
call dcopy(size_fidet,xyzzyaadj1(1,1,1,icfg),1,xyzzyaabi1(1,1,1,is),1)
xyzzyaabj1(:,is)=.true.
endif
if(xyzzyaadp1(icfg))then
call dcopy(netot*real1_complex2,xyzzyaadk1(1,1,icfg),1,xyzzyaabk1(1,1,&
&is),1)
xyzzyaabl1(:,is)=.true.
endif
if(xyzzyaadq1(icfg))then
call dcopy(size_fi_prod_det,xyzzyaadl1(1,1,1,1,1,icfg),1,xyzzyaabm1(1,&
&1,1,1,1,is),1)
xyzzyaabn1(:,is)=.true.
endif
if(xyzzyaadr1(icfg))then
call dcopy(size_prod_lapdet,xyzzyaadm1(1,1,1,1,icfg),1,xyzzyaabo1(1,1,&
&1,1,is),1)
xyzzyaabp1(:,is)=.true.
endif
if(opt_orbitals)then
xyzzyaaam1=.false.
xyzzyaaaq1=.false.
endif
end subroutine load_from_storage_slater
subroutine save_to_storage_slater(is,icfg)
implicit none
integer,intent(in) :: is,icfg
complex(dp) xyzzyaaaa56(nspin,ndet)
if(xyzzyaaam1(is).and..not.xyzzyaads1(icfg))then
xyzzyaaaa56=xyzzyaaaj1(:,:,is)
where(xyzzyaaah1(:,:,is)==xyzzyaaaf1)xyzzyaaaa56=xyzzyaaac1
where(xyzzyaaah1(:,:,is)==xyzzyaaag1)xyzzyaaaa56=xyzzyaaad1
xyzzyaadn1(:,:,icfg)=xyzzyaaaa56
xyzzyaads1(icfg)=.true.
endif
if(all(xyzzyaabj1(:,is)).and..not.xyzzyaado1(icfg))then
call dcopy(size_fidet,xyzzyaabi1(1,1,1,is),1,xyzzyaadj1(1,1,1,icfg),1)
xyzzyaado1(icfg)=.true.
endif
if(all(xyzzyaabl1(:,is)).and..not.xyzzyaadp1(icfg))then
call dcopy(netot*real1_complex2,xyzzyaabk1(1,1,is),1,xyzzyaadk1(1,1,ic&
&fg),1)
xyzzyaadp1(icfg)=.true.
endif
if(allocated(xyzzyaadl1))then
if(all(xyzzyaabn1(:,is)).and..not.xyzzyaadq1(icfg))then
call dcopy(size_fi_prod_det,xyzzyaabm1(1,1,1,1,1,is),1,xyzzyaadl1(1,1,&
&1,1,1,icfg),1)
xyzzyaadq1(icfg)=.true.
endif
endif
if(allocated(xyzzyaadm1))then
if(all(xyzzyaabp1(:,is)).and..not.xyzzyaadr1(icfg))then
call dcopy(size_prod_lapdet,xyzzyaabo1(1,1,1,1,is),1,xyzzyaadm1(1,1,1,&
&1,icfg),1)
xyzzyaadr1(icfg)=.true.
endif
endif
end subroutine save_to_storage_slater
subroutine enumerate_plot_slater(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
integer xyzzyaaaa57,xyzzyaaab57,xyzzyaaac57,xyzzyaaad57
logical xyzzyaaae57,xyzzyaaaf57
xyzzyaaaf57=.not.(present(keyword).and.present(description))
if(xyzzyaaaf57)xyzzyaaep1=0
xyzzyaaab57=0
xyzzyaaeo1=0
do xyzzyaaaa57=1,xyzzyaaek1
xyzzyaaae57=.false.
select case(xyzzyaaaa57)
case(xyzzyaael1)
xyzzyaaae57=.true.
case(xyzzyaaem1)
xyzzyaaae57=.true.
case(xyzzyaaen1)
xyzzyaaae57=.true.
end select
if(xyzzyaaae57)then
xyzzyaaab57=xyzzyaaab57+1
xyzzyaaeo1(xyzzyaaab57)=xyzzyaaaa57
endif
enddo
xyzzyaaep1(0)=xyzzyaaab57
if(.not.xyzzyaaaf57)then
do xyzzyaaaa57=1,xyzzyaaab57
keyword(xyzzyaaaa57)=xyzzyaaer1(xyzzyaaeo1(xyzzyaaaa57))
description(xyzzyaaaa57)=xyzzyaaes1(xyzzyaaeo1(xyzzyaaaa57))
enddo
endif
if(xyzzyaaaf57)then
call enumerate_plot_wfdet(xyzzyaaep1(1))
if(use_backflow)call enumerate_plot_bf(xyzzyaaep1(2))
else
xyzzyaaaa57=xyzzyaaep1(0)+1
xyzzyaaab57=sum(xyzzyaaep1(0:1))
call enumerate_plot_wfdet(xyzzyaaep1(1),keyword(xyzzyaaaa57:xyzzyaaab5&
&7),description(xyzzyaaaa57:xyzzyaaab57))
xyzzyaaaa57=sum(xyzzyaaep1(0:1))+1
xyzzyaaab57=sum(xyzzyaaep1(0:2))
if(use_backflow)call enumerate_plot_bf(xyzzyaaep1(2),keyword(xyzzyaaaa&
&57:xyzzyaaab57),description(xyzzyaaaa57:xyzzyaaab57))
endif
n=sum(xyzzyaaep1(0:2))
if(.not.xyzzyaaaf57)then
allocate(xyzzyaaeq1(sum(xyzzyaaep1(0:2))),stat=xyzzyaaac57)
call check_alloc(xyzzyaaac57,'ENUMERATE_PLOT_SLATER','which_plot_sec')
do xyzzyaaad57=0,2
xyzzyaaeq1(sum(xyzzyaaep1(0:xyzzyaaad57-1))+1:sum(xyzzyaaep1(0:xyzzyaa&
&ad57)))=xyzzyaaad57
enddo
endif
end subroutine enumerate_plot_slater
subroutine query_plot_slater(iplot,ii,rank,is_complex,has_stderr,rot_t&
&ensor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
integer xyzzyaaaa58,iplot_rel,ie,idet
logical count_only
count_only=.not.present(function_name)
if(.not.allocated(xyzzyaaeq1))call errstop_master('QUERY_PLOT_SLATER',&
&'ENUMERATE_PLOT_SLATER not called yet. Bug in calling routine.')
if(iplot<0.or.iplot>sum(xyzzyaaep1(0:2)))call errstop_master('QUERY_PL&
&OT_SLATER','IPLOT out of range. Bug in calling routine.')
xyzzyaaaa58=xyzzyaaeq1(iplot)
iplot_rel=iplot-sum(xyzzyaaep1(0:xyzzyaaaa58-1))
select case(xyzzyaaaa58)
case(0)
select case(xyzzyaaeo1(iplot_rel))
case(xyzzyaael1)
rank=0
is_complex=complex_wf
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
if(ndet==1)then
do ie=1,nele(which_spin(ii))
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Orbital #'//trim(i2s(ie))
endif
enddo
else
do idet=detstart,detstop
do ie=1,nele(which_spin(ii))
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Orbital #'//trim(i2s(ie))//' in determinant&
& #'//trim(i2s(idet))
endif
enddo
enddo
endif
case(xyzzyaaem1)
rank=1
is_complex=complex_wf
has_stderr=.false.
rot_tensor=.true.
transl_pos=.false.
nfunctions=0
if(ndet==1)then
do ie=1,nele(which_spin(ii))
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Gradient of orbital #'//trim(i2s(ie))
endif
enddo
else
do idet=detstart,detstop
do ie=1,nele(which_spin(ii))
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Gradient of orbital #'//trim(i2s(ie))//' in&
& determinant #'//trim(i2s(idet))
endif
enddo
enddo
endif
case(xyzzyaaen1)
rank=0
is_complex=complex_wf
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=0
if(ndet==1)then
do ie=1,nele(which_spin(ii))
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Laplacian of orbital #'//trim(i2s(ie))
endif
enddo
else
do idet=detstart,detstop
do ie=1,nele(which_spin(ii))
nfunctions=nfunctions+1
if(.not.count_only)then
function_name(nfunctions)='Laplacian of orbital #'//trim(i2s(ie))//' i&
&n determinant #'//trim(i2s(idet))
endif
enddo
enddo
endif
end select
case(1)
call query_plot_wfdet(iplot_rel,ii,rank,is_complex,has_stderr,rot_tens&
&or,transl_pos,nfunctions,function_name)
case(2)
call query_plot_bf(iplot_rel,ii,rank,is_complex,has_stderr,rot_tensor,&
&transl_pos,nfunctions,function_name)
end select
end subroutine query_plot_slater
subroutine get_plot_slater(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
integer xyzzyaaaa59,xyzzyaaab59,xyzzyaaac59,xyzzyaaad59,xyzzyaaae59,xy&
&zzyaaaf59,xyzzyaaag59
if(.not.allocated(xyzzyaaeq1))call errstop_master('GET_PLOT_WFN','ENUM&
&ERATE_PLOT_WFN not called yet. Bug in calling routine.')
if(iplot<0.or.iplot>sum(xyzzyaaep1(0:2)))call errstop_master('GET_PLOT&
&_WFN','IPLOT out of range. Bug in calling routine.')
xyzzyaaaa59=xyzzyaaeq1(iplot)
xyzzyaaab59=iplot-sum(xyzzyaaep1(0:xyzzyaaaa59-1))
select case(xyzzyaaaa59)
case(0)
xyzzyaaac59=which_spin(ii)
select case(xyzzyaaeo1(xyzzyaaab59))
case(xyzzyaael1)
call xyzzyaaez1(ii,xyzzyaaac59,is1,.true.,.false.)
xyzzyaaaf59=0
do xyzzyaaad59=detstart,detstop
do xyzzyaaae59=1,nele(xyzzyaaac59)
xyzzyaaaf59=xyzzyaaaf59+1
f(xyzzyaaaf59)=xyzzyaabw1(xyzzyaaae59,1,xyzzyaaad59)
if(complex_wf)then
xyzzyaaaf59=xyzzyaaaf59+1
f(xyzzyaaaf59)=xyzzyaabw1(xyzzyaaae59,1,xyzzyaaad59)
endif
enddo
enddo
case(xyzzyaaem1)
call xyzzyaaez1(ii,xyzzyaaac59,is1,.false.,.true.)
xyzzyaaaf59=0
do xyzzyaaad59=detstart,detstop
do xyzzyaaae59=1,nele(xyzzyaaac59)
do xyzzyaaag59=1,dimensionality
xyzzyaaaf59=xyzzyaaaf59+1
f(xyzzyaaaf59)=xyzzyaaca1(xyzzyaaag59,xyzzyaaae59,1,xyzzyaaad59)
if(complex_wf)then
xyzzyaaaf59=xyzzyaaaf59+1
f(xyzzyaaaf59)=xyzzyaaca1(xyzzyaaag59,xyzzyaaae59,1,xyzzyaaad59)
endif
enddo
enddo
enddo
case(xyzzyaaen1)
call xyzzyaaez1(ii,xyzzyaaac59,is1,.false.,.true.)
xyzzyaaaf59=0
do xyzzyaaad59=detstart,detstop
do xyzzyaaae59=1,nele(xyzzyaaac59)
xyzzyaaaf59=xyzzyaaaf59+1
f(xyzzyaaaf59)=xyzzyaacb1(xyzzyaaae59,1,xyzzyaaad59)
if(complex_wf)then
xyzzyaaaf59=xyzzyaaaf59+1
f(xyzzyaaaf59)=xyzzyaacb1(xyzzyaaae59,1,xyzzyaaad59)
endif
enddo
enddo
end select
case(1)
call get_plot_wfdet(xyzzyaaab59,ii,is0,is1,f)
case(2)
call get_plot_bf(xyzzyaaab59,ii,is0,is1,f)
end select
end subroutine get_plot_slater
subroutine finish_plot_slater
implicit none
xyzzyaaeo1=0
xyzzyaaep1(0:2)=0
deallocate(xyzzyaaeq1)
call finish_plot_wfdet
if(use_backflow)call finish_plot_bf
end subroutine finish_plot_slater
subroutine xyzzyaafl1(logdet,fpeinfo_det,log_pdet,lognorm_wfn,fpeinfo_&
&pdet)
implicit none
integer,intent(in) :: fpeinfo_det(nspin,ndet)
integer,intent(inout) :: fpeinfo_pdet(ndet)
complex(dp),intent(in) :: logdet(nspin,ndet)
complex(dp),intent(out) :: log_pdet(ndet),lognorm_wfn
integer xyzzyaaaa61,xyzzyaaab61,xyzzyaaac61
call timer('COMPUTE_LOG_PDET',.true.)
log_pdet(:)=czero
fpeinfo_pdet(:)=xyzzyaaae1
lognorm_wfn=czero
do xyzzyaaaa61=detstart,detstop
do xyzzyaaab61=1,nspin
if(missing_det(xyzzyaaab61,xyzzyaaaa61))cycle
if(fpeinfo_det(xyzzyaaab61,xyzzyaaaa61)/=xyzzyaaae1)then
fpeinfo_pdet(xyzzyaaaa61)=max(fpeinfo_pdet(xyzzyaaaa61),fpeinfo_det(xy&
&zzyaaab61,xyzzyaaaa61))
else
log_pdet(xyzzyaaaa61)=log_pdet(xyzzyaaaa61)+logdet(xyzzyaaab61,xyzzyaa&
&aa61)
endif
enddo
enddo
fpeinfo_pdet(detstop+1:ndet)=xyzzyaaag1
if(any(fpeinfo_pdet(detstart:detstop)==xyzzyaaae1))then
xyzzyaaac61=maxloc(real(log_pdet(detstart:detstop),dp),1,fpeinfo_pdet(&
&detstart:detstop)==xyzzyaaae1)
lognorm_wfn=log_pdet(detstart-1+xyzzyaaac61)
where(fpeinfo_pdet(detstart:detstop)==xyzzyaaae1)log_pdet(detstart:det&
&stop)=log_pdet(detstart:detstop)-lognorm_wfn
endif
call timer('COMPUTE_LOG_PDET',.false.)
end subroutine xyzzyaafl1
subroutine xyzzyaafm1(log_pdet,fpeinfo_pdet,c_pdet,logwfn_renorm,fpein&
&fo_wfn)
implicit none
integer,intent(in) :: fpeinfo_pdet(ndet)
integer,intent(out) :: fpeinfo_wfn
complex(dp),intent(in) :: log_pdet(ndet)
complex(dp),intent(out) :: c_pdet(ndet),logwfn_renorm
complex(dp) xyzzyaaaa62
call timer('COMPUTE_LOGWFN',.true.)
c_pdet=czero
logwfn_renorm=czero
fpeinfo_wfn=xyzzyaaae1
if(any(fpeinfo_pdet==xyzzyaaae1))then
where(fpeinfo_pdet==xyzzyaaae1)c_pdet=detcoef*exp(log_pdet)
xyzzyaaaa62=sum(c_pdet,fpeinfo_pdet==xyzzyaaae1)
if(xyzzyaaaa62==czero)then
fpeinfo_wfn=xyzzyaaaf1
else
logwfn_renorm=log(xyzzyaaaa62)
endif
else
fpeinfo_wfn=xyzzyaaaf1
endif
call timer('COMPUTE_LOGWFN',.false.)
end subroutine xyzzyaafm1
subroutine xyzzyaafn1(log_pdet,fpeinfo_pdet,c_pdet,logwfn_renorm,fpein&
&fo_wfn)
implicit none
integer,intent(in) :: fpeinfo_pdet(ndet)
integer,intent(out) :: fpeinfo_wfn
complex(dp),intent(in) :: log_pdet(ndet)
complex(dp),intent(out) :: c_pdet(ndet),logwfn_renorm
integer xyzzyaaaa63
complex(dp) xyzzyaaab63,xyzzyaaac63,xyzzyaaad63
call timer('COMPUTE_LOGWFN_EFF',.true.)
c_pdet=czero
xyzzyaaab63=czero
logwfn_renorm=czero
fpeinfo_wfn=xyzzyaaae1
if(any(fpeinfo_pdet==xyzzyaaae1))then
where(fpeinfo_pdet(detstart:detstop)==xyzzyaaae1)c_pdet(detstart:detst&
&op)=detcoef(detstart:detstop)*exp(log_pdet(detstart:detstop))
xyzzyaaac63=detcoef(detstart)**2
xyzzyaaad63=c_pdet(detstart)
do xyzzyaaaa63=detstart+1,detstop
if(detcoef_label(xyzzyaaaa63-1)==detcoef_label(xyzzyaaaa63))then
xyzzyaaac63=xyzzyaaac63+detcoef(xyzzyaaaa63)**2
xyzzyaaad63=xyzzyaaad63+c_pdet(xyzzyaaaa63)
else
xyzzyaaab63=xyzzyaaab63+(dble(xyzzyaaad63)**2+aimag(xyzzyaaad63)**2)/x&
&yzzyaaac63
xyzzyaaac63=detcoef(xyzzyaaaa63)**2
xyzzyaaad63=c_pdet(xyzzyaaaa63)
endif
enddo
xyzzyaaab63=xyzzyaaab63+(dble(xyzzyaaad63)**2+aimag(xyzzyaaad63)**2)/x&
&yzzyaaac63
if(xyzzyaaab63==czero)then
fpeinfo_wfn=xyzzyaaaf1
else
logwfn_renorm=0.5d0*log(xyzzyaaab63)
endif
else
fpeinfo_wfn=xyzzyaaaf1
endif
call timer('COMPUTE_LOGWFN_EFF',.false.)
end subroutine xyzzyaafn1
subroutine xyzzyaafo1(ie,ispin,idet,dbar,rpsi,q,fpeinfo,orb_m,orb_rmap&
&)
implicit none
integer,intent(in) :: ie,ispin,idet,orb_m,orb_rmap(nemax)
integer,intent(inout) :: fpeinfo
real(dp),intent(in) :: dbar(nemax,nemax,real1_complex2),rpsi(nemax,rea&
&l1_complex2)
complex(dp),intent(out) :: q
integer xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64
real(dp) xyzzyaaad64,xyzzyaaae64
call timer('CALC_Q',.true.)
fpeinfo=xyzzyaaae1
xyzzyaaab64=upd_spin(ispin,idet)
xyzzyaaac64=nele(xyzzyaaab64)
if(complex_wf)then
if(sparse)then
if(update_by_column(ispin,idet))then
xyzzyaaad64=ddot_s(orb_m,orb_rmap,dbar(1,ie,1),1,rpsi,1)-ddot_s(orb_m,&
&orb_rmap,dbar(1,ie,2),1,rpsi(1,2),1)
xyzzyaaae64=ddot_s(orb_m,orb_rmap,dbar(1,ie,2),1,rpsi,1)+ddot_s(orb_m,&
&orb_rmap,dbar(1,ie,1),1,rpsi(1,2),1)
else
xyzzyaaaa64=ie+nuc_nele(xyzzyaaab64)
xyzzyaaad64=ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa64,1,1),nemax,rpsi,1)-&
&ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa64,1,2),nemax,rpsi(1,2),1)
xyzzyaaae64=ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa64,1,2),nemax,rpsi,1)+&
&ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa64,1,1),nemax,rpsi(1,2),1)
endif
if(xyzzyaaad64*xyzzyaaad64+xyzzyaaae64*xyzzyaaae64<1.d-200)then
xyzzyaaad64=0.d0
xyzzyaaae64=0.d0
endif
else
if(update_by_column(ispin,idet))then
xyzzyaaad64=ddot(xyzzyaaac64,dbar(1,ie,1),1,rpsi,1)-ddot(xyzzyaaac64,d&
&bar(1,ie,2),1,rpsi(1,2),1)
xyzzyaaae64=ddot(xyzzyaaac64,dbar(1,ie,2),1,rpsi,1)+ddot(xyzzyaaac64,d&
&bar(1,ie,1),1,rpsi(1,2),1)
else
xyzzyaaaa64=ie+nuc_nele(xyzzyaaab64)
xyzzyaaad64=ddot(xyzzyaaac64,dbar(xyzzyaaaa64,1,1),nemax,rpsi,1)-ddot(&
&xyzzyaaac64,dbar(xyzzyaaaa64,1,2),nemax,rpsi(1,2),1)
xyzzyaaae64=ddot(xyzzyaaac64,dbar(xyzzyaaaa64,1,2),nemax,rpsi,1)+ddot(&
&xyzzyaaac64,dbar(xyzzyaaaa64,1,1),nemax,rpsi(1,2),1)
endif
endif
q=cmplx(xyzzyaaad64,xyzzyaaae64,dp)
if(xyzzyaaad64==0.d0.and.xyzzyaaae64==0.d0)then
if(all(rpsi(1:xyzzyaaac64,1:2)==0.d0))then
fpeinfo=xyzzyaaag1
else
fpeinfo=xyzzyaaaf1
endif
endif
else
if(sparse)then
if(update_by_column(ispin,idet))then
xyzzyaaad64=ddot_s(orb_m,orb_rmap,dbar(1,ie,1),1,rpsi,1)
else
xyzzyaaaa64=ie+nuc_nele(xyzzyaaab64)
xyzzyaaad64=ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa64,1,1),nemax,rpsi,1)
endif
if(xyzzyaaad64*xyzzyaaad64<1.d-200)xyzzyaaad64=0.d0
else
if(update_by_column(ispin,idet))then
xyzzyaaad64=ddot(xyzzyaaac64,dbar(1,ie,1),1,rpsi,1)
else
xyzzyaaaa64=ie+nuc_nele(xyzzyaaab64)
xyzzyaaad64=ddot(xyzzyaaac64,dbar(xyzzyaaaa64,1,1),nemax,rpsi,1)
endif
endif
q=cmplx(xyzzyaaad64,0.d0,dp)
if(xyzzyaaad64==0.d0)then
if(all(rpsi(1:xyzzyaaac64,1)==0.d0))then
fpeinfo=xyzzyaaag1
else
fpeinfo=xyzzyaaaf1
endif
endif
endif
call timer('CALC_Q',.false.)
end subroutine xyzzyaafo1
subroutine xyzzyaafp1(m,rmap,ispin,idet,dbar,schunk,sbar,sbar_piv,q,fp&
&einfo)
implicit none
integer,intent(in) :: m,rmap(m),ispin,idet
integer,intent(inout) :: sbar_piv(nemax),fpeinfo
real(dp),intent(in) :: dbar(nemax,nemax,real1_complex2),schunk(nemax,r&
&eal1_complex2,ndet,nemax,nspin)
real(dp),intent(inout) :: sbar(nemax,nemax,real1_complex2)
complex(dp),intent(out) :: q
integer xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65,xyzzyaaad65,xyzzyaaae65
logical xyzzyaaaf65
xyzzyaaae65=nele(ispin)
q=c_one
xyzzyaaaf65=.false.
fpeinfo=xyzzyaaae1
if(m==0)return
call timer('CALC_Q_BF',.true.)
if(complex_wf)then
if(update_by_column(ispin,idet))then
do xyzzyaaad65=1,m
xyzzyaaab65=rmap(xyzzyaaad65)
do xyzzyaaac65=1,m
xyzzyaaaa65=rmap(xyzzyaaac65)
sbar(xyzzyaaac65,xyzzyaaad65,1)=ddot(xyzzyaaae65,dbar(1,xyzzyaaaa65,1)&
&,1,schunk(1,1,idet,xyzzyaaab65,ispin),1)-ddot(xyzzyaaae65,dbar(1,xyzz&
&yaaaa65,2),1,schunk(1,2,idet,xyzzyaaab65,ispin),1)
sbar(xyzzyaaac65,xyzzyaaad65,2)=ddot(xyzzyaaae65,dbar(1,xyzzyaaaa65,1)&
&,1,schunk(1,2,idet,xyzzyaaab65,ispin),1)+ddot(xyzzyaaae65,dbar(1,xyzz&
&yaaaa65,2),1,schunk(1,1,idet,xyzzyaaab65,ispin),1)
enddo
enddo
else
do xyzzyaaac65=1,m
xyzzyaaaa65=rmap(xyzzyaaac65)
do xyzzyaaad65=1,m
xyzzyaaab65=rmap(xyzzyaaad65)
sbar(xyzzyaaad65,xyzzyaaac65,1)=ddot(xyzzyaaae65,dbar(xyzzyaaab65,1,1)&
&,nemax,schunk(1,1,idet,xyzzyaaaa65,ispin),1)-ddot(xyzzyaaae65,dbar(xy&
&zzyaaab65,1,2),nemax,schunk(1,2,idet,xyzzyaaaa65,ispin),1)
sbar(xyzzyaaad65,xyzzyaaac65,2)=ddot(xyzzyaaae65,dbar(xyzzyaaab65,1,1)&
&,nemax,schunk(1,2,idet,xyzzyaaaa65,ispin),1)+ddot(xyzzyaaae65,dbar(xy&
&zzyaaab65,1,2),nemax,schunk(1,1,idet,xyzzyaaaa65,ispin),1)
enddo
enddo
endif
call lu_decom_cmplx_dz(sbar,xyzzyaadi1,sbar_piv,nemax,m,xyzzyaaaf65)
fpeinfo=xyzzyaaae1
if(xyzzyaaaf65)then
q=czero
fpeinfo=xyzzyaaaf1
do xyzzyaaac65=1,m
if(all(schunk(1:xyzzyaaae65,1:2,idet,rmap(xyzzyaaac65),ispin)==0.d0))t&
&hen
fpeinfo=xyzzyaaag1
exit
endif
enddo
else
q=exp(lu_logdet_cmplx_dz(sbar,sbar_piv,nemax,m))
if(abs(q)<sqrt(tiny(1.d0)))then
q=czero
fpeinfo=xyzzyaaaf1
endif
endif
else
if(update_by_column(ispin,idet))then
do xyzzyaaad65=1,m
xyzzyaaab65=rmap(xyzzyaaad65)
do xyzzyaaac65=1,m
xyzzyaaaa65=rmap(xyzzyaaac65)
sbar(xyzzyaaac65,xyzzyaaad65,1)=ddot(xyzzyaaae65,dbar(1,xyzzyaaaa65,1)&
&,1,schunk(1,1,idet,xyzzyaaab65,ispin),1)
enddo
enddo
else
do xyzzyaaac65=1,m
xyzzyaaaa65=rmap(xyzzyaaac65)
do xyzzyaaad65=1,m
xyzzyaaab65=rmap(xyzzyaaad65)
sbar(xyzzyaaad65,xyzzyaaac65,1)=ddot(xyzzyaaae65,dbar(xyzzyaaab65,1,1)&
&,nemax,schunk(1,1,idet,xyzzyaaaa65,ispin),1)
enddo
enddo
endif
call lu_decom(sbar,sbar_piv,nemax,m,xyzzyaaaf65)
fpeinfo=xyzzyaaae1
if(xyzzyaaaf65)then
q=czero
fpeinfo=xyzzyaaaf1
do xyzzyaaac65=1,m
if(all(schunk(1:xyzzyaaae65,1:real1_complex2,idet,rmap(xyzzyaaac65),is&
&pin)==0.d0))then
fpeinfo=xyzzyaaag1
exit
endif
enddo
else
q=exp(lu_logdet(sbar,sbar_piv,nemax,m))
if(abs(q)<sqrt(tiny(1.d0)))then
q=czero
fpeinfo=xyzzyaaaf1
endif
endif
endif
call timer('CALC_Q_BF',.false.)
end subroutine xyzzyaafp1
subroutine xyzzyaafq1(ie,ispin,idet,rpsi,dbar,q,orb_m,orb_rmap)
implicit none
integer,intent(in) :: ie,ispin,idet,orb_m,orb_rmap(nemax)
real(dp),intent(in) :: rpsi(nemax,real1_complex2)
real(dp),intent(inout) :: dbar(nemax,nemax,real1_complex2)
complex(dp),intent(in) :: q
integer xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66,xyzzyaaad66
real(dp) xyzzyaaae66,xyzzyaaaf66,xyzzyaaag66,xyzzyaaah66,xyzzyaaai66
call timer('UPDATE_DBAR',.true.)
xyzzyaaab66=upd_spin(ispin,idet)
xyzzyaaad66=nele(xyzzyaaab66)
if(complex_wf)then
xyzzyaaag66=1.d0/(dble(q)**2+aimag(q)**2)
xyzzyaaae66=dble(q)*xyzzyaaag66
xyzzyaaaf66=-aimag(q)*xyzzyaaag66
if(update_by_column(ispin,idet))then
!$omp parallel do default(none) shared(ie,xyzzyaaae66,xyzzyaaaf66,xyzz&
!$omp &yaaad66,orb_m,orb_rmap,dbar,xyzzyaaab66,rpsi,sparse) private(xy&
!$omp &zzyaaaa66,xyzzyaaag66,xyzzyaaai66,xyzzyaaah66)
do xyzzyaaaa66=1,xyzzyaaad66
if(xyzzyaaaa66==ie)cycle
if(sparse)then
xyzzyaaai66=ddot_s(orb_m,orb_rmap,dbar(1,xyzzyaaaa66,1),1,rpsi,1)-ddot&
&_s(orb_m,orb_rmap,dbar(1,xyzzyaaaa66,2),1,rpsi(1,2),1)
xyzzyaaah66=ddot_s(orb_m,orb_rmap,dbar(1,xyzzyaaaa66,2),1,rpsi,1)+ddot&
&_s(orb_m,orb_rmap,dbar(1,xyzzyaaaa66,1),1,rpsi(1,2),1)
else
xyzzyaaai66=ddot(xyzzyaaad66,dbar(1,xyzzyaaaa66,1),1,rpsi,1)-ddot(xyzz&
&yaaad66,dbar(1,xyzzyaaaa66,2),1,rpsi(1,2),1)
xyzzyaaah66=ddot(xyzzyaaad66,dbar(1,xyzzyaaaa66,2),1,rpsi,1)+ddot(xyzz&
&yaaad66,dbar(1,xyzzyaaaa66,1),1,rpsi(1,2),1)
endif
xyzzyaaag66=-xyzzyaaae66*xyzzyaaai66+xyzzyaaaf66*xyzzyaaah66
xyzzyaaah66=-xyzzyaaae66*xyzzyaaah66-xyzzyaaaf66*xyzzyaaai66
call daxpy(xyzzyaaad66,xyzzyaaag66,dbar(1,ie,1),1,dbar(1,xyzzyaaaa66,1&
&),1)
call daxpy(xyzzyaaad66,-xyzzyaaah66,dbar(1,ie,2),1,dbar(1,xyzzyaaaa66,&
&1),1)
call daxpy(xyzzyaaad66,xyzzyaaag66,dbar(1,ie,2),1,dbar(1,xyzzyaaaa66,2&
&),1)
call daxpy(xyzzyaaad66,xyzzyaaah66,dbar(1,ie,1),1,dbar(1,xyzzyaaaa66,2&
&),1)
enddo
!$omp end parallel do
call dcopy(xyzzyaaad66,dbar(1,ie,1),1,xyzzyaacv1,1)
call dscal(xyzzyaaad66,xyzzyaaae66,dbar(1,ie,1),1)
call daxpy(xyzzyaaad66,-xyzzyaaaf66,dbar(1,ie,2),1,dbar(1,ie,1),1)
call dscal(xyzzyaaad66,xyzzyaaae66,dbar(1,ie,2),1)
call daxpy(xyzzyaaad66,xyzzyaaaf66,xyzzyaacv1,1,dbar(1,ie,2),1)
else
xyzzyaaac66=ie+nuc_nele(xyzzyaaab66)
do xyzzyaaaa66=1,xyzzyaaad66
if(xyzzyaaaa66==xyzzyaaac66)cycle
if(sparse)then
xyzzyaaai66=ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa66,1,1),nemax,rpsi,1)-&
&ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa66,1,2),nemax,rpsi(1,2),1)
xyzzyaaah66=ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa66,1,2),nemax,rpsi,1)+&
&ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa66,1,1),nemax,rpsi(1,2),1)
else
xyzzyaaai66=ddot(xyzzyaaad66,dbar(xyzzyaaaa66,1,1),nemax,rpsi,1)-ddot(&
&xyzzyaaad66,dbar(xyzzyaaaa66,1,2),nemax,rpsi(1,2),1)
xyzzyaaah66=ddot(xyzzyaaad66,dbar(xyzzyaaaa66,1,2),nemax,rpsi,1)+ddot(&
&xyzzyaaad66,dbar(xyzzyaaaa66,1,1),nemax,rpsi(1,2),1)
endif
xyzzyaaag66=-xyzzyaaae66*xyzzyaaai66+xyzzyaaaf66*xyzzyaaah66
xyzzyaaah66=-xyzzyaaae66*xyzzyaaah66-xyzzyaaaf66*xyzzyaaai66
call daxpy(xyzzyaaad66,xyzzyaaag66,dbar(xyzzyaaac66,1,1),nemax,dbar(xy&
&zzyaaaa66,1,1),nemax)
call daxpy(xyzzyaaad66,-xyzzyaaah66,dbar(xyzzyaaac66,1,2),nemax,dbar(x&
&yzzyaaaa66,1,1),nemax)
call daxpy(xyzzyaaad66,xyzzyaaag66,dbar(xyzzyaaac66,1,2),nemax,dbar(xy&
&zzyaaaa66,1,2),nemax)
call daxpy(xyzzyaaad66,xyzzyaaah66,dbar(xyzzyaaac66,1,1),nemax,dbar(xy&
&zzyaaaa66,1,2),nemax)
enddo
call dcopy(xyzzyaaad66,dbar(xyzzyaaac66,1,1),nemax,xyzzyaacv1,1)
call dscal(xyzzyaaad66,xyzzyaaae66,dbar(xyzzyaaac66,1,1),nemax)
call daxpy(xyzzyaaad66,-xyzzyaaaf66,dbar(xyzzyaaac66,1,2),nemax,dbar(x&
&yzzyaaac66,1,1),nemax)
call dscal(xyzzyaaad66,xyzzyaaae66,dbar(xyzzyaaac66,1,2),nemax)
call daxpy(xyzzyaaad66,xyzzyaaaf66,xyzzyaacv1,1,dbar(xyzzyaaac66,1,2),&
&nemax)
endif
else
xyzzyaaae66=1.d0/dble(q)
if(update_by_column(ispin,idet))then
!$omp parallel do default(none) shared(ie,xyzzyaaae66,xyzzyaaad66,orb_&
!$omp &m,orb_rmap,dbar,xyzzyaaab66,rpsi,sparse) private(xyzzyaaaa66,xy&
!$omp &zzyaaag66)
do xyzzyaaaa66=1,xyzzyaaad66
if(xyzzyaaaa66==ie)cycle
if(sparse)then
xyzzyaaag66=-xyzzyaaae66*ddot_s(orb_m,orb_rmap,dbar(1,xyzzyaaaa66,1),1&
&,rpsi,1)
else
xyzzyaaag66=-xyzzyaaae66*ddot(xyzzyaaad66,dbar(1,xyzzyaaaa66,1),1,rpsi&
&,1)
endif
call daxpy(xyzzyaaad66,xyzzyaaag66,dbar(1,ie,1),1,dbar(1,xyzzyaaaa66,1&
&),1)
enddo
!$omp end parallel do
call dscal(xyzzyaaad66,xyzzyaaae66,dbar(1,ie,1),1)
else
xyzzyaaac66=ie+nuc_nele(xyzzyaaab66)
do xyzzyaaaa66=1,xyzzyaaad66
if(xyzzyaaaa66==xyzzyaaac66)cycle
if(sparse)then
xyzzyaaag66=-xyzzyaaae66*ddot_s(orb_m,orb_rmap,dbar(xyzzyaaaa66,1,1),n&
&emax,rpsi,1)
else
xyzzyaaag66=-xyzzyaaae66*ddot(xyzzyaaad66,dbar(xyzzyaaaa66,1,1),nemax,&
&rpsi,1)
endif
call daxpy(xyzzyaaad66,xyzzyaaag66,dbar(xyzzyaaac66,1,1),nemax,dbar(xy&
&zzyaaaa66,1,1),nemax)
enddo
call dscal(xyzzyaaad66,xyzzyaaae66,dbar(xyzzyaaac66,1,1),nemax)
endif
endif
call timer('UPDATE_DBAR',.false.)
end subroutine xyzzyaafq1
subroutine xyzzyaafr1(m,rmap,ispin,idet,dbar,schunk,sbar,sbar_piv)
implicit none
integer,intent(in) :: m,rmap(nemax),ispin,idet,sbar_piv(nemax)
real(dp),intent(inout) :: sbar(nemax,nemax,real1_complex2),dbar(nemax,&
&nemax,real1_complex2),schunk(nemax,real1_complex2,ndet,nemax,nspin)
integer xyzzyaaaa67,xyzzyaaab67,xyzzyaaac67,xyzzyaaad67,xyzzyaaae67,xy&
&zzyaaaf67,xyzzyaaag67
real(dp) xyzzyaaah67,xyzzyaaai67
xyzzyaaag67=nele(ispin)
if(m==0)return
call timer('UPDATE_DBAR_BF',.true.)
if(complex_wf)then
if(update_by_column(ispin,idet))then
do xyzzyaaab67=1,m
xyzzyaaae67=rmap(xyzzyaaab67)
do xyzzyaaaa67=m+1,xyzzyaaag67
xyzzyaaad67=rmap(xyzzyaaaa67)
sbar(xyzzyaaaa67,xyzzyaaab67,1)=ddot(xyzzyaaag67,dbar(1,xyzzyaaad67,1)&
&,1,schunk(1,1,idet,xyzzyaaae67,ispin),1)-ddot(xyzzyaaag67,dbar(1,xyzz&
&yaaad67,2),1,schunk(1,2,idet,xyzzyaaae67,ispin),1)
sbar(xyzzyaaaa67,xyzzyaaab67,2)=ddot(xyzzyaaag67,dbar(1,xyzzyaaad67,1)&
&,1,schunk(1,2,idet,xyzzyaaae67,ispin),1)+ddot(xyzzyaaag67,dbar(1,xyzz&
&yaaad67,2),1,schunk(1,1,idet,xyzzyaaae67,ispin),1)
enddo
enddo
else
do xyzzyaaaa67=1,m
xyzzyaaad67=rmap(xyzzyaaaa67)
do xyzzyaaab67=m+1,xyzzyaaag67
xyzzyaaae67=rmap(xyzzyaaab67)
sbar(xyzzyaaab67,xyzzyaaaa67,1)=ddot(xyzzyaaag67,dbar(xyzzyaaae67,1,1)&
&,nemax,schunk(1,1,idet,xyzzyaaad67,ispin),1)-ddot(xyzzyaaag67,dbar(xy&
&zzyaaae67,1,2),nemax,schunk(1,2,idet,xyzzyaaad67,ispin),1)
sbar(xyzzyaaab67,xyzzyaaaa67,2)=ddot(xyzzyaaag67,dbar(xyzzyaaae67,1,1)&
&,nemax,schunk(1,2,idet,xyzzyaaad67,ispin),1)+ddot(xyzzyaaag67,dbar(xy&
&zzyaaae67,1,2),nemax,schunk(1,1,idet,xyzzyaaad67,ispin),1)
enddo
enddo
endif
do xyzzyaaaa67=1,m
xyzzyaacv1=0.d0
xyzzyaacv1(xyzzyaaaa67,1)=1.d0
call lu_solve_once_cmplx_dz(sbar,xyzzyaadi1,sbar_piv,xyzzyaacv1,xyzzya&
&adh1,nemax,m)
call dcopy(m,xyzzyaacv1,1,xyzzyaacw1(1,xyzzyaaaa67,1),1)
call dcopy(m,xyzzyaacv1(1,2),1,xyzzyaacw1(1,xyzzyaaaa67,2),1)
enddo
do xyzzyaaaa67=1,m
call dcopy(m,xyzzyaacw1(1,xyzzyaaaa67,1),1,sbar(1,xyzzyaaaa67,1),1)
call dcopy(m,xyzzyaacw1(1,xyzzyaaaa67,2),1,sbar(1,xyzzyaaaa67,2),1)
enddo
do xyzzyaaaa67=m+1,xyzzyaaag67
do xyzzyaaab67=1,m
xyzzyaacv1(xyzzyaaab67,1)=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,sbar(1,xy&
&zzyaaab67,1),1)-ddot(m,sbar(xyzzyaaaa67,1,2),nemax,sbar(1,xyzzyaaab67&
&,2),1)
xyzzyaacv1(xyzzyaaab67,2)=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,sbar(1,xy&
&zzyaaab67,2),1)+ddot(m,sbar(xyzzyaaaa67,1,2),nemax,sbar(1,xyzzyaaab67&
&,1),1)
enddo
call dcopy(m,xyzzyaacv1,1,sbar(xyzzyaaaa67,1,1),nemax)
call dcopy(m,xyzzyaacv1(1,2),1,sbar(xyzzyaaaa67,1,2),nemax)
enddo
if(update_by_column(ispin,idet))then
do xyzzyaaae67=1,xyzzyaaag67
do xyzzyaaac67=1,m
xyzzyaaaf67=rmap(xyzzyaaac67)
xyzzyaacv1(xyzzyaaac67,1)=dbar(xyzzyaaae67,xyzzyaaaf67,1)
xyzzyaacv1(xyzzyaaac67,2)=dbar(xyzzyaaae67,xyzzyaaaf67,2)
enddo
do xyzzyaaaa67=1,m
xyzzyaaad67=rmap(xyzzyaaaa67)
xyzzyaaah67=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,xyzzyaacv1,1)-ddot(m,sb&
&ar(xyzzyaaaa67,1,2),nemax,xyzzyaacv1(1,2),1)
xyzzyaaai67=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,xyzzyaacv1(1,2),1)+ddot&
&(m,sbar(xyzzyaaaa67,1,2),nemax,xyzzyaacv1,1)
dbar(xyzzyaaae67,xyzzyaaad67,1)=xyzzyaaah67
dbar(xyzzyaaae67,xyzzyaaad67,2)=xyzzyaaai67
enddo
do xyzzyaaaa67=m+1,xyzzyaaag67
xyzzyaaad67=rmap(xyzzyaaaa67)
xyzzyaaah67=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,xyzzyaacv1,1)-ddot(m,sb&
&ar(xyzzyaaaa67,1,2),nemax,xyzzyaacv1(1,2),1)
xyzzyaaai67=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,xyzzyaacv1(1,2),1)+ddot&
&(m,sbar(xyzzyaaaa67,1,2),nemax,xyzzyaacv1,1)
dbar(xyzzyaaae67,xyzzyaaad67,1)=dbar(xyzzyaaae67,xyzzyaaad67,1)-xyzzya&
&aah67
dbar(xyzzyaaae67,xyzzyaaad67,2)=dbar(xyzzyaaae67,xyzzyaaad67,2)-xyzzya&
&aai67
enddo
enddo
else
do xyzzyaaad67=1,xyzzyaaag67
do xyzzyaaac67=1,m
xyzzyaaaf67=rmap(xyzzyaaac67)
xyzzyaacv1(xyzzyaaac67,1)=dbar(xyzzyaaaf67,xyzzyaaad67,1)
xyzzyaacv1(xyzzyaaac67,2)=dbar(xyzzyaaaf67,xyzzyaaad67,2)
enddo
do xyzzyaaab67=1,m
xyzzyaaae67=rmap(xyzzyaaab67)
xyzzyaaah67=ddot(m,sbar(xyzzyaaab67,1,1),nemax,xyzzyaacv1,1)-ddot(m,sb&
&ar(xyzzyaaab67,1,2),nemax,xyzzyaacv1(1,2),1)
xyzzyaaai67=ddot(m,sbar(xyzzyaaab67,1,1),nemax,xyzzyaacv1(1,2),1)-ddot&
&(m,sbar(xyzzyaaab67,1,2),nemax,xyzzyaacv1,1)
dbar(xyzzyaaae67,xyzzyaaad67,1)=xyzzyaaah67
dbar(xyzzyaaae67,xyzzyaaad67,2)=xyzzyaaai67
enddo
do xyzzyaaab67=m+1,xyzzyaaag67
xyzzyaaae67=rmap(xyzzyaaab67)
xyzzyaaah67=ddot(m,sbar(xyzzyaaab67,1,1),nemax,xyzzyaacv1,1)-ddot(m,sb&
&ar(xyzzyaaab67,1,2),nemax,xyzzyaacv1(1,2),1)
xyzzyaaai67=ddot(m,sbar(xyzzyaaab67,1,1),nemax,xyzzyaacv1(1,2),1)-ddot&
&(m,sbar(xyzzyaaab67,1,2),nemax,xyzzyaacv1,1)
dbar(xyzzyaaae67,xyzzyaaad67,1)=dbar(xyzzyaaae67,xyzzyaaad67,1)-xyzzya&
&aah67
dbar(xyzzyaaae67,xyzzyaaad67,2)=dbar(xyzzyaaae67,xyzzyaaad67,2)-xyzzya&
&aai67
enddo
enddo
endif
else
if(update_by_column(ispin,idet))then
do xyzzyaaab67=1,m
xyzzyaaae67=rmap(xyzzyaaab67)
do xyzzyaaaa67=m+1,xyzzyaaag67
xyzzyaaad67=rmap(xyzzyaaaa67)
sbar(xyzzyaaaa67,xyzzyaaab67,1)=ddot(xyzzyaaag67,dbar(1,xyzzyaaad67,1)&
&,1,schunk(1,1,idet,xyzzyaaae67,ispin),1)
enddo
enddo
else
do xyzzyaaaa67=1,m
xyzzyaaad67=rmap(xyzzyaaaa67)
do xyzzyaaab67=m+1,xyzzyaaag67
xyzzyaaae67=rmap(xyzzyaaab67)
sbar(xyzzyaaab67,xyzzyaaaa67,1)=ddot(xyzzyaaag67,dbar(xyzzyaaae67,1,1)&
&,nemax,schunk(1,1,idet,xyzzyaaad67,ispin),1)
enddo
enddo
endif
xyzzyaacw1(1:m,1:m,1)=0.d0
do xyzzyaaaa67=1,m
xyzzyaacw1(xyzzyaaaa67,xyzzyaaaa67,1)=1.d0
enddo
call lu_solve_n(sbar,sbar_piv,xyzzyaacw1,nemax,m,nemax)
do xyzzyaaaa67=1,m
call dcopy(m,xyzzyaacw1(1,xyzzyaaaa67,1),1,sbar(1,xyzzyaaaa67,1),1)
enddo
do xyzzyaaaa67=m+1,xyzzyaaag67
do xyzzyaaab67=1,m
xyzzyaacv1(xyzzyaaab67,1)=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,sbar(1,xy&
&zzyaaab67,1),1)
enddo
call dcopy(m,xyzzyaacv1,1,sbar(xyzzyaaaa67,1,1),nemax)
enddo
if(update_by_column(ispin,idet))then
do xyzzyaaae67=1,xyzzyaaag67
do xyzzyaaac67=1,m
xyzzyaaaf67=rmap(xyzzyaaac67)
xyzzyaacv1(xyzzyaaac67,1)=dbar(xyzzyaaae67,xyzzyaaaf67,1)
enddo
do xyzzyaaaa67=1,m
xyzzyaaad67=rmap(xyzzyaaaa67)
xyzzyaaah67=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,xyzzyaacv1,1)
dbar(xyzzyaaae67,xyzzyaaad67,1)=xyzzyaaah67
enddo
do xyzzyaaaa67=m+1,xyzzyaaag67
xyzzyaaad67=rmap(xyzzyaaaa67)
xyzzyaaah67=ddot(m,sbar(xyzzyaaaa67,1,1),nemax,xyzzyaacv1,1)
dbar(xyzzyaaae67,xyzzyaaad67,1)=dbar(xyzzyaaae67,xyzzyaaad67,1)-xyzzya&
&aah67
enddo
enddo
else
do xyzzyaaad67=1,xyzzyaaag67
do xyzzyaaac67=1,m
xyzzyaaaf67=rmap(xyzzyaaac67)
xyzzyaacv1(xyzzyaaac67,1)=dbar(xyzzyaaaf67,xyzzyaaad67,1)
enddo
do xyzzyaaab67=1,m
xyzzyaaae67=rmap(xyzzyaaab67)
xyzzyaaah67=ddot(m,sbar(xyzzyaaab67,1,1),nemax,xyzzyaacv1,1)
dbar(xyzzyaaae67,xyzzyaaad67,1)=xyzzyaaah67
enddo
do xyzzyaaab67=m+1,xyzzyaaag67
xyzzyaaae67=rmap(xyzzyaaab67)
xyzzyaaah67=ddot(m,sbar(xyzzyaaab67,1,1),nemax,xyzzyaacv1,1)
dbar(xyzzyaaae67,xyzzyaaad67,1)=dbar(xyzzyaaae67,xyzzyaaad67,1)-xyzzya&
&aah67
enddo
enddo
endif
endif
call timer('UPDATE_DBAR_BF',.false.)
end subroutine xyzzyaafr1
subroutine xyzzyaafs1(rele,sele,val,fd,sd,eevecs,smat,dsmat,lapsmat,d2&
&smat,logdet,fpeinfo_det)
implicit none
integer,intent(in) :: sele(netot)
integer,intent(inout),optional :: fpeinfo_det(nspin,ndet)
real(dp),intent(in) :: rele(3,netot)
real(dp),intent(in),target :: eevecs(4,netot,netot)
real(dp),intent(inout),optional,target :: smat(nemax,real1_complex2,nd&
&et,nemax,nspin),dsmat(3,nemax,real1_complex2,ndet,nemax,nspin),lapsma&
&t(nemax,real1_complex2,ndet,nemax,nspin),d2smat(6,nemax,real1_complex&
&2,ndet,nemax,nspin)
complex(dp),intent(inout),optional :: logdet(nspin,ndet)
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa68,xyzzyaaab68,xyzzyaaac68,xyzzyaaad68,xyzzyaaae68
real(dp),pointer :: xyzzyaaaf68(:,:,:),xyzzyaaag68(:,:,:,:),xyzzyaaah6&
&8(:,:,:),xyzzyaaai68(:,:),xyzzyaaaj68(:,:,:,:)
logical xyzzyaaak68,xyzzyaaal68,xyzzyaaam68,xyzzyaaan68
call timer('EVAL_SMAT',.true.,collapse=.true.)
xyzzyaaak68=fd.or.sd
xyzzyaaal68=present(d2smat).and.sd
xyzzyaaam68=present(lapsmat).and.sd
if(.not.val)xyzzyaaaf68=>xyzzyaacx1
if(.not.fd)xyzzyaaag68=>xyzzyaacy1
if(.not.xyzzyaaam68)xyzzyaaah68=>xyzzyaacz1
if(.not.xyzzyaaal68)xyzzyaaaj68=>xyzzyaada1
if(.not.pairing_wf)xyzzyaaai68=>eevecs1
do xyzzyaaab68=1,nspin
xyzzyaaac68=nele(xyzzyaaab68)
xyzzyaaan68=.not.pairing_wf.or.all(update_by_column(xyzzyaaab68,:))
if(.not.xyzzyaaan68)then
if(val)xyzzyaaaf68=>xyzzyaacx1
if(fd)xyzzyaaag68=>xyzzyaacy1
if(xyzzyaaam68)xyzzyaaah68=>xyzzyaacz1
if(xyzzyaaal68)xyzzyaaaj68=>xyzzyaada1
endif
do xyzzyaaad68=1,xyzzyaaac68
xyzzyaaae68=which_ii(xyzzyaaad68,xyzzyaaab68)
if(xyzzyaaan68)then
if(val)xyzzyaaaf68=>smat(:,:,:,xyzzyaaad68,xyzzyaaab68)
if(fd)xyzzyaaag68=>dsmat(:,:,:,:,xyzzyaaad68,xyzzyaaab68)
if(xyzzyaaam68)xyzzyaaah68=>lapsmat(:,:,:,xyzzyaaad68,xyzzyaaab68)
if(xyzzyaaal68)xyzzyaaaj68=>d2smat(:,:,:,:,xyzzyaaad68,xyzzyaaab68)
endif
if(pairing_wf)xyzzyaaai68=>eevecs(:,:,xyzzyaaae68)
if(.not.xyzzyaaal68)then
call wfdet(rele(1,xyzzyaaae68),sele(xyzzyaaae68),xyzzyaaab68,wfdet_nor&
&b,wfdet_orbmask(1,xyzzyaaab68),val,xyzzyaaak68,xyzzyaadc1,xyzzyaadd1,&
&xyzzyaade1,xyzzyaaai68)
if(val)call copy_orb_to_det(1,xyzzyaaab68,wfdet_orbmap,xyzzyaadc1,xyzz&
&yaaaf68)
if(xyzzyaaak68)then
call copy_orb_to_det(3,xyzzyaaab68,wfdet_orbmap,xyzzyaadd1,xyzzyaaag68&
&)
call copy_orb_to_det(1,xyzzyaaab68,wfdet_orbmap,xyzzyaade1,xyzzyaaah68&
&)
endif
else
call wfdet(rele(1,xyzzyaaae68),sele(xyzzyaaae68),xyzzyaaab68,wfdet_nor&
&b,wfdet_orbmask(1,xyzzyaaab68),val,xyzzyaaak68,xyzzyaadc1,xyzzyaadd1,&
&xyzzyaade1,xyzzyaaai68,xyzzyaadf1)
if(val)call copy_orb_to_det(1,xyzzyaaab68,wfdet_orbmap,xyzzyaadc1,xyzz&
&yaaaf68)
if(xyzzyaaak68)then
call copy_orb_to_det(3,xyzzyaaab68,wfdet_orbmap,xyzzyaadd1,xyzzyaaag68&
&)
call copy_orb_to_det(1,xyzzyaaab68,wfdet_orbmap,xyzzyaade1,xyzzyaaah68&
&)
call copy_orb_to_det(6,xyzzyaaab68,wfdet_orbmap,xyzzyaadf1,xyzzyaaaj68&
&)
endif
endif
if(.not.xyzzyaaan68)then
if(val)call xyzzyaafu1(xyzzyaaae68,xyzzyaacx1,smat)
if(fd)call xyzzyaafv1(xyzzyaaae68,xyzzyaacy1,dsmat)
if(xyzzyaaam68)call xyzzyaafw1(xyzzyaaae68,xyzzyaacz1,lapsmat)
if(xyzzyaaal68)call xyzzyaafx1(xyzzyaaae68,xyzzyaada1,d2smat)
endif
enddo
enddo
if(present(logdet).and.val)then
logdet=czero
fpeinfo_det=xyzzyaaae1
do xyzzyaaaa68=detstart,detstop
do xyzzyaaab68=1,nspin
if(missing_det(xyzzyaaab68,xyzzyaaaa68))cycle
call xyzzyaaga1(smat,xyzzyaaab68,xyzzyaaaa68,logdet(xyzzyaaab68,xyzzya&
&aaa68),fpeinfo_det(xyzzyaaab68,xyzzyaaaa68))
enddo
enddo
endif
call timer('EVAL_SMAT',.false.)
end subroutine xyzzyaafs1
subroutine xyzzyaaft1(rele,sele,m,rmap,val,fd,eevecs,smat,dsmat,logdet&
&,fpeinfo_det,untwist)
implicit none
integer,intent(in) :: sele(netot),m(nspin),rmap(nemax,nspin)
integer,intent(inout),optional :: fpeinfo_det(nspin,ndet)
real(dp),intent(in) :: rele(3,netot)
real(dp),intent(in),target :: eevecs(4,netot,netot)
real(dp),intent(inout),target,optional :: smat(nemax,real1_complex2,nd&
&et,nemax,nspin),dsmat(3,nemax,real1_complex2,ndet,nemax,nspin)
complex(dp),intent(inout),optional :: logdet(nspin,ndet)
logical,intent(in) :: val,fd
logical,intent(in),optional :: untwist
integer xyzzyaaaa69,xyzzyaaab69,xyzzyaaac69,xyzzyaaad69,xyzzyaaae69
real(dp),pointer :: xyzzyaaaf69(:,:,:),xyzzyaaag69(:,:,:,:),xyzzyaaah6&
&9(:,:)
logical xyzzyaaai69,xyzzyaaaj69
call timer('FILL_SMAT',.true.,collapse=.true.)
xyzzyaaai69=.not.pairing_wf
if(present(untwist))xyzzyaaai69=untwist
if(.not.val)xyzzyaaaf69=>xyzzyaacx1
if(.not.fd)xyzzyaaag69=>xyzzyaacy1
if(.not.pairing_wf)xyzzyaaah69=>eevecs1
do xyzzyaaab69=1,nspin
xyzzyaaaj69=xyzzyaaai69.or.all(update_by_column(xyzzyaaab69,:))
if(.not.xyzzyaaaj69)then
if(val)xyzzyaaaf69=>xyzzyaacx1
if(fd)xyzzyaaag69=>xyzzyaacy1
endif
do xyzzyaaae69=1,m(xyzzyaaab69)
xyzzyaaac69=rmap(xyzzyaaae69,xyzzyaaab69)
xyzzyaaad69=which_ii(xyzzyaaac69,xyzzyaaab69)
if(xyzzyaaaj69)then
if(val)xyzzyaaaf69=>smat(:,:,:,xyzzyaaac69,xyzzyaaab69)
if(fd)xyzzyaaag69=>dsmat(:,:,:,:,xyzzyaaac69,xyzzyaaab69)
endif
if(pairing_wf)xyzzyaaah69=>eevecs(:,:,xyzzyaaad69)
call wfdet(rele(1,xyzzyaaad69),sele(xyzzyaaad69),xyzzyaaab69,wfdet_nor&
&b,wfdet_orbmask(1,xyzzyaaab69),val,fd,xyzzyaadc1,xyzzyaadd1,xyzzyaade&
&1,xyzzyaaah69)
if(val)call copy_orb_to_det(1,xyzzyaaab69,wfdet_orbmap,xyzzyaadc1,xyzz&
&yaaaf69)
if(fd)call copy_orb_to_det(3,xyzzyaaab69,wfdet_orbmap,xyzzyaadd1,xyzzy&
&aaag69)
if(.not.xyzzyaaaj69)then
if(val)call xyzzyaafu1(xyzzyaaad69,xyzzyaacx1,smat)
if(fd)call xyzzyaafv1(xyzzyaaad69,xyzzyaacy1,dsmat)
endif
enddo
enddo
if(present(logdet).and.val)then
logdet=czero
fpeinfo_det=xyzzyaaae1
do xyzzyaaaa69=detstart,detstop
do xyzzyaaab69=1,nspin
if(missing_det(xyzzyaaab69,xyzzyaaaa69))cycle
call xyzzyaaga1(smat,xyzzyaaab69,xyzzyaaaa69,logdet(xyzzyaaab69,xyzzya&
&aaa69),fpeinfo_det(xyzzyaaab69,xyzzyaaaa69))
enddo
enddo
endif
call timer('FILL_SMAT',.false.)
end subroutine xyzzyaaft1
subroutine xyzzyaafu1(ii,smat1,smat,back)
implicit none
integer,intent(in) :: ii
real(dp),intent(inout) :: smat1(nemax,real1_complex2,ndet),smat(nemax,&
&real1_complex2,ndet,nemax,nspin)
logical,intent(in),optional :: back
integer xyzzyaaaa70,xyzzyaaab70,xyzzyaaac70,xyzzyaaad70
logical xyzzyaaae70
call timer('TRANSFER_SMAT',.true.)
xyzzyaaae70=.false.
if(present(back))xyzzyaaae70=back
xyzzyaaaa70=which_ie(ii)
xyzzyaaab70=which_spin(ii)
if(.not.pairing_wf)then
if(.not.xyzzyaaae70)then
call dcopy(size_rpsi,smat1,1,smat(1,1,1,xyzzyaaaa70,xyzzyaaab70),1)
else
call dcopy(size_rpsi,smat(1,1,1,xyzzyaaaa70,xyzzyaaab70),1,smat1,1)
endif
else
if(.not.xyzzyaaae70)then
do xyzzyaaac70=detstart,detstop
if(update_by_column(xyzzyaaab70,xyzzyaaac70))then
call dcopy(nele(xyzzyaaab70),smat1(1,1,xyzzyaaac70),1,smat(1,1,xyzzyaa&
&ac70,xyzzyaaaa70,xyzzyaaab70),1)
if(complex_wf)call dcopy(nele(xyzzyaaab70),smat1(1,2,xyzzyaaac70),1,sm&
&at(1,2,xyzzyaaac70,xyzzyaaaa70,xyzzyaaab70),1)
else
xyzzyaaad70=upd_spin(xyzzyaaab70,xyzzyaaac70)
call dcopy(nele(xyzzyaaab70),smat1(1,1,xyzzyaaac70),1,smat(xyzzyaaaa70&
&,1,xyzzyaaac70,1,xyzzyaaad70),size_rpsi)
if(complex_wf)call dcopy(nele(xyzzyaaab70),smat1(1,2,xyzzyaaac70),1,sm&
&at(xyzzyaaaa70,2,xyzzyaaac70,1,xyzzyaaad70),size_rpsi)
endif
enddo
else
do xyzzyaaac70=detstart,detstop
if(update_by_column(xyzzyaaab70,xyzzyaaac70))then
call dcopy(nele(xyzzyaaab70),smat(1,1,xyzzyaaac70,xyzzyaaaa70,xyzzyaaa&
&b70),1,smat1(1,1,xyzzyaaac70),1)
if(complex_wf)call dcopy(nele(xyzzyaaab70),smat(1,2,xyzzyaaac70,xyzzya&
&aaa70,xyzzyaaab70),1,smat1(1,2,xyzzyaaac70),1)
else
xyzzyaaad70=upd_spin(xyzzyaaab70,xyzzyaaac70)
call dcopy(nele(xyzzyaaab70),smat(xyzzyaaaa70,1,xyzzyaaac70,1,xyzzyaaa&
&d70),size_rpsi,smat1(1,1,xyzzyaaac70),1)
if(complex_wf)call dcopy(nele(xyzzyaaab70),smat(xyzzyaaaa70,2,xyzzyaaa&
&c70,1,xyzzyaaad70),size_rpsi,smat1(1,2,xyzzyaaac70),1)
endif
enddo
endif
endif
call timer('TRANSFER_SMAT',.false.)
end subroutine xyzzyaafu1
subroutine xyzzyaafv1(ii,dsmat1,dsmat,back)
implicit none
integer,intent(in) :: ii
real(dp),intent(inout) :: dsmat1(3,nemax,real1_complex2,ndet),dsmat(3,&
&nemax,real1_complex2,ndet,nemax,nspin)
logical,intent(in),optional :: back
integer xyzzyaaaa71,xyzzyaaab71,xyzzyaaac71,xyzzyaaad71
logical xyzzyaaae71
call timer('TRANSFER_DSMAT',.true.)
xyzzyaaae71=.false.
if(present(back))xyzzyaaae71=back
xyzzyaaaa71=which_ie(ii)
xyzzyaaab71=which_spin(ii)
if(.not.pairing_wf)then
if(.not.xyzzyaaae71)then
call dcopy(size_grad,dsmat1,1,dsmat(1,1,1,1,xyzzyaaaa71,xyzzyaaab71),1&
&)
else
call dcopy(size_grad,dsmat(1,1,1,1,xyzzyaaaa71,xyzzyaaab71),1,dsmat1,1&
&)
endif
else
if(.not.xyzzyaaae71)then
do xyzzyaaac71=detstart,detstop
if(update_by_column(xyzzyaaab71,xyzzyaaac71))then
call dcopy(three_nele(xyzzyaaab71),dsmat1(1,1,1,xyzzyaaac71),1,dsmat(1&
&,1,1,xyzzyaaac71,xyzzyaaaa71,xyzzyaaab71),1)
if(complex_wf)call dcopy(three_nele(xyzzyaaab71),dsmat1(1,1,2,xyzzyaaa&
&c71),1,dsmat(1,1,2,xyzzyaaac71,xyzzyaaaa71,xyzzyaaab71),1)
else
xyzzyaaad71=upd_spin(xyzzyaaab71,xyzzyaaac71)
call dcopy3(nele(xyzzyaaab71),three_nele(xyzzyaaab71),dsmat1(1,1,1,xyz&
&zyaaac71),1,dsmat(1,xyzzyaaaa71,1,xyzzyaaac71,1,xyzzyaaad71),size_gra&
&d-2)
if(complex_wf)call dcopy3(nele(xyzzyaaab71),three_nele(xyzzyaaab71),ds&
&mat1(1,1,2,xyzzyaaac71),1,dsmat(1,xyzzyaaaa71,2,xyzzyaaac71,1,xyzzyaa&
&ad71),size_grad-2)
endif
enddo
else
do xyzzyaaac71=detstart,detstop
if(update_by_column(xyzzyaaab71,xyzzyaaac71))then
call dcopy(three_nele(xyzzyaaab71),dsmat(1,1,1,xyzzyaaac71,xyzzyaaaa71&
&,xyzzyaaab71),1,dsmat1(1,1,1,xyzzyaaac71),1)
if(complex_wf)call dcopy(three_nele(xyzzyaaab71),dsmat(1,1,2,xyzzyaaac&
&71,xyzzyaaaa71,xyzzyaaab71),1,dsmat1(1,1,2,xyzzyaaac71),1)
else
xyzzyaaad71=upd_spin(xyzzyaaab71,xyzzyaaac71)
call dcopy3(nele(xyzzyaaab71),three_nele(xyzzyaaab71),dsmat(1,xyzzyaaa&
&a71,1,xyzzyaaac71,1,xyzzyaaad71),size_grad-2,dsmat1(1,1,1,xyzzyaaac71&
&),1)
if(complex_wf)call dcopy3(nele(xyzzyaaab71),three_nele(xyzzyaaab71),ds&
&mat(1,xyzzyaaaa71,2,xyzzyaaac71,1,xyzzyaaad71),size_grad-2,dsmat1(1,1&
&,2,xyzzyaaac71),1)
endif
enddo
endif
endif
call timer('TRANSFER_DSMAT',.false.)
end subroutine xyzzyaafv1
subroutine xyzzyaafw1(ii,lapsmat1,lapsmat,back)
implicit none
integer,intent(in) :: ii
real(dp),intent(inout) :: lapsmat1(nemax,real1_complex2,ndet),lapsmat(&
&nemax,real1_complex2,ndet,nemax,nspin)
logical,intent(in),optional :: back
integer xyzzyaaaa72,xyzzyaaab72,xyzzyaaac72,xyzzyaaad72
logical xyzzyaaae72
call timer('TRANSFER_LAPSMAT',.true.)
xyzzyaaae72=.false.
if(present(back))xyzzyaaae72=back
xyzzyaaaa72=which_ie(ii)
xyzzyaaab72=which_spin(ii)
if(.not.pairing_wf)then
if(.not.xyzzyaaae72)then
call dcopy(size_rpsi,lapsmat1,1,lapsmat(1,1,1,xyzzyaaaa72,xyzzyaaab72)&
&,1)
else
call dcopy(size_rpsi,lapsmat(1,1,1,xyzzyaaaa72,xyzzyaaab72),1,lapsmat1&
&,1)
endif
else
if(.not.xyzzyaaae72)then
do xyzzyaaac72=detstart,detstop
if(update_by_column(xyzzyaaab72,xyzzyaaac72))then
call dcopy(nele(xyzzyaaab72),lapsmat1(1,1,xyzzyaaac72),1,lapsmat(1,1,x&
&yzzyaaac72,xyzzyaaaa72,xyzzyaaab72),1)
if(complex_wf)call dcopy(nele(xyzzyaaab72),lapsmat1(1,2,xyzzyaaac72),1&
&,lapsmat(1,2,xyzzyaaac72,xyzzyaaaa72,xyzzyaaab72),1)
else
xyzzyaaad72=upd_spin(xyzzyaaab72,xyzzyaaac72)
call dcopy(nele(xyzzyaaab72),lapsmat1(1,1,xyzzyaaac72),1,lapsmat(xyzzy&
&aaaa72,1,xyzzyaaac72,1,xyzzyaaad72),size_rpsi)
if(complex_wf)call dcopy(nele(xyzzyaaab72),lapsmat1(1,2,xyzzyaaac72),1&
&,lapsmat(xyzzyaaaa72,2,xyzzyaaac72,1,xyzzyaaad72),size_rpsi)
endif
enddo
else
do xyzzyaaac72=detstart,detstop
if(update_by_column(xyzzyaaab72,xyzzyaaac72))then
call dcopy(nele(xyzzyaaab72),lapsmat(1,1,xyzzyaaac72,xyzzyaaaa72,xyzzy&
&aaab72),1,lapsmat1(1,1,xyzzyaaac72),1)
if(complex_wf)call dcopy(nele(xyzzyaaab72),lapsmat(1,2,xyzzyaaac72,xyz&
&zyaaaa72,xyzzyaaab72),1,lapsmat1(1,2,xyzzyaaac72),1)
else
xyzzyaaad72=upd_spin(xyzzyaaab72,xyzzyaaac72)
call dcopy(nele(xyzzyaaab72),lapsmat(xyzzyaaaa72,1,xyzzyaaac72,1,xyzzy&
&aaad72),size_rpsi,lapsmat1(1,1,xyzzyaaac72),1)
if(complex_wf)call dcopy(nele(xyzzyaaab72),lapsmat(xyzzyaaaa72,2,xyzzy&
&aaac72,1,xyzzyaaad72),size_rpsi,lapsmat1(1,2,xyzzyaaac72),1)
endif
enddo
endif
endif
call timer('TRANSFER_LAPSMAT',.false.)
end subroutine xyzzyaafw1
subroutine xyzzyaafx1(ii,d2smat1,d2smat,back)
implicit none
integer,intent(in) :: ii
real(dp),intent(inout) :: d2smat1(6,nemax,real1_complex2,ndet),d2smat(&
&6,nemax,real1_complex2,ndet,nemax,nspin)
logical,intent(in),optional :: back
integer xyzzyaaaa73,xyzzyaaab73,xyzzyaaac73,xyzzyaaad73
logical xyzzyaaae73
call timer('TRANSFER_D2SMAT',.true.)
xyzzyaaae73=.false.
if(present(back))xyzzyaaae73=back
xyzzyaaaa73=which_ie(ii)
xyzzyaaab73=which_spin(ii)
if(.not.pairing_wf)then
if(.not.xyzzyaaae73)then
call dcopy(6*size_rpsi,d2smat1,1,d2smat(1,1,1,1,xyzzyaaaa73,xyzzyaaab7&
&3),1)
else
call dcopy(6*size_rpsi,d2smat(1,1,1,1,xyzzyaaaa73,xyzzyaaab73),1,d2sma&
&t1,1)
endif
else
if(.not.xyzzyaaae73)then
do xyzzyaaac73=detstart,detstop
if(update_by_column(xyzzyaaab73,xyzzyaaac73))then
call dcopy(six_nele(xyzzyaaab73),d2smat1(1,1,1,xyzzyaaac73),1,d2smat(1&
&,1,1,xyzzyaaac73,xyzzyaaaa73,xyzzyaaab73),1)
if(complex_wf)call dcopy(six_nele(xyzzyaaab73),d2smat1(1,1,2,xyzzyaaac&
&73),1,d2smat(1,1,2,xyzzyaaac73,xyzzyaaaa73,xyzzyaaab73),1)
else
xyzzyaaad73=upd_spin(xyzzyaaab73,xyzzyaaac73)
call dcopy6(nele(xyzzyaaab73),six_nele(xyzzyaaab73),d2smat1(1,1,1,xyzz&
&yaaac73),1,d2smat(1,xyzzyaaaa73,1,xyzzyaaac73,1,xyzzyaaad73),6*size_r&
&psi-5)
if(complex_wf)call dcopy6(nele(xyzzyaaab73),six_nele(xyzzyaaab73),d2sm&
&at1(1,1,2,xyzzyaaac73),1,d2smat(1,xyzzyaaaa73,2,xyzzyaaac73,1,xyzzyaa&
&ad73),6*size_rpsi-5)
endif
enddo
else
do xyzzyaaac73=detstart,detstop
if(update_by_column(xyzzyaaab73,xyzzyaaac73))then
call dcopy(six_nele(xyzzyaaab73),d2smat(1,1,1,xyzzyaaac73,xyzzyaaaa73,&
&xyzzyaaab73),1,d2smat1(1,1,1,xyzzyaaac73),1)
if(complex_wf)call dcopy(six_nele(xyzzyaaab73),d2smat(1,1,2,xyzzyaaac7&
&3,xyzzyaaaa73,xyzzyaaab73),1,d2smat1(1,1,2,xyzzyaaac73),1)
else
xyzzyaaad73=upd_spin(xyzzyaaab73,xyzzyaaac73)
call dcopy6(nele(xyzzyaaab73),six_nele(xyzzyaaab73),d2smat(1,xyzzyaaaa&
&73,1,xyzzyaaac73,1,xyzzyaaad73),6*size_rpsi-5,d2smat1(1,1,1,xyzzyaaac&
&73),1)
if(complex_wf)call dcopy6(nele(xyzzyaaab73),six_nele(xyzzyaaab73),d2sm&
&at(1,xyzzyaaaa73,2,xyzzyaaac73,1,xyzzyaaad73),6*size_rpsi-5,d2smat1(1&
&,1,2,xyzzyaaac73),1)
endif
enddo
endif
endif
call timer('TRANSFER_D2SMAT',.false.)
end subroutine xyzzyaafx1
subroutine xyzzyaafy1(rele,sele,dbar,logdet,fpeinfo_det,eevecs)
implicit none
integer,intent(in) :: sele(netot)
integer,intent(inout) :: fpeinfo_det(nspin,ndet)
real(dp),intent(in) :: rele(3,netot)
real(dp),intent(in),target :: eevecs(4,netot,netot)
real(dp),intent(inout) :: dbar(nemax,nemax,real1_complex2,nspin,ndet)
complex(dp),intent(inout) :: logdet(nspin,ndet)
integer xyzzyaaaa74,xyzzyaaab74,xyzzyaaac74,xyzzyaaad74
real(dp),pointer :: xyzzyaaae74(:,:)
call timer('EVAL_DBAR',.true.,collapse=.true.)
logdet=czero
fpeinfo_det=xyzzyaaae1
if(.not.pairing_wf)xyzzyaaae74=>eevecs1
do xyzzyaaad74=1,nspin
if(all(missing_det(xyzzyaaad74,:)))cycle
do xyzzyaaaa74=1,nele(xyzzyaaad74)
xyzzyaaab74=which_ii(xyzzyaaaa74,xyzzyaaad74)
if(pairing_wf)xyzzyaaae74=>eevecs(:,:,xyzzyaaab74)
call wfdet(rele(1:3,xyzzyaaab74),sele(xyzzyaaab74),xyzzyaaad74,wfdet_n&
&orb,wfdet_orbmask(1,xyzzyaaad74),.true.,.false.,xyzzyaadc1,xyzzyaadd1&
&,xyzzyaade1,xyzzyaaae74)
call copy_orb_to_det(1,xyzzyaaad74,wfdet_orbmap,xyzzyaadc1,xyzzyaacx1)
do xyzzyaaac74=detstart,detstop
if(missing_det(xyzzyaaad74,xyzzyaaac74))cycle
call dcopy(nele(xyzzyaaad74),xyzzyaacx1(1,1,xyzzyaaac74),1,xyzzyaadb1(&
&1,xyzzyaaaa74,1,xyzzyaaac74),1)
if(complex_wf)call dcopy(nele(xyzzyaaad74),xyzzyaacx1(1,2,xyzzyaaac74)&
&,1,xyzzyaadb1(1,xyzzyaaaa74,2,xyzzyaaac74),1)
enddo
enddo
do xyzzyaaac74=detstart,detstop
if(missing_det(xyzzyaaad74,xyzzyaaac74))cycle
call xyzzyaagc1(nele(xyzzyaaad74),logdet(xyzzyaaad74,xyzzyaaac74),xyzz&
&yaadb1(1,1,1,xyzzyaaac74),dbar(1,1,1,xyzzyaaad74,xyzzyaaac74),fpeinfo&
&_det(xyzzyaaad74,xyzzyaaac74))
enddo
enddo
call timer('EVAL_DBAR',.false.)
end subroutine xyzzyaafy1
subroutine xyzzyaafz1(smat,dbar,logdet,fpeinfo_det)
implicit none
integer,intent(inout) :: fpeinfo_det(nspin,ndet)
real(dp),intent(inout) :: dbar(nemax,nemax,real1_complex2,nspin,ndet),&
&smat(nemax,real1_complex2,ndet,nemax,nspin)
complex(dp),intent(inout) :: logdet(nspin,ndet)
integer xyzzyaaaa75,xyzzyaaab75,xyzzyaaac75
call timer('EVAL_DBAR_FROM_SMAT',.true.,collapse=.true.)
logdet=czero
fpeinfo_det=xyzzyaaae1
do xyzzyaaaa75=1,nspin
do xyzzyaaab75=detstart,detstop
if(missing_det(xyzzyaaaa75,xyzzyaaab75))cycle
do xyzzyaaac75=1,nele(xyzzyaaaa75)
call dcopy(nemax,smat(1,1,xyzzyaaab75,xyzzyaaac75,xyzzyaaaa75),1,xyzzy&
&aadb1(1,xyzzyaaac75,1,1),1)
if(complex_wf)call dcopy(nemax,smat(1,2,xyzzyaaab75,xyzzyaaac75,xyzzya&
&aaa75),1,xyzzyaadb1(1,xyzzyaaac75,2,1),1)
enddo
call xyzzyaagc1(nele(xyzzyaaaa75),logdet(xyzzyaaaa75,xyzzyaaab75),xyzz&
&yaadb1,dbar(1,1,1,xyzzyaaaa75,xyzzyaaab75),fpeinfo_det(xyzzyaaaa75,xy&
&zzyaaab75))
enddo
enddo
call timer('EVAL_DBAR_FROM_SMAT',.false.)
end subroutine xyzzyaafz1
subroutine xyzzyaaga1(smat,ispin,idet,logdet,fpeinfo)
implicit none
integer,intent(in) :: ispin,idet
integer,intent(inout) :: fpeinfo
real(dp),intent(in) :: smat(nemax,real1_complex2,ndet,nemax,nspin)
complex(dp),intent(out) :: logdet
integer xyzzyaaaa76
logical xyzzyaaab76
logdet=czero
fpeinfo=xyzzyaaae1
if(nele(ispin)==0)return
call timer('EVAL_ONEDET',.true.)
if(complex_wf)then
do xyzzyaaaa76=1,nele(ispin)
call dcopy(nemax,smat(1,1,idet,xyzzyaaaa76,ispin),1,xyzzyaadb1(1,xyzzy&
&aaaa76,1,1),1)
call dcopy(nemax,smat(1,2,idet,xyzzyaaaa76,ispin),1,xyzzyaadb1(1,xyzzy&
&aaaa76,2,1),1)
enddo
call lu_decom_cmplx_dz(xyzzyaadb1,xyzzyaadi1,xyzzyaact1,nemax,nele(isp&
&in),xyzzyaaab76)
if(xyzzyaaab76)then
fpeinfo=xyzzyaaaf1
do xyzzyaaaa76=1,nele(ispin)
if(all(xyzzyaadb1(1:nele(ispin),xyzzyaaaa76,1:2,1)==0.d0))then
fpeinfo=xyzzyaaag1
exit
endif
enddo
else
logdet=lu_logdet_cmplx_dz(xyzzyaadb1,xyzzyaact1,nemax,nele(ispin))
endif
else
do xyzzyaaaa76=1,nele(ispin)
call dcopy(nemax,smat(1,1,idet,xyzzyaaaa76,ispin),1,xyzzyaadb1(1,xyzzy&
&aaaa76,1,1),1)
enddo
call lu_decom(xyzzyaadb1,xyzzyaact1,nemax,nele(ispin),xyzzyaaab76)
if(xyzzyaaab76)then
fpeinfo=xyzzyaaaf1
do xyzzyaaaa76=1,nele(ispin)
if(all(xyzzyaadb1(1:nele(ispin),xyzzyaaaa76,1,1)==0.d0))then
fpeinfo=xyzzyaaag1
exit
endif
enddo
else
logdet=lu_logdet(xyzzyaadb1,xyzzyaact1,nemax,nele(ispin))
endif
endif
call timer('EVAL_ONEDET',.false.)
end subroutine xyzzyaaga1
subroutine xyzzyaagb1(ispin,idet,rele,sele,logdet,fpeinfo,eevecs)
implicit none
integer,intent(in) :: ispin,idet,sele(netot)
integer,intent(inout) :: fpeinfo
real(dp),intent(in) :: rele(3,netot)
real(dp),intent(in),target :: eevecs(4,netot,netot)
complex(dp),intent(inout) :: logdet
integer xyzzyaaaa77,xyzzyaaab77,xyzzyaaac77
real(dp),pointer :: xyzzyaaad77(:,:)
logical xyzzyaaae77
if(missing_det(ispin,idet))return
call timer('EVAL_ONEDET_FROM_SCRATCH',.true.,collapse=.true.)
logdet=czero
fpeinfo=xyzzyaaae1
xyzzyaaac77=nele(ispin)
if(.not.pairing_wf)xyzzyaaad77=>eevecs1
do xyzzyaaaa77=1,xyzzyaaac77
xyzzyaaab77=which_ii(xyzzyaaaa77,ispin)
if(pairing_wf)xyzzyaaad77=>eevecs(:,:,xyzzyaaab77)
call wfdet(rele(1:3,xyzzyaaab77),sele(xyzzyaaab77),ispin,wfdet_norb,wf&
&det_orbmask(1,ispin),.true.,.false.,xyzzyaadc1,xyzzyaadd1,xyzzyaade1,&
&xyzzyaaad77)
call copy_orb_to_det(1,ispin,wfdet_orbmap,xyzzyaadc1,xyzzyaacx1)
call dcopy(xyzzyaaac77,xyzzyaacx1(1,1,idet),1,xyzzyaadb1(1,xyzzyaaaa77&
&,1,idet),1)
if(complex_wf)call dcopy(xyzzyaaac77,xyzzyaacx1(1,2,idet),1,xyzzyaadb1&
&(1,xyzzyaaaa77,2,idet),1)
enddo
if(complex_wf)then
call lu_decom_cmplx_dz(xyzzyaadb1(1,1,1,idet),xyzzyaadi1(1,1,idet),xyz&
&zyaact1,nemax,xyzzyaaac77,xyzzyaaae77)
if(xyzzyaaae77)then
fpeinfo=xyzzyaaaf1
do xyzzyaaaa77=1,xyzzyaaac77
if(all(xyzzyaadb1(1:xyzzyaaac77,xyzzyaaaa77,1:2,idet)==0.d0))then
fpeinfo=xyzzyaaag1
exit
endif
enddo
call timer('EVAL_ONEDET_FROM_SCRATCH',.false.)
return
endif
logdet=lu_logdet_cmplx_dz(xyzzyaadb1(1,1,1,idet),xyzzyaact1,nemax,xyzz&
&yaaac77)
else
call lu_decom(xyzzyaadb1(1,1,1,idet),xyzzyaact1,nemax,xyzzyaaac77,xyzz&
&yaaae77)
if(xyzzyaaae77)then
fpeinfo=xyzzyaaaf1
do xyzzyaaaa77=1,xyzzyaaac77
if(all(xyzzyaadb1(1:xyzzyaaac77,xyzzyaaaa77,1,idet)==0.d0))then
fpeinfo=xyzzyaaag1
exit
endif
enddo
call timer('EVAL_ONEDET_FROM_SCRATCH',.false.)
return
endif
logdet=lu_logdet(xyzzyaadb1(1,1,1,idet),xyzzyaact1,nemax,xyzzyaaac77)
endif
call timer('EVAL_ONEDET_FROM_SCRATCH',.false.)
end subroutine xyzzyaagb1
subroutine xyzzyaagc1(n,logdet,smat,dbar,fpeinfo)
implicit none
integer,intent(in) :: n
integer,intent(inout) :: fpeinfo
real(dp),intent(inout) :: dbar(nemax,nemax,real1_complex2),smat(nemax,&
&nemax,real1_complex2)
complex(dp),intent(inout) :: logdet
integer xyzzyaaaa78
logical xyzzyaaab78
call timer('CALDET',.true.)
logdet=czero
fpeinfo=xyzzyaaae1
if(complex_wf)then
call lu_decom_cmplx_dz(smat,xyzzyaadi1,xyzzyaact1,nemax,n,xyzzyaaab78)
if(xyzzyaaab78)then
fpeinfo=xyzzyaaaf1
do xyzzyaaaa78=1,n
if(all(smat(1:n,xyzzyaaaa78,1:2)==0.d0))then
fpeinfo=xyzzyaaag1
exit
endif
enddo
call timer('CALDET',.false.)
return
endif
logdet=lu_logdet_cmplx_dz(smat,xyzzyaact1,nemax,n)
do xyzzyaaaa78=1,n
xyzzyaacv1(1:n,1:2)=0.d0
xyzzyaacv1(xyzzyaaaa78,1)=1.d0
call lu_solve_once_cmplx_dz(smat,xyzzyaadi1,xyzzyaact1,xyzzyaacv1,xyzz&
&yaadh1,nemax,n)
call dcopy(n,xyzzyaacv1,1,dbar(xyzzyaaaa78,1,1),nemax)
call dcopy(n,xyzzyaacv1(1,2),1,dbar(xyzzyaaaa78,1,2),nemax)
enddo
else
call lu_decom(smat,xyzzyaact1,nemax,n,xyzzyaaab78)
if(xyzzyaaab78)then
fpeinfo=xyzzyaaaf1
do xyzzyaaaa78=1,n
if(all(smat(1:n,xyzzyaaaa78,1)==0.d0))then
fpeinfo=xyzzyaaag1
exit
endif
enddo
call timer('CALDET',.false.)
return
endif
logdet=lu_logdet(smat,xyzzyaact1,nemax,n)
dbar=0.d0
do xyzzyaaaa78=1,n
dbar(xyzzyaaaa78,xyzzyaaaa78,1)=1.d0
enddo
call lu_solve_n(smat,xyzzyaact1,dbar,nemax,n,nemax)
do xyzzyaaaa78=1,n-1
call dswap(n-xyzzyaaaa78,dbar(xyzzyaaaa78+1,xyzzyaaaa78,1),1,dbar(xyzz&
&yaaaa78,xyzzyaaaa78+1,1),nemax)
enddo
endif
call timer('CALDET',.false.)
end subroutine xyzzyaagc1
subroutine xyzzyaagd1(ie,ispin,grad,dbar,q,log_pdet,fi_prod_det,fpeinf&
&o_pdet)
implicit none
integer,intent(in) :: ie,ispin,fpeinfo_pdet(ndet)
real(dp),intent(in) :: dbar(nemax,nemax,real1_complex2,nspin,ndet),gra&
&d(3,nemax,real1_complex2,ndet)
real(dp),intent(inout) :: fi_prod_det(3,ndet,nemax,real1_complex2,nspi&
&n)
complex(dp),intent(in) :: q(ndet),log_pdet(ndet)
integer xyzzyaaaa79,xyzzyaaab79,xyzzyaaac79,xyzzyaaad79
real(dp) xyzzyaaae79(3),xyzzyaaaf79(3)
complex(dp) xyzzyaaag79(3)
call timer('COMPUTE_FI_PROD_DET',.true.)
if(complex_wf)then
do xyzzyaaaa79=detstart,detstop
if(fpeinfo_pdet(xyzzyaaaa79)==xyzzyaaag1)then
fi_prod_det(1:3,xyzzyaaaa79,ie,1,ispin)=0.d0
fi_prod_det(1:3,xyzzyaaaa79,ie,2,ispin)=0.d0
cycle
elseif(fpeinfo_pdet(xyzzyaaaa79)==xyzzyaaaf1)then
fi_prod_det(1:3,xyzzyaaaa79,ie,1,ispin)=0.d0
fi_prod_det(1:3,xyzzyaaaa79,ie,2,ispin)=0.d0
cycle
endif
xyzzyaaae79=0.d0
xyzzyaaaf79=0.d0
if(update_by_column(ispin,xyzzyaaaa79))then
do xyzzyaaad79=1,dimensionality
xyzzyaaae79(xyzzyaaad79)=ddot(nele(ispin),dbar(1,ie,1,ispin,xyzzyaaaa7&
&9),1,grad(xyzzyaaad79,1,1,xyzzyaaaa79),3)-ddot(nele(ispin),dbar(1,ie,&
&2,ispin,xyzzyaaaa79),1,grad(xyzzyaaad79,1,2,xyzzyaaaa79),3)
xyzzyaaaf79(xyzzyaaad79)=ddot(nele(ispin),dbar(1,ie,2,ispin,xyzzyaaaa7&
&9),1,grad(xyzzyaaad79,1,1,xyzzyaaaa79),3)+ddot(nele(ispin),dbar(1,ie,&
&1,ispin,xyzzyaaaa79),1,grad(xyzzyaaad79,1,2,xyzzyaaaa79),3)
enddo
else
xyzzyaaab79=upd_spin(ispin,xyzzyaaaa79)
xyzzyaaac79=nuc_nele(xyzzyaaab79)+ie
do xyzzyaaad79=1,dimensionality
xyzzyaaae79(xyzzyaaad79)=-ddot(nele(xyzzyaaab79),dbar(xyzzyaaac79,1,1,&
&xyzzyaaab79,xyzzyaaaa79),nemax,grad(xyzzyaaad79,1,1,xyzzyaaaa79),3)+d&
&dot(nele(xyzzyaaab79),dbar(xyzzyaaac79,1,2,xyzzyaaab79,xyzzyaaaa79),n&
&emax,grad(xyzzyaaad79,1,2,xyzzyaaaa79),3)
xyzzyaaaf79(xyzzyaaad79)=-ddot(nele(xyzzyaaab79),dbar(xyzzyaaac79,1,1,&
&xyzzyaaab79,xyzzyaaaa79),nemax,grad(xyzzyaaad79,1,2,xyzzyaaaa79),3)-d&
&dot(nele(xyzzyaaab79),dbar(xyzzyaaac79,1,2,xyzzyaaab79,xyzzyaaaa79),n&
&emax,grad(xyzzyaaad79,1,1,xyzzyaaaa79),3)
enddo
endif
xyzzyaaag79=cmplx(xyzzyaaae79,xyzzyaaaf79,dp)*exp(log_pdet(xyzzyaaaa79&
&))/q(xyzzyaaaa79)
fi_prod_det(1:3,xyzzyaaaa79,ie,1,ispin)=dble(xyzzyaaag79)
fi_prod_det(1:3,xyzzyaaaa79,ie,2,ispin)=aimag(xyzzyaaag79)
enddo
else
do xyzzyaaaa79=detstart,detstop
if(fpeinfo_pdet(xyzzyaaaa79)==xyzzyaaag1)then
fi_prod_det(1:3,xyzzyaaaa79,ie,1,ispin)=0.d0
cycle
elseif(fpeinfo_pdet(xyzzyaaaa79)==xyzzyaaaf1)then
fi_prod_det(1:3,xyzzyaaaa79,ie,1,ispin)=0.d0
cycle
endif
xyzzyaaae79=0.d0
if(update_by_column(ispin,xyzzyaaaa79))then
do xyzzyaaad79=1,dimensionality
xyzzyaaae79(xyzzyaaad79)=ddot(nele(ispin),dbar(1,ie,1,ispin,xyzzyaaaa7&
&9),1,grad(xyzzyaaad79,1,1,xyzzyaaaa79),3)
enddo
else
xyzzyaaab79=upd_spin(ispin,xyzzyaaaa79)
xyzzyaaac79=nuc_nele(xyzzyaaab79)+ie
do xyzzyaaad79=1,dimensionality
xyzzyaaae79(xyzzyaaad79)=-ddot(nele(xyzzyaaab79),dbar(xyzzyaaac79,1,1,&
&xyzzyaaab79,xyzzyaaaa79),nemax,grad(xyzzyaaad79,1,1,xyzzyaaaa79),3)
enddo
endif
fi_prod_det(1:3,xyzzyaaaa79,ie,1,ispin)=xyzzyaaae79*dble(exp(log_pdet(&
&xyzzyaaaa79))/q(xyzzyaaaa79))
enddo
endif
call timer('COMPUTE_FI_PROD_DET',.false.)
end subroutine xyzzyaagd1
subroutine xyzzyaage1(ii,logwfn,fpeinfo_wfn,fi_prod_det,fpeinfo_pdet,f&
&idet)
implicit none
integer,intent(in) :: ii,fpeinfo_wfn,fpeinfo_pdet(ndet)
real(dp),intent(in) :: fi_prod_det(3,ndet,nemax,real1_complex2,nspin)
complex(dp),intent(in) :: logwfn
complex(dp),intent(out) :: fidet(3)
integer xyzzyaaaa80,xyzzyaaab80,xyzzyaaac80
complex(dp) xyzzyaaad80(3)
fidet=czero
if(fpeinfo_wfn/=xyzzyaaae1)return
call timer('COMPUTE_FIDET',.true.)
xyzzyaaaa80=which_ie(ii)
xyzzyaaab80=which_spin(ii)
xyzzyaaad80=czero
if(complex_wf)then
do xyzzyaaac80=detstart,detstop
if(fpeinfo_pdet(xyzzyaaac80)==xyzzyaaag1)cycle
xyzzyaaad80=xyzzyaaad80+detcoef(xyzzyaaac80)*cmplx(fi_prod_det(1:3,xyz&
&zyaaac80,xyzzyaaaa80,1,xyzzyaaab80),fi_prod_det(1:3,xyzzyaaac80,xyzzy&
&aaaa80,2,xyzzyaaab80),dp)
enddo
else
do xyzzyaaac80=detstart,detstop
if(fpeinfo_pdet(xyzzyaaac80)==xyzzyaaag1)cycle
xyzzyaaad80=xyzzyaaad80+cmplx(detcoef(xyzzyaaac80)*fi_prod_det(1:3,xyz&
&zyaaac80,xyzzyaaaa80,1,xyzzyaaab80),0.d0,dp)
enddo
endif
fidet=xyzzyaaad80/exp(logwfn)
call timer('COMPUTE_FIDET',.false.)
end subroutine xyzzyaage1
subroutine xyzzyaagf1(ie,ispin,lap,dbar,q,log_pdet,prod_lapdet,fpeinfo&
&_pdet)
implicit none
integer,intent(in) :: ie,ispin,fpeinfo_pdet(ndet)
real(dp),intent(in) :: dbar(nemax,nemax,real1_complex2,nspin,ndet),lap&
&(nemax,real1_complex2,ndet)
real(dp),intent(inout) :: prod_lapdet(ndet,nemax,real1_complex2,nspin)
complex(dp),intent(in) :: log_pdet(ndet),q(ndet)
integer xyzzyaaaa81,xyzzyaaab81,xyzzyaaac81
real(dp) xyzzyaaad81,xyzzyaaae81
complex(dp) xyzzyaaaf81
call timer('COMPUTE_PROD_LAPDET',.true.)
if(complex_wf)then
do xyzzyaaaa81=detstart,detstop
if(fpeinfo_pdet(xyzzyaaaa81)==xyzzyaaag1)then
prod_lapdet(xyzzyaaaa81,ie,1,ispin)=0.d0
prod_lapdet(xyzzyaaaa81,ie,2,ispin)=0.d0
cycle
elseif(fpeinfo_pdet(xyzzyaaaa81)==xyzzyaaaf1)then
prod_lapdet(xyzzyaaaa81,ie,1,ispin)=0.d0
prod_lapdet(xyzzyaaaa81,ie,2,ispin)=0.d0
cycle
endif
if(update_by_column(ispin,xyzzyaaaa81))then
xyzzyaaad81=ddot(nele(ispin),dbar(1,ie,1,ispin,xyzzyaaaa81),1,lap(1,1,&
&xyzzyaaaa81),1)-ddot(nele(ispin),dbar(1,ie,2,ispin,xyzzyaaaa81),1,lap&
&(1,2,xyzzyaaaa81),1)
xyzzyaaae81=ddot(nele(ispin),dbar(1,ie,2,ispin,xyzzyaaaa81),1,lap(1,1,&
&xyzzyaaaa81),1)+ddot(nele(ispin),dbar(1,ie,1,ispin,xyzzyaaaa81),1,lap&
&(1,2,xyzzyaaaa81),1)
else
xyzzyaaab81=upd_spin(ispin,xyzzyaaaa81)
xyzzyaaac81=nuc_nele(xyzzyaaab81)+ie
xyzzyaaad81=ddot(nele(xyzzyaaab81),dbar(xyzzyaaac81,1,1,xyzzyaaab81,xy&
&zzyaaaa81),nemax,lap(1,1,xyzzyaaaa81),1)-ddot(nele(xyzzyaaab81),dbar(&
&xyzzyaaac81,1,2,xyzzyaaab81,xyzzyaaaa81),nemax,lap(1,2,xyzzyaaaa81),1&
&)
xyzzyaaae81=ddot(nele(xyzzyaaab81),dbar(xyzzyaaac81,1,2,xyzzyaaab81,xy&
&zzyaaaa81),nemax,lap(1,1,xyzzyaaaa81),1)+ddot(nele(xyzzyaaab81),dbar(&
&xyzzyaaac81,1,1,xyzzyaaab81,xyzzyaaaa81),nemax,lap(1,2,xyzzyaaaa81),1&
&)
endif
xyzzyaaaf81=cmplx(xyzzyaaad81,xyzzyaaae81,dp)*exp(log_pdet(xyzzyaaaa81&
&))/q(xyzzyaaaa81)
prod_lapdet(xyzzyaaaa81,ie,1,ispin)=dble(xyzzyaaaf81)
prod_lapdet(xyzzyaaaa81,ie,2,ispin)=aimag(xyzzyaaaf81)
enddo
else
do xyzzyaaaa81=detstart,detstop
if(fpeinfo_pdet(xyzzyaaaa81)==xyzzyaaag1)then
prod_lapdet(xyzzyaaaa81,ie,1,ispin)=0.d0
cycle
elseif(fpeinfo_pdet(xyzzyaaaa81)==xyzzyaaaf1)then
prod_lapdet(xyzzyaaaa81,ie,1,ispin)=0.d0
cycle
endif
if(update_by_column(ispin,xyzzyaaaa81))then
xyzzyaaad81=ddot(nele(ispin),dbar(1,ie,1,ispin,xyzzyaaaa81),1,lap(1,1,&
&xyzzyaaaa81),1)
else
xyzzyaaab81=upd_spin(ispin,xyzzyaaaa81)
xyzzyaaac81=nuc_nele(xyzzyaaab81)+ie
xyzzyaaad81=ddot(nele(xyzzyaaab81),dbar(xyzzyaaac81,1,1,xyzzyaaab81,xy&
&zzyaaaa81),nemax,lap(1,1,xyzzyaaaa81),1)
endif
prod_lapdet(xyzzyaaaa81,ie,1,ispin)=xyzzyaaad81*dble(exp(log_pdet(xyzz&
&yaaaa81))/q(xyzzyaaaa81))
enddo
endif
call timer('COMPUTE_PROD_LAPDET',.false.)
end subroutine xyzzyaagf1
subroutine xyzzyaagg1(ii,logwfn,fpeinfo_wfn,prod_lapdet,fpeinfo_pdet,f&
&idet,lapdet)
implicit none
integer,intent(in) :: ii,fpeinfo_wfn,fpeinfo_pdet(ndet)
real(dp),intent(in) :: prod_lapdet(ndet,nemax,real1_complex2,nspin)
complex(dp),intent(in) :: fidet(3),logwfn
complex(dp),intent(out) :: lapdet
integer xyzzyaaaa82,xyzzyaaab82,xyzzyaaac82
real(dp) xyzzyaaad82
complex(dp) xyzzyaaae82
lapdet=czero
if(fpeinfo_wfn/=xyzzyaaae1)return
call timer('COMPUTE_LAPDET',.true.)
xyzzyaaaa82=which_ie(ii)
xyzzyaaab82=which_spin(ii)
xyzzyaaae82=czero
if(complex_wf)then
do xyzzyaaac82=detstart,detstop
if(fpeinfo_pdet(xyzzyaaac82)==xyzzyaaag1)cycle
xyzzyaaae82=xyzzyaaae82+detcoef(xyzzyaaac82)*cmplx(prod_lapdet(xyzzyaa&
&ac82,xyzzyaaaa82,1,xyzzyaaab82),prod_lapdet(xyzzyaaac82,xyzzyaaaa82,2&
&,xyzzyaaab82),dp)
enddo
else
do xyzzyaaac82=detstart,detstop
if(fpeinfo_pdet(xyzzyaaac82)==xyzzyaaag1)cycle
xyzzyaaae82=xyzzyaaae82+cmplx(detcoef(xyzzyaaac82)*prod_lapdet(xyzzyaa&
&ac82,xyzzyaaaa82,1,xyzzyaaab82),0.d0,dp)
enddo
endif
xyzzyaaad82=dble(fidet(1))**2+dble(fidet(2))**2+dble(fidet(3))**2
if(complex_wf)xyzzyaaad82=xyzzyaaad82+aimag(fidet(1))**2+aimag(fidet(2&
&))**2+aimag(fidet(3))**2
lapdet=xyzzyaaae82/exp(logwfn)-cmplx(xyzzyaaad82,0.d0,dp)
call timer('COMPUTE_LAPDET',.false.)
end subroutine xyzzyaagg1
subroutine xyzzyaagh1(dbar,c_pdet,logwfn,dsmat,farray_det,farray,fpein&
&fo_pdet)
implicit none
integer,intent(in) :: fpeinfo_pdet(ndet)
real(dp),intent(in) :: dbar(nemax,nemax,real1_complex2,nspin,ndet),dsm&
&at(3,nemax,real1_complex2,ndet,nemax,nspin)
real(dp),intent(inout) :: farray_det(3,real1_complex2,ndet,netot),farr&
&ay(3,real1_complex2,netot)
complex(dp),intent(in) :: c_pdet(ndet),logwfn
integer xyzzyaaaa83,xyzzyaaab83,xyzzyaaac83,xyzzyaaad83,xyzzyaaae83,xy&
&zzyaaaf83,xyzzyaaag83
real(dp) xyzzyaaah83,xyzzyaaai83,xyzzyaaaj83(ndet),xyzzyaaak83(ndet),x&
&yzzyaaal83,xyzzyaaam83
complex(dp) xyzzyaaan83
call timer('EVAL_FARRAY',.true.)
xyzzyaaaj83=dble(c_pdet)
xyzzyaaak83=aimag(c_pdet)
xyzzyaaan83=c_one/exp(logwfn)
xyzzyaaah83=dble(xyzzyaaan83)
xyzzyaaai83=aimag(xyzzyaaan83)
if(complex_wf)then
do xyzzyaaaa83=1,netot
xyzzyaaac83=which_spin(xyzzyaaaa83)
xyzzyaaab83=which_ie(xyzzyaaaa83)
xyzzyaaag83=nele(xyzzyaaac83)
do xyzzyaaad83=detstart,detstop
if(fpeinfo_pdet(xyzzyaaad83)==xyzzyaaag1)then
farray_det(1:3,1,xyzzyaaad83,xyzzyaaaa83)=0.d0
farray_det(1:3,2,xyzzyaaad83,xyzzyaaaa83)=0.d0
cycle
elseif(fpeinfo_pdet(xyzzyaaad83)==xyzzyaaaf1)then
farray_det(1:3,1,xyzzyaaad83,xyzzyaaaa83)=0.d0
farray_det(1:3,2,xyzzyaaad83,xyzzyaaaa83)=0.d0
cycle
endif
if(update_by_column(xyzzyaaac83,xyzzyaaad83))then
do xyzzyaaae83=1,dimensionality
farray_det(xyzzyaaae83,1,xyzzyaaad83,xyzzyaaaa83)=ddot(xyzzyaaag83,dba&
&r(1,xyzzyaaab83,1,xyzzyaaac83,xyzzyaaad83),1,dsmat(xyzzyaaae83,1,1,xy&
&zzyaaad83,xyzzyaaab83,xyzzyaaac83),3)-ddot(xyzzyaaag83,dbar(1,xyzzyaa&
&ab83,2,xyzzyaaac83,xyzzyaaad83),1,dsmat(xyzzyaaae83,1,2,xyzzyaaad83,x&
&yzzyaaab83,xyzzyaaac83),3)
farray_det(xyzzyaaae83,2,xyzzyaaad83,xyzzyaaaa83)=ddot(xyzzyaaag83,dba&
&r(1,xyzzyaaab83,1,xyzzyaaac83,xyzzyaaad83),1,dsmat(xyzzyaaae83,1,2,xy&
&zzyaaad83,xyzzyaaab83,xyzzyaaac83),3)+ddot(xyzzyaaag83,dbar(1,xyzzyaa&
&ab83,2,xyzzyaaac83,xyzzyaaad83),1,dsmat(xyzzyaaae83,1,1,xyzzyaaad83,x&
&yzzyaaab83,xyzzyaaac83),3)
enddo
else
xyzzyaaaf83=upd_spin(xyzzyaaac83,xyzzyaaad83)
do xyzzyaaae83=1,dimensionality
farray_det(xyzzyaaae83,1,xyzzyaaad83,xyzzyaaaa83)=-ddot(xyzzyaaag83,db&
&ar(xyzzyaaab83,1,1,xyzzyaaaf83,xyzzyaaad83),nemax,dsmat(xyzzyaaae83,x&
&yzzyaaab83,1,xyzzyaaad83,1,xyzzyaaaf83),size_grad)+ddot(xyzzyaaag83,d&
&bar(xyzzyaaab83,1,2,xyzzyaaaf83,xyzzyaaad83),nemax,dsmat(xyzzyaaae83,&
&xyzzyaaab83,2,xyzzyaaad83,1,xyzzyaaaf83),size_grad)
farray_det(xyzzyaaae83,2,xyzzyaaad83,xyzzyaaaa83)=-ddot(xyzzyaaag83,db&
&ar(xyzzyaaab83,1,1,xyzzyaaaf83,xyzzyaaad83),nemax,dsmat(xyzzyaaae83,x&
&yzzyaaab83,2,xyzzyaaad83,1,xyzzyaaaf83),size_grad)-ddot(xyzzyaaag83,d&
&bar(xyzzyaaab83,1,2,xyzzyaaaf83,xyzzyaaad83),nemax,dsmat(xyzzyaaae83,&
&xyzzyaaab83,1,xyzzyaaad83,1,xyzzyaaaf83),size_grad)
enddo
endif
enddo
do xyzzyaaae83=1,dimensionality
xyzzyaaal83=ddot(ndet,xyzzyaaaj83,1,farray_det(xyzzyaaae83,1,1,xyzzyaa&
&aa83),6)-ddot(ndet,xyzzyaaak83,1,farray_det(xyzzyaaae83,2,1,xyzzyaaaa&
&83),6)
xyzzyaaam83=ddot(ndet,xyzzyaaaj83,1,farray_det(xyzzyaaae83,2,1,xyzzyaa&
&aa83),6)+ddot(ndet,xyzzyaaak83,1,farray_det(xyzzyaaae83,1,1,xyzzyaaaa&
&83),6)
farray(xyzzyaaae83,1,xyzzyaaaa83)=xyzzyaaal83*xyzzyaaah83-xyzzyaaam83*&
&xyzzyaaai83
farray(xyzzyaaae83,2,xyzzyaaaa83)=xyzzyaaal83*xyzzyaaai83+xyzzyaaam83*&
&xyzzyaaah83
enddo
enddo
else
do xyzzyaaaa83=1,netot
xyzzyaaac83=which_spin(xyzzyaaaa83)
xyzzyaaab83=which_ie(xyzzyaaaa83)
xyzzyaaag83=nele(xyzzyaaac83)
do xyzzyaaad83=detstart,detstop
if(fpeinfo_pdet(xyzzyaaad83)==xyzzyaaag1)then
farray_det(1:3,1,xyzzyaaad83,xyzzyaaaa83)=0.d0
cycle
elseif(fpeinfo_pdet(xyzzyaaad83)==xyzzyaaaf1)then
farray_det(1:3,1,xyzzyaaad83,xyzzyaaaa83)=0.d0
cycle
endif
if(update_by_column(xyzzyaaac83,xyzzyaaad83))then
do xyzzyaaae83=1,dimensionality
farray_det(xyzzyaaae83,1,xyzzyaaad83,xyzzyaaaa83)=ddot(xyzzyaaag83,dba&
&r(1,xyzzyaaab83,1,xyzzyaaac83,xyzzyaaad83),1,dsmat(xyzzyaaae83,1,1,xy&
&zzyaaad83,xyzzyaaab83,xyzzyaaac83),3)
enddo
else
xyzzyaaaf83=upd_spin(xyzzyaaac83,xyzzyaaad83)
do xyzzyaaae83=1,dimensionality
farray_det(xyzzyaaae83,1,xyzzyaaad83,xyzzyaaaa83)=-ddot(xyzzyaaag83,db&
&ar(xyzzyaaab83,1,1,xyzzyaaaf83,xyzzyaaad83),nemax,dsmat(xyzzyaaae83,x&
&yzzyaaab83,1,xyzzyaaad83,1,xyzzyaaaf83),size_grad)
enddo
endif
enddo
do xyzzyaaae83=1,dimensionality
farray(xyzzyaaae83,1,xyzzyaaaa83)=ddot(ndet,xyzzyaaaj83,1,farray_det(x&
&yzzyaaae83,1,1,xyzzyaaaa83),3)*xyzzyaaah83
enddo
enddo
endif
call timer('EVAL_FARRAY',.false.)
end subroutine xyzzyaagh1
subroutine xyzzyaagi1(dbar,c_pdet,logwfn,dsmat,m,rmap,farray_det,farra&
&y,fpeinfo_pdet)
implicit none
integer,intent(in) :: m,rmap(m),fpeinfo_pdet(ndet)
real(dp),intent(in) :: dbar(nemax,nemax,real1_complex2,nspin,ndet),dsm&
&at(3,nemax,real1_complex2,ndet,nemax,nspin)
real(dp),intent(inout) :: farray_det(3,real1_complex2,ndet,netot),farr&
&ay(3,real1_complex2,netot)
complex(dp),intent(in) :: c_pdet(ndet),logwfn
integer xyzzyaaaa84,xyzzyaaab84,xyzzyaaac84,xyzzyaaad84,xyzzyaaae84,xy&
&zzyaaaf84,xyzzyaaag84,xyzzyaaah84
real(dp) xyzzyaaai84,xyzzyaaaj84,xyzzyaaak84(ndet),xyzzyaaal84(ndet),x&
&yzzyaaam84,xyzzyaaan84
complex(dp) xyzzyaaao84
call timer('EVAL_FARRAY_II',.true.)
xyzzyaaak84=dble(c_pdet)
xyzzyaaal84=aimag(c_pdet)
xyzzyaaao84=c_one/exp(logwfn)
xyzzyaaai84=dble(xyzzyaaao84)
xyzzyaaaj84=aimag(xyzzyaaao84)
if(complex_wf)then
do xyzzyaaaa84=1,m
xyzzyaaab84=rmap(xyzzyaaaa84)
xyzzyaaad84=which_spin(xyzzyaaab84)
xyzzyaaac84=which_ie(xyzzyaaab84)
xyzzyaaah84=nele(xyzzyaaad84)
do xyzzyaaae84=detstart,detstop
if(fpeinfo_pdet(xyzzyaaae84)==xyzzyaaag1)then
farray_det(1:3,1,xyzzyaaae84,xyzzyaaab84)=0.d0
farray_det(1:3,2,xyzzyaaae84,xyzzyaaab84)=0.d0
cycle
elseif(fpeinfo_pdet(xyzzyaaae84)==xyzzyaaaf1)then
farray_det(1:3,1,xyzzyaaae84,xyzzyaaab84)=0.d0
farray_det(1:3,2,xyzzyaaae84,xyzzyaaab84)=0.d0
cycle
endif
if(update_by_column(xyzzyaaad84,xyzzyaaae84))then
do xyzzyaaaf84=1,dimensionality
farray_det(xyzzyaaaf84,1,xyzzyaaae84,xyzzyaaab84)=ddot(xyzzyaaah84,dba&
&r(1,xyzzyaaac84,1,xyzzyaaad84,xyzzyaaae84),1,dsmat(xyzzyaaaf84,1,1,xy&
&zzyaaae84,xyzzyaaac84,xyzzyaaad84),3)-ddot(xyzzyaaah84,dbar(1,xyzzyaa&
&ac84,2,xyzzyaaad84,xyzzyaaae84),1,dsmat(xyzzyaaaf84,1,2,xyzzyaaae84,x&
&yzzyaaac84,xyzzyaaad84),3)
farray_det(xyzzyaaaf84,2,xyzzyaaae84,xyzzyaaab84)=ddot(xyzzyaaah84,dba&
&r(1,xyzzyaaac84,1,xyzzyaaad84,xyzzyaaae84),1,dsmat(xyzzyaaaf84,1,2,xy&
&zzyaaae84,xyzzyaaac84,xyzzyaaad84),3)+ddot(xyzzyaaah84,dbar(1,xyzzyaa&
&ac84,2,xyzzyaaad84,xyzzyaaae84),1,dsmat(xyzzyaaaf84,1,1,xyzzyaaae84,x&
&yzzyaaac84,xyzzyaaad84),3)
enddo
else
xyzzyaaag84=upd_spin(xyzzyaaad84,xyzzyaaae84)
do xyzzyaaaf84=1,dimensionality
farray_det(xyzzyaaaf84,1,xyzzyaaae84,xyzzyaaab84)=-ddot(xyzzyaaah84,db&
&ar(xyzzyaaac84,1,1,xyzzyaaag84,xyzzyaaae84),nemax,dsmat(xyzzyaaaf84,x&
&yzzyaaac84,1,xyzzyaaae84,1,xyzzyaaag84),size_grad)+ddot(xyzzyaaah84,d&
&bar(xyzzyaaac84,1,2,xyzzyaaag84,xyzzyaaae84),nemax,dsmat(xyzzyaaaf84,&
&xyzzyaaac84,2,xyzzyaaae84,1,xyzzyaaag84),size_grad)
farray_det(xyzzyaaaf84,2,xyzzyaaae84,xyzzyaaab84)=-ddot(xyzzyaaah84,db&
&ar(xyzzyaaac84,1,1,xyzzyaaag84,xyzzyaaae84),nemax,dsmat(xyzzyaaaf84,x&
&yzzyaaac84,2,xyzzyaaae84,1,xyzzyaaag84),size_grad)-ddot(xyzzyaaah84,d&
&bar(xyzzyaaac84,1,2,xyzzyaaag84,xyzzyaaae84),nemax,dsmat(xyzzyaaaf84,&
&xyzzyaaac84,1,xyzzyaaae84,1,xyzzyaaag84),size_grad)
enddo
endif
enddo
do xyzzyaaaf84=1,dimensionality
xyzzyaaam84=ddot(ndet,xyzzyaaak84,1,farray_det(xyzzyaaaf84,1,1,xyzzyaa&
&ab84),6)-ddot(ndet,xyzzyaaal84,1,farray_det(xyzzyaaaf84,2,1,xyzzyaaab&
&84),6)
xyzzyaaan84=ddot(ndet,xyzzyaaak84,1,farray_det(xyzzyaaaf84,2,1,xyzzyaa&
&ab84),6)+ddot(ndet,xyzzyaaal84,1,farray_det(xyzzyaaaf84,1,1,xyzzyaaab&
&84),6)
farray(xyzzyaaaf84,1,xyzzyaaab84)=xyzzyaaam84*xyzzyaaai84-xyzzyaaan84*&
&xyzzyaaaj84
farray(xyzzyaaaf84,2,xyzzyaaab84)=xyzzyaaam84*xyzzyaaaj84+xyzzyaaan84*&
&xyzzyaaai84
enddo
enddo
else
do xyzzyaaaa84=1,m
xyzzyaaab84=rmap(xyzzyaaaa84)
xyzzyaaad84=which_spin(xyzzyaaab84)
xyzzyaaac84=which_ie(xyzzyaaab84)
xyzzyaaah84=nele(xyzzyaaad84)
do xyzzyaaae84=detstart,detstop
if(fpeinfo_pdet(xyzzyaaae84)==xyzzyaaag1)then
farray_det(1:3,1,xyzzyaaae84,xyzzyaaab84)=0.d0
cycle
elseif(fpeinfo_pdet(xyzzyaaae84)==xyzzyaaaf1)then
farray_det(1:3,1,xyzzyaaae84,xyzzyaaab84)=0.d0
cycle
endif
if(update_by_column(xyzzyaaad84,xyzzyaaae84))then
do xyzzyaaaf84=1,dimensionality
farray_det(xyzzyaaaf84,1,xyzzyaaae84,xyzzyaaab84)=ddot(xyzzyaaah84,dba&
&r(1,xyzzyaaac84,1,xyzzyaaad84,xyzzyaaae84),1,dsmat(xyzzyaaaf84,1,1,xy&
&zzyaaae84,xyzzyaaac84,xyzzyaaad84),3)
enddo
else
xyzzyaaag84=upd_spin(xyzzyaaad84,xyzzyaaae84)
do xyzzyaaaf84=1,dimensionality
farray_det(xyzzyaaaf84,1,xyzzyaaae84,xyzzyaaab84)=-ddot(xyzzyaaah84,db&
&ar(xyzzyaaac84,1,1,xyzzyaaag84,xyzzyaaae84),nemax,dsmat(xyzzyaaaf84,x&
&yzzyaaac84,1,xyzzyaaae84,1,xyzzyaaag84),size_grad)
enddo
endif
enddo
do xyzzyaaaf84=1,dimensionality
farray(xyzzyaaaf84,1,xyzzyaaab84)=ddot(ndet,xyzzyaaak84,1,farray_det(x&
&yzzyaaaf84,1,1,xyzzyaaab84),3)*xyzzyaaai84
enddo
enddo
endif
call timer('EVAL_FARRAY_II',.false.)
end subroutine xyzzyaagi1
subroutine xyzzyaagj1(dbar,c_pdet,logwfn,dsmat,d2smat,bf_m2,bf_rmap2,f&
&array_det,harray,fpeinfo_pdet)
implicit none
integer,intent(in) :: bf_m2(netot),bf_rmap2(netot,netot),fpeinfo_pdet(&
&ndet)
real(dp),intent(in) :: dbar(nemax,nemax,real1_complex2,nspin,ndet),dsm&
&at(3,nemax,real1_complex2,ndet,nemax,nspin),d2smat(6,nemax,real1_comp&
&lex2,ndet,nemax,nspin),farray_det(3,real1_complex2,ndet,netot)
real(dp),intent(inout) :: harray(3,3,real1_complex2,netot,netot)
complex(dp),intent(in) :: c_pdet(ndet),logwfn
integer xyzzyaaaa85,xyzzyaaab85,xyzzyaaac85,xyzzyaaad85,xyzzyaaae85,xy&
&zzyaaaf85,xyzzyaaag85,xyzzyaaah85,xyzzyaaai85,xyzzyaaaj85,xyzzyaaak85&
&,xyzzyaaal85,xyzzyaaam85,xyzzyaaan85,xyzzyaaao85,xyzzyaaap85
real(dp) xyzzyaaaq85,xyzzyaaar85,xyzzyaaas85,xyzzyaaat85,xyzzyaaau85,x&
&yzzyaaav85,xyzzyaaaw85(ndet),xyzzyaaax85(ndet),xyzzyaaay85,xyzzyaaaz8&
&5,xyzzyaaba85,xyzzyaabb85
complex(dp) xyzzyaabc85
call timer('EVAL_HARRAY',.true.)
xyzzyaaaw85=dble(c_pdet)
xyzzyaaax85=aimag(c_pdet)
xyzzyaabc85=c_one/exp(logwfn)
xyzzyaaau85=dble(xyzzyaabc85)
xyzzyaaav85=aimag(xyzzyaabc85)
if(complex_wf)then
do xyzzyaaaa85=1,netot
xyzzyaaac85=which_spin(xyzzyaaaa85)
xyzzyaaab85=which_ie(xyzzyaaaa85)
xyzzyaaao85=nele(xyzzyaaac85)
do xyzzyaaad85=1,bf_m2(xyzzyaaaa85)
xyzzyaaae85=bf_rmap2(xyzzyaaad85,xyzzyaaaa85)
if(xyzzyaaae85<xyzzyaaaa85)cycle
xyzzyaaag85=which_spin(xyzzyaaae85)
xyzzyaaaf85=which_ie(xyzzyaaae85)
xyzzyaaap85=nele(xyzzyaaag85)
do xyzzyaaai85=1,dimensionality
do xyzzyaaaj85=1,dimensionality
xyzzyaaak85=which_d2index(xyzzyaaai85,xyzzyaaaj85)
xyzzyaaaq85=0.d0
xyzzyaaar85=0.d0
do xyzzyaaah85=detstart,detstop
if(fpeinfo_pdet(xyzzyaaah85)==xyzzyaaag1)then
cycle
elseif(fpeinfo_pdet(xyzzyaaah85)==xyzzyaaaf1)then
cycle
endif
xyzzyaaas85=0.d0
xyzzyaaat85=0.d0
if(xyzzyaaaa85==xyzzyaaae85)then
if(update_by_column(xyzzyaaac85,xyzzyaaah85))then
xyzzyaaay85=ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,d2smat(xyzzyaaak85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),6)&
&-ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,2,xyzzyaaac85,xyzzyaaah85),1,d2s&
&mat(xyzzyaaak85,1,2,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),6)
xyzzyaaaz85=ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,d2smat(xyzzyaaak85,1,2,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),6)&
&+ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,2,xyzzyaaac85,xyzzyaaah85),1,d2s&
&mat(xyzzyaaak85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),6)
else
xyzzyaaal85=upd_spin(xyzzyaaac85,xyzzyaaah85)
xyzzyaaay85=ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,1,xyzzyaaal85,xyzzyaaa&
&h85),nemax,d2smat(xyzzyaaak85,xyzzyaaab85,1,xyzzyaaah85,1,xyzzyaaal85&
&),size_sderivs)-ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,2,xyzzyaaal85,xyz&
&zyaaah85),nemax,d2smat(xyzzyaaak85,xyzzyaaab85,2,xyzzyaaah85,1,xyzzya&
&aal85),size_sderivs)
xyzzyaaaz85=ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,1,xyzzyaaal85,xyzzyaaa&
&h85),nemax,d2smat(xyzzyaaak85,xyzzyaaab85,2,xyzzyaaah85,1,xyzzyaaal85&
&),size_sderivs)+ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,2,xyzzyaaal85,xyz&
&zyaaah85),nemax,d2smat(xyzzyaaak85,xyzzyaaab85,1,xyzzyaaah85,1,xyzzya&
&aal85),size_sderivs)
endif
xyzzyaaas85=xyzzyaaas85+xyzzyaaay85
xyzzyaaat85=xyzzyaaat85+xyzzyaaaz85
else
xyzzyaaas85=xyzzyaaas85+farray_det(xyzzyaaai85,1,xyzzyaaah85,xyzzyaaaa&
&85)*farray_det(xyzzyaaaj85,1,xyzzyaaah85,xyzzyaaae85)-farray_det(xyzz&
&yaaai85,2,xyzzyaaah85,xyzzyaaaa85)*farray_det(xyzzyaaaj85,2,xyzzyaaah&
&85,xyzzyaaae85)
xyzzyaaat85=xyzzyaaat85+farray_det(xyzzyaaai85,1,xyzzyaaah85,xyzzyaaaa&
&85)*farray_det(xyzzyaaaj85,2,xyzzyaaah85,xyzzyaaae85)+farray_det(xyzz&
&yaaai85,2,xyzzyaaah85,xyzzyaaaa85)*farray_det(xyzzyaaaj85,1,xyzzyaaah&
&85,xyzzyaaae85)
if(xyzzyaaac85==xyzzyaaag85)then
if(update_by_column(xyzzyaaac85,xyzzyaaah85))then
xyzzyaaay85=ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaaj85,1,1,xyzzyaaah85,xyzzyaaaf85,xyzzyaaac85),3)-&
&ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,2,xyzzyaaac85,xyzzyaaah85),1,dsma&
&t(xyzzyaaaj85,1,2,xyzzyaaah85,xyzzyaaaf85,xyzzyaaac85),3)
xyzzyaaaz85=ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaaj85,1,2,xyzzyaaah85,xyzzyaaaf85,xyzzyaaac85),3)+&
&ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,2,xyzzyaaac85,xyzzyaaah85),1,dsma&
&t(xyzzyaaaj85,1,1,xyzzyaaah85,xyzzyaaaf85,xyzzyaaac85),3)
xyzzyaaba85=ddot(xyzzyaaao85,dbar(1,xyzzyaaaf85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaai85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)-&
&ddot(xyzzyaaao85,dbar(1,xyzzyaaaf85,2,xyzzyaaac85,xyzzyaaah85),1,dsma&
&t(xyzzyaaai85,1,2,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)
xyzzyaabb85=ddot(xyzzyaaao85,dbar(1,xyzzyaaaf85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaai85,1,2,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)+&
&ddot(xyzzyaaao85,dbar(1,xyzzyaaaf85,2,xyzzyaaac85,xyzzyaaah85),1,dsma&
&t(xyzzyaaai85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)
else
xyzzyaaal85=upd_spin(xyzzyaaac85,xyzzyaaah85)
xyzzyaaay85=-ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,1,xyzzyaaal85,xyzzyaa&
&ah85),nemax,dsmat(xyzzyaaaj85,xyzzyaaaf85,1,xyzzyaaah85,1,xyzzyaaal85&
&),size_grad)+ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,2,xyzzyaaal85,xyzzya&
&aah85),nemax,dsmat(xyzzyaaaj85,xyzzyaaaf85,2,xyzzyaaah85,1,xyzzyaaal8&
&5),size_grad)
xyzzyaaaz85=-ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,1,xyzzyaaal85,xyzzyaa&
&ah85),nemax,dsmat(xyzzyaaaj85,xyzzyaaaf85,2,xyzzyaaah85,1,xyzzyaaal85&
&),size_grad)-ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,2,xyzzyaaal85,xyzzya&
&aah85),nemax,dsmat(xyzzyaaaj85,xyzzyaaaf85,1,xyzzyaaah85,1,xyzzyaaal8&
&5),size_grad)
xyzzyaaba85=-ddot(xyzzyaaao85,dbar(xyzzyaaaf85,1,1,xyzzyaaal85,xyzzyaa&
&ah85),nemax,dsmat(xyzzyaaai85,xyzzyaaab85,1,xyzzyaaah85,1,xyzzyaaal85&
&),size_grad)+ddot(xyzzyaaao85,dbar(xyzzyaaaf85,1,2,xyzzyaaal85,xyzzya&
&aah85),nemax,dsmat(xyzzyaaai85,xyzzyaaab85,2,xyzzyaaah85,1,xyzzyaaal8&
&5),size_grad)
xyzzyaabb85=-ddot(xyzzyaaao85,dbar(xyzzyaaaf85,1,1,xyzzyaaal85,xyzzyaa&
&ah85),nemax,dsmat(xyzzyaaai85,xyzzyaaab85,2,xyzzyaaah85,1,xyzzyaaal85&
&),size_grad)-ddot(xyzzyaaao85,dbar(xyzzyaaaf85,1,2,xyzzyaaal85,xyzzya&
&aah85),nemax,dsmat(xyzzyaaai85,xyzzyaaab85,1,xyzzyaaah85,1,xyzzyaaal8&
&5),size_grad)
endif
xyzzyaaas85=xyzzyaaas85-xyzzyaaay85*xyzzyaaba85+xyzzyaaaz85*xyzzyaabb8&
&5
xyzzyaaat85=xyzzyaaat85-xyzzyaaay85*xyzzyaabb85-xyzzyaaaz85*xyzzyaaba8&
&5
elseif(pairing_wf)then
xyzzyaaal85=upd_spin(xyzzyaaac85,xyzzyaaah85)
xyzzyaaam85=upd_spin(xyzzyaaag85,xyzzyaaah85)
if(xyzzyaaam85==xyzzyaaac85)then
xyzzyaaay85=-d2smat(which_d2index(xyzzyaaai85,xyzzyaaaj85),xyzzyaaaf85&
&,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85)
xyzzyaaaz85=-d2smat(which_d2index(xyzzyaaai85,xyzzyaaaj85),xyzzyaaaf85&
&,2,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85)
do xyzzyaaan85=1,xyzzyaaao85
xyzzyaaba85=ddot(xyzzyaaao85,dbar(1,xyzzyaaan85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaai85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)-&
&ddot(xyzzyaaao85,dbar(1,xyzzyaaan85,2,xyzzyaaac85,xyzzyaaah85),1,dsma&
&t(xyzzyaaai85,1,2,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)
xyzzyaabb85=ddot(xyzzyaaao85,dbar(1,xyzzyaaan85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaai85,1,2,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)+&
&ddot(xyzzyaaao85,dbar(1,xyzzyaaan85,2,xyzzyaaac85,xyzzyaaah85),1,dsma&
&t(xyzzyaaai85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)
xyzzyaaay85=xyzzyaaay85+xyzzyaaba85*dsmat(xyzzyaaaj85,xyzzyaaaf85,1,xy&
&zzyaaah85,xyzzyaaan85,xyzzyaaac85)-xyzzyaabb85*dsmat(xyzzyaaaj85,xyzz&
&yaaaf85,2,xyzzyaaah85,xyzzyaaan85,xyzzyaaac85)
xyzzyaaaz85=xyzzyaaaz85+xyzzyaaba85*dsmat(xyzzyaaaj85,xyzzyaaaf85,2,xy&
&zzyaaah85,xyzzyaaan85,xyzzyaaac85)+xyzzyaabb85*dsmat(xyzzyaaaj85,xyzz&
&yaaaf85,1,xyzzyaaah85,xyzzyaaan85,xyzzyaaac85)
enddo
xyzzyaaba85=dbar(xyzzyaaaf85,xyzzyaaab85,1,xyzzyaaac85,xyzzyaaah85)
xyzzyaabb85=dbar(xyzzyaaaf85,xyzzyaaab85,2,xyzzyaaac85,xyzzyaaah85)
xyzzyaaas85=xyzzyaaas85+xyzzyaaay85*xyzzyaaba85-xyzzyaaaz85*xyzzyaabb8&
&5
xyzzyaaat85=xyzzyaaat85+xyzzyaaay85*xyzzyaabb85+xyzzyaaaz85*xyzzyaaba8&
&5
elseif(xyzzyaaal85==xyzzyaaag85)then
xyzzyaaay85=-d2smat(which_d2index(xyzzyaaai85,xyzzyaaaj85),xyzzyaaab85&
&,1,xyzzyaaah85,xyzzyaaaf85,xyzzyaaag85)
xyzzyaaaz85=-d2smat(which_d2index(xyzzyaaai85,xyzzyaaaj85),xyzzyaaab85&
&,2,xyzzyaaah85,xyzzyaaaf85,xyzzyaaag85)
do xyzzyaaan85=1,xyzzyaaap85
xyzzyaaba85=ddot(xyzzyaaap85,dbar(1,xyzzyaaan85,1,xyzzyaaag85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaaj85,1,1,xyzzyaaah85,xyzzyaaaf85,xyzzyaaag85),3)-&
&ddot(xyzzyaaap85,dbar(1,xyzzyaaan85,2,xyzzyaaag85,xyzzyaaah85),1,dsma&
&t(xyzzyaaaj85,1,2,xyzzyaaah85,xyzzyaaaf85,xyzzyaaag85),3)
xyzzyaabb85=ddot(xyzzyaaap85,dbar(1,xyzzyaaan85,1,xyzzyaaag85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaaj85,1,2,xyzzyaaah85,xyzzyaaaf85,xyzzyaaag85),3)+&
&ddot(xyzzyaaap85,dbar(1,xyzzyaaan85,2,xyzzyaaag85,xyzzyaaah85),1,dsma&
&t(xyzzyaaaj85,1,1,xyzzyaaah85,xyzzyaaaf85,xyzzyaaag85),3)
xyzzyaaay85=xyzzyaaay85+xyzzyaaba85*dsmat(xyzzyaaai85,xyzzyaaab85,1,xy&
&zzyaaah85,xyzzyaaan85,xyzzyaaag85)-xyzzyaabb85*dsmat(xyzzyaaai85,xyzz&
&yaaab85,2,xyzzyaaah85,xyzzyaaan85,xyzzyaaag85)
xyzzyaaaz85=xyzzyaaaz85+xyzzyaaba85*dsmat(xyzzyaaai85,xyzzyaaab85,2,xy&
&zzyaaah85,xyzzyaaan85,xyzzyaaag85)+xyzzyaabb85*dsmat(xyzzyaaai85,xyzz&
&yaaab85,1,xyzzyaaah85,xyzzyaaan85,xyzzyaaag85)
enddo
xyzzyaaba85=dbar(xyzzyaaab85,xyzzyaaaf85,1,xyzzyaaag85,xyzzyaaah85)
xyzzyaabb85=dbar(xyzzyaaab85,xyzzyaaaf85,2,xyzzyaaag85,xyzzyaaah85)
xyzzyaaas85=xyzzyaaas85+xyzzyaaay85*xyzzyaaba85-xyzzyaaaz85*xyzzyaabb8&
&5
xyzzyaaat85=xyzzyaaat85+xyzzyaaay85*xyzzyaabb85+xyzzyaaaz85*xyzzyaaba8&
&5
endif
endif
endif
xyzzyaaaq85=xyzzyaaaq85+xyzzyaaaw85(xyzzyaaah85)*xyzzyaaas85-xyzzyaaax&
&85(xyzzyaaah85)*xyzzyaaat85
xyzzyaaar85=xyzzyaaar85+xyzzyaaaw85(xyzzyaaah85)*xyzzyaaat85+xyzzyaaax&
&85(xyzzyaaah85)*xyzzyaaas85
enddo
xyzzyaaay85=xyzzyaaaq85*xyzzyaaau85-xyzzyaaar85*xyzzyaaav85
xyzzyaaaz85=xyzzyaaaq85*xyzzyaaav85+xyzzyaaar85*xyzzyaaau85
harray(xyzzyaaai85,xyzzyaaaj85,1,xyzzyaaaa85,xyzzyaaae85)=xyzzyaaay85
harray(xyzzyaaai85,xyzzyaaaj85,2,xyzzyaaaa85,xyzzyaaae85)=xyzzyaaaz85
if(xyzzyaaaa85/=xyzzyaaae85)then
harray(xyzzyaaaj85,xyzzyaaai85,1,xyzzyaaae85,xyzzyaaaa85)=xyzzyaaay85
harray(xyzzyaaaj85,xyzzyaaai85,2,xyzzyaaae85,xyzzyaaaa85)=xyzzyaaaz85
endif
enddo
enddo
enddo
enddo
else
do xyzzyaaaa85=1,netot
xyzzyaaac85=which_spin(xyzzyaaaa85)
xyzzyaaab85=which_ie(xyzzyaaaa85)
xyzzyaaao85=nele(xyzzyaaac85)
do xyzzyaaad85=1,bf_m2(xyzzyaaaa85)
xyzzyaaae85=bf_rmap2(xyzzyaaad85,xyzzyaaaa85)
if(xyzzyaaae85<xyzzyaaaa85)cycle
xyzzyaaag85=which_spin(xyzzyaaae85)
xyzzyaaaf85=which_ie(xyzzyaaae85)
xyzzyaaap85=nele(xyzzyaaag85)
do xyzzyaaai85=1,dimensionality
do xyzzyaaaj85=1,dimensionality
xyzzyaaak85=which_d2index(xyzzyaaai85,xyzzyaaaj85)
xyzzyaaaq85=0.d0
do xyzzyaaah85=detstart,detstop
if(fpeinfo_pdet(xyzzyaaah85)==xyzzyaaag1)then
cycle
elseif(fpeinfo_pdet(xyzzyaaah85)==xyzzyaaaf1)then
cycle
endif
xyzzyaaas85=0.d0
if(xyzzyaaaa85==xyzzyaaae85)then
if(update_by_column(xyzzyaaac85,xyzzyaaah85))then
xyzzyaaay85=ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,d2smat(xyzzyaaak85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),6)
else
xyzzyaaal85=upd_spin(xyzzyaaac85,xyzzyaaah85)
xyzzyaaay85=ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,1,xyzzyaaal85,xyzzyaaa&
&h85),nemax,d2smat(xyzzyaaak85,xyzzyaaab85,1,xyzzyaaah85,1,xyzzyaaal85&
&),size_sderivs)
endif
xyzzyaaas85=xyzzyaaas85+xyzzyaaay85
else
xyzzyaaas85=xyzzyaaas85+farray_det(xyzzyaaai85,1,xyzzyaaah85,xyzzyaaaa&
&85)*farray_det(xyzzyaaaj85,1,xyzzyaaah85,xyzzyaaae85)
if(xyzzyaaac85==xyzzyaaag85)then
if(update_by_column(xyzzyaaac85,xyzzyaaah85))then
xyzzyaaay85=ddot(xyzzyaaao85,dbar(1,xyzzyaaab85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaaj85,1,1,xyzzyaaah85,xyzzyaaaf85,xyzzyaaac85),3)
xyzzyaaba85=ddot(xyzzyaaao85,dbar(1,xyzzyaaaf85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaai85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)
else
xyzzyaaal85=upd_spin(xyzzyaaac85,xyzzyaaah85)
xyzzyaaay85=-ddot(xyzzyaaao85,dbar(xyzzyaaab85,1,1,xyzzyaaal85,xyzzyaa&
&ah85),nemax,dsmat(xyzzyaaaj85,xyzzyaaaf85,1,xyzzyaaah85,1,xyzzyaaal85&
&),size_grad)
xyzzyaaba85=-ddot(xyzzyaaao85,dbar(xyzzyaaaf85,1,1,xyzzyaaal85,xyzzyaa&
&ah85),nemax,dsmat(xyzzyaaai85,xyzzyaaab85,1,xyzzyaaah85,1,xyzzyaaal85&
&),size_grad)
endif
xyzzyaaas85=xyzzyaaas85-xyzzyaaay85*xyzzyaaba85
elseif(pairing_wf)then
xyzzyaaal85=upd_spin(xyzzyaaac85,xyzzyaaah85)
xyzzyaaam85=upd_spin(xyzzyaaag85,xyzzyaaah85)
if(xyzzyaaam85==xyzzyaaac85)then
xyzzyaaay85=-d2smat(which_d2index(xyzzyaaai85,xyzzyaaaj85),xyzzyaaaf85&
&,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85)
do xyzzyaaan85=1,xyzzyaaao85
xyzzyaaba85=ddot(xyzzyaaao85,dbar(1,xyzzyaaan85,1,xyzzyaaac85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaai85,1,1,xyzzyaaah85,xyzzyaaab85,xyzzyaaac85),3)
xyzzyaaay85=xyzzyaaay85+xyzzyaaba85*dsmat(xyzzyaaaj85,xyzzyaaaf85,1,xy&
&zzyaaah85,xyzzyaaan85,xyzzyaaac85)
enddo
xyzzyaaba85=dbar(xyzzyaaaf85,xyzzyaaab85,1,xyzzyaaac85,xyzzyaaah85)
xyzzyaaas85=xyzzyaaas85+xyzzyaaay85*xyzzyaaba85
elseif(xyzzyaaal85==xyzzyaaag85)then
xyzzyaaay85=-d2smat(which_d2index(xyzzyaaai85,xyzzyaaaj85),xyzzyaaab85&
&,1,xyzzyaaah85,xyzzyaaaf85,xyzzyaaag85)
do xyzzyaaan85=1,xyzzyaaap85
xyzzyaaba85=ddot(xyzzyaaap85,dbar(1,xyzzyaaan85,1,xyzzyaaag85,xyzzyaaa&
&h85),1,dsmat(xyzzyaaaj85,1,1,xyzzyaaah85,xyzzyaaaf85,xyzzyaaag85),3)
xyzzyaaay85=xyzzyaaay85+xyzzyaaba85*dsmat(xyzzyaaai85,xyzzyaaab85,1,xy&
&zzyaaah85,xyzzyaaan85,xyzzyaaag85)
enddo
xyzzyaaba85=dbar(xyzzyaaab85,xyzzyaaaf85,1,xyzzyaaag85,xyzzyaaah85)
xyzzyaaas85=xyzzyaaas85+xyzzyaaay85*xyzzyaaba85
endif
endif
endif
xyzzyaaaq85=xyzzyaaaq85+xyzzyaaaw85(xyzzyaaah85)*xyzzyaaas85
enddo
xyzzyaaay85=xyzzyaaaq85*xyzzyaaau85
harray(xyzzyaaai85,xyzzyaaaj85,1,xyzzyaaaa85,xyzzyaaae85)=xyzzyaaay85
if(xyzzyaaaa85/=xyzzyaaae85)harray(xyzzyaaaj85,xyzzyaaai85,1,xyzzyaaae&
&85,xyzzyaaaa85)=xyzzyaaay85
enddo
enddo
enddo
enddo
endif
call timer('EVAL_HARRAY',.false.)
end subroutine xyzzyaagj1
real(dp) function get_slater_rmax()
implicit none
get_slater_rmax=get_wfdet_rmax()
end function get_slater_rmax
end module slaarnacj

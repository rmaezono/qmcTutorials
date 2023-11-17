module slaarnabe
use casl
use dsp
use slaarnach
use store
use slaarnaag,   only : czero,c_one,third
use format_utils,only : i2s,wout,switch_case
use slaarnabt,   only : lu_decom,lu_decom_cmplx_dz,lu_logdet,lu_logdet&
&_cmplx_dz,lu_solve_n,lu_solve_n_cmplx_dz,ddot,resize_pointer,iswap,sw&
&ap1,dgemv,daxpy,dscal,d_or_z_gemm,d_or_z_gemv,zdotu_cc,zaxpy_cc,zscal&
&_cc,d_sparse_matrix_vector,d_sparse_vector_matrix,d2_sparse_matrix_ve&
&ctor,d2_sparse_vector_matrix,d_or_z_ger,d_or_z_scal
use parallel,    only : am_master
use run_control, only : check_alloc,errstop_master,timer
use slaarnacq, only : wfdet_norb,wfdet
use slaarnabg,    only : dimensionality
use slaarnaad,    only : setup_bf,finish_bf,accept_move_bf,clear_scrat&
&ch_bf,reset_config_bf,get_bf_x,loggrad_bf,loglap_bf,setup_bf_params,f&
&inish_bf_params,get_bf_params,put_bf_params,invalidate_param1_bf
implicit none
private
public enumerate_plot_geminal,query_plot_geminal,get_plot_geminal,fini&
&sh_plot_geminal,read_geminal,update_geminal_casl
public query_geminal_levels,query_geminal_level_details,setup_geminal,&
&finish_geminal,wfn_ratio_geminal,accept_move_geminal,reset_config_gem&
&inal,wfn_logval_geminal,wfn_loggrad_geminal,wfn_loglap_geminal,prefet&
&ch_wfn_geminal,clear_scratch_geminal,add_config_geminal_items,setup_g&
&eminal_params,finish_geminal_params,get_geminal_params,put_geminal_pa&
&rams,clone_scratch_geminal,invalidate_params_geminal,invalidate_param&
&1_geminal,setup_storage_geminal,finish_storage_geminal,load_from_stor&
&age_geminal,save_to_storage_geminal
public gen_config_geminal,delete_config_geminal,copy_config_geminal,co&
&nfig_to_pt_geminal,pt_to_config_geminal,redist_allocations_geminal,re&
&dist_load_geminal,redist_send_geminal,redist_recv_geminal,redist_save&
&_geminal,redist_deallocations_geminal,load_from_pt_geminal,save_to_pt&
&_geminal
public config_wfn_geminal
integer,parameter :: xyzzyaaaa1=1
integer,parameter :: xyzzyaaab1=2
integer,parameter :: xyzzyaaac1=0
integer :: xyzzyaaad1=xyzzyaaac1
logical :: xyzzyaaae1=.true.
real(dp),parameter :: xyzzyaaaf1=-30.d0
integer,parameter :: xyzzyaaag1=0,xyzzyaaah1=1,xyzzyaaai1=2
integer,parameter :: xyzzyaaaj1=100
integer xyzzyaaak1
real(dp),allocatable :: xyzzyaaal1(:,:,:)
real(dp),allocatable :: xyzzyaaam1(:)
logical xyzzyaaan1
integer,allocatable :: xyzzyaaao1(:,:,:),xyzzyaaap1(:,:,:)
integer,allocatable :: xyzzyaaaq1(:,:),xyzzyaaar1(:,:)
integer xyzzyaaas1,xyzzyaaat1
integer,parameter :: xyzzyaaau1=2,xyzzyaaav1=0,xyzzyaaaw1=1,xyzzyaaax1&
&=-1
integer,parameter :: xyzzyaaay1=xyzzyaaav1
integer xyzzyaaaz1,xyzzyaaba1
integer xyzzyaabb1,xyzzyaabc1,xyzzyaabd1,xyzzyaabe1
integer,allocatable :: xyzzyaabf1(:,:,:)
integer,allocatable :: xyzzyaabg1(:)
integer,allocatable :: xyzzyaabh1(:,:)
integer xyzzyaabi1,xyzzyaabj1
integer,pointer :: xyzzyaabk1(:,:,:),xyzzyaabl1(:,:)
integer,pointer :: xyzzyaabm1(:),xyzzyaabn1(:)
logical,allocatable :: xyzzyaabo1(:,:)
type config_wfn_geminal
private
logical dummy_variable
end type config_wfn_geminal
real(dp),allocatable :: xyzzyaabp1(:,:,:,:,:)
logical,allocatable :: xyzzyaabq1(:,:,:)
real(dp),allocatable :: xyzzyaabr1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaabs1(:,:,:)
real(dp),allocatable :: xyzzyaabt1(:,:,:,:,:)
logical,allocatable :: xyzzyaabu1(:,:,:)
real(dp),allocatable :: xyzzyaabv1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaabw1(:,:,:)
real(dp),allocatable :: xyzzyaabx1(:,:,:,:,:)
logical,allocatable :: xyzzyaaby1(:,:,:)
real(dp), allocatable :: xyzzyaabz1(:,:,:,:,:)
logical, allocatable :: xyzzyaaca1(:,:,:)
real(dp), allocatable :: xyzzyaacb1(:,:,:,:,:)
logical, allocatable :: xyzzyaacc1(:,:)
real(dp),allocatable :: xyzzyaacd1(:,:,:,:)
logical,allocatable :: xyzzyaace1(:,:)
real(dp),allocatable :: xyzzyaacf1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaacg1(:,:)
real(dp),allocatable :: xyzzyaach1(:,:,:,:,:)
logical,allocatable :: xyzzyaaci1(:,:)
real(dp), allocatable :: xyzzyaacj1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaack1(:,:)
real(dp),allocatable :: xyzzyaacl1(:,:,:,:,:)
logical,allocatable :: xyzzyaacm1(:,:)
real(dp),allocatable :: xyzzyaacn1(:,:,:,:,:)
logical,allocatable :: xyzzyaaco1(:,:)
real(dp),allocatable :: xyzzyaacp1(:,:,:,:)
logical,allocatable :: xyzzyaacq1(:,:)
real(dp),allocatable :: xyzzyaacr1(:,:,:,:,:)
logical,allocatable :: xyzzyaacs1(:,:)
real(dp),allocatable :: xyzzyaact1(:,:,:,:)
logical,allocatable :: xyzzyaacu1(:,:)
real(dp),allocatable :: xyzzyaacv1(:,:,:,:,:,:)
real(dp),allocatable :: xyzzyaacw1(:,:,:,:,:,:)
real(dp),allocatable :: xyzzyaacx1(:,:,:,:,:,:,:)
logical,allocatable :: xyzzyaacy1(:,:)
real(dp),allocatable :: xyzzyaacz1(:,:,:,:,:)
integer,allocatable :: xyzzyaada1(:,:)
logical,allocatable :: xyzzyaadb1(:,:)
complex(dp),allocatable :: xyzzyaadc1(:,:)
integer,allocatable :: xyzzyaadd1(:,:)
logical, allocatable :: xyzzyaade1(:,:)
complex(dp),allocatable :: xyzzyaadf1(:,:)
logical,allocatable :: xyzzyaadg1(:,:)
complex(dp),allocatable :: xyzzyaadh1(:),xyzzyaadi1(:)
logical,allocatable :: xyzzyaadj1(:)
integer,allocatable :: xyzzyaadk1(:)
real(dp), allocatable :: xyzzyaadl1(:,:,:,:)
logical, allocatable :: xyzzyaadm1(:,:)
real(dp),allocatable :: xyzzyaadn1(:,:,:,:,:,:)
logical,allocatable :: xyzzyaado1(:)
complex(dp),allocatable :: xyzzyaadp1(:,:,:)
logical,allocatable :: xyzzyaadq1(:,:)
complex(dp),allocatable :: xyzzyaadr1(:,:)
logical,allocatable :: xyzzyaads1(:,:)
complex(dp),allocatable :: xyzzyaadt1(:,:,:,:)
logical,allocatable :: xyzzyaadu1(:,:,:)
complex(dp),allocatable :: xyzzyaadv1(:,:,:)
logical,allocatable :: xyzzyaadw1(:,:,:)
real(dp),allocatable :: xyzzyaadx1(:,:,:,:,:,:,:)
logical,allocatable :: xyzzyaady1(:,:)
integer,allocatable :: xyzzyaadz1(:)
logical,allocatable :: xyzzyaaea1(:),xyzzyaaeb1(:,:,:)
real(dp),allocatable :: xyzzyaaec1(:,:),xyzzyaaed1(:,:,:),xyzzyaaee1(:&
&,:,:),xyzzyaaef1(:,:),xyzzyaaeg1(:,:)
complex(dp),allocatable :: xyzzyaaeh1(:,:),xyzzyaaei1(:,:)
real(dp),allocatable :: xyzzyaaej1(:,:,:)
complex(dp),allocatable :: xyzzyaaek1(:)
real(dp),allocatable :: xyzzyaael1(:,:,:,:)
real(dp),allocatable :: xyzzyaaem1(:,:,:,:)
contains
subroutine query_geminal_levels(nlevels)
implicit none
integer,intent(out) :: nlevels
nlevels=1
end subroutine query_geminal_levels
subroutine query_geminal_level_details(level_score,level_name)
implicit none
integer,intent(out) :: level_score(*)
character(80),intent(out) :: level_name(*)
level_score(1)=200
level_name(1)='Geminals'
end subroutine query_geminal_level_details
subroutine setup_geminal
implicit none
integer xyzzyaaaa4
integer xyzzyaaab4
allocate(xyzzyaabp1(wfdet_norb,nemax,real1_complex2,2,nscratch),xyzzya&
&abr1(3,wfdet_norb,nemax,real1_complex2,2,nscratch),xyzzyaabt1(wfdet_n&
&orb,nemax,real1_complex2,2,nscratch),xyzzyaabq1(nemax,2,nscratch),xyz&
&zyaabs1(nemax,2,nscratch),xyzzyaabu1(nemax,2,nscratch),stat=xyzzyaaab&
&4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","orb_val")
xyzzyaabp1=0.d0
xyzzyaabr1=0.d0
xyzzyaabt1=0.d0
xyzzyaabq1=.false.
xyzzyaabs1=.false.
xyzzyaabu1=.false.
allocate(xyzzyaabx1(wfdet_norb,nele(2),real1_complex2,xyzzyaaak1,nscra&
&tch),xyzzyaabz1(wfdet_norb,nele(1),real1_complex2,xyzzyaaak1,nscratch&
&),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gphi")
xyzzyaabx1=0.d0
xyzzyaabz1=0.d0
allocate(xyzzyaaby1(nele(2),xyzzyaaak1,nscratch),xyzzyaaca1(nele(1),xy&
&zzyaaak1,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gphi_valid")
xyzzyaaby1=.false.
xyzzyaaca1=.false.
allocate(xyzzyaadc1(xyzzyaaak1,nscratch),xyzzyaadd1(xyzzyaaak1,nscratc&
&h),xyzzyaade1(xyzzyaaak1,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","loggem")
xyzzyaadc1=czero
xyzzyaadd1=xyzzyaaag1
xyzzyaade1=.false.
allocate(xyzzyaadh1(nscratch),xyzzyaadi1(nscratch),xyzzyaadk1(nscratch&
&),xyzzyaadj1(nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","logwfn_scr")
xyzzyaadh1=czero
xyzzyaadi1=czero
xyzzyaadk1=xyzzyaaag1
xyzzyaadj1=.false.
allocate(xyzzyaadt1(3,xyzzyaaak1,netot,nscratch),xyzzyaadu1(xyzzyaaak1&
&,netot,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","figem")
xyzzyaadt1=0.d0
xyzzyaadu1=.false.
allocate(xyzzyaadp1(3,netot,nscratch),xyzzyaadq1(netot,nscratch),stat=&
&xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","loggrad_wfn_scr")
xyzzyaadp1=czero
xyzzyaadq1=.false.
if(.not.use_backflow)then
allocate(xyzzyaadv1(xyzzyaaak1,netot,nscratch),xyzzyaadw1(xyzzyaaak1,n&
&etot,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","lapgem")
xyzzyaadv1=0.d0
xyzzyaadw1=.false.
endif
allocate(xyzzyaadr1(netot,nscratch),xyzzyaads1(netot,nscratch),stat=xy&
&zzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","loglap_wfn_scr")
xyzzyaadr1=czero
xyzzyaads1=.false.
allocate(xyzzyaadf1(xyzzyaaak1,nscratch),xyzzyaadg1(xyzzyaaak1,nscratc&
&h),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","det_ratio_scr")
xyzzyaadf1=c_one
xyzzyaadg1=.false.
allocate(xyzzyaacb1(nele(1),nele(2),real1_complex2,xyzzyaaak1,nscratch&
&),xyzzyaacc1(xyzzyaaak1,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_matrix")
xyzzyaacb1=0.d0
xyzzyaacc1=.false.
allocate(xyzzyaacd1(nemax,real1_complex2,xyzzyaaak1,nscratch),xyzzyaac&
&e1(xyzzyaaak1,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_matrix1_chscr")
xyzzyaacd1=0.d0
xyzzyaace1=.false.
allocate(xyzzyaacf1(3,nele(1),nele(2),real1_complex2,xyzzyaaak1,nscrat&
&ch),xyzzyaacj1(3,nele(1),nele(2),real1_complex2,xyzzyaaak1,nscratch),&
&stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*_grad")
xyzzyaacf1=0.d0
xyzzyaacj1=0.d0
allocate(xyzzyaacg1(xyzzyaaak1,nscratch),xyzzyaack1(xyzzyaaak1,nscratc&
&h),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*grad_all_valid")
xyzzyaacg1=.false.
xyzzyaack1=.false.
allocate(xyzzyaach1(3,nemax,real1_complex2,xyzzyaaak1,nscratch),xyzzya&
&acl1(3,nemax,real1_complex2,xyzzyaaak1,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*grad_pair1_chscr")
xyzzyaach1=0.d0
xyzzyaacl1=0.d0
allocate(xyzzyaaci1(xyzzyaaak1,nscratch),xyzzyaacm1(xyzzyaaak1,nscratc&
&h),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*grad_pair1_chscr_val&
&id")
xyzzyaaci1=.false.
xyzzyaacm1=.false.
allocate(xyzzyaacn1(nele(1),nele(2),real1_complex2,xyzzyaaak1,nscratch&
&),xyzzyaacr1(nele(1),nele(2),real1_complex2,xyzzyaaak1,nscratch),stat&
&=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*_lap")
xyzzyaacn1=0.d0
xyzzyaacr1=0.d0
allocate(xyzzyaaco1(xyzzyaaak1,nscratch),xyzzyaacs1(xyzzyaaak1,nscratc&
&h),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*lap_all_valid")
xyzzyaaco1=.false.
xyzzyaacs1=.false.
allocate(xyzzyaacp1(nemax,real1_complex2,xyzzyaaak1,nscratch),xyzzyaac&
&t1(nemax,real1_complex2,xyzzyaaak1,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*lap_pair1_chscr")
xyzzyaacp1=0.d0
xyzzyaact1=0.d0
allocate(xyzzyaacq1(xyzzyaaak1,nscratch),xyzzyaacu1(xyzzyaaak1,nscratc&
&h),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*lap_pair1_chscr_vali&
&d")
xyzzyaacq1=.false.
xyzzyaacu1=.false.
allocate(xyzzyaacz1(nele(1),nele(2),real1_complex2,xyzzyaaak1,nscratch&
&),xyzzyaada1(xyzzyaaak1,nscratch),xyzzyaadb1(xyzzyaaak1,nscratch),sta&
&t=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_inv_matrix")
xyzzyaacz1=0.d0
xyzzyaada1=0
xyzzyaadb1=.false.
if(use_backflow)then
allocate(xyzzyaabv1(6,wfdet_norb,nemax,real1_complex2,nspin,nscratch),&
&xyzzyaabw1(nemax,real1_complex2,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","orb_sderiv")
xyzzyaabv1=0.d0
xyzzyaabw1=.false.
allocate(xyzzyaacv1(6,nemax,nemax,real1_complex2,xyzzyaaak1,nscratch),&
&xyzzyaacw1(6,nemax,nemax,real1_complex2,xyzzyaaak1,nscratch),xyzzyaac&
&x1(3,3,nemax,nemax,real1_complex2,xyzzyaaak1,nscratch),xyzzyaacy1(xyz&
&zyaaak1,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","gem_*sderiv")
xyzzyaacv1=0.d0
xyzzyaacw1=0.d0
xyzzyaacx1=0.d0
xyzzyaacy1=.false.
allocate(xyzzyaadx1(3,3,real1_complex2,netot,netot,xyzzyaaak1,nscratch&
&),xyzzyaady1(xyzzyaaak1,nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","hgem_scr")
xyzzyaadx1=0.d0
xyzzyaady1=.false.
allocate(xyzzyaadl1(3,real1_complex2,netot,nscratch),xyzzyaadm1(netot,&
&nscratch),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,"SETUP_GEMINAL","Farray")
xyzzyaadl1=0.d0
xyzzyaadm1=.false.
allocate(xyzzyaadn1(3,3,real1_complex2,netot,netot,nscratch),xyzzyaado&
&1(nscratch),stat=xyzzyaaab4)
xyzzyaadn1=0.d0
xyzzyaado1=.false.
endif
if(use_backflow)call setup_bf
do xyzzyaaaa4=1,nscratch
call clear_scratch_geminal(xyzzyaaaa4)
enddo
call xyzzyaaen1
end subroutine setup_geminal
subroutine xyzzyaaen1
implicit none
integer xyzzyaaaa5
allocate(xyzzyaadz1(nemax),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','1a')
xyzzyaadz1=0
allocate(xyzzyaaec1(4,netot),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','1b')
xyzzyaaec1=0.d0
allocate(xyzzyaaea1(wfdet_norb),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','1c')
xyzzyaaea1=.true.
allocate(xyzzyaaed1(nemax,nemax,real1_complex2),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','2')
xyzzyaaed1=0.d0
allocate(xyzzyaaee1(3,nemax,real1_complex2),xyzzyaaef1(nemax,real1_com&
&plex2),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','2b')
xyzzyaaee1=0.d0
xyzzyaaef1=0.d0
allocate(xyzzyaaeg1(3,real1_complex2),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','3')
xyzzyaaeg1=0.d0
allocate(xyzzyaaeh1(nemax,nemax),xyzzyaaei1(nemax,nemax),stat=xyzzyaaa&
&a5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','4')
xyzzyaaeh1=0.d0
allocate(xyzzyaaej1(wfdet_norb,nemax,real1_complex2),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','5')
xyzzyaaej1=0.d0
if(use_backflow)then
allocate(xyzzyaaek1(3),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','Farray')
xyzzyaaek1=czero
allocate(xyzzyaael1(3,nemax,nemax,real1_complex2),xyzzyaaem1(3,nemax,n&
&emax,real1_complex2),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'ALLOCATE_TEMPS_GEMINAL','6')
xyzzyaael1=0.d0
xyzzyaaem1=0.d0
endif
end subroutine xyzzyaaen1
subroutine finish_geminal
implicit none
deallocate(xyzzyaabp1,xyzzyaabr1,xyzzyaabt1)
deallocate(xyzzyaabq1,xyzzyaabs1,xyzzyaabu1)
deallocate(xyzzyaabx1,xyzzyaabz1)
deallocate(xyzzyaaby1,xyzzyaaca1)
deallocate(xyzzyaadc1,xyzzyaade1,xyzzyaadd1)
deallocate(xyzzyaadh1,xyzzyaadi1,xyzzyaadj1,xyzzyaadk1)
deallocate(xyzzyaadp1,xyzzyaadq1)
deallocate(xyzzyaadr1,xyzzyaads1)
deallocate(xyzzyaadt1,xyzzyaadu1)
if(.not.use_backflow)deallocate(xyzzyaadv1,xyzzyaadw1)
deallocate(xyzzyaadf1,xyzzyaadg1)
deallocate(xyzzyaacb1,xyzzyaacc1)
deallocate(xyzzyaacd1,xyzzyaace1)
deallocate(xyzzyaacf1,xyzzyaacj1)
deallocate(xyzzyaacg1,xyzzyaack1)
deallocate(xyzzyaach1,xyzzyaacl1)
deallocate(xyzzyaaci1,xyzzyaacm1)
deallocate(xyzzyaacn1,xyzzyaacr1)
deallocate(xyzzyaaco1,xyzzyaacs1)
deallocate(xyzzyaacp1,xyzzyaact1)
deallocate(xyzzyaacq1,xyzzyaacu1)
deallocate(xyzzyaacz1,xyzzyaada1,xyzzyaadb1)
if(use_backflow)then
deallocate(xyzzyaabv1,xyzzyaabw1)
deallocate(xyzzyaacy1)
deallocate(xyzzyaacv1,xyzzyaacw1)
deallocate(xyzzyaacx1)
deallocate(xyzzyaadl1,xyzzyaadm1)
deallocate(xyzzyaadx1,xyzzyaady1)
endif
if(use_backflow)call finish_bf
call xyzzyaaeo1
end subroutine finish_geminal
subroutine xyzzyaaeo1
implicit none
deallocate(xyzzyaadz1,xyzzyaaec1,xyzzyaaea1,xyzzyaaed1,xyzzyaaee1,xyzz&
&yaaef1,xyzzyaaeg1,xyzzyaaeh1,xyzzyaaei1,xyzzyaaej1)
if(use_backflow)then
deallocate(xyzzyaaek1)
deallocate(xyzzyaael1,xyzzyaaem1)
endif
end subroutine xyzzyaaeo1
subroutine wfn_ratio_geminal(is,js,ilevel,ratio,fd,sd,isnan,isinf)
implicit none
integer,intent(in) :: is,js,ilevel
complex(dp),intent(out) :: ratio
logical,intent(in) :: fd,sd
logical,intent(out) :: isnan,isinf
complex(dp) xyzzyaaaa8
call timer("WFN_RATIO_GEMINAL",.true.)
isnan=.false.
isinf=.false.
call xyzzyaaep1(is,fd,fd,sd,sd)
call xyzzyaaep1(js,fd,fd,sd,sd)
if(xyzzyaadk1(is)==xyzzyaaag1)then
if(xyzzyaadk1(js)==xyzzyaaag1)then
xyzzyaaaa8=(xyzzyaadh1(js)+xyzzyaadi1(js))-(xyzzyaadh1(is)+xyzzyaadi1(&
&is))
if(dble(xyzzyaaaa8)<xyzzyaaaf1)then
ratio=czero
else
ratio=exp(xyzzyaaaa8)
endif
else
ratio=czero
endif
else
ratio=czero
if(xyzzyaadk1(js)==xyzzyaaag1)then
isinf=.true.
else
isnan=.true.
endif
endif
call timer("WFN_RATIO_GEMINAL",.false.)
end subroutine wfn_ratio_geminal
subroutine accept_move_geminal(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,xyzzyaaae9,xyzzyaa&
&af9,xyzzyaaag9
logical xyzzyaaah9
call timer("ACCEPT_MOVE_GEMINAL",.true.)
xyzzyaaaf9=buffer_move1_from(js)
xyzzyaaah9=(xyzzyaaaf9==is.and..not.use_backflow)
if(xyzzyaaah9)then
xyzzyaaac9=buffer_move1_from_ii(js)
xyzzyaaae9=which_spin(xyzzyaaac9)
xyzzyaaad9=which_ie(xyzzyaaac9)
endif
if(xyzzyaaah9)then
if(xyzzyaabq1(xyzzyaaad9,xyzzyaaae9,js))then
xyzzyaabp1(:,xyzzyaaad9,:,xyzzyaaae9,is)=xyzzyaabp1(:,xyzzyaaad9,:,xyz&
&zyaaae9,js)
xyzzyaabq1(xyzzyaaad9,xyzzyaaae9,is)=.true.
else
xyzzyaabq1(xyzzyaaad9,xyzzyaaae9,is)=.false.
endif
if(xyzzyaabs1(xyzzyaaad9,xyzzyaaae9,js))then
xyzzyaabr1(:,:,xyzzyaaad9,:,xyzzyaaae9,is)=xyzzyaabr1(:,:,xyzzyaaad9,:&
&,xyzzyaaae9,js)
xyzzyaabs1(xyzzyaaad9,xyzzyaaae9,is)=.true.
else
xyzzyaabs1(xyzzyaaad9,xyzzyaaae9,is)=.false.
endif
if(xyzzyaabu1(xyzzyaaad9,xyzzyaaae9,js))then
xyzzyaabt1(:,xyzzyaaad9,:,xyzzyaaae9,is)=xyzzyaabt1(:,xyzzyaaad9,:,xyz&
&zyaaae9,js)
xyzzyaabu1(xyzzyaaad9,xyzzyaaae9,is)=.true.
else
xyzzyaabu1(xyzzyaaad9,xyzzyaaae9,is)=.false.
endif
else
xyzzyaabq1(:,:,is)=.false.
xyzzyaabs1(:,:,is)=.false.
xyzzyaabu1(:,:,is)=.false.
forall(xyzzyaaaa9=1:nemax,xyzzyaaab9=1:2,xyzzyaabq1(xyzzyaaaa9,xyzzyaa&
&ab9,js))
xyzzyaabp1(:,xyzzyaaaa9,:,xyzzyaaab9,is)=xyzzyaabp1(:,xyzzyaaaa9,:,xyz&
&zyaaab9,js)
xyzzyaabq1(xyzzyaaaa9,xyzzyaaab9,is)=.true.
endforall
forall(xyzzyaaaa9=1:nemax,xyzzyaaab9=1:2,xyzzyaabs1(xyzzyaaaa9,xyzzyaa&
&ab9,js))
xyzzyaabr1(:,:,xyzzyaaaa9,:,xyzzyaaab9,is)=xyzzyaabr1(:,:,xyzzyaaaa9,:&
&,xyzzyaaab9,js)
xyzzyaabs1(xyzzyaaaa9,xyzzyaaab9,is)=.true.
endforall
forall(xyzzyaaaa9=1:nemax,xyzzyaaab9=1:2,xyzzyaabu1(xyzzyaaaa9,xyzzyaa&
&ab9,js))
xyzzyaabt1(:,xyzzyaaaa9,:,xyzzyaaab9,is)=xyzzyaabt1(:,xyzzyaaaa9,:,xyz&
&zyaaab9,js)
xyzzyaabu1(xyzzyaaaa9,xyzzyaaab9,is)=.true.
endforall
endif
if(xyzzyaaah9)then
select case(xyzzyaaae9)
case(1)
xyzzyaaca1(xyzzyaaad9,:,is)=.false.
forall(xyzzyaaag9=1:xyzzyaaak1,xyzzyaaca1(xyzzyaaad9,xyzzyaaag9,js))
xyzzyaabz1(:,xyzzyaaad9,:,xyzzyaaag9,is)=xyzzyaabz1(:,xyzzyaaad9,:,xyz&
&zyaaag9,js)
xyzzyaaca1(xyzzyaaad9,xyzzyaaag9,is)=.true.
endforall
case(2)
xyzzyaaby1(xyzzyaaad9,:,is)=.false.
forall(xyzzyaaag9=1:xyzzyaaak1,xyzzyaaby1(xyzzyaaad9,xyzzyaaag9,js))
xyzzyaabx1(:,xyzzyaaad9,:,xyzzyaaag9,is)=xyzzyaabx1(:,xyzzyaaad9,:,xyz&
&zyaaag9,js)
xyzzyaaby1(xyzzyaaad9,xyzzyaaag9,is)=.true.
endforall
endselect
else
xyzzyaaby1(:,:,is)=.false.
forall(xyzzyaaaa9=1:nemax,xyzzyaaag9=1:xyzzyaaak1,xyzzyaaby1(xyzzyaaaa&
&9,xyzzyaaag9,js))
xyzzyaabx1(:,xyzzyaaaa9,:,xyzzyaaag9,is)=xyzzyaabx1(:,xyzzyaaaa9,:,xyz&
&zyaaag9,js)
xyzzyaaby1(xyzzyaaaa9,xyzzyaaag9,is)=.true.
endforall
xyzzyaaca1(:,:,is)=.false.
forall(xyzzyaaaa9=1:nemax,xyzzyaaag9=1:xyzzyaaak1,xyzzyaaca1(xyzzyaaaa&
&9,xyzzyaaag9,js))
xyzzyaabz1(:,xyzzyaaaa9,:,xyzzyaaag9,is)=xyzzyaabz1(:,xyzzyaaaa9,:,xyz&
&zyaaag9,js)
xyzzyaaca1(xyzzyaaaa9,xyzzyaaag9,is)=.true.
endforall
endif
do xyzzyaaag9=1,xyzzyaaak1
if(xyzzyaadb1(xyzzyaaag9,js))then
xyzzyaacz1(:,:,:,xyzzyaaag9,is)=xyzzyaacz1(:,:,:,xyzzyaaag9,js)
xyzzyaada1(xyzzyaaag9,is)=xyzzyaada1(xyzzyaaag9,js)
xyzzyaadb1(xyzzyaaag9,is)=.true.
elseif(xyzzyaaah9.and.xyzzyaadb1(xyzzyaaag9,is).and.xyzzyaadg1(xyzzyaa&
&ag9,js).and.xyzzyaada1(xyzzyaaag9,is)<xyzzyaaaj1)then
if(xyzzyaadd1(xyzzyaaag9,is)/=xyzzyaaag1)then
call xyzzyaafc1(xyzzyaaag9,js,is,.true.,.false.,.false.,.false.,.false&
&.)
xyzzyaaed1=xyzzyaacb1(:,:,:,xyzzyaaag9,js)
call xyzzyaagl1(nemax,xyzzyaadc1(xyzzyaaag9,js),xyzzyaaed1,xyzzyaacz1(&
&:,:,:,xyzzyaaag9,is),xyzzyaadd1(xyzzyaaag9,is))
else
call xyzzyaagn1(xyzzyaaad9,xyzzyaaae9,xyzzyaacd1(:,:,xyzzyaaag9,js),xy&
&zzyaacz1(:,:,:,xyzzyaaag9,is),xyzzyaadf1(xyzzyaaag9,js))
endif
xyzzyaadb1(xyzzyaaag9,is)=.true.
xyzzyaada1(xyzzyaaag9,is)=xyzzyaada1(xyzzyaaag9,is)+1
else
xyzzyaadb1(xyzzyaaag9,is)=.false.
endif
enddo
do xyzzyaaag9=1,xyzzyaaak1
if(xyzzyaaah9.and.xyzzyaacc1(xyzzyaaag9,is).and.xyzzyaace1(xyzzyaaag9,&
&js))then
select case(xyzzyaaae9)
case(1)
xyzzyaacb1(xyzzyaaad9,:,:,xyzzyaaag9,is)=xyzzyaacd1(:,:,xyzzyaaag9,js)
case(2)
xyzzyaacb1(:,xyzzyaaad9,:,xyzzyaaag9,is)=xyzzyaacd1(:,:,xyzzyaaag9,js)
endselect
elseif(xyzzyaacc1(xyzzyaaag9,js))then
xyzzyaacb1(:,:,:,xyzzyaaag9,is)=xyzzyaacb1(:,:,:,xyzzyaaag9,js)
xyzzyaacc1(xyzzyaaag9,is)=.true.
else
xyzzyaacc1(xyzzyaaag9,is)=.false.
endif
enddo
do xyzzyaaag9=1,xyzzyaaak1
if(xyzzyaaah9.and.xyzzyaacg1(xyzzyaaag9,is).and.xyzzyaaci1(xyzzyaaag9,&
&js))then
select case(xyzzyaaae9)
case(1)
xyzzyaacf1(:,xyzzyaaad9,:,:,xyzzyaaag9,is)=xyzzyaach1(:,:,:,xyzzyaaag9&
&,js)
case(2)
xyzzyaacf1(:,:,xyzzyaaad9,:,xyzzyaaag9,is)=xyzzyaach1(:,:,:,xyzzyaaag9&
&,js)
endselect
elseif(xyzzyaacg1(xyzzyaaag9,js))then
xyzzyaacf1(:,:,:,:,xyzzyaaag9,is)=xyzzyaacf1(:,:,:,:,xyzzyaaag9,js)
xyzzyaacg1(xyzzyaaag9,is)=.true.
else
xyzzyaacg1(xyzzyaaag9,is)=.false.
endif
enddo
do xyzzyaaag9=1,xyzzyaaak1
if(xyzzyaaah9.and.xyzzyaack1(xyzzyaaag9,is).and.xyzzyaacm1(xyzzyaaag9,&
&js))then
select case(xyzzyaaae9)
case(1)
xyzzyaacj1(:,xyzzyaaad9,:,:,xyzzyaaag9,is)=xyzzyaacl1(:,:,:,xyzzyaaag9&
&,js)
case(2)
xyzzyaacj1(:,:,xyzzyaaad9,:,xyzzyaaag9,is)=xyzzyaacl1(:,:,:,xyzzyaaag9&
&,js)
endselect
elseif(xyzzyaack1(xyzzyaaag9,js))then
xyzzyaacj1(:,:,:,:,xyzzyaaag9,is)=xyzzyaacj1(:,:,:,:,xyzzyaaag9,js)
xyzzyaack1(xyzzyaaag9,is)=.true.
else
xyzzyaack1(xyzzyaaag9,is)=.false.
endif
enddo
do xyzzyaaag9=1,xyzzyaaak1
if(xyzzyaaah9.and.xyzzyaaco1(xyzzyaaag9,is).and.xyzzyaacq1(xyzzyaaag9,&
&js))then
select case(xyzzyaaae9)
case(1)
xyzzyaacn1(xyzzyaaad9,:,:,xyzzyaaag9,is)=xyzzyaacp1(:,:,xyzzyaaag9,js)
case(2)
xyzzyaacn1(:,xyzzyaaad9,:,xyzzyaaag9,is)=xyzzyaacp1(:,:,xyzzyaaag9,js)
endselect
elseif(xyzzyaaco1(xyzzyaaag9,js))then
xyzzyaacn1(:,:,:,xyzzyaaag9,is)=xyzzyaacn1(:,:,:,xyzzyaaag9,js)
xyzzyaaco1(xyzzyaaag9,is)=.true.
else
xyzzyaaco1(xyzzyaaag9,is)=.false.
endif
enddo
do xyzzyaaag9=1,xyzzyaaak1
if(xyzzyaaah9.and.xyzzyaacs1(xyzzyaaag9,is).and.xyzzyaacu1(xyzzyaaag9,&
&js))then
select case(xyzzyaaae9)
case(1)
xyzzyaacr1(xyzzyaaad9,:,:,xyzzyaaag9,is)=xyzzyaact1(:,:,xyzzyaaag9,js)
case(2)
xyzzyaacr1(:,xyzzyaaad9,:,xyzzyaaag9,is)=xyzzyaact1(:,:,xyzzyaaag9,js)
endselect
elseif(xyzzyaacs1(xyzzyaaag9,js))then
xyzzyaacr1(:,:,:,xyzzyaaag9,is)=xyzzyaacr1(:,:,:,xyzzyaaag9,js)
xyzzyaacs1(xyzzyaaag9,is)=.true.
else
xyzzyaacs1(xyzzyaaag9,is)=.false.
endif
enddo
where(xyzzyaade1(:,js))
xyzzyaadc1(:,is)=xyzzyaadc1(:,js)
xyzzyaade1(:,is)=.true.
elsewhere
xyzzyaade1(:,is)=.false.
endwhere
xyzzyaadu1(:,:,is)=xyzzyaadu1(:,:,js)
forall(xyzzyaaag9=1:xyzzyaaak1,xyzzyaaad9=1:netot,xyzzyaadu1(xyzzyaaag&
&9,xyzzyaaad9,js))
xyzzyaadt1(:,xyzzyaaag9,xyzzyaaad9,is)=xyzzyaadt1(:,xyzzyaaag9,xyzzyaa&
&ad9,js)
endforall
if(.not.use_backflow)then
xyzzyaadw1(:,:,is)=xyzzyaadw1(:,:,js)
where(xyzzyaadw1(:,:,js))
xyzzyaadv1(:,:,is)=xyzzyaadv1(:,:,js)
endwhere
endif
if(xyzzyaadj1(js))then
xyzzyaadh1(is)=xyzzyaadh1(js)
xyzzyaadi1(is)=xyzzyaadi1(js)
xyzzyaadk1(is)=xyzzyaadk1(js)
xyzzyaadj1(is)=.true.
else
xyzzyaadj1(is)=.false.
endif
xyzzyaadq1(:,is)=xyzzyaadq1(:,js)
forall(xyzzyaaaa9=1:netot,xyzzyaadq1(xyzzyaaaa9,js))
xyzzyaadp1(:,xyzzyaaaa9,is)=xyzzyaadp1(:,xyzzyaaaa9,js)
endforall
xyzzyaads1(:,is)=xyzzyaads1(:,js)
forall(xyzzyaaaa9=1:netot,xyzzyaadq1(xyzzyaaaa9,js))
xyzzyaadr1(xyzzyaaaa9,is)=xyzzyaadr1(xyzzyaaaa9,js)
endforall
where(xyzzyaade1(:,js))
xyzzyaadd1(:,is)=xyzzyaadd1(:,js)
endwhere
if(use_backflow)then
xyzzyaabw1(:,:,is)=.false.
forall(xyzzyaaaa9=1:nemax,xyzzyaaab9=1:2,xyzzyaabw1(xyzzyaaaa9,xyzzyaa&
&ab9,js))
xyzzyaabv1(:,:,xyzzyaaaa9,:,xyzzyaaab9,is)=xyzzyaabv1(:,:,xyzzyaaaa9,:&
&,xyzzyaaab9,js)
xyzzyaabw1(xyzzyaaaa9,xyzzyaaab9,is)=.true.
endforall
xyzzyaacy1(:,is)=.false.
forall(xyzzyaaaa9=1:xyzzyaaak1,xyzzyaacy1(xyzzyaaaa9,js))
xyzzyaacv1(:,:,:,:,xyzzyaaaa9,is)=xyzzyaacv1(:,:,:,:,xyzzyaaaa9,js)
xyzzyaacw1(:,:,:,:,xyzzyaaaa9,is)=xyzzyaacw1(:,:,:,:,xyzzyaaaa9,js)
xyzzyaacx1(:,:,:,:,:,xyzzyaaaa9,is)=xyzzyaacx1(:,:,:,:,:,xyzzyaaaa9,js&
&)
xyzzyaacy1(xyzzyaaaa9,is)=.true.
endforall
xyzzyaady1(:,is)=.false.
forall(xyzzyaaaa9=1:xyzzyaaak1,xyzzyaady1(xyzzyaaaa9,js))
xyzzyaadx1(:,:,:,:,:,xyzzyaaaa9,is)=xyzzyaadx1(:,:,:,:,:,xyzzyaaaa9,js&
&)
xyzzyaady1(xyzzyaaaa9,is)=.true.
endforall
xyzzyaadm1(:,is)=.false.
forall(xyzzyaaaa9=1:netot,xyzzyaadm1(xyzzyaaaa9,js))
xyzzyaadl1(1:3,1:real1_complex2,xyzzyaaaa9,is)=xyzzyaadl1(1:3,1:real1_&
&complex2,xyzzyaaaa9,js)
xyzzyaadm1(xyzzyaaaa9,is)=.true.
endforall
if(xyzzyaado1(js))then
xyzzyaadn1(:,:,:,:,:,is)=xyzzyaadn1(:,:,:,:,:,js)
xyzzyaado1(is)=.true.
else
xyzzyaado1=.false.
endif
endif
if(use_backflow)call accept_move_bf(is,js)
call timer("ACCEPT_MOVE_GEMINAL",.false.)
end subroutine accept_move_geminal
subroutine reset_config_geminal(is,js)
implicit none
integer,intent(in) :: is,js
call clear_scratch_geminal(js)
if(use_backflow)call reset_config_bf(is,js)
end subroutine reset_config_geminal
subroutine wfn_logval_geminal(is,logwfn,iszero)
implicit none
integer,intent(in) :: is
complex(dp),intent(out) :: logwfn
logical,intent(out) :: iszero
call timer("WFN_LOGVAL_GEMINAL",.true.)
call xyzzyaaep1(is,.false.,.false.,.false.,.false.)
logwfn=xyzzyaadh1(is)+xyzzyaadi1(is)
iszero=(xyzzyaadk1(is)/=xyzzyaaag1)
call timer("WFN_LOGVAL_GEMINAL",.false.)
end subroutine wfn_logval_geminal
subroutine wfn_loggrad_geminal(ii,is,ilevel,val,sd,loggrad,isnan,isinf&
&)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: loggrad(3)
logical,intent(in) :: val,sd
logical,intent(out) :: isnan,isinf
integer xyzzyaaaa12
call timer("WFN_LOGGRAD_GEMINAL",.true.)
isnan=.false.
isinf=.false.
if(.not.xyzzyaadq1(ii,is))then
if(.not.use_backflow)then
call xyzzyaaeq1(ii,is,sd,xyzzyaadp1(1:3,ii,is),xyzzyaadq1(ii,is))
xyzzyaadq1(ii,is)=.true.
else
do xyzzyaaaa12=1,netot
if(.not.xyzzyaadm1(xyzzyaaaa12,is))then
call xyzzyaaeq1(xyzzyaaaa12,is,sd,xyzzyaaek1,xyzzyaadm1(xyzzyaaaa12,is&
&))
xyzzyaadl1(:,re,xyzzyaaaa12,is)=real(xyzzyaaek1)
xyzzyaadl1(:,im,xyzzyaaaa12,is)=aimag(xyzzyaaek1)
xyzzyaadm1(xyzzyaaaa12,is)=.true.
endif
enddo
call loggrad_bf(ii,is,sd,xyzzyaadl1(1,1,1,is),xyzzyaadp1(1:3,ii,is))
xyzzyaadq1(ii,is)=.true.
endif
endif
loggrad=xyzzyaadp1(:,ii,is)
if(xyzzyaadk1(is)/=xyzzyaaag1)then
if(all(xyzzyaadd1(:,is)==xyzzyaaag1))then
isnan=.true.
else
isinf=.true.
endif
endif
call timer("WFN_LOGGRAD_GEMINAL",.false.)
end subroutine  wfn_loggrad_geminal
subroutine wfn_loglap_geminal(ii,is,ilevel,val,fd,loglap,isnan,isinf)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: loglap
logical,intent(in) :: val,fd
logical,intent(out) :: isnan,isinf
complex(dp) xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13(3),xyzzyaaad13
call timer("WFN_LOGLAP_GEMINAL",.true.)
isnan=.false.
isinf=.false.
if(.not.xyzzyaads1(ii,is))then
call wfn_loggrad_geminal(ii,is,ilevel,.false.,.true.,xyzzyaaac13,isnan&
&,isinf)
if(.not.use_backflow)then
call xyzzyaaet1(ii,is)
if(any(xyzzyaadd1(:,is)==xyzzyaaah1))then
loglap=czero
isnan=.true.
return
endif
xyzzyaaaa13=sum(xyzzyaaam1*xyzzyaadv1(:,ii,is),xyzzyaadd1(:,is)==xyzzy&
&aaag1)
xyzzyaaab13=exp(xyzzyaadi1(is)+xyzzyaadh1(is))
xyzzyaaad13=sum(xyzzyaaac13**2)
xyzzyaadr1(ii,is)=xyzzyaaaa13/xyzzyaaab13-xyzzyaaad13
else
call xyzzyaaey1(is)
if(xyzzyaadk1(is)==xyzzyaaag1)then
call loglap_bf(ii,is,xyzzyaadl1(1,1,1,is),xyzzyaadn1(1,1,1,1,1,is),xyz&
&zyaaac13,xyzzyaadr1(ii,is))
endif
endif
xyzzyaads1(ii,is)=.true.
endif
loglap=xyzzyaadr1(ii,is)
if(xyzzyaadk1(is)/=xyzzyaaag1)then
if(all(xyzzyaadd1(:,is)==xyzzyaaag1))then
isnan=.true.
else
isinf=.true.
endif
endif
call timer("WFN_LOGLAP_GEMINAL",.false.)
end subroutine wfn_loglap_geminal
subroutine prefetch_wfn_geminal(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
end subroutine prefetch_wfn_geminal
subroutine xyzzyaaep1(is,lgrad,rgrad,llap,rlap)
implicit none
integer,intent(in) :: is
logical,intent(in) :: lgrad,rgrad,llap,rlap
integer xyzzyaaaa15
complex(dp) xyzzyaaab15
if(xyzzyaadj1(is))return
call timer("GET_LOGWFN",.true.)
call xyzzyaaer1(is,lgrad,rgrad,llap,rlap)
xyzzyaaaa15=maxloc(real(xyzzyaadc1(:,is),dp),1,xyzzyaadd1(:,is)==xyzzy&
&aaag1)
xyzzyaadi1(is)=xyzzyaadc1(xyzzyaaaa15,is)
xyzzyaaab15=sum(xyzzyaaam1*exp(xyzzyaadc1(:,is)-xyzzyaadi1(is)),xyzzya&
&add1(:,is)==xyzzyaaag1)
if(xyzzyaaab15==czero)then
xyzzyaadk1(is)=xyzzyaaah1
else
xyzzyaadh1(is)=log(xyzzyaaab15)
endif
xyzzyaadj1(is)=.true.
call timer("GET_LOGWFN",.false.)
end subroutine xyzzyaaep1
subroutine xyzzyaaeq1(ii,is,sd,dest,dest_valid)
implicit none
integer,intent(in) :: ii,is
logical,intent(in) :: sd
complex(dp),intent(inout) :: dest(3)
logical,intent(inout) :: dest_valid
integer xyzzyaaaa16
if(dest_valid)return
call timer("FILL_LOGGRAD_WFN_OR_FARRAY",.true.)
call xyzzyaaep1(is,.true.,.true.,sd,sd)
if(xyzzyaadk1(is)==xyzzyaaag1)then
call xyzzyaaes1(ii,is,sd)
dest=czero
do xyzzyaaaa16=1,xyzzyaaak1
if(xyzzyaadd1(xyzzyaaaa16,is)/=xyzzyaaag1)cycle
dest=dest+exp(xyzzyaadc1(xyzzyaaaa16,is))*xyzzyaaam1(xyzzyaaaa16)*xyzz&
&yaadt1(:,xyzzyaaaa16,ii,is)
enddo
dest=dest/exp(xyzzyaadi1(is)+xyzzyaadh1(is))
endif
dest_valid=.true.
call timer("FILL_LOGGRAD_WFN_OR_FARRAY",.false.)
end subroutine xyzzyaaeq1
recursive subroutine xyzzyaaer1(xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xy&
&zzyaaad17,xyzzyaaae17)
implicit none
integer,intent(in) :: xyzzyaaaa17
logical,intent(in) :: xyzzyaaab17,xyzzyaaac17,xyzzyaaad17,xyzzyaaae17
integer xyzzyaaaf17,xyzzyaaag17
call timer("GET_LOGGEMS",.true.)
xyzzyaaaf17=buffer_move1_from(xyzzyaaaa17)
if(xyzzyaaaf17/=0.and..not.use_backflow)then
do xyzzyaaag17=1,xyzzyaaak1
if(xyzzyaade1(xyzzyaaag17,xyzzyaaaa17))cycle
call xyzzyaafb1(xyzzyaaag17,xyzzyaaaf17)
call xyzzyaaez1(xyzzyaaag17,xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzya&
&aad17,xyzzyaaae17)
enddo
else
do xyzzyaaag17=1,xyzzyaaak1
if(xyzzyaade1(xyzzyaaag17,xyzzyaaaa17))cycle
call xyzzyaafa1(xyzzyaaag17,xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzya&
&aad17,xyzzyaaae17)
enddo
endif
call timer("GET_LOGGEMS",.false.)
end subroutine xyzzyaaer1
subroutine xyzzyaaes1(ii,is,sd)
implicit none
integer,intent(in) :: ii,is
logical,intent(in) :: sd
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xy&
&zzyaaaf18,xyzzyaaag18,xyzzyaaah18
call timer("GET_FIGEM",.true.)
xyzzyaaaa18=which_ie(ii)
xyzzyaaab18=which_spin(ii)
do xyzzyaaag18=1,xyzzyaaak1
if(.not.xyzzyaadu1(xyzzyaaag18,xyzzyaaaa18,is))then
call xyzzyaafb1(xyzzyaaag18,is)
endif
enddo
if(any(xyzzyaadd1(:,is)==xyzzyaaah1))return
xyzzyaaac18=buffer_move1_from(is)
if(xyzzyaaac18/=0)then
xyzzyaaad18=buffer_move1_from_ii(is)
xyzzyaaae18=which_ie(xyzzyaaad18)
xyzzyaaaf18=which_spin(xyzzyaaad18)
xyzzyaaah18=buffer_move1_from(xyzzyaaac18)
endif
if(xyzzyaaac18/=0)then
if(xyzzyaaaf18==xyzzyaaab18.and.xyzzyaaae18==xyzzyaaaa18)then
select case(xyzzyaaab18)
case(1)
do xyzzyaaag18=1,xyzzyaaak1
if(xyzzyaadu1(xyzzyaaag18,ii,is))cycle
call xyzzyaafd1(xyzzyaaag18,is,.false.,.true.,.false.,sd,.false.)
call xyzzyaagr1(xyzzyaach1(:,:,:,xyzzyaaag18,is),xyzzyaadd1(xyzzyaaag1&
&8,is),xyzzyaacz1(:,xyzzyaaaa18,:,xyzzyaaag18,is),xyzzyaadt1(:,xyzzyaa&
&ag18,ii,is))
xyzzyaadu1(xyzzyaaag18,ii,is)=.true.
enddo
case(2)
do xyzzyaaag18=1,xyzzyaaak1
if(xyzzyaadu1(xyzzyaaag18,ii,is))cycle
call xyzzyaafd1(xyzzyaaag18,is,.false.,.false.,.true.,.false.,sd)
call xyzzyaagr1(xyzzyaacl1(:,:,:,xyzzyaaag18,is),xyzzyaadd1(xyzzyaaag1&
&8,is),xyzzyaacz1(xyzzyaaaa18,:,:,xyzzyaaag18,is),xyzzyaadt1(:,xyzzyaa&
&ag18,ii,is))
xyzzyaadu1(xyzzyaaag18,ii,is)=.true.
enddo
endselect
else
select case(xyzzyaaab18)
case(1)
select case(xyzzyaaaf18)
case(1)
do xyzzyaaag18=1,xyzzyaaak1
if(xyzzyaadu1(xyzzyaaag18,ii,is))cycle
call xyzzyaafc1(xyzzyaaag18,xyzzyaaac18,xyzzyaaah18,.false.,.true.,.fa&
&lse.,sd,.false.)
xyzzyaaee1=xyzzyaacf1(:,xyzzyaaaa18,:,1:real1_complex2,xyzzyaaag18,xyz&
&zyaaac18)
call xyzzyaagr1(xyzzyaaee1,xyzzyaadd1(xyzzyaaag18,is),xyzzyaacz1(:,xyz&
&zyaaaa18,:,xyzzyaaag18,is),xyzzyaadt1(:,xyzzyaaag18,ii,is))
xyzzyaadu1(xyzzyaaag18,ii,is)=.true.
enddo
case(2)
do xyzzyaaag18=1,xyzzyaaak1
if(xyzzyaadu1(xyzzyaaag18,ii,is))cycle
call xyzzyaafc1(xyzzyaaag18,xyzzyaaac18,xyzzyaaah18,.false.,.true.,.fa&
&lse.,sd,.false.)
call xyzzyaafd1(xyzzyaaag18,is,.false.,.true.,.false.,sd,.false.)
xyzzyaaee1=xyzzyaacf1(:,xyzzyaaaa18,:,1:real1_complex2,xyzzyaaag18,xyz&
&zyaaac18)
xyzzyaaee1(:,xyzzyaaae18,:)=xyzzyaach1(:,xyzzyaaaa18,:,xyzzyaaag18,is)
call xyzzyaagr1(xyzzyaaee1,xyzzyaadd1(xyzzyaaag18,is),xyzzyaacz1(:,xyz&
&zyaaaa18,:,xyzzyaaag18,is),xyzzyaadt1(:,xyzzyaaag18,ii,is))
xyzzyaadu1(xyzzyaaag18,ii,is)=.true.
enddo
endselect
case(2)
select case(xyzzyaaaf18)
case(1)
do xyzzyaaag18=1,xyzzyaaak1
if(xyzzyaadu1(xyzzyaaag18,ii,is))cycle
call xyzzyaafc1(xyzzyaaag18,xyzzyaaac18,xyzzyaaah18,.false.,.false.,.t&
&rue.,.false.,sd)
call xyzzyaafd1(xyzzyaaag18,is,.false.,.false.,.true.,.false.,sd)
xyzzyaaee1=xyzzyaacj1(:,:,xyzzyaaaa18,:,xyzzyaaag18,xyzzyaaac18)
xyzzyaaee1(:,xyzzyaaae18,:)=xyzzyaacl1(:,xyzzyaaaa18,:,xyzzyaaag18,is)
call xyzzyaagr1(xyzzyaaee1,xyzzyaadd1(xyzzyaaag18,is),xyzzyaacz1(xyzzy&
&aaaa18,:,:,xyzzyaaag18,is),xyzzyaadt1(:,xyzzyaaag18,ii,is))
xyzzyaadu1(xyzzyaaag18,ii,is)=.true.
enddo
case(2)
do xyzzyaaag18=1,xyzzyaaak1
if(xyzzyaadu1(xyzzyaaag18,ii,is))cycle
call xyzzyaafc1(xyzzyaaag18,xyzzyaaac18,xyzzyaaah18,.false.,.false.,.t&
&rue.,.false.,sd)
xyzzyaaee1=xyzzyaacj1(:,:,xyzzyaaaa18,:,xyzzyaaag18,xyzzyaaac18)
call xyzzyaagr1(xyzzyaaee1,xyzzyaadd1(xyzzyaaag18,is),xyzzyaacz1(xyzzy&
&aaaa18,:,:,xyzzyaaag18,is),xyzzyaadt1(:,xyzzyaaag18,ii,is))
xyzzyaadu1(xyzzyaaag18,ii,is)=.true.
enddo
endselect
endselect
endif
else
select case(xyzzyaaab18)
case(1)
do xyzzyaaag18=1,xyzzyaaak1
if(xyzzyaadu1(xyzzyaaag18,ii,is))cycle
call xyzzyaafc1(xyzzyaaag18,is,0,.false.,.true.,.false.,sd,.false.)
call xyzzyaagr1(xyzzyaacf1(:,xyzzyaaaa18,:,:,xyzzyaaag18,is),xyzzyaadd&
&1(xyzzyaaag18,is),xyzzyaacz1(:,xyzzyaaaa18,:,xyzzyaaag18,is),xyzzyaad&
&t1(:,xyzzyaaag18,ii,is))
xyzzyaadu1(xyzzyaaag18,ii,is)=.true.
enddo
case(2)
do xyzzyaaag18=1,xyzzyaaak1
if(xyzzyaadu1(xyzzyaaag18,ii,is))cycle
call xyzzyaafc1(xyzzyaaag18,is,0,.false.,.false.,.true.,.false.,sd)
call xyzzyaagr1(xyzzyaacj1(:,:,xyzzyaaaa18,:,xyzzyaaag18,is),xyzzyaadd&
&1(xyzzyaaag18,is),xyzzyaacz1(xyzzyaaaa18,:,:,xyzzyaaag18,is),xyzzyaad&
&t1(:,xyzzyaaag18,ii,is))
xyzzyaadu1(xyzzyaaag18,ii,is)=.true.
enddo
endselect
endif
call timer("GET_FIGEM",.false.)
end subroutine xyzzyaaes1
subroutine xyzzyaaet1(ii,is)
implicit none
integer,intent(in) :: ii,is
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,xy&
&zzyaaaf19,xyzzyaaag19,xyzzyaaah19
call timer("GET_LAPGEM",.true.)
xyzzyaaaa19=which_ie(ii)
xyzzyaaab19=which_spin(ii)
do xyzzyaaag19=1,xyzzyaaak1
if(.not.xyzzyaadu1(xyzzyaaag19,xyzzyaaaa19,is))then
call xyzzyaafb1(xyzzyaaag19,is)
endif
enddo
if(any(xyzzyaadd1(:,is)==xyzzyaaah1))return
xyzzyaaac19=buffer_move1_from(is)
if(xyzzyaaac19/=0)then
xyzzyaaad19=buffer_move1_from_ii(is)
xyzzyaaae19=which_ie(xyzzyaaad19)
xyzzyaaaf19=which_spin(xyzzyaaad19)
xyzzyaaah19=buffer_move1_from(xyzzyaaac19)
endif
if(xyzzyaaac19/=0)then
if(xyzzyaaaf19==xyzzyaaab19.and.xyzzyaaae19==xyzzyaaaa19)then
select case(xyzzyaaab19)
case(1)
do xyzzyaaag19=1,xyzzyaaak1
if(xyzzyaadw1(xyzzyaaag19,ii,is))cycle
call xyzzyaafd1(xyzzyaaag19,is,.false.,.false.,.false.,.true.,.false.)
call xyzzyaags1(xyzzyaacp1(:,:,xyzzyaaag19,is),xyzzyaadc1(xyzzyaaag19,&
&is),xyzzyaadd1(xyzzyaaag19,is),xyzzyaacz1(:,xyzzyaaaa19,:,xyzzyaaag19&
&,is),xyzzyaadv1(xyzzyaaag19,ii,is))
xyzzyaadw1(xyzzyaaag19,ii,is)=.true.
enddo
case(2)
do xyzzyaaag19=1,xyzzyaaak1
if(xyzzyaadw1(xyzzyaaag19,ii,is))cycle
call xyzzyaafd1(xyzzyaaag19,is,.false.,.false.,.false.,.false.,.true.)
call xyzzyaags1(xyzzyaact1(:,:,xyzzyaaag19,is),xyzzyaadc1(xyzzyaaag19,&
&is),xyzzyaadd1(xyzzyaaag19,is),xyzzyaacz1(xyzzyaaaa19,:,:,xyzzyaaag19&
&,is),xyzzyaadv1(xyzzyaaag19,ii,is))
xyzzyaadw1(xyzzyaaag19,ii,is)=.true.
enddo
endselect
else
select case(xyzzyaaab19)
case(1)
select case(xyzzyaaaf19)
case(1)
do xyzzyaaag19=1,xyzzyaaak1
if(xyzzyaadw1(xyzzyaaag19,ii,is))cycle
call xyzzyaafc1(xyzzyaaag19,xyzzyaaac19,xyzzyaaah19,.false.,.false.,.f&
&alse.,.true.,.false.)
xyzzyaaef1=xyzzyaacn1(xyzzyaaaa19,:,1:real1_complex2,xyzzyaaag19,xyzzy&
&aaac19)
call xyzzyaags1(xyzzyaaef1,xyzzyaadc1(xyzzyaaag19,is),xyzzyaadd1(xyzzy&
&aaag19,is),xyzzyaacz1(:,xyzzyaaaa19,:,xyzzyaaag19,is),xyzzyaadv1(xyzz&
&yaaag19,ii,is))
xyzzyaadw1(xyzzyaaag19,ii,is)=.true.
enddo
case(2)
do xyzzyaaag19=1,xyzzyaaak1
if(xyzzyaadw1(xyzzyaaag19,ii,is))cycle
call xyzzyaafc1(xyzzyaaag19,xyzzyaaac19,xyzzyaaah19,.false.,.false.,.f&
&alse.,.true.,.false.)
call xyzzyaafd1(xyzzyaaag19,is,.false.,.false.,.false.,.true.,.false.)
xyzzyaaef1=xyzzyaacn1(xyzzyaaaa19,:,1:real1_complex2,xyzzyaaag19,xyzzy&
&aaac19)
xyzzyaaef1(xyzzyaaae19,:)=xyzzyaacp1(xyzzyaaaa19,:,xyzzyaaag19,is)
call xyzzyaags1(xyzzyaaef1,xyzzyaadc1(xyzzyaaag19,is),xyzzyaadd1(xyzzy&
&aaag19,is),xyzzyaacz1(:,xyzzyaaaa19,:,xyzzyaaag19,is),xyzzyaadv1(xyzz&
&yaaag19,ii,is))
xyzzyaadw1(xyzzyaaag19,ii,is)=.true.
enddo
endselect
case(2)
select case(xyzzyaaaf19)
case(1)
do xyzzyaaag19=1,xyzzyaaak1
if(xyzzyaadw1(xyzzyaaag19,ii,is))cycle
call xyzzyaafc1(xyzzyaaag19,xyzzyaaac19,xyzzyaaah19,.false.,.false.,.f&
&alse.,.false.,.true.)
call xyzzyaafd1(xyzzyaaag19,is,.false.,.false.,.false.,.false.,.true.)
xyzzyaaef1=xyzzyaacr1(:,xyzzyaaaa19,:,xyzzyaaag19,xyzzyaaac19)
xyzzyaaef1(xyzzyaaae19,:)=xyzzyaact1(xyzzyaaaa19,:,xyzzyaaag19,is)
call xyzzyaags1(xyzzyaaef1,xyzzyaadc1(xyzzyaaag19,is),xyzzyaadd1(xyzzy&
&aaag19,is),xyzzyaacz1(xyzzyaaaa19,:,:,xyzzyaaag19,is),xyzzyaadv1(xyzz&
&yaaag19,ii,is))
xyzzyaadw1(xyzzyaaag19,ii,is)=.true.
enddo
case(2)
do xyzzyaaag19=1,xyzzyaaak1
if(xyzzyaadw1(xyzzyaaag19,ii,is))cycle
call xyzzyaafc1(xyzzyaaag19,xyzzyaaac19,xyzzyaaah19,.false.,.false.,.f&
&alse.,.false.,.true.)
xyzzyaaef1=xyzzyaacr1(:,xyzzyaaaa19,:,xyzzyaaag19,xyzzyaaac19)
call xyzzyaags1(xyzzyaaef1,xyzzyaadc1(xyzzyaaag19,is),xyzzyaadd1(xyzzy&
&aaag19,is),xyzzyaacz1(xyzzyaaaa19,:,:,xyzzyaaag19,is),xyzzyaadv1(xyzz&
&yaaag19,ii,is))
xyzzyaadw1(xyzzyaaag19,ii,is)=.true.
enddo
endselect
endselect
endif
else
select case(xyzzyaaab19)
case(1)
do xyzzyaaag19=1,xyzzyaaak1
if(xyzzyaadw1(xyzzyaaag19,ii,is))cycle
call xyzzyaafc1(xyzzyaaag19,is,0,.false.,.false.,.false.,.true.,.false&
&.)
call xyzzyaags1(xyzzyaacn1(xyzzyaaaa19,:,:,xyzzyaaag19,is),xyzzyaadc1(&
&xyzzyaaag19,is),xyzzyaadd1(xyzzyaaag19,is),xyzzyaacz1(:,xyzzyaaaa19,:&
&,xyzzyaaag19,is),xyzzyaadv1(xyzzyaaag19,ii,is))
xyzzyaadw1(xyzzyaaag19,ii,is)=.true.
enddo
case(2)
do xyzzyaaag19=1,xyzzyaaak1
if(xyzzyaadw1(xyzzyaaag19,ii,is))cycle
call xyzzyaafc1(xyzzyaaag19,is,0,.false.,.false.,.false.,.false.,.true&
&.)
call xyzzyaags1(xyzzyaacr1(:,xyzzyaaaa19,:,xyzzyaaag19,is),xyzzyaadc1(&
&xyzzyaaag19,is),xyzzyaadd1(xyzzyaaag19,is),xyzzyaacz1(xyzzyaaaa19,:,:&
&,xyzzyaaag19,is),xyzzyaadv1(xyzzyaaag19,ii,is))
xyzzyaadw1(xyzzyaaag19,ii,is)=.true.
enddo
endselect
endif
call timer("GET_LAPGEM",.false.)
end subroutine xyzzyaaet1
subroutine xyzzyaaeu1(is)
implicit none
integer,intent(in) :: is
integer xyzzyaaaa20,xyzzyaaab20,xyzzyaaac20,xyzzyaaad20,xyzzyaaae20,xy&
&zzyaaaf20,xyzzyaaag20,xyzzyaaah20,xyzzyaaai20
call timer("GET_HGEMS",.true.)
do xyzzyaaaa20=1,xyzzyaaak1
if(.not.xyzzyaady1(xyzzyaaaa20,is))then
call xyzzyaafh1(xyzzyaaaa20,is,.true.,.true.,.true.,.true.)
call xyzzyaafb1(xyzzyaaaa20,is)
endif
enddo
if(any(xyzzyaadd1(:,is)==xyzzyaaah1))return
do xyzzyaaaa20=1,xyzzyaaak1
if(xyzzyaady1(xyzzyaaaa20,is))cycle
xyzzyaadx1(:,:,:,:,:,xyzzyaaaa20,is)=0.d0
do xyzzyaaaf20=1,3
call d_or_z_gemm(complex_wf,'N','N',nemax,nemax,nemax,1.d0,xyzzyaacf1(&
&xyzzyaaaf20,:,:,re,xyzzyaaaa20,is),xyzzyaacf1(xyzzyaaaf20,:,:,im,xyzz&
&yaaaa20,is),nemax,xyzzyaacz1(:,:,re,xyzzyaaaa20,is),xyzzyaacz1(:,:,im&
&,xyzzyaaaa20,is),nemax,0.d0,xyzzyaael1(xyzzyaaaf20,:,:,re),xyzzyaael1&
&(xyzzyaaaf20,:,:,im),nemax)
enddo
do xyzzyaaab20=1,nemax
call xyzzyaaev1(xyzzyaaab20,xyzzyaaaa20,is)
do xyzzyaaac20=xyzzyaaab20+1,nemax
call d_or_z_ger(complex_wf,3,3,-1.d0,xyzzyaael1(:,xyzzyaaab20,xyzzyaaa&
&c20,re),xyzzyaael1(:,xyzzyaaab20,xyzzyaaac20,im),1,xyzzyaael1(:,xyzzy&
&aaac20,xyzzyaaab20,re),xyzzyaael1(:,xyzzyaaac20,xyzzyaaab20,im),1,xyz&
&zyaadx1(:,:,re,xyzzyaaab20,xyzzyaaac20,xyzzyaaaa20,is),xyzzyaadx1(:,:&
&,im,xyzzyaaab20,xyzzyaaac20,xyzzyaaaa20,is),3)
call d_or_z_ger(complex_wf,3,3,+1.d0,xyzzyaael1(:,xyzzyaaab20,xyzzyaaa&
&b20,re),xyzzyaael1(:,xyzzyaaab20,xyzzyaaab20,im),1,xyzzyaael1(:,xyzzy&
&aaac20,xyzzyaaac20,re),xyzzyaael1(:,xyzzyaaac20,xyzzyaaac20,im),1,xyz&
&zyaadx1(:,:,re,xyzzyaaab20,xyzzyaaac20,xyzzyaaaa20,is),xyzzyaadx1(:,:&
&,im,xyzzyaaab20,xyzzyaaac20,xyzzyaaaa20,is),3)
enddo
enddo
do xyzzyaaaf20=1,3
call d_or_z_gemm(complex_wf,'N','N',nemax,nemax,nemax,1.d0,xyzzyaacz1(&
&:,:,re,xyzzyaaaa20,is),xyzzyaacz1(:,:,im,xyzzyaaaa20,is),nemax,xyzzya&
&acj1(xyzzyaaaf20,:,:,re,xyzzyaaaa20,is),xyzzyaacj1(xyzzyaaaf20,:,:,im&
&,xyzzyaaaa20,is),nemax,0.d0,xyzzyaaem1(xyzzyaaaf20,:,:,re),xyzzyaaem1&
&(xyzzyaaaf20,:,:,im),nemax)
enddo
do xyzzyaaad20=1,nemax
xyzzyaaah20=which_ii(xyzzyaaad20,2)
call xyzzyaaew1(xyzzyaaad20,xyzzyaaaa20,is)
do xyzzyaaae20=xyzzyaaad20+1,nemax
xyzzyaaai20=which_ii(xyzzyaaae20,2)
call d_or_z_ger(complex_wf,3,3,-1.d0,xyzzyaaem1(:,xyzzyaaae20,xyzzyaaa&
&d20,re),xyzzyaaem1(:,xyzzyaaae20,xyzzyaaad20,im),1,xyzzyaaem1(:,xyzzy&
&aaad20,xyzzyaaae20,re),xyzzyaaem1(:,xyzzyaaad20,xyzzyaaae20,im),1,xyz&
&zyaadx1(:,:,re,xyzzyaaah20,xyzzyaaai20,xyzzyaaaa20,is),xyzzyaadx1(:,:&
&,im,xyzzyaaah20,xyzzyaaai20,xyzzyaaaa20,is),3)
call d_or_z_ger(complex_wf,3,3,+1.d0,xyzzyaaem1(:,xyzzyaaad20,xyzzyaaa&
&d20,re),xyzzyaaem1(:,xyzzyaaad20,xyzzyaaad20,im),1,xyzzyaaem1(:,xyzzy&
&aaae20,xyzzyaaae20,re),xyzzyaaem1(:,xyzzyaaae20,xyzzyaaae20,im),1,xyz&
&zyaadx1(:,:,re,xyzzyaaah20,xyzzyaaai20,xyzzyaaaa20,is),xyzzyaadx1(:,:&
&,im,xyzzyaaah20,xyzzyaaai20,xyzzyaaaa20,is),3)
enddo
enddo
do xyzzyaaaf20=1,3
do xyzzyaaag20=1,3
call d_or_z_gemm(complex_wf,'N','N',nemax,nemax,nemax,1.d0,xyzzyaacf1(&
&xyzzyaaaf20,:,:,re,xyzzyaaaa20,is),xyzzyaacf1(xyzzyaaaf20,:,:,im,xyzz&
&yaaaa20,is),nemax,xyzzyaaem1(xyzzyaaag20,:,:,re),xyzzyaaem1(xyzzyaaag&
&20,:,:,im),nemax,0.d0,xyzzyaadx1(xyzzyaaaf20,xyzzyaaag20,re,:,nele(1)&
&+1:,xyzzyaaaa20,is),xyzzyaadx1(xyzzyaaaf20,xyzzyaaag20,im,:,nele(1)+1&
&:,xyzzyaaaa20,is),netot)
enddo
enddo
do xyzzyaaab20=1,nemax
do xyzzyaaad20=1,nemax
xyzzyaaai20=which_ii(xyzzyaaad20,2)
xyzzyaadx1(:,:,:,xyzzyaaab20,xyzzyaaai20,xyzzyaaaa20,is)=xyzzyaacx1(:,&
&:,xyzzyaaab20,xyzzyaaad20,:,xyzzyaaaa20,is)-xyzzyaadx1(:,:,:,xyzzyaaa&
&b20,xyzzyaaai20,xyzzyaaaa20,is)
call d_or_z_scal(complex_wf,9,xyzzyaacz1(xyzzyaaad20,xyzzyaaab20,re,xy&
&zzyaaaa20,is),xyzzyaacz1(xyzzyaaad20,xyzzyaaab20,im,xyzzyaaaa20,is),x&
&yzzyaadx1(1,1,re,xyzzyaaab20,xyzzyaaai20,xyzzyaaaa20,is),xyzzyaadx1(1&
&,1,im,xyzzyaaab20,xyzzyaaai20,xyzzyaaaa20,is),1)
call d_or_z_ger(complex_wf,3,3,1.d0,xyzzyaael1(:,xyzzyaaab20,xyzzyaaab&
&20,re),xyzzyaael1(:,xyzzyaaab20,xyzzyaaab20,im),1,xyzzyaaem1(:,xyzzya&
&aad20,xyzzyaaad20,re),xyzzyaaem1(:,xyzzyaaad20,xyzzyaaad20,im),1,xyzz&
&yaadx1(:,:,re,xyzzyaaab20,xyzzyaaai20,xyzzyaaaa20,is),xyzzyaadx1(:,:,&
&im,xyzzyaaab20,xyzzyaaai20,xyzzyaaaa20,is),3)
enddo
enddo
call xyzzyaaex1(xyzzyaaaa20,is)
xyzzyaady1(xyzzyaaaa20,is)=.true.
enddo
call timer("GET_HGEMS",.false.)
end subroutine xyzzyaaeu1
subroutine xyzzyaaev1(iup,igem,is)
implicit none
integer,intent(in) :: iup,igem,is
integer xyzzyaaaa21
complex(dp) xyzzyaaab21
do xyzzyaaaa21=1,3
xyzzyaaab21=zdotu_cc(nemax,xyzzyaacv1(xyzzyaaaa21,iup,1,re,igem,is),xy&
&zzyaacv1(xyzzyaaaa21,iup,1,im,igem,is),nemax*6,xyzzyaacz1(1,iup,re,ig&
&em,is),xyzzyaacz1(1,iup,im,igem,is),1)
xyzzyaadx1(xyzzyaaaa21,xyzzyaaaa21,re,iup,iup,igem,is)=real(xyzzyaaab2&
&1)
xyzzyaadx1(xyzzyaaaa21,xyzzyaaaa21,im,iup,iup,igem,is)=aimag(xyzzyaaab&
&21)
enddo
xyzzyaaab21=zdotu_cc(nemax,xyzzyaacv1(4,iup,1,re,igem,is),xyzzyaacv1(4&
&,iup,1,im,igem,is),nemax*6,xyzzyaacz1(1,iup,re,igem,is),xyzzyaacz1(1,&
&iup,im,igem,is),1)
xyzzyaadx1(1,2,re,iup,iup,igem,is)=real(xyzzyaaab21)
xyzzyaadx1(1,2,im,iup,iup,igem,is)=aimag(xyzzyaaab21)
xyzzyaaab21=zdotu_cc(nemax,xyzzyaacv1(5,iup,1,re,igem,is),xyzzyaacv1(5&
&,iup,1,im,igem,is),nemax*6,xyzzyaacz1(1,iup,re,igem,is),xyzzyaacz1(1,&
&iup,im,igem,is),1)
xyzzyaadx1(1,3,re,iup,iup,igem,is)=real(xyzzyaaab21)
xyzzyaadx1(1,3,im,iup,iup,igem,is)=aimag(xyzzyaaab21)
xyzzyaaab21=zdotu_cc(nemax,xyzzyaacv1(6,iup,1,re,igem,is),xyzzyaacv1(6&
&,iup,1,im,igem,is),nemax*6,xyzzyaacz1(1,iup,re,igem,is),xyzzyaacz1(1,&
&iup,im,igem,is),1)
xyzzyaadx1(2,3,re,iup,iup,igem,is)=real(xyzzyaaab21)
xyzzyaadx1(2,3,im,iup,iup,igem,is)=aimag(xyzzyaaab21)
end subroutine xyzzyaaev1
subroutine xyzzyaaew1(ido,igem,is)
implicit none
integer,intent(in) :: ido,igem,is
integer xyzzyaaaa22,xyzzyaaab22
complex(dp) xyzzyaaac22
xyzzyaaab22=which_ii(ido,2)
do xyzzyaaaa22=1,3
xyzzyaaac22=zdotu_cc(nemax,xyzzyaacz1(ido,1,re,igem,is),xyzzyaacz1(ido&
&,1,im,igem,is),nemax,xyzzyaacw1(xyzzyaaaa22,1,ido,re,igem,is),xyzzyaa&
&cw1(xyzzyaaaa22,1,ido,im,igem,is),6)
xyzzyaadx1(xyzzyaaaa22,xyzzyaaaa22,re,xyzzyaaab22,xyzzyaaab22,igem,is)&
&=real(xyzzyaaac22)
xyzzyaadx1(xyzzyaaaa22,xyzzyaaaa22,im,xyzzyaaab22,xyzzyaaab22,igem,is)&
&=aimag(xyzzyaaac22)
enddo
xyzzyaaac22=zdotu_cc(nemax,xyzzyaacz1(ido,1,re,igem,is),xyzzyaacz1(ido&
&,1,im,igem,is),nemax,xyzzyaacw1(4,1,ido,re,igem,is),xyzzyaacw1(4,1,id&
&o,im,igem,is),6)
xyzzyaadx1(1,2,re,xyzzyaaab22,xyzzyaaab22,igem,is)=real(xyzzyaaac22)
xyzzyaadx1(1,2,im,xyzzyaaab22,xyzzyaaab22,igem,is)=aimag(xyzzyaaac22)
xyzzyaaac22=zdotu_cc(nemax,xyzzyaacz1(ido,1,re,igem,is),xyzzyaacz1(ido&
&,1,im,igem,is),nemax,xyzzyaacw1(5,1,ido,re,igem,is),xyzzyaacw1(5,1,id&
&o,im,igem,is),6)
xyzzyaadx1(1,3,re,xyzzyaaab22,xyzzyaaab22,igem,is)=real(xyzzyaaac22)
xyzzyaadx1(1,3,im,xyzzyaaab22,xyzzyaaab22,igem,is)=aimag(xyzzyaaac22)
xyzzyaaac22=zdotu_cc(nemax,xyzzyaacz1(ido,1,re,igem,is),xyzzyaacz1(ido&
&,1,im,igem,is),nemax,xyzzyaacw1(6,1,ido,re,igem,is),xyzzyaacw1(6,1,id&
&o,im,igem,is),6)
xyzzyaadx1(2,3,re,xyzzyaaab22,xyzzyaaab22,igem,is)=real(xyzzyaaac22)
xyzzyaadx1(2,3,im,xyzzyaaab22,xyzzyaaab22,igem,is)=aimag(xyzzyaaac22)
end subroutine xyzzyaaew1
subroutine xyzzyaaex1(igem,is)
implicit none
integer,intent(in) :: igem,is
integer xyzzyaaaa23,xyzzyaaab23
do xyzzyaaaa23=1,netot
xyzzyaadx1(2,1,1:real1_complex2,xyzzyaaaa23,xyzzyaaaa23,igem,is)=xyzzy&
&aadx1(1,2,1:real1_complex2,xyzzyaaaa23,xyzzyaaaa23,igem,is)
xyzzyaadx1(3,1,1:real1_complex2,xyzzyaaaa23,xyzzyaaaa23,igem,is)=xyzzy&
&aadx1(1,3,1:real1_complex2,xyzzyaaaa23,xyzzyaaaa23,igem,is)
xyzzyaadx1(3,2,1:real1_complex2,xyzzyaaaa23,xyzzyaaaa23,igem,is)=xyzzy&
&aadx1(2,3,1:real1_complex2,xyzzyaaaa23,xyzzyaaaa23,igem,is)
enddo
do xyzzyaaaa23=1,netot
do xyzzyaaab23=1,xyzzyaaaa23-1
xyzzyaadx1(:,:,re,xyzzyaaaa23,xyzzyaaab23,igem,is)=transpose(xyzzyaadx&
&1(:,:,re,xyzzyaaab23,xyzzyaaaa23,igem,is))
xyzzyaadx1(:,:,im,xyzzyaaaa23,xyzzyaaab23,igem,is)=transpose(xyzzyaadx&
&1(:,:,im,xyzzyaaab23,xyzzyaaaa23,igem,is))
enddo
enddo
end subroutine xyzzyaaex1
subroutine xyzzyaaey1(is)
implicit none
integer,intent(in) :: is
integer xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24
complex(dp) xyzzyaaad24,xyzzyaaae24
if(xyzzyaado1(is))return
call timer("GET_HARRAY",.true.)
call xyzzyaaeu1(is)
xyzzyaadn1(1:3,1:3,1:real1_complex2,1:netot,1:netot,is)=0.d0
do xyzzyaaaa24=1,xyzzyaaak1
xyzzyaaad24=exp(xyzzyaadc1(xyzzyaaaa24,is)-xyzzyaadi1(is))*xyzzyaaam1(&
&xyzzyaaaa24)
do xyzzyaaab24=1,netot
do xyzzyaaac24=1,netot
call zaxpy_cc(9,real(xyzzyaaad24),aimag(xyzzyaaad24),xyzzyaadx1(1,1,re&
&,xyzzyaaab24,xyzzyaaac24,xyzzyaaaa24,is),xyzzyaadx1(1,1,im,xyzzyaaab2&
&4,xyzzyaaac24,xyzzyaaaa24,is),1,xyzzyaadn1(1,1,re,xyzzyaaab24,xyzzyaa&
&ac24,is),xyzzyaadn1(1,1,im,xyzzyaaab24,xyzzyaaac24,is),1)
enddo
enddo
enddo
xyzzyaaae24=exp(-xyzzyaadh1(is))
if(complex_wf)then
do xyzzyaaab24=1,netot
do xyzzyaaac24=1,netot
call zscal_cc(9,real(xyzzyaaae24),aimag(xyzzyaaae24),xyzzyaadn1(1,1,re&
&,xyzzyaaab24,xyzzyaaac24,is),xyzzyaadn1(1,1,im,xyzzyaaab24,xyzzyaaac2&
&4,is),1)
enddo
enddo
else
call dscal(9*netot**2,real(xyzzyaaae24),xyzzyaadn1(1,1,re,1,1,is),1)
endif
xyzzyaado1(is)=.true.
call timer("GET_HARRAY",.false.)
end subroutine xyzzyaaey1
subroutine xyzzyaaez1(igem,is,lgrad,rgrad,llap,rlap)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: lgrad,rgrad,llap,rlap
integer xyzzyaaaa25,xyzzyaaab25,xyzzyaaac25,xyzzyaaad25
call timer("GET_DET_RATIO_IGEM",.true.)
xyzzyaaaa25=buffer_move1_from(is)
xyzzyaaab25=buffer_move1_from_ii(is)
xyzzyaaac25=which_ie(xyzzyaaab25)
xyzzyaaad25=which_spin(xyzzyaaab25)
if(xyzzyaaae1)then
if(.not.xyzzyaade1(igem,xyzzyaaaa25))then
call errstop_master("GET_DET_RATIO_IGEM","LOGGEM in the source buffer &
&is not valid. This is unexpected... Bug.")
endif
endif
if(xyzzyaadd1(igem,xyzzyaaaa25)/=xyzzyaaag1)then
call xyzzyaafc1(igem,is,xyzzyaaaa25,.true.,lgrad,rgrad,llap,rlap)
xyzzyaaed1=xyzzyaacb1(:,:,:,igem,is)
call xyzzyaagl1(nemax,xyzzyaadc1(igem,is),xyzzyaaed1,xyzzyaacz1(:,:,:,&
&igem,is),xyzzyaadd1(igem,is))
else
call xyzzyaafd1(igem,is,.true.,lgrad,rgrad,llap,rlap)
call xyzzyaagm1(xyzzyaaac25,xyzzyaaad25,xyzzyaacz1(:,:,:,igem,xyzzyaaa&
&a25),xyzzyaacd1(:,:,igem,is),xyzzyaadf1(igem,is),xyzzyaadd1(igem,is))
if(xyzzyaadd1(igem,is)/=xyzzyaaag1)then
xyzzyaadc1(igem,is)=czero
else
xyzzyaadc1(igem,is)=xyzzyaadc1(igem,xyzzyaaaa25)+log(xyzzyaadf1(igem,i&
&s))
if(real(xyzzyaadc1(igem,is),dp)<xyzzyaaaf1)then
xyzzyaadd1(igem,is)=xyzzyaaai1
xyzzyaadc1(igem,is)=czero
endif
endif
endif
xyzzyaadg1(igem,is)=.true.
xyzzyaade1(igem,is)=.true.
call timer("GET_DET_RATIO_IGEM",.false.)
end subroutine xyzzyaaez1
subroutine xyzzyaafa1(igem,is,lgrad,rgrad,llap,rlap)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: lgrad,rgrad,llap,rlap
if(xyzzyaadb1(igem,is).and.xyzzyaade1(igem,is))return
call timer("GET_GEM_INV_ALL",.true.)
call xyzzyaafc1(igem,is,0,.true.,lgrad,rgrad,llap,rlap)
xyzzyaaed1=xyzzyaacb1(:,:,:,igem,is)
call xyzzyaagl1(nemax,xyzzyaadc1(igem,is),xyzzyaaed1,xyzzyaacz1(:,:,:,&
&igem,is),xyzzyaadd1(igem,is))
if(xyzzyaadd1(igem,is)/=xyzzyaaag1)then
xyzzyaadc1(igem,is)=czero
endif
xyzzyaade1(igem,is)=.true.
xyzzyaadb1(igem,is)=.true.
xyzzyaada1(igem,is)=0
call timer("GET_GEM_INV_ALL",.false.)
end subroutine xyzzyaafa1
subroutine xyzzyaafb1(igem,is)
implicit none
integer,intent(in) :: igem,is
integer xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xyzzyaaad27
if(xyzzyaadb1(igem,is))then
if(xyzzyaaae1.and..not.xyzzyaade1(igem,is))then
call errstop_master("GET_GEM_INV","I think the inverse is valid and th&
&e determinant is not. I'm not sure this should happen. Bug.")
else
return
endif
endif
call timer("GET_GEM_INV",.true.)
xyzzyaaaa27=buffer_move1_from(is)
if(xyzzyaaaa27==0)then
call xyzzyaafa1(igem,is,.false.,.false.,.false.,.false.)
else
if(.not.xyzzyaadb1(igem,xyzzyaaaa27).or.xyzzyaada1(igem,xyzzyaaaa27)>=&
&xyzzyaaaj1)then
call xyzzyaafa1(igem,is,.false.,.false.,.false.,.false.)
else
xyzzyaaab27=buffer_move1_from_ii(is)
xyzzyaaac27=which_ie(xyzzyaaab27)
xyzzyaaad27=which_spin(xyzzyaaab27)
call xyzzyaaez1(igem,is,.false.,.false.,.false.,.false.)
if(xyzzyaadd1(igem,is)==xyzzyaaag1.and.xyzzyaadd1(igem,xyzzyaaaa27)==x&
&yzzyaaag1)then
xyzzyaacz1(:,:,:,igem,is)=xyzzyaacz1(:,:,:,igem,xyzzyaaaa27)
call xyzzyaagn1(xyzzyaaac27,xyzzyaaad27,xyzzyaacd1(:,:,igem,is),xyzzya&
&acz1(:,:,:,igem,is),xyzzyaadf1(igem,is))
xyzzyaadb1(igem,is)=.true.
xyzzyaada1(igem,is)=xyzzyaada1(igem,is)+1
endif
endif
endif
call timer("GET_GEM_INV",.false.)
end subroutine xyzzyaafb1
recursive subroutine xyzzyaafc1(xyzzyaaaa28,xyzzyaaab28,xyzzyaaac28,xy&
&zzyaaad28,xyzzyaaae28,xyzzyaaaf28,xyzzyaaag28,xyzzyaaah28)
implicit none
integer, intent(in) :: xyzzyaaaa28,xyzzyaaab28,xyzzyaaac28
logical, intent(in) :: xyzzyaaad28,xyzzyaaae28,xyzzyaaaf28,xyzzyaaag28&
&,xyzzyaaah28
integer xyzzyaaai28
call timer("UPDATE_MATRICES_IGEM",.true.)
if(xyzzyaaac28/=0.and..not.use_backflow)then
xyzzyaaai28=buffer_move1_from(xyzzyaaac28)
call xyzzyaafc1(xyzzyaaaa28,xyzzyaaac28,xyzzyaaai28,xyzzyaaad28,xyzzya&
&aae28,xyzzyaaaf28,xyzzyaaag28,xyzzyaaah28)
call xyzzyaafd1(xyzzyaaaa28,xyzzyaaab28,xyzzyaaad28,xyzzyaaae28,xyzzya&
&aaf28,xyzzyaaag28,xyzzyaaah28)
call xyzzyaafe1(xyzzyaaaa28,xyzzyaaab28,xyzzyaaad28,xyzzyaaae28,xyzzya&
&aaf28,xyzzyaaag28,xyzzyaaah28)
else
call xyzzyaaff1(xyzzyaaaa28,xyzzyaaab28,xyzzyaaad28,xyzzyaaae28,xyzzya&
&aaf28,xyzzyaaag28,xyzzyaaah28)
endif
call timer("UPDATE_MATRICES_IGEM",.false.)
end subroutine xyzzyaafc1
subroutine xyzzyaafd1(igem,is,val,lgrad,rgrad,llap,rlap)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: val,lgrad,rgrad,llap,rlap
integer xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaae29
logical xyzzyaaaf29,xyzzyaaag29,xyzzyaaah29,xyzzyaaai29,xyzzyaaaj29
call timer("GET_CHSCR_IGEM",.true.)
xyzzyaaaf29=val.and..not.xyzzyaace1(igem,is)
xyzzyaaag29=lgrad.and..not.xyzzyaaci1(igem,is)
xyzzyaaah29=rgrad.and..not.xyzzyaacm1(igem,is)
xyzzyaaai29=llap.and..not.xyzzyaacq1(igem,is)
xyzzyaaaj29=rlap.and..not.xyzzyaacu1(igem,is)
xyzzyaaaa29=buffer_move1_from(is)
xyzzyaaab29=buffer_move1_from_ii(is)
xyzzyaaac29=which_ie(xyzzyaaab29)
xyzzyaaad29=which_spin(xyzzyaaab29)
if(xyzzyaaad29==1)then
call xyzzyaafq1(xyzzyaaac29,xyzzyaaad29,is,xyzzyaaaf29.or.xyzzyaaah29.&
&or.xyzzyaaaj29,xyzzyaaag29.or.xyzzyaaai29)
call xyzzyaafl1(igem,xyzzyaaaa29,.false.)
if(xyzzyaaaf29)then
call d_or_z_gemv(complex_wf,'T',wfdet_norb,nemax,1.d0,xyzzyaabx1(:,:,r&
&e,igem,xyzzyaaaa29),xyzzyaabx1(:,:,im,igem,xyzzyaaaa29),wfdet_norb,xy&
&zzyaabp1(1,xyzzyaaac29,re,xyzzyaaad29,is),xyzzyaabp1(1,xyzzyaaac29,im&
&,xyzzyaaad29,is),1,0.d0,xyzzyaacd1(:,re,igem,is),xyzzyaacd1(:,im,igem&
&,is),1)
xyzzyaace1(igem,is)=.true.
endif
if(xyzzyaaag29)then
call d_or_z_gemm(complex_wf,'N','N',3,nemax,wfdet_norb,1.d0,xyzzyaabr1&
&(1,1,xyzzyaaac29,re,xyzzyaaad29,is),xyzzyaabr1(1,1,xyzzyaaac29,im,xyz&
&zyaaad29,is),3,xyzzyaabx1(1,1,re,igem,xyzzyaaaa29),xyzzyaabx1(1,1,im,&
&igem,xyzzyaaaa29),wfdet_norb,0.d0,xyzzyaach1(:,:,re,igem,is),xyzzyaac&
&h1(:,:,im,igem,is),3)
xyzzyaaci1(igem,is)=.true.
endif
if(xyzzyaaah29)then
call xyzzyaafn1(xyzzyaaac29,igem,is)
do xyzzyaaae29=1,dimensionality
call d_or_z_gemv(complex_wf,'T',wfdet_norb,nemax,1.d0,xyzzyaabr1(xyzzy&
&aaae29,:,:,re,2,xyzzyaaaa29),xyzzyaabr1(xyzzyaaae29,:,:,im,2,xyzzyaaa&
&a29),wfdet_norb,xyzzyaabz1(1,xyzzyaaac29,re,igem,is),xyzzyaabz1(1,xyz&
&zyaaac29,im,igem,is),1,0.d0,xyzzyaacl1(xyzzyaaae29,:,re,igem,is),xyzz&
&yaacl1(xyzzyaaae29,:,im,igem,is),1)
enddo
xyzzyaacm1(igem,is)=.true.
endif
if(xyzzyaaai29)then
call d_or_z_gemv(complex_wf,'T',wfdet_norb,nemax,1.d0,xyzzyaabx1(:,:,r&
&e,igem,xyzzyaaaa29),xyzzyaabx1(:,:,im,igem,xyzzyaaaa29),wfdet_norb,xy&
&zzyaabt1(1,xyzzyaaac29,re,xyzzyaaad29,is),xyzzyaabt1(1,xyzzyaaac29,im&
&,xyzzyaaad29,is),1,0.d0,xyzzyaacp1(:,re,igem,is),xyzzyaacp1(:,im,igem&
&,is),1)
xyzzyaacq1(igem,is)=.true.
endif
if(xyzzyaaaj29)then
call xyzzyaafn1(xyzzyaaac29,igem,is)
call d_or_z_gemv(complex_wf,'T',wfdet_norb,nemax,1.d0,xyzzyaabt1(:,:,r&
&e,2,xyzzyaaaa29),xyzzyaabt1(:,:,im,2,xyzzyaaaa29),wfdet_norb,xyzzyaab&
&z1(1,xyzzyaaac29,re,igem,is),xyzzyaabz1(1,xyzzyaaac29,im,igem,is),1,0&
&.d0,xyzzyaact1(:,re,igem,is),xyzzyaact1(:,im,igem,is),1)
xyzzyaacu1(igem,is)=.true.
endif
else
call xyzzyaafq1(xyzzyaaac29,xyzzyaaad29,is,xyzzyaaaf29.or.xyzzyaaag29.&
&or.xyzzyaaai29,xyzzyaaah29.or.xyzzyaaaj29)
call xyzzyaafm1(igem,xyzzyaaaa29,.false.)
if(xyzzyaaaf29)then
call d_or_z_gemv(complex_wf,'T',wfdet_norb,nemax,1.d0,xyzzyaabz1(:,:,r&
&e,igem,xyzzyaaaa29),xyzzyaabz1(:,:,im,igem,xyzzyaaaa29),wfdet_norb,xy&
&zzyaabp1(1,xyzzyaaac29,re,xyzzyaaad29,is),xyzzyaabp1(1,xyzzyaaac29,im&
&,xyzzyaaad29,is),1,0.d0,xyzzyaacd1(:,re,igem,is),xyzzyaacd1(:,im,igem&
&,is),1)
xyzzyaace1(igem,is)=.true.
endif
if(xyzzyaaag29)then
call xyzzyaafo1(xyzzyaaac29,igem,is)
do xyzzyaaae29=1,dimensionality
call d_or_z_gemv(complex_wf,'T',wfdet_norb,nemax,1.d0,xyzzyaabr1(xyzzy&
&aaae29,:,:,re,1,xyzzyaaaa29),xyzzyaabr1(xyzzyaaae29,:,:,im,1,xyzzyaaa&
&a29),wfdet_norb,xyzzyaabx1(1,xyzzyaaac29,re,igem,is),xyzzyaabx1(1,xyz&
&zyaaac29,im,igem,is),1,0.d0,xyzzyaach1(xyzzyaaae29,:,re,igem,is),xyzz&
&yaach1(xyzzyaaae29,:,im,igem,is),1)
enddo
xyzzyaaci1(igem,is)=.true.
endif
if(xyzzyaaah29)then
call d_or_z_gemm(complex_wf,'N','N',3,nemax,wfdet_norb,1.d0,xyzzyaabr1&
&(1,1,xyzzyaaac29,re,xyzzyaaad29,is),xyzzyaabr1(1,1,xyzzyaaac29,im,xyz&
&zyaaad29,is),3,xyzzyaabz1(1,1,re,igem,xyzzyaaaa29),xyzzyaabz1(1,1,im,&
&igem,xyzzyaaaa29),wfdet_norb,0.d0,xyzzyaacl1(:,:,re,igem,is),xyzzyaac&
&l1(:,:,im,igem,is),3)
xyzzyaacm1(igem,is)=.true.
endif
if(xyzzyaaai29)then
call xyzzyaafo1(xyzzyaaac29,igem,is)
call d_or_z_gemv(complex_wf,'T',wfdet_norb,nemax,1.d0,xyzzyaabt1(:,:,r&
&e,1,xyzzyaaaa29),xyzzyaabt1(:,:,im,1,xyzzyaaaa29),wfdet_norb,xyzzyaab&
&x1(1,xyzzyaaac29,re,igem,is),xyzzyaabx1(1,xyzzyaaac29,im,igem,is),1,0&
&.d0,xyzzyaacp1(:,re,igem,is),xyzzyaacp1(:,im,igem,is),1)
xyzzyaacq1(igem,is)=.true.
endif
if(xyzzyaaaj29)then
call d_or_z_gemv(complex_wf,'T',wfdet_norb,nemax,1.d0,xyzzyaabz1(:,:,r&
&e,igem,xyzzyaaaa29),xyzzyaabz1(:,:,im,igem,xyzzyaaaa29),wfdet_norb,xy&
&zzyaabt1(1,xyzzyaaac29,re,xyzzyaaad29,is),xyzzyaabt1(1,xyzzyaaac29,im&
&,xyzzyaaad29,is),1,0.d0,xyzzyaact1(:,re,igem,is),xyzzyaact1(:,im,igem&
&,is),1)
xyzzyaacu1(igem,is)=.true.
endif
endif
call timer("GET_CHSCR_IGEM",.false.)
end subroutine xyzzyaafd1
subroutine xyzzyaafe1(igem,is,val,lgrad,rgrad,llap,rlap)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: val,lgrad,rgrad,llap,rlap
integer xyzzyaaaa30,xyzzyaaab30,xyzzyaaac30,xyzzyaaad30
call timer("TRANSFER_MATRICES_IGEM",.true.)
xyzzyaaaa30=buffer_move1_from(is)
xyzzyaaab30=buffer_move1_from_ii(is)
xyzzyaaac30=which_ie(xyzzyaaab30)
xyzzyaaad30=which_spin(xyzzyaaab30)
if(val.and..not.xyzzyaacc1(igem,is))then
call xyzzyaago1(xyzzyaaac30,xyzzyaaad30,nemax,xyzzyaacb1(:,:,:,igem,xy&
&zzyaaaa30),xyzzyaacd1(:,:,igem,is),xyzzyaacb1(:,:,:,igem,is))
xyzzyaacc1(igem,is)=.true.
endif
if(lgrad.and..not.xyzzyaacg1(igem,is))then
call xyzzyaagp1(xyzzyaaac30,xyzzyaaad30,nemax,xyzzyaacf1(:,:,:,:,igem,&
&xyzzyaaaa30),xyzzyaach1(:,:,:,igem,is),xyzzyaacf1(:,:,:,:,igem,is))
xyzzyaacg1(igem,is)=.true.
endif
if(rgrad.and..not.xyzzyaack1(igem,is))then
call xyzzyaagp1(xyzzyaaac30,xyzzyaaad30,nemax,xyzzyaacj1(:,:,:,:,igem,&
&xyzzyaaaa30),xyzzyaacl1(:,:,:,igem,is),xyzzyaacj1(:,:,:,:,igem,is))
xyzzyaack1(igem,is)=.true.
endif
if(llap.and..not.xyzzyaaco1(igem,is))then
call xyzzyaago1(xyzzyaaac30,xyzzyaaad30,nemax,xyzzyaacn1(:,:,:,igem,xy&
&zzyaaaa30), xyzzyaacp1(:,:,igem,is),xyzzyaacn1(:,:,:,igem,is))
xyzzyaaco1(igem,is)=.true.
endif
if(rlap.and..not.xyzzyaacs1(igem,is))then
call xyzzyaago1(xyzzyaaac30,xyzzyaaad30,nemax,xyzzyaacr1(:,:,:,igem,xy&
&zzyaaaa30), xyzzyaact1(:,:,igem,is),xyzzyaacr1(:,:,:,igem,is))
xyzzyaacs1(igem,is)=.true.
endif
call timer("TRANSFER_MATRICES_IGEM",.false.)
end subroutine xyzzyaafe1
subroutine xyzzyaaff1(igem,is,val,lgrad,rgrad,llap,rlap)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: val,lgrad,rgrad,llap,rlap
if(use_backflow)then
call xyzzyaafh1(igem,is,val,lgrad,rgrad,llap.or.rlap)
else
call xyzzyaafg1(igem,is,val,lgrad,rgrad,llap,rlap)
endif
end subroutine xyzzyaaff1
subroutine xyzzyaafg1(igem,is,val,lgrad,rgrad,llap,rlap)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: val,lgrad,rgrad,llap,rlap
integer xyzzyaaaa32
logical xyzzyaaab32,xyzzyaaac32,xyzzyaaad32,xyzzyaaae32,xyzzyaaaf32,xy&
&zzyaaag32,xyzzyaaah32,xyzzyaaai32,xyzzyaaaj32,xyzzyaaak32
call timer("GET_MATRICES_IGEM_NOBF",.true.)
xyzzyaaag32=val.and..not.xyzzyaacc1(igem,is)
xyzzyaaah32=lgrad.and..not.xyzzyaacg1(igem,is)
xyzzyaaai32=rgrad.and..not.xyzzyaack1(igem,is)
xyzzyaaaj32=llap.and..not.xyzzyaaco1(igem,is)
xyzzyaaak32=rlap.and..not.xyzzyaacs1(igem,is)
xyzzyaaac32=xyzzyaaah32.or.xyzzyaaaj32
xyzzyaaad32=xyzzyaaag32.or.xyzzyaaac32
xyzzyaaae32=xyzzyaaai32.or.xyzzyaaak32
xyzzyaaaf32=xyzzyaaag32.or.xyzzyaaae32
do xyzzyaaaa32=1,nemax
call xyzzyaafp1(xyzzyaaaa32,1,is,xyzzyaaaf32,xyzzyaaac32)
call xyzzyaafp1(xyzzyaaaa32,2,is,xyzzyaaad32,xyzzyaaae32)
enddo
if(xyzzyaaag32)then
if(.not.xyzzyaaae32.and..not.xyzzyaaac32)then
xyzzyaaab32=count(xyzzyaaby1(:,igem,is))>count(xyzzyaaca1(:,igem,is))
elseif(xyzzyaaae32)then
xyzzyaaab32=.false.
else
xyzzyaaab32=.true.
endif
call xyzzyaafi1(igem,is,xyzzyaaab32)
endif
if(xyzzyaaah32)call xyzzyaafj1(igem,is)
if(xyzzyaaai32)call xyzzyaafk1(igem,is)
if(xyzzyaaaj32)then
call xyzzyaafl1(igem,is,.false.)
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abt1(:,:,re,1,is),xyzzyaabt1(:,:,im,1,is),wfdet_norb,xyzzyaabx1(:,:,r&
&e,igem,is),xyzzyaabx1(:,:,im,igem,is),wfdet_norb,0.d0,xyzzyaacn1(:,:,&
&re,igem,is),xyzzyaacn1(:,:,im,igem,is),nemax)
xyzzyaaco1(igem,is)=.true.
endif
if(xyzzyaaak32)then
call xyzzyaafm1(igem,is,.false.)
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abz1(:,:,re,igem,is),xyzzyaabz1(:,:,im,igem,is),wfdet_norb,xyzzyaabt1&
&(:,:,re,2,is),xyzzyaabt1(:,:,im,2,is),wfdet_norb,0.d0,xyzzyaacr1(:,:,&
&re,igem,is),xyzzyaacr1(:,:,im,igem,is),nemax)
xyzzyaacs1(igem,is)=.true.
endif
call timer("GET_MATRICES_IGEM_NOBF",.false.)
end subroutine xyzzyaafg1
subroutine xyzzyaafh1(igem,is,val,lgrad,rgrad,sderivs)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: val,lgrad,rgrad,sderivs
integer xyzzyaaaa33,xyzzyaaab33,xyzzyaaac33
logical xyzzyaaad33,xyzzyaaae33,xyzzyaaaf33,xyzzyaaag33,xyzzyaaah33
call timer("GET_MATRICES_IGEM_BF",.true.)
xyzzyaaad33=val.and..not.xyzzyaacc1(igem,is)
xyzzyaaae33=lgrad.and..not.xyzzyaacg1(igem,is)
xyzzyaaaf33=rgrad.and..not.xyzzyaack1(igem,is)
xyzzyaaag33=sderivs.and..not.xyzzyaacy1(igem,is)
call get_bf_x(is,xyzzyaaad33,xyzzyaaae33.or.xyzzyaaaf33.or.xyzzyaaag33&
&,xyzzyaaag33)
do xyzzyaaaa33=1,nemax
call xyzzyaafp1(xyzzyaaaa33,1,is,xyzzyaaad33.or.xyzzyaaaf33.or.xyzzyaa&
&ag33,xyzzyaaae33.or.xyzzyaaag33)
call xyzzyaafp1(xyzzyaaaa33,2,is,xyzzyaaad33.or.xyzzyaaae33.or.xyzzyaa&
&ag33,xyzzyaaaf33.or.xyzzyaaag33)
enddo
if(xyzzyaaad33)then
if(.not.xyzzyaaaf33.and..not.xyzzyaaae33)then
xyzzyaaah33=count(xyzzyaaby1(:,igem,is))>count(xyzzyaaca1(:,igem,is))
elseif(xyzzyaaaf33)then
xyzzyaaah33=.false.
else
xyzzyaaah33=.true.
endif
call xyzzyaafi1(igem,is,xyzzyaaah33)
endif
if(xyzzyaaae33)call xyzzyaafj1(igem,is)
if(xyzzyaaaf33)call xyzzyaafk1(igem,is)
if(xyzzyaaag33)then
call xyzzyaafl1(igem,is,.true.)
do xyzzyaaab33=1,6
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abv1(xyzzyaaab33,:,:,re,1,is),xyzzyaabv1(xyzzyaaab33,:,:,im,1,is),wfd&
&et_norb,xyzzyaabx1(:,:,re,igem,is),xyzzyaabx1(:,:,im,igem,is),wfdet_n&
&orb,0.d0,xyzzyaacv1(xyzzyaaab33,:,:,re,igem,is),xyzzyaacv1(xyzzyaaab3&
&3,:,:,im,igem,is),nemax)
enddo
do xyzzyaaab33=1,6
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abz1(:,:,re,igem,is),xyzzyaabz1(:,:,im,igem,is),wfdet_norb,xyzzyaabv1&
&(xyzzyaaab33,:,:,re,2,is),xyzzyaabv1(xyzzyaaab33,:,:,im,2,is),wfdet_n&
&orb,0.d0,xyzzyaacw1(xyzzyaaab33,:,:,re,igem,is),xyzzyaacw1(xyzzyaaab3&
&3,:,:,im,igem,is),nemax)
enddo
do xyzzyaaab33=1,3
do xyzzyaaaa33=1,nemax
call xyzzyaagq1(xyzzyaaal1(:,:,igem),xyzzyaaao1(:,:,igem),xyzzyaaaq1(:&
&,igem),xyzzyaabr1(xyzzyaaab33,1,xyzzyaaaa33,re,2,is),xyzzyaabr1(xyzzy&
&aaab33,1,xyzzyaaaa33,im,2,is),3,xyzzyaaej1(1,xyzzyaaaa33,re),xyzzyaae&
&j1(1,xyzzyaaaa33,im),1)
enddo
do xyzzyaaac33=1,3
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abr1(xyzzyaaac33,:,:,re,1,is),xyzzyaabr1(xyzzyaaac33,:,:,im,1,is),wfd&
&et_norb,xyzzyaaej1(1,1,re),xyzzyaaej1(1,1,im),wfdet_norb,0.d0,xyzzyaa&
&cx1(xyzzyaaac33,xyzzyaaab33,:,:,re,igem,is),xyzzyaacx1(xyzzyaaac33,xy&
&zzyaaab33,:,:,im,igem,is),nemax)
enddo
enddo
xyzzyaacy1(igem,is)=.true.
endif
call timer("GET_MATRICES_IGEM_BF",.false.)
end subroutine xyzzyaafh1
subroutine xyzzyaafi1(igem,is,prefer_gphi)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: prefer_gphi
if(prefer_gphi)then
call xyzzyaafl1(igem,is,.false.)
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abp1(1,1,re,1,is),xyzzyaabp1(1,1,im,1,is),wfdet_norb,xyzzyaabx1(1,1,r&
&e,igem,is),xyzzyaabx1(1,1,im,igem,is),wfdet_norb,0.d0,xyzzyaacb1(1,1,&
&re,igem,is),xyzzyaacb1(1,1,im,igem,is),nemax)
else
call xyzzyaafm1(igem,is,.false.)
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abz1(1,1,re,igem,is),xyzzyaabz1(1,1,im,igem,is),wfdet_norb,xyzzyaabp1&
&(1,1,re,2,is),xyzzyaabp1(1,1,im,2,is),wfdet_norb,0.d0,xyzzyaacb1(1,1,&
&re,igem,is),xyzzyaacb1(1,1,im,igem,is),nemax)
endif
xyzzyaacc1(igem,is)=.true.
end subroutine xyzzyaafi1
subroutine xyzzyaafj1(igem,is)
implicit none
integer,intent(in) :: igem,is
integer xyzzyaaaa35
call xyzzyaafl1(igem,is,.false.)
do xyzzyaaaa35=1,dimensionality
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abr1(xyzzyaaaa35,:,:,re,1,is),xyzzyaabr1(xyzzyaaaa35,:,:,im,1,is),wfd&
&et_norb,xyzzyaabx1(:,:,re,igem,is),xyzzyaabx1(:,:,im,igem,is),wfdet_n&
&orb,0.d0,xyzzyaacf1(xyzzyaaaa35,:,:,re,igem,is),xyzzyaacf1(xyzzyaaaa3&
&5,:,:,im,igem,is),nemax)
enddo
xyzzyaacg1(igem,is)=.true.
end subroutine xyzzyaafj1
subroutine xyzzyaafk1(igem,is)
implicit none
integer,intent(in) :: igem,is
integer xyzzyaaaa36
call xyzzyaafm1(igem,is,.false.)
do xyzzyaaaa36=1,dimensionality
call d_or_z_gemm(complex_wf,'T','N',nemax,nemax,wfdet_norb,1.d0,xyzzya&
&abz1(:,:,re,igem,is),xyzzyaabz1(:,:,im,igem,is),wfdet_norb,xyzzyaabr1&
&(xyzzyaaaa36,:,:,re,2,is),xyzzyaabr1(xyzzyaaaa36,:,:,im,2,is),wfdet_n&
&orb,0.d0,xyzzyaacj1(xyzzyaaaa36,:,:,re,igem,is),xyzzyaacj1(xyzzyaaaa3&
&6,:,:,im,igem,is),nemax)
enddo
xyzzyaack1(igem,is)=.true.
end subroutine xyzzyaafk1
subroutine xyzzyaafl1(igem,is,grad_lap)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: grad_lap
integer xyzzyaaaa37
call timer("GET_GPHI",.true.)
do xyzzyaaaa37=1,nele(2)
call xyzzyaafp1(xyzzyaaaa37,2,is,.true.,grad_lap)
call xyzzyaafo1(xyzzyaaaa37, igem, is)
enddo
call timer("GET_GPHI",.false.)
end subroutine xyzzyaafl1
subroutine xyzzyaafm1(igem,is,grad_lap)
implicit none
integer,intent(in) :: igem,is
logical,intent(in) :: grad_lap
integer xyzzyaaaa38
call timer("GET_PHIG",.true.)
do xyzzyaaaa38=1,nele(1)
call xyzzyaafp1(xyzzyaaaa38,1,is,.true.,grad_lap)
call xyzzyaafn1(xyzzyaaaa38,igem,is)
enddo
call timer("GET_PHIG",.false.)
end subroutine xyzzyaafm1
subroutine xyzzyaafn1(ie,igem,is)
implicit none
integer,intent(in) :: ie,igem,is
if(xyzzyaaca1(ie,igem,is))return
call timer("GET_PHIG1",.true.)
if(xyzzyaaae1)then
if(.not.xyzzyaabq1(ie,1,is))call errstop_master("GET_PHIG1","Bug.")
endif
if(.not.xyzzyaaan1)then
call dgemv('T',wfdet_norb,wfdet_norb,1.d0,xyzzyaaal1(:,:,igem),wfdet_n&
&orb,xyzzyaabp1(:,ie,re,1,is),1,0.d0,xyzzyaabz1(:,ie,re,igem,is),1)
if(complex_wf)then
call dgemv('T',wfdet_norb,wfdet_norb,1.d0,xyzzyaaal1(:,:,igem),wfdet_n&
&orb,xyzzyaabp1(:,ie,im,1,is),1,0.d0,xyzzyaabz1(:,ie,im,igem,is),1)
endif
else
if(complex_wf)then
call d2_sparse_vector_matrix(wfdet_norb,1.d0,xyzzyaaal1(:,:,igem),wfde&
&t_norb,xyzzyaaap1(:,:,igem),xyzzyaaat1,xyzzyaaar1(:,igem),xyzzyaabp1(&
&:,ie,re,1,is),xyzzyaabp1(:,ie,im,1,is),1,xyzzyaabz1(1,ie,re,igem,is),&
&xyzzyaabz1(1,ie,im,igem,is),1)
else
call d_sparse_vector_matrix(wfdet_norb,1.d0,xyzzyaaal1(:,:,igem),wfdet&
&_norb,xyzzyaaap1(:,:,igem),xyzzyaaat1,xyzzyaaar1(:,igem),xyzzyaabp1(:&
&,ie,re,1,is),1,xyzzyaabz1(1,ie,re,igem,is),1)
endif
endif
xyzzyaaca1(ie,igem,is)=.true.
call timer("GET_PHIG1",.false.)
end subroutine xyzzyaafn1
subroutine xyzzyaafo1(ie,igem,is)
implicit none
integer,intent(in) :: ie,igem,is
if(xyzzyaaby1(ie,igem,is))return
call timer("GET_GPHI1",.true.)
if(xyzzyaaae1)then
if(.not.xyzzyaabq1(ie,2,is))call errstop_master("GET_GPHI1","Bug.")
endif
call xyzzyaagq1(xyzzyaaal1(:,:,igem),xyzzyaaao1(:,:,igem),xyzzyaaaq1(:&
&,igem),xyzzyaabp1(1,ie,re,2,is),xyzzyaabp1(1,ie,im,2,is),1,xyzzyaabx1&
&(1,ie,re,igem,is),xyzzyaabx1(1,ie,im,igem,is),1)
xyzzyaaby1(ie,igem,is)=.true.
call timer("GET_GPHI1",.false.)
end subroutine xyzzyaafo1
subroutine xyzzyaafp1(i,spin,is,val,fsd)
implicit none
integer,intent(in) :: i,spin,is
logical,intent(in) :: val,fsd
integer xyzzyaaaa41
logical xyzzyaaab41,xyzzyaaac41
if(fsd.and.(.not.xyzzyaabs1(i,spin,is).or..not.xyzzyaabu1(i,spin,is)))&
&then
xyzzyaaab41=.true.
else
xyzzyaaab41=.false.
endif
if(val.and..not.xyzzyaabq1(i,spin,is))then
xyzzyaaac41=.true.
else
xyzzyaaac41=.false.
endif
if(.not.xyzzyaaac41.and..not.xyzzyaaab41)return
call timer("GET_ORBVALS",.true.)
xyzzyaaaa41=which_ii(i,spin)
if(.not.use_backflow)then
call xyzzyaagt1(rele_scr(:,xyzzyaaaa41,is),spin,xyzzyaaac41,xyzzyaaab4&
&1,xyzzyaabp1(:,i,:,spin,is),xyzzyaabr1(:,:,i,:,spin,is),xyzzyaabt1(:,&
&i,:,spin,is))
else
if(.not.xyzzyaaab41)then
call xyzzyaagt1(bf_x_scr(:,xyzzyaaaa41,is),spin,xyzzyaaac41,xyzzyaaab4&
&1,xyzzyaabp1(:,i,:,spin,is),xyzzyaabr1(:,:,i,:,spin,is),xyzzyaabt1(:,&
&i,:,spin,is))
else
call xyzzyaagt1(bf_x_scr(:,xyzzyaaaa41,is),spin,xyzzyaaac41,xyzzyaaab4&
&1,xyzzyaabp1(:,i,:,spin,is),xyzzyaabr1(:,:,i,:,spin,is),xyzzyaabt1(:,&
&i,:,spin,is),xyzzyaabv1(:,:,i,:,spin,is))
endif
endif
if(xyzzyaaac41)xyzzyaabq1(i,spin,is)=.true.
if(xyzzyaaab41)then
xyzzyaabs1(i,spin,is)=.true.
xyzzyaabu1(i,spin,is)=.true.
if(use_backflow)xyzzyaabw1(i,spin,is)=.true.
endif
call timer("GET_ORBVALS",.false.)
end subroutine xyzzyaafp1
subroutine xyzzyaafq1(ie,ispin,is,val,fsd)
implicit none
integer,intent(in) :: ie,ispin,is
logical,intent(in) :: val,fsd
logical xyzzyaaaa42,xyzzyaaab42
xyzzyaaaa42=.false.
xyzzyaaab42=.false.
if(val.and..not.xyzzyaabq1(ie,ispin,is))xyzzyaaaa42=.true.
if(fsd.and.(.not.xyzzyaabs1(ie,ispin,is).or..not.xyzzyaabu1(ie,ispin,i&
&s)))xyzzyaaab42=.true.
if(.not.xyzzyaaaa42.and..not.xyzzyaaab42)return
call timer("GET_ORBVALS1_CHSCR",.true.)
call xyzzyaagt1(rele1_chscr(:,is),ispin,xyzzyaaaa42,xyzzyaaab42,xyzzya&
&abp1(:,ie,:,ispin,is),xyzzyaabr1(:,:,ie,:,ispin,is),xyzzyaabt1(:,ie,:&
&,ispin,is))
if(xyzzyaaaa42)xyzzyaabq1(ie,ispin,is)=.true.
if(xyzzyaaab42)then
xyzzyaabs1(ie,ispin,is)=.true.
xyzzyaabu1(ie,ispin,is)=.true.
endif
call timer("GET_ORBVALS1_CHSCR",.false.)
end subroutine xyzzyaafq1
subroutine setup_geminal_params(nparam)
implicit none
integer,intent(inout) :: nparam
nparam=0
xyzzyaabb1=0
xyzzyaabc1=0
xyzzyaabd1=0
if(opt_geminal)then
xyzzyaabb1=count(xyzzyaabf1==xyzzyaaaw1)
xyzzyaabc1=count(xyzzyaabg1==xyzzyaaaw1)
call xyzzyaagg1
call xyzzyaagh1
endif
if(use_backflow)then
call setup_bf_params(xyzzyaabd1)
endif
xyzzyaabe1=xyzzyaabb1+xyzzyaabc1+xyzzyaabd1
nparam=xyzzyaabe1
end subroutine setup_geminal_params
subroutine finish_geminal_params
implicit none
if(opt_geminal)then
deallocate(xyzzyaabh1)
deallocate(xyzzyaabo1)
endif
call finish_bf_params
end subroutine finish_geminal_params
subroutine get_geminal_params(params,has_lolim,lolim,has_hilim,hilim,i&
&s_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,l&
&abel)
implicit none
real(dp),intent(inout) :: params(xyzzyaabe1),lolim(xyzzyaabe1),hilim(x&
&yzzyaabe1)
logical,intent(inout) :: has_lolim(xyzzyaabe1),has_hilim(xyzzyaabe1),i&
&s_shallow(xyzzyaabe1),is_redundant(xyzzyaabe1),is_linear(xyzzyaabe1),&
&is_loglinear(xyzzyaabe1),has_aderiv(xyzzyaabe1),affect_map(xyzzyaabe1&
&,xyzzyaabe1)
character(2),intent(inout) :: label(xyzzyaabe1)
integer xyzzyaaaa45,xyzzyaaab45,xyzzyaaac45,xyzzyaaad45,xyzzyaaae45,xy&
&zzyaaaf45
if(.not.(opt_geminal.or.opt_backflow))return
params=0.d0
has_lolim=.false.
has_hilim=.false.
hilim=0.d0
lolim =0.d0
is_shallow=.false.
is_redundant=.false.
is_linear=.false.
is_loglinear=.false.
has_aderiv=.false.
affect_map=.true.
label='GP'
do xyzzyaaaa45=1,xyzzyaabb1
xyzzyaaab45=xyzzyaabh1(1,xyzzyaaaa45)
xyzzyaaac45=xyzzyaabh1(2,xyzzyaaaa45)
xyzzyaaad45=xyzzyaabh1(3,xyzzyaaaa45)
params(xyzzyaaaa45)=xyzzyaaal1(xyzzyaaac45,xyzzyaaad45,xyzzyaaab45)
enddo
params(xyzzyaabb1+1:xyzzyaabb1+xyzzyaabc1)=pack(xyzzyaaam1,xyzzyaabg1=&
&=xyzzyaaaw1)
if(xyzzyaabd1>0)then
xyzzyaaae45=xyzzyaabb1+xyzzyaabc1+1
xyzzyaaaf45=xyzzyaabb1+xyzzyaabc1+xyzzyaabd1
call get_bf_params(params(xyzzyaaae45:xyzzyaaaf45),has_lolim(xyzzyaaae&
&45:xyzzyaaaf45),lolim(xyzzyaaae45:xyzzyaaaf45),has_hilim(xyzzyaaae45:&
&xyzzyaaaf45),hilim(xyzzyaaae45:xyzzyaaaf45),is_shallow(xyzzyaaae45:xy&
&zzyaaaf45),is_redundant(xyzzyaaae45:xyzzyaaaf45),is_linear(xyzzyaaae4&
&5:xyzzyaaaf45),is_loglinear(xyzzyaaae45:xyzzyaaaf45),has_aderiv(xyzzy&
&aaae45:xyzzyaaaf45),affect_map(xyzzyaaae45:xyzzyaaaf45,xyzzyaaae45:xy&
&zzyaaaf45),label(xyzzyaaae45:xyzzyaaaf45))
endif
end subroutine get_geminal_params
subroutine put_geminal_params(params,ignore,iparam_buffer,prestore,bad&
&_params)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaabe1)
logical,intent(in) :: ignore(xyzzyaabe1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa46,xyzzyaaab46,xyzzyaaac46,xyzzyaaad46,xyzzyaaae46,xy&
&zzyaaaf46,xyzzyaaag46
if(.not.(opt_geminal.or.opt_backflow))return
bad_params=.false.
do xyzzyaaad46=1,xyzzyaabb1
if(ignore(xyzzyaaad46))cycle
xyzzyaaac46=xyzzyaabh1(1,xyzzyaaad46)
xyzzyaaaa46=xyzzyaabh1(2,xyzzyaaad46)
xyzzyaaab46=xyzzyaabh1(3,xyzzyaaad46)
xyzzyaaal1(xyzzyaaaa46,xyzzyaaab46,xyzzyaaac46)=params(xyzzyaaad46)
enddo
do xyzzyaaac46=1,xyzzyaaak1
if(xyzzyaabg1(xyzzyaaac46)/=xyzzyaaaw1)cycle
xyzzyaaad46=xyzzyaaad46+1
if(ignore(xyzzyaaad46))cycle
xyzzyaaam1(xyzzyaaac46)=params(xyzzyaaad46)
enddo
call xyzzyaagi1
if(xyzzyaabd1>0)then
xyzzyaaae46=xyzzyaabb1+xyzzyaabc1+1
xyzzyaaaf46=xyzzyaabb1+xyzzyaabc1+xyzzyaabd1
if(iparam_buffer>=xyzzyaaae46)then
xyzzyaaag46=iparam_buffer-xyzzyaabb1-xyzzyaabc1
else
xyzzyaaag46=0
endif
call put_bf_params(params(xyzzyaaae46:xyzzyaaaf46),ignore(xyzzyaaae46:&
&xyzzyaaaf46),xyzzyaaag46,prestore,bad_params)
endif
end subroutine put_geminal_params
subroutine invalidate_param1_geminal(is,iparam)
implicit none
integer,intent(in) :: is,iparam
integer xyzzyaaaa47
if(iparam<=xyzzyaabb1)then
do xyzzyaaaa47=1,xyzzyaaak1
if(xyzzyaabo1(xyzzyaaaa47,iparam))then
call xyzzyaafr1(is,xyzzyaaaa47,1,1)
endif
enddo
call xyzzyaafs1(is)
elseif(iparam<=xyzzyaabc1)then
call xyzzyaafs1(is)
else
call clear_scratch_geminal(is)
call invalidate_param1_bf(is,iparam)
endif
end subroutine invalidate_param1_geminal
subroutine xyzzyaafr1(is,igem,irow,icol)
implicit none
integer,intent(in) :: is,igem,irow,icol
xyzzyaaby1(:,igem,is)=.false.
xyzzyaaca1(:,igem,is)=.false.
xyzzyaacc1(igem,is)=.false.
xyzzyaacg1(igem,is)=.false.
xyzzyaack1(igem,is)=.false.
xyzzyaaco1(igem,is)=.false.
xyzzyaacs1(igem,is)=.false.
xyzzyaace1(igem,is)=.false.
xyzzyaaci1(igem,is)=.false.
xyzzyaacm1(igem,is)=.false.
xyzzyaacq1(igem,is)=.false.
xyzzyaacu1(igem,is)=.false.
xyzzyaadb1(igem,is)=.false.
xyzzyaade1(igem,is)=.false.
xyzzyaadg1(igem,is)=.false.
xyzzyaadu1(igem,:,is)=.false.
if(.not.use_backflow)then
xyzzyaadw1(igem,:,is)=.false.
else
xyzzyaadm1(:,is)=.false.
xyzzyaado1(is)=.false.
endif
end subroutine xyzzyaafr1
subroutine xyzzyaafs1(is)
implicit none
integer,intent(in) :: is
xyzzyaadj1(is)=.false.
xyzzyaadq1(:,is)=.false.
xyzzyaads1(:,is)=.false.
end subroutine xyzzyaafs1
subroutine invalidate_params_geminal(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaabe1)
end subroutine invalidate_params_geminal
subroutine read_geminal
implicit none
integer xyzzyaaaa51,xyzzyaaab51,xyzzyaaac51
logical exists,is_block,xyzzyaaad51
if(am_master)then
call wout('Geminal setup')
call wout('=============')
call wout()
endif
call query_casl_item('parameters.casl:GEMINAL',exists=exists,is_block=&
&is_block,nchildren=xyzzyaaac51)
call push_casl_context('parameters.casl:GEMINAL')
call xyzzyaaft1(xyzzyaaac51,xyzzyaaad51,xyzzyaaaz1,xyzzyaaba1)
allocate(xyzzyaaal1(wfdet_norb,wfdet_norb,xyzzyaaak1),xyzzyaabf1(wfdet&
&_norb,wfdet_norb,xyzzyaaak1),xyzzyaaeb1(wfdet_norb,wfdet_norb,xyzzyaa&
&ak1),stat=xyzzyaaab51)
call check_alloc(xyzzyaaab51,"READ_GEMINAL","gmat")
xyzzyaaal1=0.d0
xyzzyaabf1=xyzzyaaau1
xyzzyaaeb1=.false.
allocate(xyzzyaaam1(xyzzyaaak1),xyzzyaabg1(xyzzyaaak1),stat=xyzzyaaab5&
&1)
call check_alloc(xyzzyaaab51,"READ_GEMINAL","c_coeffs")
xyzzyaaam1=0.d0
xyzzyaabg1=xyzzyaaau1
do xyzzyaaaa51=1,xyzzyaaak1
call xyzzyaafu1(xyzzyaaaa51)
enddo
if(xyzzyaaad51)call xyzzyaafz1
call xyzzyaagi1
call xyzzyaafx1
call xyzzyaafy1
call xyzzyaagk1
if(xyzzyaaad1==xyzzyaaab1)then
if(am_master)then
call wout("The geminal matrices will be represented by dense matrices.&
&")
endif
xyzzyaaan1=.false.
else
call xyzzyaagf1
endif
call pop_casl_context()
end subroutine read_geminal
subroutine xyzzyaaft1(nlabels,are_constraints,xyzzyaaaz1,xyzzyaaba1)
implicit none
integer,intent(in) :: nlabels
integer,intent(out) :: xyzzyaaaz1,xyzzyaaba1
logical,intent(out) :: are_constraints
integer xyzzyaaaa52,xyzzyaaab52,xyzzyaaac52,xyzzyaaad52,xyzzyaaae52
character(casl_keysize) label,buffer
character(casl_valsize) value_buffer
are_constraints=.false.
xyzzyaaaz1=xyzzyaaau1
xyzzyaaba1=xyzzyaaau1
xyzzyaaak1=0
xyzzyaaac52=0
do xyzzyaaaa52=1,nlabels
call first_unread_child(':parameters.casl:GEMINAL',label,xyzzyaaab52,f&
&lag_as_read=.true.)
xyzzyaaae52=index(trim(label),unique_casl_string("Geminal"))
if(xyzzyaaae52/=0)then
xyzzyaaac52=xyzzyaaac52+1
buffer=trim(label(xyzzyaaae52+len("geminal"):))
read(buffer,*)xyzzyaaad52
xyzzyaaak1=max(xyzzyaaad52,xyzzyaaak1)
elseif(trim(label)==unique_casl_string("Constraints"))then
are_constraints=.true.
elseif(trim(label)==unique_casl_string("Default g optimizability"))the&
&n
call get_casl_item("Default g optimizability",value_buffer,xyzzyaaab52&
&)
xyzzyaaaz1=xyzzyaagv1(value_buffer)
elseif(trim(label)==unique_casl_string("Default c optimizability"))the&
&n
call get_casl_item("Default c optimizability",value_buffer,xyzzyaaab52&
&)
xyzzyaaba1=xyzzyaagv1(value_buffer)
else
call errstop_master("PARSE_TOP_LEVEL","Syntax error in GEMINAL block. &
&Label '"//trim(label)//"' is invalid.")
endif
enddo
if(xyzzyaaac52==0)then
call errstop_master("PARSE_TOP_LEVEL","Invalid syntax: no geminals fou&
&nd under GEMINAL block.")
endif
if(xyzzyaaak1/=xyzzyaaac52)then
call errstop_master("PARSE_TOP_LEVEL","Invalid syntax: some geminals a&
&re not specified explicitly. Please enter a geminal block for each ge&
&minal. The block may be empty.")
endif
if(xyzzyaaaz1==xyzzyaaau1)xyzzyaaaz1=xyzzyaaay1
if(xyzzyaaba1==xyzzyaaau1)xyzzyaaba1=xyzzyaaay1
end subroutine xyzzyaaft1
subroutine xyzzyaafu1(gem_index)
implicit none
integer,intent(in) :: gem_index
integer xyzzyaaaa53
logical is_block,exists
character(casl_keysize) label
call query_casl_item("Geminal "//trim(i2s(gem_index)),is_block=is_bloc&
&k)
if(.not.is_block)call errstop_master("PARSE_GEMINAL_BLOCK","Invalid sy&
&ntax: 'Geminal "//trim(i2s(gem_index))//"' should be a block. An empt&
&y block may specified using 'Geminal "//trim(i2s(gem_index))//" : [ ]&
&'.")
call query_casl_item("Geminal "//trim(i2s(gem_index))//":Parameters",e&
&xists=exists,is_block=is_block)
if(.not.exists)return
if(.not.is_block)then
call errstop_master("PARSE_GEMINAL_BLOCK","Invalid syntax in definitio&
&n of 'Parameters' section under geminal "//trim(i2s(gem_index))//". S&
&hould be a block (or absent).")
endif
call push_casl_context("Geminal "//trim(i2s(gem_index))//":Parameters"&
&)
do
call first_unread_child(':parameters.casl:GEMINAL:Geminal '//trim(i2s(&
&gem_index))//":Parameters",label,xyzzyaaaa53)
if(xyzzyaaaa53/=0)exit
if(scan(label,"g")/=0)then
call xyzzyaafv1(label,gem_index)
elseif(scan(label,"c")/=0)then
call xyzzyaafw1(gem_index)
else
call errstop_master("PARSE_GEMINAL_BLOCK","Invalid syntax in the GEMIN&
&AL block in the 'parameters.casl' input file. Keyword '"//trim(label)&
&//"' is not a valid keyword.")
endif
enddo
call pop_casl_context()
end subroutine xyzzyaafu1
subroutine xyzzyaafv1(casl_label,igem)
implicit none
integer,intent(in) :: igem
character(*),intent(in) :: casl_label
integer xyzzyaaaa54,xyzzyaaab54,xyzzyaaac54
character(casl_valsize) buffer
call xyzzyaagw1(casl_label,xyzzyaaab54,xyzzyaaac54)
if(xyzzyaaab54>wfdet_norb)then
call errstop_master('PARSE_GMAT_EL','Error while reading label '//trim&
&(casl_label)//'. The row index is larger than the total number of orb&
&itals allowed ('//trim(i2s(wfdet_norb))//'). Update the number of orb&
&itals allowed in the input file.')
endif
if(xyzzyaaab54<1)then
call errstop_master('PARSE_GMAT_EL','Error while reading label '//trim&
&(casl_label)//'. The row index is less than 1.')
endif
if(xyzzyaaac54>wfdet_norb)then
call errstop_master('PARSE_GMAT_EL','Error while reading label '//trim&
&(casl_label)//'. The column index is larger than the total number of &
&orbitals allowed ('//trim(i2s(wfdet_norb))//'). Update the number of &
&orbitals allowed in the input file.')
endif
if(xyzzyaaac54<1)then
call errstop_master('PARSE_GMAT_EL','Error while reading label '//trim&
&(casl_label)//'. The column index is less than 1.')
endif
call get_casl_item(trim(casl_label)//":%u1",xyzzyaaal1(xyzzyaaab54,xyz&
&zyaaac54,igem),xyzzyaaaa54)
if(xyzzyaaaa54/=0)then
call errstop_master('PARSE_GMAT_EL',"Parameter "//trim(casl_label)//" &
&has no specified value.")
endif
call get_casl_item(trim(casl_label)//":%u2",buffer,xyzzyaaaa54)
if(xyzzyaaaa54==0)then
xyzzyaabf1(xyzzyaaab54,xyzzyaaac54,igem)=xyzzyaagu1(buffer)
else
xyzzyaabf1(xyzzyaaab54,xyzzyaaac54,igem)=xyzzyaaaz1
endif
end subroutine xyzzyaafv1
subroutine xyzzyaafw1(gem)
implicit none
integer,intent(in) :: gem
integer xyzzyaaaa55
logical exists
character(casl_valsize) buffer
call push_casl_context("c")
call get_casl_item("%u1",xyzzyaaam1(gem),xyzzyaaaa55)
if(xyzzyaaaa55/=0)then
call errstop_master("PARSE_C_COEFF","Parameter c under geminal "//trim&
&(i2s(gem))//" has no specified value.")
endif
call query_casl_item("%u2",exists=exists)
if(exists)then
call get_casl_item("%u2",buffer,xyzzyaaaa55)
xyzzyaabg1(gem)=xyzzyaagu1(buffer)
else
xyzzyaabg1=xyzzyaaba1
endif
call pop_casl_context()
end subroutine xyzzyaafw1
subroutine xyzzyaafx1
implicit none
integer xyzzyaaaa56,xyzzyaaab56,xyzzyaaac56
do xyzzyaaaa56=1,xyzzyaaak1
do xyzzyaaac56=1,wfdet_norb
do xyzzyaaab56=1,wfdet_norb
if(xyzzyaabf1(xyzzyaaab56,xyzzyaaac56,xyzzyaaaa56)/=xyzzyaaau1)cycle
if(xyzzyaaab56==xyzzyaaac56.and.xyzzyaaab56<=nemax.and.xyzzyaaaa56==1)&
&then
xyzzyaaal1(xyzzyaaab56,xyzzyaaac56,xyzzyaaaa56)=0.d0
else
xyzzyaaal1(xyzzyaaab56,xyzzyaaac56,xyzzyaaaa56)=0.d0
endif
xyzzyaabf1(xyzzyaaab56,xyzzyaaac56,xyzzyaaaa56)=xyzzyaaaz1
enddo
enddo
enddo
end subroutine xyzzyaafx1
subroutine xyzzyaafy1
implicit none
if(xyzzyaabg1(1)==xyzzyaaau1)then
xyzzyaaam1(1)=1.d0
xyzzyaabg1(1)=xyzzyaaba1
endif
where(xyzzyaabg1(2:xyzzyaaak1)==xyzzyaaau1)
xyzzyaaam1(2:xyzzyaaak1)=0.d0
xyzzyaabg1(2:xyzzyaaak1)=xyzzyaaba1
endwhere
end subroutine xyzzyaafy1
subroutine xyzzyaafz1
implicit none
integer xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58,xyzzyaaae58,xy&
&zzyaaaf58,xyzzyaaag58,xyzzyaaah58
logical exists,is_block
character(casl_valsize) constraint
call query_casl_item("Constraints",exists=exists,is_block=is_block,nch&
&ildren=xyzzyaaaa58)
call push_casl_context("Constraints")
call query_casl_item("Equate",exists=exists)
if(exists)then
call errstop_master("PARSE_CONSTRAINTS","Syntax error in 'GEMINAL:Cons&
&traints' section of the parameters.casl file. All 'Equate' blocks sho&
&uld be followed by a numeral: eg. Equate 1 : [ ... ]")
endif
call query_casl_item("Symmetrize",exists=exists)
if(exists)then
call errstop_master("PARSE_CONSTRAINTS","Syntax error in 'GEMINAL:Cons&
&traints' section of the parameters.casl file. All 'Symmetrize' blocks&
& should be followed by a numeral: eg. Symmetrize 1 : [ ... ]")
endif
xyzzyaaab58=0
do
call query_casl_item("Equate "//trim(i2s(xyzzyaaab58+1)),exists=exists&
&,is_block=is_block)
if(exists)then
xyzzyaaab58=xyzzyaaab58+1
call xyzzyaagc1(xyzzyaaab58)
else
exit
endif
enddo
xyzzyaaac58=0
do
call query_casl_item("Symmetrize "//trim(i2s(xyzzyaaac58+1)),exists=ex&
&ists,is_block=is_block)
if(exists)then
if(.not.is_block)call errstop_master("PARSE_CONSTRAINTS","Syntax error&
& in 'Symmetrize "//trim(i2s(xyzzyaaac58+1))//"'. Should be a block. T&
&o specify an empty block, use 'Symmetrize "//trim(i2s(xyzzyaaac58+1))&
&//" : [ ]'.")
xyzzyaaac58=xyzzyaaac58+1
call xyzzyaagd1(xyzzyaaac58)
else
exit
endif
enddo
xyzzyaabi1=0
xyzzyaabj1=0
do xyzzyaaad58=1,xyzzyaaaa58-xyzzyaaab58-xyzzyaaac58
call get_casl_item("%u"//trim(i2s(xyzzyaaad58)),constraint,xyzzyaaae58&
&)
if(scan(constraint,"g")/=0)then
xyzzyaabi1=xyzzyaabi1+1
elseif(scan(constraint,"c")/=0)then
xyzzyaabj1=xyzzyaabj1+1
else
call errstop_master("PARSE_CONSTRAINTS","Invalid syntax: Constraint '"&
&//trim(constraint)//"' is not correctly specified.")
endif
enddo
call resize_pointer((/xyzzyaabi1/),xyzzyaabm1,0)
call resize_pointer((/xyzzyaabj1/),xyzzyaabn1,0)
xyzzyaaaf58=1
xyzzyaaag58=1
do xyzzyaaad58=1,xyzzyaaaa58-xyzzyaaab58-xyzzyaaac58
call get_casl_item("%u"//trim(i2s(xyzzyaaad58)),constraint,xyzzyaaae58&
&)
if(scan(constraint,"g")/=0)then
xyzzyaaah58=count_char(constraint,"=")+1
xyzzyaabm1(xyzzyaaaf58)=xyzzyaaah58
xyzzyaaaf58=xyzzyaaaf58+1
elseif(scan(constraint,"c")/=0)then
xyzzyaaah58=count_char(constraint,"=")+1
xyzzyaabn1(xyzzyaaag58)=xyzzyaaah58
xyzzyaaag58=xyzzyaaag58+1
else
call errstop_master("PARSE_CONSTRAINTS","Invalid syntax: Constraint '"&
&//trim(constraint)//"' is not correctly specified.")
endif
enddo
call resize_pointer((/3,maxval(xyzzyaabm1),xyzzyaabi1/),xyzzyaabk1)
call resize_pointer((/maxval(xyzzyaabn1),xyzzyaabj1/),xyzzyaabl1)
xyzzyaaaf58=1
xyzzyaaag58=1
do xyzzyaaad58=1,xyzzyaaaa58-xyzzyaaab58-xyzzyaaac58
call get_casl_item("%u"//trim(i2s(xyzzyaaad58)),constraint,xyzzyaaae58&
&)
if(scan(constraint,"g")/=0)then
call xyzzyaaga1(constraint,xyzzyaaaf58)
xyzzyaaaf58=xyzzyaaaf58+1
elseif(scan(constraint,"c")/=0)then
call xyzzyaagb1(constraint,xyzzyaaag58)
xyzzyaaag58=xyzzyaaag58+1
else
call errstop_master("PARSE_CONSTRAINTS","Invalid syntax: Constraint '"&
&//trim(constraint)//"' is not correctly specified.")
endif
enddo
call pop_casl_context()
end subroutine xyzzyaafz1
subroutine xyzzyaaga1(constraint,iconstraint)
implicit none
character(*),intent(in) :: constraint
integer, intent(in) :: iconstraint
integer xyzzyaaaa59,xyzzyaaab59,xyzzyaaac59,xyzzyaaad59,xyzzyaaae59,xy&
&zzyaaaf59,xyzzyaaag59
xyzzyaaaa59=1
xyzzyaaad59=xyzzyaabm1(iconstraint)
do xyzzyaaac59=1,xyzzyaaad59-1
xyzzyaaab59=scan(constraint(xyzzyaaaa59:),"=")+xyzzyaaaa59-2
call xyzzyaagx1(constraint(xyzzyaaaa59:xyzzyaaab59),xyzzyaaae59,xyzzya&
&aaf59,xyzzyaaag59)
xyzzyaabk1(:,xyzzyaaac59,iconstraint)=(/xyzzyaaae59,xyzzyaaaf59,xyzzya&
&aag59/)
xyzzyaaaa59=xyzzyaaab59+2
enddo
call xyzzyaagx1(constraint(xyzzyaaaa59:),xyzzyaaae59,xyzzyaaaf59,xyzzy&
&aaag59)
xyzzyaabk1(:,xyzzyaaad59,iconstraint)=(/xyzzyaaae59,xyzzyaaaf59,xyzzya&
&aag59/)
call xyzzyaage1(iconstraint)
end subroutine xyzzyaaga1
subroutine xyzzyaagb1(constraint,iconstraint)
implicit none
integer,intent(in) :: iconstraint
character(*),intent(in) :: constraint
integer xyzzyaaaa60,xyzzyaaab60,xyzzyaaac60,xyzzyaaad60,xyzzyaaae60,xy&
&zzyaaaf60
xyzzyaaad60=xyzzyaabn1(iconstraint)
xyzzyaaaa60=1
do xyzzyaaac60=1,xyzzyaaad60-1
xyzzyaaab60=scan(constraint(xyzzyaaaa60:),"=")+xyzzyaaaa60-2
call xyzzyaagy1(constraint(xyzzyaaaa60:xyzzyaaab60),xyzzyaaae60)
xyzzyaabl1(xyzzyaaac60,iconstraint)=xyzzyaaae60
xyzzyaaaa60=xyzzyaaab60+2
enddo
call xyzzyaagy1(constraint(xyzzyaaaa60:),xyzzyaaae60)
xyzzyaabl1(xyzzyaaad60,iconstraint)=xyzzyaaae60
xyzzyaaaf60=0
do xyzzyaaac60=1,xyzzyaaad60
xyzzyaaae60=xyzzyaabl1(xyzzyaaac60,iconstraint)
select case(xyzzyaabg1(xyzzyaaae60))
case(xyzzyaaau1)
xyzzyaabg1(xyzzyaaae60)=xyzzyaaax1
case(xyzzyaaax1)
cycle
case default
if(xyzzyaaaf60/=0)then
call errstop_master("PARSE_C_CONSTRAINTS","Contradiction in the 'Const&
&raints' section of the 'GEMINAL' block. Parameter "//trim(i2s(xyzzyaa&
&ae60))//":c is "//trim(opt2s(xyzzyaabg1(xyzzyaaae60)))//" whereas par&
&ameter "//trim(i2s(xyzzyaabl1(xyzzyaaaf60,iconstraint)))//":c is "//t&
&rim(opt2s(xyzzyaabg1(xyzzyaabl1(xyzzyaaaf60,iconstraint))))//".")
else
xyzzyaaaf60=xyzzyaaac60
endif
endselect
enddo
if(xyzzyaaaf60/=0)then
call swap1(xyzzyaabl1(1,iconstraint),xyzzyaabl1(xyzzyaaaf60,iconstrain&
&t))
else
xyzzyaabg1(xyzzyaabl1(1,iconstraint))=xyzzyaaba1
endif
end subroutine xyzzyaagb1
subroutine xyzzyaagc1(iequate)
implicit none
integer,intent(in) :: iequate
integer xyzzyaaaa61,xyzzyaaab61,xyzzyaaac61,xyzzyaaad61,xyzzyaaae61,xy&
&zzyaaaf61,xyzzyaaag61,xyzzyaaah61,xyzzyaaai61,xyzzyaaaj61,xyzzyaaak61&
&,xyzzyaaal61,xyzzyaaam61,xyzzyaaan61,xyzzyaaao61,xyzzyaaap61
logical xyzzyaaaq61,xyzzyaaar61,xyzzyaaas61,xyzzyaaat61
character(casl_valsize) buffer
call push_casl_context("Equate "//trim(i2s(iequate)))
call query_casl_item('Diagonal',exists=xyzzyaaar61)
call query_casl_item('Off-diagonal',exists=xyzzyaaas61)
call query_casl_item('All',exists=xyzzyaaat61)
if(.not.(xyzzyaaar61.or.xyzzyaaas61.or.xyzzyaaat61))call errstop_maste&
&r('PARSE_G_EQUATE','Equate '//trim(i2s(iequate))//' does not specify &
&a Diagonal/Off-diagonal/All keyword.')
if((xyzzyaaar61.and.xyzzyaaas61).or.(xyzzyaaar61.and.xyzzyaaat61).or.(&
&xyzzyaaas61.and.xyzzyaaat61))call errstop_master('PARSE_G_EQUATE','Eq&
&uate '//trim(i2s(iequate))//' specifies more than one of Diagonal/Off&
&-diagonal/All keywords.')
call get_casl_item("Geminals",buffer,xyzzyaaaj61)
if(xyzzyaaaj61==0)then
call xyzzyaaha1(buffer,xyzzyaaaa61,xyzzyaaab61,xyzzyaaaq61)
xyzzyaaac61=xyzzyaaab61-xyzzyaaaa61+1
else
xyzzyaaaa61=1
xyzzyaaab61=xyzzyaaak1
xyzzyaaac61=xyzzyaaak1
endif
xyzzyaaao61=0
xyzzyaabi1=xyzzyaabi1+1
if(xyzzyaaat61)then
call get_casl_item('All',buffer,xyzzyaaaj61)
call xyzzyaagz1(buffer,xyzzyaaad61,xyzzyaaae61,xyzzyaaag61,xyzzyaaah61&
&)
xyzzyaaaf61=xyzzyaaae61-xyzzyaaad61+1
xyzzyaaai61=xyzzyaaah61-xyzzyaaag61+1
xyzzyaaap61=xyzzyaaac61*xyzzyaaaf61*xyzzyaaai61
call resize_pointer((/xyzzyaabi1/),xyzzyaabm1,xyzzyaaap61)
call resize_pointer((/3,maxval(xyzzyaabm1),xyzzyaabi1/),xyzzyaabk1)
do xyzzyaaan61=xyzzyaaag61,xyzzyaaah61
do xyzzyaaam61=xyzzyaaad61,xyzzyaaae61
do xyzzyaaal61=xyzzyaaaa61,xyzzyaaab61
xyzzyaaao61=xyzzyaaao61+1
xyzzyaabk1(:,xyzzyaaao61,xyzzyaabi1)=(/xyzzyaaal61,xyzzyaaam61,xyzzyaa&
&an61/)
enddo
enddo
enddo
elseif(xyzzyaaar61)then
call get_casl_item('Diagonal',buffer,xyzzyaaaj61)
call xyzzyaaha1(buffer,xyzzyaaag61,xyzzyaaah61,xyzzyaaaq61)
xyzzyaaai61=xyzzyaaah61-xyzzyaaag61+1
xyzzyaaap61=xyzzyaaac61*xyzzyaaai61
call resize_pointer((/xyzzyaabi1/),xyzzyaabm1,xyzzyaaap61)
call resize_pointer((/3,maxval(xyzzyaabm1),xyzzyaabi1/),xyzzyaabk1)
do xyzzyaaal61=xyzzyaaaa61,xyzzyaaab61
do xyzzyaaan61=xyzzyaaag61,xyzzyaaah61
xyzzyaaao61=xyzzyaaao61+1
xyzzyaabk1(:,xyzzyaaao61,xyzzyaabi1)=(/xyzzyaaal61,xyzzyaaan61,xyzzyaa&
&an61/)
enddo
enddo
elseif(xyzzyaaas61)then
call get_casl_item('Off-diagonal',buffer,xyzzyaaaj61)
call xyzzyaagz1(buffer,xyzzyaaad61,xyzzyaaae61,xyzzyaaag61,xyzzyaaah61&
&)
xyzzyaaaf61=xyzzyaaae61-xyzzyaaad61+1
xyzzyaaai61=xyzzyaaah61-xyzzyaaag61+1
xyzzyaaak61=min(xyzzyaaae61,xyzzyaaah61)-max(xyzzyaaad61,xyzzyaaag61)+&
&1
xyzzyaaap61=xyzzyaaac61*(xyzzyaaaf61*xyzzyaaai61-xyzzyaaak61)
call resize_pointer((/xyzzyaabi1/),xyzzyaabm1,xyzzyaaap61)
call resize_pointer((/3,maxval(xyzzyaabm1),xyzzyaabi1/),xyzzyaabk1)
do xyzzyaaan61=xyzzyaaag61,xyzzyaaah61
do xyzzyaaam61=xyzzyaaad61,xyzzyaaae61
if(xyzzyaaam61==xyzzyaaan61)cycle
do xyzzyaaal61=xyzzyaaaa61,xyzzyaaab61
xyzzyaaao61=xyzzyaaao61+1
xyzzyaabk1(:,xyzzyaaao61,xyzzyaabi1)=(/xyzzyaaal61,xyzzyaaam61,xyzzyaa&
&an61/)
enddo
enddo
enddo
endif
call pop_casl_context()
call xyzzyaage1(xyzzyaabi1)
end subroutine xyzzyaagc1
subroutine xyzzyaagd1(isymmetrize)
implicit none
integer,intent(in) :: isymmetrize
integer xyzzyaaaa62,xyzzyaaab62,xyzzyaaac62,xyzzyaaad62,xyzzyaaae62,xy&
&zzyaaaf62,xyzzyaaag62,xyzzyaaah62,xyzzyaaai62,xyzzyaaaj62,xyzzyaaak62
logical xyzzyaaal62
character(casl_valsize) buffer
call push_casl_context("Symmetrize "//trim(i2s(isymmetrize)))
call get_casl_item("Geminals",buffer,xyzzyaaag62)
if(xyzzyaaag62==0)then
call xyzzyaaha1(buffer,xyzzyaaaa62,xyzzyaaab62,xyzzyaaal62)
xyzzyaaac62=xyzzyaaab62-xyzzyaaaa62+1
else
xyzzyaaaa62=1
xyzzyaaab62=xyzzyaaak1
xyzzyaaac62=xyzzyaaak1
endif
call get_casl_item("All",buffer,xyzzyaaag62)
if(xyzzyaaag62==0)then
call xyzzyaaha1(buffer,xyzzyaaad62,xyzzyaaae62,xyzzyaaal62)
xyzzyaaaf62=xyzzyaaae62-xyzzyaaad62+1
else
xyzzyaaad62=1
xyzzyaaae62=wfdet_norb
xyzzyaaaf62=wfdet_norb
endif
xyzzyaaak62=xyzzyaabi1
xyzzyaabi1=xyzzyaabi1+xyzzyaaac62*((xyzzyaaaf62*(xyzzyaaaf62-1))/2)
call resize_pointer((/xyzzyaabi1/),xyzzyaabm1,2)
call resize_pointer((/3,maxval(xyzzyaabm1),xyzzyaabi1/),xyzzyaabk1)
do xyzzyaaah62=xyzzyaaad62,xyzzyaaae62
do xyzzyaaai62=xyzzyaaad62,xyzzyaaah62-1
do xyzzyaaaj62=xyzzyaaaa62,xyzzyaaab62
xyzzyaaak62=xyzzyaaak62+1
xyzzyaabk1(:,1,xyzzyaaak62)=(/xyzzyaaaj62,xyzzyaaai62,xyzzyaaah62/)
xyzzyaabk1(:,2,xyzzyaaak62)=(/xyzzyaaaj62,xyzzyaaah62,xyzzyaaai62/)
call xyzzyaage1(xyzzyaaak62)
enddo
enddo
enddo
call pop_casl_context()
end subroutine xyzzyaagd1
subroutine xyzzyaage1(iconstraint)
implicit none
integer,intent(in) :: iconstraint
integer xyzzyaaaa63,xyzzyaaab63,xyzzyaaac63,xyzzyaaad63,xyzzyaaae63
xyzzyaaab63=0
do xyzzyaaaa63=1,xyzzyaabm1(iconstraint)
xyzzyaaac63=xyzzyaabk1(1,xyzzyaaaa63,iconstraint)
xyzzyaaad63=xyzzyaabk1(2,xyzzyaaaa63,iconstraint)
xyzzyaaae63=xyzzyaabk1(3,xyzzyaaaa63,iconstraint)
select case(xyzzyaabf1(xyzzyaaad63,xyzzyaaae63,xyzzyaaac63))
case(xyzzyaaau1)
xyzzyaabf1(xyzzyaaad63,xyzzyaaae63,xyzzyaaac63)=xyzzyaaax1
case(xyzzyaaax1)
cycle
case default
if(xyzzyaaab63/=0)then
call errstop_master("PARSE_G_CONSTRAINTS","Contradiction in the 'Const&
&raints' section of the 'GEMINAL' block. Parameter "//trim(inds2extgla&
&bel(xyzzyaaac63,xyzzyaaad63,xyzzyaaae63))//" is "//trim(opt2s(xyzzyaa&
&bf1(xyzzyaaad63,xyzzyaaae63,xyzzyaaac63)))//" whereas parameter "//tr&
&im(inds2extglabel(xyzzyaabk1(1,xyzzyaaab63,iconstraint),xyzzyaabk1(2,&
&xyzzyaaab63,iconstraint),xyzzyaabk1(3,xyzzyaaab63,iconstraint)))//" i&
&s "//trim(opt2s(xyzzyaabf1(xyzzyaabk1(2,xyzzyaaab63,iconstraint),xyzz&
&yaabk1(3,xyzzyaaab63,iconstraint),xyzzyaabk1(1,xyzzyaaab63,iconstrain&
&t))))//".")
else
xyzzyaaab63=xyzzyaaaa63
endif
endselect
enddo
if(xyzzyaaab63/=0)then
call iswap(3,xyzzyaabk1(:,1,iconstraint),1,xyzzyaabk1(:,xyzzyaaab63,ic&
&onstraint),1)
else
xyzzyaabf1(xyzzyaabk1(2,1,iconstraint),xyzzyaabk1(3,1,iconstraint),xyz&
&zyaabk1(1,1,iconstraint))=xyzzyaaba1
endif
end subroutine xyzzyaage1
subroutine xyzzyaagf1
implicit none
integer xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64,xyzzyaaad64,xyzzyaaae64,xy&
&zzyaaaf64,xyzzyaaag64
real(dp) xyzzyaaah64
logical,allocatable :: xyzzyaaai64(:,:,:)
allocate(xyzzyaaai64(wfdet_norb,wfdet_norb,xyzzyaaak1),stat=xyzzyaaag6&
&4)
call check_alloc(xyzzyaaag64,"PREPARE_SPARSE","gmat_is_zero")
xyzzyaaai64=.false.
where(xyzzyaabf1==xyzzyaaav1.and.xyzzyaaal1==0.d0)
xyzzyaaai64=.true.
endwhere
do xyzzyaaad64=1,xyzzyaabi1
xyzzyaaac64=xyzzyaabk1(1,1,xyzzyaaad64)
xyzzyaaaa64=xyzzyaabk1(2,1,xyzzyaaad64)
xyzzyaaab64=xyzzyaabk1(3,1,xyzzyaaad64)
if(xyzzyaabf1(xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64)==xyzzyaaav1.and.xyz&
&zyaaal1(xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64)==0.d0)then
do xyzzyaaae64=2,xyzzyaabm1(xyzzyaaad64)
xyzzyaaac64=xyzzyaabk1(1,xyzzyaaae64,xyzzyaaad64)
xyzzyaaaa64=xyzzyaabk1(2,xyzzyaaae64,xyzzyaaad64)
xyzzyaaab64=xyzzyaabk1(3,xyzzyaaae64,xyzzyaaad64)
xyzzyaaai64(xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64)=.true.
enddo
endif
enddo
xyzzyaaah64=1.d0-dble(count(xyzzyaaai64))/dble(xyzzyaaak1*wfdet_norb**&
&2)
if(am_master)then
call wout("Percentage of non-zero geminal matrix elements: "//trim(i2s&
&(int(xyzzyaaah64*100)))//"%")
endif
if(xyzzyaaah64<third.or.xyzzyaaad1==xyzzyaaaa1)then
if(am_master)then
call wout("The geminal matrices will be represented by sparse matrices&
&.")
endif
xyzzyaaan1=.true.
allocate(xyzzyaaaq1(wfdet_norb,xyzzyaaak1),xyzzyaaar1(wfdet_norb,xyzzy&
&aaak1),stat=xyzzyaaag64)
call check_alloc(xyzzyaaag64,"PREPARE_SPARSE","gmat_i*_num")
forall(xyzzyaaac64=1:xyzzyaaak1,xyzzyaaaa64=1:wfdet_norb)
xyzzyaaaq1(xyzzyaaaa64,xyzzyaaac64)=count(.not.xyzzyaaai64(xyzzyaaaa64&
&,:,xyzzyaaac64))
endforall
forall(xyzzyaaac64=1:xyzzyaaak1,xyzzyaaab64=1:wfdet_norb)
xyzzyaaar1(xyzzyaaab64,xyzzyaaac64)=count(.not.xyzzyaaai64(:,xyzzyaaab&
&64,xyzzyaaac64))
endforall
xyzzyaaas1=maxval(xyzzyaaaq1)
xyzzyaaat1=maxval(xyzzyaaar1)
allocate(xyzzyaaao1(xyzzyaaas1,wfdet_norb,xyzzyaaak1),xyzzyaaap1(xyzzy&
&aaat1,wfdet_norb,xyzzyaaak1),stat=xyzzyaaag64)
call check_alloc(xyzzyaaag64,"PREPARE_SPARSE","gmat_i*")
do xyzzyaaac64=1,xyzzyaaak1
do xyzzyaaaa64=1,wfdet_norb
xyzzyaaaf64=0
do xyzzyaaab64=1,wfdet_norb
if(.not.xyzzyaaai64(xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64))then
xyzzyaaaf64=xyzzyaaaf64+1
xyzzyaaao1(xyzzyaaaf64,xyzzyaaaa64,xyzzyaaac64)=xyzzyaaab64
endif
enddo
enddo
do xyzzyaaab64=1,wfdet_norb
xyzzyaaaf64=0
do xyzzyaaaa64=1,wfdet_norb
if(.not.xyzzyaaai64(xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64))then
xyzzyaaaf64=xyzzyaaaf64+1
xyzzyaaap1(xyzzyaaaf64,xyzzyaaab64,xyzzyaaac64)=xyzzyaaaa64
endif
enddo
enddo
enddo
else
if(am_master)then
call wout("The geminal matrices will be represented by dense matrices.&
&")
endif
xyzzyaaan1=.false.
endif
deallocate(xyzzyaaai64)
end subroutine xyzzyaagf1
subroutine xyzzyaagg1
implicit none
integer xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65,xyzzyaaad65,xyzzyaaae65,xy&
&zzyaaaf65
if(.not.opt_geminal)return
allocate(xyzzyaabh1(3,xyzzyaabb1),stat=xyzzyaaaf65)
call check_alloc(xyzzyaaaf65,"BUILD_GMAT_PARAM_LIST","gmat_params")
xyzzyaaad65=0
do xyzzyaaac65=1,xyzzyaaak1
do xyzzyaaab65=1,wfdet_norb
do xyzzyaaaa65=1,wfdet_norb
if(xyzzyaabf1(xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65)/=xyzzyaaaw1)cycle
xyzzyaaad65=xyzzyaaad65+1
xyzzyaaae65=xyzzyaagj1(xyzzyaaac65,xyzzyaaaa65,xyzzyaaab65)
if(xyzzyaaae65==0)then
xyzzyaabh1(1,xyzzyaaad65)=xyzzyaaac65
xyzzyaabh1(2,xyzzyaaad65)=xyzzyaaaa65
xyzzyaabh1(3,xyzzyaaad65)=xyzzyaaab65
else
xyzzyaabh1(:,xyzzyaaad65)=xyzzyaabk1(:,1,xyzzyaaae65)
endif
enddo
enddo
enddo
end subroutine xyzzyaagg1
subroutine xyzzyaagh1
implicit none
integer xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66,xyzzyaaad66,xyzzyaaae66,xy&
&zzyaaaf66,xyzzyaaag66
if(.not.opt_geminal)return
allocate(xyzzyaabo1(xyzzyaaak1,xyzzyaabb1),stat=xyzzyaaaa66)
call check_alloc(xyzzyaaaa66,"BUILD_G_AFFECT_MAP","g_affect_map")
xyzzyaabo1=.false.
do xyzzyaaaf66=1,xyzzyaabb1
xyzzyaaab66=xyzzyaabh1(1,xyzzyaaaf66)
xyzzyaaac66=xyzzyaabh1(2,xyzzyaaaf66)
xyzzyaaad66=xyzzyaabh1(3,xyzzyaaaf66)
xyzzyaabo1(xyzzyaaab66,xyzzyaaaf66)=.true.
xyzzyaaae66=xyzzyaagj1(xyzzyaaab66,xyzzyaaac66,xyzzyaaad66)
if(xyzzyaaae66/=0)then
do xyzzyaaag66=1,xyzzyaabm1(xyzzyaaae66)
xyzzyaaab66=xyzzyaabk1(1,xyzzyaaag66,xyzzyaaae66)
xyzzyaabo1(xyzzyaaab66,xyzzyaaaf66)=.true.
enddo
endif
enddo
end subroutine xyzzyaagh1
subroutine xyzzyaagi1
implicit none
integer xyzzyaaaa67,xyzzyaaab67,xyzzyaaac67,xyzzyaaad67,xyzzyaaae67
real(dp) xyzzyaaaf67
do xyzzyaaaa67=1,xyzzyaabi1
xyzzyaaac67=xyzzyaabk1(1,1,xyzzyaaaa67)
xyzzyaaad67=xyzzyaabk1(2,1,xyzzyaaaa67)
xyzzyaaae67=xyzzyaabk1(3,1,xyzzyaaaa67)
xyzzyaaaf67=xyzzyaaal1(xyzzyaaad67,xyzzyaaae67,xyzzyaaac67)
do xyzzyaaab67=2,xyzzyaabm1(xyzzyaaaa67)
xyzzyaaac67=xyzzyaabk1(1,xyzzyaaab67,xyzzyaaaa67)
xyzzyaaad67=xyzzyaabk1(2,xyzzyaaab67,xyzzyaaaa67)
xyzzyaaae67=xyzzyaabk1(3,xyzzyaaab67,xyzzyaaaa67)
xyzzyaaal1(xyzzyaaad67,xyzzyaaae67,xyzzyaaac67)=xyzzyaaaf67
enddo
enddo
do xyzzyaaaa67=1,xyzzyaabj1
xyzzyaaaf67=xyzzyaaam1(xyzzyaabl1(1,xyzzyaaaa67))
do xyzzyaaab67=2,xyzzyaabn1(xyzzyaaaa67)
xyzzyaaam1(xyzzyaabl1(xyzzyaaab67,xyzzyaaaa67))=xyzzyaaaf67
enddo
enddo
end subroutine xyzzyaagi1
function xyzzyaagj1(gem,row,col) result(iconstraint)
implicit none
integer,intent(in) :: gem,row,col
integer iconstraint,xyzzyaaaa68
logical xyzzyaaab68
xyzzyaaab68=.false.
do iconstraint=1,xyzzyaabi1
do xyzzyaaaa68=1,xyzzyaabm1(iconstraint)
if(gem==xyzzyaabk1(1,xyzzyaaaa68,iconstraint).and. row==xyzzyaabk1(2,x&
&yzzyaaaa68,iconstraint).and. col==xyzzyaabk1(3,xyzzyaaaa68,iconstrain&
&t))then
xyzzyaaab68=.true.
return
endif
enddo
enddo
if(.not.xyzzyaaab68)iconstraint=0
end function xyzzyaagj1
subroutine xyzzyaagk1
implicit none
integer xyzzyaaaa69,xyzzyaaab69,xyzzyaaac69
character(casl_keysize) label
character(512) errmsg
call set_casl_item("Default g optimizability",trim(opt2s(xyzzyaaav1)),&
&errmsg)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
call set_casl_item("Default c optimizability",trim(opt2s(xyzzyaaav1)),&
&errmsg)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
do xyzzyaaaa69=1,xyzzyaaak1
call push_casl_context("Geminal "//trim(i2s(xyzzyaaaa69)))
call set_casl_block("Parameters",errmsg)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
call push_casl_context("Parameters")
call set_casl_block("c",errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
call set_casl_item("c:%u1",xyzzyaaam1(xyzzyaaaa69),errmsg)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
call set_casl_item("c:%u2",trim(opt2s(xyzzyaabg1(xyzzyaaaa69))),errmsg&
&)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
do xyzzyaaab69=1,wfdet_norb
do xyzzyaaac69=1,wfdet_norb
if(xyzzyaabf1(xyzzyaaab69,xyzzyaaac69,xyzzyaaaa69)==xyzzyaaax1.or.xyzz&
&yaabf1(xyzzyaaab69,xyzzyaaac69,xyzzyaaaa69)==xyzzyaaav1.and.xyzzyaaal&
&1(xyzzyaaab69,xyzzyaaac69,xyzzyaaaa69)==0.d0)cycle
label=inds2glabel(xyzzyaaab69,xyzzyaaac69)
call set_casl_block(trim(label),errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
call set_casl_item(label//":%u1",xyzzyaaal1(xyzzyaaab69,xyzzyaaac69,xy&
&zzyaaaa69),errmsg)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
call set_casl_item(label//":%u2",trim(opt2s(xyzzyaabf1(xyzzyaaab69,xyz&
&zyaaac69,xyzzyaaaa69))),errmsg)
if(len_trim(errmsg)>0)call errstop_master('SET_GEMINAL_CASL',trim(errm&
&sg))
enddo
enddo
call pop_casl_context()
call pop_casl_context()
enddo
end subroutine xyzzyaagk1
subroutine update_geminal_casl
implicit none
integer xyzzyaaaa70,xyzzyaaab70,xyzzyaaac70
character(512) errmsg
call push_casl_context('parameters.casl:GEMINAL')
do xyzzyaaaa70=1,xyzzyaaak1
call push_casl_context("Geminal "//trim(i2s(xyzzyaaaa70)//":Parameters&
&"))
if(xyzzyaabg1(xyzzyaaaa70)==xyzzyaaaw1.or.xyzzyaabg1(xyzzyaaaa70)==xyz&
&zyaaax1)then
call set_casl_item("c:%u1",xyzzyaaam1(xyzzyaaaa70),errmsg)
if(len_trim(errmsg)>0)call errstop_master('UPDATE_GEMINAL_CASL',trim(e&
&rrmsg))
endif
do xyzzyaaab70=1,wfdet_norb
do xyzzyaaac70=1,wfdet_norb
if(xyzzyaabf1(xyzzyaaab70,xyzzyaaac70,xyzzyaaaa70)==xyzzyaaav1.or.xyzz&
&yaabf1(xyzzyaaab70,xyzzyaaac70,xyzzyaaaa70)==xyzzyaaax1)cycle
call set_casl_item(inds2glabel(xyzzyaaab70,xyzzyaaac70)//":%u1",xyzzya&
&aal1(xyzzyaaab70,xyzzyaaac70,xyzzyaaaa70),errmsg)
if(len_trim(errmsg)>0)call errstop_master('UPDATE_GEMINAL_CASL',trim(e&
&rrmsg))
enddo
enddo
call pop_casl_context()
enddo
call pop_casl_context()
end subroutine update_geminal_casl
subroutine clear_scratch_geminal(is)
implicit none
integer,intent(in) :: is
xyzzyaabq1(:,:,is)=.false.
xyzzyaabs1(:,:,is)=.false.
xyzzyaabu1(:,:,is)=.false.
xyzzyaaby1(:,:,is)=.false.
xyzzyaaca1(:,:,is)=.false.
xyzzyaade1(:,is)=.false.
xyzzyaadg1(:,is)=.false.
xyzzyaadu1(:,:,is)=.false.
if(.not.use_backflow)xyzzyaadw1(:,:,is)=.false.
xyzzyaadj1(is)=.false.
xyzzyaadq1(:,is)=.false.
xyzzyaads1(:,is)=.false.
xyzzyaacc1(:,is)=.false.
xyzzyaacg1(:,is)=.false.
xyzzyaack1(:,is)=.false.
xyzzyaaco1(:,is)=.false.
xyzzyaacs1(:,is)=.false.
xyzzyaace1(:,is)=.false.
xyzzyaaci1(:,is)=.false.
xyzzyaacm1(:,is)=.false.
xyzzyaacq1(:,is)=.false.
xyzzyaacu1(:,is)=.false.
xyzzyaadb1(:,is)=.false.
xyzzyaada1(:,is)=0
if(use_backflow)then
xyzzyaabw1(:,:,is)=.false.
xyzzyaacy1(:,is)=.false.
xyzzyaady1(:,is)=.false.
xyzzyaadm1(:,is)=.false.
xyzzyaado1(is)=.false.
endif
if(use_backflow)call clear_scratch_bf(is)
end subroutine clear_scratch_geminal
subroutine gen_config_geminal(pt_config)
implicit none
type(config_wfn_geminal),pointer :: pt_config
integer xyzzyaaaa72
allocate(pt_config,stat=xyzzyaaaa72)
call check_alloc(xyzzyaaaa72,'GEN_CONFIG_GEMINAL','container')
end subroutine gen_config_geminal
subroutine delete_config_geminal(pt_config)
implicit none
type(config_wfn_geminal),pointer :: pt_config
deallocate(pt_config)
end subroutine delete_config_geminal
subroutine copy_config_geminal(pt_from,pt_to)
implicit none
type(config_wfn_geminal),pointer :: pt_from,pt_to
end subroutine copy_config_geminal
subroutine config_to_pt_geminal(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_geminal),pointer :: pt_config
end subroutine config_to_pt_geminal
subroutine pt_to_config_geminal(pt_config)
implicit none
type(config_wfn_geminal),pointer :: pt_config
end subroutine pt_to_config_geminal
subroutine redist_allocations_geminal(kmax)
implicit none
integer,intent(in) :: kmax
end subroutine redist_allocations_geminal
subroutine redist_load_geminal(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_geminal),pointer :: pt_config
end subroutine redist_load_geminal
subroutine redist_send_geminal(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
end subroutine redist_send_geminal
subroutine redist_recv_geminal(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
end subroutine redist_recv_geminal
subroutine redist_save_geminal(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_geminal),pointer :: pt_config
end subroutine redist_save_geminal
subroutine redist_deallocations_geminal
implicit none
end subroutine redist_deallocations_geminal
subroutine load_from_pt_geminal(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_geminal),pointer :: pt_config
end subroutine load_from_pt_geminal
subroutine save_to_pt_geminal(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_geminal),pointer :: pt_config
end subroutine save_to_pt_geminal
subroutine clone_scratch_geminal(is,js)
implicit none
integer,intent(in) :: is,js
end subroutine clone_scratch_geminal
subroutine add_config_geminal_items(is)
implicit none
integer,intent(in) :: is
end subroutine add_config_geminal_items
subroutine setup_storage_geminal(nconfig,ignore)
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(xyzzyaabe1)
end subroutine setup_storage_geminal
subroutine finish_storage_geminal
implicit none
end subroutine finish_storage_geminal
subroutine load_from_storage_geminal(is,icfg)
implicit none
integer,intent(in) :: is,icfg
end subroutine load_from_storage_geminal
subroutine save_to_storage_geminal(is,icfg)
implicit none
integer,intent(in) :: is,icfg
end subroutine save_to_storage_geminal
subroutine enumerate_plot_geminal(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
n=0
end subroutine enumerate_plot_geminal
subroutine query_plot_geminal(iplot,ii,rank,is_complex,has_stderr,rot_&
&tensor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
nfunctions=0
end subroutine query_plot_geminal
subroutine get_plot_geminal(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
f(1)=0.d0
end subroutine get_plot_geminal
subroutine finish_plot_geminal
implicit none
end subroutine finish_plot_geminal
subroutine xyzzyaagl1(n,logdet,gemmat,geminv,fpeinfo)
implicit none
integer,intent(in) :: n
integer,intent(inout) :: fpeinfo
real(dp),intent(inout) :: gemmat(nemax,nemax,real1_complex2),geminv(ne&
&max,nemax,real1_complex2)
complex(dp),intent(inout) :: logdet
integer xyzzyaaaa95
logical xyzzyaaab95
logdet=czero
fpeinfo=xyzzyaaag1
if(complex_wf)then
call lu_decom_cmplx_dz(gemmat,xyzzyaaeh1,xyzzyaadz1,nemax,n,xyzzyaaab9&
&5)
if(.not.xyzzyaaab95)then
logdet=lu_logdet_cmplx_dz(gemmat,xyzzyaadz1,nemax,n)
if(dble(logdet)>xyzzyaaaf1)then
geminv=0.d0
do xyzzyaaaa95=1,n
geminv(xyzzyaaaa95,xyzzyaaaa95,1)=1.d0
enddo
call lu_solve_n_cmplx_dz(gemmat,xyzzyaaeh1,xyzzyaadz1,geminv,xyzzyaaei&
&1,nemax,n,nemax)
else
xyzzyaaab95=.true.
endif
endif
if(xyzzyaaab95)then
fpeinfo=xyzzyaaai1
endif
else
call lu_decom(gemmat,xyzzyaadz1,nemax,n,xyzzyaaab95)
if(.not.xyzzyaaab95)then
logdet=lu_logdet(gemmat,xyzzyaadz1,nemax,n)
if(dble(logdet)>xyzzyaaaf1)then
geminv=0.d0
do xyzzyaaaa95=1,n
geminv(xyzzyaaaa95,xyzzyaaaa95,1)=1.d0
enddo
call lu_solve_n(gemmat,xyzzyaadz1,geminv,nemax,n,nemax)
else
xyzzyaaab95=.true.
endif
endif
if(xyzzyaaab95)then
fpeinfo=xyzzyaaai1
endif
endif
end subroutine xyzzyaagl1
subroutine xyzzyaagm1(ie,ispin,xyzzyaacz1,gem_matrix1_ch,det_ratio,fpe&
&info)
implicit none
integer,intent(in) :: ie,ispin
integer,intent(inout) :: fpeinfo
real(dp),intent(in) :: xyzzyaacz1(nemax,nemax,real1_complex2),gem_matr&
&ix1_ch(nemax,real1_complex2)
complex(dp),intent(out) :: det_ratio
real(dp) xyzzyaaaa96
fpeinfo=xyzzyaaag1
if(complex_wf)then
if(ispin==1)then
det_ratio=zdotu_cc(nemax,gem_matrix1_ch(1,re),gem_matrix1_ch(1,im),1,x&
&yzzyaacz1(1,ie,re),xyzzyaacz1(1,ie,im),1)
else
det_ratio=zdotu_cc(nemax,gem_matrix1_ch(1,re),gem_matrix1_ch(1,im),1,x&
&yzzyaacz1(ie,1,re),xyzzyaacz1(ie,1,im),nemax)
endif
if(det_ratio==czero)then
if(all(gem_matrix1_ch(1:nemax,1:2)==0.d0))then
fpeinfo=xyzzyaaai1
else
fpeinfo=xyzzyaaai1
endif
endif
else
select case(ispin)
case(1)
xyzzyaaaa96=ddot(nemax,gem_matrix1_ch(1,1),1,xyzzyaacz1(1,ie,1),1)
case(2)
xyzzyaaaa96=ddot(nemax,gem_matrix1_ch(1,1),1,xyzzyaacz1(ie,1,1),nemax)
endselect
if(xyzzyaaaa96==0.d0)then
if(all(gem_matrix1_ch(1:nemax,1)==0.d0))then
fpeinfo=xyzzyaaai1
else
fpeinfo=xyzzyaaai1
endif
endif
det_ratio=cmplx(xyzzyaaaa96,0d0,dp)
endif
end subroutine xyzzyaagm1
subroutine xyzzyaagn1(ie,ispin,matrix1_ch,xyzzyaacz1,det_ratio)
implicit none
integer,intent(in) :: ie,ispin
real(dp),intent(in) :: matrix1_ch(nemax,real1_complex2)
real(dp),intent(inout) :: xyzzyaacz1(nemax,nemax,real1_complex2)
complex(dp),intent(in) :: det_ratio
integer xyzzyaaaa97,xyzzyaaab97
real(dp) xyzzyaaac97,xyzzyaaad97,xyzzyaaae97
complex(dp) xyzzyaaaf97,xyzzyaaag97
xyzzyaaab97=nemax
if(complex_wf)then
xyzzyaaaf97=c_one/det_ratio
xyzzyaaac97=dble(xyzzyaaaf97)
xyzzyaaad97=aimag(xyzzyaaaf97)
select case(ispin)
case(1)
do xyzzyaaaa97=1,xyzzyaaab97
if(xyzzyaaaa97==ie)cycle
xyzzyaaag97=-xyzzyaaaf97*zdotu_cc(xyzzyaaab97,matrix1_ch(1,re),matrix1&
&_ch(1,im),1,xyzzyaacz1(1,xyzzyaaaa97,re),xyzzyaacz1(1,xyzzyaaaa97,im)&
&,1)
call zaxpy_cc(xyzzyaaab97,dble(xyzzyaaag97),aimag(xyzzyaaag97),xyzzyaa&
&cz1(1,ie,re),xyzzyaacz1(1,ie,im),1,xyzzyaacz1(1,xyzzyaaaa97,re),xyzzy&
&aacz1(1,xyzzyaaaa97,im),1)
enddo
call zscal_cc(xyzzyaaab97,xyzzyaaac97,xyzzyaaad97,xyzzyaacz1(1,ie,re),&
&xyzzyaacz1(1,ie,im),1)
case(2)
do xyzzyaaaa97=1,xyzzyaaab97
if(xyzzyaaaa97==ie)cycle
xyzzyaaag97=-xyzzyaaaf97*zdotu_cc(xyzzyaaab97,xyzzyaacz1(xyzzyaaaa97,1&
&,re),xyzzyaacz1(xyzzyaaaa97,1,im),xyzzyaaab97,matrix1_ch(1,re),matrix&
&1_ch(1,im),1)
call zaxpy_cc(xyzzyaaab97,dble(xyzzyaaag97),aimag(xyzzyaaag97),xyzzyaa&
&cz1(ie,1,re),xyzzyaacz1(ie,1,im),xyzzyaaab97,xyzzyaacz1(xyzzyaaaa97,1&
&,1),xyzzyaacz1(xyzzyaaaa97,1,im),xyzzyaaab97)
enddo
call zscal_cc(xyzzyaaab97,xyzzyaaac97,xyzzyaaad97,xyzzyaacz1(ie,1,re),&
&xyzzyaacz1(ie,1,im),xyzzyaaab97)
endselect
else
xyzzyaaac97=1.d0/real(det_ratio,dp)
select case(ispin)
case(1)
do xyzzyaaaa97=1,xyzzyaaab97
if(xyzzyaaaa97==ie)cycle
xyzzyaaae97=-xyzzyaaac97*ddot(xyzzyaaab97,matrix1_ch,1,xyzzyaacz1(1,xy&
&zzyaaaa97,1),1)
call daxpy(xyzzyaaab97,xyzzyaaae97,xyzzyaacz1(1,ie,1),1,xyzzyaacz1(1,x&
&yzzyaaaa97,1),1)
enddo
call dscal(xyzzyaaab97,xyzzyaaac97,xyzzyaacz1(1,ie,1),1)
case(2)
do xyzzyaaaa97=1,xyzzyaaab97
if(xyzzyaaaa97==ie)cycle
xyzzyaaae97=-xyzzyaaac97*ddot(xyzzyaaab97,xyzzyaacz1(xyzzyaaaa97,1,1),&
&xyzzyaaab97,matrix1_ch,1)
call daxpy(xyzzyaaab97,xyzzyaaae97,xyzzyaacz1(ie,1,1),xyzzyaaab97,xyzz&
&yaacz1(xyzzyaaaa97,1,1),xyzzyaaab97)
enddo
call dscal(xyzzyaaab97,xyzzyaaac97,xyzzyaacz1(ie,1,1),xyzzyaaab97)
endselect
endif
end subroutine xyzzyaagn1
subroutine xyzzyaago1(ie,ispin,mat_size,mat_in,mat_ch,mat_out)
implicit none
integer,intent(in) :: ie,ispin,mat_size
real(dp),intent(in) :: mat_in(mat_size,mat_size,real1_complex2),mat_ch&
&(mat_size,real1_complex2)
real(dp),intent(inout) :: mat_out(mat_size,mat_size,real1_complex2)
integer xyzzyaaaa98
if(ispin==1)then
forall(xyzzyaaaa98=1:mat_size,xyzzyaaaa98/=ie)
mat_out(xyzzyaaaa98,:,:)=mat_in(xyzzyaaaa98,:,:)
endforall
mat_out(ie,:,:)=mat_ch
else
forall(xyzzyaaaa98=1:mat_size,xyzzyaaaa98/=ie)
mat_out(:,xyzzyaaaa98,:)=mat_in(:,xyzzyaaaa98,:)
endforall
mat_out(:,ie,:)=mat_ch
endif
end subroutine xyzzyaago1
subroutine xyzzyaagp1(ie,ispin,mat_size,mat_in,mat_ch,mat_out)
implicit none
integer,intent(in) :: ie,ispin,mat_size
real(dp),intent(in) :: mat_in(3,mat_size,mat_size,real1_complex2),mat_&
&ch(3,mat_size,real1_complex2)
real(dp),intent(inout) :: mat_out(3,mat_size,mat_size,real1_complex2)
integer xyzzyaaaa99
if(ispin==1)then
forall(xyzzyaaaa99=1:mat_size,xyzzyaaaa99/=ie)
mat_out(:,xyzzyaaaa99,:,:)=mat_in(:,xyzzyaaaa99,:,:)
endforall
mat_out(:,ie,:,:)=mat_ch
else
forall(xyzzyaaaa99=1:mat_size,xyzzyaaaa99/=ie)
mat_out(:,:,xyzzyaaaa99,:)=mat_in(:,:,xyzzyaaaa99,:)
endforall
mat_out(:,:,ie,:)=mat_ch
endif
end subroutine xyzzyaagp1
subroutine xyzzyaagq1(g,g_irow,g_irow_num,xvec_re,xvec_im,xstep,gx_re,&
&gx_im,gxstep)
implicit none
real(dp),intent(in) :: g(wfdet_norb,wfdet_norb),xvec_re(*),xvec_im(*)
integer,intent(in) :: g_irow(xyzzyaaas1,wfdet_norb),g_irow_num(wfdet_n&
&orb),xstep,gxstep
real(dp),intent(inout) :: gx_re(*),gx_im(*)
if(.not.xyzzyaaan1)then
call dgemv('N',wfdet_norb,wfdet_norb,1.d0,g,wfdet_norb,xvec_re(1),xste&
&p,0.d0,gx_re(1),1)
if(complex_wf)then
call dgemv('N',wfdet_norb,wfdet_norb,1.d0,g,wfdet_norb,xvec_im(1),xste&
&p,0.d0,gx_im(1),gxstep)
endif
else
if(complex_wf)then
call d2_sparse_matrix_vector(wfdet_norb,1.d0,g,wfdet_norb,g_irow,xyzzy&
&aaas1,g_irow_num,xvec_re(1),xvec_im(1),xstep,gx_re(1),gx_im(1),gxstep&
&)
else
call d_sparse_matrix_vector(wfdet_norb,1.d0,g,wfdet_norb,g_irow,xyzzya&
&aas1,g_irow_num,xvec_re(1),xstep,gx_re(1),gxstep)
endif
endif
end subroutine xyzzyaagq1
subroutine xyzzyaagr1(xyzzyaaee1,fpeinfo,gem_inv_ie,figem)
implicit none
integer,intent(in) :: fpeinfo
real(dp),intent(in) :: xyzzyaaee1(3,nemax,real1_complex2),gem_inv_ie(n&
&emax,real1_complex2)
complex(dp),intent(inout) :: figem(3)
if(fpeinfo/=xyzzyaaag1)return
xyzzyaaeg1(1:3,:)=0.d0
call d_or_z_gemv(complex_wf,'N',3,nemax,1.d0,xyzzyaaee1(:,:,re),xyzzya&
&aee1(:,:,im),3,gem_inv_ie(:,re),gem_inv_ie(:,im),1,0.d0,xyzzyaaeg1(:,&
&re),xyzzyaaeg1(:,im),1)
figem(1:3)=cmplx(xyzzyaaeg1(1:3,re),xyzzyaaeg1(1:3,im),dp)
end subroutine xyzzyaagr1
subroutine xyzzyaags1(xyzzyaaef1,loggem,fpeinfo,gem_inv_ie,lapgem)
implicit none
integer,intent(in) :: fpeinfo
real(dp),intent(in) :: xyzzyaaef1(nemax,real1_complex2),gem_inv_ie(nem&
&ax,real1_complex2)
complex(dp),intent(in) :: loggem
complex(dp),intent(out) :: lapgem
real(dp) xyzzyaaaa102
if(fpeinfo/=xyzzyaaag1)return
if(complex_wf)then
lapgem=zdotu_cc(nemax,xyzzyaaef1(1,re),xyzzyaaef1(1,im),1,gem_inv_ie(1&
&,re),gem_inv_ie(1,im),1)
lapgem=exp(loggem)*lapgem
else
xyzzyaaaa102=ddot(nemax,xyzzyaaef1(1,1),1,gem_inv_ie(1,1),1)
lapgem=exp(loggem)*cmplx(xyzzyaaaa102,0.d0,dp)
endif
end subroutine xyzzyaags1
subroutine xyzzyaagt1(rele,sele,fill_val,fill_fsd,orbval,orbgrad,orbla&
&p,orbsderivs)
integer,intent(in) :: sele
real(dp),intent(in) :: rele(3)
real(dp),intent(inout) :: orbval(wfdet_norb,real1_complex2),orbgrad(3,&
&wfdet_norb,real1_complex2),orblap(wfdet_norb,real1_complex2)
real(dp),intent(inout), optional :: orbsderivs(6,wfdet_norb,real1_comp&
&lex2)
logical,intent(in) :: fill_val,fill_fsd
if(xyzzyaaae1)then
if(present(orbsderivs).and..not.fill_fsd)then
call errstop_master("FILL_ORBS", "Unexpected Behaviour.")
endif
endif
call wfdet(rele,sele,sele,wfdet_norb,xyzzyaaea1,fill_val,fill_fsd,orbv&
&al,orbgrad,orblap,xyzzyaaec1,orbsderivs)
end subroutine xyzzyaagt1
pure integer function count_char(s,g)
implicit none
character(*),intent(in) :: s
character,intent(in) :: g
integer pos,new_pos
count_char=0
pos=0
do
new_pos=scan(s(pos+1:),g)+pos
if(new_pos/=pos)then
count_char=count_char+1
else
return
endif
pos=new_pos
enddo
end function count_char
character(casl_keysize) function inds2glabel(row,col)result(glabel)
implicit none
integer,intent(in) :: row,col
glabel="g_"//trim(i2s(row))//","//trim(i2s(col))
end function inds2glabel
character(casl_keysize) function inds2extglabel(gem,row,col)result(gla&
&bel)
implicit none
integer,intent(in) :: gem,row,col
glabel=trim(i2s(gem))//"^"//trim(inds2glabel(row,col))
end function inds2extglabel
character(casl_valsize) function opt2s(opt_code)
implicit none
integer,intent(in) :: opt_code
select case(opt_code)
case(xyzzyaaau1)
opt2s="undefined"
case(xyzzyaaav1)
opt2s="fixed"
case(xyzzyaaaw1)
opt2s="optimizable"
case(xyzzyaaax1)
opt2s="determined"
case default
call errstop_master("OPT2STRING","Unidentified opt_code. Received "//t&
&rim(i2s(opt_code))//". Bug.")
end select
end function opt2s
integer function xyzzyaagu1(s)
implicit none
character(*),intent(in) :: s
select case(trim(switch_case(s,to_lower=.true.)))
case("fixed")
xyzzyaagu1=xyzzyaaav1
case("fix")
xyzzyaagu1=xyzzyaaav1
case("optimizable")
xyzzyaagu1=xyzzyaaaw1
case("opt")
xyzzyaagu1=xyzzyaaaw1
case("determined")
xyzzyaagu1=xyzzyaaax1
case default
call errstop_master("STRING2OPT","Optimizability value '"//trim(s)//"'&
& is invalid. Allowed values are: fix[ed],opt[imizable],determined")
endselect
end function xyzzyaagu1
integer function xyzzyaagv1(string_buffer) result(opt_code)
implicit none
character(casl_valsize),intent(in) :: string_buffer
opt_code=xyzzyaagu1(string_buffer)
if(opt_code==xyzzyaaax1)call errstop_master("PARSE_DEFAULT_OPT", "The &
&default optimizability should be either fixed or optimizable.")
end function xyzzyaagv1
subroutine xyzzyaagw1(label,row,col)
implicit none
integer,intent(out) :: row,col
character(*),intent(in) :: label
integer xyzzyaaaa109,xyzzyaaab109
character(len(label)) buffer
xyzzyaaaa109=index(label,"_")
xyzzyaaab109=index(label,",")
if(xyzzyaaaa109==0.or.xyzzyaaab109==0)call errstop_master("PARSE_GLABE&
&L","Invalid syntax in "//trim(label)//".")
buffer=label(xyzzyaaaa109+1:xyzzyaaab109)
read(buffer,*)row
buffer=trim(label(xyzzyaaab109+1:))
read(buffer,*)col
end subroutine xyzzyaagw1
subroutine xyzzyaagx1(label,gem,row,col)
integer,intent(out) :: gem,row,col
character(*),intent(in) :: label
integer xyzzyaaaa110
character(len(label)) buffer
xyzzyaaaa110=scan(label,"^")
if(xyzzyaaaa110==0)call errstop_master("PARSE_EXTENDED_GLABEL","Invali&
&d syntax in '"//trim(label)//"'. The correct syntax is geminal^g_row,&
&column")
buffer=trim(label(:xyzzyaaaa110-1))
read(buffer,*)gem
call xyzzyaagw1(trim(label(xyzzyaaaa110+1:)),row,col)
end subroutine xyzzyaagx1
subroutine xyzzyaagy1(label,gem)
implicit none
integer,intent(out) :: gem
character(*),intent(in) :: label
integer xyzzyaaaa111
character(len(label)) buffer
xyzzyaaaa111=scan(label,"^")
buffer=trim(label(:xyzzyaaaa111-1))
read(buffer,*)gem
end subroutine xyzzyaagy1
subroutine xyzzyaagz1(buffer,istart1,iend1,istart2,iend2)
implicit none
integer,intent(out) :: istart1,iend1,istart2,iend2
character(*),intent(in) :: buffer
integer xyzzyaaaa112,xyzzyaaab112,xyzzyaaac112
logical xyzzyaaad112
xyzzyaaaa112=index(buffer,"(")
xyzzyaaab112=index(buffer,",")
xyzzyaaac112=index(buffer,")")
if(xyzzyaaaa112==0.or.xyzzyaaab112==0.or.xyzzyaaac112==0.or.xyzzyaaaa1&
&12>xyzzyaaab112.or.xyzzyaaab112>xyzzyaaac112)call errstop_master('PAR&
&SE_MATRIX_RANGE','Syntax error in matrix range.')
call xyzzyaaha1(buffer(xyzzyaaaa112+1:xyzzyaaab112-1),istart1,iend1,xy&
&zzyaaad112)
call xyzzyaaha1(buffer(xyzzyaaab112+1:xyzzyaaac112-1),istart2,iend2,xy&
&zzyaaad112)
end subroutine xyzzyaagz1
subroutine xyzzyaaha1(buffer,istart,iend,has_end)
implicit none
integer,intent(out) :: istart,iend
logical,intent(out) :: has_end
character(*),intent(in) :: buffer
integer xyzzyaaaa113,xyzzyaaab113,xyzzyaaac113
character(len(buffer)) :: xyzzyaaad113
xyzzyaaab113=index(buffer,":")
xyzzyaaaa113=index(buffer,"(")
xyzzyaaac113=index(buffer,")")
if(xyzzyaaac113==0)xyzzyaaac113=len(buffer)+1
if(xyzzyaaab113==0)then
has_end=.false.
xyzzyaaad113=buffer(xyzzyaaaa113+1:xyzzyaaac113-1)
read(xyzzyaaad113,*)istart
iend=istart
else
has_end=.true.
xyzzyaaad113=buffer(xyzzyaaaa113+1:xyzzyaaab113-1)
read(xyzzyaaad113,*)istart
xyzzyaaad113=buffer(xyzzyaaab113+1:xyzzyaaac113-1)
read(xyzzyaaad113,*)iend
endif
end subroutine xyzzyaaha1
end module slaarnabe

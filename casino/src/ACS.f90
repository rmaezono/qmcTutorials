module slaarnacs
use dsp
use slaarnaao
use slaarnabc
use slaarnabe
use slaarnabl
use slaarnabw
use slaarnabn
use slaarnach
use slaarnacj
use store
use slaarnaad,   only : use_gbackflow,bf_assess_check_kinetic
use slaarnaag,  only : c_one,czero
use run_control,only : timer,check_alloc,errstop,errstop_master
implicit none
private
public setup_wfn_utils,finish_wfn_utils,define_config,define_config_on&
&eelec,define_config_twoelec,accept_move,define_config_oneion,wfn_rati&
&o,wfn_logval,wfn_loggrad,wfn_loglap,prefetch_wfn,add_config_wfn_items&
&,setup_wfn_params,finish_wfn_params,get_params,put_params,put_param1,&
&which_param,any_numerical_deriv,any_analytical_deriv,param_has_aderiv&
&,clone_scratch,clear_scratch_wfn,empty_scratch_wfn,invalidate_params,&
&invalidate_param1,wfn_aderiv,setup_storage,finish_storage,load_from_s&
&torage,save_to_storage,init_wfn,update_wfn_casl,wfn_assess_check_kine&
&tic,enumerate_plot_wfn,query_plot_wfn,get_plot_wfn,finish_plot_wfn
public gen_config,delete_config,copy_config,config_to_pt,pt_to_config,&
&redist_allocations,redist_load,redist_send,redist_recv,redist_save,re&
&dist_deallocations,load_from_pt,save_to_pt
public get_wfn_rmax
public config_wfn
public nsampling_levels,sampling_level_name,etot_scr,keimag_scr,etot_v&
&alid,ecomps_scr,ke_valid,nl_valid,erest_valid,grid_angles_scr,grid_an&
&gles_valid,etot_isnan_scr,etot_isinf_scr,ke_isnan_scr,ke_isinf_scr,nl&
&_isnan_scr,nl_isinf_scr,erest_isinf_scr
public psi_s
integer nsampling_levels
integer,parameter :: xyzzyaaaa1=1,xyzzyaaab1=2
integer,allocatable :: xyzzyaaac1(:),xyzzyaaad1(:)
character(80),allocatable :: sampling_level_name(:)
integer xyzzyaaae1
logical xyzzyaaaf1
character(20) psi_s
integer,parameter :: xyzzyaaag1=2
integer xyzzyaaah1,xyzzyaaai1(xyzzyaaag1)
integer,allocatable :: xyzzyaaaj1(:)
real(dp),allocatable :: xyzzyaaak1(:),xyzzyaaal1(:),xyzzyaaam1(:)
logical xyzzyaaan1
logical,allocatable :: xyzzyaaao1(:),xyzzyaaap1(:),xyzzyaaaq1(:),xyzzy&
&aaar1(:,:)
real(dp),allocatable :: etot_scr(:),keimag_scr(:),ecomps_scr(:,:),grid&
&_angles_scr(:,:,:,:)
logical,allocatable :: etot_isnan_scr(:),etot_isinf_scr(:),ke_isnan_sc&
&r(:),ke_isinf_scr(:),nl_isnan_scr(:),nl_isinf_scr(:),erest_isinf_scr(&
&:),etot_valid(:),ke_valid(:),nl_valid(:),erest_valid(:),grid_angles_v&
&alid(:)
type config_wfn
private
type(config_wfn_slater),pointer :: pt_slater
type(config_wfn_exmol),pointer :: pt_exmol
type(config_wfn_pfaff),pointer :: pt_pfaff
type(config_wfn_mahan),pointer :: pt_mahan
type(config_wfn_geminal),pointer :: pt_geminal
type(config_wfn_jastrow),pointer :: pt_jastrow
end type config_wfn
integer,parameter :: nthings_to_plot=1,xyzzyaaas1=1
integer :: xyzzyaaat1(nthings_to_plot)=0,xyzzyaaau1(0:2)=0
integer,allocatable :: xyzzyaaav1(:)
character(64),parameter :: plot_keyword(nthings_to_plot)=(/'wfn'/),plo&
&t_description(nthings_to_plot)=(/'wave function'/)
contains
subroutine init_wfn
use slaarnabt, only : swap1
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2
integer,allocatable :: xyzzyaaag2(:)
xyzzyaaaa2=0
xyzzyaaab2=0
select case(trim(psi_s))
case('none')
xyzzyaaae1=0
xyzzyaaaa2=0
case('slater')
xyzzyaaae1=1
call query_slater_levels(xyzzyaaaa2)
case('exmol')
xyzzyaaae1=2
call query_exmol_levels(xyzzyaaaa2)
case('pfaffian')
xyzzyaaae1=3
call query_pfaff_levels(xyzzyaaaa2)
case('geminal')
xyzzyaaae1=4
call query_geminal_levels(xyzzyaaaa2)
case('mahan_ex')
xyzzyaaae1=6
call query_mahan_levels(xyzzyaaaa2)
end select
if(use_jastrow)then
call query_jastrow_levels(xyzzyaaab2)
endif
nsampling_levels=xyzzyaaaa2+xyzzyaaab2
if(nsampling_levels==0)call errstop_master('INIT_WFN','The number of s&
&ampling levels is zero.  Please make sure that you actually have a tr&
&ial wave function!  E.g., check the value of psi_s.')
if(allocated(xyzzyaaag2))deallocate(xyzzyaaag2)
if(allocated(xyzzyaaac1))deallocate(xyzzyaaac1)
if(allocated(xyzzyaaad1))deallocate(xyzzyaaad1)
if(allocated(sampling_level_name))deallocate(sampling_level_name)
allocate(xyzzyaaag2(nsampling_levels),xyzzyaaac1(0:nsampling_levels),x&
&yzzyaaad1(0:nsampling_levels),sampling_level_name(0:nsampling_levels)&
&,stat=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'INIT_WFN','level arrays')
select case(xyzzyaaae1)
case(1)
call query_slater_level_details(xyzzyaaag2(1:xyzzyaaaa2),sampling_leve&
&l_name(1:xyzzyaaaa2))
case(2)
call query_exmol_level_details(xyzzyaaag2(1:xyzzyaaaa2),sampling_level&
&_name(1:xyzzyaaaa2))
case(3)
call query_pfaff_level_details(xyzzyaaag2(1:xyzzyaaaa2),sampling_level&
&_name(1:xyzzyaaaa2))
case(4)
call query_geminal_level_details(xyzzyaaag2(1:xyzzyaaaa2),sampling_lev&
&el_name(1:xyzzyaaaa2))
case(6)
call query_mahan_level_details(xyzzyaaag2(1:xyzzyaaaa2),sampling_level&
&_name(1:xyzzyaaaa2))
end select
if(use_jastrow)then
call query_jastrow_level_details(xyzzyaaag2(xyzzyaaaa2+1:nsampling_lev&
&els),sampling_level_name(xyzzyaaaa2+1:nsampling_levels))
endif
xyzzyaaac1=0
xyzzyaaad1=0
if(xyzzyaaae1/=0)then
xyzzyaaac1(1:xyzzyaaaa2)=xyzzyaaaa1
xyzzyaaad1(1:xyzzyaaaa2)=0
endif
if(use_jastrow)then
xyzzyaaac1(xyzzyaaaa2+1:xyzzyaaaa2+xyzzyaaab2)=xyzzyaaab1
xyzzyaaad1(xyzzyaaaa2+1:xyzzyaaaa2+xyzzyaaab2)=0
endif
do xyzzyaaac2=1,nsampling_levels
xyzzyaaad2=xyzzyaaac2
do xyzzyaaae2=xyzzyaaac2+1,nsampling_levels
if(xyzzyaaag2(xyzzyaaae2)>xyzzyaaag2(xyzzyaaad2))xyzzyaaad2=xyzzyaaae2
enddo
if(xyzzyaaad2>xyzzyaaac2)then
call swap1(xyzzyaaac1(xyzzyaaac2),xyzzyaaac1(xyzzyaaad2))
call swap1(xyzzyaaad1(xyzzyaaac2),xyzzyaaad1(xyzzyaaad2))
call swap1(sampling_level_name(xyzzyaaac2),sampling_level_name(xyzzyaa&
&ad2))
call swap1(xyzzyaaag2(xyzzyaaac2),xyzzyaaag2(xyzzyaaad2))
endif
enddo
deallocate(xyzzyaaag2)
xyzzyaaaf1=(use_jastrow.and.use_gjastrow).or.(use_backflow.and.use_gba&
&ckflow)
end subroutine init_wfn
subroutine setup_wfn_utils
implicit none
integer xyzzyaaaa3
select case(xyzzyaaae1)
case(1)
call setup_slater
case(2)
call setup_exmol
case(3)
call setup_pfaff
case(4)
call setup_geminal
case(6)
call setup_mahan
end select
if(use_jastrow)call setup_jastrow
if(xyzzyaaaf1)call setup_gbasis
do xyzzyaaaa3=1,nscratch
call clear_scratch_wfn(xyzzyaaaa3)
enddo
end subroutine setup_wfn_utils
subroutine finish_wfn_utils
implicit none
select case(xyzzyaaae1)
case(1)
call finish_slater
case(2)
call finish_exmol
case(3)
call finish_pfaff
case(4)
call finish_geminal
case(6)
call finish_mahan
end select
if(use_jastrow)call finish_jastrow
if(xyzzyaaaf1)call finish_gbasis
end subroutine finish_wfn_utils
subroutine define_config(is,rele,sele)
implicit none
integer,intent(in) :: is,sele(netot)
real(dp),intent(in) :: rele(3,netot)
call timer('DEFINE_CONFIG',.true.)
call clear_scratch_wfn(is)
call xyzzyaaax1(is)
call define_config_scratch(is,rele,sele)
call timer('DEFINE_CONFIG',.false.)
end subroutine define_config
subroutine define_config_oneelec(ii,is,js,rnew,snew)
implicit none
integer,intent(in) :: ii,is,js,snew
real(dp),intent(in) :: rnew(3)
call timer('DEFINE_CONFIG_ONEELEC',.true.)
if(buffer_move1_from(js)==is)then
call xyzzyaaaw1(is,js)
else
call clear_scratch_wfn(js)
endif
call xyzzyaaax1(js)
call define_config_oneelec_scratch(ii,is,js,rnew,snew)
call timer('DEFINE_CONFIG_ONEELEC',.false.)
end subroutine define_config_oneelec
subroutine define_config_twoelec(ii,jj,is,js,riinew,siinew,rjjnew,sjjn&
&ew)
implicit none
integer,intent(in) :: ii,jj,is,js,siinew,sjjnew
real(dp),intent(in) :: riinew(3),rjjnew(3)
call timer('DEFINE_CONFIG_TWOELEC',.true.)
if(buffer_move2_from(js)==is)then
call xyzzyaaaw1(is,js)
else
call clear_scratch_wfn(js)
endif
call xyzzyaaax1(js)
call define_config_twoelec_scratch(ii,jj,is,js,riinew,siinew,rjjnew,sj&
&jnew)
call timer('DEFINE_CONFIG_TWOELEC',.false.)
end subroutine define_config_twoelec
subroutine define_config_oneion(ion,is,js,rion1_new)
implicit none
integer,intent(in) :: ion,is,js
real(dp),intent(in) :: rion1_new(3)
call timer('DEFINE_CONFIG_ONEION',.true.)
if(buffer_move_ion_from(js)==is)then
call xyzzyaaaw1(is,js)
else
call clear_scratch_wfn(js)
if(allocated(erest_valid))then
if(erest_valid(is))then
ecomps_scr(:,js)=ecomps_scr(:,is)
erest_isinf_scr(js)=erest_isinf_scr(is)
erest_valid(js)=.true.
else
erest_valid(js)=.false.
endif
if(grid_angles_valid(is))then
grid_angles_scr(:,:,:,js)=grid_angles_scr(:,:,:,is)
grid_angles_valid(js)=.true.
else
grid_angles_valid(js)=.false.
endif
endif
endif
call xyzzyaaax1(js)
call define_config_oneion_scratch(ion,is,js,rion1_new)
call timer('DEFINE_CONFIG_ONEION',.false.)
end subroutine define_config_oneion
subroutine wfn_ratio(is,js,ilevel,ratio,relprob,truncprob,isnan,isinf,&
&prefetch_fd,prefetch_sd)
implicit none
integer,intent(in) :: is,js,ilevel
real(dp),intent(out),optional :: relprob,truncprob
complex(dp),intent(out),optional :: ratio
logical,intent(in),optional :: prefetch_fd,prefetch_sd
logical,intent(out),optional :: isnan,isinf
integer xyzzyaaaa9
real(dp) xyzzyaaab9,xyzzyaaac9
complex(dp) xyzzyaaad9,xyzzyaaae9
logical xyzzyaaaf9,xyzzyaaag9,xyzzyaaah9,xyzzyaaai9,xyzzyaaaj9,xyzzyaa&
&ak9
call timer('WFN_RATIO',.true.)
xyzzyaaaf9=.false.
xyzzyaaag9=.false.
xyzzyaaah9=.false.
if(present(prefetch_fd))xyzzyaaah9=prefetch_fd
xyzzyaaai9=.false.
if(present(prefetch_sd))xyzzyaaai9=prefetch_sd
xyzzyaaaj9=xyzzyaaae1/=0.and.(ilevel==0.or.xyzzyaaac1(ilevel)==xyzzyaa&
&aa1)
xyzzyaaak9=use_jastrow.and.(ilevel==0.or.xyzzyaaac1(ilevel)==xyzzyaaab&
&1)
xyzzyaaaa9=xyzzyaaad1(ilevel)
xyzzyaaae9=c_one
if(xyzzyaaaj9)then
select case(xyzzyaaae1)
case(1)
call wfn_ratio_slater(is,js,xyzzyaaaa9,xyzzyaaad9,xyzzyaaah9,xyzzyaaai&
&9,xyzzyaaaf9,xyzzyaaag9)
case(2)
call wfn_ratio_exmol(is,js,xyzzyaaaa9,xyzzyaaad9,xyzzyaaah9,xyzzyaaai9&
&)
case(3)
call wfn_ratio_pfaff(is,js,xyzzyaaaa9,xyzzyaaad9,xyzzyaaah9,xyzzyaaai9&
&)
case(4)
call wfn_ratio_geminal(is,js,xyzzyaaaa9,xyzzyaaad9,xyzzyaaah9,xyzzyaaa&
&i9,xyzzyaaaf9,xyzzyaaag9)
case(6)
call wfn_ratio_mahan(is,js,xyzzyaaaa9,xyzzyaaad9,xyzzyaaah9,xyzzyaaai9&
&,xyzzyaaaf9,xyzzyaaag9)
end select
xyzzyaaae9=xyzzyaaae9*xyzzyaaad9
endif
if(xyzzyaaak9)then
call wfn_ratio_jastrow(is,js,xyzzyaaaa9,xyzzyaaac9,xyzzyaaah9,xyzzyaaa&
&i9)
xyzzyaaae9=xyzzyaaae9*xyzzyaaac9
endif
if(present(ratio))ratio=xyzzyaaae9
xyzzyaaab9=dble(xyzzyaaae9)**2+aimag(xyzzyaaae9)**2
if(present(relprob))relprob=xyzzyaaab9
if(present(truncprob))truncprob=min(xyzzyaaab9,1.d0)
if(present(isnan))isnan=xyzzyaaaf9
if(present(isinf))isinf=xyzzyaaag9
call timer('WFN_RATIO',.false.)
end subroutine wfn_ratio
subroutine accept_move(is,js)
implicit none
integer,intent(in) :: is,js
call timer('ACCEPT_MOVE',.true.)
select case(xyzzyaaae1)
case(1)
call accept_move_slater(is,js)
case(2)
call accept_move_exmol(is,js)
case(3)
call accept_move_pfaff(is,js)
case(4)
call accept_move_geminal(is,js)
case(6)
call accept_move_mahan(is,js)
end select
if(use_jastrow)call accept_move_jastrow(is,js)
if(xyzzyaaaf1)call accept_move_gbasis(is,js)
if(allocated(etot_valid))then
etot_scr(is)=etot_scr(js)
keimag_scr(is)=keimag_scr(js)
etot_isnan_scr(is)=etot_isnan_scr(js)
etot_isinf_scr(is)=etot_isinf_scr(js)
ecomps_scr(:,is)=ecomps_scr(:,js)
etot_valid(is)=etot_valid(js)
ke_isnan_scr(is)=ke_isnan_scr(js)
ke_isinf_scr(is)=ke_isinf_scr(js)
ke_valid(is)=ke_valid(js)
nl_isnan_scr(is)=nl_isnan_scr(js)
nl_isinf_scr(is)=nl_isinf_scr(js)
nl_valid(is)=nl_valid(js)
erest_isinf_scr(is)=erest_isinf_scr(js)
erest_valid(is)=erest_valid(js)
if(grid_angles_valid(js))grid_angles_scr(:,:,:,is)=grid_angles_scr(:,:&
&,:,js)
grid_angles_valid(is)=grid_angles_valid(js)
endif
call accept_move_scratch(is,js)
call xyzzyaaax1(is)
call timer('ACCEPT_MOVE',.false.)
end subroutine accept_move
subroutine xyzzyaaaw1(is,js)
implicit none
integer,intent(in) :: is,js
select case(xyzzyaaae1)
case(1)
call reset_config_slater(is,js)
case(2)
call reset_config_exmol(is,js)
case(3)
call reset_config_pfaff(is,js)
case(4)
call reset_config_geminal(is,js)
case(6)
call reset_config_mahan(is,js)
end select
if(use_jastrow)call reset_config_jastrow(is,js)
if(xyzzyaaaf1)call reset_config_gbasis(is,js)
if(allocated(etot_valid))then
etot_valid(js)=.false.
ke_valid(js)=.false.
nl_valid(js)=.false.
if(buffer_move_ion_from(js)==is)then
if(.not.erest_valid(js))then
if(erest_valid(is))then
ecomps_scr(:,js)=ecomps_scr(:,is)
erest_isinf_scr(js)=erest_isinf_scr(is)
erest_valid(js)=.true.
else
erest_valid(js)=.false.
endif
else
erest_valid(js)=.true.
endif
if(.not.grid_angles_valid(is))then
if(grid_angles_valid(is))then
grid_angles_scr(:,:,:,js)=grid_angles_scr(:,:,:,is)
grid_angles_valid(js)=.true.
else
grid_angles_valid(js)=.false.
endif
else
grid_angles_valid(js)=.true.
endif
else
erest_valid(js)=.false.
grid_angles_valid(js)=.false.
endif
endif
call reset_config_scratch(is,js)
end subroutine xyzzyaaaw1
subroutine clear_scratch_wfn(is,wfn_only)
implicit none
integer,intent(in) :: is
logical,intent(in),optional :: wfn_only
select case(xyzzyaaae1)
case(1)
call clear_scratch_slater(is)
case(2)
call clear_scratch_exmol(is)
case(3)
call clear_scratch_pfaff(is)
case(4)
call clear_scratch_geminal(is)
case(6)
call clear_scratch_mahan(is)
end select
if(use_jastrow)call clear_scratch_jastrow(is)
if(xyzzyaaaf1)call clear_scratch_gbasis(is)
if(allocated(etot_valid))then
etot_valid(is)=.false.
ke_valid(is)=.false.
nl_valid(is)=.false.
erest_valid(is)=.false.
grid_angles_valid(is)=.false.
endif
if(present(wfn_only))then
if(wfn_only)return
endif
call clear_scratch(is)
end subroutine clear_scratch_wfn
subroutine xyzzyaaax1(is)
implicit none
integer,intent(in) :: is
integer xyzzyaaaa13
do xyzzyaaaa13=1,nscratch
if(buffer_move1_from(xyzzyaaaa13)==is.or.buffer_move2_from(xyzzyaaaa13&
&)==is.or.buffer_move_ion_from(xyzzyaaaa13)==is)call clear_scratch_wfn&
&(xyzzyaaaa13)
enddo
end subroutine xyzzyaaax1
subroutine empty_scratch_wfn(is)
implicit none
integer,intent(in) :: is
select case(xyzzyaaae1)
case(1)
call clear_scratch_slater(is)
case(2)
call clear_scratch_exmol(is)
case(3)
call clear_scratch_pfaff(is)
case(4)
call clear_scratch_geminal(is)
case(6)
call clear_scratch_mahan(is)
end select
if(use_jastrow)call clear_scratch_jastrow(is)
if(xyzzyaaaf1)call clear_scratch_gbasis(is)
if(allocated(etot_valid))then
etot_valid(is)=.false.
ke_valid(is)=.false.
nl_valid(is)=.false.
erest_valid(is)=.false.
endif
end subroutine empty_scratch_wfn
subroutine wfn_logval(is,logwfn,iszero)
implicit none
integer,intent(in) :: is
complex(dp),intent(out) :: logwfn
logical,intent(out),optional :: iszero
real(dp) xyzzyaaaa15
complex(dp) xyzzyaaab15
logical xyzzyaaac15,xyzzyaaad15,xyzzyaaae15
call timer('WFN_LOGVAL',.true.)
if(use_altsamp.and.simplepdf==1.and.xyzzyaaae1/=1)call errstop('WFN_LO&
&GVAL','Simplified sampling not available for this wave function type.&
&')
xyzzyaaac15=(xyzzyaaae1/=0)
xyzzyaaad15=use_jastrow
logwfn=czero
xyzzyaaae15=.false.
if(xyzzyaaac15)then
select case(xyzzyaaae1)
case(1)
call wfn_logval_slater(is,xyzzyaaab15,xyzzyaaae15)
case(2)
call wfn_logval_exmol(is,xyzzyaaab15)
case(3)
call wfn_logval_pfaff(is,xyzzyaaab15)
case(4)
call wfn_logval_geminal(is,xyzzyaaab15,xyzzyaaae15)
case(6)
call wfn_logval_mahan(is,xyzzyaaab15,xyzzyaaae15)
end select
if(.not.xyzzyaaae15)logwfn=logwfn+xyzzyaaab15
endif
if(.not.xyzzyaaae15.and.xyzzyaaad15)then
call wfn_logval_jastrow(is,xyzzyaaaa15)
logwfn=logwfn+cmplx(xyzzyaaaa15,0.d0,dp)
endif
if(present(iszero))iszero=xyzzyaaae15
call timer('WFN_LOGVAL',.false.)
end subroutine wfn_logval
subroutine wfn_loggrad(ii,is,ilevel,loggrad_psi,isnan,isinf,prefetch_v&
&al,prefetch_sd,prefetch_aderiv)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: loggrad_psi(3)
logical,intent(in),optional :: prefetch_val,prefetch_sd,prefetch_aderi&
&v
logical,intent(out),optional :: isnan,isinf
integer xyzzyaaaa16
real(dp) xyzzyaaab16(3)
complex(dp) xyzzyaaac16(3)
logical xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,xyzzyaaag16,xyzzyaaah16,xy&
&zzyaaai16,xyzzyaaaj16
call timer('WFN_LOGGRAD',.true.)
xyzzyaaad16=.false.
if(present(prefetch_val))xyzzyaaad16=prefetch_val
xyzzyaaae16=.false.
if(present(prefetch_sd))xyzzyaaae16=prefetch_sd
xyzzyaaaf16=.false.
if(present(prefetch_aderiv))xyzzyaaaf16=prefetch_aderiv
xyzzyaaai16=.false.
xyzzyaaaj16=.false.
xyzzyaaag16=xyzzyaaae1/=0.and.(ilevel==0.or.xyzzyaaac1(ilevel)==xyzzya&
&aaa1)
xyzzyaaah16=use_jastrow.and.(ilevel==0.or.xyzzyaaac1(ilevel)==xyzzyaaa&
&b1)
xyzzyaaaa16=xyzzyaaad1(ilevel)
loggrad_psi=czero
if(xyzzyaaag16)then
select case(xyzzyaaae1)
case(1)
call wfn_loggrad_slater(ii,is,xyzzyaaaa16,xyzzyaaad16,xyzzyaaae16,xyzz&
&yaaac16,xyzzyaaai16,xyzzyaaaj16)
case(2)
call wfn_loggrad_exmol(ii,is,xyzzyaaaa16,xyzzyaaad16,xyzzyaaae16,xyzzy&
&aaac16)
case(3)
call wfn_loggrad_pfaff(ii,is,xyzzyaaaa16,xyzzyaaad16,xyzzyaaae16,xyzzy&
&aaac16)
case(4)
call wfn_loggrad_geminal(ii,is,xyzzyaaaa16,xyzzyaaad16,xyzzyaaae16,xyz&
&zyaaac16,xyzzyaaai16,xyzzyaaaj16)
case(6)
call wfn_loggrad_mahan(ii,is,xyzzyaaaa16,xyzzyaaad16,xyzzyaaae16,xyzzy&
&aaac16,xyzzyaaai16,xyzzyaaaj16)
end select
loggrad_psi=loggrad_psi+xyzzyaaac16
endif
if(xyzzyaaah16)then
call wfn_loggrad_jastrow(ii,is,xyzzyaaaa16,xyzzyaaad16,xyzzyaaae16,xyz&
&zyaaaf16,xyzzyaaab16)
loggrad_psi=loggrad_psi+cmplx(xyzzyaaab16(1:3),0.d0,dp)
endif
if(present(isnan))isnan=xyzzyaaai16
if(present(isinf))isinf=xyzzyaaaj16
call timer('WFN_LOGGRAD',.false.)
end subroutine wfn_loggrad
subroutine wfn_loglap(ii,is,ilevel,loglap_psi,isnan,isinf,prefetch_val&
&,prefetch_fd,prefetch_aderiv)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: loglap_psi
logical,intent(in),optional :: prefetch_val,prefetch_fd,prefetch_aderi&
&v
logical,intent(out),optional :: isnan,isinf
integer xyzzyaaaa17
real(dp) xyzzyaaab17
complex(dp) xyzzyaaac17
logical xyzzyaaad17,xyzzyaaae17,xyzzyaaaf17,xyzzyaaag17,xyzzyaaah17,xy&
&zzyaaai17,xyzzyaaaj17
call timer('WFN_LOGLAP',.true.)
xyzzyaaad17=.false.
if(present(prefetch_val))xyzzyaaad17=prefetch_val
xyzzyaaae17=.false.
if(present(prefetch_fd))xyzzyaaae17=prefetch_fd
xyzzyaaaf17=.false.
if(present(prefetch_aderiv))xyzzyaaaf17=prefetch_aderiv
xyzzyaaai17=.false.
xyzzyaaaj17=.false.
xyzzyaaag17=xyzzyaaae1/=0.and.(ilevel==0.or.xyzzyaaac1(ilevel)==xyzzya&
&aaa1)
xyzzyaaah17=use_jastrow.and.(ilevel==0.or.xyzzyaaac1(ilevel)==xyzzyaaa&
&b1)
xyzzyaaaa17=xyzzyaaad1(ilevel)
loglap_psi=czero
if(xyzzyaaag17)then
select case(xyzzyaaae1)
case(1)
call wfn_loglap_slater(ii,is,xyzzyaaaa17,xyzzyaaad17,xyzzyaaac17,xyzzy&
&aaai17,xyzzyaaaj17)
case(2)
call wfn_loglap_exmol(ii,is,xyzzyaaaa17,xyzzyaaad17,xyzzyaaae17,xyzzya&
&aac17)
case(3)
call wfn_loglap_pfaff(ii,is,xyzzyaaaa17,xyzzyaaad17,xyzzyaaae17,xyzzya&
&aac17)
case(4)
call wfn_loglap_geminal(ii,is,xyzzyaaaa17,xyzzyaaad17,xyzzyaaae17,xyzz&
&yaaac17,xyzzyaaai17,xyzzyaaaj17)
case(6)
call wfn_loglap_mahan(ii,is,xyzzyaaaa17,xyzzyaaad17,xyzzyaaae17,xyzzya&
&aac17,xyzzyaaai17,xyzzyaaaj17)
end select
loglap_psi=loglap_psi+xyzzyaaac17
endif
if(xyzzyaaah17)then
call wfn_loglap_jastrow(ii,is,xyzzyaaaa17,xyzzyaaad17,xyzzyaaae17,xyzz&
&yaaaf17,xyzzyaaab17)
loglap_psi=loglap_psi+cmplx(xyzzyaaab17,0.d0,dp)
endif
if(present(isnan))isnan=xyzzyaaai17
if(present(isinf))isinf=xyzzyaaaj17
call timer('WFN_LOGLAP',.false.)
end subroutine wfn_loglap
subroutine prefetch_wfn(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
call timer('PREFETCH_WFN',.true.)
select case(xyzzyaaae1)
case(1)
call prefetch_wfn_slater(is,fd,sd)
case(2)
call prefetch_wfn_exmol(is,fd,sd)
case(3)
call prefetch_wfn_pfaff(is,fd,sd)
case(4)
call prefetch_wfn_geminal(is,fd,sd)
case(6)
call prefetch_wfn_mahan(is,fd,sd)
end select
if(use_jastrow)call prefetch_wfn_jastrow(is,fd,sd)
call timer('PREFETCH_WFN',.false.)
end subroutine prefetch_wfn
subroutine add_config_wfn_items(is)
implicit none
integer,intent(in) :: is
select case(xyzzyaaae1)
case(1)
call add_config_slater_items(is)
case(2)
continue
case(3)
continue
case(4)
continue
case(6)
continue
end select
if(use_jastrow)call add_config_jastrow_items(is)
end subroutine add_config_wfn_items
subroutine clone_scratch(is,js,reset)
implicit none
integer,intent(in) :: is,js
logical,intent(in),optional :: reset
logical xyzzyaaaa20
call timer('CLONE_SCRATCH',.true.)
xyzzyaaaa20=.true.
if(present(reset))xyzzyaaaa20=reset
if(xyzzyaaaa20)then
call clear_scratch_wfn(js)
call xyzzyaaax1(js)
endif
call clone_scratch_geom(is,js)
select case(xyzzyaaae1)
case(1)
call clone_scratch_slater(is,js)
case(2)
call clone_scratch_exmol(is,js)
case(3)
call clone_scratch_pfaff(is,js)
case(4)
call clone_scratch_geminal(is,js)
case(6)
call clone_scratch_mahan(is,js)
end select
if(use_jastrow)call clone_scratch_jastrow(is,js)
call timer('CLONE_SCRATCH',.false.)
end subroutine clone_scratch
subroutine gen_config(pt_config_geom,pt_config_wfn)
implicit none
type(config_geom),pointer :: pt_config_geom
type(config_wfn),pointer :: pt_config_wfn
integer xyzzyaaaa21
call gen_config_geom(pt_config_geom)
allocate(pt_config_wfn,stat=xyzzyaaaa21)
call check_alloc(xyzzyaaaa21,'GEN_CONFIG_WFN','container')
nullify(pt_config_wfn%pt_slater,pt_config_wfn%pt_jastrow)
select case(xyzzyaaae1)
case(1)
call gen_config_slater(pt_config_wfn%pt_slater)
case(2)
call gen_config_exmol(pt_config_wfn%pt_exmol)
case(3)
call gen_config_pfaff(pt_config_wfn%pt_pfaff)
case(4)
call gen_config_geminal(pt_config_wfn%pt_geminal)
case(6)
call gen_config_mahan(pt_config_wfn%pt_mahan)
end select
if(use_jastrow)call gen_config_jastrow(pt_config_wfn%pt_jastrow)
end subroutine gen_config
subroutine load_from_pt(is,pt_config_geom,pt_config_wfn,rele,sele)
implicit none
integer,intent(in) :: is
integer,intent(out) :: sele(netot)
real(dp),intent(out) :: rele(3,netot)
type(config_geom),pointer :: pt_config_geom
type(config_wfn),pointer :: pt_config_wfn
call clear_scratch_wfn(is)
call xyzzyaaax1(is)
call load_from_pt_geom(is,pt_config_geom,rele,sele)
select case(xyzzyaaae1)
case(1)
call load_from_pt_slater(is,pt_config_wfn%pt_slater)
case(2)
call load_from_pt_exmol(is,pt_config_wfn%pt_exmol)
case(3)
call load_from_pt_pfaff(is,pt_config_wfn%pt_pfaff)
case(4)
call load_from_pt_geminal(is,pt_config_wfn%pt_geminal)
case(6)
call load_from_pt_mahan(is,pt_config_wfn%pt_mahan)
end select
if(use_jastrow)call load_from_pt_jastrow(is,pt_config_wfn%pt_jastrow)
end subroutine load_from_pt
subroutine save_to_pt(is,pt_config_geom,pt_config_wfn)
implicit none
integer,intent(in) :: is
type(config_geom),pointer :: pt_config_geom
type(config_wfn),pointer :: pt_config_wfn
call save_to_pt_geom(is,pt_config_geom)
select case(xyzzyaaae1)
case(1)
call save_to_pt_slater(is,pt_config_wfn%pt_slater)
case(2)
call save_to_pt_exmol(is,pt_config_wfn%pt_exmol)
case(3)
call save_to_pt_pfaff(is,pt_config_wfn%pt_pfaff)
case(4)
call save_to_pt_geminal(is,pt_config_wfn%pt_geminal)
case(6)
call save_to_pt_mahan(is,pt_config_wfn%pt_mahan)
end select
if(use_jastrow)call save_to_pt_jastrow(is,pt_config_wfn%pt_jastrow)
end subroutine save_to_pt
subroutine delete_config(pt_config_geom,pt_config_wfn)
implicit none
type(config_geom),pointer :: pt_config_geom
type(config_wfn),pointer :: pt_config_wfn
select case(xyzzyaaae1)
case(1)
call delete_config_slater(pt_config_wfn%pt_slater)
case(2)
call delete_config_exmol(pt_config_wfn%pt_exmol)
case(3)
call delete_config_pfaff(pt_config_wfn%pt_pfaff)
case(4)
call delete_config_geminal(pt_config_wfn%pt_geminal)
case(6)
call delete_config_mahan(pt_config_wfn%pt_mahan)
end select
if(use_jastrow)call delete_config_jastrow(pt_config_wfn%pt_jastrow)
deallocate(pt_config_wfn)
call delete_config_geom(pt_config_geom)
end subroutine delete_config
subroutine copy_config(pt_from_geom,pt_to_geom,pt_from_wfn,pt_to_wfn)
implicit none
type(config_geom),pointer :: pt_from_geom,pt_to_geom
type(config_wfn),pointer :: pt_from_wfn,pt_to_wfn
call copy_config_geom(pt_from_geom,pt_to_geom)
select case(xyzzyaaae1)
case(1)
call copy_config_slater(pt_from_wfn%pt_slater,pt_to_wfn%pt_slater)
case(2)
call copy_config_exmol(pt_from_wfn%pt_exmol,pt_to_wfn%pt_exmol)
case(3)
call copy_config_pfaff(pt_from_wfn%pt_pfaff,pt_to_wfn%pt_pfaff)
case(4)
call copy_config_geminal(pt_from_wfn%pt_geminal,pt_to_wfn%pt_geminal)
case(6)
call copy_config_mahan(pt_from_wfn%pt_mahan,pt_to_wfn%pt_mahan)
end select
if(use_jastrow)call copy_config_jastrow(pt_from_wfn%pt_jastrow,pt_to_w&
&fn%pt_jastrow)
end subroutine copy_config
subroutine config_to_pt(pt_config_geom,pt_config_wfn,k)
implicit none
integer,intent(in) :: k
type(config_geom),pointer :: pt_config_geom
type(config_wfn),pointer :: pt_config_wfn
call config_to_pt_geom(pt_config_geom,k)
select case(xyzzyaaae1)
case(1)
call config_to_pt_slater(pt_config_wfn%pt_slater,k)
case(2)
call config_to_pt_exmol(pt_config_wfn%pt_exmol,k)
case(3)
call config_to_pt_pfaff(pt_config_wfn%pt_pfaff,k)
case(4)
call config_to_pt_geminal(pt_config_wfn%pt_geminal,k)
case(6)
call config_to_pt_mahan(pt_config_wfn%pt_mahan,k)
end select
if(use_jastrow)call config_to_pt_jastrow(pt_config_wfn%pt_jastrow,k)
end subroutine config_to_pt
subroutine pt_to_config(pt_config_geom,pt_config_wfn)
implicit none
type(config_geom),pointer :: pt_config_geom
type(config_wfn),pointer :: pt_config_wfn
call pt_to_config_geom(pt_config_geom)
select case(xyzzyaaae1)
case(1)
call pt_to_config_slater(pt_config_wfn%pt_slater)
case(2)
call pt_to_config_exmol(pt_config_wfn%pt_exmol)
case(3)
call pt_to_config_pfaff(pt_config_wfn%pt_pfaff)
case(4)
call pt_to_config_geminal(pt_config_wfn%pt_geminal)
case(6)
call pt_to_config_mahan(pt_config_wfn%pt_mahan)
end select
if(use_jastrow)call pt_to_config_jastrow(pt_config_wfn%pt_jastrow)
end subroutine pt_to_config
subroutine redist_allocations(kmax)
implicit none
integer,intent(in) :: kmax
call redist_allocations_geom(kmax)
select case(xyzzyaaae1)
case(1)
call redist_allocations_slater(kmax)
case(2)
call redist_allocations_exmol(kmax)
case(3)
call redist_allocations_pfaff(kmax)
case(4)
call redist_allocations_geminal(kmax)
case(6)
call redist_allocations_mahan(kmax)
end select
if(use_jastrow)call redist_allocations_jastrow(kmax)
end subroutine redist_allocations
subroutine redist_load(pt_config_geom,pt_config_wfn,k,blocking)
implicit none
integer,intent(in) :: k
logical,intent(in) :: blocking
type(config_geom),pointer :: pt_config_geom
type(config_wfn),pointer :: pt_config_wfn
call redist_load_geom(pt_config_geom,k)
select case(xyzzyaaae1)
case(1)
call redist_load_slater(pt_config_wfn%pt_slater,k)
case(2)
call redist_load_exmol(pt_config_wfn%pt_exmol,k)
case(3)
call redist_load_pfaff(pt_config_wfn%pt_pfaff,k)
case(4)
call redist_load_geminal(pt_config_wfn%pt_geminal,k)
case(6)
call redist_load_mahan(pt_config_wfn%pt_mahan,k)
end select
if(use_jastrow)call redist_load_jastrow(pt_config_wfn%pt_jastrow,k,blo&
&cking)
end subroutine redist_load
subroutine redist_send(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
call redist_send_geom(jnode,kbase,k,nbt,reqbase,blocking)
select case(xyzzyaaae1)
case(1)
call redist_send_slater(jnode,kbase,k,nbt,reqbase,blocking)
case(2)
call redist_send_exmol(jnode,kbase,k,nbt,reqbase,blocking)
case(3)
call redist_send_pfaff(jnode,kbase,k,nbt,reqbase,blocking)
case(4)
call redist_send_geminal(jnode,kbase,k,nbt,reqbase,blocking)
case(6)
call redist_send_mahan(jnode,kbase,k,nbt,reqbase,blocking)
end select
if(use_jastrow)call redist_send_jastrow(jnode,kbase,k,nbt,reqbase,bloc&
&king)
end subroutine redist_send
subroutine redist_recv(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
call redist_recv_geom(jnode,kbase,k,nbt,reqbase,blocking)
select case(xyzzyaaae1)
case(1)
call redist_recv_slater(jnode,kbase,k,nbt,reqbase,blocking)
case(2)
call redist_recv_exmol(jnode,kbase,k,nbt,reqbase,blocking)
case(3)
call redist_recv_pfaff(jnode,kbase,k,nbt,reqbase,blocking)
case(4)
call redist_recv_geminal(jnode,kbase,k,nbt,reqbase,blocking)
case(6)
call redist_recv_mahan(jnode,kbase,k,nbt,reqbase,blocking)
end select
if(use_jastrow)call redist_recv_jastrow(jnode,kbase,k,nbt,reqbase,bloc&
&king)
end subroutine redist_recv
subroutine redist_save(pt_config_geom,pt_config_wfn,k)
implicit none
integer,intent(in) :: k
type(config_geom),pointer :: pt_config_geom
type(config_wfn),pointer :: pt_config_wfn
call redist_save_geom(pt_config_geom,k)
select case(xyzzyaaae1)
case(1)
call redist_save_slater(pt_config_wfn%pt_slater,k)
case(2)
call redist_save_exmol(pt_config_wfn%pt_exmol,k)
case(3)
call redist_save_pfaff(pt_config_wfn%pt_pfaff,k)
case(4)
call redist_save_geminal(pt_config_wfn%pt_geminal,k)
case(6)
call redist_save_mahan(pt_config_wfn%pt_mahan,k)
end select
if(use_jastrow)call redist_save_jastrow(pt_config_wfn%pt_jastrow,k)
end subroutine redist_save
subroutine redist_deallocations
implicit none
call redist_deallocations_geom
select case(xyzzyaaae1)
case(1)
call redist_deallocations_slater
case(2)
call redist_deallocations_exmol
case(3)
call redist_deallocations_pfaff
case(4)
call redist_deallocations_geminal
case(6)
call redist_deallocations_mahan
end select
if(use_jastrow)call redist_deallocations_jastrow
end subroutine redist_deallocations
subroutine setup_wfn_params(nparam)
implicit none
integer,intent(out) :: nparam
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34,xyzzyaaad34
select case(xyzzyaaae1)
case(1)
call setup_slater_params(xyzzyaaai1(1))
case(2)
call setup_exmol_params(xyzzyaaai1(1))
case(3)
call setup_pfaff_params(xyzzyaaai1(1))
case(4)
call setup_geminal_params(xyzzyaaai1(1))
case(6)
call setup_mahan_params(xyzzyaaai1(1))
end select
if(use_jastrow)call setup_jastrow_params(xyzzyaaai1(2))
if(xyzzyaaaf1)call setup_gbasis_params
xyzzyaaah1=sum(xyzzyaaai1,xyzzyaaai1>0)
xyzzyaaan1=.true.
allocate(xyzzyaaaj1(xyzzyaaah1),xyzzyaaao1(xyzzyaaah1),xyzzyaaak1(xyzz&
&yaaah1),xyzzyaaap1(xyzzyaaah1),xyzzyaaal1(xyzzyaaah1),xyzzyaaaq1(xyzz&
&yaaah1),xyzzyaaar1(xyzzyaaah1,xyzzyaaah1),xyzzyaaam1(xyzzyaaah1),stat&
&=xyzzyaaad34)
call check_alloc(xyzzyaaad34,'SETUP_WFN_PARAMS','limits')
xyzzyaaac34=0
do xyzzyaaaa34=1,xyzzyaaag1
xyzzyaaab34=xyzzyaaac34+1
xyzzyaaac34=xyzzyaaac34+xyzzyaaai1(xyzzyaaaa34)
xyzzyaaaj1(xyzzyaaab34:xyzzyaaac34)=xyzzyaaaa34
enddo
nparam=xyzzyaaah1
end subroutine setup_wfn_params
subroutine finish_wfn_params
implicit none
select case(xyzzyaaae1)
case(1)
call finish_slater_params
case(2)
call finish_exmol_params
case(3)
call finish_pfaff_params
case(4)
call finish_geminal_params
case(6)
call finish_mahan_params
end select
if(use_jastrow.and.opt_jastrow)call finish_jastrow_params
if(xyzzyaaaf1)call finish_gbasis_params
deallocate(xyzzyaaaj1,xyzzyaaao1,xyzzyaaak1,xyzzyaaap1,xyzzyaaal1,xyzz&
&yaaaq1,xyzzyaaar1,xyzzyaaam1)
end subroutine finish_wfn_params
subroutine get_params(params,is_shallow,is_redundant,is_linear,is_logl&
&inear,has_aderiv_out,label,rescale)
use slaarnabt, only : param_unlimiting
implicit none
real(dp),intent(inout) :: params(xyzzyaaah1)
logical,intent(in),optional :: rescale
logical,intent(inout) :: is_shallow(xyzzyaaah1),is_redundant(xyzzyaaah&
&1),is_linear(xyzzyaaah1),is_loglinear(xyzzyaaah1),has_aderiv_out(xyzz&
&yaaah1)
character(2),intent(inout) :: label(xyzzyaaah1)
integer xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36
real(dp) xyzzyaaad36,xyzzyaaae36,xyzzyaaaf36,xyzzyaaag36,xyzzyaaah36
if(present(rescale))xyzzyaaan1=rescale
xyzzyaaao1=.false.
xyzzyaaak1=0.d0
xyzzyaaap1=.false.
xyzzyaaal1=0.d0
is_shallow=.false.
is_redundant=.false.
is_linear=.false.
is_loglinear=.false.
xyzzyaaaq1=.false.
xyzzyaaar1=.false.
do xyzzyaaac36=1,xyzzyaaah1
xyzzyaaar1(xyzzyaaac36,xyzzyaaac36)=.true.
enddo
label='??'
xyzzyaaab36=0
if(xyzzyaaai1(1)>0)then
xyzzyaaaa36=xyzzyaaab36+1
xyzzyaaab36=xyzzyaaab36+xyzzyaaai1(1)
select case(xyzzyaaae1)
case(1)
call get_slater_params(params(xyzzyaaaa36:xyzzyaaab36),xyzzyaaao1(xyzz&
&yaaaa36:xyzzyaaab36),xyzzyaaak1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaap1(x&
&yzzyaaaa36:xyzzyaaab36),xyzzyaaal1(xyzzyaaaa36:xyzzyaaab36),is_shallo&
&w(xyzzyaaaa36:xyzzyaaab36),is_redundant(xyzzyaaaa36:xyzzyaaab36),is_l&
&inear(xyzzyaaaa36:xyzzyaaab36),is_loglinear(xyzzyaaaa36:xyzzyaaab36),&
&xyzzyaaaq1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaar1(xyzzyaaaa36:xyzzyaaab3&
&6,xyzzyaaaa36:xyzzyaaab36),label(xyzzyaaaa36:xyzzyaaab36))
case(2)
call get_exmol_params(params(xyzzyaaaa36:xyzzyaaab36),xyzzyaaao1(xyzzy&
&aaaa36:xyzzyaaab36),xyzzyaaak1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaap1(xy&
&zzyaaaa36:xyzzyaaab36),xyzzyaaal1(xyzzyaaaa36:xyzzyaaab36),is_shallow&
&(xyzzyaaaa36:xyzzyaaab36),is_redundant(xyzzyaaaa36:xyzzyaaab36),is_li&
&near(xyzzyaaaa36:xyzzyaaab36),is_loglinear(xyzzyaaaa36:xyzzyaaab36),x&
&yzzyaaaq1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaar1(xyzzyaaaa36:xyzzyaaab36&
&,xyzzyaaaa36:xyzzyaaab36),label(xyzzyaaaa36:xyzzyaaab36))
case(3)
call get_pfaff_params(params(xyzzyaaaa36:xyzzyaaab36),xyzzyaaao1(xyzzy&
&aaaa36:xyzzyaaab36),xyzzyaaak1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaap1(xy&
&zzyaaaa36:xyzzyaaab36),xyzzyaaal1(xyzzyaaaa36:xyzzyaaab36),is_shallow&
&(xyzzyaaaa36:xyzzyaaab36),is_redundant(xyzzyaaaa36:xyzzyaaab36),is_li&
&near(xyzzyaaaa36:xyzzyaaab36),is_loglinear(xyzzyaaaa36:xyzzyaaab36),x&
&yzzyaaaq1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaar1(xyzzyaaaa36:xyzzyaaab36&
&,xyzzyaaaa36:xyzzyaaab36),label(xyzzyaaaa36:xyzzyaaab36))
case(4)
call get_geminal_params(params(xyzzyaaaa36:xyzzyaaab36),xyzzyaaao1(xyz&
&zyaaaa36:xyzzyaaab36),xyzzyaaak1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaap1(&
&xyzzyaaaa36:xyzzyaaab36),xyzzyaaal1(xyzzyaaaa36:xyzzyaaab36),is_shall&
&ow(xyzzyaaaa36:xyzzyaaab36),is_redundant(xyzzyaaaa36:xyzzyaaab36),is_&
&linear(xyzzyaaaa36:xyzzyaaab36),is_loglinear(xyzzyaaaa36:xyzzyaaab36)&
&,xyzzyaaaq1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaar1(xyzzyaaaa36:xyzzyaaab&
&36,xyzzyaaaa36:xyzzyaaab36),label(xyzzyaaaa36:xyzzyaaab36))
case(6)
call get_mahan_params(params(xyzzyaaaa36:xyzzyaaab36),xyzzyaaao1(xyzzy&
&aaaa36:xyzzyaaab36),xyzzyaaak1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaap1(xy&
&zzyaaaa36:xyzzyaaab36),xyzzyaaal1(xyzzyaaaa36:xyzzyaaab36),is_shallow&
&(xyzzyaaaa36:xyzzyaaab36),is_redundant(xyzzyaaaa36:xyzzyaaab36),is_li&
&near(xyzzyaaaa36:xyzzyaaab36),is_loglinear(xyzzyaaaa36:xyzzyaaab36),x&
&yzzyaaaq1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaar1(xyzzyaaaa36:xyzzyaaab36&
&,xyzzyaaaa36:xyzzyaaab36),label(xyzzyaaaa36:xyzzyaaab36))
end select
endif
if(xyzzyaaai1(2)>0)then
xyzzyaaaa36=xyzzyaaab36+1
xyzzyaaab36=xyzzyaaab36+xyzzyaaai1(2)
call get_jastrow_params(params(xyzzyaaaa36:xyzzyaaab36),xyzzyaaao1(xyz&
&zyaaaa36:xyzzyaaab36),xyzzyaaak1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaap1(&
&xyzzyaaaa36:xyzzyaaab36),xyzzyaaal1(xyzzyaaaa36:xyzzyaaab36),is_shall&
&ow(xyzzyaaaa36:xyzzyaaab36),is_redundant(xyzzyaaaa36:xyzzyaaab36),is_&
&linear(xyzzyaaaa36:xyzzyaaab36),is_loglinear(xyzzyaaaa36:xyzzyaaab36)&
&,xyzzyaaaq1(xyzzyaaaa36:xyzzyaaab36),xyzzyaaar1(xyzzyaaaa36:xyzzyaaab&
&36,xyzzyaaaa36:xyzzyaaab36),label(xyzzyaaaa36:xyzzyaaab36))
endif
has_aderiv_out=xyzzyaaaq1
do xyzzyaaac36=1,xyzzyaaah1
xyzzyaaad36=params(xyzzyaaac36)
xyzzyaaaf36=xyzzyaaal1(xyzzyaaac36)
xyzzyaaag36=xyzzyaaak1(xyzzyaaac36)
xyzzyaaae36=1.d0
if(xyzzyaaan1.and.abs(xyzzyaaad36)>sqrt(tiny(1.d0)))xyzzyaaae36=abs(xy&
&zzyaaad36)
xyzzyaaah36=1.d0/xyzzyaaae36
xyzzyaaad36=xyzzyaaad36*xyzzyaaah36
xyzzyaaaf36=xyzzyaaaf36*xyzzyaaah36
xyzzyaaag36=xyzzyaaag36*xyzzyaaah36
xyzzyaaad36=param_unlimiting(xyzzyaaan1.and.vm_smooth_limits,xyzzyaaad&
&36,xyzzyaaap1(xyzzyaaac36),xyzzyaaaf36,xyzzyaaao1(xyzzyaaac36),xyzzya&
&aag36)
params(xyzzyaaac36)=xyzzyaaad36
xyzzyaaam1(xyzzyaaac36)=xyzzyaaae36
xyzzyaaal1(xyzzyaaac36)=xyzzyaaaf36
xyzzyaaak1(xyzzyaaac36)=xyzzyaaag36
enddo
end subroutine get_params
subroutine put_params(params,ignore,prestore,params_print,bad_params,r&
&estore_hint,quirky)
use slaarnabt, only : param_limiting
implicit none
integer,intent(in),optional :: restore_hint
real(dp),intent(in) :: params(:)
real(dp),intent(out),optional :: params_print(:)
logical,intent(in) :: ignore(xyzzyaaah1),prestore
logical,intent(in),optional :: quirky
logical,intent(inout),optional :: bad_params
integer xyzzyaaaa37,xyzzyaaab37,xyzzyaaac37,xyzzyaaad37,xyzzyaaae37
real(dp) xyzzyaaaf37(xyzzyaaah1),xyzzyaaag37,xyzzyaaah37
logical xyzzyaaai37,xyzzyaaaj37
call timer('PUT_PARAMS',.true.)
xyzzyaaaj37=.false.
if(present(quirky))xyzzyaaaj37=quirky
xyzzyaaac37=0
do xyzzyaaad37=1,xyzzyaaah1
if(ignore(xyzzyaaad37))cycle
xyzzyaaac37=xyzzyaaac37+1
xyzzyaaag37=params(xyzzyaaac37)
xyzzyaaah37=xyzzyaaam1(xyzzyaaad37)
xyzzyaaag37=param_limiting(xyzzyaaan1.and.vm_smooth_limits,xyzzyaaag37&
&,xyzzyaaap1(xyzzyaaad37),xyzzyaaal1(xyzzyaaad37),xyzzyaaao1(xyzzyaaad&
&37),xyzzyaaak1(xyzzyaaad37))
xyzzyaaag37=xyzzyaaag37*xyzzyaaah37
xyzzyaaaf37(xyzzyaaad37)=xyzzyaaag37
enddo
xyzzyaaab37=0
if(xyzzyaaai1(1)>0)then
xyzzyaaaa37=xyzzyaaab37+1
xyzzyaaab37=xyzzyaaab37+xyzzyaaai1(1)
xyzzyaaae37=0
if(prestore.and.present(restore_hint))xyzzyaaae37=restore_hint-xyzzyaa&
&aa37+1
if(xyzzyaaae37==0.or.(xyzzyaaae37>0.and.xyzzyaaae37<=xyzzyaaai1(1)))th&
&en
select case(xyzzyaaae1)
case(1)
call put_slater_params(xyzzyaaaf37(xyzzyaaaa37:xyzzyaaab37),ignore(xyz&
&zyaaaa37:xyzzyaaab37),0,prestore.and..not.xyzzyaaaj37,xyzzyaaai37)
case(2)
call put_exmol_params(xyzzyaaaf37(xyzzyaaaa37:xyzzyaaab37),ignore(xyzz&
&yaaaa37:xyzzyaaab37),0,prestore.and..not.xyzzyaaaj37,xyzzyaaai37)
case(3)
call put_pfaff_params(xyzzyaaaf37(xyzzyaaaa37:xyzzyaaab37),ignore(xyzz&
&yaaaa37:xyzzyaaab37),0,prestore.and..not.xyzzyaaaj37,xyzzyaaai37)
case(4)
call put_geminal_params(xyzzyaaaf37(xyzzyaaaa37:xyzzyaaab37),ignore(xy&
&zzyaaaa37:xyzzyaaab37),0,prestore.and..not.xyzzyaaaj37,xyzzyaaai37)
case(6)
call put_mahan_params(xyzzyaaaf37(xyzzyaaaa37:xyzzyaaab37),ignore(xyzz&
&yaaaa37:xyzzyaaab37),0,prestore.and..not.xyzzyaaaj37,xyzzyaaai37)
case default
xyzzyaaai37=.false.
end select
if(present(bad_params))bad_params=bad_params.or.xyzzyaaai37
endif
endif
if(xyzzyaaai1(2)>0)then
xyzzyaaaa37=xyzzyaaab37+1
xyzzyaaab37=xyzzyaaab37+xyzzyaaai1(2)
xyzzyaaae37=0
if(prestore.and.present(restore_hint))xyzzyaaae37=restore_hint-xyzzyaa&
&aa37+1
if(xyzzyaaae37==0.or.(xyzzyaaae37>0.and.xyzzyaaae37<=xyzzyaaai1(2)))th&
&en
call put_jastrow_params(xyzzyaaaf37(xyzzyaaaa37:xyzzyaaab37),ignore(xy&
&zzyaaaa37:xyzzyaaab37),0,prestore,xyzzyaaai37,xyzzyaaae37,xyzzyaaaj37&
&)
if(present(bad_params))bad_params=bad_params.or.xyzzyaaai37
endif
endif
if(present(params_print))then
xyzzyaaac37=0
do xyzzyaaad37=1,xyzzyaaah1
if(ignore(xyzzyaaad37))cycle
xyzzyaaac37=xyzzyaaac37+1
params_print(xyzzyaaac37)=xyzzyaaaf37(xyzzyaaad37)
enddo
endif
call timer('PUT_PARAMS',.false.)
end subroutine put_params
subroutine put_param1(params_in,iparam_in,param,ignore,prestore)
use slaarnabt, only : param_limiting
implicit none
integer,intent(in) :: iparam_in
real(dp),intent(in) :: param,params_in(:)
logical,intent(in) :: ignore(xyzzyaaah1),prestore
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38,xyzzyaaad38,xyzzyaaae38,xy&
&zzyaaaf38,xyzzyaaag38
real(dp) xyzzyaaah38(xyzzyaaah1),xyzzyaaai38,xyzzyaaaj38
logical xyzzyaaak38(xyzzyaaah1),xyzzyaaal38
call timer('PUT_PARAM1',.true.)
xyzzyaaaf38=which_param(iparam_in,ignore)
xyzzyaaaa38=xyzzyaaaj1(xyzzyaaaf38)
xyzzyaaac38=sum(xyzzyaaai1(1:xyzzyaaaa38))
xyzzyaaab38=xyzzyaaac38-xyzzyaaai1(xyzzyaaaa38)+1
xyzzyaaag38=xyzzyaaaf38-xyzzyaaab38+1
xyzzyaaak38(xyzzyaaab38:xyzzyaaac38)=ignore(xyzzyaaab38:xyzzyaaac38).o&
&r..not.xyzzyaaar1(xyzzyaaab38:xyzzyaaac38,xyzzyaaaf38)
xyzzyaaah38(xyzzyaaab38:xyzzyaaac38)=0.d0
do xyzzyaaae38=xyzzyaaab38,xyzzyaaac38
if(xyzzyaaak38(xyzzyaaae38))cycle
xyzzyaaad38=xyzzyaaay1(xyzzyaaae38,ignore)
xyzzyaaai38=params_in(xyzzyaaad38)
if(xyzzyaaad38==iparam_in)xyzzyaaai38=param
xyzzyaaaj38=xyzzyaaam1(xyzzyaaae38)
xyzzyaaai38=param_limiting(xyzzyaaan1.and.vm_smooth_limits,xyzzyaaai38&
&,xyzzyaaap1(xyzzyaaae38),xyzzyaaal1(xyzzyaaae38),xyzzyaaao1(xyzzyaaae&
&38),xyzzyaaak1(xyzzyaaae38))
xyzzyaaai38=xyzzyaaai38*xyzzyaaaj38
xyzzyaaah38(xyzzyaaae38)=xyzzyaaai38
enddo
select case(xyzzyaaaa38)
case(1)
select case(xyzzyaaae1)
case(1)
call put_slater_params(xyzzyaaah38(xyzzyaaab38:xyzzyaaac38),xyzzyaaak3&
&8(xyzzyaaab38:xyzzyaaac38),xyzzyaaag38,prestore,xyzzyaaal38)
case(2)
call put_exmol_params(xyzzyaaah38(xyzzyaaab38:xyzzyaaac38),xyzzyaaak38&
&(xyzzyaaab38:xyzzyaaac38),xyzzyaaag38,prestore,xyzzyaaal38)
case(3)
call put_pfaff_params(xyzzyaaah38(xyzzyaaab38:xyzzyaaac38),xyzzyaaak38&
&(xyzzyaaab38:xyzzyaaac38),xyzzyaaag38,prestore,xyzzyaaal38)
case(4)
call put_geminal_params(xyzzyaaah38(xyzzyaaab38:xyzzyaaac38),xyzzyaaak&
&38(xyzzyaaab38:xyzzyaaac38),xyzzyaaag38,prestore,xyzzyaaal38)
case(6)
call put_mahan_params(xyzzyaaah38(xyzzyaaab38:xyzzyaaac38),xyzzyaaak38&
&(xyzzyaaab38:xyzzyaaac38),xyzzyaaag38,prestore,xyzzyaaal38)
end select
case(2)
call put_jastrow_params(xyzzyaaah38(xyzzyaaab38:xyzzyaaac38),xyzzyaaak&
&38(xyzzyaaab38:xyzzyaaac38),xyzzyaaag38,prestore,xyzzyaaal38,0,.false&
&.)
end select
call timer('PUT_PARAM1',.false.)
end subroutine put_param1
subroutine invalidate_param1(is,iparam_in,ignore)
implicit none
integer,intent(in) :: is,iparam_in
logical,intent(in) :: ignore(xyzzyaaah1)
integer xyzzyaaaa39,xyzzyaaab39,xyzzyaaac39
call timer('INVALIDATE_PARAM1',.true.)
xyzzyaaac39=which_param(iparam_in,ignore)
xyzzyaaaa39=xyzzyaaaj1(xyzzyaaac39)
xyzzyaaab39=xyzzyaaac39-sum(xyzzyaaai1(1:xyzzyaaaa39-1))
select case(xyzzyaaaa39)
case(1)
select case(xyzzyaaae1)
case(1)
call invalidate_param1_slater(is,xyzzyaaab39)
case(2)
call invalidate_param1_exmol(is,xyzzyaaab39)
case(3)
call invalidate_param1_pfaff(is,xyzzyaaab39)
case(4)
call invalidate_param1_geminal(is,xyzzyaaab39)
case(6)
call invalidate_param1_mahan(is,xyzzyaaab39)
end select
case(2)
call invalidate_param1_jastrow(is,xyzzyaaab39)
end select
call timer('INVALIDATE_PARAM1',.false.)
end subroutine invalidate_param1
subroutine invalidate_params(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaaah1)
integer xyzzyaaaa40,xyzzyaaab40
call timer('INVALIDATE_PARAMS',.true.)
xyzzyaaab40=0
if(xyzzyaaai1(1)>0)then
xyzzyaaaa40=xyzzyaaab40+1
xyzzyaaab40=xyzzyaaab40+xyzzyaaai1(1)
select case(xyzzyaaae1)
case(1)
call invalidate_params_slater(ignore(xyzzyaaaa40:xyzzyaaab40))
case(2)
call invalidate_params_exmol(ignore(xyzzyaaaa40:xyzzyaaab40))
case(3)
call invalidate_params_pfaff(ignore(xyzzyaaaa40:xyzzyaaab40))
case(4)
call invalidate_params_geminal(ignore(xyzzyaaaa40:xyzzyaaab40))
case(6)
call invalidate_params_mahan(ignore(xyzzyaaaa40:xyzzyaaab40))
end select
endif
if(xyzzyaaai1(2)>0)then
xyzzyaaaa40=xyzzyaaab40+1
xyzzyaaab40=xyzzyaaab40+xyzzyaaai1(2)
call invalidate_params_jastrow(ignore(xyzzyaaaa40:xyzzyaaab40))
endif
call timer('INVALIDATE_PARAMS',.false.)
end subroutine invalidate_params
integer function which_param(iparam_in,ignore)
implicit none
integer,intent(in) :: iparam_in
logical,intent(in) :: ignore(xyzzyaaah1)
integer xyzzyaaaa41,xyzzyaaab41
xyzzyaaaa41=0
which_param=0
do xyzzyaaab41=1,xyzzyaaah1
if(ignore(xyzzyaaab41))cycle
xyzzyaaaa41=xyzzyaaaa41+1
if(xyzzyaaaa41==iparam_in)then
which_param=xyzzyaaab41
return
endif
enddo
end function which_param
integer function xyzzyaaay1(iparam,ignore)
implicit none
integer,intent(in) :: iparam
logical,intent(in) :: ignore(xyzzyaaah1)
if(iparam<1.or.iparam>xyzzyaaah1)then
xyzzyaaay1=0
elseif(ignore(iparam))then
xyzzyaaay1=0
else
xyzzyaaay1=count(.not.ignore(1:iparam))
endif
end function xyzzyaaay1
logical function param_has_aderiv(iparam_in,ignore)
implicit none
integer,intent(in) :: iparam_in
logical,intent(in) :: ignore(xyzzyaaah1)
param_has_aderiv=xyzzyaaaq1(which_param(iparam_in,ignore))
end function param_has_aderiv
logical function any_analytical_deriv(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaaah1)
any_analytical_deriv=any(.not.ignore.and.xyzzyaaaq1)
end function any_analytical_deriv
logical function any_numerical_deriv(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaaah1)
any_numerical_deriv=any(.not.ignore.and..not.xyzzyaaaq1)
end function any_numerical_deriv
subroutine update_wfn_casl
implicit none
if(use_jastrow.and.use_gjastrow)call update_jastrow_casl
if(xyzzyaaae1==4)call update_geminal_casl
end subroutine update_wfn_casl
subroutine wfn_aderiv(is,iparam_in,ignore,dloggrad_psi,dloglap_psi,dlo&
&gval_psi)
implicit none
integer,intent(in) :: is,iparam_in
complex(dp),intent(inout) :: dloggrad_psi(3,netot),dloglap_psi(netot)
complex(dp),intent(inout),optional :: dlogval_psi
logical,intent(in) :: ignore(xyzzyaaah1)
integer xyzzyaaaa47,xyzzyaaab47,xyzzyaaac47,xyzzyaaad47,xyzzyaaae47
real(dp) xyzzyaaaf47
dloggrad_psi=czero
dloglap_psi=czero
if(present(dlogval_psi))dlogval_psi=czero
xyzzyaaaa47=which_param(iparam_in,ignore)
xyzzyaaab47=xyzzyaaaj1(xyzzyaaaa47)
xyzzyaaad47=sum(xyzzyaaai1(1:xyzzyaaab47))
xyzzyaaac47=xyzzyaaad47-xyzzyaaai1(xyzzyaaab47)+1
xyzzyaaae47=xyzzyaaaa47-xyzzyaaac47+1
select case(xyzzyaaab47)
case(1)
select case(xyzzyaaae1)
case(1)
call errstop_master('WFN_ADERIV','Unsupported for Slater dets.')
case(2)
call errstop_master('WFN_ADERIV','Unsupported for EXMOL.')
case(3)
call errstop_master('WFN_ADERIV','Unsupported for Pfaffians.')
case(4)
call errstop_master('WFN_ADERIV','Unsupported for Geminals.')
case(6)
call errstop_master('WFN_ADERIV','Unsupported for Mahan excitons.')
end select
case(2)
call wfn_aderiv_jastrow(is,xyzzyaaae47,dloggrad_psi,dloglap_psi,dlogva&
&l_psi)
end select
xyzzyaaaf47=xyzzyaaam1(iparam_in)
if(present(dlogval_psi))dlogval_psi=dlogval_psi*xyzzyaaaf47
dloggrad_psi=dloggrad_psi*xyzzyaaaf47
dloglap_psi=dloglap_psi*xyzzyaaaf47
end subroutine wfn_aderiv
subroutine setup_storage(nconfig,ignore)
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(xyzzyaaah1)
integer xyzzyaaaa48,xyzzyaaab48
logical xyzzyaaac48(0)
call setup_storage_geom(nconfig)
xyzzyaaab48=0
xyzzyaaaa48=xyzzyaaab48+1
xyzzyaaab48=xyzzyaaab48+xyzzyaaai1(1)
if(xyzzyaaab48<xyzzyaaaa48)then
select case(xyzzyaaae1)
case(1)
call setup_storage_slater(nconfig,xyzzyaaac48)
case(2)
call setup_storage_exmol(nconfig,xyzzyaaac48)
case(3)
call setup_storage_pfaff(nconfig,xyzzyaaac48)
case(4)
call setup_storage_geminal(nconfig,xyzzyaaac48)
case(6)
call setup_storage_mahan(nconfig,xyzzyaaac48)
end select
else
select case(xyzzyaaae1)
case(1)
call setup_storage_slater(nconfig,ignore(xyzzyaaaa48:xyzzyaaab48))
case(2)
call setup_storage_exmol(nconfig,ignore(xyzzyaaaa48:xyzzyaaab48))
case(3)
call setup_storage_pfaff(nconfig,ignore(xyzzyaaaa48:xyzzyaaab48))
case(4)
call setup_storage_geminal(nconfig,ignore(xyzzyaaaa48:xyzzyaaab48))
case(6)
call setup_storage_mahan(nconfig,ignore(xyzzyaaaa48:xyzzyaaab48))
end select
endif
xyzzyaaaa48=xyzzyaaab48+1
xyzzyaaab48=xyzzyaaab48+xyzzyaaai1(2)
if(use_jastrow)then
if(xyzzyaaab48<xyzzyaaaa48)then
call setup_storage_jastrow(nconfig,xyzzyaaac48)
else
call setup_storage_jastrow(nconfig,ignore(xyzzyaaaa48:xyzzyaaab48))
endif
endif
end subroutine setup_storage
subroutine finish_storage
implicit none
call finish_storage_geom
select case(xyzzyaaae1)
case(1)
call finish_storage_slater
case(2)
call finish_storage_exmol
case(3)
call finish_storage_pfaff
case(4)
call finish_storage_geminal
case(6)
call finish_storage_mahan
end select
if(use_jastrow)call finish_storage_jastrow
end subroutine finish_storage
subroutine load_from_storage(is,icfg)
implicit none
integer,intent(in) :: is,icfg
call timer('LOAD_FROM_STORAGE',.true.)
call clear_scratch_wfn(is)
call xyzzyaaax1(is)
call load_from_storage_geom(is,icfg)
select case(xyzzyaaae1)
case(1)
call load_from_storage_slater(is,icfg)
case(2)
call load_from_storage_exmol(is,icfg)
case(3)
call load_from_storage_pfaff(is,icfg)
case(4)
call load_from_storage_geminal(is,icfg)
case(6)
call load_from_storage_mahan(is,icfg)
end select
if(use_jastrow)call load_from_storage_jastrow(is,icfg)
call timer('LOAD_FROM_STORAGE',.false.)
end subroutine load_from_storage
subroutine save_to_storage(is,icfg)
implicit none
integer,intent(in) :: is,icfg
call timer('SAVE_TO_STORAGE',.true.)
select case(xyzzyaaae1)
case(1)
call save_to_storage_slater(is,icfg)
case(2)
call save_to_storage_exmol(is,icfg)
case(3)
call save_to_storage_pfaff(is,icfg)
case(4)
call save_to_storage_geminal(is,icfg)
case(6)
call save_to_storage_mahan(is,icfg)
end select
if(use_jastrow)call save_to_storage_jastrow(is,icfg)
call timer('SAVE_TO_STORAGE',.false.)
end subroutine save_to_storage
subroutine enumerate_plot_wfn(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
integer xyzzyaaaa52,xyzzyaaab52,xyzzyaaac52,xyzzyaaad52
logical xyzzyaaae52,xyzzyaaaf52
xyzzyaaaf52=.not.(present(keyword).and.present(description))
if(xyzzyaaaf52)xyzzyaaau1=0
xyzzyaaab52=0
xyzzyaaat1=0
do xyzzyaaaa52=1,nthings_to_plot
xyzzyaaae52=.false.
select case(xyzzyaaaa52)
case(xyzzyaaas1)
xyzzyaaae52=.true.
end select
if(xyzzyaaae52)then
xyzzyaaab52=xyzzyaaab52+1
xyzzyaaat1(xyzzyaaab52)=xyzzyaaaa52
endif
enddo
xyzzyaaau1(0)=xyzzyaaab52
if(.not.xyzzyaaaf52)then
do xyzzyaaaa52=1,xyzzyaaab52
keyword(xyzzyaaaa52)=plot_keyword(xyzzyaaat1(xyzzyaaaa52))
description(xyzzyaaaa52)=plot_description(xyzzyaaat1(xyzzyaaaa52))
enddo
endif
if(xyzzyaaaf52)then
select case(xyzzyaaae1)
case(1)
call enumerate_plot_slater(xyzzyaaau1(1))
case(2)
call enumerate_plot_exmol(xyzzyaaau1(1))
case(3)
call enumerate_plot_pfaff(xyzzyaaau1(1))
case(4)
call enumerate_plot_geminal(xyzzyaaau1(1))
case(6)
call enumerate_plot_mahan(xyzzyaaau1(1))
end select
if(use_jastrow)call enumerate_plot_jastrow(xyzzyaaau1(2))
else
xyzzyaaaa52=xyzzyaaau1(0)+1
xyzzyaaab52=sum(xyzzyaaau1(0:1))
select case(xyzzyaaae1)
case(1)
call enumerate_plot_slater(xyzzyaaau1(1),keyword(xyzzyaaaa52:xyzzyaaab&
&52),description(xyzzyaaaa52:xyzzyaaab52))
case(2)
call enumerate_plot_exmol(xyzzyaaau1(1),keyword(xyzzyaaaa52:xyzzyaaab5&
&2),description(xyzzyaaaa52:xyzzyaaab52))
case(3)
call enumerate_plot_pfaff(xyzzyaaau1(1),keyword(xyzzyaaaa52:xyzzyaaab5&
&2),description(xyzzyaaaa52:xyzzyaaab52))
case(4)
call enumerate_plot_geminal(xyzzyaaau1(1),keyword(xyzzyaaaa52:xyzzyaaa&
&b52),description(xyzzyaaaa52:xyzzyaaab52))
case(6)
call enumerate_plot_mahan(xyzzyaaau1(1),keyword(xyzzyaaaa52:xyzzyaaab5&
&2),description(xyzzyaaaa52:xyzzyaaab52))
end select
xyzzyaaaa52=sum(xyzzyaaau1(0:1))+1
xyzzyaaab52=sum(xyzzyaaau1(0:2))
if(use_jastrow)call enumerate_plot_jastrow(xyzzyaaau1(2),keyword(xyzzy&
&aaaa52:xyzzyaaab52),description(xyzzyaaaa52:xyzzyaaab52))
endif
n=sum(xyzzyaaau1(0:2))
if(.not.xyzzyaaaf52)then
allocate(xyzzyaaav1(sum(xyzzyaaau1(0:2))),stat=xyzzyaaac52)
call check_alloc(xyzzyaaac52,'ENUMERATE_PLOT_WFN','which_plot_sec')
do xyzzyaaad52=0,2
xyzzyaaav1(sum(xyzzyaaau1(0:xyzzyaaad52-1))+1:sum(xyzzyaaau1(0:xyzzyaa&
&ad52)))=xyzzyaaad52
enddo
endif
end subroutine enumerate_plot_wfn
subroutine query_plot_wfn(iplot,ii,rank,is_complex,has_stderr,rot_tens&
&or,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
integer xyzzyaaaa53,iplot_rel
logical count_only
count_only=.not.present(function_name)
if(.not.allocated(xyzzyaaav1))call errstop_master('QUERY_PLOT_WFN','EN&
&UMERATE_PLOT_WFN not called yet. Bug in calling routine.')
if(iplot<0.or.iplot>sum(xyzzyaaau1(0:2)))call errstop_master('QUERY_PL&
&OT_WFN','IPLOT out of range. Bug in calling routine.')
xyzzyaaaa53=xyzzyaaav1(iplot)
iplot_rel=iplot-sum(xyzzyaaau1(0:xyzzyaaaa53-1))
select case(xyzzyaaaa53)
case(0)
select case(xyzzyaaat1(iplot_rel))
case(xyzzyaaas1)
rank=0
is_complex=complex_wf
has_stderr=.false.
rot_tensor=.false.
transl_pos=.false.
nfunctions=1
if(.not.count_only)then
function_name(nfunctions)='Wave function value'
endif
end select
case(1)
select case(xyzzyaaae1)
case(1)
call query_plot_slater(iplot_rel,ii,rank,is_complex,has_stderr,rot_ten&
&sor,transl_pos,nfunctions,function_name)
case(2)
call query_plot_exmol(iplot_rel,ii,rank,is_complex,has_stderr,rot_tens&
&or,transl_pos,nfunctions,function_name)
case(3)
call query_plot_pfaff(iplot_rel,ii,rank,is_complex,has_stderr,rot_tens&
&or,transl_pos,nfunctions,function_name)
case(4)
call query_plot_geminal(iplot_rel,ii,rank,is_complex,has_stderr,rot_te&
&nsor,transl_pos,nfunctions,function_name)
case(6)
call query_plot_mahan(iplot_rel,ii,rank,is_complex,has_stderr,rot_tens&
&or,transl_pos,nfunctions,function_name)
end select
case(2)
call query_plot_jastrow(iplot_rel,ii,rank,is_complex,has_stderr,rot_te&
&nsor,transl_pos,nfunctions,function_name)
end select
end subroutine query_plot_wfn
subroutine get_plot_wfn(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
integer xyzzyaaaa54,xyzzyaaab54
complex(dp) xyzzyaaac54
if(.not.allocated(xyzzyaaav1))call errstop_master('GET_PLOT_WFN','ENUM&
&ERATE_PLOT_WFN not called yet. Bug in calling routine.')
if(iplot<0.or.iplot>sum(xyzzyaaau1(0:2)))call errstop_master('GET_PLOT&
&_WFN','IPLOT out of range. Bug in calling routine.')
xyzzyaaaa54=xyzzyaaav1(iplot)
xyzzyaaab54=iplot-sum(xyzzyaaau1(0:xyzzyaaaa54-1))
select case(xyzzyaaaa54)
case(0)
select case(xyzzyaaat1(xyzzyaaab54))
case(xyzzyaaas1)
call wfn_ratio(is0,is1,0,ratio=xyzzyaaac54)
f(1)=dble(xyzzyaaac54)
if(complex_wf)f(2)=aimag(xyzzyaaac54)
end select
case(1)
select case(xyzzyaaae1)
case(1)
call get_plot_slater(xyzzyaaab54,ii,is0,is1,f)
case(2)
call get_plot_exmol(xyzzyaaab54,ii,is0,is1,f)
case(3)
call get_plot_pfaff(xyzzyaaab54,ii,is0,is1,f)
case(4)
call get_plot_geminal(xyzzyaaab54,ii,is0,is1,f)
case(6)
call get_plot_mahan(xyzzyaaab54,ii,is0,is1,f)
end select
case(2)
call get_plot_jastrow(xyzzyaaab54,ii,is0,is1,f)
end select
end subroutine get_plot_wfn
subroutine finish_plot_wfn
implicit none
xyzzyaaat1=0
xyzzyaaau1(0:2)=0
deallocate(xyzzyaaav1)
select case(xyzzyaaae1)
case(1)
call finish_plot_slater
case(2)
call finish_plot_exmol
case(3)
call finish_plot_pfaff
case(4)
call finish_plot_geminal
case(6)
call finish_plot_mahan
end select
if(use_jastrow)call finish_plot_jastrow
end subroutine finish_plot_wfn
real(dp) function get_wfn_rmax()
implicit none
select case(xyzzyaaae1)
case(1)
get_wfn_rmax=get_slater_rmax()
case(2)
get_wfn_rmax=0.d0
case(3)
get_wfn_rmax=0.d0
case(4)
get_wfn_rmax=0.d0
case(6)
get_wfn_rmax=0.d0
end select
end function get_wfn_rmax
subroutine wfn_assess_check_kinetic(is,ii,verbose)
implicit none
integer,intent(in) :: is
integer,intent(out) :: ii
logical,intent(in) :: verbose
integer xyzzyaaaa57,xyzzyaaab57
logical xyzzyaaac57
xyzzyaaac57=.true.
if(use_jastrow)then
call get_eevecs(is)
call get_eivecs(is)
call jastrow_assess_check_kinetic(eevecs_scr(1,1,1,is),eivecs_scr(1,1,&
&1,is),xyzzyaaaa57,xyzzyaaab57,xyzzyaaac57,verbose)
if(.not.xyzzyaaac57)ii=which_ii(xyzzyaaaa57,xyzzyaaab57)
endif
if(xyzzyaaac57.and.use_backflow)then
call get_eevecs(is)
call get_eivecs(is)
call bf_assess_check_kinetic(eevecs_scr(1,1,1,is),eivecs_scr(1,1,1,is)&
&,xyzzyaaaa57,xyzzyaaab57,xyzzyaaac57,verbose)
if(.not.xyzzyaaac57)ii=which_ii(xyzzyaaaa57,xyzzyaaab57)
endif
if(xyzzyaaac57)ii=0
end subroutine wfn_assess_check_kinetic
end module slaarnacs

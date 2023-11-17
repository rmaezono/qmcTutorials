module slaarnabn
use slaarnaad
use dsp
use slaarnach
use store
use slaarnaag,    only : czero
use file_utils,   only : open_units
use slaarnaas,    only : gs_kvec
use format_utils, only : wout,i2s,r2s,r2s2,write_list_int,wordwrap
use slaarnabg,     only : dimensionality,wigner_seitz_radius
use slaarnabt,    only : dcopy,ddot,lu_decom_cmplx,lu_logdet_cmplx,lu_&
&solve_n_cmplx,zdotu
use parallel,     only : am_master
use run_control,  only : errstop,errstop_master,check_alloc,timer
use slaarnabq,    only : minimum_image
implicit none
private
public enumerate_plot_mahan,query_plot_mahan,get_plot_mahan,finish_plo&
&t_mahan
public query_mahan_levels,query_mahan_level_details,setup_mahan,finish&
&_mahan,wfn_ratio_mahan,accept_move_mahan,reset_config_mahan,wfn_logva&
&l_mahan,wfn_loggrad_mahan,wfn_loglap_mahan,prefetch_wfn_mahan,clear_s&
&cratch_mahan,add_config_mahan_items,setup_mahan_params,finish_mahan_p&
&arams,get_mahan_params,put_mahan_params,clone_scratch_mahan,invalidat&
&e_params_mahan,invalidate_param1_mahan,setup_storage_mahan,finish_sto&
&rage_mahan,load_from_storage_mahan,save_to_storage_mahan,read_mahan,w&
&rite_mahan
public gen_config_mahan,delete_config_mahan,copy_config_mahan,config_t&
&o_pt_mahan,pt_to_config_mahan,redist_allocations_mahan,redist_load_ma&
&han,redist_send_mahan,redist_recv_mahan,redist_save_mahan,redist_deal&
&locations_mahan,load_from_pt_mahan,save_to_pt_mahan
public config_wfn_mahan
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1
complex(dp),allocatable :: xyzzyaaae1(:,:,:,:),xyzzyaaaf1(:,:,:,:),xyz&
&zyaaag1(:,:,:,:,:),xyzzyaaah1(:,:,:,:),xyzzyaaai1(:,:,:,:,:),xyzzyaaa&
&j1(:)
real(dp),allocatable :: xyzzyaaak1(:,:,:,:),xyzzyaaal1(:,:,:,:,:,:),xy&
&zzyaaam1(:,:,:,:,:),xyzzyaaan1(:,:,:,:,:)
complex(dp),allocatable :: xyzzyaaao1(:,:,:),xyzzyaaap1(:,:)
logical,allocatable :: xyzzyaaaq1(:),xyzzyaaar1(:),xyzzyaaas1(:),xyzzy&
&aaat1(:),xyzzyaaau1(:),xyzzyaaav1(:,:),xyzzyaaaw1(:,:),xyzzyaaax1(:),&
&xyzzyaaay1(:),xyzzyaaaz1(:),xyzzyaaba1(:),xyzzyaabb1(:)
integer xyzzyaabc1,xyzzyaabd1
integer,parameter :: xyzzyaabe1=2
integer xyzzyaabf1,xyzzyaabg1(xyzzyaabe1)
integer,allocatable :: xyzzyaabh1(:),xyzzyaabi1(:,:),xyzzyaabj1(:,:)
real(dp),allocatable :: xyzzyaabk1(:,:),xyzzyaabl1(:,:)
real(dp),allocatable :: xyzzyaabm1(:,:,:),xyzzyaabn1(:,:,:)
integer xyzzyaabo1,xyzzyaabp1
logical xyzzyaabq1
integer xyzzyaabr1
integer xyzzyaabs1
character(80) title
type config_wfn_mahan
private
logical dummy_variable
end type config_wfn_mahan
integer,allocatable :: xyzzyaabt1(:,:)
integer xyzzyaabu1,xyzzyaabv1
integer xyzzyaabw1,xyzzyaabx1
integer,allocatable :: xyzzyaaby1(:,:),xyzzyaabz1(:),xyzzyaaca1(:),xyz&
&zyaacb1(:),xyzzyaacc1(:)
real(dp),allocatable :: xyzzyaacd1(:),xyzzyaace1(:,:),xyzzyaacf1(:,:)
complex(dp),parameter :: xyzzyaacg1=(-500.d0,0.d0),xyzzyaach1=(-1000.d&
&0,0.d0)
real(dp),parameter :: xyzzyaaci1=-690.d0,xyzzyaacj1=-20.d0
integer,parameter :: xyzzyaack1=0,xyzzyaacl1=1,xyzzyaacm1=2
integer xyzzyaacn1
logical :: xyzzyaaco1=.false.
contains
subroutine query_mahan_levels(nlevels)
implicit none
integer,intent(out) :: nlevels
nlevels=1
end subroutine query_mahan_levels
subroutine query_mahan_level_details(level_score,level_name)
implicit none
integer,intent(out) :: level_score(*)
character(80),intent(out) :: level_name(*)
level_score(1)=200
level_name(1)='mahan_ex'
end subroutine query_mahan_level_details
subroutine setup_mahan
implicit none
integer xyzzyaaaa4
integer xyzzyaaab4(2),xyzzyaaac4(2),xyzzyaaad4(2),xyzzyaaae4(2),xyzzya&
&aaf4(2),xyzzyaaag4(2),xyzzyaaah4(2),xyzzyaaai4(2),xyzzyaaaj4(2),xyzzy&
&aaak4(2),xyzzyaaal4(2),xyzzyaaam4(2),xyzzyaaan4
if(dimensionality==3)then
xyzzyaabu1=6
elseif(dimensionality==2)then
xyzzyaabu1=3
else
call errstop('SETUP_MAHAN_EX','dimensionality not implemented')
endif
xyzzyaaaa1=nemax*nemax*xyzzyaacn1
xyzzyaaab1=dimensionality*xyzzyaaaa1
xyzzyaaac1=xyzzyaabu1*xyzzyaaaa1
xyzzyaaad1=3*xyzzyaabv1*nemax*xyzzyaacn1
xyzzyaaab4=0
xyzzyaaac4=0
xyzzyaaad4=0
xyzzyaaae4=0
xyzzyaaaf4=0
xyzzyaaag4=0
xyzzyaaah4=0
xyzzyaaai4=0
xyzzyaaaj4=0
xyzzyaaak4=0
xyzzyaaal4=0
xyzzyaaam4=0
call include_range((/1,nscratch/),xyzzyaaab4)
call include_range((/1,nscratch/),xyzzyaaac4)
call include_range((/1,nscratch/),xyzzyaaad4)
call include_range((/1,nscratch/),xyzzyaaak4)
call include_range((/1,nscratch/),xyzzyaaag4)
call include_range((/1,nscratch/),xyzzyaaah4)
call include_range((/1,nscratch/),xyzzyaaae4)
call include_range((/1,nscratch/),xyzzyaaaf4)
call include_range((/1,nscratch/),xyzzyaaai4)
call include_range((/1,nscratch/),xyzzyaaaj4)
call include_range((/1,nscratch/),xyzzyaaal4)
call include_range((/1,nscratch/),xyzzyaaam4)
if(xyzzyaaab4(1)/=0)then
allocate(xyzzyaaae1(nemax,nemax,xyzzyaacn1,xyzzyaaab4(1):xyzzyaaab4(2)&
&),xyzzyaaaf1(nemax,nemax,xyzzyaacn1,xyzzyaaab4(1):xyzzyaaab4(2)),stat&
&=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','orbmat')
xyzzyaaae1=0.d0
xyzzyaaaf1=0.d0
endif
allocate(xyzzyaaaq1(nscratch),xyzzyaaar1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','orbmat_valid')
if(xyzzyaaad4(1)/=0)then
allocate(xyzzyaaag1(dimensionality,nemax,nemax,xyzzyaacn1,xyzzyaaad4(1&
&):xyzzyaaad4(2)),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','gradmat')
xyzzyaaag1=0.d0
endif
allocate(xyzzyaaas1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','gradmat_valid')
if(xyzzyaaae4(1)/=0)then
allocate(xyzzyaaah1(nemax,nemax,xyzzyaacn1,xyzzyaaae4(1):xyzzyaaae4(2)&
&),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','lapmat')
xyzzyaaah1=0.d0
endif
allocate(xyzzyaaat1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','lapmat_valid')
if(xyzzyaaaf4(1)/=0)then
allocate(xyzzyaaai1(xyzzyaabu1,nemax,nemax,xyzzyaacn1,xyzzyaaaf4(1):xy&
&zzyaaaf4(2)),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','d2mat')
xyzzyaaai1=0.d0
endif
allocate(xyzzyaaau1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','d2mat_valid')
if(xyzzyaaak4(1)/=0)then
allocate(xyzzyaaaj1(xyzzyaaak4(1):xyzzyaaak4(2)),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','logdet_scr')
endif
allocate(xyzzyaaaz1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','logdet_valid')
if(xyzzyaaag4(1)/=0)then
allocate(xyzzyaaao1(3,netot,xyzzyaaag4(1):xyzzyaaag4(2)),stat=xyzzyaaa&
&n4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','fi')
xyzzyaaao1=czero
endif
allocate(xyzzyaaav1(netot,nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','fi_valid')
if(xyzzyaaah4(1)/=0)then
allocate(xyzzyaaap1(netot,xyzzyaaah4(1):xyzzyaaah4(2)),stat=xyzzyaaan4&
&)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','ti')
xyzzyaaap1=czero
endif
allocate(xyzzyaaaw1(netot,nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','ti_valid')
if(xyzzyaaai4(1)/=0)then
allocate(xyzzyaaak1(3,real1_complex2,netot,xyzzyaaai4(1):xyzzyaaai4(2)&
&),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','Farray')
xyzzyaaak1=0.d0
endif
allocate(xyzzyaaax1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','Farray_valid')
if(xyzzyaaaj4(1)/=0)then
allocate(xyzzyaaal1(3,3,real1_complex2,netot,netot,xyzzyaaaj4(1):xyzzy&
&aaaj4(2)),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','Harray')
endif
allocate(xyzzyaaay1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','Harray_valid')
if(xyzzyaaal4(1)/=0)then
allocate(xyzzyaaam1(3,xyzzyaabv1,nemax,xyzzyaacn1,xyzzyaaal4(1):xyzzya&
&aal4(2)),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','eta_G')
xyzzyaaam1=0.d0
endif
allocate(xyzzyaaba1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','eta_G_valid')
if(xyzzyaaam4(1)/=0)then
allocate(xyzzyaaan1(3,xyzzyaabv1,nemax,xyzzyaacn1,xyzzyaaam4(1):xyzzya&
&aam4(2)),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','u_G')
xyzzyaaan1=0.d0
endif
allocate(xyzzyaabb1(nscratch),stat=xyzzyaaan4)
call check_alloc(xyzzyaaan4,'SETUP_MAHAN','u_G_valid')
if(use_backflow)call setup_bf
do xyzzyaaaa4=1,nscratch
call clear_scratch_mahan(xyzzyaaaa4)
enddo
end subroutine setup_mahan
subroutine finish_mahan
implicit none
if(allocated(xyzzyaaae1))deallocate(xyzzyaaae1)
deallocate(xyzzyaaaq1)
if(allocated(xyzzyaaaf1))deallocate(xyzzyaaaf1)
deallocate(xyzzyaaar1)
if(allocated(xyzzyaaag1))deallocate(xyzzyaaag1)
deallocate(xyzzyaaas1)
if(allocated(xyzzyaaah1))deallocate(xyzzyaaah1)
deallocate(xyzzyaaat1)
if(allocated(xyzzyaaai1))deallocate(xyzzyaaai1)
deallocate(xyzzyaaau1)
if(allocated(xyzzyaaaj1))deallocate(xyzzyaaaj1)
deallocate(xyzzyaaaz1)
if(allocated(xyzzyaaao1))deallocate(xyzzyaaao1)
deallocate(xyzzyaaav1)
if(allocated(xyzzyaaap1))deallocate(xyzzyaaap1)
deallocate(xyzzyaaaw1)
if(allocated(xyzzyaaak1))deallocate(xyzzyaaak1)
deallocate(xyzzyaaax1)
if(allocated(xyzzyaaal1))deallocate(xyzzyaaal1)
deallocate(xyzzyaaay1)
if(allocated(xyzzyaaam1))deallocate(xyzzyaaam1)
deallocate(xyzzyaaba1)
if(allocated(xyzzyaaan1))deallocate(xyzzyaaan1)
deallocate(xyzzyaabb1)
if(use_backflow)call finish_bf
end subroutine finish_mahan
subroutine wfn_ratio_mahan(is,js,ilevel,ratio,fd,sd,isnan,isinf)
implicit none
integer,intent(in) :: is,js,ilevel
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: ratio
logical,intent(out) :: isnan,isinf
integer xyzzyaaaa6,xyzzyaaab6
complex(dp) xyzzyaaac6,xyzzyaaad6
logical xyzzyaaae6
xyzzyaaae6=.false.
xyzzyaaaa6=buffer_move1_from(is)
if(xyzzyaaaa6/=0)then
xyzzyaaab6=buffer_move1_from_ii(is)
xyzzyaaae6=which_spin(xyzzyaaab6)/=3
endif
xyzzyaaae6=.false.
if(xyzzyaaae6)then
else
call xyzzyaacp1(is,.true.,fd,sd)
call xyzzyaacq1(is,fd,sd,xyzzyaaac6)
call xyzzyaacp1(js,.true.,fd,sd)
call xyzzyaacq1(js,fd,sd,xyzzyaaad6)
ratio=exp(xyzzyaaad6-xyzzyaaac6)
endif
if(xyzzyaaco1)then
call wout("Old logdet is: ",xyzzyaaac6)
call wout("New logdet is: ",xyzzyaaad6)
call wout("wfn ratio is: ",ratio)
endif
isnan=.false.
isinf=.false.
end subroutine wfn_ratio_mahan
subroutine accept_move_mahan(is,js)
implicit none
integer,intent(in) :: is,js
if(xyzzyaaaq1(js))then
call zcopy(xyzzyaaaa1,xyzzyaaae1(1,1,1,js),1,xyzzyaaae1(1,1,1,is),1)
xyzzyaaaq1(is)=.true.
else
xyzzyaaaq1(is)=.false.
endif
if(xyzzyaaar1(js))then
call zcopy(xyzzyaaaa1,xyzzyaaaf1(1,1,1,js),1,xyzzyaaaf1(1,1,1,is),1)
xyzzyaaar1(is)=.true.
else
xyzzyaaar1(is)=.false.
endif
if(xyzzyaaas1(js))then
call zcopy(xyzzyaaab1,xyzzyaaag1(1,1,1,1,js),1,xyzzyaaag1(1,1,1,1,is),&
&1)
xyzzyaaas1(is)=.true.
else
xyzzyaaas1(is)=.false.
endif
if(xyzzyaaat1(js))then
call zcopy(xyzzyaaaa1,xyzzyaaah1(1,1,1,js),1,xyzzyaaah1(1,1,1,is),1)
xyzzyaaat1(is)=.true.
else
xyzzyaaat1(is)=.false.
endif
if(xyzzyaaau1(js))then
call zcopy(xyzzyaaac1,xyzzyaaai1(1,1,1,1,js),1,xyzzyaaai1(1,1,1,1,is),&
&1)
xyzzyaaau1(is)=.true.
else
xyzzyaaau1(is)=.false.
endif
if(xyzzyaaaz1(js))then
xyzzyaaaj1(is)=xyzzyaaaj1(js)
xyzzyaaaz1(is)=.true.
else
xyzzyaaaz1(is)=.false.
endif
if(xyzzyaaax1(js))then
call dcopy(size_farray,xyzzyaaak1(1,1,1,js),1,xyzzyaaak1(1,1,1,is),1)
xyzzyaaax1(is)=.true.
else
xyzzyaaax1(is)=.false.
endif
if(xyzzyaaay1(js))then
call dcopy(size_harray,xyzzyaaal1(1,1,1,1,1,js),1,xyzzyaaal1(1,1,1,1,1&
&,is),1)
xyzzyaaay1(is)=.true.
else
xyzzyaaay1(is)=.false.
endif
if(xyzzyaaba1(js))then
call dcopy(xyzzyaaad1,xyzzyaaam1(1,1,1,1,js),1,xyzzyaaam1(1,1,1,1,is),&
&1)
xyzzyaaba1(is)=.true.
else
xyzzyaaba1(is)=.false.
endif
if(xyzzyaabb1(js))then
call dcopy(xyzzyaaad1,xyzzyaaan1(1,1,1,1,js),1,xyzzyaaan1(1,1,1,1,is),&
&1)
xyzzyaabb1(is)=.true.
else
xyzzyaabb1(is)=.false.
endif
if(use_backflow)call accept_move_bf(is,js)
xyzzyaaav1(:,is)=.false.
xyzzyaaaw1(:,is)=.false.
end subroutine accept_move_mahan
subroutine reset_config_mahan(is,js)
implicit none
integer,intent(in) :: is,js
call clear_scratch_mahan(js)
end subroutine reset_config_mahan
subroutine wfn_logval_mahan(is,logwfn,iszero)
implicit none
integer,intent(in) :: is
complex(dp),intent(out) :: logwfn
logical,intent(out) :: iszero
complex(dp) xyzzyaaaa9
call xyzzyaacp1(is,.true.,.false.,.false.)
call xyzzyaacq1(is,.false.,.false.,xyzzyaaaa9)
logwfn=xyzzyaaaa9
iszero=.false.
end subroutine wfn_logval_mahan
subroutine wfn_loggrad_mahan(ii,is,ilevel,val,sd,loggrad,isnan,isinf)
implicit none
integer,intent(in) :: ii,is,ilevel
logical,intent(in) :: val,sd
complex(dp),intent(out) :: loggrad(3)
logical,intent(out) :: isnan,isinf
complex(dp) xyzzyaaaa10
call xyzzyaacp1(is,.true.,.true.,sd)
call xyzzyaacq1(is,.true.,sd,xyzzyaaaa10)
call xyzzyaacr1(is,ii,sd,xyzzyaaaa10,loggrad)
if(xyzzyaaco1)then
call wout("loggrad is: ", loggrad)
endif
isnan=.false.
isinf=.false.
end subroutine wfn_loggrad_mahan
subroutine wfn_loglap_mahan(ii,is,ilevel,val,fd,loglap,isnan,isinf)
implicit none
integer,intent(in) :: ii,is,ilevel
logical,intent(in) :: val,fd
complex(dp),intent(out) :: loglap
logical,intent(out) :: isnan,isinf
complex(dp) xyzzyaaaa11
complex(dp) xyzzyaaab11(3),xyzzyaaac11
call xyzzyaacp1(is,.true.,.true.,.true.)
call xyzzyaacq1(is,.true.,.true.,xyzzyaaaa11)
call xyzzyaacr1(is,ii,.true.,xyzzyaaaa11,xyzzyaaab11)
call xyzzyaacs1(is,ii,xyzzyaaab11,loglap)
if(xyzzyaaco1)then
call wout('Note this Laplacian info is not correct for cmplx wfns:')
call wout('lap1 is: ',-zdotu(3,xyzzyaaab11,1,xyzzyaaab11,1))
call wout('lap2 is: ',xyzzyaaac11)
call wout('loglap is: ',loglap)
endif
isnan=.false.
isinf=.false.
end subroutine wfn_loglap_mahan
subroutine prefetch_wfn_mahan(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
end subroutine prefetch_wfn_mahan
subroutine xyzzyaacp1(is,val,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa13,xyzzyaaab13
real(dp) xyzzyaaac13(dimensionality,nemax,xyzzyaacn1)
logical xyzzyaaad13,xyzzyaaae13,xyzzyaaaf13,xyzzyaaag13
xyzzyaaad13=val.and..not.xyzzyaaaq1(is)
xyzzyaaae13=fd.and..not.xyzzyaaas1(is)
if(.not.use_backflow)then
xyzzyaaaf13=sd.and..not.xyzzyaaat1(is)
else
xyzzyaaaf13=sd.and..not.xyzzyaaau1(is)
endif
if(.not.(xyzzyaaad13.or.xyzzyaaae13.or.xyzzyaaaf13))return
call timer('GET_ORBMAT',.true.)
xyzzyaaag13=.false.
xyzzyaaaa13=buffer_move1_from(is)
if(xyzzyaaaa13/=0)then
xyzzyaaag13=.true.
if(xyzzyaaad13.and..not.xyzzyaaaq1(xyzzyaaaa13))xyzzyaaag13=.false.
if(xyzzyaaae13.and..not.xyzzyaaas1(xyzzyaaaa13))xyzzyaaag13=.false.
if(xyzzyaaaf13.and..not.xyzzyaaat1(xyzzyaaaa13))xyzzyaaag13=.false.
endif
if(xyzzyaaag13.and..false.)then
call errstop("Mahan_ex","Shortcut is .true.  Not implemented.")
else
call xyzzyaacv1(is,xyzzyaaac13,xyzzyaaae13.or.xyzzyaaaf13,xyzzyaaaf13)
call xyzzyaacw1(is,xyzzyaaac13,xyzzyaaae13.or.xyzzyaaaf13,xyzzyaaaf13)
call xyzzyaacx1(is,xyzzyaaac13,xyzzyaaae13.or.xyzzyaaaf13,xyzzyaaaf13)
do xyzzyaaab13=1,xyzzyaacn1
if(heg_nele(xyzzyaaab13)==0)cycle
call xyzzyaade1(heg_nele(xyzzyaaab13),xyzzyaaac13(1,1,xyzzyaaab13),xyz&
&zyaaam1(1,1,1,xyzzyaaab13,is),xyzzyaaan1(1,1,1,xyzzyaaab13,is),xyzzya&
&aad13,xyzzyaaae13,xyzzyaaaf13,xyzzyaaae1(1,1,xyzzyaaab13,is),xyzzyaaa&
&g1(1,1,1,xyzzyaaab13,is),xyzzyaaah1(1,1,xyzzyaaab13,is),xyzzyaaai1(1,&
&1,1,xyzzyaaab13,is))
enddo
if(xyzzyaaad13)xyzzyaaaq1(is)=.true.
if(xyzzyaaae13)xyzzyaaas1(is)=.true.
if(xyzzyaaaf13)xyzzyaaat1(is)=.true.
if(xyzzyaaaf13)xyzzyaaau1(is)=.true.
endif
call timer('GET_ORBMAT',.false.)
end subroutine xyzzyaacp1
subroutine xyzzyaacq1(is,fd,sd,logdet)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: logdet
integer xyzzyaaaa14,xyzzyaaab14
complex(dp) xyzzyaaac14
logical xyzzyaaad14
if(xyzzyaaaz1(is).and.xyzzyaaar1(is))then
logdet=xyzzyaaaj1(is)
return
endif
xyzzyaaad14=fd.or.sd
xyzzyaaad14=xyzzyaaad14.and..not.xyzzyaaar1(is)
call timer('GET_CALDET',.true.)
logdet=czero
do xyzzyaaab14=1,xyzzyaacn1
if(heg_nele(xyzzyaaab14)==0)cycle
call xyzzyaadc1(heg_nele(xyzzyaaab14),xyzzyaaad14,xyzzyaaac14,xyzzyaaa&
&e1(1,1,xyzzyaaab14,is),xyzzyaaaf1(1,1,xyzzyaaab14,is),xyzzyaaaa14)
logdet=logdet+xyzzyaaac14
enddo
xyzzyaaaj1(is)=logdet
xyzzyaaaz1(is)=.true.
if(xyzzyaaad14)xyzzyaaar1(is)=.true.
call timer('GET_CALDET',.false.)
end subroutine xyzzyaacq1
subroutine xyzzyaacr1(is,ii,sd,logdet,loggrad)
implicit none
integer,intent(in) :: is,ii
complex(dp),intent(in) :: logdet
logical,intent(in) :: sd
complex(dp),intent(out) :: loggrad(3)
integer xyzzyaaaa15,xyzzyaaab15
complex(dp) xyzzyaaac15(dimensionality)
if(.not.use_backflow)then
xyzzyaaab15=which_ie(ii)
xyzzyaaaa15=which_spin(ii)
if(xyzzyaaaa15==3)then
call xyzzyaadg1(logdet,xyzzyaaag1(:,:,:,:,is),xyzzyaaaf1(:,:,:,is),xyz&
&zyaaac15)
else
call xyzzyaadf1(heg_nele(xyzzyaaaa15),logdet,xyzzyaaag1(:,:,xyzzyaaab1&
&5,xyzzyaaaa15,is),xyzzyaaaf1(:,xyzzyaaab15,xyzzyaaaa15,is),xyzzyaaac1&
&5)
endif
loggrad=czero
loggrad(1:dimensionality)=xyzzyaaac15
else
call xyzzyaact1(is)
if(sd)call xyzzyaacu1(is)
call loggrad_bf(ii,is,sd,xyzzyaaak1(1,1,1,is),loggrad)
endif
end subroutine xyzzyaacr1
subroutine xyzzyaacs1(is,ii,loggrad,lap)
implicit none
integer,intent(in) :: is,ii
complex(dp),intent(in) :: loggrad(3)
complex(dp),intent(out) :: lap
integer xyzzyaaaa16,xyzzyaaab16
xyzzyaaab16=which_ie(ii)
xyzzyaaaa16=which_spin(ii)
if(.not.use_backflow)then
if(xyzzyaaaa16==3)then
call xyzzyaadi1(xyzzyaaag1(:,:,:,:,is),xyzzyaaaf1(:,:,:,is),xyzzyaaah1&
&(:,:,:,is),loggrad,lap)
else
call xyzzyaadh1(heg_nele(xyzzyaaaa16),xyzzyaaah1(:,xyzzyaaab16,xyzzyaa&
&aa16,is),xyzzyaaaf1(:,xyzzyaaab16,xyzzyaaaa16,is),loggrad,lap)
endif
else
call xyzzyaact1(is)
call xyzzyaacu1(is)
call loglap_bf(ii,is,xyzzyaaak1(1,1,1,is),xyzzyaaal1(1,1,1,1,1,is),log&
&grad,lap)
endif
end subroutine xyzzyaacs1
subroutine xyzzyaact1(is)
implicit none
integer,intent(in) :: is
if(xyzzyaaax1(is))return
call xyzzyaadj1(xyzzyaaag1(:,:,:,:,is),xyzzyaaaf1(:,:,:,is),xyzzyaaak1&
&(:,:,:,is))
xyzzyaaax1(is)=.true.
end subroutine xyzzyaact1
subroutine xyzzyaacu1(is)
implicit none
integer,intent(in) :: is
if(xyzzyaaay1(is))return
call get_bf_x(is,.true.,.true.,.true.)
call xyzzyaadk1(xyzzyaaag1(:,:,:,:,is),xyzzyaaai1(:,:,:,:,is),xyzzyaaa&
&f1(:,:,:,is),bf_m2_scr(1,is),bf_rmap2_scr(1,1,is),xyzzyaaak1(:,:,:,is&
&),xyzzyaaal1(:,:,:,:,:,is))
xyzzyaaay1(is)=.true.
end subroutine xyzzyaacu1
subroutine xyzzyaacv1(is,args,fsd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fsd,sd
real(dp),intent(out) :: args(dimensionality,nemax,xyzzyaacn1)
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19(4)
real(dp) xyzzyaaad19(4,nemax)
args=0.d0
xyzzyaaac19(1)=1
xyzzyaaac19(2)=heg_nele(1)+1
xyzzyaaac19(3)=heg_nele(1)
xyzzyaaac19(4)=heg_nele(1)+heg_nele(2)
if(.not.use_backflow)then
call get_rsele(is)
do xyzzyaaab19=1,xyzzyaacn1
do xyzzyaaaa19=1,dimensionality
args(xyzzyaaaa19,1:heg_nele(xyzzyaaab19),xyzzyaaab19)=rele_scr(xyzzyaa&
&aa19,xyzzyaaac19(xyzzyaaab19):xyzzyaaac19(xyzzyaaab19+2),is)-rele_scr&
&(xyzzyaaaa19,netot,is)
enddo
enddo
else
if(.not.bf_x_valid(is))call get_bf_x(is,.true.,fsd,sd)
do xyzzyaaab19=1,xyzzyaacn1
do xyzzyaaaa19=1,dimensionality
args(xyzzyaaaa19,1:heg_nele(xyzzyaaab19),xyzzyaaab19)=bf_x_scr(xyzzyaa&
&aa19,xyzzyaaac19(xyzzyaaab19):xyzzyaaac19(xyzzyaaab19+2),is)-bf_x_scr&
&(xyzzyaaaa19,netot,is)
enddo
enddo
endif
do xyzzyaaab19=1,xyzzyaacn1
xyzzyaaad19=0.d0
xyzzyaaad19(1:dimensionality,1:heg_nele(xyzzyaaab19))=args(1:dimension&
&ality,1:heg_nele(xyzzyaaab19),xyzzyaaab19)
xyzzyaaad19(4,1:heg_nele(xyzzyaaab19))=sum(args(1:dimensionality,1:heg&
&_nele(xyzzyaaab19),xyzzyaaab19)**2,1)
call minimum_image(4,heg_nele(xyzzyaaab19),xyzzyaaad19)
args(1:dimensionality,1:heg_nele(xyzzyaaab19),xyzzyaaab19)=xyzzyaaad19&
&(1:dimensionality,1:heg_nele(xyzzyaaab19))
enddo
end subroutine xyzzyaacv1
subroutine xyzzyaacw1(is,args,fsd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fsd,sd
integer xyzzyaaaa20
real(dp) args(dimensionality,nemax,xyzzyaacn1)
if(xyzzyaaba1(is))then
return
else
do xyzzyaaaa20=1,xyzzyaacn1
if(heg_nele(xyzzyaaaa20)==0)cycle
call xyzzyaadq1(heg_nele(xyzzyaaaa20),args(1:dimensionality,1:heg_nele&
&(xyzzyaaaa20),xyzzyaaaa20),fsd,sd,xyzzyaaam1(1:3,1:xyzzyaabv1,1:heg_n&
&ele(xyzzyaaaa20),xyzzyaaaa20,is))
enddo
if(fsd.and.sd)xyzzyaaba1(is)=.true.
endif
end subroutine xyzzyaacw1
subroutine xyzzyaacx1(is,args,fsd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fsd,sd
integer xyzzyaaaa21
real(dp) args(dimensionality,nemax,xyzzyaacn1)
if(xyzzyaabb1(is).or..not.xyzzyaabq1)then
return
else
do xyzzyaaaa21=1,xyzzyaacn1
if(heg_nele(xyzzyaaaa21)==0)cycle
call xyzzyaadr1(heg_nele(xyzzyaaaa21),args(1:dimensionality,1:heg_nele&
&(xyzzyaaaa21),xyzzyaaaa21),fsd,sd,xyzzyaaan1(1:3,1:xyzzyaabv1,1:heg_n&
&ele(xyzzyaaaa21),xyzzyaaaa21,is))
enddo
if(fsd.and.sd)xyzzyaabb1(is)=.true.
endif
end subroutine xyzzyaacx1
subroutine setup_mahan_params(nparam)
implicit none
integer,intent(inout) :: nparam
integer xyzzyaaaa22,xyzzyaaab22,xyzzyaaac22,xyzzyaaad22,xyzzyaaae22,xy&
&zzyaaaf22,xyzzyaaag22
nparam=0
xyzzyaabf1=0
xyzzyaabg1=0
xyzzyaaaa22=0
if(opt_orbitals)then
do xyzzyaaag22=2,xyzzyaabv1
do xyzzyaaab22=1,xyzzyaabc1
if(xyzzyaabi1(xyzzyaaab22,xyzzyaaag22)==1)xyzzyaaaa22=xyzzyaaaa22+1
enddo
enddo
if(xyzzyaabq1)then
do xyzzyaaag22=2,xyzzyaabv1
do xyzzyaaab22=1,xyzzyaabd1
if(xyzzyaabj1(xyzzyaaab22,xyzzyaaag22)==1)xyzzyaaaa22=xyzzyaaaa22+1
enddo
enddo
endif
endif
xyzzyaabg1(1)=xyzzyaaaa22
call setup_bf_params(xyzzyaabg1(2))
xyzzyaabf1=sum(xyzzyaabg1,xyzzyaabg1>0)
nparam=xyzzyaabf1
allocate(xyzzyaabh1(xyzzyaabf1),stat=xyzzyaaac22)
call check_alloc(xyzzyaaac22,'SETUP_MAHAN_PARAMS','mahan_param_secs')
xyzzyaaae22=0
do xyzzyaaaf22=1,xyzzyaabe1
if(xyzzyaabg1(xyzzyaaaf22)<1)cycle
xyzzyaaad22=xyzzyaaae22+1
xyzzyaaae22=xyzzyaaae22+xyzzyaabg1(xyzzyaaaf22)
xyzzyaabh1(xyzzyaaad22:xyzzyaaae22)=xyzzyaaaf22
enddo
call xyzzyaacy1
end subroutine setup_mahan_params
subroutine finish_mahan_params
implicit none
call finish_bf_params
call xyzzyaacz1
deallocate(xyzzyaabh1)
end subroutine finish_mahan_params
subroutine xyzzyaacy1
implicit none
integer xyzzyaaaa24,xyzzyaaab24
xyzzyaaab24=xyzzyaabg1(1)
if(opt_orbitals)then
allocate(xyzzyaabm1(1:xyzzyaabc1,2:xyzzyaabv1,0:xyzzyaaab24),stat=xyzz&
&yaaaa24)
call check_alloc(xyzzyaaaa24,'SETUP_MAHAN_PBUFFER','eta_params')
xyzzyaabm1=0.d0
xyzzyaabo1=xyzzyaabc1*(xyzzyaabv1-1)
if(xyzzyaabq1)then
allocate(xyzzyaabn1(1:xyzzyaabd1,2:xyzzyaabv1,0:xyzzyaaab24),stat=xyzz&
&yaaaa24)
call check_alloc(xyzzyaaaa24,'SETUP_MAHAN_PBUFFER','u_params')
xyzzyaabn1=0.d0
xyzzyaabp1=xyzzyaabd1*(xyzzyaabv1-1)
endif
endif
end subroutine xyzzyaacy1
subroutine xyzzyaacz1
implicit none
if(opt_orbitals)then
deallocate(xyzzyaabm1)
if(xyzzyaabq1)deallocate(xyzzyaabn1)
endif
end subroutine xyzzyaacz1
subroutine get_mahan_params(params,has_lolim,lolim,has_hilim,hilim,is_&
&shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,lab&
&el)
implicit none
real(dp),intent(inout) :: params(xyzzyaabf1),lolim(xyzzyaabf1),hilim(x&
&yzzyaabf1)
logical,intent(inout) :: has_lolim(xyzzyaabf1),has_hilim(xyzzyaabf1),i&
&s_shallow(xyzzyaabf1),is_redundant(xyzzyaabf1),is_linear(xyzzyaabf1),&
&is_loglinear(xyzzyaabf1), has_aderiv(xyzzyaabf1),affect_map(xyzzyaabf&
&1,xyzzyaabf1)
character(2),intent(inout) :: label(xyzzyaabf1)
integer xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26,xyzzyaaad26,xyzzyaaae26,xy&
&zzyaaaf26
real(dp) xyzzyaaag26,xyzzyaaah26
affect_map=.false.
label='Ma'
xyzzyaaag26=1.1d-8
xyzzyaaah26=0.999999d0*wigner_seitz_radius
xyzzyaaad26=0
if(xyzzyaabg1(1)>0)then
xyzzyaaaa26=xyzzyaaad26
xyzzyaaac26=xyzzyaaad26+1
xyzzyaaad26=xyzzyaaad26+xyzzyaabg1(1)
has_lolim(xyzzyaaac26:xyzzyaaad26)=.false.
lolim(xyzzyaaac26:xyzzyaaad26)=0.d0
has_hilim(xyzzyaaac26:xyzzyaaad26)=.false.
hilim(xyzzyaaac26:xyzzyaaad26)=0.d0
is_shallow(xyzzyaaac26:xyzzyaaad26)=.false.
is_redundant(xyzzyaaac26:xyzzyaaad26)=.false.
is_linear(xyzzyaaac26:xyzzyaaad26)=.false.
is_loglinear(xyzzyaaac26:xyzzyaaad26)=.false.
has_aderiv(xyzzyaaac26:xyzzyaaad26)=.false.
do xyzzyaaae26=1,xyzzyaabf1
affect_map(xyzzyaaae26,xyzzyaaae26)=.true.
enddo
if(opt_orbitals)then
do xyzzyaaaf26=2,xyzzyaabv1
do xyzzyaaab26=1,xyzzyaabc1
if(xyzzyaabi1(xyzzyaaab26,xyzzyaaaf26)==1)then
xyzzyaaaa26=xyzzyaaaa26+1
params(xyzzyaaaa26)=xyzzyaabk1(xyzzyaaab26,xyzzyaaaf26)
if(.not.(xyzzyaaab26==1.or.xyzzyaaab26==2))params(xyzzyaaaa26)=params(&
&xyzzyaaaa26)*xyzzyaabk1(2,xyzzyaaaf26)**(xyzzyaaab26-3)
select case(1)
case(1)
if(xyzzyaaab26==2)then
has_lolim(xyzzyaaaa26)=.true.
lolim(xyzzyaaaa26)=xyzzyaaag26
has_hilim(xyzzyaaaa26)=.true.
hilim(xyzzyaaaa26)=xyzzyaaah26
is_shallow(xyzzyaaaa26)=.true.
endif
end select
endif
enddo
enddo
if(xyzzyaabq1)then
do xyzzyaaaf26=2,xyzzyaabv1
do xyzzyaaab26=1,xyzzyaabd1
if(xyzzyaabj1(xyzzyaaab26,xyzzyaaaf26)==1)then
xyzzyaaaa26=xyzzyaaaa26+1
params(xyzzyaaaa26)=xyzzyaabl1(xyzzyaaab26,xyzzyaaaf26)
if(.not.(xyzzyaaab26==1.or.xyzzyaaab26==2))params(xyzzyaaaa26)=params(&
&xyzzyaaaa26)*xyzzyaabl1(2,xyzzyaaaf26)**(xyzzyaaab26-3)
select case(1)
case(1)
if(xyzzyaaab26==2)then
has_lolim(xyzzyaaaa26)=.true.
lolim(xyzzyaaaa26)=xyzzyaaag26
has_hilim(xyzzyaaaa26)=.true.
hilim(xyzzyaaaa26)=xyzzyaaah26
is_shallow(xyzzyaaaa26)=.true.
endif
end select
endif
enddo
enddo
endif
endif
endif
if(xyzzyaabg1(2)>0)then
xyzzyaaac26=xyzzyaaad26+1
xyzzyaaad26=xyzzyaaad26+xyzzyaabg1(2)
call get_bf_params(params(xyzzyaaac26:xyzzyaaad26),has_lolim(xyzzyaaac&
&26:xyzzyaaad26),lolim(xyzzyaaac26:xyzzyaaad26),has_hilim(xyzzyaaac26:&
&xyzzyaaad26),hilim(xyzzyaaac26:xyzzyaaad26),is_shallow(xyzzyaaac26:xy&
&zzyaaad26),is_redundant(xyzzyaaac26:xyzzyaaad26),is_linear(xyzzyaaac2&
&6:xyzzyaaad26),is_loglinear(xyzzyaaac26:xyzzyaaad26),has_aderiv(xyzzy&
&aaac26:xyzzyaaad26),affect_map(xyzzyaaac26:xyzzyaaad26,xyzzyaaac26:xy&
&zzyaaad26),label(xyzzyaaac26:xyzzyaaad26))
endif
end subroutine get_mahan_params
subroutine put_mahan_params(params,ignore,iparam_buffer,prestore,bad_p&
&arams)
implicit none
integer,intent(in) :: iparam_buffer
logical,intent(in) :: ignore(xyzzyaabf1),prestore
logical,intent(out) :: bad_params
real(dp),intent(inout) :: params(xyzzyaabf1)
integer xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xyzzyaaad27,xyzzyaaae27,xy&
&zzyaaaf27,xyzzyaaag27
logical xyzzyaaah27
bad_params=.false.
xyzzyaaac27=0
xyzzyaaad27=0
if(iparam_buffer>0)then
xyzzyaaac27=xyzzyaabh1(iparam_buffer)
xyzzyaaad27=iparam_buffer-sum(xyzzyaabg1(1:xyzzyaaac27-1))
endif
xyzzyaaaf27=0
if(xyzzyaabg1(1)>0)then
xyzzyaaaa27=xyzzyaaaf27
xyzzyaaae27=xyzzyaaaf27+1
xyzzyaaaf27=xyzzyaaaf27+xyzzyaabg1(1)
if(xyzzyaaac27==1.or.xyzzyaaac27==0)then
if(prestore)then
call xyzzyaadb1(iparam_buffer)
else
if(opt_orbitals)then
do xyzzyaaag27=2,xyzzyaabv1
do xyzzyaaab27=1,xyzzyaabc1
if(xyzzyaabi1(xyzzyaaab27,xyzzyaaag27)==1)then
xyzzyaaaa27=xyzzyaaaa27+1
if(.not.ignore(xyzzyaaaa27))then
if(.not.(xyzzyaaab27==1.or.xyzzyaaab27==2))params(xyzzyaaaa27)=params(&
&xyzzyaaaa27)/xyzzyaabk1(2,xyzzyaaag27)**(xyzzyaaab27-3)
xyzzyaabk1(xyzzyaaab27,xyzzyaaag27)=params(xyzzyaaaa27)
endif
endif
enddo
enddo
select case(1)
case(1)
end select
if(xyzzyaabq1)then
do xyzzyaaag27=2,xyzzyaabv1
do xyzzyaaab27=1,xyzzyaabd1
if(xyzzyaabj1(xyzzyaaab27,xyzzyaaag27)==1)then
xyzzyaaaa27=xyzzyaaaa27+1
if(.not.ignore(xyzzyaaaa27))then
if(.not.(xyzzyaaab27==1.or.xyzzyaaab27==2))params(xyzzyaaaa27)=params(&
&xyzzyaaaa27)/xyzzyaabl1(2,xyzzyaaag27)**(xyzzyaaab27-3)
xyzzyaabl1(xyzzyaaab27,xyzzyaaag27)=params(xyzzyaaaa27)
endif
endif
enddo
select case(1)
case(1)
xyzzyaabl1(3,xyzzyaaag27)=xyzzyaabl1(4,xyzzyaaag27)*xyzzyaabl1(2,xyzzy&
&aaag27)/xyzzyaabl1(1,xyzzyaaag27)
end select
enddo
endif
endif
call xyzzyaada1(iparam_buffer)
endif
endif
endif
if(xyzzyaabg1(2)>0)then
xyzzyaaae27=xyzzyaaaf27+1
xyzzyaaaf27=xyzzyaaaf27+xyzzyaabg1(2)
if(xyzzyaaac27==2.or.xyzzyaaac27==0)then
call put_bf_params(params(xyzzyaaae27:xyzzyaaaf27),ignore(xyzzyaaae27:&
&xyzzyaaaf27),xyzzyaaad27,prestore,xyzzyaaah27)
bad_params=bad_params.or.xyzzyaaah27
endif
endif
end subroutine put_mahan_params
subroutine xyzzyaada1(indx)
implicit none
integer,intent(in) :: indx
if(opt_orbitals)then
call dcopy(xyzzyaabo1,xyzzyaabk1(1,2),1,xyzzyaabm1(1,2,indx),1)
if(xyzzyaabq1)then
call dcopy(xyzzyaabp1,xyzzyaabl1(1,2),1,xyzzyaabn1(1,2,indx),1)
endif
endif
end subroutine xyzzyaada1
subroutine xyzzyaadb1(indx)
implicit none
integer,intent(in) :: indx
if(opt_orbitals)then
call dcopy(xyzzyaabo1,xyzzyaabm1(1,2,indx),1,xyzzyaabk1(1,2),1)
if(xyzzyaabq1)then
call dcopy(xyzzyaabp1,xyzzyaabn1(1,2,indx),1,xyzzyaabl1(1,2),1)
endif
endif
end subroutine xyzzyaadb1
subroutine invalidate_param1_mahan(is,iparam)
implicit none
integer,intent(in) :: is,iparam
integer xyzzyaaaa30,xyzzyaaab30
xyzzyaaaa30=xyzzyaabh1(iparam)
xyzzyaaab30=iparam-sum(xyzzyaabg1(1:xyzzyaaaa30-1))
xyzzyaaaq1(is)=.false.
xyzzyaaar1(is)=.false.
xyzzyaaas1(is)=.false.
xyzzyaaat1(is)=.false.
xyzzyaaau1(is)=.false.
xyzzyaaaz1(is)=.false.
xyzzyaaav1(:,is)=.false.
xyzzyaaaw1(:,is)=.false.
xyzzyaaax1(is)=.false.
xyzzyaaay1(is)=.false.
xyzzyaaba1(is)=.false.
xyzzyaabb1(is)=.false.
if(xyzzyaaaa30==2)call invalidate_param1_bf(is,xyzzyaaab30)
end subroutine invalidate_param1_mahan
subroutine invalidate_params_mahan(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaabf1)
end subroutine invalidate_params_mahan
subroutine clear_scratch_mahan(is)
implicit none
integer,intent(in) :: is
xyzzyaaaq1(is)=.false.
xyzzyaaar1(is)=.false.
xyzzyaaas1(is)=.false.
xyzzyaaat1(is)=.false.
xyzzyaaau1(is)=.false.
xyzzyaaaz1(is)=.false.
xyzzyaaav1(:,is)=.false.
xyzzyaaaw1(:,is)=.false.
xyzzyaaax1(is)=.false.
xyzzyaaay1(is)=.false.
xyzzyaaba1(is)=.false.
xyzzyaabb1(is)=.false.
if(use_backflow)call clear_scratch_bf(is)
end subroutine clear_scratch_mahan
subroutine gen_config_mahan(pt_config)
implicit none
integer xyzzyaaaa33
type(config_wfn_mahan),pointer :: pt_config
allocate(pt_config,stat=xyzzyaaaa33)
call check_alloc(xyzzyaaaa33,'GEN_CONFIG_GEMINAL','container')
end subroutine gen_config_mahan
subroutine delete_config_mahan(pt_config)
implicit none
type(config_wfn_mahan),pointer :: pt_config
deallocate(pt_config)
end subroutine delete_config_mahan
subroutine copy_config_mahan(pt_from,pt_to)
implicit none
type(config_wfn_mahan),pointer :: pt_from,pt_to
end subroutine copy_config_mahan
subroutine config_to_pt_mahan(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_mahan),pointer :: pt_config
end subroutine config_to_pt_mahan
subroutine pt_to_config_mahan(pt_config)
implicit none
type(config_wfn_mahan),pointer :: pt_config
end subroutine pt_to_config_mahan
subroutine redist_allocations_mahan(kmax)
implicit none
integer,intent(in) :: kmax
end subroutine redist_allocations_mahan
subroutine redist_load_mahan(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_mahan),pointer :: pt_config
end subroutine redist_load_mahan
subroutine redist_send_mahan(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
end subroutine redist_send_mahan
subroutine redist_recv_mahan(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
end subroutine redist_recv_mahan
subroutine redist_save_mahan(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_mahan),pointer :: pt_config
end subroutine redist_save_mahan
subroutine redist_deallocations_mahan
implicit none
end subroutine redist_deallocations_mahan
subroutine load_from_pt_mahan(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_mahan),pointer :: pt_config
end subroutine load_from_pt_mahan
subroutine save_to_pt_mahan(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_mahan),pointer :: pt_config
end subroutine save_to_pt_mahan
subroutine clone_scratch_mahan(is,js)
implicit none
integer,intent(in) :: is,js
end subroutine clone_scratch_mahan
subroutine add_config_mahan_items(is)
implicit none
integer,intent(in) :: is
end subroutine add_config_mahan_items
subroutine setup_storage_mahan(nconfig,ignore)
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(xyzzyaabf1)
end subroutine setup_storage_mahan
subroutine finish_storage_mahan
implicit none
end subroutine finish_storage_mahan
subroutine load_from_storage_mahan(is,icfg)
implicit none
integer,intent(in) :: is,icfg
end subroutine load_from_storage_mahan
subroutine save_to_storage_mahan(is,icfg)
implicit none
integer,intent(in) :: is,icfg
end subroutine save_to_storage_mahan
subroutine enumerate_plot_mahan(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
n=0
end subroutine enumerate_plot_mahan
subroutine query_plot_mahan(iplot,ii,rank,is_complex,has_stderr,rot_te&
&nsor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
nfunctions=0
end subroutine query_plot_mahan
subroutine get_plot_mahan(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
f(1)=0.d0
end subroutine get_plot_mahan
subroutine finish_plot_mahan
implicit none
end subroutine finish_plot_mahan
subroutine read_mahan
implicit none
integer xyzzyaaaa56,xyzzyaaab56,xyzzyaaac56,xyzzyaaad56,xyzzyaaae56,xy&
&zzyaaaf56,xyzzyaaag56,xyzzyaaah56
character(80) char_80
real(dp) xyzzyaaai56,xyzzyaaaj56
xyzzyaaai56=0.999999d0*wigner_seitz_radius
allocate(xyzzyaabt1(dimensionality,dimensionality))
if(dimensionality==3)then
xyzzyaabt1(1:3,1)=(/1,4,5/)
xyzzyaabt1(1:3,2)=(/4,2,6/)
xyzzyaabt1(1:3,3)=(/5,6,3/)
elseif(dimensionality==2)then
xyzzyaabt1(1:2,1)=(/1,3/)
xyzzyaabt1(1:2,2)=(/3,2/)
endif
xyzzyaacn1=2
allocate(xyzzyaacb1(xyzzyaacn1))
call xyzzyaadd1()
xyzzyaacb1(:)=-1
do xyzzyaaab56=1,xyzzyaacn1
if(heg_nele(xyzzyaaab56)==1)then
xyzzyaacb1(xyzzyaaab56)=1
elseif(heg_nele(xyzzyaaab56)==0)then
xyzzyaacb1(xyzzyaaab56)=0
else
do xyzzyaaaa56=2,xyzzyaabx1
if(heg_nele(xyzzyaaab56)==xyzzyaabz1(xyzzyaaaa56+1)-1)then
xyzzyaacb1(xyzzyaaab56)=xyzzyaaaa56
exit
endif
enddo
endif
enddo
xyzzyaabv1=maxval(xyzzyaacb1(:))
if(any(xyzzyaacb1==-1))then
call errstop('READ_MAHAN','Not using a magic number of electrons.')
endif
if(xyzzyaaco1)then
call wout('Setup Mahan debugging info:')
call wout('nspin is: ',nspin)
call wout('num_spin is: ',xyzzyaacn1)
call wout('heg_nele(:) is: ',heg_nele(:))
call wout('spin_nstar(:) is: ',xyzzyaacb1(:))
call wout('nstar is: ',xyzzyaabv1)
endif
call open_units(xyzzyaabs1,xyzzyaaac56)
if(xyzzyaaac56/=0)call errstop('READ_MAHAN','Unable to find free i/o u&
&nit.')
open(unit=xyzzyaabs1,file='correlation.data',status='old',iostat=xyzzy&
&aaac56)
if(xyzzyaaac56/=0)call errstop('READ_MAHAN','Problem opening correlati&
&on.data')
if(am_master)then
call wout('Mahan exciton wave function')
call wout('===========================')
call wout('Reading correlation.data file.')
call wout()
endif
do
read(xyzzyaabs1,'(a)',iostat=xyzzyaaac56)char_80
if(trim(adjustl(char_80))=='START MAHAN')exit
if(xyzzyaaac56>0)call errstop_master('READ_MAHAN','Problem reading cor&
&relation.data.')
if(xyzzyaaac56<0)call errstop_master('READ_MAHAN','Could not find "STA&
&RT MAHAN" in correlation.data.')
enddo
read(xyzzyaabs1,*,err=100,end=101)
read(xyzzyaabs1,'(a)',err=101,end=101)title
do
read(xyzzyaabs1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='START ETA_G TERM')exit
if(trim(adjustl(char_80))=='END MAHAN')call errstop_master('READ_MAHAN&
&','Was expecting to find "START ETA_G TERM"')
enddo
read(xyzzyaabs1,*,err=100,end=101)
read(xyzzyaabs1,*,err=101,end=101)xyzzyaabc1
xyzzyaabc1=xyzzyaabc1+3
allocate(xyzzyaabk1(1:xyzzyaabc1,2:xyzzyaabv1))
allocate(xyzzyaabi1(xyzzyaabc1,2:xyzzyaabv1))
read(xyzzyaabs1,*,err=100,end=101)
read(xyzzyaabs1,*,err=101,end=101)xyzzyaaaf56
if(xyzzyaaaf56/=0)call errstop_master('READ_MAHAN','Spin-dep must be z&
&ero, currently.')
read(xyzzyaabs1,*,err=100,end=101)
read(xyzzyaabs1,*,err=101,end=101)xyzzyaaag56
if(xyzzyaaag56/=xyzzyaabv1-1)then
call wordwrap('The number of optimizable stars must be one less than t&
&he total number of occupied stars in the system.  The total number of&
& filled stars in the current system is: '//trim(adjustl(i2s(xyzzyaabv&
&1)))//'.')
call errstop_master('READ_MAHAN','Wrong number of optimizable stars.')
endif
read(xyzzyaabs1,*,err=100,end=101)
do xyzzyaaad56=2,xyzzyaabv1
read(xyzzyaabs1,*,err=101,end=101)xyzzyaabk1(1,xyzzyaaad56)
xyzzyaabi1(1,xyzzyaaad56)=0
enddo
read(xyzzyaabs1,*,err=100,end=101)
do xyzzyaaad56=2,xyzzyaabv1
read(xyzzyaabs1,*,iostat=xyzzyaaac56)xyzzyaaaj56,xyzzyaaah56
if(xyzzyaaac56/=0)then
backspace xyzzyaabs1
xyzzyaaaj56=xyzzyaaai56
xyzzyaaah56=1
endif
if(xyzzyaaah56/=0.and.xyzzyaaah56/=1)call errstop_master('READ_MAHAN',&
&'Optimizable flag should be 0 or 1.')
if(xyzzyaaaj56==0.d0)xyzzyaaaj56=xyzzyaaai56
xyzzyaabk1(2,xyzzyaaad56)=xyzzyaaaj56
xyzzyaabi1(2,xyzzyaaad56)=xyzzyaaah56
enddo
read(xyzzyaabs1,*,err=100,end=101)
do xyzzyaaad56=2,xyzzyaabv1
do xyzzyaaae56=3,xyzzyaabc1
read(xyzzyaabs1,*,iostat=xyzzyaaac56)xyzzyaabk1(xyzzyaaae56,xyzzyaaad5&
&6),xyzzyaabi1(xyzzyaaae56,xyzzyaaad56)
if(xyzzyaaac56/=0)then
backspace xyzzyaabs1
xyzzyaabk1(xyzzyaaae56,xyzzyaaad56)=0.d0
xyzzyaabi1(xyzzyaaae56,xyzzyaaad56)=1
endif
enddo
enddo
do
read(xyzzyaabs1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='END ETA_G TERM')exit
if(trim(adjustl(char_80))/='')call errstop_master('READ_MAHAN','Was ex&
&pecting to find "END ETA_G TERM"')
enddo
do
read(xyzzyaabs1,'(a)',iostat=xyzzyaaac56)char_80
if(trim(adjustl(char_80))=='END MAHAN')exit
if(xyzzyaaac56>0)call errstop_master('READ_MAHAN','Problem reading cor&
&relation.data.')
if(xyzzyaaac56<0)call errstop_master('READ_MAHAN','Could not find "END&
& MAHAN" in correlation.data.')
enddo
rewind(xyzzyaabs1)
do
read(xyzzyaabs1,'(a)',iostat=xyzzyaaac56)char_80
if(trim(adjustl(char_80))=='START MAHAN')exit
if(xyzzyaaac56>0)call errstop_master('READ_MAHAN','Problem reading cor&
&relation.data.')
if(xyzzyaaac56<0)call errstop_master('READ_MAHAN','Could not find "STA&
&RT MAHAN" in correlation.data.')
enddo
do
read(xyzzyaabs1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='START U_G TERM')then
xyzzyaabq1=.true.
xyzzyaabr1=3
exit
endif
if(trim(adjustl(char_80))=='END MAHAN')then
xyzzyaabq1=.false.
xyzzyaabr1=2
close(xyzzyaabs1)
if(am_master)call wout('No U_G term found in correlation.data file.  W&
&ill run without U_G term.')
if(am_master)call wout()
return
endif
enddo
read(xyzzyaabs1,*,err=100,end=101)
read(xyzzyaabs1,*,err=101,end=101)xyzzyaabd1
xyzzyaabd1=xyzzyaabd1+3
allocate(xyzzyaabl1(1:xyzzyaabd1,2:xyzzyaabv1))
allocate(xyzzyaabj1(xyzzyaabd1,2:xyzzyaabv1))
read(xyzzyaabs1,*,err=100,end=101)
read(xyzzyaabs1,*,err=101,end=101)xyzzyaaaf56
if(xyzzyaaaf56/=0)call errstop_master('READ_MAHAN','Spin-dep must be z&
&ero, currently.')
read(xyzzyaabs1,*,err=100,end=101)
read(xyzzyaabs1,*,err=101,end=101)xyzzyaaag56
if(xyzzyaaag56/=xyzzyaabv1-1)then
call wordwrap('The number of optimizable stars must be one less than t&
&he total number of occupied stars in the system.  The total number of&
& filled stars in the current system is: '//trim(adjustl(i2s(xyzzyaabv&
&1)))//'.')
call errstop_master('READ_MAHAN','Wrong number of optimizable stars.')
endif
read(xyzzyaabs1,*,err=100,end=101)
do xyzzyaaad56=2,xyzzyaabv1
read(xyzzyaabs1,*,err=101,end=101)xyzzyaabl1(1,xyzzyaaad56)
enddo
xyzzyaabj1(1,:)=0
read(xyzzyaabs1,*,err=100,end=101)
do xyzzyaaad56=2,xyzzyaabv1
read(xyzzyaabs1,*,iostat=xyzzyaaac56)xyzzyaaaj56,xyzzyaaah56
if(xyzzyaaac56/=0)then
backspace xyzzyaabs1
xyzzyaaaj56=xyzzyaaai56
xyzzyaaah56=1
endif
if(xyzzyaaah56/=0.and.xyzzyaaah56/=1)call errstop_master('READ_MAHAN',&
&'Optimizable flag should be 0 or 1.')
if(xyzzyaaaj56==0.d0)xyzzyaaaj56=xyzzyaaai56
xyzzyaabl1(2,xyzzyaaad56)=xyzzyaaaj56
xyzzyaabj1(2,xyzzyaaad56)=xyzzyaaah56
enddo
read(xyzzyaabs1,*,err=100,end=101)
do xyzzyaaad56=2,xyzzyaabv1
do xyzzyaaae56=4,xyzzyaabd1
read(xyzzyaabs1,*,iostat=xyzzyaaac56)xyzzyaabl1(xyzzyaaae56,xyzzyaaad5&
&6),xyzzyaabj1(xyzzyaaae56,xyzzyaaad56)
if(xyzzyaaac56/=0)then
backspace xyzzyaabs1
xyzzyaabl1(xyzzyaaae56,xyzzyaaad56)=0.d0
xyzzyaabj1(xyzzyaaae56,xyzzyaaad56)=1
endif
enddo
xyzzyaabl1(3,xyzzyaaad56)=xyzzyaabl1(4,xyzzyaaad56)*xyzzyaabl1(2,xyzzy&
&aaad56)/xyzzyaabl1(1,xyzzyaaad56)
enddo
xyzzyaabj1(3,:)=0
do
read(xyzzyaabs1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='END U_G TERM')exit
if(trim(adjustl(char_80))/='')call errstop_master('READ_MAHAN','Was ex&
&pecting to find "END U_G TERM"')
enddo
do
read(xyzzyaabs1,'(a)',iostat=xyzzyaaac56)char_80
if(trim(adjustl(char_80))=='END MAHAN')exit
if(xyzzyaaac56>0)call errstop_master('READ_MAHAN','Problem reading cor&
&relation.data.')
if(xyzzyaaac56<0)call errstop_master('READ_MAHAN','Could not find "END&
& MAHAN" in correlation.data.')
enddo
close(xyzzyaabs1)
if(am_master)then
call wout('Title: '//trim(adjustl(title)))
call wout()
call wout('Eta_G term: ')
call wout(' Expansion order: '//trim(i2s(xyzzyaabc1-3)))
call wout(' Spin dependence: '//trim(i2s(xyzzyaaaf56)))
call wout(' Total number of stars (for majority spin): '//trim(i2s(xyz&
&zyaabv1)))
call wout(' Number of optimisable stars: '//trim(i2s(xyzzyaaag56)))
do xyzzyaaad56=2,xyzzyaabv1
call wout('  Parameters for star number: '//trim(i2s(xyzzyaaad56)))
char_80=r2s2(xyzzyaabk1(1,xyzzyaaad56),'(f21.12)')
call wout('    Truncation order'//' (always fixed) :  '//trim(char_80)&
&)
char_80=r2s2(xyzzyaabk1(2,xyzzyaaad56),'(f21.12)')
if(xyzzyaabi1(2,xyzzyaaad56)==1)then
call wout('    Cutoff length'//'             (opt) :  '//trim(char_80)&
&)
else
call wout('    Cutoff length'//'           (fixed) :  '//trim(char_80)&
&)
endif
do xyzzyaaaa56=3,xyzzyaabc1
char_80=r2s2(xyzzyaabk1(xyzzyaaaa56,xyzzyaaad56),'(f21.12)')
if(xyzzyaabi1(xyzzyaaaa56,xyzzyaaad56)==1)then
call wout('    c_'//trim(i2s(xyzzyaaad56))//'_'//trim(i2s(xyzzyaaaa56-&
&3))//'                     (opt) :  '//trim(char_80))
else
call wout('    c_'//trim(i2s(xyzzyaaad56))//'_'//trim(i2s(xyzzyaaaa56-&
&3))//'                   (fixed) :  '//trim(char_80))
endif
enddo
enddo
if(xyzzyaabq1)then
call wout()
call wout('Found u_G term in correlation.data;')
call wout('will run with orbital-dependent Jastrow.')
call wout()
call wout('u_G term: ')
call wout(' Expansion order: '//trim(i2s(xyzzyaabd1-3)))
call wout(' Spin dependence: '//trim(i2s(xyzzyaaaf56)))
call wout(' Total number of stars (for majority spin): '//trim(i2s(xyz&
&zyaabv1)))
call wout(' Number of optimisable stars: '//trim(i2s(xyzzyaaag56)))
do xyzzyaaad56=2,xyzzyaabv1
call wout('  Parameters for star number: '//trim(i2s(xyzzyaaad56)))
char_80=r2s2(xyzzyaabl1(1,xyzzyaaad56),'(f21.12)')
call wout('    Truncation order'//' (always fixed) :  '//trim(char_80)&
&)
char_80=r2s2(xyzzyaabl1(2,xyzzyaaad56),'(f21.12)')
if(xyzzyaabj1(2,xyzzyaaad56)==1)then
call wout('    Cutoff length'//'             (opt) :  '//trim(char_80)&
&)
else
call wout('    Cutoff length'//'           (fixed) :  '//trim(char_80)&
&)
endif
do xyzzyaaaa56=4,xyzzyaabd1
char_80=r2s2(xyzzyaabl1(xyzzyaaaa56,xyzzyaaad56),'(f21.12)')
if(xyzzyaabj1(xyzzyaaaa56,xyzzyaaad56)==1)then
call wout('    c_'//trim(i2s(xyzzyaaad56))//'_'//trim(i2s(xyzzyaaaa56-&
&3))//'                     (opt) :  '//trim(char_80))
else
call wout('    c_'//trim(i2s(xyzzyaaad56))//'_'//trim(i2s(xyzzyaaaa56-&
&3))//'                   (fixed) :  '//trim(char_80))
endif
enddo
enddo
else
call wout('No u_G term found in correlation.data;')
call wout('will run without orbital-dependent Jastrow.')
endif
call wout()
call wout('Finished Mahan setup.')
call wout()
endif
return
100 call errstop_master('READ_MAHAN','Error reading correlation.data .&
&')
101 call errstop_master('READ_MAHAN','File correlation.data ended unex&
&pectedly.')
end subroutine read_mahan
subroutine write_mahan(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
integer xyzzyaaaa57,xyzzyaaab57,xyzzyaaac57
logical xyzzyaaad57
if(.not.am_master)return
inquire(file=trim(correlation_name),exist=xyzzyaaad57)
if(xyzzyaaad57)then
open(unit=xyzzyaabs1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaaa57)
else
open(unit=xyzzyaabs1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaaa57)
endif
if(xyzzyaaaa57/=0)call errstop('WRITE_MAHAN','Problem opening '//trim(&
&correlation_name)//'.')
write(xyzzyaabs1,*)'START MAHAN'
write(xyzzyaabs1,*)'Title'
write(xyzzyaabs1,*)trim(adjustl(title))
write(xyzzyaabs1,*)'START ETA_G TERM'
write(xyzzyaabs1,*)'Expansion order'
write(xyzzyaabs1,'(3x,a)')trim(i2s(xyzzyaabc1-3))
write(xyzzyaabs1,*)'Spin-dep (0->u=d ; 1->u/=d)'
write(xyzzyaabs1,'(3x,a)')trim(i2s(0))
write(xyzzyaabs1,*)'Number of optimizable stars'
write(xyzzyaabs1,'(3x,a)')trim(i2s(xyzzyaabv1-1))
write(xyzzyaabs1,*)'Truncation orders'
do xyzzyaaac57=2,xyzzyaabv1
write(xyzzyaabs1,*)xyzzyaabk1(1,xyzzyaaac57),'      ! C_',trim(i2s(xyz&
&zyaaac57))
enddo
write(xyzzyaabs1,*)'Cutoff radii; Optimizable (0=NO; 1=Yes)'
do xyzzyaaac57=2,xyzzyaabv1
write(xyzzyaabs1,*)xyzzyaabk1(2,xyzzyaaac57),xyzzyaabi1(2,xyzzyaaac57)&
&,'      ! L_',trim(i2s(xyzzyaaac57))
enddo
write(xyzzyaabs1,*)'Parameters; Optimizable (0=NO; 1=Yes)'
do xyzzyaaac57=2,xyzzyaabv1
do xyzzyaaab57=3,xyzzyaabc1
write(xyzzyaabs1,*)xyzzyaabk1(xyzzyaaab57,xyzzyaaac57),xyzzyaabi1(xyzz&
&yaaab57,xyzzyaaac57),'      ! c_',trim(i2s(xyzzyaaac57)),'_',trim(i2s&
&(xyzzyaaab57-3))
enddo
enddo
write(xyzzyaabs1,*)'END ETA_G TERM'
if(xyzzyaabq1)then
write(xyzzyaabs1,*)'START U_G TERM'
write(xyzzyaabs1,*)'Expansion order'
write(xyzzyaabs1,'(3x,a)')trim(i2s(xyzzyaabd1-3))
write(xyzzyaabs1,*)'Spin-dep (0->u=d ; 1->u/=d)'
write(xyzzyaabs1,'(3x,a)')trim(i2s(0))
write(xyzzyaabs1,*)'Number of optimizable stars'
write(xyzzyaabs1,'(3x,a)')trim(i2s(xyzzyaabv1-1))
write(xyzzyaabs1,*)'Truncation orders'
do xyzzyaaac57=2,xyzzyaabv1
write(xyzzyaabs1,*)xyzzyaabl1(1,xyzzyaaac57),'      ! C_',trim(i2s(xyz&
&zyaaac57))
enddo
write(xyzzyaabs1,*)'Cutoff radii; Optimizable (0=NO; 1=Yes)'
do xyzzyaaac57=2,xyzzyaabv1
write(xyzzyaabs1,*)xyzzyaabl1(2,xyzzyaaac57),xyzzyaabj1(2,xyzzyaaac57)&
&,'      ! L_',trim(i2s(xyzzyaaac57))
enddo
write(xyzzyaabs1,*)'Parameters; Optimizable (0=NO; 1=Yes)'
do xyzzyaaac57=2,xyzzyaabv1
do xyzzyaaab57=4,xyzzyaabd1
write(xyzzyaabs1,*)xyzzyaabl1(xyzzyaaab57,xyzzyaaac57),xyzzyaabj1(xyzz&
&yaaab57,xyzzyaaac57),'      ! c_',trim(i2s(xyzzyaaac57)),'_',trim(i2s&
&(xyzzyaaab57-3))
enddo
enddo
write(xyzzyaabs1,*)'END U_G TERM'
endif
write(xyzzyaabs1,*)'END MAHAN'
write(xyzzyaabs1,*)
close(xyzzyaabs1)
end subroutine write_mahan
subroutine xyzzyaadc1(n,fsd,logdet,mmat,minv,fpeinfo)
implicit none
integer,intent(in) :: n
logical,intent(in) :: fsd
integer,intent(inout) :: fpeinfo
complex(dp),intent(in) :: mmat(nemax,nemax)
complex(dp),intent(out) :: minv(nemax,nemax)
complex(dp),intent(inout) :: logdet
integer xyzzyaaaa58
logical xyzzyaaab58
complex(dp) xyzzyaaac58(n,n)
integer xyzzyaaad58(n)
complex(dp) xyzzyaaae58(n,n)
xyzzyaaac58(:,:)=czero
logdet=czero
fpeinfo=xyzzyaack1
xyzzyaaab58=.false.
xyzzyaaad58=0
xyzzyaaae58=mmat(1:n,1:n)
call lu_decom_cmplx(xyzzyaaae58,xyzzyaaad58,n,xyzzyaaab58)
if(.not.xyzzyaaab58)then
logdet=lu_logdet_cmplx(xyzzyaaae58,xyzzyaaad58,n,n)
if(dble(logdet)>xyzzyaaci1)then
if(fsd)then
xyzzyaaac58=0.d0
do xyzzyaaaa58=1,n
xyzzyaaac58(xyzzyaaaa58,xyzzyaaaa58)=1.d0
enddo
call lu_solve_n_cmplx(xyzzyaaae58,xyzzyaaad58,xyzzyaaac58,n,n,n)
endif
else
xyzzyaaab58=.true.
endif
endif
if(xyzzyaaab58)then
fpeinfo=xyzzyaacm1
endif
xyzzyaaac58=transpose(xyzzyaaac58)
minv=cmplx(0.d0,0.d0,dp)
minv(1:n,1:n)=xyzzyaaac58
end subroutine xyzzyaadc1
subroutine xyzzyaadd1()
use slaarnaan,only : lattice_generator
use slaarnabg,     only : periodicity,pb1,pb2,pb3
implicit none
integer xyzzyaaaa59,xyzzyaaab59,xyzzyaaac59,xyzzyaaad59,xyzzyaaae59,xy&
&zzyaaaf59
character(80) char80
xyzzyaabw1=nemax*100
if(mod(xyzzyaabw1,2)/=1)xyzzyaabw1=xyzzyaabw1+1
allocate(xyzzyaaby1(3,xyzzyaabw1),xyzzyaabz1((xyzzyaabw1+3)/2))
allocate(xyzzyaacd1((xyzzyaabw1+3)/2),xyzzyaace1(3,xyzzyaabw1),xyzzyaa&
&cf1(3,3))
call lattice_generator(xyzzyaabw1,periodicity,pb1,pb2,pb3,xyzzyaace1,x&
&yzzyaacd1,xyzzyaabz1,xyzzyaacf1,xyzzyaaby1,xyzzyaabx1)
allocate(xyzzyaaca1(xyzzyaabx1))
do xyzzyaaab59=1,xyzzyaabx1-1
xyzzyaaca1(xyzzyaaab59)=xyzzyaabz1(xyzzyaaab59+1)-xyzzyaabz1(xyzzyaaab&
&59)
enddo
xyzzyaaca1(xyzzyaabx1)=xyzzyaabw1-xyzzyaabz1(xyzzyaabx1)+1
do xyzzyaaaf59=1,xyzzyaacn1
if(.not.any(xyzzyaabz1(:)==heg_nele(xyzzyaaaf59)+1))then
call wordwrap('For the Mahan wave function, we currently require that &
&magic numbers (closed shells) of particles are used.  The number of p&
&articles of each type (e.g. up-spin electrons) must be equal to the n&
&umber of k-vectors in the lowest-lying (shortest) stars.  For the cur&
&rent lattice type the possible numbers of particles are:')
call write_list_int(xyzzyaabx1-1,xyzzyaabz1(2:xyzzyaabx1)-1,12,5,1)
call errstop('EVAL_LATTICE','Stopping due to unfilled shell(s) of elec&
&trons.')
endif
enddo
allocate(xyzzyaacc1(nemax))
xyzzyaaac59=0
xyzzyaaad59=0
do xyzzyaaae59=1,xyzzyaabx1
if(xyzzyaaad59<nemax)then
xyzzyaaac59=xyzzyaaad59+1
xyzzyaaad59=xyzzyaaad59+xyzzyaaca1(xyzzyaaae59)
xyzzyaacc1(xyzzyaaac59:xyzzyaaad59)=xyzzyaaae59
endif
enddo
if(xyzzyaaco1)then
call wout('Debugging info for eval_lattice')
call wout("n_vec_in_star(:) is: ",xyzzyaaca1(:))
call wout("nemax is: ",nemax)
call wout("ivec_in_star is: ",xyzzyaacc1(:))
call wout('num_g_latt is '//trim(i2s(xyzzyaabw1)))
call wout('nstar_latt is '//trim(i2s(xyzzyaabx1)))
call wout('Now printing out first_pr_in_star')
do xyzzyaaaa59=1,xyzzyaabx1
call wout(trim(i2s(xyzzyaabz1(xyzzyaaaa59))))
enddo
call wout('Now printing out pr_lattice')
do xyzzyaaaa59=1,xyzzyaabw1
call wout('pr_lattice(1:3,iii)',xyzzyaace1(1:3,xyzzyaaaa59),rfmt='(e22&
&.12)')
enddo
call wout('Now printing out pr_modsq')
do xyzzyaaaa59=1,xyzzyaabx1
char80=r2s(xyzzyaacd1(xyzzyaaaa59),'(f12.6)')
call wout(trim(char80))
enddo
endif
end subroutine xyzzyaadd1
subroutine xyzzyaade1(n,coords,eta_g_mat,u_g_mat,val,fd,sd,orbmat,grad&
&mat,lapmat,d2mat)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: coords(dimensionality,n),eta_g_mat(3,xyzzyaabv1&
&,n),u_g_mat(3,xyzzyaabv1,n)
logical,intent(in) :: val,fd,sd
complex(dp),intent(inout) :: orbmat(nemax,nemax),gradmat(dimensionalit&
&y,nemax,nemax),lapmat(nemax,nemax),d2mat(xyzzyaabu1,nemax,nemax)
integer xyzzyaaaa60,xyzzyaaab60,xyzzyaaac60
complex(dp) xyzzyaaad60,xyzzyaaae60(dimensionality),xyzzyaaaf60(xyzzya&
&abu1)
call timer('EVAL_ORBMAT',.true.)
orbmat=czero
if(fd)gradmat=czero
if(sd)then
if(.not.use_backflow)then
lapmat=czero
else
d2mat=czero
endif
endif
do xyzzyaaaa60=1,n
xyzzyaaac60=xyzzyaacc1(xyzzyaaaa60)
do xyzzyaaab60=1,n
call xyzzyaadl1(xyzzyaaaa60,coords(1,xyzzyaaab60),eta_g_mat(1,xyzzyaaa&
&c60,xyzzyaaab60),u_g_mat(1,xyzzyaaac60,xyzzyaaab60),fd,sd,xyzzyaaad60&
&,xyzzyaaae60,xyzzyaaaf60)
orbmat(xyzzyaaaa60,xyzzyaaab60)=xyzzyaaad60
if(fd)gradmat(:,xyzzyaaaa60,xyzzyaaab60)=xyzzyaaae60(:)
if(sd)then
if(.not.use_backflow)then
select case(dimensionality)
case(3)
lapmat(xyzzyaaaa60,xyzzyaaab60)=sum(xyzzyaaaf60(1:3))
case(2)
lapmat(xyzzyaaaa60,xyzzyaaab60)=sum(xyzzyaaaf60(1:2))
end select
else
d2mat(:,xyzzyaaaa60,xyzzyaaab60)=xyzzyaaaf60(:)
endif
endif
enddo
enddo
call timer('EVAL_ORBMAT',.false.)
end subroutine xyzzyaade1
subroutine xyzzyaadf1(n,logdet,gradcol,orbbarcol,loggrad)
implicit none
integer,intent(in) :: n
complex(dp),intent(in) :: logdet,gradcol(dimensionality,nemax),orbbarc&
&ol(nemax)
complex(dp),intent(out) :: loggrad(dimensionality)
integer xyzzyaaaa61
do xyzzyaaaa61=1,dimensionality
loggrad(xyzzyaaaa61)=zdotu(n,gradcol(xyzzyaaaa61,1),dimensionality,orb&
&barcol(1),1)
enddo
end subroutine xyzzyaadf1
subroutine xyzzyaadg1(logdet,gradcol,orbbarcol,loggrad)
implicit none
complex(dp),intent(in) :: logdet,gradcol(dimensionality,nemax,nemax,xy&
&zzyaacn1),orbbarcol(nemax,nemax,xyzzyaacn1)
complex(dp),intent(out) :: loggrad(dimensionality)
integer xyzzyaaaa62,xyzzyaaab62,xyzzyaaac62,xyzzyaaad62
loggrad(:)=czero
do xyzzyaaaa62=1,dimensionality
do xyzzyaaab62=1,xyzzyaacn1
do xyzzyaaac62=1,heg_nele(xyzzyaaab62)
do xyzzyaaad62=1,heg_nele(xyzzyaaab62)
loggrad(xyzzyaaaa62)=loggrad(xyzzyaaaa62)-gradcol(xyzzyaaaa62,xyzzyaaa&
&c62,xyzzyaaad62,xyzzyaaab62)*orbbarcol(xyzzyaaac62,xyzzyaaad62,xyzzya&
&aab62)
enddo
enddo
enddo
enddo
end subroutine xyzzyaadg1
subroutine xyzzyaadh1(n,lapcol,orbbarcol,loggrad,lap)
implicit none
complex(dp),intent(in) :: lapcol(nemax),orbbarcol(nemax),loggrad(3)
complex(dp),intent(out) :: lap
integer n
lap=zdotu(n,lapcol(1),1,orbbarcol(1),1)
if(dimensionality==3)then
lap=lap-(dble(loggrad(1))**2+dble(loggrad(2))**2+dble(loggrad(3))**2+a&
&imag(loggrad(1))**2+aimag(loggrad(2))**2+aimag(loggrad(3))**2)
elseif(dimensionality==2)then
lap=lap-(dble(loggrad(1))**2+dble(loggrad(2))**2+aimag(loggrad(1))**2+&
&aimag(loggrad(2))**2)
else
call errstop('EVAL_LAP','dimensionality not implemented')
endif
end subroutine xyzzyaadh1
subroutine xyzzyaadi1(gradcol,orbbarcol,lapcol,loggrad,loglap)
implicit none
complex(dp),intent(in) :: gradcol(dimensionality,nemax,nemax,xyzzyaacn&
&1),orbbarcol(nemax,nemax,xyzzyaacn1),lapcol(nemax,nemax,xyzzyaacn1),l&
&oggrad(3)
complex(dp),intent(out) :: loglap
complex(dp) xyzzyaaaa64(dimensionality,xyzzyaacn1),xyzzyaaab64(dimensi&
&onality,xyzzyaacn1),xyzzyaaac64(xyzzyaacn1),xyzzyaaad64(nemax,nemax)
integer xyzzyaaae64,xyzzyaaaf64,xyzzyaaag64,xyzzyaaah64,xyzzyaaai64
xyzzyaaaa64=czero
do xyzzyaaae64=1,dimensionality
do xyzzyaaaf64=1,xyzzyaacn1
do xyzzyaaag64=1,heg_nele(xyzzyaaaf64)
do xyzzyaaah64=1,heg_nele(xyzzyaaaf64)
xyzzyaaaa64(xyzzyaaae64,xyzzyaaaf64)=xyzzyaaaa64(xyzzyaaae64,xyzzyaaaf&
&64)-gradcol(xyzzyaaae64,xyzzyaaag64,xyzzyaaah64,xyzzyaaaf64)*orbbarco&
&l(xyzzyaaag64,xyzzyaaah64,xyzzyaaaf64)
enddo
enddo
enddo
enddo
xyzzyaaab64=czero
do xyzzyaaae64=1,dimensionality
do xyzzyaaaf64=1,xyzzyaacn1
xyzzyaaad64=czero
do xyzzyaaah64=1,heg_nele(xyzzyaaaf64)
do xyzzyaaai64=1,heg_nele(xyzzyaaaf64)
do xyzzyaaag64=1,heg_nele(xyzzyaaaf64)
xyzzyaaad64(xyzzyaaai64,xyzzyaaah64)=xyzzyaaad64(xyzzyaaai64,xyzzyaaah&
&64)+orbbarcol(xyzzyaaag64,xyzzyaaai64,xyzzyaaaf64)*gradcol(xyzzyaaae6&
&4,xyzzyaaag64,xyzzyaaah64,xyzzyaaaf64)
enddo
enddo
enddo
do xyzzyaaah64=1,heg_nele(xyzzyaaaf64)
do xyzzyaaai64=1,heg_nele(xyzzyaaaf64)
xyzzyaaab64(xyzzyaaae64,xyzzyaaaf64)=xyzzyaaab64(xyzzyaaae64,xyzzyaaaf&
&64)+xyzzyaaad64(xyzzyaaah64,xyzzyaaai64)*xyzzyaaad64(xyzzyaaai64,xyzz&
&yaaah64)
enddo
enddo
enddo
enddo
xyzzyaaac64=czero
do xyzzyaaaf64=1,xyzzyaacn1
do xyzzyaaag64=1,heg_nele(xyzzyaaaf64)
do xyzzyaaah64=1,heg_nele(xyzzyaaaf64)
xyzzyaaac64(xyzzyaaaf64)=xyzzyaaac64(xyzzyaaaf64)+lapcol(xyzzyaaag64,x&
&yzzyaaah64,xyzzyaaaf64)*orbbarcol(xyzzyaaag64,xyzzyaaah64,xyzzyaaaf64&
&)
enddo
enddo
enddo
loglap=czero
do xyzzyaaae64=1,dimensionality
loglap=loglap+2.d0*xyzzyaaaa64(xyzzyaaae64,1)*xyzzyaaaa64(xyzzyaaae64,&
&2)
loglap=loglap+sum(xyzzyaaaa64(xyzzyaaae64,:)**2)
loglap=loglap-sum(xyzzyaaab64(xyzzyaaae64,:))
enddo
loglap=loglap+sum(xyzzyaaac64(:))
if(dimensionality==3)then
loglap=loglap-(dble(loggrad(1))**2+dble(loggrad(2))**2+dble(loggrad(3)&
&)**2+aimag(loggrad(1))**2+aimag(loggrad(2))**2+aimag(loggrad(3))**2)
elseif(dimensionality==2)then
loglap=loglap-(dble(loggrad(1))**2+dble(loggrad(2))**2+aimag(loggrad(1&
&))**2+aimag(loggrad(2))**2)
else
call errstop('EVAL_LOGLAP_HOLE','Dimensionality not implemented')
endif
end subroutine xyzzyaadi1
subroutine xyzzyaadj1(gradcol,orbbarcol,farray)
implicit none
complex(dp),intent(in) :: gradcol(dimensionality,nemax,nemax,xyzzyaacn&
&1),orbbarcol(nemax,nemax,xyzzyaacn1)
real(dp),intent(out) :: farray(3,real1_complex2,netot)
complex(dp) xyzzyaaaa65(dimensionality)
integer xyzzyaaab65,xyzzyaaac65,xyzzyaaad65,xyzzyaaae65,xyzzyaaaf65,xy&
&zzyaaag65,xyzzyaaah65
farray=0.d0
do xyzzyaaac65=1,netot
xyzzyaaaf65=which_spin(xyzzyaaac65)
xyzzyaaah65=heg_nele(xyzzyaaaf65)
if(xyzzyaaaf65==3)then
xyzzyaaaa65=czero
do xyzzyaaab65=1,dimensionality
do xyzzyaaag65=1,xyzzyaacn1
do xyzzyaaae65=1,heg_nele(xyzzyaaag65)
do xyzzyaaad65=1,heg_nele(xyzzyaaag65)
xyzzyaaaa65(xyzzyaaab65)=xyzzyaaaa65(xyzzyaaab65)-gradcol(xyzzyaaab65,&
&xyzzyaaae65,xyzzyaaad65,xyzzyaaag65)*orbbarcol(xyzzyaaae65,xyzzyaaad6&
&5,xyzzyaaag65)
enddo
enddo
enddo
enddo
else
xyzzyaaad65=which_ie(xyzzyaaac65)
do xyzzyaaab65=1,dimensionality
xyzzyaaaa65(xyzzyaaab65)=zdotu(xyzzyaaah65,gradcol(xyzzyaaab65,1,xyzzy&
&aaad65,xyzzyaaaf65),dimensionality,orbbarcol(1,xyzzyaaad65,xyzzyaaaf6&
&5),1)
enddo
endif
farray(1:dimensionality,1,xyzzyaaac65)=dble(xyzzyaaaa65)
farray(1:dimensionality,2,xyzzyaaac65)=aimag(xyzzyaaaa65)
enddo
end subroutine xyzzyaadj1
subroutine xyzzyaadk1(gradcol,d2col,orbbarcol,bf_m2,bf_rmap2,farray,ha&
&rray)
implicit none
integer,intent(in) :: bf_m2(netot),bf_rmap2(netot,netot)
real(dp),intent(in) :: farray(3,real1_complex2,netot)
complex(dp),intent(in) :: d2col(xyzzyaabu1,nemax,nemax,xyzzyaacn1),orb&
&barcol(nemax,nemax,xyzzyaacn1),gradcol(dimensionality,nemax,nemax,xyz&
&zyaacn1)
real(dp),intent(out) :: harray(3,3,real1_complex2,netot,netot)
integer xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66,xyzzyaaad66,xyzzyaaae66,xy&
&zzyaaaf66,xyzzyaaag66,xyzzyaaah66,xyzzyaaai66,xyzzyaaaj66,xyzzyaaak66&
&,xyzzyaaal66,xyzzyaaam66,xyzzyaaan66
complex(dp) xyzzyaaao66
complex(dp) xyzzyaaap66(dimensionality,nemax,xyzzyaacn1),xyzzyaaaq66(d&
&imensionality,xyzzyaacn1),xyzzyaaar66(dimensionality,dimensionality,n&
&emax,xyzzyaacn1),xyzzyaaas66(dimensionality,dimensionality,xyzzyaacn1&
&),xyzzyaaat66(dimensionality,nemax,nemax,xyzzyaacn1),xyzzyaaau66(3,ne&
&tot),xyzzyaaav66(dimensionality,dimensionality,nemax,nemax,xyzzyaacn1&
&),xyzzyaaaw66(dimensionality,dimensionality,nemax,xyzzyaacn1),xyzzyaa&
&ax66(dimensionality,dimensionality,xyzzyaacn1)
harray=0.d0
xyzzyaaau66=cmplx(farray(:,1,:),farray(:,2,:),dp)
xyzzyaaav66=czero
xyzzyaaaw66=czero
xyzzyaaaw66=czero
xyzzyaaat66=czero
do xyzzyaaad66=1,xyzzyaacn1
do xyzzyaaag66=1,dimensionality
do xyzzyaaac66=1,heg_nele(xyzzyaaad66)
do xyzzyaaak66=1,heg_nele(xyzzyaaad66)
do xyzzyaaaj66=1,heg_nele(xyzzyaaad66)
xyzzyaaat66(xyzzyaaag66,xyzzyaaac66,xyzzyaaak66,xyzzyaaad66)=xyzzyaaat&
&66(xyzzyaaag66,xyzzyaaac66,xyzzyaaak66,xyzzyaaad66)+gradcol(xyzzyaaag&
&66,xyzzyaaaj66,xyzzyaaac66,xyzzyaaad66)*orbbarcol(xyzzyaaaj66,xyzzyaa&
&ak66,xyzzyaaad66)
enddo
enddo
enddo
enddo
do xyzzyaaag66=1,dimensionality
do xyzzyaaah66=1,dimensionality
do xyzzyaaac66=1,heg_nele(xyzzyaaad66)
do xyzzyaaak66=1,heg_nele(xyzzyaaad66)
xyzzyaaav66(xyzzyaaag66,xyzzyaaah66,xyzzyaaac66,xyzzyaaak66,xyzzyaaad6&
&6)=xyzzyaaat66(xyzzyaaag66,xyzzyaaac66,xyzzyaaak66,xyzzyaaad66)*xyzzy&
&aaat66(xyzzyaaah66,xyzzyaaak66,xyzzyaaac66,xyzzyaaad66)
enddo
xyzzyaaaw66(xyzzyaaag66,xyzzyaaah66,xyzzyaaac66,xyzzyaaad66)=sum(xyzzy&
&aaav66(xyzzyaaag66,xyzzyaaah66,xyzzyaaac66,:,xyzzyaaad66))
enddo
xyzzyaaax66(xyzzyaaag66,xyzzyaaah66,xyzzyaaad66)=sum(xyzzyaaaw66(xyzzy&
&aaag66,xyzzyaaah66,:,xyzzyaaad66))
enddo
enddo
enddo
xyzzyaaap66=czero
xyzzyaaaq66=czero
do xyzzyaaad66=1,xyzzyaacn1
do xyzzyaaag66=1,dimensionality
do xyzzyaaac66=1,heg_nele(xyzzyaaad66)
xyzzyaaap66(xyzzyaaag66,xyzzyaaac66,xyzzyaaad66)=xyzzyaaat66(xyzzyaaag&
&66,xyzzyaaac66,xyzzyaaac66,xyzzyaaad66)
enddo
xyzzyaaaq66(xyzzyaaag66,xyzzyaaad66)=sum(xyzzyaaap66(xyzzyaaag66,:,xyz&
&zyaaad66))
enddo
enddo
xyzzyaaar66=czero
xyzzyaaas66=czero
do xyzzyaaag66=1,dimensionality
do xyzzyaaah66=1,dimensionality
xyzzyaaai66=xyzzyaabt1(xyzzyaaag66,xyzzyaaah66)
do xyzzyaaad66=1,xyzzyaacn1
do xyzzyaaaj66=1,heg_nele(xyzzyaaad66)
do xyzzyaaac66=1,heg_nele(xyzzyaaad66)
xyzzyaaar66(xyzzyaaag66,xyzzyaaah66,xyzzyaaac66,xyzzyaaad66)=xyzzyaaar&
&66(xyzzyaaag66,xyzzyaaah66,xyzzyaaac66,xyzzyaaad66)+d2col(xyzzyaaai66&
&,xyzzyaaaj66,xyzzyaaac66,xyzzyaaad66)*orbbarcol(xyzzyaaaj66,xyzzyaaac&
&66,xyzzyaaad66)
enddo
enddo
xyzzyaaas66(xyzzyaaag66,xyzzyaaah66,xyzzyaaad66)=sum(xyzzyaaar66(xyzzy&
&aaag66,xyzzyaaah66,:,xyzzyaaad66))
enddo
enddo
enddo
do xyzzyaaaa66=1,netot-1
xyzzyaaac66=which_ie(xyzzyaaaa66)
xyzzyaaad66=which_spin(xyzzyaaaa66)
do xyzzyaaan66=1,bf_m2(xyzzyaaaa66)
xyzzyaaab66=bf_rmap2(xyzzyaaan66,xyzzyaaaa66)
if(xyzzyaaab66==netot.or.xyzzyaaaa66<xyzzyaaab66)cycle
xyzzyaaae66=which_ie(xyzzyaaab66)
xyzzyaaaf66=which_spin(xyzzyaaab66)
do xyzzyaaag66=1,dimensionality
do xyzzyaaah66=1,dimensionality
xyzzyaaai66=xyzzyaabt1(xyzzyaaag66,xyzzyaaah66)
if(xyzzyaaaa66==xyzzyaaab66)then
xyzzyaaao66=xyzzyaaar66(xyzzyaaag66,xyzzyaaah66,xyzzyaaac66,xyzzyaaad6&
&6)
else
if(xyzzyaaad66==xyzzyaaaf66)then
xyzzyaaao66=xyzzyaaau66(xyzzyaaag66,xyzzyaaaa66)*xyzzyaaau66(xyzzyaaah&
&66,xyzzyaaab66)
xyzzyaaao66=xyzzyaaao66-xyzzyaaav66(xyzzyaaag66,xyzzyaaah66,xyzzyaaac6&
&6,xyzzyaaae66,xyzzyaaad66)
else
xyzzyaaao66=xyzzyaaau66(xyzzyaaag66,xyzzyaaaa66)*xyzzyaaau66(xyzzyaaah&
&66,xyzzyaaab66)
endif
endif
harray(xyzzyaaag66,xyzzyaaah66,1,xyzzyaaaa66,xyzzyaaab66)=dble(xyzzyaa&
&ao66)
harray(xyzzyaaag66,xyzzyaaah66,2,xyzzyaaaa66,xyzzyaaab66)=aimag(xyzzya&
&aao66)
if(xyzzyaaaa66/=xyzzyaaab66)then
harray(xyzzyaaah66,xyzzyaaag66,1,xyzzyaaab66,xyzzyaaaa66)=dble(xyzzyaa&
&ao66)
harray(xyzzyaaah66,xyzzyaaag66,2,xyzzyaaab66,xyzzyaaaa66)=aimag(xyzzya&
&aao66)
endif
enddo
enddo
enddo
enddo
xyzzyaaab66=netot
xyzzyaaae66=which_ie(xyzzyaaab66)
xyzzyaaaf66=which_spin(xyzzyaaab66)
do xyzzyaaam66=1,bf_m2(xyzzyaaab66)
xyzzyaaaa66=bf_rmap2(xyzzyaaam66,xyzzyaaab66)
if(xyzzyaaaa66==xyzzyaaab66)cycle
xyzzyaaac66=which_ie(xyzzyaaaa66)
xyzzyaaad66=which_spin(xyzzyaaaa66)
if(xyzzyaaad66==1)then
xyzzyaaal66=2
else
xyzzyaaal66=1
endif
do xyzzyaaag66=1,dimensionality
do xyzzyaaah66=1,dimensionality
xyzzyaaai66=xyzzyaabt1(xyzzyaaag66,xyzzyaaah66)
xyzzyaaao66=czero
xyzzyaaao66=xyzzyaaao66-xyzzyaaap66(xyzzyaaag66,xyzzyaaac66,xyzzyaaad6&
&6)*xyzzyaaaq66(xyzzyaaah66,xyzzyaaal66)
xyzzyaaao66=xyzzyaaao66-xyzzyaaap66(xyzzyaaag66,xyzzyaaac66,xyzzyaaad6&
&6)*xyzzyaaaq66(xyzzyaaah66,xyzzyaaad66)+xyzzyaaaw66(xyzzyaaag66,xyzzy&
&aaah66,xyzzyaaac66,xyzzyaaad66)
xyzzyaaao66=xyzzyaaao66-xyzzyaaar66(xyzzyaaag66,xyzzyaaah66,xyzzyaaac6&
&6,xyzzyaaad66)
harray(xyzzyaaag66,xyzzyaaah66,1,xyzzyaaaa66,netot)=dble(xyzzyaaao66)
harray(xyzzyaaag66,xyzzyaaah66,2,xyzzyaaaa66,netot)=aimag(xyzzyaaao66)
harray(xyzzyaaah66,xyzzyaaag66,1,netot,xyzzyaaaa66)=dble(xyzzyaaao66)
harray(xyzzyaaah66,xyzzyaaag66,2,netot,xyzzyaaaa66)=aimag(xyzzyaaao66)
enddo
enddo
enddo
do xyzzyaaag66=1,dimensionality
do xyzzyaaah66=1,dimensionality
xyzzyaaao66=czero
do xyzzyaaad66=1,xyzzyaacn1
if(xyzzyaaad66==1)then
xyzzyaaal66=2
else
xyzzyaaal66=1
endif
xyzzyaaao66=xyzzyaaao66+xyzzyaaaq66(xyzzyaaag66,xyzzyaaad66)*xyzzyaaaq&
&66(xyzzyaaah66,xyzzyaaal66)
xyzzyaaao66=xyzzyaaao66+xyzzyaaaq66(xyzzyaaag66,xyzzyaaad66)*xyzzyaaaq&
&66(xyzzyaaah66,xyzzyaaad66)-xyzzyaaax66(xyzzyaaag66,xyzzyaaah66,xyzzy&
&aaad66)
xyzzyaaao66=xyzzyaaao66+xyzzyaaas66(xyzzyaaag66,xyzzyaaah66,xyzzyaaad6&
&6)
enddo
harray(xyzzyaaag66,xyzzyaaah66,1,netot,netot)=dble(xyzzyaaao66)
harray(xyzzyaaag66,xyzzyaaah66,2,netot,netot)=aimag(xyzzyaaao66)
enddo
enddo
end subroutine xyzzyaadk1
subroutine xyzzyaadl1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: r(dimensionality),eta_g_mat(3),u_g_mat(3)
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: f,df(dimensionality),d2f(xyzzyaabu1)
integer ivec
f=0.d0
if(fd.or.sd)df=0.d0
if(sd)d2f=0.d0
if(dimensionality==3)then
select case(xyzzyaabr1)
case(1)
call xyzzyaadv1(ivec,r,fd,sd,f,df,d2f)
case(2)
call xyzzyaadn1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
case(3)
call xyzzyaadm1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
end select
elseif(dimensionality==2)then
select case(xyzzyaabr1)
case(1)
call xyzzyaadw1(ivec,r,fd,sd,f,df,d2f)
case(2)
call xyzzyaadp1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
case(3)
call xyzzyaado1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
end select
endif
end subroutine xyzzyaadl1
subroutine xyzzyaadm1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: r(3),eta_g_mat(1:3),u_g_mat(1:3)
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: f,df(3),d2f(6)
integer ivec
real(dp) xyzzyaaaa68(3),xyzzyaaab68(6),xyzzyaaac68,xyzzyaaad68,xyzzyaa&
&ae68,xyzzyaaaf68,xyzzyaaag68,xyzzyaaah68,xyzzyaaai68,xyzzyaaaj68(3),x&
&yzzyaaak68,xyzzyaaal68,xyzzyaaam68(6),xyzzyaaan68(6),xyzzyaaao68(3),x&
&yzzyaaap68(6),xyzzyaaaq68,xyzzyaaar68,xyzzyaaas68,xyzzyaaat68,xyzzyaa&
&au68(3),xyzzyaaav68(6),xyzzyaaaw68
complex(dp) xyzzyaaax68,xyzzyaaay68(3),xyzzyaaaz68(6)
complex(dp),parameter :: xyzzyaaba68=(0.d0,1.d0)
logical xyzzyaabb68
xyzzyaabb68=fd.or.sd
xyzzyaaaa68=gs_kvec(:,ivec)
if(xyzzyaacc1(ivec)==1)then
xyzzyaaad68=0.d0
xyzzyaaaq68=0.d0
xyzzyaaaf68=10.d0
else
xyzzyaaad68=xyzzyaabk1(2,xyzzyaacc1(ivec))
xyzzyaaaq68=xyzzyaabl1(2,xyzzyaacc1(ivec))
xyzzyaaaf68=sqrt(ddot(3,r(1),1,r(1),1))
xyzzyaaag68=eta_g_mat(1)
xyzzyaaar68=u_g_mat(1)
if(xyzzyaabb68)then
xyzzyaaah68=eta_g_mat(2)
xyzzyaaas68=u_g_mat(2)
if(sd)then
xyzzyaaai68=eta_g_mat(3)
xyzzyaaat68=u_g_mat(3)
endif
endif
xyzzyaaae68=1.d0-xyzzyaaag68
if(xyzzyaabb68)call xyzzyaadu1(xyzzyaaaf68,r,xyzzyaaaj68,xyzzyaaak68,x&
&yzzyaaal68,xyzzyaaam68,xyzzyaaan68)
endif
xyzzyaaac68=ddot(3,xyzzyaaaa68(1),1,r(1),1)
if(sd)then
xyzzyaaab68(1)=xyzzyaaaa68(1)**2
xyzzyaaab68(2)=xyzzyaaaa68(2)**2
xyzzyaaab68(3)=xyzzyaaaa68(3)**2
xyzzyaaab68(4)=xyzzyaaaa68(1)*xyzzyaaaa68(2)
xyzzyaaab68(5)=xyzzyaaaa68(1)*xyzzyaaaa68(3)
xyzzyaaab68(6)=xyzzyaaaa68(2)*xyzzyaaaa68(3)
endif
if(xyzzyaaaf68<xyzzyaaad68)then
xyzzyaaax68=exp(xyzzyaaba68*xyzzyaaac68*xyzzyaaae68)
if(xyzzyaabb68)then
xyzzyaaao68(:)=xyzzyaaah68*xyzzyaaaj68(:)
xyzzyaaay68(:)=xyzzyaaba68*(xyzzyaaaa68(:)*xyzzyaaae68-xyzzyaaac68*xyz&
&zyaaao68(:))*xyzzyaaax68
if(sd)then
xyzzyaaap68(:)=xyzzyaaai68*xyzzyaaam68(:)+xyzzyaaah68*xyzzyaaan68(:)
xyzzyaaaz68(1:3)=xyzzyaaay68(1:3)**2/xyzzyaaax68-xyzzyaaba68*xyzzyaaax&
&68*(xyzzyaaaa68(1:3)*2.d0*xyzzyaaao68(1:3)+xyzzyaaac68*xyzzyaaap68(1:&
&3))
xyzzyaaaz68(4)=xyzzyaaay68(1)*xyzzyaaay68(2)/xyzzyaaax68-xyzzyaaba68*x&
&yzzyaaax68*(xyzzyaaaa68(1)*xyzzyaaao68(2)+xyzzyaaaa68(2)*xyzzyaaao68(&
&1)+xyzzyaaac68*xyzzyaaap68(4))
xyzzyaaaz68(5)=xyzzyaaay68(1)*xyzzyaaay68(3)/xyzzyaaax68-xyzzyaaba68*x&
&yzzyaaax68*(xyzzyaaaa68(1)*xyzzyaaao68(3)+xyzzyaaaa68(3)*xyzzyaaao68(&
&1)+xyzzyaaac68*xyzzyaaap68(5))
xyzzyaaaz68(6)=xyzzyaaay68(2)*xyzzyaaay68(3)/xyzzyaaax68-xyzzyaaba68*x&
&yzzyaaax68*(xyzzyaaaa68(2)*xyzzyaaao68(3)+xyzzyaaaa68(3)*xyzzyaaao68(&
&2)+xyzzyaaac68*xyzzyaaap68(6))
endif
endif
else
xyzzyaaax68=exp(xyzzyaaba68*xyzzyaaac68)
if(xyzzyaabb68)then
xyzzyaaay68(:)=xyzzyaaaa68(:)*xyzzyaaba68*xyzzyaaax68
if(sd)xyzzyaaaz68(:)=-xyzzyaaab68(:)*xyzzyaaax68
endif
endif
if(xyzzyaaaf68<xyzzyaaaq68)then
xyzzyaaaw68=exp(xyzzyaaar68)
f=xyzzyaaaw68*xyzzyaaax68
if(xyzzyaabb68)then
xyzzyaaau68(:)=xyzzyaaas68*xyzzyaaaj68(:)
df(:)=xyzzyaaau68(:)*f+xyzzyaaaw68*xyzzyaaay68(:)
if(sd)then
xyzzyaaav68(:)=xyzzyaaat68*xyzzyaaam68(:)+xyzzyaaas68*xyzzyaaan68(:)
d2f(1:3)=xyzzyaaav68(1:3)*f+xyzzyaaau68(1:3)*df(1:3)+xyzzyaaaw68*xyzzy&
&aaau68(1:3)*xyzzyaaay68(1:3)+xyzzyaaaw68*xyzzyaaaz68(1:3)
d2f(4:6)=xyzzyaaav68(4:6)*f+xyzzyaaaw68*xyzzyaaaz68(4:6)
d2f(4)=d2f(4)+xyzzyaaau68(1)*df(2)+xyzzyaaaw68*xyzzyaaau68(2)*xyzzyaaa&
&y68(1)
d2f(5)=d2f(5)+xyzzyaaau68(1)*df(3)+xyzzyaaaw68*xyzzyaaau68(3)*xyzzyaaa&
&y68(1)
d2f(6)=d2f(6)+xyzzyaaau68(2)*df(3)+xyzzyaaaw68*xyzzyaaau68(3)*xyzzyaaa&
&y68(2)
endif
endif
else
f=xyzzyaaax68
if(xyzzyaabb68)then
df(:)=xyzzyaaay68(:)
if(sd)d2f(:)=xyzzyaaaz68(:)
endif
endif
if(xyzzyaaco1)then
call wout("ivec_in_star(ivec) is: ",xyzzyaacc1(ivec))
call wout("L_G_eta is: ",xyzzyaaad68)
call wout("modr is: ",xyzzyaaaf68)
call wout("eta is: ",xyzzyaaag68)
call wout("onelesseta is: ",xyzzyaaae68)
call wout("L_G_u is: ",xyzzyaaaq68)
call wout("u is: ",xyzzyaaar68)
endif
end subroutine xyzzyaadm1
subroutine xyzzyaadn1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: r(3),eta_g_mat(1:3),u_g_mat(1:3)
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: f,df(3),d2f(6)
integer ivec
real(dp) xyzzyaaaa69(3),xyzzyaaab69(6),xyzzyaaac69,xyzzyaaad69,xyzzyaa&
&ae69,xyzzyaaaf69,xyzzyaaag69,xyzzyaaah69,xyzzyaaai69,xyzzyaaaj69(3),x&
&yzzyaaak69,xyzzyaaal69,xyzzyaaam69(6),xyzzyaaan69(6),xyzzyaaao69(3),x&
&yzzyaaap69(6)
complex(dp),parameter :: xyzzyaaaq69=(0.d0,1.d0)
logical xyzzyaaar69
xyzzyaaar69=fd.or.sd
xyzzyaaaa69=gs_kvec(:,ivec)
if(xyzzyaacc1(ivec)==1)then
xyzzyaaad69=0.d0
xyzzyaaaf69=10.d0
else
xyzzyaaad69=xyzzyaabk1(2,xyzzyaacc1(ivec))
xyzzyaaaf69=sqrt(ddot(3,r(1),1,r(1),1))
xyzzyaaag69=eta_g_mat(1)
if(xyzzyaaar69)then
xyzzyaaah69=eta_g_mat(2)
if(sd)then
xyzzyaaai69=eta_g_mat(3)
endif
endif
xyzzyaaae69=1.d0-xyzzyaaag69
if(xyzzyaaar69)call xyzzyaadu1(xyzzyaaaf69,r,xyzzyaaaj69,xyzzyaaak69,x&
&yzzyaaal69,xyzzyaaam69,xyzzyaaan69)
endif
xyzzyaaac69=ddot(3,xyzzyaaaa69(1),1,r(1),1)
xyzzyaaab69(1)=xyzzyaaaa69(1)**2
xyzzyaaab69(2)=xyzzyaaaa69(2)**2
xyzzyaaab69(3)=xyzzyaaaa69(3)**2
xyzzyaaab69(4)=xyzzyaaaa69(1)*xyzzyaaaa69(2)
xyzzyaaab69(5)=xyzzyaaaa69(1)*xyzzyaaaa69(3)
xyzzyaaab69(6)=xyzzyaaaa69(2)*xyzzyaaaa69(3)
if(xyzzyaaaf69<xyzzyaaad69)then
f=exp(xyzzyaaaq69*xyzzyaaac69*xyzzyaaae69)
if(xyzzyaaar69)then
xyzzyaaao69(:)=xyzzyaaah69*xyzzyaaaj69(:)
df(:)=xyzzyaaaq69*(xyzzyaaaa69(:)*xyzzyaaae69-xyzzyaaac69*xyzzyaaao69(&
&:))*f
if(sd)then
xyzzyaaap69(:)=xyzzyaaai69*xyzzyaaam69(:)+xyzzyaaah69*xyzzyaaan69(:)
d2f(1:3)=df(1:3)**2/f-xyzzyaaaq69*f*(xyzzyaaaa69(1:3)*2.d0*xyzzyaaao69&
&(1:3)+xyzzyaaac69*xyzzyaaap69(1:3))
d2f(4)=df(1)*df(2)/f-xyzzyaaaq69*f*(xyzzyaaaa69(1)*xyzzyaaao69(2)+xyzz&
&yaaaa69(2)*xyzzyaaao69(1)+xyzzyaaac69*xyzzyaaap69(4))
d2f(5)=df(1)*df(3)/f-xyzzyaaaq69*f*(xyzzyaaaa69(1)*xyzzyaaao69(3)+xyzz&
&yaaaa69(3)*xyzzyaaao69(1)+xyzzyaaac69*xyzzyaaap69(5))
d2f(6)=df(2)*df(3)/f-xyzzyaaaq69*f*(xyzzyaaaa69(2)*xyzzyaaao69(3)+xyzz&
&yaaaa69(3)*xyzzyaaao69(2)+xyzzyaaac69*xyzzyaaap69(6))
endif
endif
else
f=exp(xyzzyaaaq69*xyzzyaaac69)
if(xyzzyaaar69)then
df(:)=xyzzyaaaa69(:)*xyzzyaaaq69*f
if(sd)d2f(:)=-xyzzyaaab69(:)*f
endif
endif
if(xyzzyaaco1)then
call wout("ivec_in_star(ivec) is: ",xyzzyaacc1(ivec))
call wout("L_G is: ",xyzzyaaad69)
call wout("modr is: ",xyzzyaaaf69)
call wout("eta is: ",xyzzyaaag69)
call wout("onelesseta is: ",xyzzyaaae69)
endif
end subroutine xyzzyaadn1
subroutine xyzzyaado1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: r(2),eta_g_mat(1:3),u_g_mat(1:3)
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: f,df(2),d2f(3)
integer ivec
real(dp) xyzzyaaaa70(2),xyzzyaaab70(3),xyzzyaaac70,xyzzyaaad70,xyzzyaa&
&ae70,xyzzyaaaf70,xyzzyaaag70,xyzzyaaah70,xyzzyaaai70,xyzzyaaaj70(2),x&
&yzzyaaak70,xyzzyaaal70,xyzzyaaam70(3),xyzzyaaan70(3),xyzzyaaao70(2),x&
&yzzyaaap70(3),xyzzyaaaq70,xyzzyaaar70,xyzzyaaas70,xyzzyaaat70,xyzzyaa&
&au70(2),xyzzyaaav70(3),xyzzyaaaw70
complex(dp) xyzzyaaax70,xyzzyaaay70(2),xyzzyaaaz70(3)
complex(dp),parameter :: xyzzyaaba70=(0.d0,1.d0)
logical xyzzyaabb70
xyzzyaabb70=fd.or.sd
xyzzyaaaa70=gs_kvec(1:2,ivec)
if(xyzzyaacc1(ivec)==1)then
xyzzyaaad70=0.d0
xyzzyaaaq70=0.d0
xyzzyaaaf70=10.d0
else
xyzzyaaad70=xyzzyaabk1(2,xyzzyaacc1(ivec))
xyzzyaaaq70=xyzzyaabl1(2,xyzzyaacc1(ivec))
xyzzyaaaf70=sqrt(ddot(2,r(1),1,r(1),1))
xyzzyaaag70=eta_g_mat(1)
xyzzyaaar70=u_g_mat(1)
if(xyzzyaabb70)then
xyzzyaaah70=eta_g_mat(2)
xyzzyaaas70=u_g_mat(2)
if(sd)then
xyzzyaaai70=eta_g_mat(3)
xyzzyaaat70=u_g_mat(3)
endif
endif
xyzzyaaae70=1.d0-xyzzyaaag70
if(xyzzyaabb70)call xyzzyaadu1(xyzzyaaaf70,r,xyzzyaaaj70,xyzzyaaak70,x&
&yzzyaaal70,xyzzyaaam70,xyzzyaaan70)
endif
xyzzyaaac70=ddot(2,xyzzyaaaa70(1),1,r(1),1)
if(sd)then
xyzzyaaab70(1)=xyzzyaaaa70(1)**2
xyzzyaaab70(2)=xyzzyaaaa70(2)**2
xyzzyaaab70(3)=xyzzyaaaa70(1)*xyzzyaaaa70(2)
endif
if(xyzzyaaaf70<xyzzyaaad70)then
xyzzyaaax70=exp(xyzzyaaba70*xyzzyaaac70*xyzzyaaae70)
if(xyzzyaabb70)then
xyzzyaaao70(:)=xyzzyaaah70*xyzzyaaaj70(:)
xyzzyaaay70(:)=xyzzyaaba70*(xyzzyaaaa70(:)*xyzzyaaae70-xyzzyaaac70*xyz&
&zyaaao70(:))*xyzzyaaax70
if(sd)then
xyzzyaaap70(:)=xyzzyaaai70*xyzzyaaam70(:)+xyzzyaaah70*xyzzyaaan70(:)
xyzzyaaaz70(1:2)=xyzzyaaay70(1:2)**2/xyzzyaaax70-xyzzyaaba70*xyzzyaaax&
&70*(xyzzyaaaa70(1:2)*2.d0*xyzzyaaao70(1:2)+xyzzyaaac70*xyzzyaaap70(1:&
&2))
xyzzyaaaz70(3)=xyzzyaaay70(1)*xyzzyaaay70(2)/xyzzyaaax70-xyzzyaaba70*x&
&yzzyaaax70*(xyzzyaaaa70(1)*xyzzyaaao70(2)+xyzzyaaaa70(2)*xyzzyaaao70(&
&1)+xyzzyaaac70*xyzzyaaap70(3))
endif
endif
else
xyzzyaaax70=exp(xyzzyaaba70*xyzzyaaac70)
if(xyzzyaabb70)then
xyzzyaaay70(:)=xyzzyaaaa70(:)*xyzzyaaba70*xyzzyaaax70
if(sd)xyzzyaaaz70(:)=-xyzzyaaab70(:)*xyzzyaaax70
endif
endif
if(xyzzyaaaf70<xyzzyaaaq70)then
xyzzyaaaw70=exp(xyzzyaaar70)
f=xyzzyaaaw70*xyzzyaaax70
if(xyzzyaabb70)then
xyzzyaaau70(:)=xyzzyaaas70*xyzzyaaaj70(:)
df(:)=xyzzyaaau70(:)*f+xyzzyaaaw70*xyzzyaaay70(:)
if(sd)then
xyzzyaaav70(:)=xyzzyaaat70*xyzzyaaam70(:)+xyzzyaaas70*xyzzyaaan70(:)
d2f(1:2)=xyzzyaaav70(1:2)*f+xyzzyaaau70(1:2)*df(1:2)+xyzzyaaaw70*xyzzy&
&aaau70(1:2)*xyzzyaaay70(1:2)+xyzzyaaaw70*xyzzyaaaz70(1:2)
d2f(3)=xyzzyaaav70(3)*f+xyzzyaaaw70*xyzzyaaaz70(3)+xyzzyaaau70(1)*df(2&
&)+xyzzyaaaw70*xyzzyaaau70(2)*xyzzyaaay70(1)
endif
endif
else
f=xyzzyaaax70
if(xyzzyaabb70)then
df(:)=xyzzyaaay70(:)
if(sd)d2f(:)=xyzzyaaaz70(:)
endif
endif
if(xyzzyaaco1)then
call wout("ivec_in_star(ivec) is: ",xyzzyaacc1(ivec))
call wout("L_G_eta is: ",xyzzyaaad70)
call wout("modr is: ",xyzzyaaaf70)
call wout("eta is: ",xyzzyaaag70)
call wout("onelesseta is: ",xyzzyaaae70)
call wout("L_G_u is: ",xyzzyaaaq70)
call wout("u is: ",xyzzyaaar70)
endif
end subroutine xyzzyaado1
subroutine xyzzyaadp1(ivec,r,eta_g_mat,u_g_mat,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: r(2),eta_g_mat(1:3),u_g_mat(1:3)
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: f,df(2),d2f(3)
integer ivec
real(dp) xyzzyaaaa71(2),xyzzyaaab71(3),xyzzyaaac71,xyzzyaaad71,xyzzyaa&
&ae71,xyzzyaaaf71,xyzzyaaag71,xyzzyaaah71,xyzzyaaai71,xyzzyaaaj71(2),x&
&yzzyaaak71,xyzzyaaal71,xyzzyaaam71(3),xyzzyaaan71(3),xyzzyaaao71(2),x&
&yzzyaaap71(3)
logical xyzzyaaaq71
complex(dp),parameter :: xyzzyaaar71=(0.d0,1.d0)
xyzzyaaaq71=fd.or.sd
xyzzyaaaa71=gs_kvec(1:2,ivec)
if(xyzzyaacc1(ivec)==1)then
xyzzyaaad71=0.d0
xyzzyaaaf71=10.d0
else
xyzzyaaad71=xyzzyaabk1(2,xyzzyaacc1(ivec))
xyzzyaaaf71=sqrt(ddot(2,r(1),1,r(1),1))
xyzzyaaag71=eta_g_mat(1)
if(xyzzyaaaq71)then
xyzzyaaah71=eta_g_mat(2)
if(sd)then
xyzzyaaai71=eta_g_mat(3)
endif
endif
xyzzyaaae71=1.d0-xyzzyaaag71
if(xyzzyaaaq71)call xyzzyaadu1(xyzzyaaaf71,r,xyzzyaaaj71,xyzzyaaak71,x&
&yzzyaaal71,xyzzyaaam71,xyzzyaaan71)
endif
xyzzyaaac71=ddot(2,xyzzyaaaa71(1),1,r(1),1)
if(sd)then
xyzzyaaab71(1)=xyzzyaaaa71(1)**2
xyzzyaaab71(2)=xyzzyaaaa71(2)**2
xyzzyaaab71(3)=xyzzyaaaa71(1)*xyzzyaaaa71(2)
endif
if(xyzzyaaaf71<xyzzyaaad71)then
f=exp(xyzzyaaar71*xyzzyaaac71*xyzzyaaae71)
if(xyzzyaaaq71)then
xyzzyaaao71(:)=xyzzyaaah71*xyzzyaaaj71(:)
df(:)=xyzzyaaar71*(xyzzyaaaa71(:)*xyzzyaaae71-xyzzyaaac71*xyzzyaaao71(&
&:))*f
if(sd)then
xyzzyaaap71(:)=xyzzyaaai71*xyzzyaaam71(:)+xyzzyaaah71*xyzzyaaan71(:)
d2f(1:2)=df(1:2)**2/f-xyzzyaaar71*f*(xyzzyaaaa71(1:2)*2.d0*xyzzyaaao71&
&(1:2)+xyzzyaaac71*xyzzyaaap71(1:2))
d2f(3)=df(1)*df(2)/f-xyzzyaaar71*f*(xyzzyaaaa71(1)*xyzzyaaao71(2)+xyzz&
&yaaaa71(2)*xyzzyaaao71(1)+xyzzyaaac71*xyzzyaaap71(3))
endif
endif
else
f=exp(xyzzyaaar71*xyzzyaaac71)
if(xyzzyaaaq71)then
df(:)=xyzzyaaaa71(:)*xyzzyaaar71*f
if(sd)d2f(:)=-xyzzyaaab71(:)*f
endif
endif
if(xyzzyaaco1)then
call wout("ivec_in_star(ivec) is: ",xyzzyaacc1(ivec))
call wout("L_G_eta is: ",xyzzyaaad71)
call wout("modr is: ",xyzzyaaaf71)
call wout("eta is: ",xyzzyaaag71)
call wout("onelesseta is: ",xyzzyaaae71)
endif
end subroutine xyzzyaadp1
subroutine xyzzyaadq1(n,coords,fsd,sd,eta_g_mat)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: coords(dimensionality,n)
logical,intent(in) :: fsd,sd
real(dp),intent(inout) :: eta_g_mat(3,xyzzyaabv1,n)
integer xyzzyaaaa72,xyzzyaaab72
real(dp) xyzzyaaac72,xyzzyaaad72,xyzzyaaae72,xyzzyaaaf72
eta_g_mat=0.d0
do xyzzyaaab72=1,n
xyzzyaaac72=sqrt(ddot(dimensionality,coords(1,xyzzyaaab72),1,coords(1,&
&xyzzyaaab72),1))
do xyzzyaaaa72=2,xyzzyaabv1
call xyzzyaads1(xyzzyaaac72,xyzzyaaaa72,fsd,sd,xyzzyaaad72,xyzzyaaae72&
&,xyzzyaaaf72)
eta_g_mat(:,xyzzyaaaa72,xyzzyaaab72)=(/xyzzyaaad72,xyzzyaaae72,xyzzyaa&
&af72/)
enddo
enddo
end subroutine xyzzyaadq1
subroutine xyzzyaadr1(n,coords,fsd,sd,u_g_mat)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: coords(dimensionality,n)
logical,intent(in) :: fsd,sd
real(dp),intent(inout) :: u_g_mat(3,xyzzyaabv1,n)
integer xyzzyaaaa73,xyzzyaaab73
real(dp) xyzzyaaac73,xyzzyaaad73,xyzzyaaae73,xyzzyaaaf73
u_g_mat=0.d0
do xyzzyaaab73=1,n
xyzzyaaac73=sqrt(ddot(dimensionality,coords(1,xyzzyaaab73),1,coords(1,&
&xyzzyaaab73),1))
do xyzzyaaaa73=2,xyzzyaabv1
call xyzzyaadt1(xyzzyaaac73,xyzzyaaaa73,fsd,sd,xyzzyaaad73,xyzzyaaae73&
&,xyzzyaaaf73)
u_g_mat(:,xyzzyaaaa73,xyzzyaaab73)=(/xyzzyaaad73,xyzzyaaae73,xyzzyaaaf&
&73/)
enddo
enddo
end subroutine xyzzyaadr1
subroutine xyzzyaads1(r,istar,fsd,sd,f,df,d2f)
implicit none
integer,intent(in) :: istar
real(dp),intent(in) :: r
logical,intent(in) :: fsd,sd
real(dp),intent(out) :: f,df,d2f
integer xyzzyaaaa74,xyzzyaaab74
real(dp) xyzzyaaac74,xyzzyaaad74,xyzzyaaae74,xyzzyaaaf74,xyzzyaaag74,x&
&yzzyaaah74
real(dp) xyzzyaaai74(0:xyzzyaabc1-2),xyzzyaaaj74(0:xyzzyaabc1-2),xyzzy&
&aaak74(0:xyzzyaabc1-3),xyzzyaaal74(0:xyzzyaabc1-3)
if(r>xyzzyaabk1(2,istar))then
f=0.d0
df=0.d0
d2f=0.d0
return
endif
if(istar==1)then
f=0.d0
df=0.d0
d2f=0.d0
call errstop('ETA_G','istar should not be 1')
endif
xyzzyaaaa74=xyzzyaabc1-3
xyzzyaaae74=xyzzyaabk1(1,istar)
xyzzyaaac74=xyzzyaabk1(2,istar)
xyzzyaaaf74=1.d0/xyzzyaaac74
xyzzyaaai74(0:xyzzyaaaa74)=xyzzyaabk1(3:xyzzyaaaa74+3,istar)
xyzzyaaai74(xyzzyaaaa74+1)=0.d0
xyzzyaaad74=(1.d0-r*xyzzyaaaf74)**xyzzyaaae74
xyzzyaaag74=1.d0
do xyzzyaaab74=0,xyzzyaaaa74
xyzzyaaal74(xyzzyaaab74)=xyzzyaaag74
xyzzyaaag74=xyzzyaaag74*r
enddo
if(fsd)then
xyzzyaaaj74(0:xyzzyaaaa74)=(-xyzzyaaae74*xyzzyaaaf74)*xyzzyaaai74(0:xy&
&zzyaaaa74)
xyzzyaaah74=0.d0
do xyzzyaaab74=0,xyzzyaaaa74
xyzzyaaah74=xyzzyaaah74+(1.d0-r*xyzzyaaaf74)
xyzzyaaaj74(xyzzyaaab74)=xyzzyaaaj74(xyzzyaaab74)+xyzzyaaah74*xyzzyaaa&
&i74(xyzzyaaab74+1)
enddo
xyzzyaaaj74(xyzzyaaaa74+1)=0.d0
endif
if(sd)then
xyzzyaaak74(0:xyzzyaaaa74)=((1.d0-xyzzyaaae74)*xyzzyaaaf74)*xyzzyaaaj7&
&4(0:xyzzyaaaa74)
xyzzyaaah74=0.d0
do xyzzyaaab74=0,xyzzyaaaa74
xyzzyaaah74=xyzzyaaah74+(1.d0-r*xyzzyaaaf74)
xyzzyaaak74(xyzzyaaab74)=xyzzyaaak74(xyzzyaaab74)+xyzzyaaah74*(xyzzyaa&
&aj74(xyzzyaaab74+1)-xyzzyaaai74(xyzzyaaab74+1)*xyzzyaaaf74)
enddo
endif
f=xyzzyaaad74*ddot(xyzzyaaaa74+1,xyzzyaaai74(0),1,xyzzyaaal74(0),1)
if(fsd)then
df=(1.d0-r*xyzzyaaaf74)**(xyzzyaaae74-1.d0)*ddot(xyzzyaaaa74+1,xyzzyaa&
&aj74(0),1,xyzzyaaal74(0),1)
else
df=0.d0
endif
if(sd)then
d2f=(1.d0-r*xyzzyaaaf74)**(xyzzyaaae74-2.d0)*ddot(xyzzyaaaa74+1,xyzzya&
&aak74(0),1,xyzzyaaal74(0),1)
else
d2f=0.d0
endif
end subroutine xyzzyaads1
subroutine xyzzyaadt1(r,istar,fsd,sd,f,df,d2f)
implicit none
integer,intent(in) :: istar
real(dp),intent(in) :: r
logical,intent(in) :: fsd,sd
real(dp),intent(out) :: f,df,d2f
integer xyzzyaaaa75,xyzzyaaab75
real(dp) xyzzyaaac75,xyzzyaaad75,xyzzyaaae75,xyzzyaaaf75,xyzzyaaag75,x&
&yzzyaaah75
real(dp) xyzzyaaai75(0:xyzzyaabd1-2),xyzzyaaaj75(0:xyzzyaabd1-2),xyzzy&
&aaak75(0:xyzzyaabd1-3),xyzzyaaal75(0:xyzzyaabd1-3)
if(r>xyzzyaabl1(2,istar))then
f=0.d0
df=0.d0
d2f=0.d0
return
endif
if(istar==1)then
f=0.d0
df=0.d0
d2f=0.d0
call errstop('U_G','istar should not be 1')
endif
xyzzyaaaa75=xyzzyaabd1-3
xyzzyaaae75=xyzzyaabl1(1,istar)
xyzzyaaac75=xyzzyaabl1(2,istar)
xyzzyaaaf75=1.d0/xyzzyaaac75
xyzzyaaai75(0:xyzzyaaaa75)=xyzzyaabl1(3:xyzzyaaaa75+3,istar)
xyzzyaaai75(xyzzyaaaa75+1)=0.d0
xyzzyaaad75=(1.d0-r*xyzzyaaaf75)**xyzzyaaae75
xyzzyaaag75=1.d0
do xyzzyaaab75=0,xyzzyaaaa75
xyzzyaaal75(xyzzyaaab75)=xyzzyaaag75
xyzzyaaag75=xyzzyaaag75*r
enddo
if(fsd)then
xyzzyaaaj75(0:xyzzyaaaa75)=(-xyzzyaaae75*xyzzyaaaf75)*xyzzyaaai75(0:xy&
&zzyaaaa75)
xyzzyaaah75=0.d0
do xyzzyaaab75=0,xyzzyaaaa75
xyzzyaaah75=xyzzyaaah75+(1.d0-r*xyzzyaaaf75)
xyzzyaaaj75(xyzzyaaab75)=xyzzyaaaj75(xyzzyaaab75)+xyzzyaaah75*xyzzyaaa&
&i75(xyzzyaaab75+1)
enddo
xyzzyaaaj75(xyzzyaaaa75+1)=0.d0
endif
if(sd)then
xyzzyaaak75(0:xyzzyaaaa75)=((1.d0-xyzzyaaae75)*xyzzyaaaf75)*xyzzyaaaj7&
&5(0:xyzzyaaaa75)
xyzzyaaah75=0.d0
do xyzzyaaab75=0,xyzzyaaaa75
xyzzyaaah75=xyzzyaaah75+(1.d0-r*xyzzyaaaf75)
xyzzyaaak75(xyzzyaaab75)=xyzzyaaak75(xyzzyaaab75)+xyzzyaaah75*(xyzzyaa&
&aj75(xyzzyaaab75+1)-xyzzyaaai75(xyzzyaaab75+1)*xyzzyaaaf75)
enddo
endif
f=xyzzyaaad75*ddot(xyzzyaaaa75+1,xyzzyaaai75(0),1,xyzzyaaal75(0),1)
if(fsd)then
df=(1.d0-r*xyzzyaaaf75)**(xyzzyaaae75-1.d0)*ddot(xyzzyaaaa75+1,xyzzyaa&
&aj75(0),1,xyzzyaaal75(0),1)
else
df=0.d0
endif
if(sd)then
d2f=(1.d0-r*xyzzyaaaf75)**(xyzzyaaae75-2.d0)*ddot(xyzzyaaaa75+1,xyzzya&
&aak75(0),1,xyzzyaaal75(0),1)
else
d2f=0.d0
endif
end subroutine xyzzyaadt1
subroutine xyzzyaadu1(r,vecr,grad,grad2,lap,drdr,d2r)
implicit none
real(dp),intent(in) :: r,vecr(:)
real(dp),intent(out) :: grad(dimensionality),grad2,lap,drdr(xyzzyaabu1&
&),d2r(xyzzyaabu1)
real(dp) xyzzyaaaa76
grad2=1.d0
if(r==0.d0)then
grad(1:dimensionality)=0.d0
lap=0.d0
drdr=0.d0
d2r=0.d0
return
endif
select case(dimensionality)
case(3)
xyzzyaaaa76=1.d0/r
lap=2*xyzzyaaaa76
grad(1:3)=vecr(1:3)*xyzzyaaaa76
drdr(1:3)=grad(1:3)*grad(1:3)
drdr(4:6)=(/grad(1)*grad(2),grad(1)*grad(3),grad(2)*grad(3)/)
d2r(1:3)=xyzzyaaaa76-drdr(1:3)*xyzzyaaaa76
d2r(4:6)=-drdr(4:6)*xyzzyaaaa76
case(2)
xyzzyaaaa76=1.d0/r
lap=xyzzyaaaa76
grad(1:2)=vecr(1:2)*xyzzyaaaa76
drdr(1:2)=grad(1:2)*grad(1:2)
drdr(3)=grad(1)*grad(2)
d2r(1:2)=xyzzyaaaa76-drdr(1:2)*xyzzyaaaa76
d2r(3)=-drdr(3)*xyzzyaaaa76
case(1)
end select
end subroutine xyzzyaadu1
subroutine xyzzyaadv1(ivec,r,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: r(3)
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: f,df(3),d2f(6)
integer ivec
real(dp) xyzzyaaaa77(3),xyzzyaaab77(6),xyzzyaaac77
logical xyzzyaaad77
xyzzyaaad77=fd.or.sd
xyzzyaaaa77=gs_kvec(:,ivec)
if(xyzzyaaco1)then
call wout("ivec is: ",ivec)
call wout("Choose which orbitals to occupy by hand in mahan_ex:eval_f_&
&plane_wave.")
if(ivec==1)then
xyzzyaaaa77=xyzzyaace1(:,1)
elseif(ivec==2)then
xyzzyaaaa77=xyzzyaace1(:,7)
elseif(ivec==3)then
xyzzyaaaa77=xyzzyaace1(:,3)
elseif(ivec==4)then
xyzzyaaaa77=xyzzyaace1(:,5)
elseif(ivec==5)then
xyzzyaaaa77=xyzzyaace1(:,2)
elseif(ivec==6)then
xyzzyaaaa77=xyzzyaace1(:,4)
elseif(ivec==7)then
xyzzyaaaa77=xyzzyaace1(:,6)
elseif(ivec==8)then
xyzzyaaaa77=xyzzyaace1(:,11)
elseif(ivec==9)then
xyzzyaaaa77=xyzzyaace1(:,17)
endif
call wout("k - gs_kvec(:,ivec) is: ",xyzzyaaaa77-gs_kvec(:,ivec))
endif
xyzzyaaac77=ddot(3,xyzzyaaaa77(1),1,r(1),1)
xyzzyaaab77(1)=xyzzyaaaa77(1)**2
xyzzyaaab77(2)=xyzzyaaaa77(2)**2
xyzzyaaab77(3)=xyzzyaaaa77(3)**2
xyzzyaaab77(4)=xyzzyaaaa77(1)*xyzzyaaaa77(2)
xyzzyaaab77(5)=xyzzyaaaa77(1)*xyzzyaaaa77(3)
xyzzyaaab77(6)=xyzzyaaaa77(2)*xyzzyaaaa77(3)
f=cmplx(cos(xyzzyaaac77),sin(xyzzyaaac77),dp)
if(xyzzyaaad77)then
df(:)=xyzzyaaaa77(:)*cmplx(-sin(xyzzyaaac77),cos(xyzzyaaac77),dp)
if(sd)d2f(:)=-xyzzyaaab77(:)*cmplx(cos(xyzzyaaac77),sin(xyzzyaaac77),d&
&p)
endif
if(xyzzyaaco1)then
call wout("k is: ",xyzzyaaaa77)
call wout("r is: ",r)
call wout("f is: ",f)
endif
end subroutine xyzzyaadv1
subroutine xyzzyaadw1(ivec,r,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: r(2)
logical,intent(in) :: fd,sd
complex(dp),intent(out) :: f,df(2),d2f(3)
integer ivec
real(dp) xyzzyaaaa78(2),xyzzyaaab78(3),xyzzyaaac78
logical xyzzyaaad78
xyzzyaaad78=fd.or.sd
xyzzyaaaa78=gs_kvec(1:2,ivec)
xyzzyaaac78=ddot(2,xyzzyaaaa78(1),1,r(1),1)
xyzzyaaab78(1)=xyzzyaaaa78(1)**2
xyzzyaaab78(2)=xyzzyaaaa78(2)**2
xyzzyaaab78(3)=xyzzyaaaa78(1)*xyzzyaaaa78(2)
f=cmplx(cos(xyzzyaaac78),sin(xyzzyaaac78),dp)
if(xyzzyaaad78)then
df(:)=xyzzyaaaa78(:)*cmplx(-sin(xyzzyaaac78),cos(xyzzyaaac78),dp)
if(sd)d2f(:)=-xyzzyaaab78(:)*cmplx(cos(xyzzyaaac78),sin(xyzzyaaac78),d&
&p)
endif
if(xyzzyaaco1)then
call wout("k is: ",xyzzyaaaa78)
call wout("r is: ",r)
call wout("f is: ",f)
endif
end subroutine xyzzyaadw1
end module slaarnabn

module slaarnaao
use slaarnaad
use dsp
use slaarnach
use store
use slaarnaag,    only : pi,czero
use slaarnaan,only : ee_kato_gamma
use file_utils,   only : open_units
use format_utils, only : wout,i2s,l2s,r2s2,write_list_int
use slaarnabg,     only : dimensionality
use slaarnabt,    only : dcopy,ddot,exp_protect
use parallel,     only : am_master
use run_control,  only : errstop,errstop_master,check_alloc
implicit none
private
public read_exmol,write_exmol
public query_exmol_levels,query_exmol_level_details,setup_exmol,finish&
&_exmol,wfn_ratio_exmol,accept_move_exmol,reset_config_exmol,wfn_logva&
&l_exmol,wfn_loggrad_exmol,wfn_loglap_exmol,prefetch_wfn_exmol,clear_s&
&cratch_exmol,add_config_exmol_items,setup_exmol_params,finish_exmol_p&
&arams,get_exmol_params,put_exmol_params,clone_scratch_exmol,invalidat&
&e_params_exmol,invalidate_param1_exmol,setup_storage_exmol,finish_sto&
&rage_exmol,load_from_storage_exmol,save_to_storage_exmol,enumerate_pl&
&ot_exmol,query_plot_exmol,get_plot_exmol,finish_plot_exmol
public gen_config_exmol,delete_config_exmol,copy_config_exmol,config_t&
&o_pt_exmol,pt_to_config_exmol,redist_allocations_exmol,redist_load_ex&
&mol,redist_send_exmol,redist_recv_exmol,redist_save_exmol,redist_deal&
&locations_exmol,load_from_pt_exmol,save_to_pt_exmol
public config_wfn_exmol
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,ne,xyzzyaaad1,xyzzyaaae1,xyzz&
&yaaaf1,xyzzyaaag1(2)
integer,allocatable :: xyzzyaaah1(:),xyzzyaaai1(:,:,:),xyzzyaaaj1(:),x&
&yzzyaaak1(:),xyzzyaaal1(:,:),xyzzyaaam1(:),xyzzyaaan1(:,:),xyzzyaaao1&
&(:),xyzzyaaap1(:,:),xyzzyaaaq1(:,:),xyzzyaaar1(:),xyzzyaaas1(:),xyzzy&
&aaat1(:)
real(dp) xyzzyaaau1
real(dp),allocatable :: xyzzyaaav1(:,:),xyzzyaaaw1(:),xyzzyaaax1(:)
logical xyzzyaaay1,xyzzyaaaz1
logical,allocatable :: xyzzyaaba1(:)
character(32),allocatable :: xyzzyaabb1(:)
integer xyzzyaabc1,xyzzyaabd1,xyzzyaabe1
integer xyzzyaabf1
character(80) title
real(dp),allocatable,target :: xyzzyaabg1(:,:,:,:),xyzzyaabh1(:,:,:,:,&
&:),xyzzyaabi1(:,:,:,:),xyzzyaabj1(:,:,:,:,:)
real(dp),allocatable :: xyzzyaabk1(:,:),xyzzyaabl1(:,:),xyzzyaabm1(:),&
&xyzzyaabn1(:,:,:,:),xyzzyaabo1(:,:,:,:,:),xyzzyaabp1(:,:,:,:,:,:)
complex(dp),allocatable :: xyzzyaabq1(:,:,:),xyzzyaabr1(:,:)
logical,allocatable :: xyzzyaabs1(:),xyzzyaabt1(:),xyzzyaabu1(:),xyzzy&
&aabv1(:),xyzzyaabw1(:),xyzzyaabx1(:),xyzzyaaby1(:,:),xyzzyaabz1(:,:),&
&xyzzyaaca1(:),xyzzyaacb1(:)
integer,parameter :: xyzzyaacc1=2
integer xyzzyaacd1,xyzzyaace1(xyzzyaacc1)
integer,allocatable :: xyzzyaacf1(:)
type config_wfn_exmol
private
logical dummy_variable
end type config_wfn_exmol
real(dp),allocatable,target :: xyzzyaacg1(:,:),xyzzyaach1(:,:,:),xyzzy&
&aaci1(:,:),xyzzyaacj1(:,:,:)
real(dp),allocatable :: xyzzyaack1(:,:),xyzzyaacl1(:,:),xyzzyaacm1(:)
complex(dp),allocatable :: xyzzyaacn1(:,:,:),xyzzyaaco1(:,:)
logical,allocatable :: xyzzyaacp1(:),xyzzyaacq1(:),xyzzyaacr1(:),xyzzy&
&aacs1(:)
real(dp),allocatable :: xyzzyaact1(:,:),xyzzyaacu1(:,:,:)
integer xyzzyaacv1
contains
subroutine query_exmol_levels(nlevels)
implicit none
integer,intent(out) :: nlevels
nlevels=1
end subroutine query_exmol_levels
subroutine query_exmol_level_details(level_score,level_name)
implicit none
integer,intent(out) :: level_score(*)
character(80),intent(out) :: level_name(*)
level_score(1)=200
level_name(1)='Excitonic molecule wave function'
end subroutine query_exmol_level_details
subroutine setup_exmol
implicit none
integer xyzzyaaaa4(2),xyzzyaaab4(2),xyzzyaaac4(2),xyzzyaaad4(2),xyzzya&
&aae4(2),xyzzyaaaf4(2),xyzzyaaag4(2),xyzzyaaah4(2),xyzzyaaai4(2),xyzzy&
&aaaj4,xyzzyaaak4
xyzzyaaad4=0
xyzzyaaaa4=0
xyzzyaaab4=0
xyzzyaaac4=0
xyzzyaaae4=0
xyzzyaaaf4=0
xyzzyaaag4=0
xyzzyaaah4=0
xyzzyaaai4=0
call include_range((/1,nscratch/),xyzzyaaad4)
call include_range((/1,nscratch/),xyzzyaaaa4)
call include_range((/1,nscratch/),xyzzyaaab4)
call include_range((/1,nscratch/),xyzzyaaaf4)
call include_range((/1,nscratch/),xyzzyaaag4)
if(.not.use_backflow)then
call include_range((/1,nscratch/),xyzzyaaac4)
else
call include_range((/1,nscratch/),xyzzyaaae4)
call include_range((/1,nscratch/),xyzzyaaah4)
call include_range((/1,nscratch/),xyzzyaaai4)
endif
if(xyzzyaaad4(1)/=0)then
allocate(xyzzyaabk1(xyzzyaaaa1,xyzzyaaad4(1):xyzzyaaad4(2)),xyzzyaabl1&
&(xyzzyaaaa1,xyzzyaaad4(1):xyzzyaaad4(2)),xyzzyaabm1(xyzzyaaad4(1):xyz&
&zyaaad4(2)),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','renorm_term')
xyzzyaabk1=0.d0
xyzzyaabl1=0.d0
xyzzyaabm1=0.d0
endif
allocate(xyzzyaabw1(nscratch),xyzzyaabx1(nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','term')
if(xyzzyaaaa4(1)/=0)then
allocate(xyzzyaabg1(ne,xyzzyaaad1,xyzzyaaac1,xyzzyaaaa4(1):xyzzyaaaa4(&
&2)),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','orbmat')
xyzzyaabg1=0.d0
endif
allocate(xyzzyaabs1(nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','orbmat_valid')
if(xyzzyaaab4(1)/=0)then
allocate(xyzzyaabh1(3,ne,xyzzyaaad1,xyzzyaaac1,xyzzyaaab4(1):xyzzyaaab&
&4(2)),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','gradmat')
xyzzyaabh1=0.d0
endif
allocate(xyzzyaabt1(nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','gradmat_valid')
if(xyzzyaaac4(1)/=0)then
allocate(xyzzyaabi1(ne,xyzzyaaad1,xyzzyaaac1,xyzzyaaac4(1):xyzzyaaac4(&
&2)),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','lapmat')
xyzzyaabi1=0.d0
endif
allocate(xyzzyaabu1(nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','lapmat_valid')
if(xyzzyaaae4(1)/=0)then
allocate(xyzzyaabj1(6,ne,xyzzyaaad1,xyzzyaaac1,xyzzyaaae4(1):xyzzyaaae&
&4(2)),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','d2mat')
xyzzyaabj1=0.d0
endif
allocate(xyzzyaabv1(nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','d2mat_valid')
if(xyzzyaaaf4(1)/=0)then
allocate(xyzzyaabq1(3,xyzzyaaab1,xyzzyaaaf4(1):xyzzyaaaf4(2)),stat=xyz&
&zyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','fi')
xyzzyaabq1=czero
endif
allocate(xyzzyaaby1(xyzzyaaab1,nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','fi_valid')
if(xyzzyaaag4(1)/=0)then
allocate(xyzzyaabr1(xyzzyaaab1,xyzzyaaag4(1):xyzzyaaag4(2)),stat=xyzzy&
&aaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','ti')
xyzzyaabr1=czero
endif
allocate(xyzzyaabz1(xyzzyaaab1,nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','ti_valid')
if(xyzzyaaah4(1)/=0)then
allocate(xyzzyaabn1(3,real1_complex2,netot,xyzzyaaah4(1):xyzzyaaah4(2)&
&),xyzzyaabo1(3,real1_complex2,xyzzyaaaa1,netot,xyzzyaaah4(1):xyzzyaaa&
&h4(2)),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','Farray')
xyzzyaabn1=0.d0
xyzzyaabo1=0.d0
endif
allocate(xyzzyaaca1(nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','Farray_valid')
if(xyzzyaaai4(1)/=0)then
allocate(xyzzyaabp1(3,3,real1_complex2,netot,netot,xyzzyaaai4(1):xyzzy&
&aaai4(2)),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','Harray')
endif
allocate(xyzzyaacb1(nscratch),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','Harray_valid')
allocate(xyzzyaacg1(xyzzyaaae1,xyzzyaaac1),xyzzyaach1(3,xyzzyaaae1,xyz&
&zyaaac1),xyzzyaaci1(xyzzyaaae1,xyzzyaaac1),xyzzyaacj1(6,xyzzyaaae1,xy&
&zzyaaac1),stat=xyzzyaaaj4)
call check_alloc(xyzzyaaaj4,'SETUP_EXMOL','temps')
xyzzyaacg1=0.d0
xyzzyaach1=0.d0
xyzzyaaci1=0.d0
xyzzyaacj1=0.d0
if(use_backflow)call setup_bf
do xyzzyaaak4=1,nscratch
call clear_scratch_exmol(xyzzyaaak4)
enddo
end subroutine setup_exmol
subroutine finish_exmol
implicit none
if(allocated(xyzzyaabk1))deallocate(xyzzyaabk1,xyzzyaabl1,xyzzyaabm1)
deallocate(xyzzyaabw1,xyzzyaabx1)
if(allocated(xyzzyaabg1))deallocate(xyzzyaabg1)
deallocate(xyzzyaabs1)
if(allocated(xyzzyaabh1))deallocate(xyzzyaabh1)
deallocate(xyzzyaabt1)
if(allocated(xyzzyaabi1))deallocate(xyzzyaabi1)
deallocate(xyzzyaabu1)
if(allocated(xyzzyaabj1))deallocate(xyzzyaabj1)
deallocate(xyzzyaabv1)
if(allocated(xyzzyaabq1))deallocate(xyzzyaabq1)
deallocate(xyzzyaaby1)
if(allocated(xyzzyaabr1))deallocate(xyzzyaabr1)
deallocate(xyzzyaabz1)
if(allocated(xyzzyaabn1))deallocate(xyzzyaabn1,xyzzyaabo1)
deallocate(xyzzyaaca1)
if(allocated(xyzzyaabp1))deallocate(xyzzyaabp1)
deallocate(xyzzyaacb1)
deallocate(xyzzyaacg1,xyzzyaach1,xyzzyaaci1,xyzzyaacj1)
if(use_backflow)call finish_bf
end subroutine finish_exmol
subroutine wfn_ratio_exmol(is,js,ilevel,ratio,fd,sd)
implicit none
integer,intent(in) :: is,js,ilevel
complex(dp),intent(out) :: ratio
logical,intent(in) :: fd,sd
real(dp) xyzzyaaaa6
call xyzzyaacy1(is,fd,sd)
call xyzzyaacy1(js,fd,sd)
xyzzyaaaa6=sum(xyzzyaabk1(:,js))/sum(xyzzyaabk1(:,is))
xyzzyaaaa6=xyzzyaaaa6*exp_protect(xyzzyaabm1(js)-xyzzyaabm1(is))
ratio=cmplx(xyzzyaaaa6,0.d0,dp)
end subroutine wfn_ratio_exmol
subroutine accept_move_exmol(is,js)
implicit none
integer,intent(in) :: is,js
if(xyzzyaabx1(js))then
call dcopy(xyzzyaaaa1,xyzzyaabl1(1,js),1,xyzzyaabl1(1,is),1)
xyzzyaabm1(is)=xyzzyaabm1(js)
xyzzyaabx1(is)=.true.
else
xyzzyaabx1(is)=.false.
endif
if(xyzzyaabw1(js))then
call dcopy(xyzzyaaaa1,xyzzyaabk1(1,js),1,xyzzyaabk1(1,is),1)
xyzzyaabw1(is)=.true.
else
xyzzyaabw1(is)=.false.
endif
if(xyzzyaabs1(js))then
call dcopy(xyzzyaabc1,xyzzyaabg1(1,1,1,js),1,xyzzyaabg1(1,1,1,is),1)
xyzzyaabs1(is)=.true.
else
xyzzyaabs1(is)=.false.
endif
if(xyzzyaabt1(js))then
call dcopy(xyzzyaabd1,xyzzyaabh1(1,1,1,1,js),1,xyzzyaabh1(1,1,1,1,is),&
&1)
xyzzyaabt1(is)=.true.
else
xyzzyaabt1(is)=.false.
endif
if(xyzzyaabu1(js))then
call dcopy(xyzzyaabc1,xyzzyaabi1(1,1,1,js),1,xyzzyaabi1(1,1,1,is),1)
xyzzyaabu1(is)=.true.
else
xyzzyaabu1(is)=.false.
endif
if(xyzzyaabv1(js))then
call dcopy(xyzzyaabe1,xyzzyaabj1(1,1,1,1,js),1,xyzzyaabj1(1,1,1,1,is),&
&1)
xyzzyaabv1(is)=.true.
else
xyzzyaabv1(is)=.false.
endif
if(xyzzyaaca1(js))then
call dcopy(size_farray,xyzzyaabn1(1,1,1,js),1,xyzzyaabn1(1,1,1,is),1)
call dcopy(size_farray*xyzzyaaaa1,xyzzyaabo1(1,1,1,1,js),1,xyzzyaabo1(&
&1,1,1,1,is),1)
xyzzyaaca1(is)=.true.
else
xyzzyaaca1(is)=.false.
endif
if(xyzzyaacb1(js))then
call dcopy(size_harray,xyzzyaabp1(1,1,1,1,1,js),1,xyzzyaabp1(1,1,1,1,1&
&,is),1)
xyzzyaacb1(is)=.true.
else
xyzzyaacb1(is)=.false.
endif
if(use_backflow)call accept_move_bf(is,js)
xyzzyaaby1(:,is)=.false.
xyzzyaabz1(:,is)=.false.
end subroutine accept_move_exmol
subroutine reset_config_exmol(is,js)
implicit none
integer,intent(in) :: is,js
xyzzyaabs1(js)=.false.
xyzzyaabt1(js)=.false.
xyzzyaabu1(js)=.false.
xyzzyaabv1(js)=.false.
xyzzyaabw1(js)=.false.
xyzzyaabx1(js)=.false.
xyzzyaaby1(:,js)=.false.
xyzzyaabz1(:,js)=.false.
xyzzyaaca1(js)=.false.
xyzzyaacb1(js)=.false.
if(use_backflow)call reset_config_bf(is,js)
end subroutine reset_config_exmol
subroutine wfn_logval_exmol(is,logwfn)
implicit none
integer,intent(in) :: is
complex(dp),intent(out) :: logwfn
real(dp) xyzzyaaaa9
call xyzzyaacy1(is,.false.,.false.)
xyzzyaaaa9=sum(xyzzyaabk1(:,is))
if(xyzzyaaaa9>0.d0)then
logwfn=cmplx(log(xyzzyaaaa9)+xyzzyaabm1(is),0.d0,dp)
else
logwfn=cmplx(log(-xyzzyaaaa9)+xyzzyaabm1(is),pi,dp)
endif
end subroutine wfn_logval_exmol
subroutine wfn_loggrad_exmol(ii,is,ilevel,val,sd,loggrad)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: loggrad(3)
logical,intent(in) :: val,sd
if(.not.xyzzyaaby1(ii,is))then
if(.not.use_backflow)then
call xyzzyaacy1(is,.true.,sd)
call xyzzyaacw1(is,val,.true.,sd)
call xyzzyaadj1(ii,xyzzyaabk1(1,is),xyzzyaabh1(1,1,1,1,is),loggrad)
else
if(sd)then
call xyzzyaadb1(is)
else
call xyzzyaada1(is)
endif
call loggrad_bf(ii,is,sd,xyzzyaabn1(1,1,1,is),loggrad)
endif
xyzzyaabq1(1:3,ii,is)=loggrad(1:3)
xyzzyaaby1(ii,is)=.true.
else
loggrad(1:3)=xyzzyaabq1(1:3,ii,is)
endif
end subroutine wfn_loggrad_exmol
subroutine wfn_loglap_exmol(ii,is,ilevel,val,fd,loglap)
implicit none
integer,intent(in) :: ii,is,ilevel
complex(dp),intent(out) :: loglap
logical,intent(in) :: val,fd
complex(dp) xyzzyaaaa11(3)
if(.not.xyzzyaabz1(ii,is))then
if(.not.xyzzyaaby1(ii,is))then
call wfn_loggrad_exmol(ii,is,ilevel,.true.,.true.,xyzzyaaaa11)
else
xyzzyaaaa11=xyzzyaabq1(1:3,ii,is)
call xyzzyaacy1(is,.true.,.true.)
endif
if(.not.use_backflow)then
call xyzzyaacw1(is,val,fd,.true.)
call xyzzyaadk1(ii,xyzzyaabk1(1,is),xyzzyaabh1(1,1,1,1,is),xyzzyaabi1(&
&1,1,1,is),xyzzyaaaa11,loglap)
else
call xyzzyaadb1(is)
call loglap_bf(ii,is,xyzzyaabn1(1,1,1,is),xyzzyaabp1(1,1,1,1,1,is),xyz&
&zyaaaa11,loglap)
endif
xyzzyaabr1(ii,is)=loglap
xyzzyaabz1(ii,is)=.true.
else
loglap=xyzzyaabr1(ii,is)
endif
end subroutine wfn_loglap_exmol
subroutine prefetch_wfn_exmol(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
if(.not.use_backflow)then
call xyzzyaacw1(is,.true.,fd.or.sd,sd)
else
if(sd)then
call xyzzyaadb1(is)
elseif(fd)then
call xyzzyaada1(is)
endif
endif
end subroutine prefetch_wfn_exmol
subroutine xyzzyaacw1(is,val,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13
logical xyzzyaaad13,xyzzyaaae13,xyzzyaaaf13,xyzzyaaag13
xyzzyaaad13=val.and..not.xyzzyaabs1(is)
xyzzyaaae13=fd.and..not.xyzzyaabt1(is)
xyzzyaaaf13=sd.and..not.xyzzyaabu1(is)
if(.not.(xyzzyaaad13.or.xyzzyaaae13.or.xyzzyaaaf13))return
xyzzyaaag13=.false.
xyzzyaaac13=buffer_move1_from(is)
if(xyzzyaaac13/=0)then
xyzzyaaag13=.true.
if(xyzzyaaad13.and..not.xyzzyaabs1(xyzzyaaac13))xyzzyaaag13=.false.
if(xyzzyaaae13.and..not.xyzzyaabt1(xyzzyaaac13))xyzzyaaag13=.false.
if(xyzzyaaaf13.and..not.xyzzyaabu1(xyzzyaaac13))xyzzyaaag13=.false.
endif
if(xyzzyaaag13)then
xyzzyaaaa13=buffer_move1_from_ii(is)
xyzzyaaab13=xyzzyaaao1(xyzzyaaaa13)
call get_eevecs1_ch(xyzzyaaaa13,is)
if(xyzzyaaad13)call dcopy(xyzzyaabc1,xyzzyaabg1(1,1,1,xyzzyaaac13),1,x&
&yzzyaabg1(1,1,1,is),1)
if(xyzzyaaae13)call dcopy(xyzzyaabd1,xyzzyaabh1(1,1,1,1,xyzzyaaac13),1&
&,xyzzyaabh1(1,1,1,1,is),1)
if(xyzzyaaaf13)call dcopy(xyzzyaabc1,xyzzyaabi1(1,1,1,xyzzyaaac13),1,x&
&yzzyaabi1(1,1,1,is),1)
if(xyzzyaaba1(xyzzyaaaa13))then
call xyzzyaado1(xyzzyaaaa13,xyzzyaaat1(xyzzyaaaa13),eevecs1_chscr(1,1,&
&is),xyzzyaaad13,xyzzyaaae13,xyzzyaaaf13,orbvec=xyzzyaabg1(:,xyzzyaaab&
&13,:,is),gradvec=xyzzyaabh1(:,:,xyzzyaaab13,:,is),lapvec=xyzzyaabi1(:&
&,xyzzyaaab13,:,is))
else
call xyzzyaado1(xyzzyaaaa13,xyzzyaaat1(xyzzyaaaa13),eevecs1_chscr(1,1,&
&is),xyzzyaaad13,xyzzyaaae13,xyzzyaaaf13,orbvec=xyzzyaabg1(xyzzyaaab13&
&,:,:,is),gradvec=xyzzyaabh1(:,xyzzyaaab13,:,:,is),lapvec=xyzzyaabi1(x&
&yzzyaaab13,:,:,is))
endif
if(xyzzyaaad13)xyzzyaabs1(is)=.true.
if(xyzzyaaae13)xyzzyaabt1(is)=.true.
if(xyzzyaaaf13)xyzzyaabu1(is)=.true.
else
call get_eevecs(is)
call xyzzyaadn1(eevecs_scr(1,1,1,is),xyzzyaaad13,xyzzyaaae13,xyzzyaaaf&
&13,xyzzyaabg1(1,1,1,is),xyzzyaabh1(1,1,1,1,is),xyzzyaabi1(1,1,1,is))
if(xyzzyaaad13)xyzzyaabs1(is)=.true.
if(xyzzyaaae13)xyzzyaabt1(is)=.true.
if(xyzzyaaaf13)xyzzyaabu1(is)=.true.
endif
end subroutine xyzzyaacw1
subroutine xyzzyaacx1(is,val,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: val,fd,sd
logical xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14
xyzzyaaaa14=val.and..not.xyzzyaabs1(is)
xyzzyaaab14=fd.and..not.xyzzyaabt1(is)
xyzzyaaac14=sd.and..not.xyzzyaabv1(is)
if(.not.(xyzzyaaaa14.or.xyzzyaaab14.or.xyzzyaaac14))return
if(.not.bf_x_valid(is))call get_bf_x(is,.true.,fd.or.sd,sd)
call get_eevecs_bf(is)
call xyzzyaadn1(eevecs_bf_scr(1,1,1,is),xyzzyaaaa14,xyzzyaaab14,xyzzya&
&aac14,orbmat=xyzzyaabg1(1,1,1,is),gradmat=xyzzyaabh1(1,1,1,1,is),d2ma&
&t=xyzzyaabj1(1,1,1,1,is))
if(xyzzyaaaa14)xyzzyaabs1(is)=.true.
if(xyzzyaaab14)xyzzyaabt1(is)=.true.
if(xyzzyaaac14)xyzzyaabv1(is)=.true.
end subroutine xyzzyaacx1
subroutine xyzzyaacy1(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
integer xyzzyaaaa15,xyzzyaaab15
real(dp) xyzzyaaac15
if(xyzzyaabw1(is))return
call xyzzyaacz1(is,fd,sd)
do xyzzyaaaa15=1,xyzzyaaaa1
xyzzyaaab15=xyzzyaaah1(xyzzyaaaa15)
xyzzyaaac15=sign(1.d0,dble(xyzzyaaab15))*xyzzyaaaw1(abs(xyzzyaaab15))
xyzzyaabk1(xyzzyaaaa15,is)=exp_protect(xyzzyaabl1(xyzzyaaaa15,is)-xyzz&
&yaabm1(is))*xyzzyaaac15
enddo
xyzzyaabw1(is)=.true.
end subroutine xyzzyaacy1
subroutine xyzzyaacz1(is,fd,sd)
implicit none
integer,intent(in) :: is
logical,intent(in) :: fd,sd
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16
if(xyzzyaabx1(is))return
if(.not.use_backflow)then
call xyzzyaacw1(is,.true.,fd,sd)
else
call xyzzyaacx1(is,.true.,fd,sd)
endif
xyzzyaabl1(1:xyzzyaaaa1,is)=0.d0
do xyzzyaaab16=1,xyzzyaaad1
do xyzzyaaac16=1,ne
do xyzzyaaaa16=1,xyzzyaaaa1
xyzzyaaad16=xyzzyaaai1(xyzzyaaac16,xyzzyaaab16,xyzzyaaaa16)
if(xyzzyaaad16/=0)xyzzyaabl1(xyzzyaaaa16,is)=xyzzyaabl1(xyzzyaaaa16,is&
&)+xyzzyaabg1(xyzzyaaac16,xyzzyaaab16,xyzzyaaad16,is)
enddo
enddo
enddo
xyzzyaabm1(is)=maxval(xyzzyaabl1(1:xyzzyaaaa1,is))
xyzzyaabx1(is)=.true.
end subroutine xyzzyaacz1
subroutine xyzzyaada1(is)
implicit none
integer,intent(in) :: is
if(xyzzyaaca1(is))return
call xyzzyaacx1(is,.true.,.true.,.false.)
call xyzzyaacy1(is,.true.,.false.)
call xyzzyaadl1(xyzzyaabk1(1,is),xyzzyaabh1(1,1,1,1,is),xyzzyaabn1(1,1&
&,1,is),xyzzyaabo1(1,1,1,1,is))
xyzzyaaca1(is)=.true.
end subroutine xyzzyaada1
subroutine xyzzyaadb1(is)
implicit none
integer,intent(in) :: is
if(xyzzyaaca1(is).and.xyzzyaacb1(is))return
call get_bf_x(is,.true.,.true.,.true.)
call xyzzyaacx1(is,.true.,.true.,.true.)
call xyzzyaacy1(is,.true.,.true.)
if(.not.xyzzyaaca1(is))then
call xyzzyaadl1(xyzzyaabk1(1,is),xyzzyaabh1(1,1,1,1,is),xyzzyaabn1(1,1&
&,1,is),xyzzyaabo1(1,1,1,1,is))
xyzzyaaca1(is)=.true.
endif
call xyzzyaadm1(xyzzyaabk1(1,is),xyzzyaabj1(1,1,1,1,is),bf_m2_scr(1,is&
&),bf_rmap2_scr(1,1,is),xyzzyaabo1(1,1,1,1,is),xyzzyaabp1(1,1,1,1,1,is&
&))
xyzzyaacb1(is)=.true.
end subroutine xyzzyaadb1
subroutine setup_exmol_params(nparam)
implicit none
integer,intent(inout) :: nparam
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,xy&
&zzyaaaf19,xyzzyaaag19,xyzzyaaah19
xyzzyaace1=0
xyzzyaaab19=0
if(opt_det_coeff)then
do xyzzyaaad19=1,xyzzyaaaf1
if(xyzzyaaam1(xyzzyaaad19)==1)xyzzyaaab19=xyzzyaaab19+1
enddo
endif
if(opt_orbitals)then
do xyzzyaaaa19=1,xyzzyaaac1
do xyzzyaaac19=1,xyzzyaaak1(xyzzyaaaa19)
if(xyzzyaaal1(xyzzyaaac19,xyzzyaaaa19)==1)xyzzyaaab19=xyzzyaaab19+1
enddo
enddo
endif
xyzzyaace1(1)=xyzzyaaab19
call setup_bf_params(xyzzyaace1(2))
xyzzyaacd1=sum(xyzzyaace1,xyzzyaace1>0)
nparam=xyzzyaacd1
allocate(xyzzyaacf1(xyzzyaacd1),stat=xyzzyaaae19)
call check_alloc(xyzzyaaae19,'SETUP_EXMOL_PARAMS','exmol_param_sec')
xyzzyaaag19=0
do xyzzyaaah19=1,xyzzyaacc1
if(xyzzyaace1(xyzzyaaah19)<1)cycle
xyzzyaaaf19=xyzzyaaag19+1
xyzzyaaag19=xyzzyaaag19+xyzzyaace1(xyzzyaaah19)
xyzzyaacf1(xyzzyaaaf19:xyzzyaaag19)=xyzzyaaah19
enddo
call xyzzyaadc1
end subroutine setup_exmol_params
subroutine finish_exmol_params
implicit none
call finish_bf_params
call xyzzyaadd1
deallocate(xyzzyaacf1)
end subroutine finish_exmol_params
subroutine xyzzyaadc1
implicit none
integer xyzzyaaaa21,xyzzyaaab21
xyzzyaaab21=xyzzyaace1(1)
if(opt_det_coeff)then
allocate(xyzzyaact1(xyzzyaaaf1,0:xyzzyaaab21),stat=xyzzyaaaa21)
call check_alloc(xyzzyaaaa21,'SETUP_EXMOL_PBUFFER','termcoef')
xyzzyaact1=0.d0
endif
if(opt_orbitals)then
allocate(xyzzyaacu1(maxval(xyzzyaaak1(1:xyzzyaaac1)),xyzzyaaac1,0:xyzz&
&yaaab21),stat=xyzzyaaaa21)
call check_alloc(xyzzyaaaa21,'SETUP_EXMOL_PBUFFER','exmol_param')
xyzzyaacu1=0.d0
xyzzyaacv1=maxval(xyzzyaaak1(1:xyzzyaaac1))*xyzzyaaac1
endif
end subroutine xyzzyaadc1
subroutine xyzzyaadd1
implicit none
if(opt_det_coeff)deallocate(xyzzyaact1)
if(opt_orbitals)deallocate(xyzzyaacu1)
end subroutine xyzzyaadd1
subroutine get_exmol_params(params,has_lolim,lolim,has_hilim,hilim,is_&
&shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,lab&
&el)
implicit none
real(dp),intent(inout) :: params(xyzzyaacd1),lolim(xyzzyaacd1),hilim(x&
&yzzyaacd1)
logical,intent(inout) :: has_lolim(xyzzyaacd1),has_hilim(xyzzyaacd1),i&
&s_shallow(xyzzyaacd1),is_redundant(xyzzyaacd1),is_linear(xyzzyaacd1),&
&is_loglinear(xyzzyaacd1),has_aderiv(xyzzyaacd1),affect_map(xyzzyaacd1&
&,xyzzyaacd1)
character(2),intent(inout) :: label(xyzzyaacd1)
integer xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyzzyaaad23,xyzzyaaae23,xy&
&zzyaaaf23,xyzzyaaag23
real(dp) xyzzyaaah23
affect_map=.false.
label='Ex'
xyzzyaaah23=1.1d-8
xyzzyaaaf23=0
if(xyzzyaace1(1)>0)then
xyzzyaaaa23=xyzzyaaaf23
xyzzyaaae23=xyzzyaaaf23+1
xyzzyaaaf23=xyzzyaaaf23+xyzzyaace1(1)
has_lolim(xyzzyaaae23:xyzzyaaaf23)=.false.
lolim(xyzzyaaae23:xyzzyaaaf23)=0.d0
has_hilim(xyzzyaaae23:xyzzyaaaf23)=.false.
hilim(xyzzyaaae23:xyzzyaaaf23)=0.d0
is_shallow(xyzzyaaae23:xyzzyaaaf23)=.false.
is_redundant(xyzzyaaae23:xyzzyaaaf23)=.false.
is_linear(xyzzyaaae23:xyzzyaaaf23)=.false.
is_loglinear(xyzzyaaae23:xyzzyaaaf23)=.false.
has_aderiv(xyzzyaaae23:xyzzyaaaf23)=.false.
do xyzzyaaag23=1,xyzzyaacd1
affect_map(xyzzyaaag23,xyzzyaaag23)=.true.
enddo
if(opt_det_coeff)then
do xyzzyaaad23=1,xyzzyaaaf1
if(xyzzyaaam1(xyzzyaaad23)==1)then
xyzzyaaaa23=xyzzyaaaa23+1
params(xyzzyaaaa23)=xyzzyaaaw1(xyzzyaaad23)
is_linear(xyzzyaaaa23)=.true.
endif
enddo
xyzzyaaaw1(:)=xyzzyaaaw1(:)
endif
if(opt_orbitals)then
do xyzzyaaab23=1,xyzzyaaac1
do xyzzyaaac23=1,xyzzyaaak1(xyzzyaaab23)
if(xyzzyaaal1(xyzzyaaac23,xyzzyaaab23)==1)then
xyzzyaaaa23=xyzzyaaaa23+1
params(xyzzyaaaa23)=xyzzyaaav1(xyzzyaaac23,xyzzyaaab23)
select case(xyzzyaaaj1(xyzzyaaab23))
case(1)
if(xyzzyaaac23==2)then
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=1.d0
else
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=xyzzyaaah23
endif
case(2)
if(xyzzyaaac23==1)then
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=xyzzyaaah23
elseif(xyzzyaaac23==2)then
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=1.d0
elseif(xyzzyaaac23==4)then
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=xyzzyaaah23
is_shallow(xyzzyaaaa23)=.true.
is_redundant(xyzzyaaaa23)=all(xyzzyaaav1(4:xyzzyaaak1(xyzzyaaab23),xyz&
&zyaaab23)==0.d0)
endif
case(3)
if(xyzzyaaac23==1)then
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=xyzzyaaah23
is_shallow(xyzzyaaaa23)=.true.
is_redundant(xyzzyaaaa23)=all(xyzzyaaav1(2:xyzzyaaak1(xyzzyaaab23),xyz&
&zyaaab23)==0.d0)
endif
case(4)
continue
case(5)
if(xyzzyaaac23==2)then
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=xyzzyaaah23
else
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=xyzzyaaah23
endif
case(6)
has_lolim(xyzzyaaaa23)=.true.
lolim(xyzzyaaaa23)=xyzzyaaah23
end select
endif
enddo
enddo
endif
endif
if(xyzzyaace1(2)>0)then
xyzzyaaae23=xyzzyaaaf23+1
xyzzyaaaf23=xyzzyaaaf23+xyzzyaace1(2)
call get_bf_params(params(xyzzyaaae23:xyzzyaaaf23),has_lolim(xyzzyaaae&
&23:xyzzyaaaf23),lolim(xyzzyaaae23:xyzzyaaaf23),has_hilim(xyzzyaaae23:&
&xyzzyaaaf23),hilim(xyzzyaaae23:xyzzyaaaf23),is_shallow(xyzzyaaae23:xy&
&zzyaaaf23),is_redundant(xyzzyaaae23:xyzzyaaaf23),is_linear(xyzzyaaae2&
&3:xyzzyaaaf23),is_loglinear(xyzzyaaae23:xyzzyaaaf23),has_aderiv(xyzzy&
&aaae23:xyzzyaaaf23),affect_map(xyzzyaaae23:xyzzyaaaf23,xyzzyaaae23:xy&
&zzyaaaf23),label(xyzzyaaae23:xyzzyaaaf23))
endif
end subroutine get_exmol_params
subroutine put_exmol_params(params,ignore,iparam_buffer,prestore,bad_p&
&arams)
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaacd1)
logical,intent(in) :: ignore(xyzzyaacd1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaae24,xy&
&zzyaaaf24,xyzzyaaag24,xyzzyaaah24
real(dp) xyzzyaaai24
logical xyzzyaaaj24
bad_params=.false.
xyzzyaaae24=0
xyzzyaaaf24=0
if(iparam_buffer>0)then
xyzzyaaae24=xyzzyaacf1(iparam_buffer)
xyzzyaaaf24=iparam_buffer-sum(xyzzyaace1(1:xyzzyaaae24-1))
endif
xyzzyaaah24=0
if(xyzzyaace1(1)>0)then
xyzzyaaaa24=xyzzyaaah24
xyzzyaaag24=xyzzyaaah24+1
xyzzyaaah24=xyzzyaaah24+xyzzyaace1(1)
if(xyzzyaaae24==1.or.xyzzyaaae24==0)then
if(prestore)then
call xyzzyaadf1(iparam_buffer)
else
if(opt_det_coeff)then
do xyzzyaaad24=1,xyzzyaaaf1
if(xyzzyaaam1(xyzzyaaad24)==1)then
xyzzyaaaa24=xyzzyaaaa24+1
if(.not.ignore(xyzzyaaaa24))xyzzyaaaw1(xyzzyaaad24)=params(xyzzyaaaa24&
&)
endif
enddo
xyzzyaaai24=sqrt(sum(xyzzyaaaw1(:)**2))
if(xyzzyaaai24==0.d0)call errstop_master('PUT_EXMOL_PARAMS','All term &
&coefficients are zero.')
endif
if(opt_orbitals)then
do xyzzyaaab24=1,xyzzyaaac1
do xyzzyaaac24=1,xyzzyaaak1(xyzzyaaab24)
if(xyzzyaaal1(xyzzyaaac24,xyzzyaaab24)==1)then
xyzzyaaaa24=xyzzyaaaa24+1
if(.not.ignore(xyzzyaaaa24))xyzzyaaav1(xyzzyaaac24,xyzzyaaab24)=params&
&(xyzzyaaaa24)
endif
enddo
select case(xyzzyaaaj1(xyzzyaaab24))
case(1)
if(xyzzyaaaz1)xyzzyaaav1(3,xyzzyaaab24)=-1.d0/(xyzzyaaav1(2,xyzzyaaab2&
&4)*xyzzyaaau1)
case(2)
if(xyzzyaaaz1)xyzzyaaav1(3,xyzzyaaab24)=-1.d0/(xyzzyaaav1(2,xyzzyaaab2&
&4)*xyzzyaaau1)
xyzzyaaav1(6,xyzzyaaab24)=dble(xyzzyaaaq1(1,xyzzyaaab24))*xyzzyaaav1(5&
&,xyzzyaaab24)
case(3)
xyzzyaaav1(3,xyzzyaaab24)=dble(xyzzyaaaq1(1,xyzzyaaab24))*xyzzyaaav1(2&
&,xyzzyaaab24)
if(xyzzyaaaz1)xyzzyaaav1(3,xyzzyaaab24)=xyzzyaaav1(3,xyzzyaaab24)+xyzz&
&yaaav1(1,xyzzyaaab24)*xyzzyaaau1
case(4)
xyzzyaaav1(1,xyzzyaaab24)=-xyzzyaaau1
case(5)
continue
case(6)
continue
end select
enddo
endif
call xyzzyaade1(iparam_buffer)
endif
endif
endif
if(xyzzyaace1(2)>0)then
xyzzyaaag24=xyzzyaaah24+1
xyzzyaaah24=xyzzyaaah24+xyzzyaace1(2)
if(xyzzyaaae24==2.or.xyzzyaaae24==0)then
call put_bf_params(params(xyzzyaaag24:xyzzyaaah24),ignore(xyzzyaaag24:&
&xyzzyaaah24),xyzzyaaaf24,prestore,xyzzyaaaj24)
bad_params=bad_params.or.xyzzyaaaj24
endif
endif
end subroutine put_exmol_params
subroutine xyzzyaade1(indx)
implicit none
integer,intent(in) :: indx
if(opt_det_coeff)then
call dcopy(xyzzyaaaf1,xyzzyaaaw1(1),1,xyzzyaact1(1,indx),1)
endif
if(opt_orbitals)then
call dcopy(xyzzyaacv1,xyzzyaaav1(1,1),1,xyzzyaacu1(1,1,indx),1)
endif
end subroutine xyzzyaade1
subroutine xyzzyaadf1(indx)
implicit none
integer,intent(in) :: indx
if(opt_det_coeff)then
call dcopy(xyzzyaaaf1,xyzzyaact1(1,indx),1,xyzzyaaaw1(1),1)
endif
if(opt_orbitals)then
call dcopy(xyzzyaacv1,xyzzyaacu1(1,1,indx),1,xyzzyaaav1(1,1),1)
endif
end subroutine xyzzyaadf1
subroutine invalidate_param1_exmol(is,iparam)
implicit none
integer,intent(in) :: is,iparam
integer xyzzyaaaa27,xyzzyaaab27
xyzzyaaaa27=xyzzyaacf1(iparam)
xyzzyaaab27=iparam-sum(xyzzyaace1(1:xyzzyaaaa27-1))
xyzzyaabs1(is)=.false.
xyzzyaabt1(is)=.false.
xyzzyaabu1(is)=.false.
xyzzyaabv1(is)=.false.
xyzzyaabw1(is)=.false.
xyzzyaabx1(is)=.false.
xyzzyaaby1(:,is)=.false.
xyzzyaabz1(:,is)=.false.
xyzzyaaca1(is)=.false.
xyzzyaacb1(is)=.false.
if(xyzzyaaaa27==2)call invalidate_param1_bf(is,xyzzyaaab27)
end subroutine invalidate_param1_exmol
subroutine invalidate_params_exmol(ignore)
implicit none
logical,intent(in) :: ignore(xyzzyaacd1)
xyzzyaacr1(:)=.false.
xyzzyaacs1(:)=.false.
xyzzyaacp1(:)=.false.
xyzzyaacq1(:)=.false.
end subroutine invalidate_params_exmol
subroutine clear_scratch_exmol(is)
implicit none
integer,intent(in) :: is
xyzzyaabs1(is)=.false.
xyzzyaabt1(is)=.false.
xyzzyaabu1(is)=.false.
xyzzyaabv1(is)=.false.
xyzzyaabw1(is)=.false.
xyzzyaabx1(is)=.false.
xyzzyaaby1(:,is)=.false.
xyzzyaabz1(:,is)=.false.
xyzzyaaca1(is)=.false.
xyzzyaacb1(is)=.false.
if(use_backflow)call clear_scratch_bf(is)
end subroutine clear_scratch_exmol
subroutine gen_config_exmol(pt_config)
implicit none
integer xyzzyaaaa30
type(config_wfn_exmol),pointer :: pt_config
allocate(pt_config,stat=xyzzyaaaa30)
call check_alloc(xyzzyaaaa30,'GEN_CONFIG_EXMOL','container')
end subroutine gen_config_exmol
subroutine delete_config_exmol(pt_config)
implicit none
type(config_wfn_exmol),pointer :: pt_config
deallocate(pt_config)
end subroutine delete_config_exmol
subroutine copy_config_exmol(pt_from,pt_to)
implicit none
type(config_wfn_exmol),pointer :: pt_from,pt_to
end subroutine copy_config_exmol
subroutine config_to_pt_exmol(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_exmol),pointer :: pt_config
end subroutine config_to_pt_exmol
subroutine pt_to_config_exmol(pt_config)
implicit none
type(config_wfn_exmol),pointer :: pt_config
end subroutine pt_to_config_exmol
subroutine redist_allocations_exmol(kmax)
implicit none
integer,intent(in) :: kmax
end subroutine redist_allocations_exmol
subroutine redist_load_exmol(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_exmol),pointer :: pt_config
end subroutine redist_load_exmol
subroutine redist_send_exmol(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
end subroutine redist_send_exmol
subroutine redist_recv_exmol(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
end subroutine redist_recv_exmol
subroutine redist_save_exmol(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_wfn_exmol),pointer :: pt_config
end subroutine redist_save_exmol
subroutine redist_deallocations_exmol
implicit none
end subroutine redist_deallocations_exmol
subroutine load_from_pt_exmol(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_exmol),pointer :: pt_config
end subroutine load_from_pt_exmol
subroutine save_to_pt_exmol(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_wfn_exmol),pointer :: pt_config
end subroutine save_to_pt_exmol
subroutine clone_scratch_exmol(is,js)
implicit none
integer,intent(in) :: is,js
if(xyzzyaabx1(is).and..not.xyzzyaabx1(js))then
call dcopy(xyzzyaaaa1,xyzzyaabl1(1,is),1,xyzzyaabl1(1,js),1)
xyzzyaabm1(js)=xyzzyaabm1(is)
xyzzyaabx1(js)=.true.
endif
if(xyzzyaabw1(is).and..not.xyzzyaabw1(js))then
call dcopy(xyzzyaaaa1,xyzzyaabk1(1,is),1,xyzzyaabk1(1,js),1)
xyzzyaabw1(js)=.true.
endif
if(xyzzyaabs1(is).and..not.xyzzyaabs1(js))then
call dcopy(xyzzyaabc1,xyzzyaabg1(1,1,1,is),1,xyzzyaabg1(1,1,1,js),1)
xyzzyaabs1(js)=.true.
endif
if(xyzzyaabt1(is).and..not.xyzzyaabt1(js))then
call dcopy(xyzzyaabd1,xyzzyaabh1(1,1,1,1,is),1,xyzzyaabh1(1,1,1,1,js),&
&1)
xyzzyaabt1(js)=.true.
endif
if(xyzzyaabu1(is).and..not.xyzzyaabu1(js))then
call dcopy(xyzzyaabc1,xyzzyaabi1(1,1,1,is),1,xyzzyaabi1(1,1,1,js),1)
xyzzyaabu1(js)=.true.
endif
if(xyzzyaabv1(is).and..not.xyzzyaabv1(js))then
call dcopy(xyzzyaabe1,xyzzyaabj1(1,1,1,1,is),1,xyzzyaabj1(1,1,1,1,js),&
&1)
xyzzyaabv1(js)=.true.
endif
if(all(xyzzyaaby1(:,is)).and..not.all(xyzzyaaby1(:,js)))then
xyzzyaabq1(:,:,js)=xyzzyaabq1(:,:,is)
xyzzyaaby1(:,js)=.true.
endif
if(all(xyzzyaabz1(:,is)).and..not.all(xyzzyaabz1(:,js)))then
xyzzyaabr1(:,js)=xyzzyaabr1(:,is)
xyzzyaabz1(:,js)=.true.
endif
if(xyzzyaaca1(is).and..not.xyzzyaaca1(js))then
call dcopy(size_farray,xyzzyaabn1(1,1,1,is),1,xyzzyaabn1(1,1,1,js),1)
call dcopy(size_farray*xyzzyaaaa1,xyzzyaabo1(1,1,1,1,is),1,xyzzyaabo1(&
&1,1,1,1,js),1)
xyzzyaaca1(js)=.true.
endif
if(xyzzyaacb1(is).and..not.xyzzyaacb1(js))then
call dcopy(size_harray,xyzzyaabp1(1,1,1,1,1,is),1,xyzzyaabp1(1,1,1,1,1&
&,js),1)
xyzzyaacb1(js)=.true.
endif
end subroutine clone_scratch_exmol
subroutine add_config_exmol_items(is)
use slaarnaaf,only : add_config
implicit none
integer,intent(in) :: is
integer xyzzyaaaa44,xyzzyaaab44
real(dp) fidet(3,netot,real1_complex2),lapdet
complex(dp) xyzzyaaac44(3),xyzzyaaad44
fidet=0.d0
lapdet=0.d0
do xyzzyaaaa44=1,netot
xyzzyaaab44=which_spin(xyzzyaaaa44)
call wfn_loggrad_exmol(xyzzyaaaa44,is,0,.false.,.true.,xyzzyaaac44)
call wfn_loglap_exmol(xyzzyaaaa44,is,0,.false.,.false.,xyzzyaaad44)
fidet(:,xyzzyaaaa44,1)=dble(xyzzyaaac44(:))
if(complex_wf)fidet(:,xyzzyaaaa44,2)=aimag(xyzzyaaac44(:))
lapdet=lapdet+inv_pmass(xyzzyaaab44)*dble(xyzzyaaad44)
enddo
call add_config(modify=.true.,fidet=fidet,lapdet=lapdet)
end subroutine add_config_exmol_items
subroutine setup_storage_exmol(nconfig,ignore)
implicit none
integer,intent(in) :: nconfig
logical,intent(in) :: ignore(xyzzyaacd1)
integer xyzzyaaaa45
allocate(xyzzyaacn1(3,xyzzyaaab1,nconfig),xyzzyaaco1(xyzzyaaab1,nconfi&
&g),xyzzyaack1(xyzzyaaaa1,nconfig),xyzzyaacl1(xyzzyaaaa1,nconfig),xyzz&
&yaacm1(nconfig),stat=xyzzyaaaa45)
call check_alloc(xyzzyaaaa45,'SETUP_STORAGE_SLATER','*_store')
xyzzyaacn1=czero
xyzzyaaco1=czero
xyzzyaack1=0.d0
xyzzyaacl1=0.d0
xyzzyaacm1=0.d0
allocate(xyzzyaacp1(nconfig),xyzzyaacq1(nconfig),xyzzyaacr1(nconfig),x&
&yzzyaacs1(nconfig),stat=xyzzyaaaa45)
call check_alloc(xyzzyaaaa45,'SETUP_STORAGE_SLATER','*_svalid')
xyzzyaacp1=.false.
xyzzyaacq1=.false.
xyzzyaacr1=.false.
xyzzyaacs1=.false.
end subroutine setup_storage_exmol
subroutine finish_storage_exmol
implicit none
deallocate(xyzzyaacn1,xyzzyaaco1,xyzzyaack1,xyzzyaacl1,xyzzyaacm1)
deallocate(xyzzyaacp1,xyzzyaacq1,xyzzyaacr1,xyzzyaacs1)
end subroutine finish_storage_exmol
subroutine load_from_storage_exmol(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(xyzzyaacp1(icfg))then
xyzzyaabq1(1:3,1:xyzzyaaab1,is)=xyzzyaacn1(1:3,1:xyzzyaaab1,icfg)
xyzzyaaby1(1:xyzzyaaab1,is)=.true.
endif
if(xyzzyaacq1(icfg))then
xyzzyaabr1(1:xyzzyaaab1,is)=xyzzyaaco1(1:xyzzyaaab1,icfg)
xyzzyaabz1(1:xyzzyaaab1,is)=.true.
endif
if(xyzzyaacr1(icfg))then
call dcopy(xyzzyaaaa1,xyzzyaack1(1,icfg),1,xyzzyaabk1(1,is),1)
xyzzyaabw1(is)=.true.
endif
if(xyzzyaacs1(icfg))then
call dcopy(xyzzyaaaa1,xyzzyaacl1(1,icfg),1,xyzzyaabl1(1,is),1)
xyzzyaabm1(is)=xyzzyaacm1(icfg)
xyzzyaabx1(is)=.true.
endif
end subroutine load_from_storage_exmol
subroutine save_to_storage_exmol(is,icfg)
implicit none
integer,intent(in) :: is,icfg
if(all(xyzzyaaby1(:,is)))then
xyzzyaacn1(1:3,1:xyzzyaaab1,icfg)=xyzzyaabq1(1:3,1:xyzzyaaab1,is)
xyzzyaacp1(icfg)=.true.
endif
if(all(xyzzyaabz1(:,is)))then
xyzzyaaco1(1:xyzzyaaab1,icfg)=xyzzyaabr1(1:xyzzyaaab1,is)
xyzzyaacq1(is)=.true.
endif
if(xyzzyaabw1(is))then
call dcopy(xyzzyaaaa1,xyzzyaabk1(1,is),1,xyzzyaack1(1,icfg),1)
xyzzyaacr1(icfg)=.true.
endif
if(xyzzyaabx1(is))then
call dcopy(xyzzyaaaa1,xyzzyaabl1(1,is),1,xyzzyaacl1(1,icfg),1)
xyzzyaacm1(icfg)=xyzzyaabm1(is)
xyzzyaacs1(icfg)=.true.
endif
end subroutine save_to_storage_exmol
subroutine enumerate_plot_exmol(n,keyword,description)
implicit none
integer,intent(inout) :: n
character(64),intent(inout),optional :: keyword(:),description(:)
n=0
end subroutine enumerate_plot_exmol
subroutine query_plot_exmol(iplot,ii,rank,is_complex,has_stderr,rot_te&
&nsor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
nfunctions=0
end subroutine query_plot_exmol
subroutine get_plot_exmol(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
f(1)=0.d0
end subroutine get_plot_exmol
subroutine finish_plot_exmol
implicit none
end subroutine finish_plot_exmol
subroutine read_exmol
implicit none
integer xyzzyaaaa53,xyzzyaaab53,xyzzyaaac53,set,xyzzyaaad53,xyzzyaaae5&
&3,xyzzyaaaf53,xyzzyaaag53,xyzzyaaah53(2),xyzzyaaai53,xyzzyaaaj53,i1,i&
&2,xyzzyaaak53,xyzzyaaal53,xyzzyaaam53,xyzzyaaan53,xyzzyaaao53,xyzzyaa&
&ap53,ispin2_first
integer,allocatable :: xyzzyaaaq53(:,:,:,:),xyzzyaaar53(:,:),xyzzyaaas&
&53(:,:)
real(dp) xyzzyaaat53,xyzzyaaau53
logical xyzzyaaav53,xyzzyaaaw53
logical,parameter :: xyzzyaaax53=.true.
character(80) char_80
character(80),allocatable :: xyzzyaaay53(:)
type perm
integer,pointer :: matrix(:,:)
integer csign
type(perm),pointer :: next,next_term
end type perm
type(perm),pointer :: xyzzyaaaz53
call open_units(xyzzyaabf1,xyzzyaaab53)
if(xyzzyaaab53/=0)call errstop('READ_EXMOL','Unable to find free i/o u&
&nit.')
open(unit=xyzzyaabf1,file='correlation.data',status='old',iostat=xyzzy&
&aaab53)
if(xyzzyaaab53/=0)call errstop('READ_EXMOL','Problem opening correlati&
&on.data')
if(am_master)then
call wout('Excitonic-molecule wave function')
call wout('================================')
call wout('Reading correlation.data file.')
call wout()
endif
do
read(xyzzyaabf1,'(a)',iostat=xyzzyaaab53)char_80
if(trim(adjustl(char_80))=='START EXMOL')exit
if(xyzzyaaab53>0)call errstop_master('READ_EXMOL','Problem reading cor&
&relation.data.')
if(xyzzyaaab53<0)call errstop_master('READ_EXMOL','Could not find "STA&
&RT EXMOL" in correlation.data.')
enddo
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,'(a)',err=101,end=101)title
if(all(pfermion))then
ispin2_first=0
do xyzzyaaai53=2,nspin
if(sign(1.d0,pcharge(xyzzyaaai53))==sign(1.d0,pcharge(1)))then
if(ispin2_first>0)call errstop_master('READ_EXMOL','This wave function&
& requires that all positively- and negatively-charged particles be de&
&fined consecutively in the PARTICLES block - they are interleaved in &
&your input.')
else
if(ispin2_first==0)ispin2_first=xyzzyaaai53
endif
enddo
if(ispin2_first==0)call errstop_master('READ_EXMOL','This wave functio&
&n requires the presence of positively- and negatively-charged particl&
&es.')
xyzzyaaay1=.false.
if(mod(nspin,2)==0.and.ispin2_first==nspin/2+1)then
xyzzyaaay1=.false.
do xyzzyaaai53=1,ispin2_first-1
if(which_eqvfam(xyzzyaaai53)/=which_eqvfam(xyzzyaaai53+nspin/2))then
xyzzyaaay1=.false.
exit
endif
enddo
endif
xyzzyaaag1(1)=ispin2_first-1
xyzzyaaag1(2)=nspin-ispin2_first+1
allocate(xyzzyaaap1(maxval(xyzzyaaag1),2),stat=xyzzyaaaa53)
call check_alloc(xyzzyaaaa53,'READ_EXMOL','nele_exmol')
xyzzyaaap1=0
xyzzyaaap1(1:xyzzyaaag1(1),1)=nele(1:xyzzyaaag1(1))
xyzzyaaap1(1:xyzzyaaag1(2),2)=nele(xyzzyaaag1(1)+1:nspin)
ne=sum(xyzzyaaap1(:,1))
xyzzyaaad1=sum(xyzzyaaap1(:,2))
xyzzyaaae1=max(ne,xyzzyaaad1)
xyzzyaaab1=ne+xyzzyaaad1
if(ne==0.or.xyzzyaaad1==0)call errstop_master('READ_EXMOL','One of the&
& particle families has no particles in it.')
allocate(xyzzyaaan1(xyzzyaaae1,2),xyzzyaaao1(xyzzyaaab1),xyzzyaaba1(xy&
&zzyaaab1),xyzzyaaax1(xyzzyaaab1),xyzzyaaas1(xyzzyaaab1),xyzzyaaat1(xy&
&zzyaaab1),stat=xyzzyaaaa53)
call check_alloc(xyzzyaaaa53,'READ_EXMOL','which_part')
xyzzyaaan1=0
xyzzyaaao1=0
xyzzyaaba1=.false.
xyzzyaaax1=1.d0
do xyzzyaaae53=1,ne
xyzzyaaan1(xyzzyaaae53,1)=xyzzyaaae53
xyzzyaaao1(xyzzyaaae53)=xyzzyaaae53
xyzzyaaba1(xyzzyaaae53)=.false.
xyzzyaaax1(xyzzyaaae53)=1.d0
xyzzyaaas1(xyzzyaaae53)=2
xyzzyaaat1(xyzzyaaae53)=xyzzyaaad1
enddo
do xyzzyaaaf53=1,xyzzyaaad1
xyzzyaaan1(xyzzyaaaf53,2)=xyzzyaaaf53+ne
xyzzyaaao1(xyzzyaaaf53+ne)=xyzzyaaaf53
xyzzyaaba1(xyzzyaaaf53+ne)=.true.
xyzzyaaax1(xyzzyaaaf53+ne)=-1.d0
xyzzyaaas1(xyzzyaaaf53+ne)=1
xyzzyaaat1(xyzzyaaaf53+ne)=ne
enddo
elseif(all(.not.pfermion))then
if(no_families/=1)call errstop_master('READ_EXMOL','This wave function&
& can only be applied to systems with one bosonic ''particle family'' &
&(with its internal spin division)')
if(nspin/=2)call errstop_master('READ_EXMOL','This wave function can o&
&nly be applied to systems with two bosonic particle types.')
if(netot/=2.or.nele(1)/=nele(2))call errstop_master('READ_EXMOL','This&
& wave function can only be applied to a system of two bosons.')
ne=2
xyzzyaaad1=2
xyzzyaaae1=2
xyzzyaaab1=4
else
call errstop_master('READ_EXMOL','Bosons and fermions mixed?')
endif
call xyzzyaaba53
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,iostat=xyzzyaaab53)xyzzyaaaz1
if(xyzzyaaab53/=0)call errstop_master('READ_EXMOL','Problem reading cu&
&sp flag from correlation.data.')
xyzzyaaau1=0.d0
if(xyzzyaaaz1)then
xyzzyaaai53=which_spin(xyzzyaaan1(1,1))
xyzzyaaaj53=which_spin(xyzzyaaan1(1,2))
if(any(heg_layer(1:xyzzyaaag1(1))/=heg_layer(xyzzyaaai53)).or.any(heg_&
&layer(xyzzyaaag1(1)+1:nspin)/=heg_layer(xyzzyaaaj53)))call errstop_ma&
&ster('READ_EXMOL','Current limitation: particles in same family canno&
&t be in different layers when the cusp is imposed on the EXMOL orbita&
&ls.')
if(heg_layer(xyzzyaaai53)==heg_layer(xyzzyaaaj53))xyzzyaaau1=ee_kato_g&
&amma(xyzzyaaai53,xyzzyaaaj53)
do xyzzyaaai53=1,nspin
do xyzzyaaaj53=1,nspin
ee_cusp_in_orbital(xyzzyaaai53,xyzzyaaaj53)=xyzzyaaai53<=xyzzyaaag1(1)&
&.neqv.xyzzyaaaj53<=xyzzyaaag1(1)
enddo
enddo
endif
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,iostat=xyzzyaaab53)xyzzyaaaf1
if(xyzzyaaab53/=0)call errstop_master('READ_EXMOL','Problem reading nu&
&mber of terms from correlation.data.')
if(xyzzyaaaf1<1)call errstop_master('READ_EXMOL','Number of terms shou&
&ld be greater than zero.')
allocate(xyzzyaaaw1(xyzzyaaaf1),xyzzyaaam1(xyzzyaaaf1),stat=xyzzyaaaa5&
&3)
call check_alloc(xyzzyaaaa53,'READ_EXMOL','termcoef')
xyzzyaaaw1=1.d0
xyzzyaaam1=1
if(am_master)then
call wout('Information:')
call wout(' Title              :  '//trim(adjustl(title)))
if(ne==xyzzyaaad1)then
call wout(' System             :  Ps'//trim(i2s(ne)))
else
call wout(' System             :  Ps'//trim(i2s(min(ne,xyzzyaaad1)))//&
&repeat('+',max(0,xyzzyaaad1-ne))//repeat('-',max(0,ne-xyzzyaaad1)))
endif
call wout(' Factors per term   :  '//trim(i2s(ne*xyzzyaaad1)))
call wout(' Indep terms in sum :  '//trim(i2s(xyzzyaaaf1)))
endif
allocate(xyzzyaaas53(ne,xyzzyaaad1),stat=xyzzyaaaa53)
call check_alloc(xyzzyaaaa53,'READ_EXMOL','temp_matrix')
xyzzyaaas53=0
do xyzzyaaag53=1,xyzzyaaaf1
do
read(xyzzyaabf1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='START TERM '//trim(i2s(xyzzyaaag53)))exit
if(trim(adjustl(char_80))/='')call errstop_master('READ_EXMOL','Was ex&
&pecting to find "START TERM '//trim(i2s(xyzzyaaag53))//'".')
enddo
if(am_master)call wout(' TERM '//trim(i2s(xyzzyaaag53))//':')
read(xyzzyaabf1,*,err=100,end=101)
do xyzzyaaae53=1,ne
read(xyzzyaabf1,*,iostat=xyzzyaaab53)xyzzyaaas53(xyzzyaaae53,1:xyzzyaa&
&ad1)
enddo
if(xyzzyaaab53/=0)call errstop_master('READ_EXMOL','Problem reading th&
&e structure of term '//trim(i2s(xyzzyaaag53))//'.')
if(any(xyzzyaaas53<0))call errstop_master('READ_EXMOL','A matrix eleme&
&nt in term '//trim(i2s(xyzzyaaag53))//' is negative.')
if(all(xyzzyaaas53==0))call errstop_master('READ_EXMOL','All matrix el&
&ements in term '//trim(i2s(xyzzyaaag53))//' are zero.')
call xyzzyaabb53(xyzzyaaag53,xyzzyaaas53,xyzzyaaac53)
xyzzyaaaa1=xyzzyaaaa1+xyzzyaaac53
if(am_master)then
i1=minval(xyzzyaaas53)
i2=maxval(xyzzyaaas53)
if(i1==i2)then
call wout('  Uses function     :  '//trim(i2s(i1)))
else
call wout('  Uses functions    :  '//trim(i2s(i1))//'-'//trim(i2s(i2))&
&)
endif
call wout('  Terms gen by symm :  '//trim(i2s(xyzzyaaac53)))
endif
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,iostat=xyzzyaaab53)xyzzyaaat53,xyzzyaaac53
if(xyzzyaaab53/=0)then
backspace xyzzyaabf1
xyzzyaaat53=1.d0
xyzzyaaac53=1
endif
if(xyzzyaaac53/=0.and.xyzzyaaac53/=1)call errstop_master('READ_EXMOL',&
&'Optimizable flag should be 0 or 1 .')
xyzzyaaaw1(xyzzyaaag53)=xyzzyaaat53
xyzzyaaam1(xyzzyaaag53)=xyzzyaaac53
if(am_master)then
char_80=r2s2(xyzzyaaaw1(xyzzyaaag53),'(f21.12)')
if(xyzzyaaam1(xyzzyaaag53)==1)then
call wout('  c_'//trim(i2s(xyzzyaaag53))//'         (opt) :  '//trim(c&
&har_80))
else
call wout('  c_'//trim(i2s(xyzzyaaag53))//'       (fixed) :  '//trim(c&
&har_80))
endif
endif
do
read(xyzzyaabf1,'(a)',iostat=xyzzyaaab53)char_80
if(xyzzyaaab53/=0)call errstop_master('READ_EXMOL','Can''t find "END T&
&ERM '//trim(i2s(xyzzyaaag53))//'".')
if(trim(adjustl(char_80))=='END TERM '//trim(i2s(xyzzyaaag53)))exit
enddo
enddo
xyzzyaaau53=sqrt(sum(xyzzyaaaw1(:)**2))
if(xyzzyaaau53==0.d0)call errstop_master('READ_EXMOL','All term coeffi&
&cients are zero.')
allocate(xyzzyaaah1(xyzzyaaaa1),xyzzyaaai1(ne,xyzzyaaad1,xyzzyaaaa1),s&
&tat=xyzzyaaaa53)
call check_alloc(xyzzyaaaa53,'READ_EXMOL','fchooser')
call xyzzyaabc53(xyzzyaaaf1,xyzzyaaaa1,xyzzyaaai1,xyzzyaaah1)
xyzzyaaac1=maxval(xyzzyaaai1)
do set=1,xyzzyaaac1
if(.not.any(xyzzyaaai1==set))call errstop_master('READ_EXMOL','There i&
&s no matrix element with value '//trim(i2s(set))//', but the overall &
&range appears to be 1-'//trim(i2s(xyzzyaaac1))//'. There should be no&
& gaps.')
enddo
if(am_master)then
if(xyzzyaaax53)then
call wout(' Full wave function structure')
allocate(xyzzyaaay53(ne),stat=xyzzyaaaa53)
call check_alloc(xyzzyaaaa53,'READ_EXMOL','charvec')
xyzzyaaap53=len_trim(i2s(xyzzyaaac1))
xyzzyaaam53=(1+xyzzyaaap53)*xyzzyaaad1+2
xyzzyaaan53=ne+1-(ne+1)/2
i2=0
do xyzzyaaao53=1,xyzzyaaaf1
call wout('  Term '//trim(i2s(xyzzyaaao53))//' expands to')
xyzzyaaay53=''
i1=i2+1
i2=i2+count(abs(xyzzyaaah1)==xyzzyaaao53)
do xyzzyaaag53=i1,i2
xyzzyaaaw53=xyzzyaaah1(xyzzyaaag53)<0
do xyzzyaaae53=1,ne
char_80=''
do xyzzyaaaf53=1,xyzzyaaad1
xyzzyaaac53=xyzzyaaai1(xyzzyaaae53,xyzzyaaaf53,xyzzyaaag53)
char_80=trim(char_80)//repeat(' ',1+xyzzyaaap53-len_trim(i2s(xyzzyaaac&
&53)))//trim(i2s(xyzzyaaac53))
enddo
if(xyzzyaaan53/=xyzzyaaae53)then
xyzzyaaay53(xyzzyaaae53)=trim(xyzzyaaay53(xyzzyaaae53))//'  '//trim(ch&
&ar_80)
else
if(xyzzyaaaw53)then
xyzzyaaay53(xyzzyaaae53)=trim(xyzzyaaay53(xyzzyaaae53))//' -'//trim(ch&
&ar_80)
else
xyzzyaaay53(xyzzyaaae53)=trim(xyzzyaaay53(xyzzyaaae53))//' +'//trim(ch&
&ar_80)
endif
endif
enddo
if(len_trim(xyzzyaaay53(1))+xyzzyaaam53>76)then
do xyzzyaaae53=1,ne
call wout('  '//trim(xyzzyaaay53(xyzzyaaae53)))
enddo
call wout()
xyzzyaaay53=''
endif
enddo
if(len_trim(xyzzyaaay53(1))>0)then
do xyzzyaaae53=1,ne
call wout('  '//trim(xyzzyaaay53(xyzzyaaae53)))
enddo
call wout()
xyzzyaaay53=''
endif
enddo
deallocate(xyzzyaaay53)
endif
call wout(' Total terms in sum :  '//trim(i2s(xyzzyaaaa1)))
call wout(' Num. of functions  :  '//trim(i2s(xyzzyaaac1)))
endif
allocate(xyzzyaabb1(xyzzyaaac1),xyzzyaaaj1(xyzzyaaac1),xyzzyaaak1(xyzz&
&yaaac1),xyzzyaaar1(xyzzyaaac1),stat=xyzzyaaaa53)
xyzzyaabb1=''
xyzzyaaaj1=0
xyzzyaaak1=0
xyzzyaaar1=0
do set=1,xyzzyaaac1
do
read(xyzzyaabf1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='START FUNCTION '//trim(i2s(set)))exit
if(set==1)then
if(trim(adjustl(char_80))/='')call errstop_master('READ_EXMOL','Was ex&
&pecting to find "START FUNCTION '//trim(i2s(set))//'".')
endif
enddo
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,'(a)',err=101,end=101)xyzzyaabb1(set)
select case(trim(adjustl(xyzzyaabb1(set))))
case('exp')
xyzzyaaaj1(set)=1
xyzzyaaar1(set)=0
xyzzyaaak1(set)=3
case('exp*exppoly')
xyzzyaaaj1(set)=2
xyzzyaaar1(set)=2
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,err=100,end=101)xyzzyaaac53
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,err=100,end=101)xyzzyaaac53
xyzzyaaak1(set)=5+xyzzyaaac53
case('exppoly')
xyzzyaaaj1(set)=3
xyzzyaaar1(set)=2
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,err=100,end=101)xyzzyaaac53
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,err=100,end=101)xyzzyaaac53
xyzzyaaak1(set)=2+xyzzyaaac53
case('expproper')
if(.not.xyzzyaaaz1)call errstop_master('READ_EXMOL','In order to use t&
&he expproper pairing function, "Cusp in orbitals" needs to be set to &
&"T".')
xyzzyaaaj1(set)=4
xyzzyaaar1(set)=0
xyzzyaaak1(set)=1
case('pexp')
xyzzyaaaj1(set)=5
xyzzyaaar1(set)=0
xyzzyaaak1(set)=2
case('expproper_param')
if(.not.xyzzyaaaz1)call errstop_master('READ_EXMOL','In order to use t&
&he expproper_param pairing function, "Cusp in orbitals" needs to be s&
&et to "T".')
xyzzyaaaj1(set)=6
xyzzyaaar1(set)=0
xyzzyaaak1(set)=1
case default
call errstop_master('READ_EXMOL','Unknown function type '''//trim(adju&
&stl(xyzzyaabb1(set)))//''' in set '//trim(i2s(set)))
end select
do
read(xyzzyaabf1,'(a)',iostat=xyzzyaaab53)char_80
if(xyzzyaaab53/=0)call errstop_master('READ_EXMOL','Can''t find "END F&
&UNCTION '//trim(i2s(set))//'".')
if(trim(adjustl(char_80))=='END FUNCTION '//trim(i2s(set)))exit
enddo
enddo
xyzzyaaak53=maxval(xyzzyaaak1(1:xyzzyaaac1))
xyzzyaaal53=maxval(xyzzyaaar1(1:xyzzyaaac1))
allocate(xyzzyaaav1(xyzzyaaak53,xyzzyaaac1),xyzzyaaal1(xyzzyaaak53,xyz&
&zyaaac1),xyzzyaaaq1(xyzzyaaal53,xyzzyaaac1),stat=xyzzyaaaa53)
xyzzyaaav1(:,:)=0.d0
xyzzyaaal1(:,:)=1
rewind xyzzyaabf1
do
read(xyzzyaabf1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='START EXMOL')exit
enddo
do set=1,xyzzyaaac1
do
read(xyzzyaabf1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='START FUNCTION '//trim(i2s(set)))exit
if(set/=1)then
if(trim(adjustl(char_80))/='')call errstop_master('READ_EXMOL','Was ex&
&pecting to find "START FUNCTION '//trim(i2s(set))//'".')
endif
enddo
if(am_master)call wout(' FUNCTION '//trim(i2s(set))//':')
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,'(a)',err=101,end=101)char_80
if(am_master)then
call wout('  Function type     :  '//trim(adjustl(xyzzyaabb1(set))))
call wout('  Number of params  :  '//trim(i2s(xyzzyaaak1(set))))
endif
select case(xyzzyaaaj1(set))
case(1)
xyzzyaaav1(1,set)=2.d0
xyzzyaaav1(2,set)=1.5d0
xyzzyaaal1(3,set)=-1
case(2)
xyzzyaaav1(1,set)=2.d0
xyzzyaaav1(2,set)=1.5d0
xyzzyaaav1(4,set)=10.d0
xyzzyaaal1(6,set)=-1
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,err=100,end=101)xyzzyaaaq1(1,set)
if(xyzzyaaaq1(1,set)<0)call errstop_master('READ_EXMOL','Truncation or&
&der of polynomial cannot be negative.')
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,err=100,end=101)xyzzyaaaq1(2,set)
if(xyzzyaaaq1(2,set)<0)call errstop_master('READ_EXMOL','Expansion ord&
&er of polynomial cannot be negative.')
case(3)
xyzzyaaav1(1,set)=10.d0
xyzzyaaal1(3,set)=-1
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,err=100,end=101)xyzzyaaaq1(1,set)
if(xyzzyaaaq1(1,set)<0)call errstop_master('READ_EXMOL','Truncation or&
&der of polynomial cannot be negative.')
read(xyzzyaabf1,*,err=100,end=101)
read(xyzzyaabf1,*,err=100,end=101)xyzzyaaaq1(2,set)
if(xyzzyaaaq1(2,set)<0)call errstop_master('READ_EXMOL','Expansion ord&
&er of polynomial cannot be negative.')
case(4)
xyzzyaaav1(1,set)=-xyzzyaaau1
xyzzyaaal1(1,set)=-1
case(5)
xyzzyaaav1(1,set)=2.d0
xyzzyaaav1(2,set)=4.d0
case(6)
xyzzyaaav1(1,set)=-xyzzyaaau1
end select
if(any(xyzzyaaal1(1:xyzzyaaak1(set),set)>=0))then
read(xyzzyaabf1,*,err=100,end=101)
do xyzzyaaad53=1,xyzzyaaak1(set)
if(xyzzyaaal1(xyzzyaaad53,set)<0)cycle
read(xyzzyaabf1,*,iostat=xyzzyaaab53)xyzzyaaat53,xyzzyaaac53
if(xyzzyaaab53/=0)then
backspace xyzzyaabf1
if(am_master)call wout('  Unspecified parms :  Set to default')
exit
endif
if(xyzzyaaac53/=0.and.xyzzyaaac53/=1)call errstop_master('READ_EXMOL',&
&'Optimizable flag should be 0 or 1 .')
xyzzyaaav1(xyzzyaaad53,set)=xyzzyaaat53
xyzzyaaal1(xyzzyaaad53,set)=xyzzyaaac53
if(am_master)then
char_80=r2s2(xyzzyaaav1(xyzzyaaad53,set),'(f21.12)')
if(xyzzyaaal1(xyzzyaaad53,set)==1)then
call wout('  p_'//trim(i2s(xyzzyaaad53))//','//trim(i2s(set))//'      &
& (opt) :  '//trim(char_80))
else
call wout('  p_'//trim(i2s(xyzzyaaad53))//','//trim(i2s(set))//'     (&
&fixed) :  '//trim(char_80))
endif
endif
enddo
endif
select case(xyzzyaaaj1(set))
case(1)
if(xyzzyaaaz1)xyzzyaaav1(3,set)=-1.d0/(xyzzyaaav1(2,set)*xyzzyaaau1)
case(2)
if(xyzzyaaaz1)xyzzyaaav1(3,set)=-1.d0/(xyzzyaaav1(2,set)*xyzzyaaau1)
xyzzyaaav1(6,set)=dble(xyzzyaaaq1(1,set))*xyzzyaaav1(5,set)
case(3)
xyzzyaaav1(3,set)=dble(xyzzyaaaq1(1,set))*xyzzyaaav1(2,set)
if(xyzzyaaaz1)xyzzyaaav1(3,set)=xyzzyaaav1(3,set)+xyzzyaaav1(1,set)*xy&
&zzyaaau1
case(4)
xyzzyaaav1(1,set)=-xyzzyaaau1
case(5)
continue
case(6)
continue
end select
read(xyzzyaabf1,'(a)',iostat=xyzzyaaab53)char_80
if(xyzzyaaab53/=0)call errstop_master('READ_EXMOL','Can''t find "END F&
&UNCTION '//trim(i2s(set))//'".')
if(trim(adjustl(char_80))/='END FUNCTION '//trim(i2s(set)))call errsto&
&p_master('READ_EXMOL','Can''t find "END FUNCTION '//trim(i2s(set)) //&
&'".')
enddo
do
read(xyzzyaabf1,'(a)',err=100,end=101)char_80
if(trim(adjustl(char_80))=='END EXMOL')exit
if(trim(adjustl(char_80))/='')call errstop_master('READ_EXMOL','Was ex&
&pecting to find "END EXMOL".')
enddo
close(xyzzyaabf1)
if(am_master)call wout()
xyzzyaabc1=ne*xyzzyaaad1*xyzzyaaac1
xyzzyaabd1=3*xyzzyaabc1
xyzzyaabe1=6*xyzzyaabc1
return
100 call errstop_master('INIT_EXMOL','Error reading correlation.data .&
&')
101 call errstop_master('INIT_EXMOL','File correlation.data ended unex&
&pectedly.')
contains
subroutine xyzzyaaba53
implicit none
integer xyzzyaaaa54,xyzzyaaab54,xyzzyaaac54,xyzzyaaad54,xyzzyaaae54,xy&
&zzyaaaf54,xyzzyaaag54,xyzzyaaah54,xyzzyaaai54,xyzzyaaaj54,xyzzyaaak54&
&,xyzzyaaal54,xyzzyaaam54
nullify(xyzzyaaaz53)
xyzzyaaah53=0
do xyzzyaaaa54=1,2
xyzzyaaam54=0
if(xyzzyaaaa54==2)xyzzyaaam54=xyzzyaaag1(1)
do xyzzyaaab54=1,xyzzyaaag1(xyzzyaaaa54)
xyzzyaaae54=xyzzyaaap1(xyzzyaaab54,xyzzyaaaa54)
xyzzyaaah53(xyzzyaaaa54)=xyzzyaaah53(xyzzyaaaa54)+(xyzzyaaae54*(xyzzya&
&aae54-1))/2
if(.not.any(xyzzyaaap1(1:xyzzyaaab54-1,xyzzyaaaa54)==xyzzyaaae54).and.&
&.not.any(which_fam(xyzzyaaam54+1:xyzzyaaab54+xyzzyaaam54-1)==which_fa&
&m(xyzzyaaab54+xyzzyaaam54)))then
xyzzyaaaf54=count(xyzzyaaap1(1:xyzzyaaag1(xyzzyaaaa54),xyzzyaaaa54)==x&
&yzzyaaae54.and.which_fam(xyzzyaaam54+1:xyzzyaaag1(xyzzyaaaa54)+xyzzya&
&aam54)==which_fam(xyzzyaaab54+xyzzyaaam54))
xyzzyaaah53(xyzzyaaaa54)=xyzzyaaah53(xyzzyaaaa54)+(xyzzyaaaf54*(xyzzya&
&aaf54-1))/2
endif
enddo
enddo
allocate(xyzzyaaaq53(xyzzyaaae1,xyzzyaaae1,maxval(xyzzyaaah53),2),xyzz&
&yaaar53(maxval(xyzzyaaah53),2),stat=xyzzyaaaa53)
call check_alloc(xyzzyaaaa53,'SYMMETRIZE_TERM','trans_matrix')
xyzzyaaaq53=0
xyzzyaaar53=0
do xyzzyaaah54=1,xyzzyaaae1
xyzzyaaaq53(xyzzyaaah54,xyzzyaaah54,:,:)=1
enddo
do xyzzyaaaa54=1,2
xyzzyaaam54=0
if(xyzzyaaaa54==2)xyzzyaaam54=xyzzyaaag1(1)
xyzzyaaag54=0
do xyzzyaaab54=1,xyzzyaaag1(xyzzyaaaa54)
xyzzyaaae54=xyzzyaaap1(xyzzyaaab54,xyzzyaaaa54)
xyzzyaaaj54=sum(xyzzyaaap1(1:xyzzyaaab54-1,xyzzyaaaa54))
do xyzzyaaah54=1,xyzzyaaae54
do xyzzyaaai54=xyzzyaaah54+1,xyzzyaaae54
xyzzyaaag54=xyzzyaaag54+1
xyzzyaaar53(xyzzyaaag54,xyzzyaaaa54)=-1
xyzzyaaaq53(xyzzyaaah54+xyzzyaaaj54,xyzzyaaah54+xyzzyaaaj54,xyzzyaaag5&
&4,xyzzyaaaa54)=0
xyzzyaaaq53(xyzzyaaai54+xyzzyaaaj54,xyzzyaaai54+xyzzyaaaj54,xyzzyaaag5&
&4,xyzzyaaaa54)=0
xyzzyaaaq53(xyzzyaaah54+xyzzyaaaj54,xyzzyaaai54+xyzzyaaaj54,xyzzyaaag5&
&4,xyzzyaaaa54)=1
xyzzyaaaq53(xyzzyaaai54+xyzzyaaaj54,xyzzyaaah54+xyzzyaaaj54,xyzzyaaag5&
&4,xyzzyaaaa54)=1
enddo
enddo
if(.not.any(xyzzyaaap1(1:xyzzyaaab54-1,xyzzyaaaa54)==xyzzyaaae54).and.&
&.not.any(which_fam(xyzzyaaam54+1:xyzzyaaab54+xyzzyaaam54-1)==which_fa&
&m(xyzzyaaab54+xyzzyaaam54)))then
do xyzzyaaac54=xyzzyaaab54,xyzzyaaag1(xyzzyaaaa54)
if(xyzzyaaap1(xyzzyaaac54,xyzzyaaaa54)/=xyzzyaaae54.or.which_fam(xyzzy&
&aaac54+xyzzyaaam54)/=which_fam(xyzzyaaab54+xyzzyaaam54))cycle
xyzzyaaak54=sum(xyzzyaaap1(1:xyzzyaaac54-1,xyzzyaaaa54))
do xyzzyaaad54=xyzzyaaac54+1,xyzzyaaag1(xyzzyaaaa54)
if(xyzzyaaap1(xyzzyaaad54,xyzzyaaaa54)/=xyzzyaaae54.or.which_fam(xyzzy&
&aaad54+xyzzyaaam54)/=which_fam(xyzzyaaab54+xyzzyaaam54))cycle
xyzzyaaal54=sum(xyzzyaaap1(1:xyzzyaaad54-1,xyzzyaaaa54))
xyzzyaaag54=xyzzyaaag54+1
xyzzyaaar53(xyzzyaaag54,xyzzyaaaa54)=1
do xyzzyaaah54=1,xyzzyaaae54
xyzzyaaaq53(xyzzyaaah54+xyzzyaaak54,xyzzyaaah54+xyzzyaaak54,xyzzyaaag5&
&4,xyzzyaaaa54)=0
xyzzyaaaq53(xyzzyaaah54+xyzzyaaal54,xyzzyaaah54+xyzzyaaal54,xyzzyaaag5&
&4,xyzzyaaaa54)=0
xyzzyaaaq53(xyzzyaaah54+xyzzyaaak54,xyzzyaaah54+xyzzyaaal54,xyzzyaaag5&
&4,xyzzyaaaa54)=1
xyzzyaaaq53(xyzzyaaah54+xyzzyaaal54,xyzzyaaah54+xyzzyaaak54,xyzzyaaag5&
&4,xyzzyaaaa54)=1
enddo
enddo
enddo
endif
enddo
enddo
xyzzyaaav53=ne==xyzzyaaad1.and.xyzzyaaay1.and.all(xyzzyaaap1(:,1)==xyz&
&zyaaap1(:,2))
end subroutine xyzzyaaba53
subroutine xyzzyaabb53(xyzzyaaag53,matrix,nterm)
implicit none
integer,intent(in) :: xyzzyaaag53,matrix(ne,xyzzyaaad1)
integer,intent(out) :: nterm
integer xyzzyaaaa55,xyzzyaaab55,n,xyzzyaaac55(ne,xyzzyaaad1),xyzzyaaad&
&55,xyzzyaaae55,xyzzyaaaf55
logical xyzzyaaag55
type(perm),pointer :: xyzzyaaah55,xyzzyaaai55,xyzzyaaaj55,xyzzyaaak55
allocate(xyzzyaaai55,stat=xyzzyaaaa55)
call check_alloc(xyzzyaaaa55,'SYMMETRIZE_TERM','curr_perm')
allocate(xyzzyaaai55%matrix(ne,xyzzyaaad1),stat=xyzzyaaaa55)
call check_alloc(xyzzyaaaa55,'SYMMETRIZE_TERM','matrix')
xyzzyaaai55%matrix=matrix
xyzzyaaai55%csign=xyzzyaaag53
xyzzyaaai55%next=>xyzzyaaai55
nterm=1
if(.not.associated(xyzzyaaaz53))then
if(xyzzyaaag53/=1)call errstop('SYMMETRIZE_TERM','Not starting at 1? W&
&hy?')
xyzzyaaaz53=>xyzzyaaai55
xyzzyaaai55%next_term=>xyzzyaaaz53
else
xyzzyaaak55=>xyzzyaaaz53
xyzzyaaab55=1
do while(.not.associated(xyzzyaaak55%next_term,xyzzyaaaz53))
xyzzyaaak55=>xyzzyaaak55%next_term
xyzzyaaab55=xyzzyaaab55+1
enddo
if(xyzzyaaab55+1/=xyzzyaaag53)call errstop('SYMMETRIZE_TERM','Not in o&
&rder? Why?')
xyzzyaaak55%next_term=>xyzzyaaai55
nullify(xyzzyaaak55)
endif
xyzzyaaai55%next_term=>xyzzyaaaz53
xyzzyaaah55=>xyzzyaaai55
nullify(xyzzyaaai55)
xyzzyaaai55=>xyzzyaaah55
xyzzyaaaj55=>xyzzyaaah55
do
do xyzzyaaae55=1,3
select case(xyzzyaaae55)
case(1)
n=xyzzyaaah53(1)
case(2)
n=xyzzyaaah53(2)
case(3)
n=0
if(xyzzyaaav53)n=1
end select
do xyzzyaaaf55=1,n
select case(xyzzyaaae55)
case(1)
xyzzyaaac55=matmul(xyzzyaaaq53(1:ne,1:ne,xyzzyaaaf55,1),xyzzyaaai55%ma&
&trix)
xyzzyaaad55=xyzzyaaai55%csign*xyzzyaaar53(xyzzyaaaf55,1)
case(2)
xyzzyaaac55=matmul(xyzzyaaai55%matrix,xyzzyaaaq53(1:xyzzyaaad1,1:xyzzy&
&aaad1,xyzzyaaaf55,2))
xyzzyaaad55=xyzzyaaai55%csign*xyzzyaaar53(xyzzyaaaf55,2)
case(3)
xyzzyaaac55=transpose(xyzzyaaai55%matrix)
xyzzyaaad55=xyzzyaaai55%csign
end select
xyzzyaaak55=>xyzzyaaah55
xyzzyaaag55=.false.
do
if(all(xyzzyaaak55%matrix==xyzzyaaac55))then
if(xyzzyaaad55/=xyzzyaaak55%csign)then
if(am_master)then
call wout('Matrix with contradictory signs encountered:')
do xyzzyaaae53=1,ne
call write_list_int(xyzzyaaad1,xyzzyaaac55(xyzzyaaae53,1:xyzzyaaad1),x&
&yzzyaaad1,4,1)
enddo
endif
call errstop_master('SYMMETRIZE_TERM','Symmetry problem with term '//t&
&rim(i2s(xyzzyaaag53))//': the above matrix was generated twice via sy&
&mmetry operations, and it appears with opposite signs depending on ho&
&w it is generated.  This would suggest that this term is zero.  Check&
& that the matrix you have provided makes physical sense.')
endif
xyzzyaaag55=.true.
exit
endif
xyzzyaaak55=>xyzzyaaak55%next
if(associated(xyzzyaaak55,xyzzyaaah55))exit
enddo
nullify(xyzzyaaak55)
if(.not.xyzzyaaag55)then
nterm=nterm+1
allocate(xyzzyaaak55,stat=xyzzyaaaa55)
call check_alloc(xyzzyaaaa55,'SYMMETRIZE_TERM','tperm')
allocate(xyzzyaaak55%matrix(ne,xyzzyaaad1),stat=xyzzyaaaa55)
call check_alloc(xyzzyaaaa55,'SYMMETRIZE_TERM','matrix in tperm')
xyzzyaaak55%matrix=xyzzyaaac55
xyzzyaaak55%csign=xyzzyaaad55
xyzzyaaak55%next=>xyzzyaaah55
nullify(xyzzyaaak55%next_term)
xyzzyaaaj55%next=>xyzzyaaak55
xyzzyaaaj55=>xyzzyaaak55
nullify(xyzzyaaak55)
endif
enddo
enddo
if(associated(xyzzyaaai55,xyzzyaaaj55))exit
xyzzyaaai55=>xyzzyaaai55%next
enddo
nullify(xyzzyaaai55,xyzzyaaah55,xyzzyaaaj55)
end subroutine xyzzyaabb53
subroutine xyzzyaabc53(nterm_ind,nterm,matrices,signs)
implicit none
integer,intent(in) :: nterm_ind,nterm
integer,intent(out) :: matrices(ne,xyzzyaaad1,nterm),signs(nterm)
integer xyzzyaaaa56,xyzzyaaab56,xyzzyaaac56,xyzzyaaad56(nterm_ind)
type(perm),pointer :: xyzzyaaae56,xyzzyaaaf56,next_term,xyzzyaaag56
xyzzyaaae56=>xyzzyaaaz53
xyzzyaaaa56=0
xyzzyaaab56=0
do
xyzzyaaab56=xyzzyaaab56+1
xyzzyaaaf56=>xyzzyaaae56
xyzzyaaac56=xyzzyaaaa56
do
xyzzyaaaa56=xyzzyaaaa56+1
matrices(:,:,xyzzyaaaa56)=xyzzyaaaf56%matrix(:,:)
signs(xyzzyaaaa56)=xyzzyaaaf56%csign
if(abs(signs(xyzzyaaaa56))/=xyzzyaaab56)call errstop_master('GET_SYMME&
&TRIZED_TERMS','Bad coefficient found for term '//trim(i2s(xyzzyaaaa56&
&))//'.')
if(associated(xyzzyaaaf56%next,xyzzyaaae56))exit
xyzzyaaaf56=>xyzzyaaaf56%next
enddo
xyzzyaaad56(xyzzyaaab56)=xyzzyaaaa56-xyzzyaaac56
nullify(xyzzyaaaf56)
if(associated(xyzzyaaae56%next_term,xyzzyaaaz53))then
if(xyzzyaaaa56/=nterm)call errstop_master('GET_SYMMETRIZED_TERMS','Lis&
&t has less than NTERMS terms. Bug (1).')
if(xyzzyaaab56/=nterm_ind)call errstop_master('GET_SYMMETRIZED_TERMS',&
&'List has less than NTERMS terms. Bug (2).')
exit
endif
xyzzyaaae56=>xyzzyaaae56%next_term
enddo
xyzzyaaae56=>xyzzyaaaz53
do xyzzyaaab56=1,nterm_ind
next_term=>xyzzyaaae56%next_term
xyzzyaaaf56=>xyzzyaaae56
do xyzzyaaaa56=1,xyzzyaaad56(xyzzyaaab56)
xyzzyaaag56=>xyzzyaaaf56
xyzzyaaaf56=>xyzzyaaaf56%next
deallocate(xyzzyaaag56%matrix)
deallocate(xyzzyaaag56)
enddo
xyzzyaaae56=>next_term
enddo
nullify(xyzzyaaae56,next_term,xyzzyaaaf56,xyzzyaaag56,xyzzyaaaz53)
deallocate(xyzzyaaaq53,xyzzyaaar53)
end subroutine xyzzyaabc53
end subroutine read_exmol
subroutine write_exmol(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
integer xyzzyaaaa57,xyzzyaaab57,xyzzyaaac57,xyzzyaaad57,xyzzyaaae57,xy&
&zzyaaaf57
logical xyzzyaaag57
if(.not.am_master)return
inquire(file=trim(correlation_name),exist=xyzzyaaag57)
if(xyzzyaaag57)then
open(unit=xyzzyaabf1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaaa57)
else
open(unit=xyzzyaabf1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaaa57)
endif
if(xyzzyaaaa57/=0)call errstop('WRITE_EXMOL','Problem opening '//trim(&
&correlation_name)//'.')
write(xyzzyaabf1,*)'START EXMOL'
write(xyzzyaabf1,*)'Title'
write(xyzzyaabf1,*)trim(adjustl(title))
write(xyzzyaabf1,*)'Cusp in orbitals'
write(xyzzyaabf1,'(3x,a)')l2s(xyzzyaaaz1)
write(xyzzyaabf1,*)'Number of terms'
write(xyzzyaabf1,'(3x,a)')trim(i2s(xyzzyaaaf1))
do xyzzyaaad57=1,xyzzyaaaf1
write(xyzzyaabf1,*)'START TERM ',trim(i2s(xyzzyaaad57))
write(xyzzyaabf1,*)'Function selectors for all particle pairs'
do xyzzyaaaf57=1,xyzzyaaaa1
if(abs(xyzzyaaah1(xyzzyaaaf57))/=xyzzyaaad57)cycle
do xyzzyaaae57=1,ne
call write_list_int(xyzzyaaad1,xyzzyaaai1(xyzzyaaae57,1:xyzzyaaad1,xyz&
&zyaaaf57),xyzzyaaad1,4,1,xyzzyaabf1)
enddo
exit
enddo
write(xyzzyaabf1,*)'Coefficient ;        Optimizable (0=NO; 1=YES)'
write(xyzzyaabf1,*)xyzzyaaaw1(xyzzyaaad57)/sqrt(sum(xyzzyaaaw1(:)**2))&
&,xyzzyaaam1(xyzzyaaad57),'      ! c_',trim(i2s(xyzzyaaad57))
write(xyzzyaabf1,*)'END TERM ',trim(i2s(xyzzyaaad57))
enddo
do xyzzyaaac57=1,xyzzyaaac1
write(xyzzyaabf1,*)'START FUNCTION ',trim(i2s(xyzzyaaac57))
write(xyzzyaabf1,*)'Functional form'
write(xyzzyaabf1,'(3x,a)')trim(adjustl(xyzzyaabb1(xyzzyaaac57)))
select case(xyzzyaaaj1(xyzzyaaac57))
case(1)
continue
case(2,3)
write(xyzzyaabf1,*)'Truncation order'
write(xyzzyaabf1,'(3x,a)')trim(i2s(xyzzyaaaq1(1,xyzzyaaac57)))
write(xyzzyaabf1,*)'Expansion order'
write(xyzzyaabf1,'(3x,a)')trim(i2s(xyzzyaaaq1(2,xyzzyaaac57)))
case(4)
continue
case(5)
continue
case(6)
continue
end select
write(xyzzyaabf1,*)'Parameter ;          Optimizable (0=NO; 1=YES)'
do xyzzyaaab57=1,xyzzyaaak1(xyzzyaaac57)
if(xyzzyaaal1(xyzzyaaab57,xyzzyaaac57)/=-1)then
write(xyzzyaabf1,*)xyzzyaaav1(xyzzyaaab57,xyzzyaaac57),xyzzyaaal1(xyzz&
&yaaab57,xyzzyaaac57),'      ! p_',trim(i2s(xyzzyaaab57)),',',trim(i2s&
&(xyzzyaaac57))
endif
enddo
write(xyzzyaabf1,*)'END FUNCTION ',trim(i2s(xyzzyaaac57))
enddo
write(xyzzyaabf1,*)'END EXMOL'
write(xyzzyaabf1,*)
close(xyzzyaabf1)
end subroutine write_exmol
subroutine xyzzyaadj1(ii,renorm_term,gradmat,fi)
implicit none
integer,intent(in) :: ii
real(dp),intent(in) :: renorm_term(xyzzyaaaa1),gradmat(3,ne,xyzzyaaad1&
&,xyzzyaaac1)
complex(dp),intent(out) :: fi(3)
integer xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58,xyzzyaaae58
real(dp) xyzzyaaaf58(3),xyzzyaaag58,xyzzyaaah58(3),xyzzyaaai58
xyzzyaaab58=xyzzyaaao1(ii)
xyzzyaaai58=xyzzyaaax1(ii)
xyzzyaaaf58=0.d0
xyzzyaaag58=0.d0
do xyzzyaaaa58=1,xyzzyaaaa1
xyzzyaaah58=0.d0
if(xyzzyaaba1(ii))then
do xyzzyaaac58=1,ne
xyzzyaaae58=xyzzyaaai1(xyzzyaaac58,xyzzyaaab58,xyzzyaaaa58)
if(xyzzyaaae58/=0)xyzzyaaah58(1:3)=xyzzyaaah58(1:3)+xyzzyaaai58*gradma&
&t(1:3,xyzzyaaac58,xyzzyaaab58,xyzzyaaae58)
enddo
else
do xyzzyaaad58=1,xyzzyaaad1
xyzzyaaae58=xyzzyaaai1(xyzzyaaab58,xyzzyaaad58,xyzzyaaaa58)
if(xyzzyaaae58/=0)xyzzyaaah58(1:3)=xyzzyaaah58(1:3)+xyzzyaaai58*gradma&
&t(1:3,xyzzyaaab58,xyzzyaaad58,xyzzyaaae58)
enddo
endif
xyzzyaaaf58(1:3)=xyzzyaaaf58(1:3)+xyzzyaaah58(1:3)*renorm_term(xyzzyaa&
&aa58)
xyzzyaaag58=xyzzyaaag58+renorm_term(xyzzyaaaa58)
enddo
fi(1:3)=cmplx(xyzzyaaaf58(1:3)/xyzzyaaag58,0.d0,dp)
end subroutine xyzzyaadj1
subroutine xyzzyaadk1(ii,renorm_term,gradmat,lapmat,fi,ti)
implicit none
integer,intent(in) :: ii
real(dp),intent(in) :: renorm_term(xyzzyaaaa1),gradmat(3,ne,xyzzyaaad1&
&,xyzzyaaac1),lapmat(ne,xyzzyaaad1,xyzzyaaac1)
complex(dp),intent(in) :: fi(3)
complex(dp),intent(out) :: ti
integer xyzzyaaaa59,xyzzyaaab59,xyzzyaaac59,xyzzyaaad59,xyzzyaaae59
real(dp) xyzzyaaaf59,xyzzyaaag59,xyzzyaaah59(3),xyzzyaaai59,xyzzyaaaj5&
&9
xyzzyaaaj59=dble(fi(1))**2+dble(fi(2))**2+dble(fi(3))**2
if(complex_wf)xyzzyaaaj59=xyzzyaaaj59+aimag(fi(1))**2+aimag(fi(2))**2+&
&aimag(fi(3))**2
xyzzyaaab59=xyzzyaaao1(ii)
xyzzyaaaf59=0.d0
xyzzyaaai59=0.d0
do xyzzyaaaa59=1,xyzzyaaaa1
xyzzyaaag59=0.d0
xyzzyaaah59=0.d0
if(xyzzyaaba1(ii))then
do xyzzyaaac59=1,ne
xyzzyaaae59=xyzzyaaai1(xyzzyaaac59,xyzzyaaab59,xyzzyaaaa59)
if(xyzzyaaae59/=0)then
xyzzyaaag59=xyzzyaaag59+lapmat(xyzzyaaac59,xyzzyaaab59,xyzzyaaae59)
xyzzyaaah59(1:3)=xyzzyaaah59(1:3)+gradmat(1:3,xyzzyaaac59,xyzzyaaab59,&
&xyzzyaaae59)
endif
enddo
else
do xyzzyaaad59=1,xyzzyaaad1
xyzzyaaae59=xyzzyaaai1(xyzzyaaab59,xyzzyaaad59,xyzzyaaaa59)
if(xyzzyaaae59/=0)then
xyzzyaaag59=xyzzyaaag59+lapmat(xyzzyaaab59,xyzzyaaad59,xyzzyaaae59)
xyzzyaaah59(1:3)=xyzzyaaah59(1:3)+gradmat(1:3,xyzzyaaab59,xyzzyaaad59,&
&xyzzyaaae59)
endif
enddo
endif
xyzzyaaaf59=xyzzyaaaf59+(xyzzyaaag59+ddot(3,xyzzyaaah59(1),1,xyzzyaaah&
&59(1),1))*renorm_term(xyzzyaaaa59)
xyzzyaaai59=xyzzyaaai59+renorm_term(xyzzyaaaa59)
enddo
ti=cmplx(xyzzyaaaf59/xyzzyaaai59-xyzzyaaaj59,0.d0,dp)
end subroutine xyzzyaadk1
subroutine xyzzyaadl1(renorm_term,gradmat,farray,farray_term)
implicit none
real(dp),intent(in) :: renorm_term(xyzzyaaaa1),gradmat(3,ne,xyzzyaaad1&
&,xyzzyaaac1)
real(dp),intent(inout) :: farray(3,real1_complex2,netot),farray_term(3&
&,real1_complex2,xyzzyaaaa1,netot)
integer xyzzyaaaa60,xyzzyaaab60,xyzzyaaac60,xyzzyaaad60,xyzzyaaae60,xy&
&zzyaaaf60,xyzzyaaag60
real(dp) xyzzyaaah60(3),xyzzyaaai60,xyzzyaaaj60,xyzzyaaak60
if(.not.complex_wf)then
do xyzzyaaaa60=1,netot
xyzzyaaab60=xyzzyaaao1(xyzzyaaaa60)
xyzzyaaak60=xyzzyaaax1(xyzzyaaaa60)
do xyzzyaaaf60=1,xyzzyaaaa1
xyzzyaaah60=0.d0
if(xyzzyaaba1(xyzzyaaaa60))then
do xyzzyaaac60=1,ne
xyzzyaaae60=xyzzyaaai1(xyzzyaaac60,xyzzyaaab60,xyzzyaaaf60)
if(xyzzyaaae60/=0)xyzzyaaah60(1:3)=xyzzyaaah60(1:3)+xyzzyaaak60*gradma&
&t(1:3,xyzzyaaac60,xyzzyaaab60,xyzzyaaae60)
enddo
else
do xyzzyaaad60=1,xyzzyaaad1
xyzzyaaae60=xyzzyaaai1(xyzzyaaab60,xyzzyaaad60,xyzzyaaaf60)
if(xyzzyaaae60/=0)xyzzyaaah60(1:3)=xyzzyaaah60(1:3)+xyzzyaaak60*gradma&
&t(1:3,xyzzyaaab60,xyzzyaaad60,xyzzyaaae60)
enddo
endif
farray_term(1:3,1,xyzzyaaaf60,xyzzyaaaa60)=xyzzyaaah60(1:3)
enddo
enddo
else
call errstop('EVAL_FARRAY','Complex EXMOL wfn not supported.')
endif
if(.not.complex_wf)then
xyzzyaaai60=sum(renorm_term(1:xyzzyaaaa1))
if(xyzzyaaai60==0.d0)call errstop('EVAL_FARRAY','About to divide by ze&
&ro.')
xyzzyaaaj60=1.d0/xyzzyaaai60
else
call errstop('EVAL_FARRAY','Complex EXMOL wfn not supported.')
endif
if(.not.complex_wf)then
do xyzzyaaaa60=1,netot
do xyzzyaaag60=1,dimensionality
farray(xyzzyaaag60,1,xyzzyaaaa60)=ddot(xyzzyaaaa1,renorm_term(1),1,far&
&ray_term(xyzzyaaag60,1,1,xyzzyaaaa60),3)*xyzzyaaaj60
enddo
enddo
else
call errstop('EVAL_FARRAY','Complex EXMOL wfn not supported.')
endif
end subroutine xyzzyaadl1
subroutine xyzzyaadm1(renorm_term,d2mat,bf_m2,bf_rmap2,farray_term,har&
&ray)
implicit none
integer,intent(in) :: bf_m2(netot),bf_rmap2(netot,netot)
real(dp),intent(in) :: renorm_term(xyzzyaaaa1),d2mat(6,ne,xyzzyaaad1,x&
&yzzyaaac1),farray_term(3,real1_complex2,xyzzyaaaa1,netot)
real(dp),intent(inout) :: harray(3,3,real1_complex2,netot,netot)
integer xyzzyaaaa61,xyzzyaaab61,xyzzyaaac61,xyzzyaaad61,xyzzyaaae61,xy&
&zzyaaaf61,xyzzyaaag61,xyzzyaaah61,xyzzyaaai61,xyzzyaaaj61,xyzzyaaak61&
&,xyzzyaaal61
real(dp) xyzzyaaam61,xyzzyaaan61,xyzzyaaao61,xyzzyaaap61,xyzzyaaaq61
if(.not.complex_wf)then
xyzzyaaao61=sum(renorm_term(1:xyzzyaaaa1))
if(xyzzyaaao61==0.d0)call errstop('EVAL_HARRAY','About to divide by ze&
&ro.')
xyzzyaaap61=1.d0/xyzzyaaao61
else
call errstop('EVAL_HARRAY','Complex EXMOL wfn not supported.')
endif
if(.not.complex_wf)then
do xyzzyaaad61=1,netot
xyzzyaaae61=xyzzyaaao1(xyzzyaaad61)
do xyzzyaaaf61=1,bf_m2(xyzzyaaad61)
xyzzyaaag61=bf_rmap2(xyzzyaaaf61,xyzzyaaad61)
if(xyzzyaaag61<xyzzyaaad61)cycle
xyzzyaaah61=xyzzyaaao1(xyzzyaaag61)
do xyzzyaaaa61=1,dimensionality
do xyzzyaaab61=1,dimensionality
xyzzyaaac61=which_d2index(xyzzyaaaa61,xyzzyaaab61)
xyzzyaaam61=0.d0
do xyzzyaaal61=1,xyzzyaaaa1
xyzzyaaan61=0.d0
xyzzyaaan61=xyzzyaaan61+farray_term(xyzzyaaaa61,1,xyzzyaaal61,xyzzyaaa&
&d61)*farray_term(xyzzyaaab61,1,xyzzyaaal61,xyzzyaaag61)
if(xyzzyaaad61==xyzzyaaag61)then
if(xyzzyaaba1(xyzzyaaad61))then
do xyzzyaaai61=1,ne
xyzzyaaak61=xyzzyaaai1(xyzzyaaai61,xyzzyaaae61,xyzzyaaal61)
xyzzyaaan61=xyzzyaaan61+d2mat(xyzzyaaac61,xyzzyaaai61,xyzzyaaae61,xyzz&
&yaaak61)
enddo
else
do xyzzyaaaj61=1,xyzzyaaad1
xyzzyaaak61=xyzzyaaai1(xyzzyaaae61,xyzzyaaaj61,xyzzyaaal61)
xyzzyaaan61=xyzzyaaan61+d2mat(xyzzyaaac61,xyzzyaaae61,xyzzyaaaj61,xyzz&
&yaaak61)
enddo
endif
elseif(xyzzyaaba1(xyzzyaaad61).neqv.xyzzyaaba1(xyzzyaaag61))then
if(xyzzyaaba1(xyzzyaaad61))then
xyzzyaaak61=xyzzyaaai1(xyzzyaaah61,xyzzyaaae61,xyzzyaaal61)
xyzzyaaan61=xyzzyaaan61-d2mat(xyzzyaaac61,xyzzyaaah61,xyzzyaaae61,xyzz&
&yaaak61)
else
xyzzyaaak61=xyzzyaaai1(xyzzyaaae61,xyzzyaaah61,xyzzyaaal61)
xyzzyaaan61=xyzzyaaan61-d2mat(xyzzyaaac61,xyzzyaaae61,xyzzyaaah61,xyzz&
&yaaak61)
endif
endif
xyzzyaaam61=xyzzyaaam61+renorm_term(xyzzyaaal61)*xyzzyaaan61
enddo
xyzzyaaaq61=xyzzyaaam61*xyzzyaaap61
harray(xyzzyaaaa61,xyzzyaaab61,1,xyzzyaaad61,xyzzyaaag61)=xyzzyaaaq61
if(xyzzyaaad61/=xyzzyaaag61)harray(xyzzyaaab61,xyzzyaaaa61,1,xyzzyaaag&
&61,xyzzyaaad61)=xyzzyaaaq61
enddo
enddo
enddo
enddo
else
call errstop('EVAL_HARRAY','Complex EXMOL wave function not supported.&
&')
endif
end subroutine xyzzyaadm1
subroutine xyzzyaadn1(eevecs,val,fd,sd,orbmat,gradmat,lapmat,d2mat)
implicit none
logical,intent(in) :: val,fd,sd
real(dp),intent(in) :: eevecs(4,xyzzyaaab1,xyzzyaaab1)
real(dp),intent(inout),target,optional :: orbmat(ne,xyzzyaaad1,xyzzyaa&
&ac1),gradmat(3,ne,xyzzyaaad1,xyzzyaaac1),lapmat(ne,xyzzyaaad1,xyzzyaa&
&ac1),d2mat(6,ne,xyzzyaaad1,xyzzyaaac1)
integer xyzzyaaaa62,xyzzyaaab62
real(dp),pointer :: xyzzyaaac62(:,:),xyzzyaaad62(:,:,:),xyzzyaaae62(:,&
&:),xyzzyaaaf62(:,:,:)
logical xyzzyaaag62,xyzzyaaah62
xyzzyaaag62=present(d2mat).and.sd
xyzzyaaah62=present(lapmat).and.sd
if(.not.val)xyzzyaaac62=>xyzzyaacg1
if(.not.fd)xyzzyaaad62=>xyzzyaach1
if(.not.xyzzyaaah62)xyzzyaaae62=>xyzzyaaci1
if(.not.xyzzyaaag62)xyzzyaaaf62=>xyzzyaacj1
do xyzzyaaab62=1,xyzzyaaad1
xyzzyaaaa62=xyzzyaaan1(xyzzyaaab62,2)
if(val)xyzzyaaac62=>orbmat(:,xyzzyaaab62,:)
if(fd)xyzzyaaad62=>gradmat(:,:,xyzzyaaab62,:)
if(xyzzyaaah62)xyzzyaaae62=>lapmat(:,xyzzyaaab62,:)
if(xyzzyaaag62)xyzzyaaaf62=>d2mat(:,:,xyzzyaaab62,:)
call xyzzyaado1(xyzzyaaaa62,xyzzyaaat1(xyzzyaaaa62),eevecs(1,1,xyzzyaa&
&aa62),val,fd,sd,xyzzyaaac62,xyzzyaaad62,xyzzyaaae62,xyzzyaaaf62)
enddo
end subroutine xyzzyaadn1
subroutine xyzzyaado1(ii,norb,eevecs1,val,fd,sd,orbvec,gradvec,lapvec,&
&d2vec)
implicit none
integer,intent(in) :: ii,norb
real(dp),intent(in) :: eevecs1(4,xyzzyaaab1)
real(dp),intent(inout),optional :: orbvec(norb,xyzzyaaac1),gradvec(3,n&
&orb,xyzzyaaac1),lapvec(norb,xyzzyaaac1),d2vec(6,norb,xyzzyaaac1)
logical,intent(in) :: val,fd,sd
integer xyzzyaaaa63,xyzzyaaab63,xyzzyaaac63,xyzzyaaad63
real(dp) xyzzyaaae63,xyzzyaaaf63,xyzzyaaag63,xyzzyaaah63(3),xyzzyaaai6&
&3,xyzzyaaaj63,xyzzyaaak63,xyzzyaaal63(6),xyzzyaaam63(6),xyzzyaaan63
logical xyzzyaaao63,xyzzyaaap63,xyzzyaaaq63
xyzzyaaao63=fd.or.sd
xyzzyaaap63=present(d2vec).and.sd
xyzzyaaaq63=present(lapvec).and.sd
xyzzyaaan63=xyzzyaaax1(ii)
xyzzyaaad63=xyzzyaaas1(ii)
do xyzzyaaac63=1,norb
xyzzyaaab63=xyzzyaaan1(xyzzyaaac63,xyzzyaaad63)
xyzzyaaak63=eevecs1(4,xyzzyaaab63)
if(xyzzyaaao63)call xyzzyaadx1(xyzzyaaak63,eevecs1(1:3,xyzzyaaab63),xy&
&zzyaaah63,xyzzyaaai63,xyzzyaaaj63,xyzzyaaal63,xyzzyaaam63)
do xyzzyaaaa63=1,xyzzyaaac1
call xyzzyaadp1(xyzzyaaaa63,xyzzyaaak63,xyzzyaaao63,sd,xyzzyaaae63,xyz&
&zyaaaf63,xyzzyaaag63)
if(val)orbvec(xyzzyaaac63,xyzzyaaaa63)=xyzzyaaae63
if(fd)gradvec(1:3,xyzzyaaac63,xyzzyaaaa63)=xyzzyaaan63*xyzzyaaaf63*xyz&
&zyaaah63(1:3)
if(xyzzyaaaq63)lapvec(xyzzyaaac63,xyzzyaaaa63)=xyzzyaaag63*xyzzyaaai63&
&+xyzzyaaaf63*xyzzyaaaj63
if(xyzzyaaap63)d2vec(1:6,xyzzyaaac63,xyzzyaaaa63)=xyzzyaaag63*xyzzyaaa&
&l63(1:6)+xyzzyaaaf63*xyzzyaaam63(1:6)
enddo
enddo
end subroutine xyzzyaado1
subroutine xyzzyaadp1(set,r,fd,sd,f,df,d2f)
implicit none
integer,intent(in) :: set
real(dp),intent(in) :: r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fd,sd
integer xyzzyaaaa64
real(dp) xyzzyaaab64,xyzzyaaac64,xyzzyaaad64
f=0.d0
if(fd)df=0.d0
if(sd)d2f=0.d0
select case(xyzzyaaaj1(set))
case(1)
call xyzzyaadq1(xyzzyaaav1(1,set),xyzzyaaav1(2,set),r,fd,sd,f,df,d2f)
if(xyzzyaaaz1)then
call xyzzyaadr1(xyzzyaaav1(2,set),xyzzyaaav1(3,set),r,fd,sd,xyzzyaaab6&
&4,xyzzyaaac64,xyzzyaaad64)
f=f+xyzzyaaab64
if(fd)df=df+xyzzyaaac64
if(sd)d2f=d2f+xyzzyaaad64
endif
case(2)
call xyzzyaadq1(xyzzyaaav1(1,set),xyzzyaaav1(2,set),r,fd,sd,f,df,d2f)
if(xyzzyaaaz1)then
call xyzzyaadr1(xyzzyaaav1(2,set),xyzzyaaav1(3,set),r,fd,sd,xyzzyaaab6&
&4,xyzzyaaac64,xyzzyaaad64)
f=f+xyzzyaaab64
if(fd)df=df+xyzzyaaac64
if(sd)d2f=d2f+xyzzyaaad64
endif
xyzzyaaaa64=xyzzyaaaq1(2,set)
call xyzzyaads1(xyzzyaaaa64,xyzzyaaaq1(1,set),xyzzyaaav1(4,set),xyzzya&
&aav1(5:5+xyzzyaaaa64,set),r,fd,sd,xyzzyaaab64,xyzzyaaac64,xyzzyaaad64&
&)
f=f+xyzzyaaab64
if(fd)df=df+xyzzyaaac64
if(sd)d2f=d2f+xyzzyaaad64
case(3)
xyzzyaaaa64=xyzzyaaaq1(2,set)
call xyzzyaads1(xyzzyaaaa64,xyzzyaaaq1(1,set),xyzzyaaav1(1,set),xyzzya&
&aav1(2:2+xyzzyaaaa64,set),r,fd,sd,f,df,d2f)
case(4)
call xyzzyaadv1(xyzzyaaav1(1,set),r,fd,sd,f,df,d2f)
case(5)
call xyzzyaadw1(xyzzyaaav1(1,set),xyzzyaaav1(2,set),r,fd,sd,f,df,d2f)
case(6)
call xyzzyaadv1(xyzzyaaav1(1,set),r,fd,sd,f,df,d2f)
end select
end subroutine xyzzyaadp1
subroutine xyzzyaadq1(a,b,r,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: a,b,r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fd,sd
real(dp) xyzzyaaaa65,xyzzyaaab65,xyzzyaaac65
logical xyzzyaaad65
xyzzyaaad65=fd.or.sd
xyzzyaaaa65=1.d0/(b+r)
xyzzyaaab65=xyzzyaaaa65/a
f=-xyzzyaaab65*r*r
if(xyzzyaaad65)then
xyzzyaaac65=xyzzyaaab65*xyzzyaaaa65
if(fd)df=-xyzzyaaac65*r*(r+b+b)
if(sd)d2f=-xyzzyaaac65*xyzzyaaaa65*(b+b)*b
endif
end subroutine xyzzyaadq1
subroutine xyzzyaadr1(b,c,r,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: b,c,r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fd,sd
real(dp) xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66
logical xyzzyaaad66
xyzzyaaad66=fd.or.sd
xyzzyaaaa66=1.d0/(b+r)
xyzzyaaab66=xyzzyaaaa66/c
f=-xyzzyaaab66*r
if(xyzzyaaad66)then
xyzzyaaac66=xyzzyaaab66*xyzzyaaaa66
if(fd)df=-xyzzyaaac66*b
if(sd)d2f=xyzzyaaac66*xyzzyaaaa66*(b+b)
endif
end subroutine xyzzyaadr1
subroutine xyzzyaads1(n,c,l,coeffs,r,fd,sd,f,df,d2f)
implicit none
integer,intent(in) :: n,c
real(dp),intent(in) :: l,coeffs(0:n),r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fd,sd
real(dp) xyzzyaaaa67,xyzzyaaab67,xyzzyaaac67,xyzzyaaad67,xyzzyaaae67,x&
&yzzyaaaf67,xyzzyaaag67,xyzzyaaah67
logical xyzzyaaai67
if(r<l)then
xyzzyaaaa67=1.d0/l
xyzzyaaab67=r*xyzzyaaaa67
xyzzyaaai67=fd.or.sd
call xyzzyaadt1(n,coeffs(0:n),xyzzyaaab67,xyzzyaaai67,sd,xyzzyaaaf67,x&
&yzzyaaag67,xyzzyaaah67)
call xyzzyaadu1(c,xyzzyaaab67,xyzzyaaai67,sd,xyzzyaaac67,xyzzyaaad67,x&
&yzzyaaae67)
f=xyzzyaaaf67*xyzzyaaac67
if(fd)df=(xyzzyaaag67*xyzzyaaac67+xyzzyaaaf67*xyzzyaaad67)*xyzzyaaaa67
if(sd)d2f=(xyzzyaaah67*xyzzyaaac67+2.d0*xyzzyaaag67*xyzzyaaad67+xyzzya&
&aaf67*xyzzyaaae67)*xyzzyaaaa67*xyzzyaaaa67
else
f=0.d0
df=0.d0
d2f=0.d0
endif
end subroutine xyzzyaads1
subroutine xyzzyaadt1(n,c,r,fd,sd,f,df,d2f)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: c(0:n),r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fd,sd
integer xyzzyaaaa68
real(dp) xyzzyaaab68,xyzzyaaac68,xyzzyaaad68
if(.not.(fd.or.sd))then
f=c(0)
xyzzyaaab68=1.d0
do xyzzyaaaa68=1,n
xyzzyaaab68=xyzzyaaab68*r
f=f+c(xyzzyaaaa68)*xyzzyaaab68
enddo
elseif(.not.sd)then
f=c(0)
df=0.d0
xyzzyaaab68=1.d0
xyzzyaaac68=0.d0
do xyzzyaaaa68=1,n
xyzzyaaac68=xyzzyaaaa68*xyzzyaaab68
xyzzyaaab68=xyzzyaaab68*r
f=f+c(xyzzyaaaa68)*xyzzyaaab68
df=df+c(xyzzyaaaa68)*xyzzyaaac68
enddo
else
f=c(0)
df=0.d0
d2f=0.d0
xyzzyaaab68=1.d0
xyzzyaaac68=0.d0
xyzzyaaad68=0.d0
do xyzzyaaaa68=1,n
xyzzyaaad68=xyzzyaaaa68*xyzzyaaac68
xyzzyaaac68=xyzzyaaaa68*xyzzyaaab68
xyzzyaaab68=xyzzyaaab68*r
f=f+c(xyzzyaaaa68)*xyzzyaaab68
df=df+c(xyzzyaaaa68)*xyzzyaaac68
d2f=d2f+c(xyzzyaaaa68)*xyzzyaaad68
enddo
endif
end subroutine xyzzyaadt1
subroutine xyzzyaadu1(c,r,fd,sd,f,df,d2f)
implicit none
integer,intent(in) :: c
real(dp),intent(in) :: r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fd,sd
real(dp) xyzzyaaaa69,xyzzyaaab69,xyzzyaaac69
xyzzyaaaa69=1.d0-r
if(.not.(fd.or.sd))then
select case(c)
case(1)
f=xyzzyaaaa69
case(2)
f=xyzzyaaaa69*xyzzyaaaa69
case default
f=xyzzyaaaa69**c
end select
elseif(.not.sd)then
select case(c)
case(1)
f=xyzzyaaaa69
df=-1.d0
case(2)
f=xyzzyaaaa69*xyzzyaaaa69
df=-xyzzyaaaa69-xyzzyaaaa69
case default
xyzzyaaab69=xyzzyaaaa69**(c-1)
f=xyzzyaaab69*xyzzyaaaa69
df=-dble(c)*xyzzyaaab69
end select
else
select case(c)
case(1)
f=xyzzyaaaa69
df=-1.d0
d2f=0.d0
case(2)
f=xyzzyaaaa69*xyzzyaaaa69
df=-xyzzyaaaa69-xyzzyaaaa69
d2f=2.d0
case default
xyzzyaaac69=xyzzyaaaa69**(c-2)
xyzzyaaab69=xyzzyaaac69*xyzzyaaaa69
f=xyzzyaaab69*xyzzyaaaa69
df=-dble(c)*xyzzyaaab69
d2f=dble(c*c-c)*xyzzyaaac69
end select
endif
end subroutine xyzzyaadu1
subroutine xyzzyaadv1(a,r,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: a,r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fd,sd
logical xyzzyaaaa70
xyzzyaaaa70=fd.or.sd
f=-a*r
if(xyzzyaaaa70)then
if(fd)df=-a
if(sd)d2f=0.d0
endif
end subroutine xyzzyaadv1
subroutine xyzzyaadw1(a,b,r,fd,sd,f,df,d2f)
implicit none
real(dp),intent(in) :: a,b,r
real(dp),intent(out) :: f,df,d2f
logical,intent(in) :: fd,sd
logical xyzzyaaaa71
real(dp) xyzzyaaab71
xyzzyaaaa71=fd.or.sd
xyzzyaaab71=r*r
f=-(b+xyzzyaaab71)/(a*r)
if(xyzzyaaaa71)then
if(fd)df=-(1.d0-b/xyzzyaaab71)/a
if(sd)d2f=-2.d0*b/(a*r*xyzzyaaab71)
endif
end subroutine xyzzyaadw1
subroutine xyzzyaadx1(r,vecr,grad,grad2,lap,drdr,d2r)
implicit none
real(dp),intent(in) :: r,vecr(:)
real(dp),intent(out) :: grad(3),grad2,lap,drdr(6),d2r(6)
real(dp) xyzzyaaaa72
grad2=1.d0
if(r==0.d0)then
grad(1:3)=0.d0
lap=0.d0
drdr=0.d0
d2r=0.d0
return
endif
select case(dimensionality)
case(3)
xyzzyaaaa72=1.d0/r
lap=2*xyzzyaaaa72
grad(1:3)=vecr(1:3)*xyzzyaaaa72
drdr(1:3)=grad(1:3)*grad(1:3)
drdr(4:6)=(/grad(1)*grad(2),grad(1)*grad(3),grad(2)*grad(3)/)
d2r(1:3)=xyzzyaaaa72-drdr(1:3)*xyzzyaaaa72
d2r(4:6)=-drdr(4:6)*xyzzyaaaa72
case(2)
xyzzyaaaa72=1.d0/r
lap=xyzzyaaaa72
grad(1:2)=vecr(1:2)*xyzzyaaaa72
grad(3)=0.d0
drdr(1:2)=grad(1:2)*grad(1:2)
drdr(4)=grad(1)*grad(2)
drdr(3)=0.d0
drdr(5:6)=0.d0
d2r(1:2)=xyzzyaaaa72-drdr(1:2)*xyzzyaaaa72
d2r(4)=-drdr(4)*xyzzyaaaa72
d2r(3)=0.d0
d2r(5:6)=0.d0
case(1)
grad(1)=sign(1.d0,vecr(1))
grad(2:3)=0.d0
lap=0.d0
drdr(1)=1.d0
drdr(2:6)=0.d0
d2r(1:6)=0.d0
end select
end subroutine xyzzyaadx1
end module slaarnaao

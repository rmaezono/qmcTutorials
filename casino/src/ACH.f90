module slaarnach
use dsp
use store
use slaarnaan,only : ee_distances,ei_distances,ee_distances_all,ei_dis&
&tances_all
use format_utils, only : wout,i2s
use slaarnabg,     only : homogeneous_system,nitot,rion,ion_displaceme&
&nt,which_ion_displaced
use slaarnabt,    only : dcopy,dcopy_ee,dcopy_ee_flip,next_permutation
use run_control,  only : errstop,check_alloc
implicit none
private
public scratch_request,scratch_protect,setup_scratch,finish_scratch,wh&
&ich_scratch,define_config_scratch,define_config_oneelec_scratch,defin&
&e_config_twoelec_scratch,define_config_oneion_scratch,accept_move_scr&
&atch,reset_config_scratch,get_rsele,get_eevecs,get_eevecs1_ch,get_eiv&
&ecs,get_eivecs1_ch,clear_scratch,include_range,in_range,clone_scratch&
&_geom,setup_storage_geom,finish_storage_geom,load_from_storage_geom,p&
&ush_rion_hack
public gen_config_geom,delete_config_geom,copy_config_geom,config_to_p&
&t_geom,pt_to_config_geom,redist_allocations_geom,redist_load_geom,red&
&ist_send_geom,redist_recv_geom,redist_save_geom,redist_deallocations_&
&geom,load_from_pt_geom,save_to_pt_geom
public config_geom
public nscratch,scr_tasks,buffer_move1_from,buffer_move1_from_ii,buffe&
&r_move2_from,buffer_move2_from_ii,buffer_move2_from_jj,buffer_move_io&
&n_from,buffer_move_ion_from_ion,ratio1_from_sz,ratio1_to_sz,ratio2_fr&
&om_sz,ratio2_to_sz,ratiocfg_from_sz,ratiocfg_to_sz,ratio_ion_from_sz,&
&ratio_ion_to_sz,drift_sz,kinetic_sz,kinetic_detail_sz,wfn_detail_sz,k&
&inetic_aderiv_sz,wfn_aderiv_sz,iratio1_from,iratio1_to,iratio2_from,i&
&ratio2_to,iratiocfg_from,iratiocfg_to,iratio_ion_from,iratio_ion_to,i&
&drift,ikinetic,ikinetic_detail,iwfn_detail,ikinetic_aderiv,iwfn_aderi&
&v
public rele_scr,sele_scr,rsele_valid,eevecs_scr,eevecs_valid,eivecs_sc&
&r,eivecs_valid,rele1_chscr,sele1_chscr,rsele1_chvalid,rele2ii_chscr,r&
&ele2jj_chscr,rsele2_chvalid,rion1_chscr,rion1_chvalid,eevecs1_chscr,e&
&evecs1_chvalid,eivecs1_chscr,eivecs1_chvalid,eevecs1,eivecs1
public bf_x_scr,bf_x_valid,bf_dx_scr,bf_dx_valid,bf_d2x_scr,bf_d2x_val&
&id,eevecs_bf_scr,eevecs_bf_valid,bf_connect_scr,bf_m_scr,bf_rmap_scr,&
&bf_m2_scr,bf_rmap2_scr,bf_m_chscr,bf_rmap_chscr,bf_rmap_chvalid
integer nscratch
integer,allocatable :: buffer_move1_from(:),buffer_move1_from_ii(:),bu&
&ffer_move2_from(:),buffer_move2_from_ii(:),buffer_move2_from_jj(:),bu&
&ffer_move_ion_from(:),buffer_move_ion_from_ion(:)
integer ratio1_from_sz(2),ratio1_to_sz(2),ratio2_from_sz(2),ratio2_to_&
&sz(2),ratiocfg_from_sz(2),ratiocfg_to_sz(2),ratio_ion_from_sz(2),rati&
&o_ion_to_sz(2),drift_sz(2),kinetic_sz(2),kinetic_detail_sz(2),wfn_det&
&ail_sz(2),kinetic_aderiv_sz(2),wfn_aderiv_sz(2)
integer,allocatable :: sele_scr(:,:),sele1_chscr(:),xyzzyaaaa1(:),xyzz&
&yaaab1(:)
real(dp),allocatable,target :: rele_scr(:,:,:),rele1_chscr(:,:),rele2i&
&i_chscr(:,:),rele2jj_chscr(:,:),rion1_chscr(:,:),eevecs_scr(:,:,:,:),&
&eivecs_scr(:,:,:,:),eevecs1_chscr(:,:,:),eivecs1_chscr(:,:,:),eevecs1&
&(:,:),eivecs1(:,:)
logical,allocatable :: rsele_valid(:),eevecs_valid(:),eivecs_valid(:),&
&eevecs1_chvalid(:),eivecs1_chvalid(:),rsele1_chvalid(:),rsele2_chvali&
&d(:),rion1_chvalid(:)
integer,allocatable :: xyzzyaaac1(:),xyzzyaaad1(:),xyzzyaaae1(:)
integer,allocatable :: bf_m_scr(:,:),bf_rmap_scr(:,:,:),bf_m2_scr(:,:)&
&,bf_rmap2_scr(:,:,:),bf_m_chscr(:,:),bf_rmap_chscr(:,:,:)
logical,allocatable :: bf_connect_scr(:,:,:)
real(dp),allocatable :: bf_x_scr(:,:,:),bf_dx_scr(:,:,:,:,:),bf_d2x_sc&
&r(:,:,:,:)
real(dp),allocatable,target :: eevecs_bf_scr(:,:,:,:)
logical,allocatable :: bf_x_valid(:),bf_dx_valid(:),bf_d2x_valid(:),ee&
&vecs_bf_valid(:),bf_rmap_chvalid(:)
integer,parameter :: xyzzyaaaf1=14,xyzzyaaag1=50,iratio1_from=1,    ir&
&atio1_to=2,iratio2_from=3,    iratio2_to=4,iratiocfg_from=5,  iratioc&
&fg_to=6,iratio_ion_from=7, iratio_ion_to=8,idrift=9,          ikineti&
&c=10,ikinetic_detail=11,iwfn_detail=12,ikinetic_aderiv=13,iwfn_aderiv&
&=14
logical :: xyzzyaaah1(xyzzyaaaf1,xyzzyaaag1)=.false.,xyzzyaaai1(xyzzya&
&aag1,xyzzyaaag1)=.true.
logical,allocatable :: scr_tasks(:,:)
integer,allocatable :: xyzzyaaaj1(:)
integer :: xyzzyaaak1=0
character(30) task_name(xyzzyaaaf1)
data task_name /'ONE-ELEC. MOVE (FROM)',  'ONE-ELEC. MOVE (TO)','TWO-E&
&LEC. MOVE (FROM)',  'TWO-ELEC. MOVE (TO)','CONFIG MOVE (FROM)',     '&
&CONFIG MOVE (TO)','ONE-ION MOVE (FROM)',    'ONE-ION MOVE (TO)','DRIF&
&T VECTOR',           'KINETIC ENERGY','KINETIC ENERGY (DETAIL)','WAVE&
& FUNCTION (DETAIL)','DERIV OF KINETIC ENERGY','DERIV OF WAVE FUNCTION&
&'/
type config_geom
private
real(dp),pointer :: pt_rele(:,:)
end type config_geom
real(dp),allocatable :: xyzzyaaal1(:,:,:)
real(dp),allocatable :: xyzzyaaam1(:,:,:)
integer,allocatable :: xyzzyaaan1(:,:)
integer,allocatable :: xyzzyaaao1(:)
contains
subroutine scratch_request(ratio1_from,ratio1_to,ratio2_from,ratio2_to&
&,ratiocfg_from,ratiocfg_to,ratio_ion_from,ratio_ion_to,drift,kinetic,&
&kinetic_detail,wfn_detail,kinetic_aderiv,wfn_aderiv)
implicit none
integer,intent(inout),optional :: ratio1_from,ratio1_to,ratio2_from,ra&
&tio2_to,ratiocfg_from,ratiocfg_to,ratio_ion_from,ratio_ion_to,drift,k&
&inetic,kinetic_detail,wfn_detail,kinetic_aderiv,wfn_aderiv
if(present(ratio1_from))then
if(.not.present(ratio1_to))call errstop('SCRATCH_REQUEST','Missing arg&
&ument (ratio1_to).')
if(ratio1_from==0)then
xyzzyaaak1=xyzzyaaak1+1
ratio1_from=xyzzyaaak1
endif
xyzzyaaah1(iratio1_from,ratio1_from)=.true.
if(ratio1_to==0)then
xyzzyaaak1=xyzzyaaak1+1
ratio1_to=xyzzyaaak1
endif
xyzzyaaah1(iratio1_to,ratio1_to)=.true.
xyzzyaaai1(ratio1_from,ratio1_to)=.false.
xyzzyaaai1(ratio1_to,ratio1_from)=.false.
elseif(present(ratio1_to))then
call errstop('SCRATCH_REQUEST','Missing argument (ratio1_from).')
endif
if(present(ratio2_from))then
if(.not.present(ratio2_to))call errstop('SCRATCH_REQUEST','Missing arg&
&ument (ratio2_to).')
if(ratio2_from==0)then
xyzzyaaak1=xyzzyaaak1+1
ratio2_from=xyzzyaaak1
endif
xyzzyaaah1(iratio2_from,ratio2_from)=.true.
if(ratio2_to==0)then
xyzzyaaak1=xyzzyaaak1+1
ratio2_to=xyzzyaaak1
endif
xyzzyaaah1(iratio2_to,ratio2_to)=.true.
xyzzyaaai1(ratio2_from,ratio2_to)=.false.
xyzzyaaai1(ratio2_to,ratio2_from)=.false.
elseif(present(ratio2_to))then
call errstop('SCRATCH_REQUEST','Missing argument (ratio2_from).')
endif
if(present(ratiocfg_from))then
if(.not.present(ratiocfg_to))call errstop('SCRATCH_REQUEST','Missing a&
&rgument (ratiocfg_to).')
if(ratiocfg_from==0)then
xyzzyaaak1=xyzzyaaak1+1
ratiocfg_from=xyzzyaaak1
endif
xyzzyaaah1(iratiocfg_from,ratiocfg_from)=.true.
if(ratiocfg_to==0)then
xyzzyaaak1=xyzzyaaak1+1
ratiocfg_to=xyzzyaaak1
endif
xyzzyaaah1(iratiocfg_to,ratiocfg_to)=.true.
xyzzyaaai1(ratiocfg_from,ratiocfg_to)=.false.
xyzzyaaai1(ratiocfg_to,ratiocfg_from)=.false.
elseif(present(ratiocfg_to))then
call errstop('SCRATCH_REQUEST','Missing argument (ratiocfg_from).')
endif
if(present(ratio_ion_from))then
if(.not.present(ratio_ion_to))call errstop('SCRATCH_REQUEST','Missing &
&argument (ratio_ion_to).')
if(ratio_ion_from==0)then
xyzzyaaak1=xyzzyaaak1+1
ratio_ion_from=xyzzyaaak1
endif
xyzzyaaah1(iratio_ion_from,ratio_ion_from)=.true.
if(ratio_ion_to==0)then
xyzzyaaak1=xyzzyaaak1+1
ratio_ion_to=xyzzyaaak1
endif
xyzzyaaah1(iratio_ion_to,ratio_ion_to)=.true.
xyzzyaaai1(ratio_ion_from,ratio_ion_to)=.false.
xyzzyaaai1(ratio_ion_to,ratio_ion_from)=.false.
elseif(present(ratio_ion_to))then
call errstop('SCRATCH_REQUEST','Missing argument (ratio_ion_from).')
endif
if(present(drift))then
if(drift==0)then
xyzzyaaak1=xyzzyaaak1+1
drift=xyzzyaaak1
endif
xyzzyaaah1(idrift,drift)=.true.
endif
if(present(kinetic))then
if(kinetic==0)then
xyzzyaaak1=xyzzyaaak1+1
kinetic=xyzzyaaak1
endif
xyzzyaaah1(ikinetic,kinetic)=.true.
endif
if(present(kinetic_detail))then
if(kinetic_detail==0)then
xyzzyaaak1=xyzzyaaak1+1
kinetic_detail=xyzzyaaak1
endif
xyzzyaaah1(ikinetic_detail,kinetic_detail)=.true.
endif
if(present(wfn_detail))then
if(wfn_detail==0)then
xyzzyaaak1=xyzzyaaak1+1
wfn_detail=xyzzyaaak1
endif
xyzzyaaah1(iwfn_detail,wfn_detail)=.true.
endif
if(present(kinetic_aderiv))then
if(kinetic_aderiv==0)then
xyzzyaaak1=xyzzyaaak1+1
kinetic_aderiv=xyzzyaaak1
endif
xyzzyaaah1(ikinetic_aderiv,kinetic_aderiv)=.true.
endif
if(present(wfn_aderiv))then
if(wfn_aderiv==0)then
xyzzyaaak1=xyzzyaaak1+1
wfn_aderiv=xyzzyaaak1
endif
xyzzyaaah1(iwfn_aderiv,wfn_aderiv)=.true.
endif
end subroutine scratch_request
subroutine scratch_protect(is1,is2)
implicit none
integer,intent(inout) :: is1
integer,intent(inout),optional :: is2
if(is1==0)then
xyzzyaaak1=xyzzyaaak1+1
is1=xyzzyaaak1
endif
if(present(is2))then
if(is2==0)then
xyzzyaaak1=xyzzyaaak1+1
is2=xyzzyaaak1
endif
if(is1==is2)call errstop('SCRATCH_PROTECT','Cannot prevent a buffer fr&
&om coinciding with itself.')
xyzzyaaai1(is1,is2)=.false.
xyzzyaaai1(is2,is1)=.false.
else
xyzzyaaai1(is1,:)=.false.
xyzzyaaai1(:,is1)=.false.
endif
end subroutine scratch_protect
subroutine setup_scratch
use parallel, only : am_master
implicit none
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4,xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4,xyzzyaaak4,xyzzyaaal4&
&,xyzzyaaam4(2,xyzzyaaaf1),xyzzyaaan4,xyzzyaaao4
integer xyzzyaaap4(2),xyzzyaaaq4(2),xyzzyaaar4(2)
integer,parameter :: xyzzyaaas4=1000
integer,allocatable :: xyzzyaaat4(:),xyzzyaaau4(:),xyzzyaaav4(:)
logical xyzzyaaaw4,xyzzyaaax4,xyzzyaaay4
logical,parameter :: xyzzyaaaz4=.false.
logical,allocatable :: xyzzyaaba4(:,:)
character(1) char_1
character(80) char_80,tmpr
allocate(xyzzyaaat4(xyzzyaaak1))
do xyzzyaaab4=1,xyzzyaaak1
xyzzyaaat4(xyzzyaaab4)=xyzzyaaab4
enddo
do xyzzyaaag4=xyzzyaaaf1,0,-1
do xyzzyaaab4=1,xyzzyaaak1
if(xyzzyaaat4(xyzzyaaab4)/=xyzzyaaab4)cycle
do
xyzzyaaac4=0
xyzzyaaaf4=-1
do xyzzyaaad4=xyzzyaaab4+1,xyzzyaaak1
if(xyzzyaaat4(xyzzyaaad4)/=xyzzyaaad4.or..not.xyzzyaaai1(xyzzyaaab4,xy&
&zzyaaad4))cycle
xyzzyaaae4=count(xyzzyaaah1(:,xyzzyaaab4).and.xyzzyaaah1(:,xyzzyaaad4)&
&)
if(xyzzyaaae4>xyzzyaaaf4)then
xyzzyaaaf4=xyzzyaaae4
xyzzyaaac4=xyzzyaaad4
endif
enddo
if(xyzzyaaac4==0.or.xyzzyaaaf4<xyzzyaaag4)exit
where(xyzzyaaat4(:)==xyzzyaaac4)xyzzyaaat4(:)=xyzzyaaab4
xyzzyaaah1(:,xyzzyaaab4)=xyzzyaaah1(:,xyzzyaaab4).or.xyzzyaaah1(:,xyzz&
&yaaac4)
xyzzyaaai1(:,xyzzyaaab4)=xyzzyaaai1(:,xyzzyaaab4).and.xyzzyaaai1(:,xyz&
&zyaaac4)
xyzzyaaai1(xyzzyaaab4,:)=xyzzyaaai1(:,xyzzyaaab4)
enddo
enddo
enddo
allocate(xyzzyaaau4(xyzzyaaak1),xyzzyaaba4(xyzzyaaaf1,xyzzyaaak1))
xyzzyaaau4=0
nscratch=0
xyzzyaaba4=.false.
do xyzzyaaab4=1,xyzzyaaak1
if(xyzzyaaat4(xyzzyaaab4)/=xyzzyaaab4)cycle
nscratch=nscratch+1
where(xyzzyaaat4(:)==xyzzyaaab4)xyzzyaaau4(:)=nscratch
xyzzyaaba4(:,nscratch)=xyzzyaaah1(:,xyzzyaaab4)
enddo
if(nscratch<1)call errstop('SETUP_SCRATCH','No scratch spaces!')
deallocate(xyzzyaaat4)
xyzzyaaao4=0
allocate(xyzzyaaat4(nscratch),xyzzyaaav4(nscratch))
do xyzzyaaab4=1,nscratch
xyzzyaaat4(xyzzyaaab4)=xyzzyaaab4
enddo
xyzzyaaak4=xyzzyaaaf1*nscratch
xyzzyaaav4=xyzzyaaat4
xyzzyaaaw4=.true.
do
xyzzyaaaj4=0
do xyzzyaaal4=1,xyzzyaaaf1
xyzzyaaay4=.true.
xyzzyaaan4=0
do xyzzyaaab4=1,nscratch
xyzzyaaac4=minloc(xyzzyaaat4,1,xyzzyaaat4==xyzzyaaab4)
if(xyzzyaaba4(xyzzyaaal4,xyzzyaaac4))then
xyzzyaaay4=.false.
xyzzyaaaj4=xyzzyaaaj4+xyzzyaaan4
xyzzyaaan4=0
elseif(.not.xyzzyaaay4)then
xyzzyaaan4=xyzzyaaan4+1
endif
enddo
enddo
if(xyzzyaaaj4<xyzzyaaak4)then
xyzzyaaak4=xyzzyaaaj4
xyzzyaaav4=xyzzyaaat4
endif
if(xyzzyaaaj4==0)exit
if(xyzzyaaao4>xyzzyaaas4)exit
do
call next_permutation(nscratch,xyzzyaaat4,xyzzyaaaw4,xyzzyaaax4)
if(xyzzyaaaw4.or.xyzzyaaax4)exit
enddo
xyzzyaaao4=xyzzyaaao4+1
if(xyzzyaaax4)exit
enddo
deallocate(xyzzyaaat4)
allocate(xyzzyaaaj1(xyzzyaaak1),scr_tasks(xyzzyaaaf1,nscratch))
do xyzzyaaab4=1,xyzzyaaak1
xyzzyaaaj1(xyzzyaaab4)=xyzzyaaav4(xyzzyaaau4(xyzzyaaab4))
scr_tasks(:,xyzzyaaaj1(xyzzyaaab4))=xyzzyaaba4(:,xyzzyaaau4(xyzzyaaab4&
&))
enddo
deallocate(xyzzyaaau4,xyzzyaaav4,xyzzyaaba4)
if(am_master.and.xyzzyaaaz4)then
char_80=''
do xyzzyaaab4=1,nscratch
char_80=trim(char_80)//repeat(' ',5)//trim(i2s(xyzzyaaab4))
enddo
call wout('Scratch spaces layout            '//trim(adjustl(char_80)))
call wout(repeat('-',28+6*nscratch))
char_80=''
do xyzzyaaab4=1,nscratch
char_80=trim(char_80)//repeat(' ',5)//'X'
enddo
call wout('BASIC                            '//trim(adjustl(char_80)))
do xyzzyaaal4=1,xyzzyaaaf1
char_80=''
do xyzzyaaab4=1,nscratch
if(scr_tasks(xyzzyaaal4,xyzzyaaab4))then
char_1='X'
else
char_1='-'
endif
char_80=trim(char_80)//repeat(' ',5)//char_1
enddo
write(tmpr,'(a,t34,a)')trim(task_name(xyzzyaaal4)),trim(adjustl(char_8&
&0))
call wout(tmpr)
enddo
call wout()
endif
xyzzyaaam4=0
do xyzzyaaal4=1,xyzzyaaaf1
xyzzyaaah4=0
xyzzyaaai4=0
do xyzzyaaab4=1,nscratch
if(scr_tasks(xyzzyaaal4,xyzzyaaab4))then
if(xyzzyaaah4==0)xyzzyaaah4=xyzzyaaab4
xyzzyaaai4=xyzzyaaab4
endif
enddo
xyzzyaaam4(1:2,xyzzyaaal4)=(/xyzzyaaah4,xyzzyaaai4/)
enddo
ratio1_from_sz(1:2)=xyzzyaaam4(1:2,iratio1_from)
ratio1_to_sz(1:2)=xyzzyaaam4(1:2,iratio1_to)
ratio2_from_sz(1:2)=xyzzyaaam4(1:2,iratio2_from)
ratio2_to_sz(1:2)=xyzzyaaam4(1:2,iratio2_to)
ratiocfg_from_sz(1:2)=xyzzyaaam4(1:2,iratiocfg_from)
ratiocfg_to_sz(1:2)=xyzzyaaam4(1:2,iratiocfg_to)
ratio_ion_from_sz(1:2)=xyzzyaaam4(1:2,iratio_ion_from)
ratio_ion_to_sz(1:2)=xyzzyaaam4(1:2,iratio_ion_to)
drift_sz(1:2)=xyzzyaaam4(1:2,idrift)
kinetic_sz(1:2)=xyzzyaaam4(1:2,ikinetic)
kinetic_detail_sz(1:2)=xyzzyaaam4(1:2,ikinetic_detail)
wfn_detail_sz(1:2)=xyzzyaaam4(1:2,iwfn_detail)
kinetic_aderiv_sz(1:2)=xyzzyaaam4(1:2,ikinetic_aderiv)
wfn_aderiv_sz(1:2)=xyzzyaaam4(1:2,iwfn_aderiv)
allocate(buffer_move1_from(nscratch),buffer_move1_from_ii(nscratch),bu&
&ffer_move2_from(nscratch),buffer_move2_from_ii(nscratch),buffer_move2&
&_from_jj(nscratch),buffer_move_ion_from(nscratch),buffer_move_ion_fro&
&m_ion(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','buffer_move*')
buffer_move1_from=0
buffer_move1_from_ii=0
buffer_move2_from=0
buffer_move2_from_ii=0
buffer_move2_from_jj=0
buffer_move_ion_from=0
buffer_move_ion_from_ion=0
xyzzyaaap4=0
xyzzyaaaq4=0
xyzzyaaar4=0
call include_range(ratio1_to_sz,xyzzyaaap4)
call include_range(ratio2_to_sz,xyzzyaaaq4)
call include_range(ratio_ion_to_sz,xyzzyaaar4)
allocate(rele_scr(3,netot,nscratch),sele_scr(netot,nscratch),eevecs_sc&
&r(4,netot,netot,nscratch),eivecs_scr(4,max(nitot,1),netot,nscratch),s&
&tat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','cfg')
allocate(rsele_valid(nscratch),xyzzyaaac1(nscratch),eevecs_valid(nscra&
&tch),xyzzyaaad1(nscratch),eivecs_valid(nscratch),xyzzyaaae1(nscratch)&
&,stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','cfg_valid')
rele_scr=0.d0
sele_scr=0
eevecs_scr=0.d0
eivecs_scr=0.d0
if(xyzzyaaap4(1)/=0)then
allocate(rele1_chscr(3,xyzzyaaap4(1):xyzzyaaap4(2)),sele1_chscr(xyzzya&
&aap4(1):xyzzyaaap4(2)),eevecs1_chscr(4,netot,xyzzyaaap4(1):xyzzyaaap4&
&(2)),eivecs1_chscr(4,max(nitot,1),xyzzyaaap4(1):xyzzyaaap4(2)),stat=x&
&yzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','vecs1')
rele1_chscr=0.d0
sele1_chscr=0
eevecs1_chscr=0.d0
eivecs1_chscr=0.d0
endif
allocate(rsele1_chvalid(nscratch),eevecs1_chvalid(nscratch),eivecs1_ch&
&valid(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','vecs1_valid')
if(xyzzyaaaq4(1)/=0)then
allocate(rele2ii_chscr(3,xyzzyaaaq4(1):xyzzyaaaq4(2)),xyzzyaaaa1(xyzzy&
&aaaq4(1):xyzzyaaaq4(2)),rele2jj_chscr(3,xyzzyaaaq4(1):xyzzyaaaq4(2)),&
&xyzzyaaab1(xyzzyaaaq4(1):xyzzyaaaq4(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','vecs2')
rele2ii_chscr=0.d0
xyzzyaaaa1=0
rele2jj_chscr=0.d0
xyzzyaaab1=0
endif
allocate(rsele2_chvalid(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','vecs2_valid')
if(xyzzyaaar4(1)/=0)then
allocate(rion1_chscr(3,xyzzyaaar4(1):xyzzyaaar4(2)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','rion1_chscr')
rion1_chscr=0.d0
endif
allocate(rion1_chvalid(nscratch),stat=xyzzyaaaa4)
allocate(xyzzyaaao1(nscratch),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','rion_HACK')
xyzzyaaao1=0
allocate(eevecs1(4,netot),eivecs1(4,max(nitot,1)),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'SETUP_SCRATCH','eevecs1')
eevecs1=0.d0
eivecs1=0.d0
end subroutine setup_scratch
subroutine which_scratch(is)
implicit none
integer,intent(inout) :: is
if(is>0)is=xyzzyaaaj1(is)
end subroutine which_scratch
subroutine include_range(x,y)
implicit none
integer,intent(in) :: x(2)
integer,intent(inout) :: y(2)
if(x(1)==0)return
if(y(1)==0)then
y=x
else
y(1)=min(x(1),y(1))
y(2)=max(x(2),y(2))
endif
end subroutine include_range
logical function in_range(x,y)
implicit none
integer,intent(in) :: x,y(2)
in_range=x>=y(1).and.x<=y(2).and.x/=0
end function in_range
subroutine finish_scratch
implicit none
deallocate(buffer_move1_from,buffer_move1_from_ii,buffer_move2_from,bu&
&ffer_move2_from_ii,buffer_move2_from_jj,buffer_move_ion_from,buffer_m&
&ove_ion_from_ion)
deallocate(rele_scr,sele_scr,eevecs_scr,eivecs_scr)
deallocate(rsele_valid,xyzzyaaac1,eevecs_valid,xyzzyaaad1,eivecs_valid&
&,xyzzyaaae1)
if(allocated(eevecs1_chscr))deallocate(rele1_chscr,sele1_chscr,eevecs1&
&_chscr,eivecs1_chscr)
deallocate(rsele1_chvalid,eevecs1_chvalid,eivecs1_chvalid)
if(allocated(rele2ii_chscr))deallocate(rele2ii_chscr,xyzzyaaaa1,rele2j&
&j_chscr,xyzzyaaab1)
deallocate(rsele2_chvalid)
if(allocated(rion1_chscr))deallocate(rion1_chscr)
deallocate(rion1_chvalid)
deallocate(xyzzyaaao1)
deallocate(eevecs1,eivecs1)
xyzzyaaah1=.false.
xyzzyaaai1=.true.
deallocate(scr_tasks,xyzzyaaaj1)
xyzzyaaak1=0
end subroutine finish_scratch
subroutine define_config_scratch(is,rele,sele)
implicit none
integer,intent(in) :: is,sele(netot)
real(dp),intent(in) :: rele(3,netot)
call dcopy(three_netot,rele(1,1),1,rele_scr(1,1,is),1)
sele_scr(:,is)=sele
rsele_valid(is)=.true.
xyzzyaaac1(is)=0
xyzzyaaao1(is)=0
end subroutine define_config_scratch
subroutine define_config_oneelec_scratch(ii,is,js,rnew,snew)
implicit none
integer,intent(in) :: ii,is,js,snew
real(dp),intent(in) :: rnew(3)
buffer_move1_from(js)=is
buffer_move1_from_ii(js)=ii
rele1_chscr(:,js)=rnew
sele1_chscr(js)=snew
rsele1_chvalid(js)=.true.
if(xyzzyaaac1(js)==is)then
call dcopy(3,rnew(1),1,rele_scr(1,ii,js),1)
sele_scr(ii,js)=snew
rsele_valid(js)=.true.
xyzzyaaac1(js)=0
else
rsele_valid(js)=.false.
xyzzyaaac1(js)=0
endif
xyzzyaaao1(js)=xyzzyaaao1(is)
end subroutine define_config_oneelec_scratch
subroutine define_config_twoelec_scratch(ii,jj,is,js,riinew,siinew,rjj&
&new,sjjnew)
implicit none
integer,intent(in) :: ii,jj,is,js,siinew,sjjnew
real(dp),intent(in) :: riinew(3),rjjnew(3)
buffer_move2_from(js)=is
buffer_move2_from_ii(js)=ii
buffer_move2_from_jj(js)=jj
rele2ii_chscr(:,js)=riinew
xyzzyaaaa1(js)=siinew
rele2jj_chscr(:,js)=rjjnew
xyzzyaaab1(js)=sjjnew
rsele2_chvalid(js)=.true.
if(xyzzyaaac1(js)==is)then
call dcopy(3,riinew(1),1,rele_scr(1,ii,js),1)
call dcopy(3,rjjnew(1),1,rele_scr(1,jj,js),1)
sele_scr(ii,js)=siinew
sele_scr(jj,js)=sjjnew
rsele_valid(js)=.true.
xyzzyaaac1(js)=0
else
rsele_valid(js)=.false.
xyzzyaaac1(js)=0
endif
xyzzyaaao1(js)=xyzzyaaao1(is)
end subroutine define_config_twoelec_scratch
subroutine define_config_oneion_scratch(ion,is,js,rion1_new)
implicit none
integer,intent(in) :: ion,is,js
real(dp),intent(in) :: rion1_new(3)
buffer_move_ion_from(js)=is
buffer_move_ion_from_ion(js)=ion
call dcopy(three_netot,rele_scr(1,1,is),1,rele_scr(1,1,js),1)
sele_scr(:,js)=sele_scr(:,is)
rsele_valid(js)=.true.
xyzzyaaac1(js)=0
rion1_chscr(:,js)=rion1_new(:)
rion1_chvalid(js)=.true.
xyzzyaaao1(js)=js
end subroutine define_config_oneion_scratch
subroutine accept_move_scratch(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa13,xyzzyaaab13
if(buffer_move1_from(js)==is)then
if(rsele_valid(is))then
xyzzyaaaa13=buffer_move1_from_ii(js)
rele_scr(:,xyzzyaaaa13,is)=rele1_chscr(:,js)
sele_scr(xyzzyaaaa13,is)=sele1_chscr(js)
rsele_valid(is)=.true.
xyzzyaaac1(is)=0
elseif(rsele_valid(js))then
call dcopy(three_netot,rele_scr(1,1,js),1,rele_scr(1,1,is),1)
sele_scr(:,is)=sele_scr(:,js)
rsele_valid(is)=.true.
xyzzyaaac1(is)=0
else
rsele_valid(is)=.false.
xyzzyaaac1(is)=0
endif
elseif(buffer_move2_from(js)==is)then
if(rsele_valid(is))then
xyzzyaaaa13=buffer_move2_from_ii(js)
xyzzyaaab13=buffer_move2_from_jj(js)
rele_scr(:,xyzzyaaaa13,is)=rele2ii_chscr(:,js)
rele_scr(:,xyzzyaaab13,is)=rele2jj_chscr(:,js)
sele_scr(xyzzyaaaa13,is)=xyzzyaaaa1(js)
sele_scr(xyzzyaaab13,is)=xyzzyaaab1(js)
rsele_valid(is)=.true.
xyzzyaaac1(is)=0
elseif(rsele_valid(js))then
call dcopy(three_netot,rele_scr(1,1,js),1,rele_scr(1,1,is),1)
sele_scr(:,is)=sele_scr(:,js)
rsele_valid(is)=.true.
xyzzyaaac1(is)=0
endif
else
if(rsele_valid(js))then
call dcopy(three_netot,rele_scr(1,1,js),1,rele_scr(1,1,is),1)
sele_scr(:,is)=sele_scr(:,js)
rsele_valid(is)=.true.
xyzzyaaac1(is)=0
else
rsele_valid(is)=.false.
xyzzyaaac1(is)=0
endif
endif
if(buffer_move1_from(js)==is)then
if(eevecs_valid(is).and.eevecs1_chvalid(js))then
xyzzyaaaa13=buffer_move1_from_ii(js)
call dcopy_ee(netot,xyzzyaaaa13,eevecs1_chscr(1,1,js),eevecs_scr(1,1,1&
&,is))
eevecs_valid(is)=.true.
xyzzyaaad1(is)=0
elseif(eevecs_valid(js))then
call dcopy(four_netot_netot,eevecs_scr(1,1,1,js),1,eevecs_scr(1,1,1,is&
&),1)
eevecs_valid(is)=.true.
xyzzyaaad1(is)=0
else
eevecs_valid(is)=.false.
xyzzyaaad1(is)=0
endif
elseif(buffer_move2_from(js)==is)then
if(eevecs_valid(js))then
call dcopy(four_netot_netot,eevecs_scr(1,1,1,js),1,eevecs_scr(1,1,1,is&
&),1)
eevecs_valid(is)=.true.
xyzzyaaad1(is)=0
else
eevecs_valid(is)=.false.
xyzzyaaad1(is)=0
endif
else
if(eevecs_valid(js))then
call dcopy(four_netot_netot,eevecs_scr(1,1,1,js),1,eevecs_scr(1,1,1,is&
&),1)
eevecs_valid(is)=.true.
xyzzyaaad1(is)=0
else
eevecs_valid(is)=.false.
xyzzyaaad1(is)=0
endif
endif
if(buffer_move1_from(js)==is)then
if(eivecs_valid(is).and.eivecs1_chvalid(js))then
xyzzyaaaa13=buffer_move1_from_ii(js)
call dcopy(four_nitot,eivecs1_chscr(1,1,js),1,eivecs_scr(1,1,xyzzyaaaa&
&13,is),1)
eivecs_valid(is)=.true.
xyzzyaaae1(is)=0
elseif(eivecs_valid(js))then
call dcopy(four_nitot_netot,eivecs_scr(1,1,1,js),1,eivecs_scr(1,1,1,is&
&),1)
eivecs_valid(is)=.true.
xyzzyaaae1(is)=0
else
eivecs_valid(is)=.false.
xyzzyaaae1(is)=0
endif
elseif(buffer_move2_from(js)==is)then
if(eivecs_valid(js))then
call dcopy(four_nitot_netot,eivecs_scr(1,1,1,js),1,eivecs_scr(1,1,1,is&
&),1)
eivecs_valid(is)=.true.
xyzzyaaae1(is)=0
else
eivecs_valid(is)=.false.
xyzzyaaae1(is)=0
endif
else
if(eivecs_valid(js))then
call dcopy(four_nitot_netot,eivecs_scr(1,1,1,js),1,eivecs_scr(1,1,1,is&
&),1)
eivecs_valid(is)=.true.
xyzzyaaae1(is)=0
else
eivecs_valid(is)=.false.
xyzzyaaae1(is)=0
endif
endif
eevecs1_chvalid(is)=.false.
eivecs1_chvalid(is)=.false.
rsele1_chvalid(is)=.false.
rsele2_chvalid(is)=.false.
rion1_chvalid(is)=.false.
if(buffer_move_ion_from(js)==is)then
xyzzyaaao1(is)=0
else
xyzzyaaao1(is)=xyzzyaaao1(js)
endif
end subroutine accept_move_scratch
subroutine reset_config_scratch(is,js)
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14
if(xyzzyaaac1(js)==is)then
if(buffer_move_ion_from(js)==is)then
rsele_valid(js)=.true.
xyzzyaaac1(js)=0
else
rsele_valid(js)=.false.
xyzzyaaac1(js)=is
endif
elseif(buffer_move1_from(js)==is)then
if(rsele_valid(is).and.rsele_valid(js))then
xyzzyaaaa14=buffer_move1_from_ii(js)
rele_scr(1:3,xyzzyaaaa14,js)=rele_scr(1:3,xyzzyaaaa14,is)
sele_scr(xyzzyaaaa14,js)=sele_scr(xyzzyaaaa14,is)
rsele_valid(js)=.false.
xyzzyaaac1(js)=is
else
rsele_valid(js)=.false.
xyzzyaaac1(js)=0
endif
elseif(buffer_move2_from(js)==is)then
if(rsele_valid(is).and.rsele_valid(js))then
xyzzyaaaa14=buffer_move2_from_ii(js)
xyzzyaaab14=buffer_move2_from_jj(js)
rele_scr(1:3,xyzzyaaaa14,js)=rele_scr(1:3,xyzzyaaaa14,is)
rele_scr(1:3,xyzzyaaab14,js)=rele_scr(1:3,xyzzyaaab14,is)
sele_scr(xyzzyaaaa14,js)=sele_scr(xyzzyaaaa14,is)
sele_scr(xyzzyaaab14,js)=sele_scr(xyzzyaaab14,is)
rsele_valid(js)=.false.
xyzzyaaac1(js)=is
else
rsele_valid(js)=.false.
xyzzyaaac1(js)=0
endif
elseif(buffer_move_ion_from(js)==is)then
if(rsele_valid(js))then
rsele_valid(js)=.true.
xyzzyaaac1(js)=0
else
rsele_valid(js)=.false.
xyzzyaaac1(js)=0
endif
else
rsele_valid(js)=.false.
xyzzyaaac1(js)=0
endif
if(xyzzyaaad1(js)==is)then
if(buffer_move_ion_from(js)==is)then
eevecs_valid(js)=.true.
xyzzyaaad1(js)=0
else
eevecs_valid(js)=.false.
xyzzyaaad1(js)=is
endif
elseif(buffer_move1_from(js)==is)then
if(eevecs_valid(is).and.eevecs_valid(js))then
xyzzyaaaa14=buffer_move1_from_ii(js)
call dcopy_ee(netot,xyzzyaaaa14,eevecs_scr(1,1,xyzzyaaaa14,is),eevecs_&
&scr(1,1,1,js))
eevecs_valid(js)=.false.
xyzzyaaad1(js)=is
else
eevecs_valid(js)=.false.
xyzzyaaad1(js)=0
endif
elseif(buffer_move2_from(js)==is)then
if(eevecs_valid(is).and.eevecs_valid(js))then
xyzzyaaaa14=buffer_move2_from_ii(js)
xyzzyaaab14=buffer_move2_from_jj(js)
call dcopy_ee(netot,xyzzyaaaa14,eevecs_scr(1,1,xyzzyaaaa14,is),eevecs_&
&scr(1,1,1,js))
call dcopy_ee(netot,xyzzyaaab14,eevecs_scr(1,1,xyzzyaaab14,is),eevecs_&
&scr(1,1,1,js))
eevecs_valid(js)=.false.
xyzzyaaad1(js)=is
else
eevecs_valid(js)=.false.
xyzzyaaad1(js)=0
endif
elseif(buffer_move_ion_from(js)==is)then
if(eevecs_valid(js))then
eevecs_valid(js)=.true.
xyzzyaaad1(js)=0
else
eevecs_valid(js)=.false.
xyzzyaaad1(js)=0
endif
else
eevecs_valid(js)=.false.
xyzzyaaad1(js)=0
endif
if(xyzzyaaae1(js)==is)then
eivecs_valid(js)=.false.
xyzzyaaae1(js)=is
elseif(buffer_move1_from(js)==is)then
if(eivecs_valid(is).and.eivecs_valid(js))then
xyzzyaaaa14=buffer_move1_from_ii(js)
call dcopy(four_nitot,eivecs_scr(1,1,xyzzyaaaa14,is),1,eivecs_scr(1,1,&
&xyzzyaaaa14,js),1)
eivecs_valid(js)=.false.
xyzzyaaae1(js)=is
else
eivecs_valid(js)=.false.
xyzzyaaae1(js)=0
endif
elseif(buffer_move2_from(js)==is)then
if(eivecs_valid(is).and.eivecs_valid(js))then
xyzzyaaaa14=buffer_move2_from_ii(js)
xyzzyaaab14=buffer_move2_from_jj(js)
call dcopy(four_nitot,eivecs_scr(1,1,xyzzyaaaa14,is),1,eivecs_scr(1,1,&
&xyzzyaaaa14,js),1)
call dcopy(four_nitot,eivecs_scr(1,1,xyzzyaaab14,is),1,eivecs_scr(1,1,&
&xyzzyaaab14,js),1)
eivecs_valid(js)=.false.
xyzzyaaae1(js)=is
else
eivecs_valid(js)=.false.
xyzzyaaae1(js)=0
endif
elseif(buffer_move_ion_from(js)==is)then
if(eivecs_valid(is).and.eivecs_valid(js))then
xyzzyaaac14=buffer_move_ion_from_ion(js)
do xyzzyaaaa14=1,netot
call dcopy(4,eivecs_scr(1,xyzzyaaac14,xyzzyaaaa14,is),1,eivecs_scr(1,x&
&yzzyaaac14,xyzzyaaaa14,js),1)
enddo
eivecs_valid(js)=.false.
xyzzyaaae1(js)=is
else
eivecs_valid(js)=.false.
xyzzyaaae1(js)=0
endif
else
eivecs_valid(js)=.false.
xyzzyaaae1(js)=0
endif
eevecs1_chvalid(js)=.false.
eivecs1_chvalid(js)=.false.
rsele1_chvalid(js)=.false.
rsele2_chvalid(js)=.false.
rion1_chvalid(js)=.false.
xyzzyaaao1(js)=0
end subroutine reset_config_scratch
subroutine clear_scratch(is)
implicit none
integer,intent(in) :: is
rsele_valid(is)=.false.
xyzzyaaac1(is)=0
rsele1_chvalid(is)=.false.
rsele2_chvalid(is)=.false.
rion1_chvalid(is)=.false.
eevecs_valid(is)=.false.
xyzzyaaad1(is)=0
eevecs1_chvalid(is)=.false.
eivecs_valid(is)=.false.
xyzzyaaae1(is)=0
eivecs1_chvalid(is)=.false.
buffer_move1_from(is)=0
buffer_move1_from_ii(is)=0
buffer_move2_from(is)=0
buffer_move2_from_ii(is)=0
buffer_move2_from_jj(is)=0
buffer_move_ion_from(is)=0
buffer_move_ion_from_ion(is)=0
xyzzyaaao1(is)=0
end subroutine clear_scratch
subroutine push_rion_hack(is)
implicit none
integer,intent(in) :: is
integer xyzzyaaaa16
xyzzyaaaa16=xyzzyaaao1(is)
if(xyzzyaaaa16==0)then
which_ion_displaced=0
else
which_ion_displaced=buffer_move_ion_from_ion(xyzzyaaaa16)
ion_displacement(:)=rion1_chscr(:,xyzzyaaaa16)-rion(:,which_ion_displa&
&ced)
endif
end subroutine push_rion_hack
subroutine gen_config_geom(pt_config)
implicit none
type(config_geom),pointer :: pt_config
integer xyzzyaaaa17
allocate(pt_config,stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'GEN_CONFIG_GEOM','container')
allocate(pt_config%pt_rele(3,netot),stat=xyzzyaaaa17)
call check_alloc(xyzzyaaaa17,'GEN_CONFIG_GEOM','rele')
end subroutine gen_config_geom
subroutine delete_config_geom(pt_config)
implicit none
type(config_geom),pointer :: pt_config
deallocate(pt_config%pt_rele)
deallocate(pt_config)
end subroutine delete_config_geom
subroutine copy_config_geom(pt_from,pt_to)
implicit none
type(config_geom),pointer :: pt_from,pt_to
call dcopy(three_netot,pt_from%pt_rele,1,pt_to%pt_rele,1)
end subroutine copy_config_geom
subroutine config_to_pt_geom(pt_config,k)
use slaarnaaf, only : rele_config
implicit none
integer,intent(in) :: k
type(config_geom),pointer :: pt_config
call dcopy(three_netot,rele_config(1,1,k),1,pt_config%pt_rele,1)
end subroutine config_to_pt_geom
subroutine pt_to_config_geom(pt_config)
use slaarnaaf,only : add_config
implicit none
type(config_geom),pointer :: pt_config
call add_config(modify=.true.,rele=pt_config%pt_rele)
end subroutine pt_to_config_geom
subroutine redist_allocations_geom(kmax)
implicit none
integer,intent(in) :: kmax
integer xyzzyaaaa22
allocate(xyzzyaaal1(3,netot,kmax),stat=xyzzyaaaa22)
call check_alloc(xyzzyaaaa22,'REDIST_ALLOCATIONS_GEOM','')
end subroutine redist_allocations_geom
subroutine redist_load_geom(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_geom),pointer :: pt_config
call dcopy(three_netot,pt_config%pt_rele,1,xyzzyaaal1(1,1,k),1)
end subroutine redist_load_geom
subroutine redist_send_geom(jnode,kbase,k,nbt,reqbase,blocking)
use parallel
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa24
if(blocking)then
call mpi_ssend(xyzzyaaal1(1,1,kbase+1),k*three_netot,mpi_double_precis&
&ion,jnode,move_msg,mpi_comm_world,ierror)
else
nbt=nbt+1
call mpi_isend(xyzzyaaal1(1,1,kbase+1),k*three_netot,mpi_double_precis&
&ion,jnode,nbt,mpi_comm_world,xyzzyaaaa24,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa24
endif
call checkmpi(ierror,'ssend rele in redist_send')
end subroutine redist_send_geom
subroutine redist_recv_geom(jnode,kbase,k,nbt,reqbase,blocking)
use parallel
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa25
if(blocking)then
call mpi_recv(xyzzyaaal1(1,1,kbase+1),k*three_netot,mpi_double_precisi&
&on,jnode,move_msg,mpi_comm_world,status,ierror)
else
nbt=nbt+1
call mpi_irecv(xyzzyaaal1(1,1,kbase+1),k*three_netot,mpi_double_precis&
&ion,jnode,nbt,mpi_comm_world,xyzzyaaaa25,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa25
endif
call checkmpi(ierror,'recv rele in redist_recv')
end subroutine redist_recv_geom
subroutine redist_save_geom(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_geom),pointer :: pt_config
call dcopy(three_netot,xyzzyaaal1(1,1,k),1,pt_config%pt_rele,1)
end subroutine redist_save_geom
subroutine redist_deallocations_geom
implicit none
if(.not.allocated(xyzzyaaal1))call errstop('REDIST_DEALLOCATIONS_GEOM'&
&,'Trying to deallocate an unallocated array.')
deallocate(xyzzyaaal1)
end subroutine redist_deallocations_geom
subroutine load_from_pt_geom(is,pt_config,rele,sele)
implicit none
integer,intent(in) :: is
integer,intent(out) :: sele(netot)
real(dp),intent(out) :: rele(3,netot)
type(config_geom),pointer :: pt_config
call clear_scratch(is)
call dcopy(three_netot,pt_config%pt_rele,1,rele(1,1),1)
call dcopy(three_netot,rele(1,1),1,rele_scr(1,1,is),1)
sele=0
sele_scr(:,is)=0
rsele_valid(is)=.true.
xyzzyaaac1(is)=0
end subroutine load_from_pt_geom
subroutine save_to_pt_geom(is,pt_config)
implicit none
integer,intent(in) :: is
type(config_geom),pointer :: pt_config
call get_rsele(is)
call dcopy(three_netot,rele_scr(1,1,is),1,pt_config%pt_rele,1)
end subroutine save_to_pt_geom
subroutine clone_scratch_geom(is,js)
implicit none
integer,intent(in) :: is,js
call clear_scratch(js)
call get_rsele(is)
call dcopy(three_netot,rele_scr(1,1,is),1,rele_scr(1,1,js),1)
sele_scr(:,js)=sele_scr(:,is)
rsele_valid(js)=.true.
xyzzyaaac1(js)=0
if(eevecs_valid(is))then
call dcopy(four_netot_netot,eevecs_scr(1,1,1,is),1,eevecs_scr(1,1,1,js&
&),1)
eevecs_valid(js)=.true.
xyzzyaaad1(js)=0
else
eevecs_valid(js)=.false.
xyzzyaaad1(js)=0
endif
if(eivecs_valid(is))then
call dcopy(four_nitot_netot,eivecs_scr(1,1,1,is),1,eivecs_scr(1,1,1,js&
&),1)
eivecs_valid(js)=.true.
xyzzyaaae1(js)=0
else
eivecs_valid(js)=.false.
xyzzyaaae1(js)=0
endif
end subroutine clone_scratch_geom
recursive subroutine get_rsele(xyzzyaaaa31)
implicit none
integer,intent(in) :: xyzzyaaaa31
integer xyzzyaaab31,xyzzyaaac31,xyzzyaaad31,xyzzyaaae31
if(rsele_valid(xyzzyaaaa31))return
xyzzyaaad31=buffer_move1_from(xyzzyaaaa31)
xyzzyaaae31=buffer_move2_from(xyzzyaaaa31)
if(xyzzyaaad31/=0)then
if(.not.rsele1_chvalid(xyzzyaaaa31))call errstop('GET_RSELE','Ill-defi&
&ned config (no rsele1).')
if(xyzzyaaac1(xyzzyaaaa31)/=xyzzyaaad31)then
call get_rsele(xyzzyaaad31)
call dcopy(three_netot,rele_scr(1,1,xyzzyaaad31),1,rele_scr(1,1,xyzzya&
&aaa31),1)
sele_scr(:,xyzzyaaaa31)=sele_scr(:,xyzzyaaad31)
endif
xyzzyaaab31=buffer_move1_from_ii(xyzzyaaaa31)
rele_scr(:,xyzzyaaab31,xyzzyaaaa31)=rele1_chscr(:,xyzzyaaaa31)
sele_scr(xyzzyaaab31,xyzzyaaaa31)=sele1_chscr(xyzzyaaaa31)
rsele_valid(xyzzyaaaa31)=.true.
xyzzyaaac1(xyzzyaaaa31)=0
elseif(xyzzyaaae31/=0)then
if(.not.rsele2_chvalid(xyzzyaaaa31))call errstop('GET_RSELE','Ill-defi&
&ned config (no rsele1).')
if(xyzzyaaac1(xyzzyaaaa31)/=xyzzyaaae31)then
call get_rsele(xyzzyaaae31)
call dcopy(three_netot,rele_scr(1,1,xyzzyaaae31),1,rele_scr(1,1,xyzzya&
&aaa31),1)
sele_scr(:,xyzzyaaaa31)=sele_scr(:,xyzzyaaae31)
endif
xyzzyaaab31=buffer_move2_from_ii(xyzzyaaaa31)
xyzzyaaac31=buffer_move2_from_jj(xyzzyaaaa31)
rele_scr(:,xyzzyaaab31,xyzzyaaaa31)=rele2ii_chscr(:,xyzzyaaaa31)
rele_scr(:,xyzzyaaac31,xyzzyaaaa31)=rele2jj_chscr(:,xyzzyaaaa31)
sele_scr(xyzzyaaab31,xyzzyaaaa31)=xyzzyaaaa1(xyzzyaaaa31)
sele_scr(xyzzyaaac31,xyzzyaaaa31)=xyzzyaaab1(xyzzyaaaa31)
rsele_valid(xyzzyaaaa31)=.true.
xyzzyaaac1(xyzzyaaaa31)=0
else
call errstop('GET_RSELE','Ill-defined config (nowhere to get RSELE fro&
&m).')
endif
end subroutine get_rsele
subroutine get_eevecs1_ch(ii,is)
implicit none
integer,intent(in) :: ii,is
integer xyzzyaaaa32
if(eevecs1_chvalid(is))return
if(eevecs_valid(is))then
call dcopy(four_netot,eevecs_scr(1,1,ii,is),1,eevecs1_chscr(1,1,is),1)
else
xyzzyaaaa32=buffer_move1_from(is)
call get_rsele(xyzzyaaaa32)
call ee_distances(netot,rele1_chscr(1,is),rele_scr(1,1,xyzzyaaaa32),ee&
&vecs1_chscr(1,1,is))
endif
eevecs1_chvalid(is)=.true.
end subroutine get_eevecs1_ch
recursive subroutine get_eevecs(xyzzyaaaa33)
implicit none
integer,intent(in) :: xyzzyaaaa33
integer xyzzyaaab33,xyzzyaaac33,xyzzyaaad33,xyzzyaaae33
if(eevecs_valid(xyzzyaaaa33))return
xyzzyaaad33=buffer_move1_from(xyzzyaaaa33)
xyzzyaaae33=buffer_move2_from(xyzzyaaaa33)
if(xyzzyaaad33/=0)then
xyzzyaaab33=buffer_move1_from_ii(xyzzyaaaa33)
if(xyzzyaaad1(xyzzyaaaa33)/=xyzzyaaad33)then
call get_eevecs(xyzzyaaad33)
call dcopy(four_netot_netot,eevecs_scr(1,1,1,xyzzyaaad33),1,eevecs_scr&
&(1,1,1,xyzzyaaaa33),1)
endif
if(eevecs1_chvalid(xyzzyaaaa33))then
call dcopy_ee(netot,xyzzyaaab33,eevecs1_chscr(1,1,xyzzyaaaa33),eevecs_&
&scr(1,1,1,xyzzyaaaa33))
else
call get_rsele(xyzzyaaad33)
call ee_distances(netot,rele1_chscr(1,xyzzyaaaa33),rele_scr(1,1,xyzzya&
&aad33),eevecs_scr(1,1,xyzzyaaab33,xyzzyaaaa33))
call dcopy_ee_flip(netot,xyzzyaaab33,eevecs_scr(1,1,xyzzyaaab33,xyzzya&
&aaa33),eevecs_scr(1,1,1,xyzzyaaaa33))
endif
eevecs_valid(xyzzyaaaa33)=.true.
xyzzyaaad1(xyzzyaaaa33)=0
elseif(xyzzyaaae33/=0)then
xyzzyaaab33=buffer_move2_from_ii(xyzzyaaaa33)
xyzzyaaac33=buffer_move2_from_jj(xyzzyaaaa33)
if(xyzzyaaad1(xyzzyaaaa33)/=xyzzyaaae33)then
call get_eevecs(xyzzyaaae33)
call dcopy(four_netot_netot,eevecs_scr(1,1,1,xyzzyaaae33),1,eevecs_scr&
&(1,1,1,xyzzyaaaa33),1)
endif
call get_rsele(xyzzyaaae33)
call ee_distances(netot,rele2ii_chscr(1,xyzzyaaaa33),rele_scr(1,1,xyzz&
&yaaae33),eevecs_scr(1,1,xyzzyaaab33,xyzzyaaaa33))
call dcopy_ee_flip(netot,xyzzyaaab33,eevecs_scr(1,1,xyzzyaaab33,xyzzya&
&aaa33),eevecs_scr(1,1,1,xyzzyaaaa33))
call ee_distances(netot,rele2jj_chscr(1,xyzzyaaaa33),rele_scr(1,1,xyzz&
&yaaae33),eevecs_scr(1,1,xyzzyaaac33,xyzzyaaaa33))
call dcopy_ee_flip(netot,xyzzyaaac33,eevecs_scr(1,1,xyzzyaaac33,xyzzya&
&aaa33),eevecs_scr(1,1,1,xyzzyaaaa33))
call ee_distances(1,rele2ii_chscr(1,xyzzyaaaa33),rele2jj_chscr(1,xyzzy&
&aaaa33),eevecs_scr(1,xyzzyaaac33,xyzzyaaab33,xyzzyaaaa33))
eevecs_scr(1:3,xyzzyaaab33,xyzzyaaac33,xyzzyaaaa33)=-eevecs_scr(1:3,xy&
&zzyaaac33,xyzzyaaab33,xyzzyaaaa33)
eevecs_scr(4,xyzzyaaab33,xyzzyaaac33,xyzzyaaaa33)=eevecs_scr(4,xyzzyaa&
&ac33,xyzzyaaab33,xyzzyaaaa33)
eevecs_valid(xyzzyaaaa33)=.true.
xyzzyaaad1(xyzzyaaaa33)=0
else
call get_rsele(xyzzyaaaa33)
call ee_distances_all(netot,rele_scr(1,1,xyzzyaaaa33),eevecs_scr(1,1,1&
&,xyzzyaaaa33))
eevecs_valid(xyzzyaaaa33)=.true.
xyzzyaaad1(xyzzyaaaa33)=0
endif
end subroutine get_eevecs
subroutine get_eivecs1_ch(ii,is)
implicit none
integer,intent(in) :: ii,is
integer xyzzyaaaa34,xyzzyaaab34
real(dp) xyzzyaaac34(3)
if(homogeneous_system.or.eivecs1_chvalid(is))return
if(eivecs_valid(is))then
call dcopy(four_nitot,eivecs_scr(1,1,ii,is),1,eivecs1_chscr(1,1,is),1)
else
xyzzyaaab34=xyzzyaaao1(is)
if(xyzzyaaab34/=0)then
xyzzyaaaa34=buffer_move_ion_from_ion(xyzzyaaab34)
xyzzyaaac34(1:3)=rion(1:3,xyzzyaaaa34)
rion(1:3,xyzzyaaaa34)=rion1_chscr(1:3,xyzzyaaab34)
call ei_distances(rele1_chscr(1,is),nitot,rion,eivecs1_chscr(1,1,is))
rion(1:3,xyzzyaaaa34)=xyzzyaaac34(1:3)
else
call ei_distances(rele1_chscr(1,is),nitot,rion,eivecs1_chscr(1,1,is))
endif
endif
eivecs1_chvalid(is)=.true.
end subroutine get_eivecs1_ch
recursive subroutine get_eivecs(xyzzyaaaa35)
implicit none
integer,intent(in) :: xyzzyaaaa35
integer xyzzyaaab35,xyzzyaaac35,xyzzyaaad35,xyzzyaaae35,xyzzyaaaf35,xy&
&zzyaaag35
real(dp) xyzzyaaah35(3)
if(homogeneous_system.or.eivecs_valid(xyzzyaaaa35))return
xyzzyaaad35=buffer_move1_from(xyzzyaaaa35)
xyzzyaaae35=buffer_move2_from(xyzzyaaaa35)
if(xyzzyaaad35/=0)then
xyzzyaaab35=buffer_move1_from_ii(xyzzyaaaa35)
if(xyzzyaaae1(xyzzyaaaa35)/=xyzzyaaad35)then
call get_eivecs(xyzzyaaad35)
call dcopy(four_nitot_netot,eivecs_scr(1,1,1,xyzzyaaad35),1,eivecs_scr&
&(1,1,1,xyzzyaaaa35),1)
endif
call get_eivecs1_ch(xyzzyaaab35,xyzzyaaaa35)
call dcopy(four_nitot,eivecs1_chscr(1,1,xyzzyaaaa35),1,eivecs_scr(1,1,&
&xyzzyaaab35,xyzzyaaaa35),1)
eivecs_valid(xyzzyaaaa35)=.true.
xyzzyaaae1(xyzzyaaaa35)=0
elseif(xyzzyaaae35/=0)then
xyzzyaaab35=buffer_move2_from_ii(xyzzyaaaa35)
xyzzyaaac35=buffer_move2_from_jj(xyzzyaaaa35)
if(xyzzyaaae1(xyzzyaaaa35)/=xyzzyaaae35)then
call get_eivecs(xyzzyaaae35)
call dcopy(four_nitot_netot,eivecs_scr(1,1,1,xyzzyaaae35),1,eivecs_scr&
&(1,1,1,xyzzyaaaa35),1)
endif
xyzzyaaag35=xyzzyaaao1(xyzzyaaaa35)
if(xyzzyaaag35/=0)then
xyzzyaaaf35=buffer_move_ion_from_ion(xyzzyaaag35)
xyzzyaaah35(1:3)=rion(1:3,xyzzyaaaf35)
rion(1:3,xyzzyaaaf35)=rion1_chscr(1:3,xyzzyaaag35)
call ei_distances(rele2ii_chscr(1,xyzzyaaaa35),nitot,rion,eivecs_scr(1&
&,1,xyzzyaaab35,xyzzyaaaa35))
call ei_distances(rele2jj_chscr(1,xyzzyaaaa35),nitot,rion,eivecs_scr(1&
&,1,xyzzyaaac35,xyzzyaaaa35))
rion(1:3,xyzzyaaaf35)=xyzzyaaah35(1:3)
else
call ei_distances(rele2ii_chscr(1,xyzzyaaaa35),nitot,rion,eivecs_scr(1&
&,1,xyzzyaaab35,xyzzyaaaa35))
call ei_distances(rele2jj_chscr(1,xyzzyaaaa35),nitot,rion,eivecs_scr(1&
&,1,xyzzyaaac35,xyzzyaaaa35))
endif
eivecs_valid(xyzzyaaaa35)=.true.
xyzzyaaae1(xyzzyaaaa35)=0
else
call get_rsele(xyzzyaaaa35)
xyzzyaaag35=xyzzyaaao1(xyzzyaaaa35)
if(xyzzyaaag35/=0)then
xyzzyaaaf35=buffer_move_ion_from_ion(xyzzyaaag35)
xyzzyaaah35(1:3)=rion(1:3,xyzzyaaaf35)
rion(1:3,xyzzyaaaf35)=rion1_chscr(1:3,xyzzyaaag35)
call ei_distances_all(netot,rele_scr(1,1,xyzzyaaaa35),nitot,rion,eivec&
&s_scr(1,1,1,xyzzyaaaa35))
rion(1:3,xyzzyaaaf35)=xyzzyaaah35(1:3)
else
call ei_distances_all(netot,rele_scr(1,1,xyzzyaaaa35),nitot,rion,eivec&
&s_scr(1,1,1,xyzzyaaaa35))
endif
eivecs_valid(xyzzyaaaa35)=.true.
xyzzyaaae1(xyzzyaaaa35)=0
endif
end subroutine get_eivecs
subroutine setup_storage_geom(nconfig)
use slaarnaaf, only : rele_config,sele_config
implicit none
integer,intent(in) :: nconfig
integer xyzzyaaaa36
allocate(xyzzyaaam1(3,netot,nconfig),xyzzyaaan1(netot,nconfig),stat=xy&
&zzyaaaa36)
call check_alloc(xyzzyaaaa36,'SETUP_STORAGE_GEOM','')
if(.not.allocated(rele_config))call errstop('SETUP_STORAGE_GEOM','No c&
&onfigs to load into storage?')
if(size(rele_config,3)<nconfig)call errstop('SETUP_STORAGE_GEOM','Not &
&enough configs to load into storage?')
call dcopy(three_netot*nconfig,rele_config(1,1,1),1,xyzzyaaam1(1,1,1),&
&1)
if(allocated(sele_config))then
xyzzyaaan1(1:netot,1:nconfig)=sele_config(1:netot,1:nconfig)
else
xyzzyaaan1=0
endif
end subroutine setup_storage_geom
subroutine finish_storage_geom
implicit none
deallocate(xyzzyaaam1,xyzzyaaan1)
end subroutine finish_storage_geom
subroutine load_from_storage_geom(is,icfg)
implicit none
integer,intent(in) :: is,icfg
call clear_scratch(is)
call dcopy(three_netot,xyzzyaaam1(1,1,icfg),1,rele_scr(1,1,is),1)
sele_scr(1:netot,is)=xyzzyaaan1(1:netot,icfg)
rsele_valid(is)=.true.
xyzzyaaac1(is)=0
end subroutine load_from_storage_geom
end module slaarnach

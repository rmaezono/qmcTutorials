module slaarnaaj
use dsp
use slaarnabg
use parallel
use slaarnach
use store
use slaarnacs
use format_utils,  only : wout,i2s,r2s,wordwrap
use slaarnabt,     only : dcopy
use slaarnaca,         only : have_veep,is_ae,zion,have_ae,have_ppots
use slaarnacc,only : initialize_random,random_seed_kw,ranx,ranx_gaussi&
&an,put_random_state,get_random_state
use run_control,   only : errstop,errstop2,errstop_master,errwarn,errw&
&arn_silent,timer,tcputime,exceeds_time_limit,check_alloc
implicit none
private
public dmc_main
integer,public :: limdmc,tpdmc,max_rec_attempts,redist_period,ebest_av&
&_window,corper_dmc,ndmcave,nconfig_prelim,redist_grp_size
real(dp),public :: cerefdmc,wdmcmin,wdmcmax,targ_wt,trip_popn,dmc_init&
&_eref,alimit,dmc_equil_fixpop
logical,public :: iaccum,ibran,lwdmc,lwdmc_fixpop,growth_estimator,wri&
&teout_dmc_hist,dmc_twist_av,use_init_eref,nc_dmc,poprenorm,nucleus_gf&
&_mods,dmc_reweight_configs,dmc_spacewarping
integer,allocatable :: xyzzyaaaa1(:),xyzzyaaab1(:)
real(dp),allocatable :: xyzzyaaac1(:),xyzzyaaad1(:),xyzzyaaae1(:),xyzz&
&yaaaf1(:)
logical,parameter :: xyzzyaaag1=.true.
real(dp) xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1,xyzzyaaak1,xyzzyaaal1
type configuration
type(config_geom),pointer :: pt_geom
type(config_wfn),pointer :: pt_wfn
type(config_dmc),pointer :: pt_dmc
type(configuration),pointer :: pt_next_config,pt_prev_config
end type configuration
type config_dmc
real(dp) eloc,sloc,wdmc
integer,pointer :: pt_ageelec(:)
real(dp),pointer :: pt_pure(:,:)
end type config_dmc
integer xyzzyaaam1,xyzzyaaan1,xyzzyaaao1,xyzzyaaap1
logical xyzzyaaaq1
type(configuration),pointer :: xyzzyaaar1
real(dp),allocatable :: xyzzyaaas1(:),xyzzyaaat1(:),xyzzyaaau1(:),xyzz&
&yaaav1(:,:)
integer,allocatable :: xyzzyaaaw1(:,:)
integer :: prelim_nconfig_written=0
integer xyzzyaaax1,xyzzyaaay1,xyzzyaaaz1
integer,parameter :: xyzzyaaba1=1,xyzzyaabb1=2,xyzzyaabc1=3,xyzzyaabd1&
&=4,xyzzyaabe1=3,xyzzyaabf1=2
integer,allocatable :: xyzzyaabg1(:,:)
logical xyzzyaabh1,xyzzyaabi1,xyzzyaabj1
integer xyzzyaabk1,xyzzyaabl1,xyzzyaabm1,xyzzyaabn1,xyzzyaabo1
integer,allocatable :: xyzzyaabp1(:),xyzzyaabq1(:),xyzzyaabr1(:),xyzzy&
&aabs1(:)
logical xyzzyaabt1,xyzzyaabu1,xyzzyaabv1
integer xyzzyaabw1
integer xyzzyaabx1
character(80) tmpr
contains
subroutine xyzzyaaby1(future_in,future_nitems_in,future_length_in)
implicit none
integer,intent(in) :: future_nitems_in,future_length_in
logical,intent(in) :: future_in
xyzzyaaam1=0
nullify(xyzzyaaar1)
xyzzyaaaq1=future_in
xyzzyaaan1=future_nitems_in
xyzzyaaao1=future_length_in
xyzzyaaap1=xyzzyaaan1*xyzzyaaao1
end subroutine xyzzyaaby1
subroutine xyzzyaabz1(pt_config)
implicit none
type(configuration),pointer :: pt_config
type(configuration),pointer :: xyzzyaaaa3
call timer('GEN_CONFIG',.true.)
if(associated(pt_config))call errstop('GEN_CONFIG','Config already all&
&ocated.')
allocate(pt_config,stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'GEN_CONFIG','container')
nullify(pt_config%pt_geom,pt_config%pt_wfn,pt_config%pt_dmc)
call gen_config(pt_config%pt_geom,pt_config%pt_wfn)
call xyzzyaaca1(pt_config%pt_dmc)
xyzzyaaam1=xyzzyaaam1+1
if(associated(xyzzyaaar1))then
xyzzyaaaa3=>xyzzyaaar1%pt_prev_config
xyzzyaaaa3%pt_next_config=>pt_config
pt_config%pt_prev_config=>xyzzyaaaa3
pt_config%pt_next_config=>xyzzyaaar1
xyzzyaaar1%pt_prev_config=>pt_config
else
pt_config%pt_prev_config=>pt_config
pt_config%pt_next_config=>pt_config
xyzzyaaar1=>pt_config
endif
call timer('GEN_CONFIG',.false.)
end subroutine xyzzyaabz1
subroutine xyzzyaaca1(pt_config)
implicit none
type(config_dmc),pointer :: pt_config
allocate(pt_config,stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'GEN_CONFIG_DMC','container')
pt_config%eloc=0.d0
pt_config%sloc=0.d0
pt_config%wdmc=1.d0
allocate(pt_config%pt_ageelec(netot),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'GEN_CONFIG_DMC','ageelec')
pt_config%pt_ageelec=0
if(xyzzyaaaq1)then
allocate(pt_config%pt_pure(xyzzyaaan1,xyzzyaaao1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'GEN_CONFIG_DMC','pure')
pt_config%pt_pure=0.d0
else
nullify(pt_config%pt_pure)
endif
end subroutine xyzzyaaca1
subroutine xyzzyaacb1(is,pt_config,rele,sele,eloc,sloc,wdmc,ageelec)
implicit none
integer,intent(in) :: is
integer,intent(out) :: sele(netot),ageelec(netot)
real(dp),intent(out) :: rele(3,netot),eloc,sloc,wdmc
type(configuration),pointer :: pt_config
call timer('LOAD_FROM_PT',.true.)
call load_from_pt(is,pt_config%pt_geom,pt_config%pt_wfn,rele,sele)
call xyzzyaacc1(pt_config%pt_dmc,eloc,sloc,wdmc,ageelec)
call timer('LOAD_FROM_PT',.false.)
end subroutine xyzzyaacb1
subroutine xyzzyaacc1(pt_config,eloc,sloc,wdmc,ageelec)
implicit none
integer,intent(out) :: ageelec(netot)
real(dp),intent(out) :: eloc,sloc,wdmc
type(config_dmc),pointer :: pt_config
eloc=pt_config%eloc
sloc=pt_config%sloc
wdmc=pt_config%wdmc
ageelec(1:netot)=pt_config%pt_ageelec(1:netot)
end subroutine xyzzyaacc1
subroutine xyzzyaacd1(pt_config,pureitems,istep)
implicit none
integer,intent(in) :: istep
real(dp),intent(inout) :: pureitems(xyzzyaaan1)
real(dp),pointer :: xyzzyaaaa7(:)
type(configuration),pointer :: pt_config
xyzzyaaaa7=>pt_config%pt_dmc%pt_pure(1:xyzzyaaan1,istep)
call dcopy(xyzzyaaan1,xyzzyaaaa7,1,pureitems,1)
end subroutine xyzzyaacd1
subroutine xyzzyaace1(is,pt_config,eloc,sloc,wdmc,ageelec)
implicit none
integer,intent(in) :: is,ageelec(netot)
real(dp),intent(in) :: eloc,sloc,wdmc
type(configuration),pointer :: pt_config
call timer('SAVE_TO_PT',.true.)
call save_to_pt(is,pt_config%pt_geom,pt_config%pt_wfn)
call xyzzyaacf1(pt_config%pt_dmc,eloc,sloc,wdmc,ageelec)
call timer('SAVE_TO_PT',.false.)
end subroutine xyzzyaace1
subroutine xyzzyaacf1(pt_config,eloc,sloc,wdmc,ageelec)
implicit none
integer,intent(in) :: ageelec(netot)
real(dp),intent(in) :: eloc,sloc,wdmc
type(config_dmc),pointer :: pt_config
pt_config%eloc=eloc
pt_config%sloc=sloc
pt_config%wdmc=wdmc
pt_config%pt_ageelec(1:netot)=ageelec(1:netot)
end subroutine xyzzyaacf1
subroutine xyzzyaacg1(pt_config,pureitems,istep)
implicit none
integer,intent(in) :: istep
real(dp),intent(in) :: pureitems(xyzzyaaan1)
real(dp),pointer :: xyzzyaaaa10(:)
type(configuration),pointer :: pt_config
xyzzyaaaa10=>pt_config%pt_dmc%pt_pure(1:xyzzyaaan1,istep)
call dcopy(xyzzyaaan1,pureitems,1,xyzzyaaaa10,1)
end subroutine xyzzyaacg1
subroutine xyzzyaach1(pt_config)
implicit none
type(configuration),pointer :: pt_config,pt_prev,pt_next
logical xyzzyaaaa11
call timer('DELETE_CONFIG',.true.)
if(.not.associated(pt_config))call errstop('DELETE_CONFIG','Config not&
& allocated.')
pt_prev=>pt_config%pt_prev_config
pt_next=>pt_config%pt_next_config
xyzzyaaaa11=associated(pt_next,pt_config)
if(associated(pt_config,xyzzyaaar1))xyzzyaaar1=>pt_next
call xyzzyaaci1(pt_config%pt_dmc)
call delete_config(pt_config%pt_geom,pt_config%pt_wfn)
deallocate(pt_config)
xyzzyaaam1=xyzzyaaam1-1
if(xyzzyaaaa11)then
nullify(xyzzyaaar1,pt_config)
else
pt_prev%pt_next_config=>pt_next
pt_next%pt_prev_config=>pt_prev
pt_config=>pt_prev
endif
call timer('DELETE_CONFIG',.false.)
end subroutine xyzzyaach1
subroutine xyzzyaaci1(pt_config)
implicit none
type(config_dmc),pointer :: pt_config
if(xyzzyaaaq1)deallocate(pt_config%pt_pure)
deallocate(pt_config%pt_ageelec)
deallocate(pt_config)
end subroutine xyzzyaaci1
subroutine xyzzyaacj1(pt_from,pt_to)
implicit none
type(configuration),pointer :: pt_from,pt_to
call timer('COPY_CONFIG',.true.)
if(.not.(associated(pt_from).and.associated(pt_to)))call errstop('COPY&
&_CONFIG','Trying to copy between nonexistent configs.')
call copy_config(pt_from%pt_geom,pt_to%pt_geom,pt_from%pt_wfn,pt_to%pt&
&_wfn)
call xyzzyaack1(pt_from%pt_dmc,pt_to%pt_dmc)
call timer('COPY_CONFIG',.false.)
end subroutine xyzzyaacj1
subroutine xyzzyaack1(pt_from,pt_to)
implicit none
integer xyzzyaaaa14
type(config_dmc),pointer :: pt_from,pt_to
pt_to%eloc=pt_from%eloc
pt_to%sloc=pt_from%sloc
pt_to%wdmc=pt_from%wdmc
do xyzzyaaaa14=1,netot
pt_to%pt_ageelec(xyzzyaaaa14)=pt_from%pt_ageelec(xyzzyaaaa14)
enddo
if(xyzzyaaaq1)call dcopy(xyzzyaaap1,pt_from%pt_pure,1,pt_to%pt_pure,1)
end subroutine xyzzyaack1
subroutine xyzzyaacl1(pt_config,icfg)
implicit none
integer,intent(in) :: icfg
integer xyzzyaaaa15
type(configuration),pointer :: pt_config
if(.not.associated(xyzzyaaar1))call errstop('FIND_CONFIG','List is emp&
&ty.')
if(icfg<=0)call errstop('FIND_CONFIG','Bad CONFIG_NO index.')
pt_config=>xyzzyaaar1
do xyzzyaaaa15=2,icfg
pt_config=>pt_config%pt_next_config
if(associated(pt_config,xyzzyaaar1))call errstop('FIND_CONFIG','CONFIG&
&_NO seems to exceed the number of configs in the list.')
enddo
end subroutine xyzzyaacl1
subroutine xyzzyaacm1
implicit none
type(configuration),pointer :: xyzzyaaaa16
if(associated(xyzzyaaar1))xyzzyaaaa16=>xyzzyaaar1%pt_prev_config
do while(associated(xyzzyaaar1))
call xyzzyaach1(xyzzyaaaa16)
enddo
end subroutine xyzzyaacm1
subroutine xyzzyaacn1(pt_config,k)
implicit none
integer,intent(in) :: k
type(configuration),pointer :: pt_config
call timer('CONFIG_TO_PT',.true.)
if(.not.associated(pt_config))call errstop('CONFIG_TO_PT','Been passed&
& a nullified config pointer.')
call config_to_pt(pt_config%pt_geom,pt_config%pt_wfn,k)
call xyzzyaaco1(pt_config%pt_dmc,k)
call timer('CONFIG_TO_PT',.false.)
end subroutine xyzzyaacn1
subroutine xyzzyaaco1(pt_config,k)
use slaarnaaf,only : etot_config,stot_config,wdmc_config
implicit none
integer,intent(in) :: k
type(config_dmc),pointer :: pt_config
pt_config%eloc=etot_config(k)
pt_config%sloc=stot_config(k)
pt_config%wdmc=wdmc_config(k)
pt_config%pt_ageelec(1:netot)=0
end subroutine xyzzyaaco1
subroutine xyzzyaacp1(pt_config)
use slaarnaaf,only : add_config
implicit none
type(configuration),pointer :: pt_config
call timer('PT_TO_CONFIG',.true.)
if(.not.associated(pt_config))call errstop('PT_TO_CONFIG','Been passed&
& a nullified config pointer.')
call add_config()
call pt_to_config(pt_config%pt_geom,pt_config%pt_wfn)
call xyzzyaacq1(pt_config%pt_dmc)
call timer('PT_TO_CONFIG',.false.)
end subroutine xyzzyaacp1
subroutine xyzzyaacq1(pt_config)
use slaarnaaf,only : add_config
implicit none
type(config_dmc),pointer :: pt_config
call add_config(modify=.true.,etot=pt_config%eloc,stot=pt_config%sloc,&
&wdmc=pt_config%wdmc)
end subroutine xyzzyaacq1
subroutine xyzzyaacr1(kmax)
implicit none
integer,intent(in) :: kmax
call timer('REDIST_ALLOCATIONS',.true.)
call redist_allocations(kmax)
call xyzzyaacs1(kmax)
call timer('REDIST_ALLOCATIONS',.false.)
end subroutine xyzzyaacr1
subroutine xyzzyaacs1(kmax)
implicit none
integer,intent(in) :: kmax
allocate(xyzzyaaas1(kmax),xyzzyaaat1(kmax),xyzzyaaau1(kmax),xyzzyaaaw1&
&(netot,kmax),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'REDIST_ALLOCATIONS_DMC','')
if(xyzzyaaaq1)then
allocate(xyzzyaaav1(xyzzyaaap1,kmax),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'REDIST_ALLOCATIONS_DMC','future')
xyzzyaaav1=0.d0
endif
end subroutine xyzzyaacs1
subroutine xyzzyaact1(pt_config,k,blocking)
implicit none
integer,intent(in) :: k
logical,intent(in) :: blocking
type(configuration),pointer :: pt_config
call timer('REDIST_LOAD',.true.)
if(.not.associated(pt_config))call errstop('REDIST_LOAD','Been passed &
&a nullified config pointer.')
call redist_load(pt_config%pt_geom,pt_config%pt_wfn,k,blocking)
call xyzzyaacu1(pt_config%pt_dmc,k)
call timer('REDIST_LOAD',.false.)
end subroutine xyzzyaact1
subroutine xyzzyaacu1(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_dmc),pointer :: pt_config
xyzzyaaas1(k)=pt_config%eloc
xyzzyaaat1(k)=pt_config%sloc
xyzzyaaau1(k)=pt_config%wdmc
xyzzyaaaw1(1:netot,k)=pt_config%pt_ageelec(1:netot)
if(xyzzyaaaq1)call dcopy(xyzzyaaap1,pt_config%pt_pure,1,xyzzyaaav1(1,k&
&),1)
end subroutine xyzzyaacu1
subroutine xyzzyaacv1(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
call timer('REDIST_SEND',.true.)
call redist_send(jnode,kbase,k,nbt,reqbase,blocking)
call xyzzyaacw1(jnode,kbase,k,nbt,reqbase,blocking)
call timer('REDIST_SEND',.false.)
end subroutine xyzzyaacv1
subroutine xyzzyaacw1(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa26
if(blocking)then
call qmpi_ssend(xyzzyaaas1(kbase+1:kbase+k),jnode,move_msg,'REDIST_SEN&
&D','eloc')
call qmpi_ssend(xyzzyaaat1(kbase+1:kbase+k),jnode,move_msg,'REDIST_SEN&
&D','sloc')
call qmpi_ssend(xyzzyaaau1(kbase+1:kbase+k),jnode,move_msg,'REDIST_SEN&
&D','wdmc')
call qmpi_ssend(xyzzyaaaw1(1:netot,kbase+1:kbase+k),jnode,move_msg,'RE&
&DIST_SEND','ageelec')
if(xyzzyaaaq1)call qmpi_ssend(xyzzyaaav1(1:xyzzyaaap1,kbase+1:kbase+k)&
&,jnode,move_msg,'REDIST_SEND','pt_pure')
else
nbt=nbt+1
call mpi_isend(xyzzyaaas1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa26,ierror)
call checkmpi(ierror,'sending eloc ('//trim(i2s(nbt))//') in redist_se&
&nd_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa26
nbt=nbt+1
call mpi_isend(xyzzyaaat1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa26,ierror)
call checkmpi(ierror,'sending sloc ('//trim(i2s(nbt))//') in redist_se&
&nd_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa26
nbt=nbt+1
call mpi_isend(xyzzyaaau1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa26,ierror)
call checkmpi(ierror,'sending wdmc ('//trim(i2s(nbt))//') in redist_se&
&nd_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa26
nbt=nbt+1
call mpi_isend(xyzzyaaaw1(1,kbase+1),netot*k,mpi_integer,jnode,nbt,mpi&
&_comm_world,xyzzyaaaa26,ierror)
call checkmpi(ierror,'sending ageelec_trf ('//trim(i2s(nbt))//') in re&
&dist_send_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa26
if(xyzzyaaaq1)then
nbt=nbt+1
call mpi_isend(xyzzyaaav1(1,kbase+1),xyzzyaaap1*k,mpi_double_precision&
&,jnode,nbt,mpi_comm_world,xyzzyaaaa26,ierror)
call checkmpi(ierror,'sending pure_trf ('//trim(i2s(nbt))//') in redis&
&t_send_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa26
endif
endif
end subroutine xyzzyaacw1
subroutine xyzzyaacx1(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
call timer('REDIST_RECV',.true.)
call redist_recv(jnode,kbase,k,nbt,reqbase,blocking)
call xyzzyaacy1(jnode,kbase,k,nbt,reqbase,blocking)
call timer('REDIST_RECV',.false.)
end subroutine xyzzyaacx1
subroutine xyzzyaacy1(jnode,kbase,k,nbt,reqbase,blocking)
implicit none
integer,intent(in) :: jnode,kbase,k,reqbase
integer,intent(inout) :: nbt
logical,intent(in) :: blocking
integer xyzzyaaaa28
if(blocking)then
call qmpi_recv(xyzzyaaas1(kbase+1:kbase+k),jnode,move_msg,'REDIST_RECV&
&','eloc')
call qmpi_recv(xyzzyaaat1(kbase+1:kbase+k),jnode,move_msg,'REDIST_RECV&
&','sloc')
call qmpi_recv(xyzzyaaau1(kbase+1:kbase+k),jnode,move_msg,'REDIST_RECV&
&','wdmc')
call qmpi_recv(xyzzyaaaw1(1:netot,kbase+1:kbase+k),jnode,move_msg,'RED&
&IST_RECV','ageelec')
if(xyzzyaaaq1)call qmpi_recv(xyzzyaaav1(1:xyzzyaaap1,kbase+1:kbase+k),&
&jnode,move_msg,'REDIST_RECV','pt_pure')
else
nbt=nbt+1
call mpi_irecv(xyzzyaaas1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa28,ierror)
call checkmpi(ierror,'receiving eloc ('//trim(i2s(nbt))//') in redist_&
&recv_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa28
nbt=nbt+1
call mpi_irecv(xyzzyaaat1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa28,ierror)
call checkmpi(ierror,'receiving sloc ('//trim(i2s(nbt))//') in redist_&
&recv_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa28
nbt=nbt+1
call mpi_irecv(xyzzyaaau1(kbase+1),k,mpi_double_precision,jnode,nbt,mp&
&i_comm_world,xyzzyaaaa28,ierror)
call checkmpi(ierror,'receiving wdmc ('//trim(i2s(nbt))//') in redist_&
&recv_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa28
nbt=nbt+1
call mpi_irecv(xyzzyaaaw1(1,kbase+1),netot*k,mpi_integer,jnode,nbt,mpi&
&_comm_world,xyzzyaaaa28,ierror)
call checkmpi(ierror,'receiving ageelec_trf ('//trim(i2s(nbt))//') in &
&redist_recv_dmc')
nbreq(reqbase+nbt)=xyzzyaaaa28
if(xyzzyaaaq1)then
nbt=nbt+1
call mpi_irecv(xyzzyaaav1(1,kbase+1),xyzzyaaap1*k,mpi_double_precision&
&,jnode,nbt,mpi_comm_world,xyzzyaaaa28,ierror)
nbreq(reqbase+nbt)=xyzzyaaaa28
endif
endif
end subroutine xyzzyaacy1
subroutine xyzzyaacz1(pt_config,k)
implicit none
integer,intent(in) :: k
type(configuration),pointer :: pt_config
call timer('REDIST_SAVE',.true.)
if(.not.associated(pt_config))call errstop('REDIST_SAVE','Been passed &
&a nullified config pointer.')
call redist_save(pt_config%pt_geom,pt_config%pt_wfn,k)
call xyzzyaada1(pt_config%pt_dmc,k)
call timer('REDIST_SAVE',.false.)
end subroutine xyzzyaacz1
subroutine xyzzyaada1(pt_config,k)
implicit none
integer,intent(in) :: k
type(config_dmc),pointer :: pt_config
pt_config%eloc=xyzzyaaas1(k)
pt_config%sloc=xyzzyaaat1(k)
pt_config%wdmc=xyzzyaaau1(k)
pt_config%pt_ageelec(1:netot)=xyzzyaaaw1(1:netot,k)
if(xyzzyaaaq1)call dcopy(xyzzyaaap1,xyzzyaaav1(1,k),1,pt_config%pt_pur&
&e,1)
end subroutine xyzzyaada1
subroutine xyzzyaadb1
implicit none
call timer('REDIST_DEALLOCATIONS',.true.)
call redist_deallocations
call xyzzyaadc1
call timer('REDIST_DEALLOCATIONS',.false.)
end subroutine xyzzyaadb1
subroutine xyzzyaadc1
implicit none
deallocate(xyzzyaaas1,xyzzyaaat1,xyzzyaaau1,xyzzyaaaw1)
if(xyzzyaaaq1)deallocate(xyzzyaaav1)
end subroutine xyzzyaadc1
subroutine xyzzyaadd1(redist_type,nconfig,imult_node,ntransfers)
use slaarnaaf
implicit none
integer,intent(out),optional :: ntransfers
integer,intent(inout) :: nconfig
integer,intent(inout),optional :: imult_node(nconfig)
character(4),intent(in) :: redist_type
select case (trim(adjustl(redist_type)))
case('SEND')
call xyzzyaadg1(nconfig,imult_node,ntransfers)
case('RECV')
call xyzzyaadh1(nconfig)
case('SERC')
call xyzzyaadi1(nconfig,imult_node,ntransfers)
case default
call errstop('BRANCH_AND_REDIST','Unknown redist algorithm')
end select
end subroutine xyzzyaadd1
subroutine xyzzyaade1
implicit none
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34,xyzzyaaad34,xyzzyaaae34,xy&
&zzyaaaf34,xyzzyaaag34
integer,allocatable :: xyzzyaaah34(:),xyzzyaaai34(:)
if(allocated(xyzzyaabp1))return
if(redist_grp_size<2)call errstop('SETUP_REDIST_GROUPS','REDIST_GRP_SI&
&ZE must be a positive integer greater than 1.')
if(nnodes<redist_grp_size)redist_grp_size=nnodes
xyzzyaabk1=nnodes/redist_grp_size
allocate(xyzzyaabp1(xyzzyaabk1),xyzzyaabq1(xyzzyaabk1),xyzzyaabr1(xyzz&
&yaabk1),xyzzyaaai34(xyzzyaabk1+1),rg_comm(xyzzyaabk1+1,2),stat=xyzzya&
&aag34)
call check_alloc(xyzzyaaag34,'SETUP_REDIST_GROUPS','rg_nnodes etc.')
xyzzyaaaa34=(nnodes-xyzzyaabk1*redist_grp_size)/xyzzyaabk1
xyzzyaabp1(:)=redist_grp_size+xyzzyaaaa34
xyzzyaaab34=mod(nnodes-xyzzyaabk1*redist_grp_size,xyzzyaabk1)
do xyzzyaaac34=1,xyzzyaaab34
xyzzyaabp1(xyzzyaaac34)=xyzzyaabp1(xyzzyaaac34)+1
enddo
xyzzyaabq1(1)=0
xyzzyaabr1(1)=xyzzyaabp1(1)-1
xyzzyaabm1=1
do xyzzyaaaa34=2,xyzzyaabk1
xyzzyaabq1(xyzzyaaaa34)=xyzzyaabr1(xyzzyaaaa34-1)+1
xyzzyaabr1(xyzzyaaaa34)=xyzzyaabq1(xyzzyaaaa34)+xyzzyaabp1(xyzzyaaaa34&
&)-1
if(my_node>=xyzzyaabq1(xyzzyaaaa34))xyzzyaabm1=xyzzyaaaa34
enddo
if(any(my_node==xyzzyaabq1(:)))then
xyzzyaabt1=.true.
xyzzyaabl1=my_node
else
xyzzyaabt1=.false.
xyzzyaabl1=xyzzyaabq1(xyzzyaabm1)
endif
xyzzyaabv1=.true.
xyzzyaabn1=redist_grp_size/2
if(my_node-xyzzyaabq1(xyzzyaabm1)<xyzzyaabn1)then
xyzzyaabu1=.true.
else
xyzzyaabu1=.false.
endif
xyzzyaabo1=1
call mpi_comm_group(mpi_comm_world,xyzzyaaaf34,ierror)
call checkmpi(ierror,'mpi_comm_group in setup_redist_groups')
do xyzzyaaaa34=1,xyzzyaabk1
allocate(xyzzyaaah34(xyzzyaabp1(xyzzyaaaa34)),stat=xyzzyaaag34)
call check_alloc(xyzzyaaag34,'SETUP_REDIST_GROUPS','rg_ranks')
xyzzyaaac34=0
do xyzzyaaab34=xyzzyaabq1(xyzzyaaaa34),xyzzyaabr1(xyzzyaaaa34)
xyzzyaaac34=xyzzyaaac34+1
xyzzyaaah34(xyzzyaaac34)=xyzzyaaab34
enddo
call mpi_group_incl(xyzzyaaaf34,xyzzyaabp1(xyzzyaaaa34),xyzzyaaah34,xy&
&zzyaaai34(xyzzyaaaa34),ierror)
call checkmpi(ierror,'mpi_group_incl in setup_redist_groups')
call mpi_comm_create(mpi_comm_world,xyzzyaaai34(xyzzyaaaa34),rg_comm(x&
&yzzyaaaa34,1),ierror)
call checkmpi(ierror,'mpi_comm_create in setup_redist_groups')
deallocate(xyzzyaaah34)
enddo
do xyzzyaaaa34=1,xyzzyaabk1
allocate(xyzzyaaah34(xyzzyaabp1(xyzzyaaaa34)),stat=xyzzyaaag34)
call check_alloc(xyzzyaaag34,'SETUP_REDIST_GROUPS','rg_ranks <2>')
xyzzyaaac34=0
xyzzyaaad34=xyzzyaabq1(xyzzyaaaa34)+xyzzyaabn1
xyzzyaaae34=xyzzyaabr1(xyzzyaaaa34)+xyzzyaabn1
do xyzzyaaab34=xyzzyaaad34,xyzzyaaae34
xyzzyaaac34=xyzzyaaac34+1
xyzzyaaah34(xyzzyaaac34)=mod(xyzzyaaab34,nnodes)
enddo
call mpi_group_incl(xyzzyaaaf34,xyzzyaabp1(xyzzyaaaa34),xyzzyaaah34,xy&
&zzyaaai34(xyzzyaaaa34),ierror)
call checkmpi(ierror,'mpi_group_incl in setup_redist_groups <2>')
call mpi_comm_create(mpi_comm_world,xyzzyaaai34(xyzzyaaaa34),rg_comm(x&
&yzzyaaaa34,2),ierror)
call checkmpi(ierror,'mpi_comm_create in setup_redist_groups <2>')
deallocate(xyzzyaaah34)
enddo
allocate(xyzzyaaah34(xyzzyaabk1),stat=xyzzyaaag34)
call check_alloc(xyzzyaaag34,'SETUP_REDIST_GROUPS','rg_ranks <3>')
do xyzzyaaaa34=1,xyzzyaabk1
xyzzyaaah34(xyzzyaaaa34)=xyzzyaabq1(xyzzyaaaa34)
enddo
call mpi_group_incl(xyzzyaaaf34,xyzzyaabk1,xyzzyaaah34,xyzzyaaai34(xyz&
&zyaabk1+1),ierror)
call checkmpi(ierror,'mpi_group_incl in setup_redist_groups <3>')
call mpi_comm_create(mpi_comm_world,xyzzyaaai34(xyzzyaabk1+1),rg_comm(&
&xyzzyaabk1+1,1),ierror)
call checkmpi(ierror,'mpi_comm_create in setup_redist_groups <3>')
rg_comm(xyzzyaabk1+1,2)=rg_comm(xyzzyaabk1+1,1)
deallocate(xyzzyaaah34,xyzzyaaai34)
if(nc_dmc.or.poprenorm.or.dmc_equil_fixpop>0.d0)then
allocate(xyzzyaabs1(xyzzyaabk1),stat=xyzzyaaag34)
call check_alloc(xyzzyaaag34,'SETUP_REDIST_GROUPS','target_pop')
xyzzyaabs1(:)=nint(targ_wt)/xyzzyaabk1
do xyzzyaaaa34=1,mod(nint(targ_wt),xyzzyaabk1)
xyzzyaabs1(xyzzyaaaa34)=xyzzyaabs1(xyzzyaaaa34)+1
enddo
endif
end subroutine xyzzyaade1
subroutine xyzzyaadf1
implicit none
if(xyzzyaabv1)then
if(xyzzyaabu1)then
xyzzyaabm1=xyzzyaabm1-1
if(xyzzyaabm1==0)xyzzyaabm1=xyzzyaabk1
else
xyzzyaabl1=mod(xyzzyaabl1+xyzzyaabp1(xyzzyaabm1),nnodes)
endif
xyzzyaabq1(xyzzyaabm1)=xyzzyaabq1(xyzzyaabm1)+xyzzyaabn1
xyzzyaabv1=.false.
xyzzyaabo1=2
else
xyzzyaabq1(xyzzyaabm1)=xyzzyaabq1(xyzzyaabm1)-xyzzyaabn1
if(xyzzyaabu1)then
xyzzyaabm1=mod(xyzzyaabm1,xyzzyaabk1)+1
else
xyzzyaabl1=xyzzyaabl1-xyzzyaabp1(xyzzyaabm1)
if(xyzzyaabl1<0)xyzzyaabl1=xyzzyaabl1+nnodes
endif
xyzzyaabv1=.true.
xyzzyaabo1=1
endif
end subroutine xyzzyaadf1
subroutine xyzzyaadg1(nconfig,imult_node,ntransfers)
use slaarnaaf
use slaarnabt, only : partial_rank_int
implicit none
integer,intent(out) :: ntransfers
integer,intent(inout) :: nconfig,imult_node(nconfig)
integer xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36,xyzzyaaad36,xyzzyaaae36,xy&
&zzyaaaf36,xyzzyaaag36,xyzzyaaah36,xyzzyaaai36,xyzzyaaaj36,xyzzyaaak36&
&,xyzzyaaal36,xyzzyaaam36,xyzzyaaan36
integer,allocatable :: xyzzyaaao36(:),xyzzyaaap36(:),xyzzyaaaq36(:),xy&
&zzyaaar36(:),xyzzyaaas36(:),xyzzyaaat36(:),xyzzyaaau36(:,:),xyzzyaaav&
&36(:,:,:)
type(configuration),pointer :: xyzzyaaaw36,xyzzyaaax36
call timer('BRANCH_AND_REDIST_SEND',.true.)
xyzzyaaan36=nconfig
ntransfers=0
if(nnodes>1)then
xyzzyaaae36=xyzzyaabp1(xyzzyaabm1)-1
allocate(xyzzyaaao36(0:xyzzyaaae36),xyzzyaaap36(0:xyzzyaaae36),xyzzyaa&
&aq36(0:xyzzyaaae36),xyzzyaaar36(0:xyzzyaaae36),xyzzyaaas36(0:xyzzyaaa&
&e36),xyzzyaaat36(0:xyzzyaaae36),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SEND','iact,ncon,tpop e&
&tc..')
xyzzyaaao36=0
xyzzyaaax1=0
xyzzyaabw1=0
xyzzyaabi1=.false.
xyzzyaabh1=.false.
call timer('REDIST_GATHER_NCONFIG',.true.)
xyzzyaaad36=rg_comm(xyzzyaabm1,xyzzyaabo1)
call qmpi_allgather(nconfig,xyzzyaaap36(0:xyzzyaaae36),xyzzyaaad36,'BR&
&ANCH_AND_REDIST_SEND','ncon')
if(any(xyzzyaaap36(:)>0))then
continue
else
call errstop_master('BRANCH_AND_REDIST_SEND','Population died out in a&
& redist group. Should not happen - please report this.')
endif
call timer('REDIST_GATHER_NCONFIG',.false.)
if(xyzzyaabt1)then
if(am_master)call timer('REDIST_GATHER_IMULT',.true.)
xyzzyaaac36=maxval(xyzzyaaap36(:))
allocate(xyzzyaaau36(xyzzyaaac36,0:xyzzyaabp1(xyzzyaabm1)-1),stat=xyzz&
&yaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SEND','imult')
xyzzyaaau36=0
do xyzzyaaaa36=0,xyzzyaabp1(xyzzyaabm1)-1
xyzzyaaag36=xyzzyaabq1(xyzzyaabm1)+xyzzyaaaa36
if(xyzzyaaag36>nnodes-1)xyzzyaaag36=xyzzyaaag36-nnodes
if(xyzzyaaag36==xyzzyaabl1)then
if(xyzzyaaap36(xyzzyaaaa36)>0)then
xyzzyaaau36(1:xyzzyaaap36(xyzzyaaaa36),xyzzyaaaa36)=imult_node(1:xyzzy&
&aaap36(xyzzyaaaa36))
xyzzyaaas36(xyzzyaaaa36)=sum(xyzzyaaau36(1:xyzzyaaap36(xyzzyaaaa36),xy&
&zzyaaaa36))
else
xyzzyaaas36(xyzzyaaaa36)=0
endif
cycle
endif
if(xyzzyaaap36(xyzzyaaaa36)>0)then
call qmpi_recv(xyzzyaaau36(1:xyzzyaaap36(xyzzyaaaa36),xyzzyaaaa36),xyz&
&zyaaag36,1,'BRANCH_AND_REDIST_SEND','imult')
endif
xyzzyaaas36(xyzzyaaaa36)=sum(xyzzyaaau36(1:xyzzyaaap36(xyzzyaaaa36),xy&
&zzyaaaa36))
enddo
if(am_master)call timer('REDIST_GATHER_IMULT',.false.)
if(am_master)call timer('REDIST_CALC_INSTR',.true.)
xyzzyaaaz1=sum(xyzzyaaas36(:))
xyzzyaaaq36(:)=xyzzyaaaz1/xyzzyaabp1(xyzzyaabm1)
xyzzyaaaa36=mod(xyzzyaaaz1,xyzzyaabp1(xyzzyaabm1))
if(xyzzyaaaa36>0)then
if(xyzzyaaaa36==1)then
xyzzyaaac36=maxloc(xyzzyaaas36(:),1)-1
xyzzyaaaq36(xyzzyaaac36)=xyzzyaaaq36(xyzzyaaac36)+1
else
call partial_rank_int(xyzzyaaas36(:),xyzzyaaat36(:),xyzzyaaaa36)
do xyzzyaaab36=0,xyzzyaaaa36-1
xyzzyaaac36=xyzzyaaat36(xyzzyaaab36)-1
xyzzyaaaq36(xyzzyaaac36)=xyzzyaaaq36(xyzzyaaac36)+1
enddo
endif
endif
xyzzyaaar36(:)=xyzzyaaas36(:)-xyzzyaaaq36(:)
allocate(xyzzyaaav36(3,max_no_messages,0:xyzzyaabp1(xyzzyaabm1)-1),sta&
&t=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SEND','instructions')
xyzzyaaav36=0
loop_send: do while(any(xyzzyaaar36(:)/=0))
xyzzyaaag36=maxloc(xyzzyaaar36(:),1)-1
if(xyzzyaaar36(xyzzyaaag36)<=0)exit
if(xyzzyaaao36(xyzzyaaag36)>=max_no_messages)then
xyzzyaaaq36(xyzzyaaag36)=xyzzyaaas36(xyzzyaaag36)
xyzzyaaar36(xyzzyaaag36)=0
cycle loop_send
endif
loop_recv: do
if(any(xyzzyaaar36(:)==-xyzzyaaar36(xyzzyaaag36)))then
xyzzyaaai36=minloc(xyzzyaaar36(:),1,xyzzyaaar36(:)==-xyzzyaaar36(xyzzy&
&aaag36))-1
elseif(any(xyzzyaaar36(:)<-xyzzyaaar36(xyzzyaaag36)))then
xyzzyaaai36=minloc(xyzzyaaar36(:),1,xyzzyaaar36(:)<-xyzzyaaar36(xyzzya&
&aag36))-1
else
xyzzyaaai36=minloc(xyzzyaaar36(:),1)-1
endif
if(xyzzyaaar36(xyzzyaaai36)>=0)exit loop_send
if(-xyzzyaaao36(xyzzyaaai36)>=max_no_messages)then
xyzzyaaaq36(xyzzyaaai36)=xyzzyaaas36(xyzzyaaai36)
xyzzyaaar36(xyzzyaaai36)=0
cycle loop_recv
endif
do
xyzzyaaaf36=0
if(any(xyzzyaaau36(1:xyzzyaaap36(xyzzyaaag36),xyzzyaaag36)==xyzzyaaar3&
&6(xyzzyaaag36)+1))then
xyzzyaaaf36=maxloc(xyzzyaaau36(1:xyzzyaaap36(xyzzyaaag36),xyzzyaaag36)&
&,1,xyzzyaaau36(1:xyzzyaaap36(xyzzyaaag36),xyzzyaaag36)==xyzzyaaar36(x&
&yzzyaaag36)+1)
elseif(any(xyzzyaaau36(1:xyzzyaaap36(xyzzyaaag36),xyzzyaaag36)>xyzzyaa&
&ar36(xyzzyaaag36)+1))then
xyzzyaaaf36=maxloc(xyzzyaaau36(1:xyzzyaaap36(xyzzyaaag36),xyzzyaaag36)&
&,1,xyzzyaaau36(1:xyzzyaaap36(xyzzyaaag36),xyzzyaaag36)>xyzzyaaar36(xy&
&zzyaaag36)+1)
else
xyzzyaaaf36=maxloc(xyzzyaaau36(1:xyzzyaaap36(xyzzyaaag36),xyzzyaaag36)&
&,1)
endif
if(xyzzyaaaf36==0)call errstop('BRANCH_AND_REDIST_SEND','Bug in redist&
&ribution algorithm <1>.')
xyzzyaaaj36=min(xyzzyaaau36(xyzzyaaaf36,xyzzyaaag36),xyzzyaaar36(xyzzy&
&aaag36),-xyzzyaaar36(xyzzyaaai36))
if(xyzzyaaaj36==0)call errstop('BRANCH_AND_REDIST_SEND','Bug in redist&
&ribution algorithm <2>.')
ntransfers=ntransfers+1
xyzzyaaao36(xyzzyaaag36)=xyzzyaaao36(xyzzyaaag36)+1
xyzzyaaav36(xyzzyaabc1,xyzzyaaao36(xyzzyaaag36),xyzzyaaag36)=xyzzyaaaf&
&36
xyzzyaaav36(xyzzyaabb1,xyzzyaaao36(xyzzyaaag36),xyzzyaaag36)=xyzzyaaaj&
&36
xyzzyaaav36(xyzzyaaba1,xyzzyaaao36(xyzzyaaag36),xyzzyaaag36)=mod(xyzzy&
&aabq1(xyzzyaabm1)+xyzzyaaai36,nnodes)
xyzzyaaas36(xyzzyaaag36)=xyzzyaaas36(xyzzyaaag36)-xyzzyaaaj36
xyzzyaaau36(xyzzyaaaf36,xyzzyaaag36)=xyzzyaaau36(xyzzyaaaf36,xyzzyaaag&
&36)-xyzzyaaaj36
xyzzyaaar36(xyzzyaaag36)=xyzzyaaas36(xyzzyaaag36)-xyzzyaaaq36(xyzzyaaa&
&g36)
xyzzyaaao36(xyzzyaaai36)=xyzzyaaao36(xyzzyaaai36)-1
xyzzyaaav36(xyzzyaabb1,-xyzzyaaao36(xyzzyaaai36),xyzzyaaai36)=xyzzyaaa&
&j36
xyzzyaaav36(xyzzyaaba1,-xyzzyaaao36(xyzzyaaai36),xyzzyaaai36)=mod(xyzz&
&yaabq1(xyzzyaabm1)+xyzzyaaag36,nnodes)
xyzzyaaas36(xyzzyaaai36)=xyzzyaaas36(xyzzyaaai36)+xyzzyaaaj36
xyzzyaaar36(xyzzyaaai36)=xyzzyaaas36(xyzzyaaai36)-xyzzyaaaq36(xyzzyaaa&
&i36)
if(xyzzyaaao36(xyzzyaaag36)>=max_no_messages)then
xyzzyaaaq36(xyzzyaaag36)=xyzzyaaas36(xyzzyaaag36)
xyzzyaaar36(xyzzyaaag36)=0
endif
if(-xyzzyaaao36(xyzzyaaai36)>=max_no_messages)then
xyzzyaaaq36(xyzzyaaai36)=xyzzyaaas36(xyzzyaaai36)
xyzzyaaar36(xyzzyaaai36)=0
endif
if(xyzzyaaar36(xyzzyaaag36)==0)cycle loop_send
if(xyzzyaaar36(xyzzyaaai36)==0)cycle loop_recv
enddo
enddo loop_recv
enddo loop_send
if(am_master)call timer('REDIST_CALC_INSTR',.false.)
if(am_master)call timer('REDIST_SEND_INSTR',.true.)
do xyzzyaaaa36=0,xyzzyaabp1(xyzzyaabm1)-1
xyzzyaaag36=xyzzyaabq1(xyzzyaabm1)+xyzzyaaaa36
if(xyzzyaaag36>nnodes-1)xyzzyaaag36=xyzzyaaag36-nnodes
if(xyzzyaaag36==xyzzyaabl1)then
xyzzyaaax1=xyzzyaaao36(xyzzyaaaa36)
if(xyzzyaaax1>0)then
xyzzyaabi1=.true.
allocate(xyzzyaabg1(4,xyzzyaaax1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SEND','instructions_nod&
&e')
xyzzyaabg1=0
xyzzyaabg1(1:xyzzyaabe1,1:xyzzyaaax1)=xyzzyaaav36(1:xyzzyaabe1,1:xyzzy&
&aaax1,xyzzyaaaa36)
elseif(xyzzyaaax1<0)then
xyzzyaabh1=.true.
xyzzyaaax1=abs(xyzzyaaax1)
allocate(xyzzyaabg1(4,xyzzyaaax1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SEND','instructions_nod&
&e')
xyzzyaabg1=0
xyzzyaabg1(1:xyzzyaabf1,1:xyzzyaaax1)=xyzzyaaav36(1:xyzzyaabf1,1:xyzzy&
&aaax1,xyzzyaaaa36)
endif
cycle
endif
call qmpi_ssend(xyzzyaaao36(xyzzyaaaa36),xyzzyaaag36,2,'BRANCH_AND_RED&
&IST_SEND','instructions')
if(xyzzyaaao36(xyzzyaaaa36)>0)then
call qmpi_ssend(xyzzyaaav36(1:xyzzyaabe1,1:xyzzyaaao36(xyzzyaaaa36),xy&
&zzyaaaa36),xyzzyaaag36,3,'BRANCH_AND_REDIST_SEND','instructions')
elseif(xyzzyaaao36(xyzzyaaaa36)<0)then
call qmpi_ssend(xyzzyaaav36(1:xyzzyaabf1,1:-xyzzyaaao36(xyzzyaaaa36),x&
&yzzyaaaa36),xyzzyaaag36,4,'BRANCH_AND_REDIST_SEND','instructions')
endif
enddo
deallocate(xyzzyaaav36,xyzzyaaau36)
if(am_master)call timer('REDIST_SEND_INSTR',.false.)
else
if(nconfig>0)then
call qmpi_ssend(imult_node(1:nconfig),xyzzyaabl1,1,'BRANCH_AND_REDIST_&
&SEND','imult_node')
endif
call qmpi_recv(xyzzyaaax1,xyzzyaabl1,2,'BRANCH_AND_REDIST_SEND','instr&
&uction (1)')
if(xyzzyaaax1>0)then
xyzzyaabi1=.true.
allocate(xyzzyaabg1(4,xyzzyaaax1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SEND','instructions_nod&
&e')
xyzzyaabg1=0
call qmpi_recv(xyzzyaabg1(1:xyzzyaabe1,1:xyzzyaaax1),xyzzyaabl1,3,'BRA&
&NCH_AND_REDIST_SEND','instruction (2a)')
elseif(xyzzyaaax1<0)then
xyzzyaabh1=.true.
xyzzyaaax1=abs(xyzzyaaax1)
allocate(xyzzyaabg1(4,xyzzyaaax1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SEND','instructions_nod&
&e')
xyzzyaabg1=0
call qmpi_recv(xyzzyaabg1(1:xyzzyaabf1,1:xyzzyaaax1),xyzzyaabl1,4,'BRA&
&NCH_AND_REDIST_SEND','instruction (2b)')
endif
endif
call timer('REST_OF_REDIST',.true.)
if(xyzzyaabi1.or.xyzzyaabh1)then
xyzzyaaal36=0
do xyzzyaaaa36=1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaabd1,xyzzyaaaa36)/=0)cycle
xyzzyaabj1=.true.
xyzzyaaac36=1
xyzzyaaah36=xyzzyaaaa36
do xyzzyaaab36=xyzzyaaaa36+1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaaba1,xyzzyaaab36)==xyzzyaabg1(xyzzyaaba1,xyzzyaaaa&
&36))then
if(xyzzyaabj1)then
xyzzyaabg1(xyzzyaabd1,xyzzyaaah36)=xyzzyaaab36
xyzzyaabj1=.false.
else
xyzzyaabg1(xyzzyaabd1,xyzzyaaah36)=-xyzzyaaab36
endif
xyzzyaaac36=xyzzyaaac36+1
xyzzyaaah36=xyzzyaaab36
endif
enddo
if(xyzzyaaah36==xyzzyaaaa36)then
xyzzyaabg1(xyzzyaabd1,xyzzyaaaa36)=xyzzyaaaa36
else
xyzzyaabg1(xyzzyaabd1,xyzzyaaah36)=-xyzzyaaah36
endif
xyzzyaaal36=xyzzyaaal36+xyzzyaaac36
enddo
if(xyzzyaaal36<1)call errstop2('BRANCH_AND_REDIST_SEND','Can''t count &
&on node ',my_node)
call xyzzyaacr1(xyzzyaaal36)
endif
if(xyzzyaabh1)then
xyzzyaaak36=0
xyzzyaaay1=0
do xyzzyaaaa36=1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaabd1,xyzzyaaaa36)<0)cycle
xyzzyaaam36=0
xyzzyaaai36=xyzzyaabg1(xyzzyaaba1,xyzzyaaaa36)
xyzzyaaab36=xyzzyaaaa36
xyzzyaaac36=1
do
if(abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab36))==xyzzyaaab36)exit
xyzzyaaab36=abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab36))
xyzzyaaac36=xyzzyaaac36+1
enddo
call xyzzyaacx1(xyzzyaaai36,xyzzyaaak36,xyzzyaaac36,xyzzyaaam36,xyzzya&
&aay1,.false.)
xyzzyaaak36=xyzzyaaak36+xyzzyaaac36
xyzzyaaay1=xyzzyaaay1+xyzzyaaam36
enddo
endif
call qmc_barrier
if(xyzzyaabi1)then
xyzzyaaak36=0
xyzzyaaay1=0
do xyzzyaaaa36=1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaabd1,xyzzyaaaa36)<0)cycle
xyzzyaaam36=0
xyzzyaaai36=xyzzyaabg1(xyzzyaaba1,xyzzyaaaa36)
xyzzyaaab36=xyzzyaaaa36
xyzzyaaac36=1
do
xyzzyaaaf36=xyzzyaabg1(xyzzyaabc1,xyzzyaaab36)
xyzzyaaaj36=xyzzyaabg1(xyzzyaabb1,xyzzyaaab36)
nullify(xyzzyaaaw36)
call xyzzyaacl1(xyzzyaaaw36,xyzzyaaaf36)
call xyzzyaact1(xyzzyaaaw36,xyzzyaaak36+xyzzyaaac36,.false.)
imult_node(xyzzyaaaf36)=imult_node(xyzzyaaaf36)-xyzzyaaaj36
if(abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab36))==xyzzyaaab36)exit
xyzzyaaab36=abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab36))
xyzzyaaac36=xyzzyaaac36+1
enddo
call xyzzyaacv1(xyzzyaaai36,xyzzyaaak36,xyzzyaaac36,xyzzyaaam36,xyzzya&
&aay1,.false.)
xyzzyaaak36=xyzzyaaak36+xyzzyaaac36
xyzzyaaay1=xyzzyaaay1+xyzzyaaam36
enddo
endif
deallocate(xyzzyaaao36,xyzzyaaap36,xyzzyaaaq36,xyzzyaaar36,xyzzyaaas36&
&,xyzzyaaat36)
endif
call timer('REST_OF_REDIST',.false.)
if(.not.associated(xyzzyaaar1))return
call timer('BRANCH',.true.)
xyzzyaaaw36=>xyzzyaaar1
do xyzzyaaaf36=1,xyzzyaaan36
xyzzyaaaj36=imult_node(xyzzyaaaf36)-1
select case(xyzzyaaaj36)
case(-1)
call xyzzyaach1(xyzzyaaaw36)
case(0)
continue
case(1:)
do xyzzyaaab36=1,xyzzyaaaj36
nullify(xyzzyaaax36)
call xyzzyaabz1(xyzzyaaax36)
call xyzzyaacj1(xyzzyaaaw36,xyzzyaaax36)
enddo
case default
call errstop('BRANCH_AND_REDIST_SEND','Confused on node '//trim(i2s(my&
&_node)))
end select
nconfig=nconfig+xyzzyaaaj36
if(.not.associated(xyzzyaaaw36))exit
xyzzyaaaw36=>xyzzyaaaw36%pt_next_config
enddo
call timer('BRANCH',.false.)
if(xyzzyaaam1/=nconfig)then
if(nnodes>1)then
call wout()
call wout('In node '//trim(i2s(my_node))//':')
call wout()
endif
call wout('Number of configs (nconfig): '//trim(i2s(nconfig)))
call wout('Number of configs (tally)  : '//trim(i2s(xyzzyaaam1)))
call errstop('BRANCH_AND_REDIST_SEND','Bug: lost track of the config l&
&ist.')
endif
call timer('BRANCH_AND_REDIST_SEND',.false.)
end subroutine xyzzyaadg1
subroutine xyzzyaadh1(nconfig)
implicit none
integer,intent(inout) :: nconfig
integer xyzzyaaaa37,xyzzyaaab37,xyzzyaaac37,xyzzyaaad37,xyzzyaaae37,xy&
&zzyaaaf37
type(configuration),pointer :: xyzzyaaag37,xyzzyaaah37
call timer('BRANCH_AND_REDIST_RECV',.true.)
call timer('REDIST_MPI_WAITALL',.true.)
call mpi_waitall(xyzzyaaay1,nbreq(1:xyzzyaaay1),mpi_statuses_ignore,ie&
&rror)
call checkmpi(ierror,'MPI_WAITALL in branch_and_redist_recv')
call timer('REDIST_MPI_WAITALL',.false.)
if(xyzzyaabi1)then
call xyzzyaadb1
elseif(xyzzyaabh1)then
xyzzyaaae37=0
do xyzzyaaaa37=1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaabd1,xyzzyaaaa37)<0)cycle
xyzzyaaab37=xyzzyaaaa37
xyzzyaaac37=1
do
nullify(xyzzyaaag37)
call xyzzyaabz1(xyzzyaaag37)
call xyzzyaacz1(xyzzyaaag37,xyzzyaaae37+xyzzyaaac37)
xyzzyaaaf37=xyzzyaabg1(xyzzyaabb1,xyzzyaaab37)
if(xyzzyaaaf37>1)then
do xyzzyaaad37=1,xyzzyaaaf37-1
nullify(xyzzyaaah37)
call xyzzyaabz1(xyzzyaaah37)
call xyzzyaacj1(xyzzyaaag37,xyzzyaaah37)
enddo
endif
nullify(xyzzyaaag37)
nconfig=nconfig+xyzzyaaaf37
if(abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab37))==xyzzyaaab37)exit
xyzzyaaab37=abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab37))
xyzzyaaac37=xyzzyaaac37+1
enddo
xyzzyaaae37=xyzzyaaae37+xyzzyaaac37
enddo
call xyzzyaadb1
endif
deallocate(xyzzyaabg1)
call timer('BRANCH_AND_REDIST_RECV',.false.)
end subroutine xyzzyaadh1
subroutine xyzzyaadi1(nconfig,imult_node,ntransfers)
use slaarnaaf
use slaarnabt, only : partial_rank_int
implicit none
integer,intent(out) :: ntransfers
integer,intent(inout) :: nconfig,imult_node(nconfig)
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38,xyzzyaaad38,xyzzyaaae38,xy&
&zzyaaaf38,xyzzyaaag38,xyzzyaaah38,xyzzyaaai38,xyzzyaaaj38,xyzzyaaak38&
&,xyzzyaaal38,xyzzyaaam38,xyzzyaaan38,xyzzyaaao38
integer,allocatable :: xyzzyaaap38(:),xyzzyaaaq38(:),xyzzyaaar38(:),xy&
&zzyaaas38(:),xyzzyaaat38(:),xyzzyaaau38(:),xyzzyaaav38(:,:),xyzzyaaaw&
&38(:,:,:)
type(configuration),pointer :: xyzzyaaax38,xyzzyaaay38
call timer('BRANCH_AND_REDIST_SERC',.true.)
xyzzyaaao38=nconfig
ntransfers=0
if(nnodes>1)then
if(xyzzyaabk1>1)call xyzzyaadf1
xyzzyaaaf38=xyzzyaabp1(xyzzyaabm1)-1
allocate(xyzzyaaap38(0:xyzzyaaaf38),xyzzyaaaq38(0:xyzzyaaaf38),xyzzyaa&
&ar38(0:xyzzyaaaf38),xyzzyaaas38(0:xyzzyaaaf38),xyzzyaaat38(0:xyzzyaaa&
&f38),xyzzyaaau38(0:xyzzyaaaf38),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SENDRECV','iact,ncon,tp&
&op etc..')
xyzzyaaap38=0
xyzzyaaax1=0
xyzzyaabw1=0
xyzzyaabi1=.false.
xyzzyaabh1=.false.
call timer('REDIST_GATHER_NCONFIG',.true.)
xyzzyaaae38=rg_comm(xyzzyaabm1,xyzzyaabo1)
call qmpi_allgather(nconfig,xyzzyaaaq38(0:xyzzyaaaf38),xyzzyaaae38,'BR&
&ANCH_AND_REDIST_SENDRECV','ncon')
if(any(xyzzyaaaq38(:)>0))then
continue
else
call errstop_master('BRANCH_AND_REDIST_SENDRECV','Population died out &
&in a redist group. Should not happen - please report this.')
endif
call timer('REDIST_GATHER_NCONFIG',.false.)
if(xyzzyaabt1)then
if(am_master)call timer('REDIST_GATHER_IMULT',.true.)
xyzzyaaac38=maxval(xyzzyaaaq38(:))
allocate(xyzzyaaav38(xyzzyaaac38,0:xyzzyaabp1(xyzzyaabm1)-1),stat=xyzz&
&yaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SENDRECV','imult')
xyzzyaaav38=0
do xyzzyaaaa38=0,xyzzyaabp1(xyzzyaabm1)-1
xyzzyaaah38=xyzzyaabq1(xyzzyaabm1)+xyzzyaaaa38
if(xyzzyaaah38>nnodes-1)xyzzyaaah38=xyzzyaaah38-nnodes
if(xyzzyaaah38==xyzzyaabl1)then
if(xyzzyaaaq38(xyzzyaaaa38)>0)then
xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaaa38),xyzzyaaaa38)=imult_node(1:xyzzy&
&aaaq38(xyzzyaaaa38))
xyzzyaaat38(xyzzyaaaa38)=sum(xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaaa38),xy&
&zzyaaaa38))
else
xyzzyaaat38(xyzzyaaaa38)=0
endif
cycle
endif
if(xyzzyaaaq38(xyzzyaaaa38)>0)then
call qmpi_recv(xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaaa38),xyzzyaaaa38),xyz&
&zyaaah38,1,'BRANCH_AND_REDIST_SENDRECV','imult')
endif
xyzzyaaat38(xyzzyaaaa38)=sum(xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaaa38),xy&
&zzyaaaa38))
enddo
if(am_master)call timer('REDIST_GATHER_IMULT',.false.)
if(am_master)call timer('REDIST_CALC_INSTR',.true.)
xyzzyaaaz1=sum(xyzzyaaat38(:))
xyzzyaaar38(:)=xyzzyaaaz1/xyzzyaabp1(xyzzyaabm1)
xyzzyaaaa38=mod(xyzzyaaaz1,xyzzyaabp1(xyzzyaabm1))
if(xyzzyaaaa38>0)then
if(xyzzyaaaa38==1)then
xyzzyaaac38=maxloc(xyzzyaaat38(:),1)-1
xyzzyaaar38(xyzzyaaac38)=xyzzyaaar38(xyzzyaaac38)+1
else
call partial_rank_int(xyzzyaaat38(:),xyzzyaaau38(:),xyzzyaaaa38)
do xyzzyaaab38=0,xyzzyaaaa38-1
xyzzyaaac38=xyzzyaaau38(xyzzyaaab38)-1
xyzzyaaar38(xyzzyaaac38)=xyzzyaaar38(xyzzyaaac38)+1
enddo
endif
endif
xyzzyaaas38(:)=xyzzyaaat38(:)-xyzzyaaar38(:)
allocate(xyzzyaaaw38(3,max_no_messages,0:xyzzyaabp1(xyzzyaabm1)-1),sta&
&t=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SENDRECV','instructions&
&')
xyzzyaaaw38=0
loop_send: do while(any(xyzzyaaas38(:)/=0))
xyzzyaaah38=maxloc(xyzzyaaas38(:),1)-1
if(xyzzyaaas38(xyzzyaaah38)<=0)exit
if(xyzzyaaap38(xyzzyaaah38)>=max_no_messages)then
xyzzyaaar38(xyzzyaaah38)=xyzzyaaat38(xyzzyaaah38)
xyzzyaaas38(xyzzyaaah38)=0
cycle loop_send
endif
loop_recv: do
if(any(xyzzyaaas38(:)==-xyzzyaaas38(xyzzyaaah38)))then
xyzzyaaaj38=minloc(xyzzyaaas38(:),1,xyzzyaaas38(:)==-xyzzyaaas38(xyzzy&
&aaah38))-1
elseif(any(xyzzyaaas38(:)<-xyzzyaaas38(xyzzyaaah38)))then
xyzzyaaaj38=minloc(xyzzyaaas38(:),1,xyzzyaaas38(:)<-xyzzyaaas38(xyzzya&
&aah38))-1
else
xyzzyaaaj38=minloc(xyzzyaaas38(:),1)-1
endif
if(xyzzyaaas38(xyzzyaaaj38)>=0)exit loop_send
if(-xyzzyaaap38(xyzzyaaaj38)>=max_no_messages)then
xyzzyaaar38(xyzzyaaaj38)=xyzzyaaat38(xyzzyaaaj38)
xyzzyaaas38(xyzzyaaaj38)=0
cycle loop_recv
endif
do
xyzzyaaag38=0
if(any(xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaah38),xyzzyaaah38)==xyzzyaaas3&
&8(xyzzyaaah38)+1))then
xyzzyaaag38=maxloc(xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaah38),xyzzyaaah38)&
&,1,xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaah38),xyzzyaaah38)==xyzzyaaas38(x&
&yzzyaaah38)+1)
elseif(any(xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaah38),xyzzyaaah38)>xyzzyaa&
&as38(xyzzyaaah38)+1))then
xyzzyaaag38=maxloc(xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaah38),xyzzyaaah38)&
&,1,xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaah38),xyzzyaaah38)>xyzzyaaas38(xy&
&zzyaaah38)+1)
else
xyzzyaaag38=maxloc(xyzzyaaav38(1:xyzzyaaaq38(xyzzyaaah38),xyzzyaaah38)&
&,1)
endif
if(xyzzyaaag38==0)call errstop('BRANCH_AND_REDIST_SENDRECV','Bug in re&
&distribution algorithm <1>.')
xyzzyaaak38=min(xyzzyaaav38(xyzzyaaag38,xyzzyaaah38),xyzzyaaas38(xyzzy&
&aaah38),-xyzzyaaas38(xyzzyaaaj38))
if(xyzzyaaak38==0)call errstop('BRANCH_AND_REDIST_SENDRECV','Bug in re&
&distribution algorithm <2>.')
ntransfers=ntransfers+1
xyzzyaaap38(xyzzyaaah38)=xyzzyaaap38(xyzzyaaah38)+1
xyzzyaaaw38(xyzzyaabc1,xyzzyaaap38(xyzzyaaah38),xyzzyaaah38)=xyzzyaaag&
&38
xyzzyaaaw38(xyzzyaabb1,xyzzyaaap38(xyzzyaaah38),xyzzyaaah38)=xyzzyaaak&
&38
xyzzyaaaw38(xyzzyaaba1,xyzzyaaap38(xyzzyaaah38),xyzzyaaah38)=mod(xyzzy&
&aabq1(xyzzyaabm1)+xyzzyaaaj38,nnodes)
xyzzyaaat38(xyzzyaaah38)=xyzzyaaat38(xyzzyaaah38)-xyzzyaaak38
xyzzyaaav38(xyzzyaaag38,xyzzyaaah38)=xyzzyaaav38(xyzzyaaag38,xyzzyaaah&
&38)-xyzzyaaak38
xyzzyaaas38(xyzzyaaah38)=xyzzyaaat38(xyzzyaaah38)-xyzzyaaar38(xyzzyaaa&
&h38)
xyzzyaaap38(xyzzyaaaj38)=xyzzyaaap38(xyzzyaaaj38)-1
xyzzyaaaw38(xyzzyaabb1,-xyzzyaaap38(xyzzyaaaj38),xyzzyaaaj38)=xyzzyaaa&
&k38
xyzzyaaaw38(xyzzyaaba1,-xyzzyaaap38(xyzzyaaaj38),xyzzyaaaj38)=mod(xyzz&
&yaabq1(xyzzyaabm1)+xyzzyaaah38,nnodes)
xyzzyaaat38(xyzzyaaaj38)=xyzzyaaat38(xyzzyaaaj38)+xyzzyaaak38
xyzzyaaas38(xyzzyaaaj38)=xyzzyaaat38(xyzzyaaaj38)-xyzzyaaar38(xyzzyaaa&
&j38)
if(xyzzyaaap38(xyzzyaaah38)>=max_no_messages)then
xyzzyaaar38(xyzzyaaah38)=xyzzyaaat38(xyzzyaaah38)
xyzzyaaas38(xyzzyaaah38)=0
endif
if(-xyzzyaaap38(xyzzyaaaj38)>=max_no_messages)then
xyzzyaaar38(xyzzyaaaj38)=xyzzyaaat38(xyzzyaaaj38)
xyzzyaaas38(xyzzyaaaj38)=0
endif
if(xyzzyaaas38(xyzzyaaah38)==0)cycle loop_send
if(xyzzyaaas38(xyzzyaaaj38)==0)cycle loop_recv
enddo
enddo loop_recv
enddo loop_send
if(am_master)call timer('REDIST_CALC_INSTR',.false.)
if(am_master)call timer('REDIST_SEND_INSTR',.true.)
do xyzzyaaaa38=0,xyzzyaabp1(xyzzyaabm1)-1
xyzzyaaah38=xyzzyaabq1(xyzzyaabm1)+xyzzyaaaa38
if(xyzzyaaah38>nnodes-1)xyzzyaaah38=xyzzyaaah38-nnodes
if(xyzzyaaah38==xyzzyaabl1)then
xyzzyaaax1=xyzzyaaap38(xyzzyaaaa38)
if(xyzzyaaax1>0)then
xyzzyaabi1=.true.
allocate(xyzzyaabg1(4,xyzzyaaax1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SENDRECV','instructions&
&_node')
xyzzyaabg1=0
xyzzyaabg1(1:xyzzyaabe1,1:xyzzyaaax1)=xyzzyaaaw38(1:xyzzyaabe1,1:xyzzy&
&aaax1,xyzzyaaaa38)
elseif(xyzzyaaax1<0)then
xyzzyaabh1=.true.
xyzzyaaax1=abs(xyzzyaaax1)
allocate(xyzzyaabg1(4,xyzzyaaax1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SENDRECV','instructions&
&_node')
xyzzyaabg1=0
xyzzyaabg1(1:xyzzyaabf1,1:xyzzyaaax1)=xyzzyaaaw38(1:xyzzyaabf1,1:xyzzy&
&aaax1,xyzzyaaaa38)
endif
cycle
endif
call qmpi_ssend(xyzzyaaap38(xyzzyaaaa38),xyzzyaaah38,2,'BRANCH_AND_RED&
&IST_SENDRECV','instructions')
if(xyzzyaaap38(xyzzyaaaa38)>0)then
call qmpi_ssend(xyzzyaaaw38(1:xyzzyaabe1,1:xyzzyaaap38(xyzzyaaaa38),xy&
&zzyaaaa38),xyzzyaaah38,3,'BRANCH_AND_REDIST_SENDRECV','instructions')
elseif(xyzzyaaap38(xyzzyaaaa38)<0)then
call qmpi_ssend(xyzzyaaaw38(1:xyzzyaabf1,1:-xyzzyaaap38(xyzzyaaaa38),x&
&yzzyaaaa38),xyzzyaaah38,4,'BRANCH_AND_REDIST_SENDRECV','instructions'&
&)
endif
enddo
deallocate(xyzzyaaaw38,xyzzyaaav38)
if(am_master)call timer('REDIST_SEND_INSTR',.false.)
else
if(nconfig>0)then
call qmpi_ssend(imult_node(1:nconfig),xyzzyaabl1,1,'BRANCH_AND_REDIST_&
&SENDRECV','imult_node')
endif
call qmpi_recv(xyzzyaaax1,xyzzyaabl1,2,'BRANCH_AND_REDIST_SENDRECV','i&
&nstruction (1)')
if(xyzzyaaax1>0)then
xyzzyaabi1=.true.
allocate(xyzzyaabg1(4,xyzzyaaax1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SENDRECV','instructions&
&_node')
xyzzyaabg1=0
call qmpi_recv(xyzzyaabg1(1:xyzzyaabe1,1:xyzzyaaax1),xyzzyaabl1,3,'BRA&
&NCH_AND_REDIST_SENDRECV','instruction (2a)')
elseif(xyzzyaaax1<0)then
xyzzyaabh1=.true.
xyzzyaaax1=abs(xyzzyaaax1)
allocate(xyzzyaabg1(4,xyzzyaaax1),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'BRANCH_AND_REDIST_SENDRECV','instructions&
&_node')
xyzzyaabg1=0
call qmpi_recv(xyzzyaabg1(1:xyzzyaabf1,1:xyzzyaaax1),xyzzyaabl1,4,'BRA&
&NCH_AND_REDIST_SENDRECV','instruction (2b)')
endif
endif
call timer('REST_OF_REDIST_SENDRECV',.true.)
if(xyzzyaabi1.or.xyzzyaabh1)then
xyzzyaaam38=0
do xyzzyaaaa38=1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaabd1,xyzzyaaaa38)/=0)cycle
xyzzyaabj1=.true.
xyzzyaaac38=1
xyzzyaaai38=xyzzyaaaa38
do xyzzyaaab38=xyzzyaaaa38+1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaaba1,xyzzyaaab38)==xyzzyaabg1(xyzzyaaba1,xyzzyaaaa&
&38))then
if(xyzzyaabj1)then
xyzzyaabg1(xyzzyaabd1,xyzzyaaai38)=xyzzyaaab38
xyzzyaabj1=.false.
else
xyzzyaabg1(xyzzyaabd1,xyzzyaaai38)=-xyzzyaaab38
endif
xyzzyaaac38=xyzzyaaac38+1
xyzzyaaai38=xyzzyaaab38
endif
enddo
if(xyzzyaaai38==xyzzyaaaa38)then
xyzzyaabg1(xyzzyaabd1,xyzzyaaaa38)=xyzzyaaaa38
else
xyzzyaabg1(xyzzyaabd1,xyzzyaaai38)=-xyzzyaaai38
endif
xyzzyaaam38=xyzzyaaam38+xyzzyaaac38
enddo
if(xyzzyaaam38<1)call errstop2('BRANCH_AND_REDIST_SENDRECV','Can''t co&
&unt on node ',my_node)
call xyzzyaacr1(xyzzyaaam38)
endif
if(xyzzyaabi1)then
xyzzyaaal38=0
xyzzyaaan38=0
do xyzzyaaaa38=1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaabd1,xyzzyaaaa38)<0)cycle
xyzzyaaaj38=xyzzyaabg1(xyzzyaaba1,xyzzyaaaa38)
xyzzyaaab38=xyzzyaaaa38
xyzzyaaac38=1
do
xyzzyaaag38=xyzzyaabg1(xyzzyaabc1,xyzzyaaab38)
xyzzyaaak38=xyzzyaabg1(xyzzyaabb1,xyzzyaaab38)
nullify(xyzzyaaax38)
call xyzzyaacl1(xyzzyaaax38,xyzzyaaag38)
call xyzzyaact1(xyzzyaaax38,xyzzyaaal38+xyzzyaaac38,.true.)
imult_node(xyzzyaaag38)=imult_node(xyzzyaaag38)-xyzzyaaak38
if(abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab38))==xyzzyaaab38)exit
xyzzyaaab38=abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab38))
xyzzyaaac38=xyzzyaaac38+1
enddo
call xyzzyaacv1(xyzzyaaaj38,xyzzyaaal38,xyzzyaaac38,xyzzyaaan38,0,.tru&
&e.)
xyzzyaaal38=xyzzyaaal38+xyzzyaaac38
enddo
call xyzzyaadb1
deallocate(xyzzyaabg1)
elseif(xyzzyaabh1)then
xyzzyaaal38=0
xyzzyaaan38=0
do xyzzyaaaa38=1,xyzzyaaax1
if(xyzzyaabg1(xyzzyaabd1,xyzzyaaaa38)<0)cycle
xyzzyaaaj38=xyzzyaabg1(xyzzyaaba1,xyzzyaaaa38)
xyzzyaaab38=xyzzyaaaa38
xyzzyaaac38=1
do
if(abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab38))==xyzzyaaab38)exit
xyzzyaaab38=abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab38))
xyzzyaaac38=xyzzyaaac38+1
enddo
call xyzzyaacx1(xyzzyaaaj38,xyzzyaaal38,xyzzyaaac38,xyzzyaaan38,0,.tru&
&e.)
xyzzyaaab38=xyzzyaaaa38
xyzzyaaac38=1
do
nullify(xyzzyaaax38)
call xyzzyaabz1(xyzzyaaax38)
call xyzzyaacz1(xyzzyaaax38,xyzzyaaal38+xyzzyaaac38)
xyzzyaaak38=xyzzyaabg1(xyzzyaabb1,xyzzyaaab38)
if(xyzzyaaak38>1)then
do xyzzyaaad38=1,xyzzyaaak38-1
nullify(xyzzyaaay38)
call xyzzyaabz1(xyzzyaaay38)
call xyzzyaacj1(xyzzyaaax38,xyzzyaaay38)
enddo
endif
nullify(xyzzyaaax38)
nconfig=nconfig+xyzzyaaak38
if(abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab38))==xyzzyaaab38)exit
xyzzyaaab38=abs(xyzzyaabg1(xyzzyaabd1,xyzzyaaab38))
xyzzyaaac38=xyzzyaaac38+1
enddo
xyzzyaaal38=xyzzyaaal38+xyzzyaaac38
enddo
call xyzzyaadb1
deallocate(xyzzyaabg1)
endif
deallocate(xyzzyaaap38,xyzzyaaaq38,xyzzyaaar38,xyzzyaaas38,xyzzyaaat38&
&,xyzzyaaau38)
endif
call timer('REST_OF_REDIST_SENDRECV',.false.)
xyzzyaabh1=.false.
xyzzyaabi1=.false.
if(.not.associated(xyzzyaaar1))return
call timer('BRANCH',.true.)
xyzzyaaax38=>xyzzyaaar1
do xyzzyaaag38=1,xyzzyaaao38
xyzzyaaak38=imult_node(xyzzyaaag38)-1
select case(xyzzyaaak38)
case(-1)
call xyzzyaach1(xyzzyaaax38)
case(0)
continue
case(1:)
do xyzzyaaab38=1,xyzzyaaak38
nullify(xyzzyaaay38)
call xyzzyaabz1(xyzzyaaay38)
call xyzzyaacj1(xyzzyaaax38,xyzzyaaay38)
enddo
case default
call errstop('BRANCH_AND_REDIST_SENDRECV','Confused on node '//trim(i2&
&s(my_node)))
end select
nconfig=nconfig+xyzzyaaak38
if(.not.associated(xyzzyaaax38))exit
xyzzyaaax38=>xyzzyaaax38%pt_next_config
enddo
call timer('BRANCH',.false.)
if(xyzzyaaam1/=nconfig)then
if(nnodes>1)then
call wout()
call wout('In node '//trim(i2s(my_node))//':')
call wout()
endif
call wout('Number of configs (nconfig): '//trim(i2s(nconfig)))
call wout('Number of configs (tally)  : '//trim(i2s(xyzzyaaam1)))
call errstop('BRANCH_AND_REDIST_SENDRECV','Bug: lost track of the conf&
&ig list.')
endif
call timer('BRANCH_AND_REDIST_SERC',.false.)
end subroutine xyzzyaadi1
real(dp) function xyzzyaadj1(rold,rnew,vold,vnew,dt,diff_const)
implicit none
real(dp),intent(in) :: rold(3),rnew(3),vold(3),vnew(3),dt,diff_const
integer xyzzyaaaa39
xyzzyaaaa39=dimensionality
xyzzyaadj1=sum((vold(1:xyzzyaaaa39)+vnew(1:xyzzyaaaa39))*(rold(1:xyzzy&
&aaaa39)-rnew(1:xyzzyaaaa39)+diff_const*dt*(vold(1:xyzzyaaaa39)-vnew(1&
&:xyzzyaaaa39))))
end function xyzzyaadj1
real(dp) function xyzzyaadk1(rold,rnew,vold,vnew,dt,diff_const)
use slaarnabt,only : exp_limit
implicit none
real(dp),intent(in) :: rold(3),rnew(3),vold(3),vnew(3),dt,diff_const
xyzzyaadk1=exp_limit(xyzzyaadj1(rold,rnew,vold,vnew,dt,diff_const))
end function xyzzyaadk1
subroutine xyzzyaadl1
implicit none
if(.not.allocated(xyzzyaaae1))then
allocate(xyzzyaaac1(no_difftypes),xyzzyaaad1(no_difftypes),xyzzyaaaa1(&
&no_difftypes),xyzzyaaab1(no_difftypes),xyzzyaaae1(no_difftypes),xyzzy&
&aaaf1(no_difftypes),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'TEFF_MOVE_START','')
endif
xyzzyaaac1=0.d0
xyzzyaaad1=0.d0
if(.not.xyzzyaaag1)then
xyzzyaaaa1=0
xyzzyaaab1=0
endif
end subroutine xyzzyaadl1
subroutine xyzzyaadm1(diffsq,prob,accepted,ispin)
implicit none
integer,intent(in) :: ispin
real(dp),intent(in) :: diffsq,prob
logical,intent(in) :: accepted
integer xyzzyaaaa42
xyzzyaaaa42=which_difftype(ispin)
xyzzyaaac1(xyzzyaaaa42)=xyzzyaaac1(xyzzyaaaa42)+diffsq
if(xyzzyaaag1)then
xyzzyaaad1(xyzzyaaaa42)=xyzzyaaad1(xyzzyaaaa42)+diffsq*prob
else
if(accepted)then
xyzzyaaad1(xyzzyaaaa42)=xyzzyaaad1(xyzzyaaaa42)+diffsq
xyzzyaaab1(xyzzyaaaa42)=xyzzyaaab1(xyzzyaaaa42)+1
endif
xyzzyaaaa1(xyzzyaaaa42)=xyzzyaaaa1(xyzzyaaaa42)+1
endif
end subroutine xyzzyaadm1
subroutine xyzzyaadn1(dt,dteff)
implicit none
real(dp),intent(in) :: dt
real(dp),intent(out) :: dteff
real(dp) xyzzyaaaa43
if(xyzzyaaag1)then
xyzzyaaaa43=0.d0
if(any(xyzzyaaac1>0.d0))xyzzyaaaa43=sum(difftype_mass*xyzzyaaad1)/sum(&
&difftype_mass*xyzzyaaac1)
else
xyzzyaaaf1=0.d0
where(xyzzyaaaa1>0)xyzzyaaaf1=xyzzyaaac1/xyzzyaaaa1
xyzzyaaae1=0.d0
where(xyzzyaaab1>0)xyzzyaaae1=xyzzyaaad1/xyzzyaaab1
xyzzyaaaa43=0.d0
if(any(xyzzyaaaf1>0.d0))xyzzyaaaa43=sum(difftype_mass*xyzzyaaaa1*xyzzy&
&aaae1)/sum(difftype_mass*xyzzyaaaa1*xyzzyaaaf1)
endif
dteff=dt*xyzzyaaaa43
end subroutine xyzzyaadn1
subroutine xyzzyaado1(a,tau,dvel_in,dvel_out,diff_const)
implicit none
real(dp),intent(in) :: a,tau,dvel_in(3),diff_const
real(dp),intent(out) :: dvel_out(3)
integer xyzzyaaaa44
real(dp) xyzzyaaab44,xyzzyaaac44,xyzzyaaad44
if(limdmc==0)then
dvel_out=dvel_in
elseif(limdmc==1)then
xyzzyaaab44=1.d0/tau
do xyzzyaaaa44=1,3
if(abs(dvel_in(xyzzyaaaa44))<=xyzzyaaab44)then
dvel_out(xyzzyaaaa44)=dvel_in(xyzzyaaaa44)
else
if(dvel_in(xyzzyaaaa44)>0.d0)then
dvel_out(xyzzyaaaa44)=xyzzyaaab44
else
dvel_out(xyzzyaaaa44)=-xyzzyaaab44
endif
endif
enddo
else
xyzzyaaac44=a*2.0d0*diff_const*tau*(dvel_in(1)*dvel_in(1)+ dvel_in(2)*&
&dvel_in(2)+dvel_in(3)*dvel_in(3))
if(xyzzyaaac44>0.d0)then
xyzzyaaad44=(-1.d0+sqrt(1.d0+2.d0*xyzzyaaac44))/xyzzyaaac44
dvel_out=xyzzyaaad44*dvel_in
else
dvel_out=dvel_in
endif
endif
end subroutine xyzzyaado1
subroutine xyzzyaadp1(dmceest,inv_sqrt_tau,etrial,dvelsq,dvellimsq,eco&
&nfig_in,econfig_out)
implicit none
real(dp),intent(in) :: dmceest,etrial,inv_sqrt_tau,econfig_in,dvelsq,d&
&vellimsq
real(dp),intent(out) :: econfig_out
real(dp) xyzzyaaaa45,xyzzyaaab45
if(limdmc==0)then
econfig_out=econfig_in
elseif(limdmc==1.or.limdmc==4)then
if(limdmc==1)then
xyzzyaaaa45=2.d0*inv_sqrt_tau
else
xyzzyaaaa45=inv_sqrt_tau*sqrt_netot*0.2d0
endif
if(abs(econfig_in-etrial)>xyzzyaaaa45)then
if(econfig_in>etrial)then
econfig_out=etrial+xyzzyaaaa45
else
econfig_out=etrial-xyzzyaaaa45
endif
else
econfig_out=econfig_in
endif
else
if(dvelsq>0.d0)then
xyzzyaaab45=sqrt(dvellimsq/dvelsq)
else
xyzzyaaab45=1.d0
endif
if(limdmc==2)then
econfig_out=dmceest-(dmceest-econfig_in)*xyzzyaaab45
else
econfig_out=xyzzyaaab45*econfig_in
endif
endif
end subroutine xyzzyaadp1
subroutine xyzzyaadq1(dt)
use slaarnaag, only : twopi
implicit none
real(dp),intent(in) :: dt
if(.not.have_ae)call errstop('INIT_NUCLEUS_MODS','NUCLEUS_GF_MODS is T&
& but HAVE_AE is F.')
if(model_system)call errstop('INIT_NUCLEUS_MODS','NUCLEUS_GF_MODS is T&
& but this is a model system.')
if(.not.any(is_ae(1:nitot)))call errstop('INIT_NUCLEUS_MODS','NUCLEUS_&
&GF_MODS is T but there are no all-electron atoms.')
if(dimensionality/=3)call errstop('INIT_NUCLEUS_MODS','NUCLEUS_GF_MODS&
& is T but the dimensionality is not 3.')
xyzzyaaah1=1.d0/sqrt(2.d0*dt)
xyzzyaaai1=(twopi*dt)**(-1.5d0)
xyzzyaaaj1=1.d0/(2.d0*dt)
xyzzyaaak1=1.d0/dt
xyzzyaaal1=sqrt(dt)
end subroutine xyzzyaadq1
subroutine xyzzyaadr1(vunlim,eivecs1,a,z,mag_z,z_hat,zion_closest)
implicit none
real(dp),intent(in) :: vunlim(3),eivecs1(4,nitot)
real(dp),intent(out) :: a,z(3),mag_z,z_hat(3),zion_closest
integer xyzzyaaaa47,xyzzyaaab47
real(dp) xyzzyaaac47,xyzzyaaad47(3)
xyzzyaaab47=-1
mag_z=huge(1.d0)
do xyzzyaaaa47=1,nitot
if(is_ae(xyzzyaaaa47).and.eivecs1(4,xyzzyaaaa47)<mag_z)then
xyzzyaaab47=xyzzyaaaa47
mag_z=eivecs1(4,xyzzyaaaa47)
endif
enddo
if(xyzzyaaab47==-1)call errstop('EVAL_ALIMIT','Can''t find any bare nu&
&clei.')
z(1:3)=eivecs1(1:3,xyzzyaaab47)
if(mag_z/=0.d0)then
z_hat(1:3)=z(1:3)/mag_z
else
z_hat(1:3)=0.d0
endif
zion_closest=zion(iontype(xyzzyaaab47))
xyzzyaaac47=sqrt(vunlim(1)**2+vunlim(2)**2+vunlim(3)**2)
if(xyzzyaaac47/=0.d0)then
xyzzyaaad47(1:3)=vunlim(1:3)/xyzzyaaac47
else
xyzzyaaad47(1:3)=0.d0
endif
a=0.5d0*(1.d0+dot_product(xyzzyaaad47,z_hat)) +zion_closest**2*mag_z**&
&2/(10.d0*(4.d0+zion_closest**2*mag_z**2))
end subroutine xyzzyaadr1
subroutine xyzzyaads1(rold,rnew,vold,g_variance,sqrt_g_variance,diff_d&
&isp)
implicit none
real(dp),intent(in) :: rold(3),vold(3),g_variance,sqrt_g_variance
real(dp),intent(out) :: rnew(3),diff_disp(3)
rnew=rold
diff_disp=0.d0
diff_disp(1)=ranx_gaussian(sqrt_g_variance)
rnew(1)=rold(1)+diff_disp(1)+g_variance*vold(1)
if(dimensionality>1)then
diff_disp(2)=ranx_gaussian(sqrt_g_variance)
rnew(2)=rnew(2)+diff_disp(2)+g_variance*vold(2)
if(dimensionality>2)then
diff_disp(3)=ranx_gaussian(sqrt_g_variance)
rnew(3)=rnew(3)+diff_disp(3)+g_variance*vold(3)
endif
endif
end subroutine xyzzyaads1
subroutine xyzzyaadt1(dt,vlim,z,mag_z,z_hat,zion_closest,r_old,r_new,d&
&iff_disp,greens_fn_forwards)
use slaarnabt,     only : erfc
implicit none
real(dp),intent(in) :: dt,vlim(3),z(3),mag_z,z_hat(3),r_old(3),zion_cl&
&osest
real(dp),intent(out) :: r_new(3),diff_disp(3),greens_fn_forwards
real(dp) xyzzyaaaa49,xyzzyaaab49(3),xyzzyaaac49,xyzzyaaad49(3),xyzzyaa&
&ae49,xyzzyaaaf49,xyzzyaaag49(3),xyzzyaaah49,xyzzyaaai49(3),xyzzyaaaj4&
&9,xyzzyaaak49
xyzzyaaaj49=sqrt(zion_closest**2+xyzzyaaak1)
xyzzyaaaa49=dot_product(vlim,z_hat)
xyzzyaaab49(1:3)=vlim(1:3)-xyzzyaaaa49*z_hat(1:3)
xyzzyaaac49=sqrt(xyzzyaaab49(1)**2+xyzzyaaab49(2)**2+xyzzyaaab49(3)**2&
&)
if(xyzzyaaac49/=0.d0)then
xyzzyaaad49(1:3)=xyzzyaaab49(1:3)/xyzzyaaac49
else
xyzzyaaad49(1:3)=0.d0
endif
xyzzyaaae49=max(mag_z+xyzzyaaaa49*dt,0.d0)
if(mag_z+xyzzyaaae49/=0.d0)then
xyzzyaaaf49=2.d0*xyzzyaaac49*dt*xyzzyaaae49/(mag_z+xyzzyaaae49)
else
xyzzyaaaf49=0.d0
endif
xyzzyaaag49(1:3)=xyzzyaaae49*z_hat(1:3)+xyzzyaaaf49*xyzzyaaad49(1:3)
r_new(1:3)=r_old(1:3)+xyzzyaaag49(1:3)-z(1:3)
xyzzyaaak49=0.5d0*erfc((mag_z+xyzzyaaaa49*dt)*xyzzyaaah1)
xyzzyaaah49=ranx()
if(xyzzyaaah49<xyzzyaaak49)then
call xyzzyaadu1(xyzzyaaaj49,xyzzyaaai49)
diff_disp(1:3)=xyzzyaaai49(1:3)-xyzzyaaag49(1:3)
else
diff_disp(1)=ranx_gaussian(xyzzyaaal1)
diff_disp(2)=ranx_gaussian(xyzzyaaal1)
diff_disp(3)=ranx_gaussian(xyzzyaaal1)
xyzzyaaai49(1:3)=xyzzyaaag49(1:3)+diff_disp(1:3)
endif
r_new(1:3)=r_new(1:3)+diff_disp(1:3)
greens_fn_forwards=xyzzyaadx1(xyzzyaaaj49,xyzzyaaak49,diff_disp,xyzzya&
&aai49)
end subroutine xyzzyaadt1
subroutine xyzzyaadu1(zeta,xi)
use slaarnaag, only : twopi
implicit none
real(dp),intent(in) :: zeta
real(dp),intent(out) :: xi(3)
real(dp) xyzzyaaaa50,xyzzyaaab50,xyzzyaaac50,xyzzyaaad50,xyzzyaaae50,x&
&yzzyaaaf50,xyzzyaaag50
xyzzyaaaa50=ranx()
xyzzyaaab50=ranx()
xyzzyaaac50=ranx()
xyzzyaaad50=-log(xyzzyaaaa50*xyzzyaaab50*xyzzyaaac50)/(2.d0*zeta)
xyzzyaaae50=1.d0-2.d0*ranx()
xyzzyaaaf50=ranx()*twopi
xyzzyaaag50=xyzzyaaad50*sqrt(1.d0-xyzzyaaae50**2)
xi(1)=xyzzyaaag50*cos(xyzzyaaaf50)
xi(2)=xyzzyaaag50*sin(xyzzyaaaf50)
xi(3)=xyzzyaaad50*xyzzyaaae50
end subroutine xyzzyaadu1
real(dp) function xyzzyaadv1(xi)
implicit none
real(dp),intent(in) :: xi(3)
xyzzyaadv1=xyzzyaaai1*exp(-(xi(1)**2+xi(2)**2+xi(3)**2)*xyzzyaaaj1)
end function xyzzyaadv1
real(dp) function xyzzyaadw1(zeta,xi)
use slaarnaag, only : one_over_pi
implicit none
real(dp),intent(in) :: zeta,xi(3)
xyzzyaadw1=one_over_pi*zeta**3*exp(-2.d0*zeta*sqrt(xi(1)**2+xi(2)**2+x&
&i(3)**2))
end function xyzzyaadw1
real(dp) function xyzzyaadx1(zeta,q_tilde,diff_disp,final_z)
implicit none
real(dp),intent(in) :: zeta,q_tilde,diff_disp(3),final_z(3)
xyzzyaadx1=(1.d0-q_tilde)*xyzzyaadv1(diff_disp)+q_tilde*xyzzyaadw1(zet&
&a,final_z)
end function xyzzyaadx1
real(dp) function xyzzyaady1(dt,vlim,z,mag_z,z_hat,zion_closest,r_old,&
&r_new,greens_fn_forwards)
use slaarnabt, only : erfc
implicit none
real(dp),intent(in) :: dt,vlim(3),z(3),mag_z,z_hat(3),r_old(3),r_new(3&
&),greens_fn_forwards,zion_closest
real(dp) xyzzyaaaa54,xyzzyaaab54(3),xyzzyaaac54,xyzzyaaad54(3),xyzzyaa&
&ae54,xyzzyaaaf54,xyzzyaaag54(3),xyzzyaaah54(3),xyzzyaaai54,xyzzyaaaj5&
&4(3),xyzzyaaak54(3),xyzzyaaal54,xyzzyaaam54
xyzzyaaai54=sqrt(zion_closest**2+xyzzyaaak1)
xyzzyaaaa54=dot_product(vlim,z_hat)
xyzzyaaab54(1:3)=vlim(1:3)-xyzzyaaaa54*z_hat(1:3)
xyzzyaaac54=sqrt(xyzzyaaab54(1)**2+xyzzyaaab54(2)**2+xyzzyaaab54(3)**2&
&)
if(xyzzyaaac54/=0.d0)then
xyzzyaaad54(1:3)=xyzzyaaab54(1:3)/xyzzyaaac54
else
xyzzyaaad54(1:3)=0.d0
endif
xyzzyaaae54=max(mag_z+xyzzyaaaa54*dt,0.d0)
if(mag_z+xyzzyaaae54/=0.d0)then
xyzzyaaaf54=2.d0*xyzzyaaac54*dt*xyzzyaaae54/(mag_z+xyzzyaaae54)
else
xyzzyaaaf54=0.d0
endif
xyzzyaaag54(1:3)=xyzzyaaae54*z_hat(1:3)+xyzzyaaaf54*xyzzyaaad54(1:3)
xyzzyaaah54(1:3)=r_new(1:3)+xyzzyaaag54(1:3)-z(1:3)
xyzzyaaaj54(1:3)=r_old(1:3)-xyzzyaaah54(1:3)
xyzzyaaak54(1:3)=xyzzyaaaj54(1:3)+xyzzyaaag54(1:3)
xyzzyaaam54=0.5d0*erfc((mag_z+xyzzyaaaa54*dt)*xyzzyaaah1)
xyzzyaaal54=xyzzyaadx1(xyzzyaaai54,xyzzyaaam54,xyzzyaaaj54,xyzzyaaak54&
&)
if(greens_fn_forwards/=0.d0)then
xyzzyaady1=xyzzyaaal54/greens_fn_forwards
else
xyzzyaady1=0.d0
endif
end function xyzzyaady1
subroutine dmc_main(dt,nmove_in,nblock_in,hf_ke,hf_ex,aborted)
use slaarnaaf
use slaarnaag
use slaarnaam
use slaarnabj
use slaarnacd
use slaarnaab, only : forces_scratch_request,setup_forces_accum,finish&
&_forces_accum,eval_local_forces,forces_to_future,future_to_forces,tag&
&h_forces,nfterms,nfcomps
use slaarnaad,      only : backflow_stats
use slaarnaaq,        only : expval_scratch_request,setup_expval_accum&
&,finish_expval_accum,accumulate_expvals,expval_zero_postweight,expval&
&_postweight,write_expval,expvals,particle_is_fixed,fixed_particle,eva&
&l_dipole_moment,int_sf,eval_int_sf,finite_size_corr_xc,pair_corr,pair&
&_corr_sph,pcfs_rfix,pcf_rfix,read_expval,backup_expval_file,setup_exp&
&val,deallocate_expval,zero_expval
use file_utils,    only : open_units
use slaarnabk,  only : interaction_mpc_present,interaction_mpc_use,mpc&
&_hartree
use slaarnabs,     only : tmove_no_points,tmove_t,tmove_fixgrid,tmove_&
&points,v_non_local_tmove
use slaarnabt,     only : dscal,ddot,exp_protect
use slaarnace,    only : relativistic
implicit none
integer,intent(in) :: nmove_in,nblock_in
real(dp),intent(in) :: dt,hf_ke,hf_ex
logical,intent(out) :: aborted
integer xyzzyaaaa55,xyzzyaaab55,xyzzyaaac55
real(dp) xyzzyaaad55,xyzzyaaae55,xyzzyaaaf55,xyzzyaaag55,xyzzyaaah55,x&
&yzzyaaai55,xyzzyaaaj55,xyzzyaaak55,xyzzyaaal55
integer iblock,xyzzyaaam55,xyzzyaaan55,xyzzyaaao55,xyzzyaaap55,xyzzyaa&
&aq55,xyzzyaaar55,xyzzyaaas55,xyzzyaaat55,xyzzyaaau55,xyzzyaaav55,xyzz&
&yaaaw55,xyzzyaaax55,xyzzyaaay55,xyzzyaaaz55,xyzzyaaba55,xyzzyaabb55,x&
&yzzyaabc55
integer xyzzyaabd55,xyzzyaabe55,xyzzyaabf55,xyzzyaabg55,xyzzyaabh55,xy&
&zzyaabi55,xyzzyaabj55,xyzzyaabk55,xyzzyaabl55,xyzzyaabm55,xyzzyaabn55&
&(netot)
integer,allocatable :: xyzzyaabo55(:)
real(dp) xyzzyaabp55(3),xyzzyaabq55,prob,relprob,xyzzyaabr55,xyzzyaabs&
&55(3),xyzzyaabt55(4),xyzzyaabu55,xyzzyaabv55(2),xyzzyaabw55,xyzzyaabx&
&55,xyzzyaaby55,xyzzyaabz55,xyzzyaaca55,xyzzyaacb55,xyzzyaacc55,xyzzya&
&acd55(2),xyzzyaace55
real(dp),allocatable :: xyzzyaacf55(:,:),xyzzyaacg55(:,:),xyzzyaach55(&
&:)
logical xyzzyaaci55,xyzzyaacj55,xyzzyaack55
real(dp) xyzzyaacl55(3),xyzzyaacm55(3),xyzzyaacn55(3),xyzzyaaco55(3),x&
&yzzyaacp55(3),xyzzyaacq55,xyzzyaacr55(3),xyzzyaacs55,xyzzyaact55,xyzz&
&yaacu55,xyzzyaacv55,xyzzyaacw55,xyzzyaacx55,xyzzyaacy55,xyzzyaacz55,x&
&yzzyaada55(nspin),xyzzyaadb55(nspin),xyzzyaadc55,xyzzyaadd55(nspin)
real(dp),allocatable :: xyzzyaade55(:),xyzzyaadf55(:,:)
complex(dp) xyzzyaadg55,xyzzyaadh55(3),xyzzyaadi55
integer xyzzyaadj55,xyzzyaadk55,xyzzyaadl55
integer,allocatable :: xyzzyaadm55(:)
real(dp) wdmc,xyzzyaadn55,xyzzyaado55,xyzzyaadp55,xyzzyaadq55
real(dp),parameter :: xyzzyaadr55=5.d0
real(dp) xyzzyaads55(3),xyzzyaadt55,xyzzyaadu55(3),xyzzyaadv55,xyzzyaa&
&dw55
real(dp),allocatable :: xyzzyaadx55(:)
real(dp) xyzzyaady55,xyzzyaadz55,xyzzyaaea55,xyzzyaaeb55,xyzzyaaec55,x&
&yzzyaaed55,xyzzyaaee55,xyzzyaaef55,xyzzyaaeg55
real(dp),allocatable :: xyzzyaaeh55(:),xyzzyaaei55(:)
logical xyzzyaaej55
real(dp) xyzzyaaek55,xyzzyaael55,xyzzyaaem55
integer xyzzyaaen55
logical xyzzyaaeo55,xyzzyaaep55
integer xyzzyaaeq55,xyzzyaaer55,xyzzyaaes55
real(dp) xyzzyaaet55,xyzzyaaeu55,xyzzyaaev55,xyzzyaaew55,xyzzyaaex55,x&
&yzzyaaey55,xyzzyaaez55,xyzzyaafa55,xyzzyaafb55,xyzzyaafc55,xyzzyaafd5&
&5,xyzzyaafe55,xyzzyaaff55,xyzzyaafg55(13),xyzzyaafh55,xyzzyaafi55,xyz&
&zyaafj55
real(dp),allocatable :: xyzzyaafk55(:),xyzzyaafl55(:),xyzzyaafm55(:),x&
&yzzyaafn55(:),xyzzyaafo55(:),xyzzyaafp55(:),xyzzyaafq55(:,:)
logical isnan,isinf,xyzzyaafr55,xyzzyaafs55
integer xyzzyaaft55,xyzzyaafu55
integer xyzzyaafv55,xyzzyaafw55,xyzzyaafx55
integer xyzzyaafy55
real(dp) xyzzyaafz55,xyzzyaaga55,xyzzyaagb55,xyzzyaagc55
type(configuration),pointer :: xyzzyaagd55,xyzzyaage55,xyzzyaagf55
logical,save :: xyzzyaagg55=.false.
real(sp) xyzzyaagh55,xyzzyaagi55
integer xyzzyaagj55,xyzzyaagk55,xyzzyaagl55,xyzzyaagm55,xyzzyaagn55
integer,parameter :: xyzzyaago55=10
real(dp),parameter :: xyzzyaagp55=10.d0
real(dp),allocatable :: xyzzyaagq55(:),xyzzyaagr55(:),xyzzyaags55(:,:)&
&,xyzzyaagt55(:,:),xyzzyaagu55(:,:)
integer xyzzyaagv55,xyzzyaagw55,xyzzyaagx55
logical,parameter :: xyzzyaagy55=.false.
logical,parameter :: xyzzyaagz55=.false.
aborted=.false.
call timer('DMC',.true.)
if(am_master)then
call wout()
call wout('BEGIN DMC CALCULATION')
call wout('=====================')
call wout()
endif
call timer('SETUP',.true.,collapse=.true.)
call xyzzyaaha55
call xyzzyaade1
xyzzyaaeq55=0
xyzzyaaer55=0
xyzzyaaes55=0
call scratch_protect(xyzzyaaeq55)
if(dmc_cfg_by_cfg)then
call scratch_request(ratiocfg_from=xyzzyaaeq55,ratiocfg_to=xyzzyaaer55&
&)
call energy_scratch_request(xyzzyaaer55)
else
call scratch_request(ratio1_from=xyzzyaaeq55,ratio1_to=xyzzyaaer55)
endif
if(use_tmove)call scratch_request(ratio1_from=xyzzyaaeq55,ratio1_to=xy&
&zzyaaes55)
call scratch_request(drift=xyzzyaaeq55)
call scratch_request(drift=xyzzyaaer55)
call energy_scratch_request(xyzzyaaeq55)
if(expvals)call expval_scratch_request(xyzzyaaeq55)
if(forces)call forces_scratch_request(xyzzyaaeq55)
call setup_scratch
call which_scratch(xyzzyaaeq55)
call which_scratch(xyzzyaaer55)
if(use_tmove)call which_scratch(xyzzyaaes55)
call setup_wfn_utils
call setup_energy_utils
if(expvals)call setup_expval_accum
if(forces)call setup_forces_accum
call xyzzyaahb55
if(am_master.and.iaccum)call init_reblock(no_cols_qmc)
call timer('SETUP',.false.)
call xyzzyaaby1(use_future,xyzzyaagk55,xyzzyaagj55)
call xyzzyaahc55(.false.)
if(am_master.and..not.iaccum)xyzzyaafo55=xyzzyaafa55
if(.not.newrun.and.trim(random_seed_kw)=='timer_reset')then
if(am_master)then
call wordwrap("Re-initializing random number sequence on restart as RA&
&NDOM_SEED = 'timer_reset'.")
call wout()
endif
call initialize_random(ranluxlevel,ranprint,ranlog_unit)
endif
if(nconfig_prelim>0.and.iaccum)then
if(am_master.and.xyzzyaaav55>0)call errstop('INITIALIZE_DMC','Cannot c&
&ontinue the statistics-accumulation phase of preliminary DMC.')
call xyzzyaahh55
endif
xyzzyaaen55=0
if(use_backflow)call backflow_stats(xyzzyaabx55,xyzzyaaby55)
xyzzyaabh1=.false.
xyzzyaabi1=.false.
xyzzyaaar55=xyzzyaaaq55
iblock=1
block_loop: do
call timer('DMC_BLOCK',.true.)
if(am_master)then
xyzzyaagh55=tcputime()
xyzzyaagi55=0.
endif
if(exceeds_time_limit(iblock==1))then
if(am_master)then
call wout()
call errwarn('DMC','Time limit exceeded or about to be exceeded. Emerg&
&ency stop.')
if(am_master.and.chkpoint_level==-1.and.xyzzyaagg55)then
call wout('Writing config.out file for restart despite CHECKPOINT==-1.&
&')
call wout()
endif
call wout('CONTINUATION INFO:')
if(chkpoint_level==-1)then
call wout('Run cannot be continued since CHECKPOINT=-1 in input.')
elseif(iaccum.and.nconfig_prelim>0)then
call wout('Cannot continue the "statistics accumulation" stage of prel&
&iminary DMC.')
elseif(isdmc.and.iblock==1)then
call wout(' Suggested action: restart with more time')
call wout(' PROBLEM: first DMC block won''t fit in time slot!')
call wout(' You must re-run using a longer time slot.')
else
call wout(' Suggested action: continue run directly')
if(isvmc_dmc.and..not.iaccum)then
call wout(' Set RUNTYPE          : dmc')
elseif(iaccum.and..not.isdmc)then
call wout(' Set RUNTYPE          : dmc_stats')
endif
call wout(' Set NEWRUN           : F')
if(.not.iaccum)then
if(old_input)then
if(use_blocktime)then
call wout(' Set NMOVE_DMC_EQUIL : '//trim(i2s(xyzzyaaar55)))
else
call wout(' Set NMOVE_DMC_EQUIL : '//trim(i2s(xyzzyaaam55-iblock+1)))
endif
else
if(use_blocktime)then
call wout(' Set DMC_EQUIL_NSTEP  : '//trim(i2s(xyzzyaaar55)))
else
call wout(' Set DMC_EQUIL_NBLOCK : '//trim(i2s(xyzzyaaam55-iblock+1)))
call wout(' Set DMC_EQUIL_NSTEP  : '//trim(i2s((xyzzyaaam55-iblock+1)*&
&xyzzyaaaq55)))
endif
endif
else
if(old_input)then
if(use_blocktime)then
call wout(' Set NMOVE_DMC_STATS : '//trim(i2s(xyzzyaaar55)))
else
call wout(' Set NMOVE_DMC_STATS : '//trim(i2s(xyzzyaaam55-iblock+1)))
endif
else
if(use_blocktime)then
call wout(' Set DMC_STATS_NSTEP  : '//trim(i2s(xyzzyaaar55)))
else
call wout(' Set DMC_STATS_NBLOCK : '//trim(i2s(xyzzyaaam55-iblock+1)))
call wout(' Set DMC_STATS_NSTEP  : '//trim(i2s((xyzzyaaam55-iblock+1)*&
&xyzzyaaaq55)))
endif
endif
call wout(' Move config.out to config.in')
call wout()
call wout(' NB: The runqmc options --continue/--auto-continue will do &
&this automatically.')
endif
endif
endif
if(xyzzyaagg55)then
if(xyzzyaabh1.or.xyzzyaabi1)call xyzzyaadd1('RECV',xyzzyaaax55)
call xyzzyaahg55(.false.,.true.)
endif
aborted=.true.
call timer('DMC_BLOCK',.false.)
exit block_loop
endif
xyzzyaabj55=0
xyzzyaabh55=0
xyzzyaadl55=0
xyzzyaadp55=0.d0
xyzzyaabw55=0
xyzzyaabl55=0
if(am_master.and.writeout_dmc_hist)xyzzyaafq55=0.d0
xyzzyaaeo55=.false.
xyzzyaacb55=0.d0
xyzzyaack55=.false.
xyzzyaaan55=0
move_loop: do xyzzyaaao55=1,xyzzyaaar55
dmcave_loop: do xyzzyaaap55=1,ndmcave
call timer('DMC_STEP',.true.)
xyzzyaaan55=xyzzyaaan55+1
if(real(xyzzyaaau55,dp)+real(xyzzyaaav55,dp)+real(xyzzyaaan55-1,dp)>ma&
&x_rep_int)call errstop('DMC_MAIN','Integer overflow: ipass.')
xyzzyaabx1=xyzzyaaau55+xyzzyaaav55+xyzzyaaan55-1
xyzzyaacj55=iaccum.and.(mod(xyzzyaabx1-xyzzyaaau55,corper_dmc)==0)
if(use_future)xyzzyaagl55=mod(xyzzyaagl55,xyzzyaagj55)+1
xyzzyaafm55(:)=0.d0
if(use_future)then
xyzzyaags55=0.d0
if(xyzzyaabk1>1)xyzzyaagt55=0.d0
endif
xyzzyaabk55=0
xyzzyaabu55=0
call timer('MOVE_CONFIGS',.true.)
xyzzyaabc55=0
if(xyzzyaaax55>0)then
xyzzyaagd55=>xyzzyaaar1
if(expvals.and.iaccum.and.tpdmc>0)call expval_zero_postweight
config_loop: do
xyzzyaabc55=xyzzyaabc55+1
if(xyzzyaagy55.or.xyzzyaagz55)call wout('CONFIGURATION #'//trim(i2s(xy&
&zzyaabc55))//':')
call xyzzyaahd55
xyzzyaagd55=>xyzzyaagd55%pt_next_config
if(associated(xyzzyaagd55,xyzzyaaar1))then
if(xyzzyaabh1)xyzzyaagf55=>xyzzyaagd55%pt_prev_config
nullify(xyzzyaagd55)
exit config_loop
endif
enddo config_loop
endif
call timer('MOVE_CONFIGS',.false.)
if(xyzzyaabh1.or.xyzzyaabi1)then
xyzzyaafr55=(xyzzyaaax55==0)
call xyzzyaadd1('RECV',xyzzyaaax55)
if(xyzzyaabh1)then
call timer('MOVE_CONFIGS',.true.)
if(xyzzyaafr55)then
xyzzyaagd55=>xyzzyaaar1
else
xyzzyaagd55=>xyzzyaagf55
endif
config_loop2: do
xyzzyaabc55=xyzzyaabc55+1
if(xyzzyaagy55.or.xyzzyaagz55)call  wout('CONFIGURATION #'//trim(i2s(x&
&yzzyaabc55))//':')
if(.not.xyzzyaafr55)then
xyzzyaagd55=>xyzzyaagd55%pt_next_config
if(associated(xyzzyaagd55,xyzzyaaar1))then
nullify(xyzzyaagd55)
exit config_loop2
endif
endif
xyzzyaafr55=.false.
call xyzzyaahd55
enddo config_loop2
call timer('MOVE_CONFIGS',.false.)
endif
endif
call timer('GET_ACCUM_ENERGIES',.true.)
call qmpi_reduce(xyzzyaafm55,xyzzyaafl55,mpi_sum,'DMC_MAIN','sum_data'&
&)
if(use_future)then
call qmpi_reduce(xyzzyaags55,xyzzyaagu55,mpi_sum,'DMC_MAIN','pureave')
endif
call timer('GET_ACCUM_ENERGIES',.false.)
call timer('MPI_NCONFIG_AGEINFO',.true.)
xyzzyaabv55=(/0.d0,real(xyzzyaabk55,dp)/)
if(xyzzyaaax55>0)xyzzyaabv55(1)=real(xyzzyaabu55,dp)/real(xyzzyaaax55,&
&dp)
if(xyzzyaabt1)then
xyzzyaaay55=xyzzyaaax55
xyzzyaaba55=xyzzyaaax55
do xyzzyaaaa55=0,xyzzyaabp1(xyzzyaabm1)-1
xyzzyaaac55=xyzzyaabq1(xyzzyaabm1)+xyzzyaaaa55
if(xyzzyaaac55>nnodes-1)xyzzyaaac55=xyzzyaaac55-nnodes
if(xyzzyaaac55==xyzzyaabl1)cycle
call qmpi_recv(xyzzyaabs55,xyzzyaaac55,aveid,'DMC_MAIN','ncon_age')
xyzzyaaaz55=int(xyzzyaabs55(1))
xyzzyaaay55=xyzzyaaay55+xyzzyaaaz55
if(xyzzyaaaz55>xyzzyaaba55)xyzzyaaba55=xyzzyaaaz55
xyzzyaabv55(1)=xyzzyaabv55(1)+xyzzyaabs55(2)
if(xyzzyaabs55(3)>xyzzyaabv55(2))xyzzyaabv55(2)=xyzzyaabs55(3)
enddo
xyzzyaabt55(1)=real(xyzzyaaay55,dp)
xyzzyaabt55(2)=real(xyzzyaaba55,dp)
xyzzyaabt55(3:4)=xyzzyaabv55(1:2)
if(am_master)then
allocate(xyzzyaach55(4*xyzzyaabk1),stat=xyzzyaabw1)
else
allocate(xyzzyaach55(1),stat=xyzzyaabw1)
endif
call mpi_gather(xyzzyaabt55,4,mpi_double_precision,xyzzyaach55,4,mpi_d&
&ouble_precision,0,rg_comm(xyzzyaabk1+1,1),ierror)
if(am_master)then
xyzzyaaab55=5
do xyzzyaaaa55=2,xyzzyaabk1
xyzzyaaaz55=int(xyzzyaach55(xyzzyaaab55))
xyzzyaaay55=xyzzyaaay55+xyzzyaaaz55
xyzzyaabb55=int(xyzzyaach55(xyzzyaaab55+1))
if(xyzzyaabb55>xyzzyaaba55)xyzzyaaba55=xyzzyaabb55
xyzzyaabv55(1)=xyzzyaabv55(1)+xyzzyaach55(xyzzyaaab55+2)
if(xyzzyaach55(xyzzyaaab55+3)>xyzzyaabv55(2))xyzzyaabv55(2)=xyzzyaach5&
&5(xyzzyaaab55+3)
xyzzyaaab55=xyzzyaaab55+4
enddo
endif
deallocate(xyzzyaach55)
else
xyzzyaabs55(1)=real(xyzzyaaax55,dp)
xyzzyaabs55(2:3)=xyzzyaabv55(1:2)
call qmpi_ssend(xyzzyaabs55,xyzzyaabl1,aveid,'DMC_MAIN','ncon_age')
endif
call timer('MPI_NCONFIG_AGEINFO',.false.)
if(use_blocktime.and..not.xyzzyaack55)then
call xyzzyaahe55
if(xyzzyaack55.or.xyzzyaaao55==xyzzyaaar55)xyzzyaaaq55=xyzzyaaao55
endif
if(am_master)then
call timer('DMC_AVERAGING',.true.)
xyzzyaadn55=xyzzyaafl55(xyzzyaaft55)
if(xyzzyaadn55<=0.d0)call errstop('DMC','Total weight is zero.')
xyzzyaado55=1.d0/xyzzyaadn55
xyzzyaafk55(1:xyzzyaafv55)=xyzzyaafl55(1:xyzzyaafv55)*xyzzyaado55
if(use_future)xyzzyaagu55=xyzzyaagu55*xyzzyaado55
if(.not.xyzzyaaej55)then
xyzzyaafc55=xyzzyaafb55
xyzzyaael55=xyzzyaafk55(xyzzyaafu55)
xyzzyaaem55=xyzzyaael55*xyzzyaafb55
if(growth_estimator)xyzzyaadz55=-xyzzyaaem55
endif
if(tpdmc>0.or.growth_estimator)xyzzyaaed55=xyzzyaaem55-xyzzyaael55*xyz&
&zyaaet55
if(tpdmc>0)then
xyzzyaaas55=mod(xyzzyaabx1,tpdmc)
xyzzyaady55=xyzzyaady55+(xyzzyaaed55-xyzzyaaeh55(xyzzyaaas55))
xyzzyaaee55=exp_protect(xyzzyaady55)
xyzzyaaeh55(xyzzyaaas55)=xyzzyaaed55
endif
if(growth_estimator)then
xyzzyaaat55=mod(xyzzyaabx1,tpdmc+1)
xyzzyaadz55=xyzzyaadz55+ (xyzzyaaed55-xyzzyaaei55(xyzzyaaat55))
xyzzyaaef55=exp_protect(xyzzyaadz55)
xyzzyaaei55(xyzzyaaat55)=xyzzyaaed55
endif
if(iaccum)then
if(growth_estimator)then
xyzzyaaeb55=xyzzyaaea55
if(xyzzyaabx1>xyzzyaaau55)xyzzyaaec55=xyzzyaaec55+xyzzyaaef55*xyzzyaad&
&n55
endif
xyzzyaaea55=xyzzyaaea55+xyzzyaaee55*xyzzyaadn55
if(interaction_mpc_use)then
xyzzyaafg55(1)=xyzzyaafg55(1)+xyzzyaaee55*xyzzyaafl55(i_eloc_mpc)
else
xyzzyaafg55(1)=xyzzyaafg55(1)+xyzzyaaee55*xyzzyaafl55(i_eloc_def)
endif
xyzzyaafg55(2)=xyzzyaafg55(2)+xyzzyaaee55*xyzzyaafl55(i_ti)
xyzzyaafg55(3)=xyzzyaafg55(3)+xyzzyaaee55*xyzzyaafl55(i_kei)
xyzzyaafg55(4)=xyzzyaafg55(4)+xyzzyaaee55*xyzzyaafl55(i_fisq)
xyzzyaafg55(5)=xyzzyaafg55(5)+xyzzyaaee55*xyzzyaafl55(i_pote)
if(interaction_mpc_present)then
xyzzyaafg55(6)=xyzzyaafg55(6)+xyzzyaaee55*xyzzyaafl55(i_short)
xyzzyaafg55(7)=xyzzyaafg55(7)+xyzzyaaee55*xyzzyaafl55(i_long)
endif
xyzzyaafg55(8)=xyzzyaafg55(8)+xyzzyaaee55*xyzzyaafl55(i_potinl)
xyzzyaafg55(9)=xyzzyaafg55(9)+xyzzyaaee55*xyzzyaafl55(i_potil)
if(have_veep)then
if(isperiodic)then
xyzzyaafg55(10)=xyzzyaafg55(10)+xyzzyaaee55*xyzzyaafl55(i_ecpp_tot)
xyzzyaafg55(11:12)=0.d0
else
xyzzyaafg55(10)=xyzzyaafg55(10)+xyzzyaaee55*xyzzyaafl55(i_vcpp_ei)
xyzzyaafg55(11)=xyzzyaafg55(11)+xyzzyaaee55*xyzzyaafl55(i_vcpp_e)
xyzzyaafg55(12)=xyzzyaafg55(12)+xyzzyaaee55*xyzzyaafl55(i_vcpp_ee)
endif
endif
xyzzyaafg55(13)=xyzzyaafg55(13)+xyzzyaaee55*xyzzyaafl55(xyzzyaafu55)
if(expvals.and.tpdmc>0.and.xyzzyaacj55)call expval_postweight(xyzzyaae&
&e55)
xyzzyaaeg55=1.d0/xyzzyaaea55
xyzzyaafa55=xyzzyaafg55(1)*xyzzyaaeg55
xyzzyaaae55=real(xyzzyaabx1-xyzzyaaau55+1,dp)
if(xyzzyaaae55<xyzzyaaaf55)xyzzyaafa55=(xyzzyaafa55*xyzzyaaae55+xyzzya&
&afb55*(xyzzyaaaf55-xyzzyaaae55))/xyzzyaaaf55
xyzzyaael55=xyzzyaafg55(13)*xyzzyaaeg55
else
xyzzyaaaw55=1+mod(xyzzyaabx1,ebest_av_window)
if(interaction_mpc_use)then
xyzzyaafj55=xyzzyaafk55(i_eloc_mpc)
else
xyzzyaafj55=xyzzyaafk55(i_eloc_def)
endif
if(xyzzyaaao55==xyzzyaaaq55.and.xyzzyaaap55==ndmcave)then
xyzzyaafo55(xyzzyaaaw55)=xyzzyaafj55
xyzzyaafa55=sum(xyzzyaafo55)/xyzzyaaaf55
else
xyzzyaafa55=xyzzyaafa55+(xyzzyaafj55-xyzzyaafo55(xyzzyaaaw55))/xyzzyaa&
&af55
xyzzyaafo55(xyzzyaaaw55)=xyzzyaafj55
endif
xyzzyaael55=xyzzyaafk55(xyzzyaafu55)
endif
if(growth_estimator.and.xyzzyaaej55.and.xyzzyaabx1>xyzzyaaau55.and.iac&
&cum)xyzzyaafc55=(-1.d0/xyzzyaael55)*log(xyzzyaaec55/xyzzyaaeb55)
xyzzyaadp55=xyzzyaadp55+real(xyzzyaaay55,dp)*xyzzyaaag55/real(xyzzyaab&
&a55,dp)
xyzzyaabq55=0.d0
if(xyzzyaabh55>0)xyzzyaabq55=real(xyzzyaabj55,dp)/real(xyzzyaabh55,dp)
xyzzyaabv55(1)=xyzzyaabv55(1)*xyzzyaaag55
xyzzyaabw55=xyzzyaabw55+xyzzyaabv55(1)
if(xyzzyaabl55<nint(xyzzyaabv55(2)))xyzzyaabl55=nint(xyzzyaabv55(2))
xyzzyaafp55(:)=0.d0
xyzzyaafp55(tagh_step)=real(xyzzyaabx1+1,dp)
if(tagh_nconf>0)xyzzyaafp55(tagh_nconf)=real(xyzzyaaay55,dp)
if(tagh_eref>0)xyzzyaafp55(tagh_eref)=xyzzyaaet55*xyzzyaafe55
if(tagh_acc>0)xyzzyaafp55(tagh_acc)=xyzzyaabq55
if(tagh_ebest>0)xyzzyaafp55(tagh_ebest)=xyzzyaafa55*xyzzyaafe55
if(tagh_weight>0)xyzzyaafp55(tagh_weight)=xyzzyaaee55*xyzzyaadn55
if(tagh_energy>0)then
if(interaction_mpc_use)then
xyzzyaafp55(tagh_energy)=xyzzyaafk55(i_eloc_mpc)*xyzzyaafe55
if(tagh_etotalt>0)xyzzyaafp55(tagh_etotalt)=xyzzyaafk55(i_eloc_def)*xy&
&zzyaafe55
else
xyzzyaafp55(tagh_energy)=xyzzyaafk55(i_eloc_def)*xyzzyaafe55
if(tagh_etotalt>0)xyzzyaafp55(tagh_etotalt)=xyzzyaafk55(i_eloc_mpc)*xy&
&zzyaafe55
endif
endif
if(tagh_esqr>0)xyzzyaafp55(tagh_esqr)=xyzzyaafk55(i_esqr)*xyzzyaafe55*&
&*2
if(tagh_teff>0)xyzzyaafp55(tagh_teff)=xyzzyaafk55(xyzzyaafu55)
if(tagh_k>0)xyzzyaafp55(tagh_k)=xyzzyaafk55(i_kei)*xyzzyaafe55
if(tagh_t>0)xyzzyaafp55(tagh_t)=xyzzyaafk55(i_ti)*xyzzyaafe55
if(tagh_fisq>0)xyzzyaafp55(tagh_fisq)=xyzzyaafk55(i_fisq)*xyzzyaafe55
if(tagh_ewald>0)xyzzyaafp55(tagh_ewald)=xyzzyaafk55(i_pote)*xyzzyaafe5&
&5
if(tagh_local>0)xyzzyaafp55(tagh_local)=xyzzyaafk55(i_potil)*xyzzyaafe&
&55
if(tagh_nonlocal>0)xyzzyaafp55(tagh_nonlocal)=xyzzyaafk55(i_potinl)*xy&
&zzyaafe55
if(interaction_mpc_present)then
if(tagh_short>0)xyzzyaafp55(tagh_short)=xyzzyaafk55(i_short)*xyzzyaafe&
&55
if(tagh_long>0)xyzzyaafp55(tagh_long)=xyzzyaafk55(i_long)*xyzzyaafe55
endif
if(have_veep)then
if(tagh_cppei>0)then
if(isperiodic)then
xyzzyaafp55(tagh_cppei)=xyzzyaafk55(i_ecpp_tot)*xyzzyaafe55
else
xyzzyaafp55(tagh_cppei)=xyzzyaafk55(i_vcpp_ei)*xyzzyaafe55
endif
endif
if(tagh_cppe>0)xyzzyaafp55(tagh_cppe)=xyzzyaafk55(i_vcpp_e)*xyzzyaafe5&
&5
if(tagh_cppee>0)xyzzyaafp55(tagh_cppee)=xyzzyaafk55(i_vcpp_ee)*xyzzyaa&
&fe55
endif
if(relativistic)then
if(tagh_masspol>0)xyzzyaafp55(tagh_masspol)=xyzzyaafk55(i_emasspol)*xy&
&zzyaafe55
if(tagh_massvel>0)xyzzyaafp55(tagh_massvel)=xyzzyaafk55(i_emassvel)*xy&
&zzyaafe55
if(tagh_darwinen>0)xyzzyaafp55(tagh_darwinen)=xyzzyaafk55(i_edarwin_en&
&)*xyzzyaafe55
if(tagh_darwinee>0)xyzzyaafp55(tagh_darwinee)=xyzzyaafk55(i_edarwin_ee&
&)*xyzzyaafe55
if(tagh_retard>0)xyzzyaafp55(tagh_retard)=xyzzyaafk55(i_eretard)*xyzzy&
&aafe55
endif
if(eval_dipole_moment)then
if(tagh_dipole1>0)xyzzyaafp55(tagh_dipole1)=xyzzyaafk55(i_dipole1)
if(tagh_dipole2>0)xyzzyaafp55(tagh_dipole2)=xyzzyaafk55(i_dipole2)
if(tagh_dipole3>0)xyzzyaafp55(tagh_dipole3)=xyzzyaafk55(i_dipole3)
if(tagh_dipole_sq>0)xyzzyaafp55(tagh_dipole_sq)=xyzzyaafk55(i_dipole_s&
&q)
endif
if(forces)then
xyzzyaagx55=xyzzyaafw55-1
do xyzzyaabf55=1,nitot_forces
do xyzzyaagv55=1,naxis_forces
do xyzzyaagw55=1,nfterms
xyzzyaagx55=xyzzyaagx55+1
xyzzyaafp55(tagh_forces(xyzzyaagw55,xyzzyaagv55,xyzzyaabf55))=xyzzyaaf&
&k55(xyzzyaagx55)
enddo
enddo
enddo
endif
if(use_future)then
if(tagh_future0>0)xyzzyaafp55(tagh_future0)=xyzzyaagu55(1,0)
if(tagh_future1>0)xyzzyaafp55(tagh_future1)=xyzzyaagu55(1,1)
if(tagh_future2>0)xyzzyaafp55(tagh_future2)=xyzzyaagu55(1,2)
if(tagh_future3>0)xyzzyaafp55(tagh_future3)=xyzzyaagu55(1,3)
if(tagh_future4>0)xyzzyaafp55(tagh_future4)=xyzzyaagu55(1,4)
if(tagh_future5>0)xyzzyaafp55(tagh_future5)=xyzzyaagu55(1,5)
if(tagh_future6>0)xyzzyaafp55(tagh_future6)=xyzzyaagu55(1,6)
if(tagh_future7>0)xyzzyaafp55(tagh_future7)=xyzzyaagu55(1,7)
if(tagh_future8>0)xyzzyaafp55(tagh_future8)=xyzzyaagu55(1,8)
if(tagh_future9>0)xyzzyaafp55(tagh_future9)=xyzzyaagu55(1,9)
if(tagh_future10>0)xyzzyaafp55(tagh_future10)=xyzzyaagu55(1,10)
endif
if(mc_twist_av)then
if(tagh_hf_ke>0)xyzzyaafp55(tagh_hf_ke)=hf_ke
if(tagh_hf_ex>0)xyzzyaafp55(tagh_hf_ex)=hf_ex
endif
if(writeout_dmc_hist)then
xyzzyaafq55(:,xyzzyaaao55)=xyzzyaafq55(:,xyzzyaaao55)+xyzzyaafp55(:)
if(xyzzyaaap55==ndmcave)then
xyzzyaafq55(tagh_step,xyzzyaaao55)=xyzzyaaai55*real(xyzzyaabx1+1,dp)
endif
endif
if(iaccum)call reblock_add(xyzzyaafp55*xyzzyaadn55*xyzzyaaee55,xyzzyaa&
&dn55*xyzzyaaee55)
xyzzyaaej55=.true.
call timer('DMC_AVERAGING',.false.)
endif
if(ibran)then
call timer('DMC_IBRAN',.true.)
call timer('BCAST_EBEST_EREF',.true.)
if(nc_dmc.or.poprenorm)then
call qmpi_bcast(xyzzyaafa55,'DMC_MAIN','ebest')
xyzzyaaet55=xyzzyaafa55
else
if(am_master)then
xyzzyaaet55=xyzzyaafa55-(xyzzyaacw55/xyzzyaael55)*log(xyzzyaadn55*xyzz&
&yaaaj55)
xyzzyaacd55(1)=xyzzyaafa55
xyzzyaacd55(2)=xyzzyaaet55
endif
call qmpi_bcast(xyzzyaacd55,'DMC_MAIN','ebr')
xyzzyaafa55=xyzzyaacd55(1)
xyzzyaaet55=xyzzyaacd55(2)
endif
call timer('BCAST_EBEST_EREF',.false.)
call timer('DMC_IBRAN',.false.)
else
call timer('BCAST_EBEST_EREF',.true.)
call qmpi_bcast(xyzzyaafa55,'DMC_MAIN','ebest')
call timer('BCAST_EBEST_EREF',.false.)
endif
call timer('DMC_IBRAN',.true.)
if(ibran)then
allocate(xyzzyaadm55(xyzzyaaax55),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'DMC','mult_array')
call xyzzyaahm55(xyzzyaadm55)
if(xyzzyaaao55==xyzzyaaaq55.and.xyzzyaaap55==ndmcave)then
if(iblock>=xyzzyaaam55.and.xyzzyaaar55==xyzzyaaaq55)then
call xyzzyaadd1('SERC',xyzzyaaax55,xyzzyaadm55,xyzzyaadk55)
else
select case(chkpoint_level)
case(-1)
call xyzzyaadd1('SEND',xyzzyaaax55,xyzzyaadm55,xyzzyaadk55)
case(0)
if(((isvmc_dmc.or.isdmc_dmc).and.iaccum).or.isdmc_stats.or.isdmc_equil&
&.or.isdmc_old)then
call xyzzyaadd1('SERC',xyzzyaaax55,xyzzyaadm55,xyzzyaadk55)
else
call xyzzyaadd1('SEND',xyzzyaaax55,xyzzyaadm55,xyzzyaadk55)
endif
case(1,2)
call xyzzyaadd1('SERC',xyzzyaaax55,xyzzyaadm55,xyzzyaadk55)
end select
endif
else
call xyzzyaadd1('SEND',xyzzyaaax55,xyzzyaadm55,xyzzyaadk55)
endif
xyzzyaadl55=xyzzyaadl55+xyzzyaadk55
deallocate(xyzzyaadm55)
call timer('BCAST_TOTWEIGHT',.true.)
call qmpi_bcast(xyzzyaadn55,'DMC_MAIN','totweight')
call timer('BCAST_TOTWEIGHT',.false.)
if(trip_popn==0.d0.and.xyzzyaadn55>xyzzyaadq55)then
if(am_master)then
call wout()
tmpr=r2s(xyzzyaadr55,'(f12.3)')
call wordwrap('Population explosion encountered : the hard limit on th&
&e iteration weight ('//trim(tmpr)//' * target weight) has been exceed&
&ed.')
if(.not.iaccum)call wordwrap('Explosions during the early equilibratio&
&n phase may be avoidable by using DMC_EQUIL_FIXPOP.')
call wout()
call wordwrap('You might find it useful to use automatic catastrophe r&
&ecovery, in which case you should set the DMC_TRIP_WEIGHT keyword to &
&e.g. twice the value of DMC_TARGET_WEIGHT.')
call wout()
endif
call errstop_master('DMC','Quitting')
endif
if(trip_popn>0.d0.and.xyzzyaadn55>trip_popn)then
call errwarn('DMC','probable population explosion.')
if(am_master)then
tmpr=r2s(xyzzyaadn55,'(f12.3)')
call wout('Current weight        : '//trim(tmpr))
tmpr=r2s(trip_popn,'(f12.3)')
call wout('Catastrophe threshold : '//trim(tmpr))
endif
xyzzyaaeo55=.true.
endif
call timer('DMC_IBRAN',.false.)
endif
if(nconfig_prelim>0.and.xyzzyaacj55)call xyzzyaahi55
call timer('DMC_STEP',.false.)
if(xyzzyaaeo55)then
call xyzzyaahk55
call timer('DMC_BLOCK',.false.)
cycle block_loop
endif
enddo dmcave_loop
if(xyzzyaack55)exit move_loop
enddo move_loop
if(use_blocktime)then
xyzzyaaar55=xyzzyaaar55-xyzzyaaaq55
xyzzyaaah55=real(xyzzyaaaq55,dp)
xyzzyaaal55=1.d0/xyzzyaaah55
if(xyzzyaaar55>0)then
xyzzyaaam55=xyzzyaaam55+1
endif
endif
call timer('DMC_BLOCK_STUFF',.true.)
if(.not.xyzzyaaep55)xyzzyaaen55=0
if(am_master.and.writeout_dmc_hist)then
if(ndmcave>1)call dscal(no_cols_qmc*xyzzyaaaq55,xyzzyaaak55,xyzzyaafq5&
&5(1,1),1)
if(.not.particle_is_fixed)call write_buffered_data('DMC',xyzzyaaaq55,x&
&yzzyaafq55)
endif
if(.not.iaccum)then
xyzzyaaau55=xyzzyaaau55+xyzzyaaaq55*ndmcave
else
if(real(xyzzyaaav55,dp)+xyzzyaaah55*xyzzyaaai55>max_rep_int)call errst&
&op('DMC_MAIN','Integer overflow: stats_accum_moves.')
xyzzyaaav55=xyzzyaaav55+xyzzyaaaq55*ndmcave
endif
if(use_backflow)then
call backflow_stats(xyzzyaabx55,xyzzyaaby55)
xyzzyaabz55=xyzzyaaby55
if(dmc_cfg_by_cfg)xyzzyaabz55=xyzzyaabx55
call qmpi_reduce(xyzzyaabz55,xyzzyaaca55,mpi_sum,'DMC_MAIN','backflow &
&stats')
xyzzyaaca55=xyzzyaaca55*xyzzyaaag55
endif
if(.not.isperiodic)then
call qmpi_reduce(xyzzyaacb55,xyzzyaacc55,mpi_max,'DMC_MAIN','r2max sta&
&ts')
else
xyzzyaacc55=0.d0
endif
if(am_master)then
xyzzyaabq55=0.d0
if(xyzzyaabh55>0)xyzzyaabq55=100.d0*real(xyzzyaabj55,dp) /real(xyzzyaa&
&bh55,dp)
call wout(repeat('=',73))
call wout('In block : '//trim(i2s(iblock)))
call wout()
call wout('Number of moves in block                 : '//trim(i2s(xyzz&
&yaaaq55)))
if(nnodes>1)then
tmpr=r2s(xyzzyaadp55*100.d0*xyzzyaaal55*xyzzyaaak55,'(f12.3)')
call wout('Load-balancing efficiency (%)            : '//trim(tmpr))
call wout('Number of config transfers               : '//trim(i2s(xyzz&
&yaadl55)))
endif
tmpr=r2s(xyzzyaabq55,'(f12.3)')
call wout('Acceptance ratio (%)                     : '//trim(tmpr))
tmpr=r2s(xyzzyaafa55*xyzzyaafe55,'(f16.8)')
call wout('New best estimate of DMC energy (au)     : '//trim(tmpr))
call wout('Max no of attempts before accept move    : '//trim(i2s(xyzz&
&yaabl55)))
if(iaccum)then
tmpr=r2s(xyzzyaael55,'(f16.8)')
call wout('New best estimate of effective time step : '//trim(tmpr))
endif
if(use_backflow)then
tmpr=r2s(xyzzyaaca55,'(f8.4)')
if(.not.dmc_cfg_by_cfg)then
call wout('Particles affected per move         (%)  : '//trim(tmpr))
else
call wout('Particles within backflow range     (%)  : '//trim(tmpr))
endif
endif
if(.not.isperiodic)then
tmpr=r2s(sqrt(xyzzyaacc55),'(f16.8)')
if(homogeneous_system)then
call wout('Maximum distance from CoM (au)           : '//trim(tmpr))
else
call wout('Maximum distance from origin (au)        : '//trim(tmpr))
endif
endif
call wout()
endif
if(iblock>1.and.trip_popn>0.d0.and..not.xyzzyaaep55)then
call backup_config_file
if(expvals.and.iaccum)call backup_expval_file
endif
xyzzyaaep55=.false.
if(nconfig_prelim==0)then
call xyzzyaahg55(.not.iaccum.and.iblock==xyzzyaaam55,.false.)
else
if(.not.iaccum)call xyzzyaahg55(.false.,.false.)
endif
if(expvals.and.iaccum)call write_expval
call xyzzyaaho55(final=.false.)
if(am_master)then
call wout('Time taken in block    : : : ',dble(tcputime()-xyzzyaagh55)&
&,rfmt='(f13.4)')
call wout()
endif
call timer('DMC_BLOCK_STUFF',.false.)
call timer('DMC_BLOCK',.false.)
xyzzyaagg55=.true.
if(iblock>=xyzzyaaam55)then
if(iaccum.and.nconfig_prelim>0)then
call mpi_reduce(prelim_nconfig_written==nconfig_prelim,xyzzyaafs55,1,m&
&pi_logical,mpi_land,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Checking whether all configurations have been wr&
&itten.')
call mpi_bcast(xyzzyaafs55,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting whether all configurations have bee&
&n written.')
if(xyzzyaafs55)exit block_loop
else
exit block_loop
endif
endif
iblock=iblock+1
enddo block_loop
call timer('DMC_POSTBLOCK_STUFF',.true.)
call xyzzyaaho55(final=.true.)
if(iaccum.and.nconfig_prelim>0)call xyzzyaahj55
call xyzzyaacm1
deallocate(xyzzyaacf55,xyzzyaabo55)
if(tpdmc>0)deallocate(xyzzyaaeh55)
deallocate(xyzzyaaei55)
if(allocated(xyzzyaafq55))deallocate(xyzzyaafq55)
deallocate(xyzzyaafp55)
deallocate(xyzzyaafl55,xyzzyaafm55,xyzzyaafk55,xyzzyaafn55)
if(dmc_cfg_by_cfg)deallocate(xyzzyaacg55,xyzzyaadf55,xyzzyaade55,xyzzy&
&aadx55)
if(use_future)deallocate(xyzzyaags55,xyzzyaagt55,xyzzyaagu55)
if(am_master.and..not.iaccum)deallocate(xyzzyaafo55)
if(forces)call finish_forces_accum
if(expvals)call finish_expval_accum
if(am_master.and.iaccum)call finish_reblock
call finish_energy_utils
call finish_wfn_utils
call finish_scratch
call timer('DMC_POSTBLOCK_STUFF',.false.)
call timer('DMC',.false.)
contains
subroutine xyzzyaaha55
implicit none
if(targ_wt<=0.d0)call errstop('INITIALIZE_DMC','NCONFIG in input appea&
&rs to be zero.')
allocate(xyzzyaacf55(3,netot),xyzzyaabo55(netot),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'INITIALIZE_DMC','basic')
xyzzyaacf55=0.d0
xyzzyaabo55=0
xyzzyaafc55=0.d0
xyzzyaaed55=0.d0
xyzzyaaee55=1.d0
xyzzyaaef55=0.d0
xyzzyaaeg55=0.d0
xyzzyaada55(:)=0.5d0*inv_pmass(:)
xyzzyaadb55(:)=dt*inv_pmass(:)
xyzzyaadd55(:)=sqrt(xyzzyaadb55(:))
xyzzyaadq55=xyzzyaadr55*targ_wt
if(model_system.and.isperiodic)then
xyzzyaafe55=inv_netot
else
xyzzyaafe55=1.d0/real(npcells,dp)
endif
xyzzyaaam55=nblock_in
xyzzyaaaq55=nmove_in
xyzzyaaah55=real(xyzzyaaaq55,dp)
xyzzyaaai55=real(ndmcave,dp)
xyzzyaaal55=1.d0/xyzzyaaah55
xyzzyaaag55=1.d0/real(nnodes,dp)
xyzzyaaaj55=1.d0/targ_wt
xyzzyaaak55=1.d0/xyzzyaaai55
xyzzyaace55=1.d0/sqrt(dt)
if(poprenorm.or.nc_dmc)then
if(lwdmc)call errstop('INITIALIZE_DMC','Cannot perform weighted DMC in&
& conjunction with population renormalization.')
if(abs(targ_wt-anint(targ_wt))>1.d-7)call errstop('INITIALIZE_DMC','Ta&
&rget total weight should be an integer if population renormalization &
&is to be carried out.')
if(tpdmc>0.or.growth_estimator)call errstop('INITIALIZE_DMC','Cannot u&
&se Pi-weighting scheme or growth estimator with population renormaliz&
&ation or norm-conserving DMC.')
endif
nullify(xyzzyaagd55,xyzzyaage55)
xyzzyaacw55=min(1.d0,cerefdmc*dt)
if(nucleus_gf_mods)call xyzzyaadq1(dt)
allocate(xyzzyaaeh55(0:tpdmc-1),xyzzyaaei55(0:tpdmc),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'INITIALIZE_DMC','Pi-weight')
allocate(xyzzyaafp55(no_cols_qmc),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'INITIALIZE_DMC','hist_data')
if(am_master.and.writeout_dmc_hist)then
allocate(xyzzyaafq55(no_cols_qmc,xyzzyaaaq55),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'INITIALIZE_DMC','hist_buff')
endif
xyzzyaaaf55=real(ebest_av_window,dp)
if(am_master.and..not.iaccum)then
allocate(xyzzyaafo55(ebest_av_window),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'INITIALIZE_DMC','equil_etot')
endif
xyzzyaaep55=.false.
if(dmc_cfg_by_cfg)then
allocate(xyzzyaacg55(3,netot),xyzzyaadf55(3,netot),xyzzyaade55(netot),&
&xyzzyaadx55(netot),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'INITIALIZE_DMC','cfg-by-cfg')
endif
end subroutine xyzzyaaha55
subroutine xyzzyaahb55
implicit none
xyzzyaafv55=n_ecomp
xyzzyaaft55=xyzzyaafv55+1
xyzzyaafu55=xyzzyaafv55+2
xyzzyaafv55=xyzzyaafv55+2
xyzzyaagk55=0
xyzzyaagl55=0
xyzzyaagj55=0
if(use_future)then
xyzzyaagk55=1
xyzzyaagl55=0
xyzzyaagj55=nint(xyzzyaagp55/dt)
endif
if(forces)then
xyzzyaafw55=xyzzyaafv55+1
xyzzyaafx55=xyzzyaafv55+nfcomps
xyzzyaafv55=xyzzyaafx55
xyzzyaagk55=11*nitot_forces*naxis_forces
endif
allocate(xyzzyaafl55(xyzzyaafv55),xyzzyaafm55(xyzzyaafv55),xyzzyaafn55&
&(xyzzyaafv55),xyzzyaafk55(xyzzyaafv55),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'SETUP_ACCUMULATION_ARRAYS','sum_data arra&
&ys')
if(use_future)then
allocate(xyzzyaags55(xyzzyaagk55,0:xyzzyaago55),xyzzyaagt55(xyzzyaagk5&
&5,0:xyzzyaago55),xyzzyaagu55(xyzzyaagk55,0:xyzzyaago55),xyzzyaagq55(x&
&yzzyaagk55),xyzzyaagr55(xyzzyaagk55),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'SETUP_ACCUMULATION_ARRAYS','future_walkin&
&g')
xyzzyaags55=0.d0
xyzzyaagt55=0.d0
xyzzyaagu55=0.d0
endif
end subroutine xyzzyaahb55
subroutine xyzzyaahc55(from_backup)
implicit none
logical,intent(in) :: from_backup
integer xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58
real(dp) xyzzyaaae58,xyzzyaaaf58,xyzzyaaag58,xyzzyaaah58,xyzzyaaai58,x&
&yzzyaaaj58,etot,xyzzyaaak58,xyzzyaaal58,xyzzyaaam58,xyzzyaaan58,xyzzy&
&aaao58,xyzzyaaap58,xyzzyaaaq58,xyzzyaaar58
complex(dp) xyzzyaaas58,xyzzyaaat58
logical xyzzyaaau58(3),xyzzyaaav58(4),xyzzyaaaw58,xyzzyaaax58,xyzzyaaa&
&y58
character(20) config_item(4),opt_citem(3),extra_item(4),extra_dum(0)
call timer('READCONFIGS_DMC',.true.)
xyzzyaaay58=.false.
xyzzyaaaa58=3
config_item(1:3)=(/'RELE  ','LOGDET','ETOT  '/)
if(noncoll_spin)then
xyzzyaaaa58=xyzzyaaaa58+1
config_item(xyzzyaaaa58)='SELE  '
endif
opt_citem=(/'STOT  ','WDMC  ','VALJAS'/)
extra_item(1:4)=(/'DMC_SAVED_STATE','RANDOM         ','REBLOCK_DATA   &
&','GEOMETRY       '/)
xyzzyaaax55=-2
call load_configs(xyzzyaaax55,'VMC_DMC',config_item(1:xyzzyaaaa58),opt&
&_citem,xyzzyaaau58,extra_dum,extra_item,xyzzyaaav58,from_backup=from_&
&backup)
if(.not.xyzzyaaav58(1))then
if(dble(xyzzyaaax55)<0.5d0*targ_wt/dble(nnodes))call errwarn_silent('R&
&EADCONFIGS_DMC','Too few VMC-generated configs for current target wei&
&ght - the equilibration stage will take longer than usual.')
if(dble(xyzzyaaax55)>2.d0*targ_wt/dble(nnodes))then
call errwarn_silent('READCONFIGS_DMC','Too many VMC-generated configs &
&for current target weight - not all configs will be loaded to avoid s&
&purious population explosion reports. The equilibration stage will ta&
&ke longer than usual.')
xyzzyaaax55=int(2.d0*targ_wt/dble(nnodes))
if(xyzzyaaax55<1)xyzzyaaax55=1
endif
xyzzyaaau55=0
xyzzyaaav55=0
stot_config=etot_config
wdmc_config=1.d0
xyzzyaaag58=0.d0
if(xyzzyaaax55>0)xyzzyaaag58=sum(etot_config(1:xyzzyaaax55))
call qmpi_reduce(xyzzyaaax55,xyzzyaaac58,mpi_sum,'READCONFIGS_DMC','nc&
&onfig')
call qmpi_reduce(xyzzyaaag58,xyzzyaaaj58,mpi_sum,'READCONFIGS_DMC','el&
&oc')
if(am_master)xyzzyaafa55=xyzzyaaaj58/real(xyzzyaaac58,dp)
call qmpi_bcast(xyzzyaafa55,'READCONFIGS_DMC','ebest')
xyzzyaafb55=xyzzyaafa55
if(use_init_eref)then
xyzzyaaet55=dmc_init_eref
xyzzyaafa55=dmc_init_eref
else
xyzzyaaet55=xyzzyaafa55
endif
xyzzyaaem55=0.d0
xyzzyaael55=0.d0
xyzzyaaec55=0.d0
xyzzyaaea55=0.d0
xyzzyaaeb55=0.d0
xyzzyaady55=0.d0
xyzzyaafg55(1:13)=0.d0
if(tpdmc>0)xyzzyaaeh55=0.d0
if(growth_estimator)then
xyzzyaaei55=0.d0
xyzzyaadz55=0.d0
endif
xyzzyaaej55=.false.
if(.not.newrun)call errstop_master('READCONFIGS_DMC','VMC-generated co&
&nfigs, but NEWRUN=F.')
else
xyzzyaafa55=ebest_config
xyzzyaaet55=eref_config
xyzzyaael55=dteff_best_config
if(.not.lwdmc.or..not.lwdmc_config)wdmc_config=1.d0
if(from_backup.or..not.newrun)then
xyzzyaaau55=dmcequil_steps_config
xyzzyaaav55=dmcstats_steps_config
else
xyzzyaaau55=0
xyzzyaaav55=0
endif
if(tpdmc/=tpdmc_config)then
xyzzyaaau55=0
xyzzyaaav55=0
elseif(.not.iaccum.and.xyzzyaaav55>0)then
if(real(xyzzyaaau55,dp)+real(xyzzyaaav55,dp)>max_rep_int)call errstop_&
&master('READCONFIGS_DMC','Integer overflow: equilibration_moves.')
xyzzyaaau55=xyzzyaaau55+xyzzyaaav55
xyzzyaaav55=0
endif
if(xyzzyaaav55>0)then
xyzzyaafb55=ebest_init_config
if(growth_estimator.and..not.growth_estimator_config)then
xyzzyaaec55=0.d0
xyzzyaaeb55=0.d0
xyzzyaadz55=0.d0
xyzzyaaei55=0.d0
elseif(growth_estimator)then
xyzzyaaec55=numerator_wt2_config
xyzzyaaeb55=denominator_wt2_config
xyzzyaaei55=log_pi_wt_array2_config
xyzzyaadz55=log_pi_wt2_config
endif
xyzzyaaea55=denominator_wt_config
xyzzyaady55=log_pi_wt_config
xyzzyaafg55=numer_expect_config
if(tpdmc>0)xyzzyaaeh55=log_pi_wt_array_config
xyzzyaaej55=.true.
xyzzyaaem55=dteff_ebest_init_config
if(am_master)then
if(xyzzyaaav58(1).and.xyzzyaaav58(3))call restore_reblock_data
endif
else
xyzzyaafb55=xyzzyaafa55
xyzzyaafg55=0.d0
xyzzyaadz55=0.d0
xyzzyaaea55=0.d0
xyzzyaaec55=0.d0
xyzzyaaeb55=0.d0
xyzzyaady55=0.d0
xyzzyaaei55=0.d0
if(tpdmc>0)xyzzyaaeh55=0.d0
xyzzyaaej55=.false.
endif
endif
if(xyzzyaaav58(2))call put_random_state(random_state_config)
if(xyzzyaaax55>0)then
do xyzzyaaab58=1,xyzzyaaax55
nullify(xyzzyaagd55)
call xyzzyaabz1(xyzzyaagd55)
call xyzzyaacn1(xyzzyaagd55,xyzzyaaab58)
enddo
nullify(xyzzyaagd55)
endif
call dismantle_configs
xyzzyaaax58=.false.
if(xyzzyaaav58(1).and.(.not.from_backup))then
if(dmc_twist_av)then
xyzzyaaax58=.not.iaccum
elseif(dmc_reweight_configs.or.dmc_spacewarping)then
xyzzyaaax58=newrun
endif
endif
if(xyzzyaaax58)then
call timer('RECOMPUTE_CFGS',.true.)
if(dmc_spacewarping)then
if(.not.xyzzyaaav58(4))call errstop('READCONFIGS_DMC','DMC_SPACEWARPIN&
&G requires old GEOMETRY data. None found in config.in.')
else
if(.not.xyzzyaaau58(2))call errstop('READCONFIGS_DMC','Config reweight&
&ing requires configuration weights. None found in config.in.')
if(.not.xyzzyaaau58(3).and.use_jastrow)call errstop('READCONFIGS_DMC',&
&'Config reweighting requires VALJAS values. None found in config.in (&
&JASBUF has to be set)')
if(am_master.and.dmc_reweight_configs)then
call wout('DMC_REWEIGHT_CONFIGS activated in input.')
call wout('Updating weights of walkers read from config.in.')
call wout()
endif
endif
xyzzyaaam58=0.d0
xyzzyaaal58=0.d0
xyzzyaaah58=0.d0
xyzzyaaag58=0.d0
xyzzyaaaq58=0.d0
xyzzyaaar58=0.d0
xyzzyaaao58=1.d10
xyzzyaaap58=0.d0
xyzzyaaad58=0
if(xyzzyaaax55>0)then
xyzzyaagd55=>xyzzyaaar1
do
call xyzzyaacb1(xyzzyaaeq55,xyzzyaagd55,xyzzyaacf55,xyzzyaabo55,xyzzya&
&aaf58,xyzzyaaai58,wdmc,xyzzyaabn55)
xyzzyaaah58=xyzzyaaah58+xyzzyaaaf58*wdmc
xyzzyaaam58=xyzzyaaam58+wdmc
if(dmc_spacewarping)then
call xyzzyaahf55(rion_config(:,:),rion(:,:),xyzzyaacf55(:,:),wdmc)
call define_config(xyzzyaaeq55,xyzzyaacf55,xyzzyaabo55)
call eval_local_energy(xyzzyaaeq55,etot=xyzzyaaae58)
else
xyzzyaaaq58=xyzzyaaaq58+wdmc*wdmc
call wfn_logval(xyzzyaaeq55,xyzzyaaas58,xyzzyaaaw58)
if(xyzzyaaaw58)call errstop('READCONFIGS_DMC','Walker found on wfn=0 p&
&oint.')
call clear_scratch_wfn(xyzzyaaeq55,wfn_only=.true.)
call eval_local_energy(xyzzyaaeq55,etot=xyzzyaaae58)
call wfn_logval(xyzzyaaeq55,xyzzyaaat58,xyzzyaaaw58)
if(xyzzyaaaw58)then
wdmc=0.d0
xyzzyaaao58=0.d0
xyzzyaaad58=xyzzyaaad58+1
else
xyzzyaaan58=exp(2.d0*real(xyzzyaaat58,kind=dp)-2.d0*real(xyzzyaaas58,k&
&ind=dp))
xyzzyaaan58=min(xyzzyaaan58,10.d0)
wdmc=wdmc*xyzzyaaan58
xyzzyaaao58=min(xyzzyaaan58,xyzzyaaao58)
xyzzyaaap58=max(xyzzyaaan58,xyzzyaaap58)
xyzzyaaar58=xyzzyaaar58+wdmc*wdmc
endif
endif
xyzzyaaag58=xyzzyaaag58+wdmc*xyzzyaaae58
xyzzyaaal58=xyzzyaaal58+wdmc
call xyzzyaace1(xyzzyaaeq55,xyzzyaagd55,xyzzyaaae58,xyzzyaaai58,wdmc,x&
&yzzyaabn55)
xyzzyaagd55=>xyzzyaagd55%pt_next_config
if(associated(xyzzyaagd55,xyzzyaaar1))exit
enddo
nullify(xyzzyaagd55)
endif
if(dmc_reweight_configs.or.dmc_twist_av)then
call qmpi_reduce(xyzzyaaad58,xyzzyaaab58,mpi_sum,'READCONFIGS_DMC','wz&
&ero_counter')
call qmpi_reduce(xyzzyaaax55,xyzzyaaac58,mpi_sum,'READCONFIGS_DMC','nc&
&onfig')
if(am_master)then
if(xyzzyaaab58>0)call wout('Info: '//trim(i2s(xyzzyaaab58))//' of '//t&
&rim(i2s(xyzzyaaac58))//' walkers reweighted to zero weight.')
endif
endif
call qmpi_reduce(xyzzyaaal58,xyzzyaaaj58,mpi_sum,'READCONFIGS_DMC','ws&
&um')
if(am_master)xyzzyaaal58=xyzzyaaaj58
call qmpi_bcast(xyzzyaaal58,'READCONFIGS_DMC','wsum')
if(xyzzyaaal58<=0)call errstop_master('READCONFIGS_DMC','Configuration&
&s all dead.')
call qmpi_reduce(xyzzyaaag58,xyzzyaaaj58,mpi_sum,'READCONFIGS_DMC','el&
&oc_sum')
if(am_master)xyzzyaaag58=xyzzyaaaj58
call qmpi_reduce(xyzzyaaam58,xyzzyaaaj58,mpi_sum,'READCONFIGS_DMC','ws&
&um_old')
if(am_master)xyzzyaaam58=xyzzyaaaj58
call qmpi_bcast(xyzzyaaam58,'READCONFIGS_DMC','wsum_old')
call qmpi_reduce(xyzzyaaah58,xyzzyaaaj58,mpi_sum,'READCONFIGS_DMC','el&
&oc_old_sum')
if(am_master)xyzzyaaah58=xyzzyaaaj58
if(dmc_reweight_configs.or.dmc_twist_av)then
if(xyzzyaaax55>0)then
xyzzyaagd55=>xyzzyaaar1
do
xyzzyaagd55%pt_dmc%wdmc=xyzzyaagd55%pt_dmc%wdmc*(xyzzyaaam58/xyzzyaaal&
&58)
xyzzyaagd55=>xyzzyaagd55%pt_next_config
if(associated(xyzzyaagd55,xyzzyaaar1))exit
enddo
nullify(xyzzyaagd55)
endif
xyzzyaadn55=xyzzyaaam58
else
xyzzyaadn55=xyzzyaaal58
endif
xyzzyaado55=1.d0/xyzzyaadn55
if(dmc_reweight_configs)then
call qmpi_reduce(xyzzyaaaq58,xyzzyaaaj58,mpi_sum,'READCONFIGS_DMC','ws&
&qsum_old')
if(am_master)xyzzyaaaq58=xyzzyaaaj58
call qmpi_reduce(xyzzyaaar58,xyzzyaaaj58,mpi_sum,'READCONFIGS_DMC','ws&
&qsum_new')
if(am_master)xyzzyaaar58=xyzzyaaaj58
call qmpi_reduce(xyzzyaaao58,xyzzyaaaj58,mpi_min,'READCONFIGS_DMC','wf&
&ac_min')
if(am_master)xyzzyaaao58=xyzzyaaaj58
call qmpi_reduce(xyzzyaaap58,xyzzyaaaj58,mpi_max,'READCONFIGS_DMC','wf&
&ac_max')
if(am_master)then
xyzzyaaap58=xyzzyaaaj58
call wout('Old total weight              : ',xyzzyaaam58)
call wout('Old effective population size : ',xyzzyaaam58**2/xyzzyaaaq5&
&8)
call wout('Maximum weight factor         : ',xyzzyaaap58)
call wout('Minimum weight factor         : ',xyzzyaaao58)
call wout('New total weight (pre renorm) : ',xyzzyaaal58)
call wout('New effective population size : ',xyzzyaaal58**2/xyzzyaaar5&
&8)
tmpr=r2s((xyzzyaaal58**2*xyzzyaaaq58)/(xyzzyaaam58**2*xyzzyaaar58)*100&
&.d0,'(f24.2)')
call wout('Reweighting efficiency        : '//trim(tmpr)//'%')
call wout()
endif
endif
if(am_master)then
xyzzyaaak58=xyzzyaaah58/xyzzyaaam58
etot=xyzzyaaag58/xyzzyaaal58
if(dmc_spacewarping)then
call wout('DMC space warping activated.')
call wout()
endif
call wout('ETOT (old)                    : ',xyzzyaaak58*xyzzyaafe55)
call wout('ETOT (new)                    : ',etot*xyzzyaafe55)
call wout('EBEST (old)                   : ',xyzzyaafa55*xyzzyaafe55)
xyzzyaafa55=xyzzyaafa55-xyzzyaaak58+etot
call wout('EBEST (new)                   : ',xyzzyaafa55*xyzzyaafe55)
call wout('EREF (old)                    : ',xyzzyaaet55*xyzzyaafe55)
xyzzyaaet55=xyzzyaaet55-xyzzyaaak58+etot
call wout('EREF (new)                    : ',xyzzyaaet55*xyzzyaafe55)
call wout('Energies in au/prim_cell including N-N.')
call wout()
xyzzyaaay58=.true.
endif
call qmpi_bcast(xyzzyaafa55,'READCONFIGS_DMC','ebest')
call qmpi_bcast(xyzzyaaet55,'READCONFIGS_DMC','eref')
xyzzyaafb55=xyzzyaafa55
if(dmc_reweight_configs.or.dmc_twist_av)then
allocate(xyzzyaadm55(xyzzyaaax55),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'DMC','mult_array')
call xyzzyaahm55(xyzzyaadm55)
call xyzzyaadd1('SERC',xyzzyaaax55,xyzzyaadm55,xyzzyaadk55)
xyzzyaadl55=xyzzyaadl55+xyzzyaadk55
deallocate(xyzzyaadm55)
endif
call timer('RECOMPUTE_CFGS',.false.)
endif
if(am_master)then
if(.not.xyzzyaaay58)then
call wout('EBEST = ',xyzzyaafa55*xyzzyaafe55,' (au/prim cell inc. N-N)&
&')
call wout('EREF  = ',xyzzyaaet55*xyzzyaafe55)
call wout()
endif
if(.not.newrun.and..not.dmc_twist_av)then
call wout('Number of previous DMC stats accumulation moves : '//trim(i&
&2s(xyzzyaaav55)))
call wout()
endif
endif
call timer('READCONFIGS_DMC',.false.)
end subroutine xyzzyaahc55
subroutine xyzzyaahd55
implicit none
call xyzzyaacb1(xyzzyaaeq55,xyzzyaagd55,xyzzyaacf55,xyzzyaabo55,xyzzya&
&aeu55,xyzzyaaex55,wdmc,xyzzyaabn55)
if(particle_is_fixed)call xyzzyaahl55
xyzzyaabg55=0
xyzzyaabi55=0
call xyzzyaadl1
if(.not.dmc_cfg_by_cfg)then
call timer('DRIFT_DIFFUSE',.true.)
electron_loop: do xyzzyaabd55=1,netot
if(xyzzyaagy55)then
call wout(' PARTICLE #'//trim(i2s(xyzzyaabd55))//':')
call wout('  R_OLD: ',xyzzyaacf55(1:dimensionality,xyzzyaabd55))
endif
xyzzyaabe55=which_spin(xyzzyaabd55)
xyzzyaabg55=xyzzyaabg55+1
xyzzyaacy55=xyzzyaada55(xyzzyaabe55)
xyzzyaacx55=xyzzyaadb55(xyzzyaabe55)
xyzzyaadc55=xyzzyaadd55(xyzzyaabe55)
call wfn_loggrad(xyzzyaabd55,xyzzyaaeq55,0,xyzzyaadh55,prefetch_val=.t&
&rue.,isnan=isnan,isinf=isinf)
xyzzyaacm55=real(xyzzyaadh55,dp)
if(xyzzyaagy55)then
if(isnan)then
call wout('  DRIFT_OLD: Not a Number')
elseif(isinf)then
call wout('  DRIFT_OLD: Diverges')
else
call wout('  DRIFT_OLD: ',xyzzyaacm55(1:dimensionality))
endif
endif
if(isnan.or.isinf)call errstop('MOVE_CONFIG','Floating-point exception&
& reported when calculating old drift vector. Don''t know what to do, &
&so stopping. Please file a bug report.')
if(nucleus_gf_mods)then
call get_eivecs(xyzzyaaeq55)
call xyzzyaadr1(xyzzyaacm55,eivecs_scr(1,1,xyzzyaabd55,xyzzyaaeq55),al&
&imit,xyzzyaads55,xyzzyaadt55,xyzzyaadu55,xyzzyaadw55)
endif
call xyzzyaado1(alimit,dt,xyzzyaacm55,xyzzyaacl55,xyzzyaacy55)
if(nucleus_gf_mods)then
call xyzzyaadt1(dt,xyzzyaacl55,xyzzyaads55,xyzzyaadt55,xyzzyaadu55,xyz&
&zyaadw55,xyzzyaacf55(1,xyzzyaabd55),xyzzyaabp55,xyzzyaacp55,xyzzyaadv&
&55)
else
call xyzzyaads1(xyzzyaacf55(1,xyzzyaabd55),xyzzyaabp55,xyzzyaacl55,xyz&
&zyaacx55,xyzzyaadc55,xyzzyaacp55)
endif
if(xyzzyaagy55)call wout('  R_NEW: ',xyzzyaabp55(1:dimensionality))
xyzzyaacq55=ddot(3,xyzzyaacp55(1),1,xyzzyaacp55(1),1)
call define_config_oneelec(xyzzyaabd55,xyzzyaaeq55,xyzzyaaer55,xyzzyaa&
&bp55,xyzzyaabo55(xyzzyaabd55))
call wfn_ratio(xyzzyaaeq55,xyzzyaaer55,0,ratio=xyzzyaadg55,relprob=rel&
&prob,prefetch_fd=.true.,isnan=isnan,isinf=isinf)
if(xyzzyaagy55)then
if(isnan)then
call wout('  WFN_RATIO: Not a Number')
elseif(isinf)then
call wout('  WFN_RATIO: Diverges')
else
call wout('  WFN_RATIO: ',xyzzyaadg55)
endif
endif
if(isnan.or.isinf)then
xyzzyaadg55=czero
relprob=0.d0
endif
xyzzyaaci55=(complex_wf.or.real(xyzzyaadg55,dp)>0.d0.or..not.ibran).an&
&d.relprob>1.d-200
if(particle_is_fixed)then
if(xyzzyaabd55==fixed_particle)xyzzyaaci55=.false.
endif
if(xyzzyaaci55)then
call wfn_loggrad(xyzzyaabd55,xyzzyaaer55,0,xyzzyaadh55,isnan=isnan,isi&
&nf=isinf)
xyzzyaaco55=real(xyzzyaadh55,dp)
if(xyzzyaagy55)then
if(isnan)then
call wout('  DRIFT_NEW: Not a Number')
elseif(isinf)then
call wout('  DRIFT_NEW: Diverges')
else
call wout('  DRIFT_NEW: ',xyzzyaaco55(1:dimensionality))
endif
endif
if(isnan.or.isinf)then
xyzzyaadg55=czero
relprob=0.d0
prob=0.d0
xyzzyaaci55=.false.
else
if(nucleus_gf_mods)then
call get_eivecs(xyzzyaaer55)
call xyzzyaadr1(xyzzyaaco55,eivecs_scr(1,1,xyzzyaabd55,xyzzyaaer55),al&
&imit,xyzzyaads55,xyzzyaadt55,xyzzyaadu55,xyzzyaadw55)
endif
call xyzzyaado1(alimit,dt,xyzzyaaco55,xyzzyaacn55,xyzzyaacy55)
if(nucleus_gf_mods)then
relprob=relprob*xyzzyaady1(dt,xyzzyaacn55,xyzzyaads55,xyzzyaadt55,xyzz&
&yaadu55,xyzzyaadw55,xyzzyaacf55(1,xyzzyaabd55),xyzzyaabp55,xyzzyaadv5&
&5)
else
relprob=relprob*xyzzyaadk1(xyzzyaacf55(1,xyzzyaabd55),xyzzyaabp55,xyzz&
&yaacl55,xyzzyaacn55,dt,xyzzyaacy55)
endif
prob=min(1.d0,relprob)
xyzzyaabr55=0.d0
if(prob<1.d0)xyzzyaabr55=ranx()
xyzzyaaci55=xyzzyaabr55<prob
endif
else
xyzzyaadg55=czero
prob=0.d0
if(xyzzyaagy55)call wout('  NODE CROSSED')
endif
if(xyzzyaaci55)then
if(xyzzyaagy55)call wout('  MOVE ACCEPTED')
call accept_move(xyzzyaaeq55,xyzzyaaer55)
xyzzyaacf55(:,xyzzyaabd55)=xyzzyaabp55(:)
xyzzyaabi55=xyzzyaabi55+1
xyzzyaabn55(xyzzyaabd55)=0
call xyzzyaahp55(xyzzyaabd55,xyzzyaacf55,xyzzyaacb55)
else
if(xyzzyaagy55)call wout('  MOVE REJECTED')
xyzzyaabn55(xyzzyaabd55)=xyzzyaabn55(xyzzyaabd55)+1
endif
call xyzzyaadm1(xyzzyaacq55,prob,xyzzyaaci55,xyzzyaabe55)
enddo electron_loop
call timer('DRIFT_DIFFUSE',.false.)
else
call timer('DRIFT_DIFFUSE',.true.)
if(xyzzyaagy55)then
do xyzzyaabd55=1,netot
call wout(' R_'//trim(i2s(xyzzyaabd55))//'_OLD: ',xyzzyaacf55(1:dimens&
&ionality,xyzzyaabd55))
enddo
endif
call prefetch_wfn(xyzzyaaeq55,.true.,.false.)
if(nucleus_gf_mods)call get_eivecs(xyzzyaaeq55)
electron_loop_cfg_by_cfg: do xyzzyaabd55=1,netot
xyzzyaabg55=xyzzyaabg55+1
xyzzyaabe55=which_spin(xyzzyaabd55)
xyzzyaacy55=xyzzyaada55(xyzzyaabe55)
xyzzyaacx55=xyzzyaadb55(xyzzyaabe55)
xyzzyaadc55=xyzzyaadd55(xyzzyaabe55)
call wfn_loggrad(xyzzyaabd55,xyzzyaaeq55,0,xyzzyaadh55,isnan=isnan,isi&
&nf=isinf)
xyzzyaacm55=real(xyzzyaadh55,dp)
if(xyzzyaagy55)then
if(isnan)then
call wout(' DRIFT_'//trim(i2s(xyzzyaabd55))//'_OLD: Not a Number')
elseif(isinf)then
call wout(' DRIFT_'//trim(i2s(xyzzyaabd55))//'_OLD: Diverges')
else
call wout(' DRIFT_'//trim(i2s(xyzzyaabd55))//'_OLD: ',xyzzyaacm55(1:di&
&mensionality))
endif
endif
if(isnan.or.isinf)call errstop('MOVE_CONFIG','Floating-point exception&
& reported when calculating old drift vector. Don''t know what to do, &
&so stopping. Please file a bug report.')
if(nucleus_gf_mods)then
call xyzzyaadr1(xyzzyaacm55,eivecs_scr(1,1,xyzzyaabd55,xyzzyaaeq55),al&
&imit,xyzzyaads55,xyzzyaadt55,xyzzyaadu55,xyzzyaadw55)
endif
call xyzzyaado1(alimit,dt,xyzzyaacm55,xyzzyaadf55(1,xyzzyaabd55),xyzzy&
&aacy55)
if(nucleus_gf_mods)then
call xyzzyaadt1(dt,xyzzyaadf55(1,xyzzyaabd55),xyzzyaads55,xyzzyaadt55,&
&xyzzyaadu55,xyzzyaadw55,xyzzyaacf55(1,xyzzyaabd55),xyzzyaacg55(1,xyzz&
&yaabd55),xyzzyaacp55,xyzzyaadx55(xyzzyaabd55))
else
call xyzzyaads1(xyzzyaacf55(1,xyzzyaabd55),xyzzyaacg55(1,xyzzyaabd55),&
&xyzzyaadf55(1,xyzzyaabd55),xyzzyaacx55,xyzzyaadc55,xyzzyaacp55)
endif
xyzzyaade55(xyzzyaabd55)=ddot(3,xyzzyaacp55(1),1,xyzzyaacp55(1),1)
enddo electron_loop_cfg_by_cfg
if(xyzzyaagy55)then
do xyzzyaabd55=1,netot
call wout(' R_'//trim(i2s(xyzzyaabd55))//'_NEW: ',xyzzyaacg55(1:dimens&
&ionality,xyzzyaabd55))
enddo
endif
call define_config(xyzzyaaer55,xyzzyaacg55,xyzzyaabo55)
call wfn_ratio(xyzzyaaeq55,xyzzyaaer55,0,ratio=xyzzyaadg55,relprob=rel&
&prob,isnan=isnan,isinf=isinf)
if(xyzzyaagy55)then
if(isnan)then
call wout(' WFN_RATIO: Not a Number')
elseif(isinf)then
call wout(' WFN_RATIO: Diverges')
else
call wout(' WFN_RATIO: ',xyzzyaadg55)
endif
endif
if(isnan.or.isinf)then
xyzzyaadg55=czero
relprob=0.d0
endif
xyzzyaaci55=(complex_wf.or.real(xyzzyaadg55,dp)>0.d0.or..not.ibran).an&
&d.abs(relprob)>1.d-200
if(xyzzyaaci55)then
call prefetch_wfn(xyzzyaaer55,.true.,.false.)
if(nucleus_gf_mods)call get_eivecs(xyzzyaaer55)
xyzzyaacz55=1.d0
do xyzzyaabd55=1,netot
xyzzyaacy55=xyzzyaada55(which_spin(xyzzyaabd55))
call wfn_loggrad(xyzzyaabd55,xyzzyaaer55,0,xyzzyaadh55,isnan=isnan,isi&
&nf=isinf)
xyzzyaaco55=real(xyzzyaadh55,dp)
if(xyzzyaagy55)then
if(isnan)then
call wout(' DRIFT_'//trim(i2s(xyzzyaabd55))//'_NEW: Not a Number')
elseif(isinf)then
call wout(' DRIFT_'//trim(i2s(xyzzyaabd55))//'_NEW: Diverges')
else
call wout(' DRIFT_'//trim(i2s(xyzzyaabd55))//'_NEW: ',xyzzyaaco55(1:di&
&mensionality))
endif
endif
if(isnan.or.isinf)exit
if(nucleus_gf_mods)then
call xyzzyaadr1(xyzzyaaco55,eivecs_scr(1,1,xyzzyaabd55,xyzzyaaer55),al&
&imit,xyzzyaads55,xyzzyaadt55,xyzzyaadu55,xyzzyaadw55)
endif
call xyzzyaado1(alimit,dt,xyzzyaaco55,xyzzyaacn55,xyzzyaacy55)
if(nucleus_gf_mods)then
xyzzyaacz55=xyzzyaacz55*xyzzyaady1(dt,xyzzyaacn55,xyzzyaads55,xyzzyaad&
&t55,xyzzyaadu55,xyzzyaadw55,xyzzyaacf55(1,xyzzyaabd55),xyzzyaacg55(1,&
&xyzzyaabd55),xyzzyaadx55(xyzzyaabd55))
else
xyzzyaacz55=xyzzyaacz55*xyzzyaadk1(xyzzyaacf55(1,xyzzyaabd55),xyzzyaac&
&g55(1,xyzzyaabd55),xyzzyaadf55(1,xyzzyaabd55),xyzzyaacn55,dt,xyzzyaac&
&y55)
endif
enddo
if(isnan.or.isinf)then
relprob=0.d0
prob=0.d0
xyzzyaabr55=0.d0
xyzzyaaci55=.false.
else
relprob=relprob*xyzzyaacz55
prob=min(1.d0,relprob)
xyzzyaabr55=0.d0
if(prob<1.d0)xyzzyaabr55=ranx()
xyzzyaaci55=xyzzyaabr55<prob
endif
else
prob=0.d0
xyzzyaaci55=.false.
if(xyzzyaagy55)call wout(' NODE CROSSED')
endif
if(xyzzyaaci55)then
call accept_move(xyzzyaaeq55,xyzzyaaer55)
call dcopy(three_netot,xyzzyaacg55(1,1),1,xyzzyaacf55(1,1),1)
xyzzyaabi55=xyzzyaabi55+netot
xyzzyaabn55(:)=0
call xyzzyaahp55(0,xyzzyaacf55,xyzzyaacb55)
if(xyzzyaagy55)call wout(' MOVE ACCEPTED')
else
if(xyzzyaagy55)call wout(' MOVE REJECTED')
if(prob>0.d0)then
call timer('DRIFT_DIFFUSE',.false.)
call eval_local_energy(xyzzyaaer55,etot=xyzzyaaew55,isnan=isnan,isinf=&
&isinf)
if(xyzzyaagz55)then
if(isnan)then
call wout(' REJECTED ENERGY: Not a Number')
elseif(isinf)then
call wout(' REJECTED ENERGY: Diverges')
else
call wout(' REJECTED ENERGY: ',xyzzyaaew55)
endif
endif
if(isnan.or.isinf)then
prob=0.d0
else
xyzzyaacv55=0.d0
xyzzyaact55=0.d0
if(limdmc>=2)then
call timer('LIMIT_LOCAL_ENERGY',.true.)
do xyzzyaabd55=1,netot
xyzzyaacy55=xyzzyaada55(which_spin(xyzzyaabd55))
call wfn_loggrad(xyzzyaabd55,xyzzyaaer55,0,xyzzyaadh55,isnan=isnan,isi&
&nf=isinf)
xyzzyaacn55=real(xyzzyaadh55,dp)
if(isnan.or.isinf)exit
xyzzyaacv55=xyzzyaacv55+ddot(3,xyzzyaacn55(1),1,xyzzyaacn55(1),1)
if(nucleus_gf_mods)then
call get_eivecs(xyzzyaaer55)
call xyzzyaadr1(xyzzyaacn55,eivecs_scr(1,1,xyzzyaabd55,xyzzyaaer55),al&
&imit,xyzzyaads55,xyzzyaadt55,xyzzyaadu55,xyzzyaadw55)
endif
call xyzzyaado1(alimit,dt,xyzzyaacn55,xyzzyaacr55,xyzzyaacy55)
xyzzyaact55=xyzzyaact55+ddot(3,xyzzyaacr55(1),1,xyzzyaacr55(1),1)
enddo
if(isnan.or.isinf)prob=0.d0
call timer('LIMIT_LOCAL_ENERGY',.false.)
endif
endif
call timer('DRIFT_DIFFUSE',.true.)
endif
xyzzyaabn55(:)=xyzzyaabn55(:)+1
endif
do xyzzyaabd55=1,netot
call xyzzyaadm1(xyzzyaade55(xyzzyaabd55),prob,xyzzyaaci55,which_spin(x&
&yzzyaabd55))
enddo
call timer('DRIFT_DIFFUSE',.false.)
endif
call eval_local_energy(xyzzyaaeq55,etot=xyzzyaaev55,ecomps=xyzzyaafn55&
&,isnan=isnan,isinf=isinf)
if(xyzzyaagz55)then
if(isnan)then
call wout(' ENERGY: Not a Number')
elseif(isinf)then
call wout(' ENERGY: Diverges')
else
call wout(' ENERGY: ',xyzzyaaev55)
endif
endif
if(isnan.or.isinf)call errstop('MOVE_CONFIG','Floating-point exception&
& reported when calculating energy. Don''t know what to do, so stoppin&
&g. Please file a bug report.')
xyzzyaacs55=0.d0
xyzzyaacu55=0.d0
if(limdmc>=2)then
call timer('LIMIT_LOCAL_ENERGY',.true.)
do xyzzyaabd55=1,netot
xyzzyaacy55=xyzzyaada55(which_spin(xyzzyaabd55))
call wfn_loggrad(xyzzyaabd55,xyzzyaaeq55,0,xyzzyaadh55,isnan=isnan,isi&
&nf=isinf)
xyzzyaacn55=real(xyzzyaadh55,dp)
if(isnan.or.isinf)call errstop('MOVE_CONFIG','Floating-point exception&
& reported when calculating drift vector for local-energy limiting sch&
&eme. Don''t know what to do, so stopping. Please file a bug report.')
xyzzyaacu55=xyzzyaacu55+ddot(3,xyzzyaacn55(1),1,xyzzyaacn55(1),1)
if(nucleus_gf_mods)then
call get_eivecs(xyzzyaaeq55)
call xyzzyaadr1(xyzzyaacn55,eivecs_scr(1,1,xyzzyaabd55,xyzzyaaeq55),al&
&imit,xyzzyaads55,xyzzyaadt55,xyzzyaadu55,xyzzyaadw55)
endif
call xyzzyaado1(alimit,dt,xyzzyaacn55,xyzzyaacr55,xyzzyaacy55)
xyzzyaacs55=xyzzyaacs55+ddot(3,xyzzyaacr55(1),1,xyzzyaacr55(1),1)
enddo
call timer('LIMIT_LOCAL_ENERGY',.false.)
endif
xyzzyaabm55=sum(xyzzyaabn55(1:netot))
xyzzyaabk55=max(maxval(xyzzyaabn55(1:netot)),xyzzyaabk55)
xyzzyaabu55=xyzzyaabu55+xyzzyaabm55
if(ibran)then
call xyzzyaadp1(xyzzyaafa55,xyzzyaace55,xyzzyaafa55,xyzzyaacu55,xyzzya&
&acs55,xyzzyaaev55,xyzzyaaey55)
call xyzzyaadn1(dt,xyzzyaaek55)
if(use_tmove)then
xyzzyaaad55=-xyzzyaaek55*(0.5d0*(xyzzyaaex55+xyzzyaaey55)-xyzzyaaet55)
else
if(dmc_cfg_by_cfg)then
if(xyzzyaaci55)then
xyzzyaaad55=-xyzzyaaek55*((1.d0-0.5d0*prob)*xyzzyaaex55+0.5d0*prob*xyz&
&zyaaey55-xyzzyaaet55)
elseif(prob>0.d0)then
call xyzzyaadp1(xyzzyaafa55,xyzzyaace55,xyzzyaafa55,xyzzyaacv55,xyzzya&
&act55,xyzzyaaew55,xyzzyaaez55)
xyzzyaaad55=-xyzzyaaek55*((1.d0-0.5d0*prob)*xyzzyaaey55+0.5d0*prob*xyz&
&zyaaez55-xyzzyaaet55)
else
xyzzyaaad55=-xyzzyaaek55*(xyzzyaaey55-xyzzyaaet55)
endif
else
xyzzyaaad55=-xyzzyaaek55*(0.5d0*(xyzzyaaex55+xyzzyaaey55)-xyzzyaaet55)
endif
endif
wdmc=wdmc*exp_protect(xyzzyaaad55)
if(xyzzyaagz55)call wout(' BRANCHING FACTOR: ',wdmc)
else
wdmc=1.d0
xyzzyaaek55=dt
xyzzyaaey55=xyzzyaaev55
endif
if(expvals.and.xyzzyaacj55)call accumulate_expvals(xyzzyaaeq55,wdmc*re&
&al(corper_dmc,dp),tpdmc>0)
if(forces)then
xyzzyaafz55=0.d0
xyzzyaaga55=0.d0
if(have_ppots)then
if(.not.use_tmove)then
xyzzyaafz55=xyzzyaafn55(i_potinl)
else
call v_non_local_tmove(xyzzyaafz55,xyzzyaaga55)
endif
endif
call eval_local_forces(xyzzyaaeq55,xyzzyaaev55,xyzzyaafn55(i_kei),xyzz&
&yaafz55,xyzzyaaga55,xyzzyaafn55(xyzzyaafw55:xyzzyaafx55))
call forces_to_future(xyzzyaafn55(xyzzyaafw55:xyzzyaafx55),xyzzyaagq55&
&)
endif
if(use_future)then
if(.not.forces)xyzzyaagq55(1)=xyzzyaaev55
xyzzyaags55(1:xyzzyaagk55,0)=xyzzyaags55(1:xyzzyaagk55,0)+wdmc*xyzzyaa&
&gq55(1:xyzzyaagk55)
do xyzzyaagm55=1,xyzzyaago55
xyzzyaagn55=xyzzyaagl55-nint(real(xyzzyaagm55*xyzzyaagj55,dp)/real(xyz&
&zyaago55,dp))
if(xyzzyaagn55<=0)xyzzyaagn55=xyzzyaagn55+xyzzyaagj55
call xyzzyaacd1(xyzzyaagd55,xyzzyaagr55,xyzzyaagn55)
xyzzyaags55(1:xyzzyaagk55,xyzzyaagm55)=xyzzyaags55(1:xyzzyaagk55,xyzzy&
&aagm55)+wdmc*xyzzyaagr55(1:xyzzyaagk55)
enddo
if(forces)call future_to_forces(xyzzyaagr55,xyzzyaafn55(xyzzyaafw55:xy&
&zzyaafx55))
endif
xyzzyaafn55(xyzzyaaft55)=1.d0
xyzzyaafn55(xyzzyaafu55)=xyzzyaaek55
xyzzyaafm55=xyzzyaafm55+wdmc*xyzzyaafn55
if(use_tmove.and.have_ppots)then
call timer('T-MOVE',.true.)
do xyzzyaabd55=1,netot
if(.not.any(tmove_no_points(1:nitot,xyzzyaabd55)>=1))cycle
xyzzyaagb55=0.d0
do xyzzyaabf55=1,nitot
if(tmove_no_points(xyzzyaabf55,xyzzyaabd55)>0)then
do xyzzyaafy55=1,tmove_no_points(xyzzyaabf55,xyzzyaabd55)
if(tmove_t(xyzzyaafy55,xyzzyaabf55,xyzzyaabd55)>=0.d0)then
tmove_t(xyzzyaafy55,xyzzyaabf55,xyzzyaabd55)=0.d0
else
tmove_t(xyzzyaafy55,xyzzyaabf55,xyzzyaabd55)=-dt*tmove_t(xyzzyaafy55,x&
&yzzyaabf55,xyzzyaabd55)
xyzzyaagb55=xyzzyaagb55+tmove_t(xyzzyaafy55,xyzzyaabf55,xyzzyaabd55)
endif
enddo
endif
enddo
xyzzyaagb55=ranx()*(1.d0+xyzzyaagb55)
if(.not.tmove_fixgrid.and.xyzzyaagb55<=1.d0)then
xyzzyaafy55=0
else
xyzzyaagc55=1.d0
do xyzzyaabf55=1,nitot
do xyzzyaafy55=1,tmove_no_points(xyzzyaabf55,xyzzyaabd55)
xyzzyaagc55=xyzzyaagc55+tmove_t(xyzzyaafy55,xyzzyaabf55,xyzzyaabd55)
if(xyzzyaagb55<=xyzzyaagc55)exit
enddo
if(xyzzyaagb55<=xyzzyaagc55)exit
enddo
endif
if(xyzzyaafy55>1.or.(xyzzyaafy55==1.and..not.tmove_fixgrid))then
xyzzyaabp55=tmove_points(:,xyzzyaafy55,xyzzyaabf55,xyzzyaabd55)
call define_config_oneelec(xyzzyaabd55,xyzzyaaeq55,xyzzyaaes55,xyzzyaa&
&bp55,xyzzyaabo55(xyzzyaabd55))
call accept_move(xyzzyaaeq55,xyzzyaaes55)
xyzzyaacf55(:,xyzzyaabd55)=xyzzyaabp55(:)
xyzzyaabn55(xyzzyaabd55)=0
call xyzzyaahp55(xyzzyaabd55,xyzzyaacf55,xyzzyaacb55)
endif
enddo
call timer('T-MOVE',.false.)
endif
if(jasbuf.and.xyzzyaaao55==xyzzyaaaq55)then
call wfn_logval(xyzzyaaeq55,xyzzyaadi55)
endif
call xyzzyaace1(xyzzyaaeq55,xyzzyaagd55,xyzzyaaev55,xyzzyaaey55,wdmc,x&
&yzzyaabn55)
if(use_future)call xyzzyaacg1(xyzzyaagd55,xyzzyaagq55,xyzzyaagl55)
xyzzyaabj55=xyzzyaabj55+xyzzyaabi55
xyzzyaabh55=xyzzyaabh55+xyzzyaabg55
end subroutine xyzzyaahd55
subroutine xyzzyaahe55
implicit none
if(am_master)then
xyzzyaagi55=tcputime()-xyzzyaagh55
if(xyzzyaagi55>block_time)xyzzyaack55=.true.
endif
call qmpi_bcast(xyzzyaack55,'CHECK_BLOCKTIME_DMC','end_block')
end subroutine xyzzyaahe55
subroutine xyzzyaahf55(rion_old,rion_new,xyzzyaacf55,jacobian)
implicit none
real(dp),intent(in) :: rion_old(:,:),rion_new(:,:)
real(dp),intent(inout) :: xyzzyaacf55(:,:)
real(dp),intent(inout) :: jacobian
integer xyzzyaaaa61,xyzzyaaab61
real(dp) xyzzyaaac61(3),xyzzyaaad61(3),xyzzyaaae61,xyzzyaaaf61,xyzzyaa&
&ag61,xyzzyaaah61,xyzzyaaai61,xyzzyaaaj61,xyzzyaaak61(3),xyzzyaaal61(3&
&),xyzzyaaam61,xyzzyaaan61(3,3),xyzzyaaao61(3),xyzzyaaap61(3,3),xyzzya&
&aaq61
do xyzzyaaaa61=1,netot
xyzzyaaal61(:)=0.d0
xyzzyaaam61=0.d0
xyzzyaaan61(:,:)=0.d0
xyzzyaaao61(:)=0.d0
do xyzzyaaab61=1,nitot
xyzzyaaac61(:)=rion_new(:,xyzzyaaab61)-rion_old(:,xyzzyaaab61)
xyzzyaaad61=xyzzyaacf55(:,xyzzyaaaa61)-rion_old(:,xyzzyaaab61)
xyzzyaaaf61=1.d0/sum(xyzzyaaad61(:)**2)
xyzzyaaae61=sqrt(xyzzyaaaf61)
xyzzyaaag61=xyzzyaaaf61*xyzzyaaaf61
xyzzyaaah61=xyzzyaaag61*xyzzyaaae61
xyzzyaaai61=1.d0*xyzzyaaag61
xyzzyaaaj61=-4.d0*xyzzyaaah61
xyzzyaaak61(:)=xyzzyaaaj61*xyzzyaaad61(:)*xyzzyaaae61
xyzzyaaal61(:)=xyzzyaaal61(:)+xyzzyaaai61*xyzzyaaac61(:)
xyzzyaaam61=xyzzyaaam61+xyzzyaaai61
if(dmc_reweight_configs)then
xyzzyaaan61(:,1)=xyzzyaaan61(:,1)+xyzzyaaak61(:)*xyzzyaaac61(1)
xyzzyaaan61(:,2)=xyzzyaaan61(:,2)+xyzzyaaak61(:)*xyzzyaaac61(2)
xyzzyaaan61(:,3)=xyzzyaaan61(:,3)+xyzzyaaak61(:)*xyzzyaaac61(3)
xyzzyaaao61(:)=xyzzyaaao61(:)+xyzzyaaak61(:)
endif
enddo
xyzzyaaaq61=1.d0/xyzzyaaam61
xyzzyaacf55(:,xyzzyaaaa61)=xyzzyaacf55(:,xyzzyaaaa61)+xyzzyaaal61(:)*x&
&yzzyaaaq61
if(dmc_reweight_configs)then
xyzzyaaap61(:,1)=(/1,0,0/)+(xyzzyaaan61(:,1)-xyzzyaaao61(:)*xyzzyaaal6&
&1(1)*xyzzyaaaq61)*xyzzyaaaq61
xyzzyaaap61(:,2)=(/0,1,0/)+(xyzzyaaan61(:,2)-xyzzyaaao61(:)*xyzzyaaal6&
&1(2)*xyzzyaaaq61)*xyzzyaaaq61
xyzzyaaap61(:,3)=(/0,0,1/)+(xyzzyaaan61(:,3)-xyzzyaaao61(:)*xyzzyaaal6&
&1(3)*xyzzyaaaq61)*xyzzyaaaq61
jacobian=jacobian*(xyzzyaaap61(1,1)*(xyzzyaaap61(2,2)*xyzzyaaap61(3,3)&
&-xyzzyaaap61(2,3)*xyzzyaaap61(3,2)) +xyzzyaaap61(1,2)*(xyzzyaaap61(2,&
&3)*xyzzyaaap61(3,1)-xyzzyaaap61(2,1)*xyzzyaaap61(3,3))+xyzzyaaap61(1,&
&3)*(xyzzyaaap61(2,1)*xyzzyaaap61(3,2)-xyzzyaaap61(2,2)*xyzzyaaap61(3,&
&1)))
endif
enddo
end subroutine xyzzyaahf55
subroutine xyzzyaahg55(keep_in_memory,exceeded)
logical,intent(in) :: keep_in_memory,exceeded
integer xyzzyaaaa62,xyzzyaaab62
logical,save :: xyzzyaaac62=.false.
character(20) config_item(7),extra_item(4)
xyzzyaaaa62=5
config_item(1:5)=(/'RELE  ','LOGDET','ETOT  ','STOT  ','WDMC  '/)
if(noncoll_spin)then
xyzzyaaaa62=xyzzyaaaa62+1
config_item(xyzzyaaaa62)='SELE  '
endif
if(jasbuf.and.use_jastrow)then
xyzzyaaaa62=xyzzyaaaa62+1
config_item(xyzzyaaaa62)='VALJAS'
endif
xyzzyaaab62=3
extra_item(1:3)=(/'DMC_SAVED_STATE','RANDOM         ','GEOMETRY       &
&'/)
if(iaccum)then
xyzzyaaab62=xyzzyaaab62+1
extra_item(xyzzyaaab62)='REBLOCK_DATA   '
endif
tpdmc_config=tpdmc
lwdmc_config=lwdmc
growth_estimator_config=growth_estimator
call init_config_accumulation('DMC',xyzzyaaax55,config_item(1:xyzzyaaa&
&a62),extra_item(1:xyzzyaaab62))
if(dmc_twist_av.and..not.iaccum)then
dmcequil_steps_config=xyzzyaaau55-xyzzyaaaq55*ndmcave*iblock
else
dmcequil_steps_config=xyzzyaaau55
endif
dmcstats_steps_config=xyzzyaaav55
ebest_config=xyzzyaafa55
ebest_init_config=xyzzyaafb55
eref_config=xyzzyaaet55
dteff_ebest_init_config=xyzzyaaem55
dteff_best_config=xyzzyaael55
denominator_wt_config=xyzzyaaea55
log_pi_wt_config=xyzzyaady55
numer_expect_config=xyzzyaafg55
if(tpdmc>0)log_pi_wt_array_config=xyzzyaaeh55
if(growth_estimator)then
log_pi_wt_array2_config=xyzzyaaei55
log_pi_wt2_config=xyzzyaadz55
numerator_wt2_config=xyzzyaaec55
denominator_wt2_config=xyzzyaaeb55
endif
if(xyzzyaaax55>0)then
xyzzyaagd55=>xyzzyaaar1
do
call xyzzyaacp1(xyzzyaagd55)
xyzzyaagd55=>xyzzyaagd55%pt_next_config
if(associated(xyzzyaagd55,xyzzyaaar1))exit
enddo
nullify(xyzzyaagd55)
endif
call get_random_state(random_state_config)
if(am_master.and.iaccum)call save_reblock_data
call end_config_accumulation(exceeded)
if((iblock>=xyzzyaaam55.or.exceeded).and.rng_restart_safe)xyzzyaaac62=&
&.true.
select case(chkpoint_level)
case(-1)
if(exceeded)call write_configs(xyzzyaaac62)
case(0)
if(exceeded)then
call write_configs(xyzzyaaac62)
else
if(iblock>=xyzzyaaam55)then
if(((isvmc_dmc.or.isdmc_dmc).and.iaccum).or.isdmc_stats.or.isdmc_equil&
&.or.isdmc_old)then
call write_configs(xyzzyaaac62)
endif
endif
endif
case(1,2)
call write_configs(xyzzyaaac62)
end select
if(.not.keep_in_memory)call dismantle_configs
end subroutine xyzzyaahg55
subroutine xyzzyaahh55
integer xyzzyaaaa63,xyzzyaaab63
character(20) config_item(7),extra_item(4)
call timer('INIT_PRELIM_CONFIGS',.true.)
xyzzyaaaa63=5
config_item(1:5)=(/'RELE  ','LOGDET','ETOT  ','STOT  ','WDMC  '/)
if(noncoll_spin)then
xyzzyaaaa63=xyzzyaaaa63+1
config_item(xyzzyaaaa63)='SELE  '
endif
if(jasbuf.and.use_jastrow)then
xyzzyaaaa63=xyzzyaaaa63+1
config_item(xyzzyaaaa63)='VALJAS'
endif
xyzzyaaab63=3
extra_item(1:3)=(/'DMC_SAVED_STATE','RANDOM         ','GEOMETRY       &
&'/)
tpdmc_config=tpdmc
lwdmc_config=lwdmc
growth_estimator_config=growth_estimator
call init_config_accumulation('DMC',nconfig_prelim,config_item(1:xyzzy&
&aaaa63),extra_item(1:xyzzyaaab63))
call timer('INIT_PRELIM_CONFIGS',.false.)
end subroutine xyzzyaahh55
subroutine xyzzyaahi55
if(xyzzyaaax55>0.and.prelim_nconfig_written<nconfig_prelim)then
xyzzyaagd55=>xyzzyaaar1
do
call xyzzyaacp1(xyzzyaagd55)
xyzzyaagd55=>xyzzyaagd55%pt_next_config
prelim_nconfig_written=prelim_nconfig_written+1
if(associated(xyzzyaagd55,xyzzyaaar1) .or.prelim_nconfig_written==ncon&
&fig_prelim)exit
enddo
nullify(xyzzyaagd55)
endif
end subroutine xyzzyaahi55
subroutine xyzzyaahj55
logical :: xyzzyaaaa65=.false.
if(dmc_twist_av.and..not.iaccum)then
dmcequil_steps_config=xyzzyaaau55-xyzzyaaaq55*ndmcave*iblock
else
dmcequil_steps_config=xyzzyaaau55
endif
dmcstats_steps_config=xyzzyaaav55
ebest_config=xyzzyaafa55
ebest_init_config=xyzzyaafb55
eref_config=xyzzyaaet55
dteff_ebest_init_config=xyzzyaaem55
dteff_best_config=xyzzyaael55
denominator_wt_config=xyzzyaaea55
log_pi_wt_config=xyzzyaady55
numer_expect_config=xyzzyaafg55
if(tpdmc>0)log_pi_wt_array_config=xyzzyaaeh55
if(growth_estimator)then
log_pi_wt_array2_config=xyzzyaaei55
log_pi_wt2_config=xyzzyaadz55
numerator_wt2_config=xyzzyaaec55
denominator_wt2_config=xyzzyaaeb55
endif
call get_random_state(random_state_config)
call end_config_accumulation(.false.)
call write_configs(xyzzyaaaa65)
if(am_master)then
call wout()
call wordwrap('The initial configuration population for a large DMC ca&
&lculation has been written to config.out.')
call wout()
endif
end subroutine xyzzyaahj55
subroutine xyzzyaahk55
implicit none
real(dp) xyzzyaaaa66
logical xyzzyaaab66,xyzzyaaac66
if(xyzzyaaen55>=max_rec_attempts)call errstop_master('DMC','Have repea&
&tedly encountered persistence. Population explosion uncontrollable.')
if(xyzzyaabh1.or.xyzzyaabi1)then
call xyzzyaadd1('RECV',xyzzyaaax55)
xyzzyaabh1=.false.
xyzzyaabi1=.false.
endif
call xyzzyaacm1
if(iblock>2)then
xyzzyaaab66=.not.xyzzyaaep55
xyzzyaaac66=.true.
else
xyzzyaaab66=iblock==2
xyzzyaaac66=.false.
endif
if(am_master)then
call wout()
call wout('Catastrophe in block '//trim(i2s(iblock))//'.')
call wout('At move number       '//trim(i2s(xyzzyaaan55))//'.')
if(xyzzyaaac66)then
if(.not.expvals)then
call wout('Restarting from config.backup.')
else
call wout('Restarting from config.backup and expval.backup.')
endif
else
call wout('Restarting from config.in.')
endif
if(xyzzyaaab66)then
call wout('(Hence returning to start of previous block.)')
else
call wout('(Hence returning to start of current block.)')
endif
call wout()
endif
call xyzzyaahc55(xyzzyaaac66)
if(expvals.and.iaccum)then
if(iblock>2)then
call deallocate_expval
if(am_master)call read_expval(.true.)
call setup_expval(.true.)
else
call zero_expval
endif
endif
if(am_master.and..not.iaccum)xyzzyaafo55=xyzzyaafa55
do xyzzyaaaa55=1,xyzzyaaen55+1
xyzzyaaaa66=ranx()
xyzzyaaaa66=ranx_gaussian(0.5d0)
enddo
if(am_master.and.xyzzyaaab66.and.writeout_dmc_hist)call backtrack_hist&
&('DMC',xyzzyaaaq55)
if(xyzzyaaab66)iblock=iblock-1
xyzzyaaen55=xyzzyaaen55+1
xyzzyaaep55=.true.
end subroutine xyzzyaahk55
subroutine xyzzyaahl55
implicit none
logical xyzzyaaaa67
xyzzyaaaa67=.false.
if(pair_corr)then
if(xyzzyaacf55(1,fixed_particle)/=pcf_rfix(1).or.xyzzyaacf55(2,fixed_p&
&article)/=pcf_rfix(2).or.xyzzyaacf55(3,fixed_particle)/=pcf_rfix(3))x&
&yzzyaaaa67=.true.
elseif(pair_corr_sph)then
if(xyzzyaacf55(1,fixed_particle)/=pcfs_rfix(1).or.xyzzyaacf55(2,fixed_&
&particle)/=pcfs_rfix(2).or.xyzzyaacf55(3,fixed_particle)/=pcfs_rfix(3&
&))xyzzyaaaa67=.true.
else
call errstop('CHECK_PCF_RFIX','Called in inappropriate circumstances.'&
&)
endif
if(xyzzyaaaa67)then
call wout()
call wordwrap('Coordinates of fixed particle in rele array read from d&
&isk are not what was specified in input/expval.data.')
call errstop('CHECK_PCF_RFIX','Quitting.')
endif
end subroutine xyzzyaahl55
subroutine xyzzyaahm55(xyzzyaadm55)
implicit none
integer,intent(out) :: xyzzyaadm55(xyzzyaaax55)
integer xyzzyaaaa68,xyzzyaaab68,xyzzyaaac68
real(dp) dmult,xyzzyaaad68,xyzzyaaae68,xyzzyaaaf68(2),xyzzyaaag68(2),x&
&yzzyaaah68(2),xyzzyaaai68,xyzzyaaaj68,tmp1
real(dp),allocatable :: xyzzyaaak68(:),xyzzyaaal68(:),xyzzyaaam68(:)
type(configuration),pointer :: xyzzyaaan68
call timer('COMPUTE_MULTIPLICITIES',.true.)
if(xyzzyaabk1>1)call xyzzyaadf1
if(lwdmc)then
call qmpi_bcast(xyzzyaado55,'COMPUTE_MULTIPLICITIES','invtotweight')
xyzzyaaan68=>xyzzyaaar1
do xyzzyaaaa68=1,xyzzyaaax55
if(.not.iaccum.and.xyzzyaabx1<int(real(xyzzyaaam55,dp)*real(xyzzyaaaq5&
&5,dp)*dmc_equil_fixpop))then
dmult=xyzzyaaan68%pt_dmc%wdmc*targ_wt*xyzzyaado55
xyzzyaaan68%pt_dmc%wdmc=dmult
elseif(lwdmc_fixpop)then
dmult=xyzzyaaan68%pt_dmc%wdmc*targ_wt*xyzzyaado55
else
dmult=xyzzyaaan68%pt_dmc%wdmc
endif
if(dmult<wdmcmin)then
xyzzyaaac68=0
elseif(dmult>wdmcmax)then
if(dmult>=max_rep_int)call errstop('COMPUTE_MULTIPLICITIES','A compute&
&d multiplicity is bigger than the largest representable integer on th&
&is machine. Likely you have a really bad trial function.')
xyzzyaaac68=nint(dmult)
else
xyzzyaaac68=1
endif
xyzzyaadm55(xyzzyaaaa68)=xyzzyaaac68
if(xyzzyaaac68>1)xyzzyaaan68%pt_dmc%wdmc=xyzzyaaan68%pt_dmc%wdmc/xyzzy&
&aaac68
xyzzyaaan68=>xyzzyaaan68%pt_next_config
enddo
nullify(xyzzyaaan68)
xyzzyaaaa68=0
xyzzyaaab68=0
do while(count(xyzzyaadm55(xyzzyaaab68+1:xyzzyaaax55)==0)>1)
xyzzyaaaa68=xyzzyaaab68+1
do while(xyzzyaadm55(xyzzyaaaa68)>0)
xyzzyaaaa68=xyzzyaaaa68+1
enddo
nullify(xyzzyaaan68)
call xyzzyaacl1(xyzzyaaan68,xyzzyaaaa68)
xyzzyaaab68=xyzzyaaaa68+1
do while(xyzzyaadm55(xyzzyaaab68)>0)
xyzzyaaab68=xyzzyaaab68+1
enddo
nullify(xyzzyaage55)
call xyzzyaacl1(xyzzyaage55,xyzzyaaab68)
xyzzyaaad68=xyzzyaaan68%pt_dmc%wdmc+xyzzyaage55%pt_dmc%wdmc
if(xyzzyaaad68==0.d0)then
continue
elseif(ranx()<=xyzzyaaan68%pt_dmc%wdmc/xyzzyaaad68)then
xyzzyaadm55(xyzzyaaaa68)=1
xyzzyaaan68%pt_dmc%wdmc=xyzzyaaad68
else
xyzzyaadm55(xyzzyaaab68)=1
xyzzyaage55%pt_dmc%wdmc=xyzzyaaad68
endif
nullify(xyzzyaaan68,xyzzyaage55)
enddo
where(xyzzyaadm55(xyzzyaaab68+1:xyzzyaaax55)==0)xyzzyaadm55(xyzzyaaab6&
&8+1:xyzzyaaax55)=1
elseif(nc_dmc)then
allocate(xyzzyaaak68(xyzzyaaax55),xyzzyaaal68(xyzzyaaax55),xyzzyaaam68&
&(xyzzyaaax55),stat=xyzzyaabw1)
call check_alloc(xyzzyaabw1,'COMPUTE_MULTIPLICITIES','dmult_array')
xyzzyaaaf68(1:2)=0.d0
xyzzyaaan68=>xyzzyaaar1
do xyzzyaaaa68=1,xyzzyaaax55
xyzzyaaak68(xyzzyaaaa68)=xyzzyaaan68%pt_dmc%wdmc
xyzzyaaan68%pt_dmc%wdmc=1.d0
if(xyzzyaaak68(xyzzyaaaa68)>=1.d0)then
xyzzyaaal68(xyzzyaaaa68)=xyzzyaaak68(xyzzyaaaa68)-1.d0
xyzzyaaam68(xyzzyaaaa68)=0.d0
xyzzyaaaf68(1)=xyzzyaaaf68(1)+xyzzyaaal68(xyzzyaaaa68)
else
xyzzyaaal68(xyzzyaaaa68)=0.d0
xyzzyaaam68(xyzzyaaaa68)=1.d0-xyzzyaaak68(xyzzyaaaa68)
xyzzyaaaf68(2)=xyzzyaaaf68(2)+xyzzyaaam68(xyzzyaaaa68)
endif
xyzzyaaan68=>xyzzyaaan68%pt_next_config
enddo
nullify(xyzzyaaan68)
call qmpi_reduce(xyzzyaaaf68,xyzzyaaag68,mpi_sum,'COMPUTE_MULTIPLICITI&
&ES','sum_npm')
if(am_master)then
xyzzyaaae68=0.5d0*(xyzzyaaag68(1)+xyzzyaaag68(2))
if(xyzzyaaag68(1)>0.d0)then
xyzzyaaah68(1)=xyzzyaaae68/xyzzyaaag68(1)
else
xyzzyaaah68(1)=0.d0
endif
if(xyzzyaaag68(2)>0.d0)then
xyzzyaaah68(2)=xyzzyaaae68/xyzzyaaag68(2)
else
xyzzyaaah68(2)=0.d0
endif
endif
call qmpi_bcast(xyzzyaaah68,'COMPUTE_MULTIPLICITIES','wpm')
xyzzyaaai68=xyzzyaaah68(1)
xyzzyaaaj68=xyzzyaaah68(2)
xyzzyaadm55(:)=1
do xyzzyaaaa68=1,xyzzyaaax55
if(xyzzyaaak68(xyzzyaaaa68)>=1.d0)then
tmp1=xyzzyaaai68*xyzzyaaal68(xyzzyaaaa68)
if(tmp1>=max_rep_int)call errstop('COMPUTE_MULTIPLICITIES','A computed&
& multiplicity is bigger than the largest representable integer on thi&
&s machine. Likely you have a really bad trial function.')
if(ranx()<=mod(xyzzyaaai68*xyzzyaaal68(xyzzyaaaa68),1.d0))then
xyzzyaadm55(xyzzyaaaa68)=xyzzyaadm55(xyzzyaaaa68)+int(tmp1)+1
else
xyzzyaadm55(xyzzyaaaa68)=xyzzyaadm55(xyzzyaaaa68)+int(tmp1)
endif
else
if(ranx()<=mod(xyzzyaaaj68*xyzzyaaam68(xyzzyaaaa68),1.d0))xyzzyaadm55(&
&xyzzyaaaa68)=xyzzyaadm55(xyzzyaaaa68)-1
endif
enddo
deallocate(xyzzyaaak68,xyzzyaaal68,xyzzyaaam68)
call xyzzyaahn55
else
call qmpi_bcast(xyzzyaado55,'COMPUTE_MULTIPLICITIES','invtotweight')
xyzzyaaan68=>xyzzyaaar1
if(poprenorm.or.(.not.iaccum.and.xyzzyaabx1<int(real(xyzzyaaam55,dp)*r&
&eal(xyzzyaaaq55,dp)*dmc_equil_fixpop)))then
do xyzzyaaaa68=1,xyzzyaaax55
dmult=xyzzyaaan68%pt_dmc%wdmc*xyzzyaado55*targ_wt
tmp1=dmult+ranx()
if(tmp1>=max_rep_int)call errstop('COMPUTE_MULTIPLICITIES','A computed&
& multiplicity is bigger than the largest representable integer on thi&
&s machine. Likely you have a really bad trial function.')
xyzzyaadm55(xyzzyaaaa68)=int(tmp1)
xyzzyaaan68%pt_dmc%wdmc=1.d0
xyzzyaaan68=>xyzzyaaan68%pt_next_config
enddo
nullify(xyzzyaaan68)
call xyzzyaahn55
else
do xyzzyaaaa68=1,xyzzyaaax55
dmult=xyzzyaaan68%pt_dmc%wdmc
tmp1=dmult+ranx()
if(tmp1>=max_rep_int)call errstop('COMPUTE_MULTIPLICITIES','A computed&
& multiplicity is bigger than the largest representable integer on thi&
&s machine. Likely you have a really bad trial function.')
xyzzyaadm55(xyzzyaaaa68)=int(tmp1)
xyzzyaaan68%pt_dmc%wdmc=1.d0
xyzzyaaan68=>xyzzyaaan68%pt_next_config
enddo
nullify(xyzzyaaan68)
endif
endif
call timer('COMPUTE_MULTIPLICITIES',.false.)
end subroutine xyzzyaahm55
subroutine xyzzyaahn55
implicit none
integer xyzzyaaaa69,xyzzyaaab69,xyzzyaaac69,xyzzyaaad69,xyzzyaaae69,xy&
&zzyaaaf69,xyzzyaaag69
integer,allocatable :: xyzzyaaah69(:)
if(lwdmc)call errstop_master('RENORM','Cannot renormalize weighted DMC&
&.')
if(.not.allocated(xyzzyaabs1))call errstop_master('RENORM','RENORM cal&
&led without allocating TARGET_POP array. This is a bug.')
xyzzyaaae69=sum(xyzzyaadm55)
if(nnodes>1)then
xyzzyaaaa69=xyzzyaabl1-xyzzyaabq1(xyzzyaabm1)
if(xyzzyaaaa69<0)xyzzyaaaa69=xyzzyaaaa69+nnodes
if(xyzzyaabt1)then
allocate(xyzzyaaah69(0:xyzzyaabp1(xyzzyaabm1)-1),stat=xyzzyaabw1)
else
allocate(xyzzyaaah69(1),stat=xyzzyaabw1)
endif
call check_alloc(xyzzyaabw1,'RENORM','sum_mult_loc_array')
call mpi_gather(xyzzyaaae69,1,mpi_integer,xyzzyaaah69,1,mpi_integer,xy&
&zzyaaaa69,rg_comm(xyzzyaabm1,xyzzyaabo1),ierror)
call checkmpi(ierror,'Gathering sum_mult_loc in renorm.')
if(xyzzyaabt1)then
xyzzyaaaf69=sum(xyzzyaaah69(0:xyzzyaabp1(xyzzyaabm1)-1))
if(xyzzyaaaf69==0)call errstop('RENORM','Population died out before re&
&normalization.')
do while(xyzzyaaaf69>xyzzyaabs1(xyzzyaabm1))
xyzzyaaad69=1+int(ranx()*dble(xyzzyaaaf69))
if(xyzzyaaad69>xyzzyaaaf69)cycle
do xyzzyaaac69=0,xyzzyaabp1(xyzzyaabm1)-1
xyzzyaaad69=xyzzyaaad69-xyzzyaaah69(xyzzyaaac69)
if(xyzzyaaad69<=0)then
xyzzyaaah69(xyzzyaaac69)=xyzzyaaah69(xyzzyaaac69)-1
xyzzyaaaf69=xyzzyaaaf69-1
exit
endif
enddo
enddo
do while(xyzzyaaaf69<xyzzyaabs1(xyzzyaabm1))
xyzzyaaad69=1+int(ranx()*dble(xyzzyaaaf69))
if(xyzzyaaad69>xyzzyaaaf69)cycle
do xyzzyaaac69=0,xyzzyaabp1(xyzzyaabm1)-1
xyzzyaaad69=xyzzyaaad69-xyzzyaaah69(xyzzyaaac69)
if(xyzzyaaad69==0)then
xyzzyaaah69(xyzzyaaac69)=xyzzyaaah69(xyzzyaaac69)+1
xyzzyaaaf69=xyzzyaaaf69+1
exit
endif
enddo
enddo
xyzzyaaag69=xyzzyaaah69(xyzzyaaaa69)
do xyzzyaaac69=0,xyzzyaabp1(xyzzyaabm1)-1
if(xyzzyaaac69==xyzzyaaaa69)cycle
xyzzyaaab69=mod(xyzzyaabq1(xyzzyaabm1)+xyzzyaaac69,nnodes)
call qmpi_ssend(xyzzyaaah69(xyzzyaaac69),xyzzyaaab69,1,'RENORM','targe&
&t_pop_loc')
enddo
else
call qmpi_recv(xyzzyaaag69,xyzzyaabl1,1,'RENORM','target_pop_loc')
endif
deallocate(xyzzyaaah69)
else
xyzzyaaag69=xyzzyaabs1(1)
endif
do while(xyzzyaaae69>xyzzyaaag69)
xyzzyaaad69=1+int(ranx()*dble(xyzzyaaae69))
if(xyzzyaaad69>xyzzyaaae69)cycle
do xyzzyaadj55=1,xyzzyaaax55
xyzzyaaad69=xyzzyaaad69-xyzzyaadm55(xyzzyaadj55)
if(xyzzyaaad69<=0)then
xyzzyaadm55(xyzzyaadj55)=xyzzyaadm55(xyzzyaadj55)-1
xyzzyaaae69=xyzzyaaae69-1
exit
endif
enddo
enddo
do while(xyzzyaaae69<xyzzyaaag69)
xyzzyaaad69=1+int(ranx()*dble(xyzzyaaae69))
if(xyzzyaaad69>xyzzyaaae69)cycle
do xyzzyaadj55=1,xyzzyaaax55
xyzzyaaad69=xyzzyaaad69-xyzzyaadm55(xyzzyaadj55)
if(xyzzyaaad69<=0)then
xyzzyaadm55(xyzzyaadj55)=xyzzyaadm55(xyzzyaadj55)+1
xyzzyaaae69=xyzzyaaae69+1
exit
endif
enddo
enddo
end subroutine xyzzyaahn55
subroutine xyzzyaaho55(final)
implicit none
logical,intent(in) :: final
integer xyzzyaaaa70,xyzzyaaab70,xyzzyaaac70
real(dp) xyzzyaaad70(no_cols_qmc),xyzzyaaae70(no_cols_qmc),xyzzyaaaf70&
&(no_cols_qmc),xyzzyaaag70,xyzzyaaah70,xyzzyaaai70,xyzzyaaaj70
if(.not.am_master)return
if(.not.iaccum)return
if(dmc_twist_av)return
call timer('WRITE_DMC_STATUS',.true.)
if(final)then
xyzzyaaaa70=o
call wout(unit=xyzzyaaaa70)
call wout('Mixed estimators of the energies at the end of the run',uni&
&t=xyzzyaaaa70)
call wout('------------------------------------------------------',uni&
&t=xyzzyaaaa70)
call open_units(xyzzyaaab70,xyzzyaaac70)
if(xyzzyaaac70/=0)call errwarn('WRITE_DMC_STATUS','Unable to find free&
& i/o unitto delete dmc.status file.')
open(unit=xyzzyaaab70,file='dmc.status',iostat=xyzzyaaac70)
if(xyzzyaaac70/=0)call errwarn('WRITE_DMC_STATUS','Problem opening dmc&
&.status fordeletion after final DMC block.')
else
call open_units(xyzzyaaaa70,xyzzyaaac70)
if(xyzzyaaac70/=0)call errstop('WRITE_DMC_STATUS','Unable to find free&
& i/o unit.')
open(unit=xyzzyaaaa70,file='dmc.status',status='replace',iostat=xyzzya&
&aac70)
if(xyzzyaaac70/=0)call errstop('WRITE_DMC_STATUS','Problem opening dmc&
&.status')
call wout('Mixed estimators of the energies (current best estimate)',u&
&nit=xyzzyaaaa70)
call wout('--------------------------------------------------------',u&
&nit=xyzzyaaaa70)
endif
if(.not.isperiodic)then
xyzzyaafd55=1.d0
call wout('[All energies given in (au)]',unit=xyzzyaaaa70)
else if(model_system)then
xyzzyaafd55=1.d0
call wout('[All energies given in au per particle]',unit=xyzzyaaaa70)
else if(esupercell)then
xyzzyaafd55=real(npcells,dp)
call wout('[All energies given in au per simulation cell]',unit=xyzzya&
&aaa70)
else
xyzzyaafd55=1.d0
call wout('[All energies given in au per primitive cell]',unit=xyzzyaa&
&aa70)
endif
xyzzyaaff55=xyzzyaafe55*xyzzyaafd55
call get_reblocked_ave(xyzzyaaad70,xyzzyaaae70,xyzzyaaaf70)
if(esupercell)then
xyzzyaaad70(:)=xyzzyaaad70(:)*xyzzyaafd55
xyzzyaaae70(:)=xyzzyaaae70(:)*xyzzyaafd55
xyzzyaaaf70(:)=xyzzyaaaf70(:)*xyzzyaafd55
endif
if(xyzzyaaaf70(tagh_energy)==0.d0)then
call wordwrap('Reblocking not converged for ETOT. Too few data points?&
& Using unreblocked standard error.',xyzzyaaaa70)
elseif(xyzzyaaaf70(tagh_energy)>0.1*xyzzyaaae70(tagh_energy))then
call wordwrap('Bad reblock convergence for ETOT. Too few data points? &
&Standard error in standard error larger than 10%.',xyzzyaaaa70)
else
call wout('[All error bars obtained by reblocking]',unit=xyzzyaaaa70)
endif
call wout(unit=xyzzyaaaa70)
select case(trim(interaction))
case('ewald_mpc')
call wout('Total energy (Ewald)         = ',(/xyzzyaaad70(tagh_energy)&
&,xyzzyaaae70(tagh_energy)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaa&
&aa70)
call wout('Total energy (MPC)           = ',(/xyzzyaaad70(tagh_etotalt&
&),xyzzyaaae70(tagh_etotalt)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzy&
&aaaa70)
case('mpc_ewald')
call wout('Total energy (MPC)           = ',(/xyzzyaaad70(tagh_energy)&
&,xyzzyaaae70(tagh_energy)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaa&
&aa70)
call wout('Total energy (Ewald)         = ',(/xyzzyaaad70(tagh_etotalt&
&),xyzzyaaae70(tagh_etotalt)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzy&
&aaaa70)
case default
call wout('Total energy                 = ',(/xyzzyaaad70(tagh_energy)&
&,xyzzyaaae70(tagh_energy)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaa&
&aa70)
end select
call wout('Kinetic energy (TI)          = ',(/xyzzyaaad70(tagh_t),xyzz&
&yaaae70(tagh_t)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa70)
call wout('Kinetic energy (KEI)         = ',(/xyzzyaaad70(tagh_k),xyzz&
&yaaae70(tagh_k)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa70)
call wout('Kinetic energy (FISQ)        = ',(/xyzzyaaad70(tagh_fisq),x&
&yzzyaaae70(tagh_fisq)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa70&
&)
if(tagh_ewald>0)call wout('e-e interac. (Ewald/Coulomb) = ',(/xyzzyaaa&
&d70(tagh_ewald),xyzzyaaae70(tagh_ewald)/),rfmt='(f21.12)',rsep=' +/- &
&',unit=xyzzyaaaa70)
call wout('e-i interaction (local)      = ',(/xyzzyaaad70(tagh_local),&
&xyzzyaaae70(tagh_local)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa&
&70)
if(tagh_nonlocal>0)call wout('e-i interaction (nonlocal)   = ',(/xyzzy&
&aaad70(tagh_nonlocal),xyzzyaaae70(tagh_nonlocal)/),rfmt='(f21.12)',rs&
&ep=' +/- ',unit=xyzzyaaaa70)
if(interaction_mpc_present)then
call wout('Short-range part of MPC      = ',(/xyzzyaaad70(tagh_short),&
&xyzzyaaae70(tagh_short)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa&
&70)
call wout('Long-range part of MPC       = ',(/xyzzyaaad70(tagh_long),x&
&yzzyaaae70(tagh_long)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa70&
&)
endif
if(hartree_xc)then
call wout('Hartree energy from MPC      = ',mpc_hartree*xyzzyaafd55,rf&
&mt='(f21.12)',unit=xyzzyaaaa70)
call wout('XC energy from MPC           = ',xyzzyaaad70(tagh_short)+xy&
&zzyaaad70(tagh_long)-mpc_hartree*xyzzyaafd55,rfmt='(f21.12)',unit=xyz&
&zyaaaa70)
endif
if(int_sf)then
call eval_int_sf(xyzzyaafh55,xyzzyaafi55)
call wout('Hartree energy from SF       = ',xyzzyaafh55*xyzzyaafd55,rf&
&mt='(f21.12)',unit=xyzzyaaaa70)
call wout('XC energy from SF            = ',xyzzyaafi55*xyzzyaafd55,rf&
&mt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total e-e energy from SF     = ',(xyzzyaafh55+xyzzyaafi55)*&
&xyzzyaafd55,rfmt='(f21.12)',unit=xyzzyaaaa70)
endif
if(have_veep)then
if(isperiodic)then
call wout('Core polarization (total)    = ',(/xyzzyaaad70(tagh_cppei),&
&xyzzyaaae70(tagh_cppei)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa&
&70)
else
call wout('Core polarization (e-i term) = ',(/xyzzyaaad70(tagh_cppei),&
&xyzzyaaae70(tagh_cppei)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa&
&70)
call wout('Core polarization (e term)   = ',(/xyzzyaaad70(tagh_cppe),x&
&yzzyaaae70(tagh_cppe)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa70&
&)
call wout('Core polarization (e-e term) = ',(/xyzzyaaad70(tagh_cppee),&
&xyzzyaaae70(tagh_cppee)/),rfmt='(f21.12)',rsep=' +/- ',unit=xyzzyaaaa&
&70)
endif
endif
if(constant_energy/=0.d0)call wout('Constant energy contribs.    = ',c&
&onstant_energy*xyzzyaaff55,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout(unit=xyzzyaaaa70)
if(relativistic)then
call wout('Relativistic energy data are also present (see output of re&
&block).',unit=xyzzyaaaa70)
call wout(unit=xyzzyaaaa70)
endif
if(growth_estimator)then
call wout('Growth estimator of total energy = ',xyzzyaafc55*xyzzyaaff5&
&5,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout(unit=xyzzyaaaa70)
endif
if(esupercell)then
call wout("Dump of raw reblock data (per primitive cell)",unit=xyzzyaa&
&aa70)
call wout('---------------------------------------------',unit=xyzzyaa&
&aa70)
else
call wout("Dump of raw reblock data",unit=xyzzyaaaa70)
call wout('------------------------',unit=xyzzyaaaa70)
endif
call wout("Number of data points collected = "//trim(i2s(reblock_nstep&
&)),unit=xyzzyaaaa70)
call reblock_dump(tagh_energy,"energy",unit=xyzzyaaaa70)
call wout(unit=xyzzyaaaa70)
call popstats_dump(dt,targ_wt,unit=xyzzyaaaa70)
if(expvals.and.finite_size_corr)then
xyzzyaaag70=xyzzyaaad70(tagh_energy)
call finite_size_corr_xc
xyzzyaaai70=xc_corr*xyzzyaaff55
xyzzyaaaj70=ke_corr*xyzzyaaff55
call wout('Finite size correction data:',unit=xyzzyaaaa70)
if(chi_squared_sf>0.d0)call wout('Chi squared in structure factor fit &
&',chi_squared_sf,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total XC correction (dXC)        =  ',xyzzyaaai70,rfmt='(f2&
&1.12)',unit=xyzzyaaaa70)
call wout('Total KE correction (dKE)        =  ',xyzzyaaaj70,rfmt='(f2&
&1.12)',unit=xyzzyaaaa70)
select case(trim(interaction))
case('mpc')
call wout('Total energy (MPC)               =  ',xyzzyaaag70,rfmt='(f2&
&1.12)',unit=xyzzyaaaa70)
call wout('Total energy (MPC) + dKE         =  ',xyzzyaaag70+xyzzyaaaj&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
case('ewald_mpc')
xyzzyaaah70=xyzzyaaad70(tagh_etotalt)
call wout('Total energy (Ewald)             =  ',xyzzyaaag70,rfmt='(f2&
&1.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dXC + dKE =  ',xyzzyaaag70+xyzzyaaaj&
&70+xyzzyaaai70,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dXC       =  ',xyzzyaaag70+xyzzyaaai&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dKE       =  ',xyzzyaaag70+xyzzyaaaj&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total energy (MPC)               =  ',xyzzyaaah70,rfmt='(f2&
&1.12)',unit=xyzzyaaaa70)
call wout('Total energy (MPC) + dKE         =  ',xyzzyaaah70+xyzzyaaaj&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
case('mpc_ewald')
xyzzyaaah70=xyzzyaaad70(tagh_etotalt)
call wout('Total energy (MPC)               =  ',xyzzyaaag70,rfmt='(f2&
&1.12)',unit=xyzzyaaaa70)
call wout('Total energy (MPC) + dKE         =  ',xyzzyaaag70+xyzzyaaaj&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald)             =  ',xyzzyaaah70,rfmt='(f2&
&1.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dXC + dKE =  ',xyzzyaaah70+xyzzyaaai&
&70+xyzzyaaaj70,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dXC       =  ',xyzzyaaah70+xyzzyaaai&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dKE       =  ',xyzzyaaah70+xyzzyaaaj&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
case default
call wout('Total energy (Ewald)             =  ',xyzzyaaag70,rfmt='(f2&
&1.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dXC + dKE =  ',xyzzyaaag70+xyzzyaaai&
&70+xyzzyaaaj70,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dXC       =  ',xyzzyaaag70+xyzzyaaai&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
call wout('Total energy (Ewald) + dKE       =  ',xyzzyaaag70+xyzzyaaaj&
&70,rfmt='(f21.12)',unit=xyzzyaaaa70)
end select
call wout(unit=xyzzyaaaa70)
endif
call wout(repeat('=',73),unit=xyzzyaaaa70)
if(.not.final)then
close(xyzzyaaaa70)
open_unit(xyzzyaaaa70)=.false.
else
close(xyzzyaaab70,status='delete')
open_unit(xyzzyaaab70)=.false.
endif
call timer('WRITE_DMC_STATUS',.false.)
end subroutine xyzzyaaho55
subroutine xyzzyaahp55(xyzzyaabd55,xyzzyaacf55,xyzzyaacb55)
implicit none
integer,intent(in) :: xyzzyaabd55
real(dp),intent(in) :: xyzzyaacf55(3,netot)
real(dp),intent(inout) :: xyzzyaacb55
integer xyzzyaaaa71,xyzzyaaab71,xyzzyaaac71
real(dp) xyzzyaaad71(3),xyzzyaaae71,xyzzyaaaf71
if(isperiodic)return
xyzzyaaad71=0.d0
if(homogeneous_system)then
xyzzyaaae71=0.d0
xyzzyaaac71=0
do xyzzyaaaa71=1,nspin
xyzzyaaaf71=pmass(xyzzyaaaa71)
xyzzyaaae71=xyzzyaaae71+nele(xyzzyaaaa71)*xyzzyaaaf71
do xyzzyaaab71=1,nele(xyzzyaaaa71)
xyzzyaaac71=xyzzyaaac71+1
xyzzyaaad71(1:dimensionality)=xyzzyaaad71(1:dimensionality) +xyzzyaaaf&
&71*xyzzyaacf55(1:dimensionality,xyzzyaaac71)
enddo
enddo
xyzzyaaad71(1:dimensionality)=xyzzyaaad71(1:dimensionality)/xyzzyaaae7&
&1
endif
if(xyzzyaabd55/=0)then
xyzzyaacb55=max(xyzzyaacb55,sum((xyzzyaacf55(1:dimensionality,xyzzyaab&
&d55)-xyzzyaaad71(1:dimensionality))**2))
else
do xyzzyaaac71=1,netot
xyzzyaacb55=max(xyzzyaacb55,sum((xyzzyaacf55(1:dimensionality,xyzzyaaa&
&c71)-xyzzyaaad71(1:dimensionality))**2))
enddo
endif
end subroutine xyzzyaahp55
end subroutine dmc_main
end module slaarnaaj

module slaarnacq
use dsp
implicit none
private
public setup_wfdet,get_wfdet_rmax,wfdet,copy_orb_to_det,copy_rmap_orb_&
&to_det,setup_wfdet_params,finish_wfdet_params,get_wfdet_params,put_wf&
&det_params,enumerate_plot_wfdet,query_plot_wfdet,get_plot_wfdet,finis&
&h_plot_wfdet,get_wfdet_orbmask,wfdet_write_mdet_casl,wfdet_mdet_to_cm&
&det,wfdet_detcoef_affect_orbs
public wfdet_orbmap,wfdet_orbmask,wfdet_norb,wfdet_orbmask_scr1,wfdet_&
&orbmask_scr2
integer,parameter :: xyzzyaaaa1=4
integer xyzzyaaab1,xyzzyaaac1(xyzzyaaaa1)
integer,allocatable :: xyzzyaaad1(:)
integer :: xyzzyaaae1=0
logical :: xyzzyaaaf1=.false.,xyzzyaaag1=.false.,xyzzyaaah1=.false.,xy&
&zzyaaai1=.false.,xyzzyaaaj1=.false.,xyzzyaaak1=.false.,xyzzyaaal1=.fa&
&lse.
integer wfdet_norb,xyzzyaaam1,xyzzyaaan1
integer,allocatable :: wfdet_orbmap(:,:,:)
logical,allocatable :: wfdet_orbmask(:,:)
logical,allocatable :: wfdet_orbmask_scr1(:,:),wfdet_orbmask_scr2(:,:)
integer xyzzyaaao1,xyzzyaaap1,xyzzyaaaq1,xyzzyaaar1
integer,allocatable :: xyzzyaaas1(:,:,:),xyzzyaaat1(:),xyzzyaaau1(:,:)&
&,xyzzyaaav1(:),xyzzyaaaw1(:,:),xyzzyaaax1(:,:),xyzzyaaay1(:,:,:),xyzz&
&yaaaz1(:,:,:),xyzzyaaba1(:),xyzzyaabb1(:,:),xyzzyaabc1(:,:),xyzzyaabd&
&1(:,:,:),xyzzyaabe1(:,:)
logical xyzzyaabf1
logical,allocatable :: xyzzyaabg1(:,:)
real(dp),allocatable :: xyzzyaabh1(:),xyzzyaabi1(:,:),xyzzyaabj1(:)
contains
subroutine setup_wfdet()
use slaarnabg,    only : atom_basis_type
use run_control, only : errstop,check_alloc
use store,       only : nspin,ndet,nemax,heg_orbtype
implicit none
integer xyzzyaaaa2
select case(trim(atom_basis_type))
case('none')
xyzzyaaae1=0
case('plane-wave')
xyzzyaaae1=1
case('gaussian')
xyzzyaaae1=2
case('numerical')
xyzzyaaae1=3
case('blip')
xyzzyaaae1=4
case('non_int_he')
xyzzyaaae1=5
case('h2')
xyzzyaaae1=5
case('h3plus')
xyzzyaaae1=5
case('slater-type')
xyzzyaaae1=6
case('dimer')
xyzzyaaae1=7
case default
call errstop('SETUP_WFDET','Confused.')
end select
xyzzyaaah1=any(heg_orbtype<0)
xyzzyaaaf1=any(heg_orbtype==1)
xyzzyaaag1=any(heg_orbtype==2)
xyzzyaaaj1=any(heg_orbtype==3)
xyzzyaaak1=any(heg_orbtype==4).or.any(heg_orbtype==5).or.any(heg_orbty&
&pe==6)
xyzzyaaai1=any(heg_orbtype>100)
xyzzyaaal1=xyzzyaaah1.or.xyzzyaaaf1.or.xyzzyaaag1.or.xyzzyaaaj1.or.xyz&
&zyaaai1.or.xyzzyaaak1
allocate(wfdet_orbmap(nemax,nspin,ndet),stat=xyzzyaaaa2)
call check_alloc(xyzzyaaaa2,'SETUP_WFDET','wfdet_orbmap')
call xyzzyaabl1(wfdet_norb,wfdet_orbmap)
call xyzzyaabm1()
allocate(wfdet_orbmask(wfdet_norb,nspin),stat=xyzzyaaaa2)
call check_alloc(xyzzyaaaa2,'SETUP_WFDET','wfdet_orbmask')
call get_wfdet_orbmask(wfdet_norb,wfdet_orbmap,wfdet_orbmask)
end subroutine setup_wfdet
subroutine xyzzyaabk1(rvec,spin,jspin,norb,orbmask,val,fsd,orbval,orbg&
&rad,orblap,eevecs,orbsderivs,orb_m,orb_rmap)
use slaarnaac,           only : atomic_orb_eval
use slaarnaae,           only : blip_orb_eval
use slaarnaak,           only : dimer_orb_eval
use slaarnaap,            only : expot_is_finite
use slaarnaar,          only : expot_orb_eval
use slaarnaas,        only : free_biex1_orb_eval,free_crystal_orb_eval&
&,free_fluid_orb_eval,free_pairing_orb_eval
use slaarnaau,        only : gaussian_orb_eval
use slaarnacb,           only : pw_orb_eval
use run_control,      only : errstop,timer
use slaarnacl,         only : sto_orb_eval
use slaarnaci,              only : sdw_orb_eval
use slaarnack,  only : special_orb_eval
use store,            only : real1_complex2,use_expot,sparse,netot
implicit none
integer,intent(in) :: spin,jspin,norb
integer,optional,intent(inout) :: orb_m,orb_rmap(xyzzyaaam1)
real(dp),intent(in) :: rvec(3),eevecs(4,netot)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),optional,intent(inout) :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd,orbmask(norb)
call timer('WFDET',.true.)
if(use_expot)then
if(.not.expot_is_finite(rvec,jspin))then
if(val)orbval=0.d0
if(fsd)then
orbgrad=0.d0
orblap=0.d0
if(present(orbsderivs))orbsderivs=0.d0
endif
if(sparse)then
orb_m=0
orb_rmap=0
endif
call timer('WFDET',.false.)
return
endif
endif
select case(xyzzyaaae1)
case(0)
case(1)
call pw_orb_eval(rvec,jspin,xyzzyaaam1,orbmask,val,fsd,orbval,orbgrad,&
&orblap,orbsderivs)
case(2)
call gaussian_orb_eval(rvec,jspin,xyzzyaaam1,orbmask,val,fsd,orbval,or&
&bgrad,orblap,orbsderivs)
case(3)
call atomic_orb_eval(rvec,jspin,xyzzyaaam1,orbmask,val,fsd,orbval,orbg&
&rad,orblap,orbsderivs)
case(4)
call blip_orb_eval(rvec,jspin,xyzzyaaam1,orbmask,val,fsd,orbval,orbgra&
&d,orblap,orbsderivs,orb_m,orb_rmap)
case(5)
call special_orb_eval(rvec,jspin,xyzzyaaam1,orbmask,val,fsd,orbval,orb&
&grad,orblap,orbsderivs)
case(6)
call sto_orb_eval(rvec,jspin,xyzzyaaam1,orbmask,val,fsd,orbval,orbgrad&
&,orblap,orbsderivs)
case(7)
call dimer_orb_eval(rvec,jspin,xyzzyaaam1,orbmask,val,fsd,orbval,orbgr&
&ad,orblap,orbsderivs)
case default
call errstop('WFDET','Dazed and confused.')
end select
if(xyzzyaaal1)then
if(present(orbsderivs))then
if(xyzzyaaah1)call free_pairing_orb_eval(eevecs,jspin,norb,xyzzyaaan1,&
&orbmask(xyzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzya&
&aam1+1,1),orblap(xyzzyaaam1+1,1),orbsderivs(1,xyzzyaaam1+1,1))
if(xyzzyaaaf1)call free_fluid_orb_eval(rvec,jspin,norb,xyzzyaaan1,orbm&
&ask(xyzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaam1&
&+1,1),orblap(xyzzyaaam1+1,1),orbsderivs(1,xyzzyaaam1+1,1))
if(xyzzyaaag1)call free_crystal_orb_eval(rvec,jspin,norb,xyzzyaaan1,or&
&bmask(xyzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaa&
&m1+1,1),orblap(xyzzyaaam1+1,1),orbsderivs(1,xyzzyaaam1+1,1))
if(xyzzyaaak1)call free_biex1_orb_eval(rvec,jspin,norb,xyzzyaaan1,orbm&
&ask(xyzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaam1&
&+1,1),orblap(xyzzyaaam1+1,1),orbsderivs(1,xyzzyaaam1+1,1))
if(xyzzyaaaj1)call sdw_orb_eval(rvec,spin,norb,xyzzyaaan1,orbmask(xyzz&
&yaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaam1+1,1),or&
&blap(xyzzyaaam1+1,1),orbsderivs(1,xyzzyaaam1+1,1))
if(xyzzyaaai1)call expot_orb_eval(rvec,jspin,norb,xyzzyaaan1,orbmask(x&
&yzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaam1+1,1)&
&,orblap(xyzzyaaam1+1,1),orbsderivs(1,xyzzyaaam1+1,1))
else
if(xyzzyaaah1)call free_pairing_orb_eval(eevecs,jspin,norb,xyzzyaaan1,&
&orbmask(xyzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzya&
&aam1+1,1),orblap(xyzzyaaam1+1,1))
if(xyzzyaaaf1)call free_fluid_orb_eval(rvec,jspin,norb,xyzzyaaan1,orbm&
&ask(xyzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaam1&
&+1,1),orblap(xyzzyaaam1+1,1))
if(xyzzyaaag1)call free_crystal_orb_eval(rvec,jspin,norb,xyzzyaaan1,or&
&bmask(xyzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaa&
&m1+1,1),orblap(xyzzyaaam1+1,1))
if(xyzzyaaak1)call free_biex1_orb_eval(rvec,jspin,norb,xyzzyaaan1,orbm&
&ask(xyzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaam1&
&+1,1),orblap(xyzzyaaam1+1,1))
if(xyzzyaaaj1)call sdw_orb_eval(rvec,spin,norb,xyzzyaaan1,orbmask(xyzz&
&yaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaam1+1,1),or&
&blap(xyzzyaaam1+1,1))
if(xyzzyaaai1)call expot_orb_eval(rvec,jspin,norb,xyzzyaaan1,orbmask(x&
&yzzyaaam1+1),val,fsd,orbval(xyzzyaaam1+1,1),orbgrad(1,xyzzyaaam1+1,1)&
&,orblap(xyzzyaaam1+1,1))
endif
endif
call timer('WFDET',.false.)
end subroutine xyzzyaabk1
subroutine setup_wfdet_params(nparam)
use slaarnaas,  only : setup_freeorb_params
use slaarnabg,   only : model_system
use slaarnaau,  only : setup_gwfmolorb_params
use slaarnabu,    only : setup_orbmod_params
use store,      only : use_orbmods,opt_orbitals
use slaarnacl,   only : setup_stowf_params
use run_control,only : check_alloc
implicit none
integer,intent(out) :: nparam
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4
nparam=0
xyzzyaaab1=0
if(.not.opt_orbitals)return
xyzzyaaac1=0
if(model_system)call setup_freeorb_params(xyzzyaaac1(1))
if(use_orbmods)then
select case(xyzzyaaae1)
case(2)
call setup_gwfmolorb_params(xyzzyaaac1(3))
case(3)
call setup_orbmod_params(xyzzyaaac1(2))
case(6)
call setup_stowf_params(xyzzyaaac1(4))
end select
endif
nparam=sum(xyzzyaaac1,xyzzyaaac1>0)
xyzzyaaab1=nparam
allocate(xyzzyaaad1(xyzzyaaab1),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'SETUP_WFDET_PARAMS','wfdet_param_sec')
xyzzyaaab4=0
do xyzzyaaac4=1,xyzzyaaaa1
if(xyzzyaaac1(xyzzyaaac4)<1)cycle
xyzzyaaaa4=xyzzyaaab4+1
xyzzyaaab4=xyzzyaaab4+xyzzyaaac1(xyzzyaaac4)
xyzzyaaad1(xyzzyaaaa4:xyzzyaaab4)=xyzzyaaac4
enddo
end subroutine setup_wfdet_params
subroutine finish_wfdet_params
use slaarnaas,only : finish_freeorb_params
use slaarnabg, only : model_system
use slaarnaau,only : finish_gwfmolorb_params
use slaarnabu,  only : finish_orbmod_params
use store,    only : use_orbmods,opt_orbitals
use slaarnacl, only : finish_stowf_params
implicit none
if(.not.opt_orbitals)return
if(model_system)call finish_freeorb_params
if(use_orbmods)then
select case(xyzzyaaae1)
case(2)
call finish_gwfmolorb_params
case(3)
call finish_orbmod_params
case(6)
call finish_stowf_params
end select
endif
deallocate(xyzzyaaad1)
end subroutine finish_wfdet_params
subroutine get_wfdet_params(params,has_lolim,lolim,has_hilim,hilim,is_&
&shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map,lab&
&el)
use slaarnaas,only : get_freeorb_params
use slaarnaau,only : get_gwfmolorb_params
use slaarnabu,  only : get_orbmod_params
use slaarnacl, only : get_stowf_params
implicit none
real(dp),intent(inout) :: params(xyzzyaaab1),lolim(xyzzyaaab1),hilim(x&
&yzzyaaab1)
logical,intent(inout) :: has_lolim(xyzzyaaab1),has_hilim(xyzzyaaab1),i&
&s_shallow(xyzzyaaab1),is_redundant(xyzzyaaab1),is_linear(xyzzyaaab1),&
&is_loglinear(xyzzyaaab1),has_aderiv(xyzzyaaab1),affect_map(xyzzyaaab1&
&,xyzzyaaab1)
character(2),intent(inout) :: label(xyzzyaaab1)
integer xyzzyaaaa6,xyzzyaaab6
label='OP'
xyzzyaaab6=0
if(xyzzyaaac1(1)>0)then
xyzzyaaaa6=xyzzyaaab6+1
xyzzyaaab6=xyzzyaaab6+xyzzyaaac1(1)
call get_freeorb_params(params(xyzzyaaaa6:xyzzyaaab6),has_lolim(xyzzya&
&aaa6:xyzzyaaab6),lolim(xyzzyaaaa6:xyzzyaaab6),has_hilim(xyzzyaaaa6:xy&
&zzyaaab6),hilim(xyzzyaaaa6:xyzzyaaab6),is_shallow(xyzzyaaaa6:xyzzyaaa&
&b6),is_redundant(xyzzyaaaa6:xyzzyaaab6),is_linear(xyzzyaaaa6:xyzzyaaa&
&b6),is_loglinear(xyzzyaaaa6:xyzzyaaab6),has_aderiv(xyzzyaaaa6:xyzzyaa&
&ab6),affect_map(xyzzyaaaa6:xyzzyaaab6,xyzzyaaaa6:xyzzyaaab6),label(xy&
&zzyaaaa6:xyzzyaaab6))
endif
if(xyzzyaaac1(2)>0)then
xyzzyaaaa6=xyzzyaaab6+1
xyzzyaaab6=xyzzyaaab6+xyzzyaaac1(2)
call get_orbmod_params(params(xyzzyaaaa6:xyzzyaaab6),has_lolim(xyzzyaa&
&aa6:xyzzyaaab6),lolim(xyzzyaaaa6:xyzzyaaab6),has_hilim(xyzzyaaaa6:xyz&
&zyaaab6),hilim(xyzzyaaaa6:xyzzyaaab6),is_shallow(xyzzyaaaa6:xyzzyaaab&
&6),is_redundant(xyzzyaaaa6:xyzzyaaab6),is_linear(xyzzyaaaa6:xyzzyaaab&
&6),is_loglinear(xyzzyaaaa6:xyzzyaaab6),has_aderiv(xyzzyaaaa6:xyzzyaaa&
&b6),affect_map(xyzzyaaaa6:xyzzyaaab6,xyzzyaaaa6:xyzzyaaab6),label(xyz&
&zyaaaa6:xyzzyaaab6))
endif
if(xyzzyaaac1(3)>0)then
xyzzyaaaa6=xyzzyaaab6+1
xyzzyaaab6=xyzzyaaab6+xyzzyaaac1(3)
call get_gwfmolorb_params(params(xyzzyaaaa6:xyzzyaaab6),has_lolim(xyzz&
&yaaaa6:xyzzyaaab6),lolim(xyzzyaaaa6:xyzzyaaab6),has_hilim(xyzzyaaaa6:&
&xyzzyaaab6),hilim(xyzzyaaaa6:xyzzyaaab6),is_shallow(xyzzyaaaa6:xyzzya&
&aab6),is_redundant(xyzzyaaaa6:xyzzyaaab6),is_linear(xyzzyaaaa6:xyzzya&
&aab6),is_loglinear(xyzzyaaaa6:xyzzyaaab6),has_aderiv(xyzzyaaaa6:xyzzy&
&aaab6),affect_map(xyzzyaaaa6:xyzzyaaab6,xyzzyaaaa6:xyzzyaaab6),label(&
&xyzzyaaaa6:xyzzyaaab6))
endif
if(xyzzyaaac1(4)>0)then
xyzzyaaaa6=xyzzyaaab6+1
xyzzyaaab6=xyzzyaaab6+xyzzyaaac1(4)
call get_stowf_params(params(xyzzyaaaa6:xyzzyaaab6),has_lolim(xyzzyaaa&
&a6:xyzzyaaab6),lolim(xyzzyaaaa6:xyzzyaaab6),has_hilim(xyzzyaaaa6:xyzz&
&yaaab6),hilim(xyzzyaaaa6:xyzzyaaab6),is_shallow(xyzzyaaaa6:xyzzyaaab6&
&),is_redundant(xyzzyaaaa6:xyzzyaaab6),is_linear(xyzzyaaaa6:xyzzyaaab6&
&),is_loglinear(xyzzyaaaa6:xyzzyaaab6),has_aderiv(xyzzyaaaa6:xyzzyaaab&
&6),affect_map(xyzzyaaaa6:xyzzyaaab6,xyzzyaaaa6:xyzzyaaab6),label(xyzz&
&yaaaa6:xyzzyaaab6))
endif
end subroutine get_wfdet_params
subroutine put_wfdet_params(params,ignore,iparam_buffer,prestore,bad_p&
&arams)
use slaarnaas,only : put_freeorb_params
use slaarnaau,only : put_gwfmolorb_params
use slaarnabu,  only : put_orbmod_params
use slaarnacl, only : put_stowf_params
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaaab1)
logical,intent(in) :: ignore(xyzzyaaab1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7
logical xyzzyaaae7
bad_params=.false.
xyzzyaaac7=0
xyzzyaaad7=0
if(iparam_buffer>0)then
xyzzyaaac7=xyzzyaaad1(iparam_buffer)
xyzzyaaad7=iparam_buffer-sum(xyzzyaaac1(1:xyzzyaaac7-1))
endif
xyzzyaaab7=0
if(xyzzyaaac1(1)>0)then
xyzzyaaaa7=xyzzyaaab7+1
xyzzyaaab7=xyzzyaaab7+xyzzyaaac1(1)
if(xyzzyaaac7==1.or.xyzzyaaac7==0)then
call put_freeorb_params(params(xyzzyaaaa7:xyzzyaaab7),ignore(xyzzyaaaa&
&7:xyzzyaaab7),xyzzyaaad7,prestore,xyzzyaaae7)
bad_params=bad_params.or.xyzzyaaae7
endif
endif
if(xyzzyaaac1(2)>0)then
xyzzyaaaa7=xyzzyaaab7+1
xyzzyaaab7=xyzzyaaab7+xyzzyaaac1(2)
if(xyzzyaaac7==2.or.xyzzyaaac7==0)then
call put_orbmod_params(params(xyzzyaaaa7:xyzzyaaab7),ignore(xyzzyaaaa7&
&:xyzzyaaab7),xyzzyaaad7,prestore,xyzzyaaae7)
bad_params=bad_params.or.xyzzyaaae7
endif
endif
if(xyzzyaaac1(3)>0)then
xyzzyaaaa7=xyzzyaaab7+1
xyzzyaaab7=xyzzyaaab7+xyzzyaaac1(3)
if(xyzzyaaac7==3.or.xyzzyaaac7==0)then
call put_gwfmolorb_params(params(xyzzyaaaa7:xyzzyaaab7),ignore(xyzzyaa&
&aa7:xyzzyaaab7),xyzzyaaad7,prestore,xyzzyaaae7)
bad_params=bad_params.or.xyzzyaaae7
endif
endif
if(xyzzyaaac1(4)>0)then
xyzzyaaaa7=xyzzyaaab7+1
xyzzyaaab7=xyzzyaaab7+xyzzyaaac1(4)
if(xyzzyaaac7==4.or.xyzzyaaac7==0)then
call put_stowf_params(params(xyzzyaaaa7:xyzzyaaab7),ignore(xyzzyaaaa7:&
&xyzzyaaab7),xyzzyaaad7,prestore,xyzzyaaae7)
bad_params=bad_params.or.xyzzyaaae7
endif
endif
end subroutine put_wfdet_params
subroutine enumerate_plot_wfdet(nplot,keyword,description)
implicit none
integer,intent(inout) :: nplot
character(64),intent(inout),optional :: keyword(:),description(:)
nplot=0
end subroutine enumerate_plot_wfdet
subroutine query_plot_wfdet(iplot,ii,rank,is_complex,has_stderr,rot_te&
&nsor,transl_pos,nfunctions,function_name)
implicit none
integer,intent(in) :: iplot,ii
integer,intent(inout) :: rank,nfunctions
logical,intent(inout) :: is_complex,has_stderr,rot_tensor,transl_pos
character(64),intent(inout),optional :: function_name(:)
nfunctions=0
end subroutine query_plot_wfdet
subroutine get_plot_wfdet(iplot,ii,is0,is1,f)
implicit none
integer,intent(in) :: iplot,ii,is0,is1
real(dp),intent(out) :: f(*)
f(1)=0.d0
end subroutine get_plot_wfdet
subroutine finish_plot_wfdet
implicit none
end subroutine finish_plot_wfdet
subroutine xyzzyaabl1(norb,orbmap)
use store,           only : nemax,nspin,ndet,nele,heg_orbtype
use run_control,     only : errstop
use slaarnaac,          only : get_awfdet_orbmap
use slaarnaae,          only : get_bwfdet_orbmap
use slaarnaak,          only : get_dwfdet_orbmap
use slaarnaar,         only : get_exwfdet_orbmap
use slaarnaas,       only : get_free_orbmap
use slaarnaau,       only : get_gaussian_orbmap
use slaarnacb,          only : get_pwfdet_orbmap
use slaarnaci,             only : get_sdwfdet_orbmap
use slaarnack, only : get_special_wfn_orbmap
use slaarnacl,        only : get_stowfdet_orbmap
implicit none
integer,intent(inout) :: norb,orbmap(nemax,nspin,ndet)
integer xyzzyaaaa12(nspin),xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaa&
&ae12,xyzzyaaaf12,xyzzyaaag12
xyzzyaaaa12=0
norb=0
orbmap=0
select case(xyzzyaaae1)
case(0)
case(1)
call get_pwfdet_orbmap(xyzzyaaaa12,norb,orbmap)
case(2)
call get_gaussian_orbmap(xyzzyaaaa12,norb,orbmap)
case(3)
call get_awfdet_orbmap(xyzzyaaaa12,norb,orbmap)
case(4)
call get_bwfdet_orbmap(xyzzyaaaa12,norb,orbmap)
case(5)
call get_special_wfn_orbmap(xyzzyaaaa12,norb,orbmap)
case(6)
call get_stowfdet_orbmap(xyzzyaaaa12,norb,orbmap)
case(7)
call get_dwfdet_orbmap(xyzzyaaaa12,norb,orbmap)
case default
call errstop('GET_WFDET_ORBMAP','Dazed and confused.')
end select
xyzzyaaam1=norb
if(any(heg_orbtype/=0).and.any(heg_orbtype<=100.and.any(heg_orbtype/=3&
&)))call get_free_orbmap(xyzzyaaaa12,norb,orbmap)
if(any(heg_orbtype==3))call get_sdwfdet_orbmap(xyzzyaaaa12,norb,orbmap&
&)
if(any(heg_orbtype>100))call get_exwfdet_orbmap(xyzzyaaaa12,norb,orbma&
&p)
xyzzyaaan1=norb-xyzzyaaam1
do xyzzyaaab12=1,ndet
do xyzzyaaac12=1,nspin
do xyzzyaaad12=1,nele(xyzzyaaac12)
xyzzyaaaf12=orbmap(xyzzyaaad12,xyzzyaaac12,xyzzyaaab12)
if(xyzzyaaaf12==0)cycle
do xyzzyaaae12=xyzzyaaad12+1,nele(xyzzyaaac12)
xyzzyaaag12=orbmap(xyzzyaaae12,xyzzyaaac12,xyzzyaaab12)
if(xyzzyaaag12==0)cycle
if(xyzzyaaaf12==xyzzyaaag12)call errstop('GET_WFDET_ORBMAP','Orbital a&
&ppears more than once in a spin determinant.')
enddo
enddo
enddo
enddo
end subroutine xyzzyaabl1
subroutine get_wfdet_orbmask(norb,orbmap,orbmask)
use store,only : nemax,nspin,ndet,nele,detstart,detstop
implicit none
integer,intent(in) :: norb,orbmap(nemax,nspin,ndet)
logical,intent(out) :: orbmask(norb,nspin)
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13
orbmask=.false.
do xyzzyaaab13=1,nspin
do xyzzyaaaa13=detstart,detstop
do xyzzyaaac13=1,nele(xyzzyaaab13)
xyzzyaaad13=orbmap(xyzzyaaac13,xyzzyaaab13,xyzzyaaaa13)
if(xyzzyaaad13==0)cycle
orbmask(xyzzyaaad13,xyzzyaaab13)=.true.
enddo
enddo
enddo
end subroutine get_wfdet_orbmask
subroutine copy_orb_to_det(n,ispin,orbmap,orb_array,det_array)
use store, only : nemax,nspin,ndet,real1_complex2,detstart,detstop
implicit none
integer,intent(in) :: n,ispin,orbmap(nemax,nspin,ndet)
real(dp),intent(in) :: orb_array(n,wfdet_norb,real1_complex2)
real(dp),intent(inout) :: det_array(n,nemax,real1_complex2,ndet)
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14,xyzzyaaad14
do xyzzyaaad14=1,real1_complex2
do xyzzyaaaa14=detstart,detstop
do xyzzyaaab14=1,nemax
xyzzyaaac14=orbmap(xyzzyaaab14,ispin,xyzzyaaaa14)
if(xyzzyaaac14==0)then
det_array(1:n,xyzzyaaab14,xyzzyaaad14,xyzzyaaaa14)=0.d0
else
det_array(1:n,xyzzyaaab14,xyzzyaaad14,xyzzyaaaa14)=orb_array(1:n,xyzzy&
&aaac14,xyzzyaaad14)
endif
enddo
enddo
enddo
end subroutine copy_orb_to_det
subroutine copy_rmap_orb_to_det(ispin,orbmap,orb_m,orb_rmap,det_orb_m,&
&det_orb_rmap)
use store, only : nemax,nspin,ndet,detstart,detstop
implicit none
integer,intent(in) :: ispin,orbmap(nemax,nspin,ndet),orb_m,orb_rmap(wf&
&det_norb)
integer,intent(out) :: det_orb_m(ndet),det_orb_rmap(nemax,ndet)
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15
det_orb_m=0
det_orb_rmap=0
do xyzzyaaaa15=detstart,detstop
do xyzzyaaab15=1,nemax
xyzzyaaac15=orbmap(xyzzyaaab15,ispin,xyzzyaaaa15)
if(xyzzyaaac15==0)cycle
if(any(orb_rmap(1:orb_m)==xyzzyaaac15))then
det_orb_m(xyzzyaaaa15)=det_orb_m(xyzzyaaaa15)+1
det_orb_rmap(det_orb_m(xyzzyaaaa15),xyzzyaaaa15)=xyzzyaaab15
endif
enddo
enddo
end subroutine copy_rmap_orb_to_det
subroutine wfdet_write_mdet_casl()
use store
use casl
use format_utils, only : wout,i2s
use slaarnabp,         only : mdet_title,detcoef,detcoef_label
use parallel,     only : am_master
use run_control,  only : errwarn,errstop_master,timer
implicit none
character(512) errmsg
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16
if(xyzzyaabf1)then
call errwarn('WFDET_MDET_CASL','Already using a compressed expansion, &
&so not writing mdet.casl.')
return
endif
call timer('WFDET_WRITE_MDET_CASL',.true.)
call wout('Writing mdet.casl.')
call wout()
call set_casl_block('mdet.casl:MDET',errmsg,push=.true.)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
call set_casl_item('Title',trim(adjustl(mdet_title)),errmsg)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
call set_casl_block('Expansion',errmsg,push=.true.)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
do xyzzyaaab16=1,ndet
call set_casl_block('Term '//trim(i2s(xyzzyaaab16)),errmsg,push=.true.&
&)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
call set_casl_block('Coefficient',errmsg,prefer_inline=.true.)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
call set_casl_item('Coefficient:%u1',detcoef(xyzzyaaab16),errmsg)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
if(allocated(detcoef_label))then
call set_casl_item('Coefficient:Group',detcoef_label(xyzzyaaab16),errm&
&sg)
else
call set_casl_item('Coefficient:Group',1,errmsg)
endif
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
do xyzzyaaaa16=1,nspin
if(nele(xyzzyaaaa16)==0)cycle
call set_casl_block('Spin '//trim(i2s(xyzzyaaaa16)),errmsg,prefer_inli&
&ne=.true.,push=.true.)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
do xyzzyaaac16=1,nele(xyzzyaaaa16)
call set_casl_item('%u'//trim(i2s(xyzzyaaac16)),wfdet_orbmap(xyzzyaaac&
&16,xyzzyaaaa16,xyzzyaaab16),errmsg)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
enddo
call pop_casl_context()
enddo
call pop_casl_context()
enddo
call pop_casl_context()
call pop_casl_context()
if(am_master)then
call write_casl('mdet.casl','mdet.casl',errmsg)
if(len_trim(errmsg)>0)call errstop_master('WFDET_WRITE_MDET_CASL',trim&
&(errmsg))
endif
call timer('WFDET_WRITE_MDET_CASL',.false.)
end subroutine wfdet_write_mdet_casl
subroutine xyzzyaabm1
use casl
use slaarnaan,only : netot_nitot_products
use format_utils, only : i2s,wout
use slaarnabp,         only : detcoef,detcoef_label,orig_ndet,orig_det&
&coef
use parallel,     only : am_master
use run_control,  only : check_alloc,errstop_master
use store,        only : nele,nemax,nspin,ndet,detstart,detstop
implicit none
character(256) title
character(512) errmsg
integer ierr,xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17,xyzzyaaae&
&17,xyzzyaaaf17,xyzzyaaag17,xyzzyaaah17,xyzzyaaai17,xyzzyaaaj17,xyzzya&
&aak17
logical exists,is_block
real(dp) xyzzyaaal17
xyzzyaabf1=.false.
call read_casl('cmdet.casl',errmsg)
if(len_trim(errmsg)>0)call errstop_master('WFDET_READ_CMDET_CASL',trim&
&(errmsg))
call query_casl_item('cmdet.casl:CMDET',exists=exists)
if(.not.exists)return
xyzzyaabf1=.true.
call push_casl_context('cmdet.casl:CMDET')
call get_casl_item('Title',title,ierr=ierr)
if(ierr/=0)title='No title.'
if(am_master)then
call wout('Compressed determinant expansion')
call wout('================================')
call wout('Title: '//trim(adjustl(title)))
endif
call query_casl_item('Original expansion',exists=exists,is_block=is_bl&
&ock,nchildren=xyzzyaaaf17)
if(.not.exists.or..not.is_block)call errstop_master('WFDET_READ_CMDET_&
&CASL','Could not read "Original expansion" block.')
if(orig_ndet/=xyzzyaaaf17)call errstop_master('WFDET_READ_CMDET_CASL',&
&'File cmdet.casl contains the compression of a '//trim(i2s(orig_ndet)&
&)//'-determinant expansion, but the current expansion has '//trim(i2s&
&(xyzzyaaaf17))//' determinants.')
call push_casl_context('Original expansion')
do xyzzyaaab17=1,orig_ndet
call get_casl_item('Term '//trim(i2s(xyzzyaaab17))//':Coefficient:%u1'&
&,xyzzyaaal17,ierr=ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& coefficient #'//trim(i2s(xyzzyaaab17))//' in original expansion.')
call get_casl_item('Term '//trim(i2s(xyzzyaaab17))//':Coefficient:Grou&
&p',xyzzyaaaf17,ierr=ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& label of coefficient #'//trim(i2s(xyzzyaaab17))//' in original expan&
&sion.')
if(allocated(detcoef_label))then
if(xyzzyaaaf17/=detcoef_label(xyzzyaaab17))call errstop_master('WFDET_&
&READ_CMDET_CASL','File cmdet.casl is for an expansion whose coefficie&
&nt #'//trim(i2s(xyzzyaaab17))//' is labelled '//trim(i2s(xyzzyaaaf17)&
&)//', but in the current expansion its label is '//trim(i2s(detcoef_l&
&abel(xyzzyaaab17)))//'.')
endif
do xyzzyaaaa17=1,nspin
call query_casl_item('Term '//trim(i2s(xyzzyaaab17))//':Spin '//trim(i&
&2s(xyzzyaaaa17)),exists=exists,is_block=is_block,nchildren=xyzzyaaaf1&
&7)
if((.not.exists.or..not.is_block).and.nele(xyzzyaaaa17)>0)call errstop&
&_master('WFDET_READ_CMDET_CASL','Could not find orbital list for spin&
& channel '//trim(i2s(xyzzyaaaa17))//' for term '//trim(i2s(xyzzyaaab1&
&7))//' of the original expansion.   Perhaps this cmdet.casl file is f&
&or a different system?')
if(xyzzyaaaf17/=nele(xyzzyaaaa17))call errstop_master('WFDET_READ_CMDE&
&T_CASL','Wrong size of orbital list for spin channel '//trim(i2s(xyzz&
&yaaaa17))//' for term '//trim(i2s(xyzzyaaab17))//' of the original ex&
&pansion.  Perhaps this cmdet.casl file is for a different system?')
do xyzzyaaac17=1,nele(xyzzyaaaa17)
call get_casl_item('Term '//trim(i2s(xyzzyaaab17))//':Spin '//trim(i2s&
&(xyzzyaaaa17))//':%u'//trim(i2s(xyzzyaaac17)),xyzzyaaaf17,ierr=ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& orbital list for term #'//trim(i2s(xyzzyaaab17))//' in original expa&
&nsion.')
if(xyzzyaaaf17/=wfdet_orbmap(xyzzyaaac17,xyzzyaaaa17,xyzzyaaab17))call&
& errstop_master('WFDET_READ_CMDET_CASL','File cmdet.casl is for an ex&
&pansion whose term #'//trim(i2s(xyzzyaaab17))//' has orbital '//trim(&
&i2s(xyzzyaaaf17))//' in row '//trim(i2s(xyzzyaaac17))//' of the deter&
&minant for spin #'//trim(i2s(xyzzyaaaa17))//', but the in the current&
& expansion the orbital is '//trim(i2s(wfdet_orbmap(xyzzyaaac17,xyzzya&
&aaa17,xyzzyaaab17))))
enddo
enddo
enddo
call pop_casl_context()
xyzzyaaao1=wfdet_norb
call query_casl_item('Deduplicated coefficients',exists=exists,is_bloc&
&k=is_block,nchildren=xyzzyaaap1)
if(.not.exists.or..not.is_block)call errstop_master('WFDET_READ_CMDET_&
&CASL','Could not read "Deduplicated coefficients" block.')
call push_casl_context('Deduplicated coefficients')
allocate(xyzzyaaat1(xyzzyaaap1),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','dedup_ndetcoef')
xyzzyaaat1=0
do xyzzyaaab17=1,xyzzyaaap1
call query_casl_item('c_'//trim(i2s(xyzzyaaab17)),exists=exists,is_blo&
&ck=is_block,nchildren=xyzzyaaat1(xyzzyaaab17))
if(.not.exists.or..not.is_block)call errstop_master('WFDET_READ_CMDET_&
&CASL','Could not read block "c_'//trim(i2s(xyzzyaaab17))//'" in dedup&
&licated expansion.')
enddo
xyzzyaaah17=maxval(xyzzyaaat1)
allocate(xyzzyaaau1(xyzzyaaah17,xyzzyaaap1),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','dedup_idetcoef')
xyzzyaaau1=0
do xyzzyaaab17=1,xyzzyaaap1
do xyzzyaaad17=1,xyzzyaaat1(xyzzyaaab17)
call get_casl_item('c_'//trim(i2s(xyzzyaaab17))//':%u'//trim(i2s(xyzzy&
&aaad17)),xyzzyaaau1(xyzzyaaad17,xyzzyaaab17),ierr=ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& contribution #'//trim(i2s(xyzzyaaad17))//' to deduplicated coefficie&
&nt #'//trim(i2s(xyzzyaaab17))//'.')
enddo
enddo
call pop_casl_context()
call query_casl_item('Compressed expansion',exists=exists,is_block=is_&
&block)
if(.not.exists.or..not.is_block)call errstop_master('WFDET_READ_CMDET_&
&CASL','Could not read "Compressed expansion" block.')
call push_casl_context('Compressed expansion')
call query_casl_item('Orbital pool',exists=exists,is_block=is_block,nc&
&hildren=xyzzyaaaq1)
if(.not.exists.or..not.is_block)call errstop_master('WFDET_READ_CMDET_&
&CASL','Could not read "Orbital pool" block.')
call push_casl_context('Orbital pool')
allocate(xyzzyaaav1(xyzzyaaaq1),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','comp_nmix')
xyzzyaaav1=0
xyzzyaaaj17=0
do xyzzyaaac17=1,xyzzyaaaq1
call query_casl_item('Orbital '//trim(i2s(xyzzyaaac17)),exists=exists,&
&is_block=is_block,nchildren=xyzzyaaav1(xyzzyaaac17))
if(.not.exists.or..not.is_block)call errstop_master('WFDET_READ_CMDET_&
&CASL','Could not read block "Orbital '//trim(i2s(xyzzyaaac17))//'" in&
& compressed orbital pool.')
do xyzzyaaad17=1,xyzzyaaav1(xyzzyaaac17)
call query_casl_item('Orbital '//trim(i2s(xyzzyaaac17))//':Component '&
&//trim(i2s(xyzzyaaad17))//':Num',exists=exists,is_block=is_block,nchi&
&ldren=xyzzyaaaf17)
if(exists.and.is_block)xyzzyaaaj17=max(xyzzyaaaj17,xyzzyaaaf17)
enddo
enddo
xyzzyaaai17=maxval(xyzzyaaav1)
allocate(xyzzyaaaw1(xyzzyaaai17,xyzzyaaaq1),xyzzyaaax1(xyzzyaaai17,xyz&
&zyaaaq1),xyzzyaaay1(xyzzyaaaj17,xyzzyaaai17,xyzzyaaaq1),xyzzyaaaz1(xy&
&zzyaaaj17,xyzzyaaai17,xyzzyaaaq1),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','comp_imixcoeff')
xyzzyaaaw1=0
xyzzyaaax1=0
xyzzyaaay1=0
xyzzyaaaz1=0
do xyzzyaaac17=1,xyzzyaaaq1
call push_casl_context('Orbital '//trim(i2s(xyzzyaaac17)))
do xyzzyaaad17=1,xyzzyaaav1(xyzzyaaac17)
call push_casl_context('Component '//trim(i2s(xyzzyaaad17)))
call get_casl_item('%u1',xyzzyaaaw1(xyzzyaaad17,xyzzyaaac17),ierr=ierr&
&)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Cound not read&
& "Orbital" for component #'//trim(i2s(xyzzyaaad17))//' of compressed &
&orbital '//trim(i2s(xyzzyaaac17))//'.')
call query_casl_item('Num',exists=exists,is_block=is_block,nchildren=x&
&yzzyaaaf17)
if(exists.and.is_block)xyzzyaaax1(xyzzyaaad17,xyzzyaaac17)=xyzzyaaaf17
do xyzzyaaae17=1,xyzzyaaax1(xyzzyaaad17,xyzzyaaac17)
call get_casl_item('Num:%u'//trim(i2s(xyzzyaaae17)),xyzzyaaay1(xyzzyaa&
&ae17,xyzzyaaad17,xyzzyaaac17),ierr=ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& coefficient index #'//trim(i2s(xyzzyaaae17))//' in numerator of comp&
&onent #'//trim(i2s(xyzzyaaad17))//' of compressed orbital '//trim(i2s&
&(xyzzyaaac17))//'.')
call get_casl_item('Den:%u'//trim(i2s(xyzzyaaae17)),xyzzyaaaz1(xyzzyaa&
&ae17,xyzzyaaad17,xyzzyaaac17),ierr=ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& coefficient index #'//trim(i2s(xyzzyaaae17))//' in denominator of co&
&mponent #'//trim(i2s(xyzzyaaad17))//' of compressed orbital '//trim(i&
&2s(xyzzyaaac17))//'.')
enddo
call pop_casl_context()
enddo
call pop_casl_context()
enddo
call pop_casl_context()
call query_casl_item('Expansion',exists=exists,is_block=is_block,nchil&
&dren=xyzzyaaar1)
if(.not.exists.or..not.is_block)call errstop_master('WFDET_READ_CMDET_&
&CASL','Could not read "Expansion" block.')
call push_casl_context('Expansion')
allocate(xyzzyaaba1(xyzzyaaar1),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','comp_ndetcoef')
xyzzyaaba1=0
do xyzzyaaab17=1,xyzzyaaar1
call query_casl_item('Term '//trim(i2s(xyzzyaaab17))//':Coefficient:Nu&
&m',exists=exists,is_block=is_block,nchildren=xyzzyaaaf17)
if(exists.and.is_block)xyzzyaaba1(xyzzyaaab17)=xyzzyaaaf17
enddo
xyzzyaaak17=maxval(xyzzyaaba1)
allocate(xyzzyaabb1(xyzzyaaak17,xyzzyaaar1),xyzzyaabc1(xyzzyaaak17,xyz&
&zyaaar1),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','comp_idetcoef')
xyzzyaabb1=0
xyzzyaabc1=0
allocate(xyzzyaabd1(nemax,nspin,xyzzyaaar1),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','comp_orbmap')
xyzzyaabd1=0
do xyzzyaaab17=1,xyzzyaaar1
do xyzzyaaad17=1,xyzzyaaba1(xyzzyaaab17)
call get_casl_item('Term '//trim(i2s(xyzzyaaab17))//':Coefficient:Num:&
&%u'//trim(i2s(xyzzyaaad17)),xyzzyaabb1(xyzzyaaad17,xyzzyaaab17),ierr=&
&ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& coefficient index #'//trim(i2s(xyzzyaaad17))//' in numerator of pref&
&actor of compressed determinant #'//trim(i2s(xyzzyaaab17))//'.')
if(xyzzyaaad17<xyzzyaaba1(xyzzyaaab17))then
call get_casl_item('Term '//trim(i2s(xyzzyaaab17))//':Coefficient:Den:&
&%u'//trim(i2s(xyzzyaaad17)),xyzzyaabc1(xyzzyaaad17,xyzzyaaab17),ierr=&
&ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& coefficient index #'//trim(i2s(xyzzyaaad17))//' in denominator of pr&
&efactor of compressed determinant #'//trim(i2s(xyzzyaaab17))//'.')
endif
enddo
do xyzzyaaaa17=1,nspin
do xyzzyaaac17=1,nele(xyzzyaaaa17)
call get_casl_item('Term '//trim(i2s(xyzzyaaab17))//':Spin '//trim(i2s&
&(xyzzyaaaa17))//':%u'//trim(i2s(xyzzyaaac17)),xyzzyaabd1(xyzzyaaac17,&
&xyzzyaaaa17,xyzzyaaab17),ierr=ierr)
if(ierr/=0)call errstop_master('WFDET_READ_CMDET_CASL','Could not read&
& orbital index for row #'//trim(i2s(xyzzyaaac17))//' of the determina&
&nt for spin channel #'//trim(i2s(xyzzyaaaa17))//' of compressed term &
&#'//trim(i2s(xyzzyaaab17))//'.')
enddo
enddo
enddo
call pop_casl_context()
call pop_casl_context()
call pop_casl_context()
if(am_master)then
call wout('Compressed expansion loaded.')
call wout()
endif
allocate(xyzzyaaas1(nemax,nspin,orig_ndet),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','orig_detcoef, or&
&ig_orbmap')
orig_detcoef=detcoef
xyzzyaaas1=wfdet_orbmap
deallocate(detcoef,wfdet_orbmap)
allocate(detcoef(xyzzyaaar1),wfdet_orbmap(nemax,nspin,xyzzyaaar1),stat&
&=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','detcoef, wfdet_o&
&rbmap')
detcoef=0.d0
wfdet_orbmap=xyzzyaabd1
allocate(xyzzyaabh1(xyzzyaaap1),xyzzyaabe1(xyzzyaaai17,xyzzyaaaq1),xyz&
&zyaabi1(xyzzyaaai17,xyzzyaaaq1),xyzzyaabj1(xyzzyaaar1),xyzzyaabg1(xyz&
&zyaaao1,nspin),stat=xyzzyaaag17)
call check_alloc(xyzzyaaag17,'WFDET_READ_CMDET_CASL','dedup_detcoef, c&
&omp_mixcoeff, comp_detcoef')
xyzzyaabh1=0.d0
xyzzyaabi1=0.d0
xyzzyaabj1=0.d0
call get_wfdet_orbmask(xyzzyaaao1,xyzzyaaas1,xyzzyaabg1)
call wfdet_mdet_to_cmdet
wfdet_norb=xyzzyaaaq1
wfdet_orbmap=xyzzyaabd1
ndet=xyzzyaaar1
detstart=1
detstop=xyzzyaaar1
call netot_nitot_products
end subroutine xyzzyaabm1
subroutine wfdet_mdet_to_cmdet
use slaarnabp, only : detcoef,orig_detcoef
use slaarnabt, only : dcopy
implicit none
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xy&
&zzyaaaf18
real(dp) xyzzyaaag18,xyzzyaaah18,xyzzyaaai18
if(.not.xyzzyaabf1)then
detcoef=orig_detcoef
return
endif
do xyzzyaaaa18=1,xyzzyaaap1
xyzzyaaai18=0.d0
do xyzzyaaac18=1,xyzzyaaat1(xyzzyaaaa18)
xyzzyaaae18=xyzzyaaau1(xyzzyaaac18,xyzzyaaaa18)
xyzzyaaag18=0.d0
if(xyzzyaaae18>0)then
xyzzyaaag18=orig_detcoef(xyzzyaaae18)
elseif(xyzzyaaae18<0)then
xyzzyaaag18=-orig_detcoef(-xyzzyaaae18)
endif
xyzzyaaai18=xyzzyaaai18+xyzzyaaag18
enddo
xyzzyaabh1(xyzzyaaaa18)=xyzzyaaai18
enddo
do xyzzyaaab18=1,xyzzyaaaq1
do xyzzyaaac18=1,xyzzyaaav1(xyzzyaaab18)
xyzzyaabi1(xyzzyaaac18,xyzzyaaab18)=1.d0
xyzzyaaai18=1.d0
do xyzzyaaad18=1,xyzzyaaax1(xyzzyaaac18,xyzzyaaab18)
xyzzyaaae18=xyzzyaaay1(xyzzyaaad18,xyzzyaaac18,xyzzyaaab18)
xyzzyaaag18=1.d0
if(xyzzyaaae18>0)then
xyzzyaaag18=xyzzyaabh1(xyzzyaaae18)
elseif(xyzzyaaae18<0)then
xyzzyaaag18=-xyzzyaabh1(-xyzzyaaae18)
endif
xyzzyaaaf18=xyzzyaaaz1(xyzzyaaad18,xyzzyaaac18,xyzzyaaab18)
xyzzyaaah18=1.d0
if(xyzzyaaaf18>0)then
xyzzyaaah18=xyzzyaabh1(xyzzyaaaf18)
elseif(xyzzyaaaf18<0)then
xyzzyaaah18=-xyzzyaabh1(-xyzzyaaaf18)
endif
if(xyzzyaaah18/=0.d0)then
xyzzyaaai18=xyzzyaaai18*(xyzzyaaag18/xyzzyaaah18)
else
xyzzyaaai18=0.d0
endif
enddo
if(xyzzyaaaw1(xyzzyaaac18,xyzzyaaab18)<0)then
xyzzyaabi1(xyzzyaaac18,xyzzyaaab18)=-xyzzyaaai18
xyzzyaabe1(xyzzyaaac18,xyzzyaaab18)=-xyzzyaaaw1(xyzzyaaac18,xyzzyaaab1&
&8)
else
xyzzyaabi1(xyzzyaaac18,xyzzyaaab18)=xyzzyaaai18
xyzzyaabe1(xyzzyaaac18,xyzzyaaab18)=xyzzyaaaw1(xyzzyaaac18,xyzzyaaab18&
&)
endif
enddo
enddo
do xyzzyaaaa18=1,xyzzyaaar1
xyzzyaaai18=1.d0
do xyzzyaaac18=1,xyzzyaaba1(xyzzyaaaa18)
xyzzyaaae18=xyzzyaabb1(xyzzyaaac18,xyzzyaaaa18)
xyzzyaaag18=1.d0
if(xyzzyaaae18>0)then
xyzzyaaag18=xyzzyaabh1(xyzzyaaae18)
elseif(xyzzyaaae18<0)then
xyzzyaaag18=-xyzzyaabh1(-xyzzyaaae18)
endif
xyzzyaaaf18=xyzzyaabc1(xyzzyaaac18,xyzzyaaaa18)
xyzzyaaah18=1.d0
if(xyzzyaaaf18>0)then
xyzzyaaah18=xyzzyaabh1(xyzzyaaaf18)
elseif(xyzzyaaaf18<0)then
xyzzyaaah18=-xyzzyaabh1(-xyzzyaaaf18)
endif
if(xyzzyaaah18/=0.d0)then
xyzzyaaai18=xyzzyaaai18*(xyzzyaaag18/xyzzyaaah18)
else
xyzzyaaai18=0.d0
endif
enddo
xyzzyaabj1(xyzzyaaaa18)=xyzzyaaai18
enddo
call dcopy(xyzzyaaar1,xyzzyaabj1(1),1,detcoef(1),1)
end subroutine wfdet_mdet_to_cmdet
logical function wfdet_detcoef_affect_orbs()
implicit none
wfdet_detcoef_affect_orbs=xyzzyaabf1
end function wfdet_detcoef_affect_orbs
subroutine wfdet(rvec,spin,jspin,norb,orbmask,val,fsd,orbval,orbgrad,o&
&rblap,eevecs,orbsderivs,orb_m,orb_rmap)
use store,       only : real1_complex2,netot
use run_control, only : timer
implicit none
integer,intent(in) :: spin,jspin,norb
logical,intent(in) :: val,fsd,orbmask(norb)
real(dp),intent(in) :: rvec(3),eevecs(4,netot)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),optional,intent(inout) :: orbsderivs(6,norb,real1_complex2)
integer,optional,intent(inout) :: orb_m,orb_rmap(xyzzyaaam1)
integer xyzzyaaaa20,xyzzyaaab20,xyzzyaaac20,xyzzyaaad20
real(dp) xyzzyaaae20,xyzzyaaaf20,xyzzyaaag20(3),xyzzyaaah20,xyzzyaaai2&
&0(6)
real(dp) xyzzyaaaj20(xyzzyaaao1,real1_complex2),xyzzyaaak20(3,xyzzyaaa&
&o1,real1_complex2),xyzzyaaal20(xyzzyaaao1,real1_complex2),xyzzyaaam20&
&(6,xyzzyaaao1,real1_complex2)
if(.not.xyzzyaabf1)then
call xyzzyaabk1(rvec,spin,jspin,norb,orbmask,val,fsd,orbval,orbgrad,or&
&blap,eevecs,orbsderivs=orbsderivs,orb_m=orb_m,orb_rmap=orb_rmap)
return
endif
call timer('WFDET_CDEXP_WRAP',.true.)
if(present(orbsderivs).and.fsd)then
if(val)xyzzyaaaj20=0.d0
xyzzyaaak20=0.d0
xyzzyaaam20=0.d0
call xyzzyaabk1(rvec,spin,jspin,xyzzyaaao1,xyzzyaabg1(1,jspin),val,fsd&
&,xyzzyaaaj20,xyzzyaaak20,xyzzyaaal20,eevecs,xyzzyaaam20)
do xyzzyaaaa20=1,real1_complex2
do xyzzyaaab20=1,xyzzyaaaq1
xyzzyaaaf20=0.d0
xyzzyaaag20(1:3)=0.d0
xyzzyaaai20(1:6)=0.d0
do xyzzyaaad20=1,xyzzyaaav1(xyzzyaaab20)
xyzzyaaac20=xyzzyaabe1(xyzzyaaad20,xyzzyaaab20)
xyzzyaaae20=xyzzyaabi1(xyzzyaaad20,xyzzyaaab20)
if(val)xyzzyaaaf20=xyzzyaaaf20+xyzzyaaaj20(xyzzyaaac20,xyzzyaaaa20)*xy&
&zzyaaae20
xyzzyaaag20(1:3)=xyzzyaaag20(1:3)+xyzzyaaak20(1:3,xyzzyaaac20,xyzzyaaa&
&a20)*xyzzyaaae20
xyzzyaaai20(1:6)=xyzzyaaai20(1:6)+xyzzyaaam20(1:6,xyzzyaaac20,xyzzyaaa&
&a20)*xyzzyaaae20
enddo
if(val)orbval(xyzzyaaab20,xyzzyaaaa20)=xyzzyaaaf20
orbgrad(1:3,xyzzyaaab20,xyzzyaaaa20)=xyzzyaaag20(1:3)
orbsderivs(1:6,xyzzyaaab20,xyzzyaaaa20)=xyzzyaaai20(1:6)
enddo
enddo
else
if(val)xyzzyaaaj20=0.d0
if(fsd)then
xyzzyaaak20=0.d0
xyzzyaaal20=0.d0
endif
call xyzzyaabk1(rvec,spin,jspin,xyzzyaaao1,xyzzyaabg1(1,jspin),val,fsd&
&,xyzzyaaaj20,xyzzyaaak20,xyzzyaaal20,eevecs)
if(fsd)then
do xyzzyaaaa20=1,real1_complex2
do xyzzyaaab20=1,xyzzyaaaq1
xyzzyaaaf20=0.d0
xyzzyaaag20(1:3)=0.d0
xyzzyaaah20=0.d0
do xyzzyaaad20=1,xyzzyaaav1(xyzzyaaab20)
xyzzyaaac20=xyzzyaabe1(xyzzyaaad20,xyzzyaaab20)
xyzzyaaae20=xyzzyaabi1(xyzzyaaad20,xyzzyaaab20)
if(val)xyzzyaaaf20=xyzzyaaaf20+xyzzyaaaj20(xyzzyaaac20,xyzzyaaaa20)*xy&
&zzyaaae20
if(fsd)then
xyzzyaaag20(1:3)=xyzzyaaag20(1:3)+xyzzyaaak20(1:3,xyzzyaaac20,xyzzyaaa&
&a20)*xyzzyaaae20
xyzzyaaah20=xyzzyaaah20+xyzzyaaal20(xyzzyaaac20,xyzzyaaaa20)*xyzzyaaae&
&20
endif
enddo
if(val)orbval(xyzzyaaab20,xyzzyaaaa20)=xyzzyaaaf20
if(fsd)then
orbgrad(1:3,xyzzyaaab20,xyzzyaaaa20)=xyzzyaaag20(1:3)
orblap(xyzzyaaab20,xyzzyaaaa20)=xyzzyaaah20
endif
enddo
enddo
elseif(val)then
do xyzzyaaaa20=1,real1_complex2
do xyzzyaaab20=1,xyzzyaaaq1
xyzzyaaaf20=0.d0
do xyzzyaaad20=1,xyzzyaaav1(xyzzyaaab20)
xyzzyaaaf20=xyzzyaaaf20+xyzzyaaaj20(xyzzyaabe1(xyzzyaaad20,xyzzyaaab20&
&),xyzzyaaaa20)*xyzzyaabi1(xyzzyaaad20,xyzzyaaab20)
enddo
orbval(xyzzyaaab20,xyzzyaaaa20)=xyzzyaaaf20
enddo
enddo
endif
endif
call timer('WFDET_CDEXP_WRAP',.false.)
end subroutine wfdet
real(dp) function get_wfdet_rmax()
use slaarnaac,     only : get_awfdet_rmax
use slaarnaae,     only : get_bwfdet_rmax
use slaarnaau,  only : get_gaussian_rmax
implicit none
select case(xyzzyaaae1)
case(2)
get_wfdet_rmax=get_gaussian_rmax()
case(3)
get_wfdet_rmax=get_awfdet_rmax()
case(4)
get_wfdet_rmax=get_bwfdet_rmax()
end select
end function get_wfdet_rmax
end module slaarnacq

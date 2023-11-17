module slaarnaab
use dsp,         only : dp
use slaarnaam,only : eval_local_energy,i_kei,i_potinl,energy_scratch_r&
&equest
use format_utils,only : wout,i2s,wordwrap
use slaarnaat,   only : valence_charge
use slaarnabg,    only : iontype,rion,nitot,nitot_forces,naxis_forces,&
&nitot_forces,isperiodic,dimensionality
use slaarnabs,   only : nl_rcut,nl_nrefgrid,nl_wrkgrid_all,nl_refgridw&
&,nl_wfrat_forces,dvelf_forces,v_non_local_tmove
use slaarnabt,   only : interp_nev_with_derivs,lookup,ddot,make_repres&
&entable
use parallel,    only : am_master
use slaarnaca,       only : rcut_loc,rcut_non_loc,pp_radial_grid,fcoef&
&f_loc_full,fcoeff_non_loc_full,ncoeff,have_ae,nlang,have_ppots
use run_control, only : errstop,errwarn,errstop_master,check_alloc
use slaarnach,     only : nscratch,rele_scr,eivecs_scr,get_rsele,which&
&_scratch,scratch_request,get_eivecs
use store,       only : netot,use_tmove,isvmc,isdmc,isrmc_rmc,isrmc,is&
&dmc_dmc,isvmc_dmc
implicit none
private
public forces_scratch_request,setup_forces_accum,finish_forces_accum,e&
&val_local_forces,forces_to_future,future_to_forces,tagh_forces,nfterm&
&s,nfcomps,forces_info,initialize_forces
real(dp),allocatable :: xyzzyaaaa1(:,:),xyzzyaaab1(:,:),xyzzyaaac1(:,:&
&),xyzzyaaad1(:,:),xyzzyaaae1(:,:),xyzzyaaaf1(:,:),xyzzyaaag1(:,:),xyz&
&zyaaah1(:,:)
real(dp),allocatable :: xyzzyaaai1(:,:,:)
real(dp),allocatable :: xyzzyaaaj1(:,:)
real(dp),allocatable :: xyzzyaaak1(:,:),xyzzyaaal1(:,:),xyzzyaaam1(:,:&
&)
real(dp),allocatable :: xyzzyaaan1(:,:)
real(dp),parameter :: xyzzyaaao1=1.d-7
real(dp),allocatable :: xyzzyaaap1(:,:)
integer forces_info
integer nfcomps,nfterms
integer,allocatable :: tagh_forces(:,:,:)
integer,parameter :: xyzzyaaaq1=50
integer :: xyzzyaaar1(xyzzyaaaq1)=0,xyzzyaaas1(xyzzyaaaq1)=0,xyzzyaaat&
&1=0
integer,allocatable :: xyzzyaaau1(:)
contains
subroutine initialize_forces
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2
real(dp) xyzzyaaae2(3),xyzzyaaaf2,xyzzyaaag2(3),xyzzyaaah2(3),xyzzyaaa&
&i2,xyzzyaaaj2(3),xyzzyaaak2(3),xyzzyaaal2,xyzzyaaam2(3),xyzzyaaan2(3)&
&,xyzzyaaao2,xyzzyaaap2
real(dp),parameter :: xyzzyaaaq2=1.d-8
logical xyzzyaaar2,xyzzyaaas2,xyzzyaaat2
if(nitot==0)call errstop_master('INITIALIZE_FORCES','No atoms in this &
&system.')
if(have_ae)call errstop_master('INITIALIZE_FORCES','Calculation of ato&
&mic forces cannot be done for all-electron atoms. Pseudopotentials mu&
&st be used.')
if(any(nlang/=2))call errstop_master('INITIALIZE_FORCES','At present t&
&he force module assumes that only s, p and d angular-momentum compone&
&nts are present. Please ask someone to generalize this.')
if(isperiodic)then
naxis_forces=dimensionality
xyzzyaaar2=.false.
xyzzyaaas2=.false.
xyzzyaaat2=.false.
else
xyzzyaaar2=.true.
xyzzyaaas2=.true.
xyzzyaaat2=.true.
do xyzzyaaaa2=2,nitot
xyzzyaaae2=rion(:,1)-rion(:,xyzzyaaaa2)
xyzzyaaaf2=sqrt(sum(xyzzyaaae2(:)**2))
if(xyzzyaaaf2>xyzzyaaaq2)then
xyzzyaaag2=xyzzyaaae2/xyzzyaaaf2
xyzzyaaar2=.false.
exit
endif
enddo
if(xyzzyaaar2)call errstop_master('INITIALIZE_FORCES','Forces on a poi&
&nt-like system are zero.')
do xyzzyaaab2=xyzzyaaaa2+1,nitot
xyzzyaaah2=rion(:,1)-rion(:,xyzzyaaab2)
xyzzyaaai2=sqrt(sum(xyzzyaaah2(:)**2))
if(xyzzyaaai2<xyzzyaaaq2)continue
xyzzyaaaj2=xyzzyaaah2/xyzzyaaai2
xyzzyaaan2(1)=xyzzyaaag2(2)*xyzzyaaaj2(3)-xyzzyaaag2(3)*xyzzyaaaj2(2)
xyzzyaaan2(2)=xyzzyaaag2(3)*xyzzyaaaj2(1)-xyzzyaaag2(1)*xyzzyaaaj2(3)
xyzzyaaan2(3)=xyzzyaaag2(1)*xyzzyaaaj2(2)-xyzzyaaag2(2)*xyzzyaaaj2(1)
xyzzyaaao2=sqrt(sum(xyzzyaaan2(:)**2))
if(xyzzyaaao2>xyzzyaaaq2)then
xyzzyaaas2=.false.
exit
endif
enddo
if(.not.xyzzyaaas2)then
do xyzzyaaac2=xyzzyaaab2+1,nitot
xyzzyaaak2=rion(:,1)-rion(:,xyzzyaaac2)
xyzzyaaal2=sqrt(sum(xyzzyaaak2(:)**2))
if(xyzzyaaal2<xyzzyaaaq2)continue
xyzzyaaam2=xyzzyaaak2/xyzzyaaal2
xyzzyaaap2=sum(xyzzyaaan2(:)*xyzzyaaam2(:))
if(xyzzyaaap2>xyzzyaaaq2)then
xyzzyaaat2=.false.
exit
endif
enddo
endif
if(xyzzyaaas2)then
naxis_forces=1
if(abs(xyzzyaaag2(2))>xyzzyaaaq2)naxis_forces=2
if(abs(xyzzyaaag2(3))>xyzzyaaaq2)naxis_forces=3
elseif(xyzzyaaat2)then
naxis_forces=2
if(abs(xyzzyaaag2(3))>xyzzyaaaq2.or.abs(xyzzyaaaj2(3))>xyzzyaaaq2)naxi&
&s_forces=3
else
naxis_forces=3
endif
endif
if(am_master)then
call wout('Forces')
call wout('======')
if(xyzzyaaas2)then
call wout('Atoms in system are arranged on a line.')
elseif(xyzzyaaat2)then
call wout('Atoms in system are arranged on a plane.')
elseif(.not.isperiodic)then
call wout('Atoms in system are in a three-dimensional arrangement.')
else
call wout('System is periodic.')
endif
select case(naxis_forces)
case(1)
call wout('Forces will be evaluated along the x axis.')
case(2)
call wout('Forces will be evaluated along the x and y axes.')
case(3)
call wout('Forces will be evaluated along the x, y and z axes.')
case default
call errstop('INITIALIZE_FORCES','Problem in detecting number of axes &
&for forces evaluation.')
end select
if(xyzzyaaas2.and.naxis_forces>1)then
call errwarn('INITIALIZE_FORCES','The current system is linear, but it&
& is not aligned parallel to the x-axis. It is recommended that you ro&
&tate it in order to reduce the cost of evaluating forces.')
elseif(xyzzyaaat2.and.naxis_forces>2)then
call errwarn('INITIALIZE_FORCES','The current system is planar, but it&
& is not oriented parallel to the XY-plane. It is recommended that you&
& rotate it in order to reduce the cost of evaluating forces.')
endif
call wout()
endif
nitot_forces=nitot
call xyzzyaaav1
if(isvmc)then
nfterms=11
nfcomps=nfterms*nitot_forces*naxis_forces
elseif(isrmc.or.isrmc_rmc.or.isdmc.or.isvmc_dmc.or.isdmc_dmc)then
nfterms=22
nfcomps=nfterms*nitot_forces*naxis_forces
else
call errstop('INITIALIZE_FORCES','Problem assigning items for forces.'&
&)
endif
if(am_master.and.(forces_info>2))then
call wout('Number of forces items evaluated: '//trim(i2s(nfterms)))
call wout()
endif
allocate(tagh_forces(nfterms,3,nitot_forces),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'INITIALIZE_FORCES','tagh_forces')
end subroutine initialize_forces
subroutine forces_scratch_request(is0)
implicit none
integer,intent(inout) :: is0
integer xyzzyaaaa3
if(any(xyzzyaaar1(1:xyzzyaaat1)==is0).and.is0/=0)return
xyzzyaaat1=xyzzyaaat1+1
call energy_scratch_request(is0,nl_with_drift=.true.)
xyzzyaaar1(xyzzyaaat1)=is0
xyzzyaaaa3=0
call scratch_request(ratio_ion_from=is0,ratio_ion_to=xyzzyaaaa3)
call energy_scratch_request(xyzzyaaaa3)
xyzzyaaas1(xyzzyaaat1)=xyzzyaaaa3
end subroutine forces_scratch_request
subroutine setup_forces_accum
implicit none
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4
allocate(xyzzyaaau1(nscratch),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'SETUP_FORCES','basic0')
xyzzyaaau1=0
do xyzzyaaac4=1,xyzzyaaat1
xyzzyaaaa4=xyzzyaaar1(xyzzyaaac4)
xyzzyaaab4=xyzzyaaas1(xyzzyaaac4)
call which_scratch(xyzzyaaaa4)
call which_scratch(xyzzyaaab4)
xyzzyaaau1(xyzzyaaaa4)=xyzzyaaab4
enddo
allocate(xyzzyaaai1(3,3,nitot_forces),xyzzyaaak1(naxis_forces,nitot_fo&
&rces),xyzzyaaam1(naxis_forces,nitot_forces),xyzzyaaal1(naxis_forces,n&
&itot_forces),xyzzyaaan1(naxis_forces,nitot_forces),xyzzyaaap1(naxis_f&
&orces,nitot_forces),xyzzyaaaa1(3,nitot_forces),xyzzyaaab1(3,nitot_for&
&ces),xyzzyaaac1(3,nitot_forces),xyzzyaaad1(3,nitot_forces),xyzzyaaae1&
&(3,nitot_forces),xyzzyaaaf1(3,nitot_forces),xyzzyaaag1(3,nitot_forces&
&),xyzzyaaah1(3,nitot_forces),stat=xyzzyaaad4)
call check_alloc(xyzzyaaad4,'SETUP_FORCES','basic')
end subroutine setup_forces_accum
subroutine finish_forces_accum
implicit none
deallocate(xyzzyaaai1,xyzzyaaak1,xyzzyaaam1,xyzzyaaal1,xyzzyaaan1,xyzz&
&yaaap1,xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae1,xyzzyaa&
&af1,xyzzyaaag1,xyzzyaaah1)
deallocate(xyzzyaaau1)
xyzzyaaar1=0
xyzzyaaas1=0
xyzzyaaat1=0
end subroutine finish_forces_accum
subroutine xyzzyaaav1
implicit none
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6
real(dp) xyzzyaaad6(3),xyzzyaaae6
allocate(xyzzyaaaj1(3,nitot_forces),stat=xyzzyaaac6)
call check_alloc(xyzzyaaac6,'FORCES','neion')
do xyzzyaaaa6=1,nitot_forces
xyzzyaaad6(1:3)=0.d0
do xyzzyaaab6=1,nitot
if(xyzzyaaab6/=xyzzyaaaa6)then
xyzzyaaae6=sqrt(sum((rion(:,xyzzyaaaa6)-rion(:,xyzzyaaab6))**2))**3
xyzzyaaad6(1:3)=xyzzyaaad6(1:3)+valence_charge(xyzzyaaaa6)*valence_cha&
&rge(xyzzyaaab6)*(rion(1:3,xyzzyaaaa6)-rion(1:3,xyzzyaaab6))/xyzzyaaae&
&6
endif
enddo
xyzzyaaaj1(1:3,xyzzyaaaa6)=xyzzyaaad6(1:3)/dble(netot)
enddo
if(am_master.and.(forces_info>2))then
call wout()
do xyzzyaaaa6=1,nitot
call wout('Atom '//trim(i2s(xyzzyaaaa6))//' has valence charge ',valen&
&ce_charge(xyzzyaaaa6),rfmt='(f8.4)',adjust=.true.)
enddo
call wout()
call wordwrap('Nuclear position in a.u. for x,y,z axes for each atom i&
&n the order that atoms appear in gwfn.data file:')
do xyzzyaaaa6=1,nitot_forces
call wout('',rion(1:3,xyzzyaaaa6),rfmt='(f20.12)')
enddo
call wout()
call wordwrap('Ion-Ion forces in a.u. along x,y,z axes evaluated at ea&
&ch atom, ordered as they appear in the gwfn.data file:')
do xyzzyaaaa6=1,nitot_forces
call wout('',netot*xyzzyaaaj1(1:3,xyzzyaaaa6),rfmt='(f20.12)')
enddo
call wout()
call wout('Number of electrons: '//trim(i2s(netot)))
call wout('Number of atoms on which forces are evaluated: ' //trim(i2s&
&(nitot_forces)))
call wout('Total number of atoms in system: '//trim(i2s(nitot)))
call wout()
endif
end subroutine xyzzyaaav1
subroutine xyzzyaaaw1(ii,rele1,eivecs1)
implicit none
integer,intent(in) :: ii
real(dp),intent(in) :: rele1(3),eivecs1(4,nitot)
integer xyzzyaaaa7,xyzzyaaab7
real(dp) xyzzyaaac7(3),xyzzyaaad7(3),xyzzyaaae7
do xyzzyaaaa7=1,nitot_forces
xyzzyaaae7=eivecs1(4,xyzzyaaaa7)
if(xyzzyaaae7<=nl_rcut(xyzzyaaaa7))then
xyzzyaaad7=eivecs1(1:3,xyzzyaaaa7)
xyzzyaaac7=rele1-xyzzyaaad7
xyzzyaaab7=iontype(xyzzyaaaa7)
call xyzzyaaax1(xyzzyaaaa7,xyzzyaaad7,xyzzyaaae7,xyzzyaaac7,nl_nrefgri&
&d(xyzzyaaab7),nl_wrkgrid_all(1,1,xyzzyaaaa7,ii),nl_refgridw(1,xyzzyaa&
&ab7),nl_wfrat_forces(1,xyzzyaaaa7,ii),dvelf_forces(1,1,xyzzyaaaa7,ii)&
&)
endif
enddo
end subroutine xyzzyaaaw1
subroutine xyzzyaaax1(ion,rvec,r,ionvec,ngrid,grid,gridw,wfrat,dvelf)
implicit none
integer,intent(in) :: ngrid,ion
real(dp),intent(in) :: rvec(3),r,ionvec(3),gridw(*),wfrat(*),dvelf(3,*&
&)
real(dp),intent(inout) :: grid(3,*)
integer xyzzyaaaa8,xyzzyaaab8
real(dp) xyzzyaaac8,xyzzyaaad8(3),xyzzyaaae8,xyzzyaaaf8,xyzzyaaag8,xyz&
&zyaaah8,xyzzyaaai8,xyzzyaaaj8(3),xyzzyaaak8(3)
xyzzyaaae8=1.d0/r
xyzzyaaak8=rvec*xyzzyaaae8
xyzzyaaaa1(1:3,ion)=0.d0
xyzzyaaab1(1:3,ion)=0.d0
xyzzyaaac1(1:3,ion)=0.d0
xyzzyaaad1(1:3,ion)=0.d0
xyzzyaaae1(1:3,ion)=0.d0
xyzzyaaaf1(1:3,ion)=0.d0
xyzzyaaag1(1:3,ion)=0.d0
xyzzyaaah1(1:3,ion)=0.d0
do xyzzyaaaa8=1,ngrid
xyzzyaaaj8(1:3)=xyzzyaaae8*(grid(1:3,xyzzyaaaa8)-ionvec(1:3))
xyzzyaaac8=xyzzyaaae8*ddot(3,rvec(1),1,xyzzyaaaj8(1),1)
xyzzyaaaf8=wfrat(xyzzyaaaa8)*gridw(xyzzyaaaa8)
xyzzyaaag8=ddot(3,dvelf(1,xyzzyaaaa8),1,xyzzyaaaj8(1),1)
do xyzzyaaab8=1,naxis_forces
xyzzyaaad8(xyzzyaaab8)=gridw(xyzzyaaaa8)*(dvelf(xyzzyaaab8,xyzzyaaaa8)&
&-xyzzyaaak8(xyzzyaaab8)*xyzzyaaag8)
enddo
if(forces_info>2)then
xyzzyaaag8=3.d0*xyzzyaaac8*xyzzyaaae8*xyzzyaaaf8
xyzzyaaah8=1.5d0*xyzzyaaac8*xyzzyaaac8-0.5d0
xyzzyaaai8=xyzzyaaah8*xyzzyaaaf8
do xyzzyaaab8=1,naxis_forces
xyzzyaaaf1(xyzzyaaab8,ion)=xyzzyaaaf1(xyzzyaaab8,ion)+xyzzyaaag8*(xyzz&
&yaaac8*xyzzyaaak8(xyzzyaaab8)-xyzzyaaaj8(xyzzyaaab8))
xyzzyaaag1(xyzzyaaab8,ion)=xyzzyaaag1(xyzzyaaab8,ion)-xyzzyaaak8(xyzzy&
&aaab8)*xyzzyaaai8
xyzzyaaah1(xyzzyaaab8,ion)=xyzzyaaah1(xyzzyaaab8,ion)+xyzzyaaah8*xyzzy&
&aaad8(xyzzyaaab8)
enddo
endif
xyzzyaaag8=xyzzyaaaf8*xyzzyaaae8
xyzzyaaah8=xyzzyaaaf8*xyzzyaaac8
do xyzzyaaab8=1,naxis_forces
xyzzyaaaa1(xyzzyaaab8,ion)=xyzzyaaaa1(xyzzyaaab8,ion)+xyzzyaaad8(xyzzy&
&aaab8)
xyzzyaaab1(xyzzyaaab8,ion)=xyzzyaaab1(xyzzyaaab8,ion)-xyzzyaaak8(xyzzy&
&aaab8)*xyzzyaaaf8
xyzzyaaac1(xyzzyaaab8,ion)=xyzzyaaac1(xyzzyaaab8,ion)+(xyzzyaaac8*xyzz&
&yaaak8(xyzzyaaab8)-xyzzyaaaj8(xyzzyaaab8))*xyzzyaaag8
xyzzyaaad1(xyzzyaaab8,ion)=xyzzyaaad1(xyzzyaaab8,ion)-xyzzyaaak8(xyzzy&
&aaab8)*xyzzyaaah8
xyzzyaaae1(xyzzyaaab8,ion)=xyzzyaaae1(xyzzyaaab8,ion)+xyzzyaaac8*xyzzy&
&aaad8(xyzzyaaab8)
enddo
enddo
end subroutine xyzzyaaax1
subroutine xyzzyaaay1(eivecs1,xyzzyaaai1)
implicit none
real(dp),intent(in) :: eivecs1(4,nitot)
real(dp),intent(inout) :: xyzzyaaai1(3,3,nitot_forces)
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9
integer,parameter :: xyzzyaaae9=5
real(dp) xyzzyaaaf9,xyzzyaaag9,xyzzyaaah9(3),xyzzyaaai9(3),xyzzyaaaj9(&
&3),xyzzyaaak9(3),xyzzyaaal9(3),xyzzyaaam9,xyzzyaaan9,xyzzyaaao9,xyzzy&
&aaap9,xyzzyaaaq9,xyzzyaaar9,xyzzyaaas9(3),xyzzyaaat9(3),xyzzyaaau9(3)
do xyzzyaaaa9=1,nitot_forces
xyzzyaaac9=iontype(xyzzyaaaa9)
xyzzyaaaf9=eivecs1(4,xyzzyaaaa9)
xyzzyaaag9=1.d0/xyzzyaaaf9
xyzzyaaah9=eivecs1(1:3,xyzzyaaaa9)
xyzzyaaai9=xyzzyaaah9*xyzzyaaag9
xyzzyaaaj9=xyzzyaaai9*xyzzyaaag9*xyzzyaaag9
call lookup(pp_radial_grid(1,xyzzyaaac9),ncoeff(xyzzyaaac9),xyzzyaaaf9&
&,xyzzyaaab9)
xyzzyaaab9=min(max(xyzzyaaab9-(xyzzyaaae9-1)/2,1),ncoeff(xyzzyaaac9)+1&
&-xyzzyaaae9)
if(forces_info>2)then
xyzzyaaak9=0.d0
xyzzyaaal9=0.d0
if(xyzzyaaaf9<rcut_loc(xyzzyaaac9))then
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_loc_full(xyzzyaaab9,0,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzzyaaam9&
&,xyzzyaaan9,xyzzyaaao9)
do xyzzyaaad9=1,naxis_forces
xyzzyaaak9(xyzzyaaad9)=xyzzyaaan9*xyzzyaaai9(xyzzyaaad9)
enddo
else
do xyzzyaaad9=1,naxis_forces
xyzzyaaak9(xyzzyaaad9)=valence_charge(xyzzyaaaa9)*xyzzyaaaj9(xyzzyaaad&
&9)
enddo
endif
if(xyzzyaaaf9<rcut_non_loc(xyzzyaaac9))then
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_non_loc_full(xyzzyaaab9,1,0,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzz&
&yaaam9,xyzzyaaan9,xyzzyaaao9)
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_non_loc_full(xyzzyaaab9,2,0,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzz&
&yaaap9,xyzzyaaaq9,xyzzyaaar9)
do xyzzyaaad9=1,naxis_forces
xyzzyaaal9(xyzzyaaad9)=xyzzyaaal9(xyzzyaaad9)-3.d0*(xyzzyaaam9*(xyzzya&
&aac1(xyzzyaaad9,xyzzyaaaa9)+xyzzyaaae1(xyzzyaaad9,xyzzyaaaa9))+xyzzya&
&aad1(xyzzyaaad9,xyzzyaaaa9)*xyzzyaaan9)
xyzzyaaal9(xyzzyaaad9)=xyzzyaaal9(xyzzyaaad9)-5.d0*(xyzzyaaap9*(xyzzya&
&aaf1(xyzzyaaad9,xyzzyaaaa9)+xyzzyaaah1(xyzzyaaad9,xyzzyaaaa9))+xyzzya&
&aag1(xyzzyaaad9,xyzzyaaaa9)*xyzzyaaaq9)
enddo
endif
do xyzzyaaad9=1,naxis_forces
xyzzyaaas9(xyzzyaaad9)=xyzzyaaak9(xyzzyaaad9)+xyzzyaaal9(xyzzyaaad9)
enddo
xyzzyaaak9=0.d0
xyzzyaaal9=0.d0
if(xyzzyaaaf9<rcut_loc(xyzzyaaac9))then
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_loc_full(xyzzyaaab9,1,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzzyaaam9&
&,xyzzyaaan9,xyzzyaaao9)
do xyzzyaaad9=1,naxis_forces
xyzzyaaak9(xyzzyaaad9)=xyzzyaaan9*xyzzyaaai9(xyzzyaaad9)
enddo
else
do xyzzyaaad9=1,naxis_forces
xyzzyaaak9(xyzzyaaad9)=valence_charge(xyzzyaaaa9)*xyzzyaaaj9(xyzzyaaad&
&9)
enddo
endif
if(xyzzyaaaf9<rcut_non_loc(xyzzyaaac9))then
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_non_loc_full(xyzzyaaab9,0,1,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzz&
&yaaam9,xyzzyaaan9,xyzzyaaao9)
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_non_loc_full(xyzzyaaab9,2,1,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzz&
&yaaap9,xyzzyaaaq9,xyzzyaaar9)
do xyzzyaaad9=1,naxis_forces
xyzzyaaal9(xyzzyaaad9)=xyzzyaaal9(xyzzyaaad9)-xyzzyaaaa1(xyzzyaaad9,xy&
&zzyaaaa9)*xyzzyaaam9-xyzzyaaab1(xyzzyaaad9,xyzzyaaaa9)*xyzzyaaan9
xyzzyaaal9(xyzzyaaad9)=xyzzyaaal9(xyzzyaaad9)-5.d0*(xyzzyaaap9*(xyzzya&
&aaf1(xyzzyaaad9,xyzzyaaaa9)+xyzzyaaah1(xyzzyaaad9,xyzzyaaaa9))+xyzzya&
&aag1(xyzzyaaad9,xyzzyaaaa9)*xyzzyaaaq9)
enddo
endif
do xyzzyaaad9=1,naxis_forces
xyzzyaaat9(xyzzyaaad9)=xyzzyaaak9(xyzzyaaad9)+xyzzyaaal9(xyzzyaaad9)
enddo
endif
xyzzyaaak9=0.d0
xyzzyaaal9=0.d0
if(xyzzyaaaf9<rcut_loc(xyzzyaaac9))then
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_loc_full(xyzzyaaab9,2,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzzyaaam9&
&,xyzzyaaan9,xyzzyaaao9)
do xyzzyaaad9=1,naxis_forces
xyzzyaaak9(xyzzyaaad9)=xyzzyaaan9*xyzzyaaai9(xyzzyaaad9)
enddo
else
do xyzzyaaad9=1,naxis_forces
xyzzyaaak9(xyzzyaaad9)=valence_charge(xyzzyaaaa9)*xyzzyaaaj9(xyzzyaaad&
&9)
enddo
endif
if(xyzzyaaaf9<rcut_non_loc(xyzzyaaac9))then
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_non_loc_full(xyzzyaaab9,0,2,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzz&
&yaaam9,xyzzyaaan9,xyzzyaaao9)
call interp_nev_with_derivs(pp_radial_grid(xyzzyaaab9,xyzzyaaac9),fcoe&
&ff_non_loc_full(xyzzyaaab9,1,2,xyzzyaaac9),xyzzyaaae9,xyzzyaaaf9,xyzz&
&yaaap9,xyzzyaaaq9,xyzzyaaar9)
do xyzzyaaad9=1,naxis_forces
xyzzyaaal9(xyzzyaaad9)=xyzzyaaal9(xyzzyaaad9)-xyzzyaaaa1(xyzzyaaad9,xy&
&zzyaaaa9)*xyzzyaaam9-xyzzyaaab1(xyzzyaaad9,xyzzyaaaa9)*xyzzyaaan9
xyzzyaaal9(xyzzyaaad9)=xyzzyaaal9(xyzzyaaad9)-3.d0*(xyzzyaaap9*(xyzzya&
&aac1(xyzzyaaad9,xyzzyaaaa9)+xyzzyaaae1(xyzzyaaad9,xyzzyaaaa9))+xyzzya&
&aad1(xyzzyaaad9,xyzzyaaaa9)*xyzzyaaaq9)
enddo
endif
do xyzzyaaad9=1,naxis_forces
xyzzyaaau9(xyzzyaaad9)=xyzzyaaak9(xyzzyaaad9)+xyzzyaaal9(xyzzyaaad9)
enddo
do xyzzyaaad9=1,naxis_forces
if(forces_info>2)then
xyzzyaaai1(1,xyzzyaaad9,xyzzyaaaa9)=xyzzyaaai1(1,xyzzyaaad9,xyzzyaaaa9&
&)+xyzzyaaas9(xyzzyaaad9)+xyzzyaaaj1(xyzzyaaad9,xyzzyaaaa9)
xyzzyaaai1(2,xyzzyaaad9,xyzzyaaaa9)=xyzzyaaai1(2,xyzzyaaad9,xyzzyaaaa9&
&)+xyzzyaaat9(xyzzyaaad9)+xyzzyaaaj1(xyzzyaaad9,xyzzyaaaa9)
endif
xyzzyaaai1(3,xyzzyaaad9,xyzzyaaaa9)=xyzzyaaai1(3,xyzzyaaad9,xyzzyaaaa9&
&)+xyzzyaaau9(xyzzyaaad9)+xyzzyaaaj1(xyzzyaaad9,xyzzyaaaa9)
enddo
enddo
end subroutine xyzzyaaay1
subroutine xyzzyaaaz1(is0)
use slaarnacs,only : define_config_oneion,wfn_ratio
implicit none
integer,intent(in) :: is0
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10
real(dp) xyzzyaaad10(3),xyzzyaaae10,etot,ecomps(25),xyzzyaaaf10,xyzzya&
&aag10,xyzzyaaah10,xyzzyaaai10
complex(dp) ratio
logical isnan,isinf
xyzzyaaac10=xyzzyaaau1(is0)
do xyzzyaaaa10=1,nitot_forces
do xyzzyaaab10=1,naxis_forces
xyzzyaaad10=rion(:,xyzzyaaaa10)
xyzzyaaad10(xyzzyaaab10)=xyzzyaaad10(xyzzyaaab10)+xyzzyaaao1
call make_representable(xyzzyaaad10(xyzzyaaab10))
xyzzyaaap1(xyzzyaaab10,xyzzyaaaa10)=xyzzyaaad10(xyzzyaaab10)-rion(xyzz&
&yaaab10,xyzzyaaaa10)
call define_config_oneion(xyzzyaaaa10,is0,xyzzyaaac10,xyzzyaaad10)
call wfn_ratio(is0,xyzzyaaac10,0,ratio=ratio,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)call errstop('CALC_FORCES_WFN_ELEMENTS','Floating-po&
&int exception reported when calculating wave function ratio. Don''t k&
&now what to do, so stopping. Please file a bug report.')
xyzzyaaae10=dble(ratio)
xyzzyaaan1(xyzzyaaab10,xyzzyaaaa10)=log(abs(xyzzyaaae10))/xyzzyaaap1(x&
&yzzyaaab10,xyzzyaaaa10)
call eval_local_energy(xyzzyaaac10,etot=etot,ecomps=ecomps,fix_nl_grid&
&=.true.,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)call errstop('CALC_FORCES_WFN_ELEMENTS','Floating-po&
&int exception reported when calculating energy. Don''t know what to d&
&o, so stopping. Please file a bug report.')
xyzzyaaaf10=ecomps(i_kei)
xyzzyaaag10=ecomps(i_potinl)
xyzzyaaak1(xyzzyaaab10,xyzzyaaaa10)=xyzzyaaaf10*xyzzyaaae10
xyzzyaaai10=0.d0
xyzzyaaah10=0.d0
if(have_ppots)then
if(.not.use_tmove)then
xyzzyaaai10=xyzzyaaag10
else
call v_non_local_tmove(xyzzyaaai10,xyzzyaaah10)
endif
endif
xyzzyaaam1(xyzzyaaab10,xyzzyaaaa10)=xyzzyaaai10*xyzzyaaae10
xyzzyaaal1(xyzzyaaab10,xyzzyaaaa10)=xyzzyaaah10*xyzzyaaae10
enddo
enddo
end subroutine xyzzyaaaz1
subroutine eval_local_forces(is0,eloc,kei,poti_nl_plus,poti_nl_minus,f&
&comps)
implicit none
integer,intent(in) :: is0
real(dp),intent(in) :: eloc,poti_nl_plus,poti_nl_minus,kei
real(dp),intent(out) :: fcomps(nfcomps)
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11
real(dp) xyzzyaaae11
xyzzyaaai1=0.d0
call get_rsele(is0)
call get_eivecs(is0)
do xyzzyaaad11=1,netot
call xyzzyaaaw1(xyzzyaaad11,rele_scr(1,xyzzyaaad11,is0),eivecs_scr(1,1&
&,xyzzyaaad11,is0))
call xyzzyaaay1(eivecs_scr(1,1,xyzzyaaad11,is0),xyzzyaaai1)
enddo
call xyzzyaaaz1(is0)
xyzzyaaac11=0
fcomps=0.d0
do xyzzyaaab11=1,nitot_forces
do xyzzyaaaa11=1,naxis_forces
xyzzyaaae11=xyzzyaaan1(xyzzyaaaa11,xyzzyaaab11)
fcomps(xyzzyaaac11+1)=xyzzyaaae11
fcomps(xyzzyaaac11+2)=eloc*xyzzyaaae11
fcomps(xyzzyaaac11+3)=kei*xyzzyaaae11
fcomps(xyzzyaaac11+4)=poti_nl_plus*xyzzyaaae11
fcomps(xyzzyaaac11+5)=poti_nl_minus*xyzzyaaae11
fcomps(xyzzyaaac11+6)=(xyzzyaaak1(xyzzyaaaa11,xyzzyaaab11)-kei)/xyzzya&
&aap1(xyzzyaaaa11,xyzzyaaab11)
fcomps(xyzzyaaac11+7)=(xyzzyaaam1(xyzzyaaaa11,xyzzyaaab11)-poti_nl_plu&
&s)/xyzzyaaap1(xyzzyaaaa11,xyzzyaaab11)
fcomps(xyzzyaaac11+8)=(xyzzyaaal1(xyzzyaaaa11,xyzzyaaab11)-poti_nl_min&
&us)/xyzzyaaap1(xyzzyaaaa11,xyzzyaaab11)
fcomps(xyzzyaaac11+9)=xyzzyaaai1(3,xyzzyaaaa11,xyzzyaaab11)
if(forces_info<=2)then
fcomps(xyzzyaaac11+10)=0.d0
fcomps(xyzzyaaac11+11)=0.d0
else
fcomps(xyzzyaaac11+10)=xyzzyaaai1(2,xyzzyaaaa11,xyzzyaaab11)
fcomps(xyzzyaaac11+11)=xyzzyaaai1(1,xyzzyaaaa11,xyzzyaaab11)
endif
xyzzyaaac11=xyzzyaaac11+nfterms
enddo
enddo
end subroutine eval_local_forces
subroutine forces_to_future(fcomps,pureitems)
implicit none
real(dp),intent(in) :: fcomps(nfcomps)
real(dp),intent(inout) :: pureitems(*)
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12
xyzzyaaaa12=0
xyzzyaaab12=0
do xyzzyaaad12=1,nitot_forces
do xyzzyaaac12=1,naxis_forces
do xyzzyaaae12=1,9
pureitems(xyzzyaaab12+xyzzyaaae12)=fcomps(xyzzyaaaa12+xyzzyaaae12)
enddo
if(forces_info>2)then
pureitems(xyzzyaaab12+10)=fcomps(xyzzyaaaa12+10)
pureitems(xyzzyaaab12+11)=fcomps(xyzzyaaaa12+11)
else
pureitems(xyzzyaaab12+10:xyzzyaaab12+11)=0.d0
endif
xyzzyaaaa12=xyzzyaaaa12+nfterms
xyzzyaaab12=xyzzyaaab12+11
enddo
enddo
end subroutine forces_to_future
subroutine future_to_forces(pureitems,fcomps)
implicit none
real(dp),intent(in) :: pureitems(*)
real(dp),intent(inout) :: fcomps(nfcomps)
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13,xyzzyaaae13
xyzzyaaaa13=11
xyzzyaaab13=0
do xyzzyaaad13=1,nitot_forces
do xyzzyaaac13=1,naxis_forces
do xyzzyaaae13=1,9
fcomps(xyzzyaaaa13+xyzzyaaae13)=pureitems(xyzzyaaab13+xyzzyaaae13)
enddo
if(forces_info>2)then
fcomps(xyzzyaaaa13+10)=pureitems(xyzzyaaab13+10)
fcomps(xyzzyaaaa13+11)=pureitems(xyzzyaaab13+11)
endif
xyzzyaaaa13=xyzzyaaaa13+nfterms
xyzzyaaab13=xyzzyaaab13+11
enddo
enddo
end subroutine future_to_forces
end module slaarnaab

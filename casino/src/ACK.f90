module slaarnack
use dsp
use slaarnabg
use store
use slaarnaag,   only : root_one_over_pi,root_eight_over_pi
use format_utils,only : wout
use parallel,    only : am_master
use run_control, only : errstop,errwarn
implicit none
private
public setup_special_wfn,special_orb_eval,get_special_wfn_orbmap,get_s&
&pecial_wfn_orbdesc,get_special_wfn_ndesc
integer xyzzyaaaa1
integer,parameter :: xyzzyaaab1=0,xyzzyaaac1=0
real(dp),parameter :: xyzzyaaad1=1.4d0,xyzzyaaae1=1.65d0
contains
subroutine setup_special_wfn(eionion)
implicit none
real(dp),intent(out) :: eionion
select case(trim(atom_basis_type))
case('non_int_he')
call xyzzyaaaf1(eionion)
case('h2')
call xyzzyaaag1(eionion)
case('h3plus')
call xyzzyaaah1(eionion)
case default
if(am_master)call errstop('SETUP_SPECIAL_WFN','Special wave function n&
&ot recognized. Bug.')
end select
end subroutine setup_special_wfn
subroutine special_orb_eval(rvec,jspin,norb,orbmask,val,fsd,orbval,orb&
&grad,orblap,orbsderivs)
implicit none
integer,intent(in) :: jspin,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd,orbmask(norb)
select case(xyzzyaaaa1)
case(1)
call xyzzyaaai1(rvec,norb,val,fsd,orbval,orbgrad,orblap,orbsderivs)
case(2)
call xyzzyaaaj1(rvec,norb,val,fsd,orbval,orbgrad,orblap,orbsderivs)
case(3)
call xyzzyaaak1(rvec,norb,val,fsd,orbval,orbgrad,orblap,orbsderivs)
case default
call errstop('SPECIAL_ORB_EVAL','Special wave function code not recogn&
&ized.')
end select
end subroutine special_orb_eval
subroutine xyzzyaaaf1(eionion)
use slaarnabg
implicit none
real(dp),intent(out) :: eionion
if(am_master)then
call wout()
call wout('Exact orbitals for non-interacting helium atom.')
call wout()
endif
if(trim(interaction)/='none')call errwarn('NONINTHE_SETUP','electron-e&
&lectron interaction not deactivated. Set keyword INTERACTION to ''non&
&e'' to turn it off.')
pa1=(/500.d0,0.d0,0.d0/)
pa2=(/0.d0,500.d0,0.d0/)
pa3=(/0.d0,0.d0,500.d0/)
allocate(basis(3,1),atno(1))
nbasis=1
basis(1:3,1)=0.d0
atno(1)=2
ndet=3
xyzzyaaaa1=1
eionion=0.d0
end subroutine xyzzyaaaf1
subroutine xyzzyaaag1(eionion)
use slaarnabg
implicit none
real(dp),intent(out) :: eionion
if(am_master)then
call wout()
call wout('Parameter-less exponential orbitals for H2 molecule.')
endif
pa1=(/500.d0,0.d0,0.d0/)
pa2=(/0.d0,500.d0,0.d0/)
pa3=(/0.d0,0.d0,500.d0/)
nbasis=2
allocate(basis(3,nbasis),atno(nbasis))
basis(1:3,1)=(/0.d0,0.d0,-xyzzyaaad1*0.5d0/)
basis(1:3,2)=(/0.d0,0.d0,xyzzyaaad1*0.5d0/)
atno(1:2)=(/1,1/)
ndet=1
xyzzyaaaa1=2
eionion=1.d0/xyzzyaaad1
end subroutine xyzzyaaag1
subroutine xyzzyaaah1(eionion)
use slaarnabg
implicit none
real(dp),intent(out) :: eionion
if(am_master)then
call wout()
call wout('Parameter-less exponential orbitals for H3+ molecule.')
endif
pa1=(/500.d0,0.d0,0.d0/)
pa2=(/0.d0,500.d0,0.d0/)
pa3=(/0.d0,0.d0,500.d0/)
nbasis=3
allocate(basis(3,nbasis),atno(nbasis))
basis(1:3,1)=(/0.825d0,-0.476313972d0,0.d0/)
basis(1:3,2)=(/-0.825d0,-0.476313972d0,0.d0/)
basis(1:3,3)=(/0.d0,0.952627944d0,0.d0/)
atno(1:3)=(/1,1,1/)
ndet=1
xyzzyaaaa1=3
eionion=3.d0/xyzzyaaae1
end subroutine xyzzyaaah1
subroutine xyzzyaaai1(rvec,norb,val,fsd,orbval,orbgrad,orblap,orbsderi&
&vs)
implicit none
integer,intent(in) :: norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd
real(dp) xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7,xyzzya&
&aaf7,xyzzyaaag7
xyzzyaaaa7=sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
xyzzyaaab7=1.d0/xyzzyaaaa7
xyzzyaaac7=exp(-xyzzyaaaa7)
xyzzyaaad7=xyzzyaaac7*xyzzyaaac7
if(val)then
orbval(1,1)=root_eight_over_pi*xyzzyaaad7
orbval(2,1)=root_one_over_pi*(1.d0-xyzzyaaaa7)*xyzzyaaac7
if(complex_wf)orbval(:,2)=0.d0
endif
if(fsd)then
xyzzyaaae7=-2.d0*root_eight_over_pi*xyzzyaaad7*xyzzyaaab7
orbgrad(:,1,1)=xyzzyaaae7*rvec
xyzzyaaaf7=root_one_over_pi*(xyzzyaaaa7-2.d0)*xyzzyaaac7*xyzzyaaab7
orbgrad(:,2,1)=xyzzyaaaf7*rvec
orblap(1,1)=root_eight_over_pi*4.d0*xyzzyaaad7*(1.d0-xyzzyaaab7)
orblap(2,1)=root_one_over_pi*xyzzyaaac7*(5.d0-4.d0*xyzzyaaab7-xyzzyaaa&
&a7)
if(present(orbsderivs))then
xyzzyaaae7=2.d0*root_eight_over_pi*xyzzyaaab7*xyzzyaaad7
xyzzyaaaf7=xyzzyaaab7*(xyzzyaaab7+2.d0)
orbsderivs(1,1,1)=xyzzyaaae7*(xyzzyaaaf7*rvec(1)**2-1)
orbsderivs(2,1,1)=xyzzyaaae7*(xyzzyaaaf7*rvec(2)**2-1)
orbsderivs(3,1,1)=xyzzyaaae7*(xyzzyaaaf7*rvec(3)**2-1)
xyzzyaaae7=xyzzyaaae7*xyzzyaaaf7
orbsderivs(4,1,1)=xyzzyaaae7*rvec(1)*rvec(2)
orbsderivs(5,1,1)=xyzzyaaae7*rvec(1)*rvec(3)
orbsderivs(6,1,1)=xyzzyaaae7*rvec(2)*rvec(3)
xyzzyaaae7=root_one_over_pi*xyzzyaaac7*xyzzyaaab7
xyzzyaaaf7=-1.d0+2.d0*xyzzyaaab7*(1.d0+xyzzyaaab7)
xyzzyaaag7=xyzzyaaaa7-2.d0
orbsderivs(1,2,1)=xyzzyaaae7*(xyzzyaaaf7*rvec(1)**2+xyzzyaaag7)
orbsderivs(2,2,1)=xyzzyaaae7*(xyzzyaaaf7*rvec(2)**2+xyzzyaaag7)
orbsderivs(3,2,1)=xyzzyaaae7*(xyzzyaaaf7*rvec(3)**2+xyzzyaaag7)
xyzzyaaae7=xyzzyaaae7*xyzzyaaaf7
orbsderivs(4,2,1)=xyzzyaaae7*rvec(1)*rvec(2)
orbsderivs(5,2,1)=xyzzyaaae7*rvec(1)*rvec(3)
orbsderivs(6,2,1)=xyzzyaaae7*rvec(2)*rvec(3)
endif
if(complex_wf)then
orbgrad(:,:,2)=0.d0
orblap(:,2)=0.d0
if(present(orbsderivs))orbsderivs(:,:,2)=0.d0
endif
endif
end subroutine xyzzyaaai1
subroutine xyzzyaaaj1(rvec,norb,val,fsd,orbval,orbgrad,orblap,orbsderi&
&vs)
implicit none
integer,intent(in) :: norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd
real(dp) xyzzyaaaa8(3),xyzzyaaab8(3),xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,&
&xyzzyaaaf8,xyzzyaaag8,xyzzyaaah8
xyzzyaaaa8=rvec-(/0.d0,0.d0,-xyzzyaaad1*0.5d0/)
xyzzyaaac8=sqrt(xyzzyaaaa8(1)**2+xyzzyaaaa8(2)**2+xyzzyaaaa8(3)**2)
xyzzyaaae8=exp(-xyzzyaaac8)
xyzzyaaab8=rvec-(/0.d0,0.d0,xyzzyaaad1*0.5d0/)
xyzzyaaad8=sqrt(xyzzyaaab8(1)**2+xyzzyaaab8(2)**2+xyzzyaaab8(3)**2)
xyzzyaaaf8=exp(-xyzzyaaad8)
if(val)orbval(1,1)=xyzzyaaae8+xyzzyaaaf8
if(fsd)then
if(xyzzyaaac8/=0.d0.and.xyzzyaaad8/=0.d0)then
xyzzyaaag8=1.d0/xyzzyaaac8
xyzzyaaah8=1.d0/xyzzyaaad8
orbgrad(1:3,1,1)=-xyzzyaaae8*xyzzyaaaa8*xyzzyaaag8-xyzzyaaaf8*xyzzyaaa&
&b8*xyzzyaaah8
orblap(1,1)=xyzzyaaae8*(1.d0-2*xyzzyaaag8)+xyzzyaaaf8*(1.d0-2*xyzzyaaa&
&h8)
else
orbgrad(1:3,1,1)=0.d0
orblap(1,1)=0.d0
endif
endif
end subroutine xyzzyaaaj1
subroutine xyzzyaaak1(rvec,norb,val,fsd,orbval,orbgrad,orblap,orbsderi&
&vs)
implicit none
integer,intent(in) :: norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd
real(dp) xyzzyaaaa9(3),xyzzyaaab9(3),xyzzyaaac9(3),xyzzyaaad9,xyzzyaaa&
&e9,xyzzyaaaf9,xyzzyaaag9,xyzzyaaah9,xyzzyaaai9,xyzzyaaaj9,xyzzyaaak9,&
&xyzzyaaal9
xyzzyaaaa9=rvec-(/0.825d0,-0.476313972d0,0.d0/)
xyzzyaaad9=sqrt(xyzzyaaaa9(1)**2+xyzzyaaaa9(2)**2+xyzzyaaaa9(3)**2)
xyzzyaaag9=exp(-xyzzyaaad9)
xyzzyaaab9=rvec-(/-0.825d0,-0.476313972d0,0.d0/)
xyzzyaaae9=sqrt(xyzzyaaab9(1)**2+xyzzyaaab9(2)**2+xyzzyaaab9(3)**2)
xyzzyaaah9=exp(-xyzzyaaae9)
xyzzyaaac9=rvec-(/0.d0,0.952627944d0,0.d0/)
xyzzyaaaf9=sqrt(xyzzyaaac9(1)**2+xyzzyaaac9(2)**2+xyzzyaaac9(3)**2)
xyzzyaaai9=exp(-xyzzyaaaf9)
if(val)orbval(1,1)=xyzzyaaag9+xyzzyaaah9+xyzzyaaai9
if(fsd)then
if(xyzzyaaad9/=0.d0.and.xyzzyaaae9/=0.d0.and.xyzzyaaaf9/=0.d0)then
xyzzyaaaj9=1.d0/xyzzyaaad9
xyzzyaaak9=1.d0/xyzzyaaae9
xyzzyaaal9=1.d0/xyzzyaaaf9
orbgrad(1:3,1,1)=-xyzzyaaag9*xyzzyaaaa9*xyzzyaaaj9-xyzzyaaah9*xyzzyaaa&
&b9*xyzzyaaak9-xyzzyaaai9*xyzzyaaac9*xyzzyaaal9
orblap(1,1)=xyzzyaaag9*(1.d0-2*xyzzyaaaj9)+xyzzyaaah9*(1.d0-2*xyzzyaaa&
&k9)+xyzzyaaai9*(1.d0-2*xyzzyaaal9)
else
orbgrad(1:3,1,1)=0.d0
orblap(1,1)=0.d0
endif
endif
end subroutine xyzzyaaak1
subroutine get_special_wfn_orbmap(row_offset,norb,orbmap)
implicit none
integer,intent(inout) :: row_offset(nspin),norb,orbmap(nemax,nspin,nde&
&t)
select case(trim(atom_basis_type))
case('non_int_he')
orbmap(row_offset(1)+1,1,1)=norb+1
orbmap(row_offset(1)+1,1,2)=norb+1
orbmap(row_offset(1)+1,1,3)=norb+2
orbmap(row_offset(2)+1,2,1)=norb+1
orbmap(row_offset(2)+1,2,2)=norb+2
orbmap(row_offset(2)+1,2,3)=norb+1
row_offset(1:2)=row_offset(1:2)+1
norb=norb+2
case('h2')
orbmap(row_offset(1)+1,1,1)=norb+1
orbmap(row_offset(2)+1,2,1)=norb+1
row_offset(1:2)=row_offset(1:2)+1
norb=norb+1
case('h3plus')
orbmap(row_offset(1)+1,1,1)=norb+1
orbmap(row_offset(2)+1,2,1)=norb+1
row_offset(1:2)=row_offset(1:2)+1
norb=norb+1
case default
if(am_master)call errstop('GET_SPECIAL_WFN_ORBMAP','Special wave funct&
&ion not recognized. Bug.')
end select
end subroutine get_special_wfn_orbmap
subroutine get_special_wfn_ndesc(ndesc_int,ndesc_dp)
implicit none
integer,intent(inout) :: ndesc_int,ndesc_dp
ndesc_int=xyzzyaaab1
ndesc_dp=xyzzyaaac1
end subroutine get_special_wfn_ndesc
subroutine get_special_wfn_orbdesc(norb,ndesc_int,ndesc_dp,orbdesc_int&
&,orbdesc_dp)
implicit none
integer,intent(in) :: norb,ndesc_int,ndesc_dp
integer,intent(inout) :: orbdesc_int(ndesc_int,norb)
real(dp),intent(inout) :: orbdesc_dp(ndesc_dp,norb)
end subroutine get_special_wfn_orbdesc
end module slaarnack

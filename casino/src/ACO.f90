module slaarnaco
use dsp
use slaarnabg
use slaarnaca,      only : l_local,zion,rcut_non_loc,cppalpha,cpprbare&
&e,cpprbaree_sq,cpprbarl,cpprbarl_sq,have_veep,is_cpp,haverbarl,ion_fi&
&elds
use run_control,only : errstop,timer
implicit none
private
public compute_vcpp,compute_ecppii,compute_ii_fields
contains
subroutine compute_vcpp(elec,rele,eivecs,vcppei,vcppe,vcppee,field,ecp&
&p_nl)
use slaarnabs
use store
implicit none
integer,intent(in) :: elec
real(dp),intent(in) :: rele(3,netot),eivecs(4,nitot,netot)
real(dp),intent(out) :: vcppei,vcppe,vcppee,ecpp_nl
real(dp),intent(inout) :: field(1:3,nitot)
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2
real(dp) xyzzyaaaf2(3),xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyz&
&zyaaak2,xyzzyaaal2(3),xyzzyaaam2,xyzzyaaan2,xyzzyaaao2,xyzzyaaap2,xyz&
&zyaaaq2(3),xyzzyaaar2,xyzzyaaas2(1:3),xyzzyaaat2,xyzzyaaau2
real(dp),parameter :: xyzzyaaav2=3.809123148d0
vcppei=0.d0
vcppee=0.d0
vcppe=0.d0
ecpp_nl=0.d0
if(.not.have_veep)return
call timer('CORE_POL_EE',.true.)
if(isperiodic)then
do xyzzyaaaa2=1,nitot
xyzzyaaab2=iontype(xyzzyaaaa2)
if(is_cpp(xyzzyaaab2))then
xyzzyaaar2=xyzzyaaav2*cpprbaree(xyzzyaaab2)
xyzzyaaas2(1:3)=xyzzyaaav2*cpprbarl(1:3,xyzzyaaab2)
xyzzyaaah2=eivecs(4,xyzzyaaaa2,elec)
xyzzyaaag2=xyzzyaaah2*xyzzyaaah2
xyzzyaaap2=cpprbaree_sq(xyzzyaaab2)
if(xyzzyaaah2<xyzzyaaar2)then
xyzzyaaai2=(1.d0-exp(-xyzzyaaag2/xyzzyaaap2))**2
field(1:3,xyzzyaaaa2)=xyzzyaaai2*field(1:3,xyzzyaaaa2)
else
xyzzyaaai2=1.d0
endif
xyzzyaaau2=xyzzyaaai2*xyzzyaaai2
if(xyzzyaaah2<maxval(xyzzyaaas2))then
if(haverbarl(xyzzyaaab2))then
xyzzyaaak2=-0.5d0*cppalpha(xyzzyaaab2)/(xyzzyaaag2*xyzzyaaag2)
xyzzyaaat2=0.d0
if(xyzzyaaah2<rcut_non_loc(xyzzyaaab2))then
do xyzzyaaae2=1,3
if(xyzzyaaah2<xyzzyaaas2(xyzzyaaae2))then
xyzzyaaap2=cpprbarl_sq(xyzzyaaae2,xyzzyaaab2)
xyzzyaaai2=(1.d0-exp(-xyzzyaaag2/xyzzyaaap2))**2
else
xyzzyaaai2=1d0
endif
xyzzyaaaq2(xyzzyaaae2)=xyzzyaaai2*xyzzyaaai2
enddo
do xyzzyaaae2=1,3
xyzzyaaat2=xyzzyaaat2+(xyzzyaaaq2(xyzzyaaae2)-xyzzyaaau2)*nl_am_proj(x&
&yzzyaaae2-1,xyzzyaaaa2)
enddo
endif
ecpp_nl=ecpp_nl+xyzzyaaat2*xyzzyaaak2
endif
endif
endif
enddo
else
do xyzzyaaaa2=1,nitot
xyzzyaaab2=iontype(xyzzyaaaa2)
if(is_cpp(xyzzyaaab2))then
xyzzyaaar2=xyzzyaaav2*cpprbaree(xyzzyaaab2)
xyzzyaaas2(1:3)=xyzzyaaav2*cpprbarl(1:3,xyzzyaaab2)
xyzzyaaaf2(1:3)=eivecs(1:3,xyzzyaaaa2,elec)
xyzzyaaah2=eivecs(4,xyzzyaaaa2,elec)
xyzzyaaag2=xyzzyaaah2*xyzzyaaah2
xyzzyaaap2=cpprbaree_sq(xyzzyaaab2)
if(xyzzyaaah2<xyzzyaaar2)then
xyzzyaaai2=(1.d0-exp(-xyzzyaaag2/xyzzyaaap2))**2
else
xyzzyaaai2=1.d0
endif
xyzzyaaak2=-0.5d0*cppalpha(xyzzyaaab2)*xyzzyaaai2/(xyzzyaaag2*xyzzyaaa&
&h2)
xyzzyaaat2=0.d0
do xyzzyaaac2=1,nitot
if(xyzzyaaac2==xyzzyaaaa2)cycle
xyzzyaaal2(1:3)=rion(1:3,xyzzyaaac2)-rion(1:3,xyzzyaaaa2)
xyzzyaaam2=xyzzyaaal2(1)*xyzzyaaal2(1)+xyzzyaaal2(2)*xyzzyaaal2(2)+xyz&
&zyaaal2(3)*xyzzyaaal2(3)
xyzzyaaan2=sqrt(xyzzyaaam2)
xyzzyaaao2=xyzzyaaaf2(1)*xyzzyaaal2(1)+xyzzyaaaf2(2)*xyzzyaaal2(2)+xyz&
&zyaaaf2(3)*xyzzyaaal2(3)
xyzzyaaat2=xyzzyaaat2+zion(iontype(xyzzyaaac2))*xyzzyaaao2/(xyzzyaaam2&
&*xyzzyaaan2)
enddo
vcppei=vcppei-2.d0*xyzzyaaat2*xyzzyaaak2
xyzzyaaat2=0.d0
do xyzzyaaac2=1,netot
if(xyzzyaaac2==elec)cycle
xyzzyaaal2(1:3)=rele(1:3,xyzzyaaac2)-rion(1:3,xyzzyaaaa2)
xyzzyaaam2=xyzzyaaal2(1)*xyzzyaaal2(1)+xyzzyaaal2(2)*xyzzyaaal2(2)+xyz&
&zyaaal2(3)*xyzzyaaal2(3)
xyzzyaaan2=sqrt(xyzzyaaam2)
if(xyzzyaaan2<xyzzyaaar2)then
xyzzyaaaj2=(1.d0-exp(-xyzzyaaam2/xyzzyaaap2))**2
else
xyzzyaaaj2=1.d0
endif
xyzzyaaao2=xyzzyaaaf2(1)*xyzzyaaal2(1)+xyzzyaaaf2(2)*xyzzyaaal2(2)+xyz&
&zyaaaf2(3)*xyzzyaaal2(3)
xyzzyaaat2=xyzzyaaat2+xyzzyaaaj2*xyzzyaaao2/(xyzzyaaam2*xyzzyaaan2)
enddo
vcppee=vcppee+xyzzyaaat2*xyzzyaaak2
if(haverbarl(xyzzyaaab2))then
xyzzyaaak2=-0.5d0*cppalpha(xyzzyaaab2)/(xyzzyaaag2*xyzzyaaag2)
xyzzyaaat2=0.d0
xyzzyaaad2=l_local(xyzzyaaab2)+1
if(xyzzyaaah2<rcut_non_loc(xyzzyaaab2))then
do xyzzyaaae2=1,3
if(xyzzyaaah2<xyzzyaaas2(xyzzyaaae2))then
xyzzyaaap2=cpprbarl_sq(xyzzyaaae2,xyzzyaaab2)
xyzzyaaai2=(1.d0-exp(-xyzzyaaag2/xyzzyaaap2))**2
xyzzyaaaq2(xyzzyaaae2)=xyzzyaaai2*xyzzyaaai2
else
xyzzyaaaq2(xyzzyaaae2)=1.d0
endif
enddo
xyzzyaaat2=xyzzyaaaq2(xyzzyaaad2)
if(xyzzyaaad2==1)then
xyzzyaaat2=xyzzyaaat2+(xyzzyaaaq2(2)-xyzzyaaaq2(1))*nl_am_proj(1,xyzzy&
&aaaa2)
xyzzyaaat2=xyzzyaaat2+(xyzzyaaaq2(3)-xyzzyaaaq2(1))*nl_am_proj(2,xyzzy&
&aaaa2)
elseif(xyzzyaaad2==2)then
xyzzyaaat2=xyzzyaaat2+(xyzzyaaaq2(1)-xyzzyaaaq2(2))*nl_am_proj(0,xyzzy&
&aaaa2)
xyzzyaaat2=xyzzyaaat2+(xyzzyaaaq2(3)-xyzzyaaaq2(2))*nl_am_proj(2,xyzzy&
&aaaa2)
else
xyzzyaaat2=xyzzyaaat2+(xyzzyaaaq2(1)-xyzzyaaaq2(3))*nl_am_proj(0,xyzzy&
&aaaa2)
xyzzyaaat2=xyzzyaaat2+(xyzzyaaaq2(2)-xyzzyaaaq2(3))*nl_am_proj(1,xyzzy&
&aaaa2)
endif
vcppe=vcppe+xyzzyaaat2*xyzzyaaak2
else
if(xyzzyaaah2<xyzzyaaas2(xyzzyaaad2))then
xyzzyaaap2=cpprbarl_sq(xyzzyaaad2,xyzzyaaab2)
xyzzyaaai2=(1.d0-exp(-xyzzyaaag2/xyzzyaaap2))**2
vcppe=vcppe+xyzzyaaai2*xyzzyaaai2*xyzzyaaak2
else
vcppe=vcppe+xyzzyaaak2
endif
endif
endif
endif
enddo
endif
call timer('CORE_POL_EE',.false.)
end subroutine compute_vcpp
subroutine compute_ecppii(ecppii)
implicit none
real(dp),intent(out) :: ecppii
integer xyzzyaaaa3,xyzzyaaab3
real(dp) xyzzyaaac3,xyzzyaaad3(3),xyzzyaaae3(3),xyzzyaaaf3(4),xyzzyaaa&
&g3
ecppii=0.d0
if(have_veep)then
do xyzzyaaaa3=1,nitot
if(.not.is_cpp(iontype(xyzzyaaaa3)))cycle
xyzzyaaac3=-0.5d0*cppalpha(iontype(xyzzyaaaa3))
xyzzyaaad3(1:3)=0.d0
xyzzyaaae3(1:3)=rion(1:3,xyzzyaaaa3)
do xyzzyaaab3=1,nitot
if(xyzzyaaab3==xyzzyaaaa3)cycle
xyzzyaaaf3(1:3)=rion(1:3,xyzzyaaab3)-xyzzyaaae3(1:3)
xyzzyaaag3=(xyzzyaaaf3(1)*xyzzyaaaf3(1)+xyzzyaaaf3(2)*xyzzyaaaf3(2)+xy&
&zzyaaaf3(3)*xyzzyaaaf3(3))**1.5d0
xyzzyaaad3(1:3)=xyzzyaaad3(1:3)+zion(iontype(xyzzyaaab3))*xyzzyaaaf3(1&
&:3)/xyzzyaaag3
enddo
ecppii=ecppii+xyzzyaaac3*(xyzzyaaad3(1)*xyzzyaaad3(1)+xyzzyaaad3(2)*xy&
&zzyaaad3(2)+xyzzyaaad3(3)*xyzzyaaad3(3))
enddo
endif
end subroutine compute_ecppii
subroutine compute_ii_fields
use slaarnabq
use slaarnabk, only : ewald_2d,ewald_3d
implicit none
integer xyzzyaaaa4,xyzzyaaab4
real(dp) xyzzyaaac4,xyzzyaaad4,xyzzyaaae4(4),xyzzyaaaf4(3)
if(have_veep)then
ion_fields(:,:)=0.d0
do xyzzyaaaa4=1,nitot
if(.not.is_cpp(iontype(xyzzyaaaa4)))cycle
do xyzzyaaab4=1,nitot
if(xyzzyaaaa4==xyzzyaaab4)cycle
xyzzyaaac4=zion(iontype(xyzzyaaab4))
xyzzyaaae4(1:3)=rion(1:3,xyzzyaaab4)-rion(1:3,xyzzyaaaa4)
call minimum_image(3,1,xyzzyaaae4)
xyzzyaaae4(4)=sqrt(xyzzyaaae4(1)*xyzzyaaae4(1)+xyzzyaaae4(2)*xyzzyaaae&
&4(2)+xyzzyaaae4(3)*xyzzyaaae4(3))
select case (periodicity)
case(2)
call ewald_2d(1,xyzzyaaae4,xyzzyaaad4,xyzzyaaaf4)
case(3)
call ewald_3d(1,xyzzyaaae4,xyzzyaaad4,.false.,xyzzyaaaf4)
case default
call errstop('COMPUTE_II_FIELDS','Invalid periodicity : only 2D and 3D&
&  fields currently supported.')
end select
ion_fields(1:3,xyzzyaaaa4)=ion_fields(1:3,xyzzyaaaa4)-xyzzyaaac4*xyzzy&
&aaaf4(1:3)
enddo
enddo
endif
end subroutine compute_ii_fields
end module slaarnaco

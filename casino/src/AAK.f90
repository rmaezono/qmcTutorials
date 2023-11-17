module slaarnaak
use dsp
use slaarnabg
use parallel
use store
use slaarnaag,    only : pi
use file_utils,   only : open_units,skip
use format_utils, only : i2s,wout
use slaarnabt,    only : interp_nev_2d_with_derivs,lookup
use run_control,  only : errstop,check_alloc
implicit none
private
public dwfdet_setup,dimer_orb_eval,get_dwfdet_orbmap,get_dwfdet_orbdes&
&c,get_dwfdet_ndesc
integer,parameter :: xyzzyaaaa1=7
integer xyzzyaaab1(2),xyzzyaaac1,xyzzyaaad1,xyzzyaaae1
integer xyzzyaaaf1
integer xyzzyaaag1
real(dp) xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1
integer,allocatable :: xyzzyaaak1(:,:,:,:)
integer,allocatable :: xyzzyaaal1(:,:,:)
real(dp),allocatable :: xyzzyaaam1(:,:,:)
real(dp),allocatable :: xyzzyaaan1(:),xyzzyaaao1(:)
logical,allocatable :: xyzzyaaap1(:)
integer,parameter :: xyzzyaaaq1=3,xyzzyaaar1=0
integer,allocatable :: xyzzyaaas1(:,:)
real(dp) xyzzyaaat1(xyzzyaaaa1,xyzzyaaaa1),xyzzyaaau1(xyzzyaaaa1),xyzz&
&yaaav1(xyzzyaaaa1)
real(dp),allocatable :: xyzzyaaaw1(:),xyzzyaaax1(:),xyzzyaaay1(:),xyzz&
&yaaaz1(:),xyzzyaaba1(:),xyzzyaabb1(:)
contains
subroutine dwfdet_setup(eionion)
implicit none
real(dp),intent(out) :: eionion
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaal2&
&,xyzzyaaam2
integer,allocatable :: xyzzyaaan2(:,:)
real(dp) xyzzyaaao2,xyzzyaaap2,xyzzyaaaq2
logical xyzzyaaar2
character(30) ctmp(2)
character(80) char_80
character(1024) title
if(am_master)then
call open_units(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaab2/=0)call errstop('DWFDET_SETUP','Error getting a unit for&
& dwfn.data.')
open(xyzzyaaaa2,file='dwfn.data',status='old',action='read',iostat=xyz&
&zyaaab2)
if(xyzzyaaab2/=0)call errstop('DWFDET_SETUP','Problem opening dwfn.dat&
&a.')
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,'(a)')title
title=trim(adjustl(title))
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*)xyzzyaaao2
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*)xyzzyaaap2,xyzzyaaaq2
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*)xyzzyaaaf1
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*)xyzzyaaab1(1),xyzzyaaab1(2),xyzzyaaac1
call skip(xyzzyaaaa2,1)
call skip(xyzzyaaaa2,xyzzyaaac1*sum(xyzzyaaab1(1:2)))
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*)xyzzyaaad1,xyzzyaaae1
nelec=xyzzyaaab1(1)+xyzzyaaab1(2)
close(xyzzyaaaa2)
endif
call mpi_bcast(xyzzyaaap2,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting z1 in dwfdet_setup')
call mpi_bcast(xyzzyaaaq2,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting z2 in dwfdet_setup')
call mpi_bcast(nelec,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nelec in dwfdet_setup')
call mpi_bcast(xyzzyaaaf1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dw_norb_rad in dwfdet_setup')
call mpi_bcast(xyzzyaaab1,2,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dwnele in dwfdet_setup')
call mpi_bcast(xyzzyaaac1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dw_ndet in dwfdet_setup')
call mpi_bcast(xyzzyaaad1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dw_nnu in dwfdet_setup')
call mpi_bcast(xyzzyaaae1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dw_nmu in dwfdet_setup')
allocate(xyzzyaaak1(3,max(xyzzyaaab1(1),xyzzyaaab1(2)),2,xyzzyaaac1),x&
&yzzyaaam1(xyzzyaaad1,xyzzyaaae1,xyzzyaaaf1),xyzzyaaao1(xyzzyaaae1),xy&
&zzyaaan1(xyzzyaaad1),xyzzyaaap1(xyzzyaaaf1),xyzzyaaal1(nemax,2,xyzzya&
&aac1),xyzzyaaas1(xyzzyaaaq1,xyzzyaaac1*nelec),stat=xyzzyaaah2)
call check_alloc(xyzzyaaah2,'DWFDET_SETUP','orbital data')
xyzzyaaag1=0
xyzzyaaal1=0
xyzzyaaas1=0
xyzzyaaak1=0
if(am_master)then
open(xyzzyaaaa2,file='dwfn.data',status='old',action='read',iostat=xyz&
&zyaaab2)
if(xyzzyaaab2/=0)call errstop('DWFDET_SETUP','Problem opening dwfn.dat&
&a.')
call skip(xyzzyaaaa2,11)
do xyzzyaaad2=1,xyzzyaaac1
do xyzzyaaag2=1,2
do xyzzyaaac2=1,xyzzyaaab1(xyzzyaaag2)
read(xyzzyaaaa2,*)xyzzyaaam2,xyzzyaaak1(1,xyzzyaaac2,xyzzyaaag2,xyzzya&
&aad2),xyzzyaaak1(2,xyzzyaaac2,xyzzyaaag2,xyzzyaaad2),xyzzyaaak1(3,xyz&
&zyaaac2,xyzzyaaag2,xyzzyaaad2)
enddo
enddo
enddo
call skip(xyzzyaaaa2,2)
read(xyzzyaaaa2,*)xyzzyaaah1,xyzzyaaai1
xyzzyaaam1(:,:,:)=0.d0
do xyzzyaaaf2=1,xyzzyaaaf1
read(xyzzyaaaa2,*,iostat=xyzzyaaab2)ctmp,xyzzyaaai2
if(xyzzyaaab2/=0)exit
if(xyzzyaaai2<1.or.xyzzyaaai2>xyzzyaaaf1)call errstop('READ_DWFN','Inv&
&alid orbital index - see dwfn.data file')
do xyzzyaaak2=1,xyzzyaaae1
do xyzzyaaal2=1,xyzzyaaad1
read(xyzzyaaaa2,*)xyzzyaaam1(xyzzyaaal2,xyzzyaaak2,xyzzyaaai2)
enddo
enddo
enddo
close(xyzzyaaaa2)
do xyzzyaaaj2=1,xyzzyaaae1
xyzzyaaao1(xyzzyaaaj2)=dble(xyzzyaaaj2-1)*xyzzyaaai1
enddo
do xyzzyaaaj2=1,xyzzyaaad1
xyzzyaaan1(xyzzyaaaj2)=dble((xyzzyaaaj2-1))*xyzzyaaah1
enddo
xyzzyaaaj1=xyzzyaaao2/2.d0
call wout('Dimer configuration:')
call wout()
call wout('Det Elec. Spin iorb |m| m   idx')
call wout('-------------------------------')
do xyzzyaaad2=1,xyzzyaaac1
do xyzzyaaag2=1,2
do xyzzyaaac2=1,xyzzyaaab1(xyzzyaaag2)
xyzzyaaar2=.false.
do xyzzyaaae2=1,xyzzyaaag1
if(all(xyzzyaaak1(1:3,xyzzyaaac2,xyzzyaaag2,xyzzyaaad2)==xyzzyaaas1(1:&
&3,xyzzyaaae2)))then
xyzzyaaar2=.true.
exit
endif
enddo
if(.not.xyzzyaaar2)then
xyzzyaaag1=xyzzyaaag1+1
xyzzyaaae2=xyzzyaaag1
xyzzyaaas1(1:3,xyzzyaaae2)=xyzzyaaak1(1:3,xyzzyaaac2,xyzzyaaag2,xyzzya&
&aad2)
endif
write(char_80,'(i3,2x,i3,3x,i2,3x,i2,3x,i2,1x,i2,1x,i4)')xyzzyaaad2,xy&
&zzyaaac2,xyzzyaaag2,xyzzyaaas1(1,xyzzyaaae2),xyzzyaaas1(2,xyzzyaaae2)&
&,xyzzyaaas1(3,xyzzyaaae2),xyzzyaaae2
call wout(char_80)
xyzzyaaal1(xyzzyaaac2,xyzzyaaag2,xyzzyaaad2)=xyzzyaaae2
enddo
enddo
enddo
call wout()
endif
call mpi_bcast(xyzzyaaah1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting dw_hnu in dwfdet_setup')
call mpi_bcast(xyzzyaaai1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting dw_hmu in dwfdet_setup')
call mpi_bcast(xyzzyaaak1,3*2*max(xyzzyaaab1(1),xyzzyaaab1(2))*xyzzyaa&
&ac1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dw_ecfg in dwfdet_setup')
call mpi_bcast(xyzzyaaam1,xyzzyaaad1*xyzzyaaae1*xyzzyaaaf1,mpi_double_&
&precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dw_orb in dwfdet_setup')
call mpi_bcast(xyzzyaaag1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dw_norb in dwfdet_setup')
call mpi_bcast(xyzzyaaan1,xyzzyaaad1,mpi_double_precision,0,mpi_comm_w&
&orld,ierror)
call checkmpi(ierror,'Broadcasting dw_vnu in dwfdet_setup')
call mpi_bcast(xyzzyaaao1,xyzzyaaae1,mpi_double_precision,0,mpi_comm_w&
&orld,ierror)
call checkmpi(ierror,'Broadcasting dw_vmu in dwfdet_setup')
call mpi_bcast(xyzzyaaaj1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting dw_h in dwfdet_setup')
call mpi_bcast(xyzzyaaal1,2*max(xyzzyaaab1(1),xyzzyaaab1(2))*xyzzyaaac&
&1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting dw_orbmap in dwfdet_setup')
if(am_master)then
allocate(xyzzyaaan2(xyzzyaaaq1,xyzzyaaag1),stat=xyzzyaaah2)
call check_alloc(xyzzyaaah2,'DWFDET_SETUP','itrf2')
xyzzyaaan2(:,1:xyzzyaaag1)=xyzzyaaas1(:,1:xyzzyaaag1)
deallocate(xyzzyaaas1)
allocate(xyzzyaaas1(3,xyzzyaaag1),stat=xyzzyaaah2)
call check_alloc(xyzzyaaah2,'DWFDET_SETUP','dw_orbdesc_int')
xyzzyaaas1(:,1:xyzzyaaag1)=xyzzyaaan2(:,1:xyzzyaaag1)
deallocate(xyzzyaaan2)
else
deallocate(xyzzyaaas1)
allocate(xyzzyaaas1(xyzzyaaaq1,xyzzyaaag1),stat=xyzzyaaah2)
call check_alloc(xyzzyaaah2,'DWFDET_SETUP','dw_orbdesc_int')
endif
call mpi_bcast(xyzzyaaas1,xyzzyaaaq1*xyzzyaaag1,mpi_integer,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting dw_orbdesc_int in dwfdet_setup')
nbasis=2
allocate(basis(3,nbasis),atno(nbasis))
atno(1)=nint(xyzzyaaap2)
atno(2)=nint(xyzzyaaaq2)
basis(1,1)=0.d0
basis(2,1)=0.d0
basis(3,1)=xyzzyaaaj1
basis(1,2)=0.d0
basis(2,2)=0.d0
basis(3,2)=-xyzzyaaaj1
pa1(1)=500.d0
pa1(2)=0.d0
pa1(3)=0.d0
pa2(1)=0.d0
pa2(2)=500.d0
pa2(3)=0.d0
pa3(1)=0.d0
pa3(2)=0.d0
pa3(3)=500.d0
if(am_master)then
call wout('Dimer formed by:')
call wout(' - Atom of atomic number '//trim(i2s(atno(1)))//' located a&
&t (0, 0, ',xyzzyaaaj1,')')
call wout(' - Atom of atomic number '//trim(i2s(atno(2)))//' located a&
&t (0, 0, ',-xyzzyaaaj1,')')
call wout()
endif
eionion=atno(1)*atno(2)/(2.d0*xyzzyaaaj1)
ndet=xyzzyaaac1
allocate(xyzzyaaaw1(xyzzyaaaf1),xyzzyaaax1(xyzzyaaaf1),xyzzyaaay1(xyzz&
&yaaaf1),xyzzyaaaz1(xyzzyaaaf1),xyzzyaaba1(xyzzyaaaf1),xyzzyaabb1(xyzz&
&yaaaf1),stat=xyzzyaaah2)
call check_alloc(xyzzyaaah2,'DWFDET_SETUP','wf*')
end subroutine dwfdet_setup
subroutine dimer_orb_eval(rvec,jspin,norb,orbmask,val,fsd,orbval,orbgr&
&ad,orblap,orbsderivs)
implicit none
integer,intent(in) :: jspin,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd,orbmask(norb)
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3,xyzzyaaag3,xyzzyaaah3
real(dp) xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3,xyzzyaaal3,xyzzyaaam3,xyzzya&
&aan3,xyzzyaaao3,xyzzyaaap3,xyzzyaaaq3,xyzzyaaar3,xyzzyaaas3,xyzzyaaat&
&3,xyzzyaaau3,xyzzyaaav3,xyzzyaaaw3,xyzzyaaax3(6),xyzzyaaay3,xyzzyaaaz&
&3,xyzzyaaba3,xyzzyaabb3,xyzzyaabc3,xyzzyaabd3(3,3),xyzzyaabe3(6,3)
logical xyzzyaabf3,xyzzyaabg3
integer xyzzyaabh3,xyzzyaabi3,xyzzyaabj3,xyzzyaabk3
real(dp) xyzzyaabl3,xyzzyaabm3,xyzzyaabn3
xyzzyaaal3=rvec(1)
xyzzyaaam3=rvec(2)
xyzzyaaan3=rvec(3)
xyzzyaaao3=sqrt(xyzzyaaal3**2+xyzzyaaam3**2+(xyzzyaaan3+xyzzyaaaj1)**2&
&)
xyzzyaaap3=sqrt(xyzzyaaal3**2+xyzzyaaam3**2+(xyzzyaaan3-xyzzyaaaj1)**2&
&)
xyzzyaaai3=0.5d0*(xyzzyaaao3+xyzzyaaap3)/xyzzyaaaj1
if(xyzzyaaai3<1.d0)xyzzyaaai3=xyzzyaaai3/abs(xyzzyaaai3)
xyzzyaaaj3=0.5d0*(xyzzyaaao3-xyzzyaaap3)/xyzzyaaaj1
if(abs(xyzzyaaaj3)>1.d0)xyzzyaaaj3=xyzzyaaaj3/abs(xyzzyaaaj3)
xyzzyaaak3=atan2(xyzzyaaam3,xyzzyaaal3)
xyzzyaabm3=acos(xyzzyaaaj3)
xyzzyaabn3=dlog(xyzzyaaai3+dsqrt(xyzzyaaai3*xyzzyaaai3-1.d0) )
xyzzyaaay3=xyzzyaaak3
xyzzyaaaq3=sin(xyzzyaabm3)
xyzzyaaar3=cos(xyzzyaabm3)
xyzzyaaas3=sinh(xyzzyaabn3)
xyzzyaaat3=cosh(xyzzyaabn3)
xyzzyaaau3=sin(xyzzyaaay3)
xyzzyaaav3=cos(xyzzyaaay3)
xyzzyaaaw3=tan(xyzzyaaay3)
xyzzyaabj3=0
xyzzyaabk3=0
do
if(xyzzyaabm3>-pi)exit
xyzzyaabm3=xyzzyaabm3+2.d0*pi
enddo
do
if(xyzzyaabm3<=pi)exit
xyzzyaabm3=xyzzyaabm3-2.d0*pi
enddo
do
if(xyzzyaabm3>=0.d0)exit
xyzzyaabm3=-xyzzyaabm3
xyzzyaabj3=xyzzyaabj3+1
enddo
do
if(xyzzyaabn3>=0.d0)exit
xyzzyaabn3=-xyzzyaabn3
xyzzyaabk3=xyzzyaabk3+1
enddo
if(xyzzyaabn3>xyzzyaaao1(xyzzyaaae1))then
orbval(1:xyzzyaaag1,1)=0.d0
orbgrad(:,1:xyzzyaaag1,1)=0.d0
orblap(1:xyzzyaaag1,1)=0.d0
if(present(orbsderivs))then
orbsderivs(:,1:xyzzyaaag1,1)=0.d0
endif
return
endif
call lookup(xyzzyaaan1,xyzzyaaad1,xyzzyaabm3,xyzzyaabi3)
call lookup(xyzzyaaao1,xyzzyaaae1,xyzzyaabn3,xyzzyaabh3)
xyzzyaabi3=xyzzyaabi3-(xyzzyaaaa1-1)/2
xyzzyaabh3=min(xyzzyaabh3-(xyzzyaaaa1-1)/2,xyzzyaaae1+1-xyzzyaaaa1)
if((xyzzyaabi3>=1).and.(xyzzyaaad1>=xyzzyaabi3+xyzzyaaaa1-1))then
xyzzyaaau1=xyzzyaaan1(xyzzyaabi3:xyzzyaabi3+xyzzyaaaa1-1)
else
do xyzzyaaaa3=1,xyzzyaaaa1
xyzzyaaau1(xyzzyaaaa3)=dble(xyzzyaabi3+xyzzyaaaa3-2)*xyzzyaaah1
enddo
endif
if((xyzzyaabh3>=1).and.(xyzzyaaae1>=xyzzyaabh3+xyzzyaaaa1-1))then
xyzzyaaav1=xyzzyaaao1(xyzzyaabh3:xyzzyaabh3+xyzzyaaaa1-1)
else
do xyzzyaaaa3=1,xyzzyaaaa1
xyzzyaaav1(xyzzyaaaa3)=dble(xyzzyaabh3+xyzzyaaaa3-2)*xyzzyaaai1
enddo
endif
xyzzyaaap1=.false.
do xyzzyaaab3=1,xyzzyaaag1
if(.not.orbmask(xyzzyaaab3))cycle
xyzzyaaac3=xyzzyaaas1(1,xyzzyaaab3)
xyzzyaaag3=xyzzyaaas1(2,xyzzyaaab3)
xyzzyaaah3=xyzzyaaas1(3,xyzzyaaab3)
xyzzyaabl3=(-1.d0)**xyzzyaaag3
xyzzyaabc3=dble(xyzzyaaah3)
if(.not.xyzzyaaap1(xyzzyaaac3))then
call xyzzyaabo3(xyzzyaaac3,xyzzyaaaw1(xyzzyaaac3),xyzzyaaax1(xyzzyaaac&
&3),xyzzyaaay1(xyzzyaaac3),xyzzyaaaz1(xyzzyaaac3),xyzzyaaba1(xyzzyaaac&
&3),xyzzyaabb1(xyzzyaaac3))
xyzzyaaap1(xyzzyaaac3)=.true.
endif
if(xyzzyaaah3>0)then
xyzzyaaaz3=cos(xyzzyaabc3*xyzzyaaay3)
xyzzyaaba3=-xyzzyaabc3*sin(xyzzyaabc3*xyzzyaaay3)
xyzzyaabb3=-xyzzyaabc3**2*cos(xyzzyaabc3*xyzzyaaay3)
elseif(xyzzyaaah3<0)then
xyzzyaaaz3=sin(xyzzyaabc3*xyzzyaaay3)
xyzzyaaba3=xyzzyaabc3*cos(xyzzyaabc3*xyzzyaaay3)
xyzzyaabb3=-xyzzyaabc3**2*sin(xyzzyaabc3*xyzzyaaay3)
else
xyzzyaaaz3=1.d0
xyzzyaaba3=0.d0
xyzzyaabb3=0.d0
endif
if(val)orbval(xyzzyaaab3,1)=xyzzyaaaw1(xyzzyaaac3)*xyzzyaaaz3
if(fsd)then
xyzzyaabf3=.false.
if(abs(xyzzyaabm3)<1.d-6.or.abs(xyzzyaabm3-pi)<1.d-6)xyzzyaabf3=.true.
xyzzyaabg3=.false.
if(abs(xyzzyaabn3)<1.d-6)xyzzyaabg3=.true.
orblap(xyzzyaaab3,1)=0.d0
orblap(xyzzyaaab3,1)=orblap(xyzzyaaab3,1)+xyzzyaaaz1(xyzzyaaac3)
orblap(xyzzyaaab3,1)=orblap(xyzzyaaab3,1)+xyzzyaaba1(xyzzyaaac3)
if(xyzzyaabf3)then
orblap(xyzzyaaab3,1)=orblap(xyzzyaaab3,1)+xyzzyaaaz1(xyzzyaaac3)-xyzzy&
&aabc3**2*xyzzyaaaz1(xyzzyaaac3)
else
orblap(xyzzyaaab3,1)=orblap(xyzzyaaab3,1)+xyzzyaaar3/xyzzyaaaq3*xyzzya&
&aax1(xyzzyaaac3)-xyzzyaabc3**2*(1.0/xyzzyaaaq3**2)*xyzzyaaaw1(xyzzyaa&
&ac3)
endif
if(xyzzyaabg3)then
orblap(xyzzyaaab3,1)=orblap(xyzzyaaab3,1)+xyzzyaaba1(xyzzyaaac3)-xyzzy&
&aabc3**2*xyzzyaaaz1(xyzzyaaac3)
else
orblap(xyzzyaaab3,1)=orblap(xyzzyaaab3,1)+xyzzyaaat3/xyzzyaaas3*xyzzya&
&aay1(xyzzyaaac3)-xyzzyaabc3**2*(1.0/xyzzyaaas3**2)*xyzzyaaaw1(xyzzyaa&
&ac3)
endif
orblap(xyzzyaaab3,1)=orblap(xyzzyaaab3,1)*xyzzyaaaz3
orblap(xyzzyaaab3,1)=orblap(xyzzyaaab3,1)/(xyzzyaaaq3**2+xyzzyaaas3**2&
&)/xyzzyaaaj1**2
if(xyzzyaabf3.and.xyzzyaabg3)orblap(xyzzyaaab3,1)=0.d0
xyzzyaabd3(1,1)=-1.d0/(2.d0*xyzzyaaaj1)/xyzzyaaaq3*(xyzzyaaal3/xyzzyaa&
&ao3-xyzzyaaal3/xyzzyaaap3)
xyzzyaabd3(2,1)=-1.d0/(2.d0*xyzzyaaaj1)/xyzzyaaaq3*(xyzzyaaam3/xyzzyaa&
&ao3-xyzzyaaam3/xyzzyaaap3)
xyzzyaabd3(3,1)=-1.d0/(2.d0*xyzzyaaaj1)/xyzzyaaaq3*((xyzzyaaan3+xyzzya&
&aaj1)/xyzzyaaao3-(xyzzyaaan3-xyzzyaaaj1)/xyzzyaaap3)
xyzzyaabd3(1,2)=1.d0/(2.d0*xyzzyaaaj1)/xyzzyaaas3*(xyzzyaaal3/xyzzyaaa&
&o3+xyzzyaaal3/xyzzyaaap3)
xyzzyaabd3(2,2)=1.d0/(2.d0*xyzzyaaaj1)/xyzzyaaas3*(xyzzyaaam3/xyzzyaaa&
&o3+xyzzyaaam3/xyzzyaaap3)
xyzzyaabd3(3,2)=1.d0/(2.d0*xyzzyaaaj1)/xyzzyaaas3*((xyzzyaaan3+xyzzyaa&
&aj1)/xyzzyaaao3+(xyzzyaaan3-xyzzyaaaj1)/xyzzyaaap3)
xyzzyaabd3(1,3)=-xyzzyaaau3/(xyzzyaaas3*xyzzyaaaq3)/xyzzyaaaj1
xyzzyaabd3(2,3)=xyzzyaaav3/(xyzzyaaas3*xyzzyaaaq3)/xyzzyaaaj1
xyzzyaabd3(3,3)=0.d0
orbgrad(1,xyzzyaaab3,1)=xyzzyaabd3(1,1)*xyzzyaaax1(xyzzyaaac3)*xyzzyaa&
&az3+xyzzyaabd3(1,2)*xyzzyaaay1(xyzzyaaac3)*xyzzyaaaz3+xyzzyaabd3(1,3)&
&*xyzzyaaaw1(xyzzyaaac3)*xyzzyaaba3
orbgrad(2,xyzzyaaab3,1)=xyzzyaabd3(2,1)*xyzzyaaax1(xyzzyaaac3)*xyzzyaa&
&az3+xyzzyaabd3(2,2)*xyzzyaaay1(xyzzyaaac3)*xyzzyaaaz3+xyzzyaabd3(2,3)&
&*xyzzyaaaw1(xyzzyaaac3)*xyzzyaaba3
orbgrad(3,xyzzyaaab3,1)=xyzzyaabd3(3,1)*xyzzyaaax1(xyzzyaaac3)*xyzzyaa&
&az3+xyzzyaabd3(3,2)*xyzzyaaay1(xyzzyaaac3)*xyzzyaaaz3+xyzzyaabd3(3,3)&
&*xyzzyaaaw1(xyzzyaaac3)*xyzzyaaba3
if(xyzzyaabf3.or.xyzzyaabg3)orbgrad(:,xyzzyaaab3,1)=0.d0
if(present(orbsderivs))then
xyzzyaabe3(1,1)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaal3**2/xyzzyaaao3**3+x&
&yzzyaaal3**2/xyzzyaaap3**3+1.d0/xyzzyaaao3-1.d0/xyzzyaaap3)
xyzzyaabe3(2,1)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaam3**2/xyzzyaaao3**3+x&
&yzzyaaam3**2/xyzzyaaap3**3+1.d0/xyzzyaaao3-1.d0/xyzzyaaap3)
xyzzyaabe3(3,1)=1.d0/(2.d0*xyzzyaaaj1)*(-(xyzzyaaan3+xyzzyaaaj1)**2/xy&
&zzyaaao3**3+(xyzzyaaan3-xyzzyaaaj1)**2/xyzzyaaap3**3+1.d0/xyzzyaaao3-&
&1.d0/xyzzyaaap3)
xyzzyaabe3(4,1)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaal3*xyzzyaaam3/xyzzyaa&
&ao3**3+xyzzyaaal3*xyzzyaaam3/xyzzyaaap3**3)
xyzzyaabe3(5,1)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaal3*(xyzzyaaan3+xyzzya&
&aaj1)/xyzzyaaao3**3+xyzzyaaal3*(xyzzyaaan3-xyzzyaaaj1)/xyzzyaaap3**3)
xyzzyaabe3(6,1)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaam3*(xyzzyaaan3+xyzzya&
&aaj1)/xyzzyaaao3**3+xyzzyaaam3*(xyzzyaaan3-xyzzyaaaj1)/xyzzyaaap3**3)
xyzzyaabe3(1,1)=-(xyzzyaabe3(1,1)+xyzzyaaar3*xyzzyaabd3(1,1)**2)/xyzzy&
&aaaq3
xyzzyaabe3(2,1)=-(xyzzyaabe3(2,1)+xyzzyaaar3*xyzzyaabd3(2,1)**2)/xyzzy&
&aaaq3
xyzzyaabe3(3,1)=-(xyzzyaabe3(3,1)+xyzzyaaar3*xyzzyaabd3(3,1)**2)/xyzzy&
&aaaq3
xyzzyaabe3(4,1)=-(xyzzyaabe3(4,1)+xyzzyaaar3*xyzzyaabd3(1,1)*xyzzyaabd&
&3(2,1))/xyzzyaaaq3
xyzzyaabe3(5,1)=-(xyzzyaabe3(5,1)+xyzzyaaar3*xyzzyaabd3(1,1)*xyzzyaabd&
&3(3,1))/xyzzyaaaq3
xyzzyaabe3(6,1)=-(xyzzyaabe3(6,1)+xyzzyaaar3*xyzzyaabd3(2,1)*xyzzyaabd&
&3(3,1))/xyzzyaaaq3
xyzzyaabe3(1,2)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaal3**2/xyzzyaaao3**3-x&
&yzzyaaal3**2/xyzzyaaap3**3+1.d0/xyzzyaaao3+1.d0/xyzzyaaap3)
xyzzyaabe3(2,2)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaam3**2/xyzzyaaao3**3-x&
&yzzyaaam3**2/xyzzyaaap3**3+1.d0/xyzzyaaao3+1.d0/xyzzyaaap3)
xyzzyaabe3(3,2)=1.d0/(2.d0*xyzzyaaaj1)*(-(xyzzyaaan3+xyzzyaaaj1)**2/xy&
&zzyaaao3**3-(xyzzyaaan3-xyzzyaaaj1)**2/xyzzyaaap3**3+1.d0/xyzzyaaao3+&
&1.d0/xyzzyaaap3)
xyzzyaabe3(4,2)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaal3*xyzzyaaam3/xyzzyaa&
&ao3**3-xyzzyaaal3*xyzzyaaam3/xyzzyaaap3**3)
xyzzyaabe3(5,2)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaal3*(xyzzyaaan3+xyzzya&
&aaj1)/xyzzyaaao3**3-xyzzyaaal3*(xyzzyaaan3-xyzzyaaaj1)/xyzzyaaap3**3)
xyzzyaabe3(6,2)=1.d0/(2.d0*xyzzyaaaj1)*(-xyzzyaaam3*(xyzzyaaan3+xyzzya&
&aaj1)/xyzzyaaao3**3-xyzzyaaam3*(xyzzyaaan3-xyzzyaaaj1)/xyzzyaaap3**3)
xyzzyaabe3(1,2)=(xyzzyaabe3(1,2)-xyzzyaaat3*xyzzyaabd3(1,2)**2)/xyzzya&
&aas3
xyzzyaabe3(2,2)=(xyzzyaabe3(2,2)-xyzzyaaat3*xyzzyaabd3(2,2)**2)/xyzzya&
&aas3
xyzzyaabe3(3,2)=(xyzzyaabe3(3,2)-xyzzyaaat3*xyzzyaabd3(3,2)**2)/xyzzya&
&aas3
xyzzyaabe3(4,2)=(xyzzyaabe3(4,2)-xyzzyaaat3*xyzzyaabd3(1,2)*xyzzyaabd3&
&(2,2))/xyzzyaaas3
xyzzyaabe3(5,2)=(xyzzyaabe3(5,2)-xyzzyaaat3*xyzzyaabd3(1,2)*xyzzyaabd3&
&(3,2))/xyzzyaaas3
xyzzyaabe3(6,2)=(xyzzyaabe3(6,2)-xyzzyaaat3*xyzzyaabd3(2,2)*xyzzyaabd3&
&(3,2))/xyzzyaaas3
xyzzyaabe3(1,3)=2.d0*xyzzyaaaw3*(1.d0/(xyzzyaaaj1*xyzzyaaas3*xyzzyaaaq&
&3)**2-xyzzyaabd3(1,3)**2)
xyzzyaabe3(2,3)=-2.d0*xyzzyaaaw3*xyzzyaabd3(2,3)**2
xyzzyaabe3(3,3)=0.d0
xyzzyaabe3(4,3)=-1.d0/(xyzzyaaaj1*xyzzyaaas3*xyzzyaaaq3)**2-2.d0*xyzzy&
&aaaw3*xyzzyaabd3(1,3)*xyzzyaabd3(2,3)
xyzzyaabe3(5,3)=0.d0
xyzzyaabe3(6,3)=0.d0
xyzzyaaax3(1)=xyzzyaaaz1(xyzzyaaac3)*xyzzyaaaz3
xyzzyaaax3(2)=xyzzyaaba1(xyzzyaaac3)*xyzzyaaaz3
xyzzyaaax3(3)=xyzzyaaaw1(xyzzyaaac3)*xyzzyaabb3
xyzzyaaax3(4)=xyzzyaabb1(xyzzyaaac3)*xyzzyaaaz3
xyzzyaaax3(5)=xyzzyaaax1(xyzzyaaac3)*xyzzyaaba3
xyzzyaaax3(6)=xyzzyaaay1(xyzzyaaac3)*xyzzyaaba3
do xyzzyaaaf3=1,6
select case(xyzzyaaaf3)
case(1)
xyzzyaaad3=1
xyzzyaaae3=1
case(2)
xyzzyaaad3=2
xyzzyaaae3=2
case(3)
xyzzyaaad3=3
xyzzyaaae3=3
case(4)
xyzzyaaad3=1
xyzzyaaae3=2
case(5)
xyzzyaaad3=1
xyzzyaaae3=3
case(6)
xyzzyaaad3=2
xyzzyaaae3=3
end select
orbsderivs(xyzzyaaaf3,xyzzyaaab3,1)=xyzzyaabe3(xyzzyaaaf3,1)*xyzzyaaax&
&1(xyzzyaaac3)*xyzzyaaaz3+xyzzyaabe3(xyzzyaaaf3,2)*xyzzyaaay1(xyzzyaaa&
&c3)*xyzzyaaaz3+xyzzyaabe3(xyzzyaaaf3,3)*xyzzyaaaw1(xyzzyaaac3)*xyzzya&
&aba3+xyzzyaabd3(xyzzyaaad3,1)*xyzzyaabd3(xyzzyaaae3,1)*xyzzyaaax3(1)+&
&xyzzyaabd3(xyzzyaaad3,2)*xyzzyaabd3(xyzzyaaae3,2)*xyzzyaaax3(2)+xyzzy&
&aabd3(xyzzyaaad3,3)*xyzzyaabd3(xyzzyaaae3,3)*xyzzyaaax3(3)+(xyzzyaabd&
&3(xyzzyaaad3,1)*xyzzyaabd3(xyzzyaaae3,2)+xyzzyaabd3(xyzzyaaad3,2)*xyz&
&zyaabd3(xyzzyaaae3,1))*xyzzyaaax3(4)+(xyzzyaabd3(xyzzyaaad3,1)*xyzzya&
&abd3(xyzzyaaae3,3)+xyzzyaabd3(xyzzyaaad3,3)*xyzzyaabd3(xyzzyaaae3,1))&
&*xyzzyaaax3(5)+(xyzzyaabd3(xyzzyaaad3,2)*xyzzyaabd3(xyzzyaaae3,3)+xyz&
&zyaabd3(xyzzyaaad3,3)*xyzzyaabd3(xyzzyaaae3,2))*xyzzyaaax3(6)
enddo
if(xyzzyaabf3.or.xyzzyaabg3)orbsderivs(:,xyzzyaaab3,1)=0.d0
endif
endif
enddo
if(complex_wf)then
if(val)orbval(:,2)=0.d0
if(fsd)then
orbgrad(:,:,2)=0.d0
orblap(:,2)=0.d0
if(present(orbsderivs))orbsderivs(:,:,2)=0.d0
endif
endif
contains
subroutine xyzzyaabo3(xyzzyaaac3,wf_t,wf_nu_t,wf_mu_t,wf_nu2_t,wf_mu2_&
&t,wf_numu_t)
implicit none
integer,intent(in) :: xyzzyaaac3
real(dp),intent(out) :: wf_t,wf_nu_t,wf_mu_t,wf_nu2_t,wf_mu2_t,wf_numu&
&_t
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4
xyzzyaaae4=0
xyzzyaaaf4=0
do xyzzyaaab4=1,xyzzyaaaa1
xyzzyaaad4=xyzzyaabh3+xyzzyaaab4-1
xyzzyaaae4=0
if(xyzzyaaad4<1)then
xyzzyaaad4=(2-xyzzyaaad4)
xyzzyaaae4=1
endif
do xyzzyaaaa4=1,xyzzyaaaa1
xyzzyaaac4=xyzzyaabi3+xyzzyaaaa4-1
xyzzyaaaf4=0
if(xyzzyaaac4<1)then
xyzzyaaac4=(2-xyzzyaaac4)
xyzzyaaaf4=1
endif
if(xyzzyaaac4>xyzzyaaad1)then
xyzzyaaac4=2*xyzzyaaad1-xyzzyaaac4
xyzzyaaaf4=xyzzyaaaf4+1
endif
xyzzyaaat1(xyzzyaaaa4,xyzzyaaab4)=xyzzyaabl3**(xyzzyaaae4+xyzzyaaaf4)*&
&xyzzyaaam1(xyzzyaaac4,xyzzyaaad4,xyzzyaaac3)
enddo
enddo
call interp_nev_2d_with_derivs(xyzzyaaau1,xyzzyaaav1,xyzzyaaat1,xyzzya&
&aaa1,        xyzzyaaaa1,xyzzyaabm3,xyzzyaabn3,wf_t,wf_nu_t,wf_mu_t,wf&
&_nu2_t,wf_mu2_t,wf_numu_t)
wf_t=wf_t*xyzzyaabl3**(xyzzyaabj3+xyzzyaabk3)
wf_nu_t=wf_nu_t*xyzzyaabl3**(xyzzyaabk3)
wf_mu_t=wf_mu_t*xyzzyaabl3**(xyzzyaabj3)
wf_nu2_t=wf_nu2_t*xyzzyaabl3**(xyzzyaabj3+xyzzyaabk3)
wf_mu2_t=wf_mu2_t*xyzzyaabl3**(xyzzyaabj3+xyzzyaabk3)
wf_numu_t=wf_numu_t
end subroutine xyzzyaabo3
end subroutine dimer_orb_eval
subroutine get_dwfdet_orbmap(row_offset,norb,orbmap)
implicit none
integer,intent(inout) :: row_offset(nspin),norb,orbmap(nemax,nspin,nde&
&t)
integer xyzzyaaaa5
do xyzzyaaaa5=1,nspin
orbmap(row_offset(xyzzyaaaa5)+1:row_offset(xyzzyaaaa5)+nuc_nele(xyzzya&
&aaa5),xyzzyaaaa5,:)=norb+xyzzyaaal1(1:nuc_nele(xyzzyaaaa5),xyzzyaaaa5&
&,:)
row_offset(xyzzyaaaa5)=row_offset(xyzzyaaaa5)+nuc_nele(xyzzyaaaa5)
enddo
norb=norb+xyzzyaaag1
end subroutine get_dwfdet_orbmap
subroutine get_dwfdet_ndesc(ndesc_int,ndesc_dp)
implicit none
integer,intent(inout) :: ndesc_int,ndesc_dp
ndesc_int=xyzzyaaaq1
ndesc_dp=xyzzyaaar1
end subroutine get_dwfdet_ndesc
subroutine get_dwfdet_orbdesc(norb,ndesc_int,ndesc_dp,orbdesc_int,orbd&
&esc_dp)
implicit none
integer,intent(in) :: norb,ndesc_int,ndesc_dp
integer,intent(inout) :: orbdesc_int(ndesc_int,norb)
real(dp),intent(inout) :: orbdesc_dp(ndesc_dp,norb)
orbdesc_int(1:xyzzyaaaq1,1:xyzzyaaag1)=xyzzyaaas1(1:xyzzyaaaq1,1:xyzzy&
&aaag1)
end subroutine get_dwfdet_orbdesc
end module slaarnaak

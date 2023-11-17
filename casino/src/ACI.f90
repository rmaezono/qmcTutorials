module slaarnaci
use parallel
use slaarnaag,   only : third
use dsp,         only : dp
use format_utils,only : wout,r2s
use slaarnabg,    only : pbinv,periodicity,b1,b2,b3
use run_control, only : errstop,check_alloc
use store,       only : nele,nemax,noncoll_spin,orb_norm,use_backflow,&
&ndet,nspin,k_offset
implicit none
private
public sdw_setup,print_sdwdata,sdw_orb_eval,sdw_magvec,sdw_kvec,sdw_or&
&b_theta,get_sdwfdet_orbmap,get_sdwfdet_orbdesc,get_sdwfdet_ndesc
integer xyzzyaaaa1(3)
integer,parameter :: xyzzyaaab1=2,xyzzyaaac1=1
integer,allocatable :: xyzzyaaad1(:,:)
real(dp) sdw_magvec(3),xyzzyaaae1,xyzzyaaaf1(3),xyzzyaaag1(3)
real(dp),allocatable :: sdw_kvec(:,:),sdw_orb_theta(:),xyzzyaaah1(:,:)&
&,xyzzyaaai1(:,:,:),xyzzyaaaj1(:,:),xyzzyaaak1(:,:,:)
complex(dp),allocatable :: xyzzyaaal1(:),xyzzyaaam1(:),xyzzyaaan1(:)
contains
subroutine sdw_setup
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2
real(dp) xyzzyaaad2,xyzzyaaae2(3)
noncoll_spin=.true.
xyzzyaaad2=1.d0
do xyzzyaaaa2=2,nemax
xyzzyaaad2=xyzzyaaad2*real(xyzzyaaaa2,dp)**(-1.d0/real(nemax,dp))
enddo
if(periodicity==2)xyzzyaaad2=xyzzyaaad2**2
xyzzyaaad2=orb_norm*xyzzyaaad2
xyzzyaaae1=xyzzyaaad2**third
allocate(xyzzyaaah1(nele(1),2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'SDW_SETUP','pre-calculation arrays')
do xyzzyaaaa2=1,nele(1)
xyzzyaaah1(xyzzyaaaa2,1)=cos(0.5d0*sdw_orb_theta(xyzzyaaaa2))
xyzzyaaah1(xyzzyaaaa2,2)=sin(0.5d0*sdw_orb_theta(xyzzyaaaa2))
enddo
xyzzyaaae2=sdw_magvec*0.5d0
xyzzyaaag1=k_offset+xyzzyaaae2
xyzzyaaaf1=k_offset-xyzzyaaae2
allocate(xyzzyaaad1(3,nele(1)),xyzzyaaai1(3,nele(1),2),xyzzyaaaj1(nele&
&(1),2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'SDW_SETUP','KINT array')
if(use_backflow)then
allocate(xyzzyaaak1(6,nele(1),2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'SDW_SETUP','K_PLUS_Q_2_PROD array')
endif
do xyzzyaaaa2=1,nele(1)
xyzzyaaad1(1:3,xyzzyaaaa2)=nint(matmul(sdw_kvec(1:3,xyzzyaaaa2)-k_offs&
&et,pbinv))
xyzzyaaai1(1:3,xyzzyaaaa2,1)=sdw_kvec(1:3,xyzzyaaaa2)-xyzzyaaae2
xyzzyaaai1(1:3,xyzzyaaaa2,2)=sdw_kvec(1:3,xyzzyaaaa2)+xyzzyaaae2
do xyzzyaaab2=1,2
xyzzyaaaj1(xyzzyaaaa2,xyzzyaaab2)=dot_product(xyzzyaaai1(1:3,xyzzyaaaa&
&2,xyzzyaaab2),xyzzyaaai1(1:3,xyzzyaaaa2,xyzzyaaab2))
if(use_backflow)then
xyzzyaaak1(1,xyzzyaaaa2,xyzzyaaab2)=-xyzzyaaai1(1,xyzzyaaaa2,xyzzyaaab&
&2)*xyzzyaaai1(1,xyzzyaaaa2,xyzzyaaab2)
xyzzyaaak1(2,xyzzyaaaa2,xyzzyaaab2)=-xyzzyaaai1(2,xyzzyaaaa2,xyzzyaaab&
&2)*xyzzyaaai1(2,xyzzyaaaa2,xyzzyaaab2)
xyzzyaaak1(3,xyzzyaaaa2,xyzzyaaab2)=-xyzzyaaai1(3,xyzzyaaaa2,xyzzyaaab&
&2)*xyzzyaaai1(3,xyzzyaaaa2,xyzzyaaab2)
xyzzyaaak1(4,xyzzyaaaa2,xyzzyaaab2)=-xyzzyaaai1(1,xyzzyaaaa2,xyzzyaaab&
&2)*xyzzyaaai1(2,xyzzyaaaa2,xyzzyaaab2)
xyzzyaaak1(5,xyzzyaaaa2,xyzzyaaab2)=-xyzzyaaai1(1,xyzzyaaaa2,xyzzyaaab&
&2)*xyzzyaaai1(3,xyzzyaaaa2,xyzzyaaab2)
xyzzyaaak1(6,xyzzyaaaa2,xyzzyaaab2)=-xyzzyaaai1(2,xyzzyaaaa2,xyzzyaaab&
&2)*xyzzyaaai1(3,xyzzyaaaa2,xyzzyaaab2)
endif
enddo
enddo
do xyzzyaaaa2=1,3
xyzzyaaaa1(xyzzyaaaa2)=maxval(abs(xyzzyaaad1(xyzzyaaaa2,:)))
enddo
allocate(xyzzyaaal1(-xyzzyaaaa1(1):xyzzyaaaa1(1)),xyzzyaaam1(-xyzzyaaa&
&a1(2):xyzzyaaaa1(2)),xyzzyaaan1(-xyzzyaaaa1(3):xyzzyaaaa1(3)),stat=xy&
&zzyaaac2)
call check_alloc(xyzzyaaac2,'SDW_SETUP','EXPIKDOTR_B arrays')
end subroutine sdw_setup
subroutine print_sdwdata()
implicit none
integer xyzzyaaaa3
real(dp) xyzzyaaab3
character(80) tmpr
xyzzyaaab3=0.d0
do xyzzyaaaa3=1,nele(1)
xyzzyaaab3=xyzzyaaab3+xyzzyaaah1(xyzzyaaaa3,1)*xyzzyaaah1(xyzzyaaaa3,1&
&)*dot_product(sdw_kvec(1:3,xyzzyaaaa3)-0.5d0*sdw_magvec,sdw_kvec(1:3,&
&xyzzyaaaa3)-0.5d0*sdw_magvec)*0.5d0
xyzzyaaab3=xyzzyaaab3+xyzzyaaah1(xyzzyaaaa3,2)*xyzzyaaah1(xyzzyaaaa3,2&
&)*dot_product(sdw_kvec(1:3,xyzzyaaaa3)+0.5d0*sdw_magvec,sdw_kvec(1:3,&
&xyzzyaaaa3)+0.5d0*sdw_magvec)*0.5d0
enddo
xyzzyaaab3=xyzzyaaab3/dble(nele(1))
tmpr=r2s(xyzzyaaab3,'(f20.12)')
call wout(' HF kinetic energy         : '//trim(tmpr))
end subroutine print_sdwdata
subroutine sdw_orb_eval(rvec,spin,lnorb,norb,orbmask,val,fsd,orbval,or&
&bgrad,orblap,orbsderivs)
implicit none
integer,intent(in) :: spin,lnorb,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(lnorb,*),orbgrad(3,lnorb,*),orblap(ln&
&orb,*)
real(dp),intent(inout),optional :: orbsderivs(6,lnorb,*)
logical,intent(in) :: val,fsd,orbmask(*)
integer xyzzyaaaa4,xyzzyaaab4
real(dp) xyzzyaaac4,xyzzyaaad4,xyzzyaaae4
complex(dp) xyzzyaaaf4,xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4
if(spin==1)then
xyzzyaaae4=dot_product(xyzzyaaaf1,rvec)
elseif(spin==2)then
xyzzyaaae4=dot_product(xyzzyaaag1,rvec)
else
call errstop('SDW_ORB_EVAL','Spin is neither up nor down. This shouldn&
&''t happen.')
endif
xyzzyaaaj4=cmplx(cos(xyzzyaaae4),sin(xyzzyaaae4),dp)
xyzzyaaae4=dot_product(b1,rvec)
xyzzyaaaf4=cmplx(cos(xyzzyaaae4),sin(xyzzyaaae4),dp)
xyzzyaaae4=dot_product(b2,rvec)
xyzzyaaag4=cmplx(cos(xyzzyaaae4),sin(xyzzyaaae4),dp)
xyzzyaaae4=dot_product(b3,rvec)
xyzzyaaah4=cmplx(cos(xyzzyaaae4),sin(xyzzyaaae4),dp)
xyzzyaaal1(0)=cmplx(xyzzyaaae1,0.d0,dp)
xyzzyaaam1(0)=xyzzyaaal1(0)
xyzzyaaan1(0)=xyzzyaaal1(0)
do xyzzyaaab4=1,xyzzyaaaa1(1)
xyzzyaaal1(xyzzyaaab4)=xyzzyaaal1(xyzzyaaab4-1)*xyzzyaaaf4
xyzzyaaal1(-xyzzyaaab4)=conjg(xyzzyaaal1(xyzzyaaab4))
enddo
do xyzzyaaab4=1,xyzzyaaaa1(2)
xyzzyaaam1(xyzzyaaab4)=xyzzyaaam1(xyzzyaaab4-1)*xyzzyaaag4
xyzzyaaam1(-xyzzyaaab4)=conjg(xyzzyaaam1(xyzzyaaab4))
enddo
do xyzzyaaab4=1,xyzzyaaaa1(3)
xyzzyaaan1(xyzzyaaab4)=xyzzyaaan1(xyzzyaaab4-1)*xyzzyaaah4
xyzzyaaan1(-xyzzyaaab4)=conjg(xyzzyaaan1(xyzzyaaab4))
enddo
do xyzzyaaab4=-xyzzyaaaa1(1),xyzzyaaaa1(1)
xyzzyaaal1(xyzzyaaab4)=xyzzyaaal1(xyzzyaaab4)*xyzzyaaaj4
enddo
do xyzzyaaaa4=1,norb
if(.not.orbmask(xyzzyaaaa4))cycle
xyzzyaaai4=xyzzyaaal1(xyzzyaaad1(1,xyzzyaaaa4))*xyzzyaaam1(xyzzyaaad1(&
&2,xyzzyaaaa4))*xyzzyaaan1(xyzzyaaad1(3,xyzzyaaaa4))
xyzzyaaac4=xyzzyaaah1(xyzzyaaaa4,spin)*dble(xyzzyaaai4)
xyzzyaaad4=xyzzyaaah1(xyzzyaaaa4,spin)*aimag(xyzzyaaai4)
if(val)then
orbval(xyzzyaaaa4,1)=xyzzyaaac4
orbval(xyzzyaaaa4,2)=xyzzyaaad4
endif
if(fsd)then
orbgrad(1:3,xyzzyaaaa4,1)=-xyzzyaaai1(1:3,xyzzyaaaa4,spin)*xyzzyaaad4
orbgrad(1:3,xyzzyaaaa4,2)=xyzzyaaai1(1:3,xyzzyaaaa4,spin)*xyzzyaaac4
orblap(xyzzyaaaa4,1)=-xyzzyaaaj1(xyzzyaaaa4,spin)*xyzzyaaac4
orblap(xyzzyaaaa4,2)=-xyzzyaaaj1(xyzzyaaaa4,spin)*xyzzyaaad4
if(present(orbsderivs))then
orbsderivs(1:6,xyzzyaaaa4,1)=xyzzyaaak1(1:6,xyzzyaaaa4,spin)*xyzzyaaac&
&4
orbsderivs(1:6,xyzzyaaaa4,2)=xyzzyaaak1(1:6,xyzzyaaaa4,spin)*xyzzyaaad&
&4
endif
endif
enddo
end subroutine sdw_orb_eval
subroutine get_sdwfdet_orbmap(row_offset,norb,orbmap)
implicit none
integer,intent(inout) :: row_offset(nspin),norb,orbmap(nemax,nspin,nde&
&t)
integer xyzzyaaaa5
do xyzzyaaaa5=1,nele(1)
orbmap(row_offset(1)+xyzzyaaaa5,1,1)=norb+xyzzyaaaa5
enddo
row_offset(1)=row_offset(1)+nele(1)
norb=norb+nele(1)
end subroutine get_sdwfdet_orbmap
subroutine get_sdwfdet_ndesc(ndesc_int,ndesc_dp)
implicit none
integer,intent(inout) :: ndesc_int,ndesc_dp
ndesc_int=xyzzyaaab1
ndesc_dp=xyzzyaaac1
end subroutine get_sdwfdet_ndesc
subroutine get_sdwfdet_orbdesc(norb,ndesc_int,ndesc_dp,orbdesc_int,orb&
&desc_dp)
implicit none
integer,intent(in) :: norb,ndesc_int,ndesc_dp
integer,intent(inout) :: orbdesc_int(ndesc_int,norb)
real(dp),intent(inout) :: orbdesc_dp(ndesc_dp,norb)
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7,xyzzyaa&
&af7
logical xyzzyaaag7
xyzzyaaac7=0
xyzzyaaae7=0
do xyzzyaaaa7=1,nele(1)
xyzzyaaag7=.false.
do xyzzyaaab7=1,xyzzyaaaa7-1
if(all(sdw_kvec(1:3,xyzzyaaaa7)==sdw_kvec(1:3,xyzzyaaab7)))then
xyzzyaaag7=.true.
exit
endif
enddo
if(.not.xyzzyaaag7)then
xyzzyaaac7=xyzzyaaac7+1
xyzzyaaad7=xyzzyaaac7
orbdesc_int(1,xyzzyaaaa7)=xyzzyaaad7
else
orbdesc_int(1,xyzzyaaaa7)=orbdesc_int(1,xyzzyaaab7)
endif
xyzzyaaag7=.false.
do xyzzyaaab7=1,xyzzyaaaa7-1
if(sdw_orb_theta(xyzzyaaaa7)==sdw_orb_theta(xyzzyaaab7))then
xyzzyaaag7=.true.
exit
endif
enddo
if(.not.xyzzyaaag7)then
xyzzyaaae7=xyzzyaaae7+1
xyzzyaaaf7=xyzzyaaae7
orbdesc_int(2,xyzzyaaaa7)=xyzzyaaaf7
orbdesc_dp(1,xyzzyaaaa7)=sdw_orb_theta(xyzzyaaaa7)
else
orbdesc_int(2,xyzzyaaaa7)=orbdesc_int(2,xyzzyaaab7)
orbdesc_dp(1,xyzzyaaaa7)=orbdesc_dp(1,xyzzyaaab7)
endif
enddo
end subroutine get_sdwfdet_orbdesc
end module slaarnaci

module slaarnabs
use dsp,           only : dp
use format_utils,  only : wout,i2s
use slaarnabg,      only : nitot,nitype,iontype,nitot_forces,which_ion&
&_displaced
use slaarnaca,         only : have_ppots,rcut_non_loc,have_veep
use slaarnacc,only : ranx
use run_control,   only : errstop,timer,check_alloc
use store,         only : forces
implicit none
private
public :: init_non_local,setup_non_local_grid,print_non_local_grid,nl_&
&refgrid,nl_am_proj,v_non_local,tmove_no_points,tmove_points,tmove_t,t&
&move_fixgrid,nl_rcut,nl_wfrat_forces,dvelf_forces,nl_nrefgrid,nl_refg&
&ridw,nl_wrkgrid_all,tmove_t_moved,calc_nl_projection,v_non_local_tmov&
&e
real(dp),allocatable :: nl_rcut(:)
integer,allocatable :: nl_nrefgrid(:)
integer,allocatable :: xyzzyaaaa1(:)
integer,allocatable :: xyzzyaaab1(:)
real(dp),allocatable :: nl_refgrid(:,:,:)
real(dp),allocatable :: nl_refgridw(:,:)
real(dp),allocatable :: xyzzyaaac1(:)
real(dp),allocatable :: xyzzyaaad1(:,:),xyzzyaaae1(:)
real(dp),allocatable :: nl_wrkgrid_all(:,:,:,:),nl_wfrat_forces(:,:,:)&
&,dvelf_forces(:,:,:,:)
real(dp),allocatable :: xyzzyaaaf1(:,:,:)
real(dp),allocatable :: nl_am_proj(:,:)
real(dp),allocatable :: xyzzyaaag1(:)
integer,allocatable :: xyzzyaaah1(:)
integer,allocatable :: tmove_no_points(:,:)
real(dp),allocatable :: tmove_t(:,:,:),tmove_t_moved(:,:,:),xyzzyaaai1&
&(:,:,:),tmove_points(:,:,:,:)
logical,parameter :: tmove_fixgrid=.false.
contains
subroutine init_non_local(irule)
use slaarnaag, only : fourpi
use slaarnaca,     only : nlang,maxl_cpp
use store,     only : use_tmove,netot,forces
implicit none
integer,intent(in) :: irule(nitype)
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2
integer,parameter :: xyzzyaaag2(7)=(/1,4,6,12,18,26,50/)
xyzzyaaac2=maxval(nlang)
allocate(xyzzyaaah1(nitype),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'INIT_NON_LOCAL','0')
if(have_veep)then
xyzzyaaah1(:)=max(nlang(:),maxl_cpp)
else
xyzzyaaah1=nlang
endif
xyzzyaaad2=maxval(xyzzyaaah1)
xyzzyaaae2=maxval(irule)
xyzzyaaaf2=xyzzyaaag2(xyzzyaaae2)
allocate(xyzzyaaad1(3,xyzzyaaaf2),xyzzyaaae1(xyzzyaaaf2),xyzzyaaac1(xy&
&zzyaaaf2),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'INIT_NON_LOCAL','0')
allocate(nl_rcut(nitot),nl_am_proj(0:xyzzyaaad2,nitot),nl_nrefgrid(nit&
&ype),xyzzyaaaa1(nitype),xyzzyaaab1(nitype),nl_refgrid(3,xyzzyaaaf2,ni&
&type),nl_refgridw(xyzzyaaaf2,nitype),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'INIT_NON_LOCAL','1')
do xyzzyaaaa2=1,nitot
nl_rcut(xyzzyaaaa2)=rcut_non_loc(iontype(xyzzyaaaa2))
enddo
nl_am_proj=0.d0
nl_refgrid=0.d0
nl_refgridw=0.d0
xyzzyaaad1=0.d0
xyzzyaaac1=1.d0
xyzzyaaaa1=-1
allocate(xyzzyaaag1(0:xyzzyaaad2),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'INIT_NON_LOCAL','3')
do xyzzyaaaa2=0,xyzzyaaad2
xyzzyaaag1(xyzzyaaaa2)=fourpi*xyzzyaaal1(1.d0,xyzzyaaaa2)
enddo
call setup_non_local_grid(irule)
xyzzyaaaa2=maxval(nl_nrefgrid)
if(xyzzyaaaa2/=xyzzyaaaf2)call errstop('INIT_NON_LOCAL','Failed consis&
&tency check in number of points in spherical grid for non-local integ&
&ration. This is a bug.')
if(use_tmove)then
allocate(tmove_no_points(nitot,netot),tmove_t(xyzzyaaaf2,nitot,netot),&
&tmove_t_moved(xyzzyaaaf2,nitot,netot),xyzzyaaai1(xyzzyaaaf2,0:xyzzyaa&
&ac2,nitot),tmove_points(3,xyzzyaaaf2,nitot,netot),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'INIT_NON_LOCAL','4')
endif
if(forces)then
allocate(nl_wfrat_forces(xyzzyaaaf2,nitot,netot),dvelf_forces(3,xyzzya&
&aaf2,nitot,netot),nl_wrkgrid_all(3,xyzzyaaaf2,nitot,netot),xyzzyaaaf1&
&(xyzzyaaaf2,nitot,netot),stat=xyzzyaaab2)
call check_alloc(xyzzyaaab2,'INIT_NON_LOCAL','2')
endif
end subroutine init_non_local
subroutine setup_non_local_grid(irule)
use slaarnaag, only : one_over_root_two,twentyoneth,fifteenth,thirtiet&
&h,twelfth,sixth,eleventh,third
use store, only : use_tmove
implicit none
integer,intent(in) :: irule(nitype)
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3
real(dp) xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3,xyzzya&
&aal3,xyzzyaaam3,xyzzyaaan3,xyzzyaaao3,xyzzyaaap3(3,3),xyzzyaaaq3,xyzz&
&yaaar3,xyzzyaaas3,xyzzyaaat3,xyzzyaaau3
real(dp),parameter ::  xyzzyaaav3=4.d0/105.d0,xyzzyaaaw3=27.d0/840.d0,&
&xyzzyaaax3=4.d0/315.d0,xyzzyaaay3=64.d0/2835.d0,xyzzyaaaz3=27.d0/1280&
&.d0,xyzzyaaba3=14641.d0/725760.d0
xyzzyaaag3=sqrt(third)
xyzzyaaah3=sqrt(eleventh)
xyzzyaaai3=3.d0*xyzzyaaah3
xyzzyaaaj3=sqrt((5.d0-sqrt(5.d0))*0.1d0)
xyzzyaaak3=sqrt((5.d0+sqrt(5.d0))*0.1d0)
do xyzzyaaaa3=1,nitype
if(irule(xyzzyaaaa3)==xyzzyaaaa1(xyzzyaaaa3))cycle
select case (irule(xyzzyaaaa3))
case(1)
nl_nrefgrid(xyzzyaaaa3)=1
xyzzyaaab1(xyzzyaaaa3)=0
nl_refgrid(1,1,xyzzyaaaa3)=0
nl_refgrid(2,1,xyzzyaaaa3)=0
nl_refgrid(3,1,xyzzyaaaa3)=1.d0
nl_refgridw(1,xyzzyaaaa3)=1.d0
case(2)
nl_nrefgrid(xyzzyaaaa3)=4
xyzzyaaab1(xyzzyaaaa3)=2
nl_refgrid(1,1,xyzzyaaaa3)=xyzzyaaag3
nl_refgrid(2,1,xyzzyaaaa3)=xyzzyaaag3
nl_refgrid(3,1,xyzzyaaaa3)=xyzzyaaag3
nl_refgrid(1,2,xyzzyaaaa3)=xyzzyaaag3
nl_refgrid(2,2,xyzzyaaaa3)=-xyzzyaaag3
nl_refgrid(3,2,xyzzyaaaa3)=-xyzzyaaag3
nl_refgrid(1,3,xyzzyaaaa3)=-xyzzyaaag3
nl_refgrid(2,3,xyzzyaaaa3)=xyzzyaaag3
nl_refgrid(3,3,xyzzyaaaa3)=-xyzzyaaag3
nl_refgrid(1,4,xyzzyaaaa3)=-xyzzyaaag3
nl_refgrid(2,4,xyzzyaaaa3)=-xyzzyaaag3
nl_refgrid(3,4,xyzzyaaaa3)=xyzzyaaag3
nl_refgridw(1:4,xyzzyaaaa3)=0.25d0
case(3)
nl_nrefgrid(xyzzyaaaa3)=6
xyzzyaaab1(xyzzyaaaa3)=3
nl_refgrid(1,1,xyzzyaaaa3)=1.d0
nl_refgrid(2,1,xyzzyaaaa3)=0.d0
nl_refgrid(3,1,xyzzyaaaa3)=0.d0
nl_refgrid(1,3,xyzzyaaaa3)=0.d0
nl_refgrid(2,3,xyzzyaaaa3)=1.d0
nl_refgrid(3,3,xyzzyaaaa3)=0.d0
nl_refgrid(1,5,xyzzyaaaa3)=0.d0
nl_refgrid(2,5,xyzzyaaaa3)=0.d0
nl_refgrid(3,5,xyzzyaaaa3)=1.d0
do xyzzyaaab3=1,5,2
nl_refgrid(1:3,xyzzyaaab3+1,xyzzyaaaa3)=-nl_refgrid(1:3,xyzzyaaab3,xyz&
&zyaaaa3)
enddo
nl_refgridw(1:6,xyzzyaaaa3)=sixth
case(4)
nl_nrefgrid(xyzzyaaaa3)=12
xyzzyaaab1(xyzzyaaaa3)=5
nl_refgrid(1,1,xyzzyaaaa3)=0.d0
nl_refgrid(2,1,xyzzyaaaa3)=xyzzyaaaj3
nl_refgrid(3,1,xyzzyaaaa3)=xyzzyaaak3
nl_refgrid(1,3,xyzzyaaaa3)=0.d0
nl_refgrid(2,3,xyzzyaaaa3)=-xyzzyaaaj3
nl_refgrid(3,3,xyzzyaaaa3)=xyzzyaaak3
nl_refgrid(1,5,xyzzyaaaa3)=-xyzzyaaak3
nl_refgrid(2,5,xyzzyaaaa3)=0.d0
nl_refgrid(3,5,xyzzyaaaa3)=xyzzyaaaj3
nl_refgrid(1,7,xyzzyaaaa3)=xyzzyaaak3
nl_refgrid(2,7,xyzzyaaaa3)=0.d0
nl_refgrid(3,7,xyzzyaaaa3)=xyzzyaaaj3
nl_refgrid(1,9,xyzzyaaaa3)=-xyzzyaaaj3
nl_refgrid(2,9,xyzzyaaaa3)=xyzzyaaak3
nl_refgrid(3,9,xyzzyaaaa3)=0.d0
nl_refgrid(1,11,xyzzyaaaa3)=xyzzyaaaj3
nl_refgrid(2,11,xyzzyaaaa3)=xyzzyaaak3
nl_refgrid(3,11,xyzzyaaaa3)=0.d0
do xyzzyaaab3=1,11,2
nl_refgrid(1:3,xyzzyaaab3+1,xyzzyaaaa3)=-nl_refgrid(1:3,xyzzyaaab3,xyz&
&zyaaaa3)
enddo
nl_refgridw(1:12,xyzzyaaaa3)=twelfth
case(5)
nl_nrefgrid(xyzzyaaaa3)=18
xyzzyaaab1(xyzzyaaaa3)=5
nl_refgrid(1,1,xyzzyaaaa3)=1.d0
nl_refgrid(2,1,xyzzyaaaa3)=0.d0
nl_refgrid(3,1,xyzzyaaaa3)=0.d0
nl_refgrid(1,3,xyzzyaaaa3)=0.d0
nl_refgrid(2,3,xyzzyaaaa3)=1.d0
nl_refgrid(3,3,xyzzyaaaa3)=0.d0
nl_refgrid(1,5,xyzzyaaaa3)=0.d0
nl_refgrid(2,5,xyzzyaaaa3)=0.d0
nl_refgrid(3,5,xyzzyaaaa3)=1.d0
do xyzzyaaab3=1,5,2
nl_refgrid(1:3,xyzzyaaab3+1,xyzzyaaaa3)=-nl_refgrid(1:3,xyzzyaaab3,xyz&
&zyaaaa3)
enddo
nl_refgridw(1:6,xyzzyaaaa3)=thirtieth
do xyzzyaaab3=7,15,4
if(xyzzyaaab3==7)then
xyzzyaaad3=1
xyzzyaaae3=2
xyzzyaaaf3=3
elseif(xyzzyaaab3==11)then
xyzzyaaad3=2
xyzzyaaae3=3
xyzzyaaaf3=1
else
xyzzyaaad3=3
xyzzyaaae3=1
xyzzyaaaf3=2
endif
nl_refgrid(xyzzyaaad3,xyzzyaaab3,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+1,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+1,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+1,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+2,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+2,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+2,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+3,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+3,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+3,xyzzyaaaa3)=0.d0
enddo
nl_refgridw(7:18,xyzzyaaaa3)=fifteenth
case(6)
nl_nrefgrid(xyzzyaaaa3)=26
xyzzyaaab1(xyzzyaaaa3)=7
nl_refgrid(1,1,xyzzyaaaa3)=1.d0
nl_refgrid(2,1,xyzzyaaaa3)=0.d0
nl_refgrid(3,1,xyzzyaaaa3)=0.d0
nl_refgrid(1,3,xyzzyaaaa3)=0.d0
nl_refgrid(2,3,xyzzyaaaa3)=1.d0
nl_refgrid(3,3,xyzzyaaaa3)=0.d0
nl_refgrid(1,5,xyzzyaaaa3)=0.d0
nl_refgrid(2,5,xyzzyaaaa3)=0.d0
nl_refgrid(3,5,xyzzyaaaa3)=1.d0
do xyzzyaaab3=1,5,2
nl_refgrid(1:3,xyzzyaaab3+1,xyzzyaaaa3)=-nl_refgrid(1:3,xyzzyaaab3,xyz&
&zyaaaa3)
enddo
nl_refgridw(1:6,xyzzyaaaa3)=twentyoneth
do xyzzyaaab3=7,15,4
if(xyzzyaaab3==7)then
xyzzyaaad3=1
xyzzyaaae3=2
xyzzyaaaf3=3
elseif(xyzzyaaab3==11)then
xyzzyaaad3=2
xyzzyaaae3=3
xyzzyaaaf3=1
else
xyzzyaaad3=3
xyzzyaaae3=1
xyzzyaaaf3=2
endif
nl_refgrid(xyzzyaaad3,xyzzyaaab3,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+1,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+1,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+1,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+2,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+2,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+2,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+3,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+3,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+3,xyzzyaaaa3)=0.d0
enddo
nl_refgridw(7:18,xyzzyaaaa3)=xyzzyaaav3
do xyzzyaaab3=0,7
do xyzzyaaac3=1,3
if(btest(xyzzyaaab3,xyzzyaaac3-1))then
nl_refgrid(xyzzyaaac3,19+xyzzyaaab3,xyzzyaaaa3)=-xyzzyaaag3
else
nl_refgrid(xyzzyaaac3,19+xyzzyaaab3,xyzzyaaaa3)=xyzzyaaag3
endif
enddo
enddo
nl_refgridw(19:26,xyzzyaaaa3)=xyzzyaaaw3
case(7)
nl_nrefgrid(xyzzyaaaa3)=50
xyzzyaaab1(xyzzyaaaa3)=11
nl_refgrid(1,1,xyzzyaaaa3)=1.d0
nl_refgrid(2,1,xyzzyaaaa3)=0.d0
nl_refgrid(3,1,xyzzyaaaa3)=0.d0
nl_refgrid(1,3,xyzzyaaaa3)=0.d0
nl_refgrid(2,3,xyzzyaaaa3)=1.d0
nl_refgrid(3,3,xyzzyaaaa3)=0.d0
nl_refgrid(1,5,xyzzyaaaa3)=0.d0
nl_refgrid(2,5,xyzzyaaaa3)=0.d0
nl_refgrid(3,5,xyzzyaaaa3)=1.d0
do xyzzyaaab3=1,5,2
nl_refgrid(1:3,xyzzyaaab3+1,xyzzyaaaa3)=-nl_refgrid(1:3,xyzzyaaab3,xyz&
&zyaaaa3)
enddo
nl_refgridw(1:6,xyzzyaaaa3)=xyzzyaaax3
do xyzzyaaab3=7,15,4
if(xyzzyaaab3==7)then
xyzzyaaad3=1
xyzzyaaae3=2
xyzzyaaaf3=3
elseif(xyzzyaaab3==11)then
xyzzyaaad3=2
xyzzyaaae3=3
xyzzyaaaf3=1
else
xyzzyaaad3=3
xyzzyaaae3=1
xyzzyaaaf3=2
endif
nl_refgrid(xyzzyaaad3,xyzzyaaab3,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+1,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+1,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+1,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+2,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+2,xyzzyaaaa3)=one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+2,xyzzyaaaa3)=0.d0
nl_refgrid(xyzzyaaad3,xyzzyaaab3+3,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaae3,xyzzyaaab3+3,xyzzyaaaa3)=-one_over_root_two
nl_refgrid(xyzzyaaaf3,xyzzyaaab3+3,xyzzyaaaa3)=0.d0
enddo
nl_refgridw(7:18,xyzzyaaaa3)=xyzzyaaay3
nl_refgrid(1:3,19:26,xyzzyaaaa3)=xyzzyaaag3
do xyzzyaaab3=0,7
do xyzzyaaac3=1,3
if(btest(xyzzyaaab3,xyzzyaaac3-1))nl_refgrid(xyzzyaaac3,19+xyzzyaaab3,&
&xyzzyaaaa3)=-nl_refgrid(xyzzyaaac3,19+xyzzyaaab3,xyzzyaaaa3)
enddo
enddo
nl_refgridw(19:26,xyzzyaaaa3)=xyzzyaaaz3
nl_refgrid(1,27:34,xyzzyaaaa3)=xyzzyaaah3
nl_refgrid(2,27:34,xyzzyaaaa3)=xyzzyaaah3
nl_refgrid(3,27:34,xyzzyaaaa3)=xyzzyaaai3
nl_refgrid(1,35:42,xyzzyaaaa3)=xyzzyaaah3
nl_refgrid(2,35:42,xyzzyaaaa3)=xyzzyaaai3
nl_refgrid(3,35:42,xyzzyaaaa3)=xyzzyaaah3
nl_refgrid(1,43:50,xyzzyaaaa3)=xyzzyaaai3
nl_refgrid(2,43:50,xyzzyaaaa3)=xyzzyaaah3
nl_refgrid(3,43:50,xyzzyaaaa3)=xyzzyaaah3
do xyzzyaaab3=0,7
do xyzzyaaac3=1,3
if(btest(xyzzyaaab3,xyzzyaaac3-1))then
nl_refgrid(xyzzyaaac3,27+xyzzyaaab3,xyzzyaaaa3)=-nl_refgrid(xyzzyaaac3&
&,27+xyzzyaaab3,xyzzyaaaa3)
nl_refgrid(xyzzyaaac3,35+xyzzyaaab3,xyzzyaaaa3)=-nl_refgrid(xyzzyaaac3&
&,35+xyzzyaaab3,xyzzyaaaa3)
nl_refgrid(xyzzyaaac3,43+xyzzyaaab3,xyzzyaaaa3)=-nl_refgrid(xyzzyaaac3&
&,43+xyzzyaaab3,xyzzyaaaa3)
endif
enddo
enddo
nl_refgridw(27:50,xyzzyaaaa3)=xyzzyaaba3
case default
call errstop('SETUP_NON_LOCAL_GRID','Unsupported grid code: ' //trim(i&
&2s(irule(xyzzyaaaa3)))//'.')
end select
xyzzyaaaa1(xyzzyaaaa3)=irule(xyzzyaaaa3)
if(use_tmove.and.tmove_fixgrid)then
xyzzyaaam3=nl_refgrid(1,1,xyzzyaaaa3)
xyzzyaaan3=nl_refgrid(2,1,xyzzyaaaa3)
if(xyzzyaaam3/=0.d0.or.xyzzyaaan3/=0.d0)then
xyzzyaaao3=nl_refgrid(3,1,xyzzyaaaa3)
xyzzyaaal3=sqrt(xyzzyaaam3*xyzzyaaam3+xyzzyaaan3*xyzzyaaan3+xyzzyaaao3&
&*xyzzyaaao3)
xyzzyaaau3=sqrt(xyzzyaaam3*xyzzyaaam3+xyzzyaaan3*xyzzyaaan3)
xyzzyaaas3=xyzzyaaao3/xyzzyaaal3
xyzzyaaaq3=xyzzyaaau3/xyzzyaaal3
if(xyzzyaaau3>0.d0)then
xyzzyaaat3=xyzzyaaam3/xyzzyaaau3
xyzzyaaar3=xyzzyaaan3/xyzzyaaau3
else
xyzzyaaat3=1.d0
xyzzyaaar3=0.d0
endif
xyzzyaaap3(1,1)=xyzzyaaas3*xyzzyaaat3
xyzzyaaap3(2,1)=-xyzzyaaar3
xyzzyaaap3(3,1)=xyzzyaaaq3*xyzzyaaat3
xyzzyaaap3(1,2)=xyzzyaaas3*xyzzyaaar3
xyzzyaaap3(2,2)=xyzzyaaat3
xyzzyaaap3(3,2)=xyzzyaaaq3*xyzzyaaar3
xyzzyaaap3(1,3)=-xyzzyaaaq3
xyzzyaaap3(2,3)=0.d0
xyzzyaaap3(3,3)=xyzzyaaas3
do xyzzyaaab3=1,nl_nrefgrid(xyzzyaaaa3)
nl_refgrid(1:3,xyzzyaaab3,xyzzyaaaa3)=matmul(xyzzyaaap3,nl_refgrid(1:3&
&,xyzzyaaab3,xyzzyaaaa3))
enddo
endif
endif
enddo
end subroutine setup_non_local_grid
subroutine print_non_local_grid
implicit none
integer xyzzyaaaa4,xyzzyaaab4
real(dp) xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaaaf4
logical,parameter :: xyzzyaaag4=.false.
character(80) tmpr
if(.not.have_ppots)return
call wout('Non-local integration grids')
call wout('===========================')
do xyzzyaaab4=1,nitype
call wout('Ion type            :  '//trim(i2s(xyzzyaaab4)))
call wout('Non-local grid no.  :  '//trim(i2s(xyzzyaaaa1(xyzzyaaab4)))&
&)
call wout('Lexact              :  '//trim(i2s(xyzzyaaab1(xyzzyaaab4)))&
&)
call wout('Number of points    :  '//trim(i2s(nl_nrefgrid(xyzzyaaab4))&
&))
call wout()
if(xyzzyaaag4)then
call wout('pt.   x          y          z          weight    r')
call wout('-----------------------------------------------------------&
&')
do xyzzyaaaa4=1,nl_nrefgrid(xyzzyaaab4)
xyzzyaaac4=nl_refgrid(1,xyzzyaaaa4,xyzzyaaab4)
xyzzyaaad4=nl_refgrid(2,xyzzyaaaa4,xyzzyaaab4)
xyzzyaaae4=nl_refgrid(3,xyzzyaaaa4,xyzzyaaab4)
xyzzyaaaf4=sqrt(xyzzyaaac4*xyzzyaaac4+xyzzyaaad4*xyzzyaaad4+xyzzyaaae4&
&*xyzzyaaae4)
write(tmpr,'(i2,1x,5(f10.6,1x))')xyzzyaaaa4,xyzzyaaac4,xyzzyaaad4,xyzz&
&yaaae4,nl_refgridw(xyzzyaaaa4,xyzzyaaab4),xyzzyaaaf4
call wout(tmpr)
enddo
call wout()
endif
enddo
end subroutine print_non_local_grid
subroutine calc_nl_projection(ii,is,js,grid_angles)
use slaarnabt,  only : dcopy
use slaarnach,    only : get_rsele,buffer_move_ion_from,get_eivecs,eiv&
&ecs_scr,rele_scr,sele_scr
use store,      only : use_tmove
use slaarnacs,  only : define_config_oneelec,wfn_ratio,wfn_loggrad
implicit none
integer,intent(in) :: ii,is,js
real(dp),intent(inout) :: grid_angles(3,nitot)
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5
real(dp) xyzzyaaag5,xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5,xyzzyaaak5(3),xyz&
&zyaaal5(3)
complex(dp) xyzzyaaam5,xyzzyaaan5(3)
logical isnan,isinf
call timer('CALC_NL_PROJECTION',.true.)
call get_rsele(is)
xyzzyaaaf5=is
if(buffer_move_ion_from(is)/=0)xyzzyaaaf5=buffer_move_ion_from(is)
call get_eivecs(xyzzyaaaf5)
do xyzzyaaaa5=1,nitot
xyzzyaaag5=eivecs_scr(4,xyzzyaaaa5,ii,xyzzyaaaf5)
if(use_tmove.and.which_ion_displaced==0)tmove_no_points(xyzzyaaaa5,ii)&
&=-1
if(xyzzyaaag5<=nl_rcut(xyzzyaaaa5))then
xyzzyaaad5=iontype(xyzzyaaaa5)
if(which_ion_displaced==0)then
call timer('CONSTRUCT_GRID',.true.,collapse=.true.)
xyzzyaaah5=eivecs_scr(1,xyzzyaaaa5,ii,xyzzyaaaf5)
xyzzyaaai5=eivecs_scr(2,xyzzyaaaa5,ii,xyzzyaaaf5)
xyzzyaaaj5=eivecs_scr(3,xyzzyaaaa5,ii,xyzzyaaaf5)
xyzzyaaak5(1)=rele_scr(1,ii,is)-xyzzyaaah5
xyzzyaaak5(2)=rele_scr(2,ii,is)-xyzzyaaai5
xyzzyaaak5(3)=rele_scr(3,ii,is)-xyzzyaaaj5
call xyzzyaaaj1(xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5,xyzzyaaag5,xyzzyaaak5&
&,xyzzyaaad5,grid_angles(:,xyzzyaaaa5))
if(use_tmove)then
tmove_no_points(xyzzyaaaa5,ii)=nl_nrefgrid(xyzzyaaad5)
call dcopy(3*tmove_no_points(xyzzyaaaa5,ii),xyzzyaaad1(1,1),1,tmove_po&
&ints(1,1,xyzzyaaaa5,ii),1)
endif
call timer('CONSTRUCT_GRID',.false.)
endif
call timer('WF_RATIO',.true.)
if(use_tmove.and.tmove_fixgrid)then
xyzzyaaae1(1)=1.d0
xyzzyaaae5=2
else
xyzzyaaae5=1
endif
do xyzzyaaab5=xyzzyaaae5,nl_nrefgrid(xyzzyaaad5)
if(which_ion_displaced==0)then
xyzzyaaal5(1:3)=xyzzyaaad1(1:3,xyzzyaaab5)
else
xyzzyaaal5(1:3)=nl_wrkgrid_all(1:3,xyzzyaaab5,xyzzyaaaa5,ii)
xyzzyaaac1(xyzzyaaab5)=xyzzyaaaf1(xyzzyaaab5,xyzzyaaaa5,ii)
endif
call define_config_oneelec(ii,is,js,xyzzyaaal5,sele_scr(ii,is))
call wfn_ratio(is,js,0,ratio=xyzzyaaam5,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)xyzzyaaam5=cmplx(0.d0,0.d0,dp)
if(forces.and.which_ion_displaced==0.and.xyzzyaaaa5<=nitot_forces)then
call wfn_loggrad(ii,js,0,xyzzyaaan5)
if(isnan.or.isinf)then
xyzzyaaam5=cmplx(0.d0,0.d0,dp)
xyzzyaaan5=cmplx(0.d0,0.d0,dp)
endif
dvelf_forces(1:3,xyzzyaaab5,xyzzyaaaa5,ii)=dble(xyzzyaaan5(1:3))*dble(&
&xyzzyaaam5)
nl_wfrat_forces(xyzzyaaab5,xyzzyaaaa5,ii)=dble(xyzzyaaam5)
nl_wrkgrid_all(1:3,xyzzyaaab5,xyzzyaaaa5,ii)=xyzzyaaal5(1:3)
xyzzyaaaf1(xyzzyaaab5,xyzzyaaaa5,ii)=xyzzyaaac1(xyzzyaaab5)
endif
xyzzyaaae1(xyzzyaaab5)=dble(xyzzyaaam5)
enddo
call timer('WF_RATIO',.false.)
call timer('PROJECTION',.true.)
do xyzzyaaac5=0,xyzzyaaah1(xyzzyaaad5)
call xyzzyaaak1(xyzzyaaac5,xyzzyaaad5,xyzzyaaaa5,nl_am_proj(xyzzyaaac5&
&,xyzzyaaaa5))
enddo
call timer('PROJECTION',.false.)
endif
enddo
call timer('CALC_NL_PROJECTION',.false.)
end subroutine calc_nl_projection
subroutine xyzzyaaaj1(rx,ry,rz,r,ionvec,itype,grid_angles)
use slaarnaag, only : pi,twopi
use store,     only : use_tmove
implicit none
integer,intent(in) :: itype
real(dp),intent(in) :: r,rx,ry,rz,ionvec(3)
real(dp),intent(inout) :: grid_angles(3)
integer xyzzyaaaa6
real(dp) xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaaaf6,xyzzya&
&aag6,xyzzyaaah6,xyzzyaaai6,xyzzyaaaj6,xyzzyaaak6,xyzzyaaal6(3,3),xyzz&
&yaaam6(3),xyzzyaaan6(3),xyzzyaaao6
xyzzyaaao6=1.d0/r
if(use_tmove.and.tmove_fixgrid)then
xyzzyaaak6=sqrt(rx*rx+ry*ry)
xyzzyaaaj6=rz*xyzzyaaao6
xyzzyaaai6=xyzzyaaak6*xyzzyaaao6
if(xyzzyaaak6>0.d0)then
xyzzyaaaf6=rx/xyzzyaaak6
xyzzyaaae6=ry/xyzzyaaak6
else
xyzzyaaaf6=1.d0
xyzzyaaae6=0.d0
endif
xyzzyaaac6=twopi*ranx()
xyzzyaaah6=cos(xyzzyaaac6)
xyzzyaaag6=sin(xyzzyaaac6)
xyzzyaaal6(1,1)=xyzzyaaaj6*xyzzyaaaf6*xyzzyaaah6-xyzzyaaae6*xyzzyaaag6
xyzzyaaal6(2,1)=xyzzyaaaj6*xyzzyaaae6*xyzzyaaah6+xyzzyaaaf6*xyzzyaaag6
xyzzyaaal6(3,1)=-xyzzyaaai6*xyzzyaaah6
xyzzyaaal6(1,2)=-xyzzyaaaj6*xyzzyaaaf6*xyzzyaaag6-xyzzyaaae6*xyzzyaaah&
&6
xyzzyaaal6(2,2)=-xyzzyaaaj6*xyzzyaaae6*xyzzyaaag6+xyzzyaaaf6*xyzzyaaah&
&6
xyzzyaaal6(3,2)=xyzzyaaai6*xyzzyaaag6
xyzzyaaal6(1,3)=xyzzyaaai6*xyzzyaaaf6
xyzzyaaal6(2,3)=xyzzyaaai6*xyzzyaaae6
xyzzyaaal6(3,3)=xyzzyaaaj6
if(grid_angles(1)<0.d0)then
if(rx>=0.d0)then
grid_angles(1)=asin(xyzzyaaae6)
else
grid_angles(1)=pi-asin(xyzzyaaae6)
endif
grid_angles(2)=xyzzyaaac6
grid_angles(3)=acos(xyzzyaaaj6)
else
call errstop('CONSTRUCT_INTEGRATION_GRID','Cannot supply angles if a f&
&ixed grid is used for T moves.')
endif
else
if(grid_angles(1)<0.d0)then
xyzzyaaab6=twopi*ranx()
xyzzyaaac6=twopi*ranx()
xyzzyaaad6=acos(1.d0-2.d0*ranx())
grid_angles(1)=xyzzyaaab6
grid_angles(2)=xyzzyaaac6
grid_angles(3)=xyzzyaaad6
else
xyzzyaaab6=grid_angles(1)
xyzzyaaac6=grid_angles(2)
xyzzyaaad6=grid_angles(3)
endif
xyzzyaaae6=sin(xyzzyaaab6)
xyzzyaaag6=sin(xyzzyaaac6)
xyzzyaaai6=sin(xyzzyaaad6)
xyzzyaaaf6=cos(xyzzyaaab6)
xyzzyaaah6=cos(xyzzyaaac6)
xyzzyaaaj6=cos(xyzzyaaad6)
xyzzyaaal6(1,1)=xyzzyaaah6*xyzzyaaaf6-xyzzyaaaj6*xyzzyaaae6*xyzzyaaag6
xyzzyaaal6(2,1)=-xyzzyaaag6*xyzzyaaaf6-xyzzyaaaj6*xyzzyaaae6*xyzzyaaah&
&6
xyzzyaaal6(3,1)=xyzzyaaai6*xyzzyaaae6
xyzzyaaal6(1,2)=xyzzyaaah6*xyzzyaaae6+xyzzyaaaj6*xyzzyaaaf6*xyzzyaaag6
xyzzyaaal6(2,2)=-xyzzyaaag6*xyzzyaaae6+xyzzyaaaj6*xyzzyaaaf6*xyzzyaaah&
&6
xyzzyaaal6(3,2)=-xyzzyaaai6*xyzzyaaaf6
xyzzyaaal6(1,3)=xyzzyaaag6*xyzzyaaai6
xyzzyaaal6(2,3)=xyzzyaaah6*xyzzyaaai6
xyzzyaaal6(3,3)=xyzzyaaaj6
endif
xyzzyaaam6(1)=rx*xyzzyaaao6
xyzzyaaam6(2)=ry*xyzzyaaao6
xyzzyaaam6(3)=rz*xyzzyaaao6
do xyzzyaaaa6=1,nl_nrefgrid(itype)
xyzzyaaan6(1)=xyzzyaaal6(1,1)*nl_refgrid(1,xyzzyaaaa6,itype)+xyzzyaaal&
&6(1,2)*nl_refgrid(2,xyzzyaaaa6,itype)+xyzzyaaal6(1,3)*nl_refgrid(3,xy&
&zzyaaaa6,itype)
xyzzyaaan6(2)=xyzzyaaal6(2,1)*nl_refgrid(1,xyzzyaaaa6,itype)+xyzzyaaal&
&6(2,2)*nl_refgrid(2,xyzzyaaaa6,itype)+xyzzyaaal6(2,3)*nl_refgrid(3,xy&
&zzyaaaa6,itype)
xyzzyaaan6(3)=xyzzyaaal6(3,1)*nl_refgrid(1,xyzzyaaaa6,itype)+xyzzyaaal&
&6(3,2)*nl_refgrid(2,xyzzyaaaa6,itype)+xyzzyaaal6(3,3)*nl_refgrid(3,xy&
&zzyaaaa6,itype)
xyzzyaaad1(1,xyzzyaaaa6)=r*xyzzyaaan6(1)+ionvec(1)
xyzzyaaad1(2,xyzzyaaaa6)=r*xyzzyaaan6(2)+ionvec(2)
xyzzyaaad1(3,xyzzyaaaa6)=r*xyzzyaaan6(3)+ionvec(3)
xyzzyaaac1(xyzzyaaaa6)=xyzzyaaam6(1)*xyzzyaaan6(1)+xyzzyaaam6(2)*xyzzy&
&aaan6(2)+xyzzyaaam6(3)*xyzzyaaan6(3)
enddo
end subroutine xyzzyaaaj1
subroutine xyzzyaaak1(il,itype,ion,proj)
use slaarnaca,only : nlang
use store,only : use_tmove
implicit none
integer,intent(in) :: il,itype,ion
real(dp),intent(out) :: proj
integer xyzzyaaaa7
if(use_tmove.and.il<=nlang(itype))then
proj=0.d0
do xyzzyaaaa7=1,nl_nrefgrid(itype)
xyzzyaaai1(xyzzyaaaa7,il,ion)=xyzzyaaag1(il)*xyzzyaaal1(xyzzyaaac1(xyz&
&zyaaaa7),il)*xyzzyaaae1(xyzzyaaaa7) *nl_refgridw(xyzzyaaaa7,itype)
proj=proj+xyzzyaaai1(xyzzyaaaa7,il,ion)
enddo
else
proj=0.d0
do xyzzyaaaa7=1,nl_nrefgrid(itype)
proj=proj+xyzzyaaal1(xyzzyaaac1(xyzzyaaaa7),il)*xyzzyaaae1(xyzzyaaaa7)&
&*nl_refgridw(xyzzyaaaa7,itype)
enddo
proj=proj*xyzzyaaag1(il)
endif
end subroutine xyzzyaaak1
real(dp) function xyzzyaaal1(cos_theta,l)
use slaarnaag,only : one_over_root_fourpi,root_threeoverfourpi,root_45&
&_over_16pi,root_175_over_16pi,third,threefifths,root_11025_over_256pi&
&,sixsevenths,threethirtyfifths,root_43659_over_256pi,tenninths,fivetw&
&entyoneths
implicit none
integer,intent(in) :: l
real(dp),intent(in) :: cos_theta
real(dp) xyzzyaaaa8
select case(l)
case(0)
xyzzyaaal1=one_over_root_fourpi
case(1)
xyzzyaaal1=root_threeoverfourpi*cos_theta
case(2)
xyzzyaaal1=root_45_over_16pi*(cos_theta*cos_theta-third)
case(3)
xyzzyaaal1=root_175_over_16pi*(cos_theta*cos_theta-threefifths)*cos_th&
&eta
case(4)
xyzzyaaaa8=cos_theta*cos_theta
xyzzyaaal1=root_11025_over_256pi*((xyzzyaaaa8-sixsevenths)*xyzzyaaaa8+&
&threethirtyfifths)
case(5)
xyzzyaaaa8=cos_theta*cos_theta
xyzzyaaal1=root_43659_over_256pi*((xyzzyaaaa8-tenninths)*xyzzyaaaa8+fi&
&vetwentyoneths)*cos_theta
case default
call errstop('YL0','Spherical harmonic not implemented for l='//trim(i&
&2s(l))//'.')
end select
end function xyzzyaaal1
subroutine v_non_local(ii,eivecs,pnl_psi,isnan,isinf)
use slaarnabt, only : interp_nev,lookup,daxpy
use slaarnaca,     only : ncoeff,npoly,pp_radial_grid,nlang,l_of_non_l&
&ocal,fcoeff_non_loc
use store,     only : netot,use_tmove
implicit none
integer,intent(in) :: ii
real(dp),intent(in) :: eivecs(4,nitot,netot)
real(dp),intent(out) :: pnl_psi
logical,intent(out) :: isnan,isinf
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,xyzzyaaae9
real(dp) xyzzyaaaf9,xyzzyaaag9,xyzzyaaah9,xyzzyaaai9
pnl_psi=0.d0
isnan=.false.
isinf=.false.
do xyzzyaaac9=1,nitot
xyzzyaaae9=iontype(xyzzyaaac9)
if(ncoeff(xyzzyaaae9)/=0)then
xyzzyaaaf9=eivecs(4,xyzzyaaac9,ii)
if(use_tmove)then
if(which_ion_displaced/=0)then
tmove_t_moved(:,xyzzyaaac9,ii)=0.d0
else
tmove_t(:,xyzzyaaac9,ii)=0.d0
endif
endif
if(xyzzyaaaf9<=rcut_non_loc(xyzzyaaae9))then
call lookup(pp_radial_grid(1,xyzzyaaae9),ncoeff(xyzzyaaae9),xyzzyaaaf9&
&,xyzzyaaab9)
xyzzyaaab9=min(max(xyzzyaaab9-(npoly-1)/2,1),ncoeff(xyzzyaaae9)+1-npol&
&y)
xyzzyaaai9=0.d0
do xyzzyaaaa9=1,nlang(xyzzyaaae9)
xyzzyaaad9=l_of_non_local(xyzzyaaaa9,xyzzyaaae9)
call interp_nev(pp_radial_grid(xyzzyaaab9,xyzzyaaae9),fcoeff_non_loc(x&
&yzzyaaab9,xyzzyaaad9,xyzzyaaae9),5,xyzzyaaaf9,xyzzyaaah9,xyzzyaaag9)
xyzzyaaai9=xyzzyaaai9+nl_am_proj(xyzzyaaad9,xyzzyaaac9)*xyzzyaaah9
if(use_tmove)then
if(which_ion_displaced/=0)then
call daxpy(tmove_no_points(xyzzyaaac9,ii),xyzzyaaah9,xyzzyaaai1(1,xyzz&
&yaaad9,xyzzyaaac9),1,tmove_t_moved(1,xyzzyaaac9,ii),1)
else
call daxpy(tmove_no_points(xyzzyaaac9,ii),xyzzyaaah9,xyzzyaaai1(1,xyzz&
&yaaad9,xyzzyaaac9),1,tmove_t(1,xyzzyaaac9,ii),1)
endif
endif
enddo
pnl_psi=pnl_psi+xyzzyaaai9
endif
endif
enddo
end subroutine v_non_local
subroutine v_non_local_tmove(poti_nl_plus,poti_nl_minus)
use store, only : netot
implicit none
real(dp),intent(out) :: poti_nl_plus,poti_nl_minus
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10
poti_nl_plus=0.d0
poti_nl_minus=0.d0
do xyzzyaaaa10=1,netot
if(maxval(tmove_no_points(1:nitot,xyzzyaaaa10))<1)cycle
do xyzzyaaab10=1,nitot
if(tmove_no_points(xyzzyaaab10,xyzzyaaaa10)<1)cycle
do xyzzyaaac10=1,tmove_no_points(xyzzyaaab10,xyzzyaaaa10)
if(tmove_t(xyzzyaaac10,xyzzyaaab10,xyzzyaaaa10)>=0.d0)then
if(which_ion_displaced/=0)then
poti_nl_plus=poti_nl_plus+tmove_t_moved(xyzzyaaac10,xyzzyaaab10,xyzzya&
&aaa10)
else
poti_nl_plus=poti_nl_plus+tmove_t(xyzzyaaac10,xyzzyaaab10,xyzzyaaaa10)
endif
else
if(which_ion_displaced/=0)then
poti_nl_minus=poti_nl_minus+tmove_t_moved(xyzzyaaac10,xyzzyaaab10,xyzz&
&yaaaa10)
else
poti_nl_minus=poti_nl_minus+tmove_t(xyzzyaaac10,xyzzyaaab10,xyzzyaaaa1&
&0)
endif
endif
enddo
enddo
enddo
end subroutine v_non_local_tmove
end module slaarnabs

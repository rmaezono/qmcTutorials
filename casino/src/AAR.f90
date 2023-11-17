module slaarnaar
use slaarnaag
use dsp
use store
use file_utils,  only : open_units,skip
use format_utils,only : wout,i2s,r2s
use slaarnabg,    only : a1,a2,a3,periodicity
use parallel,    only : am_master
use run_control, only : errstop,check_alloc
implicit none
private
public init_expot_wfdet,expot_orb_eval,get_exwfdet_orbmap,get_exwfdet_&
&orbdesc,get_exwfdet_ndesc
character(80) title
integer xyzzyaaaa1
integer xyzzyaaab1,xyzzyaaac1(2)
integer,allocatable :: xyzzyaaad1(:),xyzzyaaae1(:),xyzzyaaaf1(:)
real(dp),allocatable :: xyzzyaaag1(:,:),xyzzyaaah1(:,:)
logical,allocatable :: xyzzyaaai1(:)
integer,parameter :: xyzzyaaaj1=2
integer,parameter :: xyzzyaaak1=10000
integer,allocatable :: xyzzyaaal1(:),xyzzyaaam1(:),xyzzyaaan1(:),xyzzy&
&aaao1(:,:),xyzzyaaap1(:)
real(dp),allocatable :: xyzzyaaaq1(:),xyzzyaaar1(:,:),xyzzyaaas1(:,:,:&
&),xyzzyaaat1(:,:,:),xyzzyaaau1(:,:,:),xyzzyaaav1(:,:),xyzzyaaaw1(:,:)&
&,xyzzyaaax1(:),xyzzyaaay1(:),xyzzyaaaz1(:),xyzzyaaba1(:),xyzzyaabb1(:&
&),xyzzyaabc1(:)
integer,allocatable :: xyzzyaabd1(:),xyzzyaabe1(:,:),xyzzyaabf1(:,:),x&
&yzzyaabg1(:,:)
real(dp),allocatable :: xyzzyaabh1(:,:),xyzzyaabi1(:,:),xyzzyaabj1(:,:&
&)
integer,parameter :: xyzzyaabk1=0,xyzzyaabl1=0
contains
subroutine init_expot_wfdet
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaal2&
&,xyzzyaaam2,xyzzyaaan2,xyzzyaaao2,xyzzyaaap2
real(dp) xyzzyaaaq2(3)
logical xyzzyaaar2
character(80) char_80,dumchar
if(.not.any(heg_orbtype>100))return
inquire(file='expot.data',exist=xyzzyaaar2)
if(am_master.and..not.(use_expot.and.xyzzyaaar2))call errstop('EXPOT_W&
&FDET_SETUP','To use EXPOT orbitals expot.data must be present and EXP&
&OT set to T in input.')
xyzzyaaac1=0
xyzzyaaaa1=2
xyzzyaaaj2=0
xyzzyaaak2=0
xyzzyaaal2=-1
if(am_master)then
call wout('External-potential orbitals')
call wout('===========================')
call wout('Reading external-potential orbitals from expot.data file.')
call wout()
endif
call open_units(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaab2/=0)call errstop('EXPOT_WFDET_SETUP','Unable to find free&
& i/o unit.')
open(unit=xyzzyaaaa2,file='expot.data',status='old',iostat=xyzzyaaab2)
if(xyzzyaaab2/=0)call errstop('EXPOT_WFDET_SETUP','Problem opening exp&
&ot.data.')
do
read(xyzzyaaaa2,'(a)',iostat=xyzzyaaab2)char_80
if(xyzzyaaab2/=0)exit
if(trim(adjustl(char_80))=='START VERSION')then
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaaa1
read(xyzzyaaaa2,'(a)',end=10,err=20)char_80
if(am_master.and.trim(adjustl(char_80))/='END VERSION')call errstop('E&
&XPOT_WFDET_SETUP',"Expected to find 'END VERSION' in expot.data .")
exit
endif
enddo
rewind(xyzzyaaaa2)
do
read(xyzzyaaaa2,'(a)',iostat=xyzzyaaab2)char_80
if(xyzzyaaab2/=0)call errstop('EXPOT_WFDET_SETUP',"Label 'START EXPOT_&
&WFN' not found in expot.data .")
if(trim(adjustl(char_80))=='START EXPOT_WFN')exit
enddo
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,'(a)',end=10,err=20)title
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaab1
if(am_master.and.any(heg_orbtype>100+xyzzyaaab1))call errstop('EXPOT_W&
&FDET_SETUP','Orbital assigned in free_particles block not found in ex&
&pot.data .')
if(am_master)then
call wout('Title: '//trim(title))
call wout('Number of sets : '//trim(i2s(xyzzyaaab1)))
call wout()
endif
allocate(xyzzyaaai1(xyzzyaaab1),xyzzyaaad1(xyzzyaaab1),xyzzyaaae1(xyzz&
&yaaab1),xyzzyaaaf1(xyzzyaaab1),xyzzyaaag1(3,xyzzyaaab1),xyzzyaaah1(3,&
&xyzzyaaab1),xyzzyaabc1(xyzzyaaab1),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'EXPOT_WFDET_SETUP','0')
xyzzyaaad1=0
xyzzyaaai1=.false.
xyzzyaaaf1=0
xyzzyaaag1=0.d0
xyzzyaaae1=0
xyzzyaabc1=orb_norm
do xyzzyaaad2=1,xyzzyaaab1
read(xyzzyaaaa2,'(a)',end=10,err=20)char_80
if(am_master.and.trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzya&
&aad2)))call errstop('EXPOT_WFDET_SETUP',"Expected to find 'START SET &
&"//trim(i2s(xyzzyaaad2))//"' in expot.data .")
if(.not.any(heg_orbtype==100+xyzzyaaad2))then
do
read(xyzzyaaaa2,'(a)',iostat=xyzzyaaab2)char_80
if(xyzzyaaab2/=0)call errstop('EXPOT_WFDET_SETUP',"Label 'END SET "//t&
&rim(i2s(xyzzyaaad2))//"' not found in expot.data .")
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaad2)))exit
enddo
cycle
endif
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,'(a)',end=10,err=20)char_80
select case(trim(adjustl(char_80)))
case('APERIODIC')
xyzzyaaai1(xyzzyaaad2)=.false.
case('PERIODIC')
xyzzyaaai1(xyzzyaaad2)=.true.
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP','Invalid periodicity lab&
&el; must be either PERIODIC or APERIODIC.')
end select
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,'(a)',end=10,err=20)char_80
if(.not.xyzzyaaai1(xyzzyaaad2))then
select case(trim(adjustl(char_80)))
case('GAUSSIAN')
xyzzyaaaf2=2
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP',"Aperiodic orbital type &
&'"//trim(adjustl(char_80))//"' in expot.data not recognized.")
end select
else
select case(trim(adjustl(char_80)))
case('FOURIER')
xyzzyaaaf2=1
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP',"Periodic orbital type '&
&"//trim(adjustl(char_80))//"' in expot.data not recognized.")
end select
endif
xyzzyaaad1(xyzzyaaad2)=xyzzyaaaf2
xyzzyaaac1(xyzzyaaaf2)=xyzzyaaac1(xyzzyaaaf2)+1
xyzzyaaae1(xyzzyaaad2)=xyzzyaaac1(xyzzyaaaf2)
call skip(xyzzyaaaa2,5)
select case(xyzzyaaaf2)
case(1)
call skip(xyzzyaaaa2,3)
read(xyzzyaaaa2,*,end=10,err=20)char_80
select case(trim(adjustl(char_80)))
case('ODD','odd')
xyzzyaaag2=1
case('EVEN','even')
xyzzyaaag2=2
case('NONE','none')
xyzzyaaag2=0
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP','Unrecognized symmetry o&
&f Fourier series. Must be one of EVEN, ODD or NONE.')
end select
if(xyzzyaaal2==-1.or.xyzzyaaal2==xyzzyaaag2)then
xyzzyaaal2=xyzzyaaag2
else
xyzzyaaal2=0
endif
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaah2
xyzzyaaaj2=max(xyzzyaaaj2,xyzzyaaah2)
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaai2
xyzzyaaak2=max(xyzzyaaak2,xyzzyaaai2)
if(xyzzyaaag2==1)then
call skip(xyzzyaaaa2,xyzzyaaai2*(xyzzyaaah2+5))
else
call skip(xyzzyaaaa2,xyzzyaaai2*(1+xyzzyaaah2+5))
endif
case(2)
xyzzyaaap2=0
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaah2
xyzzyaaap2=max(xyzzyaaap2,xyzzyaaah2)
call skip(xyzzyaaaa2,1+xyzzyaaah2)
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP','Bad orbital index. Bug &
&<1>.')
end select
read(xyzzyaaaa2,'(a)',iostat=xyzzyaaab2)char_80
if(am_master.and.trim(adjustl(char_80))/='END SET '//trim(i2s(xyzzyaaa&
&d2)))call errstop('EXPOT_WFDET_SETUP',"Expected to find 'END SET "//t&
&rim(i2s(xyzzyaaad2))//"' in expot.data .")
enddo
do xyzzyaaaf2=1,xyzzyaaaj1
xyzzyaaae2=xyzzyaaac1(xyzzyaaaf2)
if(xyzzyaaae2==0)cycle
select case(xyzzyaaaf2)
case(1)
allocate(xyzzyaaal1(xyzzyaaae2),xyzzyaaam1(xyzzyaaae2),xyzzyaaan1(xyzz&
&yaaae2),xyzzyaaaq1(xyzzyaaae2),xyzzyaaao1(xyzzyaaak2,xyzzyaaae2),xyzz&
&yaaau1(2,xyzzyaaak1,xyzzyaaae2),xyzzyaaav1(xyzzyaaak1,xyzzyaaae2),xyz&
&zyaaaw1(xyzzyaaaj2,xyzzyaaae2),xyzzyaaaz1(xyzzyaaak2),xyzzyaaba1(xyzz&
&yaaak2),xyzzyaabb1(xyzzyaaak2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'EXPOT_WFDET_SETUP','1a')
if(xyzzyaaal2==1.or.xyzzyaaal2==0)then
allocate(xyzzyaaat1(xyzzyaaaj2,xyzzyaaak2,xyzzyaaae2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'EXPOT_WFDET_SETUP','1b')
endif
if(xyzzyaaal2==2.or.xyzzyaaal2==0)then
allocate(xyzzyaaas1(xyzzyaaaj2,xyzzyaaak2,xyzzyaaae2),xyzzyaaar1(xyzzy&
&aaak2,xyzzyaaae2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'EXPOT_WFDET_SETUP','1c')
endif
case(2)
allocate(xyzzyaabd1(xyzzyaaae2),        xyzzyaabe1(xyzzyaaap2,xyzzyaaa&
&e2),xyzzyaabf1(xyzzyaaap2,xyzzyaaae2),xyzzyaabg1(xyzzyaaap2,xyzzyaaae&
&2),xyzzyaabh1(xyzzyaaap2,xyzzyaaae2),xyzzyaabi1(xyzzyaaap2,xyzzyaaae2&
&),xyzzyaabj1(xyzzyaaap2,xyzzyaaae2),stat=xyzzyaaac2)
call check_alloc(xyzzyaaac2,'EXPOT_WFDET_SETUP','2a')
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP','Bad orbital index. Bug &
&<2>.')
end select
enddo
rewind(xyzzyaaaa2)
do
read(xyzzyaaaa2,'(a)',iostat=xyzzyaaab2)char_80
if(xyzzyaaab2/=0)call errstop('EXPOT_WFDET_SETUP',"Label 'START EXPOT_&
&WFN' not found in expot.data on second pass (!!) .")
if(trim(adjustl(char_80))=='START EXPOT_WFN')exit
enddo
call skip(xyzzyaaaa2,4)
do xyzzyaaad2=1,xyzzyaaab1
call skip(xyzzyaaaa2,1)
if(.not.any(heg_orbtype==100+xyzzyaaad2))then
if(am_master)call wout('SET '//trim(i2s(xyzzyaaad2))//' skipped (unuse&
&d).')
do
read(xyzzyaaaa2,'(a)',end=10,err=20)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaad2)))exit
enddo
cycle
endif
xyzzyaaae2=xyzzyaaae1(xyzzyaaad2)
xyzzyaaaf2=xyzzyaaad1(xyzzyaaad2)
xyzzyaaag1(:,xyzzyaaad2)=0.d0
xyzzyaaah1(:,xyzzyaaad2)=0.d0
if(am_master)then
call wout('SET '//trim(i2s(xyzzyaaad2)))
select case(xyzzyaaaf2)
case(1)
call wout(' Orbital type                : FOURIER (periodic)')
case(2)
call wout(' Orbital type                : GAUSSIAN (aperiodic)')
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP','Orbital code '//trim(i2&
&s(xyzzyaaaf2))//' not recognized. Bug.')
end select
endif
call skip(xyzzyaaaa2,5)
read(xyzzyaaaa2,'(a)',end=10,err=20)char_80
char_80=adjustl(char_80)
xyzzyaaam2=scan(char_80,' ')
if(xyzzyaaam2/=0)char_80=char_80(:xyzzyaaam2)
select case(trim(char_80))
case('ISOTROPIC','isotropic')
xyzzyaaag2=1
case('X','x')
xyzzyaaag2=2
case('Y','y')
xyzzyaaag2=3
case('Z','z')
xyzzyaaag2=4
case('CUSTOM','custom')
xyzzyaaag2=5
case('A1','a1')
xyzzyaaag2=6
case('A2','a2')
xyzzyaaag2=7
case('A3','a3')
xyzzyaaag2=8
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP','Unrecognized direction &
&reading EXPOT_WFN in expot.data.')
end select
if(am_master.and.xyzzyaaag2>5.and.periodicity<xyzzyaaag2-5)call errsto&
&p('EXPOT_WFDET_SETUP','Direction has been defined using lattice vecto&
&r inconsistent with periodicity.')
xyzzyaaaf1(xyzzyaaad2)=xyzzyaaag2
xyzzyaaaq2=0.d0
select case(xyzzyaaag2)
case(1)
dumchar='ISOTROPIC'
case(2)
xyzzyaaaq2(1)=1.d0
dumchar='X axis'
case(3)
xyzzyaaaq2(2)=1.d0
dumchar='Y axis'
case(4)
xyzzyaaaq2(3)=1.d0
dumchar='Z axis'
case(5)
backspace(xyzzyaaaa2)
read(xyzzyaaaa2,*,end=10,err=20)dumchar,xyzzyaaaq2
dumchar='Custom: ('//trim(r2s(xyzzyaaaq2(1),'(e10.4)'))//', '//trim(r2&
&s(xyzzyaaaq2(2),'(e10.4)'))//', '//trim(r2s(xyzzyaaaq2(3),'(e10.4)'))&
&//')'
case(6)
xyzzyaaaq2=a1
dumchar='1st lattice vector'
case(7)
xyzzyaaaq2=a2
dumchar='2nd lattice vector'
case(8)
xyzzyaaaq2=a3
dumchar='3rd lattice vector'
end select
if(am_master)call wout(' Direction of change         : '//trim(dumchar&
&))
xyzzyaaag1(:,xyzzyaaad2)=xyzzyaaaq2
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaah1(1:3,xyzzyaaad2)
if(am_master)then
dumchar='('//trim(r2s(xyzzyaaah1(1,xyzzyaaad2),'(e10.4)'))//', '//trim&
&(r2s(xyzzyaaah1(2,xyzzyaaad2),'(e10.4)'))//', '//trim(r2s(xyzzyaaah1(&
&3,xyzzyaaad2),'(e10.4)'))//')'
call wout(' Origin of orbitals          : '//trim(dumchar))
endif
call skip(xyzzyaaaa2,1)
select case(xyzzyaaaf2)
case(1)
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaaq1(xyzzyaaae2)
if(am_master)then
dumchar=r2s(xyzzyaaaq1(xyzzyaaae2),'(f20.12)')
call wout(' Wavelength of series        : '//trim(dumchar))
endif
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)char_80
select case(trim(adjustl(char_80)))
case('ODD','odd')
xyzzyaaag2=1
dumchar='ODD (sine series)'
case('EVEN','even')
xyzzyaaag2=2
dumchar='EVEN (cosine series)'
case('NONE','none')
xyzzyaaag2=0
dumchar='NONE (full sine+cosine series)'
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP','Unrecognized symmetry o&
&f Fourier series on second pass (!!).')
end select
if(am_master)call wout(' Symmetry of series          : '//trim(dumchar&
&))
xyzzyaaan1(xyzzyaaae2)=xyzzyaaag2
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaam1(xyzzyaaae2)
if(am_master)call wout(' Number of terms in series   : '//trim(i2s(xyz&
&zyaaam1(xyzzyaaae2))))
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaah2
xyzzyaaal1(xyzzyaaae2)=xyzzyaaah2
if(am_master)call wout(' Number of bands in set      : '//trim(i2s(xyz&
&zyaaal1(xyzzyaaae2))))
do xyzzyaaan2=1,xyzzyaaah2
read(xyzzyaaaa2,'(a)',end=10,err=20)char_80
if(am_master.and.trim(adjustl(char_80))/='START BAND '//trim(i2s(xyzzy&
&aaan2)))call errstop('EXPOT_WFDET_SETUP',"Expected to find 'START BAN&
&D "//trim(i2s(xyzzyaaan2))//"' in expot.data .")
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaao1(xyzzyaaan2,xyzzyaaae2)
call skip(xyzzyaaaa2,1)
if(am_master)call wout(' Band '//trim(i2s(xyzzyaaan2))//' occupation  &
&         : '//trim(i2s(xyzzyaaao1(xyzzyaaan2,xyzzyaaae2))))
if(xyzzyaaag2==0.or.xyzzyaaag2==2)then
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaar1(xyzzyaaan2,xyzzyaaae2)
endif
if(xyzzyaaag2==0)then
do xyzzyaaao2=1,xyzzyaaam1(xyzzyaaae2)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaas1(xyzzyaaao2,xyzzyaaan2,xyzzy&
&aaae2),xyzzyaaat1(xyzzyaaao2,xyzzyaaan2,xyzzyaaae2)
enddo
elseif(xyzzyaaag2==1)then
do xyzzyaaao2=1,xyzzyaaam1(xyzzyaaae2)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaat1(xyzzyaaao2,xyzzyaaan2,xyzzy&
&aaae2)
enddo
elseif(xyzzyaaag2==2)then
do xyzzyaaao2=1,xyzzyaaam1(xyzzyaaae2)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaaas1(xyzzyaaao2,xyzzyaaan2,xyzzy&
&aaae2)
enddo
endif
read(xyzzyaaaa2,'(a)',end=10,err=20)char_80
if(am_master.and.trim(adjustl(char_80))/='END BAND '//trim(i2s(xyzzyaa&
&an2)))call errstop('EXPOT_WFDET_SETUP',"Expected to find 'END BAND "/&
&/trim(i2s(xyzzyaaan2))//"' in expot.data .")
enddo
case(2)
call skip(xyzzyaaaa2,1)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaabd1(xyzzyaaae2)
call skip(xyzzyaaaa2,1)
do xyzzyaaao2=1,xyzzyaabd1(xyzzyaaae2)
read(xyzzyaaaa2,*,end=10,err=20)xyzzyaabe1(xyzzyaaao2,xyzzyaaae2),xyzz&
&yaabf1(xyzzyaaao2,xyzzyaaae2),xyzzyaabg1(xyzzyaaao2,xyzzyaaae2),xyzzy&
&aabh1(xyzzyaaao2,xyzzyaaae2),xyzzyaabi1(xyzzyaaao2,xyzzyaaae2),xyzzya&
&abj1(xyzzyaaao2,xyzzyaaae2)
enddo
case default
if(am_master)call errstop('EXPOT_WFDET_SETUP',"Orbital type '"//trim(a&
&djustl(char_80))//"' in expot.data not recognized.")
end select
call skip(xyzzyaaaa2,1)
enddo
close(xyzzyaaaa2)
if(am_master)call wout()
call xyzzyaabm1
return
10 if(am_master)call errstop('EXPOT_WFDET_SETUP','Reached end-of-file &
&while reading expot.data .')
20 if(am_master)call errstop('EXPOT_WFDET_SETUP','Problem found while &
&reading expot.data .')
end subroutine init_expot_wfdet
subroutine xyzzyaabm1
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3(2),xyzzyaa&
&al3
integer,allocatable :: xyzzyaaam3(:,:),xyzzyaaan3(:)
real(dp) xyzzyaaao3,xyzzyaaap3,xyzzyaaaq3,xyzzyaaar3,xyzzyaaas3
logical xyzzyaaat3
do xyzzyaaaa3=1,xyzzyaaab1
if(xyzzyaaaf1(xyzzyaaaa3)==1)cycle
xyzzyaaao3=sqrt(sum(xyzzyaaag1(:,xyzzyaaaa3)**2))
if(am_master.and.xyzzyaaao3==0.d0)call errstop('SETUP_EXPOT_WFDET','Le&
&ngth of direction vector = 0.')
enddo
if(any(xyzzyaaad1==1))then
allocate(xyzzyaaap1(xyzzyaaab1),stat=xyzzyaaab3)
call check_alloc(xyzzyaaab3,'SETUP_EXPOT_WFDET','1')
endif
do xyzzyaaaa3=1,xyzzyaaab1
xyzzyaaac3=xyzzyaaae1(xyzzyaaaa3)
xyzzyaaal3=xyzzyaaad1(xyzzyaaaa3)
if(xyzzyaaal3==0)cycle
select case(xyzzyaaal3)
case(1)
if(am_master)then
if(periodicity<2)call errstop('SETUP_EXPOT_WFDET','Cannot use Fourier &
&expansions in systems without (at least) 2D periodicity.')
xyzzyaaat3=.false.
if(xyzzyaaaf1(xyzzyaaaa3)/=4)xyzzyaaat3=.true.
if(a1(3)/=0.d0.or.a2(3)/=0.d0)xyzzyaaat3=.true.
if(periodicity==3)then
if(a3(1)/=0.d0.or.a3(2)/=0.d0)xyzzyaaat3=.true.
endif
if(xyzzyaaat3)call errstop('SETUP_EXPOT_WFDET','Fourier expansions mus&
&t be along a3, which should be parallel to z and perpendicular to a1 &
&and a2.')
endif
xyzzyaaas3=1.d0/real(nemax*xyzzyaaam1(xyzzyaaac3),dp)
do xyzzyaaad3=2,nemax
xyzzyaaas3=xyzzyaaas3*real(xyzzyaaad3,dp)**(-1.d0/real(nemax,dp))
enddo
xyzzyaabc1(xyzzyaaaa3)=xyzzyaabc1(xyzzyaaaa3)*xyzzyaaas3
select case(xyzzyaaan1(xyzzyaaac3))
case(0)
xyzzyaaar1(:,xyzzyaaac3)=xyzzyaaar1(:,xyzzyaaac3)*xyzzyaabc1(xyzzyaaaa&
&3)
xyzzyaaas1(:,:,xyzzyaaac3)=xyzzyaaas1(:,:,xyzzyaaac3)*xyzzyaabc1(xyzzy&
&aaaa3)
xyzzyaaat1(:,:,xyzzyaaac3)=xyzzyaaat1(:,:,xyzzyaaac3)*xyzzyaabc1(xyzzy&
&aaaa3)
case(1)
xyzzyaaat1(:,:,xyzzyaaac3)=xyzzyaaat1(:,:,xyzzyaaac3)*xyzzyaabc1(xyzzy&
&aaaa3)
case(2)
xyzzyaaar1(:,xyzzyaaac3)=xyzzyaaar1(:,xyzzyaaac3)*xyzzyaabc1(xyzzyaaaa&
&3)
xyzzyaaas1(:,:,xyzzyaaac3)=xyzzyaaas1(:,:,xyzzyaaac3)*xyzzyaabc1(xyzzy&
&aaaa3)
end select
allocate(xyzzyaaam3(2,xyzzyaaak1),xyzzyaaan3(xyzzyaaak1))
xyzzyaaam3=0
xyzzyaaan3=0
xyzzyaaar3=abs(a1(1)*a2(2)-a1(2)*a2(1))
xyzzyaaaq3=twopi/sqrt(xyzzyaaar3)
xyzzyaaai3=int(sqrt(real(1+2*(xyzzyaaak1-1),dp)))-1
xyzzyaaae3=0
i_loop: do xyzzyaaad3=0,xyzzyaaai3
xyzzyaaah3=0
if(xyzzyaaad3>0)xyzzyaaah3=-xyzzyaaai3
do xyzzyaaaf3=xyzzyaaah3,xyzzyaaai3
if(xyzzyaaad3==0.and.xyzzyaaaf3==0)cycle
xyzzyaaae3=xyzzyaaae3+1
if(xyzzyaaae3>xyzzyaaak1)exit i_loop
xyzzyaaam3(1:2,xyzzyaaae3)=(/xyzzyaaad3,xyzzyaaaf3/)
xyzzyaaan3(xyzzyaaae3)=xyzzyaaad3**2+xyzzyaaaf3**2
enddo
enddo i_loop
do xyzzyaaad3=1,xyzzyaaak1-1
do xyzzyaaaf3=xyzzyaaad3+1,xyzzyaaak1
if(xyzzyaaan3(xyzzyaaaf3)<xyzzyaaan3(xyzzyaaad3))then
xyzzyaaaj3=xyzzyaaan3(xyzzyaaad3)
xyzzyaaan3(xyzzyaaad3)=xyzzyaaan3(xyzzyaaaf3)
xyzzyaaan3(xyzzyaaaf3)=xyzzyaaaj3
xyzzyaaak3=xyzzyaaam3(:,xyzzyaaad3)
xyzzyaaam3(:,xyzzyaaad3)=xyzzyaaam3(:,xyzzyaaaf3)
xyzzyaaam3(:,xyzzyaaaf3)=xyzzyaaak3
endif
enddo
enddo
xyzzyaaau1(1:2,:,xyzzyaaac3)=real(xyzzyaaam3(1:2,:),dp)*xyzzyaaaq3
xyzzyaaav1(:,xyzzyaaac3)=real(xyzzyaaan3(:),dp)*xyzzyaaaq3**2
deallocate(xyzzyaaam3,xyzzyaaan3)
xyzzyaaap1(xyzzyaaac3)=0
do xyzzyaaad3=1,nspin
if(any(heg_orbtype(xyzzyaaad3,:)==100+xyzzyaaaa3))then
xyzzyaaae3=heg_nele(xyzzyaaad3)
if(xyzzyaaae3==0)cycle
xyzzyaaaj3=0
do xyzzyaaaf3=1,xyzzyaaal1(xyzzyaaac3)
xyzzyaaag3=xyzzyaaao1(xyzzyaaaf3,xyzzyaaac3)
if(am_master.and.xyzzyaaag3>0.and.mod(xyzzyaaag3,2)/=1)call errstop('S&
&ETUP_EXPOT_WFDET','Band occupation must be an odd number.')
if(xyzzyaaae3>xyzzyaaag3)then
xyzzyaaae3=xyzzyaaae3-xyzzyaaag3
xyzzyaaaj3=max(xyzzyaaaj3,xyzzyaaag3)
else
if(am_master.and.xyzzyaaae3>0.and.mod(xyzzyaaae3,2)/=1)call errstop('S&
&ETUP_EXPOT_WFDET','Number of particles in last occupied band must be &
&an odd number.')
xyzzyaaaj3=max(xyzzyaaaj3,xyzzyaaae3)
xyzzyaaae3=0
endif
if(xyzzyaaae3==0)exit
enddo
if(am_master.and.xyzzyaaae3>0)call errstop('SETUP_EXPOT_WFDET','Number&
& of particles exceeds sum of occupations.')
if(xyzzyaaaj3>0)xyzzyaaaj3=(xyzzyaaaj3-1)/2
xyzzyaaap1(xyzzyaaac3)=max(xyzzyaaap1(xyzzyaaac3),xyzzyaaaj3)
endif
enddo
xyzzyaaap3=twopi/xyzzyaaaq1(xyzzyaaac3)
xyzzyaaao3=0.d0
do xyzzyaaad3=1,xyzzyaaam1(xyzzyaaac3)
xyzzyaaao3=xyzzyaaao3+xyzzyaaap3
xyzzyaaaw1(xyzzyaaad3,xyzzyaaac3)=xyzzyaaao3
enddo
case(2)
continue
case default
if(am_master)call errstop('SETUP_EXPOT_WFDET','Unknown orbital code. B&
&ug.')
end select
enddo
do xyzzyaaaa3=1,xyzzyaaab1
xyzzyaaal3=xyzzyaaad1(xyzzyaaaa3)
if(xyzzyaaal3==0)cycle
select case(xyzzyaaal3)
case(1)
xyzzyaaaj3=maxval(xyzzyaaap1)
allocate(xyzzyaaax1(xyzzyaaaj3),xyzzyaaay1(xyzzyaaaj3),stat=xyzzyaaab3&
&)
case(2)
continue
case default
if(am_master)call errstop('SETUP_EXPOT_WFDET','Unknown orbital code. B&
&ug.')
end select
enddo
end subroutine xyzzyaabm1
subroutine expot_orb_eval(rvec,jspin,lnorb,norb,orbmask,val,fsd,orbval&
&,orbgrad,orblap,orbsderivs)
implicit none
integer,intent(in) :: jspin,lnorb,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(lnorb,*),orbgrad(3,lnorb,*),orblap(ln&
&orb,*)
real(dp),intent(inout),optional :: orbsderivs(6,lnorb,*)
logical,intent(in) :: val,fsd,orbmask(*)
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4
real(dp) xyzzyaaad4(3)
do xyzzyaaab4=1,xyzzyaaab1
xyzzyaaaa4=xyzzyaaad1(xyzzyaaab4)
if(xyzzyaaaa4==0)cycle
xyzzyaaac4=xyzzyaaae1(xyzzyaaab4)
xyzzyaaad4=rvec-xyzzyaaah1(:,xyzzyaaab4)
select case(xyzzyaaaa4)
case(1)
call xyzzyaabr1(xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,jspin,lnorb,norb,orbm&
&ask,val,fsd,orbval,orbgrad,orblap,orbsderivs)
case(2)
call xyzzyaabn1(xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,jspin,lnorb,norb,orbm&
&ask,val,fsd,orbval,orbgrad,orblap,orbsderivs)
case default
call errstop('EXPOT_ORB_EVAL','Unknown orbital type. Bug.')
end select
enddo
end subroutine expot_orb_eval
subroutine xyzzyaabn1(iset,nset,rvec,jspin,lnorb,norb,orbmask,val,fsd,&
&orbval,orbgrad,orblap,orbsderivs)
implicit none
integer,intent(in) :: iset,nset,jspin,lnorb,norb
logical,intent(in) :: val,fsd,orbmask(*)
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(lnorb,*),orbgrad(3,lnorb,*),orblap(ln&
&orb,*)
real(dp),intent(inout),optional :: orbsderivs(6,lnorb,*)
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5
xyzzyaaac5=0
do xyzzyaaab5=1,xyzzyaaac1(2)
do xyzzyaaaa5=1,xyzzyaabd1(xyzzyaaab5)
xyzzyaaac5=xyzzyaaac5+1
if(.not.orbmask(xyzzyaaac5))cycle
if(val)then
orbval(xyzzyaaac5,1)=xyzzyaabo1(rvec(1),xyzzyaabh1(xyzzyaaaa5,nset),xy&
&zzyaabe1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(2),xyzzyaabi1(xyzzyaaaa5,n&
&set),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(3),xyzzyaabj1(xyzzy&
&aaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
endif
if(fsd)then
orbgrad(1,xyzzyaaac5,1)=xyzzyaabp1(rvec(1),xyzzyaabh1(xyzzyaaaa5,nset)&
&,xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(2),xyzzyaabi1(xyzzyaaaa&
&5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(3),xyzzyaabj1(xy&
&zzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
orbgrad(2,xyzzyaaac5,1)=xyzzyaabo1(rvec(1),xyzzyaabh1(xyzzyaaaa5,nset)&
&,xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabp1(rvec(2),xyzzyaabi1(xyzzyaaaa&
&5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(3),xyzzyaabj1(xy&
&zzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
orbgrad(3,xyzzyaaac5,1)=xyzzyaabo1(rvec(1),xyzzyaabh1(xyzzyaaaa5,nset)&
&,xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(2),xyzzyaabi1(xyzzyaaaa&
&5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabp1(rvec(3),xyzzyaabj1(xy&
&zzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
if(.not.present(orbsderivs))then
orblap(xyzzyaaac5,1)=xyzzyaabq1(rvec(1),xyzzyaabh1(xyzzyaaaa5,nset),xy&
&zzyaabe1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(2),xyzzyaabi1(xyzzyaaaa5,n&
&set),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(3),xyzzyaabj1(xyzzy&
&aaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))+xyzzyaabo1(rvec(1),xyzzyaabh&
&1(xyzzyaaaa5,nset),xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabq1(rvec(2),xy&
&zzyaabi1(xyzzyaaaa5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabo1(rve&
&c(3),xyzzyaabj1(xyzzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))+xyzzyaa&
&bo1(rvec(1),xyzzyaabh1(xyzzyaaaa5,nset),xyzzyaabe1(xyzzyaaaa5,nset))*&
&xyzzyaabo1(rvec(2),xyzzyaabi1(xyzzyaaaa5,nset),xyzzyaabf1(xyzzyaaaa5,&
&nset))*xyzzyaabq1(rvec(3),xyzzyaabj1(xyzzyaaaa5,nset),xyzzyaabg1(xyzz&
&yaaaa5,nset))
else
orbsderivs(1,xyzzyaaac5,1)=xyzzyaabq1(rvec(1),xyzzyaabh1(xyzzyaaaa5,ns&
&et),xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(2),xyzzyaabi1(xyzzya&
&aaa5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(3),xyzzyaabj1&
&(xyzzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
orbsderivs(2,xyzzyaaac5,1)=xyzzyaabo1(rvec(1),xyzzyaabh1(xyzzyaaaa5,ns&
&et),xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabq1(rvec(2),xyzzyaabi1(xyzzya&
&aaa5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(3),xyzzyaabj1&
&(xyzzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
orbsderivs(3,xyzzyaaac5,1)=xyzzyaabo1(rvec(1),xyzzyaabh1(xyzzyaaaa5,ns&
&et),xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(2),xyzzyaabi1(xyzzya&
&aaa5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabq1(rvec(3),xyzzyaabj1&
&(xyzzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
orbsderivs(4,xyzzyaaac5,1)=xyzzyaabp1(rvec(1),xyzzyaabh1(xyzzyaaaa5,ns&
&et),xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabp1(rvec(2),xyzzyaabi1(xyzzya&
&aaa5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(3),xyzzyaabj1&
&(xyzzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
orbsderivs(5,xyzzyaaac5,1)=xyzzyaabp1(rvec(1),xyzzyaabh1(xyzzyaaaa5,ns&
&et),xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabo1(rvec(2),xyzzyaabi1(xyzzya&
&aaa5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabp1(rvec(3),xyzzyaabj1&
&(xyzzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
orbsderivs(6,xyzzyaaac5,1)=xyzzyaabo1(rvec(1),xyzzyaabh1(xyzzyaaaa5,ns&
&et),xyzzyaabe1(xyzzyaaaa5,nset))*xyzzyaabp1(rvec(2),xyzzyaabi1(xyzzya&
&aaa5,nset),xyzzyaabf1(xyzzyaaaa5,nset))*xyzzyaabp1(rvec(3),xyzzyaabj1&
&(xyzzyaaaa5,nset),xyzzyaabg1(xyzzyaaaa5,nset))
endif
endif
enddo
enddo
end subroutine xyzzyaabn1
real(dp) function xyzzyaabo1(r,a,n)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r,a
real(dp) xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6
xyzzyaaaa6=r/a
xyzzyaaab6=xyzzyaaaa6*xyzzyaaaa6
xyzzyaabo1=exp(-.5d0*xyzzyaaab6)
select case(n)
case(0)
continue
case(1)
xyzzyaabo1=xyzzyaabo1*(xyzzyaaaa6+xyzzyaaaa6)
case(2)
xyzzyaabo1=xyzzyaabo1*(-2.d0+xyzzyaaab6+xyzzyaaab6+xyzzyaaab6+xyzzyaaa&
&b6)
case(3)
xyzzyaabo1=xyzzyaabo1*(xyzzyaaaa6+xyzzyaaaa6+xyzzyaaaa6+xyzzyaaaa6)*(-&
&3.d0+xyzzyaaab6+xyzzyaaab6)
case(4)
xyzzyaaac6=xyzzyaaab6*xyzzyaaab6
xyzzyaabo1=xyzzyaabo1*(12.d0-48.d0*xyzzyaaab6+16.d0*xyzzyaaac6)
case(5)
xyzzyaaac6=xyzzyaaab6*xyzzyaaab6
xyzzyaabo1=xyzzyaabo1*(xyzzyaaaa6+xyzzyaaaa6+xyzzyaaaa6+xyzzyaaaa6+xyz&
&zyaaaa6+xyzzyaaaa6+xyzzyaaaa6+xyzzyaaaa6)*(15.d0-20.d0*xyzzyaaab6+xyz&
&zyaaac6+xyzzyaaac6+xyzzyaaac6+xyzzyaaac6)
case default
call errstop('expot_Gaussian','Unknown n. Bug.')
end select
end function xyzzyaabo1
real(dp) function xyzzyaabp1(r,a,n)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r,a
real(dp) xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7
xyzzyaaae7=1.d0/a
xyzzyaaaa7=r*xyzzyaaae7
xyzzyaaab7=xyzzyaaaa7*xyzzyaaaa7
xyzzyaabp1=exp(-.5d0*xyzzyaaab7)*xyzzyaaae7
select case(n)
case(0)
xyzzyaabp1=xyzzyaabp1*(-xyzzyaaaa7)
case(1)
xyzzyaabp1=xyzzyaabp1*(2.d0-xyzzyaaab7-xyzzyaaab7)
case(2)
xyzzyaabp1=xyzzyaabp1*(xyzzyaaaa7+xyzzyaaaa7)*(5.d0-xyzzyaaab7-xyzzyaa&
&ab7)
case(3)
xyzzyaaac7=xyzzyaaab7*xyzzyaaab7
xyzzyaabp1=xyzzyaabp1*(-12.d0+36.d0*xyzzyaaab7-8.d0*xyzzyaaac7)
case(4)
xyzzyaaac7=xyzzyaaab7*xyzzyaaab7
xyzzyaabp1=xyzzyaabp1*(xyzzyaaaa7+xyzzyaaaa7+xyzzyaaaa7+xyzzyaaaa7)*(-&
&27.d0+28.d0*xyzzyaaab7-xyzzyaaac7-xyzzyaaac7-xyzzyaaac7-xyzzyaaac7)
case(5)
xyzzyaaac7=xyzzyaaab7*xyzzyaaab7
xyzzyaaad7=xyzzyaaac7*xyzzyaaab7
xyzzyaabp1=xyzzyaabp1*(120.d0-600.d0*xyzzyaaab7+320.d0*xyzzyaaac7-32.d&
&0*xyzzyaaad7)
case default
call errstop('expot_GaussianGrad','Unknown n. Bug.')
end select
end function xyzzyaabp1
real(dp) function xyzzyaabq1(r,a,n)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r,a
real(dp) xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8
xyzzyaaae8=1.d0/a
xyzzyaaaa8=r*xyzzyaaae8
xyzzyaaab8=xyzzyaaaa8*xyzzyaaaa8
xyzzyaabq1=exp(-.5d0*xyzzyaaab8)*xyzzyaaae8*xyzzyaaae8
select case(n)
case(0)
xyzzyaabq1=xyzzyaabq1*(-1.d0+xyzzyaaab8)
case(1)
xyzzyaabq1=xyzzyaabq1*(xyzzyaaaa8+xyzzyaaaa8)*(-3.d0+xyzzyaaab8)
case(2)
xyzzyaaac8=xyzzyaaab8*xyzzyaaab8
xyzzyaabq1=xyzzyaabq1*(10.d0-22.d0*xyzzyaaab8+4.d0*xyzzyaaac8)
case(3)
xyzzyaaac8=xyzzyaaab8*xyzzyaaab8
xyzzyaabq1=xyzzyaabq1*(xyzzyaaaa8+xyzzyaaaa8+xyzzyaaaa8+xyzzyaaaa8)*(2&
&1.d0-17.d0*xyzzyaaab8+xyzzyaaac8+xyzzyaaac8)
case(4)
xyzzyaaac8=xyzzyaaab8*xyzzyaaab8
xyzzyaaad8=xyzzyaaab8*xyzzyaaac8
xyzzyaabq1=xyzzyaabq1*(-108.d0+444.d0*xyzzyaaab8-192.d0*xyzzyaaac8+16.&
&d0*xyzzyaaad8)
case(5)
xyzzyaaac8=xyzzyaaab8*xyzzyaaab8
xyzzyaaad8=xyzzyaaab8*xyzzyaaac8
xyzzyaabq1=xyzzyaabq1*(xyzzyaaaa8+xyzzyaaaa8+xyzzyaaaa8+xyzzyaaaa8+xyz&
&zyaaaa8+xyzzyaaaa8+xyzzyaaaa8+xyzzyaaaa8)*(-165.d0+235.d0*xyzzyaaab8-&
&64.d0*xyzzyaaac8+xyzzyaaad8+xyzzyaaad8+xyzzyaaad8+xyzzyaaad8)
case default
call errstop('expot_GaussianLap','Unknown n. Bug.')
end select
end function xyzzyaabq1
subroutine xyzzyaabr1(iset,nset,rvec,jspin,lnorb,norb,orbmask,val,fsd,&
&orbval,orbgrad,orblap,orbsderivs)
implicit none
integer,intent(in) :: iset,nset,jspin,lnorb,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(lnorb,*),orbgrad(3,lnorb,*),orblap(ln&
&orb,*)
real(dp),intent(inout),optional :: orbsderivs(6,lnorb,*)
logical,intent(in) :: val,fsd,orbmask(*)
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,xyzzyaaae9,xyzzyaa&
&af9,xyzzyaaag9,xyzzyaaah9,xyzzyaaai9,xyzzyaaaj9,xyzzyaaak9
real(dp) xyzzyaaal9,xyzzyaaam9,xyzzyaaan9,xyzzyaaao9,xyzzyaaap9,xyzzya&
&aaq9,xyzzyaaar9,xyzzyaaas9,xyzzyaaat9,xyzzyaaau9,xyzzyaaav9,xyzzyaaaw&
&9,xyzzyaaax9,xyzzyaaay9,xyzzyaaaz9,xyzzyaaba9,xyzzyaabb9,xyzzyaabc9,x&
&yzzyaabd9,xyzzyaabe9,xyzzyaabf9,xyzzyaabg9,xyzzyaabh9,xyzzyaabi9,xyzz&
&yaabj9,xyzzyaabk9,xyzzyaabl9,xyzzyaabm9,xyzzyaabn9,xyzzyaabo9,xyzzyaa&
&bp9
xyzzyaaaa9=100+iset
xyzzyaaad9=xyzzyaaan1(nset)
xyzzyaaae9=xyzzyaaal1(nset)
xyzzyaaag9=xyzzyaaam1(nset)
xyzzyaaap9=rvec(1)
xyzzyaaaq9=rvec(2)
xyzzyaaar9=rvec(3)
do xyzzyaaai9=1,xyzzyaaap1(nset)
xyzzyaaal9=xyzzyaaap9*xyzzyaaau1(1,xyzzyaaai9,nset)+xyzzyaaaq9*xyzzyaa&
&au1(2,xyzzyaaai9,nset)
xyzzyaaax1(xyzzyaaai9)=sin(xyzzyaaal9)
xyzzyaaay1(xyzzyaaai9)=cos(xyzzyaaal9)
enddo
if(val.and..not.fsd)then
xyzzyaaaz1=0.d0
select case(xyzzyaaad9)
case(0)
do xyzzyaaaf9=1,xyzzyaaae9
xyzzyaaam9=xyzzyaaar1(xyzzyaaaf9,nset)
do xyzzyaaah9=1,xyzzyaaag9
xyzzyaaba9=xyzzyaaas1(xyzzyaaah9,xyzzyaaaf9,nset)
xyzzyaabb9=xyzzyaaat1(xyzzyaaah9,xyzzyaaaf9,nset)
xyzzyaaas9=xyzzyaaaw1(xyzzyaaah9,nset)
xyzzyaaau9=xyzzyaaas9*xyzzyaaar9
xyzzyaaam9=xyzzyaaam9+xyzzyaaba9*cos(xyzzyaaau9)+xyzzyaabb9*sin(xyzzya&
&aau9)
enddo
xyzzyaaaz1(xyzzyaaaf9)=xyzzyaaaz1(xyzzyaaaf9)+xyzzyaaam9
enddo
case(1)
do xyzzyaaaf9=1,xyzzyaaae9
xyzzyaaam9=0.d0
do xyzzyaaah9=1,xyzzyaaag9
xyzzyaabb9=xyzzyaaat1(xyzzyaaah9,xyzzyaaaf9,nset)
xyzzyaaas9=xyzzyaaaw1(xyzzyaaah9,nset)
xyzzyaaau9=xyzzyaaas9*xyzzyaaar9
xyzzyaaam9=xyzzyaaam9+xyzzyaabb9*sin(xyzzyaaau9)
enddo
xyzzyaaaz1(xyzzyaaaf9)=xyzzyaaaz1(xyzzyaaaf9)+xyzzyaaam9
enddo
case(2)
do xyzzyaaaf9=1,xyzzyaaae9
xyzzyaaam9=xyzzyaaar1(xyzzyaaaf9,nset)
do xyzzyaaah9=1,xyzzyaaag9
xyzzyaaba9=xyzzyaaas1(xyzzyaaah9,xyzzyaaaf9,nset)
xyzzyaaas9=xyzzyaaaw1(xyzzyaaah9,nset)
xyzzyaaau9=xyzzyaaas9*xyzzyaaar9
xyzzyaaam9=xyzzyaaam9+xyzzyaaba9*cos(xyzzyaaau9)
enddo
xyzzyaaaz1(xyzzyaaaf9)=xyzzyaaaz1(xyzzyaaaf9)+xyzzyaaam9
enddo
end select
else
xyzzyaaaz1=0.d0
xyzzyaaba1=0.d0
xyzzyaabb1=0.d0
select case(xyzzyaaad9)
case(0)
do xyzzyaaaf9=1,xyzzyaaae9
xyzzyaaam9=xyzzyaaar1(xyzzyaaaf9,nset)
xyzzyaaan9=0.d0
xyzzyaaao9=0.d0
do xyzzyaaah9=1,xyzzyaaag9
xyzzyaaba9=xyzzyaaas1(xyzzyaaah9,xyzzyaaaf9,nset)
xyzzyaabb9=xyzzyaaat1(xyzzyaaah9,xyzzyaaaf9,nset)
xyzzyaaas9=xyzzyaaaw1(xyzzyaaah9,nset)
xyzzyaaau9=xyzzyaaas9*xyzzyaaar9
xyzzyaaav9=cos(xyzzyaaau9)
xyzzyaaat9=sin(xyzzyaaau9)
xyzzyaaaw9=xyzzyaaba9*xyzzyaaas9
xyzzyaaay9=xyzzyaaaw9*xyzzyaaas9
xyzzyaaax9=xyzzyaabb9*xyzzyaaas9
xyzzyaaaz9=xyzzyaaax9*xyzzyaaas9
xyzzyaaam9=xyzzyaaam9+xyzzyaaba9*xyzzyaaav9
xyzzyaaam9=xyzzyaaam9+xyzzyaabb9*xyzzyaaat9
xyzzyaaan9=xyzzyaaan9-xyzzyaaaw9*xyzzyaaat9
xyzzyaaao9=xyzzyaaao9-xyzzyaaay9*xyzzyaaav9
xyzzyaaan9=xyzzyaaan9+xyzzyaaax9*xyzzyaaav9
xyzzyaaao9=xyzzyaaao9-xyzzyaaaz9*xyzzyaaat9
enddo
xyzzyaaaz1(xyzzyaaaf9)=xyzzyaaaz1(xyzzyaaaf9)+xyzzyaaam9
xyzzyaaba1(xyzzyaaaf9)=xyzzyaaba1(xyzzyaaaf9)+xyzzyaaan9
xyzzyaabb1(xyzzyaaaf9)=xyzzyaabb1(xyzzyaaaf9)+xyzzyaaao9
enddo
case(1)
do xyzzyaaaf9=1,xyzzyaaae9
xyzzyaaam9=0.d0
xyzzyaaan9=0.d0
xyzzyaaao9=0.d0
do xyzzyaaah9=1,xyzzyaaag9
xyzzyaabb9=xyzzyaaat1(xyzzyaaah9,xyzzyaaaf9,nset)
xyzzyaaas9=xyzzyaaaw1(xyzzyaaah9,nset)
xyzzyaaau9=xyzzyaaas9*xyzzyaaar9
xyzzyaaax9=xyzzyaabb9*xyzzyaaas9
xyzzyaaaz9=xyzzyaaax9*xyzzyaaas9
xyzzyaaat9=sin(xyzzyaaau9)
xyzzyaaam9=xyzzyaaam9+xyzzyaabb9*xyzzyaaat9
xyzzyaaan9=xyzzyaaan9+xyzzyaaax9*cos(xyzzyaaau9)
xyzzyaaao9=xyzzyaaao9-xyzzyaaaz9*xyzzyaaat9
enddo
xyzzyaaaz1(xyzzyaaaf9)=xyzzyaaaz1(xyzzyaaaf9)+xyzzyaaam9
xyzzyaaba1(xyzzyaaaf9)=xyzzyaaba1(xyzzyaaaf9)+xyzzyaaan9
xyzzyaabb1(xyzzyaaaf9)=xyzzyaabb1(xyzzyaaaf9)+xyzzyaaao9
enddo
case(2)
do xyzzyaaaf9=1,xyzzyaaae9
xyzzyaaam9=xyzzyaaar1(xyzzyaaaf9,nset)
xyzzyaaan9=0.d0
xyzzyaaao9=0.d0
do xyzzyaaah9=1,xyzzyaaag9
xyzzyaaba9=xyzzyaaas1(xyzzyaaah9,xyzzyaaaf9,nset)
xyzzyaaas9=xyzzyaaaw1(xyzzyaaah9,nset)
xyzzyaaau9=xyzzyaaas9*xyzzyaaar9
xyzzyaaaw9=xyzzyaaba9*xyzzyaaas9
xyzzyaaay9=xyzzyaaaw9*xyzzyaaas9
xyzzyaaav9=cos(xyzzyaaau9)
xyzzyaaam9=xyzzyaaam9+xyzzyaaba9*xyzzyaaav9
xyzzyaaan9=xyzzyaaan9-xyzzyaaaw9*sin(xyzzyaaau9)
xyzzyaaao9=xyzzyaaao9-xyzzyaaay9*xyzzyaaav9
enddo
xyzzyaaaz1(xyzzyaaaf9)=xyzzyaaaz1(xyzzyaaaf9)+xyzzyaaam9
xyzzyaaba1(xyzzyaaaf9)=xyzzyaaba1(xyzzyaaaf9)+xyzzyaaan9
xyzzyaabb1(xyzzyaaaf9)=xyzzyaabb1(xyzzyaaaf9)+xyzzyaaao9
enddo
end select
endif
if(all(heg_orbtype(jspin,:)/=xyzzyaaaa9))call errstop('EXPOT_FOURIER',&
&'Orbital type does not match for any determinant.')
if(val)orbval(1:norb,1:real1_complex2)=0.d0
if(fsd)then
orbgrad(1:3,1:norb,1:real1_complex2)=0.d0
orblap(1:norb,1:real1_complex2)=0.d0
if(present(orbsderivs))orbsderivs(1:6,1:norb,1:real1_complex2)=0.d0
endif
if(val)then
xyzzyaaab9=1
xyzzyaaak9=1
xyzzyaaaj9=0
xyzzyaaaf9=1
xyzzyaabc9=xyzzyaaaz1(1)
do
if(xyzzyaaaj9==0)then
if(orbmask(xyzzyaaab9))then
orbval(xyzzyaaab9,1)=xyzzyaabc9
xyzzyaaak9=xyzzyaaak9+1
xyzzyaaab9=xyzzyaaab9+1
else
xyzzyaaak9=xyzzyaaak9+1
xyzzyaaab9=xyzzyaaab9+1
endif
else
if(orbmask(xyzzyaaab9))then
orbval(xyzzyaaab9,1)=xyzzyaabc9*xyzzyaaay1(xyzzyaaaj9)
orbval(xyzzyaaab9+1,1)=xyzzyaabc9*xyzzyaaax1(xyzzyaaaj9)
xyzzyaaak9=xyzzyaaak9+2
xyzzyaaab9=xyzzyaaab9+2
else
xyzzyaaak9=xyzzyaaak9+1
xyzzyaaab9=xyzzyaaab9+1
endif
endif
if(xyzzyaaab9>norb)exit
if(xyzzyaaak9>xyzzyaaao1(xyzzyaaaf9,nset))then
xyzzyaaak9=1
xyzzyaaaj9=0
xyzzyaaaf9=xyzzyaaaf9+1
xyzzyaabc9=xyzzyaaaz1(xyzzyaaaf9)
else
xyzzyaaaj9=xyzzyaaaj9+1
endif
enddo
endif
if(fsd)then
xyzzyaaab9=1
xyzzyaaak9=1
xyzzyaaaj9=0
xyzzyaaaf9=1
xyzzyaabc9=xyzzyaaaz1(1)
xyzzyaabd9=xyzzyaaba1(1)
xyzzyaabe9=xyzzyaabb1(1)
do
if(xyzzyaaaj9==0)then
if(orbmask(xyzzyaaab9))then
orbgrad(1,xyzzyaaab9,1)=0.d0
orbgrad(2,xyzzyaaab9,1)=0.d0
orbgrad(3,xyzzyaaab9,1)=xyzzyaabd9
orblap(xyzzyaaab9,1)=xyzzyaabe9
if(present(orbsderivs))then
orbsderivs(1,xyzzyaaab9,1)=0.d0
orbsderivs(2,xyzzyaaab9,1)=0.d0
orbsderivs(3,xyzzyaaab9,1)=xyzzyaabe9
orbsderivs(4,xyzzyaaab9,1)=0.d0
orbsderivs(5,xyzzyaaab9,1)=0.d0
orbsderivs(6,xyzzyaaab9,1)=0.d0
endif
xyzzyaaak9=xyzzyaaak9+1
xyzzyaaab9=xyzzyaaab9+1
else
xyzzyaaak9=xyzzyaaak9+1
xyzzyaaab9=xyzzyaaab9+1
endif
else
if(orbmask(xyzzyaaab9))then
xyzzyaabf9=xyzzyaaau1(1,xyzzyaaaj9,nset)
xyzzyaabg9=xyzzyaaau1(2,xyzzyaaaj9,nset)
xyzzyaabh9=xyzzyaaav1(xyzzyaaaj9,nset)
xyzzyaabi9=xyzzyaaax1(xyzzyaaaj9)
xyzzyaabj9=xyzzyaaay1(xyzzyaaaj9)
xyzzyaabk9=xyzzyaabc9*xyzzyaabi9
xyzzyaabl9=xyzzyaabc9*xyzzyaabj9
xyzzyaabm9=xyzzyaabd9*xyzzyaabi9
xyzzyaabn9=xyzzyaabd9*xyzzyaabj9
xyzzyaabo9=xyzzyaabe9*xyzzyaabi9
xyzzyaabp9=xyzzyaabe9*xyzzyaabj9
orbgrad(1,xyzzyaaab9,1)=-xyzzyaabk9*xyzzyaabf9
orbgrad(2,xyzzyaaab9,1)=-xyzzyaabk9*xyzzyaabg9
orbgrad(3,xyzzyaaab9,1)=xyzzyaabn9
orblap(xyzzyaaab9,1)=xyzzyaabp9-xyzzyaabl9*xyzzyaabh9
if(present(orbsderivs))then
orbsderivs(1,xyzzyaaab9,1)=-xyzzyaabl9*xyzzyaabf9*xyzzyaabf9
orbsderivs(2,xyzzyaaab9,1)=-xyzzyaabl9*xyzzyaabg9*xyzzyaabg9
orbsderivs(3,xyzzyaaab9,1)=xyzzyaabp9
orbsderivs(4,xyzzyaaab9,1)=-xyzzyaabl9*xyzzyaabf9*xyzzyaabg9
orbsderivs(5,xyzzyaaab9,1)=-xyzzyaabm9*xyzzyaabf9
orbsderivs(6,xyzzyaaab9,1)=-xyzzyaabm9*xyzzyaabg9
endif
xyzzyaaac9=xyzzyaaab9+1
orbgrad(1,xyzzyaaac9,1)=xyzzyaabl9*xyzzyaabf9
orbgrad(2,xyzzyaaac9,1)=xyzzyaabl9*xyzzyaabg9
orbgrad(3,xyzzyaaac9,1)=xyzzyaabm9
orblap(xyzzyaaac9,1)=xyzzyaabo9-xyzzyaabk9*xyzzyaabh9
if(present(orbsderivs))then
orbsderivs(1,xyzzyaaac9,1)=-xyzzyaabk9*xyzzyaabf9*xyzzyaabf9
orbsderivs(2,xyzzyaaac9,1)=-xyzzyaabk9*xyzzyaabg9*xyzzyaabg9
orbsderivs(3,xyzzyaaac9,1)=xyzzyaabo9
orbsderivs(4,xyzzyaaac9,1)=-xyzzyaabk9*xyzzyaabf9*xyzzyaabg9
orbsderivs(5,xyzzyaaac9,1)=xyzzyaabn9*xyzzyaabf9
orbsderivs(6,xyzzyaaac9,1)=xyzzyaabn9*xyzzyaabg9
endif
xyzzyaaak9=xyzzyaaak9+2
xyzzyaaab9=xyzzyaaab9+2
else
xyzzyaaak9=xyzzyaaak9+2
xyzzyaaab9=xyzzyaaab9+2
endif
endif
if(xyzzyaaab9>norb)exit
if(xyzzyaaak9>xyzzyaaao1(xyzzyaaaf9,nset))then
xyzzyaaak9=1
xyzzyaaaj9=0
xyzzyaaaf9=xyzzyaaaf9+1
xyzzyaabc9=xyzzyaaaz1(xyzzyaaaf9)
xyzzyaabd9=xyzzyaaba1(xyzzyaaaf9)
xyzzyaabe9=xyzzyaabb1(xyzzyaaaf9)
else
xyzzyaaaj9=xyzzyaaaj9+1
endif
enddo
endif
end subroutine xyzzyaabr1
subroutine get_exwfdet_orbmap(row_offset,norb,orbmap)
implicit none
integer,intent(inout) :: row_offset(nspin),norb,orbmap(nemax,nspin,nde&
&t)
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10,xyzzyaaad10,xyzzyaaae10
do xyzzyaaaa10=1,xyzzyaaab1
xyzzyaaab10=xyzzyaaad1(xyzzyaaaa10)
if(xyzzyaaab10==0)cycle
select case(xyzzyaaab10)
case(1)
do xyzzyaaad10=1,nspin
do xyzzyaaac10=1,heg_nele(xyzzyaaad10)
orbmap(row_offset(xyzzyaaad10)+xyzzyaaac10,xyzzyaaad10,1)=norb+xyzzyaa&
&ac10
enddo
do xyzzyaaae10=2,ndet
orbmap(row_offset(xyzzyaaad10)+1:,xyzzyaaad10,xyzzyaaae10)=orbmap(row_&
&offset(xyzzyaaad10)+1:,xyzzyaaad10,1)
enddo
row_offset(xyzzyaaad10)=row_offset(xyzzyaaad10)+heg_nele(xyzzyaaad10)
enddo
norb=norb+maxval(heg_nele)
case(2)
do xyzzyaaad10=1,nspin
do xyzzyaaac10=1,heg_nele(xyzzyaaad10)
orbmap(row_offset(xyzzyaaad10)+xyzzyaaac10,xyzzyaaad10,1)=norb+xyzzyaa&
&ac10
enddo
do xyzzyaaae10=2,ndet
orbmap(row_offset(xyzzyaaad10)+1:,xyzzyaaad10,xyzzyaaae10)=orbmap(row_&
&offset(xyzzyaaad10)+1:,xyzzyaaad10,1)
enddo
row_offset(xyzzyaaad10)=row_offset(xyzzyaaad10)+heg_nele(xyzzyaaad10)
enddo
norb=norb+maxval(heg_nele)
case default
call errstop('GET_EXWFDET_ORBMAP','Unknown orbital type. Bug.')
end select
enddo
end subroutine get_exwfdet_orbmap
subroutine get_exwfdet_ndesc(ndesc_int,ndesc_dp)
implicit none
integer,intent(inout) :: ndesc_int,ndesc_dp
ndesc_int=xyzzyaabk1
ndesc_dp=xyzzyaabl1
end subroutine get_exwfdet_ndesc
subroutine get_exwfdet_orbdesc(norb,ndesc_int,ndesc_dp,orbdesc_int,orb&
&desc_dp)
implicit none
integer,intent(in) :: norb,ndesc_int,ndesc_dp
integer,intent(inout) :: orbdesc_int(ndesc_int,norb)
real(dp),intent(inout) :: orbdesc_dp(ndesc_dp,norb)
end subroutine get_exwfdet_orbdesc
end module slaarnaar

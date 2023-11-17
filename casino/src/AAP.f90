module slaarnaap
use dsp
use slaarnaag,   only : twopi
use file_utils,  only : open_units,skip
use format_utils,only : wout,i2s,r2s,r2s2
use slaarnabg,    only : a1,a2,a3,periodicity,dimensionality,painv
use slaarnabt,   only : lookup,interp_nev
use parallel,    only : am_master
use run_control, only : errstop,errwarn,check_alloc
use store,       only : nspin
implicit none
private
public read_expot,setup_expot,eval_expot,expot_is_finite
integer,parameter :: xyzzyaaaa1=14,xyzzyaaab1=3,xyzzyaaac1=2
integer xyzzyaaad1,xyzzyaaae1,xyzzyaaaf1(xyzzyaaaa1),xyzzyaaag1(xyzzya&
&aab1),xyzzyaaah1(xyzzyaaac1)
integer,allocatable :: xyzzyaaai1(:),xyzzyaaaj1(:),xyzzyaaak1(:),xyzzy&
&aaal1(:,:),xyzzyaaam1(:),xyzzyaaan1(:),xyzzyaaao1(:,:)
real(dp) xyzzyaaap1
real(dp),allocatable :: xyzzyaaaq1(:,:),xyzzyaaar1(:,:,:),xyzzyaaas1(:&
&,:)
logical,allocatable :: xyzzyaaat1(:)
character(80) expot_title
real(dp),allocatable :: xyzzyaaau1(:),xyzzyaaav1(:)
logical,allocatable :: xyzzyaaaw1(:),xyzzyaaax1(:),xyzzyaaay1(:)
real(dp),allocatable :: xyzzyaaaz1(:)
integer,allocatable :: xyzzyaaba1(:)
real(dp),allocatable :: xyzzyaabb1(:),xyzzyaabc1(:,:),xyzzyaabd1(:)
integer,allocatable :: xyzzyaabe1(:)
real(dp),allocatable :: xyzzyaabf1(:,:),xyzzyaabg1(:,:)
integer,allocatable :: xyzzyaabh1(:)
real(dp),allocatable :: xyzzyaabi1(:,:),xyzzyaabj1(:,:)
real(dp),allocatable :: xyzzyaabk1(:)
real(dp),allocatable :: xyzzyaabl1(:),xyzzyaabm1(:),xyzzyaabn1(:)
logical,allocatable :: xyzzyaabo1(:),xyzzyaabp1(:),xyzzyaabq1(:)
real(dp),allocatable :: xyzzyaabr1(:),xyzzyaabs1(:)
real(dp),allocatable :: xyzzyaabt1(:),xyzzyaabu1(:),xyzzyaabv1(:)
integer,allocatable :: xyzzyaabw1(:),xyzzyaabx1(:)
real(dp),allocatable :: xyzzyaaby1(:)
real(dp),allocatable :: xyzzyaabz1(:),xyzzyaaca1(:,:),xyzzyaacb1(:)
integer xyzzyaacc1
integer,allocatable :: xyzzyaacd1(:),xyzzyaace1(:)
real(dp),allocatable :: xyzzyaacf1(:),xyzzyaacg1(:)
real(dp),allocatable :: xyzzyaach1(:,:),xyzzyaaci1(:,:)
logical xyzzyaacj1,xyzzyaack1,xyzzyaacl1
integer,allocatable :: xyzzyaacm1(:),xyzzyaacn1(:)
real(dp),allocatable :: xyzzyaaco1(:)
real(dp),allocatable :: xyzzyaacp1(:,:),xyzzyaacq1(:,:)
integer,allocatable :: xyzzyaacr1(:)
real(dp),allocatable :: xyzzyaacs1(:)
real(dp),allocatable :: xyzzyaact1(:,:),xyzzyaacu1(:,:)
integer :: xyzzyaacv1(3)
real(dp),allocatable :: xyzzyaacw1(:,:,:)
real(dp) xyzzyaacx1,xyzzyaacy1,xyzzyaacz1(3)
contains
subroutine read_expot
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaal2&
&,xyzzyaaam2,xyzzyaaan2,xyzzyaaao2,xyzzyaaap2,xyzzyaaaq2,xyzzyaaar2
logical xyzzyaaas2
character(8) dumchar
character(80) char_80,tmpr,tmpr2,tmpr3
call open_units(xyzzyaaaf2,xyzzyaaae2)
if(xyzzyaaae2/=0)call errstop('READ_EXPOT','Unable to find free i/o un&
&it.')
open(unit=xyzzyaaaf2,file='expot.data',status='old',action='read',iost&
&at=xyzzyaaae2)
if(xyzzyaaae2/=0)call errstop('READ_EXPOT','Problem opening expot.data&
&.')
if(am_master)then
call wout('External potential')
call wout('==================')
call wout('Reading external potential from expot.data file.')
call wout()
endif
xyzzyaaad1=2
do
read(xyzzyaaaf2,'(a)',iostat=xyzzyaaae2)char_80
if(xyzzyaaae2/=0)exit
if(trim(adjustl(char_80))=='START VERSION')then
read(xyzzyaaaf2,*,end=10,err=10)xyzzyaaad1
read(xyzzyaaaf2,'(a)',end=10,err=10)char_80
if(am_master.and.trim(adjustl(char_80))/='END VERSION')call errstop('R&
&EAD_EXPOT',"Expected to find 'END VERSION' in expot.data .")
exit
endif
enddo
rewind(xyzzyaaaf2)
do
read(xyzzyaaaf2,'(a)',iostat=xyzzyaaae2)char_80
if(xyzzyaaae2>0.and.am_master)call errstop('READ_EXPOT','Problem readi&
&ng expot.data.')
if(xyzzyaaae2<0.and.am_master)call errstop('READ_EXPOT','Could not fin&
&d "START EXPOT" in expot.data.')
if(trim(adjustl(char_80))=='START EXPOT')exit
enddo
read(xyzzyaaaf2,*,err=10,end=10)
read(xyzzyaaaf2,'(a)',err=10,end=10)expot_title
if(am_master)call wout('Title: '//trim(adjustl(expot_title)))
read(xyzzyaaaf2,*,err=10,end=10)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaae1
if(am_master)then
call wout('Number of sets : '//trim(i2s(xyzzyaaae1)))
call wout()
endif
allocate(xyzzyaaat1(xyzzyaaae1),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','expot_periodic')
allocate(xyzzyaaaj1(xyzzyaaae1),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','expot_reptype')
allocate(xyzzyaaak1(xyzzyaaae1),xyzzyaaal1(nspin,xyzzyaaae1),stat=xyzz&
&yaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','expot_npart')
allocate(xyzzyaaai1(xyzzyaaae1),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','expot_direction')
allocate(xyzzyaaaq1(3,xyzzyaaae1),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','expot_dvec')
allocate(xyzzyaaam1(xyzzyaaae1),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','expot_ncopies')
allocate(xyzzyaaan1(xyzzyaaae1*10),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','setmap')
allocate(xyzzyaaao1(xyzzyaaae1,xyzzyaaaa1),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','invsetmap')
xyzzyaaag2=0
xyzzyaaaf1(:)=0
xyzzyaaag1(:)=0
xyzzyaaah1(:)=0
xyzzyaaai2=0
xyzzyaaaj2=0
xyzzyaaak2=0
xyzzyaaal2=0
xyzzyaaam2=0
xyzzyaacj1=.false.
xyzzyaack1=.false.
xyzzyaacl1=.false.
xyzzyaacc1=0
do xyzzyaaag2=1,xyzzyaaae1
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(am_master.and.trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzya&
&aag2)))call errstop('READ_EXPOT',"Expected to find 'START SET "//    &
&       trim(i2s(xyzzyaaag2))//"' in expot.data .")
if(xyzzyaaad1>1)then
call skip(xyzzyaaaf2,1)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
select case(trim(adjustl(char_80)))
case('all','ALL','All')
continue
case default
call skip(xyzzyaaaf2,2)
end select
endif
call skip(xyzzyaaaf2,1)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(trim(adjustl(char_80))=='PERIODIC')then
xyzzyaaas2=.true.
elseif(trim(adjustl(char_80))=='APERIODIC')then
xyzzyaaas2=.false.
else
xyzzyaaas2=.false.
call errstop('READ_EXPOT','Unrecognized periodicity in expot.data.')
endif
call skip(xyzzyaaaf2,1)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(.not.xyzzyaaas2)then
if(trim(adjustl(char_80))=='SQUARE')then
xyzzyaaaj1(xyzzyaaag2)=1
xyzzyaaaf1(1)=xyzzyaaaf1(1)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(1)
xyzzyaaao1(xyzzyaaaf1(1),1)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
call skip(xyzzyaaaf2,xyzzyaaac2+6)
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
elseif(trim(adjustl(char_80))=='LINEAR')then
xyzzyaaaj1(xyzzyaaag2)=2
xyzzyaaaf1(2)=xyzzyaaaf1(2)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(2)
xyzzyaaao1(xyzzyaaaf1(2),2)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
call skip(xyzzyaaaf2,xyzzyaaac2+4)
if(xyzzyaaac2/=1)call errstop('READ_EXPOT','Only one origin allowed fo&
&r linear potential.')
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
elseif(trim(adjustl(char_80))=='ANALYTIC')then
xyzzyaaaj1(xyzzyaaag2)=3
xyzzyaaaf1(3)=xyzzyaaaf1(3)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(3)
xyzzyaaao1(xyzzyaaaf1(3),3)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+3)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
char_80=adjustl(char_80)
xyzzyaaao2=scan(char_80,' ')
char_80=char_80(:xyzzyaaao2)
if(trim(char_80)=='GAUSSIAN')then
xyzzyaaag1(1)=xyzzyaaag1(1)+1
elseif(trim(char_80)=='HARMONIC')then
xyzzyaaag1(2)=xyzzyaaag1(2)+1
elseif(trim(char_80)=='SLAB')then
xyzzyaaag1(3)=xyzzyaaag1(3)+1
else
call errstop('READ_EXPOT','Unknown ANALYTIC potential type.')
endif
elseif(trim(adjustl(char_80))=='GAUSSIAN')then
xyzzyaaaj1(xyzzyaaag2)=4
xyzzyaaaf1(4)=xyzzyaaaf1(4)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(4)
xyzzyaaao1(xyzzyaaaf1(4),4)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+3)
read(xyzzyaaaf2,'(a)',err=10,end=10)xyzzyaaac2
xyzzyaaai2=max(xyzzyaaai2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+1)
elseif(trim(adjustl(char_80))=='NUMERICAL')then
xyzzyaaaj1(xyzzyaaag2)=5
xyzzyaaaf1(5)=xyzzyaaaf1(5)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(5)
xyzzyaaao1(xyzzyaaaf1(5),5)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaaj2=max(xyzzyaaaj2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+1)
else
call errstop('READ_EXPOT','Aperiodic representation type not recognize&
&d.)')
endif
else
if(trim(adjustl(char_80))=='CONSTANT')then
xyzzyaaaj1(xyzzyaaag2)=6
xyzzyaaaf1(6)=xyzzyaaaf1(6)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(6)
xyzzyaaao1(xyzzyaaaf1(6),6)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+4)
elseif(trim(adjustl(char_80))=='SQUARE_PERIODIC')then
xyzzyaaaj1(xyzzyaaag2)=7
xyzzyaaaf1(7)=xyzzyaaaf1(7)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(7)
xyzzyaaao1(xyzzyaaaf1(7),7)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+8)
elseif(trim(adjustl(char_80))=='SAWTOOTH')then
xyzzyaaaj1(xyzzyaaag2)=8
xyzzyaaaf1(8)=xyzzyaaaf1(8)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(8)
xyzzyaaao1(xyzzyaaaf1(8),8)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+6)
elseif(trim(adjustl(char_80))=='COSINE')then
xyzzyaaaj1(xyzzyaaag2)=9
xyzzyaaaf1(9)=xyzzyaaaf1(9)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(9)
xyzzyaaao1(xyzzyaaaf1(9),9)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+6)
elseif(trim(adjustl(char_80))=='ANALYTIC_PERIODIC')then
xyzzyaaaj1(xyzzyaaag2)=10
xyzzyaaaf1(10)=xyzzyaaaf1(10)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(10)
xyzzyaaao1(xyzzyaaaf1(10),10)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+5)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
char_80=adjustl(char_80)
xyzzyaaao2=scan(char_80,' ')
char_80=char_80(:xyzzyaaao2)
if(trim(char_80)=='GAUSSIAN')then
xyzzyaaah1(1)=xyzzyaaah1(1)+1
elseif(trim(char_80)=='SLAB')then
xyzzyaaah1(2)=xyzzyaaah1(2)+1
else
call errstop('READ_EXPOT','Unknown ANALYTIC_PERIODIC potential type.')
endif
elseif(trim(adjustl(char_80))=='FOURIER')then
xyzzyaaaj1(xyzzyaaag2)=11
xyzzyaaaf1(11)=xyzzyaaaf1(11)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(11)
xyzzyaaao1(xyzzyaaaf1(11),11)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+5)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(trim(adjustl(char_80))=='ODD')then
xyzzyaacj1=.true.
elseif(trim(adjustl(char_80))=='EVEN')then
xyzzyaack1=.true.
elseif(trim(adjustl(char_80))=='NONE')then
xyzzyaacl1=.true.
else
call errstop('READ_EXPOT','Error reading symmetry for external potenti&
&al of Fourier type.')
endif
read(xyzzyaaaf2,*,err=10,end=10)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaaa2
xyzzyaacc1=max(xyzzyaacc1,xyzzyaaaa2)
call skip(xyzzyaaaf2,xyzzyaaaa2+1)
elseif(trim(adjustl(char_80))=='GAUSSIAN_PERIODIC')then
xyzzyaaaj1(xyzzyaaag2)=12
xyzzyaaaf1(12)=xyzzyaaaf1(12)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(12)
xyzzyaaao1(xyzzyaaaf1(12),12)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+5)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaak2=max(xyzzyaaak2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+1)
elseif(trim(adjustl(char_80))=='NUMERICAL_PERIODIC')then
xyzzyaaaj1(xyzzyaaag2)=13
xyzzyaaaf1(13)=xyzzyaaaf1(13)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(13)
xyzzyaaao1(xyzzyaaaf1(13),13)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=max(xyzzyaaam2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+5)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaal2=max(xyzzyaaal2,xyzzyaaac2)
call skip(xyzzyaaaf2,xyzzyaaac2+1)
elseif(trim(adjustl(char_80))=='BLIP_PERIODIC')then
xyzzyaaaj1(xyzzyaaag2)=14
xyzzyaaaf1(14)=xyzzyaaaf1(14)+1
xyzzyaaan1(xyzzyaaag2)=xyzzyaaaf1(14)
xyzzyaaao1(xyzzyaaaf1(14),14)=xyzzyaaag2
call skip(xyzzyaaaf2,3)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaac2
xyzzyaaam2=xyzzyaaac2
if(xyzzyaaac2/=1)call errstop('READ_EXPOT','Only one potential allowed&
& for blips.')
call skip(xyzzyaaaf2,4)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaacv1
call skip(xyzzyaaaf2,1+xyzzyaacv1(1)*xyzzyaacv1(2)*xyzzyaacv1(3))
else
call errstop('READ_EXPOT','Periodic representation type not recognized&
&.')
endif
endif
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(am_master.and.trim(adjustl(char_80))/='END SET '//trim(i2s(xyzzyaaa&
&g2)))call errstop('READ_EXPOT',"Expected to find 'END SET "//        &
&   trim(i2s(xyzzyaaag2))//"' in expot.data .")
enddo
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(am_master.and.trim(adjustl(char_80))/='END EXPOT')call errstop('REA&
&D_EXPOT',"Expected to find 'END EXPOT' in expot.data .")
do xyzzyaaaa2=1,xyzzyaaaa1
xyzzyaaac2=xyzzyaaaf1(xyzzyaaaa2)
if(xyzzyaaac2==0)cycle
select case(xyzzyaaaa2)
case(1)
allocate(xyzzyaaau1(xyzzyaaac2),xyzzyaaav1(xyzzyaaac2),xyzzyaaaw1(xyzz&
&yaaac2),xyzzyaaax1(xyzzyaaac2),xyzzyaaay1(xyzzyaaac2),stat=xyzzyaaad2&
&)
if(xyzzyaaad2==0)then
xyzzyaaaw1=.false.
xyzzyaaax1=.false.
xyzzyaaay1=.false.
endif
case(2)
allocate(xyzzyaaaz1(xyzzyaaac2),stat=xyzzyaaad2)
case(3)
xyzzyaaan2=0
if(xyzzyaaag1(1)>0)xyzzyaaan2=max(xyzzyaaan2,1)
if(xyzzyaaag1(2)>0)xyzzyaaan2=max(xyzzyaaan2,0)
if(xyzzyaaag1(3)>0)xyzzyaaan2=max(xyzzyaaan2,3)
if(xyzzyaaan2>0)then
allocate(xyzzyaaba1(xyzzyaaac2),xyzzyaabb1(xyzzyaaac2),xyzzyaabd1(xyzz&
&yaaac2),xyzzyaabc1(xyzzyaaac2,xyzzyaaan2),stat=xyzzyaaad2)
else
allocate(xyzzyaaba1(xyzzyaaac2),xyzzyaabb1(xyzzyaaac2),xyzzyaabd1(xyzz&
&yaaac2),stat=xyzzyaaad2)
endif
case(4)
allocate(xyzzyaabe1(xyzzyaaac2),xyzzyaabf1(xyzzyaaai2,xyzzyaaac2),xyzz&
&yaabg1(xyzzyaaai2,xyzzyaaac2),stat=xyzzyaaad2)
case(5)
allocate(xyzzyaabh1(xyzzyaaac2),xyzzyaabi1(xyzzyaaaj2,xyzzyaaac2),xyzz&
&yaabj1(xyzzyaaaj2,xyzzyaaac2),stat=xyzzyaaad2)
case(6)
allocate(xyzzyaabk1(xyzzyaaac2),stat=xyzzyaaad2)
case(7)
allocate(xyzzyaabl1(xyzzyaaac2),xyzzyaabm1(xyzzyaaac2),xyzzyaabn1(xyzz&
&yaaac2),xyzzyaabo1(xyzzyaaac2),xyzzyaabp1(xyzzyaaac2),xyzzyaabq1(xyzz&
&yaaac2),stat=xyzzyaaad2)
if(xyzzyaaad2==0)then
xyzzyaabo1=.false.
xyzzyaabp1=.false.
xyzzyaabq1=.false.
endif
case(8)
allocate(xyzzyaabr1(xyzzyaaac2),xyzzyaabs1(xyzzyaaac2),stat=xyzzyaaad2&
&)
case(9)
allocate(xyzzyaabt1(xyzzyaaac2),xyzzyaabu1(xyzzyaaac2),xyzzyaabv1(xyzz&
&yaaac2),stat=xyzzyaaad2)
case(10)
xyzzyaaan2=0
if(xyzzyaaah1(1)>0)xyzzyaaan2=max(xyzzyaaan2,1)
if(xyzzyaaah1(2)>0)xyzzyaaan2=max(xyzzyaaan2,4)
if(xyzzyaaan2>0)then
allocate(xyzzyaabw1(xyzzyaaac2),xyzzyaabz1(xyzzyaaac2),xyzzyaacb1(xyzz&
&yaaac2),xyzzyaaca1(xyzzyaaac2,xyzzyaaan2),xyzzyaaby1(xyzzyaaac2),xyzz&
&yaabx1(xyzzyaaac2),stat=xyzzyaaad2)
else
allocate(xyzzyaabw1(xyzzyaaac2),xyzzyaabz1(xyzzyaaac2),xyzzyaacb1(xyzz&
&yaaac2),xyzzyaaby1(xyzzyaaac2),xyzzyaabx1(xyzzyaaac2),stat=xyzzyaaad2&
&)
endif
case(11)
allocate(xyzzyaacf1(xyzzyaaac2),xyzzyaace1(xyzzyaaac2),xyzzyaacd1(xyzz&
&yaaac2),xyzzyaacg1(xyzzyaaac2),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','Fourier <1>')
if(xyzzyaack1.or.xyzzyaacl1)allocate(xyzzyaach1(0:xyzzyaacc1,xyzzyaaac&
&2),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','Fourier a_n')
if(xyzzyaacj1.or.xyzzyaacl1)allocate(xyzzyaaci1(xyzzyaacc1,xyzzyaaac2)&
&,stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','Fourier b_n')
case(12)
allocate(xyzzyaacm1(xyzzyaaac2),xyzzyaaco1(xyzzyaaac2),xyzzyaacp1(xyzz&
&yaaak2,xyzzyaaac2),xyzzyaacq1(xyzzyaaak2,xyzzyaaac2),xyzzyaacn1(xyzzy&
&aaac2),stat=xyzzyaaad2)
case(13)
allocate(xyzzyaacr1(xyzzyaaac2),xyzzyaacs1(xyzzyaaac2),xyzzyaact1(xyzz&
&yaaal2,xyzzyaaac2),xyzzyaacu1(xyzzyaaal2,xyzzyaaac2),stat=xyzzyaaad2)
case(14)
allocate(xyzzyaacw1(0:xyzzyaacv1(1)-1,0:xyzzyaacv1(2)-1,0:xyzzyaacv1(3&
&)-1), stat=xyzzyaaad2)
case default
call errstop('READ_EXPOT','Error deciding which sets contain which pot&
&entials. Bug.')
end select
call check_alloc(xyzzyaaad2,'READ_EXPOT','1')
enddo
allocate(xyzzyaaar1(3,xyzzyaaam2,xyzzyaaae1),stat=xyzzyaaad2)
call check_alloc(xyzzyaaad2,'READ_EXPOT','expot_origin')
rewind(xyzzyaaaf2)
xyzzyaaah2=0
xyzzyaaaf1(:)=0
do
read(xyzzyaaaf2,'(a)',iostat=xyzzyaaae2)char_80
if(xyzzyaaae2>0.and.am_master)call errstop('READ_EXPOT','Problem readi&
&ng expot.data.')
if(xyzzyaaae2<0.and.am_master)call errstop('READ_EXPOT','Could not fin&
&d "START EXPOT" in expot.data.')
if(trim(adjustl(char_80))=='START EXPOT')exit
enddo
call skip(xyzzyaaaf2,4)
do xyzzyaaag2=1,xyzzyaaae1
call skip(xyzzyaaaf2,1)
if(am_master)call wout('SET '//trim(i2s(xyzzyaaag2)))
if(xyzzyaaad1>1)then
call skip(xyzzyaaaf2,1)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
select case(trim(adjustl(char_80)))
case('all','ALL','All')
xyzzyaaak1(xyzzyaaag2)=nspin
do xyzzyaaaa2=1,nspin
xyzzyaaal1(xyzzyaaaa2,xyzzyaaag2)=xyzzyaaaa2
enddo
case default
read(char_80,*,err=10,end=10)xyzzyaaak1(xyzzyaaag2)
if(xyzzyaaak1(xyzzyaaag2)<1.or.xyzzyaaak1(xyzzyaaag2)>nspin)call errst&
&op('READ_EXPOT','Problem with number of particle types.')
call skip(xyzzyaaaf2,1)
xyzzyaaal1(:,xyzzyaaag2)=0
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaal1(1:xyzzyaaak1(xyzzyaaag2),xy&
&zzyaaag2)
end select
else
xyzzyaaak1(xyzzyaaag2)=nspin
do xyzzyaaaa2=1,nspin
xyzzyaaal1(xyzzyaaaa2,xyzzyaaag2)=xyzzyaaaa2
enddo
endif
read(xyzzyaaaf2,*,err=10,end=10)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(trim(adjustl(char_80))=='PERIODIC')then
xyzzyaaat1(xyzzyaaag2)=.true.
elseif(trim(adjustl(char_80))=='APERIODIC')then
xyzzyaaat1(xyzzyaaag2)=.false.
else
if(am_master)call errstop('READ_EXPOT','Unrecognized periodicity')
endif
read(xyzzyaaaf2,*,err=10,end=10)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(.not.xyzzyaaat1(xyzzyaaag2))then
if(trim(adjustl(char_80))=='SQUARE')then
if(xyzzyaaaj1(xyzzyaaag2)/=1)xyzzyaaah2=1
if(am_master)call wout('Potential is aperiodic and of type SQUARE')
elseif(trim(adjustl(char_80))=='LINEAR')then
if(xyzzyaaaj1(xyzzyaaag2)/=2)xyzzyaaah2=1
if(am_master)call wout('Potential is aperiodic and of type LINEAR')
elseif(trim(adjustl(char_80))=='ANALYTIC')then
if(xyzzyaaaj1(xyzzyaaag2)/=3)xyzzyaaah2=1
if(am_master)call wout('Potential is aperiodic and of type ANALYTIC')
elseif(trim(adjustl(char_80))=='GAUSSIAN')then
if(xyzzyaaaj1(xyzzyaaag2)/=4)xyzzyaaah2=1
if(am_master)call wout('Potential is aperiodic and of type GAUSSIAN')
elseif(trim(adjustl(char_80))=='NUMERICAL')then
if(xyzzyaaaj1(xyzzyaaag2)/=5)xyzzyaaah2=1
if(am_master)call wout('Potential is aperiodic and of type NUMERICAL')
else
call errstop('READ_EXPOT','Aperiodic representation type not recognize&
&d.)')
endif
else
if(trim(adjustl(char_80))=='CONSTANT')then
if(xyzzyaaaj1(xyzzyaaag2)/=6)xyzzyaaah2=1
if(am_master)call wout('Potential is of type CONSTANT')
elseif(trim(adjustl(char_80))=='SQUARE_PERIODIC')then
if(xyzzyaaaj1(xyzzyaaag2)/=7)xyzzyaaah2=1
if(am_master)call wout('Potential is of type SQUARE_PERIODIC')
elseif(trim(adjustl(char_80))=='SAWTOOTH')then
if(xyzzyaaaj1(xyzzyaaag2)/=8)xyzzyaaah2=1
if(am_master)call wout('Potential is periodic and of type SAWTOOTH')
elseif(trim(adjustl(char_80))=='COSINE')then
if(xyzzyaaaj1(xyzzyaaag2)/=9)xyzzyaaah2=1
if(am_master)call wout('Potential is periodic and of type COSINE')
elseif(trim(adjustl(char_80))=='ANALYTIC_PERIODIC')then
if(xyzzyaaaj1(xyzzyaaag2)/=10)xyzzyaaah2=1
if(am_master)call wout('Potential is type ANALYTIC_PERIODIC')
elseif(trim(adjustl(char_80))=='FOURIER')then
if(xyzzyaaaj1(xyzzyaaag2)/=11)xyzzyaaah2=1
if(am_master)call wout('Potential is periodic and of type FOURIER')
elseif(trim(adjustl(char_80))=='GAUSSIAN_PERIODIC')then
if(xyzzyaaaj1(xyzzyaaag2)/=12)xyzzyaaah2=1
if(am_master)call wout('Potential is of type GAUSSIAN_PERIODIC')
elseif(trim(adjustl(char_80))=='NUMERICAL_PERIODIC')then
if(xyzzyaaaj1(xyzzyaaag2)/=13)xyzzyaaah2=1
if(am_master)call wout('Potential is of type NUMERICAL_PERIODIC')
elseif(trim(adjustl(char_80))=='BLIP_PERIODIC')then
if(xyzzyaaaj1(xyzzyaaag2)/=14)xyzzyaaah2=1
if(am_master)call wout('Potential is of type BLIP_PERIODIC')
else
call errstop('READ_EXPOT','Periodic representation type not recognized&
&.)')
endif
if(xyzzyaaah2/=0)call errstop('READ_EXPOT','Mismatch between first and&
& second read of expot.data')
endif
xyzzyaaaq1(:,xyzzyaaag2)=0.d0
read(xyzzyaaaf2,*,err=10,end=10)
read(xyzzyaaaf2,'(a)',err=10,end=10)char_80
if(xyzzyaaaj1(xyzzyaaag2)/=6)then
if(trim(adjustl(char_80))=='ISOTROPIC'.or.trim(adjustl(char_80))=='iso&
&tropic')then
xyzzyaaai1(xyzzyaaag2)=1
if(am_master)call wout('Direction of change         : ISOTROPIC')
elseif(trim(adjustl(char_80))=='X'.or.trim(adjustl(char_80))=='x')then
xyzzyaaai1(xyzzyaaag2)=2
xyzzyaaaq1(1,xyzzyaaag2)=1.d0
if(am_master)call wout('Direction of change         : X')
elseif(trim(adjustl(char_80))=='Y'.or.trim(adjustl(char_80))=='y')then
xyzzyaaai1(xyzzyaaag2)=3
xyzzyaaaq1(2,xyzzyaaag2)=1.d0
if(am_master)call wout('Direction of change         : Y')
elseif(trim(adjustl(char_80))=='Z'.or.trim(adjustl(char_80))=='z')then
xyzzyaaai1(xyzzyaaag2)=4
xyzzyaaaq1(3,xyzzyaaag2)=1.d0
if(am_master)call wout('Direction of change         : Z')
elseif(trim(adjustl(char_80(1:7)))=='CUSTOM'.or.trim(adjustl(char_80(1&
&:7)))=='custom')then
xyzzyaaai1(xyzzyaaag2)=5
backspace(xyzzyaaaf2)
read(xyzzyaaaf2,*,err=10,end=10)dumchar,(xyzzyaaaq1(xyzzyaaaa2,xyzzyaa&
&ag2),xyzzyaaaa2=1,3)
if(am_master)then
call wout('Direction of change         : CUSTOM')
tmpr=r2s(xyzzyaaaq1(1,xyzzyaaag2),'(f20.5)')
tmpr2=r2s(xyzzyaaaq1(2,xyzzyaaag2),'(f20.5)')
tmpr3=r2s(xyzzyaaaq1(3,xyzzyaaag2),'(f20.5)')
call wout('Custom direction vector     : '//trim(tmpr)//', '//trim(tmp&
&r2)//', '//trim(tmpr3))
endif
elseif(trim(adjustl(char_80))=='A1'.or.trim(adjustl(char_80))=='a1')th&
&en
if(.not.xyzzyaaat1(xyzzyaaag2))then
call errstop('READ_EXPOT','External potentials in *finite* systems can&
&not vary along a lattice vector.')
else
xyzzyaaai1(xyzzyaaag2)=6
xyzzyaaaq1(1:3,xyzzyaaag2)=a1(1:3)
if(am_master)call wout('Direction of change         : A1')
endif
elseif(trim(adjustl(char_80))=='A2'.or.trim(adjustl(char_80))=='a2')th&
&en
if(.not.xyzzyaaat1(xyzzyaaag2))then
call errstop('READ_EXPOT','External potentials in *finite* systems can&
&not vary along a lattice vector.')
else
xyzzyaaai1(xyzzyaaag2)=7
xyzzyaaaq1(1:3,xyzzyaaag2)=a2(1:3)
if(am_master)call wout('Direction of change         : A2')
endif
elseif(trim(adjustl(char_80))=='A3'.or.trim(adjustl(char_80))=='a3')th&
&en
if(.not.xyzzyaaat1(xyzzyaaag2))then
call errstop('READ_EXPOT','External potentials in *finite* systems can&
&not vary along a lattice vector.')
else
xyzzyaaai1(xyzzyaaag2)=8
xyzzyaaaq1(1:3,xyzzyaaag2)=a3(1:3)
if(am_master)call wout('Direction of change         : A3')
endif
else
call errstop('READ_EXPOT','Direction type not recognized.')
endif
else
xyzzyaaai1(xyzzyaaag2)=1
endif
read(xyzzyaaaf2,*,err=10,end=10)
read(xyzzyaaaf2,*,err=10,end=10)xyzzyaaam1(xyzzyaaag2)
if(am_master)call wout('Number of copies of set     : '//trim(i2s(xyzz&
&yaaam1(xyzzyaaag2))))
if(xyzzyaaam1(xyzzyaaag2)<1)then
call wout('Number of potentials to be added for set '//trim(i2s(xyzzya&
&aag2))//' must be greater than zero.')
call errstop('READ_EXPOT','Quitting.')
endif
read(xyzzyaaaf2,*,err=10,end=10)
do xyzzyaaab2=1,xyzzyaaam1(xyzzyaaag2)
read(xyzzyaaaf2,*,err=10,end=10)(xyzzyaaar1(xyzzyaaaa2,xyzzyaaab2,xyzz&
&yaaag2),xyzzyaaaa2=1,3)
if(am_master)then
tmpr=r2s(xyzzyaaar1(1,xyzzyaaab2,xyzzyaaag2),'(f20.5)')
tmpr2=r2s(xyzzyaaar1(2,xyzzyaaab2,xyzzyaaag2),'(f20.5)')
tmpr3=r2s(xyzzyaaar1(3,xyzzyaaab2,xyzzyaaag2),'(f20.5)')
if(xyzzyaaab2==1)then
if(xyzzyaaam1(xyzzyaaag2)>1)then
call wout('Origin for each copy        : '//trim(tmpr)//', '//trim(tmp&
&r2)//', '//trim(tmpr3))
else
call wout('Origin                      : '//trim(tmpr)//', '//trim(tmp&
&r2)//', '//trim(tmpr3))
endif
else
call wout('                              '//trim(tmpr)//', '//trim(tmp&
&r2)//', '//trim(tmpr3))
endif
endif
enddo
read(xyzzyaaaf2,*,err=10,end=10)
select case(xyzzyaaaj1(xyzzyaaag2))
case(1)
xyzzyaaaf1(1)=xyzzyaaaf1(1)+1
xyzzyaaac2=xyzzyaaaf1(1)
read(xyzzyaaaf2,*,err=20,end=20)
read(xyzzyaaaf2,*,err=20,end=20)xyzzyaaav1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=20,end=20)
read(xyzzyaaaf2,*,err=20,end=20)dumchar
select case(trim(dumchar))
case('inf','INF','+inf','+INF')
xyzzyaaaw1(xyzzyaaac2)=.true.
xyzzyaaau1(xyzzyaaac2)=1.d0
case('-inf','-INF')
xyzzyaaax1(xyzzyaaac2)=.true.
xyzzyaaau1(xyzzyaaac2)=-1.d0
case default
backspace(xyzzyaaaf2)
read(xyzzyaaaf2,*,err=20,end=20)xyzzyaaau1(xyzzyaaac2)
end select
xyzzyaaay1(xyzzyaaac2)=xyzzyaaaw1(xyzzyaaac2).or.xyzzyaaax1(xyzzyaaac2&
&)
if(am_master)then
tmpr=r2s(xyzzyaaau1(xyzzyaaac2),'(f20.5)')
if(xyzzyaaay1(xyzzyaaac2))tmpr='Infinity'
if(xyzzyaaau1(xyzzyaaac2)>=0.d0)then
call wout('Barrier height              : '//trim(tmpr))
else
call wout('Well depth                  : '//trim(tmpr))
endif
tmpr=r2s(xyzzyaaav1(xyzzyaaac2),'(f20.5)')
call wout('Width                       : '//trim(tmpr))
endif
case(2)
xyzzyaaaf1(2)=xyzzyaaaf1(2)+1
xyzzyaaac2=xyzzyaaaf1(2)
read(xyzzyaaaf2,*,err=25,end=25)
read(xyzzyaaaf2,*,err=25,end=25)xyzzyaaaz1(xyzzyaaac2)
if(am_master)then
tmpr=r2s(xyzzyaaaz1(xyzzyaaac2),'(f20.5)')
call wout('Electric field              : '//trim(tmpr))
endif
case(3)
xyzzyaaaf1(3)=xyzzyaaaf1(3)+1
xyzzyaaac2=xyzzyaaaf1(3)
read(xyzzyaaaf2,*,err=30,end=30)
read(xyzzyaaaf2,'(a)',err=30,end=30)char_80
char_80=adjustl(char_80)
xyzzyaaao2=scan(char_80,' ')
char_80=char_80(:xyzzyaaao2)
if(trim(char_80)=='GAUSSIAN')then
xyzzyaaba1(xyzzyaaac2)=1
backspace(xyzzyaaaf2)
read(xyzzyaaaf2,*,err=30,end=30)dumchar,xyzzyaabb1(xyzzyaaac2),xyzzyaa&
&bc1(xyzzyaaac2,1),xyzzyaabd1(xyzzyaaac2)
if(am_master)then
call wout('Analytic function           : Gaussian a*exp(-b r^2)+c')
tmpr=r2s(xyzzyaabb1(xyzzyaaac2),'(f20.5)')
tmpr2=r2s(xyzzyaabc1(xyzzyaaac2,1),'(f20.5)')
tmpr3=r2s(xyzzyaabd1(xyzzyaaac2),'(f20.5)')
call wout('                               a = '//trim(tmpr)//', b = '/&
&/trim(tmpr2)//', c = '//trim(tmpr3))
endif
elseif(trim(char_80)=='HARMONIC')then
xyzzyaaba1(xyzzyaaac2)=2
backspace(xyzzyaaaf2)
read(xyzzyaaaf2,*,err=30,end=30)dumchar,xyzzyaabb1(xyzzyaaac2),xyzzyaa&
&bd1(xyzzyaaac2)
if(am_master)then
call wout('Analytic function           : harmonic a*r^2+b')
tmpr=r2s(xyzzyaabb1(xyzzyaaac2),'(f20.5)')
tmpr2=r2s(xyzzyaabd1(xyzzyaaac2),'(f20.5)')
call wout('                               a = '//trim(tmpr)//', b = '/&
&/trim(tmpr2))
endif
elseif(trim(char_80)=='SLAB')then
xyzzyaaba1(xyzzyaaac2)=3
backspace(xyzzyaaaf2)
read(xyzzyaaaf2,*,err=30,end=30)dumchar,xyzzyaabc1(xyzzyaaac2,1),xyzzy&
&aabc1(xyzzyaaac2,2)
xyzzyaabc1(xyzzyaaac2,3)=1.d0/xyzzyaabc1(xyzzyaaac2,1)**2
xyzzyaabb1(xyzzyaaac2)=(3*xyzzyaabc1(xyzzyaaac2,1)**2)/(2*xyzzyaabc1(x&
&yzzyaaac2,2)**3)
xyzzyaabd1(xyzzyaaac2)=-(xyzzyaabc1(xyzzyaaac2,1)**2)/(4*xyzzyaabc1(xy&
&zzyaaac2,2)**3)
if(am_master)then
call wout('Analytic function           : interaction with a slab')
tmpr=r2s(xyzzyaabc1(xyzzyaaac2,1),'(f20.5)')
tmpr2=r2s(xyzzyaabc1(xyzzyaaac2,2),'(f20.5)')
call wout('                               L = '//trim(tmpr)//', r_s = &
&'//trim(tmpr2))
endif
else
call errstop('READ_EXPOT','Undefined ANALYTIC potential type.')
endif
case(4)
xyzzyaaaf1(4)=xyzzyaaaf1(4)+1
xyzzyaaac2=xyzzyaaaf1(4)
read(xyzzyaaaf2,*,err=40,end=40)
read(xyzzyaaaf2,*,err=40,end=40)xyzzyaabe1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=40,end=40)
do xyzzyaaaa2=1,xyzzyaabe1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=40,end=40)xyzzyaabf1(xyzzyaaaa2,xyzzyaaac2),xyzz&
&yaabg1(xyzzyaaaa2,xyzzyaaac2)
enddo
if(am_master)then
call wout('Number of Gaussians         : '//trim(i2s(xyzzyaabe1(xyzzya&
&aac2))))
do xyzzyaaaa2=1,xyzzyaabe1(xyzzyaaac2)
tmpr=r2s(xyzzyaabg1(xyzzyaaaa2,xyzzyaaac2),'(f20.5)')
tmpr2=r2s2(xyzzyaabf1(xyzzyaaaa2,xyzzyaaac2),'(f20.5)')
if(xyzzyaaaa2==1)then
call wout('Exponents, coefficients     : '//trim(tmpr)//', '//trim(tmp&
&r2))
else
call wout('                              '//trim(tmpr)//', '//trim(tmp&
&r2))
endif
enddo
endif
case(5)
xyzzyaaaf1(5)=xyzzyaaaf1(5)+1
xyzzyaaac2=xyzzyaaaf1(5)
read(xyzzyaaaf2,*,err=50,end=50)
read(xyzzyaaaf2,*,err=50,end=50)xyzzyaabh1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=50,end=50)
do xyzzyaaaa2=1,xyzzyaabh1(xyzzyaaag2)
read(xyzzyaaaf2,*,err=50,end=50)xyzzyaabi1(xyzzyaaaa2,xyzzyaaac2),xyzz&
&yaabj1(xyzzyaaaa2,xyzzyaaac2)
enddo
if(am_master)then
call wout('Number of grid points       : '//trim(i2s(xyzzyaabh1(xyzzya&
&aac2))))
endif
case(6)
xyzzyaaaf1(6)=xyzzyaaaf1(6)+1
xyzzyaaac2=xyzzyaaaf1(6)
read(xyzzyaaaf2,*,err=60,end=60)
read(xyzzyaaaf2,*,err=60,end=60)xyzzyaabk1(xyzzyaaac2)
if(am_master)then
tmpr=r2s(xyzzyaabk1(xyzzyaaac2),'(f20.5)')
call wout('Value of constant           : '//trim(tmpr))
endif
case(7)
xyzzyaaaf1(7)=xyzzyaaaf1(7)+1
xyzzyaaac2=xyzzyaaaf1(7)
read(xyzzyaaaf2,*,err=70,end=70)
read(xyzzyaaaf2,*,err=70,end=70)xyzzyaabl1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=70,end=70)
read(xyzzyaaaf2,*,err=70,end=70)xyzzyaabm1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=70,end=70)
read(xyzzyaaaf2,*,err=70,end=70)dumchar
select case(trim(dumchar))
case('inf','INF','+inf','+INF')
xyzzyaabo1(xyzzyaaac2)=.true.
xyzzyaabn1(xyzzyaaac2)=1.d0
case('-inf','-INF')
xyzzyaabp1(xyzzyaaac2)=.true.
xyzzyaabn1(xyzzyaaac2)=-1.d0
case default
backspace(xyzzyaaaf2)
read(xyzzyaaaf2,*,err=70,end=70)xyzzyaabn1(xyzzyaaac2)
end select
xyzzyaabq1(xyzzyaaac2)=xyzzyaabo1(xyzzyaaac2).or.xyzzyaabp1(xyzzyaaac2&
&)
if(am_master)then
tmpr=r2s(xyzzyaabn1(xyzzyaaac2),'(f20.5)')
if(xyzzyaabq1(xyzzyaaac2))tmpr='Infinity'
tmpr2=r2s(xyzzyaabm1(xyzzyaaac2),'(f20.5)')
tmpr3=r2s(xyzzyaabl1(xyzzyaaac2),'(f20.5)')
if(xyzzyaabn1(xyzzyaaac2)>=0.d0)then
call wout('Barrier height              : '//trim(tmpr))
call wout('Barrier width               : '//trim(tmpr2))
else
call wout('Well depth                  : '//trim(tmpr))
call wout('Well width                  : '//trim(tmpr2))
endif
call wout('Repeat distance             : '//trim(tmpr3))
endif
case(8)
xyzzyaaaf1(8)=xyzzyaaaf1(8)+1
xyzzyaaac2=xyzzyaaaf1(8)
read(xyzzyaaaf2,*,err=75,end=75)
read(xyzzyaaaf2,*,err=75,end=75)xyzzyaabr1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=75,end=75)
read(xyzzyaaaf2,*,err=75,end=75)xyzzyaabs1(xyzzyaaac2)
if(am_master)then
tmpr=r2s(xyzzyaabs1(xyzzyaaac2),'(f20.5)')
call wout('Electric field              : '//trim(tmpr))
tmpr=r2s(xyzzyaabr1(xyzzyaaac2),'(f20.5)')
call wout('Repeat distance (au)        : '//trim(tmpr))
endif
case(9)
xyzzyaaaf1(9)=xyzzyaaaf1(9)+1
xyzzyaaac2=xyzzyaaaf1(9)
read(xyzzyaaaf2,*,err=77,end=77)
read(xyzzyaaaf2,*,err=77,end=77)xyzzyaabv1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=77,end=77)
read(xyzzyaaaf2,*,err=77,end=77)xyzzyaabt1(xyzzyaaac2)
xyzzyaabu1(xyzzyaaac2)=twopi/xyzzyaabt1(xyzzyaaac2)
if(am_master)then
tmpr=r2s(xyzzyaabv1(xyzzyaaac2),'(f20.5)')
tmpr2=r2s(xyzzyaabt1(xyzzyaaac2),'(f20.5)')
call wout('Amplitude (au)              : '//trim(tmpr))
call wout('Wavelength (au)             : '//trim(tmpr2))
endif
case(10)
xyzzyaaaf1(10)=xyzzyaaaf1(10)+1
xyzzyaaac2=xyzzyaaaf1(10)
read(xyzzyaaaf2,*,err=80,end=80)
read(xyzzyaaaf2,*,err=80,end=80)xyzzyaaby1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=80,end=80)
read(xyzzyaaaf2,'(a)',err=80,end=80)char_80
char_80=adjustl(char_80)
xyzzyaaao2=scan(char_80,' ')
char_80=char_80(:xyzzyaaao2)
if(trim(char_80)=='GAUSSIAN')then
xyzzyaabw1(xyzzyaaac2)=1
backspace(xyzzyaaaf2)
read(xyzzyaaaf2,*,err=80,end=80)dumchar,xyzzyaabz1(xyzzyaaac2),xyzzyaa&
&ca1(xyzzyaaac2,1),xyzzyaacb1(xyzzyaaac2)
if(am_master)then
call wout('Analytic function           : Gaussian a*exp(-b r^2)+c')
tmpr=r2s(xyzzyaabz1(xyzzyaaac2),'(f20.5)')
tmpr2=r2s(xyzzyaaca1(xyzzyaaac2,1),'(f20.5)')
tmpr3=r2s(xyzzyaacb1(xyzzyaaac2),'(f20.5)')
call wout('                            : a = '//trim(tmpr)//', b = '//&
&trim(tmpr2)//', c = '//trim(tmpr3))
tmpr=r2s(xyzzyaaby1(xyzzyaaac2),'(f20.5)')
call wout('Periodic repeat distance    : '//trim(tmpr))
endif
elseif(trim(char_80)=='SLAB')then
xyzzyaabw1(xyzzyaaac2)=2
backspace(xyzzyaaaf2)
read(xyzzyaaaf2,*,err=80,end=80)dumchar,xyzzyaaca1(xyzzyaaac2,1),xyzzy&
&aaca1(xyzzyaaac2,2)
xyzzyaaca1(xyzzyaaac2,3)=xyzzyaaby1(xyzzyaaac2)/xyzzyaaca1(xyzzyaaac2,&
&1)
xyzzyaaca1(xyzzyaaac2,4)=1.d0/xyzzyaaca1(xyzzyaaac2,1)**2
xyzzyaabz1(xyzzyaaac2)=xyzzyaaca1(xyzzyaaac2,1)**2/(8*xyzzyaaca1(xyzzy&
&aaac2,3)*xyzzyaaca1(xyzzyaaac2,2)**3)
xyzzyaacb1(xyzzyaaac2)=(xyzzyaaby1(xyzzyaaac2)-xyzzyaaca1(xyzzyaaac2,1&
&))**2/(8*xyzzyaaca1(xyzzyaaac2,3)*xyzzyaaca1(xyzzyaaac2,2)**3)
if(am_master)then
call wout('Analytic function           : interaction with a slab')
tmpr=r2s(xyzzyaaca1(xyzzyaaac2,1),'(f20.5)')
tmpr2=r2s(xyzzyaaca1(xyzzyaaac2,2),'(f20.5)')
call wout('                            : L = '//trim(tmpr)//', r_s = '&
&//trim(tmpr2))
tmpr=r2s(xyzzyaaby1(xyzzyaaac2),'(f20.5)')
call wout('Periodic repeat distance    : '//trim(tmpr))
endif
else
call errstop('READ_EXPOT','Undefined ANALYTIC_PERIODIC potential type.&
&')
endif
case(11)
xyzzyaaaf1(11)=xyzzyaaaf1(11)+1
xyzzyaaac2=xyzzyaaaf1(11)
read(xyzzyaaaf2,*,err=90,end=90)
read(xyzzyaaaf2,*,err=90,end=90)xyzzyaacf1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=90,end=90)
read(xyzzyaaaf2,'(a)',err=90,end=90)char_80
if(trim(adjustl(char_80))=='ODD')then
xyzzyaace1(xyzzyaaac2)=1
elseif(trim(adjustl(char_80))=='EVEN')then
xyzzyaace1(xyzzyaaac2)=2
elseif(trim(adjustl(char_80))=='NONE')then
xyzzyaace1(xyzzyaaac2)=3
else
call errstop('READ_EXPOT','Error reading symmetry of external potentia&
&l of Fourier type')
endif
read(xyzzyaaaf2,*,err=90,end=90)
read(xyzzyaaaf2,*,err=90,end=90)xyzzyaacd1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=90,end=90)
if(xyzzyaace1(xyzzyaaac2)==2.or.xyzzyaace1(xyzzyaaac2)==3)then
read(xyzzyaaaf2,*,err=90,end=90)xyzzyaach1(0,xyzzyaaac2)
endif
if(xyzzyaace1(xyzzyaaac2)==1)then
do xyzzyaaaa2=1,xyzzyaacd1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=90,end=90)xyzzyaaci1(xyzzyaaaa2,xyzzyaaac2)
enddo
elseif(xyzzyaace1(xyzzyaaac2)==2)then
do xyzzyaaaa2=1,xyzzyaacd1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=90,end=90)xyzzyaach1(xyzzyaaaa2,xyzzyaaac2)
enddo
elseif(xyzzyaace1(xyzzyaaac2)==3)then
do xyzzyaaaa2=1,xyzzyaacd1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=90,end=90)xyzzyaach1(xyzzyaaaa2,xyzzyaaac2),xyzz&
&yaaci1(xyzzyaaaa2,xyzzyaaac2)
enddo
endif
if(am_master)then
tmpr=r2s(xyzzyaacf1(xyzzyaaac2),'(f20.5)')
call wout('Period                      : '//trim(tmpr))
call wout('Number of terms in sum      : '//trim(i2s(xyzzyaacd1(xyzzya&
&aac2))))
if(xyzzyaace1(xyzzyaaac2)==1)then
call wout('Potential is odd function   : Fourier sine series')
elseif(xyzzyaace1(xyzzyaaac2)==2)then
call wout('Potential is even function  : Fourier cosine series')
elseif(xyzzyaace1(xyzzyaaac2)==3)then
call wout('Potential has no symmetry   : Full Fourier series')
endif
endif
case(12)
xyzzyaaaf1(12)=xyzzyaaaf1(12)+1
xyzzyaaac2=xyzzyaaaf1(12)
read(xyzzyaaaf2,*,err=100,end=100)
read(xyzzyaaaf2,*,err=100,end=100)xyzzyaaco1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=100,end=100)
read(xyzzyaaaf2,*,err=100,end=100)xyzzyaacm1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=100,end=100)
do xyzzyaaaa2=1,xyzzyaacm1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=100,end=100)xyzzyaacp1(xyzzyaaaa2,xyzzyaaac2),xy&
&zzyaacq1(xyzzyaaaa2,xyzzyaaac2)
enddo
if(am_master)then
call wout('Number of Gaussians         : '//trim(i2s(xyzzyaacm1(xyzzya&
&aac2))))
tmpr=r2s(xyzzyaaco1(xyzzyaaac2),'(f20.5)')
call wout('Periodic repeat distance    : '//trim(tmpr))
endif
case(13)
xyzzyaaaf1(13)=xyzzyaaaf1(13)+1
xyzzyaaac2=xyzzyaaaf1(13)
read(xyzzyaaaf2,*,err=110,end=110)
read(xyzzyaaaf2,*,err=110,end=110)xyzzyaacs1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=110,end=110)
read(xyzzyaaaf2,*,err=110,end=110)xyzzyaacr1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=110,end=110)
do xyzzyaaaa2=1,xyzzyaacr1(xyzzyaaac2)
read(xyzzyaaaf2,*,err=110,end=110)xyzzyaact1(xyzzyaaaa2,xyzzyaaac2),xy&
&zzyaacu1(xyzzyaaaa2,xyzzyaaac2)
enddo
if(am_master)then
call wout('Number of grid points       : '//trim(i2s(xyzzyaacr1(xyzzya&
&aac2))))
tmpr=r2s(xyzzyaacs1(xyzzyaaac2),'(f20.5)')
call wout('Periodic repeat distance    : '//trim(tmpr))
endif
case(14)
xyzzyaaaf1(14)=xyzzyaaaf1(14)+1
xyzzyaaac2=xyzzyaaaf1(14)
call skip(xyzzyaaaf2,3)
call wout( 'Started reading blip potential')
do xyzzyaaap2=0,xyzzyaacv1(1)-1
do xyzzyaaaq2=0,xyzzyaacv1(2)-1
do xyzzyaaar2=0,xyzzyaacv1(3)-1
read(xyzzyaaaf2,*,end=20,err=120)xyzzyaacw1(xyzzyaaap2,xyzzyaaaq2,xyzz&
&yaaar2)
enddo
enddo
enddo
call wout( 'Finished reading blip potential')
case default
call errstop('READ_EXPOT','Bug has led to confusion about potential ty&
&pe.')
end select
call skip(xyzzyaaaf2,1)
if(am_master)call wout()
enddo
call skip(xyzzyaaaf2,1)
return
10 if(am_master)call errstop('READ_EXPOT','Problem reading expot.data &
&file.')
20 if(am_master)call errstop('READ_EXPOT','Problem reading SQUARE exte&
&rnal potential from expot.data.')
25 if(am_master)call errstop('READ_EXPOT','Problem reading LINEAR pote&
&ntial from expot.data.')
30 if(am_master)call errstop('READ_EXPOT','Problem reading ANALYTIC ex&
&ternal potential from expot.data.')
40 if(am_master)call errstop('READ_EXPOT','Problem reading GAUSSIAN ex&
&ternal potential from expot.data.')
50 if(am_master)call errstop('READ_EXPOT','Problem reading NUMERICAL p&
&otential from expot.data.')
60 if(am_master)call errstop('READ_EXPOT','Problem reading CONSTANT ex&
&ternal potential from expot.data.')
70 if(am_master)call errstop('READ_EXPOT','Problem reading SQUARE_PERI&
&ODIC external potential from expot.data.')
75 if(am_master)call errstop('READ_EXPOT','Problem reading SAWTOOTH ex&
&ternal potential from expot.data.')
77 if(am_master)call errstop('READ_EXPOT','Problem reading COSINE exte&
&rnal potential from expot.data.')
80 if(am_master)call errstop('READ_EXPOT','Problem reading ANALYTIC_PE&
&RIODIC external potential from expot.data.')
90 if(am_master)call errstop('READ_EXPOT','Problem reading FOURIER ext&
&ernalpotential from expot.data.')
100 if(am_master)call errstop('READ_EXPOT','Problem reading GAUSSIAN_P&
&ERIODIC external potential from expot.data.')
110 if(am_master)call errstop('READ_EXPOT','Problem reading NUMERICAL_&
&PERIODIC external potential from expot.data.')
120 if(am_master)call errstop('READ_EXPOT','Problem reading BLIP_PERIO&
&DIC external potential from expot.data.')
end subroutine read_expot
subroutine setup_expot
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3
real(dp) xyzzyaaag3,xyzzyaaah3
logical xyzzyaaai3
if(any(xyzzyaaai1(:)>4))then
allocate(xyzzyaaas1(3,xyzzyaaae1),stat=xyzzyaaaf3)
call check_alloc(xyzzyaaaf3,'SETUP_EXPOT','unitvector')
do xyzzyaaae3=1,xyzzyaaae1
select case(xyzzyaaai1(xyzzyaaae3))
case(5)
xyzzyaacx1=sqrt(xyzzyaaaq1(1,xyzzyaaae3)**2+xyzzyaaaq1(2,xyzzyaaae3)**&
&2+xyzzyaaaq1(3,xyzzyaaae3)**2)
if(xyzzyaacx1==0.d0)call errstop('SETUP_EXPOT','Length of custom direc&
&tion vector is zero.')
xyzzyaaas1(1:3,xyzzyaaae3)=xyzzyaaaq1(1:3,xyzzyaaae3)/xyzzyaacx1
case(6)
xyzzyaacx1=sqrt(a1(1)**2+a1(2)**2+a1(3)**2)
if(xyzzyaacx1==0.d0)call errstop('SETUP_EXPOT','Length of a1 lattice v&
&ector is zero (for some reason).')
xyzzyaaas1(1:3,xyzzyaaae3)=a1(1:3)/xyzzyaacx1
case(7)
xyzzyaacx1=sqrt(a2(1)**2+a2(2)**2+a2(3)**2)
if(xyzzyaacx1==0.d0)call errstop('SETUP_EXPOT','Length of a2 lattice v&
&ector is zero (for some reason).')
xyzzyaaas1(1:3,xyzzyaaae3)=a2(1:3)/xyzzyaacx1
case(8)
xyzzyaacx1=sqrt(a3(1)**2+a3(2)**2+a3(3)**2)
if(xyzzyaacx1==0.d0)call errstop('SETUP_EXPOT','Length of a3 lattice v&
&ector is zero (for some reason).')
xyzzyaaas1(1:3,xyzzyaaae3)=a3(1:3)/xyzzyaacx1
end select
enddo
endif
do xyzzyaaaa3=1,xyzzyaaaa1
xyzzyaaad3=xyzzyaaaf1(xyzzyaaaa3)
if(xyzzyaaad3==0)cycle
xyzzyaaai3=.false.
do xyzzyaaab3=1,xyzzyaaad3
xyzzyaaae3=xyzzyaaao1(xyzzyaaab3,xyzzyaaaa3)
select case(xyzzyaaaa3)
case(1)
xyzzyaaav1(xyzzyaaab3)=xyzzyaaav1(xyzzyaaab3)*0.5d0
case(2)
if(xyzzyaaai1(xyzzyaaae3)==1)call errstop('SETUP_EXPOT','Isotropic lin&
&ear potentials not allowed.')
case(3)
continue
case(4)
continue
case(5)
if(xyzzyaaai1(xyzzyaaae3)==1.and.any(xyzzyaabi1(:,xyzzyaaab3)<0.d0))ca&
&ll errstop('SETUP_EXPOT','For isotropic numerical potentials all grid&
& points must be positive.')
do xyzzyaaac3=1,xyzzyaabh1(xyzzyaaab3)-1
if(xyzzyaabi1(xyzzyaaac3,xyzzyaaab3)>=xyzzyaabi1(xyzzyaaac3+1,xyzzyaaa&
&b3))then
call errstop('SETUP_EXPOT','Grid values for numerical potential do not&
& increase monotonically.')
endif
enddo
case(6)
continue
case(7)
xyzzyaabm1(xyzzyaaab3)=xyzzyaabm1(xyzzyaaab3)*0.5d0
if(xyzzyaabl1(xyzzyaaab3)<xyzzyaabm1(xyzzyaaab3))then
call wout('For SQUARE_PERIODIC representation of external potential, t&
&he repeat distance')
call wout('must be greater than the width of the square well/barrier.'&
&)
call errstop('SETUP_EXPOT','Quitting.')
endif
case(8)
continue
case(9)
continue
case(10)
if(xyzzyaabw1(xyzzyaaae3)==1)then
xyzzyaaap1=abs(log(.1d0)*8.d0)
xyzzyaaah3=.5d0*xyzzyaaby1(xyzzyaaab3)
xyzzyaabx1(xyzzyaaab3)=0
do
xyzzyaaah3=xyzzyaaah3+xyzzyaaby1(xyzzyaaab3)
if(xyzzyaaca1(xyzzyaaab3,1)*xyzzyaaah3*xyzzyaaah3>xyzzyaaap1)exit
xyzzyaabx1(xyzzyaaab3)=xyzzyaabx1(xyzzyaaab3)+1
enddo
elseif(xyzzyaabw1(xyzzyaaae3)==2)then
xyzzyaabx1(xyzzyaaab3)=0
else
call errstop('SETUP_EXPOT','Unknown analytic_periodic function type.')
endif
case(11)
if(xyzzyaacf1(xyzzyaaab3)<=0.d0)then
call errstop('READ_EXPOT','Period of Fourier series expansion of exter&
&nal potential must be positive definite.')
endif
xyzzyaacg1(xyzzyaaab3)=twopi/xyzzyaacf1(xyzzyaaab3)
case(12)
xyzzyaaap1=abs(log(0.1d0)*8.d0)
xyzzyaaac3=xyzzyaacm1(xyzzyaaab3)
xyzzyaaag3=minval(xyzzyaacq1(1:xyzzyaaac3,xyzzyaaab3))
xyzzyaaah3=0.5d0*xyzzyaaco1(xyzzyaaab3)
xyzzyaacn1(xyzzyaaab3)=1
do
xyzzyaacn1(xyzzyaaab3)=xyzzyaacn1(xyzzyaaab3)+1
xyzzyaaah3=xyzzyaaah3+xyzzyaaco1(xyzzyaaab3)
if(xyzzyaaag3*xyzzyaaah3*xyzzyaaah3>xyzzyaaap1)exit
enddo
call wout('Periodic images for SET '//trim(i2s(xyzzyaaab3))//'  : '//t&
&rim(i2s(xyzzyaacn1(xyzzyaaab3))))
xyzzyaaai3=.true.
case(13)
if(xyzzyaaai1(xyzzyaaae3)==1.and.any(xyzzyaact1(:,xyzzyaaab3)<0.d0))ca&
&ll errstop('SETUP_EXPOT','For isotropic numerical potentials all grid&
& points must be positive.')
do xyzzyaaac3=1,xyzzyaacr1(xyzzyaaab3)-1
if(xyzzyaact1(xyzzyaaac3,xyzzyaaab3)>=xyzzyaact1(xyzzyaaac3+1,xyzzyaaa&
&b3))then
call errstop('SETUP_EXPOT','Grid values for numerical potential do not&
& increase monotonically.')
endif
enddo
case(14)
call wout('No additional setup for blip periodic.')
case default
call errstop('SETUP_EXPOT','No such potential representation defined.'&
&)
end select
enddo
if(xyzzyaaai3)call wout()
enddo
if(allocated(xyzzyaacf1))deallocate(xyzzyaacf1)
call xyzzyaada1
end subroutine setup_expot
subroutine xyzzyaada1
use format_utils,only : wordwrap
implicit none
integer xyzzyaaaa4
logical xyzzyaaab4,xyzzyaaac4
xyzzyaaab4=.false.
xyzzyaaac4=.false.
do xyzzyaaaa4=1,xyzzyaaae1
select case(xyzzyaaaj1(xyzzyaaaa4))
case(1:5)
if(periodicity==3)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic external potentials may not be applied to thr&
&ee-dimensionally periodic systems.')
call wout()
xyzzyaaab4=.true.
elseif(periodicity==2)then
if(dimensionality==2)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic external potentials may not be applied to 2D &
&systems with 2D periodicity.')
call wout()
xyzzyaaab4=.true.
endif
if(xyzzyaaai1(xyzzyaaaa4)==1)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic isotropic external potentials may not be used&
& in two-dimensionally periodic systems.')
call wout()
xyzzyaaab4=.true.
elseif(xyzzyaaai1(xyzzyaaaa4)<4)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic external potentials in 2D slab systems may va&
&ry only in the z direction.')
call wout()
xyzzyaaab4=.true.
elseif(xyzzyaaai1(xyzzyaaaa4)>4)then
if(.not.(xyzzyaaas1(1,xyzzyaaaa4)==0.and.xyzzyaaas1(2,xyzzyaaaa4)==0.a&
&nd.xyzzyaaas1(3,xyzzyaaaa4)/=0))then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic external potentials in 2D slab systems may va&
&ry only in the z direction.')
call wout()
xyzzyaaab4=.true.
endif
endif
elseif(periodicity==1)then
if(dimensionality==1)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic external potentials may not be applied to 1D &
&systems with 1D periodicity.')
call wout()
xyzzyaaab4=.true.
elseif(dimensionality==2)then
if(xyzzyaaai1(xyzzyaaaa4)==1)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic isotropic external potentials may not be used&
& in one-dimensionally periodic systems.')
call wout()
xyzzyaaab4=.true.
elseif(xyzzyaaai1(xyzzyaaaa4)<4)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap("Aperiodic potentials in one-dimensionally periodic 'tap&
&e' systems must be applied in the z direction.")
call wout()
xyzzyaaab4=.true.
elseif(xyzzyaaai1(xyzzyaaaa4)>4)then
if(.not.(xyzzyaaas1(1,xyzzyaaaa4)==0.and.xyzzyaaas1(2,xyzzyaaaa4)==0.a&
&nd.xyzzyaaas1(3,xyzzyaaaa4)/=0))then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap("Aperiodic potentials in one-dimensionally periodic 'tap&
&e' systems must be applied in the z direction.")
call wout()
xyzzyaaab4=.true.
endif
endif
elseif(dimensionality==3)then
if(xyzzyaaai1(xyzzyaaaa4)==1)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic isotropic external potentials may not be used&
& in one-dimensionally periodic systems.')
call wout()
xyzzyaaab4=.true.
elseif(xyzzyaaai1(xyzzyaaaa4)==2)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic external potentials in regular one-dimensiona&
&lly periodic systems may only vary in a direction within the yz plane&
&.')
call wout()
xyzzyaaab4=.true.
endif
if(xyzzyaaai1(xyzzyaaaa4)>4)then
if(xyzzyaaas1(1,xyzzyaaaa4)/=0.d0)then
call wout('SET '//trim(i2s(xyzzyaaaa4)))
call wordwrap('Aperiodic external potentials in regular one-dimensiona&
&lly periodic systems may only vary in a direction within the yz plane&
&.')
call wout()
xyzzyaaab4=.true.
endif
endif
endif
endif
case(6)
continue
case(7:13)
if(periodicity==2)then
if(xyzzyaaai1(xyzzyaaaa4)==4)then
call errwarn('READ_EXPOT','external potential periodic in z direction &
&applied to 2D system periodic in xy plane (set '//trim(i2s(xyzzyaaaa4&
&))//').')
xyzzyaaab4=.true.
elseif(xyzzyaaai1(xyzzyaaaa4)>4)then
if(xyzzyaaas1(3,xyzzyaaaa4)/=0.d0)then
call errwarn('READ_EXPOT','external potential periodic in direction wi&
&thnon-zero z component applied to 2D system periodic in xy plane (set&
& '//trim(i2s(xyzzyaaaa4))//').')
xyzzyaaab4=.true.
endif
endif
elseif(periodicity==1)then
if(xyzzyaaai1(xyzzyaaaa4)==3)then
call errwarn('READ_EXPOT','external potential periodic in y direction &
&applied to 1D system periodic in the x direction (set '//trim(i2s(xyz&
&zyaaaa4))//').')
xyzzyaaab4=.true.
elseif(xyzzyaaai1(xyzzyaaaa4)==4)then
call errwarn('READ_EXPOT','external potential periodic in z direction &
&applied to 1D system periodic in the x direction(set '//trim(i2s(xyzz&
&yaaaa4))//').')
xyzzyaaab4=.true.
elseif(xyzzyaaai1(xyzzyaaaa4)>4)then
if(xyzzyaaas1(1,xyzzyaaaa4)==0.d0.or.xyzzyaaas1(2,xyzzyaaaa4)/=0.d0.or&
&.xyzzyaaas1(3,xyzzyaaaa4)/=0.d0)then
call errwarn('READ_EXPOT','external potential periodic in direction wi&
&th non-zero y and/or z components applied to 1D system periodic in th&
&e x direction (set '//trim(i2s(xyzzyaaaa4))//').')
xyzzyaaab4=.true.
endif
endif
elseif(periodicity==0)then
call errwarn('READ_EXPOT','periodic external potential applied to fini&
&te system (set '//trim(i2s(xyzzyaaaa4))//').')
xyzzyaaab4=.true.
endif
if(xyzzyaaab4.and.xyzzyaaac4)call errstop('CHECK_EXPOT','Quitting.')
if(periodicity/=0)call xyzzyaadt1(xyzzyaaaa4)
case(14)
continue
case default
call errstop('CHECK_EXPOT','No such potential representation defined.'&
&)
end select
enddo
end subroutine xyzzyaada1
logical function expot_is_finite(rvec,ispin)
implicit none
integer,intent(in) :: ispin
real(dp),intent(in) :: rvec(3)
integer xyzzyaaaa5,xyzzyaaab5
real(dp) xyzzyaaac5
logical xyzzyaaad5
expot_is_finite=.true.
do xyzzyaaaa5=1,xyzzyaaae1
if(.not.any(xyzzyaaal1(:,xyzzyaaaa5)==ispin))cycle
xyzzyaaad5=.false.
xyzzyaaab5=xyzzyaaan1(xyzzyaaaa5)
if(xyzzyaaaj1(xyzzyaaaa5)==1)then
if(xyzzyaaay1(xyzzyaaab5))then
call xyzzyaadb1(xyzzyaaaa5,rvec,xyzzyaaac5,xyzzyaaad5)
endif
elseif(xyzzyaaaj1(xyzzyaaaa5)==7)then
if(xyzzyaabq1(xyzzyaaab5))then
call xyzzyaadh1(xyzzyaaaa5,rvec,xyzzyaaac5,xyzzyaaad5)
endif
endif
if(xyzzyaaad5)then
expot_is_finite=.false.
return
endif
enddo
end function expot_is_finite
subroutine eval_expot(rvec,ispin,pot,isinf)
implicit none
integer,intent(in) :: ispin
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: pot
logical,intent(out) :: isinf
integer xyzzyaaaa6
real(dp) xyzzyaaab6
logical xyzzyaaac6
pot=0.d0
isinf=.false.
do xyzzyaaaa6=1,xyzzyaaae1
if(.not.any(xyzzyaaal1(:,xyzzyaaaa6)==ispin))cycle
xyzzyaaac6=.false.
select case(xyzzyaaaj1(xyzzyaaaa6))
case(1)
call xyzzyaadb1(xyzzyaaaa6,rvec,xyzzyaaab6,xyzzyaaac6)
case(2)
call xyzzyaadc1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(3)
call xyzzyaadd1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(4)
call xyzzyaadf1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(5)
call xyzzyaadg1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(6)
xyzzyaaab6=xyzzyaabk1(xyzzyaaaa6)
case(7)
call xyzzyaadh1(xyzzyaaaa6,rvec,xyzzyaaab6,xyzzyaaac6)
case(8)
call xyzzyaadi1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(9)
call xyzzyaadj1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(10)
call xyzzyaadk1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(11)
call xyzzyaadm1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(12)
call xyzzyaadn1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(13)
call xyzzyaado1(xyzzyaaaa6,rvec,xyzzyaaab6)
case(14)
call xyzzyaadp1(rvec,xyzzyaaab6)
case default
call errstop('EVAL_EXPOT','No such potential representation defined.')
end select
if(xyzzyaaac6)then
xyzzyaaab6=1.d100
isinf=.true.
endif
pot=pot+xyzzyaaab6
enddo
end subroutine eval_expot
subroutine xyzzyaadb1(iset,rvec,val,infval)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
logical,intent(out) :: infval
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7
real(dp) xyzzyaaae7,xyzzyaaaf7,xyzzyaaag7,xyzzyaaah7(3)
logical xyzzyaaai7,xyzzyaaaj7
val=0.d0
xyzzyaaac7=xyzzyaaan1(iset)
xyzzyaaab7=xyzzyaaam1(iset)
infval=.false.
xyzzyaaai7=xyzzyaaaw1(xyzzyaaac7)
xyzzyaaaj7=xyzzyaaax1(xyzzyaaac7)
xyzzyaaae7=xyzzyaaav1(xyzzyaaac7)
xyzzyaaaf7=xyzzyaaau1(xyzzyaaac7)
xyzzyaaad7=xyzzyaaai1(iset)
if(.not.(xyzzyaaai7.or.xyzzyaaaj7))then
do xyzzyaaaa7=1,xyzzyaaab7
xyzzyaaah7=xyzzyaaar1(:,xyzzyaaaa7,iset)
xyzzyaaag7=abs(xyzzyaadq1(rvec,xyzzyaaah7,xyzzyaaad7,iset))
if(xyzzyaaag7<xyzzyaaae7)val=val+xyzzyaaaf7
enddo
else
do xyzzyaaaa7=1,xyzzyaaab7
xyzzyaaah7=xyzzyaaar1(:,xyzzyaaaa7,iset)
xyzzyaaag7=abs(xyzzyaadq1(rvec,xyzzyaaah7,xyzzyaaad7,iset))
if(xyzzyaaag7<xyzzyaaae7.and.xyzzyaaai7)infval=.true.
if(xyzzyaaag7>xyzzyaaae7.and.xyzzyaaaj7)infval=.true.
enddo
endif
end subroutine xyzzyaadb1
subroutine xyzzyaadc1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8
real(dp) xyzzyaaae8,xyzzyaaaf8(3)
val=0.d0
xyzzyaaac8=xyzzyaaan1(iset)
xyzzyaaab8=xyzzyaaam1(iset)
xyzzyaaad8=xyzzyaaai1(iset)
do xyzzyaaaa8=1,xyzzyaaab8
xyzzyaaaf8=xyzzyaaar1(:,xyzzyaaaa8,iset)
xyzzyaaae8=xyzzyaadq1(rvec,xyzzyaaaf8,xyzzyaaad8,iset)
val=val-xyzzyaaaz1(xyzzyaaac8)*xyzzyaaae8
enddo
end subroutine xyzzyaadc1
subroutine xyzzyaadd1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,xyzzyaaae9
real(dp) xyzzyaaaf9,xyzzyaaag9(3)
val=0.d0
xyzzyaaac9=xyzzyaaan1(iset)
xyzzyaaab9=xyzzyaaam1(iset)
xyzzyaaad9=xyzzyaaba1(xyzzyaaac9)
xyzzyaaae9=xyzzyaaai1(iset)
do xyzzyaaaa9=1,xyzzyaaab9
xyzzyaaag9=xyzzyaaar1(:,xyzzyaaaa9,iset)
xyzzyaaaf9=xyzzyaadr1(rvec,xyzzyaaag9,xyzzyaaae9,iset)
val=val+xyzzyaade1(xyzzyaaad9,xyzzyaaac9,xyzzyaaaf9)
enddo
end subroutine xyzzyaadd1
real(dp) function xyzzyaade1(itype,nset,r2)
implicit none
integer,intent(in) :: itype,nset
real(dp),intent(in) :: r2
real(dp) xyzzyaaaa10,xyzzyaaab10
xyzzyaaab10=0.d0
select case(itype)
case(1)
xyzzyaaab10=exp(-xyzzyaabc1(nset,1)*r2)
case(2)
xyzzyaaab10=r2
case(3)
xyzzyaaaa10=r2*xyzzyaabc1(nset,3)
if(xyzzyaaaa10>.25d0)then
xyzzyaaab10=sqrt(xyzzyaaaa10)
else
xyzzyaaab10=xyzzyaaaa10+.25d0
endif
case default
if(am_master)call errstop('EVAL_ANALYTIC_CORE','Unknown core function.&
&')
end select
xyzzyaade1=xyzzyaabd1(nset)+xyzzyaabb1(nset)*xyzzyaaab10
end function xyzzyaade1
subroutine xyzzyaadf1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11,xyzzyaaae11
real(dp) xyzzyaaaf11,xyzzyaaag11(3)
val=0.d0
xyzzyaaad11=xyzzyaaan1(iset)
xyzzyaaac11=xyzzyaaam1(iset)
xyzzyaaae11=xyzzyaaai1(iset)
do xyzzyaaaa11=1,xyzzyaaac11
xyzzyaaag11=xyzzyaaar1(:,xyzzyaaaa11,iset)
xyzzyaaaf11=xyzzyaadr1(rvec,xyzzyaaag11,xyzzyaaae11,iset)
do xyzzyaaab11=1,xyzzyaabe1(xyzzyaaad11)
val=val+xyzzyaabf1(xyzzyaaab11,xyzzyaaad11)*exp(-xyzzyaabg1(xyzzyaaab1&
&1,xyzzyaaad11)*xyzzyaaaf11)
enddo
enddo
end subroutine xyzzyaadf1
subroutine xyzzyaadg1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12
integer,parameter :: xyzzyaaaf12=5
real(dp) xyzzyaaag12,xyzzyaaah12,xyzzyaaai12,xyzzyaaaj12(3)
val=0.d0
xyzzyaaad12=xyzzyaaan1(iset)
xyzzyaaac12=xyzzyaaam1(iset)
xyzzyaaae12=xyzzyaaai1(iset)
do xyzzyaaaa12=1,xyzzyaaac12
xyzzyaaaj12=xyzzyaaar1(:,xyzzyaaaa12,iset)
xyzzyaaah12=xyzzyaadq1(rvec,xyzzyaaaj12,xyzzyaaae12,iset)
if(xyzzyaaah12<xyzzyaabi1(1,xyzzyaaad12))then
xyzzyaaag12=xyzzyaabj1(1,xyzzyaaad12)
elseif(xyzzyaaah12>xyzzyaabi1(xyzzyaabh1(xyzzyaaad12),xyzzyaaad12))the&
&n
xyzzyaaag12=xyzzyaabj1(xyzzyaabh1(xyzzyaaad12),xyzzyaaad12)
else
call lookup(xyzzyaabi1(1,xyzzyaaad12),xyzzyaabh1(xyzzyaaad12),xyzzyaac&
&z1(1),xyzzyaaab12)
xyzzyaaab12=min(max(xyzzyaaab12-(xyzzyaaaf12-1)/2,1),xyzzyaabh1(xyzzya&
&aad12)+1-xyzzyaaaf12)
call interp_nev(xyzzyaabi1(xyzzyaaab12,xyzzyaaad12),xyzzyaabj1(xyzzyaa&
&ab12,xyzzyaaad12),xyzzyaaaf12,xyzzyaacz1(1),xyzzyaaag12,xyzzyaaai12)
endif
val=val+xyzzyaaag12
enddo
end subroutine xyzzyaadg1
subroutine xyzzyaadh1(iset,rvec,val,infval)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
logical,intent(out) :: infval
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13
real(dp) xyzzyaaae13,xyzzyaaaf13,xyzzyaaag13,xyzzyaaah13,xyzzyaaai13(3&
&)
logical xyzzyaaaj13,xyzzyaaak13
val=0.d0
xyzzyaaac13=xyzzyaaan1(iset)
xyzzyaaab13=xyzzyaaam1(iset)
infval=.false.
xyzzyaaaj13=xyzzyaaaw1(xyzzyaaac13)
xyzzyaaak13=xyzzyaaax1(xyzzyaaac13)
xyzzyaaae13=xyzzyaaav1(xyzzyaaac13)
xyzzyaaaf13=xyzzyaaau1(xyzzyaaac13)
xyzzyaaag13=xyzzyaabl1(xyzzyaaac13)
xyzzyaaad13=xyzzyaaai1(iset)
if(.not.(xyzzyaaaj13.or.xyzzyaaak13))then
do xyzzyaaaa13=1,xyzzyaaab13
xyzzyaaai13=xyzzyaaar1(:,xyzzyaaaa13,iset)
xyzzyaaah13=abs(xyzzyaads1(rvec,xyzzyaaai13,xyzzyaaag13,xyzzyaaad13,is&
&et))
if(xyzzyaaah13<xyzzyaaae13)val=val+xyzzyaaaf13
enddo
else
do xyzzyaaaa13=1,xyzzyaaab13
xyzzyaaai13=xyzzyaaar1(:,xyzzyaaaa13,iset)
xyzzyaaah13=abs(xyzzyaads1(rvec,xyzzyaaai13,xyzzyaaag13,xyzzyaaad13,is&
&et))
if(xyzzyaaah13<xyzzyaaae13.and.xyzzyaaaj13)infval=.true.
if(xyzzyaaah13>xyzzyaaae13.and.xyzzyaaak13)infval=.true.
enddo
endif
end subroutine xyzzyaadh1
subroutine xyzzyaadi1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14,xyzzyaaad14
real(dp) xyzzyaaae14,xyzzyaaaf14(3),xyzzyaaag14
val=0.d0
xyzzyaaac14=xyzzyaaan1(iset)
xyzzyaaab14=xyzzyaaam1(iset)
xyzzyaaag14=xyzzyaabr1(xyzzyaaac14)
xyzzyaaad14=xyzzyaaai1(iset)
do xyzzyaaaa14=1,xyzzyaaab14
xyzzyaaaf14=xyzzyaaar1(:,xyzzyaaaa14,iset)
xyzzyaaae14=xyzzyaads1(rvec,xyzzyaaaf14,xyzzyaaag14,xyzzyaaad14,iset)
val=val-xyzzyaabs1(xyzzyaaac14)*xyzzyaaae14
enddo
end subroutine xyzzyaadi1
subroutine xyzzyaadj1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15,xyzzyaaad15
real(dp) xyzzyaaae15,xyzzyaaaf15(3)
val=0.d0
xyzzyaaac15=xyzzyaaan1(iset)
xyzzyaaab15=xyzzyaaam1(iset)
xyzzyaaad15=xyzzyaaai1(iset)
do xyzzyaaaa15=1,xyzzyaaab15
xyzzyaaaf15=xyzzyaaar1(:,xyzzyaaaa15,iset)
xyzzyaaae15=xyzzyaadq1(rvec,xyzzyaaaf15,xyzzyaaad15,iset)
val=val+xyzzyaabv1(xyzzyaaac15)*cos(xyzzyaaae15*xyzzyaabu1(xyzzyaaac15&
&))
enddo
end subroutine xyzzyaadj1
subroutine xyzzyaadk1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16,xyzzyaaae16,xy&
&zzyaaaf16
real(dp) xyzzyaaag16,xyzzyaaah16,xyzzyaaai16,xyzzyaaaj16,xyzzyaaak16(3&
&)
val=0.d0
xyzzyaaad16=xyzzyaaan1(iset)
xyzzyaaac16=xyzzyaaam1(iset)
xyzzyaaae16=xyzzyaabw1(iset)
xyzzyaaaf16=xyzzyaaai1(iset)
xyzzyaaaj16=xyzzyaaby1(xyzzyaaad16)
do xyzzyaaaa16=1,xyzzyaaac16
xyzzyaaak16=xyzzyaaar1(:,xyzzyaaaa16,iset)
xyzzyaaai16=xyzzyaads1(rvec,xyzzyaaak16,xyzzyaaaj16,xyzzyaaaf16,iset)
val=val+xyzzyaacb1(xyzzyaaad16)
val=val+xyzzyaadl1(xyzzyaaae16,xyzzyaaad16,xyzzyaaai16)
xyzzyaaag16=xyzzyaaai16
xyzzyaaah16=xyzzyaaai16
do xyzzyaaab16=1,xyzzyaabx1(xyzzyaaad16)
xyzzyaaag16=xyzzyaaag16+xyzzyaaaj16
val=val+xyzzyaadl1(xyzzyaaae16,xyzzyaaad16,xyzzyaaag16)
xyzzyaaah16=xyzzyaaah16-xyzzyaaaj16
val=val+xyzzyaadl1(xyzzyaaae16,xyzzyaaad16,xyzzyaaah16)
enddo
enddo
end subroutine xyzzyaadk1
real(dp) function xyzzyaadl1(itype,nset,xyzzyaacz1)
implicit none
integer,intent(in) :: itype,nset
real(dp),intent(in) :: xyzzyaacz1
real(dp) xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17,xyzzyaaae17
xyzzyaadl1=0.d0
xyzzyaaae17=xyzzyaacz1*xyzzyaacz1
select case(itype)
case(1)
xyzzyaacy1=xyzzyaaca1(nset,1)*xyzzyaaae17
if(xyzzyaacy1<xyzzyaaap1)xyzzyaadl1=xyzzyaabz1(nset)*exp(-xyzzyaacy1)
case(2)
xyzzyaaaa17=xyzzyaaca1(nset,4)*xyzzyaaae17
xyzzyaaab17=xyzzyaaca1(nset,3)
if(xyzzyaaaa17>.25d0)then
xyzzyaaac17=sqrt(xyzzyaaaa17)
xyzzyaaad17=xyzzyaaab17**2-1.d0-3*(xyzzyaaab17-2*xyzzyaaac17)**2
else
xyzzyaaad17=(1.d0-xyzzyaaab17)*(2*xyzzyaaab17-1.d0-12*xyzzyaaaa17)
endif
xyzzyaadl1=xyzzyaaad17*xyzzyaabz1(nset)
case default
if(am_master)call errstop('EVAL_ANALYTIC_PER_CORE','Unknown core funct&
&ion.')
end select
end function xyzzyaadl1
subroutine xyzzyaadm1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xy&
&zzyaaaf18
real(dp) xyzzyaaag18,xyzzyaaah18,xyzzyaaai18(3)
val=0.d0
xyzzyaaad18=xyzzyaaan1(iset)
xyzzyaaac18=xyzzyaaam1(iset)
xyzzyaaae18=xyzzyaaai1(iset)
xyzzyaaag18=xyzzyaacg1(xyzzyaaad18)
xyzzyaaaf18=xyzzyaacd1(xyzzyaaad18)
do xyzzyaaaa18=1,xyzzyaaac18
xyzzyaaai18=xyzzyaaar1(:,xyzzyaaaa18,iset)
xyzzyaaah18=xyzzyaadq1(rvec,xyzzyaaai18,xyzzyaaae18,iset)
select case(xyzzyaace1(xyzzyaaad18))
case(1)
do xyzzyaaab18=1,xyzzyaaaf18
val=val+xyzzyaaci1(xyzzyaaab18,xyzzyaaad18)*sin(xyzzyaaag18*xyzzyaaab1&
&8*xyzzyaaah18)
enddo
case(2)
val=val+0.5d0*xyzzyaach1(0,xyzzyaaad18)
do xyzzyaaab18=1,xyzzyaaaf18
val=val+xyzzyaach1(xyzzyaaab18,xyzzyaaad18)*cos(xyzzyaaag18*xyzzyaaab1&
&8*xyzzyaaah18)
enddo
case(3)
val=val+0.5d0*xyzzyaach1(0,xyzzyaaad18)
do xyzzyaaab18=1,xyzzyaaaf18
val=val+xyzzyaach1(xyzzyaaab18,xyzzyaaad18)*cos(xyzzyaaag18*xyzzyaaab1&
&8*xyzzyaaah18)
val=val+xyzzyaaci1(xyzzyaaab18,xyzzyaaad18)*sin(xyzzyaaag18*xyzzyaaab1&
&8*xyzzyaaah18)
enddo
case default
call errstop('EVAL_EXPOT_FOURIER','Incorrect symmetry of Fourier exter&
&nal potential.')
end select
enddo
end subroutine xyzzyaadm1
subroutine xyzzyaadn1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,xy&
&zzyaaaf19
real(dp) xyzzyaaag19,xyzzyaaah19,xyzzyaaai19,xyzzyaaaj19,xyzzyaaak19,x&
&yzzyaaal19(3),xyzzyaaam19,xyzzyaaan19
val=0.d0
xyzzyaaae19=xyzzyaaan1(iset)
xyzzyaaad19=xyzzyaaam1(iset)
xyzzyaaaf19=xyzzyaaai1(iset)
xyzzyaaak19=xyzzyaaco1(xyzzyaaae19)
do xyzzyaaaa19=1,xyzzyaaad19
xyzzyaaal19=xyzzyaaar1(:,xyzzyaaaa19,iset)
xyzzyaaam19=xyzzyaads1(rvec,xyzzyaaal19,xyzzyaaak19,xyzzyaaaf19,iset)
xyzzyaaan19=xyzzyaaam19**2
do xyzzyaaac19=1,xyzzyaacm1(xyzzyaaae19)
val=val+xyzzyaacp1(xyzzyaaac19,xyzzyaaae19)*exp(-xyzzyaacq1(xyzzyaaac1&
&9,xyzzyaaae19)*xyzzyaaan19)
enddo
xyzzyaaag19=xyzzyaaam19
xyzzyaaah19=xyzzyaaam19
do xyzzyaaab19=1,xyzzyaacn1(xyzzyaaae19)
xyzzyaaag19=xyzzyaaag19+xyzzyaaak19
xyzzyaaai19=xyzzyaaag19*xyzzyaaag19
xyzzyaaah19=xyzzyaaah19-xyzzyaaak19
xyzzyaaaj19=xyzzyaaah19*xyzzyaaah19
do xyzzyaaac19=1,xyzzyaacm1(xyzzyaaae19)
xyzzyaacy1=xyzzyaacq1(xyzzyaaac19,xyzzyaaae19)*xyzzyaaai19
if(xyzzyaacy1<xyzzyaaap1)val=val+xyzzyaacp1(xyzzyaaac19,xyzzyaaae19)*e&
&xp(-xyzzyaacy1)
xyzzyaacy1=xyzzyaacq1(xyzzyaaac19,xyzzyaaae19)*xyzzyaaaj19
if(xyzzyaacy1<xyzzyaaap1)val=val+xyzzyaacp1(xyzzyaaac19,xyzzyaaae19)*e&
&xp(-xyzzyaacy1)
enddo
enddo
enddo
end subroutine xyzzyaadn1
subroutine xyzzyaado1(iset,rvec,val)
implicit none
integer,intent(in) :: iset
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa20,xyzzyaaab20,xyzzyaaac20,xyzzyaaad20,xyzzyaaae20
integer,parameter :: xyzzyaaaf20=5
real(dp) xyzzyaaag20,xyzzyaaah20,xyzzyaaai20,xyzzyaaaj20(3),xyzzyaaak2&
&0
val=0.d0
xyzzyaaad20=xyzzyaaan1(iset)
xyzzyaaac20=xyzzyaaam1(iset)
xyzzyaaae20=xyzzyaaai1(iset)
xyzzyaaai20=xyzzyaacs1(xyzzyaaad20)
do xyzzyaaaa20=1,xyzzyaaac20
xyzzyaaaj20=xyzzyaaar1(:,xyzzyaaaa20,iset)
xyzzyaaak20=xyzzyaads1(rvec,xyzzyaaaj20,xyzzyaaai20,xyzzyaaae20,iset)
call lookup(xyzzyaabi1(1,xyzzyaaad20),xyzzyaabh1(xyzzyaaad20),xyzzyaaa&
&k20,xyzzyaaab20)
xyzzyaaab20=min(max(xyzzyaaab20-(xyzzyaaaf20-1)/2,1),xyzzyaabh1(xyzzya&
&aad20)+1-xyzzyaaaf20)
call interp_nev(xyzzyaabi1(xyzzyaaab20,xyzzyaaad20),xyzzyaabj1(xyzzyaa&
&ab20,xyzzyaaad20),xyzzyaaaf20,xyzzyaaak20,xyzzyaaag20,xyzzyaaah20)
val=val+xyzzyaaag20
enddo
end subroutine xyzzyaado1
subroutine xyzzyaadp1(rvec,val)
implicit none
real(dp),intent(in) :: rvec(3)
real(dp),intent(out) :: val
integer xyzzyaaaa21(4,3),xyzzyaaab21,xyzzyaaac21,xyzzyaaad21,xyzzyaaae&
&21,xyzzyaaaf21
real(dp) xyzzyaaag21(3), xyzzyaaah21(3),xyzzyaaai21(3),xyzzyaaaj21(4),&
&xyzzyaaak21(4,3),xyzzyaaal21,xyzzyaaam21(64)
xyzzyaaah21=matmul(rvec(1:3),painv)
xyzzyaaag21=dble(xyzzyaacv1)
xyzzyaaaa21(2,1)=modulo(floor(xyzzyaaah21(1)*xyzzyaaag21(1)),xyzzyaacv&
&1(1))
xyzzyaaaa21(2,2)=modulo(floor(xyzzyaaah21(2)*xyzzyaaag21(2)),xyzzyaacv&
&1(2))
xyzzyaaaa21(2,3)=modulo(floor(xyzzyaaah21(3)*xyzzyaaag21(3)),xyzzyaacv&
&1(3))
do xyzzyaaad21=1,3
xyzzyaaaa21(1,xyzzyaaad21)=modulo(xyzzyaaaa21(2,xyzzyaaad21)-1,xyzzyaa&
&cv1(xyzzyaaad21))
xyzzyaaaa21(3,xyzzyaaad21)=modulo(xyzzyaaaa21(2,xyzzyaaad21)+1,xyzzyaa&
&cv1(xyzzyaaad21))
xyzzyaaaa21(4,xyzzyaaad21)=modulo(xyzzyaaaa21(2,xyzzyaaad21)+2,xyzzyaa&
&cv1(xyzzyaaad21))
enddo
xyzzyaaai21=modulo(xyzzyaaah21*xyzzyaaag21,xyzzyaaag21)
do xyzzyaaad21=1,3
xyzzyaaaj21(1)=xyzzyaaai21(xyzzyaaad21)-dble(xyzzyaaaa21(2,xyzzyaaad21&
&)-1)
xyzzyaaaj21(2)=xyzzyaaai21(xyzzyaaad21)-dble(xyzzyaaaa21(2,xyzzyaaad21&
&))
xyzzyaaaj21(3)=xyzzyaaai21(xyzzyaaad21)-dble(xyzzyaaaa21(2,xyzzyaaad21&
&)+1)
xyzzyaaaj21(4)=xyzzyaaai21(xyzzyaaad21)-dble(xyzzyaaaa21(2,xyzzyaaad21&
&)+2)
xyzzyaaak21(1,xyzzyaaad21)=2.d0+xyzzyaaaj21(1)*(-3.d0+xyzzyaaaj21(1)*(&
&1.5d0-0.25d0*xyzzyaaaj21(1)))
xyzzyaaak21(2,xyzzyaaad21)=1.d0+xyzzyaaaj21(2)*xyzzyaaaj21(2)*(-1.5d0+&
&0.75d0*xyzzyaaaj21(2))
xyzzyaaak21(3,xyzzyaaad21)=1.d0+xyzzyaaaj21(3)*xyzzyaaaj21(3)*(-1.5d0-&
&0.75d0*xyzzyaaaj21(3))
xyzzyaaak21(4,xyzzyaaad21)=2.d0+xyzzyaaaj21(4)*(3.d0+xyzzyaaaj21(4)*(1&
&.5d0+0.25d0*xyzzyaaaj21(4)))
enddo
xyzzyaaac21=0
do xyzzyaaaf21=1,4
do xyzzyaaae21=1,4
xyzzyaaal21=xyzzyaaak21(xyzzyaaae21,2)*xyzzyaaak21(xyzzyaaaf21,3)
do xyzzyaaad21=1,4
xyzzyaaac21=xyzzyaaac21+1
xyzzyaaam21(xyzzyaaac21)=xyzzyaaak21(xyzzyaaad21,1)*xyzzyaaal21
enddo
enddo
enddo
val=0.d0
xyzzyaaac21=0
if(xyzzyaaaa21(1,1)<xyzzyaacv1(1)-3)then
xyzzyaaab21=xyzzyaaaa21(1,1)-1
do xyzzyaaaf21=1,4
do xyzzyaaae21=1,4
do xyzzyaaad21=1,4
xyzzyaaac21=xyzzyaaac21+1
val=val+xyzzyaacw1(xyzzyaaab21+xyzzyaaad21,xyzzyaaaa21(xyzzyaaae21,2),&
&xyzzyaaaa21(xyzzyaaaf21,3))*xyzzyaaam21(xyzzyaaac21)
enddo
enddo
enddo
else
do xyzzyaaaf21=1,4
do xyzzyaaae21=1,4
do xyzzyaaad21=1,4
xyzzyaaac21=xyzzyaaac21+1
val=val+xyzzyaacw1(xyzzyaaaa21(xyzzyaaad21,1),xyzzyaaaa21(xyzzyaaae21,&
&2),xyzzyaaaa21(xyzzyaaaf21,3))*xyzzyaaam21(xyzzyaaac21)
enddo
enddo
enddo
endif
end subroutine xyzzyaadp1
real(dp) function xyzzyaadq1(rvec,ro,idir,iset)
implicit none
integer,intent(in) :: idir,iset
real(dp),intent(in) :: rvec(3),ro(3)
select case(idir)
case(1)
xyzzyaadq1=sqrt(sum((rvec-ro)**2))
case(2)
xyzzyaadq1=rvec(1)-ro(1)
case(3)
xyzzyaadq1=rvec(2)-ro(2)
case(4)
xyzzyaadq1=rvec(3)-ro(3)
case(5:8)
xyzzyaadq1=sum((rvec(:)-ro(:))*xyzzyaaas1(:,iset))
case default
xyzzyaadq1=0.d0
call errstop('R_ARG','Direction not defined.')
end select
end function xyzzyaadq1
real(dp) function xyzzyaadr1(rvec,ro,idir,iset)
implicit none
integer,intent(in) :: idir,iset
real(dp),intent(in) :: rvec(3),ro(3)
select case(idir)
case(1)
xyzzyaadr1=sum(rvec-ro)**2
case(2)
xyzzyaadr1=(rvec(1)-ro(1))**2
case(3)
xyzzyaadr1=(rvec(2)-ro(2))**2
case(4)
xyzzyaadr1=(rvec(3)-ro(3))**2
case(5:8)
xyzzyaadr1=sum((rvec(:)-ro(:))*xyzzyaaas1(:,iset))**2
case default
xyzzyaadr1=0.d0
call errstop('R2_ARG','Direction not defined.')
end select
end function xyzzyaadr1
real(dp) function xyzzyaads1(rvec,ro,d,idir,iset)
implicit none
integer,intent(in) :: idir,iset
real(dp),intent(in) :: rvec(3),ro(3),d
xyzzyaacx1=xyzzyaadq1(rvec,ro,idir,iset)
xyzzyaacy1=anint(xyzzyaacx1/d)*d
xyzzyaads1=xyzzyaacx1-xyzzyaacy1
end function xyzzyaads1
subroutine xyzzyaadt1(iset)
implicit none
integer,intent(in) :: iset
real(dp) xyzzyaaaa25,xyzzyaaab25
real(dp),parameter :: xyzzyaaac25=1.d-8
logical xyzzyaaad25,xyzzyaaae25
xyzzyaaad25=.false.
xyzzyaacz1(1:3)=0.d0
xyzzyaaae25=.false.
select case(xyzzyaaaj1(iset))
case(7)
call xyzzyaadh1(iset,xyzzyaacz1,xyzzyaaaa25,xyzzyaaae25)
case(8)
call xyzzyaadi1(iset,xyzzyaacz1,xyzzyaaaa25)
case(9)
call xyzzyaadj1(iset,xyzzyaacz1,xyzzyaaaa25)
case(10)
call xyzzyaadk1(iset,xyzzyaacz1,xyzzyaaaa25)
case(11)
call xyzzyaadm1(iset,xyzzyaacz1,xyzzyaaaa25)
case(12)
call xyzzyaadn1(iset,xyzzyaacz1,xyzzyaaaa25)
case(13)
call xyzzyaado1(iset,xyzzyaacz1,xyzzyaaaa25)
case default
call errstop('CHECK_COMMENSURATE','Called with inappropriate represent&
&ation type.')
end select
if(xyzzyaaae25)xyzzyaaaa25=1.d100
xyzzyaacz1(1:3)=a1(1:3)
xyzzyaaae25=.false.
select case(xyzzyaaaj1(iset))
case(7)
call xyzzyaadh1(iset,xyzzyaacz1,xyzzyaaab25,xyzzyaaae25)
case(8)
call xyzzyaadi1(iset,xyzzyaacz1,xyzzyaaab25)
case(9)
call xyzzyaadj1(iset,xyzzyaacz1,xyzzyaaab25)
case(10)
call xyzzyaadk1(iset,xyzzyaacz1,xyzzyaaab25)
case(11)
call xyzzyaadm1(iset,xyzzyaacz1,xyzzyaaab25)
case(12)
call xyzzyaadn1(iset,xyzzyaacz1,xyzzyaaab25)
case(13)
call xyzzyaado1(iset,xyzzyaacz1,xyzzyaaab25)
end select
if(xyzzyaaae25)xyzzyaaab25=1.d100
if(abs(xyzzyaaaa25-xyzzyaaab25)>xyzzyaaac25)xyzzyaaad25=.true.
if(periodicity>1)then
xyzzyaacz1(1:3)=a2(1:3)
xyzzyaaae25=.false.
select case(xyzzyaaaj1(iset))
case(7)
call xyzzyaadh1(iset,xyzzyaacz1,xyzzyaaab25,xyzzyaaae25)
case(8)
call xyzzyaadi1(iset,xyzzyaacz1,xyzzyaaab25)
case(9)
call xyzzyaadj1(iset,xyzzyaacz1,xyzzyaaab25)
case(10)
call xyzzyaadk1(iset,xyzzyaacz1,xyzzyaaab25)
case(11)
call xyzzyaadm1(iset,xyzzyaacz1,xyzzyaaab25)
case(12)
call xyzzyaadn1(iset,xyzzyaacz1,xyzzyaaab25)
case(13)
call xyzzyaado1(iset,xyzzyaacz1,xyzzyaaab25)
end select
if(xyzzyaaae25)xyzzyaaab25=1.d100
if(abs(xyzzyaaaa25-xyzzyaaab25)>xyzzyaaac25)xyzzyaaad25=.true.
endif
if(periodicity>2)then
xyzzyaacz1(1:3)=a3(1:3)
xyzzyaaae25=.false.
select case(xyzzyaaaj1(iset))
case(7)
call xyzzyaadh1(iset,xyzzyaacz1,xyzzyaaab25,xyzzyaaae25)
case(8)
call xyzzyaadi1(iset,xyzzyaacz1,xyzzyaaab25)
case(9)
call xyzzyaadj1(iset,xyzzyaacz1,xyzzyaaab25)
case(10)
call xyzzyaadk1(iset,xyzzyaacz1,xyzzyaaab25)
case(11)
call xyzzyaadm1(iset,xyzzyaacz1,xyzzyaaab25)
case(12)
call xyzzyaadn1(iset,xyzzyaacz1,xyzzyaaab25)
case(13)
call xyzzyaado1(iset,xyzzyaacz1,xyzzyaaab25)
end select
if(xyzzyaaae25)xyzzyaaab25=1.d100
if(abs(xyzzyaaaa25-xyzzyaaab25)>xyzzyaaac25)xyzzyaaad25=.true.
endif
if(xyzzyaaad25)call errstop('CHECK_EXPOT','Periodic external potential&
& apparently not commensurate with the lattice.')
end subroutine xyzzyaadt1
end module slaarnaap

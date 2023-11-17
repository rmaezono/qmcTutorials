module slaarnaca
use dsp,only : dp
implicit none
public
integer,parameter :: npoly=5
integer,allocatable :: ncoeff(:),nlang(:),l_local(:),l_of_non_local(:,&
&:)
integer,allocatable :: nlrule1_ion(:)
real(dp) lcutofftol,nlcutofftol
real(dp),allocatable :: zion(:),atnum(:),rcut_loc(:),rcut_non_loc(:),f&
&coeff_loc(:,:),pp_radial_grid(:,:),fcoeff_non_loc(:,:,:)
integer,parameter :: maxl_cpp=2,numl_cpp=maxl_cpp+1
real(dp),allocatable :: cppalpha(:),cpprbaree(:),cpprbaree_sq(:),cpp_p&
&refac(:)
real(dp),allocatable :: cpprbarl(:,:),cpprbarl_sq(:,:)
real(dp),allocatable :: ion_fields(:,:),electron_fields(:,:)
logical have_ppots
logical have_ae
logical have_veep
logical,allocatable :: is_ae(:),is_cpp(:),haverbarl(:)
real(dp),allocatable :: fcoeff_non_loc_full(:,:,:,:),fcoeff_loc_full(:&
&,:,:)
contains
subroutine read_ppots
use slaarnaag,   only : periodic_table,periodic_table_nocap,htoev
use file_utils,  only : skip,open_units
use format_utils,only : wout,i2s,r2s,capitalize
use slaarnabg,    only : atno,nitype,nbasis,nitot,iontype,atom_basis_t&
&ype,model_system
use slaarnabt,   only : interp_nev,lookup,dcopy,dscal
use parallel,    only : am_master
use run_control, only : errstop_master,check_alloc
use store,       only : open_unit,forces
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,l,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaa&
&l2,xyzzyaaam2,xyzzyaaan2,xyzzyaaao2,xyzzyaaap2,xyzzyaaaq2,xyzzyaaar2,&
&xyzzyaaas2
integer,parameter :: xyzzyaaat2=2
integer,parameter :: xyzzyaaau2=6
real(dp) xyzzyaaav2,xyzzyaaaw2,xyzzyaaax2,xyzzyaaay2,xyzzyaaaz2,xyzzya&
&aba2
real(dp),parameter :: xyzzyaabb2=1.d-14
real(dp),allocatable :: xyzzyaabc2(:,:)
logical xyzzyaabd2
logical,allocatable :: xyzzyaabe2(:)
character(20) psp_filename,psp_filename_compress,ppchar
character(80) char80,tmpr
character(1),parameter :: xyzzyaabf2(0:xyzzyaaau2)=(/'s','p','d','f','&
&g','h','i'/)
have_ppots=.false.
have_veep=.false.
call open_units(xyzzyaaao2,xyzzyaaal2)
if(xyzzyaaal2/=0)call errstop_master('READ_PPOTS','Unable to find free&
& i/o unit')
if(am_master.and.xyzzyaaat2>=2.and..not.model_system)then
call wout('Pseudopotentials')
call wout('================')
endif
allocate(ncoeff(nitype),nlang(nitype),l_local(nitype),zion(nitype),atn&
&um(nitype),rcut_loc(nitype),rcut_non_loc(nitype),is_cpp(nitype),haver&
&barl(nitype),cppalpha(nitype),cpprbaree(nitype),cpprbarl(numl_cpp,nit&
&ype),cpprbaree_sq(nitype),cpprbarl_sq(numl_cpp,nitype),cpp_prefac(nit&
&ot),xyzzyaabe2(nitype),is_ae(nitot),stat=xyzzyaaak2)
call check_alloc(xyzzyaaak2,'READ_PPOTS','1')
do xyzzyaaae2=1,nitype
do xyzzyaaag2=1,nbasis
if(iontype(xyzzyaaag2)==xyzzyaaae2)exit
enddo
if(model_system)then
xyzzyaabd2=.false.
xyzzyaaaq2=atno(xyzzyaaag2)
xyzzyaaar2=0
else
xyzzyaaaq2=mod(atno(xyzzyaaag2),1000)
xyzzyaaar2=atno(xyzzyaaag2)/1000+1
ppchar=''
if(xyzzyaaar2>1)ppchar=trim(i2s(xyzzyaaar2))
psp_filename=trim(periodic_table_nocap(xyzzyaaaq2))//trim(adjustl(ppch&
&ar))//'_pp.data'
inquire(file=psp_filename,exist=xyzzyaabd2,err=50)
endif
if(.not.xyzzyaabd2)then
if(xyzzyaaar2>1)call errstop_master('READ_PPOTS','Cannot find addition&
&al pseudopotential file '//trim(psp_filename)//' required for atomic &
&number '//trim(i2s(atno(xyzzyaaag2)))//' specified in xwfn.data file.&
&')
select case(trim(atom_basis_type))
case('gaussian')
if(.not.have_ae)call errstop_master('READ_PPOTS','Pseudopotential file&
& '//trim(psp_filename)//' not found.')
case('plane-wave','blip')
psp_filename_compress=trim(psp_filename)//'.gz'
inquire(file=psp_filename_compress,exist=xyzzyaabd2,err=51)
if(xyzzyaabd2)then
call errstop_master('READ_PPOTS','Pseudopotential file '//trim(psp_fil&
&ename)//' appears to be gzip compressed. Please gunzip it, or, if you&
& wish to run in all-electron mode, remove or rename the gzipped file &
&as well.')
else
psp_filename_compress=trim(psp_filename)//'.bz2'
inquire(file=psp_filename_compress,exist=xyzzyaabd2,err=51)
if(xyzzyaabd2)then
call errstop_master('READ_PPOTS','Pseudopotential file '//trim(psp_fil&
&ename)//' appears to be bzip2 compressed. Please uncompress it or, if&
& you wish to run in all-electron mode, remove or rename the compresse&
&d file as well.')
else
have_ae=.true.
endif
endif
case('numerical','dimer')
have_ae=.true.
case default
continue
end select
ncoeff(xyzzyaaae2)=0
nlang(xyzzyaaae2)=0
zion(xyzzyaaae2)=xyzzyaaaq2
atnum(xyzzyaaae2)=xyzzyaaaq2
rcut_loc(xyzzyaaae2)=0.d0
rcut_non_loc(xyzzyaaae2)=0.d0
xyzzyaabe2(xyzzyaaae2)=.false.
is_cpp(xyzzyaaae2)=.false.
cppalpha(xyzzyaaae2)=0.d0
cpprbaree(xyzzyaaae2)=0.d0
if(am_master.and.xyzzyaaat2>=2.and..not.model_system)call wout('No pse&
&udopotential for '//trim(adjustl(periodic_table(xyzzyaaaq2)))//' - tr&
&eating as all-electron.')
else
open(xyzzyaaao2,file=psp_filename,status='old',action='read',iostat=xy&
&zzyaaal2)
if(xyzzyaaal2/=0)call errstop_master('READ_PPOTS','Error opening pseud&
&opotential file '//trim(psp_filename)//'.')
read(xyzzyaaao2,*,err=10,end=15)
call skip(xyzzyaaao2,9)
read(xyzzyaaao2,*,err=10,end=15)ncoeff(xyzzyaaae2)
call skip(xyzzyaaao2,1)
nlang(xyzzyaaae2)=-1
read_file : do
do xyzzyaaaa2=1,ncoeff(xyzzyaaae2)
read(xyzzyaaao2,*,iostat=xyzzyaaal2)
if(xyzzyaaal2/=0)exit read_file
enddo
read(xyzzyaaao2,'(a)',iostat=xyzzyaaal2)char80
if(xyzzyaaal2/=0)exit
char80=adjustl(char80)
call capitalize(char80,.true.)
if(char80(1:2)/='r*'.or.index(char80,'potential')==0 .or.index(char80,&
&' in ')==0)exit
nlang(xyzzyaaae2)=nlang(xyzzyaaae2)+1
if(nlang(xyzzyaaae2)>xyzzyaaau2)call errstop_master('READ_PPOTS','Plea&
&se extend the list of spectroscopic names for l=0,1,... at the top of&
& ppots.f90.')
if(index(char80,'l='//trim(i2s(nlang(xyzzyaaae2))))==0.and.index(char8&
&0,'r*'//xyzzyaabf2(nlang(xyzzyaaae2)))==0)call errstop_master('READ_P&
&POTS','Error reading '//trim(psp_filename)//': angular momentum compo&
&nents not in order of increasing l.')
enddo read_file
if(nlang(xyzzyaaae2)<0)call errstop_master('READ_PPOTS','No potential &
&data in '//trim(psp_filename)//'?  Please ensure that this file is in&
& the format described in the CASINO manual.  In particular, the angul&
&ar-momentum components should be headed by e.g. "r*potential (L=0) in&
& Ry".')
close(xyzzyaaao2)
xyzzyaabe2(xyzzyaaae2)=.true.
have_ppots=.true.
if(am_master.and.xyzzyaaat2>=4)then
call wout('Have found a pseudopotential file for '//trim(periodic_tabl&
&e(xyzzyaaaq2))//'.')
call wout('  Number of radial grid points in file         :  ' //trim(&
&i2s(ncoeff(xyzzyaaae2))))
call wout('  Highest-l angular-momentum component in file :  ' //trim(&
&i2s(nlang(xyzzyaaae2))))
endif
endif
enddo
if(have_ppots)then
xyzzyaaai2=maxval(ncoeff(1:nitype))+1
xyzzyaaap2=maxval(nlang)
allocate(pp_radial_grid(xyzzyaaai2,nitype),fcoeff_loc(xyzzyaaai2,nityp&
&e),fcoeff_non_loc(xyzzyaaai2,0:xyzzyaaap2,nitype),xyzzyaabc2(xyzzyaaa&
&i2,0:xyzzyaaap2),l_of_non_local(xyzzyaaap2,nitype),stat=xyzzyaaak2)
call check_alloc(xyzzyaaak2,'READ_PPOTS','2')
xyzzyaabc2=0.d0
if(forces) then
allocate(fcoeff_non_loc_full(xyzzyaaai2,0:xyzzyaaap2,0:xyzzyaaap2,nity&
&pe),fcoeff_loc_full(xyzzyaaai2,0:xyzzyaaap2,nitype),stat=xyzzyaaak2)
call check_alloc(xyzzyaaak2,'READ_PPOTS','2.5')
endif
do xyzzyaaae2=1,nitype
do xyzzyaaag2=1,nbasis
if(iontype(xyzzyaaag2)==xyzzyaaae2)exit
enddo
xyzzyaaar2=atno(xyzzyaaag2)/1000+1
atno(xyzzyaaag2)=mod(atno(xyzzyaaag2),1000)
ppchar=''
if(xyzzyaaar2>1)ppchar=trim(i2s(xyzzyaaar2))
psp_filename=trim(periodic_table_nocap(atno(xyzzyaaag2)))//trim(adjust&
&l(ppchar))//'_pp.data'
if(ncoeff(xyzzyaaae2)==0)then
fcoeff_loc(:,xyzzyaaae2)=0.d0
fcoeff_non_loc(:,:,xyzzyaaae2)=0.d0
if(am_master.and.xyzzyaaat2>=2)then
call wout()
call wout('Ion type '//trim(i2s(xyzzyaaae2))//', atomic no '//trim(i2s&
&(atno(xyzzyaaag2)))//', file='//trim(psp_filename))
call wout('Pseudopotential file does not exist.')
call wout()
call wout('This ion type will be treated as all-electron.')
tmpr=r2s(zion(xyzzyaaae2),'(f21.12)')
call wout('Ionic charge                            : '//trim(tmpr))
endif
else
if(am_master.and.xyzzyaaat2>=2)then
call wout()
call wout('Ion type '//trim(i2s(xyzzyaaae2))//', atomic no '//trim(i2s&
&(atno(xyzzyaaag2)))//', file='//trim(psp_filename))
endif
open(xyzzyaaao2,file=psp_filename,status='old',action='read',iostat=xy&
&zzyaaal2)
if(xyzzyaaal2/=0)call errstop_master('READ_PPOTS','Error opening '//tr&
&im(psp_filename)//' file.')
read(xyzzyaaao2,'(a80)',err=10,end=15)char80
if(am_master.and.xyzzyaaat2>=2)then
call wout('Title:  '//trim(adjustl(char80)))
call wout()
endif
read(xyzzyaaao2,*,err=10,end=15)
read(xyzzyaaao2,*,err=100,end=100)xyzzyaaaf2,xyzzyaaav2
if(xyzzyaaaf2/=atno(xyzzyaaag2).and.am_master)call errstop_master('REA&
&D_PPOTS','Atomic number in '//trim(psp_filename)//' does not agree wi&
&th value in xwfn.data.')
atnum(xyzzyaaae2)=xyzzyaaaf2
zion(xyzzyaaae2)=xyzzyaaav2
read(xyzzyaaao2,*,err=10,end=15)
read(xyzzyaaao2,'(a80)',err=10,end=15)char80
char80=adjustl(char80)
call capitalize(char80,.true.)
if(trim(char80)=='ev')then
xyzzyaaaz2=1.d0/htoev
elseif(trim(char80)=='hartree')then
xyzzyaaaz2=1.d0
elseif(trim(char80)=='rydberg')then
xyzzyaaaz2=0.5d0
else
call errstop_master('READ_PPOTS','Unrecognized energy unit in '//trim(&
&psp_filename)//'.')
endif
read(xyzzyaaao2,*,err=10,end=15)
read(xyzzyaaao2,'(a)',iostat=xyzzyaaal2)char80
if(xyzzyaaal2/=0)call errstop_master('READ_PPOTS','Error reading L_LOC&
&AL in '//trim(psp_filename)//'.')
read(char80,*,iostat=xyzzyaaal2)xyzzyaaas2,l_local(xyzzyaaae2)
if(xyzzyaaal2/=0)then
read(char80,*,iostat=xyzzyaaal2)l_local(xyzzyaaae2)
if(xyzzyaaal2/=0)call errstop_master('READ_PPOTS','Error reading L_LOC&
&AL in '//trim(psp_filename)//'.')
endif
if(l_local(xyzzyaaae2)<0.or.l_local(xyzzyaaae2)>nlang(xyzzyaaae2))call&
& errstop_master('READ_PPOTS','Invalid value of l_local in '//trim(psp&
&_filename)//'.')
read(xyzzyaaao2,*,err=10,end=15)
read(xyzzyaaao2,*,iostat=xyzzyaaal2)xyzzyaaaj2
if(xyzzyaaal2/=0)call errstop_master('READ_PPOTS','Error reading NLRUL&
&E overrides in '//trim(psp_filename)//'.')
if(xyzzyaaaj2/=0)nlrule1_ion(xyzzyaaae2)=xyzzyaaaj2
read(xyzzyaaao2,*,err=10,end=15)
read(xyzzyaaao2,*,iostat=xyzzyaaal2)ncoeff(xyzzyaaae2)
if(xyzzyaaal2/=0)call errstop_master('READ_PPOTS','Error reading numbe&
&r of grid points in '//trim(psp_filename)//'.')
if(am_master.and.xyzzyaaat2>=2)call wout('Number of grid points       &
&            :  '//trim(i2s(ncoeff(xyzzyaaae2))))
if(ncoeff(xyzzyaaae2)>xyzzyaaai2)call errstop_master('READ_PPOTS','Too&
& many grid points in '//trim(psp_filename)//'. This should not happen&
&.')
if(ncoeff(xyzzyaaae2)<npoly)call errstop_master('READ_PPOTS','Too few &
&grid points in '//trim(psp_filename)//' file.')
read(xyzzyaaao2,*,err=10,end=15)
read(xyzzyaaao2,*,err=600,end=600)xyzzyaaba2
if(xyzzyaaba2<0.d0)then
call errstop_master('READ_PPOTS','Radial grid should not be negative i&
&n '//trim(psp_filename)//'.')
elseif(xyzzyaaba2==0.d0)then
xyzzyaaad2=1
else
xyzzyaaad2=2
pp_radial_grid(2,xyzzyaaae2)=xyzzyaaba2
ncoeff(xyzzyaaae2)=ncoeff(xyzzyaaae2)+1
if(am_master.and.xyzzyaaat2>=5)call wout('r=0 grid point is missing in&
& '//trim(psp_filename)//'.  It has now been added.')
endif
pp_radial_grid(1,xyzzyaaae2)=0.d0
do xyzzyaaaa2=xyzzyaaad2+1,ncoeff(xyzzyaaae2)
read(xyzzyaaao2,*,err=600,end=600)pp_radial_grid(xyzzyaaaa2,xyzzyaaae2&
&)
if(pp_radial_grid(xyzzyaaaa2,xyzzyaaae2)<=pp_radial_grid(xyzzyaaaa2-1,&
&xyzzyaaae2))call errstop_master('READ_PPOTS','Radial grid should incr&
&ease monotonically in '//trim(psp_filename)//'.')
enddo
xyzzyaaaa2=1
do xyzzyaaam2=0,nlang(xyzzyaaae2)
if(xyzzyaaam2/=l_local(xyzzyaaae2))then
l_of_non_local(xyzzyaaaa2,xyzzyaaae2)=xyzzyaaam2
xyzzyaaaa2=xyzzyaaaa2+1
endif
enddo
do xyzzyaaam2=0,nlang(xyzzyaaae2)
read(xyzzyaaao2,*,err=10,end=15)
do xyzzyaaaa2=xyzzyaaad2,ncoeff(xyzzyaaae2)
read(xyzzyaaao2,*,iostat=xyzzyaaal2)xyzzyaabc2(xyzzyaaaa2,xyzzyaaam2)
if(xyzzyaaal2/=0)call errstop_master('READ_PPOTS','Error reading l='//&
&trim(i2s(xyzzyaaam2))//' component in '//trim(psp_filename)//'.')
if(xyzzyaaaa2>1)xyzzyaabc2(xyzzyaaaa2,xyzzyaaam2)=xyzzyaabc2(xyzzyaaaa&
&2,xyzzyaaam2)/pp_radial_grid(xyzzyaaaa2,xyzzyaaae2)
enddo
enddo
read(xyzzyaaao2,'(a)',iostat=xyzzyaaal2)char80
if(xyzzyaaal2/=0)char80=''
call capitalize(char80,.true.)
if(index(char80,'core polar')==0)then
is_cpp(xyzzyaaae2)=.false.
else
read(xyzzyaaao2,*,iostat=xyzzyaaal2)cppalpha(xyzzyaaae2),cpprbaree(xyz&
&zyaaae2)
if(xyzzyaaal2/=0)then
is_cpp(xyzzyaaae2)=.false.
else
cpprbaree_sq(xyzzyaaae2)=cpprbaree(xyzzyaaae2)*cpprbaree(xyzzyaaae2)
is_cpp(xyzzyaaae2)=.true.
haverbarl(xyzzyaaae2)=.false.
read(xyzzyaaao2,*,iostat=xyzzyaaal2)cpprbarl(1:numl_cpp,xyzzyaaae2)
if(xyzzyaaal2==0)then
cpprbarl_sq(1:numl_cpp,xyzzyaaae2)=cpprbarl(1:numl_cpp,xyzzyaaae2)*cpp&
&rbarl(1:numl_cpp,xyzzyaaae2)
haverbarl(xyzzyaaae2)=.true.
endif
endif
endif
xyzzyaabc2(1,0:nlang(xyzzyaaae2))=xyzzyaabc2(2,0:nlang(xyzzyaaae2))+(x&
&yzzyaabc2(3,0:nlang(xyzzyaaae2))-xyzzyaabc2(2,0:nlang(xyzzyaaae2)))*(&
&pp_radial_grid(2,xyzzyaaae2)/(pp_radial_grid(2,xyzzyaaae2)-pp_radial_&
&grid(3,xyzzyaaae2)))
if(trim(atom_basis_type)=='gaussian'.and.all(atno(:)==1).and.all(xyzzy&
&aabc2(:,:)==0.d0))model_system=.true.
do xyzzyaaam2=0,nlang(xyzzyaaae2)
call dscal(ncoeff(xyzzyaaae2),xyzzyaaaz2,xyzzyaabc2(1,xyzzyaaam2),1)
enddo
call dcopy(ncoeff(xyzzyaaae2),xyzzyaabc2(1,l_local(xyzzyaaae2)),1,fcoe&
&ff_loc(1,xyzzyaaae2),1)
do xyzzyaaam2=0,nlang(xyzzyaaae2)
if(xyzzyaaam2/=l_local(xyzzyaaae2))fcoeff_non_loc(1:ncoeff(xyzzyaaae2)&
&,xyzzyaaam2,xyzzyaaae2)=xyzzyaabc2(1:ncoeff(xyzzyaaae2),xyzzyaaam2)-x&
&yzzyaabc2(1:ncoeff(xyzzyaaae2),l_local(xyzzyaaae2))
enddo
do xyzzyaaam2=nlang(xyzzyaaae2),0,-1
if(xyzzyaaam2==l_local(xyzzyaaae2))exit
if(xyzzyaaam2==2.and.forces)exit
if(all(abs(fcoeff_non_loc(1:ncoeff(xyzzyaaae2),xyzzyaaam2,xyzzyaaae2))&
&<xyzzyaabb2))then
nlang(xyzzyaaae2)=nlang(xyzzyaaae2)-1
if(xyzzyaaat2>=2.and.am_master)call wout('Discarding unnecessary l='//&
&trim(i2s(xyzzyaaam2))//' angular-momentum component.')
else
exit
endif
enddo
if(nlang(xyzzyaaae2)<0)call errstop_master('READ_PPOTS','Bug.  nlang i&
&s too small.')
close(xyzzyaaao2)
open_unit(xyzzyaaao2)=.false.
if(forces)then
do xyzzyaaan2=0,nlang(xyzzyaaae2)
call dcopy(ncoeff(xyzzyaaae2),xyzzyaabc2(1,xyzzyaaan2),1,fcoeff_loc_fu&
&ll(1,xyzzyaaan2,xyzzyaaae2),1)
do xyzzyaaam2=0,nlang(xyzzyaaae2)
if(xyzzyaaam2/=xyzzyaaan2)fcoeff_non_loc_full(1:ncoeff(xyzzyaaae2),xyz&
&zyaaam2,xyzzyaaan2,xyzzyaaae2) =xyzzyaabc2(1:ncoeff(xyzzyaaae2),xyzzy&
&aaam2)-xyzzyaabc2(1:ncoeff(xyzzyaaae2),xyzzyaaan2)
enddo
enddo
endif
do xyzzyaaaa2=ncoeff(xyzzyaaae2),2,-1
if(abs(fcoeff_loc(xyzzyaaaa2,xyzzyaaae2)+zion(xyzzyaaae2)/pp_radial_gr&
&id(xyzzyaaaa2,xyzzyaaae2))>lcutofftol)exit
enddo
if(xyzzyaaaa2<ncoeff(xyzzyaaae2))xyzzyaaaa2=xyzzyaaaa2+1
rcut_loc(xyzzyaaae2)=pp_radial_grid(xyzzyaaaa2,xyzzyaaae2)
do xyzzyaaab2=ncoeff(xyzzyaaae2),2,-1
if(any(abs(fcoeff_non_loc(xyzzyaaab2,l_of_non_local(1:nlang(xyzzyaaae2&
&),xyzzyaaae2),xyzzyaaae2)) >nlcutofftol))exit
enddo
if(xyzzyaaab2<ncoeff(xyzzyaaae2))xyzzyaaab2=xyzzyaaab2+1
rcut_non_loc(xyzzyaaae2)=pp_radial_grid(xyzzyaaab2,xyzzyaaae2)
if(am_master.and.xyzzyaaat2>=2)then
tmpr=r2s(zion(xyzzyaaae2),'(f4.1)')
call wout('Pseudo-charge                           :  '//trim(tmpr))
tmpr=r2s(rcut_loc(xyzzyaaae2),'(f10.6)')
call wout('Local cutoff radius (au)                :  '//trim(tmpr))
tmpr=r2s(rcut_non_loc(xyzzyaaae2),'(f10.6)')
call wout('Non-local cutoff radius (au)            :  '//trim(tmpr))
if(xyzzyaaat2>=3)call wout('Highest-l angular-momentum component    : &
& '//trim(i2s(nlang(xyzzyaaae2))))
call lookup(pp_radial_grid(1,xyzzyaaae2),ncoeff(xyzzyaaae2),rcut_loc(x&
&yzzyaaae2),xyzzyaaah2)
xyzzyaaah2=min(max(xyzzyaaah2-(npoly-1)/2,1),ncoeff(xyzzyaaae2)+1-npol&
&y)
call interp_nev(pp_radial_grid(xyzzyaaah2,xyzzyaaae2),fcoeff_loc(xyzzy&
&aaah2,xyzzyaaae2),npoly,rcut_loc(xyzzyaaae2),xyzzyaaaw2,xyzzyaaay2)
xyzzyaaax2=(xyzzyaaaw2+zion(xyzzyaaae2)/rcut_loc(xyzzyaaae2))
call wout('Local potential angular momentum        :  ' //i2s(l_local(&
&xyzzyaaae2)))
tmpr=r2s(xyzzyaaax2,'(e16.10)')
call wout('Deviation from z/r at cutoff radius     :  '//trim(tmpr))
call lookup(pp_radial_grid(1,xyzzyaaae2),ncoeff(xyzzyaaae2),rcut_non_l&
&oc(xyzzyaaae2),xyzzyaaah2)
xyzzyaaah2=min(max(xyzzyaaah2-(npoly-1)/2,1),ncoeff(xyzzyaaae2)+1-npol&
&y)
do xyzzyaaac2=1,nlang(xyzzyaaae2)
l=l_of_non_local(xyzzyaaac2,xyzzyaaae2)
call interp_nev(pp_radial_grid(xyzzyaaah2,xyzzyaaae2),fcoeff_non_loc(x&
&yzzyaaah2,l,xyzzyaaae2),npoly,rcut_non_loc(xyzzyaaae2),xyzzyaaax2,xyz&
&zyaaay2)
call wout('Non-local potential angular momentum    :  '//trim(i2s(l)))
tmpr=r2s(abs(xyzzyaaax2),'(e16.10)')
call wout('Deviation from vlocal at cutoff radius  :  '//trim(tmpr))
enddo
if(is_cpp(xyzzyaaae2))then
call wout('Core polarization potential - ')
tmpr=r2s(cppalpha(xyzzyaaae2),'(f21.12)')
call wout('alpha                                   : '//trim(tmpr))
tmpr=r2s(cpprbaree(xyzzyaaae2),'(f21.12)')
call wout('rbaree                                  : '//trim(tmpr))
if(haverbarl(xyzzyaaae2))then
do xyzzyaaam2=1,numl_cpp
tmpr=r2s(cpprbarl(xyzzyaaam2,xyzzyaaae2),'(f21.12)')
call wout('rbar('//trim(i2s(xyzzyaaam2)) //')                         &
&        : '//trim(tmpr))
enddo
else
call wout('No rbar(l) present')
endif
endif
endif
endif
enddo
endif
if(any(is_cpp))then
have_veep=.true.
do xyzzyaaaa2=1,nitot
if(is_cpp(iontype(xyzzyaaaa2)))then
cpp_prefac(xyzzyaaaa2)=-0.5d0*cppalpha(iontype(xyzzyaaaa2))
else
cpp_prefac(xyzzyaaaa2)=0.d0
endif
enddo
else
cpp_prefac=0.d0
endif
do xyzzyaaaa2=1,nitot
is_ae(xyzzyaaaa2)=.not.xyzzyaabe2(iontype(xyzzyaaaa2))
enddo
deallocate(xyzzyaabe2)
atno(:)=mod(atno(:),1000)
if(am_master.and.xyzzyaaat2>=2)call wout()
return
10 call errstop_master('READ_PPOTS','Error reading '//trim(psp_filenam&
&e)//'.')
15 call errstop_master('READ_PPOTS','Reached end of '//trim(psp_filena&
&me)//'.')
50 call errstop_master('READ_PPOTS','Error inquiring about '//trim(psp&
&_filename)//'.')
51 call errstop_master('READ_PPOTS','Error inquiring about '//trim(psp&
&_filename_compress)//'.')
100 call errstop_master('READ_PPOTS','Error reading atno/zion in '//tr&
&im(psp_filename)//'.')
600 call errstop_master('READ_PPOTS','Error reading radial grid in '//&
&trim(psp_filename)//'.')
end subroutine read_ppots
end module slaarnaca

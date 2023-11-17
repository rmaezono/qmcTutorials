module slaarnaae
use dsp
use parallel
use slaarnaag,   only : czero,zi,one_over_twopi,pi_over_two,c_one
use file_utils,  only : open_units,skip
use format_utils,only : wout,r2s,i2s,l2s,wordwrap
use slaarnabg,    only : painv,periodicity,b1,b2,b3,pb1,pb2,pb3,basis,&
&atno,sr_lattice,pa1,pa2,pa3,num_g,a1,a2,a3,nbasis,rion,pbinv,isperiod&
&ic
use slaarnabp,        only : modified_mdet,wf_nd,wf_d,no_orb_phases,or&
&b_phase_band,orb_phase_kpoint,orb_phase_spin,orb_phase_det,orb_phase,&
&wf_nm,wf_np,broadcast_mdet,read_mdet_wfn
use run_control, only : errstop,errwarn,check_alloc,timer,errstop_mast&
&er
use shalloc,     only : shallocate_blip,deshallocate,shallocate_barrie&
&r,am_smpmaster,need_shm
use store,       only : nspin,nemax,ndet,open_unit,nuc_nele,no_ssingle&
&s,real1_complex2,complex_wf,bsmooth,localized_orbitals,orb_norm,use_b&
&ackflow
implicit none
private
public single_precision_blips,write_binary_blips,conv_binary_blips,rea&
&dbwf,blip_orb_eval,bwfdet_setup,deshalloc_bwfdet_shm,get_bwfdet_orbma&
&p,get_bwfdet_orbdesc,get_bwfdet_ndesc,get_bwfdet_rmax
integer xyzzyaaaa1(3)
real(dp) xyzzyaaab1(3),xyzzyaaac1(3)
integer xyzzyaaad1(2),xyzzyaaae1(2)
integer xyzzyaaaf1,xyzzyaaag1
integer,allocatable :: xyzzyaaah1(:,:)
integer,allocatable :: xyzzyaaai1(:,:,:,:)
real(dp) xyzzyaaaj1,xyzzyaaak1(6),xyzzyaaal1(6,6)
type ptr_sp
real(sp),pointer :: bcoefs(:,:,:)=>null()
end type ptr_sp
type ptr_dp
real(dp),pointer :: bcoefs(:,:,:)=>null()
end type ptr_dp
type(ptr_sp),allocatable :: xyzzyaaam1(:)
type(ptr_dp),allocatable :: xyzzyaaan1(:)
real(sp),pointer :: xyzzyaaao1(:,:,:,:)=>null()
real(dp),pointer :: xyzzyaaap1(:,:,:,:)=>null()
complex(sp),pointer :: xyzzyaaaq1(:,:,:,:)=>null()
complex(dp),pointer :: xyzzyaaar1(:,:,:,:)=>null()
integer xyzzyaaas1
integer xyzzyaaat1
real(dp),allocatable :: xyzzyaaau1(:,:)
complex(dp),allocatable :: xyzzyaaav1(:,:)
real(dp),allocatable :: xyzzyaaaw1(:,:,:,:)
real(dp),allocatable :: xyzzyaaax1(:,:,:)
real(dp),allocatable :: xyzzyaaay1(:,:),xyzzyaaaz1(:),xyzzyaaba1(:,:)
real(dp),allocatable,target :: xyzzyaabb1(:,:)
real(dp),allocatable :: xyzzyaabc1(:),xyzzyaabd1(:),xyzzyaabe1(:),xyzz&
&yaabf1(:)
integer,allocatable :: xyzzyaabg1(:,:),xyzzyaabh1(:),xyzzyaabi1(:,:)
logical xyzzyaabj1
logical xyzzyaabk1,xyzzyaabl1
logical xyzzyaabm1,xyzzyaabn1
logical xyzzyaabo1,single_precision_blips,write_binary_blips,conv_bina&
&ry_blips
logical,allocatable :: xyzzyaabp1(:)
logical,allocatable :: xyzzyaabq1(:)
real(dp),parameter :: xyzzyaabr1=-15.4420799369152d0,xyzzyaabs1=16.877&
&1030822967d0,xyzzyaabt1=-11.d0/6.d0,xyzzyaabu1=20.d0+10.d0*xyzzyaabt1&
&,xyzzyaabv1=-70.d0-34.d0*xyzzyaabt1,xyzzyaabw1=84.d0+39.d0*xyzzyaabt1&
&,xyzzyaabx1=-35.d0-15.d0*xyzzyaabt1
real(dp),parameter :: xyzzyaaby1(4:9)=(/xyzzyaabr1+xyzzyaabx1,xyzzyaab&
&s1-4.d0*xyzzyaabr1+xyzzyaabw1,6.d0*xyzzyaabr1-4.d0*xyzzyaabs1+xyzzyaa&
&bv1,6.d0*xyzzyaabs1-4.d0*xyzzyaabr1+xyzzyaabu1,xyzzyaabr1-4.d0*xyzzya&
&abs1,xyzzyaabs1/)
real(dp),allocatable :: xyzzyaabz1(:,:)
real(dp),allocatable :: xyzzyaaca1(:),xyzzyaacb1(:,:),xyzzyaacc1(:),xy&
&zzyaacd1(:,:)
complex(dp),allocatable :: xyzzyaace1(:),xyzzyaacf1(:,:),xyzzyaacg1(:)&
&,xyzzyaach1(:,:)
integer xyzzyaaci1,xyzzyaacj1,xyzzyaack1,xyzzyaacl1,bw_norbc
integer,allocatable :: xyzzyaacm1(:,:,:)
logical,allocatable :: xyzzyaacn1(:,:),xyzzyaaco1(:,:)
integer,allocatable :: xyzzyaacp1(:,:),xyzzyaacq1(:,:,:,:)
integer xyzzyaacr1,xyzzyaacs1,xyzzyaact1,xyzzyaacu1,xyzzyaacv1
integer,allocatable :: xyzzyaacw1(:)
integer,allocatable :: xyzzyaacx1(:),xyzzyaacy1(:),xyzzyaacz1(:)
integer,parameter :: xyzzyaada1=1,xyzzyaadb1=2,xyzzyaadc1=3
integer,allocatable :: xyzzyaadd1(:)
complex(dp),allocatable :: xyzzyaade1(:)
integer,parameter :: xyzzyaadf1=2,xyzzyaadg1=1
integer,allocatable :: xyzzyaadh1(:,:)
real(dp),allocatable :: xyzzyaadi1(:,:)
character(32),parameter :: xyzzyaadj1='CASINO v2.10-2'
contains
subroutine readbwf(eionion)
use slaarnaag,only : periodic_table_nocap
use slaarnaca,only : have_ae
implicit none
real(dp),intent(out) :: eionion
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,avc_norb,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,x&
&yzzyaaal2
integer,allocatable :: xyzzyaaam2(:),xyzzyaaan2(:)
integer xyzzyaaao2,xyzzyaaap2,xyzzyaaaq2
real(dp) xyzzyaaar2,xyzzyaaas2,xyzzyaaat2,xyzzyaaau2,xyzzyaaav2,xyzzya&
&aaw2,xyzzyaaax2
logical xyzzyaaay2,xyzzyaaaz2,xyzzyaaba2
logical xyzzyaabb2,xyzzyaabc2,xyzzyaabd2,xyzzyaabe2,xyzzyaabf2,xyzzyaa&
&bg2,xyzzyaabh2
logical,allocatable :: xyzzyaabi2(:)
character(20) psp_filename,ppchar
character(80) title,code,method,functional,pseudo_type,char80
character(80),parameter :: xyzzyaabj2='bwfn.data',xyzzyaabk2='bwfn.dat&
&a.bin',xyzzyaabl2='bwfn.data.b1',xyzzyaabm2='bwfn.data.b2',xyzzyaabn2&
&='bwfn.data.bin_mpi'
call timer("READBWF",.true.)
xyzzyaabo1=complex_wf
if(nspin/=2)call errstop_master('READBWF','Blip orbitals can only be u&
&sed for electron systems at present.')
inquire(file=trim(xyzzyaabj2),exist=xyzzyaabb2)
inquire(file=trim(xyzzyaabk2),exist=xyzzyaabc2)
inquire(file=trim(xyzzyaabl2),exist=xyzzyaabd2)
inquire(file=trim(xyzzyaabm2),exist=xyzzyaabe2)
inquire(file=trim(xyzzyaabn2),exist=xyzzyaabf2)
if(xyzzyaabf2)call errstop_master('READBWF','CASINO has dropped suppor&
&t for MPI-IO binary files.  Please remove '//trim(xyzzyaabn2)//' and &
&re-run this calculation.')
if(xyzzyaabe2)call errstop_master('READBWF','CASINO has dropped suppor&
&t for multiple bwfn.data.b* binary files. If you have a '//trim(xyzzy&
&aabj2)//' file, use that, else regenerate a single binary file.')
if(.not.(xyzzyaabb2.or.xyzzyaabc2.or.xyzzyaabd2))call errstop_master('&
&READBWF','Blip data file not found.  One of the following files shoul&
&d be present: '//trim(xyzzyaabj2)//', '//trim(xyzzyaabk2)//' or '//tr&
&im(xyzzyaabl2)//'.')
if(xyzzyaabb2.and.xyzzyaabc2)call errwarn('READBWF','Both '//trim(xyzz&
&yaabj2)//' and '//trim(xyzzyaabk2)//' are present.  Will read the lat&
&ter.')
if(xyzzyaabb2.and.xyzzyaabd2)call errwarn('READBWF','Both '//trim(xyzz&
&yaabj2)//' and '//trim(xyzzyaabl2)//' are present.  Will read the lat&
&ter.')
if(xyzzyaabd2.and.xyzzyaabc2)call errwarn('READBWF','Both '//trim(xyzz&
&yaabl2)//' and '//trim(xyzzyaabk2)//' are present.  Will read the lat&
&ter.')
call qmc_barrier
if(need_shm)then
xyzzyaaaa2=1
if(am_smpmaster)xyzzyaaaa2=0
call mpi_comm_split(mpi_comm_world,xyzzyaaaa2,my_node,xyzzyaaao2,ierro&
&r)
call checkmpi(ierror,'comm-splitting in bwfdet')
else
xyzzyaaao2=mpi_comm_world
endif
xyzzyaabg2=.false.
do
if(xyzzyaabc2)then
if(am_master)call wout('Reading binary data file '//trim(xyzzyaabk2)//&
&'.')
call timer('READ_BINARY',.true.)
call xyzzyaabr2(xyzzyaaba2)
call timer('READ_BINARY',.false.)
if(.not.xyzzyaaba2)then
call xyzzyaabt2
xyzzyaabc2=.false.
cycle
endif
call xyzzyaabv2
elseif(xyzzyaabd2)then
if(am_master)call wout('Reading old-format binary data file '//trim(xy&
&zzyaabl2)//'.')
call timer('READ_OLD_BINARY',.true.)
call xyzzyaabu2(xyzzyaaba2)
call timer('READ_OLD_BINARY',.false.)
if(.not.xyzzyaaba2)then
call xyzzyaabt2
xyzzyaabd2=.false.
cycle
endif
call xyzzyaabv2
if(conv_binary_blips.and.am_master)then
call open_units(xyzzyaaap2,xyzzyaaaq2)
if(xyzzyaaaq2/=0)call errstop('READBWF','Unable to find free i/o unit'&
&)
open(xyzzyaaap2,file=trim(xyzzyaabl2),form='unformatted',status='old',&
&action='read',iostat=xyzzyaaaq2)
if(xyzzyaaaq2==0)close(xyzzyaaap2,status='delete')
call wout('Writing binary file '//trim(xyzzyaabk2)//'.')
call xyzzyaabq2
endif
elseif(xyzzyaabb2)then
if(am_master)call wout('Reading plain data file '//trim(xyzzyaabj2)//'&
&.')
call timer('READ_FORMATTED',.true.)
call xyzzyaabp2
call timer('READ_FORMATTED',.false.)
if(write_binary_blips.and.am_master)then
call wout('Writing binary file '//trim(xyzzyaabk2)//'.')
call xyzzyaabq2
endif
else
call errstop_master('READBWF','Attempted file read(s) failed. If you h&
&ave a '//trim(xyzzyaabj2)//' file, use that.  Else, re-generate the b&
&inary file.')
endif
exit
enddo
if(am_master)call wout()
if(am_master)call xyzzyaabo2
if(xyzzyaabk1.and.(xyzzyaaaf1>1.or.any(xyzzyaaay1(1:3,1)/=0.d0)))call &
&errstop_master('READBWF','Real blip coefficients, but this is not jus&
&t a calculation at Gamma.')
if(.not.xyzzyaabk1.and.xyzzyaaaf1==1.and.all(xyzzyaaay1(1:3,1)==0.d0))&
&call errwarn('READBWF','Complex blip coefficients, but this is just a&
& calculation at Gamma.')
if(xyzzyaabo1)then
if(xyzzyaabm1)call errstop_master('READBWF','Localized orbitals should&
& not be used in conjunction with a complex wave function.')
if(xyzzyaabk1)call errstop_master('READBWF','If a complex wave functio&
&n is to be returned, the blip coefficients supplied should be complex&
&.')
if(xyzzyaacl1/=0)call errstop_master('READBWF','Complex blip wave func&
&tion, but real orbitals present.  It is likely that the value of COMP&
&LEX_WF in was changed in the input file but the bwfn.data.bin (or .b1&
&) file was not deleted.')
else
if(xyzzyaacl1/=xyzzyaack1)call errstop_master('READBWF','Real blip wav&
&e function, but real orbitals not present.  It is likely that the val&
&ue of COMPLEX_WF in was changed in the input file but the bwfn.data.b&
&in (or .b1) file was not deleted.')
endif
if(any(atno(:)>200))then
if(am_master)then
allocate(xyzzyaabi2(nbasis),xyzzyaaam2(nbasis),xyzzyaaan2(nbasis),stat&
&=xyzzyaaak2)
call check_alloc(xyzzyaaak2,'READBWF','is_pp_from_atno, etc..')
xyzzyaaam2=0
xyzzyaaan2=0
xyzzyaabi2(:)=atno(:)>200
do xyzzyaaai2=1,nbasis
atno(xyzzyaaai2)=mod(atno(xyzzyaaai2),100)
enddo
xyzzyaaaj2=0
have_ae=.false.
a:  do xyzzyaaai2=1,nbasis
do xyzzyaaac2=1,xyzzyaaai2-1
if(xyzzyaaam2(xyzzyaaac2)==atno(xyzzyaaai2))then
xyzzyaaam2(xyzzyaaai2)=xyzzyaaam2(xyzzyaaac2)
xyzzyaaan2(xyzzyaaai2)=xyzzyaaan2(xyzzyaaac2)
cycle a
endif
enddo
xyzzyaaaj2=xyzzyaaaj2+1
xyzzyaaam2(xyzzyaaai2)=atno(xyzzyaaai2)
xyzzyaaan2(xyzzyaaai2)=xyzzyaaaj2
enddo a
xyzzyaaah2=xyzzyaaaj2
do xyzzyaaal2=1,xyzzyaaah2
do xyzzyaaai2=1,nbasis
if(xyzzyaaan2(xyzzyaaai2)==xyzzyaaal2)exit
enddo
xyzzyaaaa2=mod(atno(xyzzyaaai2),1000)
xyzzyaaab2=atno(xyzzyaaai2)/1000+1
ppchar=''
if(xyzzyaaab2>1)ppchar=trim(i2s(xyzzyaaab2))
psp_filename=trim(periodic_table_nocap(xyzzyaaaa2))//trim(adjustl(ppch&
&ar))//'_pp.data'
inquire(file=psp_filename,exist=xyzzyaabh2,err=7)
if(xyzzyaabi2(xyzzyaaai2).and..not.xyzzyaabh2)call errstop('READBWF','&
&The atomic numbers in the blip file follow the +200 convention to fla&
&g pseudoatoms. However, the x_pp.data pseudopotential file for a flag&
&ged pseudoatom is missing.')
if(.not.xyzzyaabh2)have_ae=.true.
enddo
deallocate(xyzzyaabi2,xyzzyaaan2,xyzzyaaam2)
else
do xyzzyaaai2=1,nbasis
atno(xyzzyaaai2)=mod(atno(xyzzyaaai2),100)
enddo
endif
call mpi_bcast(have_ae,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting have_ae in readbwf.')
endif
call timer("READBWF",.false.)
return
7  call errstop('READBWF','Error inquiring about '//trim(psp_filename)&
&//'.')
contains
subroutine xyzzyaabo2
implicit none
character(80) tmpr
call wout('Title: '//trim(title))
call wout('Generating code                           : '//trim(code))
call wout('Method                                    : '//trim(method)&
&)
call wout('DFT functional                            : '//trim(functio&
&nal))
call wout('Pseudopotential type                      : '//trim(pseudo_&
&type))
tmpr=r2s(xyzzyaaau2,'(f12.3)')
call wout('Plane-wave cutoff (au)                    : '//trim(tmpr))
call wout()
call wout('Number of k points                        : '//trim(i2s(xyz&
&zyaaaf1)))
call wout('Max # bands per k point                   : '//trim(i2s(xyz&
&zyaaag1)))
call wout('Blip grid                                 : '//trim(i2s(xyz&
&zyaaaa1(1)))//' by '//trim(i2s(xyzzyaaaa1(2)))//' by '//trim(i2s(xyzz&
&yaaaa1(3))))
if(xyzzyaabl1)then
call wout('Spin-polarized orbital data.')
else
call wout('Non-spin-polarized orbital data.')
endif
call wout()
call wout('DFT energy and components (au per primitive cell):')
call wout('Total energy                              : ',xyzzyaaas2,ad&
&just=.true.)
call wout('Kinetic energy                            : ',xyzzyaaat2,ad&
&just=.true.)
call wout('Local potential energy                    : ',xyzzyaaav2,ad&
&just=.true.)
call wout('Non-local potential energy                : ',xyzzyaaaw2,ad&
&just=.true.)
call wout('Electron-electron energy                  : ',xyzzyaaax2,ad&
&just=.true.)
call wout('Ion-ion energy                            : ',eionion,adjus&
&t=.true.)
call wout()
if(xyzzyaabk1)then
call wout('Real blip coefficients ==> GAMMA calculation.')
else
call wout('Complex blip coefficients ==> calculation with K-POINTS.')
endif
end subroutine xyzzyaabo2
subroutine xyzzyaabp2
implicit none
integer xyzzyaaaa4,xyzzyaaab4,io,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzz&
&yaaaf4,xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4,xyzzyaaak4,xyzzyaa&
&al4,xyzzyaaam4(3)
real(dp) xyzzyaaan4(3),xyzzyaaao4
real(dp),allocatable :: xyzzyaaap4(:,:,:)
logical xyzzyaaaq4,xyzzyaaar4
call timer("READ_HDR",.true.)
call open_units(io,xyzzyaaab4)
if(xyzzyaaab4/=0)call errstop_master('READ_FORMATTED','Unable to find &
&free i/o unit')
open(io,file=trim(xyzzyaabj2),status='old',action='read',iostat=xyzzya&
&aab4)
if(xyzzyaaab4/=0)call errstop_master('READ_FORMATTED','Problem opening&
& '//trim(xyzzyaabj2)//'.')
if(am_master)then
read(io,'(a)',err=30,end=20)title
call skip(io,4)
read(io,'(a)',err=30,end=20)code
call skip(io,1)
read(io,'(a)',err=30,end=20)method
call skip(io,1)
read(io,'(a)',err=30,end=20)functional
call skip(io,1)
read(io,'(a)',err=30,end=20)pseudo_type
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaau2
call skip(io,1)
read(io,*,end=20,err=30)xyzzyaabl1
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaas2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaat2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaav2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaaw2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaax2
call skip(io,1)
read(io,*,err=30,end=20)eionion
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaae2
call skip(io,4)
read(io,*,end=20,err=30)nbasis
call skip(io,1)
allocate(atno(nbasis),basis(3,nbasis),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','atno,basis')
do xyzzyaaaj4=1,nbasis
read(io,*,err=30,end=20)atno(xyzzyaaaj4),basis(1:3,xyzzyaaaj4)
enddo
call skip(io,1)
read(io,*,end=20,err=30)pa1
read(io,*,end=20,err=30)pa2
read(io,*,end=20,err=30)pa3
call skip(io,4)
read(io,*,end=20,err=30)xyzzyaaaf2
call skip(io,xyzzyaaaf2+2)
read(io,*,end=20,err=30)xyzzyaaaa1
call skip(io,4)
read(io,*,end=20,err=30)xyzzyaaaf1
allocate(xyzzyaaay1(3,xyzzyaaaf1),xyzzyaaah1(xyzzyaaaf1,2),xyzzyaaaz1(&
&xyzzyaaaf1),xyzzyaaba1(6,xyzzyaaaf1),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','kvec,nband,ksq,kprod')
endif
call mpi_bcast(xyzzyaabl1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting spin_polarized in read_formatted.')
call mpi_bcast(nbasis,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nbasis in read_formatted.')
if(am_slave)then
allocate(atno(nbasis),basis(3,nbasis),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','atno,basis')
endif
call mpi_bcast(atno,nbasis,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting atno in read_formatted.')
call mpi_bcast(basis,3*nbasis,mpi_double_precision,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting basis in read_formatted.')
call mpi_bcast(pa1,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa1 in read_formatted.')
call mpi_bcast(pa2,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa2 in read_formatted.')
call mpi_bcast(pa3,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa3 in read_formatted.')
call mpi_bcast(xyzzyaaaa1,3,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nr in read_formatted.')
call mpi_bcast(xyzzyaaaf1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nkvec in read_formatted.')
xyzzyaaab1=dble(xyzzyaaaa1)
xyzzyaaac1(1)=xyzzyaaab1(1)**2
xyzzyaaac1(2)=xyzzyaaab1(2)**2
xyzzyaaac1(3)=xyzzyaaab1(3)**2
if(xyzzyaabl1)then
xyzzyaaad2=2
else
xyzzyaaad2=1
endif
if(am_master)then
xyzzyaaad1=0
do xyzzyaaae4=1,xyzzyaaaf1
call skip(io,1)
read(io,*,end=20,err=30)xyzzyaaal4,xyzzyaaah1(xyzzyaaae4,1:2),xyzzyaaa&
&y1(1:3,xyzzyaaae4)
xyzzyaaaz1(xyzzyaaae4)=dot_product(xyzzyaaay1(:,xyzzyaaae4),xyzzyaaay1&
&(:,xyzzyaaae4))
xyzzyaaba1(1,xyzzyaaae4)=xyzzyaaay1(1,xyzzyaaae4)*xyzzyaaay1(1,xyzzyaa&
&ae4)
xyzzyaaba1(2,xyzzyaaae4)=xyzzyaaay1(2,xyzzyaaae4)*xyzzyaaay1(2,xyzzyaa&
&ae4)
xyzzyaaba1(3,xyzzyaaae4)=xyzzyaaay1(3,xyzzyaaae4)*xyzzyaaay1(3,xyzzyaa&
&ae4)
xyzzyaaba1(4,xyzzyaaae4)=xyzzyaaay1(1,xyzzyaaae4)*xyzzyaaay1(2,xyzzyaa&
&ae4)
xyzzyaaba1(5,xyzzyaaae4)=xyzzyaaay1(1,xyzzyaaae4)*xyzzyaaay1(3,xyzzyaa&
&ae4)
xyzzyaaba1(6,xyzzyaaae4)=xyzzyaaay1(2,xyzzyaaae4)*xyzzyaaay1(3,xyzzyaa&
&ae4)
if(xyzzyaaae4==1)then
xyzzyaaad4=maxval(xyzzyaaah1(1,1:xyzzyaaad2))
allocate(xyzzyaaax1(xyzzyaaad4,xyzzyaaaf1,xyzzyaaad2),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','eigenvalue')
xyzzyaaax1=0.d0
elseif(maxval(xyzzyaaah1(xyzzyaaae4,1:xyzzyaaad2))>xyzzyaaad4)then
xyzzyaaad4=maxval(xyzzyaaah1(xyzzyaaae4,1:xyzzyaaad2))
allocate(xyzzyaaap4(xyzzyaaad4,xyzzyaaaf1,xyzzyaaad2),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','eigenvalue_temp')
xyzzyaaap4=0.d0
do xyzzyaaaf4=1,xyzzyaaae4-1
do xyzzyaaak4=1,xyzzyaaad2
do xyzzyaaag4=1,xyzzyaaah1(xyzzyaaaf4,xyzzyaaak4)
xyzzyaaap4(xyzzyaaag4,xyzzyaaaf4,xyzzyaaak4)=xyzzyaaax1(xyzzyaaag4,xyz&
&zyaaaf4,xyzzyaaak4)
enddo
enddo
enddo
deallocate(xyzzyaaax1)
allocate(xyzzyaaax1(xyzzyaaad4,xyzzyaaaf1,xyzzyaaad2),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','eigenvalue <2>')
xyzzyaaax1=xyzzyaaap4
deallocate(xyzzyaaap4)
endif
do xyzzyaaak4=1,xyzzyaaad2
xyzzyaaar4=.false.
xyzzyaaah4=0
do xyzzyaaag4=1,xyzzyaaah1(xyzzyaaae4,xyzzyaaak4)
call skip(io,1)
read(io,'(a)',end=20,err=30)char80
read(char80,*,iostat=xyzzyaaab4)xyzzyaaal4,xyzzyaaal4,xyzzyaaao4,xyzzy&
&aaaq4
if(xyzzyaaab4/=0)xyzzyaaaq4=.false.
if(xyzzyaaaq4)then
if(xyzzyaaar4)call errstop('READ_FORMATTED','Localized orbitals should&
& come before extended ones in '//trim(xyzzyaabj2)//'.')
call skip(io,1)
read(io,*,end=20,err=30)xyzzyaaan4(1:3),xyzzyaaar2,xyzzyaaar2,xyzzyaaa&
&r2,xyzzyaaac4
call skip(io,1)
read(io,*,end=20,err=30)xyzzyaaal4,xyzzyaaal4,xyzzyaaal4,xyzzyaaam4(1:&
&3)
xyzzyaaad1(xyzzyaaak4)=xyzzyaaad1(xyzzyaaak4)+1
else
xyzzyaaar4=.true.
endif
call skip(io,1)
if(xyzzyaaae4==1.and.xyzzyaaak4==1.and.xyzzyaaag4==1)then
read(io,fmt='(a)',end=20,err=30)char80
if(scan(char80,',')/=0)then
xyzzyaabk1=.false.
else
xyzzyaabk1=.true.
endif
backspace(io)
endif
if(xyzzyaaaq4)then
if(xyzzyaaak4/=xyzzyaaad2.or.xyzzyaaag4/=xyzzyaaah1(xyzzyaaae4,xyzzyaa&
&ak4))call skip(io,xyzzyaaam4(1)*xyzzyaaam4(2)*xyzzyaaam4(3))
else
xyzzyaaah4=xyzzyaaah4+1
xyzzyaaax1(xyzzyaaah4,xyzzyaaae4,xyzzyaaak4)=xyzzyaaao4
if(xyzzyaaae4/=xyzzyaaaf1.or.xyzzyaaak4/=xyzzyaaad2.or.xyzzyaaag4/=xyz&
&zyaaah1(xyzzyaaae4,xyzzyaaak4))call skip(io,xyzzyaaaa1(1)*xyzzyaaaa1(&
&2)*xyzzyaaaa1(3))
endif
enddo
enddo
enddo
xyzzyaaah1(1,1:xyzzyaaad2)=xyzzyaaah1(1,1:xyzzyaaad2)-xyzzyaaad1(1:xyz&
&zyaaad2)
xyzzyaaag1=maxval(xyzzyaaah1(:,1:xyzzyaaad2))
xyzzyaabn1=any(xyzzyaaah1(:,1:xyzzyaaad2)>0)
xyzzyaaay2=.false.
if(.not.modified_mdet)then
do
read(io,'(a)',iostat=xyzzyaaab4)char80
if(xyzzyaaab4>0)call errstop('READ_FORMATTED','Problem reading '//trim&
&(xyzzyaabj2)//'.')
if(xyzzyaaab4<0)exit
if(trim(adjustl(char80))=='WAVEFUNCTION')then
call read_mdet_wfn(io)
xyzzyaaay2=.true.
exit
endif
enddo
endif
rewind(io)
call skip(io,53+nbasis+xyzzyaaaf2)
allocate(xyzzyaaap4(xyzzyaaag1,xyzzyaaaf1,xyzzyaaad2))
xyzzyaaap4=0.d0
do xyzzyaaae4=1,xyzzyaaaf1
do xyzzyaaak4=1,xyzzyaaad2
do xyzzyaaag4=1,xyzzyaaah1(xyzzyaaae4,xyzzyaaak4)
xyzzyaaap4(xyzzyaaag4,xyzzyaaae4,xyzzyaaak4)=xyzzyaaax1(xyzzyaaag4,xyz&
&zyaaae4,xyzzyaaak4)
enddo
enddo
enddo
deallocate(xyzzyaaax1)
allocate(xyzzyaaax1(xyzzyaaag1,xyzzyaaaf1,xyzzyaaad2),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READBWF','eigenvalue')
xyzzyaaax1=xyzzyaaap4
deallocate(xyzzyaaap4)
endif
call mpi_bcast(xyzzyaaay2,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting mdet_found in read_formatted.')
if(xyzzyaaay2)call broadcast_mdet
call mpi_bcast(xyzzyaaag1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting maxband in read_formatted.')
call mpi_bcast(xyzzyaaad1,2,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_loc_orbs in read_formatted.')
call mpi_bcast(xyzzyaabk1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting gamma_only in read_formatted.')
call mpi_bcast(xyzzyaabn1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ext_orbs_present in read_formatted.&
&')
if(am_slave)then
allocate(xyzzyaaay1(3,xyzzyaaaf1),xyzzyaaaz1(xyzzyaaaf1),xyzzyaaba1(6,&
&xyzzyaaaf1),xyzzyaaah1(xyzzyaaaf1,2),xyzzyaaax1(xyzzyaaag1,xyzzyaaaf1&
&,xyzzyaaad2),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','kvec,ksq,kprod,nband')
endif
call mpi_bcast(xyzzyaaay1,3*xyzzyaaaf1,mpi_double_precision,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting kvec in read_formatted.')
call mpi_bcast(xyzzyaaaz1,xyzzyaaaf1,mpi_double_precision,0,mpi_comm_w&
&orld,ierror)
call checkmpi(ierror,'Broadcasting ksq in read_formatted.')
call mpi_bcast(xyzzyaaba1,6*xyzzyaaaf1,mpi_double_precision,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting kprod in read_formatted.')
call mpi_bcast(xyzzyaaah1,xyzzyaaaf1*2,mpi_integer,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting nband in read_formatted.')
call mpi_bcast(xyzzyaaax1,xyzzyaaag1*xyzzyaaaf1*xyzzyaaad2,mpi_double_&
&precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting eigenvalue in read_formatted.')
xyzzyaabm1=any(xyzzyaaad1(1:xyzzyaaad2)>0)
xyzzyaacj1=sum(xyzzyaaad1(1:xyzzyaaad2))
call timer("READ_HDR",.false.)
call timer("OCCUPANCY",.true.)
avc_norb=0
if(xyzzyaabn1)then
call bwfdet_setup(.true.,.true.)
do xyzzyaaak4=1,xyzzyaaad2
do xyzzyaaae4=1,xyzzyaaaf1
do xyzzyaaah4=1,xyzzyaaah1(xyzzyaaae4,xyzzyaaak4)
if((xyzzyaabl1.and.any(xyzzyaaai1(1:ndet,xyzzyaaah4,xyzzyaaae4,xyzzyaa&
&ak4)>0)).or.(.not.xyzzyaabl1.and.any(xyzzyaaai1(1:ndet,xyzzyaaah4,xyz&
&zyaaae4,1:nspin)>0)))avc_norb=avc_norb+1
enddo
enddo
enddo
endif
call timer("OCCUPANCY",.false.)
call timer("ALLOC_DATA",.true.)
if(xyzzyaabk1)then
if(xyzzyaabm1)then
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaacj1),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','AVCLOC_SP')
else
allocate(xyzzyaaan1(xyzzyaacj1),stat=xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','AVCLOC_DP')
endif
allocate(xyzzyaabi1(3,xyzzyaacj1),xyzzyaabg1(3,xyzzyaacj1),xyzzyaabb1(&
&3,xyzzyaacj1),xyzzyaabc1(xyzzyaacj1),xyzzyaabe1(xyzzyaacj1),xyzzyaabd&
&1(xyzzyaacj1),xyzzyaabf1(xyzzyaacj1),xyzzyaabh1(xyzzyaacj1),stat=xyzz&
&yaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','centres,...')
xyzzyaabi1=0
endif
if(xyzzyaabn1)then
if(single_precision_blips)then
call shallocate_blip(xyzzyaaao1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','AVC_SP')
if(am_smpmaster)xyzzyaaao1=0.
else
call shallocate_blip(xyzzyaaap1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','AVC_DP')
if(am_smpmaster)xyzzyaaap1=0.d0
endif
endif
else
if(xyzzyaabm1)call errstop_master('READ_FORMATTED','Localized orbitals&
& present for non-Gamma orbitals.')
if(single_precision_blips)then
call shallocate_blip(xyzzyaaaq1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','CAVC_SP')
if(am_smpmaster)xyzzyaaaq1=cmplx(0.0,0.0,sp)
else
call shallocate_blip(xyzzyaaar1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','CAVC_DP')
if(am_smpmaster)xyzzyaaar1=czero
endif
endif
call timer("ALLOC_DATA",.false.)
if(am_master)then
call timer("READ_DATA",.true.)
xyzzyaaai4=0
do xyzzyaaae4=1,xyzzyaaaf1
call skip(io,2)
xyzzyaaah4=0
do xyzzyaaak4=1,xyzzyaaad2
do xyzzyaaag4=1,xyzzyaaad1(xyzzyaaak4)
xyzzyaaah4=xyzzyaaah4+1
call skip(io,3)
read(io,*,err=30,end=20)xyzzyaabb1(:,xyzzyaaah4),xyzzyaabc1(xyzzyaaah4&
&),xyzzyaabe1(xyzzyaaah4),xyzzyaaar2,xyzzyaabh1(xyzzyaaah4)
xyzzyaabd1(xyzzyaaah4)=xyzzyaabc1(xyzzyaaah4)**2
xyzzyaabf1(xyzzyaaah4)=xyzzyaabe1(xyzzyaaah4)**2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaabg1(1:3,xyzzyaaah4),xyzzyaabi1(:,xyzzya&
&aah4)
call skip(io,1)
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaaah4)%bcoefs(0:xyzzyaabi1(3,xyzzyaaah4)-1,0:&
&xyzzyaabi1(2,xyzzyaaah4)-1,0:xyzzyaabi1(1,xyzzyaaah4)-1),stat=xyzzyaa&
&aa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','avcloc_sp%bcoefs')
read(io,*,end=20,err=30)xyzzyaaam1(xyzzyaaah4)%bcoefs(:,:,:)
else
allocate(xyzzyaaan1(xyzzyaaah4)%bcoefs(0:xyzzyaabi1(3,xyzzyaaah4)-1,0:&
&xyzzyaabi1(2,xyzzyaaah4)-1,0:xyzzyaabi1(1,xyzzyaaah4)-1),stat=xyzzyaa&
&aa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','avcloc_dp%bcoefs')
read(io,*,end=20,err=30)xyzzyaaan1(xyzzyaaah4)%bcoefs(:,:,:)
endif
enddo
do xyzzyaaag4=1,xyzzyaaah1(xyzzyaaae4,xyzzyaaak4)
call skip(io,3)
if((xyzzyaabl1.and.any(xyzzyaaai1(1:ndet,xyzzyaaag4,xyzzyaaae4,xyzzyaa&
&ak4)>0)).or.(.not.xyzzyaabl1.and.any(xyzzyaaai1(1:ndet,xyzzyaaag4,xyz&
&zyaaae4,1:nspin)>0)))then
xyzzyaaai4=xyzzyaaai4+1
if(single_precision_blips)then
if(xyzzyaabk1)then
read(io,*,end=20,err=30)xyzzyaaao1(xyzzyaaai4,:,:,:)
else
read(io,*,end=20,err=30)xyzzyaaaq1(xyzzyaaai4,:,:,:)
endif
else
if(xyzzyaabk1)then
read(io,*,end=20,err=30)xyzzyaaap1(xyzzyaaai4,:,:,:)
else
read(io,*,end=20,err=30)xyzzyaaar1(xyzzyaaai4,:,:,:)
endif
endif
else
call skip(io,xyzzyaaaa1(1)*xyzzyaaaa1(2)*xyzzyaaaa1(3))
endif
enddo
enddo
enddo
call timer("READ_DATA",.false.)
endif
call shallocate_barrier
do xyzzyaaae4=1,xyzzyaaaf1
do xyzzyaaag4=1,xyzzyaacj1
call mpi_bcast(xyzzyaabb1(1,xyzzyaaag4),3,mpi_double_precision,0,mpi_c&
&omm_world,ierror)
call checkmpi(ierror,'Broadcasting centrepos in read_formatted.')
call mpi_bcast(xyzzyaabc1(xyzzyaaag4),1,mpi_double_precision,0,mpi_com&
&m_world,ierror)
call checkmpi(ierror,'Broadcasting trunc_radius in read_formatted.')
call mpi_bcast(xyzzyaabe1(xyzzyaaag4),1,mpi_double_precision,0,mpi_com&
&m_world,ierror)
call checkmpi(ierror,'Broadcasting trunc_radius_p_dr in read_formatted&
&.')
call mpi_bcast(xyzzyaabh1(xyzzyaaag4),1,mpi_integer,0,mpi_comm_world,i&
&error)
call checkmpi(ierror,'Broadcasting icut in readbwf.')
call mpi_bcast(xyzzyaabg1(1,xyzzyaaag4),3,mpi_integer,0,mpi_comm_world&
&,ierror)
call checkmpi(ierror,'Broadcasting gridstartpos in read_formatted.')
call mpi_bcast(xyzzyaabi1(1,xyzzyaaag4),3,mpi_integer,0,mpi_comm_world&
&,ierror)
call checkmpi(ierror,'Broadcasting nrloc in read_formatted.')
if(.not.am_master)then
xyzzyaabd1(xyzzyaaag4)=xyzzyaabc1(xyzzyaaag4)**2
xyzzyaabf1(xyzzyaaag4)=xyzzyaabe1(xyzzyaaag4)**2
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaaag4)%bcoefs(0:xyzzyaabi1(3,xyzzyaaag4)-1,0:&
&xyzzyaabi1(2,xyzzyaaag4)-1,0:xyzzyaabi1(1,xyzzyaaag4)-1),stat=xyzzyaa&
&aa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','avcloc_sp%bcoef')
else
allocate(xyzzyaaan1(xyzzyaaag4)%bcoefs(0:xyzzyaabi1(3,xyzzyaaag4)-1,0:&
&xyzzyaabi1(2,xyzzyaaag4)-1,0:xyzzyaabi1(1,xyzzyaaag4)-1),stat=xyzzyaa&
&aa4)
call check_alloc(xyzzyaaaa4,'READ_FORMATTED','avcloc_dp%bcoef')
endif
endif
if(single_precision_blips)then
call mpi_bcast_safe(xyzzyaaam1(xyzzyaaag4)%bcoefs,int(xyzzyaabi1(1,xyz&
&zyaaag4),i64)*xyzzyaabi1(2,xyzzyaaag4)*xyzzyaabi1(3,xyzzyaaag4),mpi_r&
&eal,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting avcloc_sp%bcoefs in read_formatted.&
&')
else
call mpi_bcast_safe(xyzzyaaan1(xyzzyaaag4)%bcoefs,xyzzyaabi1(1,xyzzyaa&
&ag4)*xyzzyaabi1(2,xyzzyaaag4)*xyzzyaabi1(3,xyzzyaaag4),mpi_double_pre&
&cision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting avcloc_dp%bcoefs in read_formatted.&
&')
endif
enddo
enddo
if(xyzzyaabn1.and.am_smpmaster)then
if(single_precision_blips)then
if(xyzzyaabk1)then
call mpi_bcast_safe(xyzzyaaao1,int(avc_norb,i64)*xyzzyaaaa1(1)*xyzzyaa&
&aa1(2)*xyzzyaaaa1(3),mpi_real,0,xyzzyaaao2,ierror)
call checkmpi(ierror,'Broadcasting avc_sp in read_formatted.')
else
call mpi_bcast_safe(xyzzyaaaq1,avc_norb*xyzzyaaaa1(1)*xyzzyaaaa1(2)*xy&
&zzyaaaa1(3),mpi_complex,0,xyzzyaaao2,ierror)
call checkmpi(ierror,'Broadcasting cavc_sp in read_formatted.')
endif
else
if(xyzzyaabk1)then
call mpi_bcast_safe(xyzzyaaap1,avc_norb*xyzzyaaaa1(1)*xyzzyaaaa1(2)*xy&
&zzyaaaa1(3),mpi_double_precision,0,xyzzyaaao2,ierror)
call checkmpi(ierror,'Broadcasting avc_sp in read_formatted.')
else
call mpi_bcast_safe(xyzzyaaar1,avc_norb*xyzzyaaaa1(1)*xyzzyaaaa1(2)*xy&
&zzyaaaa1(3),mpi_double_complex,0,xyzzyaaao2,ierror)
call checkmpi(ierror,'Broadcasting cavc_dp in read_formatted.')
endif
endif
endif
call shallocate_barrier
if(am_master)then
close(io)
open_unit(io)=.false.
endif
if(xyzzyaabn1.and..not.xyzzyaabk1.and..not.xyzzyaabo1)then
deallocate(xyzzyaabp1,xyzzyaaai1)
call bwfdet_setup(.true.,.false.)
endif
call xyzzyaabw2
if(avc_norb/=bw_norbc)call errstop_master('READ_FORMATTED','Size of co&
&efficient arrays (avc_norb='//trim(i2s(avc_norb))//') does not match &
&number of unique orbitals (bw_norbc='//trim(i2s(bw_norbc))//'). This &
&indicates a bug.')
if(xyzzyaabn1)deallocate(xyzzyaabp1,xyzzyaaai1)
return
20 call errstop('READ_FORMATTED','Read past end of '//trim(xyzzyaabj2)&
&//' file.')
30 call errstop('READ_FORMATTED','Error reading '//trim(xyzzyaabj2)//'&
& file.')
end subroutine xyzzyaabp2
subroutine xyzzyaabq2
use format_utils, only : log2int
implicit none
integer xyzzyaaap2,xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaa&
&ae5
logical xyzzyaaaf5
character(16) label
if(.not.am_master)return
call timer('WRITE_BINARY',.true.)
call open_units(xyzzyaaap2,xyzzyaaaa5)
if(xyzzyaaaa5/=0)call errstop('WRITE_BINARY','Unable to find free bina&
&ry i/o unit')
inquire(file=trim(xyzzyaabk2),exist=xyzzyaaaf5)
if(xyzzyaaaf5)then
call wout('Deleting corrupt '//trim(xyzzyaabk2)//' before writing out &
&new one.')
open(xyzzyaaap2,file=trim(xyzzyaabk2),form='unformatted',status='old',&
&action='read',iostat=xyzzyaaaa5)
if(xyzzyaaaa5==0)close(xyzzyaaap2,status='delete')
endif
open(xyzzyaaap2,file=trim(xyzzyaabk2),form='unformatted',status='unkno&
&wn',iostat=xyzzyaaaa5)
if(xyzzyaaaa5/=0)call errstop('WRITE_BINARY','Problem opening '//trim(&
&xyzzyaabk2)//'.')
write(xyzzyaaap2)xyzzyaadj1
label='ENDIAN'
write(xyzzyaaap2)label
write(xyzzyaaap2)-1.d0
label='DESCRIBE'
write(xyzzyaaap2)label
write(xyzzyaaap2)title,code,method,functional,pseudo_type,xyzzyaaau2
label='SPIN_POL'
write(xyzzyaaap2)label
write(xyzzyaaap2)log2int(xyzzyaabl1)
label='ENERGIES'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaaas2,xyzzyaaat2,xyzzyaaav2,xyzzyaaaw2,xyzzyaaax&
&2
label='EIONION'
write(xyzzyaaap2)label
write(xyzzyaaap2)eionion
label='PARTICLE'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaaae2,nspin
label='BASIS'
write(xyzzyaaap2)label
write(xyzzyaaap2)nbasis
write(xyzzyaaap2)atno,basis
label='PRIM_VEC'
write(xyzzyaaap2)label
write(xyzzyaaap2)pa1,pa2,pa3
label='NBLIP'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaaaa1
label='LOC_ORB'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaaad1
label='GAMMA'
write(xyzzyaaap2)label
write(xyzzyaaap2)log2int(xyzzyaabk1)
label='EXTORB'
write(xyzzyaaap2)label
write(xyzzyaaap2)log2int(xyzzyaabn1)
label='SPBLIP'
write(xyzzyaaap2)label
write(xyzzyaaap2)log2int(single_precision_blips)
label='BANDDATA'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaaaf1
write(xyzzyaaap2)xyzzyaaay1,xyzzyaaaz1,xyzzyaaba1,xyzzyaaah1
write(xyzzyaaap2)xyzzyaaax1
label='ORBMAP'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaacm1
if(xyzzyaabn1)then
label='NORB_EXT'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaaae1,xyzzyaack1,bw_norbc,xyzzyaacl1
label='BKMAP'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaacx1,xyzzyaacy1
label='SPINMAP'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaacz1
label='ORBANGLE'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaade1,xyzzyaadd1
label='R2Z_IDX'
write(xyzzyaaap2)label
write(xyzzyaaap2)xyzzyaacw1
label='EXT_DATA'
write(xyzzyaaap2)label
write(xyzzyaaap2)avc_norb
do xyzzyaaac5=0,xyzzyaaaa1(1)-1
do xyzzyaaad5=0,xyzzyaaaa1(2)-1
do xyzzyaaae5=0,xyzzyaaaa1(3)-1
if(single_precision_blips)then
if(xyzzyaabk1)then
write(xyzzyaaap2)xyzzyaaao1(:,xyzzyaaae5,xyzzyaaad5,xyzzyaaac5)
else
write(xyzzyaaap2)xyzzyaaaq1(:,xyzzyaaae5,xyzzyaaad5,xyzzyaaac5)
endif
else
if(xyzzyaabk1)then
write(xyzzyaaap2)xyzzyaaap1(:,xyzzyaaae5,xyzzyaaad5,xyzzyaaac5)
else
write(xyzzyaaap2)xyzzyaaar1(:,xyzzyaaae5,xyzzyaaad5,xyzzyaaac5)
endif
endif
enddo
enddo
enddo
endif
if(xyzzyaabm1)then
label='LOC_DATA'
write(xyzzyaaap2)label
do xyzzyaaab5=1,xyzzyaacj1
write(xyzzyaaap2)xyzzyaabi1(:,xyzzyaaab5),xyzzyaabg1(:,xyzzyaaab5),xyz&
&zyaabb1(:,xyzzyaaab5),xyzzyaabc1(xyzzyaaab5),xyzzyaabe1(xyzzyaaab5),x&
&yzzyaabh1(xyzzyaaab5)
enddo
do xyzzyaaab5=1,xyzzyaacj1
if(single_precision_blips)then
write(xyzzyaaap2)xyzzyaaam1(xyzzyaaab5)%bcoefs
else
write(xyzzyaaap2)xyzzyaaan1(xyzzyaaab5)%bcoefs
endif
enddo
endif
close(xyzzyaaap2)
open_unit(xyzzyaaap2)=.false.
call timer('WRITE_BINARY',.false.)
end subroutine xyzzyaabq2
subroutine xyzzyaabr2(success)
use format_utils, only : int2log
implicit none
logical,intent(out) :: success
integer xyzzyaaap2,xyzzyaaaq2,xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaa&
&ad6,xyzzyaaae6,xyzzyaaaf6
real(dp) xyzzyaaag6
character(16) label
character(32) version_string
xyzzyaaaq2=0
if(am_master)then
xyzzyaaci1=0
xyzzyaabj1=.true.
call open_units(xyzzyaaap2,xyzzyaaaq2)
if(xyzzyaaaq2/=0)call errstop('READ_BINARY','Unable to find free binar&
&y i/o unit')
open(xyzzyaaap2,file=trim(xyzzyaabk2),form='unformatted',status='old',&
&iostat=xyzzyaaaq2,action='read',access='sequential')
if(xyzzyaaaq2/=0)call errstop('READ_BINARY','Problem opening '//trim(x&
&yzzyaabk2))
read(xyzzyaaap2,iostat=xyzzyaaaq2)version_string
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaaaq2=0
if(version_string/=xyzzyaadj1)then
select case(trim(version_string))
case('CASINO v2.10','CASINO v2.10-1')
call errwarn('READ_BINARY','This binary file was generated by a develo&
&pment build of CASINO 2.9.  The format changed shortly thereafter and&
& is not supported any more.')
case default
call errwarn('READ_BINARY','Unsupported binary version: '//trim(versio&
&n_string))
end select
xyzzyaaaq2=1
endif
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
do
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)label
if(xyzzyaaaq2/=0)label='exit'
endif
call mpi_bcast(label,16,mpi_character,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting label in read_binary.')
select case(trim(label))
case('ENDIAN')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaag6
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
if(abs(xyzzyaaag6+1.d0)>1.d-7)then
call errwarn('READ_BINARY','Binary file '//trim(xyzzyaabk2)//' was cre&
&ated on a machine with different endianness.')
xyzzyaaaq2=1
endif
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('DESCRIBE')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)title,code,method,functional,pseudo_&
&type,xyzzyaaau2
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('SPIN_POL')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaac6
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaabl1=int2log(xyzzyaaac6)
if(xyzzyaabl1)then
xyzzyaaad2=2
else
xyzzyaaad2=1
endif
endif
case('ENERGIES')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaas2,xyzzyaaat2,xyzzyaaav2,xyz&
&zyaaaw2,xyzzyaaax2
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('EIONION')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)eionion
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('PARTICLE')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaae2,xyzzyaaag2
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaaaq2=0
if(xyzzyaaag2/=nspin)then
call wordwrap('The number of spins specified in the input file is not &
&consistent with the number of spins stored in the binary file '//trim&
&(xyzzyaabk2)//'.')
call wout()
call wout('   nspin       = '//trim(i2s(nspin)))
call wout('   nspin_check = '//trim(i2s(xyzzyaaag2)))
call wout()
call wordwrap('Either remove the binary file '//trim(xyzzyaabk2)//' or&
& check the values of NEU and NED in the input file.')
xyzzyaaaq2=1
endif
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('BASIS')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)nbasis
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
allocate(atno(nbasis),basis(3,nbasis),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READBWF','atno,basis')
read(xyzzyaaap2,iostat=xyzzyaaaq2)atno,basis
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('PRIM_VEC')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)pa1,pa2,pa3
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('NBLIP')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaaa1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
call mpi_bcast(xyzzyaaaa1,3,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nr box in read_binary.')
xyzzyaaab1=dble(xyzzyaaaa1)
xyzzyaaac1(1:3)=xyzzyaaab1(1:3)**2
case('LOC_ORB')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaad1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaabm1=any(xyzzyaaad1(1:xyzzyaaad2)>0)
xyzzyaacj1=sum(xyzzyaaad1(1:xyzzyaaad2))
xyzzyaaci1=xyzzyaaci1+xyzzyaacj1
endif
case('GAMMA')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaac6
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)xyzzyaabk1=int2log(xyzzyaaac6)
call mpi_bcast(xyzzyaabk1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting gamma_only in read_binary.')
case('EXTORB')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaac6
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaabn1=int2log(xyzzyaaac6)
endif
case('SPBLIP')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaac6
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaaaq2=0
xyzzyaaaz2=int2log(xyzzyaaac6)
if(single_precision_blips.neqv.xyzzyaaaz2)then
call wout()
call wordwrap('The value of single_precision_blips is not consistent w&
&ith the value stored in the binary file '//trim(xyzzyaabk2)//'.')
call wout()
call wout('   value in CASINO input  = '//l2s(single_precision_blips))
call wout('   value in '//trim(xyzzyaabk2)//' = '//l2s(xyzzyaaaz2))
call wout()
xyzzyaaaq2=1
call wordwrap('Either remove the binary file '//trim(xyzzyaabk2)//' or&
& check the value of single_precision_blips in the input file.')
endif
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('BANDDATA')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaaf1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
allocate(xyzzyaaay1(3,xyzzyaaaf1),xyzzyaaaz1(xyzzyaaaf1),xyzzyaaba1(6,&
&xyzzyaaaf1),xyzzyaaah1(xyzzyaaaf1,2),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','kvec,...')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaay1,xyzzyaaaz1,xyzzyaaba1,xyz&
&zyaaah1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaaag1=maxval(xyzzyaaah1)
allocate(xyzzyaaax1(xyzzyaaag1,xyzzyaaaf1,xyzzyaaad2),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','eigenvalue')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaax1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('NORB_EXT')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaae1,xyzzyaack1,bw_norbc,xyzzy&
&aacl1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaaci1=xyzzyaaci1+xyzzyaack1
endif
case('ORBMAP')
if(am_master)then
allocate(xyzzyaacm1(nemax,nspin,ndet),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','bw_orbmap')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaacm1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('BKMAP')
if(am_master)then
allocate(xyzzyaacx1(bw_norbc),xyzzyaacy1(bw_norbc),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','corb_band_idx,...')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaacx1,xyzzyaacy1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('SPINMAP')
if(am_master)then
allocate(xyzzyaacz1(bw_norbc),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READBWF','corb_spin_idx')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaacz1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('ORBANGLE')
if(am_master)then
allocate(xyzzyaade1(xyzzyaacl1),xyzzyaadd1(xyzzyaacl1),stat=xyzzyaaaa6&
&)
call check_alloc(xyzzyaaaa6,'READ_BINARY','rorb_angle,...')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaade1,xyzzyaadd1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('R2Z_IDX')
if(am_master)then
allocate(xyzzyaacw1(xyzzyaacl1),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','Ridx2Zidx')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaacw1
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
case('LOC_DATA')
call timer("READ_LOC",.true.)
if(am_master)then
if(.not.xyzzyaabk1)call errstop('READ_BINARY','Localized orbitals pres&
&ent for non-Gamma orbitals.')
allocate(xyzzyaabi1(3,xyzzyaacj1),xyzzyaabg1(3,xyzzyaacj1),xyzzyaabb1(&
&3,xyzzyaacj1),xyzzyaabc1(xyzzyaacj1),xyzzyaabe1(xyzzyaacj1),xyzzyaabd&
&1(xyzzyaacj1),xyzzyaabf1(xyzzyaacj1),xyzzyaabh1(xyzzyaacj1),stat=xyzz&
&yaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','centres,...')
xyzzyaabi1=0
do xyzzyaaab6=1,xyzzyaacj1
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaabi1(:,xyzzyaaab6),xyzzyaabg1(&
&:,xyzzyaaab6),xyzzyaabb1(:,xyzzyaaab6),xyzzyaabc1(xyzzyaaab6),xyzzyaa&
&be1(xyzzyaaab6),xyzzyaabh1(xyzzyaaab6)
if(xyzzyaaaq2/=0)exit
xyzzyaabd1(xyzzyaaab6)=xyzzyaabc1(xyzzyaaab6)**2
xyzzyaabf1(xyzzyaaab6)=xyzzyaabe1(xyzzyaaab6)**2
if(any(xyzzyaabi1(:,xyzzyaaab6)>xyzzyaaaa1-3))xyzzyaabj1=.false.
enddo
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaacj1),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','AVCLOC_SP')
else
allocate(xyzzyaaan1(xyzzyaacj1),stat=xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','AVCLOC_DP')
endif
do xyzzyaaab6=1,xyzzyaacj1
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaaab6)%bcoefs(0:xyzzyaabi1(3,xyzzyaaab6)-1,0:&
&xyzzyaabi1(2,xyzzyaaab6)-1,0:xyzzyaabi1(1,xyzzyaaab6)-1),stat=xyzzyaa&
&aa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','avcloc_sp%bcoefs')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaam1(xyzzyaaab6)%bcoefs
if(xyzzyaaaq2/=0)exit
else
allocate(xyzzyaaan1(xyzzyaaab6)%bcoefs(0:xyzzyaabi1(3,xyzzyaaab6)-1,0:&
&xyzzyaabi1(2,xyzzyaaab6)-1,0:xyzzyaabi1(1,xyzzyaaab6)-1),stat=xyzzyaa&
&aa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','avcloc_dp%bcoefs')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaan1(xyzzyaaab6)%bcoefs
if(xyzzyaaaq2/=0)exit
endif
enddo
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
call timer("READ_LOC",.false.)
case('EXT_DATA')
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)avc_norb
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
call mpi_bcast(avc_norb,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting avc_norb in read_binary.')
call timer("READ_EXT",.true.)
if(single_precision_blips)then
if(xyzzyaabk1)then
call shallocate_blip(xyzzyaaao1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','AVC_SP')
else
call shallocate_blip(xyzzyaaaq1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','CAVC_SP')
endif
else
if(xyzzyaabk1)then
call shallocate_blip(xyzzyaaap1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','AVC_DP')
else
call shallocate_blip(xyzzyaaar1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa6)
call check_alloc(xyzzyaaaa6,'READ_BINARY','CAVC_DP')
endif
endif
xyzzyaabg2=.true.
call timer("READ_DATA",.true.)
if(am_master)then
n1_loop: do xyzzyaaad6=0,xyzzyaaaa1(1)-1
do xyzzyaaae6=0,xyzzyaaaa1(2)-1
do xyzzyaaaf6=0,xyzzyaaaa1(3)-1
if(single_precision_blips)then
if(xyzzyaabk1)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaao1(:,xyzzyaaaf6,xyzzyaaae6,x&
&yzzyaaad6)
else
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaaq1(:,xyzzyaaaf6,xyzzyaaae6,x&
&yzzyaaad6)
endif
else
if(xyzzyaabk1)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaap1(:,xyzzyaaaf6,xyzzyaaae6,x&
&yzzyaaad6)
else
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaar1(:,xyzzyaaaf6,xyzzyaaae6,x&
&yzzyaaad6)
endif
endif
if(xyzzyaaaq2/=0)exit n1_loop
enddo
enddo
enddo n1_loop
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
call timer("READ_DATA",.false.)
call timer("READ_EXT",.false.)
case('exit')
exit
case default
call errstop('READ_BINARY','Binary file label not recognized.')
end select
enddo
if(am_master)then
if(.not.xyzzyaabn1.and.any(nuc_nele(1:xyzzyaaad2)/=xyzzyaaad1(1:xyzzya&
&aad2)))then
call wordwrap('The number of electrons specified in the input file is &
&not consistent with the number of orbitals stored in the binary file &
&'//trim(xyzzyaabk2)//'.')
call wout()
call wout('   neu = '//trim(i2s(nuc_nele(1))))
call wout('   ned = '//trim(i2s(nuc_nele(2))))
call wout()
call wout('   number of up   spin orbitals = '//trim(i2s(xyzzyaaad1(1)&
&)))
call wout('   number of down spin orbitals = '//trim(i2s(xyzzyaaad1(xy&
&zzyaaad2))))
call wout()
call wordwrap('Either remove the binary file '//trim(xyzzyaabk2)//' or&
& check the values of neu and ned in the input file.')
call errstop('READ_BINARYBWF','Problem with neu/ned or '//trim(xyzzyaa&
&bk2)//' file.')
endif
endif
if(am_master)then
close(xyzzyaaap2)
open_unit(xyzzyaaap2)=.false.
endif
call shallocate_barrier
success=.true.
end subroutine xyzzyaabr2
logical function xyzzyaabs2(xyzzyaaaq2,success)
implicit none
integer,intent(in) :: xyzzyaaaq2
logical,intent(inout) :: success
if(xyzzyaaaq2/=0)then
call wordwrap('Problem reading binary file.')
call wout()
call wout('--- Data read in failed operation ---')
call wout()
call xyzzyaabo2
call wout()
call wout('--- End of data read in failed operation ---')
call wout()
success=.false.
else
success=.true.
endif
call mpi_bcast(success,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting success in successful_read.')
xyzzyaabs2=success
end function xyzzyaabs2
subroutine xyzzyaabt2
implicit none
if(allocated(atno))deallocate(atno)
if(allocated(basis))deallocate(basis)
if(allocated(xyzzyaaay1))deallocate(xyzzyaaay1)
if(allocated(xyzzyaaah1))deallocate(xyzzyaaah1)
if(allocated(xyzzyaaaz1))deallocate(xyzzyaaaz1)
if(allocated(xyzzyaaba1))deallocate(xyzzyaaba1)
if(allocated(xyzzyaaax1))deallocate(xyzzyaaax1)
if(single_precision_blips)then
if(allocated(xyzzyaaam1))deallocate(xyzzyaaam1)
else
if(allocated(xyzzyaaan1))deallocate(xyzzyaaan1)
endif
if(allocated(xyzzyaabb1))deallocate(xyzzyaabb1)
if(allocated(xyzzyaabc1))deallocate(xyzzyaabc1)
if(allocated(xyzzyaabe1))deallocate(xyzzyaabe1)
if(allocated(xyzzyaabd1))deallocate(xyzzyaabd1)
if(allocated(xyzzyaabf1))deallocate(xyzzyaabf1)
if(allocated(xyzzyaabh1))deallocate(xyzzyaabh1)
if(allocated(xyzzyaabg1))deallocate(xyzzyaabg1)
if(allocated(xyzzyaabi1))deallocate(xyzzyaabi1)
if(allocated(xyzzyaacm1))deallocate(xyzzyaacm1)
if(allocated(xyzzyaacx1))deallocate(xyzzyaacx1)
if(allocated(xyzzyaacy1))deallocate(xyzzyaacy1)
if(allocated(xyzzyaacz1))deallocate(xyzzyaacz1)
if(allocated(xyzzyaade1))deallocate(xyzzyaade1)
if(allocated(xyzzyaadd1))deallocate(xyzzyaadd1)
if(allocated(xyzzyaacw1))deallocate(xyzzyaacw1)
if(xyzzyaabg2)then
if(single_precision_blips)then
if(xyzzyaabk1)then
call deshallocate(xyzzyaaao1)
else
call deshallocate(xyzzyaaaq1)
endif
else
if(xyzzyaabk1)then
call deshallocate(xyzzyaaap1)
else
call deshallocate(xyzzyaaar1)
endif
endif
xyzzyaabg2=.false.
endif
if(allocated(xyzzyaaai1))deallocate(xyzzyaaai1)
if(allocated(xyzzyaabp1))deallocate(xyzzyaabp1)
if(allocated(xyzzyaabq1))deallocate(xyzzyaabq1)
if(allocated(xyzzyaaaw1))deallocate(xyzzyaaaw1)
if(allocated(xyzzyaabz1))deallocate(xyzzyaabz1)
if(allocated(xyzzyaaau1))deallocate(xyzzyaaau1)
if(allocated(xyzzyaaav1))deallocate(xyzzyaaav1)
end subroutine xyzzyaabt2
subroutine xyzzyaabu2(success)
implicit none
logical,intent(out) :: success
integer xyzzyaaap2,xyzzyaaaq2,xyzzyaaaa9,xyzzyaaab9,xyzzyaaac9,xyzzyaa&
&ad9,xyzzyaaae9,xyzzyaaaf9,xyzzyaaag9,xyzzyaaah9,xyzzyaaai9,num_ext_ma&
&x,xyzzyaaaj9(2),xyzzyaaak9,xyzzyaaal9,idum1,idum2,idum3,xyzzyaaam9
integer,allocatable :: xyzzyaaan9(:,:),xyzzyaaao9(:,:,:),xyzzyaaap9(:,&
&:)
logical xyzzyaaaq9
logical,allocatable :: xyzzyaaar9(:,:,:),xyzzyaaas9(:,:)
xyzzyaaaq2=0
if(am_master)then
call open_units(xyzzyaaap2,xyzzyaaaq2)
if(xyzzyaaaq2/=0)call errstop('READ_OLD_BINARY','Unable to find free b&
&inary i/o unit')
open(xyzzyaaap2,file=trim(xyzzyaabl2),form='unformatted',status='old',&
&iostat=xyzzyaaaq2,action='read',access='sequential')
if(xyzzyaaaq2/=0)call errstop('READ_OLD_BINARY','Problem opening '//tr&
&im(xyzzyaabl2)//'.')
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)title,code,method,functional,pseudo_&
&type,xyzzyaaau2,idum1,xyzzyaaas2,xyzzyaaat2,xyzzyaaav2,xyzzyaaaw2,xyz&
&zyaaax2,eionion,xyzzyaaae2,nbasis,xyzzyaaaf2,xyzzyaaaf1,xyzzyaaaa1,xy&
&zzyaaag1,idum2,idum3,xyzzyaaad1,pa1,pa2,pa3,xyzzyaaag2,num_ext_max
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaabl1=idum1/=0
xyzzyaabk1=idum2/=0
xyzzyaabn1=idum3/=0
endif
if(am_master)then
xyzzyaaaq2=0
if(xyzzyaaag2/=nspin)then
call wordwrap('The number of spins specified in the input file is not &
&consistent with the number of spins stored in the binary file '//trim&
&(xyzzyaabl2)//'.')
call wout()
call wout('   nspin       = '//trim(i2s(nspin)))
call wout('   nspin_check = '//trim(i2s(xyzzyaaag2)))
call wout()
call wordwrap('Either remove the binary file '//trim(xyzzyaabl2)//' or&
& check the values of NEU and NED in the input file.')
xyzzyaaaq2=1
endif
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
if(xyzzyaabl1)then
xyzzyaaad2=2
else
xyzzyaaad2=1
endif
xyzzyaabm1=any(xyzzyaaad1(1:xyzzyaaad2)>0)
xyzzyaacj1=sum(xyzzyaaad1(1:xyzzyaaad2))
endif
call mpi_bcast(xyzzyaaaa1,3,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nr box in readbwf.')
call mpi_bcast(xyzzyaabk1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting gamma_only in readbwf.')
call mpi_bcast(xyzzyaabn1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ext_orbs_present in readbwf.')
xyzzyaaab1=dble(xyzzyaaaa1)
xyzzyaaac1(1:3)=xyzzyaaab1(1:3)**2
if(am_master)then
allocate(atno(nbasis),basis(3,nbasis),xyzzyaaay1(3,xyzzyaaaf1),xyzzyaa&
&ah1(xyzzyaaaf1,2),xyzzyaaaz1(xyzzyaaaf1),xyzzyaaba1(6,xyzzyaaaf1),xyz&
&zyaaax1(xyzzyaaag1,xyzzyaaaf1,xyzzyaaad2),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','atno,basis,...')
allocate(xyzzyaaar9(xyzzyaacj1+xyzzyaaag1,xyzzyaaaf1,xyzzyaaad2),stat=&
&xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','templ3array')
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaay1,xyzzyaaaz1,xyzzyaaba1,atn&
&o,basis,xyzzyaaah1,xyzzyaaax1,xyzzyaaar9,xyzzyaaaj9(1:2)
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)deallocate(xyzzyaaar9)
if(am_master)read(xyzzyaaap2,iostat=xyzzyaaaq2)idum1
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)xyzzyaaaz2=idum1/=0
if(am_master)then
if(single_precision_blips.neqv.xyzzyaaaz2)then
call wout()
call wordwrap('The value of single_precision_blips is not consistent w&
&ith the value stored in the binary file '//trim(xyzzyaabl2)//'.')
call wout()
call wout('   value in CASINO input  = '//l2s(single_precision_blips))
call wout('   value in '//trim(xyzzyaabl2)//'  = '//l2s(xyzzyaaaz2))
call wout()
call wordwrap('Either remove the binary file '//trim(xyzzyaabl2)//' or&
& check the value of single_precision_blips in the input file.')
success=.false.
return
endif
endif
if(am_master)then
xyzzyaabj1=.true.
if(xyzzyaabm1)then
if(.not.xyzzyaabk1)call errstop('READ_OLD_BINARY','Localized orbitals &
&present for non-Gamma orbitals.')
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaacj1),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','AVCLOC_SP')
else
allocate(xyzzyaaan1(xyzzyaacj1),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','AVCLOC_DP')
endif
allocate(xyzzyaabi1(3,xyzzyaacj1),xyzzyaabg1(3,xyzzyaacj1),xyzzyaabb1(&
&3,xyzzyaacj1),xyzzyaabc1(xyzzyaacj1),xyzzyaabe1(xyzzyaacj1),xyzzyaabh&
&1(xyzzyaacj1),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','centres,...')
xyzzyaabi1=0
ikloop: do xyzzyaaaf9=1,xyzzyaaaf1
xyzzyaaah9=0
do xyzzyaaai9=1,xyzzyaaad2
do xyzzyaaag9=1,xyzzyaaad1(xyzzyaaai9)
xyzzyaaah9=xyzzyaaah9+1
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaabi1(:,xyzzyaaah9),xyzzyaabg1(&
&:,xyzzyaaah9),xyzzyaabb1(:,xyzzyaaah9),xyzzyaabc1(xyzzyaaah9),xyzzyaa&
&be1(xyzzyaaah9),xyzzyaaar2,xyzzyaabh1(xyzzyaaah9)
if(xyzzyaaaq2/=0)exit ikloop
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaaah9)%bcoefs(0:xyzzyaabi1(3,xyzzyaaah9)-1,0:&
&xyzzyaabi1(2,xyzzyaaah9)-1,0:xyzzyaabi1(1,xyzzyaaah9)-1),stat=xyzzyaa&
&aa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','avcloc_sp%bcoefs')
read(xyzzyaaap2,iostat=xyzzyaaaq2)(((xyzzyaaam1(xyzzyaaah9)%bcoefs(xyz&
&zyaaad9,xyzzyaaac9,xyzzyaaab9),xyzzyaaab9=0,xyzzyaabi1(1,xyzzyaaah9)-&
&1),xyzzyaaac9=0,xyzzyaabi1(2,xyzzyaaah9)-1),xyzzyaaad9=0,xyzzyaabi1(3&
&,xyzzyaaah9)-1)
if(xyzzyaaaq2/=0)exit ikloop
else
allocate(xyzzyaaan1(xyzzyaaah9)%bcoefs(0:xyzzyaabi1(3,xyzzyaaah9)-1,0:&
&xyzzyaabi1(2,xyzzyaaah9)-1,0:xyzzyaabi1(1,xyzzyaaah9)-1),stat=xyzzyaa&
&aa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','avcloc_dp%bcoefs')
read(xyzzyaaap2,iostat=xyzzyaaaq2)(((xyzzyaaan1(xyzzyaaah9)%bcoefs(xyz&
&zyaaad9,xyzzyaaac9,xyzzyaaab9),xyzzyaaab9=0,xyzzyaabi1(1,xyzzyaaah9)-&
&1),xyzzyaaac9=0,xyzzyaabi1(2,xyzzyaaah9)-1),xyzzyaaad9=0,xyzzyaabi1(3&
&,xyzzyaaah9)-1)
if(xyzzyaaaq2/=0)exit ikloop
endif
if(any(xyzzyaabi1(:,xyzzyaaah9)>xyzzyaaaa1-3))xyzzyaabj1=.false.
enddo
enddo
enddo ikloop
endif
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(xyzzyaabn1)then
if(am_master)then
allocate(xyzzyaaan9(xyzzyaaag1*xyzzyaaaf1,nspin),xyzzyaaao9(ndet,xyzzy&
&aaag1*xyzzyaaaf1,nspin),xyzzyaaas9(xyzzyaaag1*xyzzyaaaf1,nspin),xyzzy&
&aaap9(xyzzyaaag1*xyzzyaaaf1,nspin),stat=xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','tempi2array, tempi3arra&
&y, occupied')
if(trim(code)=='PWSCF')then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaan9,xyzzyaaan9,xyzzyaaao9,xyz&
&zyaaap9
else
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaan9,xyzzyaaan9,xyzzyaaao9,xyz&
&zyaaap9
endif
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
if(am_master)then
xyzzyaaas9=xyzzyaaap9/=0
deallocate(xyzzyaaap9)
xyzzyaaaq9=.not.any(xyzzyaaas9(:,:))
deallocate(xyzzyaaas9,xyzzyaaan9,xyzzyaaao9)
call bwfdet_setup(.true.,.true.)
endif
call mpi_bcast(xyzzyaaaq9,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting uncomputed_occupancy in read_old_bi&
&nary.')
endif
if(am_master)then
avc_norb=0
xyzzyaaae1=0
if(xyzzyaabl1)then
do xyzzyaaai9=1,xyzzyaaad2
do xyzzyaaaf9=1,xyzzyaaaf1
do xyzzyaaah9=1,xyzzyaaah1(xyzzyaaaf9,xyzzyaaai9)
if(any(xyzzyaaai1(1:ndet,xyzzyaaah9,xyzzyaaaf9,xyzzyaaai9)>0))then
avc_norb=avc_norb+1
xyzzyaaae1(xyzzyaaai9)=xyzzyaaae1(xyzzyaaai9)+1
endif
enddo
enddo
enddo
else
do xyzzyaaaf9=1,xyzzyaaaf1
do xyzzyaaah9=1,xyzzyaaag1
if(any(xyzzyaaai1(1:ndet,xyzzyaaah9,xyzzyaaaf9,1:nspin)>0))avc_norb=av&
&c_norb+1
do xyzzyaaam9=1,nspin
if(any(xyzzyaaai1(1:ndet,xyzzyaaah9,xyzzyaaaf9,xyzzyaaam9)>0))xyzzyaaa&
&e1(xyzzyaaam9)=xyzzyaaae1(xyzzyaaam9)+1
enddo
enddo
enddo
if(nspin==1)xyzzyaaae1(2)=xyzzyaaae1(1)
endif
endif
call mpi_bcast(avc_norb,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting avc_norb in read_old_binary.')
if(xyzzyaabn1)then
if(single_precision_blips)then
if(xyzzyaabk1)then
call shallocate_blip(xyzzyaaao1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','AVC_SP')
else
call shallocate_blip(xyzzyaaaq1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','CAVC_SP')
endif
else
if(xyzzyaabk1)then
call shallocate_blip(xyzzyaaap1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','AVC_DP')
else
call shallocate_blip(xyzzyaaar1,avc_norb,xyzzyaaaa1(3),xyzzyaaaa1(2),x&
&yzzyaaaa1(1),xyzzyaaaa9)
call check_alloc(xyzzyaaaa9,'READ_OLD_BINARY','CAVC_DP')
endif
endif
if(xyzzyaaaq9)then
if(am_master)then
xyzzyaaae9=0
ikloop2: do xyzzyaaaf9=1,xyzzyaaaf1
do xyzzyaaai9=1,xyzzyaaad2
do xyzzyaaah9=1,xyzzyaaah1(xyzzyaaaf9,xyzzyaaai9)
if((xyzzyaabl1.and.any(xyzzyaaai1(1:ndet,xyzzyaaah9-xyzzyaaad1(xyzzyaa&
&ai9),xyzzyaaaf9,xyzzyaaai9)>0)).or.(.not.xyzzyaabl1.and.any(xyzzyaaai&
&1(1:ndet,xyzzyaaah9-xyzzyaaad1(xyzzyaaai9),xyzzyaaaf9,1:nspin)>0)))th&
&en
xyzzyaaae9=xyzzyaaae9+1
if(single_precision_blips)then
if(xyzzyaabk1)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)(((xyzzyaaao1(xyzzyaaae9,xyzzyaaad9,&
&xyzzyaaac9,xyzzyaaab9),xyzzyaaab9=0,xyzzyaaaa1(1)-1),xyzzyaaac9=0,xyz&
&zyaaaa1(2)-1),xyzzyaaad9=0,xyzzyaaaa1(3)-1)
else
read(xyzzyaaap2,iostat=xyzzyaaaq2)(((xyzzyaaaq1(xyzzyaaae9,xyzzyaaad9,&
&xyzzyaaac9,xyzzyaaab9),xyzzyaaab9=0,xyzzyaaaa1(1)-1),xyzzyaaac9=0,xyz&
&zyaaaa1(2)-1),xyzzyaaad9=0,xyzzyaaaa1(3)-1)
endif
else
if(xyzzyaabk1)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)(((xyzzyaaap1(xyzzyaaae9,xyzzyaaad9,&
&xyzzyaaac9,xyzzyaaab9),xyzzyaaab9=0,xyzzyaaaa1(1)-1),xyzzyaaac9=0,xyz&
&zyaaaa1(2)-1),xyzzyaaad9=0,xyzzyaaaa1(3)-1)
else
read(xyzzyaaap2,iostat=xyzzyaaaq2)(((xyzzyaaar1(xyzzyaaae9,xyzzyaaad9,&
&xyzzyaaac9,xyzzyaaab9),xyzzyaaab9=0,xyzzyaaaa1(1)-1),xyzzyaaac9=0,xyz&
&zyaaaa1(2)-1),xyzzyaaad9=0,xyzzyaaaa1(3)-1)
endif
endif
else
read(xyzzyaaap2,iostat=xyzzyaaaq2)
endif
if(xyzzyaaaq2/=0)exit ikloop2
enddo
enddo
enddo ikloop2
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
else
if(am_master)then
spin_loop: do xyzzyaaai9=1,xyzzyaaad2
xyzzyaaak9=1
if(xyzzyaaai9>1)xyzzyaaak9=sum(xyzzyaaae1(1:xyzzyaaai9-1))+1
xyzzyaaal9=xyzzyaaak9+xyzzyaaae1(xyzzyaaai9)-1
if(xyzzyaaal9<xyzzyaaak9)cycle
do xyzzyaaad9=0,xyzzyaaaa1(3)-1
do xyzzyaaac9=0,xyzzyaaaa1(2)-1
do xyzzyaaab9=0,xyzzyaaaa1(1)-1
if(single_precision_blips)then
if(xyzzyaabk1)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaao1(xyzzyaaak9:xyzzyaaal9,xyz&
&zyaaad9,xyzzyaaac9,xyzzyaaab9)
else
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaaq1(xyzzyaaak9:xyzzyaaal9,xyz&
&zyaaad9,xyzzyaaac9,xyzzyaaab9)
endif
else
if(xyzzyaabk1)then
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaap1(xyzzyaaak9:xyzzyaaal9,xyz&
&zyaaad9,xyzzyaaac9,xyzzyaaab9)
else
read(xyzzyaaap2,iostat=xyzzyaaaq2)xyzzyaaar1(xyzzyaaak9:xyzzyaaal9,xyz&
&zyaaad9,xyzzyaaac9,xyzzyaaab9)
endif
endif
if(xyzzyaaaq2/=0)exit spin_loop
enddo
enddo
enddo
enddo spin_loop
endif
if(.not.xyzzyaabs2(xyzzyaaaq2,success))return
endif
endif
if(am_master)then
if(.not.xyzzyaabn1.and.any(nuc_nele(1:xyzzyaaad2)/=xyzzyaaad1(1:xyzzya&
&aad2)))then
call wordwrap('The number of electrons specified in the input file is &
&not consistent with the number of orbitals stored in the binary file &
&'//trim(xyzzyaabl2)//'.')
call wout()
call wout('   NEU = '//trim(i2s(nuc_nele(1))))
call wout('   NED = '//trim(i2s(nuc_nele(2))))
call wout()
call wout('   number of up   spin orbitals = '//trim(i2s(xyzzyaaad1(1)&
&)))
call wout('   number of down spin orbitals = '//trim(i2s(xyzzyaaad1(xy&
&zzyaaad2))))
call wout()
call wordwrap('Either remove the binary file '//trim(xyzzyaabl2)//' or&
& check the values of NEU and NED in the input file.')
success=.false.
return
endif
endif
if(am_master)then
if(xyzzyaabn1.and..not.xyzzyaabk1.and..not.xyzzyaabo1)then
deallocate(xyzzyaabp1,xyzzyaaai1)
call bwfdet_setup(.true.,.false.)
endif
call xyzzyaabw2
if(avc_norb/=bw_norbc)call errstop('READ_OLD_BINARY','Size of coeffici&
&ent arrays (avc_norb='//trim(i2s(avc_norb))//') does not match number&
& of unique orbitals (bw_norbc='//trim(i2s(bw_norbc))//'). This indica&
&tes a bug.')
if(xyzzyaabn1)deallocate(xyzzyaabp1,xyzzyaaai1)
endif
call shallocate_barrier
success=.true.
return
end subroutine xyzzyaabu2
subroutine xyzzyaabv2
implicit none
integer xyzzyaaaa10,xyzzyaaab10
if(am_slave)xyzzyaaci1=0
call mpi_bcast(eionion,1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting eionion in broadcast_binary.')
call mpi_bcast(nbasis,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nbasis in broadcast_binary.')
call mpi_bcast(xyzzyaaaf1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nkvec in broadcast_binary.')
call mpi_bcast(xyzzyaaag1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting maxband in broadcast_binary.')
call mpi_bcast(xyzzyaabl1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting spin_polarized in broadcast_binary.&
&')
if(.not.am_master)then
if(xyzzyaabl1)then
xyzzyaaad2=2
else
xyzzyaaad2=1
endif
endif
if(am_slave)then
allocate(xyzzyaaay1(3,xyzzyaaaf1),atno(nbasis),basis(3,nbasis),xyzzyaa&
&ah1(xyzzyaaaf1,2),xyzzyaaax1(xyzzyaaag1,xyzzyaaaf1,xyzzyaaad2),xyzzya&
&aaz1(xyzzyaaaf1),xyzzyaaba1(6,xyzzyaaaf1),stat=xyzzyaaaa10)
call check_alloc(xyzzyaaaa10,'BROADCAST_BINARY','kvec,atno,basis,nband&
&,...')
endif
call mpi_bcast(xyzzyaabn1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ext_orbs_present in broadcast_binar&
&y.')
call mpi_bcast(pa1,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa1 in broadcast_binary.')
call mpi_bcast(pa2,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa2 in broadcast_binary.')
call mpi_bcast(pa3,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa3 in broadcast_binary.')
call mpi_bcast(xyzzyaaad1,xyzzyaaad2,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting num_loc_orbs in broadcast_binary.')
call mpi_bcast(atno,nbasis,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting atno in broadcast_binary.')
call mpi_bcast(basis,3*nbasis,mpi_double_precision,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting basis in broadcast_binary.')
call mpi_bcast(xyzzyaaay1,3*xyzzyaaaf1,mpi_double_precision,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting kvec in broadcast_binary.')
call mpi_bcast(xyzzyaaaz1,xyzzyaaaf1,mpi_double_precision,0,mpi_comm_w&
&orld,ierror)
call checkmpi(ierror,'Broadcasting ksq in broadcast_binary.')
call mpi_bcast(xyzzyaaba1,6*xyzzyaaaf1,mpi_double_precision,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting kprod in broadcast_binary.')
call mpi_bcast(xyzzyaaah1,xyzzyaaaf1*2,mpi_integer,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting nband in broadcast_binary.')
call mpi_bcast(xyzzyaaax1,xyzzyaaag1*xyzzyaaaf1*xyzzyaaad2,mpi_double_&
&precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting eigenvalue in broadcast_binary.')
if(am_slave)then
xyzzyaabm1=any(xyzzyaaad1(1:xyzzyaaad2)>0)
xyzzyaacj1=sum(xyzzyaaad1(1:xyzzyaaad2))
xyzzyaaci1=xyzzyaaci1+xyzzyaacj1
if(xyzzyaabk1)then
if(xyzzyaabm1)then
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaacj1),stat=xyzzyaaaa10)
call check_alloc(xyzzyaaaa10,'BROADCAST_BINARY','AVCLOC_SP')
else
allocate(xyzzyaaan1(xyzzyaacj1),stat=xyzzyaaaa10)
call check_alloc(xyzzyaaaa10,'BROADCAST_BINARY','AVCLOC_DP')
endif
allocate(xyzzyaabi1(3,xyzzyaacj1),xyzzyaabg1(3,xyzzyaacj1),xyzzyaabb1(&
&3,xyzzyaacj1),xyzzyaabc1(xyzzyaacj1),xyzzyaabe1(xyzzyaacj1),xyzzyaabd&
&1(xyzzyaacj1),xyzzyaabf1(xyzzyaacj1),xyzzyaabh1(xyzzyaacj1),stat=xyzz&
&yaaaa10)
call check_alloc(xyzzyaaaa10,'BROADCAST_BINARY','splines,knots,nr,cent&
&res,...')
xyzzyaabi1=0
endif
else
if(xyzzyaabm1)call errstop('BROADCAST_BINARY','Localized orbitals pres&
&ent for non-Gamma orbitals.')
endif
endif
call mpi_bcast(xyzzyaabj1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nrloc_bounded in broadcast_binary.'&
&)
do xyzzyaaab10=1,xyzzyaacj1
call mpi_bcast(xyzzyaabb1(1,xyzzyaaab10),3,mpi_double_precision,0,mpi_&
&comm_world,ierror)
call checkmpi(ierror,'Broadcasting centrepos in broadcast_binary.')
call mpi_bcast(xyzzyaabc1(xyzzyaaab10),1,mpi_double_precision,0,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'Broadcasting trunc_radius in broadcast_binary.')
call mpi_bcast(xyzzyaabe1(xyzzyaaab10),1,mpi_double_precision,0,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'Broadcasting trunc_radius_p_dr in broadcast_bina&
&ry.')
call mpi_bcast(xyzzyaabh1(xyzzyaaab10),1,mpi_integer,0,mpi_comm_world,&
&ierror)
call checkmpi(ierror,'Broadcasting icut in readbwf.')
call mpi_bcast(xyzzyaabg1(1,xyzzyaaab10),3,mpi_integer,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'Broadcasting gridstartpos in broadcast_binary.')
call mpi_bcast(xyzzyaabi1(1,xyzzyaaab10),3,mpi_integer,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'Broadcasting nrloc in broadcast_binary.')
if(.not.am_master)then
xyzzyaabd1(xyzzyaaab10)=xyzzyaabc1(xyzzyaaab10)**2
xyzzyaabf1(xyzzyaaab10)=xyzzyaabe1(xyzzyaaab10)**2
if(single_precision_blips)then
allocate(xyzzyaaam1(xyzzyaaab10)%bcoefs(0:xyzzyaabi1(3,xyzzyaaab10)-1,&
&0:xyzzyaabi1(2,xyzzyaaab10)-1,0:xyzzyaabi1(1,xyzzyaaab10)-1),stat=xyz&
&zyaaaa10)
call check_alloc(xyzzyaaaa10,'BROADCAST_BINARY','avcloc_sp%bcoef')
else
allocate(xyzzyaaan1(xyzzyaaab10)%bcoefs(0:xyzzyaabi1(3,xyzzyaaab10)-1,&
&0:xyzzyaabi1(2,xyzzyaaab10)-1,0:xyzzyaabi1(1,xyzzyaaab10)-1),stat=xyz&
&zyaaaa10)
call check_alloc(xyzzyaaaa10,'BROADCAST_BINARY','avcloc_dp%bcoef')
endif
endif
if(single_precision_blips)then
call mpi_bcast_safe(xyzzyaaam1(xyzzyaaab10)%bcoefs,int(xyzzyaabi1(1,xy&
&zzyaaab10),i64)*xyzzyaabi1(2,xyzzyaaab10)*xyzzyaabi1(3,xyzzyaaab10),m&
&pi_real,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting avcloc_sp%bcoefs in broadcast_binar&
&y.')
else
call mpi_bcast_safe(xyzzyaaan1(xyzzyaaab10)%bcoefs,xyzzyaabi1(1,xyzzya&
&aab10)*xyzzyaabi1(2,xyzzyaaab10)*xyzzyaabi1(3,xyzzyaaab10),mpi_double&
&_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting avcloc_dp%bcoefs in broadcast_binar&
&y.')
endif
enddo
if(.not.am_master)then
allocate(xyzzyaacm1(nemax,nspin,ndet),stat=xyzzyaaaa10)
call check_alloc(xyzzyaaaa10,'BROADCAST_BINARY',' bw_orbmap')
endif
call mpi_bcast(xyzzyaacm1,nemax*nspin*ndet,mpi_integer,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'Broadcasting bw_orbmap in broadcast_binary.')
if(xyzzyaabn1)then
call mpi_bcast(xyzzyaaae1,2,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_ext_orbs in broadcast_binary.')
call mpi_bcast(xyzzyaack1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting bw_norb_ext in broadcast_binary.')
call mpi_bcast(bw_norbc,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting bw_norbc in broadcast_binary.')
call mpi_bcast(xyzzyaacl1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting bw_norbr in broadcast_binary.')
if(.not.am_master)then
xyzzyaaci1=xyzzyaaci1+xyzzyaack1
allocate(xyzzyaacx1(bw_norbc),xyzzyaacy1(bw_norbc),xyzzyaacz1(bw_norbc&
&),xyzzyaade1(xyzzyaacl1),xyzzyaadd1(xyzzyaacl1),xyzzyaacw1(xyzzyaacl1&
&),stat=xyzzyaaaa10)
call check_alloc(xyzzyaaaa10,'BROADCAST_BINARY',' corb_*_idx,...')
endif
call mpi_bcast(xyzzyaacx1,bw_norbc,mpi_integer,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'Broadcasting corb_band_idx in broadcast_binary.'&
&)
call mpi_bcast(xyzzyaacy1,bw_norbc,mpi_integer,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'Broadcasting corb_kvec_idx in broadcast_binary.'&
&)
call mpi_bcast(xyzzyaacz1,bw_norbc,mpi_integer,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'Broadcasting corb_spin_idx in broadcast_binary.'&
&)
call mpi_bcast(xyzzyaade1,xyzzyaacl1,mpi_double_complex,0,mpi_comm_wor&
&ld,ierror)
call checkmpi(ierror,'Broadcasting rorb_angle in broadcast_binary.')
call mpi_bcast(xyzzyaadd1,xyzzyaacl1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting rorb_angle_type in broadcast_binary&
&.')
call mpi_bcast(xyzzyaacw1,xyzzyaacl1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting Ridx2Zidx in broadcast_binary.')
if(am_smpmaster)then
if(single_precision_blips)then
if(xyzzyaabk1)then
call mpi_bcast_safe(xyzzyaaao1,int(avc_norb,i64)*xyzzyaaaa1(1)*xyzzyaa&
&aa1(2)*xyzzyaaaa1(3),mpi_real,0,xyzzyaaao2,ierror)
call checkmpi(ierror,'Broadcasting avc_sp in broadcast_binary.')
else
call mpi_bcast_safe(xyzzyaaaq1,avc_norb*xyzzyaaaa1(1)*xyzzyaaaa1(2)*xy&
&zzyaaaa1(3),mpi_complex,0,xyzzyaaao2,ierror)
call checkmpi(ierror,'Broadcasting cavc_sp in broadcast_binary.')
endif
else
if(xyzzyaabk1)then
call mpi_bcast_safe(xyzzyaaap1,avc_norb*xyzzyaaaa1(1)*xyzzyaaaa1(2)*xy&
&zzyaaaa1(3),mpi_double_precision,0,xyzzyaaao2,ierror)
call checkmpi(ierror,'Broadcasting avc_sp in broadcast_binary.')
else
call mpi_bcast_safe(xyzzyaaar1,avc_norb*xyzzyaaaa1(1)*xyzzyaaaa1(2)*xy&
&zzyaaaa1(3),mpi_double_complex,0,xyzzyaaao2,ierror)
call checkmpi(ierror,'Broadcasting cavc_dp in broadcast_binary.')
endif
endif
endif
endif
call shallocate_barrier
end subroutine xyzzyaabv2
subroutine xyzzyaabw2
implicit none
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11,xyzzyaaae11,xy&
&zzyaaaf11,xyzzyaaag11,xyzzyaaah11,xyzzyaaai11,xyzzyaaaj11,xyzzyaaak11&
&,xyzzyaaal11,xyzzyaaam11,xyzzyaaan11(nspin),xyzzyaaao11,xyzzyaaap11(n&
&spin)
integer,allocatable :: xyzzyaaaq11(:,:),xyzzyaaar11(:)
real(dp) xyzzyaaas11,xyzzyaaat11
complex(dp) xyzzyaaau11
complex(dp),allocatable :: xyzzyaaav11(:)
logical xyzzyaaaw11
allocate(xyzzyaacm1(nemax,nspin,ndet),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP',' bw_orbmap')
allocate(xyzzyaaaw1(ndet,xyzzyaaag1,xyzzyaaaf1,nspin),stat=xyzzyaaaa11&
&)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP',' bandphase')
xyzzyaacm1=0
xyzzyaaaw1=999.d0
xyzzyaaci1=0
bw_norbc=0
xyzzyaacl1=0
xyzzyaack1=0
if(xyzzyaabn1)then
xyzzyaaak11=xyzzyaaag1*xyzzyaaaf1*nspin
xyzzyaaal11=2*xyzzyaaak11*ndet
allocate(xyzzyaacx1(xyzzyaaak11),xyzzyaacy1(xyzzyaaak11),xyzzyaacz1(xy&
&zzyaaak11),xyzzyaade1(xyzzyaaal11),xyzzyaadd1(xyzzyaaal11),xyzzyaacw1&
&(xyzzyaaal11),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','corb_*_idx,...')
xyzzyaacx1=0
xyzzyaacy1=0
xyzzyaacz1=0
xyzzyaade1=czero
xyzzyaadd1=0
xyzzyaacw1=0
xyzzyaaae1=0
xyzzyaaam11=0
do xyzzyaaae11=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaae11))then
do xyzzyaaac11=1,xyzzyaaad2
do xyzzyaaaf11=1,xyzzyaaah1(xyzzyaaae11,xyzzyaaac11)
if((xyzzyaabl1.and.any(xyzzyaaai1(1:ndet,xyzzyaaaf11,xyzzyaaae11,xyzzy&
&aaac11)>0)).or.(.not.xyzzyaabl1.and.any(xyzzyaaai1(1:ndet,xyzzyaaaf11&
&,xyzzyaaae11,1:nspin)>0)))xyzzyaaam11=xyzzyaaam11+1
enddo
enddo
endif
enddo
endif
allocate(xyzzyaaaq11(ndet,nspin),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','iorb_count')
xyzzyaaaq11=0
xyzzyaaai11=0
if(xyzzyaabl1)then
xyzzyaaaj11=0
xyzzyaaao11=0
else
xyzzyaaan11=0
xyzzyaaap11=0
endif
do xyzzyaaae11=1,xyzzyaaaf1
if(xyzzyaabm1)then
do xyzzyaaab11=1,nspin
if(xyzzyaaab11==2.and.xyzzyaabl1)then
xyzzyaaac11=2
else
xyzzyaaac11=1
endif
do xyzzyaaad11=1,xyzzyaaad1(xyzzyaaac11)
if(xyzzyaabl1)then
xyzzyaaao11=xyzzyaaao11+1
else
xyzzyaaap11(xyzzyaaab11)=xyzzyaaap11(xyzzyaaab11)+1
xyzzyaaao11=xyzzyaaap11(xyzzyaaab11)
endif
xyzzyaacm1(xyzzyaaad11,xyzzyaaab11,1)=xyzzyaaao11
enddo
do xyzzyaaag11=2,ndet
xyzzyaacm1(:,xyzzyaaab11,xyzzyaaag11)=xyzzyaacm1(:,xyzzyaaab11,1)
enddo
xyzzyaaaq11(:,xyzzyaaab11)=xyzzyaaad1(xyzzyaaac11)
enddo
endif
if(xyzzyaabn1)then
if(xyzzyaabp1(xyzzyaaae11))then
do xyzzyaaab11=1,nspin
if(xyzzyaaab11==2.and.xyzzyaabl1)then
xyzzyaaac11=2
else
xyzzyaaac11=1
endif
do xyzzyaaaf11=1,xyzzyaaah1(xyzzyaaae11,xyzzyaaac11)
if((xyzzyaabl1.and.all(xyzzyaaai1(1:ndet,xyzzyaaaf11,xyzzyaaae11,xyzzy&
&aaab11)==0)).or.(.not.xyzzyaabl1.and.all(xyzzyaaai1(1:ndet,xyzzyaaaf1&
&1,xyzzyaaae11,1:nspin)==0)))cycle
xyzzyaaae1(xyzzyaaab11)=xyzzyaaae1(xyzzyaaab11)+1
if(xyzzyaabl1)then
xyzzyaaaj11=xyzzyaaaj11+1
else
xyzzyaaan11(xyzzyaaab11)=xyzzyaaan11(xyzzyaaab11)+1
xyzzyaaaj11=xyzzyaaan11(xyzzyaaab11)
xyzzyaaao11=xyzzyaaap11(xyzzyaaab11)
endif
xyzzyaacx1(xyzzyaaaj11)=xyzzyaaaf11
xyzzyaacy1(xyzzyaaaj11)=xyzzyaaae11
xyzzyaacz1(xyzzyaaaj11)=xyzzyaaac11
if(xyzzyaabo1.or.xyzzyaabk1)then
do xyzzyaaag11=1,ndet
if(xyzzyaaai1(xyzzyaaag11,xyzzyaaaf11,xyzzyaaae11,xyzzyaaab11)==1)then
xyzzyaaaq11(xyzzyaaag11,xyzzyaaab11)=xyzzyaaaq11(xyzzyaaag11,xyzzyaaab&
&11)+1
xyzzyaacm1(xyzzyaaaq11(xyzzyaaag11,xyzzyaaab11),xyzzyaaab11,xyzzyaaag1&
&1)=xyzzyaaao11+xyzzyaaaj11
endif
enddo
if(xyzzyaabk1)then
xyzzyaaai11=xyzzyaaaj11
xyzzyaacw1(xyzzyaaaj11)=xyzzyaaaj11
endif
else
do xyzzyaaag11=1,ndet
if(xyzzyaaai1(xyzzyaaag11,xyzzyaaaf11,xyzzyaaae11,xyzzyaaab11)==2)then
xyzzyaaaw11=.false.
do xyzzyaaah11=1,xyzzyaaai11
if(xyzzyaacw1(xyzzyaaah11)==xyzzyaaaj11 .and.xyzzyaadd1(xyzzyaaah11)==&
&xyzzyaada1)then
xyzzyaaaw11=.true.
exit
endif
enddo
if(.not.xyzzyaaaw11)then
xyzzyaaai11=xyzzyaaai11+1
xyzzyaaah11=xyzzyaaai11
xyzzyaade1(xyzzyaaah11)=c_one
xyzzyaadd1(xyzzyaaah11)=xyzzyaada1
xyzzyaacw1(xyzzyaaah11)=xyzzyaaaj11
endif
xyzzyaaaq11(xyzzyaaag11,xyzzyaaab11)=xyzzyaaaq11(xyzzyaaag11,xyzzyaaab&
&11)+1
xyzzyaacm1(xyzzyaaaq11(xyzzyaaag11,xyzzyaaab11),xyzzyaaab11,xyzzyaaag1&
&1)=xyzzyaaao11+xyzzyaaah11
xyzzyaaaw11=.false.
do xyzzyaaah11=1,xyzzyaaai11
if(xyzzyaacw1(xyzzyaaah11)==xyzzyaaaj11 .and.xyzzyaadd1(xyzzyaaah11)==&
&xyzzyaadb1)then
xyzzyaaaw11=.true.
exit
endif
enddo
if(.not.xyzzyaaaw11)then
xyzzyaaai11=xyzzyaaai11+1
xyzzyaaah11=xyzzyaaai11
xyzzyaade1(xyzzyaaah11)=zi
xyzzyaadd1(xyzzyaaah11)=xyzzyaadb1
xyzzyaacw1(xyzzyaaah11)=xyzzyaaaj11
endif
xyzzyaaaq11(xyzzyaaag11,xyzzyaaab11)=xyzzyaaaq11(xyzzyaaag11,xyzzyaaab&
&11)+1
xyzzyaacm1(xyzzyaaaq11(xyzzyaaag11,xyzzyaaab11),xyzzyaaab11,xyzzyaaag1&
&1)=xyzzyaaao11+xyzzyaaah11
elseif(xyzzyaaai1(xyzzyaaag11,xyzzyaaaf11,xyzzyaaae11,xyzzyaaab11)==1)&
&then
call xyzzyaadt1(xyzzyaaaj11,xyzzyaaam11,xyzzyaaag11,xyzzyaaae11,xyzzya&
&aaf11,xyzzyaaac11,xyzzyaaaw1(xyzzyaaag11,xyzzyaaaf11,xyzzyaaae11,xyzz&
&yaaac11),xyzzyaaas11,xyzzyaaat11)
xyzzyaaaw11=.false.
do xyzzyaaah11=1,xyzzyaaai11
if(xyzzyaacw1(xyzzyaaah11)==xyzzyaaaj11 .and.dble(xyzzyaade1(xyzzyaaah&
&11))==xyzzyaaas11.and.aimag(xyzzyaade1(xyzzyaaah11))==xyzzyaaat11)the&
&n
xyzzyaaaw11=.true.
exit
endif
enddo
if(.not.xyzzyaaaw11)then
xyzzyaaai11=xyzzyaaai11+1
xyzzyaaah11=xyzzyaaai11
xyzzyaaau11=cmplx(xyzzyaaas11,xyzzyaaat11,dp)
if(abs(xyzzyaaas11)<1.d3*epsilon(1.d0))xyzzyaaau11=cmplx(0.d0,xyzzyaaa&
&t11,dp)
if(abs(xyzzyaaat11)<1.d3*epsilon(1.d0))xyzzyaaau11=cmplx(xyzzyaaas11,0&
&.d0,dp)
xyzzyaade1(xyzzyaaah11)=xyzzyaaau11
if(xyzzyaade1(xyzzyaaah11)==c_one)then
xyzzyaadd1(xyzzyaaah11)=xyzzyaada1
elseif(xyzzyaade1(xyzzyaaah11)==zi)then
xyzzyaadd1(xyzzyaaah11)=xyzzyaadb1
else
xyzzyaadd1(xyzzyaaah11)=xyzzyaadc1
endif
xyzzyaacw1(xyzzyaaah11)=xyzzyaaaj11
endif
xyzzyaaaq11(xyzzyaaag11,xyzzyaaab11)=xyzzyaaaq11(xyzzyaaag11,xyzzyaaab&
&11)+1
xyzzyaacm1(xyzzyaaaq11(xyzzyaaag11,xyzzyaaab11),xyzzyaaab11,xyzzyaaag1&
&1)=xyzzyaaao11+xyzzyaaah11
endif
enddo
endif
enddo
enddo
endif
endif
enddo
deallocate(xyzzyaaaq11)
xyzzyaack1=0
if(xyzzyaabn1)then
if(xyzzyaabl1)then
bw_norbc=xyzzyaaaj11
else
bw_norbc=maxval(xyzzyaaan11)
endif
xyzzyaacl1=xyzzyaaai11
if(xyzzyaabo1)then
xyzzyaack1=bw_norbc
else
xyzzyaack1=xyzzyaacl1
endif
endif
xyzzyaaci1=xyzzyaacj1+xyzzyaack1
if(xyzzyaabn1)then
allocate(xyzzyaaar11(bw_norbc),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','itrf')
xyzzyaaar11(1:bw_norbc)=xyzzyaacx1(1:bw_norbc)
deallocate(xyzzyaacx1)
allocate(xyzzyaacx1(bw_norbc),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','corb_band_idx')
xyzzyaacx1(1:bw_norbc)=xyzzyaaar11(1:bw_norbc)
deallocate(xyzzyaaar11)
allocate(xyzzyaaar11(bw_norbc),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','itrf')
xyzzyaaar11(1:bw_norbc)=xyzzyaacy1(1:bw_norbc)
deallocate(xyzzyaacy1)
allocate(xyzzyaacy1(bw_norbc),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','corb_kvec_idx')
xyzzyaacy1(1:bw_norbc)=xyzzyaaar11(1:bw_norbc)
deallocate(xyzzyaaar11)
allocate(xyzzyaaar11(bw_norbc),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'READBWF','itrf')
xyzzyaaar11(1:bw_norbc)=xyzzyaacz1(1:bw_norbc)
deallocate(xyzzyaacz1)
allocate(xyzzyaacz1(bw_norbc),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'READBWF','corb_spin_idx')
xyzzyaacz1(1:bw_norbc)=xyzzyaaar11(1:bw_norbc)
deallocate(xyzzyaaar11)
allocate(xyzzyaaav11(xyzzyaacl1),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','ztrf')
xyzzyaaav11(1:xyzzyaacl1)=xyzzyaade1(1:xyzzyaacl1)
deallocate(xyzzyaade1)
allocate(xyzzyaade1(xyzzyaacl1),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','rorb_angle')
xyzzyaade1(1:xyzzyaacl1)=xyzzyaaav11(1:xyzzyaacl1)
deallocate(xyzzyaaav11)
allocate(xyzzyaaar11(xyzzyaacl1),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','itrf')
xyzzyaaar11(1:xyzzyaacl1)=xyzzyaadd1(1:xyzzyaacl1)
deallocate(xyzzyaadd1)
allocate(xyzzyaadd1(xyzzyaacl1),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','rorb_angle_type')
xyzzyaadd1(1:xyzzyaacl1)=xyzzyaaar11(1:xyzzyaacl1)
deallocate(xyzzyaaar11)
allocate(xyzzyaaar11(xyzzyaacl1),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','itrf')
xyzzyaaar11(1:xyzzyaacl1)=xyzzyaacw1(1:xyzzyaacl1)
deallocate(xyzzyaacw1)
allocate(xyzzyaacw1(xyzzyaacl1),stat=xyzzyaaaa11)
call check_alloc(xyzzyaaaa11,'DO_ORBMAP','Ridx2Zidx')
xyzzyaacw1(1:xyzzyaacl1)=xyzzyaaar11(1:xyzzyaacl1)
deallocate(xyzzyaaar11)
endif
if(no_orb_phases>0)deallocate(orb_phase_det,orb_phase_band,orb_phase_k&
&point,orb_phase_spin,orb_phase)
end subroutine xyzzyaabw2
end subroutine readbwf
subroutine bwfdet_setup(band_occup_only_in,skip_orbangles_in)
use slaarnaan,only : check_kpoints
use slaarnabi,         only : use_gpcc
use slaarnabt,    only : quicksort,get_numerical_orbmask,get_numerical&
&_orbrange
implicit none
logical,intent(in),optional :: band_occup_only_in,skip_orbangles_in
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12,xy&
&zzyaaaf12,i,xyzzyaaag12,xyzzyaaah12,xyzzyaaai12,xyzzyaaaj12,xyzzyaaak&
&12,xyzzyaaal12,xyzzyaaam12,xyzzyaaan12,xyzzyaaao12,xyzzyaaap12,xyzzya&
&aaq12,xyzzyaaar12,xyzzyaaas12,xyzzyaaat12,xyzzyaaau12,xyzzyaaav12
integer,allocatable :: xyzzyaaaw12(:),xyzzyaaax12(:),xyzzyaaay12(:,:),&
&xyzzyaaaz12(:,:,:)
real(dp) xyzzyaaba12(3),xyzzyaabb12,xyzzyaabc12,xyzzyaabd12,xyzzyaabe1&
&2(3),xyzzyaabf12(1:3),xyzzyaabg12,xyzzyaabh12(3),xyzzyaabi12(3),xyzzy&
&aabj12(3),xyzzyaabk12(3)
real(dp),parameter :: xyzzyaabl12=1.d-8
real(dp),parameter :: xyzzyaabm12=0.1d0
real(dp),parameter :: xyzzyaabn12=1.d0
real(dp),allocatable :: xyzzyaabo12(:),xyzzyaabp12(:,:)
logical xyzzyaabq12,band_occup_only,xyzzyaabr12,xyzzyaabs12(2),xyzzyaa&
&bt12,xyzzyaabu12
logical,parameter :: xyzzyaabv12=.true.
logical,allocatable :: ltemp(:),xyzzyaabw12(:),xyzzyaabx12(:)
character(80) char80,char80_new
character (256) tmpr
if(present(band_occup_only_in))then
band_occup_only=band_occup_only_in
else
band_occup_only=.false.
endif
if(present(skip_orbangles_in))then
xyzzyaabr12=skip_orbangles_in
else
xyzzyaabr12=.false.
endif
if(.not.band_occup_only)xyzzyaabr12=.false.
if(periodicity<=2.and.am_master.and..not.band_occup_only)then
xyzzyaaar12=sum(atno)
if(xyzzyaaar12>0)then
xyzzyaabk12(3)=xyzzyaabn12*sqrt(dot_product(painv(1:3,3),painv(1:3,3))&
&)
if(periodicity<=1)then
xyzzyaabk12(2)=xyzzyaabn12*sqrt(dot_product(painv(1:3,2),painv(1:3,2))&
&)
if(periodicity==0)xyzzyaabk12(1)=xyzzyaabn12*sqrt(dot_product(painv(1:&
&3,1),painv(1:3,1)))
endif
xyzzyaabi12=0.d0
do i=1,nbasis
xyzzyaabj12=matmul(basis(1:3,i),painv)
if(xyzzyaabj12(3)>1.d0-xyzzyaabk12(3).or.xyzzyaabj12(3)<xyzzyaabk12(3)&
&)call errstop('BWFDET_SETUP','In a reduced-periodicity blip calculati&
&on, your atom coordinates need to be in the centre of the region span&
&ned by the blip grid, which is the parallelepiped defined by the latt&
&ice vectors starting from the origin.  The problem occurs for atom ' &
&//trim(i2s(i))//' in your basis (direction of third lattice vector).'&
&)
if(periodicity<=1)then
if(xyzzyaabj12(2)>1.d0-xyzzyaabk12(2).or.xyzzyaabj12(2)<xyzzyaabk12(2)&
&)call errstop('BWFDET_SETUP','In a reduced-periodicity blip calculati&
&on, your atom coordinates need to be in the centre of the region span&
&ned by the blip grid, which is the parallelepiped defined by the latt&
&ice vectors starting from the origin.  The problem occurs for atom '/&
&/trim(i2s(i))//' in your basis (direction of second lattice vector).'&
&)
if(periodicity==0.and.xyzzyaabj12(1)>1.d0-xyzzyaabk12(1).or.xyzzyaabj1&
&2(1)<xyzzyaabk12(1))call errstop('BWFDET_SETUP','In a reduced-periodi&
&city blip calculation, your atom coordinates need to be in the centre&
& of the region spanned by the blip grid, which is the parallelepiped &
&defined by the lattice vectors starting from the origin.  The problem&
& occurs for atom '//trim(i2s(i)) //' in your basis (direction of firs&
&t lattice vectors).')
endif
xyzzyaabi12=xyzzyaabi12+xyzzyaabj12*dble(atno(i))
enddo
xyzzyaabi12=xyzzyaabi12/dble(xyzzyaaar12)
if(abs(xyzzyaabi12(3)-0.5d0)>xyzzyaabm12)call errwarn('BWFDET_SETUP','&
&The atomic coordinates in your blip file should be close to the centr&
&e of the cell defined by the lattice vectors with the origin as one o&
&f its vertices.  The problem is in the direction of the third lattice&
& vector.')
if(periodicity<=1)then
if(abs(xyzzyaabi12(2)-0.5d0)>xyzzyaabm12)call errwarn('BWFDET_SETUP','&
&The atomic coordinates in your blip file should be close to the centr&
&e of the cell defined by the lattice vectors with the origin as one o&
&f its vertices.  The problem is in the direction of the second lattic&
&e vector.')
if(periodicity==0.and.abs(xyzzyaabi12(1)-0.5d0)>xyzzyaabm12)call errwa&
&rn('BWFDET_SETUP','The atomic coordinates in your blip file should be&
& close to the centre of the cell defined by the lattice vectors with &
&the origin as one of its vertices.  The problem is in the direction o&
&f the first lattice vector.')
endif
endif
endif
allocate(xyzzyaaai1(ndet,xyzzyaaag1,xyzzyaaaf1,nspin),xyzzyaabp1(xyzzy&
&aaaf1),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','0')
xyzzyaaai1=0
xyzzyaabp1=.true.
if(xyzzyaabl1)then
xyzzyaaae12=2
else
xyzzyaaae12=1
endif
if(xyzzyaabo1)then
if(.not.band_occup_only)then
allocate(xyzzyaabx12(xyzzyaaaf1),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','TIME_REVERSED')
call check_kpoints(xyzzyaaaf1,xyzzyaaay1,xyzzyaabx12)
if(any(xyzzyaabx12))then
if(am_master)call wout('Flipping time-reversed k points in blip wave f&
&unction.')
do xyzzyaaan12=1,xyzzyaaae12
do xyzzyaaaq12=1,bw_norbc
xyzzyaaab12=xyzzyaacy1(xyzzyaaaq12)
if(xyzzyaabx12(xyzzyaaab12))then
if(single_precision_blips)then
xyzzyaaaq1(xyzzyaaaq12,:,:,:)=conjg(xyzzyaaaq1(xyzzyaaaq12,:,:,:))
else
xyzzyaaar1(xyzzyaaaq12,:,:,:)=conjg(xyzzyaaar1(xyzzyaaaq12,:,:,:))
endif
endif
enddo
enddo
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabx12(xyzzyaaab12))xyzzyaaay1(:,xyzzyaaab12)=-xyzzyaaay1(:,xy&
&zzyaaab12)
enddo
call check_kpoints(xyzzyaaaf1,xyzzyaaay1)
endif
deallocate(xyzzyaabx12)
endif
if(xyzzyaabl1)then
i=max(sum(xyzzyaaah1(:,1)),sum(xyzzyaaah1(:,2)))
else
i=sum(xyzzyaaah1(:,1))
endif
allocate(xyzzyaabo12(i),xyzzyaaax12(i),xyzzyaaaw12(i),stat=xyzzyaaaa12&
&)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','EIGTEMP')
xyzzyaabs12=.false.
do xyzzyaaan12=1,nspin
if(xyzzyaaan12==2.and.xyzzyaabl1)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
if(xyzzyaaan12==1.or.xyzzyaaad12==2)then
xyzzyaaag12=0
do xyzzyaaab12=1,xyzzyaaaf1
do xyzzyaaaf12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
xyzzyaaag12=xyzzyaaag12+1
xyzzyaabo12(xyzzyaaag12)=xyzzyaaax1(xyzzyaaaf12,xyzzyaaab12,xyzzyaaad1&
&2)
xyzzyaaax12(xyzzyaaag12)=xyzzyaaab12
enddo
enddo
if(xyzzyaaag12>0)call quicksort(xyzzyaaag12,xyzzyaabo12(1),xyzzyaaaw12&
&(1))
endif
if(nuc_nele(xyzzyaaan12)>xyzzyaaag12)then
if(am_master)call errstop('BWFDET_SETUP','An insufficient number of ei&
&genstates have been supplied.')
call qmc_barrier
endif
do i=1,nuc_nele(xyzzyaaan12)
xyzzyaaab12=xyzzyaaax12(xyzzyaaaw12(i))
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(1,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)==0)then
xyzzyaaai1(1,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)=1
exit
endif
enddo
enddo
if(am_master.and.nuc_nele(xyzzyaaan12)<xyzzyaaag12.and.nuc_nele(xyzzya&
&aan12)>0)then
if(abs(xyzzyaabo12(xyzzyaaaw12(nuc_nele(xyzzyaaan12)))-xyzzyaabo12(xyz&
&zyaaaw12(nuc_nele(xyzzyaaan12)+1)))<xyzzyaabl12*abs(xyzzyaabo12(xyzzy&
&aaaw12(nuc_nele(xyzzyaaan12)))).and.ndet==1)xyzzyaabs12(xyzzyaaan12)=&
&.true.
endif
do xyzzyaaao12=2,ndet
xyzzyaaai1(xyzzyaaao12,:,:,xyzzyaaan12)=xyzzyaaai1(1,:,:,xyzzyaaan12)
enddo
enddo
if(am_master.and..not.band_occup_only)then
if(all(xyzzyaabs12))then
call errwarn('BWFDET_SETUP','Partially occupied degenerate states at t&
&he Fermi level for both spins. Multideterminant calculation or twist-&
&averaging probably advisable.')
elseif(any(xyzzyaabs12))then
call errwarn('BWFDET_SETUP','Partially occupied degenerate states at t&
&he Fermi level for spin '//trim(i2s(maxval((/1,2/),xyzzyaabs12)))//'.&
& Multideterminant calculation or twist-averaging probably advisable.'&
&)
endif
call wout()
call wout('Blip setup')
call wout('==========')
call wout()
call wout('Periodicity : '//trim(i2s(periodicity)))
call wout('Complex blip orbitals will be used.')
call wout()
if(periodicity>0)then
allocate(xyzzyaaay12(xyzzyaaaf1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','boccband')
xyzzyaaay12=0
do xyzzyaaan12=1,nspin
if(xyzzyaaan12==2.and.xyzzyaabl1)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
do xyzzyaaab12=1,xyzzyaaaf1
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(1,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)>0)xyzzyaaay12(xyz&
&zyaaab12,xyzzyaaan12)=xyzzyaaay12(xyzzyaaab12,xyzzyaaan12)+1
enddo
enddo
enddo
if((xyzzyaabl1.or.nuc_nele(1)==nuc_nele(2)).and.all(wf_nd(1,1:nspin)==&
&0))then
xyzzyaabq12=.false.
sc: do xyzzyaaan12=1,nspin
do xyzzyaaab12=2,xyzzyaaaf1
if(xyzzyaaay12(xyzzyaaab12,xyzzyaaan12)/=xyzzyaaay12(1,xyzzyaaan12))th&
&en
xyzzyaabq12=.true.
exit sc
endif
enddo
enddo sc
if(xyzzyaabq12)then
call wout('METALLIC GROUND STATE DETECTED')
if(.not.xyzzyaabl1)then
call wout('Number of doubly-occupied bands filled at each k point in G&
&S:')
call wout('k-point    no. of filled bands')
do xyzzyaaab12=1,xyzzyaaaf1
write(tmpr,'(i4,4x,i5)')xyzzyaaab12,xyzzyaaay12(xyzzyaaab12,1)
call wout(tmpr)
enddo
else
call wout('Number of singly-occupied bands filled at each k point in G&
&S :')
call wout('Spin    k-point    no. of filled bands')
do xyzzyaaab12=1,xyzzyaaaf1
write(tmpr,'(a8,i4,4x,i5)')'up',xyzzyaaab12,xyzzyaaay12(xyzzyaaab12,1)
call wout(tmpr)
enddo
do xyzzyaaab12=1,xyzzyaaaf1
write(tmpr,'(a8,i4,4x,i5)')'down',xyzzyaaab12,xyzzyaaay12(xyzzyaaab12,&
&2)
call wout(tmpr)
enddo
endif
else
call wout('INSULATING GROUND STATE DETECTED')
if(.not.xyzzyaabl1)then
call wout('No. of doubly occupied bands at each k point: ' //trim(i2s(&
&xyzzyaaay12(1,1))))
else
call wout('No. of occupied spin-up bands at each k point  : ' //trim(i&
&2s(xyzzyaaay12(1,1))))
call wout('No. of occupied spin-down bands at each k point: ' //trim(i&
&2s(xyzzyaaay12(1,2))))
endif
endif
call wout()
endif
deallocate(xyzzyaaay12)
endif
do xyzzyaaab12=1,xyzzyaaaf1
if(periodicity<3.and.abs(dot_product(xyzzyaaay1(1:3,xyzzyaaab12),pa3(1&
&:3)))>xyzzyaabl12)call errstop('BWFDET_SETUP','k vectors inappropriat&
&e for periodicity.')
if(periodicity<2.and.abs(dot_product(xyzzyaaay1(1:3,xyzzyaaab12),pa2(1&
&:3)))>xyzzyaabl12)call errstop('BWFDET_SETUP','k vectors inappropriat&
&e for periodicity.')
if(periodicity<1.and.abs(dot_product(xyzzyaaay1(1:3,xyzzyaaab12),pa1(1&
&:3)))>xyzzyaabl12)call errstop('BWFDET_SETUP','k vectors inappropriat&
&e for periodicity.')
enddo
endif
do xyzzyaaan12=1,nspin
if(xyzzyaaan12==2.and.xyzzyaabl1)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
do xyzzyaaao12=1,ndet
if((wf_np(xyzzyaaao12,xyzzyaaan12)>0.or.wf_nm(xyzzyaaao12,xyzzyaaan12)&
&>0).and.am_master)call errstop('BWFDET_SETUP','Additions or subtracti&
&ons should not be used for blip orbitals: number of electrons is simp&
&ly determined by NEU and NED, and occupancy should only be changed us&
&ing "PR" in the MDET block.')
do xyzzyaaah12=1,wf_nd(xyzzyaaao12,xyzzyaaan12)
xyzzyaaak12=wf_d(1,xyzzyaaah12,xyzzyaaao12,xyzzyaaan12)
xyzzyaaai12=wf_d(2,xyzzyaaah12,xyzzyaaao12,xyzzyaaan12)
xyzzyaaal12=wf_d(3,xyzzyaaah12,xyzzyaaao12,xyzzyaaan12)
xyzzyaaaj12=wf_d(4,xyzzyaaah12,xyzzyaaao12,xyzzyaaan12)
if(xyzzyaaaj12<1.or.xyzzyaaaj12>xyzzyaaaf1.or.xyzzyaaai12<1.or.xyzzyaa&
&ai12>xyzzyaaaf1)then
if(am_master)call errstop('BWFDET_SETUP','k-point out of range in exci&
&tation specification.')
call qmc_barrier
endif
if(xyzzyaaal12<1.or.xyzzyaaal12>xyzzyaaah1(xyzzyaaaj12,xyzzyaaad12).or&
&.xyzzyaaak12<1.or.xyzzyaaak12>xyzzyaaah1(xyzzyaaai12,xyzzyaaad12))the&
&n
if(am_master)call errstop('BWFDET_SETUP','band index out of range in e&
&xcitation specification.')
call qmc_barrier
endif
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaal12,xyzzyaaaj12,xyzzyaaan12)/=0)then
if(am_master)call errstop('BWFDET_SETUP','Trying to promote an electro&
&n into a band that is already fully occupied.')
call qmc_barrier
endif
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaak12,xyzzyaaai12,xyzzyaaan12)/=1)then
if(am_master)call errstop('BWFDET_SETUP','Trying to promote an electro&
&n from an empty band.')
call qmc_barrier
endif
xyzzyaaai1(xyzzyaaao12,xyzzyaaal12,xyzzyaaaj12,xyzzyaaan12)=1
xyzzyaaai1(xyzzyaaao12,xyzzyaaak12,xyzzyaaai12,xyzzyaaan12)=0
enddo
i=0
do xyzzyaaab12=1,xyzzyaaaf1
i=i+sum(xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12))
enddo
if(i/=nuc_nele(xyzzyaaan12))then
if(am_master)call errstop('BWFDET_SETUP','Problem with excitations (A)&
&.')
call qmc_barrier
endif
if(any(xyzzyaaai1(xyzzyaaao12,:,:,xyzzyaaan12)<0).or.any(xyzzyaaai1(xy&
&zzyaaao12,:,:,xyzzyaaan12)>1))then
if(am_master)call errstop('BWFDET_SETUP','Problem with excitations (B)&
&.')
call qmc_barrier
endif
enddo
enddo
deallocate(xyzzyaaaw12,xyzzyaabo12,xyzzyaaax12)
if(band_occup_only)return
if(am_master.and.xyzzyaabv12)then
call wout('Detailed band-occupancy information')
call wout('-----------------------------------')
do xyzzyaaab12=1,xyzzyaaaf1
call wout('K POINT '//trim(i2s(xyzzyaaab12)))
if(isperiodic)then
write(tmpr,'("k(au):",3(1x,es23.15))')xyzzyaaay1(1:periodicity,xyzzyaa&
&ab12)
call wout(tmpr)
write(tmpr,'("kfrac:",3(1x,es23.15))')matmul(xyzzyaaay1(1:periodicity,&
&xyzzyaaab12),pbinv(1:periodicity,1:periodicity))
call wout(tmpr)
endif
do xyzzyaaan12=1,nspin
if(nuc_nele(xyzzyaaan12)==0)cycle
if(xyzzyaabl1.and.xyzzyaaan12==2)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
do xyzzyaaao12=1,ndet
if(ndet>1)then
call wout('Spin '//trim(i2s(xyzzyaaan12))//', determinant: '//trim(i2s&
&(xyzzyaaao12))//':')
else
call wout('Spin '//trim(i2s(xyzzyaaan12))//':')
endif
if(any(xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)==1))then
char80=' Bands occupied at k:'
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)==1)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam12))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam12))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Highest occupied band energy: ',maxval(xyzzyaaax1(:,xyzzyaa&
&ab12,xyzzyaaad12),xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)==&
&1))
endif
if(any(xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)==0))then
char80=' Bands unoccupied at k:'
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)==0)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam12))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam12))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Lowest unoccupied band energy: ',minval(xyzzyaaax1(:,xyzzya&
&aab12,xyzzyaaad12),xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)=&
&=0))
endif
enddo
enddo
call wout()
enddo
endif
else
allocate(xyzzyaabw12(xyzzyaaaf1),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','2')
if(.not.allocated(xyzzyaabq1))allocate(xyzzyaabq1(xyzzyaaaf1),stat=xyz&
&zyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','2')
if(.not.band_occup_only)call check_kpoints(xyzzyaaaf1,xyzzyaaay1)
xyzzyaabw12=.false.
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
do xyzzyaaac12=xyzzyaaab12+1,xyzzyaaaf1
xyzzyaaba12(1:3)=xyzzyaaay1(1:3,xyzzyaaab12)+xyzzyaaay1(1:3,xyzzyaaac1&
&2)
xyzzyaabg12=dot_product(xyzzyaaba12,pa1)*one_over_twopi
if(abs(xyzzyaabg12-anint(xyzzyaabg12))<xyzzyaabl12)then
xyzzyaabg12=dot_product(xyzzyaaba12,pa2)*one_over_twopi
if(abs(xyzzyaabg12-anint(xyzzyaabg12))<xyzzyaabl12)then
xyzzyaabg12=dot_product(xyzzyaaba12,pa3)*one_over_twopi
if(abs(xyzzyaabg12-anint(xyzzyaabg12))<xyzzyaabl12)then
xyzzyaabp1(xyzzyaaac12)=.false.
xyzzyaabw12(xyzzyaaab12)=.true.
xyzzyaabw12(xyzzyaaac12)=.true.
exit
endif
endif
endif
enddo
endif
enddo
xyzzyaabq1=.false.
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
xyzzyaaba12(1:3)=2.d0*xyzzyaaay1(1:3,xyzzyaaab12)
xyzzyaabg12=dot_product(xyzzyaaba12,pa1)*one_over_twopi
if(abs(xyzzyaabg12-anint(xyzzyaabg12))<xyzzyaabl12)then
xyzzyaabg12=dot_product(xyzzyaaba12,pa2)*one_over_twopi
if(abs(xyzzyaabg12-anint(xyzzyaabg12))<xyzzyaabl12)then
xyzzyaabg12=dot_product(xyzzyaaba12,pa3)*one_over_twopi
if(abs(xyzzyaabg12-anint(xyzzyaabg12))<xyzzyaabl12)xyzzyaabq1(xyzzyaaa&
&b12)=.true.
endif
endif
endif
enddo
if(xyzzyaabl1)then
i=2*max(sum(xyzzyaaah1(:,1)),sum(xyzzyaaah1(:,2)))
else
i=2*sum(xyzzyaaah1(:,1))
endif
allocate(xyzzyaabo12(i),xyzzyaaax12(i),xyzzyaaaw12(i),stat=xyzzyaaaa12&
&)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','eigtemp')
do xyzzyaaan12=1,nspin
if(xyzzyaaan12==2.and.xyzzyaabl1)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
if(xyzzyaaan12==1.or.xyzzyaaad12==2)then
xyzzyaaag12=0
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
do xyzzyaaaf12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
xyzzyaaag12=xyzzyaaag12+1
xyzzyaabo12(xyzzyaaag12)=xyzzyaaax1(xyzzyaaaf12,xyzzyaaab12,xyzzyaaad1&
&2)
xyzzyaaax12(xyzzyaaag12)=xyzzyaaab12
if(.not.xyzzyaabq1(xyzzyaaab12))then
xyzzyaaag12=xyzzyaaag12+1
xyzzyaabo12(xyzzyaaag12)=xyzzyaaax1(xyzzyaaaf12,xyzzyaaab12,xyzzyaaad1&
&2)
xyzzyaaax12(xyzzyaaag12)=xyzzyaaab12
endif
enddo
endif
enddo
if(xyzzyaaag12>0)call quicksort(xyzzyaaag12,xyzzyaabo12(1),xyzzyaaaw12&
&(1))
endif
if(nuc_nele(xyzzyaaan12)>xyzzyaaag12+xyzzyaaad1(xyzzyaaad12))then
if(am_master)call errstop('BWFDET_SETUP','An insufficient number of ei&
&genstates have been supplied.')
call qmc_barrier
endif
do i=1,nuc_nele(xyzzyaaan12)-xyzzyaaad1(xyzzyaaad12)
xyzzyaaab12=xyzzyaaax12(xyzzyaaaw12(i))
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if((xyzzyaabq1(xyzzyaaab12).and.xyzzyaaai1(1,xyzzyaaam12,xyzzyaaab12,x&
&yzzyaaan12)<1).or.(.not.xyzzyaabq1(xyzzyaaab12).and.xyzzyaaai1(1,xyzz&
&yaaam12,xyzzyaaab12,xyzzyaaan12)<2))then
xyzzyaaai1(1,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)=xyzzyaaai1(1,xyzzyaa&
&am12,xyzzyaaab12,xyzzyaaan12)+1
exit
endif
enddo
enddo
if(.not.band_occup_only)then
if(am_master.and.nuc_nele(xyzzyaaan12)<xyzzyaaag12+xyzzyaaad1(xyzzyaaa&
&d12).and.nuc_nele(xyzzyaaan12)>xyzzyaaad1(xyzzyaaad12))then
if(abs(xyzzyaabo12(xyzzyaaaw12(nuc_nele(xyzzyaaan12)-xyzzyaaad1(xyzzya&
&aad12)))-xyzzyaabo12(xyzzyaaaw12(nuc_nele(xyzzyaaan12)+1-xyzzyaaad1(x&
&yzzyaaad12))))<xyzzyaabl12*abs(xyzzyaabo12(xyzzyaaaw12(nuc_nele(xyzzy&
&aaan12)-xyzzyaaad1(xyzzyaaad12)))).and.ndet==1)call errwarn('BWFDET_S&
&ETUP','Partially occupied degenerate states at the Fermi level for sp&
&in '//trim(i2s(xyzzyaaan12))//'. Multideterminant calculation probabl&
&y advisable.')
endif
endif
do xyzzyaaao12=2,ndet
xyzzyaaai1(xyzzyaaao12,:,:,xyzzyaaan12)=xyzzyaaai1(1,:,:,xyzzyaaan12)
enddo
enddo
if(am_master.and..not.band_occup_only)then
call wout()
call wout('Blip setup')
call wout('==========')
call wout()
call wout('Periodicity : '//trim(i2s(periodicity)))
call wout('Real blip orbitals will be used.')
call wout()
do xyzzyaaan12=1,nspin
if(xyzzyaaan12==2.and.xyzzyaabl1)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
if(xyzzyaaad1(xyzzyaaad12)>nuc_nele(xyzzyaaan12))call errstop('BWFDET_&
&SETUP','Some localized orbitals are unoccupied.')
enddo
if(xyzzyaabm1)then
call wordwrap('Localized orbitals are present.  Note that localized or&
&bitals are not included in the list of bands.  Band indices start at &
&the first non-localized band.')
call wout()
endif
if(periodicity>0)then
allocate(xyzzyaaay12(xyzzyaaaf1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','boccband')
xyzzyaaay12=0
do xyzzyaaan12=1,nspin
if(xyzzyaaan12==2.and.xyzzyaabl1)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(1,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)>0)xyzzyaaay12(xyz&
&zyaaab12,xyzzyaaan12)=xyzzyaaay12(xyzzyaaab12,xyzzyaaan12)+1
enddo
endif
enddo
enddo
xyzzyaabu12=(xyzzyaabl1.or.nuc_nele(1)==nuc_nele(2)).and.all(wf_nd(1,1&
&:nspin)==0)
if(xyzzyaabu12)then
xyzzyaabq12=.false.
sr:  do xyzzyaaan12=1,nspin
do xyzzyaaap12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaap12))exit
enddo
do xyzzyaaab12=xyzzyaaap12+1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
if(xyzzyaaay12(xyzzyaaab12,xyzzyaaan12)/=xyzzyaaay12(xyzzyaaap12,xyzzy&
&aaan12))then
xyzzyaabq12=.true.
exit sr
endif
endif
enddo
enddo sr
endif
call wout('K POINT ANALYSIS')
call wout('  k    kx         ky         kz       use pair edge')
do xyzzyaaab12=1,xyzzyaaaf1
write(tmpr,'(i3,1x,f10.6,1x,f10.6,1x,f10.6,3x,l1,4x,l1,4x,l1)')xyzzyaa&
&ab12,xyzzyaaay1(1,xyzzyaaab12),xyzzyaaay1(2,xyzzyaaab12),xyzzyaaay1(3&
&,xyzzyaaab12),xyzzyaabp1(xyzzyaaab12),xyzzyaabw12(xyzzyaaab12),xyzzya&
&abq1(xyzzyaaab12)
call wout(tmpr)
enddo
call wout()
call wout('Any k points with edge=F give rise to independent states at&
& both k and -k.')
call wout()
allocate(xyzzyaabp12(3,xyzzyaaaf1),ltemp(xyzzyaaaf1),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','kvec2')
xyzzyaabp12=1.d5
ltemp=.true.
do xyzzyaaac12=1,xyzzyaaaf1
do xyzzyaaab12=1,num_g
xyzzyaabf12(1:3)=xyzzyaaay1(1:3,xyzzyaaac12)-0.5d0*sr_lattice(1:3,xyzz&
&yaaab12)
xyzzyaabb12=one_over_twopi*dot_product(xyzzyaabf12,a1)
xyzzyaabc12=one_over_twopi*dot_product(xyzzyaabf12,a2)
xyzzyaabd12=one_over_twopi*dot_product(xyzzyaabf12,a3)
if(abs(xyzzyaabb12-anint(xyzzyaabb12))<xyzzyaabl12.and.abs(xyzzyaabc12&
&-anint(xyzzyaabc12))<xyzzyaabl12.and.abs(xyzzyaabd12-anint(xyzzyaabd1&
&2))<xyzzyaabl12)then
ltemp(xyzzyaaac12)=.false.
xyzzyaabp12(1:3,xyzzyaaac12)=0.5d0*sr_lattice(1:3,xyzzyaaab12)
exit
endif
enddo
enddo
if(any(ltemp))call errstop('BWFDET_SETUP','Supercell k_s vector not th&
&e Gamma point or half a supercell reciprocal lattice vector. Unable t&
&o make orbitals into real Bloch functions.')
xyzzyaabe12=xyzzyaabp12(1:3,1)
do xyzzyaaab12=2,xyzzyaaaf1
if(abs(xyzzyaabe12(1)-xyzzyaabp12(1,xyzzyaaab12))>xyzzyaabl12.or.abs(x&
&yzzyaabe12(2)-xyzzyaabp12(2,xyzzyaaab12))>xyzzyaabl12.or.abs(xyzzyaab&
&e12(3)-xyzzyaabp12(3,xyzzyaaab12))>xyzzyaabl12)then
call wout('Primitive cell reciprocal lattice vectors: (au)')
call wout('pb1',pb1(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('pb2',pb2(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('pb3',pb3(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('Simulation cell reciprocal lattice vectors: (au)')
call wout('b1 ',b1(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('b2 ',b2(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('b3 ',b3(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('Input k points:')
do xyzzyaaac12=1,xyzzyaaaf1
call wout(trim(i2s(xyzzyaaac12)),xyzzyaaay1(:,xyzzyaaac12),fmt='(1x,a,&
&t5,f16.10,1x,f16.10,1x,f16.10)')
enddo
call wout('Input k point set after reduction into 1st Brillouin zone o&
&f simulation cell:')
do xyzzyaaac12=1,xyzzyaaaf1
call wout(trim(i2s(xyzzyaaac12)),xyzzyaabp12(:,xyzzyaaac12),fmt='(1x,a&
&,t5,f16.10,1x,f16.10,1x,f16.10)')
enddo
call wout()
call wordwrap('The k point set in bwfn.data does not map to a unique s&
&imulation cell k_s vector.  Are you sure the NPCELL block in input is&
& correct?')
call wout()
call errstop('BWFDET_SETUP','Quitting')
endif
enddo
call wout('MAPPING ONTO UNIQUE K_S VECTOR:')
write(tmpr,'(3f12.8,a)')xyzzyaabe12(1:3),' (Cartesian a.u.)'
call wout(tmpr)
if(any(abs(xyzzyaabe12(1:3))>xyzzyaabl12))then
xyzzyaabb12=one_over_twopi*dot_product(xyzzyaabe12,a1)
xyzzyaabc12=one_over_twopi*dot_product(xyzzyaabe12,a2)
xyzzyaabd12=one_over_twopi*dot_product(xyzzyaabe12,a3)
write(tmpr,'(f12.8,1x,f12.8,1x,f12.8,1x,a)')xyzzyaabb12,xyzzyaabc12,xy&
&zzyaabd12,' (frac supercell recip. lattice vectors)'
call wout(tmpr)
endif
call wout()
deallocate(xyzzyaabp12,ltemp)
if(xyzzyaabu12)then
if(xyzzyaabq12)then
call wout('METALLIC GROUND STATE DETECTED')
if(.not.xyzzyaabl1)then
call wout('Number of doubly-occupied bands filled at each k point in G&
&S:')
call wout('k-point    no. of filled bands')
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
write(tmpr,'(i4,4x,i5)')xyzzyaaab12,xyzzyaaay12(xyzzyaaab12,1)
call wout(tmpr)
endif
enddo
else
call wout('Number of singly-occupied bands filled at each k point in G&
&S :')
call wout('Spin    k-point    no. of filled bands')
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
write(tmpr,'(a8,i4,4x,i5)')'up',xyzzyaaab12,xyzzyaaay12(xyzzyaaab12,1)
call wout(tmpr)
endif
enddo
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
write(tmpr,'(a8,i4,4x,i5)')'down',xyzzyaaab12,xyzzyaaay12(xyzzyaaab12,&
&2)
call wout(tmpr)
endif
enddo
endif
else
call wout('INSULATING GROUND STATE DETECTED')
if(.not.xyzzyaabl1)then
call wout('No. of doubly occupied bands at each k point: '//trim(i2s(x&
&yzzyaaay12(1,1))))
else
call wout('No. of occupied spin-up bands at each k point  : '//trim(i2&
&s(xyzzyaaay12(1,1))))
call wout('No. of occupied spin-down bands at each k point: '//trim(i2&
&s(xyzzyaaay12(1,2))))
endif
endif
call wout()
endif
if(xyzzyaabm1)then
if(xyzzyaabl1)then
call wout('No. of spin-up localized bands  : '//trim(i2s(xyzzyaaad1(1)&
&)))
call wout('No. of spin-down localized bands: '//trim(i2s(xyzzyaaad1(2)&
&)))
else
call wout('No. of localized bands: '//trim(i2s(xyzzyaaad1(1))))
endif
call wout()
endif
deallocate(xyzzyaaay12)
endif
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
if(periodicity<3.and.abs(dot_product(xyzzyaaay1(1:3,xyzzyaaab12),pa3(1&
&:3)))>xyzzyaabl12)call errstop('BWFDET_SETUP','k vectors inappropriat&
&e for periodicity.')
if(periodicity<2.and.abs(dot_product(xyzzyaaay1(1:3,xyzzyaaab12),pa2(1&
&:3)))>xyzzyaabl12)call errstop('BWFDET_SETUP','k vectors inappropriat&
&e for periodicity.')
if(periodicity<1.and.abs(dot_product(xyzzyaaay1(1:3,xyzzyaaab12),pa1(1&
&:3)))>xyzzyaabl12)call errstop('BWFDET_SETUP','k vectors inappropriat&
&e for periodicity.')
endif
enddo
endif
do xyzzyaaan12=1,nspin
if(xyzzyaaan12==2.and.xyzzyaabl1)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
do xyzzyaaao12=1,ndet
do xyzzyaaah12=1,wf_nd(xyzzyaaao12,xyzzyaaan12)
xyzzyaaak12=wf_d(1,xyzzyaaah12,xyzzyaaao12,xyzzyaaan12)
xyzzyaaai12=wf_d(2,xyzzyaaah12,xyzzyaaao12,xyzzyaaan12)
xyzzyaaal12=wf_d(3,xyzzyaaah12,xyzzyaaao12,xyzzyaaan12)
xyzzyaaaj12=wf_d(4,xyzzyaaah12,xyzzyaaao12,xyzzyaaan12)
if(xyzzyaaaj12<1.or.xyzzyaaaj12>xyzzyaaaf1.or.xyzzyaaai12<1.or.xyzzyaa&
&ai12>xyzzyaaaf1)then
if(am_master)call errstop('BWFDET_SETUP','k-point out of range in exci&
&tation specification.')
call qmc_barrier
endif
if(.not.xyzzyaabp1(xyzzyaaaj12).or..not.xyzzyaabp1(xyzzyaaai12))then
if(am_master)call errstop('BWFDET_SETUP', 'Specified excitation refers&
& to a k point that is not used because it is the (-) part of a (+/-)k&
& pair.  Use the corresponding +k in your excitation specification.')
call qmc_barrier
endif
if(xyzzyaaal12<1.or.xyzzyaaal12>xyzzyaaah1(xyzzyaaaj12,xyzzyaaad12).or&
&.xyzzyaaak12<1.or.xyzzyaaak12>xyzzyaaah1(xyzzyaaai12,xyzzyaaad12))the&
&n
call errstop('BWFDET_SETUP','band index out of range in excitation spe&
&cification.')
call qmc_barrier
endif
if(xyzzyaabq1(xyzzyaaaj12))then
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaal12,xyzzyaaaj12,xyzzyaaan12)>=1)then
if(am_master)call errstop('BWFDET_SETUP', 'Trying to promote an electr&
&on into a band that is already fully occupied.')
call qmc_barrier
endif
else
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaal12,xyzzyaaaj12,xyzzyaaan12)>=2)then
if(am_master)call errstop('BWFDET_SETUP','Trying to promote an electro&
&n into a band that is already fully occupied.')
call qmc_barrier
endif
endif
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaak12,xyzzyaaai12,xyzzyaaan12)<=0)then
if(am_master)call errstop('BWFDET_SETUP', 'Trying to promote an electr&
&on from an empty band.')
call qmc_barrier
endif
xyzzyaaai1(xyzzyaaao12,xyzzyaaal12,xyzzyaaaj12,xyzzyaaan12)=xyzzyaaai1&
&(xyzzyaaao12,xyzzyaaal12,xyzzyaaaj12,xyzzyaaan12)+1
xyzzyaaai1(xyzzyaaao12,xyzzyaaak12,xyzzyaaai12,xyzzyaaan12)=xyzzyaaai1&
&(xyzzyaaao12,xyzzyaaak12,xyzzyaaai12,xyzzyaaan12)-1
enddo
i=0
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
i=i+sum(xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12))
if(xyzzyaabq1(xyzzyaaab12).and.any(xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab1&
&2,xyzzyaaan12)>1))then
if(am_master)call errstop('BWFDET_SETUP','Problem with excitations (1)&
&.')
call qmc_barrier
endif
endif
enddo
if(i/=nuc_nele(xyzzyaaan12)-xyzzyaaad1(xyzzyaaad12))then
if(am_master)call errstop('BWFDET_SETUP','Problem with excitations (2)&
&.')
call qmc_barrier
endif
if(any(xyzzyaaai1(xyzzyaaao12,:,:,xyzzyaaan12)<0).or.any(xyzzyaaai1(xy&
&zzyaaao12,:,:,xyzzyaaan12)>2))then
if(am_master)call errstop('BWFDET_SETUP','Problem with excitations (3)&
&.')
call qmc_barrier
endif
enddo
enddo
deallocate(xyzzyaaaw12,xyzzyaabo12,xyzzyaaax12)
if(band_occup_only.and..not.xyzzyaabo1)deallocate(xyzzyaabw12)
if(.not.xyzzyaabk1.and..not.xyzzyaabr12)then
if(no_orb_phases>0)then
if(am_master.and..not.band_occup_only)then
call wout('User-defined phases')
call wout('-------------------')
call wout('Number   Det   Spin   Band    k point       Phase (radians)&
&')
endif
do i=1,no_orb_phases
if(am_master.and..not.band_occup_only)then
write(tmpr,'(1x,i5,1x,i5,5x,i2,1x,i6,2x,i9,7x,f10.6)')i,orb_phase_det(&
&i),orb_phase_spin(i),orb_phase_band(i),orb_phase_kpoint(i),orb_phase(&
&i)
call wout(tmpr)
endif
if(orb_phase_det(i)<1.or.orb_phase_det(i)>ndet)then
if(am_master)call errstop('BWFDET_SETUP','Error in determinant in phas&
&e spec.')
call qmc_barrier
endif
if(orb_phase_spin(i)<1.or.orb_phase_spin(i)>nspin)then
if(am_master)call errstop('BWFDET_SETUP','Error in spin in phase spec.&
&')
call qmc_barrier
endif
if(orb_phase_spin(i)==2.and.xyzzyaabl1)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
if(orb_phase_kpoint(i)<1.or.orb_phase_kpoint(i)>xyzzyaaaf1)then
if(am_master)call errstop('BWFDET_SETUP','Error in k point in phase sp&
&ec.')
call qmc_barrier
endif
if(.not.xyzzyaabp1(orb_phase_kpoint(i)))then
if(am_master)call errstop('BWFDET_SETUP','Have specified a phase for a&
& k point that is not used because it is part of a (+/-)k pair.')
call qmc_barrier
endif
if(xyzzyaabq1(orb_phase_kpoint(i)).and..not.band_occup_only)call errwa&
&rn('BWFDET_SETUP','A phase has been specified for a k point at the BZ&
& edge.  The phase factor should only affect the normalization of the &
&wave function in this case.')
if(orb_phase_band(i)<1.or.orb_phase_band(i)>xyzzyaaah1(orb_phase_kpoin&
&t(i),xyzzyaaad12))then
if(am_master)call errstop('BWFDET_SETUP','Error in band in phase spec.&
&')
call qmc_barrier
endif
if(xyzzyaaai1(orb_phase_det(i),orb_phase_band(i),orb_phase_kpoint(i),o&
&rb_phase_spin(i))==2.and..not.band_occup_only)call errwarn('BWFDET_SE&
&TUP','A phase has been specified for a band that is occupied at both &
&k and -k.  This phase will be ignored.')
if(xyzzyaaai1(orb_phase_det(i),orb_phase_band(i),orb_phase_kpoint(i),o&
&rb_phase_spin(i))==0.and..not.band_occup_only)call errwarn('BWFDET_SE&
&TUP','A phase has been specified for a band that is unoccupied.  This&
& phase will be ignored.')
enddo
if(am_master.and..not.band_occup_only)call wout()
deallocate(orb_phase_det,orb_phase_band,orb_phase_kpoint,orb_phase_spi&
&n,orb_phase)
endif
endif
if(band_occup_only)return
if(am_master.and.xyzzyaabv12)then
call wout('Detailed band-occupancy information')
call wout('-----------------------------------')
do xyzzyaaab12=1,xyzzyaaaf1
if(xyzzyaabp1(xyzzyaaab12))then
if(xyzzyaabq1(xyzzyaaab12))then
call wout('K POINT '//trim(i2s(xyzzyaaab12))//'  [not part of a (+/-)k&
& pair]')
else
call wout('K POINT '//trim(i2s(xyzzyaaab12))//'  [part of a (+/-)k pai&
&r]')
endif
if(isperiodic)then
write(tmpr,'("k(au):",3(1x,es23.15))')xyzzyaaay1(1:periodicity,xyzzyaa&
&ab12)
call wout(tmpr)
write(tmpr,'("kfrac:",3(1x,es23.15))')matmul(xyzzyaaay1(1:periodicity,&
&xyzzyaaab12),pbinv(1:periodicity,1:periodicity))
call wout(tmpr)
endif
do xyzzyaaan12=1,nspin
if(nuc_nele(xyzzyaaan12)==0)cycle
if(xyzzyaabl1.and.xyzzyaaan12==2)then
xyzzyaaad12=2
else
xyzzyaaad12=1
endif
do xyzzyaaao12=1,ndet
if(ndet>1)then
call wout('Spin '//trim(i2s(xyzzyaaan12))//', determinant: '//trim(i2s&
&(xyzzyaaao12))//':')
else
call wout('Spin '//trim(i2s(xyzzyaaan12))//':')
endif
if(xyzzyaaad1(xyzzyaaad12)>0)call wout('Number of localized orbitals: &
&'//trim(i2s(xyzzyaaad1(xyzzyaaad12))))
if(any(xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)==2))then
char80=' Bands occupied at both k and -k:'
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)==2)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam12))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam12))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Highest occupied band energy: ',maxval(xyzzyaaax1(:,xyzzyaa&
&ab12,xyzzyaaad12),xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)==&
&2))
endif
if(any(xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)==1))then
if(xyzzyaabq1(xyzzyaaab12))then
char80=' Bands occupied at k:'
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)==1)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam12))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam12))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Highest occupied band energy: ',maxval(xyzzyaaax1(:,xyzzyaa&
&ab12,xyzzyaaad12),xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)==&
&1))
else
call wout('Single bands formed from the bands at +k and -k:')
call wout('  Band      Phase (radians)')
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)==1)then
if(allocated(xyzzyaaaw1))then
write(tmpr,'(i6,4x,f16.10)')xyzzyaaam12,xyzzyaaaw1(xyzzyaaao12,xyzzyaa&
&am12,xyzzyaaab12,xyzzyaaan12)
else
write(tmpr,'(i6,4x,a)')xyzzyaaam12,'Unknown'
endif
call wout(tmpr)
endif
enddo
endif
endif
if(any(xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)==0))then
if(xyzzyaabq1(xyzzyaaab12))then
char80=' Bands unoccupied at k:'
else
char80=' Bands unoccupied at both k and -k:'
endif
do xyzzyaaam12=1,xyzzyaaah1(xyzzyaaab12,xyzzyaaad12)
if(xyzzyaaai1(xyzzyaaao12,xyzzyaaam12,xyzzyaaab12,xyzzyaaan12)==0)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam12))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam12))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Lowest unoccupied band energy: ',minval(xyzzyaaax1(:,xyzzya&
&aab12,xyzzyaaad12),xyzzyaaai1(xyzzyaaao12,:,xyzzyaaab12,xyzzyaaan12)=&
&=0))
endif
enddo
enddo
else
write(tmpr,'(a,3(f10.6,a))')'K POINT '//trim(i2s(xyzzyaaab12))//'  [('&
&,xyzzyaaay1(1,xyzzyaaab12),',',xyzzyaaay1(2,xyzzyaaab12),',',xyzzyaaa&
&y1(3,xyzzyaaab12),')]:  NOT USED.'
call wout(tmpr)
endif
call wout()
enddo
endif
deallocate(xyzzyaabw12)
endif
xyzzyaaaj1=0.d0
do i=2,nemax
xyzzyaaaj1=xyzzyaaaj1+log(dble(i))
enddo
xyzzyaaaj1=orb_norm*exp(-xyzzyaaaj1/dble(2*nemax))
xyzzyaacv1=0
xyzzyaacr1=0
xyzzyaacs1=0
xyzzyaact1=0
xyzzyaacu1=0
if(xyzzyaabk1)then
if(xyzzyaabm1)then
xyzzyaacv1=xyzzyaacv1+1
xyzzyaacr1=xyzzyaacv1
endif
if(xyzzyaabn1)then
xyzzyaacv1=xyzzyaacv1+1
xyzzyaacs1=xyzzyaacv1
endif
allocate(xyzzyaacp1(xyzzyaacv1,nspin),xyzzyaacq1(2,xyzzyaaci1,xyzzyaac&
&v1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbrange (gamma)')
if(xyzzyaabm1)then
allocate(xyzzyaacn1(xyzzyaaci1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbmask (loc)')
call get_numerical_orbmask(xyzzyaaci1,nemax,nspin,ndet,xyzzyaacm1,xyzz&
&yaacn1)
if(xyzzyaaci1>xyzzyaacj1)xyzzyaacn1(xyzzyaacj1+1:,1:nspin)=.false.
call get_numerical_orbrange(xyzzyaaci1,nspin,xyzzyaacn1,xyzzyaacp1(xyz&
&zyaacr1,:),xyzzyaacq1(:,:,xyzzyaacr1,:))
deallocate(xyzzyaacn1)
endif
if(xyzzyaabn1)then
allocate(xyzzyaacn1(xyzzyaaci1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbmask (ext)')
call get_numerical_orbmask(xyzzyaaci1,nemax,nspin,ndet,xyzzyaacm1,xyzz&
&yaacn1)
if(xyzzyaacj1>0)xyzzyaacn1(1:xyzzyaacj1,1:nspin)=.false.
call get_numerical_orbrange(xyzzyaaci1,nspin,xyzzyaacn1,xyzzyaacp1(xyz&
&zyaacs1,:),xyzzyaacq1(:,:,xyzzyaacs1,:))
deallocate(xyzzyaacn1)
endif
else
if(xyzzyaabo1)then
xyzzyaacv1=xyzzyaacv1+1
xyzzyaact1=xyzzyaacv1
allocate(xyzzyaacp1(xyzzyaacv1,nspin),xyzzyaacq1(2,xyzzyaaci1,xyzzyaac&
&v1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbrange (non-gamma, c&
&mplx)')
allocate(xyzzyaacn1(xyzzyaaci1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbmask (non-gamma, cm&
&plx)')
call get_numerical_orbmask(xyzzyaaci1,nemax,nspin,ndet,xyzzyaacm1,xyzz&
&yaacn1)
call get_numerical_orbrange(xyzzyaaci1,nspin,xyzzyaacn1,xyzzyaacp1(xyz&
&zyaact1,:),xyzzyaacq1(:,:,xyzzyaact1,:))
deallocate(xyzzyaacn1)
else
xyzzyaacv1=xyzzyaacv1+1
xyzzyaacu1=xyzzyaacv1
xyzzyaacv1=xyzzyaacv1+1
xyzzyaact1=xyzzyaacv1
allocate(xyzzyaacp1(xyzzyaacv1,nspin),xyzzyaacq1(2,xyzzyaaci1,xyzzyaac&
&v1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbrange (non-gamma, r&
&eal)')
allocate(xyzzyaacn1(xyzzyaaci1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbmask (non-gamma, re&
&al; real)')
call get_numerical_orbmask(xyzzyaaci1,nemax,nspin,ndet,xyzzyaacm1,xyzz&
&yaacn1)
call get_numerical_orbrange(xyzzyaaci1,nspin,xyzzyaacn1,xyzzyaacp1(xyz&
&zyaacu1,:),xyzzyaacq1(:,:,xyzzyaacu1,:))
deallocate(xyzzyaacn1)
allocate(xyzzyaaaz12(nemax,nspin,ndet),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','temp_orbmap')
xyzzyaaaz12=0
do xyzzyaaao12=1,ndet
do xyzzyaaan12=1,nspin
do xyzzyaaaq12=1,nemax
if(xyzzyaacm1(xyzzyaaaq12,xyzzyaaan12,xyzzyaaao12)>0)xyzzyaaaz12(xyzzy&
&aaaq12,xyzzyaaan12,xyzzyaaao12)=xyzzyaacw1(xyzzyaacm1(xyzzyaaaq12,xyz&
&zyaaan12,xyzzyaaao12))
enddo
enddo
enddo
allocate(xyzzyaacn1(xyzzyaaci1,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbmask (non-gamma, re&
&al; real)')
call get_numerical_orbmask(xyzzyaaci1,nemax,nspin,ndet,xyzzyaaaz12,xyz&
&zyaacn1)
call get_numerical_orbrange(xyzzyaaci1,nspin,xyzzyaacn1,xyzzyaacp1(xyz&
&zyaact1,:),xyzzyaacq1(:,:,xyzzyaact1,:))
deallocate(xyzzyaaaz12,xyzzyaacn1)
endif
endif
if(xyzzyaabk1)then
localized_orbitals=xyzzyaabm1
if(xyzzyaabm1)then
if(any(xyzzyaabh1==2))then
xyzzyaabh12(1)=sqrt(dot_product(painv(:,1),painv(:,1)))
xyzzyaabh12(2)=sqrt(dot_product(painv(:,2),painv(:,2)))
xyzzyaabh12(3)=sqrt(dot_product(painv(:,3),painv(:,3)))
allocate(xyzzyaabz1(3,xyzzyaacj1),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','half_side_len')
do xyzzyaaaf12=1,xyzzyaacj1
xyzzyaabz1(1:3,xyzzyaaaf12)=xyzzyaabh12(1:3)*xyzzyaabc1(xyzzyaaaf12)
enddo
endif
endif
endif
xyzzyaaak1(1:3)=painv(1,:)**2+painv(2,:)**2+painv(3,:)**2
xyzzyaaak1(4)=2.d0*(painv(1,1)*painv(1,2)+painv(2,1)*painv(2,2)+painv(&
&3,1)*painv(3,2))
xyzzyaaak1(5)=2.d0*(painv(1,3)*painv(1,1)+painv(2,3)*painv(2,1)+painv(&
&3,3)*painv(3,1))
xyzzyaaak1(6)=2.d0*(painv(1,2)*painv(1,3)+painv(2,2)*painv(2,3)+painv(&
&3,2)*painv(3,3))
xyzzyaaal1(1,1:3)=painv(:,1)*painv(:,1)
xyzzyaaal1(2,1:3)=painv(:,2)*painv(:,2)
xyzzyaaal1(3,1:3)=painv(:,3)*painv(:,3)
xyzzyaaal1(4,1:3)=2.d0*painv(:,1)*painv(:,2)
xyzzyaaal1(5,1:3)=2.d0*painv(:,1)*painv(:,3)
xyzzyaaal1(6,1:3)=2.d0*painv(:,2)*painv(:,3)
xyzzyaaal1(1,4)=painv(1,1)*painv(2,1)
xyzzyaaal1(2,4)=painv(1,2)*painv(2,2)
xyzzyaaal1(3,4)=painv(1,3)*painv(2,3)
xyzzyaaal1(4,4)=painv(1,1)*painv(2,2)+painv(1,2)*painv(2,1)
xyzzyaaal1(5,4)=painv(1,1)*painv(2,3)+painv(1,3)*painv(2,1)
xyzzyaaal1(6,4)=painv(1,2)*painv(2,3)+painv(1,3)*painv(2,2)
xyzzyaaal1(1,5)=painv(1,1)*painv(3,1)
xyzzyaaal1(2,5)=painv(1,2)*painv(3,2)
xyzzyaaal1(3,5)=painv(1,3)*painv(3,3)
xyzzyaaal1(4,5)=painv(1,1)*painv(3,2)+painv(1,2)*painv(3,1)
xyzzyaaal1(5,5)=painv(1,1)*painv(3,3)+painv(1,3)*painv(3,1)
xyzzyaaal1(6,5)=painv(1,2)*painv(3,3)+painv(1,3)*painv(3,2)
xyzzyaaal1(1,6)=painv(2,1)*painv(3,1)
xyzzyaaal1(2,6)=painv(2,2)*painv(3,2)
xyzzyaaal1(3,6)=painv(2,3)*painv(3,3)
xyzzyaaal1(4,6)=painv(2,1)*painv(3,2)+painv(2,2)*painv(3,1)
xyzzyaaal1(5,6)=painv(2,1)*painv(3,3)+painv(2,3)*painv(3,1)
xyzzyaaal1(6,6)=painv(2,2)*painv(3,3)+painv(2,3)*painv(3,2)
if(xyzzyaabk1)then
xyzzyaaas1=bw_norbc
allocate(xyzzyaaau1(xyzzyaaas1,64),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','mavc')
xyzzyaaau1=0.d0
else
xyzzyaaat1=bw_norbc
allocate(xyzzyaaav1(xyzzyaaat1,64),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','mcavc')
xyzzyaaav1=czero
endif
allocate(xyzzyaadh1(xyzzyaadf1,xyzzyaaci1),xyzzyaadi1(xyzzyaadg1,xyzzy&
&aaci1),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbdesc_*')
xyzzyaaat12=0
do xyzzyaaaq12=1,xyzzyaack1
if(xyzzyaabo1)then
xyzzyaaau12=xyzzyaaaq12
else
xyzzyaaau12=xyzzyaacw1(xyzzyaaaq12)
endif
xyzzyaaab12=xyzzyaacy1(xyzzyaaau12)
xyzzyaaam12=xyzzyaacx1(xyzzyaaau12)
xyzzyaaad12=xyzzyaacz1(xyzzyaaau12)
xyzzyaadh1(1,xyzzyaaaq12)=xyzzyaaab12
xyzzyaabt12=.false.
do xyzzyaaav12=1,xyzzyaaaq12-1
if(xyzzyaaax1(xyzzyaaam12,xyzzyaaab12,xyzzyaaad12)==xyzzyaadi1(1,xyzzy&
&aaav12))then
xyzzyaabt12=.true.
exit
endif
enddo
if(.not.xyzzyaabt12)then
xyzzyaaat12=xyzzyaaat12+1
xyzzyaaas12=xyzzyaaat12
xyzzyaadh1(2,xyzzyaaaq12)=xyzzyaaas12
xyzzyaadi1(1,xyzzyaaaq12)=xyzzyaaax1(xyzzyaaam12,xyzzyaaab12,xyzzyaaad&
&12)
else
xyzzyaadh1(2,xyzzyaaaq12)=xyzzyaadh1(2,xyzzyaaav12)
xyzzyaadi1(1,xyzzyaaaq12)=xyzzyaadi1(1,xyzzyaaav12)
endif
enddo
deallocate(xyzzyaaax1)
allocate(xyzzyaacn1(xyzzyaaci1,2),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_orbmask')
call get_numerical_orbmask(xyzzyaaci1,nemax,nspin,ndet,xyzzyaacm1,xyzz&
&yaacn1)
if(xyzzyaabn1)then
allocate(xyzzyaaco1(bw_norbc,nspin),stat=xyzzyaaaa12)
call check_alloc(xyzzyaaaa12,'BWFDET_SETUP','bw_corbmask')
xyzzyaaco1=.false.
do xyzzyaaan12=1,nspin
do xyzzyaaaq12=1,xyzzyaack1
if(.not.xyzzyaacn1(xyzzyaaaq12,xyzzyaaan12))cycle
if(xyzzyaabo1)then
xyzzyaaco1(xyzzyaaaq12,xyzzyaaan12)=.true.
else
xyzzyaaco1(xyzzyaacw1(xyzzyaaaq12),xyzzyaaan12)=.true.
endif
enddo
enddo
endif
if(use_gpcc)call xyzzyaadu1
deallocate(xyzzyaacn1)
end subroutine bwfdet_setup
subroutine xyzzyaadt1(iorb,norb,idet,ik,ib,spin,phase,cosphase,sinphas&
&e)
implicit none
integer,intent(in) :: iorb,norb,idet,ik,ib,spin
real(dp),intent(out) :: phase,cosphase,sinphase
integer xyzzyaaaa13
real(dp) xyzzyaaab13,xyzzyaaac13
phase=0.d0
cosphase=1.d0
sinphase=0.d0
call xyzzyaaej1(iorb,norb,ik,xyzzyaaab13,xyzzyaaac13)
if(xyzzyaaab13<xyzzyaaac13)then
phase=pi_over_two
cosphase=0.d0
sinphase=1.d0
endif
if(no_orb_phases>0)then
do xyzzyaaaa13=1,no_orb_phases
if(orb_phase_det(xyzzyaaaa13)==idet.and.orb_phase_band(xyzzyaaaa13)==i&
&b.and.orb_phase_kpoint(xyzzyaaaa13)==ik.and.orb_phase_spin(xyzzyaaaa1&
&3)==spin)then
phase=orb_phase(xyzzyaaaa13)
cosphase=cos(orb_phase(xyzzyaaaa13))
sinphase=sin(orb_phase(xyzzyaaaa13))
endif
enddo
endif
end subroutine xyzzyaadt1
subroutine deshalloc_bwfdet_shm
implicit none
if(associated(xyzzyaaao1))call deshallocate(xyzzyaaao1)
if(associated(xyzzyaaap1))call deshallocate(xyzzyaaap1)
if(associated(xyzzyaaaq1))call deshallocate(xyzzyaaaq1)
if(associated(xyzzyaaar1))call deshallocate(xyzzyaaar1)
end subroutine deshalloc_bwfdet_shm
subroutine xyzzyaadu1
use slaarnabi,only : use_gpcc,initialize_gpcc,setup_gpcc,spherical_av_&
&real,spherical_av_cmplx,nradgrid,radgrid,nsphgrid,sphgrid,ae_index,na&
&eions_prim
implicit none
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15,xyzzyaaad15,xyzzyaaae15,xy&
&zzyaaaf15,xyzzyaaag15,ik,xyzzyaaah15,xyzzyaaai15,xyzzyaaaj15
real(dp) xyzzyaaak15(3)
real(dp),allocatable :: xyzzyaaal15(:,:),xyzzyaaam15(:,:,:),xyzzyaaan1&
&5(:,:),xyzzyaaao15(:,:),xyzzyaaap15(:,:)
complex(dp),allocatable :: xyzzyaaaq15(:,:)
logical xyzzyaaar15
if(xyzzyaabm1.and.am_master)call errstop('BWFDET_SETUP_GPCC','Cannot c&
&usp-correct localized orbitals at present. Ask us if you need to do t&
&his.')
if(.not.xyzzyaabo1)then
do ik=1,xyzzyaaaf1
if(xyzzyaabp1(ik).and..not.xyzzyaabq1(ik).and.any(xyzzyaaai1(:,:,ik,:)&
&==1).and.am_master)call errstop('BWFDET_SETUP_GPCC','Cannot cusp-corr&
&ect orbitals away from the BZ edge in which only one of k and -k is o&
&ccupied. Ask us if you need to do this, or use a complex wave functio&
&n.')
enddo
endif
if(xyzzyaabl1.or.nuc_nele(2)/=nuc_nele(1).or.any(xyzzyaaai1(:,:,:,2)/=&
&xyzzyaaai1(:,:,:,1)))then
xyzzyaaaa15=1
else
xyzzyaaaa15=0
endif
call initialize_gpcc(xyzzyaaaa15,bw_norbc,xyzzyaaco1,xyzzyaabk1)
allocate(xyzzyaaal15(xyzzyaaci1,real1_complex2),xyzzyaaam15(3,xyzzyaac&
&i1,real1_complex2),xyzzyaaan15(xyzzyaaci1,real1_complex2),stat=xyzzya&
&aab15)
call check_alloc(xyzzyaaab15,'BWFDET_SETUP_GPCC','valtemp')
if(xyzzyaabk1)then
allocate(xyzzyaaao15(nsphgrid,bw_norbc),stat=xyzzyaaab15)
xyzzyaaao15=0.d0
else
allocate(xyzzyaaaq15(nsphgrid,bw_norbc),stat=xyzzyaaab15)
xyzzyaaaq15=czero
endif
call check_alloc(xyzzyaaab15,'BWFDET_SETUP_GPCC','orb_sphgrid')
use_gpcc=.false.
do xyzzyaaac15=1,no_ssingles(xyzzyaaaa15)
do xyzzyaaai15=1,naeions_prim
xyzzyaaah15=ae_index(xyzzyaaai15)
do xyzzyaaad15=1,nradgrid
if(xyzzyaaad15>1)then
xyzzyaaaj15=nsphgrid
else
xyzzyaaaj15=1
endif
do xyzzyaaae15=1,xyzzyaaaj15
xyzzyaaak15=rion(1:3,xyzzyaaah15)+radgrid(xyzzyaaad15,xyzzyaaai15)*sph&
&grid(1:3,xyzzyaaae15)
call blip_orb_eval(xyzzyaaak15,xyzzyaaac15,xyzzyaaci1,xyzzyaacn1(1,xyz&
&zyaaac15),.true.,.false.,xyzzyaaal15,xyzzyaaam15,xyzzyaaan15)
xyzzyaaar15=.false.
do xyzzyaaaf15=1,xyzzyaaci1
if(.not.xyzzyaacn1(xyzzyaaaf15,xyzzyaaac15))cycle
if(xyzzyaabo1)then
xyzzyaaaq15(xyzzyaaae15,xyzzyaaaf15)=cmplx(xyzzyaaal15(xyzzyaaaf15,1),&
&xyzzyaaal15(xyzzyaaaf15,2),dp)
else
if(xyzzyaaar15)then
xyzzyaaar15=.false.
cycle
endif
xyzzyaaag15=xyzzyaacx1(xyzzyaacw1(xyzzyaaaf15))
ik=xyzzyaacy1(xyzzyaacw1(xyzzyaaaf15))
if(any(xyzzyaaai1(:,xyzzyaaag15,ik,xyzzyaaac15)==2))then
xyzzyaaaq15(xyzzyaaae15,xyzzyaacw1(xyzzyaaaf15))=cmplx(xyzzyaaal15(xyz&
&zyaaaf15,1),xyzzyaaal15(xyzzyaaaf15+1,1),dp)
xyzzyaaar15=.true.
elseif(any(xyzzyaaai1(:,xyzzyaaag15,ik,xyzzyaaac15)==1))then
if(xyzzyaabk1)then
xyzzyaaao15(xyzzyaaae15,xyzzyaacw1(xyzzyaaaf15))=xyzzyaaal15(xyzzyaaaf&
&15,1)
else
xyzzyaaaq15(xyzzyaaae15,xyzzyaacw1(xyzzyaaaf15))=cmplx(xyzzyaaal15(xyz&
&zyaaaf15,1)*dble(xyzzyaade1(xyzzyaaaf15)),xyzzyaaal15(xyzzyaaaf15,1)*&
&aimag(xyzzyaade1(xyzzyaaaf15)),dp)
endif
endif
endif
enddo
enddo
if(xyzzyaabk1)then
call spherical_av_real(xyzzyaaac15,xyzzyaaai15,xyzzyaaad15,xyzzyaaao15&
&)
else
call spherical_av_cmplx(xyzzyaaac15,xyzzyaaai15,xyzzyaaad15,xyzzyaaaq1&
&5)
endif
enddo
enddo
enddo
deallocate(xyzzyaaal15,xyzzyaaam15,xyzzyaaan15)
if(xyzzyaabk1)then
deallocate(xyzzyaaao15)
else
deallocate(xyzzyaaaq15)
endif
use_gpcc=.true.
if(xyzzyaabk1)then
call setup_gpcc
else
allocate(xyzzyaaap15(3,bw_norbc),stat=xyzzyaaab15)
call check_alloc(xyzzyaaab15,'BWFDET_SETUP_GPCC','kvec_gpcc')
xyzzyaaap15=0.d0
do xyzzyaaaf15=1,xyzzyaaci1
if(xyzzyaabo1)then
xyzzyaaap15(1:3,xyzzyaaaf15)=xyzzyaaay1(1:3,xyzzyaacy1(xyzzyaaaf15))
else
xyzzyaaap15(1:3,xyzzyaacw1(xyzzyaaaf15))=xyzzyaaay1(1:3,xyzzyaacy1(xyz&
&zyaacw1(xyzzyaaaf15)))
endif
enddo
call setup_gpcc(xyzzyaaap15)
deallocate(xyzzyaaap15)
endif
if(xyzzyaabk1)then
allocate(xyzzyaaca1(bw_norbc),xyzzyaacb1(3,bw_norbc),xyzzyaacc1(bw_nor&
&bc),stat=xyzzyaaab15)
call check_alloc(xyzzyaaab15,'BWFDET_SETUP_GPCC','cusp_val_r')
if(use_backflow)then
allocate(xyzzyaacd1(6,bw_norbc),stat=xyzzyaaab15)
call check_alloc(xyzzyaaab15,'BWFDET_SETUP_GPCC','cusp_sderivs_r')
endif
else
allocate(xyzzyaace1(bw_norbc),xyzzyaacf1(3,bw_norbc),xyzzyaacg1(bw_nor&
&bc),stat=xyzzyaaab15)
call check_alloc(xyzzyaaab15,'BWFDET_SETUP_GPCC','cusp_val_c')
if(use_backflow)then
allocate(xyzzyaach1(6,bw_norbc),stat=xyzzyaaab15)
call check_alloc(xyzzyaaab15,'BWFDET_SETUP_GPCC','cusp_sderivs_c')
endif
endif
end subroutine xyzzyaadu1
subroutine blip_orb_eval(rvec,jspin,norb,orbmask,val,fsd,orbval,orbgra&
&d,orblap,orbsderivs,orb_m,orb_rmap)
implicit none
integer,intent(in) :: jspin,norb
integer,intent(out),optional :: orb_m,orb_rmap(norb)
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd,orbmask(norb)
integer xyzzyaaaa16(3),xyzzyaaab16(3)
real(dp) xyzzyaaac16(3),xyzzyaaad16(64),xyzzyaaae16(64,3),xyzzyaaaf16(&
&64,6)
logical xyzzyaaag16,xyzzyaaah16
if(val)orbval=0.d0
if(fsd)then
orbgrad=0.d0
orblap=0.d0
if(present(orbsderivs))orbsderivs=0.d0
endif
if(present(orb_m))orb_m=0
xyzzyaaac16=matmul(rvec(1:3),painv)
call xyzzyaaek1(xyzzyaaac16,xyzzyaaaa16(1),xyzzyaaaa16(2),xyzzyaaaa16(&
&3),xyzzyaaag16)
if(.not.xyzzyaaag16)return
xyzzyaaab16(:)=xyzzyaaaa16(:)-1
where(xyzzyaaab16<0)xyzzyaaab16=xyzzyaaab16+xyzzyaaaa1
xyzzyaaah16=.false.
if(xyzzyaabk1)then
if(present(orbsderivs).and.present(orb_m))then
if(xyzzyaabm1)call xyzzyaadv1(rvec,norb,xyzzyaacv1,xyzzyaacp1(1,jspin)&
&,xyzzyaacq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzyaaab16,&
&xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap,orbsderivs,&
&orb_m,orb_rmap)
if(xyzzyaabn1)call xyzzyaadx1(rvec,jspin,norb,xyzzyaacv1,xyzzyaacp1(1,&
&jspin),xyzzyaacq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzya&
&aab16,xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap,orbsd&
&erivs,orb_m,orb_rmap)
elseif(present(orbsderivs).and..not.present(orb_m))then
if(xyzzyaabm1)call xyzzyaadv1(rvec,norb,xyzzyaacv1,xyzzyaacp1(1,jspin)&
&,xyzzyaacq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzyaaab16,&
&xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap,orbsderivs)
if(xyzzyaabn1)call xyzzyaadx1(rvec,jspin,norb,xyzzyaacv1,xyzzyaacp1(1,&
&jspin),xyzzyaacq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzya&
&aab16,xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap,orbsd&
&erivs)
elseif(.not.present(orbsderivs).and.present(orb_m))then
if(xyzzyaabm1)call xyzzyaadv1(rvec,norb,xyzzyaacv1,xyzzyaacp1(1,jspin)&
&,xyzzyaacq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzyaaab16,&
&xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap,orb_m=orb_m&
&,orb_rmap=orb_rmap)
if(xyzzyaabn1)call xyzzyaadx1(rvec,jspin,norb,xyzzyaacv1,xyzzyaacp1(1,&
&jspin),xyzzyaacq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzya&
&aab16,xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap,orb_m&
&=orb_m,orb_rmap=orb_rmap)
else
if(xyzzyaabm1)call xyzzyaadv1(rvec,norb,xyzzyaacv1,xyzzyaacp1(1,jspin)&
&,xyzzyaacq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzyaaab16,&
&xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap)
if(xyzzyaabn1)call xyzzyaadx1(rvec,jspin,norb,xyzzyaacv1,xyzzyaacp1(1,&
&jspin),xyzzyaacq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzya&
&aab16,xyzzyaaad16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap)
endif
else
if(xyzzyaabn1)then
if(present(orbsderivs))then
call xyzzyaadw1(rvec,jspin,norb,xyzzyaacv1,xyzzyaacp1(1,jspin),xyzzyaa&
&cq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzyaaab16,xyzzyaaa&
&d16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap,orbsderivs)
else
call xyzzyaadw1(rvec,jspin,norb,xyzzyaacv1,xyzzyaacp1(1,jspin),xyzzyaa&
&cq1(1,1,1,jspin),val,fsd,xyzzyaaah16,xyzzyaaac16,xyzzyaaab16,xyzzyaaa&
&d16,xyzzyaaae16,xyzzyaaaf16,orbval,orbgrad,orblap)
endif
endif
endif
end subroutine blip_orb_eval
subroutine xyzzyaadv1(rvec,norb,xyzzyaacv1,norbrange,orbrange,val,fsd,&
&blips_computed,r_frac,ir,f123,fg123,fl123,orbval,orbgrad,orblap,orbsd&
&erivs,orb_m,orb_rmap)
use slaarnabq, only : minimum_image
use slaarnabt, only : ddot
implicit none
integer,intent(in) :: norb,xyzzyaacv1,norbrange(xyzzyaacv1),orbrange(2&
&,xyzzyaaci1,xyzzyaacv1),ir(3)
integer,intent(out),optional :: orb_rmap(norb)
integer,intent(inout),optional :: orb_m
real(dp),intent(in) :: rvec(3),r_frac(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2),f123(64),fg123(64,3),fl123&
&(64,6)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd
logical,intent(inout) :: blips_computed
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17,xyzzyaaad17(3),xyzzyaaae17&
&,xyzzyaaaf17,xyzzyaaag17
real(dp) xyzzyaaah17(3),bcoef(64),xyzzyaaai17,xyzzyaaaj17(3),xyzzyaaak&
&17(6),xyzzyaaal17,xyzzyaaam17(6),xyzzyaaan17(3),xyzzyaaao17,xyzzyaaap&
&17,xyzzyaaaq17,xyzzyaaar17,xyzzyaaas17,xyzzyaaat17,xyzzyaaau17,xyzzya&
&aav17,xyzzyaaaw17,xyzzyaaax17,xyzzyaaay17,xyzzyaaaz17(1)
real(dp),pointer :: xyzzyaaba17(:,:,:)=>null()
real(sp),pointer :: xyzzyaabb17(:,:,:)=>null()
logical xyzzyaabc17,xyzzyaabd17
do xyzzyaaag17=1,norbrange(xyzzyaacr1)
xyzzyaaae17=orbrange(1,xyzzyaaag17,xyzzyaacr1)
xyzzyaaaf17=orbrange(2,xyzzyaaag17,xyzzyaacr1)
do xyzzyaaaa17=xyzzyaaae17,xyzzyaaaf17
xyzzyaabc17=.false.
xyzzyaabd17=.false.
if(xyzzyaabh1(xyzzyaaaa17)==1)then
xyzzyaaan17=rvec-xyzzyaabb1(:,xyzzyaaaa17)
if(periodicity>0)call minimum_image(3,1,xyzzyaaan17)
xyzzyaaao17=sum(xyzzyaaan17**2)
if(xyzzyaaao17<xyzzyaabf1(xyzzyaaaa17))then
xyzzyaabc17=.true.
if(bsmooth.and.xyzzyaaao17>xyzzyaabd1(xyzzyaaaa17))xyzzyaabd17=.true.
endif
else
xyzzyaaah17=matmul(rvec(1:3)-xyzzyaabb1(1:3,xyzzyaaaa17),painv)
xyzzyaaah17(1:periodicity)=modulo(xyzzyaaah17(1:periodicity)+0.5d0,1.d&
&0)-0.5d0
if(all(abs(xyzzyaaah17)<xyzzyaabz1(1:3,xyzzyaaaa17)))xyzzyaabc17=.true&
&.
endif
if(.not.xyzzyaabc17)cycle
call xyzzyaaeg1(val,fsd,blips_computed,r_frac,f123,fg123,fl123)
xyzzyaaad17(:)=modulo(ir(:)-xyzzyaabg1(:,xyzzyaaaa17),xyzzyaabi1(:,xyz&
&zyaaaa17))
if(single_precision_blips)then
xyzzyaabb17=>xyzzyaaam1(xyzzyaaaa17)%bcoefs(:,:,:)
else
xyzzyaaba17=>xyzzyaaan1(xyzzyaaaa17)%bcoefs(:,:,:)
endif
if(fsd)then
if(single_precision_blips)then
call xyzzyaady1(xyzzyaaad17,xyzzyaabi1(1,xyzzyaaaa17),1,1,1,xyzzyaabb1&
&7,bcoef)
else
call xyzzyaadz1(xyzzyaaad17,xyzzyaabi1(1,xyzzyaaaa17),1,1,1,xyzzyaaba1&
&7,bcoef)
endif
if(val.or.xyzzyaabd17)xyzzyaaai17=ddot(64,bcoef,1,f123(1),1)
do xyzzyaaab17=1,3
xyzzyaaaj17(xyzzyaaab17)=ddot(64,bcoef,1,fg123(1,xyzzyaaab17),1)
enddo
do xyzzyaaab17=1,6
xyzzyaaak17(xyzzyaaab17)=ddot(64,bcoef,1,fl123(1,xyzzyaaab17),1)
enddo
xyzzyaaaj17(1:3)=(/ddot(3,painv(1,1),3,xyzzyaaaj17(1),1),ddot(3,painv(&
&2,1),3,xyzzyaaaj17(1),1),ddot(3,painv(3,1),3,xyzzyaaaj17(1),1)/)
xyzzyaaal17=ddot(6,xyzzyaaak17(1),1,xyzzyaaak1(1),1)
if(present(orbsderivs))then
do xyzzyaaab17=1,6
xyzzyaaam17(xyzzyaaab17)=ddot(6,xyzzyaaak17(1),1,xyzzyaaal1(1,xyzzyaaa&
&b17),1)
enddo
endif
elseif(val.or.xyzzyaabd17)then
if(single_precision_blips)then
call xyzzyaaec1(xyzzyaaad17,xyzzyaabi1(1,xyzzyaaaa17),1,1,1,xyzzyaabb1&
&7,f123,xyzzyaaaz17)
else
call xyzzyaaed1(xyzzyaaad17,xyzzyaabi1(1,xyzzyaaaa17),1,1,1,xyzzyaaba1&
&7,f123,xyzzyaaaz17)
endif
xyzzyaaai17=xyzzyaaaz17(1)
endif
if(xyzzyaabd17)then
xyzzyaaap17=sqrt(xyzzyaaao17)
xyzzyaaaq17=1.d0/(xyzzyaabe1(xyzzyaaaa17)-xyzzyaabc1(xyzzyaaaa17))
xyzzyaaar17=(xyzzyaaap17-xyzzyaabc1(xyzzyaaaa17))*xyzzyaaaq17
xyzzyaaay17=1.d0/xyzzyaaap17
xyzzyaaas17=1.d0
xyzzyaaat17=0.d0
xyzzyaaaw17=0.d0
if(fsd)then
xyzzyaaau17=3.d0*xyzzyaaar17**2
xyzzyaaav17=xyzzyaaar17**3
do xyzzyaaac17=4,9
xyzzyaaax17=dble(xyzzyaaac17)*xyzzyaaau17
xyzzyaaau17=dble(xyzzyaaac17)*xyzzyaaav17
xyzzyaaav17=xyzzyaaav17*xyzzyaaar17
xyzzyaaas17=xyzzyaaas17+xyzzyaaby1(xyzzyaaac17)*xyzzyaaav17
xyzzyaaat17=xyzzyaaat17+xyzzyaaby1(xyzzyaaac17)*xyzzyaaau17
xyzzyaaaw17=xyzzyaaaw17+xyzzyaaby1(xyzzyaaac17)*xyzzyaaax17
enddo
xyzzyaaat17=xyzzyaaat17*xyzzyaaaq17
xyzzyaaaw17=xyzzyaaaw17*xyzzyaaaq17**2
xyzzyaaal17=xyzzyaaal17*xyzzyaaas17+2.d0*xyzzyaaat17*xyzzyaaay17*(dot_&
&product(xyzzyaaaj17,xyzzyaaan17)+xyzzyaaai17)+xyzzyaaaw17*xyzzyaaai17
xyzzyaaaj17=xyzzyaaaj17*xyzzyaaas17+xyzzyaaai17*xyzzyaaat17*xyzzyaaan1&
&7*xyzzyaaay17
xyzzyaaai17=xyzzyaaai17*xyzzyaaas17
else
xyzzyaaav17=xyzzyaaar17**3
do xyzzyaaac17=4,9
xyzzyaaav17=xyzzyaaav17*xyzzyaaar17
xyzzyaaas17=xyzzyaaas17+xyzzyaaby1(xyzzyaaac17)*xyzzyaaav17
enddo
xyzzyaaai17=xyzzyaaai17*xyzzyaaas17
endif
if(present(orbsderivs))call errstop('BLIP3DLOC','Second derivatives no&
&t implemented for shell region.')
endif
if(val)orbval(xyzzyaaaa17,1)=xyzzyaaai17
if(fsd)then
orbgrad(:,xyzzyaaaa17,1)=xyzzyaaaj17
orblap(xyzzyaaaa17,1)=xyzzyaaal17
if(present(orbsderivs))orbsderivs(:,xyzzyaaaa17,1)=xyzzyaaam17(:)
endif
if(present(orb_m))then
orb_m=orb_m+1
orb_rmap(orb_m)=xyzzyaaaa17
endif
enddo
enddo
end subroutine xyzzyaadv1
subroutine xyzzyaadw1(rvec,jspin,norb,xyzzyaacv1,norbrange,orbrange,va&
&l,fsd,blips_computed,r_frac,ir,f123,fg123,fl123,orbval,orbgrad,orblap&
&,orbsderivs)
use slaarnabi,      only : use_gpcc,cusp_wfdet_cmplx
use slaarnabt, only : zdot_rc
!$  use openmp_base, only : get_omp_index_range
implicit none
integer,intent(in) :: jspin,norb,xyzzyaacv1,norbrange(xyzzyaacv1),orbr&
&ange(2,xyzzyaaci1,xyzzyaacv1),ir(3)
real(dp),intent(in) :: rvec(3),r_frac(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2),f123(64),fg123(64,3),fl123&
&(64,6)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd
logical,intent(inout) :: blips_computed
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xy&
&zzyaaaf18,xyzzyaaag18,xyzzyaaah18,xyzzyaaai18
real(dp) xyzzyaaaj18,xyzzyaaak18,xyzzyaaal18
complex(dp) xyzzyaaam18(bw_norbc),xyzzyaaan18(3,bw_norbc),xyzzyaaao18(&
&6,bw_norbc),xyzzyaaap18(bw_norbc),xyzzyaaaq18(bw_norbc),xyzzyaaar18(3&
&,bw_norbc),xyzzyaaas18(6,bw_norbc),xyzzyaaat18,xyzzyaaau18(3),xyzzyaa&
&av18,xyzzyaaaw18,xyzzyaaax18(3),xyzzyaaay18(6),xyzzyaaaz18(6),xyzzyaa&
&ba18(xyzzyaaaf1)
call xyzzyaaeg1(.true.,fsd,blips_computed,r_frac,f123,fg123,fl123)
!$omp parallel default(none) shared(rvec,jspin,val,fsd,ir,f123,fg123,f&
!$omp &l123,orbval,orbgrad,orblap,orbsderivs,xyzzyaaaa1,bw_norbc,norbr&
!$omp &ange,orbrange,xyzzyaact1,xyzzyaacu1,xyzzyaacw1,xyzzyaaaq1,xyzzy&
!$omp &aaar1,xyzzyaaav1,xyzzyaaat1,xyzzyaaam18,xyzzyaaan18,xyzzyaaao18&
!$omp &,xyzzyaace1,xyzzyaacf1,xyzzyaacg1,xyzzyaach1,xyzzyaaap18,xyzzya&
!$omp &aar18,xyzzyaaaq18,xyzzyaaas18,single_precision_blips,xyzzyaaay1&
!$omp &,xyzzyaaaz1,xyzzyaaba1,xyzzyaaaf1,xyzzyaabp1,xyzzyaaba18,xyzzya&
!$omp &aaj1,use_gpcc,painv,xyzzyaaak1,xyzzyaaal1,xyzzyaabo1,xyzzyaadd1&
!$omp &,xyzzyaade1,xyzzyaacy1)private(xyzzyaaaa18,xyzzyaaac18,xyzzyaaa&
!$omp &d18,xyzzyaaaf18,xyzzyaaag18,xyzzyaaab18,xyzzyaaah18,xyzzyaaai18&
!$omp &,xyzzyaaae18,xyzzyaaaj18,xyzzyaaav18,xyzzyaaaw18,xyzzyaaax18,xy&
!$omp &zzyaaat18,xyzzyaaau18,xyzzyaaay18,xyzzyaaaz18,xyzzyaaak18,xyzzy&
!$omp &aaal18)
if(fsd)then
if(single_precision_blips)then
do xyzzyaaae18=1,norbrange(xyzzyaact1)
xyzzyaaaf18=orbrange(1,xyzzyaaae18,xyzzyaact1)
xyzzyaaag18=orbrange(2,xyzzyaaae18,xyzzyaact1)
!$  call get_omp_index_range(xyzzyaaaf18,xyzzyaaag18)
!$  if(xyzzyaaag18<xyzzyaaaf18)cycle
call xyzzyaaea1(ir,xyzzyaaaa1,xyzzyaaaf18,xyzzyaaag18,bw_norbc,xyzzyaa&
&aq1,xyzzyaaav1)
enddo
else
do xyzzyaaae18=1,norbrange(xyzzyaact1)
xyzzyaaaf18=orbrange(1,xyzzyaaae18,xyzzyaact1)
xyzzyaaag18=orbrange(2,xyzzyaaae18,xyzzyaact1)
!$  call get_omp_index_range(xyzzyaaaf18,xyzzyaaag18)
!$  if(xyzzyaaag18<xyzzyaaaf18)cycle
call xyzzyaaeb1(ir,xyzzyaaaa1,xyzzyaaaf18,xyzzyaaag18,bw_norbc,xyzzyaa&
&ar1,xyzzyaaav1)
enddo
endif
do xyzzyaaae18=1,norbrange(xyzzyaact1)
xyzzyaaaf18=orbrange(1,xyzzyaaae18,xyzzyaact1)
xyzzyaaag18=orbrange(2,xyzzyaaae18,xyzzyaact1)
!$ call get_omp_index_range(xyzzyaaaf18,xyzzyaaag18)
!$ if(xyzzyaaag18<xyzzyaaaf18)cycle
do xyzzyaaad18=xyzzyaaaf18,xyzzyaaag18
xyzzyaaam18(xyzzyaaad18)=zdot_rc(64,f123(1),1,xyzzyaaav1(xyzzyaaad18,1&
&),xyzzyaaat1)
do xyzzyaaaa18=1,3
xyzzyaaan18(xyzzyaaaa18,xyzzyaaad18)=zdot_rc(64,fg123(1,xyzzyaaaa18),1&
&,xyzzyaaav1(xyzzyaaad18,1),xyzzyaaat1)
enddo
do xyzzyaaaa18=1,6
xyzzyaaao18(xyzzyaaaa18,xyzzyaaad18)=zdot_rc(64,fl123(1,xyzzyaaaa18),1&
&,xyzzyaaav1(xyzzyaaad18,1),xyzzyaaat1)
enddo
enddo
enddo
elseif(val)then
if(single_precision_blips)then
do xyzzyaaae18=1,norbrange(xyzzyaact1)
xyzzyaaaf18=orbrange(1,xyzzyaaae18,xyzzyaact1)
xyzzyaaag18=orbrange(2,xyzzyaaae18,xyzzyaact1)
!$  call get_omp_index_range(xyzzyaaaf18,xyzzyaaag18)
!$  if(xyzzyaaag18<xyzzyaaaf18)cycle
call xyzzyaaee1(ir,xyzzyaaaa1,xyzzyaaaf18,xyzzyaaag18,bw_norbc,xyzzyaa&
&aq1,f123,xyzzyaaam18(xyzzyaaaf18))
enddo
else
do xyzzyaaae18=1,norbrange(xyzzyaact1)
xyzzyaaaf18=orbrange(1,xyzzyaaae18,xyzzyaact1)
xyzzyaaag18=orbrange(2,xyzzyaaae18,xyzzyaact1)
!$  call get_omp_index_range(xyzzyaaaf18,xyzzyaaag18)
!$  if(xyzzyaaag18<xyzzyaaaf18)cycle
call xyzzyaaef1(ir,xyzzyaaaa1,xyzzyaaaf18,xyzzyaaag18,bw_norbc,xyzzyaa&
&ar1,f123,xyzzyaaam18(xyzzyaaaf18))
enddo
endif
endif
if(use_gpcc)then
if(present(orbsderivs))then
call cusp_wfdet_cmplx(rvec,jspin,.true.,fsd,xyzzyaace1,xyzzyaacf1,xyzz&
&yaacg1,xyzzyaach1)
else
call cusp_wfdet_cmplx(rvec,jspin,.true.,fsd,xyzzyaace1,xyzzyaacf1,xyzz&
&yaacg1)
endif
!$omp  barrier ! wait until the cusp corrections are done
endif
xyzzyaaah18=1
xyzzyaaai18=xyzzyaaaf1
!$ call get_omp_index_range(xyzzyaaah18,xyzzyaaai18)
do xyzzyaaab18=xyzzyaaah18,xyzzyaaai18
if(xyzzyaabp1(xyzzyaaab18))then
xyzzyaaaj18=dot_product(rvec(1:3),xyzzyaaay1(:,xyzzyaaab18))
xyzzyaaba18(xyzzyaaab18)=cmplx(xyzzyaaaj1*cos(xyzzyaaaj18),xyzzyaaaj1*&
&sin(xyzzyaaaj18),dp)
endif
enddo
!$omp  barrier ! wait until the exponentials are available
do xyzzyaaae18=1,norbrange(xyzzyaact1)
xyzzyaaaf18=orbrange(1,xyzzyaaae18,xyzzyaact1)
xyzzyaaag18=orbrange(2,xyzzyaaae18,xyzzyaact1)
!$ call get_omp_index_range(xyzzyaaaf18,xyzzyaaag18)
!$ if(xyzzyaaag18<xyzzyaaaf18)cycle
do xyzzyaaad18=xyzzyaaaf18,xyzzyaaag18
xyzzyaaab18=xyzzyaacy1(xyzzyaaad18)
xyzzyaaav18=xyzzyaaba18(xyzzyaaab18)
xyzzyaaap18(xyzzyaaad18)=xyzzyaaav18*xyzzyaaam18(xyzzyaaad18)
if(use_gpcc)xyzzyaaap18(xyzzyaaad18)=xyzzyaaap18(xyzzyaaad18)+xyzzyaac&
&e1(xyzzyaaad18)
if(fsd)then
xyzzyaaaw18=-xyzzyaaaz1(xyzzyaaab18)*xyzzyaaav18
xyzzyaaax18=xyzzyaaay1(:,xyzzyaaab18)*zi*xyzzyaaav18
do xyzzyaaaa18=1,6
xyzzyaaat18=zdot_rc(6,xyzzyaaak1(1),1,xyzzyaaao18(1,xyzzyaaad18),1)
enddo
xyzzyaaau18(1:3)=(/zdot_rc(3,painv(1,1),3,xyzzyaaan18(1,xyzzyaaad18),1&
&),zdot_rc(3,painv(2,1),3,xyzzyaaan18(1,xyzzyaaad18),1),zdot_rc(3,pain&
&v(3,1),3,xyzzyaaan18(1,xyzzyaaad18),1)/)
xyzzyaaaq18(xyzzyaaad18)=xyzzyaaaw18*xyzzyaaam18(xyzzyaaad18)+2.d0*(xy&
&zzyaaau18(1)*xyzzyaaax18(1)+xyzzyaaau18(2)*xyzzyaaax18(2)+xyzzyaaau18&
&(3)*xyzzyaaax18(3))+xyzzyaaav18*xyzzyaaat18
xyzzyaaar18(:,xyzzyaaad18)=xyzzyaaax18*xyzzyaaam18(xyzzyaaad18)+xyzzya&
&aav18*xyzzyaaau18
if(use_gpcc)then
xyzzyaaaq18(xyzzyaaad18)=xyzzyaaaq18(xyzzyaaad18)+xyzzyaacg1(xyzzyaaad&
&18)
xyzzyaaar18(:,xyzzyaaad18)=xyzzyaaar18(:,xyzzyaaad18)+xyzzyaacf1(1:3,x&
&yzzyaaad18)
endif
if(present(orbsderivs))then
xyzzyaaay18=-xyzzyaaba1(:,xyzzyaaab18)*xyzzyaaap18(xyzzyaaad18)
xyzzyaaaz18(1:3)=2.d0*xyzzyaaax18*xyzzyaaau18
xyzzyaaaz18(4)=xyzzyaaax18(1)*xyzzyaaau18(2)+xyzzyaaax18(2)*xyzzyaaau1&
&8(1)
xyzzyaaaz18(5)=xyzzyaaax18(1)*xyzzyaaau18(3)+xyzzyaaax18(3)*xyzzyaaau1&
&8(1)
xyzzyaaaz18(6)=xyzzyaaax18(2)*xyzzyaaau18(3)+xyzzyaaax18(3)*xyzzyaaau1&
&8(2)
xyzzyaaas18(:,xyzzyaaad18)=xyzzyaaay18(:)+xyzzyaaaz18(:)+xyzzyaaav18*(&
&xyzzyaaal1(1,:)*xyzzyaaao18(1,xyzzyaaad18)+xyzzyaaal1(2,:)*xyzzyaaao1&
&8(2,xyzzyaaad18)+xyzzyaaal1(3,:)*xyzzyaaao18(3,xyzzyaaad18)+xyzzyaaal&
&1(4,:)*xyzzyaaao18(4,xyzzyaaad18)+xyzzyaaal1(5,:)*xyzzyaaao18(5,xyzzy&
&aaad18)+xyzzyaaal1(6,:)*xyzzyaaao18(6,xyzzyaaad18))
if(use_gpcc)xyzzyaaas18(:,xyzzyaaad18)=xyzzyaaas18(:,xyzzyaaad18)+xyzz&
&yaach1(:,xyzzyaaad18)
endif
endif
enddo
enddo
if(xyzzyaabo1)then
do xyzzyaaae18=1,norbrange(xyzzyaact1)
xyzzyaaaf18=orbrange(1,xyzzyaaae18,xyzzyaact1)
xyzzyaaag18=orbrange(2,xyzzyaaae18,xyzzyaact1)
!$ call get_omp_index_range(xyzzyaaaf18,xyzzyaaag18)
!$ if(xyzzyaaag18<xyzzyaaaf18)cycle
if(val)then
orbval(xyzzyaaaf18:xyzzyaaag18,1)=dble(xyzzyaaap18(xyzzyaaaf18:xyzzyaa&
&ag18))
orbval(xyzzyaaaf18:xyzzyaaag18,2)=aimag(xyzzyaaap18(xyzzyaaaf18:xyzzya&
&aag18))
endif
if(fsd)then
orbgrad(1:3,xyzzyaaaf18:xyzzyaaag18,1)=dble(xyzzyaaar18(1:3,xyzzyaaaf1&
&8:xyzzyaaag18))
orbgrad(1:3,xyzzyaaaf18:xyzzyaaag18,2)=aimag(xyzzyaaar18(1:3,xyzzyaaaf&
&18:xyzzyaaag18))
orblap(xyzzyaaaf18:xyzzyaaag18,1)=dble(xyzzyaaaq18(xyzzyaaaf18:xyzzyaa&
&ag18))
orblap(xyzzyaaaf18:xyzzyaaag18,2)=aimag(xyzzyaaaq18(xyzzyaaaf18:xyzzya&
&aag18))
if(present(orbsderivs))then
orbsderivs(1:6,xyzzyaaaf18:xyzzyaaag18,1)=dble(xyzzyaaas18(1:6,xyzzyaa&
&af18:xyzzyaaag18))
orbsderivs(1:6,xyzzyaaaf18:xyzzyaaag18,2)=aimag(xyzzyaaas18(1:6,xyzzya&
&aaf18:xyzzyaaag18))
endif
endif
enddo
else
!$omp  barrier ! wait until all complex parts are available
do xyzzyaaae18=1,norbrange(xyzzyaacu1)
xyzzyaaaf18=orbrange(1,xyzzyaaae18,xyzzyaacu1)
xyzzyaaag18=orbrange(2,xyzzyaaae18,xyzzyaacu1)
!$ call get_omp_index_range(xyzzyaaaf18,xyzzyaaag18)
!$ if(xyzzyaaag18<xyzzyaaaf18)cycle
do xyzzyaaac18=xyzzyaaaf18,xyzzyaaag18
xyzzyaaad18=xyzzyaacw1(xyzzyaaac18)
if(xyzzyaadd1(xyzzyaaac18)==xyzzyaada1)then
if(val)orbval(xyzzyaaac18,1)=dble(xyzzyaaap18(xyzzyaaad18))
if(fsd)then
orbgrad(1:3,xyzzyaaac18,1)=dble(xyzzyaaar18(1:3,xyzzyaaad18))
orblap(xyzzyaaac18,1)=dble(xyzzyaaaq18(xyzzyaaad18))
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaac18,1)=dble(xyzzyaaas18(&
&1:6,xyzzyaaad18))
endif
elseif(xyzzyaadd1(xyzzyaaac18)==xyzzyaadb1)then
if(val)orbval(xyzzyaaac18,1)=aimag(xyzzyaaap18(xyzzyaaad18))
if(fsd)then
orbgrad(1:3,xyzzyaaac18,1)=aimag(xyzzyaaar18(1:3,xyzzyaaad18))
orblap(xyzzyaaac18,1)=aimag(xyzzyaaaq18(xyzzyaaad18))
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaac18,1)=aimag(xyzzyaaas18&
&(1:6,xyzzyaaad18))
endif
else
xyzzyaaak18=dble(xyzzyaade1(xyzzyaaac18))
xyzzyaaal18=aimag(xyzzyaade1(xyzzyaaac18))
if(val)orbval(xyzzyaaac18,1)=dble(xyzzyaaap18(xyzzyaaad18))*xyzzyaaak1&
&8+aimag(xyzzyaaap18(xyzzyaaad18))*xyzzyaaal18
if(fsd)then
orbgrad(1:3,xyzzyaaac18,1)=dble(xyzzyaaar18(1:3,xyzzyaaad18))*xyzzyaaa&
&k18+aimag(xyzzyaaar18(1:3,xyzzyaaad18))*xyzzyaaal18
orblap(xyzzyaaac18,1)=dble(xyzzyaaaq18(xyzzyaaad18))*xyzzyaaak18+aimag&
&(xyzzyaaaq18(xyzzyaaad18))*xyzzyaaal18
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaac18,1)=dble(xyzzyaaas18(&
&1:6,xyzzyaaad18))*xyzzyaaak18+aimag(xyzzyaaas18(1:6,xyzzyaaad18))*xyz&
&zyaaal18
endif
endif
enddo
enddo
endif
!$omp end parallel
end subroutine xyzzyaadw1
subroutine xyzzyaadx1(rvec,jspin,norb,xyzzyaacv1,norbrange,orbrange,va&
&l,fsd,blips_computed,r_frac,ir,f123,fg123,fl123,orbval,orbgrad,orblap&
&,orbsderivs,orb_m,orb_rmap)
use slaarnabi,      only : use_gpcc,cusp_wfdet_real
use slaarnabt, only : ddot,dgemv,dadd,dcopy
!$  use openmp_base, only : get_omp_index_range
implicit none
integer,intent(in) :: jspin,norb,xyzzyaacv1,norbrange(xyzzyaacv1),orbr&
&ange(2,xyzzyaaci1,xyzzyaacv1),ir(3)
integer,intent(out),optional :: orb_rmap(norb)
integer,intent(inout),optional :: orb_m
real(dp),intent(in) :: rvec(3),r_frac(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2),f123(64),fg123(64,3),fl123&
&(64,6)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd
logical,intent(inout) :: blips_computed
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19,xyzzyaaae19,xy&
&zzyaaaf19,xyzzyaaag19
real(dp) xyzzyaaah19(bw_norbc),xyzzyaaai19(3,bw_norbc),xyzzyaaaj19(6,b&
&w_norbc),xyzzyaaak19(bw_norbc),xyzzyaaal19(6,bw_norbc)
call xyzzyaaeg1(val,fsd,blips_computed,r_frac,f123,fg123,fl123)
!$omp parallel default(none) shared(ir,bw_norbc,xyzzyaacj1,norbrange,o&
!$omp &rbrange,xyzzyaacs1,xyzzyaaao1,xyzzyaaap1,xyzzyaaaa1,xyzzyaaak1,&
!$omp &xyzzyaaad1,xyzzyaaah1,val,fsd,painv,xyzzyaaal1,f123,fg123,fl123&
!$omp &,jspin,orbval,orbgrad,orblap,orbsderivs,xyzzyaaah19,xyzzyaaai19&
!$omp &,xyzzyaaak19,xyzzyaaaj19,xyzzyaaal19,xyzzyaaau1,xyzzyaaas1,use_&
!$omp &gpcc,xyzzyaaca1,xyzzyaacb1,xyzzyaacc1,xyzzyaacd1,rvec,single_pr&
!$omp &ecision_blips) private(xyzzyaaaa19,xyzzyaaac19,xyzzyaaad19,xyzz&
!$omp &yaaae19,xyzzyaaaf19,xyzzyaaab19,xyzzyaaag19)
if(fsd)then
if(single_precision_blips)then
do xyzzyaaac19=1,norbrange(xyzzyaacs1)
xyzzyaaad19=orbrange(1,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
xyzzyaaae19=orbrange(2,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
!$  call get_omp_index_range(xyzzyaaad19,xyzzyaaae19)
!$  if(xyzzyaaae19<xyzzyaaad19)cycle
call xyzzyaady1(ir,xyzzyaaaa1,xyzzyaaad19,xyzzyaaae19,bw_norbc,xyzzyaa&
&ao1,xyzzyaaau1)
enddo
else
do xyzzyaaac19=1,norbrange(xyzzyaacs1)
xyzzyaaad19=orbrange(1,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
xyzzyaaae19=orbrange(2,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
!$  call get_omp_index_range(xyzzyaaad19,xyzzyaaae19)
!$  if(xyzzyaaae19<xyzzyaaad19)cycle
call xyzzyaadz1(ir,xyzzyaaaa1,xyzzyaaad19,xyzzyaaae19,bw_norbc,xyzzyaa&
&ap1,xyzzyaaau1)
enddo
endif
do xyzzyaaac19=1,norbrange(xyzzyaacs1)
xyzzyaaad19=orbrange(1,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
xyzzyaaae19=orbrange(2,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
!$ call get_omp_index_range(xyzzyaaad19,xyzzyaaae19)
!$ if(xyzzyaaae19<xyzzyaaad19)cycle
xyzzyaaaf19=xyzzyaaae19-xyzzyaaad19+1
if(val)then
call dgemv('N',xyzzyaaaf19,64,1.d0,xyzzyaaau1(xyzzyaaad19,1),xyzzyaaas&
&1,f123(1),1,0.d0,xyzzyaaah19(xyzzyaaad19),1)
endif
do xyzzyaaaa19=1,3
call dgemv('N',xyzzyaaaf19,64,1.d0,xyzzyaaau1(xyzzyaaad19,1),xyzzyaaas&
&1,fg123(1,xyzzyaaaa19),1,0.d0,xyzzyaaai19(xyzzyaaaa19,xyzzyaaad19),3)
call dgemv('N',xyzzyaaaf19,64,1.d0,xyzzyaaau1(xyzzyaaad19,1),xyzzyaaas&
&1,fl123(1,xyzzyaaaa19),1,0.d0,xyzzyaaaj19(xyzzyaaaa19,xyzzyaaad19),6)
call dgemv('N',xyzzyaaaf19,64,1.d0,xyzzyaaau1(xyzzyaaad19,1),xyzzyaaas&
&1,fl123(1,xyzzyaaaa19+3),1,0.d0,xyzzyaaaj19(xyzzyaaaa19+3,xyzzyaaad19&
&),6)
enddo
do xyzzyaaab19=xyzzyaaad19,xyzzyaaae19
xyzzyaaai19(1:3,xyzzyaaab19)=(/ddot(3,painv(1,1),3,xyzzyaaai19(1,xyzzy&
&aaab19),1),ddot(3,painv(2,1),3,xyzzyaaai19(1,xyzzyaaab19),1),ddot(3,p&
&ainv(3,1),3,xyzzyaaai19(1,xyzzyaaab19),1)/)
xyzzyaaak19(xyzzyaaab19)=ddot(6,xyzzyaaaj19(1,xyzzyaaab19),1,xyzzyaaak&
&1(1),1)
if(present(orbsderivs))then
do xyzzyaaaa19=1,6
xyzzyaaal19(xyzzyaaaa19,xyzzyaaab19)=ddot(6,xyzzyaaaj19(1,xyzzyaaab19)&
&,1,xyzzyaaal1(1,xyzzyaaaa19),1)
enddo
endif
enddo
enddo
elseif(val)then
if(single_precision_blips)then
do xyzzyaaac19=1,norbrange(xyzzyaacs1)
xyzzyaaad19=orbrange(1,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
xyzzyaaae19=orbrange(2,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
!$  call get_omp_index_range(xyzzyaaad19,xyzzyaaae19)
!$  if(xyzzyaaae19<xyzzyaaad19)cycle
call xyzzyaaec1(ir,xyzzyaaaa1,xyzzyaaad19,xyzzyaaae19,bw_norbc,xyzzyaa&
&ao1,f123,xyzzyaaah19(xyzzyaaad19))
enddo
else
do xyzzyaaac19=1,norbrange(xyzzyaacs1)
xyzzyaaad19=orbrange(1,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
xyzzyaaae19=orbrange(2,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
!$  call get_omp_index_range(xyzzyaaad19,xyzzyaaae19)
!$  if(xyzzyaaae19<xyzzyaaad19)cycle
call xyzzyaaed1(ir,xyzzyaaaa1,xyzzyaaad19,xyzzyaaae19,bw_norbc,xyzzyaa&
&ap1,f123,xyzzyaaah19(xyzzyaaad19))
enddo
endif
endif
if(use_gpcc)then
if(present(orbsderivs))then
call cusp_wfdet_real(rvec,jspin,val,fsd,xyzzyaaca1,xyzzyaacb1,xyzzyaac&
&c1,xyzzyaacd1)
else
call cusp_wfdet_real(rvec,jspin,val,fsd,xyzzyaaca1,xyzzyaacb1,xyzzyaac&
&c1)
endif
!$omp  barrier ! ensure that all cusp correction are done
endif
do xyzzyaaac19=1,norbrange(xyzzyaacs1)
xyzzyaaad19=orbrange(1,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
xyzzyaaae19=orbrange(2,xyzzyaaac19,xyzzyaacs1)-xyzzyaacj1
!$ call get_omp_index_range(xyzzyaaad19,xyzzyaaae19)
!$ if(xyzzyaaae19<xyzzyaaad19)cycle
xyzzyaaaf19=xyzzyaaae19-xyzzyaaad19+1
xyzzyaaag19=xyzzyaaad19+xyzzyaacj1
if(val)then
call dcopy(xyzzyaaaf19,xyzzyaaah19(xyzzyaaad19),1,orbval(xyzzyaaag19,1&
&),1)
if(use_gpcc)call dadd(xyzzyaaaf19,xyzzyaaca1(xyzzyaaad19),1,orbval(xyz&
&zyaaag19,1),1)
endif
if(fsd)then
call dcopy(3*xyzzyaaaf19,xyzzyaaai19(1,xyzzyaaad19),1,orbgrad(1,xyzzya&
&aag19,1),1)
call dcopy(xyzzyaaaf19,xyzzyaaak19(xyzzyaaad19),1,orblap(xyzzyaaag19,1&
&),1)
if(present(orbsderivs))call dcopy(6*xyzzyaaaf19,xyzzyaaal19(1,xyzzyaaa&
&d19),1,orbsderivs(1,xyzzyaaag19,1),1)
if(use_gpcc)then
call dadd(3*xyzzyaaaf19,xyzzyaacb1(1,xyzzyaaad19),1,orbgrad(1,xyzzyaaa&
&g19,1),1)
call dadd(xyzzyaaaf19,xyzzyaacc1(xyzzyaaad19),1,orblap(xyzzyaaag19,1),&
&1)
if(present(orbsderivs))then
call dadd(6*xyzzyaaaf19,xyzzyaacd1(1,xyzzyaaad19),1,orbsderivs(1,xyzzy&
&aaag19,1),1)
endif
endif
endif
enddo
!$omp end parallel
if(present(orb_m))then
do xyzzyaaac19=1,norbrange(xyzzyaacs1)
xyzzyaaad19=orbrange(1,xyzzyaaac19,xyzzyaacs1)
xyzzyaaae19=orbrange(2,xyzzyaaac19,xyzzyaacs1)
do xyzzyaaab19=xyzzyaaad19,xyzzyaaae19
orb_m=orb_m+1
orb_rmap(orb_m)=xyzzyaaab19
enddo
enddo
endif
end subroutine xyzzyaadx1
subroutine xyzzyaady1(ir,xyzzyaaaa1,iorb1,iorb2,norb,bcoef,array)
implicit none
integer,intent(in) :: ir(3),xyzzyaaaa1(3),iorb1,iorb2,norb
real(sp),intent(in) :: bcoef(norb,*)
real(dp),intent(inout) :: array(norb,64)
integer xyzzyaaaa20,xyzzyaaab20,xyzzyaaac20,xyzzyaaad20,xyzzyaaae20,xy&
&zzyaaaf20,xyzzyaaag20,xyzzyaaah20,xyzzyaaai20,xyzzyaaaj20,xyzzyaaak20&
&,xyzzyaaal20,xyzzyaaam20,xyzzyaaan20,xyzzyaaao20,xyzzyaaap20,xyzzyaaa&
&q20,xyzzyaaar20,xyzzyaaas20,xyzzyaaat20,xyzzyaaau20,xyzzyaaav20
xyzzyaaan20=ir(1)
xyzzyaaao20=ir(2)
xyzzyaaap20=ir(3)
xyzzyaaaq20=xyzzyaaaa1(1)
xyzzyaaar20=xyzzyaaaa1(2)
xyzzyaaas20=xyzzyaaaa1(3)
xyzzyaaal20=xyzzyaaas20
xyzzyaaam20=xyzzyaaal20+xyzzyaaal20+xyzzyaaal20+xyzzyaaal20
xyzzyaaaj20=xyzzyaaar20*xyzzyaaal20
xyzzyaaak20=xyzzyaaaj20+xyzzyaaaj20+xyzzyaaaj20+xyzzyaaaj20
xyzzyaaai20=xyzzyaaaq20*xyzzyaaaj20
xyzzyaaah20=(xyzzyaaan20-1)*xyzzyaaaj20+(xyzzyaaao20-1)*xyzzyaaal20+xy&
&zzyaaap20
xyzzyaaag20=0
if(xyzzyaaap20+3>=xyzzyaaas20)then
xyzzyaaad20=xyzzyaaah20
xyzzyaaat20=xyzzyaaan20-1
do xyzzyaaaa20=xyzzyaaan20,xyzzyaaan20+3
xyzzyaaad20=xyzzyaaad20+xyzzyaaaj20
xyzzyaaat20=xyzzyaaat20+1
if(xyzzyaaat20==xyzzyaaaq20)then
xyzzyaaad20=xyzzyaaad20-xyzzyaaai20
xyzzyaaat20=0
endif
xyzzyaaae20=xyzzyaaad20
xyzzyaaau20=xyzzyaaao20-1
do xyzzyaaab20=xyzzyaaao20,xyzzyaaao20+3
xyzzyaaae20=xyzzyaaae20+xyzzyaaal20
xyzzyaaau20=xyzzyaaau20+1
if(xyzzyaaau20==xyzzyaaar20)then
xyzzyaaae20=xyzzyaaae20-xyzzyaaaj20
xyzzyaaau20=0
endif
xyzzyaaaf20=xyzzyaaae20
xyzzyaaav20=xyzzyaaap20-1
do xyzzyaaac20=xyzzyaaap20,xyzzyaaap20+3
xyzzyaaaf20=xyzzyaaaf20+1
xyzzyaaav20=xyzzyaaav20+1
if(xyzzyaaav20==xyzzyaaas20)then
xyzzyaaaf20=xyzzyaaaf20-xyzzyaaal20
xyzzyaaav20=0
endif
xyzzyaaag20=xyzzyaaag20+1
array(iorb1:iorb2,xyzzyaaag20)=bcoef(iorb1:iorb2,xyzzyaaaf20)
enddo
enddo
enddo
elseif(xyzzyaaao20+3>=xyzzyaaar20)then
xyzzyaaad20=xyzzyaaah20
xyzzyaaat20=xyzzyaaan20-1
do xyzzyaaaa20=xyzzyaaan20,xyzzyaaan20+3
xyzzyaaad20=xyzzyaaad20+xyzzyaaaj20
xyzzyaaat20=xyzzyaaat20+1
if(xyzzyaaat20==xyzzyaaaq20)then
xyzzyaaad20=xyzzyaaad20-xyzzyaaai20
xyzzyaaat20=0
endif
xyzzyaaae20=xyzzyaaad20
xyzzyaaau20=xyzzyaaao20-1
do xyzzyaaab20=xyzzyaaao20,xyzzyaaao20+3
xyzzyaaae20=xyzzyaaae20+xyzzyaaal20
xyzzyaaau20=xyzzyaaau20+1
if(xyzzyaaau20==xyzzyaaar20)then
xyzzyaaae20=xyzzyaaae20-xyzzyaaaj20
xyzzyaaau20=0
endif
array(iorb1:iorb2,xyzzyaaag20+1)=bcoef(iorb1:iorb2,xyzzyaaae20+1)
array(iorb1:iorb2,xyzzyaaag20+2)=bcoef(iorb1:iorb2,xyzzyaaae20+2)
array(iorb1:iorb2,xyzzyaaag20+3)=bcoef(iorb1:iorb2,xyzzyaaae20+3)
array(iorb1:iorb2,xyzzyaaag20+4)=bcoef(iorb1:iorb2,xyzzyaaae20+4)
xyzzyaaag20=xyzzyaaag20+4
enddo
enddo
elseif(xyzzyaaan20+3>=xyzzyaaaq20)then
xyzzyaaad20=xyzzyaaah20
xyzzyaaat20=xyzzyaaan20-1
do xyzzyaaaa20=xyzzyaaan20,xyzzyaaan20+3
xyzzyaaad20=xyzzyaaad20+xyzzyaaaj20
xyzzyaaat20=xyzzyaaat20+1
if(xyzzyaaat20==xyzzyaaaq20)then
xyzzyaaad20=xyzzyaaad20-xyzzyaaai20
xyzzyaaat20=0
endif
do xyzzyaaae20=xyzzyaaad20+xyzzyaaal20,xyzzyaaad20+xyzzyaaam20,xyzzyaa&
&al20
array(iorb1:iorb2,xyzzyaaag20+1)=bcoef(iorb1:iorb2,xyzzyaaae20+1)
array(iorb1:iorb2,xyzzyaaag20+2)=bcoef(iorb1:iorb2,xyzzyaaae20+2)
array(iorb1:iorb2,xyzzyaaag20+3)=bcoef(iorb1:iorb2,xyzzyaaae20+3)
array(iorb1:iorb2,xyzzyaaag20+4)=bcoef(iorb1:iorb2,xyzzyaaae20+4)
xyzzyaaag20=xyzzyaaag20+4
enddo
enddo
else
do xyzzyaaad20=xyzzyaaah20+xyzzyaaaj20,xyzzyaaah20+xyzzyaaak20,xyzzyaa&
&aj20
do xyzzyaaae20=xyzzyaaad20+xyzzyaaal20,xyzzyaaad20+xyzzyaaam20,xyzzyaa&
&al20
array(iorb1:iorb2,xyzzyaaag20+1)=bcoef(iorb1:iorb2,xyzzyaaae20+1)
array(iorb1:iorb2,xyzzyaaag20+2)=bcoef(iorb1:iorb2,xyzzyaaae20+2)
array(iorb1:iorb2,xyzzyaaag20+3)=bcoef(iorb1:iorb2,xyzzyaaae20+3)
array(iorb1:iorb2,xyzzyaaag20+4)=bcoef(iorb1:iorb2,xyzzyaaae20+4)
xyzzyaaag20=xyzzyaaag20+4
enddo
enddo
endif
end subroutine xyzzyaady1
subroutine xyzzyaadz1(ir,xyzzyaaaa1,iorb1,iorb2,norb,bcoef,array)
implicit none
integer,intent(in) :: ir(3),xyzzyaaaa1(3),iorb1,iorb2,norb
real(dp),intent(in) :: bcoef(norb,*)
real(dp),intent(inout) :: array(norb,64)
integer xyzzyaaaa21,xyzzyaaab21,xyzzyaaac21,xyzzyaaad21,xyzzyaaae21,xy&
&zzyaaaf21,xyzzyaaag21,xyzzyaaah21,xyzzyaaai21,xyzzyaaaj21,xyzzyaaak21&
&,xyzzyaaal21,xyzzyaaam21,xyzzyaaan21,xyzzyaaao21,xyzzyaaap21,xyzzyaaa&
&q21,xyzzyaaar21,xyzzyaaas21,xyzzyaaat21,xyzzyaaau21,xyzzyaaav21
xyzzyaaan21=ir(1)
xyzzyaaao21=ir(2)
xyzzyaaap21=ir(3)
xyzzyaaaq21=xyzzyaaaa1(1)
xyzzyaaar21=xyzzyaaaa1(2)
xyzzyaaas21=xyzzyaaaa1(3)
xyzzyaaal21=xyzzyaaas21
xyzzyaaam21=xyzzyaaal21+xyzzyaaal21+xyzzyaaal21+xyzzyaaal21
xyzzyaaaj21=xyzzyaaar21*xyzzyaaal21
xyzzyaaak21=xyzzyaaaj21+xyzzyaaaj21+xyzzyaaaj21+xyzzyaaaj21
xyzzyaaai21=xyzzyaaaq21*xyzzyaaaj21
xyzzyaaah21=(xyzzyaaan21-1)*xyzzyaaaj21+(xyzzyaaao21-1)*xyzzyaaal21+xy&
&zzyaaap21
xyzzyaaag21=0
if(xyzzyaaap21+3>=xyzzyaaas21)then
xyzzyaaad21=xyzzyaaah21
xyzzyaaat21=xyzzyaaan21-1
do xyzzyaaaa21=xyzzyaaan21,xyzzyaaan21+3
xyzzyaaad21=xyzzyaaad21+xyzzyaaaj21
xyzzyaaat21=xyzzyaaat21+1
if(xyzzyaaat21==xyzzyaaaq21)then
xyzzyaaad21=xyzzyaaad21-xyzzyaaai21
xyzzyaaat21=0
endif
xyzzyaaae21=xyzzyaaad21
xyzzyaaau21=xyzzyaaao21-1
do xyzzyaaab21=xyzzyaaao21,xyzzyaaao21+3
xyzzyaaae21=xyzzyaaae21+xyzzyaaal21
xyzzyaaau21=xyzzyaaau21+1
if(xyzzyaaau21==xyzzyaaar21)then
xyzzyaaae21=xyzzyaaae21-xyzzyaaaj21
xyzzyaaau21=0
endif
xyzzyaaaf21=xyzzyaaae21
xyzzyaaav21=xyzzyaaap21-1
do xyzzyaaac21=xyzzyaaap21,xyzzyaaap21+3
xyzzyaaaf21=xyzzyaaaf21+1
xyzzyaaav21=xyzzyaaav21+1
if(xyzzyaaav21==xyzzyaaas21)then
xyzzyaaaf21=xyzzyaaaf21-xyzzyaaal21
xyzzyaaav21=0
endif
xyzzyaaag21=xyzzyaaag21+1
array(iorb1:iorb2,xyzzyaaag21)=bcoef(iorb1:iorb2,xyzzyaaaf21)
enddo
enddo
enddo
elseif(xyzzyaaao21+3>=xyzzyaaar21)then
xyzzyaaad21=xyzzyaaah21
xyzzyaaat21=xyzzyaaan21-1
do xyzzyaaaa21=xyzzyaaan21,xyzzyaaan21+3
xyzzyaaad21=xyzzyaaad21+xyzzyaaaj21
xyzzyaaat21=xyzzyaaat21+1
if(xyzzyaaat21==xyzzyaaaq21)then
xyzzyaaad21=xyzzyaaad21-xyzzyaaai21
xyzzyaaat21=0
endif
xyzzyaaae21=xyzzyaaad21
xyzzyaaau21=xyzzyaaao21-1
do xyzzyaaab21=xyzzyaaao21,xyzzyaaao21+3
xyzzyaaae21=xyzzyaaae21+xyzzyaaal21
xyzzyaaau21=xyzzyaaau21+1
if(xyzzyaaau21==xyzzyaaar21)then
xyzzyaaae21=xyzzyaaae21-xyzzyaaaj21
xyzzyaaau21=0
endif
array(iorb1:iorb2,xyzzyaaag21+1)=bcoef(iorb1:iorb2,xyzzyaaae21+1)
array(iorb1:iorb2,xyzzyaaag21+2)=bcoef(iorb1:iorb2,xyzzyaaae21+2)
array(iorb1:iorb2,xyzzyaaag21+3)=bcoef(iorb1:iorb2,xyzzyaaae21+3)
array(iorb1:iorb2,xyzzyaaag21+4)=bcoef(iorb1:iorb2,xyzzyaaae21+4)
xyzzyaaag21=xyzzyaaag21+4
enddo
enddo
elseif(xyzzyaaan21+3>=xyzzyaaaq21)then
xyzzyaaad21=xyzzyaaah21
xyzzyaaat21=xyzzyaaan21-1
do xyzzyaaaa21=xyzzyaaan21,xyzzyaaan21+3
xyzzyaaad21=xyzzyaaad21+xyzzyaaaj21
xyzzyaaat21=xyzzyaaat21+1
if(xyzzyaaat21==xyzzyaaaq21)then
xyzzyaaad21=xyzzyaaad21-xyzzyaaai21
xyzzyaaat21=0
endif
do xyzzyaaae21=xyzzyaaad21+xyzzyaaal21,xyzzyaaad21+xyzzyaaam21,xyzzyaa&
&al21
array(iorb1:iorb2,xyzzyaaag21+1)=bcoef(iorb1:iorb2,xyzzyaaae21+1)
array(iorb1:iorb2,xyzzyaaag21+2)=bcoef(iorb1:iorb2,xyzzyaaae21+2)
array(iorb1:iorb2,xyzzyaaag21+3)=bcoef(iorb1:iorb2,xyzzyaaae21+3)
array(iorb1:iorb2,xyzzyaaag21+4)=bcoef(iorb1:iorb2,xyzzyaaae21+4)
xyzzyaaag21=xyzzyaaag21+4
enddo
enddo
else
do xyzzyaaad21=xyzzyaaah21+xyzzyaaaj21,xyzzyaaah21+xyzzyaaak21,xyzzyaa&
&aj21
do xyzzyaaae21=xyzzyaaad21+xyzzyaaal21,xyzzyaaad21+xyzzyaaam21,xyzzyaa&
&al21
array(iorb1:iorb2,xyzzyaaag21+1)=bcoef(iorb1:iorb2,xyzzyaaae21+1)
array(iorb1:iorb2,xyzzyaaag21+2)=bcoef(iorb1:iorb2,xyzzyaaae21+2)
array(iorb1:iorb2,xyzzyaaag21+3)=bcoef(iorb1:iorb2,xyzzyaaae21+3)
array(iorb1:iorb2,xyzzyaaag21+4)=bcoef(iorb1:iorb2,xyzzyaaae21+4)
xyzzyaaag21=xyzzyaaag21+4
enddo
enddo
endif
end subroutine xyzzyaadz1
subroutine xyzzyaaea1(ir,xyzzyaaaa1,iorb1,iorb2,norb,bcoef,array)
implicit none
integer,intent(in) :: ir(3),xyzzyaaaa1(3),iorb1,iorb2,norb
complex(sp),intent(in) :: bcoef(norb,*)
complex(dp),intent(inout) :: array(norb,64)
integer xyzzyaaaa22,xyzzyaaab22,xyzzyaaac22,xyzzyaaad22,xyzzyaaae22,xy&
&zzyaaaf22,xyzzyaaag22,xyzzyaaah22,xyzzyaaai22,xyzzyaaaj22,xyzzyaaak22&
&,xyzzyaaal22,xyzzyaaam22,xyzzyaaan22,xyzzyaaao22,xyzzyaaap22,xyzzyaaa&
&q22,xyzzyaaar22,xyzzyaaas22,xyzzyaaat22,xyzzyaaau22,xyzzyaaav22
xyzzyaaan22=ir(1)
xyzzyaaao22=ir(2)
xyzzyaaap22=ir(3)
xyzzyaaaq22=xyzzyaaaa1(1)
xyzzyaaar22=xyzzyaaaa1(2)
xyzzyaaas22=xyzzyaaaa1(3)
xyzzyaaal22=xyzzyaaas22
xyzzyaaam22=xyzzyaaal22+xyzzyaaal22+xyzzyaaal22+xyzzyaaal22
xyzzyaaaj22=xyzzyaaar22*xyzzyaaal22
xyzzyaaak22=xyzzyaaaj22+xyzzyaaaj22+xyzzyaaaj22+xyzzyaaaj22
xyzzyaaai22=xyzzyaaaq22*xyzzyaaaj22
xyzzyaaah22=(xyzzyaaan22-1)*xyzzyaaaj22+(xyzzyaaao22-1)*xyzzyaaal22+xy&
&zzyaaap22
xyzzyaaag22=0
if(xyzzyaaap22+3>=xyzzyaaas22)then
xyzzyaaad22=xyzzyaaah22
xyzzyaaat22=xyzzyaaan22-1
do xyzzyaaaa22=xyzzyaaan22,xyzzyaaan22+3
xyzzyaaad22=xyzzyaaad22+xyzzyaaaj22
xyzzyaaat22=xyzzyaaat22+1
if(xyzzyaaat22==xyzzyaaaq22)then
xyzzyaaad22=xyzzyaaad22-xyzzyaaai22
xyzzyaaat22=0
endif
xyzzyaaae22=xyzzyaaad22
xyzzyaaau22=xyzzyaaao22-1
do xyzzyaaab22=xyzzyaaao22,xyzzyaaao22+3
xyzzyaaae22=xyzzyaaae22+xyzzyaaal22
xyzzyaaau22=xyzzyaaau22+1
if(xyzzyaaau22==xyzzyaaar22)then
xyzzyaaae22=xyzzyaaae22-xyzzyaaaj22
xyzzyaaau22=0
endif
xyzzyaaaf22=xyzzyaaae22
xyzzyaaav22=xyzzyaaap22-1
do xyzzyaaac22=xyzzyaaap22,xyzzyaaap22+3
xyzzyaaaf22=xyzzyaaaf22+1
xyzzyaaav22=xyzzyaaav22+1
if(xyzzyaaav22==xyzzyaaas22)then
xyzzyaaaf22=xyzzyaaaf22-xyzzyaaal22
xyzzyaaav22=0
endif
xyzzyaaag22=xyzzyaaag22+1
array(iorb1:iorb2,xyzzyaaag22)=bcoef(iorb1:iorb2,xyzzyaaaf22)
enddo
enddo
enddo
elseif(xyzzyaaao22+3>=xyzzyaaar22)then
xyzzyaaad22=xyzzyaaah22
xyzzyaaat22=xyzzyaaan22-1
do xyzzyaaaa22=xyzzyaaan22,xyzzyaaan22+3
xyzzyaaad22=xyzzyaaad22+xyzzyaaaj22
xyzzyaaat22=xyzzyaaat22+1
if(xyzzyaaat22==xyzzyaaaq22)then
xyzzyaaad22=xyzzyaaad22-xyzzyaaai22
xyzzyaaat22=0
endif
xyzzyaaae22=xyzzyaaad22
xyzzyaaau22=xyzzyaaao22-1
do xyzzyaaab22=xyzzyaaao22,xyzzyaaao22+3
xyzzyaaae22=xyzzyaaae22+xyzzyaaal22
xyzzyaaau22=xyzzyaaau22+1
if(xyzzyaaau22==xyzzyaaar22)then
xyzzyaaae22=xyzzyaaae22-xyzzyaaaj22
xyzzyaaau22=0
endif
array(iorb1:iorb2,xyzzyaaag22+1)=bcoef(iorb1:iorb2,xyzzyaaae22+1)
array(iorb1:iorb2,xyzzyaaag22+2)=bcoef(iorb1:iorb2,xyzzyaaae22+2)
array(iorb1:iorb2,xyzzyaaag22+3)=bcoef(iorb1:iorb2,xyzzyaaae22+3)
array(iorb1:iorb2,xyzzyaaag22+4)=bcoef(iorb1:iorb2,xyzzyaaae22+4)
xyzzyaaag22=xyzzyaaag22+4
enddo
enddo
elseif(xyzzyaaan22+3>=xyzzyaaaq22)then
xyzzyaaad22=xyzzyaaah22
xyzzyaaat22=xyzzyaaan22-1
do xyzzyaaaa22=xyzzyaaan22,xyzzyaaan22+3
xyzzyaaad22=xyzzyaaad22+xyzzyaaaj22
xyzzyaaat22=xyzzyaaat22+1
if(xyzzyaaat22==xyzzyaaaq22)then
xyzzyaaad22=xyzzyaaad22-xyzzyaaai22
xyzzyaaat22=0
endif
do xyzzyaaae22=xyzzyaaad22+xyzzyaaal22,xyzzyaaad22+xyzzyaaam22,xyzzyaa&
&al22
array(iorb1:iorb2,xyzzyaaag22+1)=bcoef(iorb1:iorb2,xyzzyaaae22+1)
array(iorb1:iorb2,xyzzyaaag22+2)=bcoef(iorb1:iorb2,xyzzyaaae22+2)
array(iorb1:iorb2,xyzzyaaag22+3)=bcoef(iorb1:iorb2,xyzzyaaae22+3)
array(iorb1:iorb2,xyzzyaaag22+4)=bcoef(iorb1:iorb2,xyzzyaaae22+4)
xyzzyaaag22=xyzzyaaag22+4
enddo
enddo
else
do xyzzyaaad22=xyzzyaaah22+xyzzyaaaj22,xyzzyaaah22+xyzzyaaak22,xyzzyaa&
&aj22
do xyzzyaaae22=xyzzyaaad22+xyzzyaaal22,xyzzyaaad22+xyzzyaaam22,xyzzyaa&
&al22
array(iorb1:iorb2,xyzzyaaag22+1)=bcoef(iorb1:iorb2,xyzzyaaae22+1)
array(iorb1:iorb2,xyzzyaaag22+2)=bcoef(iorb1:iorb2,xyzzyaaae22+2)
array(iorb1:iorb2,xyzzyaaag22+3)=bcoef(iorb1:iorb2,xyzzyaaae22+3)
array(iorb1:iorb2,xyzzyaaag22+4)=bcoef(iorb1:iorb2,xyzzyaaae22+4)
xyzzyaaag22=xyzzyaaag22+4
enddo
enddo
endif
end subroutine xyzzyaaea1
subroutine xyzzyaaeb1(ir,xyzzyaaaa1,iorb1,iorb2,norb,bcoef,array)
implicit none
integer,intent(in) :: ir(3),xyzzyaaaa1(3),iorb1,iorb2,norb
complex(dp),intent(in) :: bcoef(norb,*)
complex(dp),intent(inout) :: array(norb,64)
integer xyzzyaaaa23,xyzzyaaab23,xyzzyaaac23,xyzzyaaad23,xyzzyaaae23,xy&
&zzyaaaf23,xyzzyaaag23,xyzzyaaah23,xyzzyaaai23,xyzzyaaaj23,xyzzyaaak23&
&,xyzzyaaal23,xyzzyaaam23,xyzzyaaan23,xyzzyaaao23,xyzzyaaap23,xyzzyaaa&
&q23,xyzzyaaar23,xyzzyaaas23,xyzzyaaat23,xyzzyaaau23,xyzzyaaav23
xyzzyaaan23=ir(1)
xyzzyaaao23=ir(2)
xyzzyaaap23=ir(3)
xyzzyaaaq23=xyzzyaaaa1(1)
xyzzyaaar23=xyzzyaaaa1(2)
xyzzyaaas23=xyzzyaaaa1(3)
xyzzyaaal23=xyzzyaaas23
xyzzyaaam23=xyzzyaaal23+xyzzyaaal23+xyzzyaaal23+xyzzyaaal23
xyzzyaaaj23=xyzzyaaar23*xyzzyaaal23
xyzzyaaak23=xyzzyaaaj23+xyzzyaaaj23+xyzzyaaaj23+xyzzyaaaj23
xyzzyaaai23=xyzzyaaaq23*xyzzyaaaj23
xyzzyaaah23=(xyzzyaaan23-1)*xyzzyaaaj23+(xyzzyaaao23-1)*xyzzyaaal23+xy&
&zzyaaap23
xyzzyaaag23=0
if(xyzzyaaap23+3>=xyzzyaaas23)then
xyzzyaaad23=xyzzyaaah23
xyzzyaaat23=xyzzyaaan23-1
do xyzzyaaaa23=xyzzyaaan23,xyzzyaaan23+3
xyzzyaaad23=xyzzyaaad23+xyzzyaaaj23
xyzzyaaat23=xyzzyaaat23+1
if(xyzzyaaat23==xyzzyaaaq23)then
xyzzyaaad23=xyzzyaaad23-xyzzyaaai23
xyzzyaaat23=0
endif
xyzzyaaae23=xyzzyaaad23
xyzzyaaau23=xyzzyaaao23-1
do xyzzyaaab23=xyzzyaaao23,xyzzyaaao23+3
xyzzyaaae23=xyzzyaaae23+xyzzyaaal23
xyzzyaaau23=xyzzyaaau23+1
if(xyzzyaaau23==xyzzyaaar23)then
xyzzyaaae23=xyzzyaaae23-xyzzyaaaj23
xyzzyaaau23=0
endif
xyzzyaaaf23=xyzzyaaae23
xyzzyaaav23=xyzzyaaap23-1
do xyzzyaaac23=xyzzyaaap23,xyzzyaaap23+3
xyzzyaaaf23=xyzzyaaaf23+1
xyzzyaaav23=xyzzyaaav23+1
if(xyzzyaaav23==xyzzyaaas23)then
xyzzyaaaf23=xyzzyaaaf23-xyzzyaaal23
xyzzyaaav23=0
endif
xyzzyaaag23=xyzzyaaag23+1
array(iorb1:iorb2,xyzzyaaag23)=bcoef(iorb1:iorb2,xyzzyaaaf23)
enddo
enddo
enddo
elseif(xyzzyaaao23+3>=xyzzyaaar23)then
xyzzyaaad23=xyzzyaaah23
xyzzyaaat23=xyzzyaaan23-1
do xyzzyaaaa23=xyzzyaaan23,xyzzyaaan23+3
xyzzyaaad23=xyzzyaaad23+xyzzyaaaj23
xyzzyaaat23=xyzzyaaat23+1
if(xyzzyaaat23==xyzzyaaaq23)then
xyzzyaaad23=xyzzyaaad23-xyzzyaaai23
xyzzyaaat23=0
endif
xyzzyaaae23=xyzzyaaad23
xyzzyaaau23=xyzzyaaao23-1
do xyzzyaaab23=xyzzyaaao23,xyzzyaaao23+3
xyzzyaaae23=xyzzyaaae23+xyzzyaaal23
xyzzyaaau23=xyzzyaaau23+1
if(xyzzyaaau23==xyzzyaaar23)then
xyzzyaaae23=xyzzyaaae23-xyzzyaaaj23
xyzzyaaau23=0
endif
array(iorb1:iorb2,xyzzyaaag23+1)=bcoef(iorb1:iorb2,xyzzyaaae23+1)
array(iorb1:iorb2,xyzzyaaag23+2)=bcoef(iorb1:iorb2,xyzzyaaae23+2)
array(iorb1:iorb2,xyzzyaaag23+3)=bcoef(iorb1:iorb2,xyzzyaaae23+3)
array(iorb1:iorb2,xyzzyaaag23+4)=bcoef(iorb1:iorb2,xyzzyaaae23+4)
xyzzyaaag23=xyzzyaaag23+4
enddo
enddo
elseif(xyzzyaaan23+3>=xyzzyaaaq23)then
xyzzyaaad23=xyzzyaaah23
xyzzyaaat23=xyzzyaaan23-1
do xyzzyaaaa23=xyzzyaaan23,xyzzyaaan23+3
xyzzyaaad23=xyzzyaaad23+xyzzyaaaj23
xyzzyaaat23=xyzzyaaat23+1
if(xyzzyaaat23==xyzzyaaaq23)then
xyzzyaaad23=xyzzyaaad23-xyzzyaaai23
xyzzyaaat23=0
endif
do xyzzyaaae23=xyzzyaaad23+xyzzyaaal23,xyzzyaaad23+xyzzyaaam23,xyzzyaa&
&al23
array(iorb1:iorb2,xyzzyaaag23+1)=bcoef(iorb1:iorb2,xyzzyaaae23+1)
array(iorb1:iorb2,xyzzyaaag23+2)=bcoef(iorb1:iorb2,xyzzyaaae23+2)
array(iorb1:iorb2,xyzzyaaag23+3)=bcoef(iorb1:iorb2,xyzzyaaae23+3)
array(iorb1:iorb2,xyzzyaaag23+4)=bcoef(iorb1:iorb2,xyzzyaaae23+4)
xyzzyaaag23=xyzzyaaag23+4
enddo
enddo
else
do xyzzyaaad23=xyzzyaaah23+xyzzyaaaj23,xyzzyaaah23+xyzzyaaak23,xyzzyaa&
&aj23
do xyzzyaaae23=xyzzyaaad23+xyzzyaaal23,xyzzyaaad23+xyzzyaaam23,xyzzyaa&
&al23
array(iorb1:iorb2,xyzzyaaag23+1)=bcoef(iorb1:iorb2,xyzzyaaae23+1)
array(iorb1:iorb2,xyzzyaaag23+2)=bcoef(iorb1:iorb2,xyzzyaaae23+2)
array(iorb1:iorb2,xyzzyaaag23+3)=bcoef(iorb1:iorb2,xyzzyaaae23+3)
array(iorb1:iorb2,xyzzyaaag23+4)=bcoef(iorb1:iorb2,xyzzyaaae23+4)
xyzzyaaag23=xyzzyaaag23+4
enddo
enddo
endif
end subroutine xyzzyaaeb1
subroutine xyzzyaaec1(ir,xyzzyaaaa1,iorb1,iorb2,norb,bcoef,vector,mult&
&)
implicit none
integer,intent(in) :: ir(3),xyzzyaaaa1(3),iorb1,iorb2,norb
real(sp),intent(in) :: bcoef(norb,*)
real(dp),intent(in) :: vector(64)
real(dp),intent(inout) :: mult(1:iorb2-iorb1+1)
integer xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaae24,xy&
&zzyaaaf24,xyzzyaaag24,xyzzyaaah24,xyzzyaaai24,xyzzyaaaj24,xyzzyaaak24&
&,xyzzyaaal24,xyzzyaaam24,xyzzyaaan24,xyzzyaaao24,xyzzyaaap24,xyzzyaaa&
&q24,xyzzyaaar24,xyzzyaaas24,xyzzyaaat24,xyzzyaaau24,xyzzyaaav24,xyzzy&
&aaaw24,xyzzyaaax24,xyzzyaaay24,xyzzyaaaz24,xyzzyaaba24
real(dp) xyzzyaabb24
mult=0.d0
xyzzyaaai24=iorb1-1
xyzzyaaaj24=iorb2-xyzzyaaai24
xyzzyaaaq24=ir(1)
xyzzyaaar24=ir(2)
xyzzyaaas24=ir(3)
xyzzyaaat24=xyzzyaaaa1(1)
xyzzyaaau24=xyzzyaaaa1(2)
xyzzyaaav24=xyzzyaaaa1(3)
xyzzyaaao24=xyzzyaaav24
xyzzyaaap24=xyzzyaaao24+xyzzyaaao24+xyzzyaaao24+xyzzyaaao24
xyzzyaaam24=xyzzyaaau24*xyzzyaaao24
xyzzyaaan24=xyzzyaaam24+xyzzyaaam24+xyzzyaaam24+xyzzyaaam24
xyzzyaaal24=xyzzyaaat24*xyzzyaaam24
xyzzyaaah24=(xyzzyaaaq24-1)*xyzzyaaam24+(xyzzyaaar24-1)*xyzzyaaao24+xy&
&zzyaaas24
xyzzyaaag24=0
if(xyzzyaaas24+3>=xyzzyaaav24)then
xyzzyaaad24=xyzzyaaah24
xyzzyaaay24=xyzzyaaaq24-1
do xyzzyaaaa24=xyzzyaaaq24,xyzzyaaaq24+3
xyzzyaaad24=xyzzyaaad24+xyzzyaaam24
xyzzyaaay24=xyzzyaaay24+1
if(xyzzyaaay24==xyzzyaaat24)then
xyzzyaaad24=xyzzyaaad24-xyzzyaaal24
xyzzyaaay24=0
endif
xyzzyaaae24=xyzzyaaad24
xyzzyaaaz24=xyzzyaaar24-1
do xyzzyaaab24=xyzzyaaar24,xyzzyaaar24+3
xyzzyaaae24=xyzzyaaae24+xyzzyaaao24
xyzzyaaaz24=xyzzyaaaz24+1
if(xyzzyaaaz24==xyzzyaaau24)then
xyzzyaaae24=xyzzyaaae24-xyzzyaaam24
xyzzyaaaz24=0
endif
xyzzyaaaf24=xyzzyaaae24
xyzzyaaba24=xyzzyaaas24-1
do xyzzyaaac24=xyzzyaaas24,xyzzyaaas24+3
xyzzyaaaf24=xyzzyaaaf24+1
xyzzyaaba24=xyzzyaaba24+1
if(xyzzyaaba24==xyzzyaaav24)then
xyzzyaaaf24=xyzzyaaaf24-xyzzyaaao24
xyzzyaaba24=0
endif
xyzzyaaag24=xyzzyaaag24+1
xyzzyaabb24=vector(xyzzyaaag24)
do xyzzyaaak24=1,xyzzyaaaj24
mult(xyzzyaaak24)=mult(xyzzyaaak24)+xyzzyaabb24*bcoef(xyzzyaaai24+xyzz&
&yaaak24,xyzzyaaaf24)
enddo
enddo
enddo
enddo
elseif(xyzzyaaar24+3>=xyzzyaaau24)then
xyzzyaaad24=xyzzyaaah24
xyzzyaaay24=xyzzyaaaq24-1
do xyzzyaaaa24=xyzzyaaaq24,xyzzyaaaq24+3
xyzzyaaad24=xyzzyaaad24+xyzzyaaam24
xyzzyaaay24=xyzzyaaay24+1
if(xyzzyaaay24==xyzzyaaat24)then
xyzzyaaad24=xyzzyaaad24-xyzzyaaal24
xyzzyaaay24=0
endif
xyzzyaaae24=xyzzyaaad24
xyzzyaaaz24=xyzzyaaar24-1
do xyzzyaaab24=xyzzyaaar24,xyzzyaaar24+3
xyzzyaaae24=xyzzyaaae24+xyzzyaaao24
xyzzyaaaz24=xyzzyaaaz24+1
if(xyzzyaaaz24==xyzzyaaau24)then
xyzzyaaae24=xyzzyaaae24-xyzzyaaam24
xyzzyaaaz24=0
endif
do xyzzyaaak24=1,xyzzyaaaj24
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+1)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaae24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+2)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaae24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+3)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaae24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+4)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaae24+4)
enddo
xyzzyaaag24=xyzzyaaag24+4
enddo
enddo
elseif(xyzzyaaaq24+3>=xyzzyaaat24)then
xyzzyaaaw24=xyzzyaaao24+xyzzyaaao24
xyzzyaaax24=xyzzyaaaw24+xyzzyaaao24
xyzzyaaad24=xyzzyaaah24
xyzzyaaay24=xyzzyaaaq24-1
do xyzzyaaaa24=xyzzyaaaq24,xyzzyaaaq24+3
xyzzyaaad24=xyzzyaaad24+xyzzyaaam24
xyzzyaaay24=xyzzyaaay24+1
if(xyzzyaaay24==xyzzyaaat24)then
xyzzyaaad24=xyzzyaaad24-xyzzyaaal24
xyzzyaaay24=0
endif
do xyzzyaaak24=1,xyzzyaaaj24
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+1)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaao24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+2)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaao24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+3)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaao24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+4)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaao24+4)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+5)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaaw24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+6)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaaw24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+7)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaaw24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+8)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaaw24+4)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+9)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaax24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+10)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaax24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+11)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaax24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+12)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaax24+4)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+13)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaap24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+14)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaap24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+15)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaap24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+16)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaap24+4)
enddo
xyzzyaaag24=xyzzyaaag24+16
enddo
else
xyzzyaaaw24=xyzzyaaao24+xyzzyaaao24
xyzzyaaax24=xyzzyaaaw24+xyzzyaaao24
do xyzzyaaad24=xyzzyaaah24+xyzzyaaam24,xyzzyaaah24+xyzzyaaan24,xyzzyaa&
&am24
do xyzzyaaak24=1,xyzzyaaaj24
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+1)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaao24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+2)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaao24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+3)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaao24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+4)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaao24+4)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+5)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaaw24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+6)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaaw24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+7)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaaw24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+8)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaaw24+4)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+9)*bcoef(xyzzya&
&aai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaax24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+10)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaax24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+11)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaax24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+12)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaax24+4)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+13)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaap24+1)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+14)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaap24+2)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+15)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaap24+3)
mult(xyzzyaaak24)=mult(xyzzyaaak24)+vector(xyzzyaaag24+16)*bcoef(xyzzy&
&aaai24+xyzzyaaak24,xyzzyaaad24+xyzzyaaap24+4)
enddo
xyzzyaaag24=xyzzyaaag24+16
enddo
endif
end subroutine xyzzyaaec1
subroutine xyzzyaaed1(ir,xyzzyaaaa1,iorb1,iorb2,norb,bcoef,vector,mult&
&)
implicit none
integer,intent(in) :: ir(3),xyzzyaaaa1(3),iorb1,iorb2,norb
real(dp),intent(in) :: bcoef(norb,*),vector(64)
real(dp),intent(inout) :: mult(1:iorb2-iorb1+1)
integer xyzzyaaaa25,xyzzyaaab25,xyzzyaaac25,xyzzyaaad25,xyzzyaaae25,xy&
&zzyaaaf25,xyzzyaaag25,xyzzyaaah25,xyzzyaaai25,xyzzyaaaj25,xyzzyaaak25&
&,xyzzyaaal25,xyzzyaaam25,xyzzyaaan25,xyzzyaaao25,xyzzyaaap25,xyzzyaaa&
&q25,xyzzyaaar25,xyzzyaaas25,xyzzyaaat25,xyzzyaaau25,xyzzyaaav25,xyzzy&
&aaaw25,xyzzyaaax25,xyzzyaaay25,xyzzyaaaz25,xyzzyaaba25
real(dp) xyzzyaabb25
mult=0.d0
xyzzyaaai25=iorb1-1
xyzzyaaaj25=iorb2-xyzzyaaai25
xyzzyaaaq25=ir(1)
xyzzyaaar25=ir(2)
xyzzyaaas25=ir(3)
xyzzyaaat25=xyzzyaaaa1(1)
xyzzyaaau25=xyzzyaaaa1(2)
xyzzyaaav25=xyzzyaaaa1(3)
xyzzyaaao25=xyzzyaaav25
xyzzyaaap25=xyzzyaaao25+xyzzyaaao25+xyzzyaaao25+xyzzyaaao25
xyzzyaaam25=xyzzyaaau25*xyzzyaaao25
xyzzyaaan25=xyzzyaaam25+xyzzyaaam25+xyzzyaaam25+xyzzyaaam25
xyzzyaaal25=xyzzyaaat25*xyzzyaaam25
xyzzyaaah25=(xyzzyaaaq25-1)*xyzzyaaam25+(xyzzyaaar25-1)*xyzzyaaao25+xy&
&zzyaaas25
xyzzyaaag25=0
if(xyzzyaaas25+3>=xyzzyaaav25)then
xyzzyaaad25=xyzzyaaah25
xyzzyaaay25=xyzzyaaaq25-1
do xyzzyaaaa25=xyzzyaaaq25,xyzzyaaaq25+3
xyzzyaaad25=xyzzyaaad25+xyzzyaaam25
xyzzyaaay25=xyzzyaaay25+1
if(xyzzyaaay25==xyzzyaaat25)then
xyzzyaaad25=xyzzyaaad25-xyzzyaaal25
xyzzyaaay25=0
endif
xyzzyaaae25=xyzzyaaad25
xyzzyaaaz25=xyzzyaaar25-1
do xyzzyaaab25=xyzzyaaar25,xyzzyaaar25+3
xyzzyaaae25=xyzzyaaae25+xyzzyaaao25
xyzzyaaaz25=xyzzyaaaz25+1
if(xyzzyaaaz25==xyzzyaaau25)then
xyzzyaaae25=xyzzyaaae25-xyzzyaaam25
xyzzyaaaz25=0
endif
xyzzyaaaf25=xyzzyaaae25
xyzzyaaba25=xyzzyaaas25-1
do xyzzyaaac25=xyzzyaaas25,xyzzyaaas25+3
xyzzyaaaf25=xyzzyaaaf25+1
xyzzyaaba25=xyzzyaaba25+1
if(xyzzyaaba25==xyzzyaaav25)then
xyzzyaaaf25=xyzzyaaaf25-xyzzyaaao25
xyzzyaaba25=0
endif
xyzzyaaag25=xyzzyaaag25+1
xyzzyaabb25=vector(xyzzyaaag25)
do xyzzyaaak25=1,xyzzyaaaj25
mult(xyzzyaaak25)=mult(xyzzyaaak25)+xyzzyaabb25*bcoef(xyzzyaaai25+xyzz&
&yaaak25,xyzzyaaaf25)
enddo
enddo
enddo
enddo
elseif(xyzzyaaar25+3>=xyzzyaaau25)then
xyzzyaaad25=xyzzyaaah25
xyzzyaaay25=xyzzyaaaq25-1
do xyzzyaaaa25=xyzzyaaaq25,xyzzyaaaq25+3
xyzzyaaad25=xyzzyaaad25+xyzzyaaam25
xyzzyaaay25=xyzzyaaay25+1
if(xyzzyaaay25==xyzzyaaat25)then
xyzzyaaad25=xyzzyaaad25-xyzzyaaal25
xyzzyaaay25=0
endif
xyzzyaaae25=xyzzyaaad25
xyzzyaaaz25=xyzzyaaar25-1
do xyzzyaaab25=xyzzyaaar25,xyzzyaaar25+3
xyzzyaaae25=xyzzyaaae25+xyzzyaaao25
xyzzyaaaz25=xyzzyaaaz25+1
if(xyzzyaaaz25==xyzzyaaau25)then
xyzzyaaae25=xyzzyaaae25-xyzzyaaam25
xyzzyaaaz25=0
endif
do xyzzyaaak25=1,xyzzyaaaj25
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+1)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaae25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+2)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaae25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+3)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaae25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+4)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaae25+4)
enddo
xyzzyaaag25=xyzzyaaag25+4
enddo
enddo
elseif(xyzzyaaaq25+3>=xyzzyaaat25)then
xyzzyaaaw25=xyzzyaaao25+xyzzyaaao25
xyzzyaaax25=xyzzyaaaw25+xyzzyaaao25
xyzzyaaad25=xyzzyaaah25
xyzzyaaay25=xyzzyaaaq25-1
do xyzzyaaaa25=xyzzyaaaq25,xyzzyaaaq25+3
xyzzyaaad25=xyzzyaaad25+xyzzyaaam25
xyzzyaaay25=xyzzyaaay25+1
if(xyzzyaaay25==xyzzyaaat25)then
xyzzyaaad25=xyzzyaaad25-xyzzyaaal25
xyzzyaaay25=0
endif
do xyzzyaaak25=1,xyzzyaaaj25
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+1)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaao25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+2)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaao25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+3)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaao25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+4)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaao25+4)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+5)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaaw25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+6)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaaw25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+7)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaaw25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+8)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaaw25+4)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+9)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaax25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+10)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaax25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+11)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaax25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+12)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaax25+4)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+13)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaap25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+14)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaap25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+15)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaap25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+16)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaap25+4)
enddo
xyzzyaaag25=xyzzyaaag25+16
enddo
else
xyzzyaaaw25=xyzzyaaao25+xyzzyaaao25
xyzzyaaax25=xyzzyaaaw25+xyzzyaaao25
do xyzzyaaad25=xyzzyaaah25+xyzzyaaam25,xyzzyaaah25+xyzzyaaan25,xyzzyaa&
&am25
do xyzzyaaak25=1,xyzzyaaaj25
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+1)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaao25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+2)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaao25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+3)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaao25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+4)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaao25+4)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+5)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaaw25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+6)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaaw25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+7)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaaw25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+8)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaaw25+4)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+9)*bcoef(xyzzya&
&aai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaax25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+10)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaax25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+11)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaax25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+12)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaax25+4)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+13)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaap25+1)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+14)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaap25+2)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+15)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaap25+3)
mult(xyzzyaaak25)=mult(xyzzyaaak25)+vector(xyzzyaaag25+16)*bcoef(xyzzy&
&aaai25+xyzzyaaak25,xyzzyaaad25+xyzzyaaap25+4)
enddo
xyzzyaaag25=xyzzyaaag25+16
enddo
endif
end subroutine xyzzyaaed1
subroutine xyzzyaaee1(ir,xyzzyaaaa1,iorb1,iorb2,norb,bcoef,vector,mult&
&)
implicit none
integer,intent(in) :: ir(3),xyzzyaaaa1(3),iorb1,iorb2,norb
real(dp),intent(in) :: vector(64)
complex(sp),intent(in) :: bcoef(norb,*)
complex(dp),intent(inout) :: mult(1:iorb2-iorb1+1)
integer xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26,xyzzyaaad26,xyzzyaaae26,xy&
&zzyaaaf26,xyzzyaaag26,xyzzyaaah26,xyzzyaaai26,xyzzyaaaj26,xyzzyaaak26&
&,xyzzyaaal26,xyzzyaaam26,xyzzyaaan26,xyzzyaaao26,xyzzyaaap26,xyzzyaaa&
&q26,xyzzyaaar26,xyzzyaaas26,xyzzyaaat26,xyzzyaaau26,xyzzyaaav26,xyzzy&
&aaaw26,xyzzyaaax26,xyzzyaaay26,xyzzyaaaz26,xyzzyaaba26
real(dp) xyzzyaabb26
mult=czero
xyzzyaaai26=iorb1-1
xyzzyaaaj26=iorb2-xyzzyaaai26
xyzzyaaaq26=ir(1)
xyzzyaaar26=ir(2)
xyzzyaaas26=ir(3)
xyzzyaaat26=xyzzyaaaa1(1)
xyzzyaaau26=xyzzyaaaa1(2)
xyzzyaaav26=xyzzyaaaa1(3)
xyzzyaaao26=xyzzyaaav26
xyzzyaaap26=xyzzyaaao26+xyzzyaaao26+xyzzyaaao26+xyzzyaaao26
xyzzyaaam26=xyzzyaaau26*xyzzyaaao26
xyzzyaaan26=xyzzyaaam26+xyzzyaaam26+xyzzyaaam26+xyzzyaaam26
xyzzyaaal26=xyzzyaaat26*xyzzyaaam26
xyzzyaaah26=(xyzzyaaaq26-1)*xyzzyaaam26+(xyzzyaaar26-1)*xyzzyaaao26+xy&
&zzyaaas26
xyzzyaaag26=0
if(xyzzyaaas26+3>=xyzzyaaav26)then
xyzzyaaad26=xyzzyaaah26
xyzzyaaay26=xyzzyaaaq26-1
do xyzzyaaaa26=xyzzyaaaq26,xyzzyaaaq26+3
xyzzyaaad26=xyzzyaaad26+xyzzyaaam26
xyzzyaaay26=xyzzyaaay26+1
if(xyzzyaaay26==xyzzyaaat26)then
xyzzyaaad26=xyzzyaaad26-xyzzyaaal26
xyzzyaaay26=0
endif
xyzzyaaae26=xyzzyaaad26
xyzzyaaaz26=xyzzyaaar26-1
do xyzzyaaab26=xyzzyaaar26,xyzzyaaar26+3
xyzzyaaae26=xyzzyaaae26+xyzzyaaao26
xyzzyaaaz26=xyzzyaaaz26+1
if(xyzzyaaaz26==xyzzyaaau26)then
xyzzyaaae26=xyzzyaaae26-xyzzyaaam26
xyzzyaaaz26=0
endif
xyzzyaaaf26=xyzzyaaae26
xyzzyaaba26=xyzzyaaas26-1
do xyzzyaaac26=xyzzyaaas26,xyzzyaaas26+3
xyzzyaaaf26=xyzzyaaaf26+1
xyzzyaaba26=xyzzyaaba26+1
if(xyzzyaaba26==xyzzyaaav26)then
xyzzyaaaf26=xyzzyaaaf26-xyzzyaaao26
xyzzyaaba26=0
endif
xyzzyaaag26=xyzzyaaag26+1
xyzzyaabb26=vector(xyzzyaaag26)
do xyzzyaaak26=1,xyzzyaaaj26
mult(xyzzyaaak26)=mult(xyzzyaaak26)+xyzzyaabb26*bcoef(xyzzyaaai26+xyzz&
&yaaak26,xyzzyaaaf26)
enddo
enddo
enddo
enddo
elseif(xyzzyaaar26+3>=xyzzyaaau26)then
xyzzyaaad26=xyzzyaaah26
xyzzyaaay26=xyzzyaaaq26-1
do xyzzyaaaa26=xyzzyaaaq26,xyzzyaaaq26+3
xyzzyaaad26=xyzzyaaad26+xyzzyaaam26
xyzzyaaay26=xyzzyaaay26+1
if(xyzzyaaay26==xyzzyaaat26)then
xyzzyaaad26=xyzzyaaad26-xyzzyaaal26
xyzzyaaay26=0
endif
xyzzyaaae26=xyzzyaaad26
xyzzyaaaz26=xyzzyaaar26-1
do xyzzyaaab26=xyzzyaaar26,xyzzyaaar26+3
xyzzyaaae26=xyzzyaaae26+xyzzyaaao26
xyzzyaaaz26=xyzzyaaaz26+1
if(xyzzyaaaz26==xyzzyaaau26)then
xyzzyaaae26=xyzzyaaae26-xyzzyaaam26
xyzzyaaaz26=0
endif
do xyzzyaaak26=1,xyzzyaaaj26
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+1)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaae26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+2)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaae26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+3)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaae26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+4)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaae26+4)
enddo
xyzzyaaag26=xyzzyaaag26+4
enddo
enddo
elseif(xyzzyaaaq26+3>=xyzzyaaat26)then
xyzzyaaaw26=xyzzyaaao26+xyzzyaaao26
xyzzyaaax26=xyzzyaaaw26+xyzzyaaao26
xyzzyaaad26=xyzzyaaah26
xyzzyaaay26=xyzzyaaaq26-1
do xyzzyaaaa26=xyzzyaaaq26,xyzzyaaaq26+3
xyzzyaaad26=xyzzyaaad26+xyzzyaaam26
xyzzyaaay26=xyzzyaaay26+1
if(xyzzyaaay26==xyzzyaaat26)then
xyzzyaaad26=xyzzyaaad26-xyzzyaaal26
xyzzyaaay26=0
endif
do xyzzyaaak26=1,xyzzyaaaj26
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+1)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaao26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+2)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaao26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+3)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaao26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+4)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaao26+4)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+5)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaaw26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+6)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaaw26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+7)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaaw26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+8)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaaw26+4)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+9)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaax26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+10)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaax26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+11)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaax26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+12)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaax26+4)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+13)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaap26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+14)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaap26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+15)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaap26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+16)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaap26+4)
enddo
xyzzyaaag26=xyzzyaaag26+16
enddo
else
xyzzyaaaw26=xyzzyaaao26+xyzzyaaao26
xyzzyaaax26=xyzzyaaaw26+xyzzyaaao26
do xyzzyaaad26=xyzzyaaah26+xyzzyaaam26,xyzzyaaah26+xyzzyaaan26,xyzzyaa&
&am26
do xyzzyaaak26=1,xyzzyaaaj26
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+1)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaao26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+2)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaao26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+3)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaao26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+4)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaao26+4)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+5)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaaw26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+6)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaaw26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+7)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaaw26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+8)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaaw26+4)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+9)*bcoef(xyzzya&
&aai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaax26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+10)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaax26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+11)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaax26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+12)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaax26+4)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+13)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaap26+1)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+14)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaap26+2)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+15)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaap26+3)
mult(xyzzyaaak26)=mult(xyzzyaaak26)+vector(xyzzyaaag26+16)*bcoef(xyzzy&
&aaai26+xyzzyaaak26,xyzzyaaad26+xyzzyaaap26+4)
enddo
xyzzyaaag26=xyzzyaaag26+16
enddo
endif
end subroutine xyzzyaaee1
subroutine xyzzyaaef1(ir,xyzzyaaaa1,iorb1,iorb2,norb,bcoef,vector,mult&
&)
implicit none
integer,intent(in) :: ir(3),xyzzyaaaa1(3),iorb1,iorb2,norb
real(dp),intent(in) :: vector(64)
complex(dp),intent(in) :: bcoef(norb,*)
complex(dp),intent(inout) :: mult(1:iorb2-iorb1+1)
integer xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xyzzyaaad27,xyzzyaaae27,xy&
&zzyaaaf27,xyzzyaaag27,xyzzyaaah27,xyzzyaaai27,xyzzyaaaj27,xyzzyaaak27&
&,xyzzyaaal27,xyzzyaaam27,xyzzyaaan27,xyzzyaaao27,xyzzyaaap27,xyzzyaaa&
&q27,xyzzyaaar27,xyzzyaaas27,xyzzyaaat27,xyzzyaaau27,xyzzyaaav27,xyzzy&
&aaaw27,xyzzyaaax27,xyzzyaaay27,xyzzyaaaz27,xyzzyaaba27
real(dp) xyzzyaabb27
mult=czero
xyzzyaaai27=iorb1-1
xyzzyaaaj27=iorb2-xyzzyaaai27
xyzzyaaaq27=ir(1)
xyzzyaaar27=ir(2)
xyzzyaaas27=ir(3)
xyzzyaaat27=xyzzyaaaa1(1)
xyzzyaaau27=xyzzyaaaa1(2)
xyzzyaaav27=xyzzyaaaa1(3)
xyzzyaaao27=xyzzyaaav27
xyzzyaaap27=xyzzyaaao27+xyzzyaaao27+xyzzyaaao27+xyzzyaaao27
xyzzyaaam27=xyzzyaaau27*xyzzyaaao27
xyzzyaaan27=xyzzyaaam27+xyzzyaaam27+xyzzyaaam27+xyzzyaaam27
xyzzyaaal27=xyzzyaaat27*xyzzyaaam27
xyzzyaaah27=(xyzzyaaaq27-1)*xyzzyaaam27+(xyzzyaaar27-1)*xyzzyaaao27+xy&
&zzyaaas27
xyzzyaaag27=0
if(xyzzyaaas27+3>=xyzzyaaav27)then
xyzzyaaad27=xyzzyaaah27
xyzzyaaay27=xyzzyaaaq27-1
do xyzzyaaaa27=xyzzyaaaq27,xyzzyaaaq27+3
xyzzyaaad27=xyzzyaaad27+xyzzyaaam27
xyzzyaaay27=xyzzyaaay27+1
if(xyzzyaaay27==xyzzyaaat27)then
xyzzyaaad27=xyzzyaaad27-xyzzyaaal27
xyzzyaaay27=0
endif
xyzzyaaae27=xyzzyaaad27
xyzzyaaaz27=xyzzyaaar27-1
do xyzzyaaab27=xyzzyaaar27,xyzzyaaar27+3
xyzzyaaae27=xyzzyaaae27+xyzzyaaao27
xyzzyaaaz27=xyzzyaaaz27+1
if(xyzzyaaaz27==xyzzyaaau27)then
xyzzyaaae27=xyzzyaaae27-xyzzyaaam27
xyzzyaaaz27=0
endif
xyzzyaaaf27=xyzzyaaae27
xyzzyaaba27=xyzzyaaas27-1
do xyzzyaaac27=xyzzyaaas27,xyzzyaaas27+3
xyzzyaaaf27=xyzzyaaaf27+1
xyzzyaaba27=xyzzyaaba27+1
if(xyzzyaaba27==xyzzyaaav27)then
xyzzyaaaf27=xyzzyaaaf27-xyzzyaaao27
xyzzyaaba27=0
endif
xyzzyaaag27=xyzzyaaag27+1
xyzzyaabb27=vector(xyzzyaaag27)
do xyzzyaaak27=1,xyzzyaaaj27
mult(xyzzyaaak27)=mult(xyzzyaaak27)+xyzzyaabb27*bcoef(xyzzyaaai27+xyzz&
&yaaak27,xyzzyaaaf27)
enddo
enddo
enddo
enddo
elseif(xyzzyaaar27+3>=xyzzyaaau27)then
xyzzyaaad27=xyzzyaaah27
xyzzyaaay27=xyzzyaaaq27-1
do xyzzyaaaa27=xyzzyaaaq27,xyzzyaaaq27+3
xyzzyaaad27=xyzzyaaad27+xyzzyaaam27
xyzzyaaay27=xyzzyaaay27+1
if(xyzzyaaay27==xyzzyaaat27)then
xyzzyaaad27=xyzzyaaad27-xyzzyaaal27
xyzzyaaay27=0
endif
xyzzyaaae27=xyzzyaaad27
xyzzyaaaz27=xyzzyaaar27-1
do xyzzyaaab27=xyzzyaaar27,xyzzyaaar27+3
xyzzyaaae27=xyzzyaaae27+xyzzyaaao27
xyzzyaaaz27=xyzzyaaaz27+1
if(xyzzyaaaz27==xyzzyaaau27)then
xyzzyaaae27=xyzzyaaae27-xyzzyaaam27
xyzzyaaaz27=0
endif
do xyzzyaaak27=1,xyzzyaaaj27
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+1)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaae27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+2)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaae27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+3)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaae27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+4)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaae27+4)
enddo
xyzzyaaag27=xyzzyaaag27+4
enddo
enddo
elseif(xyzzyaaaq27+3>=xyzzyaaat27)then
xyzzyaaaw27=xyzzyaaao27+xyzzyaaao27
xyzzyaaax27=xyzzyaaaw27+xyzzyaaao27
xyzzyaaad27=xyzzyaaah27
xyzzyaaay27=xyzzyaaaq27-1
do xyzzyaaaa27=xyzzyaaaq27,xyzzyaaaq27+3
xyzzyaaad27=xyzzyaaad27+xyzzyaaam27
xyzzyaaay27=xyzzyaaay27+1
if(xyzzyaaay27==xyzzyaaat27)then
xyzzyaaad27=xyzzyaaad27-xyzzyaaal27
xyzzyaaay27=0
endif
do xyzzyaaak27=1,xyzzyaaaj27
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+1)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaao27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+2)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaao27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+3)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaao27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+4)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaao27+4)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+5)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaaw27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+6)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaaw27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+7)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaaw27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+8)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaaw27+4)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+9)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaax27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+10)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaax27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+11)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaax27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+12)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaax27+4)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+13)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaap27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+14)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaap27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+15)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaap27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+16)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaap27+4)
enddo
xyzzyaaag27=xyzzyaaag27+16
enddo
else
xyzzyaaaw27=xyzzyaaao27+xyzzyaaao27
xyzzyaaax27=xyzzyaaaw27+xyzzyaaao27
do xyzzyaaad27=xyzzyaaah27+xyzzyaaam27,xyzzyaaah27+xyzzyaaan27,xyzzyaa&
&am27
do xyzzyaaak27=1,xyzzyaaaj27
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+1)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaao27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+2)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaao27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+3)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaao27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+4)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaao27+4)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+5)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaaw27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+6)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaaw27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+7)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaaw27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+8)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaaw27+4)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+9)*bcoef(xyzzya&
&aai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaax27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+10)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaax27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+11)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaax27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+12)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaax27+4)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+13)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaap27+1)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+14)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaap27+2)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+15)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaap27+3)
mult(xyzzyaaak27)=mult(xyzzyaaak27)+vector(xyzzyaaag27+16)*bcoef(xyzzy&
&aaai27+xyzzyaaak27,xyzzyaaad27+xyzzyaaap27+4)
enddo
xyzzyaaag27=xyzzyaaag27+16
enddo
endif
end subroutine xyzzyaaef1
subroutine xyzzyaaeg1(val,fsd,blips_computed,r_frac,f123,fg123,fl123)
implicit none
real(dp),intent(in) :: r_frac(3)
real(dp),intent(inout) :: f123(64),fg123(192),fl123(384)
logical,intent(in) :: val,fsd
logical,intent(inout) :: blips_computed
integer xyzzyaaaa28,xyzzyaaab28,xyzzyaaac28,xyzzyaaad28
real(dp) xyzzyaaae28(4),xyzzyaaaf28(4),xyzzyaaag28(4),xyzzyaaah28(4),x&
&yzzyaaai28(4),xyzzyaaaj28(4),xyzzyaaak28(4),xyzzyaaal28(4),xyzzyaaam2&
&8(4),xyzzyaaan28,xyzzyaaao28,xyzzyaaap28,xyzzyaaaq28,xyzzyaaar28,xyzz&
&yaaas28,xyzzyaaat28,xyzzyaaau28,xyzzyaaav28,xyzzyaaaw28,xyzzyaaax28,x&
&yzzyaaay28,xyzzyaaaz28,xyzzyaaba28,xyzzyaabb28,xyzzyaabc28(3)
if(blips_computed)return
xyzzyaabc28(:)=modulo(r_frac(:)*xyzzyaaab1(:),1.d0)
if(.not.fsd)then
call xyzzyaaeh1(xyzzyaabc28(1),xyzzyaaae28(1))
call xyzzyaaeh1(xyzzyaabc28(2),xyzzyaaaf28(1))
call xyzzyaaeh1(xyzzyaabc28(3),xyzzyaaag28(1))
else
call xyzzyaaei1(xyzzyaabc28(1),xyzzyaaab1(1),xyzzyaaac1(1),xyzzyaaae28&
&(1),xyzzyaaah28(1),xyzzyaaak28(1))
call xyzzyaaei1(xyzzyaabc28(2),xyzzyaaab1(2),xyzzyaaac1(2),xyzzyaaaf28&
&(1),xyzzyaaai28(1),xyzzyaaal28(1))
call xyzzyaaei1(xyzzyaabc28(3),xyzzyaaab1(3),xyzzyaaac1(3),xyzzyaaag28&
&(1),xyzzyaaaj28(1),xyzzyaaam28(1))
endif
if(val.and.fsd)then
xyzzyaaad28=0
do xyzzyaaac28=1,4
xyzzyaaan28=xyzzyaaae28(xyzzyaaac28)
xyzzyaaaq28=xyzzyaaah28(xyzzyaaac28)
xyzzyaaat28=xyzzyaaak28(xyzzyaaac28)
do xyzzyaaab28=1,4
xyzzyaaao28=xyzzyaaaf28(xyzzyaaab28)
xyzzyaaar28=xyzzyaaai28(xyzzyaaab28)
xyzzyaaau28=xyzzyaaal28(xyzzyaaab28)
xyzzyaaaw28  =xyzzyaaan28*xyzzyaaao28
xyzzyaaax28 =xyzzyaaan28*xyzzyaaar28
xyzzyaaay28 =xyzzyaaaq28*xyzzyaaao28
xyzzyaaaz28=xyzzyaaan28*xyzzyaaau28
xyzzyaaba28=xyzzyaaat28*xyzzyaaao28
xyzzyaabb28=xyzzyaaaq28*xyzzyaaar28
do xyzzyaaaa28=1,4
xyzzyaaad28=xyzzyaaad28+1
xyzzyaaap28=xyzzyaaag28(xyzzyaaaa28)
xyzzyaaas28=xyzzyaaaj28(xyzzyaaaa28)
xyzzyaaav28=xyzzyaaam28(xyzzyaaaa28)
f123(xyzzyaaad28)=xyzzyaaaw28*xyzzyaaap28
fg123(xyzzyaaad28)=xyzzyaaay28*xyzzyaaap28
fg123(xyzzyaaad28+64)=xyzzyaaax28*xyzzyaaap28
fg123(xyzzyaaad28+128)=xyzzyaaaw28*xyzzyaaas28
fl123(xyzzyaaad28)=xyzzyaaba28*xyzzyaaap28
fl123(xyzzyaaad28+64)=xyzzyaaaz28*xyzzyaaap28
fl123(xyzzyaaad28+128)=xyzzyaaaw28*xyzzyaaav28
fl123(xyzzyaaad28+192)=xyzzyaabb28*xyzzyaaap28
fl123(xyzzyaaad28+256)=xyzzyaaay28*xyzzyaaas28
fl123(xyzzyaaad28+320)=xyzzyaaax28*xyzzyaaas28
enddo
enddo
enddo
elseif(fsd)then
xyzzyaaad28=0
do xyzzyaaac28=1,4
xyzzyaaan28=xyzzyaaae28(xyzzyaaac28)
xyzzyaaaq28=xyzzyaaah28(xyzzyaaac28)
xyzzyaaat28=xyzzyaaak28(xyzzyaaac28)
do xyzzyaaab28=1,4
xyzzyaaao28=xyzzyaaaf28(xyzzyaaab28)
xyzzyaaar28=xyzzyaaai28(xyzzyaaab28)
xyzzyaaau28=xyzzyaaal28(xyzzyaaab28)
xyzzyaaaw28  =xyzzyaaan28*xyzzyaaao28
xyzzyaaax28 =xyzzyaaan28*xyzzyaaar28
xyzzyaaay28 =xyzzyaaaq28*xyzzyaaao28
xyzzyaaaz28=xyzzyaaan28*xyzzyaaau28
xyzzyaaba28=xyzzyaaat28*xyzzyaaao28
xyzzyaabb28=xyzzyaaaq28*xyzzyaaar28
do xyzzyaaaa28=1,4
xyzzyaaad28=xyzzyaaad28+1
xyzzyaaap28=xyzzyaaag28(xyzzyaaaa28)
xyzzyaaas28=xyzzyaaaj28(xyzzyaaaa28)
xyzzyaaav28=xyzzyaaam28(xyzzyaaaa28)
fg123(xyzzyaaad28)=xyzzyaaay28*xyzzyaaap28
fg123(xyzzyaaad28+64)=xyzzyaaax28*xyzzyaaap28
fg123(xyzzyaaad28+128)=xyzzyaaaw28*xyzzyaaas28
fl123(xyzzyaaad28)=xyzzyaaba28*xyzzyaaap28
fl123(xyzzyaaad28+64)=xyzzyaaaz28*xyzzyaaap28
fl123(xyzzyaaad28+128)=xyzzyaaaw28*xyzzyaaav28
fl123(xyzzyaaad28+192)=xyzzyaabb28*xyzzyaaap28
fl123(xyzzyaaad28+256)=xyzzyaaay28*xyzzyaaas28
fl123(xyzzyaaad28+320)=xyzzyaaax28*xyzzyaaas28
enddo
enddo
enddo
elseif(val)then
xyzzyaaad28=0
do xyzzyaaac28=1,4
xyzzyaaan28=xyzzyaaae28(xyzzyaaac28)
do xyzzyaaab28=1,4
xyzzyaaaw28=xyzzyaaan28*xyzzyaaaf28(xyzzyaaab28)
do xyzzyaaaa28=1,4
xyzzyaaad28=xyzzyaaad28+1
f123(xyzzyaaad28)=xyzzyaaaw28*xyzzyaaag28(xyzzyaaaa28)
enddo
enddo
enddo
endif
blips_computed=.true.
end subroutine xyzzyaaeg1
subroutine xyzzyaaeh1(x,f)
implicit none
real(dp),intent(in) :: x
real(dp),intent(out) :: f(4)
real(dp) xyzzyaaaa29
xyzzyaaaa29=x+1.d0
f(1)=2.d0+xyzzyaaaa29*(-3.d0+xyzzyaaaa29*(1.5d0-0.25d0*xyzzyaaaa29))
xyzzyaaaa29=x
f(2)=1.d0+xyzzyaaaa29*xyzzyaaaa29*(-1.5d0+0.75d0*xyzzyaaaa29)
xyzzyaaaa29=x-1.d0
f(3)=1.d0+xyzzyaaaa29*xyzzyaaaa29*(-1.5d0-0.75d0*xyzzyaaaa29)
xyzzyaaaa29=x-2.d0
f(4)=2.d0+xyzzyaaaa29*(3.d0+xyzzyaaaa29*(1.5d0+0.25d0*xyzzyaaaa29))
end subroutine xyzzyaaeh1
subroutine xyzzyaaei1(x,prefac_d,prefac_d2,f,df,d2f)
implicit none
real(dp),intent(in) :: x,prefac_d,prefac_d2
real(dp),intent(out) :: f(4),df(4),d2f(4)
real(dp) xyzzyaaaa30
xyzzyaaaa30=x+1.d0
f(1)=2.d0+xyzzyaaaa30*(-3.d0+xyzzyaaaa30*(1.5d0-0.25d0*xyzzyaaaa30))
df(1)=(-3.d0+xyzzyaaaa30*(3.d0-0.75d0*xyzzyaaaa30))*prefac_d
d2f(1)=(3.d0-1.5d0*xyzzyaaaa30)*prefac_d2
xyzzyaaaa30=x
f(2)=1.d0+xyzzyaaaa30*xyzzyaaaa30*(-1.5d0+0.75d0*xyzzyaaaa30)
df(2)=xyzzyaaaa30*(-3.d0+2.25d0*xyzzyaaaa30)*prefac_d
d2f(2)=(-3.d0+4.5d0*xyzzyaaaa30)*prefac_d2
xyzzyaaaa30=x-1.d0
f(3)=1.d0+xyzzyaaaa30*xyzzyaaaa30*(-1.5d0-0.75d0*xyzzyaaaa30)
df(3)=xyzzyaaaa30*(-3.d0-2.25d0*xyzzyaaaa30)*prefac_d
d2f(3)=(-3.d0-4.5d0*xyzzyaaaa30)*prefac_d2
xyzzyaaaa30=x-2.d0
f(4)=2.d0+xyzzyaaaa30*(3.d0+xyzzyaaaa30*(1.5d0+0.25d0*xyzzyaaaa30))
df(4)=(3.d0+xyzzyaaaa30*(3.d0+0.75d0*xyzzyaaaa30))*prefac_d
d2f(4)=(3.d0+1.5d0*xyzzyaaaa30)*prefac_d2
end subroutine xyzzyaaei1
subroutine xyzzyaaej1(iorb,norb,ik,sum_orbs_real2,sum_orbs_im2)
use slaarnabt, only : inverse3
implicit none
integer,intent(in) :: iorb,norb,ik
real(dp),intent(out) :: sum_orbs_real2,sum_orbs_im2
integer xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31,xyzzyaaad31(3),xyzzyaaae31&
&(3)
real(dp) xyzzyaaaf31(3),xyzzyaaag31(3),xyzzyaaah31(3),xyzzyaaai31(3,3)&
&,xyzzyaaaj31(3,3),xyzzyaaak31,xyzzyaaal31(64),xyzzyaaam31(64,3),xyzzy&
&aaan31(64,6)
complex(dp) xyzzyaaao31,xyzzyaaap31(1)
logical xyzzyaaaq31,xyzzyaaar31
xyzzyaaai31(1,:)=pa1(:)
xyzzyaaai31(2,:)=pa2(:)
xyzzyaaai31(3,:)=pa3(:)
xyzzyaaaj31(:,:)=inverse3(xyzzyaaai31(:,:),determinant=xyzzyaaak31)
sum_orbs_real2=0.d0
sum_orbs_im2=0.d0
do xyzzyaaaa31=-3,3
xyzzyaaaf31(1)=dble(xyzzyaaaa31)*12.3456d0-0.01234d0
do xyzzyaaab31=-3,3
xyzzyaaaf31(2)=dble(xyzzyaaab31)*12.3456d0-0.04567d0
do xyzzyaaac31=-3,3
xyzzyaaaf31(3)=dble(xyzzyaaac31)*12.3456d0-0.08910d0
xyzzyaaah31=modulo(matmul(xyzzyaaaf31,xyzzyaaaj31),1.d0)
call xyzzyaaek1(xyzzyaaah31,xyzzyaaad31(1),xyzzyaaad31(2),xyzzyaaad31(&
&3),xyzzyaaaq31)
if(.not.xyzzyaaaq31)cycle
xyzzyaaae31(:)=xyzzyaaad31(:)-1
where(xyzzyaaae31<0)xyzzyaaae31=xyzzyaaae31+xyzzyaaaa1
xyzzyaaag31=xyzzyaaah31(1)*pa1(1:3)+xyzzyaaah31(2)*pa2(1:3)+xyzzyaaah3&
&1(3)*pa3(1:3)
xyzzyaaar31=.false.
call xyzzyaaeg1(.true.,.false.,xyzzyaaar31,xyzzyaaah31,xyzzyaaal31,xyz&
&zyaaam31,xyzzyaaan31)
if(single_precision_blips)then
call xyzzyaaee1(xyzzyaaae31,xyzzyaaaa1,iorb,iorb,norb,xyzzyaaaq1,xyzzy&
&aaal31,xyzzyaaap31)
else
call xyzzyaaef1(xyzzyaaae31,xyzzyaaaa1,iorb,iorb,norb,xyzzyaaar1,xyzzy&
&aaal31,xyzzyaaap31)
endif
xyzzyaaao31=exp(cmplx(0.d0,1.d0,dp)*dot_product(xyzzyaaay1(:,ik),xyzzy&
&aaag31))*xyzzyaaap31(1)
sum_orbs_real2=sum_orbs_real2+dble(xyzzyaaao31)**2
sum_orbs_im2=sum_orbs_im2+aimag(xyzzyaaao31)**2
enddo
enddo
enddo
end subroutine xyzzyaaej1
subroutine xyzzyaaek1(r,ix,iy,iz,in_grid)
implicit none
real(dp),intent(in) :: r(3)
integer,intent(out) :: ix,iy,iz
logical,intent(out) :: in_grid
if(periodicity==3)then
ix=modulo(floor(r(1)*xyzzyaaab1(1)),xyzzyaaaa1(1))
iy=modulo(floor(r(2)*xyzzyaaab1(2)),xyzzyaaaa1(2))
iz=modulo(floor(r(3)*xyzzyaaab1(3)),xyzzyaaaa1(3))
in_grid=.true.
elseif(periodicity==2)then
iz=floor(r(3)*xyzzyaaab1(3))
if(iz<0.or.iz>=xyzzyaaaa1(3))then
in_grid=.false.
else
in_grid=.true.
ix=modulo(floor(r(1)*xyzzyaaab1(1)),xyzzyaaaa1(1))
iy=modulo(floor(r(2)*xyzzyaaab1(2)),xyzzyaaaa1(2))
endif
elseif(periodicity==1)then
iy=floor(r(2)*xyzzyaaab1(2))
iz=floor(r(3)*xyzzyaaab1(3))
if(iy<0.or.iy>=xyzzyaaaa1(2).or.iz<0.or.iz>=xyzzyaaaa1(3))then
in_grid=.false.
else
in_grid=.true.
ix=modulo(floor(r(1)*xyzzyaaab1(1)),xyzzyaaaa1(1))
endif
else
ix=floor(r(1)*xyzzyaaab1(1))
iy=floor(r(2)*xyzzyaaab1(2))
iz=floor(r(3)*xyzzyaaab1(3))
in_grid=(ix>=0.and.ix<xyzzyaaaa1(1).and.iy>=0.and.iy<xyzzyaaaa1(2).and&
&.iz>=0.and.iz<xyzzyaaaa1(3))
endif
end subroutine xyzzyaaek1
subroutine get_bwfdet_orbmap(row_offset,norb,orbmap)
implicit none
integer,intent(inout) :: row_offset(nspin),norb,orbmap(nemax,nspin,nde&
&t)
integer xyzzyaaaa33
do xyzzyaaaa33=1,nspin
orbmap(row_offset(xyzzyaaaa33)+1:row_offset(xyzzyaaaa33)+nuc_nele(xyzz&
&yaaaa33),xyzzyaaaa33,:)=norb+xyzzyaacm1(1:nuc_nele(xyzzyaaaa33),xyzzy&
&aaaa33,:)
row_offset(xyzzyaaaa33)=row_offset(xyzzyaaaa33)+nuc_nele(xyzzyaaaa33)
enddo
norb=norb+xyzzyaaci1
end subroutine get_bwfdet_orbmap
subroutine get_bwfdet_ndesc(ndesc_int,ndesc_dp)
implicit none
integer,intent(inout) :: ndesc_int,ndesc_dp
ndesc_int=xyzzyaadf1
ndesc_dp=xyzzyaadg1
end subroutine get_bwfdet_ndesc
subroutine get_bwfdet_orbdesc(norb,ndesc_int,ndesc_dp,orbdesc_int,orbd&
&esc_dp)
implicit none
integer,intent(in) :: norb,ndesc_int,ndesc_dp
integer,intent(inout) :: orbdesc_int(ndesc_int,norb)
real(dp),intent(inout) :: orbdesc_dp(ndesc_dp,norb)
integer xyzzyaaaa35
if(xyzzyaabm1)then
do xyzzyaaaa35=1,xyzzyaacj1
orbdesc_int(1:xyzzyaadf1,xyzzyaaaa35)=xyzzyaaaa35
orbdesc_dp(1:xyzzyaadg1,xyzzyaaaa35)=0
enddo
endif
if(xyzzyaabn1)then
orbdesc_int(1:xyzzyaadf1,xyzzyaacj1+1:xyzzyaaci1)=xyzzyaacj1+xyzzyaadh&
&1(1:xyzzyaadf1,1:xyzzyaack1)
orbdesc_dp(1:xyzzyaadg1,xyzzyaacj1+1:xyzzyaaci1)=xyzzyaadi1(1:xyzzyaad&
&g1,1:xyzzyaack1)
endif
end subroutine get_bwfdet_orbdesc
real(dp) function get_bwfdet_rmax()
implicit none
get_bwfdet_rmax=maxval(xyzzyaabe1(:))
end function get_bwfdet_rmax
end module slaarnaae

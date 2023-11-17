module slaarnaau
use slaarnaag
use dsp
use slaarnaat
use slaarnabg
use slaarnabp
use parallel
use store
use file_utils,  only : open_units,skip
use format_utils,only : wout,i2s,r2s,display_param,wordwrap
use slaarnabt,   only : quicksort
use slaarnaca,       only : zion,have_ae,is_ae
use run_control, only : errstop,errwarn,errwarn_silent,errstop_master,&
&timer,check_alloc
use shalloc,     only : shallocate,deshallocate,shallocate_gauss,shall&
&ocate_barrier,am_smpmaster,need_shm,nnpsmp
implicit none
private
public readgw,gwfdet_setup,gaussian_orb_eval,cusp_setup,get_gaussian_r&
&max,get_gaussian_orbmap,get_gaussian_orbdesc,get_gaussian_ndesc,desha&
&lloc_gauss_shm
public read_gwfmolorb,write_gwfmolorb,setup_gwfmolorb_params,finish_gw&
&fmolorb_params,get_gwfmolorb_params,put_gwfmolorb_params
integer xyzzyaaaa1
character(80) param_title
integer xyzzyaaab1
integer xyzzyaaac1
integer,allocatable :: xyzzyaaad1(:)
integer xyzzyaaae1
integer xyzzyaaaf1
integer xyzzyaaag1
integer xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1
integer,allocatable :: xyzzyaaak1(:),xyzzyaaal1(:,:),xyzzyaaam1(:,:),x&
&yzzyaaan1(:,:),xyzzyaaao1(:,:)
integer xyzzyaaap1,xyzzyaaaq1
integer,allocatable :: xyzzyaaar1(:,:,:)
logical,allocatable :: xyzzyaaas1(:,:),xyzzyaaat1(:,:)
integer,allocatable :: xyzzyaaau1(:),xyzzyaaav1(:),xyzzyaaaw1(:)
integer,allocatable :: xyzzyaaax1(:,:,:,:)
real(dp),pointer :: xyzzyaaay1(:)=>null()
real(dp),allocatable :: xyzzyaaaz1(:,:)
real(dp),allocatable :: xyzzyaaba1(:,:,:)
real(dp),allocatable :: xyzzyaabb1(:)
real(dp),allocatable :: xyzzyaabc1(:)
real(dp),allocatable :: xyzzyaabd1(:)
integer,allocatable :: xyzzyaabe1(:),xyzzyaabf1(:),xyzzyaabg1(:)
integer xyzzyaabh1,xyzzyaabi1,xyzzyaabj1
real(dp),allocatable :: xyzzyaabk1(:)
integer,allocatable ::  xyzzyaabl1(:)
logical xyzzyaabm1,xyzzyaabn1,xyzzyaabo1
integer xyzzyaabp1
integer,parameter :: xyzzyaabq1=5
real(dp) xyzzyaabr1
real(dp),dimension(:,:),allocatable :: xyzzyaabs1,xyzzyaabt1,xyzzyaabu&
&1,xyzzyaabv1,xyzzyaabw1,xyzzyaabx1,xyzzyaaby1,xyzzyaabz1,xyzzyaaca1,x&
&yzzyaacb1,xyzzyaacc1
real(dp),dimension(:,:,:,:),allocatable :: xyzzyaacd1,xyzzyaace1,xyzzy&
&aacf1
real(dp),dimension(:,:,:,:,:),allocatable :: xyzzyaacg1
logical xyzzyaach1
integer xyzzyaaci1,xyzzyaacj1
integer,allocatable :: xyzzyaack1(:,:)
real(dp),allocatable :: xyzzyaacl1(:,:)
contains
subroutine readgw(eionion)
implicit none
real(dp),intent(out) :: eionion
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,io,xyzzyaaaj2,xyzzyaaak2,xyzzyaa&
&al2,xyzzyaaam2,xyzzyaaan2,xyzzyaaao2
integer,dimension(12) :: xyzzyaaap2,xyzzyaaaq2
integer,allocatable :: xyzzyaaar2(:)
logical xyzzyaaas2,xyzzyaaat2,xyzzyaaau2,xyzzyaaav2,xyzzyaaaw2(0:399)
logical,allocatable :: xyzzyaaax2(:),xyzzyaaay2(:)
character(2) shell_name(12)
character(6),parameter :: xyzzyaaaz2="(8i10)"
character(11) eigcheck
character(13),parameter :: xyzzyaaba2="(3(1pe20.13))",xyzzyaabb2="(4(1&
&pe20.13))"
character(20) psp_filename,ppchar
character(80) title,method,functional,tmpr
data shell_name/' s','sp',' p',' d',' f',' g',' h',' i','k','no','2s',&
&'3p'/
data xyzzyaaap2 /1,4,3,5,7,9,11,13,15,1,1,3/
data xyzzyaaaq2/1,4,3,6,10,15,21,28,36,1,1,3/
if(need_shm)then
xyzzyaaaa2=1
if(am_smpmaster)xyzzyaaaa2=0
call mpi_comm_split(mpi_comm_world,xyzzyaaaa2,my_node,xyzzyaaao2,ierro&
&r)
call checkmpi(ierror,'comm-splitting in readgw')
else
xyzzyaaao2=mpi_comm_world
endif
if(am_master)then
call open_units(io,xyzzyaaak2)
if(xyzzyaaak2/=0)call errstop('READGW','Unable to find free io unit')
open(io,file='gwfn.data',status='old',form='formatted',action='read',e&
&rr=1)
read(io,'(a)',err=2,end=2)title
call skip(io,4)
read(io,'(a)',err=2,end=2)code
call skip(io,1)
read(io,'(a)',err=2,end=2)method
call skip(io,1)
read(io,'(a)',err=2,end=2)functional
call skip(io,1)
read(io,*,err=2,end=2)periodicity
call skip(io,1)
read(io,*,err=2,end=2)xyzzyaaat2
call skip(io,1)
read(io,*,err=2,end=2)eionion
call skip(io,1)
read(io,*,err=2,end=2)xyzzyaaab1
call skip(io,4)
if(periodicity<0.or.periodicity>3)call errstop('READGW','Invalid perio&
&dicity in gwfn.data.')
if(isperiodic.and.periodicity==0) call errstop('READGW','Periodicity i&
&n gwfn.data and PERIODIC input keyword not compatible.')
if(.not.isperiodic.and.periodicity>0) call errstop('READGW','Periodici&
&ty in gwfn.data and PERIODIC input keyword not compatible.')
tmpr=adjustl(code)
if(tmpr(1:8)/='CRYSTAL9')bohr_to_angstrom=0.529177249d0
xyzzyaach1=.not.xyzzyaaat2
if(nspin==3)then
if(xyzzyaach1)then
spin_dep_in=1
else
spin_dep_in=2
endif
elseif(nspin==2.and.electron_system)then
if(xyzzyaach1)then
spin_dep_in=0
else
spin_dep_in=1
endif
else
call errstop('READGW','Need to generalize specification of spin-depend&
&ence in gwfn.data files.')
endif
xyzzyaaai1=no_ssingles(spin_dep_in)
call wout()
call wout('Title : '//trim(title))
call wout()
call wout('Generating code                           :  '//trim(adjust&
&l(code)))
call wout('Method                                    :  '//trim(adjust&
&l(method)))
call wout('DFT functional                            :  '//trim(adjust&
&l(functional)))
select case (periodicity)
case(0)
call wout('Periodicity                               :  '//trim(i2s(pe&
&riodicity))//' (molecule)')
case(1)
call wout('Periodicity                               :  '//trim(i2s(pe&
&riodicity))//' (polymer)')
case(2)
call wout('Periodicity                               :  '//trim(i2s(pe&
&riodicity))//' (slab)')
case(3)
call wout('Periodicity                               :  '//trim(i2s(pe&
&riodicity))//' (solid)')
case default
call wout('Periodicity                               :  '//trim(i2s(pe&
&riodicity)))
call errstop('READGW','Invalid periodicity. Quitting.')
end select
if(xyzzyaach1)then
call wout('Spin restricted?                          :  Yes')
else
call wout('Spin restricted?                          :  No')
endif
call wout()
if(electron_system)then
if(isperiodic)then
call wout('Electrons per cell                        :  '//trim(i2s(su&
&m(nuc_nele)/npcells)))
else
call wout('Total number of electrons                 :  '//trim(i2s(su&
&m(nuc_nele))))
endif
else
if(isperiodic)then
call wout('Particles per cell                        :  '//trim(i2s(su&
&m(nuc_nele)/npcells)))
else
call wout('Total number of particles                 :  '//trim(i2s(su&
&m(nuc_nele))))
endif
endif
read(io,*,err=2,end=2)xyzzyaaam2
call skip(io,1)
allocate(atno(xyzzyaaam2),valence_charge(xyzzyaaam2),basis(3,xyzzyaaam&
&2),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','')
read(io,xyzzyaaba2,err=2,end=2)((basis(xyzzyaaaa2,xyzzyaaab2),xyzzyaaa&
&a2=1,3),xyzzyaaab2=1,xyzzyaaam2)
call skip(io,1)
read(io,xyzzyaaaz2,err=2,end=2)(atno(xyzzyaaaa2),xyzzyaaaa2=1,xyzzyaaa&
&m2)
call skip(io,1)
read(io,xyzzyaabb2,err=2,end=2)(valence_charge(xyzzyaaaa2),xyzzyaaaa2=&
&1,xyzzyaaam2)
call skip(io,1)
if(periodicity/=0)then
read(io,xyzzyaaba2,err=2,end=2)(pa1(xyzzyaaaa2),xyzzyaaaa2=1,3)
read(io,xyzzyaaba2,err=2,end=2)(pa2(xyzzyaaaa2),xyzzyaaaa2=1,3)
read(io,xyzzyaaba2,err=2,end=2)(pa3(xyzzyaaaa2),xyzzyaaaa2=1,3)
pamat(1,1:3)=pa1(1:3)
pamat(2,1:3)=pa2(1:3)
pamat(3,1:3)=pa3(1:3)
call skip(io,4)
read(io,*,err=2,end=2)num_k
call skip(io,1)
read(io,*,err=2,end=2)num_real_k
call skip(io,1)
num_complex_k=num_k-num_real_k
num_real_k_plus_1=num_real_k+1
complex_states=.false.
if(num_complex_k>0)complex_states=.true.
num_real_k_gt_1=num_real_k>1
num_real_k_odd=mod(num_real_k,2)>0
allocate(kvec(3,num_k),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','')
read(io,xyzzyaaba2,err=2,end=2)((kvec(xyzzyaaaa2,xyzzyaaab2),xyzzyaaaa&
&2=1,3),xyzzyaaab2=1,num_k)
call skip(io,4)
call wout()
call wout('K space net information')
call wout()
call wout('No. k points                              :  '//trim(i2s(nu&
&m_k)))
call wout('No. real k points                         :  '//trim(i2s(nu&
&m_real_k)))
call wout('No. complex k points                      :  '//trim(i2s(nu&
&m_complex_k)))
call wout('k point coords')
do xyzzyaaab2=1,num_k
if(xyzzyaaab2<10)then
write(tmpr,'(i1,3x,3(1pe20.13))')xyzzyaaab2,(kvec(xyzzyaaaa2,xyzzyaaab&
&2),xyzzyaaaa2=1,3)
elseif(xyzzyaaab2<100)then
write(tmpr,'(i2,2x,3(1pe20.13))')xyzzyaaab2,(kvec(xyzzyaaaa2,xyzzyaaab&
&2),xyzzyaaaa2=1,3)
elseif(xyzzyaaab2<1000)then
write(tmpr,'(i3,1x,3(1pe20.13))')xyzzyaaab2,(kvec(xyzzyaaaa2,xyzzyaaab&
&2),xyzzyaaaa2=1,3)
else
write(tmpr,'(i4,3(1pe20.13))')xyzzyaaab2,(kvec(xyzzyaaaa2,xyzzyaaab2),&
&xyzzyaaaa2=1,3)
endif
call wout(tmpr)
enddo
else
num_k=1
num_real_k=1
num_complex_k=0
num_real_k_plus_1=2
num_real_k_odd=.true.
num_real_k_gt_1=.false.
allocate(kvec(3,1),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','')
kvec(1,1)=0.d0
kvec(2,1)=0.d0
kvec(3,1)=0.d0
pamat=0.d0
do xyzzyaaaa2=1,3
pamat(xyzzyaaaa2,xyzzyaaaa2)=500.d0
enddo
pa1(1:3)=pamat(1,1:3)
pa2(1:3)=pamat(2,1:3)
pa3(1:3)=pamat(3,1:3)
complex_states=.false.
call skip(io,3)
endif
if(num_real_k+2*num_complex_k/=npcells)call errstop('READGW','Number o&
&f k points incompatible with NCELL values in input')
read(io,*,err=2,end=2)num_centres
call skip(io,1)
read(io,*,err=2,end=2)num_shells
call skip(io,1)
read(io,*,err=2,end=2)num_ao
call skip(io,1)
read(io,*,err=2,end=2)xyzzyaaac1
call skip(io,1)
read(io,*,err=2,end=2)xyzzyaaan2
call skip(io,1)
if(num_centres/=xyzzyaaam2)then
call wordwrap('Number of atoms not equal to number of Gaussian centres&
& in gwfn.data. If you wish to incorporate floating Gaussian functions&
& not attached to a nuclear position, then define additional "ghost nu&
&clei" in the gwfn.data list of atoms, flagged by having charge zero.'&
&)
call errstop('READGW','Quitting.')
endif
call wout()
if(isperiodic)then
call wout('Basis set information (numbers per primitive cell)')
else
call wout('Basis set information')
endif
call wout()
call wout('Number of atoms                           :  '//trim(i2s(xy&
&zzyaaam2)))
call wout('Number of shells                          :  '//trim(i2s(nu&
&m_shells)))
call wout('Number of basis fns                       :  '//trim(i2s(nu&
&m_ao)))
call wout('Number of primitives                      :  '//trim(i2s(xy&
&zzyaaac1)))
call wout('Highest ang momentum                      :  '//trim(i2s(xy&
&zzyaaan2)))
allocate(shell_am(num_shells),numao_in_shell(num_shells),first_shell(n&
&um_centres+1),first_ao(num_shells+1),primitive(num_shells+1),exponent&
&(xyzzyaaac1),min_exponent(num_shells),c_prim(xyzzyaaac1),xyzzyaaaz1(3&
&,num_shells),xyzzyaaad1(num_shells),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','')
read(io,xyzzyaaaz2,err=2,end=2)(shell_am(xyzzyaaaa2),xyzzyaaaa2=1,num_&
&shells)
call skip(io,1)
read(io,xyzzyaaaz2,err=2,end=2)(xyzzyaaad1(xyzzyaaaa2),xyzzyaaaa2=1,nu&
&m_shells)
call skip(io,1)
read(io,xyzzyaaaz2,err=2,end=2)(first_shell(xyzzyaaaa2),xyzzyaaaa2=1,n&
&um_centres+1)
call skip(io,1)
read(io,xyzzyaabb2,err=2,end=2)(exponent(xyzzyaaaa2),xyzzyaaaa2=1,xyzz&
&yaaac1)
call skip(io,1)
read(io,xyzzyaabb2,err=2,end=2)(c_prim(xyzzyaaaa2),xyzzyaaaa2=1,xyzzya&
&aac1)
call skip(io,1)
if(any(shell_am(:)==2))then
allocate(c_prim2(xyzzyaaac1),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','c_prim2')
read(io,xyzzyaabb2,err=2,end=2)(c_prim2(xyzzyaaaa2),xyzzyaaaa2=1,xyzzy&
&aaac1)
call skip(io,1)
endif
read(io,xyzzyaaba2,err=2,end=2)((xyzzyaaaz1(xyzzyaaaa2,xyzzyaaab2),xyz&
&zyaaaa2=1,3),xyzzyaaab2=1,num_shells)
if(any(shell_am(:)<0))call errstop('READGW','Cartesian Gaussians not i&
&mplemented in current version of CASINO. Please convert to harmonic f&
&orm.')
tmpr=adjustl(code)
if(any(shell_am(:)==0).or.any(abs(shell_am(:))==7).or.any(abs(shell_am&
&(:))==8).or.any(abs(shell_am(:))==9).or.any(abs(shell_am(:))==10).or.&
&any(abs(shell_am(:))>12))call errstop('READGW','Invalid shell type sp&
&ecified in basis (only s,sp,p,d,f,g,2s,3p allowed.)')
if(any(abs(shell_am(:))>5).and.use_backflow)call errstop('READGW','Can&
&not use backflow with Gaussian basis functions of higher angular mome&
&ntum than f, as the required mixed second derivatives are not impleme&
&nted.')
if(tmpr(1:6)=='GAMESS')then
if(any(shell_am(:)==5).or.any(abs(shell_am(:))==6))call errstop('READG&
&W','Invalid shell type specified in basis (s,sp,p,d,2s,3p only in GAM&
&ESS).')
endif
if(tmpr(1:9)=='CRYSTAL95')then
if(any(shell_am(:)==0).or.any(abs(shell_am(:))>4))call errstop('READGW&
&','Invalid shell type specified in basis (s,sp,p,d only in CRYSTAL95)&
&.')
elseif(tmpr(1:9)=='CRYSTAL98')then
if(any(shell_am(:)==0).or.any(abs(shell_am(:))>4))call errstop('READGW&
&','Invalid shell type specified in basis (s,sp,p,d only in CRYSTAL98)&
&.')
elseif(tmpr(1:7)=='CRYSTAL')then
if(any(shell_am(:)==0).or.any(abs(shell_am(:))>5))call errstop('READGW&
&','Invalid shell type specified in basis (s,sp,p,d,f only in modern v&
&ersions of CRYSTAL.)')
endif
do xyzzyaaaa2=1,num_shells
xyzzyaaab2=shell_am(xyzzyaaaa2)
if(xyzzyaaab2>0)then
numao_in_shell(xyzzyaaaa2)=xyzzyaaap2(xyzzyaaab2)
else
if(xyzzyaaab2==-1.or.xyzzyaaab2==-2)shell_am(xyzzyaaaa2)=abs(shell_am(&
&xyzzyaaaa2))
numao_in_shell(xyzzyaaaa2)=xyzzyaaaq2(abs(xyzzyaaab2))
endif
enddo
first_ao(1)=1
primitive(1)=1
do xyzzyaaaa2=1,num_shells
first_ao(xyzzyaaaa2+1)=first_ao(xyzzyaaaa2)+numao_in_shell(xyzzyaaaa2)
primitive(xyzzyaaaa2+1)=primitive(xyzzyaaaa2)+xyzzyaaad1(xyzzyaaaa2)
enddo
xyzzyaaac2=1
do xyzzyaaaa2=1,num_shells
min_exponent(xyzzyaaaa2)=minval(exponent(xyzzyaaac2:xyzzyaaac2+xyzzyaa&
&ad1(xyzzyaaaa2)-1))
xyzzyaaac2=xyzzyaaac2+xyzzyaaad1(xyzzyaaaa2)
enddo
if(any(atno(:)>1000))call errstop('READGW','Flagging different pseudop&
&otentials by adding multiples of 1000 to the atomic number is not sup&
&ported if ATOM_BASIS_TYPE=gaussian.')
xyzzyaaaw2=.false.
call wout()
call wout('Gaussian basis sets:')
call wout(repeat('-',66))
call wout('  Atom  x(au)  y(au)  z(au)    No. Type  Exponent    Coeffi&
&cients')
call wout(repeat('-',66))
do xyzzyaaad2=1,xyzzyaaam2
xyzzyaaaa2=mod(atno(xyzzyaaad2),1000)
xyzzyaaab2=mod(xyzzyaaaa2,100)
write(tmpr,'(i3,1x,a,3f9.3)')xyzzyaaad2,periodic_table(xyzzyaaab2),(ba&
&sis(xyzzyaaac2,xyzzyaaad2),xyzzyaaac2=1,3)
call wout(tmpr)
if(xyzzyaaaw2(xyzzyaaaa2))cycle
xyzzyaaaw2(xyzzyaaaa2)=.true.
do xyzzyaaaa2=first_shell(xyzzyaaad2),first_shell(xyzzyaaad2+1)-1
xyzzyaaac2=shell_am(xyzzyaaaa2)
xyzzyaaab2=first_ao(xyzzyaaaa2)
if(xyzzyaaac2==2)then
write(tmpr,'(30x,i4,a)')xyzzyaaab2,shell_name(xyzzyaaac2)
call wout(tmpr)
do xyzzyaaab2=primitive(xyzzyaaaa2),primitive(xyzzyaaaa2+1)-1
write(tmpr,'(39x,1(1pe10.3),2(1pe10.3))')exponent(xyzzyaaab2),c_prim(x&
&yzzyaaab2),c_prim2(xyzzyaaab2)
call wout(tmpr)
enddo
else
write(tmpr,'(25x,i4,''-'',i4,1x,a)')xyzzyaaab2,first_ao(xyzzyaaaa2+1)-&
&1,shell_name(xyzzyaaac2)
call wout(tmpr)
do xyzzyaaab2=primitive(xyzzyaaaa2),primitive(xyzzyaaaa2+1)-1
write(tmpr,'(39x,1(1pe10.3),4x,1(1pe10.3))')exponent(xyzzyaaab2),c_pri&
&m(xyzzyaaab2)
call wout(tmpr)
enddo
endif
enddo
enddo
call wout(repeat('-',66))
allocate(xyzzyaaax2(xyzzyaaam2),xyzzyaaay2(xyzzyaaam2),xyzzyaaar2(xyzz&
&yaaam2),iontype_prim(xyzzyaaam2),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','is_pp_from_atno, etc..')
xyzzyaaar2=0
iontype_prim=0
xyzzyaaax2(:)=atno(:)>200
do xyzzyaaad2=1,xyzzyaaam2
atno(xyzzyaaad2)=mod(atno(xyzzyaaad2),100)
enddo
xyzzyaaay2(:)=atno(:)/=int(valence_charge(:))
xyzzyaaae2=0
have_ae=.false.
a:do xyzzyaaad2=1,xyzzyaaam2
do xyzzyaaac2=1,xyzzyaaad2-1
if(xyzzyaaar2(xyzzyaaac2)==atno(xyzzyaaad2))then
xyzzyaaar2(xyzzyaaad2)=xyzzyaaar2(xyzzyaaac2)
iontype_prim(xyzzyaaad2)=iontype_prim(xyzzyaaac2)
cycle a
endif
enddo
xyzzyaaae2=xyzzyaaae2+1
xyzzyaaar2(xyzzyaaad2)=atno(xyzzyaaad2)
iontype_prim(xyzzyaaad2)=xyzzyaaae2
enddo a
nitype=xyzzyaaae2
do xyzzyaaaf2=1,nitype
do xyzzyaaad2=1,xyzzyaaam2
if(iontype_prim(xyzzyaaad2)==xyzzyaaaf2)exit
enddo
xyzzyaaaa2=mod(atno(xyzzyaaad2),1000)
xyzzyaaab2=atno(xyzzyaaad2)/1000+1
ppchar=''
if(xyzzyaaab2>1)ppchar=trim(i2s(xyzzyaaab2))
psp_filename=trim(periodic_table_nocap(xyzzyaaaa2))//trim(adjustl(ppch&
&ar))//'_pp.data'
inquire(file=psp_filename,exist=xyzzyaaas2,err=7)
if(xyzzyaaax2(xyzzyaaad2).and..not.xyzzyaaas2)call errstop('READGW','T&
&he atomic numbers in this gwfn.data file follow the +200 convention t&
&o flag pseudoatoms. However, the x_pp.data pseudopotential file for a&
& flagged pseudoatom is missing.')
if(xyzzyaaay2(xyzzyaaad2).and..not.xyzzyaaas2)call errstop('READGW','T&
&he supplied valence charges in the gwfn.data file suggest the presenc&
&e of a pseudopotential, but the corresponding x_pp.data file does not&
& exist.')
if(.not.xyzzyaaas2)have_ae=.true.
enddo
deallocate(xyzzyaaax2,xyzzyaaay2,iontype_prim,xyzzyaaar2)
call skip(io,3)
if(any(wf_np>0).or.any(wf_nm>0))call errstop('READGW','Additions or su&
&btractions should not be used for Gaussian orbitals: number of electr&
&ons is simply determined by NEU and NED, and occupancy should only be&
& changed using "PR" in the MDET block.')
excite=(any(wf_nd(:,:)>0))
if(.not.modified_mdet)call read_mdet_wfn(io)
do
read(io,*,iostat=xyzzyaaak2)eigcheck
if(xyzzyaaak2/=0)call errstop('READGW','Problem finding orbital coeffi&
&cients.')
if(eigcheck=='EIGENVECTOR'.or.eigcheck=='ORBITAL')exit
enddo
call skip(io,1)
xyzzyaaav2=(any(wf_nd(:,:)>0))
xyzzyaaau2=.false.
if(.not.excite.and.xyzzyaaav2)then
excite=.true.
xyzzyaaau2=.true.
endif
endif
call mpi_bcast(num_k,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_k in readgw')
call mpi_bcast(num_real_k,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_real_k in readgw')
call mpi_bcast(num_complex_k,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_complex_k in readgw')
call mpi_bcast(num_ao,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_ao in readgw')
call mpi_bcast(xyzzyaaai1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting no_spin_types_in in readgw')
call mpi_bcast(detstart,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting detstart in readgw')
call mpi_bcast(detstop,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting detstop in readgw')
xyzzyaaah2=num_ao*num_ao
xyzzyaaai2=2*xyzzyaaah2
xyzzyaaag2=(num_real_k*xyzzyaaah2+num_complex_k*xyzzyaaai2)*xyzzyaaai1
call shallocate(xyzzyaaay1,xyzzyaaag2,stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','CK')
call shallocate_barrier
if(am_master)then
read(io,xyzzyaabb2,err=2,end=2)(xyzzyaaay1(xyzzyaaaa2),xyzzyaaaa2=1,xy&
&zzyaaag2)
call skip(io,1)
if(periodicity/=0)then
read(io,'(a)',err=5,end=5)eigcheck
if(scan(eigcheck,"EIGENVALUES",back=.true.)==11)then
allocate(xyzzyaaba1(num_ao,num_k,xyzzyaaai1),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','')
call skip(io,1)
do xyzzyaaac2=1,num_k
do xyzzyaaaj2=1,xyzzyaaai1
call skip(io,1)
read(io,xyzzyaabb2,err=6,end=6)(xyzzyaaba1(xyzzyaaaa2,xyzzyaaac2,xyzzy&
&aaaj2),xyzzyaaaa2=1,num_ao)
enddo
enddo
else
call errstop('READGW','Pre-2001 gwfn.data file format not supported.')
endif
endif
close(io)
open_unit(io)=.false.
endif
call mpi_bcast(xyzzyaaau2,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting rebroadcast_data in readgw')
if(xyzzyaaau2)then
call mpi_bcast(ndet,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ndet in readgw')
if(am_slave)then
deallocate(detcoef,wf_nd)
allocate(detcoef(ndet),wf_nd(ndet,2))
detcoef=1.d0
wf_nd=0
endif
call mpi_bcast(detcoef,ndet,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting detcoef in readgw')
call mpi_bcast(mods,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting mods in readgw')
if(am_slave)then
allocate(wf_d(4,mdet_max_mods,ndet,2))
wf_d=0
endif
call mpi_bcast(wf_nd,ndet*2,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting wf_nd in readgw')
call mpi_bcast(wf_d,8*mdet_max_mods*ndet,mpi_integer,0,mpi_comm_world,&
&ierror)
call checkmpi(ierror,'Broadcasting wf_d in readgw')
endif
call mpi_bcast(periodicity,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting periodicity in readgw')
call mpi_bcast(spin_dep_in,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting spin_dep_in in readgw')
call mpi_bcast(num_shells,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_shells in readgw')
call mpi_bcast(xyzzyaaac1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_prims in readgw')
call mpi_bcast(xyzzyaaam2,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_atoms in readgw')
call mpi_bcast(num_centres,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_centres in readgw')
call mpi_bcast(xyzzyaaag2,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ck_size in readgw')
call mpi_bcast(xyzzyaach1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting spin_restricted in readgw')
call mpi_bcast(complex_states,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting complex_states in readgw')
call mpi_bcast(xyzzyaaan2,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting highest_ang_mom in readgw')
call mpi_bcast(xyzzyaaab1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting gwfn_netot_per_cell in readgw')
call mpi_bcast(eionion,1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting eionion in readgw')
if(am_slave)then
num_real_k_plus_1=num_real_k+1
num_real_k_gt_1=num_real_k>1
num_real_k_odd=mod(num_real_k,2)>0
allocate(valence_charge(xyzzyaaam2),kvec(3,num_k),shell_am(num_shells)&
&,numao_in_shell(num_shells),first_shell(num_centres+1),first_ao(num_s&
&hells+1),primitive(num_shells+1),exponent(xyzzyaaac1),min_exponent(nu&
&m_shells),c_prim(xyzzyaaac1),xyzzyaaaz1(3,num_shells),xyzzyaaad1(num_&
&shells),atno(xyzzyaaam2),basis(3,xyzzyaaam2),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','')
if(periodicity>0)then
allocate(xyzzyaaba1(num_ao,num_k,xyzzyaaai1),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','')
endif
endif
call mpi_bcast(numao_in_shell,num_shells,mpi_integer,0,mpi_comm_world,&
&ierror)
call checkmpi(ierror,'Broadcasting numao_in_shell in readgw')
call mpi_bcast(primitive,num_shells+1,mpi_integer,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting primitive in readgw')
call mpi_bcast(shell_am,num_shells,mpi_integer,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'Broadcasting shell_am in readgw')
call mpi_bcast(first_ao,num_shells+1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting first_ao in readgw')
call mpi_bcast(first_shell,num_centres+1,mpi_integer,0,mpi_comm_world,&
&ierror)
call checkmpi(ierror,'Broadcasting first_shell in readgw')
call mpi_bcast(xyzzyaaad1,num_shells,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting numprims_in_shell in readgw')
call mpi_bcast(pamat,9,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pamat in readgw')
call mpi_bcast(pa1,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa1 in readgw')
call mpi_bcast(pa2,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa2 in readgw')
call mpi_bcast(pa3,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa3 in readgw')
call mpi_bcast(exponent,xyzzyaaac1,mpi_double_precision,0,mpi_comm_wor&
&ld,ierror)
call checkmpi(ierror,'Broadcasting exponent in readgw')
call mpi_bcast(c_prim,xyzzyaaac1,mpi_double_precision,0,mpi_comm_world&
&,ierror)
call checkmpi(ierror,'Broadcasting c_prim in readgw')
if(any(shell_am(:)==2))then
if(am_slave)then
allocate(c_prim2(xyzzyaaac1),stat=xyzzyaaal2)
call check_alloc(xyzzyaaal2,'READGW','c_prim2')
endif
call mpi_bcast(c_prim2,xyzzyaaac1,mpi_double_precision,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'Broadcasting c_prim2 in readgw')
endif
call mpi_bcast(basis,xyzzyaaam2*3,mpi_double_precision,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'Broadcasting basis in readgw')
call mpi_bcast(atno,xyzzyaaam2,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting atno in readgw')
call mpi_bcast(have_ae,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting have_ae in readgw')
call mpi_bcast(valence_charge,xyzzyaaam2,mpi_double_precision,0,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'Broadcasting valence_charge in readgw')
call mpi_bcast(xyzzyaaaz1,num_shells*3,mpi_double_precision,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting shell_pos in readgw')
call mpi_bcast(min_exponent,num_shells,mpi_double_precision,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting min_exponent in readgw')
call mpi_bcast(kvec,3*num_k,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting kvec in readgw')
if(periodicity>0)then
call mpi_bcast(xyzzyaaba1,num_ao*num_k*xyzzyaaai1,mpi_double_precision&
&,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting eigenvalues in readgw')
endif
if(nnodes>nnpsmp.and.am_smpmaster)then
call mpi_bcast_safe(xyzzyaaay1,xyzzyaaag2,mpi_double_precision,0,xyzzy&
&aaao2,ierror)
call checkmpi(ierror,'Broadcasting ck in readgw')
endif
nbasis=xyzzyaaam2
eionion=eionion*dble(nbasis)
if(all(atno(:)==0))model_system=.true.
ndim=periodicity
nnn(ndim+1:3)=0
call shallocate_barrier
return
1 call errstop('READGW','Cannot open gwfn.data.')
2 call errstop('READGW','Problem reading gwfn.data.')
5 call errstop('READGW','Error reading EIGENVALUE line in gwfn.data.')
6 call errstop('READGW','Error reading eigenvalues in gwfn.data.')
7 call errstop('READGW','Error inquiring about '//trim(psp_filename)//&
&'.')
end subroutine readgw
subroutine gwfdet_setup
use slaarnaan, only : check_kpoints
use slaarnabi,          only : use_gpcc
use slaarnabt,     only : get_numerical_orbmask
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3,xyzzyaaal3&
&,xyzzyaaam3,xyzzyaaan3,xyzzyaaao3,xyzzyaaap3,xyzzyaaaq3,xyzzyaaar3,xy&
&zzyaaas3,xyzzyaaat3,xyzzyaaau3,xyzzyaaav3,xyzzyaaaw3,xyzzyaaax3,xyzzy&
&aaay3,xyzzyaaaz3,xyzzyaaba3,xyzzyaabb3,xyzzyaabc3,xyzzyaabd3,xyzzyaab&
&e3,xyzzyaabf3,xyzzyaabg3,xyzzyaabh3,xyzzyaabi3,xyzzyaabj3,xyzzyaabk3,&
&xyzzyaabl3,xyzzyaabm3,xyzzyaabn3,xyzzyaabo3,xyzzyaabp3
integer,allocatable :: xyzzyaabq3(:),xyzzyaabr3(:),xyzzyaabs3(:),xyzzy&
&aabt3(:,:),xyzzyaabu3(:),xyzzyaabv3(:,:)
real(dp) xyzzyaabw3
real(dp),allocatable :: xyzzyaabx3(:)
logical xyzzyaaby3,xyzzyaabz3,xyzzyaaca3
character(80) tmpr
if(.not.isperiodic)then
xyzzyaaae1=1
else
if(.not.cusp_correction)then
xyzzyaaae1=2
else
xyzzyaaae1=3
endif
endif
rnorm=0.d0
do xyzzyaaaa3=2,nemax
rnorm=rnorm+log(dble(xyzzyaaaa3))
enddo
rnorm=exp(-rnorm/(2.d0*dble(nemax)))
rnorm=rnorm*orb_norm
one_over_rnorm=1.d0/rnorm
two_rnorm=2.d0*rnorm
xyzzyaabz3=.not.xyzzyaach1.or.(nuc_nele(1)/=nuc_nele(2))
if(nspin==3)then
if(xyzzyaabz3)then
spin_dep_gs=2
else
spin_dep_gs=1
endif
else
if(xyzzyaabz3)then
spin_dep_gs=1
else
spin_dep_gs=0
endif
endif
xyzzyaaah1=no_ssingles(spin_dep_gs)
if(spin_dep_gs<spin_dep_in)call errstop_master('GWFDET_SETUP','spin_de&
&p_gs<spin_dep_in. This is a bug.')
if((nspin/=2.and.nspin/=3).or.any(pcharge(1:2)/=-1.d0).or.any(pcharge(&
&1:2)/=-1.d0))call errstop_master('GWFDET_SETUP','First two particle t&
&ypes should be electrons. Can have one other component at present; ot&
&herwise will need to generalize the way that spin-dependence is speci&
&fied in gwfn.data.')
if(xyzzyaabz3.or.any(wf_nd(:,1)/=wf_nd(:,2)))then
spin_dep_full=1
else
spin_dep_full=0
idet_loop : do xyzzyaaav3=1,ndet
do xyzzyaaaw3=1,wf_nd(xyzzyaaav3,1)
if(.not.any(wf_d(1,1:wf_nd(xyzzyaaav3,2),xyzzyaaav3,2)==wf_d(1,xyzzyaa&
&aw3,xyzzyaaav3,1).and.wf_d(2,1:wf_nd(xyzzyaaav3,2),xyzzyaaav3,2)==wf_&
&d(2,xyzzyaaaw3,xyzzyaaav3,1)).or..not.any(wf_d(3,1:wf_nd(xyzzyaaav3,2&
&),xyzzyaaav3,2)==wf_d(3,xyzzyaaaw3,xyzzyaaav3,1).and.wf_d(4,1:wf_nd(x&
&yzzyaaav3,2),xyzzyaaav3,2)==wf_d(4,xyzzyaaaw3,xyzzyaaav3,1)))then
spin_dep_full=1
exit idet_loop
endif
enddo
enddo idet_loop
endif
if(nspin==3)spin_dep_full=spin_dep_full+1
xyzzyaaaj1=no_ssingles(spin_dep_full)
if(spin_dep_full<spin_dep_gs)call errstop_master('GWFDET_SETUP','spin_&
&dep_full<spin_dep_gs. This is a bug.')
xyzzyaaar3=0
xyzzyaaas3=0
xyzzyaaat3=0
do xyzzyaaaf3=1,nspin
xyzzyaaap3=which_ssingle(xyzzyaaaf3,spin_dep_in)
xyzzyaaag3=which_ssingle(xyzzyaaaf3,spin_dep_gs)
xyzzyaaaq3=which_ssingle(xyzzyaaaf3,spin_dep_full)
if(xyzzyaaap3<xyzzyaaar3)call errstop_master('GWFDET_SETUP','Spin '//t&
&rim(i2s(xyzzyaaaf3))//' has input spin-type '//trim(i2s(xyzzyaaap3))/&
&/' but spin '//trim(i2s(xyzzyaaaf3-1))//' has spin-type '//trim(i2s(x&
&yzzyaaar3))//'.')
if(xyzzyaaag3<xyzzyaaas3)call errstop_master('GWFDET_SETUP','Spin '//t&
&rim(i2s(xyzzyaaaf3))//' has ground-state spin-type '//trim(i2s(xyzzya&
&aag3))//' but spin '//trim(i2s(xyzzyaaaf3-1))//' has spin-type '//tri&
&m(i2s(xyzzyaaas3))//'.')
if(xyzzyaaap3<xyzzyaaat3)call errstop_master('GWFDET_SETUP','Spin '//t&
&rim(i2s(xyzzyaaaf3))//' has full spin-type '//trim(i2s(xyzzyaaaq3))//&
&' but spin '//trim(i2s(xyzzyaaaf3-1))//' has spin-type '//trim(i2s(xy&
&zzyaaat3))//'.')
if(xyzzyaaap3/=xyzzyaaar3.and.xyzzyaaag3==xyzzyaaas3)call errstop_mast&
&er('GWFDET_SETUP','Spins '//trim(i2s(xyzzyaaaf3-1))//' and '//trim(i2&
&s(xyzzyaaaf3))//' have the same ground-state spin type but different &
&input spin types.')
if(xyzzyaaag3/=xyzzyaaas3.and.xyzzyaaaq3==xyzzyaaat3)call errstop_mast&
&er('GWFDET_SETUP','Spins '//trim(i2s(xyzzyaaaf3-1))//' and '//trim(i2&
&s(xyzzyaaaf3))//' have the same full spin type but different ground-s&
&tate spin types.')
xyzzyaaar3=xyzzyaaap3
xyzzyaaas3=xyzzyaaag3
xyzzyaaat3=xyzzyaaaq3
enddo
xyzzyaaaa3=xyzzyaaab1*npcells
if(xyzzyaaaa3/=sum(nuc_nele(:)))then
if(am_master)then
call errwarn_silent('GWFDET_SETUP','Number of up/down electrons in inp&
&ut file differs from that provided in gwfn.data. You may want to chec&
&k the values of NEU and NED in the input file.')
call wout()
call wout('npcells : '//trim(i2s(npcells)))
call wout('gwfn_netot_per_cell*npcells-excitations : '//trim(i2s(xyzzy&
&aaaa3)))
if(electron_system.and.nspin==2)then
call wout('neu : '//trim(i2s(nuc_nele(1)))//' ned : '//trim(i2s(nuc_ne&
&le(2))))
else
write(tmpr,*)nele(1:nspin)
call wout('Particle numbers: '//trim(adjustl(tmpr)))
endif
call wout()
endif
endif
allocate(xyzzyaaax1(ndet,num_ao,num_k,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','band_occupancy')
xyzzyaaax1=0
xyzzyaabk3=0
xyzzyaabj3=0
xyzzyaabm3=0
xyzzyaabl3=0
xyzzyaaaq1=0
nemaxc=nemax
call check_kpoints(num_k,kvec)
if(isperiodic)then
allocate(xyzzyaabx3(num_ao*num_k),xyzzyaabr3(num_ao*num_k),xyzzyaabq3(&
&num_ao*num_k),xyzzyaabs3(num_ao*(num_real_k+num_complex_k*2)),ridx2zi&
&dx(num_ao*(num_real_k+num_complex_k*2)),nband(num_k,xyzzyaaah1),xyzzy&
&aabv3(num_k,xyzzyaaah1),xyzzyaaau1(num_ao*num_k),xyzzyaaav1(num_ao*nu&
&m_k),xyzzyaaaw1(num_ao*num_k),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','eigtemp,ktemp,indx,etc..')
xyzzyaabx3=0.d0
xyzzyaabr3=0
xyzzyaabq3=0
xyzzyaabs3=0
ridx2zidx=0
nband=0
xyzzyaabv3=0
xyzzyaaau1=0
xyzzyaaav1=0
xyzzyaaaw1=0
xyzzyaaas3=0
do xyzzyaaaf3=1,nspin
if(nuc_nele(xyzzyaaaf3)==0)cycle
xyzzyaaag3=which_ssingle(xyzzyaaaf3,spin_dep_gs)
if(xyzzyaaag3==xyzzyaaas3)cycle
xyzzyaaas3=xyzzyaaag3
xyzzyaaap3=which_ssingle(xyzzyaaaf3,spin_dep_in)
xyzzyaaai3=0
do xyzzyaaac3=1,num_k
do xyzzyaaah3=1,num_ao
xyzzyaaai3=xyzzyaaai3+1
xyzzyaabx3(xyzzyaaai3)=xyzzyaaba1(xyzzyaaah3,xyzzyaaac3,xyzzyaaap3)
xyzzyaabr3(xyzzyaaai3)=xyzzyaaac3
enddo
enddo
call quicksort(xyzzyaaai3,xyzzyaabx3(1),xyzzyaabq3(1))
xyzzyaaal3=0
do xyzzyaaaa3=1,xyzzyaaai3
xyzzyaaac3=xyzzyaabr3(xyzzyaabq3(xyzzyaaaa3))
nband(xyzzyaaac3,xyzzyaaag3)=nband(xyzzyaaac3,xyzzyaaag3)+1
xyzzyaaal3=xyzzyaaal3+1
xyzzyaabs3(xyzzyaaal3)=xyzzyaaaa3
if(xyzzyaaac3>num_real_k)then
xyzzyaaal3=xyzzyaaal3+1
xyzzyaabs3(xyzzyaaal3)=xyzzyaaaa3
endif
if(xyzzyaaal3>=nuc_nele(xyzzyaaaf3))exit
enddo
if(xyzzyaaal3>nuc_nele(xyzzyaaaf3))then
nemaxc=nemax+1
call errwarn_silent('GWFDET_SETUP','Highest occupied level is a singly&
& occupied +k/-k orbital pair.')
endif
do xyzzyaaau3=xyzzyaaaf3,nspin
if(xyzzyaaag3/=which_ssingle(xyzzyaaau3,spin_dep_gs))exit
do xyzzyaaaa3=1,nuc_nele(xyzzyaaau3)
xyzzyaaac3=xyzzyaabr3(xyzzyaabq3(xyzzyaabs3(xyzzyaaaa3)))
do xyzzyaaah3=1,nband(xyzzyaaac3,xyzzyaaag3)
if((xyzzyaaac3<=num_real_k.and.xyzzyaaax1(1,xyzzyaaah3,xyzzyaaac3,xyzz&
&yaaau3)<1) .or.(xyzzyaaac3>num_real_k.and.xyzzyaaax1(1,xyzzyaaah3,xyz&
&zyaaac3,xyzzyaaau3)<2))then
xyzzyaaax1(1,xyzzyaaah3,xyzzyaaac3,xyzzyaaau3)=xyzzyaaax1(1,xyzzyaaah3&
&,xyzzyaaac3,xyzzyaaau3)+1
exit
endif
enddo
enddo
enddo
if(am_master.and.nuc_nele(xyzzyaaaf3)<xyzzyaaai3.and.nuc_nele(xyzzyaaa&
&f3)>0)then
xyzzyaaae3=xyzzyaabs3(nuc_nele(xyzzyaaaf3))
if(abs(xyzzyaabx3(xyzzyaabq3(xyzzyaaae3))-xyzzyaabx3(xyzzyaabq3(xyzzya&
&aae3+1)))<1.d-6)then
call wout()
if(xyzzyaaah1>1)then
if(electron_system)then
if(xyzzyaaaf3==1)then
call wout('UP SPIN')
else
call wout('DOWN SPIN')
endif
else
call wout('SPIN '//trim(i2s(xyzzyaaaf3)))
endif
endif
xyzzyaabw3=abs(xyzzyaabx3(xyzzyaabq3(xyzzyaaae3)))
if(abs(xyzzyaabx3(xyzzyaabq3(xyzzyaaae3+1)))==xyzzyaabw3)then
call wout('Partially occupied strictly degenerate states at Fermi leve&
&l.')
else
call wout('Partially occupied nearly degenerate states at Fermi level.&
&')
endif
xyzzyaaam3=1
do xyzzyaaaa3=xyzzyaaae3+1,(num_ao*num_k)-1
if(abs(xyzzyaabx3(xyzzyaabq3(xyzzyaaaa3)))==xyzzyaabw3)then
xyzzyaaam3=xyzzyaaam3+1
else
exit
endif
enddo
do xyzzyaaaa3=xyzzyaaae3-1,2,-1
if(abs(xyzzyaabx3(xyzzyaabq3(xyzzyaaaa3)))==xyzzyaabw3)then
xyzzyaaam3=xyzzyaaam3+1
else
exit
endif
enddo
xyzzyaaan3=1
do xyzzyaaaa3=xyzzyaaae3+1,(num_ao*num_k)-1
if(abs(xyzzyaabx3(xyzzyaabq3(xyzzyaaaa3)))-xyzzyaabw3<1.d-6)then
xyzzyaaan3=xyzzyaaan3+1
else
exit
endif
enddo
do xyzzyaaaa3=xyzzyaaae3-1,2,-1
if(abs(abs(xyzzyaabx3(xyzzyaabq3(xyzzyaaaa3)))-xyzzyaabw3)<1.d-6)then
xyzzyaaan3=xyzzyaaan3+1
else
exit
endif
enddo
if(xyzzyaaan3>1)then
call wout()
call wout('Degeneracy of highest occupied level (exact)         :  '//&
&trim(i2s(xyzzyaaam3)))
call wout('                                     (within 10^-6)  :  '//&
&trim(i2s(xyzzyaaan3)))
endif
endif
endif
enddo
if(am_master)then
call wout('State occupation')
call wout('================')
xyzzyaaby3=.false.
do xyzzyaaag3=1,xyzzyaaah1
if(any(nband(:,xyzzyaaag3)/=nband(1,xyzzyaaag3)))then
xyzzyaaby3=.true.
exit
endif
enddo
if(xyzzyaaby3)then
call wout('METALLIC STATE DETECTED')
if(electron_system.and.nspin==2)then
if(.not.xyzzyaabz3)then
call wout('Number of doubly-occupied states filled at each k point:')
call wout('    k    nband')
do xyzzyaaac3=1,num_k
write(tmpr,'(i5,1x,i5)')xyzzyaaac3,nband(xyzzyaaac3,1)
call wout(tmpr)
enddo
else
call wout('Number of singly-occupied states filled at each k point:')
call wout('UP SPIN')
call wout('    k    nband')
do xyzzyaaac3=1,num_k
write(tmpr,'(i5,1x,i5)')xyzzyaaac3,nband(xyzzyaaac3,1)
call wout(tmpr)
enddo
call wout('DOWN SPIN')
call wout('    k    nband')
do xyzzyaaac3=1,num_k
write(tmpr,'(i5,1x,i5)')xyzzyaaac3,nband(xyzzyaaac3,2)
call wout(tmpr)
enddo
endif
else
call wout('Number of states filled at each k point:')
do xyzzyaaag3=1,xyzzyaaah1
call wout('GROUND-STATE SPIN TYPE '//trim(i2s(xyzzyaaag3))//'     (spi&
&n dependence '//trim(i2s(spin_dep_gs))//')')
call wout('    k    nband')
do xyzzyaaac3=1,num_k
write(tmpr,'(i5,1x,i5)')xyzzyaaac3,nband(xyzzyaaac3,xyzzyaaag3)
call wout(tmpr)
enddo
enddo
endif
else
call wout('INSULATING STATE DETECTED')
if(electron_system.and.nspin==2)then
if(.not.xyzzyaabz3)then
call wout('Number of doubly-occupied states filled at each k point  : &
& '//trim(i2s(nband(1,1))))
else
call wout('Number of singly-occupied states filled at each k point  : &
& '//trim(i2s(sum(nband(1,1:2)))))
endif
else
do xyzzyaaag3=1,xyzzyaaah1
call wout('Number of states for ground-state spin-type ' //trim(i2s(xy&
&zzyaaag3))//' at each k: ' //trim(i2s(nband(1,xyzzyaaag3))))
enddo
endif
endif
call wout()
endif
xyzzyaaas3=0
xyzzyaaar3=0
do xyzzyaaaf3=1,nspin
if(nuc_nele(xyzzyaaaf3)==0)cycle
xyzzyaaag3=which_ssingle(xyzzyaaaf3,spin_dep_gs)
if(xyzzyaaag3==xyzzyaaas3)cycle
xyzzyaaas3=xyzzyaaag3
xyzzyaaap3=which_ssingle(xyzzyaaaf3,spin_dep_in)
if(xyzzyaaap3==xyzzyaaar3)then
xyzzyaabj3=0
xyzzyaabk3=0
endif
xyzzyaaar3=xyzzyaaap3
xyzzyaaal3=0
k_loop: do xyzzyaaac3=1,num_k
do xyzzyaaah3=1,nband(xyzzyaaac3,xyzzyaaag3)
xyzzyaaal3=xyzzyaaal3+1
xyzzyaabj3=xyzzyaabj3+1
xyzzyaabk3=xyzzyaabk3+1
ridx2zidx(xyzzyaabj3)=xyzzyaabk3
if(xyzzyaaac3>num_real_k)then
xyzzyaaal3=xyzzyaaal3+1
if(xyzzyaaal3<=nuc_nele(xyzzyaaaf3))then
xyzzyaabj3=xyzzyaabj3+1
ridx2zidx(xyzzyaabj3)=xyzzyaabk3
endif
endif
xyzzyaaau1(xyzzyaabk3)=xyzzyaaac3
xyzzyaaav1(xyzzyaabk3)=xyzzyaaah3
xyzzyaaaw1(xyzzyaabk3)=xyzzyaaap3
if(xyzzyaaal3>=nuc_nele(xyzzyaaaf3))exit k_loop
enddo
enddo k_loop
xyzzyaabm3=max(xyzzyaabm3,xyzzyaabk3)
xyzzyaabl3=max(xyzzyaabl3,xyzzyaabj3)
enddo
xyzzyaaaq1=xyzzyaabm3
deallocate(xyzzyaabq3,xyzzyaabx3,xyzzyaabr3,xyzzyaabs3)
else
allocate(nband(1,xyzzyaaah1),xyzzyaabv3(1,xyzzyaaah1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','')
xyzzyaaar3=0
do xyzzyaaaf3=1,nspin
xyzzyaaag3=which_ssingle(xyzzyaaaf3,spin_dep_gs)
xyzzyaaap3=which_ssingle(xyzzyaaaf3,spin_dep_in)
nband(1,xyzzyaaag3)=nuc_nele(xyzzyaaaf3)
xyzzyaaax1(1,1:nband(1,xyzzyaaag3),1,xyzzyaaaf3)=1
if(xyzzyaaap3==xyzzyaaar3)then
xyzzyaabl3=max(xyzzyaabl3,nuc_nele(xyzzyaaaf3))
else
xyzzyaabl3=xyzzyaabl3+nuc_nele(xyzzyaaaf3)
endif
xyzzyaaar3=xyzzyaaap3
enddo
xyzzyaabv3=0
endif
gauss_gs_norb=xyzzyaabl3
do xyzzyaaav3=2,ndet
xyzzyaaax1(xyzzyaaav3,:,:,:)=xyzzyaaax1(1,:,:,:)
enddo
if(any(wf_nd(:,:)>0))then
excite=.true.
else
excite=.false.
endif
allocate(xyzzyaaak1(nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','num_vo')
xyzzyaaak1=0
gauss_ex_norb=0
if(excite)then
if(use_backflow.and.isperiodic)call errwarn('GWFDET_SETUP','No cusp co&
&rrections for virtual orbitals in periodic systems with backflow yet.&
&')
if(any(wf_d(2,:,:,:)>num_k).or.any(wf_d(2,:,:,:)<0).or.any(wf_d(4,:,:,&
&:)>num_k).or.any(wf_d(4,:,:,:)<0))call errstop_master('GWFDET_SETUP',&
&'Excited state k out of range.')
if(any(wf_d(3,:,:,:)>num_ao))call errstop_master('GWFDET_SETUP','Insuf&
&ficient basis functions to define requested virtual orbital(s).')
allocate(xyzzyaaal1(mdet_max_mods*ndet,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','')
xyzzyaaaa3=0
xyzzyaaab3=0
do xyzzyaaag3=1,xyzzyaaah1
xyzzyaaax3=sum(num_ao-nband(1:num_k,xyzzyaaag3))
if(xyzzyaaax3>xyzzyaaaa3)xyzzyaaaa3=xyzzyaaax3
xyzzyaaab3=xyzzyaaab3+xyzzyaaax3
enddo
allocate(xyzzyaabt3(xyzzyaaaa3,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','')
xyzzyaabt3(:,:)=0
allocate(full2vrt(xyzzyaaab3,nspin),vrt2full(xyzzyaaaa3,nspin),stat=xy&
&zzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','')
full2vrt(:,:)=0
vrt2full(:,:)=0
if(cusp_correction)then
allocate(xyzzyaaao1(xyzzyaaaa3,nspin),xyzzyaaan1(xyzzyaaaa3,nspin),sta&
&t=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','')
endif
allocate(xyzzyaaam1(mdet_max_mods*ndet,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','')
allocate(virtual_k(mdet_max_mods*ndet,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','')
virtual_k(:,:)=0
do xyzzyaaaf3=1,nspin
xyzzyaaag3=which_ssingle(xyzzyaaaf3,spin_dep_gs)
xyzzyaaap3=which_ssingle(xyzzyaaaf3,spin_dep_in)
xyzzyaaad3=0
xyzzyaaaj3=0
do xyzzyaaav3=1,ndet
do xyzzyaaay3=1,wf_nd(xyzzyaaav3,xyzzyaaaf3)
xyzzyaaaz3=wf_d(1,xyzzyaaay3,xyzzyaaav3,xyzzyaaaf3)
xyzzyaaba3=wf_d(2,xyzzyaaay3,xyzzyaaav3,xyzzyaaaf3)
xyzzyaabb3=wf_d(3,xyzzyaaay3,xyzzyaaav3,xyzzyaaaf3)
xyzzyaabc3=wf_d(4,xyzzyaaay3,xyzzyaaav3,xyzzyaaaf3)
if(xyzzyaaaz3>nband(xyzzyaaba3,xyzzyaaag3).and.am_master)then
call errstop_master('GWFDET_SETUP','Promotion problem - no electron in&
& specified orbital.')
endif
if(xyzzyaabb3<=nband(xyzzyaabc3,xyzzyaaag3).and.am_master)then
call errstop_master('GWFDET_SETUP','Specified orbital for promoted ele&
&ctron already occupied.')
endif
xyzzyaaax1(xyzzyaaav3,xyzzyaabb3,xyzzyaabc3,xyzzyaaaf3)=xyzzyaaax1(xyz&
&zyaaav3,xyzzyaabb3,xyzzyaabc3,xyzzyaaaf3)+1
xyzzyaaax1(xyzzyaaav3,xyzzyaaaz3,xyzzyaaba3,xyzzyaaaf3)=xyzzyaaax1(xyz&
&zyaaav3,xyzzyaaaz3,xyzzyaaba3,xyzzyaaaf3)-1
if(isperiodic)then
xyzzyaaca3=.false.
do xyzzyaabk3=1,xyzzyaaaq1
if(xyzzyaaav1(xyzzyaabk3)==xyzzyaabb3.and.xyzzyaaau1(xyzzyaabk3)==xyzz&
&yaabc3.and.xyzzyaaaw1(xyzzyaabk3)==xyzzyaaap3)then
xyzzyaaca3=.true.
exit
endif
enddo
if(.not.xyzzyaaca3)then
xyzzyaaaq1=xyzzyaaaq1+1
xyzzyaabk3=xyzzyaaaq1
xyzzyaaau1(xyzzyaabk3)=xyzzyaabc3
xyzzyaaav1(xyzzyaabk3)=xyzzyaabb3
xyzzyaaaw1(xyzzyaabk3)=xyzzyaaap3
endif
endif
xyzzyaabd3=0
a:   do xyzzyaaac3=1,num_k
do xyzzyaaah3=nband(xyzzyaaac3,xyzzyaaag3)+1,num_ao
xyzzyaabd3=xyzzyaabd3+1
if(xyzzyaaac3==xyzzyaabc3.and.xyzzyaaah3==xyzzyaabb3)then
xyzzyaaaj3=xyzzyaaaj3+1
if(xyzzyaabt3(xyzzyaabd3,xyzzyaaaf3)==0)then
xyzzyaaak1(xyzzyaaaf3)=xyzzyaaak1(xyzzyaaaf3)+1
gauss_ex_norb=gauss_ex_norb+1
xyzzyaabt3(xyzzyaabd3,xyzzyaaaf3)=gauss_ex_norb
xyzzyaaam1(xyzzyaaaj3,xyzzyaaaf3)=gauss_ex_norb
full2vrt(gauss_ex_norb,xyzzyaaaf3)=xyzzyaaak1(xyzzyaaaf3)
vrt2full(xyzzyaaak1(xyzzyaaaf3),xyzzyaaaf3)=gauss_ex_norb
if(isperiodic)then
ridx2zidx(gauss_gs_norb+gauss_ex_norb)=xyzzyaabk3
virtual_k(gauss_ex_norb,xyzzyaaaf3)=xyzzyaabc3
endif
else
xyzzyaaam1(xyzzyaaaj3,xyzzyaaaf3)=xyzzyaabt3(xyzzyaabd3,xyzzyaaaf3)
endif
exit a
endif
enddo
enddo a
xyzzyaaab3=0
b:   do xyzzyaaac3=1,num_k
do xyzzyaaah3=1,nband(xyzzyaaac3,xyzzyaaag3)
xyzzyaaab3=xyzzyaaab3+1
if(xyzzyaaac3==xyzzyaaba3.and.xyzzyaaah3==xyzzyaaaz3)then
xyzzyaaad3=xyzzyaaad3+1
xyzzyaaal1(xyzzyaaad3,xyzzyaaaf3)=xyzzyaaab3
if(cusp_correction)xyzzyaaan1(xyzzyaaak1(xyzzyaaaf3),xyzzyaaaf3)=xyzzy&
&aaab3
exit b
endif
if(xyzzyaaac3>num_real_k)xyzzyaaab3=xyzzyaaab3+1
enddo
enddo b
if(cusp_correction)xyzzyaaao1(xyzzyaaak1(xyzzyaaaf3),xyzzyaaaf3)=xyzzy&
&aaav3
enddo
enddo
enddo
endif
do xyzzyaaaf3=1,nspin
do xyzzyaaav3=1,ndet
if(sum(xyzzyaaax1(xyzzyaaav3,:,:,xyzzyaaaf3))/=nuc_nele(xyzzyaaaf3))ca&
&ll errstop_master('GWFDET_SETUP','Problem with excitations (1) : mism&
&atch between band_occupancy and nuc_nele.')
if(any(xyzzyaaax1(xyzzyaaav3,:,:,xyzzyaaaf3)<0).or.any(xyzzyaaax1(xyzz&
&yaaav3,:,:,xyzzyaaaf3)>2))call errstop_master('GWFDET_SETUP','Problem&
& with excitations (2) : problem with band_occupancy vector.')
enddo
enddo
allocate(xyzzyaaar1(maxval(nuc_nele),nspin,ndet),gauss_offset(nspin),g&
&auss_offsetc(nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','gauss_orbmap,gauss_offset,&
&gauss_offsetc')
xyzzyaaar1=0
xyzzyaabi3=0
xyzzyaaar3=0
do xyzzyaaaf3=1,nspin
xyzzyaaap3=which_ssingle(xyzzyaaaf3,spin_dep_in)
if(xyzzyaaap3/=xyzzyaaar3)then
do xyzzyaaaa3=1,nuc_nele(xyzzyaaaf3)
xyzzyaabi3=xyzzyaabi3+1
xyzzyaaar1(xyzzyaaaa3,xyzzyaaaf3,1)=xyzzyaabi3
enddo
else
xyzzyaaau3=sum(maxloc(nuc_nele(:xyzzyaaaf3-1),which_ssingle(:xyzzyaaaf&
&3-1,spin_dep_in)==xyzzyaaap3))
xyzzyaabh3=0
do xyzzyaaaa3=1,nuc_nele(xyzzyaaau3)
xyzzyaabh3=xyzzyaabh3+1
xyzzyaaar1(xyzzyaaaa3,xyzzyaaaf3,1)=xyzzyaaar1(xyzzyaabh3,xyzzyaaau3,1&
&)
enddo
do xyzzyaaaa3=nuc_nele(xyzzyaaau3)+1,nuc_nele(xyzzyaaaf3)
xyzzyaabi3=xyzzyaabi3+1
xyzzyaaar1(xyzzyaaaa3,xyzzyaaaf3,1)=xyzzyaabi3
enddo
endif
xyzzyaaar3=xyzzyaaap3
gauss_offset(xyzzyaaaf3)=xyzzyaaar1(1,xyzzyaaaf3,1)
if(isperiodic.and.(nuc_nele(xyzzyaaaf3)/=0))gauss_offsetc(xyzzyaaaf3)=&
&ridx2zidx(xyzzyaaar1(1,xyzzyaaaf3,1))
do xyzzyaaav3=2,ndet
xyzzyaaar1(:,xyzzyaaaf3,xyzzyaaav3)=xyzzyaaar1(:,xyzzyaaaf3,1)
enddo
enddo
if(excite)then
do xyzzyaaaf3=1,nspin
xyzzyaaaj3=0
do xyzzyaaav3=1,ndet
do xyzzyaaay3=1,wf_nd(xyzzyaaav3,xyzzyaaaf3)
xyzzyaaaj3=xyzzyaaaj3+1
xyzzyaaar1(xyzzyaaal1(xyzzyaaaj3,xyzzyaaaf3),xyzzyaaaf3,xyzzyaaav3)=ga&
&uss_gs_norb+xyzzyaaam1(xyzzyaaaj3,xyzzyaaaf3)
enddo
enddo
enddo
endif
xyzzyaaap1=gauss_gs_norb+gauss_ex_norb
if(isperiodic)then
xyzzyaaci1=2
xyzzyaacj1=1
allocate(xyzzyaack1(xyzzyaaci1,xyzzyaaap1),xyzzyaacl1(xyzzyaacj1,xyzzy&
&aaap1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','gauss_orbdesc_*')
xyzzyaack1=0
xyzzyaacl1=0.d0
xyzzyaabo3=0
do xyzzyaabh3=1,xyzzyaaap1
xyzzyaabk3=ridx2zidx(xyzzyaabh3)
xyzzyaaac3=xyzzyaaau1(xyzzyaabk3)
xyzzyaaah3=xyzzyaaav1(xyzzyaabk3)
xyzzyaaaf3=xyzzyaaaw1(xyzzyaabk3)
xyzzyaack1(1,xyzzyaabh3)=xyzzyaaac3
xyzzyaaca3=.false.
do xyzzyaabp3=1,xyzzyaabh3-1
if(xyzzyaaba1(xyzzyaaah3,xyzzyaaac3,xyzzyaaaf3)==xyzzyaacl1(1,xyzzyaab&
&p3))then
xyzzyaaca3=.true.
exit
endif
enddo
if(.not.xyzzyaaca3)then
xyzzyaabo3=xyzzyaabo3+1
xyzzyaabn3=xyzzyaabo3
xyzzyaack1(2,xyzzyaabh3)=xyzzyaabn3
xyzzyaacl1(1,xyzzyaabh3)=xyzzyaaba1(xyzzyaaah3,xyzzyaaac3,xyzzyaaaf3)
else
xyzzyaack1(2,xyzzyaabh3)=xyzzyaack1(2,xyzzyaabp3)
xyzzyaacl1(1,xyzzyaabh3)=xyzzyaacl1(1,xyzzyaabp3)
endif
enddo
deallocate(xyzzyaaba1)
else
xyzzyaaci1=0
xyzzyaacj1=0
endif
if(excite)deallocate(xyzzyaaam1)
if(isperiodic)then
allocate(xyzzyaabu3(xyzzyaaap1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','itrf <1>')
xyzzyaabu3(1:xyzzyaaap1)=ridx2zidx(1:xyzzyaaap1)
deallocate(ridx2zidx)
allocate(ridx2zidx(xyzzyaaap1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','Ridx2Zidx')
ridx2zidx(1:xyzzyaaap1)=xyzzyaabu3(1:xyzzyaaap1)
deallocate(xyzzyaabu3)
allocate(xyzzyaabu3(xyzzyaaaq1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','itrf <2>')
xyzzyaabu3(1:xyzzyaaaq1)=xyzzyaaav1(1:xyzzyaaaq1)
deallocate(xyzzyaaav1)
allocate(xyzzyaaav1(xyzzyaaaq1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','corb_band_idx')
xyzzyaaav1(1:xyzzyaaaq1)=xyzzyaabu3(1:xyzzyaaaq1)
deallocate(xyzzyaabu3)
allocate(xyzzyaabu3(xyzzyaaaq1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','itrf <3>')
xyzzyaabu3(1:xyzzyaaaq1)=xyzzyaaau1(1:xyzzyaaaq1)
deallocate(xyzzyaaau1)
allocate(xyzzyaaau1(xyzzyaaaq1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','corb_kvec_idx')
xyzzyaaau1(1:xyzzyaaaq1)=xyzzyaabu3(1:xyzzyaaaq1)
deallocate(xyzzyaabu3)
allocate(xyzzyaabu3(xyzzyaaaq1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','itrf <4>')
xyzzyaabu3(1:xyzzyaaaq1)=xyzzyaaaw1(1:xyzzyaaaq1)
deallocate(xyzzyaaaw1)
allocate(xyzzyaaaw1(xyzzyaaaq1),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','corb_spin_idx')
xyzzyaaaw1(1:xyzzyaaaq1)=xyzzyaabu3(1:xyzzyaaaq1)
deallocate(xyzzyaabu3)
endif
maxb=maxval(nband)
xyzzyaabv3(1:num_k,1:xyzzyaaah1)=num_ao-nband(1:num_k,1:xyzzyaaah1)
if(num_real_k>0)then
call shallocate(rck,maxb,num_ao,num_real_k,xyzzyaaai1,xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','rck')
size_rck=maxb*num_ao*num_real_k
endif
if(complex_states)then
call shallocate_gauss(cck,maxb,num_ao,num_real_k_plus_1,num_k,xyzzyaaa&
&i1,xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','cck')
size_cck=maxb*num_ao*num_complex_k
endif
allocate(orb1(nemaxc),orb2(nemaxc),orb3(nemaxc),orb4(nemaxc),orb5(nema&
&xc),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','orb[1-5]')
if(use_backflow)then
allocate(orb6(nemaxc),orb7(nemaxc),orb8(nemaxc),orb9(nemaxc),orb10(nem&
&axc),orb11(nemaxc),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','orb[6-11]')
endif
if(complex_states)then
allocate(corb1(maxb),corb2(maxb),corb3(maxb),corb4(maxb),corb5(maxb),s&
&tat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','corb[1-5]')
if(use_backflow)then
allocate(corb6(nemax),corb7(nemax),corb8(nemax),corb9(nemax),corb10(ne&
&max),corb11(nemax),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','corb[6-11]')
endif
endif
if(am_smpmaster)then
xyzzyaaas3=0
do xyzzyaaaf3=1,nspin
xyzzyaaag3=which_ssingle(xyzzyaaaf3,spin_dep_gs)
if(xyzzyaaag3==xyzzyaaas3)cycle
xyzzyaaas3=xyzzyaaag3
xyzzyaaap3=which_ssingle(xyzzyaaaf3,spin_dep_in)
xyzzyaaak3=num_ao*num_ao*(num_real_k+2*num_complex_k)*(xyzzyaaap3-1)
do xyzzyaaac3=1,num_real_k
do xyzzyaaah3=1,nband(xyzzyaaac3,xyzzyaaag3)
rck(xyzzyaaah3,1:num_ao,xyzzyaaac3,xyzzyaaap3)=xyzzyaaay1(xyzzyaaak3+1&
&:xyzzyaaak3+num_ao)
xyzzyaaak3=xyzzyaaak3+num_ao
enddo
xyzzyaaak3=xyzzyaaak3+num_ao*xyzzyaabv3(xyzzyaaac3,xyzzyaaag3)
enddo
do xyzzyaaac3=num_real_k_plus_1,num_k
do xyzzyaaah3=1,nband(xyzzyaaac3,xyzzyaaag3)
xyzzyaaab3=1
do xyzzyaaaa3=1,num_ao
cck(xyzzyaaah3,xyzzyaaaa3,xyzzyaaac3,xyzzyaaap3)=cmplx(xyzzyaaay1(xyzz&
&yaaak3+xyzzyaaab3),xyzzyaaay1(xyzzyaaak3+xyzzyaaab3+1),kind=dp)
xyzzyaaab3=xyzzyaaab3+2
enddo
xyzzyaaak3=xyzzyaaak3+2*num_ao
enddo
xyzzyaaak3=xyzzyaaak3+2*num_ao*xyzzyaabv3(xyzzyaaac3,xyzzyaaag3)
enddo
enddo
endif
if(allocated(xyzzyaabv3))deallocate(xyzzyaabv3)
if(num_real_k>0)then
allocate(rbf(num_ao,num_real_k),rblap(num_ao,num_real_k),rbgra1(num_ao&
&,num_real_k),rbgra2(num_ao,num_real_k),rbgra3(num_ao,num_real_k),phas&
&e(num_real_k),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','rbf etc.')
if(cusp_correction)then
allocate(rbf_s(num_ao,maxb*num_real_k),rblap_s(num_ao,maxb*num_real_k)&
&,rbgra1_s(num_ao,maxb*num_real_k),rbgra2_s(num_ao,maxb*num_real_k),rb&
&gra3_s(num_ao,maxb*num_real_k),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','rblap etc.')
if(use_backflow)then
allocate(rbsderiv1_s(num_ao,maxb*num_real_k),rbsderiv2_s(num_ao,maxb*n&
&um_real_k),rbsderiv3_s(num_ao,maxb*num_real_k),rbsderiv4_s(num_ao,max&
&b*num_real_k),rbsderiv5_s(num_ao,maxb*num_real_k),rbsderiv6_s(num_ao,&
&maxb*num_real_k),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','rbsderivs1_s etc.')
endif
endif
if(use_backflow)then
allocate(rbsderiv1(num_ao,num_real_k),rbsderiv2(num_ao,num_real_k),rbs&
&deriv3(num_ao,num_real_k),rbsderiv4(num_ao,num_real_k),rbsderiv5(num_&
&ao,num_real_k),rbsderiv6(num_ao,num_real_k),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','rbsderiv1 etc.')
endif
else
allocate(rbf(0,0),rblap(0,0),rbgra1(0,0),rbgra2(0,0),rbgra3(0,0))
if(cusp_correction)allocate(rbf_s(0,0),rblap_s(0,0),rbgra1_s(0,0),rbgr&
&a2_s(0,0),rbgra3_s(0,0))
if(use_backflow)allocate(rbsderiv1(0,0),rbsderiv2(0,0),rbsderiv3(0,0),&
&rbsderiv4(0,0),rbsderiv5(0,0),rbsderiv6(0,0))
endif
if(num_complex_k>0)then
allocate(cbf(num_ao,num_real_k_plus_1:num_k),cblap(num_ao,num_real_k_p&
&lus_1:num_k),cbgra1(num_ao,num_real_k_plus_1:num_k),cbgra2(num_ao,num&
&_real_k_plus_1:num_k),cbgra3(num_ao,num_real_k_plus_1:num_k),cphase(n&
&um_real_k_plus_1:num_k),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','cbf etc.')
if(cusp_correction)then
xyzzyaabf3=huge(1)
xyzzyaabg3=0
do xyzzyaaag3=1,xyzzyaaah1
xyzzyaaae3=sum(nband(1:num_real_k,xyzzyaaag3))
if(xyzzyaaae3<xyzzyaabf3)xyzzyaabf3=xyzzyaaae3
xyzzyaaae3=2*sum(nband(1:num_k,xyzzyaaag3))-sum(nband(1:num_real_k,xyz&
&zyaaag3))
if(xyzzyaaae3>xyzzyaabg3)xyzzyaabg3=xyzzyaaae3
enddo
allocate(cbf_s(num_ao,xyzzyaabf3:xyzzyaabg3),cblap_s(num_ao,xyzzyaabf3&
&:xyzzyaabg3),cbgra1_s(num_ao,xyzzyaabf3:xyzzyaabg3),cbgra2_s(num_ao,x&
&yzzyaabf3:xyzzyaabg3),cbgra3_s(num_ao,xyzzyaabf3:xyzzyaabg3),stat=xyz&
&zyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','cbf_s etc.')
if(use_backflow)then
allocate(cbsderiv1_s(num_ao,xyzzyaabf3:xyzzyaabg3),cbsderiv2_s(num_ao,&
&xyzzyaabf3:xyzzyaabg3),cbsderiv3_s(num_ao,xyzzyaabf3:xyzzyaabg3),cbsd&
&eriv4_s(num_ao,xyzzyaabf3:xyzzyaabg3),cbsderiv5_s(num_ao,xyzzyaabf3:x&
&yzzyaabg3),cbsderiv6_s(num_ao,xyzzyaabf3:xyzzyaabg3),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','cbsderiv1_s etc.')
endif
endif
if(use_backflow)then
allocate(cbsderiv1(num_ao,num_real_k_plus_1:num_k),cbsderiv2(num_ao,nu&
&m_real_k_plus_1:num_k),cbsderiv3(num_ao,num_real_k_plus_1:num_k),cbsd&
&eriv4(num_ao,num_real_k_plus_1:num_k),cbsderiv5(num_ao,num_real_k_plu&
&s_1:num_k),cbsderiv6(num_ao,num_real_k_plus_1:num_k),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','')
endif
else
allocate(cbf(0,0),cblap(0,0),cbgra1(0,0),cbgra2(0,0),cbgra3(0,0))
if(cusp_correction)allocate(cbf_s(0,0),cblap_s(0,0),cbgra1_s(0,0),cbgr&
&a2_s(0,0),cbgra3_s(0,0))
if(use_backflow)allocate(cbsderiv1(0,0),cbsderiv2(0,0),cbsderiv3(0,0),&
&cbsderiv4(0,0),cbsderiv5(0,0),cbsderiv6(0,0))
endif
if(excite)then
xyzzyaaaa3=maxval(xyzzyaaak1(1:nspin))
if(num_real_k>0)then
allocate(rck_ex(num_ao,xyzzyaaaa3,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','rck_ex')
endif
if(num_complex_k>0)then
allocate(cck_ex(num_ao,xyzzyaaaa3,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','cck_ex')
endif
do xyzzyaaaf3=1,nspin
xyzzyaaag3=which_ssingle(xyzzyaaaf3,spin_dep_gs)
xyzzyaaap3=which_ssingle(xyzzyaaaf3,spin_dep_in)
xyzzyaaaa3=num_ao*num_ao*(xyzzyaaap3-1)
xyzzyaaab3=2*xyzzyaaaa3
xyzzyaabe3=num_real_k*xyzzyaaaa3+num_complex_k*xyzzyaaab3
xyzzyaaaj3=0
do xyzzyaaav3=1,ndet
do xyzzyaaay3=1,wf_nd(xyzzyaaav3,xyzzyaaaf3)
xyzzyaaaz3=wf_d(1,xyzzyaaay3,xyzzyaaav3,xyzzyaaaf3)
xyzzyaaba3=wf_d(2,xyzzyaaay3,xyzzyaaav3,xyzzyaaaf3)
xyzzyaabb3=wf_d(3,xyzzyaaay3,xyzzyaaav3,xyzzyaaaf3)
xyzzyaabc3=wf_d(4,xyzzyaaay3,xyzzyaaav3,xyzzyaaaf3)
xyzzyaaak3=xyzzyaabe3
xyzzyaabd3=0
do xyzzyaaac3=1,num_real_k
do xyzzyaaah3=1,num_ao
if(xyzzyaaah3>nband(xyzzyaaac3,xyzzyaaag3))xyzzyaabd3=xyzzyaabd3+1
if(xyzzyaaac3==xyzzyaabc3.and.xyzzyaaah3==xyzzyaabb3)then
if(xyzzyaabt3(xyzzyaabd3,xyzzyaaaf3)/=0)then
xyzzyaabt3(xyzzyaabd3,xyzzyaaaf3)=0
xyzzyaaaj3=xyzzyaaaj3+1
rck_ex(1:num_ao,xyzzyaaaj3,xyzzyaaaf3)=xyzzyaaay1(xyzzyaaak3+1:xyzzyaa&
&ak3+num_ao)
endif
endif
xyzzyaaak3=xyzzyaaak3+num_ao
enddo
enddo
do xyzzyaaac3=num_real_k_plus_1,num_k
do xyzzyaaah3=1,num_ao
if(xyzzyaaah3>nband(xyzzyaaac3,xyzzyaaag3))xyzzyaabd3=xyzzyaabd3+1
if(xyzzyaaac3==xyzzyaabc3.and.xyzzyaaah3==xyzzyaabb3)then
if(xyzzyaabt3(xyzzyaabd3,xyzzyaaaf3)/=0)then
xyzzyaabt3(xyzzyaabd3,xyzzyaaaf3)=0
xyzzyaaaj3=xyzzyaaaj3+1
xyzzyaaab3=1
do xyzzyaaaa3=1,num_ao
cck_ex(xyzzyaaaa3,xyzzyaaaj3,xyzzyaaaf3)=cmplx(xyzzyaaay1(xyzzyaaak3+x&
&yzzyaaab3),xyzzyaaay1(xyzzyaaak3+xyzzyaaab3+1),kind=dp)
xyzzyaaab3=xyzzyaaab3+2
enddo
endif
endif
xyzzyaaak3=xyzzyaaak3+2*num_ao
enddo
enddo
enddo
enddo
enddo
deallocate(xyzzyaabt3)
xyzzyaaaa3=gauss_ex_norb
allocate(psi(xyzzyaaaa3),fgra1(xyzzyaaaa3),fgra2(xyzzyaaaa3),fgra3(xyz&
&zyaaaa3),flap(xyzzyaaaa3),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','psi etc.')
if(use_backflow)then
allocate(fsderiv1(xyzzyaaaa3),fsderiv2(xyzzyaaaa3),fsderiv3(xyzzyaaaa3&
&),fsderiv4(xyzzyaaaa3),fsderiv5(xyzzyaaaa3),fsderiv6(xyzzyaaaa3),stat&
&=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','fsderiv1 etc.')
endif
endif
call shallocate_barrier
call deshallocate(xyzzyaaay1)
if(.not.isperiodic.or.(isperiodic.and.cusp_correction))then
allocate(ao_m(num_ao),alap_m(num_ao),agra1_m(num_ao),agra2_m(num_ao),a&
&gra3_m(num_ao),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','ao_m etc.')
if(use_backflow)then
allocate(asderiv1_m(num_ao),asderiv2_m(num_ao),asderiv3_m(num_ao),asde&
&riv4_m(num_ao),asderiv5_m(num_ao),asderiv6_m(num_ao),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','asderiv1 etc.')
endif
endif
if(use_orbmods)call read_gwfmolorb
call xyzzyaacm1
allocate(xyzzyaaas1(xyzzyaaap1,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','gauss_orbmask')
call get_numerical_orbmask(xyzzyaaap1,nemax,nspin,ndet,xyzzyaaar1,xyzz&
&yaaas1)
if(isperiodic)then
allocate(xyzzyaaat1(xyzzyaaaq1,nspin),stat=xyzzyaaao3)
call check_alloc(xyzzyaaao3,'GWFDET_SETUP','gauss_corbmask')
xyzzyaaat1=.false.
do xyzzyaaaf3=1,nspin
do xyzzyaabh3=1,xyzzyaaap1
if(.not.xyzzyaaas1(xyzzyaabh3,xyzzyaaaf3))cycle
xyzzyaaat1(ridx2zidx(xyzzyaabh3),xyzzyaaaf3)=.true.
enddo
enddo
endif
if(have_ae)then
if(am_master)then
call wout('Gaussian cusp correction')
call wout('========================')
endif
if(cusp_correction)then
call xyzzyaaco1
call cusp_setup
if(am_master)call wout('Activated.')
else
if(am_master)call wout('Deactivated.')
endif
if(am_master)call wout()
endif
if(use_gpcc)call xyzzyaacn1
if(.not.use_gpcc.and..not.use_orbmods)deallocate(xyzzyaaas1)
end subroutine gwfdet_setup
subroutine xyzzyaacm1
implicit none
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4,xyzzyaaag4,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4,xyzzyaaak4
integer,allocatable :: xyzzyaaal4(:)
real(dp) xyzzyaaam4,xyzzyaaan4,xyzzyaaao4,xyzzyaaap4,xyzzyaaaq4(3),xyz&
&zyaaar4(3),xyzzyaaas4(3),xyzzyaaat4,xyzzyaaau4,xyzzyaaav4,xyzzyaaaw4,&
&xyzzyaaax4,xyzzyaaay4,xyzzyaaaz4,xyzzyaaba4,xyzzyaabb4,xyzzyaabc4,xyz&
&zyaabd4,xyzzyaabe4,xyzzyaabf4,xyzzyaabg4,xyzzyaabh4,xyzzyaabi4,xyzzya&
&abj4,xyzzyaabk4,xyzzyaabl4,xyzzyaabm4,xyzzyaabn4,xyzzyaabo4,xyzzyaabp&
&4,xyzzyaabq4,xyzzyaabr4,xyzzyaabs4,xyzzyaabt4,xyzzyaabu4,xyzzyaabv4,x&
&yzzyaabw4,xyzzyaabx4
real(dp) :: xyzzyaaby4=0.d0
real(dp),allocatable :: xyzzyaabz4(:),xyzzyaaca4(:),xyzzyaacb4(:),xyzz&
&yaacc4(:,:),xyzzyaacd4(:,:),xyzzyaace4(:,:),xyzzyaacf4(:,:),xyzzyaacg&
&4(:),xyzzyaach4(:,:),xyzzyaaci4(:,:,:),xyzzyaacj4(:,:,:)
logical xyzzyaack4
character(80) tmpr
screening_tolerance=abs(log(0.1d0)*dble(gautol))
if(isperiodic)then
allocate(xyzzyaacc4(3,num_g),num_centres_in_cell(num_g),num_shells_on_&
&centre(num_centres,num_g),xpos_in_cell(num_centres,num_g),ypos_in_cel&
&l(num_centres,num_g),zpos_in_cell(num_centres,num_g),shell_sequence_n&
&umber(num_g*num_shells),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','cell_coords,...(1)')
if(cusp_correction)then
allocate(which_ion(num_centres,num_g),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','which_ion')
endif
elseif(molgscreening)then
allocate(num_shells_on_centre(num_centres,ngsgrid3),shell_sequence_mol&
&(num_shells,ngsgrid3),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','num_shells_on_centre,.&
&..')
else
allocate(xyzzyaacc4(3,1),num_centres_in_cell(1),num_shells_on_centre(n&
&um_centres,1),xpos_in_cell(num_centres,1),ypos_in_cell(num_centres,1)&
&,zpos_in_cell(num_centres,1),shell_sequence_number(num_shells),stat=x&
&yzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','cell_coords,...(2)')
num_cell=1
xyzzyaacc4(1,1)=0.d0
xyzzyaacc4(2,1)=0.d0
xyzzyaacc4(3,1)=0.d0
num_centres_in_cell(1)=num_centres
do xyzzyaaab4=1,num_centres
num_shells_on_centre(xyzzyaaab4,1)=first_shell(xyzzyaaab4+1)-first_she&
&ll(xyzzyaaab4)
xyzzyaaai4=first_shell(xyzzyaaab4)
xpos_in_cell(xyzzyaaab4,1)=xyzzyaaaz1(1,xyzzyaaai4)
ypos_in_cell(xyzzyaaab4,1)=xyzzyaaaz1(2,xyzzyaaai4)
zpos_in_cell(xyzzyaaab4,1)=xyzzyaaaz1(3,xyzzyaaai4)
enddo
do xyzzyaaai4=1,num_shells
shell_sequence_number(xyzzyaaai4)=xyzzyaaai4
enddo
endif
if(isperiodic)then
if(am_master.and.printgscreening)then
call wout()
call wout('Gaussian screening')
call wout('------------------')
call wout('                  cell centre')
call wout('    cell    x(au)    y(au)    z(au)  # atoms with sig. Gaus&
&sians')
endif
allocate(xyzzyaaal4(num_centres),xyzzyaabz4(num_centres),xyzzyaaca4(nu&
&m_centres),xyzzyaacb4(num_centres),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','nshells,...')
select case(periodicity)
case(1)
xyzzyaaby4=0.5d0*pa1(1)
case(2)
allocate(xyzzyaacd4(4,2),xyzzyaacf4(4,2),xyzzyaacg4(4),stat=xyzzyaaac4&
&)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','2D')
xyzzyaaaq4(:)=0.5d0*pa1(:)
xyzzyaaar4(:)=0.5d0*pa2(:)
xyzzyaacd4(1,1:2)=-xyzzyaaaq4(1:2)-xyzzyaaar4(1:2)
xyzzyaacd4(2,1:2)=xyzzyaaaq4(1:2)-xyzzyaaar4(1:2)
xyzzyaacd4(3,1:2)=xyzzyaaaq4(1:2)+xyzzyaaar4(1:2)
xyzzyaacd4(4,1:2)=-xyzzyaaaq4(1:2)+xyzzyaaar4(1:2)
do xyzzyaaaf4=1,4
xyzzyaaag4=xyzzyaaaf4+1
if(xyzzyaaaf4==4)xyzzyaaag4=1
xyzzyaaat4=xyzzyaacd4(xyzzyaaag4,1)-xyzzyaacd4(xyzzyaaaf4,1)
xyzzyaaau4=xyzzyaacd4(xyzzyaaag4,2)-xyzzyaacd4(xyzzyaaaf4,2)
xyzzyaaam4=sqrt(xyzzyaaat4*xyzzyaaat4+xyzzyaaau4*xyzzyaaau4)
xyzzyaacg4(xyzzyaaaf4)=xyzzyaaam4
xyzzyaacf4(xyzzyaaaf4,1)=xyzzyaaat4/xyzzyaaam4
xyzzyaacf4(xyzzyaaaf4,2)=xyzzyaaau4/xyzzyaaam4
enddo
case(3)
allocate(xyzzyaace4(6,3),xyzzyaaci4(6,4,3),xyzzyaacj4(6,4,3),xyzzyaach&
&4(6,4),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','3D')
xyzzyaaaq4(:)=0.5d0*pa1(:)
xyzzyaaar4(:)=0.5d0*pa2(:)
xyzzyaaas4(:)=0.5d0*pa3(:)
xyzzyaaci4(1,1,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(1,2,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(1,3,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(1,4,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(2,1,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(2,2,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(2,3,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(2,4,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(3,1,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(3,2,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(3,3,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(3,4,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(4,1,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(4,2,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(4,3,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(4,4,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(5,1,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(5,2,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(5,3,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(5,4,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(6,1,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(6,2,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(6,3,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(6,4,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
do xyzzyaaaf4=1,6
do xyzzyaaag4=1,4
xyzzyaaah4=xyzzyaaag4+1
if(xyzzyaaag4==4)xyzzyaaah4=1
xyzzyaaat4=xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,1)-xyzzyaaci4(xyzzyaaaf4,x&
&yzzyaaag4,1)
xyzzyaaau4=xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,2)-xyzzyaaci4(xyzzyaaaf4,x&
&yzzyaaag4,2)
xyzzyaaav4=xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,3)-xyzzyaaci4(xyzzyaaaf4,x&
&yzzyaaag4,3)
xyzzyaaam4=sqrt(xyzzyaaat4*xyzzyaaat4+xyzzyaaau4*xyzzyaaau4+xyzzyaaav4&
&*xyzzyaaav4)
xyzzyaach4(xyzzyaaaf4,xyzzyaaag4)=xyzzyaaam4
xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)=xyzzyaaat4/xyzzyaaam4
xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)=xyzzyaaau4/xyzzyaaam4
xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,3)=xyzzyaaav4/xyzzyaaam4
enddo
enddo
xyzzyaace4(1,1)=a1(2)*a3(3)-a1(3)*a3(2)
xyzzyaace4(1,2)=-a1(1)*a3(3)+a1(3)*a3(1)
xyzzyaace4(1,3)=a1(1)*a3(2)-a1(2)*a3(1)
xyzzyaaam4=sqrt(xyzzyaace4(1,1)*xyzzyaace4(1,1)+xyzzyaace4(1,2)*xyzzya&
&ace4(1,2)+xyzzyaace4(1,3)*xyzzyaace4(1,3))
xyzzyaace4(1,:)=xyzzyaace4(1,:)/xyzzyaaam4
xyzzyaace4(2,1)=a1(2)*a2(3)-a1(3)*a2(2)
xyzzyaace4(2,2)=-a1(1)*a2(3)+a1(3)*a2(1)
xyzzyaace4(2,3)=a1(1)*a2(2)-a1(2)*a2(1)
xyzzyaaam4=sqrt(xyzzyaace4(2,1)*xyzzyaace4(2,1)+xyzzyaace4(2,2)*xyzzya&
&ace4(2,2)+xyzzyaace4(2,3)*xyzzyaace4(2,3))
xyzzyaace4(2,:)=xyzzyaace4(2,:)/xyzzyaaam4
xyzzyaace4(3,1)=a2(2)*a3(3)-a2(3)*a3(2)
xyzzyaace4(3,2)=-a2(1)*a3(3)+a2(3)*a3(1)
xyzzyaace4(3,3)=a2(1)*a3(2)-a2(2)*a3(1)
xyzzyaaam4=sqrt(xyzzyaace4(3,1)*xyzzyaace4(3,1)+xyzzyaace4(3,2)*xyzzya&
&ace4(3,2)+xyzzyaace4(3,3)*xyzzyaace4(3,3))
xyzzyaace4(3,:)=xyzzyaace4(3,:)/xyzzyaaam4
xyzzyaace4(4,:)=xyzzyaace4(3,:)
xyzzyaace4(5,:)=xyzzyaace4(1,:)
xyzzyaace4(6,:)=xyzzyaace4(2,:)
case default
call errstop_master('GWFDET_PRESCREEN','Confused about periodicity. Bu&
&g.')
end select
xyzzyaaaa4=0
num_cell=0
do xyzzyaaae4=1,num_g
xyzzyaaan4=p_lattice(1,xyzzyaaae4)
xyzzyaaao4=p_lattice(2,xyzzyaaae4)
xyzzyaaap4=p_lattice(3,xyzzyaaae4)
xyzzyaaaj4=0
do xyzzyaaab4=1,num_centres
xyzzyaaai4=first_shell(xyzzyaaab4)
xyzzyaaaw4=xyzzyaaaz1(1,xyzzyaaai4)
xyzzyaaax4=xyzzyaaaz1(2,xyzzyaaai4)
xyzzyaaay4=xyzzyaaaz1(3,xyzzyaaai4)
xyzzyaaaw4=xyzzyaaaw4+xyzzyaaan4
xyzzyaaax4=xyzzyaaax4+xyzzyaaao4
xyzzyaaay4=xyzzyaaay4+xyzzyaaap4
xyzzyaaba4=1.d100
select case (periodicity)
case(1)
if(xyzzyaaaw4>xyzzyaaby4)then
xyzzyaaba4=xyzzyaaaw4-xyzzyaaby4
elseif(xyzzyaaaw4<-xyzzyaaby4)then
xyzzyaaba4=-xyzzyaaby4-xyzzyaaaw4
else
xyzzyaaba4=0.d0
endif
case(2)
do xyzzyaaaf4=1,4
xyzzyaaag4=xyzzyaaaf4+1
if(xyzzyaaaf4==4)xyzzyaaag4=1
xyzzyaaat4=xyzzyaaaw4-xyzzyaacd4(xyzzyaaaf4,1)
xyzzyaaau4=xyzzyaaax4-xyzzyaacd4(xyzzyaaaf4,2)
xyzzyaaam4=abs(xyzzyaaat4*xyzzyaacf4(xyzzyaaaf4,2)-xyzzyaaau4*xyzzyaac&
&f4(xyzzyaaaf4,1))
xyzzyaabt4=abs(xyzzyaaat4*xyzzyaacf4(xyzzyaaaf4,1)+xyzzyaaau4*xyzzyaac&
&f4(xyzzyaaaf4,2))
xyzzyaaat4=xyzzyaaaw4-xyzzyaacd4(xyzzyaaag4,1)
xyzzyaaau4=xyzzyaaax4-xyzzyaacd4(xyzzyaaag4,2)
xyzzyaabu4=abs(xyzzyaaat4*xyzzyaacf4(xyzzyaaaf4,1)+xyzzyaaau4*xyzzyaac&
&f4(xyzzyaaaf4,2))
if(max(xyzzyaabt4,xyzzyaabu4)>xyzzyaacg4(xyzzyaaaf4))then
xyzzyaabv4=min(xyzzyaabt4,xyzzyaabu4)
xyzzyaabi4=sqrt(xyzzyaaam4*xyzzyaaam4+xyzzyaabv4*xyzzyaabv4)
else
xyzzyaabi4=xyzzyaaam4
endif
if(xyzzyaabi4<xyzzyaaba4)xyzzyaaba4=xyzzyaabi4
enddo
case(3)
do xyzzyaaaf4=1,6
xyzzyaaat4=xyzzyaaaw4-xyzzyaaci4(xyzzyaaaf4,1,1)
xyzzyaaau4=xyzzyaaax4-xyzzyaaci4(xyzzyaaaf4,1,2)
xyzzyaaav4=xyzzyaaay4-xyzzyaaci4(xyzzyaaaf4,1,3)
xyzzyaaam4=abs(xyzzyaaat4*xyzzyaace4(xyzzyaaaf4,1)+xyzzyaaau4*xyzzyaac&
&e4(xyzzyaaaf4,2)+xyzzyaaav4*xyzzyaace4(xyzzyaaaf4,3))
xyzzyaaat4=xyzzyaaaw4+xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,1)
xyzzyaaau4=xyzzyaaax4+xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,2)
xyzzyaaav4=xyzzyaaay4+xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,3)
xyzzyaabk4=xyzzyaaat4-xyzzyaaci4(xyzzyaaaf4,1,1)
xyzzyaabl4=xyzzyaaau4-xyzzyaaci4(xyzzyaaaf4,1,2)
xyzzyaabm4=xyzzyaaav4-xyzzyaaci4(xyzzyaaaf4,1,3)
xyzzyaabi4=abs(xyzzyaabk4*xyzzyaace4(xyzzyaaaf4,1)+xyzzyaabl4*xyzzyaac&
&e4(xyzzyaaaf4,2)+xyzzyaabm4*xyzzyaace4(xyzzyaaaf4,3))
if(xyzzyaabi4>xyzzyaaam4)then
xyzzyaaat4=xyzzyaaaw4-xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,1)
xyzzyaaau4=xyzzyaaax4-xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,2)
xyzzyaaav4=xyzzyaaay4-xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,3)
xyzzyaabk4=xyzzyaaat4-xyzzyaaci4(xyzzyaaaf4,1,1)
xyzzyaabl4=xyzzyaaau4-xyzzyaaci4(xyzzyaaaf4,1,2)
xyzzyaabm4=xyzzyaaav4-xyzzyaaci4(xyzzyaaaf4,1,3)
endif
xyzzyaabq4=xyzzyaabk4*painv(1,1)+xyzzyaabl4*painv(2,1)+xyzzyaabm4*pain&
&v(3,1)
xyzzyaabr4=xyzzyaabk4*painv(1,2)+xyzzyaabl4*painv(2,2)+xyzzyaabm4*pain&
&v(3,2)
xyzzyaabs4=xyzzyaabk4*painv(1,3)+xyzzyaabl4*painv(2,3)+xyzzyaabm4*pain&
&v(3,3)
xyzzyaabw4=1.d-8
if(xyzzyaabq4>=-xyzzyaabw4.and.xyzzyaabq4<=(1.d0+xyzzyaabw4).and.xyzzy&
&aabr4>=-xyzzyaabw4.and.xyzzyaabr4<=(1.d0+xyzzyaabw4).and.xyzzyaabs4>=&
&-xyzzyaabw4.and.xyzzyaabs4<=(1.d0+xyzzyaabw4))then
xyzzyaabi4=xyzzyaaam4
if(xyzzyaabi4<xyzzyaaba4)xyzzyaaba4=xyzzyaabi4
else
do xyzzyaaag4=1,4
xyzzyaaah4=xyzzyaaag4+1
if(xyzzyaaag4==4)xyzzyaaah4=1
xyzzyaabk4=xyzzyaaat4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaag4,1)
xyzzyaabl4=xyzzyaaau4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaag4,2)
xyzzyaabm4=xyzzyaaav4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaag4,3)
xyzzyaabn4=xyzzyaabl4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,3)-xyzzyaabm4*x&
&yzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)
xyzzyaabo4=-xyzzyaabk4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,3)+xyzzyaabm4*&
&xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)
xyzzyaabp4=xyzzyaabk4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)-xyzzyaabl4*x&
&yzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)
xyzzyaabj4=sqrt(xyzzyaabn4*xyzzyaabn4+xyzzyaabo4*xyzzyaabo4+xyzzyaabp4&
&*xyzzyaabp4)
xyzzyaabt4=abs(xyzzyaabk4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)+xyzzyaab&
&l4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)+xyzzyaabm4*xyzzyaacj4(xyzzyaaa&
&f4,xyzzyaaag4,3))
xyzzyaabk4=xyzzyaaat4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,1)
xyzzyaabl4=xyzzyaaau4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,2)
xyzzyaabm4=xyzzyaaav4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,3)
xyzzyaabu4=abs(xyzzyaabk4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)+xyzzyaab&
&l4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)+xyzzyaabm4*xyzzyaacj4(xyzzyaaa&
&f4,xyzzyaaag4,3))
if(max(xyzzyaabt4,xyzzyaabu4)>xyzzyaach4(xyzzyaaaf4,xyzzyaaag4))then
xyzzyaabv4=min(xyzzyaabt4,xyzzyaabu4)
xyzzyaabi4=sqrt(xyzzyaaam4*xyzzyaaam4+xyzzyaabj4*xyzzyaabj4+xyzzyaabv4&
&*xyzzyaabv4)
else
xyzzyaabi4=sqrt(xyzzyaaam4*xyzzyaaam4+xyzzyaabj4*xyzzyaabj4)
endif
if(xyzzyaabi4<xyzzyaaba4)xyzzyaaba4=xyzzyaabi4
enddo
endif
enddo
case default
continue
end select
xyzzyaaba4=xyzzyaaba4*xyzzyaaba4
xyzzyaaak4=0
do xyzzyaaai4=first_shell(xyzzyaaab4),first_shell(xyzzyaaab4+1)-1
xyzzyaaaz4=min_exponent(xyzzyaaai4)
xyzzyaabx4=xyzzyaaba4*xyzzyaaaz4
xyzzyaack4=.false.
if(xyzzyaabx4<=screening_tolerance)then
xyzzyaack4=.true.
else
if(periodicity==2)then
xyzzyaabq4=xyzzyaaaw4*painv(1,1)+xyzzyaaax4*painv(2,1)
xyzzyaabr4=xyzzyaaaw4*painv(1,2)+xyzzyaaax4*painv(2,2)
if(xyzzyaabq4>=-0.5d0.and.xyzzyaabq4<=0.5d0.and.xyzzyaabr4>=-0.5d0.and&
&.xyzzyaabr4<=0.5d0)xyzzyaack4=.true.
elseif(periodicity==3)then
xyzzyaabq4=xyzzyaaaw4*painv(1,1)+xyzzyaaax4*painv(2,1)+xyzzyaaay4*pain&
&v(3,1)
xyzzyaabr4=xyzzyaaaw4*painv(1,2)+xyzzyaaax4*painv(2,2)+xyzzyaaay4*pain&
&v(3,2)
xyzzyaabs4=xyzzyaaaw4*painv(1,3)+xyzzyaaax4*painv(2,3)+xyzzyaaay4*pain&
&v(3,3)
if(xyzzyaabq4>=-0.5d0.and.xyzzyaabq4<=0.5d0.and.xyzzyaabr4>=-0.5d0.and&
&.xyzzyaabr4<=0.5d0.and.xyzzyaabs4>=-0.5d0.and.xyzzyaabs4<=0.5d0)xyzzy&
&aack4=.true.
endif
endif
if(xyzzyaack4)then
xyzzyaaaa4=xyzzyaaaa4+1
xyzzyaaak4=xyzzyaaak4+1
shell_sequence_number(xyzzyaaaa4)=xyzzyaaai4
endif
enddo
if(xyzzyaaak4>0)then
xyzzyaaaj4=xyzzyaaaj4+1
xyzzyaaal4(xyzzyaaaj4)=xyzzyaaak4
xyzzyaabz4(xyzzyaaaj4)=xyzzyaaaw4
xyzzyaaca4(xyzzyaaaj4)=xyzzyaaax4
xyzzyaacb4(xyzzyaaaj4)=xyzzyaaay4
if(cusp_correction)which_ion(xyzzyaaaj4,num_cell+1)=xyzzyaaab4
endif
enddo
if(xyzzyaaaj4>0)then
num_cell=num_cell+1
xyzzyaacc4(1,num_cell)=xyzzyaaan4
xyzzyaacc4(2,num_cell)=xyzzyaaao4
xyzzyaacc4(3,num_cell)=xyzzyaaap4
num_centres_in_cell(num_cell)=xyzzyaaaj4
do xyzzyaaab4=1,xyzzyaaaj4
num_shells_on_centre(xyzzyaaab4,num_cell)=xyzzyaaal4(xyzzyaaab4)
xpos_in_cell(xyzzyaaab4,num_cell)=xyzzyaabz4(xyzzyaaab4)
ypos_in_cell(xyzzyaaab4,num_cell)=xyzzyaaca4(xyzzyaaab4)
zpos_in_cell(xyzzyaaab4,num_cell)=xyzzyaacb4(xyzzyaaab4)
enddo
if(am_master.and.printgscreening)then
write(tmpr,'(i3,i4,1x,3f9.3,1x,i5)')num_cell,xyzzyaaae4,xyzzyaaan4,xyz&
&zyaaao4,xyzzyaaap4,num_centres_in_cell(num_cell)
call wout(tmpr)
endif
endif
enddo
if(num_cell==0)call errstop_master('GWFDET_PRESCREEN','No cells contai&
&n significant shells.')
elseif(molgscreening)then
if(am_master.and.printgscreening)then
call wout()
call wout('Gaussian screening')
call wout('------------------')
endif
allocate(xyzzyaace4(6,3),xyzzyaaci4(6,4,3),xyzzyaacj4(6,4,3),xyzzyaach&
&4(6,4),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','3D <2>')
xyzzyaaaq4(:)=0.5d0*pa1(:)
xyzzyaaar4(:)=0.5d0*pa2(:)
xyzzyaaas4(:)=0.5d0*pa3(:)
xyzzyaaci4(1,1,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(1,2,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(1,3,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(1,4,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(2,1,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(2,2,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(2,3,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(2,4,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(3,1,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(3,2,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(3,3,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(3,4,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(4,1,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(4,2,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(4,3,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(4,4,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(5,1,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(5,2,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)-xyzzyaaas4(:)
xyzzyaaci4(5,3,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(5,4,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(6,1,:)=-xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(6,2,:)=xyzzyaaaq4(:)-xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(6,3,:)=xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
xyzzyaaci4(6,4,:)=-xyzzyaaaq4(:)+xyzzyaaar4(:)+xyzzyaaas4(:)
do xyzzyaaaf4=1,6
do xyzzyaaag4=1,4
xyzzyaaah4=xyzzyaaag4+1
if(xyzzyaaag4==4)xyzzyaaah4=1
xyzzyaaat4=xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,1)-xyzzyaaci4(xyzzyaaaf4,x&
&yzzyaaag4,1)
xyzzyaaau4=xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,2)-xyzzyaaci4(xyzzyaaaf4,x&
&yzzyaaag4,2)
xyzzyaaav4=xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,3)-xyzzyaaci4(xyzzyaaaf4,x&
&yzzyaaag4,3)
xyzzyaaam4=sqrt(xyzzyaaat4*xyzzyaaat4+xyzzyaaau4*xyzzyaaau4+xyzzyaaav4&
&*xyzzyaaav4)
xyzzyaach4(xyzzyaaaf4,xyzzyaaag4)=xyzzyaaam4
xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)=xyzzyaaat4/xyzzyaaam4
xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)=xyzzyaaau4/xyzzyaaam4
xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,3)=xyzzyaaav4/xyzzyaaam4
enddo
enddo
xyzzyaace4(1,1)=a1(2)*a3(3)-a1(3)*a3(2)
xyzzyaace4(1,2)=-a1(1)*a3(3)+a1(3)*a3(1)
xyzzyaace4(1,3)=a1(1)*a3(2)-a1(2)*a3(1)
xyzzyaaam4=sqrt(xyzzyaace4(1,1)*xyzzyaace4(1,1)+xyzzyaace4(1,2)*xyzzya&
&ace4(1,2)+ xyzzyaace4(1,3)*xyzzyaace4(1,3))
xyzzyaace4(1,:)=xyzzyaace4(1,:)/xyzzyaaam4
xyzzyaace4(2,1)=a1(2)*a2(3)-a1(3)*a2(2)
xyzzyaace4(2,2)=-a1(1)*a2(3)+a1(3)*a2(1)
xyzzyaace4(2,3)=a1(1)*a2(2)-a1(2)*a2(1)
xyzzyaaam4=sqrt(xyzzyaace4(2,1)*xyzzyaace4(2,1)+xyzzyaace4(2,2)*xyzzya&
&ace4(2,2)+ xyzzyaace4(2,3)*xyzzyaace4(2,3))
xyzzyaace4(2,:)=xyzzyaace4(2,:)/xyzzyaaam4
xyzzyaace4(3,1)=a2(2)*a3(3)-a2(3)*a3(2)
xyzzyaace4(3,2)=-a2(1)*a3(3)+a2(3)*a3(1)
xyzzyaace4(3,3)=a2(1)*a3(2)-a2(2)*a3(1)
xyzzyaaam4=sqrt(xyzzyaace4(3,1)*xyzzyaace4(3,1)+xyzzyaace4(3,2)*xyzzya&
&ace4(3,2)+ xyzzyaace4(3,3)*xyzzyaace4(3,3))
xyzzyaace4(3,:)=xyzzyaace4(3,:)/xyzzyaaam4
xyzzyaace4(4,:)=xyzzyaace4(3,:)
xyzzyaace4(5,:)=xyzzyaace4(1,:)
xyzzyaace4(6,:)=xyzzyaace4(2,:)
do xyzzyaaae4=1,ngsgrid3
xyzzyaaaa4=0
do xyzzyaaab4=1,num_centres
xyzzyaabe4=gsgrid(1,xyzzyaaae4)
xyzzyaabf4=gsgrid(2,xyzzyaaae4)
xyzzyaabg4=gsgrid(3,xyzzyaaae4)
xyzzyaaaw4=rion(1,xyzzyaaab4)-xyzzyaabe4
xyzzyaaax4=rion(2,xyzzyaaab4)-xyzzyaabf4
xyzzyaaay4=rion(3,xyzzyaaab4)-xyzzyaabg4
xyzzyaaba4=1.d100
do xyzzyaaaf4=1,6
xyzzyaaat4=xyzzyaaaw4-xyzzyaaci4(xyzzyaaaf4,1,1)
xyzzyaaau4=xyzzyaaax4-xyzzyaaci4(xyzzyaaaf4,1,2)
xyzzyaaav4=xyzzyaaay4-xyzzyaaci4(xyzzyaaaf4,1,3)
xyzzyaaam4=abs(xyzzyaaat4*xyzzyaace4(xyzzyaaaf4,1)+xyzzyaaau4*xyzzyaac&
&e4(xyzzyaaaf4,2)+xyzzyaaav4*xyzzyaace4(xyzzyaaaf4,3))
xyzzyaaat4=xyzzyaaaw4+xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,1)
xyzzyaaau4=xyzzyaaax4+xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,2)
xyzzyaaav4=xyzzyaaay4+xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,3)
xyzzyaabk4=xyzzyaaat4-xyzzyaaci4(xyzzyaaaf4,1,1)
xyzzyaabl4=xyzzyaaau4-xyzzyaaci4(xyzzyaaaf4,1,2)
xyzzyaabm4=xyzzyaaav4-xyzzyaaci4(xyzzyaaaf4,1,3)
xyzzyaabi4=abs(xyzzyaabk4*xyzzyaace4(xyzzyaaaf4,1)+xyzzyaabl4*xyzzyaac&
&e4(xyzzyaaaf4,2)+xyzzyaabm4*xyzzyaace4(xyzzyaaaf4,3))
if(xyzzyaabi4>xyzzyaaam4)then
xyzzyaaat4=xyzzyaaaw4-xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,1)
xyzzyaaau4=xyzzyaaax4-xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,2)
xyzzyaaav4=xyzzyaaay4-xyzzyaaam4*xyzzyaace4(xyzzyaaaf4,3)
xyzzyaabk4=xyzzyaaat4-xyzzyaaci4(xyzzyaaaf4,1,1)
xyzzyaabl4=xyzzyaaau4-xyzzyaaci4(xyzzyaaaf4,1,2)
xyzzyaabm4=xyzzyaaav4-xyzzyaaci4(xyzzyaaaf4,1,3)
endif
xyzzyaabq4=xyzzyaabk4*painv(1,1)+xyzzyaabl4*painv(2,1)+xyzzyaabm4*pain&
&v(3,1)
xyzzyaabr4=xyzzyaabk4*painv(1,2)+xyzzyaabl4*painv(2,2)+xyzzyaabm4*pain&
&v(3,2)
xyzzyaabs4=xyzzyaabk4*painv(1,3)+xyzzyaabl4*painv(2,3)+xyzzyaabm4*pain&
&v(3,3)
xyzzyaabw4=1.d-8
if(xyzzyaabq4>=-xyzzyaabw4.and.xyzzyaabq4<=(1.d0+xyzzyaabw4).and.xyzzy&
&aabr4>=-xyzzyaabw4.and.xyzzyaabr4<=(1.d0+xyzzyaabw4).and.xyzzyaabs4>=&
&-xyzzyaabw4.and.xyzzyaabs4<=(1.d0+xyzzyaabw4))then
xyzzyaabi4=xyzzyaaam4
if(xyzzyaabi4<xyzzyaaba4)xyzzyaaba4=xyzzyaabi4
else
do xyzzyaaag4=1,4
xyzzyaaah4=xyzzyaaag4+1
if(xyzzyaaag4==4)xyzzyaaah4=1
xyzzyaabk4=xyzzyaaat4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaag4,1)
xyzzyaabl4=xyzzyaaau4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaag4,2)
xyzzyaabm4=xyzzyaaav4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaag4,3)
xyzzyaabn4=xyzzyaabl4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,3)-xyzzyaabm4*x&
&yzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)
xyzzyaabo4=-xyzzyaabk4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,3)+xyzzyaabm4*&
&xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)
xyzzyaabp4=xyzzyaabk4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)-xyzzyaabl4*x&
&yzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)
xyzzyaabj4=sqrt(xyzzyaabn4*xyzzyaabn4+xyzzyaabo4*xyzzyaabo4+xyzzyaabp4&
&*xyzzyaabp4)
xyzzyaabt4=abs(xyzzyaabk4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)+xyzzyaab&
&l4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)+xyzzyaabm4*xyzzyaacj4(xyzzyaaa&
&f4,xyzzyaaag4,3))
xyzzyaabk4=xyzzyaaat4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,1)
xyzzyaabl4=xyzzyaaau4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,2)
xyzzyaabm4=xyzzyaaav4-xyzzyaaci4(xyzzyaaaf4,xyzzyaaah4,3)
xyzzyaabu4=abs(xyzzyaabk4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,1)+xyzzyaab&
&l4*xyzzyaacj4(xyzzyaaaf4,xyzzyaaag4,2)+xyzzyaabm4*xyzzyaacj4(xyzzyaaa&
&f4,xyzzyaaag4,3))
if(max(xyzzyaabt4,xyzzyaabu4)>xyzzyaach4(xyzzyaaaf4,xyzzyaaag4))then
xyzzyaabv4=min(xyzzyaabt4,xyzzyaabu4)
xyzzyaabi4=sqrt(xyzzyaaam4*xyzzyaaam4+xyzzyaabj4*xyzzyaabj4+xyzzyaabv4&
&*xyzzyaabv4)
else
xyzzyaabi4=sqrt(xyzzyaaam4*xyzzyaaam4+xyzzyaabj4*xyzzyaabj4)
endif
if(xyzzyaabi4<xyzzyaaba4)xyzzyaaba4=xyzzyaabi4
enddo
endif
enddo
xyzzyaaba4=xyzzyaaba4*xyzzyaaba4
xyzzyaaak4=0
do xyzzyaaai4=first_shell(xyzzyaaab4),first_shell(xyzzyaaab4+1)-1
xyzzyaaaz4=min_exponent(xyzzyaaai4)
xyzzyaabx4=xyzzyaaba4*xyzzyaaaz4
xyzzyaack4=.false.
if(xyzzyaabx4<=screening_tolerance)then
xyzzyaack4=.true.
else
xyzzyaabq4=xyzzyaaaw4*painv(1,1)+xyzzyaaax4*painv(2,1)+xyzzyaaay4*pain&
&v(3,1)
xyzzyaabr4=xyzzyaaaw4*painv(1,2)+xyzzyaaax4*painv(2,2)+xyzzyaaay4*pain&
&v(3,2)
xyzzyaabs4=xyzzyaaaw4*painv(1,3)+xyzzyaaax4*painv(2,3)+xyzzyaaay4*pain&
&v(3,3)
if(xyzzyaabq4>=-0.5d0.and.xyzzyaabq4<=0.5d0.and.xyzzyaabr4>=-0.5d0.and&
&.xyzzyaabr4<=0.5d0.and.xyzzyaabs4>=-0.5d0.and.xyzzyaabs4<=0.5d0)xyzzy&
&aack4=.true.
endif
if(xyzzyaack4)then
xyzzyaaaa4=xyzzyaaaa4+1
xyzzyaaak4=xyzzyaaak4+1
shell_sequence_mol(xyzzyaaaa4,xyzzyaaae4)=xyzzyaaai4
endif
enddo
num_shells_on_centre(xyzzyaaab4,xyzzyaaae4)=xyzzyaaak4
if(am_master.and.printgscreening)then
call wout(trim(i2s(xyzzyaaae4))//' '//trim(i2s(xyzzyaaab4))//' '//trim&
&(i2s(xyzzyaaak4)))
endif
enddo
enddo
deallocate(gsgrid)
endif
if(am_master.and.printgscreening)call wout()
if(periodicity==2)deallocate(xyzzyaacd4,xyzzyaacf4,xyzzyaacg4)
if(periodicity==3)deallocate(xyzzyaace4,xyzzyaaci4,xyzzyaacj4,xyzzyaac&
&h4)
if(isperiodic)then
if(num_real_k>0)then
allocate(coskdotg(num_real_k,num_cell),coskdotg_shift(num_real_k),stat&
&=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','')
endif
if(num_complex_k>0)then
allocate(expikdotg(num_real_k_plus_1:num_k,num_cell),expikdotg_shift(n&
&um_real_k_plus_1:num_k),stat=xyzzyaaac4)
call check_alloc(xyzzyaaac4,'GWFDET_PRESCREEN','')
endif
do xyzzyaaae4=1,num_cell
xyzzyaaan4=xyzzyaacc4(1,xyzzyaaae4)
xyzzyaaao4=xyzzyaacc4(2,xyzzyaaae4)
xyzzyaaap4=xyzzyaacc4(3,xyzzyaaae4)
do xyzzyaaad4=1,num_real_k
xyzzyaabb4=kvec(1,xyzzyaaad4)
xyzzyaabc4=kvec(2,xyzzyaaad4)
xyzzyaabd4=kvec(3,xyzzyaaad4)
xyzzyaabh4=xyzzyaabb4*xyzzyaaan4+xyzzyaabc4*xyzzyaaao4+xyzzyaabd4*xyzz&
&yaaap4
coskdotg(xyzzyaaad4,xyzzyaaae4)=cos(xyzzyaabh4)
enddo
do xyzzyaaad4=num_real_k_plus_1,num_k
xyzzyaabb4=kvec(1,xyzzyaaad4)
xyzzyaabc4=kvec(2,xyzzyaaad4)
xyzzyaabd4=kvec(3,xyzzyaaad4)
xyzzyaabh4=xyzzyaabb4*xyzzyaaan4+xyzzyaabc4*xyzzyaaao4+xyzzyaabd4*xyzz&
&yaaap4
expikdotg(xyzzyaaad4,xyzzyaaae4)=exp(xyzzyaabh4*zi)
enddo
enddo
endif
if(allocated(xyzzyaaal4))deallocate(xyzzyaaal4)
if(allocated(xyzzyaabz4))deallocate(xyzzyaabz4)
if(allocated(xyzzyaaca4))deallocate(xyzzyaaca4)
if(allocated(xyzzyaacb4))deallocate(xyzzyaacb4)
end subroutine xyzzyaacm1
subroutine xyzzyaacn1(iparam_buffer,prestore)
use slaarnabi,only : use_gpcc,initialize_gpcc,setup_gpcc,naeions_prim,&
&ae_index,spherical_av_real,spherical_av_cmplx,nradgrid,radgrid,nsphgr&
&id,sphgrid,gpcc_load_from_buffer,gpcc_save_to_buffer
implicit none
integer,intent(in),optional :: iparam_buffer
logical,intent(in),optional :: prestore
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5,xyzzyaaag5,xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5,xyzzyaaak5
real(dp) xyzzyaaal5(3)
real(dp),allocatable :: xyzzyaaam5(:,:),xyzzyaaan5(:,:,:),xyzzyaaao5(:&
&,:),xyzzyaaap5(:,:),xyzzyaaaq5(:,:)
complex(dp),allocatable :: xyzzyaaar5(:,:)
logical xyzzyaaas5
logical,save :: xyzzyaaat5=.true.
if(present(iparam_buffer).and.present(prestore))then
if(prestore)then
call gpcc_load_from_buffer(iparam_buffer)
return
endif
endif
if(complex_wf)call errstop_master('GWFDET_SETUP_GPCC','Cannot use GPCC&
& with Gaussian orbitals and a complex wave function.  Ask Neil if you&
& need this ability.  NB uncorrected Gaussian orbitals don''t yet work&
& with complex wave functions.')
if(xyzzyaaat5)then
if(isperiodic)then
call initialize_gpcc(spin_dep_full,xyzzyaaaq1,xyzzyaaat1,.false.)
else
call initialize_gpcc(spin_dep_full,xyzzyaaap1,xyzzyaaas1,.true.)
endif
endif
allocate(xyzzyaaam5(xyzzyaaap1,real1_complex2),xyzzyaaan5(3,xyzzyaaap1&
&,real1_complex2),xyzzyaaao5(xyzzyaaap1,real1_complex2),stat=xyzzyaaaa&
&5)
call check_alloc(xyzzyaaaa5,'GWFDET_SETUP_GPCC','valtemp')
xyzzyaaam5=0.d0
xyzzyaaan5=0.d0
xyzzyaaao5=0.d0
if(isperiodic)then
allocate(xyzzyaaar5(nsphgrid,xyzzyaaaq1),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'GWFDET_SETUP_GPCC','orb_sphgrid_c')
xyzzyaaar5=czero
else
allocate(xyzzyaaap5(nsphgrid,xyzzyaaap1),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'GWFDET_SETUP_GPCC','orb_sphgrid_r')
xyzzyaaap5=0.d0
endif
use_gpcc=.false.
do xyzzyaaab5=1,xyzzyaaaj1
do xyzzyaaae5=1,naeions_prim
xyzzyaaaf5=ae_index(xyzzyaaae5)
do xyzzyaaad5=1,nradgrid
if(xyzzyaaad5>1)then
xyzzyaaag5=nsphgrid
else
xyzzyaaag5=1
endif
do xyzzyaaac5=1,xyzzyaaag5
xyzzyaaal5=rion(1:3,xyzzyaaaf5)+radgrid(xyzzyaaad5,xyzzyaaae5)*sphgrid&
&(1:3,xyzzyaaac5)
call gaussian_orb_eval(xyzzyaaal5,xyzzyaaab5,xyzzyaaap1,xyzzyaaas1(1,x&
&yzzyaaab5),.true.,.false.,xyzzyaaam5,xyzzyaaan5,xyzzyaaao5)
do xyzzyaaaj5=1,xyzzyaaap1
if(.not.xyzzyaaas1(xyzzyaaaj5,xyzzyaaab5))cycle
if(isperiodic)then
xyzzyaaah5=xyzzyaaav1(ridx2zidx(xyzzyaaaj5))
xyzzyaaai5=xyzzyaaau1(ridx2zidx(xyzzyaaaj5))
if(any(xyzzyaaax1(:,xyzzyaaah5,xyzzyaaai5,xyzzyaaab5)==2))then
xyzzyaaas5=.false.
do xyzzyaaak5=xyzzyaaaj5+1,xyzzyaaap1
if(xyzzyaaas1(xyzzyaaak5,xyzzyaaab5).and.ridx2zidx(xyzzyaaaj5)==ridx2z&
&idx(xyzzyaaak5))then
xyzzyaaas5=.true.
exit
endif
enddo
if(xyzzyaaas5)then
xyzzyaaar5(xyzzyaaac5,ridx2zidx(xyzzyaaaj5))=cmplx(xyzzyaaam5(xyzzyaaa&
&j5,1),xyzzyaaam5(xyzzyaaak5,1),dp)
else
cycle
endif
elseif(any(xyzzyaaax1(:,xyzzyaaah5,xyzzyaaai5,xyzzyaaab5)==1))then
xyzzyaaar5(xyzzyaaac5,ridx2zidx(xyzzyaaaj5))=cmplx(xyzzyaaam5(xyzzyaaa&
&j5,1),0.d0,dp)
endif
else
xyzzyaaap5(xyzzyaaac5,xyzzyaaaj5)=xyzzyaaam5(xyzzyaaaj5,1)
endif
enddo
enddo
if(isperiodic)then
call spherical_av_cmplx(xyzzyaaab5,xyzzyaaae5,xyzzyaaad5,xyzzyaaar5)
else
call spherical_av_real(xyzzyaaab5,xyzzyaaae5,xyzzyaaad5,xyzzyaaap5)
endif
enddo
enddo
enddo
deallocate(xyzzyaaam5,xyzzyaaan5,xyzzyaaao5)
if(isperiodic)then
deallocate(xyzzyaaar5)
else
deallocate(xyzzyaaap5)
endif
use_gpcc=.true.
if(isperiodic.and.xyzzyaaat5)then
allocate(xyzzyaaaq5(3,xyzzyaaaq1),stat=xyzzyaaaa5)
call check_alloc(xyzzyaaaa5,'GWFDET_SETUP_GPCC','kvec_gpcc')
xyzzyaaaq5=0.d0
do xyzzyaaaj5=1,xyzzyaaaq1
xyzzyaaaq5(1:3,xyzzyaaaj5)=kvec(1:3,xyzzyaaau1(xyzzyaaaj5))
enddo
call setup_gpcc(xyzzyaaaq5)
deallocate(xyzzyaaaq5)
else
call setup_gpcc
endif
if(xyzzyaaat5)then
if(isperiodic)then
allocate(cusp_val_c(xyzzyaaaq1),cusp_grad_c(3,xyzzyaaaq1),cusp_lap_c(x&
&yzzyaaaq1),stat=xyzzyaaaa5)
else
allocate(cusp_val_r(xyzzyaaap1),cusp_grad_r(3,xyzzyaaap1),cusp_lap_r(x&
&yzzyaaap1),stat=xyzzyaaaa5)
endif
call check_alloc(xyzzyaaaa5,'GWFDET_SETUP_GPCC','cusp_val')
if(use_backflow)then
if(isperiodic)then
allocate(cusp_sderivs_c(6,xyzzyaaaq1),stat=xyzzyaaaa5)
else
allocate(cusp_sderivs_r(6,xyzzyaaap1),stat=xyzzyaaaa5)
endif
call check_alloc(xyzzyaaaa5,'GWFDET_SETUP_GPCC','cusp_sderivs')
endif
endif
if(present(iparam_buffer))call gpcc_save_to_buffer(iparam_buffer)
xyzzyaaat5=.false.
end subroutine xyzzyaacn1
subroutine gaussian_orb_eval(rvec,jspin,norb,orbmask,val,fsd,orbval,or&
&bgrad,orblap,orbsderivs)
use slaarnaaw,        only : gauss_mol_orb_eval
use slaarnaav,     only : gauss_mol_bf_orb_eval
use slaarnaba,        only : gauss_per_orb_eval
use slaarnaaz,   only : gauss_per_cusp_orb_eval
use slaarnaax,     only : gauss_per_bf_orb_eval
use slaarnaay,only : gauss_per_cusp_bf_orb_eval
implicit none
integer,intent(in) :: jspin,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd,orbmask(norb)
select case(xyzzyaaae1)
case(1)
if(present(orbsderivs).and.fsd)then
call gauss_mol_bf_orb_eval(rvec,jspin,norb,orbmask,val,orbval,orbgrad,&
&orbsderivs)
else
call gauss_mol_orb_eval(rvec,jspin,norb,orbmask,val,fsd,orbval,orbgrad&
&,orblap)
endif
case(2)
if(present(orbsderivs).and.fsd)then
call gauss_per_bf_orb_eval(rvec,jspin,norb,orbmask,xyzzyaaaq1,xyzzyaaa&
&t1(1,jspin),val,orbval,orbgrad,orblap,orbsderivs)
else
call gauss_per_orb_eval(rvec,jspin,norb,orbmask,xyzzyaaaq1,xyzzyaaat1(&
&1,jspin),val,fsd,orbval,orbgrad,orblap)
endif
case(3)
if(present(orbsderivs).and.fsd)then
call gauss_per_cusp_bf_orb_eval(rvec,jspin,norb,orbmask,xyzzyaaaq1,xyz&
&zyaaat1(1,jspin),val,orbval,orbgrad,orblap,orbsderivs)
else
call gauss_per_cusp_orb_eval(rvec,jspin,norb,orbmask,xyzzyaaaq1,xyzzya&
&aat1(1,jspin),val,fsd,orbval,orbgrad,orblap)
endif
case default
call errstop_master('GAUSSIAN_ORB_EVAL','Unknown gausstype.')
end select
end subroutine gaussian_orb_eval
subroutine xyzzyaaco1
implicit none
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7,xyzzyaa&
&af7
allocate(num_unique_orbitals(nspin),stat=xyzzyaaaa7)
call check_alloc(xyzzyaaaa7,'CUSP_ALLOC','num_unique_orbitals')
do xyzzyaaae7=1,nspin
if(nuc_nele(xyzzyaaae7)==0)then
num_unique_orbitals(xyzzyaaae7)=0
cycle
endif
xyzzyaaaf7=which_ssingle(xyzzyaaae7,spin_dep_gs)
xyzzyaaad7=0
do xyzzyaaac7=1,num_k
do xyzzyaaab7=1,nband(xyzzyaaac7,xyzzyaaaf7)
xyzzyaaad7=xyzzyaaad7+1
if(xyzzyaaac7>num_real_k.and.(xyzzyaaad7+1)<=nuc_nele(xyzzyaaae7))xyzz&
&yaaad7=xyzzyaaad7+1
enddo
enddo
num_unique_orbitals(xyzzyaaae7)=xyzzyaaad7+xyzzyaaak1(xyzzyaaae7)
enddo
xyzzyaabp1=maxval(num_unique_orbitals(:))
if(isperiodic.and.num_complex_k>0)then
xyzzyaaaf1=nitot
else
xyzzyaaaf1=nbasis
endif
allocate(rcusp_sq(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),acusp(xyzzyaabq1,x&
&yzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),disign(xyzzyaabp1,xyzzyaaaf1,xyzzyaa&
&aj1),pshift(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),stat=xyzzyaaaa7)
call check_alloc(xyzzyaaaa7,'CUSP_ALLOC','1')
if(isperiodic)then
allocate(exp_poly0(nemaxc),exp_poly1(nemaxc),exp_poly2(nemaxc),exp_pol&
&y3(nemaxc),exp_poly4(nemaxc),stat=xyzzyaaaa7)
call check_alloc(xyzzyaaaa7,'CUSP_ALLOC','exp_poly0')
exp_poly0=0.d0
exp_poly2=0.d0
exp_poly3=0.d0
exp_poly4=0.d0
if(use_backflow)then
allocate(exp_poly5(nemaxc),exp_poly6(nemaxc),exp_poly7(nemaxc),exp_pol&
&y8(nemaxc),exp_poly9(nemaxc),exp_poly10(nemaxc),stat=xyzzyaaaa7)
call check_alloc(xyzzyaaaa7,'CUSP_ALLOC','exp_poly5')
exp_poly5=0.d0
exp_poly6=0.d0
exp_poly7=0.d0
exp_poly8=0.d0
exp_poly9=0.d0
exp_poly10=0.d0
endif
endif
end subroutine xyzzyaaco1
subroutine cusp_setup
use slaarnaai,    only : wfdet_s
implicit none
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaa&
&af8,xyzzyaaag8,xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaak8,xyzzyaaal8&
&,xyzzyaaam8,xyzzyaaan8,xyzzyaaao8,xyzzyaaap8,xyzzyaaaq8,xyzzyaaar8,xy&
&zzyaaas8,xyzzyaaat8,xyzzyaaau8,xyzzyaaav8,xyzzyaaaw8,xyzzyaaax8,xyzzy&
&aaay8,xyzzyaaaz8,xyzzyaaba8,xyzzyaabb8,xyzzyaabc8,xyzzyaabd8,xyzzyaab&
&e8,xyzzyaabf8,xyzzyaabg8,xyzzyaabh8,xyzzyaabi8,xyzzyaabj8,xyzzyaabk8,&
&xyzzyaabl8
integer,parameter :: xyzzyaabm8=7
integer,allocatable,dimension(:) :: xyzzyaabn8,xyzzyaabo8
integer,allocatable,dimension(:,:) :: xyzzyaabp8,xyzzyaabq8
integer,allocatable,dimension(:,:,:) :: xyzzyaabr8,xyzzyaabs8
real(dp) xyzzyaabt8,xyzzyaabu8,xyzzyaabv8,xyzzyaabw8,xyzzyaabx8,xyzzya&
&aby8,xyzzyaabz8,xyzzyaaca8,xyzzyaacb8,xyzzyaacc8,xyzzyaacd8,xyzzyaace&
&8,xyzzyaacf8(xyzzyaabq1),xyzzyaacg8,xyzzyaach8,xyzzyaaci8,xyzzyaacj8,&
&xyzzyaack8,xyzzyaacl8,xyzzyaacm8,xyzzyaacn8,xyzzyaaco8,xyzzyaacp8,xyz&
&zyaacq8,xyzzyaacr8,xyzzyaacs8,xyzzyaact8,xyzzyaacu8,xyzzyaacv8,xyzzya&
&acw8,xyzzyaacx8,xyzzyaacy8(9),xyzzyaacz8(9),xyzzyaada8(0:7),xyzzyaadb&
&8,xyzzyaadc8,xyzzyaadd8,xyzzyaade8,xyzzyaadf8(0:xyzzyaabm8),xyzzyaadg&
&8(0:xyzzyaabm8),xyzzyaadh8(0:xyzzyaabm8),xyzzyaadi8(2),xyzzyaadj8(2),&
&xyzzyaadk8(2)
real(dp),parameter :: xyzzyaadl8=1.5d0,xyzzyaadm8=1.d-5,xyzzyaadn8=0.0&
&005d0,xyzzyaado8=1.d0/xyzzyaadn8,xyzzyaadp8=0.02d0,xyzzyaadq8=1.d-7
real(dp),allocatable,dimension(:) :: xyzzyaadr8,xyzzyaads8,xyzzyaadt8,&
&xyzzyaadu8,xyzzyaadv8,xyzzyaadw8,xyzzyaadx8,xyzzyaady8,xyzzyaadz8,xyz&
&zyaaea8,xyzzyaaeb8
real(dp),allocatable,dimension(:,:) :: xyzzyaaec8,xyzzyaaed8,xyzzyaaee&
&8
real(dp),allocatable,dimension(:,:,:) :: xyzzyaaef8,xyzzyaaeg8,xyzzyaa&
&eh8,xyzzyaaei8,xyzzyaaej8,xyzzyaaek8,xyzzyaael8,xyzzyaaem8,xyzzyaaen8&
&,xyzzyaaeo8
real(dp),allocatable,dimension(:,:,:,:) :: xyzzyaaep8,xyzzyaaeq8,xyzzy&
&aaer8
logical xyzzyaaes8,xyzzyaaet8,xyzzyaaeu8,xyzzyaaev8,xyzzyaaew8,xyzzyaa&
&ex8,xyzzyaaey8,xyzzyaaez8,xyzzyaafa8,xyzzyaafb8
logical,allocatable,dimension(:) :: xyzzyaafc8
character(80) tmpr,tmpr2
call timer('CUSPS',.true.)
call open_units(xyzzyaaah8,xyzzyaaal8)
if(xyzzyaaal8/=0)call errstop('CUSP_SETUP','Unable to find free i/o un&
&it.')
xyzzyaada8(0)=-0.560928d0
xyzzyaada8(1)=3.25819d0
xyzzyaada8(2)=-15.0126d0
xyzzyaada8(3)=33.7308d0
xyzzyaada8(4)=-42.8705d0
xyzzyaada8(5)=31.2276d0
xyzzyaada8(6)=-12.1316d0
xyzzyaada8(7)=1.94692d0
xyzzyaaci8=1.d0
xyzzyaaau8=int(xyzzyaaci8*xyzzyaado8)
xyzzyaabc8=nint(xyzzyaadp8/xyzzyaadn8)
xyzzyaade8=1.d0/cusp_control
xyzzyaaba8=0
xyzzyaaew8=cusp_info.and.am_master
xyzzyaaet8=.false.
xyzzyaaev8=.false.
xyzzyaaaq8=0
xyzzyaaar8=0
xyzzyaaas8=0
if(cusp_info)then
inquire(file='orbitals.in',exist=xyzzyaaev8)
if(xyzzyaaev8)then
open(xyzzyaaah8,file='orbitals.in',status='old',form='formatted',iosta&
&t=xyzzyaaal8)
if(xyzzyaaal8/=0)call errstop('CUSP_SETUP','Problem opening orbitals.i&
&n file')
read(xyzzyaaah8,*,iostat=xyzzyaaal8)xyzzyaaaq8,xyzzyaaar8,xyzzyaaas8
if(xyzzyaaas8<1.or.xyzzyaaas8>nspin.or.xyzzyaaaq8<1 .or.xyzzyaaar8<1.o&
&r.xyzzyaaaq8>num_unique_orbitals(xyzzyaaas8) .or.xyzzyaaar8>nitot.or.&
&xyzzyaaal8/=0)call errstop_master('CUSP_SETUP','Problem reading orbit&
&als.in file')
if(.not.is_ae(xyzzyaaar8))call errstop_master('CUSP_SETUP','Atom selec&
&ted for cusp plot in orbitals.in file is a pseudo-atom.')
close(xyzzyaaah8)
endif
endif
xyzzyaaeu8=.false.
xyzzyaafa8=.false.
if(xyzzyaaew8)then
call wout()
call wout('Verbose print out flagged (turn off with ''cusp_info : F'')&
&')
call wout()
if(xyzzyaaaj1==1)then
call wout('Spin restricted calculation.')
call wout()
endif
endif
if(use_orbmods)then
xyzzyaabk8=1
xyzzyaabl8=0
xyzzyaafb8=.false.
else
xyzzyaabk8=nnodes
xyzzyaabl8=my_node
xyzzyaafb8=nnodes>1
endif
allocate(xyzzyaaem8(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),xyzzyaaek8(xyzzy&
&aabp1,xyzzyaaaf1,xyzzyaaaj1),xyzzyaaer8(xyzzyaabq1,xyzzyaabp1,xyzzyaa&
&af1,xyzzyaaaj1),xyzzyaabn8(0:xyzzyaabk8),xyzzyaael8(xyzzyaabp1,xyzzya&
&aaf1,xyzzyaaaj1),xyzzyaadr8(xyzzyaaau8),xyzzyaafc8(xyzzyaaau8),xyzzya&
&ads8(xyzzyaaau8),xyzzyaaed8(xyzzyaaap1,real1_complex2),xyzzyaaeo8(3,x&
&yzzyaaap1,real1_complex2),xyzzyaaee8(xyzzyaaap1,real1_complex2),xyzzy&
&aabp8(xyzzyaabp1,xyzzyaaah1),xyzzyaabq8(xyzzyaabp1,xyzzyaaah1),stat=x&
&yzzyaabe8)
call check_alloc(xyzzyaabe8,'CUSP_SETUP','rcusp_sq etc.')
if(cusp_info)then
allocate(xyzzyaaef8(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),xyzzyaabr8(xyzzy&
&aabp1,xyzzyaaaf1,xyzzyaaaj1),xyzzyaaep8(xyzzyaabm8,xyzzyaabp1,xyzzyaa&
&af1,xyzzyaaaj1),xyzzyaaeh8(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),xyzzyaae&
&i8(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),xyzzyaaeg8(xyzzyaabp1,xyzzyaaaf1&
&,xyzzyaaaj1),xyzzyaaej8(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),xyzzyaaeq8(&
&xyzzyaabm8,xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),xyzzyaabs8(xyzzyaabp1,xy&
&zzyaaaf1,xyzzyaaaj1),xyzzyaaen8(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1),sta&
&t=xyzzyaabe8)
call check_alloc(xyzzyaabe8,'CUSP_SETUP','gauss0_s etc.')
xyzzyaabr8=0
endif
if(isperiodic)then
allocate(xyzzyaaec8(xyzzyaaaf1,xyzzyaaaj1),rcusp_sq_max(xyzzyaaaf1,xyz&
&zyaaaj1),stat=xyzzyaabe8)
call check_alloc(xyzzyaabe8,'CUSP_SETUP','rcusp_sq_max etc.')
xyzzyaaec8=0.d0
rcusp_sq_max=0.d0
endif
if(xyzzyaaev8)then
allocate(xyzzyaadw8(xyzzyaaau8),xyzzyaaea8(xyzzyaaau8),xyzzyaaeb8(xyzz&
&yaaau8),xyzzyaadx8(xyzzyaaau8),xyzzyaady8(xyzzyaaau8),xyzzyaadz8(xyzz&
&yaaau8),xyzzyaadt8(xyzzyaaau8),xyzzyaadu8(xyzzyaaau8),xyzzyaadv8(xyzz&
&yaaau8),stat=xyzzyaabe8)
endif
xyzzyaaek8=0.d0
xyzzyaaer8=0.d0
xyzzyaael8=0.d0
if(cusp_info)then
xyzzyaaef8=0.d0
xyzzyaaep8=0.d0
xyzzyaaeg8=0.d0
xyzzyaaeh8=0.d0
xyzzyaaei8=0.d0
xyzzyaaen8=0.d0
endif
if(xyzzyaabk8>1)then
acusp=0.d0
pshift=0.d0
allocate(xyzzyaabo8(xyzzyaabk8+1),stat=xyzzyaabe8)
call check_alloc(xyzzyaabe8,'CUSP_SETUP','is_per_node')
xyzzyaaag8=xyzzyaaaj1*xyzzyaaaf1
xyzzyaabo8(1:xyzzyaabk8)=xyzzyaaag8/xyzzyaabk8
do xyzzyaaaa8=1,mod(xyzzyaaag8,xyzzyaabk8)
xyzzyaabo8(xyzzyaaaa8)=xyzzyaabo8(xyzzyaaaa8)+1
enddo
xyzzyaaab8=1
xyzzyaabo8(xyzzyaabk8+1)=0
do xyzzyaaaa8=0,xyzzyaabk8
xyzzyaabn8(xyzzyaaaa8)=xyzzyaaab8
xyzzyaaab8=xyzzyaaab8+xyzzyaabo8(xyzzyaaaa8+1)
enddo
deallocate(xyzzyaabo8)
else
xyzzyaabn8(0)=1
xyzzyaabn8(1)=xyzzyaaaj1*xyzzyaaaf1+1
endif
xyzzyaabt8=0.d0
do xyzzyaaaa8=1,xyzzyaaau8
xyzzyaabt8=xyzzyaabt8+xyzzyaadn8
xyzzyaadr8(xyzzyaaaa8)=xyzzyaabt8
enddo
xyzzyaafc8(:)=.true.
do xyzzyaaao8=1,nspin
if(nuc_nele(xyzzyaaao8)==0)cycle
xyzzyaabf8=which_ssingle(xyzzyaaao8,spin_dep_gs)
xyzzyaaak8=0
do xyzzyaaac8=1,num_k
do xyzzyaaat8=1,nband(xyzzyaaac8,xyzzyaabf8)
xyzzyaaak8=xyzzyaaak8+1
xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)=xyzzyaaac8
xyzzyaabq8(xyzzyaaak8,xyzzyaabf8)=xyzzyaaat8
if(xyzzyaaac8>num_real_k.and.(xyzzyaaak8+1)<=nuc_nele(xyzzyaaao8))then
xyzzyaaak8=xyzzyaaak8+1
xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)=-xyzzyaaac8
xyzzyaabq8(xyzzyaaak8,xyzzyaabf8)=xyzzyaaat8
endif
enddo
enddo
enddo
xyzzyaaag8=0
rcusp_sq=0.d0
xyzzyaaem8=0.d0
xyzzyaabi8=0
do xyzzyaaao8=1,nspin
if(nuc_nele(xyzzyaaao8)==0)cycle
if(pmass(xyzzyaaao8)/=1.d0.or.pcharge(xyzzyaaao8)/=-1.d0)cycle
xyzzyaabh8=which_ssingle(xyzzyaaao8,spin_dep_full)
if(xyzzyaabh8==xyzzyaabi8)cycle
xyzzyaabi8=xyzzyaabh8
xyzzyaabf8=which_ssingle(xyzzyaaao8,spin_dep_gs)
xyzzyaabg8=which_ssingle(xyzzyaaao8,spin_dep_in)
do xyzzyaaaj8=1,xyzzyaaaf1
xyzzyaaai8=ion_prim(xyzzyaaaj8)
xyzzyaaag8=xyzzyaaag8+1
if(is_ae(xyzzyaaaj8).and.xyzzyaaag8>=xyzzyaabn8(xyzzyaabl8).and.xyzzya&
&aag8<xyzzyaabn8(xyzzyaabl8+1))then
xyzzyaacn8=zion(iontype_prim(xyzzyaaai8))
if(xyzzyaacn8==0.d0)cycle
xyzzyaadd8=1.d0/xyzzyaacn8
xyzzyaaez8=.false.
call gaussian_orb_eval(rion(1:3,xyzzyaaaj8),xyzzyaaao8,xyzzyaaap1,xyzz&
&yaaas1(1,xyzzyaaao8),.true.,.false.,xyzzyaaed8,xyzzyaaeo8,xyzzyaaee8)
do xyzzyaaak8=1,num_unique_orbitals(xyzzyaaao8)
if(xyzzyaaak8>nuc_nele(xyzzyaaao8))then
xyzzyaaez8=.true.
xyzzyaabd8=xyzzyaaak8-nuc_nele(xyzzyaaao8)
xyzzyaabj8=vrt2full(xyzzyaaak8-nuc_nele(xyzzyaaao8),xyzzyaaao8)
else
xyzzyaabj8=xyzzyaaak8+gauss_offset(xyzzyaaao8)-1
endif
if(.not.xyzzyaaas1(xyzzyaabj8,xyzzyaaao8))cycle
if(.not.xyzzyaaez8)then
call wfdet_s(0.d0,0,1.d0,xyzzyaaak8,xyzzyaaai8,xyzzyaaao8,xyzzyaabq8(x&
&yzzyaaak8,xyzzyaabf8),iabs(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)),xyzzyaa&
&aj8,.false.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
if(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)>0)then
xyzzyaach8=xyzzyaadi8(1)
else
xyzzyaach8=xyzzyaadi8(2)
endif
else
call wfdet_s(0.d0,0,1.d0,xyzzyaabd8,xyzzyaaai8,xyzzyaaao8,xyzzyaabd8,v&
&irtual_k(vrt2full(xyzzyaabd8,xyzzyaaao8),xyzzyaaao8),xyzzyaaaj8,.true&
&.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
xyzzyaach8=xyzzyaadi8(1)
endif
if(cusp_info)xyzzyaaef8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaach8
if(abs(xyzzyaach8)<xyzzyaadq8)cycle
xyzzyaadf8=0.d0
xyzzyaaci8=1.d0
if(xyzzyaaci8>=nearest_ion_r(xyzzyaaai8))xyzzyaaci8=0.9d0*nearest_ion_&
&r(xyzzyaaai8)
xyzzyaaav8=int(xyzzyaaci8*xyzzyaado8)
xyzzyaads8(1:xyzzyaaav8)=0.d0
xyzzyaaeu8=(xyzzyaaev8.and.xyzzyaaak8==xyzzyaaaq8.and.xyzzyaaaj8==xyzz&
&yaaar8.and.xyzzyaaao8==xyzzyaaas8)
xyzzyaaay8=int(xyzzyaaci8*xyzzyaadd8*xyzzyaado8)
xyzzyaaca8=xyzzyaadr8(xyzzyaaay8)
if(xyzzyaaba8>0)xyzzyaafc8(1:xyzzyaaay8)=.true.
xyzzyaaey8=.false.
if(.not.xyzzyaaez8)then
xyzzyaacj8=xyzzyaaed8(xyzzyaabj8,1)
else
xyzzyaacj8=xyzzyaaed8(xyzzyaaar1(xyzzyaaan1(xyzzyaabd8,xyzzyaaao8),xyz&
&zyaaao8,xyzzyaaao1(xyzzyaabd8,xyzzyaaao8)),1)
endif
if(cusp_info)xyzzyaaen8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaacj8
xyzzyaact8=xyzzyaacj8-xyzzyaach8
xyzzyaacp8=sign(1.d0,xyzzyaach8)
xyzzyaaek8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaacp8
if(xyzzyaaca8>xyzzyaadl8)then
call errwarn('CUSP_SETUP','Outside range for ideal polynomial curve fo&
&r orbital '//trim(i2s(xyzzyaaak8))//' on node '//trim(i2s(my_node))//&
&'.')
tmpr=r2s(xyzzyaaca8,'(f12.5)')
tmpr2=r2s(xyzzyaadl8,'(f12.5)')
call wout('rcmax = '//trim(tmpr)//', range = '//trim(tmpr2))
endif
xyzzyaaba8=0
xyzzyaaex8=.false.
xyzzyaacd8=0.d0
xyzzyaace8=(xyzzyaacn8**2)*xyzzyaade8
call xyzzyaafg8(xyzzyaadr8(xyzzyaaav8),xyzzyaacc8,xyzzyaacn8)
xyzzyaads8(xyzzyaaav8)=xyzzyaack8
call xyzzyaafg8(xyzzyaaca8,xyzzyaabz8,xyzzyaacn8)
do xyzzyaaaa8=xyzzyaaav8-1,1,-1
xyzzyaabt8=xyzzyaadr8(xyzzyaaaa8)
if(xyzzyaaaa8<=xyzzyaaay8.and..not.xyzzyaaex8)then
if(xyzzyaacn8/=1.d0)then
xyzzyaadb8=xyzzyaabt8
xyzzyaadc8=xyzzyaaca8
xyzzyaacb8=0.d0
do xyzzyaaab8=1,7
xyzzyaadb8=xyzzyaadb8*xyzzyaabt8
xyzzyaadc8=xyzzyaadc8*xyzzyaaca8
xyzzyaacb8=xyzzyaacb8+xyzzyaada8(xyzzyaaab8)*(xyzzyaadb8-xyzzyaadc8)
enddo
xyzzyaacb8=xyzzyaacb8*xyzzyaacn8**2+xyzzyaabz8
call xyzzyaafg8(xyzzyaabt8,xyzzyaacc8,xyzzyaacn8)
if(abs(xyzzyaacc8-xyzzyaacb8)>xyzzyaace8)then
xyzzyaaaz8=xyzzyaaaa8
xyzzyaacd8=xyzzyaabt8
xyzzyaaex8=.true.
endif
else
call xyzzyaafg8(xyzzyaabt8,xyzzyaacc8,xyzzyaacn8)
if(abs(xyzzyaacc8-xyzzyaabz8)>xyzzyaace8)then
xyzzyaaaz8=xyzzyaaaa8
xyzzyaacd8=xyzzyaabt8
xyzzyaaex8=.true.
endif
endif
else
if(.not.xyzzyaaez8)then
call wfdet_s(xyzzyaabt8,0,1.d0,xyzzyaaak8,xyzzyaaai8,xyzzyaaao8,xyzzya&
&abq8(xyzzyaaak8,xyzzyaabf8),iabs(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)),x&
&yzzyaaaj8,.false.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
if(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)>0)then
xyzzyaack8=xyzzyaadi8(1)
else
xyzzyaack8=xyzzyaadi8(2)
endif
else
call wfdet_s(xyzzyaabt8,0,1.d0,xyzzyaabd8,xyzzyaaai8,xyzzyaaao8,xyzzya&
&abd8,virtual_k(vrt2full(xyzzyaabd8,xyzzyaaao8),xyzzyaaao8),xyzzyaaaj8&
&,.true.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
xyzzyaack8=xyzzyaadi8(1)
endif
endif
xyzzyaads8(xyzzyaaaa8)=xyzzyaack8
if(xyzzyaads8(xyzzyaaaa8+1)*xyzzyaads8(xyzzyaaaa8)<0.d0)then
xyzzyaaba8=xyzzyaaba8+1
xyzzyaabu8=abs(xyzzyaabt8+xyzzyaadn8)
call xyzzyaafe8(xyzzyaabu8,xyzzyaabt8,xyzzyaadf8(xyzzyaaba8))
endif
enddo
if(xyzzyaacd8<1.d-10)then
call errstop_master('CUSP_SETUP','At least one of the orbital cusp rad&
&ii is effectively zero. Should not happen.')
endif
if(xyzzyaaba8>xyzzyaabm8-1)call errstop_master('CUSP_SETUP','Exceeded &
&maximum number of orbital nodes. Very unlikely to happen. Increase ma&
&xnodes parameter.')
if(xyzzyaaba8>0)then
xyzzyaabb8=0
xyzzyaadf8(xyzzyaaba8+1)=0.d0
xyzzyaadh8(xyzzyaaba8+1)=xyzzyaadn8+1.d-10
do xyzzyaaab8=xyzzyaaba8,1,-1
xyzzyaaaa8=int(xyzzyaadf8(xyzzyaaab8)*xyzzyaado8)
xyzzyaaae8=max(xyzzyaaaa8-xyzzyaabc8,1)
xyzzyaaaf8=min(xyzzyaaaa8+xyzzyaabc8,xyzzyaaav8)
xyzzyaadg8(xyzzyaaab8)=xyzzyaadr8(xyzzyaaae8)
xyzzyaadh8(xyzzyaaab8)=xyzzyaadr8(xyzzyaaaf8)
xyzzyaafc8(xyzzyaaae8:xyzzyaaaf8)=.false.
if(xyzzyaadg8(xyzzyaaab8)<xyzzyaaca8)xyzzyaabb8=xyzzyaaab8
enddo
if(.not.xyzzyaafc8(xyzzyaaaz8))then
xyzzyaaae8=xyzzyaaaz8
do xyzzyaaaa8=xyzzyaaae8,1,-1
if(xyzzyaafc8(xyzzyaaaa8))then
xyzzyaaaz8=xyzzyaaaa8
xyzzyaacd8=xyzzyaadr8(xyzzyaaaa8)
exit
endif
enddo
xyzzyaaey8=.true.
endif
if(xyzzyaadf8(xyzzyaabb8)>0.25d0.or.xyzzyaaey8)then
xyzzyaabt8=xyzzyaadg8(xyzzyaabb8)
call xyzzyaafg8(xyzzyaabt8,xyzzyaacw8,xyzzyaacn8)
if(xyzzyaacn8/=1.d0)then
xyzzyaadb8=xyzzyaabt8
xyzzyaadc8=xyzzyaaca8
xyzzyaacx8=0.d0
do xyzzyaaab8=1,7
xyzzyaadb8=xyzzyaadb8*xyzzyaabt8
xyzzyaadc8=xyzzyaadc8*xyzzyaaca8
xyzzyaacx8=xyzzyaacx8+xyzzyaada8(xyzzyaaab8)*(xyzzyaadb8-xyzzyaadc8)
enddo
xyzzyaacx8=xyzzyaacx8*xyzzyaacn8**2+xyzzyaabz8
else
xyzzyaacx8=xyzzyaabz8
endif
do xyzzyaaab8=xyzzyaabb8,xyzzyaaba8
xyzzyaaaf8=int(xyzzyaadg8(xyzzyaaab8)*xyzzyaado8)-1
xyzzyaaae8=int(xyzzyaadh8(xyzzyaaab8+1)*xyzzyaado8)
do xyzzyaaaa8=xyzzyaaaf8,xyzzyaaae8,-1
xyzzyaabt8=xyzzyaadr8(xyzzyaaaa8)
call xyzzyaafg8(xyzzyaabt8,xyzzyaacc8,xyzzyaacn8)
if(xyzzyaacn8/=1.d0)then
xyzzyaadb8=xyzzyaabt8
xyzzyaadc8=xyzzyaaca8
xyzzyaacb8=0.d0
do xyzzyaaad8=1,7
xyzzyaadb8=xyzzyaadb8*xyzzyaabt8
xyzzyaadc8=xyzzyaadc8*xyzzyaaca8
xyzzyaacb8=xyzzyaacb8+xyzzyaada8(xyzzyaaad8)*(xyzzyaadb8-xyzzyaadc8)
enddo
xyzzyaacb8=xyzzyaacb8*xyzzyaacn8**2+xyzzyaabz8
else
xyzzyaacb8=xyzzyaabz8
endif
if(abs(xyzzyaacc8-xyzzyaacb8)<abs(xyzzyaacw8-xyzzyaacx8))then
xyzzyaacw8=xyzzyaacc8
xyzzyaacx8=xyzzyaacb8
cycle
else
xyzzyaaca8=xyzzyaabt8
xyzzyaaay8=xyzzyaaaa8
exit
endif
enddo
call xyzzyaafg8(xyzzyaaca8,xyzzyaabz8,xyzzyaacn8)
do xyzzyaaaa8=xyzzyaaay8,xyzzyaaae8,-1
xyzzyaabt8=xyzzyaadr8(xyzzyaaaa8)
call xyzzyaafg8(xyzzyaabt8,xyzzyaacc8,xyzzyaacn8)
if(xyzzyaacn8/=1.d0)then
xyzzyaadb8=xyzzyaabt8
xyzzyaadc8=xyzzyaaca8
xyzzyaacb8=0.d0
do xyzzyaaad8=1,7
xyzzyaadb8=xyzzyaadb8*xyzzyaabt8
xyzzyaadc8=xyzzyaadc8*xyzzyaaca8
xyzzyaacb8=xyzzyaacb8+xyzzyaada8(xyzzyaaad8)*(xyzzyaadb8-xyzzyaadc8)
enddo
xyzzyaacb8=xyzzyaacb8*xyzzyaacn8**2+xyzzyaabz8
call xyzzyaafg8(xyzzyaabt8,xyzzyaacc8,xyzzyaacn8)
if(abs(xyzzyaacc8-xyzzyaacb8)>xyzzyaace8)then
xyzzyaaaz8=xyzzyaaaa8
xyzzyaacd8=xyzzyaabt8
exit
endif
else
call xyzzyaafg8(xyzzyaabt8,xyzzyaacc8,xyzzyaacn8)
if(abs(xyzzyaacc8-xyzzyaabz8)>xyzzyaace8)then
xyzzyaaaz8=xyzzyaaaa8
xyzzyaacd8=xyzzyaabt8
exit
endif
endif
enddo
enddo
endif
if(any(xyzzyaadf8(1:xyzzyaaba8)<xyzzyaacd8))then
if(xyzzyaacp8>0)then
xyzzyaaco8=2.d0*minval(xyzzyaads8(1:xyzzyaaav8))
else
xyzzyaaco8=2.d0*maxval(xyzzyaads8(1:xyzzyaaav8))
endif
else
xyzzyaaco8=0.d0
endif
else
xyzzyaaco8=0.d0
endif
xyzzyaael8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaaco8
xyzzyaacz8(1:9)=1.d100
xyzzyaaaf8=int(xyzzyaacd8*0.05d0*xyzzyaado8)
do xyzzyaaax8=1,9
xyzzyaaaa8=xyzzyaaaz8+(xyzzyaaax8-5)*xyzzyaaaf8
if(xyzzyaaaa8<0.or.xyzzyaaaa8>xyzzyaaav8)cycle
if(xyzzyaaax8>5.and..not.xyzzyaafc8(xyzzyaaaa8))exit
if(.not.xyzzyaafc8(xyzzyaaaa8))cycle
xyzzyaacd8=xyzzyaadr8(xyzzyaaaa8)
call xyzzyaafg8(xyzzyaacd8,xyzzyaaby8,xyzzyaacn8)
xyzzyaacu8=xyzzyaach8
call xyzzyaafd8(xyzzyaacu8,xyzzyaacv8)
xyzzyaacy8(xyzzyaaax8)=xyzzyaacu8
xyzzyaacz8(xyzzyaaax8)=xyzzyaacv8
enddo
if(all(xyzzyaacz8(:)==1.d100))call errstop_master('CUSP_SETUP','Proble&
&m finding optimum cusp radius : all rcusp too near node.')
xyzzyaaaw8=minloc(xyzzyaacz8(:),1)
xyzzyaaaa8=xyzzyaaaz8+(xyzzyaaaw8-5)*xyzzyaaaf8
xyzzyaaaz8=xyzzyaaaa8
xyzzyaacd8=xyzzyaadr8(xyzzyaaaa8)
xyzzyaaem8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaacd8*xyzzyaacd8
call xyzzyaafg8(xyzzyaacd8,xyzzyaaby8,xyzzyaacn8)
call xyzzyaaff8(xyzzyaacy8(xyzzyaaaw8),xyzzyaacv8)
if(cusp_info)then
xyzzyaaeh8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaacy8(xyzzyaaaw8)
xyzzyaaei8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaacz8(xyzzyaaaw8)
xyzzyaaeg8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaabr1
if(xyzzyaaba8>0)then
xyzzyaaab8=0
do xyzzyaaaa8=1,xyzzyaaba8
if(xyzzyaadf8(xyzzyaaaa8)<xyzzyaacd8)then
xyzzyaaab8=xyzzyaaab8+1
xyzzyaaep8(xyzzyaaab8,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaadf8(xyz&
&zyaaaa8)
endif
enddo
xyzzyaabr8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaaab8
endif
endif
if(xyzzyaaeu8)call xyzzyaafl8
if(xyzzyaaeu8)xyzzyaafa8=.true.
enddo
if(isperiodic)xyzzyaaec8(xyzzyaaaj8,xyzzyaabh8) =maxval(xyzzyaaem8(1:n&
&uc_nele(xyzzyaaao8),xyzzyaaaj8,xyzzyaabh8))
endif
enddo
enddo
if(xyzzyaaev8.and.xyzzyaafa8)then
xyzzyaacd8=sqrt(xyzzyaaem8(xyzzyaaaq8,xyzzyaaar8,which_ssingle(xyzzyaa&
&as8,spin_dep_full)))
xyzzyaaaz8=nint(xyzzyaacd8*xyzzyaado8)
open(xyzzyaaah8,file='local_energy.dat',status='unknown')
do xyzzyaaaa8=1,xyzzyaaaz8
write(xyzzyaaah8,'(1x,e20.10,1x,e20.10,1x,e20.10,1x,e20.10)')xyzzyaadr&
&8(xyzzyaaaa8),xyzzyaadu8(xyzzyaaaa8),xyzzyaadt8(xyzzyaaaa8),xyzzyaadv&
&8(xyzzyaaaa8)
enddo
write(xyzzyaaah8,*)'&'
close(xyzzyaaah8)
open(xyzzyaaah8,file='orbitals.dat',status='unknown',position='append'&
&)
do xyzzyaaaa8=1,xyzzyaaaz8
write(xyzzyaaah8,'(1x,e20.10,1x,e20.10,1x,e20.10)')xyzzyaadr8(xyzzyaaa&
&a8),xyzzyaadx8(xyzzyaaaa8),xyzzyaadw8(xyzzyaaaa8)
enddo
write(xyzzyaaah8,*)'&'
close(xyzzyaaah8)
open(xyzzyaaah8,file='gradients.dat',status='unknown',position='append&
&')
do xyzzyaaaa8=1,xyzzyaaaz8
write(xyzzyaaah8,'(1x,e20.10,1x,e20.10,1x,e20.10)')xyzzyaadr8(xyzzyaaa&
&a8),xyzzyaady8(xyzzyaaaa8),xyzzyaaea8(xyzzyaaaa8)
enddo
write(xyzzyaaah8,*)'&'
close(xyzzyaaah8)
open(xyzzyaaah8,file='laplacians.dat',status='unknown',position='appen&
&d')
do xyzzyaaaa8=1,xyzzyaaaz8
write(xyzzyaaah8,'(1x,e20.10,1x,e20.10,1x,e20.10)')xyzzyaadr8(xyzzyaaa&
&a8),xyzzyaadz8(xyzzyaaaa8),xyzzyaaeb8(xyzzyaaaa8)
enddo
write(xyzzyaaah8,*)'&'
close(xyzzyaaah8)
endif
if(xyzzyaafb8)then
call mpi_allreduce(xyzzyaaem8,rcusp_sq,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj&
&1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting rcusp_sq in cusp_setup')
call mpi_allreduce(xyzzyaaek8,disign,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj1,&
&mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting disign in cusp_setup')
call mpi_allreduce(xyzzyaaer8,acusp,xyzzyaabq1*xyzzyaabp1*xyzzyaaaf1*x&
&yzzyaaaj1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting acusp in cusp_setup')
call mpi_allreduce(xyzzyaael8,pshift,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj1,&
&mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pshift in cusp_setup')
deallocate(xyzzyaaem8,xyzzyaaek8,xyzzyaaer8,xyzzyaael8)
if(isperiodic)then
call mpi_allreduce(xyzzyaaec8,rcusp_sq_max,xyzzyaaaf1*xyzzyaaaj1,mpi_d&
&ouble_precision,mpi_sum,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting rcusp_sq_max in cusp_setup')
deallocate(xyzzyaaec8)
endif
else
rcusp_sq=xyzzyaaem8
disign=xyzzyaaek8
acusp=xyzzyaaer8
pshift=xyzzyaael8
deallocate(xyzzyaaem8,xyzzyaaek8,xyzzyaaer8,xyzzyaael8)
if(isperiodic)then
rcusp_sq_max=xyzzyaaec8
deallocate(xyzzyaaec8)
endif
endif
if(cusp_info)then
if(xyzzyaabk8>1.and.xyzzyaaag8>1)then
call mpi_reduce(xyzzyaaef8,xyzzyaaej8,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj1&
&,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
if(am_master)xyzzyaaef8(:,:,1:xyzzyaaaj1)=xyzzyaaej8(:,:,1:xyzzyaaaj1)
call mpi_reduce(xyzzyaaen8,xyzzyaaej8,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj1&
&,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
if(am_master)xyzzyaaen8(:,:,1:xyzzyaaaj1)=xyzzyaaej8(:,:,1:xyzzyaaaj1)
call mpi_reduce(xyzzyaaeg8,xyzzyaaej8,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj1&
&,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
if(am_master)xyzzyaaeg8(:,:,1:xyzzyaaaj1)=xyzzyaaej8(:,:,1:xyzzyaaaj1)
call mpi_reduce(xyzzyaabr8,xyzzyaabs8,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj1&
&,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
if(am_master)xyzzyaabr8(:,:,1:xyzzyaaaj1)=xyzzyaabs8(:,:,1:xyzzyaaaj1)
call mpi_reduce(xyzzyaaep8,xyzzyaaeq8,xyzzyaabm8*xyzzyaabp1*xyzzyaaaf1&
&*xyzzyaaaj1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
if(am_master)xyzzyaaep8(:,:,:,1:xyzzyaaaj1)=xyzzyaaeq8(:,:,:,1:xyzzyaa&
&aj1)
call mpi_reduce(xyzzyaaeh8,xyzzyaaej8,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj1&
&,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
if(am_master)xyzzyaaeh8(:,:,1:xyzzyaaaj1)=xyzzyaaej8(:,:,1:xyzzyaaaj1)
call mpi_reduce(xyzzyaaei8,xyzzyaaej8,xyzzyaabp1*xyzzyaaaf1*xyzzyaaaj1&
&,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
if(am_master)xyzzyaaei8(:,:,1:xyzzyaaaj1)=xyzzyaaej8(:,:,1:xyzzyaaaj1)
endif
if(am_master)then
xyzzyaaap8=0
xyzzyaacg8=0.d0
xyzzyaabi8=0
do xyzzyaaao8=1,nspin
if(nele(xyzzyaaao8)==0)cycle
if(pmass(xyzzyaaao8)/=1.d0.or.pcharge(xyzzyaaao8)/=-1.d0)cycle
xyzzyaabh8=which_ssingle(xyzzyaaao8,spin_dep_full)
if(xyzzyaabh8==xyzzyaabi8)cycle
xyzzyaabi8=xyzzyaabh8
xyzzyaabf8=which_ssingle(xyzzyaaao8,spin_dep_gs)
xyzzyaabg8=which_ssingle(xyzzyaaao8,spin_dep_in)
if(xyzzyaaaj1>1)then
if(electron_system)then
if(xyzzyaaao8==1)then
if(nuc_nele(1)>0)then
call wout('UP SPIN')
call wout()
endif
else
if(nuc_nele(2)>0)then
call wout('DOWN SPIN')
call wout()
endif
endif
else
call wout('SPIN-TYPE: '//trim(i2s(xyzzyaabh8)) //'      (spin-dependen&
&ce '//trim(i2s(spin_dep_full))//')')
call wout()
endif
endif
do xyzzyaaaj8=1,xyzzyaaaf1
xyzzyaaai8=ion_prim(xyzzyaaaj8)
if(is_ae(xyzzyaaaj8))then
xyzzyaaak8=0
do xyzzyaaac8=1,num_k
xyzzyaaes8=xyzzyaaac8>num_real_k
xyzzyaaan8=0
if(xyzzyaaes8)xyzzyaaan8=1
do xyzzyaaat8=1,nband(xyzzyaaac8,xyzzyaabf8)
do xyzzyaaam8=0,xyzzyaaan8
xyzzyaaak8=xyzzyaaak8+1
call wout('Orbital '//trim(i2s(xyzzyaaak8))//' at position of ion '//t&
&rim(i2s(xyzzyaaaj8)))
if(abs(xyzzyaaef8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8))<xyzzyaadq8)then
call wout('Orbital s component effectively zero at this nucleus.')
call wout()
cycle
endif
if(xyzzyaabr8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)>0)then
do xyzzyaaaa8=1,xyzzyaabr8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)
tmpr=r2s(xyzzyaaep8(xyzzyaaaa8,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f12&
&.5)')
call wout('Found node at r = '//trim(tmpr)//'.')
enddo
endif
if(nint(disign(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8))==1)then
call wout('Sign of orbital at nucleus                : positive')
else
call wout('Sign of orbital at nucleus                : negative')
endif
tmpr=r2s(sqrt(rcusp_sq(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)),'(f7.4)')
call wout('Cusp radius (au)                          : '//trim(tmpr))
tmpr=r2s(xyzzyaaen8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f21.12)')
call wout('Value of uncorrected orbital at nucleus   : '//trim(tmpr))
tmpr=r2s(xyzzyaaef8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f21.12)')
call wout('Value of s part of orbital at nucleus     : '//trim(tmpr))
tmpr=r2s(xyzzyaaeh8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f21.12)')
call wout('Optimum corrected s orbital at nucleus    : '//trim(tmpr))
xyzzyaacg8=xyzzyaacg8+xyzzyaaei8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)
xyzzyaaap8=xyzzyaaap8+1
tmpr=r2s(xyzzyaaei8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f16.8)')
if(xyzzyaabr8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)==1)then
tmpr2=r2s(xyzzyaadp8,'(f4.2)')
call wout('Maximum deviation from ideal local energy : '//trim(tmpr)//&
&' (EXCL. NODE +- '//trim(tmpr2)//')')
elseif(xyzzyaabr8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)>1)then
tmpr2=r2s(xyzzyaadp8,'(f4.2)')
call wout('Maximum deviation from ideal local energy : '//trim(tmpr)//&
&' (EXCL. NODES +- '//trim(tmpr2)//')')
else
call wout('Maximum deviation from ideal local energy : '//trim(tmpr))
endif
tmpr=r2s(xyzzyaaeg8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f16.8)')
call wout('Effective nuclear charge                  : '//trim(tmpr))
if(xyzzyaabr8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)==1)then
tmpr=r2s(pshift(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f16.8)')
call wout('Orbital shift to clear node               : '//trim(tmpr))
elseif(xyzzyaabr8(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)>1)then
tmpr=r2s(pshift(xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f16.8)')
call wout('Orbital shift to clear nodes              : '//trim(tmpr))
endif
call wout()
if(xyzzyaaet8)then
tmpr=r2s(acusp(1,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(f16.8)')
call wout('Polynomial parameters                    : '//trim(tmpr))
do xyzzyaaaa8=2,5
tmpr=r2s(acusp(xyzzyaaaa8,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8),'(e16.8)')
call wout('                                           '//trim(tmpr))
enddo
endif
enddo
enddo
enddo
endif
enddo
enddo
if(xyzzyaaap8>0)then
tmpr=r2s(xyzzyaacg8/dble(xyzzyaaap8),'(f16.8)')
call wout('Maximum deviation from ideal (averaged over orbitals) : '//&
&trim(tmpr))
endif
call wout()
endif
deallocate(xyzzyaabr8,xyzzyaaef8,xyzzyaaeh8,xyzzyaaei8,xyzzyaaep8,xyzz&
&yaaej8,xyzzyaaeq8,xyzzyaabs8,xyzzyaaeg8,xyzzyaaen8)
endif
open_unit(xyzzyaaah8)=.false.
if(xyzzyaaeu8)deallocate(xyzzyaadw8,xyzzyaaea8,xyzzyaaeb8,xyzzyaadx8,x&
&yzzyaady8,xyzzyaadz8,xyzzyaadt8,xyzzyaadu8,xyzzyaadv8)
deallocate(xyzzyaabn8,xyzzyaadr8,xyzzyaafc8,xyzzyaads8,xyzzyaaed8,xyzz&
&yaaeo8,xyzzyaaee8,xyzzyaabp8,xyzzyaabq8)
if(isperiodic)then
if(.not.s_plot)deallocate(ao_m,alap_m,agra1_m,agra2_m,agra3_m)
if(use_backflow)then
deallocate(asderiv1_m,asderiv2_m,asderiv3_m,asderiv4_m,asderiv5_m,asde&
&riv6_m,stat=xyzzyaabe8)
endif
endif
if(excite)deallocate(xyzzyaaao1,xyzzyaaan1)
call timer('CUSPS',.false.)
contains
subroutine xyzzyaafd8(phi_0_min,emax_min)
implicit none
real(dp),intent(inout) :: phi_0_min
real(dp),intent(out) :: emax_min
integer xyzzyaaaa9
real(dp) xyzzyaaab9,xyzzyaaac9,xyzzyaaad9,xyzzyaaae9,xyzzyaaaf9,xyzzya&
&aag9,xyzzyaaah9,xyzzyaaai9,xyzzyaaaj9,xyzzyaaak9,xyzzyaaal9,xyzzyaaam&
&9,xyzzyaaan9,xyzzyaaao9,xyzzyaaap9,xyzzyaaaq9,xyzzyaaar9,xyzzyaaas9
real(dp),parameter :: xyzzyaaat9=1.61803399d0,xyzzyaaau9=100.d0,xyzzya&
&aav9=1.d-20
xyzzyaaab9=phi_0_min
xyzzyaaac9=xyzzyaaab9*1.1d0
call xyzzyaaff8(xyzzyaaab9,xyzzyaaae9)
call xyzzyaaff8(xyzzyaaac9,xyzzyaaaf9)
if(xyzzyaaaf9>xyzzyaaae9)then
xyzzyaaah9=xyzzyaaab9
xyzzyaaab9=xyzzyaaac9
xyzzyaaac9=xyzzyaaah9
xyzzyaaah9=xyzzyaaae9
xyzzyaaae9=xyzzyaaaf9
xyzzyaaaf9=xyzzyaaah9
endif
xyzzyaaad9=xyzzyaaac9+xyzzyaaat9*(xyzzyaaac9-xyzzyaaab9)
call xyzzyaaff8(xyzzyaaad9,xyzzyaaag9)
xyzzyaaaa9=0
bracket_minimum: do
if(xyzzyaaaf9<=xyzzyaaag9)exit bracket_minimum
xyzzyaaaa9=xyzzyaaaa9+1
if(xyzzyaaaa9>100)then
call wout('WARNING: Too many iterations to bracket minimum in cusp cor&
&rection setup.')
call wout()
call wout('  Orbital : '//trim(i2s(xyzzyaaak8)))
call wout('  Ion     : '//trim(i2s(xyzzyaaai8)))
call wout('  Spin    : '//trim(i2s(xyzzyaaao8)))
call wout()
phi_0_min=xyzzyaaad9
emax_min=xyzzyaaag9
return
endif
xyzzyaaaj9=(xyzzyaaac9-xyzzyaaab9)*(xyzzyaaaf9-xyzzyaaag9)
xyzzyaaai9=(xyzzyaaad9-xyzzyaaac9)*(xyzzyaaae9-xyzzyaaaf9)
if(abs(xyzzyaaai9)-abs(xyzzyaaaj9)>xyzzyaaav9)then
xyzzyaaak9=xyzzyaaac9-((xyzzyaaac9-xyzzyaaad9)*xyzzyaaai9-(xyzzyaaac9-&
&xyzzyaaab9)*xyzzyaaaj9)/(2.d0*(xyzzyaaai9-xyzzyaaaj9))
if((xyzzyaaac9-xyzzyaaak9)*(xyzzyaaak9-xyzzyaaad9)>xyzzyaaav9)then
call xyzzyaaff8(xyzzyaaak9,xyzzyaaam9)
if(xyzzyaaam9<=xyzzyaaag9)then
xyzzyaaab9=xyzzyaaac9
xyzzyaaac9=xyzzyaaak9
xyzzyaaae9=xyzzyaaaf9
xyzzyaaaf9=xyzzyaaam9
exit bracket_minimum
elseif(xyzzyaaam9>=xyzzyaaaf9)then
xyzzyaaad9=xyzzyaaak9
xyzzyaaag9=xyzzyaaam9
exit bracket_minimum
endif
xyzzyaaak9=xyzzyaaad9+xyzzyaaat9*(xyzzyaaad9-xyzzyaaac9)
else
xyzzyaaak9=xyzzyaaac9+xyzzyaaat9*(xyzzyaaac9-xyzzyaaak9)
xyzzyaaal9=xyzzyaaac9+xyzzyaaau9*(xyzzyaaad9-xyzzyaaac9)
if((xyzzyaaad9-xyzzyaaak9)*(xyzzyaaak9-xyzzyaaal9)>0.d0)xyzzyaaak9=xyz&
&zyaaal9
endif
else
xyzzyaaak9=xyzzyaaad9+xyzzyaaat9*(xyzzyaaad9-xyzzyaaac9)
endif
call xyzzyaaff8(xyzzyaaak9,xyzzyaaam9)
xyzzyaaab9=xyzzyaaac9
xyzzyaaac9=xyzzyaaad9
xyzzyaaad9=xyzzyaaak9
xyzzyaaae9=xyzzyaaaf9
xyzzyaaaf9=xyzzyaaag9
xyzzyaaag9=xyzzyaaam9
enddo bracket_minimum
xyzzyaaap9=xyzzyaaab9
xyzzyaaas9=xyzzyaaad9
if(abs(xyzzyaaad9-xyzzyaaac9)>abs(xyzzyaaac9-xyzzyaaab9))then
xyzzyaaaq9=xyzzyaaac9
xyzzyaaan9=xyzzyaaaf9
xyzzyaaar9=xyzzyaaac9+(2.d0-xyzzyaaat9)*(xyzzyaaad9-xyzzyaaac9)
call xyzzyaaff8(xyzzyaaar9,xyzzyaaao9)
else
xyzzyaaaq9=xyzzyaaac9-(2.d0-xyzzyaaat9)*(xyzzyaaac9-xyzzyaaab9)
call xyzzyaaff8(xyzzyaaaq9,xyzzyaaan9)
xyzzyaaar9=xyzzyaaac9
xyzzyaaao9=xyzzyaaaf9
endif
find_minimum: do
if(abs(xyzzyaaas9-xyzzyaaap9)<=xyzzyaadm8*(abs(xyzzyaaaq9)+abs(xyzzyaa&
&ar9)))exit find_minimum
if(xyzzyaaao9<xyzzyaaan9)then
xyzzyaaap9=xyzzyaaaq9
xyzzyaaaq9=xyzzyaaar9
xyzzyaaar9=(xyzzyaaat9-1.d0)*xyzzyaaar9+(2.d0-xyzzyaaat9)*xyzzyaaas9
xyzzyaaan9=xyzzyaaao9
call xyzzyaaff8(xyzzyaaar9,xyzzyaaao9)
else
xyzzyaaas9=xyzzyaaar9
xyzzyaaar9=xyzzyaaaq9
xyzzyaaaq9=(xyzzyaaat9-1.d0)*xyzzyaaaq9+(2.d0-xyzzyaaat9)*xyzzyaaap9
xyzzyaaao9=xyzzyaaan9
call xyzzyaaff8(xyzzyaaaq9,xyzzyaaan9)
endif
enddo find_minimum
if(xyzzyaaan9<xyzzyaaao9)then
phi_0_min=xyzzyaaaq9
emax_min=xyzzyaaan9
else
phi_0_min=xyzzyaaar9
emax_min=xyzzyaaao9
endif
end subroutine xyzzyaafd8
subroutine xyzzyaafe8(r1,r2,xyzzyaadf8)
implicit none
real(dp),intent(in) :: r1,r2
real(dp),intent(out) :: xyzzyaadf8
integer xyzzyaaaa10
integer,parameter :: xyzzyaaab10=40
real(dp) xyzzyaaac10,xyzzyaaad10,xyzzyaaae10,xyzzyaaaf10
real(dp),parameter :: xyzzyaaag10=1.d-6
if(.not.xyzzyaaez8)then
call wfdet_s(r2,0,1.d0,xyzzyaaak8,xyzzyaaai8,xyzzyaaao8,xyzzyaabq8(xyz&
&zyaaak8,xyzzyaabf8),iabs(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)),xyzzyaaaj&
&8,.false.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
if(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)>0)then
xyzzyaaae10=xyzzyaadi8(1)
else
xyzzyaaae10=xyzzyaadi8(2)
endif
else
call wfdet_s(r2,0,1.d0,xyzzyaabd8,xyzzyaaai8,xyzzyaaao8,xyzzyaabd8,vir&
&tual_k(vrt2full(xyzzyaabd8,xyzzyaaao8),xyzzyaaao8),xyzzyaaaj8,.true.,&
&xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
xyzzyaaae10=xyzzyaadi8(1)
endif
if(.not.xyzzyaaez8)then
call wfdet_s(r1,0,1.d0,xyzzyaaak8,xyzzyaaai8,xyzzyaaao8,xyzzyaabq8(xyz&
&zyaaak8,xyzzyaabf8),iabs(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)),xyzzyaaaj&
&8,.false.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
if(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)>0)then
xyzzyaaad10=xyzzyaadi8(1)
else
xyzzyaaad10=xyzzyaadi8(2)
endif
else
call wfdet_s(r1,0,1.d0,xyzzyaabd8,xyzzyaaai8,xyzzyaaao8,xyzzyaabd8,vir&
&tual_k(vrt2full(xyzzyaabd8,xyzzyaaao8),xyzzyaaao8),xyzzyaaaj8,.true.,&
&xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
xyzzyaaad10=xyzzyaadi8(1)
endif
if(xyzzyaaad10*xyzzyaaae10>=0.d0)call errstop_master('FIND_NODE','Node&
& position not bracketed.')
if(xyzzyaaad10<0.d0)then
xyzzyaadf8=r1
xyzzyaaac10=r2-r1
else
xyzzyaadf8=r2
xyzzyaaac10=r1-r2
endif
do xyzzyaaaa10=1,xyzzyaaab10
xyzzyaaac10=xyzzyaaac10*0.5d0
xyzzyaaaf10=xyzzyaadf8+xyzzyaaac10
if(.not.xyzzyaaez8)then
call wfdet_s(xyzzyaaaf10,0,1.d0,xyzzyaaak8,xyzzyaaai8,xyzzyaaao8,xyzzy&
&aabq8(xyzzyaaak8,xyzzyaabf8),iabs(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)),&
&xyzzyaaaj8,.false.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
if(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)>0)then
xyzzyaaae10=xyzzyaadi8(1)
else
xyzzyaaae10=xyzzyaadi8(2)
endif
else
call wfdet_s(xyzzyaaaf10,0,1.d0,xyzzyaabd8,xyzzyaaai8,xyzzyaaao8,xyzzy&
&aabd8,virtual_k(vrt2full(xyzzyaabd8,xyzzyaaao8),xyzzyaaao8),xyzzyaaaj&
&8,.true.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
xyzzyaaae10=xyzzyaadi8(1)
endif
if(xyzzyaaae10<=0.d0)xyzzyaadf8=xyzzyaaaf10
if(abs(xyzzyaaac10)<xyzzyaaag10.or.xyzzyaaae10==0.d0)return
enddo
call errstop_master('FIND_NODE','Too many bisections - should not happ&
&en.')
end subroutine xyzzyaafe8
subroutine xyzzyaaff8(phi_0,emax)
implicit none
real(dp),intent(in) :: phi_0
real(dp),intent(out) :: emax
real(dp) xyzzyaaaa11
if(phi_0/=0.d0)then
xyzzyaabr1=xyzzyaacn8*(1.d0+xyzzyaact8/phi_0)
else
call errstop_master('ORB_SOLVE','Division by phi_0=0 - should not happ&
&en.')
endif
if(xyzzyaaco8==0.d0)then
xyzzyaacq8=1.d0
xyzzyaacs8=1.d0/xyzzyaack8
else
if(phi_0-xyzzyaaco8==0.d0)call errstop_master('ORB_SOLVE','Division by&
& zero in cusp setup - should not happen.')
xyzzyaacq8=phi_0/(phi_0-xyzzyaaco8)
xyzzyaacs8=1.d0/(xyzzyaack8-xyzzyaaco8)
endif
call xyzzyaafj8(phi_0)
call xyzzyaafk8
call xyzzyaafh8(xyzzyaacd8,xyzzyaaby8,xyzzyaabr1)
emax=0.d0
if(xyzzyaacn8/=1.d0)then
do xyzzyaaaa8=1,xyzzyaaaz8
xyzzyaabt8=xyzzyaadr8(xyzzyaaaa8)
xyzzyaadb8=xyzzyaabt8
xyzzyaadc8=xyzzyaacd8
xyzzyaacb8=0.d0
do xyzzyaaab8=1,7
xyzzyaadb8=xyzzyaadb8*xyzzyaabt8
xyzzyaadc8=xyzzyaadc8*xyzzyaacd8
xyzzyaacb8=xyzzyaacb8+xyzzyaada8(xyzzyaaab8)*(xyzzyaadb8-xyzzyaadc8)
enddo
xyzzyaacb8=xyzzyaacb8*xyzzyaabr1**2+xyzzyaaby8
call xyzzyaafi8(xyzzyaabt8,xyzzyaacc8)
if(xyzzyaafc8(xyzzyaaaa8))then
xyzzyaaaa11=abs(xyzzyaacc8-xyzzyaacb8)
emax=max(emax,xyzzyaaaa11)
endif
enddo
else
do xyzzyaaaa8=1,xyzzyaaaz8
xyzzyaabt8=xyzzyaadr8(xyzzyaaaa8)
call xyzzyaafi8(xyzzyaabt8,xyzzyaacc8)
if(xyzzyaafc8(xyzzyaaaa8))then
xyzzyaaaa11=abs(xyzzyaacc8-xyzzyaaby8)
emax=max(emax,xyzzyaaaa11)
endif
enddo
endif
if(emax==0.d0)call errstop_master('ORB_SOLVE','Problem in cusp correct&
&ion setup - emax is zero.')
end subroutine xyzzyaaff8
subroutine xyzzyaafg8(rrr,el,zzz)
implicit none
real(dp),intent(in) :: rrr,zzz
real(dp),intent(out) :: el
real(dp) xyzzyaaaa12,xyzzyaaab12
if(.not.xyzzyaaez8)then
call wfdet_s(rrr,1,1.d0,xyzzyaaak8,xyzzyaaai8,xyzzyaaao8,xyzzyaabq8(xy&
&zzyaaak8,xyzzyaabf8),iabs(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)),xyzzyaaa&
&j8,.false.,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
if(xyzzyaabp8(xyzzyaaak8,xyzzyaabf8)>0)then
xyzzyaack8=xyzzyaadi8(1)
xyzzyaacl8=xyzzyaadj8(1)
xyzzyaacm8=xyzzyaadk8(1)
else
xyzzyaack8=xyzzyaadi8(2)
xyzzyaacl8=xyzzyaadj8(2)
xyzzyaacm8=xyzzyaadk8(2)
endif
else
call wfdet_s(rrr,1,1.d0,xyzzyaabd8,xyzzyaaai8,xyzzyaaao8,xyzzyaabd8,vi&
&rtual_k(vrt2full(xyzzyaabd8,xyzzyaaao8),xyzzyaaao8),xyzzyaaaj8,.true.&
&,xyzzyaabg8,xyzzyaadi8,xyzzyaadj8,xyzzyaadk8)
xyzzyaack8=xyzzyaadi8(1)
xyzzyaacl8=xyzzyaadj8(1)
xyzzyaacm8=xyzzyaadk8(1)
endif
if(xyzzyaack8==0.d0)then
xyzzyaaaa12=0.d0
call errstop_master('GAUSS_ORB','Orbital divide by zero error in cusp &
&setup.')
else
xyzzyaaaa12=1.d0/xyzzyaack8
endif
xyzzyaaab12=1.d0/rrr
el=(-xyzzyaacl8*xyzzyaaab12-0.5d0*xyzzyaacm8)*xyzzyaaaa12-zzz*xyzzyaaa&
&b12
end subroutine xyzzyaafg8
subroutine xyzzyaafh8(rrr,el,zzz)
implicit none
real(dp),intent(in) :: rrr,zzz
real(dp),intent(out) :: el
real(dp) xyzzyaaaa13,xyzzyaaab13
if(xyzzyaack8==0.d0)then
xyzzyaaaa13=0.d0
call errstop_master('GAUSS_ORB_ENERGY','Orbital divide by zero error i&
&n cusp setup.')
else
xyzzyaaaa13=1.d0/xyzzyaack8
endif
xyzzyaaab13=1.d0/rrr
el=(-xyzzyaacl8*xyzzyaaab13-0.5d0*xyzzyaacm8)*xyzzyaaaa13-zzz*xyzzyaaa&
&b13
end subroutine xyzzyaafh8
subroutine xyzzyaafi8(r,el)
use slaarnabt,only : exp_protect
implicit none
real(dp),intent(in) :: r
real(dp),intent(out) :: el
real(dp) xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14,xyzzyaaad14,xyzzyaaae14,x&
&yzzyaaaf14,xyzzyaaag14
xyzzyaaaa14=xyzzyaaer8(1,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)
xyzzyaaab14=xyzzyaaer8(2,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)
xyzzyaaac14=xyzzyaaer8(3,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)
xyzzyaaad14=xyzzyaaer8(4,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)
xyzzyaaae14=xyzzyaaer8(5,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)
xyzzyaaaf14=xyzzyaaae14*r
xyzzyaabv8=xyzzyaaaa14+(xyzzyaaab14+(xyzzyaaac14+(xyzzyaaad14+xyzzyaaa&
&f14)*r)*r)*r
xyzzyaabw8=xyzzyaaab14+(2.d0*xyzzyaaac14+(3.d0*xyzzyaaad14+4.d0*xyzzya&
&aaf14)*r)*r
xyzzyaabx8=2.d0*xyzzyaaac14+(6.d0*xyzzyaaad14+12.d0*xyzzyaaaf14)*r
xyzzyaacr8=xyzzyaacp8*exp_protect(xyzzyaabv8)
if(xyzzyaacr8==0.d0.and.xyzzyaaco8==0.d0)then
xyzzyaacr8=1.d0
else
if(xyzzyaaco8+xyzzyaacr8==0.d0)then
call errstop_master('POLY_ORB','About to divide by zero in cusp setup &
&- should not happen')
endif
xyzzyaacr8=xyzzyaacr8/(xyzzyaaco8+xyzzyaacr8)
endif
xyzzyaaag14=1.d0/r
el=-(xyzzyaabw8*xyzzyaaag14+0.5d0*(xyzzyaabw8**2+xyzzyaabx8))*xyzzyaac&
&r8-xyzzyaabr1*xyzzyaaag14
end subroutine xyzzyaafi8
subroutine xyzzyaafj8(phi_0)
implicit none
real(dp),intent(in) :: phi_0
xyzzyaacf8(1)=log(abs(xyzzyaack8-xyzzyaaco8))
xyzzyaacf8(2)=xyzzyaacl8*xyzzyaacs8
xyzzyaacf8(3)=xyzzyaacm8*xyzzyaacs8
xyzzyaacf8(4)=-xyzzyaabr1*xyzzyaacq8
xyzzyaacf8(5)=log(abs(phi_0-xyzzyaaco8))
end subroutine xyzzyaafj8
subroutine xyzzyaafk8
implicit none
real(dp) xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16,xyzzyaaae16
xyzzyaaaa16=1.d0/xyzzyaacd8
xyzzyaaab16=xyzzyaaaa16*xyzzyaaaa16
xyzzyaaac16=xyzzyaaab16*(xyzzyaacf8(1)-xyzzyaacf8(5))
xyzzyaaad16=xyzzyaacf8(2)+xyzzyaacf8(4)
xyzzyaaae16=0.5d0*(xyzzyaacf8(3)-xyzzyaacf8(2)**2)
xyzzyaaer8(1,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaacf8(5)
xyzzyaaer8(2,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=xyzzyaacf8(4)
xyzzyaaer8(3,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=6.d0*xyzzyaaac16-3.d0*x&
&yzzyaaaa16*xyzzyaaad16+xyzzyaaae16
xyzzyaaer8(4,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=(-8.d0*xyzzyaaac16+xyzz&
&yaaaa16*(5.d0*xyzzyaacf8(2)+3.d0*xyzzyaacf8(4))-2.d0*xyzzyaaae16)*xyz&
&zyaaaa16
xyzzyaaer8(5,xyzzyaaak8,xyzzyaaaj8,xyzzyaabh8)=(3.d0*xyzzyaaac16-xyzzy&
&aaaa16*(xyzzyaacf8(2)+xyzzyaaad16)+xyzzyaaae16)*xyzzyaaab16
end subroutine xyzzyaafk8
subroutine xyzzyaafl8
implicit none
integer xyzzyaaaa17
real(dp) xyzzyaaab17
call xyzzyaafi8(xyzzyaacd8,xyzzyaaab17)
call xyzzyaafg8(xyzzyaacd8,xyzzyaaby8,xyzzyaabr1)
if(abs(xyzzyaabv8-xyzzyaacf8(1))>1.d-6.or.abs(xyzzyaabw8-xyzzyaacf8(2)&
&)>1.d-6.or.abs(xyzzyaabx8-xyzzyaacf8(3)+xyzzyaacf8(2)**2)>1.d-6.or.ab&
&s(xyzzyaaab17-xyzzyaaby8)>1.d-6 )then
call wout()
call wout('Orbital: '//trim(i2s(xyzzyaaak8))//'  Ion: '//trim(i2s(xyzz&
&yaaai8))//'  Spin: '//trim(i2s(xyzzyaaao8)))
call wout('The following pairs of numbers should be equal:')
tmpr=r2s(xyzzyaabv8,'(f16.8)')
tmpr2=r2s(xyzzyaacf8(1),'(f16.8)')
call wout('poly0, x(1)                              : '//trim(tmpr)//'&
& '//trim(tmpr2))
tmpr=r2s(xyzzyaabw8,'(f16.8)')
tmpr2=r2s(xyzzyaacf8(2),'(f16.8)')
call wout('poly1, x(2)                              : '//trim(tmpr)//'&
& '//trim(tmpr2))
tmpr=r2s(xyzzyaabx8,'(f16.8)')
tmpr2=r2s((xyzzyaacf8(3)-xyzzyaacf8(2)**2),'(f16.8)')
call wout('poly2, x(3)-x(2)^2                       : '//trim(tmpr)//'&
& '//trim(tmpr2))
tmpr=r2s(xyzzyaaab17,'(f16.8)')
tmpr2=r2s(xyzzyaaby8,'(f16.8)')
call wout('el_rc_p, el_rc                           : '//trim(tmpr)//'&
& '//trim(tmpr2))
call errwarn('CUSP_TEST','Error in applying constraints.')
endif
do xyzzyaaaa17=1,xyzzyaaaz8
xyzzyaabt8=xyzzyaadr8(xyzzyaaaa17)
call xyzzyaafi8(xyzzyaabt8,xyzzyaacc8)
xyzzyaadw8(xyzzyaaaa17)=xyzzyaacp8*exp(xyzzyaabv8)+xyzzyaaco8
xyzzyaaea8(xyzzyaaaa17)=xyzzyaacp8*exp(xyzzyaabv8)*xyzzyaabw8
xyzzyaaeb8(xyzzyaaaa17)=xyzzyaacp8*exp(xyzzyaabv8)*(xyzzyaabx8+xyzzyaa&
&bw8**2+(2.d0*xyzzyaabw8)/xyzzyaabt8)
xyzzyaadt8(xyzzyaaaa17)=xyzzyaacc8
call xyzzyaafg8(xyzzyaabt8,xyzzyaacc8,xyzzyaabr1)
xyzzyaadx8(xyzzyaaaa17)=xyzzyaack8
xyzzyaady8(xyzzyaaaa17)=xyzzyaacl8
xyzzyaadz8(xyzzyaaaa17)=xyzzyaacm8
xyzzyaadu8(xyzzyaaaa17)=xyzzyaacc8
if(xyzzyaacn8/=1.d0)then
xyzzyaadb8=xyzzyaabt8
xyzzyaadc8=xyzzyaacd8
xyzzyaadv8(xyzzyaaaa17)=0.d0
do xyzzyaaab8=1,7
xyzzyaadb8=xyzzyaadb8*xyzzyaabt8
xyzzyaadc8=xyzzyaadc8*xyzzyaacd8
xyzzyaadv8(xyzzyaaaa17)=xyzzyaadv8(xyzzyaaaa17)+xyzzyaada8(xyzzyaaab8)&
&*(xyzzyaadb8-xyzzyaadc8)
enddo
xyzzyaadv8(xyzzyaaaa17)=xyzzyaadv8(xyzzyaaaa17)*xyzzyaabr1**2+xyzzyaab&
&y8
else
xyzzyaadv8(xyzzyaaaa17)=xyzzyaaby8
endif
enddo
end subroutine xyzzyaafl8
end subroutine cusp_setup
subroutine read_gwfmolorb
implicit none
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xy&
&zzyaaaf18,xyzzyaaag18,xyzzyaaah18,xyzzyaaai18,xyzzyaaaj18,xyzzyaaak18&
&,xyzzyaaal18
integer :: xyzzyaaam18
character(80) char_80
real(dp),parameter :: xyzzyaaan18=1.d-4
logical xyzzyaaao18
xyzzyaabm1=.false.
xyzzyaabn1=.false.
xyzzyaabo1=.false.
call open_units(xyzzyaaaa1,xyzzyaaaa18)
if(xyzzyaaaa18/=0)call errstop('READ_GWFMOLORB','Cannot find free IO u&
&nit.')
open(unit=xyzzyaaaa1,file='correlation.data',status='old',iostat=xyzzy&
&aaaa18)
if(xyzzyaaaa18/=0)call errstop('READ_GWFMOLORB','Problem opening corre&
&lation.data.')
if(am_master)then
call wout()
call wout('Modified molecular orbitals')
call wout('===========================')
call wout('Reading molecular orbital information from correlation.data&
& file.')
call wout()
endif
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa18)char_80
if(xyzzyaaaa18>0)call errstop_master('READ_GWFMOLORB','Problem reading&
& correlation.data. Please check this file.')
if(xyzzyaaaa18<0)call errstop_master('READ_GWFMOLORB','Could not find &
&"START MOLORBMODS" in correlation.data.')
if(trim(adjustl(char_80))=='START MOLORBMODS')exit
enddo
read(xyzzyaaaa1,*,err=10,end=10)
read(xyzzyaaaa1,'(a)',err=10,end=10)param_title
if(am_master)call wout('Title: '//trim(adjustl(param_title)))
mainloop: do
read(xyzzyaaaa1,'(a)',err=10,end=10)char_80
if(trim(adjustl(char_80))=='START GAUSSIAN MO COEFFICIENTS')then
if(xyzzyaabm1)call errstop_master('READ_GWFMOLORB','Only one GAUSSIAN &
&MO COEFFICIENTS term allowed.')
if(am_master)then
call wout()
call wout()
call wout(' Gaussian molecular orbital coefficients:')
call wout()
endif
xyzzyaabm1=.true.
xyzzyaaah18=maxval(nband(1,:))*num_ao*xyzzyaaai1
allocate(xyzzyaabk1(xyzzyaaah18),stat=xyzzyaaam18)
call check_alloc(xyzzyaaam18,'READ_GWFMOLORB','rcktemp')
xyzzyaaad18=1
xyzzyaaal18=0
do xyzzyaaai18=1,nspin
xyzzyaaaj18=which_ssingle(xyzzyaaai18,spin_dep_gs)
if(xyzzyaaaj18==xyzzyaaal18)cycle
xyzzyaaal18=xyzzyaaaj18
xyzzyaaak18=which_ssingle(xyzzyaaai18,spin_dep_in)
do xyzzyaaae18=1,nband(1,xyzzyaaaj18)
xyzzyaabk1(xyzzyaaad18:xyzzyaaad18+num_ao-1)=rck(xyzzyaaae18,1:num_ao,&
&1,xyzzyaaak18)
xyzzyaaad18=xyzzyaaad18+num_ao
enddo
enddo
xyzzyaabh1=xyzzyaaah18
allocate(xyzzyaabb1(xyzzyaabh1),xyzzyaabe1(xyzzyaabh1),xyzzyaabl1(xyzz&
&yaaah18),stat=xyzzyaaam18)
call check_alloc(xyzzyaaam18,'READ_GWFMOLORB','indexrcktemp')
xyzzyaabb1=0.d0
xyzzyaabl1=0
xyzzyaabb1=xyzzyaabk1
xyzzyaabe1(:)=1
read(xyzzyaaaa1,*,err=10,end=10)
read_params_mocoeff: do xyzzyaaab18=1,xyzzyaabh1
read(xyzzyaaaa1,*,iostat=xyzzyaaaa18)xyzzyaabb1(xyzzyaaab18),xyzzyaabe&
&1(xyzzyaaab18)
if(xyzzyaaaa18/=0.and.xyzzyaaab18==1)then
backspace xyzzyaaaa1
exit read_params_mocoeff
elseif(xyzzyaaaa18/=0)then
if(am_master)then
call wout()
call wout(' A total of '//trim(i2s(xyzzyaabh1))//' coefficients is exp&
&ected.')
endif
call errstop_master('READ_GWFMOLORB','Not enough MO coefficients are s&
&upplied.')
endif
if(am_master)then
call display_param(xyzzyaabb1(xyzzyaaab18),xyzzyaabe1(xyzzyaaab18),'mo&
&c_'//trim(i2s(xyzzyaaab18)))
endif
if(xyzzyaabe1(xyzzyaaab18)/=0.and.xyzzyaabe1(xyzzyaaab18)/=1)call errs&
&top_master('READ_GWFMOLORB','Optimizable flag should be 0 or 1.')
enddo read_params_mocoeff
xyzzyaabk1=xyzzyaabb1
if(am_master)then
if(xyzzyaaaa18/=0)then
call wout(' No optimized coefficients supplied. Use the ones in gwfn.d&
&ata.')
call wout()
do xyzzyaaab18=1,xyzzyaabh1
call display_param(xyzzyaabb1(xyzzyaaab18),xyzzyaabe1(xyzzyaaab18),'mo&
&c_'//trim(i2s(xyzzyaaab18)))
enddo
endif
xyzzyaaag18=count(xyzzyaabe1>0)
call wout('  No. of free params    :  '//trim(i2s(xyzzyaaag18)))
call wout()
endif
xyzzyaaaf18=0
do xyzzyaaab18=1,xyzzyaabh1
if(xyzzyaabe1(xyzzyaaab18)==0)xyzzyaaaf18=xyzzyaaaf18+1
enddo
if(xyzzyaaaf18==0)then
if(am_master)then
call wout()
call wout(' No specification of MO coefficient optimization, so will u&
&se defaults:')
call wout('  1) Coefficients smaller than an ad hoc (hard coded) toler&
&ance set as fixed.')
call wout('  2) Largest coefficient of each orbital set as fixed.')
call wout('  3) Of remaining, smallest Norb-1 coefficients in each orb&
&ital set as fixed.')
endif
do xyzzyaaab18=1,xyzzyaabh1
if(abs(xyzzyaabb1(xyzzyaaab18))<xyzzyaaan18)then
xyzzyaabe1(xyzzyaaab18)=0
xyzzyaaaf18=xyzzyaaaf18+1
endif
enddo
xyzzyaaac18=1
xyzzyaaal18=0
do xyzzyaaai18=1,nspin
xyzzyaaaj18=which_ssingle(xyzzyaaai18,spin_dep_gs)
if(xyzzyaaaj18==xyzzyaaal18)cycle
xyzzyaaal18=xyzzyaaaj18
xyzzyaaak18=which_ssingle(xyzzyaaai18,spin_dep_in)
do xyzzyaaae18=1,nband(1,xyzzyaaaj18)
do xyzzyaaad18=1,nband(1,xyzzyaaaj18)
if(xyzzyaaad18==1)then
xyzzyaaab18=xyzzyaaac18-1+sum(maxloc(abs(xyzzyaabb1(xyzzyaaac18:xyzzya&
&aac18+num_ao-1)),xyzzyaabe1(xyzzyaaac18:xyzzyaaac18+num_ao-1)/=0))
else
xyzzyaaab18=xyzzyaaac18-1+sum(minloc(abs(xyzzyaabb1(xyzzyaaac18:xyzzya&
&aac18+num_ao-1)),xyzzyaabe1(xyzzyaaac18:xyzzyaaac18+num_ao-1)/=0))
endif
xyzzyaabe1(xyzzyaaab18)=0
xyzzyaaaf18=xyzzyaaaf18+1
enddo
xyzzyaaac18=xyzzyaaac18+num_ao
enddo
enddo
if(am_master)then
call wout(' Relaxing '//trim(i2s(xyzzyaabh1-xyzzyaaaf18))//' of '//tri&
&m(i2s(xyzzyaabh1))//'.')
call wout()
endif
endif
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa18)char_80
if(xyzzyaaaa18/=0)call errstop_master('READ_GWFMOLORB','String "END GA&
&USSIAN MO COEFFICIENTS" not found.')
if(trim(adjustl(char_80))/='END GAUSSIAN MO COEFFICIENTS')call errstop&
&_master('READ_GWFMOLORB','String "END GAUSSIAN MO COEFFICIENTS" not f&
&ound.')
elseif(trim(adjustl(char_80))=='START GAUSSIAN EXPONENTS')then
if(xyzzyaabn1)call errstop_master('READ_GWFMOLORB','Only one GAUSSIAN &
&EXPONENTS term allowed.')
if(am_master)then
call wout()
call wout(' Exponents of Gaussian primitives:')
call wout()
endif
xyzzyaabn1=.true.
xyzzyaabi1=xyzzyaaac1
allocate(xyzzyaabc1(xyzzyaabi1),xyzzyaabf1(xyzzyaabi1),stat=xyzzyaaam1&
&8)
call check_alloc(xyzzyaaam18,'READ_GWFMOLORB','opt_exponent')
xyzzyaabf1(:)=1
read(xyzzyaaaa1,*,err=10,end=10)
read_params_exponent: do xyzzyaaab18=1,xyzzyaaac1
read(xyzzyaaaa1,*,iostat=xyzzyaaaa18)xyzzyaabc1(xyzzyaaab18),xyzzyaabf&
&1(xyzzyaaab18)
if(xyzzyaaaa18/=0.and.xyzzyaaab18==1)then
backspace xyzzyaaaa1
xyzzyaabc1(:)=exponent(:)
exit read_params_exponent
elseif(xyzzyaaaa18/=0)then
if(am_master)then
call wout()
call wout(' A total of '//trim(i2s(xyzzyaabi1))//' exponents are expec&
&ted.')
endif
call errstop_master('READ_GWFMOLORB','Not enough Gaussian exponents ar&
&e supplied.')
endif
if(xyzzyaabf1(xyzzyaaab18)/=0.and.xyzzyaabf1(xyzzyaaab18)/=1)call errs&
&top_master('READ_GWFMOLORB','Optimizable flag should be 0 or 1.')
if(am_master)then
call display_param(xyzzyaabc1(xyzzyaaab18),xyzzyaabf1(xyzzyaaab18),'ex&
&po_'//trim(i2s(xyzzyaaab18)))
endif
enddo read_params_exponent
if(am_master)then
if(xyzzyaaaa18/=0)then
call wout(' No Gaussian exponents supplied. Use the ones in gwfn.data.&
&')
call wout()
do xyzzyaaab18=1,xyzzyaabi1
call display_param(xyzzyaabc1(xyzzyaaab18),xyzzyaabf1(xyzzyaaab18),'ex&
&po_'//trim(i2s(xyzzyaaab18)))
enddo
endif
xyzzyaaag18=count(xyzzyaabf1>0)
call wout('  No. of free params     :  '//trim(i2s(xyzzyaaag18)))
call wout()
endif
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa18)char_80
if(xyzzyaaaa18/=0)call errstop_master('READ_GWFMOLORB','String "END GA&
&USSIAN EXPONENTS" not found.')
if(trim(adjustl(char_80))/='END GAUSSIAN EXPONENTS')call errstop_maste&
&r('READ_GWFMOLORB','String "END GAUSSIAN EXPONENTS" not found.')
elseif(trim(adjustl(char_80))=='START CONTRACTION COEFFICIENTS')then
if(xyzzyaabo1)call errstop_master('READ_GWFMOLORB','Only one CONTRACTI&
&ON COEFFICIENTS term allowed.')
if(am_master)then
call wout()
call wout(' Contraction coefficients:')
call wout()
endif
xyzzyaabo1=.true.
if(any(shell_am(:)==2))then
xyzzyaabj1=2*xyzzyaaac1
else
xyzzyaabj1=xyzzyaaac1
endif
allocate(xyzzyaabd1(xyzzyaabj1),xyzzyaabg1(xyzzyaabj1),stat=xyzzyaaam1&
&8)
call check_alloc(xyzzyaaam18,'READ_GWFMOLORB','opt_concoeff')
xyzzyaabg1(:)=1
read(xyzzyaaaa1,*,err=10,end=10)
read_params_concoeff: do xyzzyaaab18=1,xyzzyaabj1
read(xyzzyaaaa1,*,iostat=xyzzyaaaa18)xyzzyaabd1(xyzzyaaab18),xyzzyaabg&
&1(xyzzyaaab18)
if(xyzzyaaaa18/=0.and.xyzzyaaab18==1)then
backspace xyzzyaaaa1
xyzzyaabd1(1:xyzzyaaac1)=c_prim(1:xyzzyaaac1)
if(any(shell_am(:)==2))then
xyzzyaabd1(xyzzyaaac1+1:2*xyzzyaaac1)=c_prim2(1:xyzzyaaac1)
do xyzzyaaac18=1,xyzzyaaac1
if(c_prim2(xyzzyaaac18)==0.d0)xyzzyaabg1(xyzzyaaac1+xyzzyaaac18)=0
enddo
endif
exit read_params_concoeff
elseif(xyzzyaaaa18/=0)then
if(am_master)then
call wout()
call wout(' A total of '//trim(i2s(xyzzyaabj1))//' contraction coeffic&
&ients are expected.')
endif
call errstop_master('READ_GWFMOLORB','Not enough Gaussian contraction &
&coefficients are supplied.')
endif
if(xyzzyaabg1(xyzzyaaab18)/=0.and.xyzzyaabg1(xyzzyaaab18)/=1)call errs&
&top_master('READ_GWFMOLORB','Optimizable flag should be 0 or 1.')
if(xyzzyaaab18>xyzzyaaac1)then
if(any(shell_am(:)==2))then
xyzzyaaao18=c_prim2(xyzzyaaab18-xyzzyaaac1)==0.d0.and.(xyzzyaabg1(xyzz&
&yaaab18)/=0.or.xyzzyaabd1(xyzzyaaab18)/=0)
else
xyzzyaaao18=xyzzyaabg1(xyzzyaaab18)/=0.or.xyzzyaabd1(xyzzyaaab18)/=0
endif
if(xyzzyaaao18)then
xyzzyaabg1(xyzzyaaab18)=0
xyzzyaabd1(xyzzyaaab18)=0.d0
call wout()
call errwarn('READ_GWFMOLORB','overriding values and optimizable flag &
&for conc_'//trim(i2s(xyzzyaaab18))//'. Non-sp shell contraction coeff&
&icients should have zero values and be non-optimizable.')
endif
endif
if(am_master)call display_param(xyzzyaabd1(xyzzyaaab18),xyzzyaabg1(xyz&
&zyaaab18),'conc_'//trim(i2s(xyzzyaaab18)))
enddo read_params_concoeff
if(am_master)then
if(xyzzyaaaa18/=0)then
call wout(' No Gaussian contraction coefficients supplied. Use the one&
&s in gwfn.data.')
call wout()
do xyzzyaaab18=1,xyzzyaabj1
call display_param(xyzzyaabd1(xyzzyaaab18),xyzzyaabg1(xyzzyaaab18),'co&
&nc_'//trim(i2s(xyzzyaaab18)))
enddo
endif
xyzzyaaag18=count(xyzzyaabg1>0)
call wout('  No. of free params     :  '//trim(i2s(xyzzyaaag18)))
call wout()
endif
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa18)char_80
if(xyzzyaaaa18/=0)call errstop_master('READ_GWFMOLORB','String "END CO&
&NTRACTION COEFFICIENTS" not found.')
if(trim(adjustl(char_80))/='END CONTRACTION COEFFICIENTS')call errstop&
&_master('READ_GWFMOLORB','String "END CONTRACTION COEFFICIENTS" not f&
&ound.')
elseif(trim(adjustl(char_80))=='END MOLORBMODS')then
exit mainloop
elseif(trim(adjustl(char_80))/="")then
call errstop_master('READ_GWFMOLORB','Expecting e.g. "START GAUSSIAN M&
&O COEFFICIENTS" in correlation.data.')
endif
enddo mainloop
do
read(xyzzyaaaa1,'(a)',iostat=xyzzyaaaa18)char_80
if(xyzzyaaaa18<0)exit
if(xyzzyaaaa18>0)call errstop_master('READ_GWFMOLORB','Problem reading&
& correlation.data. Please check this file.')
if(trim(adjustl(char_80))=='START MOLORBMODS')call errstop_master('REA&
&D_GWFMOLORB','There seems to be more than one set of  optimized molec&
&ular orbital coefficients in correlation.data.')
enddo
close(xyzzyaaaa1)
call xyzzyaadc1
return
10 call errstop_master('READ_GWFMOLORB','Problem reading optimized mol&
&ecular orbital coefficients in correlation.data. Please check this fi&
&le.')
end subroutine read_gwfmolorb
subroutine write_gwfmolorb(correlation_name)
implicit none
character(20),intent(in) :: correlation_name
integer xyzzyaaaa19,xyzzyaaab19
logical xyzzyaaac19
if(.not.am_master)return
inquire(file=trim(correlation_name),exist=xyzzyaaac19)
if(xyzzyaaac19)then
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='old',position&
&='append',iostat=xyzzyaaaa19)
else
open(unit=xyzzyaaaa1,file=trim(correlation_name),status='replace',iost&
&at=xyzzyaaaa19)
endif
if(xyzzyaaaa19/=0)call errstop('WRITE_GWFMOLORB','Problem opening '//t&
&rim(correlation_name)//'.')
write(xyzzyaaaa1,*)'START MOLORBMODS'
write(xyzzyaaaa1,*)'Title'
write(xyzzyaaaa1,*)trim(param_title)
if(xyzzyaabm1)then
write(xyzzyaaaa1,*)'START GAUSSIAN MO COEFFICIENTS'
write(xyzzyaaaa1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaab19=1,xyzzyaabh1
write(xyzzyaaaa1,*)xyzzyaabb1(xyzzyaaab19),xyzzyaabe1(xyzzyaaab19),'      !&
& moc_',trim(i2s(xyzzyaaab19))
enddo
write(xyzzyaaaa1,*)'END GAUSSIAN MO COEFFICIENTS'
endif
if(xyzzyaabn1)then
write(xyzzyaaaa1,*)'START GAUSSIAN EXPONENTS'
write(xyzzyaaaa1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaab19=1,xyzzyaabi1
write(xyzzyaaaa1,*)xyzzyaabc1(xyzzyaaab19),xyzzyaabf1(xyzzyaaab19),'      !&
& expo_',trim(i2s(xyzzyaaab19))
enddo
write(xyzzyaaaa1,*)'END GAUSSIAN EXPONENTS'
endif
if(xyzzyaabo1)then
write(xyzzyaaaa1,*)'START CONTRACTION COEFFICIENTS'
write(xyzzyaaaa1,*)'Parameter values  ;  Optimizable (0=NO; 1=YES)'
do xyzzyaaab19=1,xyzzyaabj1
write(xyzzyaaaa1,*)xyzzyaabd1(xyzzyaaab19),xyzzyaabg1(xyzzyaaab19),'      !&
& conc_',trim(i2s(xyzzyaaab19))
enddo
write(xyzzyaaaa1,*)'END CONTRACTION COEFFICIENTS'
endif
write(xyzzyaaaa1,*)'END MOLORBMODS'
close(xyzzyaaaa1)
end subroutine write_gwfmolorb
subroutine setup_gwfmolorb_params(nparam)
use slaarnabi,only : use_gpcc,gpcc_init_buffer
implicit none
integer,intent(out) :: nparam
integer xyzzyaaaa20
nparam=0
if(xyzzyaabm1)then
do xyzzyaaaa20=1,xyzzyaabh1
if(xyzzyaabe1(xyzzyaaaa20)==1)nparam=nparam+1
enddo
endif
if(xyzzyaabn1)then
do xyzzyaaaa20=1,xyzzyaabi1
if(xyzzyaabf1(xyzzyaaaa20)==1)nparam=nparam+1
enddo
endif
if(xyzzyaabo1)then
do xyzzyaaaa20=1,xyzzyaabj1
if(xyzzyaabg1(xyzzyaaaa20)==1)nparam=nparam+1
enddo
endif
xyzzyaaag1=nparam
if(cusp_correction)call xyzzyaacy1
if(use_gpcc)call gpcc_init_buffer(xyzzyaaag1)
end subroutine setup_gwfmolorb_params
subroutine finish_gwfmolorb_params
use slaarnabi,only : use_gpcc,gpcc_finalize_buffer
implicit none
if(cusp_correction)call xyzzyaacz1
if(use_gpcc)call gpcc_finalize_buffer
end subroutine finish_gwfmolorb_params
subroutine xyzzyaacy1
implicit none
integer xyzzyaaaa22
allocate(xyzzyaacd1(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1,0:xyzzyaaag1),xyz&
&zyaacg1(xyzzyaabq1,xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1,0:xyzzyaaag1),xyz&
&zyaace1(xyzzyaabp1,xyzzyaaaf1,xyzzyaaaj1,0:xyzzyaaag1),xyzzyaacf1(xyz&
&zyaabp1,xyzzyaaaf1,xyzzyaaaj1,0:xyzzyaaag1),stat=xyzzyaaaa22)
call check_alloc(xyzzyaaaa22,'SETUP_GWFMOLOB_PBUFFER','rcusp_sq_pbuffe&
&r etc.')
if(isperiodic)then
allocate(xyzzyaabs1(nemaxc,0:xyzzyaaag1),xyzzyaabt1(nemaxc,0:xyzzyaaag&
&1),xyzzyaabu1(nemaxc,0:xyzzyaaag1),xyzzyaabv1(nemaxc,0:xyzzyaaag1),xy&
&zzyaabw1(nemaxc,0:xyzzyaaag1),stat=xyzzyaaaa22)
call check_alloc(xyzzyaaaa22,'SETUP_GWFMOLOB_PBUFFER','exp_poly0 etc.'&
&)
if(use_backflow)then
allocate(xyzzyaabx1(nemaxc,0:xyzzyaaag1),xyzzyaaby1(nemaxc,0:xyzzyaaag&
&1),xyzzyaabz1(nemaxc,0:xyzzyaaag1),xyzzyaaca1(nemaxc,0:xyzzyaaag1),xy&
&zzyaacb1(nemaxc,0:xyzzyaaag1),xyzzyaacc1(nemaxc,0:xyzzyaaag1),stat=xy&
&zzyaaaa22)
call check_alloc(xyzzyaaaa22,'SETUP_GWFMOLOB_PBUFFER','exp_poly5 etc.'&
&)
endif
endif
end subroutine xyzzyaacy1
subroutine xyzzyaacz1
implicit none
deallocate(xyzzyaacd1,xyzzyaacg1,xyzzyaace1,xyzzyaacf1)
if(isperiodic)then
deallocate(xyzzyaabs1,xyzzyaabt1,xyzzyaabu1,xyzzyaabv1,xyzzyaabw1)
if(use_backflow)deallocate(xyzzyaabx1,xyzzyaaby1,xyzzyaabz1,xyzzyaaca1&
&,xyzzyaacb1,xyzzyaacc1)
endif
end subroutine xyzzyaacz1
subroutine get_gwfmolorb_params(params,has_lolim,lolim,has_hilim,hilim&
&,is_shallow,is_redundant,is_linear,is_loglinear,has_aderiv,affect_map&
&,label)
implicit none
character(2),intent(inout) :: label(xyzzyaaag1)
real(dp),intent(inout) :: params(xyzzyaaag1),lolim(xyzzyaaag1),hilim(x&
&yzzyaaag1)
logical,intent(inout) :: has_lolim(xyzzyaaag1),has_hilim(xyzzyaaag1),i&
&s_shallow(xyzzyaaag1),     is_redundant(xyzzyaaag1),is_linear(xyzzyaa&
&ag1),is_loglinear(xyzzyaaag1),has_aderiv(xyzzyaaag1),affect_map(xyzzy&
&aaag1,xyzzyaaag1)
integer xyzzyaaaa24,xyzzyaaab24
real(dp) xyzzyaaac24
has_lolim=.false.
lolim=0.d0
has_hilim=.false.
hilim=0.d0
is_shallow=.false.
is_redundant=.false.
is_linear=.false.
is_loglinear=.false.
has_aderiv=.false.
affect_map=.false.
do xyzzyaaaa24=1,xyzzyaaag1
affect_map(xyzzyaaaa24,xyzzyaaaa24)=.true.
enddo
label='OP'
xyzzyaaac24=1.1d-8
xyzzyaaaa24=0
if(xyzzyaabm1)then
do xyzzyaaab24=1,xyzzyaabh1
if(xyzzyaabe1(xyzzyaaab24)==1)then
xyzzyaaaa24=xyzzyaaaa24+1
params(xyzzyaaaa24)=xyzzyaabb1(xyzzyaaab24)
endif
enddo
endif
if(xyzzyaabn1)then
do xyzzyaaab24=1,xyzzyaabi1
if(xyzzyaabf1(xyzzyaaab24)==1)then
xyzzyaaaa24=xyzzyaaaa24+1
params(xyzzyaaaa24)=xyzzyaabc1(xyzzyaaab24)
has_lolim(xyzzyaaaa24)=.true.
lolim(xyzzyaaaa24)=xyzzyaaac24
endif
enddo
endif
if(xyzzyaabo1)then
do xyzzyaaab24=1,xyzzyaabj1
if(xyzzyaabg1(xyzzyaaab24)==1)then
xyzzyaaaa24=xyzzyaaaa24+1
params(xyzzyaaaa24)=xyzzyaabd1(xyzzyaaab24)
endif
enddo
endif
end subroutine get_gwfmolorb_params
subroutine put_gwfmolorb_params(params,ignore,iparam_buffer,prestore,b&
&ad_params)
use slaarnabi,only : use_gpcc,gpcc_save_to_buffer,gpcc_load_from_buffe&
&r
implicit none
integer,intent(in) :: iparam_buffer
real(dp),intent(inout) :: params(xyzzyaaag1)
logical,intent(in) :: ignore(xyzzyaaag1),prestore
logical,intent(out) :: bad_params
integer xyzzyaaaa25,xyzzyaaab25
bad_params=.false.
xyzzyaaab25=0
if(xyzzyaabm1)then
do xyzzyaaaa25=1,xyzzyaabh1
if(xyzzyaabe1(xyzzyaaaa25)==1)then
xyzzyaaab25=xyzzyaaab25+1
if(.not.ignore(xyzzyaaab25))xyzzyaabb1(xyzzyaaaa25)=params(xyzzyaaab25&
&)
endif
enddo
xyzzyaabk1=xyzzyaabb1
endif
if(xyzzyaabn1)then
do xyzzyaaaa25=1,xyzzyaabi1
if(xyzzyaabf1(xyzzyaaaa25)==1)then
xyzzyaaab25=xyzzyaaab25+1
if(.not.ignore(xyzzyaaab25))xyzzyaabc1(xyzzyaaaa25)=params(xyzzyaaab25&
&)
endif
enddo
endif
if(xyzzyaabo1)then
do xyzzyaaaa25=1,xyzzyaabj1
if(xyzzyaabg1(xyzzyaaaa25)==1)then
xyzzyaaab25=xyzzyaaab25+1
if(.not.ignore(xyzzyaaab25))xyzzyaabd1(xyzzyaaaa25)=params(xyzzyaaab25&
&)
endif
enddo
endif
call xyzzyaadc1
if(cusp_correction)then
if(prestore)then
call xyzzyaadb1(iparam_buffer)
return
endif
if(have_ae)call cusp_setup
call xyzzyaada1(iparam_buffer)
endif
if(use_gpcc)then
if(prestore)then
call gpcc_load_from_buffer(iparam_buffer)
else
call xyzzyaacn1(iparam_buffer,prestore)
call gpcc_save_to_buffer(iparam_buffer)
endif
endif
end subroutine put_gwfmolorb_params
subroutine xyzzyaada1(indx)
use slaarnaat, only: exp_poly0,exp_poly1,exp_poly2,exp_poly3,exp_poly4&
&,exp_poly5,exp_poly6,exp_poly7,exp_poly8,exp_poly9,exp_poly10,rcusp_s&
&q,acusp,disign,pshift
implicit none
integer,intent(in) :: indx
xyzzyaacd1(:,:,:,indx)=rcusp_sq(:,:,:)
xyzzyaacg1(:,:,:,:,indx)=acusp(:,:,:,:)
xyzzyaace1(:,:,:,indx)=disign(:,:,:)
xyzzyaacf1(:,:,:,indx)=pshift(:,:,:)
if(isperiodic)then
xyzzyaabs1(:,indx)=exp_poly0(:)
xyzzyaabt1(:,indx)=exp_poly1(:)
xyzzyaabu1(:,indx)=exp_poly2(:)
xyzzyaabv1(:,indx)=exp_poly3(:)
xyzzyaabw1(:,indx)=exp_poly4(:)
if(use_backflow)then
xyzzyaabx1(:,indx)=exp_poly5(:)
xyzzyaaby1(:,indx)=exp_poly6(:)
xyzzyaabz1(:,indx)=exp_poly7(:)
xyzzyaaca1(:,indx)=exp_poly8(:)
xyzzyaacb1(:,indx)=exp_poly9(:)
xyzzyaacc1(:,indx)=exp_poly10(:)
endif
endif
end subroutine xyzzyaada1
subroutine xyzzyaadb1(indx)
use slaarnaat, only: exp_poly0,exp_poly1,exp_poly2,exp_poly3,exp_poly4&
&,exp_poly5,exp_poly6,exp_poly7,exp_poly8,exp_poly9,exp_poly10,rcusp_s&
&q,acusp,disign,pshift
implicit none
integer,intent(in) :: indx
rcusp_sq(:,:,:)=xyzzyaacd1(:,:,:,indx)
acusp(:,:,:,:)=xyzzyaacg1(:,:,:,:,indx)
disign(:,:,:)=xyzzyaace1(:,:,:,indx)
pshift(:,:,:)=xyzzyaacf1(:,:,:,indx)
if(isperiodic)then
exp_poly0(:)=xyzzyaabs1(:,indx)
exp_poly1(:)=xyzzyaabt1(:,indx)
exp_poly2(:)=xyzzyaabu1(:,indx)
exp_poly3(:)=xyzzyaabv1(:,indx)
exp_poly4(:)=xyzzyaabw1(:,indx)
if(use_backflow)then
exp_poly5(:)=xyzzyaabx1(:,indx)
exp_poly6(:)=xyzzyaaby1(:,indx)
exp_poly7(:)=xyzzyaabz1(:,indx)
exp_poly8(:)=xyzzyaaca1(:,indx)
exp_poly9(:)=xyzzyaacb1(:,indx)
exp_poly10(:)=xyzzyaacc1(:,indx)
endif
endif
end subroutine xyzzyaadb1
subroutine xyzzyaadc1
implicit none
integer xyzzyaaaa28,xyzzyaaab28,xyzzyaaac28,xyzzyaaad28,xyzzyaaae28,xy&
&zzyaaaf28,xyzzyaaag28,xyzzyaaah28
if(xyzzyaabm1)then
xyzzyaaad28=1
xyzzyaaah28=0
do xyzzyaaae28=1,nspin
xyzzyaaaf28=which_ssingle(xyzzyaaae28,spin_dep_gs)
if(xyzzyaaaf28==xyzzyaaah28)cycle
xyzzyaaah28=xyzzyaaaf28
xyzzyaaag28=which_ssingle(xyzzyaaae28,spin_dep_in)
do xyzzyaaac28=1,nband(1,xyzzyaaaf28)
rck(xyzzyaaac28,1:num_ao,1,xyzzyaaag28)=xyzzyaabk1(xyzzyaaad28:xyzzyaa&
&ad28+num_ao-1)
xyzzyaaad28=xyzzyaaad28+num_ao
enddo
enddo
endif
if(xyzzyaabn1)then
exponent(:)=xyzzyaabc1(:)
xyzzyaaab28=1
do xyzzyaaaa28=1,num_shells
min_exponent(xyzzyaaaa28)=minval(exponent(xyzzyaaab28:xyzzyaaab28+xyzz&
&yaaad1(xyzzyaaaa28)-1))
xyzzyaaab28=xyzzyaaab28+xyzzyaaad1(xyzzyaaaa28)
enddo
endif
if(xyzzyaabo1)then
c_prim(:)=xyzzyaabd1(1:xyzzyaaac1)
if(any(shell_am(:)==2))c_prim2(:)=xyzzyaabd1(xyzzyaaac1+1:2*xyzzyaaac1&
&)
endif
end subroutine xyzzyaadc1
subroutine get_gaussian_orbmap(row_offset,norb,orbmap)
implicit none
integer,intent(inout) :: row_offset(nspin),norb,orbmap(nemax,nspin,nde&
&t)
integer xyzzyaaaa29
do xyzzyaaaa29=1,nspin
orbmap(row_offset(xyzzyaaaa29)+1:row_offset(xyzzyaaaa29)+nuc_nele(xyzz&
&yaaaa29),xyzzyaaaa29,:)=xyzzyaaar1(1:nuc_nele(xyzzyaaaa29),xyzzyaaaa2&
&9,:)
row_offset(xyzzyaaaa29)=row_offset(xyzzyaaaa29)+nuc_nele(xyzzyaaaa29)
enddo
norb=norb+xyzzyaaap1
end subroutine get_gaussian_orbmap
subroutine get_gaussian_ndesc(ndesc_int,ndesc_dp)
implicit none
integer,intent(inout) :: ndesc_int,ndesc_dp
ndesc_int=xyzzyaaci1
ndesc_dp=xyzzyaacj1
end subroutine get_gaussian_ndesc
subroutine get_gaussian_orbdesc(norb,ndesc_int,ndesc_dp,orbdesc_int,or&
&bdesc_dp)
implicit none
integer,intent(in) :: norb,ndesc_int,ndesc_dp
integer,intent(inout) :: orbdesc_int(ndesc_int,norb)
real(dp),intent(inout) :: orbdesc_dp(ndesc_dp,norb)
if(isperiodic)then
orbdesc_int(1:xyzzyaaci1,1:xyzzyaaap1)=xyzzyaack1(1:xyzzyaaci1,1:xyzzy&
&aaap1)
orbdesc_dp(1:xyzzyaacj1,1:xyzzyaaap1)=xyzzyaacl1(1:xyzzyaacj1,1:xyzzya&
&aap1)
endif
end subroutine get_gaussian_orbdesc
real(dp) function get_gaussian_rmax()
use slaarnabg,only : nitot,isperiodic
use run_control,only : errstop_master
implicit none
real(dp) xyzzyaaaa32
if(nitot/=1.or.isperiodic)call errstop_master('GET_GAUSSIAN_RMAX','Thi&
&s is not a single-atom system.')
xyzzyaaaa32=maxval(min_exponent)
get_gaussian_rmax=sqrt(screening_tolerance/xyzzyaaaa32)
end function get_gaussian_rmax
subroutine deshalloc_gauss_shm
implicit none
if(associated(xyzzyaaay1))call deshallocate(xyzzyaaay1)
if(associated(rck))call deshallocate(rck)
if(associated(cck))call deshallocate(cck)
end subroutine deshalloc_gauss_shm
end module slaarnaau

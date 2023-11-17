module slaarnacb
use dsp
use slaarnaag,   only : czero,zi,one_over_twopi,c_one
use file_utils,  only : skip,open_units
use format_utils,only : wout,d2s,i2s,r2s,wordwrap
use slaarnabg,    only : nbasis,atno,basis,pa1,pa2,pa3,periodicity,pb1&
&1,pb12,pb13,pb22,pb23,pb33,pb1,pb2,pb3,a1,a2,a3,npcells,model_system,&
&sr_lattice,num_g,rion,painv,b1,b2,b3,pbinv
use run_control, only : errstop,errwarn,errwarn_silent,check_alloc
use store,       only : nemax,ndet,real1_complex2,complex_wf,nspin,nuc&
&_nele,open_unit,use_backflow,orb_norm,no_ssingles
implicit none
private
public readpwf,pwfdet_setup,pw_orb_eval,get_pwfdet_orbmap,get_pwfdet_o&
&rbdesc,get_pwfdet_ndesc
integer xyzzyaaaa1
integer,allocatable :: xyzzyaaab1(:,:)
complex(dp),allocatable :: xyzzyaaac1(:,:),xyzzyaaad1(:),xyzzyaaae1(:)
integer xyzzyaaaf1
real(dp),allocatable :: xyzzyaaag1(:,:)
real(dp),allocatable :: xyzzyaaah1(:,:,:),xyzzyaaai1(:,:,:)
complex(dp),allocatable :: xyzzyaaaj1(:,:,:),xyzzyaaak1(:,:,:)
real(dp),allocatable :: xyzzyaaal1(:,:,:),xyzzyaaam1(:,:,:),xyzzyaaan1&
&(:,:,:),xyzzyaaao1(:,:,:),xyzzyaaap1(:,:,:),xyzzyaaaq1(:,:,:),xyzzyaa&
&ar1(:,:,:),xyzzyaaas1(:,:,:),xyzzyaaat1(:,:,:),xyzzyaaau1(:,:,:),xyzz&
&yaaav1(:,:,:),xyzzyaaaw1(:,:,:),xyzzyaaax1(:,:,:),xyzzyaaay1(:,:,:),x&
&yzzyaaaz1(:,:,:),xyzzyaaba1(:,:,:),xyzzyaabb1(:,:,:),xyzzyaabc1(:,:,:&
&),xyzzyaabd1(:,:,:),xyzzyaabe1(:,:,:)
complex(dp),allocatable :: xyzzyaabf1(:,:,:),xyzzyaabg1(:,:,:),xyzzyaa&
&bh1(:,:,:),xyzzyaabi1(:,:,:),xyzzyaabj1(:,:,:),xyzzyaabk1(:,:,:),xyzz&
&yaabl1(:,:,:),xyzzyaabm1(:,:,:),xyzzyaabn1(:,:,:),xyzzyaabo1(:,:,:),x&
&yzzyaabp1(:,:,:),xyzzyaabq1(:,:,:),xyzzyaabr1(:,:,:),xyzzyaabs1(:,:,:&
&),xyzzyaabt1(:,:,:),xyzzyaabu1(:,:,:),xyzzyaabv1(:,:,:),xyzzyaabw1(:,&
&:,:),xyzzyaabx1(:,:,:),xyzzyaaby1(:,:,:)
logical xyzzyaabz1
logical xyzzyaaca1
real(dp) xyzzyaacb1
integer xyzzyaacc1,xyzzyaacd1
integer,allocatable :: xyzzyaace1(:,:)
integer,allocatable :: xyzzyaacf1(:,:,:,:)
real(dp) xyzzyaacg1
real(dp),allocatable :: xyzzyaach1(:,:,:,:),xyzzyaaci1(:,:,:,:),xyzzya&
&acj1(:,:,:,:)
real(dp),allocatable :: xyzzyaack1(:,:,:)
real(dp),allocatable :: xyzzyaacl1(:,:)
logical xyzzyaacm1
logical,allocatable :: xyzzyaacn1(:)
logical,allocatable :: xyzzyaaco1(:)
complex(dp),allocatable :: xyzzyaacp1(:),xyzzyaacq1(:,:),xyzzyaacr1(:)&
&,xyzzyaacs1(:,:)
complex(dp),allocatable :: xyzzyaact1(:),xyzzyaacu1(:,:),xyzzyaacv1(:)&
&,xyzzyaacw1(:,:)
integer xyzzyaacx1,xyzzyaacy1,xyzzyaacz1
integer,allocatable :: xyzzyaada1(:,:,:)
integer,allocatable :: xyzzyaadb1(:)
integer,allocatable :: xyzzyaadc1(:),xyzzyaadd1(:),xyzzyaade1(:)
complex(dp),allocatable :: xyzzyaadf1(:)
integer,allocatable :: xyzzyaadg1(:)
integer,parameter :: xyzzyaadh1=1,xyzzyaadi1=2,xyzzyaadj1=3
logical,allocatable :: xyzzyaadk1(:,:),xyzzyaadl1(:,:)
integer,parameter :: xyzzyaadm1=2,xyzzyaadn1=1
integer,allocatable :: xyzzyaado1(:,:)
real(dp),allocatable :: xyzzyaadp1(:,:)
contains
subroutine readpwf(eionion)
use parallel
use slaarnaag, only : periodic_table_nocap
use slaarnaca, only : have_ae
implicit none
real(dp),intent(out) :: eionion
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,io,xyzzyaaad2,xyzzyaaae2,xyzz&
&yaaaf2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaa&
&al2,xyzzyaaam2,xyzzyaaan2,xyzzyaaao2
integer,allocatable :: xyzzyaaap2(:),xyzzyaaaq2(:)
logical xyzzyaaar2
logical,allocatable :: xyzzyaaas2(:)
real(dp) xyzzyaaat2,xyzzyaaau2,xyzzyaaav2,xyzzyaaaw2,xyzzyaaax2
character(80) title,code,method,functional,pseudo_type,sline,tmpr
character(20) psp_filename,ppchar
xyzzyaaca1=complex_wf
if(nspin/=2.and.am_master)call errstop('READPWF','Plane-wave orbitals &
&can only be used for electron systems at present.')
if(am_master)then
call open_units(io,xyzzyaaak2)
if(xyzzyaaak2/=0)call errstop('READPWF','Unable to find free i/o unit.&
&')
open(io,file='pwfn.data',status='old',action='read',err=10)
call skip(io,15)
read(io,*,end=20,err=30)xyzzyaacm1
call skip(io,18)
read(io,*,end=20,err=30)nbasis
call skip(io,nbasis+9)
read(io,*,end=20,err=30)xyzzyaaaf1
call skip(io,xyzzyaaaf1+5)
read(io,*,end=20,err=30)xyzzyaacc1
allocate(atno(nbasis),basis(3,nbasis),xyzzyaaag1(3,xyzzyaaaf1),xyzzyaa&
&cl1(3,xyzzyaacc1),xyzzyaace1(xyzzyaacc1,2),stat=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'READPWF','')
if(xyzzyaacm1)then
xyzzyaaai2=2
else
xyzzyaaai2=1
endif
do xyzzyaaac2=1,xyzzyaacc1
call skip(io,1)
read(io,*,end=20,err=30)xyzzyaaag2,xyzzyaace1(xyzzyaaac2,1:2)
do xyzzyaaah2=1,xyzzyaaai2
do xyzzyaaad2=1,xyzzyaace1(xyzzyaaac2,xyzzyaaah2)
call skip(io,3)
if(xyzzyaaac2==1.and.xyzzyaaad2==1)then
read(io,fmt='(a)',end=20,err=30)sline
if(scan(sline,",")/=0)then
xyzzyaabz1=.false.
else
xyzzyaabz1=.true.
endif
backspace(io)
endif
call skip(io,xyzzyaaaf1)
enddo
enddo
enddo
xyzzyaacd1=maxval(xyzzyaace1)
allocate(xyzzyaack1(xyzzyaacd1,xyzzyaacc1,2),stat=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'READPWF','Eigenvalue')
if(xyzzyaabz1)then
allocate(xyzzyaaah1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'READPWF','CKG')
xyzzyaaah1=0.d0
if(xyzzyaacm1)then
allocate(xyzzyaaai1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'READPWF','CKG2')
xyzzyaaai1=0.d0
endif
else
allocate(xyzzyaaaj1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'READPWF','CCKG')
xyzzyaaaj1=czero
if(xyzzyaacm1)then
allocate(xyzzyaaak1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'READPWF','CCKG2')
xyzzyaaak1=czero
endif
endif
rewind(io)
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
read(io,*,err=30,end=20)xyzzyaaat2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaacm1
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaau2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaacb1
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaav2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaaw2
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaax2
call skip(io,1)
read(io,*,err=30,end=20)eionion
call skip(io,1)
read(io,*,err=30,end=20)xyzzyaaaj2
call skip(io,4)
read(io,*,end=20,err=30)nbasis
call skip(io,1)
do xyzzyaaaa2=1,nbasis
read(io,*,end=20,err=30)atno(xyzzyaaaa2),basis(1:3,xyzzyaaaa2)
enddo
call skip(io,1)
read(io,*,end=20,err=30)pa1
read(io,*,end=20,err=30)pa2
read(io,*,end=20,err=30)pa3
call skip(io,4)
read(io,*,end=20,err=30)xyzzyaaaf1
call skip(io,1)
do xyzzyaaae2=1,xyzzyaaaf1
read(io,*,end=20,err=30)xyzzyaaag1(1:3,xyzzyaaae2)
enddo
call skip(io,4)
read(io,*,end=20,err=30)xyzzyaacc1
do xyzzyaaac2=1,xyzzyaacc1
call skip(io,1)
read(io,*,end=20,err=30)xyzzyaaag2,xyzzyaace1(xyzzyaaac2,1:2),xyzzyaac&
&l1(1:3,xyzzyaaac2)
do xyzzyaaah2=1,xyzzyaaai2
do xyzzyaaad2=1,xyzzyaace1(xyzzyaaac2,xyzzyaaah2)
call skip(io,1)
read(io,*,end=20,err=30)xyzzyaaaa2,xyzzyaaab2,xyzzyaack1(xyzzyaaad2,xy&
&zzyaaac2,xyzzyaaah2)
call skip(io,1)
if(xyzzyaabz1)then
if(xyzzyaaah2==1)then
do xyzzyaaae2=1,xyzzyaaaf1
read(io,*,end=20,err=30)xyzzyaaah1(xyzzyaaae2,xyzzyaaad2,xyzzyaaac2)
enddo
else
do xyzzyaaae2=1,xyzzyaaaf1
read(io,*,end=20,err=30)xyzzyaaai1(xyzzyaaae2,xyzzyaaad2,xyzzyaaac2)
enddo
endif
else
if(xyzzyaaah2==1)then
do xyzzyaaae2=1,xyzzyaaaf1
read(io,*,end=20,err=30)xyzzyaaaj1(xyzzyaaae2,xyzzyaaad2,xyzzyaaac2)
enddo
else
do xyzzyaaae2=1,xyzzyaaaf1
read(io,*,end=20,err=30)xyzzyaaak1(xyzzyaaae2,xyzzyaaad2,xyzzyaaac2)
enddo
endif
endif
enddo
enddo
enddo
if(.not.xyzzyaacm1)then
xyzzyaack1(:,:,2)=xyzzyaack1(:,:,1)
xyzzyaace1(:,2)=xyzzyaace1(:,1)
endif
close(io)
open_unit(io)=.false.
if(any(xyzzyaaag1(3,:)/=0.d0))then
periodicity=3
else
if(any(xyzzyaaag1(2,:)/=0.d0))then
periodicity=2
else
if(any(xyzzyaaag1(1,:)/=0.d0))then
periodicity=1
else
if(xyzzyaaaf1>0)call errstop('READPWF','All G appear to be zero.')
endif
endif
endif
if(periodicity==1)call errstop('READPWF','Deduced one-dimensional peri&
&odicity from G vectors in pwfn.data file. Not currently allowed.')
call wout()
call wout('Title: '//trim(title))
call wout('Generating code                           :  '//trim(adjust&
&l(code)))
call wout('Method                                    :  '//trim(adjust&
&l(method)))
call wout('DFT functional                            :  '//trim(adjust&
&l(functional)))
call wout('Pseudopotential type                      :  '//trim(adjust&
&l(pseudo_type)))
tmpr=r2s(xyzzyaaat2,'(f12.3)')
call wout('Plane-wave cutoff (au)                    :  '//trim(tmpr))
call wout('Deduced periodicity                       :  '//trim(i2s(pe&
&riodicity)))
call wout()
call wout('Number of k points                        :  '//trim(i2s(xy&
&zzyaacc1)))
call wout('Max # bands per k point                   :  '//trim(i2s(xy&
&zzyaacd1)))
call wout('Number of G vectors                       :  '//trim(i2s(xy&
&zzyaaaf1)))
call wout()
if(xyzzyaaau2/=0.d0.or.xyzzyaacb1/=0.d0)then
call wout('DFT energy and components (au per primitive cell):')
tmpr=r2s(xyzzyaaau2,'(f16.10)')
call wout('Total energy                              :  '//trim(tmpr))
tmpr=r2s(xyzzyaacb1,'(f16.10)')
call wout('Kinetic energy                            :  '//trim(tmpr))
tmpr=r2s(xyzzyaaav2,'(f16.10)')
call wout('Local potential energy                    :  '//trim(tmpr))
tmpr=r2s(xyzzyaaaw2,'(f16.10)')
call wout('Non-local potential energy                :  '//trim(tmpr))
tmpr=r2s(xyzzyaaax2,'(f16.10)')
call wout('Electron-electron energy                  :  '//trim(tmpr))
tmpr=r2s(eionion,'(f16.10)')
call wout('Ion-ion energy                            :  '//trim(tmpr))
call wout()
else
call wout('Information about DFT energies apparently not supplied.')
call wout()
endif
if(xyzzyaabz1)then
call wout('Real orbital coefficients ==> inversion symmetry')
else
call wout('Complex orbital coefficients ==> no inversion symmetry')
endif
if(xyzzyaaca1)then
call wout('Complex plane-wave orbitals will be evaluated.')
if(.not.complex_wf)call errstop('READPWF','Complex orbitals to be retu&
&rned, but trial wave function is real.')
else
call wout('Real plane-wave orbitals will be evaluated.')
if(complex_wf)call errwarn('READPWF','Real orbitals to be returned,but&
& trial wave function is complex.')
endif
if(any(atno(:)>200))then
allocate(xyzzyaaas2(nbasis),xyzzyaaap2(nbasis),xyzzyaaaq2(nbasis),stat&
&=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'READPWF','is_pp_from_atno, etc..')
xyzzyaaap2=0
xyzzyaaaq2=0
xyzzyaaas2(:)=atno(:)>200
do xyzzyaaam2=1,nbasis
atno(xyzzyaaam2)=mod(atno(xyzzyaaam2),100)
enddo
xyzzyaaan2=0
have_ae=.false.
a: do xyzzyaaam2=1,nbasis
do xyzzyaaac2=1,xyzzyaaam2-1
if(xyzzyaaap2(xyzzyaaac2)==atno(xyzzyaaam2))then
xyzzyaaap2(xyzzyaaam2)=xyzzyaaap2(xyzzyaaac2)
xyzzyaaaq2(xyzzyaaam2)=xyzzyaaaq2(xyzzyaaac2)
cycle a
endif
enddo
xyzzyaaan2=xyzzyaaan2+1
xyzzyaaap2(xyzzyaaam2)=atno(xyzzyaaam2)
xyzzyaaaq2(xyzzyaaam2)=xyzzyaaan2
enddo a
xyzzyaaao2=xyzzyaaan2
do xyzzyaaal2=1,xyzzyaaao2
do xyzzyaaam2=1,nbasis
if(xyzzyaaaq2(xyzzyaaam2)==xyzzyaaal2)exit
enddo
xyzzyaaaa2=mod(atno(xyzzyaaam2),1000)
xyzzyaaab2=atno(xyzzyaaam2)/1000+1
ppchar=''
if(xyzzyaaab2>1)ppchar=trim(i2s(xyzzyaaab2))
psp_filename=trim(periodic_table_nocap(xyzzyaaaa2))//trim(adjustl(ppch&
&ar))//'_pp.data'
inquire(file=psp_filename,exist=xyzzyaaar2,err=7)
if(xyzzyaaas2(xyzzyaaam2).and..not.xyzzyaaar2)call errstop('READPWF','&
&The atomic numbers in this pwfn.data file follow the +200 convention &
&to flag pseudoatoms. However, the x_pp.data pseudopotential file for &
&a flagged pseudoatom is missing.')
if(.not.xyzzyaaar2)have_ae=.true.
enddo
deallocate(xyzzyaaas2,xyzzyaaaq2,xyzzyaaap2)
endif
endif
call mpi_bcast(eionion,1,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting eionion in readpwf.')
call mpi_bcast(pa1,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa1 in readpwf.')
call mpi_bcast(pa2,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa2 in readpwf.')
call mpi_bcast(pa3,3,mpi_double_precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting pa3 in readpwf.')
call mpi_bcast(periodicity,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting periodicity in readpwf.')
call mpi_bcast(xyzzyaaai2,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting num_spins in readpwf.')
call mpi_bcast(nbasis,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nbasis in readpwf.')
call mpi_bcast(xyzzyaaaf1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ngvec in readpwf.')
call mpi_bcast(xyzzyaacc1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting nkvec in readpwf.')
call mpi_bcast(xyzzyaacd1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting maxband in readpwf.')
call mpi_bcast(xyzzyaabz1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting real_pw_coeffs in readpwf.')
call mpi_bcast(xyzzyaacm1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting spin_polarized in readpwf.')
call mpi_bcast(have_ae,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting have_ae in readpwf.')
if(am_slave)then
allocate(atno(nbasis),basis(3,nbasis),stat=xyzzyaaaf2)
call check_alloc(xyzzyaaaf2,'READPWF','basis [slaves]')
if(xyzzyaabz1)then
allocate(xyzzyaaag1(3,xyzzyaaaf1),xyzzyaacl1(3,xyzzyaacc1),xyzzyaace1(&
&xyzzyaacc1,2),xyzzyaack1(xyzzyaacd1,xyzzyaacc1,2),xyzzyaaah1(xyzzyaaa&
&f1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaf2)
else
allocate(xyzzyaaag1(3,xyzzyaaaf1),xyzzyaacl1(3,xyzzyaacc1),xyzzyaace1(&
&xyzzyaacc1,2),xyzzyaack1(xyzzyaacd1,xyzzyaacc1,2),xyzzyaaaj1(xyzzyaaa&
&f1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaf2)
endif
call check_alloc(xyzzyaaaf2,'READPWF','wavefun [slaves]')
if(xyzzyaacm1)then
if(xyzzyaabz1)then
allocate(xyzzyaaai1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaf2)
else
allocate(xyzzyaaak1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaf2)
endif
call check_alloc(xyzzyaaaf2,'READPWF','wavefun <2> [slaves]')
endif
endif
call mpi_bcast(basis,nbasis*3,mpi_double_precision,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting basis in readpwf.')
call mpi_bcast(atno,nbasis,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting atno in readpwf.')
call mpi_bcast(xyzzyaaag1,xyzzyaaaf1*3,mpi_double_precision,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting gvecwf in readpwf.')
call mpi_bcast(xyzzyaacl1,xyzzyaacc1*3,mpi_double_precision,0,mpi_comm&
&_world,ierror)
call checkmpi(ierror,'Broadcasting kvec in readpwf.')
call mpi_bcast(xyzzyaace1,xyzzyaacc1*2,mpi_integer,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting nband in readpwf.')
call mpi_bcast(xyzzyaack1,2*xyzzyaacc1*xyzzyaacd1,mpi_double_precision&
&,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting eigenvalue in readpwf.')
if(xyzzyaabz1)then
call mpi_bcast(xyzzyaaah1,xyzzyaacc1*xyzzyaacd1*xyzzyaaaf1,mpi_double_&
&precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ckg in readpwf.')
if(xyzzyaacm1)then
call mpi_bcast(xyzzyaaai1,xyzzyaacc1*xyzzyaacd1*xyzzyaaaf1,mpi_double_&
&precision,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ckg2 in readpwf.')
endif
else
call mpi_bcast(xyzzyaaaj1,xyzzyaacc1*xyzzyaacd1*xyzzyaaaf1,mpi_double_&
&complex,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting cckg in readpwf.')
if(xyzzyaacm1)then
call mpi_bcast(xyzzyaaak1,xyzzyaacc1*xyzzyaacd1*xyzzyaaaf1,mpi_double_&
&complex,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting cckg2 in readpwf.')
endif
endif
model_system=all(atno==0)
return
7  call errstop('READPWF','Error inquiring about '//trim(psp_filename)&
&//'.')
10 call errstop('READPWF','Cannot find pwfn.data file.')
20 call errstop('READPWF','Read past end of pwfn.data file.')
30 call errstop('READPWF','Error reading pwfn.data file.')
end subroutine readpwf
subroutine pwfdet_setup(kwarn)
use slaarnaag,only : pi_over_two
use slaarnaan, only : check_kpoints
use slaarnabi, only : use_gpcc
use slaarnabp, only : no_orb_phases,orb_phase_band,orb_phase_kpoint,or&
&b_phase_spin,orb_phase_det,orb_phase,wf_np,wf_nm,wf_nd,wf_d
use slaarnabt, only : inverse3,quicksort,get_numerical_orbmask
use parallel, only  : am_master
implicit none
logical,intent(in) :: kwarn
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3,i,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3,xyzzyaaa&
&l3,xyzzyaaam3,xyzzyaaan3,xyzzyaaao3,xyzzyaaap3,xyzzyaaaq3,xyzzyaaar3,&
&xyzzyaaas3,xyzzyaaat3,xyzzyaaau3,xyzzyaaav3,xyzzyaaaw3,xyzzyaaax3,xyz&
&zyaaay3,xyzzyaaaz3,xyzzyaaba3,xyzzyaabb3,xyzzyaabc3
integer,allocatable :: xyzzyaabd3(:),xyzzyaabe3(:),xyzzyaabf3(:,:),xyz&
&zyaabg3(:,:),xyzzyaabh3(:)
real(dp) xyzzyaabi3(3),xyzzyaabj3,xyzzyaabk3,xyzzyaabl3,xyzzyaabm3(3),&
&xyzzyaabn3(1:3),xyzzyaabo3,xyzzyaabp3(3,3),xyzzyaabq3(3,3),xyzzyaabr3&
&,xyzzyaabs3,xyzzyaabt3,xyzzyaabu3,xyzzyaabv3,xyzzyaabw3,xyzzyaabx3,xy&
&zzyaaby3,xyzzyaabz3
real(dp),parameter :: xyzzyaaca3=1.d-8
real(dp),allocatable :: xyzzyaacb3(:),xyzzyaacc3(:,:),xyzzyaacd3(:,:),&
&xyzzyaace3(:,:)
complex(dp) xyzzyaacf3,xyzzyaacg3
complex(dp),allocatable :: xyzzyaach3(:)
logical xyzzyaaci3,xyzzyaacj3,xyzzyaack3
logical,parameter :: xyzzyaacl3=.true.
logical,allocatable :: ltemp(:),xyzzyaacm3(:)
character(80) char80,char80_new
character (256) tmpr
if(xyzzyaaca1)then
allocate(xyzzyaacf1(ndet,xyzzyaacd1,xyzzyaacc1,nspin),xyzzyaacn1(xyzzy&
&aacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','0')
xyzzyaacf1=0
xyzzyaacn1=.true.
if(xyzzyaacm1)then
xyzzyaaae3=2
else
xyzzyaaae3=1
endif
call check_kpoints(xyzzyaacc1,xyzzyaacl1)
if(xyzzyaacm1)then
i=max(sum(xyzzyaace1(:,1)),sum(xyzzyaace1(:,2)))
else
i=sum(xyzzyaace1(:,1))
endif
allocate(xyzzyaacb3(i),xyzzyaabe3(i),xyzzyaabd3(i),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','EIGTEMP')
do xyzzyaaao3=1,nspin
if(xyzzyaaao3==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
if(xyzzyaaao3==1.or.xyzzyaaad3==2)then
xyzzyaaag3=0
do xyzzyaaab3=1,xyzzyaacc1
do xyzzyaaaf3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
xyzzyaaag3=xyzzyaaag3+1
xyzzyaacb3(xyzzyaaag3)=xyzzyaack1(xyzzyaaaf3,xyzzyaaab3,xyzzyaaad3)
xyzzyaabe3(xyzzyaaag3)=xyzzyaaab3
enddo
enddo
if(xyzzyaaag3>0)call quicksort(xyzzyaaag3,xyzzyaacb3(1),xyzzyaabd3(1))
endif
if(nuc_nele(xyzzyaaao3)>xyzzyaaag3)call errstop('PWFDET_SETUP','An ins&
&ufficient number of eigenstates have been supplied.')
do i=1,nuc_nele(xyzzyaaao3)
xyzzyaaab3=xyzzyaabe3(xyzzyaabd3(i))
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(1,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==0)then
xyzzyaacf1(1,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)=1
exit
endif
enddo
enddo
if(am_master.and.nuc_nele(xyzzyaaao3)<xyzzyaaag3.and.nuc_nele(xyzzyaaa&
&o3)>0)then
if(abs(xyzzyaacb3(xyzzyaabd3(nuc_nele(xyzzyaaao3)))-xyzzyaacb3(xyzzyaa&
&bd3(nuc_nele(xyzzyaaao3)+1)))<xyzzyaaca3*abs(xyzzyaacb3(xyzzyaabd3(nu&
&c_nele(xyzzyaaao3)))).and.ndet==1)call errwarn('PWFDET_SETUP','Partia&
&lly occupied degenerate states at the Fermi level for spin '//trim(i2&
&s(xyzzyaaao3))//'. Multideterminant calculation probably advisable.')
endif
do xyzzyaaap3=2,ndet
xyzzyaacf1(xyzzyaaap3,:,:,xyzzyaaao3)=xyzzyaacf1(1,:,:,xyzzyaaao3)
enddo
enddo
if(am_master)then
call wout()
call wout('Plane-wave setup')
call wout('================')
call wout()
call wout('Periodicity : '//trim(i2s(periodicity)))
call wout('Complex plane-wave orbitals will be used.')
call wout()
if(.not.complex_wf)call errstop('PWFDET_SETUP','PW routine requested t&
&o return complex orbitals, but CASINO is not using a complex wave fun&
&ction.')
if(periodicity>0)then
allocate(xyzzyaabf3(xyzzyaacc1,nspin),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','boccband')
xyzzyaabf3=0
do xyzzyaaao3=1,nspin
if(xyzzyaaao3==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
do xyzzyaaab3=1,xyzzyaacc1
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(1,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)>0)xyzzyaabf3(xyzzyaa&
&ab3,xyzzyaaao3)=xyzzyaabf3(xyzzyaaab3,xyzzyaaao3)+1
enddo
enddo
enddo
if((xyzzyaacm1.or.nuc_nele(1)==nuc_nele(2)).and.all(wf_nd(1,1:nspin)==&
&0))then
xyzzyaaci3=.false.
sc:  do xyzzyaaao3=1,nspin
do xyzzyaaab3=2,xyzzyaacc1
if(xyzzyaabf3(xyzzyaaab3,xyzzyaaao3)/=xyzzyaabf3(1,xyzzyaaao3))then
xyzzyaaci3=.true.
exit sc
endif
enddo
enddo sc
if(xyzzyaaci3)then
call wout('METALLIC GROUND STATE DETECTED')
if(.not.xyzzyaacm1)then
call wout('Number of doubly-occupied bands filled at each k point in G&
&S:')
call wout('k-point    no. of filled bands')
do xyzzyaaab3=1,xyzzyaacc1
write(tmpr,'(i4,4x,i5)')xyzzyaaab3,xyzzyaabf3(xyzzyaaab3,1)
call wout(tmpr)
enddo
else
call wout('Number of singly-occupied bands filled at each k point in G&
&S :')
call wout('Spin    k-point    no. of filled bands')
do xyzzyaaab3=1,xyzzyaacc1
write(tmpr,'(a8,i4,4x,i5)')'up',xyzzyaaab3,xyzzyaabf3(xyzzyaaab3,1)
call wout(tmpr)
enddo
do xyzzyaaab3=1,xyzzyaacc1
write(tmpr,'(a8,i4,4x,i5)')'down',xyzzyaaab3,xyzzyaabf3(xyzzyaaab3,2)
call wout(tmpr)
enddo
endif
else
call wout('INSULATING GROUND STATE DETECTED')
if(.not.xyzzyaacm1)then
call wout('No. of doubly occupied bands at each k point: '//trim(i2s(x&
&yzzyaabf3(1,1))))
else
call wout('No. of occupied spin-up bands at each k point  : '//trim(i2&
&s(xyzzyaabf3(1,1))))
call wout('No. of occupied spin-down bands at each k point: '//trim(i2&
&s(xyzzyaabf3(1,2))))
endif
endif
call wout()
endif
deallocate(xyzzyaabf3)
endif
do xyzzyaaab3=1,xyzzyaacc1
if(periodicity<3.and.abs(dot_product(xyzzyaacl1(1:3,xyzzyaaab3),pa3(1:&
&3)))>xyzzyaaca3)call errstop('PWFDET_SETUP','k vectors inappropriate &
&for periodicity.')
if(periodicity<2.and.abs(dot_product(xyzzyaacl1(1:3,xyzzyaaab3),pa2(1:&
&3)))>xyzzyaaca3)call errstop('PWFDET_SETUP','k vectors inappropriate &
&for periodicity.')
if(periodicity<1.and.abs(dot_product(xyzzyaacl1(1:3,xyzzyaaab3),pa1(1:&
&3)))>xyzzyaaca3)call errstop('PWFDET_SETUP','k vectors inappropriate &
&for periodicity.')
enddo
endif
do xyzzyaaao3=1,nspin
if(xyzzyaaao3==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
do xyzzyaaap3=1,ndet
if((wf_np(xyzzyaaap3,xyzzyaaao3)>0.or.wf_nm(xyzzyaaap3,xyzzyaaao3)>0).&
&and.am_master)call errstop('PWFDET_SETUP','Additions or subtractions &
&should not be used for plane-wave orbitals: number of electrons is si&
&mply determined by NEU and NED, and occupancy should only be changed &
&using "PR" in the MDET block.')
do xyzzyaaah3=1,wf_nd(xyzzyaaap3,xyzzyaaao3)
xyzzyaaak3=wf_d(1,xyzzyaaah3,xyzzyaaap3,xyzzyaaao3)
xyzzyaaai3=wf_d(2,xyzzyaaah3,xyzzyaaap3,xyzzyaaao3)
xyzzyaaal3=wf_d(3,xyzzyaaah3,xyzzyaaap3,xyzzyaaao3)
xyzzyaaaj3=wf_d(4,xyzzyaaah3,xyzzyaaap3,xyzzyaaao3)
if(xyzzyaaaj3<1.or.xyzzyaaaj3>xyzzyaacc1.or.xyzzyaaai3<1.or.xyzzyaaai3&
&>xyzzyaacc1)call errstop('PWFDET_SETUP','k-point out of range in exci&
&tation specification.')
if(xyzzyaaal3<1.or.xyzzyaaal3>xyzzyaace1(xyzzyaaaj3,xyzzyaaad3).or.xyz&
&zyaaak3<1.or.xyzzyaaak3>xyzzyaace1(xyzzyaaai3,xyzzyaaad3))call errsto&
&p('PWFDET_SETUP','band index out of range in excitation specification&
&.')
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaal3,xyzzyaaaj3,xyzzyaaao3)/=0)call err&
&stop('PWFDET_SETUP','Trying to promote an electron into a band that i&
&s already fully occupied.')
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaak3,xyzzyaaai3,xyzzyaaao3)/=1)call err&
&stop('PWFDET_SETUP','Trying to promote an electron from an empty band&
&.')
xyzzyaacf1(xyzzyaaap3,xyzzyaaal3,xyzzyaaaj3,xyzzyaaao3)=1
xyzzyaacf1(xyzzyaaap3,xyzzyaaak3,xyzzyaaai3,xyzzyaaao3)=0
enddo
i=0
do xyzzyaaab3=1,xyzzyaacc1
i=i+sum(xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3))
enddo
if(i/=nuc_nele(xyzzyaaao3))call errstop('PWFDET_SETUP','Problem with e&
&xcitations (A).')
if(any(xyzzyaacf1(xyzzyaaap3,:,:,xyzzyaaao3)<0).or.any(xyzzyaacf1(xyzz&
&yaaap3,:,:,xyzzyaaao3)>1))call errstop('PWFDET_SETUP','Problem with e&
&xcitations (B).')
enddo
enddo
deallocate(xyzzyaabd3,xyzzyaacb3,xyzzyaabe3)
if(am_master.and.xyzzyaacl3)then
call wout('Detailed band-occupancy information')
call wout('-----------------------------------')
do xyzzyaaab3=1,xyzzyaacc1
call wout('K POINT '//trim(i2s(xyzzyaaab3)))
write(tmpr,'("k(au):",3(1x,es23.15))')xyzzyaacl1(1:3,xyzzyaaab3)
call wout(tmpr)
write(tmpr,'("kfrac:",3(1x,es23.15))')matmul(xyzzyaacl1(1:3,xyzzyaaab3&
&),pbinv)
call wout(tmpr)
do xyzzyaaao3=1,nspin
if(nuc_nele(xyzzyaaao3)==0)cycle
if(xyzzyaacm1.and.xyzzyaaao3==2)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
do xyzzyaaap3=1,ndet
if(ndet>1)then
call wout('Spin '//trim(i2s(xyzzyaaao3))//', determinant: '//trim(i2s(&
&xyzzyaaap3))//':')
else
call wout('Spin '//trim(i2s(xyzzyaaao3))//':')
endif
if(any(xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==1))then
char80=' Bands occupied at k:'
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==1)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam3))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam3))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Highest occupied band energy: ',maxval(xyzzyaack1(:,xyzzyaa&
&ab3,xyzzyaaad3),xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==1))
endif
if(any(xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==0))then
char80=' Bands unoccupied at k:'
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==0)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam3))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam3))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Lowest unoccupied band energy: ',minval(xyzzyaack1(:,xyzzya&
&aab3,xyzzyaaad3),xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==0))
endif
enddo
enddo
call wout()
enddo
endif
else
allocate(xyzzyaacn1(xyzzyaacc1),xyzzyaacm3(xyzzyaacc1),xyzzyaaco1(xyzz&
&yaacc1),xyzzyaacf1(ndet,xyzzyaacd1,xyzzyaacc1,nspin),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','2')
xyzzyaacf1=0
call check_kpoints(xyzzyaacc1,xyzzyaacl1)
if(xyzzyaacm1)then
xyzzyaaae3=2
else
xyzzyaaae3=1
endif
xyzzyaacn1=.true.
xyzzyaacm3=.false.
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
do xyzzyaaac3=xyzzyaaab3+1,xyzzyaacc1
xyzzyaabi3(1:3)=xyzzyaacl1(1:3,xyzzyaaab3)+xyzzyaacl1(1:3,xyzzyaaac3)
xyzzyaabo3=dot_product(xyzzyaabi3,pa1)*one_over_twopi
if(abs(xyzzyaabo3-anint(xyzzyaabo3))<xyzzyaaca3)then
xyzzyaabo3=dot_product(xyzzyaabi3,pa2)*one_over_twopi
if(abs(xyzzyaabo3-anint(xyzzyaabo3))<xyzzyaaca3)then
xyzzyaabo3=dot_product(xyzzyaabi3,pa3)*one_over_twopi
if(abs(xyzzyaabo3-anint(xyzzyaabo3))<xyzzyaaca3)then
xyzzyaacn1(xyzzyaaac3)=.false.
xyzzyaacm3(xyzzyaaab3)=.true.
xyzzyaacm3(xyzzyaaac3)=.true.
exit
endif
endif
endif
enddo
endif
enddo
if(any(.not.xyzzyaacn1))call errwarn_silent('PWFDET_SETUP','Redundant &
&paired k points in pwfn.data file.')
xyzzyaaco1=.false.
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
xyzzyaabi3(1:3)=2.d0*xyzzyaacl1(1:3,xyzzyaaab3)
xyzzyaabo3=dot_product(xyzzyaabi3,pa1)*one_over_twopi
if(abs(xyzzyaabo3-anint(xyzzyaabo3))<xyzzyaaca3)then
xyzzyaabo3=dot_product(xyzzyaabi3,pa2)*one_over_twopi
if(abs(xyzzyaabo3-anint(xyzzyaabo3))<xyzzyaaca3)then
xyzzyaabo3=dot_product(xyzzyaabi3,pa3)*one_over_twopi
if(abs(xyzzyaabo3-anint(xyzzyaabo3))<xyzzyaaca3)xyzzyaaco1(xyzzyaaab3)&
&=.true.
endif
endif
endif
enddo
if(xyzzyaacm1)then
i=2*max(sum(xyzzyaace1(:,1)),sum(xyzzyaace1(:,2)))
else
i=2*sum(xyzzyaace1(:,1))
endif
allocate(xyzzyaacb3(i),xyzzyaabe3(i),xyzzyaabd3(i),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','EIGTEMP')
do xyzzyaaao3=1,nspin
if(xyzzyaaao3==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
if(xyzzyaaao3==1.or.xyzzyaaad3==2)then
xyzzyaaag3=0
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
do xyzzyaaaf3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
xyzzyaaag3=xyzzyaaag3+1
xyzzyaacb3(xyzzyaaag3)=xyzzyaack1(xyzzyaaaf3,xyzzyaaab3,xyzzyaaad3)
xyzzyaabe3(xyzzyaaag3)=xyzzyaaab3
if(.not.xyzzyaaco1(xyzzyaaab3))then
xyzzyaaag3=xyzzyaaag3+1
xyzzyaacb3(xyzzyaaag3)=xyzzyaack1(xyzzyaaaf3,xyzzyaaab3,xyzzyaaad3)
xyzzyaabe3(xyzzyaaag3)=xyzzyaaab3
endif
enddo
endif
enddo
if(xyzzyaaag3>0)call quicksort(xyzzyaaag3,xyzzyaacb3(1),xyzzyaabd3(1))
endif
if(nuc_nele(xyzzyaaao3)>xyzzyaaag3)call errstop('PWFDET_SETUP','An ins&
&ufficient number of eigenstates have been supplied.')
do i=1,nuc_nele(xyzzyaaao3)
xyzzyaaab3=xyzzyaabe3(xyzzyaabd3(i))
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if((xyzzyaaco1(xyzzyaaab3).and.xyzzyaacf1(1,xyzzyaaam3,xyzzyaaab3,xyzz&
&yaaao3)<1).or.(.not.xyzzyaaco1(xyzzyaaab3).and.xyzzyaacf1(1,xyzzyaaam&
&3,xyzzyaaab3,xyzzyaaao3)<2))then
xyzzyaacf1(1,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)=xyzzyaacf1(1,xyzzyaaam3&
&,xyzzyaaab3,xyzzyaaao3)+1
exit
endif
enddo
enddo
if(am_master.and.nuc_nele(xyzzyaaao3)<xyzzyaaag3.and.nuc_nele(xyzzyaaa&
&o3)>0)then
if(abs(xyzzyaacb3(xyzzyaabd3(nuc_nele(xyzzyaaao3)))-xyzzyaacb3(xyzzyaa&
&bd3(nuc_nele(xyzzyaaao3)+1)))<xyzzyaaca3*abs(xyzzyaacb3(xyzzyaabd3(nu&
&c_nele(xyzzyaaao3)))).and.ndet==1)call errwarn('PWFDET_SETUP','Partia&
&lly occupied degenerate states at the Fermi level for spin '//trim(i2&
&s(xyzzyaaao3))//'. Multideterminant calculation probably advisable.')
endif
do xyzzyaaap3=2,ndet
xyzzyaacf1(xyzzyaaap3,:,:,xyzzyaaao3)=xyzzyaacf1(1,:,:,xyzzyaaao3)
enddo
enddo
if(am_master)then
call wout()
call wout('Plane-wave setup')
call wout('================')
call wout()
call wout('Periodicity : '//trim(i2s(periodicity)))
call wout('Real plane-wave orbitals will be used.')
call wout()
if(complex_wf)call errwarn('PWFDET_SETUP','PW routine requested to ret&
&urn real orbitals, but CASINO is using a complex wave function.')
if(periodicity>0)then
allocate(xyzzyaabf3(xyzzyaacc1,nspin),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','boccband')
xyzzyaabf3=0
do xyzzyaaao3=1,nspin
if(xyzzyaaao3==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(1,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)>0)xyzzyaabf3(xyzzyaa&
&ab3,xyzzyaaao3)=xyzzyaabf3(xyzzyaaab3,xyzzyaaao3)+1
enddo
endif
enddo
enddo
xyzzyaack3=(xyzzyaacm1.or.nuc_nele(1)==nuc_nele(2)).and.all(wf_nd(1,1:&
&nspin)==0)
if(xyzzyaack3)then
xyzzyaaci3=.false.
sr:  do xyzzyaaao3=1,nspin
do xyzzyaaaq3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaaq3))exit
enddo
do xyzzyaaab3=xyzzyaaaq3+1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
if(xyzzyaabf3(xyzzyaaab3,xyzzyaaao3)/=xyzzyaabf3(xyzzyaaaq3,xyzzyaaao3&
&))then
xyzzyaaci3=.true.
exit sr
endif
endif
enddo
enddo sr
endif
call wout('K POINT ANALYSIS')
call wout('  k    kx         ky         kz       use pair edge')
do xyzzyaaab3=1,xyzzyaacc1
write(tmpr,'(i3,1x,f10.6,1x,f10.6,1x,f10.6,3x,l1,4x,l1,4x,l1)')xyzzyaa&
&ab3,xyzzyaacl1(1,xyzzyaaab3),xyzzyaacl1(2,xyzzyaaab3),xyzzyaacl1(3,xy&
&zzyaaab3),xyzzyaacn1(xyzzyaaab3),xyzzyaacm3(xyzzyaaab3),xyzzyaaco1(xy&
&zzyaaab3)
call wout(tmpr)
enddo
call wout()
call wout('Any k points with edge=F give rise to independent states at&
& both k and -k.')
call wout()
allocate(xyzzyaacc3(3,xyzzyaacc1),ltemp(xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','kvec2')
xyzzyaacc3=1.d5
ltemp=.true.
do xyzzyaaac3=1,xyzzyaacc1
do xyzzyaaab3=1,num_g
xyzzyaabn3(1:3)=xyzzyaacl1(1:3,xyzzyaaac3)-0.5d0*sr_lattice(1:3,xyzzya&
&aab3)
xyzzyaabj3=one_over_twopi*dot_product(xyzzyaabn3,a1)
xyzzyaabk3=one_over_twopi*dot_product(xyzzyaabn3,a2)
xyzzyaabl3=one_over_twopi*dot_product(xyzzyaabn3,a3)
if(abs(xyzzyaabj3-anint(xyzzyaabj3))<xyzzyaaca3.and.abs(xyzzyaabk3-ani&
&nt(xyzzyaabk3))<xyzzyaaca3.and.abs(xyzzyaabl3-anint(xyzzyaabl3))<xyzz&
&yaaca3)then
ltemp(xyzzyaaac3)=.false.
xyzzyaacc3(1:3,xyzzyaaac3)=0.5d0*sr_lattice(1:3,xyzzyaaab3)
exit
endif
enddo
enddo
if(any(ltemp))call errstop('PWFDET_SETUP','Supercell k_s vector not th&
&e Gamma point or half a supercell reciprocal lattice vector. Unable t&
&o make orbitals into real Bloch functions.')
xyzzyaabm3=xyzzyaacc3(1:3,1)
do xyzzyaaab3=2,xyzzyaacc1
if(abs(xyzzyaabm3(1)-xyzzyaacc3(1,xyzzyaaab3))>xyzzyaaca3.or.abs(xyzzy&
&aabm3(2)-xyzzyaacc3(2,xyzzyaaab3))>xyzzyaaca3.or.abs(xyzzyaabm3(3)-xy&
&zzyaacc3(3,xyzzyaaab3))>xyzzyaaca3)then
call wout('Primitive cell reciprocal lattice vectors: (au)')
call wout('pb1',pb1(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('pb2',pb2(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('pb3',pb3(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('Simulation cell reciprocal lattice vectors: (au)')
call wout('b1 ',b1(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('b2 ',b2(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('b3 ',b3(1:3),fmt='(1x,a,f16.10,1x,f16.10,1x,f16.10)')
call wout('Input k points:')
do xyzzyaaac3=1,xyzzyaacc1
call wout(trim(i2s(xyzzyaaac3)),xyzzyaacl1(:,xyzzyaaac3),fmt='(1x,a,t5&
&,f16.10,1x,f16.10,1x,f16.10)')
enddo
call wout('Input k point set after reduction into 1st Brillouin zone o&
&f simulation cell:')
do xyzzyaaac3=1,xyzzyaacc1
call wout(trim(i2s(xyzzyaaac3)),xyzzyaacc3(:,xyzzyaaac3),fmt='(1x,a,t5&
&,f16.10,1x,f16.10,1x,f16.10)')
enddo
call wout()
call wordwrap('The k point set in pwfn.data does not map to a unique s&
&imulation cell k_s vector.  Are you sure the NPCELL block in input is&
& correct?')
call wout()
call errstop('PWFDET_SETUP','Quitting')
endif
enddo
call wout('MAPPING ONTO UNIQUE K_S VECTOR:')
write(tmpr,'(3f12.8,a)')xyzzyaabm3(1:3),' (Cartesian a.u.)'
call wout(tmpr)
if(any(abs(xyzzyaabm3(1:3))>xyzzyaaca3))then
xyzzyaabj3=one_over_twopi*dot_product(xyzzyaabm3,a1)
xyzzyaabk3=one_over_twopi*dot_product(xyzzyaabm3,a2)
xyzzyaabl3=one_over_twopi*dot_product(xyzzyaabm3,a3)
write(tmpr,'(f12.8,1x,f12.8,1x,f12.8,1x,a)')xyzzyaabj3,xyzzyaabk3,xyzz&
&yaabl3,' (frac supercell recip. lattice vectors)'
call wout(tmpr)
endif
call wout()
deallocate(xyzzyaacc3,ltemp)
if(xyzzyaack3)then
if(xyzzyaaci3)then
call wout('METALLIC GROUND STATE DETECTED')
if(.not.xyzzyaacm1)then
call wout('Number of doubly-occupied bands filled at each k point in G&
&S:')
call wout('k-point    no. of filled bands')
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
write(tmpr,'(i4,4x,i5)')xyzzyaaab3,xyzzyaabf3(xyzzyaaab3,1)
call wout(tmpr)
endif
enddo
else
call wout('Number of singly-occupied bands filled at each k point in G&
&S :')
call wout('Spin    k-point    no. of filled bands')
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
write(tmpr,'(a8,i4,4x,i5)')'up',xyzzyaaab3,xyzzyaabf3(xyzzyaaab3,1)
call wout(tmpr)
endif
enddo
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
write(tmpr,'(a8,i4,4x,i5)')'down',xyzzyaaab3,xyzzyaabf3(xyzzyaaab3,2)
call wout(tmpr)
endif
enddo
endif
else
call wout('INSULATING GROUND STATE DETECTED')
if(.not.xyzzyaacm1)then
call wout('No. of doubly occupied bands at each k point: '//trim(i2s(x&
&yzzyaabf3(1,1))))
else
call wout('No. of occupied spin-up bands at each k point  : '//trim(i2&
&s(xyzzyaabf3(1,1))))
call wout('No. of occupied spin-down bands at each k point: '//trim(i2&
&s(xyzzyaabf3(1,2))))
endif
endif
call wout()
endif
deallocate(xyzzyaabf3)
endif
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
if(periodicity<3.and.abs(dot_product(xyzzyaacl1(1:3,xyzzyaaab3),pa3(1:&
&3)))>xyzzyaaca3)call errstop('PWFDET_SETUP','k vectors inappropriate &
&for periodicity.')
if(periodicity<2.and.abs(dot_product(xyzzyaacl1(1:3,xyzzyaaab3),pa2(1:&
&3)))>xyzzyaaca3)call errstop('PWFDET_SETUP','k vectors inappropriate &
&for periodicity.')
if(periodicity<1.and.abs(dot_product(xyzzyaacl1(1:3,xyzzyaaab3),pa1(1:&
&3)))>xyzzyaaca3)call errstop('PWFDET_SETUP','k vectors inappropriate &
&for periodicity.')
endif
enddo
endif
do xyzzyaaao3=1,nspin
if(xyzzyaaao3==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
do xyzzyaaap3=1,ndet
do xyzzyaaah3=1,wf_nd(xyzzyaaap3,xyzzyaaao3)
xyzzyaaak3=wf_d(1,xyzzyaaah3,xyzzyaaap3,xyzzyaaao3)
xyzzyaaai3=wf_d(2,xyzzyaaah3,xyzzyaaap3,xyzzyaaao3)
xyzzyaaal3=wf_d(3,xyzzyaaah3,xyzzyaaap3,xyzzyaaao3)
xyzzyaaaj3=wf_d(4,xyzzyaaah3,xyzzyaaap3,xyzzyaaao3)
if(xyzzyaaaj3<1.or.xyzzyaaaj3>xyzzyaacc1.or.xyzzyaaai3<1.or.xyzzyaaai3&
&>xyzzyaacc1)call errstop('PWFDET_SETUP','k-point out of range in exci&
&tation specification.')
if(.not.xyzzyaacn1(xyzzyaaaj3).or..not.xyzzyaacn1(xyzzyaaai3))call err&
&stop('PWFDET_SETUP','Specified excitation refers to a k point that is&
& not used because it is the (-) part of a (+/-)k pair.  Use the corre&
&sponding +k in your excitation specification.')
if(xyzzyaaal3<1.or.xyzzyaaal3>xyzzyaace1(xyzzyaaaj3,xyzzyaaad3).or.xyz&
&zyaaak3<1.or.xyzzyaaak3>xyzzyaace1(xyzzyaaai3,xyzzyaaad3))call errsto&
&p('PWFDET_SETUP','band index out of range in excitation specification&
&.')
if(xyzzyaaco1(xyzzyaaaj3))then
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaal3,xyzzyaaaj3,xyzzyaaao3)>=1)call err&
&stop('PWFDET_SETUP','Trying to promote an electron into a band that i&
&s already fully occupied.')
else
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaal3,xyzzyaaaj3,xyzzyaaao3)>=2)call err&
&stop('PWFDET_SETUP','Trying to promote an electron into a band that i&
&s already fully occupied.')
endif
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaak3,xyzzyaaai3,xyzzyaaao3)<=0)call err&
&stop('PWFDET_SETUP','Trying to promote an electron from an empty band&
&.')
xyzzyaacf1(xyzzyaaap3,xyzzyaaal3,xyzzyaaaj3,xyzzyaaao3)=xyzzyaacf1(xyz&
&zyaaap3,xyzzyaaal3,xyzzyaaaj3,xyzzyaaao3)+1
xyzzyaacf1(xyzzyaaap3,xyzzyaaak3,xyzzyaaai3,xyzzyaaao3)=xyzzyaacf1(xyz&
&zyaaap3,xyzzyaaak3,xyzzyaaai3,xyzzyaaao3)-1
enddo
i=0
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
i=i+sum(xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3))
if(xyzzyaaco1(xyzzyaaab3).and.any(xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,x&
&yzzyaaao3)>1))call errstop('PWFDET_SETUP','Problem with excitations (&
&1).')
endif
enddo
if(i/=nuc_nele(xyzzyaaao3))call errstop('PWFDET_SETUP','Problem with e&
&xcitations (2).')
if(any(xyzzyaacf1(xyzzyaaap3,:,:,xyzzyaaao3)<0) .or.any(xyzzyaacf1(xyz&
&zyaaap3,:,:,xyzzyaaao3)>2))call errstop('PWFDET_SETUP','Problem with &
&excitations (3).')
enddo
enddo
deallocate(xyzzyaabd3,xyzzyaacb3,xyzzyaabe3)
allocate(xyzzyaach1(ndet,xyzzyaacd1,xyzzyaacc1,nspin),xyzzyaaci1(ndet,&
&xyzzyaacd1,xyzzyaacc1,nspin),xyzzyaacj1(ndet,xyzzyaacd1,xyzzyaacc1,ns&
&pin),xyzzyaacd3(xyzzyaacd1,xyzzyaacc1),xyzzyaace3(xyzzyaacd1,xyzzyaac&
&c1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','2.5')
xyzzyaach1=0.d0
xyzzyaaci1=1.d0
xyzzyaacj1=0.d0
do xyzzyaaad3=1,xyzzyaaae3
call xyzzyaads1(xyzzyaacd3,xyzzyaace3,xyzzyaaad3)
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
do xyzzyaaaf3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacd3(xyzzyaaaf3,xyzzyaaab3)<xyzzyaace3(xyzzyaaaf3,xyzzyaaab3)&
&)then
xyzzyaach1(:,xyzzyaaaf3,xyzzyaaab3,xyzzyaaad3)=pi_over_two
xyzzyaaci1(:,xyzzyaaaf3,xyzzyaaab3,xyzzyaaad3)=0.d0
xyzzyaacj1(:,xyzzyaaaf3,xyzzyaaab3,xyzzyaaad3)=1.d0
endif
enddo
endif
enddo
enddo
deallocate(xyzzyaacd3,xyzzyaace3)
if(.not.xyzzyaacm1)then
xyzzyaach1(:,:,:,2)=xyzzyaach1(:,:,:,1)
xyzzyaaci1(:,:,:,2)=xyzzyaaci1(:,:,:,1)
xyzzyaacj1(:,:,:,2)=xyzzyaacj1(:,:,:,1)
endif
if(no_orb_phases>0)then
if(am_master)then
call wout('User-defined phases')
call wout('-------------------')
call wout('Number   Det   Spin   Band    k point       Phase (radians)&
&')
endif
do i=1,no_orb_phases
if(am_master)then
write(tmpr,'(1x,i5,1x,i5,5x,i2,1x,i6,2x,i9,7x,f10.6)')i,orb_phase_det(&
&i),orb_phase_spin(i),orb_phase_band(i),orb_phase_kpoint(i),orb_phase(&
&i)
call wout(tmpr)
endif
if(orb_phase_det(i)<1.or.orb_phase_det(i)>ndet)call errstop('PWFDET_SE&
&TUP','Error in determinant in phase spec.')
if(orb_phase_spin(i)<1.or.orb_phase_spin(i)>nspin)call errstop('PWFDET&
&_SETUP','Error in spin in phase spec.')
if(orb_phase_spin(i)==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
if(orb_phase_kpoint(i)<1.or.orb_phase_kpoint(i)>xyzzyaacc1)call errsto&
&p('PWFDET_SETUP','Error in k point in phase spec.')
if(.not.xyzzyaacn1(orb_phase_kpoint(i)))call errstop('PWFDET_SETUP','H&
&ave specified a phase for a k point that is not used because it is pa&
&rt of a (+/-)k pair.')
if(xyzzyaaco1(orb_phase_kpoint(i)).and.am_master)call errwarn('PWFDET_&
&SETUP','a phase has been specified for a k point at the BZ edge. The &
&phase factor should only affect the normalization of the wave functio&
&n in this case.')
if(orb_phase_band(i)<1.or.orb_phase_band(i)>xyzzyaace1(orb_phase_kpoin&
&t(i),xyzzyaaad3))call errstop('PWFDET_SETUP','Error in band in phase &
&spec.')
xyzzyaach1(orb_phase_det(i),orb_phase_band(i),orb_phase_kpoint(i),orb_&
&phase_spin(i))=orb_phase(i)
xyzzyaaci1(orb_phase_det(i),orb_phase_band(i),orb_phase_kpoint(i),orb_&
&phase_spin(i))=cos(orb_phase(i))
xyzzyaacj1(orb_phase_det(i),orb_phase_band(i),orb_phase_kpoint(i),orb_&
&phase_spin(i))=sin(orb_phase(i))
if(xyzzyaacf1(orb_phase_det(i),orb_phase_band(i),orb_phase_kpoint(i),o&
&rb_phase_spin(i))==2.and.am_master)call errwarn('PWFDET_SETUP','a pha&
&se has been specified for a band that is occupied at both k and -k.  &
&This phase will be ignored.')
if(xyzzyaacf1(orb_phase_det(i),orb_phase_band(i),orb_phase_kpoint(i),o&
&rb_phase_spin(i))==0.and.am_master)call errwarn('PWFDET_SETUP','a pha&
&se has been specified for a band that is unoccupied.  This phase will&
& be ignored.')
enddo
if(am_master)call wout()
deallocate(orb_phase_det,orb_phase_band,orb_phase_kpoint,orb_phase_spi&
&n,orb_phase)
endif
if(am_master.and.xyzzyaacl3)then
call wout('Detailed band-occupancy information')
call wout('-----------------------------------')
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
if(xyzzyaaco1(xyzzyaaab3))then
call wout('K POINT '//trim(i2s(xyzzyaaab3))//'  [not part of a (+/-)k &
&pair]')
else
call wout('K POINT '//trim(i2s(xyzzyaaab3))//'  [part of a (+/-)k pair&
&]')
endif
write(tmpr,'("k(au):",3(1x,es23.15))')xyzzyaacl1(1:3,xyzzyaaab3)
call wout(tmpr)
write(tmpr,'("kfrac:",3(1x,es23.15))')matmul(xyzzyaacl1(1:3,xyzzyaaab3&
&),pbinv)
call wout(tmpr)
do xyzzyaaao3=1,nspin
if(nuc_nele(xyzzyaaao3)==0)cycle
if(xyzzyaacm1.and.xyzzyaaao3==2)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
do xyzzyaaap3=1,ndet
if(ndet>1)then
call wout('Spin '//trim(i2s(xyzzyaaao3))//', determinant: '//trim(i2s(&
&xyzzyaaap3))//':')
else
call wout('Spin '//trim(i2s(xyzzyaaao3))//':')
endif
if(any(xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==2))then
char80=' Bands occupied at both k and -k:'
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==2)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam3))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam3))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Highest occupied band energy: ',maxval(xyzzyaack1(:,xyzzyaa&
&ab3,xyzzyaaad3),xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==2))
endif
if(any(xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==1))then
if(xyzzyaaco1(xyzzyaaab3))then
char80=' Bands occupied at k:'
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==1)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam3))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam3))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Highest occupied band energy: ',maxval(xyzzyaack1(:,xyzzyaa&
&ab3,xyzzyaaad3),xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==1))
else
call wout('Single bands formed from the bands at +k and -k:')
call wout('  Band      Phase (radians)')
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==1)then
write(tmpr,'(i6,4x,f16.10)')xyzzyaaam3,xyzzyaach1(xyzzyaaap3,xyzzyaaam&
&3,xyzzyaaab3,xyzzyaaao3)
call wout(tmpr)
endif
enddo
endif
endif
if(any(xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==0))then
if(xyzzyaaco1(xyzzyaaab3))then
char80=' Bands unoccupied at k:'
else
char80=' Bands unoccupied at both k and -k:'
endif
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==0)then
char80_new=trim(char80)//' '//trim(i2s(xyzzyaaam3))
if(len_trim(char80_new)>72)then
call wout(trim(char80),fmt='(a)')
char80=' '//trim(i2s(xyzzyaaam3))
else
char80=char80_new
endif
endif
enddo
call wout(trim(char80),fmt='(a)')
call wout('Lowest unoccupied band energy: ',minval(xyzzyaack1(:,xyzzya&
&aab3,xyzzyaaad3),xyzzyaacf1(xyzzyaaap3,:,xyzzyaaab3,xyzzyaaao3)==0))
endif
enddo
enddo
else
write(tmpr,'(a,3(f10.6,a))')'K POINT '//trim(i2s(xyzzyaaab3))//'  [(',&
&xyzzyaacl1(1,xyzzyaaab3),',',xyzzyaacl1(2,xyzzyaaab3),',',xyzzyaacl1(&
&3,xyzzyaaab3),')]:  NOT USED.'
call wout(tmpr)
endif
call wout()
enddo
endif
deallocate(xyzzyaacm3)
endif
if(.not.xyzzyaacm1)then
xyzzyaaaz3=xyzzyaacd1*xyzzyaacc1
else
xyzzyaaaz3=sum(xyzzyaace1)
endif
xyzzyaaay3=2*xyzzyaaaz3*ndet
allocate(xyzzyaada1(nemax,nspin,ndet),xyzzyaabg3(ndet,nspin),xyzzyaadc&
&1(xyzzyaaaz3),xyzzyaadd1(xyzzyaaaz3),xyzzyaade1(xyzzyaaaz3),xyzzyaadf&
&1(xyzzyaaay3),xyzzyaadg1(xyzzyaaay3),xyzzyaadb1(xyzzyaaay3),stat=xyzz&
&yaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','pw_orbmap,...')
xyzzyaada1=0
xyzzyaabg3=0
xyzzyaadc1=0
xyzzyaadd1=0
xyzzyaade1=0
xyzzyaadf1=czero
xyzzyaadg1=0
xyzzyaadb1=0
xyzzyaaaw3=0
xyzzyaaax3=0
do xyzzyaaao3=1,nspin
if(xyzzyaaao3==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
if(.not.xyzzyaacm1)xyzzyaaax3=0
do xyzzyaaab3=1,xyzzyaacc1
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaad3)
if((xyzzyaacm1.and.all(xyzzyaacf1(1:ndet,xyzzyaaam3,xyzzyaaab3,xyzzyaa&
&ao3)==0)).or.(.not.xyzzyaacm1.and.all(xyzzyaacf1(1:ndet,xyzzyaaam3,xy&
&zzyaaab3,1:nspin)==0)))cycle
xyzzyaaax3=xyzzyaaax3+1
xyzzyaadc1(xyzzyaaax3)=xyzzyaaam3
xyzzyaadd1(xyzzyaaax3)=xyzzyaaab3
xyzzyaade1(xyzzyaaax3)=xyzzyaaad3
if(xyzzyaaca1)then
do xyzzyaaap3=1,ndet
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==1)then
xyzzyaabg3(xyzzyaaap3,xyzzyaaao3)=xyzzyaabg3(xyzzyaaap3,xyzzyaaao3)+1
xyzzyaada1(xyzzyaabg3(xyzzyaaap3,xyzzyaaao3),xyzzyaaao3,xyzzyaaap3)=xy&
&zzyaaax3
endif
enddo
else
do xyzzyaaap3=1,ndet
if(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==2)then
xyzzyaacj3=.false.
do xyzzyaaav3=1,xyzzyaaaw3
if(xyzzyaadb1(xyzzyaaav3)/=xyzzyaaax3)cycle
if(xyzzyaadg1(xyzzyaaav3)==xyzzyaadh1)then
xyzzyaacj3=.true.
exit
endif
enddo
if(.not.xyzzyaacj3)then
xyzzyaaaw3=xyzzyaaaw3+1
xyzzyaaav3=xyzzyaaaw3
xyzzyaadf1(xyzzyaaav3)=c_one
xyzzyaadg1(xyzzyaaav3)=xyzzyaadh1
xyzzyaadb1(xyzzyaaav3)=xyzzyaaax3
endif
xyzzyaabg3(xyzzyaaap3,xyzzyaaao3)=xyzzyaabg3(xyzzyaaap3,xyzzyaaao3)+1
xyzzyaada1(xyzzyaabg3(xyzzyaaap3,xyzzyaaao3),xyzzyaaao3,xyzzyaaap3)=xy&
&zzyaaav3
xyzzyaacj3=.false.
do xyzzyaaav3=1,xyzzyaaaw3
if(xyzzyaadb1(xyzzyaaav3)/=xyzzyaaax3)cycle
if(xyzzyaadg1(xyzzyaaav3)==xyzzyaadi1)then
xyzzyaacj3=.true.
exit
endif
enddo
if(.not.xyzzyaacj3)then
xyzzyaaaw3=xyzzyaaaw3+1
xyzzyaaav3=xyzzyaaaw3
xyzzyaadf1(xyzzyaaav3)=zi
xyzzyaadg1(xyzzyaaav3)=xyzzyaadi1
xyzzyaadb1(xyzzyaaav3)=xyzzyaaax3
endif
xyzzyaabg3(xyzzyaaap3,xyzzyaaao3)=xyzzyaabg3(xyzzyaaap3,xyzzyaaao3)+1
xyzzyaada1(xyzzyaabg3(xyzzyaaap3,xyzzyaaao3),xyzzyaaao3,xyzzyaaap3)=xy&
&zzyaaav3
elseif(xyzzyaacf1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)==1)then
xyzzyaaby3=xyzzyaaci1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)
xyzzyaabz3=xyzzyaacj1(xyzzyaaap3,xyzzyaaam3,xyzzyaaab3,xyzzyaaao3)
xyzzyaacj3=.false.
do xyzzyaaav3=1,xyzzyaaaw3
if(xyzzyaadb1(xyzzyaaav3)/=xyzzyaaax3)cycle
if(dble(xyzzyaadf1(xyzzyaaav3))==xyzzyaaby3.and.aimag(xyzzyaadf1(xyzzy&
&aaav3))==xyzzyaabz3)then
xyzzyaacj3=.true.
exit
endif
enddo
if(.not.xyzzyaacj3)then
xyzzyaaaw3=xyzzyaaaw3+1
xyzzyaaav3=xyzzyaaaw3
xyzzyaacf3=cmplx(xyzzyaaby3,xyzzyaabz3,dp)
if(abs(xyzzyaaby3)<1.d3*epsilon(1.d0))xyzzyaacf3=cmplx(0.d0,xyzzyaabz3&
&,dp)
if(abs(xyzzyaabz3)<1.d3*epsilon(1.d0))xyzzyaacf3=cmplx(xyzzyaaby3,0.d0&
&,dp)
xyzzyaadf1(xyzzyaaav3)=xyzzyaacf3
if(xyzzyaadf1(xyzzyaaav3)==c_one)then
xyzzyaadg1(xyzzyaaav3)=xyzzyaadh1
elseif(xyzzyaadf1(xyzzyaaav3)==zi)then
xyzzyaadg1(xyzzyaaav3)=xyzzyaadi1
else
xyzzyaadg1(xyzzyaaav3)=xyzzyaadj1
endif
xyzzyaadb1(xyzzyaaav3)=xyzzyaaax3
endif
xyzzyaabg3(xyzzyaaap3,xyzzyaaao3)=xyzzyaabg3(xyzzyaaap3,xyzzyaaao3)+1
xyzzyaada1(xyzzyaabg3(xyzzyaaap3,xyzzyaaao3),xyzzyaaao3,xyzzyaaap3)=xy&
&zzyaaav3
endif
enddo
endif
enddo
enddo
enddo
deallocate(xyzzyaabg3)
xyzzyaacz1=xyzzyaaax3
xyzzyaacy1=xyzzyaaaw3
if(xyzzyaaca1)then
xyzzyaacx1=xyzzyaacz1
else
xyzzyaacx1=xyzzyaacy1
endif
allocate(xyzzyaado1(xyzzyaadm1,xyzzyaacx1),xyzzyaadp1(xyzzyaadn1,xyzzy&
&aacx1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'BWFDET_SETUP','pw_orbdesc_*')
xyzzyaabb3=0
do xyzzyaaan3=1,xyzzyaacx1
if(xyzzyaaca1)then
xyzzyaaax3=xyzzyaaan3
else
xyzzyaaax3=xyzzyaadb1(xyzzyaaan3)
endif
xyzzyaaab3=xyzzyaadd1(xyzzyaaax3)
xyzzyaaam3=xyzzyaadc1(xyzzyaaax3)
xyzzyaaad3=xyzzyaade1(xyzzyaaax3)
xyzzyaado1(1,xyzzyaaan3)=xyzzyaaab3
xyzzyaacj3=.false.
do xyzzyaabc3=1,xyzzyaaan3-1
if(xyzzyaack1(xyzzyaaam3,xyzzyaaab3,xyzzyaaad3)==xyzzyaadp1(1,xyzzyaab&
&c3))then
xyzzyaacj3=.true.
exit
endif
enddo
if(.not.xyzzyaacj3)then
xyzzyaabb3=xyzzyaabb3+1
xyzzyaaba3=xyzzyaabb3
xyzzyaado1(2,xyzzyaaan3)=xyzzyaaba3
xyzzyaadp1(1,xyzzyaaan3)=xyzzyaack1(xyzzyaaam3,xyzzyaaab3,xyzzyaaad3)
else
xyzzyaado1(2,xyzzyaaan3)=xyzzyaado1(2,xyzzyaabc3)
xyzzyaadp1(1,xyzzyaaan3)=xyzzyaadp1(1,xyzzyaabc3)
endif
enddo
deallocate(xyzzyaack1)
allocate(xyzzyaach3(xyzzyaacy1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','ztrf')
xyzzyaach3(1:xyzzyaacy1)=xyzzyaadf1(1:xyzzyaacy1)
deallocate(xyzzyaadf1)
allocate(xyzzyaadf1(xyzzyaacy1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','rorb_angle')
xyzzyaadf1(1:xyzzyaacy1)=xyzzyaach3(1:xyzzyaacy1)
deallocate(xyzzyaach3)
allocate(xyzzyaabh3(xyzzyaacy1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','itrf')
xyzzyaabh3(1:xyzzyaacy1)=xyzzyaadg1(1:xyzzyaacy1)
deallocate(xyzzyaadg1)
allocate(xyzzyaadg1(xyzzyaacy1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','rorb_angle_type')
xyzzyaadg1(1:xyzzyaacy1)=xyzzyaabh3(1:xyzzyaacy1)
deallocate(xyzzyaabh3)
allocate(xyzzyaabh3(xyzzyaacy1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','itrf')
xyzzyaabh3(1:xyzzyaacy1)=xyzzyaadb1(1:xyzzyaacy1)
deallocate(xyzzyaadb1)
allocate(xyzzyaadb1(xyzzyaacy1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','Ridx2Zidx')
xyzzyaadb1(1:xyzzyaacy1)=xyzzyaabh3(1:xyzzyaacy1)
deallocate(xyzzyaabh3)
allocate(xyzzyaabh3(xyzzyaacz1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','itrf')
xyzzyaabh3(1:xyzzyaacz1)=xyzzyaadc1(1:xyzzyaacz1)
deallocate(xyzzyaadc1)
allocate(xyzzyaadc1(xyzzyaacz1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','corb_band_idx')
xyzzyaadc1(1:xyzzyaacz1)=xyzzyaabh3(1:xyzzyaacz1)
deallocate(xyzzyaabh3)
allocate(xyzzyaabh3(xyzzyaacz1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','itrf')
xyzzyaabh3(1:xyzzyaacz1)=xyzzyaadd1(1:xyzzyaacz1)
deallocate(xyzzyaadd1)
allocate(xyzzyaadd1(xyzzyaacz1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','corb_kvec_idx')
xyzzyaadd1(1:xyzzyaacz1)=xyzzyaabh3(1:xyzzyaacz1)
deallocate(xyzzyaabh3)
allocate(xyzzyaabh3(xyzzyaacz1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','itrf')
xyzzyaabh3(1:xyzzyaacz1)=xyzzyaade1(1:xyzzyaacz1)
deallocate(xyzzyaade1)
allocate(xyzzyaade1(xyzzyaacz1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','corb_spin_idx')
xyzzyaade1(1:xyzzyaacz1)=xyzzyaabh3(1:xyzzyaacz1)
deallocate(xyzzyaabh3)
xyzzyaacg1=0.d0
do i=2,nemax
xyzzyaacg1=xyzzyaacg1+log(dble(i))
enddo
xyzzyaacg1=orb_norm*exp(-xyzzyaacg1/dble(2*nemax))
if(am_master)call xyzzyaadr1(kwarn)
if(use_backflow)then
if(xyzzyaabz1)then
allocate(xyzzyaaal1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaam1(xyzzy&
&aaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaan1(xyzzyaaaf1,xyzzyaacd1,xyzzyaa&
&cc1),xyzzyaaar1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaat1(xyzzyaaa&
&f1,xyzzyaacd1,xyzzyaacc1),xyzzyaaau1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1&
&),xyzzyaaav1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaaw1(xyzzyaaaf1,&
&xyzzyaacd1,xyzzyaacc1),xyzzyaaax1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),x&
&yzzyaaay1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','1a')
if(xyzzyaacm1)then
allocate(xyzzyaaao1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaap1(xyzzy&
&aaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaaq1(xyzzyaaaf1,xyzzyaacd1,xyzzyaa&
&cc1),xyzzyaaas1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaaz1(xyzzyaaa&
&f1,xyzzyaacd1,xyzzyaacc1),xyzzyaaba1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1&
&),xyzzyaabb1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabc1(xyzzyaaaf1,&
&xyzzyaacd1,xyzzyaacc1),xyzzyaabd1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),x&
&yzzyaabe1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','1b')
endif
else
allocate(xyzzyaabf1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabg1(xyzzy&
&aaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabh1(xyzzyaaaf1,xyzzyaacd1,xyzzyaa&
&cc1),xyzzyaabl1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabn1(xyzzyaaa&
&f1,xyzzyaacd1,xyzzyaacc1),xyzzyaabo1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1&
&),xyzzyaabp1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabq1(xyzzyaaaf1,&
&xyzzyaacd1,xyzzyaacc1),xyzzyaabr1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),x&
&yzzyaabs1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','1c')
if(xyzzyaacm1)then
allocate(xyzzyaabi1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabj1(xyzzy&
&aaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabk1(xyzzyaaaf1,xyzzyaacd1,xyzzyaa&
&cc1),xyzzyaabm1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabt1(xyzzyaaa&
&f1,xyzzyaacd1,xyzzyaacc1),xyzzyaabu1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1&
&),xyzzyaabv1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabw1(xyzzyaaaf1,&
&xyzzyaacd1,xyzzyaacc1),xyzzyaabx1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),x&
&yzzyaaby1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','1d')
endif
endif
else
if(xyzzyaabz1)then
allocate(xyzzyaaal1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaam1(xyzzy&
&aaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaan1(xyzzyaaaf1,xyzzyaacd1,xyzzyaa&
&cc1),xyzzyaaar1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','1a2')
if(xyzzyaacm1)then
allocate(xyzzyaaao1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaap1(xyzzy&
&aaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaaaq1(xyzzyaaaf1,xyzzyaacd1,xyzzyaa&
&cc1),xyzzyaaas1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','1b2')
endif
else
allocate(xyzzyaabf1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabg1(xyzzy&
&aaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabh1(xyzzyaaaf1,xyzzyaacd1,xyzzyaa&
&cc1),xyzzyaabl1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','1c2')
if(xyzzyaacm1)then
allocate(xyzzyaabi1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabj1(xyzzy&
&aaaf1,xyzzyaacd1,xyzzyaacc1),xyzzyaabk1(xyzzyaaaf1,xyzzyaacd1,xyzzyaa&
&cc1),xyzzyaabm1(xyzzyaaaf1,xyzzyaacd1,xyzzyaacc1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','1d2')
endif
endif
endif
allocate(xyzzyaaab1(3,xyzzyaaaf1),xyzzyaaad1(xyzzyaacc1),xyzzyaaae1(xy&
&zzyaaaf1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','gmap,...')
xyzzyaabp3(1,1)=pb11
xyzzyaabp3(1,2)=pb12
xyzzyaabp3(1,3)=pb13
xyzzyaabp3(2,1)=pb12
xyzzyaabp3(2,2)=pb22
xyzzyaabp3(2,3)=pb23
xyzzyaabp3(3,1)=pb13
xyzzyaabp3(3,2)=pb23
xyzzyaabp3(3,3)=pb33
xyzzyaabq3=inverse3(xyzzyaabp3)
xyzzyaaaa1=1
do xyzzyaaar3=1,xyzzyaaaf1
xyzzyaabr3=xyzzyaaag1(1,xyzzyaaar3)*pb1(1)+xyzzyaaag1(2,xyzzyaaar3)*pb&
&1(2)+xyzzyaaag1(3,xyzzyaaar3)*pb1(3)
xyzzyaabs3=xyzzyaaag1(1,xyzzyaaar3)*pb2(1)+xyzzyaaag1(2,xyzzyaaar3)*pb&
&2(2)+xyzzyaaag1(3,xyzzyaaar3)*pb2(3)
xyzzyaabt3=xyzzyaaag1(1,xyzzyaaar3)*pb3(1)+xyzzyaaag1(2,xyzzyaaar3)*pb&
&3(2)+xyzzyaaag1(3,xyzzyaaar3)*pb3(3)
xyzzyaaas3=nint(xyzzyaabq3(1,1)*xyzzyaabr3+xyzzyaabq3(1,2)*xyzzyaabs3+&
&xyzzyaabq3(1,3)*xyzzyaabt3)
xyzzyaaat3=nint(xyzzyaabq3(2,1)*xyzzyaabr3+xyzzyaabq3(2,2)*xyzzyaabs3+&
&xyzzyaabq3(2,3)*xyzzyaabt3)
xyzzyaaau3=nint(xyzzyaabq3(3,1)*xyzzyaabr3+xyzzyaabq3(3,2)*xyzzyaabs3+&
&xyzzyaabq3(3,3)*xyzzyaabt3)
xyzzyaaab1(1,xyzzyaaar3)=xyzzyaaas3
xyzzyaaab1(2,xyzzyaaar3)=xyzzyaaat3
xyzzyaaab1(3,xyzzyaaar3)=xyzzyaaau3
xyzzyaaaa1=max(abs(xyzzyaaas3),abs(xyzzyaaat3),abs(xyzzyaaau3),xyzzyaa&
&aa1)
enddo
allocate(xyzzyaaac1(3,-xyzzyaaaa1:xyzzyaaaa1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','etemp array')
xyzzyaaac1(1:3,0)=c_one
do xyzzyaaao3=1,xyzzyaaae3
if(xyzzyaaao3==2.and.xyzzyaacm1)then
xyzzyaaad3=2
else
xyzzyaaad3=1
endif
do xyzzyaaab3=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaab3))then
do xyzzyaaam3=1,xyzzyaace1(xyzzyaaab3,xyzzyaaao3)
xyzzyaabj3=xyzzyaacl1(1,xyzzyaaab3)
xyzzyaabk3=xyzzyaacl1(2,xyzzyaaab3)
xyzzyaabl3=xyzzyaacl1(3,xyzzyaaab3)
do xyzzyaaar3=1,xyzzyaaaf1
xyzzyaabu3=xyzzyaaag1(1,xyzzyaaar3)+xyzzyaabj3
xyzzyaabw3=xyzzyaaag1(2,xyzzyaaar3)+xyzzyaabk3
xyzzyaabv3=xyzzyaaag1(3,xyzzyaaar3)+xyzzyaabl3
xyzzyaabx3=xyzzyaabu3*xyzzyaabu3+xyzzyaabw3*xyzzyaabw3+xyzzyaabv3*xyzz&
&yaabv3
if(xyzzyaabz1)then
if(xyzzyaaad3==1)then
xyzzyaabo3=xyzzyaaah1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)
xyzzyaaal1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabu3*xyzzyaabo3
xyzzyaaam1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabw3*xyzzyaabo3
xyzzyaaan1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabv3*xyzzyaabo3
xyzzyaaar1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabx3*xyzzyaabo3
if(use_backflow)then
xyzzyaaat1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabu3*xy&
&zzyaabo3
xyzzyaaau1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabw3*xyzzyaabw3*xy&
&zzyaabo3
xyzzyaaav1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabv3*xyzzyaabv3*xy&
&zzyaabo3
xyzzyaaaw1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabw3*xy&
&zzyaabo3
xyzzyaaax1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabv3*xy&
&zzyaabo3
xyzzyaaay1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabw3*xyzzyaabv3*xy&
&zzyaabo3
endif
else
xyzzyaabo3=xyzzyaaai1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)
xyzzyaaao1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabu3*xyzzyaabo3
xyzzyaaap1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabw3*xyzzyaabo3
xyzzyaaaq1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabv3*xyzzyaabo3
xyzzyaaas1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabx3*xyzzyaabo3
if(use_backflow)then
xyzzyaaaz1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabu3*xy&
&zzyaabo3
xyzzyaaba1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabw3*xyzzyaabw3*xy&
&zzyaabo3
xyzzyaabb1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabv3*xyzzyaabv3*xy&
&zzyaabo3
xyzzyaabc1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabw3*xy&
&zzyaabo3
xyzzyaabd1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabv3*xy&
&zzyaabo3
xyzzyaabe1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabw3*xyzzyaabv3*xy&
&zzyaabo3
endif
endif
else
if(xyzzyaaad3==1)then
xyzzyaacf3=xyzzyaaaj1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)
xyzzyaacg3=cmplx(-aimag(xyzzyaacf3),dble(xyzzyaacf3),dp)
xyzzyaabf1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabu3*xyzzyaacg3
xyzzyaabg1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabw3*xyzzyaacg3
xyzzyaabh1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabv3*xyzzyaacg3
xyzzyaabl1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabx3*xyzzyaacf3
if(use_backflow)then
xyzzyaabn1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabu3*xy&
&zzyaacf3
xyzzyaabo1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabw3*xyzzyaabw3*xy&
&zzyaacf3
xyzzyaabp1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabv3*xyzzyaabv3*xy&
&zzyaacf3
xyzzyaabq1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabw3*xy&
&zzyaacf3
xyzzyaabr1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabv3*xy&
&zzyaacf3
xyzzyaabs1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabw3*xyzzyaabv3*xy&
&zzyaacf3
endif
else
xyzzyaacf3=xyzzyaaak1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)
xyzzyaacg3=cmplx(-aimag(xyzzyaacf3),dble(xyzzyaacf3),dp)
xyzzyaabi1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabu3*xyzzyaacg3
xyzzyaabj1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabw3*xyzzyaacg3
xyzzyaabk1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=xyzzyaabv3*xyzzyaacg3
xyzzyaabm1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabx3*xyzzyaacf3
if(use_backflow)then
xyzzyaabt1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabu3*xy&
&zzyaacf3
xyzzyaabu1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabw3*xyzzyaabw3*xy&
&zzyaacf3
xyzzyaabv1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabv3*xyzzyaabv3*xy&
&zzyaacf3
xyzzyaabw1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabw3*xy&
&zzyaacf3
xyzzyaabx1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabu3*xyzzyaabv3*xy&
&zzyaacf3
xyzzyaaby1(xyzzyaaar3,xyzzyaaam3,xyzzyaaab3)=-xyzzyaabw3*xyzzyaabv3*xy&
&zzyaacf3
endif
endif
endif
enddo
enddo
endif
enddo
enddo
allocate(xyzzyaacp1(xyzzyaacz1),xyzzyaacq1(3,xyzzyaacz1),xyzzyaacr1(xy&
&zzyaacz1),xyzzyaacs1(6,xyzzyaacz1),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','val_c')
allocate(xyzzyaadk1(xyzzyaacx1,2),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','pw_orbmask')
call get_numerical_orbmask(xyzzyaacx1,nemax,nspin,ndet,xyzzyaada1,xyzz&
&yaadk1)
allocate(xyzzyaadl1(xyzzyaacz1,nspin),stat=xyzzyaaaa3)
call check_alloc(xyzzyaaaa3,'PWFDET_SETUP','pw_corbmask')
xyzzyaadl1=.false.
do xyzzyaaao3=1,nspin
do xyzzyaaan3=1,xyzzyaacx1
if(.not.xyzzyaadk1(xyzzyaaan3,xyzzyaaao3))cycle
if(xyzzyaaca1)then
xyzzyaadl1(xyzzyaaan3,xyzzyaaao3)=.true.
else
xyzzyaadl1(xyzzyaadb1(xyzzyaaan3),xyzzyaaao3)=.true.
endif
enddo
enddo
if(use_gpcc)call xyzzyaadq1
deallocate(xyzzyaadk1)
end subroutine pwfdet_setup
subroutine xyzzyaadq1
use slaarnabi, only : use_gpcc,initialize_gpcc,setup_gpcc,spherical_av&
&_cmplx,nradgrid,radgrid,nsphgrid,sphgrid,ae_index,naeions_prim
implicit none
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4,xyzzyaaae4,xyzzyaa&
&af4,xyzzyaaag4,ik,xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4
real(dp) xyzzyaaak4(3)
real(dp),allocatable :: xyzzyaaal4(:,:),xyzzyaaam4(:,:,:),xyzzyaaan4(:&
&,:),xyzzyaaao4(:,:)
complex(dp),allocatable :: xyzzyaaap4(:,:)
logical xyzzyaaaq4
if(.not.complex_wf)then
do ik=1,xyzzyaacc1
if(xyzzyaacn1(ik).and..not.xyzzyaaco1(ik).and.any(xyzzyaacf1(:,:,ik,:)&
&==1))call errstop('PWFDET_SETUP_GPCC','Cannot cusp-correct orbitals a&
&way from the BZ edge in which only one of k and -k is occupied.  Ask &
&Neil if you need to do this, or use a complex wave function.')
enddo
endif
if(xyzzyaacm1.or.nuc_nele(2)/=nuc_nele(1).or.any(xyzzyaacf1(:,:,:,2)/=&
&xyzzyaacf1(:,:,:,1)))then
xyzzyaaaa4=1
else
xyzzyaaaa4=0
endif
call initialize_gpcc(xyzzyaaaa4,xyzzyaacz1,xyzzyaadl1,.false.)
allocate(xyzzyaaal4(xyzzyaacx1,real1_complex2),xyzzyaaam4(3,xyzzyaacx1&
&,real1_complex2),xyzzyaaan4(xyzzyaacx1,real1_complex2),stat=xyzzyaaab&
&4)
call check_alloc(xyzzyaaab4,'PWFDET_SETUP_GPCC','valtemp')
allocate(xyzzyaaap4(nsphgrid,xyzzyaacz1),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,'PWFDET_SETUP_GPCC','orb_sphgrid')
xyzzyaaap4=czero
use_gpcc=.false.
do xyzzyaaac4=1,no_ssingles(xyzzyaaaa4)
do xyzzyaaai4=1,naeions_prim
xyzzyaaah4=ae_index(xyzzyaaai4)
do xyzzyaaad4=1,nradgrid
if(xyzzyaaad4>1)then
xyzzyaaaj4=nsphgrid
else
xyzzyaaaj4=1
endif
do xyzzyaaae4=1,xyzzyaaaj4
xyzzyaaak4=rion(1:3,xyzzyaaah4)+radgrid(xyzzyaaad4,xyzzyaaai4)*sphgrid&
&(1:3,xyzzyaaae4)
call pw_orb_eval(xyzzyaaak4,xyzzyaaac4,xyzzyaacx1,xyzzyaadk1(1,xyzzyaa&
&ac4),.true.,.false.,xyzzyaaal4,xyzzyaaam4,xyzzyaaan4)
xyzzyaaaq4=.false.
do xyzzyaaaf4=1,xyzzyaacx1
if(.not.xyzzyaadk1(xyzzyaaaf4,xyzzyaaac4))cycle
if(xyzzyaaca1)then
xyzzyaaap4(xyzzyaaae4,xyzzyaaaf4)=cmplx(xyzzyaaal4(xyzzyaaaf4,1),xyzzy&
&aaal4(xyzzyaaaf4,2),dp)
else
if(xyzzyaaaq4)then
xyzzyaaaq4=.false.
cycle
endif
xyzzyaaag4=xyzzyaadc1(xyzzyaadb1(xyzzyaaaf4))
ik=xyzzyaadd1(xyzzyaadb1(xyzzyaaaf4))
if(any(xyzzyaacf1(:,xyzzyaaag4,ik,xyzzyaaac4)==2))then
xyzzyaaap4(xyzzyaaae4,xyzzyaadb1(xyzzyaaaf4))=cmplx(xyzzyaaal4(xyzzyaa&
&af4,1),xyzzyaaal4(xyzzyaaaf4+1,1),dp)
xyzzyaaaq4=.true.
elseif(any(xyzzyaacf1(:,xyzzyaaag4,ik,xyzzyaaac4)==1))then
xyzzyaaap4(xyzzyaaae4,xyzzyaadb1(xyzzyaaaf4))=cmplx(xyzzyaaal4(xyzzyaa&
&af4,1)*dble(xyzzyaadf1(xyzzyaaaf4)),xyzzyaaal4(xyzzyaaaf4,1)*aimag(xy&
&zzyaadf1(xyzzyaaaf4)),dp)
endif
endif
enddo
enddo
call spherical_av_cmplx(xyzzyaaac4,xyzzyaaai4,xyzzyaaad4,xyzzyaaap4)
enddo
enddo
enddo
deallocate(xyzzyaaal4,xyzzyaaam4,xyzzyaaan4,xyzzyaaap4)
use_gpcc=.true.
allocate(xyzzyaaao4(3,xyzzyaacz1),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,'PWFDET_SETUP_GPCC','kvec_gpcc')
xyzzyaaao4=0.d0
do xyzzyaaaf4=1,xyzzyaacx1
if(xyzzyaaca1)then
xyzzyaaao4(1:3,xyzzyaaaf4)=xyzzyaacl1(1:3,xyzzyaadd1(xyzzyaaaf4))
else
xyzzyaaao4(1:3,xyzzyaadb1(xyzzyaaaf4))=xyzzyaacl1(1:3,xyzzyaadd1(xyzzy&
&aadb1(xyzzyaaaf4)))
endif
enddo
call setup_gpcc(xyzzyaaao4)
deallocate(xyzzyaaao4)
allocate(xyzzyaact1(xyzzyaacz1),xyzzyaacu1(3,xyzzyaacz1),xyzzyaacv1(xy&
&zzyaacz1),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,'PWFDET_SETUP_GPCC','cusp_val_c')
if(use_backflow)then
allocate(xyzzyaacw1(6,xyzzyaacz1),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,'PWFDET_SETUP_GPCC','cusp_sderivs_c')
endif
end subroutine xyzzyaadq1
subroutine xyzzyaadr1(kwarn)
use slaarnabp,only : modified_mdet
use slaarnabt, only : ddot
implicit none
logical,intent(in) :: kwarn
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5,xyzzyaaad5,xyzzyaaae5,xyzzyaa&
&af5,xyzzyaaag5,xyzzyaaah5,xyzzyaaai5
real(dp) xyzzyaaaj5,xyzzyaaak5,xyzzyaaal5(3),xyzzyaaam5,xyzzyaaan5
real(dp),parameter :: xyzzyaaao5=1.d-8,xyzzyaaap5=xyzzyaaao5*xyzzyaaao&
&5
real(dp),allocatable :: xyzzyaaaq5(:)
logical xyzzyaaar5,xyzzyaaas5,xyzzyaaat5
xyzzyaaar5=.false.
if(xyzzyaacb1/=0.d0)then
call wout('Kinetic energy from pwfn.data                 : '//trim(d2s&
&(xyzzyaacb1)))
xyzzyaaas5=.true.
else
call wout('DFT kinetic energy not given in pwfn.data. No check made.')
xyzzyaaas5=.false.
endif
if(ndet>1.or.modified_mdet)xyzzyaaas5=.false.
if(xyzzyaacm1)then
xyzzyaaai5=2
else
xyzzyaaai5=1
endif
xyzzyaaat5=.true.
spin_loop: do xyzzyaaag5=1,xyzzyaaai5
do xyzzyaaac5=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaac5))then
allocate(xyzzyaaaq5(xyzzyaace1(xyzzyaaac5,xyzzyaaag5)),stat=xyzzyaaah5&
&)
call check_alloc(xyzzyaaah5,'PWFDET_CF_DFT_KE','orbnorm2')
do xyzzyaaad5=1,xyzzyaace1(xyzzyaaac5,xyzzyaaag5)
if(xyzzyaabz1)then
if(xyzzyaaag5==1)then
xyzzyaaaq5(xyzzyaaad5)=ddot(xyzzyaaaf1,xyzzyaaah1(1,xyzzyaaad5,xyzzyaa&
&ac5),1,xyzzyaaah1(1,xyzzyaaad5,xyzzyaaac5),1)
else
xyzzyaaaq5(xyzzyaaad5)=ddot(xyzzyaaaf1,xyzzyaaai1(1,xyzzyaaad5,xyzzyaa&
&ac5),1,xyzzyaaai1(1,xyzzyaaad5,xyzzyaaac5),1)
endif
else
if(xyzzyaaag5==1)then
xyzzyaaaq5(xyzzyaaad5)=abs(dot_product(xyzzyaaaj1(1:xyzzyaaaf1,xyzzyaa&
&ad5,xyzzyaaac5),xyzzyaaaj1(1:xyzzyaaaf1,xyzzyaaad5,xyzzyaaac5)))
else
xyzzyaaaq5(xyzzyaaad5)=abs(dot_product(xyzzyaaak1(1:xyzzyaaaf1,xyzzyaa&
&ad5,xyzzyaaac5),xyzzyaaak1(1:xyzzyaaaf1,xyzzyaaad5,xyzzyaaac5)))
endif
endif
enddo
do xyzzyaaad5=1,xyzzyaace1(xyzzyaaac5,xyzzyaaag5)-1
do xyzzyaaae5=xyzzyaaad5+1,xyzzyaace1(xyzzyaaac5,xyzzyaaag5)
if(xyzzyaabz1)then
if(xyzzyaaag5==1)then
xyzzyaaan5=ddot(xyzzyaaaf1,xyzzyaaah1(1,xyzzyaaad5,xyzzyaaac5),1,xyzzy&
&aaah1(1,xyzzyaaae5,xyzzyaaac5),1)**2
else
xyzzyaaan5=ddot(xyzzyaaaf1,xyzzyaaai1(1,xyzzyaaad5,xyzzyaaac5),1,xyzzy&
&aaai1(1,xyzzyaaae5,xyzzyaaac5),1)**2
endif
else
if(xyzzyaaag5==1)then
xyzzyaaan5=abs(dot_product(xyzzyaaaj1(1:xyzzyaaaf1,xyzzyaaad5,xyzzyaaa&
&c5),xyzzyaaaj1(1:xyzzyaaaf1,xyzzyaaae5,xyzzyaaac5)))**2
else
xyzzyaaan5=abs(dot_product(xyzzyaaak1(1:xyzzyaaaf1,xyzzyaaad5,xyzzyaaa&
&c5),xyzzyaaak1(1:xyzzyaaaf1,xyzzyaaae5,xyzzyaaac5)))**2
endif
endif
if(xyzzyaaan5>xyzzyaaap5*xyzzyaaaq5(xyzzyaaad5)*xyzzyaaaq5(xyzzyaaae5)&
&)then
xyzzyaaat5=.false.
exit spin_loop
endif
enddo
enddo
deallocate(xyzzyaaaq5)
endif
enddo
enddo spin_loop
if(.not.xyzzyaaat5)then
call wout()
call errwarn('PWFDET_CF_DFT_KE','Your plane-wave orbitals are nonortho&
&gonal.')
xyzzyaaas5=.false.
endif
do xyzzyaaaa5=1,ndet
xyzzyaaaj5=0.d0
do xyzzyaaab5=1,nspin
if(xyzzyaaab5==2.and.xyzzyaacm1)then
xyzzyaaag5=2
else
xyzzyaaag5=1
endif
do xyzzyaaac5=1,xyzzyaacc1
do xyzzyaaad5=1,xyzzyaace1(xyzzyaaac5,xyzzyaaab5)
if(xyzzyaacf1(xyzzyaaaa5,xyzzyaaad5,xyzzyaaac5,xyzzyaaab5)>0)then
xyzzyaaak5=0.d0
if(xyzzyaabz1)then
if(xyzzyaaag5==1)then
do xyzzyaaaf5=1,xyzzyaaaf1
xyzzyaaal5=xyzzyaacl1(:,xyzzyaaac5)+xyzzyaaag1(:,xyzzyaaaf5)
xyzzyaaam5=xyzzyaaal5(1)*xyzzyaaal5(1)+xyzzyaaal5(2)*xyzzyaaal5(2)+xyz&
&zyaaal5(3)*xyzzyaaal5(3)
xyzzyaaak5=xyzzyaaak5+xyzzyaaah1(xyzzyaaaf5,xyzzyaaad5,xyzzyaaac5)**2*&
&xyzzyaaam5
enddo
else
do xyzzyaaaf5=1,xyzzyaaaf1
xyzzyaaal5=xyzzyaacl1(:,xyzzyaaac5)+xyzzyaaag1(:,xyzzyaaaf5)
xyzzyaaam5=xyzzyaaal5(1)*xyzzyaaal5(1)+xyzzyaaal5(2)*xyzzyaaal5(2)+xyz&
&zyaaal5(3)*xyzzyaaal5(3)
xyzzyaaak5=xyzzyaaak5+xyzzyaaai1(xyzzyaaaf5,xyzzyaaad5,xyzzyaaac5)**2*&
&xyzzyaaam5
enddo
endif
else
if(xyzzyaaag5==1)then
do xyzzyaaaf5=1,xyzzyaaaf1
xyzzyaaal5=xyzzyaacl1(:,xyzzyaaac5)+xyzzyaaag1(:,xyzzyaaaf5)
xyzzyaaam5=xyzzyaaal5(1)*xyzzyaaal5(1)+xyzzyaaal5(2)*xyzzyaaal5(2)+xyz&
&zyaaal5(3)*xyzzyaaal5(3)
xyzzyaaak5=xyzzyaaak5+dble(xyzzyaaaj1(xyzzyaaaf5,xyzzyaaad5,xyzzyaaac5&
&)*conjg(xyzzyaaaj1(xyzzyaaaf5,xyzzyaaad5,xyzzyaaac5)))*xyzzyaaam5
enddo
else
do xyzzyaaaf5=1,xyzzyaaaf1
xyzzyaaal5=xyzzyaacl1(:,xyzzyaaac5)+xyzzyaaag1(:,xyzzyaaaf5)
xyzzyaaam5=xyzzyaaal5(1)*xyzzyaaal5(1)+xyzzyaaal5(2)*xyzzyaaal5(2)+xyz&
&zyaaal5(3)*xyzzyaaal5(3)
xyzzyaaak5=xyzzyaaak5+dble(xyzzyaaak1(xyzzyaaaf5,xyzzyaaad5,xyzzyaaac5&
&)*conjg(xyzzyaaak1(xyzzyaaaf5,xyzzyaaad5,xyzzyaaac5)))*xyzzyaaam5
enddo
endif
endif
xyzzyaaaj5=xyzzyaaaj5+dble(xyzzyaacf1(xyzzyaaaa5,xyzzyaaad5,xyzzyaaac5&
&,xyzzyaaab5))*xyzzyaaak5
endif
enddo
enddo
enddo
xyzzyaaaj5=0.5d0*xyzzyaaaj5/dble(npcells)
if(ndet>1)then
call wout('Calculated kinetic energy of determinant '//trim(i2s(xyzzya&
&aaa5))//' : '//trim(d2s(xyzzyaaaj5)))
else
call wout('Calculated kinetic energy                     : '//trim(d2s&
&(xyzzyaaaj5)))
if(abs(xyzzyaacb1-xyzzyaaaj5)>xyzzyaaao5*abs(xyzzyaacb1))xyzzyaaar5=.t&
&rue.
endif
enddo
if(xyzzyaaas5)then
if(xyzzyaaar5)then
if(kwarn)then
call errwarn_silent('PWFDET_CF_DFT_KE','kinetic energy mismatch betwee&
&n calculated KE and the value in pwfn.data.  Sometimes this can''t be&
& helped e.g. fractional occupation numbers in the DFT calculation. Th&
&e program will continue since you have activated the KWARN flag in th&
&e input file.')
else
call errstop('PWFDET_CF_DFT_KE','Kinetic energy mismatch. Probably som&
&e kind of state occupation problem (inevitable if using fractional oc&
&cupation in the DFT). Activate the KWARN flag in input to disable thi&
&s check.')
endif
else
call wout('KE check passed.')
endif
endif
call wout()
end subroutine xyzzyaadr1
subroutine xyzzyaads1(sum_orbs_real2,sum_orbs_im2,spin)
use slaarnabt, only : zdotu,zdot_rc
implicit none
integer,intent(in) :: spin
real(dp),intent(out) :: sum_orbs_real2(xyzzyaacd1,xyzzyaacc1),sum_orbs&
&_im2(xyzzyaacd1,xyzzyaacc1)
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaa&
&af6,xyzzyaaag6
real(dp) xyzzyaaah6,xyzzyaaai6,xyzzyaaaj6(3),xyzzyaaak6(3)
complex(dp) xyzzyaaal6,xyzzyaaam6
complex(dp),allocatable :: xyzzyaaan6(:)
allocate(xyzzyaaan6(xyzzyaaaf1),stat=xyzzyaaad6)
call check_alloc(xyzzyaaad6,'SUM_ORBS_OVER_GRID_WITH_PW','orbegdr')
sum_orbs_real2=0.d0
sum_orbs_im2=0.d0
do xyzzyaaae6=-3,3
xyzzyaaaj6(1)=dble(xyzzyaaae6)*12.3456d0-0.01234d0
do xyzzyaaaf6=-3,3
xyzzyaaaj6(2)=dble(xyzzyaaaf6)*12.3456d0-0.04567d0
do xyzzyaaag6=-3,3
xyzzyaaaj6(3)=dble(xyzzyaaag6)*12.3456d0-0.08910d0
xyzzyaaak6=matmul(xyzzyaaaj6,painv)
do xyzzyaaab6=1,xyzzyaaaf1
xyzzyaaah6=xyzzyaaag1(1,xyzzyaaab6)*xyzzyaaak6(1)+xyzzyaaag1(2,xyzzyaa&
&ab6)*xyzzyaaak6(2)+xyzzyaaag1(3,xyzzyaaab6)*xyzzyaaak6(3)
xyzzyaaan6(xyzzyaaab6)=cmplx(cos(xyzzyaaah6),sin(xyzzyaaah6),dp)
enddo
do xyzzyaaac6=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaac6).and.xyzzyaaco1(xyzzyaaac6))then
xyzzyaaai6=xyzzyaacl1(1,xyzzyaaac6)*xyzzyaaak6(1)+xyzzyaacl1(2,xyzzyaa&
&ac6)*xyzzyaaak6(2)+xyzzyaacl1(3,xyzzyaaac6)*xyzzyaaak6(3)
xyzzyaaal6=cmplx(cos(xyzzyaaai6),sin(xyzzyaaai6),dp)
do xyzzyaaaa6=1,xyzzyaace1(xyzzyaaac6,spin)
if(spin==1)then
if(xyzzyaabz1)then
xyzzyaaam6=xyzzyaaal6*zdot_rc(xyzzyaaaf1,xyzzyaaah1(1,xyzzyaaaa6,xyzzy&
&aaac6),1,xyzzyaaan6(1),1)
else
xyzzyaaam6=xyzzyaaal6*zdotu(xyzzyaaaf1,xyzzyaaaj1(1,xyzzyaaaa6,xyzzyaa&
&ac6),1,xyzzyaaan6(1),1)
endif
else
if(xyzzyaabz1)then
xyzzyaaam6=xyzzyaaal6*zdot_rc(xyzzyaaaf1,xyzzyaaai1(1,xyzzyaaaa6,xyzzy&
&aaac6),1,xyzzyaaan6(1),1)
else
xyzzyaaam6=xyzzyaaal6*zdotu(xyzzyaaaf1,xyzzyaaak1(1,xyzzyaaaa6,xyzzyaa&
&ac6),1,xyzzyaaan6(1),1)
endif
endif
sum_orbs_real2(xyzzyaaaa6,xyzzyaaac6)=sum_orbs_real2(xyzzyaaaa6,xyzzya&
&aac6)+dble(xyzzyaaam6)**2
sum_orbs_im2(xyzzyaaaa6,xyzzyaaac6)=sum_orbs_im2(xyzzyaaaa6,xyzzyaaac6&
&)+aimag(xyzzyaaam6)**2
enddo
endif
enddo
enddo
enddo
enddo
deallocate(xyzzyaaan6)
end subroutine xyzzyaads1
subroutine pw_orb_eval(rvec,jspin,norb,orbmask,val,fsd,orbval,orbgrad,&
&orblap,orbsderivs)
use slaarnabi, only : use_gpcc,cusp_wfdet_cmplx
use slaarnabt, only : zdotu,zdot_rc
implicit none
integer,intent(in) :: jspin,norb
real(dp),intent(in) :: rvec(3)
real(dp),intent(inout) :: orbval(norb,real1_complex2),orbgrad(3,norb,r&
&eal1_complex2),orblap(norb,real1_complex2)
real(dp),intent(inout),optional :: orbsderivs(6,norb,real1_complex2)
logical,intent(in) :: val,fsd,orbmask(norb)
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7,xyzzyaa&
&af7,xyzzyaaag7
real(dp) xyzzyaaah7,xyzzyaaai7,xyzzyaaaj7,xyzzyaaak7,xyzzyaaal7,xyzzya&
&aam7,xyzzyaaan7
if(jspin==2.and.xyzzyaacm1)then
xyzzyaaaa7=2
else
xyzzyaaaa7=1
endif
xyzzyaaah7=rvec(1)
xyzzyaaai7=rvec(2)
xyzzyaaaj7=rvec(3)
xyzzyaaak7=xyzzyaaah7*pb1(1)+xyzzyaaai7*pb1(2)+xyzzyaaaj7*pb1(3)
xyzzyaaac1(1,1)=cmplx(cos(xyzzyaaak7),sin(xyzzyaaak7),dp)
xyzzyaaac1(1,-1)=conjg(xyzzyaaac1(1,1))
xyzzyaaak7=xyzzyaaah7*pb2(1)+xyzzyaaai7*pb2(2)+xyzzyaaaj7*pb2(3)
xyzzyaaac1(2,1)=cmplx(cos(xyzzyaaak7),sin(xyzzyaaak7),dp)
xyzzyaaac1(2,-1)=conjg(xyzzyaaac1(2,1))
xyzzyaaak7=xyzzyaaah7*pb3(1)+xyzzyaaai7*pb3(2)+xyzzyaaaj7*pb3(3)
xyzzyaaac1(3,1)=cmplx(cos(xyzzyaaak7),sin(xyzzyaaak7),dp)
xyzzyaaac1(3,-1)=conjg(xyzzyaaac1(3,1))
do xyzzyaaag7=2,xyzzyaaaa1
xyzzyaaac1(1,xyzzyaaag7)=xyzzyaaac1(1,xyzzyaaag7-1)*xyzzyaaac1(1,1)
xyzzyaaac1(2,xyzzyaaag7)=xyzzyaaac1(2,xyzzyaaag7-1)*xyzzyaaac1(2,1)
xyzzyaaac1(3,xyzzyaaag7)=xyzzyaaac1(3,xyzzyaaag7-1)*xyzzyaaac1(3,1)
xyzzyaaac1(1,-xyzzyaaag7)=conjg(xyzzyaaac1(1,xyzzyaaag7))
xyzzyaaac1(2,-xyzzyaaag7)=conjg(xyzzyaaac1(2,xyzzyaaag7))
xyzzyaaac1(3,-xyzzyaaag7)=conjg(xyzzyaaac1(3,xyzzyaaag7))
enddo
do xyzzyaaac7=1,xyzzyaaaf1
xyzzyaaae1(xyzzyaaac7)=xyzzyaaac1(1,xyzzyaaab1(1,xyzzyaaac7))*xyzzyaaa&
&c1(2,xyzzyaaab1(2,xyzzyaaac7))*xyzzyaaac1(3,xyzzyaaab1(3,xyzzyaaac7))
enddo
do xyzzyaaad7=1,xyzzyaacc1
if(xyzzyaacn1(xyzzyaaad7))then
xyzzyaaal7=xyzzyaaah7*xyzzyaacl1(1,xyzzyaaad7)+xyzzyaaai7*xyzzyaacl1(2&
&,xyzzyaaad7)+xyzzyaaaj7*xyzzyaacl1(3,xyzzyaaad7)
xyzzyaaad1(xyzzyaaad7)=cmplx(xyzzyaacg1*cos(xyzzyaaal7),xyzzyaacg1*sin&
&(xyzzyaaal7),dp)
endif
enddo
if(use_gpcc)then
if(present(orbsderivs))then
call cusp_wfdet_cmplx(rvec,jspin,val,fsd,xyzzyaact1,xyzzyaacu1,xyzzyaa&
&cv1,xyzzyaacw1)
else
call cusp_wfdet_cmplx(rvec,jspin,val,fsd,xyzzyaact1,xyzzyaacu1,xyzzyaa&
&cv1)
endif
endif
xyzzyaacp1=czero
xyzzyaacq1=czero
xyzzyaacr1=czero
xyzzyaacs1=czero
do xyzzyaaae7=1,xyzzyaacz1
if(.not.xyzzyaadl1(xyzzyaaae7,jspin))cycle
xyzzyaaab7=xyzzyaadc1(xyzzyaaae7)
xyzzyaaad7=xyzzyaadd1(xyzzyaaae7)
if(val)then
if(xyzzyaaaa7==1)then
if(xyzzyaabz1)then
xyzzyaacp1(xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaah1(1,xyzzyaaab7,xyzz&
&yaaad7),1,xyzzyaaae1(1),1)
else
xyzzyaacp1(xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaaaj1(1,x&
&yzzyaaab7,xyzzyaaad7),1)
endif
else
if(xyzzyaabz1)then
xyzzyaacp1(xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaai1(1,xyzzyaaab7,xyzz&
&yaaad7),1,xyzzyaaae1(1),1)
else
xyzzyaacp1(xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaaak1(1,x&
&yzzyaaab7,xyzzyaaad7),1)
endif
endif
if(use_gpcc)then
xyzzyaacp1(xyzzyaaae7)=xyzzyaaad1(xyzzyaaad7)*xyzzyaacp1(xyzzyaaae7)+x&
&yzzyaact1(xyzzyaaae7)
else
xyzzyaacp1(xyzzyaaae7)=xyzzyaaad1(xyzzyaaad7)*xyzzyaacp1(xyzzyaaae7)
endif
endif
if(fsd)then
if(xyzzyaaaa7==1)then
if(xyzzyaabz1)then
xyzzyaacq1(1,xyzzyaaae7)=zi*zdot_rc(xyzzyaaaf1,xyzzyaaal1(1,xyzzyaaab7&
&,xyzzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacq1(2,xyzzyaaae7)=zi*zdot_rc(xyzzyaaaf1,xyzzyaaam1(1,xyzzyaaab7&
&,xyzzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacq1(3,xyzzyaaae7)=zi*zdot_rc(xyzzyaaaf1,xyzzyaaan1(1,xyzzyaaab7&
&,xyzzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacr1(xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaar1(1,xyzzyaaab7,xyzz&
&yaaad7),1,xyzzyaaae1(1),1)
if(present(orbsderivs))then
xyzzyaacs1(1,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaat1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(2,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaau1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(3,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaav1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(4,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaaw1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(5,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaax1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(6,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaay1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
endif
else
xyzzyaacq1(1,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabf1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacq1(2,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabg1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacq1(3,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabh1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacr1(xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabl1(1,x&
&yzzyaaab7,xyzzyaaad7),1)
if(present(orbsderivs))then
xyzzyaacs1(1,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabn1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(2,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabo1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(3,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabp1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(4,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabq1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(5,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabr1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(6,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabs1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
endif
endif
else
if(xyzzyaabz1)then
xyzzyaacq1(1,xyzzyaaae7)=zi*zdot_rc(xyzzyaaaf1,xyzzyaaao1(1,xyzzyaaab7&
&,xyzzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacq1(2,xyzzyaaae7)=zi*zdot_rc(xyzzyaaaf1,xyzzyaaap1(1,xyzzyaaab7&
&,xyzzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacq1(3,xyzzyaaae7)=zi*zdot_rc(xyzzyaaaf1,xyzzyaaaq1(1,xyzzyaaab7&
&,xyzzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacr1(xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaas1(1,xyzzyaaab7,xyzz&
&yaaad7),1,xyzzyaaae1(1),1)
if(present(orbsderivs))then
xyzzyaacs1(1,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaaz1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(2,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaaba1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(3,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaabb1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(4,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaabc1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(5,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaabd1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
xyzzyaacs1(6,xyzzyaaae7)=zdot_rc(xyzzyaaaf1,xyzzyaabe1(1,xyzzyaaab7,xy&
&zzyaaad7),1,xyzzyaaae1(1),1)
endif
else
xyzzyaacq1(1,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabi1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacq1(2,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabj1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacq1(3,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabk1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacr1(xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabm1(1,x&
&yzzyaaab7,xyzzyaaad7),1)
if(present(orbsderivs))then
xyzzyaacs1(1,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabt1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(2,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabu1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(3,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabv1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(4,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabw1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(5,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaabx1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
xyzzyaacs1(6,xyzzyaaae7)=zdotu(xyzzyaaaf1,xyzzyaaae1(1),1,xyzzyaaby1(1&
&,xyzzyaaab7,xyzzyaaad7),1)
endif
endif
endif
if(use_gpcc)then
xyzzyaacq1(1:3,xyzzyaaae7)=xyzzyaaad1(xyzzyaaad7)*xyzzyaacq1(1:3,xyzzy&
&aaae7)+xyzzyaacu1(1:3,xyzzyaaae7)
xyzzyaacr1(xyzzyaaae7)=xyzzyaaad1(xyzzyaaad7)*xyzzyaacr1(xyzzyaaae7)+x&
&yzzyaacv1(xyzzyaaae7)
if(present(orbsderivs))xyzzyaacs1(1:6,xyzzyaaae7)=xyzzyaaad1(xyzzyaaad&
&7)*xyzzyaacs1(1:6,xyzzyaaae7)+xyzzyaacw1(1:6,xyzzyaaae7)
else
xyzzyaacq1(1:3,xyzzyaaae7)=xyzzyaaad1(xyzzyaaad7)*xyzzyaacq1(1:3,xyzzy&
&aaae7)
xyzzyaacr1(xyzzyaaae7)=xyzzyaaad1(xyzzyaaad7)*xyzzyaacr1(xyzzyaaae7)
if(present(orbsderivs))xyzzyaacs1(1:6,xyzzyaaae7)=xyzzyaaad1(xyzzyaaad&
&7)*xyzzyaacs1(1:6,xyzzyaaae7)
endif
endif
enddo
if(xyzzyaaca1)then
do xyzzyaaaf7=1,xyzzyaacx1
if(.not.orbmask(xyzzyaaaf7))cycle
if(val)then
orbval(xyzzyaaaf7,1)=dble(xyzzyaacp1(xyzzyaaaf7))
orbval(xyzzyaaaf7,2)=aimag(xyzzyaacp1(xyzzyaaaf7))
endif
if(fsd)then
orbgrad(1:3,xyzzyaaaf7,1)=dble(xyzzyaacq1(1:3,xyzzyaaaf7))
orbgrad(1:3,xyzzyaaaf7,2)=aimag(xyzzyaacq1(1:3,xyzzyaaaf7))
orblap(xyzzyaaaf7,1)=dble(xyzzyaacr1(xyzzyaaaf7))
orblap(xyzzyaaaf7,2)=aimag(xyzzyaacr1(xyzzyaaaf7))
if(present(orbsderivs))then
orbsderivs(1:6,xyzzyaaaf7,1)=dble(xyzzyaacs1(1:6,xyzzyaaaf7))
orbsderivs(1:6,xyzzyaaaf7,2)=aimag(xyzzyaacs1(1:6,xyzzyaaaf7))
endif
endif
enddo
else
do xyzzyaaaf7=1,xyzzyaacx1
if(.not.orbmask(xyzzyaaaf7))cycle
xyzzyaaae7=xyzzyaadb1(xyzzyaaaf7)
if(xyzzyaadg1(xyzzyaaaf7)==xyzzyaadh1)then
if(val)orbval(xyzzyaaaf7,1)=dble(xyzzyaacp1(xyzzyaaae7))
if(fsd)then
orbgrad(1:3,xyzzyaaaf7,1)=dble(xyzzyaacq1(1:3,xyzzyaaae7))
orblap(xyzzyaaaf7,1)=dble(xyzzyaacr1(xyzzyaaae7))
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaaf7,1)=dble(xyzzyaacs1(1:&
&6,xyzzyaaae7))
endif
elseif(xyzzyaadg1(xyzzyaaaf7)==xyzzyaadi1)then
if(val)orbval(xyzzyaaaf7,1)=aimag(xyzzyaacp1(xyzzyaaae7))
if(fsd)then
orbgrad(1:3,xyzzyaaaf7,1)=aimag(xyzzyaacq1(1:3,xyzzyaaae7))
orblap(xyzzyaaaf7,1)=aimag(xyzzyaacr1(xyzzyaaae7))
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaaf7,1)=aimag(xyzzyaacs1(1&
&:6,xyzzyaaae7))
endif
else
xyzzyaaam7=dble(xyzzyaadf1(xyzzyaaaf7))
xyzzyaaan7=aimag(xyzzyaadf1(xyzzyaaaf7))
if(val)orbval(xyzzyaaaf7,1)=dble(xyzzyaacp1(xyzzyaaae7))*xyzzyaaam7+ai&
&mag(xyzzyaacp1(xyzzyaaae7))*xyzzyaaan7
if(fsd)then
orbgrad(1:3,xyzzyaaaf7,1)=dble(xyzzyaacq1(1:3,xyzzyaaae7))*xyzzyaaam7+&
&aimag(xyzzyaacq1(1:3,xyzzyaaae7))*xyzzyaaan7
orblap(xyzzyaaaf7,1)=dble(xyzzyaacr1(xyzzyaaae7))*xyzzyaaam7+aimag(xyz&
&zyaacr1(xyzzyaaae7))*xyzzyaaan7
if(present(orbsderivs))orbsderivs(1:6,xyzzyaaaf7,1)=dble(xyzzyaacs1(1:&
&6,xyzzyaaae7))*xyzzyaaam7+aimag(xyzzyaacs1(1:6,xyzzyaaae7))*xyzzyaaan&
&7
endif
endif
enddo
endif
if(complex_wf.and..not.xyzzyaaca1)then
if(val)orbval(:,2)=0.d0
if(fsd)then
orbgrad(:,:,2)=0.d0
orblap(:,2)=0.d0
if(present(orbsderivs))orbsderivs(:,:,2)=0.d0
endif
endif
end subroutine pw_orb_eval
subroutine get_pwfdet_orbmap(row_offset,norb,orbmap)
implicit none
integer,intent(inout) :: row_offset(nspin),norb,orbmap(nemax,nspin,nde&
&t)
integer xyzzyaaaa8
do xyzzyaaaa8=1,nspin
orbmap(row_offset(xyzzyaaaa8)+1:row_offset(xyzzyaaaa8)+nuc_nele(xyzzya&
&aaa8),xyzzyaaaa8,:)=xyzzyaada1(1:nuc_nele(xyzzyaaaa8),xyzzyaaaa8,:)
row_offset(xyzzyaaaa8)=row_offset(xyzzyaaaa8)+nuc_nele(xyzzyaaaa8)
enddo
norb=norb+xyzzyaacx1
end subroutine get_pwfdet_orbmap
subroutine get_pwfdet_ndesc(ndesc_int,ndesc_dp)
implicit none
integer,intent(inout) :: ndesc_int,ndesc_dp
ndesc_int=xyzzyaadm1
ndesc_dp=xyzzyaadn1
end subroutine get_pwfdet_ndesc
subroutine get_pwfdet_orbdesc(norb,ndesc_int,ndesc_dp,orbdesc_int,orbd&
&esc_dp)
implicit none
integer,intent(in) :: norb,ndesc_int,ndesc_dp
integer,intent(inout) :: orbdesc_int(ndesc_int,norb)
real(dp),intent(inout) :: orbdesc_dp(ndesc_dp,norb)
orbdesc_int(1:xyzzyaadm1,1:xyzzyaacx1)=xyzzyaado1(1:xyzzyaadm1,1:xyzzy&
&aacx1)
orbdesc_dp(1:xyzzyaadn1,1:xyzzyaacx1)=xyzzyaadp1(1:xyzzyaadn1,1:xyzzya&
&acx1)
end subroutine get_pwfdet_orbdesc
end module slaarnacb

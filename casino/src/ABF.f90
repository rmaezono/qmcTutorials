module slaarnabf
use dsp
use parallel,    only : am_master
use slaarnaag,   only : pi,twopi,fourpi,one_over_twopi,twothirds,root_&
&two,czero
use file_utils,  only : open_units
use format_utils,only : wout,i2s,r2s,r2s2
use slaarnabg,    only : a1,a2,a3,pa1,pa2,pa3,pb1,pb2,pb3,wigner_seitz&
&_radius,periodicity,scell_matrix,npcells,nbasis,atno,basis,dimensiona&
&lity,atom_basis_type,amat,ainv,binv,volume,area
use slaarnabt,   only : quicksort,ddot,zscal,dscal
use slaarnaca,       only : have_veep
use run_control, only : errstop,errstop_master,errwarn,timer,check_all&
&oc
use singleton,   only : fftn
use store,       only : nele,nemax,ndet,netot,nspin,complex_wf,noncoll&
&_spin,real1_complex2,localized_orbitals,self_term,open_unit
implicit none
private
public run_mpc_generation
real(dp) xyzzyaaaa1
contains
subroutine run_mpc_generation(mpc_cutoff)
implicit none
real(dp),intent(inout) :: mpc_cutoff
integer xyzzyaaaa2,xyzzyaaab2
integer,parameter :: xyzzyaaac2=1
integer,parameter :: xyzzyaaad2=30000,xyzzyaaae2=50
character(80) tmpr
if(mpc_cutoff<=0.d0)then
if(periodicity==3)then
mpc_cutoff=(3.d0*pi**2*dble(xyzzyaaad2)/(root_two*volume))**twothirds
elseif(periodicity==2)then
mpc_cutoff=twopi*dble(xyzzyaaae2)/area
else
call errstop('RUN_MPC_GENERATION','MPC interaction is not available fo&
&r 1D-periodic systems.')
endif
endif
if(am_master)then
call open_units(xyzzyaaaa2,xyzzyaaab2)
if(xyzzyaaab2/=0)call errstop('RUN_MPC_GENERATION','Unable to find fre&
&e io unit')
open(xyzzyaaaa2,file='mpc.data',form='formatted',status='unknown',iost&
&at=xyzzyaaab2)
if(xyzzyaaab2/=0)call errstop('RUN_MPC_GENERATION','Unable to open mpc&
&.data.')
call wout('Generation of mpc.data file')
call wout('===========================')
call wout()
tmpr=r2s(mpc_cutoff,'(f16.4)')
call wout('Using energy cutoff of '//trim(tmpr)//' au.')
call wout()
write(xyzzyaaaa2,'(a)')'START MPC DATA'
write(xyzzyaaaa2,'(a)')'Title'
write(xyzzyaaaa2,'(a)')'  No title given'
write(xyzzyaaaa2,'(a)')'File version'
write(xyzzyaaaa2,'(1x,a)')trim(i2s(xyzzyaaac2))
write(xyzzyaaaa2,*)
endif
call timer('DENSITY_FFT',.true.)
if(am_master)call wout('Fourier transforming the density...')
call xyzzyaaab1(xyzzyaaaa2,mpc_cutoff)
if(am_master)call wout('Done.')
call timer('DENSITY_FFT',.false.)
call timer('EEPOT_FFT',.true.)
if(am_master)call wout('Fourier transforming 1/r...')
call xyzzyaaac1(xyzzyaaaa2,mpc_cutoff)
if(am_master)call wout('Done.')
call timer('EEPOT_FFT',.false.)
if(am_master)then
write(xyzzyaaaa2,'(a)')
write(xyzzyaaaa2,'(a)')'END MPC DATA'
call wout()
call wout('Written mpc.data file.')
call wout()
have_veep=.false.
close(xyzzyaaaa2)
open_unit(xyzzyaaaa2)=.false.
endif
end subroutine run_mpc_generation
subroutine xyzzyaaab1(io_mpc,mpc_cutoff)
use parallel
use slaarnacq, only : wfdet,wfdet_norb,wfdet_orbmap,wfdet_orbmask,copy&
&_orb_to_det
implicit none
integer,intent(in) :: io_mpc
real(dp),intent(in) :: mpc_cutoff
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3(3),xyzzyaaak3,xyzzyaa&
&al3,xyzzyaaam3,xyzzyaaan3,xyzzyaaao3(3),xyzzyaaap3
integer,allocatable :: xyzzyaaaq3(:),xyzzyaaar3(:),xyzzyaaas3(:)
real(dp) xyzzyaaat3(3),xyzzyaaau3,xyzzyaaav3,xyzzyaaaw3,xyzzyaaax3(3),&
&xyzzyaaay3,xyzzyaaaz3,xyzzyaaba3,xyzzyaabb3,xyzzyaabc3,xyzzyaabd3,xyz&
&zyaabe3,xyzzyaabf3,xyzzyaabg3
real(dp),allocatable :: xyzzyaabh3(:),xyzzyaabi3(:,:),xyzzyaabj3(:,:),&
&xyzzyaabk3(:,:,:),xyzzyaabl3(:,:,:,:),xyzzyaabm3(:,:,:),xyzzyaabn3(:,&
&:),xyzzyaabo3(:,:,:),xyzzyaabp3(:,:)
complex(dp) xyzzyaabq3
complex(dp),allocatable :: xyzzyaabr3(:,:,:),xyzzyaabs3(:,:,:),xyzzyaa&
&bt3(:)
character(80) tmpr,tmpr2,tmpr3
if(ndet>1)call errstop_master('GENERATE_DENSITY','At present the singl&
&e-particle charge density can only be calculated using a single-deter&
&minant wave function.')
if(localized_orbitals)call errstop_master('GENERATE_DENSITY','The sing&
&le-particle charge density can only be calculated using orthogonal or&
&bitals.')
xyzzyaabd3=sqrt(dot_product(pa1,pa1))
xyzzyaabe3=sqrt(dot_product(pa2,pa2))
xyzzyaabf3=sqrt(dot_product(pa3,pa3))
call xyzzyaaah1(xyzzyaabd3,xyzzyaabe3,xyzzyaabf3,mpc_cutoff,xyzzyaaao3&
&)
if(periodicity==2)then
if(dimensionality==2)then
xyzzyaaao3(3)=1
elseif(dimensionality==3)then
xyzzyaaay3=sqrt(pa1(1)**2+pa1(2)**2)
xyzzyaaaz3=sqrt(pa2(1)**2+pa2(2)**2)
if(pa3(3)>10.d0*max(xyzzyaaay3,xyzzyaaaz3))call errwarn('GENERATE_DENS&
&ITY','Artificial periodicity in z-direction is very large; reduce it &
&if possible.')
else
call errstop_master('GENERATE_DENSITY','Single-particle density cannot&
& yet be calculated for 1D systems.')
endif
endif
if(periodicity/=2.and.periodicity/=3)then
call errstop_master('GENERATE_DENSITY','The density generator was call&
&ed with an unsupported periodicity.')
endif
xyzzyaaap3=xyzzyaaao3(1)*xyzzyaaao3(2)*xyzzyaaao3(3)
if(am_master)call wout('Density generator using '//trim(i2s(xyzzyaaao3&
&(1)))//'x'//trim(i2s(xyzzyaaao3(2)))//'x'//trim(i2s(xyzzyaaao3(3)))//&
&' grid.')
if(xyzzyaaap3>2**19)call errwarn('GENERATE_DENSITY','The density gener&
&ator grid is extremely large. Consider reducing MPC_CUTOFF.')
xyzzyaaai3=xyzzyaaaf1(xyzzyaabd3,xyzzyaabe3,xyzzyaabf3,pb1,pb2,pb3,mpc&
&_cutoff)
allocate(xyzzyaaas3(0:nnodes-1),xyzzyaaar3(0:nnodes-1),stat=xyzzyaaah3&
&)
call check_alloc(xyzzyaaah3,'GENERATE_DENSITY','limits')
xyzzyaaak3=xyzzyaaap3/nnodes
xyzzyaaal3=mod(xyzzyaaap3,nnodes)
xyzzyaaas3(0)=1
if(xyzzyaaal3>0)then
xyzzyaaar3(0)=xyzzyaaak3+1
else
xyzzyaaar3(0)=xyzzyaaak3
endif
do xyzzyaaam3=1,nnodes-1
xyzzyaaas3(xyzzyaaam3)=xyzzyaaar3(xyzzyaaam3-1)+1
if(xyzzyaaam3<xyzzyaaal3)then
xyzzyaaar3(xyzzyaaam3)=xyzzyaaas3(xyzzyaaam3)+xyzzyaaak3
else
xyzzyaaar3(xyzzyaaam3)=xyzzyaaas3(xyzzyaaam3)+xyzzyaaak3-1
endif
enddo
xyzzyaaay3=1.d0/dble(xyzzyaaao3(1))
xyzzyaaaz3=1.d0/dble(xyzzyaaao3(2))
xyzzyaaba3=1.d0/dble(xyzzyaaao3(3))
allocate(xyzzyaabn3(wfdet_norb,real1_complex2),xyzzyaabo3(3,wfdet_norb&
&,real1_complex2),xyzzyaabp3(wfdet_norb,real1_complex2),stat=xyzzyaaah&
&3)
call check_alloc(xyzzyaaah3,'GENERATE_DENSITY','orbval_t,...')
allocate(xyzzyaabk3(nemax,real1_complex2,ndet),xyzzyaabl3(3,nemax,real&
&1_complex2,ndet),xyzzyaabm3(nemax,real1_complex2,ndet),xyzzyaabs3(0:x&
&yzzyaaao3(1)-1,0:xyzzyaaao3(2)-1,0:xyzzyaaao3(3)-1),xyzzyaabj3(4,neto&
&t),stat=xyzzyaaah3)
call check_alloc(xyzzyaaah3,'GENERATE_DENSITY','rpsi,...')
xyzzyaabj3=0.d0
xyzzyaaaa3=1
do xyzzyaaad3=0,xyzzyaaao3(3)-1
xyzzyaaaw3=xyzzyaaba3*dble(xyzzyaaad3)
do xyzzyaaac3=0,xyzzyaaao3(2)-1
xyzzyaaav3=xyzzyaaaz3*dble(xyzzyaaac3)
do xyzzyaaab3=0,xyzzyaaao3(1)-1
if(xyzzyaaaa3>=xyzzyaaas3(my_node).and.xyzzyaaaa3<=xyzzyaaar3(my_node)&
&)then
xyzzyaabc3=0.d0
xyzzyaaau3=xyzzyaaay3*dble(xyzzyaaab3)
xyzzyaaax3=xyzzyaaau3*pa1+xyzzyaaav3*pa2+xyzzyaaaw3*pa3
if(trim(atom_basis_type)=='gaussian'.and.periodicity==2.and.xyzzyaaax3&
&(3)>0.75d0*pa3(3))xyzzyaaax3(3)=xyzzyaaax3(3)-pa3(3)
if(trim(atom_basis_type)/='gaussian'.or.periodicity/=2.or.abs(xyzzyaaa&
&x3(3))<0.25d0*pa3(3))then
do xyzzyaaan3=1,nspin
if(nele(xyzzyaaan3)>0)then
call wfdet(xyzzyaaax3,1,xyzzyaaan3,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&n3),.true.,.false.,xyzzyaabn3,xyzzyaabo3,xyzzyaabp3,xyzzyaabj3)
call copy_orb_to_det(1,xyzzyaaan3,wfdet_orbmap,xyzzyaabn3,xyzzyaabk3)
xyzzyaabc3=xyzzyaabc3+ddot(nele(xyzzyaaan3),xyzzyaabk3(1,1,1),1,xyzzya&
&abk3(1,1,1),1)
if(complex_wf)xyzzyaabc3=xyzzyaabc3+ddot(nele(xyzzyaaan3),xyzzyaabk3(1&
&,2,1),1,xyzzyaabk3(1,2,1),1)
if(noncoll_spin)then
call wfdet(xyzzyaaax3,2,xyzzyaaan3,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&n3),.true.,.false.,xyzzyaabn3,xyzzyaabo3,xyzzyaabp3,xyzzyaabj3)
call copy_orb_to_det(1,xyzzyaaan3,wfdet_orbmap,xyzzyaabn3,xyzzyaabk3)
xyzzyaabc3=xyzzyaabc3+ddot(nele(xyzzyaaan3),xyzzyaabk3(1,1,1),1,xyzzya&
&abk3(1,1,1),1)
if(complex_wf)xyzzyaabc3=xyzzyaabc3+ddot(nele(xyzzyaaan3),xyzzyaabk3(1&
&,2,1),1,xyzzyaabk3(1,2,1),1)
endif
endif
enddo
endif
xyzzyaabs3(xyzzyaaab3,xyzzyaaac3,xyzzyaaad3)=cmplx(xyzzyaabc3,0.d0,dp)
else
xyzzyaabs3(xyzzyaaab3,xyzzyaaac3,xyzzyaaad3)=czero
endif
xyzzyaaaa3=xyzzyaaaa3+1
enddo
enddo
enddo
deallocate(xyzzyaaas3,xyzzyaaar3,xyzzyaabj3,xyzzyaabk3,xyzzyaabm3,xyzz&
&yaabl3)
allocate(xyzzyaabr3(0:xyzzyaaao3(1)-1,0:xyzzyaaao3(2)-1,0:xyzzyaaao3(3&
&)-1),stat=xyzzyaaah3)
call check_alloc(xyzzyaaah3,'GENERATE_DENSITY','rhor')
call mpi_reduce(xyzzyaabs3,xyzzyaabr3,xyzzyaaap3,mpi_double_complex,mp&
&i_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Collecting density data in denft.')
deallocate(xyzzyaabs3)
if(am_master)then
call timer('DENSITY_FFT',.true.)
call fftn(xyzzyaabr3,xyzzyaaao3,inv=.true.)
call timer('DENSITY_FFT',.false.)
allocate(xyzzyaabi3(3,xyzzyaaai3),xyzzyaabt3(xyzzyaaai3),xyzzyaabh3(xy&
&zzyaaai3),stat=xyzzyaaah3)
call check_alloc(xyzzyaaah3,'GENERATE_DENSITY','gvec')
xyzzyaabg3=2.d0*mpc_cutoff
xyzzyaaae3=0
xyzzyaaaj3=xyzzyaaao3/2
do xyzzyaaad3=0,xyzzyaaao3(3)-1
if(xyzzyaaad3>xyzzyaaaj3(3))then
xyzzyaaaw3=dble(xyzzyaaad3-xyzzyaaao3(3))
else
xyzzyaaaw3=dble(xyzzyaaad3)
endif
do xyzzyaaac3=0,xyzzyaaao3(2)-1
if(xyzzyaaac3>xyzzyaaaj3(2))then
xyzzyaaav3=dble(xyzzyaaac3-xyzzyaaao3(2))
else
xyzzyaaav3=dble(xyzzyaaac3)
endif
do xyzzyaaab3=0,xyzzyaaao3(1)-1
if(xyzzyaaab3>xyzzyaaaj3(1))then
xyzzyaaau3=dble(xyzzyaaab3-xyzzyaaao3(1))
else
xyzzyaaau3=dble(xyzzyaaab3)
endif
xyzzyaaat3=xyzzyaaau3*pb1+xyzzyaaav3*pb2+xyzzyaaaw3*pb3
xyzzyaabb3=dot_product(xyzzyaaat3,xyzzyaaat3)
if(xyzzyaabb3<xyzzyaabg3)then
xyzzyaaae3=xyzzyaaae3+1
xyzzyaabt3(xyzzyaaae3)=xyzzyaabr3(xyzzyaaab3,xyzzyaaac3,xyzzyaaad3)
xyzzyaabi3(1:3,xyzzyaaae3)=xyzzyaaat3
xyzzyaabh3(xyzzyaaae3)=xyzzyaabb3
endif
enddo
enddo
enddo
deallocate(xyzzyaabr3)
xyzzyaabq3=cmplx(dble(sum(nele))/dble(npcells),0.d0,dp)/xyzzyaabt3(1)
call zscal(xyzzyaaai3,xyzzyaabq3,xyzzyaabt3(1),1)
allocate(xyzzyaaaq3(xyzzyaaai3),stat=xyzzyaaah3)
call check_alloc(xyzzyaaah3,'GENERATE_DENSITY','gindx')
call quicksort(xyzzyaaai3,xyzzyaabh3(1),xyzzyaaaq3(1))
write(io_mpc,'(a)')'START DENSITY DATA'
write(io_mpc,'(a)')'Self-consistent charge density in reciprocal space&
&'
write(io_mpc,'(a)')'Real space primitive cell translation vectors (au)&
&'
tmpr=r2s2(pa1(1),'(es24.16)')
tmpr2=r2s2(pa1(2),'(es24.16)')
tmpr3=r2s2(pa1(3),'(es24.16)')
write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
tmpr=r2s2(pa2(1),'(es24.16)')
tmpr2=r2s2(pa2(2),'(es24.16)')
tmpr3=r2s2(pa2(3),'(es24.16)')
write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
tmpr=r2s2(pa3(1),'(es24.16)')
tmpr2=r2s2(pa3(2),'(es24.16)')
tmpr3=r2s2(pa3(3),'(es24.16)')
write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
write(io_mpc,'(a)')'Number of atoms in the primitive cell'
write(io_mpc,'(1x,a)')trim(i2s(nbasis))
write(io_mpc,'(a)')'Positions of atoms (au)'
do xyzzyaaaf3=1,nbasis
tmpr=r2s2(basis(1,xyzzyaaaf3),'(es24.16)')
tmpr2=r2s2(basis(2,xyzzyaaaf3),'(es24.16)')
tmpr3=r2s2(basis(3,xyzzyaaaf3),'(es24.16)')
write(io_mpc,'(1x,4(1x,a))')trim(i2s(atno(xyzzyaaaf3))),trim(tmpr),tri&
&m(tmpr2),trim(tmpr3)
enddo
write(io_mpc,'(a)')'Energy cutoff used for G-vectors (au)'
tmpr=r2s(mpc_cutoff,'(es24.16)')
write(io_mpc,'(2x,a)')trim(tmpr)
write(io_mpc,'(a)')'Number of particle types (1=electrons,2=electrons/&
&holes)'
write(io_mpc,*)' 1'
write(io_mpc,'(a)')'Number of G-vectors'
write(io_mpc,'(1x,a)')trim(i2s(xyzzyaaai3))
write(io_mpc,'(a)')'Primitive cell G-vectors (au)'
do xyzzyaaae3=1,xyzzyaaai3
xyzzyaaag3=xyzzyaaaq3(xyzzyaaae3)
tmpr=r2s2(xyzzyaabi3(1,xyzzyaaag3),'(es24.16)')
tmpr2=r2s2(xyzzyaabi3(2,xyzzyaaag3),'(es24.16)')
tmpr3=r2s2(xyzzyaabi3(3,xyzzyaaag3),'(es24.16)')
write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
enddo
write(io_mpc,'(a)')'START SET 1'
write(io_mpc,'(a)')'Particle type'
write(io_mpc,'(a)')'  1'
write(io_mpc,'(a)')'Complex charge density (real part, imaginary part)&
&'
do xyzzyaaae3=1,xyzzyaaai3
xyzzyaaag3=xyzzyaaaq3(xyzzyaaae3)
tmpr=r2s2(dble(xyzzyaabt3(xyzzyaaag3)),'(es24.16)')
tmpr2=r2s2(aimag(xyzzyaabt3(xyzzyaaag3)),'(es24.16)')
write(io_mpc,'(2(1x,a))')trim(tmpr),trim(tmpr2)
enddo
write(io_mpc,'(a)')'END SET 1'
write(io_mpc,'(a)')'END DENSITY DATA'
write(io_mpc,*)
deallocate(xyzzyaabt3,xyzzyaabi3,xyzzyaabh3,xyzzyaaaq3)
else
deallocate(xyzzyaabr3)
endif
end subroutine xyzzyaaab1
subroutine xyzzyaaac1(io_mpc,mpc_cutoff)
use slaarnabk, only : init_interactions
use parallel,     only : am_master
implicit none
integer,intent(in) :: io_mpc
real(dp),intent(in) :: mpc_cutoff
integer xyzzyaaaa4,xyzzyaaab4,xyzzyaaac4,xyzzyaaad4
real(dp) xyzzyaaae4,xyzzyaaaf4,xyzzyaaag4
real(dp),allocatable :: xyzzyaaah4(:),xyzzyaaai4(:),xyzzyaaaj4(:,:)
character(80) tmpr,tmpr2,tmpr3,tmpr4
if(.not.am_master)return
if(periodicity==2)then
if(trim(atom_basis_type)/='gaussian'.and.trim(atom_basis_type)/='blip'&
&.and.trim(atom_basis_type)/='none')call errstop('GENERATE_MPC','This &
&routine works in 2D only for finite-width slabs. Somebody please fix &
&this.')
call wout('Need Ewald interaction for 2D case.')
call wout()
call init_interactions
elseif(periodicity/=3)then
call errstop('GENERATE_MPC','This routine was called with an unsupport&
&ed periodicity')
endif
xyzzyaaaa1=wigner_seitz_radius**2
xyzzyaaae4=sqrt(dot_product(pa1,pa1))
xyzzyaaaf4=sqrt(dot_product(pa2,pa2))
xyzzyaaag4=sqrt(dot_product(pa3,pa3))
xyzzyaaac4=xyzzyaaaf1(xyzzyaaae4,xyzzyaaaf4,xyzzyaaag4,pb1,pb2,pb3,mpc&
&_cutoff)
allocate(xyzzyaaah4(xyzzyaaac4),xyzzyaaaj4(3,xyzzyaaac4),xyzzyaaai4(xy&
&zzyaaac4),stat=xyzzyaaab4)
call check_alloc(xyzzyaaab4,'GENERATE_MPC','')
call xyzzyaaag1(xyzzyaaae4,xyzzyaaaf4,xyzzyaaag4,pb1,pb2,pb3,mpc_cutof&
&f,xyzzyaaac4,xyzzyaaaj4,xyzzyaaai4)
call xyzzyaaad1(mpc_cutoff,xyzzyaaac4,xyzzyaaaj4,xyzzyaaah4)
if(periodicity==3)then
xyzzyaaah4(1)=xyzzyaaah4(1)*volume+twopi/5.d0*xyzzyaaaa1
do xyzzyaaaa4=2,xyzzyaaac4
xyzzyaaah4(xyzzyaaaa4)=xyzzyaaah4(xyzzyaaaa4)*volume+fourpi/xyzzyaaai4&
&(xyzzyaaaa4)**2+3.d0*fourpi/(xyzzyaaai4(xyzzyaaaa4)**4*xyzzyaaaa1)*(c&
&os(xyzzyaaai4(xyzzyaaaa4)*wigner_seitz_radius)-sin(xyzzyaaai4(xyzzyaa&
&aa4)*wigner_seitz_radius)/(xyzzyaaai4(xyzzyaaaa4)*wigner_seitz_radius&
&))
enddo
else
xyzzyaaah4=xyzzyaaah4*dble(npcells)
endif
call wout()
call wout('Total number of G vectors                              : '/&
&/trim(i2s(xyzzyaaac4)))
write(io_mpc,'(a)',iostat=xyzzyaaad4)'START EEPOT DATA'
if(xyzzyaaad4/=0)call errstop('GENERATE_MPC','Error writing to mpc.dat&
&a.')
write(io_mpc,'(a)')'Real space simulation cell translation vectors (au&
&)'
tmpr=r2s2(a1(1),'(es24.16)')
tmpr2=r2s2(a1(2),'(es24.16)')
tmpr3=r2s2(a1(3),'(es24.16)')
write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
tmpr=r2s2(a2(1),'(es24.16)')
tmpr2=r2s2(a2(2),'(es24.16)')
tmpr3=r2s2(a2(3),'(es24.16)')
write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
tmpr=r2s2(a3(1),'(es24.16)')
tmpr2=r2s2(a3(2),'(es24.16)')
tmpr3=r2s2(a3(3),'(es24.16)')
write(io_mpc,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
write(io_mpc,'(a)')'Diagonal elements of supercell matrix'
write(io_mpc,'(1x,3(1x,a))')trim(i2s(scell_matrix(1,1))),trim(i2s(scel&
&l_matrix(2,2))),trim(i2s(scell_matrix(3,3)))
write(io_mpc,'(a)')'Energy cutoff used for G-vectors (au)'
tmpr=r2s(mpc_cutoff,'(es24.16)')
write(io_mpc,'(2x,a)')trim(tmpr)
write(io_mpc,'(a)')'Number of G-vectors'
write(io_mpc,'(1x,a)')trim(i2s(xyzzyaaac4))
write(io_mpc,'(a)')'Primitive-cell G-vectors (au) and Fourier coeffici&
&ents'
do xyzzyaaaa4=1,xyzzyaaac4
tmpr=r2s2(xyzzyaaaj4(1,xyzzyaaaa4),'(es24.16)')
tmpr2=r2s2(xyzzyaaaj4(2,xyzzyaaaa4),'(es24.16)')
tmpr3=r2s2(xyzzyaaaj4(3,xyzzyaaaa4),'(es24.16)')
tmpr4=r2s2(xyzzyaaah4(xyzzyaaaa4),'(es24.16)')
write(io_mpc,'(4(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3),trim(tmpr4&
&)
enddo
write(io_mpc,'(a)')'END EEPOT DATA'
deallocate(xyzzyaaah4,xyzzyaaai4,xyzzyaaaj4)
end subroutine xyzzyaaac1
subroutine xyzzyaaad1(mpc_cutoff,no_gvec,gvec,gg)
implicit none
integer,intent(in) :: no_gvec
real(dp),intent(in) :: mpc_cutoff,gvec(3,no_gvec)
real(dp),intent(out) :: gg(no_gvec)
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5(3),xyzzyaaad5,xyzzyaaae5(3)
real(dp) xyzzyaaaf5,xyzzyaaag5,xyzzyaaah5,xyzzyaaai5,xyzzyaaaj5,xyzzya&
&aak5,xyzzyaaal5,xyzzyaaam5,xyzzyaaan5,xyzzyaaao5,xyzzyaaap5,xyzzyaaaq&
&5,xyzzyaaar5
real(dp),allocatable :: xyzzyaaas5(:),xyzzyaaat5(:),xyzzyaaau5(:)
xyzzyaaaf5=sqrt(dot_product(a1,a1))
xyzzyaaag5=sqrt(dot_product(a2,a2))
xyzzyaaah5=sqrt(dot_product(a3,a3))
call xyzzyaaah1(xyzzyaaaf5,xyzzyaaag5,xyzzyaaah5,mpc_cutoff,xyzzyaaae5&
&)
allocate(xyzzyaaas5(no_gvec),xyzzyaaat5(no_gvec),xyzzyaaau5(no_gvec),s&
&tat=xyzzyaaab5)
call check_alloc(xyzzyaaab5,'EXTRAP_GG','')
xyzzyaaac5=xyzzyaaae5*4
call wout('MPC generator using '//trim(i2s(xyzzyaaac5(1)))//'x'//trim(&
&i2s(xyzzyaaac5(2)))//'x'//trim(i2s(xyzzyaaac5(3)))//' grid.')
if(xyzzyaaac5(1)*xyzzyaaac5(2)*xyzzyaaac5(3)>2**24)call errwarn('EXTRA&
&P_GG','The MPC grid is extremely large. Consider reducing MPC_CUTOFF.&
&')
xyzzyaaak5=1.d0/dble(xyzzyaaac5(1)**2)
call xyzzyaaae1(xyzzyaaac5,no_gvec,gvec,xyzzyaaau5)
xyzzyaaac5=xyzzyaaac5/2
xyzzyaaaj5=1.d0/dble(xyzzyaaac5(1)**2)
call xyzzyaaae1(xyzzyaaac5,no_gvec,gvec,xyzzyaaat5)
xyzzyaaac5=xyzzyaaac5/2
xyzzyaaai5=1.d0/dble(xyzzyaaac5(1)**2)
call xyzzyaaae1(xyzzyaaac5,no_gvec,gvec,xyzzyaaas5)
xyzzyaaal5=xyzzyaaaj5*xyzzyaaak5/((xyzzyaaai5-xyzzyaaaj5)*(xyzzyaaai5-&
&xyzzyaaak5))
xyzzyaaam5=xyzzyaaai5*xyzzyaaak5/((xyzzyaaaj5-xyzzyaaai5)*(xyzzyaaaj5-&
&xyzzyaaak5))
xyzzyaaan5=xyzzyaaai5*xyzzyaaaj5/((xyzzyaaak5-xyzzyaaai5)*(xyzzyaaak5-&
&xyzzyaaaj5))
xyzzyaaao5=xyzzyaaak5/(xyzzyaaak5-xyzzyaaaj5)
xyzzyaaap5=xyzzyaaaj5/(xyzzyaaaj5-xyzzyaaak5)
xyzzyaaad5=-1
xyzzyaaar5=-1.d0
do xyzzyaaaa5=1,no_gvec
xyzzyaaaq5=xyzzyaaao5*xyzzyaaat5(xyzzyaaaa5)+xyzzyaaap5*xyzzyaaau5(xyz&
&zyaaaa5)
gg(xyzzyaaaa5)=xyzzyaaal5*xyzzyaaas5(xyzzyaaaa5)+xyzzyaaam5*xyzzyaaat5&
&(xyzzyaaaa5)+xyzzyaaan5*xyzzyaaau5(xyzzyaaaa5)
if(abs(gg(xyzzyaaaa5)-xyzzyaaaq5)>xyzzyaaar5)then
xyzzyaaad5=xyzzyaaaa5
xyzzyaaar5=abs(gg(xyzzyaaaa5)-xyzzyaaaq5)
endif
enddo
call wout()
call wout('Analysis of extrapolation of Fourier coefficients of g(r):'&
&)
call wout()
call wout('    1/[FFT grid(1)]^2             g_{G=0}              g_{G&
&=worst case}')
call wout('-----------------------------------------------------------&
&-----------------')
call wout('',(/xyzzyaaai5,xyzzyaaas5(1),xyzzyaaas5(xyzzyaaad5)/),fmt='&
&(a,3(1x,es24.16))')
call wout('',(/xyzzyaaaj5,xyzzyaaat5(1),xyzzyaaat5(xyzzyaaad5)/),fmt='&
&(a,3(1x,es24.16))')
call wout('',(/xyzzyaaak5,xyzzyaaau5(1),xyzzyaaau5(xyzzyaaad5)/),fmt='&
&(a,3(1x,es24.16))')
call wout('      0 (lin. extrap.)   ',(/xyzzyaaao5*xyzzyaaat5(1)+xyzzy&
&aaap5*xyzzyaaau5(1),xyzzyaaao5*xyzzyaaat5(xyzzyaaad5)+xyzzyaaap5*xyzz&
&yaaau5(xyzzyaaad5)/),fmt='(a,2(1x,es24.16))')
call wout('      0 (quad. extrap.)  ',(/gg(1),gg(xyzzyaaad5)/),fmt='(a&
&,2(1x,es24.16))')
call wout()
call wout('Please verify that the extrapolated coefficients look corre&
&ct.')
deallocate(xyzzyaaas5,xyzzyaaat5,xyzzyaaau5)
end subroutine xyzzyaaad1
subroutine xyzzyaaae1(ngrid,no_gvec,gvec,gg)
use slaarnabk, only : ewald_2d
use slaarnabq, only : minimum_image,min_image_brute_force
implicit none
integer,intent(in) :: ngrid(3),no_gvec
real(dp),intent(in) :: gvec(3,no_gvec)
real(dp),intent(out) :: gg(no_gvec)
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6,xyzzyaa&
&af6(3)
real(dp) xyzzyaaag6(3,3),xyzzyaaah6,xyzzyaaai6,xyzzyaaaj6(4),xyzzyaaak&
&6,xyzzyaaal6,xyzzyaaam6(3),xyzzyaaan6(3),xyzzyaaao6(3,3)
complex(dp),allocatable :: xyzzyaaap6(:,:,:)
xyzzyaaak6=-0.5d0/wigner_seitz_radius**3
xyzzyaaal6=1.5d0/wigner_seitz_radius
allocate(xyzzyaaap6(0:ngrid(1)-1,0:ngrid(2)-1,0:ngrid(3)-1),stat=xyzzy&
&aaae6)
call check_alloc(xyzzyaaae6,'FFT_G','FFT grid')
if(periodicity==2)then
xyzzyaaam6(1)=sqrt(dot_product(ainv(1:3,1),ainv(1:3,1)))
xyzzyaaam6(2)=sqrt(dot_product(ainv(1:3,2),ainv(1:3,2)))
xyzzyaaam6(3)=sqrt(dot_product(ainv(1:3,3),ainv(1:3,3)))
xyzzyaaao6=transpose(amat)
endif
xyzzyaaag6(1:3,1)=a1(1:3)/dble(ngrid(1))
xyzzyaaag6(1:3,2)=a2(1:3)/dble(ngrid(2))
xyzzyaaag6(1:3,3)=a3(1:3)/dble(ngrid(3))
do xyzzyaaac6=0,ngrid(3)-1
do xyzzyaaab6=0,ngrid(2)-1
do xyzzyaaaa6=0,ngrid(1)-1
if(periodicity==3)then
xyzzyaaaj6(1:3)=dble(xyzzyaaaa6)*xyzzyaaag6(1:3,1)+dble(xyzzyaaab6)*xy&
&zzyaaag6(1:3,2)+dble(xyzzyaaac6)*xyzzyaaag6(1:3,3)
call minimum_image(4,1,xyzzyaaaj6)
xyzzyaaah6=xyzzyaaaj6(1)**2+xyzzyaaaj6(2)**2+xyzzyaaaj6(3)**2
if(xyzzyaaah6<xyzzyaaaa1)then
xyzzyaaai6=xyzzyaaak6*xyzzyaaah6+xyzzyaaal6
else
xyzzyaaai6=1.d0/sqrt(xyzzyaaah6)
endif
else
if(xyzzyaaaa6>0.or.xyzzyaaab6>0.or.xyzzyaaac6>0)then
xyzzyaaan6=dble(xyzzyaaaa6)*xyzzyaaag6(1:3,1)+dble(xyzzyaaab6)*xyzzyaa&
&ag6(1:3,2)+dble(xyzzyaaac6)*xyzzyaaag6(1:3,3)
call min_image_brute_force(3,xyzzyaaan6,xyzzyaaao6,ainv,xyzzyaaaj6,xyz&
&zyaaam6)
xyzzyaaaj6(4)=sqrt(xyzzyaaaj6(1)**2+xyzzyaaaj6(2)**2+xyzzyaaaj6(3)**2)
call ewald_2d(1,xyzzyaaaj6,xyzzyaaai6)
xyzzyaaai6=xyzzyaaai6-1.d0/xyzzyaaaj6(4)
else
xyzzyaaai6=2.d0*self_term
endif
endif
xyzzyaaap6(xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6)=cmplx(xyzzyaaai6,0.d0,dp)
enddo
enddo
enddo
call timer('DENSITY_FFT',.true.)
call fftn(xyzzyaaap6,ngrid,inv=.true.)
call timer('DENSITY_FFT',.false.)
do xyzzyaaad6=1,no_gvec
xyzzyaaaf6=nint(matmul(gvec(1:3,xyzzyaaad6),binv))
if(xyzzyaaaf6(1)<0)xyzzyaaaf6(1)=xyzzyaaaf6(1)+ngrid(1)
if(xyzzyaaaf6(2)<0)xyzzyaaaf6(2)=xyzzyaaaf6(2)+ngrid(2)
if(xyzzyaaaf6(3)<0)xyzzyaaaf6(3)=xyzzyaaaf6(3)+ngrid(3)
gg(xyzzyaaad6)=dble(xyzzyaaap6(xyzzyaaaf6(1),xyzzyaaaf6(2),xyzzyaaaf6(&
&3)))
enddo
deallocate(xyzzyaaap6)
call dscal(no_gvec,1.d0/sqrt(dble(product(ngrid))),gg(1),1)
end subroutine xyzzyaaae1
integer function xyzzyaaaf1(maglv1,maglv2,maglv3,rlv1,rlv2,rlv3,mpc_cu&
&toff)
implicit none
real(dp),intent(in) :: mpc_cutoff,maglv1,maglv2,maglv3,rlv1(3),rlv2(3)&
&,rlv3(3)
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7(3)
real(dp) xyzzyaaae7(3),xyzzyaaaf7
xyzzyaaaf7=2.d0*mpc_cutoff
xyzzyaaad7(1)=floor(one_over_twopi*sqrt(xyzzyaaaf7)*maglv1)
xyzzyaaad7(2)=floor(one_over_twopi*sqrt(xyzzyaaaf7)*maglv2)
xyzzyaaad7(3)=floor(one_over_twopi*sqrt(xyzzyaaaf7)*maglv3)
xyzzyaaaf1=0
do xyzzyaaac7=-xyzzyaaad7(3),xyzzyaaad7(3)
do xyzzyaaab7=-xyzzyaaad7(2),xyzzyaaad7(2)
do xyzzyaaaa7=-xyzzyaaad7(1),xyzzyaaad7(1)
xyzzyaaae7=dble(xyzzyaaaa7)*rlv1+dble(xyzzyaaab7)*rlv2+dble(xyzzyaaac7&
&)*rlv3
if(dot_product(xyzzyaaae7,xyzzyaaae7)<xyzzyaaaf7)xyzzyaaaf1=xyzzyaaaf1&
&+1
enddo
enddo
enddo
end function xyzzyaaaf1
subroutine xyzzyaaag1(maglv1,maglv2,maglv3,rlv1,rlv2,rlv3,mpc_cutoff,n&
&o_gvec,gvec,gmag)
implicit none
integer,intent(in) :: no_gvec
real(dp),intent(in) :: mpc_cutoff,maglv1,maglv2,maglv3,rlv1(3),rlv2(3)&
&,rlv3(3)
real(dp),intent(out) :: gvec(3,no_gvec),gmag(no_gvec)
integer xyzzyaaaa8,xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8(3),xyzz&
&yaaaf8,xyzzyaaag8
integer,allocatable :: xyzzyaaah8(:)
real(dp) xyzzyaaai8(3),xyzzyaaaj8,xyzzyaaak8
real(dp),allocatable :: xyzzyaaal8(:,:),xyzzyaaam8(:)
xyzzyaaak8=2.d0*mpc_cutoff
xyzzyaaae8(1)=floor(one_over_twopi*sqrt(xyzzyaaak8)*maglv1)
xyzzyaaae8(2)=floor(one_over_twopi*sqrt(xyzzyaaak8)*maglv2)
xyzzyaaae8(3)=floor(one_over_twopi*sqrt(xyzzyaaak8)*maglv3)
xyzzyaaad8=0
do xyzzyaaac8=-xyzzyaaae8(3),xyzzyaaae8(3)
do xyzzyaaab8=-xyzzyaaae8(2),xyzzyaaae8(2)
do xyzzyaaaa8=-xyzzyaaae8(1),xyzzyaaae8(1)
xyzzyaaai8=dble(xyzzyaaaa8)*rlv1+dble(xyzzyaaab8)*rlv2+dble(xyzzyaaac8&
&)*rlv3
xyzzyaaaj8=dot_product(xyzzyaaai8,xyzzyaaai8)
if(xyzzyaaaj8<xyzzyaaak8)then
xyzzyaaad8=xyzzyaaad8+1
gvec(1:3,xyzzyaaad8)=xyzzyaaai8
gmag(xyzzyaaad8)=sqrt(xyzzyaaaj8)
endif
enddo
enddo
enddo
if(xyzzyaaad8/=no_gvec)call errstop('GET_G_VECTORS','Bug: number of G &
&vectors is incorrect.')
allocate(xyzzyaaah8(no_gvec),xyzzyaaam8(no_gvec),xyzzyaaal8(3,no_gvec)&
&,stat=xyzzyaaaf8)
call check_alloc(xyzzyaaaf8,'GET_G_VECTORS','')
call quicksort(no_gvec,gmag(1),xyzzyaaah8(1))
xyzzyaaam8=gmag
xyzzyaaal8=gvec
do xyzzyaaad8=1,no_gvec
xyzzyaaag8=xyzzyaaah8(xyzzyaaad8)
gmag(xyzzyaaad8)=xyzzyaaam8(xyzzyaaag8)
gvec(1:3,xyzzyaaad8)=xyzzyaaal8(1:3,xyzzyaaag8)
enddo
deallocate(xyzzyaaah8,xyzzyaaam8,xyzzyaaal8)
end subroutine xyzzyaaag1
subroutine xyzzyaaah1(maglv1,maglv2,maglv3,mpc_cutoff,ngrid)
implicit none
integer,intent(out) :: ngrid(3)
real(dp),intent(in) :: maglv1,maglv2,maglv3,mpc_cutoff
integer xyzzyaaaa9(3)
real(dp) xyzzyaaab9
integer,parameter :: xyzzyaaac9=16
xyzzyaaab9=1+2*floor(sqrt(2.d0*mpc_cutoff)*maglv1*one_over_twopi)
xyzzyaaaa9(1)=ceiling(log(dble(xyzzyaaab9))/log(2.d0))
xyzzyaaab9=1+2*floor(sqrt(2.d0*mpc_cutoff)*maglv2*one_over_twopi)
xyzzyaaaa9(2)=ceiling(log(dble(xyzzyaaab9))/log(2.d0))
xyzzyaaab9=1+2*floor(sqrt(2.d0*mpc_cutoff)*maglv3*one_over_twopi)
xyzzyaaaa9(3)=ceiling(log(dble(xyzzyaaab9))/log(2.d0))
do while(sum(xyzzyaaaa9)<xyzzyaaac9)
xyzzyaaaa9=xyzzyaaaa9+1
enddo
ngrid=2**xyzzyaaaa9
end subroutine xyzzyaaah1
end module slaarnabf

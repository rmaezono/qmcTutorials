module slaarnacr
use dsp
use store
use format_utils,only : wout,i2s,r2s,wordwrap
use run_control, only : errstop,errstop_master,errwarn,check_alloc
implicit none
private
public check_kinetic,orbital_check,ke_verbose,ke_forgive
logical ke_verbose,ke_forgive
real(dp),parameter :: xyzzyaaaa1=1.d-3,xyzzyaaab1=1.d-5,xyzzyaaac1=1.d&
&-7
character(7),parameter :: xyzzyaaad1(-1:3)=(/'N/A    ','bad    ','poor&
&   ','good   ','optimal'/)
real(dp),parameter :: xyzzyaaae1=0.1d0
real(dp),parameter :: xyzzyaaaf1=.59d0
integer,parameter  :: xyzzyaaag1=20
contains
subroutine orbital_check(rvec_in,rele,rpsi,grad,lap,sderivs)
use parallel
use slaarnacq
use slaarnaag,     only : czero
use slaarnaan, only : ee_distances
use slaarnabg,      only : dimensionality
implicit none
real(dp),intent(inout) :: rvec_in(3),rele(3,netot),grad(3,nemax,real1_&
&complex2,ndet),lap(nemax,real1_complex2,ndet),rpsi(nemax,real1_comple&
&x2,ndet)
real(dp),intent(inout),optional:: sderivs(6,nemax,real1_complex2,ndet)
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaal2&
&,xyzzyaaam2(6),xyzzyaaan2
real(dp) xyzzyaaao2,xyzzyaaap2(3),xyzzyaaaq2(3),xyzzyaaar2(4,netot),xy&
&zzyaaas2,xyzzyaaat2(6),xyzzyaaau2,xyzzyaaav2,xyzzyaaaw2,xyzzyaaax2,xy&
&zzyaaay2,xyzzyaaaz2,xyzzyaaba2(6),xyzzyaabb2(wfdet_norb,real1_complex&
&2),xyzzyaabc2(3,wfdet_norb,real1_complex2),xyzzyaabd2(wfdet_norb,real&
&1_complex2),xyzzyaabe2(6,wfdet_norb,real1_complex2),xyzzyaabf2(3,wfde&
&t_norb,real1_complex2),xyzzyaabg2(wfdet_norb,real1_complex2)
real(dp),allocatable :: xyzzyaabh2(:,:,:,:,:),xyzzyaabi2(:,:,:,:),xyzz&
&yaabj2(:,:,:),xyzzyaabk2(:,:,:,:)
complex(dp) xyzzyaabl2,xyzzyaabm2,xyzzyaabn2(6),xyzzyaabo2,xyzzyaabp2,&
&xyzzyaabq2(6)
logical xyzzyaabr2,xyzzyaabs2,xyzzyaabt2,xyzzyaabu2,xyzzyaabv2,xyzzyaa&
&bw2,xyzzyaabx2,xyzzyaaby2,xyzzyaabz2
character(80) tmpr,tmpr2,tmpr3
real(dp),parameter :: xyzzyaaca2=1.d-2
real(dp),parameter :: xyzzyaacb2=1.d-2
real(dp),parameter :: xyzzyaacc2=1.d-2
real(dp),parameter :: xyzzyaacd2=1.d-9
real(dp),parameter :: xyzzyaace2=1.d-4
integer :: xyzzyaacf2=2
if(am_master)then
if(present(sderivs))then
call wout('Test gradient, Laplacian and second derivatives of orbitals&
&')
call wout('===========================================================&
&')
else
call wout('Test gradient and Laplacian of orbitals')
call wout('=======================================')
endif
call wout('[Turn off by setting "checkwfn : F" in the input file]')
call wout()
endif
xyzzyaaad2=which_spin(1)
if(present(sderivs))then
allocate(xyzzyaabh2(3,dimensionality,nele(xyzzyaaad2),real1_complex2,n&
&det),xyzzyaabi2(3,nemax,real1_complex2,ndet),xyzzyaabj2(nemax,real1_c&
&omplex2,ndet),xyzzyaabk2(6,nele(xyzzyaaad2),real1_complex2,ndet),stat&
&=xyzzyaaae2)
else
allocate(xyzzyaabh2(3,dimensionality,nele(xyzzyaaad2),real1_complex2,n&
&det),xyzzyaabi2(3,nemax,real1_complex2,ndet),xyzzyaabj2(nemax,real1_c&
&omplex2,ndet),stat=xyzzyaaae2)
endif
call check_alloc(xyzzyaaae2,'ORBITAL_CHECK','')
xyzzyaaap2(1:3)=rvec_in(1:3)
xyzzyaaaq2(1:3)=xyzzyaaap2(1:3)
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
if(xyzzyaacf2==1)then
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.true.,xyzzyaabb2,xyzzyaabc2,xyzzyaabd2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
call copy_orb_to_det(3,xyzzyaaad2,wfdet_orbmap,xyzzyaabc2,grad)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabd2,lap)
elseif(xyzzyaacf2==2)then
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabc2,xyzzyaabd2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.false.,.true.,xyzzyaabb2,xyzzyaabc2,xyzzyaabd2,xyzzyaaar2)
call copy_orb_to_det(3,xyzzyaaad2,wfdet_orbmap,xyzzyaabc2,grad)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabd2,lap)
else
call errstop_master('ORBITAL_CHECK','Unknown value of internal paramet&
&er SUB_SELECT.')
endif
if(present(sderivs))then
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.false.,.true.,xyzzyaabb2,xyzzyaabc2,xyzzyaabd2,xyzzyaaar2,xyzzya&
&abe2)
call copy_orb_to_det(6,xyzzyaaad2,wfdet_orbmap,xyzzyaabe2,sderivs)
endif
do xyzzyaaaf2=1,ndet
do xyzzyaaan2=1,real1_complex2
do xyzzyaaag2=1,nele(xyzzyaaad2)
do xyzzyaaaa2=1,dimensionality
xyzzyaabh2(2,xyzzyaaaa2,xyzzyaaag2,xyzzyaaan2,xyzzyaaaf2) =rpsi(xyzzya&
&aag2,xyzzyaaan2,xyzzyaaaf2)
enddo
enddo
enddo
enddo
xyzzyaaao2=xyzzyaaae1
xyzzyaaac2=0
xyzzyaabv2=.false.
xyzzyaabu2=.false.
xyzzyaabw2=.false.
xyzzyaabr2=.false.
xyzzyaabs2=.false.
xyzzyaabt2=.false.
xyzzyaaak2=0
xyzzyaaal2=0
do
do xyzzyaaaa2=1,dimensionality
xyzzyaaaq2(1:3)=xyzzyaaap2(1:3)
xyzzyaaaq2(xyzzyaaaa2)=xyzzyaaap2(xyzzyaaaa2)-xyzzyaaao2
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabf2,xyzzyaabg2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
xyzzyaabh2(1,xyzzyaaaa2,1:nele(xyzzyaaad2),1:real1_complex2,1:ndet) =r&
&psi(1:nele(xyzzyaaad2),1:real1_complex2,1:ndet)
xyzzyaaaq2(xyzzyaaaa2)=xyzzyaaap2(xyzzyaaaa2)+xyzzyaaao2
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabf2,xyzzyaabg2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
xyzzyaabh2(3,xyzzyaaaa2,1:nele(xyzzyaaad2),1:real1_complex2,1:ndet)=rp&
&si(1:nele(xyzzyaaad2),1:real1_complex2,1:ndet)
enddo
if(present(sderivs).and.dimensionality>1)then
xyzzyaaaq2(1)=xyzzyaaap2(1)-xyzzyaaao2
xyzzyaaaq2(2)=xyzzyaaap2(2)-xyzzyaaao2
xyzzyaaaq2(3)=xyzzyaaap2(3)
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabf2,xyzzyaabg2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
xyzzyaabk2(1,1:nele(xyzzyaaad2),1:real1_complex2,1:ndet) =rpsi(1:nele(&
&xyzzyaaad2),1:real1_complex2,1:ndet)
xyzzyaaaq2(1)=xyzzyaaap2(1)+xyzzyaaao2
xyzzyaaaq2(2)=xyzzyaaap2(2)+xyzzyaaao2
xyzzyaaaq2(3)=xyzzyaaap2(3)
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabf2,xyzzyaabg2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
xyzzyaabk2(2,1:nele(xyzzyaaad2),1:real1_complex2,1:ndet) =rpsi(1:nele(&
&xyzzyaaad2),1:real1_complex2,1:ndet)
if(dimensionality==3)then
xyzzyaaaq2(1)=xyzzyaaap2(1)
xyzzyaaaq2(2)=xyzzyaaap2(2)-xyzzyaaao2
xyzzyaaaq2(3)=xyzzyaaap2(3)-xyzzyaaao2
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabf2,xyzzyaabg2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
xyzzyaabk2(3,1:nele(xyzzyaaad2),1:real1_complex2,1:ndet)=rpsi(1:nele(x&
&yzzyaaad2),1:real1_complex2,1:ndet)
xyzzyaaaq2(1)=xyzzyaaap2(1)
xyzzyaaaq2(2)=xyzzyaaap2(2)+xyzzyaaao2
xyzzyaaaq2(3)=xyzzyaaap2(3)+xyzzyaaao2
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabf2,xyzzyaabg2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
xyzzyaabk2(4,1:nele(xyzzyaaad2),1:real1_complex2,1:ndet)=rpsi(1:nele(x&
&yzzyaaad2),1:real1_complex2,1:ndet)
xyzzyaaaq2(1)=xyzzyaaap2(1)-xyzzyaaao2
xyzzyaaaq2(2)=xyzzyaaap2(2)
xyzzyaaaq2(3)=xyzzyaaap2(3)-xyzzyaaao2
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabf2,xyzzyaabg2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
xyzzyaabk2(5,1:nele(xyzzyaaad2),1:real1_complex2,1:ndet)=rpsi(1:nele(x&
&yzzyaaad2),1:real1_complex2,1:ndet)
xyzzyaaaq2(1)=xyzzyaaap2(1)+xyzzyaaao2
xyzzyaaaq2(2)=xyzzyaaap2(2)
xyzzyaaaq2(3)=xyzzyaaap2(3)+xyzzyaaao2
rpsi=666.661d0
if(pairing_wf)call ee_distances(netot,xyzzyaaaq2,rele,xyzzyaaar2)
call wfdet(xyzzyaaaq2,1,xyzzyaaad2,wfdet_norb,wfdet_orbmask(1,xyzzyaaa&
&d2),.true.,.false.,xyzzyaabb2,xyzzyaabf2,xyzzyaabg2,xyzzyaaar2)
call copy_orb_to_det(1,xyzzyaaad2,wfdet_orbmap,xyzzyaabb2,rpsi)
xyzzyaabk2(6,1:nele(xyzzyaaad2),1:real1_complex2,1:ndet)=rpsi(1:nele(x&
&yzzyaaad2),1:real1_complex2,1:ndet)
endif
endif
if(xyzzyaaac2==0)then
do xyzzyaaaf2=1,ndet
do xyzzyaaag2=1,nele(xyzzyaaad2)
if(all(abs(xyzzyaabh2(1:3,1:dimensionality,xyzzyaaag2,1:real1_complex2&
&,xyzzyaaaf2))<xyzzyaacd2))then
tmpr=r2s(xyzzyaacd2,'(g18.1)')
call wout('Orbital '//trim(i2s(xyzzyaaag2))//', spin '//trim(i2s(xyzzy&
&aaad2))//', det '//trim(i2s(xyzzyaaaf2))//' appears to be zero (< '//&
&trim(tmpr)//').')
endif
enddo
enddo
endif
if(.not.xyzzyaabu2)then
xyzzyaabr2=.false.
xyzzyaaau2=0.d0
xyzzyaaak2=0
xyzzyaaay2=0.d0
endif
if(.not.xyzzyaabv2)then
xyzzyaabs2=.false.
xyzzyaaas2=0.d0
xyzzyaaal2=0
xyzzyaaaz2=0.d0
endif
if(present(sderivs))then
if(.not.xyzzyaabw2)then
xyzzyaabt2=.false.
xyzzyaaat2(1:6)=0.d0
xyzzyaaam2(1:6)=0
xyzzyaaba2(1:6)=0.d0
endif
endif
do xyzzyaaaf2=1,ndet
do xyzzyaaag2=1,nele(xyzzyaaad2)
xyzzyaabx2=.false.
xyzzyaaby2=.false.
xyzzyaabz2=.false.
if(.not.xyzzyaabu2)then
do xyzzyaaaa2=1,dimensionality
if(complex_wf)then
xyzzyaabl2=cmplx(xyzzyaabh2(3,xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2)-xyzz&
&yaabh2(1,xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2),xyzzyaabh2(3,xyzzyaaaa2,&
&xyzzyaaag2,2,xyzzyaaaf2)-xyzzyaabh2(1,xyzzyaaaa2,xyzzyaaag2,2,xyzzyaa&
&af2),dp)/(2.d0*xyzzyaaao2)
xyzzyaabo2=cmplx(grad(xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2),grad(xyzzyaa&
&aa2,xyzzyaaag2,2,xyzzyaaaf2),dp)
else
xyzzyaabl2=cmplx((xyzzyaabh2(3,xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2) -xy&
&zzyaabh2(1,xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2))/(2.d0*xyzzyaaao2),0.d&
&0,dp)
xyzzyaabo2=cmplx(grad(xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2),0.d0,dp)
endif
if(abs(xyzzyaabo2)<xyzzyaacd2)then
if(abs(xyzzyaabl2-xyzzyaabo2)<xyzzyaace2)then
xyzzyaabx2=.false.
else
xyzzyaabr2=.true.
xyzzyaabx2=.true.
endif
else
if(abs(xyzzyaabl2-xyzzyaabo2)<xyzzyaaca2*abs(xyzzyaabo2))then
xyzzyaabx2=.false.
else
xyzzyaabr2=.true.
xyzzyaabx2=.true.
endif
xyzzyaaak2=xyzzyaaak2+1
xyzzyaaay2=xyzzyaaay2+abs((xyzzyaabl2-xyzzyaabo2)/xyzzyaabo2)
endif
if(xyzzyaabx2)then
if(xyzzyaaac2>=xyzzyaaag1.or.ke_verbose)then
call wout('Failed gradient test.')
call wout('Spin                        : '//trim(i2s(xyzzyaaad2)))
call wout('Orbital                     : '//trim(i2s(xyzzyaaag2)))
call wout('Determinant                 : '//trim(i2s(xyzzyaaaf2)))
tmpr=r2s(dble(xyzzyaaap2(1)),'(es12.4)')
tmpr2=r2s(dble(xyzzyaaap2(2)),'(es12.4)')
tmpr3=r2s(dble(xyzzyaaap2(3)),'(es12.4)')
call wout('Position                    : ('//trim(tmpr)//', '//trim(tm&
&pr2)//', '//trim(tmpr3)//')')
call wout('Cartesian component         : '//trim(i2s(xyzzyaaaa2)))
tmpr=r2s(dble(xyzzyaabo2),'(f16.8)')
call wout('Analytic gradient (real)    : '//trim(tmpr))
tmpr=r2s(dble(xyzzyaabl2),'(f16.8)')
call wout('Numerical gradient (real)   : '//trim(tmpr))
if(complex_wf)then
tmpr=r2s(aimag(xyzzyaabo2),'(f16.8)')
call wout('Analytic gradient (imag)    : '//trim(tmpr))
tmpr=r2s(aimag(xyzzyaabl2),'(f16.8)')
call wout('Numerical gradient (imag)   : '//trim(tmpr))
endif
tmpr=r2s(xyzzyaaca2,'(f16.8)')
call wout('Fract. difference tolerance : '//trim(tmpr))
tmpr=r2s(xyzzyaaae1,'(f16.8)')
call wout('Initial stepsize            : '//trim(tmpr))
tmpr=r2s(xyzzyaaao2,'(f16.8)')
call wout('Current stepsize            : '//trim(tmpr))
call wout()
if(xyzzyaaac2>=xyzzyaaag1)then
if(nnodes>1)then
call errstop('ORBITAL_CHECK','Failed gradient test on node '//trim(i2s&
&(my_node))//'.')
else
call errstop('ORBITAL_CHECK','Failed gradient test.')
endif
endif
endif
endif
xyzzyaaau2=xyzzyaaau2+abs(xyzzyaabo2-xyzzyaabl2)
enddo
endif
if(.not.xyzzyaabv2)then
xyzzyaabm2=czero
if(complex_wf)then
do xyzzyaaaa2=1,dimensionality
xyzzyaabm2=xyzzyaabm2+cmplx(xyzzyaabh2(3,xyzzyaaaa2,xyzzyaaag2,1,xyzzy&
&aaaf2)+xyzzyaabh2(1,xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaab&
&h2(2,xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2),xyzzyaabh2(3,xyzzyaaaa2,xyzz&
&yaaag2,2,xyzzyaaaf2)+xyzzyaabh2(1,xyzzyaaaa2,xyzzyaaag2,2,xyzzyaaaf2)&
&-2.d0*xyzzyaabh2(2,xyzzyaaaa2,xyzzyaaag2,2,xyzzyaaaf2),dp)/(xyzzyaaao&
&2*xyzzyaaao2)
enddo
xyzzyaabp2=cmplx(lap(xyzzyaaag2,1,xyzzyaaaf2),lap(xyzzyaaag2,2,xyzzyaa&
&af2),dp)
else
do xyzzyaaaa2=1,dimensionality
xyzzyaabm2=xyzzyaabm2+cmplx((xyzzyaabh2(3,xyzzyaaaa2,xyzzyaaag2,1,xyzz&
&yaaaf2)+xyzzyaabh2(1,xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaa&
&bh2(2,xyzzyaaaa2,xyzzyaaag2,1,xyzzyaaaf2))/(xyzzyaaao2*xyzzyaaao2),0.&
&d0,dp)
enddo
xyzzyaabp2=cmplx(lap(xyzzyaaag2,1,xyzzyaaaf2),0.d0,dp)
endif
if(abs(xyzzyaabp2)<xyzzyaacd2)then
if(abs(xyzzyaabm2-xyzzyaabp2)<xyzzyaace2)then
xyzzyaaby2=.false.
else
xyzzyaabs2=.true.
xyzzyaaby2=.true.
endif
else
if(abs(xyzzyaabm2-xyzzyaabp2)<xyzzyaacb2*abs(xyzzyaabp2))then
xyzzyaaby2=.false.
else
xyzzyaabs2=.true.
xyzzyaaby2=.true.
endif
xyzzyaaal2=xyzzyaaal2+1
xyzzyaaaz2=xyzzyaaaz2+abs((xyzzyaabm2-xyzzyaabp2)/xyzzyaabp2)
endif
if(xyzzyaaby2)then
if(xyzzyaaac2>=xyzzyaaag1.or.ke_verbose)then
call wout('Failed Laplacian test.')
call wout('Spin                        : '//trim(i2s(xyzzyaaad2)))
call wout('Orbital                     : '//trim(i2s(xyzzyaaag2)))
call wout('Determinant                 : '//trim(i2s(xyzzyaaaf2)))
tmpr=r2s(dble(xyzzyaaap2(1)),'(es12.4)')
tmpr2=r2s(dble(xyzzyaaap2(2)),'(es12.4)')
tmpr3=r2s(dble(xyzzyaaap2(3)),'(es12.4)')
call wout('Position                    : ('//trim(tmpr)//', '//trim(tm&
&pr2)//', '//trim(tmpr3)//')')
tmpr=r2s(dble(xyzzyaabp2),'(f16.8)')
call wout('Analytic Laplacian (real)   : '//trim(tmpr))
tmpr=r2s(dble(xyzzyaabm2),'(f16.8)')
call wout('Numerical Laplacian (real)  : '//trim(tmpr))
if(complex_wf)then
tmpr=r2s(aimag(xyzzyaabp2),'(f16.8)')
call wout('Analytic Laplacian (imag)   : '//trim(tmpr))
tmpr=r2s(aimag(xyzzyaabm2),'(f16.8)')
call wout('Numerical Laplacian (imag)  : '//trim(tmpr))
endif
tmpr=r2s(xyzzyaacb2,'(f16.8)')
call wout('Fract. difference tolerance : '//trim(tmpr))
tmpr=r2s(xyzzyaaae1,'(f16.8)')
call wout('Initial stepsize            : '//trim(tmpr))
tmpr=r2s(xyzzyaaao2,'(f16.8)')
call wout('Current stepsize            : '//trim(tmpr))
call wout()
if(xyzzyaaac2>=xyzzyaaag1)then
if(nnodes>1)then
call errstop('ORBITAL_CHECK','Failed Laplacian test on node ' //trim(i&
&2s(my_node))//'.')
else
call errstop('ORBITAL_CHECK','Failed Laplacian test.')
endif
endif
endif
endif
xyzzyaaas2=xyzzyaaas2+abs(xyzzyaabm2-xyzzyaabp2)
endif
if(present(sderivs))then
if(.not.xyzzyaabw2)then
xyzzyaabn2=czero
if(complex_wf)then
do xyzzyaaab2=1,dimensionality
xyzzyaabn2(xyzzyaaab2)=cmplx(xyzzyaabh2(3,xyzzyaaab2,xyzzyaaag2,1,xyzz&
&yaaaf2)+xyzzyaabh2(1,xyzzyaaab2,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaa&
&bh2(2,xyzzyaaab2,xyzzyaaag2,1,xyzzyaaaf2),xyzzyaabh2(3,xyzzyaaab2,xyz&
&zyaaag2,2,xyzzyaaaf2)+xyzzyaabh2(1,xyzzyaaab2,xyzzyaaag2,2,xyzzyaaaf2&
&)-2.d0*xyzzyaabh2(2,xyzzyaaab2,xyzzyaaag2,2,xyzzyaaaf2),dp)/(xyzzyaaa&
&o2*xyzzyaaao2)
enddo
if(dimensionality>1)then
xyzzyaabn2(4)=cmplx((xyzzyaabk2(1,xyzzyaaag2,1,xyzzyaaaf2)+xyzzyaabk2(&
&2,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,2,xyzzyaaag2,1,xyzzyaaaf&
&2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sderivs(1,xyzzyaaag2,1,xyzzy&
&aaaf2)+sderivs(2,xyzzyaaag2,1,xyzzyaaaf2)),(xyzzyaabk2(1,xyzzyaaag2,2&
&,xyzzyaaaf2)+xyzzyaabk2(2,xyzzyaaag2,2,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,&
&2,xyzzyaaag2,2,xyzzyaaaf2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sder&
&ivs(1,xyzzyaaag2,2,xyzzyaaaf2)+sderivs(2,xyzzyaaag2,2,xyzzyaaaf2)),dp&
&)
if(dimensionality==3)then
xyzzyaabn2(5)=cmplx((xyzzyaabk2(5,xyzzyaaag2,1,xyzzyaaaf2)+xyzzyaabk2(&
&6,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,2,xyzzyaaag2,1,xyzzyaaaf&
&2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sderivs(1,xyzzyaaag2,1,xyzzy&
&aaaf2)+sderivs(3,xyzzyaaag2,1,xyzzyaaaf2)),(xyzzyaabk2(5,xyzzyaaag2,2&
&,xyzzyaaaf2)+xyzzyaabk2(6,xyzzyaaag2,2,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,&
&2,xyzzyaaag2,2,xyzzyaaaf2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sder&
&ivs(1,xyzzyaaag2,2,xyzzyaaaf2)+sderivs(3,xyzzyaaag2,2,xyzzyaaaf2)),dp&
&)
xyzzyaabn2(6)=cmplx((xyzzyaabk2(3,xyzzyaaag2,1,xyzzyaaaf2)+xyzzyaabk2(&
&4,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,2,xyzzyaaag2,1,xyzzyaaaf&
&2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sderivs(2,xyzzyaaag2,1,xyzzy&
&aaaf2)+sderivs(3,xyzzyaaag2,1,xyzzyaaaf2)),(xyzzyaabk2(3,xyzzyaaag2,2&
&,xyzzyaaaf2)+xyzzyaabk2(4,xyzzyaaag2,2,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,&
&2,xyzzyaaag2,2,xyzzyaaaf2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sder&
&ivs(2,xyzzyaaag2,2,xyzzyaaaf2)+sderivs(3,xyzzyaaag2,2,xyzzyaaaf2)),dp&
&)
endif
endif
do xyzzyaaab2=1,6
xyzzyaabq2(xyzzyaaab2)=cmplx(sderivs(xyzzyaaab2,xyzzyaaag2,1,xyzzyaaaf&
&2),sderivs(xyzzyaaab2,xyzzyaaag2,2,xyzzyaaaf2),dp)
enddo
else
do xyzzyaaab2=1,dimensionality
xyzzyaabn2(xyzzyaaab2)=cmplx((xyzzyaabh2(3,xyzzyaaab2,xyzzyaaag2,1,xyz&
&zyaaaf2)+xyzzyaabh2(1,xyzzyaaab2,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzya&
&abh2(2,xyzzyaaab2,xyzzyaaag2,1,xyzzyaaaf2))/(xyzzyaaao2*xyzzyaaao2),0&
&.d0,dp)
enddo
if(dimensionality>1)then
xyzzyaabn2(4)=cmplx((xyzzyaabk2(1,xyzzyaaag2,1,xyzzyaaaf2)+xyzzyaabk2(&
&2,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,2,xyzzyaaag2,1,xyzzyaaaf&
&2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sderivs(1,xyzzyaaag2,1,xyzzy&
&aaaf2)+sderivs(2,xyzzyaaag2,1,xyzzyaaaf2)),0.d0,dp)
if(dimensionality==3)then
xyzzyaabn2(5)=cmplx((xyzzyaabk2(5,xyzzyaaag2,1,xyzzyaaaf2)+xyzzyaabk2(&
&6,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,2,xyzzyaaag2,1,xyzzyaaaf&
&2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sderivs(1,xyzzyaaag2,1,xyzzy&
&aaaf2)+sderivs(3,xyzzyaaag2,1,xyzzyaaaf2)),0.d0,dp)
xyzzyaabn2(6)=cmplx((xyzzyaabk2(3,xyzzyaaag2,1,xyzzyaaaf2)+xyzzyaabk2(&
&4,xyzzyaaag2,1,xyzzyaaaf2)-2.d0*xyzzyaabh2(2,2,xyzzyaaag2,1,xyzzyaaaf&
&2))*0.5d0/(xyzzyaaao2*xyzzyaaao2)-0.5d0*(sderivs(2,xyzzyaaag2,1,xyzzy&
&aaaf2)+sderivs(3,xyzzyaaag2,1,xyzzyaaaf2)),0.d0,dp)
endif
endif
do xyzzyaaab2=1,6
xyzzyaabq2(xyzzyaaab2)=cmplx(sderivs(xyzzyaaab2,xyzzyaaag2,1,xyzzyaaaf&
&2),0.d0,dp)
enddo
endif
xyzzyaabz2=.false.
do xyzzyaaab2=1,6
if((dimensionality==2.and.(xyzzyaaab2==3.or.xyzzyaaab2>4)).or.(dimensi&
&onality==1.and.xyzzyaaab2/=1))cycle
if(abs(xyzzyaabq2(xyzzyaaab2))<xyzzyaacd2)then
if(.not.abs(xyzzyaabn2(xyzzyaaab2)-xyzzyaabq2(xyzzyaaab2))<xyzzyaace2)&
&then
xyzzyaabt2=.true.
xyzzyaabz2=.true.
endif
else
if(.not.abs(xyzzyaabn2(xyzzyaaab2)-xyzzyaabq2(xyzzyaaab2))<xyzzyaacb2*&
&abs(xyzzyaabq2(xyzzyaaab2)))then
xyzzyaabt2=.true.
xyzzyaabz2=.true.
endif
xyzzyaaam2(xyzzyaaab2)=xyzzyaaam2(xyzzyaaab2)+1
xyzzyaaba2(xyzzyaaab2)=xyzzyaaba2(xyzzyaaab2)+abs((xyzzyaabn2(xyzzyaaa&
&b2)-xyzzyaabq2(xyzzyaaab2))/xyzzyaabq2(xyzzyaaab2))
endif
enddo
if(xyzzyaabz2)then
if(xyzzyaaac2>=xyzzyaaag1.or.ke_verbose)then
call wout('Failed second derivatives test.')
call wout('Spin                        : '//trim(i2s(xyzzyaaad2)))
call wout('Orbital                     : '//trim(i2s(xyzzyaaag2)))
call wout('Determinant                 : '//trim(i2s(xyzzyaaaf2)))
tmpr=r2s(dble(xyzzyaaap2(1)),'(es12.4)')
tmpr2=r2s(dble(xyzzyaaap2(2)),'(es12.4)')
tmpr3=r2s(dble(xyzzyaaap2(3)),'(es12.4)')
call wout('Position                    : ('//trim(tmpr)//', '//trim(tm&
&pr2)//', '//trim(tmpr3)//')')
call wout('Analytic derivatives        : ')
tmpr=r2s(dble(xyzzyaabq2(1)),'(f16.8)')
call wout(' d^2 Phi/dx dx  (real)        '//trim(tmpr))
if(dimensionality>1)then
tmpr=r2s(dble(xyzzyaabq2(2)),'(f16.8)')
call wout(' d^2 Phi/dy dy  (real)        '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(dble(xyzzyaabq2(3)),'(f16.8)')
call wout(' d^2 Phi/dz dz  (real)        '//trim(tmpr))
endif
tmpr=r2s(dble(xyzzyaabq2(4)),'(f16.8)')
call wout(' d^2 Phi/dx dy  (real)        '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(dble(xyzzyaabq2(5)),'(f16.8)')
call wout(' d^2 Phi/dx dz  (real)        '//trim(tmpr))
tmpr=r2s(dble(xyzzyaabq2(6)),'(f16.8)')
call wout(' d^2 Phi/dy dz  (real)        '//trim(tmpr))
endif
endif
call wout('Numerical derivatives       : ')
tmpr=r2s(dble(xyzzyaabn2(1)),'(f16.8)')
call wout(' d^2 Phi/dx dx  (real)        '//trim(tmpr))
if(dimensionality>1)then
tmpr=r2s(dble(xyzzyaabn2(2)),'(f16.8)')
call wout(' d^2 Phi/dy dy  (real)        '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(dble(xyzzyaabn2(3)),'(f16.8)')
call wout(' d^2 Phi/dz dz  (real)        '//trim(tmpr))
endif
tmpr=r2s(dble(xyzzyaabn2(4)),'(f16.8)')
call wout(' d^2 Phi/dx dy  (real)        '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(dble(xyzzyaabn2(5)),'(f16.8)')
call wout(' d^2 Phi/dx dz  (real)        '//trim(tmpr))
tmpr=r2s(dble(xyzzyaabn2(6)),'(f16.8)')
call wout(' d^2 Phi/dy dz  (real)        '//trim(tmpr))
endif
endif
if(complex_wf)then
tmpr=r2s(aimag(xyzzyaabq2(1)),'(f16.8)')
call wout(' d^2 Phi/dx dx  (imag)        '//trim(tmpr))
if(dimensionality>1)then
tmpr=r2s(aimag(xyzzyaabq2(2)),'(f16.8)')
call wout(' d^2 Phi/dy dy  (imag)        '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(aimag(xyzzyaabq2(3)),'(f16.8)')
call wout(' d^2 Phi/dz dz  (imag)        '//trim(tmpr))
endif
tmpr=r2s(aimag(xyzzyaabq2(4)),'(f16.8)')
call wout(' d^2 Phi/dx dy  (imag)        '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(aimag(xyzzyaabq2(5)),'(f16.8)')
call wout(' d^2 Phi/dx dz  (imag)        '//trim(tmpr))
tmpr=r2s(aimag(xyzzyaabq2(6)),'(f16.8)')
call wout(' d^2 Phi/dy dz  (imag)        '//trim(tmpr))
endif
endif
call wout('Numerical derivatives       : ')
tmpr=r2s(aimag(xyzzyaabn2(1)),'(f16.8)')
call wout(' d^2 Phi/dx dx  (imag)        '//trim(tmpr))
if(dimensionality>1)then
tmpr=r2s(aimag(xyzzyaabn2(2)),'(f16.8)')
call wout(' d^2 Phi/dy dy  (imag)        '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(aimag(xyzzyaabn2(3)),'(f16.8)')
call wout(' d^2 Phi/dz dz  (imag)        '//trim(tmpr))
endif
tmpr=r2s(aimag(xyzzyaabn2(4)),'(f16.8)')
call wout(' d^2 Phi/dx dy  (imag)        '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(aimag(xyzzyaabn2(5)),'(f16.8)')
call wout(' d^2 Phi/dx dz  (imag)        '//trim(tmpr))
tmpr=r2s(aimag(xyzzyaabn2(6)),'(f16.8)')
call wout(' d^2 Phi/dy dz  (imag)        '//trim(tmpr))
endif
endif
endif
tmpr=r2s(xyzzyaacc2,'(f16.8)')
call wout('Fract. difference tolerance : '//trim(tmpr))
tmpr=r2s(xyzzyaaae1,'(f16.8)')
call wout('Initial stepsize            : '//trim(tmpr))
tmpr=r2s(xyzzyaaao2,'(f16.8)')
call wout('Current stepsize            : '//trim(tmpr))
call wout()
if(xyzzyaaac2>=xyzzyaaag1)then
if(nnodes>1)then
call errstop('ORBITAL_CHECK','Failed second derivatives test on node '&
&//trim(i2s(my_node))//'.')
else
call errstop('ORBITAL_CHECK','Failed second derivatives test.')
endif
endif
endif
endif
do xyzzyaaab2=1,6
if((dimensionality==2.and.(xyzzyaaab2==3.or.xyzzyaaab2>4)) .or.(dimens&
&ionality==1.and.xyzzyaaab2/=1))cycle
xyzzyaaat2(xyzzyaaab2)=xyzzyaaat2(xyzzyaaab2)+abs(xyzzyaabn2(xyzzyaaab&
&2)-xyzzyaabq2(xyzzyaaab2))
enddo
endif
endif
enddo
enddo
if(.not.xyzzyaabu2)then
if(.not.xyzzyaabr2)then
xyzzyaabu2=.true.
xyzzyaaaw2=xyzzyaaao2
xyzzyaaah2=xyzzyaaac2
xyzzyaaau2=xyzzyaaau2/dble(dimensionality*nele(xyzzyaaad2)*ndet)
if(ke_verbose)then
call wout('Passed gradient.')
call wout()
endif
if(xyzzyaaak2>0)then
xyzzyaaay2=xyzzyaaay2/dble(xyzzyaaak2)
else
xyzzyaaay2=0.d0
endif
endif
endif
if(.not.xyzzyaabs2.and..not.xyzzyaabv2)then
xyzzyaabv2=.true.
xyzzyaaav2=xyzzyaaao2
xyzzyaaai2=xyzzyaaac2
xyzzyaaas2=xyzzyaaas2/dble(nele(xyzzyaaad2)*ndet)
if(ke_verbose)then
call wout('Passed Laplacian.')
call wout()
endif
if(xyzzyaaal2>0)then
xyzzyaaaz2=xyzzyaaaz2/dble(xyzzyaaal2)
else
xyzzyaaaz2=0.d0
endif
endif
if(present(sderivs))then
if(.not.xyzzyaabt2.and..not.xyzzyaabw2)then
xyzzyaabw2=.true.
xyzzyaaax2=xyzzyaaao2
xyzzyaaaj2=xyzzyaaac2
do xyzzyaaab2=1,6
if((dimensionality==2.and.(xyzzyaaab2==3.or.xyzzyaaab2>4)) .or.(dimens&
&ionality==1.and.xyzzyaaab2/=1))cycle
xyzzyaaat2(xyzzyaaab2)=xyzzyaaat2(xyzzyaaab2)/dble(nele(xyzzyaaad2)*nd&
&et)
enddo
if(ke_verbose)then
call wout('Passed second derivatives test.')
call wout()
endif
do xyzzyaaab2=1,6
if((dimensionality==2.and.(xyzzyaaab2==3.or.xyzzyaaab2>4)) .or.(dimens&
&ionality==1.and.xyzzyaaab2/=1))cycle
if(xyzzyaaam2(xyzzyaaab2)>0)then
xyzzyaaba2(xyzzyaaab2)=xyzzyaaba2(xyzzyaaab2)/dble(xyzzyaaal2)
else
xyzzyaaba2(xyzzyaaab2)=0.d0
endif
enddo
endif
endif
if(present(sderivs))then
if(xyzzyaabv2.and.xyzzyaabu2.and.xyzzyaabw2)exit
else
if(xyzzyaabv2.and.xyzzyaabu2)exit
endif
xyzzyaaao2=xyzzyaaao2*xyzzyaaaf1
xyzzyaaac2=xyzzyaaac2+1
enddo
if(my_node<12)then
call wout()
call wout('The orbitals have passed the numerical derivatives test.')
call wout()
endif
if(am_master.and.nnodes>1)then
if(nnodes<=12)then
call wordwrap('(The test results on slave nodes are being written to t&
&heir respective output files, which will be appended to the main outp&
&ut file at the end of the run by the "runqmc" script.)')
else
call wordwrap('(The test results on the first twelve slave nodes are b&
&eing written to their respective output files, which will be appended&
& to the main output file at the end of the run by the "runqmc" script&
&.)')
endif
call wout()
endif
if(my_node<12)then
tmpr=r2s(xyzzyaaae1,'(f16.8)')
call wout('Initial stepsize for gradients      : '//trim(tmpr))
tmpr=r2s(xyzzyaaaw2,'(f16.8)')
call wout('Final stepsize for gradients        : '//trim(tmpr))
call wout('Number of adjustments to stepsize   : '//trim(i2s(xyzzyaaah&
&2)))
tmpr=r2s(xyzzyaaca2,'(f16.8)')
call wout('Fractional tolerance for gradients  : '//trim(tmpr))
tmpr=r2s(xyzzyaaau2,'(f16.8)')
call wout('Average abs. diff. of gradients     : '//trim(tmpr))
tmpr=r2s(xyzzyaaay2,'(f16.8)')
call wout('Average abs. fractional difference  : '//trim(tmpr))
call wout('Number of zeroes encountered        : '//trim(i2s(dimension&
&ality*nele(xyzzyaaad2)*ndet-xyzzyaaak2)))
call wout()
tmpr=r2s(xyzzyaaae1,'(f16.8)')
call wout('Initial stepsize for Laplacians     : '//trim(tmpr))
tmpr=r2s(xyzzyaaav2,'(f16.8)')
call wout('Final stepsize for Laplacians       : '//trim(tmpr))
call wout('Number of adjustments to stepsize   : '//trim(i2s(xyzzyaaai&
&2)))
tmpr=r2s(xyzzyaacb2,'(f16.8)')
call wout('Fractional tolerance for Laplacians : '//trim(tmpr))
tmpr=r2s(xyzzyaaas2,'(f16.8)')
call wout('Average abs. diff. of Laplacians    : '//trim(tmpr))
tmpr=r2s(xyzzyaaaz2,'(f16.8)')
call wout('Average abs. fractional difference  : '//trim(tmpr))
call wout('Number of zeroes encountered        : '//trim(i2s(nele(xyzz&
&yaaad2)*ndet-xyzzyaaal2)))
call wout()
if(present(sderivs))then
tmpr=r2s(xyzzyaaae1,'(f16.8)')
call wout('Initial stepsize for 2nd derivatives: '//trim(tmpr))
tmpr=r2s(xyzzyaaax2,'(f16.8)')
call wout('Final stepsize for derivatives was  : '//trim(tmpr))
call wout('Number of adjustments to stepsize   : '//trim(i2s(xyzzyaaaj&
&2)))
tmpr=r2s(xyzzyaacc2,'(f16.8)')
call wout('Frac. tolerance for 2nd derivatives : '//trim(tmpr))
call wout('Average abs. diff. of -')
tmpr=r2s(xyzzyaaat2(1),'(f16.8)')
call wout(' d^2 Phi/dx dx                      : '//trim(tmpr))
if(dimensionality>1)then
tmpr=r2s(xyzzyaaat2(2),'(f16.8)')
call wout(' d^2 Phi/dy dy                      : '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(xyzzyaaat2(3),'(f16.8)')
call wout(' d^2 Phi/dz dz                      : '//trim(tmpr))
endif
tmpr=r2s(xyzzyaaat2(4),'(f16.8)')
call wout(' d^2 Phi/dx dy                      : '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(xyzzyaaat2(5),'(f16.8)')
call wout(' d^2 Phi/dx dz                      : '//trim(tmpr))
tmpr=r2s(xyzzyaaat2(6),'(f16.8)')
call wout(' d^2 Phi/dy dz                      : '//trim(tmpr))
endif
endif
call wout('Average abs. fractional difference -  ')
tmpr=r2s(xyzzyaaba2(1),'(f16.8)')
call wout(' d^2 Phi/dx dx                      : '//trim(tmpr))
if(dimensionality>1)then
tmpr=r2s(xyzzyaaba2(2),'(f16.8)')
call wout(' d^2 Phi/dy dy                      : '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(xyzzyaaba2(3),'(f16.8)')
call wout(' d^2 Phi/dz dz                      : '//trim(tmpr))
endif
tmpr=r2s(xyzzyaaba2(4),'(f16.8)')
call wout(' d^2 Phi/dx dy                      : '//trim(tmpr))
if(dimensionality==3)then
tmpr=r2s(xyzzyaaba2(5),'(f16.8)')
call wout(' d^2 Phi/dx dz                      : '//trim(tmpr))
tmpr=r2s(xyzzyaaba2(6),'(f16.8)')
call wout(' d^2 Phi/dy dz                      : '//trim(tmpr))
endif
endif
call wout('Number of zeroes encountered - ')
call wout(' d^2 Phi/dx dx                      : '//trim(i2s(nele(xyzz&
&yaaad2)*ndet-xyzzyaaam2(1))))
if(dimensionality>1)then
call wout(' d^2 Phi/dy dy                      : '//trim(i2s(nele(xyzz&
&yaaad2)*ndet-xyzzyaaam2(2))))
if(dimensionality==3)then
call wout(' d^2 Phi/dz dz                      : '//trim(i2s(nele(xyzz&
&yaaad2)*ndet-xyzzyaaam2(3))))
endif
call wout(' d^2 Phi/dx dy                      : '//trim(i2s(nele(xyzz&
&yaaad2)*ndet-xyzzyaaam2(4))))
if(dimensionality==3)then
call wout(' d^2 Phi/dx dz                      : '//trim(i2s(nele(xyzz&
&yaaad2)*ndet-xyzzyaaam2(5))))
call wout(' d^2 Phi/dy dz                      : '//trim(i2s(nele(xyzz&
&yaaad2)*ndet-xyzzyaaam2(6))))
endif
endif
call wout()
endif
endif
deallocate(xyzzyaabh2,xyzzyaabi2,xyzzyaabj2)
if(present(sderivs))deallocate(xyzzyaabk2)
end subroutine orbital_check
subroutine check_kinetic(is,js)
use slaarnach
use slaarnacs
use slaarnaag, only : czero
use slaarnabg,  only : dimensionality
use slaarnabt, only : ddot
implicit none
integer,intent(in) :: is,js
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3(nsampling_levels),&
&xyzzyaaae3(nsampling_levels),xyzzyaaaf3(nsampling_levels),xyzzyaaag3(&
&nsampling_levels),xyzzyaaah3,xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3
real(dp) xyzzyaaal3(nsampling_levels),xyzzyaaam3(nsampling_levels),xyz&
&zyaaan3(nsampling_levels),xyzzyaaao3(nsampling_levels),xyzzyaaap3(nsa&
&mpling_levels),xyzzyaaaq3(nsampling_levels),xyzzyaaar3(nsampling_leve&
&ls),xyzzyaaas3(nsampling_levels),xyzzyaaat3(3),xyzzyaaau3,xyzzyaaav3,&
&xyzzyaaaw3(3),xyzzyaaax3(3),xyzzyaaay3,xyzzyaaaz3,xyzzyaaba3
real(dp),parameter :: xyzzyaabb3=epsilon(1.d0)
real(dp),parameter :: xyzzyaabc3=1.d150
complex(dp) xyzzyaabd3(3,2,nsampling_levels),xyzzyaabe3(3),xyzzyaabf3
logical xyzzyaabg3(nsampling_levels),xyzzyaabh3(nsampling_levels),xyzz&
&yaabi3(nsampling_levels),xyzzyaabj3(nsampling_levels),xyzzyaabk3(nsam&
&pling_levels),xyzzyaabl3(nsampling_levels),isnan,isinf,nan,inf
character(80) tmpr
if(ke_verbose)then
call wout('Kinetic energy check')
call wout('====================')
endif
call get_rsele(is)
call wfn_assess_check_kinetic(is,xyzzyaaaa3,ke_verbose)
if(xyzzyaaaa3==0)xyzzyaaaa3=1
if(xyzzyaaaa3==on_top_ii.or.xyzzyaaaa3==on_top_jj)then
do xyzzyaaaa3=1,netot
if(xyzzyaaaa3==on_top_ii.or.xyzzyaaaa3==on_top_jj)cycle
exit
enddo
endif
xyzzyaaab3=which_ie(xyzzyaaaa3)
xyzzyaaac3=which_spin(xyzzyaaaa3)
if(ke_verbose)call wout('Checking KE for particle '//trim(i2s(xyzzyaaa&
&b3))//' of spin '//trim(i2s(xyzzyaaac3)))
isnan=.false.
isinf=.false.
do xyzzyaaah3=1,nsampling_levels
call wfn_loggrad(xyzzyaaaa3,is,xyzzyaaah3,xyzzyaabe3,prefetch_sd=.true&
&.,isnan=nan,isinf=inf)
xyzzyaaaw3=dble(xyzzyaabe3)
xyzzyaaal3(xyzzyaaah3)=ddot(3,xyzzyaaaw3(1),1,xyzzyaaaw3(1),1)
if(complex_wf)then
xyzzyaaax3=aimag(xyzzyaabe3)
xyzzyaaal3(xyzzyaaah3)=xyzzyaaal3(xyzzyaaah3)+ddot(3,xyzzyaaax3(1),1,x&
&yzzyaaax3(1),1)
endif
xyzzyaaal3(xyzzyaaah3)=0.5d0*xyzzyaaal3(xyzzyaaah3)*inv_pmass(xyzzyaaa&
&c3)
xyzzyaabi3(xyzzyaaah3)=abs(xyzzyaaal3(xyzzyaaah3))<xyzzyaabb3
call wfn_loglap(xyzzyaaaa3,is,xyzzyaaah3,xyzzyaabf3,isnan=nan,isinf=in&
&f)
xyzzyaaam3(xyzzyaaah3)=-0.25d0*dble(xyzzyaabf3)*inv_pmass(xyzzyaaac3)
xyzzyaabh3(xyzzyaaah3)=abs(xyzzyaaam3(xyzzyaaah3))<xyzzyaabb3
isnan=isnan.or.nan
isinf=isinf.or.inf
enddo
if(isnan.or.isinf)then
if(ke_verbose)then
call wout('Analytical derivatives misbehave. Test failed.')
call wout()
else
call wout('Kinetic energy test failed: analytical derivatives misbehav&
&e.')
call wout()
endif
return
endif
if(ke_verbose)then
call wout('------------------------------------')
do xyzzyaaah3=1,nsampling_levels
call wout(trim(sampling_level_name(xyzzyaaah3))//' :')
tmpr=r2s(xyzzyaaal3(xyzzyaaah3),'(e16.8)')
if(xyzzyaabi3(xyzzyaaah3))tmpr=trim(tmpr)//' [assumed zero]'
call wout(' Analytical Fisq    : '//trim(tmpr))
tmpr=r2s(xyzzyaaam3(xyzzyaaah3),'(e16.8)')
if(xyzzyaabh3(xyzzyaaah3))tmpr=trim(tmpr)//' [assumed zero]'
call wout(' Analytical Ti      : '//trim(tmpr))
enddo
endif
xyzzyaaad3=0
xyzzyaaae3=0
xyzzyaaap3=0.d0
xyzzyaaaq3=0.d0
xyzzyaabg3=.false.
xyzzyaaai3=0
xyzzyaaau3=xyzzyaaae1
do while(.not.all(xyzzyaabg3).and.xyzzyaaai3<xyzzyaaag1)
xyzzyaaai3=xyzzyaaai3+1
xyzzyaaav3=1.d0/xyzzyaaau3
isnan=.false.
isinf=.false.
do xyzzyaaaj3=1,dimensionality
do xyzzyaaak3=1,2
xyzzyaaat3=rele_scr(:,xyzzyaaaa3,is)
if(xyzzyaaak3==1)then
xyzzyaaat3(xyzzyaaaj3)=xyzzyaaat3(xyzzyaaaj3)-xyzzyaaau3
else
xyzzyaaat3(xyzzyaaaj3)=xyzzyaaat3(xyzzyaaaj3)+xyzzyaaau3
endif
call define_config_oneelec(xyzzyaaaa3,is,js,xyzzyaaat3,sele_scr(xyzzya&
&aaa3,is))
do xyzzyaaah3=1,nsampling_levels
if(xyzzyaabg3(xyzzyaaah3))cycle
call wfn_ratio(is,js,xyzzyaaah3,ratio=xyzzyaabd3(xyzzyaaaj3,xyzzyaaak3&
&,xyzzyaaah3),isnan=nan,isinf=inf)
isnan=isnan.or.nan
isinf=isinf.or.inf
enddo
enddo
enddo
do xyzzyaaah3=1,nsampling_levels
xyzzyaabj3(xyzzyaaah3)=.false.
if(xyzzyaabg3(xyzzyaaah3))cycle
if(any(abs(dble(xyzzyaabd3(1:dimensionality,1:2,xyzzyaaah3)))>xyzzyaab&
&c3).or.isnan.or.isinf)then
xyzzyaaan3(xyzzyaaah3)=sqrt(huge(1.d0))
xyzzyaaao3(xyzzyaaah3)=sqrt(huge(1.d0))
xyzzyaabj3(xyzzyaaah3)=.true.
cycle
endif
xyzzyaabe3=czero
xyzzyaaay3=0.d0
do xyzzyaaaj3=1,dimensionality
xyzzyaabe3(xyzzyaaaj3)=(xyzzyaabd3(xyzzyaaaj3,2,xyzzyaaah3)-xyzzyaabd3&
&(xyzzyaaaj3,1,xyzzyaaah3))*0.5d0*xyzzyaaav3
xyzzyaaay3=xyzzyaaay3+(dble(xyzzyaabd3(xyzzyaaaj3,2,xyzzyaaah3))+dble(&
&xyzzyaabd3(xyzzyaaaj3,1,xyzzyaaah3))-2.d0)*xyzzyaaav3*xyzzyaaav3-dble&
&(xyzzyaabe3(xyzzyaaaj3))**2-aimag(xyzzyaabe3(xyzzyaaaj3))**2
enddo
xyzzyaaaw3=dble(xyzzyaabe3)
xyzzyaaan3(xyzzyaaah3)=ddot(3,xyzzyaaaw3(1),1,xyzzyaaaw3(1),1)
if(complex_wf)then
xyzzyaaax3=aimag(xyzzyaabe3)
xyzzyaaan3(xyzzyaaah3)=xyzzyaaan3(xyzzyaaah3)+ddot(3,xyzzyaaax3(1),1,x&
&yzzyaaax3(1),1)
endif
xyzzyaaan3(xyzzyaaah3)=0.5d0*xyzzyaaan3(xyzzyaaah3)*inv_pmass(xyzzyaaa&
&c3)
xyzzyaaao3(xyzzyaaah3)=-xyzzyaaay3*0.25d0*inv_pmass(xyzzyaaac3)
enddo
if(ke_verbose)then
call wout('------------------------------------')
call wout('Iteration           : '//trim(i2s(xyzzyaaai3)))
tmpr=r2s(xyzzyaaau3,'(e16.8)')
call wout('Stepsize            : '//trim(tmpr))
endif
do xyzzyaaah3=1,nsampling_levels
if(xyzzyaabg3(xyzzyaaah3))then
if(ke_verbose)call wout(trim(sampling_level_name(xyzzyaaah3))//' : don&
&e')
cycle
endif
if(xyzzyaabj3(xyzzyaaah3))then
xyzzyaaaz3=sqrt(huge(1.d0))
elseif(.not.xyzzyaabi3(xyzzyaaah3))then
xyzzyaaaz3=abs(1.d0-xyzzyaaan3(xyzzyaaah3)/xyzzyaaal3(xyzzyaaah3))
else
xyzzyaaaz3=abs(xyzzyaaan3(xyzzyaaah3))
endif
if(xyzzyaabj3(xyzzyaaah3))then
xyzzyaaba3=sqrt(huge(1.d0))
elseif(.not.xyzzyaabh3(xyzzyaaah3))then
xyzzyaaba3=abs(1.d0-xyzzyaaao3(xyzzyaaah3)/xyzzyaaam3(xyzzyaaah3))
else
xyzzyaaba3=abs(xyzzyaaao3(xyzzyaaah3))
endif
if(xyzzyaaaz3<xyzzyaaac1)then
xyzzyaaad3(xyzzyaaah3)=max(xyzzyaaad3(xyzzyaaah3),3)
elseif(xyzzyaaaz3<xyzzyaaab1)then
xyzzyaaad3(xyzzyaaah3)=max(xyzzyaaad3(xyzzyaaah3),2)
elseif(xyzzyaaaz3<xyzzyaaaa1)then
xyzzyaaad3(xyzzyaaah3)=max(xyzzyaaad3(xyzzyaaah3),1)
endif
if(xyzzyaaba3<xyzzyaaac1)then
xyzzyaaae3(xyzzyaaah3)=max(xyzzyaaae3(xyzzyaaah3),3)
elseif(xyzzyaaba3<xyzzyaaab1)then
xyzzyaaae3(xyzzyaaah3)=max(xyzzyaaae3(xyzzyaaah3),2)
elseif(xyzzyaaba3<xyzzyaaaa1)then
xyzzyaaae3(xyzzyaaah3)=max(xyzzyaaae3(xyzzyaaah3),1)
endif
if(xyzzyaaaz3<xyzzyaaap3(xyzzyaaah3).or.xyzzyaaai3==1)then
xyzzyaaaf3(xyzzyaaah3)=xyzzyaaai3
xyzzyaaap3(xyzzyaaah3)=xyzzyaaaz3
xyzzyaaar3(xyzzyaaah3)=xyzzyaaan3(xyzzyaaah3)
xyzzyaabk3(xyzzyaaah3)=xyzzyaabj3(xyzzyaaah3)
endif
if(xyzzyaaba3<xyzzyaaaq3(xyzzyaaah3).or.xyzzyaaai3==1)then
xyzzyaaag3(xyzzyaaah3)=xyzzyaaai3
xyzzyaaaq3(xyzzyaaah3)=xyzzyaaba3
xyzzyaaas3(xyzzyaaah3)=xyzzyaaao3(xyzzyaaah3)
xyzzyaabl3(xyzzyaaah3)=xyzzyaabj3(xyzzyaaah3)
endif
xyzzyaabg3(xyzzyaaah3)=xyzzyaaad3(xyzzyaaah3)==3.and.xyzzyaaae3(xyzzya&
&aah3)==3
if(ke_verbose)then
call wout(trim(sampling_level_name(xyzzyaaah3))//' :')
if(xyzzyaabj3(xyzzyaaah3))then
call wout(' Numerical Fisq     : diverges')
call wout(' Rel. diff. in Fisq : undefined')
else
tmpr=r2s(xyzzyaaan3(xyzzyaaah3),'(e16.8)')
call wout(' Numerical Fisq     : '//trim(tmpr))
tmpr=r2s(xyzzyaaaz3,'(e16.8)')
call wout(' Rel. diff. in Fisq : '//trim(tmpr))
endif
call wout(' Fisq quality       : '//trim(xyzzyaaad1(xyzzyaaad3(xyzzyaa&
&ah3))))
if(xyzzyaabj3(xyzzyaaah3))then
call wout(' Numerical Ti       : diverges')
call wout(' Rel. diff. in Ti   : undefined')
else
tmpr=r2s(xyzzyaaao3(xyzzyaaah3),'(e16.8)')
call wout(' Numerical Ti       : '//trim(tmpr))
tmpr=r2s(xyzzyaaba3,'(e16.8)')
call wout(' Rel. diff. in Ti   : '//trim(tmpr))
endif
call wout(' Ti quality         : '//trim(xyzzyaaad1(xyzzyaaae3(xyzzyaa&
&ah3))))
endif
enddo
xyzzyaaau3=xyzzyaaau3*xyzzyaaaf1
enddo
if(ke_verbose)then
call wout('====================================')
call wout()
endif
if(any(xyzzyaaad3==0).or.any(xyzzyaaae3==0))then
call wout('Wave function failed kinetic energy check')
call wout('-----------------------------------------')
do xyzzyaaah3=1,nsampling_levels
if(xyzzyaaad3(xyzzyaaah3)==0.or.xyzzyaaae3(xyzzyaaah3)==0)then
call wout('Test failed for '//trim(sampling_level_name(xyzzyaaah3))//'&
&:')
tmpr=r2s(xyzzyaaal3(xyzzyaaah3),'(e16.8)')
if(xyzzyaabi3(xyzzyaaah3))tmpr=trim(tmpr)//' [assumed zero]'
call wout(' Analytical Fisq     : '//trim(tmpr))
if(xyzzyaabk3(xyzzyaaah3))then
tmpr='diverges'
else
tmpr=r2s(xyzzyaaar3(xyzzyaaah3),'(e16.8)')
endif
call wout(' Numerical Fisq      : '//trim(tmpr)//' (for stepsize #'//t&
&rim(i2s(xyzzyaaaf3(xyzzyaaah3)))//')')
tmpr=r2s(xyzzyaaam3(xyzzyaaah3),'(e16.8)')
if(xyzzyaabh3(xyzzyaaah3))tmpr=trim(tmpr)//' [assumed zero]'
call wout(' Analytical Ti       : '//trim(tmpr))
if(xyzzyaabl3(xyzzyaaah3))then
tmpr='diverges'
else
tmpr=r2s(xyzzyaaas3(xyzzyaaah3),'(e16.8)')
endif
call wout(' Numerical Ti        : '//trim(tmpr)//' (for stepsize #'//t&
&rim(i2s(xyzzyaaag3(xyzzyaaah3)))//')')
endif
enddo
if(ke_forgive)then
call wout()
call errwarn('CHECK_KE','This may indicate a problem.  Continuing run,&
& however, because the KE_FORGIVE flag is set to T.')
call wout()
else
call errstop('CHECK_KE','Stopping...')
endif
endif
call wout('Kinetic energy check performed.')
do xyzzyaaah3=1,nsampling_levels
call wout(' '//trim(sampling_level_name(xyzzyaaah3))//' - '//'gradient&
&: '//trim(xyzzyaaad1(xyzzyaaad3(xyzzyaaah3)))//', Laplacian: '//trim(&
&xyzzyaaad1(xyzzyaaae3(xyzzyaaah3)))//'.')
enddo
call wout('End of report.')
call wout()
end subroutine check_kinetic
end module slaarnacr

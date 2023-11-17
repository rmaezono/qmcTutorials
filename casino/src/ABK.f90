module slaarnabk
use dsp
use parallel
use slaarnaag,   only : two_over_root_pi,one_over_fourpisquared,root_p&
&i,c_one,pi,pi_over_four,third,fourpi,czero,four_over_root_pi,pi_over_&
&two,sixth,one_over_pi,one_over_twopi,twopi,one_over_root_pi,euler
use slaarnaaq,      only : qmc_density_mpc,permit_den_symm,den_nsets,d&
&en_sc,den_gvec,expval_den_mpc,expval_ngvec_truncated,expval_gvec,expv&
&al_complex_den
use file_utils,  only : skip,open_units
use format_utils,only : wout,r2s,i2s,wordwrap
use slaarnabg,    only : periodicity,a1,a2,a3,area,first_sr_in_star,nu&
&m_g,s_modsq,sr_modsq,first_s_in_star,sr_lattice,sr_lmng,volume,b1,b2,&
&b3,s_lattice,pa1,pa2,pa3,basis,atno,nbasis,ngpri,pvolume,pb1,pb2,pb3,&
&npcells,atom_basis_type,wigner_seitz_radius,painv,lsize,rionion,hard_&
&sphere,hard_diam,hard_op_spins,hard_diam_sq,dimensionality
use run_control, only : errstop,errstop_master,timer,check_alloc
use store,       only : nspin,nele,ewald_control,self_term,netot,pchar&
&ge,no_families,calc_field,fam_charge,interaction,heg_orbtype,which_fa&
&m,electron_system,which_spin,inv_pmass,heg_zlayer,heg_ylayer,heg_laye&
&r,heg_nlayers,fix_holes,isgen_mpc,hartree_xc,int_name,man_int_params,&
&man_int_op_spins,rstar
implicit none
private
public init_interactions,ee_potential,en_potential,nuclear_repulsion_e&
&nergy,ewald_2d,ewald_3d,mpc_correction,interaction_mpc_present,biex3p&
&ot,compute_fourier_basis_mpc,interaction_mpc_use,interaction1_present&
&,homogeneous_mpc,mpc_hartree,man_int_params,man_int_op_spins,blip_mpc&
&,man_int_poly_order,man_int_poly_coeffs,man_int_poly_cutoff,eval_fsco&
&rr_consts,man_int_tilt_poly_order,man_int_tilt_poly_coeffs,iinteracti&
&on,icoulomb,inone,ilogarithmic,pot_2d_int,i2d_int,pseudo_ewald
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae1,xyzzyaa&
&af1,xyzzyaaag1,xyzzyaaah1,xyzzyaaai1,xyzzyaaaj1,xyzzyaaak1,xyzzyaaal1
integer,allocatable :: xyzzyaaam1(:),xyzzyaaan1(:),xyzzyaaao1(:,:)
real(dp) xyzzyaaap1,xyzzyaaaq1,xyzzyaaar1,xyzzyaaas1,xyzzyaaat1,xyzzya&
&aau1,xyzzyaaav1,xyzzyaaaw1,xyzzyaaax1,xyzzyaaay1,xyzzyaaaz1,xyzzyaaba&
&1,xyzzyaabb1,xyzzyaabc1,xyzzyaabd1,xyzzyaabe1,xyzzyaabf1,xyzzyaabg1,x&
&yzzyaabh1,xyzzyaabi1,xyzzyaabj1,xyzzyaabk1,xyzzyaabl1,xyzzyaabm1(0:40&
&03),mpc_correction,mpc_hartree,xyzzyaabn1,xyzzyaabo1,xyzzyaabp1,xyzzy&
&aabq1,xyzzyaabr1
real(dp),allocatable :: xyzzyaabs1(:),xyzzyaabt1(:),xyzzyaabu1(:),xyzz&
&yaabv1(:),xyzzyaabw1(:),xyzzyaabx1(:),xyzzyaaby1(:),xyzzyaabz1(:),xyz&
&zyaaca1(:),xyzzyaacb1(:,:),xyzzyaacc1(:,:),xyzzyaacd1(:,:)
complex(dp),allocatable :: xyzzyaace1(:),xyzzyaacf1(:),xyzzyaacg1(:)
!$omp threadprivate(xyzzyaace1,xyzzyaacf1,xyzzyaacg1)
logical interaction1_present,homogeneous_mpc,interaction_mpc_present,i&
&nteraction_mpc_use,pseudo_ewald
integer xyzzyaach1,xyzzyaaci1,xyzzyaacj1,xyzzyaack1
integer,allocatable :: xyzzyaacl1(:,:)
real(dp),allocatable :: xyzzyaacm1(:,:),xyzzyaacn1(:,:)
complex(dp),allocatable :: xyzzyaaco1(:,:)
integer,parameter :: inone=0,icoulomb=1,xyzzyaacp1=2,xyzzyaacq1=3,xyzz&
&yaacr1=4,xyzzyaacs1=5,xyzzyaact1=6,xyzzyaacu1=7,xyzzyaacv1=8,ilogarit&
&hmic=9,i2d_int=10,xyzzyaacw1=11,xyzzyaacx1=12,xyzzyaacy1=13,xyzzyaacz&
&1=14
integer iinteraction
integer xyzzyaada1(3)
real(dp) xyzzyaadb1(3),xyzzyaadc1(3)
real(dp) xyzzyaadd1
real(dp),allocatable :: xyzzyaade1(:,:,:)
logical blip_mpc
integer man_int_poly_order
integer man_int_tilt_poly_order
real(dp),allocatable :: man_int_poly_coeffs(:)
real(dp),allocatable :: man_int_tilt_poly_coeffs(:)
real(dp) man_int_poly_cutoff
logical,parameter :: xyzzyaadf1=.false.
contains
subroutine init_interactions
use slaarnaas, only : heg_mpc_long,half_mpc_constant,excite_heg,harmwi&
&re_b
implicit none
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaal2&
&,xyzzyaaam2,xyzzyaaan2
real(dp) xyzzyaaao2,xyzzyaaap2,xyzzyaaaq2,xyzzyaaar2,xyzzyaaas2,xyzzya&
&aat2,xyzzyaaau2,xyzzyaaav2,xyzzyaaaw2,xyzzyaaax2,xyzzyaaay2,xyzzyaaaz&
&2,xyzzyaaba2,xyzzyaabb2,xyzzyaabc2,xyzzyaabd2,xyzzyaabe2,xyzzyaabf2,x&
&yzzyaabg2
logical xyzzyaabh2
character(80) tmpr
call open_units(xyzzyaaai1,xyzzyaaal2)
if(xyzzyaaal2/=0)call errstop('INIT_INTERACTIONS','Unable to find free&
& i/o unit.')
if(am_master)then
call wout('Interactions')
call wout('============')
select case(trim(interaction))
case('none')
call wout('Non-interacting system.')
call wout()
case('coulomb','ewald','mpc','ewald_mpc','mpc_ewald')
call wout('Interaction type: Coulomb')
case('manual')
continue
case default
call errstop('INIT_INTERACTIONS','Unknown interaction type '//trim(int&
&eraction))
end select
endif
homogeneous_mpc=.false.
self_term=0.d0
select case(trim(interaction))
case('coulomb')
if(am_master)then
call wout('System is aperiodic. Using 1/r only.')
call wout()
endif
iinteraction=icoulomb
case('manual')
call xyzzyaadg1(iinteraction)
case default
select case(trim(interaction))
case('ewald')
iinteraction=xyzzyaacp1
case('mpc')
iinteraction=xyzzyaacq1
case('ewald_mpc','mpc_ewald')
iinteraction=xyzzyaacr1
end select
if(pseudo_ewald)call xyzzyaadg1(xyzzyaaan2)
if(interaction_mpc_present)then
if(am_master)then
if(periodicity==2)then
if(qmc_density_mpc)then
call wout('Setup 2D MPC interaction (using QMC density).')
else
call wout('Setup 2D MPC interaction.')
endif
else
if(qmc_density_mpc)then
call wout('Setup 3D MPC interaction (using QMC density).')
else
call wout('Setup 3D MPC interaction.')
endif
endif
endif
if(trim(atom_basis_type)=='none')then
if(.not.any(heg_orbtype<0).and..not.any(heg_orbtype==2).and.electron_s&
&ystem)then
inquire(file='mpc.data',exist=xyzzyaabh2)
if(.not.xyzzyaabh2)then
if(excite_heg)call errstop('INIT_INTERACTIONS','Cannot use homogeneous&
& MPC for HEG excited states. Please supply mpc.data file.')
if(am_master)call wout('Calculating long-range MPC interaction.')
homogeneous_mpc=.true.
if(blip_mpc)call errstop('INIT_INTERACTIONS','You should set blip_mpc=&
&F for homogeneous systems.')
call heg_mpc_long
else
homogeneous_mpc=.false.
if(am_master)then
call wout('Reading in FT of f(r) and density from mpc.data.')
call xyzzyaaed1
endif
if(nnodes>1)call xyzzyaaee1
endif
else
call errstop('INIT_INTERACTIONS','Need to check whether the MPC intera&
&ction works for the system requested in the FREE_PARTICLES block in i&
&nput.')
endif
else
if(periodicity==1)call errstop('INIT COULOMB','MPC interaction not imp&
&lemented for one-dimensional periodic systems.')
if(am_master.and.(.not.isgen_mpc))call xyzzyaaed1
if(nnodes>1)call xyzzyaaee1
endif
if(am_master.and.homogeneous_mpc)call wout('Long-range MPC constant (a&
&u): ',2.d0*half_mpc_constant)
if(periodicity==3.and.nspin<3.and.blip_mpc)call xyzzyaaef1
endif
select case(periodicity)
case(1)
xyzzyaabg1=a1(1)
xyzzyaabj1=1.d0/xyzzyaabg1
xyzzyaabf1=-xyzzyaabj1
xyzzyaabh1=xyzzyaabg1/24.d0
xyzzyaaap1=-xyzzyaabg1*xyzzyaabg1*xyzzyaabh1*7.d0/240.d0
xyzzyaaaa2=int(ewald_control*0.036d0+3.6d0)
xyzzyaaaa1=max(xyzzyaaaa2+xyzzyaaaa2+1,3)
xyzzyaaaq1=xyzzyaaaa1*xyzzyaabg1*0.5d0
xyzzyaabe2=0.d0
do xyzzyaaac2=2,xyzzyaaaa1
xyzzyaabe2=xyzzyaabe2+1.d0/abs(s_lattice(1,xyzzyaaac2))
enddo
xyzzyaabf2=1.d0/(xyzzyaaaq1*xyzzyaaaq1)
xyzzyaabg2=xyzzyaabf2*xyzzyaabf2
if(harmwire_b<0.d0)then
self_term=(2*xyzzyaabf1*log(2*xyzzyaaaq1)-2*xyzzyaabh1*xyzzyaabf2-12*x&
&yzzyaaap1*xyzzyaabg2+xyzzyaabe2)*0.5d0
else
self_term=0.5d0*xyzzyaadj1(harmwire_b)
endif
if(am_master)then
if(ngpri<xyzzyaaaa1)then
call wout('Too few direct lattice vectors for 1D Madelung sums')
call wout('ngpri: '//trim(i2s(ngpri)))
call wout('nreal:'//trim(i2s(xyzzyaaaa1)))
call wout()
call errstop('INIT_INTERACTIONS','Quitting.')
endif
call wout()
call wout('Setup 1D periodic Coulomb interaction.')
call wout()
tmpr=r2s(xyzzyaabg1,'(es20.8)')
call wout('Length of unit cell          (au) :  '//trim(tmpr))
tmpr=r2s(xyzzyaaaq1,'(es20.8)')
call wout('Length of quantum zone       (au) :  '//trim(tmpr))
call wout('No. real space lattice vectors    :  '//trim(i2s(xyzzyaaaa1&
&)))
call wout('Electron self-image term v_M (au) :  ',2.d0*self_term)
call wout()
endif
case(2)
call xyzzyaady1
xyzzyaaar2=1.d0/sqrt(area)
xyzzyaaao2=(100.d0+ewald_control)*xyzzyaaar2*0.01d0
xyzzyaabe1=10.d0*xyzzyaaao2*xyzzyaaao2
xyzzyaabd1=2.4d0*xyzzyaaar2
xyzzyaaax1=xyzzyaabd1*xyzzyaabd1
xyzzyaaar1=1.d0/xyzzyaabd1
xyzzyaaat2=0.25d0/xyzzyaaax1
xyzzyaabi1=xyzzyaabd1*two_over_root_pi
xyzzyaaav1=24.d0/xyzzyaaax1
xyzzyaaaw1=4.d0*xyzzyaaav1
xyzzyaabk1=-(2.d0*root_pi*xyzzyaaar1)/area
xyzzyaabl1=2.d0*xyzzyaabk1
xyzzyaaad2=2
xyzzyaaab1=1
xyzzyaaas2=sr_modsq(xyzzyaaad2)
do
do xyzzyaaac2=first_sr_in_star(xyzzyaaad2),first_sr_in_star(xyzzyaaad2&
&+1)-1,2
xyzzyaaab1=xyzzyaaab1+1
if(xyzzyaaab1>=(num_g+1)/2)call errstop('INIT_INTERACTIONS:','Ran out &
&of predefined reciprocal lattice vectors.')
enddo
xyzzyaaad2=xyzzyaaad2+1
xyzzyaaas2=sr_modsq(xyzzyaaad2)
if(xyzzyaaas2*one_over_fourpisquared>=xyzzyaabe1)exit
enddo
xyzzyaaae2=xyzzyaaad2-1
allocate(xyzzyaaao1(2,xyzzyaaab1),xyzzyaaam1(2*xyzzyaaab1),xyzzyaabs1(&
&xyzzyaaab1),xyzzyaabt1(xyzzyaaab1),xyzzyaaca1(xyzzyaaab1),xyzzyaabw1(&
&xyzzyaaab1),xyzzyaabx1(xyzzyaaab1),xyzzyaabz1(xyzzyaaab1),stat=xyzzya&
&aam2)
call check_alloc(xyzzyaaam2,'INIT_INTERACTIONS','1')
xyzzyaaab1=1
xyzzyaaax2=0.d0
xyzzyaaac1=0
xyzzyaaad1=0
xyzzyaaaf1=0
xyzzyaaag1=0
do xyzzyaaad2=2,xyzzyaaae2
xyzzyaaas2=sr_modsq(xyzzyaaad2)
do xyzzyaaac2=first_sr_in_star(xyzzyaaad2),first_sr_in_star(xyzzyaaad2&
&+1)-1,2
xyzzyaaab1=xyzzyaaab1+1
xyzzyaaao1(1,xyzzyaaab1)=sr_lmng(1,xyzzyaaac2)
xyzzyaaao1(2,xyzzyaaab1)=sr_lmng(2,xyzzyaaac2)
xyzzyaaac1=min(xyzzyaaac1,xyzzyaaao1(1,xyzzyaaab1))
xyzzyaaad1=min(xyzzyaaad1,xyzzyaaao1(2,xyzzyaaab1))
xyzzyaaaf1=max(xyzzyaaaf1,xyzzyaaao1(1,xyzzyaaab1))
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaao1(2,xyzzyaaab1))
xyzzyaabc2=sqrt(xyzzyaaas2*one_over_fourpisquared)
xyzzyaabs1(xyzzyaaab1)=sr_lattice(1,xyzzyaaac2)
xyzzyaabt1(xyzzyaaab1)=sr_lattice(2,xyzzyaaac2)
xyzzyaabw1(xyzzyaaab1)=1.d0/(area*xyzzyaabc2)
xyzzyaabx1(xyzzyaaab1)=pi*xyzzyaaar1*xyzzyaabc2
xyzzyaabz1(xyzzyaaab1)=xyzzyaaas2*xyzzyaaat2
xyzzyaaca1(xyzzyaaab1)=sqrt(xyzzyaabs1(xyzzyaaab1)**2+xyzzyaabt1(xyzzy&
&aaab1)**2)
enddo
enddo
xyzzyaaai2=max(xyzzyaaaf1,1)
xyzzyaaaj2=max(xyzzyaaag1,1)
xyzzyaaaf2=min(xyzzyaaac1,-1)
xyzzyaaag2=min(xyzzyaaad1,-1)
allocate(xyzzyaace1(xyzzyaaaf2:xyzzyaaai2),xyzzyaacf1(xyzzyaaag2:xyzzy&
&aaaj2),stat=xyzzyaaam2)
call check_alloc(xyzzyaaam2,'INIT_INTERACTIONS','2')
xyzzyaace1(0)=c_one
xyzzyaacf1(0)=c_one
xyzzyaaad2=1
xyzzyaaaa1=0
xyzzyaaaw2=0.d0
do
xyzzyaaav2=s_modsq(xyzzyaaad2)*xyzzyaaax1
if(xyzzyaaav2>=16.d0)then
do xyzzyaaac2=first_s_in_star(xyzzyaaad2),first_s_in_star(xyzzyaaad2+1&
&)-1
xyzzyaaaa1=xyzzyaaaa1+1
if(xyzzyaaaa1>=num_g)call errstop('INIT_INTERACTIONS:','Ran out of pre&
&defined real space lattice vectors.')
xyzzyaaau2=((xyzzyaaav2+xyzzyaabb1)*exp(-xyzzyaaav2)/((xyzzyaaav2+xyzz&
&yaabc1)*(xyzzyaaav2+xyzzyaaav2)+xyzzyaaba1))
xyzzyaaaw2=xyzzyaaaw2-xyzzyaaau2
enddo
else
do xyzzyaaac2=first_s_in_star(xyzzyaaad2),first_s_in_star(xyzzyaaad2+1&
&)-1
xyzzyaaaa1=xyzzyaaaa1+1
if(xyzzyaaaa1>=num_g)call errstop('INIT_INTERACTIONS:','Ran out of pre&
&defined real space lattice vectors.')
xyzzyaaaa2=int(xyzzyaaaz1*xyzzyaaav2+0.5d0)
xyzzyaaap2=xyzzyaaay1*xyzzyaaaa2-xyzzyaaav2
xyzzyaaau2=((xyzzyaabm1(4*xyzzyaaaa2)*xyzzyaaap2+xyzzyaabm1(4*xyzzyaaa&
&a2+1))*xyzzyaaap2+xyzzyaabm1(4*xyzzyaaaa2+2))*xyzzyaaap2+xyzzyaabm1(4&
&*xyzzyaaaa2+3)
xyzzyaaaw2=xyzzyaaaw2+xyzzyaaau2
if(xyzzyaaav2>=1.d-16)xyzzyaaaw2=xyzzyaaaw2-sqrt(pi_over_four/xyzzyaaa&
&v2)
enddo
endif
xyzzyaaad2=xyzzyaaad2+1
if(s_modsq(xyzzyaaad2)>xyzzyaaav1)exit
enddo
xyzzyaaat1=0.d0
xyzzyaaas1=0.d0
xyzzyaaax2=0.d0
do xyzzyaaac2=2,xyzzyaaab1
call xyzzyaaea1(xyzzyaabx1(xyzzyaaac2),xyzzyaaaq2)
xyzzyaaax2=xyzzyaaaq2*xyzzyaabw1(xyzzyaaac2)+xyzzyaaax2
enddo
call xyzzyaadz1(xyzzyaaaq2,xyzzyaabd2)
xyzzyaaaw2=-xyzzyaaaw2*xyzzyaabi1
self_term=0.5d0*(xyzzyaaax2+xyzzyaaaw2+xyzzyaaaq2)
if(am_master)then
call wout()
call wout('Setup 2D Ewald interaction.')
call wout()
call wout('No. of G in reciprocal space sum  :  '//trim(i2s(xyzzyaaab1&
&-1)))
call wout('No. of vectors in real space sum  :  '//trim(i2s(xyzzyaaaa1&
&)))
tmpr=r2s(xyzzyaabe1,'(1e20.8)')
call wout('Maximum |G|^2 (au)                :  '//trim(tmpr))
tmpr=r2s(xyzzyaabd1,'(es20.8)')
call wout('Gamma^(1/2) = Gaussian half-width :  '//trim(tmpr))
call wout('Electron self-image term v_M (au) :  ',2.d0*self_term)
call wout()
endif
case(3)
call xyzzyaady1
xyzzyaabb2=fourpi/volume
xyzzyaaar2=(1.d0/volume)**third
xyzzyaaao2=(100.d0+ewald_control)*xyzzyaaar2*0.01d0
xyzzyaabe1=12.8d0*xyzzyaaao2*xyzzyaaao2
xyzzyaabd1=2.8d0*xyzzyaaar2
xyzzyaaax1=xyzzyaabd1*xyzzyaabd1
xyzzyaabi1=xyzzyaabd1*two_over_root_pi
xyzzyaaat2=0.25d0/xyzzyaaax1
xyzzyaaav1=24.d0/xyzzyaaax1
xyzzyaaaw1=4.d0*xyzzyaaav1
xyzzyaaad2=2
xyzzyaaab1=1
xyzzyaaas2=sr_modsq(xyzzyaaad2)
do
do xyzzyaaac2=first_sr_in_star(xyzzyaaad2),first_sr_in_star(xyzzyaaad2&
&+1)-1,2
xyzzyaaab1=xyzzyaaab1+1
if(xyzzyaaab1>=(num_g+1)/2)call errstop('INIT_INTERACTIONS:','Ran out &
&of predefined reciprocal lattice vectors.')
enddo
xyzzyaaad2=xyzzyaaad2+1
xyzzyaaas2=sr_modsq(xyzzyaaad2)
if(xyzzyaaas2*one_over_fourpisquared>=xyzzyaabe1)exit
enddo
allocate(xyzzyaaao1(3,xyzzyaaab1),xyzzyaabs1(xyzzyaaab1),xyzzyaabt1(xy&
&zzyaaab1),xyzzyaabu1(xyzzyaaab1),xyzzyaabv1(xyzzyaaab1),xyzzyaaam1(2*&
&xyzzyaaab1),xyzzyaaby1(xyzzyaaab1),stat=xyzzyaaam2)
call check_alloc(xyzzyaaam2,'INIT_INTERACTIONS','3')
xyzzyaaae2=xyzzyaaad2-1
xyzzyaaab1=1
xyzzyaaax2=0.d0
xyzzyaaac1=0
xyzzyaaad1=0
xyzzyaaae1=0
xyzzyaaaf1=0
xyzzyaaag1=0
xyzzyaaah1=0
do xyzzyaaad2=2,xyzzyaaae2
xyzzyaaas2=sr_modsq(xyzzyaaad2)
xyzzyaaba2=xyzzyaaas2*xyzzyaaat2
xyzzyaaaz2=xyzzyaabb2/xyzzyaaas2
do xyzzyaaac2=first_sr_in_star(xyzzyaaad2),first_sr_in_star(xyzzyaaad2&
&+1)-1,2
xyzzyaaab1=xyzzyaaab1+1
xyzzyaaao1(1,xyzzyaaab1)=sr_lmng(1,xyzzyaaac2)
xyzzyaaao1(2,xyzzyaaab1)=sr_lmng(2,xyzzyaaac2)
xyzzyaaao1(3,xyzzyaaab1)=sr_lmng(3,xyzzyaaac2)
xyzzyaaac1=min(xyzzyaaac1,xyzzyaaao1(1,xyzzyaaab1))
xyzzyaaad1=min(xyzzyaaad1,xyzzyaaao1(2,xyzzyaaab1))
xyzzyaaae1=min(xyzzyaaae1,xyzzyaaao1(3,xyzzyaaab1))
xyzzyaaaf1=max(xyzzyaaaf1,xyzzyaaao1(1,xyzzyaaab1))
xyzzyaaag1=max(xyzzyaaag1,xyzzyaaao1(2,xyzzyaaab1))
xyzzyaaah1=max(xyzzyaaah1,xyzzyaaao1(3,xyzzyaaab1))
xyzzyaabs1(xyzzyaaab1)=sr_lattice(1,xyzzyaaac2)
xyzzyaabt1(xyzzyaaab1)=sr_lattice(2,xyzzyaaac2)
xyzzyaabu1(xyzzyaaab1)=sr_lattice(3,xyzzyaaac2)
xyzzyaaby1(xyzzyaaab1)=2.d0*exp(-xyzzyaaba2)*xyzzyaaaz2
xyzzyaaax2=xyzzyaaax2+xyzzyaaby1(xyzzyaaab1)
enddo
enddo
xyzzyaaai2=max(xyzzyaaaf1,1)
xyzzyaaaj2=max(xyzzyaaag1,1)
xyzzyaaak2=max(xyzzyaaah1,1)
xyzzyaaaf2=min(xyzzyaaac1,-1)
xyzzyaaag2=min(xyzzyaaad1,-1)
xyzzyaaah2=min(xyzzyaaae1,-1)
!$omp parallel default(none) shared(xyzzyaaaf2,xyzzyaaai2,xyzzyaaag2,x&
!$omp &yzzyaaaj2,xyzzyaaah2,xyzzyaaak2) private(xyzzyaaam2)
allocate(xyzzyaace1(xyzzyaaaf2:xyzzyaaai2),xyzzyaacf1(xyzzyaaag2:xyzzy&
&aaaj2),xyzzyaacg1(xyzzyaaah2:xyzzyaaak2),stat=xyzzyaaam2)
call check_alloc(xyzzyaaam2,'INIT_INTERACTIONS','4')
xyzzyaace1(0)=c_one
xyzzyaacf1(0)=c_one
xyzzyaacg1(0)=c_one
!$omp end parallel
if(xyzzyaaaf1/=0)then
call wordwrap('CASINO requires the list of supercell reciprocal lattic&
&e vectors to have the following property (for i even and > 0 ): G_i+1&
&(1:3) = -G_i(1:3) and G_i(1) < 0.  For some reason this is not the ca&
&se. This should not happen.')
call errstop('INIT_INTERACTIONS','maxng1 /= 0 in Ewald setup.')
endif
allocate(xyzzyaaan1(xyzzyaaac1:0),stat=xyzzyaaam2)
call check_alloc(xyzzyaaam2,'INIT_INTERACTIONS','5')
xyzzyaaab2=1
do xyzzyaaaa2=2,xyzzyaaab1
if(xyzzyaaao1(1,xyzzyaaaa2)==0)then
xyzzyaaam1(2*(xyzzyaaab2-1)+1)=xyzzyaaao1(2,xyzzyaaaa2)
xyzzyaaam1(2*(xyzzyaaab2-1)+2)=xyzzyaaao1(3,xyzzyaaaa2)
xyzzyaabv1(xyzzyaaab2)=xyzzyaaby1(xyzzyaaaa2)
xyzzyaaab2=xyzzyaaab2+1
endif
enddo
xyzzyaaan1(0)=xyzzyaaab2-1
do xyzzyaaac2=-1,xyzzyaaac1,-1
do xyzzyaaaa2=2,xyzzyaaab1
if(xyzzyaaao1(1,xyzzyaaaa2)==xyzzyaaac2)then
xyzzyaaam1(2*(xyzzyaaab2-1)+1)=xyzzyaaao1(2,xyzzyaaaa2)
xyzzyaaam1(2*(xyzzyaaab2-1)+2)=xyzzyaaao1(3,xyzzyaaaa2)
xyzzyaabv1(xyzzyaaab2)=xyzzyaaby1(xyzzyaaaa2)
xyzzyaaab2=xyzzyaaab2+1
endif
enddo
xyzzyaaan1(xyzzyaaac2)=xyzzyaaab2-1-sum(xyzzyaaan1(xyzzyaaac2+1:0))
enddo
if(xyzzyaaab2/=xyzzyaaab1)call errstop('INIT_INTERACTIONS','j /= nreci&
&p in Ewald setup.')
xyzzyaaay2=two_over_root_pi*sqrt(xyzzyaaax1)
xyzzyaaad2=2
xyzzyaaaa1=1
xyzzyaaaw2=0.d0
do
xyzzyaaav2=s_modsq(xyzzyaaad2)*xyzzyaaax1
if(xyzzyaaav2>=16.d0)then
do xyzzyaaac2=first_s_in_star(xyzzyaaad2),first_s_in_star(xyzzyaaad2+1&
&)-1
xyzzyaaaa1=xyzzyaaaa1+1
if(xyzzyaaaa1>=num_g)call errstop('INIT_INTERACTIONS:','Ran out of pre&
&defined real space lattice vectors.')
xyzzyaaau2=(xyzzyaaav2+xyzzyaabb1)*exp(-xyzzyaaav2)/((xyzzyaaav2+xyzzy&
&aabc1)*(xyzzyaaav2+xyzzyaaav2)+xyzzyaaba1)
xyzzyaaaw2=xyzzyaaaw2-xyzzyaaau2
enddo
else
do xyzzyaaac2=first_s_in_star(xyzzyaaad2),first_s_in_star(xyzzyaaad2+1&
&)-1
xyzzyaaaa1=xyzzyaaaa1+1
if(xyzzyaaaa1>=num_g)call errstop('INIT_INTERACTIONS:','Ran out of pre&
&defined real space lattice vectors.')
xyzzyaaaa2=int(xyzzyaaaz1*xyzzyaaav2+0.5d0)
xyzzyaaap2=xyzzyaaay1*xyzzyaaaa2-xyzzyaaav2
xyzzyaaau2=((xyzzyaabm1(4*xyzzyaaaa2)*xyzzyaaap2+xyzzyaabm1(4*xyzzyaaa&
&a2+1))*xyzzyaaap2+xyzzyaabm1(4*xyzzyaaaa2+2))*xyzzyaaap2+xyzzyaabm1(4&
&*xyzzyaaaa2+3)
if(xyzzyaaav2>=(1.d-16))xyzzyaaau2=xyzzyaaau2-sqrt(pi_over_four/xyzzya&
&aav2)
xyzzyaaaw2=xyzzyaaaw2+xyzzyaaau2
enddo
endif
xyzzyaaad2=xyzzyaaad2+1
if(s_modsq(xyzzyaaad2)>xyzzyaaav1)exit
enddo
xyzzyaaau1=pi/(xyzzyaaax1*volume)
xyzzyaaaw2=-xyzzyaaaw2*xyzzyaabi1
self_term=0.5d0*(xyzzyaaaw2+xyzzyaaax2-xyzzyaaay2-xyzzyaaau1)
if(am_master)then
call wout()
call wout('Setup 3D Ewald interaction.')
call wout()
call wout('No. of G in reciprocal space sum   :  '//trim(i2s(xyzzyaaab&
&1-1)))
call wout('No. of vectors in real space sum   :  '//trim(i2s(xyzzyaaaa&
&1)))
tmpr=r2s(xyzzyaabe1,'(1e20.8)')
call wout('Maximum |G|^2 (au)                 :  '//trim(tmpr))
tmpr=r2s(xyzzyaabd1,'(es20.8)')
call wout('Gamma^(1/2) = Gaussian half-width  :  '//trim(tmpr))
call wout('Electron self-image term v_M (au)  :  ',2.d0*self_term)
call wout()
endif
end select
end select
if(.not.electron_system)then
allocate(xyzzyaacd1(nspin,nspin),stat=xyzzyaaam2)
call check_alloc(xyzzyaaam2,'INIT_INTERACTIONS','pcharge_pcharge')
do xyzzyaaab2=1,nspin
do xyzzyaaaa2=1,nspin
xyzzyaacd1(xyzzyaaaa2,xyzzyaaab2)=pcharge(xyzzyaaaa2)*pcharge(xyzzyaaa&
&b2)
enddo
enddo
endif
end subroutine init_interactions
subroutine xyzzyaadg1(iinteraction)
implicit none
integer,intent(inout) :: iinteraction
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3
real(dp) xyzzyaaaf3,xyzzyaaag3,xyzzyaaah3,xyzzyaaai3
real(dp) xyzzyaaaj3
real(dp),parameter :: xyzzyaaak3=1.d-8
self_term=0.d0
select case(trim(int_name))
case('square_well')
if(am_master)then
call wout('Interaction type: square well potential (manual)')
call wout('Height          : ',man_int_params(2))
call wout('Width           : ',man_int_params(1))
endif
if(man_int_params(1)<0.d0)call errstop_master('INIT_MAN_INT','Well_wid&
&th must be positive.')
iinteraction=xyzzyaacs1
case('poschl_teller')
xyzzyaadd1=2.d0*man_int_params(2)*man_int_params(1)**2*inv_pmass(1)
if(am_master)then
call wout('Interaction type: Poschl-Teller potential (manual)')
call wout('v_0             : ',man_int_params(2))
call wout('mu              : ',man_int_params(1))
endif
if(man_int_params(1)<0.d0)call errstop_master('INIT_MAN_INT','Well_wid&
&th must be positive.')
if(.not.allocated(inv_pmass))call errstop_master('INIT_MAN_INT','Array&
& ''inv_pmass'' not allocated. Should''t happen.')
do xyzzyaaaa3=2,nspin
if(inv_pmass(xyzzyaaaa3)/=inv_pmass(1))call errstop_master('INIT_MAN_I&
&NT','Masses of all particles should be equal when using the Poschl-Te&
&ller potential.')
enddo
xyzzyaaaj3=xyzzyaadn1(wigner_seitz_radius)
if(am_master)call wout('Potential at WSR: ',xyzzyaaaj3)
if(abs(xyzzyaaaj3)>1.d-5)call errstop_master('INIT_MAN_INT','Poschl-Te&
&ller interaction > 1.d-5 at radius of sphere inscribed in Wigner-Seit&
&z cell.')
iinteraction=xyzzyaact1
case('hard_sphere')
hard_sphere=.true.
hard_op_spins=man_int_op_spins
hard_diam=2.d0*man_int_params(1)
hard_diam_sq=hard_diam**2
xyzzyaabn1=man_int_params(2)
if(xyzzyaabn1/=0.d0)then
if(dimensionality==2)then
xyzzyaaag3=sqrt(sum(a1(1:2)**2))
xyzzyaaah3=sqrt(sum(a2(1:2)**2))
xyzzyaaai3=abs(dot_product(a1(1:2),a2(1:2)))/(xyzzyaaag3*xyzzyaaah3)
if(abs(xyzzyaaag3-xyzzyaaah3)<xyzzyaaak3*xyzzyaaag3)then
if(abs(xyzzyaaai3-0.5d0)<xyzzyaaak3)then
xyzzyaaaj1=86
elseif(abs(xyzzyaaai3)<xyzzyaaak3)then
xyzzyaaaj1=119
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
if(xyzzyaaaj1+1>size(first_s_in_star))call errstop('INIT_MAN_INT','NST&
&AR_R_TO_MINUS3 is too large (1).')
if(first_s_in_star(xyzzyaaaj1+1)==0)call errstop('INIT_MAN_INT','NSTAR&
&_R_TO_MINUS3 is too large (2).')
xyzzyaaab3=1
do xyzzyaaac3=2,xyzzyaaaj1
xyzzyaaad3=first_s_in_star(xyzzyaaac3)
xyzzyaaae3=first_s_in_star(xyzzyaaac3+1)-xyzzyaaad3
xyzzyaaab3=xyzzyaaab3+xyzzyaaae3
self_term=self_term+dble(xyzzyaaae3)*sum(s_lattice(1:dimensionality,xy&
&zzyaaad3)**2)**(-1.5d0)
enddo
xyzzyaaaf3=sqrt(dble(xyzzyaaab3)*area*one_over_pi)
xyzzyaabo1=twopi/(area*xyzzyaaaf3)
xyzzyaabp1=0.75/xyzzyaaaf3**2
self_term=0.5d0*xyzzyaabn1*(self_term+xyzzyaabo1)
else
call errstop_master('INIT_MAN_INT','r^(-3) tail is only implemented fo&
&r 2D case at present.')
endif
endif
if(am_master)then
call wout('Interaction type: hard-sphere potential (manual)')
call wout('Hard-sphere diameter                     : ',hard_diam)
if(hard_op_spins)call wout('Opposite spins only')
endif
if(hard_diam<0.d0)call errstop_master('INIT_MAN_INT','Hard-sphere radi&
&us should be non-negative.')
if(am_master)then
if(xyzzyaabn1/=0.d0)then
call wout('Weight of r^(-3) tail                    : ',xyzzyaabn1)
call wout('Number of stars of lattice vectors in sum: '//trim(i2s(xyzz&
&yaaaj1)))
call wout('Radius of cutoff sphere                  : ',xyzzyaaaf3)
call wout('Self-term v_M                            : ',2.d0*self_term&
&)
endif
call wout
endif
iinteraction=xyzzyaacu1
case('polynomial')
if(am_master)then
call wout("Interaction type: polynomial potential (manual)")
call wout("Cutoff : ",man_int_poly_cutoff)
call wout("Non-zero coefficients : ")
do xyzzyaaaa3=0,man_int_poly_order
if(man_int_poly_coeffs(xyzzyaaaa3)/=0.d0)then
call wout("  r^"//trim(i2s(xyzzyaaaa3))//" : ",man_int_poly_coeffs(xyz&
&zyaaaa3))
endif
enddo
call wout
endif
iinteraction=xyzzyaacv1
case("logarithmic")
rstar=man_int_params(1)
if(am_master)then
call wout("Interaction type: logarithmic (manual)")
call wout("rstar : ",rstar)
call wout
endif
if(rstar<=0.d0)call errstop_master('INIT_MAN_INT','Should have rstar>0&
& in logarithmic interaction.')
if(dimensionality/=2)call errstop_master('INIT_MAN_INT','System should&
& be 2D for logarithmic interaction.')
if(periodicity/=0)call errstop_master('INIT_MAN_INT','System should be&
& aperiodic for logarithmic interaction.')
iinteraction=ilogarithmic
case("2D_int")
rstar=man_int_params(1)
if(am_master)then
call wout("Interaction type: 2D Keldysh interaction (manual)")
call wout("rstar : ",rstar)
if(xyzzyaadf1)call wordwrap("A simple approximation to the 2D interact&
&ion will be used, because the SIMPLE_APPROX_2D_INT parameter is set t&
&o T in INTERACTIONS.")
call wout
endif
if(rstar<=0.d0)call errstop_master('INIT_MAN_INT','Should have rstar>0&
& in 2D interaction.')
if(dimensionality/=2)call errstop_master('INIT_MAN_INT','System should&
& be, er, 2D for 2D interaction.')
if(periodicity/=0)call errstop_master('INIT_MAN_INT','System should be&
& aperiodic for 2D interaction.')
iinteraction=i2d_int
case('dipole')
xyzzyaabq1=man_int_params(1)
if(dimensionality==2)then
xyzzyaaag3=sqrt(sum(a1(1:2)**2))
xyzzyaaah3=sqrt(sum(a2(1:2)**2))
xyzzyaaai3=abs(dot_product(a1(1:2),a2(1:2)))/(xyzzyaaag3*xyzzyaaah3)
if(abs(xyzzyaaag3-xyzzyaaah3)<xyzzyaaak3*xyzzyaaag3)then
if(abs(xyzzyaaai3-0.5d0)<xyzzyaaak3)then
xyzzyaaal1=86
elseif(abs(xyzzyaaai3)<xyzzyaaak3)then
xyzzyaaal1=119
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
if(xyzzyaaal1+1>size(first_s_in_star))call errstop_master('INIT_MAN_IN&
&T','NSTAR_DIPOLE is too large (1).')
if(first_s_in_star(xyzzyaaal1+1)==0)call errstop_master('INIT_MAN_INT'&
&,'NSTAR_DIPOLE is too large (2).')
xyzzyaaab3=1
do xyzzyaaac3=2,xyzzyaaal1
xyzzyaaad3=first_s_in_star(xyzzyaaac3)
xyzzyaaae3=first_s_in_star(xyzzyaaac3+1)-xyzzyaaad3
xyzzyaaab3=xyzzyaaab3+xyzzyaaae3
self_term=self_term+dble(xyzzyaaae3)*sum(s_lattice(1:dimensionality,xy&
&zzyaaad3)**2)**(-1.5d0)
enddo
xyzzyaaaf3=sqrt(dble(xyzzyaaab3)*area*one_over_pi)
xyzzyaabo1=twopi/(area*xyzzyaaaf3)
xyzzyaabp1=0.75/xyzzyaaaf3**2
self_term=0.5d0*xyzzyaabq1*(self_term+xyzzyaabo1)
else
call errstop_master('INIT_MAN_INT','Dipole interaction is only impleme&
&nted for 2D case at present.')
endif
if(am_master)call wout('Interaction type: dipolar potential (manual)')
if(xyzzyaabq1<0.d0)call errstop_master('INIT_MAN_INT','Give dipole str&
&ength as d^2.')
if(am_master)then
call wout('Dipolar interaction strength             : ',xyzzyaabq1)
call wout('Number of stars of lattice vectors in sum: '//trim(i2s(xyzz&
&yaaal1)))
call wout('Self-term v_M                            : ',2.d0*self_term&
&)
endif
iinteraction=xyzzyaacw1
case('pseudodipole')
if(man_int_poly_cutoff>wigner_seitz_radius)call errstop_master('INIT_M&
&AN_INT','Cutoff radius should be less than the Wigner-Seitz radius.')
xyzzyaabq1=man_int_params(1)
if(dimensionality==2)then
xyzzyaaag3=sqrt(sum(a1(1:2)**2))
xyzzyaaah3=sqrt(sum(a2(1:2)**2))
xyzzyaaai3=abs(dot_product(a1(1:2),a2(1:2)))/(xyzzyaaag3*xyzzyaaah3)
if(abs(xyzzyaaag3-xyzzyaaah3)<xyzzyaaak3*xyzzyaaag3)then
if(abs(xyzzyaaai3-0.5d0)<xyzzyaaak3)then
xyzzyaaal1=86
elseif(abs(xyzzyaaai3)<xyzzyaaak3)then
xyzzyaaal1=119
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
if(xyzzyaaal1+1>size(first_s_in_star))call errstop_master('INIT_MAN_IN&
&T','NSTAR_DIPOLE is too large (1).')
if(first_s_in_star(xyzzyaaal1+1)==0)call errstop_master('INIT_MAN_INT'&
&,'NSTAR_DIPOLE is too large (2).')
xyzzyaaab3=1
do xyzzyaaac3=2,xyzzyaaal1
xyzzyaaad3=first_s_in_star(xyzzyaaac3)
xyzzyaaae3=first_s_in_star(xyzzyaaac3+1)-xyzzyaaad3
xyzzyaaab3=xyzzyaaab3+xyzzyaaae3
self_term=self_term+dble(xyzzyaaae3)*sum(s_lattice(1:dimensionality,xy&
&zzyaaad3)**2)**(-1.5d0)
enddo
xyzzyaaaf3=sqrt(dble(xyzzyaaab3)*area*one_over_pi)
xyzzyaabo1=twopi/(area*xyzzyaaaf3)
xyzzyaabp1=0.75/xyzzyaaaf3**2
self_term=0.5d0*xyzzyaabq1*(self_term+xyzzyaabo1)
else
call errstop_master('INIT_MAN_INT','Pseudodipole interaction is only i&
&mplemented for 2D case at present.')
endif
if(am_master)call wout('Interaction type: pseudodipolar potential (man&
&ual)')
if(xyzzyaabq1<0.d0)call errstop_master('INIT_MAN_INT','Give dipole str&
&ength as d^2.')
if(am_master)then
call wout('Dipolar interaction strength             : ',xyzzyaabq1)
call wout("Cutoff : ",man_int_poly_cutoff)
call wout("Non-zero coefficients : ")
do xyzzyaaaa3=0,man_int_poly_order
if(man_int_poly_coeffs(xyzzyaaaa3)/=0.d0)call wout("  r^"//trim(i2s(xy&
&zzyaaaa3))//" : ",man_int_poly_coeffs(xyzzyaaaa3))
enddo
call wout('Number of stars of lattice vectors in sum: '//trim(i2s(xyzz&
&yaaal1)))
call wout('Self-term v_M                            : ',2.d0*self_term&
&)
endif
iinteraction=xyzzyaacx1
case('2D_tilted_dipole')
xyzzyaabq1=man_int_params(1)
xyzzyaabr1=man_int_params(2)
if(dimensionality==2)then
xyzzyaaag3=sqrt(sum(a1(1:2)**2))
xyzzyaaah3=sqrt(sum(a2(1:2)**2))
xyzzyaaai3=abs(dot_product(a1(1:2),a2(1:2)))/(xyzzyaaag3*xyzzyaaah3)
if(abs(xyzzyaaag3-xyzzyaaah3)<xyzzyaaak3*xyzzyaaag3)then
if(abs(xyzzyaaai3-0.5d0)<xyzzyaaak3)then
xyzzyaaal1=86
elseif(abs(xyzzyaaai3)<xyzzyaaak3)then
xyzzyaaal1=119
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
if(xyzzyaaal1+1>size(first_s_in_star))call errstop_master('INIT_MAN_IN&
&T','NSTAR_DIPOLE is too large (1).')
if(first_s_in_star(xyzzyaaal1+1)==0)call errstop_master('INIT_MAN_INT'&
&,'NSTAR_DIPOLE is too large (2).')
self_term=0.d0
else
call errstop_master('INIT_MAN_INT','Tilted dipole interaction is only &
&implemented for 2D case at present.')
endif
if(am_master)call wout('Interaction type: tilted dipolar potential (ma&
&nual)')
if(xyzzyaabq1<0.d0)call errstop_master('INIT_MAN_INT','Give dipole str&
&ength as d^2.')
if(xyzzyaabr1<0.d0)call errstop_master('INIT_MAN_INT','The tilt angle &
&should be between 0 and pi.')
if(xyzzyaabr1>pi)call errstop_master('INIT_MAN_INT','The tilt angle sh&
&ould be between 0 and pi.')
if(am_master)then
call wout('Dipolar interaction strength             : ',xyzzyaabq1)
call wout('Tilt angle away from vertical            : ',xyzzyaabr1)
call wout('Number of stars of lattice vectors in sum: '//trim(i2s(xyzz&
&yaaal1)))
endif
iinteraction=xyzzyaacy1
case('2D_tilt_pseudodipole')
if(man_int_poly_cutoff>wigner_seitz_radius)call errstop_master('INIT_M&
&AN_INT','Cutoff radius should be less than the Wigner-Seitz radius.')
xyzzyaabq1=man_int_params(1)
xyzzyaabr1=man_int_params(2)
if(dimensionality==2)then
xyzzyaaag3=sqrt(sum(a1(1:2)**2))
xyzzyaaah3=sqrt(sum(a2(1:2)**2))
xyzzyaaai3=abs(dot_product(a1(1:2),a2(1:2)))/(xyzzyaaag3*xyzzyaaah3)
if(abs(xyzzyaaag3-xyzzyaaah3)<xyzzyaaak3*xyzzyaaag3)then
if(abs(xyzzyaaai3-0.5d0)<xyzzyaaak3)then
xyzzyaaal1=86
elseif(abs(xyzzyaaai3)<xyzzyaaak3)then
xyzzyaaal1=119
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
else
call errstop_master('INIT_MAN_INT','You might want to set up a new def&
&ault behaviour for the number of stars of lattice vectors to use when&
& evaluating the interaction.')
endif
if(xyzzyaaal1+1>size(first_s_in_star))call errstop_master('INIT_MAN_IN&
&T','NSTAR_DIPOLE is too large (1).')
if(first_s_in_star(xyzzyaaal1+1)==0)call errstop_master('INIT_MAN_INT'&
&,'NSTAR_DIPOLE is too large (2).')
else
call errstop_master('INIT_MAN_INT','Tilted pseudodipole interaction is&
& only implemented for 2D case at present.')
endif
if(am_master)call wout('Interaction type: tilted pseudodipolar potenti&
&al (manual)')
if(xyzzyaabq1<0.d0)call errstop_master('INIT_MAN_INT','Give dipole str&
&ength as d^2.')
if(xyzzyaabr1<0.d0)call errstop_master('INIT_MAN_INT','The tilt angle &
&should be between 0 and pi.')
if(xyzzyaabr1>pi)call errstop_master('INIT_MAN_INT','The tilt angle sh&
&ould be between 0 and pi.')
if(am_master)then
call wout('Dipolar interaction strength             : ',xyzzyaabq1)
call wout('Tilt angle away from vertical            : ',xyzzyaabr1)
call wout("Cutoff : ",man_int_poly_cutoff)
call wout("Non-zero coefficients : ")
do xyzzyaaaa3=0,man_int_poly_order
if(man_int_poly_coeffs(xyzzyaaaa3)/=0.d0)call wout("  r^"//trim(i2s(xy&
&zzyaaaa3))//" : ",man_int_poly_coeffs(xyzzyaaaa3))
enddo
call wout("Non-zero coefficients of cos(2 phi) : ")
do xyzzyaaaa3=0,man_int_tilt_poly_order
if(man_int_tilt_poly_coeffs(xyzzyaaaa3)/=0.d0)call wout("  r^"//trim(i&
&2s(xyzzyaaaa3))//" : ",man_int_tilt_poly_coeffs(xyzzyaaaa3))
enddo
call wout('Number of stars of lattice vectors in sum: '//trim(i2s(xyzz&
&yaaal1)))
endif
iinteraction=xyzzyaacz1
case default
call errstop_master('INIT_MAN_INT','Unknown manual interaction '//trim&
&(int_name))
end select
end subroutine xyzzyaadg1
subroutine en_potential(n,z,vecji,enpot,field_at_ion)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: z(n)
real(dp),intent(in) :: vecji(4,n)
real(dp),intent(out) :: enpot,field_at_ion(3,n)
integer xyzzyaaaa4
real(dp) xyzzyaaab4,xyzzyaaac4(4)
enpot=0.d0
select case(periodicity)
case(3)
if(calc_field)then
do xyzzyaaaa4=1,n
xyzzyaaac4=vecji(1:4,xyzzyaaaa4)
call ewald_3d(1,xyzzyaaac4,xyzzyaaab4,.false.,field_at_ion(1,xyzzyaaaa&
&4))
enpot=enpot-z(xyzzyaaaa4)*xyzzyaaab4
enddo
else
do xyzzyaaaa4=1,n
xyzzyaaac4=vecji(1:4,xyzzyaaaa4)
call ewald_3d(1,xyzzyaaac4,xyzzyaaab4,.false.)
enpot=enpot-z(xyzzyaaaa4)*xyzzyaaab4
enddo
endif
case(2)
if(calc_field)then
do xyzzyaaaa4=1,n
xyzzyaaac4=vecji(1:4,xyzzyaaaa4)
call ewald_2d(1,xyzzyaaac4,xyzzyaaab4,field_at_ion(1,xyzzyaaaa4))
enpot=enpot-z(xyzzyaaaa4)*xyzzyaaab4
enddo
else
do xyzzyaaaa4=1,n
xyzzyaaac4=vecji(1:4,xyzzyaaaa4)
call ewald_2d(1,xyzzyaaac4,xyzzyaaab4)
enpot=enpot-z(xyzzyaaaa4)*xyzzyaaab4
enddo
endif
case(1)
do xyzzyaaaa4=1,n
xyzzyaaac4=vecji(1:4,xyzzyaaaa4)
call xyzzyaadi1(0,1,xyzzyaaac4,xyzzyaaab4)
enpot=enpot-z(xyzzyaaaa4)*xyzzyaaab4
enddo
case default
call errstop('EN_POTENTIAL','Invalid periodicity - this is a bug.')
end select
end subroutine en_potential
subroutine ee_potential(ie,b,vecji,eepot,just_want_j_gt_i,isinf)
implicit none
integer,intent(in) :: ie,b
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: eepot(3)
logical,intent(in) :: just_want_j_gt_i
logical,intent(out) :: isinf
call timer('EE_POTENTIAL',.true.)
isinf=.false.
select case(iinteraction)
case(inone)
eepot=0.d0
case(icoulomb)
call xyzzyaaec1(ie,vecji,eepot(1),just_want_j_gt_i,isinf)
eepot(2:3)=0.d0
case(xyzzyaacp1)
call xyzzyaadh1(ie,vecji,eepot(1),just_want_j_gt_i,pseudo_ewald)
eepot(2:3)=0.d0
case(xyzzyaacq1)
call xyzzyaaeb1(ie,b,vecji,eepot(2),eepot(3),just_want_j_gt_i,isinf)
eepot(1)=0.d0
case(xyzzyaacr1)
call xyzzyaadh1(ie,vecji,eepot(1),just_want_j_gt_i,pseudo_ewald)
call xyzzyaaeb1(ie,b,vecji,eepot(2),eepot(3),just_want_j_gt_i,isinf)
case(xyzzyaacs1)
call xyzzyaadl1(ie,vecji,man_int_params(1),man_int_params(2),eepot(1),&
&just_want_j_gt_i)
eepot(2:3)=0.d0
case(xyzzyaact1)
call xyzzyaadm1(ie,vecji,eepot(1),just_want_j_gt_i)
eepot(2:3)=0.d0
case(xyzzyaacu1)
call xyzzyaado1(ie,vecji,eepot(1),just_want_j_gt_i)
eepot(2:3)=0.d0
case(xyzzyaacv1)
call xyzzyaadr1(ie,vecji,eepot(1),just_want_j_gt_i)
eepot(2:3)=0.d0
case(ilogarithmic)
call xyzzyaads1(ie,vecji,eepot(1),just_want_j_gt_i,isinf)
eepot(2:3)=0.d0
case(i2d_int)
call xyzzyaadt1(ie,vecji,eepot(1),just_want_j_gt_i,isinf)
eepot(2:3)=0.d0
case(xyzzyaacw1)
call xyzzyaadu1(ie,vecji,eepot(1),just_want_j_gt_i)
eepot(2:3)=0.d0
case(xyzzyaacx1)
call xyzzyaadv1(ie,vecji,eepot(1),just_want_j_gt_i)
eepot(2:3)=0.d0
case(xyzzyaacy1)
call xyzzyaadx1(ie,vecji,eepot(1),just_want_j_gt_i)
eepot(2:3)=0.d0
case(xyzzyaacz1)
call xyzzyaadw1(ie,vecji,eepot(1),just_want_j_gt_i)
eepot(2:3)=0.d0
case default
call errstop('EE_POTENTIAL','Puzzled.')
end select
call timer('EE_POTENTIAL',.false.)
end subroutine ee_potential
subroutine xyzzyaadh1(ie,vecji,ewald,just_want_j_gt_i,pseudo)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: ewald
logical,intent(in) :: just_want_j_gt_i,pseudo
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6,xyzzyaaad6,xyzzyaaae6
real(dp) xyzzyaaaf6,xyzzyaaag6
call timer('EWALD',.true.)
select case(periodicity)
case(3)
if(electron_system)then
if(just_want_j_gt_i)then
if(ie<netot)then
call ewald_3d(netot-ie,vecji(1,ie+1),xyzzyaaag6,pseudo)
ewald=xyzzyaaag6+self_term
else
ewald=self_term
endif
else
if(ie>1)then
call ewald_3d(ie-1,vecji(1,1),xyzzyaaaf6,pseudo)
else
xyzzyaaaf6=0.d0
endif
if(ie<netot)then
call ewald_3d(netot-ie,vecji(1,ie+1),xyzzyaaag6,pseudo)
else
xyzzyaaag6=0.d0
endif
ewald=0.5d0*(xyzzyaaaf6+xyzzyaaag6)+self_term
endif
else
xyzzyaaad6=0
do xyzzyaaaa6=1,nspin
xyzzyaaad6=xyzzyaaad6+nele(xyzzyaaaa6)
if(ie<=xyzzyaaad6)exit
enddo
xyzzyaaac6=ie-xyzzyaaad6+nele(xyzzyaaaa6)
if(just_want_j_gt_i)then
if(xyzzyaaac6<nele(xyzzyaaaa6))then
call ewald_3d(nele(xyzzyaaaa6)-xyzzyaaac6,vecji(1,ie+1),xyzzyaaag6,pse&
&udo)
ewald=(xyzzyaaag6+self_term)*xyzzyaacd1(xyzzyaaaa6,xyzzyaaaa6)
else
ewald=self_term*xyzzyaacd1(xyzzyaaaa6,xyzzyaaaa6)
endif
xyzzyaaae6=xyzzyaaad6+1
do xyzzyaaab6=xyzzyaaaa6+1,nspin
call ewald_3d(nele(xyzzyaaab6),vecji(1,xyzzyaaae6),xyzzyaaag6,pseudo)
ewald=ewald+xyzzyaaag6*xyzzyaacd1(xyzzyaaab6,xyzzyaaaa6)
xyzzyaaae6=xyzzyaaae6+nele(xyzzyaaab6)
enddo
else
ewald=0.d0
xyzzyaaae6=1
do xyzzyaaab6=1,nspin
if(xyzzyaaab6==xyzzyaaaa6)then
if(xyzzyaaac6>1)then
call ewald_3d(xyzzyaaac6-1,vecji(1,xyzzyaaae6),xyzzyaaaf6,pseudo)
else
xyzzyaaaf6=0.d0
endif
if(xyzzyaaac6<nele(xyzzyaaaa6))then
call ewald_3d(nele(xyzzyaaaa6)-xyzzyaaac6,vecji(1,ie+1),xyzzyaaag6,pse&
&udo)
else
xyzzyaaag6=0.d0
endif
ewald=ewald+(0.5d0*(xyzzyaaaf6+xyzzyaaag6)+self_term)*xyzzyaacd1(xyzzy&
&aaaa6,xyzzyaaaa6)
else
call ewald_3d(nele(xyzzyaaab6),vecji(1,xyzzyaaae6),xyzzyaaag6,pseudo)
ewald=ewald+0.5d0*xyzzyaaag6*xyzzyaacd1(xyzzyaaab6,xyzzyaaaa6)
endif
xyzzyaaae6=xyzzyaaae6+nele(xyzzyaaab6)
enddo
endif
endif
case(2)
if(electron_system)then
if(just_want_j_gt_i)then
if(ie<netot)then
call ewald_2d(netot-ie,vecji(1,ie+1),xyzzyaaag6)
ewald=xyzzyaaag6+self_term
else
ewald=self_term
endif
else
if(ie>1)then
call ewald_2d(ie-1,vecji(1,1),xyzzyaaaf6)
else
xyzzyaaaf6=0.d0
endif
if(ie<netot)then
call ewald_2d(netot-ie,vecji(1,ie+1),xyzzyaaag6)
else
xyzzyaaag6=0.d0
endif
ewald=0.5d0*(xyzzyaaaf6+xyzzyaaag6)+self_term
endif
else
xyzzyaaad6=0
do xyzzyaaaa6=1,nspin
xyzzyaaad6=xyzzyaaad6+nele(xyzzyaaaa6)
if(ie<=xyzzyaaad6)exit
enddo
xyzzyaaac6=ie-xyzzyaaad6+nele(xyzzyaaaa6)
if(just_want_j_gt_i)then
if(xyzzyaaac6<nele(xyzzyaaaa6))then
call ewald_2d(nele(xyzzyaaaa6)-xyzzyaaac6,vecji(1,ie+1),xyzzyaaag6)
ewald=(xyzzyaaag6+self_term)*xyzzyaacd1(xyzzyaaaa6,xyzzyaaaa6)
else
ewald=self_term*xyzzyaacd1(xyzzyaaaa6,xyzzyaaaa6)
endif
xyzzyaaae6=xyzzyaaad6+1
do xyzzyaaab6=xyzzyaaaa6+1,nspin
call ewald_2d(nele(xyzzyaaab6),vecji(1,xyzzyaaae6),xyzzyaaag6)
ewald=ewald+xyzzyaaag6*xyzzyaacd1(xyzzyaaab6,xyzzyaaaa6)
xyzzyaaae6=xyzzyaaae6+nele(xyzzyaaab6)
enddo
else
ewald=0.d0
xyzzyaaae6=1
do xyzzyaaab6=1,nspin
if(xyzzyaaab6==xyzzyaaaa6)then
if(xyzzyaaac6>1)then
call ewald_2d(xyzzyaaac6-1,vecji(1,xyzzyaaae6),xyzzyaaaf6)
else
xyzzyaaaf6=0.d0
endif
if(xyzzyaaac6<nele(xyzzyaaaa6))then
call ewald_2d(nele(xyzzyaaaa6)-xyzzyaaac6,vecji(1,ie+1),xyzzyaaag6)
else
xyzzyaaag6=0.d0
endif
ewald=ewald+(0.5d0*(xyzzyaaaf6+xyzzyaaag6)+self_term)*xyzzyaacd1(xyzzy&
&aaaa6,xyzzyaaaa6)
else
call ewald_2d(nele(xyzzyaaab6),vecji(1,xyzzyaaae6),xyzzyaaag6)
ewald=ewald+0.5d0*xyzzyaaag6*xyzzyaacd1(xyzzyaaab6,xyzzyaaaa6)
endif
xyzzyaaae6=xyzzyaaae6+nele(xyzzyaaab6)
enddo
endif
endif
case(1)
if(electron_system)then
if(just_want_j_gt_i)then
if(ie<netot)then
call xyzzyaadi1(0,netot-ie,vecji(1,ie+1),xyzzyaaag6)
ewald=xyzzyaaag6+self_term
else
ewald=self_term
endif
else
call xyzzyaadi1(ie,netot,vecji,xyzzyaaag6)
ewald=0.5d0*xyzzyaaag6+self_term
endif
else
xyzzyaaad6=0
do xyzzyaaaa6=1,nspin
xyzzyaaad6=xyzzyaaad6+nele(xyzzyaaaa6)
if(ie<=xyzzyaaad6)exit
enddo
xyzzyaaac6=ie-xyzzyaaad6+nele(xyzzyaaaa6)
if(just_want_j_gt_i)then
if(xyzzyaaac6<nele(xyzzyaaaa6))then
call xyzzyaadi1(0,nele(xyzzyaaaa6)-xyzzyaaac6,vecji(1,ie+1),xyzzyaaag6&
&)
ewald=(xyzzyaaag6+self_term)*xyzzyaacd1(xyzzyaaaa6,xyzzyaaaa6)
else
ewald=self_term*xyzzyaacd1(xyzzyaaaa6,xyzzyaaaa6)
endif
xyzzyaaae6=xyzzyaaad6+1
do xyzzyaaab6=xyzzyaaaa6+1,nspin
call xyzzyaadi1(0,nele(xyzzyaaab6),vecji(1,xyzzyaaae6),xyzzyaaag6)
ewald=ewald+xyzzyaaag6*xyzzyaacd1(xyzzyaaab6,xyzzyaaaa6)
xyzzyaaae6=xyzzyaaae6+nele(xyzzyaaab6)
enddo
else
ewald=0.d0
xyzzyaaae6=1
do xyzzyaaab6=1,nspin
if(xyzzyaaab6==xyzzyaaaa6)then
call xyzzyaadi1(xyzzyaaac6,nele(xyzzyaaab6),vecji(1,xyzzyaaae6),xyzzya&
&aag6)
ewald=ewald+(0.5d0*xyzzyaaag6+self_term)*xyzzyaacd1(xyzzyaaab6,xyzzyaa&
&ab6)
else
call xyzzyaadi1(0,nele(xyzzyaaab6),vecji(1,xyzzyaaae6),xyzzyaaag6)
ewald=ewald+0.5d0*xyzzyaaag6*xyzzyaacd1(xyzzyaaab6,xyzzyaaaa6)
endif
xyzzyaaae6=xyzzyaaae6+nele(xyzzyaaab6)
enddo
endif
endif
case default
call errstop('EE_EVAL_EWALD','Invalid periodicity - this is a bug.')
end select
call timer('EWALD',.false.)
end subroutine xyzzyaadh1
subroutine xyzzyaadi1(ie,n,vecji,pote)
use slaarnaag,only : root_pi
use slaarnaas,only : harmwire_b
use slaarnabi,only : setup_spline_reg,eval_spline_reg
use slaarnabt,only : erfc,exp_int
implicit none
integer,intent(in) :: ie,n
real(dp),intent(in) :: vecji(4,n)
real(dp),intent(out) :: pote
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7,xyzzyaaad7,xyzzyaaae7,xyzzyaa&
&af7
integer,save :: xyzzyaaag7
real(dp) xyzzyaaah7(3),xyzzyaaai7,xyzzyaaaj7,xyzzyaaak7,xyzzyaaal7,xyz&
&zyaaam7,xyzzyaaan7,xyzzyaaao7,xyzzyaaap7,xyzzyaaaq7,xyzzyaaar7,xyzzya&
&aas7,xyzzyaaat7,xyzzyaaau7,xyzzyaaav7,xyzzyaaaw7,xyzzyaaax7,xyzzyaaay&
&7,xyzzyaaaz7,xyzzyaaba7,xyzzyaabb7,xyzzyaabc7,xyzzyaabd7,xyzzyaabe7,x&
&yzzyaabf7,xyzzyaabg7,xyzzyaabh7,xyzzyaabi7
real(dp),save :: xyzzyaabj7
real(dp),parameter :: xyzzyaabk7=1.d-7
real(dp),allocatable,save :: xyzzyaabl7(:),xyzzyaabm7(:),xyzzyaabn7(:)
logical,save :: xyzzyaabo7=.true.
if(harmwire_b>0.d0)then
if(xyzzyaabo7)then
xyzzyaabo7=.false.
xyzzyaabg7=harmwire_b*harmwire_b
xyzzyaabh7=0.5d0/harmwire_b
xyzzyaabi7=root_pi*2.d0*xyzzyaabh7
xyzzyaaag7=netot*50
xyzzyaabj7=wigner_seitz_radius*2.d0
allocate(xyzzyaabl7(xyzzyaaag7),xyzzyaabm7(xyzzyaaag7),xyzzyaabn7(xyzz&
&yaaag7),stat=xyzzyaaaf7)
call check_alloc(xyzzyaaaf7,'EWALD_1D','')
xyzzyaabn7(1)=0.d0
xyzzyaabl7(1)=2.d0*xyzzyaadj1(harmwire_b)+xyzzyaabi7
do xyzzyaaad7=2,xyzzyaaag7
xyzzyaaay7=dble(xyzzyaaad7-1)*xyzzyaabj7/dble(2*(xyzzyaaag7-1))
xyzzyaabn7(xyzzyaaad7)=xyzzyaaay7
xyzzyaaav7=xyzzyaabi7
xyzzyaaaw7=0.d0
do xyzzyaaae7=-1,1,2
xyzzyaaab7=xyzzyaaae7
xyzzyaaba7=0.d0
do
xyzzyaaaz7=abs(xyzzyaaay7-dble(xyzzyaaab7)*xyzzyaabj7)*xyzzyaabh7
xyzzyaaaz7=xyzzyaadk1(xyzzyaaaz7)
xyzzyaaba7=xyzzyaaba7+xyzzyaaaz7
if(xyzzyaaav7*abs(xyzzyaaaz7)<xyzzyaabk7)exit
xyzzyaaab7=xyzzyaaab7+xyzzyaaae7
enddo
xyzzyaaaw7=xyzzyaaaw7+xyzzyaaav7*xyzzyaaba7
enddo
xyzzyaaaz7=abs(xyzzyaaay7)*xyzzyaabh7
xyzzyaaaz7=xyzzyaadk1(xyzzyaaaz7)
xyzzyaaaw7=xyzzyaaaw7+xyzzyaaaz7*xyzzyaaav7
xyzzyaaba7=0.d0
xyzzyaaaz7=0.d0
xyzzyaaav7=-2.d0
do xyzzyaaae7=-1,1,2
xyzzyaaab7=xyzzyaaae7
do
xyzzyaaaz7=(1.d0/abs(xyzzyaaay7-dble(xyzzyaaab7)*xyzzyaabj7))*(1.d0-er&
&fc((abs(xyzzyaaay7-dble(xyzzyaaab7)*xyzzyaabj7)*xyzzyaabh7)))
xyzzyaaba7=xyzzyaaba7+xyzzyaaaz7
if(abs(xyzzyaaav7*xyzzyaaaz7)<xyzzyaabk7)exit
xyzzyaaab7=xyzzyaaab7+xyzzyaaae7
enddo
enddo
xyzzyaaaw7=xyzzyaaaw7+xyzzyaaav7*xyzzyaaba7
if(xyzzyaaay7==0.d0)then
xyzzyaaaz7=one_over_root_pi/harmwire_b
else
xyzzyaaaz7=(1.d0/abs(xyzzyaaay7))*(1.d0-erfc(abs(xyzzyaaay7)*xyzzyaabh&
&7))
endif
xyzzyaaaw7=xyzzyaaaw7+xyzzyaaav7*xyzzyaaaz7
xyzzyaaax7=0.d0
xyzzyaaab7=1
xyzzyaabe7=twopi/xyzzyaabj7
xyzzyaaav7=4.d0/xyzzyaabj7
do
xyzzyaabf7=xyzzyaabe7*dble(xyzzyaaab7)
xyzzyaaaz7=cos(xyzzyaabf7*xyzzyaaay7)*exp_int(xyzzyaabg7*xyzzyaabf7*xy&
&zzyaabf7)
xyzzyaaax7=xyzzyaaax7+xyzzyaaaz7
if(xyzzyaabf7*harmwire_b>12.d0)exit
xyzzyaaab7=xyzzyaaab7+1
enddo
xyzzyaaax7=xyzzyaaax7*xyzzyaaav7
xyzzyaabl7(xyzzyaaad7)=xyzzyaaaw7+xyzzyaaax7
enddo
xyzzyaabl7=0.5d0*xyzzyaabl7
call setup_spline_reg(xyzzyaaag7,xyzzyaabn7,xyzzyaabl7,xyzzyaabm7)
endif
pote=0.d0
do xyzzyaaac7=1,n
xyzzyaaay7=abs(vecji(1,xyzzyaaac7))
if(xyzzyaaay7>0.5d0*xyzzyaabj7)call errstop('EWALD_1D','Attempting to &
&evaluate the quasi-1d interaction beyond the extent of the spline gri&
&d. Bug.')
call eval_spline_reg(xyzzyaaag7,xyzzyaabn7,xyzzyaabl7,xyzzyaabm7,xyzzy&
&aaay7,.true.,.false.,xyzzyaabd7,xyzzyaabb7,xyzzyaabc7)
pote=pote+xyzzyaabd7
enddo
else
pote=0.d0
do xyzzyaaaa7=1,ie-1
xyzzyaaah7=vecji(1:3,xyzzyaaaa7)
xyzzyaaai7=xyzzyaaah7(2)*xyzzyaaah7(2)+xyzzyaaah7(3)*xyzzyaaah7(3)
xyzzyaaak7=xyzzyaaah7(1)
xyzzyaaaj7=0.d0
do xyzzyaaab7=1,xyzzyaaaa1
xyzzyaaaj7=1.d0/sqrt((xyzzyaaak7+s_lattice(1,xyzzyaaab7))**2+xyzzyaaai&
&7)+xyzzyaaaj7
enddo
xyzzyaaal7=xyzzyaaaq1+xyzzyaaak7
xyzzyaaam7=xyzzyaaaq1-xyzzyaaak7
xyzzyaaan7=xyzzyaaal7*xyzzyaaal7
xyzzyaaao7=xyzzyaaam7*xyzzyaaam7
xyzzyaaap7=sqrt(xyzzyaaan7+xyzzyaaai7)
xyzzyaaaq7=sqrt(xyzzyaaao7+xyzzyaaai7)
xyzzyaaar7=1.d0/xyzzyaaap7
xyzzyaaas7=1.d0/xyzzyaaaq7
xyzzyaaat7=xyzzyaaar7*xyzzyaaar7
xyzzyaaau7=xyzzyaaas7*xyzzyaaas7
xyzzyaaar7=xyzzyaaal7*xyzzyaaar7*xyzzyaaat7
xyzzyaaas7=xyzzyaaam7*xyzzyaaas7*xyzzyaaau7
pote=(log((xyzzyaaap7+xyzzyaaal7)*(xyzzyaaaq7+xyzzyaaam7))*xyzzyaabf1-&
&(xyzzyaaar7+xyzzyaaas7)*xyzzyaabh1+((9.d0-xyzzyaaat7*xyzzyaaan7*15.d0&
&)*xyzzyaaar7*xyzzyaaat7+(9.d0-xyzzyaaau7*xyzzyaaao7*15.d0)*xyzzyaaas7&
&*xyzzyaaau7)*xyzzyaaap1+xyzzyaaaj7)+pote
enddo
do xyzzyaaaa7=ie+1,n
xyzzyaaah7=vecji(1:3,xyzzyaaaa7)
xyzzyaaai7=xyzzyaaah7(2)*xyzzyaaah7(2)+xyzzyaaah7(3)*xyzzyaaah7(3)
xyzzyaaak7=xyzzyaaah7(1)
xyzzyaaaj7=0.d0
do xyzzyaaab7=1,xyzzyaaaa1
xyzzyaaaj7=1.d0/sqrt((xyzzyaaak7+s_lattice(1,xyzzyaaab7))**2+xyzzyaaai&
&7)+xyzzyaaaj7
enddo
xyzzyaaal7=xyzzyaaaq1+xyzzyaaak7
xyzzyaaam7=xyzzyaaaq1-xyzzyaaak7
xyzzyaaan7=xyzzyaaal7*xyzzyaaal7
xyzzyaaao7=xyzzyaaam7*xyzzyaaam7
xyzzyaaap7=sqrt(xyzzyaaan7+xyzzyaaai7)
xyzzyaaaq7=sqrt(xyzzyaaao7+xyzzyaaai7)
xyzzyaaar7=1.d0/xyzzyaaap7
xyzzyaaas7=1.d0/xyzzyaaaq7
xyzzyaaat7=xyzzyaaar7*xyzzyaaar7
xyzzyaaau7=xyzzyaaas7*xyzzyaaas7
xyzzyaaar7=xyzzyaaal7*xyzzyaaar7*xyzzyaaat7
xyzzyaaas7=xyzzyaaam7*xyzzyaaas7*xyzzyaaau7
pote=(log((xyzzyaaap7+xyzzyaaal7)*(xyzzyaaaq7+xyzzyaaam7))*xyzzyaabf1-&
&(xyzzyaaar7+xyzzyaaas7)*xyzzyaabh1+((9.d0-xyzzyaaat7*xyzzyaaan7*15.d0&
&)*xyzzyaaar7*xyzzyaaat7+(9.d0-xyzzyaaau7*xyzzyaaao7*15.d0)*xyzzyaaas7&
&*xyzzyaaau7)*xyzzyaaap1+xyzzyaaaj7)+pote
enddo
endif
end subroutine xyzzyaadi1
real(dp) function xyzzyaadj1(b)
use slaarnaag,only : root_pi
use slaarnabt,only : erfc,exp_int
implicit none
real(dp),intent(in) :: b
integer xyzzyaaaa8
real(dp) xyzzyaaab8,xyzzyaaac8,xyzzyaaad8,xyzzyaaae8,xyzzyaaaf8,xyzzya&
&aag8,xyzzyaaah8,xyzzyaaai8,xyzzyaaaj8,xyzzyaaak8
real(dp),parameter :: xyzzyaaal8=1.d-7
xyzzyaaae8=2.d0*wigner_seitz_radius
xyzzyaaab8=2.d0*root_pi/b
xyzzyaaaj8=0.5d0/b
xyzzyaaak8=b*b
xyzzyaaaa8=1
xyzzyaaac8=0.d0
do
xyzzyaaad8=abs(dble(xyzzyaaaa8)*xyzzyaaae8)*xyzzyaaaj8
xyzzyaaad8=xyzzyaadk1(xyzzyaaad8)
xyzzyaaac8=xyzzyaaac8+xyzzyaaad8
if(abs(xyzzyaaab8*xyzzyaaad8)<xyzzyaaal8)exit
xyzzyaaaa8=xyzzyaaaa8+1
enddo
xyzzyaaaf8=xyzzyaaab8*xyzzyaaac8
xyzzyaaac8=-2.d0/(b*root_pi)
xyzzyaaad8=0.d0
xyzzyaaaa8=1
do
xyzzyaaad8=-(2.d0/abs(-dble(xyzzyaaaa8)*xyzzyaaae8))*(1.d0-erfc((abs(-&
&dble(xyzzyaaaa8)*xyzzyaaae8)*xyzzyaaaj8)))
xyzzyaaac8=xyzzyaaac8+2.d0*xyzzyaaad8
if(abs(2.d0*xyzzyaaad8)<xyzzyaaal8)exit
xyzzyaaaa8=xyzzyaaaa8+1
enddo
xyzzyaaaf8=xyzzyaaaf8+xyzzyaaac8
xyzzyaaag8=0.d0
xyzzyaaaa8=1
xyzzyaaah8=twopi/xyzzyaaae8
xyzzyaaab8=4.d0/xyzzyaaae8
do
xyzzyaaai8=xyzzyaaah8*xyzzyaaaa8
xyzzyaaag8=xyzzyaaag8+exp_int(xyzzyaaak8*xyzzyaaai8*xyzzyaaai8)
if(xyzzyaaai8*b>12.d0)exit
xyzzyaaaa8=xyzzyaaaa8+1
enddo
xyzzyaaag8=xyzzyaaab8*xyzzyaaag8
xyzzyaadj1=0.5d0*(xyzzyaaaf8+xyzzyaaag8)
end function xyzzyaadj1
real(dp) function xyzzyaadk1(x)
use slaarnaag,only : one_over_root_pi
use slaarnabt,only : gamma_ser
implicit none
real(dp),intent(in) :: x
integer xyzzyaaaa9
integer,parameter :: xyzzyaaab9=100
real(dp) xyzzyaaac9,xyzzyaaad9,xyzzyaaae9,xyzzyaaaf9,xyzzyaaag9,xyzzya&
&aah9,xyzzyaaai9
real(dp),parameter :: xyzzyaaaj9=1.d-12,xyzzyaaak9=1.d-30
if(x<0.d0)call errstop('VSR_TERM_1DHEG','Bad argments in short-ranged &
&part of quasi-1d potential.')
xyzzyaaac9=x*x
if(xyzzyaaac9<1.5d0)then
xyzzyaadk1=exp(xyzzyaaac9)*(1.d0-gamma_ser(xyzzyaaac9))
else
xyzzyaaae9=xyzzyaaac9+0.5d0
xyzzyaaaf9=1.d0/xyzzyaaak9
xyzzyaaag9=1.d0/xyzzyaaae9
xyzzyaaai9=xyzzyaaag9
do xyzzyaaaa9=1,xyzzyaaab9
xyzzyaaad9=dble(xyzzyaaaa9)*(0.5d0-dble(xyzzyaaaa9))
xyzzyaaae9=xyzzyaaae9+2
xyzzyaaag9=xyzzyaaad9*xyzzyaaag9+xyzzyaaae9
if(abs(xyzzyaaag9)<xyzzyaaak9)xyzzyaaag9=xyzzyaaak9
xyzzyaaaf9=xyzzyaaae9+xyzzyaaad9/xyzzyaaaf9
if(abs(xyzzyaaaf9)<xyzzyaaak9)xyzzyaaaf9=xyzzyaaak9
xyzzyaaag9=1.d0/xyzzyaaag9
xyzzyaaah9=xyzzyaaag9*xyzzyaaaf9
xyzzyaaai9=xyzzyaaai9*xyzzyaaah9
if(abs(xyzzyaaah9-1.d0)<xyzzyaaaj9)exit
enddo
if(xyzzyaaaa9==xyzzyaaab9+1)call errstop('VSR_TERM_1DHEG','Short-range&
&d part of quasi-1d interaction failed to converge.')
xyzzyaadk1=xyzzyaaai9*x*one_over_root_pi
endif
end function xyzzyaadk1
subroutine ewald_2d(n,vecji,pote,field)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: vecji(4,n)
real(dp),intent(out) :: pote
real(dp),intent(out),optional :: field(3,n)
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10,xyzzyaaad10,xyzzyaaae10
real(dp) xyzzyaaaf10,xyzzyaaag10,xyzzyaaah10,xyzzyaaai10,xyzzyaaaj10(2&
&),xyzzyaaak10,xyzzyaaal10,xyzzyaaam10,xyzzyaaan10,xyzzyaaao10,xyzzyaa&
&ap10,xyzzyaaaq10,xyzzyaaar10,xyzzyaaas10,xyzzyaaat10,xyzzyaaau10,xyzz&
&yaaav10,xyzzyaaaw10(2),xyzzyaaax10,xyzzyaaay10,xyzzyaaaz10,xyzzyaaba1&
&0,xyzzyaabb10,xyzzyaabc10,xyzzyaabd10,xyzzyaabe10,xyzzyaabf10,xyzzyaa&
&bg10
complex(dp) xyzzyaabh10
pote=0.d0
do xyzzyaaad10=1,n
xyzzyaaaj10=vecji(1:2,xyzzyaaad10)
xyzzyaaap10=vecji(3,xyzzyaaad10)
xyzzyaaag10=xyzzyaaaj10(1)*xyzzyaaaj10(1)+xyzzyaaaj10(2)*xyzzyaaaj10(2&
&)+xyzzyaaap10*xyzzyaaap10
xyzzyaaag10=sqrt(xyzzyaaag10*xyzzyaaaw1)+xyzzyaaag10+xyzzyaaav1
xyzzyaaah10=0.d0
xyzzyaaax10=0.d0
xyzzyaaay10=0.d0
xyzzyaaaz10=0.d0
do xyzzyaaae10=1,lsize
if(s_modsq(xyzzyaaae10)>=xyzzyaaag10)exit
do xyzzyaaac10=first_s_in_star(xyzzyaaae10),first_s_in_star(xyzzyaaae1&
&0+1)-1
xyzzyaaan10=s_lattice(1,xyzzyaaac10)+xyzzyaaaj10(1)
xyzzyaaao10=s_lattice(2,xyzzyaaac10)+xyzzyaaaj10(2)
xyzzyaaar10=xyzzyaaan10*xyzzyaaan10+xyzzyaaao10*xyzzyaaao10+xyzzyaaap1&
&0*xyzzyaaap10
if(xyzzyaaar10>=xyzzyaaav1)cycle
xyzzyaaak10=xyzzyaaar10*xyzzyaaax1
if(xyzzyaaak10>=16.d0)then
xyzzyaaam10=exp(-xyzzyaaak10)
xyzzyaaaq10=((xyzzyaaak10+xyzzyaabb1)*xyzzyaaam10/((xyzzyaaak10+xyzzya&
&abc1)*(xyzzyaaak10+xyzzyaaak10)+xyzzyaaba1))
xyzzyaaah10=xyzzyaaah10-xyzzyaaaq10
if(present(field).and.xyzzyaaar10>0.d0)then
xyzzyaabe10=1.d0/xyzzyaaar10*(xyzzyaaam10-xyzzyaaaq10)
xyzzyaaax10=xyzzyaaax10-xyzzyaaan10*xyzzyaabe10
xyzzyaaay10=xyzzyaaay10-xyzzyaaao10*xyzzyaabe10
xyzzyaaaz10=xyzzyaaaz10-xyzzyaaap10*xyzzyaabe10
endif
else
xyzzyaaaa10=int(xyzzyaaaz1*xyzzyaaak10+0.5d0)
xyzzyaaal10=xyzzyaaay1*xyzzyaaaa10-xyzzyaaak10
xyzzyaaaq10=((xyzzyaabm1(4*xyzzyaaaa10)*xyzzyaaal10+xyzzyaabm1(4*xyzzy&
&aaaa10+1))*xyzzyaaal10+xyzzyaabm1(4*xyzzyaaaa10+2))*xyzzyaaal10+xyzzy&
&aabm1(4*xyzzyaaaa10+3)
if(xyzzyaaak10>=(1.d-16))xyzzyaaaq10=xyzzyaaaq10-sqrt(pi_over_four/xyz&
&zyaaak10)
xyzzyaaah10=xyzzyaaah10+xyzzyaaaq10
if(present(field).and.xyzzyaaar10>0.d0)then
xyzzyaaam10=exp(-xyzzyaaak10)
xyzzyaabe10=1.d0/xyzzyaaar10*(xyzzyaaam10-xyzzyaaaq10)
xyzzyaaax10=xyzzyaaax10-xyzzyaaan10*xyzzyaabe10
xyzzyaaay10=xyzzyaaay10-xyzzyaaao10*xyzzyaabe10
xyzzyaaaz10=xyzzyaaaz10-xyzzyaaap10*xyzzyaabe10
endif
endif
enddo
enddo
xyzzyaaah10=xyzzyaaah10*xyzzyaabi1
xyzzyaaai10=0.d0
xyzzyaaba10=0.d0
xyzzyaabb10=0.d0
xyzzyaabc10=0.d0
xyzzyaaaw10(1)=b1(1)*xyzzyaaaj10(1)+b1(2)*xyzzyaaaj10(2)
xyzzyaaaw10(2)=b2(1)*xyzzyaaaj10(1)+b2(2)*xyzzyaaaj10(2)
xyzzyaaas10=cos(xyzzyaaaw10(1))
xyzzyaaat10=cos(xyzzyaaaw10(2))
xyzzyaaau10=sin(xyzzyaaaw10(1))
xyzzyaaav10=sin(xyzzyaaaw10(2))
xyzzyaace1(1)=cmplx(xyzzyaaas10,xyzzyaaau10,dp)
xyzzyaacf1(1)=cmplx(xyzzyaaat10,xyzzyaaav10,dp)
xyzzyaace1(-1)=cmplx(xyzzyaaas10,-xyzzyaaau10,dp)
xyzzyaacf1(-1)=cmplx(xyzzyaaat10,-xyzzyaaav10,dp)
do xyzzyaaab10=2,xyzzyaaaf1
xyzzyaace1(xyzzyaaab10)=xyzzyaace1(1)*xyzzyaace1(xyzzyaaab10-1)
enddo
do xyzzyaaab10=-2,xyzzyaaac1,-1
xyzzyaace1(xyzzyaaab10)=xyzzyaace1(-1)*xyzzyaace1(xyzzyaaab10+1)
enddo
do xyzzyaaab10=2,xyzzyaaag1
xyzzyaacf1(xyzzyaaab10)=xyzzyaacf1(1)*xyzzyaacf1(xyzzyaaab10-1)
enddo
do xyzzyaaab10=-2,xyzzyaaad1,-1
xyzzyaacf1(xyzzyaaab10)=xyzzyaacf1(-1)*xyzzyaacf1(xyzzyaaab10+1)
enddo
xyzzyaaak10=xyzzyaaap10*xyzzyaabd1
xyzzyaaat1=xyzzyaaak10*xyzzyaaak10
if(.not.present(field))then
if(xyzzyaaap10/=0.d0)then
xyzzyaaas1=abs(xyzzyaaak10)
do xyzzyaaac10=2,xyzzyaaab1
call xyzzyaaea1(xyzzyaabx1(xyzzyaaac10),xyzzyaaaf10)
xyzzyaaai10=real(xyzzyaace1(xyzzyaaao1(1,xyzzyaaac10))*xyzzyaacf1(xyzz&
&yaaao1(2,xyzzyaaac10)),dp)*xyzzyaabw1(xyzzyaaac10)*xyzzyaaaf10+xyzzya&
&aai10
enddo
else
do xyzzyaaac10=2,xyzzyaaab1
xyzzyaaak10=xyzzyaabz1(xyzzyaaac10)
if(xyzzyaaak10>=16.d0)then
xyzzyaaam10=exp(-xyzzyaaak10)
xyzzyaaaf10=((xyzzyaaak10+xyzzyaabb1)*xyzzyaaam10/((xyzzyaaak10+xyzzya&
&abc1)*(xyzzyaaak10+xyzzyaaak10)+xyzzyaaba1))
else
xyzzyaaaa10=int(xyzzyaaaz1*xyzzyaaak10+0.5d0)
xyzzyaaal10=xyzzyaaay1*xyzzyaaaa10-xyzzyaaak10
xyzzyaaaf10=((xyzzyaabm1(4*xyzzyaaaa10)*xyzzyaaal10+xyzzyaabm1(4*xyzzy&
&aaaa10+1))*xyzzyaaal10+xyzzyaabm1(4*xyzzyaaaa10+2))*xyzzyaaal10+xyzzy&
&aabm1(4*xyzzyaaaa10+3)
if(xyzzyaaak10>=(1.d-16))xyzzyaaaf10=xyzzyaaaf10-sqrt(pi_over_four/xyz&
&zyaaak10)
endif
xyzzyaaai10=real(xyzzyaace1(xyzzyaaao1(1,xyzzyaaac10))*xyzzyaacf1(xyzz&
&yaaao1(2,xyzzyaaac10)),dp)*xyzzyaaaf10*xyzzyaabl1+xyzzyaaai10
enddo
endif
call xyzzyaadz1(xyzzyaaaf10,xyzzyaabd10)
xyzzyaaai10=xyzzyaaaf10+xyzzyaaai10
pote=pote+xyzzyaaai10-xyzzyaaah10
else
if(xyzzyaaap10/=0.d0)then
xyzzyaaas1=abs(xyzzyaaak10)
do xyzzyaaac10=2,xyzzyaaab1
call xyzzyaaea1(xyzzyaabx1(xyzzyaaac10),xyzzyaaaf10)
xyzzyaabh10=xyzzyaace1(xyzzyaaao1(1,xyzzyaaac10))*xyzzyaacf1(xyzzyaaao&
&1(2,xyzzyaaac10))
xyzzyaaaf10=xyzzyaabw1(xyzzyaaac10)*xyzzyaaaf10
xyzzyaabf10=real(xyzzyaabh10,dp)*xyzzyaaaf10
xyzzyaabg10=aimag(xyzzyaabh10)*xyzzyaaaf10
xyzzyaaai10=xyzzyaabf10+xyzzyaaai10
xyzzyaabc10=xyzzyaabc10-xyzzyaabf10*xyzzyaaca1(xyzzyaaac10)
xyzzyaaba10=xyzzyaaba10-xyzzyaabg10*xyzzyaabs1(xyzzyaaac10)
xyzzyaabb10=xyzzyaabb10-xyzzyaabg10*xyzzyaabt1(xyzzyaaac10)
enddo
else
do xyzzyaaac10=2,xyzzyaaab1
xyzzyaaak10=xyzzyaabz1(xyzzyaaac10)
if(xyzzyaaak10>=16.d0)then
xyzzyaaam10=exp(-xyzzyaaak10)
xyzzyaaaf10=((xyzzyaaak10+xyzzyaabb1)*xyzzyaaam10/((xyzzyaaak10+xyzzya&
&abc1)*(xyzzyaaak10+xyzzyaaak10)+xyzzyaaba1))
else
xyzzyaaaa10=int(xyzzyaaaz1*xyzzyaaak10+0.5d0)
xyzzyaaal10=xyzzyaaay1*xyzzyaaaa10-xyzzyaaak10
xyzzyaaaf10=((xyzzyaabm1(4*xyzzyaaaa10)*xyzzyaaal10+xyzzyaabm1(4*xyzzy&
&aaaa10+1))*xyzzyaaal10+xyzzyaabm1(4*xyzzyaaaa10+2))*xyzzyaaal10+xyzzy&
&aabm1(4*xyzzyaaaa10+3)
if(xyzzyaaak10>=(1.d-16))xyzzyaaaf10=xyzzyaaaf10-sqrt(pi_over_four/xyz&
&zyaaak10)
endif
xyzzyaabh10=xyzzyaace1(xyzzyaaao1(1,xyzzyaaac10))*xyzzyaacf1(xyzzyaaao&
&1(2,xyzzyaaac10))
xyzzyaaaf10=xyzzyaaaf10*xyzzyaabl1
xyzzyaabf10=real(xyzzyaabh10,dp)*xyzzyaaaf10
xyzzyaabg10=aimag(xyzzyaabh10)*xyzzyaaaf10
xyzzyaaai10=xyzzyaabf10+xyzzyaaai10
xyzzyaabc10=xyzzyaabc10-xyzzyaabf10
xyzzyaaba10=xyzzyaaba10-xyzzyaabg10*xyzzyaabs1(xyzzyaaac10)
xyzzyaabb10=xyzzyaabb10-xyzzyaabg10*xyzzyaabt1(xyzzyaaac10)
enddo
endif
call xyzzyaadz1(xyzzyaaaf10,xyzzyaabd10)
xyzzyaaai10=xyzzyaaaf10+xyzzyaaai10
if(xyzzyaaap10/=0.d0)xyzzyaabc10=xyzzyaabc10+xyzzyaabd10/xyzzyaaap10
pote=pote+xyzzyaaai10-xyzzyaaah10
field(1,xyzzyaaad10)=xyzzyaaax10*xyzzyaabi1+xyzzyaaba10
field(2,xyzzyaaad10)=xyzzyaaay10*xyzzyaabi1+xyzzyaabb10
field(3,xyzzyaaad10)=xyzzyaaaz10*xyzzyaabi1+xyzzyaabc10
endif
enddo
end subroutine ewald_2d
subroutine ewald_3d(n,vecji,pote,pseudo,field)
!$ use openmp_base, only : newald_min
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: vecji(4,n)
real(dp),intent(out) :: pote
logical,intent(in) :: pseudo
real(dp),intent(out),optional :: field(3,n)
integer xyzzyaaaa11,xyzzyaaab11,xyzzyaaac11,xyzzyaaad11,xyzzyaaae11,xy&
&zzyaaaf11
real(dp) xyzzyaaag11,xyzzyaaah11,xyzzyaaai11,xyzzyaaaj11,xyzzyaaak11(3&
&),xyzzyaaal11(3),xyzzyaaam11,xyzzyaaan11,xyzzyaaao11,xyzzyaaap11,xyzz&
&yaaaq11,xyzzyaaar11,xyzzyaaas11,xyzzyaaat11,xyzzyaaau11,xyzzyaaav11,x&
&yzzyaaaw11,xyzzyaaax11,xyzzyaaay11,xyzzyaaaz11,xyzzyaaba11,xyzzyaabb1&
&1,xyzzyaabc11,xyzzyaabd11,xyzzyaabe11,xyzzyaabf11,xyzzyaabg11,xyzzyaa&
&bh11,xyzzyaabi11,xyzzyaabj11
complex(dp) xyzzyaabk11,xyzzyaabl11
logical xyzzyaabm11
pote=0.d0
xyzzyaabm11=.false.
if(present(field))xyzzyaabm11=.true.
!$omp parallel do default(none) shared(xyzzyaabm11,n,vecji,xyzzyaaay1,&
!$omp &field,xyzzyaaaw1,xyzzyaaav1,first_s_in_star,s_lattice,xyzzyaaax&
!$omp &1,xyzzyaaba1,xyzzyaabb1,xyzzyaabc1,xyzzyaaaz1,xyzzyaabm1,s_mods&
!$omp &q,calc_field,xyzzyaabi1,b1,b2,b3,xyzzyaaag1,xyzzyaaah1,xyzzyaaa&
!$omp &c1,xyzzyaaad1,xyzzyaaae1,xyzzyaaan1,xyzzyaaam1,xyzzyaabv1,xyzzy&
!$omp &aaau1,xyzzyaaao1,xyzzyaaby1,xyzzyaabs1,xyzzyaabt1,xyzzyaabu1,xy&
!$omp &zzyaaab1,pseudo,man_int_poly_cutoff,man_int_poly_order,man_int_&
!$omp &poly_coeffs) private(xyzzyaaae11,xyzzyaaal11,xyzzyaaag11,xyzzya&
!$omp &aah11,xyzzyaaaf11,xyzzyaaai11,xyzzyaaav11,xyzzyaaaw11,xyzzyaaax&
!$omp &11,xyzzyaaad11,xyzzyaabb11,xyzzyaabc11,xyzzyaabd11,xyzzyaabg11,&
!$omp &xyzzyaaam11,xyzzyaabf11,xyzzyaabe11,xyzzyaaau11,xyzzyaaaa11,xyz&
!$omp &zyaaan11,xyzzyaaaj11,xyzzyaaay11,xyzzyaaaz11,xyzzyaaba11,xyzzya&
!$omp &aak11,xyzzyaaao11,xyzzyaaap11,xyzzyaaaq11,xyzzyaaar11,xyzzyaaas&
!$omp &11,xyzzyaaat11,xyzzyaaab11,xyzzyaaac11,xyzzyaabl11,xyzzyaabk11,&
!$omp &xyzzyaabh11,xyzzyaabi11,xyzzyaabj11) reduction(+:pote) if(n>=ne&
!$omp &wald_min)
do xyzzyaaae11=1,n
if(pseudo)then
if(vecji(4,xyzzyaaae11)<=man_int_poly_cutoff)then
xyzzyaabh11=vecji(4,xyzzyaaae11)/man_int_poly_cutoff
xyzzyaabi11=0.d0
pote=pote-1.d0/vecji(4,xyzzyaaae11)+1.d0/man_int_poly_cutoff+1.d0/man_&
&int_poly_cutoff*(1.d0-xyzzyaabh11)*xyzzyaabh11*xyzzyaabh11
xyzzyaabi11=(0.5d0+xyzzyaabh11)*man_int_poly_coeffs(1)
xyzzyaabj11=xyzzyaabh11
do xyzzyaaaa11=2,man_int_poly_order
xyzzyaabj11=xyzzyaabj11*xyzzyaabh11
xyzzyaabi11=xyzzyaabi11+xyzzyaabj11*man_int_poly_coeffs(xyzzyaaaa11)
enddo
pote=pote+(1.d0-xyzzyaabh11)*(1.d0-xyzzyaabh11)*xyzzyaabi11
endif
endif
xyzzyaaal11=vecji(1:3,xyzzyaaae11)
xyzzyaaag11=vecji(4,xyzzyaaae11)**2
xyzzyaaah11=sqrt(xyzzyaaag11*xyzzyaaaw1)+xyzzyaaag11+xyzzyaaav1
xyzzyaaaf11=1
xyzzyaaai11=0.d0
xyzzyaaav11=0.d0
xyzzyaaaw11=0.d0
xyzzyaaax11=0.d0
do xyzzyaaaf11=1,lsize
if(s_modsq(xyzzyaaaf11)>=xyzzyaaah11)exit
do xyzzyaaad11=first_s_in_star(xyzzyaaaf11),first_s_in_star(xyzzyaaaf1&
&1+1)-1
xyzzyaabb11=s_lattice(1,xyzzyaaad11)+xyzzyaaal11(1)
xyzzyaabc11=s_lattice(2,xyzzyaaad11)+xyzzyaaal11(2)
xyzzyaabd11=s_lattice(3,xyzzyaaad11)+xyzzyaaal11(3)
xyzzyaabg11=xyzzyaabb11*xyzzyaabb11+xyzzyaabc11*xyzzyaabc11+xyzzyaabd1&
&1*xyzzyaabd11
if(xyzzyaabg11>=xyzzyaaav1)cycle
xyzzyaaam11=xyzzyaabg11*xyzzyaaax1
if(xyzzyaaam11>=16.d0)then
xyzzyaabf11=exp(-xyzzyaaam11)
xyzzyaabe11=((xyzzyaaam11+xyzzyaabb1)*xyzzyaabf11/((xyzzyaaam11+xyzzya&
&abc1)*(xyzzyaaam11+xyzzyaaam11)+xyzzyaaba1))
xyzzyaaai11=xyzzyaaai11-xyzzyaabe11
if(xyzzyaabm11.and.xyzzyaabg11>0)then
xyzzyaaau11=(1.d0/xyzzyaabg11)*(xyzzyaabf11-xyzzyaabe11)
xyzzyaaav11=xyzzyaaav11-xyzzyaabb11*xyzzyaaau11
xyzzyaaaw11=xyzzyaaaw11-xyzzyaabc11*xyzzyaaau11
xyzzyaaax11=xyzzyaaax11-xyzzyaabd11*xyzzyaaau11
endif
else
xyzzyaaaa11=int(xyzzyaaaz1*xyzzyaaam11+0.5d0)
xyzzyaaan11=xyzzyaaay1*xyzzyaaaa11-xyzzyaaam11
xyzzyaabe11=((xyzzyaabm1(4*xyzzyaaaa11)*xyzzyaaan11+xyzzyaabm1(4*xyzzy&
&aaaa11+1))*xyzzyaaan11+xyzzyaabm1(4*xyzzyaaaa11+2))*xyzzyaaan11+xyzzy&
&aabm1(4*xyzzyaaaa11+3)
if(xyzzyaaam11>=(1.d-16))xyzzyaabe11=xyzzyaabe11-sqrt(pi_over_four/xyz&
&zyaaam11)
xyzzyaaai11=xyzzyaaai11+xyzzyaabe11
if(xyzzyaabm11.and.xyzzyaabg11>0)then
xyzzyaabf11=exp(-xyzzyaaam11)
xyzzyaaau11=(1.d0/xyzzyaabg11)*(xyzzyaabf11-xyzzyaabe11)
xyzzyaaav11=xyzzyaaav11-xyzzyaabb11*xyzzyaaau11
xyzzyaaaw11=xyzzyaaaw11-xyzzyaabc11*xyzzyaaau11
xyzzyaaax11=xyzzyaaax11-xyzzyaabd11*xyzzyaaau11
endif
endif
enddo
enddo
xyzzyaaai11=xyzzyaaai11*xyzzyaabi1
xyzzyaaaj11=0.d0
xyzzyaaay11=0.d0
xyzzyaaaz11=0.d0
xyzzyaaba11=0.d0
xyzzyaaak11(1)=b1(1)*xyzzyaaal11(1)+b1(2)*xyzzyaaal11(2)+b1(3)*xyzzyaa&
&al11(3)
xyzzyaaak11(2)=b2(1)*xyzzyaaal11(1)+b2(2)*xyzzyaaal11(2)+b2(3)*xyzzyaa&
&al11(3)
xyzzyaaak11(3)=b3(1)*xyzzyaaal11(1)+b3(2)*xyzzyaaal11(2)+b3(3)*xyzzyaa&
&al11(3)
xyzzyaaao11=cos(xyzzyaaak11(1))
xyzzyaaap11=cos(xyzzyaaak11(2))
xyzzyaaaq11=cos(xyzzyaaak11(3))
xyzzyaaar11=sin(xyzzyaaak11(1))
xyzzyaaas11=sin(xyzzyaaak11(2))
xyzzyaaat11=sin(xyzzyaaak11(3))
xyzzyaacf1(1)=cmplx(xyzzyaaap11,xyzzyaaas11,dp)
xyzzyaacg1(1)=cmplx(xyzzyaaaq11,xyzzyaaat11,dp)
xyzzyaacf1(-1)=cmplx(xyzzyaaap11,-xyzzyaaas11,dp)
xyzzyaacg1(-1)=cmplx(xyzzyaaaq11,-xyzzyaaat11,dp)
do xyzzyaaab11=2,xyzzyaaag1
xyzzyaacf1(xyzzyaaab11)=xyzzyaacf1(1)*xyzzyaacf1(xyzzyaaab11-1)
enddo
do xyzzyaaac11=2,xyzzyaaah1
xyzzyaacg1(xyzzyaaac11)=xyzzyaacg1(1)*xyzzyaacg1(xyzzyaaac11-1)
enddo
do xyzzyaaab11=-2,xyzzyaaad1,-1
xyzzyaacf1(xyzzyaaab11)=xyzzyaacf1(-1)*xyzzyaacf1(xyzzyaaab11+1)
enddo
do xyzzyaaac11=-2,xyzzyaaae1,-1
xyzzyaacg1(xyzzyaaac11)=xyzzyaacg1(-1)*xyzzyaacg1(xyzzyaaac11+1)
enddo
if(.not.xyzzyaabm11)then
do xyzzyaaab11=1,xyzzyaaan1(0)
xyzzyaaaj11=real(xyzzyaacf1(xyzzyaaam1(2*xyzzyaaab11-1))*xyzzyaacg1(xy&
&zzyaaam1(2*xyzzyaaab11)),dp)*xyzzyaabv1(xyzzyaaab11)+xyzzyaaaj11
enddo
xyzzyaabl11=cmplx(xyzzyaaao11,-xyzzyaaar11,dp)
xyzzyaabk11=xyzzyaabl11
xyzzyaaac11=xyzzyaaan1(0)+1
do xyzzyaaaa11=-1,xyzzyaaac1,-1
do xyzzyaaab11=1,xyzzyaaan1(xyzzyaaaa11)
xyzzyaaaj11=real(xyzzyaabk11*xyzzyaacf1(xyzzyaaam1(2*xyzzyaaac11-1))*x&
&yzzyaacg1(xyzzyaaam1(2*xyzzyaaac11)),dp)*xyzzyaabv1(xyzzyaaac11)+xyzz&
&yaaaj11
xyzzyaaac11=xyzzyaaac11+1
enddo
xyzzyaabk11=xyzzyaabk11*xyzzyaabl11
enddo
pote=pote+(xyzzyaaaj11-xyzzyaaai11-xyzzyaaau1)
else
xyzzyaace1(-1)=cmplx(xyzzyaaao11,-xyzzyaaar11,dp)
do xyzzyaaaa11=-2,xyzzyaaac1,-1
xyzzyaace1(xyzzyaaaa11)=xyzzyaace1(-1)*xyzzyaace1(xyzzyaaaa11+1)
enddo
do xyzzyaaaa11=2,xyzzyaaab1
xyzzyaabk11=xyzzyaace1(xyzzyaaao1(1,xyzzyaaaa11))*xyzzyaacf1(xyzzyaaao&
&1(2,xyzzyaaaa11))*xyzzyaacg1(xyzzyaaao1(3,xyzzyaaaa11))
xyzzyaaaj11=xyzzyaaaj11+real(xyzzyaabk11,dp)*xyzzyaaby1(xyzzyaaaa11)
xyzzyaaau11=aimag(xyzzyaabk11)*xyzzyaaby1(xyzzyaaaa11)
xyzzyaaay11=xyzzyaaay11-xyzzyaaau11*xyzzyaabs1(xyzzyaaaa11)
xyzzyaaaz11=xyzzyaaaz11-xyzzyaaau11*xyzzyaabt1(xyzzyaaaa11)
xyzzyaaba11=xyzzyaaba11-xyzzyaaau11*xyzzyaabu1(xyzzyaaaa11)
enddo
pote=pote+(xyzzyaaaj11-xyzzyaaai11-xyzzyaaau1)
field(1,xyzzyaaae11)=xyzzyaaav11*xyzzyaabi1+xyzzyaaay11
field(2,xyzzyaaae11)=xyzzyaaaw11*xyzzyaabi1+xyzzyaaaz11
field(3,xyzzyaaae11)=xyzzyaaax11*xyzzyaabi1+xyzzyaaba11
endif
enddo
!$omp end parallel do
end subroutine ewald_3d
subroutine xyzzyaadl1(ie,vecji,well_width,well_height,direct,just_want&
&_j_gt_i)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot),well_width,well_height
real(dp),intent(out) :: direct
logical,intent(in) :: just_want_j_gt_i
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12,xy&
&zzyaaaf12
real(dp),allocatable,save :: xyzzyaaag12(:,:)
logical,save :: xyzzyaaah12=.true.
if(xyzzyaaah12)then
if(heg_nlayers>1)then
allocate(xyzzyaaag12(nspin,nspin),stat=xyzzyaaaf12)
call check_alloc(xyzzyaaaf12,'EE_EVAL_SQUARE_WELL','sqrt_wsq_zsq')
do xyzzyaaaa12=1,nspin
xyzzyaaag12(xyzzyaaaa12,xyzzyaaaa12)=well_width
do xyzzyaaae12=xyzzyaaaa12+1,nspin
xyzzyaaag12(xyzzyaaae12,xyzzyaaaa12)=well_width**2-(heg_zlayer(heg_lay&
&er(xyzzyaaae12))-heg_zlayer(heg_layer(xyzzyaaaa12)))**2-(heg_ylayer(h&
&eg_layer(xyzzyaaae12))-heg_ylayer(heg_layer(xyzzyaaaa12)))**2
if(xyzzyaaag12(xyzzyaaae12,xyzzyaaaa12)<0)call errstop('EE_EVAL_SQUARE&
&_WELL','well_width < layer separation.')
xyzzyaaag12(xyzzyaaae12,xyzzyaaaa12)=sqrt(xyzzyaaag12(xyzzyaaae12,xyzz&
&yaaaa12))
xyzzyaaag12(xyzzyaaaa12,xyzzyaaae12)=xyzzyaaag12(xyzzyaaae12,xyzzyaaaa&
&12)
enddo
enddo
endif
xyzzyaaah12=.false.
endif
direct=0.d0
xyzzyaaae12=which_spin(ie)
if(heg_nlayers==1)then
if(just_want_j_gt_i)then
xyzzyaaac12=nele(xyzzyaaae12)
do xyzzyaaaa12=xyzzyaaae12+1,nspin
xyzzyaaad12=xyzzyaaac12+1
xyzzyaaac12=xyzzyaaac12+nele(xyzzyaaaa12)
do xyzzyaaab12=xyzzyaaad12,xyzzyaaac12
if(vecji(4,xyzzyaaab12)<well_width)direct=direct+well_height
enddo
enddo
else
xyzzyaaac12=0
do xyzzyaaaa12=1,xyzzyaaae12-1
xyzzyaaad12=xyzzyaaac12+1
xyzzyaaac12=xyzzyaaac12+nele(xyzzyaaaa12)
do xyzzyaaab12=xyzzyaaad12,xyzzyaaac12
if(vecji(4,xyzzyaaab12)<well_width)direct=direct+well_height
enddo
enddo
xyzzyaaad12=xyzzyaaac12+1
xyzzyaaac12=xyzzyaaac12+nele(xyzzyaaae12)
do xyzzyaaaa12=xyzzyaaae12+1,nspin
xyzzyaaad12=xyzzyaaac12+1
xyzzyaaac12=xyzzyaaac12+nele(xyzzyaaaa12)
do xyzzyaaab12=xyzzyaaad12,xyzzyaaac12
if(vecji(4,xyzzyaaab12)<well_width)direct=direct+well_height
enddo
enddo
direct=0.5d0*direct
endif
else
if(just_want_j_gt_i)then
xyzzyaaac12=nele(xyzzyaaae12)
do xyzzyaaaa12=xyzzyaaae12+1,nspin
xyzzyaaad12=xyzzyaaac12+1
xyzzyaaac12=xyzzyaaac12+nele(xyzzyaaaa12)
do xyzzyaaab12=xyzzyaaad12,xyzzyaaac12
if(vecji(4,xyzzyaaab12)<xyzzyaaag12(xyzzyaaaa12,xyzzyaaae12))direct=di&
&rect+well_height
enddo
enddo
else
xyzzyaaac12=0
do xyzzyaaaa12=1,xyzzyaaae12-1
xyzzyaaad12=xyzzyaaac12+1
xyzzyaaac12=xyzzyaaac12+nele(xyzzyaaaa12)
do xyzzyaaab12=xyzzyaaad12,xyzzyaaac12
if(vecji(4,xyzzyaaab12)<xyzzyaaag12(xyzzyaaaa12,xyzzyaaae12))direct=di&
&rect+well_height
enddo
enddo
xyzzyaaad12=xyzzyaaac12+1
xyzzyaaac12=xyzzyaaac12+nele(xyzzyaaae12)
do xyzzyaaaa12=xyzzyaaae12+1,nspin
xyzzyaaad12=xyzzyaaac12+1
xyzzyaaac12=xyzzyaaac12+nele(xyzzyaaaa12)
do xyzzyaaab12=xyzzyaaad12,xyzzyaaac12
if(vecji(4,xyzzyaaab12)<xyzzyaaag12(xyzzyaaaa12,xyzzyaaae12))direct=di&
&rect+well_height
enddo
enddo
direct=0.5d0*direct
endif
endif
end subroutine xyzzyaadl1
subroutine xyzzyaadm1(ie,vecji,direct,just_want_j_gt_i)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: direct
logical,intent(in) :: just_want_j_gt_i
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13,xyzzyaaae13,xy&
&zzyaaaf13
real(dp) xyzzyaaag13
real(dp),allocatable,save :: xyzzyaaah13(:,:)
logical,save :: xyzzyaaai13=.true.
if(xyzzyaaai13)then
if(heg_nlayers>1)then
allocate(xyzzyaaah13(nspin,nspin),stat=xyzzyaaaf13)
call check_alloc(xyzzyaaaf13,'EE_EVAL_POSCHL_TELLER','yzsq')
do xyzzyaaae13=1,nspin
xyzzyaaah13(xyzzyaaae13,xyzzyaaae13)=0.d0
do xyzzyaaaa13=xyzzyaaae13+1,nspin
xyzzyaaah13(xyzzyaaaa13,xyzzyaaae13)=(heg_ylayer(heg_layer(xyzzyaaaa13&
&))-heg_ylayer(heg_layer(xyzzyaaae13)))**2+(heg_zlayer(heg_layer(xyzzy&
&aaaa13))-heg_zlayer(heg_layer(xyzzyaaae13)))**2
xyzzyaaah13(xyzzyaaae13,xyzzyaaaa13)=xyzzyaaah13(xyzzyaaaa13,xyzzyaaae&
&13)
enddo
enddo
endif
xyzzyaaai13=.false.
endif
direct=0.d0
xyzzyaaae13=which_spin(ie)
if(heg_nlayers==1)then
if(just_want_j_gt_i)then
xyzzyaaac13=nele(xyzzyaaae13)
do xyzzyaaaa13=xyzzyaaae13+1,nspin
xyzzyaaad13=xyzzyaaac13+1
xyzzyaaac13=xyzzyaaac13+nele(xyzzyaaaa13)
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
if(vecji(4,xyzzyaaab13)<wigner_seitz_radius)direct=direct+xyzzyaadn1(v&
&ecji(4,xyzzyaaab13))
enddo
enddo
else
xyzzyaaac13=0
do xyzzyaaaa13=1,xyzzyaaae13-1
xyzzyaaad13=xyzzyaaac13+1
xyzzyaaac13=xyzzyaaac13+nele(xyzzyaaaa13)
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
if(vecji(4,xyzzyaaab13)<wigner_seitz_radius)direct=direct+xyzzyaadn1(v&
&ecji(4,xyzzyaaab13))
enddo
enddo
xyzzyaaad13=xyzzyaaac13+1
xyzzyaaac13=xyzzyaaac13+nele(xyzzyaaae13)
do xyzzyaaaa13=xyzzyaaae13+1,nspin
xyzzyaaad13=xyzzyaaac13+1
xyzzyaaac13=xyzzyaaac13+nele(xyzzyaaaa13)
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
if(vecji(4,xyzzyaaab13)<wigner_seitz_radius)direct=direct+xyzzyaadn1(v&
&ecji(4,xyzzyaaab13))
enddo
enddo
direct=0.5d0*direct
endif
else
if(just_want_j_gt_i)then
xyzzyaaac13=nele(xyzzyaaae13)
do xyzzyaaaa13=xyzzyaaae13+1,nspin
xyzzyaaad13=xyzzyaaac13+1
xyzzyaaac13=xyzzyaaac13+nele(xyzzyaaaa13)
if(heg_layer(xyzzyaaae13)==heg_layer(xyzzyaaaa13))then
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
if(vecji(4,xyzzyaaab13)<wigner_seitz_radius)direct=direct+xyzzyaadn1(v&
&ecji(4,xyzzyaaab13))
enddo
else
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
xyzzyaaag13=sqrt(vecji(4,xyzzyaaab13)**2+xyzzyaaah13(xyzzyaaaa13,xyzzy&
&aaae13))
if(xyzzyaaag13<wigner_seitz_radius)direct=direct+xyzzyaadn1(xyzzyaaag1&
&3)
enddo
endif
enddo
else
xyzzyaaac13=0
do xyzzyaaaa13=1,xyzzyaaae13-1
xyzzyaaad13=xyzzyaaac13+1
xyzzyaaac13=xyzzyaaac13+nele(xyzzyaaaa13)
if(heg_layer(xyzzyaaae13)==heg_layer(xyzzyaaaa13))then
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
if(vecji(4,xyzzyaaab13)<wigner_seitz_radius)direct=direct+xyzzyaadn1(v&
&ecji(4,xyzzyaaab13))
enddo
else
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
xyzzyaaag13=sqrt(vecji(4,xyzzyaaab13)**2+xyzzyaaah13(xyzzyaaaa13,xyzzy&
&aaae13))
if(xyzzyaaag13<wigner_seitz_radius)direct=direct+xyzzyaadn1(xyzzyaaag1&
&3)
enddo
endif
enddo
xyzzyaaad13=xyzzyaaac13+1
xyzzyaaac13=xyzzyaaac13+nele(xyzzyaaae13)
do xyzzyaaaa13=xyzzyaaae13+1,nspin
xyzzyaaad13=xyzzyaaac13+1
xyzzyaaac13=xyzzyaaac13+nele(xyzzyaaaa13)
if(heg_layer(xyzzyaaae13)==heg_layer(xyzzyaaaa13))then
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
if(vecji(4,xyzzyaaab13)<wigner_seitz_radius)direct=direct+xyzzyaadn1(v&
&ecji(4,xyzzyaaab13))
enddo
else
do xyzzyaaab13=xyzzyaaad13,xyzzyaaac13
xyzzyaaag13=sqrt(vecji(4,xyzzyaaab13)**2+xyzzyaaah13(xyzzyaaaa13,xyzzy&
&aaae13))
if(xyzzyaaag13<wigner_seitz_radius)direct=direct+xyzzyaadn1(xyzzyaaag1&
&3)
enddo
endif
enddo
direct=0.5d0*direct
endif
endif
end subroutine xyzzyaadm1
real(dp) function xyzzyaadn1(r)
implicit none
real(dp),intent(in) :: r
real(dp) xyzzyaaaa14
xyzzyaaaa14=cosh(r*man_int_params(1))
xyzzyaadn1=xyzzyaadd1/(xyzzyaaaa14*xyzzyaaaa14)
end function xyzzyaadn1
subroutine xyzzyaado1(ie,vecji,int_pot,just_want_j_gt_i)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: int_pot
logical,intent(in) :: just_want_j_gt_i
logical xyzzyaaaa15
integer xyzzyaaab15,xyzzyaaac15,xyzzyaaad15
xyzzyaaaa15=periodicity==2
xyzzyaaac15=which_spin(ie)
do xyzzyaaab15=1,netot
if(ie/=xyzzyaaab15)then
if(xyzzyaaab15>ie.or..not.just_want_j_gt_i)then
xyzzyaaad15=which_spin(xyzzyaaab15)
if(.not.hard_op_spins.or.xyzzyaaac15/=xyzzyaaad15)then
if(vecji(4,xyzzyaaab15)<=hard_diam)call errstop('EE_EVAL_HARD_SPHERE',&
&'Two hard spheres are overlapping.')
endif
endif
endif
enddo
int_pot=0.d0
if(xyzzyaabn1/=0.d0)then
if(.not.just_want_j_gt_i)then
do xyzzyaaab15=1,ie-1
int_pot=int_pot+xyzzyaadp1(vecji(1:3,xyzzyaaab15),vecji(4,xyzzyaaab15)&
&,xyzzyaaaa15)
enddo
endif
do xyzzyaaab15=ie+1,netot
int_pot=int_pot+xyzzyaadp1(vecji(1:3,xyzzyaaab15),vecji(4,xyzzyaaab15)&
&,xyzzyaaaa15)
enddo
int_pot=int_pot*xyzzyaabn1
if(.not.just_want_j_gt_i)int_pot=int_pot*0.5d0
endif
int_pot=int_pot+self_term
end subroutine xyzzyaado1
real(dp) function xyzzyaadp1(rvec,r,want_periodic_images)
implicit none
real(dp),intent(in) :: rvec(3),r
logical,intent(in) :: want_periodic_images
integer xyzzyaaaa16,xyzzyaaab16
real(dp) xyzzyaaac16(dimensionality)
if(dimensionality==2)then
if(want_periodic_images)then
xyzzyaadp1=1.d0/r**3+xyzzyaabo1*(1.d0+xyzzyaabp1*r*r)
else
xyzzyaadp1=1.d0/r**3
endif
else
call errstop('PERIODIC_R_TO_MINUS3','Only the 2D case is coded up at p&
&resent.')
endif
if(want_periodic_images)then
do xyzzyaaaa16=2,xyzzyaaaj1
do xyzzyaaab16=first_s_in_star(xyzzyaaaa16),first_s_in_star(xyzzyaaaa1&
&6+1)-1
xyzzyaaac16(1:dimensionality)=rvec(1:dimensionality)+s_lattice(1:dimen&
&sionality,xyzzyaaab16)
xyzzyaadp1=xyzzyaadp1+dot_product(xyzzyaaac16(1:dimensionality),xyzzya&
&aac16(1:dimensionality))**(-1.5d0)
enddo
enddo
endif
end function xyzzyaadp1
real(dp) function xyzzyaadq1(rvec,r,theta,want_periodic_images)
implicit none
real(dp),intent(in) :: rvec(3),r,theta
logical,intent(in) :: want_periodic_images
integer xyzzyaaaa17,xyzzyaaab17
real(dp) xyzzyaaac17(dimensionality)
if(dimensionality==2)then
if(want_periodic_images)then
xyzzyaadq1=(1.d0-3.d0*sin(theta)*sin(theta)*(rvec(1)*rvec(1))/(r*r))/r&
&**3+xyzzyaabo1*(1.d0+xyzzyaabp1*r*r)
else
xyzzyaadq1=(1.d0-3.d0*sin(theta)*sin(theta)*(rvec(1)*rvec(1))/(r*r))/r&
&**3
endif
else
call errstop('PERIODIC_R_TO_MINUS3','Only the 2D case is coded up at p&
&resent.')
endif
if(want_periodic_images)then
do xyzzyaaaa17=2,xyzzyaaaj1
do xyzzyaaab17=first_s_in_star(xyzzyaaaa17),first_s_in_star(xyzzyaaaa1&
&7+1)-1
xyzzyaaac17(1:dimensionality)=rvec(1:dimensionality)+s_lattice(1:dimen&
&sionality,xyzzyaaab17)
xyzzyaadq1=xyzzyaadq1+(1.d0-3.d0*sin(theta)*sin(theta)*(xyzzyaaac17(1)&
&*xyzzyaaac17(1))/(dot_product(xyzzyaaac17(1:dimensionality),xyzzyaaac&
&17(1:dimensionality))))*dot_product(xyzzyaaac17(1:dimensionality),xyz&
&zyaaac17(1:dimensionality))**(-1.5d0)
enddo
enddo
endif
end function xyzzyaadq1
subroutine xyzzyaadr1(ie,vecji,int_pot,just_want_j_gt_i)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: int_pot
logical,intent(in) :: just_want_j_gt_i
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18
xyzzyaaab18=which_spin(ie)
int_pot=0.d0
xyzzyaaac18=nele(xyzzyaaab18)
do xyzzyaaaa18=xyzzyaaab18+1,nspin
xyzzyaaad18=xyzzyaaac18+1
xyzzyaaac18=xyzzyaaac18+nele(xyzzyaaaa18)
int_pot=int_pot+sum(polynomial_interaction(vecji(4,xyzzyaaad18:xyzzyaa&
&ac18)),vecji(4,xyzzyaaad18:xyzzyaaac18)<man_int_poly_cutoff)
enddo
if(.not.just_want_j_gt_i)then
xyzzyaaac18=0
do xyzzyaaaa18=1,xyzzyaaab18-1
xyzzyaaad18=xyzzyaaac18+1
xyzzyaaac18=xyzzyaaac18+nele(xyzzyaaab18)
int_pot=int_pot+sum(polynomial_interaction(vecji(4,xyzzyaaad18:xyzzyaa&
&ac18)),vecji(4,xyzzyaaad18:xyzzyaaac18)<man_int_poly_cutoff)
enddo
int_pot=0.5d0*int_pot
endif
end subroutine xyzzyaadr1
real(dp) elemental function polynomial_interaction(r) result(vij)
implicit none
real(dp),intent(in) :: r
integer icoeff
vij=man_int_poly_coeffs(man_int_poly_order)
do icoeff=man_int_poly_order-1,0,-1
vij=vij*r+man_int_poly_coeffs(icoeff)
enddo
end function polynomial_interaction
real(dp) elemental function tilted_polynomial_interaction(x,r) result(&
&vij)
implicit none
real(dp),intent(in) :: r,x
integer icoeff
vij=man_int_tilt_poly_coeffs(man_int_tilt_poly_order)
do icoeff=man_int_tilt_poly_order-1,0,-1
vij=vij*r+man_int_tilt_poly_coeffs(icoeff)
enddo
vij=vij*(2.d0*x*x/(r*r)-1.d0)
end function tilted_polynomial_interaction
subroutine xyzzyaads1(ie,vecji,v,just_want_j_gt_i,isinf)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: v
logical,intent(in) :: just_want_j_gt_i
logical,intent(out) :: isinf
integer xyzzyaaaa19,xyzzyaaab19
real(dp) xyzzyaaac19
v=0.d0
isinf=.false.
xyzzyaaab19=which_spin(ie)
if(just_want_j_gt_i)then
do xyzzyaaaa19=ie+1,netot
xyzzyaaac19=rstar+abs(vecji(3,xyzzyaaaa19))
if(vecji(4,xyzzyaaaa19)==0.d0)then
isinf=.true.
return
endif
v=v+xyzzyaacd1(which_spin(xyzzyaaaa19),xyzzyaaab19)*(log(2.d0*xyzzyaaa&
&c19/vecji(4,xyzzyaaaa19))-euler)/xyzzyaaac19
enddo
else
do xyzzyaaaa19=1,netot
if(xyzzyaaaa19==ie)cycle
xyzzyaaac19=rstar+abs(vecji(3,xyzzyaaaa19))
if(vecji(4,xyzzyaaaa19)==0.d0)then
isinf=.true.
return
endif
v=v+xyzzyaacd1(which_spin(xyzzyaaaa19),xyzzyaaab19)*(log(2.d0*xyzzyaaa&
&c19/vecji(4,xyzzyaaaa19))-euler)/xyzzyaaac19
enddo
v=0.5d0*v
endif
end subroutine xyzzyaads1
subroutine xyzzyaadt1(ie,vecji,v,just_want_j_gt_i,isinf)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: v
logical,intent(in) :: just_want_j_gt_i
logical,intent(out) :: isinf
integer xyzzyaaaa20,xyzzyaaab20
v=0.d0
isinf=.false.
xyzzyaaab20=which_spin(ie)
if(just_want_j_gt_i)then
do xyzzyaaaa20=ie+1,netot
if(vecji(4,xyzzyaaaa20)==0.d0)then
isinf=.true.
return
endif
v=v+xyzzyaacd1(which_spin(xyzzyaaaa20),xyzzyaaab20)*pot_2d_int(vecji(4&
&,xyzzyaaaa20),vecji(3,xyzzyaaaa20))
enddo
else
do xyzzyaaaa20=1,netot
if(xyzzyaaaa20==ie)cycle
if(vecji(4,xyzzyaaaa20)==0.d0)then
isinf=.true.
return
endif
v=v+xyzzyaacd1(which_spin(xyzzyaaaa20),xyzzyaaab20)*pot_2d_int(vecji(4&
&,xyzzyaaaa20),vecji(3,xyzzyaaaa20))
enddo
v=0.5d0*v
endif
end subroutine xyzzyaadt1
real(dp) function pot_2d_int(rp,z)
use slaarnaag,only : log2
implicit none
real(dp),intent(in) :: rp,z
integer xyzzyaaaa21
real(dp) xyzzyaaab21,xyzzyaaac21,xyzzyaaad21,xyzzyaaae21,xyzzyaaaf21,x&
&yzzyaaag21,xyzzyaaah21,xyzzyaaai21,xyzzyaaaj21,xyzzyaaak21,xyzzyaaal2&
&1,xyzzyaaam21,xyzzyaaan21
real(dp),parameter :: xyzzyaaao21(60)=(/0.5772156649015328606065120900&
&82d0,-0.422784335098467139393487909918d0,-0.9227843350984671393934879&
&09918d0,-1.25611766843180047272682124325d0,-1.50611766843180047272682&
&124325d0,-1.70611766843180047272682124325d0,-1.8727843350984671393934&
&8790992d0,-2.01564147795560999653634505277d0,-2.140641477955609996536&
&34505277d0,-2.25175258906672110764745616389d0,-2.35175258906672110764&
&745616389d0,-2.44266167997581201673836525479d0,-2.5259950133091453500&
&7169858813d0,-2.60291809023222227314862166505d0,-2.674346661660793701&
&72005023648d0,-2.74101332832746036838671690315d0,-2.80351332832746036&
&838671690315d0,-2.86233685773922507426906984432d0,-2.9178924132947806&
&2982462539988d0,-2.97052399224214905087725697883d0,-3.020523992242149&
&05087725697883d0,-3.06814303986119666992487602645d0,-3.11359758531574&
&212447033057190d0,-3.15707584618530734186163491973d0,-3.1987425128519&
&7400852830158639d0,-3.23874251285197400852830158639d0,-3.277204051313&
&51247006676312485d0,-3.31424108835054950710380016189d0,-3.34995537406&
&483522138951444761d0,-3.38443813268552487656192824071d0,-3.4177714660&
&1885820989526157404d0,-3.45002953053498724215332609017d0,-3.481279530&
&53498724215332609017d0,-3.51158256083801754518362912047d0,-3.54099432&
&554389989812480559106d0,-3.56956575411532846955337701963d0,-3.5973435&
&3189310624733115479741d0,-3.62437055892013327435818182444d0,-3.650686&
&34839381748488449761391d0,-3.67632737403484312591013863955d0,-3.70132&
&737403484312591013863955d0,-3.72571761793728215030038254199d0,-3.7495&
&2714174680595982419206580d0,-3.77278295570029433191721532162d0,-3.795&
&51022842756705918994259434d0,-3.81773245064978928141216481657d0,-3.83&
&947158108457189010781699048d0,-3.86074817682925274117164677771d0,-3.8&
&8158151016258607450498011105d0,-3.90198967342789219695395970288d0,-3.&
&92198967342789219695395970288d0,-3.94159751656514709891474401661d0,-3&
&.96082828579591632968397478584d0,-3.97969621032421821647642761603d0,-&
&3.99821472884273673499494613455d0,-4.01639654702455491681312795273d0,&
&-4.03425368988169777395598509558d0,-4.05179754953082058097352895523d0&
&,-4.06903892884116540855973585179d0,-4.08598808138353828991566805518d&
&0/)
real(dp),parameter :: xyzzyaaap21=euler-log2
xyzzyaaan21=1.d0/(rstar+abs(z))
xyzzyaaab21=rp*xyzzyaaan21
if(xyzzyaadf1)then
pot_2d_int=-(log(xyzzyaaab21/(1.d0+xyzzyaaab21))+xyzzyaaap21*exp(-xyzz&
&yaaab21))*xyzzyaaan21
else
if(xyzzyaaab21>18.d0)then
xyzzyaaac21=0.d0
xyzzyaaad21=100.d0
xyzzyaaaf21=1.d0
xyzzyaaaj21=1.d0/xyzzyaaab21
xyzzyaaal21=xyzzyaaaj21*xyzzyaaaj21
do xyzzyaaaa21=0,200
xyzzyaaae21=xyzzyaaad21
if(mod(xyzzyaaaa21,2)==0)then
xyzzyaaad21=xyzzyaaaf21*xyzzyaaaf21*xyzzyaaaj21
else
xyzzyaaad21=-xyzzyaaaf21*xyzzyaaaf21*xyzzyaaaj21
endif
if(abs(xyzzyaaad21)>abs(xyzzyaaae21).or.abs(xyzzyaaad21/xyzzyaaac21)<1&
&.d-18)exit
xyzzyaaac21=xyzzyaaac21+xyzzyaaad21
xyzzyaaaf21=xyzzyaaaf21*dble(2*xyzzyaaaa21+1)
xyzzyaaaj21=xyzzyaaaj21*xyzzyaaal21
enddo
pot_2d_int=xyzzyaaac21*xyzzyaaan21
else
xyzzyaaah21=0.d0
xyzzyaaaf21=1.d0
xyzzyaaag21=1.d0
xyzzyaaak21=xyzzyaaab21*xyzzyaaab21
xyzzyaaam21=log(0.5d0*xyzzyaaab21)
xyzzyaaaj21=1.d0
do xyzzyaaaa21=1,60
if(mod(xyzzyaaaa21,2)==0)then
xyzzyaaai21=((xyzzyaaam21+xyzzyaaao21(xyzzyaaaa21))/(xyzzyaaaf21*xyzzy&
&aaaf21)-xyzzyaaab21/(xyzzyaaag21*xyzzyaaag21))*xyzzyaaaj21
else
xyzzyaaai21=-((xyzzyaaam21+xyzzyaaao21(xyzzyaaaa21))/(xyzzyaaaf21*xyzz&
&yaaaf21)-xyzzyaaab21/(xyzzyaaag21*xyzzyaaag21))*xyzzyaaaj21
endif
xyzzyaaah21=xyzzyaaah21+xyzzyaaai21
if(abs(xyzzyaaai21/xyzzyaaah21)<1.d-18)exit
xyzzyaaaf21=xyzzyaaaf21*dble(2*xyzzyaaaa21)
xyzzyaaag21=xyzzyaaag21*dble(2*xyzzyaaaa21+1)
xyzzyaaaj21=xyzzyaaaj21*xyzzyaaak21
enddo
pot_2d_int=xyzzyaaah21*xyzzyaaan21
endif
endif
end function pot_2d_int
subroutine xyzzyaadu1(ie,vecji,int_pot,just_want_j_gt_i)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: int_pot
logical,intent(in) :: just_want_j_gt_i
logical xyzzyaaaa22
integer xyzzyaaab22
int_pot=0.d0
xyzzyaaaa22=periodicity==2
if(.not.just_want_j_gt_i)then
do xyzzyaaab22=1,ie-1
int_pot=int_pot+xyzzyaadp1(vecji(1:3,xyzzyaaab22),vecji(4,xyzzyaaab22)&
&,xyzzyaaaa22)
enddo
endif
do xyzzyaaab22=ie+1,netot
int_pot=int_pot+xyzzyaadp1(vecji(1:3,xyzzyaaab22),vecji(4,xyzzyaaab22)&
&,xyzzyaaaa22)
enddo
int_pot=int_pot*xyzzyaabq1
if(.not.just_want_j_gt_i)int_pot=int_pot*0.5d0
int_pot=int_pot+self_term
end subroutine xyzzyaadu1
subroutine xyzzyaadv1(ie,vecji,int_pot,just_want_j_gt_i)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: int_pot
logical,intent(in) :: just_want_j_gt_i
logical xyzzyaaaa23
integer xyzzyaaab23
int_pot=0.d0
xyzzyaaaa23=periodicity==2
if(.not.just_want_j_gt_i)then
do xyzzyaaab23=1,ie-1
int_pot=int_pot+xyzzyaabq1*xyzzyaadp1(vecji(1:3,xyzzyaaab23),vecji(4,x&
&yzzyaaab23),xyzzyaaaa23)
if(vecji(4,xyzzyaaab23)<man_int_poly_cutoff)then
int_pot=int_pot+polynomial_interaction(vecji(4,xyzzyaaab23))
int_pot=int_pot-xyzzyaabq1*xyzzyaadp1(vecji(1:3,xyzzyaaab23),vecji(4,x&
&yzzyaaab23),.false.)
endif
enddo
endif
do xyzzyaaab23=ie+1,netot
int_pot=int_pot+xyzzyaabq1*xyzzyaadp1(vecji(1:3,xyzzyaaab23),vecji(4,x&
&yzzyaaab23),xyzzyaaaa23)
if(vecji(4,xyzzyaaab23)<man_int_poly_cutoff)then
int_pot=int_pot+polynomial_interaction(vecji(4,xyzzyaaab23))
int_pot=int_pot-xyzzyaabq1*xyzzyaadp1(vecji(1:3,xyzzyaaab23),vecji(4,x&
&yzzyaaab23),.false.)
endif
enddo
if(.not.just_want_j_gt_i)int_pot=int_pot*0.5d0
int_pot=int_pot+self_term
end subroutine xyzzyaadv1
subroutine xyzzyaadw1(ie,vecji,int_pot,just_want_j_gt_i)
implicit none
integer,intent(in) :: ie
integer :: xyzzyaaaa24,xyzzyaaab24,xyzzyaaac24,xyzzyaaad24,xyzzyaaae24
real(dp) :: xyzzyaaaf24,xyzzyaaag24,xyzzyaaah24
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: int_pot
logical,intent(in) :: just_want_j_gt_i
logical xyzzyaaai24
integer xyzzyaaaj24
int_pot=0.d0
xyzzyaaah24=0.d0
xyzzyaaai24=periodicity==2
xyzzyaaad24=1
if(xyzzyaaai24)then
do xyzzyaaae24=2,xyzzyaaal1
xyzzyaaab24=first_s_in_star(xyzzyaaae24)
xyzzyaaac24=first_s_in_star(xyzzyaaae24+1)-xyzzyaaab24
xyzzyaaad24=xyzzyaaad24+xyzzyaaac24
xyzzyaaag24=0.d0
do xyzzyaaaa24=xyzzyaaab24,xyzzyaaab24+xyzzyaaac24-1
xyzzyaaag24=xyzzyaaag24+(1.d0-3.d0*sin(xyzzyaabr1)*sin(xyzzyaabr1)*((s&
&_lattice(1,xyzzyaaaa24)*s_lattice(1,xyzzyaaaa24))/(sum(s_lattice(1:di&
&mensionality,xyzzyaaab24)**2))))*sum(s_lattice(1:dimensionality,xyzzy&
&aaab24)**2)**(-1.5d0)
enddo
xyzzyaaah24=xyzzyaaah24+xyzzyaaag24
enddo
xyzzyaaaf24=sqrt(dble(xyzzyaaad24)*area*one_over_pi)
xyzzyaabo1=pi*(2.d0-3.d0*sin(xyzzyaabr1)*sin(xyzzyaabr1))/(area*xyzzya&
&aaf24)
xyzzyaabp1=0.75/xyzzyaaaf24**2
endif
if(.not.just_want_j_gt_i)then
do xyzzyaaaj24=1,ie-1
int_pot=int_pot+xyzzyaabq1*xyzzyaadq1(vecji(1:3,xyzzyaaaj24),vecji(4,x&
&yzzyaaaj24),xyzzyaabr1,xyzzyaaai24)
if(vecji(4,xyzzyaaaj24)<man_int_poly_cutoff)then
int_pot=int_pot+polynomial_interaction(vecji(4,xyzzyaaaj24))+tilted_po&
&lynomial_interaction(vecji(1,xyzzyaaaj24),vecji(4,xyzzyaaaj24))
int_pot=int_pot-xyzzyaabq1*xyzzyaadq1(vecji(1:3,xyzzyaaaj24),vecji(4,x&
&yzzyaaaj24),xyzzyaabr1,.false.)
endif
enddo
endif
do xyzzyaaaj24=ie+1,netot
int_pot=int_pot+xyzzyaabq1*xyzzyaadq1(vecji(1:3,xyzzyaaaj24),vecji(4,x&
&yzzyaaaj24),xyzzyaabr1,xyzzyaaai24)
if(vecji(4,xyzzyaaaj24)<man_int_poly_cutoff)then
int_pot=int_pot+polynomial_interaction(vecji(4,xyzzyaaaj24))+tilted_po&
&lynomial_interaction(vecji(1,xyzzyaaaj24),vecji(4,xyzzyaaaj24))
int_pot=int_pot-xyzzyaabq1*xyzzyaadq1(vecji(1:3,xyzzyaaaj24),vecji(4,x&
&yzzyaaaj24),xyzzyaabr1,.false.)
endif
enddo
if(.not.just_want_j_gt_i)int_pot=int_pot*0.5d0
xyzzyaaah24=0.5d0*xyzzyaabq1*(xyzzyaaah24+xyzzyaabo1)
if(xyzzyaaai24)int_pot=int_pot+xyzzyaaah24
end subroutine xyzzyaadw1
subroutine xyzzyaadx1(ie,vecji,int_pot,just_want_j_gt_i)
implicit none
integer,intent(in) :: ie
integer :: xyzzyaaaa25,xyzzyaaab25,xyzzyaaac25,xyzzyaaad25,xyzzyaaae25
real(dp) :: xyzzyaaaf25,xyzzyaaag25,xyzzyaaah25
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: int_pot
logical,intent(in) :: just_want_j_gt_i
logical xyzzyaaai25
integer xyzzyaaaj25
int_pot=0.d0
xyzzyaaah25=0.d0
xyzzyaaai25=periodicity==2
xyzzyaaad25=1
if(xyzzyaaai25)then
do xyzzyaaae25=2,xyzzyaaal1
xyzzyaaab25=first_s_in_star(xyzzyaaae25)
xyzzyaaac25=first_s_in_star(xyzzyaaae25+1)-xyzzyaaab25
xyzzyaaad25=xyzzyaaad25+xyzzyaaac25
xyzzyaaag25=0.d0
do xyzzyaaaa25=xyzzyaaab25,xyzzyaaab25+xyzzyaaac25-1
xyzzyaaag25=xyzzyaaag25+(1.d0-3.d0*sin(xyzzyaabr1)*sin(xyzzyaabr1)*((s&
&_lattice(1,xyzzyaaaa25)*s_lattice(1,xyzzyaaaa25))/(sum(s_lattice(1:di&
&mensionality,xyzzyaaab25)**2))))*sum(s_lattice(1:dimensionality,xyzzy&
&aaab25)**2)**(-1.5d0)
enddo
xyzzyaaah25=xyzzyaaah25+xyzzyaaag25
enddo
xyzzyaaaf25=sqrt(dble(xyzzyaaad25)*area*one_over_pi)
xyzzyaabo1=pi*(2.d0-3.d0*sin(xyzzyaabr1)*sin(xyzzyaabr1))/(area*xyzzya&
&aaf25)
xyzzyaabp1=0.75/xyzzyaaaf25**2
endif
if(.not.just_want_j_gt_i)then
do xyzzyaaaj25=1,ie-1
int_pot=int_pot+xyzzyaadq1(vecji(1:3,xyzzyaaaj25),vecji(4,xyzzyaaaj25)&
&,xyzzyaabr1,xyzzyaaai25)
enddo
endif
do xyzzyaaaj25=ie+1,netot
int_pot=int_pot+xyzzyaadq1(vecji(1:3,xyzzyaaaj25),vecji(4,xyzzyaaaj25)&
&,xyzzyaabr1,xyzzyaaai25)
enddo
xyzzyaaah25=0.5d0*xyzzyaabq1*(xyzzyaaah25+xyzzyaabo1)
int_pot=int_pot*xyzzyaabq1
if(.not.just_want_j_gt_i)int_pot=int_pot*0.5d0
if(xyzzyaaai25)int_pot=int_pot+xyzzyaaah25
end subroutine xyzzyaadx1
subroutine biex3pot(rele,biexpot)
use slaarnaas,only : me_biex3,mh_biex3,mu_biex3,xx_sep
implicit none
real(dp),intent(in) :: rele(3,netot)
real(dp),intent(out) :: biexpot
real(dp) xyzzyaaaa26,xyzzyaaab26,xyzzyaaac26,xyzzyaaad26,xyzzyaaae26,x&
&yzzyaaaf26
real(dp) xyzzyaaag26,xyzzyaaah26,xyzzyaaai26,xyzzyaaaj26,xyzzyaaak26,x&
&yzzyaaal26
real(dp) xyzzyaaam26,xyzzyaaan26,xyzzyaaao26,xyzzyaaap26,xyzzyaaaq26,x&
&yzzyaaar26,xyzzyaaas26,xyzzyaaat26
real(dp),save :: xyzzyaaau26,xyzzyaaav26,xyzzyaaaw26
logical,save :: xyzzyaaax26=.true.
if(xyzzyaaax26)then
xyzzyaaau26=mu_biex3/me_biex3
xyzzyaaav26=mu_biex3/mh_biex3
xyzzyaaaw26=abs(heg_zlayer(2)-heg_zlayer(1))
xyzzyaaax26=.false.
endif
xyzzyaaaa26=rele(1,1)
xyzzyaaad26=rele(2,1)
xyzzyaaai26=sqrt(xyzzyaaaa26*xyzzyaaaa26+xyzzyaaad26*xyzzyaaad26+xyzzy&
&aaaw26*xyzzyaaaw26)
xyzzyaaab26=rele(1,2)
xyzzyaaae26=rele(2,2)
xyzzyaaal26=sqrt(xyzzyaaab26*xyzzyaaab26+xyzzyaaae26*xyzzyaaae26+xyzzy&
&aaaw26*xyzzyaaaw26)
xyzzyaaac26=xyzzyaaaa26-xyzzyaaab26
xyzzyaaaf26=xyzzyaaad26-xyzzyaaae26
if(fix_holes)then
xyzzyaaam26=xx_sep-xyzzyaaac26
xyzzyaaan26=-xyzzyaaaf26
xyzzyaaao26=xx_sep
xyzzyaaap26=0.d0
xyzzyaaaq26=xx_sep-xyzzyaaaa26
xyzzyaaar26=xyzzyaaad26
xyzzyaaas26=xx_sep+xyzzyaaab26
xyzzyaaat26=xyzzyaaae26
else
xyzzyaaam26=xx_sep+xyzzyaaau26*xyzzyaaac26
xyzzyaaan26=xyzzyaaau26*xyzzyaaaf26
xyzzyaaao26=xx_sep+xyzzyaaav26*(-xyzzyaaac26)
xyzzyaaap26=xyzzyaaav26*(-xyzzyaaaf26)
xyzzyaaaq26=xx_sep+xyzzyaaau26*xyzzyaaaa26+xyzzyaaav26*xyzzyaaab26
xyzzyaaar26=xyzzyaaau26*xyzzyaaad26+xyzzyaaav26*xyzzyaaae26
xyzzyaaas26=xx_sep-xyzzyaaav26*xyzzyaaaa26-xyzzyaaau26*xyzzyaaab26
xyzzyaaat26=-xyzzyaaav26*xyzzyaaad26-xyzzyaaau26*xyzzyaaae26
endif
xyzzyaaag26=sqrt(xyzzyaaam26*xyzzyaaam26+xyzzyaaan26*xyzzyaaan26)
xyzzyaaah26=sqrt(xyzzyaaao26*xyzzyaaao26+xyzzyaaap26*xyzzyaaap26)
xyzzyaaaj26=sqrt(xyzzyaaaq26*xyzzyaaaq26+xyzzyaaar26*xyzzyaaar26+xyzzy&
&aaaw26*xyzzyaaaw26)
xyzzyaaak26=sqrt(xyzzyaaas26*xyzzyaaas26+xyzzyaaat26*xyzzyaaat26+xyzzy&
&aaaw26*xyzzyaaaw26)
biexpot=1.d0/xyzzyaaag26+1.d0/xyzzyaaah26-1.d0/xyzzyaaai26-1.d0/xyzzya&
&aal26-1.d0/xyzzyaaaj26-1.d0/xyzzyaaak26
end subroutine biex3pot
subroutine xyzzyaady1
implicit none
integer xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xyzzyaaad27,xyzzyaaae27,xy&
&zzyaaaf27
integer,parameter :: xyzzyaaag27=13
real(dp) xyzzyaaah27,xyzzyaaai27
real(dp),parameter :: xyzzyaaaj27=16.d0,xyzzyaaak27=0.016d0
real(dp),allocatable :: xyzzyaaal27(:),xyzzyaaam27(:),xyzzyaaan27(:),x&
&yzzyaaao27(:)
allocate(xyzzyaaal27(191),xyzzyaaam27(1001),xyzzyaaan27(1001),xyzzyaaa&
&o27(14014),stat=xyzzyaaaf27)
call check_alloc(xyzzyaaaf27,'FTABLE','')
do xyzzyaaab27=1,191
xyzzyaaal27(xyzzyaaab27)=1.d0/(xyzzyaaab27+xyzzyaaab27-1)
enddo
xyzzyaaay1=xyzzyaaak27
xyzzyaaaz1=1.d0/xyzzyaaay1
xyzzyaaba1=(xyzzyaaaj27+2.5d0)*1.5d0/(xyzzyaaaj27+4.5d0)
xyzzyaabb1=xyzzyaaba1+1.d0
xyzzyaabc1=xyzzyaaba1+1.5d0
xyzzyaaad27=xyzzyaaag27*1001
xyzzyaaab27=xyzzyaaag27+100
xyzzyaaah27=xyzzyaaal27(xyzzyaaab27+1)
do xyzzyaaaa27=1,1001
xyzzyaaai27=(xyzzyaaaa27-1)*xyzzyaaay1
xyzzyaaam27(xyzzyaaaa27)=exp(-xyzzyaaai27)
xyzzyaaan27(xyzzyaaaa27)=xyzzyaaai27+xyzzyaaai27
xyzzyaaao27(xyzzyaaad27+xyzzyaaaa27)=xyzzyaaan27(xyzzyaaaa27)*xyzzyaaa&
&h27
enddo
3 xyzzyaaah27=xyzzyaaal27(xyzzyaaab27)
do xyzzyaaaa27=1,1001
xyzzyaaao27(xyzzyaaad27+xyzzyaaaa27)=(xyzzyaaao27(xyzzyaaad27+xyzzyaaa&
&a27)*xyzzyaaan27(xyzzyaaaa27)+xyzzyaaam27(xyzzyaaaa27))*xyzzyaaah27
enddo
xyzzyaaab27=xyzzyaaab27-1
if(xyzzyaaab27/=xyzzyaaag27)goto 3
6 xyzzyaaac27=xyzzyaaad27-1001
xyzzyaaah27=xyzzyaaal27(xyzzyaaab27)
do xyzzyaaaa27=1,1001
xyzzyaaao27(xyzzyaaac27+xyzzyaaaa27)=(xyzzyaaao27(xyzzyaaad27+xyzzyaaa&
&a27)*xyzzyaaan27(xyzzyaaaa27)+xyzzyaaam27(xyzzyaaaa27))*xyzzyaaah27
enddo
xyzzyaaad27=xyzzyaaac27
xyzzyaaab27=xyzzyaaab27-1
if(xyzzyaaab27/=0)goto 6
do xyzzyaaae27=0,1000
xyzzyaabm1(3+4*xyzzyaaae27)=xyzzyaaao27(xyzzyaaae27+1)
xyzzyaabm1(2+4*xyzzyaaae27)=xyzzyaaao27(xyzzyaaae27+1002)
xyzzyaabm1(1+4*xyzzyaaae27)=xyzzyaaao27(xyzzyaaae27+2003)*0.5d0
xyzzyaabm1(4*xyzzyaaae27)=xyzzyaaao27(xyzzyaaae27+3004)*sixth
enddo
deallocate(xyzzyaaal27,xyzzyaaam27,xyzzyaaan27,xyzzyaaao27)
end subroutine xyzzyaady1
subroutine xyzzyaadz1(f,fz)
implicit none
real(dp),intent(out) :: f,fz
integer xyzzyaaaa28
real(dp) xyzzyaaab28,xyzzyaaac28,xyzzyaaad28,xyzzyaaae28
xyzzyaaab28=exp(-xyzzyaaat1)*xyzzyaabk1
xyzzyaaac28=xyzzyaaat1+xyzzyaaat1
if(xyzzyaaat1<16.d0)then
xyzzyaaaa28=int(xyzzyaaaz1*xyzzyaaat1+0.5d0)
xyzzyaaae28=xyzzyaaay1*xyzzyaaaa28-xyzzyaaat1
xyzzyaaad28=(((xyzzyaabm1(4*xyzzyaaaa28)*xyzzyaaae28+xyzzyaabm1(4*xyzz&
&yaaaa28+1))*xyzzyaaae28+xyzzyaabm1(4*xyzzyaaaa28+2))*xyzzyaaae28+xyzz&
&yaabm1(4*xyzzyaaaa28+3))
xyzzyaaad28=xyzzyaaad28*xyzzyaabk1
else
xyzzyaaad28=sqrt(pi_over_two/xyzzyaaac28)*xyzzyaabk1-((xyzzyaaat1+xyzz&
&yaabb1)*xyzzyaaab28/((xyzzyaaat1+xyzzyaabc1)*xyzzyaaac28+xyzzyaaba1))
endif
fz=xyzzyaaad28*xyzzyaaac28
f=fz+xyzzyaaab28
end subroutine xyzzyaadz1
subroutine xyzzyaaea1(phi,f)
implicit none
real(dp),intent(in) :: phi
real(dp),intent(out) :: f
integer xyzzyaaaa29
real(dp) xyzzyaaab29,xyzzyaaac29,xyzzyaaad29,xyzzyaaae29,xyzzyaaaf29,x&
&yzzyaaag29,xyzzyaaah29,xyzzyaaai29
xyzzyaaab29=phi*phi+xyzzyaaat1
xyzzyaaac29=exp(-xyzzyaaab29)*0.5d0
xyzzyaaad29=phi+xyzzyaaas1
xyzzyaaae29=xyzzyaaad29*xyzzyaaad29
xyzzyaaaf29=exp(xyzzyaaab29-xyzzyaaae29)
if(xyzzyaaae29<16.d0)then
xyzzyaaaa29=int(xyzzyaaaz1*xyzzyaaae29+0.5d0)
xyzzyaaai29=xyzzyaaay1*xyzzyaaaa29-xyzzyaaae29
xyzzyaaag29=(1.d0-(((xyzzyaabm1(4*xyzzyaaaa29)*xyzzyaaai29+xyzzyaabm1(&
&4*xyzzyaaaa29+1))*xyzzyaaai29+xyzzyaabm1(4*xyzzyaaaa29+2))*xyzzyaaai2&
&9+xyzzyaabm1(4*xyzzyaaaa29+3))*xyzzyaaad29*two_over_root_pi)/xyzzyaaa&
&f29
else
xyzzyaaag29=(xyzzyaaae29+xyzzyaabb1)*xyzzyaaac29*xyzzyaaad29*four_over&
&_root_pi/((xyzzyaaae29+xyzzyaabc1)*(xyzzyaaae29+xyzzyaaae29)+xyzzyaab&
&a1)
endif
xyzzyaaad29=phi-xyzzyaaas1
xyzzyaaae29=xyzzyaaad29*xyzzyaaad29
if(xyzzyaaae29<16.d0)then
xyzzyaaaa29=int(xyzzyaaaz1*xyzzyaaae29+0.5d0)
xyzzyaaai29=xyzzyaaay1*xyzzyaaaa29-xyzzyaaae29
xyzzyaaah29=(1.d0-(((xyzzyaabm1(4*xyzzyaaaa29)*xyzzyaaai29+xyzzyaabm1(&
&4*xyzzyaaaa29+1))*xyzzyaaai29+xyzzyaabm1(4*xyzzyaaaa29+2))*xyzzyaaai2&
&9+xyzzyaabm1(4*xyzzyaaaa29+3))*xyzzyaaad29*two_over_root_pi)*xyzzyaaa&
&f29
else
xyzzyaaah29=(xyzzyaaae29+xyzzyaabb1)*xyzzyaaac29*xyzzyaaaf29*xyzzyaaad&
&29*four_over_root_pi/((xyzzyaaae29+xyzzyaabc1)*(xyzzyaaae29+xyzzyaaae&
&29)+xyzzyaaba1)
if(xyzzyaaad29<(0.d0))xyzzyaaah29=xyzzyaaah29+2.d0
xyzzyaaah29=xyzzyaaah29*xyzzyaaaf29
endif
f=xyzzyaaag29+xyzzyaaah29
end subroutine xyzzyaaea1
subroutine xyzzyaaeb1(ie,b,vecji,mpc_short,mpc_long,just_want_j_gt_i,i&
&sinf)
use slaarnaas, only : half_mpc_constant
use slaarnabt, only : ddot
implicit none
integer,intent(in) :: ie,b
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: mpc_short,mpc_long
logical,intent(in) :: just_want_j_gt_i
logical,intent(out) :: isinf
integer xyzzyaaaa30,xyzzyaaab30,xyzzyaaac30,xyzzyaaad30
call timer('MPC',.true.)
call xyzzyaaec1(ie,vecji,mpc_short,just_want_j_gt_i,isinf)
if(.not.homogeneous_mpc)then
if(electron_system)then
if(blip_mpc.and.periodicity==3)then
mpc_long=xyzzyaaeg1(xyzzyaadc1)
else
mpc_long=ddot(xyzzyaaak1,xyzzyaacb1(1,1),1,xyzzyaacm1(1,b),1)-ddot(xyz&
&zyaaak1,xyzzyaacc1(1,1),1,xyzzyaacn1(1,b),1)
endif
else
xyzzyaaad30=0
do xyzzyaaaa30=1,nspin
xyzzyaaad30=xyzzyaaad30+nele(xyzzyaaaa30)
if(ie<=xyzzyaaad30)exit
enddo
xyzzyaaab30=which_fam(xyzzyaaaa30)
mpc_long=0.d0
do xyzzyaaac30=1,no_families
mpc_long=mpc_long+(ddot(xyzzyaaak1,xyzzyaacb1(1,xyzzyaaac30),1,xyzzyaa&
&cm1(1,b),1)-ddot(xyzzyaaak1,xyzzyaacc1(1,xyzzyaaac30),1,xyzzyaacn1(1,&
&b),1))*fam_charge(xyzzyaaac30)
enddo
mpc_long=mpc_long*fam_charge(xyzzyaaab30)
endif
else
mpc_long=half_mpc_constant
endif
call timer('MPC',.false.)
end subroutine xyzzyaaeb1
subroutine xyzzyaaec1(ie,vecji,direct,just_want_j_gt_i,isinf)
implicit none
integer,intent(in) :: ie
real(dp),intent(in) :: vecji(4,netot)
real(dp),intent(out) :: direct
logical,intent(in) :: just_want_j_gt_i
logical,intent(out) :: isinf
integer xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31,xyzzyaaad31,xyzzyaaae31,xy&
&zzyaaaf31
real(dp) xyzzyaaag31,xyzzyaaah31
real(dp),allocatable,save :: xyzzyaaai31(:,:)
logical,save :: xyzzyaaaj31=.true.
if(xyzzyaaaj31)then
if(heg_nlayers>1)then
allocate(xyzzyaaai31(nspin,nspin),stat=xyzzyaaaf31)
call check_alloc(xyzzyaaaf31,'COUL_EVAL_DIRECT','yzsq')
do xyzzyaaab31=1,nspin
xyzzyaaai31(xyzzyaaab31,xyzzyaaab31)=0.d0
do xyzzyaaac31=xyzzyaaab31+1,nspin
xyzzyaaai31(xyzzyaaac31,xyzzyaaab31)=(heg_ylayer(heg_layer(xyzzyaaac31&
&))-heg_ylayer(heg_layer(xyzzyaaab31)))**2+(heg_zlayer(heg_layer(xyzzy&
&aaac31))-heg_zlayer(heg_layer(xyzzyaaab31)))**2
xyzzyaaai31(xyzzyaaab31,xyzzyaaac31)=xyzzyaaai31(xyzzyaaac31,xyzzyaaab&
&31)
enddo
enddo
endif
xyzzyaaaj31=.false.
endif
direct=0.d0
isinf=.false.
if(heg_nlayers==1)then
if(electron_system)then
if(just_want_j_gt_i)then
do xyzzyaaaa31=ie+1,netot
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
direct=direct+1.d0/xyzzyaaah31
endif
enddo
else
do xyzzyaaaa31=1,ie-1
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
direct=direct+0.5d0/xyzzyaaah31
endif
enddo
do xyzzyaaaa31=ie+1,netot
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
direct=direct+0.5d0/xyzzyaaah31
endif
enddo
endif
else
xyzzyaaad31=0
do xyzzyaaab31=1,nspin
xyzzyaaad31=xyzzyaaad31+nele(xyzzyaaab31)
if(ie<=xyzzyaaad31)exit
enddo
if(just_want_j_gt_i)then
xyzzyaaag31=0.d0
do xyzzyaaaa31=ie+1,xyzzyaaad31
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
direct=direct+xyzzyaaag31*xyzzyaacd1(xyzzyaaab31,xyzzyaaab31)
do xyzzyaaac31=xyzzyaaab31+1,nspin
xyzzyaaae31=xyzzyaaad31+1
xyzzyaaad31=xyzzyaaad31+nele(xyzzyaaac31)
xyzzyaaag31=0.d0
do xyzzyaaaa31=xyzzyaaae31,xyzzyaaad31
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
direct=direct+xyzzyaaag31*xyzzyaacd1(xyzzyaaac31,xyzzyaaab31)
enddo
else
xyzzyaaad31=0
do xyzzyaaac31=1,nspin
xyzzyaaae31=xyzzyaaad31+1
xyzzyaaad31=xyzzyaaad31+nele(xyzzyaaac31)
xyzzyaaag31=0.d0
if(xyzzyaaac31==xyzzyaaab31)then
do xyzzyaaaa31=xyzzyaaae31,ie-1
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
do xyzzyaaaa31=ie+1,xyzzyaaad31
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
else
do xyzzyaaaa31=xyzzyaaae31,xyzzyaaad31
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
endif
direct=direct+xyzzyaaag31*xyzzyaacd1(xyzzyaaac31,xyzzyaaab31)
enddo
direct=0.5d0*direct
endif
endif
else
xyzzyaaad31=0
do xyzzyaaab31=1,nspin
xyzzyaaad31=xyzzyaaad31+nele(xyzzyaaab31)
if(ie<=xyzzyaaad31)exit
enddo
if(just_want_j_gt_i)then
xyzzyaaag31=0.d0
do xyzzyaaaa31=ie+1,xyzzyaaad31
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
direct=direct+xyzzyaaag31*xyzzyaacd1(xyzzyaaab31,xyzzyaaab31)
do xyzzyaaac31=xyzzyaaab31+1,nspin
xyzzyaaae31=xyzzyaaad31+1
xyzzyaaad31=xyzzyaaad31+nele(xyzzyaaac31)
xyzzyaaag31=0.d0
if(heg_layer(xyzzyaaab31)==heg_layer(xyzzyaaac31))then
do xyzzyaaaa31=xyzzyaaae31,xyzzyaaad31
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
else
do xyzzyaaaa31=xyzzyaaae31,xyzzyaaad31
xyzzyaaah31=sqrt(vecji(4,xyzzyaaaa31)**2+xyzzyaaai31(xyzzyaaac31,xyzzy&
&aaab31))
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
endif
direct=direct+xyzzyaaag31*xyzzyaacd1(xyzzyaaac31,xyzzyaaab31)
enddo
else
xyzzyaaad31=0
do xyzzyaaac31=1,nspin
xyzzyaaae31=xyzzyaaad31+1
xyzzyaaad31=xyzzyaaad31+nele(xyzzyaaac31)
xyzzyaaag31=0.d0
if(xyzzyaaac31==xyzzyaaab31)then
do xyzzyaaaa31=xyzzyaaae31,ie-1
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
do xyzzyaaaa31=ie+1,xyzzyaaad31
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
else
if(heg_layer(xyzzyaaab31)==heg_layer(xyzzyaaac31))then
do xyzzyaaaa31=xyzzyaaae31,xyzzyaaad31
xyzzyaaah31=vecji(4,xyzzyaaaa31)
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
else
do xyzzyaaaa31=xyzzyaaae31,xyzzyaaad31
xyzzyaaah31=sqrt(vecji(4,xyzzyaaaa31)**2+xyzzyaaai31(xyzzyaaac31,xyzzy&
&aaab31))
if(xyzzyaaah31==0.d0)then
isinf=.true.
else
xyzzyaaag31=xyzzyaaag31+1.d0/xyzzyaaah31
endif
enddo
endif
endif
direct=direct+xyzzyaaag31*xyzzyaacd1(xyzzyaaac31,xyzzyaaab31)
enddo
direct=0.5d0*direct
endif
endif
end subroutine xyzzyaaec1
subroutine xyzzyaaed1
use slaarnabt,only : dscal,approx_equal
implicit none
integer xyzzyaaaa32,xyzzyaaab32,xyzzyaaac32,xyzzyaaad32,xyzzyaaae32,xy&
&zzyaaaf32,xyzzyaaag32,xyzzyaaah32,xyzzyaaai32,xyzzyaaaj32,xyzzyaaak32&
&,xyzzyaaal32,xyzzyaaam32
real(dp) xyzzyaaan32,xyzzyaaao32,xyzzyaaap32,xyzzyaaaq32,xyzzyaaar32,x&
&yzzyaaas32,xyzzyaaat32,xyzzyaaau32,xyzzyaaav32,xyzzyaaaw32,xyzzyaaax3&
&2,xyzzyaaay32,xyzzyaaaz32(3),xyzzyaaba32(3),xyzzyaabb32(3),xyzzyaabc3&
&2,xyzzyaabd32,xyzzyaabe32,xyzzyaabf32,xyzzyaabg32,xyzzyaabh32,xyzzyaa&
&bi32,xyzzyaabj32,xyzzyaabk32(3),xyzzyaabl32(3),xyzzyaabm32(3),xyzzyaa&
&bn32(3),xyzzyaabo32,xyzzyaabp32,xyzzyaabq32,xyzzyaabr32,xyzzyaabs32,x&
&yzzyaabt32
real(dp),parameter :: xyzzyaabu32=1.d-8,xyzzyaabv32=1.d-6,xyzzyaabw32=&
&1.d-8,xyzzyaabx32=1.d-6
real(dp),allocatable :: xyzzyaaby32(:),xyzzyaabz32(:),xyzzyaaca32(:),x&
&yzzyaacb32(:),xyzzyaacc32(:)
complex(dp) xyzzyaacd32
character(80) tmpr
open(xyzzyaaai1,file='mpc.data',status='old',iostat=xyzzyaaae32,action&
&='read')
if(xyzzyaaae32/=0)call errstop('READ_MPC','Problem opening mpc.data fi&
&le.')
call skip(xyzzyaaai1,9)
read(xyzzyaaai1,*,err=6,end=6)xyzzyaabk32
read(xyzzyaaai1,*,err=6,end=6)xyzzyaabl32
read(xyzzyaaai1,*,err=6,end=6)xyzzyaabm32
if(any(abs(pa1-xyzzyaabk32)>xyzzyaabu32).or.any(abs(pa2-xyzzyaabl32)>x&
&yzzyaabu32).or.any(abs(pa3-xyzzyaabm32)>xyzzyaabu32))call errstop('RE&
&AD_MPC','Primitive lattice vectors in DENSITY section of mpc.data inc&
&ompatible with those in xwfn.data.')
read(xyzzyaaai1,*,err=6,end=6)
read(xyzzyaaai1,*,err=6,end=6)xyzzyaaaa32
if(xyzzyaaaa32/=nbasis)call errstop('READ_MPC','Mismatch in number of &
&atoms per primitive cell between mpc.data/xwfn.data.')
read(xyzzyaaai1,*,err=6,end=6)
do xyzzyaaaa32=1,nbasis
read(xyzzyaaai1,*,err=6,end=6)xyzzyaaac32,xyzzyaabn32
if(xyzzyaaac32/=atno(xyzzyaaaa32))call errstop('READ_MPC','Mismatch in&
& atomic numbers between mpc.data/xwfn.data.')
if(any(abs(basis(:,xyzzyaaaa32)-xyzzyaabn32)>xyzzyaabu32))call errstop&
&('READ_MPC','Mismatch between atomic positions in xwfn.data and those&
& in mpc.data.')
enddo
call skip(xyzzyaaai1,3)
read(xyzzyaaai1,*,err=6,end=6)xyzzyaacj1
if(xyzzyaacj1/=no_families)call errstop('READ_MPC','Number of particle&
& types in mpc.data does not match computed number of spin-families.')
if(qmc_density_mpc.and.den_nsets/=xyzzyaacj1)call errstop('READ_DENSIT&
&Y','No. of density sets for different particle types not equal in den&
&sity.data and expval.data.')
read(xyzzyaaai1,*,err=6,end=6)
read(xyzzyaaai1,*,err=6,end=6)xyzzyaaci1
allocate(den_gvec(3,xyzzyaaci1),den_sc(xyzzyaaci1,xyzzyaacj1),xyzzyaac&
&l1(3,xyzzyaaci1),stat=xyzzyaaad32)
call check_alloc(xyzzyaaad32,'READ_MPC','density arrays')
allocate(xyzzyaacm1(xyzzyaaci1,1),xyzzyaacn1(xyzzyaaci1,1),stat=xyzzya&
&aad32)
call check_alloc(xyzzyaaad32,'READ_MPC','den_cosgdotr array')
xyzzyaacm1=0.d0
xyzzyaacn1=0.d0
read(xyzzyaaai1,*,err=6,end=6)
do xyzzyaaaa32=1,xyzzyaaci1
read(xyzzyaaai1,*,err=6,end=6)den_gvec(1:3,xyzzyaaaa32)
enddo
do xyzzyaaab32=1,xyzzyaacj1
call skip(xyzzyaaai1,2)
read(xyzzyaaai1,*,err=6,end=6)xyzzyaaaa32
if(xyzzyaaaa32/=xyzzyaaab32)call errstop('READ_MPC','Error reading par&
&ticle type of set '//trim(i2s(xyzzyaaab32))//'.')
read(xyzzyaaai1,*,err=6,end=6)
do xyzzyaaaa32=1,xyzzyaaci1
read(xyzzyaaai1,*,err=6,end=6)xyzzyaaaw32,xyzzyaaax32
den_sc(xyzzyaaaa32,xyzzyaaab32)=cmplx(xyzzyaaaw32,xyzzyaaax32,dp)
enddo
enddo
call skip(xyzzyaaai1,3)
call skip(xyzzyaaai1,2)
read(xyzzyaaai1,*,err=7,end=7)xyzzyaaaz32
read(xyzzyaaai1,*,err=7,end=7)xyzzyaaba32
read(xyzzyaaai1,*,err=7,end=7)xyzzyaabb32
xyzzyaaan32=sqrt(a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3))
xyzzyaaaq32=sqrt(xyzzyaaaz32(1)*xyzzyaaaz32(1)+xyzzyaaaz32(2)*xyzzyaaa&
&z32(2)+xyzzyaaaz32(3)*xyzzyaaaz32(3))
xyzzyaaao32=sqrt(a2(1)*a2(1)+a2(2)*a2(2)+a2(3)*a2(3))
xyzzyaaar32=sqrt(xyzzyaaba32(1)*xyzzyaaba32(1)+xyzzyaaba32(2)*xyzzyaab&
&a32(2)+xyzzyaaba32(3)*xyzzyaaba32(3))
xyzzyaaap32=sqrt(a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3))
xyzzyaaas32=sqrt(xyzzyaabb32(1)*xyzzyaabb32(1)+xyzzyaabb32(2)*xyzzyaab&
&b32(2)+xyzzyaabb32(3)*xyzzyaabb32(3))
xyzzyaabc32=dot_product(a1,xyzzyaaaz32)
xyzzyaabd32=dot_product(a2,xyzzyaaba32)
xyzzyaabe32=dot_product(a3,xyzzyaabb32)
if(abs(xyzzyaabc32-xyzzyaaan32*xyzzyaaaq32)>xyzzyaabu32.or.abs(xyzzyaa&
&bd32-xyzzyaaao32*xyzzyaaar32)>xyzzyaabu32.or.abs(xyzzyaabe32-xyzzyaaa&
&p32*xyzzyaaas32)>xyzzyaabu32)then
call wordwrap('Lattice vectors in EEPOT DATA section of mpc.data don''&
&t point in same direction as those in xwfn.data. Regenerate mpc.data &
&file.')
call wout('(x)wfn.data:')
call wout('          ',a1,rfmt='(f18.12)')
call wout('          ',a2,rfmt='(f18.12)')
call wout('          ',a3,rfmt='(f18.12)')
call wout('mpc.data:')
call wout('          ',xyzzyaaaz32,rfmt='(f18.12)')
call wout('          ',xyzzyaaba32,rfmt='(f18.12)')
call wout('          ',xyzzyaabb32,rfmt='(f18.12)')
call wout()
call errstop('READ_MPC','Quitting')
endif
if(abs(xyzzyaaaq32)<xyzzyaabu32.or.abs(xyzzyaaar32)<xyzzyaabu32.or.abs&
&(xyzzyaaas32)<xyzzyaabu32)then
call wordwrap('Length of one or more lattice vectors in EEPOT DATA sec&
&tion of mpc.data file is zero. This will lead to a division by zero e&
&rror.')
call errstop('READ_MPC','Quitting.')
endif
xyzzyaabf32=xyzzyaaan32/xyzzyaaaq32
xyzzyaabg32=xyzzyaaao32/xyzzyaaar32
xyzzyaabh32=xyzzyaaap32/xyzzyaaas32
if(periodicity==2)then
if(abs(xyzzyaabf32-xyzzyaabg32)>xyzzyaabu32)then
call wout()
call wordwrap('Lattice vectors in EEPOT section of mpc.data need to be&
& obtainable from lattice vectors in (x)wfn.data through scaling by a &
&*single* parameter.  This isn''t the case.  You need to regenerate mp&
&c.data.')
call wout('(x)wfn.data:')
call wout('          ',a1(1:2),rfmt='(f18.12)')
call wout('          ',a2(1:2),rfmt='(f18.12)')
call wout('mpc.data:')
call wout('          ',xyzzyaaaz32(1:2),rfmt='(f18.12)')
call wout('          ',xyzzyaaba32(1:2),rfmt='(f18.12)')
call wout()
call errstop('READ_MPC','Quitting.')
endif
if(abs(xyzzyaabh32-1.d0)>xyzzyaabu32)then
call wordwrap('The third lattice vector in quasi-2D systems must be th&
&e same in the EEPOT section of mpc.data and (x)wfn.data.  It isn''t.'&
&)
call wout('(x)wfn.data:')
call wout('          ',a3(3),rfmt='(f18.12)')
call wout('mpc.data:')
call wout('          ',xyzzyaabb32(3),rfmt='(f18.12)')
call wout()
call errstop('READ_MPC','Quitting.')
endif
else
if(abs(xyzzyaabf32-xyzzyaabg32)>xyzzyaabu32.or.abs(xyzzyaabf32-xyzzyaa&
&bh32)>xyzzyaabu32)then
call wout()
call wordwrap('Lattice vectors in EEPOT section of mpc.data need to be&
& obtainable from lattice vectors in (x)wfn.data through scaling by a &
&*single* parameter.  This isn''t the case.  You need to regenerate mp&
&c.data.')
call wout('(x)wfn.data:')
call wout('          ',a1,rfmt='(f18.12)')
call wout('          ',a2,rfmt='(f18.12)')
call wout('          ',a3,rfmt='(f18.12)')
call wout('mpc.data:')
call wout('          ',xyzzyaaaz32,rfmt='(f18.12)')
call wout('          ',xyzzyaaba32,rfmt='(f18.12)')
call wout('          ',xyzzyaabb32,rfmt='(f18.12)')
call wout()
call errstop('READ_MPC','Quitting.')
endif
endif
call skip(xyzzyaaai1,5)
read(xyzzyaaai1,*,err=7,end=7)xyzzyaaaf32
allocate(xyzzyaabz32(xyzzyaaaf32),xyzzyaaca32(xyzzyaaaf32),xyzzyaacb32&
&(xyzzyaaaf32),xyzzyaacc32(xyzzyaaaf32),stat=xyzzyaaad32)
call check_alloc(xyzzyaaad32,'READ_MPC','EEPOT G and f(G)')
call skip(xyzzyaaai1,1)
do xyzzyaaaa32=1,xyzzyaaaf32
read(xyzzyaaai1,*,err=7,end=7)xyzzyaabz32(xyzzyaaaa32),xyzzyaaca32(xyz&
&zyaaaa32),xyzzyaacb32(xyzzyaaaa32),xyzzyaacc32(xyzzyaaaa32)
enddo
xyzzyaaan32=a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3)
xyzzyaaao32=xyzzyaaaz32(1)*xyzzyaaaz32(1)+xyzzyaaaz32(2)*xyzzyaaaz32(2&
&)+xyzzyaaaz32(3)*xyzzyaaaz32(3)
xyzzyaaat32=xyzzyaaan32/xyzzyaaao32
xyzzyaaau32=sqrt(xyzzyaaao32/xyzzyaaan32)
if(xyzzyaaat32-1.d0>xyzzyaabu32)then
call wout('Read mpc.data file')
tmpr=r2s(xyzzyaaat32,'(es24.16)')
call wout('Rescaling EEPOT G vectors by            : '//trim(tmpr))
tmpr=r2s(xyzzyaaau32,'(es24.16)')
call wout('Rescaling EEPOT Fourier coefficients by : '//trim(tmpr))
call wout()
call dscal(xyzzyaaaf32,xyzzyaaau32,xyzzyaabz32(1),1)
call dscal(xyzzyaaaf32,xyzzyaaau32,xyzzyaaca32(1),1)
call dscal(xyzzyaaaf32,xyzzyaaau32,xyzzyaacb32(1),1)
call dscal(xyzzyaaaf32,xyzzyaaat32,xyzzyaacc32(1),1)
else
call wout()
call wout('Successfully read mpc.data file. Rescaling of data not requ&
&ired.')
endif
close(xyzzyaaai1)
xyzzyaaaa32=-1
do xyzzyaaac32=1,xyzzyaacj1
do xyzzyaaab32=xyzzyaaci1,1,-1
if(abs(den_sc(xyzzyaaab32,xyzzyaaac32))>xyzzyaabv32)exit
enddo
xyzzyaaaa32=max(xyzzyaaaa32,xyzzyaaab32)
enddo
xyzzyaabi32=den_gvec(1,xyzzyaaaa32)*den_gvec(1,xyzzyaaaa32)+den_gvec(2&
&,xyzzyaaaa32)*den_gvec(2,xyzzyaaaa32)+den_gvec(3,xyzzyaaaa32)*den_gve&
&c(3,xyzzyaaaa32)
xyzzyaack1=xyzzyaaci1
do xyzzyaaab32=xyzzyaaaa32+1,xyzzyaaci1
xyzzyaabj32=den_gvec(1,xyzzyaaab32)*den_gvec(1,xyzzyaaab32)+den_gvec(2&
&,xyzzyaaab32)*den_gvec(2,xyzzyaaab32)+den_gvec(3,xyzzyaaab32)*den_gve&
&c(3,xyzzyaaab32)
if(abs(xyzzyaabj32-xyzzyaabi32)>xyzzyaabw32)then
xyzzyaack1=xyzzyaaab32-1
exit
endif
enddo
xyzzyaach1=0
do xyzzyaaaa32=1,xyzzyaaci1
xyzzyaaag32=nint(one_over_twopi*(den_gvec(1,xyzzyaaaa32)*pa1(1)+den_gv&
&ec(2,xyzzyaaaa32)*pa1(2)+den_gvec(3,xyzzyaaaa32)*pa1(3)))
if(abs(xyzzyaaag32)>xyzzyaach1)xyzzyaach1=abs(xyzzyaaag32)
xyzzyaacl1(1,xyzzyaaaa32)=xyzzyaaag32
xyzzyaaah32=nint(one_over_twopi*(den_gvec(1,xyzzyaaaa32)*pa2(1)+den_gv&
&ec(2,xyzzyaaaa32)*pa2(2)+den_gvec(3,xyzzyaaaa32)*pa2(3)))
if(abs(xyzzyaaah32)>xyzzyaach1)xyzzyaach1=abs(xyzzyaaah32)
xyzzyaacl1(2,xyzzyaaaa32)=xyzzyaaah32
xyzzyaaai32=nint(one_over_twopi*(den_gvec(1,xyzzyaaaa32)*pa3(1)+den_gv&
&ec(2,xyzzyaaaa32)*pa3(2)+den_gvec(3,xyzzyaaaa32)*pa3(3)))
if(abs(xyzzyaaai32)>xyzzyaach1)xyzzyaach1=abs(xyzzyaaai32)
xyzzyaacl1(3,xyzzyaaaa32)=xyzzyaaai32
enddo
if(xyzzyaach1<1)call errstop('READ_MPC','G vector sets for density in &
&mpc.data must include G/=0')
if(permit_den_symm.and.qmc_density_mpc)then
if(expval_complex_den)then
do xyzzyaaak32=1,den_nsets
xyzzyaaac32=1
do
xyzzyaacd32=czero
xyzzyaaal32=0
xyzzyaabi32=den_gvec(1,xyzzyaaac32)*den_gvec(1,xyzzyaaac32)+den_gvec(2&
&,xyzzyaaac32)*den_gvec(2,xyzzyaaac32)+den_gvec(3,xyzzyaaac32)*den_gve&
&c(3,xyzzyaaac32)
do
xyzzyaabo32=abs(real(expval_den_mpc(xyzzyaaac32+xyzzyaaal32,xyzzyaaak3&
&2),dp))
xyzzyaabp32=abs(aimag(expval_den_mpc(xyzzyaaac32+xyzzyaaal32,xyzzyaaak&
&32)))
xyzzyaacd32=xyzzyaacd32+cmplx(xyzzyaabo32,xyzzyaabp32,dp)
xyzzyaaal32=xyzzyaaal32+1
if(xyzzyaaac32+xyzzyaaal32>xyzzyaaci1)exit
xyzzyaabj32=den_gvec(1,xyzzyaaac32+xyzzyaaal32)*den_gvec(1,xyzzyaaac32&
&+xyzzyaaal32)+den_gvec(2,xyzzyaaac32+xyzzyaaal32)*den_gvec(2,xyzzyaaa&
&c32+xyzzyaaal32)+den_gvec(3,xyzzyaaac32+xyzzyaaal32)*den_gvec(3,xyzzy&
&aaac32+xyzzyaaal32)
if(.not.approx_equal(xyzzyaabi32,xyzzyaabj32,xyzzyaabx32))exit
xyzzyaabq32=dble(den_sc(xyzzyaaac32+xyzzyaaal32,xyzzyaaak32))
xyzzyaabr32=dble(den_sc(xyzzyaaac32,xyzzyaaak32))
xyzzyaabs32=aimag(den_sc(xyzzyaaac32+xyzzyaaal32,xyzzyaaak32))
xyzzyaabt32=aimag(den_sc(xyzzyaaac32,xyzzyaaak32))
if(.not.approx_equal(abs(xyzzyaabq32),abs(xyzzyaabr32),xyzzyaabx32).or&
&..not.approx_equal(abs(xyzzyaabs32),abs(xyzzyaabt32),xyzzyaabx32))exi&
&t
enddo
if(xyzzyaaal32>1)xyzzyaacd32=xyzzyaacd32/real(xyzzyaaal32,dp)
do xyzzyaaab32=xyzzyaaac32,xyzzyaaac32+xyzzyaaal32-1
xyzzyaabo32=sign(1.d0,real(expval_den_mpc(xyzzyaaab32,xyzzyaaak32),dp)&
&)*real(xyzzyaacd32,dp)
xyzzyaabp32=sign(1.d0,aimag(expval_den_mpc(xyzzyaaab32,xyzzyaaak32)))*&
&aimag(xyzzyaacd32)
expval_den_mpc(xyzzyaaab32,xyzzyaaak32)=cmplx(xyzzyaabo32,xyzzyaabp32,&
&dp)
enddo
xyzzyaaac32=xyzzyaaac32+xyzzyaaal32
if(xyzzyaaac32>xyzzyaack1.or.xyzzyaaac32>expval_ngvec_truncated(xyzzya&
&aak32))exit
enddo
enddo
else
do xyzzyaaak32=1,den_nsets
xyzzyaaac32=1
do
xyzzyaacd32=czero
xyzzyaaal32=0
xyzzyaabi32=den_gvec(1,xyzzyaaac32)*den_gvec(1,xyzzyaaac32)+den_gvec(2&
&,xyzzyaaac32)*den_gvec(2,xyzzyaaac32)+den_gvec(3,xyzzyaaac32)*den_gve&
&c(3,xyzzyaaac32)
do
xyzzyaabo32=abs(real(expval_den_mpc(xyzzyaaac32+xyzzyaaal32,xyzzyaaak3&
&2),dp))
xyzzyaabp32=0.d0
xyzzyaacd32=xyzzyaacd32+cmplx(xyzzyaabo32,xyzzyaabp32,dp)
xyzzyaaal32=xyzzyaaal32+1
if(xyzzyaaac32+xyzzyaaal32>xyzzyaaci1)exit
xyzzyaabj32=den_gvec(1,xyzzyaaac32+xyzzyaaal32)*den_gvec(1,xyzzyaaac32&
&+xyzzyaaal32)+den_gvec(2,xyzzyaaac32+xyzzyaaal32)*den_gvec(2,xyzzyaaa&
&c32+xyzzyaaal32)+den_gvec(3,xyzzyaaac32+xyzzyaaal32)*den_gvec(3,xyzzy&
&aaac32+xyzzyaaal32)
if(.not.approx_equal(xyzzyaabi32,xyzzyaabj32,xyzzyaabx32))exit
xyzzyaabq32=dble(den_sc(xyzzyaaac32+xyzzyaaal32,xyzzyaaak32))
xyzzyaabr32=dble(den_sc(xyzzyaaac32,xyzzyaaak32))
xyzzyaabs32=0.d0
xyzzyaabt32=0.d0
if(.not.approx_equal(abs(xyzzyaabq32),abs(xyzzyaabr32),xyzzyaabx32))ex&
&it
enddo
if(xyzzyaaal32>1)xyzzyaacd32=xyzzyaacd32/real(xyzzyaaal32,dp)
do xyzzyaaab32=xyzzyaaac32,xyzzyaaac32+xyzzyaaal32-1
xyzzyaabo32=sign(1.d0,real(expval_den_mpc(xyzzyaaab32,xyzzyaaak32),dp)&
&)*real(xyzzyaacd32,dp)
xyzzyaabp32=0.d0
expval_den_mpc(xyzzyaaab32,xyzzyaaak32)=cmplx(xyzzyaabo32,xyzzyaabp32,&
&dp)
enddo
xyzzyaaac32=xyzzyaaac32+xyzzyaaal32
if(xyzzyaaac32>xyzzyaack1.or.xyzzyaaac32>expval_ngvec_truncated(xyzzya&
&aak32))exit
enddo
enddo
endif
endif
allocate(xyzzyaaco1(3,-xyzzyaach1:xyzzyaach1),stat=xyzzyaaad32)
call check_alloc(xyzzyaaad32,'READ_MPC','den_mwork')
if(qmc_density_mpc)then
xyzzyaaak1=expval_ngvec_truncated(1)
else
xyzzyaaak1=xyzzyaack1
endif
allocate(xyzzyaacb1(xyzzyaaak1,xyzzyaacj1),xyzzyaacc1(xyzzyaaak1,xyzzy&
&aacj1),stat=xyzzyaaad32)
call check_alloc(xyzzyaaad32,'READ_MPC','mpc_gspace')
allocate(xyzzyaaby32(xyzzyaaak1),stat=xyzzyaaad32)
call check_alloc(xyzzyaaad32,'READ_MPC','g2')
xyzzyaaay32=1.d0/pvolume
mpc_correction=0.d0
mpc_hartree=0.d0
do xyzzyaaaj32=1,xyzzyaacj1
select case (periodicity)
case(2)
if(qmc_density_mpc)then
xyzzyaacb1(1,xyzzyaaaj32)=xyzzyaacc32(1)*dble(expval_den_mpc(1,xyzzyaa&
&aj32))
xyzzyaacc1(1,xyzzyaaaj32)=-xyzzyaacc32(1)*aimag(expval_den_mpc(1,xyzzy&
&aaaj32))
mpc_correction=mpc_correction+0.5d0*(xyzzyaacb1(1,xyzzyaaaj32)*dble(ex&
&pval_den_mpc(1,xyzzyaaaj32))-xyzzyaacc1(1,xyzzyaaaj32)*aimag(expval_d&
&en_mpc(1,xyzzyaaaj32)))
else
xyzzyaacb1(1,xyzzyaaaj32)=xyzzyaacc32(1)*dble(den_sc(1,xyzzyaaaj32))
xyzzyaacc1(1,xyzzyaaaj32)=-xyzzyaacc32(1)*aimag(den_sc(1,xyzzyaaaj32))
mpc_correction=mpc_correction+0.5d0*(xyzzyaacb1(1,xyzzyaaaj32)*dble(de&
&n_sc(1,xyzzyaaaj32))-xyzzyaacc1(1,xyzzyaaaj32)*aimag(den_sc(1,xyzzyaa&
&aj32)))
endif
case(3)
if(qmc_density_mpc)then
xyzzyaacb1(1,xyzzyaaaj32)=-xyzzyaaay32*xyzzyaacc32(1)*dble(expval_den_&
&mpc(1,xyzzyaaaj32))
xyzzyaacc1(1,xyzzyaaaj32)=xyzzyaaay32*xyzzyaacc32(1)*aimag(expval_den_&
&mpc(1,xyzzyaaaj32))
mpc_correction=mpc_correction+0.5d0*(xyzzyaacb1(1,xyzzyaaaj32)*dble(ex&
&pval_den_mpc(1,xyzzyaaaj32))-xyzzyaacc1(1,xyzzyaaaj32)*aimag(expval_d&
&en_mpc(1,xyzzyaaaj32)))
else
xyzzyaacb1(1,xyzzyaaaj32)=-xyzzyaaay32*xyzzyaacc32(1)*dble(den_sc(1,xy&
&zzyaaaj32))
xyzzyaacc1(1,xyzzyaaaj32)=xyzzyaaay32*xyzzyaacc32(1)*aimag(den_sc(1,xy&
&zzyaaaj32))
mpc_correction=mpc_correction+0.5d0*(xyzzyaacb1(1,xyzzyaaaj32)*dble(de&
&n_sc(1,xyzzyaaaj32))-xyzzyaacc1(1,xyzzyaaaj32)*aimag(den_sc(1,xyzzyaa&
&aj32)))
endif
case default
call errstop('READ_MPC','Dimensionality error.')
end select
if(qmc_density_mpc)then
do xyzzyaaaa32=1,xyzzyaaak1
xyzzyaaby32(xyzzyaaaa32)=expval_gvec(1,xyzzyaaaa32,1)*expval_gvec(1,xy&
&zzyaaaa32,1)+expval_gvec(2,xyzzyaaaa32,1)*expval_gvec(2,xyzzyaaaa32,1&
&)+expval_gvec(3,xyzzyaaaa32,1)*expval_gvec(3,xyzzyaaaa32,1)
enddo
else
do xyzzyaaaa32=1,xyzzyaaak1
xyzzyaaby32(xyzzyaaaa32)=den_gvec(1,xyzzyaaaa32)*den_gvec(1,xyzzyaaaa3&
&2)+den_gvec(2,xyzzyaaaa32)*den_gvec(2,xyzzyaaaa32)+den_gvec(3,xyzzyaa&
&aa32)*den_gvec(3,xyzzyaaaa32)
enddo
endif
xyzzyaaam32=0
select case (periodicity)
case(2)
a : do xyzzyaaaa32=2,xyzzyaaak1
do xyzzyaaab32=2,xyzzyaaaf32
if(qmc_density_mpc)then
if(abs(expval_gvec(1,xyzzyaaaa32,1)-xyzzyaabz32(xyzzyaaab32))<xyzzyaab&
&u32.and.abs(expval_gvec(2,xyzzyaaaa32,1)-xyzzyaaca32(xyzzyaaab32))<xy&
&zzyaabu32.and.abs(expval_gvec(3,xyzzyaaaa32,1)-xyzzyaacb32(xyzzyaaab3&
&2))<xyzzyaabu32)then
xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj32)=xyzzyaacc32(xyzzyaaab32)*dble(expv&
&al_den_mpc(xyzzyaaaa32,xyzzyaaaj32))
xyzzyaacc1(xyzzyaaaa32,xyzzyaaaj32)=-xyzzyaacc32(xyzzyaaab32)*aimag(ex&
&pval_den_mpc(xyzzyaaaa32,xyzzyaaaj32))
mpc_correction=mpc_correction+0.5d0*(xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj3&
&2)*dble(expval_den_mpc(xyzzyaaaa32,xyzzyaaaj32))-xyzzyaacc1(xyzzyaaaa&
&32,xyzzyaaaj32)*aimag(expval_den_mpc(xyzzyaaaa32,xyzzyaaaj32)))
cycle a
endif
else
if(abs(den_gvec(1,xyzzyaaaa32)-xyzzyaabz32(xyzzyaaab32))<xyzzyaabu32.a&
&nd.abs(den_gvec(2,xyzzyaaaa32)-xyzzyaaca32(xyzzyaaab32))<xyzzyaabu32.&
&and.abs(den_gvec(3,xyzzyaaaa32)-xyzzyaacb32(xyzzyaaab32))<xyzzyaabu32&
&)then
xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj32)=xyzzyaacc32(xyzzyaaab32)*dble(den_&
&sc(xyzzyaaaa32,xyzzyaaaj32))
xyzzyaacc1(xyzzyaaaa32,xyzzyaaaj32)=-xyzzyaacc32(xyzzyaaab32)*aimag(de&
&n_sc(xyzzyaaaa32,xyzzyaaaj32))
mpc_correction=mpc_correction+0.5d0*(xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj3&
&2)*dble(den_sc(xyzzyaaaa32,xyzzyaaaj32))-xyzzyaacc1(xyzzyaaaa32,xyzzy&
&aaaj32)*aimag(den_sc(xyzzyaaaa32,xyzzyaaaj32)))
cycle a
endif
endif
enddo
xyzzyaaam32=xyzzyaaam32+1
xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj32)=0.d0
xyzzyaacc1(xyzzyaaaa32,xyzzyaaaj32)=0.d0
enddo a
if(xyzzyaaam32>0)call wout('Number of primitive G vectors with no matc&
&hing supercell G-vector: '//trim(i2s(xyzzyaaam32)))
case(3)
b : do xyzzyaaaa32=2,xyzzyaaak1
do xyzzyaaab32=2,xyzzyaaaf32
if(qmc_density_mpc)then
if(abs(expval_gvec(1,xyzzyaaaa32,1)-xyzzyaabz32(xyzzyaaab32))<xyzzyaab&
&u32.and.abs(expval_gvec(2,xyzzyaaaa32,1)-xyzzyaaca32(xyzzyaaab32))<xy&
&zzyaabu32.and.abs(expval_gvec(3,xyzzyaaaa32,1)-xyzzyaacb32(xyzzyaaab3&
&2))<xyzzyaabu32)then
xyzzyaaav32=(fourpi/xyzzyaaby32(xyzzyaaaa32)-xyzzyaacc32(xyzzyaaab32))&
&*xyzzyaaay32
xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj32)=xyzzyaaav32*dble(expval_den_mpc(xy&
&zzyaaaa32,xyzzyaaaj32))
xyzzyaacc1(xyzzyaaaa32,xyzzyaaaj32)=-xyzzyaaav32*aimag(expval_den_mpc(&
&xyzzyaaaa32,xyzzyaaaj32))
mpc_correction=mpc_correction+0.5d0*(xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj3&
&2)*dble(expval_den_mpc(xyzzyaaaa32,xyzzyaaaj32))-xyzzyaacc1(xyzzyaaaa&
&32,xyzzyaaaj32)*aimag(expval_den_mpc(xyzzyaaaa32,xyzzyaaaj32)))
cycle b
endif
else
if(abs(den_gvec(1,xyzzyaaaa32)-xyzzyaabz32(xyzzyaaab32))<xyzzyaabu32.a&
&nd.abs(den_gvec(2,xyzzyaaaa32)-xyzzyaaca32(xyzzyaaab32))<xyzzyaabu32.&
&and.abs(den_gvec(3,xyzzyaaaa32)-xyzzyaacb32(xyzzyaaab32))<xyzzyaabu32&
&)then
xyzzyaaav32=(fourpi/xyzzyaaby32(xyzzyaaaa32)-xyzzyaacc32(xyzzyaaab32))&
&*xyzzyaaay32
xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj32)=xyzzyaaav32*dble(den_sc(xyzzyaaaa3&
&2,xyzzyaaaj32))
xyzzyaacc1(xyzzyaaaa32,xyzzyaaaj32)=-xyzzyaaav32*aimag(den_sc(xyzzyaaa&
&a32,xyzzyaaaj32))
mpc_correction=mpc_correction+0.5d0*(xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj3&
&2)*dble(den_sc(xyzzyaaaa32,xyzzyaaaj32))-xyzzyaacc1(xyzzyaaaa32,xyzzy&
&aaaj32)*aimag(den_sc(xyzzyaaaa32,xyzzyaaaj32)))
cycle b
endif
endif
enddo
xyzzyaaam32=xyzzyaaam32+1
xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj32)=0.d0
xyzzyaacc1(xyzzyaaaa32,xyzzyaaaj32)=0.d0
enddo b
if(xyzzyaaam32>0)call wout('Number of primitive G vectors with no matc&
&hing supercell G-vector: '//trim(i2s(xyzzyaaam32)))
if(hartree_xc)then
c :  do xyzzyaaaa32=2,xyzzyaaak1
do xyzzyaaab32=2,xyzzyaaaf32
if(qmc_density_mpc)then
if(abs(expval_gvec(1,xyzzyaaaa32,1)-xyzzyaabz32(xyzzyaaab32))<xyzzyaab&
&u32.and.abs(expval_gvec(2,xyzzyaaaa32,1)-xyzzyaaca32(xyzzyaaab32))<xy&
&zzyaabu32.and.abs(expval_gvec(3,xyzzyaaaa32,1)-xyzzyaacb32(xyzzyaaab3&
&2))<xyzzyaabu32)then
xyzzyaaav32=fourpi*xyzzyaaay32/xyzzyaaby32(xyzzyaaaa32)
mpc_hartree=mpc_hartree+xyzzyaaav32*(dble(expval_den_mpc(xyzzyaaaa32,x&
&yzzyaaaj32))**2+aimag(expval_den_mpc(xyzzyaaaa32,xyzzyaaaj32))**2)
cycle c
endif
else
if(abs(den_gvec(1,xyzzyaaaa32)-xyzzyaabz32(xyzzyaaab32))<xyzzyaabu32.a&
&nd.abs(den_gvec(2,xyzzyaaaa32)-xyzzyaaca32(xyzzyaaab32))<xyzzyaabu32.&
&and.abs(den_gvec(3,xyzzyaaaa32)-xyzzyaacb32(xyzzyaaab32))<xyzzyaabu32&
&)then
xyzzyaaav32=fourpi*xyzzyaaay32/xyzzyaaby32(xyzzyaaaa32)
mpc_hartree=mpc_hartree+xyzzyaaav32*(dble(den_sc(xyzzyaaaa32,xyzzyaaaj&
&32))**2+aimag(den_sc(xyzzyaaaa32,xyzzyaaaj32))**2)
cycle c
endif
endif
enddo
xyzzyaacb1(xyzzyaaaa32,xyzzyaaaj32)=0.d0
xyzzyaacc1(xyzzyaaaa32,xyzzyaaaj32)=0.d0
enddo c
endif
end select
enddo
mpc_correction=mpc_correction*real(npcells,dp)
if(hartree_xc)mpc_hartree=mpc_hartree*0.5d0
deallocate(xyzzyaaby32,xyzzyaabz32,xyzzyaaca32,xyzzyaacb32,xyzzyaacc32&
&)
return
6 call errstop('READ_MPC','Problem reading DENSITY DATA section of mpc&
&.data file.')
7 call errstop('READ_MPC','Problem reading EEPOT DATA section of mpc.d&
&ata file.')
end subroutine xyzzyaaed1
subroutine xyzzyaaee1
implicit none
integer xyzzyaaaa33
call mpi_bcast(xyzzyaaci1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting den_ngvec in setup_parallel_mpc.')
call mpi_bcast(xyzzyaack1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting den_ngvec_truncated in setup_parall&
&el_mpc.')
call mpi_bcast(xyzzyaaak1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ngtrunc in setup_parallel_mpc.')
call mpi_bcast(xyzzyaacj1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting ntypes in setup_parallel_mpc.')
call mpi_bcast(xyzzyaach1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting den_grange in setup_parallel_mpc.')
if(am_slave)then
allocate(xyzzyaacb1(xyzzyaaak1,xyzzyaacj1),xyzzyaacc1(xyzzyaaak1,xyzzy&
&aacj1),stat=xyzzyaaaa33)
call check_alloc(xyzzyaaaa33,'SETUP_PARALLEL_MPC','mpc_gspace')
allocate(xyzzyaaco1(3,-xyzzyaach1:xyzzyaach1),stat=xyzzyaaaa33)
call check_alloc(xyzzyaaaa33,'SETUP_PARALLEL_MPC','den_mwork')
allocate(xyzzyaacm1(xyzzyaaci1,1),xyzzyaacn1(xyzzyaaci1,1),stat=xyzzya&
&aaa33)
call check_alloc(xyzzyaaaa33,'SETUP_PARALLEL_MPC','den_cosgdotr')
xyzzyaacm1=0.d0
xyzzyaacn1=0.d0
allocate(xyzzyaacl1(3,xyzzyaaci1),stat=xyzzyaaaa33)
call check_alloc(xyzzyaaaa33,'SETUP_PARALLEL_MPC','den_pgmap')
endif
call mpi_bcast(xyzzyaacb1,xyzzyaaak1*xyzzyaacj1,mpi_double_precision,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting mpc_gspace_real in setup_parallel_m&
&pc.')
call mpi_bcast(xyzzyaacc1,xyzzyaaak1*xyzzyaacj1,mpi_double_precision,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting mpc_gspace_imag in setup_parallel_m&
&pc.')
call mpi_bcast(mpc_correction,1,mpi_double_precision,0,mpi_comm_world,&
&ierror)
call checkmpi(ierror,'Broadcasting mpc_correction in setup_parallel_mp&
&c.')
if(hartree_xc)then
call mpi_bcast(mpc_hartree,xyzzyaacj1,mpi_double_precision,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'Broadcasting mpc_hartree in setup_parallel_mpc.'&
&)
endif
call mpi_bcast(xyzzyaacl1,3*xyzzyaaci1,mpi_integer,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting den_pgmap in setup_parallel_mpc.')
end subroutine xyzzyaaee1
subroutine nuclear_repulsion_energy(n,z,nn)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: z(n)
real(dp),intent(out) :: nn
integer xyzzyaaaa34,xyzzyaaab34
real(dp) xyzzyaaac34
nn=0.d0
select case(periodicity)
case(3)
do xyzzyaaaa34=1,n
do xyzzyaaab34=xyzzyaaaa34+1,n
call ewald_3d(1,rionion(:,xyzzyaaab34,xyzzyaaaa34),xyzzyaaac34,.false.&
&)
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaab34)*xyzzyaaac34
enddo
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaaa34)*self_term
enddo
case(2)
do xyzzyaaaa34=1,n
do xyzzyaaab34=xyzzyaaaa34+1,n
call ewald_2d(1,rionion(:,xyzzyaaab34,xyzzyaaaa34),xyzzyaaac34)
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaab34)*xyzzyaaac34
enddo
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaaa34)*self_term
enddo
case(1)
do xyzzyaaaa34=1,n
do xyzzyaaab34=xyzzyaaaa34+1,n
call xyzzyaadi1(0,1,rionion(:,xyzzyaaab34,xyzzyaaaa34),xyzzyaaac34)
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaab34)*xyzzyaaac34
enddo
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaaa34)*self_term
enddo
case(0)
select case(iinteraction)
case(icoulomb,inone)
do xyzzyaaaa34=1,n-1
do xyzzyaaab34=xyzzyaaaa34+1,n
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaab34)/rionion(4,xyzzyaaab34,xyzzyaaaa34)
enddo
enddo
case(ilogarithmic)
do xyzzyaaaa34=1,n-1
do xyzzyaaab34=xyzzyaaaa34+1,n
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaab34)*(log(2.d0*(rstar+abs(rionion(3,xyz&
&zyaaab34,xyzzyaaaa34)))/rionion(4,xyzzyaaab34,xyzzyaaaa34)-euler)/(rs&
&tar+abs(rionion(3,xyzzyaaab34,xyzzyaaaa34))))
enddo
enddo
case(i2d_int)
do xyzzyaaaa34=1,n-1
do xyzzyaaab34=xyzzyaaaa34+1,n
nn=nn+z(xyzzyaaaa34)*z(xyzzyaaab34)*pot_2d_int(rionion(4,xyzzyaaab34,x&
&yzzyaaaa34),rionion(3,xyzzyaaab34,xyzzyaaaa34))
enddo
enddo
case default
call errstop('NUCLEAR_REPULSION_ENERGY','Currently only the logarithmi&
&c and 2D_int manual interactions are available for fixed charges.')
end select
case default
call errstop('NUCLEAR_REPULSION_ENERGY','Invalid periodicity ==> bug.'&
&)
end select
nn=nn/real(npcells,dp)
end subroutine nuclear_repulsion_energy
subroutine compute_fourier_basis_mpc(rvec,b)
implicit none
integer,intent(in) :: b
real(dp),intent(in) :: rvec(3)
integer xyzzyaaaa35
real(dp) xyzzyaaab35(3),xyzzyaaac35,xyzzyaaad35,xyzzyaaae35,xyzzyaaaf3&
&5,xyzzyaaag35,xyzzyaaah35
complex(dp) xyzzyaaai35
if(homogeneous_mpc)return
call timer('MPC_BASIS',.true.)
if(blip_mpc.and.electron_system.and.periodicity==3)then
xyzzyaadc1=rvec
else
xyzzyaaab35(1)=pb1(1)*rvec(1)+pb1(2)*rvec(2)+pb1(3)*rvec(3)
xyzzyaaab35(2)=pb2(1)*rvec(1)+pb2(2)*rvec(2)+pb2(3)*rvec(3)
xyzzyaaab35(3)=pb3(1)*rvec(1)+pb3(2)*rvec(2)+pb3(3)*rvec(3)
xyzzyaaac35=cos(xyzzyaaab35(1))
xyzzyaaad35=cos(xyzzyaaab35(2))
xyzzyaaae35=cos(xyzzyaaab35(3))
xyzzyaaaf35=sin(xyzzyaaab35(1))
xyzzyaaag35=sin(xyzzyaaab35(2))
xyzzyaaah35=sin(xyzzyaaab35(3))
xyzzyaaco1(1,0)=c_one
xyzzyaaco1(2,0)=c_one
xyzzyaaco1(3,0)=c_one
xyzzyaaco1(1,1)=cmplx(xyzzyaaac35,xyzzyaaaf35,dp)
xyzzyaaco1(2,1)=cmplx(xyzzyaaad35,xyzzyaaag35,dp)
xyzzyaaco1(3,1)=cmplx(xyzzyaaae35,xyzzyaaah35,dp)
xyzzyaaco1(1,-1)=cmplx(xyzzyaaac35,-xyzzyaaaf35,dp)
xyzzyaaco1(2,-1)=cmplx(xyzzyaaad35,-xyzzyaaag35,dp)
xyzzyaaco1(3,-1)=cmplx(xyzzyaaae35,-xyzzyaaah35,dp)
do xyzzyaaaa35=2,xyzzyaach1
xyzzyaaco1(1,xyzzyaaaa35)=xyzzyaaco1(1,1)*xyzzyaaco1(1,xyzzyaaaa35-1)
xyzzyaaco1(2,xyzzyaaaa35)=xyzzyaaco1(2,1)*xyzzyaaco1(2,xyzzyaaaa35-1)
xyzzyaaco1(3,xyzzyaaaa35)=xyzzyaaco1(3,1)*xyzzyaaco1(3,xyzzyaaaa35-1)
xyzzyaaco1(1,-xyzzyaaaa35)=conjg(xyzzyaaco1(1,xyzzyaaaa35))
xyzzyaaco1(2,-xyzzyaaaa35)=conjg(xyzzyaaco1(2,xyzzyaaaa35))
xyzzyaaco1(3,-xyzzyaaaa35)=conjg(xyzzyaaco1(3,xyzzyaaaa35))
enddo
do xyzzyaaaa35=1,xyzzyaaci1
xyzzyaaai35=xyzzyaaco1(1,xyzzyaacl1(1,xyzzyaaaa35))*xyzzyaaco1(2,xyzzy&
&aacl1(2,xyzzyaaaa35))*xyzzyaaco1(3,xyzzyaacl1(3,xyzzyaaaa35))
xyzzyaacm1(xyzzyaaaa35,b)=dble(xyzzyaaai35)
xyzzyaacn1(xyzzyaaaa35,b)=aimag(xyzzyaaai35)
enddo
endif
call timer('MPC_BASIS',.false.)
end subroutine compute_fourier_basis_mpc
subroutine xyzzyaaef1
use singleton, only : fftn
implicit none
integer xyzzyaaaa36,xyzzyaaab36,xyzzyaaac36,xyzzyaaad36,xyzzyaaae36(3)&
&,xyzzyaaaf36,xyzzyaaag36
real(dp) xyzzyaaah36(3),xyzzyaaai36
real(dp),allocatable :: xyzzyaaaj36(:,:)
complex(dp),allocatable :: xyzzyaaak36(:,:,:)
xyzzyaaaa36=3
do xyzzyaaab36=1,3
xyzzyaada1(xyzzyaaab36)=2*(maxval(abs(xyzzyaacl1(xyzzyaaab36,:)))*xyzz&
&yaaaa36+1)
enddo
xyzzyaadb1=real(xyzzyaada1,dp)
allocate(xyzzyaaak36(xyzzyaada1(1),xyzzyaada1(2),xyzzyaada1(3)),xyzzya&
&ade1(0:xyzzyaada1(1)-1,0:xyzzyaada1(2)-1,0:xyzzyaada1(3)-1),xyzzyaaaj&
&36(maxval(xyzzyaada1),3),stat=xyzzyaaag36)
call check_alloc(xyzzyaaag36,'CREATE_BLIP_FIT_MPC','grid_in')
xyzzyaaak36=0d0
do xyzzyaaaf36=1,xyzzyaack1
xyzzyaaae36=xyzzyaacl1(:,xyzzyaaaf36)+1
do xyzzyaaab36=1,3
if(xyzzyaaae36(xyzzyaaab36)<=0)xyzzyaaae36(xyzzyaaab36)=xyzzyaaae36(xy&
&zzyaaab36)+xyzzyaada1(xyzzyaaab36)
enddo
xyzzyaaak36(xyzzyaaae36(1),xyzzyaaae36(2),xyzzyaaae36(3))=cmplx(xyzzya&
&acb1(xyzzyaaaf36,1),xyzzyaacc1(xyzzyaaaf36,1),dp)
enddo
xyzzyaaah36=2d0*pi/xyzzyaadb1
do xyzzyaaab36=1,3
do xyzzyaaac36=-xyzzyaada1(xyzzyaaab36)/2+1,xyzzyaada1(xyzzyaaab36)/2
xyzzyaaad36=xyzzyaaac36+1
if(xyzzyaaad36<=0)xyzzyaaad36=xyzzyaaad36+xyzzyaada1(xyzzyaaab36)
xyzzyaaaj36(xyzzyaaad36,xyzzyaaab36)=1d0/(1d0+0.5d0*cos(xyzzyaaah36(xy&
&zzyaaab36)*real(xyzzyaaac36,dp)))
enddo
enddo
do xyzzyaaad36=1,xyzzyaada1(3)
do xyzzyaaac36=1,xyzzyaada1(2)
xyzzyaaai36=xyzzyaaaj36(xyzzyaaac36,2)*xyzzyaaaj36(xyzzyaaad36,3)
do xyzzyaaab36=1,xyzzyaada1(1)
xyzzyaaak36(xyzzyaaab36,xyzzyaaac36,xyzzyaaad36)=xyzzyaaak36(xyzzyaaab&
&36,xyzzyaaac36,xyzzyaaad36)*xyzzyaaaj36(xyzzyaaab36,1)*xyzzyaaai36
enddo
enddo
enddo
call fftn(xyzzyaaak36,xyzzyaada1,inv=.true.)
xyzzyaade1=real(xyzzyaaak36,dp)*sqrt(real(xyzzyaada1(1)*xyzzyaada1(2)*&
&xyzzyaada1(3),dp))
deallocate(xyzzyaaak36,xyzzyaaaj36)
end subroutine xyzzyaaef1
real(dp) function xyzzyaaeg1(rvec)
implicit none
real(dp),intent(in) :: rvec(3)
integer xyzzyaaaa37,xyzzyaaab37,xyzzyaaac37,xyzzyaaad37(4,3)
real(dp) xyzzyaaae37(3),xyzzyaaaf37(3),xyzzyaaag37(4),xyzzyaaah37(4,3)&
&,xyzzyaaai37
xyzzyaaae37=matmul(rvec(1:3),painv)
xyzzyaaad37(2,1)=modulo(floor(xyzzyaaae37(1)*xyzzyaadb1(1)),xyzzyaada1&
&(1))
xyzzyaaad37(2,2)=modulo(floor(xyzzyaaae37(2)*xyzzyaadb1(2)),xyzzyaada1&
&(2))
xyzzyaaad37(2,3)=modulo(floor(xyzzyaaae37(3)*xyzzyaadb1(3)),xyzzyaada1&
&(3))
do xyzzyaaaa37=1,3
xyzzyaaad37(1,xyzzyaaaa37)=modulo(xyzzyaaad37(2,xyzzyaaaa37)-1,xyzzyaa&
&da1(xyzzyaaaa37))
xyzzyaaad37(3,xyzzyaaaa37)=modulo(xyzzyaaad37(2,xyzzyaaaa37)+1,xyzzyaa&
&da1(xyzzyaaaa37))
xyzzyaaad37(4,xyzzyaaaa37)=modulo(xyzzyaaad37(2,xyzzyaaaa37)+2,xyzzyaa&
&da1(xyzzyaaaa37))
enddo
xyzzyaaaf37=modulo(xyzzyaaae37*xyzzyaadb1,xyzzyaadb1)
do xyzzyaaaa37=1,3
xyzzyaaag37(1)=xyzzyaaaf37(xyzzyaaaa37)-dble(xyzzyaaad37(2,xyzzyaaaa37&
&)-1)
xyzzyaaag37(2)=xyzzyaaaf37(xyzzyaaaa37)-dble(xyzzyaaad37(2,xyzzyaaaa37&
&))
xyzzyaaag37(3)=xyzzyaaaf37(xyzzyaaaa37)-dble(xyzzyaaad37(2,xyzzyaaaa37&
&)+1)
xyzzyaaag37(4)=xyzzyaaaf37(xyzzyaaaa37)-dble(xyzzyaaad37(2,xyzzyaaaa37&
&)+2)
xyzzyaaah37(1,xyzzyaaaa37)=2d0+xyzzyaaag37(1)*(-3d0+xyzzyaaag37(1)*(1.&
&5d0-0.25d0*xyzzyaaag37(1)))
xyzzyaaah37(2,xyzzyaaaa37)=1d0+xyzzyaaag37(2)*xyzzyaaag37(2)*(-1.5d0+0&
&.75d0*xyzzyaaag37(2))
xyzzyaaah37(3,xyzzyaaaa37)=1d0+xyzzyaaag37(3)*xyzzyaaag37(3)*(-1.5d0-0&
&.75d0*xyzzyaaag37(3))
xyzzyaaah37(4,xyzzyaaaa37)=2d0+xyzzyaaag37(4)*(3d0+xyzzyaaag37(4)*(1.5&
&d0+0.25d0*xyzzyaaag37(4)))
enddo
xyzzyaaeg1=0d0
do xyzzyaaac37=1,4
do xyzzyaaab37=1,4
xyzzyaaai37=xyzzyaaah37(xyzzyaaab37,2)*xyzzyaaah37(xyzzyaaac37,3)
do xyzzyaaaa37=1,4
xyzzyaaeg1=xyzzyaaeg1+xyzzyaade1(xyzzyaaad37(xyzzyaaaa37,1),xyzzyaaad3&
&7(xyzzyaaab37,2),xyzzyaaad37(xyzzyaaac37,3))*xyzzyaaah37(xyzzyaaaa37,&
&1)*xyzzyaaai37
enddo
enddo
enddo
end function xyzzyaaeg1
subroutine eval_fscorr_consts(c)
use slaarnaag,only : twothirds,fourthirds
use slaarnabg,only : binv,b11,b22,b33
use slaarnabt,only : quicksort,ddot
implicit none
real(dp),intent(out) :: c
integer xyzzyaaaa38,xyzzyaaab38,xyzzyaaac38,xyzzyaaad38,xyzzyaaae38,xy&
&zzyaaaf38,xyzzyaaag38,xyzzyaaah38,xyzzyaaai38,xyzzyaaaj38, xyzzyaaak3&
&8,xyzzyaaal38,xyzzyaaam38,xyzzyaaan38
integer,allocatable :: xyzzyaaao38(:),xyzzyaaap38(:)
integer,parameter :: xyzzyaaaq38=50000000
integer,parameter :: xyzzyaaar38=3
integer,parameter :: xyzzyaaas38=2
real(dp) alpha,xyzzyaaat38,xyzzyaaau38(3),xyzzyaaav38(3),xyzzyaaaw38(3&
&),xyzzyaaax38(3),xyzzyaaay38(3),xyzzyaaaz38, xyzzyaaba38,xyzzyaabb38,&
&xyzzyaabc38,xyzzyaabd38,xyzzyaabe38,xyzzyaabf38(xyzzyaaar38), xyzzyaa&
&bg38(xyzzyaaar38),xyzzyaabh38(xyzzyaaar38),xyzzyaabi38
real(dp),allocatable :: xyzzyaabj38(:),xyzzyaabk38(:),xyzzyaabl38(:),x&
&yzzyaabm38(:), xyzzyaabn38(:)
real(dp),parameter :: xyzzyaabo38=1.d-10
logical,parameter :: xyzzyaabp38=.true.,xyzzyaabq38=.false.
character(80) char80
if(.not.am_master)return
call timer('EVAL_FSCORR_CONSTS',.true.)
if(periodicity/=2.and.periodicity/=3)call errstop('EVAL_FSCORR_CONSTS'&
&, 'Finite-size corrections are only implemented for 3D- or 2D-periodi&
&c systems.')
if(xyzzyaaas38>1)then
call wout('Evaluation of lattice-dependent constants for finite-size c&
&orrection')
call wout('===========================================================&
&=========')
endif
call timer('FSCORR_RAW_GRID_SETUP',.true.)
if(xyzzyaaas38>2)call wout('Target number of G vectors in sum    : ' /&
&/trim(i2s(xyzzyaaaq38)))
if(periodicity==3)then
xyzzyaaaz38=(6.d0*pi**2*dble(xyzzyaaaq38)/volume)**third
else
xyzzyaaaz38=sqrt(fourpi*dble(xyzzyaaaq38)/area)
endif
if(xyzzyaaas38>2)call wout('Cutoff radius in reciprocal space    : ',x&
&yzzyaaaz38)
xyzzyaaba38=xyzzyaaaz38**2
xyzzyaaae38=int(xyzzyaaaz38*sqrt(dot_product(binv(1:3,1),binv(1:3,1)))&
&)
xyzzyaaaf38=int(xyzzyaaaz38*sqrt(dot_product(binv(1:3,2),binv(1:3,2)))&
&)
if(periodicity==3)then
xyzzyaaag38=int(xyzzyaaaz38*sqrt(dot_product(binv(1:3,3),binv(1:3,3)))&
&)
else
xyzzyaaag38=0
endif
if(periodicity==3)then
xyzzyaabi38=min(b11,b22,b33)
else
xyzzyaabi38=min(b11,b22)
endif
xyzzyaabb38=xyzzyaabo38*xyzzyaabi38
xyzzyaaad38=xyzzyaaae38+1+xyzzyaaaf38*(2*xyzzyaaae38+1)+xyzzyaaag38*(2&
&*xyzzyaaaf38+1) *(2*xyzzyaaae38+1)
if(xyzzyaaas38>2)call wout('Number of G vectors in full grid     : ' /&
&/trim(i2s(2*xyzzyaaad38-1)))
allocate(xyzzyaabj38(xyzzyaaad38),stat=xyzzyaaaa38)
call check_alloc(xyzzyaaaa38,'EVAL_FSCORR_CONSTS','gsq_grid')
xyzzyaaax38=-dble(xyzzyaaae38+1)*b1
xyzzyaaay38=-dble(xyzzyaaaf38+1)*b2
xyzzyaaah38=1
xyzzyaabj38(1)=0.d0
xyzzyaaau38=0.d0
do xyzzyaaai38=1,xyzzyaaae38
xyzzyaaau38=xyzzyaaau38+b1
xyzzyaabe38=sum(xyzzyaaau38(1:periodicity)**2)
if(xyzzyaabe38<xyzzyaaba38)then
xyzzyaaah38=xyzzyaaah38+1
xyzzyaabj38(xyzzyaaah38)=xyzzyaabe38
endif
enddo
xyzzyaaaw38=0.d0
do xyzzyaaaj38=1,xyzzyaaaf38
xyzzyaaaw38=xyzzyaaaw38+b2
xyzzyaaau38=xyzzyaaaw38+xyzzyaaax38
do xyzzyaaai38=-xyzzyaaae38,xyzzyaaae38
xyzzyaaau38=xyzzyaaau38+b1
xyzzyaabe38=sum(xyzzyaaau38(1:periodicity)**2)
if(xyzzyaabe38<xyzzyaaba38)then
xyzzyaaah38=xyzzyaaah38+1
xyzzyaabj38(xyzzyaaah38)=xyzzyaabe38
endif
enddo
enddo
xyzzyaaav38=0.d0
do xyzzyaaak38=1,xyzzyaaag38
xyzzyaaav38=xyzzyaaav38+b3
xyzzyaaaw38=xyzzyaaav38+xyzzyaaay38
do xyzzyaaaj38=-xyzzyaaaf38,xyzzyaaaf38
xyzzyaaaw38=xyzzyaaaw38+b2
xyzzyaaau38=xyzzyaaaw38+xyzzyaaax38
do xyzzyaaai38=-xyzzyaaae38,xyzzyaaae38
xyzzyaaau38=xyzzyaaau38+b1
xyzzyaabe38=sum(xyzzyaaau38(1:periodicity)**2)
if(xyzzyaabe38<xyzzyaaba38)then
xyzzyaaah38=xyzzyaaah38+1
xyzzyaabj38(xyzzyaaah38)=xyzzyaabe38
endif
enddo
enddo
enddo
xyzzyaaad38=xyzzyaaah38
call timer('FSCORR_RAW_GRID_SETUP',.false.)
call timer('FSCORR_SORT',.true.)
allocate(xyzzyaaao38(xyzzyaaad38),stat=xyzzyaaaa38)
call check_alloc(xyzzyaaaa38,'EVAL_FSCORR_CONSTS','ira')
call quicksort(xyzzyaaad38,xyzzyaabj38(1),xyzzyaaao38(1))
call timer('FSCORR_SORT',.false.)
call timer('FSCORR_STARS',.true.)
allocate(xyzzyaabk38(xyzzyaaad38),xyzzyaaap38(xyzzyaaad38),stat=xyzzya&
&aaa38)
call check_alloc(xyzzyaaaa38,'EVAL_FSCORR_CONSTS','gsq_star')
xyzzyaaac38=1
xyzzyaabk38(1)=0.d0
xyzzyaaap38(1)=1
do xyzzyaaah38=2,xyzzyaaad38
if(xyzzyaabj38(xyzzyaaao38(xyzzyaaah38))>xyzzyaabk38(xyzzyaaac38)+xyzz&
&yaabb38)then
xyzzyaaac38=xyzzyaaac38+1
xyzzyaabk38(xyzzyaaac38)=xyzzyaabj38(xyzzyaaao38(xyzzyaaah38))
xyzzyaaap38(xyzzyaaac38)=2
else
xyzzyaaap38(xyzzyaaac38)=xyzzyaaap38(xyzzyaaac38)+2
endif
enddo
deallocate(xyzzyaabj38,xyzzyaaao38)
xyzzyaaam38=xyzzyaaac38-1
if(xyzzyaaas38>2)call wout('Number of stars in sum               : ' /&
&/trim(i2s(xyzzyaaac38)))
xyzzyaaal38=sum(xyzzyaaap38(1:xyzzyaaac38))
if(xyzzyaaas38>2)call wout('Actual number of G vectors in sum    : ' /&
&/trim(i2s(xyzzyaaal38)))
if(periodicity==3)then
xyzzyaaaz38=(6.d0*pi**2*dble(xyzzyaaal38)/volume)**third
else
xyzzyaaaz38=sqrt(fourpi*dble(xyzzyaaal38)/area)
endif
if(xyzzyaaas38>2)call wout('Effective cutoff radius              : ',x&
&yzzyaaaz38)
call timer('FSCORR_STARS',.false.)
call timer('FSCORR_BUFFERS',.true.)
if(xyzzyaabp38)then
allocate(xyzzyaabm38(2:xyzzyaaac38),stat=xyzzyaaaa38)
call check_alloc(xyzzyaaaa38,'EVAL_FSCORR_CONSTS','ninstar_g')
endif
if(xyzzyaabq38)then
if(periodicity/=3)call errstop('EVAL_FSCORR_CONSTS', 'Please work out &
&the expression for the HF finite-size correction in 2D.')
allocate(xyzzyaabn38(2:xyzzyaaac38),stat=xyzzyaaaa38)
call check_alloc(xyzzyaaaa38,'EVAL_FSCORR_CONSTS','ninstar_rec_g')
endif
if(xyzzyaabp38)then
if(xyzzyaabq38)then
do xyzzyaaab38=2,xyzzyaaac38
xyzzyaabn38(xyzzyaaab38)=dble(xyzzyaaap38(xyzzyaaab38))/sqrt(xyzzyaabk&
&38(xyzzyaaab38))
xyzzyaabm38(xyzzyaaab38)=xyzzyaabn38(xyzzyaaab38)*xyzzyaabk38(xyzzyaaa&
&b38)
enddo
else
if(periodicity==3)then
do xyzzyaaab38=2,xyzzyaaac38
xyzzyaabm38(xyzzyaaab38)=dble(xyzzyaaap38(xyzzyaaab38))*sqrt(xyzzyaabk&
&38(xyzzyaaab38))
enddo
else
do xyzzyaaab38=2,xyzzyaaac38
xyzzyaabm38(xyzzyaaab38)=dble(xyzzyaaap38(xyzzyaaab38))*xyzzyaabk38(xy&
&zzyaaab38)**0.25d0
enddo
endif
endif
elseif(xyzzyaabq38)then
do xyzzyaaab38=2,xyzzyaaac38
xyzzyaabn38(xyzzyaaab38)=dble(xyzzyaaap38(xyzzyaaab38))/sqrt(xyzzyaabk&
&38(xyzzyaaab38))
enddo
else
call errstop('EVAL_FSCORR_CONSTS','Bug.')
endif
deallocate(xyzzyaaap38)
call timer('FSCORR_BUFFERS',.false.)
allocate(xyzzyaabl38(2:xyzzyaaac38),stat=xyzzyaaaa38)
call check_alloc(xyzzyaaaa38,'EVAL_FSCORR_CONSTS','exp_magsq')
call timer('FSCORR_CALCULATIONS',.true.)
if(xyzzyaaar38==3)then
xyzzyaabf38(3)=250.d0/xyzzyaaba38
xyzzyaabf38(2)=4.d0*xyzzyaabf38(3)
xyzzyaabf38(1)=25.d0*xyzzyaabf38(3)
if(xyzzyaabf38(1)>0.5d0/xyzzyaabi38)call errstop('EVAL_FSCORR_CONSTS',&
&'Need more G vectors.')
endif
do xyzzyaaan38=1,xyzzyaaar38
alpha=xyzzyaabf38(xyzzyaaan38)
do xyzzyaaab38=2,xyzzyaaac38
xyzzyaabl38(xyzzyaaab38)=exp(-alpha*xyzzyaabk38(xyzzyaaab38))
enddo
if(xyzzyaabp38)then
if(periodicity==3)then
xyzzyaabc38=volume*exp(-alpha*xyzzyaaba38)*(1.d0+alpha*xyzzyaaba38)/(t&
&wopi*alpha)**2
else
xyzzyaabc38=0.d0
endif
xyzzyaabd38=ddot(xyzzyaaam38,xyzzyaabl38(2),1,xyzzyaabm38(2),1)+xyzzya&
&abc38
if(periodicity==3)then
xyzzyaabg38(xyzzyaaan38)=volume**fourthirds*(0.25d0/(pi*alpha**2) -pi/&
&volume*xyzzyaabd38)
else
xyzzyaabg38(xyzzyaaan38)=area**1.25d0*(0.226600619263869d0/alpha**1.25&
&d0 -pi/area*xyzzyaabd38)
endif
endif
if(xyzzyaabq38)then
xyzzyaabc38=volume*exp(-alpha*xyzzyaaba38)/(twopi**2*alpha)
xyzzyaabd38=ddot(xyzzyaaam38,xyzzyaabl38(2),1,xyzzyaabn38(2),1) +xyzzy&
&aabc38
xyzzyaabh38(xyzzyaaan38)=volume**twothirds*(one_over_twopi/alpha -twop&
&i/volume*xyzzyaabd38)
endif
enddo
deallocate(xyzzyaabk38,xyzzyaabl38)
if(xyzzyaabp38)deallocate(xyzzyaabm38)
if(xyzzyaabq38)deallocate(xyzzyaabn38)
if(xyzzyaabp38)then
if(xyzzyaaas38>1)then
call wordwrap('Please check that the following quadratic extrapolation&
& of C to alpha=0 looks reasonable.  [See Eqs. (56) and (59) of Phys. &
&Rev. B 78, 125106 (2008) for the definition of C in 3D and 2D, resect&
&ively.]')
call wout('===================================================')
call wout('         alpha                     C')
call wout('---------------------------------------------------')
do xyzzyaaan38=1,xyzzyaaar38
write(char80,'(1x,es24.16,2x,es24.16)')xyzzyaabf38(xyzzyaaan38),xyzzya&
&abg38(xyzzyaaan38)
call wout(trim(char80))
enddo
endif
if(xyzzyaaar38==3)then
if(xyzzyaaas38>1)call wout('------------------------------------------&
&---------')
c=xyzzyaabg38(1)*xyzzyaabf38(2)*xyzzyaabf38(3) /((xyzzyaabf38(1)-xyzzy&
&aabf38(2))*(xyzzyaabf38(1)-xyzzyaabf38(3))) +xyzzyaabg38(2)*xyzzyaabf&
&38(1)*xyzzyaabf38(3) /((xyzzyaabf38(2)-xyzzyaabf38(1))*(xyzzyaabf38(2&
&)-xyzzyaabf38(3))) +xyzzyaabg38(3)*xyzzyaabf38(1)*xyzzyaabf38(2) /((x&
&yzzyaabf38(3)-xyzzyaabf38(1))*(xyzzyaabf38(3)-xyzzyaabf38(2)))
if(xyzzyaaas38>1)then
write(char80,'(1x,es24.16,2x,es24.16)')0.d0,c
call wout(trim(char80))
endif
endif
if(xyzzyaaas38>1)then
call wout('===================================================')
call wout
endif
endif
if(xyzzyaabq38)then
if(xyzzyaaas38>1)then
call wordwrap('Please check that the following quadratic extrapolation&
& of C_HF to alpha=0 looks reasonable.  [See Eq. (40) of Phys. Rev. B &
&78, 125106 (2008) for the definition of C_HF.]')
call wout('===================================================')
call wout('         alpha                     C_HF')
call wout('---------------------------------------------------')
do xyzzyaaan38=1,xyzzyaaar38
write(char80,'(1x,es24.16,2x,es24.16)')xyzzyaabf38(xyzzyaaan38),xyzzya&
&abh38(xyzzyaaan38)
call wout(trim(char80))
enddo
endif
if(xyzzyaaar38==3)then
if(xyzzyaaas38>1)call wout('------------------------------------------&
&---------')
xyzzyaaat38=xyzzyaabh38(1)*xyzzyaabf38(2)*xyzzyaabf38(3) /((xyzzyaabf3&
&8(1)-xyzzyaabf38(2))*(xyzzyaabf38(1)-xyzzyaabf38(3))) +xyzzyaabh38(2)&
&*xyzzyaabf38(1)*xyzzyaabf38(3) /((xyzzyaabf38(2)-xyzzyaabf38(1))*(xyz&
&zyaabf38(2)-xyzzyaabf38(3))) +xyzzyaabh38(3)*xyzzyaabf38(1)*xyzzyaabf&
&38(2) /((xyzzyaabf38(3)-xyzzyaabf38(1))*(xyzzyaabf38(3)-xyzzyaabf38(2&
&)))
if(xyzzyaaas38>1)then
write(char80,'(1x,es24.16,2x,es24.16)')0.d0,xyzzyaaat38
call wout(trim(char80))
endif
endif
if(xyzzyaaas38>1)then
call wout('===================================================')
call wout
endif
endif
call timer('FSCORR_CALCULATIONS',.false.)
call timer('EVAL_FSCORR_CONSTS',.false.)
end subroutine eval_fscorr_consts
end module slaarnabk

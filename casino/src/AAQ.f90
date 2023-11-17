module slaarnaaq
use dsp
use slaarnabg
use parallel
use slaarnach
use store
use slaarnacq
use slaarnacs
use slaarnaag,     only : czero,one_over_twopi,twopi,fourpi,one_over_p&
&i,third,pi
use slaarnaan, only : mol_system_size
use file_utils,    only : open_units
use format_utils,  only : wout,i2s,r2s,r2s2,wordwrap
use slaarnabt,     only : approx_equal,quicksort,gauss_quadrature,dsca&
&l,daxpy,bessel_j0,bessel_j1,ddot,exp_protect,sph_bessel,chebyshev,exp&
&onential
use slaarnacc,only : ranx
use run_control,   only : errstop,errstop2,errstop_master,errwarn,time&
&r,check_alloc
implicit none
private
public read_expval,setup_expval,write_expval,eval_int_sf,setup_emt,dip&
&ole_calc,eval_residuals_xc_corr,eval_jacobian_xc_corr,accumulate_expv&
&als,finite_size_corr_xc,expval_zero_postweight,expval_postweight,swit&
&ch_spin_density_to_pcf,switch_spin_density_to_pcf_sph,write_density,e&
&xpval_scratch_request,setup_expval_accum,finish_expval_accum,contact_&
&den_calc,backup_expval_file,deallocate_expval,zero_expval,expval_pcf
public expvals,density,spin_density,spin_density_mat,pair_corr,pair_co&
&rr_sph,loc_tensor,structure_factor,structure_factor_sph,onep_density_&
&mat,twop_density_mat,cond_fraction,mom_den,int_sf,eval_dipole_moment,&
&eval_contact_den,finite_density,population,mol_density,mol_spin_densi&
&ty,twop_dm_mom,cond_frac_mom,pcfs_nbins_in,pcfs_rcutoff_in
public expval_error_bars,pcf_rfix,pcfs_rfix,particle_is_fixed,input_rf&
&ix,input_type_fixed,fixed_particle,expval_ka,expval_kb,expval_kc,expv&
&al_kd,expval_nk,expval_cutoff,expval_nkdim,expval_nkgrids,expval_kgri&
&d_id,nexpvals_need_kgrid,xc_corr_method,turn_on_pair_corr,turn_on_pai&
&r_corr_sph
public den_nsets,permit_den_symm,den_sc,tol_expval_diff,expval_den_mpc&
&,qmc_density_mpc,expval_ngvec,expval_ngvec_truncated,expval_gvec,den_&
&gvec,homogeneous_density,expval_complex_den
public a,nparam1,nparam2
integer xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1,xyzzyaaad1,xyzzyaaae1,xyzzyaa&
&af1,xyzzyaaag1,xyzzyaaah1,xyzzyaaai1
integer,parameter :: xyzzyaaaj1=1,xyzzyaaak1=10,xyzzyaaal1=17,nexpvals&
&_need_kgrid=3
real(dp) expval_cutoff,xyzzyaaam1,xyzzyaaan1,xyzzyaaao1,xyzzyaaap1(3),&
&xyzzyaaaq1(3),xyzzyaaar1(3)
real(dp) :: tol_expval_diff=1.d-6
logical expvals,xyzzyaaas1,xyzzyaaat1,xyzzyaaau1,xyzzyaaav1,xyzzyaaaw1&
&,xyzzyaaax1(xyzzyaaal1),int_sf,xyzzyaaay1,xyzzyaaaz1,xyzzyaaba1,xyzzy&
&aabb1,qmc_density_mpc,density,xyzzyaabc1,spin_density,xyzzyaabd1,spin&
&_density_mat,xyzzyaabe1,loc_tensor,xyzzyaabf1,xyzzyaabg1,structure_fa&
&ctor,xyzzyaabh1,xyzzyaabi1,structure_factor_sph,xyzzyaabj1,xyzzyaabk1&
&,pair_corr,xyzzyaabl1,xyzzyaabm1,pair_corr_sph,xyzzyaabn1,xyzzyaabo1,&
&onep_density_mat,xyzzyaabp1,xyzzyaabq1,twop_density_mat,xyzzyaabr1,xy&
&zzyaabs1,cond_fraction,xyzzyaabt1,xyzzyaabu1,expval_error_bars,mom_de&
&n,xyzzyaabv1,xyzzyaabw1,finite_density,xyzzyaabx1,xyzzyaaby1,turn_on_&
&pair_corr,turn_on_pair_corr_sph,population,xyzzyaabz1,xyzzyaaca1,mol_&
&density,xyzzyaacb1,xyzzyaacc1,twop_dm_mom,xyzzyaacd1,xyzzyaace1,cond_&
&frac_mom,xyzzyaacf1,xyzzyaacg1
character(80),allocatable :: xyzzyaach1(:)
integer xyzzyaaci1,xyzzyaacj1(3,3),xyzzyaack1,xyzzyaacl1,xyzzyaacm1
integer,allocatable :: xyzzyaacn1(:),xyzzyaaco1(:)
real(dp) xyzzyaacp1,xyzzyaacq1,xyzzyaacr1
logical permit_den_symm
character(80) title
integer xyzzyaacs1(3),xyzzyaact1(3)
integer xyzzyaacu1,xyzzyaacv1,xyzzyaacw1,xyzzyaacx1(9)
integer,allocatable :: expval_ngvec(:),expval_ngvec_truncated(:),xyzzy&
&aacy1(:,:,:),xyzzyaacz1(:)
real(dp) xyzzyaada1(9)
real(dp),allocatable :: xyzzyaadb1(:),expval_gvec(:,:,:),xyzzyaadc1(:,&
&:,:),xyzzyaadd1(:,:,:)
complex(dp),allocatable :: xyzzyaade1(:,:,:),xyzzyaadf1(:,:,:)
integer expval_nkgrids,expval_nkdim(nexpvals_need_kgrid),expval_kgrid_&
&id(nexpvals_need_kgrid),expval_nk(3,nexpvals_need_kgrid)
real(dp) expval_ka(3,nexpvals_need_kgrid),expval_kb(3,nexpvals_need_kg&
&rid),expval_kc(3,nexpvals_need_kgrid),expval_kd(3,nexpvals_need_kgrid&
&)
real(dp),allocatable :: den_gvec(:,:)
complex(dp),allocatable :: den_sc(:,:)
integer xyzzyaadg1,den_nsets
integer,allocatable :: xyzzyaadh1(:)
real(dp),allocatable :: xyzzyaadi1(:),xyzzyaadj1(:)
complex(dp),allocatable :: xyzzyaadk1(:,:),xyzzyaadl1(:,:),xyzzyaadm1(&
&:,:),expval_den_mpc(:,:),xyzzyaadn1(:,:)
logical expval_complex_den
character(3) accumulation_method_den
integer xyzzyaado1,xyzzyaadp1
integer,allocatable :: xyzzyaadq1(:)
real(dp),allocatable :: xyzzyaadr1(:),xyzzyaads1(:)
complex(dp),allocatable :: xyzzyaadt1(:,:),xyzzyaadu1(:,:),xyzzyaadv1(&
&:,:),xyzzyaadw1(:,:)
logical xyzzyaadx1
character(3) accumulation_method_sden
integer xyzzyaady1,xyzzyaadz1
integer,allocatable :: xyzzyaaea1(:)
real(dp),allocatable :: xyzzyaaeb1(:),xyzzyaaec1(:)
complex(dp),allocatable :: xyzzyaaed1(:,:,:),xyzzyaaee1(:,:,:),xyzzyaa&
&ef1(:,:,:)
logical xyzzyaaeg1
character(3) accumulation_method_sdenm
integer input_type_fixed,fixed_particle
real(dp) input_rfix(3)
logical particle_is_fixed
integer xyzzyaaeh1,xyzzyaaei1,xyzzyaaej1
integer,allocatable :: xyzzyaaek1(:)
real(dp) :: pcf_rfix(3)
real(dp),allocatable :: xyzzyaael1(:),xyzzyaaem1(:)
complex(dp),allocatable :: expval_pcf(:,:),xyzzyaaen1(:,:),xyzzyaaeo1(&
&:,:),xyzzyaaep1(:,:)
logical xyzzyaaeq1
character(3) accumulation_method_pcf
integer xyzzyaaer1,xyzzyaaes1,xyzzyaaet1,xyzzyaaeu1,pcfs_nbins_in
integer,allocatable :: xyzzyaaev1(:,:)
real(dp) xyzzyaaew1,pcfs_rfix(3),xyzzyaaex1,xyzzyaaey1,xyzzyaaez1,xyzz&
&yaafa1,pcfs_rcutoff_in
real(dp),allocatable :: xyzzyaafb1(:),xyzzyaafc1(:,:,:),xyzzyaafd1(:,:&
&,:),xyzzyaafe1(:,:,:),xyzzyaaff1(:,:,:)
character(3) accumulation_method_pcfs
logical homogeneous_density
integer xyzzyaafg1
real(dp), allocatable :: xyzzyaafh1(:),xyzzyaafi1(:),xyzzyaafj1(:),xyz&
&zyaafk1(:)
complex(dp),allocatable :: xyzzyaafl1(:,:,:),xyzzyaafm1(:,:,:),xyzzyaa&
&fn1(:,:,:),xyzzyaafo1(:,:,:)
logical xyzzyaafp1
character(3) accumulation_method_lt
integer xyzzyaafq1,xyzzyaafr1,xyzzyaafs1
integer,allocatable :: xyzzyaaft1(:,:),xyzzyaafu1(:)
real(dp) xyzzyaafv1,xyzzyaafw1,xyzzyaafx1
real(dp),allocatable :: xyzzyaafy1(:),xyzzyaafz1(:,:),xyzzyaaga1(:,:),&
&xyzzyaagb1(:,:),xyzzyaagc1(:,:),xyzzyaagd1(:),xyzzyaage1(:)
complex(dp),allocatable :: xyzzyaagf1(:,:),xyzzyaagg1(:,:),xyzzyaagh1(&
&:,:),xyzzyaagi1(:,:),xyzzyaagj1(:,:)
logical xyzzyaagk1
character(3) accumulation_method_sf
integer xyzzyaagl1,xyzzyaagm1
integer,allocatable :: xyzzyaagn1(:,:)
real(dp) xyzzyaago1,xyzzyaagp1
real(dp),allocatable :: xyzzyaagq1(:),xyzzyaagr1(:),xyzzyaags1(:),xyzz&
&yaagt1(:),xyzzyaagu1(:,:),xyzzyaagv1(:,:),xyzzyaagw1(:,:),xyzzyaagx1(&
&:,:),xyzzyaagy1(:)
character(3) accumulation_method_sf_sph
integer xyzzyaagz1,xyzzyaaha1,xyzzyaahb1
integer,allocatable :: xyzzyaahc1(:),xyzzyaahd1(:,:)
real(dp) xyzzyaahe1
real(dp),allocatable :: xyzzyaahf1(:,:),xyzzyaahg1(:,:),xyzzyaahh1(:,:&
&),xyzzyaahi1(:,:),xyzzyaahj1(:,:),xyzzyaahk1(:,:),xyzzyaahl1(:,:),xyz&
&zyaahm1(:,:),xyzzyaahn1(:,:),xyzzyaaho1(:),xyzzyaahp1(:,:),xyzzyaahq1&
&(:,:),xyzzyaahr1(:,:)
character(3) accumulation_method_onep_dm
integer xyzzyaahs1,xyzzyaaht1
integer,allocatable :: xyzzyaahu1(:),xyzzyaahv1(:,:)
real(dp),allocatable :: xyzzyaahw1(:),xyzzyaahx1(:),xyzzyaahy1(:),xyzz&
&yaahz1(:,:),xyzzyaaia1(:,:),xyzzyaaib1(:,:),xyzzyaaic1(:,:),xyzzyaaid&
&1(:,:),xyzzyaaie1(:,:),xyzzyaaif1(:),xyzzyaaig1(:,:),xyzzyaaih1(:,:),&
&xyzzyaaii1(:),xyzzyaaij1(:,:)
character(3) accumulation_method_mom_den
integer xyzzyaaik1,xyzzyaail1,xyzzyaaim1
integer,allocatable :: xyzzyaain1(:),xyzzyaaio1(:,:,:),xyzzyaaip1(:),x&
&yzzyaaiq1(:,:,:),xyzzyaair1(:,:),xyzzyaais1(:),xyzzyaait1(:)
real(dp) xyzzyaaiu1,xyzzyaaiv1
real(dp),allocatable :: xyzzyaaiw1(:,:),xyzzyaaix1(:,:),xyzzyaaiy1(:,:&
&),xyzzyaaiz1(:,:),xyzzyaaja1(:,:),xyzzyaajb1(:,:),xyzzyaajc1(:,:),xyz&
&zyaajd1(:,:),xyzzyaaje1(:,:),xyzzyaajf1(:,:),xyzzyaajg1(:,:),xyzzyaaj&
&h1(:,:),xyzzyaaji1(:)
character(3) accumulation_method_twop_dm
integer xyzzyaajj1,xyzzyaajk1,xyzzyaajl1
integer,allocatable :: xyzzyaajm1(:),xyzzyaajn1(:,:,:),xyzzyaajo1(:),x&
&yzzyaajp1(:,:,:),xyzzyaajq1(:,:),xyzzyaajr1(:),xyzzyaajs1(:)
real(dp) xyzzyaajt1,xyzzyaaju1
real(dp),allocatable :: xyzzyaajv1(:,:),xyzzyaajw1(:,:),xyzzyaajx1(:,:&
&),xyzzyaajy1(:,:),xyzzyaajz1(:,:),xyzzyaaka1(:,:),xyzzyaakb1(:,:),xyz&
&zyaakc1(:,:),xyzzyaakd1(:,:),xyzzyaake1(:,:),xyzzyaakf1(:,:),xyzzyaak&
&g1(:,:),xyzzyaakh1(:)
character(3) accum_method_cond_frac
integer xyzzyaaki1,xyzzyaakj1,xyzzyaakk1,xyzzyaakl1
integer,parameter :: xyzzyaakm1=1,xyzzyaakn1=2,xyzzyaako1=3,xyzzyaakp1&
&=4
real(dp) xyzzyaakq1,xyzzyaakr1,xyzzyaaks1,xyzzyaakt1,    xyzzyaaku1,xy&
&zzyaakv1,xyzzyaakw1,xyzzyaakx1,    xyzzyaaky1,xyzzyaakz1,xyzzyaala1,x&
&yzzyaalb1,  xyzzyaalc1,xyzzyaald1
real(dp),allocatable :: xyzzyaale1(:,:),xyzzyaalf1(:,:),xyzzyaalg1(:,:&
&),xyzzyaalh1(:,:),xyzzyaali1(:,:),xyzzyaalj1(:,:),xyzzyaalk1(:,:),xyz&
&zyaall1(:,:),xyzzyaalm1(:,:),xyzzyaaln1(:,:),xyzzyaalo1(:,:),xyzzyaal&
&p1(:,:),xyzzyaalq1(:,:),xyzzyaalr1(:,:),xyzzyaals1(:,:),xyzzyaalt1(:,&
&:),xyzzyaalu1(:)
character(3) accumulation_method_fin_den
character(20) fin_den_basis
integer xyzzyaalv1,xyzzyaalw1(3)
real(dp) xyzzyaalx1(3),xyzzyaaly1(3),xyzzyaalz1,xyzzyaama1,xyzzyaamb1,&
&xyzzyaamc1,xyzzyaamd1,xyzzyaame1,xyzzyaamf1,xyzzyaamg1
logical mol_spin_density
real(dp),allocatable :: xyzzyaamh1(:,:,:,:),xyzzyaami1(:,:,:,:),xyzzya&
&amj1(:,:,:,:),xyzzyaamk1(:,:,:,:),xyzzyaaml1(:,:,:)
character(3) accumulation_method_mol_den
real(dp) xyzzyaamm1(3)
logical :: eval_dipole_moment=.false.
integer xyzzyaamn1
integer,allocatable :: xyzzyaamo1(:),xyzzyaamp1(:,:)
real(dp),allocatable :: xyzzyaamq1(:),xyzzyaamr1(:),xyzzyaams1(:),xyzz&
&yaamt1(:,:),xyzzyaamu1(:,:),xyzzyaamv1(:,:),xyzzyaamw1(:,:),xyzzyaamx&
&1(:,:),xyzzyaamy1(:,:),xyzzyaamz1(:,:),xyzzyaana1(:,:),xyzzyaanb1(:)
character(3) accumulation_method_population
logical :: eval_contact_den=.false.
integer xyzzyaanc1,xyzzyaand1
integer,allocatable :: xyzzyaane1(:),xyzzyaanf1(:,:,:),xyzzyaang1(:)
real(dp) xyzzyaanh1
real(dp),allocatable :: xyzzyaani1(:),xyzzyaanj1(:),xyzzyaank1(:),xyzz&
&yaanl1(:,:),xyzzyaanm1(:,:),xyzzyaann1(:,:),xyzzyaano1(:,:),xyzzyaanp&
&1(:,:),xyzzyaanq1(:,:),xyzzyaanr1(:,:),xyzzyaans1(:,:),xyzzyaant1(:),&
&xyzzyaanu1(:,:),xyzzyaanv1(:)
character(3) accumulation_method_twop_dm_mom
integer xyzzyaanw1,xyzzyaanx1
integer,allocatable :: xyzzyaany1(:),xyzzyaanz1(:,:,:),xyzzyaaoa1(:)
real(dp) xyzzyaaob1
real(dp),allocatable :: xyzzyaaoc1(:),xyzzyaaod1(:),xyzzyaaoe1(:),xyzz&
&yaaof1(:,:),xyzzyaaog1(:,:),xyzzyaaoh1(:,:),xyzzyaaoi1(:,:),xyzzyaaoj&
&1(:,:),xyzzyaaok1(:,:),xyzzyaaol1(:,:),xyzzyaaom1(:,:),xyzzyaaon1(:),&
&xyzzyaaoo1(:,:),xyzzyaaop1(:)
character(3) accum_method_cond_frac_mom
integer xyzzyaaoq1,xc_corr_method,xyzzyaaor1,xyzzyaaos1,nparam1,nparam&
&2
real(dp),allocatable :: xyzzyaaot1(:),xyzzyaaou1(:),xyzzyaaov1(:),a(:)&
&,xyzzyaaow1(:),xyzzyaaox1(:),xyzzyaaoy1(:),weight(:)
logical,parameter :: xyzzyaaoz1=.false.
integer,parameter :: xyzzyaapa1=50
integer :: xyzzyaapb1(xyzzyaapa1)=0,xyzzyaapc1(xyzzyaapa1)=0,xyzzyaapd&
&1(xyzzyaapa1)=0,xyzzyaape1=0
integer,allocatable :: xyzzyaapf1(:),xyzzyaapg1(:)
contains
subroutine read_expval(from_backup)
implicit none
logical,intent(in) :: from_backup
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2
logical xyzzyaaad2
xyzzyaaam1=2.d0*expval_cutoff
homogeneous_density=(trim(atom_basis_type)=='none')
if(homogeneous_density)then
do xyzzyaaab2=1,nspin
if(nele(xyzzyaaab2)>0.and.any(heg_orbtype(xyzzyaaab2,:)==0.or.heg_orbt&
&ype(xyzzyaaab2,:)>=100.or.heg_orbtype(xyzzyaaab2,:)==2))homogeneous_d&
&ensity=.false.
enddo
endif
xyzzyaabb1=periodicity/=0.and.(trim(interaction)=='mpc'.or.trim(intera&
&ction)=='ewald_mpc'.or.trim(interaction)=='mpc_ewald').and.qmc_densit&
&y_mpc
call open_units(xyzzyaaaf1,xyzzyaaae1)
if(xyzzyaaae1/=0)call errstop('READ_EXPVAL','Cannot find free io unit.&
&')
if(.not.from_backup)then
if(expvals)then
call wout('Requested expectation values:')
if(density)call wout('- density')
if(spin_density)call wout('- spin density')
if(spin_density_mat)call wout('- spin density matrix')
if(pair_corr)call wout('- reciprocal space pair correlation function')
if(pair_corr_sph)call wout('- spherical real space pair correlation fu&
&nction')
if(loc_tensor)call wout('- localization tensor')
if(structure_factor)call wout('- structure factor')
if(structure_factor_sph)call wout('- spherically-averaged structure fa&
&ctor')
if(onep_density_mat)call wout('- one-particle density matrix')
if(twop_density_mat)call wout('- two-particle density matrix')
if(cond_fraction)call wout('- condensate fraction estimator')
if(mom_den)call wout('- momentum density')
if(finite_density)call wout('- finite density')
if(eval_dipole_moment)call wout('- dipole moment (accumulated in hist &
&file)')
if(population)call wout('- population')
if(mol_density)call wout('- molecular density')
if(twop_dm_mom)call wout('- two-particle density matrix (momentum spac&
&e)')
if(cond_frac_mom)call wout('- condensate fraction estimator (momentum &
&space)')
else
call wout('None requested.')
call wout()
call wout('Reading expval.data anyway since QMC density requested for &
&MPC interaction.')
endif
endif
call xyzzyaapi1
allocate(xyzzyaaco1(no_families))
xyzzyaaco1=0
do xyzzyaaaa2=1,no_families
do xyzzyaaab2=1,nspin
if(which_fam(xyzzyaaab2)==xyzzyaaaa2)xyzzyaaco1(xyzzyaaaa2)=xyzzyaaco1&
&(xyzzyaaaa2)+nele(xyzzyaaab2)
enddo
enddo
if(.not.from_backup)then
inquire(file='expval.data',exist=xyzzyaaav1)
else
inquire(file='expval.backup',exist=xyzzyaaav1)
endif
if(xyzzyaaav1)then
call xyzzyaaph1(from_backup)
xyzzyaaaw1=.false.
if(xyzzyaaav1)then
if(density.and..not.xyzzyaabc1)xyzzyaaaw1=.true.
if(spin_density.and..not.xyzzyaabd1)xyzzyaaaw1=.true.
if(spin_density_mat.and..not.xyzzyaabe1)xyzzyaaaw1=.true.
if(pair_corr.and..not.xyzzyaabl1)xyzzyaaaw1=.true.
if(pair_corr_sph.and..not.xyzzyaabn1)xyzzyaaaw1=.true.
if(loc_tensor.and..not.xyzzyaabf1)xyzzyaaaw1=.true.
if(structure_factor.and..not.xyzzyaabh1)xyzzyaaaw1=.true.
if(structure_factor_sph.and..not.xyzzyaabj1)xyzzyaaaw1=.true.
if(onep_density_mat.and..not.xyzzyaabp1)xyzzyaaaw1=.true.
if(twop_density_mat.and..not.xyzzyaabr1)xyzzyaaaw1=.true.
if(cond_fraction.and..not.xyzzyaabt1)xyzzyaaaw1=.true.
if(mom_den.and..not.xyzzyaabv1)xyzzyaaaw1=.true.
if(finite_density.and..not.xyzzyaabx1)xyzzyaaaw1=.true.
if(mol_density.and..not.xyzzyaacb1)xyzzyaaaw1=.true.
if(population.and..not.xyzzyaabz1)xyzzyaaaw1=.true.
if(twop_dm_mom.and..not.xyzzyaacd1)xyzzyaaaw1=.true.
if(cond_frac_mom.and..not.xyzzyaacf1)xyzzyaaaw1=.true.
endif
call xyzzyaapk1
xyzzyaaac2=0
xyzzyaacw1=xyzzyaacv1
if(xyzzyaaaw1)then
if(density.or.spin_density.or.spin_density_mat.or.pair_corr.or.structu&
&re_factor)then
if(xyzzyaacv1>0)then
if(periodicity<=3.and.periodicity>=0)then
call xyzzyaaqv1(xyzzyaacs1)
else
call errstop('READ_EXPVAL','Periodicity error.')
endif
xyzzyaact1(:)=xyzzyaacs1(:)/2
if(density.or.spin_density.or.spin_density_mat)then
xyzzyaaap1=pb1
xyzzyaaaq1=pb2
xyzzyaaar1=pb3
xyzzyaaac2=xyzzyaasb1(xyzzyaact1,xyzzyaaam1)
xyzzyaaay1=.true.
do xyzzyaaaa2=1,xyzzyaacv1
if(abs(xyzzyaada1(xyzzyaaaa2)-xyzzyaaam1)<1.d-7.and.xyzzyaacx1(xyzzyaa&
&aa2)==xyzzyaaac2)then
xyzzyaaay1=.false.
endif
enddo
if(xyzzyaaay1)then
xyzzyaacw1=xyzzyaacw1+1
if(.not.from_backup)then
call wout()
call wout('Required primitive G set not in pre-existing expval.data.')
endif
else
if(.not.from_backup)then
call wout()
call wout('All required primitive G sets present in pre-existing expva&
&l.data.')
call wout('Skipping primitive G set generation step.')
endif
endif
else
xyzzyaaay1=.false.
endif
if(pair_corr.or.structure_factor)then
xyzzyaaap1=b1
xyzzyaaaq1=b2
xyzzyaaar1=b3
xyzzyaaac2=xyzzyaasb1(xyzzyaact1,xyzzyaaam1)
xyzzyaaaz1=.true.
do xyzzyaaaa2=1,xyzzyaacv1
if(abs(xyzzyaada1(xyzzyaaaa2)-xyzzyaaam1)<1.d-7.and.xyzzyaacx1(xyzzyaa&
&aa2)==xyzzyaaac2)then
xyzzyaaaz1=.false.
endif
enddo
if(xyzzyaaaz1)then
xyzzyaacw1=xyzzyaacw1+1
call wout()
call wout('Required supercell G set not in pre-existing expval.data.')
else
call wout()
call wout('All required supercell G sets present in pre-existing expva&
&l.data.')
call wout('Skipping supercell G set generation step.')
endif
else
xyzzyaaaz1=.false.
endif
else
call wout()
call wout('Required G set not already in pre-existing expval.data.')
endif
endif
endif
if(xyzzyaacv1>0)then
call xyzzyaapl1(xyzzyaaac2)
endif
if(xyzzyaabc1)call xyzzyaapm1
if(xyzzyaabd1)call xyzzyaapn1
if(xyzzyaabe1)call xyzzyaapo1
if(xyzzyaabl1)call xyzzyaapp1
if(xyzzyaabn1)call xyzzyaapq1
if(xyzzyaabf1)call xyzzyaapr1
if(xyzzyaabh1)call xyzzyaaps1
if(xyzzyaabj1)call xyzzyaapt1
if(xyzzyaabp1)call xyzzyaapu1
if(xyzzyaabr1)call xyzzyaapv1
if(xyzzyaabt1)call xyzzyaapw1
if(xyzzyaabv1)call xyzzyaapx1
if(xyzzyaabx1)call xyzzyaapy1
if(population)call xyzzyaapz1
if(xyzzyaacb1)call xyzzyaaqa1
if(xyzzyaacd1)call xyzzyaaqb1
if(xyzzyaacf1)call xyzzyaaqc1
close(xyzzyaaaf1)
if(newrun)then
xyzzyaaad2=.false.
if(xyzzyaabc1.and.density)xyzzyaaad2=.true.
if(xyzzyaabd1.and.spin_density)xyzzyaaad2=.true.
if(xyzzyaabe1.and.spin_density_mat)xyzzyaaad2=.true.
if(xyzzyaabl1.and.pair_corr.and.xyzzyaabm1)xyzzyaaad2=.true.
if(xyzzyaabn1.and.pair_corr_sph.and.xyzzyaabo1)xyzzyaaad2=.true.
if(xyzzyaabf1.and.loc_tensor.and.xyzzyaabg1)xyzzyaaad2=.true.
if(xyzzyaabh1.and.structure_factor.and.xyzzyaabi1)xyzzyaaad2=.true.
if(xyzzyaabj1.and.structure_factor_sph.and.xyzzyaabk1)xyzzyaaad2=.true&
&.
if(xyzzyaabp1.and.onep_density_mat.and.xyzzyaabq1)xyzzyaaad2=.true.
if(xyzzyaabr1.and.twop_density_mat.and.xyzzyaabs1)xyzzyaaad2=.true.
if(xyzzyaabt1.and.cond_fraction.and.xyzzyaabu1)xyzzyaaad2=.true.
if(xyzzyaabv1.and.mom_den.and.xyzzyaabw1)xyzzyaaad2=.true.
if(xyzzyaabx1.and.finite_density.and.xyzzyaaby1)xyzzyaaad2=.true.
if(xyzzyaabz1.and.population.and.xyzzyaaca1)xyzzyaaad2=.true.
if(xyzzyaacb1.and.mol_density.and.xyzzyaacc1)xyzzyaaad2=.true.
if(xyzzyaacd1.and.twop_dm_mom.and.xyzzyaace1)xyzzyaaad2=.true.
if(xyzzyaacf1.and.cond_frac_mom.and.xyzzyaacg1)xyzzyaaad2=.true.
if(xyzzyaaad2)then
call wout()
call errstop('READ_EXPVAL','Requested continued accumulation of expect&
&ation value already present in expval.data. Must restart with NEWRUN &
&set to F in input to ensure a different random number sequence.')
endif
endif
else
if(.not.from_backup)then
call wout()
call wout('No pre-existing expval.data file found.')
allocate(xyzzyaach1(1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_EXPVAL','header_text')
xyzzyaaau1=.true.
xyzzyaaag1=1
xyzzyaach1(1)='No header text'
else
call errstop('READ_EXPVAL','Attempting to restore expval.data from bac&
&kup but the backup file does not exist.')
endif
endif
call xyzzyaapj1
end subroutine read_expval
subroutine xyzzyaaph1(from_backup)
implicit none
logical,intent(in) :: from_backup
integer xyzzyaaaa3
character(80) char_80
if(.not.from_backup)then
call wout()
call wout('Pre-existing expval.data file found. Reading data.')
endif
open(unit=xyzzyaaaf1,file='expval.data',status='old',iostat=xyzzyaaae1&
&)
if(xyzzyaaae1/=0)call errstop('FIND_EXPVAL_SETS','Problem opening expv&
&al.data.')
xyzzyaaaa3=1
xyzzyaaag1=0
xyzzyaaau1=.false.
xyzzyaaas1=.false.
xyzzyaaat1=.false.
xyzzyaabc1=.false.
xyzzyaabd1=.false.
xyzzyaabe1=.false.
xyzzyaabl1=.false.
xyzzyaabn1=.false.
xyzzyaabf1=.false.
xyzzyaabh1=.false.
xyzzyaabj1=.false.
xyzzyaabp1=.false.
xyzzyaabr1=.false.
xyzzyaabt1=.false.
xyzzyaabv1=.false.
xyzzyaabx1=.false.
xyzzyaabz1=.false.
xyzzyaacb1=.false.
xyzzyaacd1=.false.
xyzzyaacf1=.false.
xyzzyaabm1=.false.
xyzzyaabo1=.false.
xyzzyaabg1=.false.
xyzzyaabi1=.false.
xyzzyaabk1=.false.
xyzzyaabq1=.false.
xyzzyaabs1=.false.
xyzzyaabu1=.false.
xyzzyaabw1=.false.
xyzzyaaby1=.false.
xyzzyaaca1=.false.
xyzzyaacc1=.false.
xyzzyaace1=.false.
xyzzyaacg1=.false.
allocate(xyzzyaach1(xyzzyaaak1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'FIND_EXPVAL_SETS','header_text')
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)then
call errstop('FIND_EXPVAL_SETS','Problem reading expval.data.')
elseif(xyzzyaaae1<0)then
if(xyzzyaaag1==0)call wout('No header is present in this file.')
exit
endif
if(trim(adjustl(char_80))=='START HEADER')then
if(xyzzyaaag1>0)call errstop('FIND_EXPVAL_SETS','There is more than on&
&e header in expval.data.')
if(.not.from_backup)call wout('HEADER:')
do xyzzyaaah1=1,xyzzyaaak1
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)xyzzyaach1(xyzzyaaah1)
if(xyzzyaaae1/=0)call errstop('FIND_EXPVAL_SETS','Problem reading expv&
&al.data.')
xyzzyaach1(xyzzyaaah1)=adjustl(xyzzyaach1(xyzzyaaah1))
if(trim(adjustl(xyzzyaach1(xyzzyaaah1)))=='END HEADER')exit
if(.not.from_backup)call wout(" "//trim(xyzzyaach1(xyzzyaaah1)))
xyzzyaaag1=xyzzyaaag1+1
if(xyzzyaaag1>=xyzzyaaak1)call errstop('FIND_EXPVAL_SETS','Header is m&
&ore than ' //trim(i2s(xyzzyaaak1-1))//' long. Try to be more concise.&
&')
enddo
xyzzyaaau1=.true.
elseif(trim(adjustl(char_80))=='START EXPVAL')then
if(.not.from_backup)call wout('EXPECTATION VALUES ALREADY PRESENT:')
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)then
call errstop('FIND_EXPVAL_SETS','Problem reading expval.data <2>.')
elseif(xyzzyaaae1<0)then
exit
endif
select case(trim(adjustl(char_80)))
case('END EXPVAL')
exit
case('START DENSITY')
if(xyzzyaabc1)call errstop('FIND_EXPVAL_SETS','More than one DENSITY b&
&lock in expval.data.')
xyzzyaabc1=.true.
if(.not.from_backup)call wout('- density')
case('START SPIN DENSITY')
if(xyzzyaabd1)call errstop('FIND_EXPVAL_SETS','More than one SPIN DENS&
&ITY block in expval.data.')
xyzzyaabd1=.true.
if(.not.from_backup)call wout('- spin density')
case('START SPIN-DENSITY MATRIX')
if(xyzzyaabe1)call errstop('FIND_EXPVAL_SETS','More than one SPIN-DENS&
&ITY MATRIX block in expval.data.')
xyzzyaabe1=.true.
if(.not.from_backup)call wout('- spin density matrix')
case('START RECIPROCAL-SPACE PCF')
if(xyzzyaabl1)call errstop('FIND_EXPVAL_SETS','More than one RECIPROCA&
&L-SPACE PCF block in expval.data.')
xyzzyaabl1=.true.
if(.not.from_backup)call wout('- reciprocal space pair correlation fun&
&ction')
case('START SPHERICAL PCF')
if(xyzzyaabn1)call errstop('FIND_EXPVAL_SETS','More than one SPHERICAL&
& PCF block in expval.data.')
xyzzyaabn1=.true.
if(.not.from_backup)call wout('- spherical real space pair correlation&
& function')
case('START LOCALIZATION TENSOR')
if(xyzzyaabf1)call errstop('FIND_EXPVAL_SETS','More than one LOCALIZAT&
&ION TENSOR block in expval.data.')
xyzzyaabf1=.true.
if(.not.from_backup)call wout('- localization tensor')
case('START STRUCTURE FACTOR')
if(xyzzyaabh1)call errstop('FIND_EXPVAL_SETS','More than one STRUCTURE&
& FACTOR block in expval.data.')
xyzzyaabh1=.true.
if(.not.from_backup)call wout('- structure factor')
case('START SPHERICAL STRUCTURE FACTOR')
if(xyzzyaabj1)call errstop('FIND_EXPVAL_SETS','More than one SPHERICAL&
& STRUCTURE FACTOR block in expval.data.')
xyzzyaabj1=.true.
if(.not.from_backup)call wout('- spherically-averaged structure factor&
&')
case('START ONE-PARTICLE DENSITY MATRIX')
if(xyzzyaabp1)call errstop('FIND_EXPVAL_SETS','More than one ONE-PARTI&
&CLE DENSITY MATRIX block in expval.data.')
xyzzyaabp1=.true.
if(.not.from_backup)call wout('- one-particle density matrix')
case('START TWO-PARTICLE DENSITY MATRIX')
if(xyzzyaabr1)call errstop('FIND_EXPVAL_SETS','More than one TWO-PARTI&
&CLE DENSITY MATRIX block in expval.data.')
xyzzyaabr1=.true.
if(.not.from_backup)call wout('- two-particle density matrix')
case('START CONDENSATE FRACTION')
if(xyzzyaabt1)call errstop('FIND_EXPVAL_SETS','More than one CONDENSAT&
&E FRACTION block in expval.data.')
xyzzyaabt1=.true.
if(.not.from_backup)call wout('- condensate fraction estimator')
case('START MOMENTUM DENSITY')
if(xyzzyaabv1)call errstop('FIND_EXPVAL_SETS','More than one MOMENTUM &
&DENSITY block in expval.data.')
xyzzyaabv1=.true.
if(.not.from_backup)call wout('- momentum density')
case('START FINITE DENSITY')
if(xyzzyaabx1)call errstop('FIND_EXPVAL_SETS','More than one FINITE DE&
&NSITY block in expval.data.')
xyzzyaabx1=.true.
if(.not.from_backup)call wout('- finite_density')
case('START POPULATION')
if(xyzzyaabz1)call errstop('FIND_EXPVAL_SETS','More than one POPULATIO&
&N block in expval.data.')
xyzzyaabz1=.true.
if(.not.from_backup)call wout('- population')
case('START MOLECULAR DENSITY')
if(xyzzyaacb1)call errstop('FIND_EXPVAL_SETS','More than one MOLECULAR&
& DENSITY block in expval.data.')
xyzzyaacb1=.true.
if(.not.from_backup)call wout('- mol_density')
case('START TWO-PARTICLE DENSITY MATRIX (MOMENTUM SPACE)')
if(xyzzyaacd1)call errstop('FIND_EXPVAL_SETS','More than one TWO-PARTI&
&CLE DENSITY MATRIX (MOMENTUM SPACE) block in expval.data.')
xyzzyaacd1=.true.
if(.not.from_backup)call wout('- two-particle density matrix (momentum&
& space)')
case('START CONDENSATE FRACTION (MOMENTUM SPACE)')
if(xyzzyaacf1)call errstop('FIND_EXPVAL_SETS','More than one CONDENSAT&
&E FRACTION (MOMENTUM SPACE) block in expval.data.')
xyzzyaacf1=.true.
if(.not.from_backup)call wout('- condensate fraction estimator (moment&
&um space)')
case default
if(trim(adjustl(char_80))=='START GVECTOR SET '//trim(i2s(xyzzyaaaa3))&
&)then
read(xyzzyaaaf1,*,iostat=xyzzyaaae1)
if(xyzzyaaae1/=0)call errstop('FIND_EXPVAL_SETS','Problem reading expv&
&al.data.')
read(xyzzyaaaf1,*,iostat=xyzzyaaae1)xyzzyaada1(xyzzyaaaa3)
if(xyzzyaaae1/=0)call errstop('FIND_EXPVAL_SETS','Problem reading expv&
&al.data.')
read(xyzzyaaaf1,*,iostat=xyzzyaaae1)
if(xyzzyaaae1/=0)call errstop('FIND_EXPVAL_SETS','Problem reading expv&
&al.data.')
read(xyzzyaaaf1,*,iostat=xyzzyaaae1)xyzzyaacx1(xyzzyaaaa3)
if(xyzzyaaae1/=0)call errstop('FIND_EXPVAL_SETS','Problem reading expv&
&al.data.')
xyzzyaaaa3=xyzzyaaaa3+1
endif
end select
enddo
xyzzyaaas1=.true.
endif
enddo
rewind(xyzzyaaaf1)
if(.not.xyzzyaaas1)then
call errstop('FIND_EXPVAL_SETS','The expval.data file must contain a"S&
&TART EXPVAL" block.')
endif
if(.not.xyzzyaaau1)then
if(.not.allocated(xyzzyaach1))allocate(xyzzyaach1(1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'FIND_EXPVAL_SETS','header_text')
xyzzyaaau1=.true.
xyzzyaaag1=1
xyzzyaach1(1)='No title given.'
endif
end subroutine xyzzyaaph1
subroutine xyzzyaapi1
implicit none
if(density.and.homogeneous_density)then
call errstop('CHECK_EXPVAL1','System has homogeneous density - no need&
& to accumulate it.')
endif
if(spin_density.and.homogeneous_density)then
call errstop('CHECK_EXPVAL1','System has homogeneous spin densities - &
&no need to accumulate them.')
endif
if(spin_density_mat.and..not.noncoll_spin) call errstop('CHECK_EXPVAL1&
&','Can accumulate spin density matrix only in systems with non-collin&
&ear spins.')
if(spin_density.and.noncoll_spin) call errstop('CHECK_EXPVAL1','In non&
&-collinear system should accumulate spin density matrix rather than s&
&pin density.')
if(spin_density_mat.and.(isdmc.or.isvmc_dmc.or.isdmc_dmc))call errstop&
&('CHECK_EXPVAL1','DMC accumulation of spin density matrix not current&
&ly possible (no DMC algorithm for non-collinear spins).')
if(pair_corr.and.homogeneous_density)then
call wout()
call errwarn('CHECK_EXPVAL1','Calculating the reciprocal-space pair-co&
&rrelation function for a system with homogeneous density, but you pro&
&bably want to use the real-space spherical algorithm (activated with &
&the PAIR_CORR_SPH keyword) instead.')
endif
if(pair_corr.and.pair_corr_sph)call errstop('CHECK_EXPVAL1','Cannot si&
&multaneously accumulate spherical and reciprocal space PCFs.')
if(pair_corr.and.(density.or.spin_density.or.loc_tensor.or.structure_f&
&actor.or.structure_factor_sph))then
call wout()
call wordwrap('Reciprocal-space pair correlation function must not be &
&accumulated simultaneously with other expectation values, as algorith&
&m requires particles to be fixed in position.')
call errstop('CHECK_EXPVAL1','Quitting.')
endif
if(pair_corr_sph.and..not.homogeneous_density)then
if(input_type_fixed/=-999)then
call errstop('CHECK_EXPVAL1','Accumulation of spherical pair correlati&
&on function with fixed particles has not yet been implemented.')
endif
if(density.or.spin_density.or.loc_tensor.or.structure_factor.or.struct&
&ure_factor_sph)then
call wout()
call wordwrap('Spherical pair correlation function for inhomogeneous s&
&ystems must not be accumulated simultaneously with other expectation &
&values, as algorithm requires particles to be fixed in position.')
call errstop('CHECK_EXPVAL1','Quitting.')
endif
endif
if(isperiodic)then
if((pair_corr.or.(pair_corr_sph.and..not.homogeneous_density)).and.inp&
&ut_type_fixed/=-999)then
xyzzyaaai1=size(nele)
if(input_type_fixed<1.or.input_type_fixed>xyzzyaaai1)then
call wout()
call wout('Type of fixed particle in PCF_RFIX block in input file is o&
&ut of range.')
if(xyzzyaaai1==1)then
call wout('For this system you should set it to 1.')
elseif(xyzzyaaai1==2)then
call wout('Must be either 1 or '//trim(i2s(xyzzyaaai1))//'.')
else
call wout('Must be between 1 and '//trim(i2s(xyzzyaaai1))//'.')
endif
call errstop('CHECK_EXPVAL1','Quitting.')
endif
if(nele(input_type_fixed)==0)then
call wout()
call wordwrap('Requested pair-correlation around fixed particle of typ&
&e ' //trim(i2s(input_type_fixed))//', but there are no such particles&
&.')
call errstop('CHECK_EXPVAL1','Quitting.')
endif
endif
endif
if(.not.isperiodic)then
if(density)call errstop('CHECK_EXPVAL1','Density accumulation for fini&
&te systems such as atoms and ions done through finite_density. Please&
& file a bug report.')
if(spin_density)call errstop('CHECK_EXPVAL1','Spin-density accumulatio&
&n for finite systems currently non-functional.')
if(pair_corr)call errstop('CHECK_EXPVAL1','Pair-correlation function m&
&ay not be calculated for finite systems.')
if(loc_tensor)call errstop('CHECK_EXPVAL1','Localization-tensoraccumul&
&ation for finite systems currently non-functional.')
if(structure_factor)call errstop('CHECK_EXPVAL1','Structure factor may&
& not be calculated for finite systems.')
if(spin_density_mat)call errstop('CHECK_EXPVAL1','Spin-density matrix &
&accumulation for finite systems currently non-functional.')
if(onep_density_mat)call errstop('CHECK_EXPVAL1','One-particle density&
& matrix accumulation for finite systems currently non-functional.')
if(twop_density_mat)call errstop('CHECK_EXPVAL1','Two-particle density&
& matrix accumulation for finite systems currently non-functional.')
if(cond_fraction)call errstop('CHECK_EXPVAL1','Condensate fraction est&
&imator accumulation for finite systems currently non-functional.')
if(twop_dm_mom)call errstop('CHECK_EXPVAL1','Two-particle density matr&
&ix (momentum space) accumulation for finite systems currently non-fun&
&ctional.')
if(cond_frac_mom)call errstop('CHECK_EXPVAL1','Condensate fraction (mo&
&mentum space) estimator accumulation for finite systems currently non&
&-functional.')
endif
if(onep_density_mat.and..not.homogeneous_density)then
call errstop('CHECK_EXPVAL1','One-particle density matrix for inhomoge&
&neous systems not implemented.')
endif
if(twop_density_mat.and..not.homogeneous_density)then
call errstop('CHECK_EXPVAL1','Two-particle density matrix for inhomoge&
&neous systems not implemented.')
endif
if(cond_fraction.and..not.homogeneous_density)then
call errstop('CHECK_EXPVAL1','Condensate fraction estimator for inhomo&
&geneous systems not implemented.')
endif
if(structure_factor_sph.and..not.homogeneous_density)call errstop('CHE&
&CK_EXPVAL1','Spherical structure factor for inhomogeneous systems not&
& implemented.')
if(mom_den.and..not.homogeneous_density)call errstop('CHECK_EXPVAL1','&
&Momentum density for inhomogeneous systems not implemented.')
if(finite_density.and.isperiodic)call errstop('CHECK_EXPVAL1','Finite &
&density mistakenly flagged for periodic system. This is a bug.')
if(finite_density.and.nitot/=1)call errstop('CHECK_EXPVAL1','Finite de&
&nsity accumulation flagged for molecule: this facility is only availa&
&ble for single atoms.')
if(finite_density.and.nspin/=2)call errstop('CHECK_EXPVAL1','Finite de&
&nsity accumulation only available for systems with electrons only.')
if(mol_density.and.isperiodic)call errstop('CHECK_EXPVAL1','Molecular &
&density mistakenly flagged for periodic system. This is a bug.')
if(finite_density.and.nspin/=2)call errstop('CHECK_EXPVAL1','Molecular&
& density accumulation only available for systems with electrons only.&
&')
if(twop_dm_mom.and..not.homogeneous_density)call errstop('CHECK_EXPVAL&
&1','Two-particle density matrix (momentum space) for inhomogeneous sy&
&stems not implemented.')
if(cond_frac_mom.and..not.homogeneous_density)call errstop('CHECK_EXPV&
&AL1','Condensate fraction estimator (momentum space) for inhomogeneou&
&s systems not implemented.')
end subroutine xyzzyaapi1
subroutine xyzzyaapj1
implicit none
integer xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5
character(80) tmpr,tmpr2,tmpr3
if(xyzzyaaav1)then
if(xyzzyaabb1.and..not.xyzzyaabc1)then
call errstop('CHECK_EXPVAL2','Requested QMC density to be used for MPC&
& interaction, but no density block in expval.data.')
endif
if(int_sf.and..not.structure_factor.and..not.xyzzyaabh1)then
call errstop('CHECK_EXPVAL2','INT_SF is flagged in input so either an &
&expval.data file containing the accumulated structure factor should e&
&xist, or STRUCTURE_FACTOR must be set to T in input.')
endif
if(pair_corr.and..not.xyzzyaabl1.and.input_type_fixed==-999)then
call errstop('CHECK_EXPVAL2','Fixed particle coords for PCF not given &
&in input or in expval.data. Use the PCF_RFIX input block to do this.'&
&)
endif
if(isperiodic)then
if(pair_corr_sph.and..not.xyzzyaabn1.and.input_type_fixed==-999.and..n&
&ot.homogeneous_density)then
call errstop('CHECK_EXPVAL2','Accumulation of spherical PCF not yetimp&
&lemented for inhomogeneous periodic systems.')
call errstop('CHECK_EXPVAL2','Fixed particle coords for PCF not given &
&in input or in expval.data. Use the PCF_RFIX input block to do this.'&
&)
endif
endif
if(structure_factor_sph.and..not.xyzzyaabj1.and.expval_nkdim(2)==-999)&
&then
call errstop('CHECK_EXPVAL2','Radial k point grid for spherical struct&
&ure factor is given neither in input (keyword EXPVAL_KGRID) nor in ex&
&pval.data.')
endif
if(mom_den.and.xyzzyaabv1.and.onep_density_mat.and.xyzzyaabp1)then
if(xyzzyaagz1/=xyzzyaahs1.or.xyzzyaahb1/=xyzzyaaht1)call errstop('CHEC&
&K_EXPVAL2','When the one-particle density matrix and the momentum den&
&sity are both present in expval.data, the number of sets and of rando&
&m points must be equal in both of them.')
do xyzzyaaac5=1,xyzzyaagz1
do xyzzyaaaa5=1,nspin
if(xyzzyaahd1(xyzzyaaaa5,xyzzyaaac5)/=xyzzyaahv1(xyzzyaaaa5,xyzzyaaac5&
&))call errstop('CHECK_EXPVAL2','When the one-particle density matrix &
&and the momentum density are both present in expval.data, the particl&
&e types associated with each set must be the same, and must be given &
&in the same order.')
enddo
enddo
endif
if(twop_density_mat.and.xyzzyaabr1.and.twop_dm_mom.and.xyzzyaacd1)then
if(xyzzyaaik1/=xyzzyaaik1.or.xyzzyaaim1/=xyzzyaand1.or.xyzzyaaiv1/=xyz&
&zyaanh1)call errstop('CHECK_EXPVAL2','When the two-particle density m&
&atrix and the two-particle density matrix (momentum space) are both p&
&resent in expval.data, the number of sets, random points and fraction&
& of pairs to smaple must be equal in both of them.')
do xyzzyaaac5=1,xyzzyaaik1
do xyzzyaaab5=1,max_spin_pairs
if(xyzzyaaio1(1,xyzzyaaab5,xyzzyaaac5)/=xyzzyaanf1(1,xyzzyaaab5,xyzzya&
&aac5).or.xyzzyaaio1(2,xyzzyaaab5,xyzzyaaac5)/=xyzzyaanf1(2,xyzzyaaab5&
&,xyzzyaaac5))call errstop('CHECK_EXPVAL2','When the two-particle dens&
&ity matrix and the two-particle density matrix (momentum space) are b&
&oth present in expval.data, the particle pair types associated with e&
&ach set must be the same, and must be given in the same order.')
enddo
enddo
endif
else
if(pair_corr)then
if(input_type_fixed/=-999)then
call wout()
call wout('Coordinates of fixed particle given in input for reciprocal&
&-space PCF.')
tmpr=r2s(input_rfix(1),'(f18.6)')
tmpr2=r2s(input_rfix(2),'(f18.6)')
tmpr3=r2s(input_rfix(3),'(f18.6)')
call wout('Particle of type '//trim(i2s(input_type_fixed))//' fixed at&
& '//trim(tmpr)//', '//trim(tmpr2)//', '//trim(tmpr3)//'.')
else
call errstop('CHECK_EXPVAL2','Fixed particle coords for PCF not given &
&in input or in expval.data. Use the PCF_RFIX input block to do this.'&
&)
endif
endif
if(pair_corr_sph)then
call wout()
if(isperiodic)then
if(homogeneous_density)then
call wout('Using homogeneous density algorithm for spherical PCF.')
else
call errstop('CHECK_EXPVAL2','Accumulation of spherical PCF not yetimp&
&lemented for inhomogeneous periodic systems.')
call wout('Using fixed particle algorithm for spherical PCF (inhomogen&
&eous density).')
if(input_type_fixed/=-999)then
call wout('Coordinates given in input.')
tmpr=r2s(input_rfix(1),'(f18.6)')
tmpr2=r2s(input_rfix(2),'(f18.6)')
tmpr3=r2s(input_rfix(3),'(f18.6)')
call wout('Particle of type '//trim(i2s(input_type_fixed))//' fixed at&
& '//trim(tmpr)//', '//trim(tmpr2)//', '//trim(tmpr3)//'.')
else
call errstop('CHECK_EXPVAL2','Fixed particle coords for PCF not given &
&in input or in expval.data. Use the PCF_RFIX input block to do this.'&
&)
endif
endif
else
if(trim(atom_basis_type)=='none')then
call wout('Using homogeneous density algorithm in finite system for sp&
&herical PCF.')
if(input_type_fixed/=-999)then
call wout()
call wout('Ignoring fixed particle information in PCF_RFIX block in in&
&put.')
endif
else
call errstop('CHECK_EXPVAL2','Accumulation of spherical PCF not implem&
&ented for finite systems containing atoms.')
endif
endif
endif
if(xyzzyaabb1)then
call errstop('CHECK_EXPVAL2','Requested QMC density to be used for MPC&
& interation, but expval.data does not exist.')
endif
if(int_sf.and..not.structure_factor)then
call errstop('CHECK_EXPVAL2','INT_SF is flagged in input so either an &
&expval.data file containing the accumulated structure factor should e&
&xist, or STRUCTURE_FACTOR must be set to T in input.')
endif
if(structure_factor_sph)then
if(expval_nkdim(2)/=-999)then
if(expval_nkdim(2)>1)call errstop('CHECK_EXPVAL2','One-dimensional k g&
&rid required for spherical structure factor. Check EXPVAL_KGRID input&
&.')
call wout()
call wout('Radial k point grid given in input for spherical structure &
&factor:')
xyzzyaago1=expval_ka(1,2)
xyzzyaagp1=expval_kb(1,2)
xyzzyaagl1=expval_nk(1,2)
tmpr=r2s(xyzzyaago1,'(f18.6)')
tmpr2=r2s(xyzzyaagp1,'(f18.6)')
call wout('Evenly-spaced grid of '//trim(i2s(xyzzyaagl1))//' k-points &
&from '//trim(tmpr)//' to '//trim(tmpr2)//'.')
else
call errstop('CHECK_EXPVAL2','Radial k point grid for spherical struct&
&ure factor not given in input (keyword EXPVAL_KGRID) or in expval.dat&
&a.')
endif
endif
endif
if(isperiodic)then
if(pair_corr_sph.and.input_type_fixed/=-999.and.homogeneous_density)ca&
&ll errwarn('CHECK_EXPVAL2','Fixed particle coords for spherical PCF i&
&n input PCF_RFIX block will be ignored as system has homogeneous dens&
&ity.')
endif
end subroutine xyzzyaapj1
subroutine xyzzyaapk1
implicit none
integer xyzzyaaaa6,xyzzyaaab6,xyzzyaaac6
real(dp) xyzzyaaad6(3),xyzzyaaae6(3),xyzzyaaaf6(3),xyzzyaaag6
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaac6)char_80
if(xyzzyaaac6>0)call errstop('READ_EXPVAL_BLOCK','Problem reading expv&
&al.data.')
if(xyzzyaaac6<0)call errstop('READ_EXPVAL_BLOCK','Could not find "STAR&
&T EXPVAL" in expval.data.')
if(trim(adjustl(char_80))=='START EXPVAL')exit
enddo
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
char_80=trim(adjustl(char_80))
read(xyzzyaaaf1,'(a)',err=1,end=1)title
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaack1
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,iostat=xyzzyaaac6)xyzzyaacl1
if(xyzzyaaac6/=0)then
rewind(xyzzyaaaf1)
return
endif
xyzzyaaat1=.true.
if(xyzzyaacl1>nspin)call errstop2('READ_EXPVAL_BLOCK','Too many partic&
&le types in expval.data: ',xyzzyaacl1)
allocate(xyzzyaacn1(xyzzyaacl1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_EXPVAL_BLOCK','nele_expval')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)(xyzzyaacn1(xyzzyaaaa6),xyzzyaaaa6=1,xyz&
&zyaacl1)
if(any(xyzzyaacn1(:)/=nele(:)))call errstop('READ_EXPVAL_BLOCK','Misma&
&tch in number of each type of particle in expval.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaacm1
if(xyzzyaacm1/=dimensionality)call errstop('READ_EXPVAL_BLOCK','Dimens&
&ionality mismatch between wave function and expval.data file')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaci1
if(xyzzyaaci1/=periodicity)call errstop('READ_EXPVAL_BLOCK','Periodici&
&ty mismatch between wave function and expval.data file')
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab6=1,3
do xyzzyaaaa6=1,3
if(xyzzyaaaa6==xyzzyaaab6)then
xyzzyaacj1(xyzzyaaaa6,xyzzyaaab6)=1
else
xyzzyaacj1(xyzzyaaaa6,xyzzyaaab6)=0
endif
enddo
enddo
select case(xyzzyaaci1)
case(0)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)
case(1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaad6(1)
if(abs(xyzzyaaad6(1)-pa1(1))>1.d-8)call errstop('READ_EXPVAL_BLOCK','P&
&rimitive translation vector mismatch in expval.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaacj1(1,1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaacq1
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaag6
if(abs(xyzzyaaag6-wigner_seitz_radius)>1.d-8)call errstop('READ_EXPVAL&
&_BLOCK','Mismatch of radius of line inscribed in Wigner-Seitz cell in&
& expval.data.')
case(2)
read(xyzzyaaaf1,*,err=1,end=1)(xyzzyaaad6(xyzzyaaaa6),xyzzyaaaa6=1,2)
read(xyzzyaaaf1,*,err=1,end=1)(xyzzyaaae6(xyzzyaaaa6),xyzzyaaaa6=1,2)
if(any(abs(xyzzyaaad6(1:2)-pa1(1:2))>1.d-8).or.any(abs(xyzzyaaae6(1:2)&
&-pa2(1:2))>1.d-8))call errstop('READ_EXPVAL_BLOCK','Primitive transla&
&tion vector mismatch in expval.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
read(char_80,*,iostat=xyzzyaaac6)xyzzyaacj1(1,1),xyzzyaacj1(2,2),xyzzy&
&aacj1(1,2),xyzzyaacj1(2,1)
if(xyzzyaaac6/=0)then
xyzzyaacj1(1,2)=0
xyzzyaacj1(2,1)=0
read(char_80,*,err=1,end=1)xyzzyaacj1(1,1),xyzzyaacj1(2,2)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaacp1
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaag6
if(abs(xyzzyaaag6-wigner_seitz_radius)>1.d-8)call errstop('READ_EXPVAL&
&_BLOCK','Mismatch of radius of circle inscribed in Wigner-Seitz cell &
&in expval.data.')
case(3)
read(xyzzyaaaf1,*,err=1,end=1)(xyzzyaaad6(xyzzyaaaa6),xyzzyaaaa6=1,3)
read(xyzzyaaaf1,*,err=1,end=1)(xyzzyaaae6(xyzzyaaaa6),xyzzyaaaa6=1,3)
read(xyzzyaaaf1,*,err=1,end=1)(xyzzyaaaf6(xyzzyaaaa6),xyzzyaaaa6=1,3)
if(any(abs(xyzzyaaad6(:)-pa1(:))>1.d-8).or.any(abs(xyzzyaaae6(:)-pa2(:&
&))>1.d-8).or.any(abs(xyzzyaaaf6(:)-pa3(:))>1.d-8))call errstop('READ_&
&EXPVAL_BLOCK','Primitive translation vector mismatch in expval.data.'&
&)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
read(char_80,*,iostat=xyzzyaaac6)xyzzyaacj1(1,1),xyzzyaacj1(2,2),xyzzy&
&aacj1(3,3),xyzzyaacj1(1,2),xyzzyaacj1(1,3),xyzzyaacj1(2,1),xyzzyaacj1&
&(2,3),xyzzyaacj1(3,1),xyzzyaacj1(3,2)
if(xyzzyaaac6/=0)then
xyzzyaacj1=0
read(char_80,*,err=1,end=1)xyzzyaacj1(1,1),xyzzyaacj1(2,2),xyzzyaacj1(&
&3,3)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaacr1
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaag6
if(abs(xyzzyaaag6-wigner_seitz_radius)>1.d-8)call errstop('READ_EXPVAL&
&_BLOCK','Mismatch of radius of sphere inscribed in Wigner-Seitz cell &
&in expval.data.')
case default
call errstop('READ_EXPVAL_BLOCK','Invalid periodicity.')
end select
if(any(scell_matrix/=xyzzyaacj1))call errstop('READ_EXPVAL_BLOCK','SCE&
&LL_MATRIX mismatch between input and expval.data files.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaacv1
if(xyzzyaacv1>9)call errstop('READ_EXPVAL_BLOCK','No more than nine G &
&vector sets may be defined in expval.data.')
rewind(xyzzyaaaf1)
return
1 call errstop('READ_EXPVAL_BLOCK','Problem reading START EXPVAL secti&
&on of expval.data file.')
end subroutine xyzzyaapk1
subroutine xyzzyaapl1(ng)
implicit none
integer,intent(in) :: ng
integer xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7
character(80) char_80
allocate(xyzzyaadb1(xyzzyaacw1),expval_ngvec(xyzzyaacw1),xyzzyaadc1(3,&
&3,xyzzyaacw1),xyzzyaadd1(3,3,xyzzyaacw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_G_VECTOR_SETS','e_cutoff,expval_ngve&
&c,expval_latvec,expval_rlatvec')
if(xyzzyaabb1)then
allocate(expval_ngvec_truncated(xyzzyaacw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_G_VECTOR_SETS','expval_ngvec_truncat&
&ed')
endif
do xyzzyaaaa7=1,xyzzyaacv1
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_G_VECTOR_SETS','Problem reading exp&
&val.data.')
if(xyzzyaaae1<0)call errstop2('READ_G_VECTOR_SETS','Could not find "ST&
&ART GVECTOR SET x" in expval.data for x = ',xyzzyaaaa7)
if(trim(adjustl(char_80))=='START GVECTOR SET '//trim(i2s(xyzzyaaaa7))&
&)exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)expval_ngvec(xyzzyaaaa7)
enddo
xyzzyaacu1=max(ng,maxval(expval_ngvec(1:xyzzyaacv1)))
allocate(expval_gvec(3,xyzzyaacu1,xyzzyaacw1),xyzzyaacy1(3,xyzzyaacu1,&
&xyzzyaacw1),xyzzyaacz1(xyzzyaacw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_G_VECTOR_SETS','expval_gvec,pgmap,gr&
&ange')
rewind(xyzzyaaaf1)
do xyzzyaaaa7=1,xyzzyaacv1
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(trim(adjustl(char_80))=='START GVECTOR SET '//trim(i2s(xyzzyaaaa7))&
&)exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaadb1(xyzzyaaaa7)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)expval_ngvec(xyzzyaaaa7)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab7=1,3
read(xyzzyaaaf1,*,err=1,end=1)(xyzzyaadd1(xyzzyaaac7,xyzzyaaab7,xyzzya&
&aaa7),xyzzyaaac7=1,3)
enddo
if(abs(xyzzyaadd1(1,1,xyzzyaaaa7)-pb1(1))<1.d-8.and.abs(xyzzyaadd1(2,1&
&,xyzzyaaaa7)-pb1(2))<1.d-8.and.abs(xyzzyaadd1(3,1,xyzzyaaaa7)-pb1(3))&
&<1.d-8)then
xyzzyaadc1(1:3,1,xyzzyaaaa7)=pa1(1:3)
xyzzyaadc1(1:3,2,xyzzyaaaa7)=pa2(1:3)
xyzzyaadc1(1:3,3,xyzzyaaaa7)=pa3(1:3)
elseif(abs(xyzzyaadd1(1,1,xyzzyaaaa7)-b1(1))<1.d-8.and.abs(xyzzyaadd1(&
&2,1,xyzzyaaaa7)-b1(2))<1.d-8.and.abs(xyzzyaadd1(3,1,xyzzyaaaa7)-b1(3)&
&)<1.d-8)then
xyzzyaadc1(1:3,1,xyzzyaaaa7)=a1(1:3)
xyzzyaadc1(1:3,2,xyzzyaaaa7)=a2(1:3)
xyzzyaadc1(1:3,3,xyzzyaaaa7)=a3(1:3)
else
call errstop('READ_G_VECTOR_SETS','Set of reciprocal lattice vectors i&
&n expval.data do not correspond to the primitive cell or supercell.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab7=1,expval_ngvec(xyzzyaaaa7)
read(xyzzyaaaf1,*,err=1,end=1)(expval_gvec(xyzzyaaac7,xyzzyaaab7,xyzzy&
&aaaa7),xyzzyaaac7=1,3)
enddo
expval_gvec(:,expval_ngvec(xyzzyaaaa7)+1:xyzzyaacu1,xyzzyaaaa7)=0.d0
enddo
allocate(xyzzyaadf1(xyzzyaacu1,xyzzyaacw1,1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_G_VECTOR_SETS','expval_expigdotr')
xyzzyaadf1(:,:,:)=czero
rewind(xyzzyaaaf1)
return
1 call errstop('READ_G_VECTOR_SETS','Problem reading G vector sets in &
&expval.data file.')
end subroutine xyzzyaapl1
subroutine xyzzyaapm1
implicit none
integer xyzzyaaaa8,xyzzyaaab8
character(2) nset
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_DENSITY','Problem reading expval.da&
&ta file.')
if(xyzzyaaae1<0)call errstop('READ_DENSITY','Could not find START DENS&
&ITY in expval.data file.')
if(trim(adjustl(char_80))=='START DENSITY')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_den=trim(adjustl(char_80))
if(accumulation_method_den/='VMC'.and.accumulation_method_den/='DMC')c&
&all errstop('READ_DENSITY','Density accumulation method (VMC or DMC) &
&in expval.data not recognized.')
if(density)then
if(isvmc.and.accumulation_method_den=='DMC')call errstop('READ_DENSITY&
&','DMC density accumulation on top of VMC-accumulated density in expv&
&al.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_den=='VMC&
&')call errstop('READ_DENSITY','VMC density accumulation on top of DMC&
&-accumulated density in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaadg1
if(xyzzyaadg1<1.or.xyzzyaadg1>xyzzyaacv1)call errstop('READ_DENSITY','&
&Invalid specification of which G vector set to use.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)den_nsets
if(den_nsets>nspin)then
call errstop('READ_DENSITY','More density sets in expval.data than def&
&ined types of particle.')
endif
allocate(xyzzyaadj1(den_nsets),xyzzyaadi1(den_nsets),xyzzyaadh1(den_ns&
&ets),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_DENSITY','expval_den_weight')
allocate(xyzzyaadl1(expval_ngvec(xyzzyaadg1),den_nsets),xyzzyaadk1(exp&
&val_ngvec(xyzzyaadg1),den_nsets),xyzzyaadm1(expval_ngvec(xyzzyaadg1),&
&den_nsets),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_DENSITY','expval_den')
if(xyzzyaabb1)then
allocate(expval_den_mpc(expval_ngvec(xyzzyaadg1),den_nsets),stat=xyzzy&
&aaad1)
call check_alloc(xyzzyaaad1,'READ_DENSITY','expval_den_mpc')
endif
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaadn1(expval_ngvec(xyzzyaadg1),den_nsets),stat=xyzzyaaad&
&1)
call check_alloc(xyzzyaaad1,'READ_DENSITY','expval_den_dmc')
endif
xyzzyaadi1(:)=0.d0
xyzzyaadk1(:,:)=czero
expval_complex_den=.false.
if(.not.inversion_symmetry)expval_complex_den=.true.
do xyzzyaaaa8=1,den_nsets
write(nset,'(i2)')xyzzyaaaa8
nset=trim(adjustl(nset))
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//nset)then
call errstop2('READ_DENSITY','Error reading START SET x line for x = '&
&,xyzzyaaaa8)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaadh1(xyzzyaaaa8)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaadj1(xyzzyaaaa8)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab8=1,expval_ngvec(xyzzyaadg1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaao1,xyzzyaaan1
xyzzyaadl1(xyzzyaaab8,xyzzyaaaa8)=cmplx(xyzzyaaao1,xyzzyaaan1,kind=dp)
enddo
if(xyzzyaabb1)then
do xyzzyaaab8=1,expval_ngvec(xyzzyaadg1)
expval_den_mpc(xyzzyaaab8,xyzzyaaaa8)=xyzzyaadl1(xyzzyaaab8,xyzzyaaaa8&
&)
enddo
endif
if(any(aimag(xyzzyaadl1(:,xyzzyaaaa8))>tol_expval_diff).and..not.expva&
&l_complex_den)call errstop('READ_DENSITY','Structure has inversion sy&
&mmetry but density Fourier coeffs in expval.data are complex.')
read(xyzzyaaaf1,*,err=1,end=1)
enddo
return
1 call errstop('READ_DENSITY','Problem reading DENSITY section of expv&
&al.data file.')
end subroutine xyzzyaapm1
subroutine xyzzyaapn1
implicit none
integer xyzzyaaaa9,xyzzyaaab9
character(2) nset
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_SPIN_DENSITY','Problem reading expv&
&al.data file.')
if(xyzzyaaae1<0)call errstop('READ_SPIN_DENSITY','Could not find START&
& SPIN DENSITY in expval.data file.')
if(trim(adjustl(char_80))=='START SPIN DENSITY')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_sden=trim(adjustl(char_80))
if(accumulation_method_sden/='VMC'.and.accumulation_method_sden/='DMC'&
&)call errstop('READ_SPIN_DENSITY','Spin density accumulation method (&
&VMC or DMC) in expval.data not recognized.')
if(spin_density)then
if(isvmc.and.accumulation_method_sden=='DMC')call errstop('READ_SPIN_D&
&ENSITY','DMC spin density accumulation on top of VMC-accumulated spin&
& density in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_sden=='VM&
&C')call errstop('READ_SPIN_DENSITY','VMC spin density accumulation on&
& top of DMC-accumulated spin density in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaado1
if(xyzzyaado1<1.or.xyzzyaado1>xyzzyaacv1)call errstop('READ_SPIN_DENSI&
&TY','Invalid specification of which G vector set to use.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaadp1
if(xyzzyaadp1>nspin)then
call errstop('READ_SPIN_DENSITY','More spin density sets in expval.dat&
&a than defined types of particle.')
endif
allocate(xyzzyaads1(xyzzyaadp1),xyzzyaadr1(xyzzyaadp1),xyzzyaadq1(xyzz&
&yaadp1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_SPIN_DENSITY','expval_sden_weight ar&
&rays')
allocate(xyzzyaadu1(expval_ngvec(xyzzyaado1),xyzzyaadp1),xyzzyaadt1(ex&
&pval_ngvec(xyzzyaado1),xyzzyaadp1),xyzzyaadv1(expval_ngvec(xyzzyaado1&
&),xyzzyaadp1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_SPIN_DENSITY','expval_sden')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaadw1(expval_ngvec(xyzzyaado1),xyzzyaadp1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'READ_SPIN_DENSITY','expval_sden_dmc')
endif
xyzzyaadr1(:)=0.d0
xyzzyaadt1(:,:)=czero
if(.not.inversion_symmetry)then
xyzzyaadx1=.true.
else
xyzzyaadx1=.false.
endif
do xyzzyaaaa9=1,xyzzyaadp1
write(nset,'(i2)')xyzzyaaaa9
nset=trim(adjustl(nset))
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//nset)call errstop2('READ_SPIN&
&_DENSITY','Error reading START SET x line for x = ',xyzzyaaaa9)
read(xyzzyaaaf1,*,err=1,end=1)
if(nele(xyzzyaaaa9)>0)then
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaadq1(xyzzyaaaa9)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaads1(xyzzyaaaa9)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab9=1,expval_ngvec(xyzzyaado1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaao1,xyzzyaaan1
xyzzyaadu1(xyzzyaaab9,xyzzyaaaa9)=cmplx(xyzzyaaao1,xyzzyaaan1,kind=dp)
enddo
if(any(aimag(xyzzyaadu1(:,xyzzyaaaa9))>tol_expval_diff).and..not.xyzzy&
&aadx1)call errstop('READ_SPIN_DENSITY','Structure has inversion symme&
&try but spin density Fourier coefficients in expval.data are complex.&
&')
read(xyzzyaaaf1,*,err=1,end=1)
else
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaads1(xyzzyaaaa9)=0.d0
xyzzyaadu1(1:expval_ngvec(xyzzyaado1),xyzzyaaaa9)=czero
endif
enddo
return
1 call errstop('READ_SPIN_DENSITY','Problem reading SPIN DENSITY secti&
&on of expval.data file.')
end subroutine xyzzyaapn1
subroutine xyzzyaapo1
implicit none
integer xyzzyaaaa10,xyzzyaaab10,xyzzyaaac10
character(2) nset
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_SPIN_DENSITY_MATRIX','Problem readi&
&ng expval.data file.')
if(xyzzyaaae1<0)call errstop('READ_SPIN_DENSITY_MATRIX','Could not fin&
&d START SPIN-DENSITY MATRIX in expval.data file.')
if(trim(adjustl(char_80))=='START SPIN-DENSITY MATRIX')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_sdenm=trim(adjustl(char_80))
if(accumulation_method_sdenm/='VMC'.and.accumulation_method_sdenm/='DM&
&C')call errstop('READ_SPIN_DENSITY_MATRIX','Spin density matrix accum&
&ulation method (VMC or DMC) in expval.data not recognized.')
if(spin_density_mat)then
if(isvmc.and.accumulation_method_sdenm=='DMC')call errstop('READ_SPIN_&
&DENSITY_MATRIX','DMC spin density matrix accumulation on top of VMC-a&
&ccumulated spin density matrix in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_sdenm=='V&
&MC') call errstop('READ_SPIN_DENSITY_MATRIX','VMC spin density matrix&
& accumulation on top of DMC-accumulated spin density matrix in expval&
&.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaady1
if(xyzzyaady1<1.or.xyzzyaady1>xyzzyaacv1)call errstop('READ_SPIN_DENSI&
&TY_MATRIX','Invalid specification of which G vector set to use.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaadz1
if(xyzzyaadz1>nspin)then
call errstop('READ_SPIN_DENSITY_MATRIX','More spin density matrix sets&
& in expval.data than defined types of particle.')
endif
allocate(xyzzyaaec1(xyzzyaadz1),xyzzyaaeb1(xyzzyaadz1),xyzzyaaea1(xyzz&
&yaadz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_SPIN_DENSITY_MATRIX','expval_sdenm_w&
&eight')
allocate(xyzzyaaee1(expval_ngvec(xyzzyaady1),4,xyzzyaadz1),xyzzyaaed1(&
&expval_ngvec(xyzzyaady1),4,xyzzyaadz1),xyzzyaaef1(expval_ngvec(xyzzya&
&ady1),4,xyzzyaadz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_SPIN_DENSITY_MATRIX','expval_sdenm')
xyzzyaaeb1(:)=0.d0
xyzzyaaed1(:,:,:)=czero
xyzzyaaeg1=.true.
do xyzzyaaaa10=1,xyzzyaadz1
write(nset,'(i2)')xyzzyaaaa10
nset=trim(adjustl(nset))
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//nset)then
call errstop2('READ_SPIN_DENSITY_MATRIX','Error reading START SET x li&
&ne for x = ',xyzzyaaaa10)
endif
read(xyzzyaaaf1,*,err=1,end=1)
if(nele(xyzzyaaaa10)>0)then
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaea1(xyzzyaaaa10)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaec1(xyzzyaaaa10)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaac10=1,4
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
select case(xyzzyaaac10)
case(1)
if(trim(adjustl(char_80))/='1 1')call errstop('READ_SPIN_DENSITY_MATRI&
&X','Error reading component 1 1.')
case(2)
if(trim(adjustl(char_80))/='1 2')call errstop('READ_SPIN_DENSITY_MATRI&
&X','Error reading component 1 2.')
case(3)
if(trim(adjustl(char_80))/='2 1')call errstop('READ_SPIN_DENSITY_MATRI&
&X','Error reading component 2 1.')
case(4)
if(trim(adjustl(char_80))/='2 2')call errstop('READ_SPIN_DENSITY_MATRI&
&X','Error reading component 2 2.')
end select
do xyzzyaaab10=1,expval_ngvec(xyzzyaady1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaao1,xyzzyaaan1
xyzzyaaee1(xyzzyaaab10,xyzzyaaac10,xyzzyaaaa10)=cmplx(xyzzyaaao1,xyzzy&
&aaan1,kind=dp)
enddo
enddo
if(any(aimag(xyzzyaaee1(:,:,xyzzyaaaa10))>tol_expval_diff).and. .not.x&
&yzzyaaeg1)call errstop('READ_SPIN_DENSITY_MATRIX','Structure has inve&
&rsion symmetry but spin density mat Fourier coeffs in expval.data are&
& complex.')
read(xyzzyaaaf1,*,err=1,end=1)
else
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaaec1(xyzzyaaaa10)=0.d0
xyzzyaaee1(1:expval_ngvec(xyzzyaady1),1:4,xyzzyaaaa10)=czero
endif
enddo
return
1 call errstop('READ_SPIN_DENSITY_MATRIX','Problem reading SPIN-DENSIT&
&Y MATRIX section of expval.data file.')
end subroutine xyzzyaapo1
subroutine xyzzyaapp1
implicit none
integer xyzzyaaaa11,xyzzyaaab11
character(2) nset
character(80) char_80,tmpr,tmpr2,tmpr3
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_PAIR_CORR','Problem reading expval.&
&data file.')
if(xyzzyaaae1<0)call errstop('READ_PAIR_CORR','Could not find START RE&
&CIPROCAL-SPACE PCF in expval.data file.')
if(trim(adjustl(char_80))=='START RECIPROCAL-SPACE PCF')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_pcf=trim(adjustl(char_80))
if(accumulation_method_pcf/='VMC'.and.accumulation_method_pcf/='DMC')c&
&all errstop('READ_PAIR_CORR','Reciprocal-space PCF accumulation metho&
&d (VMC or DMC) in expval.data not recognized.')
if(pair_corr)then
if(isvmc.and.accumulation_method_pcf=='DMC')call errstop('READ_PAIR_CO&
&RR','DMC PCF accumulation on top of VMC-accumulated PCF in expval.dat&
&a.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_pcf=='VMC&
&')call errstop('READ_PAIR_CORR','VMC PCF accumulation on top of DMC-a&
&ccumulated PCF in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaej1
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)(pcf_rfix(xyzzyaaaa11),xyzzyaaaa11=1,3)
xyzzyaaai1=size(nele)
if(xyzzyaaej1<1.or.xyzzyaaej1>xyzzyaaai1)then
call wout()
call wout('Type of fixed particle in expval.data PAIR_CORR block is ou&
&t of range.')
if(xyzzyaaai1==1)then
call wout('For this system it must be 1.')
elseif(xyzzyaaai1==2)then
call wout('Must be either 1 or '//trim(i2s(xyzzyaaai1))//'.')
else
call wout('Must be between 1 and '//trim(i2s(xyzzyaaai1))//'.')
endif
call errstop('READ_PAIR_CORR','Quitting.')
endif
if(nele(xyzzyaaej1)==0)then
call wout()
call wout('Existing pair-correlation function around fixed particle of&
& type '//trim(i2s(xyzzyaaej1))//', but')
call wout('there are no such particles.')
call errstop('READ_PAIR_COOR','Quitting.')
endif
if(input_type_fixed==-999)then
call wout()
call wout('Coordinates of fixed particle given in expval.data for reci&
&procal-space PCF.')
tmpr=r2s(pcf_rfix(1),'(f18.6)')
tmpr2=r2s(pcf_rfix(2),'(f18.6)')
tmpr3=r2s(pcf_rfix(3),'(f18.6)')
call wout('Particle of type '//trim(i2s(xyzzyaaej1))//' fixed at '//tr&
&im(tmpr)//', '//trim(tmpr2)//', '//trim(tmpr3)//'.')
else
if(any(pcf_rfix(:)/=input_rfix(:)))then
call errstop('READ_PAIR_CORR','Coordinates of fixed point given in inp&
&ut and in pre-existing expval.data do not match.')
endif
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaeh1
if(xyzzyaaeh1<1.or.xyzzyaaeh1>xyzzyaacv1)call errstop('READ_PAIR_CORR'&
&,'Invalid specification of which G vector set to use for reciprocal s&
&pace PCF.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaei1
if(xyzzyaaei1>nspin)then
call errstop('READ_PAIR_CORR','More RECIPROCAL-SPACE PCF sets in expva&
&l.data than defined types of particle.')
endif
allocate(xyzzyaaem1(xyzzyaaei1),xyzzyaael1(xyzzyaaei1),xyzzyaaek1(xyzz&
&yaaei1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_PAIR_CORR','expval_pcf_weight')
allocate(xyzzyaaen1(expval_ngvec(xyzzyaaeh1),xyzzyaaei1),expval_pcf(ex&
&pval_ngvec(xyzzyaaeh1),xyzzyaaei1),xyzzyaaeo1(expval_ngvec(xyzzyaaeh1&
&),xyzzyaaei1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_PAIR_CORR','expval_pcf')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaaep1(expval_ngvec(xyzzyaaeh1),xyzzyaaei1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'READ_PAIR_CORR','expval_pcf_dmc')
endif
xyzzyaael1(:)=0.d0
expval_pcf(:,:)=czero
xyzzyaaeq1=.false.
if(.not.inversion_symmetry)xyzzyaaeq1=.true.
loop_sets: do xyzzyaaaa11=1,xyzzyaaei1
write(nset,'(i2)')xyzzyaaaa11
nset=trim(adjustl(nset))
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//nset)then
call errstop2('READ_PAIR_CORR','Error reading START SET x line for x =&
& ',xyzzyaaaa11)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaek1(xyzzyaaaa11)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaem1(xyzzyaaaa11)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab11=1,expval_ngvec(xyzzyaaeh1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaaa11)))then
xyzzyaaem1(xyzzyaaaa11)=0.d0
xyzzyaaen1(:,xyzzyaaaa11)=czero
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaaen1(xyzzyaaab11,xyzzyaaaa11)
enddo
if(.not.xyzzyaaeq1)then
if(any(aimag(xyzzyaaen1(:,xyzzyaaaa11))>tol_expval_diff).and..not.xyzz&
&yaaeq1)call errstop('READ_PAIR_CORR','Structure has inversion symmetr&
&y but PCF Fourier coeffs in expval.data are complex.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaabm1=.true.
enddo loop_sets
return
1 call errstop('READ_PAIR_CORR','Problem reading RECIPROCAL-SPACE PCF &
&section of expval.data file.')
end subroutine xyzzyaapp1
subroutine xyzzyaapq1
implicit none
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12
character(2) nset
character(80) char_80,tmpr,tmpr2,tmpr3
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_PAIR_CORR_SPH','Problem reading exp&
&val.data file.')
if(xyzzyaaae1<0)call errstop('READ_PAIR_CORR_SPH','Could not find STAR&
&T SPHERICAL PCF in expval.data file.')
if(trim(adjustl(char_80))=='START SPHERICAL PCF')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_pcfs=trim(adjustl(char_80))
if(accumulation_method_pcfs/='VMC'.and.accumulation_method_pcfs/='DMC'&
&)call errstop('READ_PAIR_CORR_SPH','Spherical PCF accumulation method&
& (VMC or DMC) in expval.data not recognized.')
if(pair_corr_sph)then
if(isvmc.and.accumulation_method_pcfs=='DMC')call errstop('READ_PAIR_C&
&ORR_SPH','DMC spherical PCF accumulation on top of VMC-accumulated sp&
&herical PCF in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_pcfs=='VM&
&C') call errstop('READ_PAIR_CORR_SPH','VMC spherical PCF accumulation&
& on top of DMC-accumulated spherical PCF in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaes1
if(xyzzyaaes1<1)call errstop('READ_PAIR_CORR_SPH','Need at least one b&
&in for accumulation of spherical PCF in expval.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaew1
if(xyzzyaaew1<=0.d0)call errstop('READ_PAIR_CORR_SPH','Invalid cutoff &
&radius for spherical PCF in expval.data.')
if(isperiodic.and.xyzzyaaew1>wigner_seitz_radius)call errstop ('READ_P&
&AIR_CORR_SPH','Cutoff radius for spherical PCF in expval.data > radiu&
&s of sphere inscribed in Wigner-Seitz cell.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaet1
if(.not.isperiodic)then
if(xyzzyaaet1/=2)call errstop('READ_PAIR_CORR_SPH','Accumulation mode &
&must be 2 (no fixed particles) for finite systems.')
else
if(xyzzyaaet1==2.and..not.homogeneous_density)call errstop ('READ_PAIR&
&_CORR_SPH','Accumulation mode set to homogeneous in expval.data but c&
&urrent system is inhomogeneous.')
endif
if(xyzzyaaet1/=1.and.xyzzyaaet1/=2)call errstop('READ_PAIR_CORR_SPH','&
&Invalid accumulation mode for spherical PCF in expval.data.')
if(xyzzyaaet1==1)then
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaeu1
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)(pcfs_rfix(xyzzyaaaa12),xyzzyaaaa12=1,3)
if(input_type_fixed==-999)then
call wout()
call wout('Coordinates of fixed particle given in expval.data for sphe&
&rical PCF.')
tmpr=r2s(pcfs_rfix(1),'(f18.6)')
tmpr2=r2s(pcfs_rfix(2),'(f18.6)')
tmpr3=r2s(pcfs_rfix(3),'(f18.6)')
call wout('Particle of type '//trim(i2s(xyzzyaaeu1))//' fixed at '//tr&
&im(tmpr)//', '//trim(tmpr2)//', '//trim(tmpr3)//'.')
if(nele(xyzzyaaeu1)==0)then
call wout()
call wout('But there are no defined particles of this type..')
call errstop('READ_EXPVAL','Quitting.')
endif
else
if(any(pcfs_rfix(:)/=input_rfix(:)))then
call errstop('READ_PAIR_CORR_SPH','Fixed particle coords given in inpu&
&t and in expval.data do not match.')
endif
endif
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaer1
if(xyzzyaaet1==1)then
if(xyzzyaaer1>nspin)then
call errstop('READ_PAIR_CORR_SPH','More SPHERICAL PCF sets in expval.d&
&ata than number of pairs for fixed particle.')
endif
else
xyzzyaaab12=0
do xyzzyaaaa12=1,nspin
xyzzyaaab12=xyzzyaaab12+xyzzyaaaa12
enddo
if(xyzzyaaer1>xyzzyaaab12)then
call errstop('READ_PAIR_CORR_SPH','More SPHERICAL PCF sets in expval.d&
&ata than types of particle pairs.')
endif
endif
allocate(xyzzyaafb1(xyzzyaaer1),xyzzyaaev1(2,xyzzyaaer1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'READ_PAIR_CORR_SPH','expval_pcfs_weight_s&
&tart,pcfs_ptype')
if(xyzzyaaet1==1)then
allocate(xyzzyaafd1(xyzzyaaes1,1,nspin),xyzzyaafc1(xyzzyaaes1,1,nspin)&
&,xyzzyaafe1(xyzzyaaes1,1,nspin),stat=xyzzyaaad1)
else
allocate(xyzzyaafd1(xyzzyaaes1,nspin,nspin),xyzzyaafc1(xyzzyaaes1,nspi&
&n,nspin),xyzzyaafe1(xyzzyaaes1,nspin,nspin),stat=xyzzyaaad1)
endif
call check_alloc(xyzzyaaad1,'READ_PAIR_CORR_SPH','expval_pcfs')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
if(xyzzyaaet1==1)then
allocate(xyzzyaaff1(xyzzyaaes1,1,nspin),stat=xyzzyaaad1)
else
allocate(xyzzyaaff1(xyzzyaaes1,nspin,nspin),stat=xyzzyaaad1)
endif
call check_alloc(xyzzyaaad1,'READ_PAIR_CORR_SPH','expval_pcfs_dmc')
endif
xyzzyaafb1(:)=0.d0
xyzzyaaey1=0.d0
xyzzyaaez1=0.d0
xyzzyaafd1(:,:,:)=0.d0
xyzzyaafc1(:,:,:)=0.d0
loop_sets: do xyzzyaaad12=1,xyzzyaaer1
write(nset,'(i2)')xyzzyaaad12
nset=trim(adjustl(nset))
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//nset)then
call errstop2('READ_PAIR_CORR_SPH','Error reading START SET x line for&
& x = ',xyzzyaaad12)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaev1(1,xyzzyaaad12),xyzzyaaev1(2,x&
&yzzyaaad12)
if(xyzzyaaet1==1)then
if(xyzzyaaev1(1,xyzzyaaad12)/=xyzzyaaeu1)then
call wout()
call wordwrap('First particle type in pair must equal the type of the &
&fixed particle when accumulating the spherical PCF in inhomogeneous s&
&ystems.')
call errstop('READ_PAIR_CORR_SPH','Quitting.')
endif
endif
if(any(xyzzyaaev1(:,xyzzyaaad12)<1).or.any(xyzzyaaev1(:,xyzzyaaad12)>n&
&spin))then
call errstop('READ_PAIR_CORR_SPH','Types of particle in pair are out o&
&f range in spherical PCF block in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaafb1(xyzzyaaad12)
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaaab12=xyzzyaaev1(1,xyzzyaaad12)
xyzzyaaac12=xyzzyaaev1(2,xyzzyaaad12)
do xyzzyaaaa12=1,xyzzyaaes1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaad12)))then
xyzzyaafb1(xyzzyaaad12)=0.d0
xyzzyaafd1(:,xyzzyaaab12,xyzzyaaac12)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaafd1(xyzzyaaaa12,xyzzyaaab12,xyzzyaaa&
&c12)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaabo1=.true.
enddo loop_sets
return
1 call errstop('READ_PAIR_CORR_SPH','Problem reading SPHERICAL PCF sec&
&tion of expval.data file.')
end subroutine xyzzyaapq1
subroutine xyzzyaapr1
implicit none
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13
character(80) char_80
xyzzyaafg1=nspin
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_LOC_TENSOR','Problem readingexpval.&
&data file.')
if(xyzzyaaae1<0)call errstop('READ_LOC_TENSOR','Could not find START L&
&OCALIZATION TENSOR in expval.data file.')
if(trim(adjustl(char_80))=='START LOCALIZATION TENSOR')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)' ,err=1,end=1)char_80
accumulation_method_lt=trim(adjustl(char_80))
if(accumulation_method_lt/='VMC' .and. accumulation_method_lt/='DMC')c&
&all errstop('READ_LOC_TENSOR','Localization tensor accumulation metho&
&d (VMC or DMC) in expval.data not recognized.' )
if(loc_tensor)then
if(isvmc.and.accumulation_method_lt/='VMC')call errstop('READ_LOC_TENS&
&OR','VMC localization tensor accumulation on top of DMC-accumulated l&
&ocalization tensor in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_lt/='DMC'&
&)call errstop('READ_LOC_TENSOR','DMC localization tensor accumulation&
& on top of VMC-accumulated localization tensor in expval.data.')
endif
allocate(xyzzyaafm1(periodicity,periodicity,xyzzyaafg1),xyzzyaafi1(xyz&
&zyaafg1),xyzzyaafl1(periodicity,periodicity,xyzzyaafg1),xyzzyaafh1(xy&
&zzyaafg1),xyzzyaafn1(periodicity,periodicity,xyzzyaafg1),xyzzyaafj1(x&
&yzzyaafg1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_LOC_TENSOR','expval_lt_*')
xyzzyaafi1(:)=0.d0
xyzzyaafh1(:)=0.d0
xyzzyaafj1(:)=0.d0
xyzzyaafm1(:,:,:)=czero
xyzzyaafl1(:,:,:)=czero
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaafo1(periodicity,periodicity,xyzzyaafg1),xyzzyaafk1(xyz&
&zyaafg1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_LOC_TENSOR','expval_lt_dmc')
xyzzyaafo1=czero
xyzzyaafk1=0.d0
endif
xyzzyaafp1=.not.inversion_symmetry
do xyzzyaaaa13=1,xyzzyaafg1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaaa13)))call e&
&rrstop2('READ_LOC_TENSOR','Error reading START SET x line for x = ',x&
&yzzyaaaa13)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaafi1(xyzzyaaaa13)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab13=1,periodicity
do xyzzyaaac13=1,periodicity
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaao1,xyzzyaaan1
xyzzyaafm1(xyzzyaaac13,xyzzyaaab13,xyzzyaaaa13)=cmplx(xyzzyaaao1,xyzzy&
&aaan1,dp)
enddo
enddo
if(any(aimag(xyzzyaafm1(1:periodicity,1:periodicity,xyzzyaaaa13))>tol_&
&expval_diff).and..not.xyzzyaafp1)call errstop('READ_LOC_TENSOR','Stru&
&cture has inversion symmetrybut <ei(k.r)> has Fourier coeffs in expva&
&l.data are complex.')
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab13=1,periodicity
read(xyzzyaaaf1,*,err=1,end=1)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='END SET '//trim(i2s(xyzzyaaaa13)))call err&
&stop('READ_LOC_TENSOR','Error reading END SET '//trim(i2s(xyzzyaaaa13&
&))//' line.')
enddo
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='END LOCALIZATION TENSOR')call errstop('REA&
&D_LOC_TENSOR','Error reading END LOCALIZATION TENSOR line.')
return
1 call errstop('READ_LOC_TENSOR','Problem reading LOCALIZATION TENSOR &
&section of expval.data file.')
end subroutine xyzzyaapr1
subroutine xyzzyaaps1
implicit none
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14
character(2) nset
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_STRUCTURE_FACTOR','Problem reading &
&expval.data file.')
if(xyzzyaaae1<0)call errstop('READ_STRUCTURE_FACTOR','Could not find S&
&TART STRUCTURE FACTOR in expval.data file.')
if(trim(adjustl(char_80))=='START STRUCTURE FACTOR')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_sf=trim(adjustl(char_80))
if(accumulation_method_sf/='VMC'.and.accumulation_method_sf/='DMC') ca&
&ll errstop('READ_STRUCTURE_FACTOR','Structure factor accumulation met&
&hod (VMC or DMC) in expval.data not recognized.')
if(structure_factor)then
if(isvmc.and.accumulation_method_sf=='DMC')call errstop('READ_STRUCTUR&
&E_FACTOR','DMC structure factor accumulation on top of VMC-accumulate&
&d structure factor in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_sf=='VMC'&
&)call errstop('READ_STRUCTURE_FACTOR','VMC structure factor accumulat&
&ion on top of DMC-accumulated structure factor in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaafr1
if(xyzzyaafr1<1.or.xyzzyaafr1>xyzzyaacv1) call errstop('READ_STRUCTURE&
&_FACTOR','Invalid specification of which G vector set to use.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaafq1
xyzzyaaab14=0
do xyzzyaaaa14=1,nspin
xyzzyaaab14=xyzzyaaab14+xyzzyaaaa14
enddo
if(xyzzyaafq1>xyzzyaaab14)then
call errstop('READ_STRUCTURE_FACTOR','More structure factor sets in ex&
&pval.data than defined types of particle pairs.')
endif
allocate(xyzzyaafy1(xyzzyaafq1),xyzzyaaft1(2,xyzzyaafq1),xyzzyaafz1(ex&
&pval_ngvec(xyzzyaafr1),xyzzyaafq1),xyzzyaaga1(expval_ngvec(xyzzyaafr1&
&),xyzzyaafq1),xyzzyaagf1(expval_ngvec(xyzzyaafr1),nspin),xyzzyaagb1(e&
&xpval_ngvec(xyzzyaafr1),xyzzyaafq1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_STRUCTURE_FACTOR','expval_sf arrays'&
&)
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagc1(expval_ngvec(xyzzyaafr1),xyzzyaafq1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'READ_STRUCTURE_FACTOR','expval_sf_dmc')
endif
xyzzyaafy1(:)=0.d0
xyzzyaafv1=0.d0
xyzzyaafw1=0.d0
xyzzyaafz1(:,:)=0.d0
xyzzyaaga1(:,:)=0.d0
xyzzyaagf1(:,:)=czero
loop_sets: do xyzzyaaac14=1,xyzzyaafq1
write(nset,'(i2)')xyzzyaaac14
nset=trim(adjustl(nset))
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//nset)then
call errstop2('READ_STRUCTURE_FACTOR','Error reading START SET x line &
&for x = ',xyzzyaaac14)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaft1(1,xyzzyaaac14),xyzzyaaft1(2,x&
&yzzyaaac14)
read(xyzzyaaaf1,*,err=1,end=1)
if(nele(xyzzyaaft1(1,xyzzyaaac14))/=0.and.nele(xyzzyaaft1(2,xyzzyaaac1&
&4))/=0)then
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaafy1(xyzzyaaac14)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaab14=1,expval_ngvec(xyzzyaafr1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaac14)))then
xyzzyaafy1(xyzzyaaac14)=0.d0
xyzzyaafz1(:,xyzzyaaac14)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaafz1(xyzzyaaab14,xyzzyaaac14)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaabi1=.true.
else
read(xyzzyaaaf1,*,err=1,end=1)
endif
enddo loop_sets
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaafs1
if(homogeneous_density.and.xyzzyaafs1/=0)then
call errstop('READ_STRUCTURE_FACTOR','Density is homogeneous but the S&
&TRUCTURE FACTOR block in expval.data contains an accumulated spin den&
&sity.')
endif
if(.not.homogeneous_density.and.xyzzyaafs1==0)then
call errstop('READ_STRUCTURE_FACTOR','Density is inhomogeneous but the&
& STRUCTURE FACTOR block in expval.data does not contain an accumulate&
&d spin density.')
endif
if(xyzzyaafs1>nspin)then
call errstop('READ_STRUCTURE_FACTOR','More spin density sets (for stru&
&cture factor) in expval.data than defined types of particle.')
endif
if(.not.homogeneous_density)then
allocate(xyzzyaagd1(xyzzyaafs1),xyzzyaage1(xyzzyaafs1),xyzzyaafu1(xyzz&
&yaafs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_STRUCTURE_FACTOR','expval_sfsden_wei&
&ght,sfsden_ptype')
allocate(xyzzyaagh1(expval_ngvec(xyzzyaafr1),xyzzyaafs1),xyzzyaagg1(ex&
&pval_ngvec(xyzzyaafr1),xyzzyaafs1),xyzzyaagi1(expval_ngvec(xyzzyaafr1&
&),xyzzyaafs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_STRUCTURE_FACTOR','expval_sfsden')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagj1(expval_ngvec(xyzzyaafr1),xyzzyaafs1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'READ_STRUCTURE_FACTOR','expval_sfsden_dmc&
&')
endif
xyzzyaagd1(:)=0.d0
xyzzyaage1(:)=0.d0
xyzzyaagh1(:,:)=czero
xyzzyaagg1(:,:)=czero
if(.not.inversion_symmetry)then
xyzzyaagk1=.true.
else
xyzzyaagk1=.false.
endif
loop_sets2: do xyzzyaaac14=1,xyzzyaafs1
write(nset,'(i2)')xyzzyaaac14
nset=trim(adjustl(nset))
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//nset)then
call errstop2('READ_STRUCTURE_FACTOR','Error reading START SET x line &
&for spin density part of structure factor for x = ',xyzzyaaac14)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaafu1(xyzzyaaac14)
if(nele(xyzzyaafu1(xyzzyaaac14))>0)then
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaagd1(xyzzyaaac14)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaaa14=1,expval_ngvec(xyzzyaafr1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaac14)))then
xyzzyaagd1(xyzzyaaac14)=0.d0
xyzzyaagh1(:,xyzzyaaac14)=czero
cycle loop_sets2
endif
read(char_80,*,err=1,end=1)xyzzyaaao1,xyzzyaaan1
xyzzyaagh1(xyzzyaaaa14,xyzzyaaac14)=cmplx(xyzzyaaao1,xyzzyaaan1,kind=d&
&p)
enddo
if(.not.xyzzyaagk1)then
if(any(aimag(xyzzyaagh1(:,xyzzyaaac14))>tol_expval_diff))call errstop(&
&'READ_STRUCTURE_FACTOR','Structure has inversion symmetry but spin de&
&nsity (for structure factor) Fourier coefficients in expval.data are &
&complex.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaabi1=.true.
else
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)
endif
enddo loop_sets2
endif
return
1 call errstop('READ_STRUCTURE_FACTOR','Problem reading STRUCTURE FACT&
&OR section of expval.data file.')
end subroutine xyzzyaaps1
subroutine xyzzyaapt1
implicit none
integer xyzzyaaaa15,xyzzyaaab15
character(80) char_80,tmpr,tmpr2
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_STRUCTURE_FACTOR_SPH','Problem read&
&ing expval.data file.')
if(xyzzyaaae1<0)call errstop('READ_STRUCTURE_FACTOR_SPH','Could not fi&
&nd START SPHERICAL STRUCTURE FACTOR in expval.data file.')
if(trim(adjustl(char_80))=='START SPHERICAL STRUCTURE FACTOR')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_sf_sph=trim(adjustl(char_80))
if(accumulation_method_sf_sph/='VMC'.and.accumulation_method_sf_sph/='&
&DMC') call errstop('READ_STRUCTURE_FACTOR_SPH','Spherical structure f&
&actor accumulation method (VMC or DMC) in expval.data not recognized.&
&')
if(structure_factor_sph)then
if(isvmc.and.accumulation_method_sf_sph=='DMC')call errstop('READ_STRU&
&CTURE_FACTOR_SPH','DMC spherical structure factor accumulation on top&
& of VMC-accumulated spherical structure factor in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_sf_sph=='&
&VMC') call errstop('READ_STRUCTURE_FACTOR_SPH','VMC spherical structu&
&re factor accumulation on top of DMC-accumulated spherical structure &
&factor in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaago1,xyzzyaagp1,xyzzyaagl1
if(expval_nk(1,2)==-999)then
call wout()
call wout('Radial k point grid given in expval.data for spherical stru&
&cture factor:')
tmpr=r2s(xyzzyaago1,'(f18.6)')
tmpr2=r2s(xyzzyaagp1,'(f18.6)')
call wout('Grid of '//trim(i2s(xyzzyaagl1))//' k points from '//trim(t&
&mpr)//' to '//trim(tmpr2)//'.')
else
if(abs(expval_ka(1,2)-xyzzyaago1)>1.d-10.or.abs(expval_kb(1,2)-xyzzyaa&
&gp1)>1.d-10.or.expval_nk(1,2)/=xyzzyaagl1)then
call errstop('READ_STRUCTURE_FACTOR_SPH','Radial k point grids given i&
&n input and in expval.data do not match.')
endif
endif
allocate(xyzzyaagt1(xyzzyaagl1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_STRUCTURE_FACTOR_SPH','expval_sf_sph&
&_k')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaagm1
if(xyzzyaagm1/=no_spairs(levels_spairs))call errstop('READ_STRUCTURE_F&
&ACTOR_SPH','At present the number of spherical structure factor sets &
&must equal the number of spin-pairs.')
allocate(xyzzyaagq1(xyzzyaagm1),xyzzyaags1(xyzzyaagm1),xyzzyaagn1(2,xy&
&zzyaagm1),xyzzyaagu1(xyzzyaagl1,xyzzyaagm1),xyzzyaagv1(xyzzyaagl1,xyz&
&zyaagm1),xyzzyaagy1(xyzzyaagm1),xyzzyaagx1(xyzzyaagl1,xyzzyaagm1),sta&
&t=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_STRUCTURE_FACTOR_SPH','expval_sf_sph&
&')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagw1(xyzzyaagl1,xyzzyaagm1),xyzzyaagr1(xyzzyaagm1),stat&
&=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_STRUCTURE_FACTOR_SPH','expval_sf_sph&
&_dmc')
endif
xyzzyaagx1(:,:)=0.d0
xyzzyaagy1(:)=0.d0
xyzzyaagu1(:,:)=0.d0
xyzzyaagq1(:)=0.d0
loop_sets: do xyzzyaaab15=1,xyzzyaagm1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaab15)))call e&
&rrstop('READ_STRUCTURE_FACTOR_SPH','Error reading START SET ' //trim(&
&i2s(xyzzyaaab15))//'.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaagn1(1:2,xyzzyaaab15)
if(xyzzyaaab15/=which_spair(xyzzyaagn1(1,xyzzyaaab15),xyzzyaagn1(2,xyz&
&zyaaab15),levels_spairs)) call errstop('READ_STRUCTURE_FACTOR_SPH','P&
&article pair not as expected in set '//trim(i2s(xyzzyaaab15))//'.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaagy1(xyzzyaaab15)
if(xyzzyaagy1(xyzzyaaab15)<0.d0)call errstop('READ_STRUCTURE_FACTOR_SP&
&H','Negative weight.')
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaaa15=1,xyzzyaagl1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaab15)))then
xyzzyaagy1(xyzzyaaab15)=0.d0
xyzzyaagx1(:,xyzzyaaab15)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaagt1(xyzzyaaaa15),xyzzyaagx1(xyzzyaaa&
&a15,xyzzyaaab15)
enddo
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='END SET '//trim(i2s(xyzzyaaab15)))call err&
&stop('READ_STRUCTURE_FACTOR_SPH','Error reading END SET ' //trim(i2s(&
&xyzzyaaab15))//'.')
xyzzyaabk1=.true.
enddo loop_sets
return
1 call errstop('READ_STRUCTURE_FACTOR_SPH','Problem reading SPHERICAL &
&STRUCTURE FACTOR section of expval.data file.')
end subroutine xyzzyaapt1
subroutine xyzzyaapu1
implicit none
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_ONEP_DENSITY_MATRIX','Problem readi&
&ng expval.data file.')
if(xyzzyaaae1<0)call errstop('READ_ONEP_DENSITY_MATRIX','Could not fin&
&d START ONE-PARTICLE DENSITY MATRIX in expval.data file.')
if(trim(adjustl(char_80))=='START ONE-PARTICLE DENSITY MATRIX')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_onep_dm=trim(adjustl(char_80))
if(accumulation_method_onep_dm/='VMC'.and.accumulation_method_onep_dm/&
&='DMC')call errstop('READ_ONEP_DENSITY_MATRIX','One-particle density &
&matrix accumulation method (VMC or DMC) in expval.data not recognized&
&.')
if(onep_density_mat)then
if(isvmc.and.accumulation_method_onep_dm=='DMC')call errstop('READ_ONE&
&P_DENSITY_MATRIX','DMC one-particle density matrix accumulation on to&
&p of VMC-accumulated one-particle density matrix in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_onep_dm==&
&'VMC') call errstop('READ_ONEP_DENSITY_MATRIX','VMC one-particle dens&
&ity matrix accumulation on top of DMC-accumulated one-particle densit&
&y matrix in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaagz1
if(on_top_ii==0)then
if(xyzzyaagz1>nspin)call errstop('READ_ONEP_DENSITY_MATRIX','More one-&
&particle density matrix sets in expval.data than defined types of par&
&ticle.')
else
if(xyzzyaagz1>1)call errstop('READ_ONEP_DENSITY_MATRIX','More than one&
& one-particle density matrix sets in expval.data, but this is an ON_T&
&OP_PAIR calculation.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaha1
if(xyzzyaaha1<2)call errstop('READ_ONEP_DENSITY_MATRIX','Need at least&
& two bins for accumulation of one-particle density matrix in expval.d&
&ata.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaahb1
if(xyzzyaahb1<1)call errstop('READ_ONEP_DENSITY_MATRIX','Need at least&
& one random point for accumulation of one-particle density matrix in &
&expval.data.')
allocate(xyzzyaahg1(xyzzyaaha1,xyzzyaagz1),xyzzyaahk1(xyzzyaaha1,xyzzy&
&aagz1),xyzzyaahl1(xyzzyaaha1,xyzzyaagz1),xyzzyaahf1(xyzzyaaha1,xyzzya&
&agz1),xyzzyaahi1(xyzzyaaha1,xyzzyaagz1),xyzzyaahj1(xyzzyaaha1,xyzzyaa&
&gz1),xyzzyaahh1(xyzzyaaha1,xyzzyaagz1),xyzzyaahm1(xyzzyaaha1,xyzzyaag&
&z1),xyzzyaahn1(xyzzyaaha1,xyzzyaagz1),xyzzyaahc1(xyzzyaagz1),xyzzyaah&
&d1(nspin,xyzzyaagz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_ONEP_DENSITY_MATRIX','expval_onep_dm&
&')
xyzzyaahg1=0.d0
xyzzyaahk1=0.d0
xyzzyaahl1=0.d0
xyzzyaahf1=0.d0
xyzzyaahi1=0.d0
xyzzyaahj1=0.d0
xyzzyaahh1=0.d0
xyzzyaahm1=0.d0
xyzzyaahn1=0.d0
xyzzyaahc1=0
xyzzyaahd1=0
if(accumulation_method_onep_dm=='DMC')then
allocate(xyzzyaahr1(xyzzyaaha1,xyzzyaagz1),xyzzyaahp1(xyzzyaaha1,xyzzy&
&aagz1),xyzzyaahq1(xyzzyaaha1,xyzzyaagz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_ONEP_DENSITY_MATRIX','expval_onep_dm&
&_*_dmc')
xyzzyaahr1=0.d0
xyzzyaahp1=0.d0
xyzzyaahq1=0.d0
endif
loop_sets: do xyzzyaaab16=1,xyzzyaagz1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaab16)))then
call errstop2('READ_ONEP_DENSITY_MATRIX','Error reading START SET x li&
&ne for x = ',xyzzyaaab16)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='on-top pair')then
if(on_top_ii==0)call errstop('READ_ONEP_DENSITY_MATRIX','Attempting to&
& continue an on-top density matrix accumulation, but ON_TOP_PAIR not &
&set.')
xyzzyaahc1(xyzzyaaab16)=1
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaac16,xyzzyaaad16
if(xyzzyaaac16/=which_spin(on_top_ii).or.xyzzyaaad16/=which_spin(on_to&
&p_jj))then
if(xyzzyaaad16==which_spin(on_top_ii).and.xyzzyaaac16==which_spin(on_t&
&op_jj))call errstop('READ_ONEP_DENSITY_MATRIX','Attempting to continu&
&e on-top density matrix accumulation, but ON_TOP_PAIR declaration has&
& lines in reverse order, which could cause issues.')
call errstop('READ_ONEP_DENSITY_MATRIX','Attempting to continue on-top&
& density matrix accumulation, but ON_TOP_PAIR specifies different par&
&ticle types from original run.')
endif
xyzzyaahd1(1,xyzzyaaab16)=which_spin(on_top_ii)
else
read(char_80,*,err=1)xyzzyaahc1(xyzzyaaab16)
if(xyzzyaahc1(xyzzyaaab16)<1.or.xyzzyaahc1(xyzzyaaab16)>nspin)call err&
&stop2('READ_ONEP_DENSITY_MATRIX','Problem with number of particle typ&
&es in set ',xyzzyaaab16)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaahd1(1:xyzzyaahc1(xyzzyaaab16),xyz&
&zyaaab16)
if(any(xyzzyaahd1(1:xyzzyaahc1(xyzzyaaab16),xyzzyaaab16)<1).or.any(xyz&
&zyaahd1(1:xyzzyaahc1(xyzzyaaab16),xyzzyaaab16)>nspin))call errstop2('&
&READ_ONEP_DENSITY_MATRIX','Particle type out of range in particle lis&
&t for set ',xyzzyaaab16)
endif
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaaa16=1,xyzzyaaha1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaab16)))then
xyzzyaahg1(:,xyzzyaaab16)=0.d0
xyzzyaahk1(:,xyzzyaaab16)=0.d0
xyzzyaahl1(:,xyzzyaaab16)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaahg1(xyzzyaaaa16,xyzzyaaab16),xyzzyaa&
&hl1(xyzzyaaaa16,xyzzyaaab16),xyzzyaahk1(xyzzyaaaa16,xyzzyaaab16)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaabq1=.true.
enddo loop_sets
return
1 call errstop('READ_ONEP_DENSITY_MATRIX','Problem reading ONE-PARTICL&
&E DENSITY MATRIX section of expval.data file.')
end subroutine xyzzyaapu1
subroutine xyzzyaapv1
implicit none
integer xyzzyaaaa17,xyzzyaaab17
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_TWOP_DENSITY_MATRIX','Problem readi&
&ng expval.data file.')
if(xyzzyaaae1<0)call errstop('READ_TWOP_DENSITY_MATRIX','Could not fin&
&d START TWO-PARTICLE DENSITY MATRIX in expval.data file.')
if(trim(adjustl(char_80))=='START TWO-PARTICLE DENSITY MATRIX')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_twop_dm=trim(adjustl(char_80))
if(accumulation_method_twop_dm/='VMC'.and.accumulation_method_twop_dm/&
&='DMC')call errstop('READ_TWOP_DENSITY_MATRIX','Two-particle density &
&matrix accumulation method (VMC or DMC) in expval.data not recognized&
&.')
if(twop_density_mat)then
if(isvmc.and.accumulation_method_twop_dm=='DMC')call errstop('READ_TWO&
&P_DENSITY_MATRIX','DMC two-particle density matrix accumulation on to&
&p of VMC-accumulated two-particle density matrix in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_twop_dm==&
&'VMC') call errstop('READ_TWOP_DENSITY_MATRIX','VMC two-particle dens&
&ity matrix accumulation on top of DMC-accumulated two-particle densit&
&y matrix in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaik1
if(xyzzyaaik1>max_spin_pairs)call errstop('READ_TWOP_DENSITY_MATRIX','&
&More two-particle density matrix sets in expval.data than possible sp&
&in-pair types.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaail1
if(xyzzyaail1<2)call errstop('READ_TWOP_DENSITY_MATRIX','Need at least&
& two bins for accumulation of two-particle density matrix in expval.d&
&ata.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaim1
if(xyzzyaaim1<1)call errstop('READ_TWOP_DENSITY_MATRIX','Need at least&
& one random point for accumulation of two-particle density matrix in &
&expval.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaiv1
if(xyzzyaaiv1<=0.d0.or.xyzzyaaiv1>1.d0)call errstop('READ_TWOP_DENSITY&
&_MATRIX','Fraction of number of pairs must be greater than zero and l&
&ess than or equal to one in expval.data.')
allocate(xyzzyaaix1(xyzzyaail1,xyzzyaaik1),xyzzyaajb1(xyzzyaail1,xyzzy&
&aaik1),xyzzyaajc1(xyzzyaail1,xyzzyaaik1),xyzzyaaiw1(xyzzyaail1,xyzzya&
&aik1),xyzzyaaiz1(xyzzyaail1,xyzzyaaik1),xyzzyaaja1(xyzzyaail1,xyzzyaa&
&ik1),xyzzyaaiy1(xyzzyaail1,xyzzyaaik1),xyzzyaajd1(xyzzyaail1,xyzzyaai&
&k1),xyzzyaaje1(xyzzyaail1,xyzzyaaik1),xyzzyaain1(xyzzyaaik1),xyzzyaai&
&o1(2,max_spin_pairs,xyzzyaaik1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_TWOP_DENSITY_MATRIX','expval_twop_dm&
&')
xyzzyaaix1=0.d0
xyzzyaajb1=0.d0
xyzzyaajc1=0.d0
xyzzyaaiw1=0.d0
xyzzyaaiz1=0.d0
xyzzyaaja1=0.d0
xyzzyaaiy1=0.d0
xyzzyaajd1=0.d0
xyzzyaaje1=0.d0
xyzzyaain1=0
xyzzyaaio1=0
if(accumulation_method_twop_dm=='DMC')then
allocate(xyzzyaajh1(xyzzyaail1,xyzzyaaik1),xyzzyaajf1(xyzzyaail1,xyzzy&
&aaik1),xyzzyaajg1(xyzzyaail1,xyzzyaaik1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_TWOP_DENSITY_MATRIX','expval_twop_dm&
&_*_dmc')
xyzzyaajh1=0.d0
xyzzyaajf1=0.d0
xyzzyaajg1=0.d0
endif
loop_sets: do xyzzyaaab17=1,xyzzyaaik1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaab17)))then
call errstop2('READ_TWOP_DENSITY_MATRIX','Error reading START SET x li&
&ne for x = ',xyzzyaaab17)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaain1(xyzzyaaab17)
if(xyzzyaain1(xyzzyaaab17)<1.or.xyzzyaain1(xyzzyaaab17)>max_spin_pairs&
&)call errstop2('READ_TWOP_DENSITY_MATRIX','Problem with number of par&
&ticles-pair types in set ',xyzzyaaab17)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaio1(1:2,1:xyzzyaain1(xyzzyaaab17)&
&,xyzzyaaab17)
read(xyzzyaaaf1,*,err=1,end=1)
if(any(xyzzyaaio1(1:2,1:xyzzyaain1(xyzzyaaab17),xyzzyaaab17)<1).or.any&
&(xyzzyaaio1(1:2,1:xyzzyaain1(xyzzyaaab17),xyzzyaaab17)>nspin))call er&
&rstop2('READ_TWOP_DENSITY_MATRIX','Particle type out of range in part&
&icle-pair list for set ',xyzzyaaab17)
do xyzzyaaaa17=1,xyzzyaail1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaab17)))then
xyzzyaaix1(:,xyzzyaaab17)=0.d0
xyzzyaajb1(:,xyzzyaaab17)=0.d0
xyzzyaajc1(:,xyzzyaaab17)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaaix1(xyzzyaaaa17,xyzzyaaab17),xyzzyaa&
&jc1(xyzzyaaaa17,xyzzyaaab17),xyzzyaajb1(xyzzyaaaa17,xyzzyaaab17)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaabs1=.true.
enddo loop_sets
return
1 call errstop('READ_TWOP_DENSITY_MATRIX','Problem reading TWO-PARTICL&
&E DENSITY MATRIX section of expval.data file.')
end subroutine xyzzyaapv1
subroutine xyzzyaapw1
implicit none
integer xyzzyaaaa18,xyzzyaaab18
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_COND_FRACTION','Problem reading exp&
&val.data file.')
if(xyzzyaaae1<0)call errstop('READ_COND_FRACTION','Could not find STAR&
&T CONDENSATE FRACTION in expval.data file.')
if(trim(adjustl(char_80))=='START CONDENSATE FRACTION')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accum_method_cond_frac=trim(adjustl(char_80))
if(accum_method_cond_frac/='VMC'.and.accum_method_cond_frac/='DMC')cal&
&l errstop('READ_COND_FRACTION','Condensate fraction estimator accumul&
&ation method in expval.data not recognized, should be VMC or DMC.')
if(cond_fraction)then
if(isvmc.and.accum_method_cond_frac=='DMC')call errstop('READ_COND_FRA&
&CTION','DMC condensate fraction estimator accumulation on top of VMC-&
&accumulated data in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accum_method_cond_frac =='VMC&
&')call errstop('READ_COND_FRACTION','VMC condensate fraction estimato&
&r accumulation on top of DMC-accumulated data in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaajj1
if(xyzzyaajj1>max_spin_pairs)call errstop('READ_COND_FRACTION','More c&
&ondensate-fraction estimator sets in expval.data than possible spin-p&
&air types.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaajk1
if(xyzzyaajk1<2)call errstop('READ_COND_FRACTION','Need at least two b&
&ins for accumulation of condensate fraction estimator in expval.data.&
&')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaajl1
if(xyzzyaajl1<1)call errstop('READ_COND_FRACTION','Need at least one r&
&andom point for accumulation of condensate fraction estimator in expv&
&al.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaju1
if(xyzzyaaju1<=0.d0.or.xyzzyaaju1>1.d0)call errstop('READ_COND_FRACTIO&
&N','Fraction of number of pairs must be greater than zero and less th&
&an or equal to one in expval.data.')
allocate(xyzzyaajw1(xyzzyaajk1,xyzzyaajj1),xyzzyaaka1(xyzzyaajk1,xyzzy&
&aajj1),xyzzyaakb1(xyzzyaajk1,xyzzyaajj1),xyzzyaajv1(xyzzyaajk1,xyzzya&
&ajj1),xyzzyaajy1(xyzzyaajk1,xyzzyaajj1),xyzzyaajz1(xyzzyaajk1,xyzzyaa&
&jj1),xyzzyaajx1(xyzzyaajk1,xyzzyaajj1),xyzzyaakc1(xyzzyaajk1,xyzzyaaj&
&j1),xyzzyaakd1(xyzzyaajk1,xyzzyaajj1),xyzzyaajm1(xyzzyaajj1),xyzzyaaj&
&n1(2,max_spin_pairs,xyzzyaajj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_COND_FRACTION','expval_cond_frac')
xyzzyaajw1=0.d0
xyzzyaaka1=0.d0
xyzzyaakb1=0.d0
xyzzyaajv1=0.d0
xyzzyaajy1=0.d0
xyzzyaajz1=0.d0
xyzzyaajx1=0.d0
xyzzyaakc1=0.d0
xyzzyaakd1=0.d0
xyzzyaajm1=0
xyzzyaajn1=0
if(accum_method_cond_frac=='DMC')then
allocate(xyzzyaakg1(xyzzyaajk1,xyzzyaajj1),xyzzyaake1(xyzzyaajk1,xyzzy&
&aajj1),xyzzyaakf1(xyzzyaajk1,xyzzyaajj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_COND_FRACTION','expval_cond_frac_*_d&
&mc')
xyzzyaakg1=0.d0
xyzzyaake1=0.d0
xyzzyaakf1=0.d0
endif
loop_sets: do xyzzyaaab18=1,xyzzyaajj1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaab18)))then
call errstop2('READ_COND_FRACTION','Error reading START SET x line for&
& x = ',xyzzyaaab18)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaajm1(xyzzyaaab18)
if(xyzzyaajm1(xyzzyaaab18)<1.or.xyzzyaajm1(xyzzyaaab18)>max_spin_pairs&
&)call errstop2('READ_COND_FRACTION','Problem with number of particles&
&-pair types in set ',xyzzyaaab18)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaajn1(1:2,1:xyzzyaajm1(xyzzyaaab18)&
&,xyzzyaaab18)
read(xyzzyaaaf1,*,err=1,end=1)
if(any(xyzzyaajn1(1:2,1:xyzzyaajm1(xyzzyaaab18),xyzzyaaab18)<1).or.any&
&(xyzzyaajn1(1:2,1:xyzzyaajm1(xyzzyaaab18),xyzzyaaab18)>nspin))call er&
&rstop2('READ_COND_FRACTION','Particle type out of range in particle-p&
&air list for set ',xyzzyaaab18)
do xyzzyaaaa18=1,xyzzyaajk1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaab18)))then
xyzzyaajw1(:,xyzzyaaab18)=0.d0
xyzzyaaka1(:,xyzzyaaab18)=0.d0
xyzzyaakb1(:,xyzzyaaab18)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaajw1(xyzzyaaaa18,xyzzyaaab18),xyzzyaa&
&kb1(xyzzyaaaa18,xyzzyaaab18),xyzzyaaka1(xyzzyaaaa18,xyzzyaaab18)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaabu1=.true.
enddo loop_sets
return
1 call errstop('READ_COND_FRACTION','Problem reading CONDENSATE FRACTI&
&ON section of expval.data file.')
end subroutine xyzzyaapw1
subroutine xyzzyaapx1
implicit none
integer xyzzyaaaa19,xyzzyaaab19,xyzzyaaac19,xyzzyaaad19
character(80) char_80
real(dp) xyzzyaaae19
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_MOM_DEN','Problem reading expval.da&
&ta file.')
if(xyzzyaaae1<0)call errstop('READ_MOM_DEN','Could not find START MOME&
&NTUM DENSITY in expval.data file.')
if(trim(adjustl(char_80))=='START MOMENTUM DENSITY')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_mom_den=trim(adjustl(char_80))
if(accumulation_method_mom_den/='VMC'.and.accumulation_method_mom_den/&
&='DMC')call errstop('READ_MOM_DEN','Momentum density accumulation met&
&hod (VMC or DMC) in expval.data not recognized.')
if(mom_den)then
if(isvmc.and.accumulation_method_mom_den=='DMC')call errstop('READ_MOM&
&_DEN','DMC momentum density accumulation on top of VMC-accumulated mo&
&mentum density in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_mom_den==&
&'VMC')call errstop('READ_MOM_DEN','VMC momentum density accumulation &
&on top of DMC-accumulated momentum density in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaahs1
if(on_top_ii==0)then
if(xyzzyaahs1>nspin)call errstop('READ_MOM_DEN','More momentum density&
& sets in expval.data than defined types of particle.')
else
if(xyzzyaahs1>1)call errstop('READ_MOM_DEN','More than one momentum de&
&nsity sets in expval.data, but this is an ON_TOP_PAIR calculation.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaht1
if(xyzzyaaht1<1)call errstop('READ_MOM_DEN','Need at least one random &
&point for accumulation of momentum density in expval.data.')
allocate(xyzzyaahx1(xyzzyaahs1),xyzzyaaib1(expval_ngvec(1),xyzzyaahs1)&
&,xyzzyaaic1(expval_ngvec(1),xyzzyaahs1),xyzzyaahw1(xyzzyaahs1),xyzzya&
&ahz1(expval_ngvec(1),xyzzyaahs1),xyzzyaaia1(expval_ngvec(1),xyzzyaahs&
&1),xyzzyaahy1(xyzzyaahs1),xyzzyaaid1(expval_ngvec(1),xyzzyaahs1),xyzz&
&yaaie1(expval_ngvec(1),xyzzyaahs1),xyzzyaahu1(xyzzyaahs1),xyzzyaahv1(&
&nspin,xyzzyaahs1),xyzzyaaij1(3,expval_ngvec(1)),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOM_DEN','expval_mom_den')
xyzzyaahx1=0.d0
xyzzyaaib1=0.d0
xyzzyaaic1=0.d0
xyzzyaahw1=0.d0
xyzzyaahz1=0.d0
xyzzyaaia1=0.d0
xyzzyaahy1=0.d0
xyzzyaaid1=0.d0
xyzzyaaie1=0.d0
xyzzyaahu1=0
xyzzyaahv1=0
do xyzzyaaaa19=1,expval_ngvec(1)
xyzzyaaij1(:,xyzzyaaaa19)=expval_gvec(:,xyzzyaaaa19,1)-k_offset(:)
enddo
if(accumulation_method_mom_den=='DMC')then
allocate(xyzzyaaii1(xyzzyaahs1),xyzzyaaig1(expval_ngvec(1),xyzzyaahs1)&
&,xyzzyaaih1(expval_ngvec(1),xyzzyaahs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_MOM_DEN','expval_mom_den_*_dmc')
xyzzyaaii1=0.d0
xyzzyaaig1=0.d0
xyzzyaaih1=0.d0
endif
loop_sets: do xyzzyaaab19=1,xyzzyaahs1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaab19)))then
call errstop2('READ_MOM_DEN','Error reading START SET x line for x = '&
&,xyzzyaaab19)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='on-top pair')then
if(on_top_ii==0)call errstop('READ_MOM_DEN','Attempting to continue an&
& on-top momentum density accumulation, but ON_TOP_PAIR not set.')
xyzzyaahu1(xyzzyaaab19)=1
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaac19,xyzzyaaad19
if(xyzzyaaac19/=which_spin(on_top_ii).or.xyzzyaaad19/=which_spin(on_to&
&p_jj))then
if(xyzzyaaad19==which_spin(on_top_ii).and.xyzzyaaac19==which_spin(on_t&
&op_jj))call errstop('READ_MOM_DEN','Attempting to continue on-top mom&
&entum density accumulation, but ON_TOP_PAIR declaration has lines in &
&reverse order, which could cause issues.')
call errstop('READ_MOM_DEN','Attempting to continue on-top momentum de&
&nsity accumulation, but ON_TOP_PAIR specifies different particle type&
&s from original run.')
endif
xyzzyaahv1(1,xyzzyaaab19)=which_spin(on_top_ii)
else
read(char_80,*,err=1)xyzzyaahu1(xyzzyaaab19)
if(xyzzyaahu1(xyzzyaaab19)<1.or.xyzzyaahu1(xyzzyaaab19)>nspin)call err&
&stop2('READ_MOM_DEN','Problem with number of particle types in set ',&
&xyzzyaaab19)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaahv1(1:xyzzyaahu1(xyzzyaaab19),xyz&
&zyaaab19)
if(any(xyzzyaahv1(1:xyzzyaahu1(xyzzyaaab19),xyzzyaaab19)<1).or.any(xyz&
&zyaahv1(1:xyzzyaahu1(xyzzyaaab19),xyzzyaaab19)>nspin))call errstop2('&
&READ_MOM_DEN','Particle type out of range in particle list for set ',&
&xyzzyaaab19)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaahx1(xyzzyaaab19)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaaa19=1,expval_ngvec(1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaab19)))then
xyzzyaahx1(xyzzyaaab19)=0.d0
xyzzyaaib1(:,xyzzyaaab19)=0.d0
xyzzyaaic1(:,xyzzyaaab19)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaaae19,xyzzyaaic1(xyzzyaaaa19,xyzzyaaa&
&b19),xyzzyaaib1(xyzzyaaaa19,xyzzyaaab19)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaabw1=.true.
enddo loop_sets
return
1 call errstop('READ_MOM_DEN','Problem reading MOMENTUM DENSITY sectio&
&n of expval.data file.')
end subroutine xyzzyaapx1
subroutine xyzzyaapy1
implicit none
integer xyzzyaaaa20,xyzzyaaab20
character(80) char_80
xyzzyaaki1=3
xyzzyaakj1=2
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_FINITE_DENSITY','Problem reading ex&
&pval.data file.')
if(xyzzyaaae1<0)call errstop('READ_FINITE_DENSITY','Could not find STA&
&RT FINITE DENSITY in expval.data file.')
if(trim(adjustl(char_80))=='START FINITE DENSITY')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_fin_den=trim(adjustl(char_80))
if(accumulation_method_fin_den/='VMC'.and.accumulation_method_fin_den/&
&='DMC')call errstop('READ_FINITE_DENSITY','Finite density accumulatio&
&n method (VMC or DMC) in expval.data not recognized.')
if(finite_density)then
if(isvmc.and.accumulation_method_fin_den=='DMC')call errstop('READ_FIN&
&ITE_DENSITY','DMC finite density accumulation on top of VMC-accumulat&
&ed finite density in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_fin_den==&
&'VMC')call errstop('READ_FINITE_DENSITY','VMC finite density accumula&
&tion on top of DMC-accumulated finite density in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)fin_den_basis
select case(trim(adjustl(fin_den_basis)))
case('step')
xyzzyaakl1=xyzzyaakm1
case('sph_bessel')
xyzzyaakl1=xyzzyaakn1
call errstop('READ_FINITE_DENSITY','Full functionality is not availabl&
&e for a Bessel function basis.')
case('chebyshev')
xyzzyaakl1=xyzzyaako1
call errstop('READ_FINITE_DENSITY','Full functionality is not availabl&
&e for a Chebyshev polynomial basis.')
case('exponential')
xyzzyaakl1=xyzzyaakp1
call errstop('READ_FINITE_DENSITY','Full functionality is not availabl&
&e for an exponential function basis.')
case default
call errstop('READ_FINITE_DENSITY','Finite density basis function in e&
&xpval.data not recognized.')
end select
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaakk1
if(xyzzyaakk1<2)call errstop('READ_FINITE_DENSITY','Need at least two &
&bins for accumulation of finite density in expval.data.')
select case(xyzzyaakl1)
case(xyzzyaakm1,xyzzyaako1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaakq1
xyzzyaakr1=xyzzyaakq1
case(xyzzyaakp1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaakr1
xyzzyaakq1=xyzzyaakr1
case(xyzzyaakn1)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaakq1,xyzzyaakr1
end select
allocate(xyzzyaalg1(xyzzyaakk1,xyzzyaaki1),xyzzyaalh1(xyzzyaakk1,xyzzy&
&aaki1),xyzzyaale1(xyzzyaakk1,xyzzyaaki1),xyzzyaalf1(xyzzyaakk1,xyzzya&
&aki1),xyzzyaali1(xyzzyaakk1,xyzzyaaki1),xyzzyaalj1(xyzzyaakk1,xyzzyaa&
&ki1),xyzzyaalo1(xyzzyaakk1,xyzzyaakj1),xyzzyaalp1(xyzzyaakk1,xyzzyaak&
&j1),xyzzyaalm1(xyzzyaakk1,xyzzyaakj1),xyzzyaaln1(xyzzyaakk1,xyzzyaakj&
&1),xyzzyaalq1(xyzzyaakk1,xyzzyaakj1),xyzzyaalr1(xyzzyaakk1,xyzzyaakj1&
&),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_FINITE_DENSITY','expval_fin_den arra&
&ys')
xyzzyaalg1=0.d0
xyzzyaalh1=0.d0
xyzzyaale1=0.d0
xyzzyaalf1=0.d0
xyzzyaali1=0.d0
xyzzyaalj1=0.d0
xyzzyaalo1=0.d0
xyzzyaalp1=0.d0
xyzzyaalm1=0.d0
xyzzyaaln1=0.d0
xyzzyaalq1=0.d0
xyzzyaalr1=0.d0
xyzzyaaks1=0.d0
xyzzyaakt1=0.d0
xyzzyaaku1=0.d0
xyzzyaakw1=0.d0
xyzzyaakx1=0.d0
xyzzyaaky1=0.d0
xyzzyaala1=0.d0
xyzzyaalb1=0.d0
xyzzyaalc1=0.d0
if(accumulation_method_fin_den=='DMC')then
allocate(xyzzyaalk1(xyzzyaakk1,xyzzyaaki1),xyzzyaall1(xyzzyaakk1,xyzzy&
&aaki1),xyzzyaals1(xyzzyaakk1,xyzzyaakj1),xyzzyaalt1(xyzzyaakk1,xyzzya&
&akj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_FINITE_DENSITY','expval_fin_den DMC &
&arrays')
xyzzyaalk1=0.d0
xyzzyaall1=0.d0
xyzzyaals1=0.d0
xyzzyaalt1=0.d0
xyzzyaakv1=0.d0
xyzzyaakz1=0.d0
xyzzyaald1=0.d0
endif
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaakt1,xyzzyaakx1,xyzzyaalb1
loop_sets_ij: do xyzzyaaab20=1,xyzzyaaki1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START e-e SET '//trim(i2s(xyzzyaaab20)))th&
&en
if(trim(adjustl(char_80))=='END FINITE DENSITY')then
backspace(xyzzyaaaf1)
exit loop_sets_ij
endif
call errstop2('READ_FINITE_DENSITY','Error reading START e-e SET x lin&
&e for x = ',xyzzyaaab20)
endif
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END e-e SET '//trim(i2s(xyzzyaaab20)))cycl&
&e loop_sets_ij
do xyzzyaaaa20=1,xyzzyaakk1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
read(char_80,*,err=1,end=1)xyzzyaalg1(xyzzyaaaa20,xyzzyaaab20),xyzzyaa&
&lh1(xyzzyaaaa20,xyzzyaaab20)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
enddo loop_sets_ij
loop_sets_i: do xyzzyaaab20=1,xyzzyaakj1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START e-nucleus SET '//trim(i2s(xyzzyaaab2&
&0)))then
if(trim(adjustl(char_80))=='END FINITE DENSITY')then
backspace(xyzzyaaaf1)
exit loop_sets_i
endif
call errstop2('READ_FINITE_DENSITY','Error reading START e-nucleus SET&
& x line for x = ',xyzzyaaab20)
endif
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END e-nucleus SET '//trim(i2s(xyzzyaaab20)&
&))cycle loop_sets_i
do xyzzyaaaa20=1,xyzzyaakk1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
read(char_80,*,err=1,end=1)xyzzyaalo1(xyzzyaaaa20,xyzzyaaab20),xyzzyaa&
&lp1(xyzzyaaaa20,xyzzyaaab20)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
enddo loop_sets_i
return
1 call errstop('READ_FINITE_DENSITY','Problem reading FINITE DENSITY s&
&ection of expval.data file.')
end subroutine xyzzyaapy1
subroutine xyzzyaapz1
implicit none
integer xyzzyaaaa21,xyzzyaaab21,xyzzyaaac21,xyzzyaaad21
character(80) char_80
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_POPULATION','Problem reading expval&
&.data file.')
if(xyzzyaaae1<0)call errstop('READ_POPULATION','Could not find START P&
&OPULATION in expval.data file.')
if(trim(adjustl(char_80))=='START POPULATION')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_population=trim(adjustl(char_80))
if(accumulation_method_population/='VMC'.and.accumulation_method_popul&
&ation /='DMC')call errstop('READ_POPULATION','Population accumulation&
& method (VMC or DMC) in expval.data not recognized.')
if(population)then
if(isvmc.and.accumulation_method_population=='DMC')call errstop('READ_&
&POPULATION','DMC population accumulation on top of VMC-accumulated po&
&pulation in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_populatio&
&n =='VMC')call errstop('READ_POPULATION','VMC population accumulation&
& on top of DMC-accumulated population in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaamn1
if(xyzzyaamn1>nspin)call errstop('READ_POPULATION','More population se&
&ts in expval.data than defined types of particle.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaac21
if(xyzzyaaac21/=nitot)call errstop('READ_POPULATION','Number of ions i&
&n POPULATION block in expval.data differs from the number of ions in &
&xwfn.data.')
allocate(xyzzyaamr1(xyzzyaamn1),xyzzyaamv1(nitot,xyzzyaamn1),xyzzyaamw&
&1(nitot,xyzzyaamn1),xyzzyaamq1(xyzzyaamn1),xyzzyaamt1(nitot,xyzzyaamn&
&1),xyzzyaamu1(nitot,xyzzyaamn1),xyzzyaams1(xyzzyaamn1),xyzzyaamx1(nit&
&ot,xyzzyaamn1),xyzzyaamy1(nitot,xyzzyaamn1),xyzzyaamo1(xyzzyaamn1),xy&
&zzyaamp1(nspin,xyzzyaamn1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_POPULATION','expval_population')
xyzzyaamr1=0.d0
xyzzyaamv1=0.d0
xyzzyaamw1=0.d0
xyzzyaamq1=0.d0
xyzzyaamt1=0.d0
xyzzyaamu1=0.d0
xyzzyaams1=0.d0
xyzzyaamx1=0.d0
xyzzyaamy1=0.d0
xyzzyaamo1=0
xyzzyaamp1=0
if(accumulation_method_population=='DMC')then
allocate(xyzzyaanb1(xyzzyaamn1),xyzzyaamz1(nitot,xyzzyaamn1),xyzzyaana&
&1(nitot,xyzzyaamn1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_POPULATION','expval_population_*_dmc&
&')
xyzzyaanb1=0.d0
xyzzyaamz1=0.d0
xyzzyaana1=0.d0
endif
loop_sets: do xyzzyaaab21=1,xyzzyaamn1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaab21)))call e&
&rrstop2('READ_POPULATION','Error reading START SET x line for x = ',x&
&yzzyaaab21)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaamo1(xyzzyaaab21)
if(xyzzyaamo1(xyzzyaaab21)<1.or.xyzzyaamo1(xyzzyaaab21)>nspin) call er&
&rstop2('READ_POPULATION','Problem with number of particles types in s&
&et ',xyzzyaaab21)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaamp1(1:xyzzyaamo1(xyzzyaaab21),xyz&
&zyaaab21)
if(any(xyzzyaamp1(1:xyzzyaamo1(xyzzyaaab21),xyzzyaaab21)<1).or.any(xyz&
&zyaamp1(1:xyzzyaamo1(xyzzyaaab21),xyzzyaaab21)>nspin))call errstop2('&
&READ_POPULATION','Particle type out of range in particle list for set&
& ',xyzzyaaab21)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaamr1(xyzzyaaab21)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaaa21=1,nitot
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaab21)))then
xyzzyaamr1(xyzzyaaab21)=0.d0
xyzzyaamv1(:,xyzzyaaab21)=0.d0
xyzzyaamw1(:,xyzzyaaab21)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaaad21,xyzzyaamw1(xyzzyaaaa21,xyzzyaaa&
&b21),xyzzyaamv1(xyzzyaaaa21,xyzzyaaab21)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaaca1=.true.
enddo loop_sets
return
1 call errstop('READ_POPULATION','Problem reading POPULATION section o&
&f expval.data file.')
end subroutine xyzzyaapz1
subroutine xyzzyaaqa1
implicit none
integer xyzzyaaaa22,xyzzyaaab22,xyzzyaaac22,xyzzyaaad22
character(80) char_80
if(mol_spin_density)then
xyzzyaalv1=nspin
else
xyzzyaalv1=1
endif
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_MOL_DENSITY','Problem reading expva&
&l.data file.')
if(xyzzyaaae1<0)call errstop('READ_MOL_DENSITY','Could not find START &
&MOLECULAR DENSITY in expval.data file.')
if(trim(adjustl(char_80))=='START MOLECULAR DENSITY')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_mol_den=trim(adjustl(char_80))
if(accumulation_method_mol_den/='VMC'.and.accumulation_method_mol_den/&
&='DMC')call errstop('READ_MOL_DENSITY','Finite density accumulation m&
&ethod (VMC or DMC) in expval.data not recognized.')
if(mol_density)then
if(isvmc.and.accumulation_method_mol_den=='DMC')call errstop('READ_MOL&
&_DENSITY','DMC finite density accumulation on top of VMC-accumulated &
&finite density in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_mol_den==&
&'VMC')call errstop('READ_MOL_DENSITY','VMC mol density accumulation o&
&n top of DMC-accumulated mol density in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaalw1(1:3)
if(product(xyzzyaalw1)<2)call errstop('READ_MOL_DENSITY','Need at leas&
&t two bins for accumulation of finite density in expval.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaly1(1:3)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaalx1(1:3)
if(.not.all(xyzzyaalx1<=xyzzyaaly1))call errstop('READ_MOL_DENSITY','A&
&ll components of point B must be <= those of point A.')
allocate(xyzzyaami1(xyzzyaalw1(1),xyzzyaalw1(2),                     x&
&yzzyaalw1(3),xyzzyaalv1),xyzzyaamh1(xyzzyaalw1(1),xyzzyaalw1(2),     &
&          xyzzyaalw1(3),xyzzyaalv1),xyzzyaamj1(xyzzyaalw1(1),xyzzyaal&
&w1(2),                     xyzzyaalw1(3),xyzzyaalv1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_MOL_DENSITY','expval_mol_den arrays'&
&)
xyzzyaamh1=0.d0
xyzzyaami1=0.d0
xyzzyaamj1=0.d0
xyzzyaalz1=0.d0
xyzzyaama1=0.d0
xyzzyaamb1=0.d0
xyzzyaamd1=0.d0
xyzzyaame1=0.d0
xyzzyaamf1=0.d0
if(accumulation_method_mol_den=='DMC')then
allocate(xyzzyaamk1(xyzzyaalw1(1),xyzzyaalw1(2),                   xyz&
&zyaalw1(3),xyzzyaalv1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_MOL_DENSITY','expval_mol_den DMC arr&
&ays')
xyzzyaamk1=0.d0
xyzzyaamc1=0.d0
xyzzyaamg1=0.d0
endif
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaama1,xyzzyaame1
loop_sets_i: do xyzzyaaad22=1,xyzzyaalv1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaad22)))then
if(trim(adjustl(char_80))=='END MOLECULAR DENSITY')then
backspace(xyzzyaaaf1)
exit loop_sets_i
endif
call errstop2('READ_MOL_DENSITY','Error reading START SET x line for x&
& = ',xyzzyaaad22)
endif
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaad22)))cycle lo&
&op_sets_i
do xyzzyaaac22=1,xyzzyaalw1(3)
do xyzzyaaab22=1,xyzzyaalw1(2)
do xyzzyaaaa22=1,xyzzyaalw1(1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
read(char_80,*,err=1,end=1)xyzzyaami1(xyzzyaaaa22,xyzzyaaab22,xyzzyaaa&
&c22,xyzzyaaad22)
enddo
enddo
enddo
read(xyzzyaaaf1,*,err=1,end=1)
enddo loop_sets_i
return
1 call errstop('READ_MOL_DENSITY','Problem reading MOLECULAR DENSITY s&
&ection of expval.data file.')
end subroutine xyzzyaaqa1
subroutine xyzzyaaqb1
implicit none
integer xyzzyaaaa23,xyzzyaaab23
character(80) char_80
real(dp) xyzzyaaac23
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_TWOP_DM_MOM','Problem reading expva&
&l.data file.')
if(xyzzyaaae1<0)call errstop('READ_TWOP_DM_MOM','Could not find TWO-PA&
&RTICLE DENSITY MATRIX (MOMENTUM SPACE) in expval.data file.')
if(trim(adjustl(char_80))=='START TWO-PARTICLE DENSITY MATRIX (MOMENTU&
&M SPACE)')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accumulation_method_twop_dm_mom=trim(adjustl(char_80))
if(accumulation_method_twop_dm_mom/='VMC'.and.accumulation_method_twop&
&_dm_mom/='DMC')call errstop('READ_TWOP_DM_MOM','Momentum density accu&
&mulation method (VMC or DMC) in expval.data not recognized.')
if(twop_dm_mom)then
if(isvmc.and.accumulation_method_twop_dm_mom=='DMC')call errstop('READ&
&_TWOP_DM_MOM','DMC momentum density accumulation on top of VMC-accumu&
&lated momentum density in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accumulation_method_twop_dm_m&
&om=='VMC')call errstop('READ_TWOP_DM_MOM','VMC momentum density accum&
&ulation on top of DMC-accumulated momentum density in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaanc1
if(xyzzyaanc1>max_spin_pairs)call errstop('READ_TWOP_DM_MOM','More mom&
&entum density sets in expval.data than defined types of particle.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaand1
if(xyzzyaand1<1)call errstop('READ_TWOP_DM_MOM','Need at least one ran&
&dom point for accumulation of two-particle density matrix (momentum s&
&pace) in expval.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaanh1
if(xyzzyaanh1<=0.d0.or.xyzzyaanh1>1.d0)call errstop('READ_TWOP_DM_MOM'&
&,'Fraction of number of pairs must be greater than zero and less than&
& or equal to one in expval.data.')
allocate(xyzzyaanj1(xyzzyaanc1),xyzzyaann1(expval_ngvec(1),xyzzyaanc1)&
&,xyzzyaano1(expval_ngvec(1),xyzzyaanc1),xyzzyaani1(xyzzyaanc1),xyzzya&
&anl1(expval_ngvec(1),xyzzyaanc1),xyzzyaanm1(expval_ngvec(1),xyzzyaanc&
&1),xyzzyaank1(xyzzyaanc1),xyzzyaanp1(expval_ngvec(1),xyzzyaanc1),xyzz&
&yaanq1(expval_ngvec(1),xyzzyaanc1),xyzzyaane1(xyzzyaanc1),xyzzyaanf1(&
&2,max_spin_pairs,xyzzyaanc1),xyzzyaanu1(3,expval_ngvec(1)),stat=xyzzy&
&aaad1)
call check_alloc(xyzzyaaad1,'READ_TWOP_DM_MOM','expval_twop_dm_mom')
xyzzyaanj1=0.d0
xyzzyaann1=0.d0
xyzzyaano1=0.d0
xyzzyaani1=0.d0
xyzzyaanl1=0.d0
xyzzyaanm1=0.d0
xyzzyaank1=0.d0
xyzzyaanp1=0.d0
xyzzyaanq1=0.d0
xyzzyaane1=0
xyzzyaanf1=0
do xyzzyaaaa23=1,expval_ngvec(1)
xyzzyaanu1(:,xyzzyaaaa23)=expval_gvec(:,xyzzyaaaa23,1)-k_offset(:)
enddo
if(accumulation_method_twop_dm_mom=='DMC')then
allocate(xyzzyaant1(xyzzyaanc1),xyzzyaanr1(expval_ngvec(1),xyzzyaanc1)&
&,xyzzyaans1(expval_ngvec(1),xyzzyaanc1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_TWOP_DM_MOM','expval_twop_dm_mom_*_d&
&mc')
xyzzyaant1=0.d0
xyzzyaanr1=0.d0
xyzzyaans1=0.d0
endif
loop_sets: do xyzzyaaab23=1,xyzzyaanc1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaab23)))then
call errstop2('READ_TWOP_DM_MOM','Error reading START SET x line for x&
& = ',xyzzyaaab23)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaane1(xyzzyaaab23)
if(xyzzyaane1(xyzzyaaab23)<1.or.xyzzyaane1(xyzzyaaab23)>max_spin_pairs&
&)call errstop2('READ_TWOP_DM_MOM','Problem with number of particles-p&
&air types in set ',xyzzyaaab23)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaanf1(1:2,1:xyzzyaane1(xyzzyaaab23)&
&,xyzzyaaab23)
if(any(xyzzyaanf1(1:2,1:xyzzyaane1(xyzzyaaab23),xyzzyaaab23)<1).or.any&
&(xyzzyaanf1(1:2,1:xyzzyaane1(xyzzyaaab23),xyzzyaaab23)>nspin))call er&
&rstop2('READ_TWOP_DM_MOM','Particle type out of range in particle-pai&
&r list for set ',xyzzyaaab23)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaanj1(xyzzyaaab23)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaaa23=1,expval_ngvec(1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaab23)))then
xyzzyaanj1(xyzzyaaab23)=0.d0
xyzzyaann1(:,xyzzyaaab23)=0.d0
xyzzyaano1(:,xyzzyaaab23)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaaac23,xyzzyaano1(xyzzyaaaa23,xyzzyaaa&
&b23),xyzzyaann1(xyzzyaaaa23,xyzzyaaab23)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaace1=.true.
enddo loop_sets
return
1 call errstop('READ_TWOP_DM_MOM','Problem reading TWO-PARTICLE DENSIT&
&Y MATRIX (MOMENTUM SPACE) section of expval.data file.')
end subroutine xyzzyaaqb1
subroutine xyzzyaaqc1
implicit none
integer xyzzyaaaa24,xyzzyaaab24
character(80) char_80
real(dp) xyzzyaaac24
do
read(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)char_80
if(xyzzyaaae1>0)call errstop('READ_COND_FRAC_MOM','Problem reading exp&
&val.data file.')
if(xyzzyaaae1<0)call errstop('READ_COND_FRAC_MOM','Could not find STAR&
&T CONDENSATE FRACTION in expval.data file.')
if(trim(adjustl(char_80))=='START CONDENSATE FRACTION (MOMENTUM SPACE)&
&')exit
enddo
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
accum_method_cond_frac_mom=trim(adjustl(char_80))
if(accum_method_cond_frac_mom/='VMC'.and.accum_method_cond_frac_mom/='&
&DMC')call errstop('READ_COND_FRAC_MOM','Condensate fraction estimator&
& accumulation method in expval.data not recognized, should be VMC or &
&DMC.')
if(cond_frac_mom)then
if(isvmc.and.accum_method_cond_frac_mom=='DMC')call errstop('READ_COND&
&_FRAC_MOM','DMC condensate fraction estimator accumulation on top of &
&VMC-accumulated data in expval.data.')
if((isdmc.or.isvmc_dmc.or.isdmc_dmc).and.accum_method_cond_frac_mom ==&
&'VMC')call errstop('READ_COND_FRAC_MOM','VMC condensate fraction esti&
&mator accumulation on top of DMC-accumulated data in expval.data.')
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaanw1
if(xyzzyaanw1>max_spin_pairs)call errstop('READ_COND_FRAC_MOM','More c&
&ondensate-fraction estimator sets in expval.data than possible spin-p&
&air types.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaanx1
if(xyzzyaanx1<1)call errstop('READ_COND_FRAC_MOM','Need at least one r&
&andom point for accumulation of condensate fraction estimator in expv&
&al.data.')
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaob1
if(xyzzyaaob1<=0.d0.or.xyzzyaaob1>1.d0)call errstop('READ_COND_FRAC_MO&
&M','Fraction of number of pairs must be greater than zero and less th&
&an or equal to one in expval.data.')
allocate(xyzzyaaod1(xyzzyaanw1),xyzzyaaoh1(expval_ngvec(1),xyzzyaanw1)&
&,xyzzyaaoi1(expval_ngvec(1),xyzzyaanw1),xyzzyaaoc1(xyzzyaanw1),xyzzya&
&aof1(expval_ngvec(1),xyzzyaanw1),xyzzyaaog1(expval_ngvec(1),xyzzyaanw&
&1),xyzzyaaoe1(xyzzyaanw1),xyzzyaaoj1(expval_ngvec(1),xyzzyaanw1),xyzz&
&yaaok1(expval_ngvec(1),xyzzyaanw1),xyzzyaany1(xyzzyaanw1),xyzzyaanz1(&
&2,max_spin_pairs,xyzzyaanw1),xyzzyaaoo1(3,expval_ngvec(1)),stat=xyzzy&
&aaad1)
call check_alloc(xyzzyaaad1,'READ_COND_FRAC_MOM','expval_cond_frac_mom&
&')
xyzzyaaod1=0.d0
xyzzyaaoh1=0.d0
xyzzyaaoi1=0.d0
xyzzyaaoc1=0.d0
xyzzyaaof1=0.d0
xyzzyaaog1=0.d0
xyzzyaaoe1=0.d0
xyzzyaaoj1=0.d0
xyzzyaaok1=0.d0
xyzzyaany1=0
xyzzyaanz1=0
do xyzzyaaaa24=1,expval_ngvec(1)
xyzzyaaoo1(:,xyzzyaaaa24)=expval_gvec(:,xyzzyaaaa24,1)-k_offset(:)
enddo
if(accum_method_cond_frac_mom=='DMC')then
allocate(xyzzyaaon1(xyzzyaanw1),xyzzyaaol1(expval_ngvec(1),xyzzyaanw1)&
&,xyzzyaaom1(expval_ngvec(1),xyzzyaanw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_COND_FRAC_MOM','expval_cond_frac_mom&
&_*_dmc')
xyzzyaaon1=0.d0
xyzzyaaol1=0.d0
xyzzyaaom1=0.d0
endif
loop_sets: do xyzzyaaab24=1,xyzzyaanw1
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))/='START SET '//trim(i2s(xyzzyaaab24)))then
call errstop2('READ_COND_FRAC_MOM','Error reading START SET x line for&
& x = ',xyzzyaaab24)
endif
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaany1(xyzzyaaab24)
if(xyzzyaany1(xyzzyaaab24)<1.or.xyzzyaany1(xyzzyaaab24)>max_spin_pairs&
&)call errstop2('READ_COND_FRAC_MOM','Problem with number of particles&
&-pair types in set ',xyzzyaaab24)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaanz1(1:2,1:xyzzyaany1(xyzzyaaab24)&
&,xyzzyaaab24)
if(any(xyzzyaanz1(1:2,1:xyzzyaany1(xyzzyaaab24),xyzzyaaab24)<1).or.any&
&(xyzzyaanz1(1:2,1:xyzzyaany1(xyzzyaaab24),xyzzyaaab24)>nspin))call er&
&rstop2('READ_COND_FRAC_MOM','Particle type out of range in particle-p&
&air list for set ',xyzzyaaab24)
read(xyzzyaaaf1,*,err=1,end=1)
read(xyzzyaaaf1,*,err=1,end=1)xyzzyaaod1(xyzzyaaab24)
read(xyzzyaaaf1,*,err=1,end=1)
do xyzzyaaaa24=1,expval_ngvec(1)
read(xyzzyaaaf1,'(a)',err=1,end=1)char_80
if(trim(adjustl(char_80))=='END SET '//trim(i2s(xyzzyaaab24)))then
xyzzyaaod1(xyzzyaaab24)=0.d0
xyzzyaaoh1(:,xyzzyaaab24)=0.d0
xyzzyaaoi1(:,xyzzyaaab24)=0.d0
cycle loop_sets
endif
read(char_80,*,err=1,end=1)xyzzyaaac24,xyzzyaaoi1(xyzzyaaaa24,xyzzyaaa&
&b24),xyzzyaaoh1(xyzzyaaaa24,xyzzyaaab24)
enddo
read(xyzzyaaaf1,*,err=1,end=1)
xyzzyaacg1=.true.
enddo loop_sets
return
1 call errstop('READ_COND_FRAC_MOM','Problem reading CONDENSATE FRACTI&
&ON section of expval.data file.')
end subroutine xyzzyaaqc1
subroutine write_expval
implicit none
integer xyzzyaaaa25
character(80) tmpr,tmpr2,tmpr3
if(am_master)then
open(xyzzyaaaf1,file='expval.data',status='unknown',iostat=xyzzyaaae1)
if(xyzzyaaae1/=0)call errstop('WRITE_EXPVAL','Error opening expval.dat&
&a file for writing.')
write(xyzzyaaaf1,'(a)',iostat=xyzzyaaae1)'START HEADER'
if(xyzzyaaae1/=0)call errstop('WRITE_EXPVAL','Error writing to expval.&
&data file.')
do xyzzyaaah1=1,xyzzyaaag1
write(xyzzyaaaf1,'(1x,a)')trim(xyzzyaach1(xyzzyaaah1))
enddo
write(xyzzyaaaf1,'(a)')'END HEADER'
write(xyzzyaaaf1,*)
write(xyzzyaaaf1,'(a)')'START EXPVAL'
write(xyzzyaaaf1,'(a)')'Title'
write(xyzzyaaaf1,'(1x,a)')trim(adjustl(title))
write(xyzzyaaaf1,'(a)')'File version'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaack1))
write(xyzzyaaaf1,'(a)')'Number of particle types (e.g. 2=electrons, 4=&
&electrons+holes)'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaacl1))
write(xyzzyaaaf1,'(a)')'Number of each type of particle'
write(xyzzyaaaf1,'(20(1x,a))')(trim(i2s(xyzzyaacn1(xyzzyaaaa25))),xyzz&
&yaaaa25=1,xyzzyaacl1)
write(xyzzyaaaf1,'(a)')'Dimensionality'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaacm1))
write(xyzzyaaaf1,'(a)')'Periodicity'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaci1))
write(xyzzyaaaf1,'(a)')'Primitive translation vectors (au)'
select case(xyzzyaaci1)
case(0)
write(xyzzyaaaf1,'(a)')' Finite system'
write(xyzzyaaaf1,'(a)')'Size of simulation cell'
write(xyzzyaaaf1,'(a)')' Finite system'
case(1)
tmpr=r2s2(pa1(1),'(es24.16)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr)
write(xyzzyaaaf1,'(a)')'Supercell matrix'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaacj1(1,1)))
write(xyzzyaaaf1,'(a)')'Length of simulation cell'
tmpr=r2s(xyzzyaacq1,'(es24.16)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
write(xyzzyaaaf1,'(a)')'"Radius" of simulation Wigner-Seitz cell'
tmpr=r2s(wigner_seitz_radius,'(es24.16)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
case(2)
tmpr=r2s2(pa1(1),'(es24.16)')
tmpr2=r2s2(pa1(2),'(es24.16)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2)
tmpr=r2s2(pa2(1),'(es24.16)')
tmpr2=r2s2(pa2(2),'(es24.16)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2)
write(xyzzyaaaf1,'(a)')'Supercell matrix (11),(22),(12),(21)'
write(xyzzyaaaf1,'(4(1x,a))')trim(i2s(xyzzyaacj1(1,1))),trim(i2s(xyzzy&
&aacj1(2,2))),trim(i2s(xyzzyaacj1(1,2))),trim(i2s(xyzzyaacj1(2,1)))
write(xyzzyaaaf1,'(a)')'Area of simulation cell'
tmpr=r2s(xyzzyaacp1,'(es24.16)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
write(xyzzyaaaf1,'(a)')'Radius of circle inscribed in Wigner-Seitz cel&
&l of simulation cell'
tmpr=r2s(wigner_seitz_radius,'(es24.16)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
case(3)
tmpr=r2s2(pa1(1),'(es24.16)')
tmpr2=r2s2(pa1(2),'(es24.16)')
tmpr3=r2s2(pa1(3),'(es24.16)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
tmpr=r2s2(pa2(1),'(es24.16)')
tmpr2=r2s2(pa2(2),'(es24.16)')
tmpr3=r2s2(pa2(3),'(es24.16)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
tmpr=r2s2(pa3(1),'(es24.16)')
tmpr2=r2s2(pa3(2),'(es24.16)')
tmpr3=r2s2(pa3(3),'(es24.16)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
write(xyzzyaaaf1,'(a)')'Supercell matrix (11),(22),(33),(12),(13),(21)&
&,(23),(31),(32)'
write(xyzzyaaaf1,'(9(1x,a))')trim(i2s(xyzzyaacj1(1,1))),trim(i2s(xyzzy&
&aacj1(2,2))),trim(i2s(xyzzyaacj1(3,3))),trim(i2s(xyzzyaacj1(1,2))),tr&
&im(i2s(xyzzyaacj1(1,3))),trim(i2s(xyzzyaacj1(2,1))),trim(i2s(xyzzyaac&
&j1(2,3))),trim(i2s(xyzzyaacj1(3,1))),trim(i2s(xyzzyaacj1(3,2)))
write(xyzzyaaaf1,'(a)')'Volume of simulation cell'
tmpr=r2s(xyzzyaacr1,'(es24.16)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
write(xyzzyaaaf1,'(a)')'Radius of sphere inscribed in Wigner-Seitz cel&
&l of simulation cell'
tmpr=r2s(wigner_seitz_radius,'(es24.16)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
case default
call errstop('WRITE_EXPVAL','Invalid periodicity.')
end select
write(xyzzyaaaf1,'(a)')'Number of available G-vector sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaacw1))
write(xyzzyaaaf1,*)
do xyzzyaaaa25=1,xyzzyaacw1
call xyzzyaaqd1(xyzzyaaaa25)
write(xyzzyaaaf1,*)
enddo
endif
if(xyzzyaaax1(1))call write_density
if(xyzzyaaax1(2))call xyzzyaaqe1
if(xyzzyaaax1(3))call xyzzyaaqf1
if(xyzzyaaax1(4))call xyzzyaaqg1
if(xyzzyaaax1(5))call xyzzyaaqh1
if(xyzzyaaax1(6))call xyzzyaaqi1
if(xyzzyaaax1(7))call xyzzyaaqj1
if(xyzzyaaax1(8))call xyzzyaaqk1
if(xyzzyaaax1(9))call xyzzyaaql1
if(xyzzyaaax1(10))call xyzzyaaqm1
if(xyzzyaaax1(11))call xyzzyaaqn1
if(xyzzyaaax1(12))call xyzzyaaqp1
if(xyzzyaaax1(13))call xyzzyaaqr1
if(xyzzyaaax1(14))call xyzzyaaqs1
if(xyzzyaaax1(15))call xyzzyaaqt1
if(xyzzyaaax1(16))call xyzzyaaqq1
if(xyzzyaaax1(17))call xyzzyaaqo1
if(am_master)then
write(xyzzyaaaf1,'(a)')'END EXPVAL'
close(xyzzyaaaf1)
endif
end subroutine write_expval
subroutine xyzzyaaqd1(set)
implicit none
integer,intent(in) :: set
integer xyzzyaaaa26
character(80) tmpr,tmpr2,tmpr3
write(xyzzyaaaf1,'(a)')'START GVECTOR SET '//trim(i2s(set))
write(xyzzyaaaf1,'(a)')'Energy cutoff (au) used to generate set'
tmpr=r2s(xyzzyaadb1(set),'(es24.16)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
write(xyzzyaaaf1,'(a)')'Number of G-vectors in set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(expval_ngvec(set)))
if(abs(xyzzyaadd1(1,1,set)-pb1(1))<1.d-8.and.abs(xyzzyaadd1(2,1,set)-p&
&b1(2))<1.d-8.and.abs(xyzzyaadd1(3,1,set)-pb1(3))<1.d-8)then
write(xyzzyaaaf1,'(a)')'Primitive reciprocal lattice vectors (au)'
else
write(xyzzyaaaf1,'(a)')'Supercell reciprocal lattice vectors (au)'
endif
do xyzzyaaaa26=1,3
tmpr=r2s2(xyzzyaadd1(1,xyzzyaaaa26,set),'(es24.16)')
tmpr2=r2s2(xyzzyaadd1(2,xyzzyaaaa26,set),'(es24.16)')
tmpr3=r2s2(xyzzyaadd1(3,xyzzyaaaa26,set),'(es24.16)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
enddo
write(xyzzyaaaf1,'(a)')'G-vector components Gx, Gy, Gz (au)'
do xyzzyaaaa26=1,expval_ngvec(set)
tmpr=r2s2(expval_gvec(1,xyzzyaaaa26,set),'(es24.16)')
tmpr2=r2s2(expval_gvec(2,xyzzyaaaa26,set),'(es24.16)')
tmpr3=r2s2(expval_gvec(3,xyzzyaaaa26,set),'(es24.16)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
enddo
write(xyzzyaaaf1,'(a)')'END GVECTOR SET '//trim(i2s(set))
end subroutine xyzzyaaqd1
subroutine write_density
implicit none
integer xyzzyaaaa27,xyzzyaaab27,xyzzyaaac27,xyzzyaaad27,xyzzyaaae27
real(dp) xyzzyaaaf27,xyzzyaaag27,xyzzyaaah27,xyzzyaaai27,xyzzyaaaj27,x&
&yzzyaaak27,xyzzyaaal27,xyzzyaaam27,xyzzyaaan27
real(dp) :: xyzzyaaao27=1.d-6
complex(dp) xyzzyaaap27
if(density)then
xyzzyaadm1=czero
call mpi_reduce(xyzzyaadk1,xyzzyaadm1,expval_ngvec(xyzzyaadg1)*den_nse&
&ts,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_den in write_density.')
endif
if(am_master)then
if(density)then
do xyzzyaaad27=1,den_nsets
xyzzyaadi1(xyzzyaaad27)=real(xyzzyaadm1(1,xyzzyaaad27),dp)/real(xyzzya&
&aco1(xyzzyaaad27),dp)+xyzzyaadj1(xyzzyaaad27)
if(xyzzyaadi1(xyzzyaaad27)==0.d0)call errstop('WRITE_DENSITY','Zero we&
&ight.')
xyzzyaaaf27=1.d0/(xyzzyaadi1(xyzzyaaad27)*real(npcells,dp))
xyzzyaaah27=xyzzyaadj1(xyzzyaaad27)*real(npcells,dp)
do xyzzyaaaa27=1,expval_ngvec(xyzzyaadg1)
xyzzyaaap27=xyzzyaaah27*xyzzyaadl1(xyzzyaaaa27,xyzzyaaad27)+xyzzyaadm1&
&(xyzzyaaaa27,xyzzyaaad27)
xyzzyaadm1(xyzzyaaaa27,xyzzyaaad27)=xyzzyaaaf27*xyzzyaaap27
enddo
enddo
else
xyzzyaadi1(:)=xyzzyaadj1(:)
endif
if(permit_den_symm)then
if(.not.allocated(den_sc))call errstop('WRITE_DENSITY','NO SCF density&
& to perform QMC density symmetrization. Please use RUNTYPE=gen_mpc to&
& generate an mpc.data file first.')
xyzzyaaae27=size(den_sc,1)
do xyzzyaaad27=1,den_nsets
xyzzyaaaa27=1
do
xyzzyaaap27=czero
xyzzyaaac27=0
xyzzyaaam27=den_gvec(1,xyzzyaaaa27)*den_gvec(1,xyzzyaaaa27)+den_gvec(2&
&,xyzzyaaaa27)*den_gvec(2,xyzzyaaaa27)+den_gvec(3,xyzzyaaaa27)*den_gve&
&c(3,xyzzyaaaa27)
do
xyzzyaaaf27=abs(real(xyzzyaadm1(xyzzyaaaa27+xyzzyaaac27,xyzzyaaad27),d&
&p))
xyzzyaaag27=abs(aimag(xyzzyaadm1(xyzzyaaaa27+xyzzyaaac27,xyzzyaaad27))&
&)
xyzzyaaap27=xyzzyaaap27+cmplx(xyzzyaaaf27,xyzzyaaag27,dp)
xyzzyaaac27=xyzzyaaac27+1
if(xyzzyaaaa27+xyzzyaaac27>xyzzyaaae27.or.xyzzyaaaa27+xyzzyaaac27>expv&
&al_ngvec(xyzzyaaad27))exit
xyzzyaaan27=den_gvec(1,xyzzyaaaa27+xyzzyaaac27)*den_gvec(1,xyzzyaaaa27&
&+xyzzyaaac27)+den_gvec(2,xyzzyaaaa27+xyzzyaaac27)*den_gvec(2,xyzzyaaa&
&a27+xyzzyaaac27)+den_gvec(3,xyzzyaaaa27+xyzzyaaac27)*den_gvec(3,xyzzy&
&aaaa27+xyzzyaaac27)
if(.not.approx_equal(xyzzyaaam27,xyzzyaaan27,xyzzyaaao27))exit
xyzzyaaai27=dble(den_sc(xyzzyaaaa27+xyzzyaaac27,xyzzyaaad27))
xyzzyaaaj27=dble(den_sc(xyzzyaaaa27,xyzzyaaad27))
xyzzyaaak27=aimag(den_sc(xyzzyaaaa27+xyzzyaaac27,xyzzyaaad27))
xyzzyaaal27=aimag(den_sc(xyzzyaaaa27,xyzzyaaad27))
if(.not.approx_equal(abs(xyzzyaaai27),abs(xyzzyaaaj27),xyzzyaaao27).or&
&..not.approx_equal(abs(xyzzyaaak27),abs(xyzzyaaal27),xyzzyaaao27))exi&
&t
enddo
if(xyzzyaaac27>1)xyzzyaaap27=xyzzyaaap27/real(xyzzyaaac27,dp)
do xyzzyaaab27=xyzzyaaaa27,xyzzyaaaa27+xyzzyaaac27-1
xyzzyaaaf27=sign(1.d0,real(xyzzyaadm1(xyzzyaaab27,xyzzyaaad27),dp))*re&
&al(xyzzyaaap27,dp)
xyzzyaaag27=sign(1.d0,aimag(xyzzyaadm1(xyzzyaaab27,xyzzyaaad27)))*aima&
&g(xyzzyaaap27)
xyzzyaadm1(xyzzyaaab27,xyzzyaaad27)=cmplx(xyzzyaaaf27,xyzzyaaag27,dp)
enddo
xyzzyaaaa27=xyzzyaaaa27+xyzzyaaac27
if(xyzzyaaaa27>xyzzyaaae27.or.xyzzyaaaa27+xyzzyaaac27>expval_ngvec(xyz&
&zyaaad27))exit
enddo
enddo
endif
write(xyzzyaaaf1,'(a)')'START DENSITY'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_den
write(xyzzyaaaf1,'(a)')'Use G-vector set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaadg1))
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(den_nsets))
do xyzzyaaad27=1,den_nsets
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaad27))
write(xyzzyaaaf1,'(a)')'Particle type'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaadh1(xyzzyaaad27)))
write(xyzzyaaaf1,'(a)')'Total weight'
write(xyzzyaaaf1,*)xyzzyaadi1(xyzzyaaad27)
write(xyzzyaaaf1,'(a)')'Complex charge-density coefficients (real part&
&, imaginary part)'
if(density)then
if(expval_complex_den)then
do xyzzyaaaa27=1,expval_ngvec(xyzzyaadg1)
xyzzyaaap27=xyzzyaadm1(xyzzyaaaa27,xyzzyaaad27)
write(xyzzyaaaf1,*)real(xyzzyaaap27,dp),aimag(xyzzyaaap27)
enddo
else
do xyzzyaaaa27=1,expval_ngvec(xyzzyaadg1)
xyzzyaaaf27=real(xyzzyaadm1(xyzzyaaaa27,xyzzyaaad27),dp)
write(xyzzyaaaf1,*)xyzzyaaaf27,' 0.'
enddo
endif
else
do xyzzyaaaa27=1,expval_ngvec(xyzzyaadg1)
xyzzyaaap27=xyzzyaadl1(xyzzyaaaa27,xyzzyaaad27)
write(xyzzyaaaf1,*)real(xyzzyaaap27,dp),aimag(xyzzyaaap27)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaad27))
enddo
write(xyzzyaaaf1,'(a)')'END DENSITY'
write(xyzzyaaaf1,*)
endif
end subroutine write_density
subroutine xyzzyaaqe1
implicit none
integer xyzzyaaaa28,xyzzyaaab28
real(dp) xyzzyaaac28,xyzzyaaad28
complex(dp) xyzzyaaae28
if(spin_density)then
xyzzyaadv1=czero
call mpi_reduce(xyzzyaadt1,xyzzyaadv1,expval_ngvec(xyzzyaado1)*xyzzyaa&
&dp1,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_sden in write_spin_density.&
&')
endif
if(am_master)then
if(spin_density)then
do xyzzyaaab28=1,xyzzyaadp1
if(nele(xyzzyaaab28)/=0)then
xyzzyaadr1(xyzzyaaab28)=real(xyzzyaadv1(1,xyzzyaaab28),dp)/real(nele(x&
&yzzyaaab28),dp)+xyzzyaads1(xyzzyaaab28)
if(xyzzyaadr1(xyzzyaaab28)==0.d0)call errstop('WRITE_SPIN_DENSITY','Ze&
&ro weight.')
xyzzyaaac28=1.d0/(xyzzyaadr1(xyzzyaaab28)*real(npcells,dp))
xyzzyaaad28=xyzzyaads1(xyzzyaaab28)*real(npcells,dp)
do xyzzyaaaa28=1,expval_ngvec(xyzzyaado1)
xyzzyaaae28=xyzzyaaad28*xyzzyaadu1(xyzzyaaaa28,xyzzyaaab28)+xyzzyaadv1&
&(xyzzyaaaa28,xyzzyaaab28)
xyzzyaadv1(xyzzyaaaa28,xyzzyaaab28)=xyzzyaaac28*xyzzyaaae28
enddo
endif
enddo
else
if(.not.pair_corr)xyzzyaadr1(:)=xyzzyaads1(:)
endif
write(xyzzyaaaf1,'(a)')'START SPIN DENSITY'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_sden
write(xyzzyaaaf1,'(a)')'Use G-vector set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaado1))
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaadp1))
do xyzzyaaab28=1,xyzzyaadp1
if(nele(xyzzyaaab28)/=0)then
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab28))
write(xyzzyaaaf1,'(a)')'Particle type'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaadq1(xyzzyaaab28)))
write(xyzzyaaaf1,'(a)')'Total weight'
write(xyzzyaaaf1,*)xyzzyaadr1(xyzzyaaab28)
write(xyzzyaaaf1,'(a)')'Complex charge-density coefficients (real part&
&, imaginary part)'
if(spin_density)then
if(xyzzyaadx1)then
do xyzzyaaaa28=1,expval_ngvec(xyzzyaado1)
xyzzyaaae28=xyzzyaadv1(xyzzyaaaa28,xyzzyaaab28)
write(xyzzyaaaf1,*)real(xyzzyaaae28,dp),aimag(xyzzyaaae28)
enddo
else
do xyzzyaaaa28=1,expval_ngvec(xyzzyaado1)
xyzzyaaac28=real(xyzzyaadv1(xyzzyaaaa28,xyzzyaaab28),dp)
write(xyzzyaaaf1,*)xyzzyaaac28,' 0.'
enddo
endif
else
do xyzzyaaaa28=1,expval_ngvec(xyzzyaado1)
xyzzyaaae28=xyzzyaadu1(xyzzyaaaa28,xyzzyaaab28)
write(xyzzyaaaf1,*)real(xyzzyaaae28,dp),aimag(xyzzyaaae28)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab28))
else
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab28))
write(xyzzyaaaf1,'(1x,a)')'No particles with this spin.'
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab28))
endif
enddo
write(xyzzyaaaf1,'(a)')'END SPIN DENSITY'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqe1
subroutine xyzzyaaqf1
implicit none
integer xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29
real(dp) xyzzyaaad29,xyzzyaaae29
complex(dp) xyzzyaaaf29
if(spin_density_mat)then
xyzzyaaef1=czero
call mpi_reduce(xyzzyaaed1,xyzzyaaef1,4*expval_ngvec(xyzzyaady1)*xyzzy&
&aadz1,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_sdenm in write_spin_density&
&_matrix.')
endif
if(am_master)then
if(spin_density_mat)then
do xyzzyaaac29=1,xyzzyaadz1
if(nele(xyzzyaaac29)/=0)then
xyzzyaaeb1(xyzzyaaac29)=real(xyzzyaaef1(1,1,xyzzyaaac29)+xyzzyaaef1(1,&
&4,xyzzyaaac29),dp)/real(nele(xyzzyaaac29),dp)+xyzzyaaec1(xyzzyaaac29)
if(xyzzyaaeb1(xyzzyaaac29)==0.d0)call errstop('WRITE_SPIN_DENSITY_MATR&
&IX','Zero weight.')
xyzzyaaad29=1.d0/(xyzzyaaeb1(xyzzyaaac29)*real(npcells,dp))
xyzzyaaae29=xyzzyaaec1(xyzzyaaac29)*real(npcells,dp)
do xyzzyaaab29=1,4
do xyzzyaaaa29=1,expval_ngvec(xyzzyaady1)
xyzzyaaaf29=xyzzyaaae29*xyzzyaaee1(xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29&
&)+xyzzyaaef1(xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29)
xyzzyaaef1(xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29)=xyzzyaaad29*xyzzyaaaf2&
&9
enddo
enddo
endif
enddo
else
xyzzyaaeb1(:)=xyzzyaaec1(:)
endif
write(xyzzyaaaf1,'(a)')'START SPIN-DENSITY MATRIX'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_sdenm
write(xyzzyaaaf1,'(a)')'Use G-vector set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaady1))
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaadz1))
do xyzzyaaac29=1,xyzzyaadz1
if(nele(xyzzyaaac29)/=0)then
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaac29))
write(xyzzyaaaf1,'(a)')'Particle type'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaea1(xyzzyaaac29)))
write(xyzzyaaaf1,'(a)')'Total weight'
write(xyzzyaaaf1,*)xyzzyaaeb1(xyzzyaaac29)
write(xyzzyaaaf1,'(a)')'Complex charge-density coefficients (real part&
&, imaginary part)'
if(spin_density_mat)then
if(xyzzyaaeg1)then
do xyzzyaaab29=1,4
write(xyzzyaaaf1,'(a)')'Matrix Component'
select case(xyzzyaaab29)
case(1)
write(xyzzyaaaf1,'(a)')' 1 1'
case(2)
write(xyzzyaaaf1,'(a)')' 1 2'
case(3)
write(xyzzyaaaf1,'(a)')' 2 1'
case(4)
write(xyzzyaaaf1,'(a)')' 2 2'
end select
do xyzzyaaaa29=1,expval_ngvec(xyzzyaady1)
xyzzyaaaf29=xyzzyaaef1(xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29)
write(xyzzyaaaf1,*)real(xyzzyaaaf29,dp),aimag(xyzzyaaaf29)
enddo
enddo
else
do xyzzyaaab29=1,4
write(xyzzyaaaf1,'(a)')'Matrix Component'
select case(xyzzyaaab29)
case(1)
write(xyzzyaaaf1,'(a)')' 1 1'
case(2)
write(xyzzyaaaf1,'(a)')' 1 2'
case(3)
write(xyzzyaaaf1,'(a)')' 2 1'
case(4)
write(xyzzyaaaf1,'(a)')' 2 2'
end select
do xyzzyaaaa29=1,expval_ngvec(xyzzyaady1)
xyzzyaaad29=real(xyzzyaaef1(xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29),dp)
write(xyzzyaaaf1,*)xyzzyaaad29,' 0.'
enddo
enddo
endif
else
do xyzzyaaab29=1,4
write(xyzzyaaaf1,'(a)')'Matrix Component'
select case(xyzzyaaab29)
case(1)
write(xyzzyaaaf1,'(a)')' 1 1'
case(2)
write(xyzzyaaaf1,'(a)')' 1 2'
case(3)
write(xyzzyaaaf1,'(a)')' 2 1'
case(4)
write(xyzzyaaaf1,'(a)')' 2 2'
end select
do xyzzyaaaa29=1,expval_ngvec(xyzzyaady1)
xyzzyaaaf29=xyzzyaaee1(xyzzyaaaa29,xyzzyaaab29,xyzzyaaac29)
write(xyzzyaaaf1,*)real(xyzzyaaaf29,dp),aimag(xyzzyaaaf29)
enddo
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaac29))
else
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaac29))
write(xyzzyaaaf1,'(1x,a)')'No particles of this type.'
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaac29))
endif
enddo
write(xyzzyaaaf1,'(a)')'END SPIN-DENSITY MATRIX'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqf1
subroutine xyzzyaaqg1
implicit none
character(80) tmpr,tmpr2,tmpr3
integer xyzzyaaaa30,xyzzyaaab30
real(dp) xyzzyaaac30,xyzzyaaad30
complex(dp) xyzzyaaae30
if(pair_corr)then
xyzzyaaeo1=czero
call mpi_reduce(expval_pcf,xyzzyaaeo1,expval_ngvec(xyzzyaaeh1)*xyzzyaa&
&ei1,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_pcf in write_pair_corr.')
endif
if(am_master)then
if(pair_corr)then
do xyzzyaaab30=1,xyzzyaaei1
xyzzyaael1(xyzzyaaab30)=real(xyzzyaaeo1(1,xyzzyaaab30),dp)/real(nele(x&
&yzzyaaab30),dp)+xyzzyaaem1(xyzzyaaab30)
if(xyzzyaael1(xyzzyaaab30)==0.d0)call errstop('WRITE_PAIR_CORR','Zero &
&weight.')
xyzzyaaac30=1.d0/(xyzzyaael1(xyzzyaaab30)*real(npcells,dp))
xyzzyaaad30=xyzzyaaem1(xyzzyaaab30)*real(npcells,dp)
do xyzzyaaaa30=1,expval_ngvec(xyzzyaaeh1)
xyzzyaaae30=xyzzyaaad30*xyzzyaaen1(xyzzyaaaa30,xyzzyaaab30)+xyzzyaaeo1&
&(xyzzyaaaa30,xyzzyaaab30)
xyzzyaaeo1(xyzzyaaaa30,xyzzyaaab30)=xyzzyaaac30*xyzzyaaae30
enddo
enddo
else
xyzzyaael1(:)=xyzzyaaem1(:)
endif
write(xyzzyaaaf1,'(a)')'START RECIPROCAL-SPACE PCF'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_pcf
write(xyzzyaaaf1,'(a)')'Fixed particle type'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaej1))
write(xyzzyaaaf1,'(a)')'Fixed particle position (au)'
tmpr=r2s(pcf_rfix(1),'(f18.12)')
tmpr2=r2s(pcf_rfix(2),'(f18.12)')
tmpr3=r2s(pcf_rfix(3),'(f18.12)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
write(xyzzyaaaf1,'(a)')'Use G-vector set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaeh1))
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaei1))
do xyzzyaaab30=1,xyzzyaaei1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab30))
write(xyzzyaaaf1,'(a)')'Particle type'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaek1(xyzzyaaab30)))
write(xyzzyaaaf1,'(a)')'Total weight'
write(xyzzyaaaf1,*)xyzzyaael1(xyzzyaaab30)
write(xyzzyaaaf1,'(a)')'Complex charge-density coefficients (real part&
&, imaginary part)'
if(pair_corr)then
if(xyzzyaaeq1)then
do xyzzyaaaa30=1,expval_ngvec(xyzzyaaeh1)
xyzzyaaae30=xyzzyaaeo1(xyzzyaaaa30,xyzzyaaab30)
write(xyzzyaaaf1,*)real(xyzzyaaae30,dp),aimag(xyzzyaaae30)
enddo
else
do xyzzyaaaa30=1,expval_ngvec(xyzzyaaeh1)
xyzzyaaac30=real(xyzzyaaeo1(xyzzyaaaa30,xyzzyaaab30),dp)
write(xyzzyaaaf1,*)xyzzyaaac30,' 0.'
enddo
endif
else
do xyzzyaaaa30=1,expval_ngvec(xyzzyaaeh1)
xyzzyaaae30=xyzzyaaen1(xyzzyaaaa30,xyzzyaaab30)
write(xyzzyaaaf1,*)real(xyzzyaaae30,dp),aimag(xyzzyaaae30)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab30))
enddo
write(xyzzyaaaf1,'(a)')'END RECIPROCAL-SPACE PCF'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqg1
subroutine xyzzyaaqh1
implicit none
integer xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31,xyzzyaaad31
real(dp) xyzzyaaae31,xyzzyaaaf31,xyzzyaaag31
character(80) tmpr,tmpr2,tmpr3
if(pair_corr_sph)then
xyzzyaafe1=0.d0
xyzzyaaez1=0.d0
if(xyzzyaaet1==1)then
call mpi_reduce(xyzzyaafc1,xyzzyaafe1,xyzzyaaes1*nspin,mpi_double_prec&
&ision,mpi_sum,0,mpi_comm_world,ierror)
else
call mpi_reduce(xyzzyaafc1,xyzzyaafe1,xyzzyaaes1*nspin*nspin,mpi_doubl&
&e_precision,mpi_sum,0,mpi_comm_world,ierror)
endif
call checkmpi(ierror,'mpi_reduce of expval_pcfs in write_pair_corr_sph&
&.')
call mpi_reduce(xyzzyaaey1,xyzzyaaez1,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_pcfs_weight in write_pair_c&
&orr_sph.')
endif
if(am_master)then
if(pair_corr_sph)then
do xyzzyaaad31=1,xyzzyaaer1
xyzzyaaab31=xyzzyaaev1(1,xyzzyaaad31)
xyzzyaaac31=xyzzyaaev1(2,xyzzyaaad31)
xyzzyaaae31=xyzzyaaez1+xyzzyaafb1(xyzzyaaad31)
if(xyzzyaaae31==0.d0)call errstop('WRITE_PAIR_CORR_SPH','Zero weight f&
&or set '//trim(i2s(xyzzyaaad31))//'.')
do xyzzyaaaa31=1,xyzzyaaes1
xyzzyaaaf31=xyzzyaafb1(xyzzyaaad31)*xyzzyaafd1(xyzzyaaaa31,xyzzyaaab31&
&,xyzzyaaac31)
xyzzyaaag31=xyzzyaafe1(xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31)
xyzzyaafe1(xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31)=(xyzzyaaaf31+xyzzyaaag&
&31)/xyzzyaaae31
enddo
enddo
endif
write(xyzzyaaaf1,'(a)')'START SPHERICAL PCF'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_pcfs
write(xyzzyaaaf1,'(a)')'Number of bins'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaes1))
write(xyzzyaaaf1,'(a)')'Cutoff radius (au)'
tmpr=r2s(xyzzyaaew1,'(f18.12)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
write(xyzzyaaaf1,'(a)')'Accumulation mode (1=fixed particle,2=homogene&
&ous)'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaet1))
if(xyzzyaaet1==1)then
write(xyzzyaaaf1,'(a)')'Fixed particle type'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaeu1))
write(xyzzyaaaf1,'(a)')'Fixed particle position (au)'
tmpr=r2s(pcfs_rfix(1),'(f18.12)')
tmpr2=r2s(pcfs_rfix(2),'(f18.12)')
tmpr3=r2s(pcfs_rfix(3),'(f18.12)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(tmpr3)
endif
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaer1))
do xyzzyaaad31=1,xyzzyaaer1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaad31))
write(xyzzyaaaf1,'(a)')'Types of particles in pair'
xyzzyaaab31=xyzzyaaev1(1,xyzzyaaad31)
xyzzyaaac31=xyzzyaaev1(2,xyzzyaaad31)
write(xyzzyaaaf1,'(2(1x,a))')trim(i2s(xyzzyaaab31)),trim(i2s(xyzzyaaac&
&31))
write(xyzzyaaaf1,'(a)')'Total weight'
if(pair_corr_sph)then
if(isvmc)then
write(xyzzyaaaf1,*)anint(xyzzyaaez1+xyzzyaafb1(xyzzyaaad31))
else
write(xyzzyaaaf1,*)xyzzyaaez1+xyzzyaafb1(xyzzyaaad31)
endif
else
write(xyzzyaaaf1,*)xyzzyaafb1(xyzzyaaad31)
endif
write(xyzzyaaaf1,'(a)')'Bin contents'
if(pair_corr_sph)then
do xyzzyaaaa31=1,xyzzyaaes1
write(xyzzyaaaf1,*)xyzzyaafe1(xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31)
enddo
else
do xyzzyaaaa31=1,xyzzyaaes1
write(xyzzyaaaf1,*)xyzzyaafd1(xyzzyaaaa31,xyzzyaaab31,xyzzyaaac31)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaad31))
enddo
write(xyzzyaaaf1,'(a)')'END SPHERICAL PCF'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqh1
subroutine xyzzyaaqi1
implicit none
integer xyzzyaaaa32,xyzzyaaab32,xyzzyaaac32
real(dp) xyzzyaaad32,xyzzyaaae32,xyzzyaaaf32,xyzzyaaag32,xyzzyaaah32,x&
&yzzyaaai32,xyzzyaaaj32(periodicity,nspin),xyzzyaaak32,xyzzyaaal32,xyz&
&zyaaam32,xyzzyaaan32(periodicity,nspin),xyzzyaaao32(periodicity,perio&
&dicity,nspin)
complex(dp) xyzzyaaap32
logical xyzzyaaaq32(periodicity,periodicity,nspin)
character(80) tmpr
character(256) xyzzyaaah1
xyzzyaaaj32=0.d0
xyzzyaaan32=0.d0
xyzzyaaao32=0.d0
xyzzyaaaq32=.false.
select case(periodicity)
case(1)
xyzzyaaam32=xyzzyaacq1
case(2)
xyzzyaaam32=xyzzyaacp1
case(3)
xyzzyaaam32=xyzzyaacr1
end select
if(loc_tensor)then
xyzzyaafn1(:,:,:)=czero
call mpi_reduce(xyzzyaafl1,xyzzyaafn1,periodicity*periodicity*xyzzyaaf&
&g1,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_lt in write_localization_te&
&nsor.')
call mpi_reduce(xyzzyaafh1,xyzzyaafj1,xyzzyaafg1,mpi_double_precision,&
&mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_lt_weight in write_localiza&
&tion_tensor.')
endif
if(am_master)then
if(loc_tensor)then
do xyzzyaaac32=1,xyzzyaafg1
if(xyzzyaafj1(xyzzyaaac32)==0.d0)call errstop('WRITE_LOC_TENSOR','Zero&
& weight.')
xyzzyaaag32=xyzzyaafi1(xyzzyaaac32)
xyzzyaaak32=1.d0/(xyzzyaafj1(xyzzyaaac32)+xyzzyaaag32)
do xyzzyaaaa32=1,periodicity
do xyzzyaaab32=1,periodicity
xyzzyaaap32=xyzzyaaag32*xyzzyaafm1(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32&
&)+xyzzyaafn1(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32)
xyzzyaafn1(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32)=xyzzyaaak32*xyzzyaaap3&
&2
enddo
enddo
enddo
else
xyzzyaafn1(:,:,:)=xyzzyaafm1(:,:,:)
xyzzyaafj1(:)=xyzzyaafi1(:)
endif
do xyzzyaaac32=1,xyzzyaafg1
do xyzzyaaaa32=1,periodicity
xyzzyaaah32=sqrt(dot_product(bmat(1:periodicity,xyzzyaaaa32),bmat(1:pe&
&riodicity,xyzzyaaaa32)))
xyzzyaaak32=sqrt(dble(xyzzyaafn1(xyzzyaaaa32,xyzzyaaaa32,xyzzyaaac32)*&
&conjg(xyzzyaafn1(xyzzyaaaa32,xyzzyaaaa32,xyzzyaaac32))))
xyzzyaaae32=dble(nele(xyzzyaaac32))*xyzzyaaah32
if(xyzzyaaae32/=0.d0.and.xyzzyaaak32/=0.d0)then
xyzzyaaad32=-xyzzyaaae32*xyzzyaaah32
xyzzyaaao32(xyzzyaaaa32,xyzzyaaaa32,xyzzyaaac32)=(2.d0/xyzzyaaad32)*lo&
&g(xyzzyaaak32)
xyzzyaaaj32(xyzzyaaaa32,xyzzyaaac32)=(1.d0/xyzzyaaae32)*aimag(log(xyzz&
&yaafn1(xyzzyaaaa32,xyzzyaaaa32,xyzzyaaac32)))
xyzzyaaan32(xyzzyaaaa32,xyzzyaaac32)=(dble(nele(xyzzyaaac32))/dble(xyz&
&zyaaam32))*xyzzyaaaj32(xyzzyaaaa32,xyzzyaaac32)
else
xyzzyaaaq32(xyzzyaaaa32,xyzzyaaaa32,xyzzyaaac32)=.true.
endif
do xyzzyaaab32=xyzzyaaaa32+1,periodicity
xyzzyaaai32=sqrt(dot_product(bmat(1:periodicity,xyzzyaaab32),bmat(1:pe&
&riodicity,xyzzyaaab32)))
xyzzyaaad32=-dble(nele(xyzzyaaac32))*xyzzyaaah32*xyzzyaaai32
xyzzyaaal32=sqrt(dble(xyzzyaafn1(xyzzyaaab32,xyzzyaaab32,xyzzyaaac32)*&
&conjg(xyzzyaafn1(xyzzyaaab32,xyzzyaaab32,xyzzyaaac32))))
xyzzyaaaf32=sqrt(dble(xyzzyaafn1(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32)*&
&conjg(xyzzyaafn1(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32))))
xyzzyaaao32(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32)=0.d0
if(xyzzyaaad32/=0.d0.and.xyzzyaaak32/=0.d0.and.xyzzyaaal32/=0.d0.and.x&
&yzzyaaaf32/=0.d0)then
xyzzyaaao32(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32)=(log(xyzzyaaak32)+log&
&(xyzzyaaal32)-log(xyzzyaaaf32))/xyzzyaaad32
xyzzyaaao32(xyzzyaaaa32,xyzzyaaab32,xyzzyaaac32)=xyzzyaaao32(xyzzyaaab&
&32,xyzzyaaaa32,xyzzyaaac32)
else
xyzzyaaaq32(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32)=.true.
xyzzyaaaq32(xyzzyaaaa32,xyzzyaaab32,xyzzyaaac32)=.true.
endif
enddo
enddo
enddo
write(xyzzyaaaf1,'(a)')'START LOCALIZATION TENSOR'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_lt
do xyzzyaaac32=1,xyzzyaafg1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaac32))
write(xyzzyaaaf1,'(a)')'Total weight'
write(xyzzyaaaf1,*)xyzzyaafj1(xyzzyaaac32)
write(xyzzyaaaf1,'(a)')'Expectation value  <exp(ik.r)>c'
if(xyzzyaafp1)then
do xyzzyaaaa32=1,periodicity
do xyzzyaaab32=1,periodicity
write(xyzzyaaaf1,*)real(xyzzyaafn1(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32&
&),dp),aimag(xyzzyaafn1(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32))
enddo
enddo
else
do xyzzyaaaa32=1,periodicity
do xyzzyaaab32=1,periodicity
write(xyzzyaaaf1,*)real(xyzzyaafn1(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32&
&),dp),0.d0
enddo
enddo
endif
write(xyzzyaaaf1,'(a)')'LOCALIZATION TENSOR ELEMENTS <r^2>c'
do xyzzyaaaa32=1,periodicity
xyzzyaaah1=''
do xyzzyaaab32=1,periodicity
tmpr=''
if(xyzzyaaaq32(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32))then
tmpr='Infinity'
else
tmpr=r2s2(xyzzyaaao32(xyzzyaaab32,xyzzyaaaa32,xyzzyaaac32),'(es)')
endif
xyzzyaaah1=trim(xyzzyaaah1)//' '//trim(tmpr)
enddo
write(xyzzyaaaf1,'(a)')trim(xyzzyaaah1)
enddo
write(xyzzyaaaf1,'(a)')'DIPOLE MOMENT <r>c'
xyzzyaaah1=''
do xyzzyaaaa32=1,periodicity
tmpr=''
if(xyzzyaaaq32(xyzzyaaaa32,xyzzyaaaa32,xyzzyaaac32))then
tmpr='Infinity'
else
tmpr=r2s2(xyzzyaaaj32(xyzzyaaaa32,xyzzyaaac32),'(es)')
endif
xyzzyaaah1=trim(xyzzyaaah1)//' '//trim(tmpr)
enddo
write(xyzzyaaaf1,'(a)')trim(xyzzyaaah1)
write(xyzzyaaaf1,'(a)')'POLARIZATION'
xyzzyaaah1=''
do xyzzyaaaa32=1,periodicity
tmpr=''
if(xyzzyaaaq32(xyzzyaaaa32,xyzzyaaaa32,xyzzyaaac32))then
tmpr='Infinity'
else
tmpr=r2s2(xyzzyaaan32(xyzzyaaaa32,xyzzyaaac32),'(es)')
endif
xyzzyaaah1=trim(xyzzyaaah1)//' '//trim(tmpr)
enddo
write(xyzzyaaaf1,'(a)')trim(xyzzyaaah1)
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaac32))
enddo
write(xyzzyaaaf1,'(a)')'END LOCALIZATION TENSOR'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqi1
subroutine xyzzyaaqj1
implicit none
integer xyzzyaaaa33,xyzzyaaab33,xyzzyaaac33
real(dp) xyzzyaaad33,xyzzyaaae33,xyzzyaaaf33
complex(dp) xyzzyaaag33
if(structure_factor)then
xyzzyaagb1(:,:)=0.d0
call mpi_reduce(xyzzyaaga1,xyzzyaagb1,expval_ngvec(xyzzyaafr1)*xyzzyaa&
&fq1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_sf in write_structure_facto&
&r.')
call mpi_reduce(xyzzyaafv1,xyzzyaafw1,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_sf_weight in write_structur&
&e_factor.')
if(.not.homogeneous_density)then
xyzzyaagi1(:,:)=czero
call mpi_reduce(xyzzyaagg1,xyzzyaagi1,expval_ngvec(xyzzyaafr1)*xyzzyaa&
&fs1,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_sfsden in write_structure_f&
&actor')
endif
endif
if(am_master)then
if(structure_factor)then
do xyzzyaaac33=1,xyzzyaafq1
xyzzyaaad33=(xyzzyaafw1+xyzzyaafy1(xyzzyaaac33))*(real(npcells,dp)**2)
if(xyzzyaaad33==0.d0)call errstop('WRITE_STRUCTURE_FACTOR','Zero weigh&
&t for structure factor.')
do xyzzyaaab33=1,expval_ngvec(xyzzyaafr1)
xyzzyaaae33=real(npcells**2,dp)*xyzzyaafy1(xyzzyaaac33)*xyzzyaafz1(xyz&
&zyaaab33,xyzzyaaac33)+xyzzyaagb1(xyzzyaaab33,xyzzyaaac33)
xyzzyaagb1(xyzzyaaab33,xyzzyaaac33)=xyzzyaaae33/xyzzyaaad33
enddo
enddo
do xyzzyaaac33=1,xyzzyaafs1
if(nele(xyzzyaaac33)/=0)then
xyzzyaage1(xyzzyaaac33)=real(xyzzyaagi1(1,xyzzyaaac33),dp)/real(nele(x&
&yzzyaaac33),dp)+xyzzyaagd1(xyzzyaaac33)
if(xyzzyaage1(xyzzyaaac33)==0.d0)call errstop('WRITE_STRUCTURE_FACTOR'&
&,'Zero weight for spin-density part of structure factor.')
xyzzyaaad33=1.d0/(xyzzyaage1(xyzzyaaac33)*real(npcells,dp))
xyzzyaaaf33=real(npcells,dp)*xyzzyaagd1(xyzzyaaac33)
do xyzzyaaaa33=1,expval_ngvec(xyzzyaafr1)
xyzzyaaag33=xyzzyaaaf33*xyzzyaagh1(xyzzyaaaa33,xyzzyaaac33)+xyzzyaagi1&
&(xyzzyaaaa33,xyzzyaaac33)
xyzzyaagi1(xyzzyaaaa33,xyzzyaaac33)=xyzzyaaag33*xyzzyaaad33
enddo
endif
enddo
endif
write(xyzzyaaaf1,'(a)')'START STRUCTURE FACTOR'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_sf
write(xyzzyaaaf1,'(a)')'Use G-vector set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaafr1))
write(xyzzyaaaf1,'(a)')'Number of sets for rho_a(G)*rho_b(-G)'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaafq1))
do xyzzyaaac33=1,xyzzyaafq1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaac33))
write(xyzzyaaaf1,'(a)')'Types of particles in pair'
write(xyzzyaaaf1,'(2(1x,a))')trim(i2s(xyzzyaaft1(1,xyzzyaaac33))),trim&
&(i2s(xyzzyaaft1(2,xyzzyaaac33)))
if(nele(xyzzyaaft1(1,xyzzyaaac33))==0.or.nele(xyzzyaaft1(2,xyzzyaaac33&
&))==0)then
write(xyzzyaaaf1,'(1x,a)')'No such pairs.'
else
write(xyzzyaaaf1,'(a)')'Total weight'
if(structure_factor)then
if(isvmc)then
write(xyzzyaaaf1,*)anint(xyzzyaafw1+xyzzyaafy1(xyzzyaaac33))
else
write(xyzzyaaaf1,*)xyzzyaafw1+xyzzyaafy1(xyzzyaaac33)
endif
else
write(xyzzyaaaf1,*)xyzzyaafy1(xyzzyaaac33)
endif
write(xyzzyaaaf1,'(a)')'rho_a(G)*rho_b(-G)'
if(structure_factor)then
do xyzzyaaab33=1,expval_ngvec(xyzzyaafr1)
write(xyzzyaaaf1,*)xyzzyaagb1(xyzzyaaab33,xyzzyaaac33)
enddo
else
do xyzzyaaab33=1,expval_ngvec(xyzzyaafr1)
write(xyzzyaaaf1,*)xyzzyaafz1(xyzzyaaab33,xyzzyaaac33)
enddo
endif
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaac33))
enddo
write(xyzzyaaaf1,'(a)')'Number of sets for spin density part'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaafs1))
do xyzzyaaac33=1,xyzzyaafs1
if(nele(xyzzyaaac33)/=0)then
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaac33))
write(xyzzyaaaf1,'(a)')'Particle type'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaafu1(xyzzyaaac33)))
write(xyzzyaaaf1,'(a)')'Total weight'
if(structure_factor)then
write(xyzzyaaaf1,*)xyzzyaage1(xyzzyaaac33)
else
write(xyzzyaaaf1,*)xyzzyaagd1(xyzzyaaac33)
endif
write(xyzzyaaaf1,'(a)')'Complex charge-density coefficients (real part&
&, imaginary part)'
if(structure_factor)then
if(xyzzyaagk1)then
do xyzzyaaaa33=1,expval_ngvec(xyzzyaafr1)
xyzzyaaag33=xyzzyaagi1(xyzzyaaaa33,xyzzyaaac33)
write(xyzzyaaaf1,*)real(xyzzyaaag33,dp),aimag(xyzzyaaag33)
enddo
else
do xyzzyaaaa33=1,expval_ngvec(xyzzyaafr1)
xyzzyaaad33=real(xyzzyaagi1(xyzzyaaaa33,xyzzyaaac33),dp)
write(xyzzyaaaf1,*)xyzzyaaad33,' 0.'
enddo
endif
else
do xyzzyaaaa33=1,expval_ngvec(xyzzyaafr1)
xyzzyaaag33=xyzzyaagh1(xyzzyaaaa33,xyzzyaaac33)
write(xyzzyaaaf1,*)real(xyzzyaaag33,dp),aimag(xyzzyaaag33)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaac33))
else
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaac33))
write(xyzzyaaaf1,'(a)')'Particle type'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaafu1(xyzzyaaac33)))
write(xyzzyaaaf1,'(1x,a)')'No such particles.'
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaac33))
endif
enddo
write(xyzzyaaaf1,'(a)')'END STRUCTURE FACTOR'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqj1
subroutine xyzzyaaqk1
implicit none
integer xyzzyaaaa34,xyzzyaaab34,xyzzyaaac34,xyzzyaaad34
real(dp) xyzzyaaae34,xyzzyaaaf34,xyzzyaaag34
character(80) tmpr,tmpr2
if(structure_factor_sph)then
xyzzyaagv1(:,:)=0.d0
xyzzyaags1(:)=0.d0
call mpi_reduce(xyzzyaagu1,xyzzyaagv1,xyzzyaagl1*xyzzyaagm1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_sf_sph in write_structure_f&
&actor_sph.')
call mpi_reduce(xyzzyaagq1,xyzzyaags1,xyzzyaagm1,mpi_double_precision,&
&mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_sf_sph_weight in write_stru&
&cture_factor_sph.')
endif
if(am_master)then
if(structure_factor_sph)then
do xyzzyaaaa34=1,nspin
do xyzzyaaab34=xyzzyaaaa34,nspin
xyzzyaaad34=which_spair(xyzzyaaab34,xyzzyaaaa34,levels_spairs)
if(xyzzyaags1(xyzzyaaad34)==0.d0)call errstop('WRITE_STRUCTURE_FACTOR_&
&SPH','Zero weight for set ' //trim(i2s(xyzzyaaad34))//'.')
call dscal(xyzzyaagl1,2.d0/xyzzyaags1(xyzzyaaad34),xyzzyaagv1(1,xyzzya&
&aad34),1)
if(xyzzyaaaa34==xyzzyaaab34)xyzzyaagv1(1:xyzzyaagl1,xyzzyaaad34)=xyzzy&
&aagv1(1:xyzzyaagl1,xyzzyaaad34)+dble(nele(xyzzyaaaa34))/dble(netot)
enddo
enddo
if(dimensionality==3)then
xyzzyaaaf34=-fourpi*wigner_seitz_radius**3/(dble(netot)*volume)
elseif(dimensionality==2)then
xyzzyaaaf34=-twopi*wigner_seitz_radius**2/(dble(netot)*area)
else
xyzzyaaaf34=-2.d0*wigner_seitz_radius/(dble(netot)*abs(a1(1)))
endif
do xyzzyaaac34=1,xyzzyaagl1
if(xyzzyaagt1(xyzzyaaac34)/=0.d0)then
xyzzyaaae34=xyzzyaagt1(xyzzyaaac34)*wigner_seitz_radius
if(dimensionality==3)then
xyzzyaaag34=xyzzyaaaf34*(sin(xyzzyaaae34)-xyzzyaaae34*cos(xyzzyaaae34)&
&)/xyzzyaaae34**3
elseif(dimensionality==2)then
xyzzyaaag34=xyzzyaaaf34*bessel_j1(xyzzyaaae34)/xyzzyaaae34
else
xyzzyaaag34=xyzzyaaaf34*sin(xyzzyaaae34)/xyzzyaaae34
endif
else
if(dimensionality==3)then
xyzzyaaag34=xyzzyaaaf34*third
elseif(dimensionality==2)then
xyzzyaaag34=xyzzyaaaf34*0.5d0
else
xyzzyaaag34=xyzzyaaaf34
endif
endif
do xyzzyaaaa34=1,nspin
xyzzyaaad34=which_spair(xyzzyaaaa34,xyzzyaaaa34,levels_spairs)
xyzzyaagv1(xyzzyaaac34,xyzzyaaad34)=xyzzyaagv1(xyzzyaaac34,xyzzyaaad34&
&)+dble(nele(xyzzyaaaa34)**2)*xyzzyaaag34
do xyzzyaaab34=xyzzyaaaa34+1,nspin
xyzzyaaad34=which_spair(xyzzyaaab34,xyzzyaaaa34,levels_spairs)
xyzzyaagv1(xyzzyaaac34,xyzzyaaad34)=xyzzyaagv1(xyzzyaaac34,xyzzyaaad34&
&)+dble(2*nele(xyzzyaaaa34)*nele(xyzzyaaab34))*xyzzyaaag34
enddo
enddo
enddo
endif
do xyzzyaaad34=1,xyzzyaagm1
if(xyzzyaagy1(xyzzyaaad34)>0.d0)then
call dscal(xyzzyaagl1,xyzzyaags1(xyzzyaaad34)/(xyzzyaags1(xyzzyaaad34)&
&+xyzzyaagy1(xyzzyaaad34)),xyzzyaagv1(1,xyzzyaaad34),1)
xyzzyaags1(xyzzyaaad34)=xyzzyaags1(xyzzyaaad34)+xyzzyaagy1(xyzzyaaad34&
&)
call daxpy(xyzzyaagl1,xyzzyaagy1(xyzzyaaad34)/xyzzyaags1(xyzzyaaad34),&
&xyzzyaagx1(1,xyzzyaaad34),1,xyzzyaagv1(1,xyzzyaaad34),1)
endif
enddo
write(xyzzyaaaf1,'(a)')'START SPHERICAL STRUCTURE FACTOR'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_sf_sph
write(xyzzyaaaf1,'(a)')'Radial k point grid'
tmpr=r2s(xyzzyaago1,'(f18.12)')
tmpr2=r2s(xyzzyaagp1,'(f18.12)')
write(xyzzyaaaf1,'(3(1x,a))')trim(tmpr),trim(tmpr2),trim(i2s(xyzzyaagl&
&1))
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaagm1))
do xyzzyaaad34=1,xyzzyaagm1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaad34))
write(xyzzyaaaf1,'(a)')'Particle types'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaagn1(1,xyzzyaaad34)))//' ' //&
&trim(i2s(xyzzyaagn1(2,xyzzyaaad34)))
write(xyzzyaaaf1,'(a)')'Total weight'
write(xyzzyaaaf1,*)xyzzyaags1(xyzzyaaad34)
write(xyzzyaaaf1,'(a)')'k,SF(k)'
do xyzzyaaac34=1,xyzzyaagl1
write(xyzzyaaaf1,*)xyzzyaagt1(xyzzyaaac34),xyzzyaagv1(xyzzyaaac34,xyzz&
&yaaad34)
enddo
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaad34))
enddo
write(xyzzyaaaf1,'(a)')'END SPHERICAL STRUCTURE FACTOR'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqk1
subroutine xyzzyaaql1
implicit none
integer xyzzyaaaa35,xyzzyaaab35
real(dp) xyzzyaaac35,xyzzyaaad35
character(80) char_80
if(onep_density_mat)then
xyzzyaahm1=0.d0
xyzzyaahn1=0.d0
call mpi_reduce(xyzzyaahi1,xyzzyaahm1,xyzzyaaha1*xyzzyaagz1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_onep_dm in write_onep_densi&
&ty_mat.')
call mpi_reduce(xyzzyaahj1,xyzzyaahn1,xyzzyaaha1*xyzzyaagz1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_onep_dm2 in write_onep_dens&
&ity_mat.')
call mpi_reduce(xyzzyaahf1,xyzzyaahh1,xyzzyaaha1*xyzzyaagz1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_onep_dm_weight in write_one&
&p_density_mat.')
endif
if(am_master)then
if(onep_density_mat)then
do xyzzyaaab35=1,xyzzyaagz1
do xyzzyaaaa35=1,xyzzyaaha1
xyzzyaahh1(xyzzyaaaa35,xyzzyaaab35)=xyzzyaahg1(xyzzyaaaa35,xyzzyaaab35&
&)+xyzzyaahh1(xyzzyaaaa35,xyzzyaaab35)
xyzzyaaac35=xyzzyaahm1(xyzzyaaaa35,xyzzyaaab35)+xyzzyaahk1(xyzzyaaaa35&
&,xyzzyaaab35)*xyzzyaahg1(xyzzyaaaa35,xyzzyaaab35)
xyzzyaaad35=xyzzyaahn1(xyzzyaaaa35,xyzzyaaab35)+xyzzyaahl1(xyzzyaaaa35&
&,xyzzyaaab35)*xyzzyaahg1(xyzzyaaaa35,xyzzyaaab35)
if(xyzzyaahh1(xyzzyaaaa35,xyzzyaaab35)==0.d0)then
xyzzyaahm1(xyzzyaaaa35,xyzzyaaab35)=0.d0
xyzzyaahl1(xyzzyaaaa35,xyzzyaaab35)=0.d0
else
xyzzyaahm1(xyzzyaaaa35,xyzzyaaab35)=xyzzyaaac35/xyzzyaahh1(xyzzyaaaa35&
&,xyzzyaaab35)
xyzzyaahn1(xyzzyaaaa35,xyzzyaaab35)=xyzzyaaad35/xyzzyaahh1(xyzzyaaaa35&
&,xyzzyaaab35)
endif
enddo
enddo
else
xyzzyaahh1(:,:)=xyzzyaahg1(:,:)
endif
write(xyzzyaaaf1,'(a)')'START ONE-PARTICLE DENSITY MATRIX'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_onep_dm
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaagz1))
write(xyzzyaaaf1,'(a)')'Number of bins'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaha1))
write(xyzzyaaaf1,'(a)')'Number of random points to sample'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaahb1))
do xyzzyaaab35=1,xyzzyaagz1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab35))
write(xyzzyaaaf1,'(a)')'Number of particle types in set'
if(on_top_ii==0)then
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaahc1(xyzzyaaab35)))
char_80=''
do xyzzyaaaa35=1,xyzzyaahc1(xyzzyaaab35)
char_80=trim(adjustl(char_80))//' '//trim(i2s(xyzzyaahd1(xyzzyaaaa35,x&
&yzzyaaab35)))
enddo
else
write(xyzzyaaaf1,'(1x,a)')'on-top pair'
char_80=' '//trim(i2s(which_spin(on_top_ii)))//' '//trim(i2s(which_spi&
&n(on_top_jj)))
endif
write(xyzzyaaaf1,'(a)')'Particle types'
write(xyzzyaaaf1,'(1x,a)')trim(char_80)
write(xyzzyaaaf1,'(a)')'Weight(r),OBDM(r)**2,OBDM(r)'
if(onep_density_mat)then
do xyzzyaaaa35=1,xyzzyaaha1
write(xyzzyaaaf1,*)xyzzyaahh1(xyzzyaaaa35,xyzzyaaab35),xyzzyaahn1(xyzz&
&yaaaa35,xyzzyaaab35),xyzzyaahm1(xyzzyaaaa35,xyzzyaaab35)
enddo
else
do xyzzyaaaa35=1,xyzzyaaha1
write(xyzzyaaaf1,*)xyzzyaahg1(xyzzyaaaa35,xyzzyaaab35),xyzzyaahl1(xyzz&
&yaaaa35,xyzzyaaab35),xyzzyaahk1(xyzzyaaaa35,xyzzyaaab35)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab35))
enddo
write(xyzzyaaaf1,'(a)')'END ONE-PARTICLE DENSITY MATRIX'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaql1
subroutine xyzzyaaqm1
implicit none
integer xyzzyaaaa36,xyzzyaaab36
real(dp) xyzzyaaac36,xyzzyaaad36
character(80) char_80,tmpr
if(twop_density_mat)then
xyzzyaajd1=0.d0
xyzzyaaje1=0.d0
call mpi_reduce(xyzzyaaiz1,xyzzyaajd1,xyzzyaail1*xyzzyaaik1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_twop_dm in write_twop_densi&
&ty_mat.')
call mpi_reduce(xyzzyaaja1,xyzzyaaje1,xyzzyaail1*xyzzyaaik1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_twop_dm2 in write_twop_dens&
&ity_mat.')
call mpi_reduce(xyzzyaaiw1,xyzzyaaiy1,xyzzyaail1*xyzzyaaik1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_twop_dm_weight in write_two&
&p_density_mat.')
endif
if(am_master)then
if(twop_density_mat)then
do xyzzyaaab36=1,xyzzyaaik1
do xyzzyaaaa36=1,xyzzyaail1
xyzzyaaiy1(xyzzyaaaa36,xyzzyaaab36)=xyzzyaaix1(xyzzyaaaa36,xyzzyaaab36&
&)+xyzzyaaiy1(xyzzyaaaa36,xyzzyaaab36)
xyzzyaaac36=xyzzyaajd1(xyzzyaaaa36,xyzzyaaab36)+xyzzyaajb1(xyzzyaaaa36&
&,xyzzyaaab36)*xyzzyaaix1(xyzzyaaaa36,xyzzyaaab36)
xyzzyaaad36=xyzzyaaje1(xyzzyaaaa36,xyzzyaaab36)+xyzzyaajc1(xyzzyaaaa36&
&,xyzzyaaab36)*xyzzyaaix1(xyzzyaaaa36,xyzzyaaab36)
if(xyzzyaaiy1(xyzzyaaaa36,xyzzyaaab36)==0.d0)then
xyzzyaajd1(xyzzyaaaa36,xyzzyaaab36)=0.d0
xyzzyaaje1(xyzzyaaaa36,xyzzyaaab36)=0.d0
else
xyzzyaajd1(xyzzyaaaa36,xyzzyaaab36)=xyzzyaaac36/xyzzyaaiy1(xyzzyaaaa36&
&,xyzzyaaab36)
xyzzyaaje1(xyzzyaaaa36,xyzzyaaab36)=xyzzyaaad36/xyzzyaaiy1(xyzzyaaaa36&
&,xyzzyaaab36)
endif
enddo
enddo
else
xyzzyaaiy1(:,:)=xyzzyaaix1(:,:)
endif
write(xyzzyaaaf1,'(a)')'START TWO-PARTICLE DENSITY MATRIX'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_twop_dm
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaik1))
write(xyzzyaaaf1,'(a)')'Number of bins'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaail1))
write(xyzzyaaaf1,'(a)')'Number of random points to sample'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaim1))
write(xyzzyaaaf1,'(a)')'Fraction of particle pairs to sample at random&
&'
tmpr=r2s(xyzzyaaiv1,'(f16.4)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
do xyzzyaaab36=1,xyzzyaaik1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab36))
write(xyzzyaaaf1,'(a)')'Number of particle-pair types in set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaain1(xyzzyaaab36)))
char_80=''
do xyzzyaaaa36=1,xyzzyaain1(xyzzyaaab36)
char_80=trim(adjustl(char_80))//'  '//trim(i2s(xyzzyaaio1(1,xyzzyaaaa3&
&6,xyzzyaaab36)))//' '//trim(i2s(xyzzyaaio1(2,xyzzyaaaa36,xyzzyaaab36)&
&))
enddo
write(xyzzyaaaf1,'(a)')'Particle-pair types'
write(xyzzyaaaf1,'(1x,a)')trim(adjustl(char_80))
write(xyzzyaaaf1,'(a)')'Weight(r),TBDM(r)**2,TBDM(r)'
if(twop_density_mat)then
do xyzzyaaaa36=1,xyzzyaail1
write(xyzzyaaaf1,*)xyzzyaaiy1(xyzzyaaaa36,xyzzyaaab36),xyzzyaaje1(xyzz&
&yaaaa36,xyzzyaaab36),xyzzyaajd1(xyzzyaaaa36,xyzzyaaab36)
enddo
else
do xyzzyaaaa36=1,xyzzyaail1
write(xyzzyaaaf1,*)xyzzyaaix1(xyzzyaaaa36,xyzzyaaab36),xyzzyaajc1(xyzz&
&yaaaa36,xyzzyaaab36),xyzzyaajb1(xyzzyaaaa36,xyzzyaaab36)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab36))
enddo
write(xyzzyaaaf1,'(a)')'END TWO-PARTICLE DENSITY MATRIX'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqm1
subroutine xyzzyaaqn1
implicit none
integer xyzzyaaaa37,xyzzyaaab37
real(dp) xyzzyaaac37,xyzzyaaad37
character(80) char_80,tmpr
if(cond_fraction)then
xyzzyaakc1=0.d0
xyzzyaakd1=0.d0
call mpi_reduce(xyzzyaajy1,xyzzyaakc1,xyzzyaajk1*xyzzyaajj1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_cond_frac in write_cond_fra&
&ction.')
call mpi_reduce(xyzzyaajz1,xyzzyaakd1,xyzzyaajk1*xyzzyaajj1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_cond_frac2 in write_cond_fr&
&action.')
call mpi_reduce(xyzzyaajv1,xyzzyaajx1,xyzzyaajk1*xyzzyaajj1,mpi_double&
&_precision,mpi_sum,0,     mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_cond_frac_weight in write_c&
&ond_fraction.')
endif
if(am_master)then
if(cond_fraction)then
do xyzzyaaab37=1,xyzzyaajj1
do xyzzyaaaa37=1,xyzzyaajk1
xyzzyaajx1(xyzzyaaaa37,xyzzyaaab37)=xyzzyaajw1(xyzzyaaaa37,xyzzyaaab37&
&)+xyzzyaajx1(xyzzyaaaa37,xyzzyaaab37)
xyzzyaaac37=xyzzyaakc1(xyzzyaaaa37,xyzzyaaab37)+xyzzyaaka1(xyzzyaaaa37&
&,xyzzyaaab37)*xyzzyaajw1(xyzzyaaaa37,xyzzyaaab37)
xyzzyaaad37=xyzzyaakd1(xyzzyaaaa37,xyzzyaaab37)+xyzzyaakb1(xyzzyaaaa37&
&,xyzzyaaab37)*xyzzyaajw1(xyzzyaaaa37,xyzzyaaab37)
if(xyzzyaajx1(xyzzyaaaa37,xyzzyaaab37)==0.d0)then
xyzzyaakc1(xyzzyaaaa37,xyzzyaaab37)=0.d0
xyzzyaakd1(xyzzyaaaa37,xyzzyaaab37)=0.d0
else
xyzzyaakc1(xyzzyaaaa37,xyzzyaaab37)=xyzzyaaac37/xyzzyaajx1(xyzzyaaaa37&
&,xyzzyaaab37)
xyzzyaakd1(xyzzyaaaa37,xyzzyaaab37)=xyzzyaaad37/xyzzyaajx1(xyzzyaaaa37&
&,xyzzyaaab37)
endif
enddo
enddo
else
xyzzyaajx1(:,:)=xyzzyaajw1(:,:)
endif
write(xyzzyaaaf1,'(a)')'START CONDENSATE FRACTION'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accum_method_cond_frac
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaajj1))
write(xyzzyaaaf1,'(a)')'Number of bins'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaajk1))
write(xyzzyaaaf1,'(a)')'Number of random points to sample'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaajl1))
write(xyzzyaaaf1,'(a)')'Fraction of particle pairs to sample at random&
&'
tmpr=r2s(xyzzyaaju1,'(f16.4)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
do xyzzyaaab37=1,xyzzyaajj1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab37))
write(xyzzyaaaf1,'(a)')'Number of particle-pair types in set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaajm1(xyzzyaaab37)))
char_80=''
do xyzzyaaaa37=1,xyzzyaajm1(xyzzyaaab37)
char_80=trim(adjustl(char_80))//'  '//trim(i2s(xyzzyaajn1(1,xyzzyaaaa3&
&7,xyzzyaaab37)))//' '//trim(i2s(xyzzyaajn1(2,xyzzyaaaa37,xyzzyaaab37)&
&))
enddo
write(xyzzyaaaf1,'(a)')'Particle-pair types'
write(xyzzyaaaf1,'(1x,a)')trim(adjustl(char_80))
write(xyzzyaaaf1,'(a)')'Weight(r),CF(r)**2,CF(r)'
if(cond_fraction)then
do xyzzyaaaa37=1,xyzzyaajk1
write(xyzzyaaaf1,*)xyzzyaajx1(xyzzyaaaa37,xyzzyaaab37),xyzzyaakd1(xyzz&
&yaaaa37,xyzzyaaab37),xyzzyaakc1(xyzzyaaaa37,xyzzyaaab37)
enddo
else
do xyzzyaaaa37=1,xyzzyaajk1
write(xyzzyaaaf1,*)xyzzyaajw1(xyzzyaaaa37,xyzzyaaab37),xyzzyaakb1(xyzz&
&yaaaa37,xyzzyaaab37),xyzzyaaka1(xyzzyaaaa37,xyzzyaaab37)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab37))
enddo
write(xyzzyaaaf1,'(a)')'END CONDENSATE FRACTION'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqn1
subroutine xyzzyaaqo1
implicit none
integer xyzzyaaaa38,xyzzyaaab38
real(dp) xyzzyaaac38,xyzzyaaad38,xyzzyaaae38,xyzzyaaaf38
character(80) char_80,tmpr
xyzzyaaaf38=1.d0
if(cond_frac_mom)then
xyzzyaaoj1=0.d0
xyzzyaaok1=0.d0
call mpi_reduce(xyzzyaaof1,xyzzyaaoj1,expval_ngvec(1)*xyzzyaanw1,mpi_d&
&ouble_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_cond_frac_mom in write_cond&
&_frac_mom.')
call mpi_reduce(xyzzyaaog1,xyzzyaaok1,expval_ngvec(1)*xyzzyaanw1,mpi_d&
&ouble_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_cond_frac_mom2 in write_con&
&d_frac_mom.')
call mpi_reduce(xyzzyaaoc1,xyzzyaaoe1,xyzzyaanw1,mpi_double_precision,&
&mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_cond_frac_mom_wt in write_c&
&ond_frac_mom.')
endif
if(am_master)then
if(cond_frac_mom)then
do xyzzyaaab38=1,xyzzyaanw1
xyzzyaaoe1(xyzzyaaab38)=xyzzyaaod1(xyzzyaaab38)+xyzzyaaoe1(xyzzyaaab38&
&)
do xyzzyaaaa38=1,expval_ngvec(1)
xyzzyaaac38=xyzzyaaoj1(xyzzyaaaa38,xyzzyaaab38)+xyzzyaaoh1(xyzzyaaaa38&
&,xyzzyaaab38)*xyzzyaaod1(xyzzyaaab38)
xyzzyaaad38=xyzzyaaok1(xyzzyaaaa38,xyzzyaaab38)+xyzzyaaoi1(xyzzyaaaa38&
&,xyzzyaaab38)*xyzzyaaod1(xyzzyaaab38)
if(xyzzyaaoe1(xyzzyaaab38)==0.d0)then
xyzzyaaoj1(xyzzyaaaa38,xyzzyaaab38)=0.d0
xyzzyaaoi1(xyzzyaaaa38,xyzzyaaab38)=0.d0
else
xyzzyaaoj1(xyzzyaaaa38,xyzzyaaab38)=xyzzyaaac38/xyzzyaaoe1(xyzzyaaab38&
&)
xyzzyaaok1(xyzzyaaaa38,xyzzyaaab38)=xyzzyaaad38/xyzzyaaoe1(xyzzyaaab38&
&)
endif
enddo
enddo
else
xyzzyaaoe1(:)=xyzzyaaod1(:)
endif
write(xyzzyaaaf1,'(a)')'START CONDENSATE FRACTION (MOMENTUM SPACE)'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accum_method_cond_frac_mom
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaanw1))
write(xyzzyaaaf1,'(a)')'Number of random points to sample'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaanx1))
write(xyzzyaaaf1,'(a)')'Fraction of particle pairs to sample at random&
&'
tmpr=r2s(xyzzyaaob1,'(f16.4)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
do xyzzyaaab38=1,xyzzyaanw1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab38))
write(xyzzyaaaf1,'(a)')'Number of particle-pair types in set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaany1(xyzzyaaab38)))
char_80=''
do xyzzyaaaa38=1,xyzzyaany1(xyzzyaaab38)
char_80=trim(adjustl(char_80))//'  '//trim(i2s(xyzzyaanz1(1,xyzzyaaaa3&
&8,xyzzyaaab38)))//' '//trim(i2s(xyzzyaanz1(2,xyzzyaaaa38,xyzzyaaab38)&
&))
enddo
write(xyzzyaaaf1,'(a)')'Particle-pair types'
write(xyzzyaaaf1,'(1x,a)')trim(adjustl(char_80))
write(xyzzyaaaf1,'(a)')'Weight for this set'
write(xyzzyaaaf1,*)xyzzyaaoe1(xyzzyaaab38)
write(xyzzyaaaf1,'(a)')'|k|/k_F,cond_frac_mom(k)**2,cond_frac_mom(k)'
if(cond_frac_mom)then
do xyzzyaaaa38=1,expval_ngvec(1)
xyzzyaaae38=sqrt(ddot(3,xyzzyaaoo1(1,xyzzyaaaa38),1,xyzzyaaoo1(1,xyzzy&
&aaaa38),1))*xyzzyaaaf38
write(xyzzyaaaf1,*)xyzzyaaae38,xyzzyaaok1(xyzzyaaaa38,xyzzyaaab38),xyz&
&zyaaoj1(xyzzyaaaa38,xyzzyaaab38)
enddo
else
do xyzzyaaaa38=1,expval_ngvec(1)
xyzzyaaae38=sqrt(ddot(3,xyzzyaaoo1(1,xyzzyaaaa38),1,xyzzyaaoo1(1,xyzzy&
&aaaa38),1))*xyzzyaaaf38
write(xyzzyaaaf1,*)xyzzyaaae38,xyzzyaaoi1(xyzzyaaaa38,xyzzyaaab38),xyz&
&zyaaoh1(xyzzyaaaa38,xyzzyaaab38)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab38))
enddo
write(xyzzyaaaf1,'(a)')'END CONDENSATE FRACTION (MOMENTUM SPACE)'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqo1
subroutine xyzzyaaqp1
implicit none
integer xyzzyaaaa39,xyzzyaaab39
real(dp) xyzzyaaac39,xyzzyaaad39,xyzzyaaae39,xyzzyaaaf39
character(80) char_80
xyzzyaaaf39=1.d0
if(mom_den)then
xyzzyaaid1=0.d0
xyzzyaaie1=0.d0
call mpi_reduce(xyzzyaahz1,xyzzyaaid1,expval_ngvec(1)*xyzzyaahs1,mpi_d&
&ouble_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_mom_den in write_mom_den.')
call mpi_reduce(xyzzyaaia1,xyzzyaaie1,expval_ngvec(1)*xyzzyaahs1,mpi_d&
&ouble_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_mom_den2 in write_mom_den.'&
&)
call mpi_reduce(xyzzyaahw1,xyzzyaahy1,xyzzyaahs1,mpi_double_precision,&
&mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_mom_den_weight in write_mom&
&_den.')
endif
if(am_master)then
if(mom_den)then
do xyzzyaaab39=1,xyzzyaahs1
xyzzyaahy1(xyzzyaaab39)=xyzzyaahx1(xyzzyaaab39)+xyzzyaahy1(xyzzyaaab39&
&)
do xyzzyaaaa39=1,expval_ngvec(1)
xyzzyaaac39=xyzzyaaid1(xyzzyaaaa39,xyzzyaaab39)+xyzzyaaib1(xyzzyaaaa39&
&,xyzzyaaab39)*xyzzyaahx1(xyzzyaaab39)
xyzzyaaad39=xyzzyaaie1(xyzzyaaaa39,xyzzyaaab39)+xyzzyaaic1(xyzzyaaaa39&
&,xyzzyaaab39)*xyzzyaahx1(xyzzyaaab39)
if(xyzzyaahy1(xyzzyaaab39)==0.d0)then
xyzzyaaid1(xyzzyaaaa39,xyzzyaaab39)=0.d0
xyzzyaaic1(xyzzyaaaa39,xyzzyaaab39)=0.d0
else
xyzzyaaid1(xyzzyaaaa39,xyzzyaaab39)=xyzzyaaac39/xyzzyaahy1(xyzzyaaab39&
&)
xyzzyaaie1(xyzzyaaaa39,xyzzyaaab39)=xyzzyaaad39/xyzzyaahy1(xyzzyaaab39&
&)
endif
enddo
enddo
else
xyzzyaahy1(:)=xyzzyaahx1(:)
endif
write(xyzzyaaaf1,'(a)')'START MOMENTUM DENSITY'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_mom_den
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaahs1))
write(xyzzyaaaf1,'(a)')'Number of random points to sample'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaaht1))
do xyzzyaaab39=1,xyzzyaahs1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab39))
write(xyzzyaaaf1,'(a)')'Number of particle types in set'
if(on_top_ii==0)then
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaahu1(xyzzyaaab39)))
char_80=''
do xyzzyaaaa39=1,xyzzyaahu1(xyzzyaaab39)
char_80=trim(adjustl(char_80))//' '//trim(i2s(xyzzyaahv1(xyzzyaaaa39,x&
&yzzyaaab39)))
enddo
else
write(xyzzyaaaf1,'(1x,a)')'on-top pair'
char_80=' '//trim(i2s(which_spin(on_top_ii)))//' '//trim(i2s(which_spi&
&n(on_top_jj)))
endif
write(xyzzyaaaf1,'(a)')'Particle types'
write(xyzzyaaaf1,'(1x,a)')trim(char_80)
write(xyzzyaaaf1,'(a)')'Weight for this set'
write(xyzzyaaaf1,*)xyzzyaahy1(xyzzyaaab39)
write(xyzzyaaaf1,'(a)')'|k|/k_F,MOM_DEN(k)**2,MOM_DEN(k)'
if(mom_den)then
do xyzzyaaaa39=1,expval_ngvec(1)
xyzzyaaae39=sqrt(ddot(3,xyzzyaaij1(1,xyzzyaaaa39),1,xyzzyaaij1(1,xyzzy&
&aaaa39),1))*xyzzyaaaf39
write(xyzzyaaaf1,*)xyzzyaaae39,xyzzyaaie1(xyzzyaaaa39,xyzzyaaab39),xyz&
&zyaaid1(xyzzyaaaa39,xyzzyaaab39)
enddo
else
do xyzzyaaaa39=1,expval_ngvec(1)
xyzzyaaae39=sqrt(ddot(3,xyzzyaaij1(1,xyzzyaaaa39),1,xyzzyaaij1(1,xyzzy&
&aaaa39),1))*xyzzyaaaf39
write(xyzzyaaaf1,*)xyzzyaaae39,xyzzyaaic1(xyzzyaaaa39,xyzzyaaab39),xyz&
&zyaaib1(xyzzyaaaa39,xyzzyaaab39)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab39))
enddo
write(xyzzyaaaf1,'(a)')'END MOMENTUM DENSITY'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqp1
subroutine xyzzyaaqq1
implicit none
integer xyzzyaaaa40,xyzzyaaab40
real(dp) xyzzyaaac40,xyzzyaaad40,xyzzyaaae40,xyzzyaaaf40
character(80) char_80,tmpr
xyzzyaaaf40=1.d0
if(twop_dm_mom)then
xyzzyaanp1=0.d0
xyzzyaanq1=0.d0
call mpi_reduce(xyzzyaanl1,xyzzyaanp1,expval_ngvec(1)*xyzzyaanc1,mpi_d&
&ouble_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_twop_dm_mom in write_twop_d&
&m_mom.')
call mpi_reduce(xyzzyaanm1,xyzzyaanq1,expval_ngvec(1)*xyzzyaanc1,mpi_d&
&ouble_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_twop_dm_mom2 in write_twop_&
&dm_mom.')
call mpi_reduce(xyzzyaani1,xyzzyaank1,xyzzyaanc1,mpi_double_precision,&
&mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_twop_dm_mom_weight in write&
&_twop_dm_mom.')
endif
if(am_master)then
if(twop_dm_mom)then
do xyzzyaaab40=1,xyzzyaanc1
xyzzyaank1(xyzzyaaab40)=xyzzyaanj1(xyzzyaaab40)+xyzzyaank1(xyzzyaaab40&
&)
do xyzzyaaaa40=1,expval_ngvec(1)
xyzzyaaac40=xyzzyaanp1(xyzzyaaaa40,xyzzyaaab40)+xyzzyaann1(xyzzyaaaa40&
&,xyzzyaaab40)*xyzzyaanj1(xyzzyaaab40)
xyzzyaaad40=xyzzyaanq1(xyzzyaaaa40,xyzzyaaab40)+xyzzyaano1(xyzzyaaaa40&
&,xyzzyaaab40)*xyzzyaanj1(xyzzyaaab40)
if(xyzzyaank1(xyzzyaaab40)==0.d0)then
xyzzyaanp1(xyzzyaaaa40,xyzzyaaab40)=0.d0
xyzzyaano1(xyzzyaaaa40,xyzzyaaab40)=0.d0
else
xyzzyaanp1(xyzzyaaaa40,xyzzyaaab40)=xyzzyaaac40/xyzzyaank1(xyzzyaaab40&
&)
xyzzyaanq1(xyzzyaaaa40,xyzzyaaab40)=xyzzyaaad40/xyzzyaank1(xyzzyaaab40&
&)
endif
enddo
enddo
else
xyzzyaank1(:)=xyzzyaanj1(:)
endif
write(xyzzyaaaf1,'(a)')'START TWO-PARTICLE DENSITY MATRIX (MOMENTUM SP&
&ACE)'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_twop_dm_mom
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaanc1))
write(xyzzyaaaf1,'(a)')'Number of random points to sample'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaand1))
write(xyzzyaaaf1,'(a)')'Fraction of particle pairs to sample at random&
&'
tmpr=r2s(xyzzyaanh1,'(f16.4)')
write(xyzzyaaaf1,'(1x,a)')trim(tmpr)
do xyzzyaaab40=1,xyzzyaanc1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab40))
write(xyzzyaaaf1,'(a)')'Number of particle-pair types in set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaane1(xyzzyaaab40)))
char_80=''
do xyzzyaaaa40=1,xyzzyaane1(xyzzyaaab40)
char_80=trim(adjustl(char_80))//'  '//trim(i2s(xyzzyaanf1(1,xyzzyaaaa4&
&0,xyzzyaaab40)))//' '//trim(i2s(xyzzyaanf1(2,xyzzyaaaa40,xyzzyaaab40)&
&))
enddo
write(xyzzyaaaf1,'(a)')'Particle-pair types'
write(xyzzyaaaf1,'(1x,a)')trim(adjustl(char_80))
write(xyzzyaaaf1,'(a)')'Weight for this set'
write(xyzzyaaaf1,*)xyzzyaank1(xyzzyaaab40)
write(xyzzyaaaf1,'(a)')'|k|/k_F,twop_dm_mom(k)**2,twop_dm_mom(k)'
if(twop_dm_mom)then
do xyzzyaaaa40=1,expval_ngvec(1)
xyzzyaaae40=sqrt(ddot(3,xyzzyaanu1(1,xyzzyaaaa40),1,xyzzyaanu1(1,xyzzy&
&aaaa40),1))*xyzzyaaaf40
write(xyzzyaaaf1,*)xyzzyaaae40,xyzzyaanq1(xyzzyaaaa40,xyzzyaaab40),xyz&
&zyaanp1(xyzzyaaaa40,xyzzyaaab40)
enddo
else
do xyzzyaaaa40=1,expval_ngvec(1)
xyzzyaaae40=sqrt(ddot(3,xyzzyaanu1(1,xyzzyaaaa40),1,xyzzyaanu1(1,xyzzy&
&aaaa40),1))*xyzzyaaaf40
write(xyzzyaaaf1,*)xyzzyaaae40,xyzzyaano1(xyzzyaaaa40,xyzzyaaab40),xyz&
&zyaann1(xyzzyaaaa40,xyzzyaaab40)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab40))
enddo
write(xyzzyaaaf1,'(a)')'END TWO-PARTICLE DENSITY MATRIX (MOMENTUM SPAC&
&E)'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqq1
subroutine xyzzyaaqr1
implicit none
integer xyzzyaaaa41,xyzzyaaab41
real(dp) xyzzyaaac41,xyzzyaaad41
if(finite_density)then
xyzzyaali1=0.d0
xyzzyaalj1=0.d0
xyzzyaalq1=0.d0
xyzzyaalr1=0.d0
call mpi_reduce(xyzzyaale1,xyzzyaali1,xyzzyaakk1*xyzzyaaki1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_fin_den_ij in write_finite_&
&density.')
call mpi_reduce(xyzzyaalf1,xyzzyaalj1,xyzzyaakk1*xyzzyaaki1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_fin_den2_ij in write_finite&
&_density.')
call mpi_reduce(xyzzyaalm1,xyzzyaalq1,xyzzyaakk1*xyzzyaakj1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_fin_den_i in write_finite_d&
&ensity.')
call mpi_reduce(xyzzyaaln1,xyzzyaalr1,xyzzyaakk1*xyzzyaakj1,mpi_double&
&_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_fin_den2_i in write_finite_&
&density.')
call mpi_reduce(xyzzyaaks1,xyzzyaaku1,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_fin_den_nstep in write_fini&
&te_density.')
call mpi_reduce(xyzzyaakw1,xyzzyaaky1,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_fin_den_weight in write_fin&
&ite_density.')
call mpi_reduce(xyzzyaala1,xyzzyaalc1,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_fin_den_weight2 in write_fi&
&nite_density.')
endif
if(am_master)then
if(finite_density)then
xyzzyaaku1=xyzzyaakt1+xyzzyaaku1
xyzzyaaky1=xyzzyaakx1+xyzzyaaky1
xyzzyaalc1=xyzzyaalb1+xyzzyaalc1
do xyzzyaaab41=1,xyzzyaaki1
do xyzzyaaaa41=1,xyzzyaakk1
xyzzyaaac41=xyzzyaali1(xyzzyaaaa41,xyzzyaaab41)+xyzzyaalg1(xyzzyaaaa41&
&,xyzzyaaab41)*xyzzyaakx1
xyzzyaaad41=xyzzyaalj1(xyzzyaaaa41,xyzzyaaab41)+xyzzyaalh1(xyzzyaaaa41&
&,xyzzyaaab41)*xyzzyaakx1
if(xyzzyaaky1==0.d0)then
xyzzyaali1(xyzzyaaaa41,xyzzyaaab41)=0.d0
xyzzyaalj1(xyzzyaaaa41,xyzzyaaab41)=0.d0
else
xyzzyaali1(xyzzyaaaa41,xyzzyaaab41)=xyzzyaaac41/xyzzyaaky1
xyzzyaalj1(xyzzyaaaa41,xyzzyaaab41)=xyzzyaaad41/xyzzyaaky1
endif
enddo
enddo
do xyzzyaaab41=1,xyzzyaakj1
do xyzzyaaaa41=1,xyzzyaakk1
xyzzyaaac41=xyzzyaalq1(xyzzyaaaa41,xyzzyaaab41)+xyzzyaalo1(xyzzyaaaa41&
&,xyzzyaaab41)*xyzzyaakx1
xyzzyaaad41=xyzzyaalr1(xyzzyaaaa41,xyzzyaaab41)+xyzzyaalp1(xyzzyaaaa41&
&,xyzzyaaab41)*xyzzyaakx1
if(xyzzyaaky1==0.d0)then
xyzzyaalq1(xyzzyaaaa41,xyzzyaaab41)=0.d0
xyzzyaalr1(xyzzyaaaa41,xyzzyaaab41)=0.d0
else
xyzzyaalq1(xyzzyaaaa41,xyzzyaaab41)=xyzzyaaac41/xyzzyaaky1
xyzzyaalr1(xyzzyaaaa41,xyzzyaaab41)=xyzzyaaad41/xyzzyaaky1
endif
enddo
enddo
endif
write(xyzzyaaaf1,'(a)')'START FINITE DENSITY'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_fin_den
write(xyzzyaaaf1,'(a)')'Basis functions'
write(xyzzyaaaf1,'(1x,a)')fin_den_basis
if(xyzzyaakl1==xyzzyaakn1.or.xyzzyaakl1==xyzzyaako1.or.xyzzyaakl1==xyz&
&zyaakp1)then
call errstop('WRITE_FINITE_DENSITY','Full functionality is not availab&
&le for the Bessel function, Chebyshev polynomial or exponential funct&
&ion bases.')
endif
write(xyzzyaaaf1,'(a)')'Expansion order'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaakk1))
select case(xyzzyaakl1)
case(xyzzyaakm1,xyzzyaako1)
write(xyzzyaaaf1,'(a)')'Cutoff'
write(xyzzyaaaf1,*)xyzzyaakq1
case(xyzzyaakp1)
write(xyzzyaaaf1,'(a)')'Tail Offset'
write(xyzzyaaaf1,*)xyzzyaakr1
case(xyzzyaakn1)
write(xyzzyaaaf1,'(a)')'Cutoff, Tail offset'
write(xyzzyaaaf1,*)xyzzyaakq1,xyzzyaakr1
end select
write(xyzzyaaaf1,'(a)')'Nstep, Total weight, Total weight^2'
write(xyzzyaaaf1,*)xyzzyaaku1,xyzzyaaky1,xyzzyaalc1
loop_sets_ij: do xyzzyaaab41=1,xyzzyaaki1
write(xyzzyaaaf1,'(a)')'START e-e SET '//trim(i2s(xyzzyaaab41))
write(xyzzyaaaf1,'(a)')'n_ij,(n_ij)**2'
if(finite_density)then
do xyzzyaaaa41=1,xyzzyaakk1
write(xyzzyaaaf1,*)xyzzyaali1(xyzzyaaaa41,xyzzyaaab41),xyzzyaalj1(xyzz&
&yaaaa41,xyzzyaaab41)
enddo
else
do xyzzyaaaa41=1,xyzzyaakk1
write(xyzzyaaaf1,*)xyzzyaalg1(xyzzyaaaa41,xyzzyaaab41),xyzzyaalh1(xyzz&
&yaaaa41,xyzzyaaab41)
enddo
endif
write(xyzzyaaaf1,'(a)')'END e-e SET '//trim(i2s(xyzzyaaab41))
enddo loop_sets_ij
loop_sets_i: do xyzzyaaab41=1,xyzzyaakj1
write(xyzzyaaaf1,'(a)')'START e-nucleus SET '//trim(i2s(xyzzyaaab41))
write(xyzzyaaaf1,'(a)')'n_i,(n_i)**2'
if(finite_density)then
do xyzzyaaaa41=1,xyzzyaakk1
write(xyzzyaaaf1,*)xyzzyaalq1(xyzzyaaaa41,xyzzyaaab41),xyzzyaalr1(xyzz&
&yaaaa41,xyzzyaaab41)
enddo
else
do xyzzyaaaa41=1,xyzzyaakk1
write(xyzzyaaaf1,*)xyzzyaalo1(xyzzyaaaa41,xyzzyaaab41),xyzzyaalp1(xyzz&
&yaaaa41,xyzzyaaab41)
enddo
endif
write(xyzzyaaaf1,'(a)')'END e-nucleus SET '//trim(i2s(xyzzyaaab41))
enddo loop_sets_i
write(xyzzyaaaf1,'(a)')'END FINITE DENSITY'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqr1
subroutine xyzzyaaqs1
implicit none
integer xyzzyaaaa42,xyzzyaaab42
real(dp) xyzzyaaac42,xyzzyaaad42
character(80) char_80
if(population)then
xyzzyaamx1=0.d0
xyzzyaamy1=0.d0
call mpi_reduce(xyzzyaamt1,xyzzyaamx1,nitot*xyzzyaamn1,mpi_double_prec&
&ision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_population in WRITE_POPULAT&
&ION.')
call mpi_reduce(xyzzyaamu1,xyzzyaamy1,nitot*xyzzyaamn1,mpi_double_prec&
&ision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_population2 in WRITE_POPULA&
&TION.')
call mpi_reduce(xyzzyaamq1,xyzzyaams1,xyzzyaamn1,mpi_double_precision,&
&mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_population_weight in WRITE_&
&POPULATION.')
endif
if(am_master)then
if(population)then
do xyzzyaaab42=1,xyzzyaamn1
xyzzyaams1(xyzzyaaab42)=xyzzyaamr1(xyzzyaaab42)+xyzzyaams1(xyzzyaaab42&
&)
do xyzzyaaaa42=1,nitot
xyzzyaaac42=xyzzyaamx1(xyzzyaaaa42,xyzzyaaab42)+xyzzyaamv1(xyzzyaaaa42&
&,xyzzyaaab42)*xyzzyaamr1(xyzzyaaab42)
xyzzyaaad42=xyzzyaamy1(xyzzyaaaa42,xyzzyaaab42)+xyzzyaamw1(xyzzyaaaa42&
&,xyzzyaaab42)*xyzzyaamr1(xyzzyaaab42)
if(xyzzyaams1(xyzzyaaab42)==0.d0)then
xyzzyaamx1(xyzzyaaaa42,xyzzyaaab42)=0.d0
xyzzyaamw1(xyzzyaaaa42,xyzzyaaab42)=0.d0
else
xyzzyaamx1(xyzzyaaaa42,xyzzyaaab42)=xyzzyaaac42/xyzzyaams1(xyzzyaaab42&
&)
xyzzyaamy1(xyzzyaaaa42,xyzzyaaab42)=xyzzyaaad42/xyzzyaams1(xyzzyaaab42&
&)
endif
enddo
enddo
else
xyzzyaams1(:)=xyzzyaamr1(:)
endif
write(xyzzyaaaf1,'(a)')'START POPULATION'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_population
write(xyzzyaaaf1,'(a)')'Number of sets'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaamn1))
write(xyzzyaaaf1,'(a)')'Number of ions'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(nitot))
do xyzzyaaab42=1,xyzzyaamn1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaab42))
write(xyzzyaaaf1,'(a)')'Number of particle types in set'
write(xyzzyaaaf1,'(1x,a)')trim(i2s(xyzzyaamo1(xyzzyaaab42)))
char_80=''
do xyzzyaaaa42=1,xyzzyaamo1(xyzzyaaab42)
char_80=trim(adjustl(char_80))//' ' //trim(i2s(xyzzyaamp1(xyzzyaaaa42,&
&xyzzyaaab42)))
enddo
write(xyzzyaaaf1,'(a)')'Particle types'
write(xyzzyaaaf1,'(1x,a)')trim(char_80)
write(xyzzyaaaf1,'(a)')'Weight for this set'
write(xyzzyaaaf1,*)xyzzyaams1(xyzzyaaab42)
write(xyzzyaaaf1,'(a)')'ion, population^2, population'
if(population)then
do xyzzyaaaa42=1,nitot
write(xyzzyaaaf1,*)xyzzyaaaa42,xyzzyaamy1(xyzzyaaaa42,xyzzyaaab42),xyz&
&zyaamx1(xyzzyaaaa42,xyzzyaaab42)
enddo
else
do xyzzyaaaa42=1,nitot
write(xyzzyaaaf1,*)xyzzyaaaa42,xyzzyaamw1(xyzzyaaaa42,xyzzyaaab42),xyz&
&zyaamv1(xyzzyaaaa42,xyzzyaaab42)
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaab42))
enddo
write(xyzzyaaaf1,'(a)')'END POPULATION'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqs1
subroutine xyzzyaaqt1
implicit none
integer xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43,xyzzyaaad43
real(dp) xyzzyaaae43
if(mol_density)then
xyzzyaamj1=0.d0
call mpi_reduce(xyzzyaamh1,xyzzyaamj1,product(xyzzyaalw1)*xyzzyaalv1,m&
&pi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_mol_den in write_mol_densit&
&y.')
call mpi_reduce(xyzzyaalz1,xyzzyaamb1,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_mol_den_nstep in write_mol_&
&density.')
call mpi_reduce(xyzzyaamd1,xyzzyaamf1,1,mpi_double_precision,mpi_sum,0&
&,mpi_comm_world,ierror)
call checkmpi(ierror,'mpi_reduce of expval_mol_den_weight in write_mol&
&_density.')
endif
if(nnodes==1)then
xyzzyaamj1=xyzzyaamh1
xyzzyaamf1=xyzzyaamd1
endif
if(am_master)then
if(mol_density)then
xyzzyaamb1=xyzzyaama1+xyzzyaamb1
xyzzyaamf1=xyzzyaame1+xyzzyaamf1
do xyzzyaaad43=1,xyzzyaalv1
do xyzzyaaaa43=1,xyzzyaalw1(1)
do xyzzyaaab43=1,xyzzyaalw1(2)
do xyzzyaaac43=1,xyzzyaalw1(3)
xyzzyaaae43=xyzzyaamj1(xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43,xyzzyaaad43&
&)+xyzzyaami1(xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43,xyzzyaaad43)*xyzzyaa&
&me1
if(xyzzyaamf1==0.d0)then
xyzzyaamj1(xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43,xyzzyaaad43)=0.d0
else
xyzzyaamj1(xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43,xyzzyaaad43)=xyzzyaaae4&
&3/xyzzyaamf1
endif
enddo
enddo
enddo
enddo
endif
write(xyzzyaaaf1,'(a)')'START MOLECULAR DENSITY'
write(xyzzyaaaf1,'(a)')'Accumulation carried out using'
write(xyzzyaaaf1,'(1x,a)')accumulation_method_mol_den
write(xyzzyaaaf1,'(a)')'Grid size'
write(xyzzyaaaf1,*)trim(i2s(xyzzyaalw1(1))),' ',trim(i2s(xyzzyaalw1(2)&
&)),' ',trim(i2s(xyzzyaalw1(3)))
write(xyzzyaaaf1,'(a)')'Coordinates of A (max corner)'
write(xyzzyaaaf1,*)xyzzyaaly1(1:3)
write(xyzzyaaaf1,'(a)')'Coordinates of B (min corner)'
write(xyzzyaaaf1,*)xyzzyaalx1(1:3)
write(xyzzyaaaf1,'(a)')'Nstep, Total weight'
write(xyzzyaaaf1,*)xyzzyaamb1,xyzzyaamf1
loop_sets_i: do xyzzyaaad43=1,xyzzyaalv1
write(xyzzyaaaf1,'(a)')'START SET '//trim(i2s(xyzzyaaad43))
write(xyzzyaaaf1,'(a)')'n_i'
if(mol_density)then
do xyzzyaaac43=1,xyzzyaalw1(3)
do xyzzyaaab43=1,xyzzyaalw1(2)
do xyzzyaaaa43=1,xyzzyaalw1(1)
write(xyzzyaaaf1,*)xyzzyaamj1(xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43,xyzz&
&yaaad43)
enddo
enddo
enddo
else
do xyzzyaaac43=1,xyzzyaalw1(3)
do xyzzyaaab43=1,xyzzyaalw1(2)
do xyzzyaaaa43=1,xyzzyaalw1(1)
write(xyzzyaaaf1,*)xyzzyaami1(xyzzyaaaa43,xyzzyaaab43,xyzzyaaac43,xyzz&
&yaaad43)
enddo
enddo
enddo
endif
write(xyzzyaaaf1,'(a)')'END SET '//trim(i2s(xyzzyaaad43))
enddo loop_sets_i
write(xyzzyaaaf1,'(a)')'END MOLECULAR DENSITY'
write(xyzzyaaaf1,*)
endif
end subroutine xyzzyaaqt1
subroutine setup_expval(from_backup)
implicit none
logical,intent(in) :: from_backup
integer xyzzyaaaa44
xyzzyaaam1=2.d0*expval_cutoff
homogeneous_density=(trim(atom_basis_type)=='none')
if(homogeneous_density)then
do xyzzyaaaa44=1,nspin
if(nele(xyzzyaaaa44)>0.and.any(heg_orbtype(xyzzyaaaa44,:)==0.or.heg_or&
&btype(xyzzyaaaa44,:)>=100.or.heg_orbtype(xyzzyaaaa44,:)==2))homogeneo&
&us_density=.false.
enddo
endif
xyzzyaabb1=periodicity/=0.and.(trim(interaction)=='mpc'.or.trim(intera&
&ction)=='ewald_mpc'.or.trim(interaction)=='mpc_ewald').and.qmc_densit&
&y_mpc
if(density.or.spin_density.or.spin_density_mat.or.pair_corr.or.structu&
&re_factor.or.mom_den.or.xyzzyaabb1.or.twop_dm_mom.or.cond_frac_mom)th&
&en
xyzzyaaba1=.true.
else
xyzzyaaba1=.false.
endif
if(am_master)then
xyzzyaaax1(:)=.false.
if(xyzzyaabc1.or.density)xyzzyaaax1(1)=.true.
if(xyzzyaabd1.or.spin_density)xyzzyaaax1(2)=.true.
if(xyzzyaabe1.or.spin_density_mat)xyzzyaaax1(3)=.true.
if(pair_corr)xyzzyaaax1(2)=.true.
if(xyzzyaabl1)xyzzyaaax1(4)=.true.
if(xyzzyaabn1.or.pair_corr_sph)xyzzyaaax1(5)=.true.
if(xyzzyaabf1.or.loc_tensor)xyzzyaaax1(6)=.true.
if(xyzzyaabh1.or.structure_factor)xyzzyaaax1(7)=.true.
if(xyzzyaabj1.or.structure_factor_sph)xyzzyaaax1(8)=.true.
if(xyzzyaabp1.or.onep_density_mat)xyzzyaaax1(9)=.true.
if(xyzzyaabr1.or.twop_density_mat)xyzzyaaax1(10)=.true.
if(xyzzyaabt1.or.cond_fraction)xyzzyaaax1(11)=.true.
if(xyzzyaabv1.or.mom_den)xyzzyaaax1(12)=.true.
if(xyzzyaabx1.or.finite_density)xyzzyaaax1(13)=.true.
if(xyzzyaabz1.or.population)xyzzyaaax1(14)=.true.
if(xyzzyaacb1.or.mol_density)xyzzyaaax1(15)=.true.
if(xyzzyaacd1.or.twop_dm_mom)xyzzyaaax1(16)=.true.
if(xyzzyaacf1.or.cond_frac_mom)xyzzyaaax1(17)=.true.
if(all(.not.xyzzyaaax1(:)))call errstop('SETUP_EXPVAL','No expectation&
& values flagged for writing. Should not happen.')
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaaax1,xyzzyaaal1,mpi_logical,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting expval_write in setup_expval.')
call mpi_bcast(xyzzyaaaw1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting new_expvals in setup_expval.')
call mpi_bcast(xyzzyaaav1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_present in setup_expval.')
call mpi_bcast(xyzzyaaat1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_nonempty in setup_expval.')
if(am_slave)allocate(xyzzyaaco1(no_families),stat=xyzzyaaad1)
call mpi_bcast(xyzzyaaco1,no_families,mpi_integer,0,mpi_comm_world,ier&
&ror)
call checkmpi(ierror,'Broadcasting nparticles in setup_expval.')
endif
if(.not.xyzzyaaav1.or..not.xyzzyaaat1)then
title='No title given'
xyzzyaack1=xyzzyaaaj1
xyzzyaacl1=nspin
allocate(xyzzyaacn1(xyzzyaacl1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_EXPVAL','nele_expval')
xyzzyaacn1(:)=nele(:)
xyzzyaacm1=dimensionality
xyzzyaaci1=periodicity
select case(xyzzyaaci1)
case(1)
xyzzyaacq1=a1(1)
case(2)
xyzzyaacp1=area
case(3)
xyzzyaacr1=volume
end select
xyzzyaacj1=scell_matrix
xyzzyaacw1=0
if(xyzzyaaba1)then
if(density.or.spin_density.or.spin_density_mat)then
if(structure_factor.or.pair_corr.or.(from_backup.and.turn_on_pair_corr&
&))then
xyzzyaacw1=2
else
xyzzyaacw1=1
endif
else
xyzzyaacw1=1
endif
allocate(xyzzyaadc1(3,3,xyzzyaacw1),xyzzyaadd1(3,3,xyzzyaacw1),stat=xy&
&zzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_EXPVAL','expval_latvec,expval_rlatv&
&ec')
if(structure_factor.or.mom_den.or.pair_corr.or.twop_dm_mom.or.cond_fra&
&c_mom)then
xyzzyaadc1(:,1,1)=a1(:)
xyzzyaadc1(:,2,1)=a2(:)
xyzzyaadc1(:,3,1)=a3(:)
xyzzyaadd1(:,1,1)=b1(:)
xyzzyaadd1(:,2,1)=b2(:)
xyzzyaadd1(:,3,1)=b3(:)
call xyzzyaaqu1(1,from_backup)
endif
if(density.or.spin_density.or.spin_density_mat)then
xyzzyaadc1(:,1,xyzzyaacw1)=pa1(:)
xyzzyaadc1(:,2,xyzzyaacw1)=pa2(:)
xyzzyaadc1(:,3,xyzzyaacw1)=pa3(:)
xyzzyaadd1(:,1,xyzzyaacw1)=pb1(:)
xyzzyaadd1(:,2,xyzzyaacw1)=pb2(:)
xyzzyaadd1(:,3,xyzzyaacw1)=pb3(:)
call xyzzyaaqu1(xyzzyaacw1,from_backup)
endif
endif
else
if(nnodes>1)then
call mpi_bcast(xyzzyaaci1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_periodicity in setup_expval.&
&')
select case(xyzzyaaci1)
case(1)
call mpi_bcast(xyzzyaacq1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting simcell_length in setup_expval.')
case(2)
call mpi_bcast(xyzzyaacp1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting simcell_area in setup_expval.')
case(3)
call mpi_bcast(xyzzyaacr1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'Broadcasting simcell_volume in setup_expval.')
end select
endif
if(xyzzyaaba1)then
if(nnodes>1)then
call mpi_bcast(xyzzyaacv1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_ngsets_in_file in setup_expv&
&al.')
call mpi_bcast(xyzzyaacw1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_ngsets in setup_expval.')
endif
if(.not.xyzzyaaaw1)then
if(nnodes>1)then
if(am_slave)then
allocate(expval_ngvec(xyzzyaacw1),xyzzyaadc1(3,3,xyzzyaacw1),xyzzyaadd&
&1(3,3,xyzzyaacw1),stat=xyzzyaaad1)
if(xyzzyaabb1)allocate(expval_ngvec_truncated(xyzzyaacw1),stat=xyzzyaa&
&ad1)
endif
call mpi_bcast(expval_ngvec,xyzzyaacw1,mpi_integer,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting expval_ngvec in setup_expval.')
call mpi_bcast(xyzzyaacu1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_ngvec_max in setup_expval.')
if(am_slave)then
allocate(expval_gvec(3,xyzzyaacu1,xyzzyaacw1),xyzzyaadc1(3,3,xyzzyaacw&
&1),xyzzyaadd1(3,3,xyzzyaacw1),stat=xyzzyaaad1)
endif
call mpi_bcast(xyzzyaadc1,3*3*xyzzyaacw1,mpi_double_precision,0,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_latvec in setup_expval.')
call mpi_bcast(xyzzyaadd1,3*3*xyzzyaacw1,mpi_double_precision,0,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_rlatvec in setup_expval.')
call mpi_bcast(expval_gvec,3*xyzzyaacu1*xyzzyaacw1,mpi_double_precisio&
&n,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_gvec in setup_expval.')
endif
else
if(xyzzyaacv1>0)then
if(am_master)then
if(xyzzyaaaz1)then
xyzzyaadc1(:,1,xyzzyaacv1+1)=a1(:)
xyzzyaadc1(:,2,xyzzyaacv1+1)=a2(:)
xyzzyaadc1(:,3,xyzzyaacv1+1)=a3(:)
xyzzyaadd1(:,1,xyzzyaacv1+1)=b1(:)
xyzzyaadd1(:,2,xyzzyaacv1+1)=b2(:)
xyzzyaadd1(:,3,xyzzyaacv1+1)=b3(:)
call xyzzyaaqu1(xyzzyaacv1+1,from_backup)
xyzzyaadb1(xyzzyaacv1+1)=expval_cutoff
endif
if(xyzzyaaay1)then
xyzzyaadc1(:,1,xyzzyaacw1)=pa1(:)
xyzzyaadc1(:,2,xyzzyaacw1)=pa2(:)
xyzzyaadc1(:,3,xyzzyaacw1)=pa3(:)
xyzzyaadd1(:,1,xyzzyaacw1)=pb1(:)
xyzzyaadd1(:,2,xyzzyaacw1)=pb2(:)
xyzzyaadd1(:,3,xyzzyaacw1)=pb3(:)
call xyzzyaaqu1(xyzzyaacw1,from_backup)
xyzzyaadb1(xyzzyaacw1)=expval_cutoff
endif
xyzzyaacu1=maxval(expval_ngvec(:))
endif
if(nnodes>1)then
if(am_slave)then
allocate(expval_ngvec(xyzzyaacw1),xyzzyaadc1(3,3,xyzzyaacw1),xyzzyaadd&
&1(3,3,xyzzyaacw1),stat=xyzzyaaad1)
if(xyzzyaabb1)allocate(expval_ngvec_truncated(xyzzyaacw1),stat=xyzzyaa&
&ad1)
endif
call mpi_bcast(expval_ngvec,xyzzyaacw1,mpi_integer,0,mpi_comm_world,ie&
&rror)
call checkmpi(ierror,'Broadcasting expval_ngvec in setup_expval.')
if(am_slave)then
xyzzyaacu1=maxval(expval_ngvec(:))
allocate(expval_gvec(3,xyzzyaacu1,xyzzyaacw1),xyzzyaadc1(3,3,xyzzyaacw&
&1),xyzzyaadd1(3,3,xyzzyaacw1),stat=xyzzyaaad1)
endif
call mpi_bcast(xyzzyaadc1,3*3*xyzzyaacw1,mpi_double_precision,0,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_latvec in setup_expval.')
call mpi_bcast(xyzzyaadd1,3*3*xyzzyaacw1,mpi_double_precision,0,mpi_co&
&mm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_rlatvec in setup_expval.')
call mpi_bcast(expval_gvec,3*xyzzyaacu1*xyzzyaacw1,mpi_double_precisio&
&n,0,mpi_comm_world,ierror)
call checkmpi(ierror,'Broadcasting expval_gvec in setup_expval.')
endif
else
if(density.or.spin_density.or.spin_density_mat)then
if(structure_factor.or.pair_corr)then
xyzzyaacw1=2
else
xyzzyaacw1=1
endif
else
xyzzyaacw1=1
endif
allocate(xyzzyaadc1(3,3,xyzzyaacw1),xyzzyaadd1(3,3,xyzzyaacw1),stat=xy&
&zzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_EXPVAL','expval_latvec,expval_rlatv&
&ec')
if(structure_factor.or.mom_den.or.pair_corr.or.twop_dm_mom.or.cond_fra&
&c_mom)then
xyzzyaadc1(:,1,1)=a1(:)
xyzzyaadc1(:,2,1)=a2(:)
xyzzyaadc1(:,3,1)=a3(:)
xyzzyaadd1(:,1,1)=b1(:)
xyzzyaadd1(:,2,1)=b2(:)
xyzzyaadd1(:,3,1)=b3(:)
call xyzzyaaqu1(1,from_backup)
endif
if(density.or.spin_density.or.spin_density_mat)then
xyzzyaadc1(:,1,xyzzyaacw1)=pa1(:)
xyzzyaadc1(:,2,xyzzyaacw1)=pa2(:)
xyzzyaadc1(:,3,xyzzyaacw1)=pa3(:)
xyzzyaadd1(:,1,xyzzyaacw1)=pb1(:)
xyzzyaadd1(:,2,xyzzyaacw1)=pb2(:)
xyzzyaadd1(:,3,xyzzyaacw1)=pb3(:)
call xyzzyaaqu1(xyzzyaacw1,from_backup)
endif
endif
endif
endif
endif
if(xyzzyaaba1)call xyzzyaaqw1
if(density)call xyzzyaaqx1
if(spin_density.and..not.(from_backup.and.turn_on_pair_corr))call xyzz&
&yaaqy1
if(spin_density_mat)call xyzzyaaqz1
if(pair_corr.or.(from_backup.and.turn_on_pair_corr))then
spin_density=.true.
if(turn_on_pair_corr.and.pair_corr)spin_density=.false.
call xyzzyaaqy1
call xyzzyaara1(from_backup)
endif
if(pair_corr_sph)call xyzzyaarb1
if(loc_tensor)call xyzzyaarc1
if(structure_factor)call xyzzyaard1
if(structure_factor_sph)call xyzzyaare1
if(onep_density_mat)call xyzzyaarf1
if(twop_density_mat)call xyzzyaarg1
if(cond_fraction)call xyzzyaarh1
if(mom_den)call xyzzyaari1
if(finite_density)call xyzzyaarj1
if(population)call xyzzyaark1
if(mol_density)call xyzzyaarl1
if(twop_dm_mom)call xyzzyaarm1
if(cond_frac_mom)call xyzzyaarn1
if(am_master.and..not.from_backup)then
if(isperiodic.and.xyzzyaaba1.and.inversion_symmetry)then
call wout()
call wout('Assuming real Fourier coefficients since structure has inve&
&rsion symmetry.')
endif
call wout()
endif
end subroutine setup_expval
subroutine xyzzyaaqu1(set,from_backup)
implicit none
integer,intent(in) :: set
logical,intent(in) :: from_backup
integer xyzzyaaaa45,xyzzyaaab45,xyzzyaaac45,xyzzyaaad45,xyzzyaaae45,xy&
&zzyaaaf45,xyzzyaaag45,xyzzyaaah45
integer,allocatable :: xyzzyaaai45(:)
real(dp) xyzzyaaaj45,xyzzyaaak45,xyzzyaaal45,xyzzyaaam45,xyzzyaaan45,x&
&yzzyaaao45,xyzzyaaap45,xyzzyaaaq45,xyzzyaaar45,xyzzyaaas45
real(dp),allocatable :: xyzzyaaat45(:,:),xyzzyaaau45(:)
character(80) tmpr
if(am_master.and..not.from_backup)then
call wout()
tmpr=r2s(expval_cutoff,'(f15.4)')
call wout('Generating new set of G vectors using cutoff of '//trim(tmp&
&r)//' au.')
endif
if(periodicity<=3.and.periodicity>=0)then
call xyzzyaaqv1(xyzzyaacs1)
if(am_master.and..not.from_backup)then
call wout()
call wout('Size of raw grid of G vectors')
call wout('(only the vectors within the cutoff will be used)')
call wout(' n1 = '//trim(i2s(xyzzyaacs1(1))))
if(periodicity>1)call wout(' n2 = '//trim(i2s(xyzzyaacs1(2))))
if(periodicity==3)call wout(' n3 = '//trim(i2s(xyzzyaacs1(3))))
endif
else
call errstop('GENERATE_EXPVAL_G_VECTORS','Periodicity error.')
endif
xyzzyaact1(:)=xyzzyaacs1(:)/2
if(set==1)then
allocate(expval_ngvec(xyzzyaacw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'GENERATE_EXPVAL_G_VECTORS','expval_ngvec'&
&)
endif
xyzzyaaap1(1:3)=xyzzyaadd1(1:3,1,set)
xyzzyaaaq1(1:3)=xyzzyaadd1(1:3,2,set)
xyzzyaaar1(1:3)=xyzzyaadd1(1:3,3,set)
expval_ngvec(set)=xyzzyaasb1(xyzzyaact1,xyzzyaaam1)
xyzzyaaad45=xyzzyaasb1(xyzzyaacs1,xyzzyaaam1)
if(expval_ngvec(set)/=xyzzyaaad45)then
xyzzyaaaj45=0.d0
xyzzyaaak45=max(abs(xyzzyaaam1)*2.d0,0.d0)
if(xyzzyaaak45<=0.d0)then
xyzzyaaac45=-4
else
xyzzyaaac45=nint(log(xyzzyaaak45)/log(10.d0))
endif
xyzzyaaaa45=xyzzyaaac45
do
expval_ngvec(set)=xyzzyaasb1(xyzzyaact1,xyzzyaaak45)
xyzzyaaad45=xyzzyaasb1(xyzzyaacs1,xyzzyaaak45)
if(expval_ngvec(set)/=xyzzyaaad45)exit
xyzzyaaak45=xyzzyaaak45+10.d0**xyzzyaaaa45
xyzzyaaaa45=xyzzyaaaa45+1
enddo
do
if(xyzzyaaak45-xyzzyaaaj45<1.d-7)exit
xyzzyaaal45=0.5d0*(xyzzyaaaj45+xyzzyaaak45)
expval_ngvec(set)=xyzzyaasb1(xyzzyaact1,xyzzyaaal45)
xyzzyaaad45=xyzzyaasb1(xyzzyaacs1,xyzzyaaal45)
if(expval_ngvec(set)/=xyzzyaaad45)then
xyzzyaaak45=xyzzyaaal45
else
xyzzyaaaj45=xyzzyaaal45
endif
enddo
if(am_master)then
if(expval_cutoff>0.d0)then
call wout()
call wout('ERROR : Values chosen for following parameters in expval.f9&
&0 module:')
call wout(' n1 = '//trim(i2s(xyzzyaacs1(1))))
call wout(' n2 = '//trim(i2s(xyzzyaacs1(2))))
call wout(' n3 = '//trim(i2s(xyzzyaacs1(3))))
tmpr=r2s(expval_cutoff,'(f16.4)')
call wout('are too small for cutoff EXPVAL_CUTOFF = '//trim(tmpr)//' a&
&u.')
call wout()
call wordwrap('There are G vectors which are below the cutoff that are&
& not included in the n1*n2*n3 grid. BUG')
call wout()
endif
endif
endif
if(expval_cutoff<=0.d0)then
call wout()
call wout('EXPVAL_CUTOFF parameter set to zero or less in input (this &
&parameter now has no upper limit).')
call wout()
call errstop_master('GENERATE_EXPVAL_G_VECTORS','Quitting.')
endif
if(set==1)then
allocate(expval_gvec(3,expval_ngvec(1),xyzzyaacw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'GENERATE_EXPVAL_G_VECTORS','expval_gvec')
endif
allocate(xyzzyaaat45(3,expval_ngvec(set)),xyzzyaaau45(expval_ngvec(set&
&)),xyzzyaaai45(expval_ngvec(set)),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'GENERATE_EXPVAL_G_VECTORS','gtemp,gmagn,i&
&ndx')
xyzzyaaah45=0
xyzzyaaap1(1:3)=xyzzyaadd1(1:3,1,set)
xyzzyaaaq1(1:3)=xyzzyaadd1(1:3,2,set)
xyzzyaaar1(1:3)=xyzzyaadd1(1:3,3,set)
do xyzzyaaag45=0,xyzzyaacs1(3)-1
if(xyzzyaaag45>xyzzyaact1(3))then
xyzzyaaas45=real(xyzzyaaag45-xyzzyaacs1(3),dp)
else
xyzzyaaas45=real(xyzzyaaag45,dp)
endif
do xyzzyaaaf45=0,xyzzyaacs1(2)-1
if(xyzzyaaaf45>xyzzyaact1(2))then
xyzzyaaar45=real(xyzzyaaaf45-xyzzyaacs1(2),dp)
else
xyzzyaaar45=real(xyzzyaaaf45,dp)
endif
do xyzzyaaae45=0,xyzzyaacs1(1)-1
if(xyzzyaaae45>xyzzyaact1(1))then
xyzzyaaaq45=real(xyzzyaaae45-xyzzyaacs1(1),dp)
else
xyzzyaaaq45=real(xyzzyaaae45,dp)
endif
xyzzyaaam45=xyzzyaaaq45*xyzzyaaap1(1)+xyzzyaaar45*xyzzyaaaq1(1)+xyzzya&
&aas45*xyzzyaaar1(1)
xyzzyaaan45=xyzzyaaaq45*xyzzyaaap1(2)+xyzzyaaar45*xyzzyaaaq1(2)+xyzzya&
&aas45*xyzzyaaar1(2)
xyzzyaaao45=xyzzyaaaq45*xyzzyaaap1(3)+xyzzyaaar45*xyzzyaaaq1(3)+xyzzya&
&aas45*xyzzyaaar1(3)
xyzzyaaap45=xyzzyaaam45**2+xyzzyaaan45**2+xyzzyaaao45**2
if(xyzzyaaap45<=xyzzyaaam1)then
xyzzyaaah45=xyzzyaaah45+1
xyzzyaaat45(1,xyzzyaaah45)=xyzzyaaam45
xyzzyaaat45(2,xyzzyaaah45)=xyzzyaaan45
xyzzyaaat45(3,xyzzyaaah45)=xyzzyaaao45
xyzzyaaau45(xyzzyaaah45)=xyzzyaaap45
endif
enddo
enddo
enddo
call quicksort(expval_ngvec(set),xyzzyaaau45(1),xyzzyaaai45(1))
do xyzzyaaaa45=1,expval_ngvec(set)
xyzzyaaab45=xyzzyaaai45(xyzzyaaaa45)
expval_gvec(1:3,xyzzyaaaa45,set)=xyzzyaaat45(1:3,xyzzyaaab45)
enddo
deallocate(xyzzyaaai45)
do xyzzyaaab45=2,expval_ngvec(set),2
xyzzyaaap1(:)=expval_gvec(1:3,xyzzyaaab45,set)
do xyzzyaaaa45=xyzzyaaab45+1,expval_ngvec(set)
if(abs(xyzzyaaap1(1)+expval_gvec(1,xyzzyaaaa45,set))<1.d-8.and. abs(xy&
&zzyaaap1(2)+expval_gvec(2,xyzzyaaaa45,set))<1.d-8.and. abs(xyzzyaaap1&
&(3)+expval_gvec(3,xyzzyaaaa45,set))<1.d-8)then
xyzzyaaaq1(:)=expval_gvec(1:3,xyzzyaaab45+1,set)
expval_gvec(:,xyzzyaaab45+1,set)=expval_gvec(:,xyzzyaaaa45,set)
expval_gvec(:,xyzzyaaaa45,set)=xyzzyaaaq1(:)
endif
enddo
enddo
deallocate(xyzzyaaat45,xyzzyaaau45)
end subroutine xyzzyaaqu1
subroutine xyzzyaaqv1(ngrid_pts)
implicit none
integer,intent(out) :: ngrid_pts(3)
real(dp) xyzzyaaaa46(3)
xyzzyaaaa46(1)=(sqrt(xyzzyaaam1)*sqrt(a1(1)**2+a1(2)**2+a1(3)**2))
xyzzyaaaa46(2)=(sqrt(xyzzyaaam1)*sqrt(a2(1)**2+a2(2)**2+a2(3)**2))
xyzzyaaaa46(3)=(sqrt(xyzzyaaam1)*sqrt(a3(1)**2+a3(2)**2+a3(3)**2))
xyzzyaaaa46=ceiling(xyzzyaaaa46/pi)
ngrid_pts=int(xyzzyaaaa46)+1
if(dimensionality==2)ngrid_pts(3)=1
if(dimensionality==1)ngrid_pts(2:3)=1
end subroutine
subroutine xyzzyaaqw1
implicit none
integer xyzzyaaaa47,xyzzyaaab47,xyzzyaaac47
real(dp) xyzzyaaad47,xyzzyaaae47
if(am_master.and..not.xyzzyaaav1)then
xyzzyaacu1=maxval(expval_ngvec(:))
allocate(xyzzyaacy1(3,xyzzyaacu1,xyzzyaacw1),xyzzyaacz1(xyzzyaacw1),xy&
&zzyaadb1(xyzzyaacw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FOURIER_BASIS','expval_pgmap,grange&
&')
if(xyzzyaabb1)then
allocate(expval_ngvec_truncated(xyzzyaacw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FOURIER_BASIS','expval_ngvec_trunca&
&ted')
endif
allocate(xyzzyaadf1(xyzzyaacu1,xyzzyaacw1,1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FOURIER_BASIS','expval_expigdotr')
xyzzyaadf1(:,:,:)=czero
xyzzyaadb1(:)=expval_cutoff
endif
if(am_slave)then
xyzzyaacu1=maxval(expval_ngvec(:))
allocate(xyzzyaacy1(3,xyzzyaacu1,xyzzyaacw1),xyzzyaacz1(xyzzyaacw1),st&
&at=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FOURIER_BASIS','1')
if(xyzzyaabb1)then
if(.not.allocated(expval_ngvec_truncated))then
allocate(expval_ngvec_truncated(xyzzyaacw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FOURIER_BASIS','2')
endif
endif
allocate(xyzzyaadf1(xyzzyaacu1,xyzzyaacw1,1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FOURIER_BASIS','3')
xyzzyaadf1(:,:,:)=czero
endif
if(xyzzyaabb1.and.am_master)then
do xyzzyaaac47=1,xyzzyaacw1
do xyzzyaaaa47=expval_ngvec(xyzzyaaac47),1,-1
if(abs(expval_den_mpc(xyzzyaaaa47,xyzzyaaac47))>tol_expval_diff)exit
enddo
xyzzyaaad47=sqrt(expval_gvec(1,xyzzyaaaa47,xyzzyaaac47)*expval_gvec(1,&
&xyzzyaaaa47,xyzzyaaac47)+expval_gvec(2,xyzzyaaaa47,xyzzyaaac47)*expva&
&l_gvec(2,xyzzyaaaa47,xyzzyaaac47)+expval_gvec(3,xyzzyaaaa47,xyzzyaaac&
&47)*expval_gvec(3,xyzzyaaaa47,xyzzyaaac47))
expval_ngvec_truncated(xyzzyaaac47)=expval_ngvec(xyzzyaaac47)
do xyzzyaaab47=xyzzyaaaa47+1,expval_ngvec(xyzzyaaac47)
xyzzyaaae47=sqrt(expval_gvec(1,xyzzyaaab47,xyzzyaaac47)*expval_gvec(1,&
&xyzzyaaab47,xyzzyaaac47)+expval_gvec(2,xyzzyaaab47,xyzzyaaac47)*expva&
&l_gvec(2,xyzzyaaab47,xyzzyaaac47)+expval_gvec(3,xyzzyaaab47,xyzzyaaac&
&47)*expval_gvec(3,xyzzyaaab47,xyzzyaaac47))
if(abs(xyzzyaaae47-xyzzyaaad47)>1.d-8)then
expval_ngvec_truncated(xyzzyaaac47)=xyzzyaaab47-1
exit
endif
enddo
enddo
endif
xyzzyaacz1(:)=0
do xyzzyaaab47=1,xyzzyaacw1
xyzzyaaap1(1:3)=xyzzyaadc1(1:3,1,xyzzyaaab47)
xyzzyaaaq1(1:3)=xyzzyaadc1(1:3,2,xyzzyaaab47)
xyzzyaaar1(1:3)=xyzzyaadc1(1:3,3,xyzzyaaab47)
do xyzzyaaaa47=1,expval_ngvec(xyzzyaaab47)
xyzzyaaaa1=nint(one_over_twopi*(expval_gvec(1,xyzzyaaaa47,xyzzyaaab47)&
&*xyzzyaaap1(1)+expval_gvec(2,xyzzyaaaa47,xyzzyaaab47)*xyzzyaaap1(2)+e&
&xpval_gvec(3,xyzzyaaaa47,xyzzyaaab47)*xyzzyaaap1(3)))
if(abs(xyzzyaaaa1)>xyzzyaacz1(xyzzyaaab47))xyzzyaacz1(xyzzyaaab47)=abs&
&(xyzzyaaaa1)
xyzzyaacy1(1,xyzzyaaaa47,xyzzyaaab47)=xyzzyaaaa1
xyzzyaaab1=nint(one_over_twopi*(expval_gvec(1,xyzzyaaaa47,xyzzyaaab47)&
&*xyzzyaaaq1(1)+expval_gvec(2,xyzzyaaaa47,xyzzyaaab47)*xyzzyaaaq1(2)+e&
&xpval_gvec(3,xyzzyaaaa47,xyzzyaaab47)*xyzzyaaaq1(3)))
if(abs(xyzzyaaab1)>xyzzyaacz1(xyzzyaaab47))xyzzyaacz1(xyzzyaaab47)=abs&
&(xyzzyaaab1)
xyzzyaacy1(2,xyzzyaaaa47,xyzzyaaab47)=xyzzyaaab1
xyzzyaaac1=nint(one_over_twopi*(expval_gvec(1,xyzzyaaaa47,xyzzyaaab47)&
&*xyzzyaaar1(1)+expval_gvec(2,xyzzyaaaa47,xyzzyaaab47)*xyzzyaaar1(2)+e&
&xpval_gvec(3,xyzzyaaaa47,xyzzyaaab47)*xyzzyaaar1(3)))
if(abs(xyzzyaaac1)>xyzzyaacz1(xyzzyaaab47))xyzzyaacz1(xyzzyaaab47)=abs&
&(xyzzyaaac1)
xyzzyaacy1(3,xyzzyaaaa47,xyzzyaaab47)=xyzzyaaac1
enddo
enddo
if(any(xyzzyaacz1(:)<1))then
call errstop('SETUP_FOURIER_BASIS','G vector sets in expval.data must &
&include non-zero G. Cutoff not large enough?')
endif
xyzzyaaaa47=maxval(xyzzyaacz1(:))
allocate(xyzzyaade1(3,-xyzzyaaaa47:xyzzyaaaa47,xyzzyaacw1),stat=xyzzya&
&aad1)
call check_alloc(xyzzyaaad1,'SETUP_FOURIER_BASIS','expval_mwork')
end subroutine xyzzyaaqw1
subroutine xyzzyaaqx1
implicit none
integer xyzzyaaaa48
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabc1)then
xyzzyaadg1=xyzzyaacw1
den_nsets=no_families
allocate(xyzzyaadj1(den_nsets),xyzzyaadi1(den_nsets),xyzzyaadh1(den_ns&
&ets),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_DENSITY','expval_den_weight')
do xyzzyaaaa48=1,den_nsets
xyzzyaadh1(xyzzyaaaa48)=xyzzyaaaa48
enddo
allocate(xyzzyaadl1(expval_ngvec(xyzzyaadg1),den_nsets),xyzzyaadk1(exp&
&val_ngvec(xyzzyaadg1),den_nsets),xyzzyaadm1(expval_ngvec(xyzzyaadg1),&
&den_nsets),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_DENSITY','expval_den')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaadn1(expval_ngvec(xyzzyaadg1),den_nsets),stat=xyzzyaaad&
&1)
call check_alloc(xyzzyaaad1,'SETUP_DENSITY','expval_den_dmc')
endif
xyzzyaadj1(:)=0.d0
xyzzyaadi1(:)=0.d0
xyzzyaadl1(:,:)=czero
xyzzyaadk1(:,:)=czero
if(inversion_symmetry)then
expval_complex_den=.false.
else
expval_complex_den=.true.
endif
if(isvmc)then
accumulation_method_den='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_den='DMC'
else
call errstop('SETUP_DENSITY','Confusion over what accumulation method &
&I''m using. Should not happen.')
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaadg1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast in setup_density')
call mpi_bcast(den_nsets,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast in setup_density')
if(am_slave)then
allocate(xyzzyaadk1(expval_ngvec(xyzzyaadg1),den_nsets),xyzzyaadm1(exp&
&val_ngvec(xyzzyaadg1),den_nsets),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_DENSITY','expval_den')
if(xyzzyaabb1)then
allocate(expval_den_mpc(expval_ngvec(xyzzyaadg1),den_nsets),stat=xyzzy&
&aaad1)
call check_alloc(xyzzyaaad1,'SETUP_DENSITY','expval_den_mpc')
endif
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaadn1(expval_ngvec(xyzzyaadg1),den_nsets),stat=xyzzyaaad&
&1)
call check_alloc(xyzzyaaad1,'SETUP_DENSITY','expval_den_dmc')
endif
xyzzyaadk1(:,:)=czero
if(xyzzyaabb1)expval_den_mpc=czero
endif
call mpi_bcast(expval_complex_den,1,mpi_logical,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'broadcast in setup_density')
endif
end subroutine xyzzyaaqx1
subroutine xyzzyaaqy1
implicit none
integer xyzzyaaaa49
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabd1)then
xyzzyaado1=xyzzyaacw1
xyzzyaadp1=nspin
allocate(xyzzyaads1(xyzzyaadp1),xyzzyaadr1(xyzzyaadp1),xyzzyaadq1(xyzz&
&yaadp1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_SPIN_DENSITY','expval_sden_weight')
do xyzzyaaaa49=1,xyzzyaadp1
xyzzyaadq1(xyzzyaaaa49)=xyzzyaaaa49
enddo
allocate(xyzzyaadu1(expval_ngvec(xyzzyaado1),xyzzyaadp1),xyzzyaadt1(ex&
&pval_ngvec(xyzzyaado1),xyzzyaadp1),xyzzyaadv1(expval_ngvec(xyzzyaado1&
&),xyzzyaadp1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_SPIN_DENSITY','expval_sden')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaadw1(expval_ngvec(xyzzyaado1),xyzzyaadp1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_SPIN_DENSITY','expval_sden_dmc')
endif
xyzzyaads1(:)=0.d0
xyzzyaadr1(:)=0.d0
xyzzyaadu1(:,:)=czero
xyzzyaadt1(:,:)=czero
if(inversion_symmetry)then
xyzzyaadx1=.false.
else
xyzzyaadx1=.true.
endif
if(isvmc)then
accumulation_method_sden='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_sden='DMC'
else
call errstop('SETUP_SPIN_DENSITY','Confusion over what accumulation me&
&thod I''m using. Should not happen.')
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaado1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sden_gset in setup_spin_density')
call mpi_bcast(xyzzyaadp1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sden_nsets in setup_spin_density')
if(am_slave)then
allocate(xyzzyaadt1(expval_ngvec(xyzzyaado1),xyzzyaadp1),xyzzyaadv1(ex&
&pval_ngvec(xyzzyaado1),xyzzyaadp1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_SPIN_DENSITY','expval_sden')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaadw1(expval_ngvec(xyzzyaado1),xyzzyaadp1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_SPIN_DENSITY','expval_sden_dmc')
endif
xyzzyaadt1(:,:)=czero
endif
call mpi_bcast(xyzzyaadx1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast in setup_spin_density')
endif
end subroutine xyzzyaaqy1
subroutine xyzzyaaqz1
implicit none
integer xyzzyaaaa50
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabe1)then
xyzzyaady1=xyzzyaacw1
xyzzyaadz1=nspin
allocate(xyzzyaaec1(xyzzyaadz1),xyzzyaaeb1(xyzzyaadz1),xyzzyaaea1(xyzz&
&yaadz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_SPIN_DENSITY_MATRIX','expval_sdenm_&
&weight')
do xyzzyaaaa50=1,xyzzyaadz1
xyzzyaaea1(xyzzyaaaa50)=xyzzyaaaa50
enddo
allocate(xyzzyaaee1(expval_ngvec(xyzzyaady1),4,xyzzyaadz1),xyzzyaaed1(&
&expval_ngvec(xyzzyaady1),4,xyzzyaadz1),xyzzyaaef1(expval_ngvec(xyzzya&
&ady1),4,xyzzyaadz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_SPIN_DENSITY_MATRIX','expval_sdenm'&
&)
xyzzyaaec1(:)=0.d0
xyzzyaaeb1(:)=0.d0
xyzzyaaee1(:,:,:)=czero
xyzzyaaed1(:,:,:)=czero
xyzzyaaeg1=.true.
if(isvmc)then
accumulation_method_sdenm='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_sdenm='DMC'
else
call errstop('SETUP_SPIN_DENSITY_MATRIX','Confusion over what accumula&
&tion method I''m using. Should not happen.')
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaady1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sdenm_gset in setup_spin_density')
call mpi_bcast(xyzzyaadz1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sdenm_nsets in setup_spin_density')
if(am_slave)then
allocate(xyzzyaaee1(expval_ngvec(xyzzyaady1),4,xyzzyaadz1),xyzzyaaed1(&
&expval_ngvec(xyzzyaady1),4,xyzzyaadz1),xyzzyaaef1(expval_ngvec(xyzzya&
&ady1),4,xyzzyaadz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_SPIN_DENSITY_MATRIX','expval_sdenm'&
&)
xyzzyaaed1(:,:,:)=czero
endif
call mpi_bcast(xyzzyaaeg1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast in setup_spin_density_matrix')
endif
end subroutine xyzzyaaqz1
subroutine xyzzyaara1(from_backup)
implicit none
logical,intent(in) :: from_backup
integer xyzzyaaaa51,xyzzyaaab51,xyzzyaaac51
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabl1)then
if(.not.from_backup)xyzzyaaeh1=xyzzyaacv1+1
xyzzyaaei1=nspin
allocate(xyzzyaaem1(xyzzyaaei1),xyzzyaael1(xyzzyaaei1),xyzzyaaek1(xyzz&
&yaaei1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR','expval_pcf_weight')
do xyzzyaaaa51=1,xyzzyaaei1
xyzzyaaek1(xyzzyaaaa51)=xyzzyaaaa51
enddo
allocate(xyzzyaaen1(expval_ngvec(xyzzyaaeh1),xyzzyaaei1),expval_pcf(ex&
&pval_ngvec(xyzzyaaeh1),xyzzyaaei1),xyzzyaaeo1(expval_ngvec(xyzzyaaeh1&
&),xyzzyaaei1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR','expval_pcf')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaaep1(expval_ngvec(xyzzyaaeh1),xyzzyaaei1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR','expval_pcf_dmc')
endif
xyzzyaaem1(:)=0.d0
xyzzyaael1(:)=0.d0
xyzzyaaen1(:,:)=czero
expval_pcf(:,:)=czero
if(inversion_symmetry)then
xyzzyaaeq1=.false.
else
xyzzyaaeq1=.true.
endif
if(isvmc)then
accumulation_method_pcf='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_pcf='DMC'
else
call errstop('SETUP_PAIR_CORR','Confusion over what accumulation metho&
&d I''m using. Should not happen.')
endif
xyzzyaaej1=input_type_fixed
pcf_rfix(1:3)=input_rfix(1:3)
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaaeh1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast pcf_gset in setup_pair_corr')
call mpi_bcast(xyzzyaaei1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast pcf_nsets in setup_pair_corr')
if(am_slave)then
allocate(xyzzyaaen1(expval_ngvec(xyzzyaaeh1),xyzzyaaei1),expval_pcf(ex&
&pval_ngvec(xyzzyaaeh1),xyzzyaaei1),xyzzyaaeo1(expval_ngvec(xyzzyaaeh1&
&),xyzzyaaei1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR','expval_pcf')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaaep1(expval_ngvec(xyzzyaaeh1),xyzzyaaei1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR','expval_pcf_dmc')
endif
expval_pcf(:,:)=czero
endif
call mpi_bcast(xyzzyaaeq1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast in setup_pair_corr')
call mpi_bcast(xyzzyaaej1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast pcf_type_fixed in setup_pair_corr')
call mpi_bcast(pcf_rfix,3,mpi_double_precision,0,mpi_comm_world,ierror&
&)
call checkmpi(ierror,'broadcast pcf_rfix in setup_pair_corr')
endif
xyzzyaaaa51=0
do xyzzyaaac51=1,nspin
do xyzzyaaab51=1,nele(xyzzyaaac51)
xyzzyaaaa51=xyzzyaaaa51+1
if(xyzzyaaac51==xyzzyaaej1.and.xyzzyaaab51==1)fixed_particle=xyzzyaaaa&
&51
enddo
enddo
particle_is_fixed=.false.
end subroutine xyzzyaara1
subroutine xyzzyaarb1
implicit none
integer xyzzyaaaa52,xyzzyaaab52,xyzzyaaac52
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabn1)then
if(isvmc)then
accumulation_method_pcfs='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_pcfs='DMC'
else
call errstop('SETUP_PAIR_CORR_SPH','Confusion over what accumulation m&
&ethod I''m using. Should not happen.')
endif
if(homogeneous_density)then
xyzzyaaet1=2
else
xyzzyaaet1=1
endif
if(pcfs_nbins_in>0)then
xyzzyaaes1=pcfs_nbins_in
else
xyzzyaaes1=5000
endif
if(pcfs_rcutoff_in>0.d0)then
if(isperiodic)then
xyzzyaaew1=min(wigner_seitz_radius-1.d-6,pcfs_rcutoff_in)
else
xyzzyaaet1=2
xyzzyaaew1=pcfs_rcutoff_in
endif
else
if(isperiodic)then
xyzzyaaew1=wigner_seitz_radius-1.d-6
else
xyzzyaaet1=2
xyzzyaaew1=max(mol_system_size(),10.d0)
endif
endif
if(xyzzyaaet1==1)then
xyzzyaaer1=nspin
xyzzyaaeu1=input_type_fixed
pcfs_rfix(1:3)=input_rfix(1:3)
else
xyzzyaaab52=0
do xyzzyaaaa52=1,nspin
xyzzyaaab52=xyzzyaaab52+xyzzyaaaa52
enddo
xyzzyaaer1=xyzzyaaab52
endif
allocate(xyzzyaafb1(xyzzyaaer1),xyzzyaaev1(2,xyzzyaaer1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR_SPH','expval_pcfs_weight_&
&start,pcfs_ptype')
if(xyzzyaaet1==1)then
do xyzzyaaaa52=1,xyzzyaaer1
xyzzyaaev1(1,xyzzyaaaa52)=input_type_fixed
xyzzyaaev1(2,xyzzyaaaa52)=xyzzyaaaa52
enddo
allocate(xyzzyaafd1(xyzzyaaes1,1,nspin),xyzzyaafc1(xyzzyaaes1,1,nspin)&
&,xyzzyaafe1(xyzzyaaes1,1,nspin),stat=xyzzyaaad1)
else
xyzzyaaac52=0
do xyzzyaaaa52=1,nspin
do xyzzyaaab52=xyzzyaaaa52,nspin
xyzzyaaac52=xyzzyaaac52+1
xyzzyaaev1(1,xyzzyaaac52)=xyzzyaaaa52
xyzzyaaev1(2,xyzzyaaac52)=xyzzyaaab52
enddo
enddo
allocate(xyzzyaafd1(xyzzyaaes1,nspin,nspin),xyzzyaafc1(xyzzyaaes1,nspi&
&n,nspin),xyzzyaafe1(xyzzyaaes1,nspin,nspin),stat=xyzzyaaad1)
endif
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR_SPH','expval_pcfs')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
if(xyzzyaaet1==1)then
allocate(xyzzyaaff1(xyzzyaaes1,1,nspin),stat=xyzzyaaad1)
else
allocate(xyzzyaaff1(xyzzyaaes1,nspin,nspin),stat=xyzzyaaad1)
endif
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR_SPH','expval_pcfs_dmc')
endif
xyzzyaafb1(:)=0.d0
xyzzyaaey1=0.d0
xyzzyaaez1=0.d0
xyzzyaafd1(:,:,:)=0.d0
xyzzyaafc1(:,:,:)=0.d0
xyzzyaafe1(:,:,:)=0.d0
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaaes1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast pcfs_nbins in setup_pair_corr_sph')
call mpi_bcast(xyzzyaaer1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast pcfs_nsets in setup_pair_corr_sph')
call mpi_bcast(xyzzyaaet1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast pcfs_accmode in setup_pair_corr_sph')
call mpi_bcast(xyzzyaaew1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast pcfs_rcutoff in setup_pair_corr_sph')
call mpi_bcast(xyzzyaaeu1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast pcfs_type_fixed in setup_pair_corr_sph&
&')
call mpi_bcast(pcfs_rfix,3,mpi_double_precision,0,mpi_comm_world,ierro&
&r)
call checkmpi(ierror,'broadcast pcfs_rfix in setup_pair_corr_sph')
if(am_slave)then
allocate(xyzzyaafb1(xyzzyaaer1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR_SPH','expval_pcfs_weight_&
&start')
if(xyzzyaaet1==1)then
allocate(xyzzyaafd1(xyzzyaaes1,1,nspin),xyzzyaafc1(xyzzyaaes1,1,nspin)&
&,xyzzyaafe1(xyzzyaaes1,1,nspin),stat=xyzzyaaad1)
else
allocate(xyzzyaafd1(xyzzyaaes1,nspin,nspin),xyzzyaafc1(xyzzyaaes1,nspi&
&n,nspin),xyzzyaafe1(xyzzyaaes1,nspin,nspin),stat=xyzzyaaad1)
endif
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR_SPH','expval_pcfs')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
if(xyzzyaaet1==1)then
allocate(xyzzyaaff1(xyzzyaaes1,1,nspin),stat=xyzzyaaad1)
else
allocate(xyzzyaaff1(xyzzyaaes1,nspin,nspin),stat=xyzzyaaad1)
endif
call check_alloc(xyzzyaaad1,'SETUP_PAIR_CORR_SPH','expval_pcfs_dmc')
endif
xyzzyaafb1(:)=0.d0
xyzzyaaey1=0.d0
xyzzyaaez1=0.d0
xyzzyaafd1(:,:,:)=0.d0
xyzzyaafc1(:,:,:)=0.d0
xyzzyaafe1(:,:,:)=0.d0
endif
endif
if(xyzzyaaet1==1)then
xyzzyaaaa52=0
do xyzzyaaac52=1,nspin
do xyzzyaaab52=1,nele(xyzzyaaac52)
xyzzyaaaa52=xyzzyaaaa52+1
if(xyzzyaaac52==xyzzyaaeu1.and.xyzzyaaab52==1)fixed_particle=xyzzyaaaa&
&52
enddo
enddo
endif
if(isperiodic)then
xyzzyaaex1=real(xyzzyaaes1,dp)/wigner_seitz_radius
else
xyzzyaaex1=real(xyzzyaaes1,dp)/xyzzyaaew1
endif
end subroutine xyzzyaarb1
subroutine xyzzyaarc1
implicit none
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabf1)then
xyzzyaafg1=nspin
allocate(xyzzyaafm1(periodicity,periodicity,xyzzyaafg1),xyzzyaafi1(xyz&
&zyaafg1),xyzzyaafl1(periodicity,periodicity,xyzzyaafg1),xyzzyaafh1(xy&
&zzyaafg1),xyzzyaafn1(periodicity,periodicity,xyzzyaafg1),xyzzyaafj1(x&
&yzzyaafg1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_LOC_TENSOR','expval_lt_*')
xyzzyaafi1(:)=0.d0
xyzzyaafh1(:)=0.d0
xyzzyaafj1(:)=0.d0
xyzzyaafm1(:,:,:)=czero
xyzzyaafl1(:,:,:)=czero
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaafo1(periodicity,periodicity,xyzzyaafg1),xyzzyaafk1(xyz&
&zyaafg1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'READ_LOC_TENSOR','expval_lt_dmc')
xyzzyaafo1=czero
xyzzyaafk1=0.d0
endif
xyzzyaafp1=.not.inversion_symmetry
if(isvmc)then
accumulation_method_lt='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_lt='DMC'
else
call errstop('SETUP_LOC_TENSOR','Confusion over what accumulation meth&
&odI''m using. Should not happen.')
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaafg1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast in setup_loc_tensor')
if(am_slave)then
allocate(xyzzyaafl1(periodicity,periodicity,xyzzyaafg1),xyzzyaafh1(xyz&
&zyaafg1),xyzzyaafn1(periodicity,periodicity,xyzzyaafg1),xyzzyaafj1(xy&
&zzyaafg1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_LOC_TENSOR','expval_lt')
xyzzyaafl1(:,:,:)=czero
xyzzyaafh1(:)=0.d0
xyzzyaafn1(:,:,:)=czero
xyzzyaafj1(:)=0.d0
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaafo1(periodicity,periodicity,xyzzyaafg1),xyzzyaafk1(xyz&
&zyaafg1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_LOC_TENSOR','expval_lt_dmc array')
xyzzyaafo1(:,:,:)=czero
xyzzyaafk1=0.d0
endif
endif
call mpi_bcast(xyzzyaafp1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast in setup_loc_tensor')
endif
end subroutine xyzzyaarc1
subroutine xyzzyaard1
implicit none
integer xyzzyaaaa54,xyzzyaaab54,xyzzyaaac54
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabh1)then
if(isvmc)then
accumulation_method_sf='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_sf='DMC'
else
call errstop('SETUP_STRUCTURE_FACTOR','Confusion over what accumulatio&
&n method I''m using. Should not happen.')
endif
xyzzyaafq1=0
do xyzzyaaaa54=1,nspin
xyzzyaafq1=xyzzyaafq1+xyzzyaaaa54
enddo
allocate(xyzzyaafy1(xyzzyaafq1),xyzzyaaft1(2,xyzzyaafq1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sf_weight&
&_start,sf_ptype')
if(homogeneous_density)then
xyzzyaafs1=0
else
xyzzyaafs1=nspin
allocate(xyzzyaagd1(xyzzyaafs1),xyzzyaage1(xyzzyaafs1),xyzzyaafu1(xyzz&
&yaafs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sfsden_we&
&ight/sfsden_ptype')
endif
xyzzyaafr1=xyzzyaacv1+1
xyzzyaaac54=0
do xyzzyaaaa54=1,nspin
do xyzzyaaab54=xyzzyaaaa54,nspin
xyzzyaaac54=xyzzyaaac54+1
xyzzyaaft1(1,xyzzyaaac54)=xyzzyaaaa54
xyzzyaaft1(2,xyzzyaaac54)=xyzzyaaab54
enddo
enddo
do xyzzyaaaa54=1,xyzzyaafs1
xyzzyaafu1(xyzzyaaaa54)=xyzzyaaaa54
enddo
allocate(xyzzyaafz1(expval_ngvec(xyzzyaafr1),xyzzyaafq1),xyzzyaaga1(ex&
&pval_ngvec(xyzzyaafr1),xyzzyaafq1),xyzzyaagf1(expval_ngvec(xyzzyaafr1&
&),nspin),xyzzyaagb1(expval_ngvec(xyzzyaafr1),xyzzyaafq1),stat=xyzzyaa&
&ad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sf')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagc1(expval_ngvec(xyzzyaafr1),xyzzyaafq1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sf_dmc')
endif
if(.not.homogeneous_density)then
allocate(xyzzyaagh1(expval_ngvec(xyzzyaafr1),xyzzyaafs1),xyzzyaagg1(ex&
&pval_ngvec(xyzzyaafr1),xyzzyaafs1),xyzzyaagi1(expval_ngvec(xyzzyaafr1&
&),xyzzyaafs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sf,expval&
&_sfsden')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagj1(expval_ngvec(xyzzyaafr1),xyzzyaafs1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sf_dmc,ex&
&pval_sfsden_dmc')
endif
endif
xyzzyaafy1(:)=0.d0
xyzzyaafv1=0.d0
xyzzyaafw1=0.d0
xyzzyaafz1(:,:)=0.d0
xyzzyaaga1(:,:)=0.d0
xyzzyaagf1(:,:)=czero
if(.not.homogeneous_density)then
xyzzyaagd1(:)=0.d0
xyzzyaage1(:)=0.d0
xyzzyaagh1(:,:)=czero
xyzzyaagg1(:,:)=czero
if(inversion_symmetry)then
xyzzyaagk1=.false.
else
xyzzyaagk1=.true.
endif
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaafr1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sf_gset in setup_structure_factor')
call mpi_bcast(xyzzyaafq1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sf_nsets in setup_structure_factor')
if(.not.homogeneous_density)then
call mpi_bcast(xyzzyaafs1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sfsden_nsets in setup_structure_factor&
&')
endif
if(am_slave)then
allocate(xyzzyaaga1(expval_ngvec(xyzzyaafr1),xyzzyaafq1),xyzzyaagf1(ex&
&pval_ngvec(xyzzyaafr1),nspin),xyzzyaagb1(expval_ngvec(xyzzyaafr1),xyz&
&zyaafq1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sf')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagc1(expval_ngvec(xyzzyaafr1),xyzzyaafq1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sf_dmc')
endif
xyzzyaagf1(:,:)=czero
xyzzyaaga1(:,:)=0.d0
if(.not.homogeneous_density)then
allocate(xyzzyaagg1(expval_ngvec(xyzzyaafr1),xyzzyaafs1),xyzzyaagi1(ex&
&pval_ngvec(xyzzyaafr1),xyzzyaafs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sfsden')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagj1(expval_ngvec(xyzzyaafr1),xyzzyaafs1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR','expval_sfsden_dm&
&c')
endif
xyzzyaagg1(:,:)=czero
endif
endif
call mpi_bcast(xyzzyaagk1,1,mpi_logical,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast in setup_structure_factor')
endif
end subroutine xyzzyaard1
subroutine xyzzyaare1
implicit none
integer xyzzyaaaa55,xyzzyaaab55,xyzzyaaac55
real(dp) xyzzyaaad55
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabj1)then
xyzzyaago1=expval_ka(1,2)
xyzzyaagp1=expval_kb(1,2)
xyzzyaagl1=expval_nk(1,2)
if(xyzzyaago1<0.d0.or.xyzzyaagp1<0.d0)call errstop('SETUP_STRUCTURE_FA&
&CTOR_SPH','Spherical structure factor: the coordinates of the k point&
&s in the radial grid must not be negative. Change the EXPVAL_KGRID bl&
&ock in input.')
if(xyzzyaago1>xyzzyaagp1.or.xyzzyaagl1<1.or.((xyzzyaago1==xyzzyaagp1).&
&and.xyzzyaagl1>1))call errstop('SETUP_STRUCTURE_FACTOR_SPH','Spherica&
&l structure factor: error in specification of radial k point grid. Ch&
&ange the EXPVAL_KGRID block in input.')
xyzzyaagm1=no_spairs(levels_spairs)
allocate(xyzzyaagq1(xyzzyaagm1),xyzzyaagn1(2,xyzzyaagm1),xyzzyaags1(xy&
&zzyaagm1),xyzzyaagt1(xyzzyaagl1),xyzzyaagy1(xyzzyaagm1),stat=xyzzyaaa&
&d1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR_SPH','expval_sf_sp&
&h_weight,sf_sph_ptype,expval_sf_sph_k')
xyzzyaaac55=0
do xyzzyaaaa55=1,nspin
do xyzzyaaab55=xyzzyaaaa55,nspin
xyzzyaaac55=xyzzyaaac55+1
xyzzyaagn1(1,xyzzyaaac55)=xyzzyaaaa55
xyzzyaagn1(2,xyzzyaaac55)=xyzzyaaab55
enddo
enddo
xyzzyaagt1(1)=xyzzyaago1
xyzzyaagt1(xyzzyaagl1)=xyzzyaagp1
if(xyzzyaagl1>2)then
xyzzyaaad55=(xyzzyaagp1-xyzzyaago1)/real(xyzzyaagl1-1,dp)
else
xyzzyaaad55=0.d0
endif
do xyzzyaaaa55=1,xyzzyaagl1-2
xyzzyaagt1(xyzzyaaaa55+1)=xyzzyaago1+real(xyzzyaaaa55,dp)*xyzzyaaad55
enddo
allocate(xyzzyaagu1(xyzzyaagl1,xyzzyaagm1),xyzzyaagv1(xyzzyaagl1,xyzzy&
&aagm1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR_SPH','expval_sf_sp&
&h')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagw1(xyzzyaagl1,xyzzyaagm1),xyzzyaagr1(xyzzyaagm1),stat&
&=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR_SPH','expval_sf_sp&
&h_dmc')
endif
xyzzyaagu1(:,:)=0.d0
xyzzyaagq1(:)=0.d0
xyzzyaagy1(:)=0.d0
if(isvmc)then
accumulation_method_sf_sph='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_sf_sph='DMC'
else
call errstop('SETUP_STRUCTURE_FACTOR_SPH','Confusion over what accumul&
&ation method is being used. Should not happen.')
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaagm1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sf_sph_nsets in setup_structure_factor&
&_sph')
call mpi_bcast(xyzzyaagl1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast sf_sph_nk in setup_structure_factor_sp&
&h')
call mpi_bcast(xyzzyaago1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast sf_sph_ka in setup_structure_factor_sp&
&h')
call mpi_bcast(xyzzyaagp1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast sf_sph_kb in setup_structure_factor_sp&
&h')
if(am_slave)then
allocate(xyzzyaagt1(xyzzyaagl1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR_SPH','expval_sf_sp&
&h_k')
xyzzyaagt1(1)=xyzzyaago1
xyzzyaagt1(xyzzyaagl1)=xyzzyaagp1
if(xyzzyaagl1>2)then
xyzzyaaad55=(xyzzyaagp1-xyzzyaago1)/real(xyzzyaagl1-1,dp)
else
xyzzyaaad55=0.d0
endif
do xyzzyaaaa55=1,xyzzyaagl1-2
xyzzyaagt1(xyzzyaaaa55+1)=xyzzyaago1+real(xyzzyaaaa55,dp)*xyzzyaaad55
enddo
allocate(xyzzyaagu1(xyzzyaagl1,xyzzyaagm1),xyzzyaagv1(xyzzyaagl1,xyzzy&
&aagm1),xyzzyaagq1(xyzzyaagm1),xyzzyaags1(xyzzyaagm1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR_SPH','expval_sf_sp&
&h')
if(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
allocate(xyzzyaagw1(xyzzyaagl1,xyzzyaagm1),xyzzyaagr1(xyzzyaagm1),stat&
&=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_STRUCTURE_FACTOR_SPH','expval_sf_sp&
&h_dmc')
endif
xyzzyaagu1(:,:)=0.d0
xyzzyaagq1(:)=0.d0
endif
endif
end subroutine xyzzyaare1
subroutine xyzzyaarf1
implicit none
integer xyzzyaaaa56,xyzzyaaab56,xyzzyaaac56,xyzzyaaad56
real(dp) xyzzyaaae56
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabp1)then
xyzzyaaha1=100
if(.not.xyzzyaabv1)then
xyzzyaahb1=20
if(on_top_ii==0)then
xyzzyaagz1=no_families
else
xyzzyaagz1=1
endif
else
xyzzyaahb1=xyzzyaaht1
xyzzyaagz1=xyzzyaahs1
endif
allocate(xyzzyaahg1(xyzzyaaha1,xyzzyaagz1),xyzzyaahk1(xyzzyaaha1,xyzzy&
&aagz1),xyzzyaahl1(xyzzyaaha1,xyzzyaagz1),xyzzyaahf1(xyzzyaaha1,xyzzya&
&agz1),xyzzyaahi1(xyzzyaaha1,xyzzyaagz1),xyzzyaahj1(xyzzyaaha1,xyzzyaa&
&gz1),xyzzyaahh1(xyzzyaaha1,xyzzyaagz1),xyzzyaahm1(xyzzyaaha1,xyzzyaag&
&z1),xyzzyaahn1(xyzzyaaha1,xyzzyaagz1),xyzzyaahc1(xyzzyaagz1),xyzzyaah&
&d1(nspin,xyzzyaagz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_ONEP_DENSITY_MATRIX','expval_onep_d&
&m')
xyzzyaahg1=0.d0
xyzzyaahk1=0.d0
xyzzyaahl1=0.d0
xyzzyaahf1=0.d0
xyzzyaahi1=0.d0
xyzzyaahj1=0.d0
xyzzyaahh1=0.d0
xyzzyaahm1=0.d0
xyzzyaahn1=0.d0
xyzzyaahc1=0
xyzzyaahd1=0
if(.not.xyzzyaabv1)then
if(on_top_ii==0)then
do xyzzyaaaa56=1,no_families
xyzzyaaac56=0
do xyzzyaaab56=1,nspin
if(which_fam(xyzzyaaab56)==xyzzyaaaa56)then
xyzzyaaac56=xyzzyaaac56+1
xyzzyaahd1(xyzzyaaac56,xyzzyaaaa56)=xyzzyaaab56
endif
enddo
xyzzyaahc1(xyzzyaaaa56)=xyzzyaaac56
enddo
else
xyzzyaahd1(1,1)=which_spin(on_top_ii)
xyzzyaahc1(1)=1
endif
if(isvmc)then
accumulation_method_onep_dm='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_onep_dm='DMC'
else
call errstop('SETUP_ONEP_DENSITY_MATRIX','Confusion over what accumula&
&tion method I''m using. Should not happen.')
endif
else
xyzzyaahd1=xyzzyaahv1
xyzzyaahc1=xyzzyaahu1
accumulation_method_onep_dm=accumulation_method_mom_den
endif
if(accumulation_method_onep_dm=='DMC')then
allocate(xyzzyaahr1(xyzzyaaha1,xyzzyaagz1),xyzzyaahp1(xyzzyaaha1,xyzzy&
&aagz1),              xyzzyaahq1(xyzzyaaha1,xyzzyaagz1),stat=xyzzyaaad&
&1)
call check_alloc(xyzzyaaad1,'SETUP_ONEP_DENSITY_MATRIX','expval_onep_d&
&m_*_dmc')
xyzzyaahr1=0.d0
xyzzyaahp1=0.d0
xyzzyaahq1=0.d0
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaagz1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast onep_dm_nsets in setup_onep_density_ma&
&trix')
call mpi_bcast(xyzzyaaha1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast onep_dm_nbins in setup_onep_density_ma&
&trix')
call mpi_bcast(xyzzyaahb1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast onep_dm_nrandom in setup_onep_density_&
&matrix')
call mpi_bcast(accumulation_method_onep_dm,3,mpi_character,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'broadcast accumulation_method_onep_dm in setup_o&
&nep_density_matrix')
if(am_slave)then
allocate(xyzzyaahg1(xyzzyaaha1,xyzzyaagz1),xyzzyaahk1(xyzzyaaha1,xyzzy&
&aagz1),xyzzyaahl1(xyzzyaaha1,xyzzyaagz1),xyzzyaahf1(xyzzyaaha1,xyzzya&
&agz1),xyzzyaahi1(xyzzyaaha1,xyzzyaagz1),xyzzyaahj1(xyzzyaaha1,xyzzyaa&
&gz1),xyzzyaahh1(xyzzyaaha1,xyzzyaagz1),xyzzyaahm1(xyzzyaaha1,xyzzyaag&
&z1),xyzzyaahn1(xyzzyaaha1,xyzzyaagz1),xyzzyaahc1(xyzzyaagz1),xyzzyaah&
&d1(nspin,xyzzyaagz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_ONEP_DENSITY_MATRIX','expval_onep_d&
&m')
xyzzyaahg1=0.d0
xyzzyaahk1=0.d0
xyzzyaahl1=0.d0
xyzzyaahf1=0.d0
xyzzyaahi1=0.d0
xyzzyaahj1=0.d0
xyzzyaahh1=0.d0
xyzzyaahm1=0.d0
xyzzyaahn1=0.d0
xyzzyaahc1=0
xyzzyaahd1=0
if(accumulation_method_onep_dm=='DMC')then
allocate(xyzzyaahr1(xyzzyaaha1,xyzzyaagz1),xyzzyaahp1(xyzzyaaha1,xyzzy&
&aagz1),              xyzzyaahq1(xyzzyaaha1,xyzzyaagz1),stat=xyzzyaaad&
&1)
call check_alloc(xyzzyaaad1,'SETUP_ONEP_DENSITY_MATRIX','expval_onep_d&
&m_*_dmc')
xyzzyaahr1=0.d0
xyzzyaahp1=0.d0
xyzzyaahq1=0.d0
endif
endif
call mpi_bcast(xyzzyaahc1,xyzzyaagz1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast onep_dm_nptypes_in_set in setup_onep_d&
&ensity_matrix')
call mpi_bcast(xyzzyaahd1,xyzzyaagz1*nspin,mpi_integer,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'broadcast onep_dm_ptype_in_set in setup_onep_den&
&sity_matrix')
endif
xyzzyaahe1=real(xyzzyaaha1,dp)/wigner_seitz_radius
select case(periodicity)
case(1)
xyzzyaaae56=1.d0/xyzzyaacq1
case(2)
xyzzyaaae56=1.d0/xyzzyaacp1
case(3)
xyzzyaaae56=1.d0/xyzzyaacr1
end select
allocate(xyzzyaaho1(xyzzyaagz1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_ONEP_DENSITY_MATRIX','onep_dm_facto&
&r')
xyzzyaaho1=0.d0
if(on_top_ii==0)then
do xyzzyaaad56=1,xyzzyaagz1
do xyzzyaaaa56=1,xyzzyaahc1(xyzzyaaad56)
xyzzyaaho1(xyzzyaaad56)=xyzzyaaho1(xyzzyaaad56)+xyzzyaaae56*real(nele(&
&xyzzyaahd1(xyzzyaaaa56,xyzzyaaad56)),dp)
enddo
enddo
else
xyzzyaaho1(1)=xyzzyaaae56
endif
end subroutine xyzzyaarf1
subroutine xyzzyaarg1
implicit none
integer xyzzyaaaa57,xyzzyaaab57,xyzzyaaac57,xyzzyaaad57,xyzzyaaae57,xy&
&zzyaaaf57,xyzzyaaag57,xyzzyaaah57,xyzzyaaai57,xyzzyaaaj57,xyzzyaaak57&
&,xyzzyaaal57,xyzzyaaam57
logical xyzzyaaan57
real(dp) xyzzyaaao57
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabr1)then
xyzzyaail1=100
xyzzyaaik1=no_spairs(min(1,levels_spairs))
xyzzyaaim1=20
xyzzyaaiv1=min(1.d0,2.d0/dble(maxval(nele)))
allocate(xyzzyaaix1(xyzzyaail1,xyzzyaaik1),xyzzyaajb1(xyzzyaail1,xyzzy&
&aaik1),xyzzyaajc1(xyzzyaail1,xyzzyaaik1),xyzzyaaiw1(xyzzyaail1,xyzzya&
&aik1),xyzzyaaiz1(xyzzyaail1,xyzzyaaik1),xyzzyaaja1(xyzzyaail1,xyzzyaa&
&ik1),xyzzyaaiy1(xyzzyaail1,xyzzyaaik1),xyzzyaajd1(xyzzyaail1,xyzzyaai&
&k1),xyzzyaaje1(xyzzyaail1,xyzzyaaik1),xyzzyaain1(xyzzyaaik1),xyzzyaai&
&o1(2,max_spin_pairs,xyzzyaaik1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DENSITY_MATRIX','expval_twop_d&
&m')
xyzzyaaix1=0.d0
xyzzyaajb1=0.d0
xyzzyaajc1=0.d0
xyzzyaaiw1=0.d0
xyzzyaaiz1=0.d0
xyzzyaaja1=0.d0
xyzzyaaiy1=0.d0
xyzzyaajd1=0.d0
xyzzyaaje1=0.d0
xyzzyaain1=0
xyzzyaaio1=0
do xyzzyaaac57=1,no_spairs(min(1,levels_spairs))
xyzzyaaad57=0
do xyzzyaaaa57=1,nspin
do xyzzyaaab57=xyzzyaaaa57,nspin
if(which_spair(xyzzyaaaa57,xyzzyaaab57,min(1,levels_spairs))==xyzzyaaa&
&c57)then
xyzzyaaad57=xyzzyaaad57+1
xyzzyaaio1(1,xyzzyaaad57,xyzzyaaac57)=xyzzyaaaa57
xyzzyaaio1(2,xyzzyaaad57,xyzzyaaac57)=xyzzyaaab57
endif
enddo
enddo
xyzzyaain1(xyzzyaaac57)=xyzzyaaad57
enddo
if(isvmc)then
accumulation_method_twop_dm='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_twop_dm='DMC'
allocate(xyzzyaajh1(xyzzyaail1,xyzzyaaik1),xyzzyaajf1(xyzzyaail1,xyzzy&
&aaik1),xyzzyaajg1(xyzzyaail1,xyzzyaaik1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DENSITY_MATRIX','expval_twop_d&
&m_*_dmc')
xyzzyaajh1=0.d0
xyzzyaajf1=0.d0
xyzzyaajg1=0.d0
else
call errstop('SETUP_TWOP_DENSITY_MATRIX','Confusion over what accumula&
&tion method I''m using. Should not happen.')
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaaik1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast twop_dm_nsets in setup_twop_density_ma&
&trix')
call mpi_bcast(xyzzyaail1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast twop_dm_nbins in setup_twop_density_ma&
&trix')
call mpi_bcast(xyzzyaaim1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast twop_dm_nrandom in setup_twop_density_&
&matrix')
call mpi_bcast(xyzzyaaiv1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast twop_dm_frac_pairs in setup_twop_densi&
&ty_matrix')
call mpi_bcast(accumulation_method_twop_dm,3,mpi_character,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'broadcast accumulation_method_twop_dm in setup_t&
&wop_density_matrix')
if(am_slave)then
allocate(xyzzyaaix1(xyzzyaail1,xyzzyaaik1),  xyzzyaajb1(xyzzyaail1,xyz&
&zyaaik1),xyzzyaajc1(xyzzyaail1,xyzzyaaik1),xyzzyaaiw1(xyzzyaail1,xyzz&
&yaaik1),xyzzyaaiz1(xyzzyaail1,xyzzyaaik1),xyzzyaaja1(xyzzyaail1,xyzzy&
&aaik1),xyzzyaaiy1(xyzzyaail1,xyzzyaaik1),xyzzyaajd1(xyzzyaail1,xyzzya&
&aik1),xyzzyaaje1(xyzzyaail1,xyzzyaaik1),xyzzyaain1(xyzzyaaik1),xyzzya&
&aio1(2,max_spin_pairs,xyzzyaaik1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DENSITY_MATRIX','expval_twop_d&
&m')
xyzzyaaix1=0.d0
xyzzyaajb1=0.d0
xyzzyaajc1=0.d0
xyzzyaaiw1=0.d0
xyzzyaaiz1=0.d0
xyzzyaaja1=0.d0
xyzzyaaiy1=0.d0
xyzzyaajd1=0.d0
xyzzyaaje1=0.d0
xyzzyaain1=0
xyzzyaaio1=0
if(accumulation_method_twop_dm=='DMC')then
allocate(xyzzyaajh1(xyzzyaail1,xyzzyaaik1),xyzzyaajf1(xyzzyaail1,xyzzy&
&aaik1),              xyzzyaajg1(xyzzyaail1,xyzzyaaik1),stat=xyzzyaaad&
&1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DENSITY_MATRIX','expval_twop_d&
&m_*_dmc')
xyzzyaajh1=0.d0
xyzzyaajf1=0.d0
xyzzyaajg1=0.d0
endif
endif
call mpi_bcast(xyzzyaain1,xyzzyaaik1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast twop_dm_nptypes_in_set in setup_twop_d&
&ensity_matrix')
call mpi_bcast(xyzzyaaio1,xyzzyaaik1*2*max_spin_pairs,mpi_integer,0,mp&
&i_comm_world,ierror)
call checkmpi(ierror,'broadcast twop_dm_ptype_in_set in setup_twop_den&
&sity_matrix')
endif
xyzzyaaiu1=real(xyzzyaail1,dp)/wigner_seitz_radius
select case(periodicity)
case(1)
xyzzyaaao57=1.d0/xyzzyaacq1**2
case(2)
xyzzyaaao57=1.d0/xyzzyaacp1**2
case(3)
xyzzyaaao57=1.d0/xyzzyaacr1**2
end select
allocate(xyzzyaaip1(xyzzyaaik1),xyzzyaais1(xyzzyaaik1),xyzzyaaji1(xyzz&
&yaaik1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DENSITY_MATRIX','twop_dm_npair&
&_in_set')
if(twop_dm_mom)then
allocate(xyzzyaanv1(xyzzyaaik1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DENSITY_MATRIX_mom','twop_dm_m&
&om_factor')
endif
do xyzzyaaae57=1,xyzzyaaik1
xyzzyaaaf57=0
xyzzyaaam57=0
do xyzzyaaaa57=1,xyzzyaain1(xyzzyaaae57)
xyzzyaaai57=xyzzyaaio1(1,xyzzyaaaa57,xyzzyaaae57)
xyzzyaaaj57=xyzzyaaio1(2,xyzzyaaaa57,xyzzyaaae57)
if(xyzzyaaai57/=xyzzyaaaj57)then
xyzzyaaaf57=xyzzyaaaf57+nele(xyzzyaaai57)*nele(xyzzyaaaj57)
xyzzyaaam57=xyzzyaaam57+nele(xyzzyaaai57)*nele(xyzzyaaaj57)
else
xyzzyaaaf57=xyzzyaaaf57+(nele(xyzzyaaai57)*(nele(xyzzyaaaj57)-1))/2
xyzzyaaam57=xyzzyaaam57+nele(xyzzyaaai57)*(nele(xyzzyaaaj57)-1)
endif
enddo
xyzzyaaip1(xyzzyaaae57)=xyzzyaaaf57
xyzzyaais1(xyzzyaaae57)=nint(xyzzyaaaf57*xyzzyaaiv1)
xyzzyaaji1(xyzzyaaae57)=xyzzyaaao57*dble(xyzzyaaam57)
if(twop_dm_mom)xyzzyaanv1(xyzzyaaae57)=xyzzyaaao57*dble(xyzzyaaam57)
if(xyzzyaaaf57/=0.and.xyzzyaais1(xyzzyaaae57)<1)call errstop('SETUP_TW&
&OP_DENSITY_MATRIX','Fraction of particle-pairs too small.')
enddo
allocate(xyzzyaaiq1(2,maxval(xyzzyaaip1),xyzzyaaik1),xyzzyaair1(maxval&
&(xyzzyaais1),xyzzyaaik1),xyzzyaait1(maxval(xyzzyaais1)),stat=xyzzyaaa&
&d1)
xyzzyaair1=0
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DENSITY_MATRIX','twop_dm_pair'&
&)
do xyzzyaaae57=1,xyzzyaaik1
xyzzyaaac57=0
do xyzzyaaaa57=1,xyzzyaain1(xyzzyaaae57)
xyzzyaaai57=xyzzyaaio1(1,xyzzyaaaa57,xyzzyaaae57)
xyzzyaaaj57=xyzzyaaio1(2,xyzzyaaaa57,xyzzyaaae57)
xyzzyaaan57=(xyzzyaaai57==xyzzyaaaj57)
do xyzzyaaak57=1,nele(xyzzyaaai57)
do xyzzyaaal57=1,nele(xyzzyaaaj57)
if(xyzzyaaan57.and.xyzzyaaak57>=xyzzyaaal57)cycle
xyzzyaaac57=xyzzyaaac57+1
xyzzyaaag57=which_ii(xyzzyaaak57,xyzzyaaai57)
xyzzyaaah57=which_ii(xyzzyaaal57,xyzzyaaaj57)
xyzzyaaiq1(1,xyzzyaaac57,xyzzyaaae57)=xyzzyaaag57
xyzzyaaiq1(2,xyzzyaaac57,xyzzyaaae57)=xyzzyaaah57
enddo
enddo
enddo
enddo
if(xyzzyaaiv1==1.d0)then
do xyzzyaaae57=1,xyzzyaaik1
do xyzzyaaaa57=1,xyzzyaaip1(xyzzyaaae57)
xyzzyaair1(xyzzyaaaa57,xyzzyaaae57)=xyzzyaaaa57
enddo
enddo
endif
end subroutine xyzzyaarg1
subroutine xyzzyaarh1
implicit none
integer xyzzyaaaa58,xyzzyaaab58,xyzzyaaac58,xyzzyaaad58,xyzzyaaae58,xy&
&zzyaaaf58,xyzzyaaag58,xyzzyaaah58,xyzzyaaai58,xyzzyaaaj58,xyzzyaaak58&
&,xyzzyaaal58,xyzzyaaam58
logical xyzzyaaan58
real(dp) xyzzyaaao58
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabt1)then
xyzzyaajk1=100
xyzzyaajj1=no_spairs(min(1,levels_spairs))
xyzzyaajl1=20
xyzzyaaju1=min(1.d0,2.d0/dble(maxval(nele)))
allocate(xyzzyaajw1(xyzzyaajk1,xyzzyaajj1),xyzzyaaka1(xyzzyaajk1,xyzzy&
&aajj1),xyzzyaakb1(xyzzyaajk1,xyzzyaajj1),xyzzyaajv1(xyzzyaajk1,xyzzya&
&ajj1),xyzzyaajy1(xyzzyaajk1,xyzzyaajj1),xyzzyaajz1(xyzzyaajk1,xyzzyaa&
&jj1),xyzzyaajx1(xyzzyaajk1,xyzzyaajj1),xyzzyaakc1(xyzzyaajk1,xyzzyaaj&
&j1),xyzzyaakd1(xyzzyaajk1,xyzzyaajj1),xyzzyaajm1(xyzzyaajj1),xyzzyaaj&
&n1(2,max_spin_pairs,xyzzyaajj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRACTION','expval_cond_frac')
xyzzyaajw1=0.d0
xyzzyaaka1=0.d0
xyzzyaakb1=0.d0
xyzzyaajv1=0.d0
xyzzyaajy1=0.d0
xyzzyaajz1=0.d0
xyzzyaajx1=0.d0
xyzzyaakc1=0.d0
xyzzyaakd1=0.d0
xyzzyaajm1=0
xyzzyaajn1=0
do xyzzyaaac58=1,no_spairs(min(1,levels_spairs))
xyzzyaaad58=0
do xyzzyaaaa58=1,nspin
do xyzzyaaab58=xyzzyaaaa58,nspin
if(which_spair(xyzzyaaaa58,xyzzyaaab58,min(1,levels_spairs))==xyzzyaaa&
&c58)then
xyzzyaaad58=xyzzyaaad58+1
xyzzyaajn1(1,xyzzyaaad58,xyzzyaaac58)=xyzzyaaaa58
xyzzyaajn1(2,xyzzyaaad58,xyzzyaaac58)=xyzzyaaab58
endif
enddo
enddo
xyzzyaajm1(xyzzyaaac58)=xyzzyaaad58
enddo
if(isvmc)then
accum_method_cond_frac='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accum_method_cond_frac='DMC'
allocate(xyzzyaakg1(xyzzyaajk1,xyzzyaajj1),xyzzyaake1(xyzzyaajk1,xyzzy&
&aajj1),xyzzyaakf1(xyzzyaajk1,xyzzyaajj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRACTION','expval_cond_frac_*_&
&dmc')
xyzzyaakg1=0.d0
xyzzyaake1=0.d0
xyzzyaakf1=0.d0
else
call errstop('SETUP_COND_FRACTION','Confusion over what accumulation m&
&ethod I''m using. Should not happen.')
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaajj1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast cond_frac_nsets in setup_cond_fraction&
&')
call mpi_bcast(xyzzyaajk1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast cond_frac_nbins in setup_cond_fraction&
&')
call mpi_bcast(xyzzyaajl1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast cond_frac_nrandom in setup_cond_fracti&
&on')
call mpi_bcast(xyzzyaaju1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast cond_frac_frac_pairs in setup_cond_fra&
&ction')
call mpi_bcast(accum_method_cond_frac,3,mpi_character,0,mpi_comm_world&
&,ierror)
call checkmpi(ierror,'broadcast accum_method_cond_frac in setup_cond_f&
&raction')
if(am_slave)then
allocate(xyzzyaajw1(xyzzyaajk1,xyzzyaajj1),  xyzzyaaka1(xyzzyaajk1,xyz&
&zyaajj1),xyzzyaakb1(xyzzyaajk1,xyzzyaajj1),xyzzyaajv1(xyzzyaajk1,xyzz&
&yaajj1),xyzzyaajy1(xyzzyaajk1,xyzzyaajj1),xyzzyaajz1(xyzzyaajk1,xyzzy&
&aajj1),xyzzyaajx1(xyzzyaajk1,xyzzyaajj1),xyzzyaakc1(xyzzyaajk1,xyzzya&
&ajj1),xyzzyaakd1(xyzzyaajk1,xyzzyaajj1),xyzzyaajm1(xyzzyaajj1),xyzzya&
&ajn1(2,max_spin_pairs,xyzzyaajj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRACTION','expval_cond_frac')
xyzzyaajw1=0.d0
xyzzyaaka1=0.d0
xyzzyaakb1=0.d0
xyzzyaajv1=0.d0
xyzzyaajy1=0.d0
xyzzyaajz1=0.d0
xyzzyaajx1=0.d0
xyzzyaakc1=0.d0
xyzzyaakd1=0.d0
xyzzyaajm1=0
xyzzyaajn1=0
if(accum_method_cond_frac=='DMC')then
allocate(xyzzyaakg1(xyzzyaajk1,xyzzyaajj1),xyzzyaake1(xyzzyaajk1,xyzzy&
&aajj1),xyzzyaakf1(xyzzyaajk1,xyzzyaajj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRACTION','expval_cond_frac_*_&
&dmc')
xyzzyaakg1=0.d0
xyzzyaake1=0.d0
xyzzyaakf1=0.d0
endif
endif
call mpi_bcast(xyzzyaajm1,xyzzyaajj1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast cond_frac_nptypes_in_set in setup_cond&
&_fraction')
call mpi_bcast(xyzzyaajn1,xyzzyaajj1*2*max_spin_pairs,mpi_integer,0,mp&
&i_comm_world,ierror)
call checkmpi(ierror,'broadcast cond_frac_ptype_in_set in setup_cond_f&
&raction')
endif
xyzzyaajt1=real(xyzzyaajk1,dp)/wigner_seitz_radius
select case(periodicity)
case(1)
xyzzyaaao58=1.d0/xyzzyaacq1**2
case(2)
xyzzyaaao58=1.d0/xyzzyaacp1**2
case(3)
xyzzyaaao58=1.d0/xyzzyaacr1**2
end select
allocate(xyzzyaajo1(xyzzyaajj1),xyzzyaajr1(xyzzyaajj1),xyzzyaakh1(xyzz&
&yaajj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRACTION','cond_frac_npair_in_&
&set')
if(cond_frac_mom)then
allocate(xyzzyaaop1(xyzzyaajj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP__COND_FRACTION','cond_frac_mom_fact&
&or')
endif
do xyzzyaaae58=1,xyzzyaajj1
xyzzyaaaf58=0
xyzzyaaam58=0
do xyzzyaaaa58=1,xyzzyaajm1(xyzzyaaae58)
xyzzyaaai58=xyzzyaajn1(1,xyzzyaaaa58,xyzzyaaae58)
xyzzyaaaj58=xyzzyaajn1(2,xyzzyaaaa58,xyzzyaaae58)
if(xyzzyaaai58/=xyzzyaaaj58)then
xyzzyaaaf58=xyzzyaaaf58+nele(xyzzyaaai58)*nele(xyzzyaaaj58)
xyzzyaaam58=xyzzyaaam58+nele(xyzzyaaai58)*nele(xyzzyaaaj58)
else
xyzzyaaaf58=xyzzyaaaf58+(nele(xyzzyaaai58)*(nele(xyzzyaaaj58)-1))/2
xyzzyaaam58=xyzzyaaam58+nele(xyzzyaaai58)*(nele(xyzzyaaaj58)-1)
endif
enddo
xyzzyaajo1(xyzzyaaae58)=xyzzyaaaf58
xyzzyaajr1(xyzzyaaae58)=nint(xyzzyaaaf58*xyzzyaaju1)
xyzzyaakh1(xyzzyaaae58)=xyzzyaaao58*dble(xyzzyaaam58)
if(cond_frac_mom)xyzzyaaop1(xyzzyaaae58)=xyzzyaaao58*dble(xyzzyaaam58)
if(xyzzyaaaf58/=0.and.xyzzyaajr1(xyzzyaaae58)<1)call errstop('SETUP_CO&
&ND_FRACTION','Fraction of particle-pairs too small.')
enddo
allocate(xyzzyaajp1(2,maxval(xyzzyaajo1),xyzzyaajj1),xyzzyaajq1(maxval&
&(xyzzyaajr1),xyzzyaajj1),xyzzyaajs1(maxval(xyzzyaajr1)),stat=xyzzyaaa&
&d1)
xyzzyaajq1=0
call check_alloc(xyzzyaaad1,'SETUP_COND_FRACTION','cond_frac_pair')
do xyzzyaaae58=1,xyzzyaajj1
xyzzyaaac58=0
do xyzzyaaaa58=1,xyzzyaajm1(xyzzyaaae58)
xyzzyaaai58=xyzzyaajn1(1,xyzzyaaaa58,xyzzyaaae58)
xyzzyaaaj58=xyzzyaajn1(2,xyzzyaaaa58,xyzzyaaae58)
xyzzyaaan58=(xyzzyaaai58==xyzzyaaaj58)
do xyzzyaaak58=1,nele(xyzzyaaai58)
do xyzzyaaal58=1,nele(xyzzyaaaj58)
if(xyzzyaaan58.and.xyzzyaaak58>=xyzzyaaal58)cycle
xyzzyaaac58=xyzzyaaac58+1
xyzzyaaag58=which_ii(xyzzyaaak58,xyzzyaaai58)
xyzzyaaah58=which_ii(xyzzyaaal58,xyzzyaaaj58)
xyzzyaajp1(1,xyzzyaaac58,xyzzyaaae58)=xyzzyaaag58
xyzzyaajp1(2,xyzzyaaac58,xyzzyaaae58)=xyzzyaaah58
enddo
enddo
enddo
enddo
if(xyzzyaaju1==1.d0)then
do xyzzyaaae58=1,xyzzyaajj1
do xyzzyaaaa58=1,xyzzyaajo1(xyzzyaaae58)
xyzzyaajq1(xyzzyaaaa58,xyzzyaaae58)=xyzzyaaaa58
enddo
enddo
endif
end subroutine xyzzyaarh1
subroutine xyzzyaari1
implicit none
integer xyzzyaaaa59,xyzzyaaab59,xyzzyaaac59,xyzzyaaad59
real(dp) xyzzyaaae59
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabv1)then
if(.not.xyzzyaabp1)then
xyzzyaaht1=20
if(on_top_ii==0)then
xyzzyaahs1=no_families
else
xyzzyaahs1=1
endif
else
xyzzyaaht1=xyzzyaahb1
xyzzyaahs1=xyzzyaagz1
endif
allocate(xyzzyaahx1(xyzzyaahs1),xyzzyaaib1(expval_ngvec(1),xyzzyaahs1)&
&,xyzzyaaic1(expval_ngvec(1),xyzzyaahs1),xyzzyaahw1(xyzzyaahs1),xyzzya&
&ahz1(expval_ngvec(1),xyzzyaahs1),xyzzyaaia1(expval_ngvec(1),xyzzyaahs&
&1),xyzzyaahy1(xyzzyaahs1),xyzzyaaid1(expval_ngvec(1),xyzzyaahs1),xyzz&
&yaaie1(expval_ngvec(1),xyzzyaahs1),xyzzyaahu1(xyzzyaahs1),xyzzyaahv1(&
&nspin,xyzzyaahs1),xyzzyaaij1(3,expval_ngvec(1)),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOM_DEN','expval_mom_den')
xyzzyaahx1=0.d0
xyzzyaaib1=0.d0
xyzzyaaic1=0.d0
xyzzyaahw1=0.d0
xyzzyaahz1=0.d0
xyzzyaaia1=0.d0
xyzzyaahy1=0.d0
xyzzyaaid1=0.d0
xyzzyaaie1=0.d0
xyzzyaahu1=0
xyzzyaahv1=0
do xyzzyaaaa59=1,expval_ngvec(1)
xyzzyaaij1(:,xyzzyaaaa59)=expval_gvec(:,xyzzyaaaa59,1)-k_offset(:)
enddo
if(.not.xyzzyaabp1)then
if(on_top_ii==0)then
do xyzzyaaaa59=1,no_families
xyzzyaaac59=0
do xyzzyaaab59=1,nspin
if(which_fam(xyzzyaaab59)==xyzzyaaaa59)then
xyzzyaaac59=xyzzyaaac59+1
xyzzyaahv1(xyzzyaaac59,xyzzyaaaa59)=xyzzyaaab59
endif
enddo
xyzzyaahu1(xyzzyaaaa59)=xyzzyaaac59
enddo
else
xyzzyaahv1(1,1)=which_spin(on_top_ii)
xyzzyaahu1(1)=1
endif
if(isvmc)then
accumulation_method_mom_den='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_mom_den='DMC'
else
call errstop('SETUP_MOM_DEN','Confusion over what accumulation method &
&I''m using. Should not happen.')
endif
else
xyzzyaahv1=xyzzyaahd1
xyzzyaahu1=xyzzyaahc1
accumulation_method_mom_den=accumulation_method_onep_dm
endif
if(accumulation_method_mom_den=='DMC')then
allocate(xyzzyaaii1(xyzzyaahs1),xyzzyaaig1(expval_ngvec(1),xyzzyaahs1)&
&,xyzzyaaih1(expval_ngvec(1),xyzzyaahs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOM_DEN','expval_mom_den_*_dmc')
xyzzyaaii1=0.d0
xyzzyaaig1=0.d0
xyzzyaaih1=0.d0
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaahs1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast mom_den_nsets in setup_mom_den')
call mpi_bcast(expval_ngvec(1),1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast expval_ngvec(1) in setup_mom_den')
call mpi_bcast(xyzzyaaht1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast mom_den_nrandom in setup_mom_den')
call mpi_bcast(accumulation_method_mom_den,3,mpi_character,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'broadcast accumulation_method_mom_den in setup_m&
&om_den')
if(am_slave)then
allocate(xyzzyaahx1(xyzzyaahs1),xyzzyaaib1(expval_ngvec(1),xyzzyaahs1)&
&,xyzzyaaic1(expval_ngvec(1),xyzzyaahs1),xyzzyaahw1(xyzzyaahs1),xyzzya&
&ahz1(expval_ngvec(1),xyzzyaahs1),xyzzyaaia1(expval_ngvec(1),xyzzyaahs&
&1),xyzzyaahy1(xyzzyaahs1),xyzzyaaid1(expval_ngvec(1),xyzzyaahs1),xyzz&
&yaaie1(expval_ngvec(1),xyzzyaahs1),xyzzyaahu1(xyzzyaahs1),xyzzyaahv1(&
&nspin,xyzzyaahs1),xyzzyaaij1(3,expval_ngvec(1)),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOM_DEN','expval_mom_den')
xyzzyaahx1=0.d0
xyzzyaaib1=0.d0
xyzzyaaic1=0.d0
xyzzyaahw1=0.d0
xyzzyaahz1=0.d0
xyzzyaaia1=0.d0
xyzzyaahy1=0.d0
xyzzyaaid1=0.d0
xyzzyaaie1=0.d0
xyzzyaahu1=0
xyzzyaahv1=0
do xyzzyaaaa59=1,expval_ngvec(1)
xyzzyaaij1(:,xyzzyaaaa59)=expval_gvec(:,xyzzyaaaa59,1)-k_offset(:)
enddo
if(accumulation_method_mom_den=='DMC')then
allocate(xyzzyaaii1(xyzzyaahs1),xyzzyaaig1(expval_ngvec(1),xyzzyaahs1)&
&,xyzzyaaih1(expval_ngvec(1),xyzzyaahs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOM_DEN','expval_mom_den_*_dmc')
xyzzyaaii1=0.d0
xyzzyaaig1=0.d0
xyzzyaaih1=0.d0
endif
endif
call mpi_bcast(xyzzyaahu1,xyzzyaahs1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast mom_den_nptypes_in_set in setup_moment&
&um_density')
call mpi_bcast(xyzzyaahv1,xyzzyaahs1*nspin,mpi_integer,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'broadcast mom_den_ptype_in_set in setup_momentum&
&_density')
endif
allocate(xyzzyaaif1(xyzzyaahs1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOM_DEN','mom_den_factor')
xyzzyaaif1=0.d0
if(on_top_ii==0)then
select case(periodicity)
case(1)
xyzzyaaae59=xyzzyaacq1/(dble(netot)*twopi)
case(2)
xyzzyaaae59=xyzzyaacp1/(dble(netot)*twopi**2)
case(3)
xyzzyaaae59=xyzzyaacr1/(dble(netot)*twopi**3)
end select
do xyzzyaaad59=1,xyzzyaahs1
do xyzzyaaaa59=1,xyzzyaahu1(xyzzyaaad59)
xyzzyaaif1(xyzzyaaad59)=xyzzyaaif1(xyzzyaaad59)+xyzzyaaae59*real(nele(&
&xyzzyaahv1(xyzzyaaaa59,xyzzyaaad59)),dp)
enddo
enddo
else
select case(periodicity)
case(1)
xyzzyaaif1(1)=xyzzyaacq1/twopi
case(2)
xyzzyaaif1(1)=xyzzyaacp1/twopi**2
case(3)
xyzzyaaif1(1)=xyzzyaacr1/twopi**3
end select
endif
end subroutine xyzzyaari1
subroutine xyzzyaarj1
implicit none
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabx1)then
xyzzyaaki1=3
xyzzyaakj1=2
xyzzyaakq1=20.0d0
xyzzyaakr1=xyzzyaakq1
fin_den_basis='step'
xyzzyaakl1=xyzzyaakm1
xyzzyaakk1=10000
allocate(xyzzyaalg1(xyzzyaakk1,xyzzyaaki1),xyzzyaalh1(xyzzyaakk1,xyzzy&
&aaki1),xyzzyaale1(xyzzyaakk1,xyzzyaaki1),xyzzyaalf1(xyzzyaakk1,xyzzya&
&aki1),xyzzyaali1(xyzzyaakk1,xyzzyaaki1),xyzzyaalj1(xyzzyaakk1,xyzzyaa&
&ki1),xyzzyaalo1(xyzzyaakk1,xyzzyaakj1),xyzzyaalp1(xyzzyaakk1,xyzzyaak&
&j1),xyzzyaalm1(xyzzyaakk1,xyzzyaakj1),xyzzyaaln1(xyzzyaakk1,xyzzyaakj&
&1),xyzzyaalq1(xyzzyaakk1,xyzzyaakj1),xyzzyaalr1(xyzzyaakk1,xyzzyaakj1&
&),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FINITE_DENSITY','expval_fin_den arr&
&ays')
xyzzyaalg1=0.d0
xyzzyaalh1=0.d0
xyzzyaale1=0.d0
xyzzyaalf1=0.d0
xyzzyaali1=0.d0
xyzzyaalj1=0.d0
xyzzyaalo1=0.d0
xyzzyaalp1=0.d0
xyzzyaalm1=0.d0
xyzzyaaln1=0.d0
xyzzyaalq1=0.d0
xyzzyaalr1=0.d0
if(isvmc)then
accumulation_method_fin_den='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_fin_den='DMC'
else
call errstop('SETUP_FINITE_DENSITY','Confusion over what accumulation &
&method I''m using. Should not happen.')
endif
if(accumulation_method_fin_den=='DMC')then
allocate(xyzzyaalk1(xyzzyaakk1,xyzzyaaki1),xyzzyaall1(xyzzyaakk1,xyzzy&
&aaki1),xyzzyaals1(xyzzyaakk1,xyzzyaakj1),xyzzyaalt1(xyzzyaakk1,xyzzya&
&akj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FINITE_DENSITY','expval_fin_den DMC&
& arrays')
xyzzyaalk1=0.d0
xyzzyaall1=0.d0
xyzzyaals1=0.d0
xyzzyaalt1=0.d0
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaaki1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast fin_den_nsets_ij in setup_finite_densi&
&ty')
call mpi_bcast(xyzzyaakj1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast fin_den_nsets_i in setup_finite_densit&
&y')
call mpi_bcast(xyzzyaakk1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast fin_den_norder in setup_finite_density&
&')
call mpi_bcast(xyzzyaakq1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast fin_den_cutoff in setup_finite_density&
&')
call mpi_bcast(xyzzyaakl1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast fin_den_ibasis in setup_finite_density&
&')
call mpi_bcast(accumulation_method_fin_den,3,mpi_character,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'broadcast accumulation_method_fin_den in setup_f&
&inite_density')
if(am_slave)then
allocate(xyzzyaalg1(xyzzyaakk1,xyzzyaaki1),xyzzyaalh1(xyzzyaakk1,xyzzy&
&aaki1),xyzzyaale1(xyzzyaakk1,xyzzyaaki1),xyzzyaalf1(xyzzyaakk1,xyzzya&
&aki1),xyzzyaali1(xyzzyaakk1,xyzzyaaki1),xyzzyaalj1(xyzzyaakk1,xyzzyaa&
&ki1),xyzzyaalo1(xyzzyaakk1,xyzzyaakj1),xyzzyaalp1(xyzzyaakk1,xyzzyaak&
&j1),xyzzyaalm1(xyzzyaakk1,xyzzyaakj1),xyzzyaaln1(xyzzyaakk1,xyzzyaakj&
&1),xyzzyaalq1(xyzzyaakk1,xyzzyaakj1),xyzzyaalr1(xyzzyaakk1,xyzzyaakj1&
&),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FINITE_DENSITY','expval_fin_den arr&
&ays')
xyzzyaalg1=0.d0
xyzzyaalh1=0.d0
xyzzyaale1=0.d0
xyzzyaalf1=0.d0
xyzzyaali1=0.d0
xyzzyaalj1=0.d0
xyzzyaalo1=0.d0
xyzzyaalp1=0.d0
xyzzyaalm1=0.d0
xyzzyaaln1=0.d0
xyzzyaalq1=0.d0
xyzzyaalr1=0.d0
if(accumulation_method_fin_den=='DMC')then
allocate(xyzzyaalk1(xyzzyaakk1,xyzzyaaki1),xyzzyaall1(xyzzyaakk1,xyzzy&
&aaki1),xyzzyaals1(xyzzyaakk1,xyzzyaakj1),xyzzyaalt1(xyzzyaakk1,xyzzya&
&akj1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FINITE_DENSITY','expval_fin_den DMC&
& arrays')
xyzzyaalk1=0.d0
xyzzyaall1=0.d0
xyzzyaals1=0.d0
xyzzyaalt1=0.d0
endif
endif
endif
if(.not.xyzzyaaav1.or..not.xyzzyaabx1)then
xyzzyaakt1=0.d0
xyzzyaakx1=0.d0
xyzzyaalb1=0.d0
endif
xyzzyaaks1=0.d0
xyzzyaaku1=0.d0
xyzzyaakv1=0.d0
xyzzyaakw1=0.d0
xyzzyaaky1=0.d0
xyzzyaakz1=0.d0
xyzzyaala1=0.d0
xyzzyaalc1=0.d0
xyzzyaald1=0.d0
allocate(xyzzyaalu1(xyzzyaakk1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_FINITE_DENSITY','fin_den_phi')
end subroutine xyzzyaarj1
subroutine xyzzyaark1
implicit none
integer xyzzyaaaa61,xyzzyaaab61,xyzzyaaac61
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaabz1)then
xyzzyaamn1=no_families
allocate(xyzzyaamr1(xyzzyaamn1),xyzzyaamv1(nitot,xyzzyaamn1),xyzzyaamw&
&1(nitot,xyzzyaamn1),xyzzyaamq1(xyzzyaamn1),xyzzyaamt1(nitot,xyzzyaamn&
&1),xyzzyaamu1(nitot,xyzzyaamn1),xyzzyaams1(xyzzyaamn1),xyzzyaamx1(nit&
&ot,xyzzyaamn1),xyzzyaamy1(nitot,xyzzyaamn1),xyzzyaamo1(xyzzyaamn1),xy&
&zzyaamp1(nspin,xyzzyaamn1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_POPULATION','expval_population')
xyzzyaamr1=0.d0
xyzzyaamv1=0.d0
xyzzyaamw1=0.d0
xyzzyaamq1=0.d0
xyzzyaamt1=0.d0
xyzzyaamu1=0.d0
xyzzyaams1=0.d0
xyzzyaamx1=0.d0
xyzzyaamy1=0.d0
xyzzyaamo1=0
xyzzyaamp1=0
do xyzzyaaaa61=1,no_families
xyzzyaaac61=0
do xyzzyaaab61=1,nspin
if(which_fam(xyzzyaaab61)==xyzzyaaaa61)then
xyzzyaaac61=xyzzyaaac61+1
xyzzyaamp1(xyzzyaaac61,xyzzyaaaa61)=xyzzyaaab61
endif
enddo
xyzzyaamo1(xyzzyaaaa61)=xyzzyaaac61
enddo
if(isvmc)then
accumulation_method_population='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_population='DMC'
else
call errstop('SETUP_POPULATION','Confusion over what accumulation meth&
&od I''m using. Should not happen.')
endif
if(accumulation_method_population=='DMC')then
allocate(xyzzyaanb1(xyzzyaamn1),xyzzyaamz1(nitot,xyzzyaamn1),xyzzyaana&
&1(nitot,xyzzyaamn1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_POPULATION','expval_population_*_dm&
&c')
xyzzyaanb1=0.d0
xyzzyaamz1=0.d0
xyzzyaana1=0.d0
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaamn1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast population_nsets in SETUP_POPULATION')
call mpi_bcast(nitot,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast nitot in SETUP_POPULATION')
call mpi_bcast(accumulation_method_population,3,mpi_character,0, mpi_c&
&omm_world,ierror)
call checkmpi(ierror,'broadcast accumulation_method_population in SETU&
&P_POPULATION')
if(am_slave)then
allocate(xyzzyaamr1(xyzzyaamn1),xyzzyaamv1(nitot,xyzzyaamn1),xyzzyaamw&
&1(nitot,xyzzyaamn1),xyzzyaamq1(xyzzyaamn1),xyzzyaamt1(nitot,xyzzyaamn&
&1),xyzzyaamu1(nitot,xyzzyaamn1),xyzzyaams1(xyzzyaamn1),xyzzyaamx1(nit&
&ot,xyzzyaamn1),xyzzyaamy1(nitot,xyzzyaamn1),xyzzyaamo1(xyzzyaamn1),xy&
&zzyaamp1(nspin,xyzzyaamn1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_POPULATION','expval_population')
xyzzyaamr1=0.d0
xyzzyaamv1=0.d0
xyzzyaamw1=0.d0
xyzzyaamq1=0.d0
xyzzyaamt1=0.d0
xyzzyaamu1=0.d0
xyzzyaams1=0.d0
xyzzyaamx1=0.d0
xyzzyaamy1=0.d0
xyzzyaamo1=0
xyzzyaamp1=0
if(accumulation_method_population=='DMC')then
allocate(xyzzyaanb1(xyzzyaamn1),xyzzyaamz1(nitot,xyzzyaamn1),xyzzyaana&
&1(nitot,xyzzyaamn1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_POPULATION','expval_population_*_dm&
&c')
xyzzyaanb1=0.d0
xyzzyaamz1=0.d0
xyzzyaana1=0.d0
endif
endif
call mpi_bcast(xyzzyaamo1,xyzzyaamn1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast population_nptypes_in_set in SETUP_POP&
&ULATION')
call mpi_bcast(xyzzyaamp1,xyzzyaamn1*nspin,mpi_integer,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'broadcast population_ptype_in_set in SETUP_POPUL&
&ATION')
endif
end subroutine xyzzyaark1
subroutine xyzzyaarl1
implicit none
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaacb1)then
if(mol_spin_density)then
xyzzyaalv1=nspin
else
xyzzyaalv1=1
endif
xyzzyaalw1=(/20,20,20/)
xyzzyaaly1=(/5.d0,5.d0,5.d0/)
xyzzyaalx1=(/-5.d0,-5.d0,-5.d0/)
allocate(xyzzyaami1(xyzzyaalw1(1),xyzzyaalw1(2),                     x&
&yzzyaalw1(3),xyzzyaalv1),xyzzyaamh1(xyzzyaalw1(1),xyzzyaalw1(2),     &
&          xyzzyaalw1(3),xyzzyaalv1),xyzzyaamj1(xyzzyaalw1(1),xyzzyaal&
&w1(2),                     xyzzyaalw1(3),xyzzyaalv1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOL_DENSITY','expval_mol_den arrays&
&')
xyzzyaami1=0.d0
xyzzyaamh1=0.d0
xyzzyaamj1=0.d0
if(isvmc)then
accumulation_method_mol_den='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_mol_den='DMC'
else
call errstop('SETUP_MOL_DENSITY','Confusion over what accumulation met&
&hod I''m using. Should not happen.')
endif
if(accumulation_method_mol_den=='DMC')then
allocate(xyzzyaamk1(xyzzyaalw1(1),xyzzyaalw1(2),                   xyz&
&zyaalw1(3),xyzzyaalv1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOL_DENSITY','expval_mol_den DMC ar&
&rays')
xyzzyaamk1=0.d0
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaalv1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast mol_den_nsets in setup_mol_density')
call mpi_bcast(xyzzyaalw1,3,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast mol_den_norder in setup_mol_density')
call mpi_bcast(xyzzyaalx1,3,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast mol_den_mincoord in setup_mol_density'&
&)
call mpi_bcast(xyzzyaaly1,3,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast mol_den_maxcoord in setup_mol_density'&
&)
call mpi_bcast(accumulation_method_mol_den,3,mpi_character,0,mpi_comm_&
&world,ierror)
call checkmpi(ierror,'broadcast accumulation_method_mol_den in setup_m&
&ol_density')
if(am_slave)then
allocate(xyzzyaami1(xyzzyaalw1(1),xyzzyaalw1(2),                     x&
&yzzyaalw1(3),xyzzyaalv1),xyzzyaamh1(xyzzyaalw1(1),xyzzyaalw1(2),     &
&          xyzzyaalw1(3),xyzzyaalv1),xyzzyaamj1(xyzzyaalw1(1),xyzzyaal&
&w1(2),                     xyzzyaalw1(3),xyzzyaalv1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOL_DENSITY','expval_mol_den arrays&
&')
xyzzyaami1=0.d0
xyzzyaamh1=0.d0
xyzzyaamj1=0.d0
if(accumulation_method_mol_den=='DMC')then
allocate(xyzzyaamk1(xyzzyaalw1(1),xyzzyaalw1(2),                   xyz&
&zyaalw1(3),xyzzyaalv1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOL_DENSITY','expval_mol_den DMC ar&
&rays')
xyzzyaamk1=0.d0
endif
endif
endif
if(.not.xyzzyaaav1.or..not.xyzzyaacb1)then
xyzzyaama1=0.d0
xyzzyaame1=0.d0
endif
xyzzyaalz1=0.d0
xyzzyaamb1=0.d0
xyzzyaamc1=0.d0
xyzzyaamd1=0.d0
xyzzyaamf1=0.d0
xyzzyaamg1=0.d0
allocate(xyzzyaaml1(xyzzyaalw1(1),xyzzyaalw1(2),xyzzyaalw1(3)),stat=xy&
&zzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_MOL_DENSITY','mol_den_phi')
end subroutine xyzzyaarl1
subroutine xyzzyaarm1
implicit none
integer xyzzyaaaa63,xyzzyaaab63,xyzzyaaac63,xyzzyaaad63,xyzzyaaae63,xy&
&zzyaaaf63,xyzzyaaag63,xyzzyaaah63,xyzzyaaai63,xyzzyaaaj63,xyzzyaaak63&
&,xyzzyaaal63,xyzzyaaam63
logical xyzzyaaan63
real(dp) xyzzyaaao63
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaacd1)then
if(.not.xyzzyaabr1)then
xyzzyaanc1=no_spairs(min(1,levels_spairs))
xyzzyaand1=20
xyzzyaanh1=min(1.d0,2.d0/dble(maxval(nele)))
else
xyzzyaanc1=xyzzyaaik1
xyzzyaand1=xyzzyaaim1
xyzzyaanh1=xyzzyaaiv1
endif
allocate(xyzzyaanj1(xyzzyaanc1),xyzzyaann1(expval_ngvec(1),xyzzyaanc1)&
&,xyzzyaano1(expval_ngvec(1),xyzzyaanc1),xyzzyaani1(xyzzyaanc1),xyzzya&
&anl1(expval_ngvec(1),xyzzyaanc1),xyzzyaanm1(expval_ngvec(1),xyzzyaanc&
&1),xyzzyaank1(xyzzyaanc1),xyzzyaanp1(expval_ngvec(1),xyzzyaanc1),xyzz&
&yaanq1(expval_ngvec(1),xyzzyaanc1),xyzzyaane1(xyzzyaanc1),xyzzyaanf1(&
&2,max_spin_pairs,xyzzyaanc1),xyzzyaanu1(3,expval_ngvec(1)),stat=xyzzy&
&aaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DM_MOM','expval_twop_dm_mom')
xyzzyaanj1=0.d0
xyzzyaann1=0.d0
xyzzyaano1=0.d0
xyzzyaani1=0.d0
xyzzyaanl1=0.d0
xyzzyaanm1=0.d0
xyzzyaank1=0.d0
xyzzyaanp1=0.d0
xyzzyaanq1=0.d0
xyzzyaane1=0
xyzzyaanf1=0
do xyzzyaaaa63=1,expval_ngvec(1)
xyzzyaanu1(:,xyzzyaaaa63)=expval_gvec(:,xyzzyaaaa63,1)-k_offset(:)
enddo
if(.not.xyzzyaabr1)then
do xyzzyaaac63=1,no_spairs(min(1,levels_spairs))
xyzzyaaad63=0
do xyzzyaaaa63=1,nspin
do xyzzyaaab63=xyzzyaaaa63,nspin
if(which_spair(xyzzyaaaa63,xyzzyaaab63,min(1,levels_spairs))==xyzzyaaa&
&c63)then
xyzzyaaad63=xyzzyaaad63+1
xyzzyaanf1(1,xyzzyaaad63,xyzzyaaac63)=xyzzyaaaa63
xyzzyaanf1(2,xyzzyaaad63,xyzzyaaac63)=xyzzyaaab63
endif
enddo
enddo
xyzzyaane1(xyzzyaaac63)=xyzzyaaad63
enddo
if(isvmc)then
accumulation_method_twop_dm_mom='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accumulation_method_twop_dm_mom='DMC'
else
call errstop('SETUP_TWOP_DM_MOM','Confusion over what accumulation met&
&hod I''m using. Should not happen.')
endif
else
xyzzyaanf1=xyzzyaanf1
xyzzyaane1=xyzzyaane1
accumulation_method_twop_dm_mom=accumulation_method_twop_dm_mom
endif
if(accumulation_method_twop_dm_mom=='DMC')then
allocate(xyzzyaant1(xyzzyaanc1),xyzzyaanr1(expval_ngvec(1),xyzzyaanc1)&
&,xyzzyaans1(expval_ngvec(1),xyzzyaanc1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DM_MOM','expval_twop_dm_mom_*_&
&dmc')
xyzzyaant1=0.d0
xyzzyaanr1=0.d0
xyzzyaans1=0.d0
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaanc1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast twop_dm_mom_nsets in setup_twop_dm_mom&
&')
call mpi_bcast(expval_ngvec(1),1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast expval_ngvec(1) in setup_twop_dm_mom')
call mpi_bcast(xyzzyaand1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast twop_dm_mom_nrandom in setup_twop_dm_m&
&om')
call mpi_bcast(xyzzyaanh1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast twop_dm_mom_frac_pairs in setup_twop_d&
&m_mom')
call mpi_bcast(accumulation_method_twop_dm_mom,3,mpi_character,0,mpi_c&
&omm_world,ierror)
call checkmpi(ierror,'broadcast accumulation_method_twop_dm_mom in set&
&up_twop_dm_mom')
if(am_slave)then
allocate(xyzzyaanj1(xyzzyaanc1),xyzzyaann1(expval_ngvec(1),xyzzyaanc1)&
&,xyzzyaano1(expval_ngvec(1),xyzzyaanc1),xyzzyaani1(xyzzyaanc1),xyzzya&
&anl1(expval_ngvec(1),xyzzyaanc1),xyzzyaanm1(expval_ngvec(1),xyzzyaanc&
&1),xyzzyaank1(xyzzyaanc1),xyzzyaanp1(expval_ngvec(1),xyzzyaanc1),xyzz&
&yaanq1(expval_ngvec(1),xyzzyaanc1),xyzzyaane1(xyzzyaanc1),xyzzyaanf1(&
&2,max_spin_pairs,xyzzyaanc1),xyzzyaanu1(3,expval_ngvec(1)),stat=xyzzy&
&aaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DM_MOM','expval_twop_dm_mom')
xyzzyaanj1=0.d0
xyzzyaann1=0.d0
xyzzyaano1=0.d0
xyzzyaani1=0.d0
xyzzyaanl1=0.d0
xyzzyaanm1=0.d0
xyzzyaank1=0.d0
xyzzyaanp1=0.d0
xyzzyaanq1=0.d0
xyzzyaane1=0
xyzzyaanf1=0
do xyzzyaaaa63=1,expval_ngvec(1)
xyzzyaanu1(:,xyzzyaaaa63)=expval_gvec(:,xyzzyaaaa63,1)-k_offset(:)
enddo
if(accumulation_method_twop_dm_mom=='DMC')then
allocate(xyzzyaant1(xyzzyaanc1),xyzzyaanr1(expval_ngvec(1),xyzzyaanc1)&
&,xyzzyaans1(expval_ngvec(1),xyzzyaanc1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DM_MOM','expval_twop_dm_mom_*_&
&dmc')
xyzzyaant1=0.d0
xyzzyaanr1=0.d0
xyzzyaans1=0.d0
endif
endif
call mpi_bcast(xyzzyaane1,xyzzyaanc1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast twop_dm_mom_nptypes_in_set in setup_mo&
&mentum_density')
call mpi_bcast(xyzzyaanf1,xyzzyaanc1*nspin,mpi_integer,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'broadcast twop_dm_mom_ptype_in_set in setup_mome&
&ntum_density')
endif
select case(periodicity)
case(1)
xyzzyaaao63=1.d0/xyzzyaacq1**2
case(2)
xyzzyaaao63=1.d0/xyzzyaacp1**2
case(3)
xyzzyaaao63=1.d0/xyzzyaacr1**2
end select
if(.not.twop_density_mat)then
allocate(xyzzyaaip1(xyzzyaanc1),xyzzyaais1(xyzzyaanc1),xyzzyaanv1(xyzz&
&yaanc1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DM_MOM','twop_dm_npair_in_set'&
&)
do xyzzyaaae63=1,xyzzyaanc1
xyzzyaaaf63=0
xyzzyaaam63=0
do xyzzyaaaa63=1,xyzzyaane1(xyzzyaaae63)
xyzzyaaai63=xyzzyaanf1(1,xyzzyaaaa63,xyzzyaaae63)
xyzzyaaaj63=xyzzyaanf1(2,xyzzyaaaa63,xyzzyaaae63)
if(xyzzyaaai63/=xyzzyaaaj63)then
xyzzyaaaf63=xyzzyaaaf63+nele(xyzzyaaai63)*nele(xyzzyaaaj63)
xyzzyaaam63=xyzzyaaam63+nele(xyzzyaaai63)*nele(xyzzyaaaj63)
else
xyzzyaaaf63=xyzzyaaaf63+(nele(xyzzyaaai63)*(nele(xyzzyaaaj63)-1))/2
xyzzyaaam63=xyzzyaaam63+nele(xyzzyaaai63)*(nele(xyzzyaaaj63)-1)
endif
enddo
xyzzyaaip1(xyzzyaaae63)=xyzzyaaaf63
xyzzyaais1(xyzzyaaae63)=nint(xyzzyaaaf63*xyzzyaanh1)
xyzzyaanv1(xyzzyaaae63)=xyzzyaaao63*dble(xyzzyaaam63)
if(xyzzyaaaf63/=0.and.xyzzyaais1(xyzzyaaae63)<1)call errstop('SETUP_TW&
&OP_DM_MOM','Fraction of particle-pairs too small.')
enddo
allocate(xyzzyaaiq1(2,maxval(xyzzyaaip1),xyzzyaanc1),xyzzyaair1(maxval&
&(xyzzyaais1),xyzzyaanc1),xyzzyaait1(maxval(xyzzyaais1)),stat=xyzzyaaa&
&d1)
xyzzyaair1=0
call check_alloc(xyzzyaaad1,'SETUP_TWOP_DM_MOM','twop_dm_mom_pair')
do xyzzyaaae63=1,xyzzyaanc1
xyzzyaaac63=0
do xyzzyaaaa63=1,xyzzyaane1(xyzzyaaae63)
xyzzyaaai63=xyzzyaanf1(1,xyzzyaaaa63,xyzzyaaae63)
xyzzyaaaj63=xyzzyaanf1(2,xyzzyaaaa63,xyzzyaaae63)
xyzzyaaan63=(xyzzyaaai63==xyzzyaaaj63)
do xyzzyaaak63=1,nele(xyzzyaaai63)
do xyzzyaaal63=1,nele(xyzzyaaaj63)
if(xyzzyaaan63.and.xyzzyaaak63>=xyzzyaaal63)cycle
xyzzyaaac63=xyzzyaaac63+1
xyzzyaaag63=which_ii(xyzzyaaak63,xyzzyaaai63)
xyzzyaaah63=which_ii(xyzzyaaal63,xyzzyaaaj63)
xyzzyaaiq1(1,xyzzyaaac63,xyzzyaaae63)=xyzzyaaag63
xyzzyaaiq1(2,xyzzyaaac63,xyzzyaaae63)=xyzzyaaah63
enddo
enddo
enddo
enddo
if(xyzzyaanh1==1.d0)then
do xyzzyaaae63=1,xyzzyaanc1
do xyzzyaaaa63=1,xyzzyaaip1(xyzzyaaae63)
xyzzyaair1(xyzzyaaaa63,xyzzyaaae63)=xyzzyaaaa63
enddo
enddo
endif
endif
end subroutine xyzzyaarm1
subroutine xyzzyaarn1
implicit none
integer xyzzyaaaa64,xyzzyaaab64,xyzzyaaac64,xyzzyaaad64,xyzzyaaae64,xy&
&zzyaaaf64,xyzzyaaag64,xyzzyaaah64,xyzzyaaai64,xyzzyaaaj64,xyzzyaaak64&
&,xyzzyaaal64,xyzzyaaam64
logical xyzzyaaan64
real(dp) xyzzyaaao64
if(am_master)then
if(.not.xyzzyaaav1.or..not.xyzzyaacf1)then
if(.not.xyzzyaabt1)then
xyzzyaanw1=no_spairs(min(1,levels_spairs))
xyzzyaanx1=20
xyzzyaaob1=min(1.d0,2.d0/dble(maxval(nele)))
else
xyzzyaanw1=xyzzyaajj1
xyzzyaanx1=xyzzyaajl1
xyzzyaaob1=xyzzyaaju1
endif
allocate(xyzzyaaod1(xyzzyaanw1),xyzzyaaoh1(expval_ngvec(1),xyzzyaanw1)&
&,xyzzyaaoi1(expval_ngvec(1),xyzzyaanw1),xyzzyaaoc1(xyzzyaanw1),xyzzya&
&aof1(expval_ngvec(1),xyzzyaanw1),xyzzyaaog1(expval_ngvec(1),xyzzyaanw&
&1),xyzzyaaoe1(xyzzyaanw1),xyzzyaaoj1(expval_ngvec(1),xyzzyaanw1),xyzz&
&yaaok1(expval_ngvec(1),xyzzyaanw1),xyzzyaany1(xyzzyaanw1),xyzzyaanz1(&
&2,max_spin_pairs,xyzzyaanw1),xyzzyaaoo1(3,expval_ngvec(1)),stat=xyzzy&
&aaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRAC_MOM','expval_cond_frac_mo&
&m')
xyzzyaaod1=0.d0
xyzzyaaoh1=0.d0
xyzzyaaoi1=0.d0
xyzzyaaoc1=0.d0
xyzzyaaof1=0.d0
xyzzyaaog1=0.d0
xyzzyaaoe1=0.d0
xyzzyaaoj1=0.d0
xyzzyaaok1=0.d0
xyzzyaany1=0
xyzzyaanz1=0
do xyzzyaaaa64=1,expval_ngvec(1)
xyzzyaaoo1(:,xyzzyaaaa64)=expval_gvec(:,xyzzyaaaa64,1)-k_offset(:)
enddo
if(.not.xyzzyaabt1)then
do xyzzyaaac64=1,no_spairs(min(1,levels_spairs))
xyzzyaaad64=0
do xyzzyaaaa64=1,nspin
do xyzzyaaab64=xyzzyaaaa64,nspin
if(which_spair(xyzzyaaaa64,xyzzyaaab64,min(1,levels_spairs))==xyzzyaaa&
&c64)then
xyzzyaaad64=xyzzyaaad64+1
xyzzyaanz1(1,xyzzyaaad64,xyzzyaaac64)=xyzzyaaaa64
xyzzyaanz1(2,xyzzyaaad64,xyzzyaaac64)=xyzzyaaab64
endif
enddo
enddo
xyzzyaany1(xyzzyaaac64)=xyzzyaaad64
enddo
if(isvmc)then
accum_method_cond_frac_mom='VMC'
elseif(isdmc.or.isvmc_dmc.or.isdmc_dmc)then
accum_method_cond_frac_mom='DMC'
else
call errstop('SETUP_COND_FRAC_MOM','Confusion over what accumulation m&
&ethod I''m using. Should not happen.')
endif
else
xyzzyaanz1=xyzzyaanz1
xyzzyaany1=xyzzyaany1
accum_method_cond_frac_mom=accum_method_cond_frac_mom
endif
if(accum_method_cond_frac_mom=='DMC')then
allocate(xyzzyaaon1(xyzzyaanw1),xyzzyaaol1(expval_ngvec(1),xyzzyaanw1)&
&,xyzzyaaom1(expval_ngvec(1),xyzzyaanw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRAC_MOM','expval_cond_frac_mo&
&m_*_dmc')
xyzzyaaon1=0.d0
xyzzyaaol1=0.d0
xyzzyaaom1=0.d0
endif
endif
endif
if(nnodes>1)then
call mpi_bcast(xyzzyaanw1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast cond_frac_mom_nsets in setup_cond_frac&
&_mom')
call mpi_bcast(expval_ngvec(1),1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast expval_ngvec(1) in setup_cond_frac_mom&
&')
call mpi_bcast(xyzzyaanx1,1,mpi_integer,0,mpi_comm_world,ierror)
call checkmpi(ierror,'broadcast cond_frac_mom_nrandom in setup_cond_fr&
&ac_mom')
call mpi_bcast(xyzzyaaob1,1,mpi_double_precision,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast cond_frac_mom_frac_pairs in setup_cond&
&_frac_mom')
call mpi_bcast(accum_method_cond_frac_mom,3,mpi_character,0,mpi_comm_w&
&orld,ierror)
call checkmpi(ierror,'broadcast accum_method_cond_frac_mom in setup_co&
&nd_frac_mom')
if(am_slave)then
allocate(xyzzyaaod1(xyzzyaanw1),xyzzyaaoh1(expval_ngvec(1),xyzzyaanw1)&
&,xyzzyaaoi1(expval_ngvec(1),xyzzyaanw1),xyzzyaaoc1(xyzzyaanw1),xyzzya&
&aof1(expval_ngvec(1),xyzzyaanw1),xyzzyaaog1(expval_ngvec(1),xyzzyaanw&
&1),xyzzyaaoe1(xyzzyaanw1),xyzzyaaoj1(expval_ngvec(1),xyzzyaanw1),xyzz&
&yaaok1(expval_ngvec(1),xyzzyaanw1),xyzzyaany1(xyzzyaanw1),xyzzyaanz1(&
&2,max_spin_pairs,xyzzyaanw1),xyzzyaaoo1(3,expval_ngvec(1)),stat=xyzzy&
&aaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRAC_MOM','expval_cond_frac_mo&
&m')
xyzzyaaod1=0.d0
xyzzyaaoh1=0.d0
xyzzyaaoi1=0.d0
xyzzyaaoc1=0.d0
xyzzyaaof1=0.d0
xyzzyaaog1=0.d0
xyzzyaaoe1=0.d0
xyzzyaaoj1=0.d0
xyzzyaaok1=0.d0
xyzzyaany1=0
xyzzyaanz1=0
do xyzzyaaaa64=1,expval_ngvec(1)
xyzzyaaoo1(:,xyzzyaaaa64)=expval_gvec(:,xyzzyaaaa64,1)-k_offset(:)
enddo
if(accum_method_cond_frac_mom=='DMC')then
allocate(xyzzyaaon1(xyzzyaanw1),xyzzyaaol1(expval_ngvec(1),xyzzyaanw1)&
&,xyzzyaaom1(expval_ngvec(1),xyzzyaanw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRAC_MOM','expval_cond_frac_mo&
&m_*_dmc')
xyzzyaaon1=0.d0
xyzzyaaol1=0.d0
xyzzyaaom1=0.d0
endif
endif
call mpi_bcast(xyzzyaany1,xyzzyaanw1,mpi_integer,0,mpi_comm_world,ierr&
&or)
call checkmpi(ierror,'broadcast cond_frac_mom_nptypes_in_set in setup_&
&momentum_density')
call mpi_bcast(xyzzyaanz1,xyzzyaanw1*nspin,mpi_integer,0,mpi_comm_worl&
&d,ierror)
call checkmpi(ierror,'broadcast cond_frac_mom_ptype_in_set in setup_mo&
&mentum_density')
endif
select case(periodicity)
case(1)
xyzzyaaao64=1.d0/xyzzyaacq1**2
case(2)
xyzzyaaao64=1.d0/xyzzyaacp1**2
case(3)
xyzzyaaao64=1.d0/xyzzyaacr1**2
end select
if(.not.cond_fraction)then
allocate(xyzzyaajo1(xyzzyaanw1),xyzzyaajr1(xyzzyaanw1),xyzzyaaop1(xyzz&
&yaanw1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'SETUP_COND_FRAC_MOM','cond_frac_npair_in_&
&set')
do xyzzyaaae64=1,xyzzyaanw1
xyzzyaaaf64=0
xyzzyaaam64=0
do xyzzyaaaa64=1,xyzzyaany1(xyzzyaaae64)
xyzzyaaai64=xyzzyaanz1(1,xyzzyaaaa64,xyzzyaaae64)
xyzzyaaaj64=xyzzyaanz1(2,xyzzyaaaa64,xyzzyaaae64)
if(xyzzyaaai64/=xyzzyaaaj64)then
xyzzyaaaf64=xyzzyaaaf64+nele(xyzzyaaai64)*nele(xyzzyaaaj64)
xyzzyaaam64=xyzzyaaam64+nele(xyzzyaaai64)*nele(xyzzyaaaj64)
else
xyzzyaaaf64=xyzzyaaaf64+(nele(xyzzyaaai64)*(nele(xyzzyaaaj64)-1))/2
xyzzyaaam64=xyzzyaaam64+nele(xyzzyaaai64)*(nele(xyzzyaaaj64)-1)
endif
enddo
xyzzyaajo1(xyzzyaaae64)=xyzzyaaaf64
xyzzyaajr1(xyzzyaaae64)=nint(xyzzyaaaf64*xyzzyaaob1)
xyzzyaaop1(xyzzyaaae64)=xyzzyaaao64*dble(xyzzyaaam64)
if(xyzzyaaaf64/=0.and.xyzzyaajr1(xyzzyaaae64)<1)call errstop('SETUP_CO&
&ND_FRAC_MOM','Fraction of particle-pairs too small.')
enddo
allocate(xyzzyaajp1(2,maxval(xyzzyaajo1),xyzzyaanw1),xyzzyaajq1(maxval&
&(xyzzyaajr1),xyzzyaanw1),xyzzyaajs1(maxval(xyzzyaajr1)),stat=xyzzyaaa&
&d1)
xyzzyaajq1=0
call check_alloc(xyzzyaaad1,'SETUP_COND_FRAC_MOM','cond_frac_mom_pair'&
&)
do xyzzyaaae64=1,xyzzyaanw1
xyzzyaaac64=0
do xyzzyaaaa64=1,xyzzyaany1(xyzzyaaae64)
xyzzyaaai64=xyzzyaanz1(1,xyzzyaaaa64,xyzzyaaae64)
xyzzyaaaj64=xyzzyaanz1(2,xyzzyaaaa64,xyzzyaaae64)
xyzzyaaan64=(xyzzyaaai64==xyzzyaaaj64)
do xyzzyaaak64=1,nele(xyzzyaaai64)
do xyzzyaaal64=1,nele(xyzzyaaaj64)
if(xyzzyaaan64.and.xyzzyaaak64>=xyzzyaaal64)cycle
xyzzyaaac64=xyzzyaaac64+1
xyzzyaaag64=which_ii(xyzzyaaak64,xyzzyaaai64)
xyzzyaaah64=which_ii(xyzzyaaal64,xyzzyaaaj64)
xyzzyaajp1(1,xyzzyaaac64,xyzzyaaae64)=xyzzyaaag64
xyzzyaajp1(2,xyzzyaaac64,xyzzyaaae64)=xyzzyaaah64
enddo
enddo
enddo
enddo
if(xyzzyaaob1==1.d0)then
do xyzzyaaae64=1,xyzzyaanw1
do xyzzyaaaa64=1,xyzzyaajo1(xyzzyaaae64)
xyzzyaajq1(xyzzyaaaa64,xyzzyaaae64)=xyzzyaaaa64
enddo
enddo
endif
endif
end subroutine xyzzyaarn1
subroutine expval_scratch_request(is0)
implicit none
integer,intent(inout) :: is0
integer xyzzyaaaa65
if(any(xyzzyaapb1(1:xyzzyaape1)==is0).and.is0/=0)return
xyzzyaape1=xyzzyaape1+1
if(spin_density_mat.or.((onep_density_mat.or.mom_den).and.on_top_ii==0&
&).or.(cond_fraction.or.cond_frac_mom))then
xyzzyaaaa65=0
call scratch_request(ratio1_from=is0,ratio1_to=xyzzyaaaa65)
xyzzyaapc1(xyzzyaape1)=xyzzyaaaa65
xyzzyaapb1(xyzzyaape1)=is0
endif
if((twop_density_mat.or.twop_dm_mom).or.(cond_fraction.or.cond_frac_mo&
&m).or.((onep_density_mat.or.mom_den).and.on_top_ii/=0))then
xyzzyaaaa65=0
call scratch_request(ratio2_from=is0,ratio2_to=xyzzyaaaa65)
xyzzyaapd1(xyzzyaape1)=xyzzyaaaa65
xyzzyaapb1(xyzzyaape1)=is0
endif
end subroutine expval_scratch_request
subroutine setup_expval_accum
implicit none
integer xyzzyaaaa66,xyzzyaaab66,xyzzyaaac66
if(any(xyzzyaapc1/=0))then
allocate(xyzzyaapf1(nscratch))
xyzzyaapf1=0
do xyzzyaaaa66=1,xyzzyaape1
xyzzyaaab66=xyzzyaapb1(xyzzyaaaa66)
call which_scratch(xyzzyaaab66)
xyzzyaaac66=xyzzyaapc1(xyzzyaaaa66)
call which_scratch(xyzzyaaac66)
xyzzyaapf1(xyzzyaaab66)=xyzzyaaac66
enddo
endif
if(any(xyzzyaapd1/=0))then
allocate(xyzzyaapg1(nscratch))
xyzzyaapg1=0
do xyzzyaaaa66=1,xyzzyaape1
xyzzyaaab66=xyzzyaapb1(xyzzyaaaa66)
call which_scratch(xyzzyaaab66)
xyzzyaaac66=xyzzyaapd1(xyzzyaaaa66)
call which_scratch(xyzzyaaac66)
xyzzyaapg1(xyzzyaaab66)=xyzzyaaac66
enddo
endif
end subroutine setup_expval_accum
subroutine finish_expval_accum
implicit none
if(allocated(xyzzyaapf1))deallocate(xyzzyaapf1)
if(allocated(xyzzyaapg1))deallocate(xyzzyaapg1)
xyzzyaapb1=0
xyzzyaapc1=0
xyzzyaapd1=0
xyzzyaape1=0
end subroutine finish_expval_accum
subroutine accumulate_expvals(is0,weight,postweight)
implicit none
integer,intent(in) :: is0
logical,intent(in) :: postweight
real(dp),intent(in) :: weight
integer xyzzyaaaa68,xyzzyaaab68,xyzzyaaac68,xyzzyaaad68,xyzzyaaae68
real(dp) xyzzyaaaf68(3)
call timer('ACCUMULATE_EXPVALS',.true.)
xyzzyaaad68=0
if(allocated(xyzzyaapf1))xyzzyaaad68=xyzzyaapf1(is0)
xyzzyaaae68=0
if(allocated(xyzzyaapg1))xyzzyaaae68=xyzzyaapg1(is0)
call get_rsele(is0)
if(structure_factor)xyzzyaagf1=czero
if(pair_corr_sph.or.structure_factor_sph)call get_eevecs(is0)
if(population)call get_eivecs(is0)
if(xyzzyaaba1)then
if(postweight)then
do xyzzyaaaa68=1,netot
xyzzyaaab68=which_spin(xyzzyaaaa68)
xyzzyaaac68=which_fam(xyzzyaaab68)
xyzzyaaaf68=rele_scr(:,xyzzyaaaa68,is0)
call xyzzyaarp1(xyzzyaaaf68,1)
if(density)call xyzzyaarq1(xyzzyaadg1,expval_complex_den,xyzzyaadn1(1,&
&xyzzyaaac68),weight)
if(spin_density)call xyzzyaarq1(xyzzyaado1,xyzzyaadx1,xyzzyaadw1(1,xyz&
&zyaaab68),weight)
if(spin_density_mat)continue
if(pair_corr)call xyzzyaarq1(xyzzyaaeh1,xyzzyaaeq1,xyzzyaaep1(1,xyzzya&
&aab68),weight)
if(structure_factor)then
call xyzzyaarq1(xyzzyaafr1,.true.,xyzzyaagf1(1,xyzzyaaab68),1.d0)
if(.not.homogeneous_density)call xyzzyaarq1(xyzzyaafr1,xyzzyaagk1,xyzz&
&yaagg1(1,xyzzyaaab68),weight)
endif
enddo
else
do xyzzyaaaa68=1,netot
xyzzyaaab68=which_spin(xyzzyaaaa68)
xyzzyaaac68=which_fam(xyzzyaaab68)
xyzzyaaaf68=rele_scr(:,xyzzyaaaa68,is0)
call xyzzyaarp1(xyzzyaaaf68,1)
if(density)call xyzzyaarq1(xyzzyaadg1,expval_complex_den,xyzzyaadk1(1,&
&xyzzyaaac68),weight)
if(spin_density)call xyzzyaarq1(xyzzyaado1,xyzzyaadx1,xyzzyaadt1(1,xyz&
&zyaaab68),weight)
if(spin_density_mat)call xyzzyaarr1(xyzzyaaaa68,xyzzyaaab68,is0,xyzzya&
&aad68,weight)
if(pair_corr)call xyzzyaarq1(xyzzyaaeh1,xyzzyaaeq1,expval_pcf(1,xyzzya&
&aab68),weight)
if(structure_factor)then
call xyzzyaarq1(xyzzyaafr1,.true.,xyzzyaagf1(1,xyzzyaaab68),1.d0)
if(.not.homogeneous_density)call xyzzyaarq1(xyzzyaafr1,xyzzyaagk1,xyzz&
&yaagg1(1,xyzzyaaab68),weight)
endif
enddo
endif
endif
if(pair_corr_sph)call xyzzyaars1(is0,weight,postweight)
if(structure_factor)call xyzzyaart1(weight,postweight)
if(structure_factor_sph)call xyzzyaaru1(is0,weight,postweight)
if(loc_tensor)call xyzzyaaro1(is0,weight,postweight)
if(population)call xyzzyaarz1(is0,weight,postweight)
if(onep_density_mat.or.mom_den)call xyzzyaarv1(is0,xyzzyaaad68,xyzzyaa&
&ae68,weight,postweight)
if(twop_density_mat.or.twop_dm_mom)call xyzzyaarw1(is0,xyzzyaaae68,wei&
&ght,postweight)
if(cond_fraction.or.cond_frac_mom)call xyzzyaarx1(is0,xyzzyaaad68,xyzz&
&yaaae68,weight,postweight)
if(finite_density)call xyzzyaary1(is0,weight,postweight)
if(mol_density)call xyzzyaasa1(is0,weight,postweight)
call timer('ACCUMULATE_EXPVALS',.false.)
end subroutine accumulate_expvals
subroutine xyzzyaaro1(is0,weight,postweight)
implicit none
integer,intent(in) :: is0
real(dp),intent(in),optional :: weight
logical,intent(in) :: postweight
complex(dp) xyzzyaaaa69(periodicity,periodicity)
integer xyzzyaaab69,xyzzyaaac69,xyzzyaaad69,xyzzyaaae69,xyzzyaaaf69
real(dp) xyzzyaaag69(3),xyzzyaaah69(3)
do xyzzyaaad69=1,nspin
xyzzyaaag69=0.d0
do xyzzyaaae69=1,nele(xyzzyaaad69)
xyzzyaaaf69=which_ii(xyzzyaaae69,xyzzyaaad69)
xyzzyaaag69(:)=xyzzyaaag69(:)+rele_scr(:,xyzzyaaaf69,is0)
enddo
do xyzzyaaab69=1,periodicity
xyzzyaaah69(xyzzyaaab69)=ddot(periodicity,xyzzyaaag69(1),1,bmat(1,xyzz&
&yaaab69),1)
enddo
xyzzyaaaa69=czero
do xyzzyaaab69=1,periodicity
xyzzyaaaa69(xyzzyaaab69,xyzzyaaab69)=exp(cmplx(0.d0,xyzzyaaah69(xyzzya&
&aab69),dp))
do xyzzyaaac69=xyzzyaaab69+1,periodicity
xyzzyaaaa69(xyzzyaaac69,xyzzyaaab69)=exp(cmplx(0.d0,xyzzyaaah69(xyzzya&
&aac69)-xyzzyaaah69(xyzzyaaab69),dp))
xyzzyaaaa69(xyzzyaaab69,xyzzyaaac69)=conjg(xyzzyaaaa69(xyzzyaaac69,xyz&
&zyaaab69))
enddo
enddo
if(postweight)then
do xyzzyaaab69=1,periodicity
xyzzyaafo1(xyzzyaaab69,xyzzyaaab69,xyzzyaaad69)=xyzzyaafo1(xyzzyaaab69&
&,xyzzyaaab69,xyzzyaaad69)+xyzzyaaaa69(xyzzyaaab69,xyzzyaaab69)*weight
do xyzzyaaac69=xyzzyaaab69+1,periodicity
xyzzyaafo1(xyzzyaaac69,xyzzyaaab69,xyzzyaaad69)=xyzzyaafo1(xyzzyaaac69&
&,xyzzyaaab69,xyzzyaaad69)+xyzzyaaaa69(xyzzyaaac69,xyzzyaaab69)*weight
xyzzyaafo1(xyzzyaaab69,xyzzyaaac69,xyzzyaaad69)=xyzzyaafo1(xyzzyaaac69&
&,xyzzyaaab69,xyzzyaaad69)
enddo
enddo
xyzzyaafk1(xyzzyaaad69)=xyzzyaafk1(xyzzyaaad69)+weight
else
do xyzzyaaab69=1,periodicity
xyzzyaafl1(xyzzyaaab69,xyzzyaaab69,xyzzyaaad69)=xyzzyaafl1(xyzzyaaab69&
&,xyzzyaaab69,xyzzyaaad69)+xyzzyaaaa69(xyzzyaaab69,xyzzyaaab69)*weight
do xyzzyaaac69=xyzzyaaab69+1,periodicity
xyzzyaafl1(xyzzyaaac69,xyzzyaaab69,xyzzyaaad69)=xyzzyaafl1(xyzzyaaac69&
&,xyzzyaaab69,xyzzyaaad69)+xyzzyaaaa69(xyzzyaaac69,xyzzyaaab69)*weight
xyzzyaafl1(xyzzyaaab69,xyzzyaaac69,xyzzyaaad69)=xyzzyaafl1(xyzzyaaac69&
&,xyzzyaaab69,xyzzyaaad69)
enddo
enddo
xyzzyaafh1(xyzzyaaad69)=xyzzyaafh1(xyzzyaaad69)+weight
endif
enddo
end subroutine xyzzyaaro1
subroutine xyzzyaarp1(rvec,b)
implicit none
integer,intent(in) :: b
real(dp),intent(in) :: rvec(3)
integer xyzzyaaaa70,xyzzyaaab70
real(dp) xyzzyaaac70(3),xyzzyaaad70,xyzzyaaae70,xyzzyaaaf70,xyzzyaaag7&
&0,xyzzyaaah70,xyzzyaaai70
call timer('EXPVAL_FOURIER_BASIS',.true.)
do xyzzyaaab70=1,xyzzyaacw1
xyzzyaaap1(1:3)=xyzzyaadd1(1:3,1,xyzzyaaab70)
xyzzyaaaq1(1:3)=xyzzyaadd1(1:3,2,xyzzyaaab70)
xyzzyaaar1(1:3)=xyzzyaadd1(1:3,3,xyzzyaaab70)
xyzzyaaac70(1)=xyzzyaaap1(1)*rvec(1)+xyzzyaaap1(2)*rvec(2)+xyzzyaaap1(&
&3)*rvec(3)
xyzzyaaac70(2)=xyzzyaaaq1(1)*rvec(1)+xyzzyaaaq1(2)*rvec(2)+xyzzyaaaq1(&
&3)*rvec(3)
xyzzyaaac70(3)=xyzzyaaar1(1)*rvec(1)+xyzzyaaar1(2)*rvec(2)+xyzzyaaar1(&
&3)*rvec(3)
xyzzyaaad70=cos(xyzzyaaac70(1))
xyzzyaaae70=cos(xyzzyaaac70(2))
xyzzyaaaf70=cos(xyzzyaaac70(3))
xyzzyaaag70=sin(xyzzyaaac70(1))
xyzzyaaah70=sin(xyzzyaaac70(2))
xyzzyaaai70=sin(xyzzyaaac70(3))
xyzzyaade1(1:3,0,xyzzyaaab70)=cmplx(1.d0,0.d0,kind=dp)
xyzzyaade1(1,1,xyzzyaaab70)=cmplx(xyzzyaaad70,xyzzyaaag70,kind=dp)
xyzzyaade1(2,1,xyzzyaaab70)=cmplx(xyzzyaaae70,xyzzyaaah70,kind=dp)
xyzzyaade1(3,1,xyzzyaaab70)=cmplx(xyzzyaaaf70,xyzzyaaai70,kind=dp)
xyzzyaade1(1,-1,xyzzyaaab70)=cmplx(xyzzyaaad70,-xyzzyaaag70,kind=dp)
xyzzyaade1(2,-1,xyzzyaaab70)=cmplx(xyzzyaaae70,-xyzzyaaah70,kind=dp)
xyzzyaade1(3,-1,xyzzyaaab70)=cmplx(xyzzyaaaf70,-xyzzyaaai70,kind=dp)
do xyzzyaaaa70=2,xyzzyaacz1(xyzzyaaab70)
xyzzyaade1(1,xyzzyaaaa70,xyzzyaaab70)=xyzzyaade1(1,1,xyzzyaaab70)*xyzz&
&yaade1(1,xyzzyaaaa70-1,xyzzyaaab70)
xyzzyaade1(2,xyzzyaaaa70,xyzzyaaab70)=xyzzyaade1(2,1,xyzzyaaab70)*xyzz&
&yaade1(2,xyzzyaaaa70-1,xyzzyaaab70)
xyzzyaade1(3,xyzzyaaaa70,xyzzyaaab70)=xyzzyaade1(3,1,xyzzyaaab70)*xyzz&
&yaade1(3,xyzzyaaaa70-1,xyzzyaaab70)
xyzzyaade1(1,-xyzzyaaaa70,xyzzyaaab70)=conjg(xyzzyaade1(1,xyzzyaaaa70,&
&xyzzyaaab70))
xyzzyaade1(2,-xyzzyaaaa70,xyzzyaaab70)=conjg(xyzzyaade1(2,xyzzyaaaa70,&
&xyzzyaaab70))
xyzzyaade1(3,-xyzzyaaaa70,xyzzyaaab70)=conjg(xyzzyaade1(3,xyzzyaaaa70,&
&xyzzyaaab70))
enddo
do xyzzyaaaa70=1,expval_ngvec(xyzzyaaab70)
xyzzyaaaa1=xyzzyaacy1(1,xyzzyaaaa70,xyzzyaaab70)
xyzzyaaab1=xyzzyaacy1(2,xyzzyaaaa70,xyzzyaaab70)
xyzzyaaac1=xyzzyaacy1(3,xyzzyaaaa70,xyzzyaaab70)
xyzzyaadf1(xyzzyaaaa70,xyzzyaaab70,b)=xyzzyaade1(1,xyzzyaaaa1,xyzzyaaa&
&b70)*xyzzyaade1(2,xyzzyaaab1,xyzzyaaab70)*xyzzyaade1(3,xyzzyaaac1,xyz&
&zyaaab70)
enddo
enddo
call timer('EXPVAL_FOURIER_BASIS',.false.)
end subroutine xyzzyaarp1
subroutine expval_zero_postweight
implicit none
call timer('EXPVAL_ZERO_POSTWEIGHT',.true.)
if(density)xyzzyaadn1(:,:)=czero
if(spin_density)xyzzyaadw1(:,:)=czero
if(pair_corr)xyzzyaaep1(:,:)=czero
if(pair_corr_sph)then
xyzzyaaff1(:,:,:)=0.d0
xyzzyaafa1=0.d0
endif
if(loc_tensor)then
xyzzyaafo1(:,:,:)=czero
xyzzyaafk1(:)=0.d0
endif
if(structure_factor)then
xyzzyaagc1(:,:)=0.d0
xyzzyaafx1=0.d0
if(.not.homogeneous_density)xyzzyaagj1(:,:)=czero
endif
if(structure_factor_sph)then
xyzzyaagw1(:,:)=0.d0
xyzzyaagr1(:)=0.d0
endif
if(onep_density_mat)then
xyzzyaahp1(:,:)=0.d0
xyzzyaahq1(:,:)=0.d0
xyzzyaahr1(:,:)=0.d0
endif
if(twop_density_mat)then
xyzzyaajf1(:,:)=0.d0
xyzzyaajg1(:,:)=0.d0
xyzzyaajh1(:,:)=0.d0
endif
if(cond_fraction)then
xyzzyaake1(:,:)=0.d0
xyzzyaakf1(:,:)=0.d0
xyzzyaakg1(:,:)=0.d0
endif
if(mom_den)then
xyzzyaaig1(:,:)=0.d0
xyzzyaaih1(:,:)=0.d0
xyzzyaaii1(:)=0.d0
endif
if(finite_density)then
xyzzyaalk1(:,:)=0.d0
xyzzyaall1(:,:)=0.d0
xyzzyaals1(:,:)=0.d0
xyzzyaalt1(:,:)=0.d0
xyzzyaakv1=0.d0
xyzzyaakz1=0.d0
xyzzyaald1=0.d0
endif
if(population)then
xyzzyaamz1(:,:)=0.d0
xyzzyaana1(:,:)=0.d0
xyzzyaanb1(:)=0.d0
endif
if(mol_density)then
xyzzyaamk1(:,:,:,:)=0.d0
xyzzyaamc1=0.d0
xyzzyaamg1=0.d0
endif
if(twop_dm_mom)then
xyzzyaanr1(:,:)=0.d0
xyzzyaans1(:,:)=0.d0
xyzzyaant1(:)=0.d0
endif
if(cond_frac_mom)then
xyzzyaaol1(:,:)=0.d0
xyzzyaaom1(:,:)=0.d0
xyzzyaaon1(:)=0.d0
endif
call timer('EXPVAL_ZERO_POSTWEIGHT',.false.)
end subroutine expval_zero_postweight
subroutine expval_postweight(weight)
implicit none
real(dp),intent(in) :: weight
call timer('EXPVAL_POSTWEIGHT',.true.)
if(density)xyzzyaadk1(:,:)=xyzzyaadk1(:,:)+weight*xyzzyaadn1(:,:)
if(spin_density)xyzzyaadt1(:,:)=xyzzyaadt1(:,:)+weight*xyzzyaadw1(:,:)
if(pair_corr)expval_pcf(:,:)=expval_pcf(:,:)+weight*xyzzyaaep1(:,:)
if(pair_corr_sph)then
xyzzyaafc1(:,:,:)=xyzzyaafc1(:,:,:)+weight*xyzzyaaff1(:,:,:)
xyzzyaaey1=xyzzyaaey1+weight*xyzzyaafa1
endif
if(loc_tensor)then
xyzzyaafl1(:,:,:)=xyzzyaafl1(:,:,:)+weight*xyzzyaafo1(:,:,:)
xyzzyaafh1(:)=xyzzyaafh1(:)+weight*xyzzyaafk1(:)
endif
if(structure_factor)then
xyzzyaaga1(:,:)=xyzzyaaga1(:,:)+weight*xyzzyaagc1(:,:)
xyzzyaafv1=xyzzyaafv1+weight*xyzzyaafx1
if(.not.homogeneous_density)xyzzyaagg1(:,:)=xyzzyaagg1(:,:)+weight*xyz&
&zyaagj1(:,:)
endif
if(structure_factor_sph)then
xyzzyaagu1(:,:)=xyzzyaagu1(:,:)+weight*xyzzyaagw1(:,:)
xyzzyaagq1(:)=xyzzyaagq1(:)+weight*xyzzyaagr1(:)
endif
if(onep_density_mat)then
xyzzyaahi1(:,:)=xyzzyaahi1(:,:)+weight*xyzzyaahp1(:,:)
xyzzyaahj1(:,:)=xyzzyaahj1(:,:)+weight*xyzzyaahq1(:,:)
xyzzyaahf1(:,:)=xyzzyaahf1(:,:)+weight*xyzzyaahr1(:,:)
endif
if(twop_density_mat)then
xyzzyaaiz1(:,:)=xyzzyaaiz1(:,:)+weight*xyzzyaajf1(:,:)
xyzzyaaja1(:,:)=xyzzyaaja1(:,:)+weight*xyzzyaajg1(:,:)
xyzzyaaiw1(:,:)=xyzzyaaiw1(:,:)+weight*xyzzyaajh1(:,:)
endif
if(cond_fraction)then
xyzzyaajy1(:,:)=xyzzyaajy1(:,:)+weight*xyzzyaake1(:,:)
xyzzyaajz1(:,:)=xyzzyaajz1(:,:)+weight*xyzzyaakf1(:,:)
xyzzyaajv1(:,:)=xyzzyaajv1(:,:)+weight*xyzzyaakg1(:,:)
endif
if(mom_den)then
xyzzyaahz1(:,:)=xyzzyaahz1(:,:)+weight*xyzzyaaig1(:,:)
xyzzyaaia1(:,:)=xyzzyaaia1(:,:)+weight*xyzzyaaih1(:,:)
xyzzyaahw1(:)=xyzzyaahw1(:)+weight*xyzzyaaii1(:)
endif
if(finite_density)then
xyzzyaale1(:,:)=xyzzyaale1(:,:)+weight*xyzzyaalk1(:,:)
xyzzyaalf1(:,:)=xyzzyaalf1(:,:)+weight*xyzzyaall1(:,:)
xyzzyaalm1(:,:)=xyzzyaalm1(:,:)+weight*xyzzyaals1(:,:)
xyzzyaaln1(:,:)=xyzzyaaln1(:,:)+weight*xyzzyaalt1(:,:)
xyzzyaaks1=xyzzyaaks1+xyzzyaakv1
xyzzyaakw1=xyzzyaakw1+weight*xyzzyaakz1
xyzzyaala1=xyzzyaala1+weight*xyzzyaald1
endif
if(population)then
call daxpy(nitot*xyzzyaamn1,weight,xyzzyaamz1(1,1),1, xyzzyaamt1(1,1),&
&1)
call daxpy(nitot*xyzzyaamn1,weight,xyzzyaana1(1,1),1, xyzzyaamu1(1,1),&
&1)
call daxpy(xyzzyaamn1,weight,xyzzyaanb1(1),1, xyzzyaamq1(1),1)
endif
if(mol_density)then
xyzzyaamh1(:,:,:,:)=xyzzyaamh1(:,:,:,:)+weight*xyzzyaamk1(:,:,:,:)
xyzzyaalz1=xyzzyaalz1+xyzzyaamc1
xyzzyaamd1=xyzzyaamd1+weight*xyzzyaamg1
endif
if(twop_dm_mom)then
xyzzyaanl1(:,:)=xyzzyaanl1(:,:)+weight*xyzzyaanr1(:,:)
xyzzyaanm1(:,:)=xyzzyaanm1(:,:)+weight*xyzzyaans1(:,:)
xyzzyaani1(:)=xyzzyaani1(:)+weight*xyzzyaant1(:)
endif
if(cond_frac_mom)then
xyzzyaaof1(:,:)=xyzzyaaof1(:,:)+weight*xyzzyaaol1(:,:)
xyzzyaaog1(:,:)=xyzzyaaog1(:,:)+weight*xyzzyaaom1(:,:)
xyzzyaaoc1(:)=xyzzyaaoc1(:)+weight*xyzzyaaon1(:)
endif
call timer('EXPVAL_POSTWEIGHT',.false.)
end subroutine expval_postweight
subroutine xyzzyaarq1(gset,complex_mode,cdata,rweight,cweight)
implicit none
integer,intent(in) :: gset
logical,intent(in) :: complex_mode
real(dp),intent(in),optional :: rweight
complex(dp),intent(in),optional :: cweight
complex(dp),intent(inout) :: cdata(*)
integer xyzzyaaaa73
real(dp) xyzzyaaab73
if(.not.present(cweight))then
if(.not.present(rweight))then
if(complex_mode)then
do xyzzyaaaa73=1,expval_ngvec(gset)
cdata(xyzzyaaaa73)=cdata(xyzzyaaaa73)+xyzzyaadf1(xyzzyaaaa73,gset,1)
enddo
else
do xyzzyaaaa73=1,expval_ngvec(gset)
xyzzyaaab73=real(cdata(xyzzyaaaa73),dp)+real(xyzzyaadf1(xyzzyaaaa73,gs&
&et,1),dp)
cdata(xyzzyaaaa73)=cmplx(xyzzyaaab73,0.d0,dp)
enddo
endif
else
if(complex_mode)then
do xyzzyaaaa73=1,expval_ngvec(gset)
cdata(xyzzyaaaa73)=cdata(xyzzyaaaa73)+rweight*xyzzyaadf1(xyzzyaaaa73,g&
&set,1)
enddo
else
do xyzzyaaaa73=1,expval_ngvec(gset)
xyzzyaaab73=real(cdata(xyzzyaaaa73),dp)+rweight*real(xyzzyaadf1(xyzzya&
&aaa73,gset,1),dp)
cdata(xyzzyaaaa73)=cmplx(xyzzyaaab73,0.d0,dp)
enddo
endif
endif
else
if(complex_mode)then
do xyzzyaaaa73=1,expval_ngvec(gset)
cdata(xyzzyaaaa73)=cdata(xyzzyaaaa73)+cweight*xyzzyaadf1(xyzzyaaaa73,g&
&set,1)
enddo
else
do xyzzyaaaa73=1,expval_ngvec(gset)
xyzzyaaab73=real(cdata(xyzzyaaaa73),dp)+real(cweight*xyzzyaadf1(xyzzya&
&aaa73,gset,1),dp)
cdata(xyzzyaaaa73)=cmplx(xyzzyaaab73,0.d0,dp)
enddo
endif
endif
end subroutine xyzzyaarq1
subroutine xyzzyaarr1(ii,ispin,is,js,weight)
implicit none
integer,intent(in) :: ii,ispin,is,js
real(dp),intent(in) :: weight
integer xyzzyaaaa74
logical isnan,isinf
complex(dp) xyzzyaaab74
xyzzyaaaa74=sele_scr(ii,is)
call define_config_oneelec(ii,is,js,rele_scr(1,ii,is),3-xyzzyaaaa74)
call wfn_ratio(is,js,0,ratio=xyzzyaaab74,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)xyzzyaaab74=0.d0
xyzzyaaab74=xyzzyaaab74*weight
if(xyzzyaaaa74==1)then
call xyzzyaarq1(xyzzyaady1,xyzzyaaeg1,xyzzyaaed1(1,1,ispin),weight)
call xyzzyaarq1(xyzzyaady1,xyzzyaaeg1,xyzzyaaed1(1,3,ispin),cweight=xy&
&zzyaaab74)
else
call xyzzyaarq1(xyzzyaady1,xyzzyaaeg1,xyzzyaaed1(1,4,ispin),weight)
call xyzzyaarq1(xyzzyaady1,xyzzyaaeg1,xyzzyaaed1(1,2,ispin),cweight=xy&
&zzyaaab74)
endif
end subroutine xyzzyaarr1
subroutine xyzzyaars1(is,weight,postweight)
implicit none
integer,intent(in) :: is
real(dp),intent(in) :: weight
logical,intent(in) :: postweight
integer xyzzyaaaa75,xyzzyaaab75,xyzzyaaac75,xyzzyaaad75,xyzzyaaae75,xy&
&zzyaaaf75,xyzzyaaag75,xyzzyaaah75
real(dp) xyzzyaaai75
if(xyzzyaaet1==1)call errstop_master('ACCUMULATE_PCF_SPH','Routine not&
& yet coded for inhomogeneous case.')
if(.not.postweight)then
do xyzzyaaaa75=1,netot
xyzzyaaab75=which_spin(xyzzyaaaa75)
xyzzyaaac75=0
do xyzzyaaae75=1,nspin
xyzzyaaag75=min(xyzzyaaab75,xyzzyaaae75)
xyzzyaaah75=max(xyzzyaaab75,xyzzyaaae75)
do xyzzyaaad75=1,nele(xyzzyaaae75)
xyzzyaaac75=xyzzyaaac75+1
if(xyzzyaaac75==xyzzyaaaa75)cycle
xyzzyaaai75=eevecs_scr(4,xyzzyaaac75,xyzzyaaaa75,is)*xyzzyaaex1
xyzzyaaaf75=int(xyzzyaaai75)+1
if(xyzzyaaaf75>xyzzyaaes1)cycle
xyzzyaafc1(xyzzyaaaf75,xyzzyaaag75,xyzzyaaah75)=xyzzyaafc1(xyzzyaaaf75&
&,xyzzyaaag75,xyzzyaaah75)+weight
enddo
enddo
enddo
xyzzyaaey1=xyzzyaaey1+weight
else
do xyzzyaaaa75=1,netot
xyzzyaaab75=which_spin(xyzzyaaaa75)
xyzzyaaac75=0
do xyzzyaaae75=1,nspin
xyzzyaaag75=min(xyzzyaaab75,xyzzyaaae75)
xyzzyaaah75=max(xyzzyaaab75,xyzzyaaae75)
do xyzzyaaad75=1,nele(xyzzyaaae75)
xyzzyaaac75=xyzzyaaac75+1
if(xyzzyaaac75==xyzzyaaaa75)cycle
xyzzyaaai75=eevecs_scr(4,xyzzyaaac75,xyzzyaaaa75,is)*xyzzyaaex1
xyzzyaaaf75=int(xyzzyaaai75)+1
if(xyzzyaaaf75>xyzzyaaes1)cycle
xyzzyaaff1(xyzzyaaaf75,xyzzyaaag75,xyzzyaaah75)=xyzzyaaff1(xyzzyaaaf75&
&,xyzzyaaag75,xyzzyaaah75)+weight
enddo
enddo
enddo
xyzzyaafa1=xyzzyaafa1+weight
endif
end subroutine xyzzyaars1
subroutine xyzzyaart1(weight,postweight)
implicit none
logical,intent(in) :: postweight
real(dp),intent(in) :: weight
integer xyzzyaaaa76,xyzzyaaab76,xyzzyaaac76,xyzzyaaad76
if(.not.postweight)then
xyzzyaaad76=0
do xyzzyaaaa76=1,nspin
do xyzzyaaab76=xyzzyaaaa76,nspin
xyzzyaaad76=xyzzyaaad76+1
xyzzyaaga1(1,xyzzyaaad76)=xyzzyaaga1(1,xyzzyaaad76)+weight*real(xyzzya&
&agf1(1,xyzzyaaaa76)*xyzzyaagf1(1,xyzzyaaab76),dp)
do xyzzyaaac76=2,expval_ngvec(xyzzyaafr1),2
xyzzyaaga1(xyzzyaaac76,xyzzyaaad76)=xyzzyaaga1(xyzzyaaac76,xyzzyaaad76&
&)+weight*real(xyzzyaagf1(xyzzyaaac76,xyzzyaaaa76)*xyzzyaagf1(xyzzyaaa&
&c76+1,xyzzyaaab76),dp)
xyzzyaaga1(xyzzyaaac76+1,xyzzyaaad76)=xyzzyaaga1(xyzzyaaac76+1,xyzzyaa&
&ad76)+weight*real(xyzzyaagf1(xyzzyaaac76+1,xyzzyaaaa76)*xyzzyaagf1(xy&
&zzyaaac76,xyzzyaaab76),dp)
enddo
enddo
enddo
xyzzyaafv1=xyzzyaafv1+weight
else
xyzzyaaad76=0
do xyzzyaaaa76=1,nspin
do xyzzyaaab76=xyzzyaaaa76,nspin
xyzzyaaad76=xyzzyaaad76+1
xyzzyaagc1(1,xyzzyaaad76)=xyzzyaagc1(1,xyzzyaaad76)+weight*real(xyzzya&
&agf1(1,xyzzyaaaa76)*xyzzyaagf1(1,xyzzyaaab76),dp)
do xyzzyaaac76=2,expval_ngvec(xyzzyaafr1),2
xyzzyaagc1(xyzzyaaac76,xyzzyaaad76)=xyzzyaagc1(xyzzyaaac76,xyzzyaaad76&
&)+weight*real(xyzzyaagf1(xyzzyaaac76,xyzzyaaaa76)*xyzzyaagf1(xyzzyaaa&
&c76+1,xyzzyaaab76),dp)
xyzzyaagc1(xyzzyaaac76+1,xyzzyaaad76)=xyzzyaagc1(xyzzyaaac76+1,xyzzyaa&
&ad76)+weight*real(xyzzyaagf1(xyzzyaaac76+1,xyzzyaaaa76)*xyzzyaagf1(xy&
&zzyaaac76,xyzzyaaab76),dp)
enddo
enddo
enddo
xyzzyaafx1=xyzzyaafx1+weight
endif
end subroutine xyzzyaart1
subroutine xyzzyaaru1(is,weight,postweight)
implicit none
integer,intent(in) :: is
real(dp),intent(in) :: weight
logical,intent(in) :: postweight
integer xyzzyaaaa77,xyzzyaaab77,xyzzyaaac77,xyzzyaaad77,xyzzyaaae77,xy&
&zzyaaaf77
real(dp) xyzzyaaag77,xyzzyaaah77,xyzzyaaai77
do xyzzyaaaa77=1,netot
xyzzyaaab77=which_spin(xyzzyaaaa77)
do xyzzyaaac77=xyzzyaaaa77+1,netot
xyzzyaaad77=which_spin(xyzzyaaac77)
xyzzyaaaf77=which_spair(xyzzyaaad77,xyzzyaaab77,levels_spairs)
xyzzyaaag77=eevecs_scr(4,xyzzyaaac77,xyzzyaaaa77,is)
if(xyzzyaaag77<=wigner_seitz_radius)then
do xyzzyaaae77=1,xyzzyaagl1
xyzzyaaah77=xyzzyaagt1(xyzzyaaae77)*xyzzyaaag77
if(xyzzyaaah77==0.d0)then
xyzzyaaai77=1.d0
else
if(dimensionality==3)then
xyzzyaaai77=sin(xyzzyaaah77)/xyzzyaaah77
elseif(dimensionality==2)then
xyzzyaaai77=bessel_j0(xyzzyaaah77)
else
xyzzyaaai77=cos(xyzzyaaah77)
endif
endif
if(postweight)then
xyzzyaagw1(xyzzyaaae77,xyzzyaaaf77)=xyzzyaagw1(xyzzyaaae77,xyzzyaaaf77&
&)+xyzzyaaai77*weight
else
xyzzyaagu1(xyzzyaaae77,xyzzyaaaf77)=xyzzyaagu1(xyzzyaaae77,xyzzyaaaf77&
&)+xyzzyaaai77*weight
endif
enddo
endif
enddo
enddo
xyzzyaagq1(1:xyzzyaagm1)=xyzzyaagq1(1:xyzzyaagm1)+weight*netot
end subroutine xyzzyaaru1
subroutine xyzzyaarv1(is,js,ks,weight,postweight)
implicit none
integer,intent(in) :: is,js,ks
real(dp),intent(in) :: weight
logical,intent(in) :: postweight
integer xyzzyaaaa78,xyzzyaaab78,xyzzyaaac78,xyzzyaaad78,xyzzyaaae78,xy&
&zzyaaaf78,xyzzyaaag78,xyzzyaaah78
integer,save :: xyzzyaaai78,xyzzyaaaj78
integer,allocatable,save :: xyzzyaaak78(:),xyzzyaaal78(:,:)
real(dp) xyzzyaaam78,xyzzyaaan78(3),xyzzyaaao78,xyzzyaaap78,xyzzyaaaq7&
&8
real(dp),allocatable,save :: xyzzyaaar78(:,:)
complex(dp) xyzzyaaas78,xyzzyaaat78,xyzzyaaau78
logical isnan,isinf
logical,allocatable,save :: xyzzyaaav78(:)
logical,save :: xyzzyaaaw78=.true.
if(xyzzyaaaw78)then
if(onep_density_mat)then
xyzzyaaai78=xyzzyaahb1
xyzzyaaaj78=xyzzyaagz1
else
xyzzyaaai78=xyzzyaaht1
xyzzyaaaj78=xyzzyaahs1
endif
allocate(xyzzyaaar78(3,2*xyzzyaaai78),xyzzyaaak78(2*xyzzyaaai78),xyzzy&
&aaal78(nspin,xyzzyaaaj78),xyzzyaaav78(2*xyzzyaaai78),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'ACCUMULATE_ONEP_DM','0.7')
xyzzyaaal78=0
if(onep_density_mat)then
xyzzyaaal78=xyzzyaahd1
else
xyzzyaaal78=xyzzyaahv1
endif
xyzzyaaaw78=.false.
endif
if(.not.mom_den)then
xyzzyaaav78=.true.
xyzzyaaah78=0
do xyzzyaaac78=1,xyzzyaaai78
call xyzzyaasc1(xyzzyaaar78(1:3,xyzzyaaac78),xyzzyaaam78)
xyzzyaaak78(xyzzyaaac78)=int(xyzzyaaam78*xyzzyaahe1)+1
enddo
elseif(.not.onep_density_mat)then
xyzzyaaav78=.false.
xyzzyaaah78=0
do xyzzyaaac78=1,xyzzyaaai78
call xyzzyaasd1(xyzzyaaar78(1:3,xyzzyaaac78),xyzzyaaam78,xyzzyaaav78(x&
&yzzyaaac78))
enddo
else
xyzzyaaav78=.false.
do xyzzyaaac78=1,xyzzyaaai78
call xyzzyaasd1(xyzzyaaar78(1:3,xyzzyaaac78),xyzzyaaam78,xyzzyaaav78(x&
&yzzyaaac78))
xyzzyaaak78(xyzzyaaac78)=int(xyzzyaaam78*xyzzyaahe1)+1
enddo
xyzzyaaah78=count(.not.xyzzyaaav78(1:xyzzyaaai78))
do xyzzyaaac78=xyzzyaaai78+1,xyzzyaaai78+xyzzyaaah78
call xyzzyaasc1(xyzzyaaar78(1:3,xyzzyaaac78),xyzzyaaam78)
xyzzyaaak78(xyzzyaaac78)=int(xyzzyaaam78*xyzzyaahe1)+1
enddo
endif
do xyzzyaaaa78=1,xyzzyaaaj78
do xyzzyaaaf78=1,nspin
if(all(xyzzyaaal78(:,xyzzyaaaa78)/=xyzzyaaaf78))cycle
do xyzzyaaab78=1,nele(xyzzyaaaf78)
xyzzyaaae78=which_ii(xyzzyaaab78,xyzzyaaaf78)
if(on_top_ii/=0.and.xyzzyaaae78/=on_top_ii)cycle
do xyzzyaaad78=1,xyzzyaaai78+xyzzyaaah78
xyzzyaaan78=rele_scr(:,xyzzyaaae78,is)+xyzzyaaar78(:,xyzzyaaad78)
if(on_top_ii==0)then
call define_config_oneelec(xyzzyaaae78,is,js,xyzzyaaan78,sele_scr(xyzz&
&yaaae78,is))
call wfn_ratio(is,js,0,ratio=xyzzyaaas78,isnan=isnan,isinf=isinf)
else
call define_config_twoelec(xyzzyaaae78,on_top_jj,is,ks,xyzzyaaan78,sel&
&e_scr(xyzzyaaae78,is),xyzzyaaan78,sele_scr(on_top_jj,is))
call wfn_ratio(is,ks,0,ratio=xyzzyaaas78,isnan=isnan,isinf=isinf)
endif
if(isnan.or.isinf)xyzzyaaas78=0.d0
if(onep_density_mat.and.xyzzyaaav78(xyzzyaaad78))then
xyzzyaaag78=xyzzyaaak78(xyzzyaaad78)
xyzzyaaao78=xyzzyaaho1(xyzzyaaaa78)*real(xyzzyaaas78,dp)
xyzzyaaap78=xyzzyaaao78*xyzzyaaao78
xyzzyaaao78=weight*xyzzyaaao78
xyzzyaaap78=weight*xyzzyaaap78
if(.not.postweight)then
xyzzyaahi1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaahi1(xyzzyaaag78,xyzzyaaaa78&
&)+xyzzyaaao78
xyzzyaahj1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaahj1(xyzzyaaag78,xyzzyaaaa78&
&)+xyzzyaaap78
xyzzyaahf1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaahf1(xyzzyaaag78,xyzzyaaaa78&
&)+weight
else
xyzzyaahp1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaahp1(xyzzyaaag78,xyzzyaaaa78&
&)+xyzzyaaao78
xyzzyaahq1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaahq1(xyzzyaaag78,xyzzyaaaa78&
&)+xyzzyaaap78
xyzzyaahr1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaahr1(xyzzyaaag78,xyzzyaaaa78&
&)+weight
endif
endif
if(mom_den.and.xyzzyaaad78<=xyzzyaaai78)then
xyzzyaaau78=xyzzyaaif1(xyzzyaaaa78)*xyzzyaaas78
if(.not.postweight)then
xyzzyaahw1(xyzzyaaaa78)=xyzzyaahw1(xyzzyaaaa78)+weight
else
xyzzyaaii1(xyzzyaaaa78)=xyzzyaaii1(xyzzyaaaa78)+weight
endif
do xyzzyaaag78=1,expval_ngvec(1)
xyzzyaaaq78=ddot(3,xyzzyaaij1(1,xyzzyaaag78),1,xyzzyaaar78(1,xyzzyaaad&
&78),1)
xyzzyaaat78=cmplx(cos(xyzzyaaaq78),sin(xyzzyaaaq78),dp)
xyzzyaaao78=dble(xyzzyaaau78*xyzzyaaat78)
xyzzyaaap78=xyzzyaaao78*xyzzyaaao78
if(.not.postweight)then
xyzzyaahz1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaahz1(xyzzyaaag78,xyzzyaaaa78&
&)+xyzzyaaao78
xyzzyaaia1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaaia1(xyzzyaaag78,xyzzyaaaa78&
&)+xyzzyaaap78
else
xyzzyaaao78=weight*xyzzyaaao78
xyzzyaaap78=weight*xyzzyaaap78
xyzzyaaig1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaaig1(xyzzyaaag78,xyzzyaaaa78&
&)+xyzzyaaao78
xyzzyaaih1(xyzzyaaag78,xyzzyaaaa78)=xyzzyaaih1(xyzzyaaag78,xyzzyaaaa78&
&)+xyzzyaaap78
endif
enddo
endif
enddo
enddo
enddo
enddo
end subroutine xyzzyaarv1
subroutine xyzzyaarw1(is,js,weight,postweight)
implicit none
integer,intent(in) :: is,js
real(dp),intent(in) :: weight
logical,intent(in) :: postweight
integer xyzzyaaaa79,xyzzyaaab79,xyzzyaaac79,xyzzyaaad79,xyzzyaaae79,xy&
&zzyaaaf79,xyzzyaaag79,xyzzyaaah79
integer,save :: xyzzyaaai79,xyzzyaaaj79
real(dp),save :: xyzzyaaak79
integer,allocatable,save :: xyzzyaaal79(:),xyzzyaaam79(:),xyzzyaaan79(&
&:,:,:)
real(dp) xyzzyaaao79,xyzzyaaap79(3),xyzzyaaaq79(3),xyzzyaaar79,xyzzyaa&
&as79,xyzzyaaat79
real(dp),allocatable,save :: xyzzyaaau79(:,:)
complex(dp) xyzzyaaav79,xyzzyaaaw79,xyzzyaaax79
logical isnan,isinf
logical,allocatable,save :: xyzzyaaay79(:)
logical,save :: xyzzyaaaz79=.true.
if(xyzzyaaaz79)then
if(twop_density_mat)then
xyzzyaaai79=xyzzyaaim1
xyzzyaaaj79=xyzzyaaik1
xyzzyaaak79=xyzzyaaiv1
else
xyzzyaaai79=xyzzyaand1
xyzzyaaaj79=xyzzyaanc1
xyzzyaaak79=xyzzyaanh1
endif
allocate(xyzzyaaau79(3,2*xyzzyaaai79),xyzzyaaal79(2*xyzzyaaai79),xyzzy&
&aaam79(xyzzyaaaj79),xyzzyaaan79(2,max_spin_pairs,xyzzyaaaj79),xyzzyaa&
&ay79(2*xyzzyaaai79),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'ACCUMULATE_COND_FRAC','0.7')
xyzzyaaam79=0
xyzzyaaan79=0
if(twop_density_mat)then
xyzzyaaam79=xyzzyaain1
xyzzyaaan79=xyzzyaaio1
else
xyzzyaaam79=xyzzyaane1
xyzzyaaan79=xyzzyaanf1
endif
xyzzyaaaz79=.false.
endif
if(.not.twop_dm_mom)then
xyzzyaaay79=.true.
xyzzyaaah79=0
do xyzzyaaab79=1,xyzzyaaai79
call xyzzyaasc1(xyzzyaaau79(1:3,xyzzyaaab79),xyzzyaaao79)
xyzzyaaal79(xyzzyaaab79)=int(xyzzyaaao79*xyzzyaaiu1)+1
enddo
elseif(.not.twop_density_mat)then
xyzzyaaay79=.false.
xyzzyaaah79=0
do xyzzyaaab79=1,xyzzyaaai79
call xyzzyaasd1(xyzzyaaau79(1:3,xyzzyaaab79),xyzzyaaao79,xyzzyaaay79(x&
&yzzyaaab79))
enddo
else
xyzzyaaay79=.false.
do xyzzyaaab79=1,xyzzyaaai79
call xyzzyaasd1(xyzzyaaau79(1:3,xyzzyaaab79),xyzzyaaao79,xyzzyaaay79(x&
&yzzyaaab79))
xyzzyaaal79(xyzzyaaab79)=int(xyzzyaaao79*xyzzyaaiu1)+1
enddo
xyzzyaaah79=count(.not.xyzzyaaay79(1:xyzzyaaai79))
do xyzzyaaab79=xyzzyaaai79+1,xyzzyaaai79+xyzzyaaah79
call xyzzyaasc1(xyzzyaaau79(1:3,xyzzyaaab79),xyzzyaaao79)
xyzzyaaal79(xyzzyaaab79)=int(xyzzyaaao79*xyzzyaaiu1)+1
enddo
endif
if(xyzzyaaak79<=.5d0)then
do xyzzyaaaa79=1,xyzzyaaaj79
do xyzzyaaab79=1,xyzzyaais1(xyzzyaaaa79)
do
xyzzyaaac79=int(ranx()*xyzzyaaip1(xyzzyaaaa79))+1
if(xyzzyaaab79>1.and.any(xyzzyaair1(1:xyzzyaaab79,xyzzyaaaa79)==xyzzya&
&aac79))cycle
xyzzyaair1(xyzzyaaab79,xyzzyaaaa79)=xyzzyaaac79
exit
enddo
enddo
enddo
elseif(xyzzyaaak79<1.d0)then
do xyzzyaaaa79=1,xyzzyaaaj79
do xyzzyaaab79=1,xyzzyaaip1(xyzzyaaaa79)-xyzzyaais1(xyzzyaaaa79)
do
xyzzyaaac79=int(ranx()*xyzzyaaip1(xyzzyaaaa79))+1
if(xyzzyaaab79>1.and.any(xyzzyaait1(1:xyzzyaaab79)==xyzzyaaac79))cycle
xyzzyaait1(xyzzyaaab79)=xyzzyaaac79
exit
enddo
enddo
xyzzyaaac79=xyzzyaaip1(xyzzyaaaa79)-xyzzyaais1(xyzzyaaaa79)
xyzzyaaad79=0
do xyzzyaaab79=1,xyzzyaaip1(xyzzyaaaa79)
if(all(xyzzyaait1(1:xyzzyaaac79)/=xyzzyaaab79))then
xyzzyaaad79=xyzzyaaad79+1
xyzzyaair1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaaab79
endif
enddo
enddo
endif
do xyzzyaaaa79=1,xyzzyaaaj79
do xyzzyaaab79=1,xyzzyaais1(xyzzyaaaa79)
xyzzyaaae79=xyzzyaair1(xyzzyaaab79,xyzzyaaaa79)
xyzzyaaaf79=xyzzyaaiq1(1,xyzzyaaae79,xyzzyaaaa79)
xyzzyaaag79=xyzzyaaiq1(2,xyzzyaaae79,xyzzyaaaa79)
do xyzzyaaac79=1,xyzzyaaai79
xyzzyaaap79=rele_scr(:,xyzzyaaaf79,is)+xyzzyaaau79(:,xyzzyaaac79)
xyzzyaaaq79=rele_scr(:,xyzzyaaag79,is)+xyzzyaaau79(:,xyzzyaaac79)
call define_config_twoelec(xyzzyaaaf79,xyzzyaaag79,is,js,xyzzyaaap79,s&
&ele_scr(xyzzyaaaf79,is),xyzzyaaaq79,sele_scr(xyzzyaaag79,is))
call wfn_ratio(is,js,0,ratio=xyzzyaaav79,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)xyzzyaaav79=0.d0
if(twop_density_mat.and.xyzzyaaay79(xyzzyaaac79))then
xyzzyaaad79=xyzzyaaal79(xyzzyaaac79)
xyzzyaaar79=xyzzyaaji1(xyzzyaaaa79)*real(xyzzyaaav79,dp)
xyzzyaaas79=xyzzyaaar79*xyzzyaaar79
xyzzyaaar79=xyzzyaaar79*weight
xyzzyaaas79=xyzzyaaas79*weight
if(.not.postweight)then
xyzzyaaiz1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaaiz1(xyzzyaaad79,xyzzyaaaa79&
&)+xyzzyaaar79
xyzzyaaja1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaaja1(xyzzyaaad79,xyzzyaaaa79&
&)+xyzzyaaas79
xyzzyaaiw1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaaiw1(xyzzyaaad79,xyzzyaaaa79&
&)+weight
else
xyzzyaajf1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaajf1(xyzzyaaad79,xyzzyaaaa79&
&)+xyzzyaaar79
xyzzyaajg1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaajg1(xyzzyaaad79,xyzzyaaaa79&
&)+xyzzyaaas79
xyzzyaajh1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaajh1(xyzzyaaad79,xyzzyaaaa79&
&)+weight
endif
endif
if(twop_dm_mom.and.xyzzyaaac79<=xyzzyaaai79)then
xyzzyaaax79=xyzzyaanv1(xyzzyaaaa79)*xyzzyaaav79
if(.not.postweight)then
xyzzyaani1(xyzzyaaaa79)=xyzzyaani1(xyzzyaaaa79)+weight
else
xyzzyaant1(xyzzyaaaa79)=xyzzyaant1(xyzzyaaaa79)+weight
endif
do xyzzyaaad79=1,expval_ngvec(1)
xyzzyaaat79=ddot(3,xyzzyaanu1(1,xyzzyaaad79),1,xyzzyaaau79(1,xyzzyaaac&
&79),1)
xyzzyaaaw79=cmplx(cos(xyzzyaaat79),sin(xyzzyaaat79),dp)
xyzzyaaar79=dble(xyzzyaaax79*xyzzyaaaw79)
xyzzyaaas79=xyzzyaaar79*xyzzyaaar79
if(.not.postweight)then
xyzzyaanl1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaanl1(xyzzyaaad79,xyzzyaaaa79&
&)+xyzzyaaar79
xyzzyaanm1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaanm1(xyzzyaaad79,xyzzyaaaa79&
&)+xyzzyaaas79
else
xyzzyaaar79=weight*xyzzyaaar79
xyzzyaaas79=weight*xyzzyaaas79
xyzzyaanr1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaanr1(xyzzyaaad79,xyzzyaaaa79&
&)+xyzzyaaar79
xyzzyaans1(xyzzyaaad79,xyzzyaaaa79)=xyzzyaans1(xyzzyaaad79,xyzzyaaaa79&
&)+xyzzyaaas79
endif
enddo
endif
enddo
enddo
enddo
end subroutine xyzzyaarw1
subroutine xyzzyaarx1(is,js,ks,weight,postweight)
implicit none
integer,intent(in) :: is,js,ks
real(dp),intent(in) :: weight
logical,intent(in) :: postweight
integer xyzzyaaaa80,xyzzyaaab80,xyzzyaaac80,xyzzyaaad80,xyzzyaaae80,xy&
&zzyaaaf80,xyzzyaaag80
integer,save :: xyzzyaaah80,xyzzyaaai80,xyzzyaaaj80
real(dp),save :: xyzzyaaak80
integer,allocatable,save :: xyzzyaaal80(:),xyzzyaaam80(:),xyzzyaaan80(&
&:,:,:)
real(dp) xyzzyaaao80,xyzzyaaap80(3),xyzzyaaaq80(3),xyzzyaaar80(3),xyzz&
&yaaas80(3),xyzzyaaat80,xyzzyaaau80,xyzzyaaav80
real(dp),allocatable,save :: xyzzyaaaw80(:,:)
logical isnan,isinf,xyzzyaaax80
logical,allocatable,save :: xyzzyaaay80(:)
logical,save :: xyzzyaaaz80=.true.
complex(dp) xyzzyaaba80,xyzzyaabb80,xyzzyaabc80,xyzzyaabd80,xyzzyaabe8&
&0,xyzzyaabf80,xyzzyaabg80
if(xyzzyaaaz80)then
if(cond_fraction)then
xyzzyaaah80=xyzzyaajl1
xyzzyaaaj80=xyzzyaajj1
xyzzyaaak80=xyzzyaaju1
else
xyzzyaaah80=xyzzyaanx1
xyzzyaaaj80=xyzzyaanw1
xyzzyaaak80=xyzzyaaob1
endif
allocate(xyzzyaaaw80(3,2*xyzzyaaah80),xyzzyaaal80(2*xyzzyaaah80),xyzzy&
&aaam80(xyzzyaaaj80),xyzzyaaan80(2,max_spin_pairs,xyzzyaaaj80),xyzzyaa&
&ay80(2*xyzzyaaah80),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'ACCUMULATE_COND_FRAC','0.7')
xyzzyaaam80=0
xyzzyaaan80=0
if(cond_fraction)then
xyzzyaaam80=xyzzyaajm1
xyzzyaaan80=xyzzyaajn1
else
xyzzyaaam80=xyzzyaany1
xyzzyaaan80=xyzzyaanz1
endif
xyzzyaaaz80=.false.
endif
if(.not.cond_frac_mom)then
xyzzyaaay80=.true.
xyzzyaaai80=0
do xyzzyaaab80=1,xyzzyaaah80
call xyzzyaasc1(xyzzyaaaw80(1:3,xyzzyaaab80),xyzzyaaao80)
xyzzyaaal80(xyzzyaaab80)=int(xyzzyaaao80*xyzzyaajt1)+1
enddo
elseif(.not.cond_fraction)then
xyzzyaaay80=.false.
xyzzyaaai80=0
do xyzzyaaab80=1,xyzzyaaah80
call xyzzyaasd1(xyzzyaaaw80(1:3,xyzzyaaab80),xyzzyaaao80,xyzzyaaay80(x&
&yzzyaaab80))
enddo
else
xyzzyaaay80=.false.
do xyzzyaaab80=1,xyzzyaaah80
call xyzzyaasd1(xyzzyaaaw80(1:3,xyzzyaaab80),xyzzyaaao80,xyzzyaaay80(x&
&yzzyaaab80))
xyzzyaaal80(xyzzyaaab80)=int(xyzzyaaao80*xyzzyaajt1)+1
enddo
xyzzyaaai80=count(.not.xyzzyaaay80(1:xyzzyaaah80))
do xyzzyaaab80=xyzzyaaah80+1,xyzzyaaah80+xyzzyaaai80
call xyzzyaasc1(xyzzyaaaw80(1:3,xyzzyaaab80),xyzzyaaao80)
xyzzyaaal80(xyzzyaaab80)=int(xyzzyaaao80*xyzzyaajt1)+1
enddo
endif
if(xyzzyaaak80<=.5d0)then
do xyzzyaaaa80=1,xyzzyaaaj80
do xyzzyaaab80=1,xyzzyaajr1(xyzzyaaaa80)
do
xyzzyaaac80=int(ranx()*xyzzyaajo1(xyzzyaaaa80))+1
if(xyzzyaaab80>1.and.any(xyzzyaajq1(1:xyzzyaaab80,xyzzyaaaa80)==xyzzya&
&aac80))cycle
xyzzyaajq1(xyzzyaaab80,xyzzyaaaa80)=xyzzyaaac80
exit
enddo
enddo
enddo
elseif(xyzzyaaak80<1.d0)then
do xyzzyaaaa80=1,xyzzyaaaj80
do xyzzyaaab80=1,xyzzyaajo1(xyzzyaaaa80)-xyzzyaajr1(xyzzyaaaa80)
do
xyzzyaaac80=int(ranx()*xyzzyaajo1(xyzzyaaaa80))+1
if(xyzzyaaab80>1.and.any(xyzzyaajs1(1:xyzzyaaab80)==xyzzyaaac80))cycle
xyzzyaajs1(xyzzyaaab80)=xyzzyaaac80
exit
enddo
enddo
xyzzyaaac80=xyzzyaajo1(xyzzyaaaa80)-xyzzyaajr1(xyzzyaaaa80)
xyzzyaaad80=0
do xyzzyaaab80=1,xyzzyaajo1(xyzzyaaaa80)
if(all(xyzzyaajs1(1:xyzzyaaac80)/=xyzzyaaab80))then
xyzzyaaad80=xyzzyaaad80+1
xyzzyaajq1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaaab80
endif
enddo
enddo
endif
do xyzzyaaaa80=1,xyzzyaaaj80
xyzzyaaax80=xyzzyaaan80(1,1,xyzzyaaaa80)==xyzzyaaan80(2,1,xyzzyaaaa80)
do xyzzyaaab80=1,xyzzyaajr1(xyzzyaaaa80)
xyzzyaaae80=xyzzyaajq1(xyzzyaaab80,xyzzyaaaa80)
xyzzyaaaf80=xyzzyaajp1(1,xyzzyaaae80,xyzzyaaaa80)
xyzzyaaag80=xyzzyaajp1(2,xyzzyaaae80,xyzzyaaaa80)
do xyzzyaaac80=1,xyzzyaaah80
xyzzyaaap80=rele_scr(:,xyzzyaaaf80,is)+xyzzyaaaw80(:,xyzzyaaac80)
xyzzyaaaq80=rele_scr(:,xyzzyaaag80,is)+xyzzyaaaw80(:,xyzzyaaac80)
call define_config_oneelec(xyzzyaaag80,is,js,xyzzyaaaq80,sele_scr(xyzz&
&yaaag80,is))
call wfn_ratio(is,js,0,ratio=xyzzyaabb80,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)xyzzyaabb80=0.d0
call define_config_oneelec(xyzzyaaaf80,is,js,xyzzyaaap80,sele_scr(xyzz&
&yaaaf80,is))
call wfn_ratio(is,js,0,ratio=xyzzyaaba80,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)xyzzyaaba80=0.d0
if(xyzzyaaax80)then
xyzzyaaar80=rele_scr(:,xyzzyaaag80,is)+xyzzyaaaw80(:,xyzzyaaac80)
xyzzyaaas80=rele_scr(:,xyzzyaaaf80,is)+xyzzyaaaw80(:,xyzzyaaac80)
call define_config_oneelec(xyzzyaaag80,is,js,xyzzyaaas80,sele_scr(xyzz&
&yaaag80,is))
call wfn_ratio(is,js,0,ratio=xyzzyaabd80,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)xyzzyaabd80=0.d0
call define_config_oneelec(xyzzyaaaf80,is,js,xyzzyaaar80,sele_scr(xyzz&
&yaaaf80,is))
call wfn_ratio(is,js,0,ratio=xyzzyaabc80,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)xyzzyaabc80=0.d0
endif
call define_config_twoelec(xyzzyaaaf80,xyzzyaaag80,is,ks,xyzzyaaap80,s&
&ele_scr(xyzzyaaaf80,is),xyzzyaaaq80,sele_scr(xyzzyaaag80,is))
call wfn_ratio(is,ks,0,ratio=xyzzyaabe80,isnan=isnan,isinf=isinf)
if(isnan.or.isinf)xyzzyaabe80=0.d0
xyzzyaabe80=xyzzyaabe80-xyzzyaaba80*xyzzyaabb80
if(xyzzyaaax80)xyzzyaabe80=xyzzyaabe80+xyzzyaabc80*xyzzyaabd80
if(cond_fraction.and.xyzzyaaay80(xyzzyaaac80))then
xyzzyaaad80=xyzzyaaal80(xyzzyaaac80)
xyzzyaaat80=xyzzyaakh1(xyzzyaaaa80)*real(xyzzyaabe80,dp)
xyzzyaaau80=xyzzyaaat80*xyzzyaaat80
xyzzyaaat80=xyzzyaaat80*weight
xyzzyaaau80=xyzzyaaau80*weight
if(.not.postweight)then
xyzzyaajy1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaajy1(xyzzyaaad80,xyzzyaaaa80&
&)+xyzzyaaat80
xyzzyaajz1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaajz1(xyzzyaaad80,xyzzyaaaa80&
&)+xyzzyaaau80
xyzzyaajv1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaajv1(xyzzyaaad80,xyzzyaaaa80&
&)+weight
else
xyzzyaake1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaake1(xyzzyaaad80,xyzzyaaaa80&
&)+xyzzyaaat80
xyzzyaakf1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaakf1(xyzzyaaad80,xyzzyaaaa80&
&)+xyzzyaaau80
xyzzyaakg1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaakg1(xyzzyaaad80,xyzzyaaaa80&
&)+weight
endif
endif
if(cond_frac_mom.and.xyzzyaaac80<=xyzzyaaah80)then
xyzzyaabg80=xyzzyaaop1(xyzzyaaaa80)*xyzzyaabe80
if(.not.postweight)then
xyzzyaaoc1(xyzzyaaaa80)=xyzzyaaoc1(xyzzyaaaa80)+weight
else
xyzzyaaon1(xyzzyaaaa80)=xyzzyaaon1(xyzzyaaaa80)+weight
endif
do xyzzyaaad80=1,expval_ngvec(1)
xyzzyaaav80=ddot(3,xyzzyaaoo1(1,xyzzyaaad80),1,xyzzyaaaw80(1,xyzzyaaac&
&80),1)
xyzzyaabf80=cmplx(cos(xyzzyaaav80),sin(xyzzyaaav80),dp)
xyzzyaaat80=dble(xyzzyaabg80*xyzzyaabf80)
xyzzyaaau80=xyzzyaaat80*xyzzyaaat80
if(.not.postweight)then
xyzzyaaof1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaaof1(xyzzyaaad80,xyzzyaaaa80&
&)+xyzzyaaat80
xyzzyaaog1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaaog1(xyzzyaaad80,xyzzyaaaa80&
&)+xyzzyaaau80
else
xyzzyaaat80=weight*xyzzyaaat80
xyzzyaaau80=weight*xyzzyaaau80
xyzzyaaol1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaaol1(xyzzyaaad80,xyzzyaaaa80&
&)+xyzzyaaat80
xyzzyaaom1(xyzzyaaad80,xyzzyaaaa80)=xyzzyaaom1(xyzzyaaad80,xyzzyaaaa80&
&)+xyzzyaaau80
endif
enddo
endif
enddo
enddo
enddo
end subroutine xyzzyaarx1
subroutine xyzzyaary1(is,weight,postweight)
implicit none
integer,intent(in) :: is
logical,intent(in) :: postweight
real(dp),intent(in) :: weight
integer xyzzyaaaa81,xyzzyaaab81,xyzzyaaac81,xyzzyaaad81,xyzzyaaae81,xy&
&zzyaaaf81,xyzzyaaag81,xyzzyaaah81,xyzzyaaai81
real(dp) xyzzyaaaj81,xyzzyaaak81,xyzzyaaal81,xyzzyaaam81,xyzzyaaan81,x&
&yzzyaaao81,xyzzyaaap81,xyzzyaaaq81,xyzzyaaar81,xyzzyaaas81
call get_eivecs(is)
xyzzyaaal81=xyzzyaakq1
xyzzyaaaq81=1.6d0
xyzzyaaam81=(dble(netot)-5.d0/16.d0)*xyzzyaaaq81**(-10)
if(xyzzyaakl1==xyzzyaakm1)then
xyzzyaaas81=dble(xyzzyaakk1)/(2.d0*xyzzyaaal81)
xyzzyaaar81=dble(xyzzyaakk1)/xyzzyaaal81
endif
if(.not.postweight)then
xyzzyaaks1=xyzzyaaks1+1.d0
xyzzyaakw1=xyzzyaakw1+weight
xyzzyaala1=xyzzyaala1+weight*weight
else
xyzzyaakv1=xyzzyaakv1+1.d0
xyzzyaakz1=xyzzyaakz1+weight
xyzzyaald1=xyzzyaald1+weight*weight
endif
call get_eevecs(is)
do xyzzyaaae81=1,nspin
do xyzzyaaaf81=xyzzyaaae81,nspin
if((xyzzyaaae81==1).and.(xyzzyaaaf81==1))xyzzyaaai81=1
if((xyzzyaaae81==2).and.(xyzzyaaaf81==2))xyzzyaaai81=3
if(xyzzyaaae81/=xyzzyaaaf81)xyzzyaaai81=2
do xyzzyaaaa81=1,nele(xyzzyaaae81)
xyzzyaaag81=1
if(xyzzyaaae81==xyzzyaaaf81)xyzzyaaag81=xyzzyaaaa81+1
do xyzzyaaab81=xyzzyaaag81,nele(xyzzyaaaf81)
xyzzyaaac81=which_ii(xyzzyaaaa81,xyzzyaaae81)
xyzzyaaad81=which_ii(xyzzyaaab81,xyzzyaaaf81)
xyzzyaaak81=eevecs_scr(4,xyzzyaaac81,xyzzyaaad81,is)
if(xyzzyaakl1/=xyzzyaakp1.and.xyzzyaaak81>=2.d0*xyzzyaaal81)cycle
if(xyzzyaakl1==xyzzyaakm1)then
xyzzyaaah81=int(xyzzyaaak81*xyzzyaaas81)+1
xyzzyaaan81=1.d0
xyzzyaaao81=xyzzyaaan81*xyzzyaaan81
xyzzyaaan81=weight*xyzzyaaan81
xyzzyaaao81=weight*xyzzyaaao81
xyzzyaale1(xyzzyaaah81,xyzzyaaai81)=xyzzyaale1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaan81
xyzzyaalf1(xyzzyaaah81,xyzzyaaai81)=xyzzyaalf1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaao81
else
call xyzzyaasf1(xyzzyaaak81,xyzzyaakk1,2.d0*xyzzyaaal81,xyzzyaaam81,xy&
&zzyaaaq81,xyzzyaaap81,xyzzyaalu1)
do xyzzyaaah81=1,xyzzyaakk1
xyzzyaaan81=xyzzyaalu1(xyzzyaaah81)*xyzzyaaap81
xyzzyaaao81=xyzzyaaan81*xyzzyaaan81
xyzzyaaan81=weight*xyzzyaaan81
xyzzyaaao81=weight*xyzzyaaao81
if(.not.postweight)then
xyzzyaale1(xyzzyaaah81,xyzzyaaai81)=xyzzyaale1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaan81
xyzzyaalf1(xyzzyaaah81,xyzzyaaai81)=xyzzyaalf1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaao81
else
xyzzyaalk1(xyzzyaaah81,xyzzyaaai81)=xyzzyaalk1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaan81
xyzzyaall1(xyzzyaaah81,xyzzyaaai81)=xyzzyaall1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaao81
endif
enddo
endif
enddo
enddo
enddo
enddo
call get_eivecs(is)
do xyzzyaaae81=1,nspin
xyzzyaaai81=xyzzyaaae81
do xyzzyaaaa81=1,nele(xyzzyaaae81)
xyzzyaaac81=which_ii(xyzzyaaaa81,xyzzyaaae81)
xyzzyaaaj81=eivecs_scr(4,1,xyzzyaaac81,is)
if(xyzzyaakl1/=xyzzyaakp1.and.xyzzyaaaj81>=xyzzyaaal81)cycle
if(xyzzyaakl1==xyzzyaakm1)then
xyzzyaaah81=int(xyzzyaaaj81*xyzzyaaar81)+1
xyzzyaaan81=1.d0
xyzzyaaao81=xyzzyaaan81*xyzzyaaan81
xyzzyaaan81=weight*xyzzyaaan81
xyzzyaaao81=weight*xyzzyaaao81
xyzzyaalm1(xyzzyaaah81,xyzzyaaai81)=xyzzyaalm1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaan81
xyzzyaaln1(xyzzyaaah81,xyzzyaaai81)=xyzzyaaln1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaao81
else
call xyzzyaasf1(xyzzyaaaj81,xyzzyaakk1,xyzzyaaal81,xyzzyaaam81,xyzzyaa&
&aq81,xyzzyaaap81,xyzzyaalu1)
do xyzzyaaah81=1,xyzzyaakk1
xyzzyaaan81=xyzzyaalu1(xyzzyaaah81)*xyzzyaaap81
xyzzyaaao81=xyzzyaaan81*xyzzyaaan81
xyzzyaaan81=weight*xyzzyaaan81
xyzzyaaao81=weight*xyzzyaaao81
if(.not.postweight)then
xyzzyaalm1(xyzzyaaah81,xyzzyaaai81)=xyzzyaalm1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaan81
xyzzyaaln1(xyzzyaaah81,xyzzyaaai81)=xyzzyaaln1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaao81
else
xyzzyaals1(xyzzyaaah81,xyzzyaaai81)=xyzzyaals1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaan81
xyzzyaalt1(xyzzyaaah81,xyzzyaaai81)=xyzzyaalt1(xyzzyaaah81,xyzzyaaai81&
&)+xyzzyaaao81
endif
enddo
endif
enddo
enddo
end subroutine xyzzyaary1
subroutine xyzzyaarz1(is,weight,postweight)
implicit none
integer,intent(in) :: is
real(dp),intent(in) :: weight
logical,intent(in) :: postweight
integer xyzzyaaaa82,xyzzyaaab82,xyzzyaaac82,xyzzyaaad82,xyzzyaaae82,xy&
&zzyaaaf82
real(dp) xyzzyaaag82,xyzzyaaah82(nitot)
if(nitot==0)return
do xyzzyaaaa82=1,xyzzyaamn1
xyzzyaaah82=0.d0
do xyzzyaaad82=1,nspin
if(all(xyzzyaamp1(:,xyzzyaaaa82)/=xyzzyaaad82))cycle
do xyzzyaaab82=1,nele(xyzzyaaad82)
xyzzyaaac82=which_ii(xyzzyaaab82,xyzzyaaad82)
xyzzyaaaf82=1
xyzzyaaag82=eivecs_scr(4,1,xyzzyaaac82,is)
do xyzzyaaae82=2,nitot
if(eivecs_scr(4,xyzzyaaae82,xyzzyaaac82,is)<xyzzyaaag82)then
xyzzyaaaf82=xyzzyaaae82
xyzzyaaag82=eivecs_scr(4,xyzzyaaae82,xyzzyaaac82,is)
endif
enddo
xyzzyaaah82(xyzzyaaaf82)=xyzzyaaah82(xyzzyaaaf82)+1.d0
enddo
enddo
if(.not.postweight)then
call daxpy(nitot,weight,xyzzyaaah82(1),1,xyzzyaamt1(1,xyzzyaaaa82),1)
xyzzyaamu1(:,xyzzyaaaa82)=xyzzyaamu1(:,xyzzyaaaa82)+xyzzyaaah82(:)**2*&
&weight
else
call daxpy(nitot,weight,xyzzyaaah82(1),1,xyzzyaamz1(1,xyzzyaaaa82),1)
xyzzyaana1(:,xyzzyaaaa82)=xyzzyaana1(:,xyzzyaaaa82) +xyzzyaaah82(:)**2&
&*weight
endif
enddo
if(.not.postweight)then
xyzzyaamq1(:)=xyzzyaamq1(:)+weight
else
xyzzyaanb1(:)=xyzzyaanb1(:)+weight
endif
end subroutine xyzzyaarz1
subroutine xyzzyaasa1(is,weight,postweight)
implicit none
integer,intent(in) :: is
logical,intent(in) :: postweight
real(dp),intent(in) :: weight
integer xyzzyaaaa83,xyzzyaaab83,xyzzyaaac83,xyzzyaaad83,xyzzyaaae83(3)
real(dp) xyzzyaaaf83(3),xyzzyaaag83(3),xyzzyaaah83(3,netot)
call get_rsele(is)
xyzzyaaah83=rele_scr(:,:,is)
xyzzyaaag83=dble(xyzzyaalw1)/(xyzzyaaly1-xyzzyaalx1)
if(.not.postweight)then
xyzzyaalz1=xyzzyaalz1+1.d0
xyzzyaamd1=xyzzyaamd1+weight
else
xyzzyaamc1=xyzzyaamc1+1.d0
xyzzyaamg1=xyzzyaamg1+weight
endif
call get_eivecs(is)
do xyzzyaaac83=1,nspin
if(mol_spin_density)then
xyzzyaaad83=xyzzyaaac83
else
xyzzyaaad83=1
endif
do xyzzyaaaa83=1,nele(xyzzyaaac83)
xyzzyaaab83=which_ii(xyzzyaaaa83,xyzzyaaac83)
xyzzyaaaf83=xyzzyaaah83(1:3,xyzzyaaab83)
xyzzyaaae83=1+floor((xyzzyaaaf83-xyzzyaalx1)*xyzzyaaag83)
if(any(xyzzyaaae83<1).or.any(xyzzyaaae83>xyzzyaalw1))cycle
xyzzyaamh1(xyzzyaaae83(1),xyzzyaaae83(2),xyzzyaaae83(3),xyzzyaaad83)=x&
&yzzyaamh1(xyzzyaaae83(1),xyzzyaaae83(2),xyzzyaaae83(3),xyzzyaaad83)+w&
&eight
enddo
enddo
end subroutine xyzzyaasa1
integer function xyzzyaasb1(pointrange,cutoff)
implicit none
integer,intent(in) :: pointrange(3)
real(dp),intent(in) :: cutoff
integer xyzzyaaaa84,xyzzyaaab84,xyzzyaaac84,xyzzyaaad84,xyzzyaaae84,xy&
&zzyaaaf84,xyzzyaaag84,xyzzyaaah84,xyzzyaaai84,xyzzyaaaj84
real(dp) xyzzyaaak84,xyzzyaaal84,xyzzyaaam84,xyzzyaaan84,xyzzyaaao84,x&
&yzzyaaap84,xyzzyaaaq84
xyzzyaaae84=-pointrange(1)+1
xyzzyaaaf84=pointrange(1)
if(dimensionality==3)then
xyzzyaaag84=-pointrange(2)+1
xyzzyaaah84=pointrange(2)
xyzzyaaai84=-pointrange(3)+1
xyzzyaaaj84=pointrange(3)
elseif(dimensionality==2)then
xyzzyaaag84=-pointrange(2)+1
xyzzyaaah84=pointrange(2)
xyzzyaaai84=0
xyzzyaaaj84=0
else
xyzzyaaag84=0
xyzzyaaah84=0
xyzzyaaai84=0
xyzzyaaaj84=0
endif
xyzzyaaad84=0
do xyzzyaaac84=xyzzyaaai84,xyzzyaaaj84
xyzzyaaam84=real(xyzzyaaac84,dp)
do xyzzyaaab84=xyzzyaaag84,xyzzyaaah84
xyzzyaaal84=real(xyzzyaaab84,dp)
do xyzzyaaaa84=xyzzyaaae84,xyzzyaaaf84
xyzzyaaak84=real(xyzzyaaaa84,dp)
xyzzyaaan84=xyzzyaaak84*xyzzyaaap1(1)+xyzzyaaal84*xyzzyaaaq1(1)+xyzzya&
&aam84*xyzzyaaar1(1)
xyzzyaaao84=xyzzyaaak84*xyzzyaaap1(2)+xyzzyaaal84*xyzzyaaaq1(2)+xyzzya&
&aam84*xyzzyaaar1(2)
xyzzyaaap84=xyzzyaaak84*xyzzyaaap1(3)+xyzzyaaal84*xyzzyaaaq1(3)+xyzzya&
&aam84*xyzzyaaar1(3)
xyzzyaaaq84=xyzzyaaan84**2+xyzzyaaao84**2+xyzzyaaap84**2
if(xyzzyaaaq84<=cutoff)xyzzyaaad84=xyzzyaaad84+1
enddo
enddo
enddo
xyzzyaasb1=xyzzyaaad84
end function xyzzyaasb1
subroutine xyzzyaasc1(rvec,r)
implicit none
real(dp),intent(out) :: rvec(3),r
real(dp) xyzzyaaaa85,xyzzyaaab85,xyzzyaaac85,xyzzyaaad85,xyzzyaaae85,x&
&yzzyaaaf85,xyzzyaaag85
if(dimensionality==3)then
r=wigner_seitz_radius*ranx()**third
xyzzyaaac85=1.d0-2.d0*ranx()
xyzzyaaaa85=twopi*ranx()
xyzzyaaab85=sqrt(1.d0-xyzzyaaac85*xyzzyaaac85)
xyzzyaaae85=cos(xyzzyaaaa85)
xyzzyaaad85=sqrt(1.d0-xyzzyaaae85*xyzzyaaae85)
if(xyzzyaaaa85>pi)xyzzyaaad85=-xyzzyaaad85
xyzzyaaaf85=r*xyzzyaaab85
rvec(1:3)=(/xyzzyaaaf85*xyzzyaaae85,xyzzyaaaf85*xyzzyaaad85,r*xyzzyaaa&
&c85/)
elseif(dimensionality==2)then
r=wigner_seitz_radius*sqrt(ranx())
xyzzyaaaa85=twopi*ranx()
xyzzyaaae85=cos(xyzzyaaaa85)
xyzzyaaad85=sqrt(1.d0-xyzzyaaae85*xyzzyaaae85)
if(xyzzyaaaa85>pi)xyzzyaaad85=-xyzzyaaad85
rvec(1:3)=(/r*xyzzyaaae85,r*xyzzyaaad85,0.d0/)
else
xyzzyaaag85=a1(1)*(ranx()-0.5d0)
r=abs(xyzzyaaag85)
rvec(1:3)=(/xyzzyaaag85,0.d0,0.d0/)
endif
end subroutine xyzzyaasc1
subroutine xyzzyaasd1(rvec,r,inside_wsr)
implicit none
real(dp),intent(out) :: rvec(3),r
logical,intent(out) :: inside_wsr
if(dimensionality==3)then
rvec(1:3)=(ranx()-0.5d0)*a1+(ranx()-0.5d0)*a2+(ranx()-0.5d0)*a3
r=sqrt(rvec(1)*rvec(1)+rvec(2)*rvec(2)+rvec(3)*rvec(3))
inside_wsr=r<=wigner_seitz_radius
elseif(dimensionality==2)then
rvec(1:2)=(ranx()-0.5d0)*a1(1:2)+(ranx()-0.5d0)*a2(1:2)
rvec(3)=0.d0
r=sqrt(rvec(1)*rvec(1)+rvec(2)*rvec(2))
inside_wsr=r<=wigner_seitz_radius
else
rvec(1)=a1(1)*(ranx()-0.5d0)
rvec(2:3)=0.d0
r=abs(rvec(1))
inside_wsr=.true.
endif
end subroutine xyzzyaasd1
subroutine switch_spin_density_to_pcf
implicit none
pair_corr=.true.
spin_density=.false.
particle_is_fixed=.true.
xyzzyaaax1(4)=.true.
if(am_master)xyzzyaadu1(:,:)=xyzzyaadv1(:,:)
end subroutine switch_spin_density_to_pcf
subroutine switch_spin_density_to_pcf_sph
implicit none
pair_corr_sph=.true.
spin_density=.false.
particle_is_fixed=.true.
xyzzyaaax1(5)=.true.
if(am_master)xyzzyaadu1(:,:)=xyzzyaadv1(:,:)
end subroutine switch_spin_density_to_pcf_sph
subroutine setup_emt
use slaarnaca,only : zion
implicit none
integer xyzzyaaaa89
xyzzyaamm1=0.d0
do xyzzyaaaa89=1,nitot
xyzzyaamm1=xyzzyaamm1+zion(iontype(xyzzyaaaa89))*rion(1:3,xyzzyaaaa89)
enddo
end subroutine setup_emt
subroutine dipole_calc(rele,dipole_moment)
implicit none
real(dp),intent(in) :: rele(3,netot)
real(dp),intent(out) :: dipole_moment(3)
integer xyzzyaaaa90
dipole_moment=xyzzyaamm1
do xyzzyaaaa90=1,netot
dipole_moment=dipole_moment-rele(1:3,xyzzyaaaa90)
enddo
end subroutine dipole_calc
subroutine contact_den_calc(dvel,relkei,eevecs,contact_den)
implicit none
real(dp),intent(in) :: dvel(3,netot,real1_complex2),relkei(netot),eeve&
&cs(4,netot,netot)
real(dp),intent(out) :: contact_den
integer xyzzyaaaa91,xyzzyaaab91,xyzzyaaac91,xyzzyaaad91,xyzzyaaae91,xy&
&zzyaaaf91
real(dp) xyzzyaaag91,xyzzyaaah91
xyzzyaaag91=0.d0
xyzzyaaaf91=3
do xyzzyaaab91=1,nele(xyzzyaaaf91)
xyzzyaaad91=which_ii(xyzzyaaab91,xyzzyaaaf91)
do xyzzyaaae91=1,2
do xyzzyaaaa91=1,nele(xyzzyaaae91)
xyzzyaaac91=which_ii(xyzzyaaaa91,xyzzyaaae91)
xyzzyaaag91=xyzzyaaag91+1.d0/eevecs(4,xyzzyaaac91,xyzzyaaad91)
enddo
enddo
enddo
xyzzyaaah91=0.d0
do xyzzyaaac91=1,netot
xyzzyaaah91=xyzzyaaah91+(relkei(xyzzyaaac91)*pmass(which_spin(xyzzyaaa&
&c91))-0.5d0*dot_product(dvel(1:3,xyzzyaaac91,1),dvel(1:3,xyzzyaaac91,&
&1)))
enddo
contact_den=one_over_twopi*xyzzyaaag91*xyzzyaaah91
end subroutine contact_den_calc
subroutine eval_int_sf(hartree,xc)
implicit none
real(dp),intent(out) :: hartree,xc
integer xyzzyaaaa92,xyzzyaaab92,xyzzyaaac92,xyzzyaaad92,xyzzyaaae92
real(dp) xyzzyaaaf92,xyzzyaaag92,xyzzyaaah92,xyzzyaaai92,xyzzyaaaj92
complex(dp) xyzzyaaak92
if(.not.structure_factor)then
xyzzyaagb1(:,:)=xyzzyaafz1(:,:)
xyzzyaagi1(:,:)=xyzzyaagh1(:,:)
endif
allocate(xyzzyaaot1(expval_ngvec(xyzzyaafr1)),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'EVAL_INT_SF','sf_sum')
xyzzyaaaf92=real(npcells**2,dp)/real(netot,dp)
xyzzyaaot1=0.d0
if(homogeneous_density)then
do xyzzyaaae92=1,xyzzyaafq1
xyzzyaaab92=xyzzyaaft1(1,xyzzyaaae92)
xyzzyaaac92=xyzzyaaft1(2,xyzzyaaae92)
if(xyzzyaaab92/=xyzzyaaac92)then
xyzzyaaag92=2.d0
else
xyzzyaaag92=1.d0
endif
xyzzyaaot1(1)=xyzzyaaot1(1)+xyzzyaaag92*(xyzzyaagb1(1,xyzzyaaae92)-nel&
&e(xyzzyaaab92)*nele(xyzzyaaac92))
do xyzzyaaaa92=2,expval_ngvec(xyzzyaafr1)
xyzzyaaot1(xyzzyaaaa92)=xyzzyaaot1(xyzzyaaaa92)+xyzzyaaag92*xyzzyaagb1&
&(xyzzyaaaa92,xyzzyaaae92)
enddo
enddo
xyzzyaaot1(:)=xyzzyaaot1(:)*xyzzyaaaf92
else
do xyzzyaaae92=1,xyzzyaafq1
xyzzyaaab92=xyzzyaaft1(1,xyzzyaaae92)
xyzzyaaac92=xyzzyaaft1(2,xyzzyaaae92)
if(xyzzyaaab92/=xyzzyaaac92)then
xyzzyaaag92=2.d0
else
xyzzyaaag92=1.d0
endif
xyzzyaaot1(1)=xyzzyaaot1(1)+xyzzyaaag92*(xyzzyaagb1(1,xyzzyaaae92)-abs&
&(xyzzyaagi1(1,xyzzyaaab92)*xyzzyaagi1(1,xyzzyaaac92)))
xyzzyaaaa92=1
do xyzzyaaad92=2,expval_ngvec(xyzzyaafr1),2
xyzzyaaaa92=xyzzyaaaa92+1
xyzzyaaot1(xyzzyaaaa92)=xyzzyaaot1(xyzzyaaaa92)+xyzzyaaag92*(xyzzyaagb&
&1(xyzzyaaaa92,xyzzyaaae92)-abs(xyzzyaagi1(xyzzyaaaa92,xyzzyaaab92)*xy&
&zzyaagi1(xyzzyaaaa92+1,xyzzyaaac92)))
xyzzyaaaa92=xyzzyaaaa92+1
xyzzyaaot1(xyzzyaaaa92)=xyzzyaaot1(xyzzyaaaa92)+xyzzyaaag92*(xyzzyaagb&
&1(xyzzyaaaa92,xyzzyaaae92)-abs(xyzzyaagi1(xyzzyaaaa92,xyzzyaaab92)*xy&
&zzyaagi1(xyzzyaaaa92-1,xyzzyaaac92)))
enddo
enddo
xyzzyaaot1(:)=xyzzyaaot1(:)*xyzzyaaaf92
endif
if(homogeneous_density)then
xyzzyaaai92=0.d0
xyzzyaaaj92=0.d0
do xyzzyaaaa92=2,expval_ngvec(xyzzyaafr1),2
xyzzyaaaf92=expval_gvec(1,xyzzyaaaa92,xyzzyaafr1)**2+expval_gvec(2,xyz&
&zyaaaa92,xyzzyaafr1)**2+expval_gvec(3,xyzzyaaaa92,xyzzyaafr1)**2
xyzzyaaah92=1.d0/xyzzyaaaf92
xyzzyaaai92=xyzzyaaai92+xyzzyaaah92*(xyzzyaaot1(xyzzyaaaa92)+xyzzyaaot&
&1(xyzzyaaaa92+1)-2.d0)
enddo
else
xyzzyaaai92=0.d0
xyzzyaaaj92=0.d0
do xyzzyaaaa92=2,expval_ngvec(xyzzyaafr1),2
xyzzyaaaf92=expval_gvec(1,xyzzyaaaa92,xyzzyaafr1)**2+expval_gvec(2,xyz&
&zyaaaa92,xyzzyaafr1)**2+expval_gvec(3,xyzzyaaaa92,xyzzyaafr1)**2
xyzzyaaah92=1.d0/xyzzyaaaf92
xyzzyaaai92=xyzzyaaai92+xyzzyaaah92*(xyzzyaaot1(xyzzyaaaa92)+xyzzyaaot&
&1(xyzzyaaaa92+1)-2.d0)
xyzzyaaak92=sum(xyzzyaagi1(xyzzyaaaa92,1:xyzzyaafs1))
xyzzyaaak92=xyzzyaaak92*sum(xyzzyaagi1(xyzzyaaaa92+1,1:xyzzyaafs1))
xyzzyaaaj92=xyzzyaaaj92+xyzzyaaah92*abs(xyzzyaaak92)
enddo
endif
xc=(real(netot,dp)/real(npcells,dp)**2)*twopi*xyzzyaaai92/pvolume+ (re&
&al(netot,dp)/real(npcells,dp))*self_term
hartree=fourpi*xyzzyaaaj92/pvolume
deallocate(xyzzyaaot1)
end subroutine eval_int_sf
subroutine finite_size_corr_xc
implicit none
integer xyzzyaaaa93,xyzzyaaab93,xyzzyaaac93,xyzzyaaad93,xyzzyaaae93,xy&
&zzyaaaf93,xyzzyaaag93
integer,allocatable :: xyzzyaaah93(:)
real(dp) xyzzyaaai93,xyzzyaaaj93,xyzzyaaak93,xyzzyaaal93,xyzzyaaam93,x&
&yzzyaaan93
real(dp),external :: calcf_xc_corr,calcf_xc_corr2
logical,parameter :: xyzzyaaao93=.true.
if(periodicity/=2.and.periodicity/=3)call errstop('FINITE_SIZE_CORR_XC&
&', 'Exercise: derive the theory of finite-size corrections in 1D.')
allocate(xyzzyaaou1(expval_ngvec(xyzzyaafr1)),xyzzyaaov1(expval_ngvec(&
&xyzzyaafr1)),xyzzyaaah93(expval_ngvec(xyzzyaafr1)+1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'FINITE_SIZE_CORR_XC','star information.')
xyzzyaaaf93=1
xyzzyaaoq1=1
xyzzyaaou1(1)=0.d0
xyzzyaaov1(1)=0.d0
xyzzyaaah93(1)=1
do xyzzyaaaa93=2,expval_ngvec(xyzzyaafr1)
xyzzyaaai93=expval_gvec(1,xyzzyaaaa93,xyzzyaafr1)**2+expval_gvec(2,xyz&
&zyaaaa93,xyzzyaafr1)**2+expval_gvec(3,xyzzyaaaa93,xyzzyaafr1)**2
if(xyzzyaaai93-xyzzyaaov1(xyzzyaaoq1)>1.d-8)then
xyzzyaaoq1=xyzzyaaoq1+1
xyzzyaaah93(xyzzyaaoq1)=xyzzyaaaa93
xyzzyaaov1(xyzzyaaoq1)=xyzzyaaai93
xyzzyaaou1(xyzzyaaoq1)=sqrt(xyzzyaaai93)
xyzzyaaaf93=0
endif
xyzzyaaaf93=xyzzyaaaf93+1
enddo
xyzzyaaah93(xyzzyaaoq1+1)=xyzzyaaah93(xyzzyaaoq1)+xyzzyaaaf93
allocate(xyzzyaaot1(xyzzyaaoq1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'FINITE_SIZE_CORR_XC','sf_sum')
xyzzyaaot1(1:xyzzyaaoq1)=0.d0
if(xyzzyaafs1==0)then
do xyzzyaaag93=1,xyzzyaafq1
xyzzyaaac93=xyzzyaaft1(1,xyzzyaaag93)
xyzzyaaad93=xyzzyaaft1(2,xyzzyaaag93)
if(xyzzyaaac93/=xyzzyaaad93)then
xyzzyaaal93=2.d0
else
xyzzyaaal93=1.d0
endif
xyzzyaaot1(1)=xyzzyaaot1(1)+xyzzyaaal93*(xyzzyaagb1(1,xyzzyaaag93)-xyz&
&zyaacn1(xyzzyaaac93)*xyzzyaacn1(xyzzyaaad93))
xyzzyaaaa93=1
do xyzzyaaae93=2,xyzzyaaoq1
xyzzyaaaf93=0
do xyzzyaaab93=xyzzyaaah93(xyzzyaaae93),xyzzyaaah93(xyzzyaaae93+1)-1,2
xyzzyaaaa93=xyzzyaaaa93+1
xyzzyaaot1(xyzzyaaae93)=xyzzyaaot1(xyzzyaaae93)+xyzzyaaal93*(xyzzyaagb&
&1(xyzzyaaaa93,xyzzyaaag93)+xyzzyaagb1(xyzzyaaaa93+1,xyzzyaaag93))
xyzzyaaaa93=xyzzyaaaa93+1
xyzzyaaaf93=xyzzyaaaf93+2
enddo
enddo
enddo
else
do xyzzyaaag93=1,xyzzyaafq1
xyzzyaaac93=xyzzyaaft1(1,xyzzyaaag93)
xyzzyaaad93=xyzzyaaft1(2,xyzzyaaag93)
if(xyzzyaaac93/=xyzzyaaad93)then
xyzzyaaal93=2.d0
else
xyzzyaaal93=1.d0
endif
xyzzyaaot1(1)=xyzzyaaot1(1)+xyzzyaaal93*(xyzzyaagb1(1,xyzzyaaag93)-rea&
&l(xyzzyaagi1(1,xyzzyaaac93)*xyzzyaagi1(1,xyzzyaaad93),dp))
xyzzyaaaa93=1
do xyzzyaaae93=2,xyzzyaaoq1
xyzzyaaaf93=0
do xyzzyaaab93=xyzzyaaah93(xyzzyaaae93),xyzzyaaah93(xyzzyaaae93+1)-1,2
xyzzyaaaa93=xyzzyaaaa93+1
xyzzyaaot1(xyzzyaaae93)=xyzzyaaot1(xyzzyaaae93)+xyzzyaaal93*(xyzzyaagb&
&1(xyzzyaaaa93,xyzzyaaag93)+xyzzyaagb1(xyzzyaaaa93+1,xyzzyaaag93)-2.d0&
&*real(xyzzyaagi1(xyzzyaaaa93,xyzzyaaac93)*xyzzyaagi1(xyzzyaaaa93+1,xy&
&zzyaaad93),dp))
xyzzyaaaa93=xyzzyaaaa93+1
xyzzyaaaf93=xyzzyaaaf93+2
enddo
enddo
enddo
endif
xyzzyaaal93=real(npcells**2,dp)/real(netot,dp)
do xyzzyaaae93=2,xyzzyaaoq1
xyzzyaaaf93=xyzzyaaah93(xyzzyaaae93+1)-xyzzyaaah93(xyzzyaaae93)
if(xyzzyaaaf93/=0)then
xyzzyaaot1(xyzzyaaae93)=xyzzyaaot1(xyzzyaaae93)*xyzzyaaal93/real(xyzzy&
&aaaf93,dp)
else
call errstop('FINITE_SIZE_CORR_XC','Zero G vectors in star - should no&
&t happen.')
endif
enddo
select case(xc_corr_method)
case(1)
if(periodicity==3)then
xyzzyaaak93=xyzzyaaot1(2)/xyzzyaaou1(2)**2
xc_corr=twopi*xyzzyaaak93*real(netot,dp)/volume
else
xyzzyaaak93=xyzzyaaot1(2)/xyzzyaaou1(2)**1.5d0
xc_corr=finite_size_const_c*dble(netot)*xyzzyaaak93/area**1.25d0
endif
chi_squared_sf=0.d0
case(2)
if(periodicity/=3)call errstop('FINITE_SIZE_CORR_XC', 'Need periodicit&
&y=3 if XC_CORR_METHOD=2.')
call xyzzyaase1
if(xyzzyaaao93)then
call gauss_quadrature(0.d0,xyzzyaaou1(xyzzyaaor1),100,calcf_xc_corr,xy&
&zzyaaaj93)
else
call gauss_quadrature(0.d0,xyzzyaaou1(xyzzyaaor1),100,calcf_xc_corr2,x&
&yzzyaaaj93)
endif
xyzzyaaam93=xyzzyaaaj93*one_over_pi
if(xyzzyaaao93)then
xyzzyaaal93=0.d0
do xyzzyaaab93=2,xyzzyaaor1
xyzzyaaaf93=xyzzyaaah93(xyzzyaaab93+1)-xyzzyaaah93(xyzzyaaab93)
xyzzyaaak93=1.d0/xyzzyaaov1(xyzzyaaab93)
xyzzyaaal93=xyzzyaaal93+real(xyzzyaaaf93,dp)*xyzzyaaak93*xyzzyaaoy1(xy&
&zzyaaab93)
enddo
xyzzyaaan93=xyzzyaaal93*twopi/volume
else
xyzzyaaal93=0.d0
do xyzzyaaab93=2,xyzzyaaor1
xyzzyaaaf93=xyzzyaaah93(xyzzyaaab93+1)-xyzzyaaah93(xyzzyaaab93)
xyzzyaaak93=1.d0/xyzzyaaov1(xyzzyaaab93)
xyzzyaaal93=xyzzyaaal93+real(xyzzyaaaf93,dp)*xyzzyaaak93*(xyzzyaaoy1(x&
&yzzyaaab93)-1.d0)
enddo
xyzzyaaan93=xyzzyaaal93*twopi/volume
endif
if(xyzzyaaao93)then
xc_corr=xyzzyaaam93-xyzzyaaan93
else
xc_corr=xyzzyaaam93-xyzzyaaan93-self_term
endif
case default
call errstop('FINITE_SIZE_CORR_XC','Unknown method type XC_CORR_METHOD&
&')
end select
deallocate(xyzzyaaah93,xyzzyaaou1,xyzzyaaov1,xyzzyaaot1)
if(allocated(xyzzyaaoy1))deallocate(xyzzyaaoy1)
end subroutine finite_size_corr_xc
subroutine xyzzyaase1
use toms573,only : nl2sol,dfault
implicit none
integer xyzzyaaaa94,xyzzyaaab94
integer,allocatable :: xyzzyaaac94(:)
real(dp) xyzzyaaad94,xyzzyaaae94,xyzzyaaaf94
real(dp),allocatable :: xyzzyaaag94(:),xyzzyaaah94(:)
logical,parameter :: xyzzyaaai94=.false.
external calcr_xc_corr,calcj_xc_corr
nparam1=8
nparam2=6
xyzzyaaos1=nparam1+nparam2
xyzzyaaor1=xyzzyaaoq1
if(xyzzyaaoz1)then
allocate(weight(xyzzyaaor1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'FIT_STRUCTURE_FACTOR','weights')
weight(1)=xyzzyaaou1(2)
do xyzzyaaaa94=2,xyzzyaaor1-1
weight(xyzzyaaaa94)=abs(xyzzyaaou1(xyzzyaaaa94+1)-xyzzyaaou1(xyzzyaaaa&
&94-1))
enddo
weight(xyzzyaaor1)=abs(xyzzyaaou1(xyzzyaaor1)-xyzzyaaou1(xyzzyaaor1-1)&
&)
xyzzyaaaf94=maxval(weight(:))
weight(:)=weight(:)/xyzzyaaaf94
endif
allocate(xyzzyaaac94(60+xyzzyaaos1),xyzzyaaag94(93+xyzzyaaor1*(xyzzyaa&
&os1+3)+xyzzyaaos1*(3*xyzzyaaos1+33)/2),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'FIT_STRUCTURE_FACTOR','NL2SOL arrays')
call dfault(xyzzyaaac94,xyzzyaaag94)
xyzzyaaac94(14:15)=0
xyzzyaaac94(19:24)=0
xyzzyaaac94(21)=-1
xyzzyaaac94(17)=500
xyzzyaaac94(18)=100
allocate(a(xyzzyaaos1),xyzzyaaow1(xyzzyaaor1),xyzzyaaox1(xyzzyaaor1),x&
&yzzyaaoy1(xyzzyaaor1),stat=xyzzyaaad1)
call check_alloc(xyzzyaaad1,'FIT_STRUCTURE_FACTOR','work arrays')
a(:)=0.d0
call wout()
call wout('Evaluating XC finite size correction using XC_CORR_METHOD==&
&2')
call nl2sol(xyzzyaaor1,xyzzyaaos1,a,calcr_xc_corr,calcj_xc_corr,xyzzya&
&aac94,xyzzyaaag94)
if(xyzzyaaac94(1)==3)then
call wout('Parameter convergence. Structure factor fit complete.')
elseif(xyzzyaaac94(1)==4)then
call wout('Relative-function convergence. Structure factor fit complet&
&e.')
elseif(xyzzyaaac94(1)==5)then
call wout('Parameter and relative-function convergence. Structure fact&
&or fit complete.')
elseif(xyzzyaaac94(1)==6)then
call wout('Absolute function convergence. Structure factor fit complet&
&e.')
elseif(xyzzyaaac94(1)==7)then
call wout('Singular convergence. Structure factor fit complete.')
elseif(xyzzyaaac94(1)==8)then
call wout('False convergence. Structure factor fit complete.')
elseif(xyzzyaaac94(1)==9.or.xyzzyaaac94(1)==10)then
call wout('Maximum number of iterations reached. Structure factor fit &
&complete.')
else
call wout('NL2SOL return code: '//trim(i2s(xyzzyaaac94(1))))
call wout('Structure factor fitting procedure was unsuccessful.')
endif
chi_squared_sf=0.d0
do xyzzyaaab94=1,xyzzyaaor1
xyzzyaaad94=0.d0
do xyzzyaaaa94=1,nparam1
xyzzyaaad94=xyzzyaaad94+a(xyzzyaaaa94)*xyzzyaaou1(xyzzyaaab94)**(xyzzy&
&aaaa94+1)
enddo
xyzzyaaoy1(xyzzyaaab94)=1.d0-exp(xyzzyaaad94)
xyzzyaaad94=0.d0
do xyzzyaaaa94=1,nparam2
xyzzyaaad94=xyzzyaaad94+a(nparam1+xyzzyaaaa94)*xyzzyaaou1(xyzzyaaab94)&
&**(xyzzyaaaa94+1)
enddo
xyzzyaaoy1(xyzzyaaab94)=xyzzyaaoy1(xyzzyaaab94)*(xyzzyaaad94+1.d0)
chi_squared_sf=chi_squared_sf+(xyzzyaaot1(xyzzyaaab94)-xyzzyaaoy1(xyzz&
&yaaab94))**2
enddo
if(xyzzyaaai94)then
allocate(xyzzyaaah94(200))
xyzzyaaae94=xyzzyaaou1(xyzzyaaor1)/200.d0
xyzzyaaah94=0.d0
xyzzyaaaf94=0.d0
do xyzzyaaab94=1,200
xyzzyaaad94=0.d0
do xyzzyaaaa94=1,nparam1
xyzzyaaad94=xyzzyaaad94+a(xyzzyaaaa94)*xyzzyaaaf94**(xyzzyaaaa94+1)
enddo
xyzzyaaah94(xyzzyaaab94)=1.d0-exp(xyzzyaaad94)
xyzzyaaad94=0.d0
do xyzzyaaaa94=1,nparam2
xyzzyaaad94=xyzzyaaad94+a(nparam1+xyzzyaaaa94)*xyzzyaaaf94**(xyzzyaaaa&
&94+1)
enddo
xyzzyaaah94(xyzzyaaab94)=xyzzyaaah94(xyzzyaaab94)*(xyzzyaaad94+1.d0)
xyzzyaaaf94=xyzzyaaaf94+xyzzyaaae94
enddo
call wout('DATA')
do xyzzyaaaa94=1,xyzzyaaoq1
call wout('',(/xyzzyaaou1(xyzzyaaaa94),xyzzyaaot1(xyzzyaaaa94)/))
enddo
call wout('&')
xyzzyaaaf94=0.d0
do xyzzyaaaa94=1,200
call wout('',(/xyzzyaaaf94,xyzzyaaah94(xyzzyaaaa94)/))
xyzzyaaaf94=xyzzyaaaf94+xyzzyaaae94
enddo
deallocate(xyzzyaaah94)
endif
deallocate(xyzzyaaac94,xyzzyaaag94,xyzzyaaow1,xyzzyaaox1)
if(xyzzyaaoz1)deallocate(weight)
end subroutine xyzzyaase1
subroutine eval_residuals_xc_corr(ai,residual)
implicit none
real(dp),intent(in) :: ai(:)
real(dp),intent(out) :: residual(:)
integer xyzzyaaaa95,xyzzyaaab95
real(dp) xyzzyaaac95,xyzzyaaad95,xyzzyaaae95
do xyzzyaaaa95=1,xyzzyaaor1
xyzzyaaac95=0.d0
do xyzzyaaab95=1,nparam1
xyzzyaaac95=xyzzyaaac95+ai(xyzzyaaab95)*xyzzyaaou1(xyzzyaaaa95)**(xyzz&
&yaaab95+1)
enddo
xyzzyaaad95=0.d0
do xyzzyaaab95=1,nparam2
xyzzyaaad95=xyzzyaaad95+ai(nparam1+xyzzyaaab95)*xyzzyaaou1(xyzzyaaaa95&
&)**(xyzzyaaab95+1)
enddo
xyzzyaaae95=exp_protect(xyzzyaaac95)
residual(xyzzyaaaa95)=xyzzyaaot1(xyzzyaaaa95)-((1.d0-xyzzyaaae95)*(1.d&
&0+xyzzyaaad95))
enddo
if(xyzzyaaoz1)residual(1:xyzzyaaor1)=residual(1:xyzzyaaor1)*weight(1:x&
&yzzyaaor1)
end subroutine eval_residuals_xc_corr
subroutine eval_jacobian_xc_corr(ai,jacobian)
implicit none
real(dp),intent(in) :: ai(:)
real(dp),intent(out) :: jacobian(:)
integer xyzzyaaaa96,xyzzyaaab96,xyzzyaaac96,xyzzyaaad96
xyzzyaaow1(:)=0.d0
do xyzzyaaaa96=2,xyzzyaaor1
do xyzzyaaab96=1,nparam1
xyzzyaaow1(xyzzyaaaa96)=xyzzyaaow1(xyzzyaaaa96)+ai(xyzzyaaab96)*xyzzya&
&aou1(xyzzyaaaa96)**(xyzzyaaab96+1)
enddo
enddo
do xyzzyaaaa96=1,xyzzyaaor1
xyzzyaaow1(xyzzyaaaa96)=exp_protect(xyzzyaaow1(xyzzyaaaa96))
enddo
xyzzyaaox1(:)=0.d0
do xyzzyaaaa96=2,xyzzyaaor1
do xyzzyaaab96=1,nparam2
xyzzyaaox1(xyzzyaaaa96)=xyzzyaaox1(xyzzyaaaa96)+ai(nparam1+xyzzyaaab96&
&)*xyzzyaaou1(xyzzyaaaa96)**(xyzzyaaab96+1)
enddo
enddo
xyzzyaaac96=0
do xyzzyaaaa96=1,nparam1
do xyzzyaaab96=1,xyzzyaaor1
xyzzyaaac96=xyzzyaaac96+1
jacobian(xyzzyaaac96)=(1.d0+xyzzyaaox1(xyzzyaaab96))*xyzzyaaow1(xyzzya&
&aab96)*xyzzyaaou1(xyzzyaaab96)**(xyzzyaaaa96+1)
enddo
enddo
xyzzyaaad96=xyzzyaaac96
if(xyzzyaaoz1)jacobian(1:xyzzyaaac96)=jacobian(1:xyzzyaaac96)*weight(x&
&yzzyaaab96)
do xyzzyaaaa96=1,nparam2
do xyzzyaaab96=1,xyzzyaaor1
xyzzyaaac96=xyzzyaaac96+1
jacobian(xyzzyaaac96)=(xyzzyaaow1(xyzzyaaab96)-1.d0)*xyzzyaaou1(xyzzya&
&aab96)**(xyzzyaaaa96+1)
enddo
enddo
if(xyzzyaaoz1)jacobian(xyzzyaaad96:xyzzyaaac96)=jacobian(xyzzyaaad96:x&
&yzzyaaac96)*weight(xyzzyaaab96)
end subroutine eval_jacobian_xc_corr
subroutine xyzzyaasf1(r,n,r_c,k_0,exp_step,weight,f)
implicit none
integer,intent(in) :: n
real(dp),intent(in) :: r,r_c,k_0,exp_step
real(dp),intent(inout) :: f(n),weight
real(dp) xyzzyaaaa97
integer xyzzyaaab97
select case(xyzzyaakl1)
case(xyzzyaakm1)
weight=1.d0
xyzzyaaaa97=r_c/dble(n)
f=0.d0
xyzzyaaab97=int(r/xyzzyaaaa97)+1
f(xyzzyaaab97)=1.d0
case(xyzzyaakn1)
weight=1.d0
call sph_bessel(r,n,r_c,f)
case(xyzzyaako1)
weight=1.d0
call chebyshev(r,n,r_c,f)
case(xyzzyaakp1)
weight=1.d0
call exponential(r,n,k_0,exp_step,f)
case default
call errstop('FIN_DEN_BASIS_FNC','Basis type not recognised.')
end select
end subroutine xyzzyaasf1
subroutine backup_expval_file
integer xyzzyaaaa98,xyzzyaaab98
logical xyzzyaaac98
character(80) expval_line
if(am_master)then
call timer('BACKUP_EXPVAL_FILE',.true.)
inquire(file='expval.data',exist=xyzzyaaac98)
if(.not.xyzzyaaac98)return
call open_units(xyzzyaaaa98,xyzzyaaae1)
if(xyzzyaaae1/=0)call errstop('BACKUP_EXPVAL_FILE','Could not find fre&
&e i/o unit for expval.data file.')
open(unit=xyzzyaaaa98,file="expval.data",form="formatted",status="old"&
&,action="read")
call open_units(xyzzyaaab98,xyzzyaaae1)
if(xyzzyaaae1/=0)call errstop('BACKUP_EXPVAL_FILE','Could not find fre&
&e i/o unit for expval.backup file.')
open(unit=xyzzyaaab98,file="expval.backup",form="formatted",status="re&
&place",action="write")
do
read(unit=xyzzyaaaa98,fmt="(a)",end=10)expval_line
write(unit=xyzzyaaab98,fmt="(a)")trim(expval_line)
enddo
10 close(unit=xyzzyaaaa98)
close(unit=xyzzyaaab98)
open_unit(xyzzyaaaa98)=.false.
open_unit(xyzzyaaab98)=.false.
call qmc_barrier
call timer('BACKUP_EXPVAL_FILE',.false.)
else
call qmc_barrier
endif
end subroutine backup_expval_file
subroutine zero_expval
implicit none
if(xyzzyaaba1)then
if(density)then
xyzzyaadk1(:,:)=czero
endif
if(spin_density)then
xyzzyaadt1(:,:)=czero
endif
if(spin_density_mat)then
xyzzyaaed1(:,:,:)=czero
endif
if(pair_corr)then
if(.not.spin_density)xyzzyaadt1(:,:)=czero
expval_pcf(:,:)=czero
endif
if(structure_factor)then
xyzzyaaga1(:,:)=0.d0
xyzzyaafv1=0.d0
if(.not.homogeneous_density)xyzzyaagg1(:,:)=czero
endif
endif
if(pair_corr_sph)then
xyzzyaafc1(:,:,:)=0.d0
xyzzyaaey1=0.d0
endif
if(structure_factor_sph)then
xyzzyaagu1(:,:)=0.d0
xyzzyaagq1(:)=0.d0
endif
if(loc_tensor)then
xyzzyaafl1(:,:,:)=czero
xyzzyaafh1(:)=0.d0
endif
if(onep_density_mat.or.mom_den)then
xyzzyaahi1(:,:)=0.d0
xyzzyaahj1(:,:)=0.d0
xyzzyaahf1(:,:)=0.d0
if(mom_den)then
xyzzyaahz1(:,:)=0.d0
xyzzyaaia1(:,:)=0.d0
xyzzyaahw1(:)=0.d0
endif
endif
if(twop_density_mat)then
xyzzyaaiz1(:,:)=0.d0
xyzzyaaja1(:,:)=0.d0
xyzzyaaiw1(:,:)=0.d0
endif
if(cond_fraction)then
xyzzyaajy1(:,:)=0.d0
xyzzyaajz1(:,:)=0.d0
xyzzyaajv1(:,:)=0.d0
endif
if(finite_density)then
xyzzyaaks1=0.d0
xyzzyaakw1=0.d0
xyzzyaala1=0.d0
xyzzyaale1(:,:)=0.d0
xyzzyaalf1(:,:)=0.d0
xyzzyaalm1(:,:)=0.d0
xyzzyaaln1(:,:)=0.d0
endif
if(population)then
xyzzyaamt1(:,:)=0.d0
xyzzyaamu1(:,:)=0.d0
xyzzyaamq1(:)=0.d0
endif
if(mol_density)then
xyzzyaalz1=0.d0
xyzzyaamd1=0.d0
xyzzyaamh1(:,:,:,:)=0.d0
endif
if(twop_dm_mom)then
xyzzyaanl1(:,:)=0.d0
xyzzyaanm1(:,:)=0.d0
xyzzyaani1(:)=0.d0
endif
if(cond_frac_mom)then
xyzzyaaof1(:,:)=0.d0
xyzzyaaog1(:,:)=0.d0
xyzzyaaoc1(:)=0.d0
endif
end subroutine zero_expval
subroutine deallocate_expval
implicit none
if(allocated(xyzzyaaco1))deallocate(xyzzyaaco1)
if(allocated(xyzzyaach1))deallocate(xyzzyaach1)
if(allocated(xyzzyaacn1))deallocate(xyzzyaacn1)
if(allocated(xyzzyaadb1))deallocate(xyzzyaadb1)
if(allocated(expval_ngvec))deallocate(expval_ngvec)
if(allocated(xyzzyaadc1))deallocate(xyzzyaadc1)
if(allocated(xyzzyaadc1))deallocate(xyzzyaadc1)
if(allocated(xyzzyaadd1))deallocate(xyzzyaadd1)
if(allocated(expval_ngvec_truncated))deallocate(expval_ngvec_truncated&
&)
if(allocated(expval_gvec))deallocate(expval_gvec)
if(allocated(xyzzyaacy1))deallocate(xyzzyaacy1)
if(allocated(xyzzyaacz1))deallocate(xyzzyaacz1)
if(allocated(xyzzyaadf1))deallocate(xyzzyaadf1)
if(allocated(xyzzyaade1))deallocate(xyzzyaade1)
if(allocated(xyzzyaadj1))deallocate(xyzzyaadj1)
if(allocated(xyzzyaadi1))deallocate(xyzzyaadi1)
if(allocated(xyzzyaadh1))deallocate(xyzzyaadh1)
if(allocated(xyzzyaadl1))deallocate(xyzzyaadl1)
if(allocated(xyzzyaadk1))deallocate(xyzzyaadk1)
if(allocated(xyzzyaadm1))deallocate(xyzzyaadm1)
if(allocated(expval_den_mpc))deallocate(expval_den_mpc)
if(allocated(xyzzyaadn1))deallocate(xyzzyaadn1)
if(allocated(xyzzyaads1))deallocate(xyzzyaads1)
if(allocated(xyzzyaadr1))deallocate(xyzzyaadr1)
if(allocated(xyzzyaadq1))deallocate(xyzzyaadq1)
if(allocated(xyzzyaadu1))deallocate(xyzzyaadu1)
if(allocated(xyzzyaadt1))deallocate(xyzzyaadt1)
if(allocated(xyzzyaadv1))deallocate(xyzzyaadv1)
if(allocated(xyzzyaadw1))deallocate(xyzzyaadw1)
if(allocated(xyzzyaaec1))deallocate(xyzzyaaec1)
if(allocated(xyzzyaaeb1))deallocate(xyzzyaaeb1)
if(allocated(xyzzyaaee1))deallocate(xyzzyaaee1)
if(allocated(xyzzyaaed1))deallocate(xyzzyaaed1)
if(allocated(xyzzyaaef1))deallocate(xyzzyaaef1)
if(allocated(xyzzyaaem1))deallocate(xyzzyaaem1)
if(allocated(xyzzyaael1))deallocate(xyzzyaael1)
if(allocated(xyzzyaaek1))deallocate(xyzzyaaek1)
if(allocated(xyzzyaaen1))deallocate(xyzzyaaen1)
if(allocated(expval_pcf))deallocate(expval_pcf)
if(allocated(xyzzyaaeo1))deallocate(xyzzyaaeo1)
if(allocated(xyzzyaaep1))deallocate(xyzzyaaep1)
if(allocated(xyzzyaafb1))deallocate(xyzzyaafb1)
if(allocated(xyzzyaaev1))deallocate(xyzzyaaev1)
if(allocated(xyzzyaafd1))deallocate(xyzzyaafd1)
if(allocated(xyzzyaafc1))deallocate(xyzzyaafc1)
if(allocated(xyzzyaafe1))deallocate(xyzzyaafe1)
if(allocated(xyzzyaaff1))deallocate(xyzzyaaff1)
if(allocated(xyzzyaafm1))deallocate(xyzzyaafm1)
if(allocated(xyzzyaafi1))deallocate(xyzzyaafi1)
if(allocated(xyzzyaafl1))deallocate(xyzzyaafl1)
if(allocated(xyzzyaafh1))deallocate(xyzzyaafh1)
if(allocated(xyzzyaafn1))deallocate(xyzzyaafn1)
if(allocated(xyzzyaafj1))deallocate(xyzzyaafj1)
if(allocated(xyzzyaafo1))deallocate(xyzzyaafo1)
if(allocated(xyzzyaafk1))deallocate(xyzzyaafk1)
if(allocated(xyzzyaafy1))deallocate(xyzzyaafy1)
if(allocated(xyzzyaaft1))deallocate(xyzzyaaft1)
if(allocated(xyzzyaafz1))deallocate(xyzzyaafz1)
if(allocated(xyzzyaaga1))deallocate(xyzzyaaga1)
if(allocated(xyzzyaagf1))deallocate(xyzzyaagf1)
if(allocated(xyzzyaagb1))deallocate(xyzzyaagb1)
if(allocated(xyzzyaagc1))deallocate(xyzzyaagc1)
if(allocated(xyzzyaagd1))deallocate(xyzzyaagd1)
if(allocated(xyzzyaafu1))deallocate(xyzzyaafu1)
if(allocated(xyzzyaage1))deallocate(xyzzyaage1)
if(allocated(xyzzyaagh1))deallocate(xyzzyaagh1)
if(allocated(xyzzyaagg1))deallocate(xyzzyaagg1)
if(allocated(xyzzyaagi1))deallocate(xyzzyaagi1)
if(allocated(xyzzyaagj1))deallocate(xyzzyaagj1)
if(allocated(xyzzyaagq1))deallocate(xyzzyaagq1)
if(allocated(xyzzyaags1))deallocate(xyzzyaags1)
if(allocated(xyzzyaagt1))deallocate(xyzzyaagt1)
if(allocated(xyzzyaagn1))deallocate(xyzzyaagn1)
if(allocated(xyzzyaagu1))deallocate(xyzzyaagu1)
if(allocated(xyzzyaagv1))deallocate(xyzzyaagv1)
if(allocated(xyzzyaagy1))deallocate(xyzzyaagy1)
if(allocated(xyzzyaagx1))deallocate(xyzzyaagx1)
if(allocated(xyzzyaagw1))deallocate(xyzzyaagw1)
if(allocated(xyzzyaagr1))deallocate(xyzzyaagr1)
if(allocated(xyzzyaahg1))deallocate(xyzzyaahg1)
if(allocated(xyzzyaahk1))deallocate(xyzzyaahk1)
if(allocated(xyzzyaahl1))deallocate(xyzzyaahl1)
if(allocated(xyzzyaahf1))deallocate(xyzzyaahf1)
if(allocated(xyzzyaaho1))deallocate(xyzzyaaho1)
if(allocated(xyzzyaahi1))deallocate(xyzzyaahi1)
if(allocated(xyzzyaahj1))deallocate(xyzzyaahj1)
if(allocated(xyzzyaahh1))deallocate(xyzzyaahh1)
if(allocated(xyzzyaahm1))deallocate(xyzzyaahm1)
if(allocated(xyzzyaahn1))deallocate(xyzzyaahn1)
if(allocated(xyzzyaahc1))deallocate(xyzzyaahc1)
if(allocated(xyzzyaahd1))deallocate(xyzzyaahd1)
if(allocated(xyzzyaahr1))deallocate(xyzzyaahr1)
if(allocated(xyzzyaahp1))deallocate(xyzzyaahp1)
if(allocated(xyzzyaahq1))deallocate(xyzzyaahq1)
if(allocated(xyzzyaaix1))deallocate(xyzzyaaix1)
if(allocated(xyzzyaajb1))deallocate(xyzzyaajb1)
if(allocated(xyzzyaajc1))deallocate(xyzzyaajc1)
if(allocated(xyzzyaaiw1))deallocate(xyzzyaaiw1)
if(allocated(xyzzyaaiz1))deallocate(xyzzyaaiz1)
if(allocated(xyzzyaaja1))deallocate(xyzzyaaja1)
if(allocated(xyzzyaaiy1))deallocate(xyzzyaaiy1)
if(allocated(xyzzyaajd1))deallocate(xyzzyaajd1)
if(allocated(xyzzyaaje1))deallocate(xyzzyaaje1)
if(allocated(xyzzyaain1))deallocate(xyzzyaain1)
if(allocated(xyzzyaaio1))deallocate(xyzzyaaio1)
if(allocated(xyzzyaaip1))deallocate(xyzzyaaip1)
if(allocated(xyzzyaais1))deallocate(xyzzyaais1)
if(allocated(xyzzyaaiq1))deallocate(xyzzyaaiq1)
if(allocated(xyzzyaair1))deallocate(xyzzyaair1)
if(allocated(xyzzyaait1))deallocate(xyzzyaait1)
if(allocated(xyzzyaaji1))deallocate(xyzzyaaji1)
if(allocated(xyzzyaajh1))deallocate(xyzzyaajh1)
if(allocated(xyzzyaajf1))deallocate(xyzzyaajf1)
if(allocated(xyzzyaajg1))deallocate(xyzzyaajg1)
if(allocated(xyzzyaajw1))deallocate(xyzzyaajw1)
if(allocated(xyzzyaaka1))deallocate(xyzzyaaka1)
if(allocated(xyzzyaakb1))deallocate(xyzzyaakb1)
if(allocated(xyzzyaajv1))deallocate(xyzzyaajv1)
if(allocated(xyzzyaajy1))deallocate(xyzzyaajy1)
if(allocated(xyzzyaajz1))deallocate(xyzzyaajz1)
if(allocated(xyzzyaajx1))deallocate(xyzzyaajx1)
if(allocated(xyzzyaakc1))deallocate(xyzzyaakc1)
if(allocated(xyzzyaakd1))deallocate(xyzzyaakd1)
if(allocated(xyzzyaajm1))deallocate(xyzzyaajm1)
if(allocated(xyzzyaajn1))deallocate(xyzzyaajn1)
if(allocated(xyzzyaajo1))deallocate(xyzzyaajo1)
if(allocated(xyzzyaajr1))deallocate(xyzzyaajr1)
if(allocated(xyzzyaajp1))deallocate(xyzzyaajp1)
if(allocated(xyzzyaajq1))deallocate(xyzzyaajq1)
if(allocated(xyzzyaajs1))deallocate(xyzzyaajs1)
if(allocated(xyzzyaakh1))deallocate(xyzzyaakh1)
if(allocated(xyzzyaakg1))deallocate(xyzzyaakg1)
if(allocated(xyzzyaake1))deallocate(xyzzyaake1)
if(allocated(xyzzyaakf1))deallocate(xyzzyaakf1)
if(allocated(xyzzyaahx1))deallocate(xyzzyaahx1)
if(allocated(xyzzyaaib1))deallocate(xyzzyaaib1)
if(allocated(xyzzyaaic1))deallocate(xyzzyaaic1)
if(allocated(xyzzyaahw1))deallocate(xyzzyaahw1)
if(allocated(xyzzyaahz1))deallocate(xyzzyaahz1)
if(allocated(xyzzyaaia1))deallocate(xyzzyaaia1)
if(allocated(xyzzyaahy1))deallocate(xyzzyaahy1)
if(allocated(xyzzyaaid1))deallocate(xyzzyaaid1)
if(allocated(xyzzyaaie1))deallocate(xyzzyaaie1)
if(allocated(xyzzyaahu1))deallocate(xyzzyaahu1)
if(allocated(xyzzyaahv1))deallocate(xyzzyaahv1)
if(allocated(xyzzyaaii1))deallocate(xyzzyaaii1)
if(allocated(xyzzyaaig1))deallocate(xyzzyaaig1)
if(allocated(xyzzyaaih1))deallocate(xyzzyaaih1)
if(allocated(xyzzyaaif1))deallocate(xyzzyaaif1)
if(allocated(xyzzyaaij1))deallocate(xyzzyaaij1)
if(allocated(xyzzyaalg1))deallocate(xyzzyaalg1)
if(allocated(xyzzyaalh1))deallocate(xyzzyaalh1)
if(allocated(xyzzyaale1))deallocate(xyzzyaale1)
if(allocated(xyzzyaalf1))deallocate(xyzzyaalf1)
if(allocated(xyzzyaali1))deallocate(xyzzyaali1)
if(allocated(xyzzyaalj1))deallocate(xyzzyaalj1)
if(allocated(xyzzyaalo1))deallocate(xyzzyaalo1)
if(allocated(xyzzyaalp1))deallocate(xyzzyaalp1)
if(allocated(xyzzyaalm1))deallocate(xyzzyaalm1)
if(allocated(xyzzyaaln1))deallocate(xyzzyaaln1)
if(allocated(xyzzyaalq1))deallocate(xyzzyaalq1)
if(allocated(xyzzyaalr1))deallocate(xyzzyaalr1)
if(allocated(xyzzyaalk1))deallocate(xyzzyaalk1)
if(allocated(xyzzyaall1))deallocate(xyzzyaall1)
if(allocated(xyzzyaals1))deallocate(xyzzyaals1)
if(allocated(xyzzyaalt1))deallocate(xyzzyaalt1)
if(allocated(xyzzyaalu1))deallocate(xyzzyaalu1)
if(allocated(xyzzyaamr1))deallocate(xyzzyaamr1)
if(allocated(xyzzyaamv1))deallocate(xyzzyaamv1)
if(allocated(xyzzyaamw1))deallocate(xyzzyaamw1)
if(allocated(xyzzyaamq1))deallocate(xyzzyaamq1)
if(allocated(xyzzyaamt1))deallocate(xyzzyaamt1)
if(allocated(xyzzyaamu1))deallocate(xyzzyaamu1)
if(allocated(xyzzyaams1))deallocate(xyzzyaams1)
if(allocated(xyzzyaamx1))deallocate(xyzzyaamx1)
if(allocated(xyzzyaamy1))deallocate(xyzzyaamy1)
if(allocated(xyzzyaamo1))deallocate(xyzzyaamo1)
if(allocated(xyzzyaamp1))deallocate(xyzzyaamp1)
if(allocated(xyzzyaanb1)) deallocate(xyzzyaanb1)
if(allocated(xyzzyaamz1))deallocate(xyzzyaamz1)
if(allocated(xyzzyaana1))deallocate(xyzzyaana1)
if(allocated(xyzzyaami1))deallocate(xyzzyaami1)
if(allocated(xyzzyaamh1))deallocate(xyzzyaamh1)
if(allocated(xyzzyaamj1))deallocate(xyzzyaamj1)
if(allocated(xyzzyaamk1))deallocate(xyzzyaamk1)
if(allocated(xyzzyaaml1))deallocate(xyzzyaaml1)
if(allocated(xyzzyaanj1))deallocate(xyzzyaanj1)
if(allocated(xyzzyaann1))deallocate(xyzzyaann1)
if(allocated(xyzzyaano1))deallocate(xyzzyaano1)
if(allocated(xyzzyaani1))deallocate(xyzzyaani1)
if(allocated(xyzzyaanl1))deallocate(xyzzyaanl1)
if(allocated(xyzzyaanm1))deallocate(xyzzyaanm1)
if(allocated(xyzzyaank1))deallocate(xyzzyaank1)
if(allocated(xyzzyaanp1))deallocate(xyzzyaanp1)
if(allocated(xyzzyaanq1))deallocate(xyzzyaanq1)
if(allocated(xyzzyaane1))deallocate(xyzzyaane1)
if(allocated(xyzzyaanf1))deallocate(xyzzyaanf1)
if(allocated(xyzzyaang1))deallocate(xyzzyaang1)
if(allocated(xyzzyaanv1))deallocate(xyzzyaanv1)
if(allocated(xyzzyaant1))deallocate(xyzzyaant1)
if(allocated(xyzzyaanr1))deallocate(xyzzyaanr1)
if(allocated(xyzzyaans1))deallocate(xyzzyaans1)
if(allocated(xyzzyaaod1))deallocate(xyzzyaaod1)
if(allocated(xyzzyaaoh1))deallocate(xyzzyaaoh1)
if(allocated(xyzzyaaoi1))deallocate(xyzzyaaoi1)
if(allocated(xyzzyaaoc1))deallocate(xyzzyaaoc1)
if(allocated(xyzzyaaof1))deallocate(xyzzyaaof1)
if(allocated(xyzzyaaog1))deallocate(xyzzyaaog1)
if(allocated(xyzzyaaoe1))deallocate(xyzzyaaoe1)
if(allocated(xyzzyaaoj1))deallocate(xyzzyaaoj1)
if(allocated(xyzzyaaok1))deallocate(xyzzyaaok1)
if(allocated(xyzzyaany1))deallocate(xyzzyaany1)
if(allocated(xyzzyaanz1))deallocate(xyzzyaanz1)
if(allocated(xyzzyaaoa1))deallocate(xyzzyaaoa1)
if(allocated(xyzzyaaop1))deallocate(xyzzyaaop1)
if(allocated(xyzzyaaon1))deallocate(xyzzyaaon1)
if(allocated(xyzzyaaol1))deallocate(xyzzyaaol1)
if(allocated(xyzzyaaom1))deallocate(xyzzyaaom1)
end subroutine deallocate_expval
end module slaarnaaq
subroutine calcr_xc_corr(ndata,nparam,ai,nf,res,uiparm,urparm,ufparm)
use dsp
use slaarnaaq,only : eval_residuals_xc_corr
implicit none
integer,intent(in) :: ndata,nparam
integer,intent(inout) :: nf
integer,intent(inout),optional :: uiparm(:)
real(dp),intent(in) :: ai(:)
real(dp),intent(out) :: res(:)
real(dp),intent(inout),optional :: urparm(:),ufparm
call eval_residuals_xc_corr(ai,res)
end subroutine calcr_xc_corr
subroutine calcj_xc_corr(ndata,nparam,ai,nf,jac,uiparm,urparm,ufparm)
use dsp
use slaarnaaq,only : eval_jacobian_xc_corr
implicit none
integer,intent(in) :: ndata,nparam
integer,intent(inout) :: nf
integer,intent(inout),optional :: uiparm(:)
real(dp),intent(in) :: ai(:)
real(dp),intent(out) :: jac(:)
real(dp),intent(inout),optional :: urparm(:),ufparm
call eval_jacobian_xc_corr(ai,jac)
end subroutine calcj_xc_corr
function calcf_xc_corr(xx)
use dsp
use slaarnaaq,only : a,nparam1,nparam2
implicit none
real(dp),intent(in) :: xx
integer xyzzyaaaa103
real(dp) xyzzyaaab103,xyzzyaaac103,calcf_xc_corr
xyzzyaaab103=0.d0
do xyzzyaaaa103=1,nparam1
xyzzyaaab103=xyzzyaaab103+a(xyzzyaaaa103)*xx**(xyzzyaaaa103+1)
enddo
xyzzyaaac103=0.d0
do xyzzyaaaa103=1,nparam2
xyzzyaaac103=xyzzyaaac103+a(nparam1+xyzzyaaaa103)*xx**(xyzzyaaaa103+1)
enddo
calcf_xc_corr=(1.d0-exp(xyzzyaaab103))*(1.d0+xyzzyaaac103)
end function calcf_xc_corr
function calcf_xc_corr2(xx)
use dsp
use slaarnaaq,only : a,nparam1,nparam2
implicit none
real(dp),intent(in) :: xx
integer xyzzyaaaa104
real(dp) xyzzyaaab104,xyzzyaaac104,calcf_xc_corr2
xyzzyaaab104=0.d0
do xyzzyaaaa104=1,nparam1
xyzzyaaab104=xyzzyaaab104+a(xyzzyaaaa104)*xx**(xyzzyaaaa104+1)
enddo
xyzzyaaac104=0.d0
do xyzzyaaaa104=1,nparam2
xyzzyaaac104=xyzzyaaac104+a(nparam1+xyzzyaaaa104)*xx**(xyzzyaaaa104+1)
enddo
calcf_xc_corr2=(1.d0-exp(xyzzyaaab104))*(1.d0+xyzzyaaac104)-1.d0
end function calcf_xc_corr2

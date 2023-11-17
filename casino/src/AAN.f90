module slaarnaan
use slaarnaag
use dsp
use slaarnabg
use slaarnabq
use store
use format_utils,only : wout,i2s,wordwrap
use slaarnaat,   only : cell_map,gsgrid,gs_origin,molgscreening,cusp_c&
&orrection
use run_control, only : errstop,errstop_master,timer,errwarn,check_all&
&oc
use slaarnabt, only : inverse3,mxmb
implicit none
private
public init_geometry,print_geometry,ei_distances,ee_distances,ee_dista&
&nces_all,ei_distances_all,gvsort,lattice_generator,neighbour_analysis&
&,mol_system_size,make_movie,netot_nitot_products,ee_kato_gamma,en_kat&
&o_gamma,check_kpoints,points,vector_difference,atoms_label_species,at&
&oms_label_pcell,map_to_simcell
contains
subroutine points(rele,sele,initial_rele,initial_rele_set,printout)
use slaarnaap,         only : expot_is_finite
use parallel,      only : am_master
use slaarnaca,         only : atnum,zion
use slaarnacc,only : ranx,ranx_gaussian
use slaarnaci,           only : sdw_orb_theta
implicit none
integer, intent(inout) :: sele(netot)
real(dp),intent(in) :: initial_rele(3,netot)
real(dp),intent(inout) :: rele(3,netot)
logical,intent(in) :: initial_rele_set(netot)
logical,intent(in),optional :: printout
integer xyzzyaaaa2,xyzzyaaab2,xyzzyaaac2,xyzzyaaad2,xyzzyaaae2,xyzzyaa&
&af2,xyzzyaaag2,xyzzyaaah2,xyzzyaaai2,xyzzyaaaj2,xyzzyaaak2,xyzzyaaal2&
&,xyzzyaaam2,xyzzyaaan2
integer,parameter :: xyzzyaaao2=1000
integer,allocatable :: xyzzyaaap2(:)
real(dp) xyzzyaaaq2,xyzzyaaar2,xyzzyaaas2,xyzzyaaat2,xyzzyaaau2,xyzzya&
&aav2,xyzzyaaaw2(3),xyzzyaaax2(3),xyzzyaaay2(3),xyzzyaaaz2(3),xyzzyaab&
&a2,xyzzyaabb2,xyzzyaabc2,xyzzyaabd2,xyzzyaabe2(4,1),xyzzyaabf2,xyzzya&
&abg2
real(dp),parameter :: xyzzyaabh2=2.d0
logical xyzzyaabi2,xyzzyaabj2
logical,allocatable :: xyzzyaabk2(:)
character(80) tmpr
xyzzyaabi2=.false.
if(present(printout))xyzzyaabi2=printout
if(.not.allocated(neion))then
allocate(neion(nitot,nspin),stat=xyzzyaaai2)
call check_alloc(xyzzyaaai2,'POINTS','neion')
endif
xyzzyaaaf2=1
if(nitot>0.and..not.ignore_ionic_interactions)then
if(.not.init_by_ion.and..not.init_by_iontype)then
neion(:,:)=0
allocate(xyzzyaaap2(nspin),stat=xyzzyaaai2)
call check_alloc(xyzzyaaai2,'POINTS','e_remaining')
xyzzyaaap2(1:nspin)=nele(1:nspin)
if(nspin==2.and.electron_system)then
do
do xyzzyaaaa2=1,nitot
do xyzzyaaae2=1,2
if(sum(neion(xyzzyaaaa2,:))<zion(iontype(xyzzyaaaa2)))then
xyzzyaaaf2=xyzzyaaae2
if(mod(xyzzyaaaa2,2)==0)xyzzyaaaf2=abs(xyzzyaaae2-3)
if(xyzzyaaap2(xyzzyaaaf2)>0)then
neion(xyzzyaaaa2,xyzzyaaaf2)=neion(xyzzyaaaa2,xyzzyaaaf2)+1
xyzzyaaap2(xyzzyaaaf2)=xyzzyaaap2(xyzzyaaaf2)-1
endif
endif
enddo
enddo
if(all(xyzzyaaap2(:)==0))exit
if(all(real((neion(:,1)+neion(:,2)),dp)==zion(iontype(:))))exit
enddo
do xyzzyaaae2=1,2
do
if(xyzzyaaap2(xyzzyaaae2)==0)exit
do xyzzyaaaa2=1,nitot
neion(xyzzyaaaa2,xyzzyaaae2)=neion(xyzzyaaaa2,xyzzyaaae2)+1
xyzzyaaap2(xyzzyaaae2)=xyzzyaaap2(xyzzyaaae2)-1
if(xyzzyaaap2(xyzzyaaae2)==0)exit
enddo
enddo
enddo
else
xyzzyaaba2=0.d0
xyzzyaabb2=0.d0
do xyzzyaaan2=1,nitype
if(zion(xyzzyaaan2)>0.d0)then
xyzzyaaba2=xyzzyaaba2+zion(xyzzyaaan2)*dble(nion(xyzzyaaan2))
else
xyzzyaabb2=xyzzyaabb2+zion(xyzzyaaan2)*dble(nion(xyzzyaaan2))
endif
enddo
do xyzzyaaae2=1,nspin
do xyzzyaaaa2=1,nitot
if(pcharge(xyzzyaaae2)<0.d0.and.zion(iontype(xyzzyaaaa2))>0.d0)then
neion(xyzzyaaaa2,xyzzyaaae2)=int(zion(iontype(xyzzyaaaa2))/xyzzyaaba2*&
&dble(nele(xyzzyaaae2)))
elseif(pcharge(xyzzyaaae2)>0.d0.and.zion(iontype(xyzzyaaaa2))<0.d0)the&
&n
neion(xyzzyaaaa2,xyzzyaaae2)=int(zion(iontype(xyzzyaaaa2))/xyzzyaabb2*&
&dble(nele(xyzzyaaae2)))
endif
xyzzyaaap2(xyzzyaaae2)=xyzzyaaap2(xyzzyaaae2)-neion(xyzzyaaaa2,xyzzyaa&
&ae2)
enddo
enddo
do xyzzyaaae2=1,nspin
do while(xyzzyaaap2(xyzzyaaae2)>0)
xyzzyaabc2=(zion(iontype(1))+sum(neion(1,1:nspin)*pcharge(1:nspin))) *&
&pcharge(xyzzyaaae2)
xyzzyaaal2=1
do xyzzyaaaa2=2,nitot
xyzzyaabd2=(zion(iontype(xyzzyaaaa2))+sum(neion(xyzzyaaaa2,1:nspin)*pc&
&harge(1:nspin))) *pcharge(xyzzyaaae2)
if(xyzzyaabd2<xyzzyaabc2)then
xyzzyaabc2=xyzzyaabd2
xyzzyaaal2=xyzzyaaaa2
endif
enddo
neion(xyzzyaaal2,xyzzyaaae2)=neion(xyzzyaaal2,xyzzyaaae2)+1
xyzzyaaap2(xyzzyaaae2)=xyzzyaaap2(xyzzyaaae2)-1
enddo
enddo
endif
deallocate(xyzzyaaap2)
if(any(neion<0))call errstop('POINTS','Error in points algorithm for g&
&enerating initial distribution. Please file bug report.')
do xyzzyaaae2=1,nspin
if(sum(neion(:,xyzzyaaae2))/=nele(xyzzyaaae2))call errstop('POINTS','E&
&rror in points algorithm for generating initial electron distribution&
&. Please file bug report.')
enddo
endif
do xyzzyaaae2=1,nspin
xyzzyaaab2=0
if(xyzzyaaae2>1)xyzzyaaab2=sum(nele(1:xyzzyaaae2-1))
do xyzzyaaaa2=1,nitot
do xyzzyaaac2=1,neion(xyzzyaaaa2,xyzzyaaae2)
xyzzyaaab2=xyzzyaaab2+1
if(dimensionality==3)then
xyzzyaaaq2=ranx()
xyzzyaaar2=xyzzyaaaq2**third
xyzzyaaas2=ranx()
xyzzyaaat2=twopi*ranx()
xyzzyaaau2=2.d0*xyzzyaabh2*xyzzyaaar2*sqrt(xyzzyaaas2*(1.d0-xyzzyaaas2&
&))
rele(1,xyzzyaaab2)=rion(1,xyzzyaaaa2)+xyzzyaaau2*sin(xyzzyaaat2)
rele(2,xyzzyaaab2)=rion(2,xyzzyaaaa2)+xyzzyaaau2*cos(xyzzyaaat2)
rele(3,xyzzyaaab2)=rion(3,xyzzyaaaa2)+xyzzyaabh2*xyzzyaaar2*(1.d0-2.d0&
&*xyzzyaaas2)
elseif(dimensionality==2)then
xyzzyaaaq2=ranx()
xyzzyaaat2=twopi*ranx()
xyzzyaaau2=sqrt(xyzzyaaaq2)*xyzzyaabh2
rele(1,xyzzyaaab2)=rion(1,xyzzyaaaa2)+xyzzyaaau2*cos(xyzzyaaat2)
rele(2,xyzzyaaab2)=rion(2,xyzzyaaaa2)+xyzzyaaau2*sin(xyzzyaaat2)
rele(3,xyzzyaaab2)=heg_zlayer(heg_layer(which_spin(xyzzyaaab2)))
else
rele(1,xyzzyaaab2)=rion(1,xyzzyaaaa2)+2.d0*xyzzyaabh2*(ranx()-0.5d0)
rele(2,xyzzyaaab2)=heg_ylayer(heg_layer(which_spin(xyzzyaaab2)))
rele(3,xyzzyaaab2)=heg_zlayer(heg_layer(which_spin(xyzzyaaab2)))
endif
enddo
enddo
enddo
elseif(model_system)then
do xyzzyaaae2=1,nspin
if(heg_orbtype(xyzzyaaae2,1)==2)then
xyzzyaaag2=which_ssingle(xyzzyaaae2,spin_dep_wc)
xyzzyaaac2=0
xyzzyaaad2=sum(nele(1:xyzzyaaae2-1))+nuc_nele(xyzzyaaae2)
do xyzzyaaab2=xyzzyaaad2+1,xyzzyaaad2+heg_nele(xyzzyaaae2)
xyzzyaaac2=xyzzyaaac2+1
crystal_attempt: do xyzzyaaam2=0,xyzzyaaao2
if(xyzzyaaam2==xyzzyaaao2)call errstop('POINTS','Cannot place hard sph&
&eres.  Need either a lower density or a larger Gaussian exponent or a&
& better points algorithm.')
if(type_wc==1.or.type_wc==3)then
rele(1,xyzzyaaab2)=site_pos_in_cell(1,xyzzyaaac2,1,heg_slatt(xyzzyaaae&
&2,1)) +ranx_gaussian(sqrt(0.25d0/wc_gauss_exp(xyzzyaaag2)))
else
rele(1,xyzzyaaab2)=site_pos_in_cell(1,xyzzyaaac2,1,heg_slatt(xyzzyaaae&
&2,1))
endif
if(dimensionality==3)then
if(type_wc==1.or.type_wc==3)then
rele(2,xyzzyaaab2)=site_pos_in_cell(2,xyzzyaaac2,1,heg_slatt(xyzzyaaae&
&2,1)) +ranx_gaussian(sqrt(0.25d0/wc_gauss_exp(xyzzyaaag2)))
rele(3,xyzzyaaab2)=site_pos_in_cell(3,xyzzyaaac2,1,heg_slatt(xyzzyaaae&
&2,1)) +ranx_gaussian(sqrt(0.25d0/wc_gauss_exp(xyzzyaaag2)))
else
rele(2,xyzzyaaab2)=site_pos_in_cell(2,xyzzyaaac2,1,heg_slatt(xyzzyaaae&
&2,1))
rele(3,xyzzyaaab2)=site_pos_in_cell(3,xyzzyaaac2,1,heg_slatt(xyzzyaaae&
&2,1))
endif
elseif(dimensionality==2)then
if(type_wc==1.or.type_wc==3)then
rele(2,xyzzyaaab2)=site_pos_in_cell(2,xyzzyaaac2,1,heg_slatt(xyzzyaaae&
&2,1)) +ranx_gaussian(sqrt(0.25d0/wc_gauss_exp(xyzzyaaag2)))
else
rele(2,xyzzyaaab2)=site_pos_in_cell(2,xyzzyaaac2,1,heg_slatt(xyzzyaaae&
&2,1))
endif
rele(3,xyzzyaaab2)=heg_zlayer(heg_layer(xyzzyaaae2))
else
rele(2,xyzzyaaab2)=heg_ylayer(heg_layer(xyzzyaaae2))
rele(3,xyzzyaaab2)=heg_zlayer(heg_layer(xyzzyaaae2))
endif
if(hard_sphere)then
do xyzzyaaak2=1,xyzzyaaab2-1
if(hard_op_spins.and.which_spin(xyzzyaaak2)==which_spin(xyzzyaaab2))cy&
&cle
xyzzyaabe2(1:3,1)=rele(1:3,xyzzyaaak2)-rele(1:3,xyzzyaaab2)
call minimum_image(4,1,xyzzyaabe2)
if(sum(xyzzyaabe2(1:dimensionality,1)**2)<=hard_diam_sq)cycle crystal_&
&attempt
enddo
endif
exit
enddo crystal_attempt
enddo
elseif(heg_orbtype(xyzzyaaae2,1)<0.and.xyzzyaaae2>-heg_orbtype(xyzzyaa&
&ae2,1))then
xyzzyaaad2=0
if(xyzzyaaae2>1)xyzzyaaad2=sum(nele(1:xyzzyaaae2-1))+nuc_nele(xyzzyaaa&
&e2)
xyzzyaaaf2=-heg_orbtype(xyzzyaaae2,1)
xyzzyaaac2=0
if(xyzzyaaaf2>1)xyzzyaaac2=sum(nele(1:xyzzyaaaf2-1))+nuc_nele(xyzzyaaa&
&f2)
do xyzzyaaab2=xyzzyaaad2+1,xyzzyaaad2+heg_nele(xyzzyaaae2)
xyzzyaaac2=xyzzyaaac2+1
rele(1,xyzzyaaab2)=rele(1,xyzzyaaac2)+2.d0*(ranx()-0.5d0)
if(dimensionality==3)then
rele(2,xyzzyaaab2)=rele(2,xyzzyaaac2)+2.d0*(ranx()-0.5d0)
rele(3,xyzzyaaab2)=rele(3,xyzzyaaac2)+2.d0*(ranx()-0.5d0)
elseif(dimensionality==2)then
rele(2,xyzzyaaab2)=rele(2,xyzzyaaac2)+2.d0*(ranx()-0.5d0)
rele(3,xyzzyaaab2)=heg_zlayer(heg_layer(xyzzyaaae2))
else
rele(2,xyzzyaaab2)=heg_ylayer(heg_layer(xyzzyaaae2))
rele(3,xyzzyaaab2)=heg_zlayer(heg_layer(xyzzyaaae2))
endif
enddo
elseif(heg_orbtype(xyzzyaaae2,1)==4)then
xyzzyaaav2=abs((pmass(1)+pmass(3)) /(2*pcharge(1)*pcharge(3)*pmass(1)*&
&pmass(3)))
rele(1,xyzzyaaae2)=ranx_gaussian(xyzzyaaav2)
rele(2,xyzzyaaae2)=ranx_gaussian(xyzzyaaav2)
if(xyzzyaaae2==2.or.xyzzyaaae2==4)rele(1,xyzzyaaae2)=rele(1,xyzzyaaae2&
&)+0.5d0*heg_zlayer(2)
rele(3,xyzzyaaae2)=heg_zlayer(heg_layer(xyzzyaaae2))
elseif(heg_orbtype(xyzzyaaae2,1)==6)then
xyzzyaaav2=0.25d0
rele(1,xyzzyaaae2)=ranx_gaussian(xyzzyaaav2)
rele(2,xyzzyaaae2)=ranx_gaussian(xyzzyaaav2)
rele(3,xyzzyaaae2)=0.d0
else
xyzzyaaaw2=a1
xyzzyaaax2=a2
xyzzyaaay2=a3
xyzzyaaaz2=0.d0
if(periodicity<3)then
xyzzyaaay2(3)=2.d0
xyzzyaaaz2(3)=-1.d0
if(periodicity<2)then
xyzzyaaax2(2)=2.d0
xyzzyaaaz2(2)=-1.d0
if(periodicity==0)then
xyzzyaaaw2(1)=2.d0
xyzzyaaaz2(1)=-1.d0
endif
endif
endif
if(dimensionality==1.and.periodicity==1.and..not.use_expot)xyzzyaaaz2(&
&1)=xyzzyaaaw2(1)*ranx()
xyzzyaaad2=0
if(xyzzyaaae2>1)xyzzyaaad2=sum(nele(1:xyzzyaaae2-1))+nuc_nele(xyzzyaaa&
&e2)
do xyzzyaaab2=xyzzyaaad2+1,xyzzyaaad2+heg_nele(xyzzyaaae2)
fluid_attempt: do xyzzyaaam2=0,xyzzyaaao2
if(xyzzyaaam2==xyzzyaaao2)call errstop('POINTS','Cannot place hard sph&
&eres at this density.  Need a more sophisticated points algorithm or &
&a lower density.')
if(dimensionality==3)then
rele(1:3,xyzzyaaab2)=ranx()*xyzzyaaaw2+ranx()*xyzzyaaax2+ranx()*xyzzya&
&aay2+xyzzyaaaz2
elseif(dimensionality==2)then
rele(1:2,xyzzyaaab2)=ranx()*xyzzyaaaw2(1:2)+ranx()*xyzzyaaax2(1:2)+xyz&
&zyaaaz2(1:2)
rele(3,xyzzyaaab2)=heg_zlayer(heg_layer(xyzzyaaae2))
else
if(periodicity==0.or.use_expot)then
rele(1,xyzzyaaab2)=ranx()*xyzzyaaaw2(1)+xyzzyaaaz2(1)
else
rele(1,xyzzyaaab2)=xyzzyaaaw2(1)*real(xyzzyaaab2-xyzzyaaad2-1,dp)/real&
&(heg_nele(xyzzyaaae2),dp) +xyzzyaaaz2(1)
endif
rele(2,xyzzyaaab2)=heg_ylayer(heg_layer(xyzzyaaae2))
rele(3,xyzzyaaab2)=heg_zlayer(heg_layer(xyzzyaaae2))
endif
if(hard_sphere)then
do xyzzyaaak2=1,xyzzyaaab2-1
if(hard_op_spins.and.which_spin(xyzzyaaak2)==which_spin(xyzzyaaab2))cy&
&cle
xyzzyaabe2(1:3,1)=rele(1:3,xyzzyaaak2)-rele(1:3,xyzzyaaab2)
call minimum_image(4,1,xyzzyaabe2)
if(sum(xyzzyaabe2(1:dimensionality,1)**2)<=hard_diam_sq)cycle fluid_at&
&tempt
enddo
endif
if(expot_is_finite(rele(:,xyzzyaaab2),xyzzyaaae2))exit
enddo fluid_attempt
enddo
endif
enddo
else
call errstop_master('POINTS','There do not seem to be any nuclei and t&
&his is not a model system.  Not sure what to do.')
endif
do xyzzyaaah2=1,dimensionality
where(initial_rele_set)rele(xyzzyaaah2,:)=initial_rele(xyzzyaaah2,:)
enddo
if(on_top_ii/=0)rele(1:dimensionality,on_top_jj)=rele(1:dimensionality&
&,on_top_ii)
if(noncoll_spin)then
if(any(heg_orbtype(:,1)==3))then
do xyzzyaaaj2=1,netot
xyzzyaabf2=cos(0.5d0*sdw_orb_theta(xyzzyaaaj2))**2
xyzzyaabg2=sin(0.5d0*sdw_orb_theta(xyzzyaaaj2))**2
if(xyzzyaabf2>xyzzyaabg2)then
sele(xyzzyaaaj2)=1
else
sele(xyzzyaaaj2)=2
endif
enddo
else
call errwarn('POINTS','experimental POINTS algorithm in use for non-co&
&llinear spins.')
do xyzzyaaaj2=1,netot
if(xyzzyaaaj2>int(dble(netot)*(1.d0-sdw_init_polarization)*0.5d0))then
sele(xyzzyaaaj2)=1
else
sele(xyzzyaaaj2)=2
endif
enddo
endif
endif
if(xyzzyaabi2.and.nitot>1.and.am_master.and..not.ignore_ionic_interact&
&ions)then
if(model_system)then
xyzzyaabj2=.false.
else
xyzzyaabj2=.true.
il : do xyzzyaaaj2=1,nitot
do xyzzyaaak2=1,nitot
if(iontype(xyzzyaaaj2)==iontype(xyzzyaaak2).and.any(neion(xyzzyaaaj2,1&
&:nspin)/=neion(xyzzyaaak2,1:nspin)))then
xyzzyaabj2=.false.
exit il
endif
enddo
enddo il
endif
if(xyzzyaabj2)then
allocate(xyzzyaabk2(0:92))
xyzzyaabk2=.false.
if(nspin==2.and.electron_system)then
call wout('No of up/down spin electrons initially associated with each&
& ion:')
call wout('-----------------------------------------------------------&
&-----')
call wout('Element, # up-spin, # down-spin  ')
else
call wout('No of up/down spin particles initially associated with each&
& ion:')
call wout('-----------------------------------------------------------&
&-----')
call wout('Element, # particles of each type')
endif
do xyzzyaaaj2=1,nitot
xyzzyaaak2=nint(atnum(iontype(xyzzyaaaj2)))
if(.not.xyzzyaabk2(xyzzyaaak2))then
if(nspin==2.and.electron_system)then
if(zion(iontype(xyzzyaaaj2))>1.d-10)then
write(tmpr,'(a2,5x,i8,1x,i10)')periodic_table(xyzzyaaak2),neion(xyzzya&
&aaj2,1:nspin)
else
write(tmpr,'(a4,3x,i8,1x,i10)')'Null',neion(xyzzyaaaj2,1:nspin)
endif
else
if(zion(iontype(xyzzyaaaj2))>1.d-10)then
write(tmpr,'(1x,a2,5x,13(1x,i4))')periodic_table(xyzzyaaak2),neion(xyz&
&zyaaaj2,1:nspin)
else
write(tmpr,'(1x,a4,3x,13(1x,i4))')'Null',neion(xyzzyaaaj2,1:nspin)
endif
endif
call wout(tmpr)
xyzzyaabk2(xyzzyaaak2)=.true.
endif
enddo
deallocate(xyzzyaabk2)
else
if(nitot<=64)then
if(nspin==2.and.electron_system)then
call wout('No of up/down spin electrons initially associated with each&
& ion:')
call wout('-----------------------------------------------------------&
&-----')
if(all(zion(iontype(:))>1.d-10))then
call wout('Ion,    # up-spin, # down-spin')
else
call wout('Site,   # up-spin, # down-spin')
endif
do xyzzyaaaj2=1,nitot
if(xyzzyaaaj2<10)then
write(tmpr,'(i1,5x,i8,1x,i10)')xyzzyaaaj2,neion(xyzzyaaaj2,1:nspin)
elseif(xyzzyaaaj2>=10.and.xyzzyaaaj2<100)then
write(tmpr,'(i2,4x,i8,1x,i10)')xyzzyaaaj2,neion(xyzzyaaaj2,1:nspin)
elseif(xyzzyaaaj2>=100.and.xyzzyaaaj2<1000)then
write(tmpr,'(i3,3x,i8,1x,i10)')xyzzyaaaj2,neion(xyzzyaaaj2,1:nspin)
elseif(xyzzyaaaj2>=1000.and.xyzzyaaaj2<10000)then
write(tmpr,'(i4,2x,i8,1x,i10)')xyzzyaaaj2,neion(xyzzyaaaj2,1:nspin)
else
write(tmpr,'(i6,i10,1x,i10)')xyzzyaaaj2,neion(xyzzyaaaj2,1:nspin)
endif
call wout(tmpr)
enddo
else
call wout('No of up/down spin particles initially associated with each&
& ion:')
call wout('-----------------------------------------------------------&
&-----')
if(all(zion(iontype(:))>1.d-10))then
call wout('Ion,    # particles of each type')
else
call wout('Site,   # particles of each type')
endif
do xyzzyaaaj2=1,nitot
write(tmpr,'(1x,i4,5x,13(1x,i4))')xyzzyaaaj2,neion(xyzzyaaaj2,1:nspin)
call wout(tmpr)
enddo
endif
else
call wordwrap('Suppressing printout of the number of electrons associa&
&ted with each ion, since there are too many ions. Read the input file&
& instead.')
endif
endif
if(.not.init_by_ion.and..not.init_by_iontype)then
if(xyzzyaabj2.or.(.not.xyzzyaabj2.and.nitot<=64))then
call wout('[Use edist_by_ion/edist_by_iontype in input to change this]&
&')
call wout()
endif
endif
call wout()
endif
end subroutine points
subroutine init_geometry
implicit none
integer xyzzyaaaa3,xyzzyaaab3,xyzzyaaac3,xyzzyaaad3,xyzzyaaae3,xyzzyaa&
&af3,xyzzyaaag3,xyzzyaaah3(3),xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3,xyzzyaa&
&al3,xyzzyaaam3,xyzzyaaan3,xyzzyaaao3,xyzzyaaap3(3),xyzzyaaaq3(3)
integer,allocatable :: xyzzyaaar3(:),xyzzyaaas3(:,:)
real(dp) xyzzyaaat3,xyzzyaaau3,xyzzyaaav3,xyzzyaaaw3,xyzzyaaax3,xyzzya&
&aay3,xyzzyaaaz3(3),xyzzyaaba3(3),xyzzyaabb3(3),xyzzyaabc3(3),xyzzyaab&
&d3(3),xyzzyaabe3(3)
real(dp),parameter :: xyzzyaabf3=1.d-12
logical xyzzyaabg3
pamat(1,:)=pa1(:)
pamat(2,:)=pa2(:)
pamat(3,:)=pa3(:)
painv(:,:)=inverse3(pamat(:,:),determinant=pvolume)
pvolume=abs(pvolume)
if(pvolume<1.d-16)call errstop('INIT_GEOMETRY','Primitive lattice vect&
&ors apparently linearly dependent.')
pbmat(:,:)=transpose(painv(:,:))*twopi
pbinv(:,:)=transpose(pamat(:,:))*one_over_twopi
pb1(:)=pbmat(1,:)
pb2(:)=pbmat(2,:)
pb3(:)=pbmat(3,:)
pb11=dot_product(pb1,pb1)
pb22=dot_product(pb2,pb2)
pb33=dot_product(pb3,pb3)
pb12=dot_product(pb1,pb2)
pb13=dot_product(pb1,pb3)
pb23=dot_product(pb2,pb3)
scell_mat_inv=inverse3(dble(scell_matrix),xyzzyaabg3)
if(xyzzyaabg3)call errstop('INIT_GEOMETRY','SCELL_MATRIX is singular.'&
&)
allocate(xyzzyaaas3(3,npcells),stat=xyzzyaaal3)
if(xyzzyaaal3/=0)call errstop('INIT_GEOMETRY','Allocation error: pcell&
&_int.')
do xyzzyaaab3=1,3
xyzzyaaap3(xyzzyaaab3)=0
xyzzyaaaq3(xyzzyaaab3)=0
do xyzzyaaac3=1,3
if(scell_matrix(xyzzyaaac3,xyzzyaaab3)>0)then
xyzzyaaap3(xyzzyaaab3)=xyzzyaaap3(xyzzyaaab3)+scell_matrix(xyzzyaaac3,&
&xyzzyaaab3)
else
xyzzyaaaq3(xyzzyaaab3)=xyzzyaaaq3(xyzzyaaab3)+scell_matrix(xyzzyaaac3,&
&xyzzyaaab3)
endif
enddo
enddo
xyzzyaaao3=0
outer: do xyzzyaaai3=xyzzyaaaq3(1),xyzzyaaap3(1)
xyzzyaabc3=dble(xyzzyaaai3)*scell_mat_inv(1,1:3)
do xyzzyaaaj3=xyzzyaaaq3(2),xyzzyaaap3(2)
xyzzyaabd3=xyzzyaabc3+dble(xyzzyaaaj3)*scell_mat_inv(2,1:3)
do xyzzyaaak3=xyzzyaaaq3(3),xyzzyaaap3(3)
xyzzyaabe3=xyzzyaabd3+dble(xyzzyaaak3)*scell_mat_inv(3,1:3)
if(all(xyzzyaabe3>=-xyzzyaabf3).and.all(xyzzyaabe3<1.d0-xyzzyaabf3))th&
&en
xyzzyaaao3=xyzzyaaao3+1
xyzzyaaas3(1:3,xyzzyaaao3)=(/xyzzyaaai3,xyzzyaaaj3,xyzzyaaak3/)
if(xyzzyaaao3==npcells)exit outer
endif
enddo
enddo
enddo outer
if(xyzzyaaao3/=npcells)call errstop('INIT_GEOMETRY','Have been unable &
&to find all the primitive cells.')
amat=matmul(dble(scell_matrix),pamat)
a1(:)=amat(1,:)
a2(:)=amat(2,:)
a3(:)=amat(3,:)
ainv(:,:)=inverse3(amat,xyzzyaabg3)
if(xyzzyaabg3)call errstop('INIT_GEOMETRY','SCELL_MATRIX is singular.'&
&)
bmat=transpose(ainv)*twopi
binv=transpose(amat)*one_over_twopi
b1(:)=bmat(1,:)
b2(:)=bmat(2,:)
b3(:)=bmat(3,:)
b11=dot_product(b1,b1)
b22=dot_product(b2,b2)
b33=dot_product(b3,b3)
b12=dot_product(b1,b2)
b13=dot_product(b1,b3)
b23=dot_product(b2,b3)
volume=pvolume*real(npcells,dp)
if(periodicity==1)then
if(a1(2)/=0.d0.or.a1(3)/=0.d0)call errstop('INIT_GEOMETRY','For 1D sys&
&tems the first lattice vector should point in the x direction, which &
&should be the direction of translational symmetry.')
if(dot_product(a1,a2)/=0.d0.or.dot_product(a1,a3)/=0.d0)call errstop('&
&INIT_GEOMETRY','For 1D systems the second and third "lattice vectors"&
& should be orthogonal to the first lattice vector.')
endif
if(periodicity==2)then
if(a1(3)/=0.d0.or.a2(3)/=0.d0)call errstop('INIT_GEOMETRY','For 2D sys&
&tems the first two lattice vectors must lie in the (x,y) plane. The s&
&ystem should have translational symmetry in the direction of the firs&
&t two lattice vectors.')
if(dot_product(a1,a3)/=0.d0.or.dot_product(a2,a3)/=0.d0)call errstop('&
&INIT_GEOMETRY','For 2D systems the first and second lattice vectors s&
&hould be orthogonal to the third "lattice vector".')
area=abs(amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1))
endif
if(periodicity==2.or.periodicity==3)then
xyzzyaaat3=sqrt(min(sum(a1(:)**2),sum(a2(:)**2),sum(a3(:)**2)))
xyzzyaaah3(1)=int(xyzzyaaat3*sqrt(sum(b1**2))*one_over_twopi+0.01d0)
xyzzyaaah3(2)=int(xyzzyaaat3*sqrt(sum(b2**2))*one_over_twopi+0.01d0)
xyzzyaaah3(3)=int(xyzzyaaat3*sqrt(sum(b3**2))*one_over_twopi+0.01d0)
xyzzyaaav3=xyzzyaaat3**2
do xyzzyaaai3=1,xyzzyaaah3(1)
do xyzzyaaaj3=-xyzzyaaah3(2),xyzzyaaah3(2)
do xyzzyaaak3=-xyzzyaaah3(3),xyzzyaaah3(3)
xyzzyaaaw3=dble(xyzzyaaai3)
xyzzyaaax3=dble(xyzzyaaaj3)
xyzzyaaay3=dble(xyzzyaaak3)
xyzzyaaau3=(xyzzyaaaw3*a1(1)+xyzzyaaax3*a2(1)+xyzzyaaay3*a3(1))**2+(xy&
&zzyaaaw3*a1(2)+xyzzyaaax3*a2(2)+xyzzyaaay3*a3(2))**2+(xyzzyaaaw3*a1(3&
&)+xyzzyaaax3*a2(3)+xyzzyaaay3*a3(3))**2
if(xyzzyaaau3<xyzzyaaav3)xyzzyaaav3=xyzzyaaau3
enddo
enddo
enddo
do xyzzyaaaj3=1,xyzzyaaah3(2)
do xyzzyaaak3=-xyzzyaaah3(3),xyzzyaaah3(3)
xyzzyaaax3=dble(xyzzyaaaj3)
xyzzyaaay3=dble(xyzzyaaak3)
xyzzyaaau3=(xyzzyaaax3*a2(1)+xyzzyaaay3*a3(1))**2+(xyzzyaaax3*a2(2)+xy&
&zzyaaay3*a3(2))**2+(xyzzyaaax3*a2(3)+xyzzyaaay3*a3(3))**2
if(xyzzyaaau3<xyzzyaaav3)xyzzyaaav3=xyzzyaaau3
enddo
enddo
do xyzzyaaak3=1,xyzzyaaah3(3)
xyzzyaaay3=dble(xyzzyaaak3)
xyzzyaaau3=(xyzzyaaay3*a3(1))**2+(xyzzyaaay3*a3(2))**2+(xyzzyaaay3*a3(&
&3))**2
if(xyzzyaaau3<xyzzyaaav3)xyzzyaaav3=xyzzyaaau3
enddo
wigner_seitz_radius=0.5d0*sqrt(xyzzyaaav3)
elseif(periodicity==1)then
wigner_seitz_radius=0.5d0*a1(1)
endif
nitot=npcells*nbasis
if(nbasis>0)then
allocate(iontype(nitot),rion(3,nitot),xyzzyaaar3(nbasis),rionion(4,nit&
&ot,nitot),iontype_prim(nbasis),ion_prim(nitot),nbasis_to_nitot(xyzzya&
&aaq3(1):xyzzyaaap3(1),xyzzyaaaq3(2):xyzzyaaap3(2),xyzzyaaaq3(3):xyzzy&
&aaap3(3),nbasis),gpcell(3,nitot),stat=xyzzyaaal3)
call check_alloc(xyzzyaaal3,'INIT_GEOMETRY','1')
xyzzyaaar3=0
iontype_prim=0
which_ion_displaced=0
ion_displacement=0.d0
xyzzyaaaa3=0
a:do xyzzyaaab3=1,nbasis
do xyzzyaaad3=1,xyzzyaaab3-1
if(xyzzyaaar3(xyzzyaaad3)==atno(xyzzyaaab3))then
xyzzyaaar3(xyzzyaaab3)=xyzzyaaar3(xyzzyaaad3)
iontype_prim(xyzzyaaab3)=iontype_prim(xyzzyaaad3)
cycle a
endif
enddo
xyzzyaaaa3=xyzzyaaaa3+1
xyzzyaaar3(xyzzyaaab3)=atno(xyzzyaaab3)
iontype_prim(xyzzyaaab3)=xyzzyaaaa3
enddo a
nitype=xyzzyaaaa3
allocate(nion(nitype),stat=xyzzyaaal3)
call check_alloc(xyzzyaaal3,'INIT_GEOMETRY','2')
nion=0
xyzzyaaaa3=1
do xyzzyaaag3=1,npcells
xyzzyaabb3=dble(xyzzyaaas3(1,xyzzyaaag3))*pa1+dble(xyzzyaaas3(2,xyzzya&
&aag3))*pa2+dble(xyzzyaaas3(3,xyzzyaaag3))*pa3
do xyzzyaaaf3=1,nbasis
gpcell(1:3,xyzzyaaaa3)=xyzzyaabb3
rion(1:3,xyzzyaaaa3)=xyzzyaabb3+basis(1:3,xyzzyaaaf3)
xyzzyaaae3=iontype_prim(xyzzyaaaf3)
iontype(xyzzyaaaa3)=xyzzyaaae3
nion(xyzzyaaae3)=nion(xyzzyaaae3)+1
ion_prim(xyzzyaaaa3)=xyzzyaaaf3
xyzzyaaaa3=xyzzyaaaa3+1
enddo
enddo
deallocate(xyzzyaaas3)
rionion(:,:,:)=0.d0
do xyzzyaaaf3=1,nitot
call ei_distances(rion(:,xyzzyaaaf3),nitot,rion,rionion(:,xyzzyaaaf3,:&
&))
enddo
deallocate(xyzzyaaar3)
endif
if(isperiodic)then
call lattice_generator(num_g,periodicity,pa1,pa2,pa3,p_lattice,p_modsq&
&,first_p_in_star,painv,p_lmng,xyzzyaaan3)
call lattice_generator(num_g,periodicity,pb1,pb2,pb3,pr_lattice,pr_mod&
&sq,first_pr_in_star,pbinv,pr_lmng,xyzzyaaan3)
call lattice_generator(num_g,periodicity,a1,a2,a3,s_lattice,s_modsq,fi&
&rst_s_in_star,ainv,s_lmng,xyzzyaaan3)
call lattice_generator(num_g,periodicity,b1,b2,b3,sr_lattice,sr_modsq,&
&first_sr_in_star,binv,sr_lmng,xyzzyaaan3)
if(nbasis>0)then
xyzzyaaaa3=1
xyzzyaaah3=0
xyzzyaaam3=periodicity
do xyzzyaaag3=1,npcells
do xyzzyaaaf3=1,nbasis
xyzzyaaaz3(1:xyzzyaaam3)=matmul(rion(1:xyzzyaaam3,xyzzyaaaa3),ainv(1:x&
&yzzyaaam3,1:xyzzyaaam3))+1.d-5
xyzzyaaba3(1:xyzzyaaam3)=dble(floor(xyzzyaaaz3(1:xyzzyaaam3)))
xyzzyaaaz3(1:xyzzyaaam3)=matmul(xyzzyaaba3(1:xyzzyaaam3),amat(1:xyzzya&
&aam3,1:xyzzyaaam3))
xyzzyaaba3(1:xyzzyaaam3)=rion(1:xyzzyaaam3,xyzzyaaaa3)-xyzzyaaaz3(1:xy&
&zzyaaam3)
xyzzyaaaz3(1:xyzzyaaam3)=matmul(xyzzyaaba3(1:xyzzyaaam3),painv(1:xyzzy&
&aaam3,1:xyzzyaaam3))+1.d-5
xyzzyaaah3(1:xyzzyaaam3)=floor(xyzzyaaaz3(1:xyzzyaaam3))
nbasis_to_nitot(xyzzyaaah3(1),xyzzyaaah3(2),xyzzyaaah3(3),xyzzyaaaf3)=&
&xyzzyaaaa3
xyzzyaaaa3=xyzzyaaaa3+1
enddo
enddo
endif
else
p_lattice(1:3,1)=0.d0
p_lmng(1:3,1)=0
if(trim(atom_basis_type)=='gaussian'.and.molgscreening)call xyzzyaaaa1
endif
end subroutine init_geometry
subroutine netot_nitot_products
implicit none
nitot_netot=nitot*netot
netot_netot=netot*netot
four_nitot_netot=4*nitot*netot
four_netot_netot=4*netot*netot
four_nitot=4*nitot
four_netot=4*netot
three_nitot=3*nitot
three_nitot_netot=3*nitot*netot
three_nemax=3*nemax
six_nemax=6*nemax
four_netot_3=4*netot-3
three_nemax_2=3*nemax-2
three_netot=3*netot
three_nitot_2=3*nitot-2
three_netot_2=3*netot-2
three_nitot_netot_2=3*nitot*netot-2
nine_netot=9*netot
nine_nitot=9*nitot
nine_netot_8=9*netot-8
size_det=nspin*ndet
three_netot_netot=3*netot*netot
nine_netot_netot=9*netot*netot
nemax2_nspin_ndet=nemax*nemax*nspin*ndet
size_dbar=nemax2_nspin_ndet*real1_complex2
size_onedbar=nemax*nemax*real1_complex2
size_fidet=three_netot*real1_complex2
size_prod_lapdet=size_det*nemax*real1_complex2
size_fi_prod_det=3*size_prod_lapdet
nspin_ndet=nspin*ndet
size_rpsi=nemax*ndet*real1_complex2
size_grad=3*size_rpsi
size_sderivs=6*size_rpsi
size_dsmat=3*size_dbar
size_d2smat=6*size_dbar
size_farray=3*netot*real1_complex2
size_harray=9*netot*netot*real1_complex2
sqrt_netot=sqrt(dble(netot))
end subroutine netot_nitot_products
subroutine ei_distances(rvec,ni,rion,eivecs)
use slaarnabt,only : ddot
!$ use openmp_base, only : nitot_min,get_omp_index_range
implicit none
integer,intent(in) :: ni
real(dp),intent(in) :: rvec(3),rion(3,ni)
real(dp),intent(out) :: eivecs(4,ni)
integer xyzzyaaaa5
!$ integer xyzzyaaab5,xyzzyaaac5
call timer('EI_DISTANCES',.true.)
if(isperiodic)then
!$omp parallel default(none) shared(ni,eivecs,rvec,rion,dimensionality&
!$omp &) private(xyzzyaaaa5,xyzzyaaab5,xyzzyaaac5) if(ni>=nitot_min)
!$omp  do
do xyzzyaaaa5=1,ni
eivecs(1,xyzzyaaaa5)=rvec(1)-rion(1,xyzzyaaaa5)
eivecs(2,xyzzyaaaa5)=rvec(2)-rion(2,xyzzyaaaa5)
eivecs(3,xyzzyaaaa5)=rvec(3)-rion(3,xyzzyaaaa5)
enddo
!$omp  enddo
!$  if(ni>=nitot_min)then
!$   xyzzyaaab5=1
!$ xyzzyaaac5=ni
!$   call get_omp_index_range(xyzzyaaab5,xyzzyaaac5)
!$   if(xyzzyaaac5>=xyzzyaaab5)call minimum_image(4,xyzzyaaac5-xyzzyaa&
!$ &ab5+1,eivecs(1,xyzzyaaab5))
!$omp  barrier
!$  else
call minimum_image(4,ni,eivecs)
!$  endif
!$omp  do
do xyzzyaaaa5=1,ni
eivecs(4,xyzzyaaaa5)=sqrt(ddot(dimensionality,eivecs(1,xyzzyaaaa5),1,e&
&ivecs(1,xyzzyaaaa5),1))
enddo
!$omp  enddo
!$omp end parallel
else
!$omp parallel do default(none) shared(ni,eivecs,rvec,rion,dimensional&
!$omp &ity) private(xyzzyaaaa5) if(ni>=nitot_min)
do xyzzyaaaa5=1,ni
eivecs(1,xyzzyaaaa5)=rvec(1)-rion(1,xyzzyaaaa5)
eivecs(2,xyzzyaaaa5)=rvec(2)-rion(2,xyzzyaaaa5)
eivecs(3,xyzzyaaaa5)=rvec(3)-rion(3,xyzzyaaaa5)
eivecs(4,xyzzyaaaa5)=sqrt(ddot(dimensionality,eivecs(1,xyzzyaaaa5),1,e&
&ivecs(1,xyzzyaaaa5),1))
enddo
!$omp end parallel do
endif
call timer('EI_DISTANCES',.false.)
end subroutine ei_distances
subroutine ei_distances_all(ne,rele,ni,rion,eivecs)
use slaarnabt,only : ddot
!$ use openmp_base, only : ne_min_ee
implicit none
integer,intent(in) :: ne,ni
real(dp),intent(in) :: rele(3,ne),rion(3,ni)
real(dp),intent(out) :: eivecs(4,ni,ne)
integer xyzzyaaaa6,xyzzyaaab6
call timer('EI_DISTANCES_ALL',.true.)
if(isperiodic)then
!$omp parallel default(none) shared(netot,ni,eivecs,rele,rion,dimensio&
!$omp &nality) private(xyzzyaaaa6,xyzzyaaab6) if(netot>=ne_min_ee)
!$omp  do
do xyzzyaaaa6=1,netot
do xyzzyaaab6=1,ni
eivecs(1,xyzzyaaab6,xyzzyaaaa6)=rele(1,xyzzyaaaa6)-rion(1,xyzzyaaab6)
eivecs(2,xyzzyaaab6,xyzzyaaaa6)=rele(2,xyzzyaaaa6)-rion(2,xyzzyaaab6)
eivecs(3,xyzzyaaab6,xyzzyaaaa6)=rele(3,xyzzyaaaa6)-rion(3,xyzzyaaab6)
enddo
call minimum_image(4,ni,eivecs(1,1,xyzzyaaaa6))
do xyzzyaaab6=1,ni
eivecs(4,xyzzyaaab6,xyzzyaaaa6)=sqrt(ddot(dimensionality,eivecs(1,xyzz&
&yaaab6,xyzzyaaaa6),1,eivecs(1,xyzzyaaab6,xyzzyaaaa6),1))
enddo
enddo
!$omp  enddo
!$omp end parallel
else
!$omp parallel do default(none) shared(netot,ni,eivecs,rele,rion,dimen&
!$omp &sionality) private(xyzzyaaaa6,xyzzyaaab6) if(netot>ne_min_ee)
do xyzzyaaaa6=1,netot
do xyzzyaaab6=1,ni
eivecs(1,xyzzyaaab6,xyzzyaaaa6)=rele(1,xyzzyaaaa6)-rion(1,xyzzyaaab6)
eivecs(2,xyzzyaaab6,xyzzyaaaa6)=rele(2,xyzzyaaaa6)-rion(2,xyzzyaaab6)
eivecs(3,xyzzyaaab6,xyzzyaaaa6)=rele(3,xyzzyaaaa6)-rion(3,xyzzyaaab6)
eivecs(4,xyzzyaaab6,xyzzyaaaa6)=sqrt(ddot(dimensionality,eivecs(1,xyzz&
&yaaab6,xyzzyaaaa6),1,eivecs(1,xyzzyaaab6,xyzzyaaaa6),1))
enddo
enddo
!$omp end parallel do
endif
call timer('EI_DISTANCES_ALL',.false.)
end subroutine ei_distances_all
subroutine ee_distances(ne,rvec,rele,eevecs1)
!$ use openmp_base, only : ne_min_ee,get_omp_index_range
use slaarnabt,       only : ddot
implicit none
integer,intent(in) :: ne
real(dp),intent(in) :: rvec(3),rele(3,ne)
real(dp),intent(out) :: eevecs1(4,ne)
integer xyzzyaaaa7
!$ integer xyzzyaaab7,xyzzyaaac7
call timer('EE_DISTANCES',.true.)
if(isperiodic)then
!$omp parallel default(none) shared(ne,eevecs1,rvec,rele,dimensionalit&
!$omp &y) private(xyzzyaaaa7,xyzzyaaab7,xyzzyaaac7) if(ne>=ne_min_ee)
!$omp  do
do xyzzyaaaa7=1,ne
eevecs1(1,xyzzyaaaa7)=rvec(1)-rele(1,xyzzyaaaa7)
eevecs1(2,xyzzyaaaa7)=rvec(2)-rele(2,xyzzyaaaa7)
eevecs1(3,xyzzyaaaa7)=rvec(3)-rele(3,xyzzyaaaa7)
enddo
!$omp  enddo
!$  if(ne>=ne_min_ee)then
!$   xyzzyaaab7=1
!$ xyzzyaaac7=ne
!$   call get_omp_index_range(xyzzyaaab7,xyzzyaaac7)
!$   if(xyzzyaaac7>=xyzzyaaab7)call minimum_image(4,xyzzyaaac7-xyzzyaa&
!$ &ab7+1,eevecs1(1,xyzzyaaab7))
!$omp  barrier
!$  else
call minimum_image(4,ne,eevecs1)
!$  endif
!$omp  do
do xyzzyaaaa7=1,ne
eevecs1(4,xyzzyaaaa7)=sqrt(ddot(dimensionality,eevecs1(1,xyzzyaaaa7),1&
&,eevecs1(1,xyzzyaaaa7),1))
enddo
!$omp  enddo
!$omp end parallel
else
!$omp parallel do default(none) shared(ne,eevecs1,rvec,rele,dimensiona&
!$omp &lity) private(xyzzyaaaa7) if(ne>=ne_min_ee)
do xyzzyaaaa7=1,ne
eevecs1(1,xyzzyaaaa7)=rvec(1)-rele(1,xyzzyaaaa7)
eevecs1(2,xyzzyaaaa7)=rvec(2)-rele(2,xyzzyaaaa7)
eevecs1(3,xyzzyaaaa7)=rvec(3)-rele(3,xyzzyaaaa7)
eevecs1(4,xyzzyaaaa7)=sqrt(ddot(dimensionality,eevecs1(1,xyzzyaaaa7),1&
&,eevecs1(1,xyzzyaaaa7),1))
enddo
!$omp end parallel do
endif
call timer('EE_DISTANCES',.false.)
end subroutine ee_distances
subroutine ee_distances_all(ne,rele,eevecs)
use slaarnabt,only : ddot
!$ use openmp_base, only : ne_min_ee
implicit none
integer,intent(in) :: ne
real(dp),intent(in) :: rele(3,ne)
real(dp),intent(out) :: eevecs(4,ne,ne)
integer xyzzyaaaa8,xyzzyaaab8
real(dp) xyzzyaaac8(3),xyzzyaaad8,xyzzyaaae8(3)
call timer('EE_DISTANCES_ALL',.true.)
if(isperiodic)then
!$omp parallel default(none) shared(ne,rele,eevecs,dimensionality) pri&
!$omp &vate(xyzzyaaac8,xyzzyaaae8,xyzzyaaad8,xyzzyaaaa8,xyzzyaaab8) if&
!$omp &(ne>=ne_min_ee)
!$omp  do
do xyzzyaaaa8=1,ne-1
xyzzyaaac8(1:3)=rele(1:3,xyzzyaaaa8)
eevecs(1:4,xyzzyaaaa8,xyzzyaaaa8)=0.d0
do xyzzyaaab8=xyzzyaaaa8+1,ne
eevecs(1,xyzzyaaab8,xyzzyaaaa8)=xyzzyaaac8(1)-rele(1,xyzzyaaab8)
eevecs(2,xyzzyaaab8,xyzzyaaaa8)=xyzzyaaac8(2)-rele(2,xyzzyaaab8)
eevecs(3,xyzzyaaab8,xyzzyaaaa8)=xyzzyaaac8(3)-rele(3,xyzzyaaab8)
enddo
call minimum_image(4,ne-xyzzyaaaa8,eevecs(1,xyzzyaaaa8+1,xyzzyaaaa8))
do xyzzyaaab8=xyzzyaaaa8+1,ne
xyzzyaaae8(1:3)=eevecs(1:3,xyzzyaaab8,xyzzyaaaa8)
xyzzyaaad8=sqrt(ddot(dimensionality,xyzzyaaae8(1),1,xyzzyaaae8(1),1))
eevecs(4,xyzzyaaab8,xyzzyaaaa8)=xyzzyaaad8
eevecs(1,xyzzyaaaa8,xyzzyaaab8)=-xyzzyaaae8(1)
eevecs(2,xyzzyaaaa8,xyzzyaaab8)=-xyzzyaaae8(2)
eevecs(3,xyzzyaaaa8,xyzzyaaab8)=-xyzzyaaae8(3)
eevecs(4,xyzzyaaaa8,xyzzyaaab8)=xyzzyaaad8
enddo
enddo
!$omp  enddo
!$omp end parallel
eevecs(1:4,ne,ne)=0.d0
else
!$omp parallel do default(none) shared(ne,rele,eevecs,dimensionality) &
!$omp &private(xyzzyaaac8,xyzzyaaaa8,xyzzyaaab8,xyzzyaaae8,xyzzyaaad8)&
!$omp & if(ne>=ne_min_ee)
do xyzzyaaaa8=1,ne-1
xyzzyaaac8(1:3)=rele(1:3,xyzzyaaaa8)
eevecs(1:4,xyzzyaaaa8,xyzzyaaaa8)=0.d0
do xyzzyaaab8=xyzzyaaaa8+1,ne
xyzzyaaae8(1)=xyzzyaaac8(1)-rele(1,xyzzyaaab8)
xyzzyaaae8(2)=xyzzyaaac8(2)-rele(2,xyzzyaaab8)
xyzzyaaae8(3)=xyzzyaaac8(3)-rele(3,xyzzyaaab8)
xyzzyaaad8=sqrt(ddot(dimensionality,xyzzyaaae8(1),1,xyzzyaaae8(1),1))
eevecs(1,xyzzyaaab8,xyzzyaaaa8)=xyzzyaaae8(1)
eevecs(2,xyzzyaaab8,xyzzyaaaa8)=xyzzyaaae8(2)
eevecs(3,xyzzyaaab8,xyzzyaaaa8)=xyzzyaaae8(3)
eevecs(4,xyzzyaaab8,xyzzyaaaa8)=xyzzyaaad8
eevecs(1,xyzzyaaaa8,xyzzyaaab8)=-xyzzyaaae8(1)
eevecs(2,xyzzyaaaa8,xyzzyaaab8)=-xyzzyaaae8(2)
eevecs(3,xyzzyaaaa8,xyzzyaaab8)=-xyzzyaaae8(3)
eevecs(4,xyzzyaaaa8,xyzzyaaab8)=xyzzyaaad8
enddo
enddo
!$omp end parallel do
eevecs(1:4,ne,ne)=0.d0
endif
call timer('EE_DISTANCES_ALL',.false.)
end subroutine ee_distances_all
subroutine vector_difference(rvec1,rvec2,rdiff)
implicit none
real(dp),intent(in) :: rvec1(3),rvec2(3)
real(dp),intent(out) :: rdiff(3)
real(dp) xyzzyaaaa9(4)
xyzzyaaaa9(1:3)=rvec2(1:3)-rvec1(1:3)
xyzzyaaaa9(4)=0.d0
if(isperiodic)call minimum_image(4,1,xyzzyaaaa9)
rdiff(1:3)=xyzzyaaaa9(1:3)
end subroutine vector_difference
real(dp) function ee_kato_gamma(ispin,jspin,indist_in,f_zero)
implicit none
integer,intent(in) :: ispin,jspin
logical,intent(in),optional :: indist_in
character(8),intent(out),optional :: f_zero
real(dp) xyzzyaaaa10,xyzzyaaab10
logical xyzzyaaac10
xyzzyaaac10=.true.
if(present(indist_in))xyzzyaaac10=indist_in
xyzzyaaac10=ispin==jspin.and.xyzzyaaac10
ee_kato_gamma=0.d0
if(present(f_zero))f_zero='r'
select case(trim(interaction))
case('none')
continue
case('coulomb','ewald','mpc','ewald_mpc','mpc_ewald')
if(heg_layer(ispin)/=heg_layer(jspin))then
if(present(f_zero))then
f_zero='quasi'
else
return
endif
endif
if(ispin==jspin)then
if(pinfmass(ispin))then
xyzzyaaab10=0.d0
else
xyzzyaaab10=pmass(ispin)
endif
else
if(pinfmass(ispin).and.pinfmass(jspin))then
xyzzyaaab10=0.d0
elseif(pinfmass(ispin))then
xyzzyaaab10=2*pmass(jspin)
elseif(pinfmass(jspin))then
xyzzyaaab10=2*pmass(ispin)
else
xyzzyaaab10=2.d0*pmass(ispin)*pmass(jspin)/(pmass(ispin)+pmass(jspin))
endif
endif
xyzzyaaaa10=0.d0
if(xyzzyaaac10)then
select case(dimensionality)
case(3)
xyzzyaaaa10=0.25d0
case(2)
xyzzyaaaa10=third
case(1)
xyzzyaaaa10=0.5d0
end select
else
select case(dimensionality)
case(3)
xyzzyaaaa10=0.5d0
case(2)
xyzzyaaaa10=1.d0
case(1)
xyzzyaaaa10=0.5d0
end select
endif
ee_kato_gamma=xyzzyaaaa10*xyzzyaaab10*pcharge(ispin)*pcharge(jspin)
case('manual')
select case(trim(int_name))
case('square_well')
continue
case('poschl_teller')
continue
case('hard_sphere')
continue
case('polynomial')
continue
case('logarithmic','2D_int')
if(present(f_zero))then
f_zero='r2log'
else
ee_kato_gamma=0.d0
return
endif
if(ispin==jspin)then
if(pinfmass(ispin))then
xyzzyaaab10=0.d0
else
xyzzyaaab10=pmass(ispin)
endif
else
if(pinfmass(ispin).and.pinfmass(jspin))then
xyzzyaaab10=0.d0
elseif(pinfmass(ispin))then
xyzzyaaab10=2*pmass(jspin)
elseif(pinfmass(jspin))then
xyzzyaaab10=2*pmass(ispin)
else
xyzzyaaab10=2.d0*pmass(ispin)*pmass(jspin)/(pmass(ispin)+pmass(jspin))
endif
endif
xyzzyaaaa10=0.d0
if(xyzzyaaac10)then
select case(dimensionality)
case(3)
xyzzyaaaa10=0.2d0
case(2)
xyzzyaaaa10=0.25d0
case(1)
xyzzyaaaa10=third
end select
else
select case(dimensionality)
case(3)
xyzzyaaaa10=third
case(2)
xyzzyaaaa10=0.5d0
case(1)
xyzzyaaaa10=1.d0
end select
endif
ee_kato_gamma=-xyzzyaaaa10*xyzzyaaab10*pcharge(ispin)*pcharge(jspin)/(&
&2.d0*man_int_params(1))
case('dipole')
if(present(f_zero))then
f_zero='1/sqrtr'
else
ee_kato_gamma=0.d0
return
endif
if(ispin==jspin)then
if(pinfmass(ispin))then
xyzzyaaab10=0.d0
else
xyzzyaaab10=pmass(ispin)
endif
else
if(pinfmass(ispin).and.pinfmass(jspin))then
xyzzyaaab10=0.d0
elseif(pinfmass(ispin))then
xyzzyaaab10=2*pmass(jspin)
elseif(pinfmass(jspin))then
xyzzyaaab10=2*pmass(ispin)
else
xyzzyaaab10=2.d0*pmass(ispin)*pmass(jspin)/(pmass(ispin)+pmass(jspin))
endif
endif
ee_kato_gamma=-sqrt(xyzzyaaab10*abs(man_int_params(1)))
case('2D_tilted_dipole')
if(present(f_zero))then
f_zero='1/sqrtr'
else
ee_kato_gamma=0.d0
return
endif
if(ispin==jspin)then
if(pinfmass(ispin))then
xyzzyaaab10=0.d0
else
xyzzyaaab10=pmass(ispin)
endif
else
if(pinfmass(ispin).and.pinfmass(jspin))then
xyzzyaaab10=0.d0
elseif(pinfmass(ispin))then
xyzzyaaab10=2*pmass(jspin)
elseif(pinfmass(jspin))then
xyzzyaaab10=2*pmass(ispin)
else
xyzzyaaab10=2.d0*pmass(ispin)*pmass(jspin)/(pmass(ispin)+pmass(jspin))
endif
endif
if(man_int_params(2)>0.95531661812450927817d0)then
call errstop('EE_KATO_GAMMA','2D_tilted_dipole: theta must be less tha&
&n 1/2 arccos(-1/sqrt(3))')
else
ee_kato_gamma=-0.5d0*sqrt(xyzzyaaab10*abs(man_int_params(1))*(4.d0-6.d&
&0*sin(man_int_params(2))*sin(man_int_params(2))))
endif
case('pseudodipole')
continue
case('2D_tilt_pseudodipole')
continue
case default
continue
end select
case default
continue
end select
end function ee_kato_gamma
real(dp) function en_kato_gamma(ispin,z)
implicit none
integer,intent(in) :: ispin
real(dp),intent(in) :: z
real(dp) xyzzyaaaa11,xyzzyaaab11
xyzzyaaab11=2*pmass(ispin)
xyzzyaaaa11=0.d0
select case(dimensionality)
case(3)
xyzzyaaaa11=0.5d0
case(2)
xyzzyaaaa11=1.d0
case(1)
xyzzyaaaa11=0.5d0
end select
en_kato_gamma=xyzzyaaaa11*xyzzyaaab11*pcharge(ispin)*z
end function en_kato_gamma
subroutine lattice_generator(ngvec,periodicity,v1,v2,v3,lattice_point,&
&squared_length,first_in_star,inverse_latvec,lmng,nstar)
implicit none
integer,intent(in) :: ngvec,periodicity
integer,intent(out) :: lmng(3,ngvec),first_in_star((ngvec+3)/2),nstar
real(dp),intent(in) :: v1(3),v2(3),v3(3)
real(dp),intent(out) :: squared_length((ngvec+3)/2),lattice_point(3,ng&
&vec),inverse_latvec(3,3)
integer xyzzyaaaa12,xyzzyaaab12,xyzzyaaac12,xyzzyaaad12,xyzzyaaae12,xy&
&zzyaaaf12,xyzzyaaag12,xyzzyaaah12,xyzzyaaai12,xyzzyaaaj12,xyzzyaaak12&
&,xyzzyaaal12,xyzzyaaam12,xyzzyaaan12,xyzzyaaao12,xyzzyaaap12,xyzzyaaa&
&q12,xyzzyaaar12,xyzzyaaas12,xyzzyaaat12,xyzzyaaau12,xyzzyaaav12
integer,allocatable :: xyzzyaaaw12(:),xyzzyaaax12(:),xyzzyaaay12(:),xy&
&zzyaaaz12(:),xyzzyaaba12(:),xyzzyaabb12(:)
real(dp) xyzzyaabc12,xyzzyaabd12,xyzzyaabe12,xyzzyaabf12,r,xyzzyaabg12&
&,xyzzyaabh12(3,3)
real(dp),parameter :: xyzzyaabi12=1.d-6
real(dp),allocatable :: xyzzyaabj12(:),xyzzyaabk12(:),xyzzyaabl12(:),x&
&yzzyaabm12(:)
logical xyzzyaabn12
logical,parameter :: xyzzyaabo12=.false.
character(80) tmpr
if(mod(ngvec,2)/=1)call errstop('LATTICE_GENERATOR','ngvec should be o&
&dd.')
xyzzyaaar12=(ngvec+1)/2
xyzzyaaas12=ngvec*2
xyzzyaaat12=xyzzyaaas12+xyzzyaaas12+2
xyzzyaaau12=3*xyzzyaaas12
allocate(xyzzyaabb12(ngvec),xyzzyaaaw12(xyzzyaaas12),xyzzyaaax12(xyzzy&
&aaas12),xyzzyaaay12(xyzzyaaau12),xyzzyaaaz12(xyzzyaaas12),xyzzyaaba12&
&(xyzzyaaas12),xyzzyaabj12(xyzzyaaau12),xyzzyaabk12(xyzzyaaas12),xyzzy&
&aabl12(xyzzyaaas12),xyzzyaabm12(xyzzyaaau12),stat=xyzzyaaav12)
call check_alloc(xyzzyaaav12,'LATTICE_GENERATOR','')
xyzzyaabh12(1,1:3)=v1(1:3)
xyzzyaabh12(2,1:3)=v2(1:3)
xyzzyaabh12(3,1:3)=v3(1:3)
lattice_point(1:3,1)=0.d0
squared_length(1)=0.d0
lmng(1:3,1)=0
xyzzyaabb12(1)=1
first_in_star(1)=1
nstar=1
ngpri=1
xyzzyaabe12=0.d0
xyzzyaabd12=0.d0
xyzzyaabc12=0.d0
inverse_latvec=inverse3(xyzzyaabh12,xyzzyaabn12)
if(xyzzyaabn12)call errstop('LATTICE_GENERATOR','Lattice vectors are l&
&inearly dependent.')
xyzzyaabg12=0.d0
select case (periodicity)
case(0)
return
case(1)
xyzzyaabc12=sqrt(inverse_latvec(1,1)**2+inverse_latvec(2,1)**2+inverse&
&_latvec(3,1)**2)
case(2)
xyzzyaabd12=sqrt(inverse_latvec(1,2)**2+inverse_latvec(2,2)**2+inverse&
&_latvec(3,2)**2)
xyzzyaabg12=max(xyzzyaabg12,xyzzyaabd12)
xyzzyaabc12=sqrt(inverse_latvec(1,1)**2+inverse_latvec(2,1)**2+inverse&
&_latvec(3,1)**2)
case(3)
xyzzyaabe12=sqrt(inverse_latvec(1,3)**2+inverse_latvec(2,3)**2+inverse&
&_latvec(3,3)**2)
xyzzyaabg12=xyzzyaabe12
xyzzyaabd12=sqrt(inverse_latvec(1,2)**2+inverse_latvec(2,2)**2+inverse&
&_latvec(3,2)**2)
xyzzyaabg12=max(xyzzyaabg12,xyzzyaabd12)
xyzzyaabc12=sqrt(inverse_latvec(1,1)**2+inverse_latvec(2,1)**2+inverse&
&_latvec(3,1)**2)
case default
call errstop('LATTICE_GENERATOR','Subroutine has been called with peri&
&odicity '//trim(i2s(periodicity))//'.')
end select
xyzzyaabg12=max(xyzzyaabg12,xyzzyaabc12)
xyzzyaabf12=1.d0/dble(periodicity)
xyzzyaaan12=(int(dble(xyzzyaaas12)**xyzzyaabf12)-1)/2
r=dble(xyzzyaaan12)/xyzzyaabg12
xyzzyaabg12=r*(1.d0+xyzzyaabi12)
do
do
xyzzyaaaa12=int(xyzzyaabg12*xyzzyaabc12)
xyzzyaaab12=int(xyzzyaabg12*xyzzyaabd12)
xyzzyaaac12=int(xyzzyaabg12*xyzzyaabe12)
if(((xyzzyaaaa12+xyzzyaaaa12+1)*(xyzzyaaab12+xyzzyaaab12+1)*(xyzzyaaac&
&12+xyzzyaaac12+1))>=xyzzyaaat12)exit
xyzzyaaad12=xyzzyaaaa12
xyzzyaaae12=xyzzyaaab12
xyzzyaaaf12=xyzzyaaac12
xyzzyaabg12=xyzzyaabg12+r
enddo
r=r*0.5d0
xyzzyaabg12=xyzzyaabg12-r
if(r<0.2d0)exit
enddo
xyzzyaaab12=-xyzzyaaae12
xyzzyaaac12=-xyzzyaaaf12
xyzzyaaal12=0
xyzzyaaaj12=xyzzyaaae12
do xyzzyaaag12=-xyzzyaaad12,0
if(xyzzyaaag12==0)xyzzyaaaj12=0
xyzzyaaak12=xyzzyaaaf12
do xyzzyaaah12=xyzzyaaab12,xyzzyaaaj12
if(abs(xyzzyaaah12)==xyzzyaaag12)xyzzyaaak12=-1
do xyzzyaaai12=xyzzyaaac12,xyzzyaaak12
xyzzyaaal12=xyzzyaaal12+3
xyzzyaaay12(xyzzyaaal12-2)=xyzzyaaag12
xyzzyaaay12(xyzzyaaal12-1)=xyzzyaaah12
xyzzyaaay12(xyzzyaaal12)=xyzzyaaai12
enddo
enddo
enddo
do xyzzyaaap12=1,xyzzyaaal12
xyzzyaabm12(xyzzyaaap12)=xyzzyaaay12(xyzzyaaap12)
xyzzyaabj12(xyzzyaaap12)=0.d0
enddo
xyzzyaaal12=xyzzyaaal12/3
call mxmb(xyzzyaabm12,size(xyzzyaabm12),3,1,xyzzyaabh12,1,3,xyzzyaabj1&
&2,3,1,xyzzyaaal12,periodicity,periodicity)
xyzzyaaaq12=0
do xyzzyaaap12=1,xyzzyaaal12
xyzzyaaaq12=xyzzyaaaq12+3
xyzzyaabm12(xyzzyaaap12)=xyzzyaabj12(xyzzyaaaq12-2)**2+xyzzyaabj12(xyz&
&zyaaaq12-1)**2+xyzzyaabj12(xyzzyaaaq12)**2
enddo
xyzzyaaaj12=min(xyzzyaaal12,xyzzyaaar12-1)
call gvsort(xyzzyaabm12,xyzzyaaaz12,xyzzyaabk12,xyzzyaabl12,xyzzyaaaw1&
&2,xyzzyaaax12,xyzzyaaba12,xyzzyaaal12,xyzzyaaaj12,xyzzyaabm12(xyzzyaa&
&al12)*xyzzyaabi12)
r=0.d0
do xyzzyaaaq12=1,xyzzyaaaj12
xyzzyaaao12=xyzzyaaaz12(xyzzyaaaq12)
xyzzyaabf12=xyzzyaabm12(xyzzyaaao12)
xyzzyaaal12=ngpri+1
ngpri=ngpri+2
xyzzyaabb12(xyzzyaaal12)=ngpri
xyzzyaabb12(ngpri)=xyzzyaaal12
if(xyzzyaabf12>=r)then
nstar=nstar+1
first_in_star(nstar)=xyzzyaaal12
squared_length(nstar)=xyzzyaabf12
r=xyzzyaabf12*(1.d0+xyzzyaabi12)
endif
xyzzyaaao12=xyzzyaaao12*3-3
do xyzzyaaan12=1,3
xyzzyaaao12=xyzzyaaao12+1
lmng(xyzzyaaan12,xyzzyaaal12)=xyzzyaaay12(xyzzyaaao12)
lmng(xyzzyaaan12,ngpri)=-xyzzyaaay12(xyzzyaaao12)
lattice_point(xyzzyaaan12,xyzzyaaal12)=xyzzyaabj12(xyzzyaaao12)
lattice_point(xyzzyaaan12,ngpri)=-xyzzyaabj12(xyzzyaaao12)
enddo
enddo
xyzzyaaaj12=nstar+1
first_in_star(xyzzyaaaj12)=ngpri+1
squared_length(xyzzyaaaj12)=1.d25
if(xyzzyaabo12)then
call wout()
call wout('Lattice generator verbose output:')
call wout('  Max lattice vector indices : [ '//trim(i2s(xyzzyaaad12))/&
&/', '//trim(i2s(xyzzyaaae12))//', '//trim(i2s(xyzzyaaaf12))//' ]')
call wout('  No. of vectors created     : '//trim(i2s(ngpri)))
call wout('  No. of stars               : '//trim(i2s(nstar)))
call wout('  Rmax                       : ',sqrt(squared_length(nstar)&
&))
xyzzyaaaj12=nstar
do xyzzyaaam12=1,xyzzyaaaj12
xyzzyaaal12=first_in_star(xyzzyaaam12)
xyzzyaaao12=first_in_star(xyzzyaaam12+1)
xyzzyaaaq12=xyzzyaaao12-xyzzyaaal12
call wout()
call wout('  Star #'//trim(i2s(xyzzyaaam12))//' has '//trim(i2s(xyzzya&
&aaq12))//' vectors, r=',sqrt(squared_length(xyzzyaaam12)))
call wout('   lg i(lg) lx   ly   lz         x           y          z')
do xyzzyaaaq12=xyzzyaaal12,xyzzyaaao12-1
write(tmpr,'(5i5,4x,3f12.6,2i6)')xyzzyaaaq12,xyzzyaabb12(xyzzyaaaq12),&
&lmng(1:3,xyzzyaaaq12),lattice_point(1:3,xyzzyaaaq12)
call wout(tmpr)
enddo
enddo
endif
deallocate(xyzzyaabb12,xyzzyaaaw12,xyzzyaaax12,xyzzyaaay12,xyzzyaaaz12&
&,xyzzyaaba12,xyzzyaabj12,xyzzyaabk12,xyzzyaabl12,xyzzyaabm12)
if(periodicity==1)then
inverse_latvec(1,1)=1.d0/xyzzyaabh12(1,1)
elseif(periodicity==2)then
xyzzyaabg12=1.d0/(v1(1)*v2(2)-v1(2)*v2(1))
inverse_latvec(1,1)=v2(2)*xyzzyaabg12
inverse_latvec(2,1)=-v2(1)*xyzzyaabg12
inverse_latvec(1,2)=-v1(2)*xyzzyaabg12
inverse_latvec(2,2)=v1(1)*xyzzyaabg12
endif
end subroutine lattice_generator
subroutine xyzzyaaaa1
implicit none
integer xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13,xyzzyaaad13,xyzzyaaae13
real(dp) xyzzyaaaf13,xyzzyaaag13,xyzzyaaah13,xyzzyaaai13,xyzzyaaaj13,x&
&yzzyaaak13,xyzzyaaal13,xyzzyaaam13,xyzzyaaan13,xyzzyaaao13,xyzzyaaap1&
&3,xyzzyaaaq13
real(dp),parameter :: xyzzyaaar13=1.d0/ngsgrid_minus_1,xyzzyaaas13=4.d&
&0
logical xyzzyaaat13
allocate(cell_map(0:ngsgrid_minus_1,0:ngsgrid_minus_1,0:ngsgrid_minus_&
&1),gsgrid(3,ngsgrid3),stat=xyzzyaaae13)
call check_alloc(xyzzyaaae13,'MOLGSCREEN_LATTICE','')
xyzzyaaaf13=minval(rion(1,:))-xyzzyaaas13
xyzzyaaag13=maxval(rion(1,:))+xyzzyaaas13
xyzzyaaah13=minval(rion(2,:))-xyzzyaaas13
xyzzyaaai13=maxval(rion(2,:))+xyzzyaaas13
xyzzyaaaj13=minval(rion(3,:))-xyzzyaaas13
xyzzyaaak13=maxval(rion(3,:))+xyzzyaaas13
a1(1)=(xyzzyaaag13-xyzzyaaaf13)
a2(2)=(xyzzyaaai13-xyzzyaaah13)
a3(3)=(xyzzyaaak13-xyzzyaaaj13)
xyzzyaaao13=a1(1)*xyzzyaaar13
pa1(1)=xyzzyaaao13
pa1(2)=0.d0
pa1(3)=0.d0
xyzzyaaap13=a2(2)*xyzzyaaar13
pa2(1)=0.d0
pa2(2)=xyzzyaaap13
pa2(3)=0.d0
xyzzyaaaq13=a3(3)*xyzzyaaar13
pa3(1)=0.d0
pa3(2)=0.d0
pa3(3)=xyzzyaaaq13
pamat(1,1:3)=pa1
pamat(2,1:3)=pa2
pamat(3,1:3)=pa3
painv=inverse3(pamat,xyzzyaaat13)
if(xyzzyaaat13)then
call wout('pa1: ',pa1)
call wout('pa2: ',pa2)
call wout('pa3: ',pa3)
call errstop('MOLGSCREEN_LATTICE ','Lattice vector matrix singular.')
endif
gs_origin(1)=xyzzyaaaf13
gs_origin(2)=xyzzyaaah13
gs_origin(3)=xyzzyaaaj13
xyzzyaaad13=0
do xyzzyaaac13=0,ngsgrid_minus_1
do xyzzyaaab13=0,ngsgrid_minus_1
do xyzzyaaaa13=0,ngsgrid_minus_1
xyzzyaaad13=xyzzyaaad13+1
cell_map(xyzzyaaaa13,xyzzyaaab13,xyzzyaaac13)=xyzzyaaad13
xyzzyaaal13=dble(xyzzyaaaa13)*xyzzyaaar13
gsgrid(1,xyzzyaaad13)=xyzzyaaal13*a1(1)+gs_origin(1)
xyzzyaaam13=dble(xyzzyaaab13)*xyzzyaaar13
gsgrid(2,xyzzyaaad13)=xyzzyaaam13*a2(2)+gs_origin(2)
xyzzyaaan13=dble(xyzzyaaac13)*xyzzyaaar13
gsgrid(3,xyzzyaaad13)=xyzzyaaan13*a3(3)+gs_origin(3)
enddo
enddo
enddo
end subroutine xyzzyaaaa1
subroutine gvsort(r,ip,alow,ahi,istart,nentry,jp,n,nlim,toll)
implicit none
integer,intent(in) :: n,nlim
integer,intent(inout) :: ip(*),istart(*),nentry(*),jp(*)
real(dp),intent(in) :: r(*),toll
real(dp),intent(inout) :: alow(*),ahi(*)
integer xyzzyaaaa14,xyzzyaaab14,xyzzyaaac14,xyzzyaaad14,xyzzyaaae14,xy&
&zzyaaaf14,xyzzyaaag14,xyzzyaaah14,xyzzyaaai14
real(dp) xyzzyaaaj14,xyzzyaaak14,xyzzyaaal14,xyzzyaaam14,xyzzyaaan14
xyzzyaaak14=r(1)
ip(1)=1
xyzzyaaaj14=xyzzyaaak14
do xyzzyaaad14=2,n
ip(xyzzyaaad14)=xyzzyaaad14
xyzzyaaam14=r(xyzzyaaad14)
xyzzyaaak14=min(xyzzyaaak14,xyzzyaaam14)
xyzzyaaaj14=max(xyzzyaaaj14,xyzzyaaam14)
enddo
xyzzyaaai14=1
istart(1)=0
nentry(1)=n
alow(1)=xyzzyaaak14
ahi(1)=xyzzyaaaj14
1 xyzzyaaak14=alow(xyzzyaaai14)
xyzzyaaaj14=ahi(xyzzyaaai14)
if((xyzzyaaaj14-xyzzyaaak14)<toll)goto 2
xyzzyaaag14=nentry(xyzzyaaai14)
xyzzyaaaa14=istart(xyzzyaaai14)
if(xyzzyaaag14>=9)then
xyzzyaaan14=(xyzzyaaak14+xyzzyaaaj14)*0.5d0
ahi(xyzzyaaai14+1)=xyzzyaaaj14
xyzzyaaal14=xyzzyaaaj14
xyzzyaaaj14=xyzzyaaak14
xyzzyaaaf14=0
xyzzyaaah14=0
do xyzzyaaad14=1,xyzzyaaag14
xyzzyaaac14=ip(xyzzyaaaa14+xyzzyaaad14)
xyzzyaaam14=r(xyzzyaaac14)
if(xyzzyaaam14>=xyzzyaaan14)then
xyzzyaaah14=xyzzyaaah14+1
jp(xyzzyaaah14)=xyzzyaaac14
xyzzyaaal14=min(xyzzyaaal14,xyzzyaaam14)
cycle
endif
xyzzyaaaf14=xyzzyaaaf14+1
ip(xyzzyaaaa14+xyzzyaaaf14)=xyzzyaaac14
xyzzyaaaj14=max(xyzzyaaaj14,xyzzyaaam14)
enddo
ahi(xyzzyaaai14)=xyzzyaaaj14
nentry(xyzzyaaai14)=xyzzyaaaf14
xyzzyaaaa14=xyzzyaaaa14+xyzzyaaaf14
if(xyzzyaaaa14>=nlim)goto 1
xyzzyaaai14=xyzzyaaai14+1
istart(xyzzyaaai14)=xyzzyaaaa14
nentry(xyzzyaaai14)=xyzzyaaah14
alow(xyzzyaaai14)=xyzzyaaal14
do xyzzyaaad14=1,xyzzyaaah14
ip(xyzzyaaad14+xyzzyaaaa14)=jp(xyzzyaaad14)
enddo
goto 1
endif
xyzzyaaaf14=xyzzyaaag14-1
do xyzzyaaad14=1,xyzzyaaaf14
xyzzyaaac14=ip(xyzzyaaad14+xyzzyaaaa14)
xyzzyaaam14=r(xyzzyaaac14)
xyzzyaaah14=xyzzyaaad14+1
do xyzzyaaae14=xyzzyaaah14,xyzzyaaag14
xyzzyaaab14=ip(xyzzyaaae14+xyzzyaaaa14)
xyzzyaaan14=r(xyzzyaaab14)
if(xyzzyaaan14>=xyzzyaaam14)cycle
ip(xyzzyaaae14+xyzzyaaaa14)=xyzzyaaac14
xyzzyaaac14=xyzzyaaab14
xyzzyaaam14=xyzzyaaan14
enddo
ip(xyzzyaaad14+xyzzyaaaa14)=xyzzyaaac14
enddo
2 xyzzyaaai14=xyzzyaaai14-1
if(xyzzyaaai14/=0)goto 1
end subroutine gvsort
subroutine print_geometry
use slaarnabt,only : isdiag_int
implicit none
integer xyzzyaaaa15,xyzzyaaab15,xyzzyaaac15
real(dp) xyzzyaaad15(3)
character(80) tmpr
if(isperiodic)then
select case (periodicity)
case(1)
call wout('Simulation cell : 1D periodic (supercell geometry)')
case(2)
call wout('Simulation cell : 2D periodic (supercell geometry)')
case(3)
call wout('Simulation cell : 3D periodic (supercell geometry)')
case default
call errstop('PRINT_GEOMETRY','Periodicity problem.  Bug.')
end select
else
call wout('Simulation cell : Finite (atomic/molecular geometry)')
endif
if(isperiodic)then
select case(periodicity)
case(1)
call wout('Primitive cells : '//trim(i2s(npcells)))
call wout()
call wout('Primitive cell lattice vector (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pa1
call wout(tmpr)
call wout()
call wout('Primitive cell reciprocal lattice vector (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pb1
call wout(tmpr)
call wout()
call wout('Simulation cell lattice vector (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')a1
call wout(tmpr)
call wout()
call wout('Simulation cell reciprocal lattice vector (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')b1
call wout(tmpr)
call wout()
call wout('Radius of line inscribed in Wigner-Seitz cell   ',wigner_se&
&itz_radius,rfmt='(e18.12)')
call wout()
case(2)
if(isdiag_int(3,scell_matrix))then
call wout('Primitive cells : '//trim(i2s(scell_matrix(1,1)))//' x '//t&
&rim(i2s(scell_matrix(2,2)))//' = '//trim(i2s(npcells)))
call wout()
endif
call wout('Primitive cell lattice vectors (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pa1
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pa2
call wout(tmpr)
call wout()
call wout('Primitive cell reciprocal lattice vectors (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pb1
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pb2
call wout(tmpr)
call wout()
call wout('Simulation cell lattice vectors (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')a1
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')a2
call wout(tmpr)
call wout()
call wout('Simulation cell reciprocal lattice vectors (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')b1
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')b2
call wout(tmpr)
call wout()
call wout('Radius of circle inscribed in Wigner-Seitz cell ',wigner_se&
&itz_radius,rfmt='(e18.12)')
call wout()
case(3)
if(isdiag_int(3,scell_matrix))then
call wout('Primitive cells : '//trim(i2s(scell_matrix(1,1)))//' x '//t&
&rim(i2s(scell_matrix(2,2)))//' x '//trim(i2s(scell_matrix(3,3))) //' &
&= '//trim(i2s(npcells)))
call wout()
endif
call wout('Primitive cell lattice vectors (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pa1
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pa2
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pa3
call wout(tmpr)
call wout()
call wout('Primitive cell reciprocal lattice vectors (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pb1
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pb2
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')pb3
call wout(tmpr)
call wout()
call wout('Primitive cell volume                           ',pvolume,r&
&fmt='(e18.12)')
call wout()
call wout('Simulation cell lattice vectors (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')a1
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')a2
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')a3
call wout(tmpr)
call wout()
call wout('Simulation cell reciprocal lattice vectors (au)')
call wout(repeat('-',66))
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')b1
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')b2
call wout(tmpr)
write(tmpr,'(10x,f18.12,1x,f18.12,1x,f18.12)')b3
call wout(tmpr)
call wout()
call wout('Simulation cell volume                          ',volume,rf&
&mt='(e18.12)')
if(periodicity==1)then
call wout('Radius of line inscribed in Wigner-Seitz cell   ',wigner_se&
&itz_radius,rfmt='(e18.12)')
elseif(periodicity==2)then
call wout('Radius of circle inscribed in Wigner-Seitz cell ',wigner_se&
&itz_radius,rfmt='(e18.12)')
else
call wout('Radius of sphere inscribed in Wigner-Seitz cell ',wigner_se&
&itz_radius,rfmt='(e18.12)')
endif
call wout()
end select
endif
if(nbasis>0)then
if(isperiodic)then
call wout()
call wout('Atom Atno  Type  Position (fractional)            Position &
&(Cartesian au)')
call wout(repeat('-',78))
xyzzyaaaa15=0
do xyzzyaaab15=1,npcells
do xyzzyaaac15=1,nbasis
xyzzyaaaa15=xyzzyaaaa15+1
select case (periodicity)
case(1)
xyzzyaaad15(1)=rion(1,xyzzyaaaa15)/pa1(1)
xyzzyaaad15(2:3)=0.d0
case(2)
xyzzyaaad15(1:2)=matmul(rion(1:2,xyzzyaaaa15),painv(1:2,1:2))
xyzzyaaad15(3)=0.d0
case(3)
xyzzyaaad15(1:3)=matmul(rion(1:3,xyzzyaaaa15),painv)
end select
write(tmpr,'(i4,1x,i3,1x,i3,6(1x,f10.4))')xyzzyaaaa15,atno(xyzzyaaac15&
&),iontype(xyzzyaaaa15),xyzzyaaad15(1),xyzzyaaad15(2),xyzzyaaad15(3),r&
&ion(1,xyzzyaaaa15),rion(2,xyzzyaaaa15),rion(3,xyzzyaaaa15)
call wout(tmpr)
enddo
enddo
else
call wout()
call wout('Atom Atno  Type              Position (Cartesian au)')
call wout(repeat('-',61))
do xyzzyaaaa15=1,nbasis
write(tmpr,'(i4,1x,i3,1x,i3,4x,3(1x,f14.8))')xyzzyaaaa15,atno(xyzzyaaa&
&a15),iontype(xyzzyaaaa15),rion(1,xyzzyaaaa15),rion(2,xyzzyaaaa15),rio&
&n(3,xyzzyaaaa15)
call wout(tmpr)
enddo
endif
call xyzzyaaab1
call wout()
if(neighprint==0.and.(trim(atom_basis_type)=='gaussian'.and.cusp_corre&
&ction))neighprint=1
if(neighprint>0)then
call neighbour_analysis(.true.)
call wout()
else
if(nbasis>1)then
call wout('Neighbour analysis turned off (NEIGHPRINT==0).')
call wout()
endif
endif
endif
end subroutine print_geometry
subroutine neighbour_analysis(pr)
implicit none
logical,intent(in) :: pr
integer xyzzyaaaa16,xyzzyaaab16,xyzzyaaac16,xyzzyaaad16,xyzzyaaae16,xy&
&zzyaaaf16,xyzzyaaag16,xyzzyaaah16,xyzzyaaai16,xyzzyaaaj16,xyzzyaaak16&
&,xyzzyaaal16,xyzzyaaam16,xyzzyaaan16,xyzzyaaao16,xyzzyaaap16,xyzzyaaa&
&q16,xyzzyaaar16,xyzzyaaas16,xyzzyaaat16,xyzzyaaau16,xyzzyaaav16
integer,parameter :: xyzzyaaaw16=30,xyzzyaaax16=1000
integer,allocatable,dimension(:) :: xyzzyaaay16,xyzzyaaaz16,xyzzyaaba1&
&6,xyzzyaabb16,xyzzyaabc16,xyzzyaabd16
real(dp) xyzzyaabe16,xyzzyaabf16,xyzzyaabg16,xyzzyaabh16,xyzzyaabi16
real(dp),parameter :: xyzzyaabj16=1.d-8
real(dp),allocatable :: xyzzyaabk16(:)
character(80) tmpr,tmpr2,tmpr3
character(2),allocatable :: xyzzyaabl16(:)
if(.not.pr)call errstop('NEIGHBOUR_ANALYSIS','This routine should only&
& be called by the master node.')
xyzzyaaag16=0
xyzzyaaat16=min(neighprint,xyzzyaaaw16)
xyzzyaaam16=xyzzyaaat16*xyzzyaaax16
xyzzyaaan16=xyzzyaaam16+1
neighprint=neighprint+1
if(isperiodic)then
xyzzyaaau16=num_g
else
xyzzyaaau16=1
endif
allocate(xyzzyaabl16(nbasis),xyzzyaaay16(nbasis),xyzzyaabk16(xyzzyaaan&
&16),xyzzyaaaz16(xyzzyaaan16),xyzzyaaba16(xyzzyaaan16),xyzzyaabb16(xyz&
&zyaaan16),nearest_ion_d(nitype),nearest_ion_r(nbasis),stat=xyzzyaaag1&
&6)
call check_alloc(xyzzyaaag16,'NEIGHBOUR_ANALYSIS','1')
nearest_ion_d(:)=0.d0
nearest_ion_r(:)=1.d100
if(isperiodic)then
call wout('Neighbour analysis of atoms in zero primitive cell')
call wout(repeat('-',50))
if(any(atno(:)==0))call wout('Floating Gaussians not included.')
call wout('  Atom   N     R(Ang)     R(au)                    Neighbou&
&rs')
else
if(count(atno(:)/=0)>1)then
call wout('Neighbour analysis of the atoms in the molecule')
call wout(repeat('-',47))
if(any(atno(:)==0))call wout('Floating Gaussians not included.')
call wout('  Atom   N     R(Ang)     R(au)    Neighbours')
else
call wout('Neighbour analysis of single atom not required.')
return
endif
endif
call timer('NEIGHBOUR_ANALYSIS',.true.)
xyzzyaaah16=nbasis
do xyzzyaaai16=1,xyzzyaaah16
xyzzyaaay16(xyzzyaaai16)=xyzzyaaai16
enddo
xyzzyaaaz16(1)=1
if(model_system)then
do xyzzyaaab16=1,nbasis
if(atno(xyzzyaaab16)>0)then
xyzzyaabl16(xyzzyaaab16)='D+'
elseif(atno(xyzzyaaab16)<0)then
xyzzyaabl16(xyzzyaaab16)='A-'
else
xyzzyaabl16(xyzzyaaab16)=''
endif
enddo
else
xyzzyaabl16(1:nbasis)=periodic_table(atno(1:nbasis))
endif
do xyzzyaaaf16=1,xyzzyaaah16
if(atno(xyzzyaaaf16)==0)cycle
xyzzyaaai16=xyzzyaaay16(xyzzyaaaf16)
xyzzyaaak16=0
do xyzzyaaaj16=1,nbasis
if(atno(xyzzyaaaj16)==0)cycle
xyzzyaabe16=basis(1,xyzzyaaaj16)-basis(1,xyzzyaaai16)
xyzzyaabf16=basis(2,xyzzyaaaj16)-basis(2,xyzzyaaai16)
xyzzyaabg16=basis(3,xyzzyaaaj16)-basis(3,xyzzyaaai16)
do xyzzyaaaa16=1,xyzzyaaau16
if(isperiodic)then
xyzzyaabi16=(xyzzyaabe16+p_lattice(1,xyzzyaaaa16))**2+(xyzzyaabf16+p_l&
&attice(2,xyzzyaaaa16))**2+(xyzzyaabg16+p_lattice(3,xyzzyaaaa16))**2
else
xyzzyaabi16=xyzzyaabe16**2+xyzzyaabf16**2+xyzzyaabg16**2
endif
do xyzzyaaal16=1,xyzzyaaak16
if(xyzzyaabi16>xyzzyaabk16(xyzzyaaal16))cycle
do xyzzyaaab16=xyzzyaaak16,xyzzyaaal16,-1
xyzzyaabk16(xyzzyaaab16+1)=xyzzyaabk16(xyzzyaaab16)
xyzzyaaba16(xyzzyaaab16+1)=xyzzyaaba16(xyzzyaaab16)
xyzzyaabb16(xyzzyaaab16+1)=xyzzyaabb16(xyzzyaaab16)
enddo
goto 2
enddo
xyzzyaaal16=xyzzyaaak16+1
2   if(xyzzyaaak16/=xyzzyaaam16)xyzzyaaak16=xyzzyaaak16+1
xyzzyaabk16(xyzzyaaal16)=xyzzyaabi16
xyzzyaaba16(xyzzyaaal16)=xyzzyaaaj16
xyzzyaabb16(xyzzyaaal16)=xyzzyaaaa16
enddo
enddo
xyzzyaaad16=1
a:do xyzzyaaal16=2,xyzzyaaak16
xyzzyaabi16=xyzzyaabk16(xyzzyaaal16)
if(xyzzyaabi16-xyzzyaabk16(xyzzyaaal16-1)<xyzzyaabj16)then
xyzzyaaao16=atno(xyzzyaaba16(xyzzyaaal16-1))
if(atno(xyzzyaaba16(xyzzyaaal16))==xyzzyaaao16)cycle
xyzzyaaap16=xyzzyaaba16(xyzzyaaal16)
xyzzyaaae16=xyzzyaabb16(xyzzyaaal16)
do xyzzyaaab16=xyzzyaaal16+1,xyzzyaaak16
if(xyzzyaabk16(xyzzyaaab16)-xyzzyaabi16>=xyzzyaabj16)exit
if(atno(xyzzyaaba16(xyzzyaaab16))/=xyzzyaaao16)cycle
xyzzyaabk16(xyzzyaaal16)=xyzzyaabk16(xyzzyaaab16)
xyzzyaaba16(xyzzyaaal16)=xyzzyaaba16(xyzzyaaab16)
xyzzyaabb16(xyzzyaaal16)=xyzzyaabb16(xyzzyaaab16)
xyzzyaabk16(xyzzyaaab16)=xyzzyaabi16
xyzzyaaba16(xyzzyaaab16)=xyzzyaaap16
xyzzyaabb16(xyzzyaaab16)=xyzzyaaae16
cycle a
enddo
endif
xyzzyaaad16=xyzzyaaad16+1
xyzzyaaaz16(xyzzyaaad16)=xyzzyaaal16
enddo a
xyzzyaaaz16(xyzzyaaad16+1)=xyzzyaaak16+1
xyzzyaaas16=min(xyzzyaaad16,neighprint)
xyzzyaaar16=0
do xyzzyaaad16=1,xyzzyaaas16
xyzzyaaab16=xyzzyaaaz16(xyzzyaaad16+1)-xyzzyaaaz16(xyzzyaaad16)
xyzzyaaar16=max(xyzzyaaab16,xyzzyaaar16)
enddo
allocate(xyzzyaabc16(xyzzyaaar16),xyzzyaabd16(xyzzyaaar16),stat=xyzzya&
&aag16)
call check_alloc(xyzzyaaag16,'NEIGHBOUR_ANALYSIS','2')
do xyzzyaaad16=1,xyzzyaaas16
xyzzyaaaq16=0
do xyzzyaaab16=xyzzyaaaz16(xyzzyaaad16),xyzzyaaaz16(xyzzyaaad16+1)-1
xyzzyaaaa16=xyzzyaabb16(xyzzyaaab16)
xyzzyaaaj16=xyzzyaaba16(xyzzyaaab16)
xyzzyaaaq16=xyzzyaaaq16+1
xyzzyaabd16(xyzzyaaaq16)=xyzzyaaaa16
xyzzyaabc16(xyzzyaaaq16)=xyzzyaaaj16
enddo
xyzzyaaab16=xyzzyaaaz16(xyzzyaaad16+1)-1
xyzzyaabh16=sqrt(xyzzyaabk16(xyzzyaaab16))
if(xyzzyaaad16==2)nearest_ion_r(xyzzyaaaf16)=xyzzyaabh16
if(xyzzyaaad16==2)then
xyzzyaaav16=iontype(xyzzyaaaf16)
if(xyzzyaabh16>nearest_ion_d(xyzzyaaav16))nearest_ion_d(xyzzyaaav16)=x&
&yzzyaabh16
endif
if(xyzzyaabh16>xyzzyaabj16.and.pr)then
tmpr2=''
write(tmpr,'(i3,1x,a2,i4,2f11.4,t34)')xyzzyaaaf16,xyzzyaabl16(xyzzyaaa&
&f16),xyzzyaaaq16,xyzzyaabh16*bohr_to_angstrom,xyzzyaabh16
do xyzzyaaac16=1,xyzzyaaaq16
write(tmpr3,'(i4,1x,a2,2x,3i2)')xyzzyaabc16(xyzzyaaac16),xyzzyaabl16(x&
&yzzyaabc16(xyzzyaaac16)),(p_lmng(xyzzyaaab16,xyzzyaabd16(xyzzyaaac16)&
&),xyzzyaaab16=1,3)
tmpr2=trim(tmpr2)//tmpr3
if(mod(xyzzyaaac16,3)==0.or.xyzzyaaac16==xyzzyaaaq16)then
write(tmpr3,'(a,t33,a)')trim(tmpr),trim(tmpr2)
call wout(tmpr3)
tmpr=''
tmpr2=''
endif
enddo
endif
enddo
deallocate(xyzzyaabc16,xyzzyaabd16)
if(xyzzyaaas16>1.and.xyzzyaaaf16/=xyzzyaaah16.and.pr)call wout()
enddo
deallocate(xyzzyaabl16,xyzzyaaay16,xyzzyaaaz16,xyzzyaabk16,xyzzyaaba16&
&,xyzzyaabb16)
call timer('NEIGHBOUR_ANALYSIS',.false.)
end subroutine neighbour_analysis
subroutine map_to_simcell(ldvec,nvec,vec)
implicit none
integer,intent(in) :: ldvec,nvec
real(dp),intent(inout) :: vec(ldvec,nvec)
integer xyzzyaaaa17,xyzzyaaab17,xyzzyaaac17
if(periodicity<1)return
do xyzzyaaaa17=1,nvec
do xyzzyaaab17=1,periodicity
xyzzyaaac17=floor(sum(vec(1:periodicity,xyzzyaaaa17)*ainv(1:periodicit&
&y,xyzzyaaab17)))
vec(1:periodicity,xyzzyaaaa17)=vec(1:periodicity,xyzzyaaaa17)-xyzzyaaa&
&c17*amat(xyzzyaaab17,1:periodicity)
enddo
enddo
end subroutine map_to_simcell
subroutine make_movie(rele)
use slaarnaca, only : zion
implicit none
real(dp),intent(in) :: rele(3,netot)
character(2) symbol
integer xyzzyaaaa18,xyzzyaaab18,xyzzyaaac18,xyzzyaaad18,xyzzyaaae18,xy&
&zzyaaaf18,xyzzyaaag18,xyzzyaaah18
real(dp),allocatable,save :: xyzzyaaai18(:,:),xyzzyaaaj18(:,:)
logical,save :: xyzzyaaak18=.true.
if(xyzzyaaak18)then
allocate(xyzzyaaai18(3,netot),stat=xyzzyaaah18)
call check_alloc(xyzzyaaah18,'MAKE_MOVIE','temp_rele')
if(nitot/=0)then
allocate(xyzzyaaaj18(3,nitot),stat=xyzzyaaah18)
call check_alloc(xyzzyaaah18,'MAKE_MOVIE','temp_rion')
xyzzyaaaj18=rion
call map_to_simcell(3,nitot,xyzzyaaaj18)
endif
xyzzyaaak18=.false.
endif
xyzzyaaai18=rele
call map_to_simcell(3,netot,xyzzyaaai18)
if(movie_supercells)then
xyzzyaaaf18=3**min(periodicity,2)
xyzzyaaag18=3**periodicity
else
xyzzyaaaf18=1
xyzzyaaag18=1
endif
write(io_mov,*)(nitot+netot)*xyzzyaaaf18
write(io_mov,*)'Input geometry'
if(movie_supercells)then
if(nitot/=0)then
do xyzzyaaad18=1,nitot
symbol=periodic_table(mod(atno(mod(xyzzyaaad18-1,nbasis)+1),1000))
do xyzzyaaae18=1,xyzzyaaag18
if(s_lmng(3,xyzzyaaae18)/=0)cycle
write(io_mov,'(a,4f12.6)')symbol,xyzzyaaaj18(1:2,xyzzyaaad18)+s_lattic&
&e(1:2,xyzzyaaae18),xyzzyaaaj18(3,xyzzyaaad18),zion(iontype(xyzzyaaad1&
&8))
enddo
enddo
endif
xyzzyaaac18=0
do xyzzyaaaa18=1,nspin
do xyzzyaaab18=1,nele(xyzzyaaaa18)
xyzzyaaac18=xyzzyaaac18+1
do xyzzyaaae18=1,xyzzyaaag18
if(s_lmng(3,xyzzyaaae18)/=0)cycle
write(io_mov,'(a,4f12.6)')'p'//trim(i2s(xyzzyaaaa18)),xyzzyaaai18(1:2,&
&xyzzyaaac18)+s_lattice(1:2,xyzzyaaae18),xyzzyaaai18(3,xyzzyaaac18),pc&
&harge(xyzzyaaaa18)
enddo
enddo
enddo
else
if(nitot/=0)then
do xyzzyaaad18=1,nitot
write(io_mov,'(a,4f12.6)')periodic_table(mod(atno(mod(xyzzyaaad18-1,nb&
&asis)+1),1000)),xyzzyaaaj18(1:3,xyzzyaaad18),zion(iontype(xyzzyaaad18&
&))
enddo
endif
xyzzyaaac18=0
do xyzzyaaaa18=1,nspin
do xyzzyaaab18=1,nele(xyzzyaaaa18)
xyzzyaaac18=xyzzyaaac18+1
write(io_mov,'(a,4f12.6)')'p'//trim(i2s(xyzzyaaaa18)),xyzzyaaai18(1:3,&
&xyzzyaaac18),pcharge(xyzzyaaaa18)
enddo
enddo
endif
end subroutine make_movie
real(dp) function mol_system_size()
implicit none
integer xyzzyaaaa19,xyzzyaaab19
real(dp) xyzzyaaac19,xyzzyaaad19
real(dp),parameter :: xyzzyaaae19=1.25d0
real(dp),parameter :: xyzzyaaaf19=2.5d0
if(periodicity>0)call errstop('MOL_SYSTEM_SIZE','Called for periodic s&
&ystem.')
mol_system_size=0.d0
if(nitot>1)then
xyzzyaaac19=0.d0
do xyzzyaaaa19=1,nitot-1
do xyzzyaaab19=xyzzyaaaa19+1,nitot
xyzzyaaad19=(rion(1,xyzzyaaaa19)-rion(1,xyzzyaaab19))**2+(rion(2,xyzzy&
&aaaa19)-rion(2,xyzzyaaab19))**2+(rion(3,xyzzyaaaa19)-rion(3,xyzzyaaab&
&19))**2
if(xyzzyaaad19>xyzzyaaac19)xyzzyaaac19=xyzzyaaad19
enddo
enddo
mol_system_size=sqrt(xyzzyaaac19)*xyzzyaaae19
if(mol_system_size<xyzzyaaaf19)mol_system_size=xyzzyaaaf19
elseif(nitot==1)then
mol_system_size=xyzzyaaaf19
endif
end function mol_system_size
subroutine xyzzyaaab1
implicit none
integer xyzzyaaaa20,xyzzyaaab20,xyzzyaaac20,xyzzyaaad20
real(dp) xyzzyaaae20,xyzzyaaaf20,xyzzyaaag20,xyzzyaaah20
logical,allocatable :: xyzzyaaai20(:)
allocate(xyzzyaaai20(nbasis),stat=xyzzyaaaa20)
call check_alloc(xyzzyaaaa20,'CHECK_INVERSION_SYMMETRY','')
xyzzyaaai20(:)=.false.
inversion_symmetry=.false.
a: do xyzzyaaab20=1,nbasis
if(atno(xyzzyaaab20)==0)then
xyzzyaaai20(xyzzyaaab20)=.true.
cycle
endif
do xyzzyaaac20=xyzzyaaab20,nbasis
if(xyzzyaaai20(xyzzyaaac20))cycle
if(atno(xyzzyaaac20)==0)then
xyzzyaaai20(xyzzyaaac20)=.true.
cycle
endif
if(atno(xyzzyaaac20)/=atno(xyzzyaaab20))cycle
xyzzyaaae20=basis(1,xyzzyaaab20)+basis(1,xyzzyaaac20)
xyzzyaaaf20=basis(2,xyzzyaaab20)+basis(2,xyzzyaaac20)
xyzzyaaag20=basis(3,xyzzyaaab20)+basis(3,xyzzyaaac20)
select case(periodicity)
case(0)
xyzzyaaah20=xyzzyaaae20**2+xyzzyaaaf20*2+xyzzyaaag20**2
if(xyzzyaaah20<1.d-10)then
xyzzyaaai20(xyzzyaaab20)=.true.
xyzzyaaai20(xyzzyaaac20)=.true.
cycle a
endif
case(1)
do xyzzyaaad20=1,num_g
xyzzyaaah20=(xyzzyaaae20+p_lattice(1,xyzzyaaad20))**2
if(xyzzyaaah20<1.d-10)then
xyzzyaaai20(xyzzyaaab20)=.true.
xyzzyaaai20(xyzzyaaac20)=.true.
cycle a
endif
enddo
case(2)
do xyzzyaaad20=1,num_g
xyzzyaaah20=(xyzzyaaae20+p_lattice(1,xyzzyaaad20))**2+(xyzzyaaaf20+p_l&
&attice(2,xyzzyaaad20))**2
if(xyzzyaaah20<1.d-10)then
xyzzyaaai20(xyzzyaaab20)=.true.
xyzzyaaai20(xyzzyaaac20)=.true.
cycle a
endif
enddo
case(3)
do xyzzyaaad20=1,num_g
xyzzyaaah20=(xyzzyaaae20+p_lattice(1,xyzzyaaad20))**2+(xyzzyaaaf20+p_l&
&attice(2,xyzzyaaad20))**2+(xyzzyaaag20+p_lattice(3,xyzzyaaad20))**2
if(xyzzyaaah20<1.d-10)then
xyzzyaaai20(xyzzyaaab20)=.true.
xyzzyaaai20(xyzzyaaac20)=.true.
cycle a
endif
enddo
case default
call errstop('CHECK_INVERSION_SYMMETRY','Invalid periodicity.')
end select
enddo
enddo a
call wout()
if(all(xyzzyaaai20(:)))then
inversion_symmetry=.true.
if(isperiodic)then
call wout('Crystal appears to have inversion symmetry.')
else
if(nitot/=1)call wout('Molecule appears to have inversion symmetry.')
endif
else
if(isperiodic)then
call wout('Crystal does not appear to have inversion symmetry.')
else
if(nitot/=1)call wout('Molecule does not appear to have inversion symm&
&etry.')
endif
endif
deallocate(xyzzyaaai20)
end subroutine xyzzyaaab1
subroutine check_kpoints(nk,kvec,time_reversed)
implicit none
integer,intent(in) :: nk
real(dp),intent(in) :: kvec(3,nk)
logical,intent(out),optional :: time_reversed(nk)
integer xyzzyaaaa21,xyzzyaaab21,xyzzyaaac21
real(dp) xyzzyaaad21(3,nk),xyzzyaaae21(3),two_k_o_rlv(3),xyzzyaaaf21(3&
&)
real(dp),parameter :: tol=1.d-8
if(present(time_reversed))time_reversed=.false.
if(isperiodic.and.nk>1)then
xyzzyaaac21=periodicity
xyzzyaaad21(1:xyzzyaaac21,1)=0.d0
do xyzzyaaaa21=2,nk
xyzzyaaad21(1:xyzzyaaac21,xyzzyaaaa21)=matmul(amat(1:xyzzyaaac21,1:xyz&
&zyaaac21),kvec(1:xyzzyaaac21,1)-kvec(1:xyzzyaaac21,xyzzyaaaa21))*one_&
&over_twopi
if(present(time_reversed).and.any(abs(xyzzyaaad21(1:xyzzyaaac21,xyzzya&
&aaa21)-anint(xyzzyaaad21(1:xyzzyaaac21,xyzzyaaaa21)))>tol))then
xyzzyaaad21(1:xyzzyaaac21,xyzzyaaaa21)=matmul(amat(1:xyzzyaaac21,1:xyz&
&zyaaac21),kvec(1:xyzzyaaac21,1)+kvec(1:xyzzyaaac21,xyzzyaaaa21))*one_&
&over_twopi
time_reversed(xyzzyaaaa21)=.true.
endif
if(any(abs(xyzzyaaad21(1:xyzzyaaac21,xyzzyaaaa21)-anint(xyzzyaaad21(1:&
&xyzzyaaac21,xyzzyaaaa21)))>tol))call errstop_master('CHECK_KPOINTS','&
&The k vectors should lie on a grid of simulation-cell reciprocal latt&
&ice vectors, where the grid may be offset from the origin.  But k vec&
&tors 1 and '//trim(i2s(xyzzyaaaa21)) //' are not separated by a simul&
&ation-cell reciprocal-lattice vector.  Please check that the NPCELL b&
&lock in the input file is set to specify the appropriate simulation s&
&upercell.')
do xyzzyaaab21=1,xyzzyaaaa21-1
xyzzyaaaf21(1:xyzzyaaac21)=matmul(scell_mat_inv(1:xyzzyaaac21,1:xyzzya&
&aac21),xyzzyaaad21(1:xyzzyaaac21,xyzzyaaab21)-xyzzyaaad21(1:xyzzyaaac&
&21,xyzzyaaaa21))
if(all(abs(xyzzyaaaf21(1:xyzzyaaac21)-anint(xyzzyaaaf21(1:xyzzyaaac21)&
&))<tol))call errstop_master('CHECK_KPOINTS','Two k points ('//trim(i2&
&s(xyzzyaaab21))//' and ' //trim(i2s(xyzzyaaaa21))//') are equivalent.&
&')
enddo
enddo
endif
if(isperiodic.and.nk>=1)then
xyzzyaaae21(1)=dot_product(binv(1:3,1),binv(1:3,1))
xyzzyaaae21(2)=dot_product(binv(1:3,2),binv(1:3,2))
xyzzyaaae21(3)=dot_product(binv(1:3,3),binv(1:3,3))
call min_image_brute_force(periodicity,kvec(1:3,1),transpose(bmat),bin&
&v,k_offset,xyzzyaaae21)
two_k_o_rlv=2.d0*matmul(k_offset,binv)
if(.not.complex_wf.and.any(abs(two_k_o_rlv -anint(two_k_o_rlv))>tol*sq&
&rt(dot_product(b1,b1))))call errstop_master('CHECK_KPOINTS','If the w&
&ave function is real then the simulation-cell Bloch vector ought to b&
&e half a simulation-cell reciprocal lattice vector.')
endif
end subroutine check_kpoints
subroutine atoms_label_pcell(natoms_listed_pcell,set_atom_list_pcell,n&
&atoms_listed_scell,set_atom_list_scell)
implicit none
integer,intent(in) :: natoms_listed_pcell, set_atom_list_pcell(natoms_&
&listed_pcell)
integer,intent(out) :: natoms_listed_scell,set_atom_list_scell(nitot)
integer xyzzyaaaa22,xyzzyaaab22
set_atom_list_scell(1:natoms_listed_pcell) =set_atom_list_pcell(1:nato&
&ms_listed_pcell)
natoms_listed_scell=natoms_listed_pcell
do xyzzyaaab22=2,npcells
do xyzzyaaaa22=1,natoms_listed_pcell
natoms_listed_scell=natoms_listed_scell+1
set_atom_list_scell(natoms_listed_scell) =set_atom_list_scell(natoms_l&
&isted_scell-natoms_listed_pcell)+nbasis
enddo
enddo
end subroutine atoms_label_pcell
subroutine atoms_label_species(nspecies_listed,label_species, natoms_l&
&isted_scell,set_atom_list_scell)
implicit none
integer,intent(in) :: nspecies_listed,label_species(nspecies_listed)
integer,intent(out) :: natoms_listed_scell,set_atom_list_scell(nitot)
integer xyzzyaaaa23,xyzzyaaab23
natoms_listed_scell=0
do xyzzyaaaa23=1,nitot
do xyzzyaaab23=1,nspecies_listed
if(iontype(xyzzyaaaa23)==label_species(xyzzyaaab23))then
natoms_listed_scell=natoms_listed_scell+1
set_atom_list_scell(natoms_listed_scell)=xyzzyaaaa23
endif
enddo
enddo
end subroutine atoms_label_species
end module slaarnaan

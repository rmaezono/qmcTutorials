module slaarnabg
use dsp
implicit none
integer scell_matrix(3,3),npcells
real(dp) scell_mat_inv(3,3)
character(20) atom_basis_type
logical isperiodic,model_system,homogeneous_system
integer dimensionality
integer nbasis,periodicity
integer,dimension(:),allocatable :: atno
real(dp) pa1(3),pa2(3),pa3(3)
real(dp),dimension(:,:),allocatable :: basis
real(dp) pb1(3),pb2(3),pb3(3),pvolume
integer nitot,nitype
integer,allocatable,dimension(:) :: iontype,nion,ion_prim,iontype_prim
integer,allocatable,dimension(:,:,:,:) :: nbasis_to_nitot
real(dp) a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),volume,b11,b22,b33,b12,b1&
&3,b23,pb11,pb22,pb33,pb12,pb13,pb23,area,wigner_seitz_radius
real(dp),allocatable,dimension(:,:) :: rion,gpcell
real(dp),allocatable,dimension(:,:,:) :: rionion
real(dp) pamat(3,3),pbmat(3,3),amat(3,3),bmat(3,3)
real(dp) ainv(3,3),painv(3,3),binv(3,3),pbinv(3,3)
integer,parameter :: lsize=500,num_g=lsize+lsize-1
integer p_lmng(3,num_g),pr_lmng(3,num_g),s_lmng(3,num_g),sr_lmng(3,num&
&_g),first_p_in_star(lsize+1),first_pr_in_star(lsize+1),first_s_in_sta&
&r(lsize+1),first_sr_in_star(lsize+1),ngpri
integer,allocatable :: spin_lattice(:),lat_point(:,:,:)
real(dp) p_lattice(3,num_g),pr_lattice(3,num_g),s_lattice(3,num_g),sr_&
&lattice(3,num_g),p_modsq(lsize+1),pr_modsq(lsize+1),s_modsq(lsize+1),&
&sr_modsq(lsize+1)
logical orthogonal,axis_aligned,fcc,bcc
real(dp) scubelen,inv_scubelen,ma1(3),ma2(3),ma3(3),tb(3,3)
integer neighprint
integer,parameter :: ngsgrid=10,ngsgrid_minus_1=9,ngsgrid3=1000
real(dp),allocatable :: nearest_ion_d(:),nearest_ion_r(:)
logical inversion_symmetry
integer nitot_forces,naxis_forces
integer which_ion_displaced
real(dp) ion_displacement(3)
integer setting
character(12) crystal_system
logical :: hard_sphere=.false.
logical :: hard_op_spins=.false.
real(dp) :: hard_diam=-1.d0,hard_diam_sq=-1.d0
logical :: ignore_ionic_interactions=.false.
end module slaarnabg

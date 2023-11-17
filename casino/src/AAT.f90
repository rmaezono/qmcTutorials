module slaarnaat
use dsp
implicit none
integer num_shells,num_ao,num_k,num_real_k,num_complex_k,num_real_k_pl&
&us_1,num_centres,num_cell,maxb,nemaxc,ndim,nnn(3),spin_dep_gs,spin_de&
&p_in,spin_dep_full,size_rck,size_cck
integer,dimension(:),allocatable :: shell_am,numao_in_shell,first_shel&
&l,first_ao,primitive,num_centres_in_cell,shell_sequence_number,num_un&
&ique_orbitals,gauss_offset,gauss_offsetc
integer,dimension(:,:),allocatable :: num_shells_on_centre,virtual_k,n&
&band,which_ion,shell_sequence_mol
integer,dimension(:,:,:),allocatable :: cell_map
real(dp) screening_tolerance,rnorm,new_origin(3),new_origin_to_rvec(3)&
&,one_over_rnorm,two_rnorm,gs_origin(3),gautol,cusp_threshold,cusp_con&
&trol
real(dp),dimension(:),allocatable :: min_exponent,exponent,c_prim,c_pr&
&im2,phase,coskdotg_shift,ao_m,alap_m,agra1_m,agra2_m,agra3_m,asderiv1&
&_m,asderiv2_m,asderiv3_m,asderiv4_m,asderiv5_m,asderiv6_m,psi,fgra1,f&
&gra2,fgra3,flap,orb1,orb2,orb3,orb4,orb5,orb6,orb7,orb8,orb9,orb10,or&
&b11,fsderiv1,fsderiv2,fsderiv3,fsderiv4,fsderiv5,fsderiv6,valence_cha&
&rge
real(dp),dimension(:,:),allocatable :: kvec,xpos_in_cell,ypos_in_cell,&
&zpos_in_cell,coskdotg,rbf,rblap,rbgra1,rbgra2,rbgra3,rbf_s,rblap_s,rb&
&gra1_s,rbgra2_s,rbgra3_s,rcusp_sq_max,rbsderiv1,rbsderiv2,rbsderiv3,r&
&bsderiv4,rbsderiv5,rbsderiv6,rbsderiv1_s,rbsderiv2_s,rbsderiv3_s,rbsd&
&eriv4_s,rbsderiv5_s,rbsderiv6_s,gsgrid
real(dp),dimension(:,:,:),allocatable :: rck_ex
real(dp),dimension(:,:,:,:),pointer :: rck=>null()
complex(dp),dimension(:),allocatable :: expikdotg_shift,cphase,corb1,c&
&orb2,corb3,corb4,corb5,corb6,corb7,corb8,corb9,corb10,corb11
complex(dp),dimension(:,:),allocatable :: expikdotg,cbf,cblap,cbgra1,c&
&bgra2,cbgra3,cbf_s,cblap_s,cbgra1_s,cbgra2_s,cbgra3_s,cbsderiv1,cbsde&
&riv2,cbsderiv3,cbsderiv4,cbsderiv5,cbsderiv6,cbsderiv1_s,cbsderiv2_s,&
&cbsderiv3_s,cbsderiv4_s,cbsderiv5_s,cbsderiv6_s
complex(dp),dimension(:,:,:),allocatable :: cck_ex
complex(dp),dimension(:,:,:,:),pointer :: cck=>null()
logical excite,complex_states,cusp_correction,molgscreening,num_real_k&
&_gt_1,num_real_k_odd,s_plot,printgscreening,cusp_info
character(80) code
real(dp),dimension(:),allocatable :: exp_poly0,exp_poly1,exp_poly2,exp&
&_poly3,exp_poly4,exp_poly5,exp_poly6,exp_poly7,exp_poly8,exp_poly9,ex&
&p_poly10
real(dp),dimension(:,:,:),allocatable :: rcusp_sq,disign,pshift
real(dp),dimension(:,:,:,:),allocatable :: acusp
real(dp),allocatable :: cusp_val_r(:),cusp_grad_r(:,:),cusp_lap_r(:),c&
&usp_sderivs_r(:,:)
complex(dp),allocatable :: cusp_val_c(:),cusp_grad_c(:,:),cusp_lap_c(:&
&),cusp_sderivs_c(:,:)
integer gauss_gs_norb,gauss_ex_norb
integer,allocatable :: ridx2zidx(:),full2vrt(:,:),vrt2full(:,:)
end module slaarnaat

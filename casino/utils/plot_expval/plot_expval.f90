 PROGRAM plot_expval
!-----------------------------------------------------------------------------!
! PLOT_EXPVAL                                                                 !
! Expectation value plotting program                                          !
!                                                                             !
! Reads expval.data file produced by CASINO, performs appropriate conversions !
! on the data therein, and writes the results to files visualizable by        !
! xmgr, gnuplot, etc..                                                        !
!                                                                             !
! Currently supports :                                                        !
!                                                                             !
! (1)  density                                                                !
! (2)  spin density                                                           !
! (3)  spin density matrix                                                    !
! (4)  reciprocal space pair-correlation function                             !
! (5)  spherical pair correlation function                                    !
! (6)  localization tensor (NOT YET)                                          !
! (7)  structure factor                                                       !
! (8)  spherically-averaged structure factor                                  !
! (9)  one-particle density matrix                                            !
! (10) two-particle density matrix                                            !
! (11) condensate fraction estimator (unbiased TBDM, goes as TBDM-OBDM**2)    !
! (12) momentum density                                                       !
! (13) finite density                                                         !
! (14) molecular density                                                      !
! (15) two-particle density matrix (momentum space)                           !
! (16) condensate fraction estimator (momentum space)                         !
!                                                                             !
! MDT 11.2005                                                                 !
!                                                                             !
! CHANGES                                                                     !
! =======                                                                     !
! NDD 04.2007  Added ability to plot expvals for 1D systems.                  !
! NDD 11.2007  Fixed bug affecting charge density in systems w/o inv. sym.    !
! PS  10.2009  Added functionality for finite density.                        !
! NDD 05.2012  Added ability to plot the spin density matrix.                 !
! NDD 11.2013  Enabled nondiagonal supercell matrices.                        !
!-----------------------------------------------------------------------------!
 USE dsp
 USE esdf
 USE run_control
 USE format_utils, ONLY : write_list_int,i2s,wordwrap
 IMPLICIT NONE

 INTEGER,PARAMETER :: num_expvals=16
 INTEGER ialloc,ierr,header_lines,line,npoint1,npoint2,npoint3,nex,netot,&
  &inex(num_expvals),gset(num_expvals),iplot,iset,n1,n2,n3
 INTEGER,PARAMETER :: max_header_lines=10,io_expval=10
 REAL(dp) timag,treal,vec1(3),vec2(3),vec3(3)
 REAL(dp),PARAMETER :: pi=3.14159265358979324d0,twopi=2.d0*pi,&
  &one_over_twopi=1.d0/twopi
 LOGICAL expval_in_expval,density_in_expval,&
  &spin_density_in_expval,spin_density_mat_in_expval,pair_corr_in_expval,&
  &pair_corr_sph_in_expval,loc_tensor_in_expval,structure_factor_in_expval,&
  &structure_factor_sph_in_expval,onep_density_mat_in_expval,&
  &twop_density_mat_in_expval,cond_fraction_in_expval,mom_den_in_expval,&
  &file_present,twop_dm_with_onep_dm,twop_dm_rescale,finite_density_in_expval,&
  &mol_density_in_expval,twop_dm_mom_in_expval,twop_dm_mom_rescale,&
  &cond_fraction_mom_in_expval
 CHARACTER(80) :: char_80
 CHARACTER(80),ALLOCATABLE :: header_text(:)

! Basic info
 INTEGER expval_periodicity,expval_scell_matrix(3,3),npcells,expval_ngsets,&
  &file_version,num_particle_types,expval_dimensionality
 INTEGER,ALLOCATABLE :: nele_expval(:)
 REAL(dp) simcell_area,simcell_length,simcell_volume,pa1(3),pa2(3),pa3(3),&
  &pamat(3,3),pbmat(3,3),amat(3,3),bmat(3,3),wigner_seitz_radius
 CHARACTER(80) title

! Plot description
 INTEGER nlines,plot_dim
 REAL(dp) apoint(3),bpoint(3),cpoint(3),dpoint(3),xprod(3),r_ab(3),r_ac(3),&
  &r_ad(3),step_ab(3),step_ac(3),step_ad(3),r(3),r1(3),step_ab_length,expval,&
  &expval2

! Parameters for G vector set generation
! INTEGER,PARAMETER :: ngrid_0d=128
! INTEGER,PARAMETER :: ngrid_1d(3)=(/64,128,128/)
! INTEGER,PARAMETER :: ngrid_2d(3)=(/64,64,128/)
! INTEGER,PARAMETER :: ngrid_3d=64

! G vector sets
 INTEGER expval_grange
 INTEGER,ALLOCATABLE :: expval_ngvec(:),expval_pgmap(:,:)
 REAL(dp),ALLOCATABLE :: e_cutoff(:),expval_gvec(:,:,:),expval_amat(:,:,:),&
  &expval_bmat(:,:,:)
 COMPLEX(dp),ALLOCATABLE :: expval_expigdotr(:),expval_mwork(:,:)

! QMC density
 INTEGER den_gset,den_nsets
 INTEGER,ALLOCATABLE :: den_ptype(:)
 REAL(dp),ALLOCATABLE :: expval_den_weight(:)
 COMPLEX(dp),ALLOCATABLE :: expval_den(:,:)
 CHARACTER(3) accum_method_den

! Spin density
 INTEGER sden_gset,sden_nsets
 INTEGER,ALLOCATABLE :: sden_ptype(:)
 REAL(dp),ALLOCATABLE :: expval_sden_weight(:)
 COMPLEX(dp),ALLOCATABLE :: expval_sden(:,:)
 CHARACTER(3) accum_method_sden

! Spin density matrix
 INTEGER sdenm_gset,sdenm_nsets
 INTEGER,ALLOCATABLE :: sdenm_ptype(:)
 REAL(dp),ALLOCATABLE :: expval_sdenm_weight(:)
 COMPLEX(dp),ALLOCATABLE :: expval_sdenm(:,:,:)
 CHARACTER(3) accum_method_sdenm

! Reciprocal space pair correlation function
 INTEGER pcf_gset,pcf_nsets,pcf_type_fixed
 INTEGER,ALLOCATABLE :: pcf_ptype(:)
 REAL(dp) :: pcf_rfix(3)
 REAL(dp),ALLOCATABLE :: expval_pcf_weight(:)
 COMPLEX(dp),ALLOCATABLE :: expval_pcf(:,:)
 CHARACTER(3) accum_method_pcf

! Spherical pair correlation function
 INTEGER pcfs_nsets,pcfs_nbins,pcfs_nbins_new,pcfs_accmode,pcfs_type_fixed
 INTEGER,ALLOCATABLE :: pcfs_ptype(:,:)
 REAL(dp) pcfs_rcutoff,pcfs_rfix(3)
 REAL(dp),ALLOCATABLE :: expval_pcfs_weight(:),expval_pcfs(:,:,:),&
  &expval_pcfs_reblocked(:,:,:)
 LOGICAL reblocked
 CHARACTER(3) accum_method_pcfs

! Localization tensor
 INTEGER lt_nsets

! Structure factor
 INTEGER sf_nsets,sf_gset,sfsden_nsets,sf_option
 INTEGER,ALLOCATABLE :: sf_ptype(:,:),sfsden_ptype(:)
 REAL(dp) expval_sf_weight,expval_sfsden_weight
 REAL(dp),ALLOCATABLE :: expval_sf(:,:)
 COMPLEX(dp),ALLOCATABLE :: expval_sfsden(:,:)
 LOGICAL homogeneous_structure_factor
 CHARACTER(3) accum_method_sf

! Spherical structure factor
 INTEGER sf_sph_nsets,sf_sph_k1,sf_sph_k2,sf_sph_nk
 INTEGER,ALLOCATABLE :: sf_sph_ptype(:,:)
 REAL(dp),ALLOCATABLE :: expval_sf_sph_k(:),expval_sf_sph_weight(:),&
  &expval_sf_sph(:,:)
 CHARACTER(3) accum_method_sf_sph

! On-top OBDM / momentum density particle types.
 INTEGER :: on_top_ispin=0,on_top_jspin=0

! One-particle density matrix
 INTEGER onep_dm_nsets,onep_dm_nbins
 INTEGER,ALLOCATABLE :: onep_dm_nptypes_in_set(:),onep_dm_ptype_in_set(:,:)
 REAL(dp),ALLOCATABLE :: expval_onep_dm_weight(:,:),&
  &expval_onep_dm(:,:),expval_onep_dm2(:,:)
 CHARACTER(3) accum_method_onep_dm

! Two-particle density matrix
 INTEGER twop_dm_nsets,twop_dm_nbins
 INTEGER,ALLOCATABLE :: twop_dm_nptypes_in_set(:),twop_dm_ptype_in_set(:,:,:)
 REAL(dp),ALLOCATABLE :: expval_twop_dm_weight(:,:),expval_twop_dm(:,:),&
  &expval_twop_dm2(:,:)
 CHARACTER(3) accum_method_twop_dm

! Condensate fraction
 INTEGER cond_frac_nsets,cond_frac_nbins
 INTEGER,ALLOCATABLE :: cond_frac_nptypes_in_set(:),&
  &cond_frac_ptype_in_set(:,:,:)
 REAL(dp),ALLOCATABLE :: expval_cond_frac_weight(:,:),expval_cond_frac(:,:),&
  &expval_cond_frac2(:,:)
 CHARACTER(3) accum_method_cond_frac

! Momentum density
 INTEGER mom_den_nsets
 INTEGER,ALLOCATABLE :: mom_den_nptypes_in_set(:),mom_den_ptype_in_set(:,:)
 REAL(dp),ALLOCATABLE :: expval_mom_den_weight(:),expval_mom_den(:,:),&
  &expval_mom_den2(:,:),expval_mom_den_k(:,:)
 CHARACTER(3) accum_method_mom_den

! Finite Density
 INTEGER,PARAMETER :: fin_den_f_indx_step=1,fin_den_f_indx_bessel=2,&
  &fin_den_f_indx_cheb=3,fin_den_f_indx_exp=4,&
  &fin_den_power(6)=(/-2,-1,0,1,2,3/)
 INTEGER :: fin_den_nsets_ij,fin_den_nsets_i,fin_den_norder,fin_den_ibasis,&
  &fin_den_z,io_rad_mom=11
 REAL(dp),ALLOCATABLE :: expval_fin_den_ij(:,:),expval_fin_den2_ij(:,:),&
  &expval_fin_den_i(:,:),expval_fin_den2_i(:,:)
 REAL(dp) fin_den_cutoff,fin_den_tail_offset,expval_fin_den_nstep,&
  &expval_fin_den_weight,expval_fin_den_weight2,fin_den_ip
 CHARACTER(3) accum_method_fin_den
 CHARACTER(15) fin_den_basis

! Molecular Density
 INTEGER mol_den_nsets,mol_den_norder(3)
 REAL(dp),ALLOCATABLE :: expval_mol_den(:,:,:,:)
 REAL(dp) expval_mol_den_nstep,expval_mol_den_weight,mol_den_mincoord(3),&
  &mol_den_maxcoord(3)
 CHARACTER(3) accum_method_mol_den

! Two-particle density matrix (momentum space)
 INTEGER twop_dm_mom_nsets
 INTEGER,ALLOCATABLE :: twop_dm_mom_nptypes_in_set(:),&
  &twop_dm_mom_ptype_in_set(:,:,:)
 REAL(dp),ALLOCATABLE :: expval_twop_dm_mom_weight(:),expval_twop_dm_mom(:,:),&
  &expval_twop_dm_mom2(:,:),expval_twop_dm_mom_k(:,:)
 CHARACTER(3) accum_method_twop_dm_mom

! Condensate fraction (momentum space)
 INTEGER cond_frac_mom_nsets
 INTEGER,ALLOCATABLE :: cond_frac_mom_nptypes_in_set(:),&
  &cond_frac_mom_ptype_in_set(:,:,:)
 REAL(dp),ALLOCATABLE :: expval_cond_frac_mom_weight(:),&
  &expval_cond_frac_mom(:,:),expval_cond_frac_mom2(:,:),&
  &expval_cond_frac_mom_k(:,:)
 CHARACTER(3) accum_method_cond_frac_mom

 gset(:)=0

 write(6,*)
 write(6,*)'PLOT_EXPVAL : visualization of CASINO expectation values'
 write(6,*)'========================================================'
 write(6,*)

 inquire(file='expval.data',exist=file_present)
 if(.not.file_present)then
  write(6,*)"No file 'expval.data' containing the expectation value data."
  write(6,*)
  call errstop('PLOT_EXPVAL','Quitting.')
 endif

 call read_expval

 if(density_in_expval.or.spin_density_in_expval.or.spin_density_mat_in_expval&
  &.or.pair_corr_in_expval)then
  inquire(file='input',exist=file_present)
  if(.not.file_present)then
   write(6,*)"No file 'input' containing the controlling plot_expval block."
   write(6,*)
   call errstop('PLOT_EXPVAL','Quitting.')
  endif
  call read_input
 else
  write(6,*)'Skipping unnecessary input read.'
  write(6,*)
 endif

 call talk_to_user

 call write_expval

 stop


CONTAINS


 SUBROUTINE read_expval
!-----------------------------------------------------------------------!
! Call routine to determine what data blocks expval.data contains, then !
! call the relevant additional routines to read these blocks.           !
!-----------------------------------------------------------------------!
 IMPLICIT NONE

 write(6,*)'Found expval.data file. Reading data.'
 write(6,*)

 open(unit=io_expval,file='expval.data',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_EXPVAL','Problem opening expval.data.')

! Find out which data sets are present in expval.data.
 call find_expval_sets

! Detailed read of 'START EXPVAL' block
 call read_expval_block

! If there are G vector sets given, then read them.
 if(expval_ngsets>0)then
  call read_g_vector_sets
 endif

! Detailed read of whatever sets are present
 if(density_in_expval)call read_density
 if(spin_density_in_expval)call read_spin_density
 if(spin_density_mat_in_expval)call read_spin_density_matrix
 if(pair_corr_in_expval)call read_pair_corr
 if(pair_corr_sph_in_expval)call read_pair_corr_sph
 if(loc_tensor_in_expval)call read_loc_tensor
 if(structure_factor_in_expval)call read_structure_factor
 if(structure_factor_sph_in_expval)call read_structure_factor_sph
 if(onep_density_mat_in_expval)call read_onep_density_matrix
 if(twop_density_mat_in_expval)call read_twop_density_matrix
 if(cond_fraction_in_expval)call read_cond_fraction
 if(mom_den_in_expval)call read_mom_den
 if(finite_density_in_expval)call read_finite_density
 if(mol_density_in_expval)call read_mol_density
 if(twop_dm_mom_in_expval)call read_twop_dm_mom
 if(cond_fraction_mom_in_expval)call read_cond_fraction_mom

 close(io_expval)

 END SUBROUTINE read_expval


 SUBROUTINE read_input
 IMPLICIT NONE
 INTEGER j
 CHARACTER(16) plot_dim_string
 REAL(dp),PARAMETER :: tol=1.d-13

 call esdf_init('input') ! i.e. read the input file.
 call esdf_warnout

! Read the plot_expval block.
 npoint1=0 ; npoint2=0 ; npoint3=0
 apoint=0.d0 ; bpoint=0.d0 ; cpoint=0.d0 ; dpoint=0.d0

 if(esdf_block('plot_expval',nlines))then

  write(6,*)'Found input file. Reading plot_expval block.'
  write(6,*)

! Item to plot
! Dimensionality
  read(block_data(1),*,err=2,end=3)plot_dim_string
  select case(trim(adjustl(plot_dim_string)))
  case('line','LINE','1','1D','1d','1-D','1-d')
   plot_dim=1
  case('plane','PLANE','2','2D','2d','2-D','2-d')
   plot_dim=2
  case('volume','VOLUME','3','3D','3d','3-D','3-d')
   plot_dim=3
  case default
   write(6,*)'The second line in the plot_expval block must be one of:'
   write(6,*)' line, plane, volume (or similar).'
   call errstop('READ_INPUT','Stopping.')
  end select
! No. of points to plot, starting point and ending point A & B [& C [& D]]
  select case(plot_dim)
  case(1)
   read(block_data(2),*,err=2,end=3)npoint1
   read(block_data(3),*,err=2,end=3)(apoint(j),j=1,3)
   read(block_data(4),*,err=2,end=3)(bpoint(j),j=1,3)
  case(2)
   read(block_data(2),*,err=2,end=3)npoint1,npoint2
   read(block_data(3),*,err=2,end=3)(apoint(j),j=1,3)
   read(block_data(4),*,err=2,end=3)(bpoint(j),j=1,3)
   read(block_data(5),*,err=2,end=3)(cpoint(j),j=1,3)
  case(3)
   read(block_data(2),*,err=2,end=3)npoint1,npoint2,npoint3
   read(block_data(3),*,err=2,end=3)(apoint(j),j=1,3)
   read(block_data(4),*,err=2,end=3)(bpoint(j),j=1,3)
   read(block_data(5),*,err=2,end=3)(cpoint(j),j=1,3)
   read(block_data(6),*,err=2,end=3)(dpoint(j),j=1,3)
  end select

! Check line lengths and collinearity.
  r_ab=bpoint-apoint ; treal=sqrt(sum(r_ab(:)**2))
  if(treal==0.d0)call errstop('READ_INPUT','Length AB zero in plot_expval.')
  if(plot_dim>1)then
   r_ac=cpoint-apoint ; treal=sqrt(sum(r_ac(:)**2))
   if(treal==0.d0)call errstop('READ_INPUT','Length AC zero in plot_expval.')
   xprod(1)=r_ab(2)*r_ac(3)-r_ab(3)*r_ac(2)
   xprod(2)=r_ab(3)*r_ac(1)-r_ab(1)*r_ac(3)
   xprod(3)=r_ab(1)*r_ac(2)-r_ab(2)*r_ac(1)
   if(all(abs(xprod)<tol))call errstop('READ_INPUT',&
    &'Collinear points A, B and C in plot_expval.')
  endif
  if(plot_dim>2)then
   r_ad=dpoint-apoint ; treal=sqrt(sum(r_ad(:)**2))
   if(treal==0.d0)call errstop('READ_INPUT','Length AD zero in plot_expval.')
   xprod(1)=r_ab(2)*r_ad(3)-r_ab(3)*r_ad(2)
   xprod(2)=r_ab(3)*r_ad(1)-r_ab(1)*r_ad(3)
   xprod(3)=r_ab(1)*r_ad(2)-r_ab(2)*r_ad(1)
   if(all(abs(xprod)<tol))call errstop('READ_INPUT',&
    &'Collinear points A, B and D in qmc_plot.')
   xprod(1)=r_ac(2)*r_ad(3)-r_ac(3)*r_ad(2)
   xprod(2)=r_ac(3)*r_ad(1)-r_ac(1)*r_ad(3)
   xprod(3)=r_ac(1)*r_ad(2)-r_ac(2)*r_ad(1)
   if(all(abs(xprod)<tol))call errstop('READ_INPUT',&
    &'Collinear points A, C and D in plot_expval.')
  endif

! Define steps along each line to get the required number of points.
  step_ab=r_ab/dble(npoint1-1)
  if(plot_dim>1)step_ac=r_ac/dble(npoint2-1)
  if(plot_dim>2)step_ad=r_ad/dble(npoint3-1)

  select case(plot_dim)
  case(1) ; step_ab_length=sqrt(step_ab(1)**2+step_ab(2)**2+step_ab(3)**2)
  case default ; step_ab_length=0.d0
  end select

 else

  call errstop('READ_INPUT','The input file does not contain a plot_expval &
   &block')

 endif

 return

2 call errstop('READ_INPUT','Error reading plot_expval block.')
3 call errstop('READ_INPUT','Reached end of line reading plot_expval block.')

 END SUBROUTINE read_input


 SUBROUTINE find_expval_sets
!------------------------------------------------------------------!
! If expval.data exists at the start of the run, do a cursory read !
! and note which expectation values are already present.           !
!------------------------------------------------------------------!
 IMPLICIT NONE
 header_lines=0
 density_in_expval=.false.          ; spin_density_in_expval=.false.
 spin_density_mat_in_expval=.false. ; pair_corr_in_expval=.false.
 pair_corr_sph_in_expval=.false.    ; loc_tensor_in_expval=.false.
 structure_factor_in_expval=.false. ; structure_factor_sph_in_expval=.false.
 onep_density_mat_in_expval=.false. ; twop_density_mat_in_expval=.false.
 cond_fraction_in_expval=.false.    ; mom_den_in_expval=.false.
 finite_density_in_expval=.false.   ; mol_density_in_expval=.false.
 twop_dm_mom_in_expval=.false.      ; cond_fraction_mom_in_expval=.false.
 expval_in_expval=.false.

 allocate(header_text(max_header_lines),stat=ierr)
 if(ierr/=0)call errstop('FIND_EXPVAL_SETS','Allocation problem: &
  &header_text.')

! Look for headers
 do

  read(io_expval,'(a)',iostat=ierr)char_80

  if(ierr>0)then
   call errstop('FIND_EXPVAL_SETS','Problem reading expval.data.')
  elseif(ierr<0)then
   if(header_lines==0)write(6,*)'No header is present in this file.'
   exit
  endif

  if(trim(adjustl(char_80))=='START HEADER')then

   if(header_lines>0)call errstop('FIND_EXPVAL_SETS',&
    &'There is more than one header in expval.data.')
   write(6,*)'HEADER:'
   do line=1,max_header_lines
    read(io_expval,'(a)',iostat=ierr)header_text(line)
    if(ierr/=0)call errstop('FIND_EXPVAL_SETS','Problem reading expval.data.')
    header_text(line)=adjustl(header_text(line))
    if(trim(adjustl(header_text(line)))=='END HEADER')exit
    write(6,*)" "//trim(header_text(line))
    header_lines=header_lines+1
    if(header_lines>=max_header_lines)call errstop&
     &('FIND_EXPVAL_SETS','Header is more than '&
     &//trim(i2s(max_header_lines-1))//' long. &
     &Please try to be more concise.')
   enddo ! line

  elseif(trim(adjustl(char_80))=='START EXPVAL')then

   write(6,*)
   write(6,*)'EXPECTATION VALUES PRESENT:'
   nex=0
   do ! headers in the EXPVAL block
    read(io_expval,'(a)',iostat=ierr)char_80
    if(trim(adjustl(char_80))=='END EXPVAL')exit
    if(ierr>0)then
     call errstop('FIND_EXPVAL_SETS','Problem reading expval.data <2>.')
    elseif(ierr<0)then
     exit
    endif
    if(trim(adjustl(char_80))=='START DENSITY')then
     if(density_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one DENSITY block in expval.data.')
     nex=nex+1 ; inex(nex)=1 ; density_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') density'
    elseif(trim(adjustl(char_80))=='START SPIN DENSITY')then
     if(spin_density_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one SPIN DENSITY block in expval.data.')
     nex=nex+1 ; inex(nex)=2 ; spin_density_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') spin density'
    elseif(trim(adjustl(char_80))=='START SPIN-DENSITY MATRIX')then
     if(spin_density_mat_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one SPIN-DENSITY MATRIX block in expval.data.')
     nex=nex+1 ; inex(nex)=3 ; spin_density_mat_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') spin density matrix'
    elseif(trim(adjustl(char_80))=='START RECIPROCAL-SPACE PCF')then
     if(pair_corr_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one RECIPROCAL-SPACE PCF block in expval.data.')
     nex=nex+1 ; inex(nex)=4 ; pair_corr_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') reciprocal space pair correlation function'
    elseif(trim(adjustl(char_80))=='START SPHERICAL PCF')then
     if(pair_corr_sph_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one SPHERICAL PCF block in expval.data.')
     nex=nex+1 ; inex(nex)=5 ; pair_corr_sph_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') spherical pair correlation function'
    elseif(trim(adjustl(char_80))=='START LOCALIZATION TENSOR')then
     if(loc_tensor_in_expval)call errstop('FIND_EXPVAL_SETS', &
      &'More than one LOCALIZATION TENSOR block in expval.data.')
     nex=nex+1 ; inex(nex)=6 ; loc_tensor_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') localization tensor'
    elseif(trim(adjustl(char_80))=='START STRUCTURE FACTOR')then
     if(structure_factor_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one STRUCTURE FACTOR block in expval.data.')
     nex=nex+1 ; inex(nex)=7 ; structure_factor_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') structure factor'
    elseif(trim(adjustl(char_80))=='START SPHERICAL STRUCTURE FACTOR')then
     if(structure_factor_sph_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one SPHERICAL STRUCTURE FACTOR block in expval.data.')
     nex=nex+1 ; inex(nex)=8 ; structure_factor_sph_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') spherical structure factor'
    elseif(trim(adjustl(char_80))=='START ONE-PARTICLE DENSITY MATRIX')then
     if(onep_density_mat_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one ONE-PARTICLE DENSITY MATRIX block in expval.data.')
     nex=nex+1 ; inex(nex)=9 ; onep_density_mat_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') one-particle density matrix'
    elseif(trim(adjustl(char_80))=='START TWO-PARTICLE DENSITY MATRIX')then
     if(twop_density_mat_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one TWO-PARTICLE DENSITY MATRIX block in expval.data.')
     nex=nex+1 ; inex(nex)=10 ; twop_density_mat_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') two-particle density matrix'
    elseif(trim(adjustl(char_80))=='START CONDENSATE FRACTION')then
     if(cond_fraction_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one CONDENSATE FRACTION block in expval.data.')
     nex=nex+1 ; inex(nex)=11 ; cond_fraction_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') condensate fraction'
    elseif(trim(adjustl(char_80))=='START MOMENTUM DENSITY')then
     if(mom_den_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one MOMENTUM DENSITY block in expval.data.')
     nex=nex+1 ; inex(nex)=12 ; mom_den_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') momentum density'
    elseif(trim(adjustl(char_80))=='START FINITE DENSITY')then
     if(finite_density_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one FINITE DENSITY block in expval.data.')
     nex=nex+1 ; inex(nex)=13 ; finite_density_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') finite density'
    elseif(trim(adjustl(char_80))=='START MOLECULAR DENSITY')then
     if(mol_density_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one MOLECULAR DENSITY block in expval.data.')
     nex=nex+1 ; inex(nex)=14 ; mol_density_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') molecular density'
    elseif(trim(adjustl(char_80))=='START TWO-PARTICLE DENSITY MATRIX &
     &(MOMENTUM SPACE)')then
     if(twop_dm_mom_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one TWO-PARTICLE DENSITY MATRIX (MOMENTUM SPACE) block in &
      &expval.data.')
     nex=nex+1 ; inex(nex)=15 ; twop_dm_mom_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') two-particle density matrix &
      &(momentum space)'
    elseif(trim(adjustl(char_80))=='START CONDENSATE FRACTION &
     &(MOMENTUM SPACE)')then
     if(cond_fraction_mom_in_expval)call errstop('FIND_EXPVAL_SETS',&
      &'More than one CONDENSATE FRACTION (MOMENTUM SPACE) block in &
      &expval.data.')
     nex=nex+1 ; inex(nex)=16 ; cond_fraction_mom_in_expval=.true.
     write(6,*)'(',trim(i2s(nex)),') condensate fraction (momentum space)'
    endif
   enddo ! headers in the expval block
   write(6,*)

   expval_in_expval=.true.

  endif ! char_80

 enddo ! Look for header.

 rewind(io_expval)

 if(.not.expval_in_expval)then
  call errstop('FIND_EXPVAL_SETS','The expval.data file must contain a&
   &"START EXPVAL" block.')
 endif

 END SUBROUTINE find_expval_sets


 SUBROUTINE read_expval_block
!--------------------------------------------------------!
! Read basic info following START EXPVAL in expval.data. !
!--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,ierr

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_EXPVAL_BLOCK','Problem reading expval.data.')
  if(ierr<0)call errstop('READ_EXPVAL_BLOCK','Could not find &
   &"START EXPVAL" in expval.data.')
  if(trim(adjustl(char_80))=='START EXPVAL')exit
 enddo
 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)title
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)file_version
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)num_particle_types
 allocate(nele_expval(num_particle_types),stat=ialloc)
 if(ialloc/=0)call errstop('READ_EXPVAL_BLOCK','Error allocating nele_expval &
  &array')
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)(nele_expval(i),i=1,num_particle_types)
 netot=sum(nele_expval(:))
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)expval_dimensionality
 if(expval_dimensionality<1.or.expval_dimensionality>3) &
  &call errstop('READ_EXPVAL_BLOCK','Invalid dimensionality in expval.data. &
  &Must be 1, 2 or 3.')
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)expval_periodicity
 if(expval_periodicity<0.or.expval_periodicity>expval_dimensionality) &
  &call errstop('READ_EXPVAL_BLOCK','Invalid periodicity in expval.data.')
 read(io_expval,*,err=1,end=1)
 pa1=0.d0 ; pa2=0.d0 ; pa3=0.d0
 do j=1,3
  do i=1,3
   if(i==j)then
    expval_scell_matrix(i,j)=1
   else
    expval_scell_matrix(i,j)=0
   endif ! i=j
  enddo ! i
 enddo ! j
 select case(expval_periodicity)
 case(0)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)
 case(1)
  read(io_expval,*,err=1,end=1)pa1(1)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_scell_matrix(1,1)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)simcell_length
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)wigner_seitz_radius
  pa2=(/0.d0,1.d0,0.d0/)
  pa3=(/0.d0,0.d0,1.d0/)
 case(2)
  read(io_expval,*,err=1,end=1)(pa1(i),i=1,2)
  read(io_expval,*,err=1,end=1)(pa2(i),i=1,2)
  read(io_expval,*,err=1,end=1)
  read(io_expval,'(a)',err=1,end=1)char_80
! Strange order of supercell matrix elements for backwards compatibility
! - diagonal elements must be read first.
  read(char_80,*,iostat=ierr)expval_scell_matrix(1,1),expval_scell_matrix(2,2),&
   &expval_scell_matrix(1,2),expval_scell_matrix(2,1)
  if(ierr/=0)then
   expval_scell_matrix(1,2)=0 ; expval_scell_matrix(2,1)=0
   read(char_80,*,err=1,end=1)expval_scell_matrix(1,1),expval_scell_matrix(2,2)
  endif ! ierr/=0
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)simcell_area
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)wigner_seitz_radius
  pa3=(/0.d0,0.d0,1.d0/)
 case(3)
  read(io_expval,*,err=1,end=1)(pa1(i),i=1,3)
  read(io_expval,*,err=1,end=1)(pa2(i),i=1,3)
  read(io_expval,*,err=1,end=1)(pa3(i),i=1,3)
  read(io_expval,*,err=1,end=1)
  read(io_expval,'(a)',err=1,end=1)char_80
! Strange order of supercell matrix elements for backwards compatibility
! - diagonal elements must be read first.
  read(char_80,*,iostat=ierr)expval_scell_matrix(1,1),expval_scell_matrix(2,2),&
   &expval_scell_matrix(3,3),expval_scell_matrix(1,2),expval_scell_matrix(1,3),&
   &expval_scell_matrix(2,1),expval_scell_matrix(2,3),expval_scell_matrix(3,1),&
   &expval_scell_matrix(3,2)
  if(ierr/=0)then
   expval_scell_matrix=0
   read(char_80,*,err=1,end=1)expval_scell_matrix(1,1),&
    &expval_scell_matrix(2,2),expval_scell_matrix(3,3)
  endif ! ierr/=0
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)simcell_volume
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)wigner_seitz_radius
 case default
  call errstop('READ_EXPVAL_BLOCK','Invalid periodicity.')
 end select
 if(expval_periodicity>0)then
  pamat(1,:)=pa1(:) ; pamat(2,:)=pa2(:) ; pamat(3,:)=pa3(:)
  call inv_33(pamat,pbmat)
  pbmat(:,:)=transpose(pbmat(:,:))*twopi
  amat=matmul(dble(expval_scell_matrix),pamat)
  call inv_33(amat,bmat)
  bmat=transpose(bmat)*twopi
 endif
 npcells=abs(expval_scell_matrix(1,1)*(expval_scell_matrix(2,2)&
  &*expval_scell_matrix(3,3)-expval_scell_matrix(2,3)*expval_scell_matrix(3,2))&
  &+expval_scell_matrix(2,1)*(expval_scell_matrix(3,2)*expval_scell_matrix(1,3)&
  &-expval_scell_matrix(1,2)*expval_scell_matrix(3,3))+expval_scell_matrix(3,1)&
  &*(expval_scell_matrix(1,2)*expval_scell_matrix(2,3)-expval_scell_matrix(1,3)&
  &*expval_scell_matrix(2,2)))
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)expval_ngsets
 if(expval_ngsets>9)call errstop('READ_EXPVAL_BLOCK','No more than nine &
  &G vector sets may be defined in expval.data.')

 rewind(io_expval)

 return

1 call errstop('READ_EXPVAL_BLOCK','Problem reading START EXPVAL &
  &section of expval.data file.')
 END SUBROUTINE read_expval_block


 SUBROUTINE read_g_vector_sets
!------------------------------------------------!
! Read in GVECTOR SET n blocks from expval.data. !
!------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k
 CHARACTER(1) set_g
 REAL(dp),PARAMETER :: tol=1.d-8

 allocate(e_cutoff(expval_ngsets),expval_ngvec(expval_ngsets),&
  &expval_amat(3,3,expval_ngsets),expval_bmat(3,3,expval_ngsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_G_VECTOR_SETS','Allocation problem : &
  &e_cutoff,expval_ngvec,expval_amat,expval_bmat')

! Skip read to get expval_gvec array size.
 do i=1,expval_ngsets

  write(set_g,'(i1)')i

  do
   read(io_expval,'(a)',iostat=ierr)char_80
   if(ierr>0)call errstop('READ_G_VECTOR_SETS','Problem reading expval.data.')
   if(ierr<0)call errstop('READ_G_VECTOR_SETS','Could not find &
    &"START GVECTOR SET x" in expval.data for x = '//trim(i2s(i)))
   if(trim(adjustl(char_80))=='START GVECTOR SET '//set_g)exit
  enddo

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_ngvec(i)

 enddo ! G vector sets

 i=maxval(expval_ngvec(:))
 allocate(expval_gvec(3,i,expval_ngsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_G_VECTOR_SETS','Allocation problem : &
  &expval_gvec array.')

 rewind(io_expval)

! Detailed read
 do i=1,expval_ngsets

  write(set_g,'(i1)')i

  do
   read(io_expval,'(a)',iostat=ierr)char_80
   if(trim(adjustl(char_80))=='START GVECTOR SET '//set_g)exit
  enddo

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)e_cutoff(i)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_ngvec(i)
  read(io_expval,*,err=1,end=1)

  do j=1,3
   read(io_expval,*,err=1,end=1)(expval_bmat(j,k,i),k=1,3)
  enddo
  if(all(abs(expval_bmat(1:expval_periodicity,1:expval_periodicity,i)&
   &-pbmat(1:expval_periodicity,1:expval_periodicity))<tol))then
   expval_amat(1:3,1:3,i)=pamat
  elseif(all(abs(expval_bmat(1:expval_periodicity,1:expval_periodicity,i)&
   &-bmat(1:expval_periodicity,1:expval_periodicity))<tol))then
   expval_amat(1:3,1:3,i)=amat
  else
   call errstop('READ_G_VECTOR_SETS','Set of reciprocal lattice vectors in &
    &expval.data do not correspond to the primitive cell or supercell.')
  endif

  read(io_expval,*,err=1,end=1)
  do j=1,expval_ngvec(i)
   read(io_expval,*,err=1,end=1)(expval_gvec(k,j,i),k=1,3)
  enddo

 enddo ! G vector sets

 rewind(io_expval)

 return

1 call errstop('READ_G_VECTOR_SETS','Problem reading G vector &
  &sets in expval.data file.')
 END SUBROUTINE read_g_vector_sets


 SUBROUTINE read_density
!-------------------------------------------!
! Read in a DENSITY block from expval.data. !
!-------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j
 CHARACTER(2) nset
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_DENSITY','Problem reading expval.data file.')
  if(ierr<0)call errstop('READ_DENSITY','Could not find START DENSITY in &
   &expval.data file.')
  if(trim(adjustl(char_80))=='START DENSITY')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_den=trim(adjustl(char_80))
 if(accum_method_den/='VMC'.and.accum_method_den/='DMC') &
  &call errstop('READ_DENSITY','Density accumulation method (VMC or DMC) in &
  &expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)den_gset
 if(den_gset<1.or.den_gset>expval_ngsets)call errstop('READ_DENSITY',&
  &'Invalid specification of which G vector set was used.')
 gset(1)=den_gset

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)den_nsets

 allocate(expval_den_weight(den_nsets),den_ptype(den_nsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_DENSITY','Allocation error: &
  &expval_den_weight/den_ptype arrays')
 allocate(expval_den(expval_ngvec(den_gset),den_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_DENSITY','Allocation error: expval_den arrays')

 do i=1,den_nsets

  write(nset,'(i2)')i
  nset=trim(adjustl(nset))
  read(io_expval,'(a)',err=1,end=1)char_80

  if(trim(adjustl(char_80))/='START SET '//nset)then
   call errstop('READ_DENSITY','Error reading START SET x line for x = '//&
    &trim(i2s(i)))
  endif

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)den_ptype(i)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_den_weight(i)
  read(io_expval,*,err=1,end=1)

  do j=1,expval_ngvec(den_gset)
   read(io_expval,*,err=1,end=1)treal,timag
   expval_den(j,i)=cmplx(treal,timag,kind=dp)
  enddo

  read(io_expval,*,err=1,end=1)

 enddo ! density sets

 return

1 call errstop('READ_DENSITY','Problem reading DENSITY &
  &section of expval.data file')
 END SUBROUTINE read_density


 SUBROUTINE read_spin_density
!------------------------------------------------!
! Read in a SPIN DENSITY block from expval.data. !
!------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j
 CHARACTER(2) nset
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_SPIN_DENSITY','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_SPIN_DENSITY','Could not find &
   &START SPIN DENSITY in expval.data file.')
  if(trim(adjustl(char_80))=='START SPIN DENSITY')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_sden=trim(adjustl(char_80))
 if(accum_method_sden/='VMC'.and.accum_method_sden/='DMC')&
  &call errstop('READ_SPIN_DENSITY','Spin density accumulation method (VMC or &
  &DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sden_gset
 if(sden_gset<1.or.sden_gset>expval_ngsets)call errstop('READ_SPIN_DENSITY', &
  &'Invalid specification of which G vector set was used.')
 gset(2)=sden_gset

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sden_nsets

 allocate(expval_sden_weight(sden_nsets),sden_ptype(sden_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_SPIN_DENSITY','Allocation error: &
  &expval_sden_weight/sden_ptype arrays')
 allocate(expval_sden(expval_ngvec(sden_gset),sden_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_SPIN_DENSITY','Allocation error: &
  &expval_sden arrays')

 do i=1,sden_nsets

  write(nset,'(i2)')i
  nset=trim(adjustl(nset))
  read(io_expval,'(a)',err=1,end=1)char_80

  if(trim(adjustl(char_80))/='START SET '//nset)then
   call errstop('READ_SPIN_DENSITY','Error reading START SET x line for x = '&
    &//trim(i2s(i)))
  endif

  read(io_expval,*,err=1,end=1)

  if(nele_expval(i)>0)then

   read(io_expval,*,err=1,end=1)sden_ptype(i)
   read(io_expval,*,err=1,end=1)
   read(io_expval,*,err=1,end=1)expval_sden_weight(i)
   read(io_expval,*,err=1,end=1)

   do j=1,expval_ngvec(sden_gset)
    read(io_expval,*,err=1,end=1)treal,timag
    expval_sden(j,i)=cmplx(treal,timag,kind=dp)
   enddo

   read(io_expval,*,err=1,end=1)

  else

   read(io_expval,*,err=1,end=1)

  endif

 enddo ! spin density sets

 return

1 call errstop('READ_SPIN_DENSITY','Problem reading SPIN DENSITY &
  &section of expval.data file')
 END SUBROUTINE read_spin_density


 SUBROUTINE read_spin_density_matrix
!-------------------------------------------------------!
! Read in a SPIN-DENSITY MATRIX block from expval.data. !
!-------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k
 CHARACTER(2) nset
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_SPIN_DENSITY_MATRIX','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_SPIN_DENSITY_MATRIX','Could not find &
   &START SPIN-DENSITY MATRIX in expval.data file.')
  if(trim(adjustl(char_80))=='START SPIN-DENSITY MATRIX')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_sdenm=trim(adjustl(char_80))
 if(accum_method_sdenm/='VMC'.and.accum_method_sdenm/='DMC')&
  &call errstop('READ_SPIN_DENSITY_MATRIX','Spin density matrix accumulation &
  &method (VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sdenm_gset
 if(sdenm_gset<1.or.sdenm_gset>expval_ngsets)call errstop&
  &('READ_SPIN_DENSITY_MATRIX','Invalid specification of which G vector set to&
  & use.')
 gset(3)=sdenm_gset

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sdenm_nsets

 allocate(expval_sdenm_weight(sdenm_nsets),sdenm_ptype(sdenm_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_SPIN_DENSITY_MATRIX','Allocation error: &
  &expval_sdenm_weight/sdenm_ptype arrays')
 allocate(expval_sdenm(expval_ngvec(sdenm_gset),4,sdenm_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_SPIN_DENSITY_MATRIX',&
  &'Allocation error: expval_sdenm array')

 do i=1,sdenm_nsets

  write(nset,'(i2)')i
  nset=trim(adjustl(nset))
  read(io_expval,'(a)',err=1,end=1)char_80

  if(trim(adjustl(char_80))/='START SET '//nset)then
   call errstop('READ_SPIN_DENSITY_MATRIX','Error reading START SET x line &
    &for x = '//trim(i2s(i)))
  endif
  read(io_expval,'(a)',err=1,end=1)char_80

  if(nele_expval(i)>0)then

   if(trim(adjustl(char_80))=="No particles of this type.")call &
    &errstop('READ_SPIN_DENSITY_MATRIX','Error: no data provided for a set &
    &that should be non-empty.')

   read(io_expval,*,err=1,end=1)sdenm_ptype(i)
   read(io_expval,*,err=1,end=1)
   read(io_expval,*,err=1,end=1)expval_sdenm_weight(i)
   read(io_expval,*,err=1,end=1)

   do k=1,4
    read(io_expval,*,err=1,end=1)
    read(io_expval,'(a)',err=1,end=1)char_80
    select case(k)
    case(1)
     if(trim(adjustl(char_80))/='1 1')call errstop('READ_SPIN_DENSITY_MATRIX',&
      &'Error reading component 1 1.')
    case(2)
     if(trim(adjustl(char_80))/='1 2')call errstop('READ_SPIN_DENSITY_MATRIX',&
      &'Error reading component 1 2.')
    case(3)
     if(trim(adjustl(char_80))/='2 1')call errstop('READ_SPIN_DENSITY_MATRIX',&
      &'Error reading component 2 1.')
    case(4)
     if(trim(adjustl(char_80))/='2 2')call errstop('READ_SPIN_DENSITY_MATRIX',&
      &'Error reading component 2 2.')
    end select
    do j=1,expval_ngvec(sdenm_gset)
     read(io_expval,*,err=1,end=1)treal,timag
     expval_sdenm(j,k,i)=cmplx(treal,timag,kind=dp)
    enddo ! j
   enddo ! k

   read(io_expval,*,err=1,end=1)

  else

   if(trim(adjustl(char_80))/="No particles of this type.")call &
    &errstop('READ_SPIN_DENSITY_MATRIX','Error: data provided for a set &
    &that should be empty.')

   read(io_expval,*,err=1,end=1)

  endif

 enddo ! spin density matrix sets

 return

1 call errstop('READ_SPIN_DENSITY_MATRIX','Problem reading SPIN-DENSITY MATRIX&
  & section of expval.data file')
 END SUBROUTINE read_spin_density_matrix


 SUBROUTINE read_pair_corr
!--------------------------------------------------------!
! Read in a RECIPROCAL-SPACE PCF block from expval.data. !
!--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j
 CHARACTER(2) nset
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_PAIR_CORR','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_PAIR_CORR','Could not find &
   &START RECIPROCAL-SPACE PCF in expval.data file.')
  if(trim(adjustl(char_80))=='START RECIPROCAL-SPACE PCF')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_pcf=trim(adjustl(char_80))
 if(accum_method_pcf/='VMC'.and.accum_method_pcf/='DMC') &
  &call errstop('READ_PAIR_CORR','Reciprocal-space PCF accumulation method &
  &(VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)pcf_type_fixed
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)(pcf_rfix(i),i=1,3)

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)pcf_gset
 if(pcf_gset<1.or.pcf_gset>expval_ngsets)call errstop('READ_PAIR_CORR', &
  &'Invalid specification of which G vector set to use for rec. space PCF.')
 gset(4)=pcf_gset

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)pcf_nsets

 allocate(expval_pcf_weight(pcf_nsets),pcf_ptype(pcf_nsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_PAIR_CORR','Allocation error: &
  &expval_pcf_weight/pcf_ptype arrays')
 allocate(expval_pcf(expval_ngvec(pcf_gset),pcf_nsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_PAIR_CORR','Allocation error: expval_pcf&
  & array')

 do i=1,pcf_nsets

  write(nset,'(i2)')i
  nset=trim(adjustl(nset))
  read(io_expval,'(a)',err=1,end=1)char_80

  if(trim(adjustl(char_80))/='START SET '//nset)call errstop('READ_PAIR_CORR',&
   &'Error reading START SET x line for x = '//trim(i2s(i)))

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)pcf_ptype(i)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_pcf_weight(i)
  read(io_expval,*,err=1,end=1)

  do j=1,expval_ngvec(pcf_gset)
   read(io_expval,*,err=1,end=1)treal,timag
   expval_pcf(j,i)=cmplx(treal,timag,kind=dp)
  enddo

  read(io_expval,*,err=1,end=1)

 enddo ! spin density sets

 return

1 call errstop('READ_PAIR_CORR','Problem reading RECIPROCAL-SPACE PCF &
  &section of expval.data file')

 END SUBROUTINE read_pair_corr


 SUBROUTINE read_pair_corr_sph
!-------------------------------------------------!
! Read in a SPHERICAL PCF block from expval.data. !
!-------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,set
 CHARACTER(2) nset
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_PAIR_CORR_SPH','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_PAIR_CORR_SPH','Could not find &
   &START SPHERICAL PCF in expval.data file.')
  if(trim(adjustl(char_80))=='START SPHERICAL PCF')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_pcfs=trim(adjustl(char_80))
 if(accum_method_pcfs/='VMC'.and.accum_method_pcfs/='DMC') &
  &call errstop('READ_PAIR_CORR_SPH','Spherical PCF accumulation method (VMC&
  & or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)pcfs_nbins
 if(pcfs_nbins<1)call errstop('READ_PAIR_CORR_SPH',&
  &'Need at least 1 bin for accumulation of spherical PCF in expval.data.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)pcfs_rcutoff
 if(pcfs_rcutoff<=0.d0)call errstop('READ_PAIR_CORR_SPH',&
  &'Invalid cutoff radius for spherical PCF in expval.data.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)pcfs_accmode
 if(pcfs_accmode/=1.and.pcfs_accmode/=2)call errstop('READ_PAIR_CORR_SPH',&
  &'Invalid accumulation mode for spherical PCF in expval.data.')

 if(pcfs_accmode==1)then
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)pcfs_type_fixed
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)(pcfs_rfix(i),i=1,3)
 endif

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)pcfs_nsets

 allocate(expval_pcfs_weight(pcfs_nsets),pcfs_ptype(2,pcfs_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_PAIR_CORR_SPH','Allocation error: &
  &expval_pcfs_weight/pcfs_ptype')
 if(pcfs_accmode==1)then
  allocate(expval_pcfs(pcfs_nbins,1,num_particle_types),stat=ialloc)
 else
  allocate(expval_pcfs(pcfs_nbins,num_particle_types,num_particle_types),&
   &stat=ialloc)
 endif
 if(ialloc/=0)call errstop('READ_PAIR_CORR_SPH','Allocation error: expval_pcfs &
  &arrays')

 do set=1,pcfs_nsets

  write(nset,'(i2)')set
  nset=trim(adjustl(nset))
  read(io_expval,'(a)',err=1,end=1)char_80

  if(trim(adjustl(char_80))/='START SET '//nset)call errstop&
   &('READ_PAIR_CORR_SPH','Error reading START SET x line for x = '//&
   &trim(i2s(set)))

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)pcfs_ptype(1,set),pcfs_ptype(2,set)
  if(pcfs_accmode==1)then
   if(pcfs_ptype(1,set)/=pcfs_type_fixed)then
    write(6,*)
    write(6,*)'First particle type in pair must equal the type of the fixed'
    write(6,*)'particle when accumulating the spherical PCF in inhomogeneous'
    write(6,*)'systems.'
    call errstop('READ_PAIR_CORR_SPH','Quitting.')
   endif
  endif
  if(any(pcfs_ptype(:,set)<1))then
   call errstop('READ_PAIR_CORR_SPH','Types of particle in pair are out of&
    & range in spherical PCF block in expval.data.')
  endif

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_pcfs_weight(set)
  read(io_expval,*,err=1,end=1)

  j=pcfs_ptype(1,set) ; k=pcfs_ptype(2,set)
  do i=1,pcfs_nbins
   read(io_expval,*,err=1,end=1)expval_pcfs(i,j,k)
  enddo

  read(io_expval,*,err=1,end=1)

 enddo ! spherical PCF sets

 return

1 call errstop('READ_PAIR_CORR_SPH','Problem reading SPHERICAL PCF &
  &section of expval.data file')
 END SUBROUTINE read_pair_corr_sph


 SUBROUTINE read_loc_tensor
!-------------------------------------------------------!
! Read in a LOCALIZATION TENSOR block from expval.data. !
!-------------------------------------------------------!
 IMPLICIT NONE
 call errstop('READ_LOC_TENSOR','Routine not yet coded')
 END SUBROUTINE read_loc_tensor


 SUBROUTINE read_structure_factor
!----------------------------------------------------!
! Read in a STRUCTURE FACTOR block from expval.data. !
!----------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,k,set
 CHARACTER(2) nset
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_STRUCTURE_FACTOR','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_STRUCTURE_FACTOR','Could not find &
   &START STRUCTURE FACTOR in expval.data file.')
  if(trim(adjustl(char_80))=='START STRUCTURE FACTOR')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_sf=trim(adjustl(char_80))
 if(accum_method_sf/='VMC'.and.accum_method_sf/='DMC')&
  &call errstop('READ_STRUCTURE_FACTOR','Structure factor accumulation method &
   &(VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sf_gset

 if(sf_gset<1.or.sf_gset>expval_ngsets)call errstop('READ_STRUCTURE_FACTOR',&
  &'Invalid specification of which G vector set to use.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sf_nsets

 allocate(sf_ptype(2,sf_nsets),expval_sf(expval_ngvec(sf_gset),sf_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_STRUCTURE_FACTOR','Allocation error: &
  &sf_ptype/expval_sf arrays')

 do set=1,sf_nsets

  write(nset,'(i2)')set
  nset=trim(adjustl(nset))
  read(io_expval,'(a)',err=1,end=1)char_80

  if(trim(adjustl(char_80))/='START SET '//nset)call errstop&
   &('READ_STRUCTURE_FACTOR','Error reading START SET x line for x = '//&
   &trim(i2s(set)))

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)sf_ptype(1,set),sf_ptype(2,set)

  if(nele_expval(sf_ptype(1,set))/=0.and.nele_expval(sf_ptype(2,set))/=0)then
   read(io_expval,*,err=1,end=1)
   read(io_expval,*,err=1,end=1)expval_sf_weight
   read(io_expval,*,err=1,end=1)
   do k=1,expval_ngvec(sf_gset)
    read(io_expval,*,err=1,end=1)expval_sf(k,set)
   enddo
  else
   read(io_expval,*,err=1,end=1)
  endif

  read(io_expval,*,err=1,end=1)

 enddo ! sets

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sfsden_nsets

 if(sfsden_nsets==0)then
  homogeneous_structure_factor=.true.
 else
  homogeneous_structure_factor=.false.
 endif

 allocate(sfsden_ptype(sfsden_nsets),expval_sfsden(expval_ngvec(sf_gset),&
  &sfsden_nsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_STRUCTURE_FACTOR','Allocation error: &
  &sfsden_ptype/expval_sfsden arrays')

 do set=1,sfsden_nsets

  write(nset,'(i2)')set
  nset=trim(adjustl(nset))
  read(io_expval,'(a)',err=1,end=1)char_80

  if(trim(adjustl(char_80))/='START SET '//nset)then
   call errstop('READ_STRUCTURE_FACTOR','Error reading START SET x line &
    &for spin density part of structure factor for x = '//trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1)

  if(nele_expval(set)>0)then

   read(io_expval,*,err=1,end=1)sfsden_ptype(set)
   read(io_expval,*,err=1,end=1)
   read(io_expval,*,err=1,end=1)expval_sfsden_weight
   read(io_expval,*,err=1,end=1)

   do i=1,expval_ngvec(sf_gset)
    read(io_expval,*,err=1,end=1)treal,timag
    expval_sfsden(i,set)=cmplx(treal,timag,kind=dp)
   enddo

   read(io_expval,*,err=1,end=1)

  else

   read(io_expval,*,err=1,end=1)

  endif

 enddo ! spin density for structure factor sets

 return

1 call errstop('READ_STRUCTURE_FACTOR','Problem reading STRUCTURE FACTOR &
   &section of expval.data file.')

 END SUBROUTINE read_structure_factor


 SUBROUTINE read_structure_factor_sph
!--------------------------------------------------------------!
! Read in a SPHERICAL STRUCTURE FACTOR block from expval.data. !
!--------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER k,set
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_STRUCTURE_FACTOR_SPH','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_STRUCTURE_FACTOR_SPH','Could not find &
   &START SPHERICAL STRUCTURE FACTOR in expval.data file.')
  if(trim(adjustl(char_80))=='START SPHERICAL STRUCTURE FACTOR')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_sf_sph=trim(adjustl(char_80))
 if(accum_method_sf_sph/='VMC'.and.accum_method_sf_sph/='DMC')&
  &call errstop('READ_STRUCTURE_FACTOR_SPH','Spherical structure factor &
  &accumulation method (VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sf_sph_k1,sf_sph_k2,sf_sph_nk

 allocate(expval_sf_sph_k(sf_sph_nk),stat=ialloc)
 if(ialloc/=0)call errstop('READ_STRUCTURE_FACTOR_SPH','Allocation error: &
  &expval_sf_sph_k')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)sf_sph_nsets

 allocate(expval_sf_sph_weight(sf_sph_nsets),sf_sph_ptype(2,sf_sph_nsets),&
  &expval_sf_sph(sf_sph_nk,sf_sph_nsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_STRUCTURE_FACTOR_SPH','Allocation error: &
  &expval_sf_sph arrays')

 do set=1,sf_sph_nsets

  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START SET '//trim(i2s(set)))call &
   &errstop('READ_STRUCTURE_FACTOR_SPH','Error reading START SET '&
   &//trim(i2s(set))//'.')

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)sf_sph_ptype(1:2,set)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_sf_sph_weight(set)
  read(io_expval,*,err=1,end=1)

  do k=1,sf_sph_nk
   read(io_expval,*,err=1,end=1)expval_sf_sph_k(k),expval_sf_sph(k,set)
  enddo ! k

  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='END SET '//trim(i2s(set)))call &
   &errstop('READ_STRUCTURE_FACTOR_SPH','Error reading END SET '&
   &//trim(i2s(set))//'.')

 enddo ! sets

 return

1 call errstop('READ_STRUCTURE_FACTOR_SPH','Problem reading SPHERICAL &
   &STRUCTURE FACTOR section of expval.data file.')

 END SUBROUTINE read_structure_factor_sph


 SUBROUTINE read_onep_density_matrix
 IMPLICIT NONE
 INTEGER i,set,idum,idum1,idum2
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_ONEP_DENSITY_MATRIX','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_ONEP_DENSITY_MATRIX','Could not find &
   &START ONE-PARTICLE DENSITY MATRIX in expval.data file.')
  if(trim(adjustl(char_80))=='START ONE-PARTICLE DENSITY MATRIX')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_onep_dm=trim(adjustl(char_80))
 if(accum_method_onep_dm/='VMC'.and.accum_method_onep_dm/='DMC')&
  &call errstop('READ_ONEP_DENSITY_MATRIX','One-particle density matrix &
  &accumulation method (VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)onep_dm_nsets
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)onep_dm_nbins
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)idum

 allocate(expval_onep_dm_weight(onep_dm_nbins,onep_dm_nsets),&
  &expval_onep_dm(onep_dm_nbins,onep_dm_nsets),&
  &expval_onep_dm2(onep_dm_nbins,onep_dm_nsets),&
  &onep_dm_nptypes_in_set(onep_dm_nsets),&
  &onep_dm_ptype_in_set(num_particle_types,onep_dm_nsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_ONEP_DENSITY_MATRIX','Allocation error: &
  &expval_onep_dm arrays')
 expval_onep_dm_weight=0.d0 ; expval_onep_dm=0.d0
 expval_onep_dm2=0.d0 ; onep_dm_nptypes_in_set=0
 onep_dm_ptype_in_set=0

 do set=1,onep_dm_nsets

  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START SET '//trim(i2s(set)))then
   call errstop('READ_ONEP_DENSITY_MATRIX','Error reading START SET x line &
    &for x = '//trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1)
  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))=='on-top pair')then
   onep_dm_nptypes_in_set(set)=1
   read(io_expval,*,err=1,end=1)
   read(io_expval,*,err=1,end=1)idum1,idum2
   if(idum1<1.or.idum1>num_particle_types.or.idum2<1.or.&
    &idum2>num_particle_types)call errstop('READ_ONEP_DENSITY_MATRIX',&
    &'Particle types for on-top density matrix out of range.')
   if(on_top_ispin/=0.and.(idum1/=on_top_ispin.or.idum2/=on_top_jspin))&
    &call errstop('READ_ONEP_DENSITY_MATRIX','Particle types for on-top &
    &density matrix differ from previously-read on-top particle types.')
   on_top_ispin=idum1
   on_top_jspin=idum2
   onep_dm_ptype_in_set(1,set)=on_top_ispin
  else
   read(char_80,*,err=1,end=1)onep_dm_nptypes_in_set(set)
   read(io_expval,*,err=1,end=1)
   read(io_expval,*,err=1,end=1)&
    &onep_dm_ptype_in_set(1:onep_dm_nptypes_in_set(set),set)
  endif

  read(io_expval,*,err=1,end=1)
  do i=1,onep_dm_nbins
   read(io_expval,*,err=1,end=1)expval_onep_dm_weight(i,set),&
    &expval_onep_dm2(i,set),expval_onep_dm(i,set)
  enddo

  read(io_expval,*,err=1,end=1) ! END SET x

 enddo ! sets

 return

1 call errstop('READ_ONEP_DENSITY_MATRIX','Problem reading ONE-PARTICLE &
   &DENSITY MATRIX section of expval.data file.')

 END SUBROUTINE read_onep_density_matrix


 SUBROUTINE read_twop_density_matrix
 IMPLICIT NONE
 INTEGER i,set
 REAL(dp) dumr
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_TWOP_DENSITY_MATRIX','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_TWOP_DENSITY_MATRIX','Could not find &
   &START TWO-PARTICLE DENSITY MATRIX in expval.data file.')
  if(trim(adjustl(char_80))=='START TWO-PARTICLE DENSITY MATRIX')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_twop_dm=trim(adjustl(char_80))
 if(accum_method_twop_dm/='VMC'.and.accum_method_twop_dm/='DMC')&
  &call errstop('READ_TWOP_DENSITY_MATRIX','One-particle density matrix &
  &accumulation method (VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)twop_dm_nsets
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)twop_dm_nbins
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)dumr
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)dumr

 allocate(expval_twop_dm_weight(twop_dm_nbins,twop_dm_nsets),&
  &expval_twop_dm(twop_dm_nbins,twop_dm_nsets),&
  &expval_twop_dm2(twop_dm_nbins,twop_dm_nsets),&
  &twop_dm_nptypes_in_set(twop_dm_nsets),&
  &twop_dm_ptype_in_set(2,((num_particle_types+1)*num_particle_types)/2,&
  &twop_dm_nsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_TWOP_DENSITY_MATRIX','Allocation error: &
  &expval_twop_dm arrays')
 expval_twop_dm_weight=0.d0 ; expval_twop_dm=0.d0
 expval_twop_dm2=0.d0 ; twop_dm_nptypes_in_set=0
 twop_dm_ptype_in_set=0

 do set=1,twop_dm_nsets

  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START SET '//trim(i2s(set)))then
   call errstop('READ_TWOP_DENSITY_MATRIX','Error reading START SET x line &
    &for x = '//trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)twop_dm_nptypes_in_set(set)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)&
   &twop_dm_ptype_in_set(1:2,1:twop_dm_nptypes_in_set(set),set)
  read(io_expval,*,err=1,end=1)

  do i=1,twop_dm_nbins
   read(io_expval,*,err=1,end=1)expval_twop_dm_weight(i,set),&
    &expval_twop_dm2(i,set),expval_twop_dm(i,set)
  enddo

  read(io_expval,*,err=1,end=1) ! END SET x

 enddo ! sets

 return

1 call errstop('READ_TWOP_DENSITY_MATRIX','Problem reading TWO-PARTICLE &
   &DENSITY MATRIX section of expval.data file.')

 END SUBROUTINE read_twop_density_matrix


 SUBROUTINE read_cond_fraction
 IMPLICIT NONE
 INTEGER i,set
 REAL(dp) dumr
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_COND_FRACTION','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_COND_FRACTION','Could not find &
   &START CONDENSATE FRACTION in expval.data file.')
  if(trim(adjustl(char_80))=='START CONDENSATE FRACTION')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_cond_frac=trim(adjustl(char_80))
 if(accum_method_cond_frac/='VMC'.and.&
  &accum_method_cond_frac/='DMC')call errstop('READ_COND_FRACTION',&
  &'Condensate fraction estimator accumulation method (VMC or DMC) in &
  &expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)cond_frac_nsets
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)cond_frac_nbins
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)dumr
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)dumr

 allocate(expval_cond_frac_weight(cond_frac_nbins,cond_frac_nsets),&
  &expval_cond_frac(cond_frac_nbins,cond_frac_nsets),&
  &expval_cond_frac2(cond_frac_nbins,cond_frac_nsets),&
  &cond_frac_nptypes_in_set(cond_frac_nsets),&
  &cond_frac_ptype_in_set(2,((num_particle_types+1)*num_particle_types)/2,&
  &cond_frac_nsets),stat=ialloc)
 if(ialloc/=0)call errstop('READ_COND_FRACTION','Allocation error: &
  &expval_cond_frac arrays')
 expval_cond_frac_weight=0.d0 ; expval_cond_frac=0.d0
 cond_frac_nptypes_in_set=0 ; cond_frac_ptype_in_set=0

 do set=1,cond_frac_nsets

  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START SET '//trim(i2s(set)))then
   call errstop('READ_COND_FRACTION','Error reading START SET x line &
    &for x = '//trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)cond_frac_nptypes_in_set(set)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)&
   &cond_frac_ptype_in_set(1:2,1:cond_frac_nptypes_in_set(set),set)
  read(io_expval,*,err=1,end=1)

  do i=1,cond_frac_nbins
   read(io_expval,*,err=1,end=1)expval_cond_frac_weight(i,set),&
    &expval_cond_frac2(i,set),expval_cond_frac(i,set)
  enddo

  read(io_expval,*,err=1,end=1) ! END SET x

 enddo ! sets

 return

1 call errstop('READ_COND_FRACTION','Problem reading CONDENSATE &
   &FRACTION section of expval.data file.')

 END SUBROUTINE read_cond_fraction


 SUBROUTINE read_mom_den
 IMPLICIT NONE
 INTEGER i,set,idum,idum1,idum2
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_MOM_DEN','Problem reading expval.data file.')
  if(ierr<0)call errstop('READ_MOM_DEN','Could not find START MOMENTUM &
   &DENSITY in expval.data file.')
  if(trim(adjustl(char_80))=='START MOMENTUM DENSITY')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_mom_den=trim(adjustl(char_80))
 if(accum_method_mom_den/='VMC'.and.accum_method_mom_den/='DMC')&
  &call errstop('READ_MOM_DEN','Momentum density accumulation method (VMC or &
  &DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)mom_den_nsets
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)idum

 allocate(expval_mom_den_weight(mom_den_nsets),&
  &expval_mom_den(expval_ngvec(1),mom_den_nsets),&
  &expval_mom_den2(expval_ngvec(1),mom_den_nsets),&
  &expval_mom_den_k(expval_ngvec(1),mom_den_nsets),&
  &mom_den_nptypes_in_set(mom_den_nsets),&
  &mom_den_ptype_in_set(num_particle_types,mom_den_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_MOM_DEN','Allocation error: expval_mom_den &
  &arrays')
 expval_mom_den_weight=0.d0 ; expval_mom_den=0.d0
 expval_mom_den2=0.d0 ; expval_mom_den_k=0.d0
 mom_den_nptypes_in_set=0 ; mom_den_ptype_in_set=0

 do set=1,mom_den_nsets

  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START SET '//trim(i2s(set)))then
   call errstop('READ_MOM_DEN','Error reading START SET x line for x = '//&
    &trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1)
  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))=='on-top pair')then
   mom_den_nptypes_in_set(set)=1
   read(io_expval,*,err=1,end=1)
   read(io_expval,*,err=1,end=1)idum1,idum2
   if(idum1<1.or.idum1>num_particle_types.or.idum2<1.or.&
    &idum2>num_particle_types)call errstop('READ_MOM_DEN','Particle types for &
    &on-top density matrix out of range.')
   if(on_top_ispin/=0.and.(idum1/=on_top_ispin.or.idum2/=on_top_jspin))&
    &call errstop('READ_MOM_DEN','Particle types for on-top density matrix &
    &differ from previously-read on-top particle types.')
   on_top_ispin=idum1
   on_top_jspin=idum2
   mom_den_ptype_in_set(1,set)=on_top_ispin
  else
   read(char_80,*,err=1)mom_den_nptypes_in_set(set)
   read(io_expval,*,err=1,end=1)
   read(io_expval,*,err=1,end=1)&
    &mom_den_ptype_in_set(1:mom_den_nptypes_in_set(set),set)
  endif
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_mom_den_weight(set)
  read(io_expval,*,err=1,end=1)

  do i=1,expval_ngvec(1)
   read(io_expval,*,err=1,end=1)expval_mom_den_k(i,set),&
    &expval_mom_den2(i,set),expval_mom_den(i,set)
  enddo

  read(io_expval,*,err=1,end=1) ! END SET x

 enddo ! sets

 return

1 call errstop('READ_MOM_DEN','Problem reading MOMENTUM DENSITY section of &
  &expval.data file.')

 END SUBROUTINE read_mom_den


 SUBROUTINE read_finite_density
!--------------------------------------------------!
! Read in a FINITE DENSITY block from expval.data. !
!--------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,set
 CHARACTER(80) char_80

! Set defaults.
 fin_den_nsets_ij=3 ! uu,ud,dd
 fin_den_nsets_i=2 ! u,d
 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_FINITE_DENSITY','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_FINITE_DENSITY','Could not find &
   &START FINITE DENSITY in expval.data file.')
  if(trim(adjustl(char_80))=='START FINITE DENSITY')exit
 enddo

 read(io_expval,*,err=1,end=1) ! accumulation method
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_fin_den=trim(adjustl(char_80))
 if(accum_method_fin_den/='VMC'.and.accum_method_fin_den/='DMC')&
  &call errstop('READ_FINITE_DENSITY','Finite density &
  &accumulation method (VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1) ! basis
 read(io_expval,*,err=1,end=1)fin_den_basis
 select case(trim(adjustl(fin_den_basis)))
 case('step')
  fin_den_ibasis=fin_den_f_indx_step
 case('sph_bessel')
  fin_den_ibasis=fin_den_f_indx_bessel
  call errstop('READ_FINITE_DENSITY','Full functionality &
   &is not available for a Bessel function basis.')
 case('chebyshev')
  fin_den_ibasis=fin_den_f_indx_cheb
  call errstop('READ_FINITE_DENSITY','Full functionality &
   &is not available for a Chebyshev polynomial basis.')
 case('exponential')
  fin_den_ibasis=fin_den_f_indx_exp
  call errstop('READ_FINITE_DENSITY','Full functionality &
   &is not available for an exponential function basis.')
 case default
  call errstop('READ_FINITE_DENSITY','Finite density &
   &basis function in expval.data not recognized.')
 end select
 read(io_expval,*,err=1,end=1) ! expansion order
 read(io_expval,*,err=1,end=1)fin_den_norder
 if(fin_den_norder<2)call errstop('READ_FINITE_DENSITY','Need at least &
  &two bins for accumulation of finite density in expval.data.')
 select case(fin_den_ibasis)
 case(fin_den_f_indx_step,fin_den_f_indx_cheb)
  read(io_expval,*,err=1,end=1) ! cutoff
  read(io_expval,*,err=1,end=1)fin_den_cutoff
  fin_den_tail_offset=fin_den_cutoff
 case(fin_den_f_indx_exp)
  read(io_expval,*,err=1,end=1) ! tail offset
  read(io_expval,*,err=1,end=1)fin_den_tail_offset
  fin_den_cutoff=fin_den_tail_offset
 case(fin_den_f_indx_bessel)
  read(io_expval,*,err=1,end=1) ! cutoff, tail offset
  read(io_expval,*,err=1,end=1)fin_den_cutoff,fin_den_tail_offset
 end select

 allocate(&
  &expval_fin_den_ij(fin_den_norder,fin_den_nsets_ij),&
  &expval_fin_den2_ij(fin_den_norder,fin_den_nsets_ij),&
  &expval_fin_den_i(fin_den_norder,fin_den_nsets_i),&
  &expval_fin_den2_i(fin_den_norder,fin_den_nsets_i),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_FINITE_DENSITY','Allocation error: &
  &expval_fin_den arrays')
 expval_fin_den_ij=0.d0 ; expval_fin_den_i=0.d0
 expval_fin_den2_ij=0.d0 ; expval_fin_den2_i=0.d0
 expval_fin_den_nstep=0.d0 ; expval_fin_den_weight=0.d0
 expval_fin_den_weight2=0.d0

 read(io_expval,'(a)',err=1,end=1)char_80 ! Nstep, Total weight, Total weight^2
 read(io_expval,*,err=1,end=1)expval_fin_den_nstep,expval_fin_den_weight,&
  &expval_fin_den_weight2

! e-e READ
loop_sets_ij: do set=1,fin_den_nsets_ij
  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START e-e SET '//trim(i2s(set)))then
   call errstop('READ_FINITE_DENSITY','Error reading START e-e SET x line &
    &for x = '//trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1) ! n_ij, n_ij^2
  do i=1,fin_den_norder
   read(io_expval,'(a)',err=1,end=1)char_80
   read(char_80,*,err=1,end=1)&
    &expval_fin_den_ij(i,set),expval_fin_den2_ij(i,set)
  enddo ! order
  read(io_expval,*,err=1,end=1) ! END e-e SET x
 enddo loop_sets_ij ! sets_ij

! e-nucleus READ
loop_sets_i: do set=1,fin_den_nsets_i
  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START e-nucleus SET '//trim(i2s(set)))then
   call errstop('READ_FINITE_DENSITY','Error reading START e-nucleus SET x &
    &line for x = '//trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1) ! n_i, n_i^2
  do i=1,fin_den_norder
   read(io_expval,'(a)',err=1,end=1)char_80
   read(char_80,*,err=1,end=1)&
    &expval_fin_den_i(i,set),expval_fin_den2_i(i,set)
  enddo ! order
  read(io_expval,*,err=1,end=1) ! END e-nucleus SET x
 enddo loop_sets_i ! sets_i

 return

1 call errstop('READ_FINITE_DENSITY','Problem reading FINITE DENSITY &
   &section of expval.data file.')

 END SUBROUTINE read_finite_density


 SUBROUTINE read_mol_density
!--------------------------------------------------!
! Read in a FINITE DENSITY block from expval.data. !
!--------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,set
 CHARACTER(80) char_80

! Set defaults.
 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_MOL_DENSITY','Problem reading &
   &expval.data file.')
  if(ierr<0)call errstop('READ_MOL_DENSITY','Could not find &
   &START MOLECULAR DENSITY in expval.data file.')
  if(trim(adjustl(char_80))=='START MOLECULAR DENSITY')exit
 enddo

 read(io_expval,*,err=1,end=1) ! accumulation method
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_mol_den=trim(adjustl(char_80))
 if(accum_method_mol_den/='VMC'.and.accum_method_mol_den/='DMC')&
  &call errstop('READ_MOL_DENSITY','Molecular density &
  &accumulation method (VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1) ! expansion order
 read(io_expval,*,err=1,end=1)mol_den_norder(1:3)
 if(product(mol_den_norder)<2)call errstop('READ_MOL_DENSITY','Need at least &
  &two bins for accumulation of molecular density in expval.data.')
 read(io_expval,*,err=1,end=1) ! coordinates of A (min corner)
 read(io_expval,*,err=1,end=1)mol_den_maxcoord(1:3)
 read(io_expval,*,err=1,end=1) ! coordinates of B (max corner)
 read(io_expval,*,err=1,end=1)mol_den_mincoord(1:3)
 if(.not.all(mol_den_mincoord<=mol_den_maxcoord))call errstop&
  &('READ_MOL_DENSITY','All components of point A must be <= those of &
  &point B.')
 read(io_expval,*,err=1,end=1) ! Nstep, Total weight
 read(io_expval,*,err=1,end=1)expval_mol_den_nstep,expval_mol_den_weight

! Count number of sets.
 mol_den_nsets=0
 set=0
 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr/=0)call errstop('READ_MOL_DENSITY','Error reading expval.data.')
  if(trim(adjustl(char_80))=='START SET '//trim(i2s(set+1)))then
   if(mol_den_nsets/=set)call errstop('READ_MOL_DENSITY','START/END SET &
    &mismatch in expval.data (1).')
   set=set+1
  elseif(trim(adjustl(char_80))=='END SET '//trim(i2s(set)))then
   if(mol_den_nsets+1/=set)call errstop('READ_MOL_DENSITY','START/END SET &
    &mismatch in expval.data (2).')
   mol_den_nsets=set
  elseif(trim(adjustl(char_80))=='END MOLECULAR DENSITY')then
   if(mol_den_nsets/=set)call errstop('READ_MOL_DENSITY','START/END SET &
    &mismatch in expval.data (3).')
   exit
  endif
 enddo
 if(mol_den_nsets<1)call errstop('READ_MOL_DENSITY','No data found in &
  &expval.data.')

! Go back to where we were.
 rewind(io_expval)
 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr/=0)call errstop('READ_MOL_DENSITY','Rewind paradox.')
  if(trim(adjustl(char_80))=='START MOLECULAR DENSITY')exit
 enddo
 read(io_expval,*,err=1,end=1) ! accumulation method
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1) ! expansion order
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1) ! coordinates of A (max corner)
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1) ! coordinates of B (min corner)
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1) ! Nstep, Total weight
 read(io_expval,*,err=1,end=1)

! Allocate work arrays.
 allocate(&
  &expval_mol_den(mol_den_norder(1),mol_den_norder(2),mol_den_norder(3),&
  &               mol_den_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_MOL_DENSITY','Allocation error: &
  &expval_mol_den arrays')
 expval_mol_den=0.d0
 expval_mol_den_nstep=0.d0
 expval_mol_den_weight=0.d0

! e-nucleus READ
 loop_sets_i: do set=1,mol_den_nsets
  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START SET '//trim(i2s(set)))then
   call errstop('READ_MOL_DENSITY','Error reading START SET x &
    &line for x = '//trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1) ! n_i
  do k=1,mol_den_norder(3)
   do j=1,mol_den_norder(2)
    do i=1,mol_den_norder(1)
     read(io_expval,'(a)',err=1,end=1)char_80
     read(char_80,*,err=1,end=1)expval_mol_den(i,j,k,set)
    enddo ! i
   enddo ! j
  enddo ! k
  read(io_expval,*,err=1,end=1) ! END e-nucleus SET x
 enddo loop_sets_i ! sets_i

 return

1 call errstop('READ_MOL_DENSITY','Problem reading MOLECULAR DENSITY &
   &section of expval.data file.')

 END SUBROUTINE read_mol_density


 SUBROUTINE read_twop_dm_mom
 IMPLICIT NONE
 INTEGER i,set
 REAL (dp) dumr
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_TWOP_DM_MOM','Problem reading expval.data file.')
  if(ierr<0)call errstop('READ_TWOP_DM_MOM','Could not find START TWO-PARTICLE &
   &DENSITY MATRIX (MOMENTUM SPACE) in expval.data file.')
  if(trim(adjustl(char_80))=='START TWO-PARTICLE DENSITY MATRIX &
   &(MOMENTUM SPACE)')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_twop_dm_mom=trim(adjustl(char_80))
 if(accum_method_twop_dm_mom/='VMC'.and.&
  &accum_method_twop_dm_mom/='DMC')&
  &call errstop('READ_TWOP_DM_MOM','Momentum density accumulation method &
  &(VMC or DMC) in expval.data not recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)twop_dm_mom_nsets
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)dumr
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)dumr

 allocate(expval_twop_dm_mom_weight(twop_dm_mom_nsets),&
  &expval_twop_dm_mom(expval_ngvec(1),twop_dm_mom_nsets),&
  &expval_twop_dm_mom2(expval_ngvec(1),twop_dm_mom_nsets),&
  &expval_twop_dm_mom_k(expval_ngvec(1),twop_dm_mom_nsets),&
  &twop_dm_mom_nptypes_in_set(twop_dm_mom_nsets),&
  &twop_dm_mom_ptype_in_set(2,((num_particle_types+1)*num_particle_types)/2,&
  &twop_dm_mom_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_TWOP_DM_MOM','Allocation error: &
  &expval_twop_dm_mom arrays')
 expval_twop_dm_mom_weight=0.d0 ; expval_twop_dm_mom=0.d0
 expval_twop_dm_mom2=0.d0 ; expval_twop_dm_mom_k=0.d0
 twop_dm_mom_nptypes_in_set=0 ; twop_dm_mom_ptype_in_set=0

 do set=1,twop_dm_mom_nsets

  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START SET '//trim(i2s(set)))then
   call errstop('READ_TWOP_DM_MOM','Error reading START SET x line for x = '//&
    &trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)twop_dm_mom_nptypes_in_set(set)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)&
   &twop_dm_mom_ptype_in_set(1:2,1:twop_dm_mom_nptypes_in_set(set),set)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_twop_dm_mom_weight(set)
  read(io_expval,*,err=1,end=1)

  do i=1,expval_ngvec(1)
   read(io_expval,*,err=1,end=1)expval_twop_dm_mom_k(i,set),&
    &expval_twop_dm_mom2(i,set),expval_twop_dm_mom(i,set)
  enddo

  read(io_expval,*,err=1,end=1) ! END SET x

 enddo ! sets

 return

1 call errstop('READ_TWOP_DM_MOM','Problem reading TWO-PARTICLE DENSITY &
  &MATRIX (MOMENTUM SPACE) section of expval.data file.')

 END SUBROUTINE read_twop_dm_mom


 SUBROUTINE read_cond_fraction_mom
 IMPLICIT NONE
 INTEGER i,set
 REAL (dp) dumr
 CHARACTER(80) char_80

 do
  read(io_expval,'(a)',iostat=ierr)char_80
  if(ierr>0)call errstop('READ_COND_FRAC_MOM','Problem reading expval.data &
   &file.')
  if(ierr<0)call errstop('READ_COND_FRAC_MOM','Could not find START &
   &CONDENSATE FRACTION (MOMENTUM SPACE) in expval.data file.')
  if(trim(adjustl(char_80))=='START CONDENSATE FRACTION (MOMENTUM SPACE)')exit
 enddo

 read(io_expval,*,err=1,end=1)
 read(io_expval,'(a)',err=1,end=1)char_80
 accum_method_cond_frac_mom=trim(adjustl(char_80))
 if(accum_method_cond_frac_mom/='VMC'.and.&
  &accum_method_cond_frac_mom/='DMC')call errstop('READ_COND_FRAC_MOM',&
  &'Momentum density accumulation method (VMC or DMC) in expval.data not &
  &recognized.')

 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)cond_frac_mom_nsets
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)dumr
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)dumr

 allocate(expval_cond_frac_mom_weight(cond_frac_mom_nsets),&
  &expval_cond_frac_mom(expval_ngvec(1),cond_frac_mom_nsets),&
  &expval_cond_frac_mom2(expval_ngvec(1),cond_frac_mom_nsets),&
  &expval_cond_frac_mom_k(expval_ngvec(1),cond_frac_mom_nsets),&
  &cond_frac_mom_nptypes_in_set(cond_frac_mom_nsets),&
  &cond_frac_mom_ptype_in_set(2,((num_particle_types+1)*num_particle_types)/2,&
  &cond_frac_mom_nsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_COND_FRAC_MOM','Allocation error: &
  &expval_cond_frac_mom arrays')
 expval_cond_frac_mom_weight=0.d0 ; expval_cond_frac_mom=0.d0
 expval_cond_frac_mom2=0.d0 ; expval_cond_frac_mom_k=0.d0
 cond_frac_mom_nptypes_in_set=0 ; cond_frac_mom_ptype_in_set=0

 do set=1,cond_frac_mom_nsets

  read(io_expval,'(a)',err=1,end=1)char_80
  if(trim(adjustl(char_80))/='START SET '//trim(i2s(set)))then
   call errstop('READ_COND_FRAC_MOM','Error reading START SET x line for x = '&
    &//trim(i2s(set)))
  endif

  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)cond_frac_mom_nptypes_in_set(set)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)&
   &cond_frac_mom_ptype_in_set(1:2,1:cond_frac_mom_nptypes_in_set(set),set)
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)expval_cond_frac_mom_weight(set)
  read(io_expval,*,err=1,end=1)

  do i=1,expval_ngvec(1)
   read(io_expval,*,err=1,end=1)expval_cond_frac_mom_k(i,set),&
    &expval_cond_frac_mom2(i,set),expval_cond_frac_mom(i,set)
  enddo

  read(io_expval,*,err=1,end=1) ! END SET x

 enddo ! sets

 return

1 call errstop('READ_COND_FRAC_MOM','Problem reading CONDENSATE FRACTION &
  &(MOMENTUM SPACE) section of expval.data file.')

 END SUBROUTINE read_cond_fraction_mom


 SUBROUTINE talk_to_user
!------------------------------------------------------!
! Discuss with the user precisely what he wants to do. !
!------------------------------------------------------!
 IMPLICIT NONE
 CHARACTER(2048) char_2048
 INTEGER i,j,k,k2,blocksize,reblock(100),kstart,kend,jset
 REAL(dp) temp

! Ask what to plot.
 if(nex>1)then
  i=0
  do while(i<1.or.i>nex)
   write(6,*)
   write(6,*)'Which of these do you want to plot?'
   write(6,*)
   read(5,*,iostat=ierr)i
   if(ierr>0)i=0
   write(6,*)
  enddo
  iplot=inex(i)
 else
  iplot=inex(1)
 endif

! Ask additional expval-specific questions.
 select case(iplot)

 case(1) ! Density
  if(den_nsets>1)then
   iset=0
   do while(iset<1.or.iset>den_nsets)
    write(6,*)'The DENSITY block contains ',trim(i2s(den_nsets)),' sets'
    write(6,*)'corresponding to different particle types.'
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif

 case(2) ! Spin density
  if(sden_nsets>1)then
   iset=0
   do while(iset<1.or.iset>sden_nsets)
    write(6,*)'The SPIN DENSITY block contains ',trim(i2s(sden_nsets)),' sets'
    write(6,*)'corresponding to different spins and/or particle types.'
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif

 case(3) ! Spin density matrix
  if(sdenm_nsets>1)then
   if(count(nele_expval>0)/=1)then
    call errstop('TALK_TO_USER','Only one spin-density matrix set should have &
     &a nonzero number of particles at present.')
    iset=0
    do while(iset<1.or.iset>sdenm_nsets)
     write(6,*)'The SPIN-DENSITY MATRIX block contains ',&
      &trim(i2s(sdenm_nsets)),' sets'
     write(6,*)'corresponding to different spins and/or particle types.'
     write(6,*)'Which do you want to plot?'
     write(6,*)
     read(5,*,iostat=ierr)iset
     if(ierr>0)iset=0
     write(6,*)
    enddo
   else
    iset=-1
    do jset=1,sdenm_nsets
     if(nele_expval(jset)>0)then
      iset=jset
      exit
     endif
    enddo ! jset
    if(iset==-1)call errstop('TALK_TO_USER','Bug.')
   endif ! Just one set has a nonzero number of particles.
  else
   iset=1
  endif

 case(4) ! Reciprocal space PCF
  if(pcf_nsets>1)then
   iset=-1
   do while(iset<0.or.iset>pcf_nsets)
    write(6,*)'The RECIPROCAL-SPACE PCF block contains ',trim(i2s(pcf_nsets)),&
     &' sets'
    write(6,*)'corresponding to the following particle pairs'
    write(6,*)
    do j=1,pcf_nsets
     write(6,*)'(',trim(i2s(j)),') ',trim(i2s(pcf_type_fixed)),'-',&
      &trim(i2s(pcf_ptype(j)))
    enddo
    write(6,*)
    write(6,*)'Which do you want to plot? (0 for all)'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=-1
    write(6,*)
   enddo
  else
   iset=1
  endif

 case(5) ! Spherical PCF
  if(pcfs_nsets>1)then
   iset=-1
   do while(iset<0.or.iset>pcfs_nsets)
    write(6,*)'The SPHERICAL PCF block contains ',trim(i2s(pcfs_nsets)),' sets'
    write(6,*)'corresponding to the following particle pairs :'
    do j=1,pcfs_nsets
     if(nele_expval(pcfs_ptype(1,j))==0.or.nele_expval(pcfs_ptype(2,j))==0&
      &.or.(pcfs_ptype(1,j)==pcfs_ptype(2,j)&
      &.and.nele_expval(pcfs_ptype(1,j))==1))then
      write(6,*)'(',trim(i2s(j)),') ',trim(i2s(pcfs_ptype(1,j))),'-',&
       &trim(i2s(pcfs_ptype(2,j)))//'   <- NO DATA'
     else
      write(6,*)'(',trim(i2s(j)),') ',trim(i2s(pcfs_ptype(1,j))),'-',&
       &trim(i2s(pcfs_ptype(2,j)))
     endif
    enddo ! j
    write(6,*)
    write(6,*)'Which do you want to plot? (0 for all)'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=-1
    write(6,*)
   enddo

   write(6,*)
   write(6,*)'These data were accumulated in ',trim(i2s(pcfs_nbins)),' bins,&
    & and can be smoothed by reblocking.'
   j=1
   do i=1,100
    if(mod(pcfs_nbins,i)==0)then
     reblock(j)=i
     j=j+1
    endif
   enddo
   write(6,*)'Suggested block sizes are :'
   call write_list_int(j-1,reblock,16,4,1,6)
   write(6,*)'or any other divisor of the total number of bins.'
   write(6,*)'Input a block size (0 for no reblocking) :'
   blocksize=-1
   do while(blocksize<0)
    read(5,*,iostat=ierr)blocksize
    if(ierr/=0)blocksize=-1
    write(6,*)
   enddo
   if(blocksize>0)then
    if(mod(pcfs_nbins,blocksize)/=0)then
     write(6,*)'Warning: block size not a divisor of number of bins.'
     write(6,*)'         Last '//trim(i2s(mod(pcfs_nbins,blocksize)))//&
      &'bins will be ignored.'
    endif
    write(6,*)
    pcfs_nbins_new=pcfs_nbins/blocksize
    allocate(expval_pcfs_reblocked(pcfs_nbins_new,num_particle_types,&
     &num_particle_types),stat=ialloc)
    if(ialloc/=0)call errstop('TALK_TO_USER','Allocation error : &
     &expval_pcfs_reblocked.')
    do i=1,num_particle_types
     do j=1,i
      do k=1,pcfs_nbins_new
       kstart=(k-1)*blocksize+1
       kend=kstart+blocksize-1
       temp=0.d0
       do k2=kstart,kend
        temp=temp+expval_pcfs(k2,j,i)
       enddo
       expval_pcfs_reblocked(k,j,i)=temp
      enddo
     enddo
    enddo
    deallocate(expval_pcfs)
    reblocked=.true.
   else
    reblocked=.false.
    write(6,*)'No reblocking transformation will be performed.'
    write(6,*)
   endif
  else ! pcfs_nsets=1
   iset=1
  endif
  plot_dim=1

 case(6) ! Localization tensor
  if(lt_nsets>1)then
   iset=0
   do while(iset<1.or.iset>lt_nsets)
    write(6,*)'The LOCALIZATION TENSOR block contains ',trim(i2s(lt_nsets)),&
     &' sets.'
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif

 case(7) ! Structure factor

  if(sf_nsets>1)then
   iset=-1
   do while(iset<0.or.iset>sf_nsets+1)
    write(6,*)'The STRUCTURE FACTOR block contains ',trim(i2s(sf_nsets)),&
     &' sets corresponding to'
    write(6,*)'different particle pairs.'
    write(6,*)
    write(6,*)'You can plot :'
    write(6,*)'(0) All spin pairs.'
    do j=1,sf_nsets
     write(6,*)'(',trim(i2s(j)),') ',trim(i2s(sf_ptype(1,j))),'-',&
      &trim(i2s(sf_ptype(2,j)))
    enddo
    write(6,*)'(',trim(i2s(sf_nsets+1)),') The total summed over spin pairs &
     &(usual).'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=-1
    write(6,*)
   enddo
  else
   iset=1
  endif

  write(6,*)'These data were accumulated on a grid of points {G} in reciprocal &
   &space.'
  sf_option=0
  do while(sf_option<1.or.sf_option>2)
   write(6,*)'You can:'
   write(6,*)
   write(6,*)'(1) Compute the spherical average (over stars of G vectors).'
   write(6,*)'(2) Dump the whole grid into a file (NOT YET FUNCTIONAL).'
   write(6,*)
   read(5,*,iostat=ierr)sf_option
   write(6,*)
  enddo
  if(sf_option==1)plot_dim=1
  if(sf_option==2)then
   call errstop('TALK_TO_USER','This option is not yet implemented.')
  endif

 case(8) ! Spherical structure factor

  if(sf_sph_nsets>1)then
   iset=-1
   do while(iset<0.or.iset>sf_sph_nsets+1)
    call wordwrap('The SPHERICAL STRUCTURE FACTOR block contains '&
     &//trim(i2s(sf_sph_nsets))//' sets.')
    write(6,*)'You can plot :'
    write(6,*)'(0) All spin pairs.'
    do j=1,sf_sph_nsets
     write(6,*)'(',trim(i2s(j)),') ',trim(i2s(sf_sph_ptype(1,j))),'-',&
      &trim(i2s(sf_sph_ptype(2,j)))
    enddo
    write(6,*)'(',trim(i2s(sf_sph_nsets+1)),') The total summed over spin &
     &pairs (usual).'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)i=-1
    write(6,*)
   enddo
  else
   iset=1
  endif
  plot_dim=1

 case(9) ! One-particle density matrix
  if(onep_dm_nsets>1)then
   iset=0
   do while(iset<1.or.iset>onep_dm_nsets)
    write(6,*)'The ONE-PARTICLE DENSITY MATRIX block contains ',&
     &trim(i2s(onep_dm_nsets)),' sets corresponding to different &
     &particle types:'
    do i=1,onep_dm_nsets
     char_2048=trim(i2s(onep_dm_ptype_in_set(1,i)))
     do j=2,onep_dm_nptypes_in_set(i)
      char_2048=trim(char_2048)//', '//trim(i2s(onep_dm_ptype_in_set(j,i)))
     enddo ! j
     write(6,*)trim(i2s(i))//': '//trim(char_2048)
    enddo ! i
    write(6,*)
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif
  plot_dim=1

 case(10) ! Two-particle density matrix
  if(twop_dm_nsets>1)then
   iset=0
   do while(iset<1.or.iset>twop_dm_nsets)
    write(6,*)'The TWO-PARTICLE DENSITY MATRIX block contains ',&
     &trim(i2s(twop_dm_nsets)),' sets corresponding to different &
     &particle types:'
    do i=1,twop_dm_nsets
     char_2048=trim(i2s(twop_dm_ptype_in_set(1,1,i)))//'-'//&
      &trim(i2s(twop_dm_ptype_in_set(2,1,i)))
     do j=2,twop_dm_nptypes_in_set(i)
      char_2048=trim(char_2048)//', '//trim(i2s(twop_dm_ptype_in_set(1,j,i)))&
       &//'-'//trim(i2s(twop_dm_ptype_in_set(2,j,i)))
     enddo ! j
     write(6,*)trim(i2s(i))//': '//trim(char_2048)
    enddo ! i
    write(6,*)
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif
  plot_dim=1 ; twop_dm_with_onep_dm=.false. ; twop_dm_rescale=.false.
  if(onep_density_mat_in_expval)then
   i=-1
   do while(i<0.or.i>2)
    write(6,*)'Perform post-processing?'
    write(6,*)'0: no post-processing'
    write(6,*)'1: renormalize so lim r->inf is condensate fraction'
    write(6,*)'2: cancel one-body contributions and renormalize like option 1'
    write(6,*)
    read(5,*,iostat=ierr)i
    if(ierr>0)i=0
    write(6,*)
   enddo
   twop_dm_rescale=i>0
   twop_dm_with_onep_dm=i==2
   if(twop_dm_with_onep_dm)then
    do j=1,twop_dm_nptypes_in_set(iset)
     if(twop_dm_ptype_in_set(1,j,iset)==twop_dm_ptype_in_set(2,j,iset))then
      write(6,*)
      write(6,*)'WARNING: This set contains identical-particle pairs.  &
       &Subtracting one-body'
      write(6,*)'         contributions for identical-particle pairs requires &
       &information which'
      write(6,*)'         is not currently accumulated in expval.data.  The &
       &resulting plot is'
      write(6,*)'         likely to be incorrect.'
      write(6,*)'         You should either plot the two-particle density &
       &matrix without'
      write(6,*)'         subtracting one-body contributions or plot the &
       &"condensate fraction'
      write(6,*)'         estimator" expectation value which does not have &
       &this problem.'
      write(6,*)
     endif
    enddo ! j
   endif
  endif

 case(11) ! Condensate fraction
  if(cond_frac_nsets>1)then
   iset=0
   do while(iset<1.or.iset>cond_frac_nsets)
    write(6,*)'The CONDENSATE FRACTION block contains ',&
     &trim(i2s(cond_frac_nsets)),' sets corresponding to different &
     &particle types:'
    do i=1,cond_frac_nsets
     char_2048=trim(i2s(cond_frac_ptype_in_set(1,1,i)))//'-'//&
      &trim(i2s(cond_frac_ptype_in_set(2,1,i)))
     do j=2,cond_frac_nptypes_in_set(i)
      char_2048=trim(char_2048)//', '//&
       &trim(i2s(cond_frac_ptype_in_set(1,j,i)))//'-'//&
       &trim(i2s(cond_frac_ptype_in_set(2,j,i)))
     enddo ! j
     write(6,*)trim(i2s(i))//': '//trim(char_2048)
    enddo ! i
    write(6,*)
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif
  plot_dim=1

 case(12) ! Momentum density
  if(mom_den_nsets>1)then
   iset=0
   do while(iset<1.or.iset>mom_den_nsets)
    write(6,*)'The MOMENTUM DENSITY block contains ',&
     &trim(i2s(mom_den_nsets)),' sets corresponding to different &
     &particle types.'
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif
  plot_dim=1

 case(13) ! Finite density
  if((fin_den_nsets_ij+fin_den_nsets_i)>1)then
   iset=0
   do while(iset<1.or.iset>2)
    call wordwrap('The FINITE DENSITY block contains (1) electron-electron &
     &moments "ij" and (2) electron-nucleus moments "i". Radial moments &
     &<r^n> for n=-2,-1,1,2,3 will be calculated.')
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif
  plot_dim=1
  write(6,*)'Enter Z, the charge on the nucleus.'
  write(6,*)
  read(5,*,iostat=ierr)fin_den_z
  if(ierr>0)write(6,*)'Value of Z not recognized.'
  write(6,*)
  write(6,*)'Enter ionization potential of atom/ion. This is only neccessary &
   &if fitting with constraints.'
  write(6,*)
  read(5,*,iostat=ierr)fin_den_ip
  if(ierr>0)write(6,*)'Value of ionization potential not recognized.'
  write(6,*)

 case(14) ! Molecular density
  if(mol_den_nsets>1)then
   iset=0
   do while(iset<1.or.iset>mol_den_nsets)
    write(6,*)'The MOLECULAR DENSITY block contains ',&
     &trim(i2s(mol_den_nsets)),' sets'
    write(6,*)'corresponding to different particle types.'
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif
  select case(count(mol_den_norder==1))
  case(0)
   plot_dim=3
  case(1)
   plot_dim=2
  case(2,3)
   plot_dim=1
  end select

 case(15) ! Two-particle density matrix (momentum space)
  if(twop_dm_mom_nsets>1)then
   iset=0
   do while(iset<1.or.iset>twop_dm_mom_nsets)
    write(6,*)'The TWO-PARTICLE DENSITY MATRIX block contains ',&
     &trim(i2s(twop_dm_mom_nsets)),' sets corresponding to different &
     &particle types:'
    do i=1,twop_dm_mom_nsets
     char_2048=trim(i2s(twop_dm_mom_ptype_in_set(1,1,i)))//'-'//&
      &trim(i2s(twop_dm_mom_ptype_in_set(2,1,i)))
     do j=2,twop_dm_mom_nptypes_in_set(i)
      char_2048=trim(char_2048)//', '//&
       &trim(i2s(twop_dm_mom_ptype_in_set(1,j,i)))&
       &//'-'//trim(i2s(twop_dm_mom_ptype_in_set(2,j,i)))
     enddo ! j
     write(6,*)trim(i2s(i))//': '//trim(char_2048)
    enddo ! i
    write(6,*)
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif
  plot_dim=1 ; ; twop_dm_mom_rescale=.false.
  i=-1
  do while(i<0.or.i>1)
   write(6,*)'Perform post-processing?'
   write(6,*)'0: no post-processing'
   write(6,*)'1: renormalize so lim r->inf is condensate fraction'
   write(6,*)
   read(5,*,iostat=ierr)i
   if(ierr>0)i=0
   write(6,*)
  enddo
  twop_dm_mom_rescale=i>0

 case(16) ! Condensate fraction (momentum space)
  if(cond_frac_mom_nsets>1)then
   iset=0
   do while(iset<1.or.iset>cond_frac_mom_nsets)
    write(6,*)'The CONDENSATE FRACTION block contains ',&
     &trim(i2s(cond_frac_mom_nsets)),' sets corresponding to different &
     &particle types:'
    do i=1,cond_frac_mom_nsets
     char_2048=trim(i2s(cond_frac_mom_ptype_in_set(1,1,i)))//'-'//&
      &trim(i2s(cond_frac_mom_ptype_in_set(2,1,i)))
     do j=2,cond_frac_mom_nptypes_in_set(i)
      char_2048=trim(char_2048)//', '//&
       &trim(i2s(cond_frac_mom_ptype_in_set(1,j,i)))//'-'//&
       &trim(i2s(cond_frac_mom_ptype_in_set(2,j,i)))
     enddo ! j
     write(6,*)trim(i2s(i))//': '//trim(char_2048)
    enddo ! i
    write(6,*)
    write(6,*)'Which do you want to plot?'
    write(6,*)
    read(5,*,iostat=ierr)iset
    if(ierr>0)iset=0
    write(6,*)
   enddo
  else
   iset=1
  endif
  plot_dim=1

 case default
  call errstop('TALK_TO_USER','Invalid code for expval plot selection - bug.')

 end select

! Do required setup based on user responses.
 select case(iplot)
  case(5,7,8,9,10,11,12,13,14,15,16)
   continue ! No setup
  case default
   call setup_fourier_basis
 end select

 END SUBROUTINE talk_to_user


 SUBROUTINE write_expval
 IMPLICIT NONE
 select case(plot_dim)
 case(1)
  open(10,file='lineplot.dat',status='unknown',action='write',iostat=ierr)
  write(6,*)'Plotting expval to lineplot.dat. View with xmgr/grace.'
 case(2)
  open(10,file='2Dplot.dat',status='unknown',action='write',iostat=ierr)
  write(6,*)'Plotting expval to 2Dplot.dat. View with gnuplot (via plot_2D &
   &script)'
 case(3)
  open(10,file='3Dplot.dat',status='unknown',action='write',iostat=ierr)
  write(6,*)'Plotting expval to 3Dplot.dat. View with gnuplot (via plot_2D &
   &script)'
 case default
  call errstop('WRITE_EXPVAL','Invalid plot_dim code - bug.')
 end select
 write(6,*)
 if(ierr/=0)call errstop('WRITE_EXPVAL','Error opening plot file for writing.')
 select case(iplot)
 case(1) ; call plot_density
 case(2) ; call plot_spin_density
 case(3) ; call plot_spin_density_matrix
 case(4) ; call plot_pair_corr
 case(5) ; call plot_pair_corr_sph
 case(6) ; call plot_localization_tensor
 case(7) ; call plot_structure_factor
 case(8) ; call plot_structure_factor_sph
 case(9) ; call plot_onep_density_matrix
 case(10); call plot_twop_density_matrix
 case(11); call plot_cond_fraction
 case(12); call plot_mom_den
 case(13); call plot_finite_density
 case(14); call plot_mol_density
 case(15); call plot_twop_dm_mom
 case(16); call plot_cond_fraction_mom
 case default
  call errstop('WRITE_EXPVAL','Invalid code for expval plot selection - bug.')
 end select

 close(10)

 END SUBROUTINE write_expval


 SUBROUTINE plot_density
!------------------------------------------!
! Plot the DENSITY block from expval.data. !
!------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,g
 REAL(dp) x

 select case (plot_dim)
 case(1) ! line

  if(r_ab(2)==0.d0.and.r_ab(3)==0.d0)then
   x=apoint(1)
  elseif(r_ab(1)==0.d0.and.r_ab(3)==0.d0)then
   x=apoint(2)
  elseif(r_ab(1)==0.d0.and.r_ab(2)==0.d0)then
   x=apoint(3)
  else
   x=0.d0
  endif
  r(:)=apoint(:)
  do i=1,npoint1
   expval=0.d0
   call compute_fourier_basis(r)
   do g=1,expval_ngvec(den_gset)
    expval=expval+real(expval_den(g,iset)*expval_expigdotr(g),dp)
   enddo
   write(10,'(e18.10,1x,e18.10)')x,expval
   r(:)=r(:)+step_ab(:)
   x=x+step_ab_length
  enddo

 case(2) ! plane

  r(:)=apoint(:) ; r1(:)=apoint(:)
  do j=1,npoint2
   do i=1,npoint1
    expval=0.d0
    call compute_fourier_basis(r)
    do g=1,expval_ngvec(den_gset)
     expval=expval+real(expval_den(g,iset)*expval_expigdotr(g),dp)
    enddo
    write(10,'(e18.10,3(1x,e18.10))')r(:),expval
    r(:)=r(:)+step_ab(:)
   enddo
   write(10,*)
   r1(:)=r1(:)+step_ac(:) ; r(:)=r1(:)
  enddo

 case(3) ! volume

  do k=0,npoint3-1
   do j=0,npoint2-1
    do i=1,npoint1
     expval=0.d0
     r(:)=apoint(:)+dble(i-1)*step_ab(:)+dble(j)*step_ac(:)+dble(k)*step_ad(:)
     call compute_fourier_basis(r)
     do g=1,expval_ngvec(den_gset)
      expval=expval+real(expval_den(g,iset)*expval_expigdotr(g),dp)
     enddo
     write(10,'(e18.10,3(1x,e18.10))')r(:),expval
    enddo
    write(10,*)
   enddo
  enddo

 case default
  call errstop('PLOT_DENSITY','Invalid dimensionality for plot.')

 end select

 END SUBROUTINE plot_density


 SUBROUTINE plot_spin_density
!-----------------------------------------------!
! Plot the SPIN DENSITY block from expval.data. !
!-----------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,g
 REAL(dp) x

 select case (plot_dim)
 case(1) ! line

  if(r_ab(2)==0.d0.and.r_ab(3)==0.d0)then
   x=apoint(1)
  elseif(r_ab(1)==0.d0.and.r_ab(3)==0.d0)then
   x=apoint(2)
  elseif(r_ab(1)==0.d0.and.r_ab(2)==0.d0)then
   x=apoint(3)
  else
   x=0.d0
  endif
  r(:)=apoint(:)
  do i=1,npoint1
   expval=0.d0
   call compute_fourier_basis(r)
   do g=1,expval_ngvec(sden_gset)
    expval=expval+real(expval_sden(g,iset)*expval_expigdotr(g),dp)
   enddo
   write(10,'(e18.10,1x,e18.10)')x,expval
   r(:)=r(:)+step_ab(:)
   x=x+step_ab_length
  enddo

 case(2) ! plane

  r(:)=apoint(:) ; r1(:)=apoint(:)
  do j=1,npoint2
   do i=1,npoint1
    expval=0.d0
    call compute_fourier_basis(r)
    do g=1,expval_ngvec(sden_gset)
     expval=expval+real(expval_sden(g,iset)*expval_expigdotr(g),dp)
    enddo
    write(10,'(e18.10,3(1x,e18.10))')r(:),expval
    r(:)=r(:)+step_ab(:)
   enddo
   write(10,*)
   r1(:)=r1(:)+step_ac(:) ; r(:)=r1(:)
  enddo

 case(3) ! volume

  do k=0,npoint3-1
   do j=0,npoint2-1
    do i=1,npoint1
     expval=0.d0
     r(:)=apoint(:)+dble(i-1)*step_ab(:)+dble(j)*step_ac(:)+dble(k)*step_ad(:)
     call compute_fourier_basis(r)
     do g=1,expval_ngvec(sden_gset)
      expval=expval+real(expval_sden(g,iset)*expval_expigdotr(g),dp)
     enddo
     write(10,'(e18.10,3(1x,e18.10))')r(:),expval
    enddo
    write(10,*)
   enddo
  enddo

 case default
  call errstop('PLOT_SPIN_DENSITY','Invalid dimensionality for plot.')
 end select

 END SUBROUTINE plot_spin_density


 SUBROUTINE plot_spin_density_matrix
!------------------------------------------------!
! Plot the SPIN-DENSITY MATRIX from expval.data. !
!------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,g,sdm_el,ierr,sdm_or_mag
 REAL(dp) x,xarr(npoint1),mag(npoint1,3),den(npoint1)
 REAL(dp),PARAMETER :: tol=1.d-3
 COMPLEX(dp) sdm(npoint1,4),rec_vol

! SDM data from CASINO are normalised such that SDM_11+SDM_22=N, where N is the
! electron number.  Divide by the cell volume to get the usual definition of the
! SDM, etc.
 rec_vol=cmplx(1.d0/simcell_volume,0.d0,dp)

! Does the user want SDM itself or the magnetisation density?
 do
   write(6,*)'Would you like to plot (1) the spin-density matrix itself or &
     &(2) the'
   write(6,*)'magnetization density?'
   read(5,*,iostat=ierr)sdm_or_mag
   if(ierr/=0)sdm_or_mag=-1
   if(sdm_or_mag==1.or.sdm_or_mag==2)exit
   write(6,*)'Please try again.  Enter "1" or "2".'
 enddo
 write(6,*)

 select case (plot_dim)
 case(1) ! line

! Evaluate the set of x points.
  if(r_ab(2)==0.d0.and.r_ab(3)==0.d0)then
   x=apoint(1)
  elseif(r_ab(1)==0.d0.and.r_ab(3)==0.d0)then
   x=apoint(2)
  elseif(r_ab(1)==0.d0.and.r_ab(2)==0.d0)then
   x=apoint(3)
  else
   x=0.d0
  endif
  do i=1,npoint1
   xarr(i)=x
   x=x+step_ab_length
  enddo ! i

! Evaluate the spin density matrix at each point along the line.
  do sdm_el=1,4
   r(:)=apoint(:)
   do i=1,npoint1
    call compute_fourier_basis(r)
    sdm(i,sdm_el)=cmplx(0.d0,0.d0,dp)
    do g=1,expval_ngvec(sdenm_gset)
     sdm(i,sdm_el)=sdm(i,sdm_el)+expval_sdenm(g,sdm_el,iset)*expval_expigdotr(g)
    enddo ! g
    sdm(i,sdm_el)=sdm(i,sdm_el)*rec_vol
    r(:)=r(:)+step_ab(:)
   enddo ! i
  enddo ! sdm_el

! Check the SDM has the required properties.
! This keeps failing, so just issue a warning.
  do i=1,npoint1
   if(abs(aimag(sdm(i,1)))>tol*abs(dble(sdm(i,1))))then
    call errwarn('PLOT_SPIN_DENSITY_MATRIX','The (1,1) element of the &
     &spin-density matrix is complex.  It should be real.')
    write(6,*)'(1,1) element of SDM at point '//trim(i2s(i))//': ',sdm(i,1)
    write(6,*)
    exit
   endif ! SDM_11 not real.
   if(abs(dble(sdm(i,2))-dble(sdm(i,3)))>tol*abs(sdm(i,2)))then
    call errwarn('PLOT_SPIN_DENSITY_MATRIX','The (1,2) and (2,1) elements of &
     &the spin-density matrix should be a complex-conjugate pair.')
    write(6,*)'(1,2) element of SDM at point '//trim(i2s(i))//': ',sdm(i,2)
    write(6,*)'(2,1) element of SDM at point '//trim(i2s(i))//': ',sdm(i,3)
    write(6,*)
    exit
   endif ! SDM_12 and SDM_22 not a complex conjugate pair.
   if(abs(aimag(sdm(i,2))+aimag(sdm(i,3)))>tol*abs(sdm(i,2)))then
    call errwarn('PLOT_SPIN_DENSITY_MATRIX','The (1,2) and (2,1) elements of &
     &the spin-density matrix should be a complex-conjugate pair.')
    write(6,*)'(1,2) element of SDM at point '//trim(i2s(i))//': ',sdm(i,2)
    write(6,*)'(2,1) element of SDM at point '//trim(i2s(i))//': ',sdm(i,3)
    write(6,*)
    exit
   endif ! SDM_12 and SDM_22 not a complex conjugate pair.
   if(abs(aimag(sdm(i,4)))>tol*abs(dble(sdm(i,4))))then
    call errwarn('PLOT_SPIN_DENSITY_MATRIX','The (2,2) element of the &
     &spin-density matrix is complex.  It should be real.')
    write(6,*)'(2,2) element of SDM at point '//trim(i2s(i))//': ',sdm(i,4)
    write(6,*)
    exit
   endif ! SDM_22 not real.
  enddo ! i
! Average the off-diagonal elements of the SDM (which are supposed to be a
! complex-conjugate pair).
  do i=1,npoint1
   sdm(i,2)=0.5d0*(sdm(i,2)+conjg(sdm(i,3)))
   sdm(i,3)=conjg(sdm(i,2))
  enddo ! i

! Write out the SDM.
  if(sdm_or_mag==1)then
    do i=1,npoint1
      write(10,'(e18.10,1x,e18.10)')xarr(i),dble(sdm(i,1))
    enddo ! i
    write(10,'(a)')'&'
    do i=1,npoint1
      write(10,'(e18.10,1x,e18.10)')xarr(i),dble(sdm(i,2))
    enddo ! i
    write(10,'(a)')'&'
    do i=1,npoint1
      write(10,'(e18.10,1x,e18.10)')xarr(i),aimag(sdm(i,2))
    enddo ! i
    write(10,'(a)')'&'
    do i=1,npoint1
      write(10,'(e18.10,1x,e18.10)')xarr(i),dble(sdm(i,4))
    enddo ! i
    write(6,*)'Let (SDM) denote spin-density maxtrix and let M be &
      &magnetization density.'
    write(6,*)'The 1st data set in lineplot.dat is (SDM)_11.'
    write(6,*)'The 2nd data set in lineplot.dat is Re[(SDM)_12] = M_x/2.'
    write(6,*)'The 3rd data set in lineplot.dat is Im[(SDM)_12] = -M_y/2.'
    write(6,*)'The 4th data set in lineplot.dat is (SDM)_22.'
    write(6,*)
    write(6,*)'The density is n=(SDM)_11+(SDM)_22.'
    write(6,*)'M_z=[(SDM)_11-(SDM)_22].'
    write(6,*)
  else
    ! Calculate the magnetisation density.
    do i=1,npoint1
      mag(i,1)=2.d0*dble(sdm(i,2))             ! M_x
      mag(i,2)=-2.d0*aimag(sdm(i,2))           ! M_y
      mag(i,3)=dble(sdm(i,1))-dble(sdm(i,4))   ! M_z
      den(i)=dble(sdm(i,1))+dble(sdm(i,4))     ! n
    enddo ! i
    do j=1,3
      do i=1,npoint1
        write(10,'(e18.10,1x,e18.10)')xarr(i),mag(i,j)
      enddo ! i
      write(10,'(a)')'&'
    enddo ! j
    do i=1,npoint1
      write(10,'(e18.10,1x,e18.10)')xarr(i),den(i)
    enddo ! i
    write(6,*)'Let M be the magnetization density and n be the density.'
    write(6,*)'The 1st data set in lineplot.dat is M_x.'
    write(6,*)'The 2nd data set in lineplot.dat is M_y.'
    write(6,*)'The 3rd data set in lineplot.dat is M_z.'
    write(6,*)'The 4th data set in lineplot.dat is n.'
    write(6,*)
  endif ! sdm_or_mag

 case(2) ! plane
  call errstop('PLOT_SPIN_DENSITY_MATRIX','The next person who needs a 2D &
   &plot of the spin-density matrix is very welcome to code it up.')
 case(3) ! volume
  call errstop('PLOT_SPIN_DENSITY_MATRIX','The next person who needs a 3D &
   &plot of the spin-density matrix is very welcome to code it up.')
 case default
  call errstop('PLOT_SPIN_DENSITY_MATRIX','Invalid dimensionality for plot.')
 end select

 END SUBROUTINE plot_spin_density_matrix


 SUBROUTINE plot_pair_corr
!--------------------------------------------------------!
! Plot the RECIPROCAL-SPACE PCF block from expval.data.  !
!--------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,g,jset,set
 REAL(dp) x

 if(iset/=0)then
  pcf_nsets=1 ; jset=iset
 endif

 select case (plot_dim)
 case(1) ! line

  do set=1,pcf_nsets
   if(iset==0)jset=set

   if(r_ab(2)==0.d0.and.r_ab(3)==0.d0)then
    x=apoint(1)
   elseif(r_ab(1)==0.d0.and.r_ab(3)==0.d0)then
    x=apoint(2)
   elseif(r_ab(1)==0.d0.and.r_ab(2)==0.d0)then
    x=apoint(3)
   else
    x=0.d0
   endif
   r(:)=apoint(:)
   do i=1,npoint1
    expval=0.d0 ; expval2=0.d0
    call compute_fourier_basis(r)
    do g=1,expval_ngvec(pcf_gset)
     expval=expval+real(expval_pcf(g,jset)*expval_expigdotr(g),dp)
     expval2=expval2+real(expval_sden(g,jset)*expval_expigdotr(g),dp)
    enddo
    write(10,'(e18.10,1x,e18.10)')x,expval/expval2
    r(:)=r(:)+step_ab(:)
    x=x+step_ab_length
   enddo

   if(iset==0.and.set<pcf_nsets)write(10,*)'&'
  enddo ! sets

 case(2) ! plane

  do set=1,pcf_nsets
   if(iset==0)jset=set

   r(:)=apoint(:) ; r1(:)=apoint(:)
   do j=1,npoint2
    do i=1,npoint1
     expval=0.d0 ; expval2=0.d0
     call compute_fourier_basis(r)
     do g=1,expval_ngvec(sden_gset)
      expval=expval+real(expval_pcf(g,jset)*expval_expigdotr(g),dp)
      expval2=expval2+real(expval_sden(g,jset)*expval_expigdotr(g),dp)
     enddo
     write(10,'(e18.10,3(1x,e18.10))')r(:),expval/expval2
     r(:)=r(:)+step_ab(:)
    enddo
    write(10,*)
    r1(:)=r1(:)+step_ac(:) ; r(:)=r1(:)
   enddo

   if(iset==0.and.set<pcf_nsets)then
    write(10,*) ; write(10,*)
   endif
  enddo ! sets

 case(3) ! volume

  do set=1,pcf_nsets
   if(iset==0)jset=set

   do k=0,npoint3-1
    do j=0,npoint2-1
     do i=1,npoint1
      r(:)=apoint(:)+dble(i-1)*step_ab(:)+dble(j)*step_ac(:)+dble(k)*step_ad(:)
      expval=0.d0 ; expval2=0.d0
      call compute_fourier_basis(r)
      do g=1,expval_ngvec(pcf_gset)
       expval=expval+real(expval_pcf(g,jset)*expval_expigdotr(g),dp)
       expval2=expval2+real(expval_sden(g,jset)*expval_expigdotr(g),dp)
      enddo
      write(10,'(e18.10,3(1x,e18.10))')r(:),expval/expval2
     enddo
     write(10,*)
    enddo
   enddo

   if(iset==0.and.set<pcf_nsets)then
    write(10,*) ; write(10,*)
   endif
  enddo ! sets


 case default
  call errstop('PLOT_PAIR_CORR','Invalid dimensionality for plot.')
 end select

 END SUBROUTINE plot_pair_corr


 SUBROUTINE plot_pair_corr_sph
!-----------------------------------------------------------------------------!
! Plot the SPHERICAL PCF block from expval.data.                              !
!                                                                             !
! The value of the pair correlation function is collected into bins of width  !
! dr. At a point rn its value is calculated (for any type of particles) as    !
! the average over a simulation of the quantity                               !
!    V*N(n)                                                                   !
!   --------                                                                  !
!   N1*N2*Vn                                                                  !
! where N(n) is the number of pairs of particles in the nth bin (radius       !
! between n*dr and (n-1)*dr, Vn is the volume of the nth bin and N1,N2 are    !
! the number of particles of each type.                                       !
!                                                                             !
! For 3D:                                                                     !
!   Vn = 4*pi*dr^3*(n^2-n+1/3)                                                !
!   rn = 3*((n*dr)^4-((n-1)*dr)^4)/4/((n*dr)^3-((n-1)*dr)^3)                  !
! For 2D:                                                                     !
!   Vn = pi*dr^2*(2n-1)                                                       !
!   rn = 2*((n*dr)^3-((n-1)*dr)^3)/3/((n*dr)^2-((n-1)*dr)^2)                  !
! For 1D:                                                                     !
!   Vn = 2dr                                                                  !
!   rn = (n-1/2)dr                                                            !
!                                                                             !
! This routine takes the data N(n) in each bin, multiplies by the appropriate !
! factors and prints the result so that it can be plotted.                    !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,set,jset
 REAL(dp) const,dr,dn,rn,vn
 REAL(dp),PARAMETER :: third=1.d0/3.d0

!
! Define normalization (or not)
!
! Not sure what normalization to use in systems without a well-defined
! defined volume/area. Ordinarily we e.g divide by N/V so the PCF goes to
! 1 at long range at the value of the 'average density', but what to do
! in a molecule/polymer/slab geometry?
!
! Possibly choose the normalization of g so that the XC interaction
! energy comes out as
!
!                 n(r)n(r')
! \int dr \int dr' --------  [g(r,r')-1]
!                  |r-r'|
!
! Do later when someone actually wants to do this.

 if(expval_periodicity/=0)then
  if(expval_periodicity/=expval_dimensionality)call &
   &errstop('PLOT_PAIR_CORR_SPH','This subroutine requires the periodicity to &
   &equal the dimensionality.')
 else
  simcell_length=1.d0
  simcell_area=1.d0
  simcell_volume=1.d0
 endif

 if(reblocked)then

  dr=pcfs_rcutoff/real(pcfs_nbins_new,dp)
  if(iset/=0)then
   pcfs_nsets=1 ; jset=iset
  endif

  do set=1,pcfs_nsets

   if(iset==0)jset=set
   j=pcfs_ptype(1,jset) ; k=pcfs_ptype(2,jset)

   if(iset==0.and.(nele_expval(j)==0.or.nele_expval(k)==0&
    &.or.(j==k.and.nele_expval(j)==1)))cycle

   if(nele_expval(j)<=0.or.nele_expval(k)<=0)call &
    &errstop('PLOT_PAIR_CORR_SPH','No data in particle-pair set ' &
    &//trim(i2s(j))//'-'//trim(i2s(k))//'. You should only plot sets that &
    &contain data.')

   if(j==k)then
    const=1.d0/dble(nele_expval(j)*nele_expval(j))
   else
    const=0.5d0/dble(nele_expval(j)*nele_expval(k))
   endif
   select case(expval_dimensionality)
    case(1) ; const=const*simcell_length*0.5d0/dr
    case(2) ; const=const*simcell_area*0.5d0/pi/(dr*dr)
    case(3) ; const=const*simcell_volume*0.25d0/pi/(dr*dr*dr)
   end select

   select case (plot_dim)
    case(1)
     do i=1,pcfs_nbins_new
      dn=real(i,dp)
      select case(expval_dimensionality)
      case(1)
       vn=1.d0
       rn=(dn-0.5d0)*dr
      case(2)
       vn=dn-0.5d0
       rn=dr*(dn**2-dn+third)/vn
      case(3)
       vn=dn**2-dn+third
       rn=dr*(dn**3-1.5d0*dn**2+dn-0.25d0)/vn
      end select
      write(10,'(e18.10,1x,e18.10)')rn,expval_pcfs_reblocked(i,j,k)*const/vn
     enddo
    case default
     call errstop('PLOT_PAIR_CORR_SPH','Spherical PCF plots must be &
      &one-dimensional.')
   end select

   if(iset==0.and.set<pcfs_nsets)write(10,*)'&'

  enddo ! sets

 else ! no reblocking

  dr=pcfs_rcutoff/real(pcfs_nbins,dp)
  if(iset/=0)then
   pcfs_nsets=1 ; jset=iset
  endif

  do set=1,pcfs_nsets

   if(iset==0)jset=set
   j=pcfs_ptype(1,jset) ; k=pcfs_ptype(2,jset)

   if(iset==0.and.(nele_expval(j)==0.or.nele_expval(k)==0&
    &.or.(j==k.and.nele_expval(j)==1)))cycle

   if(j==k)then
    const=1.d0/dble(nele_expval(j)*nele_expval(j))
   else
    const=0.5d0/dble(nele_expval(j)*nele_expval(k))
   endif
   select case(expval_dimensionality)
   case(1) ; const=const*simcell_length*0.5d0/dr
   case(2) ; const=const*simcell_area*0.5d0/pi/(dr*dr)
   case(3) ; const=const*simcell_volume*0.25d0/pi/(dr*dr*dr)
   case default ; call errstop('PLOT_PAIR_CORR_SPH','Bug.')
   end select

   select case (plot_dim)
    case(1)
     do i=1,pcfs_nbins
      dn=real(i,dp)
      select case(expval_dimensionality)
       case(1)
        vn=1.d0
        rn=(dn-0.5d0)*dr
       case(2)
        vn=dn-0.5d0
        rn=dr*(dn**2-dn+third)/vn
       case(3)
        vn=dn**2-dn+third
        rn=dr*(dn**3-1.5d0*dn**2+dn-0.25d0)/vn
      end select
      write(10,'(e18.10,1x,e18.10)')rn,expval_pcfs(i,j,k)*const/vn
     enddo
    case default
     call errstop('PLOT_PAIR_CORR_SPH','Spherical PCF plots must be &
      &one-dimensional.')
    end select

   if(iset==0.and.set<pcfs_nsets)write(10,*)'&'

  enddo ! sets

 endif

 END SUBROUTINE plot_pair_corr_sph


 SUBROUTINE plot_localization_tensor
!------------------------------------------------------!
! Plot the LOCALIZATION TENSOR block from expval.data. !
!------------------------------------------------------!
 IMPLICIT NONE
 call errstop('PLOT_LOCALIZATION_TENSOR','Routine not yet coded')
 select case (plot_dim)
 case(1) ! line
  continue
 case(2) ! plane
  continue
 case(3) ! volume
  continue
 case default
  call errstop('PLOT_LOCALIZATION_TENSOR','Invalid dimensionality for plot.')
 end select
 END SUBROUTINE plot_localization_tensor


 SUBROUTINE plot_structure_factor
!---------------------------------------------------!
! Plot the STRUCTURE FACTOR block from expval.data. !
!---------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,jset,k,g,n,ng,nstar,set
 INTEGER,ALLOCATABLE :: first_g_in_star(:)
 REAL(dp) gmagsq,t1,t2
 REAL(dp) :: zero=0.d0
 REAL(dp),ALLOCATABLE :: gmagsq_star(:),sf_sum(:)
 LOGICAL print_this,single_pair,all_pairs,total

 select case(sf_option)

 case(1) ! spherical average

! Work out star structure of G vector set
  allocate(gmagsq_star(expval_ngvec(sf_gset)),&
   &first_g_in_star(expval_ngvec(sf_gset)+1),stat=ialloc)
  if(ialloc/=0)call errstop('PLOT_STRUCTURE_FACTOR','Allocation error : &
   &star information.')

  ng=1
  nstar=1
  gmagsq_star(1)=0.d0
  first_g_in_star(1)=1
  do g=2,expval_ngvec(sf_gset)
   gmagsq=expval_gvec(1,g,sf_gset)**2+expval_gvec(2,g,sf_gset)**2+&
    &expval_gvec(3,g,sf_gset)**2
   if(gmagsq-gmagsq_star(nstar)>1.d-8)then ! we have discovered a new star
    nstar=nstar+1
    first_g_in_star(nstar)=g
    gmagsq_star(nstar)=gmagsq
    ng=0
   endif
   ng=ng+1
  enddo ! G
  first_g_in_star(nstar+1)=first_g_in_star(nstar)+ng
  allocate(sf_sum(nstar),stat=ialloc)
  if(ialloc/=0)call errstop('PLOT_STRUCTURE_FACTOR','Allocation error : &
   &sf_sum array.')

! Spherically average the structure factor data over the G in each star,
! and write result to plot file.

  single_pair=(iset/=0.and.iset/=sf_nsets+1)
  all_pairs=(iset==0)
  total=(iset==sf_nsets+1)
  if(single_pair)then
   sf_nsets=1 ; jset=iset
  endif
  if(total)sf_sum(1:nstar)=0.d0

  if(homogeneous_structure_factor)then

   do set=1,sf_nsets

    if(all_pairs.or.total)jset=set
    j=sf_ptype(1,jset) ; k=sf_ptype(2,jset)
    if(j/=k)then
     t2=2.d0
    else
     t2=1.d0
    endif
    print_this=(.not.total.or.(total.and.set==sf_nsets))
    if(.not.total)sf_sum(1:nstar)=0.d0

    sf_sum(1)=sf_sum(1)+t2*(expval_sf(1,jset)-nele_expval(j)*nele_expval(k))
    if(print_this)write(10,'(e18.10,1x,e18.10)')zero,sf_sum(1)

    g=1
    do n=2,nstar
     ng=0
     do i=first_g_in_star(n),first_g_in_star(n+1)-1,2
      g=g+1
      sf_sum(n)=sf_sum(n)+t2*(expval_sf(g,jset)+expval_sf(g+1,jset))
      g=g+1
      ng=ng+2
     enddo ! G in star
     if(print_this)then
      if(ng/=0)then
        t1=(sf_sum(n)*real(npcells**2,dp))/(real(ng,dp)*real(netot,dp))
        write(10,'(e18.10,1x,e18.10)')sqrt(gmagsq_star(n)),t1
      else
       call errstop('PLOT_STRUCTURE_FACTOR','Zero G vectors in star - should &
        &not happen.')
      endif
     endif
    enddo ! stars

    if(all_pairs.and.set<sf_nsets)write(10,*)'&'

   enddo ! SF sets

  else ! inhomogeneous

   do set=1,sf_nsets

    if(all_pairs.or.total)jset=set
    j=sf_ptype(1,jset) ; k=sf_ptype(2,jset)
    if(j/=k)then
     t2=2.d0
    else
     t2=1.d0
    endif
    print_this=(.not.total.or.(total.and.set==sf_nsets))
    if(.not.total)sf_sum(1:nstar)=0.d0

    sf_sum(1)=sf_sum(1)+t2*(expval_sf(1,jset)- &
     &real(expval_sfsden(1,j)*expval_sfsden(1,k),dp)) ! i.e. zero
    if(print_this)write(10,'(e18.10,1x,e18.10)')zero,sf_sum(1)

    g=1
    do n=2,nstar
     ng=0
     do i=first_g_in_star(n),first_g_in_star(n+1)-1,2
      g=g+1
      sf_sum(n)=sf_sum(n)+t2*(expval_sf(g,jset)+expval_sf(g+1,jset)-&
       &2.d0*real(expval_sfsden(g,j)*expval_sfsden(g+1,k),dp))
      g=g+1
      ng=ng+2
     enddo ! G in star

     if(print_this)then
      if(ng/=0)then
       t1=(sf_sum(n)*real(npcells**2,dp))/(real(ng,dp)*real(netot,dp))
      else
       call errstop('PLOT_STRUCTURE_FACTOR','Zero G vectors in star - should &
        &not happen.')
      endif
      write(10,'(e18.10,1x,e18.10)')sqrt(gmagsq_star(n)),t1
     endif

    enddo ! stars

    if(all_pairs.and.set<sf_nsets)write(10,*)'&'

   enddo ! SF sets

  endif ! homogeneous density or not

 case(2) ! raw data

  call errstop('PLOT_STRUCTURE_FACTOR','Plotting option 2 for structure factor&
   & not yet implemented.')

  select case (plot_dim)
  case(1) ! line
   continue
  case(2) ! plane
   continue
  case(3) ! volume
   continue
  case default
   call errstop('PLOT_STRUCTURE_FACTOR','Invalid dimensionality for plot.')
  end select

 end select

 deallocate(gmagsq_star,first_g_in_star)

 END SUBROUTINE plot_structure_factor


 SUBROUTINE plot_structure_factor_sph
!-------------------------------------------------------------!
! Plot the SPHERICAL STRUCTURE FACTOR block from expval.data. !
!-------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER k,set

 select case(plot_dim)

 case(1) ! line

  if(iset==0)then
   do set=1,sf_sph_nsets
    if(set>1)write(10,*)'&'
    do k=1,sf_sph_nk
     write(10,'(e18.10,1x,e18.10)')expval_sf_sph_k(k),expval_sf_sph(k,set)
    enddo ! k
   enddo ! set
  elseif(iset==sf_sph_nsets+1)then
   do k=1,sf_sph_nk
    write(10,'(e18.10,1x,e18.10)')expval_sf_sph_k(k),sum(expval_sf_sph(k,:))
   enddo ! k
  else
   do k=1,sf_sph_nk
    write(10,'(e18.10,1x,e18.10)')expval_sf_sph_k(k),expval_sf_sph(k,iset)
   enddo ! k
  endif ! iset

 case default
  call errstop('PLOT_STRUCTURE_FACTOR_SPH','Spherical structure factor plots &
   &must be one-dimensional.')
 end select

 END SUBROUTINE plot_structure_factor_sph


 SUBROUTINE plot_onep_density_matrix
!--------------------------------------------------------------!
! Plot the ONE-PARTICLE DENSITY MATRIX block from expval.data. !
!--------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i
 REAL(dp) r,delta,x,dx

 delta=wigner_seitz_radius/real(onep_dm_nbins,dp)

 select case(plot_dim)
 case(1) ! line
  do i=1,onep_dm_nbins
   select case(expval_periodicity)
   case(1) ; r=real(2*i-1,dp)*delta*.5d0
   case(2) ; r=real(6*i*(i-1)+2,dp)*delta/real(6*i-3,dp)
   case(3) ; r=real(6*i*(i*(2*i-3)+2)-3,dp)*delta/real(12*i*(i-1)+4,dp)
   end select
   x=expval_onep_dm(i,iset) ; dx=0.d0
   if(expval_onep_dm_weight(i,iset)>1.d0)then
    dx=sqrt(abs(expval_onep_dm2(i,iset)-expval_onep_dm(i,iset)**2)/&
     &(expval_onep_dm_weight(i,iset)-1.d0))
   endif
   write(10,'(e18.10,1x,e18.10,1x,e18.10)')r,x,dx
  enddo
 case default
  call errstop('PLOT_ONEP_DENSITY_MATRIX','One-particle density matrix plots &
   &must be one-dimensional.')
 end select

 END SUBROUTINE plot_onep_density_matrix


 SUBROUTINE plot_twop_density_matrix
!--------------------------------------------------------------!
! Plot the TWO-PARTICLE DENSITY MATRIX block from expval.data. !
!--------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,n,spin1,set1,oset(2,twop_dm_nptypes_in_set(iset))
 REAL(dp) r,delta,x,dx,ox,ox1,ox2,odx_2,odx1_2,odx2_2,t1,&
  &ofactor(2,twop_dm_nptypes_in_set(iset))
 LOGICAL pairmat(num_particle_types,num_particle_types)

! Bin width.
 delta=wigner_seitz_radius/real(twop_dm_nbins,dp)

! Prepare for subtracting one-body contributions using OBDM data.
 if(twop_dm_with_onep_dm)then
! Check grid consistency.
  if(twop_dm_nbins/=onep_dm_nbins)call errstop('PLOT_TWOP_DENSITY_MATRIX',&
   &'Must use same number of bins for TBDM and OBDM.')
! Find sets matching each particle pair in the TBDM set, and construct
! rescaling factors for each pair.
  oset=0
  ofactor=0.d0
  do i=1,twop_dm_nptypes_in_set(iset)
   do j=1,2
    do set1=1,onep_dm_nsets
     n=0
     do k=1,onep_dm_nptypes_in_set(set1)
      spin1=onep_dm_ptype_in_set(k,set1)
      n=n+nele_expval(spin1)
      if(spin1==twop_dm_ptype_in_set(j,i,iset))oset(j,i)=set1
     enddo ! k
     if(oset(j,i)/=0)exit
    enddo ! set1
    if(oset(j,i)==0)call errstop('PLOT_TWOP_DENSITY_MATRIX',&
     &'Could not find matching OBDM sets to subtract from TBDM.')
    if(n<1)call errstop('PLOT_TWOP_DENSITY_MATRIX',&
     &'Problem counting particles in OBDM set.')
    ofactor(j,i)=dble(nele_expval(twop_dm_ptype_in_set(j,i,iset)))/dble(n)
   enddo ! j
   if(twop_dm_ptype_in_set(1,i,iset)==twop_dm_ptype_in_set(2,i,iset))&
    &ofactor(2,i)=ofactor(2,i)*&
    &dble(nele_expval(twop_dm_ptype_in_set(1,i,iset))-1)/&
    &dble(nele_expval(twop_dm_ptype_in_set(1,i,iset)))
  enddo ! i
 endif ! twop_dm_with_onep_dm

! Calculate rescaling factor so lim r->inf is condensate fraction.
 t1=1.d0
 if(twop_dm_rescale)then
  ! Flag which spins are paired.
  pairmat=.false.
  do i=1,twop_dm_nptypes_in_set(iset)
   pairmat(twop_dm_ptype_in_set(1,i,iset),twop_dm_ptype_in_set(2,i,iset))=.true.
   pairmat(twop_dm_ptype_in_set(2,i,iset),twop_dm_ptype_in_set(1,i,iset))=.true.
  enddo ! i
  ! Number of pairs given by max number of particles paired to any particle
  ! type.
  n=0
  do i=1,num_particle_types
   n=max(n,sum(nele_expval,pairmat(:,i)))
  enddo ! i
  ! Evaluate prefactor.
  select case(expval_periodicity)
  case(1) ; t1=(simcell_length**2)/real(n,dp)
  case(2) ; t1=(simcell_area**2)/real(n,dp)
  case(3) ; t1=(simcell_volume**2)/real(n,dp)
  end select
 endif ! twop_dm_rescale

 select case(plot_dim)
 case(1) ! line
  do i=1,twop_dm_nbins
   select case(expval_periodicity)
   case(1) ; r=real(2*i-1,dp)*delta*.5d0
   case(2) ; r=real(6*i*(i-1)+2,dp)*delta/real(6*i-3,dp)
   case(3) ; r=real(6*i*(i*(2*i-3)+2)-3,dp)*delta/real(12*i*(i-1)+4,dp)
   end select
   x=expval_twop_dm(i,iset)
   dx=0.d0
   if(expval_twop_dm_weight(i,iset)>1.d0)dx=&
    &sqrt(abs(expval_twop_dm2(i,iset)-expval_twop_dm(i,iset)**2)/&
    &(expval_twop_dm_weight(i,iset)-1.d0))
   if(twop_dm_with_onep_dm)then
    ox=0.d0
    odx_2=0.d0
    do j=1,twop_dm_nptypes_in_set(iset)
     ox1=ofactor(1,j)*expval_onep_dm(i,oset(1,j))
     ox2=ofactor(2,j)*expval_onep_dm(i,oset(2,j))
     ox=ox+ox1*ox2
     if(expval_onep_dm_weight(i,oset(1,j))>1.d0.and.&
      &expval_onep_dm_weight(i,oset(2,j))>1.d0)then
      odx1_2=(expval_onep_dm2(i,oset(1,j))-&
       &expval_onep_dm(i,oset(1,j))**2)/&
       &(expval_onep_dm_weight(i,oset(1,j))-1.d0)
      odx2_2=(expval_onep_dm2(i,oset(2,j))-&
       &expval_onep_dm(i,oset(2,j))**2)/&
       &(expval_onep_dm_weight(i,oset(2,j))-1.d0)
      odx_2=odx_2+ox1**2*odx2_2+ox2**2*odx1_2
     endif
    enddo ! j
    x=x-ox
    dx=sqrt(abs(dx**2+odx_2))
   endif
   if(twop_dm_rescale)then
    x=t1*x
    dx=t1*dx
   endif
   write(10,'(e18.10,1x,e18.10,1x,e18.10)')r,x,dx
  enddo
 case default
  call errstop('PLOT_TWOP_DENSITY_MATRIX','Two-particle density matrix plots &
   &must be one-dimensional.')
 end select

 END SUBROUTINE plot_twop_density_matrix


 SUBROUTINE plot_cond_fraction
!------------------------------------------------------!
! Plot the CONDENSATE FRACTION block from expval.data. !
!------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,n
 REAL(dp) r,delta,x,dx,t1
 LOGICAL pairmat(num_particle_types,num_particle_types)

 delta=wigner_seitz_radius/real(cond_frac_nbins,dp)
 ! Flag which spins are paired.
 pairmat=.false.
 do i=1,twop_dm_nptypes_in_set(iset)
  pairmat(twop_dm_ptype_in_set(1,i,iset),twop_dm_ptype_in_set(2,i,iset))=.true.
  pairmat(twop_dm_ptype_in_set(2,i,iset),twop_dm_ptype_in_set(1,i,iset))=.true.
 enddo ! i
 ! Number of pairs given by max number of particles paired to any particle
 ! type.
 n=0
 do i=1,num_particle_types
  n=max(n,sum(nele_expval,pairmat(:,i)))
 enddo ! i
 ! Evaluate prefactor.
 select case(expval_periodicity)
 case(1) ; t1=(simcell_length**2)/real(n,dp)
 case(2) ; t1=(simcell_area**2)/real(n,dp)
 case(3) ; t1=(simcell_volume**2)/real(n,dp)
 end select

 select case(plot_dim)
 case(1) ! line
  do i=1,cond_frac_nbins
   select case(expval_periodicity)
   case(1) ; r=real(2*i-1,dp)*delta*.5d0
   case(2) ; r=real(6*i*(i-1)+2,dp)*delta/real(6*i-3,dp)
   case(3) ; r=real(6*i*(i*(2*i-3)+2)-3,dp)*delta/real(12*i*(i-1)+4,dp)
   end select
   dx=0.d0 ; x=expval_cond_frac(i,iset)*t1
   if(expval_cond_frac_weight(i,iset)>1.d0)then
    dx=sqrt(abs(expval_cond_frac2(i,iset)-expval_cond_frac(i,iset)**2)/&
     &(expval_cond_frac_weight(i,iset)-1.d0))
    dx=dx*t1
   endif
   write(10,'(e18.10,1x,e18.10,1x,e18.10)')r,x,dx
  enddo
 case default
  call errstop('PLOT_COND_FRACTION','Condensate fraction estimator plots &
   &must be one-dimensional.')
 end select

 END SUBROUTINE plot_cond_fraction


 SUBROUTINE plot_mom_den
!---------------------------------------------------!
! Plot the MOMENTUM DENSITY block from expval.data. !
!---------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,nsame
 REAL(dp) x,x2,w,k,dx,inv_nsame,rs,kf,k_factor,md_factor
 REAL(dp),PARAMETER :: tol_same=1.d-12

 select case(plot_dim)
 case(1) ! line
! Sort k points by size (only necessary if k_offset/=0).
  do i=1,expval_ngvec(1)
   k=expval_mom_den_k(i,iset)
   do j=i+1,expval_ngvec(1)
    if(k>expval_mom_den_k(j,iset))then
     x=expval_mom_den(j,iset) ; x2=expval_mom_den2(j,iset)
     k=expval_mom_den_k(j,iset)
     expval_mom_den(j,iset)=expval_mom_den(i,iset)
     expval_mom_den2(j,iset)=expval_mom_den2(i,iset)
     expval_mom_den_k(j,iset)=expval_mom_den_k(i,iset)
     expval_mom_den(i,iset)=x
     expval_mom_den2(i,iset)=x2
     expval_mom_den_k(i,iset)=k
    endif
   enddo
  enddo
! Rescaling factors for on-top momentum density.
  k_factor=1.d0
  md_factor=1.d0
  if(on_top_ispin/=0)then
   if(nele_expval(on_top_ispin)==1.or.nele_expval(on_top_jspin)==1)then
    write(6,*)'On-top momentum density detected.'
    write(6,*)'Rescaling so that momentum density integrates to pi * kF^2'
    write(6,*)'and plotting against k/kF .'
    write(6,*)
    select case(expval_periodicity)
    case(1)
     rs=simcell_length/(2.d0*dble(netot-1))
     kf=pi/(4.d0*rs)
     k_factor=1.d0/kf
     md_factor=pi/(4.d0*rs)
    case(2)
     rs=sqrt(simcell_area/(pi*dble(netot-1)))
     kf=sqrt(2.d0)/rs
     k_factor=1.d0/kf
     md_factor=2.d0*pi/rs**2
    case(3)
     rs=(3.d0*simcell_volume/(4.d0*pi*dble(netot-1)))**(1.d0/3.d0)
     kf=(9.d0*pi/4.d0)**(1.d0/3.d0)/rs
     k_factor=1.d0/kf
     md_factor=3.d0*pi**2/rs**3
    end select
   endif
  endif
  x=0.d0 ; x2=0.d0 ; w=0.d0 ; nsame=0
  do i=1,expval_ngvec(1)
   x=x+expval_mom_den(i,iset) ; x2=x2+expval_mom_den2(i,iset)
   w=w+expval_mom_den_weight(iset) ; nsame=nsame+1
   if(i<expval_ngvec(1))then
    if(abs(expval_mom_den_k(i+1,iset)-expval_mom_den_k(i,iset))<tol_same)cycle
   endif
   inv_nsame=1.d0/dble(nsame) ; x=x*inv_nsame ; x2=x2*inv_nsame ; w=w*inv_nsame
   dx=0.d0 ; if(w>1.d0)dx=sqrt(abs(x2-x*x)/(w-1.d0))
   write(10,'(e18.10,1x,e18.10,1x,e18.10)')&
    &expval_mom_den_k(i,iset)*k_factor,x*md_factor,dx*md_factor
   x=0.d0 ; x2=0.d0 ; w=0.d0 ; nsame=0
  enddo
 case default
  call errstop('PLOT_MOM_DEN','Momentum density plots must be one-dimensional.')
 end select

 END SUBROUTINE plot_mom_den


 SUBROUTINE plot_finite_density
!-------------------------------------------------!
! Plot the FINITE DENSITY block from expval.data. !
!-------------------------------------------------!
 IMPLICIT NONE
 REAL(dp) alpha,beta,exp_step,gamma_cusp,k_0
 INTEGER npoints,i

 if(fin_den_ibasis==fin_den_f_indx_step)then
  npoints=fin_den_norder
 else
  npoints=1000
 endif
 alpha=sqrt(2.d0*fin_den_ip)
 beta=dble(fin_den_z-netot+1)/alpha-1.d0
 exp_step=1.6d0
 k_0=exp_step**(-10)*(dble(netot)-5.d0/16.d0)

 open(unit=io_rad_mom,file='rad_mom.dat',status='new',iostat=ierr)
 if(ierr/=0)call errstop('PLOT_FINITE_DENSITY','Problem opening rad_mom.data.')

 if(iset==1)then ! Intracule density

  do i=1,fin_den_nsets_ij
   if(i==1.or.i==3)gamma_cusp=0.5d0
   if(i==2)gamma_cusp=0.25d0
   write(10,*)'# h(r_ij) set ',i
   call fit_fin_den(expval_fin_den_ij(:,i),expval_fin_den2_ij(:,i),gamma_cusp,&
    &2.d0*fin_den_cutoff,2.d0*fin_den_tail_offset,fin_den_ip,alpha,beta,&
    &exp_step,k_0,expval_fin_den_nstep,expval_fin_den_weight,&
    &expval_fin_den_weight2,npoints)
  enddo
  write(10,*)'# h(r_ij) total'
  gamma_cusp=0.25d0
  call fit_fin_den(expval_fin_den_ij(:,1)+expval_fin_den_ij(:,2)+&
   &expval_fin_den_ij(:,3),expval_fin_den2_ij(:,1)+expval_fin_den2_ij(:,2)+&
   &expval_fin_den2_ij(:,3),gamma_cusp,2.d0*fin_den_cutoff,&
   &2.d0*fin_den_tail_offset,fin_den_ip,alpha,beta,exp_step,k_0,&
   &expval_fin_den_nstep,expval_fin_den_weight,expval_fin_den_weight2,npoints)

 elseif(iset==2)then ! Finite density

  gamma_cusp=-dble(fin_den_z)
  do i=1,fin_den_nsets_i
   write(10,*)'# rho(r_i) set ',i
   call fit_fin_den(expval_fin_den_i(:,i),expval_fin_den2_i(:,i),gamma_cusp,&
    &fin_den_cutoff,fin_den_tail_offset,fin_den_ip,alpha,beta,exp_step,k_0,&
    &expval_fin_den_nstep,expval_fin_den_weight,expval_fin_den_weight2,npoints)
  enddo
  write(10,*)'# rho(r_i) total'
  call fit_fin_den(expval_fin_den_i(:,1)+expval_fin_den_i(:,2),&
   &expval_fin_den2_i(:,1)+expval_fin_den2_i(:,2),gamma_cusp,&
   &fin_den_cutoff,fin_den_tail_offset,fin_den_ip,alpha,beta,exp_step,k_0,&
   &expval_fin_den_nstep,expval_fin_den_weight,expval_fin_den_weight2,npoints)
 endif

 close(io_rad_mom)

 END SUBROUTINE plot_finite_density


 SUBROUTINE fit_fin_den(accum,accum2,gamma_cusp,r_c,L,fin_den_ip,alpha,beta,&
  &exp_step,k_0,nstep,w,w2,npoints)
!-----------------------------------------------------------!
! Fit FINITE DENSITY for a single e-e or e-n set and write. !
!-----------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: npoints
 REAL(dp),INTENT(in) :: accum(:),accum2(:),gamma_cusp,r_c,L,fin_den_ip,alpha,&
  &beta,exp_step,k_0,nstep,w,w2
 INTEGER n_min,n_max,i,j,k,n,ierr,ialloc
 INTEGER,PARAMETER :: nconstraints=4
 INTEGER,ALLOCATABLE :: piv_h(:),piv_g(:)
 REAL(dp),ALLOCATABLE :: h(:,:,:),h_0(:,:),h_0_inv(:,:),h_pn(:,:),coeff(:),&
  &dcoeff2(:),daccum2(:),coeff_const(:),dcoeff2_const(:),phi(:),dphi(:),&
  &constraint1(:,:),constraint2(:),gam(:,:),gam_inv(:,:),lambda(:),&
  &dlambda2(:),f0(:),df0(:),fL(:),dfL(:),d2fL(:),rad_mom(:),drad_mom2(:)
 REAL(dp) omega,t1,t2,t3,f,df2,dfdr,d_dfdr2,r
 CHARACTER(80) char_80
 CHARACTER(256) char_256
 EXTERNAL dgetrf,dgetrs
 ! Fitting with constraints and radial moments evaluation is not implemented
 ! for the step function basis. Use nl2sol to fit binned data.
 LOGICAL,PARAMETER :: constraint=.false.,calc_rad_mom=.false.

! Range over which radial moments are to be calculated.
 n_min=minval(fin_den_power)
 n_max=maxval(fin_den_power)

 allocate(&
  &h(fin_den_norder,fin_den_norder,n_min:n_max),&
  &h_0(fin_den_norder,fin_den_norder),&
  &h_0_inv(fin_den_norder,fin_den_norder),piv_h(fin_den_norder),&
  &h_pn(fin_den_norder,n_min:n_max),&
  &coeff(fin_den_norder),dcoeff2(fin_den_norder),daccum2(fin_den_norder),&
  &coeff_const(fin_den_norder),dcoeff2_const(fin_den_norder),&
  &phi(fin_den_norder),dphi(fin_den_norder),&
  &constraint1(nconstraints,fin_den_norder),&
  &constraint2(nconstraints),&
  &gam(nconstraints,nconstraints),&
  &gam_inv(nconstraints,nconstraints),piv_g(nconstraints),&
  &lambda(nconstraints),dlambda2(nconstraints),&
  &f0(fin_den_norder),df0(fin_den_norder),&
  &fL(fin_den_norder),dfL(fin_den_norder),d2fL(fin_den_norder),&
  &rad_mom(n_min:n_max),drad_mom2(n_min:n_max),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('FIT_FIN_DEN','Allocation error.')

! Evaluate overlap matrix H (basis specific).
 h_0=0.d0
 select case(fin_den_ibasis)
 case(fin_den_f_indx_step)
  do i=1,fin_den_norder
   h_0(i,i)=(4.d0*pi/3.d0)*(dble(i)**3-dble(i-1)**3)*r_c**3/&
    &dble(fin_den_norder)**3 ! normalise volume
  enddo
 case(fin_den_f_indx_bessel)
  do i=1,fin_den_norder
   h_0(i,i)=1.d0
  enddo
 case(fin_den_f_indx_cheb)
  call eval_cheb_overlap(fin_den_norder,n_min,n_max,L,h)
  h_0(1:fin_den_norder,1:fin_den_norder)=h(1:fin_den_norder,1:fin_den_norder,0)
 case(fin_den_f_indx_exp)
  call eval_exp_overlap(fin_den_norder,k_0,exp_step,L,h_0)
 case default
  call errstop('FIT_FIN_DEN','Invalid basis set.')
 end select

! Evaluate inverse of H.
 h_0_inv=0.d0
 do i=1,fin_den_norder
  h_0_inv(i,i)=1.d0
 enddo
 call dgetrf(fin_den_norder,fin_den_norder,h_0,fin_den_norder,piv_h,ierr)
 if(ierr/=0)call errstop('FIT_FIN_DEN','h_0_inv: DGETRF says ierr = '//&
  &trim(i2s(ierr))//'.')
 call dgetrs('N',fin_den_norder,fin_den_norder,h_0,fin_den_norder,piv_h,&
  &h_0_inv,fin_den_norder,ierr)
 if(ierr/=0)call errstop('FIT_FIN_DEN','h_0_inv: DGETRS says ierr = '//&
  &trim(i2s(ierr))//'.')

! Get daccum2.
 do i=1,fin_den_norder
  daccum2(i)=abs(accum2(i)-accum(i)**2)/((1.d0-w2/w**2)*dble(nstep))
 enddo

! Calculate coeffs and dcoeffs of each basis function.
 coeff=0.d0 ; dcoeff2=0.d0
 do i=1,fin_den_norder
  coeff(i)=sum(h_0_inv(i,:)*accum(:))
  dcoeff2(i)=sum(h_0_inv(i,:)**2*daccum2(:))
 enddo

 if(constraint)then
  write(6,*)"Constraints are being applied."

! Calculate constraints.
  constraint1=0.d0 ; constraint2=0.d0
  constraint2=(/0.d0,0.d0,0.d0,1.d0/)
  call fin_den_basis_fnc(0.d0,fin_den_norder,r_c,k_0,exp_step,f=f0,df=df0)
  call fin_den_basis_fnc(L,fin_den_norder,r_c,k_0,exp_step,f=fL,df=dfL,d2f=d2fL)
  do i=1,fin_den_norder
   if(fin_den_ibasis==fin_den_f_indx_step)then
    call errstop('PLOT_FINITE_DENSITY','Constraints cannot be applied in the &
     &case of a step function basis.')
   elseif(fin_den_ibasis/=fin_den_f_indx_step)then
    constraint1(1,i)=df0(i)-2.d0*gamma_cusp*f0(i)
    constraint1(2,i)=dfL(i)-fL(i)*(2.d0*beta/L-2.d0*alpha)
    constraint1(3,i)=d2fL(i)-fL(i)*(2.d0*beta*(2.d0*beta-1.d0)/(L*L)-&
     &8.d0*beta*alpha/L+8.d0*fin_den_ip)
    select case(fin_den_ibasis)
    case(fin_den_f_indx_bessel)
     omega=dble(i)*pi/r_c
     t1=2.d0*twopi*(sin(omega*L)/(omega*omega)-L*cos(omega*L)/omega)/&
      &sqrt(twopi*r_c)
    case(fin_den_f_indx_cheb)
     t1=h(i,1,0)
    case(fin_den_f_indx_exp)
     call eval_exp_sph_ave(fin_den_norder,n_min,n_max,k_0,exp_step,L,h_pn)
     t1=h_pn(i,0)
    case default
     call errstop('PLOT_FINITE_DENSITY','Invalid basis set.')
    end select
    t2=fL(i)*exp(2.d0*alpha*L)/((2.d0*alpha)**(2.d0*beta+3.d0)*L**&
     &(2.d0*beta))
    t3=t2*2.d0*twopi*gamma_inc_upper(2.d0*beta+3.d0,2.d0*alpha*L)
    constraint1(4,i)=t1+t3
   endif
  enddo

! Calculate Lagrange multipliers.
  lambda=0.d0 ; dlambda2=0.d0
  gam=0.d0 ; gam_inv=0.d0
  do i=1,nconstraints
   do j=1,nconstraints
    do k=1,fin_den_norder
     gam(i,j)=gam(i,j)+constraint1(i,k)*sum(h_0_inv(k,:)*constraint1(j,:))
    enddo
   enddo
  enddo
  do i=1,nconstraints
   gam_inv(i,i)=1.d0
  enddo
  call dgetrf(nconstraints,nconstraints,gam,nconstraints,piv_g,ierr)
  if(ierr/=0)call errstop('PLOT_FINITE_DENSITY','gam_inv: DGETRF says ierr = '&
   &//trim(i2s(ierr))//'.')
  call dgetrs('N',nconstraints,nconstraints,gam,nconstraints,piv_g,&
   &gam_inv,nconstraints,ierr)
  if(ierr/=0)call errstop('PLOT_FINITE_DENSITY','gam_inv: DGETRS says ierr = '&
   &//trim(i2s(ierr))//'.')
  do i=1,nconstraints
   do j=1,nconstraints
    lambda(i)=lambda(i)+gam_inv(i,j)*(sum(constraint1(j,:)*coeff(:))-&
     &constraint2(j))
    dlambda2(i)=dlambda2(i)+gam_inv(i,j)**2*&
     &(sum(constraint1(j,:)**2*dcoeff2(:)))
   enddo
  enddo

! Calculate constrained coefficients.
  coeff_const=coeff ; dcoeff2_const=dcoeff2
  do i=1,fin_den_norder
   do j=1,fin_den_norder
    coeff_const(i)=coeff_const(i)-h_0_inv(i,j)*sum(constraint1(:,j)*lambda(:))
    dcoeff2_const(i)=dcoeff2_const(i)+h_0_inv(i,j)**2*sum(constraint1(:,j)**2*&
     &dlambda2(:))
   enddo
  enddo

 else ! no constraints

  coeff_const=coeff
  dcoeff2_const=dcoeff2

 endif

! Calculate radial moments.
! Note: not implemented for step functions (use fitted form) and spherical
! bessel functions.
 if(calc_rad_mom)then
  rad_mom=0.d0 ; drad_mom2=0.d0
  do n=n_min,n_max,1
   t1=exp(-2.d0*alpha*L)*((2.d0*alpha)**(2.d0*beta+dble(n)+3.d0)*L**&
    &(2.d0*beta))
   t2=2.d0*twopi*gamma_inc_upper(2.d0*beta+dble(n)+3.d0,2.d0*alpha*L)/t1
   select case(fin_den_ibasis)
   case(fin_den_f_indx_step)
    call errstop('PLOT_FINITE_DENSITY','Radial moments cannot be calculated &
     &for step function basis.')
   case(fin_den_f_indx_bessel)
    call errstop('PLOT_FINITE_DENSITY','Radial moments cannot be calculated &
     &for sph bessel function basis.')
   case(fin_den_f_indx_cheb)
    rad_mom(n)=sum(h(:,1,n)*coeff_const(:))+sum(fL(:)*coeff_const(:))*t2
    drad_mom2(n)=sum(h(:,1,n)**2*dcoeff2_const(:))+sum(fL(:)**2*&
     &dcoeff2_const(:))*t2**2
   case(fin_den_f_indx_exp)
    rad_mom(n)=sum(h_pn(:,n)*coeff_const(:))+sum(fL(:)*coeff_const(:))*t2
    drad_mom2(n)=sum(h_pn(:,n)**2*dcoeff2_const(:))+sum(fL(:)**2*&
     &dcoeff2_const(:))*t2**2
   end select
  enddo
!  Write radial moments.
  if((fin_den_ibasis==fin_den_f_indx_cheb).or.&
    &(fin_den_ibasis==fin_den_f_indx_exp))then
   do n=n_min,n_max,1
    char_256=''
    write(char_80,*)n
    char_256=trim(char_256)//' '//trim(adjustl(char_80))
    write(char_80,*)rad_mom(n)
    char_256=trim(char_256)//' '//trim(adjustl(char_80))
    write(char_80,*)sqrt(drad_mom2(n))
    char_256=trim(char_256)//' '//trim(adjustl(char_80))
    write(io_rad_mom,'(a)')trim(char_256)
   enddo
  endif
 endif

! Calculate binned data and radial moments for step function data.
 rad_mom=0.d0
 drad_mom2=0.d0
 do i=0,npoints
  if(fin_den_ibasis==fin_den_f_indx_step)then
   if(i==0)cycle
   t1=dble(i)
   r=(6*t1*(t1*(2*t1-3)+2)-3)/(12*t1*(t1-1)+4)*r_c/dble(npoints) ! bin "centre"
  else
   r=dble(i)*r_c/dble(npoints)
  endif
  f=0.d0 ; df2=0.d0 ; dfdr=0.d0 ; d_dfdr2=0.d0
  call fin_den_basis_fnc(r,fin_den_norder,r_c,k_0,exp_step,f=phi,df=dphi)
  do j=1,fin_den_norder
   f=f+coeff_const(j)*phi(j)
   df2=df2+dcoeff2_const(j)*phi(j)*phi(j)
   dfdr=dfdr+coeff_const(j)*dphi(j)
   d_dfdr2=d_dfdr2+dcoeff2_const(j)*dphi(j)*dphi(j)
  enddo ! norder
! Radial moments
  if(i==0.or.i==npoints)then
   do n=n_min,n_max,1
    rad_mom(n)=rad_mom(n)+r**(n+2)*f/2.d0
    drad_mom2(n)=drad_mom2(n)+r**(2*n+4)*df2/4.d0
   enddo
  else
   do n=n_min,n_max,1
    rad_mom(n)=rad_mom(n)+r**(n+2)*f
    drad_mom2(n)=drad_mom2(n)+r**(2*n+4)*df2
   enddo
  endif
! Write data.
  char_256=''
  write(char_80,*)r
  char_256=trim(char_256)//' '//trim(adjustl(char_80))
  write(char_80,*)f
  char_256=trim(char_256)//' '//trim(adjustl(char_80))
  write(char_80,*)sqrt(df2)
  char_256=trim(char_256)//' '//trim(adjustl(char_80))
  write(char_80,*)dfdr
  char_256=trim(char_256)//' '//trim(adjustl(char_80))
  write(char_80,*)sqrt(d_dfdr2)
  char_256=trim(char_256)//' '//trim(adjustl(char_80))
  write(10,'(a)')trim(char_256)
 enddo ! r values
 write(10,*)'&'
 char_256=''
 write(char_80,*)netot
 char_256=trim(char_256)//' '//trim(adjustl(char_80))
 do n=n_min,n_max,1
  rad_mom(n)=4.d0*pi*r_c/dble(npoints)*rad_mom(n)
  drad_mom2(n)=4.d0*pi*r_c/dble(npoints)*sqrt(drad_mom2(n)) ! d_rad_mom
  write(char_80,*)rad_mom(n) ! rad_mom
  char_256=trim(char_256)//' '//trim(adjustl(char_80))
  write(char_80,*)drad_mom2(n) ! d_rad_mom
  char_256=trim(char_256)//' '//trim(adjustl(char_80))
 enddo
 write(io_rad_mom,'(a)')trim(char_256)

 deallocate(h,h_0,h_0_inv,piv_h,h_pn,coeff,dcoeff2,daccum2,coeff_const,&
  &dcoeff2_const,phi,dphi,constraint1,constraint2,gam,gam_inv,&
  &lambda,dlambda2,f0,df0,fL,dfL,d2fL,rad_mom,drad_mom2,stat=ialloc)
 if(ialloc/=0)call errstop('FIT_FIN_DEN','Deallocation error.')

 END SUBROUTINE fit_fin_den


 SUBROUTINE fin_den_basis_fnc(r,n,r_c,k_0,exp_step,f,df,d2f)
!-------------------------------------------------------------------------!
! This subroutine acts as a wrapper around the different basis functions. !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 REAL(dp),INTENT(in) :: r,r_c,k_0,exp_step
 REAL(dp),INTENT(inout),OPTIONAL :: f(n),df(n),d2f(n)
 INTEGER i
 REAL(dp) binsize
 LOGICAL val,fd,sd

 val=present(f)
 fd=present(df)
 sd=present(d2f)

 select case(fin_den_ibasis)
 case(fin_den_f_indx_step) ! step
  do i=1,n
   if(val)then
    binsize=r_c/dble(n)
    if(r>=dble(i-1)*binsize.and.r<dble(i)*binsize)then
     f(i)=1.d0
    else
     f(i)=0.d0
    endif
   endif
  enddo
  if(fd.or.sd)df=0.d0
  if(sd)d2f=0.d0
 case(fin_den_f_indx_bessel) ! bessel
  if(sd)then
   call sph_bessel(r,n,r_c,f,df,d2f)
  elseif(fd)then
   call sph_bessel(r,n,r_c,f,df)
  elseif(val)then
   call sph_bessel(r,n,r_c,f)
  endif
 case(fin_den_f_indx_cheb) ! chebyshev
  if(sd)then
   call chebyshev(r,n,r_c,f,df,d2f)
  elseif(fd)then
   call chebyshev(r,n,r_c,f,df)
  elseif(val)then
   call chebyshev(r,n,r_c,f)
  endif
 case(fin_den_f_indx_exp) ! exponential
  if(sd)then
   call exponential(r,n,k_0,exp_step,f,df,d2f)
  elseif(fd)then
   call exponential(r,n,k_0,exp_step,f,df)
  elseif(val)then
   call exponential(r,n,k_0,exp_step,f)
  endif
 case default
  call errstop('CHG_DEN_BASIS_FNC','Basis type not recognised.')
 end select

 END SUBROUTINE fin_den_basis_fnc


 REAL(dp) FUNCTION rho_tail(r,L,alpha,beta,rho_L)
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: r,L,alpha,beta,rho_L
 rho_tail=rho_L*(r/L)**(2*beta)*exp(-2.d0*alpha*(r-L))
 END FUNCTION rho_tail


 SUBROUTINE setup_fourier_basis
!----------------------------------------------------------!
! Perform required setup for calculation of Fourier basis. !
!----------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,n

 n=gset(iplot)

 allocate(expval_pgmap(3,expval_ngvec(n)),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_FOURIER_BASIS','Allocation problem : &
  &expval_pgmap/grange arrays.')
 allocate(expval_expigdotr(expval_ngvec(n)),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_FOURIER_BASIS','Allocation error: &
  &expval_expigdotr array.')

! Map G vectors onto repeats of primitive lattice.
 expval_grange=0

 vec1(1:3)=expval_amat(1,1:3,n)
 vec2(1:3)=expval_amat(2,1:3,n)
 vec3(1:3)=expval_amat(3,1:3,n)

 do i=1,expval_ngvec(n)

  n1=nint(one_over_twopi*(expval_gvec(1,i,n)*vec1(1)+&
   &expval_gvec(2,i,n)*vec1(2)+expval_gvec(3,i,n)*vec1(3)))
  if(abs(n1)>expval_grange)expval_grange=abs(n1)
  expval_pgmap(1,i)=n1

  n2=nint(one_over_twopi*(expval_gvec(1,i,n)*vec2(1)+&
   &expval_gvec(2,i,n)*vec2(2)+expval_gvec(3,i,n)*vec2(3)))
  if(abs(n2)>expval_grange)expval_grange=abs(n2)
  expval_pgmap(2,i)=n2

  n3=nint(one_over_twopi*(expval_gvec(1,i,n)*vec3(1)+&
   &expval_gvec(2,i,n)*vec3(2)+expval_gvec(3,i,n)*vec3(3)))
  if(abs(n3)>expval_grange)expval_grange=abs(n3)
  expval_pgmap(3,i)=n3

 enddo ! G vectors in set

 if(expval_grange<1)then
  call errstop('SETUP_FOURIER_BASIS','G vector sets in expval.data must &
   &include non-zero G.')
 endif

 allocate(expval_mwork(3,-expval_grange:expval_grange),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_FOURIER_BASIS','Allocation problem : &
  &expval_mwork')

 END SUBROUTINE setup_fourier_basis


 SUBROUTINE compute_fourier_basis(rvec)
!--------------------------------------------------------------------!
! Compute the required exp(-iG.r) terms used in the reciprocal space !
! representations of expectation values.                             !
!--------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,n
 REAL(dp),INTENT(in) :: rvec(3)
 REAL(dp) p(3),cos1,cos2,cos3,sin1,sin2,sin3

 n=gset(iplot)

 vec1(1:3)=expval_bmat(1,1:3,n)
 vec2(1:3)=expval_bmat(2,1:3,n)
 vec3(1:3)=expval_bmat(3,1:3,n)

 p(1)=-(vec1(1)*rvec(1)+vec1(2)*rvec(2)+vec1(3)*rvec(3))
 p(2)=-(vec2(1)*rvec(1)+vec2(2)*rvec(2)+vec2(3)*rvec(3))
 p(3)=-(vec3(1)*rvec(1)+vec3(2)*rvec(2)+vec3(3)*rvec(3))

 cos1=cos(p(1)) ; cos2=cos(p(2)) ; cos3=cos(p(3))
 sin1=sin(p(1)) ; sin2=sin(p(2)) ; sin3=sin(p(3))

 expval_mwork(1:3,0)=cmplx(1.d0,0.d0,dp)
! Must have G/=0 in basis or will crash here.
 expval_mwork(1,1)=cmplx(cos1,sin1,kind=dp)
 expval_mwork(2,1)=cmplx(cos2,sin2,kind=dp)
 expval_mwork(3,1)=cmplx(cos3,sin3,kind=dp)
 expval_mwork(1,-1)=cmplx(cos1,-sin1,kind=dp)
 expval_mwork(2,-1)=cmplx(cos2,-sin2,kind=dp)
 expval_mwork(3,-1)=cmplx(cos3,-sin3,kind=dp)

 do i=2,expval_grange
  expval_mwork(1,i)=expval_mwork(1,1)*expval_mwork(1,i-1)
  expval_mwork(2,i)=expval_mwork(2,1)*expval_mwork(2,i-1)
  expval_mwork(3,i)=expval_mwork(3,1)*expval_mwork(3,i-1)
  expval_mwork(1,-i)=conjg(expval_mwork(1,i))
  expval_mwork(2,-i)=conjg(expval_mwork(2,i))
  expval_mwork(3,-i)=conjg(expval_mwork(3,i))
 enddo

! Insert alternative faster version for real arrays - xxxx.
 do i=1,expval_ngvec(gset(iplot))
  n1=expval_pgmap(1,i); n2=expval_pgmap(2,i); n3=expval_pgmap(3,i)
  expval_expigdotr(i)=expval_mwork(1,n1)*expval_mwork(2,n2)*expval_mwork(3,n3)
 enddo

 END SUBROUTINE compute_fourier_basis


 SUBROUTINE sph_bessel(r,n,r_c,bessel,dbessel,d2bessel)
!--------------------------------------------------------------------!
! This function returns the spherical bessel function of order n and !
! cutoff r_c evaluated at r. It is normalised between 0 and r_c      !
!                                                                    !
! phi_n(r) = 1/sqrt(2pi*r_c)(sin(n*pi*r/r_c)/r)                      !
! int(phi_n(r)*phi_m(r)d3r)=del_n,m                                  !
!--------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 REAL(dp),INTENT(in) :: r,r_c
 REAL(dp),INTENT(inout) :: bessel(n)
 REAL(dp),INTENT(inout),OPTIONAL :: dbessel(n),d2bessel(n)
 INTEGER i
 REAL(dp) prefactor,pi_over_rc,alpha,arg,cos_arg,sin_arg,inv_r,inv_r2,inv_r3
 LOGICAL fd,sd

 fd=present(dbessel)
 sd=fd.and.present(d2bessel)

 prefactor=1.d0/sqrt(twopi*r_c)
 pi_over_rc=pi/r_c
 alpha=0.d0

 if(r/=0.d0)then
  inv_r=1.d0/r
  if(fd)inv_r2=inv_r*inv_r
  if(sd)inv_r3=inv_r2*inv_r
 endif

 do i=1,n
  alpha=alpha+pi_over_rc
  arg=alpha*r
  sin_arg=sin(arg)
  if(r==0.d0)then
   bessel(i)=prefactor*alpha
   if(fd)dbessel(i)=0.d0
   if(sd)d2bessel(i)=-2.d0*bessel(i)*alpha*alpha
  else
   bessel(i)=prefactor*inv_r*sin_arg
   if(fd)then
    cos_arg=cos(arg)
    dbessel(i)=prefactor*inv_r2*(arg*cos_arg-sin_arg)
    if(sd)d2bessel(i)=prefactor*inv_r3*((2.d0-arg*arg)*sin_arg-2.d0*arg*cos_arg)
   endif
  endif
 enddo

 END SUBROUTINE sph_bessel


 SUBROUTINE chebyshev(r,n,r_c,cheb,dcheb,d2cheb)
!-----------------------------------------------------------------------!
! This function returns the chebyshev functions of the first kind up to !
! order n between 0 and r_c.                                            !
!-----------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 REAL(dp),INTENT(in) :: r,r_c
 REAL(dp),INTENT(inout) :: cheb(n)
 REAL(dp),INTENT(inout),OPTIONAL :: dcheb(n),d2cheb(n)
 INTEGER i
 REAL(dp) x
 LOGICAL fd,sd

 if(n==0)return

 x=2.d0*r/r_c-1.d0
 fd=present(dcheb)
 sd=fd.and.present(d2cheb)

! Initial values
 cheb(1)=1.d0
 if(fd)dcheb(1)=0.d0
 if(sd)d2cheb(1)=0.d0
 if(n>1)then
  cheb(2)=x
  if(fd)dcheb(2)=1.d0
  if(sd)d2cheb(2)=0.d0
 endif

! Recurrence relation
 do i=3,n
  cheb(i)=2*x*cheb(i-1)-cheb(i-2)
  if(fd)dcheb(i)=2.d0*(x*dcheb(i-1)+cheb(i-1))-dcheb(i-2)
  if(sd)d2cheb(i)=2.d0*(x*d2cheb(i-1)+2.d0*dcheb(i-1))-d2cheb(i-2)
 enddo

! Chain rule
 do i=1,n
  if(fd)dcheb(i)=2.d0/r_c*dcheb(i)
  if(sd)d2cheb(i)=4.d0/r_c**2*d2cheb(i)
 enddo

 END SUBROUTINE chebyshev


 SUBROUTINE exponential(r,n,k_0,exp_step,expn,dexpn,d2expn)
!------------------------------------------------------!
! This function returns the n exponentials of the form !
! exp(-k0*exp_step^(p-p_o)r).                          !
!------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 REAL(dp),INTENT(in) :: r,k_0,exp_step
 REAL(dp),INTENT(inout) :: expn(n)
 REAL(dp),INTENT(inout),OPTIONAL :: dexpn(n),d2expn(n)
 REAL(dp) power
 INTEGER i,p_min
 LOGICAL fd,sd

 fd=present(dexpn)
 sd=fd.and.present(d2expn)
 p_min=-((n-1)/2)

 power=-k_0*exp_step**p_min
 do i=1,n
  expn(i)=exp(power*r)
  if(fd)dexpn(i)=power*expn(i)
  if(sd)d2expn(i)=power*dexpn(i)
  power=exp_step*power
 enddo

 END SUBROUTINE exponential


 REAL(dp) FUNCTION gamma_ln(x_plus_one)
!--------------------------------------------------!
! This function returns the value of ln(gamma(x)). !
!--------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: x_plus_one
 INTEGER i
 REAL(dp),PARAMETER :: cof(6)=(/76.18009173d0,-86.50532033d0,24.0109822d0,&
  &-1.231739516d0,0.120858003d-2,-0.536382d-5/),factor=2.50662827465d0
 REAL(dp) x,zgh,tot
 x=x_plus_one-1.d0
 zgh=x+5.5d0
 zgh=(x+0.5d0)*log(zgh)-zgh
 tot=1.d0
 do i=1,6
  x=x+1.0d0
  tot=tot+cof(i)/x
 enddo
 gamma_ln=zgh+log(factor*tot)
 END FUNCTION gamma_ln


 REAL(dp) FUNCTION gamma_inc_upper(a,x)
!-------------------------------------------------------------------------!
! This function returns the value of the upper incomplete gamma function. !
!-------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: a,x
 REAL(dp) g_ln,gam_ser,gam_cf
 if(x<0.d0.or.a<=0.d0)call errstop('GAMMA_INC_UPPER','Out of bound arguments.')
 if(x<(a+1.d0))then
  call gamma_ser(a,x,g_ln,gam_ser)
  gamma_inc_upper=exp(g_ln)*(1.d0-gam_ser)
 else
  call gamma_cf(a,x,g_ln,gam_cf)
  gamma_inc_upper=exp(g_ln)*gam_cf
 endif
 END FUNCTION gamma_inc_upper


 SUBROUTINE gamma_ser(a,x,g_ln,gam_ser)
!------------------------------------------------------------!
! This function returns the incomplete gamma function Q(a,x) !
! evaluated by its continued fractions representation.       !
!------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: a,x
 REAL(dp),INTENT(out) :: g_ln,gam_ser
 INTEGER i
 INTEGER,PARAMETER :: maxiter=100
 REAL(dp) b,tot,del
 REAL(dp),PARAMETER :: eps=3.d-12
 g_ln=gamma_ln(a)
 if(x<=0.d0)then
  if(x<0.d0)call errstop('GAMMA_SER','Out of bound arguments.')
  gam_ser=0.d0
  return
 endif
 b=a
 tot=1.d0/a
 del=tot
 do i=1,maxiter
  b=b+1.d0
  del=del*x/b
  tot=tot+del
  if(abs(del)<(abs(tot)*eps))then
   gam_ser=tot*exp(-x+a*log(x)-g_ln)
   return
  endif
 enddo
 call wordwrap('GAMMA_SER: Not converged. a is too large or maxiter is too &
  &small.')
 END SUBROUTINE gamma_ser


 SUBROUTINE gamma_cf(a,x,g_ln,gam_cf)
!------------------------------------------------------------!
! This function returns the incomplete gamma function Q(a,x) !
! evaluated by its continued fractions representation.       !
!------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: a,x
 REAL(dp),INTENT(out) :: g_ln,gam_cf
 INTEGER i
 INTEGER,PARAMETER :: maxiter=100
 REAL(dp) g0,a0,a1,b0,b1,factor,ad,d,adf,g
 REAL(dp),PARAMETER :: eps=3.d-12
 g_ln=gamma_ln(a)
 g0=0.d0
 a0=1.d0
 a1=x
 b0=0.d0
 b1=1.d0
 factor=1.d0
 do i=1,maxiter
  ad=dble(i)
  d=ad-a
  a0=(a1+a0*d)*factor
  b0=(b1+b0*d)*factor
  adf=ad*factor
  a1=x*a0+adf*a1
  b1=x*b0+adf*b1
  if(a1/=0.d0)then
   factor=1.d0/a1
   g=b1*factor
   if(abs((g-g0)/g)<eps)then
    gam_cf=exp(-x+a*log(x)-g_ln)*g
    return
   else
    g0=g
   endif
  endif
 enddo
 call wordwrap('GAMMA_CF: Not converged. a is too large or maxiter is too &
  &small.')
 END SUBROUTINE gamma_cf


 SUBROUTINE eval_cheb_overlap(p_max,n_min,n_max,L,H)
!-------------------------------------------------------------------!
! Evaluates the overlap integral of spherically symmetric Chebyshev !
! polynomials with respect to a weight r^n.                         !
! H_p,q,n = 4*pi*[integral(r^(n+2)*h_p(r)*h_q(r)dr) from O to L]    !
!-------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: p_max,n_min,n_max
 REAL(dp),INTENT(in) :: L
 REAL(dp),INTENT(out) :: H(p_max,p_max,n_min:n_max)
 INTEGER n,p,q,ialloc,n_aux
 REAL(dp) fourpi,L_n3
 REAL(dp),ALLOCATABLE :: H_aux(:,:,:)

! Initialize.
 fourpi=2.d0*twopi
 n_aux=n_max+2*p_max
 allocate(H_aux(p_max,p_max,n_min:n_aux),stat=ialloc)
 if(ialloc/=0)call errstop('EVAL_CHEB_OVERLAP','Allocation error (H_aux).')
 H_aux=0.d0

! Compute (1:2,1:2,:) elements.
 do n=n_aux,n_min,-1
  H_aux(1,1,n)=fourpi/dble(n+3)
  H_aux(2,1,n)=fourpi*dble(n+2)/dble((n+3)*(n+4))
  H_aux(1,2,n)=H_aux(2,1,n)
  H_aux(2,2,n)=fourpi*dble(n*(n+5)+8)/dble((n+3)*(n+4)*(n+5))
 enddo

! Compute (3:,1,:) elements.
 do n=n_aux-1,n_min,-1
  do p=3,p_max
   H_aux(p,1,n)=4.d0*H_aux(p-1,1,n+1)-2.d0*H_aux(p-1,1,n)-H_aux(p-2,1,n)
   H_aux(1,p,n)=H_aux(p,1,n)
  enddo ! p
 enddo ! n

! Compute (3:,2,:) elements.
 do n=n_aux-2,n_min,-1
  do p=3,p_max
   H_aux(p,2,n)=8.d0*H_aux(p-1,1,n+2)-8.d0*H_aux(p-1,1,n+1)+&
    &2.d0*H_aux(p-1,1,n)-2.d0*H_aux(p-2,1,n+1)+H_aux(p-2,1,n)
   H_aux(2,p,n)=H_aux(p,2,n)
  enddo ! p
 enddo ! n

! Compute rest of matrix.
 do n=n_aux-3,n_min,-1
  do q=3,p_max
   do p=q,p_max
    H_aux(p,q,n)=16.d0*H_aux(p-1,q-1,n+2)-&
     &16.d0*H_aux(p-1,q-1,n+1)+4.d0*H_aux(p-1,q-1,n)-&
     &4.d0*H_aux(p-1,q-2,n+1)+2.d0*H_aux(p-1,q-2,n)-&
     &4.d0*H_aux(p-2,q-1,n+1)+2.d0*H_aux(p-2,q-1,n)+&
     &H_aux(p-2,q-2,n)
    H_aux(q,p,n)=H_aux(p,q,n)
   enddo
  enddo
 enddo

! Copy to output matrix.
 do n=n_min,n_max
  L_n3=L**(n+3)
  H(1:p_max,1:p_max,n)=H_aux(1:p_max,1:p_max,n)*L_n3
 enddo ! n

 deallocate(H_aux)

 END SUBROUTINE eval_cheb_overlap


 SUBROUTINE eval_exp_overlap(norder,k_0,exp_step,L,H)
!-----------------------------------------------------------!
! Evaluates the overlap integral of spherically symmetric   !
! exponentials with respect to a weight r^n.                !
! H_p,q = 4*pi*[integral(r^2*h_p(r)*h_q(r)dr) from O to L], !
! where h_p(r)=exp(-k_0*2^(p-p_o)*r).                       !
!-----------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: norder
 REAL(dp),INTENT(in) :: L,k_0,exp_step
 REAL(dp),INTENT(out) :: H(norder,norder)
 INTEGER p_min,p,q
 REAL(dp) fourpi,k_p,k_q,pref,k_pq,t1

 fourpi=2.d0*twopi
 p_min=-((norder-1)/2)
 k_p=k_0*exp_step**p_min

 do p=1,norder
  k_q=k_p
  do q=p,norder
   k_pq=k_p+k_q
   pref=-fourpi/(k_pq**3)
   t1=exp(-k_pq*L)*(k_pq*L*(k_pq*L+2.d0)+2.d0)-2.d0
   H(p,q)=pref*t1
   if(q/=p)H(q,p)=H(p,q)
   k_q=k_q*exp_step
  enddo ! q
  k_p=k_p*exp_step
 enddo ! p

 END SUBROUTINE eval_exp_overlap


 SUBROUTINE eval_exp_sph_ave(norder,n_min,n_max,k_0,exp_step,L,H)
!----------------------------------------------------------!
! Evaluates the spherical average of spherically symmetric !
! exponentials with respect to a weight r^n.               !
! H_p,n = 4*pi*[integral(r^(n+2)*h_p(r)dr) from O to L],   !
! where h_p(r)=exp(-k_0*2^(p-p_o)*r).                      !
!----------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: norder,n_min,n_max
 REAL(dp),INTENT(in) :: L,k_0,exp_step
 REAL(dp),INTENT(out) :: H(norder,n_min:n_max)
 INTEGER p_min,n,m,p,k
 REAL(dp) fourpi,k_p,k_pL,pref,num,den,L_mk,t1,t2,t3,t4

 fourpi=2.d0*twopi
 p_min=-((norder-1)/2)
 k_p=k_0*exp_step**p_min

 do p=1,norder

  if(k_p*L>=40.d0)then ! large p
   pref=-fourpi/k_p
   t1=L**(n_min+2)
   H(p,-2)=pref*(exp(-k_p*L)-1.d0)
   do n=n_min+1,n_max
    m=n+2
    t2=0.d0
    L_mk=L**(m-1)
    num=dble(m)
    den=k_p
    do k=1,m
     t2=t2+(num*L_mk/den)
     if(k==m)then
      t3=num
     else
      num=num*dble(m-k)
      den=den*k_p
      L_mk=L_mk/L
     endif
    enddo ! k
    t3=t3/(k_p**m)
    H(p,n)=pref*(exp(-k_p*L)*(t1+t2)-t3)
    t1=t1*L
   enddo ! n

  else ! small p - exponential "remainder" method
   k_pL=k_p*L
   pref=fourpi*exp(-k_pL)
   t1=k_p**2
   H(p,-2)=-fourpi/k_p*(exp(-k_p*L)-1.d0)
   do n=n_min+1,n_max
    t3=0.d0
    t4=k_pL**(n+2)
    k=n+2
    do
     k=k+1
     t4=t4*k_pL/dble(k)
     if(abs(t4)<=abs(t3)*epsilon(1.d0))exit
     t3=t3+t4
    enddo ! k
    H(p,n)=pref*(t3/t1)
    t1=t1*k_p
   enddo ! n
  endif
  k_p=k_p*exp_step
 enddo ! p

 END SUBROUTINE eval_exp_sph_ave


 SUBROUTINE plot_mol_density
!----------------------------------------------------!
! Plot the MOLECULAR DENSITY block from expval.data. !
!----------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k
 REAL(dp) rvec(3),bin_delta(3)

! Compute bin sizes.
 bin_delta=(mol_den_maxcoord-mol_den_mincoord)/dble(mol_den_norder)

 select case(plot_dim)
 case(1) ! lineplot
  do k=1,mol_den_norder(3)
   do j=1,mol_den_norder(2)
    do i=1,mol_den_norder(1)
     rvec=((/dble(i)-0.5d0,dble(j)-0.5d0,dble(k)-0.5d0/))*bin_delta
     r=sqrt(sum(rvec**2))
     write(10,*)r,expval_mol_den(i,j,k,iset)
    enddo ! i
   enddo ! j
  enddo ! k

 case(2) ! 2D_plot - insert gnuplot blank lines where appropriate
  if(mol_den_norder(3)==1)then ! xy plane
   k=1
   do j=1,mol_den_norder(2)
    do i=1,mol_den_norder(1)
     rvec=mol_den_mincoord+((/dble(i)-0.5d0,dble(j)-0.5d0,dble(k)-0.5d0/))*&
      &bin_delta
     write(10,*)rvec,expval_mol_den(i,j,k,iset)
    enddo ! i
    write(10,*)
   enddo ! j
  elseif(mol_den_norder(2)==1)then ! xz plane
   j=1
   do k=1,mol_den_norder(3)
    do i=1,mol_den_norder(1)
     rvec=mol_den_mincoord+((/dble(i)-0.5d0,dble(j)-0.5d0,dble(k)-0.5d0/))*&
      &bin_delta
     write(10,*)rvec,expval_mol_den(i,j,k,iset)
    enddo ! i
    write(10,*)
   enddo ! k
  else ! yz plane
   i=1
   do k=1,mol_den_norder(3)
    do j=1,mol_den_norder(2)
     rvec=mol_den_mincoord+((/dble(i)-0.5d0,dble(j)-0.5d0,dble(k)-0.5d0/))*&
      &bin_delta
     write(10,*)rvec,expval_mol_den(i,j,k,iset)
    enddo ! j
    write(10,*)
   enddo ! k
  endif ! which plane

 case(3) ! 3D_plot
  do k=1,mol_den_norder(3)
   do j=1,mol_den_norder(2)
    do i=1,mol_den_norder(1)
     rvec=mol_den_mincoord+((/dble(i)-0.5d0,dble(j)-0.5d0,dble(k)-0.5d0/))*&
      &bin_delta
     write(10,*)rvec,expval_mol_den(i,j,k,iset)
    enddo ! i
   enddo ! j
  enddo ! k

 end select

 END SUBROUTINE plot_mol_density


 SUBROUTINE plot_twop_dm_mom
!---------------------------------------------------------------------------!
! Plot the TWO-BODY DENSITY MATRIX (MOMENTUM SPACE) block from expval.data. !
!---------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,n,nsame
 REAL(dp) x,x2,w,k,dx,inv_nsame,t1
 REAL(dp),PARAMETER :: tol_same=1.d-12

! Calculate rescaling factor so lim r->inf is condensate fraction.
 t1=1.d0
 if(twop_dm_mom_rescale)then
  n=0
  do i=1,twop_dm_mom_nptypes_in_set(iset)
   n=n+min(nele_expval(twop_dm_mom_ptype_in_set(1,i,iset)),&
    &nele_expval(twop_dm_mom_ptype_in_set(2,i,iset)))
  enddo ! i
  select case(expval_periodicity)
  case(1) ; t1=(simcell_length**2)/real(n,dp)
  case(2) ; t1=(simcell_area**2)/real(n,dp)
  case(3) ; t1=(simcell_volume**2)/real(n,dp)
  end select
 endif ! twop_dm_mom_rescale

 select case(plot_dim)
 case(1) ! line
! Sort k points by size (only necessary if k_offset/=0).
  do i=1,expval_ngvec(1)
   k=expval_twop_dm_mom_k(i,iset)
   do j=i+1,expval_ngvec(1)
    if(k>expval_twop_dm_mom_k(j,iset))then
     x=expval_twop_dm_mom(j,iset) ; x2=expval_twop_dm_mom2(j,iset)
     k=expval_twop_dm_mom_k(j,iset)
     expval_twop_dm_mom(j,iset)=expval_twop_dm_mom(i,iset)
     expval_twop_dm_mom2(j,iset)=expval_twop_dm_mom2(i,iset)
     expval_twop_dm_mom_k(j,iset)=expval_twop_dm_mom_k(i,iset)
     expval_twop_dm_mom(i,iset)=x
     expval_twop_dm_mom2(i,iset)=x2
     expval_twop_dm_mom_k(i,iset)=k
    endif
   enddo
  enddo
  x=0.d0 ; x2=0.d0 ; w=0.d0 ; nsame=0
  do i=1,expval_ngvec(1)
   x=x+expval_twop_dm_mom(i,iset) ; x2=x2+expval_twop_dm_mom2(i,iset)
   w=w+expval_twop_dm_mom_weight(iset) ; nsame=nsame+1
   if(i<expval_ngvec(1))then
    if(abs(expval_twop_dm_mom_k(i+1,iset)-expval_twop_dm_mom_k(i,iset))<&
     &tol_same)cycle
   endif
   inv_nsame=1.d0/dble(nsame) ; x=x*inv_nsame ; x2=x2*inv_nsame ; w=w*inv_nsame
   dx=0.d0 ; if(w>1.d0)dx=sqrt(abs(x2-x*x)/(w-1.d0))
   if(twop_dm_mom_rescale)then
    x=t1*x
    dx=t1*dx
   endif
   write(10,'(e18.10,1x,e18.10,1x,e18.10)')&
    &expval_twop_dm_mom_k(i,iset),x,dx
   x=0.d0 ; x2=0.d0 ; w=0.d0 ; nsame=0
  enddo
 case default
  call errstop('PLOT_TWOP_DM_MOM','Momentum density plots must be &
   &one-dimensional.')
 end select

 END SUBROUTINE plot_twop_dm_mom


 SUBROUTINE plot_cond_fraction_mom
!------------------------------------------------------!
! Plot the CONDENSATE FRACTION block from expval.data. !
!------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,nsame,n
 REAL(dp) x,x2,w,k,dx,t1,inv_nsame
 REAL(dp),PARAMETER :: tol_same=1.d-12

  n=0
  do i=1,cond_frac_mom_nptypes_in_set(iset)
   n=n+min(nele_expval(cond_frac_mom_ptype_in_set(1,i,iset)),&
    &nele_expval(cond_frac_mom_ptype_in_set(2,i,iset)))
  enddo ! i
  select case(expval_periodicity)
  case(1) ; t1=(simcell_length**2)/real(n,dp)
  case(2) ; t1=(simcell_area**2)/real(n,dp)
  case(3) ; t1=(simcell_volume**2)/real(n,dp)
  end select

 select case(plot_dim)
 case(1) ! line
  ! Sort k points by size (only necessary if k_offset/=0).
  do i=1,expval_ngvec(1)
   k=expval_cond_frac_mom_k(i,iset)
   do j=i+1,expval_ngvec(1)
    if(k>expval_cond_frac_mom_k(j,iset))then
     x=expval_cond_frac_mom(j,iset) ; x2=expval_cond_frac_mom2(j,iset)
     k=expval_cond_frac_mom_k(j,iset)
     expval_cond_frac_mom(j,iset)=expval_cond_frac_mom(i,iset)
     expval_cond_frac_mom2(j,iset)=expval_cond_frac_mom2(i,iset)
     expval_cond_frac_mom_k(j,iset)=expval_cond_frac_mom_k(i,iset)
     expval_cond_frac_mom(i,iset)=x
     expval_cond_frac_mom2(i,iset)=x2
     expval_cond_frac_mom_k(i,iset)=k
    endif
   enddo
  enddo
  x=0.d0 ; x2=0.d0 ; w=0.d0 ; nsame=0
  do i=1,expval_ngvec(1)
   x=x+expval_cond_frac_mom(i,iset) ; x2=x2+expval_cond_frac_mom2(i,iset)
   w=w+expval_cond_frac_mom_weight(iset) ; nsame=nsame+1
   if(i<expval_ngvec(1))then
    if(abs(expval_cond_frac_mom_k(i+1,iset)-expval_cond_frac_mom_k(i,iset))<&
     &tol_same)cycle
   endif
   inv_nsame=1.d0/dble(nsame) ; x=x*inv_nsame ; x2=x2*inv_nsame ; w=w*inv_nsame
   dx=0.d0 ; if(w>1.d0)dx=sqrt(abs(x2-x*x)/(w-1.d0))
   write(10,'(e18.10,1x,e18.10,1x,e18.10)')&
    &expval_cond_frac_mom_k(i,iset),x*t1,dx*t1
   x=0.d0 ; x2=0.d0 ; w=0.d0 ; nsame=0
  enddo

 case default
  call errstop('PLOT_COND_FRACTION_MOM','Condensate fraction &
   &(momentum space) estimator plots must be one-dimensional.')
 end select

 END SUBROUTINE plot_cond_fraction_mom


 SUBROUTINE inv_33(A,B,det)
!------------------------------------------------------------------------!
! This subroutine calculates the inverse B of matrix A.  Optionally, the !
! determinant of A is returned.  A and B are real, 3x3 matrices.         !
!------------------------------------------------------------------------!
  IMPLICIT NONE
  REAL(dp),INTENT(in) :: A(3,3)
  REAL(dp),INTENT(out) :: B(3,3)
  REAL(dp),INTENT(out),OPTIONAL :: det
  REAL(dp) d
  d=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(2,1)*(A(3,2)*A(1,3)-A(1,2)*A(3,3))+ &
   &A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
  if(d==0.d0)call errstop('INV_33','Singular matrix')
  if(present(det))det=d
  d=1.d0/d
  B(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))*d
  B(1,2)=(A(3,2)*A(1,3)-A(1,2)*A(3,3))*d
  B(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))*d
  B(2,1)=(A(3,1)*A(2,3)-A(2,1)*A(3,3))*d
  B(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))*d
  B(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))*d
  B(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))*d
  B(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))*d
  B(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))*d
 END SUBROUTINE inv_33


END PROGRAM plot_expval

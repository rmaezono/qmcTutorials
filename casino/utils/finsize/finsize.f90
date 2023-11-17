 MODULE finsize
!-----------------------------------------------------------------------------!
! FINSIZE                                                                     !
! Calculates post-calculation table of KE and XC finite-size corrections and  !
! finite-size-corrected Ewald and/or MPC energies.                            !
!                                                                             !
! Previous VMC or DMC calculation should have been run with the flag          !
! FINITE_SIZE_CORR activated in the input file.                               !
!                                                                             !
! This utility requires as input some structure factor data accumulated in    !
! the expval.data file, as well as the files 'input' and 'out'.               !
!                                                                             !
! Note that the KE correction is not recomputed from the current              !
! correlation.data file, but is taken from the value printed in the out file. !
!                                                                             !
! MDT 4.2007                                                                  !
!-----------------------------------------------------------------------------!
 USE dsp,          ONLY : dp
 USE esdf,         ONLY : esdf_init,esdf_defined,esdf_string,esdf_integer, &
  &esdf_boolean,esdf_warnout
 USE format_utils, ONLY : i2s
 USE run_control,  ONLY : errstop
 USE store,        ONLY : o
 IMPLICIT NONE

 INTEGER ialloc,ierr,netot
 INTEGER,PARAMETER :: current_file_version=1,max_header_lines=10,io_expval=10,&
  &io_out=11
 REAL(dp) timag,treal
 REAL(dp),PARAMETER :: pi=3.14159265358979324d0,twopi=2.d0*pi,&
  &one_over_twopi=1.d0/twopi,fourpi=2.d0*twopi,one_over_pi=1.d0/pi
 LOGICAL expval_in_expval,structure_factor_in_expval,file_present
 CHARACTER(80) :: char_80

! Basic info
 INTEGER expval_periodicity,npcells,expval_scell_matrix(3,3),&
  &expval_ngsets,file_version,num_particle_types,expval_dimensionality
 INTEGER,ALLOCATABLE :: nele_expval(:)
 REAL(dp) simcell_area,simcell_length,simcell_volume,pa1(3),pa2(3),pa3(3),&
  &pamat(3,3),pbmat(3,3),amat(3,3),bmat(3,3),wigner_seitz_radius,self_term
 CHARACTER(80) title

! G vector sets
 INTEGER,ALLOCATABLE :: expval_ngvec(:)
 REAL(dp),ALLOCATABLE :: e_cutoff(:),expval_gvec(:,:,:)

! Structure factor
 INTEGER sf_nsets,sf_gset,sfsden_nsets
 INTEGER,ALLOCATABLE :: sf_ptype(:,:),sfsden_ptype(:)
 REAL(dp) expval_sf_weight,expval_sfsden_weight
 REAL(dp),ALLOCATABLE :: expval_sf(:,:)
 COMPLEX(dp),ALLOCATABLE :: expval_sfsden(:,:)
 LOGICAL homogeneous_structure_factor
 CHARACTER(3) accumulation_method_sf

! Finite_size_corrections
 INTEGER nstar,xc_corr_method,ndata,nparam1,nparam2,nparam
 REAL(dp) xc_corr,ke_corr,chi_squared_sf,ewald_energy,mpc_energy, &
  &finite_size_const_c
 REAL(dp),ALLOCATABLE :: sf_sum(:),gs(:),g2s(:),a(:),weight(:),sf_fit(:),&
  &epoly(:),polyb(:)
 LOGICAL esupercell,electron_gas
 CHARACTER(20) interaction,atom_basis_type
 LOGICAL,PARAMETER :: USE_WEIGHTS=.false.

! Overflow protection
 REAL(dp),PARAMETER :: min_exp=1.d-150,max_exp=1.d150
 REAL(dp) minimum_exp_arg_actual,maximum_exp_arg_actual,&
  &minimum_exp_arg,maximum_exp_arg,largest_representable_number,&
  &smallest_representable_number

! ELIMINATE THIS
  INTEGER ndet


CONTAINS


  SUBROUTINE finsize_main
  IMPLICIT NONE

! ELIMINATE THIS
  ndet=1


!----------------------------------------------------------------------------

  write(6,*)
  write(6,*)'FINSIZE: calculation of finite-size corrections from structure &
   &factor'
  write(6,*)'===============================================================&
   &======'
  write(6,*)

  inquire(file='expval.data',exist=file_present)
  if(.not.file_present)then
   write(6,*)"No file 'expval.data' containing the structure factor data."
   write(6,*)
   call errstop('FINSIZE_MAIN','Quitting.')
  endif
  write(6,*)'Found expval.data file. Attempting to read data.'
  call read_expval
  write(6,*)'Successfully read structure factor data.'
  write(6,*)

  inquire(file='input',exist=file_present)
  if(.not.file_present)then
   write(6,*)"Need to read file 'input' but it is not present."
   write(6,*)
   call errstop('FINSIZE_MAIN','Quitting.')
  endif
  write(6,*)'Found input file. Attempting to read data.'
  call read_input
  write(6,*)'Done.'
  write(6,*)

  inquire(file='out',exist=file_present)
  if(.not.file_present)then
   write(6,*)"Need to read file 'out' but it is not present."
   write(6,*)
   call errstop('FINSIZE_MAIN','Quitting.')
  endif
  write(6,*)'Found out file. Attempting to read data.'
  call read_output
  write(6,*)'Electron self-image term (au)  : ',2.d0*self_term
  write(6,*)'KE finite size correction (au) : ',ke_corr
  write(6,*)
  write(6,*)'===============================================================&
   &======'
  write(6,*)

  call talk_to_user

  call eval_finite_size_corr

  stop

 END SUBROUTINE finsize_main


 SUBROUTINE read_expval
!-----------------------------------------------------------------------------!
! Call routine to determine what data blocks expval.data contains, then       !
! call the relevant additional routines to read these blocks.                 !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE

 open(unit=io_expval,file='expval.data',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_EXPVAL','Problem opening expval.data.')

! Find out which data sets are present in expval.data
 call find_expval_sets

! Detailed read of 'START EXPVAL' block
 call read_expval_block

! If there are G vector sets given, then read them
 if(expval_ngsets>0)then
  call read_g_vector_sets
 endif

! Detailed read of whatever sets are present.
 call read_structure_factor
 if(structure_factor_in_expval)then
  write(6,*)'File contains structure factor block.'
 else
  call errstop('READ_EXPVAL','No structure factor blocks were found in this &
   &file.')
 endif

 close(io_expval)

 END SUBROUTINE read_expval


 SUBROUTINE find_expval_sets
!---------------------------------------------------------------------!
! If expval.data exists do a cursory read to check the structure of   !
! the file is correct and that it contains a structure factor block.  !
!---------------------------------------------------------------------!
 IMPLICIT NONE

 structure_factor_in_expval=.false. ; expval_in_expval=.false.

 do

  read(io_expval,'(a)',iostat=ierr)char_80

  if(ierr>0)then
   call errstop('FIND_EXPVAL_SETS','Problem reading expval.data.')
  elseif(ierr<0)then
   exit
  endif

  if(trim(adjustl(char_80))=='START EXPVAL')then

   do
    read(io_expval,'(a)',iostat=ierr)char_80
    if(trim(adjustl(char_80))=='END EXPVAL')exit
    if(ierr>0)then
     call errstop('FIND_EXPVAL_SETS','Problem reading expval.data <2>.')
    elseif(ierr<0)then
     exit
    endif
    if(trim(adjustl(char_80))=='START STRUCTURE FACTOR')then
     if(structure_factor_in_expval)call errstop('FIND_EXPVAL_SETS', &
      &'More than one STRUCTURE FACTOR block in expval.data.')
     structure_factor_in_expval=.true.
    endif
   enddo ! headers in the expval block

   expval_in_expval=.true.

  endif ! char_80

 enddo

 rewind(io_expval)

 if(.not.expval_in_expval)then
  call errstop('FIND_EXPVAL_SETS','The expval.data file must contain a&
   &"START EXPVAL" block.')
 endif

 END SUBROUTINE find_expval_sets


 SUBROUTINE read_expval_block
!-----------------------------------------------------------------------------!
! Read basic info following START EXPVAL in expval.data.                      !
!-----------------------------------------------------------------------------!
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
 if(expval_dimensionality/=2.and.expval_dimensionality/=3) &
  &call errstop('READ_EXPVAL_BLOCK','Invalid dimensionality in expval.data. &
  &Must be 2 or 3.')
 read(io_expval,*,err=1,end=1)
 read(io_expval,*,err=1,end=1)expval_periodicity
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
   read(char_80,*,err=1,end=1)expval_scell_matrix(1,1),expval_scell_matrix(2,2),&
    &expval_scell_matrix(3,3)
  endif ! ierr/=0
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)simcell_volume
  read(io_expval,*,err=1,end=1)
  read(io_expval,*,err=1,end=1)wigner_seitz_radius
 case default
  call errstop('READ_EXPVAL_BLOCK','Invalid periodicity.')
 end select
 pamat(1,:)=pa1(:) ; pamat(2,:)=pa2(:) ; pamat(3,:)=pa3(:)
 call inv_33(pamat,pbmat)
 pbmat(:,:)=transpose(pbmat(:,:))*twopi
 amat=matmul(dble(expval_scell_matrix),pamat)
 call inv_33(amat,bmat)
 bmat=transpose(bmat)*twopi
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
!-----------------------------------------------------------------------------!
! Read in GVECTOR SET n blocks from expval.data.                              !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k
 REAL(dp) temp_bmat(3,3)
 REAL(dp),PARAMETER :: tol=1.d-8
 CHARACTER(1) set_g

 allocate(e_cutoff(expval_ngsets),expval_ngvec(expval_ngsets),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('READ_G_VECTOR_SETS','Allocation problem : &
  &e_cutoff,expval_ngvec')

! Skip read to get expval_gvec array size
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
   read(io_expval,*,err=1,end=1)temp_bmat(j,1:3)
  enddo ! j
  if(.not.all(abs(temp_bmat(1:expval_periodicity,1:expval_periodicity) &
   &-pbmat(1:expval_periodicity,1:expval_periodicity))<tol).and.&
   &.not.all(abs(temp_bmat(1:expval_periodicity,1:expval_periodicity)&
   &-bmat(1:expval_periodicity,1:expval_periodicity))<tol))call &
   &errstop('READ_G_VECTOR_SETS','Set of reciprocal lattice vectors in &
   &expval.data do not correspond to the primitive cell or supercell.')

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


 SUBROUTINE read_structure_factor
!-----------------------------------------------------------------------------!
! Read in a STRUCTURE FACTOR block from expval.data.                          !
!-----------------------------------------------------------------------------!
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
 accumulation_method_sf=trim(adjustl(char_80))
 if(accumulation_method_sf/='VMC'.and.accumulation_method_sf/='DMC') &
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

  if(trim(adjustl(char_80))/='START SET '//nset)then
   call errstop('READ_STRUCTURE_FACTOR','Error reading START SET x line &
    &for x = '//trim(i2s(set)))
  endif

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


 SUBROUTINE read_input
!-----------------------------------------------------------------------------!
! Read the input file for required keywords.                                  !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 LOGICAL check

 call esdf_init('input') ! i.e. read the input file
 call esdf_warnout

 check=esdf_defined('interaction','T')
 if(.not.check)then
  call errstop('READ_INPUT','INTERACTION must be defined in input file.')
 endif
 check=esdf_defined('atom_basis_type','T')
 if(.not.check)then
  call errstop('READ_INPUT','ATOM_BASIS_TYPE must be defined in input file.')
 endif
 check=esdf_defined('xc_corr_method','I')
 if(.not.check)then
  xc_corr_method=1
 else
  xc_corr_method=esdf_integer('xc_corr_method',1)
 endif

 interaction=esdf_string('interaction','default')
 atom_basis_type=esdf_string('atom_basis_type','default')
 esupercell=esdf_boolean('esupercell',.false.)

 if(trim(adjustl(atom_basis_type))=='none')then
  electron_gas=.true.
 else
  electron_gas=.false.
 endif

 if(trim(interaction)=='default')then
  call errstop('READ_INPUT','INTERACTION keyword may not take the value &
   &"default"')
 endif

 END SUBROUTINE read_input


 SUBROUTINE read_output
!-----------------------------------------------------------------------------!
! Read the output file to get the self-image term.                            !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i
 REAL(dp) tempr
 LOGICAL found_self_image,found_ke_corr,found_c
 CHARACTER(1) a

 found_self_image=.false. ;found_ke_corr=.false. ; found_c=.false.

 open(unit=io_out,file='out',status='old',iostat=ierr)
 if(ierr/=0)call errstop('READ_OUTPUT','Problem opening out file.')

 do

  read(io_out,'(a)',iostat=ierr)char_80
  if(ierr>0)then
   call errstop('READ_OUTPUT','Problem reading out file.')
  elseif(ierr<0)then
   exit
  endif

  if(index(char_80,'Electron self-image term v_M (au)')/=0)then
   read(char_80,*,err=1,end=1)(a,i=1,6),self_term
   self_term=0.5d0*self_term
   found_self_image=.true.
  endif

  if(electron_gas)then
   if(trim(adjustl(char_80(1:40)))=='Total KE correction (a.u. per particle)')&
    &then
    read(char_80,*,err=1,end=1)(a,i=1,7),ke_corr
    if(esupercell)ke_corr=ke_corr*real(netot,dp)
    found_ke_corr=.true.
   endif
  else
   if(trim(adjustl(char_80(1:42)))=='Total KE correction (a.u. per prim. &
    &cell)')then
    read(char_80,*,err=1,end=1)(a,i=1,8),ke_corr
    if(esupercell)ke_corr=ke_corr*real(npcells,dp)
    found_ke_corr=.true.
   endif
  endif ! electron_gas

  if(trim(adjustl(char_80))=='alpha                     C')then
   do i=1,5
    read(io_out,'(a)',iostat=ierr)char_80
    if(ierr/=0)call errstop('READ_OUTPUT','Error reading out.')
   enddo ! i
   read(io_out,*,err=13,end=13)tempr,finite_size_const_c
   found_c=.true.
  endif ! C

  if(found_self_image.and.found_ke_corr.and.found_c)exit

 enddo

 if(.not.found_self_image)call errstop('READ_OUTPUT','Reached end of out file &
  &without finding electron self-image term.')

 if(.not.found_ke_corr)then
  write(6,*)'Reached end of out file without finding total KE correction.'
  write(6,*)'Setting KE correction to zero.'
  write(6,*)
  ke_corr=0.d0
 endif ! found_ke_corr

 if(.not.found_c)then
  write(6,*)'Reached end of out file without finding finite-size corr. const.'
  write(6,*)'Setting finite-size correction constant C to zero.'
  write(6,*)
  finite_size_const_c=0.d0
 endif ! found_ke_corr

 close(io_out)

 return

1 call errstop('READ_OUTPUT','Error reading electron self-image line of out &
   &file.')

13 call errstop('READ_OUTPUT','Error reading finite-size constant C from out &
   &file.')

 END SUBROUTINE read_output


 SUBROUTINE talk_to_user
!-----------------------------------------------------------------------------!
! Discuss with the user precisely what he wants to do.                        !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE

! Ask required questions

 if(esupercell)then
  write(6,*)'ESUPERCELL set to T in input file.'
  write(6,*)'All total energies and corrections are thus PER SIMULATION CELL.'
  write(6,*)
 else
  if(electron_gas)then
   write(6,*)'All total energies and corrections are PER PARTICLE.'
  else
   write(6,*)'All total energies and corrections are PER PRIMITIVE CELL.'
  endif
  write(6,*)
 endif

 ewald_energy=1.d10 ; mpc_energy=1.d10
 select case(trim(interaction))
  case('ewald')
   write(6,*)'Please input the value of the total Ewald energy:'
   do while(ewald_energy==1.d10)
    read(5,*,iostat=ierr)ewald_energy
    if(ierr>0)ewald_energy=1.d10
   enddo
  case('mpc')
   write(6,*)'Please input the value of the total MPC energy:'
   do while(mpc_energy==1.d10)
    read(5,*,iostat=ierr)mpc_energy
    if(ierr>0)mpc_energy=1.d10
   enddo
  case('ewald_mpc')
   write(6,*)'Please input the value of the total Ewald energy:'
   do while(ewald_energy==1.d10)
    read(5,*,iostat=ierr)ewald_energy
    if(ierr>0)ewald_energy=1.d10
   enddo
   write(6,*)'Please input the value of the total MPC energy:'
   do while(mpc_energy==1.d10)
    read(5,*,iostat=ierr)mpc_energy
    if(ierr>0)mpc_energy=1.d10
   enddo
  case('mpc_ewald') ! mpc_ewald
   write(6,*)'Please input the value of the total MPC energy:'
   do while(mpc_energy==1.d10)
    read(5,*,iostat=ierr)mpc_energy
    if(ierr>0)mpc_energy=1.d10
   enddo
   write(6,*)'Please input the value of the total Ewald energy:'
   do while(ewald_energy==1.d10)
    read(5,*,iostat=ierr)ewald_energy
    if(ierr>0)ewald_energy=1.d10
   enddo
  case default
   call errstop('TALK_TO_USER','Error in select case for interaction keyword &
    &: bug.')
 end select

 END SUBROUTINE talk_to_user


 SUBROUTINE eval_finite_size_corr
!-----------------------------------------------------------------------------!
! Call the relevant routines to compute the finite size corrections, add      !
! them to the input total energies, and report.                               !
!-----------------------------------------------------------------------------!
 IMPLICIT NONE

 call overflow_protection

 call finite_size_corr_xc

 write(6,*)
 select case(xc_corr_method)
  case(1)
   write(6,'(t2,a)')'Finite size correction data:'
  case(2)
   write(6,'(t2,a)')'Finite size correction data (METHOD 2):'
  case default
   call errstop('EVAL_FINITE_SIZE_CORR','Unknown XC_CORR_METHOD')
 end select
 if(trim(interaction)/='mpc')then
  if(chi_squared_sf>0.d0)write(6,1)'Chi squared for structure factor fit', &
   &chi_squared_sf
  write(6,1)'Total XC correction (dXC)',xc_corr
 endif
 write(6,1)'Total KE correction (dKE)',ke_corr
 select case(trim(interaction))
  case ('ewald')
   write(6,1)'Total energy (Ewald)',ewald_energy
   write(6,1)'Total energy (Ewald) + dXC + dKE',ewald_energy+xc_corr+ke_corr
   write(6,1)'Total energy (Ewald) + dXC',ewald_energy+xc_corr
   write(6,1)'Total energy (Ewald) + dKE',ewald_energy+ke_corr
  case('mpc')
   write(6,1)'Total energy (MPC)',mpc_energy
   write(6,1)'Total energy (MPC) + dKE',mpc_energy+ke_corr
  case('ewald_mpc')
   write(6,1)'Total energy (Ewald)',ewald_energy
   write(6,1)'Total energy (Ewald) + dXC + dKE',ewald_energy+xc_corr+ke_corr
   write(6,1)'Total energy (Ewald) + dXC',ewald_energy+xc_corr
   write(6,1)'Total energy (Ewald) + dKE',ewald_energy+ke_corr
   write(6,1)'Total energy (MPC)',mpc_energy
   write(6,1)'Total energy (MPC) + dKE',mpc_energy+ke_corr
  case('mpc_ewald')
   write(6,1)'Total energy (MPC)',mpc_energy
   write(6,1)'Total energy (MPC) + dKE',mpc_energy+ke_corr
   write(6,1)'Total energy (Ewald)',ewald_energy
   write(6,1)'Total energy (Ewald) + dXC + dKE',ewald_energy+ke_corr+xc_corr
   write(6,1)'Total energy (Ewald) + dXC',ewald_energy+xc_corr
   write(6,1)'Total energy (Ewald) + dKE',ewald_energy+ke_corr
 end select
 write(6,*)

1 format(t2,a,t38,'(au) =',t48,f21.12)

 END SUBROUTINE eval_finite_size_corr


 SUBROUTINE finite_size_corr_xc
!--------------------------------------------------------------------!
! Compute and report XC finite size correction to the total energy.  !
!                                                                    !
! MDT 11.2006                                                        !
!--------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER g,i,j,k,n,ng,set
 INTEGER,ALLOCATABLE :: first_g_in_star(:)
 REAL(dp) g2,integral,t1,t2,term1,term2
! REAL(dp),EXTERNAL :: calcf_test
 LOGICAL,PARAMETER :: integrate_s=.true. ! if .false. integrate s-1
 INTERFACE
  FUNCTION calcf_xc_corr(xx)
   USE dsp,ONLY : dp
   REAL(dp),INTENT(in) :: xx
   REAL(dp) calcf_xc_corr
  END FUNCTION calcf_xc_corr
  FUNCTION calcf_xc_corr2(xx)
   USE dsp,ONLY : dp
   REAL(dp),INTENT(in) :: xx
   REAL(dp) calcf_xc_corr2
  END FUNCTION calcf_xc_corr2
  SUBROUTINE gauss_quadrature(x1,x2,n,func,integral)
   USE dsp,ONLY : dp
   INTEGER,INTENT(in) :: n
   REAL(dp),INTENT(in) :: x1,x2
   REAL(dp),INTENT(out) :: integral
   REAL(dp),EXTERNAL :: func
  END SUBROUTINE gauss_quadrature
 END INTERFACE

! Work out star structure of G vector set

 allocate(gs(expval_ngvec(sf_gset)),g2s(expval_ngvec(sf_gset)),&
  &first_g_in_star(expval_ngvec(sf_gset)+1),stat=ialloc)

 if(ialloc/=0)call errstop('FINITE_SIZE_CORR_XC','Allocation error : &
  &star information.')

 ng=1
 nstar=1
 gs(1)=0.d0 ; g2s(1)=0.d0
 first_g_in_star(1)=1
 do g=2,expval_ngvec(sf_gset)
  g2=expval_gvec(1,g,sf_gset)**2+expval_gvec(2,g,sf_gset)**2+&
   &expval_gvec(3,g,sf_gset)**2
  if(g2-g2s(nstar)>1.d-8)then ! we have discovered a new star
   nstar=nstar+1
   first_g_in_star(nstar)=g
   g2s(nstar)=g2
   gs(nstar)=sqrt(g2)
   ng=0
  endif
  ng=ng+1
 enddo ! G
 first_g_in_star(nstar+1)=first_g_in_star(nstar)+ng

 allocate(sf_sum(nstar),stat=ialloc)
 if(ialloc/=0)call errstop('FINITE_SIZE_CORR_XC','Allocation error : &
  &sf_sum array.')

! Spherically average the structure factor data over the G in each star.

 sf_sum(1:nstar)=0.d0

 if(sfsden_nsets==0)then ! homogeneous structure factor

  do set=1,sf_nsets

   j=sf_ptype(1,set) ; k=sf_ptype(2,set)
   if(j/=k)then
    t2=2.d0
   else
    t2=1.d0
   endif

   sf_sum(1)=sf_sum(1)+t2*(expval_sf(1,set)-nele_expval(j)*nele_expval(k))

   g=1
   do n=2,nstar
    ng=0
    do i=first_g_in_star(n),first_g_in_star(n+1)-1,2
     g=g+1
     sf_sum(n)=sf_sum(n)+t2*(expval_sf(g,set)+expval_sf(g+1,set))
     g=g+1
     ng=ng+2
    enddo ! G in star
   enddo ! stars

  enddo ! SF sets

 else ! inhomogeneous structure factor

  do set=1,sf_nsets

   j=sf_ptype(1,set) ; k=sf_ptype(2,set)
   if(j/=k)then
    t2=2.d0
   else
    t2=1.d0
   endif

   sf_sum(1)=sf_sum(1)+t2*(expval_sf(1,set)- &
    &real(expval_sfsden(1,j)*expval_sfsden(1,k),dp)) ! i.e. zero

   g=1
   do n=2,nstar
    ng=0
    do i=first_g_in_star(n),first_g_in_star(n+1)-1,2
     g=g+1
     sf_sum(n)=sf_sum(n)+t2*(expval_sf(g,set)+expval_sf(g+1,set)-&
      &2.d0*real(expval_sfsden(g,j)*expval_sfsden(g+1,k),dp))
     g=g+1
     ng=ng+2
    enddo ! G in star
   enddo ! stars

  enddo ! SF sets

 endif ! homogeneous density or not

 t2=real(npcells**2,dp)/real(netot,dp)
 do n=2,nstar
  ng=first_g_in_star(n+1)-first_g_in_star(n)
  if(ng/=0)then
   sf_sum(n)=sf_sum(n)*t2/real(ng,dp)
  else
   call errstop('FINITE_SIZE_CORR_XC','Zero G vectors in star - should &
    &not happen.')
  endif
 enddo

 select case(xc_corr_method)

  case(1) ! simple(!) method i.e. what Chiesa does.

   if(expval_periodicity==3)then
    t1=sf_sum(2)/gs(2)**2 ! Approximation to lim_k->0 S(k)/k^2.
    xc_corr=twopi*t1*real(netot,dp)/simcell_volume
   elseif(expval_periodicity==2)then ! 2D
! Use Eq. (59) of PRB 78, 125106 (2008).
    t1=sf_sum(2)/gs(2)**1.5d0 ! Approximation to lim_k->0 S(k)/k^(3/2).
    xc_corr=finite_size_const_c*dble(netot)*t1/simcell_area**1.25d0
   else
    call errstop('FINITE_SIZE_CORR_XC','Finite-size corrections only &
     &available for 3D- or 2D-periodic systems.')
   endif ! expval_periodicity

  case(2)

! Fit S = [1-exp(\sum_{n=2}^n1 a_n*G^n)] * (\sum_{m=2}^n2 b_m*G^m + 1)

   call fit_structure_factor

! Evaluate term 1

   if(integrate_s)then
    call gauss_quadrature(0.d0,gs(ndata),100,calcf_xc_corr,integral)
   else ! integrate s-1
    call gauss_quadrature(0.d0,gs(ndata),100,calcf_xc_corr2,integral)
   endif

   term1=integral*one_over_pi

! Evaluate term 2 [ twopi * e^2 / \Gamma \sum_{Gshells} n(G) * S(G) / G^2 ]

   if(integrate_s)then

    t2=0.d0
    do i=2,ndata ! shells of G
     ng=first_g_in_star(i+1)-first_g_in_star(i)
     t1=1.d0/g2s(i)
     t2=t2+real(ng,dp)*t1*sf_fit(i)
    enddo
    term2=t2*twopi/simcell_volume

   else ! integrate s-1

    t2=0.d0
    do i=2,ndata ! shells of G
     ng=first_g_in_star(i+1)-first_g_in_star(i)
     t1=1.d0/g2s(i)
     t2=t2+real(ng,dp)*t1*(sf_fit(i)-1.d0)
    enddo
    term2=t2*twopi/simcell_volume

   endif

! Work out final XC correction per supercell

   if(integrate_s)then
    xc_corr=term1-term2
   else
    xc_corr=term1-term2-self_term
   endif

  case default

   call errstop('FINITE_SIZE_CORR_XC','Unknown method type XC_CORR_METHOD')

 end select

 if(.not.esupercell)then
  if(electron_gas)then
   xc_corr=xc_corr/real(netot,dp)
  else
   xc_corr=xc_corr/real(npcells,dp)
  endif
 endif

 deallocate(first_g_in_star,g2s,sf_sum)

 END SUBROUTINE finite_size_corr_xc


 SUBROUTINE fit_structure_factor
!-----------------------------------------------------------------------------!
! Fit non-linear function to spherically-averaged structure factor as follows:!
!                                                                             !
! S(k) = [1 - exp (\sum_{n=2}^N a_n * G^n) ] * [ \sum_{m=2}^M b_m G^m ]       !
!                                                                             !
! MDT 6.2007                                                                  !
!-----------------------------------------------------------------------------!
 USE toms573,ONLY : nl2sol,dfault
 IMPLICIT NONE
 INTEGER i,j
 INTEGER,ALLOCATABLE :: iv(:)
 REAL(dp) poly,step,t1
 REAL(dp),ALLOCATABLE :: v(:),sf_fit2(:)
 LOGICAL,PARAMETER :: DEBUG=.false.
 EXTERNAL calcr_xc_corr,calcj_xc_corr

 nparam1=8 ; nparam2=6
 nparam=nparam1+nparam2
 ndata=nstar

 if(USE_WEIGHTS)then
! Construct weights (which seem to make things worse, so don't use by default..)
  allocate(weight(ndata),stat=ialloc)
  if(ialloc/=0)call errstop('FIT_STRUCTURE_FACTOR','Allocation problem : &
   &weights.')
  weight(1)=gs(2)
  do i=2,ndata-1
   weight(i)=abs(gs(i+1)-gs(i-1))
  enddo
  weight(ndata)=abs(gs(ndata)-gs(ndata-1))
  t1=maxval(weight(:))
  weight(:)=weight(:)/t1
 endif

! Allocate input/output/working arrays for NL2SOL.
 allocate(iv(60+nparam),v(93+ndata*(nparam+3)+nparam*(3*nparam+33)/2),&
  &stat=ialloc)
 if(ialloc/=0)call errstop('FIT_STRUCTURE_FACTOR','Allocation problem: &
  &NL2SOL arrays.')

! Set NL2SOL defaults for control parameters

 call dfault(iv,v)
 iv(14:15)=0 ; iv(19:24)=0 ; iv(21)=-1 ! suppress output
 iv(17)=500  ! Maximum number of function evaluations.
 iv(18)=100 ! Maximum number of iterations.

! Allocation
 allocate(a(nparam),epoly(ndata),polyb(ndata),sf_fit(ndata),stat=ialloc)
 if(ialloc/=0)call errstop('FIT_STRUCTURE_FACTOR','Allocation problem: &
  &work arrays.')

! Zero the parameter vector
 a(:)=0.d0

 write(o,*)
 write(o,*)'Evaluating XC finite size correction using XC_CORR_METHOD==2'

! Optimize
 call nl2sol(ndata,nparam,a,calcr_xc_corr,calcj_xc_corr,iv,v)

 if(iv(1)==3)then
  write(o,*)'Parameter convergence. Structure factor fit complete.'
 elseif(iv(1)==4)then
  write(o,*)'Relative-function convergence. Structure factor fit complete.'
 elseif(iv(1)==5)then
  write(o,*)'Parameter and relative-function convergence. Structure factor&
   & fit complete.'
 elseif(iv(1)==6)then
  write(o,*)'Absolute function convergence. Structure factor fit complete.'
 elseif(iv(1)==7)then
  write(o,*)'Singular convergence. Structure factor fit complete.'
 elseif(iv(1)==8)then
  write(o,*)'False convergence. Structure factor fit complete.'
 elseif(iv(1)==9.or.iv(1)==10)then
  write(o,*)'Maximum number of iterations reached. Structure factor fit &
   &complete.'
 else
  write(o,*)'NL2SOL return code: '//trim(i2s(iv(1)))
  write(o,*)'Structure factor fitting procedure was unsuccessful.'
 endif ! iv value.

! Evaluate fitted structure factor and hence chi-squared function

 chi_squared_sf=0.d0
 do j=1,ndata
  poly=0.d0
  do i=1,nparam1
   poly=poly+a(i)*gs(j)**(i+1)
  enddo
  sf_fit(j)=1.d0-exp(poly)
  poly=0.d0
  do i=1,nparam2
   poly=poly+a(nparam1+i)*gs(j)**(i+1)
  enddo
  sf_fit(j)=sf_fit(j)*(poly+1.d0)
  chi_squared_sf=chi_squared_sf+(sf_sum(j)-sf_fit(j))**2
 enddo

 if(DEBUG)then
  allocate(sf_fit2(200))
  step=gs(ndata)/200.d0
  sf_fit2=0.d0
  t1=0.d0
  do j=1,200
   poly=0.d0
   do i=1,nparam1
    poly=poly+a(i)*t1**(i+1)
   enddo
   sf_fit2(j)=1.d0-exp(poly)
   poly=0.d0
   do i=1,nparam2
    poly=poly+a(nparam1+i)*t1**(i+1)
   enddo
   sf_fit2(j)=sf_fit2(j)*(poly+1.d0)
   t1=t1+step
  enddo
  write(o,*)'DATA'
  do i=1,nstar
   write(o,*)gs(i),sf_sum(i)
  enddo
  write(o,*)'&'
  t1=0.d0
  do i=1,200
   write(o,*)t1,sf_fit2(i)
   t1=t1+step
  enddo
  deallocate(sf_fit2)
 endif ! DEBUG

 deallocate(iv,v,epoly,polyb)
 if(USE_WEIGHTS)deallocate(weight)

 END SUBROUTINE fit_structure_factor


 SUBROUTINE eval_residuals_xc_corr(ai,residual)
!------------------------------------------------------------------------!
! Evaluate vector of residuals for NL2SOL/fit_structure_factor_exp.      !
!------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: ai(:)
 REAL(dp),INTENT(out) :: residual(:)
 INTEGER i,j
 REAL(dp) poly1,poly2,expfn

 do i=1,ndata
  poly1=0.d0
  do j=1,nparam1
   poly1=poly1+ai(j)*gs(i)**(j+1)
  enddo
  poly2=0.d0
  do j=1,nparam2
   poly2=poly2+ai(nparam1+j)*gs(i)**(j+1)
  enddo
  if(poly1<maximum_exp_arg)then
   expfn=exp(poly1)
  else
   expfn=max_exp
  endif
  residual(i)=sf_sum(i)-((1.d0-expfn)*(1.d0+poly2))
 enddo
 if(USE_WEIGHTS)residual(1:ndata)=residual(1:ndata)*weight(1:ndata)

 END SUBROUTINE eval_residuals_xc_corr


 SUBROUTINE eval_jacobian_xc_corr(ai,jacobian)
!------------------------------------------------------------------------!
! Evaluate Jacobian matrix for NL2SOL/fit_structure_factor_exp.          !
!------------------------------------------------------------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: ai(:)
 REAL(dp),INTENT(out) :: jacobian(:)
 INTEGER i,j,ii,ii1

 epoly(:)=0.d0
 do i=2,ndata
  do j=1,nparam1
   epoly(i)=epoly(i)+ai(j)*gs(i)**(j+1)
  enddo
 enddo
 do i=1,ndata
  if(epoly(i)<maximum_exp_arg)then
   epoly(i)=exp(epoly(i))
  else
   epoly(i)=max_exp
  endif
 enddo
 polyb(:)=0.d0
 do i=2,ndata
  do j=1,nparam2
   polyb(i)=polyb(i)+ai(nparam1+j)*gs(i)**(j+1)
  enddo
 enddo
 ii=0
 do i=1,nparam1
  do j=1,ndata
   ii=ii+1
   jacobian(ii)=(1.d0+polyb(j))*epoly(j)*gs(j)**(i+1)
  enddo
 enddo
 ii1=ii
 if(USE_WEIGHTS)jacobian(1:ii)=jacobian(1:ii)*weight(j)
 do i=1,nparam2
  do j=1,ndata
   ii=ii+1
   jacobian(ii)=(epoly(j)-1.d0)*gs(j)**(i+1)
  enddo
 enddo
 if(USE_WEIGHTS)jacobian(ii1:ii)=jacobian(ii1:ii)*weight(j)

 END SUBROUTINE eval_jacobian_xc_corr


 SUBROUTINE overflow_protection
!--------------------------------------------------------------------------!
! Very basic overflow protection.                                          !
! Note that most machines handle exp() underflow by setting the result to  !
! zero, but define the relevant quantities here anyway.                    !
!                                                                          !
! MDT 2000                                                                 !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 largest_representable_number=huge(1.d0)  ! ~1.797693134862316E+308
 smallest_representable_number=tiny(1.d0) ! ~2.225073858507201E-308
 maximum_exp_arg_actual=log(largest_representable_number) ! ~709.782712893384
 minimum_exp_arg_actual=log(smallest_representable_number) ! ~-708.396418532264
! Define working limits to the size of arguments of the exp function based on
! an assumption of a largest representable number of 1.d300 (or whatever
! max_exp is set to in constants.f90). This is so we don't calculate
! f=exp(maximum_exp_arg_actual) then go on to multiply it by a number greater
! than 1 .
 maximum_exp_arg=log(max_exp)
 minimum_exp_arg=log(min_exp)
 if(maximum_exp_arg>=maximum_exp_arg_actual)call errstop('OVERFLOW_PROTECTION',&
  &'Largest representable double precision number on this machine is &
  & too small for the inbuilt exp() overflow protection to work.')
 if(minimum_exp_arg<=minimum_exp_arg_actual)call errstop('OVERFLOW_PROTECTION',&
  &'Smallest representable double precision number on this machine is &
  & too large for the inbuilt exp() underflow protection to work.')
 END SUBROUTINE overflow_protection


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
  d=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(2,1)*(A(3,2)*A(1,3)-A(1,2)*A(3,3))+&
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


END MODULE finsize


SUBROUTINE calcr_xc_corr(ndata,nparam,ai,nf,res,uiparm,urparm,ufparm)
!------------------------------------------------------------------------!
! Returns vector of residuals for NL2SOL/fit_structure_factor_exp.       !
!------------------------------------------------------------------------!
 USE dsp
 USE finsize,ONLY : eval_residuals_xc_corr
 IMPLICIT NONE
 INTEGER,INTENT(in) :: ndata,nparam
 INTEGER,INTENT(inout) :: nf
 INTEGER,INTENT(inout),OPTIONAL :: uiparm(:)
 REAL(dp),INTENT(in) :: ai(:)
 REAL(dp),INTENT(out) :: res(:)
 REAL(dp),INTENT(inout),OPTIONAL :: urparm(:),ufparm

 call eval_residuals_xc_corr(ai,res)

END SUBROUTINE calcr_xc_corr


SUBROUTINE calcj_xc_corr(ndata,nparam,ai,nf,jac,uiparm,urparm,ufparm)
!------------------------------------------------------------------------!
! Returns Jacobian matrix required by NL2SOL/fit_structure_factor_exp.   !
!------------------------------------------------------------------------!
 USE dsp
 USE finsize,ONLY : eval_jacobian_xc_corr
 IMPLICIT NONE
 INTEGER,INTENT(in) :: ndata,nparam
 INTEGER,INTENT(inout) :: nf
 INTEGER,INTENT(inout),OPTIONAL :: uiparm(:)
 REAL(dp),INTENT(in) :: ai(:)
 REAL(dp),INTENT(out) :: jac(:)
 REAL(dp),INTENT(inout),OPTIONAL :: urparm(:),ufparm

 call eval_jacobian_xc_corr(ai,jac)

END SUBROUTINE calcj_xc_corr


FUNCTION calcf_xc_corr(xx)
!------------------------------------------------------------------------!
! Evaluates fitted function for numerical integration following call of  !
! fit_structure_factor_exp (integrating S rather than S-1)               !
!------------------------------------------------------------------------!
 USE dsp
 USE finsize,ONLY : a,nparam1,nparam2
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: xx
 INTEGER i
 REAL(dp) p1,p2,calcf_xc_corr

 p1=0.d0
 do i=1,nparam1
  p1=p1+a(i)*xx**(i+1)
 enddo
 p2=0.d0
 do i=1,nparam2
  p2=p2+a(nparam1+i)*xx**(i+1)
 enddo
 calcf_xc_corr=(1.d0-exp(p1))*(1.d0+p2)

END FUNCTION calcf_xc_corr


FUNCTION calcf_xc_corr2(xx)
!------------------------------------------------------------------------!
! Evaluates fitted function for numerical integration following call of  !
! fit_structure_factor_exp (integrating S-1 rather than S)               !
!------------------------------------------------------------------------!
 USE dsp
 USE finsize,ONLY : a,nparam1,nparam2
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: xx
 INTEGER i
 REAL(dp) p1,p2,calcf_xc_corr2

 p1=0.d0
 do i=1,nparam1
  p1=p1+a(i)*xx**(i+1)
 enddo
 p2=0.d0
 do i=1,nparam2
  p2=p2+a(nparam1+i)*xx**(i+1)
 enddo
 calcf_xc_corr2=(1.d0-exp(p1))*(1.d0+p2)-1.d0

END FUNCTION calcf_xc_corr2


!FUNCTION calcf_test(xx)
!------------------------------------------------------------------------!
! Evaluates simple known function for Gaussian quadrature testing.       !
!------------------------------------------------------------------------!
! USE dsp
! IMPLICIT NONE
! REAL(dp),INTENT(in) :: xx
! REAL(dp) calcf_test
!
! e.g. integrate x^2 from 0 to 4 - should give [x^3/3] = 21.3333333333333333
!
! calcf_test=xx**2
!
!END FUNCTION calcf_test

PROGRAM format_configs
!--------------------------------------------------------------------------!
! Program for formatting a config.in file to produce a config.in_formatted !
! file or vice versa.                                                      !
!--------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0),sp=kind(1.0)
 LOGICAL lcon_in,lcon_out,lcon_backup,lcon_inf,lcon_outf,lcon_backupf,&
  &from_bin(3)
 CHARACTER(32) :: con_in='config.in',con_out='config.out',&
  &con_backup='config.backup',con_inf='config.in_formatted',&
  &con_outf='config.out_formatted',con_backupf='config.backup_formatted'
 CHARACTER(32),DIMENSION(3,2) :: con_choice
 INTEGER ichoice,nchoice,ierr

 inquire(file=trim(con_in),exist=lcon_in)
 inquire(file=trim(con_inf),exist=lcon_inf)
 inquire(file=trim(con_out),exist=lcon_out)
 inquire(file=trim(con_outf),exist=lcon_outf)
 inquire(file=trim(con_backup),exist=lcon_backup)
 inquire(file=trim(con_backupf),exist=lcon_backupf)

 if(lcon_in.and.lcon_inf)write(6,*)'Found both ',trim(con_in),' and ',&
  &trim(con_inf),'. Ignored.'
 if(lcon_out.and.lcon_outf)write(6,*)'Found both ',trim(con_out),' and ',&
  &trim(con_outf),'. Ignored.'
 if(lcon_backup.and.lcon_backupf)write(6,*)'Found both ',trim(con_backup),&
  &' and ',trim(con_backupf),'. Ignored.'

 nchoice=0
 if(lcon_in.neqv.lcon_inf)then
  nchoice=nchoice+1
  if(lcon_in)then
   con_choice(nchoice,1)=con_in ; con_choice(nchoice,2)=con_inf
   from_bin(nchoice)=.true.
  else
   con_choice(nchoice,1)=con_inf ; con_choice(nchoice,2)=con_in
   from_bin(nchoice)=.false.
  endif
 endif
 if(lcon_out.neqv.lcon_outf)then
  nchoice=nchoice+1
  if(lcon_out)then
   con_choice(nchoice,1)=con_out ; con_choice(nchoice,2)=con_outf
   from_bin(nchoice)=.true.
  else
   con_choice(nchoice,1)=con_outf ; con_choice(nchoice,2)=con_out
   from_bin(nchoice)=.false.
  endif
 endif
 if(lcon_backup.neqv.lcon_backupf)then
  nchoice=nchoice+1
  if(lcon_backup)then
   con_choice(nchoice,1)=con_backup ; con_choice(nchoice,2)=con_backupf
   from_bin(nchoice)=.true.
  else
   con_choice(nchoice,1)=con_backupf ; con_choice(nchoice,2)=con_backup
   from_bin(nchoice)=.false.
  endif
 endif

 select case(nchoice)
 case(0) ; write(6,*)'No config files to format.' ; stop
 case(1) ; ichoice=1
 case default
  do ichoice=1,nchoice
   write(6,*)trim(i2s(ichoice)),': ',trim(con_choice(ichoice,1)),' -> ',&
    &trim(con_choice(ichoice,2))
  enddo
  write(6,*)'Choose one of the above options:'
  read(5,*,iostat=ierr)ichoice
  if(ierr/=0)then
   write(6,*)'Not a valid option. Quitting.' ; stop
  endif
  if(ichoice<1.or.ichoice>nchoice)then
   write(6,*)'Not a valid option. Quitting.' ; stop
  endif
 end select
 write(6,*)'Performing conversion ',trim(con_choice(ichoice,1)),' -> ',&
  &trim(con_choice(ichoice,2))
 if(from_bin(ichoice))then
  call b2a_config_file(con_choice(ichoice,1),con_choice(ichoice,2))
 else
  call a2b_config_file(con_choice(ichoice,1),con_choice(ichoice,2))
 endif
 write(6,*)'Done.'


CONTAINS


 SUBROUTINE check_alloc(ialloc,rout,var)
  INTEGER,INTENT(in) :: ialloc
  CHARACTER(*),INTENT(in) :: rout,var
  if(ialloc/=0)call errstop(rout,"allocation problem ("//var//")")
 END SUBROUTINE check_alloc


 SUBROUTINE b2a_config_file(con_from,con_to)
!---------------------------------------------------!
! Clone of copy_config_file, with formatted output. !
!---------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: con_from,con_to
  CHARACTER(20) label
  INTEGER io_from,io_to,ierr,ialloc,inode
! Data in config file
  INTEGER nspin_in,ndet_in,nemax_in,netot_in,nitem_in,item,nextra_in,&
   &no_difftypes_in,vmc_steps_in,r1_c2_in,nconfig_in,iconfig,complex_wf_in, &
   &noncoll_spin_in
  INTEGER,ALLOCATABLE :: nele_in(:),sele_vmc_in(:),sele_in(:)
  COMPLEX(dp),ALLOCATABLE :: logdet_in(:,:)
  CHARACTER(20) gen_by_in,interaction_in,atom_basis_type_in
  CHARACTER(20),ALLOCATABLE :: config_item_in(:)
  REAL(sp) total_time_in
  REAL(dp) lapdet_in,etot_in,local_potential_in,nltot_in,stot_in,wdmc_in,&
   &valjas_in,logp_in,twist_in(5)
  REAL(dp),ALLOCATABLE :: dtvmc_array_in(:),rele_vmc_in(:,:),rele_in(:,:),&
   &fidet_in(:,:,:),fi_prod_det_in(:,:,:,:,:),prod_lapdet_in(:,:,:,:)
  INTEGER dmcequil_steps_in,dmcstats_steps_in,tpdmc_in,nnodes_in,&
   &random_state_in(25),nitot_in
  INTEGER ranx_max,ranx_gauss_max,ranx_indx,ranx_gauss_indx
  REAL(sp),ALLOCATABLE :: ranx_buffer(:),gauss_buffer(:)
  INTEGER lwdmc_in,growth_estimator_in
  REAL(dp) ebest_in,ebest_init_in,eref_in,dteff_ebest_init_in,dteff_best_in,&
   &numerator_wt2_in,denominator_wt_in,denominator_wt2_in,log_pi_wt_in,&
   &log_pi_wt2_in,numer_expect_in(13),final_vmcE_in,final_vmcdE_in,&
   &final_vmcdEu_in,final_vmcE2_in,final_vmcdE2_in,final_vmcdE2u_in,&
   &final_vmcvar_in,final_vmcmove_in,final_vmctau_in
  REAL(dp),ALLOCATABLE :: log_pi_wt_array_in(:),log_pi_wt_array2_in(:),&
   &rion_in(:,:)
! On-the-fly reblocking data
  INTEGER reblock_nbs_in,reblock_nstep_in,reblock_nobs_in,popstats_in
  REAL(dp) reblock_sum_w_in,reblock_sum_w2_in
  INTEGER,ALLOCATABLE :: reblock_block_length_in(:),reblock_s_u_f_in(:)
  REAL(dp),ALLOCATABLE :: reblock_sum_w2_closed_in(:),reblock_sum_w_open_in(:),&
   &reblock_sum_o2_closed_in(:,:),reblock_sum_o_open_in(:,:),&
   &reblock_sum_o_in(:),reblock_sum_ow0_in(:),reblock_sum_ow2_in(:),&
   &reblock_sum_o2w0_in(:),reblock_sum_o2w2_in(:)

! Initialize data
  gen_by_in='NONE' ; nitem_in=0 ; nconfig_in=-1
  nspin_in=-1 ; ndet_in=-1 ; nnodes_in=0
  complex_wf_in=0 ; r1_c2_in=1
  noncoll_spin_in=0
  no_difftypes_in=-1
  io_from=10 ; io_to=11

! Open files
  open(unit=io_from,file=trim(con_from),form='unformatted',status='old',&
   &iostat=ierr)
  if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem opening '//&
   &trim(con_from)//' for copy.')
  open(unit=io_to,file=trim(con_to),form='formatted',status='new',&
   &iostat=ierr)
  if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem opening '//trim(con_to)//&
   &'.')

! Read first label to check for correct file format.
  read(io_from,iostat=ierr)label
  if(ierr<0)call errstop('B2A_CONFIG_FILE','Config file appears to be empty?')
  if(ierr>0)call errstop('B2A_CONFIG_FILE','Problem reading first label in &
   &config file. Possibly it was produced on a different computer with an &
   &incompatible binary format?')
  if(trim(label)/='INFO')call errstop('B2A_CONFIG_FILE','First label in &
   &config file is wrong. Perhaps this config file is in old format? Try &
   &using the update_config utility.')
  rewind(io_from)

! Start copy

! Read each section
  do
   read(io_from,iostat=ierr)label
   if(ierr<0)then
    if(nconfig_in<0)then
     nconfig_in=0 ; exit
    endif
    call errstop('B2A_CONFIG_FILE','End-of-file reached while reading section &
     &name.')
   endif
   if(ierr>0)call errstop('B2A_CONFIG_FILE','Problem reading section name.')
   write(io_to,*)label
   select case(trim(label))
   case('INFO') ! system/run info, mostly for error-checking
    do
     read(io_from,iostat=ierr)label
     if(ierr<0)call errstop('B2A_CONFIG_FILE','End-of-file reached while &
      &reading info label.')
     if(ierr>0)call errstop('B2A_CONFIG_FILE','Problem reading info label.')
     write(io_to,*)label
     if(trim(label)=='END INFO')exit
     select case(trim(label))
     case('NEXTRA') ! number of 'extra' sections in this file
      read(io_from,iostat=ierr)nextra_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NEXTRA.')
      write(io_to,*)nextra_in
     case('GEN_BY')
      read(io_from,iostat=ierr)gen_by_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading GEN_BY.')
      write(io_to,*)gen_by_in
     case('NSPIN')
      read(io_from,iostat=ierr)nspin_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NSPIN.')
      write(io_to,*)nspin_in
     case('NELE')
      if(nspin_in<0)call errstop('B2A_CONFIG_FILE','NELE found before NSPIN &
       &in config file. This is a bug.')
      allocate(nele_in(nspin_in),stat=ialloc)
      if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Problem allocating &
       &NELE_IN.')
      read(io_from,iostat=ierr)nele_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NELE.')
      write(io_to,*)nele_in
      nemax_in=maxval(nele_in) ; netot_in=sum(nele_in)
      deallocate(nele_in)
     case('NDET')
      read(io_from,iostat=ierr)ndet_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NDET.')
      write(io_to,*)ndet_in
     case('COMPLEX_WF')
      read(io_from,iostat=ierr)complex_wf_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading COMPLEX_WF.')
      if(complex_wf_in/=0)r1_c2_in=2
      write(io_to,*)complex_wf_in
     case('NONCOLL_SPIN')
      read(io_from,iostat=ierr)noncoll_spin_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NONCOLL_SPIN.')
      write(io_to,*)noncoll_spin_in
     case('INTERACTION')
      read(io_from,iostat=ierr)interaction_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading INTERACTION.')
      write(io_to,*)interaction_in
     case('ATOM_BASIS_TYPE')
      read(io_from,iostat=ierr)atom_basis_type_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &ATOM_BASIS_TYPE.')
      write(io_to,*)atom_basis_type_in
     case('TOTAL_TIME')
      read(io_from,iostat=ierr)total_time_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &TOTAL_TIME.')
      write(io_to,*)total_time_in
     case default ! unknown
      call errstop('B2A_CONFIG_FILE','Found unknown info '//trim(label)//'.')
     end select
    enddo
   case('DEFINE_CONFIGS')
    read(io_from,iostat=ierr)nconfig_in
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading number of &
     &configs.')
    write(io_to,*)nconfig_in
    read(io_from,iostat=ierr)nitem_in
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading number of &
     &items per config.')
    write(io_to,*)nitem_in
    allocate(config_item_in(nitem_in),stat=ialloc)
    if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Allocation problem &
     &(config_item).')
    do item=1,nitem_in
     read(io_from,iostat=ierr)label
     if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading config &
      &content label.')
     write(io_to,*)label
     config_item_in(item)=trim(label)
     select case(trim(label))
     case('RELE','SELE','ETOT','LOGDET','FIDET','FI_PROD_DET','LAPDET',&
      &'PROD_LAPDET','LOCAL_POTENTIAL','NLTOT','STOT','WDMC','VALJAS','LOGP',&
      &'TWIST')
      continue
     case default
      call errstop('B2A_CONFIG_FILE',"Unrecognized item '"//trim(label)//"'.")
     end select
    enddo
    read(io_from,iostat=ierr)label
    if(ierr/=0)call errstop('B2A_CONFIG_FILE',"Problem reading 'END &
     &DEFINE_CONFIGS'.")
    if(trim(label)/='END DEFINE_CONFIGS')call errstop('B2A_CONFIG_FILE',&
     &"Expected to find 'END DEFINE_CONFIGS' but didn't.")
    write(io_to,*)label
   case('VMC_SAVED_STATE')
    do
     read(io_from,iostat=ierr)label
     if(ierr/=0)call errstop('B2A_CONFIG_FILE',"Problem reading label in &
      &VMC_SAVED_STATE section in config file.")
     write(io_to,*)label
     select case(trim(label))
     case('NO_DIFFTYPES')
      read(io_from,iostat=ierr)no_difftypes_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &NO_DIFFTYPES.')
      write(io_to,*)no_difftypes_in
     case('DTVMC')
      if(no_difftypes_in<1)call errstop('B2A_CONFIG_FILE','DTVMC found &
       &before NO_DIFFTYPES in config file.')
      allocate(dtvmc_array_in(no_difftypes_in),stat=ialloc)
      if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Allocation problem &
       &(DTVMC_ARRAY_IN).')
      read(io_from,iostat=ierr)dtvmc_array_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading DTVMC.')
      write(io_to,*)dtvmc_array_in
      deallocate(dtvmc_array_in)
     case('VMC_STEPS')
      read(io_from,iostat=ierr)vmc_steps_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading VMC_STEPS.')
      write(io_to,*)vmc_steps_in
     case('NNODES')
      read(io_from,iostat=ierr)nnodes_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NNODES.')
      write(io_to,*)nnodes_in
     case('RELE')
      if(netot_in<1)call errstop('B2A_CONFIG_FILE','NETOT undefined when &
       &reading RELE VMC_SAVED_STATE block. This is a bug.')
      allocate(rele_vmc_in(3,netot_in),stat=ialloc)
      if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Allocation problem &
       &(RELE_VMC_IN).')
      do inode=1,nnodes_in
       read(io_from,iostat=ierr)rele_vmc_in
       if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading RELE_VMC.')
       write(io_to,*)rele_vmc_in
      enddo ! inode
      deallocate(rele_vmc_in)
     case('SELE')
      if(netot_in<1)call errstop('B2A_CONFIG_FILE','NETOT undefined when &
       &reading SELE in VMC_SAVED_STATE block. This is a bug.')
      allocate(sele_vmc_in(netot_in),stat=ialloc)
      if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Allocation problem &
       &(SELE_VMC_CONFIG).')
      do inode=1,nnodes_in
       read(io_from,iostat=ierr)sele_vmc_in
       if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading SELE.')
       write(io_to,*)sele_vmc_in
      enddo ! inode
      deallocate(sele_vmc_in)
     case('END VMC_SAVED_STATE')
      exit
     case default
      call errstop('B2A_CONFIG_FILE','Unrecognized information '//&
       &trim(label)//' in VMC_SAVED_STATE section of config file.')
     end select
    enddo
   case('FINAL_VMC_RESULT')
    do
     read(io_from,iostat=ierr)label
     if(ierr/=0)call errstop('B2A_CONFIG_FILE',"Problem reading label &
      &in FINAL_VMC_RESULT section in config file.")
     write(io_to,*)label
     select case(trim(label))
     case('ENERGY')
      read(io_from,iostat=ierr)final_vmcE_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading ENERGY')
      write(io_to,*)final_vmcE_in
     case('ERRORBAR')
      read(io_from,iostat=ierr)final_vmcdE_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading ERRORBAR')
      write(io_to,*)final_vmcdE_in
     case('ERRORBARU')
      read(io_from,iostat=ierr)final_vmcdEu_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading ERRORBARU')
      write(io_to,*)final_vmcdEu_in
     case('VARIANCE')
      read(io_from,iostat=ierr)final_vmcvar_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading VARIANCE')
      write(io_to,*)final_vmcvar_in
     case('CORRTIME')
      read(io_from,iostat=ierr)final_vmctau_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading CORRTIME')
      write(io_to,*)final_vmctau_in
     case('TOTMOVE')
      read(io_from,iostat=ierr)final_vmcmove_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading TOTMOVE')
      write(io_to,*)final_vmcmove_in
     case('ENERGY2')
      read(io_from,iostat=ierr)final_vmcE2_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading ENERGY2')
      write(io_to,*)final_vmcE2_in
     case('ERRORBAR2')
      read(io_from,iostat=ierr)final_vmcdE2_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading ERRORBAR2')
      write(io_to,*)final_vmcdE2_in
     case('ERRORBAR2U')
      read(io_from,iostat=ierr)final_vmcdE2u_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading ERRORBAR2U')
      write(io_to,*)final_vmcdE2u_in
     case('END FINAL_VMC_RESULT')
      exit
     case default
      call errstop('B2A_CONFIG_FILE','Unrecognized information '//&
       &trim(label)//' in FINAL_VMC_RESULT section of config file.')
     end select
    enddo
   case('DMC_SAVED_STATE')
    do
     read(io_from,iostat=ierr)label
     if(ierr/=0)call errstop('B2A_CONFIG_FILE',"Problem reading label &
      &in VMC_SAVED_STATE section in config file.")
     write(io_to,*)label
     select case(trim(label))
     case('EBEST')
      read(io_from,iostat=ierr)ebest_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading EBEST.')
      write(io_to,*)ebest_in
     case('EBEST_INIT')
      read(io_from,iostat=ierr)ebest_init_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading EBEST_INIT.')
      write(io_to,*)ebest_init_in
     case('EREF')
      read(io_from,iostat=ierr)eref_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading EREF.')
      write(io_to,*)eref_in
     case('DTEFF_EBEST_INIT')
      read(io_from,iostat=ierr)dteff_ebest_init_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &DTEFF_EBEST_INIT.')
      write(io_to,*)dteff_ebest_init_in
     case('DTEFF_BEST')
      read(io_from,iostat=ierr)dteff_best_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading DTEFF_BEST.')
      write(io_to,*)dteff_best_in
     case('DMCEQUIL_STEPS')
      read(io_from,iostat=ierr)dmcequil_steps_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &DMCEQUIL_STEPS.')
      write(io_to,*)dmcequil_steps_in
     case('DMCSTATS_STEPS')
      read(io_from,iostat=ierr)dmcstats_steps_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &DMCSTATS_STEPS.')
      write(io_to,*)dmcstats_steps_in
     case('TPDMC')
      read(io_from,iostat=ierr)tpdmc_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading TPDMC.')
      write(io_to,*)tpdmc_in
     case('LWDMC')
      read(io_from,iostat=ierr)lwdmc_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading LWDMC.')
      write(io_to,*)lwdmc_in
     case('GROWTH_ESTIMATOR')
      read(io_from,iostat=ierr)growth_estimator_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &GROWTH_ESTIMATOR.')
      write(io_to,*)growth_estimator_in
     case('NUMER_EXPECT')
      read(io_from,iostat=ierr)numer_expect_in(1:13)
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &NUMER_EXPECT.')
      write(io_to,*)numer_expect_in(1:13)
     case('DENOMINATOR_WT')
      read(io_from,iostat=ierr)denominator_wt_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &DENOMINATOR_WT.')
      write(io_to,*)denominator_wt_in
     case('LOG_PI_WT')
      read(io_from,iostat=ierr)log_pi_wt_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading LOG_PI_WT.')
      write(io_to,*)log_pi_wt_in
     case('LOG_PI_WT2')
      read(io_from,iostat=ierr)log_pi_wt2_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading LOG_PI_WT2.')
      write(io_to,*)log_pi_wt2_in
     case('LOG_PI_WT_ARRAY')
      if(tpdmc_in<0)call errstop('B2A_CONFIG_FILE','LOG_PI_WT_ARRAY found &
       &before TPDMC in config file. This is a bug.')
      if(tpdmc_in==0)call errstop('B2A_CONFIG_FILE','LOG_PI_WT_ARRAY found &
       &in config file even though the stored TPDMC is zero. This is a bug.')
      allocate(log_pi_wt_array_in(0:tpdmc_in-1),stat=ialloc)
      if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Allocation problem &
       &(LOG_PI_WT_ARRAY)')
      read(io_from,iostat=ierr)log_pi_wt_array_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &LOG_PI_WT_ARRAY.')
      write(io_to,*)log_pi_wt_array_in
      deallocate(log_pi_wt_array_in)
     case('LOG_PI_WT_ARRAY2')
      if(tpdmc_in<0)call errstop('B2A_CONFIG_FILE',&
       &'LOG_PI_WT_ARRAY2 found before TPDMC in config file. This is a bug.')
      allocate(log_pi_wt_array2_in(0:tpdmc_in),stat=ialloc)
      if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Allocation problem &
       &(LOG_PI_WT_ARRAY2)')
      read(io_from,iostat=ierr)log_pi_wt_array2_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &LOG_PI_WT_ARRAY2.')
      write(io_to,*)log_pi_wt_array2_in
      deallocate(log_pi_wt_array2_in)
     case('NUMERATOR_WT2')
      read(io_from,iostat=ierr)numerator_wt2_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &NUMERATOR_WT2.')
      write(io_to,*)numerator_wt2_in
     case('DENOMINATOR_WT2')
      read(io_from,iostat=ierr)denominator_wt2_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &DENOMINATOR_WT2.')
      write(io_to,*)denominator_wt2_in
     case('END DMC_SAVED_STATE')
      exit
     case default
      call errstop('B2A_CONFIG_FILE','Unrecognized information '//&
       &trim(label)//' in DMC_SAVED_STATE section of config file.')
     end select
    enddo
   case('REBLOCK_DATA')
    read(io_from,iostat=ierr)reblock_nbs_in,reblock_nstep_in,&
     &reblock_nobs_in,reblock_sum_w_in,popstats_in
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading number of &
     &reblock blocks.')
    write(io_to,*)reblock_nbs_in,reblock_nstep_in,&
     &reblock_nobs_in,reblock_sum_w_in,popstats_in
    allocate(reblock_block_length_in(reblock_nbs_in),&
     &reblock_s_u_f_in(reblock_nbs_in),&
     &reblock_sum_w2_closed_in(reblock_nbs_in),&
     &reblock_sum_w_open_in(reblock_nbs_in),&
     &reblock_sum_o2_closed_in(reblock_nobs_in,reblock_nbs_in),&
     &reblock_sum_o_open_in(reblock_nobs_in,reblock_nbs_in),&
     &reblock_sum_o_in(reblock_nobs_in),&
     &stat=ialloc)
    call check_alloc(ialloc,'B2A_CONFIG_FILE','reblock_block_length_config')
    read(io_from,iostat=ierr)reblock_block_length_in,&
     &reblock_s_u_f_in,&
     &reblock_sum_w2_closed_in,&
     &reblock_sum_w_open_in,&
     &reblock_sum_o2_closed_in,&
     &reblock_sum_o_open_in,&
     &reblock_sum_o_in
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading reblock blocks.')
    write(io_to,*)reblock_block_length_in,&
     &reblock_s_u_f_in,&
     &reblock_sum_w2_closed_in,&
     &reblock_sum_w_open_in,&
     &reblock_sum_o2_closed_in,&
     &reblock_sum_o_open_in,&
     &reblock_sum_o_in
    deallocate(reblock_block_length_in,reblock_s_u_f_in,&
     &reblock_sum_w2_closed_in,reblock_sum_w_open_in,reblock_sum_o2_closed_in,&
     &reblock_sum_o_open_in,reblock_sum_o_in)
    if(popstats_in/=0)then
     allocate( reblock_sum_ow0_in(reblock_nobs_in),&
      &reblock_sum_ow2_in(reblock_nobs_in),&
      &reblock_sum_o2w0_in(reblock_nobs_in),&
      &reblock_sum_o2w2_in(reblock_nobs_in),&
      &stat=ialloc)
     call check_alloc(ialloc,'B2A_CONFIG_FILE','reblock_sum_ow0_in')
     read(io_from,iostat=ierr)reblock_sum_w2_in,&
      &reblock_sum_ow0_in,&
      &reblock_sum_ow2_in,&
      &reblock_sum_o2w0_in,&
      &reblock_sum_o2w2_in
     if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading popstats.')
     write(io_to,*)reblock_sum_w2_in,&
      &reblock_sum_ow0_in,&
      &reblock_sum_ow2_in,&
      &reblock_sum_o2w0_in,&
      &reblock_sum_o2w2_in
     deallocate( reblock_sum_ow0_in,reblock_sum_ow2_in,&
      &reblock_sum_o2w0_in,reblock_sum_o2w2_in)
    endif
    read(io_from,iostat=ierr)label
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading label END &
     &REBLOCK_DATA in config file.')
    if(trim(label)/='END REBLOCK_DATA')call errstop('B2A_CONFIG_FILE',&
     &'Wrong label "'//trim(label)//'" instead of END REBLOCK_DATA in config &
     &file.')
    write(io_to,*)label
   case('RANDOM')
    read(io_from,iostat=ierr)label
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading label NNODES &
     &in RANDOM section of config file.')
    if(trim(label)/='NNODES')call errstop('B2A_CONFIG_FILE','Wrong label &
     &"'//trim(label)//'" instead of NNODES in RANDOM section of config file.')
    write(io_to,*)label
    read(io_from,iostat=ierr)nnodes_in
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NNODES in &
     &RANDOM section of config file.')
    write(io_to,*)nnodes_in
    do inode=1,nnodes_in
     read(io_from,iostat=ierr)random_state_in
     if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading RANDOM_STATE &
      &for node '//trim(i2s(inode))//' in RANDOM section of config file.')
     write(io_to,*)random_state_in
    enddo
    read(io_from,iostat=ierr)label
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading RANDOM &
     &section of config file.')
    if(trim(label)/='BUFFER'.and.trim(label)/='END RANDOM')call &
     &errstop('B2A_CONFIG_FILE','Wrong label "'//trim(label)//'" instead of&
     & BUFFER or END RANDOM in config file.')
    if(trim(label)=='BUFFER')then
     write(io_to,*)label
     read(io_from,iostat=ierr)ranx_max,ranx_gauss_max
     if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading BUFFER sizes &
      & in RANDOM section of config file.')
     write(io_to,*)ranx_max,ranx_gauss_max
     do inode=1,nnodes_in
      read(io_from,iostat=ierr)ranx_indx,ranx_gauss_indx
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading BUFFER &
       & positions in RANDOM section of config file.')
      write(io_to,*)ranx_indx,ranx_gauss_indx
     enddo
     allocate(ranx_buffer(ranx_max),stat=ierr)
     if(ierr/=0)call errstop('B2A_CONFIG_FILE','Allocation problem for &
      &ranx_buffer.')
     do inode=1,nnodes_in
      read(io_from,iostat=ierr)ranx_buffer
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading RANX_BUFFER &
       & in RANDOM section of config file.')
      write(io_to,*)ranx_buffer
     enddo
     deallocate(ranx_buffer)
     allocate(gauss_buffer(ranx_gauss_max),stat=ierr)
     if(ierr/=0)call errstop('B2A_CONFIG_FILE','Allocation problem for &
      &gauss_buffer.')
     do inode=1,nnodes_in
      read(io_from,iostat=ierr)gauss_buffer
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading GAUSS_BUFFER &
       &in RANDOM section of config file.')
      write(io_to,*)gauss_buffer
     enddo
     deallocate(gauss_buffer)
     read(io_from,iostat=ierr)label
     if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading END RANDOM &
      &in config file.')
     if(trim(label)/='END RANDOM')call errstop('B2A_CONFIG_FILE',&
      &'Wrong label "'//trim(label)//'" instead of END RANDOM in config file.')
    endif
    write(io_to,*)label
   case('GEOMETRY')
    read(io_from,iostat=ierr)nitot_in
    if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NITOT.')
    write(io_to,*)nitot_in
    if(nitot_in>0)then
     allocate(rion_in(3,nitot_in),stat=ialloc)
     if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Allocation problem RION_IN')
     read(io_from,iostat=ierr)rion_in
     if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading RION.')
     write(io_to,*)rion_in
     deallocate(rion_in)
    endif ! nitot_in>0
    read(io_from,iostat=ierr)label
    if(ierr/=0)call errstop('B2A_CONFIG_FILE',"Problem reading 'END &
     &GEOMETRY' label in config file.")
    if(trim(label)/='END GEOMETRY')call errstop('B2A_CONFIG_FILE','Wrong &
     &label "'//trim(label)//'" read in when expecting "END GEOMETRY".')
    write(io_to,*)label
   case('CONFIGS')
    exit
   case default
    call errstop('B2A_CONFIG_FILE','Unrecognized section '//trim(label)//&
     &' in config file.')
   end select
  enddo

! Configs, only if config label found
  if(nconfig_in<0)call errstop('B2A_CONFIG_FILE','Number of configs undefined &
   &when about to read configs.')

  if(nconfig_in>0)then

! Check that configs are defined before loading them
   if(nitem_in==0)call errstop('B2A_CONFIG_FILE','Config contents undefined &
    &when about to read configs.')

! Allocate arrays.
   do item=1,nitem_in
    select case(trim(config_item_in(item)))
    case('RELE') ; allocate(rele_in(3,netot_in),stat=ialloc)
    case('SELE') ; allocate(sele_in(netot_in),stat=ialloc)
    case('LOGDET') ; allocate(logdet_in(nspin_in,ndet_in),stat=ialloc)
    case('FIDET') ; allocate(fidet_in(3,netot_in,r1_c2_in),stat=ialloc)
    case('FI_PROD_DET') ; allocate(fi_prod_det_in&
     &(3,ndet_in,nemax_in,r1_c2_in,nspin_in),stat=ialloc)
    case('PROD_LAPDET') ; allocate(prod_lapdet_in&
     &(ndet_in,nemax_in,r1_c2_in,nspin_in),stat=ialloc)
    case default ; ialloc=0
    end select
    if(ialloc/=0)call errstop('B2A_CONFIG_FILE','Allocation problem ('//&
     &trim(config_item_in(item))//').')
   enddo

! Copy configs.
   do iconfig=1,nconfig_in
    do item=1,nitem_in
     select case(trim(config_item_in(item)))
     case('RELE')
      read(io_from,iostat=ierr)rele_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading RELE for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)rele_in
     case('SELE')
      read(io_from,iostat=ierr)sele_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading SELE for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)sele_in
     case('ETOT')
      read(io_from,iostat=ierr)etot_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading ETOT for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)etot_in
     case('LOGDET')
      read(io_from,iostat=ierr)logdet_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading LOGDET for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)logdet_in
     case('FIDET')
      read(io_from,iostat=ierr)fidet_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading FIDET for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)fidet_in
     case('FI_PROD_DET')
      read(io_from,iostat=ierr)fi_prod_det_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading FI_PROD_DET &
       &for config '//trim(i2s(iconfig))//'.')
      write(io_to,*)fi_prod_det_in
     case('LAPDET')
      read(io_from,iostat=ierr)lapdet_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading LAPDET for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)lapdet_in
     case('PROD_LAPDET')
      read(io_from,iostat=ierr)prod_lapdet_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading PROD_LAPDET &
       &for config '//trim(i2s(iconfig))//'.')
      write(io_to,*)prod_lapdet_in
     case('LOCAL_POTENTIAL')
      read(io_from,iostat=ierr)local_potential_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading &
       &LOCAL_POTENTIAL for config '//trim(i2s(iconfig))//'.')
      write(io_to,*)local_potential_in
     case('NLTOT')
      read(io_from,iostat=ierr)nltot_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading NLTOT for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)nltot_in
     case('STOT')
      read(io_from,iostat=ierr)stot_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading STOT for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)stot_in
     case('WDMC')
      read(io_from,iostat=ierr)wdmc_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading WDMC for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)wdmc_in
     case('VALJAS')
      read(io_from,iostat=ierr)valjas_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading VALJAS for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)valjas_in
     case('LOGP')
      read(io_from,iostat=ierr)logp_in
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading LOGP for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)logp_in
     case('TWIST')
      read(io_from,iostat=ierr)twist_in(:)
      if(ierr/=0)call errstop('B2A_CONFIG_FILE','Problem reading TWIST for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to,*)twist_in(:)
     end select
    enddo ! item
   enddo ! iconfig

! Deallocate temp arrays
   do item=1,nitem_in
    select case(trim(config_item_in(item)))
    case('RELE') ; deallocate(rele_in)
    case('SELE') ; deallocate(sele_in)
    case('LOGDET') ; deallocate(logdet_in)
    case('FIDET') ; deallocate(fidet_in)
    case('FI_PROD_DET') ; deallocate(fi_prod_det_in)
    case('PROD_LAPDET') ; deallocate(prod_lapdet_in)
    end select
   enddo

  endif ! nconfig_in>0

  if(allocated(config_item_in))deallocate(config_item_in)

! Close files
  close(io_from)
  close(io_to)

 END SUBROUTINE b2a_config_file


 SUBROUTINE a2b_config_file(con_from,con_to)
!--------------------------------------------------!
! Clone of copy_config_file, with formatted input. !
!--------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: con_from,con_to
  CHARACTER(20) label
  INTEGER io_from,io_to,ierr,ialloc,inode
! Data in config file
  INTEGER nspin_in,ndet_in,nemax_in,netot_in,nitem_in,item,nextra_in,&
   &no_difftypes_in,vmc_steps_in,r1_c2_in,nconfig_in,iconfig,complex_wf_in, &
   &noncoll_spin_in
  INTEGER,ALLOCATABLE :: nele_in(:),sele_vmc_in(:),sele_in(:)
  COMPLEX(dp),ALLOCATABLE :: logdet_in(:,:)
  CHARACTER(20) gen_by_in,interaction_in,atom_basis_type_in
  CHARACTER(20),ALLOCATABLE :: config_item_in(:)
  REAL(sp) total_time_in
  REAL(dp) lapdet_in,etot_in,local_potential_in,nltot_in,stot_in,wdmc_in,&
   &valjas_in,logp_in,twist_in(5)
  REAL(dp),ALLOCATABLE :: dtvmc_array_in(:),rele_vmc_in(:,:),rele_in(:,:),&
   &fidet_in(:,:,:),fi_prod_det_in(:,:,:,:,:),prod_lapdet_in(:,:,:,:)
  INTEGER dmcequil_steps_in,dmcstats_steps_in,tpdmc_in,nnodes_in,&
   &random_state_in(25),nitot_in
  INTEGER ranx_max,ranx_gauss_max,ranx_indx,ranx_gauss_indx
  REAL(sp),ALLOCATABLE :: ranx_buffer(:),gauss_buffer(:)
  INTEGER lwdmc_in,growth_estimator_in
  REAL(dp) ebest_in,ebest_init_in,eref_in,dteff_ebest_init_in,dteff_best_in,&
   &numerator_wt2_in,denominator_wt_in,denominator_wt2_in,log_pi_wt_in,&
   &log_pi_wt2_in,numer_expect_in(13),final_vmcE_in,final_vmcdE_in,&
   &final_vmcdEu_in,final_vmcE2_in,final_vmcdE2_in,final_vmcdE2u_in,&
   &final_vmcvar_in,final_vmcmove_in,final_vmctau_in
  REAL(dp),ALLOCATABLE :: log_pi_wt_array_in(:),log_pi_wt_array2_in(:),&
   &rion_in(:,:)
! On-the-fly reblocking data
  INTEGER reblock_nbs_in,reblock_nstep_in,reblock_nobs_in
  REAL(dp) reblock_sum_w_in,reblock_sum_w2_in
  INTEGER popstats_in
  INTEGER,ALLOCATABLE :: reblock_block_length_in(:),reblock_s_u_f_in(:)
  REAL(dp),ALLOCATABLE :: reblock_sum_w2_closed_in(:),reblock_sum_w_open_in(:),&
   &reblock_sum_o2_closed_in(:,:),reblock_sum_o_open_in(:,:),&
   &reblock_sum_o_in(:),reblock_sum_ow0_in(:),reblock_sum_ow2_in(:),&
   &reblock_sum_o2w0_in(:),reblock_sum_o2w2_in(:)

! Initialize data
  gen_by_in='NONE' ; nitem_in=0 ; nconfig_in=-1
  nspin_in=-1 ; ndet_in=-1 ; nnodes_in=0
  complex_wf_in=0 ; r1_c2_in=1
  noncoll_spin_in=0
  no_difftypes_in=-1
  io_from=10 ; io_to=11

! Open files
  open(unit=io_from,file=trim(con_from),form='formatted',status='old',&
   &iostat=ierr)
  if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem opening '//&
   &trim(con_from)//' for copy.')
  open(unit=io_to,file=trim(con_to),form='unformatted',status='new',&
   &iostat=ierr)
  if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem opening '//trim(con_to)//&
   &'.')

! Read first label to check for correct file format.
  read(io_from,'(a)',iostat=ierr)label
  if(ierr<0)call errstop('A2B_CONFIG_FILE','Config file appears to be empty?')
  if(ierr>0)call errstop('A2B_CONFIG_FILE','Problem reading first label in &
   &config file. Possibly it was produced on a different computer with an &
   &incompatible binary format?')
  label=adjustl(label)
  if(trim(label)/='INFO')call errstop('A2B_CONFIG_FILE','First label in &
   &config file is wrong. Perhaps this config file is in old format? Try &
   &using the update_config utility.')
  rewind(io_from)

! Start copy

! Read each section
  do
   read(io_from,'(a)',iostat=ierr)label
   if(ierr<0)then
    if(nconfig_in<0)then
     nconfig_in=0 ; exit
    endif
    call errstop('A2B_CONFIG_FILE','End-of-file reached while reading section &
     &name.')
   endif
   if(ierr>0)call errstop('A2B_CONFIG_FILE','Problem reading section name.')
   label=adjustl(label)
   write(io_to)label
   select case(trim(label))
   case('INFO') ! system/run info, mostly for error-checking
    do
     read(io_from,'(a)',iostat=ierr)label
     if(ierr<0)call errstop('A2B_CONFIG_FILE','End-of-file reached while &
      &reading info label.')
     if(ierr>0)call errstop('A2B_CONFIG_FILE','Problem reading info label.')
     label=adjustl(label)
     write(io_to)label
     if(trim(label)=='END INFO')exit
     select case(trim(label))
     case('NEXTRA') ! number of 'extra' sections in this file
      read(io_from,*,iostat=ierr)nextra_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NEXTRA.')
      write(io_to)nextra_in
     case('GEN_BY')
      read(io_from,'(a)',iostat=ierr)gen_by_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading GEN_BY.')
      gen_by_in=trim(adjustl(gen_by_in))
      write(io_to)gen_by_in
     case('NSPIN')
      read(io_from,*,iostat=ierr)nspin_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NSPIN.')
      write(io_to)nspin_in
     case('NELE')
      if(nspin_in<0)call errstop('A2B_CONFIG_FILE','NELE found before NSPIN &
       &in config file. This is a bug.')
      allocate(nele_in(nspin_in),stat=ialloc)
      if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Problem allocating &
       &NELE_IN.')
      read(io_from,*,iostat=ierr)nele_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NELE.')
      write(io_to)nele_in
      nemax_in=maxval(nele_in) ; netot_in=sum(nele_in)
      deallocate(nele_in)
     case('NDET')
      read(io_from,*,iostat=ierr)ndet_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NDET.')
      write(io_to)ndet_in
     case('COMPLEX_WF')
      read(io_from,*,iostat=ierr)complex_wf_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading COMPLEX_WF.')
      if(complex_wf_in/=0)r1_c2_in=2
      write(io_to)complex_wf_in
     case('NONCOLL_SPIN')
      read(io_from,*,iostat=ierr)noncoll_spin_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NONCOLL_SPIN.')
      write(io_to)noncoll_spin_in
     case('INTERACTION')
      read(io_from,'(a)',iostat=ierr)interaction_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading INTERACTION.')
      interaction_in=trim(adjustl(interaction_in))
      write(io_to)interaction_in
     case('ATOM_BASIS_TYPE')
      read(io_from,'(a)',iostat=ierr)atom_basis_type_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &ATOM_BASIS_TYPE.')
      atom_basis_type_in=trim(adjustl(atom_basis_type_in))
      write(io_to)atom_basis_type_in
     case('TOTAL_TIME')
      read(io_from,*,iostat=ierr)total_time_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &TOTAL_TIME.')
      write(io_to)total_time_in
     case default ! unknown
      call errstop('A2B_CONFIG_FILE','Found unknown info '//trim(label)//'.')
     end select
    enddo
   case('DEFINE_CONFIGS')
    read(io_from,*,iostat=ierr)nconfig_in
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading number of &
     &configs.')
    write(io_to)nconfig_in
    read(io_from,*,iostat=ierr)nitem_in
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading number of &
     &items per config.')
    write(io_to)nitem_in
    allocate(config_item_in(nitem_in),stat=ialloc)
    if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Allocation problem &
     &(config_item).')
    do item=1,nitem_in
     read(io_from,'(a)',iostat=ierr)label
     if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading config &
      &content label.')
     label=adjustl(label)
     write(io_to)label
     config_item_in(item)=trim(label)
     select case(trim(label))
     case('RELE','SELE','ETOT','LOGDET','FIDET','FI_PROD_DET','LAPDET',&
      &'PROD_LAPDET','LOCAL_POTENTIAL','NLTOT','STOT','WDMC','VALJAS','LOGP',&
      &'TWIST')
      continue
     case default
      call errstop('A2B_CONFIG_FILE',"Unrecognized item '"//trim(label)//"'.")
     end select
    enddo
    read(io_from,'(a)',iostat=ierr)label
    if(ierr/=0)call errstop('A2B_CONFIG_FILE',"Problem reading 'END &
     &DEFINE_CONFIGS'.")
    label=adjustl(label)
    if(trim(label)/='END DEFINE_CONFIGS')call errstop('A2B_CONFIG_FILE',&
     &"Expected to find 'END DEFINE_CONFIGS' but didn't.")
    write(io_to)label
   case('VMC_SAVED_STATE')
    do
     read(io_from,'(a)',iostat=ierr)label
     if(ierr/=0)call errstop('A2B_CONFIG_FILE',"Problem reading label in &
      &VMC_SAVED_STATE section in config file.")
     label=adjustl(label)
     write(io_to)label
     select case(trim(label))
     case('NO_DIFFTYPES')
      read(io_from,*,iostat=ierr)no_difftypes_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &NO_DIFFTYPES.')
      write(io_to)no_difftypes_in
     case('DTVMC')
      if(no_difftypes_in<1)call errstop('A2B_CONFIG_FILE','DTVMC found &
       &before NO_DIFFTYPES in config file.')
      allocate(dtvmc_array_in(no_difftypes_in),stat=ialloc)
      if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Allocation problem &
       &(DTVMC_ARRAY_IN).')
      read(io_from,*,iostat=ierr)dtvmc_array_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading DTVMC.')
      write(io_to)dtvmc_array_in
      deallocate(dtvmc_array_in)
     case('VMC_STEPS')
      read(io_from,*,iostat=ierr)vmc_steps_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading VMC_STEPS.')
      write(io_to)vmc_steps_in
     case('NNODES')
      read(io_from,*,iostat=ierr)nnodes_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NNODES.')
      write(io_to)nnodes_in
     case('RELE')
      if(netot_in<1)call errstop('A2B_CONFIG_FILE','NETOT undefined when &
       &reading RELE VMC_SAVED_STATE block. This is a bug.')
      allocate(rele_vmc_in(3,netot_in),stat=ialloc)
      if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Allocation problem &
       &(RELE_VMC_IN).')
      do inode=1,nnodes_in
       read(io_from,*,iostat=ierr)rele_vmc_in
       if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading RELE_VMC.')
       write(io_to)rele_vmc_in
      enddo ! inode
      deallocate(rele_vmc_in)
     case('SELE')
      if(netot_in<1)call errstop('A2B_CONFIG_FILE','NETOT undefined when &
       &reading SELE in VMC_SAVED_STATE block. This is a bug.')
      allocate(sele_vmc_in(netot_in),stat=ialloc)
      if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Allocation problem &
       &(SELE_VMC_CONFIG).')
      do inode=1,nnodes_in
       read(io_from,*,iostat=ierr)sele_vmc_in
       if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading SELE.')
       write(io_to)sele_vmc_in
      enddo ! inode
      deallocate(sele_vmc_in)
     case('END VMC_SAVED_STATE')
      exit
     case default
      call errstop('A2B_CONFIG_FILE','Unrecognized information '//&
       &trim(label)//' in VMC_SAVED_STATE section of config file.')
     end select
    enddo
   case('FINAL_VMC_RESULT')
    do
     read(io_from,'(a)',iostat=ierr)label
     if(ierr/=0)call errstop('A2B_CONFIG_FILE',"Problem reading label &
      &in FINAL_VMC_RESULT section in config file.")
     label=adjustl(label)
     if(trim(label)/='END FINAL_VMC_RESUL')then
      write(io_to)label
     else
      write(io_to)'END FINAL_VMC_RESULT' ! too long for character format
     endif
     select case(trim(label))
     case('ENERGY')
      read(io_from,*,iostat=ierr)final_vmcE_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading ENERGY')
      write(io_to)final_vmcE_in
     case('ERRORBAR')
      read(io_from,*,iostat=ierr)final_vmcdE_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading ERRORBAR')
      write(io_to)final_vmcdE_in
     case('ERRORBARU')
      read(io_from,*,iostat=ierr)final_vmcdEu_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading ERRORBARU')
      write(io_to)final_vmcdEu_in
     case('VARIANCE')
      read(io_from,*,iostat=ierr)final_vmcvar_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading VARIANCE')
      write(io_to)final_vmcvar_in
     case('CORRTIME')
      read(io_from,*,iostat=ierr)final_vmctau_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading CORRTIME')
      write(io_to)final_vmctau_in
     case('TOTMOVE')
      read(io_from,*,iostat=ierr)final_vmcmove_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading TOTMOVE')
      write(io_to)final_vmcmove_in
     case('ENERGY2')
      read(io_from,*,iostat=ierr)final_vmcE2_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading ENERGY2')
      write(io_to)final_vmcE2_in
     case('ERRORBAR2')
      read(io_from,*,iostat=ierr)final_vmcdE2_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading ERRORBAR2')
      write(io_to)final_vmcdE2_in
     case('ERRORBAR2U')
      read(io_from,*,iostat=ierr)final_vmcdE2u_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading ERRORBAR2U')
      write(io_to)final_vmcdE2u_in
!    case('END FINAL_VMC_RESULT') ! too long for character length - can't change
     case('END FINAL_VMC_RESUL')
      exit
     case default
      call errstop('A2B_CONFIG_FILE','Unrecognized information '//&
       &trim(label)//' in FINAL_VMC_RESULT section of config file.')
     end select
    enddo
   case('DMC_SAVED_STATE')
    do
     read(io_from,'(a)',iostat=ierr)label
     if(ierr/=0)call errstop('A2B_CONFIG_FILE',"Problem reading label &
      &in VMC_SAVED_STATE section in config file.")
     label=adjustl(label)
     write(io_to)label
     select case(trim(label))
     case('EBEST')
      read(io_from,*,iostat=ierr)ebest_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading EBEST.')
      write(io_to)ebest_in
     case('EBEST_INIT')
      read(io_from,*,iostat=ierr)ebest_init_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading EBEST_INIT.')
      write(io_to)ebest_init_in
     case('EREF')
      read(io_from,*,iostat=ierr)eref_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading EREF.')
      write(io_to)eref_in
     case('DTEFF_EBEST_INIT')
      read(io_from,*,iostat=ierr)dteff_ebest_init_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &DTEFF_EBEST_INIT.')
      write(io_to)dteff_ebest_init_in
     case('DTEFF_BEST')
      read(io_from,*,iostat=ierr)dteff_best_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading DTEFF_BEST.')
      write(io_to)dteff_best_in
     case('DMCEQUIL_STEPS')
      read(io_from,*,iostat=ierr)dmcequil_steps_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &DMCEQUIL_STEPS.')
      write(io_to)dmcequil_steps_in
     case('DMCSTATS_STEPS')
      read(io_from,*,iostat=ierr)dmcstats_steps_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &DMCSTATS_STEPS.')
      write(io_to)dmcstats_steps_in
     case('TPDMC')
      read(io_from,*,iostat=ierr)tpdmc_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading TPDMC.')
      write(io_to)tpdmc_in
     case('LWDMC')
      read(io_from,*,iostat=ierr)lwdmc_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading LWDMC.')
      write(io_to)lwdmc_in
     case('GROWTH_ESTIMATOR')
      read(io_from,*,iostat=ierr)growth_estimator_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &GROWTH_ESTIMATOR.')
      write(io_to)growth_estimator_in
     case('NUMER_EXPECT')
      read(io_from,*,iostat=ierr)numer_expect_in(1:13)
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &NUMER_EXPECT.')
      write(io_to)numer_expect_in(1:13)
     case('DENOMINATOR_WT')
      read(io_from,*,iostat=ierr)denominator_wt_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &DENOMINATOR_WT.')
      write(io_to)denominator_wt_in
     case('LOG_PI_WT')
      read(io_from,*,iostat=ierr)log_pi_wt_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading LOG_PI_WT.')
      write(io_to)log_pi_wt_in
     case('LOG_PI_WT2')
      read(io_from,*,iostat=ierr)log_pi_wt2_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading LOG_PI_WT2.')
      write(io_to)log_pi_wt2_in
     case('LOG_PI_WT_ARRAY')
      if(tpdmc_in<0)call errstop('A2B_CONFIG_FILE','LOG_PI_WT_ARRAY found &
       &before TPDMC in config file. This is a bug.')
      if(tpdmc_in==0)call errstop('A2B_CONFIG_FILE','LOG_PI_WT_ARRAY found &
       &in config file even though the stored TPDMC is zero. This is a bug.')
      allocate(log_pi_wt_array_in(0:tpdmc_in-1),stat=ialloc)
      if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Allocation problem &
       &(LOG_PI_WT_ARRAY)')
      read(io_from,*,iostat=ierr)log_pi_wt_array_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &LOG_PI_WT_ARRAY.')
      write(io_to)log_pi_wt_array_in
      deallocate(log_pi_wt_array_in)
     case('LOG_PI_WT_ARRAY2')
      if(tpdmc_in<0)call errstop('A2B_CONFIG_FILE',&
       &'LOG_PI_WT_ARRAY2 found before TPDMC in config file. This is a bug.')
      allocate(log_pi_wt_array2_in(0:tpdmc_in),stat=ialloc)
      if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Allocation problem &
       &(LOG_PI_WT_ARRAY2)')
      read(io_from,*,iostat=ierr)log_pi_wt_array2_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &LOG_PI_WT_ARRAY2.')
      write(io_to)log_pi_wt_array2_in
      deallocate(log_pi_wt_array2_in)
     case('NUMERATOR_WT2')
      read(io_from,*,iostat=ierr)numerator_wt2_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &NUMERATOR_WT2.')
      write(io_to)numerator_wt2_in
     case('DENOMINATOR_WT2')
      read(io_from,*,iostat=ierr)denominator_wt2_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &DENOMINATOR_WT2.')
      write(io_to)denominator_wt2_in
     case('END DMC_SAVED_STATE')
      exit
     case default
      call errstop('A2B_CONFIG_FILE','Unrecognized information '//&
       &trim(label)//' in DMC_SAVED_STATE section of config file.')
     end select
    enddo
   case('REBLOCK_DATA')
    read(io_from,*,iostat=ierr)reblock_nbs_in,reblock_nstep_in,&
     &reblock_nobs_in,reblock_sum_w_in,popstats_in
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading number of &
     &reblock blocks.')
    write(io_to)reblock_nbs_in,reblock_nstep_in,&
     &reblock_nobs_in,reblock_sum_w_in,popstats_in
    allocate(reblock_block_length_in(reblock_nbs_in),&
     &reblock_s_u_f_in(reblock_nbs_in),&
     &reblock_sum_w2_closed_in(reblock_nbs_in),&
     &reblock_sum_w_open_in(reblock_nbs_in),&
     &reblock_sum_o2_closed_in(reblock_nobs_in,reblock_nbs_in),&
     &reblock_sum_o_open_in(reblock_nobs_in,reblock_nbs_in),&
     &reblock_sum_o_in(reblock_nobs_in),&
     &stat=ialloc)
    call check_alloc(ialloc,'A2B_CONFIG_FILE','reblock_block_length_config')
    read(io_from,*,iostat=ierr)reblock_block_length_in,&
     &reblock_s_u_f_in,&
     &reblock_sum_w2_closed_in,&
     &reblock_sum_w_open_in,&
     &reblock_sum_o2_closed_in,&
     &reblock_sum_o_open_in,&
     &reblock_sum_o_in
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading reblock blocks.')
    write(io_to)reblock_block_length_in,&
     &reblock_s_u_f_in,&
     &reblock_sum_w2_closed_in,&
     &reblock_sum_w_open_in,&
     &reblock_sum_o2_closed_in,&
     &reblock_sum_o_open_in,&
     &reblock_sum_o_in
    deallocate(reblock_block_length_in,reblock_s_u_f_in,&
     &reblock_sum_w2_closed_in,reblock_sum_w_open_in,reblock_sum_o2_closed_in,&
     &reblock_sum_o_open_in,reblock_sum_o_in)
    if(popstats_in/=0)then
     allocate( reblock_sum_ow0_in(reblock_nobs_in),&
      &reblock_sum_ow2_in(reblock_nobs_in),&
      &reblock_sum_o2w0_in(reblock_nobs_in),&
      &reblock_sum_o2w2_in(reblock_nobs_in),&
      &stat=ialloc)
     call check_alloc(ialloc,'A2B_CONFIG_FILE','reblock_sum_ow0_in')
     read(io_from,*,iostat=ierr)reblock_sum_w2_in,&
      &reblock_sum_ow0_in,&
      &reblock_sum_ow2_in,&
      &reblock_sum_o2w0_in,&
      &reblock_sum_o2w2_in
     if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading popstats.')
     write(io_to)reblock_sum_w2_in,&
      &reblock_sum_ow0_in,&
      &reblock_sum_ow2_in,&
      &reblock_sum_o2w0_in,&
      &reblock_sum_o2w2_in
     deallocate( reblock_sum_ow0_in,reblock_sum_ow2_in,&
      &reblock_sum_o2w0_in,reblock_sum_o2w2_in)
    endif
    read(io_from,'(a)',iostat=ierr)label
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading label END &
     &REBLOCK_DATA in config file.')
    label=adjustl(label)
    if(trim(label)/='END REBLOCK_DATA')call errstop('A2B_CONFIG_FILE',&
     &'Wrong label "'//trim(label)//'" instead of END REBLOCK_DATA in &
     &config file.')
    write(io_to)label
   case('RANDOM')
    read(io_from,'(a)',iostat=ierr)label
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading label NNODES &
     &in RANDOM section of config file.')
    label=adjustl(label)
    if(trim(label)/='NNODES')call errstop('A2B_CONFIG_FILE','Wrong label &
     &"'//trim(label)//'" instead of NNODES in RANDOM section of config file.')
    write(io_to)label
    read(io_from,*,iostat=ierr)nnodes_in
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NNODES in &
     &RANDOM section of config file.')
    write(io_to)nnodes_in
    do inode=1,nnodes_in
     read(io_from,*,iostat=ierr)random_state_in
     if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading RANDOM_STATE &
      &for node '//trim(i2s(inode))//' in RANDOM section of config file.')
     write(io_to)random_state_in
    enddo
    read(io_from,'(a)',iostat=ierr)label
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading RANDOM &
     &section of config file.')
    label=adjustl(label)
    if(trim(label)/='BUFFER'.and.trim(label)/='END RANDOM')call &
     &errstop('A2B_CONFIG_FILE','Wrong label "'//trim(label)//'" instead of&
     & BUFFER or END RANDOM in config file.')
    if(trim(label)=='BUFFER')then
     write(io_to)label
     read(io_from,*,iostat=ierr)ranx_max,ranx_gauss_max
     if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading BUFFER sizes &
      & in RANDOM section of config file.')
     write(io_to)ranx_max,ranx_gauss_max
     do inode=1,nnodes_in
      read(io_from,*,iostat=ierr)ranx_indx,ranx_gauss_indx
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading BUFFER data &
       & in RANDOM section of config file.')
      write(io_to)ranx_indx,ranx_gauss_indx
     enddo
     allocate(ranx_buffer(ranx_max),gauss_buffer(ranx_gauss_max),stat=ierr)
     if(ierr/=0)call errstop('A2B_CONFIG_FILE','Allocation problem for &
      &ranx_buffer/gauss_buffer.')
     do inode=1,nnodes_in
      read(io_from,*,iostat=ierr)ranx_buffer
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading RANX_BUFFER &
       &in RANDOM section of config file.')
      write(io_to)ranx_buffer
     enddo
     do inode=1,nnodes_in
      read(io_from,*,iostat=ierr)gauss_buffer
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading GAUSS_BUFFER &
       &in RANDOM section of config file.')
      write(io_to)gauss_buffer
     enddo
     deallocate(ranx_buffer,gauss_buffer)
     read(io_from,'(a)',iostat=ierr)label
     if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading END RANDOM &
      &in config file.')
     label=adjustl(label)
     if(trim(label)/='END RANDOM')call errstop('A2B_CONFIG_FILE',&
      &'Wrong label "'//trim(label)//'" instead of END RANDOM in config file.')
    endif
    write(io_to)label
   case('GEOMETRY')
    read(io_from,*,iostat=ierr)nitot_in
    if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NITOT.')
    write(io_to)nitot_in
    if(nitot_in>0)then
     allocate(rion_in(3,nitot_in),stat=ialloc)
     if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Allocation problem (RION_IN)')
     read(io_from,*,iostat=ierr)rion_in
     if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading RION.')
     write(io_to)rion_in
     deallocate(rion_in)
    endif ! nitot_in>0
    read(io_from,'(a)',iostat=ierr)label
    if(ierr/=0)call errstop('A2V_CONFIG_FILE',"Problem reading 'END &
     &GEOMETRY' label in config file.")
    label=adjustl(label)
    if(trim(label)/='END GEOMETRY')call errstop('A2B_CONFIG_FILE','Wrong &
     &label "'//trim(label)//'" read in when expecting "END GEOMETRY".')
    write(io_to)label
   case('CONFIGS')
    exit
   case default
    call errstop('A2B_CONFIG_FILE','Unrecognized section '//trim(label)//&
     &' in config file.')
   end select
  enddo

! Configs, only if config label found
  if(nconfig_in<0)call errstop('A2B_CONFIG_FILE','Number of configs undefined &
   &when about to read configs.')

  if(nconfig_in>0)then

! Check that configs are defined before loading them
   if(nitem_in==0)call errstop('A2B_CONFIG_FILE','Config contents undefined &
    &when about to read configs.')

! Allocate arrays.
   do item=1,nitem_in
    select case(trim(config_item_in(item)))
    case('RELE') ; allocate(rele_in(3,netot_in),stat=ialloc)
    case('SELE') ; allocate(sele_in(netot_in),stat=ialloc)
    case('LOGDET') ; allocate(logdet_in(nspin_in,ndet_in),stat=ialloc)
    case('FIDET') ; allocate(fidet_in(3,netot_in,r1_c2_in),stat=ialloc)
    case('FI_PROD_DET') ; allocate(fi_prod_det_in&
     &(3,ndet_in,nemax_in,r1_c2_in,nspin_in),stat=ialloc)
    case('PROD_LAPDET') ; allocate(prod_lapdet_in&
     &(ndet_in,nemax_in,r1_c2_in,nspin_in),stat=ialloc)
    case default ; ialloc=0
    end select
    if(ialloc/=0)call errstop('A2B_CONFIG_FILE','Allocation problem ('//&
     &trim(config_item_in(item))//').')
   enddo

! Copy configs.
   do iconfig=1,nconfig_in
    do item=1,nitem_in
     select case(trim(config_item_in(item)))
     case('RELE')
      read(io_from,*,iostat=ierr)rele_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading RELE for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)rele_in
     case('SELE')
      read(io_from,*,iostat=ierr)sele_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading SELE for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)sele_in
     case('ETOT')
      read(io_from,*,iostat=ierr)etot_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading ETOT for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)etot_in
     case('LOGDET')
      read(io_from,*,iostat=ierr)logdet_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading LOGDET for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)logdet_in
     case('FIDET')
      read(io_from,*,iostat=ierr)fidet_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading FIDET for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)fidet_in
     case('FI_PROD_DET')
      read(io_from,*,iostat=ierr)fi_prod_det_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading FI_PROD_DET &
       &for config '//trim(i2s(iconfig))//'.')
      write(io_to)fi_prod_det_in
     case('LAPDET')
      read(io_from,*,iostat=ierr)lapdet_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading LAPDET for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)lapdet_in
     case('PROD_LAPDET')
      read(io_from,*,iostat=ierr)prod_lapdet_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading PROD_LAPDET &
       &for config '//trim(i2s(iconfig))//'.')
      write(io_to)prod_lapdet_in
     case('LOCAL_POTENTIAL')
      read(io_from,*,iostat=ierr)local_potential_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading &
       &LOCAL_POTENTIAL for config '//trim(i2s(iconfig))//'.')
      write(io_to)local_potential_in
     case('NLTOT')
      read(io_from,*,iostat=ierr)nltot_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading NLTOT for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)nltot_in
     case('STOT')
      read(io_from,*,iostat=ierr)stot_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading STOT for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)stot_in
     case('WDMC')
      read(io_from,*,iostat=ierr)wdmc_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading WDMC for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)wdmc_in
     case('VALJAS')
      read(io_from,*,iostat=ierr)valjas_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading VALJAS for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)valjas_in
     case('LOGP')
      read(io_from,*,iostat=ierr)logp_in
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading LOGP for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)logp_in
     case('TWIST')
      read(io_from,*,iostat=ierr)twist_in(:)
      if(ierr/=0)call errstop('A2B_CONFIG_FILE','Problem reading TWIST for &
       &config '//trim(i2s(iconfig))//'.')
      write(io_to)twist_in(:)
     end select
    enddo ! item
   enddo ! iconfig

! Deallocate temp arrays
   do item=1,nitem_in
    select case(trim(config_item_in(item)))
    case('RELE') ; deallocate(rele_in)
    case('SELE') ; deallocate(sele_in)
    case('LOGDET') ; deallocate(logdet_in)
    case('FIDET') ; deallocate(fidet_in)
    case('FI_PROD_DET') ; deallocate(fi_prod_det_in)
    case('PROD_LAPDET') ; deallocate(prod_lapdet_in)
    end select
   enddo

  endif ! nconfig_in>0

  if(allocated(config_item_in))deallocate(config_item_in)

! Close files
  close(io_from)
  close(io_to)

 END SUBROUTINE a2b_config_file


 SUBROUTINE errstop(routine,error)
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: routine,error
 write(6,*)
 write(6,*)'ERROR : '//trim(routine)
 write(6,*)trim(error)
 write(6,*)
 stop
 END SUBROUTINE errstop


 CHARACTER(20) FUNCTION i2s(n)
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 INTEGER i,j
 CHARACTER tmp,sign
 if(n==0)then
  i2s='0' ; return
 endif
 sign=' ' ; if(n<0)sign='-'
 do i=1,len(i2s)
  i2s(i:i)=' '
 enddo
 i=abs(n)
 do j=1,len(i2s)
  if(i==0)exit
  i2s(j:j)=achar(ichar('0')+mod(i,10))
  i=i/10
 enddo
 i=1 ; j=len_trim(i2s)
 do
  if(i>=j)exit
  tmp=i2s(j:j)
  i2s(j:j)=i2s(i:i)
  i2s(i:i)=tmp
  i=i+1
  j=j-1
 enddo
 i2s=trim(sign)//i2s
 END FUNCTION i2s


END PROGRAM format_configs

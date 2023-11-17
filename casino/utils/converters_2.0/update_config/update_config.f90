PROGRAM update_config
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0),io=10
 CHARACTER(20) gen_by,label
 CHARACTER(160) con_name,pos_name
 CHARACTER(6) weight_type,chk_random,chk_stepno
 CHARACTER(5) chk_noncoll_spin,chk_dtvmc
 INTEGER i,j,nnodes,nspin,nemax,netot,ndet,tot_ncon,idum,idum2,idum3,icon,&
  &real1_complex2,inode,ispin,idet,tpdmc,mx,equilibration_moves,ie,       &
  &stats_accum_moves,dumnnodes,ierr,ialloc,ioption,no_eqvfamilies,        &
  &step_number,nextra,nitem
 LOGICAL complex_wf,noncoll_spin,lwdmc,growth_estimator,present_in,&
  &present_out,present_posin,present_posout,con_exists,pos_exists,in_deck,&
  &out_deck,ran_exists,use_mpc
 REAL(dp) dteff_ebest_init,dteff_best,denominator_wt,log_pi_wt,&
  &numerator_wt2,denominator_wt2,log_pi_wt2,tempr,tempr3(3),tempi,tempi3(3),&
  &eref,ebest,constant_energy
 COMPLEX(dp) maxlogdet,tempc,tempc3(3)
 INTEGER,ALLOCATABLE :: ncon(:),dumncon(:),nele(:),dumnele(:),sele(:,:),&
  &dumsele(:),all_isdext(:,:),sele_vmc_config(:,:)
 REAL(dp),ALLOCATABLE :: rele(:,:,:),detr(:,:,:),etot(:),etotmpc(:),stot(:),&
  &stotmpc(:),petot(:,:),nltot(:),laptot(:),fitot(:,:,:,:,:),&
  &fi_prod_det(:,:,:,:,:,:),prod_lapdet(:,:,:,:,:),wdmc(:),numer_expect(:),&
  &log_pi_wt_term_array(:),log_pi_wt_term_array2(:),&
  &fi_prod_det_new(:,:,:,:,:,:),fitot_new(:,:,:,:,:),&
  &dumrele(:,:),rele_vmc_config(:,:,:),dtvmc_config(:)
 COMPLEX(dp),ALLOCATABLE :: detc(:,:,:),logdet(:,:,:),logprod_det(:,:)

 write(6,*)'UPDATE_CONFIG'
 write(6,*)'============='
 write(6,*)

 inquire(file='config.in',exist=present_in)
 inquire(file='vmc.posin',exist=present_posin)
 inquire(file='config.out',exist=present_out)
 inquire(file='vmc.posout',exist=present_posout)
 inquire(file='.varmin_ran_save',exist=ran_exists)

 in_deck=present_in.or.present_posin
 out_deck=present_out.or.present_posout

 if(.not.(in_deck.or.out_deck))call errstop('No IN or OUT config files found &
  &for update.')

 if(in_deck.and.out_deck)then
  write(6,*)'1) config.in and/or vmc.posin'
  write(6,*)'2) config.out and/or vmc.posout'
  write(6,*)
  write(6,*)'Please select which input deck to convert:'
  read(5,*)ioption
  if(ioption<1.or.ioption>2)call errstop('Quitting.')
  select case(ioption)
  case(1) ; out_deck=.false.
  case(2) ; in_deck=.false.
  end select
 endif

 if(in_deck)then
  con_name='config.in' ; con_exists=present_in
  pos_name='vmc.posin' ; pos_exists=present_posin
 else
  con_name='config.out' ; con_exists=present_out
  pos_name='vmc.posout' ; pos_exists=present_posout
 endif

 write(6,*)'Converting '//trim(con_name)//' and/or '//trim(pos_name)//' to &
  &new format.'
 write(6,*)

! Set some things
 nnodes=-1 ; netot=-1 ; gen_by='NONE' ; tot_ncon=0

! Read config.x if present
 if(con_exists)then

! Ask user for data
  write(6,*)"Please provide the following data:"
  write(6,*)
  write(6,*)"Enter value for COMPLEX_WF (T/F) [usually 'F']:"
  read(5,*)complex_wf
  write(6,*)
  real1_complex2=1 ; if(complex_wf)real1_complex2=2
  write(6,*)"Enter value for NONCOLL_SPIN (T/F) [usually 'F']:"
  read(5,*)noncoll_spin
  write(6,*)
  write(6,*)"Using MPC (T) as main interaction, or not (F)? [usually 'F']"
  read(5,*)use_mpc
  write(6,*)
  write(6,*)"Constant energy? [value on 8th line of 2.0-style .hist file, plus"
  write(6,*)"                  ion-ion core-polarization potential if you have"
  write(6,*)"                  one, and the result multiplied by the number of"
  write(6,*)"                  primitive cells in the simulation cell]"
  read(5,*)constant_energy
  write(6,*)
  write(6,*)'Reading ',trim(con_name),'...'

! Open file
  open(io,file=trim(con_name),status='old',form='unformatted')

! Read header and check this is not a new-style config file
  read(io,iostat=ierr)label
  if(ierr==0)then
   if(trim(label)=='INFO')call errstop('This file is already in new format.')
  endif
  rewind(io)

! Read nnodes, ncon.
  read(io,iostat=ierr)nnodes
  if(ierr/=0)call errstop('Problem reading NNODES.')
  allocate(ncon(nnodes),dumncon(nnodes),stat=ialloc)
  if(ialloc/=0)call errstop('Problem allocating NCON.')
  read(io,iostat=ierr)ncon(:)
  if(ierr/=0)call errstop('Problem reading NCON.')
  if(nnodes>1)then
   if(any(ncon(2:nnodes)/=ncon(1)))gen_by='DMC'
  endif
  tot_ncon=sum(ncon)

! Read nspin, nele, ndet and find out whether sele is present and whether
! complex_wf is .true.
  read(io,iostat=ierr)nspin
  if(ierr/=0)call errstop('Problem reading NSPIN.')
  allocate(nele(nspin),dumnele(nspin),stat=ialloc)
  if(ialloc/=0)call errstop('Problem allocating NELE.')
  rewind(io) ; read(io)dumnnodes ; read(io)dumncon
  read(io,iostat=ierr)idum,idum2,nele
  if(ierr/=0)call errstop('Problem reading NELE.')
  nemax=maxval(nele) ; netot=sum(nele) ; ndet=0
  allocate(dumrele(3,netot),dumsele(netot),stat=ialloc)
  if(ialloc/=0)call errstop('Problem allocating DUMRELE.')
  dumrele=0.d0 ; dumsele=0.d0
  rewind(io) ; read(io)dumnnodes ; read(io)dumncon
  if(noncoll_spin)then
   read(io,iostat=ierr)idum,idum2,dumnele,dumrele,dumsele,ndet
   if(ierr/=0)call errstop('Problem reading NDET.')
  else
   read(io,iostat=ierr)idum,idum2,dumnele,dumrele,ndet
   if(ierr/=0)call errstop('Problem reading NDET.')
  endif
  rewind(io) ; read(io)dumnnodes ; read(io)dumncon
  deallocate(dumncon,dumrele,dumsele)

  if(complex_wf)then
   allocate(detc(nspin,ndet,tot_ncon),stat=ialloc)
  else
   allocate(detr(nspin,ndet,tot_ncon),stat=ialloc)
  endif
  if(ialloc/=0)call errstop('Cannot allocate DET arrays.')
  allocate(logdet(nspin,ndet,tot_ncon),rele(3,netot,tot_ncon),&
   &sele(netot,tot_ncon),etot(tot_ncon),etotmpc(tot_ncon),stot(tot_ncon),&
   &stotmpc(tot_ncon),petot(4,tot_ncon),nltot(tot_ncon),laptot(tot_ncon),&
   &fitot(nemax,3,real1_complex2,nspin,tot_ncon),&
   &fi_prod_det(ndet,nemax,3,real1_complex2,nspin,tot_ncon),&
   &prod_lapdet(ndet,nemax,real1_complex2,nspin,tot_ncon),stat=ialloc)
  if(ialloc/=0)call errstop('Cannot allocate config arrays.')
  rele=0.d0 ; sele=0 ; fitot=0.d0 ; fi_prod_det=0.d0 ; prod_lapdet=0.d0
  etot=0.d0 ; etotmpc=0.d0 ; stot=0.d0 ; stotmpc=0.d0 ; petot=0.d0
  nltot=0.d0 ; laptot=0.d0 ; logdet=cmplx(0.d0,0.d0,dp)

! Read all configs
  do icon=1,tot_ncon
   if(noncoll_spin)then
    if(complex_wf)then
     read(io,iostat=ierr)idum,idum2,dumnele,rele(:,:,icon),sele,idum3,&
      &detc(:,:,icon),etot(icon),etotmpc(icon),stot(icon),stotmpc(icon),&
      &petot(:,icon),nltot(icon),laptot(icon),fitot(:,:,:,:,icon),&
      &fi_prod_det(:,:,:,:,:,icon),prod_lapdet(:,:,:,:,icon)
    else
     read(io,iostat=ierr)idum,idum2,dumnele,rele(:,:,icon),sele,idum3,&
      &detr(:,:,icon),etot(icon),etotmpc(icon),stot(icon),stotmpc(icon),&
      &petot(:,icon),nltot(icon),laptot(icon),fitot(:,:,:,:,icon),&
      &fi_prod_det(:,:,:,:,:,icon),prod_lapdet(:,:,:,:,icon)
    endif
   else
    if(complex_wf)then
     read(io,iostat=ierr)idum,idum2,dumnele,rele(:,:,icon),idum3,&
      &detc(:,:,icon),etot(icon),etotmpc(icon),stot(icon),stotmpc(icon),&
      &petot(:,icon),nltot(icon),laptot(icon),fitot(:,:,:,:,icon),&
      &fi_prod_det(:,:,:,:,:,icon),prod_lapdet(:,:,:,:,icon)
    else
     read(io,iostat=ierr)idum,idum2,dumnele,rele(:,:,icon),idum3,&
      &detr(:,:,icon),etot(icon),etotmpc(icon),stot(icon),stotmpc(icon),&
      &petot(:,icon),nltot(icon),laptot(icon),fitot(:,:,:,:,icon),&
      &fi_prod_det(:,:,:,:,:,icon),prod_lapdet(:,:,:,:,icon)
    endif
   endif
   if(ierr/=0)call errstop('Problem reading config '//trim(i2s(icon))//'.')
  enddo ! icon

  if(any(fitot/=0.d0).or.any(fi_prod_det/=0.d0).or.any(prod_lapdet/=0.d0).or.&
   &any(petot/=0.d0).or.any(nltot/=0.d0).or.any(laptot/=0.d0))then
   if(trim(gen_by)=='DMC')call errstop('File appeared to be from DMC, but has &
    &nonzero energy components.')
   gen_by='VMC'
  endif

  if(.not.noncoll_spin)sele(:,:)=0

! Read additional data
  read(io,iostat=ierr)ebest,eref,stats_accum_moves
  if(ierr/=0)backspace(io)
  read(io,iostat=ierr)weight_type
  if(ierr==0.and.trim(weight_type)/='RANDOM')then
   if(trim(gen_by)=='VMC')call errstop('File appeared to be from VMC, but has &
    &a final weights section.')
   gen_by='DMC'
   select case(trim(weight_type))
   case('NEW_WG') ; lwdmc=.true.  ; growth_estimator=.true.
   case('NEW_WN') ; lwdmc=.true.  ; growth_estimator=.false.
   case('NEW_UG') ; lwdmc=.false. ; growth_estimator=.true.
   case('NEW_UN') ; lwdmc=.false. ; growth_estimator=.false.
   case default ; call errstop('Error reading weights label.')
   end select
   if(lwdmc)then
    allocate(wdmc(tot_ncon),stat=ialloc)
    if(ialloc/=0)call errstop('Problem allocating WDMC.')
    wdmc=1.d0 ; icon=0
    do inode=1,nnodes
     read(io,iostat=ierr)wdmc(icon+1:icon+ncon(inode))
     if(ierr/=0)call errstop('Problem reading WDMC.')
     icon=icon+ncon(inode)
    enddo
   endif
   allocate(numer_expect(13),stat=ialloc)
   if(ialloc/=0)call errstop('Problem allocating NUMER_EXPECT.')
   read(io,iostat=ierr)tpdmc,equilibration_moves
   if(ierr/=0)call errstop('Problem reading TPDMC.')
   read(io,iostat=ierr)dteff_ebest_init,dteff_best
   if(ierr/=0)call errstop('Problem reading DTEFF.')
   if(tpdmc>0)then
    allocate(log_pi_wt_term_array(0:tpdmc-1),log_pi_wt_term_array2(0:tpdmc),&
     &stat=ialloc)
    if(ialloc/=0)call errstop('Problem allocating LOG_PI_WT_TERM_ARRAY.')
    read(io,iostat=ierr)numer_expect(1:13),denominator_wt,log_pi_wt,&
     &log_pi_wt_term_array(0:tpdmc-1)
   else
    read(io,iostat=ierr)numer_expect(1:13),denominator_wt,log_pi_wt
   endif
   if(ierr/=0)call errstop('Problem reading NUMER_EXPECT.')
   if(growth_estimator)then
    read(io,iostat=ierr)numerator_wt2,denominator_wt2,log_pi_wt2,&
     &log_pi_wt_term_array2(0:tpdmc)
    if(ierr/=0)call errstop('Problem reading growth-estimator data.')
   else
    numerator_wt2=0.d0 ; denominator_wt2=0.d0
    log_pi_wt2=0.d0    ; log_pi_wt_term_array2=0.d0
   endif
  else
   backspace(io)
   if(trim(gen_by)=='DMC')call errstop('File appeared to be from DMC, but no &
    &weights label is present.')
   gen_by='VMC'
  endif
  read(io,iostat=ierr)weight_type
  if(ierr==0.and.trim(weight_type)/='RANDOM')then
   write(6,*)'Copying random seeds to qmc.ran_new.'
   read(io,iostat=ierr)dumnnodes
   if(ierr/=0)call errstop('Problem reading NNODES (random seed).')
   if(dumnnodes/=nnodes)call errstop("NNODES in random seed section doesn't &
    &match value at top.")
   allocate(all_isdext(25,nnodes),stat=ialloc)
   if(ialloc/=0)call errstop('Problem allocating ALL_ISDEXT.')
   read(io,iostat=ierr)((all_isdext(i,j),i=1,25),j=1,nnodes)
   if(ierr/=0)call errstop('Problem reading ALL_ISDEXT.')
  endif

  close(io)
  write(6,*)'Done.'
  write(6,*)'Converting data...'

! Transform dets into logdet
  if(complex_wf)then
   do icon=1,tot_ncon
    do idet=1,ndet
     do ispin=1,nspin
      if(detc(ispin,idet,icon)/=cmplx(0.d0,0.d0,dp))&
       &logdet(ispin,idet,icon)=log(detc(ispin,idet,icon))
     enddo
    enddo
   enddo
  else ! .not.complex_wf
   do icon=1,tot_ncon
    do idet=1,ndet
     do ispin=1,nspin
      if(detr(ispin,idet,icon)/=0.d0)&
       &logdet(ispin,idet,icon)=log(cmplx(detr(ispin,idet,icon),0.d0,dp))
     enddo
    enddo
   enddo
  endif

  if(trim(gen_by)=='VMC')then
   allocate(logprod_det(idet,icon),stat=ialloc)
   if(ialloc/=0)call errstop('Problem allocating LOGPROD_DET.')
   do icon=1,tot_ncon
    do idet=1,ndet
     logprod_det(idet,icon)=sum(logdet(:,idet,icon))
    enddo
   enddo
! Transform fi_prod_det: multiply by prod_det/maxval(prod_det)
! Transform prod_lapdet: divide by maxval(prod_det)
! Transform fitot: reindex
   allocate(fi_prod_det_new(3,ndet,nemax,real1_complex2,nspin,tot_ncon),&
    &fitot_new(3,nemax,real1_complex2,nspin,tot_ncon),stat=ialloc)
   if(ialloc/=0)call errstop('Problem allocating FI_PROD_DET_NEW')
   fi_prod_det_new=0.d0 ; fitot_new=0.d0 ; tempi3=0.d0 ; tempi=0.d0
   do icon=1,tot_ncon
    mx=maxval(maxloc(dble(logprod_det(:,icon))))
    maxlogdet=logprod_det(mx,icon)
    do ispin=1,nspin
     do ie=1,nemax
      fitot_new(1:3,ie,:,ispin,icon)=fitot(ie,1:3,:,ispin,icon)
      do idet=1,ndet
       tempr3=fi_prod_det(idet,ie,1:3,1,ispin,icon)
       if(complex_wf)tempi3=fi_prod_det(idet,ie,1:3,2,ispin,icon)
       tempc3=cmplx(tempr3,tempi3,dp)*exp(logdet(ispin,idet,icon)-maxlogdet)
       fi_prod_det_new(1:3,idet,ie,1,ispin,icon)=dble(tempc3)
       if(complex_wf)fi_prod_det_new(1:3,idet,ie,2,ispin,icon)=aimag(tempc3)
       tempr=prod_lapdet(idet,ie,1,ispin,icon)
       if(complex_wf)tempi=prod_lapdet(idet,ie,2,ispin,icon)
       tempc=cmplx(tempr,tempi,dp)/exp(maxlogdet)
       prod_lapdet(idet,ie,1,ispin,icon)=dble(tempc)
       if(complex_wf)prod_lapdet(idet,ie,2,ispin,icon)=aimag(tempc)
      enddo
     enddo
    enddo
   enddo
   deallocate(logprod_det)
  endif
  deallocate(fitot,fi_prod_det)

  write(6,*)'Done.'
  write(6,*)

 endif ! con_exists

! Read vmc.posx if present
 if(pos_exists)then

! Ask user for data if not known yet.
  if(.not.con_exists)then
   write(6,*)"Please provide the following data:"
   write(6,*)
   write(6,*)"Enter value for COMPLEX_WF (T/F) [usually 'F']:"
   read(5,*)complex_wf
   write(6,*)
   write(6,*)"Enter value for NONCOLL_SPIN (T/F) [usually 'F']:"
   read(5,*)noncoll_spin
   write(6,*)
   write(6,*)'Enter value for NSPIN [usually 2]:'
   read(5,*)nspin
   write(6,*)
   allocate(nele(nspin))
   write(6,*)'Enter number of particles per spin channel, separated by &
    &spaces:'
   read(5,*)nele
   write(6,*)
   netot=sum(nele) ; nemax=maxval(nele)
   write(6,*)'Enter number of terms in multideterminant expansion [usually 1]:'
   read(5,*)ndet
   write(6,*)
   write(6,*)"Enter number of nodes used in this calculation:"
   read(5,*)nnodes
   write(6,*)
   write(6,*)'Using MPC (T) as main interaction, or other interaction (F)?'
   read(5,*)use_mpc
   write(6,*)
  endif
  write(6,*)'Enter number of inequivalent particle families in this &
   &calculation'
  write(6,*)' [1 if electron-only, 2 if electron-hole, etc.]'
  read(5,*)no_eqvfamilies
  write(6,*)

  write(6,*)'Reading ',trim(pos_name),'...'

! Open file
  open(io,file=trim(pos_name),status='old',form='unformatted')

! Read header and check this is not a new-style config file
  read(io,iostat=ierr)chk_noncoll_spin
  if(ierr==0)then
   if(trim(chk_noncoll_spin)=='NCSPN')then
    noncoll_spin=.true.
   else
    if(noncoll_spin)call errstop('Info mismatch: NONCOLL_SPIN appears to be &
     &F in file '//trim(pos_name)//'.')
    backspace(io)
   endif
  else
   if(noncoll_spin)call errstop('Info mismatch: NONCOLL_SPIN appears to be &
    &F in file '//trim(pos_name)//'.')
   backspace(io)
  endif

  if(.not.noncoll_spin)then
   allocate(rele_vmc_config(3,netot,nnodes),stat=ialloc)
   if(ialloc/=0)call errstop('Allocation problem (RELE_VMC_CONFIG).')
   do inode=1,nnodes
    read(io,iostat=ierr)rele_vmc_config(:,:,inode)
    if(ierr/=0)call errstop('Problem reading RELE_VMC_CONFIG for node '//&
     &trim(i2s(inode))//'.')
   enddo ! inode
  else
   allocate(rele_vmc_config(3,netot,nnodes),sele_vmc_config(netot,nnodes),&
    &stat=ialloc)
   if(ialloc/=0)call errstop('Allocation problem ([R/S]ELE_VMC_CONFIG).')
   do inode=1,nnodes
    read(io,iostat=ierr)rele_vmc_config(:,:,inode),sele_vmc_config(:,inode)
    if(ierr/=0)call errstop('Problem reading [R/S]ELE_VMC_CONFIG for node '//&
     &trim(i2s(inode))//'.')
   enddo ! inode
  endif

  read(io,iostat=ierr)chk_dtvmc
  if(ierr==0)then
   if(trim(chk_dtvmc)=='DTVMC')then
    allocate(dtvmc_config(no_eqvfamilies),stat=ialloc)
    if(ialloc/=0)call errstop('Problem allocating DTVMC_CONFIG.')
    read(io,iostat=ierr)dtvmc_config
    if(ierr/=0)call errstop('Problem reading DTVMC.')
   else
    backspace(io)
   endif
  else
   backspace(io)
  endif

  step_number=-1
  read(io,iostat=ierr)chk_stepno
  if(ierr==0)then
   if(trim(chk_stepno)=='STEPNO')then
    read(io,iostat=ierr)step_number
    if(ierr/=0)call errstop('Problem reading STEP_NUMBER.')
   else
    backspace(io)
   endif
  else
   backspace(io)
  endif

  read(io,iostat=ierr)chk_random
  if(ierr==0)then
   if(trim(chk_random)=='RANDOM')then
    read(io,iostat=ierr)dumnnodes
    if(ierr/=0)call errstop('Problem reading NNODES.')
    if(nnodes/=-1.and.dumnnodes/=nnodes)call errstop('File inconsistent with &
     &previously gathered value of NNODES.')
    nnodes=dumnnodes
    if(.not.allocated(all_isdext))then
     allocate(all_isdext(25,nnodes),stat=ialloc)
     if(ialloc/=0)call errstop('Problem allocating ALL_ISDEXT.')
    endif
    read(io,iostat=ierr)((all_isdext(i,j),i=1,25),j=1,nnodes)
    if(ierr/=0)call errstop('Problem reading ALL_ISDEXT.')
   else
    backspace(io)
   endif
  else
   backspace(io)
  endif

  close(io)
  write(6,*)'Done.'
  write(6,*)

 endif ! pos_exists

! Read from .varmin_ran_save if present
 if(ran_exists)then
  write(6,*)'Reading .varmin_ran_save...'
! Open file
  open(io,file='.varmin_ran_save',status='old',form='unformatted')
  read(io,iostat=ierr)dumnnodes
  if(ierr/=0)call errstop('Problem reading NNODES.')
  if(nnodes/=-1.and.dumnnodes/=nnodes)call errstop('File inconsistent with &
   &previously gathered value of NNODES.')
  nnodes=dumnnodes
  if(.not.allocated(all_isdext))then
   allocate(all_isdext(25,nnodes),stat=ialloc)
   if(ialloc/=0)call errstop('Problem allocating ALL_ISDEXT.')
  endif
  read(io,iostat=ierr)((all_isdext(i,j),i=1,25),j=1,nnodes)
  if(ierr/=0)call errstop('Problem reading ALL_ISDEXT.')
  close(io)
  write(6,*)'Done.'
  write(6,*)
 endif ! ran_exists

 write(6,*)'Writing new-style config file...'
 write(6,*)

! Now dump to new file
 if(trim(gen_by)=='NONE')gen_by='VMC'
 open(io,file=trim(con_name)//'_new',status='replace',form='unformatted')
 call write_label(io,'INFO')
 nextra=0
 if(pos_exists)nextra=nextra+1
 if(allocated(all_isdext))nextra=nextra+1
 if(con_exists)then
  if(trim(gen_by)=='DMC')nextra=nextra+1
 endif
 call write_label(io,'NEXTRA') ; write(io)nextra
 call write_label(io,'GEN_BY') ; write(io)gen_by
 call write_label(io,'NSPIN') ; write(io)nspin
 call write_label(io,'NELE') ; write(io)nele
 call write_label(io,'NDET') ; write(io)ndet
 call write_label(io,'COMPLEX_WF') ; write(io)complex_wf
 call write_label(io,'NONCOLL_SPIN') ; write(io)noncoll_spin
 call write_label(io,'END INFO')

 if(tot_ncon>0)then
  call write_label(io,'DEFINE_CONFIGS')
  write(io)tot_ncon
  if(trim(gen_by)=='DMC')then
   nitem=4
   if(noncoll_spin)nitem=nitem+1
   if(lwdmc)nitem=nitem+1
   write(io)nitem
   call write_label(io,'RELE')
   if(noncoll_spin)call write_label(io,'SELE')
   call write_label(io,'ETOT')
   call write_label(io,'LOGDET')
   call write_label(io,'STOT')
   if(lwdmc)call write_label(io,'WDMC')
  else
   nitem=7
   if(noncoll_spin)nitem=nitem+1
   if(ndet>1)nitem=nitem+2
   write(io)nitem
   call write_label(io,'RELE')
   if(noncoll_spin)call write_label(io,'SELE')
   call write_label(io,'ETOT')
   call write_label(io,'LOGDET')
   call write_label(io,'FIDET')
   if(ndet>1)call write_label(io,'FI_PROD_DET')
   call write_label(io,'LAPDET')
   if(ndet>1)call write_label(io,'PROD_LAPDET')
   call write_label(io,'LOCAL_POTENTIAL')
   call write_label(io,'NLTOT')
  endif
  call write_label(io,'END DEFINE_CONFIGS')
 endif

 if(trim(gen_by)=='DMC')then
! DMC saved state
  call write_label(io,'DMC_SAVED_STATE')
  call write_label(io,'EBEST') ; write(io)ebest+constant_energy
  call write_label(io,'EREF') ; write(io)eref+constant_energy
  call write_label(io,'DTEFF_EBEST_INIT') ; write(io)dteff_ebest_init
  call write_label(io,'DTEFF_BEST') ; write(io)dteff_best
  call write_label(io,'DMCEQUIL_STEPS') ; write(io)equilibration_moves
  call write_label(io,'DMCSTATS_STEPS') ; write(io)stats_accum_moves
  call write_label(io,'TPDMC') ; write(io)tpdmc
  call write_label(io,'LWDMC') ; write(io)lwdmc
  call write_label(io,'GROWTH_ESTIMATOR') ; write(io)growth_estimator
  numer_expect(1)=numer_expect(1)+constant_energy*denominator_wt
  call write_label(io,'NUMER_EXPECT') ; write(io)numer_expect
  call write_label(io,'DENOMINATOR_WT') ; write(io)denominator_wt
  call write_label(io,'LOG_PI_WT') ; write(io)log_pi_wt
  if(tpdmc>0)then
   call write_label(io,'LOG_PI_WT_ARRAY') ; write(io)log_pi_wt_term_array
  endif
  if(growth_estimator)then
   call write_label(io,'NUMERATOR_WT2') ; write(io)numerator_wt2
   call write_label(io,'DENOMINATOR_WT2') ; write(io)denominator_wt2
   call write_label(io,'LOG_PI_WT2') ; write(io)log_pi_wt2
   call write_label(io,'LOG_PI_WT_ARRAY2') ; write(io)log_pi_wt_term_array2
  endif
  call write_label(io,'END DMC_SAVED_STATE')
 elseif(pos_exists)then
! VMC saved state
  call write_label(io,'VMC_SAVED_STATE')
  call write_label(io,'NO_EQVFAMILIES') ; write(io)no_eqvfamilies
  call write_label(io,'DTVMC') ; write(io)dtvmc_config
  call write_label(io,'VMC_STEPS') ; write(io)step_number
  call write_label(io,'NNODES') ; write(io)nnodes
  call write_label(io,'RELE')
  do inode=1,nnodes
   write(io)rele_vmc_config(:,:,inode)
  enddo
  if(noncoll_spin)then
   call write_label(io,'SELE')
   do inode=1,nnodes
    write(io)sele_vmc_config(:,inode)
   enddo
  endif
  call write_label(io,'END VMC_SAVED_STATE')
 endif

! Random state
 if(allocated(all_isdext))then
  call write_label(io,'RANDOM')
  call write_label(io,'NNODES') ; write(io)nnodes
  do inode=1,nnodes
   write(io)all_isdext(:,inode)
  enddo
  call write_label(io,'END RANDOM')
 endif

! Configs
 if(tot_ncon>0)then
  call write_label(io,'CONFIGS')
  if(trim(gen_by)=='DMC')then
   do icon=1,tot_ncon
    write(io)rele(:,:,icon)
    if(noncoll_spin)write(io)sele(:,icon)
    if(use_mpc)then
     write(io)etotmpc(icon)+constant_energy
    else
     write(io)etot(icon)+constant_energy
    endif
    write(io)logdet(:,:,icon)
    if(use_mpc)then
     write(io)stotmpc(icon)+constant_energy
    else
     write(io)stot(icon)+constant_energy
    endif
    if(lwdmc)write(io)wdmc(icon)
   enddo
  else ! VMC
   do icon=1,tot_ncon
    write(io)rele(:,:,icon)
    if(noncoll_spin)write(io)sele(:,icon)
    if(use_mpc)then
     write(io)etotmpc(icon)+constant_energy
    else
     write(io)etot(icon)+constant_energy
    endif
    write(io)logdet(:,:,icon)
    write(io)fitot_new(:,:,:,:,icon)
    if(ndet>1)write(io)fi_prod_det_new(:,:,:,:,:,icon)
    write(io)laptot(icon)
    if(ndet>1)write(io)prod_lapdet(:,:,:,:,icon)
    if(use_mpc)then
     write(io)petot(1,icon)+petot(2,icon)+petot(4,icon)
    else
     write(io)petot(3,icon)+petot(4,icon)
    endif
    write(io)nltot(icon)
   enddo
  endif
 endif ! tot_ncon>0

 close(io)
 write(6,*)'Done.'
 write(6,*)'Program finished.'


CONTAINS


 SUBROUTINE write_label(io,label)
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 INTEGER,INTENT(in) :: io
 CHARACTER(20) label20
 INTEGER ierr
 label20=label
 write(io,iostat=ierr)label20
 if(ierr/=0)call errstop('Problem writing label: '//trim(label)//'.')
 END SUBROUTINE


 SUBROUTINE errstop(error)
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: error
 write(6,*)
 write(6,*)'ERROR:'
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


END PROGRAM update_config

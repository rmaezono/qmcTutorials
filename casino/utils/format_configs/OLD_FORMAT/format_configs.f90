PROGRAM format_configs
!----------------------------------------------------------------------------!
! FORMAT CONFIGS                                                             !
! --------------                                                             !
! Format a binary config.in file and make a readable config.in_formatted     !
! (or vice versa)                                                            !
!                                                                            !
! 3.2001 MDT                                                                 !
!                                                                            !
! Changes                                                                    !
! -------                                                                    !
! 3.2002 NDD - Look for info about weights, pi-weights and eff timestep.     !
! 8.2002 MDT - Fixed writing pi weights if tpdmc=0.                          !
! 8.2002 MDT - Look for info about saved state of random number generator.   !
! 1.2005 AM  - Changed format of config.in_formatted so it's more readable.  !
!----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER nspin,ndet,nemax,icon,nslater,netot,nnodes,i,j,n,ialloc,io_status,    &
  &total_moves_so_far,dummy_int,tpdmc,equilibration_moves,rnodes,is,ie,idet,   &
  &islater
 INTEGER,ALLOCATABLE :: nele(:),ncon(:),all_isdext(:,:)
 DOUBLE PRECISION etot,etotmpc,etotlim,etotmpclim,petot(4),nltot,laptot,ebest, &
  &eref,wdmcin,numer_expect(1:13),denominator_wt,numerator_wt2,denominator_wt2,&
  &log_pi_wt,log_pi_wt2,dteff_ebest_init,dteff_best
 DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: log_pi_wt_term_array,            &
  &log_pi_wt_term_array2
 DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: det,rele
 DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: fitot,prod_lapdet
 DOUBLE PRECISION,DIMENSION(:,:,:,:),ALLOCATABLE :: fi_prod
 CHARACTER(6)more_to_read,dummy_char
 LOGICAL utof,posin_exists,posinf_exists,found_random

 ialloc=0
 found_random=.false.

 inquire(file='config.in',exist=posin_exists)
 inquire(file='config.in_formatted',exist=posinf_exists)

 if(posin_exists.and.posinf_exists)then
  write(6,*)
  write(6,*)'Both config.in and config.in_formatted already exist. Delete &
   &one of'
  write(6,*)'them and FORMAT_CONFIGS will convert the other.'
  write(6,*)
  stop
 endif

 if(.not.posin_exists.and..not.posinf_exists)then
  write(6,*)
  write(6,*)'FORMAT_CONFIGS requires either a config.in or a config.in_forma&
   &tted'
  write(6,*)'file to convert and it can''t find either.'
  write(6,*)
  stop
 endif

 if(posin_exists)utof=.true.
 if(posinf_exists)utof=.false.

 write(6,*)
 if(utof)then
  write(6,*)'Formatting unformatted configs'
  write(6,*)'------------------------------'
  open(unit=8,file='config.in',form='unformatted',status='old',err=10)
  open(unit=9,file='config.in_formatted',form='formatted',status='new',err=20)
 else
  write(6,*)'Unformatting formatted configs'
  write(6,*)'------------------------------'
  open(unit=8,file='config.in_formatted',form='formatted',status='old',err=30)
  open(unit=9,file='config.in',form='unformatted',status='new',err=40)
 endif

 if(utof)then
  read(8)nnodes
  read(8)
  read(8)nspin
 else
  read(8,*)
  read(8,*)nnodes
  read(8,*)
  read(8,*)(i,j=1,nnodes)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)nspin
 endif
 allocate(nele(nspin),ncon(nnodes),stat=ialloc)

 if(ialloc/=0)call errstop('FORMAT_CONFIGS','Allocation problem <1>')
 rewind(8)

 if(utof)then
  read(8)nnodes
  read(8)ncon
  read(8)nspin,nslater,nele
 else
  read(8,*)
  read(8,*)nnodes
  read(8,*)
  read(8,*)(i,j=1,nnodes)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)nspin
  read(8,*)
  read(8,*)nslater
  read(8,*)
  read(8,*)
  do is=1,nspin
   read(8,*)dummy_int,nele(is)
  enddo
 endif
 nemax=maxval(nele)
 netot=sum(nele)
 allocate(rele(3,netot),fitot(nemax,3,nspin),stat=ialloc)
 if(ialloc/=0)call errstop('FORMAT_CONFIGS','Allocation problem <2>')

 rewind(8)

 if(utof)then
  read(8)
  read(8)
  read(8)nspin,nslater,nele,rele,ndet
 else
  read(8,*)
  read(8,*)nnodes
  read(8,*)
  read(8,*)ncon
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)nspin
  read(8,*)
  read(8,*)nslater
  read(8,*)
  read(8,*)
  do is=1,nspin
   read(8,*)dummy_int,nele(is)
  enddo
  read(8,*)
  read(8,*)
  read(8,*)
  read(8,*)
  do ie=1,netot
   read(8,*)dummy_int,rele(1:3,ie)
  enddo
  read(8,*)
  read(8,*)
  read(8,*)ndet
 endif
 allocate(det(nslater,ndet),fi_prod(ndet,nemax,3,nspin), &
  &prod_lapdet(ndet,nemax,nspin),stat=ialloc)
 if(ialloc/=0)call errstop('FORMAT_CONFIGS','Allocation problem <3>')
 rewind(8)

 write(6,*)
 if(nnodes==1)then
  write(6,*)'Data produced on a single node.'
 else
  write(6,'(3a)')' Data produced by ',trim(i2s(nnodes)),' nodes.'
 endif

 if(nnodes<=0)call errstop('FORMAT_CONFIGS','That can''t be right')
 rewind(8)

 if(utof)then
  read(8)nnodes
  read(8)(ncon(i),i=1,nnodes)
 else
  read(8,*)
  read(8,*)nnodes
  read(8,*)
  read(8,*)(ncon(i),i=1,nnodes)
 endif
 if(utof)then
  write(9,*)'Number of nodes'
  write(9,'(3x,a)')trim(i2s(nnodes))
  write(9,*)'No. of configs on node'
  do i=1,nnodes
   write(9,'(3x,a)')trim(i2s(ncon(i)))
  enddo
 else
  write(9)nnodes
  write(9)(ncon(i),i=1,nnodes)
 endif
 do n=1,nnodes
  if(utof)then
   write(9,*)
   write(9,*)'Config information on node ',trim(i2s(n))
   write(9,*)'============================='
   write(9,*)
  else
   read(8,*)
   read(8,*)
   read(8,*)
   read(8,*)
  endif
  write(6,'(4a)')' Number of configs on node ',trim(i2s(n)),' = ', &
   &trim(i2s(ncon(n)))
  do icon=1,ncon(n)
   if(utof)then
    write(9,'(a,a,a)')' --------------------------- Config no. ',&
     &trim(i2s(icon)),' -------------------------------'
    write(9,*)
    read(8,err=70,end=90)nspin,nslater,nele,rele,ndet,det,etot,etotmpc,       &
     &etotlim,etotmpclim,petot(1),petot(2),petot(3),petot(4),nltot,laptot,     &
     &fitot,fi_prod,prod_lapdet
    write(9,*)'No. of spins'
    write(9,'(3x,a)')trim(i2s(nspin))
    write(9,*)'No. of Slater determinants'
    write(9,'(3x,a)')trim(i2s(nslater))
    write(9,*)'No. of electrons'
    write(9,*)' spin   no.'
    do is=1,nspin
     write(9,'(3x,3a)')trim(i2s(is)),'     ',trim(i2s(nele(is)))
    enddo
    write(9,*)
    write(9,*)'Positions of electrons'
    write(9,*)'----------------------'
    write(9,*)'el           x                      y                       z'
    do ie=1,netot
     write(9,'(1x,a,3e24.15)')trim(i2s(ie)),rele(1:3,ie)
    enddo
    write(9,*)
    write(9,*)'No. of determinants'
    write(9,'(3x,a)')trim(i2s(ndet))
    write(9,*)
    write(9,*)'Values of Determinants '
    write(9,*)'-----------------------'
    write(9,*)'nslater ndet      det        '
    do islater=1,nslater
     do idet=1,ndet
      write(9,'(3x,4a,e24.15)')trim(i2s(islater)),'     ',trim(i2s(idet)),' ',&
       &det(islater,idet)
     enddo
    enddo
    write(9,*)
    write(9,*)'Total Energies, total potentials and total Laplacians'
    write(9,*)'-----------------------------------------------------'
    write(9,'(a30,es24.16)')'Ewald local energy         : ',etot
    write(9,'(a30,es24.16)')'MPC local energy           : ',etotmpc
    write(9,'(a30,es24.16)')'Ewald local energy lim     : ',etotlim
    write(9,'(a30,es24.16)')'MPC local energy lim       : ',etotmpclim
    write(9,'(a30,es24.16)')'Short MPC pot.             : ',petot(1)
    write(9,'(a30,es24.16)')'Long MPC pot.              : ',petot(2)
    write(9,'(a30,es24.16)')'Ewald e-e pot.             : ',petot(3)
    write(9,'(a30,es24.16)')'Local ionic pot.           : ',petot(4)
    write(9,'(a30,es24.16)')'Non-local energy           : ',nltot
    write(9,'(a30,es24.16)')'Laplacians of lndet        : ',laptot
    write(9,*)
    write(9,*)'Drift vector (fitot) '
    write(9,*)'------------------'
    write(9,*)'el sp           x                     y                      z'
    do ie=1,nemax
     do is=1,nspin
      write(9,'(1x,3a,3e23.15,i3,i3)')trim(i2s(ie)),'  ',trim(i2s(is)),&
       &fitot(ie,1:3,is)
     enddo
    enddo
    write(9,*)
    write(9,*)'Product of the drift vector (fi_prod)'
    write(9,*)'-------------------------------------'
    write(9,*)'id el sp          x                    y                       z'
    do idet=1,ndet
     do ie=1,nemax
      do is=1,nspin
       write(9,'(1x,5a,3e23.15)')trim(i2s(idet)),'  ',trim(i2s(ie)),'  ',&
        &trim(i2s(is)),fi_prod(idet,ie,1:3,is)
      enddo
     enddo
    enddo
    write(9,*)
    write(9,*)'Product of Laplacian (prod_lapdet)'
    write(9,*)'----------------------------------'
    write(9,*)'id el sp      prod_lapdet'
    do idet=1,ndet
     do ie=1,nemax
      do is=1,nspin
       write(9,'(1x,5a,e23.15)')trim(i2s(idet)),'  ',trim(i2s(ie)),'  ',&
        &trim(i2s(is)),prod_lapdet(idet,ie,is)
      enddo
     enddo
    enddo
    write(9,*)
    if(icon==ncon(n))then
     write(9,'(a50,i5,a20)')'====================== End of Config info on node'&
      &,n,' ============================'
     write(9,*)
    endif
   else
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)nspin
    read(8,*)
    read(8,*)nslater
    read(8,*)
    read(8,*)
    do is=1,nspin
     read(8,*)dummy_int,nele(is)
    enddo
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)
    do ie=1,netot
     read(8,*)dummy_int,rele(1:3,ie)
    enddo
    read(8,*)
    read(8,*)
    read(8,*)ndet
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)
     do islater=1,nslater
      do idet=1,ndet
       read(8,*)dummy_int,dummy_int,det(islater,idet)
      enddo
     enddo
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,'(a30,es24.16)')dummy_char,etot
    read(8,'(a30,es24.16)')dummy_char,etotmpc
    read(8,'(a30,es24.16)')dummy_char,etotlim
    read(8,'(a30,es24.16)')dummy_char,etotmpclim
    read(8,'(a30,es24.16)')dummy_char,petot(1)
    read(8,'(a30,es24.16)')dummy_char,petot(2)
    read(8,'(a30,es24.16)')dummy_char,petot(3)
    read(8,'(a30,es24.16)')dummy_char,petot(4)
    read(8,'(a30,es24.16)')dummy_char,nltot
    read(8,'(a30,es24.16)')dummy_char,laptot
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)
     do ie=1,nemax
      do is=1,nspin
       read(8,*)dummy_int,dummy_int,fitot(ie,1:3,is)
      enddo
     enddo
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)
     do idet=1,ndet
      do ie=1,nemax
       do is=1,nspin
        read(8,*)dummy_int,dummy_int,dummy_int,fi_prod(idet,ie,1:3,is)
       enddo
      enddo
     enddo
    read(8,*)
    read(8,*)
    read(8,*)
    read(8,*)
     do idet=1,ndet
      do ie=1,nemax
       do is=1,nspin
        read(8,*)dummy_int,dummy_int,dummy_int,prod_lapdet(idet,ie,is)
       enddo
      enddo
     enddo
    read(8,*)
    if(icon==ncon(n))then
    read(8,*)
    read(8,*)
    endif
    write(9)nspin,nslater,nele,rele,ndet,det,etot,etotmpc,etotlim,etotmpclim, &
     &petot(1),petot(2),petot(3),petot(4),nltot,laptot,fitot,fi_prod,prod_lapdet
   endif
  enddo
 enddo

 deallocate(det,fi_prod,rele,fitot)

 if(utof)then
  read(8,end=97)ebest,eref,total_moves_so_far
  write(9,'(a36,es24.15)')'Best estimate of the energy       :',ebest
  write(9,'(a36,es24.15)')'Reference  energy                 :',eref
  write(9,'(a36,i24)')'Total number of moves so far      :',total_moves_so_far
  write(9,*)
 else
  read(8,'(a36,es24.15)')dummy_char,ebest
  read(8,'(a36,es24.15)')dummy_char,eref
  read(8,'(a36,i24)')dummy_char,total_moves_so_far
  read(8,*)
  write(9)ebest,eref,total_moves_so_far
 endif

 if(utof)then
  read(8,iostat=io_status)more_to_read
 else
  read(8,*,iostat=io_status)more_to_read
 endif

 write(6,*)

 if(io_status==0)then

  if(more_to_read/='NEW_WG'.and.more_to_read/='NEW_UG' &
   &.and.more_to_read/='NEW_WN'.and.more_to_read/='NEW_UN'&
   &.and.more_to_read/='RANDOM') &
   &call errstop('FORMAT_CONFIGS','Should have read "NEW_W" or "NEW_U" or &
   &"RANDOM" at end of configs file. Found something else.')

  if(utof)then
   write(9,*)more_to_read
  else
   write(9)more_to_read
  endif

  if(more_to_read/='RANDOM')then

   write(6,*)'Information about effective timestep and Pi-weights &
    &included in files.'


   if(more_to_read=='NEW_WN'.or.more_to_read=='NEW_WG')then
    write(6,*)'Weighted DMC run: configuration weights supplied.'
   if(utof)then
    write(9,*)'Config. No.            Node           Weight of Config.'
   else
    read(8,*)dummy_char
   endif
    do i=1,nnodes
     do j=1,ncon(i)
      if(utof)then
       read(8,end=102)wdmcin
       write(9,'(i6,i20,a10,es24.15)')j,i,':',wdmcin
      else
       read(8,'(i6,i20,a10,es24.15)',end=105)dummy_int,dummy_int,dummy_char,&
        &wdmcin
       write(9)wdmcin
      endif
     enddo
    enddo
   endif

   if(utof)then
    read(8,end=102)tpdmc,equilibration_moves
    write(9,*)
    write(9,'(a36,i24)')' DMC T_p parameter                 : ',tpdmc
    write(9,'(a36,i24)')' Number of equilibration moves     : ',&
     &equilibration_moves
   else
    read(8,*)
    read(8,'(a36,i24)',end=105)dummy_char,tpdmc
    read(8,'(a36,i24)',end=105)dummy_char,equilibration_moves
    write(9)tpdmc,equilibration_moves
   endif

   allocate(log_pi_wt_term_array(0:tpdmc-1),log_pi_wt_term_array2(0:tpdmc), &
    &stat=ialloc)
   if(ialloc/=0)call errstop('FORMAT_CONFIGS','Allocatiopn problem <4>.')

   if(utof)then

    read(8,end=102)dteff_ebest_init,dteff_best
    if(tpdmc>0)then
     read(8,end=102)numer_expect(1:13),denominator_wt,log_pi_wt,         &
     &log_pi_wt_term_array(0:tpdmc-1)
    else
     read(8,end=102)numer_expect(1:13),denominator_wt,log_pi_wt
    endif
    if(more_to_read=='NEW_WG'.or.more_to_read=='NEW_UG')then
     read(8,end=102)numerator_wt2,denominator_wt2,log_pi_wt2,           &
     &log_pi_wt_term_array2(0:tpdmc)
    endif
    write(9,*)
    write(9,*)'Numerators in expectation values'
    write(9,*)'--------------------------------'
    write(9,'(a36,es24.15)')'   EBEST                           : ',      &
     &numer_expect(1)
    write(9,'(a36,es24.15)')'   TI                              : ',      &
     &numer_expect(2)
    write(9,'(a36,es24.15)')'   KEI                             : ',      &
     &numer_expect(3)
    write(9,'(a36,es24.15)')'   FISQ                            : ',      &
     &numer_expect(4)
    write(9,'(a36,es24.15)')'   POTE                            : ',      &
     &numer_expect(5)
    write(9,'(a36,es24.15)')'   SHORT                           : ',      &
     &numer_expect(6)
    write(9,'(a36,es24.15)')'   LONG                            : ',      &
     &numer_expect(7)
    write(9,'(a36,es24.15)')'   NL                              : ',      &
     &numer_expect(8)
    write(9,'(a36,es24.15)')'   LOC                             : ',      &
     &numer_expect(9)
    write(9,'(a36,es24.15)')'   VCPPEI                          : ',      &
     &numer_expect(10)
    write(9,'(a36,es24.15)')'   VCPPE                           : ',      &
     &numer_expect(11)
    write(9,'(a36,es24.15)')'   VCPPEE                          : ',      &
     &numer_expect(12)
    write(9,'(a36,es24.15)')'   TAU_EFFECTIVE                   : ',      &
     &numer_expect(13)
    write(9,*)
    write(9,'(a36,es24.15)')' Denominator in expectation values : ', &
     &denominator_wt

    write(9,'(a36,es24.15)')' Logarithm of Pi-weight            : ', &
     &log_pi_wt

    if(more_to_read=='NEW_WG'.or.more_to_read=='NEW_UG')then
     write(9,'(a36,es24.15)')' Numerator for growth estimator    : ', &
      &numerator_wt2
     write(9,'(a36,es24.15)')' Denominator for growth estimator  : ', &
      &denominator_wt2
     write(9,'(a36,es24.15)')' Logarithm of Pi-weight for g. e   : ', &
      &log_pi_wt2
    endif

    write(9,'(a36,es24.15)')' Best estimate of effective tau    : ', &
     &dteff_best
    write(9,'(a36,es24.15)')' Initial effective tau times ebest : ', &
     &dteff_ebest_init
    if(tpdmc>0)then
     write(9,*)
     write(9,*)'Arrays of factors of Pi-weights (tpdmc>0)'
     write(9,'(es24.15)')log_pi_wt_term_array(0:tpdmc-1)
     write(9,*)
    endif
    if(more_to_read=='NEW_WG'.or.more_to_read=='NEW_UG')then
     write(9,*)'Arrays of factors of Pi-weights'
     write(9,'(es24.15)')log_pi_wt_term_array2(0:tpdmc)
     write(9,*)
    endif

   else ! ftou
    read(8,*)
    read(8,*)
    read(8,*)
    do i=1,13
    read(8,'(a36,es24.15)',end=105)dummy_char,numer_expect(i)
    enddo
    read(8,*)
    read(8,'(a36,es24.15)',end=105)dummy_char,denominator_wt
    read(8,'(a36,es24.15)',end=105)dummy_char,log_pi_wt
    if(more_to_read=='NEW_WG'.or.more_to_read=='NEW_UG')then
     read(8,'(a36,es24.15)',end=105)dummy_char,numerator_wt2
     read(8,'(a36,es24.15)',end=105)dummy_char,denominator_wt2
     read(8,'(a36,es24.15)',end=105)dummy_char,log_pi_wt2
    endif
    read(8,'(a36,es24.15)',end=105)dummy_char,dteff_best
    read(8,'(a36,es24.15)',end=105)dummy_char,dteff_ebest_init
    if(tpdmc>0)then
     read(8,*)
     read(8,*)
     read(8,*)log_pi_wt_term_array(0:tpdmc-1)
     read(8,*)
    endif
    if(more_to_read=='NEW_WG'.or.more_to_read=='NEW_UG')then
     read(8,*)
     read(8,*)log_pi_wt_term_array2(0:tpdmc)
     read(8,*)
    endif
    write(9)dteff_ebest_init,dteff_best
    if(tpdmc>0)then
     write(9)numer_expect(1:13),denominator_wt,log_pi_wt,         &
      &log_pi_wt_term_array(0:tpdmc-1)
    else
     write(9)numer_expect(1:13),denominator_wt,log_pi_wt
    endif
    if(more_to_read=='NEW_WG'.or.more_to_read=='NEW_UG')then
     write(9)numerator_wt2,denominator_wt2,log_pi_wt2,           &
      &log_pi_wt_term_array2(0:tpdmc)
    endif

   endif ! utof or ftou

   deallocate(log_pi_wt_term_array,log_pi_wt_term_array2)

  else ! found RANDOM

   if(utof)then
    read(8)rnodes
   else
    read(8,*)rnodes
   endif
   if(rnodes/=nnodes)call errstop('FORMAT_CONFIGS','Number of nodes in &
    &saved state of random number generator different to the declared &
    &number of nnodes in the config data. Should not happen.')
   allocate(all_isdext(25,rnodes),stat=ialloc)
   if(ialloc/=0)call errstop('FORMAT_CONFIGS','Allocation problem <4>')
   if(utof)then
    write(9,*)rnodes
   else
    write(9)rnodes
   endif
   if(utof)then
    read(8)((all_isdext(i,j),i=1,25),j=1,rnodes)
   else
    read(8,*)((all_isdext(i,j),i=1,25),j=1,rnodes)
   endif
   if(utof)then
    write(9,*)((all_isdext(i,j),i=1,25),j=1,rnodes)
   else
    write(9)((all_isdext(i,j),i=1,25),j=1,rnodes)
   endif

   write(6,*)'Written saved state of random number generator.'
   found_random=.true.

  endif

 else ! iostat/=0

  write(6,*)'No information about effective timestep/Pi-weights &
   &or saved state of random number generator included in file.'

 endif

 deallocate(nele,ncon)

 if(.not.found_random)then

  if(utof)then
   read(8,iostat=io_status)more_to_read
  else
   read(8,*,iostat=io_status)more_to_read
  endif

  if(io_status==0.and.more_to_read=='RANDOM')then

   if(utof)then
    write(9,*)more_to_read
   else
    write(9)more_to_read
   endif

   if(utof)then
    read(8)rnodes
   else
    read(8,*)rnodes
   endif
   if(rnodes/=nnodes)call errstop('FORMAT_CONFIGS','Number of nodes in &
    &saved state of random number generator different to the declared &
    &number of nnodes in the config data. Should not happen.')
   allocate(all_isdext(25,rnodes),stat=ialloc)
   if(ialloc/=0)call errstop('FORMAT_CONFIGS','Allocation problem <5>')
   if(utof)then
    write(9,*)rnodes
   else
    write(9)rnodes
   endif
   if(utof)then
    read(8)((all_isdext(i,j),i=1,25),j=1,rnodes)
   else
    read(8,*)((all_isdext(i,j),i=1,25),j=1,rnodes)
   endif
   if(utof)then
    write(9,*)((all_isdext(i,j),i=1,25),j=1,rnodes)
   else
    write(9)((all_isdext(i,j),i=1,25),j=1,rnodes)
   endif

   write(6,*)'Written saved state of random number generator.'
   found_random=.true.

  endif

 endif ! not found_random

 write(6,*)
 if(utof)then
  write(6,*)'Formatted config data written to config.in_formatted.'
 else
  write(6,*)'Unformatted config data written to config.in.'
 endif

 close(8)
 close(9)
 stop

10 call errstop('FORMAT_CONFIGS','Error opening config.in file.')
20 call errstop('FORMAT_CONFIGS','Error opening config.in_formatted file.')
30 call errstop('FORMAT_CONFIGS','Error opening config.in_formatted file <2>.')
40 call errstop('FORMAT_CONFIGS','Error opening config.in file <2>.')
70 call errstop('FORMAT_CONFIGS','Error reading config.in')
90 call errstop('FORMAT_CONFIGS','Unexpectedly reached end of file when reading&
   & config.in <2>.')
97 call errstop('FORMAT_CONFIGS','Unexpectedly reached end of file when reading&
   & config.in <3>.')
102 call errstop('FORMAT_CONFIGS','Unexpectedly reached end of file when&
   & reading config.in <4>.')
105 call errstop('FORMAT_CONFIGS','Unexpectedly reached end of file when&
   & reading config.in_formatted <4>.')

CONTAINS


 SUBROUTINE errstop(subroutine,error)
! Write out routine name and error message then stop
  IMPLICIT NONE
  CHARACTER(*),INTENT(in) :: subroutine,error
  write(6,1)subroutine,error
1 format(/1x,'ERROR : ',a,/1x,a/)
  stop
 END SUBROUTINE errstop


 CHARACTER(20) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER i,j,n
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

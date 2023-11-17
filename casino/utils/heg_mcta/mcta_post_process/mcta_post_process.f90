PROGRAM mcta_post_process
  !----------------------------------------------------!
  ! Program to post-process MCTA-QMC data for the HEG. !
  !----------------------------------------------------!
  USE heg_tools
  IMPLICIT NONE
  ! System variables.
  INTEGER nspin,nemax,netot,dimensionality
  DOUBLE PRECISION r_s,A(3,3),B(3,3),Ainv(3,3),Binv(3,3),volume,self_term,&
     &inv_netot,inv_netot_volume
  INTEGER,ALLOCATABLE :: nele(:)
  DOUBLE PRECISION,ALLOCATABLE :: pmass(:),inv_pmass(:),pcharge(:),pcharge2(:)
  ! QMC run parameters and data.
  INTEGER qmc_nproc,nstep_per_twist,ntwist,nproc_per_twist,&
     &nstep
  LOGICAL isvmc,isdmc,is_ewald
  DOUBLE PRECISION,ALLOCATABLE :: qmc_E(:),qmc_W(:),qmc_K(:),qmc_X(:),&
     &qmc_Epp(:)
  ! Misc variables.
  CHARACTER(128) line,word,runtype,char1,char2,hist_file
  INTEGER i,j,itwist,idum,istep,ispin,ierr,ialloc
  LOGICAL exists,have_fit,have_mfit
  DOUBLE PRECISION r1,r2,r3,a_ctt,k_offset(3),sum_k,sum_x
  ! MCTA-HF energy components.
  DOUBLE PRECISION hf_ke,err_hf_ke,hf_xc,err_hf_xc
  ! Fit parameters.
  DOUBLE PRECISION fit_a,fit_b,mfit_a,err_mfit_a,mfit_b,err_mfit_b
  ! Results at various stages.
  DOUBLE PRECISION naive_E,err_naive_E,final_Epp,err_final_Epp,&
     &final_E,err_final_E,final_CE,err_final_CE
  ! k-lattices.
  INTEGER n1,n2,n3,n,kindx_min,kindx_max,ikmin,ikmax
  DOUBLE PRECISION maxk,maxk_offset
  INTEGER,ALLOCATABLE :: gamma_kindx_lattice(:),kindx_lattice(:)
  DOUBLE PRECISION,ALLOCATABLE :: gamma_k_lattice(:,:),gamma_k2_lattice(:),&
     &k_lattice(:,:),k2_lattice(:)
  ! Parameters (I/O units).
  INTEGER,PARAMETER :: io1=7,io2=8
  ! Constants (local to main routine).
  DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0
  DOUBLE PRECISION,PARAMETER :: two_pi=2.d0*pi
  DOUBLE PRECISION,PARAMETER :: four_pi_thirds=4.d0*pi/3.d0
  DOUBLE PRECISION,PARAMETER :: one_third=1.d0/3.d0

  ! Open file to read system and run parameters from.
  inquire(file='out',exist=exists)
  if(.not.exists)call quit('File "out" not found.')
  open(unit=io1,file='out',status='old',iostat=ierr)
  if(ierr/=0)call quit('Problem opening out file.')

  ! Initialize parameters to defaults.
  qmc_nproc=1
  runtype=''
  hist_file=''
  is_ewald=.false.
  nspin=2
  nstep_per_twist=0
  r_s=-1.d0
  dimensionality=3
  A=0.d0

  ! Read system and run parameters.
  do

    ! Read a line.
    read(io1,'(a)',iostat=ierr)line
    if(ierr<0)call quit('End of "out" file before reaching "Setup complete." &
       &line.')
    if(ierr>0)call quit('Problem reading "out" file.')

    ! Set things depending on which line we've read.
    if(match_line(line,'Running in parallel'))then
      word=extract_word(line,5)
      read(word,*,iostat=ierr)qmc_nproc
      if(ierr/=0)call quit('Problem reading number of MPI processes from &
         &"out" file.')
      if(qmc_nproc<1)call quit('Bad number of MPI processes read from &
         &"out" file.')

    elseif(match_line(line,'RUNTYPE'))then
      runtype=extract_word(line,-1)
      select case(trim(runtype))
      case('dmc','vmc_dmc','vmc_dmc_equil','dmc_stats','dmc_equil','dmc_dmc')
        runtype='DMC'
        hist_file='dmc.hist'
        isdmc=.true.
      case('vmc')
        runtype='VMC'
        hist_file='vmc.hist'
        isvmc=.true.
      case default
        call quit('Unsupported RUNTYPE read from "out" file.')
      end select

    elseif(match_line(line,'INTERACTION'))then
      is_ewald=trim(extract_word(line,-1))=='ewald'

    elseif(match_line(line,'VMC_NSTEP'))then
      if(isvmc)then
        word=extract_word(line,-1)
        read(word,*,iostat=ierr)nstep_per_twist
        if(ierr/=0)call quit('Problem reading VMC_NSTEP from "out" file.')
        if(nstep_per_twist<qmc_nproc)call quit('Bad value of VMC_NSTEP read &
           &from "out" file.')
        if(mod(nstep_per_twist,qmc_nproc)==0)then
          nstep_per_twist=nstep_per_twist/qmc_nproc
        else
          nstep_per_twist=nstep_per_twist/qmc_nproc+1
        endif
      endif

    elseif(match_line(line,'DMC_STATS_NSTEP'))then
      if(isdmc)then
        word=extract_word(line,-1)
        read(word,*)nstep_per_twist
        if(ierr/=0)call quit('Problem reading DMC_STATS_NSTEP from "out" &
           &file.')
        if(nstep_per_twist<1)call quit('Bad value of DMC_STATS_NSTEP read &
           &from "out" file.')
      endif

    elseif(match_line(line,'Particle name Charge Mass'))then
      read(io1,'(a)',iostat=ierr)line
      if(ierr/=0)call quit('Problem reading particle info in "out" file.')
      nspin=0
      do
        read(io1,'(a)',iostat=ierr)line
        if(ierr/=0)call quit('Problem reading particle info in "out" file.')
        write(char1,*)nspin+1
        if(.not.match_line(line,trim(adjustl(char1))//':'))exit
        nspin=nspin+1
        if(trim(extract_word(line,-1))/='Fermion')call quit('Non-fermion &
           &particles declared in "out" file.')
      enddo
      allocate(nele(nspin),pmass(nspin),inv_pmass(nspin),pcharge(nspin),&
         &pcharge2(nspin),stat=ialloc)
      if(ialloc/=0)call quit('Allocation problem (nele).')
      nele=0
      ! Get back to top of table, assuming no read errors on what's already
      ! been read, and read charges and masses.
      rewind(io1)
      do
        read(io1,'(a)')line
        if(match_line(line,'Particle name Charge Mass'))exit
      enddo
      read(io1,'(a)')line
      do ispin=1,nspin
        read(io1,'(a)')line
        char1=extract_word(line,-4)
        read(char1,*,iostat=ierr)pcharge(ispin)
        if(ierr/=0)call quit('Problem reading particle charge from "out" &
           &file.')
        char1=extract_word(line,-3)
        read(char1,*,iostat=ierr)pmass(ispin)
        if(ierr/=0)call quit('Problem reading particle mass from "out" &
           &file.')
      enddo
      inv_pmass(:)=1.d0/pmass(:)
      pcharge2(:)=pcharge(:)*pcharge(:)

    elseif(match_line(line,'MD term Det Particle'))then
      read(io1,'(a)',iostat=ierr)line
      if(ierr/=0)call quit('Problem reading free-particle info in "out" file.')
      do ispin=1,nspin
        read(io1,'(a)',iostat=ierr)line
        if(ierr/=0)call quit('Problem reading free-particle info in "out" &
           &file.')
        if(ispin==1)then
          word=extract_word(line,4)
        else
          word=extract_word(line,3)
        endif
        read(word,*,iostat=ierr)nele(ispin)
        if(ierr/=0)call quit('Problem parsing free-particle info in "out" &
           &file.')
      enddo ! ispin

    elseif(match_line(line,'Dimensionality'))then
      word=extract_word(line,-1)
      read(word,*,iostat=ierr)dimensionality
      if(ierr/=0)call quit('Problem reading dimensionality from "out" file.')
      if(dimensionality<1.or.dimensionality>3)call quit('Bad dimensionality &
         &value read from "out" file.')

    elseif(match_line(line,'r_s parameter'))then
      word=extract_word(line,-1)
      read(word,*,iostat=ierr)r_s
      if(ierr/=0)call quit('Problem reading r_s from "out" file.')
      if(r_s<0.d0)call quit('Bad r_s value read from "out" file.')

    elseif(match_line(line,'A1 ='))then
      do j=1,dimensionality
        word=extract_word(line,-dimensionality-1+j)
        read(word,*,iostat=ierr)A(j,1)
        if(ierr/=0)call quit('Problem reading lattice vector from "out" file.')
      enddo ! j
      do i=2,dimensionality
        read(io1,'(a)')line
        do j=1,dimensionality
          word=extract_word(line,-dimensionality-1+j)
          read(word,*,iostat=ierr)A(j,i)
          if(ierr/=0)call quit('Problem reading lattice vector from "out" &
             &file.')
        enddo ! j
      enddo ! i

    elseif(match_line(line,'Setup complete.'))then
      exit

    endif

  enddo ! loop over lines

  ! Check we have read everything we need.
  if(len_trim(runtype)==0.or.len_trim(hist_file)==0)call quit('Did not find &
     &RUNTYPE value in "out" file.')
  if(.not.allocated(nele))call quit('Did not find particle definitions in &
     &"out" file.')
  if(all(nele==0))call quit('Did not find particle numbers in "out" file.')
  if(r_s<0.d0)call quit('Did not find r_s value in "out" file.')
  if(nstep_per_twist<1)call quit('Did not find *_NSTEP value in "out" file.')

  ! Set things depending on particle numbers.
  nemax=maxval(nele)
  netot=sum(nele)
  inv_netot=1.d0/dble(netot)

  ! Rescale simulation cell vectors.
  Ainv=transpose(inverse(dimensionality,A,volume))
  if(abs(volume)<1.d-10*maxval(abs(A))**dimensionality)call quit&
     &('Simulation cell has zero volume.')
  select case(dimensionality)
  case(3)
    a_ctt=(volume/(dble(netot)*four_pi_thirds*r_s**3))**one_third
  case(2)
    a_ctt=(volume/(dble(netot)*pi*r_s**2))**0.5d0
  case(1)
    a_ctt=volume/(dble(netot)*2.d0*r_s)
  end select
  A=A/a_ctt

  ! Invert simulation cell matrix.
  Ainv=inverse(dimensionality,A,volume)
  B=transpose(Ainv)*two_pi
  Binv=transpose(A)/two_pi
  inv_netot_volume=inv_netot/volume

  ! Get the self-interaction term.
  call compute_self_term(dimensionality,volume,A,Ainv,B,Binv,self_term,ierr)
  if(ierr/=0)call quit('Problem in compute_self_term.')

  ! Construct lattice.
  maxk=lattice_max_r(dimensionality,B,nemax)
  maxk_offset=overestimate_wsr(dimensionality,B)
  call lattice_size(dimensionality,Binv,maxk+maxk_offset,n1,n2,n3,n)
  allocate(gamma_k_lattice(3,n),gamma_k2_lattice(n),gamma_kindx_lattice(n),&
     &k_lattice(3,n),k2_lattice(n),kindx_lattice(n),stat=ialloc)
  if(ialloc/=0)call quit('Allocation problem (*_lattice).')
  call gen_lattice_unshifted(dimensionality,B,n1,n2,n3,gamma_k_lattice,&
     &gamma_k2_lattice,gamma_kindx_lattice)
  kindx_min=n
  kindx_max=0
  do ispin=1,nspin
    call locate_reoccup_range(n,nemax,gamma_k2_lattice,gamma_kindx_lattice,&
       &maxk_offset,ikmin,ikmax)
    if(kindx_min>ikmin)kindx_min=ikmin
    if(kindx_max<ikmax)kindx_max=ikmax
  enddo ! ispin

  ! Report system parameters.
  write(char1,*)netot
  write(char2,'(f16.5)')r_s
  write(6,*)
  write(6,*)trim(adjustl(char1))//'-particle gas in '//&
     &trim(lattice_name(dimensionality,A))//' lattice at r_s [a.u.] = '//&
     &trim(adjustl(char2))

  ! Go over rest of 'out' file once to get the number of k offsets.
  ntwist=0
  nproc_per_twist=0
  do
    read(io1,'(a)',iostat=ierr)line
    if(ierr<0)exit
    if(ierr>0)call quit('Problem reading "out" file.')
    if(match_line(line,'NEW K-VECTOR'))then
      ntwist=ntwist+1
      if(nproc_per_twist==0)then
        nproc_per_twist=1
      elseif(nproc_per_twist/=1)then
        call quit('Found global K_OFFSET after blocks with per-process &
           &K_OFFSETs.')
      endif
    elseif(match_line(line,'PROCESSOR NEW K-VECTOR'))then
      ntwist=ntwist+1
      if(nproc_per_twist==0)then
        nproc_per_twist=qmc_nproc
      elseif(nproc_per_twist/=qmc_nproc)then
        call quit('Found per-process K_OFFSETs after blocks with global &
           &K_OFFSET.')
      endif
    endif
  enddo ! loop over lines
  if(ntwist==0)call quit('No twist-averaged data found.')
  nstep=ntwist*nstep_per_twist
  ! Report.
  write(char1,*)ntwist
  write(6,*)trim(adjustl(char1))//' k-vector offsets used in '//&
     &trim(runtype)//' run.'
  write(6,*)

  ! Go over 'out' file again to get the k offsets, and compute Hartree-Fock
  ! energy components while at it.
  rewind(io1)
  allocate(qmc_W(nstep),qmc_E(nstep),qmc_K(nstep),qmc_X(nstep),&
     &qmc_Epp(nstep),stat=ialloc)
  if(ialloc/=0)call quit('Allocation problem (qmc_*).')
  qmc_W=0.d0
  qmc_E=0.d0
  qmc_K=0.d0
  qmc_X=0.d0
  itwist=0
  istep=0
  k_offset=0.d0
  do
    read(io1,'(a)',iostat=ierr)line
    if(ierr/=0)exit
    if(match_line(line,'NEW K-VECTOR'))then
      itwist=itwist+1
      read(io1,*,iostat=ierr)k_offset(1:dimensionality)
      if(ierr/=0)call quit('Problem reading global K_OFFSET from "out" file.')
      call shift_lattice(n,kindx_min,kindx_max,gamma_k_lattice,&
         &gamma_kindx_lattice,k_offset,k_lattice,k2_lattice,kindx_lattice)
      qmc_K(istep+1:istep+nstep_per_twist)=hf_kinetic_energy(nspin,nele,&
         &inv_pmass,inv_netot,k2_lattice,kindx_lattice)
      if(is_ewald)qmc_X(istep+1:istep+nstep_per_twist)=hf_exchange_energy&
         &(dimensionality,nspin,nele,pcharge2,inv_netot_volume,self_term,&
         &k_lattice,kindx_lattice)
      istep=istep+nstep_per_twist
    elseif(match_line(line,'PROCESSOR NEW K-VECTOR'))then
      itwist=itwist+1
      sum_k=0.d0
      if(is_ewald)sum_x=0.d0
      do i=1,nproc_per_twist
        read(io1,*,iostat=ierr)idum,k_offset(1:dimensionality)
        if(ierr/=0)call quit('Problem reading per-process K_OFFSET from &
           &"out" file.')
        call shift_lattice(n,kindx_min,kindx_max,gamma_k_lattice,&
           &gamma_kindx_lattice,k_offset,k_lattice,k2_lattice,kindx_lattice)
        sum_k=sum_k+hf_kinetic_energy(nspin,nele,inv_pmass,inv_netot,&
           &k2_lattice,kindx_lattice)
        if(is_ewald)sum_x=sum_x+hf_exchange_energy(dimensionality,nspin,nele,&
           &pcharge2,inv_netot_volume,self_term,k_lattice,kindx_lattice)
      enddo ! i
      qmc_K(istep+1:istep+nstep_per_twist)=sum_k/dble(nproc_per_twist)
      if(is_ewald)qmc_X(istep+1:istep+nstep_per_twist)=&
         &sum_x/dble(nproc_per_twist)
      istep=istep+nstep_per_twist
    endif
  enddo ! loop over lines

  ! Done reading out file.
  close(io1)

  ! Open file to read raw history from.
  inquire(file=trim(hist_file),exist=exists)
  if(.not.exists)call quit('File '//trim(hist_file)//' not found.')
  open(unit=io2,file=trim(hist_file),status='old',iostat=ierr)
  if(ierr/=0)call quit('Problem opening '//trim(hist_file)//' file.')

  ! Skip header and/or equilibration.
  do
    read(io2,'(a)',iostat=ierr)line
    if(ierr/=0)call quit('Problem skipping header and/or equilibration in "'//&
       &trim(hist_file)//'" file.')
    if(isvmc)then
      if(match_line(line,'# Raw QMC data'))exit
    elseif(isdmc)then
      if(match_line(line,'#### START STATS'))exit
    endif
  enddo ! loop over lines

  ! Read energies and weights.
  do istep=1,nstep
    read(io2,*,iostat=ierr)idum,r1,r2,r3,line
    if(ierr<0)call quit('Not enough data in "'//trim(hist_file)//'" file.')
    if(ierr>0)call quit('Problem reading "'//trim(hist_file)//'" file.')
    if(isvmc)then
      qmc_W(istep)=1.d0
      qmc_E(istep)=r1
    elseif(isdmc)then
      qmc_W(istep)=r1
      qmc_E(istep)=r3
    endif
  enddo ! loop over lines

  ! Check that file ends here.
  read(io2,*,iostat=ierr)idum,r1,r2,r3,line
  if(ierr>=0)write(6,*)'WARNING: "'//trim(hist_file)//'" file contains more &
     &data than indicated by "out" file.'

  ! Done reading history file.
  close(io2)

  ! Reblock unprocessed data and report.
  call reblock_weighted(nstep,qmc_E,qmc_W,naive_E,err_naive_E)
  write(6,*)'Reblocked unprocessed energy [a.u./part.]:'
  write(char1,*)naive_E
  write(char2,*)err_naive_E
  write(6,*)'  E_unprocessed = '//trim(adjustl(char1))//' +/- '//&
     &trim(adjustl(char2))
  write(6,*)

  ! Perform single-block fit, process energy and report.
  call fit_twist_data(nstep,qmc_W,qmc_E,qmc_K,qmc_X,is_ewald,fit_a,fit_b,ierr)
  write(6,*)'Single-block post-processing:'
  if(ierr/=0)then
    write(6,*)'  Could not fit.'
    have_fit=.false.
  else
    qmc_Epp=qmc_E+fit_a*qmc_K+fit_b*qmc_X
    call reblock_weighted(nstep,qmc_Epp,qmc_W,final_Epp,err_final_Epp)
    write(char1,*)fit_a
    write(6,*)'  a = '//trim(adjustl(char1))
    write(char1,*)fit_b
    write(6,*)'  b = '//trim(adjustl(char1))
    if(is_ewald)then
      write(6,*)'  Energy (offset by a*HF_KE + b*HF_XE) [a.u./part.]:'
    else
      write(6,*)'  Energy (offset by a*HF_KE) [a.u./part.]:'
    endif
    write(char1,*)final_Epp
    write(char2,*)err_final_Epp
    write(6,*)'    E_pp_single = '//trim(adjustl(char1))//' +/- '//&
       &trim(adjustl(char2))
    have_fit=.true.
  endif
  write(6,*)

  ! Perform multi-block fit, process energy and report.
  call multi_fit_twist_data(ntwist,nstep,qmc_W,qmc_E,qmc_K,qmc_X,is_ewald,&
     &mfit_a,err_mfit_a,mfit_b,err_mfit_b,ierr)
  write(6,*)'Multi-block post-processing:'
  if(ierr/=0)then
    write(6,*)'  Could not fit'
    have_mfit=.false.
  else
    qmc_Epp=qmc_E+mfit_a*qmc_K+mfit_b*qmc_X
    call reblock_weighted(nstep,qmc_Epp,qmc_W,final_Epp,err_final_Epp)
    write(char1,*)mfit_a
    write(char2,*)err_mfit_a
    write(6,*)'  a = '//trim(adjustl(char1))//' +/- '//trim(adjustl(char2))
    write(char1,*)mfit_b
    write(char2,*)err_mfit_b
    write(6,*)'  b = '//trim(adjustl(char1))//' +/- '//trim(adjustl(char2))
    if(is_ewald)then
      write(6,*)'  Energy (offset by a*HF_KE + b*HF_XE) [a.u./part.]:'
    else
      write(6,*)'  Energy (offset by a*HF_KE) [a.u./part.]:'
    endif
    write(char1,*)final_Epp
    write(char2,*)err_final_Epp
    write(6,*)'    E_pp_multi = '//trim(adjustl(char1))//' +/- '//&
       &trim(adjustl(char2))
    have_mfit=.true.
  endif
  write(6,*)

  ! Check that at least one of the above has succeeded.
  if(.not.have_fit.and..not.have_mfit)call quit('Cannot continue.')
  if(have_fit.and..not.have_mfit)then
    mfit_a=fit_a
    mfit_b=fit_b
  endif

  ! Fetch from stdin.
  write(6,*)'Enter the twist-averaged HF kinetic energy and error bar &
     &[a.u./part.]:'
  read(5,*,iostat=ierr)hf_ke,err_hf_ke
  if(ierr/=0)call quit('Problem reading HF kinetic energy from stdin.')
  if(is_ewald)then
    write(6,*)'Enter the twist-averaged HF exchange energy and error bar &
       &[a.u./part.]:'
    read(5,*,iostat=ierr)hf_xc,err_hf_xc
    if(ierr/=0)call quit('Problem reading HF exchange energy from stdin.')
  else
    hf_xc=0.d0
    err_hf_xc=0.d0
  endif
  write(6,*)

  ! Compute final total energy and report.
  final_E=final_Epp-mfit_a*hf_ke-mfit_b*hf_xc
  err_final_E=sqrt(err_final_Epp**2+(mfit_a*err_hf_ke)**2+&
     &(mfit_b*err_hf_xc)**2)
  write(6,*)'Final corrected total energy [a.u./part.]:'
  write(char1,*)final_E
  write(char2,*)err_final_E
  write(6,*)'  E_final = '//trim(adjustl(char1))//' +/- '//trim(adjustl(char2))
  write(6,*)
  final_CE=final_Epp+(-1.d0-mfit_a)*hf_ke+(-1.d0-mfit_b)*hf_xc
  err_final_CE=sqrt(err_final_Epp**2+((-1.d0-mfit_a)*err_hf_ke)**2+&
     &((-1.d0-mfit_b)*err_hf_xc)**2)
  write(6,*)'Final corrected correlation energy [a.u./part.]:'
  write(char1,*)final_CE
  write(char2,*)err_final_CE
  write(6,*)'  CE_final = '//trim(adjustl(char1))//' +/- '//&
     &trim(adjustl(char2))
  write(6,*)
  write(6,*)'Relative reduction in standard error of total energy from &
     &post-processing:'
  if(err_final_E>0.d0)then
    write(char1,'(f16.5)')err_naive_E/err_final_E
    write(6,*)'  dE_unprocessed / dE_final = '//trim(adjustl(char1))
  elseif(err_naive_E>0.d0)then
    write(6,*)'  dE_unprocessed = 0 < dE_final'
  else
    write(6,*)'  dE_unprocessed = dE_final = 0'
  endif
  write(6,*)



CONTAINS


  SUBROUTINE fit_twist_data(n,W,E,K,X,use_X,a,b,ierr)
    !------------------------------------------------------------!
    ! Minimize the variance of the modified local energy E' with !
    ! respect to the fit parameters a and b:                     !
    !   E' = E + a*K + b*X                                       !
    ! This reduces to the very simple problem of minimizing      !
    !   F = c0 + c1*a**2 + c2*b**2 + 2*c3*a + 2*c4*b + 2*c5*a*b  !
    ! where                                                      !
    !   c0 = <E**2> - <E>**2  (unneeded)                         !
    !   c1 = <K**2> - <K>**2                                     !
    !   c2 = <X**2> - <X>**2                                     !
    !   c3 = <E*K> - <E>*<K>                                     !
    !   c4 = <E*X> - <E>*<X>                                     !
    !   c5 = <K*X> - <K>*<X>                                     !
    !------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    LOGICAL,INTENT(in) :: use_X
    DOUBLE PRECISION,INTENT(in) :: W(n),E(n),K(n),X(n)
    INTEGER,INTENT(out) :: ierr
    DOUBLE PRECISION,INTENT(out) :: a,b
    DOUBLE PRECISION sum_W,ave_E,ave_K,ave_K2,ave_EK,ave_X,ave_X2,&
       &ave_EX,ave_KX,c1,c2,c3,c4,c5,denom

    ! Compute statistics.
    sum_W=sum(W)
    ave_E=sum(W*E)/sum_W
    ave_K=sum(W*K)/sum_W
    ave_K2=sum(W*K**2)/sum_W
    ave_EK=sum(W*E*K)/sum_W
    if(use_X)then
      ave_X=sum(W*X)/sum_W
      ave_X2=sum(W*X**2)/sum_W
      ave_EX=sum(W*E*X)/sum_W
      ave_KX=sum(W*K*X)/sum_W
    endif

    ! Solve minimization problem.
    c1=ave_K2-ave_K**2
    c3=ave_EK-ave_E*ave_K
    if(use_X)then
      c2=ave_X2-ave_X**2
      c4=ave_EX-ave_E*ave_X
      c5=ave_KX-ave_K*ave_X
      denom=c5**2-c1*c2
      if(are_equal(denom,0.d0))then
        ierr=1
        a=0.d0
        b=0.d0
      else
        ierr=0
        a=(c2*c3-c4*c5)/denom
        b=(c1*c4-c3*c5)/denom
      endif
    else
      if(are_equal(c1,0.d0))then
        ierr=1
        a=0.d0
        b=0.d0
      else
        ierr=0
        a=-c3/c1
        b=0.d0
      endif
    endif

  END SUBROUTINE fit_twist_data


  SUBROUTINE multi_fit_twist_data(ntwist,nstep,W,E,K,X,use_X,&
     &a,best_err_a,b,best_err_b,ierr)
    !--------------------------------------------------------------!
    ! Use a blocking method to find unbiased values of the a and b !
    ! fit parameters with error bars.                              !
    !--------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ntwist,nstep
    LOGICAL,INTENT(in) :: use_X
    DOUBLE PRECISION,INTENT(in) :: W(nstep),E(nstep),K(nstep),X(nstep)
    INTEGER,INTENT(out) :: ierr
    DOUBLE PRECISION,INTENT(out) :: a,best_err_a,b,best_err_b
    INTEGER nstep_per_twist,ntwist_per_block,nstep_per_block,nblock,iblock,&
       &jblock,i1,i2,n12,new_ntwist_per_block
    DOUBLE PRECISION best_err_err_a,best_err_err_b,fit_a,fit_b,sum_w,sum_w2,&
       &var_factor,ave_a,ave_a2,err_a,err_err_a,ave_b,ave_b2,err_b,err_err_b
    DOUBLE PRECISION a_vector(ntwist),b_vector(ntwist),w_vector(ntwist)
    ! Parameters.
    INTEGER,PARAMETER :: MIN_NTWIST_IN_BLOCK=4

    ! Intialize.
    a=0.d0
    best_err_a=-1.d0
    best_err_err_a=-1.d0
    b=0.d0
    best_err_b=-1.d0
    best_err_err_b=-1.d0
    nstep_per_twist=nstep/ntwist

    ! Loop over block lengths.
    ntwist_per_block=MIN_NTWIST_IN_BLOCK
    do while(ntwist_per_block<=ntwist)

      nstep_per_block=nstep_per_twist*ntwist_per_block
      nblock=int(ntwist/ntwist_per_block)
      if(mod(ntwist,ntwist_per_block)>=MIN_NTWIST_IN_BLOCK)nblock=nblock+1

      ! Compute the set of values of a and b for each block of this size.
      a_vector=0.d0
      w_vector=0.d0
      if(use_X)b_vector=0.d0
      ierr=0
      jblock=0
      do iblock=1,nblock
        jblock=jblock+1
        i1=(iblock-1)*nstep_per_block+1
        i2=min(iblock*nstep_per_block,nstep)
        n12=i2-i1+1
        call fit_twist_data(n12,W(i1),E(i1),K(i1),X(i1),use_X,fit_a,fit_b,ierr)
        if(ierr==0)then
          a_vector(jblock)=fit_a
          if(use_X)b_vector(jblock)=fit_b
          w_vector(jblock)=dble(n12)/dble(nstep_per_block)
        else
          jblock=jblock-1
        endif
      enddo ! iblock
      nblock=jblock

      if(nblock>0)then

        ! Get mean values of a and b.
        sum_w=sum(w_vector(1:nblock))
        ave_a=sum(w_vector(1:nblock)*a_vector(1:nblock))/sum_w
        if(use_X)ave_b=sum(w_vector(1:nblock)*b_vector(1:nblock))/sum_w

        if(nblock>1)then

          ! Get statistical error bars on a and b.
          sum_w2=sum(w_vector(1:nblock)**2)
          var_factor=sum_w/(sum_w-sum_w2/sum_w)
          ave_a2=sum(w_vector(1:nblock)*a_vector(1:nblock)**2)/sum_w
          err_a=sqrt((ave_a2-ave_a**2)*var_factor/sum_w)
          err_err_a=err_a*2*var_factor/sum_w
          if(use_X)then
            ave_b2=sum(w_vector(1:nblock)*b_vector(1:nblock)**2)/sum_w
            err_b=sqrt((ave_b2-ave_b**2)*var_factor/sum_w)
            err_err_b=err_b*2*var_factor/sum_w
          endif

          ! Replace current best values of a and/or b if unset, or if the errors
          ! have smaller three-sigma intervals than current best.
          if(best_err_a<0.0.or.err_a+3.d0*err_err_a<&
             &best_err_a+3.d0*best_err_err_a)then
            a=ave_a
            best_err_a=err_a
            best_err_err_a=err_err_a
          endif
          if(use_X)then
            if(best_err_b<0.0.or.err_b+3.d0*err_err_b<&
               &best_err_b+3.d0*best_err_err_b)then
              b=ave_b
              best_err_b=err_b
              best_err_err_b=err_err_b
            endif
          endif

        else ! nblock==1

          ! No errorbar, so replace best values only if unset - if the best
          ! fit comes from here, an error will be flagged, but with this one
          ! could check if the values of a and b are sensible.
          if(best_err_a<0.0)then
            a=ave_a
            best_err_a=-1.d0
            best_err_err_a=-1.d0
          endif
          if(use_X)then
            if(best_err_b<0.0)then
              b=ave_b
              best_err_b=-1.d0
              best_err_err_b=-1.d0
            endif
          endif

        endif ! more than one valid block

      endif ! any valid blocks

      ! Increase block size.
      new_ntwist_per_block=int(1.2d0*dble(ntwist_per_block))
      if(new_ntwist_per_block>ntwist)new_ntwist_per_block=ntwist
      if(new_ntwist_per_block==ntwist_per_block)new_ntwist_per_block=&
         &new_ntwist_per_block+1
      ntwist_per_block=new_ntwist_per_block

    enddo ! loop over block sizes

    ! Set error flag.
    if(best_err_a<0.d0.or.(use_X.and.best_err_b<0.d0))then
      ierr=1
    else
      ierr=0
    endif

  END SUBROUTINE multi_fit_twist_data


  SUBROUTINE quit(msg)
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: msg
    write(6,*)msg
    stop
  END SUBROUTINE quit


END PROGRAM mcta_post_process

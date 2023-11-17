PROGRAM mcta_hf
  !--------------------------------------------------------!
  ! Program to calculate the MCTA-HF energy of a system to !
  ! a given precision.                                     !
  !--------------------------------------------------------!
  USE heg_tools
  IMPLICIT NONE
  ! MPI stuff (global).
  INCLUDE "mpif.h"
  INTEGER ierror,nproc,my_proc
  LOGICAL am_master
  ! Data types.
  INTEGER,PARAMETER :: i64=selected_int_kind(15)
  ! System variables.
  INTEGER nspin,nemax,netot,dimensionality
  DOUBLE PRECISION r_s,A(3,3),B(3,3),Ainv(3,3),Binv(3,3),volume,self_term,&
     &inv_netot,inv_netot_volume
  INTEGER,ALLOCATABLE :: nele(:)
  DOUBLE PRECISION,ALLOCATABLE :: pmass(:),inv_pmass(:),pcharge(:),pcharge2(:)
  ! Misc variables.
  CHARACTER(256) line
  CHARACTER(80) char1,char2
  INTEGER ie,je,i,j,ispin,ierr,ialloc
  INTEGER(i64) ntwist
  DOUBLE PRECISION a_ctt,target_error,hf_ke,err_hf_ke,hf_xc,err_hf_xc,k,x,&
     &ntwist_dble,Avec(3)
  ! k-lattice.
  INTEGER n1,n2,n3,n,kindx_min,kindx_max,ikmin,ikmax
  DOUBLE PRECISION maxk,maxk_offset
  INTEGER,ALLOCATABLE :: gamma_kindx_lattice(:)
  DOUBLE PRECISION,ALLOCATABLE :: gamma_k_lattice(:,:),gamma_k2_lattice(:)
  ! Arrays for communication.
  DOUBLE PRECISION,ALLOCATABLE :: hf_ke_gather(:),err_hf_ke_gather(:),&
     &hf_xc_gather(:),err_hf_xc_gather(:),ntwist_dble_gather(:)
  ! Constants.
  DOUBLE PRECISION,PARAMETER :: pi=3.14159265358979324d0
  DOUBLE PRECISION,PARAMETER :: two_pi=2.d0*pi
  DOUBLE PRECISION,PARAMETER :: four_pi_thirds=4.d0*pi/3.d0
  DOUBLE PRECISION,PARAMETER :: one_third=1.d0/3.d0

  ! Initialize MPI.
  call mpi_init(ierror)
  call mpi_comm_size(mpi_comm_world,nproc,ierror)
  call mpi_comm_rank(mpi_comm_world,my_proc,ierror)
  am_master=my_proc==0
  allocate(hf_ke_gather(nproc),err_hf_ke_gather(nproc),hf_xc_gather(nproc),&
     &err_hf_xc_gather(nproc),ntwist_dble_gather(nproc),stat=ialloc)
  if(ialloc/=0)call quit('Allocation problem (*_gather).')

  ! Loop over systems.
  do

    ! Read dimensionality (or quit signal).
    if(am_master)then
      write(6,*)
      write(6,*)'Enter dimensionality (0 to exit):'
      read(5,*,iostat=ierr)dimensionality
    endif
    call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
    if(ierr<0)call quit('Done.')
    if(ierr>0)call quit('Problem reading dimensionality.')
    call mpi_bcast(dimensionality,1,mpi_integer,0,mpi_comm_world,ierror)
    if(dimensionality==0)call quit('Done.')
    if(dimensionality<1.or.dimensionality>3)call quit('Invalid &
       &dimensionality.')

    ! Read number of particles.
    if(am_master)then
      write(6,*)'Enter particle numbers (one integer per particle type in one &
         &line):'
      read(5,'(a)',iostat=ierr)line
      nspin=0
      if(ierr==0)then
        do
          char1=extract_word(line,nspin+1)
          if(len_trim(char1)==0)exit
          read(char1,*,iostat=ierr)i
          if(ierr/=0)exit
          nspin=nspin+1
        enddo ! loop over words in line
      endif ! ierr/=0
    endif ! am_master
    call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
    if(ierr/=0)call quit('Problem reading number of particle types.')
    call mpi_bcast(nspin,1,mpi_integer,0,mpi_comm_world,ierror)
    if(nspin<1)call quit('Problem reading number of particle types.')
    allocate(nele(nspin),pmass(nspin),inv_pmass(nspin),pcharge(nspin),&
       &pcharge2(nspin),stat=ialloc)
    if(ialloc/=0)call quit('Allocation problem (nele).')
    if(am_master)then
      do ispin=1,nspin
        char1=extract_word(line,ispin)
        read(char1,*)nele(ispin)
      enddo ! ispin
    endif
    call mpi_bcast(nele,nspin,mpi_integer,0,mpi_comm_world,ierror)
    if(any(nele<=0))call quit('Negative/zero number of particles.')

    ! Read particle masses.
    if(am_master)then
      write(6,*)'Enter particle masses (one number per particle type in one &
         &line) [a.u.]:'
      read(5,'(a)',iostat=ierr)line
      pmass=0.d0
      if(ierr==0)then
        do ispin=1,nspin
          char1=extract_word(line,ispin)
          read(char1,*,iostat=ierr)pmass(ispin)
          if(ierr/=0)exit
        enddo ! ispin
      endif ! ierr/=0
    endif ! am_master
    call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
    if(ierr/=0)call quit('Problem reading particle masses.')
    call mpi_bcast(pmass,nspin,mpi_double_precision,0,mpi_comm_world,ierror)
    if(any(pmass<=0.d0))call quit('Negative/zero particle masses.')
    inv_pmass(:)=1.d0/pmass(:)

    ! Read particle masses.
    if(am_master)then
      write(6,*)'Enter particle charges (one number per particle type in one &
         &line) [a.u.]:'
      read(5,'(a)',iostat=ierr)line
      pcharge=0.d0
      if(ierr==0)then
        do ispin=1,nspin
          char1=extract_word(line,ispin)
          read(char1,*,iostat=ierr)pcharge(ispin)
          if(ierr/=0)exit
        enddo ! ispin
      endif ! ierr/=0
    endif ! am_master
    call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
    if(ierr/=0)call quit('Problem reading particle charges.')
    call mpi_bcast(pcharge,nspin,mpi_double_precision,0,mpi_comm_world,ierror)
    pcharge2(:)=pcharge(:)*pcharge(:)

    ! Read density.
    if(am_master)then
      write(6,*)'Enter r_s [a.u.]:'
      read(5,*,iostat=ierr)r_s
    endif
    call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
    if(ierr/=0)call quit('Problem reading r_s.')
    call mpi_bcast(r_s,1,mpi_double_precision,0,mpi_comm_world,ierror)
    if(r_s<=0.d0)call quit('Invalid r_s value.')

    ! Read cell geometry.
    A=0.d0
    if(dimensionality>1)then
      do i=1,dimensionality
        if(am_master)then
          write(char1,*)i
          write(6,*)'Enter unscaled simulation cell vector #'//&
             &trim(adjustl(char1))//':'
          read(5,*,iostat=ierr)Avec(1:dimensionality)
        endif
        call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
        if(ierr/=0)call quit('Problem reading cell vector.')
        call mpi_bcast(Avec,dimensionality,mpi_double_precision,0,&
           &mpi_comm_world,ierror)
        A(1:dimensionality,i)=Avec(1:dimensionality)
      enddo ! i
    else
      A(1,1)=1.d0
    endif

    ! Read target error bar.
    if(am_master)then
      write(6,*)'Enter target error bar [a.u./part.]:'
      read(5,*,iostat=ierr)target_error
    endif
    call mpi_bcast(ierr,1,mpi_integer,0,mpi_comm_world,ierror)
    if(ierr/=0)call quit('Problem reading target error bar.')
    call mpi_bcast(target_error,1,mpi_double_precision,0,mpi_comm_world,ierror)
    if(target_error<=0.d0)call quit('Invalid target error value.')
    if(am_master)write(6,*)

    ! Set things depending on particle numbers.
    nemax=maxval(nele)
    netot=sum(nele)
    inv_netot=1.d0/dble(netot)

    ! Report system parameters.
    if(am_master)then
      write(char1,*)netot
      write(char2,'(f16.5)')r_s
      write(6,*)trim(adjustl(char1))//'-particle gas in '//&
         &trim(lattice_name(dimensionality,A))//' lattice at r_s [a.u.] = '//&
         &trim(adjustl(char2))
    endif

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

    ! Invert matrix.
    Ainv=inverse(dimensionality,A,volume)
    B=transpose(Ainv)*two_pi
    Binv=transpose(A)/two_pi
    inv_netot_volume=inv_netot/volume

    ! Get the self-interaction term.
    call compute_self_term(dimensionality,volume,A,Ainv,B,Binv,self_term,ierr)
    if(ierr/=0)call quit('Problem in compute_self_term.')
    if(am_master)then
      write(char1,*)self_term*2.d0
      write(6,*)'  Electron self-image energy [a.u./part.]: '//&
         &trim(adjustl(char1))
    endif

    ! Construct lattice.
    maxk=lattice_max_r(dimensionality,B,nemax)
    maxk_offset=overestimate_wsr(dimensionality,B)
    call lattice_size(dimensionality,Binv,maxk+maxk_offset,n1,n2,n3,n)
    allocate(gamma_k_lattice(3,n),gamma_k2_lattice(n),gamma_kindx_lattice(n),&
       &stat=ialloc)
    if(ialloc/=0)call quit('Allocation problem (gamma_*_lattice).')
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
    if(am_master)then
      write(6,*)'  Reciprocal lattice:'
      write(char1,*)n
      write(6,*)'    k vectors generated  : '//trim(adjustl(char1))
      write(char1,*)kindx_max
      write(6,*)'    Active k vectors     : '//trim(adjustl(char1))
      write(char1,*)kindx_max-kindx_min+1
      write(6,*)'    Reoccupable k vectors: '//trim(adjustl(char1))
    endif

    ! Sanity check: is -k occupied when k is?
    if(am_master)then
      do ie=1,nemax
        i=gamma_kindx_lattice(ie)
        do je=1,nemax
          j=gamma_kindx_lattice(je)
          if(are_equal(gamma_k_lattice(1,j),-gamma_k_lattice(1,i)).and.&
             &are_equal(gamma_k_lattice(2,j),-gamma_k_lattice(2,i)).and.&
             &are_equal(gamma_k_lattice(3,j),-gamma_k_lattice(3,i)))exit
        enddo ! je
        if(je>nemax)then
          write(char1,*)i
          write(6,*)'  WARNING: vector #',trim(adjustl(char1)),' has no &
             &reciprocal in set'
        endif
      enddo ! ie
    endif ! am_master

    ! Report HF energies at Gamma.
    if(am_master)then
      k=hf_kinetic_energy(nspin,nele,inv_pmass,inv_netot,gamma_k2_lattice,&
         &gamma_kindx_lattice)
      x=hf_exchange_energy(dimensionality,nspin,nele,pcharge2,inv_netot_volume,&
         &self_term,gamma_k_lattice,gamma_kindx_lattice)
      write(6,*)'  HF energy at Gamma point [a.u./part.]:'
      write(char1,*)k
      write(6,*)'    K = ',trim(adjustl(char1))
      write(char1,*)x
      write(6,*)'    X = ',trim(adjustl(char1))
      write(char1,*)k+x
      write(6,*)'    E = ',trim(adjustl(char1))
    endif ! am_master

    ! Wait here for master above.
    call mpi_barrier(mpi_comm_world,ierror)

    ! Compute MCTA-HF energy.
    if(am_master)write(6,*)'  MCTA-HF energy [a.u./part.]:'
    call twist_averaged_hf_energy(dimensionality,nspin,nele,inv_pmass,&
       &pcharge2,self_term,inv_netot,inv_netot_volume,B,Binv,&
       &n,kindx_min,kindx_max,gamma_k_lattice,gamma_kindx_lattice,&
       &target_error*sqrt(dble(nproc)),hf_ke,err_hf_ke,hf_xc,err_hf_xc,&
       &ntwist,ierr)
    if(ierr/=0)call quit('Problem in twist_averaged_hf_energy.')

    ! Wait here for sync.
    call mpi_barrier(mpi_comm_world,ierror)

    ! Gather all data on master.
    ntwist_dble=dble(ntwist)
    call mpi_gather(ntwist_dble,1,mpi_double_precision,ntwist_dble_gather,1,&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    call mpi_gather(hf_ke,1,mpi_double_precision,hf_ke_gather,1,&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    call mpi_gather(err_hf_ke,1,mpi_double_precision,err_hf_ke_gather,1,&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    call mpi_gather(hf_xc,1,mpi_double_precision,hf_xc_gather,1,&
       &mpi_double_precision,0,mpi_comm_world,ierror)
    call mpi_gather(err_hf_xc,1,mpi_double_precision,err_hf_xc_gather,1,&
       &mpi_double_precision,0,mpi_comm_world,ierror)

    ! Process the data and report.
    if(am_master)then
      ntwist_dble=sum(ntwist_dble_gather)
      hf_ke=sum(hf_ke_gather*ntwist_dble_gather)/ntwist_dble
      err_hf_ke=sqrt(sum(ntwist_dble_gather*&
         &((ntwist_dble_gather-1.d0)*err_hf_ke_gather**2 + hf_ke_gather**2))/&
         &(ntwist_dble*(ntwist_dble-1.d0)) - (hf_ke**2)/(ntwist_dble-1.d0))
      hf_xc=sum(hf_xc_gather*ntwist_dble_gather)/ntwist_dble
      err_hf_xc=sqrt(sum(ntwist_dble_gather*&
         &((ntwist_dble_gather-1.d0)*err_hf_xc_gather**2 + hf_xc_gather**2))/&
         &(ntwist_dble*(ntwist_dble-1.d0)) - (hf_xc**2)/(ntwist_dble-1.d0))
      write(char1,*)hf_ke ; write(char2,*)err_hf_ke
      write(6,*)'    K = ',trim(adjustl(char1)),' +/- ',trim(adjustl(char2))
      write(char1,*)hf_xc ; write(char2,*)err_hf_xc
      write(6,*)'    X = ',trim(adjustl(char1)),' +/- ',trim(adjustl(char2))
      write(char1,*)hf_ke+hf_xc ; write(char2,*)sqrt(err_hf_ke**2+err_hf_xc**2)
      write(6,*)'    E = ',trim(adjustl(char1)),' +/- ',trim(adjustl(char2))
      write(char1,*)nint(ntwist_dble,i64)
      write(6,*)'    Info: used a total of ',trim(adjustl(char1)),' twists.'
      write(6,*)
    endif ! am_master

    ! Clean up for next system.
    deallocate(gamma_k_lattice,gamma_k2_lattice,gamma_kindx_lattice)
    deallocate(nele,pmass,inv_pmass,pcharge,pcharge2)

  enddo ! loop over systems

  ! NB, we never get here.

  ! Clean up a bit.
  deallocate(hf_ke_gather,err_hf_ke_gather,hf_xc_gather,err_hf_xc_gather,&
     &ntwist_dble_gather)

  ! Finish MPI.
  call mpi_finalize(ierror)


CONTAINS


  SUBROUTINE twist_averaged_hf_energy(dimensionality,nspin,nele,inv_pmass,&
     &pcharge2,self_term,inv_netot,inv_netot_volume,B,Binv,&
     &n,kindx_min,kindx_max,gamma_k_lattice,gamma_kindx_lattice,&
     &target_error,hf_ke,err_hf_ke,hf_xc,err_hf_xc,ntwist,ierr)
    !------------------------------------------------------------!
    ! Compute the Monte Carlo twist-averaged Hartree-Fock energy !
    ! components for the curren system to the given target       !
    ! accuracy.                                                  !
    !------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: dimensionality,nspin,nele(nspin),n,kindx_min,&
       &kindx_max,gamma_kindx_lattice(n)
    DOUBLE PRECISION,INTENT(in) :: inv_pmass(nspin),pcharge2(nspin),self_term,&
       &inv_netot,inv_netot_volume,B(3,3),Binv(3,3),gamma_k_lattice(3,n),&
       &target_error
    INTEGER,INTENT(out) :: ierr
    DOUBLE PRECISION,INTENT(out) :: hf_ke,err_hf_ke,hf_xc,err_hf_xc
    INTEGER(i64),INTENT(out) :: ntwist
    INTEGER icheck
    DOUBLE PRECISION target_error2,xivec(3),k_offset(3),k,x,sum_k,sum_k2,&
       &sum_x,sum_x2,hf_ke2,err2_hf_ke,hf_xc2,err2_hf_xc,&
       &inv_ntwist,inv_ntwist_1
    DOUBLE PRECISION k_lattice(3,n),k2_lattice(n)
    INTEGER kindx_lattice(n)
    ! Number of steps between checks for convergence.
    INTEGER,PARAMETER :: nstep_check=1000

    ierr=0

    ! Initialize output variables.
    ntwist=0_i64
    hf_ke=0.d0
    err_hf_ke=0.d0
    hf_xc=0.d0
    err_hf_xc=0.d0
    target_error2=target_error*target_error

    ! Initialize accumulators.
    sum_k=0.d0
    sum_k2=0.d0
    sum_x=0.d0
    sum_x2=0.d0

    ! Loop over random twists.
    icheck=0
    do
      ntwist=ntwist+1_i64
      icheck=icheck+1

      ! Generate random offset (=twist)
      call random_number(xivec(1:dimensionality))
      do i=1,dimensionality
        k_offset(i)=sum(xivec(1:dimensionality)*B(i,1:dimensionality))
      enddo ! i
      call min_image(dimensionality,B,Binv,k_offset,ierr)
      if(ierr/=0)return

      ! Shift the lattice at Gamma by k_offset.
      call shift_lattice(n,kindx_min,kindx_max,gamma_k_lattice,&
         &gamma_kindx_lattice,k_offset,k_lattice,k2_lattice,kindx_lattice)

      ! Compute HF kinetic energy.
      k=hf_kinetic_energy(nspin,nele,inv_pmass,inv_netot,k2_lattice,&
         &kindx_lattice)
      sum_k=sum_k+k
      sum_k2=sum_k2+k*k

      ! Compute HF exchange energy.
      x=hf_exchange_energy(dimensionality,nspin,nele,pcharge2,&
         &inv_netot_volume,self_term,k_lattice,kindx_lattice)
      sum_x=sum_x+x
      sum_x2=sum_x2+x*x

      ! Check for convergence.
      if(icheck>=nstep_check)then

        inv_ntwist=1.d0/dble(ntwist)
        inv_ntwist_1=1.d0/dble(ntwist-1)
        hf_ke=sum_k*inv_ntwist
        hf_ke2=sum_k2*inv_ntwist
        err2_hf_ke=(hf_ke2-hf_ke*hf_ke)*inv_ntwist_1
        hf_xc=sum_x*inv_ntwist
        hf_xc2=sum_x2*inv_ntwist
        err2_hf_xc=(hf_xc2-hf_xc*hf_xc)*inv_ntwist_1
        if(err2_hf_ke+err2_hf_xc<target_error2)then
          err_hf_ke=sqrt(max(err2_hf_ke,0.d0))
          err_hf_xc=sqrt(max(err2_hf_xc,0.d0))
          exit
        endif ! done
        icheck=0

      endif ! time to check for convergence

    enddo ! loop over random twists

  END SUBROUTINE twist_averaged_hf_energy


  SUBROUTINE quit(msg)
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: msg
    if(am_master)write(6,*)msg
    call mpi_finalize(ierror)
    stop
  END SUBROUTINE quit


END PROGRAM mcta_hf

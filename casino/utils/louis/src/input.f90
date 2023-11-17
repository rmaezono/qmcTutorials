 MODULE input
  USE comms
  USE esdf
  USE run_control, ONLY : errstop
  USE parallel,ONLY : am_master,nnodes
  USE store
  IMPLICIT NONE
  PRIVATE
  PUBLIC set_input_parameters,check_input_parameters
  INTEGER i,k,ij,ik,il,nlines
  REAL(dp) t1


 CONTAINS


  SUBROUTINE set_input_parameters
  IMPLICIT NONE
 
  calc_type=esdf_string('calc_type','trajectory')
  dim=esdf_integer('dimensionality',2)
  int_algorithm=esdf_string('int_algorithm','runge-kutta')
  timing_info=esdf_boolean('timing_info',.false.)
  verbose=esdf_boolean('verbose',.false.)
  vel_type=esdf_string('vel_type','deBB')
  curlweight=esdf_double('curlweight',1.d0)
  maxstp=esdf_integer('maxstp',1000000)
  init_eps=esdf_double('init_eps',0.0001d0)
  init_h=esdf_double('init_h',0.0001d0)
  converge_maxdiff=esdf_double('converge_maxdiff',0.0001d0)
  time_direction=esdf_string('time_direction','forward')
  ntrajectories=esdf_integer('ntrajectories',1)
  if(trim(adjustl(time_direction))=='forward'.or.trim(adjustl(time_direction))&
   &=='forwards')then
   traj_time_start=esdf_double('traj_time_start',0.d0)
   traj_time_end=esdf_double('traj_time_end',-999.d0)
  else
   traj_time_start=esdf_double('traj_time_start',-999.d0)
   traj_time_end=esdf_double('traj_time_end',0.d0)
  endif
  plot_all_traj=esdf_boolean('plot_all_traj',.false.)
  dtsave=esdf_double('dtsave',0.d0)
  phase_noise=esdf_integer('phase_noise',0)
  den_ntimes=esdf_integer('den_ntimes',1)
  dentype=esdf_integer('dentype',0)
  cell_x=esdf_double('cell_x',pi)
  if(dim>1)then
   cell_y=esdf_double('cell_y',pi)
  else
   cell_y=esdf_double('cell_y',0.d0)
  endif
  if(dim==3)then
   cell_z=esdf_double('cell_z',pi)
  else
   cell_z=esdf_double('cell_z',0.d0)
  endif
  nlattice=esdf_integer('nlattice',1024)
  nscgrain=esdf_integer('nscgrain',64)
  nsmoothstep=esdf_integer('nsmoothstep',8)
  den_time_start=esdf_double('den_time_start',0.d0)
  den_time_end=esdf_double('den_time_end',-999.d0)
  read_backtracked=esdf_boolean('read_backtracked',.false.)
  save_backtracked=esdf_boolean('save_backtracked',.true.)
  plot_raw=esdf_boolean('plot_raw',.true.)
  plot_cg=esdf_boolean('plot_cg',.true.)
  plot_smooth=esdf_boolean('plot_smooth',.true.)
  hfunction=esdf_boolean('hfunction',.true.)
  plot_h_integrand=esdf_boolean('plot_h_integrand',.false.)
  writefail=esdf_boolean('writefail',.false.)
  nstep_histogram=esdf_boolean('nstep_histogram',.false.)
  if(dim==3)then
   fastmode=esdf_boolean('fastmode',.true.)
  else
   fastmode=esdf_boolean('fastmode',.false.)
  endif
  wfn_type=esdf_string('wfn_type','sine_wave')
  phase_format=esdf_string('phase_format','input')
  phase_preset=esdf_integer('phase_preset',0)
  weight_format=esdf_string('weight_format','input')
  weight_preset=esdf_integer('weight_preset',0)
  negphase=esdf_boolean('negphase',.false.)
  transposephase=esdf_boolean('transposephase',.false.)
  testrun=esdf_boolean('testrun',.false.)
  nmodes=esdf_integer('nmodes',0)

! Block input

  if(esdf_block('trajectory_start',nlines))then
   if(nlines/=ntrajectories)call errstop('SET_INPUT_PARAMETERS','There should&
    & be NTRAJECTORIES lines in the TRAJECTORY_START block in input.')
   allocate(xstart(ntrajectories),ystart(ntrajectories),zstart(ntrajectories),&
    &stat=ialloc)
   if(ialloc/=0)call errstop('SET_INPUT_PARAMETERS','Unable to allocate &
    &vectors for trajectory starting positions.')
   select case(dim)
    case(1)
     do i=1,ntrajectories
      read(block_data(i),*)xstart(i)
     enddo
     ystart(:)=0.d0 ; zstart(:)=0.d0
    case(2)
     do i=1,ntrajectories
      read(block_data(i),*)xstart(i),ystart(i)
     enddo
     zstart(:)=0.d0
    case(3) 
     do i=1,ntrajectories
      read(block_data(i),*)xstart(i),ystart(i),zstart(i)
     enddo
   end select
  endif

  if(esdf_block('coarse_graining_lengths',nlines))then
   if(nlines/=2)call errstop('SET_INPUT_PARAMETERS','There should be two&
    & lines in the COARSE_GRAINING_LENGTHS block in input.')
   read(block_data(1),*)num_ncgrain
   if(num_ncgrain<1.or.num_ncgrain>10)call errstop('SET_INPUT_PARAMETERS',&
    &'Invalid value for NUM_NCGRAIN on first line of COARSE_GRAINING_LENGTHS&
    & block in input. Must be between 1 and 10.')
   if(num_ncgrain/=1.and.dim==3.and.trim(adjustl(calc_type))=='density')call&
    & errstop('SET_INPUT_PARAMETERS','Because of the way the algorithm works,&
    & only one coarse-graining length may currently be defined in input for&
    & 3-dimensional density calculations.')
   allocate(ncgrain(num_ncgrain),stat=ialloc)
   if(ialloc/=0)call errstop('SET_INPUT_PARAMETERS','Allocation error for&
    & coarse-graining lengths.')
   read(block_data(2),*,err=1,end=1)ncgrain(1:num_ncgrain)
  endif

  if(esdf_block('3dplanes',nlines))then
   if(nlines/=2)call errstop('SET_INPUT_PARAMETERS','There should be two&
    & lines in the 3DPLANES block in input.')
   read(block_data(1),*)n3dplanes
   if(n3dplanes<1.or.n3dplanes>10)call errstop('SET_INPUT_PARAMETERS',&
    &'Invalid value for N3DPLANES on first line of 3DPLANES&
    & block in input. Must be between 1 and 10.')
   allocate(plane_number(n3dplanes),stat=ialloc)
   if(ialloc/=0)call errstop('SET_INPUT_PARAMETERS','Allocation error for&
    & PLANE_NUMBER.')
   read(block_data(2),*,err=2,end=2)plane_number(1:n3dplanes)
  else
   allocate(plane_number(1),stat=ialloc)
   if(ialloc/=0)call errstop('SET_INPUT_PARAMETERS','Allocation error for&
    & PLANE_NUMBER <2>.')
   n3dplanes=1
   plane_number(1)=nlattice/2
  endif

! Set phases and/or weights depending on choice of wavefunction.
  select case(trim(adjustl(wfn_type)))
   case('sine_wave','sin_wave','sine_wfn','sin_wfn')
    wfn_type='sine_wave'
    call set_phases
   case('scaled_sine','scaled_sin','scaled_sin_wfn','scaled_sine_wfn')
    wfn_type='scaled_sine'
    call set_phases
    call set_weights
  end select

  return

1 call errstop('SET_INPUT_PARAMETERS','Error reading COARSE_GRAINING_LENGTHS&
   & block in input.')
2 call errstop('SET_INPUT_PARAMETERS','Error reading 3DPLANES block in input.')

  END SUBROUTINE set_input_parameters


  SUBROUTINE set_phases
  USE random_numbers,ONLY : default_seed,randomseed,ranx_twopi,initialize_random
  IMPLICIT NONE
  INTEGER randomseed_in,idum1,idum2,ierr

  ierr=0
  if(nmodes<=0)call errstop('SET_PHASES','Number of wave function modes in &
   &superposition must be set in input through keyword NMODES.')
  
  if(phase_preset<0)call errstop('SET_PHASES',&
   &'PHASE_PRESET must be a positive integer or zero.')

! Set random seed.
  select case (trim(adjustl(phase_format)))
   case('input','default')
    randomseed=default_seed
   case('random')
    if(am_master)call system_clock(randomseed,idum1,idum2)
    call mpi_bcast(randomseed,1,mpi_integer,0,mpi_comm_world,ierr)
   case('preset')
    randomseed_in=default_seed
    call get_preset_randomseed
   case default
    call errstop('SET_PHASES','Unknown value for PHASE_FORMAT in input.')
  end select

  if(trim(adjustl(phase_format))/='input')call initialize_random

! Define the phases.
  select case(dim)

   case(1)
    nmodes_dim=nmodes
    allocate(theta_1d(nmodes),stat=ialloc)
    if(ialloc/=0)call errstop('SET_PHASES','Allocation &
     &problem with 1D phase array.')

    select case (trim(adjustl(phase_format)))  

     case('input')
      if(esdf_block('phases',nlines))then
       if(nlines/=nmodes)call errstop('SET_PHASES','Should be &
        &nmodes lines in the phases block in input.')
       do i=1,nmodes_dim
        read(block_data(i),*)theta_1d(i)
       enddo
      endif

     case('random','default','preset')
 
      do i=1,nmodes_dim
       theta_1d(i)=ranx_twopi()
      enddo

    end select

   case(2)
    nmodes_dim=nint(sqrt(real(nmodes,dp)))
    if(nmodes_dim*nmodes_dim/=nmodes)call errstop('SET_PHASES',&
     &'NMODES input parameter must have an integer square root in two&
     & dimensions.')
    allocate(theta_2d(nmodes_dim,nmodes_dim),stat=ialloc)
    if(ialloc/=0)call errstop('SET_PHASE','Allocation &
     &problem with 2D phase array.')

    select case (trim(adjustl(phase_format)))  

     case('input')
      if(esdf_block('phases',nlines))then
       if(nlines/=nmodes)call errstop('SET_PHASES','Should be &
        &NMODES lines in the PHASES block in input.')
       i=0
       do il=1,nmodes_dim
        do ik=1,nmodes_dim
         i=i+1
         read(block_data(i),*)theta_2d(ik,il)
        enddo
       enddo
      endif

     case('random','default','preset')

      do il=1,nmodes_dim
       do ik=1,nmodes_dim
        theta_2d(ik,il)=ranx_twopi()
       enddo
      enddo

    end select

   case(3)
    nmodes_dim=nint((real(nmodes,dp)**(1.d0/3.d0)))
    if(nmodes_dim*nmodes_dim*nmodes_dim/=nmodes)call &
     &errstop('SET_PHASES','NMODES input parameter must have &
     &an integer cube root in three dimensions.')
    allocate(theta_3d(nmodes_dim,nmodes_dim,nmodes_dim),stat=ialloc)
    if(ialloc/=0)call errstop('SET_PHASES','Allocation &
     &problem with 3D phase array.')

    select case (trim(adjustl(phase_format)))  

     case('input')
      if(esdf_block('phases',nlines))then
       if(nlines/=nmodes)call errstop('SET_PHASES','Should be &
        &NMODES lines in the PHASES block in input.')
       i=0
       do il=1,nmodes_dim
        do ik=1,nmodes_dim
         do ij=1,nmodes_dim
          i=i+1
          read(block_data(i),*)theta_3d(ij,ik,il)
         enddo
        enddo
       enddo
      endif

     case('random','default','preset')
      do il=1,nmodes_dim
       do ik=1,nmodes_dim
        do ij=1,nmodes_dim
         theta_3d(ij,ik,il)=ranx_twopi()
        enddo
       enddo
      enddo

    end select

  end select


  CONTAINS


   SUBROUTINE get_preset_randomseed
   IMPLICIT NONE
   INTEGER,PARAMETER :: ia=16807,im=2147483647,iq=127773,ir=2836,mask=123459786

   if(phase_preset==0)then
    randomseed=randomseed_in
    return
   endif
   randomseed=randomseed_in
   randomseed=ieor(randomseed,mask)
   do i=1,phase_preset
    k=randomseed/iq
    randomseed=ia*(randomseed-k*iq)-ir*k
    if(randomseed<0)randomseed=randomseed+im
   enddo
   randomseed=ieor(randomseed,mask)

   END SUBROUTINE get_preset_randomseed

 
  END SUBROUTINE set_phases


  SUBROUTINE set_weights
  USE random_numbers,ONLY : default_seed,randomseed,ranx_twopi,initialize_random
  IMPLICIT NONE
  INTEGER randomseed_in,idum1,idum2,ierr
  REAL(dp) sumsq_weights

  sumsq_weights=0
  ierr=0
  if(nmodes<=0)call errstop('SET_WEIGHTS','Number of wave function modes in & 
   &superposition must be set in input through keyword NMODES.') 
  if(weight_preset<0)call errstop('SET_WEIGHTS',&
   &'WEIGHT_PRESET must be a positive integer or zero.')

  if(trim(adjustl(phase_format))=='input')then
! Set random seed.
   select case (trim(adjustl(weight_format)))
    case('input','default')
     randomseed=default_seed
    case('random')
     if(am_master)call system_clock(randomseed,idum1,idum2)
     randomseed=default_seed
     call mpi_bcast(randomseed,1,mpi_integer,0,mpi_comm_world,ierr)
    case('preset')
     randomseed_in=default_seed
     call get_preset_randomseed
    case default
     call errstop('SET_WEIGHTS','Unknown value for WEIGHT_FORMAT in input.')
   end select

   if(trim(adjustl(weight_format))/='input')call initialize_random
  endif

! Define the weights.
  select case(dim)

   case(1)
    nmodes_dim=nmodes
    allocate(weights_1d(nmodes),stat=ialloc)
    if(ialloc/=0)call errstop('SET_WEIGHTS','Allocation &
     &problem with 1D weights array.')

    select case (trim(adjustl(weight_format)))  

     case('input')
      if(esdf_block('weights',nlines))then
       if(nlines/=nmodes)call errstop('SET_WEIGHTS','Should be &
        &nmodes lines in the WEIGHTS block in input.')
       do i=1,nmodes_dim
        read(block_data(i),*)weights_1d(i)
        sumsq_weights=sumsq_weights+weights_1d(i)**2.d0
       enddo
      endif

     case('random','default','preset')
 
      do i=1,nmodes_dim
       weights_1d(i)=ranx_twopi()
       sumsq_weights=sumsq_weights+weights_1d(i)**2.d0
      enddo

    end select

    weights_1d=weights_1d/sqrt(sumsq_weights*cell_x/2.d0)

   case(2)
    nmodes_dim=nint(sqrt(real(nmodes,dp)))
    if(nmodes_dim*nmodes_dim/=nmodes)call errstop('SET_WEIGHTS',&
     &'NMODES input parameter must have an integer square root in two&
     & dimensions.')
    allocate(weights_2d(nmodes_dim,nmodes_dim),stat=ialloc)
    if(ialloc/=0)call errstop('SET_WEIGHTS','Allocation &
     &problem with 2D weights array.')

    select case (trim(adjustl(weight_format)))  

     case('input')
      if(esdf_block('weights',nlines))then
       if(nlines/=nmodes)call errstop('SET_WEIGHTS','Should be &
        &NMODES lines in the WEIGHT block in input.')
       i=0
       do il=1,nmodes_dim
        do ik=1,nmodes_dim
         i=i+1
         read(block_data(i),*)weights_2d(ik,il)
         sumsq_weights=sumsq_weights+weights_2d(ik,il)**2.d0
        enddo
       enddo
      endif

     case('random','default','preset')

      do il=1,nmodes_dim
       do ik=1,nmodes_dim
        weights_2d(ik,il)=ranx_twopi()
        sumsq_weights=sumsq_weights+weights_2d(ik,il)**2.d0
       enddo
      enddo

    end select

    weights_2d=weights_2d/sqrt(sumsq_weights*cell_x*cell_y/4.d0)

   case(3)
    nmodes_dim=nint((real(nmodes,dp)**(1.d0/3.d0)))
    if(nmodes_dim*nmodes_dim*nmodes_dim/=nmodes)call &
     &errstop('SET_WEIGHTS','NMODES input parameter must have &
     &an integer cube root in three dimensions.')
    allocate(weights_3d(nmodes_dim,nmodes_dim,nmodes_dim),stat=ialloc)
    if(ialloc/=0)call errstop('SET_WEIGHTS','Allocation &
     &problem with 3D weights array.')

    select case (trim(adjustl(weight_format)))  

     case('input')
      if(esdf_block('weights',nlines))then
       if(nlines/=nmodes)call errstop('SET_WEIGHTS','Should be &
        &NMODES lines in the WEIGHTS block in input.')
       i=0
       do il=1,nmodes_dim
        do ik=1,nmodes_dim
         do ij=1,nmodes_dim
          i=i+1
          read(block_data(i),*)weights_3d(ij,ik,il)
          sumsq_weights=sumsq_weights+weights_3d(ij,ik,il)**2.d0
         enddo
        enddo
       enddo
      endif

     case('random','default','preset')
      do il=1,nmodes_dim
       do ik=1,nmodes_dim
        do ij=1,nmodes_dim
         weights_3d(ij,ik,il)=ranx_twopi()
         sumsq_weights=sumsq_weights+weights_3d(ij,ik,il)**2.d0
        enddo
       enddo
      enddo

    end select

    weights_3d=weights_3d/sqrt(sumsq_weights*cell_x*cell_y*cell_z/8.d0)

  end select


  CONTAINS


   SUBROUTINE get_preset_randomseed
   IMPLICIT NONE
   INTEGER,PARAMETER :: ia=16807,im=2147483647,iq=127773,ir=2836,mask=123459786

   if(phase_preset==0)then
    randomseed=randomseed_in
    return
   endif
   randomseed=randomseed_in
   randomseed=ieor(randomseed,mask)
   do i=1,phase_preset
    k=randomseed/iq
    randomseed=ia*(randomseed-k*iq)-ir*k
    if(randomseed<0)randomseed=randomseed+im
   enddo
   randomseed=ieor(randomseed,mask)

   END SUBROUTINE get_preset_randomseed

 
  END SUBROUTINE set_weights

 
  SUBROUTINE check_input_parameters
  IMPLICIT NONE

  if(trim(adjustl(calc_type))/='trajectory'.and.&
   &trim(adjustl(calc_type))/='density')call errstop&
   &('CHECK_INPUT_PARAMETERS','Input keyword CALC_TYPE invalid. Currently &
   &allowed values: "trajectory", "density".')
 
  if(dim/=1.and.dim/=2.and.dim/=3)call errstop&
   &('CHECK_INPUT_PARAMETERS','DIMENSIONALITY must be 1, 2 or 3')

  if(trim(adjustl(int_algorithm))/='runge-kutta'.and.&
   &trim(adjustl(int_algorithm))/='runge_kutta'.and.&
   &trim(adjustl(int_algorithm))/='Runge-Kutta'.and.&
   &trim(adjustl(int_algorithm))/='Runge_Kutta'.and.&
   &trim(adjustl(int_algorithm))/='bulirsch-stoer'.and.&
   &trim(adjustl(int_algorithm))/='bulirsch_stoer'.and.&
   &trim(adjustl(int_algorithm))/='Bulirsch-Stoer'.and.&
   &trim(adjustl(int_algorithm))/='Bulirsch_Stoer')call errstop&
   &('CHECK_INPUT_PARAMETERS','Keyword INT_ALGORITHM currently restricted to &
   &"runge-kutta" or "bulirsch-stoer".')

   if(trim(adjustl(int_algorithm))=='runge_kutta')int_algorithm='runge-kutta'
   if(trim(adjustl(int_algorithm))=='Runge-Kutta')int_algorithm='runge-kutta'
   if(trim(adjustl(int_algorithm))=='Runge_Kutta')int_algorithm='runge-kutta'
   if(trim(adjustl(int_algorithm))=='bulirsch_stoer')int_algorithm=&
    &'bulirsch-stoer'
   if(trim(adjustl(int_algorithm))=='Bulirsch-Stoer')int_algorithm=&
    &'bulirsch-stoer'
   if(trim(adjustl(int_algorithm))=='Bulirsch_Stoer')int_algorithm=&
    &'bulirsch-stoer'

  if(maxstp<1)call errstop('CHECK_INPUT_PARAMETERS',&
   &'Invalid value of MAXSTP.')

  if(init_eps<100.d0*mineps)call errstop('CHECK_INPUT_PARAMETERS',&
   &'INIT_EPS input parameter must be at least 100 times the internal &
   &MINEPS parameter (which should be 10^-12).')

  if(trim(adjustl(vel_type))/='deBB'.and.trim(adjustl(vel_type))/='curl1')&
   &call errstop('CHECK_INPUT_PARAMETERS','Unknown value for keyword &
   &VEL_TYPE. Currently allowed values: "deBB", "curl1".')

  if(dim==1.and.trim(adjustl(vel_type))=='curl1')call errstop&
   &('CHECK_INPUT_PARAMETERS','Velocity-type CURL1 not possible for 1D&
   & calculations.')

  if(trim(adjustl(vel_type))=='curl1'.and.curlweight==0.d0)vel_type='deBB'

  if(trim(adjustl(time_direction))/='forward'.and.trim(adjustl(time_direction))&
   &/='forwards'.and.trim(adjustl(time_direction))/='back'.and.&
   &trim(adjustl(time_direction))/='backward'.and.trim(adjustl(time_direction))&
   &/='backwards')call errstop('CHECK_INPUT_PARAMETERS',&
   &'Keyword TIME_DIRECTION currently restricted to "forward" or "back".')
  if(trim(adjustl(time_direction))=='forwards')time_direction='forward'
  if(trim(adjustl(time_direction))=='backward')time_direction='back'
  if(trim(adjustl(time_direction))=='backwards')time_direction='back'
 
  if(trim(adjustl(calc_type))=='trajectory')then

   if(trim(adjustl(time_direction))=='forward')then
    if(traj_time_end==-999.d0)then
     call errstop('CHECK_INPUT_PARAMETERS','For forwards trajectory &
      &calculations a value for the TRAJ_TIME_END keyword must be &
      &supplied in input, since no sensible default can be assumed.')
    endif
   else ! trim(adjustl(time_direction)=='back'
    if(traj_time_start==-999.d0)then
     call errstop('CHECK_INPUT_PARAMETERS','For backwards trajectory &
      &calculations a value for the TRAJ_TIME_START keyword must be &
      &supplied in input, since no sensible default can be assumed.')
    endif
   endif ! time_direction

   if(any(xstart(:)<=0.d0).or.any(xstart(:)>=cell_x))call errstop&
    &('CHECK_INPUT_PARAMETERS','All initial x components of position in the &
    &TRAJECTORY_START block must be greater than 0 and less than CELL_X.')
   if(dim>1.and.(any(ystart(:)<=0.d0).or.any(ystart(:)>=cell_y)))call &
    &errstop('CHECK_INPUT_PARAMETERS','All initial y components of &
    &position in the TRAJECTORY_START block must be greater than 0 and less &
    &than CELL_Y.')
   if(dim>2.and.(any(zstart(:)<=0.d0).or.any(zstart(:)>=cell_z)))call &
    &errstop('CHECK_INPUT_PARAMETERS','All initial z components of &
    &position in the TRAJECTORY_START block must be greater than 0 and less &
    &than CELL_Z.')

   if(traj_time_end==traj_time_start)call errstop&
    &('CHECK_INPUT_PARAMETERS','TRAJ_TIME_START equals TRAJ_TIME_END. No &
    &propagation required,')

   if(trim(adjustl(time_direction))=='forward')then
    if(traj_time_end<traj_time_start)then
     call errstop('CHECK_INPUT_PARAMETERS','Forward-in-time trajectory &
      &calculation but TRAJ_TIME_END is earlier than TRAJ_TIME_START.')
    endif
   else ! trim(adjustl(time_direction)=='back'
    if(traj_time_start<traj_time_end)then
     call errstop('CHECK_INPUT_PARAMETERS','Backward-in-time trajectory&
      & calculation but TRAJ_TIME_END is earlier than TRAJ_TIME_START.')
    endif
   endif

   if(plot_all_traj.and.ntrajectories>1)call errstop&
    &('CHECK_INPUT_PARAMETERS','The PLOT_ALL_TRAJ option must be set to false&
    & if there is more than one trajectory.')

   if(dtsave<0.d0)call errstop('CHECK_INPUT_PARAMETERS','Invalid&
    & value of DTSAVE keyword.')

   if(phase_noise<0.or.phase_noise>14)call errstop('CHECK_INPUT_PARAMETERS',&
    &'PHASE_NOISE keyword out of range (0:14).')

  else ! calc_type==density 

   if(den_ntimes<1)call errstop('CHECK_INPUT_PARAMETERS','Invalid value &
    &for DEN_NTIMES in input.')

   if(den_ntimes==1.and.(den_time_start/=den_time_end))call errstop&
    &('CHECK_INPUT_PARAMETERS','If number of density plots DEN_NTIMES equals 1,&
    & then DEN_TIME_START must equal DEN_TIME_END in input.')

   if(den_ntimes>1.and.(den_time_start==den_time_end))call errstop&
    &('CHECK_INPUT_PARAMETERS','If number of density plots DEN_NTIMES > 1,&
    & then DEN_TIME_START must not equal DEN_TIME_END in input.')

   select case(dim)
    case(1)
     if(cell_x<=0.d0)call errstop('CHECK_INPUT_PARAMETERS',&
      &'Invalid value for keyword CELL_X in input.')
     if(cell_y/=0.d0.or.cell_z/=0.d0)call errstop&
      &('CHECK_INPUT_PARAMETERS','Input parameters CELL_Y and &
      &CELL_Z must be zero for a 1D calculation.')
    case(2) 
     if(cell_x<=0.d0.or.cell_y<=0.d0)call errstop&
      &('CHECK_INPUT_PARAMETERS','Invalid value for keyword CELL_X and/or&
      & CELL_Y in input.')
     if(cell_z/=0.d0)call errstop('CHECK_INPUT_PARAMETERS',&
     &'Input parameter CELL_Z must be zero for a 2D calculation.')
     if(cell_x/=cell_y)call errstop('CHECK_INPUT_PARAMETERS',&
      &'Input parameters CELL_X and CELL_Y must be equal (cell for the moment&
      & assumed to be square in 2D).')
    case(3)
     if(cell_x<=0.d0.or.cell_y<=0.d0.or.cell_z<=0.d0)call &
      &errstop('CHECK_INPUT_PARAMETERS','Invalid value for keywords&
      & CELL_X, CELL_Y, CELL_Z in input.')
     if(cell_x/=cell_y.or.cell_z/=cell_x)call errstop&
      &('CHECK_INPUT_PARAMETERS','Input parameters CELL_X and CELL_Y must be &
      &equal (cell for the moment assumed to be a cube in 3D).')
   end select

   if(cell_x/=pi)then
    write(o,*)'Definition of pi: ',pi
    call errstop('CHECK_INPUT_PARAMETERS','For the moment the &
     &CELL_X input parameter must equal pi. If you would like this changed &
     &please ask Mike.')
   endif
 
   if(nlattice<1)call errstop('CHECK_INPUT_PARAMETERS','Invalid value&
    & for keyword NLATTICE in input.')

   if(any(ncgrain(:)<1))call errstop('CHECK_INPUT_PARAMETERS','Invalid value&
    & for one or more COARSE_GRAINING_LENGTHS in input.')

   if(any(ncgrain(:)==1))call errstop('CHECK_INPUT_PARAMETERS','Invalid value&
    & for one or more COARSE_GRAINING_LENGTHS in input. Minimum value is 2.')

   if(any(ncgrain(:)>nlattice))call errstop('CHECK_INPUT_PARAMETERS','&
    &COARSE_GRAINING_LENGTHS may not be larger than keyword NLATTICE in input.')
 
   if(any(mod(nlattice,ncgrain(:))/=0))call errstop('CHECK_INPUT_PARAMETERS',&
    &'NLATTICE keyword must be an integer multiple of each of the&
    & COARSE_GRAINING LENGTHS (e.g. 1024,32).')

   if(nscgrain<1)call errstop('CHECK_INPUT_PARAMETERS','Invalid value&
    & for keyword NSCGRAIN in input.')

   if(nsmoothstep<1)call errstop('CHECK_INPUT_PARAMETERS','Invalid value&
    & for keyword NSMOOTHSTEP in input.')

   if(nscgrain>nlattice)call errstop('CHECK_INPUT_PARAMETERS','Keyword&
    & NSCGRAIN may not be larger than keyword NLATTICE in input.')
 
   if(mod(nlattice,nscgrain)/=0)call errstop('CHECK_INPUT_PARAMETERS',&
    &'NLATTICE keyword must be an integer multiple of NSCGRAIN (e.g. 1024,64).')

   t1=real((nlattice-nscgrain)/nsmoothstep+1,dp)
   if(anint(t1)/=t1)call errstop('CHECK_INPUT_PARAMETERS',&
    &'Invalid value for NSMOOTHSTEP keyword. You must therefore ensure that&
    & (NLATTICE-NSCGRAIN)/NSMOOTHSTEP + 1 is a whole number to make sure the&
    & translated cells fit in the box.')

   if(den_time_start>den_time_end)call errstop('CHECK_INPUT_PARAMETERS',&
    &'Keyword DEN_TIME_START is later than DEN_TIME_END in input.')

   if(read_backtracked.and.save_backtracked)call errstop&
    &('CHECK_INPUT_PARAMETERS','READ_BACKTRACKED and SAVE_BACKTRACKED keywords&
    & may not both be T.')

   if(dim==3.and.read_backtracked)call errstop('CHECK_INPUT_PARAMETERS',&
    &'The READ_BACKTRACKED facility is not available in 3 dimensions.')

   if(dim==3.and.save_backtracked)call errstop('CHECK_INPUT_PARAMETERS',&
    &'The SAVE_BACKTRACKED facility is not available in 3 dimensions.')
 
   if(dim==3.and.plot_smooth)call errstop('CHECK_INPUT_PARAMAETERS',&
    &'The PLOT_SMOOTH facility is not available in 3 dimensions, &
    &though it could easily be coded up if anyone is interested.')

   if(dim==3.and.nnodes<2)call errstop('CHECK_INPUT_PARAMETERS','3D density &
    &computation must be run on more than 1 node, since only slave nodes &
    &do any real work with the chosen algorithm. With realistic &
    &system sizes, you can make that "considerably more than 1 node" if &
    &you want to avoid blowing the memory or dying of boredom.')

   if(.not.plot_raw.and..not.plot_cg.and..not.plot_smooth.and..not.hfunction)&
    &call errstop('CHECK_INPUT_PARAMETERS','Density mode flagged but &
    &PLOT_RAW, PLOT_CG, PLOT_SMOOTH, and HFUNCTION are all turned off, which &
    &means that the code will chug away for ages and then not produce &
    &anything. This has the same net result as not doing the calculation &
    &at all so I''ll stop now in order to save electricity.')

  endif ! trajectory or density

  if(trim(adjustl(wfn_type))/='sine_wave'.and.trim(adjustl(wfn_type))/=&
   &'sin_wave'.and.trim(adjustl(wfn_type))/='scaled_sin'.and.&
   &trim(adjustl(wfn_type))/='scaled_sine'&
   &)call errstop('CHECK_INPUT_PARAMETERS',&
   &'Currently only "sine_wave" and "scaled_sine" are allowed for keyword &
   &WFN_TYPE in input.')

  if(trim(adjustl(calc_type))=='trajectory'.and..not.&
   &esdf_block('trajectory_start',nlines))call errstop&
   &('CHECK_INPUT_PARAMETERS','TRAJECTORY_START block must exist in input for&
   & CALC_TYPE="trajectory".')

  if(trim(adjustl(phase_format))=='input'.and.&
   &(trim(adjustl(wfn_type))=='sine_wave'.or.trim(adjustl(wfn_type))==&
   &'sin_wave'.or.trim(adjustl(wfn_type))=='scaled_sin'.or.&
   &trim(adjustl(wfn_type))=='scaled_sine').and..not.&
   &esdf_block('phases',nlines))call errstop&
   &('CHECK_INPUT_PARAMETERS','Phases block must exist in input for&
   & WFN_TYPE="sine_wave" or WFN_TYPE="scaled_sine".')

  if(trim(adjustl(weight_format))=='input'.and.&
   &(trim(adjustl(wfn_type))=='scaled_sine'.or.trim(adjustl(wfn_type))==&
   &'scaled_sin').and..not.esdf_block('weights',nlines))call errstop&
   &('CHECK_INPUT_PARAMETERS','Weights block must exist in input for&
   & WFN_TYPE="scaled_sine".')

  END SUBROUTINE check_input_parameters

 END MODULE input

 PROGRAM louis
!-----------------------------------------------------------------------------!
! LOUIS                                                                       !
! =====                                                                       !
!                                                                             !
! Version : 1.1                                                               !
!                                                                             !
! Program to propagate quantum particles moving according to de Broglie-Bohm  !
! dynamics (and appropriate generalizations).                                 !
!                                                                             !
! Program will plot trajectories, calculate time-dependent densities,         !
! psi squared, and H functions.                                               !
!                                                                             !
! Mike Towler (Feb 2010)                                                      !
!                                                                             !
! CHANGES                                                                     !
! -------                                                                     !
! None yet.                                                                   !
!-----------------------------------------------------------------------------!
 USE bulirsch_stoer
 USE dsp
 USE store
 USE esdf
 USE numerical
 USE parallel
 USE runge_kutta
 USE choose_wfn
 USE eval_density
 USE sin_wfn
 USE scaled_sin_wfn
 USE format_utils,ONLY : i2s,i2s64,r2s,d2s0,global_time_heading
 USE input,ONLY : set_input_parameters, check_input_parameters
 USE random_numbers,ONLY : randomseed,ranx_pm,initialize_random
 USE run_control,ONLY: errstop,errstop_master,errstop2,timer,timer_start,&
  &timer_end,reallocate
 IMPLICIT NONE

! Local variables
 INTEGER i,j,k,ii,jj,kk,ij,ik,n,nn,nnn,ncg,nc,ng,nplot,nfail,nstep,noscgcell,&
  &noscgcellsq,ngoodcg,ncgrainsq,nscgrainsq,ncgrain3,nscgrain3,ncgrow,nscgrow,&
  &nbase,is,ierr,nzero,nzero_total,nextra,mdttype1,mdttype2,nbins,maxstp_in,&
  &epsit,ncgcellsq,ncgcell3,nminh,nminh_total,kbase
 INTEGER(i64) noktot,nbadtot,noktot_total,nbadtot_total
 INTEGER(i64),ALLOCATABLE :: ntemp(:)
 INTEGER,ALLOCATABLE :: is_per_node(:),nfail_total(:),nhist(:),nhist_total(:),&
  &ncgcell(:)
 REAL(dp) t1,t2,t3,tt,den_time,den_timestep,cell_over_nlat,xdiff,rho,rho_0,&
  &psi2_0,psi2_t,one_over_ncgrain,one_over_nscgrain,one_over_ncgrainsq,&
  &one_over_nscgrainsq,one_over_ncgrain3,one_over_nscgrain3,nhist_scale,&
  &fita,fitb,ua,ub,chisq,emergency_maxdiff,val_interp,err_interp,points(10),&
  &repair(10)
 REAL(dp),ALLOCATABLE :: x0(:),x1(:),xlast(:),density_raw_n(:),density_raw(:),&
  &density_cg(:),density_smooth(:),psi_squared_raw_n(:),&
  &psi_squared_raw(:),psi_squared_cg(:),psi_squared_smooth(:),H(:,:),lnH(:),&
  &ttt(:),htemp(:),htemp2(:),xback1_n(:),xback1(:),xback2_n(:),xback2(:),&
  &xsave(:,:)
 LOGICAL failed,file_present,interpolation_required,nminhfail
 CHARACTER(24) tmpr,tmpr2,tmpr3
 CHARACTER(40) filename
! 3D comms groups
 INTEGER nredistgrps,rgmymaster,rggroup,my_rgnode,my_rgnode2,rgnnodes
 INTEGER,ALLOCATABLE :: rg_nnodes(:),rgstart(:),rgstop(:),rg_comm(:),&
  &rgrequest(:),rgrequest2(:),rgrequest3(:),which_rgnodes(:)
 LOGICAL rgmaster
 LOGICAL,ALLOCATABLE :: good_traj_save(:),stopped(:)

! Initialize MPI.
 ierr=0
 call init_parallel

! Open output file(s).
 o=7
 if(am_master)then
  output_file='out'
 else
  output_file='.out_node'//trim(i2s(my_node))
 endif
 open(unit=o,file=trim(output_file),status='unknown',position='append')

! Initialize timing routines.
 call timer_start ! proper timer
 
 call timer('SETUP',.true.)

 if(am_master)then
  call global_time_heading(.true.)
  write(o,'(1x,72(''=''))')
  write(o,*)' _                _'
  write(o,*)'| |    ___  _   _(_)___'
  write(o,*)'| |   / _ *| | | | / __|'
  write(o,*)'| |__| (_) | |_| | *__ *'
  write(o,*)'|_____*___/ *__,_|_|___*  1.2'
  write(o,*)
  write(o,*)'M.D. Towler, University of Cambridge (2010)'
  write(o,*)
 endif

! Get input parameters, check and do relevant setup.

 open_unit=.false.
 call esdf_init('input') ! read the input file (on each processor)
 if(am_master)call esdf_warnout
 call set_input_parameters
 if(am_master)call check_input_parameters ! for detectable errors
 call initial_setup
 call esdf_close
 if(trim(adjustl(phase_format))=='input'.and.phase_noise>0)&
  &call initialize_random

 if(am_master)then
  if(nnodes==1)then
   write(o,*)'Running on single processor.'
  else
   if(trim(adjustl(calc_type))=='density')then
    write(o,*)'Running in parallel on ',trim(i2s(nnodes)),' processors.'
   else
    write(o,*)'Parallel calculation requested - not necessary in &
     &trajectory mode.'
    write(o,*)'Only the master node will do any work.'
   endif
  endif
  write(o,*)
 endif

! Specify defaults for integration parameters.

 eps=init_eps             ! Initial value of overall tolerance level   
 h1=init_h                ! Guessed first stepsize
 maxdiff=converge_maxdiff ! Distance tolerance for convergence
 hmin=1.d-14              ! Min allowed stepsize
 maxstp_in=maxstp 
 emergency_maxdiff=10*maxdiff

 call timer('SETUP',.false.)

! Call driver routine with appropriate velocity, integrator and dimensionality.

 select case(trim(adjustl(int_algorithm)))

  case ('runge-kutta')

   select case(trim(adjustl(vel_type)))

    case('deBB')
     select case(dim)
      case(1)
       select case(trim(adjustl(wfn_type)))
        case('sine_wave')
         call driver(sw_velocity_deBB_1d,stepper_rk)
        case('scaled_sine')
         call driver(ss_velocity_deBB_1d,stepper_rk)
       end select
      case(2)
       select case(trim(adjustl(wfn_type)))
        case('sine_wave')
         call driver(sw_velocity_deBB_2d,stepper_rk)
        case('scaled_sine')
         call driver(ss_velocity_deBB_2d,stepper_rk)
       end select
      case(3)
       select case(trim(adjustl(wfn_type)))
        case('sine_wave')
         call driver(sw_velocity_deBB_3d,stepper_rk)
        case('scaled_sine')
         call driver(ss_velocity_deBB_3d,stepper_rk)
       end select
     end select

    case('curl1')
     select case(dim)
      case(1)
       call errstop('LOUIS','Velocity-type CURL1 not possible in 1D.')
      case(2)
       select case(trim(adjustl(wfn_type)))
        case('sine_wave')
         call driver(sw_velocity_curl1_2d,stepper_rk)
        case('scaled_sine')
         call driver(ss_velocity_curl1_2d,stepper_rk)
       end select
      case(3)
       call errstop('LOUIS','Curl1 velocity routine not yet coded for 3D.')
!      call driver(velocity_curl1_3d,stepper_rk)
     end select

    case default
     call errstop_master('LOUIS','Unknown velocity type.')

   end select

  case ('bulirsch-stoer')

   select case(trim(adjustl(vel_type)))

    case('deBB')
     select case(dim)
      case(1)
       select case(trim(adjustl(wfn_type)))
        case('sine_wave')
         call driver(sw_velocity_deBB_1d,stepper_bs)
        case('scaled_sine')
         call driver(ss_velocity_deBB_1d,stepper_bs)
       end select
      case(2)
       select case(trim(adjustl(wfn_type)))
        case('sine_wave')
         call driver(sw_velocity_deBB_2d,stepper_bs)
        case('scaled_sine')
         call driver(ss_velocity_deBB_2d,stepper_bs)
       end select
      case(3)
       select case(trim(adjustl(wfn_type)))
        case('sine_wave')
         call driver(sw_velocity_deBB_3d,stepper_bs)
        case('scaled_sine')
         call driver(ss_velocity_deBB_3d,stepper_bs)
       end select
     end select

    case('curl1')
     select case(dim)
      case(1)
       call errstop('LOUIS','Velocity-type CURL1 not possible in 1D.')
      case(2)
       select case(trim(adjustl(wfn_type)))
        case('sine_wave')
         call driver(sw_velocity_curl1_2d,stepper_bs)
        case('scaled_sine')
         call driver(ss_velocity_curl1_2d,stepper_bs)
       end select
      case(3)
       call errstop('LOUIS','Curl1 velocity routine not yet coded for 3D.')
!      call driver(velocity_curl1_3d,stepper_bs)
     end select

    case default
     call errstop_master('LOUIS','Unknown velocity type.')

   end select

  case default
   call errstop('LOUIS','Unknown integration algorithm.')

 end select

 if(am_master)then
  call timer_end
  call global_time_heading(.false.)
 endif

 call end_parallel
 
 stop


 CONTAINS


  SUBROUTINE initial_setup
!------------------------------------------------------------!
! Basic setup (allocation, parallel stuff, basic constants). !
!------------------------------------------------------------!
  IMPLICIT NONE

  if(negphase)then
    select case(dim)
     case(1)
      theta_1d=-theta_1d
     case(2)
      theta_2d=-theta_2d
     case(3)
      theta_3d=-theta_3d
    end select
  endif

  if(transposephase)then
   if(dim==2)then
    theta_2d=transpose(theta_2d)
   else
    call errstop_master('INITIAL_SETUP','The TRANSPOSEPHASE keyword only&
     & works in two dimensions.') 
   endif
  endif

  two_over_pi_n=two_over_pi/(real(nmodes_dim,dp))
  four_over_pi2_n2=four_over_pi_squared/(real(nmodes_dim**2,dp))
  eight_over_pi3_n3=eight_over_pi_cubed/(real(nmodes_dim**3,dp))
  half_cell_x=0.5d0*cell_x
  half_cell_y=0.5d0*cell_y
  half_cell_z=0.5d0*cell_z

! Allocate stuff.
  allocate(x0(dim),x1(dim),xlast(dim),stat=ialloc)
  if(ialloc/=0)call errstop_master('INITIAL_SETUP','Allocation error <1>.')

  if(trim(adjustl(calc_type))=='density')then

   time_direction='backward'
   nlatticesq=nlattice*nlattice
   allocate(ncgcell(num_ncgrain),stat=ialloc)
   if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to&
    & allocate NCGCELL vector.')
   do i=1,num_ncgrain
    ncgcell(i)=nlattice/ncgrain(i)  ! no of coarse-graining cells in the box
   enddo
   noscgcell=(nlattice-nscgrain)/nsmoothstep+1 ! no of overlapping smoothed 
                                         ! coarse-graining cells in the box
   noscgcellsq=noscgcell*noscgcell

   select case(dim)
    case(1)
     allocate(latpos1(nlattice),stat=ialloc)
     if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
      &allocate LATPOS1 vector. ')
     if(plot_cg.or.hfunction)then
      allocate(latpos1_cg(maxval(ncgcell(:))),stat=ialloc)
      if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
       &allocate LATPOS1_CG vector.')
     endif
     if(plot_smooth)then
      allocate(latpos1_smooth(noscgcell),stat=ialloc)
      if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
       &allocate LATPOS1_SMOOTH vector.')
     endif
    case(2)
     allocate(latpos2(2,nlatticesq),stat=ialloc)
     if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
      &allocate LATPOS2 vector.')
     if(plot_cg.or.hfunction)then
      allocate(latpos2_cg(2,maxval(ncgcell(:))**2),stat=ialloc)
      if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
       &allocate LATPOS2_CG vector.')
     endif
     if(plot_smooth)then
      allocate(latpos2_smooth(2,noscgcellsq),stat=ialloc)
      if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
       &allocate LATPOS2_SMOOTH vector.')
     endif
    case(3)
     allocate(latpos2(3,nlatticesq*ncgrain(1)),stat=ialloc)
     if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
      &allocate LATPOS2 vector.')
     if(plot_cg.or.hfunction)then
      allocate(latpos2_cg(3,ncgcell(1)**2),stat=ialloc)
      if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
       &allocate LATPOS2_CG vector.')
     endif
     if(plot_smooth)then
      allocate(latpos2_smooth(3,noscgcellsq),stat=ialloc)
      if(ialloc/=0)call errstop_master('INITIAL_SETUP','Insufficient memory to &
       &allocate LATPOS2_SMOOTH vector.')
     endif
    case default
     call errstop_master('INITIAL_SETUP','Unknown dimensionality.')
   end select

   if(nnodes>1)then

    if(dim==3)then
     if(am_master.and.nnodes<ncgcell(1)*2)then
      call errstop('INITIAL_SETUP','Number of processors insufficient for 3d&
       & calculation. Current algorithm requires a number of processors at &
       & least equal to twice the number of coarse-graining cells per dimension&
       & (and preferably considerably more).')
     endif
     call setup_comm_groups! Setup communicator groups for parallel 3D algorithm
     if(rgnnodes-1>nlattice*ncgrain(1).and.am_master)&
      &call errstop2('INITIAL_SETUP','For this system size LOUIS cannot &
      &exploit more than the following number of processors: ',&
      &nlatticesq+ncgcell(1))
    else
     allocate(is_per_node(0:nnodes-1),stat=ialloc)
     if(ialloc/=0)call errstop_master('INITIAL_SETUP','Cannot allocate &
      &IS_PER_NODE vector.')
     is=nlattice**dim
     is_per_node(0:nnodes-1)=is/nnodes
     nextra=mod(is,nnodes)
     do i=0,nextra-1
      is_per_node(i)=is_per_node(i)+1
     enddo
     n=is_per_node(0) ! equivalent to n=maxval(is_per_node(0:nnodes-1))
     allocate(density_raw_n(n),psi_squared_raw_n(n),stat=ialloc)
     if(ialloc/=0)call errstop('INITIAL_SETUP','Unable to allocate &
      &DENSITY_RAW_N and PSI_SQUARED_RAW_N vectors.')
     if(save_backtracked)then
      allocate(xback1_n(n),stat=ialloc)
      if(ialloc/=0)call errstop('INITIAL_SETUP','Unable to allocate XBACK1_N &
       &vector.')
      if(dim==2)allocate(xback2_n(n),stat=ialloc)
       if(ialloc/=0)call errstop('INITIAL_SETUP','Unable to allocate XBACK2_N &
      &vector.')
     endif
    endif

   else ! 1 node - not allowed for 3D case

    is_per_node(0)=nlattice**dim
    allocate(density_raw_n(is_per_node(0)),psi_squared_raw_n(is_per_node(0)),&
     &stat=ialloc)
    if(ialloc/=0)call errstop('INITIAL_SETUP','Unable to allocate &
     &DENSITY_RAW_N and PSI_SQUARED_RAW_N vectors.')
    if(save_backtracked)then
     allocate(xback1_n(is_per_node(0)),stat=ialloc)
     if(ialloc/=0)call errstop('INITIAL_SETUP','Unable to allocate XBACK1_N &
      &vector.')
     if(dim==2)allocate(xback2_n(is_per_node(0)),stat=ialloc)
     if(ialloc/=0)call errstop('INITIAL_SETUP','Unable to allocate XBACK2_N &
      &vector.')
    endif

   endif ! parallel or not

   allocate(nfail_total(den_ntimes),stat=ialloc)
   if(ialloc/=0)call errstop('INITIAL_SETUP','Unable to allocate NFAIL_TOTAL&
    & vector.')

   if(hfunction)then
    allocate(H(num_ncgrain,den_ntimes),lnH(den_ntimes),&
     &ttt(den_ntimes),stat=ialloc)
    if(ialloc/=0)call errstop('INITIAL_SETUP','Unable to allocate H vector.')
   endif

   if(nstep_histogram)then
    nbins=1000
    nhist_scale=real(nbins,dp)/real(maxstp,dp)
    allocate(nhist(nbins),stat=ialloc)
    if(ialloc/=0)call errstop_master('INITIAL_SETUP','Unable to allocate NHIST&
     & vector.')
    nhist(:)=0
   endif

  endif ! density

  if(trim(adjustl(calc_type))=='trajectory')then
   if(trim(adjustl(time_direction))=='forward')then
    n=int((traj_time_end-traj_time_start)/(fourpi+1.d-6))+1
   else
    n=int((traj_time_start-traj_time_end)/(fourpi+1.d-6))+1
   endif
   maxstp=n*maxstp
  endif

  END SUBROUTINE initial_setup


  SUBROUTINE driver(velocity,stepper)
!-------------------------------------------------------------------------!
! Print information and select calculation type.                          !
!-------------------------------------------------------------------------!
  IMPLICIT NONE

  INTERFACE
   SUBROUTINE velocity(t,x,vel)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
   SUBROUTINE stepper(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
    REAL(dp),INTENT(inout) :: x(dim)
    REAL(dp),INTENT(in) :: htry,eps
    REAL(dp),INTENT(out) :: hdid,hnext
    REAL(dp),INTENT(inout) :: t
    INTERFACE
     SUBROUTINE velocity(t,x,vel)
      USE dsp
      USE store, ONLY : dim
      IMPLICIT NONE
      REAL(dp),INTENT(in) :: t
      REAL(dp),INTENT(in) :: x(dim)
      REAL(dp),INTENT(out) :: vel(dim)
     END SUBROUTINE velocity
    END INTERFACE
   END SUBROUTINE stepper
  END INTERFACE
  if(am_master)then
   select case(dim)
    case(1)
     write(o,*)'DIMENSIONALITY         : 1d'
    case(2)
     write(o,*)'DIMENSIONALITY         : 2d'
    case(3)
     write(o,*)'DIMENSIONALITY         : 3d'
    case default
     call errstop('DRIVER','Unknown dimensionality.')
   end select
   select case(trim(adjustl(int_algorithm)))
    case ('runge-kutta')
     write(o,*)'INTEGRATION ALGORITHM  : Runge-Kutta'
    case ('bulirsch-stoer')
     write(o,*)'INTEGRATION ALGORITHM  : Bulirsch-Stoer'
     write(o,*)' WARNING: large stepsizes - plots of trajectories&
      & will not look smooth.'
    case default
     call errstop('DRIVER','Unknown integration algorithm.')
   end select
   select case(trim(adjustl(vel_type)))
    case('deBB')
     write(o,*)'VELOCITY TYPE          : Standard de Broglie-Bohm'
    case('curl1')
     write(o,*)'VELOCITY TYPE          : De Broglie-Bohm plus curl |Psi^2|'
    case default
     call errstop('DRIVER','Unknown velocity type.')
   end select
   select case(trim(adjustl(wfn_type)))
    case('sine_wave')
     write(o,*)'WAVEFUNCTION TYPE      : sine wave'
    case('scaled_sine')
     write(o,*)'WAVEFUNCTION TYPE      : scaled sine wave'
    case default
     call errstop('DRIVER','Unknown wavefunction type,')
   end select
 
   if(trim(adjustl(calc_type))=='density')then
    if(num_ncgrain==1)then
     t1=(real(ncgrain(1),dp)*cell_x)/real(nlattice,dp)
     tmpr=r2s(t1,'(f12.6)')
     write(o,*)'COARSE-GRAINING LENGTH : ',trim(tmpr)
    endif
   endif

   call energy_and_variance
   tmpr=r2s(mean_energy,'(f12.6)')
   write(o,*)'MEAN ENERGY            : ',trim(tmpr)
   tmpr=r2s(variance,'(f12.6)')
   write(o,*)'VARIANCE               : ',trim(tmpr)
   tmpr2=r2s(cvariance,'(f12.6)')
   if(trim(tmpr)/=trim(tmpr2))then
    write(o,*)'VARIANCE (METHOD 2)    : ',trim(tmpr)
   endif

   if(fastmode)then
    write(o,*)'FASTMODE               : T'
   else 
    write(o,*)'FASTMODE               : F'
   endif

   write(o,*)'PHASE_NOISE            : ',trim(i2s(phase_noise))

   select case(trim(adjustl(phase_format)))
    case('input')
     write(o,*)'PHASE_FORMAT           : input'
    case('random')
     write(o,*)'PHASE_FORMAT           : random'
    case('default')
     write(o,*)'PHASE_FORMAT           : default'
    case('preset')
     write(o,*)'PHASE_FORMAT           : preset'
     write(o,*)'PHASE_PRESET           : ',trim(i2s(phase_preset))
   end select
   write(o,*)'NUMBER OF MODES        : ',trim(i2s(nmodes))
   if(trim(adjustl(phase_format))/='input')then
    write(o,*)'RANDOM SEED            : ',trim(i2s(randomseed))
   endif
   write(o,*)'INITIAL PHASES'
   select case(dim)
    case(1)
     do i=1,nmodes_dim
      write(o,'(2x,a)')d2s0(theta_1d(i))
     enddo
    case(2)
     do j=1,nmodes_dim
      do i=1,nmodes_dim
       write(o,'(2x,a)')d2s0(theta_2d(i,j))
      enddo
     enddo
    case(3)
     do k=1,nmodes_dim
      do j=1,nmodes_dim
       do i=1,nmodes_dim
        write(o,'(2x,a)')d2s0(theta_3d(i,j,k))
       enddo
      enddo
     enddo
   end select
   if(trim(adjustl(wfn_type))=='scaled_sine')then
    select case(trim(adjustl(weight_format)))
     case('input')
      write(o,*)'WEIGHTS_FORMAT          : input'
     case('random')
      write(o,*)'WEIGHTS_FORMAT          : random'
     case('default')
      write(o,*)'WEIGHTS_FORMAT          : default'
     case('preset')
      write(o,*)'WEIGHTS_FORMAT          : preset'
      write(o,*)'WEIGHTS_PRESET          : ',trim(i2s(weight_preset))
    end select
    if(trim(adjustl(phase_format))=='input'.and.&
     &trim(adjustl(weight_format))/='input')then
     write(o,*)'RANDOM SEED            : ',trim(i2s(randomseed))
    endif
    write(o,*)'WEIGHTS'  
    select case(dim)
     case(1)
      do i=1,nmodes_dim
       write(o,'(2x,a)')d2s0(weights_1d(i))
      enddo
     case(2)
      do j=1,nmodes_dim
       do i=1,nmodes_dim
        write(o,'(2x,a)')d2s0(weights_2d(i,j))
       enddo
      enddo
     case(3)
      do k=1,nmodes_dim
       do j=1,nmodes_dim
        do i=1,nmodes_dim
         write(o,'(2x,a)')d2s0(weights_3d(i,j,k))
        enddo
       enddo
      enddo
    end select
   endif
  endif

  select case(trim(adjustl(calc_type)))
   case('trajectory')
    if(am_master)call trajectory_driver(velocity,stepper)
   case('density')
    call density_driver(velocity,stepper)
   case default
    call errstop_master('DRIVER', 'Type of calculation not recognized. &
     &Check CALC_TYPE value in the input file.')
  end select
  
  END SUBROUTINE driver


  SUBROUTINE trajectory_driver(velocity,stepper)
!--------------------------------!
! Drive trajectory calculations. !
!--------------------------------!
  IMPLICIT NONE
  INTEGER ntraj
  LOGICAL op,warned_fastmode

  INTERFACE
   SUBROUTINE velocity(t,x,vel)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
   SUBROUTINE stepper(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
    REAL(dp),INTENT(inout) :: x(dim)
    REAL(dp),INTENT(in) :: htry,eps
    REAL(dp),INTENT(out) :: hdid,hnext
    REAL(dp),INTENT(inout) :: t
    INTERFACE
     SUBROUTINE velocity(t,x,vel)
      USE dsp
      USE store, ONLY : dim
      IMPLICIT NONE
      REAL(dp),INTENT(in) :: t
      REAL(dp),INTENT(in) :: x(dim)
      REAL(dp),INTENT(out) :: vel(dim)
     END SUBROUTINE velocity
    END INTERFACE
   END SUBROUTINE stepper
  END INTERFACE

  call timer('TRAJECTORY PROPAGATION',.true.)

  write(o,*)
  write(o,*)'Trajectory mode'
  write(o,*)'==============='
  write(o,*)
  write(o,*)'Number of trajectories : ',trim(i2s(ntrajectories))
  if(traj_time_start==0.d0)then
   write(o,*)'TIME START             : 0.0'
  else
   tmpr=r2s(traj_time_start,'(f25.16)')
   write(o,*)'TIME START             : ',trim(tmpr)
  endif
  tmpr=r2s(traj_time_end,'(f25.16)')
  write(o,*)'TIME END               : ',trim(tmpr)

  do ntraj=1,ntrajectories

   write(o,*)
   if(ntrajectories>1)then
    write(o,*)'TRAJECTORY ',trim(i2s(ntraj))
    if(ntraj<10)then
     write(o,*)'------------'
    elseif(ntraj<100)then
     write(o,*)'-------------'
    else
     write(o,*)'--------------'
    endif
   endif

   if(ntraj>1.and.phase_noise>0)call add_noise_to_phases

   if(plot_all_traj)then
    save_steps=.true.
   else
    save_steps=.false. ! don't save 1st attempt unless plot of all requested
   endif

   x0(1)=xstart(ntraj) 
   if(dim>1)x0(2)=ystart(ntraj) 
   if(dim==3)x0(3)=zstart(ntraj)

   write(o,*)'INITIAL POSITION'
   tmpr=r2s(x0(1),'(f25.16)')
   select case(dim)
    case(1)
     write(o,"(2x,a)")trim(tmpr)
     call eval_psi_squared_1d(traj_time_start,x0,rho)
    case(2)
     tmpr2=r2s(x0(2),'(f25.16)')
     write(o,"(1x,2(1x,a))")trim(tmpr),trim(tmpr2)
     call eval_psi_squared_2d(traj_time_start,x0,rho)
    case(3)
     tmpr2=r2s(x0(2),'(f25.16)')
     tmpr3=r2s(x0(3),'(f25.16)')
     write(o,"(1x,3(1x,a))")trim(tmpr),trim(tmpr2),trim(tmpr3)
     call eval_psi_squared_3d(traj_time_start,x0,rho)
   end select

   write(o,*)
   write(o,*)'SQUARE OF WAVE FUNCTION AT THIS POINT'
   tmpr=r2s(rho,'(f25.16)')
   write(o,'(2x,a)')trim(tmpr)
   if(rho==0.d0)then
    if(ntrajectories>1)then
     write(o,*)
     write(o,*)'Selected starting point lies on a wave function node.'
     write(o,*)'Cannot start here so will skip current trajectory.'
     write(o,*)
     cycle
    else
     call errstop('TRAJECTORY_DRIVER','Selected starting point lies on a wave &
      &function node. Cannot start trajectory here.')
    endif
   endif

   xlast=1.d10
   xdiff=1.d10
   warned_fastmode=.false.
   failed=.false. ; epsit=1 

   do while(eps>mineps.and.xdiff>maxdiff)
    x1=x0
    nminh=0
    select case(dim)
     case(1)
      call propagate(x1,traj_time_start,traj_time_end,velocity,stepper)
      if(good_traj)xdiff=abs(x1(1)-xlast(1))
     case(2)
      call propagate(x1,traj_time_start,traj_time_end,velocity,stepper)
      if(good_traj)xdiff=sqrt((x1(1)-xlast(1))**2+(x1(2)-xlast(2))**2)
     case(3)
      call propagate(x1,traj_time_start,traj_time_end,velocity,stepper)
      if(good_traj)xdiff=sqrt((x1(1)-xlast(1))**2+(x1(2)-xlast(2))**2+&
       &(x1(3)-xlast(3))**2)
    end select

    if(plot_all_traj)then
     if(ntraj==1)open(14,file='trajectory'//trim(i2s(epsit))//'.dat',&
      &status='unknown')
     if(dim==1)then
      do n=1,nsteps
       write(14,"(2(1x,f25.16))")tp(n),xp(1,n)
      enddo
     else
      do n=1,nsteps
       write(14,"(3(1x,f25.16))")xp(1:dim,n)
      enddo
     endif
     if(ntraj==ntrajectories)close(14)
    endif

    if(good_traj)then
     xlast=x1
     if(epsit==1)then
      write(o,*)
      write(o,*)'POSSIBLE FINAL POSITIONS'
     endif
     write(o,*)
     tmpr=r2s(x1(1),'(f25.16)')
     tmpr2=r2s(eps,'(f15.12)')
     select case(dim)
      case(1)
       if(nsteps>0)then
        write(o,"(2x,a)")trim(tmpr)//'[EPS: '//trim(tmpr2)//&
         &' NSTEPS: '//trim(i2s(nsteps))//']'
       else
        write(o,"(2x,a)")trim(tmpr)//' '//trim(tmpr2)//' [EPS: '//trim(tmpr)//&
         &']'
       endif
      case(2)
       tmpr2=r2s(x1(2),'(f25.16)')
       tmpr3=r2s(eps,'(f15.12)')
       if(nsteps>0)then
        write(o,"(2x,a)")trim(tmpr)//' '//trim(tmpr2)//' [EPS: '//trim(tmpr3)//&
         &' NSTEPS: '//trim(i2s(nsteps))//']'
       else
        write(o,"(2x,a)")trim(tmpr)//' '//trim(tmpr2)//' [EPS: '//trim(tmpr3)//&
         &']'
       endif
      case(3)
       tmpr2=r2s(x1(2),'(f25.16)')
       tmpr3=r2s(x1(3),'(f25.16)')
       write(o,"(1x,3(1x,a))")trim(tmpr),trim(tmpr2),trim(tmpr3)
       tmpr=r2s(eps,'(f15.12)')
       if(nsteps>0)then
        write(o,'(2x,a)')'[EPS: '//trim(tmpr)//' NSTEPS: '//&
         &trim(i2s(nsteps))//']'
       else
        write(o,'(2x,a)')'[EPS: '//trim(tmpr)//']'
       endif
     end select
     failed=.false.
    else 
     write(o,*)
     write(o,'(2x,a,f15.12,a)')'FAILED TO COMPLETE PROPAGATION FOR EPS =',eps,&
      &'.'
     if(nminh/=0)then
      write(o,*)'STEPSIZE SMALLER THAN MINIMUM.'
     else
      write(o,'(2x,a)')'REACHED MAXIMUM NUMBER OF ITERATIONS MAXSTP = '//&
       &trim(i2s(maxstp))//'.' 
     endif
     failed=.true.
     if(fastmode)then
      if(.not.warned_fastmode)write(o,'(2x,a)')'FASTMODE ACTIVATED IN input - &
       &WOULD EXIT HERE IN DENSITY MODE.'
      warned_fastmode=.true.
     endif
     maxstp=min(maxstp*10,max_maxstp)
     tmpr=r2s(emergency_maxdiff,'(f6.4)')
     if(maxdiff/=emergency_maxdiff)write(o,'(2x,a)')'SWITCHING TO LESS&
      & TIGHT MAXDIFF TOLERANCE   = '//trim(tmpr)//'.'
     maxdiff=emergency_maxdiff

    endif

    eps=0.1d0*eps
    epsit=epsit+1
    save_steps=.true.

   enddo

   if(.not.failed)then
    write(o,*)
    if(eps<=mineps)then
     if(xdiff>maxdiff)then
      write(o,*)'REACHED MIN VALUE OF EPS WITHOUT CONVERGING.'
      write(o,*)'DISTANCE BETWEEN LAST TWO FINAL POSITIONS : ',xdiff
      write(o,*)'DISTANCE TOLERANCE                        : ',maxdiff
      failed=.true.
     else
      write(o,*)'CONVERGENCE ACHIEVED ON MINIMUM VALUE OF EPS.'
      write(o,*)'DISTANCE BETWEEN LAST TWO FINAL POSITIONS : ',xdiff
      write(o,*)'DISTANCE TOLERANCE                        : ',maxdiff
     endif
    else
     write(o,*)'CONVERGED'
     write(o,*)'DISTANCE BETWEEN LAST TWO FINAL POSITIONS : ',xdiff
     write(o,*)'DISTANCE TOLERANCE                        : ',maxdiff
    endif
   endif

   if(.not.plot_all_traj.and..not.failed)then
    if(ntraj==1)then
     open(14,file='trajectory.dat',status='unknown')
    else
     inquire(14,opened=op)
     if(.not.op)open(14,file='trajectory.dat',status='unknown')
     write(14,*)
    endif
    if(dim==1)then
     do n=1,nsteps
      write(14,"(2(1x,f25.16))")tp(n),xp(1,n)
     enddo
    else
     do n=1,nsteps
      write(14,"(3(1x,f25.16))")xp(1:dim,n)
     enddo
    endif
    if(ntraj==ntrajectories)then
     inquire(14,opened=op) ; if(op)close(14)
    endif
   endif

   if(.not.failed)then
    write(o,*)
    write(o,*)'FINAL POSITION'
    tmpr=r2s(x1(1),'(f25.16)')
    select case(dim)
     case(1)
      write(o,"(2x,a)")trim(tmpr)
     case(2)
      tmpr2=r2s(x1(2),'(f25.16)')
      write(o,"(1x,2(1x,a))")trim(tmpr),trim(tmpr2)
     case(3)
      tmpr2=r2s(x1(2),'(f25.16)')
      tmpr3=r2s(x1(3),'(f25.16)')
      write(o,"(1x,3(1x,a))")trim(tmpr),trim(tmpr2),trim(tmpr3)
    end select
   endif

   write(o,*)
   if(plot_all_traj)then

    if(.not.failed)then
     write(o,*)'PLOT_ALL_TRAJ keyword active.'
     if(ntraj==1)then
      write(o,*)"Final trajectories with different EPS tolerances written to &
       &'trajectoryX.dat'." 
     else
      write(o,*)"Final trajectories with different EPS tolerances appended to &
       &'trajectoryX.dat'." 
     endif
    else
     write(o,*)'PLOT_ALL_TRAJ keyword active.'
     if(ntraj==1)then
      write(o,*)"Final partial trajectories with different EPS &
       &written to 'trajectoryX.dat'"
     else
      write(o,*)"Final partial trajectories with different EPS &
       &appended to 'trajectoryX.dat'"
     endif
    endif

   else

    if(.not.failed)then
     if(ntraj==1)then
      write(o,*)"Final trajectory written to file 'trajectory.dat'."
     else
      write(o,*)"Final trajectory appended to file 'trajectory.dat'."
     endif
    else
     if(epsit>1)then
      if(ntraj==1)then
       write(o,*)"Final partial trajectory written to file 'trajectory.dat'."
      else
       write(o,*)"Final partial trajectory appended to file 'trajectory.dat'."
      endif
     else
      write(o,*)"Failed on first iteration - trajectory not saved. Set &
       &PLOT_ALL_TRAJ"
      write(o,*)'to T to plot the failed partial trajectory.'
     endif
    endif

   endif

   write(o,*)
   write(o,*)'Number of good steps                        : ',trim(i2s(nok))
   write(o,*)'Number of bad (but retried and fixed) steps : ',trim(i2s(nbad))

   eps=init_eps

  enddo ! ntrajectories

  if(ntrajectories>1)then
   write(o,*)
   write(o,*)'FINISHED ALL TRAJECTORIES'
  endif
  write(o,*)
  if(plot_all_traj)then
   if(dim==3)then
    write(o,*)'View trajectoryX.dat files with e.g. xmgrace.'
   else
    write(o,*)'View trajectoryX.dat files with e.g. plot_louis.'
   endif
  else
   if(dim==3)then
    write(o,*)'View trajectory.dat file with e.g. plot_louis.'
   else
    write(o,*)'View trajectory.dat file with e.g. xmgrace.'
   endif
  endif

  call TIMER('TRAJECTORY PROPAGATION',.false.)
 
  END SUBROUTINE trajectory_driver


  SUBROUTINE add_noise_to_phases
  IMPLICIT NONE
  REAL(dp),PARAMETER :: scale(14)=(/0.1d0,1.d-2,1.d-3,1.d-4,1.d-5,1.d-6,1.d-7,&
   &1.d-8,1.d-9,1.d-10,1.d-11,1.d-12,1.d-13,1.d-14/)
  REAL(dp),ALLOCATABLE,SAVE :: theta_1d_init(:),theta_2d_init(:,:),&
   &theta_3d_init(:,:,:)
  INTEGER ipn,jpn,kpn
  LOGICAL,SAVE :: first_call=.true.,print_noise=.true.

  if(first_call)then
   select case(dim)
    case(1)
     allocate(theta_1d_init(nmodes_dim),stat=ialloc)
     if(ialloc/=0)call errstop('ADD_NOISE_TO_PHASES','Phase allocation error.')
     theta_1d_init(:)=theta_1d(:)
    case(2)
     allocate(theta_2d_init(nmodes_dim,nmodes_dim),stat=ialloc)
     if(ialloc/=0)call errstop('ADD_NOISE_TO_PHASES','Phase allocation error.')
     theta_2d_init(:,:)=theta_2d(:,:)
    case(3)
     allocate(theta_3d_init(nmodes_dim,nmodes_dim,nmodes_dim),stat=ialloc)
     if(ialloc/=0)call errstop('ADD_NOISE_TO_PHASES','Phase allocation error.')
     theta_3d_init(:,:,:)=theta_3d(:,:,:)
   end select
   if(trim(adjustl(calc_type))=='density')print_noise=.false.
   first_call=.false.
  endif

  if(print_noise)then
   write(o,*)
   write(o,*)'Random noise added to phases.'
   write(o,*)'NEW PHASES'
  endif

  select case(dim)
   case(1)
    do ipn=1,nmodes_dim
     theta_1d(ipn)=theta_1d_init(ipn)+ranx_pm()*scale(phase_noise)
     if(print_noise)write(o,'(2x,a)')d2s0(theta_1d(ipn))
    enddo
   case(2)
    do jpn=1,nmodes_dim
     do ipn=1,nmodes_dim
      theta_2d(ipn,jpn)=theta_2d_init(ipn,jpn)+ranx_pm()*scale(phase_noise)
      if(print_noise)write(o,'(2x,a)')d2s0(theta_2d(ipn,jpn))
     enddo
    enddo
   case(3)
    do kpn=1,nmodes_dim
     do jpn=1,nmodes_dim
      do ipn=1,nmodes_dim
       theta_3d(ipn,jpn,k)=theta_3d_init(ipn,jpn,k)+ranx_pm()*scale(phase_noise)
       if(print_noise)write(o,'(2x,a)')d2s0(theta_3d(ipn,jpn,kpn))
      enddo
     enddo
    enddo

  end select
  if(print_noise)write(o,*)

  END SUBROUTINE add_noise_to_phases
 

  SUBROUTINE density_driver(velocity,stepper)
!-------------------------------------------------------------------------!
! Driver for density and H function calculations.                         !
!-------------------------------------------------------------------------!
   IMPLICIT NONE

   INTERFACE
   SUBROUTINE velocity(t,x,vel)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
   SUBROUTINE stepper(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
    REAL(dp),INTENT(inout) :: x(dim)
    REAL(dp),INTENT(in) :: htry,eps
    REAL(dp),INTENT(out) :: hdid,hnext
    REAL(dp),INTENT(inout) :: t
    INTERFACE
     SUBROUTINE velocity(t,x,vel)
      USE dsp
      USE store, ONLY : dim
      IMPLICIT NONE
      REAL(dp),INTENT(in) :: t
      REAL(dp),INTENT(in) :: x(dim)
      REAL(dp),INTENT(out) :: vel(dim)
     END SUBROUTINE velocity
    END INTERFACE
   END SUBROUTINE stepper
  END INTERFACE

  if(am_master)then
   write(o,*)
   write(o,*)'Density mode'
   write(o,*)'============'
   write(o,*)
   write(o,*)'Total no of points in the lattice : ',trim(i2s(nlattice**dim))
  endif

  save_steps=.false.

! Following 8 lines necessarily done in INITIAL_SETUP routine (for allocation)
! nlatticesq=nlattice*nlattice
! do i=1,num_ncgrain
!  ncgcell(i)=nlattice/ncgrain(i)  ! no of coarse-graining cells in the box
!  ncgcellsq(i)=ncgcell(i)**2
! enddo
! noscgcell=(nlattice-nscgrain)/nsmoothstep+1 ! no of overlapping smoothed 
!                                       ! coarse-graining cells in the box
! noscgcellsq=noscgcell*noscgcell
  one_over_nscgrain=1.d0/real(nscgrain,dp)
  nscgrainsq=nscgrain**2
  one_over_nscgrainsq=1.d0/real(nscgrainsq,dp)
  nscgrain3=nscgrain**3
  one_over_nscgrain3=1.d0/real(nscgrain3,dp)
  nscgrow=nlattice*nsmoothstep

  cell_over_nlat=cell_x/real(nlattice,dp)

! Begin the density calculation.
  if(den_ntimes==1)then
   den_timestep=0.d0 ! no propagation required to produce only 1 density
  else
   den_timestep=(den_time_end-den_time_start)/real(den_ntimes-1,dp)
  endif

  noktot=0_i64
  nbadtot=0_i64

  den_time=den_time_start

  select case(dim)
   case(1)
    call density_1d(velocity,stepper)
   case(2)
    call density_2d(velocity,stepper)
   case(3)
    call density_3d(velocity,stepper)
  end select

  if(.not.read_backtracked)then

   if(verbose.and.nnodes>1.and..not.(den_time_start==0.d0.and.den_ntimes==1))&
    &then
    if(am_master)then
     write(o,*)
     write(o,*)'------------VERBOSE------------'
     if(ialloc/=0)call errstop('DENSITY_DRIVER','Cannot allocate NTEMP.')
    endif
    allocate(ntemp(1:nnodes),stat=ialloc) ! on all nodes because of compiler
    call mpi_gather(noktot,1,mpi_integer8,ntemp,1,mpi_integer8,0,&
     &mpi_comm_world,ierr)
    call checkmpi(ierr,'Error in MPI_GATHER operation for noktot counter.')
    if(am_master)then
     write(o,*)'Total number of good +h steps:'
     do n=1,nnodes
      write(o,*)'Node ',trim(i2s(n)),': ',trim(i2s64(ntemp(n)))
     enddo
    endif
    call mpi_gather(nbadtot,1,mpi_integer8,ntemp,1,mpi_integer8,0,&
     &mpi_comm_world,ierr)
    call checkmpi(ierr,'Error in MPI_GATHER operation for nbadtot counter.')
    if(am_master)then
     write(o,*)'Total number of bad (but retried and fixed) steps:'
     do n=1,nnodes
      write(o,*)'Node ',trim(i2s(n)),': ',trim(i2s64(ntemp(n)))
     enddo
     deallocate(ntemp)
     write(o,*)'----------END-VERBOSE----------'
    endif
   endif
  
   noktot_total=0_i64 ; nbadtot_total=0_i64
   call mpi_reduce(noktot,noktot_total,1,mpi_integer8,mpi_sum,0,&
     &mpi_comm_world,ierr)
   call checkmpi(ierr,'Error in MPI_REDUCE operation for noktot counter.')
   call mpi_reduce(nbadtot,nbadtot_total,1,mpi_integer8,mpi_sum,0,&
     &mpi_comm_world,ierr)
   call checkmpi(ierr,'Error in MPI_REDUCE operation for nbadtot counter.')
   if(nstep_histogram)then
    if(am_master)then
     allocate(nhist_total(nbins),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_DRIVER','Unable to allocate&
      & NHIST_TOTAL vector.')
     nhist_total(:)=0
    endif
    call mpi_reduce(nhist,nhist_total,nbins,mpi_integer,&
     &mpi_sum,0,mpi_comm_world,ierr)
    call checkmpi(ierr,'Error in MPI_REDUCE operation for nhist counter.')
   endif

  endif

  if(am_master)then

   if(den_ntimes>1)then

    write(o,*)
    write(o,*)'All timesteps done.'

    if(hfunction)then

     do j=1,num_ncgrain
      if(den_ntimes>2)then
       write(o,*)
       t1=(real(ncgrain(j),dp)*cell_x)/real(nlattice,dp)
       tmpr=r2s(t1,'(f12.6)')
       if(num_ncgrain>1)write(o,*)'COARSE-GRAINING-LENGTH '//&
        &trim(i2s(ncgrain(j))) //' ('//trim(tmpr)//')'
      endif
      tmpr=''
      if(num_ncgrain>1)tmpr='_cg='//trim(i2s(ncgrain(j)))
      filename='H_vs_t'//trim(tmpr)//'.dat'
      open(file=filename,unit=14,status='unknown')
      den_time=den_time_start
      n=0
      do i=1,den_ntimes
       write(14,'(2f25.16)')den_time,H(j,i)
       if(H(j,i)>0.d0)then
        n=n+1
        ttt(n)=den_time ; lnH(n)=log(H(j,i)) 
       endif
       den_time=den_time+den_timestep
      enddo
      close(14)
      if(den_ntimes>2)write(o,*)'Written time-dependent H function.'
      filename='lnH_vs_t'//trim(tmpr)//'.dat'
      open(file=filename,unit=14,status='unknown')
      do i=1,den_ntimes
       if(H(j,i)>0.d0)write(14,'(2f25.16)')ttt(i),lnH(i)
      enddo
      if(den_ntimes>2)then
       write(o,*)'Written time-dependent log of H function.'
       call straight_line_fit(ttt(1:n),lnH(1:n),n,fita,fitb,ua,ub,chisq)
       write(o,'(1x,a)')'Fit lnH vs t to straight line y=a+bt.'
       tmpr=r2s(fita,'(f12.6)') ; tmpr2=r2s(ua,'(f12.6)')
       write(o,'(1x,5a)')'Best a parameter: ',trim(tmpr),' with uncertainty ',&
        &trim(tmpr2),'.'
       tmpr=r2s(fitb,'(f12.6)') ; tmpr2=r2s(ub,'(f12.6)')
       write(o,'(1x,5a)')'Best b parameter: ',trim(tmpr),' with uncertainty ',&
        &trim(tmpr2),'.'
       tmpr=r2s(chisq,'(f12.6)')
       write(o,'(1x,2a)')'Value of chi^2  :  ',trim(tmpr)
       tmpr=r2s(-1.d0/fitb,'(f12.6)')
       write(o,'(1x,2a)')'Assuming H(t)=H(0)exp(-t/t_c) then RELAXATION TIME&
        & t_c = ',trim(tmpr)
       write(14,*)
       write(14,*)'# Straight-line fit'
       do i=1,den_ntimes
        write(14,*)ttt(i),fita+fitb*ttt(i)
       enddo
      endif
      close(14)

     enddo ! num_ncgrain

     if(den_ntimes==1)then
      write(o,*)'Written time-dependent H function.'
      write(o,*)'Written time-dependent log of H function.'
      write(o,*)'Not enough data for fit to estimate relaxation time&
       & (increase DEN_NTIMES).'
     endif
     if(den_ntimes==2)then
      write(o,*)'Written time-dependent H functions.'
      write(o,*)'Written time-dependent log of H functions.'
      write(o,*)'Not enough data for fit to estimate relaxation time&
       & (increase DEN_NTIMES).'
     endif
 
    endif ! hfunction

   endif ! den_ntimes>1

   if(.not.read_backtracked)then
    write(o,*)
    if(.not.(den_ntimes==1.and.den_time_start==0.d0))then
     write(o,*)'Total number of good +h steps                        : ',&
     &trim(i2s64(noktot_total))
     write(o,*)'Total number of bad (but retried and fixed) +h steps : ',&
      &trim(i2s64(nbadtot_total))
     if(noktot_total==0.and.nbadtot_total==0)call errstop('DENSITY_DRIVER',&
      &'Error working out percentage of good steps. Should not happen.')
     tmpr=r2s(100.d0*real(noktot_total,dp)/real(noktot_total+nbadtot_total,dp),&
      &'(f12.2)')
     write(o,'(1x,a)')'Percentage of good +h steps                          : '&
      &//trim(tmpr)//'%'
     tmpr=r2s(real(noktot_total+nbadtot_total,dp)/real(nlattice,dp)/&
      &real(den_ntimes,dp),'(f25.2)')
     write(o,'(1x,a)')'Average number of steps to complete one trajectory&
      &   : '//trim(tmpr)
     if(nstep_histogram)then
      open(file='nstep_histogram.dat',unit=22,status='unknown')
      do i=1,nbins
       write(22,*)int(real(i-1,dp)/nhist_scale),nhist_total(i)
      enddo
      write(o,*)'Written histogram of number of steps per trajectory.'
      write(o,*)
     endif
     write(o,*)'Number of trajectories for which backtracking failed : ',&
      &trim(i2s(sum(nfail_total)))
     if(sum(nfail_total)>0)then
      if(nlattice**dim==0)call errstop('DENSITY_DRIVER','Error working out &
       &percentage of failed trajectories. Should not happen.')
      t1=100.d0*real(sum(nfail_total),dp)/real((nlattice**dim),dp)
      if(t1>1.d0)then
       tmpr=r2s(t1,'(f12.2)')
      else
       if(t1>0.1d0)then
        tmpr=r2s(t1,'(f12.3)')
       else
        if(t1>0.01d0)then
         tmpr=r2s(t1,'(f12.4)')
        else
         if(t1>0.001d0)then
          tmpr=r2s(t1,'(f12.5)')
         else
          if(t1>0.001d0)then
           tmpr=r2s(t1,'(f12.5)')
          else
           tmpr=r2s(t1,'(f12.8)')
          endif
         endif
        endif
       endif
      endif
      write(o,'(1x,a)')'Percentage failed points                             :&
       & '//trim(tmpr)//'%'
     endif
    endif
   endif

  endif ! am_master
 
  END SUBROUTINE density_driver


  SUBROUTINE propagate(xx,time1,time2,velocity,stepper)
  IMPLICIT NONE
  REAL(dp),INTENT(inout) :: xx(dim)
  REAL(dp),INTENT(in) :: time1,time2
  REAL(dp),PARAMETER :: tiny=1.d-30

  INTERFACE
   SUBROUTINE velocity(t,x,vel)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
   SUBROUTINE stepper(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
    REAL(dp),INTENT(inout) :: x(dim)
    REAL(dp),INTENT(in) :: htry,eps
    REAL(dp),INTENT(out) :: hdid,hnext
    REAL(dp),INTENT(inout) :: t
    INTERFACE
     SUBROUTINE velocity(t,x,vel)
      USE dsp
      USE store, ONLY : dim
      IMPLICIT NONE
      REAL(dp),INTENT(in) :: t
      REAL(dp),INTENT(in) :: x(dim)
      REAL(dp),INTENT(out) :: vel(dim)
     END SUBROUTINE velocity
    END INTERFACE
   END SUBROUTINE stepper
  END INTERFACE

  INTEGER nstp
  REAL(dp) h,hdid,hnext,t,tsav
  REAL(dp),DIMENSION(dim) :: vel,x,xscal

  t=time1
  h=sign(h1,time2-time1)
  nok=0
  nbad=0
  nminhfail=.false.
  nsteps=0
  x(:)=xx(:)

  if(save_steps)then
   tsav=t-2.d0*dtsave ! assures storage of first step
   nullify(tp,xp)
   allocate(tp(256))
   allocate(xp(size(xx),size(tp)))
  endif

  do nstp=1,maxstp ! take at most maxstp steps

   call velocity(t,x,vel)
   xscal(:)=abs(x(:))+abs(h*vel(:))+tiny

   if(save_steps.and.(abs(t-tsav)>abs(dtsave)))then
    nsteps=nsteps+1
    if(nsteps>size(tp))then
     tp=>reallocate(tp,2*size(tp))
     xp=>reallocate(xp,size(xp,1),size(tp))
    endif
    tp(nsteps)=t
    xp(:,nsteps)=x(:)
    tsav=t
   endif

   if((t+h-time2)*(t+h-time1)>0.d0)h=time2-t!if stepsize can overshoot, decrease

   call stepper(x,vel,t,h,eps,xscal,hdid,hnext,velocity)

   if(hdid==h)then
    nok=nok+1
   else
    nbad=nbad+1
   endif

   if((t-time2)*(time2-time1)>=0.d0)then ! Are we done?
    xx(:)=x(:)

    if(save_steps)then
     nsteps=nsteps+1
     if(nsteps>size(tp))then
      tp=>reallocate(tp,2*size(tp))
      xp=>reallocate(xp,size(xp,1),size(tp))
     endif
     tp(nsteps)=t
     xp(:,nsteps)=x(:)
     tsav=t
    endif

    good_traj=.true.
    return ! normal exit
   endif

   if(abs(hnext)<hmin)then ! stepsize smaller than minimum
    nminhfail=.true. ; good_traj=.false. ; exit
   endif

   h=hnext

  enddo ! nstp

  good_traj=.false.
  xx(:)=x(:)

  END SUBROUTINE propagate


  SUBROUTINE density_1d(velocity,stepper)
!-------------------------------------------------------------------------!
! Calculate the raw, coarse-grained, and smoothed 1D densities and psi^2. !
!-------------------------------------------------------------------------!
  IMPLICIT NONE

  INTERFACE
   SUBROUTINE velocity(t,x,vel)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
   SUBROUTINE stepper(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
    REAL(dp),INTENT(inout) :: x(dim)
    REAL(dp),INTENT(in) :: htry,eps
    REAL(dp),INTENT(out) :: hdid,hnext
    REAL(dp),INTENT(inout) :: t
    INTERFACE
     SUBROUTINE velocity(t,x,vel)
      USE dsp
      USE store, ONLY : dim
      IMPLICIT NONE
      REAL(dp),INTENT(in) :: t
      REAL(dp),INTENT(in) :: x(dim)
      REAL(dp),INTENT(out) :: vel(dim)
     END SUBROUTINE velocity
    END INTERFACE
   END SUBROUTINE stepper
  END INTERFACE

! Generate 1d raw lattice.
  do i=1,nlattice
   latpos1(i)=(real(i,dp)-0.5d0)*cell_over_nlat
  enddo

  if(plot_smooth)then
! Generate 1d smoothing lattice.
   n=nscgrain/2
   do i=1,noscgcell
    latpos1_smooth(i)=(real(n,dp)-0.5d0)*cell_over_nlat
    n=n+nsmoothstep
   enddo
  endif

  if(read_backtracked)writefail=.false.

  if(writefail)then
   if(nnodes>1)then
    open(file='failed_trajectories'//trim(i2s(my_node))//'.dat',unit=21,&
     &status='unknown')
    write(21,*)'# NODE ',trim(i2s(my_node))
   else
    open(file='failed_trajectories.dat',unit=21,status='unknown')
   endif
   write(21,*)'# START TIME, START POSITION OF FAILED TRAJECTORIES'
  endif

  nfail_total(:)=0

! Loop over time sequence.
  do nplot=1,den_ntimes

   if(am_master)then
    if(den_time==0.d0)then
     tmpr='0.0'
     tmpr2='0.0'
    else
     tmpr=r2s(den_time,'(f25.16)')
     tmpr2=r2s(den_time,'(f12.2)')
    endif
   endif

   if(.not.read_backtracked)then

    call timer('BACKTRACK',.true.)
 
    if(am_master)then
     write(o,*)
     write(o,*)'Computing data at t = ',trim(adjustl(tmpr)),'.'
     if(den_time/=0.d0)then
      write(o,*)'Backtracking trajectories to zero for each point.'
     else
      write(o,*)'No backtracking required at time zero.'
     endif
    endif

    n=int((den_time)/(fourpi+1.d-6))+1
    maxstp_in=n*maxstp_in
    max_maxstp=n*max_maxstp
 
    nfail=0 ; nzero=0 ; nzero_total=0 ; nminh=0 ; nminh_total=0 

    n=my_node+1

    if(den_time==0.d0)then ! no propagation required

     do i=1,is_per_node(my_node) ! i.e. i=1,nlattice on 1 node
      x0(1)=latpos1(n)
      call eval_psi_squared_t0_1d(x0,rho)
      psi_squared_raw_n(i)=rho
      call eval_density_t0_1d(dentype,x0,rho)
      density_raw_n(i)=rho
      n=n+nnodes
     enddo

    else ! t > 0.0

     do i=1,is_per_node(my_node) ! i.e. i=1,nlattice on 1 node
      x0(1)=latpos1(n)
      if(phase_noise>0)call add_noise_to_phases
! psi^2 at lattice x, t
      call eval_psi_squared_1d(den_time,x0,psi2_t)
      if(psi2_t==0.d0)then
       nzero=nzero+1 ; good_traj=.false.
      else
       psi_squared_raw_n(i)=psi2_t
! Find backtracked positions.
       eps=init_eps ; xlast=1.d10 ; xdiff=1.d10 
       maxdiff=converge_maxdiff ; maxstp=maxstp_in
       do while(eps>mineps.and.xdiff>maxdiff)
        x1=x0
        call propagate(x1,den_time,0.d0,velocity,stepper)
        if(good_traj)then
         xdiff=abs(x1(1)-xlast(1))
         xlast=x1
        else
         if(fastmode)then
          exit
         else
          maxstp=min(maxstp*10,max_maxstp)
          maxdiff=emergency_maxdiff
         endif
        endif
        eps=0.1d0*eps
       enddo
       noktot=noktot+nok
       nbadtot=nbadtot+nbad
       if(nstep_histogram)then
        ii=int(real(nok+nbad,dp)*nhist_scale)
        nhist(ii+1)=nhist(ii+1)+1
       endif
      endif
      if(good_traj)then
       if(save_backtracked)xback1_n(i)=x1(1)
! rho at backtracked x, t0
       call eval_density_t0_1d(dentype,x1,rho_0)
! psi^2 at backtracked x, t0
       call eval_psi_squared_t0_1d(x1,psi2_0)
! rho at lattice x, t
       density_raw_n(i)=rho_0*psi2_t/psi2_0
      else
       nfail=nfail+1
       if(nminhfail)nminh=nminh+1
       if(writefail)write(21,*)den_time,latpos1(n),' 0.d0 0.d0'
       latpos1(n)=-latpos1(n)
       if(save_backtracked)xback1_n(i)=latpos1(n)
       psi_squared_raw_n(i)=0.d0
       density_raw_n(i)=0.d0
      endif
      n=n+nnodes ! load balancing
     enddo

    endif ! time > 0.d0 or not

    call timer('BACKTRACK',.false.)

    if(nnodes>1)then
     call timer('BARRIER',.true.)
     call barrier
     call timer('BARRIER',.false.)
    endif
  
    if(am_master)then
     allocate(density_raw(nlattice),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_1D','Insufficient memory to &
      &allocate raw density vector.')
     allocate(psi_squared_raw(nlattice),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_1D','Insufficient memory to &
      &allocate raw psi_squared vector.')
     if(save_backtracked)then
      allocate(xback1(nlattice),stat=ialloc)
      if(ialloc/=0)call errstop_master('DENSITY_1D','Insufficient memory to &
       &allocate XBACK vector.')
     endif   
    endif

    if(nnodes>1)then
     call timer('MPI',.true.)
     n=is_per_node(0)
     call mpi_type_vector(n,1,nnodes,mpi_double_precision,mdttype1,ierr) 
     call mpi_type_commit(mdttype1,ierr)

! mpi_type_create_resized followed by mpi_gather or mpi_gatherv doesn't appear 
! to work as advertised for interlaced receive vectors..(?) Use hideous set of 
! send/receives instead.
!  call mpi_type_create_resized(mdttype1,lb,extent,mdttype2,ierr)
!  call mpi_type_commit(mdttype2,ierr)
!  call mpi_gather(density_raw_n,n,mpi_double_precision,density_raw,1,mdttype2,&
!   &0,mpi_comm_world,ierr)
!  call mpi_gather(psi_squared_raw_n,n,mpi_double_precision,psi_squared_raw,&
!   &1,mdttype2,0,mpi_comm_world,ierr)
!  call mpi_type_free(mdttype2,ierr)

     if(nextra>0)then
      call mpi_type_vector(n-1,1,nnodes,mpi_double_precision,mdttype2,ierr) 
      call mpi_type_commit(mdttype2,ierr)
     endif
     if(am_slave)then
      if(nextra==0.or.my_node<nextra)then
       call mpi_ssend(density_raw_n,n,mpi_double_precision,0,my_node,&
        &mpi_comm_world,ierr)
      else
       call mpi_ssend(density_raw_n,n-1,mpi_double_precision,0,my_node,&
        &mpi_comm_world,ierr)
      endif
      call checkmpi(ierr,'MPI_SSEND error on slave <1>')
     else ! am_master
      j=1
      do i=1,n
       density_raw(j)=density_raw_n(i) ; j=j+nnodes
      enddo
      if(nextra==0)then
       do i=1,nnodes-1
        call mpi_recv(density_raw(i+1),1,mdttype1,i,i,mpi_comm_world,status,&
         &ierr)
       enddo
      else
       do i=1,nextra-1
        call mpi_recv(density_raw(i+1),1,mdttype1,i,i,mpi_comm_world,status,&
         &ierr)
       enddo
       do i=nextra,nnodes-1
        call mpi_recv(density_raw(i+1),1,mdttype2,i,i,mpi_comm_world,status,&
         &ierr)
       enddo
      endif
      call checkmpi(ierr,'MPI_RECV error on master <1>')
     endif ! master/slave
     if(am_slave)then
      if(nextra==0.or.my_node<nextra)then
       call mpi_ssend(psi_squared_raw_n,n,mpi_double_precision,0,my_node,&
        &mpi_comm_world,ierr)
      else 
       call mpi_ssend(psi_squared_raw_n,n-1,mpi_double_precision,0,my_node,&
        &mpi_comm_world,ierr)
      endif
      call checkmpi(ierr,'MPI_SSEND error on slave <2>')
     else ! am_master
      j=1
      do i=1,n
       psi_squared_raw(j)=psi_squared_raw_n(i) ; j=j+nnodes
      enddo
      if(nextra==0)then
       do i=1,nnodes-1
        call mpi_recv(psi_squared_raw(i+1),1,mdttype1,i,i,mpi_comm_world,&
         &status,ierr) 
       enddo
      else
       do i=1,nextra-1
        call mpi_recv(psi_squared_raw(i+1),1,mdttype1,i,i,mpi_comm_world,&
         &status,ierr) 
       enddo
       do i=nextra,nnodes-1
        call mpi_recv(psi_squared_raw(i+1),1,mdttype2,i,i,mpi_comm_world,&
         &status,ierr) 
       enddo
      endif
      call checkmpi(ierr,'MPI_RECV error on master <2>')
     endif ! master/slave
    if(save_backtracked)then
     if(am_slave)then
      if(nextra==0.or.my_node<nextra)then
        call mpi_ssend(xback1_n,n,mpi_double_precision,0,my_node,&
         &mpi_comm_world,ierr)
       else
        call mpi_ssend(xback1_n,n-1,mpi_double_precision,0,my_node,&
         &mpi_comm_world,ierr)
       endif
       call checkmpi(ierr,'MPI_SSEND error on slave <3>')
      else ! am_master
       j=1
       do i=1,n
        xback1(j)=xback1_n(i) ; j=j+nnodes
       enddo
       if(nextra==0)then
        do i=1,nnodes-1
         call mpi_recv(xback1(i+1),1,mdttype1,i,i,mpi_comm_world,status,ierr)
        enddo
       else
        do i=1,nextra-1
         call mpi_recv(xback1(i+1),1,mdttype1,i,i,mpi_comm_world,status,ierr)
        enddo
        do i=nextra,nnodes-1
         call mpi_recv(xback1(i+1),1,mdttype2,i,i,mpi_comm_world,status,ierr)
        enddo
       endif
       call checkmpi(ierr,'MPI_RECV error on master <3>')
      endif ! master/slave
     endif ! save_backtracked
 
     call mpi_type_free(mdttype1,ierr)
     if(nextra>0)call mpi_type_free(mdttype2,ierr)

     call timer('MPI',.false.)
 
    else ! nnodes==1

     density_raw(:)=density_raw_n(:)
     psi_squared_raw(:)=psi_squared_raw_n(:)
     if(save_backtracked)xback1(:)=xback1_n(:)

    endif

    call mpi_reduce(nfail,nfail_total(nplot),1,mpi_integer,mpi_sum,0,&
     &mpi_comm_world,ierr)
    call mpi_reduce(nzero,nzero_total,1,mpi_integer,mpi_sum,0,&
     &mpi_comm_world,ierr)
    call mpi_reduce(nminh,nminh_total,1,mpi_integer,mpi_sum,0,&
     &mpi_comm_world,ierr)

    if(am_master)then
     write(o,*)'Done.'
     if(den_time/=0.d0)then
      if(nfail_total(nplot)>0)then
       write(o,'(1x,a)')'Number of failed trajectories          : '&
        &//trim(i2s(nfail_total(nplot)))
       write(o,'(1x,a)')'Number failed for too many steps       : '&
        &//trim(i2s(nfail_total(nplot)-nzero_total-nminh_total))
       write(o,'(1x,a)')'Number failed for too small stepsize   : '&
        &//trim(i2s(nminh_total))
       write(o,'(1x,a)')'Number failed because started on nodes : '&
        &//trim(i2s(nzero_total))
      else
       write(o,'(1x,a)')'Number of failed trajectories: '//&
        &trim(i2s(nfail_total(nplot)))
      endif
     endif
     if(plot_raw.or.plot_cg.or.plot_smooth.or.plot_h_integrand)then
      write(o,*)'Writing requested data to disk.'
     endif
    endif
 
   else ! read_backtracked=T

    if(am_master)then

     allocate(density_raw(nlattice),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_1D','Insufficient memory to &
      &allocate raw density vector.')
     allocate(psi_squared_raw(nlattice),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_1D','Insufficient memory to &
      &allocate raw psi_squared vector.')
     allocate(xback1(nlattice),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_1D','Insufficient memory to &
      &allocate XBACK vector.')

     write(o,*)
     write(o,*)'Computing data at t = ',trim(adjustl(tmpr)),'.'
     if(den_time/=0.d0)then
      write(o,*)'Reading previously backtracked positions from disk.'
     else
      write(o,*)'No backtracking required at time zero.'
     endif

     if(den_time==0.d0)then

      do i=1,nlattice
       x0(1)=latpos1(i)
       call eval_psi_squared_t0_1d(x0,rho)
       psi_squared_raw(i)=rho
       call eval_density_t0_1d(dentype,x0,rho)
       density_raw(i)=rho
      enddo

     else ! den_time/=0.d0

      filename='backtracked_positions_t='//trim(adjustl(tmpr2))//'.dat'
      inquire(file=filename,exist=file_present)
      if(.not.file_present)then
       write(o,'(1x,a)')'Required filename: '//filename
       call errstop('DENSITY_1D','READ_BACKTRACKED=T but&
       & file containing backtracked positions is missing.')
      endif
      open(file=filename,unit=14,status='old')
      read(14,*,iostat=ierr)
      if(ierr/=0)call errstop('DENSITY_1D','Error reading backtracked&
       & positions file.')
      read(14,*,iostat=ierr)t1
      if(ierr/=0)call errstop('DENSITY_1D','Error reading backtracked&
       & positions file (second line).')
      if(abs(t1-den_time)>1.d-6)call errstop('DENSITY_1D','File containing&
       & backtracked positions has wrong timestamp.')
      do i=1,nlattice
       read(14,*,err=1,end=2)xback1(i)
      enddo
      write(o,*)'Done.'
      write(o,*)
      do i=1,nlattice
       if(xback1(i)>0.d0)then
        call eval_psi_squared_1d(den_time,latpos1(i),psi2_t)
        psi_squared_raw(i)=psi2_t
        call eval_density_t0_1d(dentype,xback1(i),rho_0)
        call eval_psi_squared_t0_1d(xback1(i),psi2_0)
        density_raw(i)=rho_0*psi2_t/psi2_0
       else
        nfail_total(nplot)=nfail_total(nplot)+1
        latpos1(i)=-latpos1(i)
        psi_squared_raw(i)=0.d0
        density_raw(i)=0.d0
       endif
      enddo

     endif ! den_time==0.d0 or not

    endif ! am_master

   endif ! read_backtracked or not

   if(am_master)then

! Write out backtracked positions, if requested.
    if(save_backtracked.and.den_time/=0.d0)then

     call timer('RAW WRITE',.true.)

     filename='backtracked_positions_t='//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,*)'BACKTRACKED POSITIONS from time on second line.'
     write(14,*)den_time
     do i=1,nlattice
      write(14,'(f25.16)')xback1(i)
     enddo
     close(14)
     write(o,*)'Written backtracked positions.'

     call timer('RAW WRITE',.false.)

    endif ! save_backtracked

! Write out raw densities and psi^2 if requested.
    if(plot_raw)then

     call timer('RAW WRITE',.true.)

     filename='density_raw_t='//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
     do i=1,nlattice
      write(14,'(2f25.16)')abs(latpos1(i)),density_raw(i)
     enddo
     close(14)
     write(o,*)'Written raw density.'

     filename='psisq_raw_t='//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
     do i=1,nlattice
      write(14,'(2f25.16)')abs(latpos1(i)),psi_squared_raw(i)
     enddo
     close(14)
     write(o,*)'Written raw psi squared.'

     call timer('RAW WRITE',.false.)

    endif ! plot_raw

! Compute and write out non-overlapping coarse-grained densities and
! psi^2 if requested.
    if(plot_cg.or.hfunction)then

     call timer('COARSE-GRAINING',.true.)
 
     do ncg=1,num_ncgrain 

      nc=ncgcell(ncg) ; ng=ncgrain(ncg)
      one_over_ncgrain=1.d0/real(ng,dp)
      interpolation_required=.false.

! Generate 1d coarse-grained lattice.
      t1=real(ng,dp)*cell_over_nlat
      do i=1,nc
       latpos1_cg(i)=(real(i,dp)-0.5d0)*t1
      enddo

      allocate(density_cg(nc),psi_squared_cg(nc),stat=ialloc)
      if(ialloc/=0)call errstop_master('DENSITY_1D','Insufficient memory to &
       &allocate coarse-graining vectors.')

      if(nfail_total(nplot)==0)then
       n=1
       do i=1,nc
        t1=0.d0 ; t2=0.d0
        do ii=1,ng
         t1=t1+density_raw(n)
         t2=t2+psi_squared_raw(n)
         n=n+1
        enddo
        density_cg(i)=t1
        psi_squared_cg(i)=t2
       enddo
       density_cg(:)=density_cg(:)*one_over_ncgrain
       psi_squared_cg(:)=psi_squared_cg(:)*one_over_ncgrain
      else ! nfail_total > 0
       n=1
       do i=1,nc
        ngoodcg=ng ; t1=0.d0 ; t2=0.d0
        do ii=1,ng
         if(latpos1(n)<0.d0)then
          ngoodcg=ngoodcg-1
         else
          t1=t1+density_raw(n)
          t2=t2+psi_squared_raw(n)
         endif
         n=n+1
        enddo
        if(ngoodcg==0)then
         density_cg(i)=-1.d0 ; psi_squared_cg(i)=-1.d0
         interpolation_required=.true.
        else
         density_cg(i)=t1/real(ngoodcg,dp)
         psi_squared_cg(i)=t2/real(ngoodcg,dp)
        endif
       enddo
      endif ! nfail_total

! Interpolate to repair (very unlikely in 1D) gaps in coarse-grained density 
! and psi_squared due to failed trajectories.
      if(interpolation_required)then
       k=1
       do i=1,nc
        if(density_cg(k)==-1.d0)then
         ii=max(1,k-5) ; jj=min(nc,k+5)
         n=jj-ii ; nn=-1 ; nnn=0
         do kk=1,n
          if(ii+kk-1==k)nn=0
          if(density_cg(ii+kk+nn)/=-1.d0)then
           nnn=nnn+1
           points(nnn)=latpos1_cg(ii+kk+nn)
           repair(nnn)=density_cg(ii+kk+nn)
          endif
         enddo
         if(nnn>1)then
          call interp_nev(points(1:nnn),repair(1:nnn),nnn,latpos1_cg(k),&
           &val_interp,err_interp)
          density_cg(k)=val_interp
         endif
        endif
        if(density_cg(k)<0.d0)density_cg(k)=0.d0
        k=k+1
       enddo
       k=1
       do i=1,nc
        if(psi_squared_cg(k)==-1.d0)then
         ii=max(1,k-5) ; jj=min(nc,k+5)
         n=jj-ii ; nn=-1 ; nnn=0
         do kk=1,n
          if(ii+kk-1==k)nn=0
          if(psi_squared_cg(ii+kk+nn)/=-1.d0)then
           nnn=nnn+1
           points(nnn)=latpos1_cg(ii+kk+nn)
           repair(nnn)=psi_squared_cg(ii+kk+nn)
          endif
         enddo
         if(nnn>1)then
          call interp_nev(points(1:nnn),repair(1:nnn),nnn,latpos1_cg(k),&
           &val_interp,err_interp)
          psi_squared_cg(k)=val_interp
         endif
        endif
        if(psi_squared_cg(k)<=0.d0)psi_squared_cg(k)=1.d-5
        k=k+1
       enddo
      endif ! interpolation_required

! Write coarse-grained quantities.
      if(plot_cg)then

       tmpr=''
       if(num_ncgrain>1)tmpr='_cg='//trim(i2s(ng))

       filename='density_cg_t='//trim(adjustl(tmpr2))//trim(tmpr)//'.dat'
       open(file=filename,unit=14,status='unknown')
       write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
       do i=1,nc
        write(14,'(2f25.16)')latpos1_cg(i),density_cg(i)
       enddo
       close(14)
       if(num_ncgrain==1)write(o,*)'Written coarse-grained density.'

       filename='psisq_cg_t='//trim(adjustl(tmpr2))//trim(tmpr)//'.dat'
       open(file=filename,unit=14,status='unknown')
       write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
       do i=1,nc
        write(14,'(2f25.16)')latpos1_cg(i),psi_squared_cg(i)
       enddo
       close(14)
       if(num_ncgrain==1)write(o,*)'Written coarse-grained psi squared.'

      endif

! Calculate H function.
      if(hfunction.or.plot_h_integrand)then
        do i=1,nc
        if(density_cg(i)/=0.d0.and.psi_squared_cg(i)/=0.d0)then
         t1=log(density_cg(i)/psi_squared_cg(i)) 
         density_cg(i)=density_cg(i)*t1 
        else
         density_cg(i)=0.d0 ! limit of x ln ax = 0 as x-->0+
        endif
       enddo
 
       if(plot_h_integrand)then
        filename='h_integrand_t='//trim(adjustl(tmpr2))//trim(tmpr)//'.dat'
        open(file=filename,unit=14,status='unknown')
        write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
        write(14,'(2f25.16)')0.d0,0.d0
        do i=1,nc
         write(14,'(2f25.16)')latpos1_cg(i),density_cg(i)
        enddo
        write(14,'(2f25.16)')cell_x,0.d0
        close(14)
       endif
 
       if(hfunction)then
        if(nc<6)call errstop('DENSITY_1D','Not enough coarse-graining&
         & cells to calculate the H function - need at least 6.')
        t1=latpos1_cg(1)*2.d0
        H(ncg,nplot)=0.5d0*t1*(density_cg(1)+density_cg(nc)) ! linear to 0
        density_cg(1)=density_cg(1)*0.375d0   
        density_cg(2)=density_cg(2)*seven_sixths
        density_cg(3)=density_cg(3)*twentythree_over_twentyfour
        density_cg(nc-2)=density_cg(nc-2)*twentythree_over_twentyfour
        density_cg(nc-1)=density_cg(nc-1)*seven_sixths
        density_cg(nc)=density_cg(nc)*0.375d0
        H(ncg,nplot)=H(ncg,nplot)+sum(density_cg(:))*t1
       endif

      endif ! hfunction or plot_h_integrand

      deallocate(density_cg,psi_squared_cg)

     enddo ! num_ncgrain
     if(num_ncgrain>1)write(o,*)'Written coarse-grained densities.'
     if(num_ncgrain>1)write(o,*)'Written coarse-grained psi-squareds.'

     call timer('COARSE-GRAINING',.false.)

    endif ! plot_cg or hfunction

! Compute and write out smoothed densities and psi^2 if requested.
    if(plot_smooth)then
     call timer('SMOOTHING',.true.)

     allocate(density_smooth(noscgcell),psi_squared_smooth(noscgcell),&
      &stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_1D','Insufficient memory to &
      &allocate smoothing vectors.')

     nn=0
     if(nfail_total(nplot)==0)then
      do i=1,noscgcell
       t1=0.d0 ; n=nn*nsmoothstep+1
       do ii=1,nscgrain
        t1=t1+density_raw(n)
        n=n+1
       enddo
       density_smooth(i)=t1
       nn=nn+1
      enddo
      density_smooth(:)=density_smooth(:)*one_over_nscgrain
     else ! nfail_total > 0
      do i=1,noscgcell
       ngoodcg=nscgrain ; t1=0.d0 ; n=nn*nsmoothstep+1
       do ii=1,nscgrain
        if(latpos1(n)<0.d0)then
         ngoodcg=ngoodcg-1
        else
         t1=t1+density_raw(n)
        endif
        n=n+1
       enddo
       if(ngoodcg==0)then
        density_smooth(i)=0.d0
       else
        density_smooth(i)=t1/real(ngoodcg,dp)
       endif
       nn=nn+1
      enddo
     endif ! nfail_total

     nn=0
     do i=1,noscgcell
      t1=0.d0 ; n=nn*nsmoothstep+1
      do ii=1,nscgrain
       t1=t1+psi_squared_raw(n)
       n=n+1
      enddo
      psi_squared_smooth(i)=t1
      nn=nn+1
     enddo
     psi_squared_smooth(:)=psi_squared_smooth(:)*one_over_nscgrain

     filename='density_smoothed_t='//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
     do i=1,noscgcell
      write(14,'(2f25.16)')latpos1_smooth(i),density_smooth(i)
     enddo
     close(14)
     write(o,*)'Written smoothed coarse-grained density.'

     filename='psisq_smoothed_t='//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
     do i=1,noscgcell
      write(14,'(2f25.16)')latpos1_smooth(i),psi_squared_smooth(i)
     enddo
     close(14)
     write(o,*)'Written smoothed coarse-grained psi squared.'

     deallocate(density_smooth,psi_squared_smooth)

     call timer('SMOOTHING',.false.)

    endif ! plot_smooth

    if(plot_h_integrand)then
     if(num_ncgrain==1)then
      write(o,*)'Written integrand of coarse-grained H function formula.'
     else
      write(o,*)'Written integrands of coarse-grained H function formula.'
     endif
    endif
    if(hfunction)then
     if(num_ncgrain==1)then
      tmpr=r2s(H(1,nplot),'(f12.4)')
      write(o,*)'H function at this timestep: ',trim(tmpr)
     else
      write(o,*)'H function at this timestep: '
      do i=1,num_ncgrain
       tmpr=r2s(H(i,nplot),'(f12.4)')
       t1=(real(ncgrain(i),dp)*cell_x)/real(nlattice,dp)
       tmpr2=r2s(t1,'(f12.6)')
       if(ncgrain(i)<10)then
        write(o,*)trim(tmpr)//' (CG length '//trim(i2s(ncgrain(i)))//'  = ',&
         &trim(tmpr2)//')'
       else
        write(o,*)trim(tmpr)//' (CG length '//trim(i2s(ncgrain(i)))//' = ',&
         &trim(tmpr2)//')'
       endif
      enddo
     endif
    endif

    deallocate(density_raw,psi_squared_raw)
    if(save_backtracked.or.read_backtracked)deallocate(xback1)

   endif ! am_master

   den_time=den_time+den_timestep

   if(nfail_total(nplot)>0)latpos1=abs(latpos1)

  enddo ! nplot

  if(writefail)then
   if(nnodes==1.and.sum(nfail_total)==0)then
    close(21,status='delete')
   else
    close(21)
   endif
  endif

  return

1 call errstop('DENSITY_1D','Error reading file of backtracked positions.')
2 call errstop('DENSITY_1D','Premature end-of-file reading backtracked&
   & positions.')

  END SUBROUTINE density_1d


  SUBROUTINE density_2d(velocity,stepper)
!-------------------------------------------------------------------------!
! Calculate the raw, coarse-grained, and smoothed 2D densities and psi^2. !
!-------------------------------------------------------------------------!
  IMPLICIT NONE

  INTERFACE
   SUBROUTINE velocity(t,x,vel)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
   SUBROUTINE stepper(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
    REAL(dp),INTENT(inout) :: x(dim)
    REAL(dp),INTENT(in) :: htry,eps
    REAL(dp),INTENT(out) :: hdid,hnext
    REAL(dp),INTENT(inout) :: t
    INTERFACE
     SUBROUTINE velocity(t,x,vel)
      USE dsp
      USE store, ONLY : dim
      IMPLICIT NONE
      REAL(dp),INTENT(in) :: t
      REAL(dp),INTENT(in) :: x(dim)
      REAL(dp),INTENT(out) :: vel(dim)
     END SUBROUTINE velocity
    END INTERFACE
   END SUBROUTINE stepper
  END INTERFACE

  call timer('LATTICE GENERATION',.true.)

! Generate a raw 2d lattice.
  i=0
  do ik=1,nlattice
   do ij=1,nlattice
    i=i+1
    latpos2(1,i)=real(ij,dp)-0.5d0
    latpos2(2,i)=real(ik,dp)-0.5d0
   enddo
  enddo
  latpos2(:,:)=latpos2(:,:)*cell_over_nlat

  if(plot_smooth)then
! Generate a smoothing 2d lattice.
   jj=nscgrain/2 ; k=jj ; i=1
   do ik=1,noscgcell
    j=jj
    do ij=1,noscgcell
     latpos2_smooth(1,i)=(real(j,dp)-0.5d0)*cell_over_nlat
     latpos2_smooth(2,i)=(real(k,dp)-0.5d0)*cell_over_nlat
     j=j+nsmoothstep ; i=i+1
    enddo
    k=k+nsmoothstep
   enddo
  endif

  call timer('LATTICE GENERATION',.false.)

  if(read_backtracked)writefail=.false.

  if(writefail)then
   if(nnodes>1)then
    open(file='failed_trajectories'//trim(i2s(my_node))//'.dat',unit=21,&
     &status='unknown')
    write(21,*)'# NODE ',trim(i2s(my_node))
   else
    open(file='failed_trajectories.dat',unit=21,status='unknown')
   endif
   write(21,*)'# START TIME, START POSITION OF FAILED TRAJECTORIES'
  endif

  nfail_total(:)=0

! Loop over time sequence.
  do nplot=1,den_ntimes

   if(am_master)then
    if(den_time==0.d0)then
     tmpr='0.0'
     tmpr2='0.0'
    else
     tmpr=r2s(den_time,'(f25.16)')
     tmpr2=r2s(den_time,'(f12.2)')
    endif
   endif

   if(.not.read_backtracked)then

    call timer('BACKTRACK',.true.)

    if(am_master)then
     write(o,*)
     write(o,*)'Computing data at t = ',trim(adjustl(tmpr)),'.'
     if(den_time/=0.d0)then
      write(o,*)'Backtracking trajectories to zero for each point.'
     else
      write(o,*)'No backtracking required at time zero.'
     endif
    endif

    n=int((den_time)/(fourpi+1.d-6))+1
    maxstp_in=n*maxstp_in
    max_maxstp=n*max_maxstp

    nfail=0 ; nzero=0 ; nzero_total=0 ; nminh=0 ; nminh_total=0 ; n=my_node+1

    if(den_time==0.d0)then ! no propagation required

     do i=1,is_per_node(my_node) ! i.e. i=1,nlatticesq on 1 node
      x0(1)=latpos2(1,n)
      x0(2)=latpos2(2,n)
      call eval_psi_squared_t0_2d(x0,rho)
      psi_squared_raw_n(i)=rho
      call eval_density_t0_2d(dentype,x0,rho)
      density_raw_n(i)=rho
      n=n+nnodes
     enddo

    else ! t > 0.0

     do i=1,is_per_node(my_node) ! i.e. i=1,nlatticesq on 1 node
      x0(1)=latpos2(1,n)
      x0(2)=latpos2(2,n)
      if(phase_noise>0)call add_noise_to_phases
! psi^2 at lattice x, t
      call eval_psi_squared_2d(den_time,x0,psi2_t)
      if(psi2_t==0.d0)then
       nzero=nzero+1 ; good_traj=.false.
      else
       psi_squared_raw_n(i)=psi2_t
! Find backtracked positions.
       eps=init_eps ; xlast=1.d10 ; xdiff=1.d10 
       maxdiff=converge_maxdiff ; maxstp=maxstp_in
       do while(eps>mineps.and.xdiff>maxdiff)
        x1=x0
        call propagate(x1,den_time,0.d0,velocity,stepper)
        if(good_traj)then
         xdiff=sqrt((x1(1)-xlast(1))**2+(x1(2)-xlast(2))**2)
         xlast=x1
        else
         if(fastmode)then
          exit
         else
          maxstp=min(maxstp*10,max_maxstp)
          maxdiff=emergency_maxdiff
         endif
        endif
        eps=0.1d0*eps
       enddo
       noktot=noktot+nok
       nbadtot=nbadtot+nbad
       if(nstep_histogram)then
        ii=int(real(nok+nbad,dp)*nhist_scale)
        nhist(ii+1)=nhist(ii+1)+1
       endif
      endif
      if(good_traj)then
       if(save_backtracked)then
        xback1_n(i)=x1(1) ; xback2_n(i)=x1(2)
       endif
! rho at backtracked x, t0
       call eval_density_t0_2d(dentype,x1,rho_0)
! psi^2 at backtracked x, t0
       call eval_psi_squared_t0_2d(x1,psi2_0)
! rho at lattice x, t
       density_raw_n(i)=rho_0*psi2_t/psi2_0
      else
       nfail=nfail+1
       if(nminhfail)nminh=nminh+1
       if(writefail)write(21,*)den_time,latpos2(1:2,n),' 0.d0'
       latpos2(1,n)=-latpos2(1,n)
       psi_squared_raw_n(i)=0.d0
       density_raw_n(i)=0.d0
      endif
      n=n+nnodes
     enddo

    endif ! time > 0.d0 or not

    call timer('BACKTRACK',.false.)

    if(nnodes>1)then
     call timer('BARRIER',.true.)
     call barrier
     call timer('BARRIER',.false.)
    endif

    if(am_master)then
     allocate(density_raw(nlatticesq),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_2D','Insufficient memory to &
      &allocate raw density vector.')
     allocate(psi_squared_raw(nlatticesq),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_2D','Insufficient memory to &
      &allocate raw psi_squared vector.')
     if(save_backtracked)then
      allocate(xback1(nlatticesq),stat=ialloc)
      if(ialloc/=0)call errstop_master('DENSITY_2D','Insufficient memory to &
       &allocate XBACK1 vector.')
      allocate(xback2(nlatticesq),stat=ialloc)
      if(ialloc/=0)call errstop_master('DENSITY_2D','Insufficient memory to &
       &allocate XBACK2 vector.')
     endif   
    endif

    if(nnodes>1)then

     call timer('MPI',.true.)
     n=is_per_node(0)
     call mpi_type_vector(n,1,nnodes,mpi_double_precision,mdttype1,ierr)
     call mpi_type_commit(mdttype1,ierr)

! mpi_type_create_resized followed by mpi_gather or mpi_gatherv doesn't appear
! to work as advertised for interlaced receive vectors..(?) Use hideous set of 
! send/receives instead.
!  call mpi_type_create_resized(mdttype1,lb,extent,mdttype2,ierr)
!  call mpi_type_commit(mdttype2,ierr)
!  call mpi_gather(density_raw_n,n,mpi_double_precision,density_raw,1,mdttype2,&
!   &0,mpi_comm_world,ierr)
!  call mpi_gather(psi_squared_raw_n,n,mpi_double_precision,psi_squared_raw,&
!   &1,mdttype2,0,mpi_comm_world,ierr)
!  call mpi_type_free(mdttype2,ierr)

     if(nextra>0)then
      call mpi_type_vector(n-1,1,nnodes,mpi_double_precision,mdttype2,ierr)
      call mpi_type_commit(mdttype2,ierr)
     endif
     if(am_slave)then
      if(nextra==0.or.my_node<nextra)then
       call mpi_ssend(density_raw_n,n,mpi_double_precision,0,my_node,&
        &mpi_comm_world,ierr)
      else
       call mpi_ssend(density_raw_n,n-1,mpi_double_precision,0,my_node,&
        &mpi_comm_world,ierr)
      endif
      call checkmpi(ierr,'MPI_SSEND error on slave <1>')
     else ! am_master
      j=1
      do i=1,n
       density_raw(j)=density_raw_n(i) ; j=j+nnodes
      enddo
      if(nextra==0)then
       do i=1,nnodes-1
        call mpi_recv(density_raw(i+1),1,mdttype1,i,i,mpi_comm_world,status,&
         &ierr)
       enddo
      else
       do i=1,nextra-1
        call mpi_recv(density_raw(i+1),1,mdttype1,i,i,mpi_comm_world,status,&
         &ierr)
       enddo
       do i=nextra,nnodes-1
        call mpi_recv(density_raw(i+1),1,mdttype2,i,i,mpi_comm_world,status,&
         &ierr)
       enddo
      endif
      call checkmpi(ierr,'MPI_RECV error on master <1>')
     endif
     if(am_slave)then
      if(nextra==0.or.my_node<nextra)then
       call mpi_ssend(psi_squared_raw_n,n,mpi_double_precision,0,my_node,&
        &mpi_comm_world,ierr)
      else
       call mpi_ssend(psi_squared_raw_n,n-1,mpi_double_precision,0,my_node,&
        &mpi_comm_world,ierr)
      endif
      call checkmpi(ierr,'MPI_SSEND error on slave <2>')
     else ! am_master
      j=1
      do i=1,n
       psi_squared_raw(j)=psi_squared_raw_n(i) ; j=j+nnodes
      enddo
      if(nextra==0)then
       do i=1,nnodes-1
        call mpi_recv(psi_squared_raw(i+1),1,mdttype1,i,i,mpi_comm_world,&
         &status,ierr)
       enddo
      else
       do i=1,nextra-1
        call mpi_recv(psi_squared_raw(i+1),1,mdttype1,i,i,mpi_comm_world,&
         &status,ierr)
       enddo
       do i=nextra,nnodes-1
        call mpi_recv(psi_squared_raw(i+1),1,mdttype2,i,i,mpi_comm_world,&
         &status,ierr)
       enddo
      endif
      call checkmpi(ierr,'MPI_RECV error on master <2>')
     endif
     if(save_backtracked)then
      if(am_slave)then
       if(nextra==0.or.my_node<nextra)then
        call mpi_ssend(xback1_n,n,mpi_double_precision,0,my_node,&
         &mpi_comm_world,ierr)
       else
        call mpi_ssend(xback1_n,n-1,mpi_double_precision,0,my_node,&
         &mpi_comm_world,ierr)
       endif
       call checkmpi(ierr,'MPI_SSEND error on slave <3>')
      else ! am_master
       j=1
       do i=1,n
        xback1(j)=xback1_n(i) ; j=j+nnodes
       enddo
       if(nextra==0)then
        do i=1,nnodes-1
         call mpi_recv(xback1(i+1),1,mdttype1,i,i,mpi_comm_world,status,ierr)
        enddo
       else
        do i=1,nextra-1
         call mpi_recv(xback1(i+1),1,mdttype1,i,i,mpi_comm_world,status,ierr)
        enddo
        do i=nextra,nnodes-1
         call mpi_recv(xback1(i+1),1,mdttype2,i,i,mpi_comm_world,status,ierr)
        enddo
       endif
       call checkmpi(ierr,'MPI_RECV error on master <3>')
      endif ! master/slave
      if(am_slave)then
       if(nextra==0.or.my_node<nextra)then
        call mpi_ssend(xback2_n,n,mpi_double_precision,0,my_node,&
         &mpi_comm_world,ierr)
       else
        call mpi_ssend(xback2_n,n-1,mpi_double_precision,0,my_node,&
         &mpi_comm_world,ierr)
       endif
       call checkmpi(ierr,'MPI_SSEND error on slave <3>')
      else ! am_master
       j=1
       do i=1,n
        xback2(j)=xback2_n(i) ; j=j+nnodes
       enddo
       if(nextra==0)then
        do i=1,nnodes-1
         call mpi_recv(xback2(i+1),1,mdttype1,i,i,mpi_comm_world,status,ierr)
        enddo
       else
        do i=1,nextra-1
         call mpi_recv(xback2(i+1),1,mdttype1,i,i,mpi_comm_world,status,ierr)
        enddo
        do i=nextra,nnodes-1
         call mpi_recv(xback2(i+1),1,mdttype2,i,i,mpi_comm_world,status,ierr)
        enddo
       endif
       call checkmpi(ierr,'MPI_RECV error on master <3>')
      endif ! master/slave
     endif ! save_backtracked
     call mpi_type_free(mdttype1,ierr)
     if(nextra>0)call mpi_type_free(mdttype2,ierr)

    else

     density_raw(:)=density_raw_n(:)
     psi_squared_raw(:)=psi_squared_raw_n(:)
     if(save_backtracked)then
      xback1(:)=xback1_n(:) ; xback2(:)=xback2_n(:)
     endif

    endif

    call mpi_reduce(nfail,nfail_total(nplot),1,mpi_integer,mpi_sum,0,&
     &mpi_comm_world,ierr)
    call mpi_reduce(nzero,nzero_total,1,mpi_integer,mpi_sum,0,&
     &mpi_comm_world,ierr)
    call mpi_reduce(nminh,nminh_total,1,mpi_integer,mpi_sum,0,&
     &mpi_comm_world,ierr)

    if(nnodes>1)call timer('MPI',.false.)

    if(am_master)then
     write(o,*)'Done.'
     if(den_time/=0.d0)then
      if(nfail_total(nplot)>0)then
       write(o,'(1x,a)')'Total number of failed trajectories    : '&
        &//trim(i2s(nfail_total(nplot)))
       write(o,'(1x,a)')'Number failed for too many steps       : '&
        &//trim(i2s(nfail_total(nplot)-nzero_total-nminh_total))
       write(o,'(1x,a)')'Number failed for too small stepsize   : '&
        &//trim(i2s(nminh_total))
       write(o,'(1x,a)')'Number failed because started on nodes : '&
        &//trim(i2s(nzero_total))
      else
       write(o,'(1x,a)')'Number of failed trajectories: '//&
        &trim(i2s(nfail_total(nplot)))
      endif
     endif
     if(plot_raw.or.plot_cg.or.plot_smooth.or.plot_h_integrand)then
      write(o,*)'Writing requested data to disk.'
     endif
    endif

   else ! read_backtracked=T

    if(am_master)then

     allocate(density_raw(nlatticesq),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_2D','Cannot &
      &allocate raw density vector.')
     allocate(psi_squared_raw(nlatticesq),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_2D','Cannot &
      &allocate raw psi_squared vector.')
     allocate(xback1(nlatticesq),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_2D','Cannot &
      &allocate XBACK1 vector.')
     allocate(xback2(nlatticesq),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_2D','Cannot &
      &allocate XBACK2 vector.')

     write(o,*)
     write(o,*)'Computing data at t = ',trim(adjustl(tmpr)),'.'
     if(den_time/=0.d0)then
      write(o,*)'Reading previously backtracked positions from disk.'
     else
      write(o,*)'No backtracking required at time zero.'
     endif

     if(den_time==0.d0)then

      do i=1,nlatticesq
       x0(1)=latpos2(1,i)
       x0(2)=latpos2(2,i)
       call eval_psi_squared_t0_2d(x0,rho)
       psi_squared_raw(i)=rho
       call eval_density_t0_2d(dentype,x0,rho)
       density_raw(i)=rho
      enddo

     else ! den_time/=0.d0
 
      filename='backtracked_positions_t'//trim(adjustl(tmpr2))//'.dat'
      inquire(file=filename,exist=file_present)
      if(.not.file_present)then
       write(o,'(1x,a)')'Required filename: '//filename
       call errstop('DENSITY_2D','READ_BACKTRACKED=T but&
       & file containing backtracked positions is missing.')
      endif
      open(file=filename,unit=14,status='old')
      read(14,*,iostat=ierr)
      if(ierr/=0)call errstop('DENSITY_2D','Error reading backtracked&
       & positions file.')
      read(14,*,iostat=ierr)t1
      if(ierr/=0)call errstop('DENSITY_2D','Error reading backtracked&
       & positions file (second line).')
      if(abs(t1-den_time)>1.d-6)call errstop('DENSITY_2D','File containing&
       & backtracked positions has wrong timestamp.')
      do i=1,nlatticesq
       read(14,*,err=1,end=2)xback1(i),xback2(i)
      enddo
      write(o,*)'Done.'
      do i=1,nlatticesq
       if(xback1(i)>0.d0)then
        call eval_psi_squared_2d(den_time,latpos2(1:2,i),psi2_t)
        psi_squared_raw(i)=psi2_t
        x0(1)=xback1(i) ; x0(2)=xback2(i)
        call eval_density_t0_2d(dentype,x0,rho_0)
        call eval_psi_squared_t0_2d(x0,psi2_0)
        density_raw(i)=rho_0*psi2_t/psi2_0
       else
        nfail_total(nplot)=nfail_total(nplot)+1
        latpos2(1,i)=-latpos2(1,i)
        psi_squared_raw(i)=0.d0
        density_raw(i)=0.d0
       endif
      enddo

     endif ! den_time==0.d0 or not
 
    endif ! am_master

   endif ! read_backtracked or not

   if(am_master)then

! Write out backtracked positions if requested.
    if(save_backtracked.and.den_time/=0.d0)then

     call timer('RAW WRITE',.true.)

     filename='backtracked_positions_t'//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,*)'BACKTRACKED POSITIONS from time on second line.'
     write(14,*)den_time
     i=0
     do ik=1,nlattice
      do ij=1,nlattice
       i=i+1
       write(14,'(2f25.16)')xback1(i),xback2(i)
      enddo
     enddo
     close(14)
     write(o,*)'Written backtracked positions.'

     call timer('RAW WRITE',.false.)

    endif ! save_backtracked

    if(plot_raw)then
! Write out raw density and psi^2 if requested.

     call timer('RAW WRITE',.true.)
 
     filename='density_raw_t'//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
     write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
     i=0
     do ik=1,nlattice
      if(ik>1)write(14,*) ! gnuplot line separator
      do ij=1,nlattice
       i=i+1
       write(14,'(2f25.16,a,f25.16)')abs(latpos2(1,i)),latpos2(2,i),&
        &' 0.',density_raw(i)
      enddo
     enddo
     close(14)
     write(o,*)'Written raw density.'

     filename='psisq_raw_t'//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
     write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
     i=0
     do ik=1,nlattice
      if(ik>1)write(14,*) ! gnuplot line separator
      do ij=1,nlattice
       i=i+1
       write(14,'(2f25.16,a,f25.16)')abs(latpos2(1,i)),latpos2(2,i),&
        &' 0.',psi_squared_raw(i)
      enddo
     enddo
     close(14)
     write(o,*)'Written raw psi squared.'

     call timer('RAW WRITE',.false.)

    endif ! plot_raw

! Compute and write out non-overlapping coarse-grained densities if requested.
    if(plot_cg.or.hfunction)then
 
     call timer('COARSE-GRAINING',.true.)

     do ncg=1,num_ncgrain

      nc=ncgcell(ncg) ; ng=ncgrain(ncg)
      ncgcellsq=nc**2 ; ncgrainsq=ng**2
      one_over_ncgrainsq=1.d0/real(ncgrainsq,dp)
      ncgrow=nlattice*ng
      interpolation_required=.false.

! Generate a 2d coarse-graining lattice.
      i=0
      t1=real(ng,dp)*cell_over_nlat
      do ik=1,nc
       do ij=1,nc
        i=i+1
        latpos2_cg(1,i)=(real(ij,dp)-0.5d0)*t1
        latpos2_cg(2,i)=(real(ik,dp)-0.5d0)*t1
       enddo
      enddo

      allocate(density_cg(ncgcellsq),psi_squared_cg(ncgcellsq),stat=ialloc)
      if(ialloc/=0)call errstop_master('DENSITY_2D','Insufficient memory to &
       &allocate coarse-graining vectors.')

      k=0 ; nbase=1
      nstep=nlattice-ng
      if(nfail_total(nplot)==0)then
       do j=1,nc
        nn=0
        do i=1,nc
         t1=0.d0 ; t2=0.d0 ; n=nbase+nn*ng
         do jj=1,ng
          do ii=1,ng
           t1=t1+density_raw(n)
           t2=t2+psi_squared_raw(n)
           n=n+1
          enddo
          n=n+nstep
         enddo
         k=k+1 ; nn=nn+1
         density_cg(k)=t1 
         psi_squared_cg(k)=t2
        enddo
        nbase=ncgrow*j+1
       enddo
       density_cg(:)=density_cg(:)*one_over_ncgrainsq
       psi_squared_cg(:)=psi_squared_cg(:)*one_over_ncgrainsq
      else ! nfail_total > 0
       do j=1,nc
        nn=0
        do i=1,nc
         t1=0.d0 ; t2=0.d0 ; ngoodcg=ncgrainsq ; n=nbase+nn*ng
         do jj=1,ng
          do ii=1,ng
           if(latpos2(1,n)<0.d0)then
            ngoodcg=ngoodcg-1
           else
            t1=t1+density_raw(n)
            t2=t2+psi_squared_raw(n)
           endif
           n=n+1
          enddo
          n=n+nstep
         enddo
         k=k+1 ; nn=nn+1
         if(ngoodcg==0)then
          density_cg(k)=-1.d0 ; psi_squared_cg(k)=-1.d0
          interpolation_required=.true.
         else
          density_cg(k)=t1/real(ngoodcg,dp)
          psi_squared_cg(k)=t2/real(ngoodcg,dp)
         endif
        enddo
        nbase=ncgrow*j+1
       enddo
      endif ! nfail_total

! Interpolate to repair gaps in coarse-grained density and psi^2 due to 
! failed trajectories.

      if(interpolation_required)then
       k=1 ; kbase=1
       do j=1,nc
        do i=1,nc
         if(density_cg(k)==-1.d0)then
          ii=max(kbase,k-5) ; jj=min(kbase+nc-1,k+5)
          n=jj-ii ; nn=-1 ; nnn=0
          do kk=1,n
           if(ii+kk-1==k)nn=0
           if(density_cg(ii+kk+nn)/=-1.d0)then
            nnn=nnn+1
            points(nnn)=latpos2_cg(1,ii+kk+nn)
            repair(nnn)=density_cg(ii+kk+nn)
           endif
          enddo
          if(nnn>1)then
           call interp_nev(points(1:nnn),repair(1:nnn),nnn,latpos2_cg(1,k),&
            &val_interp,err_interp)
           density_cg(k)=val_interp
          endif
         endif
         if(density_cg(k)<0.d0)density_cg(k)=0.d0       
         k=k+1
        enddo
        kbase=kbase+nc
       enddo
       k=1 ; kbase=1
       do j=1,nc
        do i=1,nc
         if(psi_squared_cg(k)==-1.d0)then
          ii=max(kbase,k-5) ; jj=min(kbase+nc-1,k+5)
          n=jj-ii ; nn=-1 ; nnn=0
          do kk=1,n
           if(ii+kk-1==k)nn=0
           if(psi_squared_cg(ii+kk+nn)/=-1.d0)then
            nnn=nnn+1
            points(nnn)=latpos2_cg(1,ii+kk+nn)
            repair(nnn)=psi_squared_cg(ii+kk+nn)
           endif
          enddo
          if(nnn>1)then
           call interp_nev(points(1:nnn),repair(1:nnn),nnn,latpos2_cg(1,k),&
            &val_interp,err_interp)
           psi_squared_cg(k)=val_interp
          endif
         endif
         if(psi_squared_cg(k)<=0.d0)psi_squared_cg(k)=1.d-5
         k=k+1
        enddo
        kbase=kbase+nc
       enddo
      endif ! interpolation_required

      if(plot_cg)then

       tmpr=''
       if(num_ncgrain>1)tmpr='_cg'//trim(i2s(ng))
 
       filename='density_cg_t'//trim(adjustl(tmpr2))//trim(tmpr)//'.dat'
       open(file=filename,unit=14,status='unknown')
       write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
       write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
       n=1
       do j=1,nc
        if(j>1)write(14,*) ! gnuplot line separator
        do i=1,nc
         write(14,'(2f25.16,a,f25.16)')latpos2_cg(1,n),latpos2_cg(2,n),' 0.',&
          &density_cg(n)
         n=n+1
        enddo
       enddo
       close(14)
       if(num_ncgrain==1)write(o,*)'Written coarse-grained density.'
 
       filename='psisq_cg_t'//trim(adjustl(tmpr2))//trim(tmpr)//'.dat'
       open(file=filename,unit=14,status='unknown')
       write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
       write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
       n=1
       do j=1,nc
        if(j>1)write(14,*) ! gnuplot line separator
        do i=1,nc
         write(14,'(2f25.16,a,f25.16)')latpos2_cg(1,n),latpos2_cg(2,n),' 0.',&
          &psi_squared_cg(n)
         n=n+1
        enddo
       enddo
       close(14)
       if(num_ncgrain==1)write(o,*)'Written coarse-grained psi squared.'

      endif ! plot_cg

      if(hfunction.or.plot_h_integrand)then

       call timer('H FUNCTION',.true.)

       do i=1,ncgcellsq
        if(density_cg(i)/=0.d0.and.psi_squared_cg(i)/=0.d0)then
         density_cg(i)=density_cg(i)*log(density_cg(i)/psi_squared_cg(i))
        else
         density_cg(i)=0.d0 ! limit of x ln ax = 0 as x-->0+
        endif
       enddo

       if(plot_h_integrand)then
        filename='h_integrand_t'//trim(adjustl(tmpr2))//trim(tmpr)//'.dat'
        open(file=filename,unit=14,status='unknown')
        write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
        write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
        write(14,'(f25.16,a)')0.d0,' 0. 0. 0.'
        do i=1,nc
         write(14,'(f25.16,a)')latpos2_cg(1,i),' 0. 0. 0.'
        enddo
        write(14,'(f25.16,a)')cell_x,' 0. 0. 0.'
        n=1
        do j=1,nc
         t1=latpos2_cg(2,n)
         write(14,*) ! gnuplot line separator
         write(14,'(2f25.16,a,f25.16)')0.d0,t1,' 0.',0.d0
         do i=1,nc
          write(14,'(2f25.16,a,f25.16)')latpos2_cg(1,i),t1,' 0.',density_cg(n)
          n=n+1
         enddo
         write(14,'(2f25.16,a,f25.16)')cell_x,t1,' 0.',0.d0
        enddo
        write(14,*) ! gnuplot line separator
        write(14,'(2f25.16,a,f25.16)')0.d0,cell_y,' 0.',0.d0
        do i=1,nc
         write(14,'(2f25.16,a,f25.16)')latpos2_cg(1,i),cell_y,' 0.',0.d0
        enddo
        write(14,'(2f25.16,a,f25.16)')cell_x,cell_y,' 0.',0.d0
        close(14)
       endif

       if(hfunction)then
        allocate(htemp(nc),stat=ialloc)
        if(ialloc/=0)call errstop('DENSITY_2D','Unable to allocate HTEMP &
         &vector.')
        if(nc<6)call errstop('DENSITY_2D','Not enough coarse-graining&
         & cells to calculate the H function - need at least 6 per dimension.')

! Main integral
        t1=2.d0*latpos2_cg(1,1)
        n=0
        do i=1,nc
         htemp(i)=0.5d0*t1*(density_cg(n+1)+density_cg(n+nc)) ! linear to 0
         tt=density_cg(n+1)*0.375d0+density_cg(n+2)*seven_sixths+&
          &density_cg(n+3)*twentythree_over_twentyfour
         n=n+nc
         tt=tt+density_cg(n-2)*twentythree_over_twentyfour+&
          &density_cg(n-1)*seven_sixths+density_cg(n)*0.375d0
         htemp(i)=htemp(i)+(tt+sum(density_cg(n-nc+4:n-3)))*t1
        enddo
        htemp(1)=htemp(1)*0.375d0
        htemp(2)=htemp(2)*seven_sixths
        htemp(3)=htemp(3)*twentythree_over_twentyfour
        htemp(nc-2)=htemp(nc-2)*twentythree_over_twentyfour
        htemp(nc-1)=htemp(nc-1)*seven_sixths
        htemp(nc)=htemp(nc)*0.375d0
        H(ncg,nplot)=sum(htemp(:))*latpos2_cg(2,1)*2.d0
  
! Two missing rows
        t2=0.5d0*latpos2_cg(2,1)
        n=ncgcellsq-nc
        do i=1,nc
         htemp(i)=t2*(density_cg(i)+density_cg(n+i))
        enddo
        htemp(1)=htemp(1)*0.375d0
        htemp(2)=htemp(2)*seven_sixths
        htemp(3)=htemp(3)*twentythree_over_twentyfour
        htemp(nc-2)=htemp(nc-2)*twentythree_over_twentyfour
        htemp(nc-1)=htemp(nc-1)*seven_sixths
        htemp(nc)=htemp(nc)*0.375d0
        H(ncg,nplot)=H(ncg,nplot)+sum(htemp(:))*t1

! Four missing corners
        t3=t1*t2/6.d0
        H(ncg,nplot)=H(ncg,nplot)+t3*(density_cg(1)+density_cg(nc)+ &
         &density_cg(ncgcellsq-nc+1)+density_cg(ncgcellsq))

        deallocate(htemp)
       endif

       call timer('H FUNCTION',.false.)

      endif ! hfunction or plot_h_integrand

      deallocate(density_cg,psi_squared_cg)

     enddo ! num_ncgrain
     if(num_ncgrain>1)then
      write(o,*)'Written coarse-grained densities.'
      write(o,*)'Written coarse-grained psi-squareds.'
     endif

     call timer('COARSE-GRAINING',.false.)

    endif ! plot_cg or hfunction

! Compute and write out smoothed densities, if plot_smooth.
    if(plot_smooth)then
     call timer('SMOOTHING',.true.)

     allocate(density_smooth(noscgcellsq),psi_squared_smooth(noscgcellsq),&
      &stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_2D','Insufficient memory to &
      &allocate smoothing vectors.')

     nstep=nlattice-nscgrain
     k=0 ; nbase=1
     if(nfail_total(nplot)==0)then
      do j=1,noscgcell
       nn=0
       do i=1,noscgcell
        t1=0.d0 ; n=nbase+nn*nsmoothstep
        do jj=1,nscgrain
         do ii=1,nscgrain
          t1=t1+density_raw(n)
          n=n+1
         enddo
         n=n+nstep
        enddo
        k=k+1 ; nn=nn+1
        density_smooth(k)=t1
       enddo
       nbase=nscgrow*j+1
      enddo
      density_smooth(:)=density_smooth(:)*one_over_nscgrainsq
     else ! nfail_total > 0
      do j=1,noscgcell
       nn=0
       do i=1,noscgcell
        ngoodcg=nscgrainsq ; t1=0.d0 ; n=nbase+nn*nsmoothstep
        do jj=1,nscgrain
         do ii=1,nscgrain
          if(latpos2(1,n)<0.d0)then
           ngoodcg=ngoodcg-1
          else
           t1=t1+density_raw(n)
          endif
          n=n+1
         enddo
         n=n+nstep
        enddo
        k=k+1 ; nn=nn+1
        if(ngoodcg==0)then
         density_smooth(k)=0.d0
        else
         density_smooth(k)=t1/real(ngoodcg,dp)
        endif
       enddo
       nbase=nscgrow*j+1
      enddo
     endif ! nfail_total

     filename='density_smoothed_t'//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
     write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
     n=1
     do j=1,noscgcell
      if(j>1)write(14,*) ! gnuplot line separator
      do i=1,noscgcell
       write(14,'(2f25.16,a,f25.16)')latpos2_smooth(1,n),latpos2_smooth(2,n),&
        &' 0.',density_smooth(n)
       n=n+1
      enddo
     enddo
     close(14)
     write(o,*)'Written smoothed coarse-grained density.'

     k=0 ; nbase=1
     do j=1,noscgcell
      nn=0
      do i=1,noscgcell
       t1=0.d0 ; n=nbase+nn*nsmoothstep
       do jj=1,nscgrain
        do ii=1,nscgrain
         t1=t1+psi_squared_raw(n)
         n=n+1
        enddo
        n=n+nstep
       enddo
       k=k+1 ; nn=nn+1
       psi_squared_smooth(k)=t1
      enddo
      nbase=nscgrow*j+1
     enddo
     psi_squared_smooth(:)=psi_squared_smooth(:)*one_over_nscgrainsq

     filename='psisq_smoothed_t'//trim(adjustl(tmpr2))//'.dat'
     open(file=filename,unit=14,status='unknown')
     write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
     write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
     n=1
     do j=1,noscgcell
      if(j>1)write(14,*) ! gnuplot line separator
      do i=1,noscgcell
       write(14,'(2f25.16,a,f25.16)')latpos2_smooth(1,n),latpos2_smooth(2,n),&
        &' 0.',psi_squared_smooth(n)
       n=n+1
      enddo
     enddo

     close(14)
     write(o,*)'Written smoothed coarse-grained psi squared.'

     deallocate(density_smooth,psi_squared_smooth)

     call timer('SMOOTHING',.false.)

    endif ! plot_smooth

    if(plot_h_integrand)write(o,*)&
     &'Written integrand of coarse-grained H function formula.'
    if(hfunction)then
     if(num_ncgrain==1)then
      tmpr=r2s(H(1,nplot),'(f12.4)')
      write(o,*)'H function at this timestep: ',trim(tmpr)
     else
      write(o,*)'H function at this timestep: '
      do i=1,num_ncgrain
       tmpr=r2s(H(i,nplot),'(f12.4)')
       t1=(real(ncgrain(i),dp)*cell_x)/real(nlattice,dp)
       tmpr2=r2s(t1,'(f12.6)')
       if(ncgrain(i)<10)then
        write(o,*)trim(tmpr)//' (CG length '//trim(i2s(ncgrain(i)))//'  = '//&
         &trim(tmpr2)//')'
       else
        write(o,*)trim(tmpr)//' (CG length '//trim(i2s(ncgrain(i)))//' = '//&
         &trim(tmpr2)//')'
       endif
      enddo
     endif
    endif

    deallocate(density_raw,psi_squared_raw)
    if(save_backtracked.or.read_backtracked)deallocate(xback1,xback2)

   endif ! am_master

   den_time=den_time+den_timestep

   if(nfail_total(nplot)>0.and..not.read_backtracked)latpos2(1,:)=&
    &abs(latpos2(1,:))

  enddo ! nplot

  if(writefail)then
   if(nnodes==1.and.sum(nfail_total)==0)then
    close(21,status='delete')
   else
    close(21)
   endif
  endif

  return

1 call errstop('DENSITY_2D','Error reading file of backtracked positions.')
2 call errstop('DENSITY_2D','Premature end-of-file reading backtracked&
   & positions.')

  END SUBROUTINE density_2d


  SUBROUTINE density_3d(velocity,stepper)
!-------------------------------------------------------------------------!
! Calculate the coarse-grained 3D densities and psi^2 and the H function. !
!-------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER plane_start,plane_stop,nplane,signal,outcount,string,stringpoint,&
   &nstrings,stringbase,rgrequestn,rgrequest2n,rgrequest3n,inode
  REAL(dp) zcoord
  REAL(dp),ALLOCATABLE :: h_integrand(:)
  LOGICAL finished

  INTERFACE
   SUBROUTINE velocity(t,x,vel)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: t
    REAL(dp),INTENT(in) :: x(dim)
    REAL(dp),INTENT(out) :: vel(dim)
   END SUBROUTINE velocity
   SUBROUTINE stepper(x,vel,t,htry,eps,xscal,hdid,hnext,velocity)
    USE dsp
    USE store, ONLY : dim
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: vel(dim),xscal(dim)
    REAL(dp),INTENT(inout) :: x(dim)
    REAL(dp),INTENT(in) :: htry,eps
    REAL(dp),INTENT(out) :: hdid,hnext
    REAL(dp),INTENT(inout) :: t
    INTERFACE
     SUBROUTINE velocity(t,x,vel)
      USE dsp
      USE store, ONLY : dim
      IMPLICIT NONE
      REAL(dp),INTENT(in) :: t
      REAL(dp),INTENT(in) :: x(dim)
      REAL(dp),INTENT(out) :: vel(dim)
     END SUBROUTINE velocity
    END INTERFACE
   END SUBROUTINE stepper
  END INTERFACE

  call timer('LATTICE GENERATION',.true.)

  nstrings=nlattice*ncgrain(1)
  nlattice3=nlatticesq*ncgrain(1)
  ncgcell3=ncgcell(1)**3

  if(rgmaster)then
   allocate(density_raw(nlattice3),psi_squared_raw(nlattice3),stat=ialloc)
   if(ialloc/=0)call errstop_master('DENSITY_3D','Allocation problem: &
    &density_raw etc.')
   stopped(:)=.false.
  endif

! Which planes are assigned to my group?
  plane_start=1+(rggroup-1)*ncgrain(1)
  plane_stop=ncgrain(1)+(rggroup-1)*ncgrain(1)

! Generate a raw 3d lattice for the slab of planes assigned to this group.
  n=0 ; ik=0
group: do i=1,ncgcell(1)
   if(i==rggroup)then
    do j=1,ncgrain(1)
     ik=ik+1
     do ij=1,nlattice
      do ii=1,nlattice
       n=n+1
       latpos2(1,n)=(real(ii,dp)-0.5d0)
       latpos2(2,n)=(real(ij,dp)-0.5d0)
       latpos2(3,n)=(real(ik,dp)-0.5d0)
      enddo
     enddo
    enddo
    exit group
   else
    ik=ik+ncgrain(1)
   endif
  enddo group
  latpos2(:,:)=latpos2(:,:)*cell_over_nlat


  call timer('LATTICE GENERATION',.false.)

! Initialize failed trajectory stuff.
  if(writefail)then
   if(nnodes>1)then
    open(file='failed_trajectories'//trim(i2s(my_node))//'.dat',unit=21,&
     &status='unknown')
    write(21,*)'# NODE ',trim(i2s(my_node))
   else
    open(file='failed_trajectories.dat',unit=21,status='unknown')
   endif
   write(21,*)'# START TIME, START POSITION OF FAILED TRAJECTORIES'
  endif

  nfail_total(:)=0

! Loop over time sequence.
  do nplot=1,den_ntimes

   if(rgmaster)then
    if(den_time==0.d0)then
     tmpr='0.0'
     tmpr2='0.0'
    else
     tmpr=r2s(den_time,'(f25.16)')
     tmpr2=r2s(den_time,'(f12.2)')
    endif
   endif

   if(am_master)then
    write(o,*)
    write(o,*)'Computing data at t = ',trim(adjustl(tmpr)),'.'
    if(den_time/=0.d0)then
     write(o,*)'Backtracking trajectories to zero for each point.'
    else
     write(o,*)'No backtracking required at time zero.'
    endif
   endif

   if(rgmaster)then
    if(den_time==0.d0)then
     tmpr='0.0'
     tmpr2='0.0'
    else
     tmpr=r2s(den_time,'(f25.16)')
     tmpr2=r2s(den_time,'(f12.2)')
    endif
   endif

   n=int((den_time)/(fourpi+1.d-6))+1
   maxstp_in=n*maxstp_in
   max_maxstp=n*max_maxstp

   nfail=0 ; nzero=0 ; nzero_total=0 ; nminh=0 ; nminh_total=0

   if(den_time==0.d0)then ! no propagation required

!    n=1
!    if(.not.rgmaster)n=sum(is_per_node(0:my_rgnode-1))+1
!    do i=1,is_per_node(my_rgnode)
!     x0(1:3)=latpos2(1:3,n)
!     call eval_psi_squared_t0_3d(x0,rho)
!     psi_squared_raw_n(i)=rho
!     call eval_density_t0_3d(dentype,x0,rho)
!     density_raw_n(i)=rho
!     n=n+1
!    enddo
! Is it worth parallelizing this? If so need to gather 
! density_raw_n-->density_raw etc. on rgmaster later

    if(rgmaster)then
     do n=1,nlattice3
      x0(1:3)=latpos2(1:3,n)
      call eval_psi_squared_t0_3d(x0,rho)
      psi_squared_raw(n)=rho
      call eval_density_t0_3d(dentype,x0,rho)
      density_raw(n)=rho
     enddo
    endif

   else ! t > 0.0

    call timer('BACKTRACK',.true.)

    if(rgmaster)then
 
     finished=.false.

    call timer('MPI1',.true.)
     do inode=1,rgnnodes-1 ! all rgslaves
      call mpi_irecv(signal,1,mpi_integer,inode,1,rg_comm(rggroup),&
       &rgrequest(inode),ierror)
      call checkmpi(ierror,'mpi_irecv <1> in density_3d')
      n=(inode-1)*nlattice+1
      call mpi_irecv(density_raw(n),nlattice,mpi_double_precision,inode,2,&
       &rg_comm(rggroup),rgrequest2(inode),ierror)
      call checkmpi(ierror,'mpi_irecv <2> in density_3d')
      call mpi_irecv(psi_squared_raw(n),nlattice,mpi_double_precision,inode,3,&
       &rg_comm(rggroup),rgrequest3(inode),ierror)
      call checkmpi(ierror,'mpi_irecv <3> in density_3d')
     enddo ! rgslaves
    call timer('MPI1',.false.)

     stringbase=rgnnodes-1
     string=stringbase   

strings: do
     
! Wait for finished signals of any responding nodes
    call timer('MPI2',.true.)
      call mpi_waitsome(rgnnodes-1,rgrequest,outcount,which_rgnodes,&
       &mpi_statuses_ignore,ierror)
      call checkmpi(ierror,'mpi_waitsome in density_3d')
    call timer('MPI2',.false.)

! Flag last pass through string loop if now have enough assigned
! nodes to do all remaining strings.
      if(stringbase+outcount>=nstrings)finished=.true.

! Some or all of the now free nodes are required to do some more work
! and some might not be. Tell them which string to do in the former case,
! and tell them to quit in the latter.
    call timer('MPI3',.true.)
      do i=1,outcount
       string=string+1
       stringpoint=(string-1)*nlattice+1
       j=which_rgnodes(i)
       call mpi_isend(stringpoint,1,mpi_integer,j,4,rg_comm(rggroup),&
        &rgrequest(j),ierror)
       call checkmpi(ierror,'mpi_isend <1> in density_3d')
       if(string==nstrings)then
        do k=i+1,outcount
         stringpoint=0
         j=which_rgnodes(k)
         call mpi_isend(stringpoint,1,mpi_integer,j,4,rg_comm(rggroup),&
          &rgrequest(j),ierror)
         call checkmpi(ierror,'mpi_isend <2> in density_3d')
         stopped(j)=.true.
        enddo
        exit
       endif
      enddo
    call timer('MPI3',.false.)

! Wait for the density and psi_squared of all current responders.
    call timer('MPI4',.true.)
      do i=1,outcount
       j=which_rgnodes(i)
       call mpi_wait(rgrequest2(j),mpi_status_ignore,ierror)
       call checkmpi(ierror,'mpi_wait <1> in density_3d')
       call mpi_wait(rgrequest3(j),mpi_status_ignore,ierror)
       call checkmpi(ierror,'mpi_wait <2> in density_3d')
      enddo
    call timer('MPI4',.false.)

! Wait for instruction signals to be received.

    call timer('MPI5',.true.)
      do i=1,outcount
       j=which_rgnodes(i)
       call mpi_wait(rgrequest(j),mpi_status_ignore,ierror)
       call checkmpi(ierror,'mpi_wait <3> in density_3d')
      enddo
    call timer('MPI5',.false.)

! Post receives for signal, density, psi_squared for any
! newly initiated strings.
    call timer('MPI6',.true.)
      string=stringbase
      do i=1,outcount
       string=string+1
       j=which_rgnodes(i)
       call mpi_irecv(signal,1,mpi_integer,j,1,rg_comm(rggroup),&
        &rgrequest(j),ierror)
       call checkmpi(ierror,'mpi_irecv <4> in density_3d')
       n=(string-1)*nlattice+1
       call mpi_irecv(density_raw(n),nlattice,mpi_double_precision,j,2,&
        &rg_comm(rggroup),rgrequest2(j),ierror)
       call checkmpi(ierror,'mpi_irecv <5> in density_3d')
       call mpi_irecv(psi_squared_raw(n),nlattice,mpi_double_precision,j,3,&
        &rg_comm(rggroup),rgrequest3(j),ierror)
       call checkmpi(ierror,'mpi_irecv <6> in density_3d')
       if(string==nstrings)exit
      enddo
    call timer('MPI6',.false.)

! In final pass through string loop, wait for receipt of final 
! density and psi_squareds, and kill any node that hasn't already
! been stopped.
      if(finished)then
    call timer('MPI7',.true.)
       call mpi_waitall(rgnnodes-1,rgrequest,mpi_statuses_ignore,ierror)
       call checkmpi(ierror,'mpi_wait <4> in density_3d')
    call timer('MPI7',.false.)
    call timer('MPI8',.true.)
       call mpi_waitall(rgnnodes-1,rgrequest2,mpi_statuses_ignore,ierror)
       call checkmpi(ierror,'mpi_wait <5> in density_3d')
    call timer('MPI8',.false.)
       call mpi_waitall(rgnnodes-1,rgrequest3,mpi_statuses_ignore,ierror)
       call checkmpi(ierror,'mpi_wait <6> in density_3d')
    call timer('MPI9',.true.)
       do j=1,rgnnodes-1
        if(.not.stopped(j))then
         call mpi_isend(0,1,mpi_integer,j,4,rg_comm(rggroup),&
          &rgrequest(j),ierror)
         call checkmpi(ierror,'mpi_isend <3> in density_3d')
        endif 
       enddo
    call timer('MPI9',.false.)
       exit strings
      endif

      stringbase=string

     enddo strings

    else ! rgslave

     allocate(density_raw_n(nlattice),psi_squared_raw_n(nlattice),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_3D','Allocation problem: &
      &density_raw_n etc.')

     rgrequest2n=mpi_request_null ; rgrequest3n=mpi_request_null

     n=(my_rgnode-1)*nlattice+1

strings2: do ! strings
 
      do i=1,nlattice
       x0(1:3)=latpos2(1:3,n)
       if(phase_noise>0)call add_noise_to_phases
! psi^2 at lattice x, t
       call eval_psi_squared_3d(den_time,x0,psi2_t)
       if(psi2_t==0.d0)then
        nzero=nzero+1 ; good_traj=.false.
       else
        psi_squared_raw_n(i)=psi2_t
! Find backtracked positions.
        eps=init_eps ; xlast=1.d10 ; xdiff=1.d10 
        maxdiff=converge_maxdiff ; maxstp=maxstp_in
        do while(eps>mineps.and.xdiff>maxdiff)
         x1=x0
         call propagate(x1,den_time,0.d0,velocity,stepper)
         if(good_traj)then
          xdiff=sqrt((x1(1)-xlast(1))**2+(x1(2)-xlast(2))**2+(x1(3)-&
           &xlast(3))**2)
          xlast=x1
         else
          if(fastmode)then
           exit
          else
           maxstp=min(maxstp*10,max_maxstp)
           maxdiff=emergency_maxdiff
          endif
         endif
         eps=0.1d0*eps
        enddo
        noktot=noktot+nok
        nbadtot=nbadtot+nbad
        if(nstep_histogram)then
         ii=int(real(nok+nbad,dp)*nhist_scale)
         nhist(ii+1)=nhist(ii+1)+1
        endif
       endif
       n=n+1
       good_traj_save(i)=good_traj
       xsave(:,i)=x1(:)
      enddo ! i=1,nlattice

      call mpi_wait(rgrequest2n,mpi_status_ignore,ierror)
      call checkmpi(ierror,'mpi_wait <3> in density_3d')
      call mpi_wait(rgrequest3n,mpi_status_ignore,ierror)
      call checkmpi(ierror,'mpi_wait <4> in density_3d')

      do i=1,nlattice
       if(good_traj_save(i))then
! rho at backtracked x, t0
        x1(:)=xsave(:,i)
        call eval_density_t0_3d(dentype,x1,rho_0)
! psi^2 at backtracked x, t0
        call eval_psi_squared_t0_3d(x1,psi2_0)
! rho at lattice x, t
        density_raw_n(i)=rho_0*psi_squared_raw_n(i)/psi2_0
       else
        nfail=nfail+1
        if(nminhfail)nminh=nminh+1
        if(writefail)write(21,*)den_time,latpos2(1:3,i)
        latpos2(1,i)=-latpos2(1,i)
        psi_squared_raw_n(i)=0.d0
        density_raw_n(i)=0.d0
       endif
      enddo

      call mpi_isend(1,1,mpi_integer,0,1,rg_comm(rggroup),rgrequestn,ierror)
      call checkmpi(ierror,'mpi_isend <3> in density_3d')
      call mpi_isend(density_raw_n,nlattice,mpi_double_precision,0,2,&
       &rg_comm(rggroup),rgrequest2n,ierror)
      call checkmpi(ierror,'mpi_isend <4> in density_3d')
      call mpi_isend(psi_squared_raw_n,nlattice,mpi_double_precision,0,&
       &3,rg_comm(rggroup),rgrequest3n,ierror)
      call checkmpi(ierror,'mpi_isend <5> in density_3d')
      call mpi_recv(stringpoint,1,mpi_integer,0,4,rg_comm(rggroup),status,&
       &ierror)
      call checkmpi(ierror,'mpi_recv <1> in density_3d')
      call mpi_wait(rgrequestn,mpi_status_ignore,ierror)
      call checkmpi(ierror,'mpi_wait <5> in density_3d')

      if(stringpoint==0)exit strings2
      n=stringpoint

     enddo strings2 

     deallocate(density_raw_n,psi_squared_raw_n)
   
    endif

    call timer('BACKTRACK',.false.)

   endif ! time > 0.d0 or not

   call barrier

   call timer('MPI',.true.)
   call mpi_reduce(nfail,nfail_total(nplot),1,mpi_integer,mpi_sum,0,&
    &mpi_comm_world,ierr)
   call mpi_reduce(nzero,nzero_total,1,mpi_integer,mpi_sum,0,&
    &mpi_comm_world,ierr)
   call mpi_reduce(nminh,nminh_total,1,mpi_integer,mpi_sum,0,&
    &mpi_comm_world,ierr)
   call timer('MPI',.false.)

   if(rgmaster)then

    if(am_master)then
     write(o,*)'Done.'
     if(den_time/=0.d0)then
      if(nfail_total(nplot)>0)then
       write(o,'(1x,a)')'Total number of failed trajectories    : '&
        &//trim(i2s(nfail_total(nplot)))
       write(o,'(1x,a)')'Number failed for too many steps       : '&
        &//trim(i2s(nfail_total(nplot)-nzero_total-nminh_total))
       write(o,'(1x,a)')'Number failed for too small stepsize   : '&
        &//trim(i2s(nminh_total))
       write(o,'(1x,a)')'Number failed because started on nodes : '&
        &//trim(i2s(nzero_total))
      else
       write(o,'(1x,a)')'Number of failed trajectories: '//&
        &trim(i2s(nfail_total(nplot)))
      endif
     endif
     if(plot_raw.or.plot_cg.or.plot_smooth.or.plot_h_integrand)then
      write(o,*)'Writing requested data to disk.'
     endif
    endif ! am_master

    if(plot_raw)then
! Write out raw density and psi^2 if requested for planes selected in input.

     call timer('RAW WRITE',.true.)

     do nplane=1,n3dplanes

      if(plane_number(nplane)<plane_start.or.plane_number(nplane)>&
       &plane_stop)cycle  ! if this plane is not in my comm group
  
      filename='density_raw_t'//trim(adjustl(tmpr2))//'_p'//&
       &trim(adjustl(i2s(plane_number(nplane))))//'.dat'
      open(file=filename,unit=14,status='unknown')
      write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
      write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
!      write(14,'(a,f25.16)')'#ZRANGE 0.0 ',cell_z
      i=(plane_number(nplane)-plane_start)*nlatticesq
      do ik=1,nlattice
       if(ik>1)write(14,*) ! gnuplot line separator
       do ij=1,nlattice
        i=i+1
        write(14,'(4f25.16)')abs(latpos2(1,i)),latpos2(2,i),&
         &latpos2(3,i),density_raw(i)
       enddo
      enddo
      close(14)
      if(am_master)write(o,*)'Written raw density.'

      filename='psisq_raw_t'//trim(adjustl(tmpr2))//'_p'//&
       &trim(adjustl(i2s(plane_number(nplane))))//'.dat'
      open(file=filename,unit=14,status='unknown')
      write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
      write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
!      write(14,'(a,f25.16)')'#ZRANGE 0.0 ',cell_z
      i=(plane_number(nplane)-plane_start)*nlatticesq
      do ik=1,nlattice
       if(ik>1)write(14,*) ! gnuplot line separator
       do ij=1,nlattice
        i=i+1
        write(14,'(4f25.16)')abs(latpos2(1,i)),latpos2(2,i),&
         &latpos2(3,i),psi_squared_raw(i)
       enddo
      enddo
      close(14)
      if(am_master)write(o,*)'Written raw psi squared.'

     enddo ! planes

     call timer('RAW WRITE',.false.)

    endif ! plot_raw

! Compute and write out non-overlapping coarse-grained densities if requested.
    if(plot_cg.or.hfunction)then

     call timer('COARSE_GRAINING',.true.)

     nc=ncgcell(1) ; ng=ncgrain(1) 
     ncgcellsq=nc*nc ; ncgrainsq=ng*ng
     ncgrain3=ncgrain(1)**3
     one_over_ncgrain3=1.d0/real(ncgrain3,dp)
     ncgrow=nlattice*ng
     interpolation_required=.false.

! Generate a coarse-graining lattice for the current 3D slab.
     i=0
     do ik=1,nc
      do ij=1,nc
       i=i+1
       latpos2_cg(1,i)=real(ij,dp)-0.5d0
       latpos2_cg(2,i)=real(ik,dp)-0.5d0
      enddo
     enddo
     latpos2_cg(3,:)=real(rggroup,dp)-0.5d0
     t1=real(ng,dp)*cell_over_nlat
     latpos2_cg(:,:)=latpos2_cg(:,:)*t1

     allocate(density_cg(ncgcellsq),psi_squared_cg(ncgcellsq),stat=ialloc)
     if(ialloc/=0)call errstop_master('DENSITY_3D','Insufficient memory to &
      &allocate coarse-graining vectors.')

     k=0 ; nbase=1
     nstep=nlattice-ng

     if(nfail_total(nplot)==0)then

      do j=1,nc
       nn=0
       do i=1,nc
        t1=0.d0 ; t2=0.d0 
        do kk=0,ng-1
         n=nbase+nn*ng+kk*nlatticesq
         do jj=1,ng
          do ii=1,ng
           t1=t1+density_raw(n)
           t2=t2+psi_squared_raw(n)
           n=n+1
          enddo ! ii
          n=n+nstep
         enddo ! jj
        enddo ! kk
        k=k+1 ; nn=nn+1
        density_cg(k)=t1 
        psi_squared_cg(k)=t2
       enddo ! i
       nbase=ncgrow*j+1
      enddo ! j
      density_cg(:)=density_cg(:)*one_over_ncgrain3
      psi_squared_cg(:)=psi_squared_cg(:)*one_over_ncgrain3

     else ! nfail_total > 0

      do j=1,nc
       nn=0
       do i=1,nc
        t1=0.d0 ; t2=0.d0 ; ngoodcg=ncgrainsq 
        do kk=0,ng-1
         n=nbase+nn*ng+kk*nlatticesq
         do jj=1,ng
          do ii=1,ng
           if(latpos2(1,n)<0.d0)then
            ngoodcg=ngoodcg-1
           else
            t1=t1+density_raw(n)
            t2=t2+psi_squared_raw(n)
           endif
           n=n+1
          enddo ! ii
          n=n+nstep
         enddo ! jj
        enddo ! kk
        k=k+1 ; nn=nn+1
        if(ngoodcg==0)then
         density_cg(k)=-1.d0 ; psi_squared_cg(k)=-1.d0
         interpolation_required=.true.
        else
         density_cg(k)=t1/real(ngoodcg,dp)
         psi_squared_cg(k)=t2/real(ngoodcg,dp)
        endif
       enddo
       nbase=ncgrow*j+1
      enddo

     endif ! nfail_total

     if(interpolation_required)call errstop('DENSITY_3D','Interpolation of &
      &3D coarse-grained density required, but Mike hasn''t coded it up yet.')

     if(plot_cg)then

      if(any(plane_number(:)>=plane_start.and.plane_number(:)<=plane_stop))then
      
       filename='density_cg_t'//trim(adjustl(tmpr2))//'_pc'//&
        &trim(adjustl(i2s(rggroup)))//'.dat'
       open(file=filename,unit=16,status='unknown')
       write(16,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
       write(16,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
!       write(16,'(a,f25.16)')'#ZRANGE 0.0 ',cell_z
       n=0
!      do k=1,ncgcell(1)
!       if(k>1)write(16,*) ! gnuplot plane separator
        do j=1,ncgcell(1)
         if(j>1)write(16,*) ! gnuplot line separator
         do i=1,ncgcell(1)
          n=n+1
          write(16,'(4f25.16)')latpos2_cg(1,n),latpos2_cg(2,n),latpos2_cg(3,n),&
           &density_cg(n)
         enddo
        enddo
!      enddo

       filename='psisq_cg_t'//trim(adjustl(tmpr2))//'_pc'//&
        &trim(adjustl(i2s(rggroup)))//'.dat'
       open(file=filename,unit=17,status='unknown')
       write(17,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
       write(17,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
!       write(17,'(a,f25.16)')'#ZRANGE 0.0 ',cell_z
       n=0
!      do k=1,ncgcell(1)
!       if(k>1)write(17,*) ! gnuplot plane separator
        do j=1,ncgcell(1)
         if(j>1)write(17,*) ! gnuplot line separator
         do i=1,ncgcell(1)
          n=n+1
          write(17,'(4f25.16)')latpos2_cg(1,n),latpos2_cg(2,n),latpos2_cg(3,n),&
           &psi_squared_cg(n)
         enddo
        enddo
!      enddo
       if(am_master)then
        write(o,*)'Written coarse-grained density.'
        write(o,*)'Written coarse-grained psi squared.'
       endif
 
       close(16) ; close(17)

      endif ! this cg plane is to be plotted.
     endif ! plot_cg

     if(hfunction.or.plot_h_integrand)then

      call timer('H FUNCTION',.true.)

      do i=1,ncgcellsq 
       if(density_cg(i)/=0.d0.and.psi_squared_cg(i)/=0.d0)then
        density_cg(i)=density_cg(i)*log(density_cg(i)/psi_squared_cg(i))
       else
        density_cg(i)=0.d0 ! limit of x ln ax = 0 as x-->0+
       endif
      enddo

      if(plot_h_integrand)then
       filename='h_integrand_t='//trim(adjustl(tmpr2))//'.dat'
       open(file=filename,unit=14,status='unknown')
       write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
       write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
       write(14,'(a,f25.16)')'#ZRANGE 0.0 ',cell_z

       write(14,'(f25.16,a)')0.d0,' 0. 0. 0.'
       do i=1,ncgcell(1)
        write(14,'(f25.16,a)')latpos2_cg(1,i),' 0. 0. 0.'
       enddo
       write(14,'(f25.16,a)')cell_x,' 0. 0. 0.'

       n=1
       do j=1,ncgcell(1)
        t1=latpos2_cg(2,n)
        write(14,*) ! gnuplot line separator
        write(14,'(2f25.16,a,f25.16)')0.d0,t1,' 0.',0.d0
        do i=1,ncgcell(1)
         write(14,'(2f25.16,a,f25.16)')latpos2_cg(1,i),t1,' 0.',density_cg(n)
         n=n+1
        enddo
        write(14,'(2f25.16,a,f25.16)')cell_x,t1,' 0.',0.d0
       enddo
       write(14,*) ! gnuplot line separator
       write(14,'(2f25.16,a,f25.16)')0.d0,cell_y,' 0.',0.d0
       do i=1,ncgcell(1)
        write(14,'(2f25.16,a,f25.16)')latpos2_cg(1,i),cell_y,' 0.',0.d0
       enddo
       write(14,'(2f25.16,a,f25.16)')cell_x,cell_y,' 0.',0.d0

       do k=1,ncgcell(1)
        zcoord=(real(k,dp)-0.5d0)*real(ncgrain(1),dp)*cell_over_nlat
        write(14,*) ! gnuplot line separator
        write(14,'(2(f25.16,a))')0.d0,' 0. ',zcoord,' 0.'
        do i=1,ncgcell(1)
         write(14,'(2(f25.16,a))')latpos2_cg(1,i),' 0. ',zcoord,' 0.'
        enddo
        write(14,'(2(f25.16,a))')cell_x,' 0. ',zcoord,' 0.'
        n=1
        do j=1,ncgcell(1)
         t1=latpos2_cg(2,n)
         write(14,*) ! gnuplot line separator
         write(14,'(4f25.16)')0.d0,t1,zcoord,0.d0
         do i=1,ncgcell(1)
          write(14,'(4f25.16)')latpos2_cg(1,i),t1,zcoord,density_cg(n)
          n=n+1
         enddo
         write(14,'(4f25.16)')cell_x,t1,zcoord,0.d0
        enddo
        write(14,*) ! gnuplot line separator
        write(14,'(4f25.16)')0.d0,cell_y,zcoord,0.d0
        do i=1,ncgcell(1)
         write(14,'(4f25.16)')latpos2_cg(1,i),cell_y,zcoord,0.d0
        enddo
        write(14,'(4f25.16)')cell_x,cell_y,zcoord,0.d0
       enddo ! k
 
       write(14,'(2(f25.16,a))')0.d0,' 0. ',cell_z,' 0.'
       do i=1,ncgcell(1)
        write(14,'(2(f25.16,a))')latpos2_cg(1,i),' 0. ',cell_z,' 0.'
       enddo
       write(14,'(2(f25.16,a))')cell_x,' 0. ',cell_z,' 0.'
       n=1
       do j=1,ncgcell(1)
        t1=latpos2_cg(2,n)
        write(14,*) ! gnuplot line separator
        write(14,'(4f25.16)')0.d0,t1,cell_z,0.d0
        do i=1,ncgcell(1)
         write(14,'(4f25.16)')latpos2_cg(1,i),t1,cell_z,density_cg(n)
         n=n+1
        enddo
        write(14,'(4f25.16)')cell_x,t1,cell_z,0.d0
       enddo
       write(14,*) ! gnuplot line separator
       write(14,'(4f25.16)')0.d0,cell_y,cell_z,0.d0
       do i=1,ncgcell(1)
        write(14,'(4f25.16)')latpos2_cg(1,i),cell_y,cell_z,0.d0
       enddo
       write(14,'(4f25.16)')cell_x,cell_y,cell_z,0.d0
       close(14)
 
       if(am_master)then
        write(o,*)'Written integrand of coarse-grained H function formula.'
       endif
      endif ! h_integrand

! H function calculation
      if(hfunction)then

! Gather all slabs H integrands on the master.
       call mpi_barrier(rg_comm(nredistgrps+1),ierror)
       if(am_master)allocate(h_integrand(ncgcell3),stat=ialloc)
       if(ialloc/=0)call errstop('DENSITY_3D','Unable to allocate H_INTEGRAND &
        &vector.')

       call mpi_gather(density_cg,ncgcellsq,mpi_double_precision,h_integrand,&
        &ncgcellsq,mpi_double_precision,0,rg_comm(nredistgrps+1),ierror)

       if(am_master)then

        allocate(htemp(ncgcell(1)),htemp2(ncgcell(1)),stat=ialloc)
        if(ialloc/=0)call errstop('DENSITY_3D','Unable to allocate &
         &HTEMP/HTEMP2 vectors.')

        if(ncgcell(1)<6)call errstop('DENSITY_3D','Not enough coarse-graining&
         & cells to calculate the H function - need at least 6 per dimension.')
   
! Main integral
        n=0
        t1=2.d0*latpos2_cg(1,1) 
        t2=2.d0*latpos2_cg(2,1) 
        t3=2.d0*latpos2_cg(3,1)

        do j=1,ncgcell(1) ! loop over xy planes at increasing z
         do i=1,ncgcell(1)
          htemp(i)=0.25d0*t1*(h_integrand(n+1)+h_integrand(n+ncgcell(1)))!lin->0
          tt=h_integrand(n+1)*0.375d0+h_integrand(n+2)*seven_sixths+&
           &h_integrand(n+3)*twentythree_over_twentyfour
          n=n+ncgcell(1)
          tt=tt+h_integrand(n-2)*twentythree_over_twentyfour+&
           &h_integrand(n-1)*seven_sixths+h_integrand(n)*0.375d0
          htemp(i)=htemp(i)+(tt+sum(h_integrand(n-ncgcell(1)+4:n-3)))*t1
         enddo ! i
         htemp(1)=htemp(1)*0.375d0
         htemp(2)=htemp(2)*seven_sixths
         htemp(3)=htemp(3)*twentythree_over_twentyfour
         htemp(ncgcell(1)-2)=htemp(ncgcell(1)-2)*twentythree_over_twentyfour
         htemp(ncgcell(1)-1)=htemp(ncgcell(1)-1)*seven_sixths
         htemp(ncgcell(1))=htemp(ncgcell(1))*0.375d0
         htemp2(j)=sum(htemp(:))*t2
        enddo ! j
        htemp2(1)=htemp2(1)*0.375d0
        htemp2(2)=htemp2(2)*seven_sixths
        htemp2(3)=htemp2(3)*twentythree_over_twentyfour
        htemp2(ncgcell(1)-2)=htemp2(ncgcell(1)-2)*twentythree_over_twentyfour
        htemp2(ncgcell(1)-1)=htemp2(ncgcell(1)-1)*seven_sixths
        htemp2(ncgcell(1))=htemp2(ncgcell(1))*0.375d0
        H(1,nplot)=sum(htemp2(:))*t3

! Two side slabs including one pair of edges
        nnn=ncgcell(1)**2-ncgcell(1)
        do j=1,ncgcell(1) ! z
         k=(j-1)*ncgcell(1)**2
         n=k+nnn
         htemp2(j)=0.25d0*t1*(h_integrand(k+1)+h_integrand(k+ncgcell(1))+&
          h_integrand(n+1)+h_integrand(n+ncgcell(1))) ! lin->0
         do i=1,ncgcell(1) ! x
          k=k+1 ; n=n+1
          htemp(i)=h_integrand(k)+h_integrand(n)
         enddo
         htemp(1)=htemp(1)*0.375d0
         htemp(2)=htemp(2)*seven_sixths
         htemp(3)=htemp(3)*twentythree_over_twentyfour
         htemp(ncgcell(1)-2)=htemp(ncgcell(1)-2)*twentythree_over_twentyfour
         htemp(ncgcell(1)-1)=htemp(ncgcell(1)-1)*seven_sixths
         htemp(ncgcell(1))=htemp(ncgcell(1))*0.375d0
         htemp2(j)=htemp2(j)+t1*sum(htemp(:))
         htemp2(j)=0.25d0*t2*htemp2(j) ! lin->0
        enddo
        htemp2(1)=htemp2(1)*0.375d0
        htemp2(2)=htemp2(2)*seven_sixths
        htemp2(3)=htemp2(3)*twentythree_over_twentyfour
        htemp2(ncgcell(1)-2)=htemp2(ncgcell(1)-2)*twentythree_over_twentyfour
        htemp2(ncgcell(1)-1)=htemp2(ncgcell(1)-1)*seven_sixths
        htemp2(ncgcell(1))=htemp2(ncgcell(1))*0.375d0
        H(1,nplot)=H(1,nplot)+sum(htemp2(:))*t3
 
! Bottom and top slab including one pair of edges
        n=0 ; k=ncgcell3-ncgcell(1)**2
        do i=1,ncgcell(1)
         htemp(i)=0.25d0*t1*(h_integrand(n+1)+h_integrand(n+ncgcell(1))+&
          &h_integrand(k+1)+h_integrand(k+ncgcell(1))) ! lin->0

         tt=h_integrand(n+1)*0.375d0+h_integrand(n+2)*seven_sixths+&
          &h_integrand(n+3)*twentythree_over_twentyfour
         n=n+ncgcell(1)
         tt=tt+h_integrand(n-2)*twentythree_over_twentyfour+&
          &h_integrand(n-1)*seven_sixths+h_integrand(n)*0.375d0
         tt=tt+sum(h_integrand(n-ncgcell(1)+4:n-3))

         tt=tt+h_integrand(k+1)*0.375d0+h_integrand(k+2)*seven_sixths+&
          &h_integrand(k+3)*twentythree_over_twentyfour
         k=k+ncgcell(1)
         tt=tt+h_integrand(k-2)*twentythree_over_twentyfour+&
          &h_integrand(k-1)*seven_sixths+h_integrand(k)*0.375d0
         tt=(tt+sum(h_integrand(k-ncgcell(1)+4:k-3)))*t1
         htemp(i)=0.25d0*t3*(htemp(i)+tt) ! lin->0
        enddo
        htemp(1)=htemp(1)*0.375d0
        htemp(2)=htemp(2)*seven_sixths
        htemp(3)=htemp(3)*twentythree_over_twentyfour
        htemp(ncgcell(1)-2)=htemp(ncgcell(1)-2)*twentythree_over_twentyfour
        htemp(ncgcell(1)-1)=htemp(ncgcell(1)-1)*seven_sixths
        htemp(ncgcell(1))=htemp(ncgcell(1))*0.375d0
        H(1,nplot)=H(1,nplot)+sum(htemp(:))*t2

! Four missing edges
        do i=1,4
         select case(i)
          case(1) ; n=0
          case(2) ; n=ncgcell(1)**2-ncgcell(1)
          case(3) ; n=ncgcell3-ncgcell(1)**2
          case(4) ; n=ncgcell3-ncgcell(1)
         end select
         htemp(1)=0.25d0*t1*(h_integrand(n+1)+h_integrand(n+ncgcell(1)))!lin->0
         tt=h_integrand(n+1)*0.375d0+h_integrand(n+2)*seven_sixths+&
          &h_integrand(n+3)*twentythree_over_twentyfour
         n=n+ncgcell(1)
         tt=tt+h_integrand(n-2)*twentythree_over_twentyfour+&
          &h_integrand(n-1)*seven_sixths+h_integrand(n)*0.375d0
         tt=(tt+sum(h_integrand(n-ncgcell(1)+4:n-3)))*t1
         htemp(1)=0.0625d0*t3*t2*(htemp(1)+tt) ! lin->0
         H(1,nplot)=H(1,nplot)+htemp(1)
        enddo

        deallocate(htemp,htemp2)

        tmpr=r2s(H(1,nplot),'(f12.4)')
        write(o,*)'H function at this timestep: ',trim(tmpr)

       endif ! am_master

      endif ! hfunction

      call timer('H FUNCTION',.false.)
      call timer('COARSE_GRAINING',.false.)

     endif ! hfunction or plot_h_integrand

     deallocate(density_cg,psi_squared_cg)
     if(am_master)deallocate(h_integrand)
  
    endif ! plot_cg or h function

   endif ! rgmaster

   if(nfail>0)latpos2(1,:)=abs(latpos2(1,:))

   den_time=den_time+den_timestep
 
   if(rgmaster)stopped(:)=.false.
   call barrier
  enddo ! nplot

  if(writefail)then
   if(nnodes==1.and.sum(nfail_total)==0)then
    close(21,status='delete')
   else
    close(21)
   endif
  endif

  END SUBROUTINE density_3d

 
  SUBROUTINE setup_comm_groups
!----------------------------------------------------------------------------!
! SETUP_COMM_GROUPS                                                          !
! -------------------                                                        !
! Setup the structure of the subgroups of processors ('comm groups') that    !
! are used in the parallel 3D algorithm.                                     !
!                                                                            !
! MDT 1.2011                                                                 !
!----------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER i,j,k,world_group,ialloc
 INTEGER,ALLOCATABLE :: rg_ranks(:),rg_grp(:)

 nredistgrps=ncgcell(1)
 allocate(rg_nnodes(nredistgrps),rgstart(nredistgrps),rgstop(nredistgrps),&
  &rg_grp(nredistgrps+1),rg_comm(nredistgrps+1),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_COMM_GROUPS','Allocation problem: &
  &rg_nnodes etc.')
 rg_nnodes(:)=nnodes/nredistgrps
 j=mod(nnodes,nredistgrps)
 do k=1,j
  rg_nnodes(k)=rg_nnodes(k)+1
 enddo
 rgstart(1)=0 ; rgstop(1)=rg_nnodes(1)-1 ; rggroup=1
 do i=2,nredistgrps
  rgstart(i)=rgstop(i-1)+1
  rgstop(i)=rgstart(i)+rg_nnodes(i)-1
  if(my_node>=rgstart(i))rggroup=i
 enddo
 rgnnodes=rg_nnodes(rggroup)
 if(any(my_node==rgstart(:)))then
  rgmaster=.true.
  rgmymaster=my_node
 else
  rgmaster=.false.
  rgmymaster=rgstart(rggroup)
 endif

 call mpi_comm_group(mpi_comm_world,world_group,ierror)
 call checkmpi(ierror,'mpi_comm_group in setup_comm_groups')

 do i=1,nredistgrps
  allocate(rg_ranks(rg_nnodes(i)),stat=ialloc)
  if(ialloc/=0)call errstop('SETUP_COMM_GROUPS','Allocation problem: &
   &rg_ranks')
  k=0
  do j=rgstart(i),rgstop(i)
   k=k+1
   rg_ranks(k)=j
  enddo
  call mpi_group_incl(world_group,rg_nnodes(i),rg_ranks,rg_grp(i),&
   &ierror)
  call checkmpi(ierror,'mpi_group_incl in setup_comm_groups')
  call mpi_comm_create(mpi_comm_world,rg_grp(i),rg_comm(i),ierror)
  call checkmpi(ierror,'mpi_comm_create in setup_comm_groups')
  deallocate(rg_ranks)
 enddo

 allocate(rg_ranks(nredistgrps),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_COMM_GROUPS','Allocation problem: &
  &rg_ranks')

 do i=1,nredistgrps
  rg_ranks(i)=rgstart(i)
 enddo

 call mpi_group_incl(world_group,nredistgrps,rg_ranks,rg_grp(nredistgrps+1),&
  &ierror)
 call checkmpi(ierror,'mpi_group_incl in setup_comm_groups <2>')
 call mpi_comm_create(mpi_comm_world,rg_grp(nredistgrps+1),&
  &rg_comm(nredistgrps+1),ierror)
 call checkmpi(ierror,'mpi_comm_create in setup_comm_groups <2>')
 deallocate(rg_ranks,rg_grp)

 call mpi_comm_rank(rg_comm(rggroup),my_rgnode,ierror)
 call checkmpi(ierror,'mpi_comm_rank in setup_comm_groups <2>')

 if(rgmaster)then
  call mpi_comm_rank(rg_comm(nredistgrps+1),my_rgnode2,ierror)
  call checkmpi(ierror,'mpi_comm_rank in setup_comm_groups <2>')
  allocate(stopped(rgnnodes-1),stat=ialloc)
  if(ialloc/=0)call errstop('SETUP_COMM_GROUPS','Allocation problem: &
   & stopped array')
 endif

 allocate(rgrequest(rgnnodes-1),rgrequest2(rgnnodes-1),rgrequest3(rgnnodes-1),&
  &which_rgnodes(rgnnodes-1),good_traj_save(nlattice),&
  &xsave(3,nlattice),stat=ialloc)
 if(ialloc/=0)call errstop('SETUP_COMM_GROUPS','Allocation problem: &
  &rgrequest')

 END SUBROUTINE setup_comm_groups


END PROGRAM louis

! Interpolate to repair gaps in coarse-grained density and psi^2 due to
! failed trajectories
!    if(interpolation_required)then
!     m=1 ; mbase=1
!     do k=1,ncgcell(1)
!      do j=1,ncgcell(1)
!       do i=1,ncgcell(1)
!        if(density_cg(m)==-1.d0)then
!         ii=max(mbase,m-5) ; jj=min(mbase+ncgcell(1)-1,m+5)
!         n=jj-ii ; nn=-1 ; nnn=0
!         do kk=1,n
!          if(ii+kk-1==m)nn=0
!          if(density_cg(ii+kk+nn)/=-1.d0)then
!           nnn=nnn+1
!           points(nnn)=latpos2_cg(1,ii+kk+nn)
!           repair(nnn)=density_cg(ii+kk+nn)
!          endif
!         enddo
!         if(nnn>1)then
!          call interp_nev(points(1:nnn),repair(1:nnn),nnn,latpos2_cg(1,m),&
!           &val_interp,err_interp)
!          density_cg(m)=val_interp
!         endif
!        endif
!        if(density_cg(m)<0.d0)density_cg(m)=0.d0
!        m=m+1
!       enddo ! i
!       mbase=mbase+ncgcell(1)
!      enddo ! j 
!     enddo ! k 
!     m=1 ; mbase=1
!     do k=1,ncgcell(1)
!      do j=1,ncgcell(1)
!       do i=1,ncgcell(1)
!        if(psi_squared_cg(m)==-1.d0)then
!         ii=max(mbase,m-5) ; jj=min(mbase+ncgcell(1)-1,m+5)
!         n=jj-ii ; nn=-1 ; nnn=0
!         do kk=1,n
!          if(ii+kk-1==m)nn=0
!          if(psi_squared_cg(ii+kk+nn)/=-1.d0)then
!           nnn=nnn+1
!           points(nnn)=latpos2_cg(1,ii+kk+nn)
!           repair(nnn)=psi_squared_cg(ii+kk+nn)
!          endif
!         enddo
!         if(nnn>1)then
!          call interp_nev(points(1:nnn),repair(1:nnn),nnn,latpos2_cg(1,m),&
!           &val_interp,err_interp)
!          psi_squared_cg(m)=val_interp
!         endif
!        endif
!        if(psi_squared_cg(m)<=0.d0)psi_squared_cg(m)=1.d-5
!        m=m+1
!       enddo ! i
!       mbase=mbase+nc
!      enddo ! j 
!     enddo ! k 
!    endif ! interpolation_required


!      if(plot_h_integrand)then
!       filename='h_integrand_t='//trim(adjustl(tmpr2))//'.dat'
!       open(file=filename,unit=14,status='unknown')
!       write(14,'(a,f25.16)')'#XRANGE 0.0 ',cell_x
!       write(14,'(a,f25.16)')'#YRANGE 0.0 ',cell_y
!       write(14,'(a,f25.16)')'#ZRANGE 0.0 ',cell_z
!
!       write(14,'(f25.16,a)')0.d0,' 0. 0. 0.'
!       do i=1,ncgcell(1)
!        write(14,'(f25.16,a)')latpos2_cg(1,i),' 0. 0. 0.'
!       enddo
!       write(14,'(f25.16,a)')cell_x,' 0. 0. 0.'
!
!       n=1
!       do j=1,ncgcell(1)
!        t1=latpos2_cg(2,n)
!        write(14,*) ! gnuplot line separator
!        write(14,'(2f25.16,a,f25.16)')0.d0,t1,' 0.',0.d0
!        do i=1,ncgcell(1)
!         write(14,'(2f25.16,a,f25.16)')latpos2_cg(1,i),t1,' 0.',density_cg(n)
!         n=n+1
!        enddo
!        write(14,'(2f25.16,a,f25.16)')cell_x,t1,' 0.',0.d0
!       enddo
!       write(14,*) ! gnuplot line separator
!       write(14,'(2f25.16,a,f25.16)')0.d0,cell_y,' 0.',0.d0
!       do i=1,ncgcell(1)
!        write(14,'(2f25.16,a,f25.16)')latpos2_cg(1,i),cell_y,' 0.',0.d0
!       enddo
!       write(14,'(2f25.16,a,f25.16)')cell_x,cell_y,' 0.',0.d0
!
!       do k=1,ncgcell(1)
!        zcoord=(real(k,dp)-0.5d0)*real(ncgrain(1),dp)*cell_over_nlat
!        write(14,*) ! gnuplot line separator
!        write(14,'(2(f25.16,a))')0.d0,' 0. ',zcoord,' 0.'
!        do i=1,ncgcell(1)
!         write(14,'(2(f25.16,a))')latpos2_cg(1,i),' 0. ',zcoord,' 0.'
!        enddo
!        write(14,'(2(f25.16,a))')cell_x,' 0. ',zcoord,' 0.'
!        n=1
!        do j=1,ncgcell(1)
!         t1=latpos2_cg(2,n)
!         write(14,*) ! gnuplot line separator
!         write(14,'(4f25.16)')0.d0,t1,zcoord,0.d0
!         do i=1,ncgcell(1)
!          write(14,'(4f25.16)')latpos2_cg(1,i),t1,zcoord,density_cg(n)
!          n=n+1
!         enddo
!         write(14,'(4f25.16)')cell_x,t1,zcoord,0.d0
!        enddo
!        write(14,*) ! gnuplot line separator
!        write(14,'(4f25.16)')0.d0,cell_y,zcoord,0.d0
!        do i=1,ncgcell(1)
!         write(14,'(4f25.16)')latpos2_cg(1,i),cell_y,zcoord,0.d0
!        enddo
!        write(14,'(4f25.16)')cell_x,cell_y,zcoord,0.d0
!       enddo ! k
! 
!       write(14,'(2(f25.16,a))')0.d0,' 0. ',cell_z,' 0.'
!       do i=1,ncgcell(1)
!        write(14,'(2(f25.16,a))')latpos2_cg(1,i),' 0. ',cell_z,' 0.'
!       enddo
!       write(14,'(2(f25.16,a))')cell_x,' 0. ',cell_z,' 0.'
!       n=1
!       do j=1,ncgcell(1)
!        t1=latpos2_cg(2,n)
!        write(14,*) ! gnuplot line separator
!        write(14,'(4f25.16)')0.d0,t1,cell_z,0.d0
!        do i=1,ncgcell(1)
!         write(14,'(4f25.16)')latpos2_cg(1,i),t1,cell_z,density_cg(n)
!         n=n+1
!        enddo
!        write(14,'(4f25.16)')cell_x,t1,cell_z,0.d0
!       enddo
!       write(14,*) ! gnuplot line separator
!       write(14,'(4f25.16)')0.d0,cell_y,cell_z,0.d0
!       do i=1,ncgcell(1)
!        write(14,'(4f25.16)')latpos2_cg(1,i),cell_y,cell_z,0.d0
!       enddo
!       write(14,'(4f25.16)')cell_x,cell_y,cell_z,0.d0
!       close(14)
! 
!       if(am_master)then
!        write(o,*)'Written integrand of coarse-grained H function formula.'
!       endif
!      endif ! h_integrand

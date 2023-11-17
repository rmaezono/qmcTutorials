MODULE store
!--------------------------------------------!
! Module to store frequently used variables. !
!--------------------------------------------!
 USE dsp
 IMPLICIT NONE

 INTEGER ialloc
 INTEGER o                ! Output unit (7 by default)
 REAL(sp) max_cpu_time,max_real_time ! CPU/real time limits
 LOGICAL errstop_skip     ! controls skipping of errstop calls
 LOGICAL open_unit(99)    ! Are the corresponding IO units free?
 LOGICAL wout_inhibit_node ! 'wout(io=o)' can be inhibited on this node
 CHARACTER(20) output_file

! Constants
 REAL(dp),PARAMETER :: pi=3.141592653589793238d0,twopi=2.d0*pi,fourpi=4.d0*pi,&
  &one_over_pi=1.d0/pi,two_over_pi=2.d0/pi,one_over_pi_squared=1.d0/pi**2,&
  &pi_over_two=pi/2.d0,four_over_pi=4.d0/pi,four_over_pi_squared=4.d0/pi**2,&
  &eight_over_pi_cubed=8.d0/pi**3,& 
  &sixteen_over_pi_squared=four_over_pi_squared*4.d0,&
  &twentythree_over_twentyfour=23.d0/24.d0,seven_sixths=7.d0/6.d0
 REAL(dp),PARAMETER :: mineps=1.d-12 ! Min value of overall tolerance level

! Input parameters
 INTEGER dim,maxstp,nmodes,nmodes_dim,den_ntimes,nscgrain,nlattice,&
  &nlattice3,nsmoothstep,nlatticesq,dentype,phase_preset,weight_preset,&
  &ntrajectories,phase_noise,num_ncgrain,n3dplanes
 INTEGER,ALLOCATABLE :: ncgrain(:),plane_number(:)
 REAL(dp) traj_time_start,traj_time_end,cell_x,cell_y,cell_z,den_time_start,&
  &den_time_end,init_eps,init_h,converge_maxdiff,curlweight
 REAL(dp),ALLOCATABLE,SAVE :: theta_1d(:),theta_2d(:,:),theta_3d(:,:,:),&
  &weights_1d(:),weights_2d(:,:),weights_3d(:,:,:),&
  &xstart(:),ystart(:),zstart(:)
 LOGICAL timing_info,verbose,negphase,transposephase,plot_all_traj,hfunction,&
  &testrun,plot_raw,plot_cg,plot_smooth,plot_h_integrand,writefail,&
  &nstep_histogram,read_backtracked,save_backtracked,fastmode
 CHARACTER(20) :: calc_type,int_algorithm,time_direction,wfn_type,vel_type,&
  &phase_format,weight_format

! Working parameters
 INTEGER nok,nbad,nsteps
 INTEGER :: max_maxstp=1000000000
 REAL(dp) eps,h1,hmin,dtsave,maxdiff,two_over_pi_n,four_over_pi2_n2,&
  &eight_over_pi3_n3,half_cell_x,half_cell_y,half_cell_z,mean_energy,&
  &variance,cvariance
 REAL(dp),POINTER :: tp(:),xp(:,:)
 REAL(dp),ALLOCATABLE :: latpos1(:),latpos1_cg(:),latpos1_smooth(:),&
  &latpos2(:,:),latpos2_cg(:,:),latpos2_smooth(:,:)
 LOGICAL good_traj,save_steps

END MODULE store

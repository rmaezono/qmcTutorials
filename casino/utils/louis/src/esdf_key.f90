MODULE esdf_key
!-------------------------------------------------------------------------!
! Module to hold keyword list for the 'input' file. This must be updated  !
! as new keywords are brought into existence or removed.                  !
!                                                                         !
! The 'label' is the label as used in calling the esdf routines           !
! 'typ' defines the type, with the following syntax. It is 3 characters   !
! long.                                                                   !
! The first indicates:                                                    !
!  I - integer                                                            !
!  S - single                                                             !
!  D - double                                                             !
!  P - physical                                                           !
!  T - string (text)                                                      !
!  E - defined (exists)                                                   !
!  L - boolean (logical)                                                  !
!  B - block                                                              !
! The second is always a colon (:)                                        !
! The third indicates the "level" of the keyword                          !
!  B - Basic                                                              !
!  I - Intermediate                                                       !
!  E - Expert                                                             !
!  D - Dummy                                                              !
!                                                                         !
! 'dscrpt' is a description of the variable. It should contain a (short)  !
! title enclosed between *! ... !*, and then a more detailed description  !
! of the variable.                                                        !
!                                                                         !
! MDT 12.2009                                                             !
!                                                                         !
! Changes                                                                 !
! -------                                                                 !
! None yet.                                                               !
!-------------------------------------------------------------------------!

IMPLICIT NONE
PRIVATE
PUBLIC load_keywords,kw_type,kw,numkw

TYPE kw_type
 CHARACTER(30) :: label
 CHARACTER(3) :: typ
 CHARACTER(2000) :: dscrpt
END TYPE kw_type

TYPE kwlist_item
 TYPE(kw_type) kw
 TYPE(kwlist_item),POINTER :: next
END TYPE kwlist_item

INTEGER :: numkw=0
TYPE(kw_type),ALLOCATABLE :: kw(:)
TYPE(kwlist_item),POINTER :: first_kw,current_kw,new_kw

CONTAINS


SUBROUTINE load_keywords
!---------------------------------------------------------------!
! Now define the keywords and load them into the keyword arrays !
!---------------------------------------------------------------!
 IMPLICIT NONE
 LOGICAL,SAVE :: first_call=.true.

 if(.not.first_call)return
 first_call=.false.

! Syntax:
!
! call add_kw("label","typ",&
!  &"*! Short title !*&
!  & Description in as many lines as required (up to 2000 characters), and&
!  & <-this space goes here, not at the end of the previous line.  This&
!  & symbol->$ &
!  & forces a newline, which is good for lists, etc.; spaces after the&
!  & newline are ignored.")
!
! --- START KEYWORD LIST ---

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! MDT NOTE: The keyword descriptions should be considered approximate/provisional 
! until this note is removed..
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 call add_kw("calc_type","T:B",&
  &"*! Single_trajectory/density !*&
  & Select type of calculation. Currently this can take the value: ;$&
  & 'single_trajectory' (propagate a single trajectory and write the coords&
  & to disk);$ & 
  & 'density' (calculate the density on a coarse-grained lattice at one or&
  & more times.)")

 call add_kw("dimensionality","I:B",&
  &"*! Number of space dimensions !*&
  & Select 1D/2D/3D trajectories. Currently the code can only do 2D.")

 call add_kw("int_algorithm","T:B",&
  &"*! Integration algorithm !*&
  & The code can use various numerical integration algorithms to propagate&
  & the trajectories. Currently the options are: ;$& 
  & 'runge-kutta' (The Runge-Kutta algorithm) ;$&
  & 'bulirsch-stoer (The Bulirsch-Stoer algorithm)")

 call add_kw("timing_info","L:B",&
  &"*! Timers on/off !*&
  & Toggle the production of timing info for the code. Note that activation&
  & of this option could conceivable slow down program execution.")

 call add_kw("verbose","L:B",&
  &"*! Toggle verbose output !*&
  & Produce more verbose output (whatever that means).")

 call add_kw("vel_type","T:B",&
  &"*! Velocity formula !*&
  & String to select the velocity formula to be used in computing the&
  & dynamics of the particle. Currently the options are 'deBB' - &
  & the standard de Broglie-bohm dynamics - or 'curl1'&
  & which corresponds to deBB plus the curl of psi^2 (i.e. the curl term&
  & in equation 19 of Colin and Struyve's paper (arXiv:0911.2823) with&
  & f=|\Psi^2|). Alternative velocity forumulae could easily be added.")

 call add_kw("curlweight","D:B",&
  &"*! Weight of vel eq. curl term !*&
  & Weight of curl term in velocity formula i.e. mu in Colin and Struye's&
  & paper arXiv:0911.2823. Meaningful only if vel_type='curl1'.")

 call add_kw("maxstp","I:B",&
  &"*! Max num integration iterations !*&
  & Maximum number of steps to be taken in the computation of a single&
  & trajectory over any period of 4pi seconds or fraction thereof. This is& 
  & likely to be exceeded for very long times, or for very tight tolerances&
  & (small EPS). Increasing this from its default value of 100000 will lead&
  & to fewer failed trajectories but a (possibly hugely) increased total CPU&
  & time and crappy load balancing on parallel machines. It has a maximum&
  & value MAX_MAXSTP which is defined in the store.f90 module to be 1000000000&
  & (as of 3/2010, MDT) - increase this at your own risk.")

 call add_kw("init_eps","D:B",&
  &"*! Initial value of EPS !*&
  & Initial value of EPS, the overall tolerance level eps for the Runge-Kutta&
  & integration algorithm. This value will be repeatedly decreased by a&
  & factor of ten and the trajectory recalculated until the separation of&
  & the final positions in successive iterations differs by less than&
  & INIT_MAXDIFF, or EPS gets reduced to 1.d-12, whichever happens sooner.&
  & A default value of 0.0001 was found to give the fastest overall&
  & calculations in limited tests.")

 call add_kw("init_h","D:B",&
  &"*! Guessed first stepsize h !*&
  & Initial value for the stepsize h along the computed trajectory. This&
  & value will be changed according to the local behaviour of the trajectory&
  & by the Runge-Kutta adaptive stepsize control algorithm.")

 call add_kw("converge_maxdiff","D:B",&
  &"*! Is trajectory converged? !*&
  & The overall tolerance level EPS for the Runge-Kutta integration algorithm&
  & is progressively decreased from its initial value of INIT_EPS by a&
  & successive factors of ten and the trajectory recalculated in each&
  & iteration. A trajectory is regarded as converged if the separation of&
  & its final positions from the final position in the previous iteration&
  & differs by less than CONVERGE_MAXDIFF, or EPS gets reduced to 1.d-12,&
  &  whichever happens sooner.")

 call add_kw("time_direction","T:B",&
  &"*! Forward or backward in time !*&
  & Indicate whether it is intended for a single trajectory calculation to&
  & travel forwards or backwards in time. Must match the order of &
  & TRAJ_TIME_START and TRAJ_TIME_END.")

 call add_kw("ntrajectories","I:B",&
  &"*! No. of trajectories to plot !*&
  & In trajectory mode, LOUIS will calculate full trajectories from&
  & TRAJ_TIME_START to TRAJ_TIME_END for each of NTRAJECTORIES different&
  & starting points. These starting points should be collected into&
  & NTRAJECTORIES lines of the TRAJECTORY_START block.")

 call add_kw("traj_time_start","D:B",&
  &"*! Initial time !*&
  & Initial time for single trajectory calculation.")

 call add_kw("traj_time_end","D:B",&
  &"*! Final time !*&
  & Final time for single trajectory calculation.")

 call add_kw("plot_all_traj","L:B",&
  &"*! Plot all trajectories !*&
  & If this keyword is set to T, then a trajectory plot is produced for&
  & all trajectories calculated with the different values of the tolerance&
  & parameter EPS - even the non-converged ones. ")

 call add_kw("dtsave","D:B",&
  &"*! Min delta t store results !*&
  & When plotting single trajectories, results are stored only at&
  & time intervals greater than DTSAVE. This saves time and memory, while&
  & still producing essentially the same plot visually. Note that this&
  & is basically pointless so I will remove this eventually.")

 call add_kw("phase_noise","I:B",&
  &"*! Phase noise in nth dp !*&
  & When doing multiple trajectory calculations (NTRAJECTORIES>1) one may&
  & optionally add random noise in the PHASE_NOISEth decimal place of the&
  & set of wave function phases on every trajectory after the first. Having&
  & PHASE_NOISE=0 is equivalent to turning this off. Clearly this is intended&
  & for the case when each of the trajectories in the TRAJECTORY_START list has&
  & the same starting position. Maximum value is 14.")

 call add_kw("trajectory_start","B:B",&
  &"*! Starting positions block !*&
  & In trajectory mode, LOUIS will calculate full trajectories from&
  & TRAJ_TIME_START to TRAJ_TIME_END for each of NTRAJECTORIES different&
  & starting points. These starting points should be collected into&
  & NTRAJECTORIES lines of the TRAJECTORY_START block. In dimensions less&
  & than three the redundant coordinates will be ignored, irrespective&
  & of whether you include them.")

 call add_kw("den_ntimes","I:B",&
  &"*! No. of density plots !*&
  & In density mode, LOUIS will calculate a time-ordered sequence of plots&
  & of the particle density and the wave function, starting at DEN_TIME_START&
  & and finishing at DEN_TIME_END. The DEN_NTIMES keyword gives the number of&
  & such plots.")

 call add_kw("dentype","I:B",&
  &"*! Selection of initial density !*&
  & The integer keyword DENTYPE selects a particular initial density. The&
  & different types available should be documented in the manual, but if not&
  & see the eval_density_to_Xd routine in SIN_WFN.f90 (where X is the&
  & dimensionality of the system).")

 call add_kw("cell_x","D:B",&
  &"*! Size of cell: x !*&
  & Gives the length of the simulation cell along the x direction.")

 call add_kw("cell_y","D:B",&
  &"*! Size of cell: y !*&
  & Gives the length of the simulation cell along the y direction. &
  & Currently, in 2D the cell is assumed to be a square and in 3D the&
  & cell is assumed to be a cube, and thus CELL_Y must equal&
  & CELL_X. This restriction will be removed in future.")

 call add_kw("cell_z","D:B",&
  &"*! Size of cell: z !*&
  & Currently, in 2D the cell is assumed to be a square and in 3D the&
  & cell is assumed to be a cube, and thus CELL_Z must equal&
  & CELL_X and CELL Y. This restriction may be removed in future.")

 call add_kw("nlattice","I:B",&
  &"*! No. of lattice points !*&
  & Number of lattice points per dimension in the box on which the raw&
  & density/wave function is explicitly calculated. NLATTICE is constrained to&
  & be an integer multiple of NCGRAIN and NSCGRAIN, and has a default value of&
  & 1024. ")

 call add_kw("nscgrain","I:B",&
  &"*! No. of points in smoothed cg cell !*&
  & The number of points per dimension in an overlapping smooth coarse-graining&
  & cell. The total number of points per dimension in the entire box (NLATTICE)&
  & must be an integer multiple of this. The smoothed coarse-grained density&
  & is the average of the densities at all points within such a cell for which&
  & backtracking could be carried out. Each subsequent point at which the&
  & smoothed density is calculated is found by translating the smoothing cell&
  & by NSMOOTHSTEP points. In 3D where smoothed calculations are not&
  & possible, this keyword is ignored. The default value of NSCGRAIN is 64.")

 call add_kw("nsmoothstep","I:B",&
  &"*! Step in smoothed coarse-graining !*&
  & In the smooth coarse-graining procedure, the smoothed coarse-grained&
  & density is the average of the densities at all NSCGRAIN**dim points within&
  & an overlapping smooth coarse-graining cell for which backtracking&
  & could be carried out. Each subsequent point at which the smoothed density&
  & is calculated is found by translating the smoothing cell by NSMOOTHSTEP&
  & points (to generate all possible cells which fit in the box in all DIM&
  & dimensions). You must therefore ensure that&
  & (NLATTICE-NSCGRAIN)/NSMOOTHSTEP + 1 is a whole number to make sure the&
  & translated cells fit in the box. The keyword is ignored in 3D where&
  & smoothed calculations are not possible. The default value of NSCGRAIN&
  & is 8.")

 call add_kw("den_time_start","D:B",&
  &"*! Time for first density plot !*&
  & In density mode a sequence of DEN_NTIMES density plots at specific times&
  & will be created. DEN_TIME_START is the initial time in this sequence.")

 call add_kw("den_time_end","D:B",&
  &"*! Time for last density plot !*&
  & In density mode a sequence of DEN_NTIMES density plots at specific times&
  & will be created. DEN_TIME_END is the final time in this sequence.")

 call add_kw("hfunction","L:B",&
  &"*! Compute H-function !*&
  & Toggle the computation of the H function in density mode. This keyword is&
  & ignored in single-trajectory mode. Default F.")

 call add_kw("read_backtracked","L:B",&
  &"*! Read backtrack positions !*&
  & If READ_BACKTRACKED is T, then the backtracked positions at t=0 of all&
  & particles starting at each lattice point at some time t1 will be read&
  & from a file 'backtracked_positions.dat'. Such files may be produced by&
  & setting the SAVE_BACKTRACKED keyword in a previous run - the files will&
  & include the time in their filename (e.g. backtracked_positions_t=t1.dat)&
  & and must therefore be renamed before being read in. The purpose of this&
  & functionality is to be able to calculate density plots, H functions,&
  & relaxation times etc. with e.g. different coarse-graining lengths without&
  & having to recompute the backtracked positions (which is 95% of the cost of&
  & a calculation. The input file must be set up for a density calculation&
  & at a single time (DEN_NTIMES=1,DEN_TIME_START=DEN_TIME_END/=0) and this&
  & time must correspond to the 2nd line in the backtracked_positions.dat file.&
  & Note this facility is disabled in 3D because of the huge size of the file.")

 call add_kw("save_backtracked","L:B",&
  &"*! Save backtrack positions !*&
  & If SAVE_BACKTRACKED is T, then the backtracked positions at t=0 of all&
  & particles starting at each lattice point at some time t1 will be saved to a&
  & file 'backtracked_positions_t=t1.dat' (one file for each of the DEN_NTIMES&
  & starting points except t=0). This is not done by default. Note this&
  & facility is disabled in 3D because of the huge size of the file.")

 call add_kw("plot_raw","L:B",&
  &"*! Plot raw densities !*&
  & Write the raw density and psi squared at each time step to&
  & density_raw_t=xx.dat and psi_squared_t=xx.dat files. Note these files can&
  & be very large. Note that for 3D calculations, this data is written only&
  & for selected 2D planes (which may be specified in the 3DPLANES input&
  & block) in the files density_raw_t=xx_p=yy.dat etc.")

 call add_kw("plot_cg","L:B",&
  &"*! Plot coarse-grained densities !*&
  & Write the coarse-grained density and psi squared at each time step to&
  & density_cg_t=xx.dat and psi_squared_cg_t=xx.dat files. In 3D, the&
  & coarse_grained density is written out for all xy slabs of coarse-graining&
  & cells that contain a raw data plane specified in the 3DPLANE block (using&
  & the file name density_cg_t=xx_pc=yy.data etc.")

 call add_kw("plot_smooth","L:B",&
  &"*! Plot smoothed densities !*&
  & Write the smoothed coarse-grained density and psi squared at each time&
  & step to a density_smooth_t=xx.dat and psi_squared_smooth_t=xx.dat files.")

 call add_kw("plot_h_integrand","L:B",&
  &"*! Plot integrand of H function !*&
  & Write the integrand of the formula used to calculate the H function& 
  & for each coarse-graining cell. i.e. rho_cg * log ( rho_cg / psi^2_cg ).")

 call add_kw("writefail","L:B",&
  &"*! Write initpos failed traj !*&
  & Write the initial time and positions of all failed trajectories in a&
  & density calculation to a file 'failed_trajectories.dat' (in parallel&
  & separate files will be created for each node. These strange trajectories&
  & may then be studied in single-trajectory mode.")

 call add_kw("nstep_histogram","L:B",&
  &"*! Bin number of steps !*&
  & Collect into 1000 bins - scaled onto 1:MAXSTP - the total number of steps&
  & taken for each single trajectory in a density calculation. This can be&
  & plotted as a histogram.")

 call add_kw("fastmode","L:B",&
  &"*! Hard trajectories method !*&
  & Controls behaviour when LOUIS encounters 'hard trajectories'. These are&
  & trajectories where the number of steps required exceeds the maximum&
  & value MAXSTP, usually because the path is winding at high velocity&
  & round and round a moving node. If fastmode is set to T, no extra effort&
  & is expended in trying to converge these trajectories and problem&
  & trajectories are simply labelled as 'unable to be backtracked'. Any&
  & resulting gaps in the coarse-grained density and psi^2 data are filled in&
  & by an approximate interpolation procedure. If fastmode is set to F, then&
  & the code will repeatedly attempt to converge problem trajectories by&
  & using tricks, such as increasing maxstp, tightening the integration&
  & accuracy controls, loosening the accuracy tolerance on the final position,&
  & and other things. These additional measures can add considerably to the&
  & cost of the calculation, and it is still likely that a number of points&
  & will be impossible to backtrack without extremely long run times (the&
  & hopefully smaller number of gaps in the coarse-grained density data will &
  & still be filled in using approximate interpolation). Note the above&
  & assumes LOUIS is running in density mode, but the keyword also has&
  & the same significance in trajectory mode, except for the stuff&
  & about interpolation. Note also that for FASTMODE to work properly, then&
  & the MAXSTP parameter should be set somewhere near its default value of&
  & 10^5. To gauge the effect on the accuracy, look at the computed H-function&
  & values with and without FASTMODE = T.") 

 call add_kw("coarse_graining_lengths","B:B",&
  &"*! No. of points in cg cell !*&
  & A block defining a set of coarse-graining lengths to use in computing&
  & coarse-grained quantities such as the H function. The block contains: $&
  & LINE 1: The number of coarse-graining lengths in the set (NUM_NCGRAIN) ;$&
  & LINE 2: NUM_NCGRAIN values for the coarse-graining length NCGRAIN. ;$&
  & Each NCGRAIN is the number of points per dimension in a non-overlapping&
  & coarse-graining cell. The total number of points per dimension in the&
  & entire box (NLATTICE) must be an integer multiple of each of the set of&
  & values of NCGRAIN. The coarse-grained density (for example) is the&
  & average of the densities at all points within the coarse-graining cell for&
  & which backtracking could be carried out. A typical value for NCGRAIN is&
  & 32.")

 call add_kw("3dplanes","B:B",&
  &"*! Which planes to plot in 3D !*&
  & In 3 dimensions, one can choose a set of xy planes at different z values&
  & for plotting quantities such as the density. The planes can be manually&
  & selected in this input block. The block contains: $&
  & LINE 1: the number of planes to plot (N3DPLANES) ;$&
  & LINE 2: N3DPLANES selected planes (which have integer sequence numbers&
  & from 1:NLATTICE). $&
  & Note that because of the way the algorithm is parallelized it is not&
  & currently possible to plot data from planes other than xy. The maximum&
  & number of planes is hardwired to 10 - clearly this can be easily modified.")

 call add_kw("wfn_type","T:B",&
  &"*! Type of wave function !*&
  & Currently only 'sine_wave' allowed. Proper keyword description to be&
  & inserted.")

 call add_kw("nmodes","I:B",&
  &"*! Number of modes !*&
  & Number of eigenstates in the superposition making up the wave function.")

 call add_kw("phase_format","T:B",&
  &"*! How to choose phases !*&
  & This keyword determines how to choose the set of NMODES initial t=0 phases&
  & for the time-dependent complex exponentials in the superposition making&
  & up the wave function. The various options are ;$&
  & 'input' - read the NMODES phases directly from the PHASES block in the &
  & input file. The PHASES block is ignored if PHASE_FORMAT has a value&
  & other than 'input'.;$&
  & 'random' - use a truly random set of phases with the seed derived from the&
  & system timer.;$&
  & 'default' - use a fixed set of phases generated from a predefined seed&
  & 31415927.$&
  & 'preset' - LOUIS will use the Nth&
  & random number in a preset sequence as the seed for the generation of&
  & a set of random phases; N is the value of PHASE_PRESET. The point of&
  & this is supposed to be that you can do a *reproducible* set of (say) 10&
  & calculations which differ only in the set of phases used. Each set of&
  & NMODES phases is defined uniquely by the single integer PHASE_PRESET,&
  & which might run from 1-10 in this case to define the 10 sets of&
  & calculations. Note that PHASE_FORMAT='preset'/PHASE_PRESET=0 is&
  & entirely equivalent to PHASE_FORMAT='default'.")

 call add_kw("phase_preset","I:B",&
  &"*! Which set of preset phases !*&
  & If PHASE_FORMAT has the value 'preset' then LOUIS will use the Nth&
  & random number in a preset sequence as the seed for the generation of&
  & a set of random phases; N is the value of PHASE_PRESET. The point of&
  & this is supposed to be that you can do a *reproducible* set of (say) 10&
  & calculations which differ only in the set of phases used. Each set of&
  & NMODES phases is defined uniquely by the single integer PHASE_PRESET,&
  & which might run from 1-10 in this case to define the 10 sets of&
  & calculations.")

 call add_kw("weight_format","T:B",&
  &"*! How to choose weights !*&
  & See 'PHASE_FORMAT' for a detailed description.")

 call add_kw("weight_preset","I:B",&
  &"*! Set of preset amplitudes !*&
  & Choose a reproducible set of random amplitudes to use for the wave&
  & function. See PHASE_PRESET for a detailed description.")

 call add_kw("negphase","L:B",&
  &"*! Reverse sign of phases!*&
  & Reverse sign of all the phases in the phases block. Option included&
  & since both Westman and Colin coded up the wrong sign of the phases in&
  & their programs.")

 call add_kw("transposephase","L:B",&
  &"*! Transpose phase matrix!*&
  & Transpose the phase matrix read in from the PHASES block. Only relevant in&
  & two dimensions. Should flip x and y in a trajectory plot. Used as a&
  & consistency check.")

 call add_kw("phases","B:B",&
  &"*! Number of modes !*&
  & Block containing NMODES phases to define the starting wave function; these&
  & are the initial t=0 phases of the complex exponentials in each of the&
  & eigenfunctions in the superposition making up the wave function. This&
  & block of phases is ignored unless PHASE_FORMAT has the value 'input'. The&
  & numbers should be one per line and in Fortran column major order i.e. in&
  & 2D if you have a conventional matrix, write it into the phases block by&
  & proceeding down the first column, then the second column etc.")

 call add_kw("weights","B:B",&
  &"*! Amplitudes of modes !*&
  & Block containing NMODES amplitudes to define the starting wave function;&
  & these are the coefficients of the eigenfunctions making up the wave&
  & function. This block is ignored unless PHASE_FORMAT has the value 'input'&
  & and WFN_TYPE has the value 'scaled_sin'. The numbers should be one per&
  & line and ordered in the same way as the PHASES block")

 call add_kw("testrun","L:B",&
  &"*! Test run!*&
  & Do a test run, rather than a full production run. Intended for debugging&
  & or profiling. What a test run actually involves is up to the programmer.&
  & It is highly likely that when you try it, it won't do anything at all.")

! PLOT KEYWORDS

 call add_kw("plot_output","T:B",&
  &"*! Type of output for the plot !*&
  & Select the output type for the plot. For popout windows, choose 'wxt' or&
  & 'x11'; 'wxt' allowing interactive zoom level. For files, choose from &
  &'gif', 'jpg', 'png', 'svg', 'postscript', 'eps', 'pdf'. Option 'agif' will&
  & produce an animated gif, with each datafile being plotted as a frame in&
  & the animation. Note that availability of output types depends on your&
  & build of gnuplot.")

 call add_kw("plot_title","T:B",&
  &"*! Set a title for the plot !*&
  & Set the title for the plots. Allows three special options: 'auto' will&
  & automatically generate a title if the datafile is an output from LOUIS;&
  & 'none' will not show a title; 'prompt' will prompt the user for a title for&
  & each datafile. Anything else will be interpreted as a string to use as the&
  & title for all datafiles.")

 call add_kw("plot_type","T:B",&
  &"*! Type of plot to produce !*&
  & Select the style of plot to produce: $&
  & 'trajectory' shows the position of a particle in time for 1D (or 2D&
  & with 'spacetime' on) or the path of a single particle for 3D or 2D with&
  & 'spacetime' off;$& 
  &'density' shows rho (or |psi|^2) at all positions in the box (as a line&
  & in 1D or a surface in 2D);$&
  & 'bar' is for plotting coarse-grained data in 1D or 2D and shows flat bars&
  & indicating the value of rho (or |psi|^2) in each coarse-graining cell;$&
  & 'map' can only be used in 2D and is for showing rho (or |psi|^2) on&
  & a 2D plane either as contours or represented by the colour at each point.$&
  & 'sweep' and 'sweep3d' are to be used for plotting 3D density graphs. They &
  & both require plot output to be one of 'agif','postscript','pdf'")

 call add_kw("axislabels","L:B",&
  &"*! Add labels to the axes !*&
  & If T, labels will be added to x,y,z and t axes. No label will be added&
  & to the |psi|^2 or rho axes, as either could be plotted.")

 call add_kw("spacetime","L:B",&
  &"*! 2D trajectory in spacetime !*&
  & Display 2D trajectories with time as the third axis.$&
  & NOT CURRENTLY IMPLEMENTED.")

 call add_kw("framerate","I:B",&
  &"*! in fps (agif only) !*&
  & Set the framerate (in frames per second) for animated gifs.")

 call add_kw("plane","T:B",&
  &"*! 3D density plane constant !*& 
  & When plotting 3D densities using the 'sweep' plot type, take cross-&
  &sections with this axis co-ordinate constant, e.g. plane 'x' will plot&
  & cross-sections of constant x")

 call add_kw("autoscale_x","L:B",&
  &"*! Autoscale the x axis !*&
  & Automatically scale the x axis based on the range of data in the datafile.&
  & If multiple datafiles are supplied, the user will be prompted as to&
  & whether to use the same scaling for all datafiles, or scale them&
  & independently.&
  & Note that this doesn't necessarily autoscale the tick placing. See 'xtics'&
  & for more detail.")

 call add_kw("xmin","D:B",&
  &"*! Minimum for the x-axis !*&
  & Lowest x value to be plotted. Not used if autoscaling.")

 call add_kw("xmax","D:B",&
  &"*! Maximum for the x-axis !*&
  & Highest x value to be plotted. Not used if autoscaling.")

 call add_kw("xticks","D:B",&
  &"*! Spacing of ticks on x-axis !*&
  & Spacing of the ticks (and labels) on the x axis. 3 options for this:$&
  & positive value - tick spacing always set to this value;$&
  & zero - ticks always off;$&
  & negative value - tick spacing set automatically by gnuplot.$&
  & Note that this is independent of whether the axes are autoscaled.")

 call add_kw("autoscale_y","L:B",&
  &"*! Autoscale the y axis !*&
  & Automatically scale the y axis based on the range of data in the datafile.&
  & If multiple datafiles are supplied, the user will be prompted as to&
  & whether to use the same scaling for all datafiles, or scale them&
  & independently.&
  & Note that this doesn't necessarily autoscale the tick placing. See 'xtics'&
  & for more detail.")

 call add_kw("ymin","D:B",&
  &"*! Minimum for the y-axis !*&
  & Lowest y value to be plotted. Not used if autoscaling.")

 call add_kw("ymax","D:B",&
  &"*! Maximum for the y-axis !*&
  & Highest y value to be plotted. Not used if autoscaling.")

 call add_kw("yticks","D:B",&
  &"*! Spacing of ticks on y-axis !*&
  & Spacing of the ticks (and labels) on the y axis. 3 options for this:$&
  & Positive value - tick spacing always set to this value;$&
  & zero - ticks always off;$&
  & negative value - tick spacing set automatically by gnuplot.$&
  & Note that this is independent of whether the axes are autoscaled.")

 call add_kw("autoscale_z","L:B",&
  &"*! Autoscale the z axis !*&
  & Automatically scale the z axis based on the range of data in the datafile.&
  & If multiple datafiles are supplied, the user will be prompted as to&
  & whether to use the same scaling for all datafiles, or scale them&
  & independently.&
  & Note that this doesn't necessarily autoscale the tick placing. See 'xtics'&
  & for more detail.")

 call add_kw("zmin","D:B",&
  &"*! Minimum for the z-axis !*&
  & Lowest z value to be plotted. Not used if autoscaling.")

 call add_kw("zmax","D:B",&
  &"*! Maximum for the z-axis !*&
  & Highest z value to be plotted. Not used if autoscaling.")

 call add_kw("zticks","D:B",&
  &"*! Spacing of ticks on z-axis !*&
  & Spacing of the ticks (and labels) on the z axis. 3 options for this:$&
  & Positive value - tick spacing always set to this value;$&
  & zero - ticks always off;$&
  & negative value - tick spacing set automatically by gnuplot.$&
  & Note that this is independent of whether the axes are autoscaled.")

 call add_kw("autoscale_t","L:B",&
  &"*! Autoscale the t axis !*&
  & Automatically scale the t axis based on the range of data in the datafile.&
  & If multiple datafiles are supplied, the user will be prompted as to&
  & whether to use the same scaling for all datafiles, or scale them&
  & independently.&
  & Note that this doesn't necessarily autoscale the tick placing. See 'xtics'&
  & for more detail.")

 call add_kw("tmin","D:B",&
  &"*! Minimum for the t-axis !*&
  & Lowest t value to be plotted. Not used if autoscaling.")

 call add_kw("tmax","D:B",&
  &"*! Maximum for the t-axis !*&
  & Highest t value to be plotted. Not used if autoscaling.")

 call add_kw("tticks","D:B",&
  &"*! Spacing of ticks on t-axis !*&
  & Spacing of the ticks (and labels) on the t axis. 3 options for this:$&
  & Positive value - tick spacing always set to this value;$&
  & zero - ticks always off;$&
  & negative value - tick spacing set automatically by gnuplot.$&
  & Note that this is independent of whether the axes are autoscaled.")

 call add_kw("colour","L:B",&
  &"*! Use colour for plots !*&
  & If true, plots will have coloured lines and surfaces. If false,&
  & greyscale is used where possible, but some output types will still&
  & have coloured lines.")

 call add_kw("color","L:B",&
  &"*! Use color for plots !*&
  & If true, plots will have colored lines and surfaces. If false,&
  & greyscale is used where possible, but some output types will still&
  & have colored lines.")

 call add_kw("linegradient","L:B",&
  &"*! Line colour indicates z value !*&
  & If true, the contours in 2D map plots and the mesh in 2D density plots&
  & will be coloured according to the value of rho (or |psi|^2) at that&
  & point. Lines will default to black if the shading is also on so that&
  & they are visible.")

 call add_kw("tickslevel","D:B",&
  &"*! Gap xy plane -> z-axis zero !*&
  & Gap between the xy plane and the zero of the z-axis.")

 call add_kw("autorotate","L:B",&
  &"*! Use gnuplot default angles !*&
  & If true, the gnuplot default viewing angle is used. This corresponds to&
  & 60,30 in the manual settings (see 'rot_x' and 'rot_z').")

 call add_kw("rot_x","D:B",&
  &"*! Rotation about x-axis !*&
  & Looking directly onto the x-y plane with the x-axis horizontal,&
  & this is the degree of rotation about the x-axis. This is applied&
  & before 'rot_z'.")

 call add_kw("rot_z","D:B",&
  &"*! Rotation about x-axis !*&
  & After applying 'rot_x', this is the degree of rotation about the z-axis.")

 call add_kw("colour_key","L:B",&
  &"*! Display a key to colours !*&
  & In plots where colour of a line or surface encodes the value of rho&
  & (or |psi|^2), this will toggle a key relating the colour to a numerical&
  & value.")

 call add_kw("color_key","L:B",&
  &"*! Display a key to colors !*&
  & In plots where color of a line or surface encodes the value of rho&
  & (or |psi|^2), this will toggle a key relating the color to a numerical&
  & value.")

 call add_kw("shading","L:B",&
  &"*! Toggle shading !*&
  & For 2D density plots, this will shade the surface according to the value&
  & of rho; for a 2D bar graph the bars will be shaded, and for a 2D map&
  & plot the shading will be on the x-y plane. No effect on any other type&
  & of plot.")

 call add_kw("mesh","L:B",&
  &"*! Show a mesh on a surface !*&
  & For 2D density plots, show a mesh on top of the surface; for a 2D bar&
  & graph, this will be an outline around the bars. No effect on any other&
  & type of plot.")

 call add_kw("contours","L:B",&
  &"*! Show contours on a map !*&
  & Show contours on a 2D map plot. No effect on any other type of plot.")

 call add_kw("contour_no","I:B",&
  &"*! Number of contours !*&
  & Tell gnuplot how many contours to display. It sometimes actually uses&
  & that number.")

 call add_kw("contour_labels","L:B",&
  &"*! Display contour labels !*&
  & Show contour labels on a 2D map plot. No effect on any other type of plot.")

! --- END OF KEYWORD LIST ---

! Build the kw(:) vector that esdf.f90 uses.
 call build_kw_vector

END SUBROUTINE load_keywords


SUBROUTINE add_kw(label,typ,dscrpt)
 IMPLICIT NONE
 CHARACTER(*),INTENT(in) :: label
 CHARACTER(*),INTENT(in) :: typ
 CHARACTER(*),INTENT(in) :: dscrpt
 allocate(new_kw)
 new_kw%kw%label=trim(adjustl(label))
 new_kw%kw%typ=trim(adjustl(typ))
 new_kw%kw%dscrpt=trim(adjustl(dscrpt))
 if(numkw>0)then
  current_kw%next=>new_kw
 else
  first_kw=>new_kw
 endif
 current_kw=>new_kw
 current_kw%next=>first_kw
 numkw=numkw+1
END SUBROUTINE add_kw


SUBROUTINE build_kw_vector
 IMPLICIT NONE
 INTEGER ikw
 allocate(kw(numkw))
 current_kw=>first_kw
 do ikw=1,numkw-1
  kw(ikw)=current_kw%kw
  new_kw=>current_kw%next
  deallocate(current_kw)
  current_kw=>new_kw
 enddo
 kw(numkw)=current_kw%kw
 deallocate(current_kw)
END SUBROUTINE build_kw_vector


END MODULE esdf_key

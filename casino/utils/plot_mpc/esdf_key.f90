module esdf_key
implicit none
private
public load_keywords,kw_type,kw,numkw
type kw_type
character(30) :: label
character(3) :: typ
character(2100) :: dscrpt
end type kw_type
type kwlist_item
type(kw_type) kw
type(kwlist_item),pointer :: next
end type kwlist_item
integer :: numkw=0
type(kw_type),allocatable :: kw(:)
type(kwlist_item),pointer :: xyzzyaaaa1,xyzzyaaab1,xyzzyaaac1
contains
subroutine load_keywords
implicit none
logical,save :: xyzzyaaaa2=.true.
if(.not.xyzzyaaaa2)return
xyzzyaaaa2=.false.
call xyzzyaaad1("neu","I:B","*! Number of up electrons !* For real sys&
&tems containing atoms, NEU is the total number of spin-up electrons r&
&eferenced by the many-body wave function (for periodic systems, this &
&is the number of spin-up electrons in the simulation cell, rather tha&
&n the underlying primitive cell). The number of spin-down electrons i&
&s given by the keyword NED.$$ Note that in the presence of addition o&
&r subtraction excitations, NEU refers to the state of the system AFTE&
&R the required number of electrons have been added or removed.  For m&
&odel electron(-hole) phases such as the HEG, set NEU to zero and use &
&the FREE_PARTICLES block to define the number of spin-up electrons.")
call xyzzyaaad1("ned","I:B","*! Number of down electrons !* For real s&
&ystems containing atoms, NED is the total number of spin-down electro&
&ns referenced by the many-body wave function (for periodic systems, t&
&his is the number of spin-down electrons in the simulation cell, rath&
&er than the underlying primitive cell). The number of spin-up electro&
&ns is given by the keyword NEU.$$ Note that in the presence of additi&
&on or subtraction excitations, NED  refers to the state of the system&
& AFTER the required number of electrons have been added or removed. F&
&or model electron(-hole) phases such as the HEG, set NED to zero and &
&use the FREE_PARTICLES block to define the number of spin-down electr&
&ons.")
call xyzzyaaad1("nhu","I:I","*! Number of spin-up fermions other than &
&electrons in real systems !* NHU is the number of spin-up fermions ot&
&her than electrons in real systems.  For example, if you are interest&
&ed in positronic molecules then NHU should be set to 1, and the up-sp&
&in positron (""spin"" 3) should be defined appropriately in the PARTI&
&CLES block.  Likewise for muonic systems.  At present only the Gaussi&
&an routines can be used to return orbitals for species other than ele&
&ctrons.  (There exists a modified version of the GAUSSIAN code that c&
&an be used to generate orbitals for positronic molecules.)  If you wa&
&nt to use a different basis or study muons or something then you will&
& need to make the appropriate changes first.  For model systems such &
&as electron-hole gases, etc., please use the FREE_PARTICLES block to &
&define the number of spin-up holes.")
call xyzzyaaad1("nhd","I:I","*! Number of spin-down fermions other tha&
&n electrons in real systems !* Like NHU, but for spin-down particles.&
&")
call xyzzyaaad1("runtype","T:B","*! Type of calculation !* This keywor&
&d specifies the type of QMC run to be carried out.  It can take the f&
&ollowing values:$ 'vmc' (perform a single VMC simulation);$ 'dmc_equi&
&l' (perform DMC equilibration);$ 'dmc_stats' (perform DMC statistics &
&accumulation);$ 'dmc_dmc' (perform DMC equilibration, then statistics&
& accumulation);$ 'vmc_dmc' (perform VMC, then DMC equilibration, then&
& DMC statistics accumulation);$ 'vmc_dmc_equil' (perform VMC, then DM&
&C equilibration);$ 'opt' (perform a single optimization run);$ 'vmc_o&
&pt' (perform OPT_CYCLES cycles of VMC and optimization, alternately);&
&$ 'opt_vmc' (perform OPT_CYCLES cycles of optimization and VMC, alter&
&nately);$ 'gen_mdet_casl' (generate an mdet.casl file for post-proces&
&sing by det_compress);$ 'gen_mpc' (generate an mpc.data file enabling&
& the use of the MPC interaction);$ 'rmc' (perform a single reptation &
&simulation: EITHER equilibration OR statistics accumulation -- EXPERI&
&MENTAL);$ 'rmc_rmc' (perform reptation equilibration, then reptation &
&statistics accumulation -- EXPERIMENTAL);$ 'plot' (perform plot speci&
&fied by either block PLOT or block QMC_PLOT).;$ NOTE: in earlier vers&
&ions of the code, we used 'runtype : dmc' with the now-redundant keyw&
&ord 'iaccumulate : T or F' to indicate whether stats accumulation was&
& activated or not. This usage is now deprecated and - unless iaccumul&
&ate is specifically defined in input then 'runtype : dmc' is just a s&
&ynonym for 'runtype : dmc_dmc'.")
call xyzzyaaad1("atom_basis_type","T:B","*! Basis set type !* ATOM_BAS&
&IS_TYPE selects the basis set in which the atom-centred orbitals are &
&expanded (thus choosing which file to read the orbitals from), or mor&
&e generally, the 'type of orbital' to be used.  Possible values are:$&
&$ * 'none'       : [default] no atoms are present, and therefore no a&
&tomic orbitals are read in;$$ * 'plane-wave' : use a plane-wave basis&
& set; the orbitals are read in from pwfn.data;$$ * 'gaussian'   : use&
& a Gaussian basis set; the orbitals are read in from gwfn.data;$$ * '&
&slater-type': use a Slater-type orbital basis set; the orbitals are r&
&ead in from stowfn.data;$$ * 'numerical'  : use orbitals tabulated on&
& a grid (atomic systems only); the orbitals are read in from awfn.dat&
&a;$$ * 'dimer    '  : use orbitals tabulated on a grid (molecular dim&
&ers only); the orbitals are read in from dwfn.data;$$ * 'blip'       &
&: use a blip basis set; the orbitals are read in from bwfn.data;$$ So&
&me special wave function types are also available:$$ * 'nonint_he'  :&
& use exact orbitals for a non-interacting Helium atom.$$ * 'h2' : wav&
&e function for the H2 molecule where each orbital is the sum over hyd&
&rogen nuclei of a parameter-less exponential centred at each nucleus.&
&$$ * 'h3plus' : wave function for the H3+ molecular ion where each or&
&bital is the sum over hydrogen nuclei of a parameter-less exponential&
& centred at each nucleus.$$ For free-particle and external-potential-&
&related orbitals, set ATOM_BASIS_TYPE to 'none' and use the input blo&
&ck FREE_PARTICLES.")
call xyzzyaaad1("dmc_method","I:I","*! DMC method !* DMC_METHOD select&
&s which version of DMC to use: $ 1 = electron-by-electron algorithm;$&
& 2 = configuration-by-configuration algorithm.$ Method 1 is the defau&
&lt.")
call xyzzyaaad1("vmc_equil_nstep","I:B","*! Number of VMC equilibratio&
&n steps !* Total number of equilibration steps in VMC; should normall&
&y be at least a few thousand. Equilibration is only performed if NEWR&
&UN=T, thus the value of VMC_EQUIL_NSTEP is ignored on restarts. This &
&is a single-processor quantity, that is, all processors run VMC_EQUIL&
&_NSTEP equilibration steps. The value of VMC_DECORR_PERIOD is ignored&
& during equilibration.")
call xyzzyaaad1("vmc_nstep","I:B","*! Number of (main) VMC steps !* To&
&tal number of VMC steps summed over all processors; this corresponds &
&to the total number of particle configurations for which the energy a&
&nd other quantities to be averaged are calculated (required only if S&
&TOP_METHOD = 'nstep').$$ Note that because adjacent moves are likely &
&to be serially correlated, there is also an inner decorrelation loop &
&of length VMC_DECORR_PERIOD, so the total number of configuration mov&
&es attempted in a VMC run following equilibration is  VMC_NSTEP*VMC_D&
&ECORR_PERIOD.$$ On parallel machines, each core will do the same numb&
&er of steps and for each step the energy is averaged over the cores a&
&nd written to the vmc.hist file (which will ultimately contain VMC_NS&
&TEP/NCORES lines - though VMC_AVE_PERIOD adjacent lines may be averag&
&ed over to reduce the file size). This means that if VMC_NSTEP is not&
& divisible by the number of cores then it will internally be rounded &
&up to the nearest multiple of the number of cores (example: on a 12-c&
&ore machine, given VMC_NSTEP=20 in input, CASINO will round up VMC_NS&
&TEP to 24; each core will then do two steps and a total of two record&
&s will be written to vmc.hist, each of which is an average of 12 ener&
&gies). On a single-core machine with VMC_NSTEP=20, CASINO will move t&
&he single config 20 times, and 20 records will be written to vmc.hist&
&.$$ Note the VMC_NBLOCK or BLOCK_TIME keywords may be used to vary th&
&e frequency with which checkpointing is done i.e. how often we write &
&the data to disk; this does not affect the total number of VMC steps &
&and expectation values such as average energy should be independent o&
&f it.$$ Instead of executing a fixed number of VMC steps, CASINO may &
&be requested to do as many steps as are required to attain a target e&
&rror bar or to make the error bar as small as possible - see keyword &
&STOP_METHOD.")
call xyzzyaaad1("vmc_nblock","I:I","*! Number of blocks in VMC !* Sett&
&ing VMC_NBLOCK is one of two ways of specifying the number of blocks &
&into which the total VMC run is divided post-equilibration (the other&
& way being to specify BLOCK_TIME). The number of blocks determines ho&
&w often the output, history and configuration/checkpoint files are wr&
&itten to disk. More specifically, at the end of each block:$$ (1) the&
& processor- and block-averaged energies and a short 'report' are writ&
&ten to out.$$ (2) the processor-averaged energies for each step in th&
&e current block are appended to vmc.hist (and other quantities to exp&
&val.data).$$ (3) the current VMC state plus any accumulated configs a&
&re written to the config.out file (this latter only if the CHECKPOINT&
& input keyword is increased to 2 from its default value of 1 - otherw&
&ise config.out is only written after the end of the final block).$$ N&
&ote that the total energy and error bar should be effectively indepen&
&dent of VMC_NBLOCK (provided it is ensured that the random number seq&
&uence is independent of the number of blocks - which it has not been &
&at various periods in CASINO's history, though it should be now). Not&
&e that the value of VMC_NBLOCK is ignored if VMC_NTWIST>0, or if BLOC&
&K_TIME>0.0. The default value of VMC_NBLOCK is 1.")
call xyzzyaaad1("periodic","L:B","*! Periodic boundary conditions !* P&
&ERIODIC should be T if the system is periodic in 1, 2 or 3 dimensions&
&.  It should be F if system is finite without periodic boundary condi&
&tions.")
call xyzzyaaad1("testrun","L:B","*! Test run flag !* If TESTRUN is T, &
&CASINO will just read in all input files, print out the input informa&
&tion, then stop without performing any QMC calculation.")
call xyzzyaaad1("newrun","L:B","*! New run or continue old !* NEWRUN d&
&etermines whether this is a new run or a continuation of a previous o&
&ne.$$ * VMC:$ NEWRUN=T: VMC_EQUIL_NSTEP Metropolis steps are performe&
&d on a set of  randomly generated configs before accumulation of stat&
&istics begins.$ NEWRUN=F: a set of old (and presumably equilibrated) &
&electron positions are read from a config.in file, and no Metropolis &
&equilibration steps are performed before accumulation of statistics.$&
& * DMC:$ NEWRUN=T: a set of VMC configs is read from a config.in file&
& and the initial best estimate of the energy EBEST is calculated as t&
&he mean energy of these configs. In the special case of DMC_REWEIGHT_&
&CONFIGS=T, config.in may also contain DMC configs which are reweighte&
&d to a new wave function and EBEST shifted by the current energy diff&
& between the wave functions.$ NEWRUN=F: a set of DMC configs is read &
&from a config.in file and EBEST is not recomputed but taken to have t&
&he value written on the end of the config.in file (presumably by a pr&
&evious DMC run - either equilibration or accumulation).")
call xyzzyaaad1("density","L:B","*! Accumulate density !* IF DENSITY i&
&s set to T the charge density is accumulated and written to the expva&
&l.data file.")
call xyzzyaaad1("dtvmc","D:B","*! VMC time step !* DTVMC is the time s&
&tep for VMC runs in atomic units. In systems with more than one 'fami&
&ly' of particles, the time step used is the value of DTVMC divided by&
& the particle mass. If OPT_DTVMC is set to 1, the time step is optimi&
&zed to give acceptance ratios of about 50% (this is done for each 'fa&
&mily' individually). Both DTVMC and OPT_DTVMC are ignored if input bl&
&ock DTVMCS is present.")
call xyzzyaaad1("dtdmc","D:B","*! DMC time step !* DTDMC is the time s&
&tep for DMC runs in atomic units.")
call xyzzyaaad1("tpdmc","I:E","*! DMC number of factors in Pi-wts !* T&
&PDMC (T_p) is the number of time steps for which the effects of chang&
&es in the (theoretically constant) reference energy should be undone &
&in order to estimate the DMC energy at a given point. It is assumed t&
&hat the best estimates of the DMC energy separated by an amount great&
&er than this are not correlated by fluctuations in the reference ener&
&gy.  Thus T_p should exceed the timescale of fluctuations in the refe&
&rence energy. Umrigar suggests using T_p=10/tau where tau is the time&
& step.  If you set it to 9999 in the input, then the code will automa&
&tically use this value; if you set it to 0 then the reweighting schem&
&e for population control biasing will not be used. The latter is the &
&default.$$ Note that this procedure is not really all that useful.  I&
&f you suffer from population control biasing, then the scheme will he&
&lp to correct it, but you have to run for longer because the Pi-weigh&
&ts fluctuate, increasing the variance of the energy estimate. In fact&
&, it is not clear that you gain anything, because the extra run time &
&turns out to be about the same as you would need if you simply used a&
& sufficiently large population.  So use more configs to avoid populat&
&ion control biasing, basically.")
call xyzzyaaad1("vmc_decorr_period","I:I","*! VMC decorrelation period !&
&* Length of inner decorrelation loop in VMC. The code will do VMC_DEC&
&ORR_PERIOD configuration moves between successive evaluations of the &
&local energy (and other quantities to be averaged) in order to ensure&
& the configurations used are not significantly correlated. Setting VM&
&C_DECORR_PERIOD to a value greater than 1 should reduce serial correl&
&ation, but the length of the run will be increased.$$ It is normally &
&stated that typical values might be 3 or 4 for a pure VMC calculation&
&, and >= 10 during a config-generation run for optimization or DMC (t&
&hough this depends on the system, and clearly setting the decorrelati&
&on period to a higher value in DMC config generation is less importan&
&t than in wave function optimization, since the correlations will dis&
&appear as the DMC calculation evolves). The defaults are 3 for VMC ca&
&lculations and 15 for config generation runs. A slight complication i&
&s that if VMC_NSTEP is greater than VMC_NCONFIG_WRITE (as it might be&
& if you need more moves than the number of desired configs to calcula&
&te a VMC energy with a small enough error bar), then CASINO is able t&
&o exploit the extra moves to space the config writes further apart, a&
&nd it is no longer necessary for the inner decorrelation loop to be s&
&o long in config generation runs. In such a case, VMC_DECORR_PERIOD s&
&hould be taken to represent the *minimum* number of steps separating &
&config writes; internally the length of the decorrelation loop will b&
&e reduced as far as practical without going below the VMC default of &
&3.$$ Note that VMC_NSTEP refers to the number of moves at which energ&
&ies are evaluated or configurations are (sometimes) written out; this&
& is not affected by the value of VMC_DECORR_PERIOD. The value of VMC_&
&DECORR_PERIOD is ignored during equilibration.$$ Note that CASINO can&
& automatically determine VMC_DECORR_PERIOD to maximize run efficiency&
&. This is done by adding an extra set of moves after equilibration to&
& compute an estimation of the correlation time of the local energies.&
& The feature is enabled by setting VMC_DECORR_PERIOD=0.")
call xyzzyaaad1("vmc_nconfig_write","I:B","*! Number of configs to wri&
&te in VMC !* Total number of configurations to be written out in VMC &
&for later use (wave-function optimization or DMC). This number must b&
&e <= VMC_NSTEP (though you may want to set VMC_NSTEP to be significan&
&tly greater than VMC_NCONFIG_WRITE to get an acceptable error bar on &
&the energy; this is useful for e.g. judging the success of an optimiz&
&ation after each stage). Since each processor always does the same nu&
&mber of steps, then VMC_NCONFIG_WRITE (and VMC_NSTEP) will be rounded&
& up to the nearest multiple of the number of processors (e.g. VMC_NCO&
&NFIG_WRITE=20 will be rounded up internally to 24 on 12 cores, and 24&
& configs will be written to config.out - 2 from each core). Note that&
& the config.out file will still be written even if VMC_NCONFIG_WRITE &
&is zero, since this file is used to store the current state of the sy&
&stem at the end of every VMC block (equivalent to writing one config,&
& though of course multiple cores write multiple configs to save the s&
&tate). Writing of config.out may be suppressed completely with an app&
&ropriate value for the CHECKPOINT keyword, and the data will be held &
&in memory between different stages of the calculation.")
call xyzzyaaad1("con_loc","T:B","*! Configs directory !* Directory in &
&which config files are kept to be written/read. Default: './'.")
call xyzzyaaad1("dmc_target_weight","D:B","*! Total target weight in D&
&MC !* Total target weight in DMC, summed over all processors. This is&
& synonymous with the ""target population"" of configs, except that DM&
&C_TARGET_WEIGHT is allowed to be non-integer. Typically DMC_TARGET_WE&
&IGHT will be the same as VMC_NCONFIG_WRITE when RUNTYPE=vmc_dmc, thou&
&gh it does not have to be.$$ It may seem bizarre to allow non-integer&
& total target weights, but a possible use for this is in increasing t&
&he parallel efficiency when you have a very small population per core&
&.$$ Suppose your target weight is 1 config per core and you are runni&
&ng on 100000 cores. Half the time your total population will be a bit&
& higher than 100000, and when this happens nearly all of your cores w&
&ill spend half their time twiddling their thumbs waiting for the smal&
&l number of cores that have two configurations to finish the iteratio&
&n.  So around 25% of the computer time is wasted.$$ If instead you se&
&t your target weight to 0.98 configs per core then it is very unlikel&
&y that any cores will have two configurations.  Instead you have on a&
&verage 2% of your cores sitting idle, which is sad, but still more ef&
&ficient than having large numbers of cores wait for a small number of&
& over-burdened cores.")
call xyzzyaaad1("growth_estimator","L:I","*! Calculate DMC growth esti&
&mator !* If this flag is set to T then the growth estimator of the DM&
&C energy is evaluated in addition to the usual mixed estimator.  A st&
&atistically significant difference between the mixed estimator and th&
&e growth estimator for the energy normally implies the presence of ti&
&me-step bias.  Other than that, the growth estimator is not generally&
& useful because it has a significantly greater statistical error than&
& the mixed estimator.")
call xyzzyaaad1("gautol","D:I","*! Gaussian tolerance !* GAUTOL contro&
&ls the accuracy of evaluating orbitals expanded in Gaussian basis set&
&s.  Roughly speaking, any contribution from a Gaussian-type function &
&to an orbital at a point in space is neglected if exp(-a * r^2) is le&
&ss than 10^-GAUTOL.  Typical values in solids are around 6-7.  This p&
&arameter can have a large effect on the cost of a calculation in peri&
&odic systems - you can wind it up much higher in molecules.")
call xyzzyaaad1("dbarrc","I:I","*! DBAR recalculation period !* DBARRC&
& is the number of moves between full recalculation of the cofactor (=&
&'DBAR') matrices.  Basically every time an electron move is accepted &
&in equilibration/VMC/DMC, the update_dbar routine is called which upd&
&ates these matrices using the efficient Eq. (26) of Fahy et al., PRB &
&42, 3503 (1990).  As a numerical precaution that the (unstable!) upda&
&te procedure is working, every DBARRC accepted moves the DBAR matrice&
&s and determinant are recomputed from scratch from the orbitals in th&
&e Slater matrix. If the new DBAR differs by too much from the old upd&
&ated DBAR, then the program ought to be stopped (but in fact it is no&
&t, since it happens at least once per simulation and is very irritati&
&ng). In principle one can boost the value of DBARRC up to a fairly la&
&rge number before this happens, and this is a good idea since reevalu&
&ation of the matrix costs quite a lot.  The default value is 100,000.&
& Tests show that for systems with 1024 particles per spin channel it &
&may be safe to do up to 1,000,000 updates with accuracy better than s&
&ingle precision.$$ NOTE: the DBARRC keyword had a different meaning v&
&ery early in the life of CASINO where it was appropriate to set it to&
& a value of e.g. 10. Re-using an old input file with such a setting n&
&ow can cause the code to slow down by 1-2 orders of magnitude without&
& the user necessarily understanding why (real world examples have bee&
&n observed of people doing this). It is therefore now forbidden to se&
&t DBARRC to a value lower than the default; if the user has a genuine&
& reason for wishing to do this he/she may search for the error trap i&
&n the source code and in the runqmc script and comment it out.")
call xyzzyaaad1("interaction","T:I","*! Interaction type !* INTERACTIO&
&N determines the type of interaction between the simulated particles.&
&  INTERACTION can take the following values:$ 'none'       : non-inte&
&racting particles;$ 'coulomb'    : Coulomb interaction;$ 'ewald'     &
& : periodic Coulomb interaction computed using Ewald summation;$ 'mpc&
&'        : periodic Coulomb interaction computed using the model peri&
&odic Coulomb (MPC) method;$ 'ewald_mpc'  : compute and report both Ew&
&ald and MPC results, but use Ewald in DMC propagation;$ 'mpc_ewald'  &
&: compute and report both Ewald and MPC results, but use MPC in DMC p&
&ropagation;$ 'manual'     : compute a user-defined interaction (see t&
&he  MANUAL_INTERACTION block input keyword);$ 'ewaldpp'$ 'ewaldpp_mpc&
&'$ 'mpc_ewaldpp': as their above counterparts, but using an electron-&
&electron pseudopotential for the Ewald interaction, whose parameters &
&(see Eq. (4) of https://doi.org/10.1103/PhysRevB.92.075106) must be s&
&pecified in the MANUAL_INTERACTION block.$$ The values 'coulomb' and &
&'ewald' can be used interchangeably, although 'coulomb' should strict&
&ly refer to aperiodic systems and 'ewald' to periodic systems.$$ The &
&MPC interaction is generally significantly faster than the Ewald inte&
&raction and should give smaller finite size effects.  Note that curre&
&ntly the MPC interaction is not implemented for 1D systems.")
call xyzzyaaad1("manual_interaction","B:I","*! Manual particle interac&
&tions !* This block defines the interaction form and parameters when &
&INTERACTION = 'manual'.  The first line specifies the interaction typ&
&e, and subsequent lines define the parameters, to be given as 'name :&
& value' (real and int) or 'name' (boolean).$$ The possible 'manual' i&
&nteractions and their corresponding parameters (real-valued [a.u.] an&
&d mandatory unless otherwise stated; m and q are taken from PARTICLES&
& block) are:$$ * 'square_well' $| V = height for r<width, 0 otherwise&
&$ | params: height (<0 for well), width$$ * 'poschl_teller' $| V = 2 &
&v_0 mu^2 / [m cosh^2(mu*r)]$ | params: mu, v_0 (<0 for well)$$ * 'har&
&d_sphere' $| V = Infinity for r<=D, lambda/r^3 otherwise$ | params: D&
& or R (=D/2), lambda (opt, default=0),$| op_spins (opt bool, restrict&
& to opp-spin)$$ * 'polynomial' $| V = sum_{k=0}^{order-1} c_k r^k for&
& r<cutoff, else 0$ | params: order (integer), cutoff, c_0..c_{order-1&
&} (opt, default=0)$$ * 'logarithmic' $| V = qi*qj * [log(2*rstar/r) -&
& Euler]$ | params: rstar (length scale)$$ * '2D_int' $| V = qi*qj * F&
&(r,rstar)$ | params: rstar (length scale)$ | note: F -> 1/r for r>>rs&
&tar and F -> log(2*rstar/r)-Euler for r<<rstar$$ * 'dipole'$ | V = d^&
&2/r^3$ | params: d^2$$ * 'pseudodipole'$ | (like polynomial, with ext&
&ra param d^2)$$ * '2D_tilted_dipole'$| (like dipole, with extra param&
& theta)$$ * '2D_tilt_pseudodipole'$| (like pseudodipole, with extra p&
&aram theta)")
call xyzzyaaad1("npcell","B:B","*! # primcells/axis !* NPCELL is a vec&
&tor of length three giving the number of primitive cells in each dime&
&nsion that make up the simulation cell. NB, for the polymer case, npc&
&ell(2) and npcell(3) must be 1; for the 2D slab case npcell(3) must b&
&e 1. If you wish to use a more general 3x3x3 supercell matrix to defi&
&ne the simulation cell, then you can define this with the SCELL_MATRI&
&X keyword.")
call xyzzyaaad1("scell_matrix","B:E","*! # Supercell matrix array !* S&
&CELL_MATRIX is a 3x3 integer matrix giving the supercell lattice vect&
&ors in terms of the primitive-cell lattice vectors.  This is a genera&
&lization of NPCELL.")
call xyzzyaaad1("input_example","E:B","*! Print example input !* If IN&
&PUT_EXAMPLE is T, then an example of a QMC input file will be written&
& out with all currently known keywords and their default values.  A m&
&odified version of this can be used as an input file in future runs."&
&)
call xyzzyaaad1("ibran","L:E","*! DMC branching !* IBRAN enables weigh&
&ting and branching in DMC.  Normally IBRAN = T for a DMC run; a value&
& of F is allowed for DMC checking purposes.")
call xyzzyaaad1("limdmc","I:E","*! Green function mods !* LIMDMC is us&
&ed to set the type of modifications to the Green function (see manual&
&).  Possible values are:$ 0 = no modifications applied;$ 1 = Depasqua&
&le limits to drift velocity and energy;$ 2 = Umrigar mods (need ALIMI&
&T);$ 3 = Rothstein-Vrbik mods (need ALIMIT).$ 4 = Umrigar mods to dri&
&ft velocity (need ALIMIT); Zen-Sorella-Alfe` mods to energy.$$ We str&
&ongly recommend you use LIMDMC=4.")
call xyzzyaaad1("alimit","D:E","*! DMC limit parameter !* ALIMIT is a &
&parameter required when LIMDMC is 2, 3 or 4.  A value of 0.25 was sug&
&gested by Umrigar et al. for all-electron calculations, but a value o&
&f 0.5 may be more appropriate for pseudopotential calculations.  The &
&answers are essentially insensitive to the precise value.  ALIMIT is &
&not required if NUCLEUS_GF_MODS is set to true.")
call xyzzyaaad1("nucleus_gf_mods","L:E","*! Green fn mods for bare nuc&
&lei !* NUCLEUS_GF_MODS is the switch for enabling the use of the modi&
&fications to the DMC Green's function in the presence of bare nuclei,&
& suggested in J. Chem Phys. 99, 2865 (1993).")
call xyzzyaaad1("mpc_cutoff","P:E","*! G vector cutoff for MPC !* MPC_&
&CUTOFF is the energy cutoff for G vectors used in (a) the FFT of the &
&MPC interaction, and (b) the FFT of the one-particle density required&
& when generating the mpc.data file. We use the convention that we inc&
&lude the set of G vectors such that (1/2)|G|^2<mpc_cutoff. The progra&
&m will suggest a value for MPC_CUTOFF if the existing value is unsuit&
&able, or if the user inputs a value of zero. This is a keyword of typ&
&e 'Physical' hence you need to supply units, such as 'ev', 'ry', 'har&
&tree', 'kcal/mol' etc. The default is 30.d0 hartree.")
call xyzzyaaad1("non_local_grid","I:I","*! Non-local integration rule !&
&* NON_LOCAL_GRID selects the grid for non-local integration, ranging &
&from coarse (low NON_LOCAL_GRID value) to fine (high NON_LOCAL_GRID v&
&alue) to finer grids.  The value is assumed to be the same for all at&
&oms if it is controlled through this keyword; you can provide an over&
&ride value of NON_LOCAL_GRID for particular atoms at the top of the c&
&orresponding pseudopotential file, where it is called NLRULE1.$$ The &
&following table gives the grid details:$$ +--------------------------&
&-------------------------------+$ | NON_LOCAL_GRID    Exactly integra&
&tes l=...   No. points |$ +------------------------------------------&
&---------------+$ |       1                     0                    &
&  1    |$ |       2                     2                      4    |&
&$ |       3                     3                      6    |$ |     &
&  4                     5                     12    |$ |       5     &
&                5                     18    |$ |       6             &
&        7                     26    |$ |       7                    1&
&1                     50    |$ +-------------------------------------&
&--------------------+$$ Notice that NON_LOCAL_GRID=5 offers no theore&
&tical advantage over NON_LOCAL_GRID=4, and is significantly more expe&
&nsive (+50% points). We recommend that NON_LOCAL_GRID=5 not be used.$&
&$ The default value is NON_LOCAL_GRID=4 (this is also adopted if NON_&
&LOCAL_GRID is given a negative value).")
call xyzzyaaad1("opt_info","I:E","*! Optimizer information level !* Co&
&ntrols amount of information displayed and/or written out during opti&
&mization.$ In variance minimization:$ 1 = display no information;$ 2 &
&= display variance and energies at each iteration;$ 3 = also print pa&
&rameters, derivatives and intermediate evaluations;$ 4 = calculate an&
&d print weights as well;$ 5 = write out configs and their energies et&
&c as they are read in (lots of data).$ In energy minimization:$ 1 = v&
&ery little information, no extra output files;$ 2 = basic information&
&, no extra output files;$ 3 = full information, write matrix algebra &
&log file;$ 4 = also write full B~ and BH~ matrix files;$ 5 = also wri&
&te SVD components files (if SVD used).")
call xyzzyaaad1("vm_reweight","L:E","*! Varmin reweighting !* In a var&
&iance minimization calculation, if VM_REWEIGHT is F then all weights &
&are set to unity.  This is the default.  What does this mean? When us&
&ing weights in correlated sampling, the procedure may become numerica&
&lly unstable, particularly for large system sizes. The characteristic&
& of these instabilities is that during the minimization procedure a f&
&ew configurations (often only one) acquire a very large weight.  The &
&estimate of the variance is then reduced almost to zero by a set of p&
&arameters which are found to give extremely poor results in a subsequ&
&ent QMC calculation.  One can overcome this instability by using more&
& configurations. Various alternative ways of dealing with this instab&
&ility have been devised.  One method is to limit the upper value of t&
&he weights or to set the weights equal to unity.  In our calculations&
& for large systems we normally set the weights to unity, which is ach&
&ieved by setting VM_REWEIGHT to F.  It can be verified that the effec&
&t of this approximation is in general negligible by regenerating conf&
&igurations and carrying out a new variance minimization process.  On &
&the other hand, if the initial trial wave function is poor, then some&
& of the configuration generation/variance minimization cycles can be &
&bypassed by using the weights.  There is also some evidence that the &
&optimization of 'difficult' parameters works better when VM_REWEIGHT &
&is T.")
call xyzzyaaad1("opt_fixnl","L:E","*! Fix nonlocal energies in optimiz&
&ation !* Fix non-local contribution to local energy in optimization. &
&In VARMIN, this gives a large speedup for systems where pseudopotenti&
&als are used, and yields good results in most cases. In EMIN the spee&
&dup is the same, but tends to destabilize the optimization. Defaults &
&are T for VARMIN, and F for EMIN. For both methods, when this is F, t&
&he non-local integration grids are fixed for each electron in each co&
&nfiguration to ensure consistent, smooth variance/energy surfaces.")
call xyzzyaaad1("opt_maxeval","I:I","*! Varmin max evaluations !* Maxi&
&mum number of evaluations during optimization.")
call xyzzyaaad1("vm_forgiving","L:E","*! Varmin no whinge !* If VM_FOR&
&GIVING is set to T [default], CASINO won't consider it an error if th&
&e config energies from VMC do not agree with the initial energies in &
&variance minimization.")
call xyzzyaaad1("opt_jastrow","L:I","*! Optimize Jastrow factor !* Dur&
&ing optimization, allow the parameters in the Jastrow factor to be op&
&timized.")
call xyzzyaaad1("opt_det_coeff","L:I","*! Optimize determinant coeffs !&
&* During optimization, allow the determinant coefficients to be optim&
&ized.")
call xyzzyaaad1("opt_geminal","L:I","*! Optimize geminal coeffs !* Dur&
&ing optimization, allow the optimization of the geminal coefficient  &
&matrix.")
call xyzzyaaad1("postfit_vmc","L:B","*! Perform post-fit VMC !* If POS&
&TFIT_VMC is set to T then an extra VMC calculation will be performed &
&with the final optimized wave function when RUNTYPE=vmc_opt or opt_vm&
&c, to enable one to see the effect of the final optimization on the e&
&nergy etc.. This is done by default. Unless POSTFIT_KEEP_CFG is set t&
&o T, this final VMC run will not generate any configurations.")
call xyzzyaaad1("postfit_keep_cfg","L:B","*! Keep postfit VMC configur&
&ations !* If POSTFIT_KEEP_CFG is set to T then the configurations gen&
&erated in  the post-fit VMC calculation will be written to config.out&
& (the default is F.")
call xyzzyaaad1("lcutofftol","D:E","*! Local PP cutoff tol !* LCUTOFFT&
&OL is used to define the cutoff radius for the local part of the pseu&
&dopotential.  It is the maximum deviation of the local potential from&
& -z/r at the local cutoff radius.  The default of 1.d-5 is normally a&
&dequate.")
call xyzzyaaad1("nlcutofftol","D:E","*! Nonlocal PP cutoff tol !* NLCU&
&TOFFTOL is used to define the cutoff radius for the non-local parts o&
&f the pseudopotential.  It is defined as the maximum deviation of the&
& non_local potentials from the local_potential at the non-local cutof&
&f radius.  The default of 1.d-5 is normally adequate.")
call xyzzyaaad1("orbbuf","L:I","*! DMC orbital buffering !* Setting OR&
&BBUF=T turns on DMC orbital buffering.  This is an efficiency device &
&in which buffered copies of orbitals/gradients/Laplacians are kept fo&
&r later reuse.  This has a high memory cost.  Orbital buffering shoul&
&d always be used unless you start running out of memory, hence the ab&
&ility to turn it off.")
call xyzzyaaad1("blip_periodicity","I:I","*! Periodicity with blip bas&
&is !* Orbitals expanded in a blip basis can be periodic in zero, one,&
& two or three dimensions.  BLIP_PERIODICITY specifies the number of d&
&imensions in which the orbitals are periodic (-1,0,1,2,3).  Note that&
& if BLIP_PERIODICITY is 1 then the system is assumed to be periodic i&
&n the direction of lattice vector 1, while if BLIP_PERIODICITY is 2 t&
&hen the system is periodic in the directions of lattice vectors 1 and&
& 2.  If BLIP_PERIODICITY is -1 (the default) then the periodicity is &
&deduced from the value of the PERIODIC keyword (F-->0, T-->3).  In al&
&l cases, the simulation cell is the parallelepiped defined by the lat&
&tice vectors placed at the origin so in non-3D cases make sure you pu&
&t your atoms somewhere near the middle of it. Note that k points may &
&only be used in periodic directions.")
call xyzzyaaad1("expot","L:B","*! Use external potential !* If EXPOT i&
&s set to T, then an external potential is read from the file 'expot.d&
&ata' and included as a summed contribution to the total energy.")
call xyzzyaaad1("magnetic_field","L:B","*! Use external magnetic field !&
&* If MAGNETIC_FIELD is set to T, then an external magnetic field is r&
&ead from the file 'expot.data'.  The magnetic field contributes to th&
&e kinetic energy.")
call xyzzyaaad1("opt_orbitals","L:I","*! Optimize orbital parameters !&
&* During optimization, allow the parameters in the orbitals to be opt&
&imized.")
call xyzzyaaad1("structure_factor","L:E","*! Accumulate structure fact&
&or !* If STRUCTURE_FACTOR is set to true then the structure factor wi&
&ll be accumulated in the expval.data file (periodic systems only).")
call xyzzyaaad1("cerefdmc","D:E","*! DMC EREF update const !* Constant&
& used in updating the reference energy EREF during initial DMC diffus&
&ion to ground state using the latest population control algorithm. A &
&value of 1.d0 is usually appropriate.")
call xyzzyaaad1("makemovie","L:I","*! Make a movie !* Plot the electro&
&n positions every movieplot moves.")
call xyzzyaaad1("movieplot","I:I","*! Frame length for movie !* Plot t&
&he electron positions every movieplot moves.")
call xyzzyaaad1("movienode","I:I","*! Movie processor !* Plot the elec&
&tron positions on processor: movienode.")
call xyzzyaaad1("moviecells","L:I","*! Plot n.n. supercells in movie !&
&* If false makemovie will plot the unit cell, if true, the n.n. cells&
& in the xy plane will also be written.")
call xyzzyaaad1("edist_by_ion","B:I","*! Initial electron distribution !&
&* The optional EDIST_BY_ION block allows fine control of the initial &
&distribution of the electrons before equilibration starts.  The stand&
&ard algorithm shares out the electrons amongst the various ions weigh&
&ted by the pseudo-charge/atomic number of the ion.  Each electron is &
&placed randomly on the surface of a sphere surrounding its parent ion&
&.  There are certain situations, for example a simple crystal with a &
&very large lattice constant, where the standard algorithm in the POIN&
&TS routine may give a bad initial distribution which cannot be undone&
& by equilibrating for a reasonable amount of time.  This keyword allo&
&ws a user-defined set of electron/ion associations to be supplied.  T&
&he syntax is to supply N_ion lines within the block which look like e&
&.g. 1 4 4, where the three numbers are: the ion sequence number; the &
&number of up-spin electrons associated with this ion ; the number of &
&down-spin electrons associated with this ion. Alternatively one may u&
&se the EDIST_BY_IONTYPE keyword block where you replace the ion seque&
&nce number with the ion type  sequence number and the information is &
&supplied only for each particular  type of ion. This generally saves &
&typing.")
call xyzzyaaad1("edist_by_iontype","B:I","*! Initial electron distribu&
&tion !* The optional EDIST_BY_IONTYPE block allows fine control of th&
&e initial distribution of the electrons before equilibration starts. &
& The standard algorithm shares out the electrons amongst the various &
&ions weighted by the pseudo-charge/atomic number of the ion.  Each el&
&ectron is placed randomly on the surface of a sphere surrounding its &
&parent ion.  There are certain situations, for example a simple cryst&
&al with a very large lattice constant, where the standard algorithm i&
&n the POINTS routine may give a bad initial distribution which cannot&
& be undone by equilibrating for a reasonable amount of time.  This ke&
&yword allows a user-defined set of electron/ion associations to be su&
&pplied.  The syntax is to supply N_ion lines within the block which l&
&ook like e.g. 1 4 4, where the three numbers are: the ion type sequen&
&ce number; the number of up-spin electrons associated with each atom &
&of this type ; the number of down-spin electrons associated with each&
& atom of this type.  Alternatively one may use the EDIST_BY_ION keywo&
&rd block where you replace the sequence number of each type of atom w&
&ith the sequence number of each ion.  This will obviously be useful i&
&n magnetic systems, and in putting zigzag stripes of oxygen holes in &
&those manganites and high-Tc superconductors which CASINO is so good &
&at.")
call xyzzyaaad1("ewald_control","D:E","*! Ewald accuracy control !* EW&
&ALD_CONTROL is the percentage increase (from the default) of the cuto&
&ff radius for the reciprocal space sum in the Ewald interaction - use&
&d for calculating electrostatic interactions between particles in per&
&iodic systems.  Its default value is zero.  Increasing this will caus&
&e more vectors to be included in the sum, the effect of which is to i&
&ncrease the range of the Ewald gamma parameter over which the energy &
&is constant (the default gamma should lie somewhere in the middle of &
&this range).  This need only be done in exceptional circumstances and&
& the default should be fine for the general user.")
call xyzzyaaad1("orb_norm","D:I","*! Orbital normalization !* Allows u&
&ser to change normalization of orbitals by multiplying all of them by&
& this constant.  Of course this should have no effect on the energy b&
&ut it can be useful if the Slater determinant starts going singular, &
&as it might for example for some very dilute/low density systems.")
call xyzzyaaad1("printgscreening","L:E","*! Print Gaussian screening i&
&nfo !* Before doing a periodic Gaussian calculation, CASINO prepares &
&lists of potentially significant (primitive) cells and sites in each &
&such cell which could contain Gaussians having a non-zero value in a &
&reference primitive cell centred on the origin.  Zero is defined as 1&
&0^-GAUTOL. Turning on the PRINTGSCREENING flag prints out the importa&
&nt information about this screening - it is turned off by default.")
call xyzzyaaad1("qmc_plot","B:I","*! Plot orbitals etc. in line/plane/&
&volume !* This block allows you to plot the value of certain quantiti&
&es along a segment A-B / parallelogram with sides A-B, A-C / parallel&
&epiped with edges A-B, A-C, A-D.  The data will be plotted in a forma&
&t suitable for xmgr/grace in the file 'lineplot.dat' or in a format s&
&uitable for gnuplot in '2Dplot.dat' or '3Dplot.dat'.$$ In order to pr&
&oduce the plot, RUNTYPE must be set to 'plot'.$$ The block has the fo&
&llowing format: $ LINE 1: what_to_plot (orb, orb_gradx, orb_grady, or&
&b_gradz, orb_lap, wfn, nodes, energy, eipot, expot);$ LINE 2: dimensi&
&onality ndim (1,2,3);$ LINE 3: no of points along the ndim directions&
&;$ LINES 4-: xyz coords of point A; of point B; of point C (if reqd.)&
&;$ of point D (if reqd.).$$ Three additional lines have to be added f&
&or orbital plots:$ LINE A1: number of orbitals (norb);$ LINE A2: norb&
& integers identifying the orbitals to be plotted;$ LINE A3: norb inte&
&gers identifying the spin/species for each of the orbitals.$$ For wav&
&e-function and local-energy plots, the coordinate defined by LINES 4-&
& refer to electron 1. The coordinates of all remaining electrons are &
&taken from a configuration obtained by VMC equilibration (without Jas&
&trow factor). Alternatively the positions of several electrons may be&
& fixed by the following lines:$ LINE B1: number of fixed electrons (>&
&= 0);$ LINE B2 and onwards (if the number of fixed electrons is nonze&
&ro): spin, number and (x, y, z) coordinates for each of the fixed ele&
&ctrons.$$ For electron-ion potential, external potential and node plo&
&ts, no lines have to be added.")
call xyzzyaaad1("plot","B:E","*! Plot quantities in line/plane/volume !&
&* This block allows you to plot the value of certain quantities along&
& a segment A-B / parallelogram with sides A-B, A-C / parallelepiped w&
&ith edges A-B, A-C, A-D.  The data will be plotted in a format suitab&
&le for post-processing by a CASINO utility (yet to be written), to a &
&file named '1Dplot.dat', '2Dplot.dat' or '3Dplot.dat'.  This facility&
& will replace QMC_PLOT eventually.$$ In order to produce the plot, RU&
&NTYPE must be set to 'plot'.$$ The block has the following format:$$ &
&<what-to-plot>$ electron <ie> spin <ispin>$ <dimensionality>$ grid <n&
&grid>$ A <A-coords>$ B <B-coords>$ C <C-coords>$ D <D-coords>$ fix el&
&ectron <ie_fix1> spin <ispin_fix1> @ <fix1-coords>$$ where:$ - <what-&
&to-plot> indicates what to plot, e.g. orb_grad, wfn, energy, etc. The&
& list of available plot subjects is dynamic, so you should set this t&
&o 'help' and run CASINO, which will display the list for the current &
&system;$ - <ie> and <ispin> are the particle and particle-type indice&
&s of the particle that is moved to generate the plot.$ - <dimensional&
&ity> is the dimensionality of the plot.$ - <ngrid> is a set of <dimen&
&sionality> integers defining the density of the plot grid$ - <A-coord&
&s>, <B-coords>, <C-coords> and <D-coords> are the coordinates of A, B&
&, C and D (only specify those required according to <dimensionality>)&
&$ - <ie_fix1> and <ispin_fix1> are the particle and particle-type ind&
&ices of a particle that ought to be fixed at <fix1-coords>; any numbe&
&r of particles can be fixed with additional 'fix' lines.")
call xyzzyaaad1("expval_cutoff","P:E","*! G vector cutoff for expval !&
&* EXPVAL_CUTOFF is the energy cutoff for G vectors used in the evalua&
&tion of expectation values accumulated in reciprocal space (e.g. the &
&density, spin density, pair-correlation function etc.). We use the co&
&nvention that we include the set of G vectors such that (1/2)|G|^2<mp&
&c_cutoff (or |G|^2<mpc_cutoff before September 2013). The value of EX&
&PVAL_CUTOFF is ignored if an expval.data file is already present, in &
&which case the G vector set(s) given therein are used instead.  If yo&
&u set it to zero, then the program will suggest a value. This is a ke&
&yword of type 'Physical' hence you need to supply units, such as 'ev'&
&, 'ry', 'hartree', 'kcal/mol' etc. The default is 30.0 hartree.")
call xyzzyaaad1("lwdmc","L:E","*! Enable weighted DMC !* Enable weight&
&ed DMC, where each configuration carries a weight that is simply mult&
&iplied by the branching factor after each move; only if the weight of&
& a configuration goes outside certain bounds (above WDMCMAX or below &
&WDMCMIN) is it allowed to branch or be combined with another configur&
&ation. This should reduce excessive population fluctuations, which is&
& generally held to be a good thing. Note that setting LWDMC=T means t&
&hat your population will generally fluctuate around a value other tha&
&n DMC_TARGET_WEIGHT (after an initial transient);  the chances of bei&
&ng killed if your weight is below 1 or duplicated if your weight is a&
&bove 1 depend on the values of WDMCMIN and WDMCMAX, and in general th&
&is is not symmetrical.")
call xyzzyaaad1("dmc_equil_fixpop","D:E","*! Fix pop and total weight &
&during initial DMC equilibration !* If the VMC and DMC energies are v&
&ery different, the population increases during the initial phase of e&
&quilibration before the reference energy can counteract. This paramet&
&er (between 0.0 and 1.0) specifies an initial fraction of the equilib&
&ration phase during which the population and total weight are fixed t&
&o the target weight. Setting this parameter to e.g. 0.5 will prevent &
&such explosions and should have a negligible effect on the equilibrat&
&ion time.")
call xyzzyaaad1("lwdmc_fixpop","L:E","*! Fix population in lwdmc !* Th&
&is flag activates the LWDMC variant with fixed population. By interpr&
&eting WDMCMIN and WDMCMAX relative towards the current population the&
& population and the total weight are decoupled. The population is nea&
&rly fixed while the total weight fluctuates as usual. While this gene&
&rally reduces the statistical efficency of the DMC algorithm, it is a&
& simple way to eliminate population explosions or extinction in cases&
& of small population and large population fluctuation. WARNING: this &
&is *not* a solution for walkers trapped in singular points of the wav&
&e functions, nor is it a solution for populations that get trapped in&
& high-energy states. Be careful about this option when you do not kno&
&w the reason for the population problems in the first place.")
call xyzzyaaad1("dmc_norm_conserve","L:E","*! Enable norm-conserving D&
&MC !* Use the norm-conserving DMC algorithm.  This eliminates populat&
&ion fluctuations.  Experimental algorithm: use with caution.")
call xyzzyaaad1("dmc_poprenorm","L:E","*! Enable DMC pop-renormalizati&
&on !* Control the DMC configuration population by randomly deleting o&
&r copying configurations after branching with the reference energy se&
&t equal to the best estimate of the ground-state energy. This can be &
&used to maintain a constant population of configs per core, provided &
&the value of DMC_TARGET_WEIGHT is an integer multiple of the number o&
&f cores. Note that non-integer values of DMC_TARGET_WEIGHT are not al&
&lowed when using DMC_POPRENORM. Note also that DMC_POPRENORM is not i&
&n general recommended because of the population control errors it can&
& theoretically introduce, though in general these are likely to be sm&
&all.")
call xyzzyaaad1("wdmcmin","D:E","*! Minimum weight !* IF LWDMC=T then &
&WDMCMIN is the minimum weight.  Now type 'casinohelp lwdmc'.")
call xyzzyaaad1("wdmcmax","D:E","*! Maximum weight !* IF LWDMC=T then &
&WDMCMAX is the maximum weight.  Now type  'casinohelp lwdmc.'")
call xyzzyaaad1("checkwfn","L:E","*! Numcheck analytic wfn derivs !* E&
&nable a numerical check of the analytic orbital derivatives coded in &
&the various routines such as gauss_per/gauss_mol/bwfdet/pwfdet etc.")
call xyzzyaaad1("vmc_ave_period","I:I","*! Energy-averaging period in &
&VMC !* Number of consecutive local energies that are averaged togethe&
&r in VMC before writing them to the vmc.hist file. The only effect of&
& this keyword is to reduce the size of vmc.hist: the number of lines &
&written in a VMC calculation is VMC_NSTEP/VMC_AVE_PERIOD.")
call xyzzyaaad1("pair_corr_sph","L:E","*! Accumulate real-space PCF !*&
& If PAIR_CORR_SPH is set to true then the spherically-averaged real-s&
&pace pair-correlation function will be accumulated in the expval.data&
& file (via a process of 'binning' the electron-electron separations).&
& This currently works for periodic homogeneous systems and finite iso&
&tropic systems such as electron-hole biexcitons. For periodic systems&
& with atoms you can use the PAIR_CORR keyword instead which gives you&
& the full (non-spherically averaged) pair-correlation function accumu&
&lated in reciprocal space.")
call xyzzyaaad1("pcfs_nbins","I:E","*! Number of bins for real-space P&
&CF !* Number of bins to be used when accumulating the real-space pair&
&-correlation function.  Enter a negative value or omit the keyword to&
& use the default value.  If an expval.data file is present then the v&
&alue given in expval.data will be used and the input keyword will be &
&ignored.")
call xyzzyaaad1("pcfs_rcutoff","P:E","*! Radius of binned region for r&
&eal-space PCF !* Radius of region to be considered when accumulating &
&the pair-correlation function.  The default value in periodic systems&
& (the Wigner-Seitz cell radius) is generally appropriate; however, fo&
&r finite systems such as excitonic complexes, PCFS_RCUTOFF should be &
&set to something rather larger than the size of the complex.  Note th&
&at units (e.g., bohr) should be supplied after the value of PCFS_RCUT&
&OFF.  Enter a negative value or omit the keyword to use the default v&
&alue.  If an expval.data file is present then the value given in expv&
&al.data will be used and the input keyword will be ignored.")
call xyzzyaaad1("kwarn","L:I","*! Disable KE check in PW basis !* If t&
&he kwarn flag is set to T, then the routine PWFDET_SETUP will issue a&
& warning whenever the kinetic energy calculated from the supplied orb&
&itals differs from the DFT kinetic energy given in the pwfn.data file&
& by more than an internal tolerance (usually set to 10^-6). If the fl&
&ag is F, then CASINO will stop with an error message on detecting thi&
&s condition. Note that in cases where the DFT calculation which gener&
&ated the orbitals used fractional occupation numbers, the kinetic ene&
&rgy mismatch is very likely to occur since QMC deals in principle onl&
&y with integer occupation numbers, hence the existence of this flag."&
&)
call xyzzyaaad1("writeout_vmc_hist","L:I","*! Write vmc.hist file in V&
&MC !* If WRITEOUT_VMC_HIST is set to T then the energy components, et&
&c., are written to the file vmc.hist during a VMC simulation.  This w&
&ill not occur if WRITEOUT_VMC_HIST is set to F.  Furthermore, the con&
&fig.out file required to continue a VMC calculation will only be prod&
&uced if WRITEOUT_VMC_HIST is T.")
call xyzzyaaad1("writeout_dmc_hist","L:I","*! Write dmc.hist file in D&
&MC !* If WRITEOUT_DMC_HIST is set to T then the energy components, et&
&c., are written to the file dmc.hist during a DMC simulation.  This w&
&ill not occur if WRITEOUT_DMC_HIST is set to F.")
call xyzzyaaad1("vmc_method","I:I","*! VMC method !* VMC_METHOD select&
&s which version of VMC to use:$ 1 = Electron-by-electron algorithm, e&
&valuating configuration energies at the end of the configuration move&
&; or$ 3 = Configuration-by-configuration algorithm, evaluating config&
&uration energies before and after the move and adding a weighted sum &
&of these to the accumulation arrays.$$ Method 1 is the default and is&
& recommended.$$ There used to be a Method 2 but after extensive testi&
&ng we concluded that it did not offer any advantage over the other me&
&thods and was hard to support.")
call xyzzyaaad1("vmc_ionjump","D:E","*! Probability for trying jumps b&
&etween ions in VMC !* In the special case of nearly separated molecul&
&e fragments (e.g., for intermolecular forces) electrons may get trapp&
&ed in the energetically less favourable fragment for a long time duri&
&ng a VMC run due to the large distance between the fragments. Setting&
& vmc_ionjump to a small, nonzero probability will cause the VMC routi&
&ne to try a long-distance jump from one ion to another once in a whil&
&e, improving the sampling of disjoint areas of the configuration spac&
&e. Detailed balance is preserved.")
call xyzzyaaad1("jasbuf","L:I","*! Jastrow buffering !* If JASBUF is s&
&et to T then the one-body terms in the Jastrow factor for each electr&
&on in each configuration is buffered: saves time at the expense of me&
&mory. Clearly this will have no effect in systems without one-body te&
&rms in the Jastrow (these are the Chi and Q terms at present).")
call xyzzyaaad1("neighprint","I:I","*! Neighbour analysis !* NEIGHPRIN&
&T = n will generate a printout of the first n stars of neighbours of &
&each atom in the primitive cell, with the relevant interatomic distan&
&ces given in both Angstrom and a.u. If n=0 or if you are an atom-free&
& electron or electron-hole fluid phase, then the keyword has no effec&
&t. Note that activating cusp corrections when using a Gaussian basis &
&(the default) will trigger a neighbour analysis irrespective of the v&
&alue of this keyword.")
call xyzzyaaad1("ranluxlevel","I:E","*! Quality/cost of random numbers !&
&* To generate parallel streams of pseudo-random numbers for its stoch&
&astic algorithms, CASINO uses an implementation of RANLUX. This is an&
& advanced pseudo-random number generator based on the RCARRY algorith&
&m proposed in 1991 by Marsaglia and Zaman. RCARRY used a subtract-and&
&-borrow algorithm with a period on the order of 10^171 but still had &
&detectable correlations between numbers. Martin Luescher proposed the&
& RANLUX algorithm in 1993; RANLUX generates pseudo-random numbers usi&
&ng RCARRY but throws away numbers to destroy correlations. RANLUX tra&
&des execution speed for quality through the choice of a 'luxury level&
&' given in CASINO by the RANLUXLEVEL input keyword. By choosing a lar&
&ger luxury setting one gets better random numbers slower. By the test&
&s available at the time it was proposed, RANLUX at its higher setting&
&s appears to give a significant advance in quality over previous gene&
&rators.$$ The luxury setting RANLUXLEVEL must be in the range 0-4.$ L&
&evel 0: equivalent to the original RCARRY of Marsaglia and Zaman, ver&
&y long period, but fails many tests.$ Level 1: considerable improveme&
&nt in quality over level 0, now passes the gap test, but still fails &
&spectral test.$ Level 2: passes all known tests, but theoretically st&
&ill defective.$ Level 3 [DEFAULT]: any theoretically possible correla&
&tions have very small chance of being observed.$ Level 4: highest pos&
&sible luxury, all 24 bits chaotic.")
call xyzzyaaad1("ranprint","I:E","*! Print random numbers !* Setting t&
&his keyword to a value greater than zero will cause the first RANPRIN&
&T numbers generated by the CASINO random number generator to be print&
&ed to a file 'random.log'. On parallel machines the numbers generated&
& on all processors are printed. The run script should pick out random&
&.log files from different stages of a calculation (e.g. VMC config ge&
&n/DMC equil/DMC stats accumulation) and rename them appropriately.")
call xyzzyaaad1("sparse","L:E","*! Activate sparse algorithms !* CASIN&
&O is capable of using sparse matrix algebra in some algorithms for ef&
&ficiency purposes. For systems which are definitely not sparse (orbit&
&als not well localized) then attempting to use sparse algorithms migh&
&t actually slow things down. Thus until we work out a better way you &
&can toggle this behaviour with the SPARSE flag. In these algorithms a&
& matrix element is considered to be zero if it less than the value of&
& the input keyword SPARSE_THRESHOLD.")
call xyzzyaaad1("sparse_threshold","D:E","*! Sparsity threshold !* CAS&
&INO sometimes uses sparse matrix algebra for efficiency purposes. Thi&
&s keyword defines a threshold such that a matrix element is taken to &
&be zero if it is less than this threshold. Changing this quantity mig&
&ht be used to trade speed for accuracy in some cases.")
call xyzzyaaad1("cusp_correction","L:I","*! Cusp-correct AE Gaussians &
&and STOs !* When expanded in a basis set of Gaussian functions, the e&
&lectron-nucleus cusp present in all-electron calculations is not repr&
&esented correctly, since the gradient of an atom-centred Gaussian is &
&necessarily zero there. Clearly this only matters for s-type GTFs, si&
&nce all functions of higher angular momentum vanish at the nucleus.  &
&When the CUSP_CORRECTION flag is activated, the s-type GTFs centred o&
&n each atom are replaced within a sphere of radius r_c by a function &
&which ensures that the electron-nucleus cusp condition is obeyed. Thi&
&s procedure can be expected to greatly reduce fluctuations in the loc&
&al energy in all-electron Gaussian calculations.$$ For STO wave funct&
&ions, the cusp condition can be exactly satisfied by a linear constra&
&int. This flag determines whether this constraint is enforced when re&
&ading in or optimizing a wave function.")
call xyzzyaaad1("dmc_equil_nstep","I:B","*! No of steps in DMC equil !&
&* Number of DMC steps performed on each processor in the DMC equilibr&
&ation phase, and consequently, the total number of local energy sampl&
&es (averaged over configs and processors) written to the dmc.hist fil&
&e. The equilibration phase may be partitioned into DMC_EQUIL_NBLOCK b&
&locks, but this does not affect the total number of steps (just how f&
&requently stuff is written out). However, if DMC_EQUIL_NSTEP is not d&
&ivisible by the number of blocks, then it will be rounded up to the n&
&earest multiple of DMC_EQUIL_NBLOCK. Furthermore, DMC_AVE_PERIOD cons&
&ecutive local energies may be averaged together in DMC before writing&
& them to the dmc.hist file (hence reducing its size), but again, if D&
&MC_EQUIL_NSTEP is not divisible by DMC_AVE_PERIOD, it will be rounded&
& up to the nearest multiple of it. Note the difference in parallel be&
&haviour compared to VMC_NSTEP, which is not a per processor quantity;&
& this is because the DMC phase is parallelized over configs.")
call xyzzyaaad1("dmc_equil_nblock","I:I","*! Number of blocks in DMC e&
&quilibration !* In cases when the BLOCK_TIME keyword is not used this&
& keyword defines the number of blocks into which the DMC equilibratio&
&n phase is divided (if DMC_EQUIL_NSTEP is not divisible by DMC_EQUIL_&
&NBLOCK, then the number of steps will be increased to the nearest mul&
&tiple of the number of blocks). Note that having multiple blocks does&
& not increase the amount of data collected, merely the frequency with&
& which data is written to files; the final answer should be essential&
&ly the same, irrespective of the number of blocks. Specifically, at t&
&he end of each equilibration block, the following significant actions&
& are performed:$$ (1) Write processor- and config-averaged data to dm&
&c.hist (one line for each step in the current block).$$ (2) Print mon&
&itoring data to the output file (block-averaged quantities).$$ (3) Ma&
&ke a backup copy of the config.out file (if catastrophe protection is&
& turned on with the DMC_TRIP_WEIGHT keyword).$$ (4) Write the dmc.sta&
&tus file.$$ (5) Write the current state of the system, and all config&
&s in the current population to the config.out file (note that by sett&
&ing the CHECKPOINT keyword to 0, this step can be skipped until the e&
&nd of the final block, or skipped completely if CHECKPOINT=-1, but th&
&is is not the default).$$ Note that if accumulating expectation value&
&s other than the energy, data is not written to the expval.data file &
&after each block, as it would be during the statistics accumulation p&
&hase.$$ Note that having too many blocks will make the code slower, a&
&nd if the run is not massively long it is perfectly in order to have &
&only one DMC block (which is the default).")
call xyzzyaaad1("dmc_stats_nstep","I:B","*! No of steps in DMC stats a&
&ccum !* If the value of STOP_METHOD='nstep', then this is the number &
&of DMC steps performed on each processor in the statistics accumulati&
&on phase, and consequently, the total number of local energy samples &
&(averaged over configs and processors) written to the dmc.hist file. &
&The accumulation phase may be partitioned into DMC_STATS_NBLOCK block&
&s, but this does not affect the total number of steps (just how frequ&
&ently stuff is written out). However, if DMC_STATS_NSTEP is not divis&
&ible by the number of blocks, then it will be rounded up to the neare&
&st multiple of DMC_STATS_NBLOCK. Furthermore, DMC_AVE_PERIOD consecut&
&ive local energies may be averaged together in DMC before writing the&
&m to the dmc.hist file (hence reducing its size), but again, if DMC_S&
&TATS_NSTEP is not divisible by DMC_AVE_PERIOD, it will be rounded up &
&to the nearest multiple of it. Note the difference in parallel behavi&
&our compared to VMC_NSTEP, which is not a per processor quantity; thi&
&s is because the DMC phase is parallelized over configs.$$ Note that &
&instead of executing a fixed number of DMC steps, CASINO may be reque&
&sted to do as many steps as are required to attain a target error bar&
& or to make the error bar as small as possible - see keyword STOP_MET&
&HOD.")
call xyzzyaaad1("dmc_stats_nblock","I:I","*! Number of blocks in DMC s&
&tatistics accumulation !* In cases when the BLOCK_TIME keyword is not&
& used this keyword defines the number of blocks into which the DMC st&
&atistics accumulation phase is divided (if DMC_STATS_NSTEP is not div&
&isible by DMC_STATS_NBLOCK, then the number of steps will be increase&
&d to the nearest multiple of the number of blocks). Note that having &
&multiple blocks does not increase the amount of data collected, merel&
&y the frequency with which data is written to files; the final answer&
& should be essentially the same, irrespective of the number of blocks&
&. Specifically, at the end of each accumulation block, the following &
&significant actions are performed:$$ (1) Write processor- and config-&
&averaged data to dmc.hist (one line for each step in the current bloc&
&k).$$ (2) Write processor- and config-averaged data to the expval.dat&
&a file (if accumulating expectation values other than the energy).$$ &
&(3) Print monitoring data to the output file (block-averaged quantiti&
&es).$$ (4) Make a backup copy of the config.out file (if catastrophe &
&protection is turned on with the DMC_TRIP_WEIGHT keyword).$$ (5) Make&
& a backup copy of the expval.data file (if it exists, and if catastro&
&phe protection is turned on with the DMC_TRIP_WEIGHT keyword).$$ (6) &
&Write the dmc.status file.$$ (7) Write the current state of the syste&
&m, and all configs in the current population to the config.out file (&
&note that by setting the CHECKPOINT keyword to 0, this step can be sk&
&ipped until the end of the final block, or skipped completely if CHEC&
&KPOINT=-1, but this is not the default).$$ Note that having too many &
&blocks will make the code slower, and if the run is not massively lon&
&g it is perfectly in order to have only one DMC block (which is the d&
&efault).")
call xyzzyaaad1("cusp_info","L:E","*! Print Gaussian cusp info !* If C&
&USP_CORRECTION is set to TRUE for an all-electron Gaussian basis set &
&calculation, then CASINO will alter the orbitals inside a small radiu&
&s around each nucleus in such a way that they obey the electron-nucle&
&ar cusp condition. If CUSP_INFO is set to true, then information abou&
&t precisely how this is being done will be printed to the output file&
&. Be aware that in large systems, this may produce a lot of output. F&
&urthermore, if you create a file called 'orbitals.in' containing an i&
&nteger triplet specifying which orbital/ion/spin you want, the code w&
&ill additionally print graphs of the specified orbital, radial gradie&
&nt, Laplacian and 'one-electron local energy' to the files orbitals.d&
&at, gradients.dat, laplacians.dat and local_energy.dat.  These graphs&
& may be viewed using xmgr/grace or similar plotting programs.")
call xyzzyaaad1("opt_maxiter","I:E","*! Max iterations in optimization !&
&* OPT_MAXITER specifies the largest permitted number of iterations of&
& the minimizer in both VARMIN and EMIN. Default number: 10.")
call xyzzyaaad1("ewald_check","L:E","*! Perform Ewald accuracy check !&
&* CASINO and the wave function generating program should be able to c&
&alculate the same value for the nuclear repulsion energy, given the s&
&ame crystal structure. By default CASINO therefore computes the Ewald&
& interaction and compares it with the value given in the wave functio&
&n file. If they differ by more than 10^-5, then CASINO will stop and &
&complain. If you have a justifiable reason for doing so, you may turn&
& off this check by setting EWALD_CHECK to F.")
call xyzzyaaad1("permit_den_symm","L:I","*! Symmetrize QMC charge data !&
&* If this flag is set to T then the symmetry of the SCF charge densit&
&y (in mpc.data) will be imposed upon the QMC charge-density data for &
&use in the MPC interaction (with the justification that imposing an e&
&xact condition on the charge density can''t hurt) and also when writi&
&ng the QMC density to expval.data. It is possible however that DMC wi&
&ll break the symmetry of the SCF calculation; in this case the user s&
&hould turn off PERMIT_DEN_SYMM.")
call xyzzyaaad1("qmc_density_mpc","L:I","*! Use QMC density in MPC int !&
&* If this flag is set to T then the QMC charge density data at the en&
&d of the expval.data file will be used to compute the MPC interaction&
&, rather than the default SCF density in the mpc.data file. This is l&
&ikely to be useful in cases such as the Wigner crystal where the Hart&
&ree-Fock charge density is very different to the true charge density &
&(it is too localized) as opposed to, say, the Fermi fluid where the H&
&artree-Fock charge density is exact. Note that using this option is l&
&ikely to increase the time taken to evaluate the MPC interaction; in &
&both DFT and QMC cases, the code counts backwards from the end of the&
& list of G vectors and discards all those before the first non-zero o&
&ne (where zero is defined by some threshold like 1.d-6). In the HF/DF&
&T case this tends to give a large reduction in the size of the vector&
& to be evaluated. However, the random noise in the QMC density coeffi&
&cients is likely to exceed the zero threshold for all G, and the vect&
&or will likely be untruncated. It is important therefore to use a val&
&ue for EXPVAL_CUTOFF which is not too large when using this facility,&
& in order that the total number of G vectors in the expansion is not &
&too large.")
call xyzzyaaad1("esupercell","L:I","*! Energy/per supercell in output !&
&* By default total energies and their components are printed as energ&
&ies per primitive cell.  Switching this flag to T forces printing of &
&energies per simulation cell in the output file.")
call xyzzyaaad1("dmc_trip_weight","D:I","*! DMC catastrophe threshold !&
&* In the course of a DMC simulation, it is possible for a configurati&
&on ""population explosion"" to occur. If DMC_TRIP_WEIGHT is set to 0 &
&then nothing will be done about this. If DMC_TRIP_WEIGHT>0 then it wi&
&ll attempt to restart with a different random number sequence from th&
&e beginning of the previous block if the iteration weight exceeds DMC&
&_TRIP_WEIGHT. A general suggestion for its value would be 2-3 times D&
&MC_TARGET_WEIGHT (but see the discussion in the manual about this).")
call xyzzyaaad1("max_rec_attempts","I:I","*! Maximum number DMC recove&
&ries !* This is the maximum number of times DMC will attempt to resta&
&rt a block if it continues to encounter catastrophes. Relevant only i&
&f the DMC_TRIP_WEIGHT keyword is set to a nonzero value.")
call xyzzyaaad1("molgscreening","L:I","*! Screen molecular Gaussians !&
&* Toggle on and off the use of screening in Gaussian basis set calcul&
&ations of molecules i.e. the division of space into boxes and the pre&
&paration of lists of which Gaussian basis functions have a significan&
&t weight in each box. Should speed up the calculation of large molecu&
&les. The screening information can take up a reasonable amount of mem&
&ory, hence this keyword.")
call xyzzyaaad1("bsmooth","L:E","*! Smooth truncate loc orbs !* If bsm&
&ooth is set to true then localized orbitals are interpolated smoothly&
& to zero beyond their cutoff radius.  Otherwise, they are truncated a&
&bruptly.  The latter is the default as it generally produces better r&
&esults.")
call xyzzyaaad1("relativistic","L:E","*! Relativistic correction !* If&
& RELATIVISTIC is T, then calculate relativistic corrections to the en&
&ergy using perturbation theory.  Note that for the moment this can on&
&ly be done for closed-shell systems.")
call xyzzyaaad1("isotope_mass","D:E","*! Nuclear mass override !* This&
& keyword can be used to define the nuclear mass (in amu) if one wishe&
&s to override the default value (which is averaged over isotopes acco&
&rding to their abundances).  The default is used if ISOTOPE_MASS is s&
&et to zero. The atomic mass unit (amu) in this sense means 'the ratio&
& of the average mass per atom of the element to 1/12 of the mass of 1&
&2C.")
call xyzzyaaad1("vm_w_min","D:E","*! VM minimum weight !* Minimum valu&
&e that a configuration weight may take during weighted variance minim&
&ization. This parameter should have a value between zero and one. Not&
&e that the limiting is not applied if VM_W_MAX = 0.")
call xyzzyaaad1("vm_w_max","D:E","*! VM maximum weight !* Maximum valu&
&e that a configuration weight may take during weighted variance minim&
&ization. Set this to zero if you do not wish to limit the weights; ot&
&herwise it should be greater than 1.")
call xyzzyaaad1("opt_dtvmc","I:B","*! Optimize VMC time step !* IF OPT&
&_DTVMC is set to 1 the VMC time step (initially the value of keyword &
&DTVMC) is optimized so that the VMC acceptance ratio is roughly 50%. &
&This is the default. To prevent this, set OPT_DTVMC to 0. The value o&
&f OPT_DTVMC is ignored if the input block DTVMCS is supplied.$$ CASIN&
&O can also maximize the diffusion constant with respect to DTVMC. Thi&
&s can be enabled by setting OPT_DTVMC=2. In a first stage, DTVMC is v&
&aried to get an acceptance ratio of 50%, so as to have decent statist&
&ics to perform the diffusion-constant maximization stage. This is onl&
&y useful for VMC_METHOD=3, where it is the default.")
call xyzzyaaad1("pcf_rfix","B:E","*! Fixed particle type and position &
&in PCF calc. !* This block contains two lines. The first line gives t&
&he type of particle to be fixed during accumulation of the pair corre&
&lation function g(r,r') (e.g. 1-4 typically up/down spin electron, up&
&/down spin hole); the second line gives the coordinates of the positi&
&on r at which to fix it (in au). This applies to the reciprocal-space&
& PCF activated with the PAIR_CORR input keyword. It also applies *in &
&principle* to the spherical real space PCF activated with the PAIR_CO&
&RR_SPH input keyword, in the sense that the format of expval.data all&
&ows it, but the accumulation of the spherical PCF with fixed particle&
&s has not yet been implemented.")
call xyzzyaaad1("on_top_pair","B:E","*! Particles to place on top of e&
&ach other in RPMD calc. !* This block contains two lines consisting o&
&f the particle type and index (integers) of each of two particles to &
&be forced to stay on top of each other throughout a VMC calculation. &
&This is intended for the evaluation of recombining-pair momentum dens&
&ities and 'one-body' density matrices in hole-in-HEG systems -- wrong&
& values will be reported for other expectation values, including the &
&energy. Note that the time step of the first-specified particle appli&
&es to the pair. This block must not be used in DMC or optimization ru&
&ns.")
call xyzzyaaad1("spin_density","L:E","*! Accumulate spin densities !* &
&Setting SPIN_DENSITY to T will activate the accumulation of separate &
&up- and down-spin densities in the expval.data file.")
call xyzzyaaad1("pair_corr","L:E","*! Accumulate rec. space PCF !* Set&
& PAIR_CORR to T to accumulate the reciprocal-space pair-correlation f&
&unction in the expval.data file. Currently restricted to periodic sys&
&tems. Note you also need to give the position and type of a fixed par&
&ticle using the PCF_RFIX block. Note that if the density is homogeneo&
&us, or if you want only the spherically-averaged pair-correlation fun&
&ction (for a very restricted class of systems), you should use the PA&
&IR_CORR_SPH keyword.")
call xyzzyaaad1("opt_cycles","I:B","*! Number of optimization cycles !&
&* Number of cycles of configuration generation and optimization runs &
&to be carried out if RUNTYPE=vmc_opt or opt_vmc. For variance minimiz&
&ation, 3--6 cycles is typical; for energy minimization, 5--10 cycles &
&is usual unless only determinant coefficients are being optimized, in&
& which case 1--2 cycles will be enough.")
call xyzzyaaad1("loc_tensor","L:E","*! Accumulate localization tensor !&
&* If LOC_TENSOR is set to true then the localization tensor will be a&
&ccumulated in the expval.data file (periodic systems only).")
call xyzzyaaad1("jastrow_plot","B:I","*! Plot components of Jastrow fa&
&ctor !* This block allows you to plot the u(rij), w(rij), chi(ri), f(&
&ri,rj,rij), p(rij) and q(ri) terms in the Jastrow factor.$$ Either TH&
&REE or SIX lines should be given in the input block:$ LINE 1: flag fo&
&r whether the Jastrow factor is to be plotted (0=NO, 1=YES);$ LINE 2:&
& spin of particle i (=1,2,..);$ LINE 3: spin of particle j (=1,2,..);&
&$ The following three lines are optional:$ LINE 4: (x,y,z)-position o&
&f particle j (used in plots of f and p);$ LINE 5: vector with the dir&
&ection in which i is moved (used in plots of f, p and q);$ LINE 6: po&
&sition vector of a point on the straight line along which electron i &
&moves (used in plots of f, p and q).$$ Note that the nucleus is assum&
&ed to lie at the origin in plots of f. The jastrow_value_f_?.dat file&
&s contain the value of f against the distance from the point given in&
& line 6.  Likewise for plots of p and q.$$ For u and chi we just plot&
& u(r) and chi(r) and the derivatives thereof.")
call xyzzyaaad1("popstats","L:I","*! Collect population statistics !* &
&To measure the statistical efficiency of a DMC calculation, the varia&
&nce of the local energy needs to be obtained. This measurement is nee&
&ded to estimate the effect of the exponentially scaling statistical i&
&nefficiency for large systems. See PRB 81, 035119 (2010).")
call xyzzyaaad1("timing_info","L:I","*! Activate subroutine timers !* &
&Setting TIMING_INFO to F (the default) will turn off the collection o&
&f subroutine timings.  You might want to do this as the the timing ro&
&utines can adversely affect system performance on certain computers (&
&such as Alpha or PC clusters) especially for small systems.")
call xyzzyaaad1("cusp_control","D:E","*! All-electron cusp radius cont&
&rol !* This undocumented parameter is required in the *OLD* procedure&
& (which is in fact much better than the NEW one, and is the default) &
&used to make all-electron orbitals expanded in Gaussian functions sat&
&isfy the electron-nuclear cusp conditions (activated by the CUSP_CORR&
&ECTION input keyword).  To activate the old algorithm for doing this,&
& set OLD_CUSP_RADIUS=.true. in the routine cusp_setup in gaussians.f9&
&0. The radius inside which the form of the orbitals is modified is de&
&termined by looking at fluctuations in the local energy. This radius &
&('rcusp') is set to the largest distance from the nucleus at which th&
&e deviation from the 'ideal' curve has a magnitude of greater then (z&
&ion^2/CUSP_CONTROL), where zion is the nuclear charge.  The default v&
&alue of CUSP_CONTROL is 50.0.  Note that this keyword will have no ef&
&fect on hydrogen atoms, which are treated as a special case.")
call xyzzyaaad1("vm_use_E_guess","L:E","*! Use guess of GS energy !* S&
&etting the flag to true will cause the 'variance' in VARMIN to be eva&
&luated using VM_E_GUESS in place of the average energy of the config &
&set, in an attempt to combine energy minimization with VARMIN.  Other&
&wise the least-squares function will simply be the variance of the co&
&nfiguration local energies.")
call xyzzyaaad1("vm_E_guess","P:E","*! Estimate of GS energy !* If VM_&
&USE_E_GUESS is true then VM_E_GUESS should be supplied as an estimate&
& of the ground-state energy. This is a keyword of type 'Physical' hen&
&ce you need to supply units, such as 'ev', 'ry', 'hartree', 'kcal/mol&
&' etc.")
call xyzzyaaad1("e_offset","P:I","*! Energy offset !* The E_OFFSET key&
&word gives a constant shift in the total energy per electron such tha&
&t the final results will be E = E_calc - E_offset. The default is zer&
&o. This is a keyword of type 'Physical' hence you need to supply unit&
&s, such as 'ev', 'ry', 'hartree', 'kcal/mol' etc.")
call xyzzyaaad1("vm_smooth_limits","L:I","*! Smooth Jastrow cutoff lim&
&its !* By setting this keyword to T (default), the optimizing routine&
& used in variance minimization is sent a smoothed version of the set &
&of parameters.  This only affects those which are to remain bounded s&
&uch as Jastrow cutoffs.  The result is a set of parameters which can &
&vary in the range (-Inf,+Inf), which can be more convenient than igno&
&ring out-of-range values without the minimizer knowing.  A suitable h&
&yperbolic function is used for mapping 'limited' values into 'extende&
&d' ones and vice versa.")
call xyzzyaaad1("backflow","L:B","*! Use backflow corrections !* Turns&
& on backflow corrections. Backflow parameters are read from correlati&
&on.data and, if optimized (OPT_BACKFLOW = T), written to correlation.&
&out(.x) .")
call xyzzyaaad1("opt_backflow","L:B","*! Optimize backflow parameters !&
&* During optimization, allow the backflow parameters to be optimized.&
&")
call xyzzyaaad1("use_jastrow","L:B","*! Use a Jastrow function !* Use &
&a wave function of the Slater-Jastrow form, where the Jastrow factor &
&exp(J) is an optimizable object that multiplies the determinant part &
&in order to introduce correlation in the system.  The Jastrow factor &
&must be provided in the correlation.data file --see the CASINO manual&
& for information on the format.")
call xyzzyaaad1("plot_backflow","B:I","*! Generate plot of backflow tr&
&ansformation !* This block allows plotting the backflow transformatio&
&n after VMC equilibration. The block should contain 2 lines, plus an &
&optional line for plotting the Phi term:$ LINE 1: 0 or 1 to (de-)acti&
&vate this facility;$ LINE 2: kspin, knumber and zposition;$ LINE 3: v&
&alue of fixed electron-ion distance riI.$$ This will produce various &
&files:$ * bfconfig.dat (reference config),$ * bfconfigx.dat (associat&
&ed quasi-particle config),$ * bfions.dat (coordinates of nuclei for w&
&hich backflow terms exist),$ * bfeta_<s>.dat (eta vs rij for each spi&
&n-pair type <s>),$ * bfmu_<s>_<set>.dat (mu vs ri for each spin type &
&<s> in each set <set>),$ * bfphi.dat (contribution of Phi to the 3D b&
&ackflow displacement on electron j vs 2D projection of rjI on the pla&
&ne defined by electron j, ion I in set <set>, and electron i at dista&
&nce riI from the nucleus with spin such that <s> is the spin-pair typ&
&e of i and j),$ * bffield.dat (3D backflow displacement on electron (&
&kspin,knumber) vs its 2D position on the plane z=zposition).")
call xyzzyaaad1("vm_linjas_method","T:I","*! Optimization method for a&
&cc varmin !* VM_LINJAS_METHOD specifies the method used to minimize t&
&he quartic LSF.  VM_LINJAS_METHOD should be one of:$ 'CG' (conjugate &
&gradients),$ 'MC' (Monte Carlo),$ 'LM' (line minimization),$ 'SD' (st&
&eepest descents),$ 'BFGS' (Broyden-Fletcher-Goldfarb-Shanno),$ 'BFGS_&
&MC' (BFGS and Monte Carlo),$ 'CG_MC' (conjugate gradients and Monte C&
&arlo),$ 'GN' (Gauss-Newton) or$ 'GN_MC' (Gauss-Newton and Monte Carlo&
&).")
call xyzzyaaad1("vm_linjas_its","I:E","*! Max iterations in acc varmin !&
&* VM_LINJAS_ITS specifies the maximum number of conjugate-gradients, &
&steepest-descent or BFGS iterations to be performed if VM_LINJAS_ITS &
&is 'CG', 'SD', 'BFGS', 'CG_MC' or 'BFGS_MC'.  If VM_LINJAS_ITS is 'MC&
&', 'LM', 'CG_MC' or 'BFGS_MC' then it specifies the number of line mi&
&nimizations to be performed.")
call xyzzyaaad1("splot","L:E","*! Line plotter s components !* Tell li&
&ne plotter to do s component of orbitals rather than full plot. Usefu&
&l for analysing Gaussian cusp corrections.")
call xyzzyaaad1("cusp_threshold","D:E","*! Zero orbital threshold !* I&
&f the magnitude of the s component of a Gaussian orbital is less than&
& this threshold, then it will not be cusp corrected.")
call xyzzyaaad1("vm_filter","L:I","*! Filter outlying configs !* This &
&keyword activates filtering of configurations in VARMIN by making the&
& weights (artificially) energy-dependent, i.e., W_i = W(|E_i-E_ave|).&
&  This method uses two parameters: VM_FILTER_THRES and VM_FILTER_WIDT&
&H.")
call xyzzyaaad1("vm_filter_thres","D:I","*! Filter threshold !* When l&
&imiting outlying configs in VARMIN (by setting the VM_FILTER flag to &
&T), the maximum deviation from the average energy at which the (artif&
&icial) weight of a configuration W_i =  W(|E_i-E_ave|) is kept equal &
&to unity is VM_FILTER_THRES times the square root of the unreweighted&
& variance.  After such limit, the weight is brought to zero using a g&
&aussian of width VM_FILTER_WIDTH times the square root of the unrewei&
&ghted variance.  By default, VM_FILTER_THRES is 4.0 and VM_FILTER_WID&
&TH is 2.0.")
call xyzzyaaad1("vm_filter_width","D:I","*! Gaussian filter width !* W&
&hen limiting outlying configs in VARMIN (by setting the VM_FILTER fla&
&g to T), the maximum deviation from the average energy at which the (&
&artificial) weight of a configuration W_i =  W(|E_i-E_ave|) is kept e&
&qual to unity is VM_FILTER_THRES times the square root of the unrewei&
&ghted variance.  After such limit, the weight is brought to zero usin&
&g a gaussian of width VM_FILTER_WIDTH times the square root of the un&
&reweighted variance.  By default, VM_FILTER_THRES is 4.0 and VM_FILTE&
&R_WIDTH is 2.0.")
call xyzzyaaad1("max_cpu_time","P:I","*! Maximum CPU time !* If the CP&
&U time elapsed since the start of a QMC simulation exceeds MAX_CPU_TI&
&ME and a suitable point in the algorithm is reached, then CASINO will&
& halt gracefully. This should make it easier to carry out e.g. multip&
&le DMC runs on a computer with a queueing system, particularly when u&
&sed with the --continue or --auto-continue runqmc options.$$ The way &
&this works is that at the end of each block of moves, CASINO will che&
&ck whether doing one more block will exceed the time limit. If so, it&
& will perform an emergency stop, writing to the output file any chang&
&es to the input file that must be made in order to restart the job (i&
&n a form readable both by humans and by runqmc). The user must theref&
&ore define the block length appropriately -- most usefully via the BL&
&OCK_TIME keyword -- such that the time taken per block is a sufficien&
&tly small fraction of MAX_CPU_TIME. $$Note that in DMC, if CHECKPOINT&
&=0 in input and there is a job time limit, it is strongly recommended&
& that MAX_CPU_TIME is used to ensure the config.out file is written o&
&ut if the required CPU time is longer than the time limit. MAX_CPU_TI&
&ME is a physical parameter with dimensions of time, and the units mus&
&t be specified as e.g. 1 day, 24 hr, 1440 min, or 86400 s. See also t&
&he MAX_REAL_TIME keyword.")
call xyzzyaaad1("max_real_time","P:I","*! Maximum real time !* If the &
&wall-clock time elapsed since the start of a QMC simulation exceeds M&
&AX_REAL_TIME and a suitable point in the algorithm is reached, then C&
&ASINO will halt gracefully.  This should make it easier to carry out &
&e.g. multiple DMC runs on a computer with a queueing system, particul&
&arly when used with the --continue or --auto-continue runqmc options.&
&$$ The way this works is that at the end of each block of moves, CASI&
&NO will check whether doing one more block will exceed the time limit&
&. If so, it will perform an emergency stop, writing to the output fil&
&e any changes to the input file that must be made in order to restart&
& the job (in a form readable both by humans and by runqmc). The user &
&must therefore define the block length appropriately -- most usefully&
& via the BLOCK_TIME keyword -- such that the time taken per block is &
&a sufficiently small fraction of MAX_REAL_TIME. $$ Note that in DMC, &
&if CHECKPOINT=0 in input and there is a job time limit, it is strongl&
&y recommended that MAX_REAL_TIME is used to ensure the config.out fil&
&e is written out if the required time is longer than the time limit. &
&MAX_REAL_TIME is a physical parameter with dimensions of time, and th&
&e units must be specified as e.g. 1 day, 24 hr, 1440 min, or 86400 s.&
& See also the MAX_CPU_TIME keyword.")
call xyzzyaaad1("use_orbmods","L:B","*! Single particle orbital modifi&
&cations !* If USE_ORBMODS is set to T then the orbital-modification b&
&lock in correlation.data will be read, and the modifications will be &
&applied to the single particle orbitals orbitals.  This only applies &
&to ATOM_BASIS_TYPEs 'numerical', 'gaussian' or 'slater-type.")
call xyzzyaaad1("particles","B:E","*! Define custom quantum particles !&
&* Using the PARTICLES block the user can define quantum particles (ot&
&her than electrons, which can be introduced using NEU/NED) to be used&
& in the QMC calculation.  The format of each line is '<i> <charge/|e|&
&> <mass/m_e> <spin/hbar> <name>'.  A negative value of the mass indic&
&ates that the following three lines give an (anisotropic) 3x3 mass te&
&nsor.  CASINO should decide whether each particle type is a fermion o&
&r a boson (based on the spin), and select the appropriate way to comb&
&ine the one-particle orbitals (symmetric combination [not currently i&
&mplemented] or antisymmetric Slater determinants).  The particles def&
&ined here can be assigned orbitals using the FREE_PARTICLES block.")
call xyzzyaaad1("free_particles","B:E","*! Parameters for free particl&
&es !* This block sets the parameters defining the behaviour of the or&
&bitals that are not atom-related in the system.  The geometry of the &
&system can be given using$ 'r_s <rs>',$ 'dimensionality <d>',$ 'cell_&
&geometry' (followed by d lines with d reals corresponding to the unsc&
&aled cell vectors), $ 'z-separation <z>' (assign particles to the top&
& layer (2D) or wire (1D) using 'top-layer <p1> [<p2> [...]]'), $ 'heg&
&_nlayers <n>'     (define number n of layers),$ 'heg_ylayer  <l> <y>'&
& (define y-coordinate of ""layer"" l),$ 'heg_zlayer  <l> <z>' (define&
& z-coordinate of layer l), and $ 'heg_layer   <p> <l>' (assign partic&
&le type p to layer l).$ These parameters are only required if ATOM_BA&
&SIS_TYPE is 'none'.$$ The number and type of the orbitals can be give&
&n using lines with the syntax 'particle <i> : <n> orbitals <orb> [orb&
&-options]' (if all determinants contain the same orbital type), or 'p&
&article <i> det <det> : <n> orbitals <orb> [orb-options]', where <det&
&> is the term in the multideterminant expansion, <i> must be 1, 2 or &
&a number given in the PARTICLES block (1/2 are up/down electrons), <n&
&> is the number of free particles/orbitals belonging in the <det>-th &
&determinant, and '<orb> [orb-options]' is one of the following: 'free&
&', 'crystal sublattice <s>', 'biex1', 'biex2', 'biex3' or 'pairing <j&
&>', <j> being the particle type which <i> is paired with.  Should the&
& orbitals have optimizable parameters, these must be provided in corr&
&elation.data.  Wigner crystal geometry is specified using the keyword&
&s 'crystal_type <type> <n> sublattice[s] [repeat <r>]' (type = 'cubic&
&', 'fcc', 'bcc', 'rectangular', 'hexagonal' or 'triangular', which mu&
&st match 'dimensionality' and 'cell_geometry', or 'manual'), and 'sub&
&lattice <s> [antiferro[magnetic]] offset <x y z>' for predefined latt&
&ices, and 'sublattice <s> manual <n> site[s]' followed by <n> lines o&
&f the form <x y z> defining the sites for manual lattices.")
call xyzzyaaad1("fixed_particles","B:E","*! Create a set of fixed, cha&
&rged particles !* When setting up a model system, one can place fixed&
&, charged particles within the simulation cell by using the fixed_par&
&ticles block. This can be used to study, e.g., electron-hole complexe&
&s in the presence of fixed donor and acceptor ions.  The block consis&
&ts of one line for each fixed particle, where the lines are of the fo&
&rm '<charge> <x> <y> <z>', where x, y and z are the Cartesian compone&
&nts of the fixed charge's position.  The charge must be an integer.")
call xyzzyaaad1("primitive_cell","B:E","*! Override prim cell latvec !&
&* Sometimes the 'primitive lattice vectors' in the xwfn.data file do &
&not correspond to the true primitive cell.  If these are needed (e.g.&
& for accumulation of the density) then override values can be supplie&
&d here.")
call xyzzyaaad1("plot_expval","B:I","*! Plot exp values in line/plane/&
&volume !* The utility 'plot_expval' allows you to plot expectation va&
&lues calculated by CASINO and stored in the file 'expval.data'.  It t&
&akes its instructions from this block in input which tells it about t&
&he geometrical region over which the data will be plotted.  Where the&
& geometry is clear, this input block is not required (e.g. spherical &
&PCF/structure factor).  The geometrical region may be a line AB / pla&
&ne AB-AC / volume AB-AC-AD. The data will be plotted in a format suit&
&able for xmgr/grace in the file 'lineplot.dat' or in a format suitabl&
&e for gnuplot in '2Dplot.dat' or '3Dplot.dat'.  These latter two file&
&s can be quickly visualized with the 'plot_2D' utility.$$ The block, &
&which is ignored by CASINO itself, has the following format:$ LINE 1:&
& dimensionality of plot ndim 1/2/3, OR EQUIVALENTLY, line/plane/volum&
&e;$ LINE 2: No. of points along each of the ndim directions;$ LINES 3&
&-: xyz coords of point A ; of point B ; of point C (if reqd.) ; of po&
&int D (if reqd.).")
call xyzzyaaad1("ke_verbose","L:E","*! Detailed info KE tests !* CASIN&
&O performs numerical tests to determine whether the kinetic energies &
&computed during the run will be correct.  Such tests are carried out &
&after VMC equilibration, and will only produce concise output about t&
&he outcome.  However, if the flag KE_VERBOSE is set to T, CASINO will&
& print out information throughout the process.  The default is F.  Se&
&e also KE_FORGIVE.")
call xyzzyaaad1("ke_forgive","L:E","*! Allow failing KE tests !* CASIN&
&O performs numerical tests to determine whether the kinetic energies &
&computed during the run will be correct.  If KE_FORGIVE is set to F, &
&CASINO will regard this as an error and stop.  The default is T.  Not&
&e that although the procedure is rather stable, there may be cases in&
& which poor numerics causes narrow fails.  See also KE_VERBOSE.")
call xyzzyaaad1("dipole_moment","L:B","*! Accumulate elec. dipole mome&
&nt !* If this flag is set to T then CASINO will accumulate the expect&
&ation value of the electric dipole moment p. It will also evaluate th&
&e expectation of p^2. The data p_x, p_y, p_z and p^2 are written to (&
&v/d)mc.hist like energy components, rather than into expval.data, and&
& their value and error bars are determined by reblocking. This can on&
&ly be done for finite systems.$$ Note that the CASINO reblock utility&
& reports only the components and not the magnitude of the dipole mome&
&nt in order to allow the user to decide how to deal with the symmetry&
&. Suppose that symmetry dictates the dipole moment will point in the &
&x direction. The y and z components should be zero, but there will be&
& some noise when they are evaluated in QMC. If you work out p=sqrt(p_&
&x^2+p_y^2+p_z^2) then you will get something larger than p_x, tending&
& to p_x in the limit of perfect sampling (i.e. a biased estimate with&
& finite sampling). You will also get larger error bars on p than p_x.&
&")
call xyzzyaaad1("struc_factor_sph","L:E","*! Accumulate sph. struc. fa&
&ctor !* If STRUC_FACTOR_SPH is set to T then the spherically-averaged&
& structure factor will be accumulated in the expval.data file.  You s&
&hould also define a one-dimensional k point grid on which to calculat&
&e it using the EXPVAL_KGRID keyword block. Note this is implemented o&
&nly for homogeneous systems.")
call xyzzyaaad1("expval_kgrid","B:E","*! k point grid for exp values !&
&* This block contains a specification of one or more k point grids de&
&fined in 1, 2 or 3 dimensions (NOTE: this is more general than CASINO&
& actually requires at the moment). One-dimensional grid defined by li&
&ne AB, two dimensional grid by plane AB-AC, three-dimensional grid by&
& parallelepiped AB-AC-AD, all with an appropriate number of k points &
&along each direction. These grids may be used in the calculation of v&
&arious expectation values. Only used if appropriate expval keywords a&
&re set to T in input.$$ The block consists of the following lines:$ L&
&INE 1: No of k grids defined in this block;$ LINE 2: Which expectatio&
&n value uses this grid? (Currently 1=BLANK, 2=spherical structure fac&
&tor, 3=BLANK);$ LINE 3: dimensionality of current k grid NKDIM (1-3);&
&$ LINE 4: k point A coordinates - origin (au);$ LINE 5: k point B coo&
&rdinates (au), number of k along AB;$ LINE 6: [IF NKDIM=2 or 3] k poi&
&nt C coordinates (au), number of k along AC;$ LINE 7: [IF NKDIM=3] k &
&point D coordinates (au), number of k along AD ; Then repeat lines 2 &
&to 7 for each additional k grid.$$ Note that spherically averaged qua&
&ntities require only a one-dimensional k grid independently of the di&
&mensionality of the system, and that this is a radial coordinate so o&
&nly one number is required to specify each k point coordinate.")
call xyzzyaaad1("opt_method","T:I","*! Optimization method !* There ar&
&e currently four optimization methods implemented in CASINO:$ 'varmin&
&': variance minimization$ 'varmin_linjas': an accelerated variance mi&
&nimization technique for parameters that appear linearly in the Jastr&
&ow;$ 'emin': linear least-squares energy minimization$ 'madmin': mini&
&mization of the mean absolute deviation of the set of local energies &
&from the median$$ All these methods can be used to optimize all types&
& of parameters (Jastrow, orbitals, backflow and determinant coefficie&
&nts), except 'varmin_linjas' which can only be applied to optimize Ja&
&strow factors.$$ There are other keywords that affect the behaviour o&
&f each of these methods.  The default value of OPT_METHOD is 'varmin'&
&.")
call xyzzyaaad1("onep_density_mat","L:I","*! Accumulate 1p density mat&
&rix !* If ONEP_DENSITY_MAT is set to T, then the spherically averaged&
& one-particle density matrix will be computed.  This is only possible&
& if the system is homogeneous for the moment.")
call xyzzyaaad1("twop_density_mat","L:I","*! Accumulate 2p density mat&
&rix !* If TWOP_DENSITY_MAT is set to T, then the spherically averaged&
& two-particle density matrix will be computed.  This is only possible&
& if the system is homogeneous for the moment, and will increase the c&
&ost of the calculation significantly.")
call xyzzyaaad1("int_sf","L:E","*! Calc e-e int from strucfac !* If IN&
&T_SF is set to T, then the electron-electron interaction energy for a&
& periodic system will be calculated in terms of the structure factor.&
&  The structure factor should either have been accumulated in a previ&
&ous run and stored in an available expval.data file, or its accumulat&
&ion should be flagged for the current run.  Using this method the tot&
&al interaction energy can be separated into Hartree and exchange-corr&
&elation terms. [NB: this is a deprecated keyword, and should be repla&
&ced by HARTREE_XC.]")
call xyzzyaaad1("hartree_xc","L:E","*! Calc e-e int from strucfac/MPC !&
&* Flag the computation of separate Hartree and exchange-correlation (&
&XC) parts of the electron-electron interaction energy for a periodic &
&system. This may be done in two different ways, namely the structure &
&factor method and the MPC method. The computation thus requires eithe&
&r (1) structure factor information from a previously accumulated expv&
&al.data file or from setting STRUCTURE_FACTOR=T, or (2) the MPC inter&
&action to be active (through INTERACTION=mpc, mpc_ewald or ewald_mpc)&
&. If both these things are true then both methods will be used to com&
&pute the hartree/XC energies (the resulting numbers should agree reas&
&onably closely). If neither are true, then this keyword has no effect&
&. The default is T. Note that the MPC version only works with 3D peri&
&odicity.")
call xyzzyaaad1("future_walking","L:I","*! Enable future walking !* If&
& this flag is set to T, then future walking will be used to evaluate &
&pure estimators in DMC.  This is currently implemented only for sampl&
&ing the total energy.")
call xyzzyaaad1("cond_fraction","L:I","*! Accumulate cond fraction !* &
&If COND_FRACTION is set to T, then an improved estimator of the spher&
&ically averaged two-particle density matrix, from which one-body cont&
&ributions are subtracted, will be computed.  This is only available i&
&f the system is homogeneous, for the moment.")
call xyzzyaaad1("complex_wf","L:I","*! Use complex Slater wfn !* If CO&
&MPLEX_WF is set to T then CASINO will use a complex wave function and&
& the fixed-phase (rather than fixed-node) approximation will be appli&
&ed in DMC.  Using complex arithmetic is slower, and is only necessary&
& if the Hamiltonian is complex (e.g. if a magnetic field is present) &
&or if the boundary conditions on the wave function force it to be com&
&plex (e.g. for a periodic system with a general set of twist angles).&
&")
call xyzzyaaad1("virtual_node","I:E","*! Processor no in virtual paral&
&lel vm !* This parameter is not to be set manually.")
call xyzzyaaad1("virtual_nconfig","I:E","*! nconfigs in virtual parall&
&el vm !* This parameter is not to be set manually.")
call xyzzyaaad1("virtual_nnodes","I:E","*! nprocessors during virtual &
&parallel vm !* This parameter is not to be set manually.")
call xyzzyaaad1("use_tmove","L:I","*! Casula nl PP scheme in DMC !* If&
& USE_TMOVE is T then the Casula nonlocal pseudopotential scheme will &
&be used in DMC. So-called 'T-moves' will be performed in order to giv&
&e a DMC energy that is greater than or equal to the ground-state ener&
&gy. This  violates the detailed-balance principle at finite time step&
&s, but greatly improves the stability of the DMC algorithm when nonlo&
&cal pseudopotentials are used. The advantages of T-moves are that the&
&y restore the variational principle and help to prevent population ex&
&plosions; the disadvantages of T-moves are that the magnitude of the &
&error due to the locality approximation is generally larger, although&
& always positive, and the time-step bias is generally worse. [This la&
&tter problem is alleviated, to some extent, by using a symmetric bran&
&ching factor (Casula 2010) as opposed to the asymmetric one suggested&
& in his 2006 paper. This advice was implemented in CASINO in June 201&
&4.]. A further disadvantage is that this option requires a truly enor&
&mous amount of memory in systems with large numbers of particles (see&
&ing if this can be reduced remains a project). The default of USE_TMO&
&VE is F for the moment but unless memory issues are encountered we no&
&w actively recommend its use.")
call xyzzyaaad1("finite_size_corr","L:I","*! Eval. finite size correct&
&ions !* Calculate finite size corrections to kinetic energy and elect&
&ron-electron interaction energy using the Chiesa-Ceperley-Martin-Holz&
&mann / Drummond-Needs-Sorouri-Foulkes scheme.")
call xyzzyaaad1("random_seed","T:E","*! Random seed !* This keyword de&
&termines which random seed to use for the RANLUX random-number genera&
&tor. The default value of RANDOM_SEED is 'timer', which causes the sy&
&stem timer to be used as the seed. If RANDOM_SEED is set to 'standard&
&', the seed 314159265 is used. If the value of RANDOM_SEED is an inte&
&ger, that integer will be used as the random seed. The seed is printe&
&d to the output file so that calculations using RANDOM_SEED='timer' c&
&an be reproduced afterwards.$$ Note that, if RANDOM_SEED is an intege&
&r or 'standard' or 'timer' then, when restarting from a previous calc&
&ulation the value of RANDOM_SEED is ignored (except for any initial s&
&etup, such as evaluation of the twist-averaged Hartree-Fock energy of&
& a homogeneous electron gas), and the random-number sequence will gen&
&erally be continued from the saved state of the random-number generat&
&or stored in the config.in file. However, if RANDOM_SEED is 'timer_re&
&set' then the generator will be re-initialized from the system clock &
&after the config.in file is read. This might be useful, for example, &
&if a prior test has revealed that the standard sequence will lead to &
&a configuration giving rise to a population explosion.")
call xyzzyaaad1("custom_striplet_dep","B:E","*! Custom spin-triplet de&
&pendence !* This input block can be used to create spin-triplet group&
&ings (for the Jastrow H term only so far). CASINO does not currently &
&know how to generate non-trivial groupings automatically, and this bl&
&ock is the way to define them.  The format of the block is as in the &
&following example:$$ %block custom_striplet_dep$ no_striplet_deps 1$ &
&striplet_dep -1 4$ 1=2-3,1=2-4,1=1-3,1=1-4,2=2-3,2=2-4$ 1-3=4,2-3=4,1&
&-3=3,1-4=4,2-3=3,2-4=4$ 1=1=1,1=1=2,1=2=2,2=2=2$ 3=3=3,3=3=4,3=4=4,4=&
&4=4$ %endblock custom_striplet_dep$$ The equal signs denote that the &
&two particle types are equivalent within the group, whereas the hyphe&
&ns denote they are different. All spin-triplets must be included in a&
& group. See also CUSTOM_SPAIR_DEP and CUSTOM_SSINGLE_DEP.")
call xyzzyaaad1("custom_spair_dep","B:E","*! Custom spin-pair dependen&
&ce !* This input block can be used to create new spin-pair groupings &
&for the Jastrow factor, etc. For example, if one were studying a para&
&magnetic fluid bilayer, with spin-up and spin-down electrons in one p&
&lane (spins 1 and 2) and spin-up and spin-down electrons on the other&
& plane (spins 3 and 4) then one would want sets of u(r_ij) terms for &
&same-plane, same-spin pairs, same-plane, opposite-spin pairs and oppo&
&site-plane pairs. Here is an example:$ %block custom_spair_dep$ no_sp&
&air_deps 1    # Number of custom spin dependences$ spair_dep -1 3    &
& # Label (-1,-2,-3,...) and number of spin groups$ 1-1,2-2,3-3,4-4$ 1&
&-2,2-4$ 1-3,1-4,2-3,2-4$ %endblock custom_spair_dep$ All spin-pairs m&
&ust be included in a group. See also CUSTOM_SSINGLE_DEP.")
call xyzzyaaad1("custom_ssingle_dep","B:E","*! Custom spin-single depe&
&ndence !* This input block can be used to create new spin-single grou&
&pings for the Jastrow factor, etc.  Here is an example:$ %block custo&
&m_ssingle_dep$ no_ssingle_deps 1    # Number of custom spin dependenc&
&es$ ssingle_dep -1 2     # Label (-1,-2,-3,...) and number of spin gr&
&oups$ 1,2$ 3,4$ %endblock custom_ssingle_dep$ All spins must be inclu&
&ded in a group. See also CUSTOM_SPAIR_DEP.")
call xyzzyaaad1("use_gpcc","L:I","*! General-purpose cusp correction !&
&* If use_gpcc is set to T then short-ranged functions will be added t&
&o the orbitals to ensure that the Kato cusp conditions are satisfied.&
&")
call xyzzyaaad1("vmc_ntwist","I:I","*! Number of twist angles in VMC !&
&* Number of different 'twists' or offsets to the grid of k vectors to&
& be applied during a VMC twist-averaging run. Note that the usual key&
&words define the run length for a single twist angle, thus the run le&
&ngth is increased by a factor of VMC_NTWIST. Note also that if VMC_NT&
&WIST > 0, the values of VMC_NBLOCK and BLOCK_TIME are ignored. Settin&
&g this keyword to a value greater than zero requires the use of a com&
&plex wave function (COMPLEX_WF : T). $$ Note twist-averaging wholly w&
&ithin CASINO can currently be done only for electron(-hole) fluid pha&
&ses; for real systems with atoms one needs to couple with an external&
& code to regenerate the wave function after each twist. The various t&
&wistav_xxx scripts in CASINO/utils/twist can help with this - see the&
& manual.")
call xyzzyaaad1("vmc_reequil_nstep","I:I","*! Number of post-twist equ&
&il steps !* Total number of steps to take for a re-equilibration when&
& doing a twist-averaged VMC run. Currently, a re-equilibration only t&
&akes place when the twist angle is changed. Note this is a single-pro&
&cessor quantity; all cores run VMC_REEQUIL_NSTEP reequilibration step&
&s. Note also that the value of VMC_DECORR_PERIOD is ignored in reequi&
&librations. Electron(-hole) fluid phases only. For more details, see &
&the VMC_NTWIST keyword.")
call xyzzyaaad1("dmc_ntwist","I:I","*! Number of twist angles in DMC !&
&* Number of different 'twists' or offsets to the grid of k vectors to&
& be applied during DMC statistics accumulation. If DMC_NTWIST is 0 th&
&en the twist angle is not changed during DMC. After each change of tw&
&ist angle, the set of configurations needs to be re-equilibrated: hen&
&ce a value needs to be specified for DMC_REEQUIL_NSTEP. Setting DMC_N&
&TWIST>0 requires the use of a complex wave function (COMPLEX_WF : T).&
& Note that the usual keywords define the run length for a single twis&
&t angle, thus the run length is increased by a factor of DMC_NTWIST. &
&(Note twist-averaging wholly within CASINO can currently be done only&
& for electron(-hole) fluid phases; for real systems with atoms one ne&
&eds to couple with an external code to regenerate the wave function a&
&fter each twist. The various twistav_xxx scripts in CASINO/utils/twis&
&t can help with this - see the manual.)")
call xyzzyaaad1("dmc_reequil_nstep","I:I","*! No of post-twist equil m&
&oves !* Number of steps to take for a re-equilibration when doing a t&
&wist-averaged DMC run. Currently, a re-equilibration only takes place&
& when the twist angle is changed. Electron(-hole) fluid phases only. &
&For more details, see the DMC_NTWIST keyword.")
call xyzzyaaad1("dmc_reequil_nblock","I:I","*! No. of post-twist equil&
& blocks !* Number of blocks in which to divide a re-equilibration whe&
&n doing a twist-averaged DMC run. Currently, a re-equilibration only &
&takes place when the twist angle is changed.  Electron(-hole) fluid p&
&hases only. For more details, see the DMC_NTWIST keyword.")
call xyzzyaaad1("dmc_reweight_conf","L:I","*! Update walker weights re&
&ad in from config.in !* Weights of walkers are recomputed after readi&
&ng config.in to correct for a modified wave function. This allows con&
&tinuous QMC-MD computations as described in PhysRevLett.94.056403.")
call xyzzyaaad1("dmc_spacewarping","L:I","*! Adjust electron coordinat&
&es to new wave function !* Electronic positions are adjusted to follo&
&w the ionic positions when adapting an existing population to a new w&
&ave function. The method follows the description in PhysRevB.61.R1629&
&1 and is typically combined with DMC_REWEIGHT_CONF.")
call xyzzyaaad1("emin_xi_value","D:E","*! Emin semi-orthog param !* Du&
&ring energy minimization, this parameter determines the wave function&
& with respect to which the linear basis of first derivatives is semi-&
&orthogonalized. Please see the manual for more details.")
call xyzzyaaad1("ebest_av_window","I:E","*! Av. window for DMC equil !&
&* During DMC equilibration the best estimate of the ground-state ener&
&gy is taken to be the average local energy over the last EBEST_AV_WIN&
&DOW moves.")
call xyzzyaaad1("xc_corr_method","I:E","*! XC finite size corr algorit&
&hm !* If set to 1, the XC finite-size correction will be evaluated by&
& determining the coefficient of k^2 in the structure factor; if set t&
&o 2, the XC correction will be evaluated by fitting the whole structu&
&re factor.  Method 1 (the default) is recommended.")
call xyzzyaaad1("emin_min_energy","P:I","*! Emin min. E threshold !* T&
&he stability of energy minimization can be improved by supplying a th&
&reshold below which energies will be ignored (difficult optimizations&
& can produce spurious too-low energies). If this keyword is not prese&
&nt in the input file, no threshold will be applied. This is a keyword&
& of type 'Physical' hence you need to supply units, such as 'ev', 'ry&
&', 'hartree', 'kcal/mol' etc.")
call xyzzyaaad1("expval_error_bars","L:I","*! Error bars for exp value&
&s !* If this flag is set, CASINO will, where practicable, accumulate &
&the additional quantities required to evaluate error bars on requeste&
&d expectation values. This will increase the size of the expval.data &
&file and slow down the calculation slightly. At present this function&
&ality is limited to : structure_factor.")
call xyzzyaaad1("bf_sparse","L:E","*! Woodbury formula det update !* B&
&F_SPARSE can be used to speed up calculations on large systems by usi&
&ng the Woodbury update formula instead of recalculating the determina&
&nts. Let N=number of electrons per determinant, and M=number of elect&
&rons within the backflow range of a given electron (on average). The &
&Woodbury update scales as M * N**2, whereas recomputing the determina&
&nt scales as N**3. It is advantageous to turn BF_SPARSE to T when M/N&
& < 1/3 for acceptance ratios of about 1/2. For calculating non-local &
&energies using the Woodbury formula is advantageous when M/N < 1/2. N&
&ote that when pairing wave functions are used scaling is rather worse&
& than this. If for a given system M is a constant (backflow range is &
&independent of N), using this feature makes backflow calculations sca&
&le as N**3, like Slater-Jastrow calculations. The default value for B&
&F_SPARSE is F.")
call xyzzyaaad1("opt_strict","L:I","*! Halt opt on apparent divergence !&
&* Setting OPT_STRICT=T will cause CASINO to stop a vmc_opt or opt_vmc&
& run if the VMC energies are incremented within a 99.7% confidence in&
&terval during two consecutive cycles. Intended for not wasting CPU ti&
&me in times of scarcity. Default value is F.")
call xyzzyaaad1("opt_noctf_cycles","I:E","*! Fix cutoffs for some cycl&
&es !* Supplying a positive integer X for this keyword will cause all &
&'shallow' parameters (cut-off lengths in the Jastrow factor, backflow&
& transformation and orbitals) to remain fixed for the final X cycles &
&of a multi-cycle optimization run.  This is potentially useful for en&
&ergy minimization, which can be adversely affected by the presence of&
& optimizable cut-off parameters.  OPT_NOCTF_CYCLES defaults to 0 (i.e&
&., cut-offs are never fixed).")
call xyzzyaaad1("mom_den","L:I","*! Accumulate momentum density !* If &
&set to T the momentum density will be accumulated. Exclusively for HE&
&Gs at the moment.")
call xyzzyaaad1("dtvmcs","B:E","*! List of VMC time steps !* Use this &
&keyword to specify a VMC time step for each particle family explicitl&
&y, as well as to determine whether to optimize each of them individua&
&lly. The contents of this block override the values of DTVMC and OPT_&
&DTVMC. One line is to be written for each 'family' of particles, the &
&format of each line being:$ <dtvmc> <opt_dtvmc>$ where <dtvmc> is the&
& value of the time step, and <opt_dtvmc> can be 0 or 1, indicating wh&
&ether to optimize the corresponding <dtvmc> or not.")
call xyzzyaaad1("initial_config","B:E","*! Initial VMC positions !* Us&
&e this keyword if you want to specify the initial VMC configuration t&
&o use instead of the random one generated by the POINTS routine. It i&
&s possible to specify the positions of only some of the particles. Th&
&e format of each line in this block is:$ <spin> <number> <x> <y> <z>$&
& where <spin> is the spin index of the particle, <number> is the inde&
&x of the particle within its spin channel, and <x> <y> <z> is the pos&
&ition of the particle.")
call xyzzyaaad1("dmc_decorr_period","I:I","*! DMC decorrelation period !&
&* Length of the inner decorrelation loop in DMC. The algorithm will p&
&erform DMC_DECORR_PERIOD configuration moves between successive evalu&
&ations of expectation values other than the energy. Setting DMC_DECOR&
&R_PERIOD to a value greater than 1 should reduce the serial correlati&
&on of the data.$$ Notice that DMC_DECORR_PERIOD differs from its VMC &
&counterpart in that in DMC local energies are calculated at intermedi&
&ate steps (they must), and these additional values are averaged into &
&the energy data. Therefore, for calculations which do not require exp&
&ectation values other than the energy, changing DMC_DECORR_PERIOD fro&
&m 1 to some value x is equivalent to multiplying both DMC_[EQUIL|STAT&
&S]_NSTEP and DMC_AVE_PERIOD by x.  In a preliminary DMC calculation D&
&MC_DECORR_PERIOD specifies the frequency with which configurations ar&
&e written out.")
call xyzzyaaad1("dmc_ave_period","I:I","*! Energy-averaging period in &
&DMC !* Number of consecutive local energies that are averaged togethe&
&r in DMC before writing them to the dmc.hist file. The only effect of&
& this keyword is reduce the number of lines in dmc.hist by a factor o&
&f 1/DMC_AVE_PERIOD. Note that if DMC_EQUIL_NSTEP or DMC_STATS_NSTEP a&
&re not divisible by DMC_AVE_PERIOD, they will be rounded up to the ne&
&arest integer multiple of it.")
call xyzzyaaad1("forces","L:I","*! Calculate atomic forces in VMC or D&
&MC !* Forces are only implemented for the Gaussian basis set and are &
&only supposed to work with pseudopotentials in order to eliminate the&
& electron-nucleus singularity. The keyword forces_info can be chosen &
&to vary the level of output.")
call xyzzyaaad1("forces_info","I:E","*! Forces information level !* Co&
&ntrols the amount of information calculated/displayed during calculat&
&ions:$ 2 = displays no extra information; the Hellmann-Feynman force &
&is evaluated with the d-channel of the pseudopotential chosen local a&
&nd s-d and p-d channels nonlocal$ 5 = calculates/displays two additio&
&nal Hellmann-Feynman force estimators where the s- or p- channels of &
&the pseudopotential components are chosen local.")
call xyzzyaaad1("checkpoint","I:E","*! Checkpointing level !* This int&
&eger-valued keyword determines how much CASINO should worry about sav&
&ing checkpoint data to config.* files (which can take a  significant &
&amount of time, especially with large systems done on many cores and &
&can reduce the parallel efficiency - since the slower blocking redist&
&ribution algorithm must be used at the end of every block when we wri&
&te out a config file). CHECKPOINT can take four values:$$ '2' : save &
&data after every block in both VMC and DMC, and save the state  of th&
&e random number generator in OPT runs.$$ '1' [default] : as '2', but &
& save data in VMC only after the last block when RUNTYPE=vmc_opt, opt&
&_vmc or vmc_dmc (still after every block if RUNTYPE=vmc).$$ '0' : onl&
&y save data at the end of the run, for continuation purposes. This is&
& safe only if used in conjunction with the MAX_CPU_TIME keyword (sinc&
&e then the config file will be automatically written out if CASINO se&
&es the job is about to run into an imposed time limit, even if we hav&
&e not completed the full number of requested blocks). $$ '-1' : do no&
&t write config file at all, ever. Note this value should be chosen on&
&ly if you *know* that the job will fit in any imposed time limit , an&
&d that such a run will be long enough to give an acceptably small err&
&or bar, since it will be impossible to subsequently continue the run.&
&$$ CHECKPOINT=0 or -1 clashes with the DMC catastrophe-recovery facil&
&ity, for which each DMC block needs to be checkpointed. The value of &
&CHECKPOINT is thus set to 1 regardless of the input value if DMC_TRIP&
&_WEIGHT > 0 .")
call xyzzyaaad1("fix_holes","L:E","*! Choose constraint for BIEX3 !* T&
&his keyword is used to define the reference points for the exciton-ex&
&citon separation when using BIEX3. Setting fix_holes to T means that &
&the two holes are fixed at a distance xx_sep apart. The default is F,&
& in which case the centres of mass of the two excitons are fixed inst&
&ead. If BIEX3 is not being used this keyword is ignored.")
call xyzzyaaad1("allow_nochi_atoms","L:E","*! Permit atoms no chi (etc&
&) terms !* If this keyword is set to T then CASINO will issue a warni&
&ng message when some atoms are not included in any sets of chi or f t&
&erms in the Jastrow factor and mu and Phi terms in the backflow funct&
&ion.  Otherwise, CASINO halts with an error if some atoms are not inc&
&luded in these terms.")
call xyzzyaaad1("psi_s","T:I","*! Wave function form for Psi_S !* This&
& keyword can take one of the following values:$ 'none': sets Psi_S to&
& one$ 'slater': uses Slater determinants or multi-determinant expansi&
&on (default)$ 'exmol': uses a specially-crafted wave function for exc&
&itonic and positronic molecules$ 'geminal': uses a single geminal$ 'm&
&ahan': uses custom wave function for impurity-in-HEG calculations.")
call xyzzyaaad1("sp_blips","L:E","*! Blip orbitals single prec coeffs !&
&*  Single particle orbitals that appear in the determinants can take &
&a great deal of memory when expanded in a blip basis. With SP_BLIPS=T&
& CASINO represents the blip coefficients using single precision real/&
&complex numbers, which will halve the memory required. This parameter&
& is only relevant when ATOM_BASIS_TYPE=blip. Default value is F.")
call xyzzyaaad1("write_binary_blips","L:E","*! Write formatted bwfn.da&
&ta as binary bwfn.data.bin!*  The formatted blip data file bwfn.data &
&can be very large. Consequently reading this file can be slow. Settin&
&g WRITE_BINARY_BLIPS to T will cause CASINO to write the binary file &
&bwfn.data.bin provided there are no pre-existing bwfn.data.bin files &
&in the run directory. When CASINO is run again it will first attempt &
&to read bwfn.data.bin rather then read bwfn.data, so start up will be&
& faster. If one wishes to change the occupied orbitals delete the exi&
&sting binary file bwfn.data.bin in the run directory. This parameter &
&is only relevant when the ATOM_BASIS_TYPE=blip. Default value is T.")
call xyzzyaaad1("conv_binary_blips","L:E","*! Convert old bwfn.data.b1&
& to new bwfn.data.bin. !*  In November 2011, a new format binary blip&
& file (bwfn.data.bin) was introduced, which is now written out by def&
&ault when CASINO reads in formatted bwfn.data files. The previous bin&
&ary format - bwfn.data.b1 - is still supported, not least because at &
&the time of the introduction of the new format, DFT codes such as PWS&
&CF still produced the old-format b1 file natively (without the interm&
&ediate formatted file ever having existed). By default, b1 files are &
&treated exactly as bin files. If the value of CONV_BINARY_BLIPS is se&
&t to T, then after reading in a b1 file, the data will be converted a&
&nd written out as bwfn.data.bin and the old b1 file will be deleted, &
&prior to continuing the calculation as normal. This can save disk spa&
&ce since bin files are generally smaller than b1 files, and they can &
&be read in somewhat faster (which is advantageous if the same bin fil&
&e is to be used in multiple calculations). The bin files are also mor&
&e portable. Default value is F.")
call xyzzyaaad1("blip_mpc","L:E","*! Blip long-range part of MPC !*  I&
&f BLIP_MPC is set to T and one is using the MPC interaction in a syst&
&em that is periodic in all three dimensions and consists of only elec&
&trons, then the long-range portion of the MPC potential will be evalu&
&ated using three-dimensional B-splines, 'Blips'. In some systems sett&
&ing this to T can greatly speed up the calculation. The default is F.&
&")
call xyzzyaaad1("small_transfer","L:E","*! Prevent transfer of large a&
&rrays !*  If SMALL_TRANSFER is set to T, the DBAR matrices and any po&
&tentially large optional data are not transferred between processors &
&in DMC configuration redistribution. The default is F. Set to T if yo&
&u run into problems with parallel transfers.")
call xyzzyaaad1("opt_small_buffers","L:E","*! Prevent buffering large &
&arrays in opt !*  If OPT_SMALL_BUFFERS is set to T, large intermediat&
&e arrays will not be buffered during optimization.  This only affects&
& optimizations involving determinant coefficients in the absence of b&
&ackflow.")
call xyzzyaaad1("opt_plan","B:I","*! Multi-cycle optimization plan !* &
& This block allows specifying different parameters for each optimizat&
&ion cycle for RUNTYPE = 'vmc_opt', 'opt_vmc' or 'opt'.  The block has&
& one line per optimization cycle (the block length overrides the valu&
&e of OPT_CYCLES), each containing the cycle index followed by any num&
&ber of blank-separated '<keyword>=<value>' assignments.$$ Valid keywo&
&rds are:$ * method: sets OPT_METHOD to <value> for the cycle (string)&
&$ * jastrow: sets OPT_JASTROW to <value> for the cycle (Boolean)$ * b&
&ackflow: sets OPT_BACKFLOW to <value> for the cycle (Boolean)$ * det_&
&coeff: sets OPT_DET_COEFF to <value> for the cycle (Boolean)$ * orbit&
&als: sets OPT_ORBITALS to <value> for the cycle (Boolean)$ * geminal:&
& sets OPT_GEMINAL to <value> for the cycle (Boolean)$ * maxiter: sets&
& OPT_MAXITER to <value> for the cycle (Boolean)$ * fix_cutoffs: deter&
&mines whether to fix cut-offs (T) or not (F) for the cycle (Boolean; &
&analogous to OPT_NOCTF_CYCLES)$$ Input keywords will remain at their &
&provided/default values for all cycles for which they are not modifie&
&d by the corresponding OPT_PLAN line.")
call xyzzyaaad1("rmc_rep_length","I:B","*! No of configs in reptile !*&
& Determines the number of configurations in a reptile for reptation."&
&)
call xyzzyaaad1("dtrmc","D:B","*! RMC time step !* Reptation time step&
&.")
call xyzzyaaad1("rmc_move_length","I:I","*! No of configs in RMC move !&
&* Determines the number of configurations in an RMC move (should be 1&
& if RMC_BOUNCE=T). Default is 1.")
call xyzzyaaad1("rmc_bounce","L:I","*! Use bounce algorithm for RMC !*&
& Determines if bounce algorithm is used.  Default is T.")
call xyzzyaaad1("rmc_meas_pos","L:I","*! Measure electron pos for RMC !&
&* Determines if electron positions are measure in RMC. Default is F."&
&)
call xyzzyaaad1("rmc_decorr_period","I:I","*! RMC decorrelation period !&
&* Length of the inner decorrelation loop in RMC. The algorithm will p&
&erform RMC_DECORR_PERIOD configuration moves between successive evalu&
&ations of the local energy and other expectation values. Setting RMC_&
&DECORR_PERIOD to a value greater than 1 should reduce the serial corr&
&elation of the data, but notice that the length of the run will be in&
&creased.")
call xyzzyaaad1("rmc_ave_period","I:I","*! Energy-averaging period in &
&RMC !* Number of consecutive local energies that are averaged togethe&
&r in RMC before writing them to the rmc.hist file. The only effect of&
& this keyword is reduce the number of lines in rmc.hist by a factor o&
&f 1/RMC_AVE_PERIOD. NOTE: THIS FUNCTIONALITY IS NOT YET IMPLEMENTED -&
& RMC_AVE_PERIOD MUST BE SET TO ONE AT PRESENT.")
call xyzzyaaad1("rmc_equil_nstep","I:B","*! No of steps in RMC equil !&
&* Total number of RMC steps performed in the RMC equilibration stage.&
& Notice that this number will be rounded up to the nearest multiple o&
&f RMC_EQUIL_NBLOCK times RMC_AVE_PERIOD.")
call xyzzyaaad1("rmc_equil_nblock","I:I","*! No of blocks in RMC equil !&
&* Number of blocks into which the total RMC equilibration run length &
&is divided. The value of RMC_EQUIL_NBLOCK determines how often the ou&
&tput file is written to.$$ Default: 1.")
call xyzzyaaad1("rmc_stats_nstep","I:B","*! Number of steps in RMC sta&
&ts accum !* Total number of RMC steps performed in the RMC statistics&
&-accumulation stage. Notice that this number will be rounded up to th&
&e nearest multiple of RMC_STATS_NBLOCK times RMC_AVE_PERIOD.")
call xyzzyaaad1("rmc_stats_nblock","I:I","*! No of blocks in RMC stats&
& accum !* Number of blocks into which the total RMC statistics-accumu&
&lation run is divided. The value of RMC_STATS_NBLOCK determines how o&
&ften the output file is written to.$$ Default: 1.")
call xyzzyaaad1("dmc_init_eref","P:E","*! Initial reference energy !* &
& If set, DMC_INIT_EREF defines the initial reference energy for a DMC&
& calculation. If unset, the VMC energy is used instead (default). Thi&
&s keyword is ignored if the initial configurations come from DMC, in &
&which case the previous DMC best estimate of the energy is used inste&
&ad. This is a keyword of type 'Physical' hence you need to supply uni&
&ts, such as 'ev', 'ry', 'hartree', 'kcal/mol' etc.")
call xyzzyaaad1("use_gjastrow","L:E","*! Use 'gjastrow' Jastrow functi&
&on !*  Use the 'gjastrow' Jastrow factor (T) or the Drummond-Towler-N&
&eeds Jastrow factor (F).  CASINO automatically detects the presence o&
&f JASTROW blocks in the parameters.casl and correlation.data files to&
& initialize the value of this keyword.  Explicitly setting the value &
&of USE_GJASTROW is useful when both correlation.data and parameters.c&
&asl are present.")
call xyzzyaaad1("use_gbackflow","L:E","*! DEV: use gbackflow !*  RESER&
&VED KEYWORD FOR DEVELOPMENT, NO EFFECT.")
call xyzzyaaad1("gen_gjastrow","L:E","*! Convert to gjastrow !*  Setti&
&ng this flag to T triggers the conversion of a Jastrow factor read fr&
&om correlation.data into a gjastrow, which is written to parameters.c&
&asl_converted .  Renaming this file to parameters.casl will cause CAS&
&INO to use the resulting gjastrow in subsequent runs. The default val&
&ue of this keyword is F.")
call xyzzyaaad1("checkpoint_ncpu","I:E","*! Num cores checkpt read gro&
&ups !*  This keyword can be used to specify how to group CPUs for rea&
&ding 'config.in' checkpoint files. Having many CPUs access the same f&
&ile at the same time is not a good idea; therefore we form groups of &
&CHECKPOINT_NCPU cores in which only one of them accesses data. The de&
&fault value is the total number of cores (NNODES internally), but dep&
&ending on the hardware you run on you may want to set CHECKPOINT_NCPU&
& to a different value (between 1 and NNODES). Note that in the case t&
&hat NNODES is not exactly divisible by CHECKPOINT_NCPU, then the rema&
&inder will be distributed over the existing groups, and some of the g&
&roups will therefore contain CHECKPOINT_NCPU+1 cores.")
call xyzzyaaad1("emin_auto_varmin","L:E","*! Min variance in 1st emin &
&cycle !*  EMIN is known to have trouble in the first cycle of optimiz&
&ing a Jastrow factor, where the initial configurations are generated &
&using HF-VMC and the Jastrow (which only contains a cusp initially) i&
&s switched on. Setting EMIN_AUTO_VARMIN to T (default) allows the EMI&
&N module to minimize the variance in the first cycle under such condi&
&tions and then switch back to optimizing the energy in later cycles [&
&FUNCTIONALITY CURRENTLY DISABLED].")
call xyzzyaaad1("contact_den","L:I","*! Accumulate spatial overlap !* &
& If this flag is set then CASINO will accumulate the spatial overlap &
&between electrons and a positron.")
call xyzzyaaad1("dtvmc_shift","D:E","*! VMC transition probability shi&
&ft !* DTVMC_SHIFT is an optional shift in the VMC transition probabil&
&ity which can be used to 'encourage' electrons to be more mobile. DTV&
&MC_SHIFT is expressed in units of the square root of DTVMC.")
call xyzzyaaad1("dmc_nconf_prelim","I:E","*! # of configs to generate &
&in prelim DMC calc !* This is the approximate number of configuration&
&s to generate in a preliminary DMC calculation.")
call xyzzyaaad1("redist_grp_size","I:E","*! Number of processors in re&
&dist group !* In the branch_and_redist algorithm (which does redistri&
&bution of configs across cores in DMC) we must decide which pairs of &
&cores are involved in config transfers, and how many configs are to b&
&e transferred in each operation. There is an optimal algorithm for do&
&ing this (involving looking at individual config multiplicities and t&
&he exact excess or deficit of configs relative to a target on each co&
&re). If we consider *all* the cores, then this algorithm scales linea&
&rly with the number of cores, eventually becoming so expensive that f&
&or a fixed number of configs the code actually becomes slower if we i&
&ncrease the number of cores. We therefore parallelize the algorithm; &
&to do with this we form groups of processors ('redist groups') of siz&
&e REDIST_GRP_SIZE (plus some remainder). When calculating the vector &
&of instructions, only transfers within these groups are contemplated,&
& and the cost for working out what to send where no longer increases &
&with the number of cores (above a certain size).")
call xyzzyaaad1("allow_slave_write","L:I","*! Toggle slave write to ou&
&tput !* The ALLOW_SLAVE_WRITE flag can be used to allow/disallow slav&
&e processor output to the main output file. The default is to allow i&
&t. The ability to turn this off can be useful when you''re running on&
& a million cores.")
call xyzzyaaad1("dmc_md","L:I","*! DMC molecular dynamics !* If DMC_MD&
& is T then in a DMC calculation we assume we are doing molecular dyna&
&mics and that we are restarting from a converged wave function for a &
&slightly different nuclear configuration. In practice, all this means&
& is that the number of steps performed are given by DMCMD_EQUIL_NSTEP&
& and DMCMD_STATS_NSTEP, rather than DMC_EQUIL_NSTEP/DMC_STATS_NSTEP (&
&the number of blocks is assumed to be 1 in the MD case, and the value&
& of BLOCK_TIME is ignored). The number of moves necessary will be gre&
&atly reduced from the normal case. See also DMC_REWEIGHT_CONF and DMC&
&_SPACEWARPING. The necessary manipulations are automated by the runqm&
&cmd script.")
call xyzzyaaad1("dmcmd_equil_nstep","I:B","*! No of steps in DMC-MD eq&
&uil !* Total number of DMC steps performed in the DMC equilibration s&
&tage when we are doing a non-initial step in a DMC molecular dynamics&
& calculation (we already have a quasi-converged wave function for a s&
&lightly different nuclear configuration). The number of blocks is ass&
&umed to be 1.")
call xyzzyaaad1("dmcmd_stats_nstep","I:B","*! No of steps in DMC-MD st&
&ats accum !* Total number of DMC steps performed in the DMC statistic&
&s-accumulation stage when we are doing a non-initial step in a DMC mo&
&lecular dynamics calculation (we already have a quasi-converged wave &
&function for a slightly different nuclear configuration). The number &
&of blocks is assumed to be 1.")
call xyzzyaaad1("rng_restart_safe","L:E","*! Continuity of RNG through&
& restart !* We would like e.g. a 1000 move VMC run, and two 500 move &
&VMC runs linked together by a restart, to give the same answer (in th&
&e sense that we end up with the same vmc.hist file). Unfortunately th&
&ey do not in general since the pseudorandom number sequence is affect&
&ed by the restart. This is because random numbers are generated somet&
&hing like 63 at a time and stored in a buffer until needed (this buff&
&er being refilled when necessary). In the normal way of saving a poin&
&t in the random number sequence, any unused numbers in the buffer are&
& discarded, which means the final answer will be different to the unr&
&estarted case. If the keyword RNG_RESTART_SAFE is T (which is actuall&
&y now the default) then the whole current buffer is saved in the fina&
&l config.out file as well as the current state of the random sequence&
& (necessarily fixed at the end of the current buffer). This allows mu&
&ltiple step runs to give the same answer as single step runs, at the &
&expense of slightly larger config files.")
call xyzzyaaad1("shm_size_nproc","I:E","*! No of assumed MPI processes&
& in Shm test !* In Shm calculations on Blue Gene machines one needs t&
&o know in advance the number of MB of shared memory required, so that&
& one may set the BG_SHAREDMEMSIZE environment variable (which can be &
&done by means of the --user.shemsize argument to runqmc). CASINO will&
& calculate this number and print it to output at the end of the setup&
& process (within the scope of TESTRUN=T). However, the amount of shar&
&ed memory required depends on the the number of MPI processes per nod&
&e (or per shared memory partition). If one ultimately wishes to run o&
&n, say, half a million cores, it may be desirable to execute a test r&
&un on just a few cores on your personal laptop, rather than waiting a&
& week for the full job to sit in a queue. For the purposes of computi&
&ng the size of the shared memory partition, one may therefore set the&
& number of desired processes/node by setting SHM_SIZE_NPROC, and this&
& value will be used in computation of the shared memory size rather t&
&han the actual number of processes/node being used in the test run (u&
&nless SHM_SIZE_NPROC=0 - which is the default). Note that the CASINO &
&test run must be done in Shm mode.")
call xyzzyaaad1("vmc_sampling","T:E","*! Type of VMC sampling distribu&
&tion !* This keyword allows using alternative sampling distributions &
&instead of the square of the trial wave function in VMC and wave func&
&tion optimization.  Each sampling distribution has its own set of adv&
&antages and disadvantages.$$ Possible values of this keyword are:$$ -&
& 'standard': use the square of the wave function (default).  This is &
&the usual way of running VMC calculations.  This form has the drawbac&
&k of poor sampling near the nodes, which negatively impacts optimizat&
&ion. This form should be used when generating configurations for DMC.&
&$$ - 'optimum': use the optimum sampling distribution, which achieves&
& the smallest possible variance of the local energies.  This form enh&
&ances sampling near the nodes of the trial wave function, potentially&
& improving optimization.  This form is expensive because it requires &
&the local energy to be evaluated at every move.  With this form you s&
&hould set the additional input parameters VMC_OPTIMUM_E0 and VMC_OPTI&
&MUM_EW.$$ - 'HF optimum': use the optimum sampling distribution for t&
&he HF wave function.  This does not necessarily reduce the variance w&
&ith respect to 'standard', but otherwise offers the same advantages a&
&s 'optimum' at a much reduced cost.  Again, you should set the additi&
&onal input parameters VMC_OPTIMUM_E0 and VMC_OPTIMUM_EW.$$ - 'efficie&
&nt': use a probability distribution designed to be inexpensive to eva&
&luate and has just the correct properties to enhance sampling near th&
&e nodes.  This form tends to yield the best performance, but it is on&
&ly applicable to multideterminant wave functions.")
call xyzzyaaad1("vmc_optimum_e0","P:E","*! Centre parameter for optimu&
&m VMC sampling !* This keyword controls the centre parameter used for&
& optimum VMC sampling, which is enabled by setting VMC_SAMPLING to 'o&
&ptimum' or 'HF optimum'. By default this is 0.d0 hartree. It should b&
&e set to an estimate of the ground-state energy of the system under c&
&onsideration. This is a keyword of type 'Physical' hence you need to &
&supply units, such as 'ev', 'ry', 'hartree', 'kcal/mol' etc.")
call xyzzyaaad1("vmc_optimum_ew","P:E","*! Width parameter for optimum&
& VMC sampling !* This keyword controls the width parameter used for o&
&ptimum VMC sampling, which is enabled by setting VMC_SAMPLING to 'opt&
&imum' or 'HF optimum'. By default this is 100.d0 hartree. It should b&
&e set to an estimate of the expected width of the local energy distri&
&bution. This is a keyword of type 'Physical' hence you need to supply&
& units, such as 'ev', 'ry', 'hartree', 'kcal/mol' etc.")
call xyzzyaaad1("population","L:I","*! Accumulate ionic populations !*&
& If set to T, ionic populations will be evaluated by Voronoi partitio&
&ning of the charge density.")
call xyzzyaaad1("block_time","P:I","*! Approx CPU time per block !* If&
& BLOCK_TIME is greater than 0.d0, then the number of blocks of moves &
&implied by VMC_NBLOCK, DMC_EQUIL_NBLOCK, or  DMC_STATS_NBLOCK will be&
& ignored. Instead, CASINO will do everything it normally does at the &
&end of a block approximately every BLOCK_TIME seconds of CPU time.$$ &
&For VMC, the actions performed after a block are: (1) write data to o&
&ut, vmc.hist, and possibly expval.data; (2) write current VMC state p&
&lus any accumulated configs to config.out (this latter only if CHECKP&
&OINT is increased to 2 from its default of 1 - otherwise config.out i&
&s only written at the end of the last block).$$ For DMC, the actions &
&performed after a block are: (1) write data to out, dmc.hist, and pos&
&sibly expval.data (the latter not during equilibration); (2) make a b&
&ackup copy of the config.out/expval.data file (if catastrophe protect&
&ion is turned on via DMC_TRIP_WEIGHT); (3) Write the dmc.status file &
&(except after the last block); (4) Write current state of the system,&
& and all configs in the current population to config.out (note that b&
&y setting CHECKPOINT to 0, this step can be skipped until the end of &
&the last block, or skipped completely if CHECKPOINT=-1, but this is n&
&ot the default).$$ Note the above actions can take a long time, espec&
&ially if they involve writing to disk, so it is better to do them as &
&infrequently as possible (i.e. large value of BLOCK_TIME). Obviously &
&if the stopping criterion (number of moves, target error bar..) impli&
&es that the run will stop before BLOCK_TIME minutes have elapsed, the&
&n the total run time can be shorter than BLOCK_TIME. Note using BLOCK&
&_TIME implies that multiple repetitions of the same run will not nece&
&ssarily lead to the same answer in parallel calculations (as the numb&
&er of runs done in BLOCK_TIME seconds is defined by the master and wh&
&at happens on the other slaves can mess around with their random numb&
&er sequences in an unpredictable way).$$ Note finally that BLOCK_TIME&
& is a physical parameter with dimensions of time; the units must be s&
&pecified as e.g. 1 day, 24 hr, 1440 min, or 86400 s.")
call xyzzyaaad1("stop_method","T:B","*! How to terminate VMC/DMC run !&
&* The STOP_METHOD keyword defines how a VMC/DMC run is to be terminat&
&ed. It may take the values 'nstep', 'target_error', or 'small_error'.&
&$$ The classic method is 'nstep' which means simply: perform the numb&
&er of VMC/DMC steps implied by the input keywords VMC_NSTEP or DMC_ST&
&ATS_NSTEP then stop. The error bar then is what it is (it may be too &
&large or smaller than required).$$ If STOP_METHOD = 'target_error' th&
&en the run will continue until the error bar on the total energy (cor&
&rected on the fly for serial correlation) is approximately equal to t&
&hat defined by the TARGET_ERROR input keyword, subject to the constra&
&int that the *estimated* CPU time required on the master process (sum&
&med over restarts if necessary) will not exceed STOP_TIME. CASINO is &
&able to approximately estimate the required time by analysing how the&
& error bar decreases as a function of the number of moves, and as soo&
&n as it is reasonably confident that the desired target_error is too &
&small and cannot be reached, then the code will stop (in a restartabl&
&e condition). On halting in this manner, an estimate of the CPU time &
&required to get a range of error bars will be written to the output f&
&ile. Note that the method used to estimate the required time assumes &
&the validity of the central limit theorem, which is only approximatel&
&y valid in most cases.$$ If STOP_METHOD = 'small_error', CASINO will &
&attempt to make the error bar as small as possible in a 'reasonable t&
&ime' defined by the value of the STOP_TIME keyword. 'As small as poss&
&ible' means what it says, but taking account of the fact that there i&
&s an error bar on the error bar and it is somewhat pointless to reduc&
&e the error bar below its significant precision.$$ Note in both the l&
&ast two cases CASINO has a minimum run length needed to get a reasona&
&ble estimate of the variance.")
call xyzzyaaad1("target_error","D:B","*! Target error bar on the energ&
&y !* If STOP_METHOD = 'target_error' then the run will continue until&
& the error bar on the total energy (corrected on the fly for serial c&
&orrelation) is approximately equal to that defined by the TARGET_ERRO&
&R keyword, subject to the constraint that the *estimated* CPU time re&
&quired on the master process (summed over restarts if necessary) will&
& not exceed STOP_TIME. CASINO is able to approximately estimate the r&
&equired time by analysing how the error bar decreases as a function o&
&f the number of moves, and as soon as it is reasonably confident that&
& the desired target_error is too small and cannot be reached, then th&
&e code will stop (in a restartable condition). On halting in this man&
&ner, an estimate of the CPU time required to get a range of error bar&
&s will be written to the output file. Note that the method used to es&
&timate the required time assumes the validity of the central limit th&
&eorem, which is only approximately valid in most cases.")
call xyzzyaaad1("stop_time","P:B","*! Stop method definition of 'reaso&
&nable time'  !* If STOP_METHOD = 'target_error' then the run will con&
&tinue until the error bar on the total energy (corrected on the fly f&
&or serial correlation) is approximately equal to that defined by the &
&TARGET_ERROR input keyword, subject to the constraint that the *estim&
&ated* CPU time required on the master process (summed over restarts i&
&f necessary) will not exceed STOP_TIME. CASINO is able to approximate&
&ly estimate the required time by analysing how the error bar decrease&
&s as a function of the number of moves, and as soon as it is reasonab&
&ly confident that the desired target_error is too small and cannot be&
& reached, then the code will stop (in a restartable condition). On ha&
&lting in this manner, an estimate of the CPU time required to get a r&
&ange of error bars will be written to the output file. Note that the &
&method used to estimate the required time assumes the validity of the&
& central limit theorem, which is only approximately valid in most cas&
&es.$$ If STOP_METHOD = 'small_error', CASINO will attempt to make the&
& error bar as small as possible in a 'reasonable time' defined by the&
& value of STOP_TIME. 'As small as possible' means what it says, but t&
&aking account the fact that there is an error bar on the error bar an&
&d it is somewhat pointless to reduce the error bar below its signific&
&ant precision.")
call xyzzyaaad1("twop_dm_mom","L:I","*! Accum 2p momentum density !* I&
&f TWOP_DM_MOM is set to T, then the Fourier transform of the two-part&
&icle density matrix will be computed.  This is only possible if the s&
&ystem is homogeneous for the moment, and will increase the cost of th&
&e calculation significantly.")
call xyzzyaaad1("cond_fraction_mom","L:I","*! Accum strict 2p momentum&
& density !* If COND_FRACTION_MOM is set to T, then an improved estima&
&tor of the Fourier transform of the two-particle density matrix, from&
& which one-body contributions are subtracted, will be computed.  This&
& is only available if the system is homogeneous, for the moment.")
call xyzzyaaad1("nequil","I:B","*! (OLD) Number of equilibration steps !&
&* NEQUIL is the number of Metropolis equilibration steps.  Note that &
&CORPER is not accounted for, i.e., NEQUIL configuration move attempts&
& are made.$$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT&
& FOR IT WILL BE REMOVED IN EARLY 2015. USE vmc_equil_nstep INSTEAD.]"&
&)
call xyzzyaaad1("nmove","I:B","*! (OLD) Number of VMC moves !* Number &
&of moves per block in the main VMC loop, which corresponds to the num&
&ber of configurations in a block for which the energy and other expec&
&tation values are stored. Note that there is an inner decorrelation l&
&oop of length CORPER, and an additional expectation-value averaging l&
&oop of length NVMCAVE, so the total number of configuration moves att&
&empted in a block is NMOVE*CORPER*NVMCAVE. If PARALLEL_KEYWORDS is se&
&t to 'per_node', then NMOVE is a per-processor quantity, and the tota&
&l number of moves per block in the main VMC loop is NMOVE times the n&
&umber of processors. If PARALLEL_KEYWORDS is 'total', then NMOVE is a&
& total quantity.$$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATED - S&
&UPPORT FOR IT WILL BE REMOVED IN EARLY 2015. USE vmc_nstep INSTEAD.]"&
&)
call xyzzyaaad1("nblock","I:B","*! (OLD) Number of VMC blocks !* NBLOC&
&K is the total number of blocks of NMOVE moves in a VMC run.  NB, you&
& should use the reblock utility to investigate the effect of varying &
&the block size on the variance after the calculation is completed. In&
& VMC, NBLOCK just determines how often block averaged quantities are &
&written to the output file.$$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEP&
&RECATED - SUPPORT FOR IT WILL BE REMOVED IN EARLY 2015. USE vmc_nbloc&
&k INSTEAD.]")
call xyzzyaaad1("corper","I:B","*! (OLD) VMC energy evaluation period !&
&* VMC only. The local energy is calculated only once every CORPER con&
&figuration moves. If CORPER is sufficiently large then the VMC energy&
& data are uncorrelated, so that the naive error bars displayed in the&
& out file are accurate, and the VMC-generated configurations are unco&
&rrelated, which is very important when performing wave-function optim&
&ization or generating the initial configuration population for a DMC &
&calculation.  If one is simply interested in obtaining a VMC energy t&
&hen CORPER should be 3 or 4; if one is using VMC to generate configur&
&ations for DMC or variance minimization then CORPER should usually be&
& in excess of 10 (e.g. 15 is typical).$$ [NOTE: THIS KEYWORD IS NOW S&
&EVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMOVED IN EARLY 2015. US&
&E vmc_decorr_period INSTEAD.]")
call xyzzyaaad1("nwrcon","I:B","*! (OLD) Number of configs to write !*&
& NWRCON is the number of configurations to be written out in VMC, for&
& later use (wave-function optimization or DMC). If PARALLEL_KEYWORDS &
&is set to 'per_node' (the default) then NWRCON is a per-processor qua&
&ntity, i.e. the total number of configurations written is NWRCON mult&
&iplied by the number of processors.  If PARALLEL_KEYWORDS is set to '&
&total' then NWRCON is the total number of configurations written.$$ [&
&NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WILL B&
&E REMOVED IN EARLY 2015. USE vmc_nconfig_write INSTEAD.]")
call xyzzyaaad1("nconfig","D:B","*! (OLD) DMC target weight !* Target &
&weight in DMC.  This is synonymous with ""target population"", except&
& that NCONFIG is allowed to be a non-integer. If PARALLEL_KEYWORDS is&
& set to 'per_node' (the default) then NCONFIG is a per-processor quan&
&tity, i.e. the total weight is NCONFIG multiplied by the number of pr&
&ocessors.  If PARALLEL_KEYWORDS is set to 'total' then NCONFIG is the&
& total target weight.$$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATE&
&D - SUPPORT FOR IT WILL BE REMOVED IN EARLY 2015. USE dmc_target_weig&
&ht INSTEAD.]")
call xyzzyaaad1("nvmcave","I:I","*! (OLD) Average successive points !*&
& Instead of writing out VMC energies etc every time they are calculat&
&ed, we average over NVMCAVE evaluations before writing to the vmc.his&
&t file.  Note that the total number of moves of all electrons in a bl&
&ock is given by NMOVE*CORPER*NVMCAVE so you will need to reduce NMOVE&
& proportionately if you increase NVMCAVE.$$ [NOTE: THIS KEYWORD IS NO&
&W SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMOVED IN EARLY 2015.&
& USE vmc_ave_period INSTEAD.]")
call xyzzyaaad1("nmove_dmc_equil","I:B","*! (OLD) Number of moves DMC &
&equil !* NMOVE_DMC_EQUIL is the number of moves of all electrons in a&
& block during DMC equilibration.$$ [NOTE: THIS KEYWORD IS NOW SEVEREL&
&Y DEPRECATED - SUPPORT FOR IT WILL BE REMOVED IN EARLY 2015. USE dmc_&
&equil_nstep INSTEAD.]")
call xyzzyaaad1("nblock_dmc_equil","I:B","*! (OLD) Number of blocks DM&
&C equil !* NBLOCK_DMC_EQUIL is the total number of blocks of NMOVE_DM&
&C_EQUIL moves during the DMC equilibration phase.$$ [NOTE: THIS KEYWO&
&RD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMOVED IN EAR&
&LY 2015. USE dmc_equil_nblock INSTEAD.]")
call xyzzyaaad1("nmove_dmc_stats","I:B","*! (OLD) Number of moves DMC &
&stats accum !* NMOVE_DMC_STATS is the number of moves of all electron&
&s in a block during the DMC statistics accumulation phase.$$ [NOTE: T&
&HIS KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMOV&
&ED IN EARLY 2015. USE dmc_stats_nstep INSTEAD.]")
call xyzzyaaad1("nblock_dmc_stats","I:B","*! (OLD) Number of blocks DM&
&C stats accum !* NBLOCK_DMC_STATS is the total number of blocks of NM&
&OVE_DMC_STATS moves during the DMC statistics accumulation phase.  NB&
&, you should use the reblock utility in the utils directory to see th&
&e effect of varying the block size on the variance after the calculat&
&ion is completed.$$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATED - &
&SUPPORT FOR IT WILL BE REMOVED IN EARLY 2015. USE dmc_stats_nblock IN&
&STEAD.]")
call xyzzyaaad1("trip_popn","D:I","*! (OLD) DMC recovery population !*&
& In the course of a DMC simulation, it is possible for a configuratio&
&n ""population explosion"" to occur.  If TRIP_POPN is set to 0 then n&
&othing will be done about this. If TRIP_POPN>0 then it will attempt t&
&o restart the block if the iteration weight exceeds TRIP_POPN.  A gen&
&eral suggestion for its value would be three times NCONFIG (but see t&
&he discussion in the manual about this). If PARALLEL_KEYWORDS is set &
&to 'per_node', then TRIP_POPN is a per-processor quantity, and the th&
&reshold for catastrophe is TRIP_POPN times the number of processors. &
&If PARALLEL_KEYWORDS is 'total', then TRIP_POPN is a total quantity.$&
&$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WIL&
&L BE REMOVED IN EARLY 2015. USE dmc_trip_weight INSTEAD.]")
call xyzzyaaad1("vmc_twist_av","L:I","*! (OLD) Perform VMC twist avera&
&ging !* Perform random changes of twist angle during a VMC simulation&
&.  This can only be done for electron(-hole) fluid phases at present.&
&  NEQUIL_TA must be given a positive value in this case.  The k-vecto&
&r offset is changed at the start of each block in VMC.$$ [NOTE: THIS &
&KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMOVED I&
&N EARLY 2015. USE vmc_ntwist INSTEAD.]")
call xyzzyaaad1("nequil_ta","I:I","*! (OLD) Number of VMC post-twist-c&
&hange equilibration moves !* Number of equilibration VMC moves to mak&
&e after each change of k-vector offset. HEG only.$$ [NOTE: THIS KEYWO&
&RD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMOVED IN EAR&
&LY 2015. USE vmc_reequil_nstep INSTEAD.]")
call xyzzyaaad1("num_dmc_twists","I:I","*! (OLD) Number of DMC twist a&
&ngles !* Number of different offsets to the grid of k vectors to be a&
&pplied during DMC statistics accumulation. HEG only.$$ [NOTE: THIS KE&
&YWORD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMOVED IN &
&EARLY 2015. USE num_dmc_twists INSTEAD.]")
call xyzzyaaad1("nmove_dmct_equil","I:I","*! (OLD) Number of DMC moves&
& per block during post-twist-change equil !* Number of equilibration &
&DMC moves per block to make after each change of k-vector offset. Ele&
&ctron(-hole) fluid phases only.$$ [NOTE: THIS KEYWORD IS NOW SEVERELY&
& DEPRECATED - SUPPORT FOR IT WILL BE REMOVED IN EARLY 2015. USE dmc_r&
&eequil_nstep INSTEAD.]")
call xyzzyaaad1("nblock_dmct_equil","I:I","*! (OLD) Number of blocks o&
&f DMC moves during post-twist-change equil !* Number of blocks of equ&
&ilibration DMC moves after each change of k-vector offset. HEG only.$&
&$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WIL&
&L BE REMOVED IN EARLY 2015. USE dmc_reequil_nblock INSTEAD.]")
call xyzzyaaad1("corper_dmc","I:I","*! (OLD) DMC correlation period !*&
& When gathering expectation values in DMC, it is inefficient to compu&
&te the expectation values at every iteration.  Instead one can calcul&
&ate the expectation values every CORPER_DMCth iteration.  Note that, &
&unlike its VMC counterpart, CORPER_DMC does not affect the number of &
&moves carried out.  In a preliminary DMC calculation, CORPER_DMC spec&
&ifies the frequency with which configurations are written out.$$ [NOT&
&E: THIS KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE R&
&EMOVED IN EARLY 2015. USE dmc_decorr_period INSTEAD.]")
call xyzzyaaad1("ndmcave","I:I","*! (OLD) Average successive points in&
& DMC !* Instead of writing out DMC energies etc every time they are c&
&alculated, we average over NDMCAVE evaluations before writing to the &
&dmc.hist file. Note that the total number of moves of all electrons i&
&n a block is given by NMOVE_DMC_[EQUIL,STATS]*NDMCAVE so you will nee&
&d to reduce NMOVE_DMC_[EQUIL,STATS] proportionately if you increase N&
&DMCAVE.$$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT FO&
&R IT WILL BE REMOVED IN EARLY 2015. USE dmc_ave_period INSTEAD.]")
call xyzzyaaad1("nmove_dmcmd_equil","I:B","*! (OLD) Number of moves DM&
&C-MD equil !* NMOVE_DMC_EQUIL_MD is the number of moves of all electr&
&ons in a block during DMC equilibration, when we are doing a non-init&
&ial step in a DMC molecular dynamics calculation (i.e. we already hav&
&e a quasi-converged wave function for a slightly different nuclear co&
&nfiguration)..$$ [NOTE: THIS KEYWORD IS NOW SEVERELY DEPRECATED - SUP&
&PORT FOR IT WILL BE REMOVED IN EARLY 2015. USE dmcmd_equil_nstep INST&
&EAD.]")
call xyzzyaaad1("nmove_dmcmd_stats","I:B","*! (OLD) Number of moves DM&
&C-MD stats accum !* NMOVE_DMC_STATS_MD is the number of moves of all &
&electrons in a block during the DMC statistics accumulation phase, wh&
&en we are doing a non-initial step in a DMC molecular dynamics calcul&
&ation (i.e. we already have a quasi-converged wave function for a sli&
&ghtly different nuclear configuration)..$$ [NOTE: THIS KEYWORD IS NOW&
& SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMOVED IN EARLY 2015. &
&USE dmcmd_stats_nstep INSTEAD.]")
call xyzzyaaad1("nconfig_prelim","I:E","*! # of configs to generate in&
& prelim DMC calc !* This is the approximate number of configurations &
&per processor to generate in a preliminary DMC calculation.$$ [NOTE: &
&THIS KEYWORD IS NOW SEVERELY DEPRECATED - SUPPORT FOR IT WILL BE REMO&
&VED IN EARLY 2015. USE dmc_nconf_prelim INSTEAD.]")
call xyzzyaaad1("single_precision_blips","L:E","*! REDUNDANT Blip orbi&
&tal single prec coeffs !*  SINGLE_PRECISION_BLIPS is redundant.  Use &
&SP_BLIPS instead.")
call xyzzyaaad1("btype","I:B","*! REDUNDANT: Basis set type !* BTYPE i&
&s redundant.  Use ATOM_BASIS_TYPE instead.")
call xyzzyaaad1("special_wfn","T:E","*! REDUNDANT: Special wave functi&
&on !* SPECIAL_WFN is redundant.  Use ATOM_BASIS_TYPE instead.")
call xyzzyaaad1("iterac","I:I","*! REDUNDANT: ee interaction type !* I&
&TERAC is redundant.  Use INTERACTION instead.")
call xyzzyaaad1("no_ee_int","L:E","*! REDUNDANT: Turn off e-e interact&
&ion !* NO_EE_INT is redundant.  Use INTERACTION instead.")
call xyzzyaaad1("nlrule1","I:I","*! REDUNDANT: NL int rule (VMC/DMC) !&
&* NLRULE1 is redundant. Use NON_LOCAL_GRID instead.")
call xyzzyaaad1("iaccumulate","L:I","*! REDUNDANT: DMC stage !* IACCUM&
&ULATE is redundant. Use RUNTYPE = dmc_equil or dmc_stats instead.")
call xyzzyaaad1("use_molorbmods","L:B","*! REDUNDANT: Molecular-orbita&
&l mods !* USE_MOLORBMODS is redundant. Use USE_ORBMODS instead.")
call xyzzyaaad1("nlrule2","I:I","*! REDUNDANT: NL int rule (configs) !&
&* NLRULE2 is redundant and its value is ignored. See NON_LOCAL_GRID."&
&)
call xyzzyaaad1("calc_variance","L:I","*! REDUNDANT: Calculate VMC var&
&iance !* CALC_VARIANCE is redundant and its value is ignored.")
call xyzzyaaad1("vm_deriv_buffer","L:I","*! REDUNDANT: Buffer WF secti&
&ons !* VM_DERIV_BUFFER is redundant and its value is ignored.")
call xyzzyaaad1("vm_dist_buffer","L:E","*! REDUNDANT: Buffer distances !&
&* VM_DIST_BUFFER is redundant and its value is ignored.")
call xyzzyaaad1("bf_save_memory","L:E","*! REDUNDANT: Disable backflow&
& buffers !* BF_SAVE_MEMORY is redundant and its value is ignored.")
call xyzzyaaad1("emin_sampling","T:E","*! REDUNDANT: Type of sampling &
&in EMIN !* EMIN_SAMPLING is redundant and its value is ignored.")
call xyzzyaaad1("dmc_npops","I:E","*! REDUNDANT: Number of independent&
& populations in DMC !* Removed because it doesn''t help.")
call xyzzyaaad1("spindensitymat","L:I","*! REDUNDANT: Keyword !* SPIND&
&ENSITYMAT is redundant. Implied by DENSITY/SPIN_DENSITY=T in non-coll&
&inear system.")
call xyzzyaaad1("redist_period","I:E","*! REDUNDANT: Redistribution pe&
&riod !* REDIST_PERIOD is redundant and its value is ignored.")
call xyzzyaaad1("num_cpus_in_group","I:I","*! REDUNDANT: Num CPUs to s&
&hare blip orbs !* NUM_CPUS_IN_GROUP is redundant and its value is ign&
&ored.")
call xyzzyaaad1("use_mpiio","L:E","*! REDUNDANT: MPI IO to r/w binary &
&blip file !* USE_MPIIO is redundant and its value is ignored.")
call xyzzyaaad1("have_ae","L:B","*! REDUNDANT: System contains all-ele&
&ctron nuclei !*  HAVE_AE is redundant and its value is ignored.")
call xyzzyaaad1("allow_ae_ppots","L:I","*! REDUNDANT: Allow mixed ae/p&
&p !*  ALLOW_AE_PPOTS is redundant and its value is ignored.")
call xyzzyaaad1("parallel_keywords","T:E","*! REDUNDANT: Parallel keyw&
&ord interpretation !* When PARALLEL_KEYWORDS was set to 'per_node' [d&
&efault], keywords NMOVE, NWRCON, NCONFIG and TRIP_POPN were interpret&
&ed as numbers per node, as in previous versions of CASINO. When PARAL&
&LEL_KEYWORDS was set to 'total', the value of these keywords was divi&
&ded by the number of processors, so that they represent total quantit&
&ies.$$ This was an initial attempt at providing simplified input for &
&parallel runs. However, the preferred approach is now to replace:$ * &
&NMOVE, NWRCON, CORPER, NVMCAVE, NBLOCK, NEQUIL, NCONFIG, NMOVE_DMC_[E&
&QUIL|STATS], NBLOCK_DMC_[EQUIL|STATS], NDMCAVE, CORPER_DMC, TRIP_POPN&
&, VMC_TWIST_AV, NEQUIL_TA, NUM_DMC_TWISTS, NMOVE_DMCT_EQUIL, and NBLO&
&CK_DMCT_EQUIL$ with$ * VMC_NSTEP, VMC_NCONFIG_WRITE, VMC_DECORR_PERIO&
&D, VMC_AVE_PERIOD, VMC_NBLOCK, VMC_EQUIL_NSTEP, DMC_TARGET_WEIGHT, DM&
&C_[EQUIL|STATS]_NSTEP, DMC_[EQUIL|STATS]_NBLOCK, DMC_AVE_PERIOD, DMC_&
&DECORR_PERIOD, DMC_TRIP_WEIGHT, VMC_NTWIST, VMC_REEQUIL_NSTEP, DMC_NT&
&WIST, DMC_REEQUIL_NSTEP, and DMC_REEQUIL_NBLOCK,$ respectively, to ac&
&hieve the same effect (but notice the slightly different meanings of &
&some of the keywords in the two sets!).")
call xyzzyaaae1
end subroutine load_keywords
subroutine xyzzyaaad1(label,typ,dscrpt)
implicit none
character(*),intent(in) :: label
character(*),intent(in) :: typ
character(*),intent(in) :: dscrpt
allocate(xyzzyaaac1)
xyzzyaaac1%kw%label=trim(adjustl(label))
xyzzyaaac1%kw%typ=trim(adjustl(typ))
xyzzyaaac1%kw%dscrpt=trim(adjustl(dscrpt))
if(numkw>0)then
xyzzyaaab1%next=>xyzzyaaac1
else
xyzzyaaaa1=>xyzzyaaac1
endif
xyzzyaaab1=>xyzzyaaac1
xyzzyaaab1%next=>xyzzyaaaa1
numkw=numkw+1
end subroutine xyzzyaaad1
subroutine xyzzyaaae1
implicit none
integer xyzzyaaaa4
allocate(kw(numkw))
xyzzyaaab1=>xyzzyaaaa1
do xyzzyaaaa4=1,numkw-1
kw(xyzzyaaaa4)=xyzzyaaab1%kw
xyzzyaaac1=>xyzzyaaab1%next
deallocate(xyzzyaaab1)
xyzzyaaab1=>xyzzyaaac1
enddo
kw(numkw)=xyzzyaaab1%kw
deallocate(xyzzyaaab1)
end subroutine xyzzyaaae1
end module esdf_key

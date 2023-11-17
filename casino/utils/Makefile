#------------------------------------------------------------------#
# CASINO utils Makefile                                            #
#                                                                  #
# Allows easy maintenance and installation of all tools and        #
# utilities supplied with the CASINO distribution. The Makefile    #
# will go through all the relevant utility directories, compiling  #
# and copying everything necessary into the bin directory          #
# CASINO/bin_qmc which you should include in your path.            #
#                                                                  #
# All system dependent stuff is derived from the include           #
# files in the arch/data directory.                                #
#                                                                  #
# Targets:                                                         #
#  default: Make all utils and stick in relevant places.           #
#  clean  : Removes all object files.                              #
#  vclean : Removes all object and executable files.               #
#                                                                  #
# MDT 1999                                                         #
#------------------------------------------------------------------#
# Basic stuff
UTILDIR = $(PWD)
TOOLDIR = $(UTILDIR)/../build-tools
INCBASE = $(UTILDIR)/../arch/data
BINBASE = $(UTILDIR)/../bin_qmc
LIBBASE = $(UTILDIR)/../lib/zlib
OBJBASE = $(UTILDIR)/zlib
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

start: default
.PHONY: start

# System-independent variables
UTILS_MODE     =
UNATIVE_       =
UNATIVE_native = _NATIVE
BINDIR         = $(BINBASE)/utils/$(CASINO_ARCH)
BINDIR_FEATURE = $(BINDIR)/opt
LIBDIR         = $(LIBBASE)/$(CASINO_ARCH)
OBJDIR         = $(OBJBASE)/$(CASINO_ARCH)
MKBINDIR       = [ -d $(BINDIR) ] || mkdir -p $(BINDIR)
MKDUMMY_nopath = { $(MKBINDIR) ;\
                 bash -c '$(BASH_PRETTY) ; compile_info DUMMY utils $(@F)' ;\
                 rm -f $@ 2> /dev/null ;\
                 touch $(UTILDIR)/dummy_script/dummy_script ;\
                 ln -s $(UTILDIR)/dummy_script/dummy_script $@ ; }

# Defaults for system-dependent variables
ENVIRONMENT_COMMAND        = :
ENVIRONMENT_COMMAND_NATIVE = :

# System-dependent include file
include $(INCBASE)/$(CASINO_ARCH).arch

# Use the same compilers as for src, unless overrides defined
UNATIVE                   = $(UNATIVE_$(UTILS_MODE))
ENVIRONMENT_COMMAND_UTILS = $(ENVIRONMENT_COMMAND$(UNATIVE))
F90_UTILS                 = $(F90$(UNATIVE))
FFLAGS_opt_UTILS          = $(FFLAGS_opt$(UNATIVE))
FFLAGS_all_UTILS          = $(FFLAGS_all$(UNATIVE))
INCLUDE_DIR_UTILS         = $(INCLUDE_DIR$(UNATIVE))
# NB, following only used to determine whether to compile some utils
CC_UTILS                  = $(CC$(UNATIVE))
CXX_UTILS                 = $(CXX$(UNATIVE))

# System-dependent variables from data in include files
FFLAGS_UTILS = $(FFLAGS_opt_UTILS) $(FFLAGS_all_UTILS) $(INCLUDE_DIR_UTILS)

# Main target
.PHONY: default
default: info wrap

# Info
.PHONY: info
info:
	@bash -c '$(BASH_PRETTY) ;\
        info_title "CASINO utils MAKE information" ;\
        info_equal "LN" "ln -s" ;\
        [ "$(ENVIRONMENT_COMMAND_UTILS)" != : ]\
         && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND_UTILS)" ;\
        info_equal "F90" "$(F90_UTILS) $(FFLAGS_UTILS)" ;\
        info_end'

# Types of utilities according to how they are made
UTILS = $(SCRIPT_UTILS) $(SINGLE_F90_UTILS) $(SUBMAKE_UTILS)

# Utilities which are scripts
SCRIPT_UTILS = $(SCRIPT_UTILS_all) $(SCRIPT_UTILS_py2)
# For all CASINO_ARCHs (bash scripts)
SCRIPT_UTILS_all = \
$(BINDIR)/billy \
$(BINDIR)/casinohelp \
$(BINDIR)/clearup \
$(BINDIR)/clearup_twistav \
$(BINDIR)/crysgen98 \
$(BINDIR)/rna \
$(BINDIR)/extr_pwscf \
$(BINDIR)/extr_casino \
$(BINDIR)/gamess2qmc \
$(BINDIR)/getdet \
$(BINDIR)/get_exciton_binding \
$(BINDIR)/get_qp_gap \
$(BINDIR)/graphdmc \
$(BINDIR)/graphit \
$(BINDIR)/killvp \
$(BINDIR)/louishelp \
$(BINDIR)/louisplot2d \
$(BINDIR)/louisplot3d \
$(BINDIR)/menugrep \
$(BINDIR)/modify_inputs \
$(BINDIR)/modify_louis_inputs \
$(BINDIR)/movie2avi \
$(BINDIR)/mpr \
$(BINDIR)/multirun \
$(BINDIR)/opt_crystal \
$(BINDIR)/plot_2D \
$(BINDIR)/plot_bffield \
$(BINDIR)/plot_bfphi \
$(BINDIR)/plot_louis \
$(BINDIR)/plot_reblock \
$(BINDIR)/runcrystal \
$(BINDIR)/runlouis \
$(BINDIR)/runpwscf \
$(BINDIR)/runqmc \
$(BINDIR)/runqmcmd \
$(BINDIR)/runvp \
$(BINDIR)/syncvp \
$(BINDIR)/tailvp \
$(BINDIR)/tcm_comps \
$(BINDIR)/twistanalysis_castep \
$(BINDIR)/twistav_castep \
$(BINDIR)/twistav_pwscf \
$(BINDIR)/update_input \
$(BINDIR)/update_src \
$(BINDIR)/ve \
$(BINDIR)/envmc \
$(BINDIR)/vevp \
$(BINDIR)/envmcvp
# For some systems only
SCRIPT_UTILS_py2 = \
$(BINDIR)/adf2stowf \
$(BINDIR)/molden2qmc

# Dependencies
$(BINDIR)/adf2stowf: wfn_converters/adf/adf2stowf.py
$(BINDIR)/billy: billy/billy
$(BINDIR)/casinohelp: help/casinohelp
$(BINDIR)/clearup: clearup/clearup
$(BINDIR)/clearup_twistav: twist/clearup_twistav
$(BINDIR)/crysgen98: wfn_converters/crystal9x/crysgen98
$(BINDIR)/rna : rna/rna
$(BINDIR)/extr_pwscf : runqmcmd/extr_pwscf
$(BINDIR)/extr_casino : runqmcmd/extr_casino
$(BINDIR)/gamess2qmc : wfn_converters/gamess/gamess2qmc
$(BINDIR)/getdet : wfn_converters/gamess/getdet
$(BINDIR)/get_exciton_binding : get_gaps/get_exciton_binding
$(BINDIR)/get_qp_gap : get_gaps/get_qp_gap
$(BINDIR)/graphdmc: graphdmc/graphdmc
$(BINDIR)/graphit: graphdmc/graphit
$(BINDIR)/killvp : runvp/killvp
$(BINDIR)/louishelp: louis/utils/help/louishelp
$(BINDIR)/louisplot2d: louis/utils/plot_other/louisplot2d
$(BINDIR)/louisplot3d: louis/utils/plot_other/louisplot3d
$(BINDIR)/menugrep : menugrep/menugrep
$(BINDIR)/modify_inputs: modify_inputs/modify_inputs
$(BINDIR)/modify_louis_inputs: louis/utils/change_inputs/modify_louis_inputs
$(BINDIR)/molden2qmc: wfn_converters/molden/molden2qmc.py
$(BINDIR)/movie2avi: movie2avi/movie2avi
$(BINDIR)/mpr: pretty_printer/mpr
$(BINDIR)/multirun: multirun/multirun
$(BINDIR)/opt_crystal: opt_crystal/opt_crystal
$(BINDIR)/plot_2D: plot_2D/plot_2D
$(BINDIR)/plot_bffield: plot_backflow/plot_bffield
$(BINDIR)/plot_bfphi: plot_backflow/plot_bfphi
$(BINDIR)/plot_louis: louis/utils/plot/plot_louis
$(BINDIR)/plot_reblock: reblock/plot_reblock
$(BINDIR)/runcrystal: wfn_converters/crystal14/crystal_run_script/runcrystal
$(BINDIR)/runlouis: louis/utils/runlouis/runlouis
$(BINDIR)/runpwscf: wfn_converters/pwscf/runpwscf
$(BINDIR)/runqmc: runqmc/runqmc
$(BINDIR)/runqmcmd: runqmcmd/runqmcmd
$(BINDIR)/runvp : runvp/runvp
$(BINDIR)/syncvp : runvp/syncvp
$(BINDIR)/tailvp : runvp/tailvp
$(BINDIR)/tcm_comps: tcm_comps/tcm_comps
$(BINDIR)/twistanalysis_castep: twist/CASTEP/twistanalysis_castep
$(BINDIR)/twistav_castep: twist/CASTEP/twistav_castep
$(BINDIR)/twistav_pwscf: twist/PWSCF/twistav_pwscf
$(BINDIR)/update_input: converters_2.0/update_input/update_input
$(BINDIR)/update_src: update_src/update_src
$(BINDIR)/ve: envmc/ve
$(BINDIR)/envmc: envmc/envmc
$(BINDIR)/vevp : runvp/vevp
$(BINDIR)/envmcvp : runvp/envmcvp

# Shortcut for scripts group and individual utils
.PHONY: scripts
scripts: $(SCRIPT_UTILS)
.PHONY: adf2stowf
adf2stowf: $(BINDIR)/adf2stowf
.PHONY: billy
billy: $(BINDIR)/billy
.PHONY: casinohelp
casinohelp: $(BINDIR)/casinohelp
.PHONY: clearup
clearup: $(BINDIR)/clearup
.PHONY: clearup_twistav
clearup_twistav: $(BINDIR)/clearup_twistav
.PHONY: crysgen98
crysgen98: $(BINDIR)/crysgen98
.PHONY: rna
rna : $(BINDIR)/rna
.PHONY: gamess2qmc
gamess2qmc : $(BINDIR)/gamess2qmc
.PHONY: getdet
getdet : $(BINDIR)/getdet
.PHONY: get_exciton_binding
get_exciton_binding : $(BINDIR)/get_exciton_binding
.PHONY: get_qp_gap
get_qp_gap : $(BINDIR)/get_qp_gap
.PHONY: graphdmc
graphdmc: $(BINDIR)/graphdmc
.PHONY: graphit
graphit: $(BINDIR)/graphit
.PHONY: killvp
killvp : $(BINDIR)/killvp
.PHONY: louishelp
louishelp: $(BINDIR)/louishelp
.PHONY: louisplot2d
louisplot2d: $(BINDIR)/louisplot2d
.PHONY: louisplot3d
louisplot3d: $(BINDIR)/louisplot3d
.PHONY: menugrep
menugrep : $(BINDIR)/menugrep
.PHONY: modify_inputs
modify_inputs: $(BINDIR)/modify_inputs
.PHONY: modify_louis_inputs
modify_louis_inputs: $(BINDIR)/modify_louis_inputs
.PHONY: molden2qmc
molden2qmc: $(BINDIR)/molden2qmc
.PHONY: movie2avi
movie2avi: $(BINDIR)/movie2avi
.PHONY: mpr
mpr: $(BINDIR)/mpr
.PHONY: multirun
multirun: $(BINDIR)/multirun
.PHONY: opt_crystal
opt_crystal: $(BINDIR)/opt_crystal
.PHONY: plot_2D
plot_2D: $(BINDIR)/plot_2D
.PHONY: plot_bffield
plot_bffield: $(BINDIR)/plot_bffield
.PHONY: plot_bfphi
plot_bfphi: $(BINDIR)/plot_bfphi
.PHONY: plot_louis
plot_louis: $(BINDIR)/plot_louis
.PHONY: plot_reblock
plot_reblock: $(BINDIR)/plot_reblock
.PHONY: runcrystal
runcrystal: $(BINDIR)/runcrystal
.PHONY: runlouis
runlouis: $(BINDIR)/runlouis
.PHONY: runpwscf
runpwscf: $(BINDIR)/runpwscf
.PHONY: runqmc
runqmc: $(BINDIR)/runqmc
.PHONY: runqmcmd
runqmcmd: $(BINDIR)/runqmcmd
.PHONY: runvp
runvp : $(BINDIR)/runvp
.PHONY: syncvp
syncvp : $(BINDIR)/syncvp
.PHONY: tailvp
tailvp : $(BINDIR)/tailvp
.PHONY: tcm_comps
tcm_comps: $(BINDIR)/tcm_comps
.PHONY: twistanalysis_castep
twistanalysis: $(BINDIR)/twistanalysis_castep
.PHONY: twistav_castep
twistav: $(BINDIR)/twistav_castep
.PHONY: twistav_pwscf
twistav: $(BINDIR)/twistav_pwscf
.PHONY: update_input
update_input: $(BINDIR)/update_input
.PHONY: update_src
update_src: $(BINDIR)/update_src
.PHONY: ve
ve: $(BINDIR)/ve
.PHONY: envmc
envmc: $(BINDIR)/envmc
.PHONY: vevp
vevp : $(BINDIR)/vevp
.PHONY: envmcvp
envmcvp : $(BINDIR)/envmcvp

# How to make them
LINKSCRIPT = { $(MKBINDIR) ;\
           bash -c '$(BASH_PRETTY) ; compile_info LN utils $(@F)' ;\
           rm -f $@ 2> /dev/null ; ln -s $(UTILDIR)/$? $@ ; }
$(SCRIPT_UTILS_all):
	@$(LINKSCRIPT)
$(SCRIPT_UTILS_py2):
	@[ ! -z "`bash -c 'type -P python2 2> /dev/null'`" ]\
        && $(LINKSCRIPT) || $(MKDUMMY_nopath)


# Utilities with an independent Makefile. NB, the dependencies are
# directly handled by the relevant Makefiles.
SUBMAKE_UTILS = $(SUBMAKE_UTILS_all) $(SUBMAKE_UTILS_c) $(SUBMAKE_UTILS_cxx)
# For all CASINO_ARCHs
SUBMAKE_UTILS_all = \
 $(BINDIR)/blip \
 $(BINDIR)/casinohelp_helper \
 $(BINDIR)/crysgen06 \
 $(BINDIR)/crysgen09 \
 $(BINDIR)/crysgen14 \
 $(BINDIR)/dfit \
 $(BINDIR)/extrapolate_tau \
 $(BINDIR)/finsize \
 $(BINDIR)/gaussiantoqmc \
 $(BINDIR_FEATURE)/louis \
 $(BINDIR)/louishelp_helper \
 $(BINDIR)/localizer \
 $(BINDIR)/mcta_post_process \
 $(BINDIR)/plot_expval \
 $(BINDIR)/plot_mpc
# $(BINDIR)/mcta_hf
# For some CASINO_ARCHs only
SUBMAKE_UTILS_c = \
 $(BINDIR)/det_compress \
 $(BINDIR)/jeep_to_pwfn
SUBMAKE_UTILS_cxx = \
 $(BINDIR)/ppconvert

# Directory and makefile list
MAKEDIR_blip              = wfn_converters/pw_to_blips
MAKEDIR_casinohelp_helper = help
MAKEDIR_crysgen06         = wfn_converters/crystal06
MAKEDIR_crysgen09         = wfn_converters/crystal09
MAKEDIR_crysgen14         = wfn_converters/crystal14
MAKEDIR_det_compress      = det_compress
MAKEDIR_dfit              = billy/fortran/dfit
MAKEDIR_extrapolate_tau   = extrapolate_tau
MAKEDIR_finsize           = finsize
MAKEDIR_gaussiantoqmc     = wfn_converters/gaussian9x-03-09
MAKEDIR_jeep_to_pwfn      = wfn_converters/gp-jeep-qbox
MAKEDIR_localizer         = wfn_converters/pw_to_blips/localizer
MAKEDIR_louis             = louis/src
MAKEDIR_louishelp_helper  = louis/utils/help
MAKEDIR_mcta_hf           = heg_mcta/mcta_hf
MAKEDIR_mcta_post_process = heg_mcta/mcta_post_process
MAKEDIR_plot_expval       = plot_expval
MAKEDIR_plot_mpc          = plot_mpc
MAKEDIR_ppconvert = pseudo_converters/abinit/casino_to_abinit/ppconvert_esler

# Dependencies
$(BINDIR)/blip: \
 wfn_converters/pw_to_blips/blip.f90 \
 wfn_converters/pw_to_blips/blipp.f90 \
 wfn_converters/pw_to_blips/blipk.f90 \
 wfn_converters/pw_to_blips/helpers.f90 \
 wfn_converters/pw_to_blips/rng.f90 \
 wfn_converters/pw_to_blips/singleton.f90
$(BINDIR)/casinohelp_helper: \
 help/casinohelp_helper.f90 \
 help/dsp.f90 \
 help/esdf.f90 \
 help/esdf_key.f90 \
 help/format_utils.f90 \
 help/run_control.f90 \
 help/store.f90
$(BINDIR)/crysgen06: \
 wfn_converters/crystal06/casino_interface.f90 \
 wfn_converters/crystal06/condft.f90 \
 wfn_converters/crystal06/convrs.f90 \
 wfn_converters/crystal06/cryapi_inp_c06.f90 \
 wfn_converters/crystal06/crysgen06.f90 \
 wfn_converters/crystal06/estros.f90 \
 wfn_converters/crystal06/expo_module.f90 \
 wfn_converters/crystal06/expt.f90 \
 wfn_converters/crystal06/expu.f90 \
 wfn_converters/crystal06/lmaxxx.f90 \
 wfn_converters/crystal06/matout.f90 \
 wfn_converters/crystal06/mcm.f90 \
 wfn_converters/crystal06/memory_screen.f90 \
 wfn_converters/crystal06/paral1_module.f90 \
 wfn_converters/crystal06/rotation_matrices.f90 \
 wfn_converters/crystal06/rotmatrix.f90 \
 wfn_converters/crystal06/shell_info.f90 \
 wfn_converters/crystal06/sphfac_module.f90 \
 wfn_converters/crystal06/vrslat.f90
$(BINDIR)/crysgen09: \
 wfn_converters/crystal09/casino_interface.f90 \
 wfn_converters/crystal09/condft.f90 \
 wfn_converters/crystal09/convrs.f90 \
 wfn_converters/crystal09/cryapi_inp_c09.f90 \
 wfn_converters/crystal09/crysgen09.f90 \
 wfn_converters/crystal09/estros.f90 \
 wfn_converters/crystal09/expo_module.f90 \
 wfn_converters/crystal09/expt.f90 \
 wfn_converters/crystal09/expu.f90 \
 wfn_converters/crystal09/lmaxxx.f90 \
 wfn_converters/crystal09/matout.f90 \
 wfn_converters/crystal09/mcm.f90 \
 wfn_converters/crystal09/memory_screen.f90 \
 wfn_converters/crystal09/paral1_module.f90 \
 wfn_converters/crystal09/rotation_matrices.f90 \
 wfn_converters/crystal09/rotmatrix.f90 \
 wfn_converters/crystal09/shell_info.f90 \
 wfn_converters/crystal09/sphfac_module.f90 \
 wfn_converters/crystal09/vrslat.f90
$(BINDIR)/crysgen14: \
 wfn_converters/crystal14/casino_interface.f90 \
 wfn_converters/crystal14/condft.f90 \
 wfn_converters/crystal14/convrs.f90 \
 wfn_converters/crystal14/cryapi_inp_c09.f90 \
 wfn_converters/crystal14/crysgen14.f90 \
 wfn_converters/crystal14/estros.f90 \
 wfn_converters/crystal14/expo_module.f90 \
 wfn_converters/crystal14/expt.f90 \
 wfn_converters/crystal14/expu.f90 \
 wfn_converters/crystal14/lmaxxx.f90 \
 wfn_converters/crystal14/matout.f90 \
 wfn_converters/crystal14/mcm.f90 \
 wfn_converters/crystal14/memory_screen.f90 \
 wfn_converters/crystal14/paral1_module.f90 \
 wfn_converters/crystal14/rotation_matrices.f90 \
 wfn_converters/crystal14/rotmatrix.f90 \
 wfn_converters/crystal14/shell_info.f90 \
 wfn_converters/crystal14/sphfac_module.f90 \
 wfn_converters/crystal14/vrslat.f90
$(BINDIR)/det_compress: \
 det_compress/casl.f90 \
 det_compress/det_compress.f90 \
 det_compress/lp2fortran.c \
 det_compress/tools.f90
$(BINDIR)/dfit: \
 billy/fortran/dfit/dfit.f90 \
 billy/fortran/dfit/dsp.f90 \
 billy/fortran/dfit/errvrs.f90 \
 billy/fortran/dfit/lesolv.f90 \
 billy/fortran/dfit/mxmb.f90
$(BINDIR)/extrapolate_tau: \
 extrapolate_tau/extrapolate_tau.f90
$(BINDIR)/finsize: \
 finsize/dsp.f90 \
 finsize/esdf.f90 \
 finsize/esdf_key.f90 \
 finsize/finsize.f90 \
 finsize/format_utils.f90 \
 finsize/main.f90 \
 finsize/nl2sol.f90 \
 finsize/quadrature.f90 \
 finsize/run_control.f90 \
 finsize/store.f90
$(BINDIR)/gaussiantoqmc: \
 wfn_converters/gaussian9x-03-09/analyze_cis_state.f90 \
 wfn_converters/gaussian9x-03-09/awk_like.f90 \
 wfn_converters/gaussian9x-03-09/cas_wfn.f90 \
 wfn_converters/gaussian9x-03-09/cas_write.f90 \
 wfn_converters/gaussian9x-03-09/cis_data.f90 \
 wfn_converters/gaussian9x-03-09/con_coeffs.f90 \
 wfn_converters/gaussian9x-03-09/fatal.f90 \
 wfn_converters/gaussian9x-03-09/g94_wavefunction.f90 \
 wfn_converters/gaussian9x-03-09/g_d_type.f90 \
 wfn_converters/gaussian9x-03-09/g_s_type.f90 \
 wfn_converters/gaussian9x-03-09/gaussiantoqmc.f90 \
 wfn_converters/gaussian9x-03-09/get_gauss_version.f90 \
 wfn_converters/gaussian9x-03-09/integ_params.f90 \
 wfn_converters/gaussian9x-03-09/max_coincidence.f90 \
 wfn_converters/gaussian9x-03-09/normalization_check.f90 \
 wfn_converters/gaussian9x-03-09/normalize_ci.f90 \
 wfn_converters/gaussian9x-03-09/numsrt_2way.f90 \
 wfn_converters/gaussian9x-03-09/pack_evcoeffs.f90 \
 wfn_converters/gaussian9x-03-09/paramfile.f90 \
 wfn_converters/gaussian9x-03-09/psi.f90 \
 wfn_converters/gaussian9x-03-09/qmc_write.f90 \
 wfn_converters/gaussian9x-03-09/re_sum.f90 \
 wfn_converters/gaussian9x-03-09/read_G9xout.f90 \
 wfn_converters/gaussian9x-03-09/read_fchk.f90 \
 wfn_converters/gaussian9x-03-09/rejig.f90 \
 wfn_converters/gaussian9x-03-09/resum_cas.f90 \
 wfn_converters/gaussian9x-03-09/shell_centres.f90 \
 wfn_converters/gaussian9x-03-09/sum_degen_excite.f90 \
 wfn_converters/gaussian9x-03-09/user_control.f90 \
 wfn_converters/gaussian9x-03-09/wfn_construct.f90 \
 wfn_converters/gaussian9x-03-09/wfn_test.f90
$(BINDIR)/jeep_to_pwfn: \
 wfn_converters/gp-jeep-qbox/gen_gvec.f90 \
 wfn_converters/gp-jeep-qbox/jeep_to_pwfn.f90 \
 wfn_converters/gp-jeep-qbox/jeep_io.c
$(BINDIR)/localizer: \
 wfn_converters/pw_to_blips/localizer/localizer.f90 \
 wfn_converters/pw_to_blips/localizer/singleton.f90
$(BINDIR)/louis: \
 louis/src/bulirsch_stoer.f90 \
 louis/src/choose_wfn.f90 \
 louis/src/comms_serial.f90 \
 louis/src/comms_parallel_mpi1.f90 \
 louis/src/comms_parallel_mpi2.f90 \
 louis/src/dsp.f90 \
 louis/src/esdf.f90 \
 louis/src/esdf_key.f90 \
 louis/src/etime.c \
 louis/src/eval_density.f90 \
 louis/src/file_utils.f90 \
 louis/src/format_utils.f90 \
 louis/src/input.f90 \
 louis/src/louis.f90 \
 louis/src/numerical.f90 \
 louis/src/parallel.f90 \
 louis/src/random_numbers.f90 \
 louis/src/run_control.f90 \
 louis/src/runge_kutta.f90 \
 louis/src/scaled_sin_wfn.f90 \
 louis/src/sin_wfn.f90 \
 louis/src/store.f90
$(BINDIR)/louishelp_helper: \
 louis/utils/help/louishelp_helper.f90 \
 louis/utils/help/dsp.f90 \
 louis/utils/help/esdf.f90 \
 louis/utils/help/esdf_key.f90 \
 louis/utils/help/run_control.f90 \
 louis/utils/help/store.f90
$(BINDIR)/mcta_hf: \
 heg_mcta/mcta_hf/heg_tools.f90 \
 heg_mcta/mcta_hf/mcta_hf.f90
$(BINDIR)/mcta_post_process: \
 heg_mcta/mcta_post_process/heg_tools.f90 \
 heg_mcta/mcta_post_process/mcta_post_process.f90
$(BINDIR)/plot_expval: \
 plot_expval/dsp.f90 \
 plot_expval/esdf.f90 \
 plot_expval/esdf_key.f90 \
 plot_expval/format_utils.f90 \
 plot_expval/plot_expval.f90 \
 plot_expval/run_control.f90 \
 plot_expval/store.f90
$(BINDIR)/plot_mpc: \
 plot_mpc/dsp.f90 \
 plot_mpc/esdf.f90 \
 plot_mpc/esdf_key.f90 \
 plot_mpc/format_utils.f90 \
 plot_mpc/plot_mpc.f90 \
 plot_mpc/run_control.f90 \
 plot_mpc/store.f90
$(BINDIR)/ppconvert: \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/CubicSpline.cc \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/CubicSpline.h \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/GeneralGrid.h \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/NLPPClass.cc \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/NLPPClass.h \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/ParseCommand.cc \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/ParseCommand.h \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/ParserClass.cc \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/ParserClass.h \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/XMLWriterClass2.cc \
 pseudo_converters/abinit/casino_to_abinit/ppconvert_esler/XMLWriterClass2.h

# Shortcut for submake group and individual utils
.PHONY: submake
submake: $(SUBMAKE_UTILS)
.PHONY: blip
blip: $(BINDIR)/blip
.PHONY: casinohelp_helper
casinohelp_helper: $(BINDIR)/casinohelp_helper
.PHONY: crysgen06
crysgen06: $(BINDIR)/crysgen06
.PHONY: crysgen09
crysgen09: $(BINDIR)/crysgen09
.PHONY: crysgen14
crysgen14: $(BINDIR)/crysgen14
.PHONY: det_compress
det_compress: $(BINDIR)/det_compress
.PHONY: dfit
dfit: $(BINDIR)/dfit
.PHONY: extrapolate_tau
extrapolate_tau: $(BINDIR)/extrapolate_tau
.PHONY: finsize
finsize: $(BINDIR)/finsize
.PHONY: gaussiantoqmc
gaussiantoqmc: $(BINDIR)/gaussiantoqmc
.PHONY: jeep_to_pwfn
jeep_to_pwfn: $(BINDIR)/jeep_to_pwfn
.PHONY: localizer
localizer: $(BINDIR)/localizer
.PHONY: louis
louis: $(BINDIR)/louis
.PHONY: louishelp_helper
louishelp_helper: $(BINDIR)/louishelp_helper
.PHONY: mcta_hf
mcta_hf: $(BINDIR)/mcta_hf
.PHONY: mcta_post_process
mcta_post_process: $(BINDIR)/mcta_post_process
.PHONY: plot_expval
plot_expval: $(BINDIR)/plot_expval
.PHONY: plot_mpc
plot_mpc: $(BINDIR)/plot_mpc
.PHONY: ppconvert
ppconvert: $(BINDIR)/ppconvert

# How to make them
# NB, GNU make requires '$(MAKE)' to be explicitly specified in the
# rule in order for the jobserver to be available to sub-makes. So the
# following statement:
#  { [ -L "$@" ] && rm -f "$@" ; $(MKBINDIR) ; cd $(MAKEDIR_$(@F))\
#    && $(MAKE) 'PWD=$(UTILDIR)/$(MAKEDIR_$(@F))' || exit 1 ; }
# is split around $(MAKE):
SUBMAKE1 = { [ -L "$@" ] && rm -f "$@" ; $(MKBINDIR) ; cd $(MAKEDIR_$(@F)) &&
SUBMAKE2 = 'PWD=$(UTILDIR)/$(MAKEDIR_$(@F))' || exit 1 ; }
$(SUBMAKE_UTILS_all):
	@$(SUBMAKE1) $(MAKE) $(SUBMAKE2)
$(SUBMAKE_UTILS_c):
	@[ ! -z "$(CC_UTILS)" ] && $(SUBMAKE1) $(MAKE) $(SUBMAKE2)\
        || $(MKDUMMY_nopath)
$(SUBMAKE_UTILS_cxx):
	@[ ! -z "$(CXX_UTILS)" ] && $(SUBMAKE1) $(MAKE) $(SUBMAKE2)\
        || $(MKDUMMY_nopath)

# Utilities that only require compiling a single F90 file
SINGLE_F90_UTILS = $(SINGLE_F90_UTILS_all)
# For all CASINO_ARCHs
SINGLE_F90_UTILS_all = \
$(BINDIR)/abinit_to_casino_pp \
$(BINDIR)/champ_to_casino_pp \
$(BINDIR)/casino2gon \
$(BINDIR)/combine_cdata \
$(BINDIR)/crystaltoqmc \
$(BINDIR)/egaussian \
$(BINDIR)/envmc_helper \
$(BINDIR)/extrapolated_est \
$(BINDIR)/expval_average \
$(BINDIR)/findmax \
$(BINDIR)/format_configs \
$(BINDIR)/helpbilly \
$(BINDIR)/helpbilly2 \
$(BINDIR)/input_kw_conv \
$(BINDIR)/ion_dist \
$(BINDIR)/jeep_pp \
$(BINDIR)/make_new_mpc \
$(BINDIR)/make_p_stars \
$(BINDIR)/nstring \
$(BINDIR)/percent_calculator \
$(BINDIR)/plot_hist \
$(BINDIR)/ptm \
$(BINDIR)/quad_fit \
$(BINDIR)/quickblock \
$(BINDIR)/reblock \
$(BINDIR)/supercell \
$(BINDIR)/tahelper \
$(BINDIR)/twistoffset_castep \
$(BINDIR)/twistoffset_pwscf \
$(BINDIR)/update_config \
$(BINDIR)/update_hist

# Dependencies
$(BINDIR)/abinit_to_casino_pp: pseudo_converters/abinit/abinit_to_casino/abinit_to_casino_pp.f90
$(BINDIR)/champ_to_casino_pp: pseudo_converters/champ/champ_to_casino_pp.f90
$(BINDIR)/casino2gon: pseudo_converters/pwscf/casino2gon/casino2gon.f90
$(BINDIR)/combine_cdata: combine/combine_cdata.f90
$(BINDIR)/crystaltoqmc: wfn_converters/crystal9x/crystaltoqmc.f90
$(BINDIR)/egaussian: wfn_converters/gaussian9x-03-09/UTILITIES/egaussian.f90
$(BINDIR)/envmc_helper: envmc/envmc_helper.f90
$(BINDIR)/expval_average: plot_expval/expval_average/expval_average.f90
$(BINDIR)/extrapolated_est: extrapolated_est/extrapolated_est.f90
$(BINDIR)/findmax: plot_2D/findmax.f90
$(BINDIR)/format_configs: format_configs/format_configs.f90
$(BINDIR)/helpbilly: billy/fortran/helpbilly.f90
$(BINDIR)/helpbilly2: billy/fortran/helpbilly2.f90
$(BINDIR)/input_kw_conv: input_kw_conv/input_kw_conv.f90
$(BINDIR)/ion_dist: ion_dist/ion_dist.f90
$(BINDIR)/jeep_pp: pseudo_converters/gp-jeep-qbox/jeep_pp.f90
$(BINDIR)/make_p_stars: make_p_stars/make_p_stars.f90
$(BINDIR)/make_new_mpc: converters_2.0/make_new_mpc/make_new_mpc.f90
$(BINDIR)/nstring: nstring/nstring.f90
$(BINDIR)/percent_calculator: percent_calculator/percent_calculator.f90
$(BINDIR)/plot_hist: plot_hist/plot_hist.f90
$(BINDIR)/ptm: ptm/ptm.f90
$(BINDIR)/quad_fit: quad_fit/quad_fit.f90
$(BINDIR)/quickblock: reblock/quickblock.f90
$(BINDIR)/reblock: reblock/reblock.f90
$(BINDIR)/supercell: supercell/supercell.f90
$(BINDIR)/tahelper: twist/tahelper.f90
$(BINDIR)/twistoffset_castep: twist/CASTEP/twistoffset_castep.f90
$(BINDIR)/twistoffset_pwscf: twist/PWSCF/twistoffset_pwscf.f90
$(BINDIR)/update_config: converters_2.0/update_config/update_config.f90
$(BINDIR)/update_hist: converters_2.0/update_hist/update_hist.f90

# Shortcut for f90 group and individual utils
.PHONY: f90
f90: $(SINGLE_F90_UTILS_all)
.PHONY: abinit_to_casino_pp
abinit_to_casino_pp: $(BINDIR)/abinit_to_casino_pp
.PHONY: champ_to_casino_pp
champ_to_casino_pp: $(BINDIR)/champ_to_casino_pp
.PHONY: casino2gon
casino2gon: $(BINDIR)/casino2gon
.PHONY: combine_cdata
combine_cdata: $(BINDIR)/combine_cdata
.PHONY: crystaltoqmc
crystaltoqmc: $(BINDIR)/crystaltoqmc
.PHONY: egaussian
egaussian: $(BINDIR)/egaussian
.PHONY: envmc_helper
ve_helper: $(BINDIR)/envmc_helper
.PHONY: extrapolated_est
extrapolated_est: $(BINDIR)/extrapolated_est
.PHONY: expval_average
expval_average: $(BINDIR)/expval_average
.PHONY: findmax
findmax: $(BINDIR)/findmax
.PHONY: format_configs
format_configs: $(BINDIR)/format_configs
.PHONY: helpbilly
helpbilly: $(BINDIR)/helpbilly
.PHONY: helpbilly2
helpbilly2: $(BINDIR)/helpbilly2
.PHONY: input_kw_conv
input_kw_conv: $(BINDIR)/input_kw_conv
.PHONY: ion_dist
ion_dist: $(BINDIR)/ion_dist
.PHONY: jeep_pp
jeep_pp: $(BINDIR)/jeep_pp
.PHONY: make_new_mpc
make_new_mpc: $(BINDIR)/make_new_mpc
.PHONY: make_p_stars
make_p_stars: $(BINDIR)/make_p_stars
.PHONY: nstring
nstring: $(BINDIR)/nstring
.PHONY: percent_calculator
percent_calculator: $(BINDIR)/percent_calculator
.PHONY: plot_hist
plot_hist: $(BINDIR)/plot_hist
.PHONY: ptm
ptm: $(BINDIR)/ptm
.PHONY: quad_fit
quad_fit: $(BINDIR)/quad_fit
.PHONY: quickblock
quickblock: $(BINDIR)/quickblock
.PHONY: reblock
reblock: $(BINDIR)/reblock
.PHONY: supercell
supercell: $(BINDIR)/supercell
.PHONY: tahelper
tahelper: $(BINDIR)/tahelper
.PHONY: twistoffset_castep
twistoffset: $(BINDIR)/twistoffset_castep
.PHONY: twistoffset_pwscf
twistoffset: $(BINDIR)/twistoffset_pwscf
.PHONY: update_config
update_config: $(BINDIR)/update_config
.PHONY: update_hist
update_hist: $(BINDIR)/update_hist

# How to make them
F90_SINGLE = {\
           bash -c '$(BASH_PRETTY) ; compile_info F90 utils $(@F)' ;\
           rm -f $@ 2> /dev/null ;\
           [ -d $(OBJDIR)/$(@F) ] || mkdir -p $(OBJDIR)/$(@F) ;\
           $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)/$(@F)\
           && $(F90_UTILS) $(FFLAGS_UTILS) -o $@ $(UTILDIR)/$? || exit 1 ; }
$(SINGLE_F90_UTILS_all):
	@$(F90_SINGLE)

# Target for creating/updating wrapper scripts
.PHONY: wrap utils

utils: $(UTILS)

wrap: info utils
	@bash -c '$(BASH_PRETTY) ; compile_info WRAP utils' ;\
        $(MKBINDIR) ; cd $(BINBASE) && $(TOOLDIR)/make-wrappers\
        --environment-command='$(ENVIRONMENT_COMMAND_UTILS)'

# Target for deleting object files
.PHONY: clean
clean:
	@echo "Cleaning object files [$(CASINO_ARCH)]." ;\
        rm -rf $(OBJDIR) ; :

# Target for deleting object files and compiled utils
.PHONY: vclean
vclean:
	@echo "Cleaning object files and compiled utilities [$(CASINO_ARCH)]."\
        ; rm -rf $(OBJDIR) $(BINDIR) ; :

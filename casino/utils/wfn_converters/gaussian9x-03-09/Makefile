#----------------------------------------------#
# General Makefile for gaussiantoqmc program.  #
#                                              #
# MDT 4.2003                                   #
#----------------------------------------------#
# Basic stuff
UTILEXE = gaussiantoqmc
SRCDIR  = $(PWD)
UTILDIR = $(SRCDIR)/../..
TOOLDIR = $(UTILDIR)/../build-tools
INCBASE = $(UTILDIR)/../arch/data
BINDIR  = $(UTILDIR)/../bin_qmc/utils/$(CASINO_ARCH)
OBJDIR  = $(UTILDIR)/zlib/$(CASINO_ARCH)/$(UTILEXE)
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

start: default
.PHONY: start

# System-independent variables
UTILS_MODE     =
UNATIVE_       =
UNATIVE_native = _NATIVE

# Defaults for system-dependent variables
LDF90                      = $(F90)
LDFLAGS_opt                = $(FFLAGS_opt)
LDFLAGS_all                = $(FFLAGS_all)
ENVIRONMENT_COMMAND        = :
LDF90_NATIVE               = $(F90_NATIVE)
LDFLAGS_opt_NATIVE         = $(FFLAGS_opt_NATIVE)
LDFLAGS_all_NATIVE         = $(FFLAGS_all_NATIVE)
ENVIRONMENT_COMMAND_NATIVE = :

# System dependent include files
include $(INCBASE)/$(CASINO_ARCH).arch

# Pick which set of variables to use
UNATIVE                   = $(UNATIVE_$(UTILS_MODE))
ENVIRONMENT_COMMAND_UTILS = $(ENVIRONMENT_COMMAND$(UNATIVE))
F90_UTILS                 = $(F90$(UNATIVE))
FFLAGS_opt_UTILS          = $(FFLAGS_opt$(UNATIVE))
FFLAGS_all_UTILS          = $(FFLAGS_all$(UNATIVE))
INCLUDE_DIR_UTILS         = $(INCLUDE_DIR$(UNATIVE))
LDF90_UTILS               = $(LDF90$(UNATIVE))
LDFLAGS_opt_UTILS         = $(LDFLAGS_opt$(UNATIVE))
LDFLAGS_all_UTILS         = $(LDFLAGS_all$(UNATIVE))
LIB_PATH_UTILS            = $(LIB_PATH$(UNATIVE))
LDLIBS_opt_UTILS          = $(LDLIBS_opt$(UNATIVE))
LDLIBS_all_UTILS          = $(LDLIBS_all$(UNATIVE))

# System-dependent variables from data in include files
FFLAGS_UTILS  = $(FFLAGS_opt_UTILS) $(FFLAGS_all_UTILS) $(INCLUDE_DIR_UTILS)
LDFLAGS_UTILS = $(LDFLAGS_opt_UTILS) $(LDFLAGS_all_UTILS) $(LIB_PATH_UTILS)
LDLIBS_UTILS  = $(LDLIBS_opt_UTILS) $(LDLIBS_all_UTILS)

# Main target
.PHONY: default
default: $(BINDIR)/$(UTILEXE)

# 'Clean' target
.PHONY: clean
clean:
	@rm -rf $(OBJDIR) $(BINDIR)/$(UTILEXE)

# Object-file list
OBJ = \
$(OBJDIR)/analyze_cis_state.o \
$(OBJDIR)/awk_like.o \
$(OBJDIR)/cas_wfn.o \
$(OBJDIR)/cas_write.o \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/con_coeffs.o \
$(OBJDIR)/fatal.o \
$(OBJDIR)/g94_wavefunction.o \
$(OBJDIR)/g_d_type.o \
$(OBJDIR)/g_s_type.o \
$(OBJDIR)/gaussiantoqmc.o \
$(OBJDIR)/get_gauss_version.o \
$(OBJDIR)/integ_params.o \
$(OBJDIR)/max_coincidence.o \
$(OBJDIR)/normalization_check.o \
$(OBJDIR)/normalize_ci.o \
$(OBJDIR)/numsrt_2way.o \
$(OBJDIR)/pack_evcoeffs.o \
$(OBJDIR)/paramfile.o \
$(OBJDIR)/psi.o \
$(OBJDIR)/qmc_write.o \
$(OBJDIR)/re_sum.o \
$(OBJDIR)/read_G9xout.o \
$(OBJDIR)/read_fchk.o \
$(OBJDIR)/rejig.o \
$(OBJDIR)/resum_cas.o \
$(OBJDIR)/shell_centres.o \
$(OBJDIR)/sum_degen_excite.o \
$(OBJDIR)/user_control.o \
$(OBJDIR)/wfn_construct.o \
$(OBJDIR)/wfn_test.o
# End object-file list

# Module dependencies
ANALYZE_CIS_STATE_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/analyze_cis_state.o : $(SRCDIR)/analyze_cis_state.f90 $(ANALYZE_CIS_STATE_DEP)
analyze_cis_state.o_name = analyze_cis_state

AWK_LIKE_DEP = \
$(OBJDIR)/paramfile.o
$(OBJDIR)/awk_like.o : $(SRCDIR)/awk_like.f90 $(AWK_LIKE_DEP)
awk_like.o_name = awk_like

CAS_WFN_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/cas_wfn.o : $(SRCDIR)/cas_wfn.f90 $(CAS_WFN_DEP)
cas_wfn.o_name = cas_wfn

CAS_WRITE_DEP = \
$(OBJDIR)/cis_data.o
$(OBJDIR)/cas_write.o : $(SRCDIR)/cas_write.f90 $(CAS_WRITE_DEP)
cas_write.o_name = cas_write

CIS_DATA_DEP = \
$(OBJDIR)/paramfile.o
$(OBJDIR)/cis_data.o : $(SRCDIR)/cis_data.f90 $(CIS_DATA_DEP)
cis_data.o_name = cis_data

CON_COEFFS_DEP = \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/con_coeffs.o : $(SRCDIR)/con_coeffs.f90 $(CON_COEFFS_DEP)
con_coeffs.o_name = con_coeffs

$(OBJDIR)/fatal.o : $(SRCDIR)/fatal.f90
fatal.o_name = fatal

G94_WAVEFUNCTION_DEP = \
$(OBJDIR)/paramfile.o
$(OBJDIR)/g94_wavefunction.o : $(SRCDIR)/g94_wavefunction.f90 $(G94_WAVEFUNCTION_DEP)
g94_wavefunction.o_name = g94_wavefunction

G_D_TYPE_DEP = \
$(OBJDIR)/paramfile.o
$(OBJDIR)/g_d_type.o : $(SRCDIR)/g_d_type.f90 $(G_D_TYPE_DEP)
g_d_type.o_name = g_d_type

G_S_TYPE_DEP = \
$(OBJDIR)/paramfile.o
$(OBJDIR)/g_s_type.o : $(SRCDIR)/g_s_type.f90 $(G_S_TYPE_DEP)
g_s_type.o_name = g_s_type

$(OBJDIR)/gaussiantoqmc.o : $(SRCDIR)/gaussiantoqmc.f90
gaussiantoqmc.o_name = gaussiantoqmc

GET_GAUSS_VERSION_DEP = \
$(OBJDIR)/awk_like.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/get_gauss_version.o : $(SRCDIR)/get_gauss_version.f90 $(GET_GAUSS_VERSION_DEP)
get_gauss_version.o_name = get_gauss_version

INTEG_PARAMS_DEP = \
$(OBJDIR)/paramfile.o
$(OBJDIR)/integ_params.o : $(SRCDIR)/integ_params.f90 $(INTEG_PARAMS_DEP)
integ_params.o_name = integ_params

$(OBJDIR)/max_coincidence.o : $(SRCDIR)/max_coincidence.f90
max_coincidence.o_name = max_coincidence

NORMALIZATION_CHECK_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o \
$(OBJDIR)/integ_params.o
$(OBJDIR)/normalization_check.o : $(SRCDIR)/normalization_check.f90 $(NORMALIZATION_CHECK_DEP)
normalization_check.o_name = normalization_check

NORMALIZE_CI_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/normalize_ci.o : $(SRCDIR)/normalize_ci.f90 $(NORMALIZE_CI_DEP)
normalize_ci.o_name = normalize_ci

$(OBJDIR)/numsrt_2way.o : $(SRCDIR)/numsrt_2way.f90
numsrt_2way.o_name = numsrt_2way

PACK_EVCOEFFS_DEP = \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/pack_evcoeffs.o : $(SRCDIR)/pack_evcoeffs.f90 $(PACK_EVCOEFFS_DEP)
pack_evcoeffs.o_name = pack_evcoeffs

$(OBJDIR)/paramfile.o : $(SRCDIR)/paramfile.f90
paramfile.o_name = paramfile

PSI_DEP = \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/psi.o : $(SRCDIR)/psi.f90 $(PSI_DEP)
psi.o_name = psi

QMC_WRITE_DEP = \
$(OBJDIR)/awk_like.o \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/qmc_write.o : $(SRCDIR)/qmc_write.f90 $(QMC_WRITE_DEP)
qmc_write.o_name = qmc_write

RE_SUM_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/re_sum.o : $(SRCDIR)/re_sum.f90 $(RE_SUM_DEP)
re_sum.o_name = re_sum

READ_G9XOUT_DEP = \
$(OBJDIR)/awk_like.o \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/read_G9xout.o : $(SRCDIR)/read_G9xout.f90 $(READ_G9XOUT_DEP)
read_G9xout.o_name = read_G9xout

READ_FCHK_DEP = \
$(OBJDIR)/awk_like.o \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/read_fchk.o : $(SRCDIR)/read_fchk.f90 $(READ_FCHK_DEP)
read_fchk.o_name = read_fchk

$(OBJDIR)/rejig.o : $(SRCDIR)/rejig.f90
rejig.o_name = rejig

RESUM_CAS_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o \
$(OBJDIR)/rejig.o
$(OBJDIR)/resum_cas.o : $(SRCDIR)/resum_cas.f90 $(RESUM_CAS_DEP)
resum_cas.o_name = resum_cas

SHELL_CENTRES_DEP = \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/shell_centres.o : $(SRCDIR)/shell_centres.f90 $(SHELL_CENTRES_DEP)
shell_centres.o_name = shell_centres

SUM_DEGEN_EXCITE_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/sum_degen_excite.o : $(SRCDIR)/sum_degen_excite.f90 $(SUM_DEGEN_EXCITE_DEP)
sum_degen_excite.o_name = sum_degen_excite

USER_CONTROL_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o
$(OBJDIR)/user_control.o : $(SRCDIR)/user_control.f90 $(USER_CONTROL_DEP)
user_control.o_name = user_control

WFN_CONSTRUCT_DEP = \
$(OBJDIR)/cis_data.o \
$(OBJDIR)/g94_wavefunction.o \
$(OBJDIR)/integ_params.o
$(OBJDIR)/wfn_construct.o : $(SRCDIR)/wfn_construct.f90 $(WFN_CONSTRUCT_DEP)
wfn_construct.o_name = wfn_construct

WFN_TEST_DEP = \
$(OBJDIR)/g94_wavefunction.o \
$(OBJDIR)/integ_params.o
$(OBJDIR)/wfn_test.o : $(SRCDIR)/wfn_test.f90 $(WFN_TEST_DEP)
wfn_test.o_name = wfn_test
# End module dependencies

# How to make the objects
$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info F90 utils $(UTILEXE)\
        $($(@F)_name).f90' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)\
        && $(F90_UTILS) -c $(FFLAGS_UTILS) $(SRCDIR)/$($(@F)_name).f90
# How to link the binary
$(BINDIR)/$(UTILEXE) : $(OBJ)
	@[ -d $(BINDIR) ] || mkdir -p $(BINDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info LDF90 utils $(UTILEXE)' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ;\
        $(LDF90_UTILS) $(LDFLAGS_UTILS) -o $(@) $(OBJ) $(LDLIBS_UTILS)

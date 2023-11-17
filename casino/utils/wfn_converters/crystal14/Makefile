#------------------------------------------------------------------#
# Makefile for crysgen14 programme.                                #
#                                                                  #
# MDT 9.2008                                                       #
#------------------------------------------------------------------#
# Basic stuff
UTILEXE = crysgen14
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
$(OBJDIR)/casino_interface.o \
$(OBJDIR)/condft.o \
$(OBJDIR)/convrs.o \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/crysgen14.o \
$(OBJDIR)/estros.o \
$(OBJDIR)/expo_module.o \
$(OBJDIR)/expt.o \
$(OBJDIR)/expu.o \
$(OBJDIR)/lmaxxx.o \
$(OBJDIR)/matout.o \
$(OBJDIR)/mcm.o \
$(OBJDIR)/memory_screen.o \
$(OBJDIR)/paral1_module.o \
$(OBJDIR)/rotation_matrices.o \
$(OBJDIR)/rotmatrix.o \
$(OBJDIR)/shell_info.o \
$(OBJDIR)/sphfac_module.o \
$(OBJDIR)/vrslat.o
# End object-file list

# Module dependencies
CASINO_INTERFACE_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/expo_module.o
$(OBJDIR)/casino_interface.o : $(SRCDIR)/casino_interface.f90 $(CASINO_INTERFACE_DEP)
casino_interface.o_name = casino_interface

CONDFT_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/lmaxxx.o \
$(OBJDIR)/shell_info.o \
$(OBJDIR)/sphfac_module.o
$(OBJDIR)/condft.o : $(SRCDIR)/condft.f90 $(CONDFT_DEP)
condft.o_name = condft

CONVRS_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/lmaxxx.o
$(OBJDIR)/convrs.o : $(SRCDIR)/convrs.f90 $(CONVRS_DEP)
convrs.o_name = convrs

CRYAPI_INP_C09_DEP = \
$(OBJDIR)/memory_screen.o
$(OBJDIR)/cryapi_inp_c09.o : $(SRCDIR)/cryapi_inp_c09.f90 $(CRYAPI_INP_C09_DEP)
cryapi_inp_c09.o_name = cryapi_inp_c09

CRYSGEN14_DEP = \
$(OBJDIR)/casino_interface.o \
$(OBJDIR)/cryapi_inp_c09.o
$(OBJDIR)/crysgen14.o : $(SRCDIR)/crysgen14.f90 $(CRYSGEN14_DEP)
crysgen14.o_name = crysgen14

ESTROS_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/expo_module.o \
$(OBJDIR)/memory_screen.o \
$(OBJDIR)/rotmatrix.o
$(OBJDIR)/estros.o : $(SRCDIR)/estros.f90 $(ESTROS_DEP)
estros.o_name = estros

EXPO_MODULE_DEP = \
$(OBJDIR)/cryapi_inp_c09.o
$(OBJDIR)/expo_module.o : $(SRCDIR)/expo_module.f90 $(EXPO_MODULE_DEP)
expo_module.o_name = expo_module

EXPT_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/expo_module.o
$(OBJDIR)/expt.o : $(SRCDIR)/expt.f90 $(EXPT_DEP)
expt.o_name = expt

EXPU_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/expo_module.o
$(OBJDIR)/expu.o : $(SRCDIR)/expu.f90 $(EXPU_DEP)
expu.o_name = expu

$(OBJDIR)/lmaxxx.o : $(SRCDIR)/lmaxxx.f90
lmaxxx.o_name = lmaxxx

MATOUT_DEP = \
$(OBJDIR)/cryapi_inp_c09.o
$(OBJDIR)/matout.o : $(SRCDIR)/matout.f90 $(MATOUT_DEP)
matout.o_name = matout

MCM_DEP = \
$(OBJDIR)/cryapi_inp_c09.o
$(OBJDIR)/mcm.o : $(SRCDIR)/mcm.f90 $(MCM_DEP)
mcm.o_name = mcm

$(OBJDIR)/memory_screen.o : $(SRCDIR)/memory_screen.f90
memory_screen.o_name = memory_screen

$(OBJDIR)/paral1_module.o : $(SRCDIR)/paral1_module.f90
paral1_module.o_name = paral1_module

ROTATION_MATRICES_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/lmaxxx.o \
$(OBJDIR)/paral1_module.o \
$(OBJDIR)/rotmatrix.o \
$(OBJDIR)/shell_info.o \
$(OBJDIR)/sphfac_module.o
$(OBJDIR)/rotation_matrices.o : $(SRCDIR)/rotation_matrices.f90 $(ROTATION_MATRICES_DEP)
rotation_matrices.o_name = rotation_matrices

ROTMATRIX_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/lmaxxx.o
$(OBJDIR)/rotmatrix.o : $(SRCDIR)/rotmatrix.f90 $(ROTMATRIX_DEP)
rotmatrix.o_name = rotmatrix

SHELL_INFO_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/lmaxxx.o
$(OBJDIR)/shell_info.o : $(SRCDIR)/shell_info.f90 $(SHELL_INFO_DEP)
shell_info.o_name = shell_info

SPHFAC_MODULE_DEP = \
$(OBJDIR)/cryapi_inp_c09.o \
$(OBJDIR)/lmaxxx.o
$(OBJDIR)/sphfac_module.o : $(SRCDIR)/sphfac_module.f90 $(SPHFAC_MODULE_DEP)
sphfac_module.o_name = sphfac_module

VRSLAT_DEP = \
$(OBJDIR)/cryapi_inp_c09.o
$(OBJDIR)/vrslat.o : $(SRCDIR)/vrslat.f90 $(VRSLAT_DEP)
vrslat.o_name = vrslat
# End module dependencies
# How to make the objects
$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info F90 utils $(UTILEXE)\
        $($(@F)_name).f90' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)\
        && $(F90_UTILS) -c $(FFLAGS_UTILS) $(SRCDIR)/$($(@F)_name).f90
$(BINDIR)/$(UTILEXE) : $(OBJ)
	@[ -d $(BINDIR) ] || mkdir -p $(BINDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info LDF90 utils $(UTILEXE)' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ;\
        $(LDF90_UTILS) $(LDFLAGS_UTILS) -o $(@) $(OBJ) $(LDLIBS_UTILS)

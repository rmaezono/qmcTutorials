#------------------------------------------------------------------#
# Makefile for plot_mpc program.                                   #
#                                                                  #
# MDT 5.2006                                                       #
#------------------------------------------------------------------#
# Basic stuff
UTILEXE = plot_mpc
SRCDIR  = $(PWD)
UTILDIR = $(SRCDIR)/..
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
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf.o \
$(OBJDIR)/esdf_key.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/plot_mpc.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
# End object-file list

# Module dependencies
$(OBJDIR)/dsp.o : $(SRCDIR)/dsp.f90
dsp.o_name = dsp

ESDF_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf_key.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/run_control.o
$(OBJDIR)/esdf.o : $(SRCDIR)/esdf.f90 $(ESDF_DEP)
esdf.o_name = esdf

$(OBJDIR)/esdf_key.o : $(SRCDIR)/esdf_key.f90
esdf_key.o_name = esdf_key

FORMAT_UTILS_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/format_utils.o : $(SRCDIR)/format_utils.f90 $(FORMAT_UTILS_DEP)
format_utils.o_name = format_utils

PLOT_MPC_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/run_control.o
$(OBJDIR)/plot_mpc.o : $(SRCDIR)/plot_mpc.f90 $(PLOT_MPC_DEP)
plot_mpc.o_name = plot_mpc

RUN_CONTROL_DEP = \
$(OBJDIR)/format_utils.o
$(OBJDIR)/run_control.o : $(SRCDIR)/run_control.f90 $(RUN_CONTROL_DEP)
run_control.o_name = run_control

$(OBJDIR)/store.o : $(SRCDIR)/store.f90
store.o_name = store
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

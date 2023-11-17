#----------------------------------------------------------------------#
# Makefile for LOUIS program                                           #
# ==========================                                           #
#                                                                      #
# Note this Makefile is written so that LOUIS can be compiled with     #
# features (opt/debug/prof/dev etc.) but this will only work if        #
# you specifically type 'make debug' etc. in this directory. Compiling #
# the utilities through the utils Makefile will always result in       #
# compilation of the opt version.                                      #
#                                                                      #
# MDT 12.2009                                                          #
#----------------------------------------------------------------------#

# Basic stuff
UTILEXE = louis
SRCDIR  = $(PWD)
UTILDIR = $(SRCDIR)/../..
TOOLDIR = $(UTILDIR)/../build-tools
F90WRAP = $(TOOLDIR)/f90wrap
INCBASE = $(UTILDIR)/../arch/data
SILENT  = 0
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

start: default
.PHONY: start

# Default feature and version.
VERSION    = opt
FEATURES   =
USE_OPENMP = no
USE_SHM    = no

# Feature and version selection targets.
.PHONY: NoFeatures NoFeatures/opt NoFeatures/dev NoFeatures/debug\
 NoFeatures/prof nofeatures nofeatures/opt nofeatures/dev nofeatures/debug\
 nofeatures/prof opt dev debug prof
NoFeatures NoFeatures/opt nofeatures nofeatures/opt opt:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=opt
NoFeatures/dev nofeatures/dev dev:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=dev
NoFeatures/debug nofeatures/debug debug:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=debug
NoFeatures/prof nofeatures/prof prof:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=prof

# System-independent variables
BINDIR       = $(UTILDIR)/../bin_qmc/utils/$(CASINO_ARCH)/$(FEATURES)/$(VERSION)
OBJDIR       = $(UTILDIR)/zlib/$(CASINO_ARCH)/$(FEATURES)/$(VERSION)/$(UTILEXE)/

MKOBJDIR            = [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
MKBINDIR            = [ -d $(BINDIR) ] || mkdir -p $(BINDIR)
ETIME_yes           = $(OBJDIR)/etime.o
ETIME_no            =
ETIMESRC_yes        = $(SRCDIR)/etime.c
ETIMESRC_no         =
MPI_VERSION         = 2
COMMS_TYPE_single   = serial
COMMS_TYPE_parallel = parallel_mpi$(MPI_VERSION)
COMMS_TYPE_cluster  = parallel_mpi$(MPI_VERSION)

# System-independent feature variables
WINWRAP_yes        = $(TOOLDIR)/winwrap
WINWRAP_no         =

# Defaults for system-dependent variables
TYPE                       = single
NATIVE_WINDOWS             = no
NEED_ETIME                 = no
MODNAME_BUG                = 0
LDF90                      = $(F90)
LDFLAGS_opt                = $(FFLAGS_opt)
LDFLAGS_all                = $(FFLAGS_all)
ENVIRONMENT_COMMAND        = :

# System dependent include files
include $(INCBASE)/$(CASINO_ARCH).arch

# System-dependent variables from data in include files
COMMS_FILE     = comms_$(COMMS_TYPE_$(TYPE)).o
ETIME          = $(ETIME_$(NEED_ETIME))
ETIMESRC       = $(ETIMESRC_$(NEED_ETIME))
WINWRAP        = $(WINWRAP_$(NATIVE_WINDOWS))
FFLAGS         = $(FFLAGS_$(VERSION)) $(FFLAGS_all) $(INCLUDE_DIR)
CFLAGS         = $(CFLAGS_$(VERSION)) $(CFLAGS_all) $(CFLAGS_F90_INTERFACE)\
 $(INCLUDE_DIR)
LDFLAGS        = $(LDFLAGS_$(VERSION)) $(LDFLAGS_all) $(LIB_PATH)

# Main target
.PHONY: default
default: $(BINDIR)/$(UTILEXE)

# 'Clean' target
.PHONY: clean
clean:
	@rm -rf $(OBJDIR) $(BINDIR)/$(UTILEXE)

# Object-file list
OBJ = \
$(OBJDIR)/bulirsch_stoer.o \
$(OBJDIR)/choose_wfn.o \
$(OBJDIR)/$(COMMS_FILE) \
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf.o \
$(OBJDIR)/esdf_key.o \
$(OBJDIR)/eval_density.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/input.o \
$(OBJDIR)/louis.o \
$(OBJDIR)/numerical.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/random_numbers.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/runge_kutta.o \
$(OBJDIR)/scaled_sin_wfn.o \
$(OBJDIR)/sin_wfn.o \
$(OBJDIR)/store.o
# End object-file list

# Module dependencies
BULIRSCH_STOER_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/bulirsch_stoer.o : $(SRCDIR)/bulirsch_stoer.f90 $(BULIRSCH_STOER_DEP)
bulirsch_stoer.o_name = bulirsch_stoer

CHOOSE_WFN_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/scaled_sin_wfn.o \
$(OBJDIR)/sin_wfn.o \
$(OBJDIR)/store.o
$(OBJDIR)/choose_wfn.o : $(SRCDIR)/choose_wfn.f90 $(CHOOSE_WFN_DEP)
choose_wfn.o_name = choose_wfn

$(OBJDIR)/comms_parallel_mpi1.o : $(SRCDIR)/comms_parallel_mpi1.f90
comms_parallel_mpi1.o_name = comms_parallel_mpi1

$(OBJDIR)/comms_parallel_mpi2.o : $(SRCDIR)/comms_parallel_mpi2.f90
comms_parallel_mpi2.o_name = comms_parallel_mpi2

COMMS_SERIAL_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/comms_serial.o : $(SRCDIR)/comms_serial.f90 $(COMMS_SERIAL_DEP)
comms_serial.o_name = comms_serial

$(OBJDIR)/dsp.o : $(SRCDIR)/dsp.f90
dsp.o_name = dsp

ESDF_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf_key.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/esdf.o : $(SRCDIR)/esdf.f90 $(ESDF_DEP)
esdf.o_name = esdf

$(OBJDIR)/esdf_key.o : $(SRCDIR)/esdf_key.f90
esdf_key.o_name = esdf_key

EVAL_DENSITY_DEP = \
$(OBJDIR)/choose_wfn.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/eval_density.o : $(SRCDIR)/eval_density.f90 $(EVAL_DENSITY_DEP)
eval_density.o_name = eval_density

FILE_UTILS_DEP = \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/file_utils.o : $(SRCDIR)/file_utils.f90 $(FILE_UTILS_DEP)
file_utils.o_name = file_utils

FORMAT_UTILS_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/format_utils.o : $(SRCDIR)/format_utils.f90 $(FORMAT_UTILS_DEP)
format_utils.o_name = format_utils

INPUT_DEP = \
$(OBJDIR)/$(COMMS_FILE) \
$(OBJDIR)/esdf.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/random_numbers.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/input.o : $(SRCDIR)/input.f90 $(INPUT_DEP)
input.o_name = input

LOUIS_DEP = \
$(OBJDIR)/bulirsch_stoer.o \
$(OBJDIR)/choose_wfn.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf.o \
$(OBJDIR)/eval_density.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/input.o \
$(OBJDIR)/numerical.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/random_numbers.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/runge_kutta.o \
$(OBJDIR)/scaled_sin_wfn.o \
$(OBJDIR)/sin_wfn.o \
$(OBJDIR)/store.o
$(OBJDIR)/louis.o : $(SRCDIR)/louis.f90 $(LOUIS_DEP)
louis.o_name = louis

NUMERICAL_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/numerical.o : $(SRCDIR)/numerical.f90 $(NUMERICAL_DEP)
numerical.o_name = numerical

PARALLEL_DEP = \
$(OBJDIR)/$(COMMS_FILE) \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/parallel.o : $(SRCDIR)/parallel.f90 $(PARALLEL_DEP)
parallel.o_name = parallel

RANDOM_NUMBERS_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/random_numbers.o : $(SRCDIR)/random_numbers.f90 $(RANDOM_NUMBERS_DEP)
random_numbers.o_name = random_numbers

RUN_CONTROL_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/store.o
$(OBJDIR)/run_control.o : $(SRCDIR)/run_control.f90 $(RUN_CONTROL_DEP)
run_control.o_name = run_control

RUNGE_KUTTA_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/runge_kutta.o : $(SRCDIR)/runge_kutta.f90 $(RUNGE_KUTTA_DEP)
runge_kutta.o_name = runge_kutta

SCALED_SIN_WFN_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/scaled_sin_wfn.o : $(SRCDIR)/scaled_sin_wfn.f90 $(SCALED_SIN_WFN_DEP)
scaled_sin_wfn.o_name = scaled_sin_wfn

SIN_WFN_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/sin_wfn.o : $(SRCDIR)/sin_wfn.f90 $(SIN_WFN_DEP)
sin_wfn.o_name = sin_wfn

STORE_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/store.o : $(SRCDIR)/store.f90 $(STORE_DEP)
store.o_name = store
# End module dependencies

# How to compile things
$(OBJ):
	@$(MKOBJDIR) ;\
        action=`bash -c '$(BASH_PRETTY) ; compile_info F90 utils louis\
         $($(@F)_name).f90'` ;\
        inaction=`bash -c '$(BASH_PRETTY) ; nocompile_info F90 utils louis\
         $($(@F)_name).f90'` ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(F90WRAP) $($(@F)_name) $(SRCDIR)/$($(@F)_name).f90 "$(WINWRAP)\
        $(F90) -c $(FFLAGS) $(SRCDIR)/$($(@F)_name).f90 -I$(SRCDIR)/.."\
        --mods $($($(@F)_name)_mods) --deps $($($(@F)_name)_depmod) --silent\
        $(SILENT) --action "$$action" --inaction "$$inaction"\
        --modname_bug $(MODNAME_BUG)

$(ETIME): $(ETIMESRC)
	@$(MKOBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info CC utils louis etime.c' ;\
        $(ENVIRONMENT_COMMAND) ;\
        cd $(OBJDIR) && $(WINWRAP) $(CC) $(CFLAGS) $(CFLAGS_ETIME)\
        -c $(ETIMESRC) $(INCLUDE_DIR)

$(BINDIR)/$(UTILEXE): $(OBJ) $(ETIME)
	@$(MKBINDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info LDF90 utils louis $(UTILEXE)' ;\
        $(ENVIRONMENT_COMMAND) ;\
        $(WINWRAP) $(LDF90) $(LDFLAGS) -o $(BINDIR)/$(UTILEXE) $(OBJ)\
        $(ETIME)


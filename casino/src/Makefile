#----------------------------------------------------------------------#
# CASINO sources Makefile                                              #
#                                                                      #
# All system dependent stuff is derived from the relevant include      #
# file in the CASINO/arch/data directory. All object files are         #
# put in an appropriate subdirectory of CASINO/src/zlib to reduce      #
# clutter and to allow different sets of object files to be maintained #
# for different machines and compilation levels.                       #
#                                                                      #
# The machine architecture should be defined before compilation by     #
# setting the CASINO_ARCH environment variable.                        #
#                                                                      #
# The code is compiled with full optimization by default; debugging,   #
# profiled and development versions can also be compiled (see below).  #
# The binaries and object files for different versions are kept        #
# separate. The different versions can be accessed using appropriate   #
# command-line flags for runqmc.                                       #
#                                                                      #
# IMPORTANT: the module dependencies are generated automatically with  #
# the 'update_makefile.sh' script. DO NOT MODIFY THE DEPENDENCIES BY   #
# HAND!!                                                               #
#                                                                      #
# Generic form of target:                                              #
#    <features>/<version>.<action>                                     #
# where <version> can be:                                              #
#  - opt  : use full optimization (default).                           #
#  - debug: use debugging compiler flags.                              #
#  - dev  : use full optimization, but keep object files and binary    #
#           separate from 'opt' version.                               #
#  - prof : use profiling compiler flags (if a profiler is available). #
# <features> can be:                                                   #
#  - NoFeatures / nofeatures : no special features (default).          #
#  - Openmp / openmp         : build OpenMP support.                   #
#  - Shm / shm               : build shared memory (SHM) support       #
#  - OpenmpShm / openmpshm   : build OpenMP and SHM support (broken?). #
# and <action> can be:                                                 #
#  - <empty>: perform compilation (default).                           #
#  - clean  : remove object files.                                     #
#  - vclean : remove object files and binaries.                        #
#                                                                      #
# E.g., one can type:                                                  #
#  make                                                                #
#  make opt                                                            #
#  make NoFeatures                                                     #
#  make nofeatures                                                     #
#  make NoFeatures/opt                                                 #
# to the same effect.                                                  #
#                                                                      #
# For the 'clean' and 'vclean' targets, 'all' can be used to refer     #
# to all features or versions, e.g.:                                   #
#  make all.clean          # same as 'make NoFeatures/all.vclean'      #
#  make all/dev.vclean                                                 #
#  make Openmp/all.vclean                                              #
#  make all/all.vclean                                                 #
#                                                                      #
# There are also advanced options that can be modified at compile time #
# should you need to. They are especially useful when trying to come   #
# up with an include file for your CASINO_ARCH, or to build a          #
# different binary when developing or debugging. The main options are  #
# (specify as an argument of 'make', e.g. 'make SILENT=1'):            #
#                                                                      #
#  VERSION     : set to 'opt', 'dev', 'debug' or 'prof' for different  #
#                compiler flags and target directories.                #
#  USE_OPENMP  : set to 'yes' for OpenMP build, else to 'no'.          #
#  USE_SHM     : set to 'yes' for SHM build, else to 'no'.             #
#  EXECUTABLE  : sets the name of the binary. Default: 'casino'. This  #
#                can be used in combination with runqmc's '--binary'   #
#                option to maintain a set of executables during        #
#                developemnt.                                          #
#  F90         : name of the Fortran compiler or MPI-Fortran wrapper.  #
#  LDF90       : linker (F90 by default).                              #
#  CFLAGS_x    : C compiler flags (x=opt|dev|debug|prof)               #
#  FFLAGS_x    : Fortran compiler flags (x=opt|dev|debug|prof)         #
#  LDFLAGS_x   : linker flags (FFLAGS_x by default; x=opt|dev|debug|   #
#                prof).                                                #
#  LDLIBS_x    : library flags to pass to the linker ('' by default;   #
#                x=opt|dev|debug|prof)                                 #
#  SILENT      : set to 1 to inhibit compiler verbosity. Default: 0.   #
#  TYPE        : set to 'single' to compile CASINO against the fake    #
#                MPI library comms_serial. Else, it is assumed that    #
#                the MPI library will be used instead (specified as    #
#                an include, or implicitly via an 'mpif90'-type        #
#                wrapper.                                              #
#  NEED_ETIME  : 'yes' if your Fortran compiler does not support the   #
#                ETIME extension, 'no' otherwise (default).            #
#  HAVE_BLAS   : if you set HAVE_BLAS to 'yes' and BLAS_yes to the     #
#  BLAS_yes      appropriate -L and/or -l flags (e.g., '-lmkl -lguide' #
#                for the Intel MKL libraries), CASINO will be compiled #
#                against the specified libraries. If HAVE_BLAS is      #
#                'no', BLAS_yes is ignored and the BLAS library under  #
#                CASINO/lib/BLAS is used instead.                      #
#  HAVE_LAPACK : if you set HAVE_LAPACK to 'yes' and LAPACK_yes to the #
#  LAPACK_yes    appropriate -L and/or -l flags (e.g., '-lmkl_lapack'  #
#                for the Intel MKL LAPACK libraries), CASINO will be   #
#                compiled against the specified libraries. If          #
#                HAVE_LAPACK is 'no', LAPACK_yes is ignored and the    #
#                LAPACK library under CASINO/lib/LAPACK is used        #
#                instead.                                              #
#  SUPPORT_OPENMP   : set to 'yes' if your compiler supports OpenMP.   #
#  CFLAGS_OPENMP_yes: set to the additional flags required for OpenMP  #
#                     compilation (e.g. '-fopenmp').                   #
#  FFLAGS_OPENMP_yes: set to the additional flags required for OpenMP  #
#                     compilation (e.g. '-fopenmp').                   #
#  SUPPORT_SHM      : set to 'yes' if your architecture supports SYSV  #
#                     IPC for shared-memory calculations.              #
#                                                                      #
# The CASINO Team 1999-2013                                            #
#----------------------------------------------------------------------#
# Basic stuff
SRCDIR         = $(PWD)
TOOLDIR        = $(SRCDIR)/../build-tools
BINBASE        = $(SRCDIR)/../bin_qmc
LIBBASE        = $(SRCDIR)/../lib/zlib
OBJBASE        = $(SRCDIR)/zlib
INCBASE        = $(SRCDIR)/../arch/data
EXECUTABLE     = casino
F90WRAP        = $(TOOLDIR)/f90wrap
VERSION_FILE   = $(SRCDIR)/../VERSION
VERSION_SCRIPT = $(TOOLDIR)/update-version
SILENT         = 0
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

.PHONY: start
start: default

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
.PHONY: Openmp Openmp/opt Openmp/dev Openmp/debug Openmp/prof\
 openmp openmp/opt openmp/dev openmp/debug openmp/prof
Openmp Openmp/opt openmp openmp/opt:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=opt
Openmp/dev openmp/dev:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=dev
Openmp/debug openmp/debug:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=debug
Openmp/prof openmp/prof:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=prof
.PHONY: Shm Shm/opt Shm/dev Shm/debug Shm/prof shm shm/opt shm/dev shm/debug \
 shm/prof
Shm Shm/opt shm shm/opt:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=opt
Shm/dev shm/dev:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=dev
Shm/debug shm/debug:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=debug
Shm/prof shm/prof:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=prof
.PHONY: OpenmpShm OpenmpShm/opt OpenmpShm/dev OpenmpShm/debug OpenmpShm/prof\
 openmpshm openmpshm/opt openmpshm/dev openmpshm/debug openmpshm/prof
OpenmpShm OpenmpShm/opt openmpshm openmpshm/opt:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=opt
OpenmpShm/dev openmpshm/dev:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=dev
OpenmpShm/debug openmpshm/debug:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=debug
OpenmpShm/prof openmpshm/prof:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=prof

# Extract FEATURES flag from USE_* variables.
# IMPORTANT: should new features be added, continue using wiki-like
# capitalization and concatenating the strings in alphabetical order to form
# the FEATURES variable.
FEATURE_OPENMP_yes = Openmp
FEATURE_OPENMP_no  =
FEATURE_SHM_yes    = Shm
FEATURE_SHM_no     =
FEATURES = $(FEATURE_OPENMP_$(USE_OPENMP))$(FEATURE_SHM_$(USE_SHM))

# System-independent variables
BINDIR              = $(BINBASE)/$(CASINO_ARCH)/$(FEATURES)/$(VERSION)
LIBDIR              = $(LIBBASE)/$(CASINO_ARCH)
OBJDIR              = $(OBJBASE)/$(CASINO_ARCH)/$(FEATURES)/$(VERSION)
LAPACK_no           = $(LIBDIR)/liblapack.a
LDLAPACK_no         = -L$(LIBDIR) -llapack
BLAS_no             = $(LIBDIR)/libblas.a
BLAS_yes            =
LDBLAS_no           = -L$(LIBDIR) -lblas
LDBLAS_yes          =
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
FFLAGS_OPENMP_no   =
CFLAGS_OPENMP_no   =
#CXXFLAGS_OPENMP_no =
LDFLAGS_OPENMP_no  =
ALLOC_yes          = $(OBJDIR)/alloc_shm.o
ALLOC_no           =
ALLOCSRC_yes       = $(SRCDIR)/alloc_shm.c
ALLOCSRC_no        =
SHM_TYPE_yes       = smp
SHM_TYPE_no        = nonsmp
WINWRAP_yes        = $(TOOLDIR)/winwrap
WINWRAP_no         =
ALLOC              = $(ALLOC_$(USE_SHM))
ALLOCSRC           = $(ALLOCSRC_$(USE_SHM))
SHM_TYPE           = $(SHM_TYPE_$(USE_SHM))
SHM_FILE           = shalloc_$(SHM_TYPE).o

# Defaults for system-dependent variables
TYPE                = single
NATIVE_WINDOWS      = no
NEED_ETIME          = no
HAVE_LAPACK         = no
LAPACK_yes          =
HAVE_BLAS           = no
LDLAPACK_yes        =
LDF90               = $(F90)
MODNAME_BUG         = 0
SUPPORT_OPENMP      = no
SUPPORT_SHM         = yes
FFLAGS_NOOPT        = -O0
CFLAGS_SHM          = -DSHM_SYSV
FFLAGS_dev          = $(FFLAGS_opt)
CFLAGS_dev          = $(CFLAGS_opt)
#CXXFLAGS_dev        = $(CXXFLAGS_opt)
LDLIBS_dev          = $(LDLIBS_opt)
LDFLAGS_opt         = $(FFLAGS_opt)
LDFLAGS_dev         = $(FFLAGS_dev)
LDFLAGS_debug       = $(FFLAGS_debug)
LDFLAGS_prof        = $(FFLAGS_prof)
LDFLAGS_OPENMP_yes  = $(FFLAGS_OPENMP_yes)
ENVIRONMENT_COMMAND = :

# System-dependent include files
include $(INCBASE)/$(CASINO_ARCH).arch

# System-dependent variables from data in include files
COMMS_FILE     = comms_$(COMMS_TYPE_$(TYPE)).o
ETIME          = $(ETIME_$(NEED_ETIME))
ETIMESRC       = $(ETIMESRC_$(NEED_ETIME))
LAPACK         = $(LAPACK_$(HAVE_LAPACK))
LDLAPACK       = $(LDLAPACK_$(HAVE_LAPACK))
BLAS           = $(BLAS_$(HAVE_BLAS))
LDBLAS         = $(LDBLAS_$(HAVE_BLAS))
WINWRAP        = $(WINWRAP_$(NATIVE_WINDOWS))
FFLAGS_OPENMP  = $(FFLAGS_OPENMP_$(USE_OPENMP))
CFLAGS_OPENMP  = $(CFLAGS_OPENMP_$(USE_OPENMP))
LDFLAGS_OPENMP = $(LDFLAGS_OPENMP_$(USE_OPENMP))
FFLAGS         = $(FFLAGS_$(VERSION)) $(FFLAGS_all) $(FFLAGS_OPENMP)\
 $(INCLUDE_DIR)
FFLAGS0        = $(FFLAGS_NOOPT) $(FFLAGS_all) $(FFLAGS_OPENMP) $(INCLUDE_DIR)
CFLAGS         = $(CFLAGS_$(VERSION)) $(CFLAGS_all) $(CFLAGS_F90_INTERFACE)\
 $(CFLAGS_OPENMP) $(INCLUDE_DIR)
#CXXFLAGS       = $(CXXFLAGS_$(VERSION)) $(CXXFLAGS_all) $(CXXFLAGS_OPENMP)\
# $(INCLUDE_DIR)
LDFLAGS        = $(LDFLAGS_$(VERSION)) $(LDFLAGS_all) $(LDFLAGS_OPENMP)\
 $(LIB_PATH)
LDLIBS         = $(LDLAPACK) $(LDBLAS) $(LDLIBS_$(VERSION)) $(LDLIBS_all)

# Main target
.PHONY: default
default: $(BLAS) $(LAPACK) info check_features $(BINDIR)/$(EXECUTABLE)

# Info.
.PHONY: info
info:
	@bash -c '$(BASH_PRETTY) ;\
        info_title "CASINO sources MAKE information" ;\
        info_colon "Compiling version" "$(VERSION)" ;\
        [ -z "$(FEATURES)" ] \
         && info_colon "Compiling features" "none" \
         || info_colon "Compiling features" "$(FEATURES)" ;\
        [ "$(NEED_ETIME)" = yes ] \
         && info_line "ETIME extension to be compiled separately." ;\
        [ "$(USE_OPENMP)" = yes ] && info_line "OpenMP enabled." ;\
        [ "$(USE_SHM)" = yes ] \
         && info_line "Shared memory allocation (SHM) enabled." ;\
        [ "$(NEED_ETIME)$(USE_SHM)" != nono ] \
         && info_equal "CC" "$(WINWRAP) $(CC) $(CFLAGS) $(INCLUDE_DIR)";\
        [ no != no ] \
         && info_equal "CXX" "$(WINWRAP) $(CXX) $(CXXFLAGS) $(INCLUDE_DIR)" ;\
        [ "$(MODNAME_BUG)" = 1 ] \
         && info_line "Compiler suffers from MOD name bug." ;\
        [ "$(ENVIRONMENT_COMMAND)" != : ]\
         && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
        info_equal "F90" "$(WINWRAP) $(F90) $(FFLAGS)" ;\
        info_equal "F90-O0" "$(WINWRAP) $(F90) $(FFLAGS0)" ;\
        info_equal "LDF90"\
         "$(WINWRAP) $(LDF90) $(LDFLAGS) $(ETIME) $(LDLIBS)" ;\
        info_end'

# Check features are available.
.PHONY: check_features
check_features:
	@die=0 ;\
        [ "$(USE_OPENMP)$(SUPPORT_OPENMP)" = yesno ] \
         && echo "OpenMP build not supported for this architecture."\
         && die=1 ; \
        [ "$(USE_SHM)$(SUPPORT_SHM)" = yesno ] \
         && echo "SHM build not supported for this architecture."\
         && die=1 ; \
        [ $$die = 1 ] && exit 1 || :

# Targets for deleting object files for one version and one feature.
.PHONY: clean
clean:
	@[ -z "$(FEATURES)" ]\
         && echo "Cleaning object files [$(CASINO_ARCH), $(VERSION)]."\
         || echo "Cleaning object files [$(CASINO_ARCH), $(VERSION),\
 $(FEATURES)]." ;\
        rm -rf $(OBJDIR)
.PHONY: NoFeatures.clean NoFeatures/opt.clean NoFeatures/dev.clean\
 NoFeatures/debug.clean NoFeatures/prof.clean nofeatures.clean\
 nofeatures/opt.clean nofeatures/dev.clean nofeatures/debug.clean\
 nofeatures/prof.clean opt.clean dev.clean debug.clean prof.clean 
NoFeatures.clean NoFeatures/opt.clean nofeatures.clean nofeatures/opt.clean\
 opt.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=opt clean
NoFeatures/dev.clean nofeatures/dev.clean dev.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=dev clean
NoFeatures/debug.clean nofeatures/debug.clean debug.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=debug clean
NoFeatures/prof.clean nofeatures/prof.clean prof.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=prof clean
.PHONY: Openmp.clean Openmp/opt.clean Openmp/dev.clean Openmp/debug.clean\
 Openmp/prof.clean openmp.clean openmp/opt.clean openmp/dev.clean\
 openmp/debug.clean openmp/prof.clean
Openmp.clean Openmp/opt.clean openmp.clean openmp/opt.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=opt clean
Openmp/dev.clean openmp/dev.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=dev clean
Openmp/debug.clean openmp/debug.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=debug clean
Openmp/prof.clean openmp/prof.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=prof clean
.PHONY: Shm.clean Shm/opt.clean Shm/dev.clean Shm/debug prof.clean\
 Shm/prof.clean shm.clean shm/opt.clean shm/dev.clean shm/debug prof.clean\
 shm/prof.clean
Shm.clean Shm/opt.clean shm.clean shm/opt.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=opt clean
Shm/dev.clean shm/dev.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=dev clean
Shm/debug.clean shm/debug.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=debug clean
Shm/prof.clean shm/prof.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=prof clean
.PHONY: OpenmpShm.clean OpenmpShm/opt.clean OpenmpShm/dev.clean\
 OpenmpShm/debug.clean OpenmpShm/prof.clean openmpshm.clean\
 openmpshm/opt.clean openmpshm/dev.clean openmpshm/debug.clean\
 openmpshm/prof.clean
OpenmpShm.clean OpenmpShm/opt.clean openmpshm.clean openmpshm/opt.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=opt clean
OpenmpShm/dev.clean openmpshm/dev.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=dev clean
OpenmpShm/debug.clean openmpshm/debug.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=debug clean
OpenmpShm/prof.clean openmpshm/prof.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=prof clean

# Targets for deleting object files for all versions and one feature.
.PHONY: all.clean NoFeatures/all.clean Openmp/all.clean Shm/all.clean\
 OpenmpShm/all.clean nofeature/all.clean openmp/all.cleanshm/all.clean\
 openmpshm/all.clean
all.clean NoFeatures/all.clean nofeatures/all.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=opt clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=no VERSION=dev clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=no VERSION=debug clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=no VERSION=prof clean
Openmp/all.clean openmp/all.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=opt clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=dev clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=debug clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=prof clean
Shm/all.clean shm/all.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=opt clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=dev clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=debug clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=prof clean
OpenmpShm/all.clean openmpshm/all.clean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=opt clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=dev clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=debug clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=prof clean

# Targets for deleting object files for one version and all features.
.PHONY: all/opt.clean all/debug.clean all/dev.clean all/prof.clean
all/opt.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=opt clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=opt clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=opt clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=opt clean
all/debug.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=debug clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=debug clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=debug clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=debug clean
all/dev.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=dev clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=dev clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=dev clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=dev clean
all/prof.clean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=prof clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=prof clean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=prof clean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=prof clean

# Targets for deleting object files for all versions and all features (along 
# with their base directory).
.PHONY: all/all.clean
all/all.clean:
	@echo "Cleaning object files [$(CASINO_ARCH), all versions,\
 all features]" ;\
        rm -rf $(OBJBASE)/$(CASINO_ARCH)

# Targets for deleting object files and binary for one version and one feature.
.PHONY: vclean
vclean:
	@[ -z "$(FEATURES)" ]\
         && echo "Cleaning object files and binary [$(CASINO_ARCH),\
 $(VERSION)]."\
         || echo "Cleaning object files and binary [$(CASINO_ARCH),\
 $(VERSION), $(FEATURES)]." ;\
        rm -rf $(OBJDIR) $(BINDIR)/$(EXECUTABLE)
.PHONY: NoFeatures.vclean NoFeatures/opt.vclean NoFeatures/dev.vclean\
 NoFeatures/debug.vclean NoFeatures/prof.vclean nofeatures.vclean\
 nofeatures/opt.vclean nofeatures/dev.vclean nofeatures/debug.vclean\
 nofeatures/prof.vclean opt.vclean dev.vclean debug.vclean prof.vclean 
NoFeatures.vclean NoFeatures/opt.vclean nofeatures.vclean nofeatures/opt.vclean\ opt.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=opt vclean
NoFeatures/dev.vclean nofeatures/dev.vclean dev.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=dev vclean
NoFeatures/debug.vclean nofeatures/debug.vclean debug.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=debug vclean
NoFeatures/prof.vclean nofeatures/prof.vclean prof.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=prof vclean
.PHONY: Openmp.vclean Openmp/opt.vclean Openmp/dev.vclean Openmp/debug.vclean\
 Openmp/prof.vclean openmp.vclean openmp/opt.vclean openmp/dev.vclean\
 openmp/debug.vclean openmp/prof.vclean
Openmp.vclean Openmp/opt.vclean openmp.vclean openmp/opt.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=opt vclean
Openmp/dev.vclean openmp/dev.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=dev vclean
Openmp/debug.vclean openmp/debug.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=debug vclean
Openmp/prof.vclean openmp/prof.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=prof vclean
.PHONY: Shm.vclean Shm/opt.vclean Shm/dev.vclean Shm/debug prof.vclean\
 Shm/prof.vclean shm.vclean shm/opt.vclean shm/dev.vclean shm/debug prof.vclean\
 shm/prof.vclean
Shm.vclean Shm/opt.vclean shm.vclean shm/opt.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=opt vclean
Shm/dev.vclean shm/dev.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=dev vclean
Shm/debug.vclean shm/debug.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=debug vclean
Shm/prof.vclean shm/prof.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=prof vclean
.PHONY: OpenmpShm.vclean OpenmpShm/opt.vclean OpenmpShm/dev.vclean\
 OpenmpShm/debug.vclean OpenmpShm/prof.vclean
OpenmpShm.vclean OpenmpShm/opt.vclean openmpshm.vclean openmpshm/opt.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=opt vclean
OpenmpShm/dev.vclean openmpshm/dev.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=dev vclean
OpenmpShm/debug.vclean openmpshm/debug.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=debug vclean
OpenmpShm/prof.vclean openmpshm/prof.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=prof vclean

# Targets for deleting object files and binaries for all versions and one
# feature.
.PHONY: all.vclean NoFeatures/all.vclean Openmp/all.vclean Shm/all.vclean\
 OpenmpShm/all.vclean noFeatures/all.vclean openmp/all.vclean shm/all.vclean\
 openmpshm/all.vclean
all.vclean NoFeatures/all.vclean nofeatures/all.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=opt vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=no VERSION=dev vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=no VERSION=debug vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=no VERSION=prof vclean
Openmp/all.vclean openmp/all.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=opt vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=dev vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=debug vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=prof vclean
Shm/all.vclean shm/all.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=opt vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=dev vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=debug vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=prof vclean
OpenmpShm/all.vclean openmpshm/all.vclean:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=opt vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=dev vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=debug vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=prof vclean

# Targets for deleting object files and binaries for one version and all
# features.
.PHONY: all/opt.vclean all/debug.vclean all/dev.vclean all/prof.vclean
all/opt.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=opt vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=opt vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=opt vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=opt vclean
all/debug.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=debug vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=debug vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=debug vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=debug vclean
all/dev.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=dev vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=dev vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=dev vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=dev vclean
all/prof.vclean:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=prof vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=prof vclean &&\
        $(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=prof vclean &&\
        $(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=prof vclean

# Targets for deleting object files and binaries for all versions and all
# features (along with their base directories).
.PHONY: all/all.vclean
all/all.vclean:
	@echo "Cleaning object files and binary [$(CASINO_ARCH), all versions,\
 all features]" ;\
        rm -rf $(OBJBASE)/$(CASINO_ARCH) $(BINBASE)/$(CASINO_ARCH)

# Object-file list
OBJ = \
$(OBJDIR)/AAA.o \
$(OBJDIR)/AAB.o \
$(OBJDIR)/AAC.o \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAE.o \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAH.o \
$(OBJDIR)/AAI.o \
$(OBJDIR)/AAJ.o \
$(OBJDIR)/AAK.o \
$(OBJDIR)/AAL.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAO.o \
$(OBJDIR)/AAP.o \
$(OBJDIR)/AAQ.o \
$(OBJDIR)/AAR.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/AAU.o \
$(OBJDIR)/AAV.o \
$(OBJDIR)/AAW.o \
$(OBJDIR)/AAX.o \
$(OBJDIR)/AAY.o \
$(OBJDIR)/AAZ.o \
$(OBJDIR)/ABA.o \
$(OBJDIR)/ABB.o \
$(OBJDIR)/ABC.o \
$(OBJDIR)/ABD.o \
$(OBJDIR)/ABE.o \
$(OBJDIR)/ABF.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABH.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABJ.o \
$(OBJDIR)/ABK.o \
$(OBJDIR)/ABL.o \
$(OBJDIR)/ABM.o \
$(OBJDIR)/ABN.o \
$(OBJDIR)/ABO.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABR.o \
$(OBJDIR)/ABS.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ABU.o \
$(OBJDIR)/ABV.o \
$(OBJDIR)/ABW.o \
$(OBJDIR)/ABX.o \
$(OBJDIR)/ABY.o \
$(OBJDIR)/ABZ.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACB.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACD.o \
$(OBJDIR)/ACE.o \
$(OBJDIR)/ACF.o \
$(OBJDIR)/ACG.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACI.o \
$(OBJDIR)/ACJ.o \
$(OBJDIR)/ACK.o \
$(OBJDIR)/ACL.o \
$(OBJDIR)/ACM.o \
$(OBJDIR)/ACN.o \
$(OBJDIR)/ACO.o \
$(OBJDIR)/ACP.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/ACR.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/$(COMMS_FILE) \
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf.o \
$(OBJDIR)/esdf_key.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/nl2sol.o \
$(OBJDIR)/openmp_base.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/$(SHM_FILE) \
$(OBJDIR)/singleton.o \
$(OBJDIR)/store.o
# End object-file list

# Module structure
AAA_mods = slaarnAAA
AAB_mods = slaarnAAB
AAC_mods = slaarnAAC
AAD_mods = slaarnAAD
AAE_mods = slaarnAAE
AAF_mods = slaarnAAF
AAG_mods = slaarnAAG
AAH_mods = slaarnAAH
AAI_mods = slaarnAAI
AAJ_mods = slaarnAAJ
AAK_mods = slaarnAAK
AAL_mods = slaarnAAL
AAM_mods = slaarnAAM
AAN_mods = slaarnAAN
AAO_mods = slaarnAAO
AAP_mods = slaarnAAP
AAQ_mods = slaarnAAQ calcr_xc_corr calcj_xc_corr calcf_xc_corr \
calcf_xc_corr2
AAR_mods = slaarnAAR
AAS_mods = slaarnAAS
AAT_mods = slaarnAAT
AAU_mods = slaarnAAU
AAV_mods = slaarnAAV
AAW_mods = slaarnAAW
AAX_mods = slaarnAAX
AAY_mods = slaarnAAY
AAZ_mods = slaarnAAZ
ABA_mods = slaarnABA
ABB_mods = slaarnABB
ABC_mods = slaarnABC
ABD_mods = slaarnABC_noopt
ABE_mods = slaarnABE
ABF_mods = slaarnABF
ABG_mods = slaarnABG
ABH_mods = slaarnABH
ABI_mods = slaarnABI objfn_inner objfn_outer
ABJ_mods = slaarnABJ
ABK_mods = slaarnABK
ABL_mods = slaarnABL
ABM_mods = slaarnABM
ABN_mods = slaarnABN
ABO_mods = main
ABP_mods = slaarnABP
ABQ_mods = slaarnABQ
ABR_mods = monte_carlo
ABS_mods = slaarnABS
ABT_mods = slaarnABT
ABU_mods = slaarnABU
ABV_mods = slaarnABV
ABW_mods = slaarnABW
ABX_mods = slaarnABX
ABY_mods = slaarnABY
ABZ_mods = slaarnABYter
ACA_mods = slaarnACA
ACB_mods = slaarnACB
ACC_mods = slaarnACC
ACD_mods = slaarnACD
ACE_mods = slaarnACE
ACF_mods = slaarnACF
ACG_mods = slaarnACG
ACH_mods = slaarnACH
ACI_mods = slaarnACI
ACJ_mods = slaarnACJ
ACK_mods = slaarnACK
ACL_mods = slaarnACL
ACM_mods = slaarnACM madr
ACN_mods = slaarnACM_linjas calcr_varmin_linjas calcj_varmin_linjas
ACO_mods = slaarnACO
ACP_mods = slaarnACP
ACQ_mods = slaarnACQ
ACR_mods = slaarnACR
ACS_mods = slaarnACS
casl_mods = casl
comms_parallel_mpi1_mods = comms
comms_parallel_mpi2_mods = comms
comms_serial_mods = comms
dsp_mods = dsp
esdf_mods = esdf
esdf_key_mods = esdf_key
file_utils_mods = file_utils
format_utils_mods = format_utils
nl2sol_mods = machine_constants toms573
openmp_base_mods = openmp_base
parallel_mods = parallel
run_control_mods = run_control
shalloc_nonsmp_mods = shalloc
shalloc_smp_mods = shalloc
singleton_mods = singleton
store_mods = store
# End module structure

# Module dependencies
AAA_depmod = dsp format_utils parallel run_control slaarnAAG slaarnAAM \
slaarnACG slaarnACH slaarnACQ slaarnACS store
AAA_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/ACG.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAA.o : $(SRCDIR)/AAA.f90 $(AAA_DEP)
AAA.o_name = AAA

AAB_depmod = dsp format_utils parallel run_control slaarnAAM slaarnAAT \
slaarnABG slaarnABS slaarnABT slaarnACA slaarnACH slaarnACS store
AAB_DEP = \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABS.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAB.o : $(SRCDIR)/AAB.f90 $(AAB_DEP)
AAB.o_name = AAB

AAC_depmod = dsp file_utils format_utils parallel run_control slaarnABG \
slaarnABT slaarnABU slaarnACA store
AAC_DEP = \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ABU.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAC.o : $(SRCDIR)/AAC.f90 $(AAC_DEP)
AAC.o_name = AAC

AAD_depmod = dsp slaarnABB slaarnABV
AAD_DEP = \
$(OBJDIR)/ABB.o \
$(OBJDIR)/ABV.o \
$(OBJDIR)/dsp.o
$(OBJDIR)/AAD.o : $(SRCDIR)/AAD.f90 $(AAD_DEP)
AAD.o_name = AAD

AAE_depmod = dsp file_utils format_utils openmp_base parallel run_control \
shalloc slaarnAAG slaarnAAN slaarnABG slaarnABI slaarnABP slaarnABQ \
slaarnABT slaarnACA store
AAE_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/openmp_base.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/$(SHM_FILE) \
$(OBJDIR)/store.o
$(OBJDIR)/AAE.o : $(SRCDIR)/AAE.f90 $(AAE_DEP)
AAE.o_name = AAE

AAF_depmod = dsp file_utils format_utils parallel run_control shalloc \
slaarnABG slaarnABT slaarnACC store
AAF_DEP = \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/$(SHM_FILE) \
$(OBJDIR)/store.o
$(OBJDIR)/AAF.o : $(SRCDIR)/AAF.f90 $(AAF_DEP)
AAF.o_name = AAF

AAG_depmod = dsp
AAG_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/AAG.o : $(SRCDIR)/AAG.f90 $(AAG_DEP)
AAG.o_name = AAG

AAH_depmod = file_utils format_utils parallel run_control slaarnAAD \
slaarnAAO slaarnAAS slaarnAAU slaarnABG slaarnABL slaarnABN slaarnABP \
slaarnABU slaarnACL slaarnACS store
AAH_DEP = \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAO.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/AAU.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABL.o \
$(OBJDIR)/ABN.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABU.o \
$(OBJDIR)/ACL.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAH.o : $(SRCDIR)/AAH.f90 $(AAH_DEP)
AAH.o_name = AAH

AAI_depmod = dsp run_control slaarnAAG slaarnAAT slaarnABG slaarnABP \
slaarnABT store
AAI_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAI.o : $(SRCDIR)/AAI.f90 $(AAI_DEP)
AAI.o_name = AAI

AAJ_depmod = dsp file_utils format_utils parallel run_control slaarnAAB \
slaarnAAD slaarnAAF slaarnAAG slaarnAAM slaarnAAQ slaarnABG slaarnABJ \
slaarnABK slaarnABS slaarnABT slaarnACA slaarnACC slaarnACD slaarnACE \
slaarnACH slaarnACS store
AAJ_DEP = \
$(OBJDIR)/AAB.o \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAQ.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABJ.o \
$(OBJDIR)/ABK.o \
$(OBJDIR)/ABS.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACD.o \
$(OBJDIR)/ACE.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAJ.o : $(SRCDIR)/AAJ.f90 $(AAJ_DEP)
AAJ.o_name = AAJ

AAK_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnABG slaarnABT store
AAK_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAK.o : $(SRCDIR)/AAK.f90 $(AAK_DEP)
AAK.o_name = AAK

AAL_depmod = casl dsp format_utils parallel run_control slaarnAAA slaarnAAF \
slaarnAAG slaarnAAH slaarnAAM slaarnAAS slaarnABG slaarnABT slaarnACA \
slaarnACC slaarnACH slaarnACS store
AAL_DEP = \
$(OBJDIR)/AAA.o \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAH.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAL.o : $(SRCDIR)/AAL.f90 $(AAL_DEP)
AAL.o_name = AAL

AAM_depmod = dsp run_control slaarnAAF slaarnAAG slaarnAAP slaarnAAQ \
slaarnABG slaarnABK slaarnABM slaarnABS slaarnABT slaarnACA slaarnACE \
slaarnACH slaarnACO slaarnACS store
AAM_DEP = \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAP.o \
$(OBJDIR)/AAQ.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABK.o \
$(OBJDIR)/ABM.o \
$(OBJDIR)/ABS.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACE.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACO.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAM.o : $(SRCDIR)/AAM.f90 $(AAM_DEP)
AAM.o_name = AAM

AAN_depmod = dsp format_utils openmp_base parallel run_control slaarnAAG \
slaarnAAP slaarnAAT slaarnABG slaarnABQ slaarnABT slaarnACA slaarnACC \
slaarnACI store
AAN_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAP.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACI.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/openmp_base.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAN.o : $(SRCDIR)/AAN.f90 $(AAN_DEP)
AAN.o_name = AAN

AAO_depmod = dsp file_utils format_utils parallel run_control slaarnAAD \
slaarnAAF slaarnAAG slaarnAAN slaarnABG slaarnABT slaarnACH store
AAO_DEP = \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAO.o : $(SRCDIR)/AAO.f90 $(AAO_DEP)
AAO.o_name = AAO

AAP_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnABG slaarnABT store
AAP_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAP.o : $(SRCDIR)/AAP.f90 $(AAP_DEP)
AAP.o_name = AAP

AAQ_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnAAN slaarnABG slaarnABT slaarnACA slaarnACC slaarnACH slaarnACQ \
slaarnACS store toms573
AAQ_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/nl2sol.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAQ.o : $(SRCDIR)/AAQ.f90 $(AAQ_DEP)
AAQ.o_name = AAQ

AAR_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnABG store
AAR_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAR.o : $(SRCDIR)/AAR.f90 $(AAR_DEP)
AAR.o_name = AAR

AAS_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnAAN slaarnAAT slaarnABG slaarnABP slaarnABQ slaarnABT slaarnACA \
slaarnACC slaarnACI store
AAS_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACI.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAS.o : $(SRCDIR)/AAS.f90 $(AAS_DEP)
AAS.o_name = AAS

AAT_depmod = dsp
AAT_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/AAT.o : $(SRCDIR)/AAT.f90 $(AAT_DEP)
AAT.o_name = AAT

AAU_depmod = dsp file_utils format_utils parallel run_control shalloc \
slaarnAAG slaarnAAI slaarnAAN slaarnAAT slaarnAAV slaarnAAW slaarnAAX \
slaarnAAY slaarnAAZ slaarnABA slaarnABG slaarnABI slaarnABP slaarnABT \
slaarnACA store
AAU_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAI.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/AAV.o \
$(OBJDIR)/AAW.o \
$(OBJDIR)/AAX.o \
$(OBJDIR)/AAY.o \
$(OBJDIR)/AAZ.o \
$(OBJDIR)/ABA.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/$(SHM_FILE) \
$(OBJDIR)/store.o
$(OBJDIR)/AAU.o : $(SRCDIR)/AAU.f90 $(AAU_DEP)
AAU.o_name = AAU

AAV_depmod = dsp slaarnAAI slaarnAAT slaarnABG slaarnABI slaarnABP \
slaarnABT slaarnACA store
AAV_DEP = \
$(OBJDIR)/AAI.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAV.o : $(SRCDIR)/AAV.f90 $(AAV_DEP)
AAV.o_name = AAV

AAW_depmod = dsp slaarnAAI slaarnAAT slaarnABG slaarnABI slaarnABP \
slaarnABT slaarnACA store
AAW_DEP = \
$(OBJDIR)/AAI.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAW.o : $(SRCDIR)/AAW.f90 $(AAW_DEP)
AAW.o_name = AAW

AAX_depmod = dsp slaarnAAG slaarnAAT slaarnABG slaarnABI slaarnABP \
slaarnABT store
AAX_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAX.o : $(SRCDIR)/AAX.f90 $(AAX_DEP)
AAX.o_name = AAX

AAY_depmod = dsp slaarnAAG slaarnAAT slaarnABG slaarnABP slaarnABT \
slaarnACA store
AAY_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAY.o : $(SRCDIR)/AAY.f90 $(AAY_DEP)
AAY.o_name = AAY

AAZ_depmod = dsp slaarnAAG slaarnAAT slaarnABG slaarnABP slaarnABT \
slaarnACA store
AAZ_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/AAZ.o : $(SRCDIR)/AAZ.f90 $(AAZ_DEP)
AAZ.o_name = AAZ

ABA_depmod = dsp slaarnAAG slaarnAAT slaarnABG slaarnABI slaarnABP \
slaarnABT store
ABA_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABA.o : $(SRCDIR)/ABA.f90 $(ABA_DEP)
ABA.o_name = ABA

ABB_depmod = casl dsp slaarnAAG slaarnABG store
ABB_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABB.o : $(SRCDIR)/ABB.f90 $(ABB_DEP)
ABB.o_name = ABB

ABC_depmod = casl dsp format_utils parallel run_control slaarnAAG slaarnABG \
slaarnABT slaarnACH store
ABC_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABC.o : $(SRCDIR)/ABC.f90 $(ABC_DEP)
ABC.o_name = ABC

ABD_depmod = dsp format_utils parallel run_control slaarnABT
ABD_DEP = \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o
$(OBJDIR)/ABD.o : $(SRCDIR)/ABD.f90 $(ABD_DEP)
ABD.o_name = ABD

ABE_depmod = casl dsp format_utils parallel run_control slaarnAAD slaarnAAG \
slaarnABG slaarnABT slaarnACH slaarnACQ store
ABE_DEP = \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABE.o : $(SRCDIR)/ABE.f90 $(ABE_DEP)
ABE.o_name = ABE

ABF_depmod = dsp file_utils format_utils parallel run_control singleton \
slaarnAAG slaarnABG slaarnABK slaarnABQ slaarnABT slaarnACA slaarnACQ store
ABF_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABK.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/singleton.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABF.o : $(SRCDIR)/ABF.f90 $(ABF_DEP)
ABF.o_name = ABF

ABG_depmod = dsp
ABG_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/ABG.o : $(SRCDIR)/ABG.f90 $(ABG_DEP)
ABG.o_name = ABG

ABH_depmod = casl dsp format_utils parallel run_control slaarnAAF slaarnAAN \
slaarnABC slaarnABC_noopt slaarnABG slaarnABT slaarnACA slaarnACH store
ABH_DEP = \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABC.o \
$(OBJDIR)/ABD.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABH.o : $(SRCDIR)/ABH.f90 $(ABH_DEP)
ABH.o_name = ABH

ABI_depmod = dsp format_utils parallel run_control slaarnAAG slaarnABG \
slaarnABQ slaarnABT slaarnACA store
ABI_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABI.o : $(SRCDIR)/ABI.f90 $(ABI_DEP)
ABI.o_name = ABI

ABJ_depmod = dsp file_utils format_utils parallel run_control slaarnAAB \
slaarnAAQ slaarnABG slaarnACE store
ABJ_DEP = \
$(OBJDIR)/AAB.o \
$(OBJDIR)/AAQ.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ACE.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABJ.o : $(SRCDIR)/ABJ.f90 $(ABJ_DEP)
ABJ.o_name = ABJ

ABK_depmod = dsp file_utils format_utils openmp_base parallel run_control \
singleton slaarnAAG slaarnAAQ slaarnAAS slaarnABG slaarnABI slaarnABT store
ABK_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAQ.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/openmp_base.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/singleton.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABK.o : $(SRCDIR)/ABK.f90 $(ABK_DEP)
ABK.o_name = ABK

ABL_depmod = dsp run_control slaarnABH slaarnABX
ABL_DEP = \
$(OBJDIR)/ABH.o \
$(OBJDIR)/ABX.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o
$(OBJDIR)/ABL.o : $(SRCDIR)/ABL.f90 $(ABL_DEP)
ABL.o_name = ABL

ABM_depmod = dsp file_utils format_utils parallel run_control store
ABM_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABM.o : $(SRCDIR)/ABM.f90 $(ABM_DEP)
ABM.o_name = ABM

ABN_depmod = dsp file_utils format_utils parallel run_control slaarnAAD \
slaarnAAG slaarnAAN slaarnAAS slaarnABG slaarnABQ slaarnABT slaarnACH store
ABN_DEP = \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABN.o : $(SRCDIR)/ABN.f90 $(ABN_DEP)
ABN.o_name = ABN

ABO_depmod = format_utils openmp_base parallel run_control shalloc store
ABO_DEP = \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/openmp_base.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/$(SHM_FILE) \
$(OBJDIR)/store.o
$(OBJDIR)/ABO.o : $(SRCDIR)/ABO.f90 $(ABO_DEP)
ABO.o_name = ABO

ABP_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnABG slaarnABT store
ABP_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABP.o : $(SRCDIR)/ABP.f90 $(ABP_DEP)
ABP.o_name = ABP

ABQ_depmod = dsp slaarnABG
ABQ_DEP = \
$(OBJDIR)/ABG.o \
$(OBJDIR)/dsp.o
$(OBJDIR)/ABQ.o : $(SRCDIR)/ABQ.f90 $(ABQ_DEP)
ABQ.o_name = ABQ

ABR_depmod = casl dsp esdf file_utils format_utils openmp_base parallel \
run_control shalloc slaarnAAB slaarnAAC slaarnAAD slaarnAAE slaarnAAF \
slaarnAAG slaarnAAH slaarnAAJ slaarnAAK slaarnAAL slaarnAAN slaarnAAO \
slaarnAAP slaarnAAQ slaarnAAR slaarnAAS slaarnAAT slaarnAAU slaarnABE \
slaarnABF slaarnABG slaarnABI slaarnABJ slaarnABK slaarnABL slaarnABM \
slaarnABN slaarnABP slaarnABS slaarnABT slaarnABY slaarnABYter slaarnACA \
slaarnACB slaarnACC slaarnACE slaarnACF slaarnACI slaarnACJ slaarnACK \
slaarnACL slaarnACM slaarnACM_linjas slaarnACO slaarnACP slaarnACQ \
slaarnACR slaarnACS store
ABR_DEP = \
$(OBJDIR)/AAB.o \
$(OBJDIR)/AAC.o \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAE.o \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAH.o \
$(OBJDIR)/AAJ.o \
$(OBJDIR)/AAK.o \
$(OBJDIR)/AAL.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAO.o \
$(OBJDIR)/AAP.o \
$(OBJDIR)/AAQ.o \
$(OBJDIR)/AAR.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/AAU.o \
$(OBJDIR)/ABE.o \
$(OBJDIR)/ABF.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABJ.o \
$(OBJDIR)/ABK.o \
$(OBJDIR)/ABL.o \
$(OBJDIR)/ABM.o \
$(OBJDIR)/ABN.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABS.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ABY.o \
$(OBJDIR)/ABZ.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACB.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACE.o \
$(OBJDIR)/ACF.o \
$(OBJDIR)/ACI.o \
$(OBJDIR)/ACJ.o \
$(OBJDIR)/ACK.o \
$(OBJDIR)/ACL.o \
$(OBJDIR)/ACM.o \
$(OBJDIR)/ACN.o \
$(OBJDIR)/ACO.o \
$(OBJDIR)/ACP.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/ACR.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/openmp_base.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/$(SHM_FILE) \
$(OBJDIR)/store.o
$(OBJDIR)/ABR.o : $(SRCDIR)/ABR.f90 $(ABR_DEP)
ABR.o_name = ABR

ABS_depmod = dsp format_utils run_control slaarnAAG slaarnABG slaarnABT \
slaarnACA slaarnACC slaarnACH slaarnACS store
ABS_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABS.o : $(SRCDIR)/ABS.f90 $(ABS_DEP)
ABS.o_name = ABS

ABT_depmod = dsp format_utils run_control slaarnAAG
ABT_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/run_control.o
$(OBJDIR)/ABT.o : $(SRCDIR)/ABT.f90 $(ABT_DEP)
ABT.o_name = ABT

ABU_depmod = dsp file_utils format_utils parallel run_control slaarnABT \
slaarnACA store
ABU_DEP = \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABU.o : $(SRCDIR)/ABU.f90 $(ABU_DEP)
ABU.o_name = ABU

ABV_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnAAN slaarnABG slaarnABT slaarnACA slaarnACH store
ABV_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABV.o : $(SRCDIR)/ABV.f90 $(ABV_DEP)
ABV.o_name = ABV

ABW_depmod = dsp run_control slaarnAAG slaarnACH store
ABW_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABW.o : $(SRCDIR)/ABW.f90 $(ABW_DEP)
ABW.o_name = ABW

ABX_depmod = dsp file_utils format_utils parallel run_control singleton \
slaarnAAF slaarnAAG slaarnAAN slaarnAAS slaarnAAT slaarnABG slaarnABI \
slaarnABQ slaarnABT slaarnACA slaarnACC slaarnACH store
ABX_DEP = \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/singleton.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABX.o : $(SRCDIR)/ABX.f90 $(ABX_DEP)
ABX.o_name = ABX

ABY_depmod = dsp file_utils format_utils parallel run_control slaarnAAM \
slaarnAAN slaarnABG slaarnACH slaarnACP slaarnACS store
ABY_DEP = \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACP.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABY.o : $(SRCDIR)/ABY.f90 $(ABY_DEP)
ABY.o_name = ABY

ABZ_depmod = dsp file_utils format_utils parallel run_control slaarnAAI \
slaarnAAM slaarnAAN slaarnAAP slaarnAAT slaarnABG slaarnACH slaarnACP \
slaarnACQ slaarnACS store
ABZ_DEP = \
$(OBJDIR)/AAI.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAP.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACP.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ABZ.o : $(SRCDIR)/ABZ.f90 $(ABZ_DEP)
ABZ.o_name = ABZ

ACA_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnABG slaarnABT store
ACA_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACA.o : $(SRCDIR)/ACA.f90 $(ACA_DEP)
ACA.o_name = ACA

ACB_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnAAN slaarnABG slaarnABI slaarnABP slaarnABT slaarnACA store
ACB_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABI.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACB.o : $(SRCDIR)/ACB.f90 $(ACB_DEP)
ACB.o_name = ACB

ACC_depmod = dsp format_utils parallel run_control store
ACC_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACC.o : $(SRCDIR)/ACC.f90 $(ACC_DEP)
ACC.o_name = ACC

ACD_depmod = dsp format_utils run_control slaarnAAF slaarnAAG slaarnABG \
slaarnABJ store
ACD_DEP = \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABJ.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACD.o : $(SRCDIR)/ACD.f90 $(ACD_DEP)
ACD.o_name = ACD

ACE_depmod = dsp format_utils parallel run_control slaarnAAG slaarnABG \
slaarnABK slaarnABT slaarnACA slaarnACS store
ACE_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABK.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACE.o : $(SRCDIR)/ACE.f90 $(ACE_DEP)
ACE.o_name = ACE

ACF_depmod = dsp file_utils format_utils parallel run_control slaarnAAB \
slaarnAAM slaarnAAN slaarnABG slaarnABS slaarnABT slaarnACA slaarnACC \
slaarnACH slaarnACS store
ACF_DEP = \
$(OBJDIR)/AAB.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABS.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACF.o : $(SRCDIR)/ACF.f90 $(ACF_DEP)
ACF.o_name = ACF

ACG_depmod = dsp run_control slaarnACH slaarnACQ store
ACG_DEP = \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACG.o : $(SRCDIR)/ACG.f90 $(ACG_DEP)
ACG.o_name = ACG

ACH_depmod = dsp format_utils parallel run_control slaarnAAF slaarnAAN \
slaarnABG slaarnABT store
ACH_DEP = \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACH.o : $(SRCDIR)/ACH.f90 $(ACH_DEP)
ACH.o_name = ACH

ACI_depmod = dsp format_utils parallel run_control slaarnAAG slaarnABG \
store
ACI_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACI.o : $(SRCDIR)/ACI.f90 $(ACI_DEP)
ACI.o_name = ACI

ACJ_depmod = dsp format_utils parallel run_control slaarnAAD slaarnAAF \
slaarnAAG slaarnABG slaarnABP slaarnABT slaarnACH slaarnACQ store
ACJ_DEP = \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACJ.o : $(SRCDIR)/ACJ.f90 $(ACJ_DEP)
ACJ.o_name = ACJ

ACK_depmod = dsp format_utils parallel run_control slaarnAAG slaarnABG \
store
ACK_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACK.o : $(SRCDIR)/ACK.f90 $(ACK_DEP)
ACK.o_name = ACK

ACL_depmod = dsp file_utils format_utils parallel run_control slaarnAAG \
slaarnAAT slaarnABG slaarnABP slaarnABT store
ACL_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAT.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACL.o : $(SRCDIR)/ACL.f90 $(ACL_DEP)
ACL.o_name = ACL

ACM_depmod = casl dsp format_utils machine_constants parallel run_control \
slaarnAAA slaarnAAF slaarnAAG slaarnAAH slaarnAAM slaarnAAS slaarnABG \
slaarnABT slaarnACC slaarnACH slaarnACS store toms573
ACM_DEP = \
$(OBJDIR)/AAA.o \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAH.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/nl2sol.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACM.o : $(SRCDIR)/ACM.f90 $(ACM_DEP)
ACM.o_name = ACM

ACN_depmod = casl dsp format_utils parallel run_control slaarnAAF slaarnAAG \
slaarnAAH slaarnAAM slaarnABL slaarnABT slaarnACC slaarnACH slaarnACS store \
toms573
ACN_DEP = \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAH.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/ABL.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/nl2sol.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACN.o : $(SRCDIR)/ACN.f90 $(ACN_DEP)
ACN.o_name = ACN

ACO_depmod = dsp run_control slaarnABG slaarnABK slaarnABQ slaarnABS \
slaarnACA store
ACO_DEP = \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABK.o \
$(OBJDIR)/ABQ.o \
$(OBJDIR)/ABS.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACO.o : $(SRCDIR)/ACO.f90 $(ACO_DEP)
ACO.o_name = ACO

ACP_depmod = dsp file_utils format_utils parallel run_control slaarnAAA \
slaarnAAB slaarnAAD slaarnAAF slaarnAAG slaarnAAM slaarnAAN slaarnAAQ \
slaarnAAS slaarnABG slaarnABJ slaarnABK slaarnABT slaarnACA slaarnACC \
slaarnACD slaarnACE slaarnACH slaarnACR slaarnACS store
ACP_DEP = \
$(OBJDIR)/AAA.o \
$(OBJDIR)/AAB.o \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAF.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAM.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAQ.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABJ.o \
$(OBJDIR)/ABK.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACA.o \
$(OBJDIR)/ACC.o \
$(OBJDIR)/ACD.o \
$(OBJDIR)/ACE.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACR.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/file_utils.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACP.o : $(SRCDIR)/ACP.f90 $(ACP_DEP)
ACP.o_name = ACP

ACQ_depmod = casl dsp format_utils parallel run_control slaarnAAC slaarnAAE \
slaarnAAK slaarnAAN slaarnAAP slaarnAAR slaarnAAS slaarnAAU slaarnABG \
slaarnABP slaarnABT slaarnABU slaarnACB slaarnACI slaarnACK slaarnACL store
ACQ_DEP = \
$(OBJDIR)/AAC.o \
$(OBJDIR)/AAE.o \
$(OBJDIR)/AAK.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/AAP.o \
$(OBJDIR)/AAR.o \
$(OBJDIR)/AAS.o \
$(OBJDIR)/AAU.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABP.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ABU.o \
$(OBJDIR)/ACB.o \
$(OBJDIR)/ACI.o \
$(OBJDIR)/ACK.o \
$(OBJDIR)/ACL.o \
$(OBJDIR)/casl.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACQ.o : $(SRCDIR)/ACQ.f90 $(ACQ_DEP)
ACQ.o_name = ACQ

ACR_depmod = dsp format_utils parallel run_control slaarnAAG slaarnAAN \
slaarnABG slaarnABT slaarnACH slaarnACQ slaarnACS store
ACR_DEP = \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAN.o \
$(OBJDIR)/ABG.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACQ.o \
$(OBJDIR)/ACS.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACR.o : $(SRCDIR)/ACR.f90 $(ACR_DEP)
ACR.o_name = ACR

ACS_depmod = dsp run_control slaarnAAD slaarnAAG slaarnAAO slaarnABC \
slaarnABE slaarnABL slaarnABN slaarnABT slaarnABW slaarnACH slaarnACJ store
ACS_DEP = \
$(OBJDIR)/AAD.o \
$(OBJDIR)/AAG.o \
$(OBJDIR)/AAO.o \
$(OBJDIR)/ABC.o \
$(OBJDIR)/ABE.o \
$(OBJDIR)/ABL.o \
$(OBJDIR)/ABN.o \
$(OBJDIR)/ABT.o \
$(OBJDIR)/ABW.o \
$(OBJDIR)/ACH.o \
$(OBJDIR)/ACJ.o \
$(OBJDIR)/dsp.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/ACS.o : $(SRCDIR)/ACS.f90 $(ACS_DEP)
ACS.o_name = ACS

$(OBJDIR)/casl.o : $(SRCDIR)/casl.f90
casl.o_name = casl

$(OBJDIR)/comms_parallel_mpi1.o : $(SRCDIR)/comms_parallel_mpi1.f90
comms_parallel_mpi1.o_name = comms_parallel_mpi1

$(OBJDIR)/comms_parallel_mpi2.o : $(SRCDIR)/comms_parallel_mpi2.f90
comms_parallel_mpi2.o_name = comms_parallel_mpi2

comms_serial_depmod = dsp
COMMS_SERIAL_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/comms_serial.o : $(SRCDIR)/comms_serial.f90 $(COMMS_SERIAL_DEP)
comms_serial.o_name = comms_serial

$(OBJDIR)/dsp.o : $(SRCDIR)/dsp.f90
dsp.o_name = dsp

esdf_depmod = dsp esdf_key format_utils run_control
ESDF_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/esdf_key.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/run_control.o
$(OBJDIR)/esdf.o : $(SRCDIR)/esdf.f90 $(ESDF_DEP)
esdf.o_name = esdf

$(OBJDIR)/esdf_key.o : $(SRCDIR)/esdf_key.f90
esdf_key.o_name = esdf_key

file_utils_depmod = run_control store
FILE_UTILS_DEP = \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/file_utils.o : $(SRCDIR)/file_utils.f90 $(FILE_UTILS_DEP)
file_utils.o_name = file_utils

format_utils_depmod = dsp store
FORMAT_UTILS_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/store.o
$(OBJDIR)/format_utils.o : $(SRCDIR)/format_utils.f90 $(FORMAT_UTILS_DEP)
format_utils.o_name = format_utils

nl2sol_depmod = run_control
NL2SOL_DEP = \
$(OBJDIR)/run_control.o
$(OBJDIR)/nl2sol.o : $(SRCDIR)/nl2sol.f90 $(NL2SOL_DEP)
nl2sol.o_name = nl2sol

openmp_base_depmod = parallel
OPENMP_BASE_DEP = \
$(OBJDIR)/parallel.o
$(OBJDIR)/openmp_base.o : $(SRCDIR)/openmp_base.f90 $(OPENMP_BASE_DEP)
openmp_base.o_name = openmp_base

parallel_depmod = comms dsp format_utils store
PARALLEL_DEP = \
$(OBJDIR)/$(COMMS_FILE) \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/store.o
$(OBJDIR)/parallel.o : $(SRCDIR)/parallel.f90 $(PARALLEL_DEP)
parallel.o_name = parallel

run_control_depmod = dsp format_utils parallel store
RUN_CONTROL_DEP = \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/store.o
$(OBJDIR)/run_control.o : $(SRCDIR)/run_control.f90 $(RUN_CONTROL_DEP)
run_control.o_name = run_control

shalloc_nonsmp_depmod = dsp
SHALLOC_NONSMP_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/shalloc_nonsmp.o : $(SRCDIR)/shalloc_nonsmp.f90 $(SHALLOC_NONSMP_DEP)
shalloc_nonsmp.o_name = shalloc_nonsmp

shalloc_smp_depmod = comms dsp format_utils parallel run_control slaarnACC \
store
SHALLOC_SMP_DEP = \
$(OBJDIR)/ACC.o \
$(OBJDIR)/$(COMMS_FILE) \
$(OBJDIR)/dsp.o \
$(OBJDIR)/format_utils.o \
$(OBJDIR)/parallel.o \
$(OBJDIR)/run_control.o \
$(OBJDIR)/store.o
$(OBJDIR)/shalloc_smp.o : $(SRCDIR)/shalloc_smp.f90 $(SHALLOC_SMP_DEP)
shalloc_smp.o_name = shalloc_smp

$(OBJDIR)/singleton.o : $(SRCDIR)/singleton.f90
singleton.o_name = singleton

store_depmod = dsp
STORE_DEP = \
$(OBJDIR)/dsp.o
$(OBJDIR)/store.o : $(SRCDIR)/store.f90 $(STORE_DEP)
store.o_name = store
# End module dependencies

# How to compile things
$(BLAS_no):
	@echo "BLAS library missing" ;\
        echo "Type 'make' under the CASINO directory to build the entire\
 distribution." ;\
        exit 1

$(LAPACK_no):
	@echo "LAPACK library missing" ;\
        echo "Type 'make' under the CASINO directory to build the entire\
 distribution." ;\
        exit 1

$(VERSION_FILE):
	@cd $(SRCDIR)/.. && $(VERSION_SCRIPT)

$(OBJDIR)/main.o: $(VERSION_FILE)

$(OBJ):
	@$(MKOBJDIR) ;\
        action=`bash -c '$(BASH_PRETTY) ; compile_info F90 src\
         $($(@F)_name).f90'` ;\
        inaction=`bash -c '$(BASH_PRETTY) ; nocompile_info F90 src\
         $($(@F)_name).f90'` ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(F90WRAP) $($(@F)_name) $(SRCDIR)/$($(@F)_name).f90 "$(WINWRAP)\
        $(F90) -c $(FFLAGS) $(SRCDIR)/$($(@F)_name).f90 -I$(SRCDIR)/.."\
        --mods $($($(@F)_name)_mods) --deps $($($(@F)_name)_depmod) --silent\
        $(SILENT) --action "$$action" --inaction "$$inaction"\
        --modname_bug $(MODNAME_BUG)

$(OBJ_NOOPT):
	@$(MKOBJDIR) ;\
        action=`bash -c '$(BASH_PRETTY) ; compile_info F90-O0 src\
         $($(@F)_name).f90'` ;\
        inaction=`bash -c '$(BASH_PRETTY) ; nocompile_info F90-O0 src\
         $($(@F)_name).f90'` ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(F90WRAP) $($(@F)_name) $(SRCDIR)/$($(@F)_name).f90 "$(WINWRAP)\
        $(F90) -c $(FFLAGS0) $(SRCDIR)/$($(@F)_name).f90 -I$(SRCDIR)/.."\
        --mods $($($(@F)_name)_mods) --deps $($($(@F)_name)_depmod) --silent\
        $(SILENT) --action "$$action" --inaction "$$inaction"\
        --modname_bug $(MODNAME_BUG)

$(ETIME): $(ETIMESRC)
	@$(MKOBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info CC src etime.c' ;\
        $(ENVIRONMENT_COMMAND) ;\
        cd $(OBJDIR) && $(WINWRAP) $(CC) $(CFLAGS) $(CFLAGS_ETIME)\
        -c $(ETIMESRC) $(INCLUDE_DIR)

$(ALLOC): $(ALLOCSRC)
	@$(MKOBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info CC src alloc_shm.c' ;\
        $(ENVIRONMENT_COMMAND) ;\
        cd $(OBJDIR) && $(WINWRAP) $(CC) $(CFLAGS) $(CFLAGS_SHM)\
        -c $(ALLOCSRC) $(INCLUDE_DIR)

$(BINDIR)/$(EXECUTABLE): $(OBJ) $(OBJ_NOOPT) $(ETIME) $(ALLOC) $(BLAS) $(LAPACK)
	@$(MKBINDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info LDF90 src $(EXECUTABLE)' ;\
        $(ENVIRONMENT_COMMAND) ;\
        $(WINWRAP) $(LDF90) $(LDFLAGS) -o $(BINDIR)/$(EXECUTABLE) $(OBJ)\
        $(OBJ_NOOPT) $(ETIME) $(ALLOC) $(LDLIBS)

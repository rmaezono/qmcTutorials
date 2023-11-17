#------------------------------------------------------------------#
# LOUIS utils Makefile                                             #
#                                                                  #
# Allows easy maintenance and installation of all tools and        #
# utilities supplied with the LOUIS distribution. The Makefile     #
# will go through all the relevant utility directories, compiling  #
# and copying everything necessary into the bin directory          #
# LOUIS/bin which you should include in your path.                 #
#                                                                  #
# All system dependent stuff is derived from the include           #
# files in the zmakes directory.                                   #
#                                                                  #
# Targets:                                                         #
#  default: Make all utils and stick in relevant places.           #
#  clean  : Removes all object files.                              #
#  vclean : Removes all object and executable files.               #
#                                                                  #
# MDT 2010                                                         #
#------------------------------------------------------------------#
# Basic stuff
UTILDIR = $(PWD)
INCBASE = $(UTILDIR)/../arch/make.inc
BINBASE = $(UTILDIR)/../bin
OBJBASE = $(UTILDIR)/zlib
.SILENT:

start: default
.PHONY: start

# System-independent variables
BINDIR     = $(BINBASE)/utils/$(LOUIS_ARCH)
OBJDIR     = $(OBJBASE)/$(LOUIS_ARCH)
MKBINDIR   = [ -d $(BINDIR) ] || mkdir -p $(BINDIR)
MKDUMMY_nopath = { $(MKBINDIR) ; echo "[DUMMY] $(@F)" ;\
                 rm -f $(BINDIR)/$@ 2> /dev/null ;\
                 touch $(UTILDIR)/dummy_script/dummy_script ;\
                 ln -s $(UTILDIR)/dummy_script/dummy_script $(BINDIR)/$@ ; }

# Defaults for system-dependent variables
LDF90 = $(F90)

# System-dependent include file
include $(INCBASE)/$(LOUIS_ARCH).inc
include $(INCBASE)/users/$(LOUIS_ARCH)/$(LOUIS_ID).inc

# System-dependent variables from data in include files
FFLAGS   = $(FFLAGS_opt) $(INCLUDE_DIR)
LDFLAGS  = $(LDFLAGS_opt)

# Main target
.PHONY: default
default: info wrap

# Info
.PHONY: info
info:
	@echo "LOUIS utils MAKE information:" ;\
        echo " * LN = ln -s" ;\
        echo " * F90 = $(F90) $(FFLAGS)" ;\
        echo " * CC = $(CC) $(CCFLAGS)" ;\
        [ -z "$(CXX)" ] || echo " * CXX = $(CXX) $(CXXFLAGS)" ;\
        echo " * LDF90 = $(LDF90) $(LDFLAGS) $(LDLIBS)" ;\
        echo "---"

# Types of utilities according to how they are made
UTILS = $(SCRIPT_UTILS) $(SINGLE_F90_UTILS) $(SUBMAKE_UTILS)

# Utilities which are scripts
SCRIPT_UTILS = $(SCRIPT_UTILS_all)
# For all LOUIS_ARCHs
SCRIPT_UTILS_all = \
$(BINDIR)/change_inputs \
$(BINDIR)/louishelp \
$(BINDIR)/louis_help \
$(BINDIR)/plot_louis \
$(BINDIR)/plot2d \
$(BINDIR)/plot3d \
$(BINDIR)/runlouis

# Dependencies
$(BINDIR)/change_inputs : change_inputs/change_inputs
$(BINDIR)/louishelp : help/louishelp
$(BINDIR)/louis_help : help/louis_help
$(BINDIR)/plot_louis : plot/plot_louis
$(BINDIR)/plot2d : plot_other/plot2d
$(BINDIR)/plot3d : plot_other/plot3d
$(BINDIR)/runlouis : runlouis/runlouis

# How to make them
LINKSCRIPT = { $(MKBINDIR) ; echo "[LN] $(@F)" ; rm -f $@ 2> /dev/null ;\
             ln -s $(UTILDIR)/$? $@ ; }
$(SCRIPT_UTILS_all):
	@$(LINKSCRIPT)

# Utilities with an independent Makefile. NB, the dependencies are
# directly handled by the relevant Makefiles.
SUBMAKE_UTILS = $(SUBMAKE_UTILS_all) $(SUBMAKE_UTILS_cxx)
# For all LOUIS_ARCHs
SUBMAKE_UTILS_all = louis_help_helper
# For some LOUIS_ARCHs only
SUBMAKE_UTILS_cxx  = 
.PHONY: $(SUBMAKE_UTILS)

# Directory and makefile list
MAKEDIR_louis_help_helper = help

# How to make them
# NB, GNU make requires '$(MAKE)' to be explicitly specified in the
# rule in order for the jobserver to be available to sub-makes. So the
# following statement:
#  { $(MKBINDIR) ; echo "[MAKE] $(@)" ; cd $(MAKEDIR_$(@F)) &&\
#    $(MAKE) 'PWD=$(UTILDIR)/$(MAKEDIR_$(@F))' || exit 1 ; }
# is split around $(MAKE):
SUBMAKE1 = { $(MKBINDIR) ; echo "[MAKE] $@" ; cd $(MAKEDIR_$@) &&
SUBMAKE2 = 'PWD=$(UTILDIR)/$(MAKEDIR_$@)' || exit 1 ; }
$(SUBMAKE_UTILS_all):
	@$(SUBMAKE1) $(MAKE) $(SUBMAKE2)
$(SUBMAKE_UTILS_cxx):
	@[ ! -z "$(CXX)" ] && $(SUBMAKE1) $(MAKE) $(SUBMAKE2)\
        || $(MKDUMMY_nopath)

# Utilities that only require compiling a single F90 file
SINGLE_F90_UTILS = $(SINGLE_F90_UTILS_all)
# For all LOUIS_ARCHs
SINGLE_F90_UTILS_all = 
#$(BINDIR)/example_fortran_binary

# Dependencies
#$(BINDIR)/examlpe_fortran_binary: arse/example_fortran_binary.f90

# How to make them
F90_SINGLE = { echo "[F90] $(@F)" ; rm -f $@ 2> /dev/null ;\
             [ -d $(OBJDIR)/$(@F) ] || mkdir -p $(OBJDIR)/$(@F) ;\
             cd $(OBJDIR)/$(@F) && $(F90) $(FFLAGS) -o $@ $(UTILDIR)/$? ||\
             exit 1 ; }
$(SINGLE_F90_UTILS_all):
	@$(F90_SINGLE)

# Target for creating/updating wrapper scripts
.PHONY: wrap utils build

utils: $(UTILS)

build: $(SINGLE_F90_UTILS) $(SUBMAKE_UTILS)

wrap: utils
	@echo "[WRAP] utils" ; $(MKBINDIR) ;\
        cd $(BINBASE) && $(UTILDIR)/wrappers/make_wrappers

# Target for deleting object files
.PHONY: clean
clean:
	@echo "Cleaning object files [$(LOUIS_ARCH)]." ;\
        rm -rf $(OBJDIR) ; :

# Target for deleting object files and compiled utils
.PHONY: vclean
vclean:
	@echo "Cleaning object files and compiled utilities [$(LOUIS_ARCH)]." ;\
        rm -rf $(OBJDIR) $(BINDIR) ; :

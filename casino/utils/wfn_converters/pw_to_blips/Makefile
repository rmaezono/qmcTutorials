#------------------------------------------------------------#
# Makefile for the blip converter.                           #
#                                                            #
# MDT 9.2003                                                 #
#------------------------------------------------------------#
# Basic stuff
UTILEXE = blip
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
$(OBJDIR)/blip.o \
$(OBJDIR)/helpers.o \
$(OBJDIR)/rng.o \
$(OBJDIR)/singleton.o
# End object-file list

# Module dependencies
BLIP_DEP = \
$(OBJDIR)/helpers.o \
$(OBJDIR)/rng.o \
$(OBJDIR)/singleton.o
$(OBJDIR)/blip.o : $(SRCDIR)/blip.f90 $(SRCDIR)/blipp.f90 $(SRCDIR)/blipk.f90 $(BLIP_DEP)
blip.o_name = blip

$(OBJDIR)/helpers.o : $(SRCDIR)/helpers.f90
helpers.o_name = helpers

$(OBJDIR)/rng.o : $(SRCDIR)/rng.f90
rng.o_name = rng

$(OBJDIR)/singleton.o : $(SRCDIR)/singleton.f90
singleton.o_name = singleton
# End module dependencies

# How to make the objects
$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info F90 utils $(UTILEXE)\
        $($(@F)_name).f90' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)\
        && $(F90_UTILS) -c $(FFLAGS_UTILS) $(SRCDIR)/$($(@F)_name).f90\
        -I$(SRCDIR)
$(BINDIR)/$(UTILEXE) : $(OBJ)
	@[ -d $(BINDIR) ] || mkdir -p $(BINDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info LDF90 utils $(UTILEXE)' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ;\
        $(LDF90_UTILS) $(LDFLAGS_UTILS) -o $(@) $(OBJ) $(LDLIBS_UTILS)

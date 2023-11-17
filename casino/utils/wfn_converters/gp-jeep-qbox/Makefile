#------------------------------------------------------------#
# Makefile for the jeep_to_pwfn program.                     #
#                                                            #
# MDT/NDD 10.2003                                            #
#------------------------------------------------------------#
# Basic stuff
UTILEXE = jeep_to_pwfn
SRCDIR  = $(PWD)
UTILDIR = $(SRCDIR)/../..
TOOLDIR = $(UTILDIR)/../build-tools
INCBASE = $(UTILDIR)/../arch/data
BINDIR  = $(UTILDIR)/../bin_qmc/utils/$(CASINO_ARCH)
OBJDIR  = $(UTILDIR)/zlib/$(CASINO_ARCH)/$(UTILEXE)
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

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

# Use the same compilers as for src, unless overrides defined
UNATIVE                    = $(UNATIVE_$(UTILS_MODE))
ENVIRONMENT_COMMAND_UTILS  = $(ENVIRONMENT_COMMAND$(UNATIVE))
F90_UTILS                  = $(F90$(UNATIVE))
FFLAGS_opt_UTILS           = $(FFLAGS_opt$(UNATIVE))
FFLAGS_all_UTILS           = $(FFLAGS_all$(UNATIVE))
INCLUDE_DIR_UTILS          = $(INCLUDE_DIR$(UNATIVE))
LDF90_UTILS                = $(LDF90$(UNATIVE))
LDFLAGS_opt_UTILS          = $(LDFLAGS_opt$(UNATIVE))
LDFLAGS_all_UTILS          = $(LDFLAGS_all$(UNATIVE))
LIB_PATH_UTILS             = $(LIB_PATH$(UNATIVE))
LDLIBS_opt_UTILS           = $(LDLIBS_opt$(UNATIVE))
LDLIBS_all_UTILS           = $(LDLIBS_all$(UNATIVE))
CC_UTILS                   = $(CC$(UNATIVE))
CFLAGS_opt_UTILS           = $(CFLAGS_opt$(UNATIVE))
CFLAGS_all_UTILS           = $(CFLAGS_all$(UNATIVE))
CFLAGS_F90_INTERFACE_UTILS = $(CFLAGS_F90_INTERFACE$(UNATIVE))

# System-dependent variables from data in include files
FFLAGS_UTILS  = $(FFLAGS_opt_UTILS) $(FFLAGS_all_UTILS) $(INCLUDE_DIR_UTILS)
LDFLAGS_UTILS = $(LDFLAGS_opt_UTILS) $(LDFLAGS_all_UTILS) $(LIB_PATH_UTILS)
CFLAGS_UTILS  = $(CFLAGS_opt_UTILS) $(CFLAGS_all_UTILS)\
 $(CFLAGS_F90_INTERFACE_UTILS)
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
$(OBJDIR)/gen_gvec.o \
$(OBJDIR)/jeep_to_pwfn.o
# End object-file list

# Module dependencies
$(OBJDIR)/gen_gvec.o : $(SRCDIR)/gen_gvec.f90
gen_gvec.o_name = gen_gvec

JEEP_TO_PWFN_DEP = \
$(OBJDIR)/gen_gvec.o
$(OBJDIR)/jeep_to_pwfn.o : $(SRCDIR)/jeep_to_pwfn.f90 $(JEEP_TO_PWFN_DEP)
jeep_to_pwfn.o_name = jeep_to_pwfn
# End module dependencies

# How to make the objects
$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info F90 utils $(UTILEXE)\
        $($(@F)_name).f90' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)\
        && $(F90_UTILS) -c $(FFLAGS_UTILS) $(SRCDIR)/$($(@F)_name).f90
$(OBJDIR)/jeep_io.o : $(SRCDIR)/jeep_io.c
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info CC utils $(UTILEXE)\
        jeep_io.c' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)\
        && $(CC_UTILS) $(CFLAGS_UTILS) -c $(SRCDIR)/jeep_io.c
$(BINDIR)/$(UTILEXE) : $(OBJDIR)/jeep_io.o $(OBJ)
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info LDF90 utils $(UTILEXE)' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ;\
        $(LDF90_UTILS) $(LDFLAGS_UTILS) -o $(@) $(OBJDIR)/jeep_io.o $(OBJ)\
        $(LDLIBS_UTILS)

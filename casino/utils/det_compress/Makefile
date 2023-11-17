#----------------------------------------#
# Makefile for the det_compress utility. #
#                                        #
# PLR 01.2013                            #
#----------------------------------------#
# Basic stuff
UTILEXE = det_compress
SRCDIR  = $(PWD)
UTILDIR = $(SRCDIR)/..
TOOLDIR = $(UTILDIR)/../build-tools
INCBASE = $(UTILDIR)/../arch/data
LIBDIR  = $(UTILDIR)/../lib/zlib/$(CASINO_ARCH)
HDRDIR  = $(UTILDIR)/../lib
BINDIR  = $(UTILDIR)/../bin_qmc/utils/$(CASINO_ARCH)
OBJDIR  = $(UTILDIR)/zlib/$(CASINO_ARCH)/$(UTILEXE)
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

# System-independent variables
UTILS_MODE          =
UNATIVE_            =
UNATIVE_native      = _NATIVE
LPSOLVE_no          = $(LIBDIR)/liblpsolve_native.a
LPSOLVE_yes         =
COLAMD_no           = $(LIBDIR)/libcolamd_native.a
COLAMD_yes          =
LDLPSOLVE_no        = -L$(LIBDIR) -llpsolve_native -ldl
LDLPSOLVE_yes       =
LDCOLAMD_no         = -L$(LIBDIR) -lcolamd_native -ldl
LDCOLAMD_yes        =
INCLUDE_LPSOLVE_no  = -I$(HDRDIR)/lpsolve
INCLUDE_LPSOLVE_yes =
INCLUDE_COLAMD_no   = -I$(HDRDIR)/colamd
INCLUDE_COLAMD_yes  =
ETIME_yes           = $(OBJDIR)/etime.o
ETIME_no            =
ETIMESRC_yes        = $(SRCDIR)/etime.c
ETIMESRC_no         =

# Defaults for system-dependent variables
LDF90                      = $(F90)
LDFLAGS_opt                = $(FFLAGS_opt)
LDFLAGS_all                = $(FFLAGS_all)
ENVIRONMENT_COMMAND        = :
HAVE_LPSOLVE_NATIVE        = no
HAVE_COLAMD_NATIVE         = no
LDF90_NATIVE               = $(F90_NATIVE)
LDFLAGS_opt_NATIVE         = $(FFLAGS_opt_NATIVE)
LDFLAGS_all_NATIVE         = $(FFLAGS_all_NATIVE)
ENVIRONMENT_COMMAND_NATIVE = :
NEED_ETIME                 = no
HAVE_LPSOLVE               = no
HAVE_COLAMD                = no
ISOVAR_VERSION             = 1

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
HAVE_LPSOLVE_UTILS         = $(HAVE_LPSOLVE$(UNATIVE))
HAVE_COLAMD_UTILS          = $(HAVE_COLAMD$(UNATIVE))
CC_UTILS                   = $(CC$(UNATIVE))
CFLAGS_opt_UTILS           = $(CFLAGS_opt$(UNATIVE))
CFLAGS_all_UTILS           = $(CFLAGS_all$(UNATIVE))
CFLAGS_F90_INTERFACE_UTILS = $(CFLAGS_F90_INTERFACE$(UNATIVE))

# System-dependent variables from data in include files
LPSOLVE_UTILS         = $(LPSOLVE_$(HAVE_LPSOLVE_UTILS))
LDLPSOLVE_UTILS       = $(LDLPSOLVE_$(HAVE_LPSOLVE_UTILS))
INCLUDE_LPSOLVE_UTILS = $(INCLUDE_LPSOLVE_$(HAVE_LPSOLVE_UTILS))
COLAMD_UTILS          = $(COLAMD_$(HAVE_COLAMD_UTILS))
LDCOLAMD_UTILS        = $(LDCOLAMD_$(HAVE_COLAMD_UTILS))
INCLUDE_COLAMD_UTILS  = $(INCLUDE_COLAMD_$(HAVE_COLAMD_UTILS))
FFLAGS_UTILS          = $(FFLAGS_opt_UTILS) $(FFLAGS_all_UTILS)\
 $(INCLUDE_DIR_UTILS)
LDFLAGS_UTILS         = $(LDFLAGS_opt_UTILS) $(LDFLAGS_all_UTILS)\
 $(LIB_PATH_UTILS)
CFLAGS_UTILS          = $(CFLAGS_opt_UTILS) $(CFLAGS_all_UTILS)\
 $(CFLAGS_F90_INTERFACE_UTILS)
LDLIBS_UTILS          = $(LDLPSOLVE_UTILS) $(LDCOLAMD_UTILS)\
 $(LDLIBS_opt_UTILS) $(LDLIBS_all_UTILS)

# Main target
.PHONY: default
default: $(BINDIR)/$(UTILEXE)

# 'Clean' target
.PHONY: clean
clean:
	@rm -rf $(OBJDIR) $(BINDIR)/$(UTILEXE)

# Object-file list
OBJ = \
$(OBJDIR)/casl.o \
$(OBJDIR)/det_compress.o \
$(OBJDIR)/tools.o
# End object-file list

# Module dependencies
$(OBJDIR)/casl.o : $(SRCDIR)/casl.f90
casl.o_name = casl

DET_COMPRESS_DEP = \
$(OBJDIR)/casl.o \
$(OBJDIR)/tools.o
$(OBJDIR)/det_compress.o : $(SRCDIR)/det_compress.f90 $(DET_COMPRESS_DEP)
det_compress.o_name = det_compress

$(OBJDIR)/tools.o : $(SRCDIR)/tools.f90
tools.o_name = tools
# End module dependencies

# How to make the objects
$(LPSOLVE_no):
	@echo "LP_SOLVE library missing" ;\
        echo "Type 'make' under the CASINO directory to build the entire\
 distribution." ;\
        exit 1
$(COLAMD_no):
	@echo "COLAMD library missing" ;\
        echo "Type 'make' under the CASINO directory to build the entire\
 distribution." ;\
        exit 1
$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info F90 utils $(UTILEXE)\
        $($(@F)_name).f90' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)\
        && $(F90_UTILS) -c $(FFLAGS_UTILS) $(SRCDIR)/$($(@F)_name).f90
$(OBJDIR)/lp2fortran.o : $(SRCDIR)/lp2fortran.c
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info CC utils $(UTILEXE)\
        lp2fortran.c' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)\
        && $(CC_UTILS) $(CFLAGS_UTILS) -c $(SRCDIR)/lp2fortran.c\
        $(INCLUDE_LPSOLVE_UTILS) $(INCLUDE_COLAMD_UTILS)
$(BINDIR)/$(UTILEXE) : $(OBJ) $(OBJDIR)/lp2fortran.o $(LPSOLVE_UTILS)\
  $(COLAMD_UTILS)
	@[ -d $(BINDIR) ] || mkdir -p $(BINDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info LDF90 utils $(UTILEXE)' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ;\
        $(LDF90_UTILS) $(LDFLAGS_UTILS) -o $(@) $(OBJ) $(OBJDIR)/lp2fortran.o\
        $(LDLIBS_UTILS)

#-----------------------------------------------#
# Makefile to create static a library for BLAS  #
# (levels 1, 2 and 3). BLAS as distributed with #
# LAPACK 3.2, only double-precision routines.   #
#                                               #
# PLR 12.2008                                   #
#-----------------------------------------------#
# Basic stuff
SHELL          = /bin/sh
BLAS_SRCDIR    = $(PWD)
INCBASE        = $(BLAS_SRCDIR)/../../arch/data
OBJBASE        = $(BLAS_SRCDIR)/../zlib
TOOLDIR        = $(BLAS_SRCDIR)/../../build-tools
BLAS_AR        = libblas.a
BLAS_AR_NATIVE = libblas_native.a
SILENT         = 0
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

# System-independent variables
LIBDIR          = $(OBJBASE)/$(CASINO_ARCH)
OBJDIR          = $(LIBDIR)/blas
MKOBJDIR        = [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
BLAS_yes        =
BLAS_no         = $(LIBDIR)/$(BLAS_AR)
OBJDIR_NATIVE   = $(LIBDIR)/blas_native
MKOBJDIR_NATIVE = [ -d $(OBJDIR_NATIVE) ] || mkdir -p $(OBJDIR_NATIVE)
BLAS_NATIVE_yes =
BLAS_NATIVE_no  = $(LIBDIR)/$(BLAS_AR_NATIVE)

# Defaults for system-dependent variables
AR                         = ar crs
REDIR_0                    =
REDIR_1                    = >& /dev/null
FFLAGS_libs                = $(FFLAGS_opt) $(FFLAGS_all)
FFLAGS0_libs               = -O0 $(FFLAGS_all)
HAVE_BLAS                  = no
ENVIRONMENT_COMMAND        = :
FFLAGS_libs_NATIVE         = $(FFLAGS_opt_NATIVE) $(FFLAGS_all_NATIVE)
FFLAGS0_libs_NATIVE        = -O0 $(FFLAGS_all_NATIVE)
HAVE_BLAS_NATIVE           = no
ENVIRONMENT_COMMAND_NATIVE = :

# System-dependent include files
include $(INCBASE)/$(CASINO_ARCH).arch

# System-dependent variables from data in include files
REDIR                   = $(REDIR_$(SILENT))
BLAS                    = $(BLAS_$(HAVE_BLAS))
HAVE_BLAS_NATIVE_native = $(HAVE_BLAS_NATIVE)
HAVE_BLAS_NATIVE_       = $(HAVE_BLAS)
BLAS_NATIVE             = $(BLAS_NATIVE_$(HAVE_BLAS_NATIVE_$(UTILS_MODE)))

# Main target
.PHONY: default
default: info $(BLAS) $(BLAS_NATIVE)

# Info
.PHONY: info
info:
	@bash -c '$(BASH_PRETTY) ;\
        info_title "CASINO BLAS MAKE information" ;\
        if [ "$(UTILS_MODE)" != native ] ; then\
         if [ "$(HAVE_BLAS)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
          info_equal "F90" "$(F90) -c $(FFLAGS_libs)" ;\
          info_equal "F90-O0" "$(F90) -c $(FFLAGS0_libs)" ;\
          info_equal "AR" "$(AR)" ;\
         else\
          info_line "Using external BLAS library." ;\
         fi ;\
        else\
         if [ "$(HAVE_BLAS)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
          info_equal "F90" "$(F90) -c $(FFLAGS_libs)" ;\
          info_equal "F90-O0" "$(F90) -c $(FFLAGS0_libs)" ;\
          info_equal "AR" "$(AR)" ;\
         else\
          info_line "Using external BLAS library for main binary." ;\
         fi ;\
         if [ "$(HAVE_BLAS_NATIVE)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND_NATIVE)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND_NATIVE"\
           "$(ENVIRONMENT_COMMAND_NATIVE)" ;\
          info_equal "F90n" "$(F90_NATIVE) -c $(FFLAGS_libs_NATIVE)" ;\
          info_equal "F90n-O0" "$(F90_NATIVE) -c $(FFLAGS0_libs_NATIVE)" ;\
          info_equal "ARn" "$(AR)" ;\
         else\
          info_line "Using external BLAS library for utilities." ;\
         fi ;\
        fi ;\
        info_end'

# Target for deleting object files
.PHONY: clean
clean:
	@echo "Cleaning object files [$(CASINO_ARCH)]." ;\
        rm -rf $(OBJDIR) $(OBJDIR_NATIVE)

# Target for deleting object files
.PHONY: vclean
vclean:
	@echo "Cleaning object files and libraries [$(CASINO_ARCH)]." ;\
        rm -rf $(OBJDIR) $(OBJDIR_NATIVE) $(LIBDIR)/$(BLAS_AR)\
        $(LIBDIR)/$(BLAS_AR_NATIVE)

# Object-file list
OBJ_OPT = $(OBJDIR)/daxpy.o $(OBJDIR)/dcabs1.o $(OBJDIR)/dcopy.o \
 $(OBJDIR)/ddot.o $(OBJDIR)/dgemm.o $(OBJDIR)/dgemv.o $(OBJDIR)/dger.o \
 $(OBJDIR)/dnrm2.o $(OBJDIR)/drot.o $(OBJDIR)/dscal.o $(OBJDIR)/dswap.o \
 $(OBJDIR)/dsymm.o $(OBJDIR)/dsymv.o $(OBJDIR)/dtrmm.o $(OBJDIR)/dtrmv.o \
 $(OBJDIR)/dtrsm.o $(OBJDIR)/dzasum.o $(OBJDIR)/dznrm2.o $(OBJDIR)/idamax.o \
 $(OBJDIR)/izamax.o $(OBJDIR)/lsame.o $(OBJDIR)/xerbla.o $(OBJDIR)/zaxpy.o \
 $(OBJDIR)/zcopy.o $(OBJDIR)/zdotc.o $(OBJDIR)/zdotu.o $(OBJDIR)/zdscal.o \
 $(OBJDIR)/zgemm.o $(OBJDIR)/zgemv.o $(OBJDIR)/zgerc.o $(OBJDIR)/zgeru.o \
 $(OBJDIR)/zhemm.o $(OBJDIR)/zhemv.o $(OBJDIR)/zher2.o $(OBJDIR)/zher2k.o \
 $(OBJDIR)/zherk.o $(OBJDIR)/zscal.o $(OBJDIR)/zswap.o $(OBJDIR)/ztrmm.o \
 $(OBJDIR)/ztrmv.o $(OBJDIR)/ztrsm.o $(OBJDIR)/ztrsv.o
OBJ_NOOPT =
OBJ_OPT_NATIVE = $(OBJDIR_NATIVE)/daxpy.o $(OBJDIR_NATIVE)/dcabs1.o \
 $(OBJDIR_NATIVE)/dcopy.o $(OBJDIR_NATIVE)/ddot.o $(OBJDIR_NATIVE)/dgemm.o \
 $(OBJDIR_NATIVE)/dgemv.o $(OBJDIR_NATIVE)/dger.o $(OBJDIR_NATIVE)/dnrm2.o \
 $(OBJDIR_NATIVE)/drot.o $(OBJDIR_NATIVE)/dscal.o $(OBJDIR_NATIVE)/dswap.o \
 $(OBJDIR_NATIVE)/dsymm.o $(OBJDIR_NATIVE)/dsymv.o $(OBJDIR_NATIVE)/dtrmm.o \
 $(OBJDIR_NATIVE)/dtrmv.o $(OBJDIR_NATIVE)/dtrsm.o $(OBJDIR_NATIVE)/dzasum.o \
 $(OBJDIR_NATIVE)/dznrm2.o $(OBJDIR_NATIVE)/idamax.o $(OBJDIR_NATIVE)/izamax.o \
 $(OBJDIR_NATIVE)/lsame.o $(OBJDIR_NATIVE)/xerbla.o $(OBJDIR_NATIVE)/zaxpy.o \
 $(OBJDIR_NATIVE)/zcopy.o $(OBJDIR_NATIVE)/zdotc.o $(OBJDIR_NATIVE)/zdotu.o \
 $(OBJDIR_NATIVE)/zdscal.o $(OBJDIR_NATIVE)/zgemm.o $(OBJDIR_NATIVE)/zgemv.o \
 $(OBJDIR_NATIVE)/zgerc.o $(OBJDIR_NATIVE)/zgeru.o $(OBJDIR_NATIVE)/zhemm.o \
 $(OBJDIR_NATIVE)/zhemv.o $(OBJDIR_NATIVE)/zher2.o $(OBJDIR_NATIVE)/zher2k.o \
 $(OBJDIR_NATIVE)/zherk.o $(OBJDIR_NATIVE)/zscal.o $(OBJDIR_NATIVE)/zswap.o \
 $(OBJDIR_NATIVE)/ztrmm.o $(OBJDIR_NATIVE)/ztrmv.o $(OBJDIR_NATIVE)/ztrsm.o \
 $(OBJDIR_NATIVE)/ztrsv.o
OBJ_NOOPT_NATIVE =
# End object-file list

# All objects
OBJ               = $(OBJ_OPT) $(OBJ_NOOPT)
OBJ_NATIVE_       =
OBJ_NATIVE_native = $(OBJ_OPT_NATIVE) $(OBJ_NOOPT_NATIVE)
OBJ_NATIVE        = $(OBJ_NATIVE_$(UTILS_MODE))

# How to compile things
$(BLAS): $(OBJ)
	@bash -c '$(MKOBJDIR) ; a=$(@F) ; b=$${a%.a} ;\
        $(BASH_PRETTY) ; compile_info AR lib BLAS ;\
        rm -f $(BLAS) ; $(AR) $(BLAS) $(OBJ) $(REDIR)'

$(OBJ_OPT):
	@bash -c '$(MKOBJDIR) ; o=$(@F) ; b=$${o%.o} ; f=$$b.f ;\
        $(BASH_PRETTY) ; compile_info F90 lib BLAS $$f ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(F90) $(FFLAGS_libs) -c $(BLAS_SRCDIR)/$$f -o $@ $(REDIR)'

$(OBJ_NOOPT):
	@bash -c '$(MKOBJDIR) ; o=$(@F) ; b=$${o%.o} ; f=$$b.f ;\
        $(BASH_PRETTY) ; compile_info F90-O0 lib BLAS $$f ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(F90) $(FFLAGS0_libs) -c $(BLAS_SRCDIR)/$$f -o $@ $(REDIR)'

$(BLAS_NATIVE): $(OBJ_NATIVE)
	@bash -c '$(BASH_PRETTY) ;\
        if [ "$(UTILS_MODE)" = native ] ; then \
         $(MKOBJDIR_NATIVE) ; a=$(@F) ; b=$${a%.a} ;\
         compile_info ARn lib BLAS_NATIVE ;\
         rm -f $(BLAS_NATIVE) ; $(AR) $(BLAS_NATIVE) $(OBJ_NATIVE) $(REDIR) ;\
        else \
         compile_info LN lib BLAS_NATIVE ;\
         rm -f $(BLAS_NATIVE) ; ln -s $(LIBDIR)/$(BLAS_AR) $(BLAS_NATIVE) ;\
        fi'

$(OBJ_OPT_NATIVE):
	@bash -c '$(MKOBJDIR_NATIVE) ; o=$(@F) ; b=$${o%.o} ; f=$$b.f ;\
        $(BASH_PRETTY) ; compile_info F90n lib BLAS_NATIVE $$f ;\
        $(ENVIRONMENT_COMMAND_NATIVE) ; cd $(OBJDIR_NATIVE) &&\
        $(F90_NATIVE) $(FFLAGS_libs_NATIVE) -c $(BLAS_SRCDIR)/$$f\
        -o $@ $(REDIR)'

$(OBJ_NOOPT_NATIVE):
	@bash -c '$(MKOBJDIR_NATIVE) ; o=$(@F) ; b=$${o%.o} ; f=$$b.f ;\
        $(BASH_PRETTY) ; compile_info F90n-O0 lib BLAS_NATIVE $$f ;\
        $(ENVIRONMENT_COMMAND_NATIVE) ; cd $(OBJDIR_NATIVE) &&\
        $(F90_NATIVE) $(FFLAGS0_libs_NATIVE) -c $(BLAS_SRCDIR)/$$f\
        -o $@ $(REDIR)'

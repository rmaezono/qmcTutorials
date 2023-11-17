#----------------------------------------------------------#
# Makefile to create a static library for LAPACK (v3.1.1). #
# LAPACK as distributed from NETLIB, only double-precision #
# routines.                                                #
#                                                          #
# PLR 12.2008                                              #
#----------------------------------------------------------#
# Basic stuff
SHELL            = /bin/sh
LAPACK_SRCDIR    = $(PWD)
INCBASE          = $(LAPACK_SRCDIR)/../../arch/data
OBJBASE          = $(LAPACK_SRCDIR)/../zlib
TOOLDIR          = $(LAPACK_SRCDIR)/../../build-tools
LAPACK_AR        = liblapack.a
LAPACK_AR_NATIVE = liblapack_native.a
SILENT           = 0
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

# System-independent variables
LIBDIR            = $(OBJBASE)/$(CASINO_ARCH)
OBJDIR            = $(LIBDIR)/lapack
MKOBJDIR          = [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
LAPACK_yes        =
LAPACK_no         = $(LIBDIR)/$(LAPACK_AR)
OBJDIR_NATIVE     = $(LIBDIR)/lapack_native
MKOBJDIR_NATIVE   = [ -d $(OBJDIR_NATIVE) ] || mkdir -p $(OBJDIR_NATIVE)
LAPACK_NATIVE_yes =
LAPACK_NATIVE_no  = $(LIBDIR)/$(LAPACK_AR_NATIVE)

# Defaults for system-dependent variables
AR                         = ar crs
REDIR_0                    =
REDIR_1                    = >& /dev/null
FFLAGS_libs                = $(FFLAGS_opt) $(FFLAGS_all)
FFLAGS0_libs               = -O0 $(FFLAGS_all)
HAVE_LAPACK                = no
ENVIRONMENT_COMMAND        = :
FFLAGS_libs_NATIVE         = $(FFLAGS_opt_NATIVE) $(FFLAGS_all_NATIVE)
FFLAGS0_libs_NATIVE        = -O0 $(FFLAGS_all_NATIVE)
HAVE_LAPACK_NATIVE         = no
ENVIRONMENT_COMMAND_NATIVE = :

# System-dependent include files
include $(INCBASE)/$(CASINO_ARCH).arch

# System-dependent variables from data in include files
REDIR                     = $(REDIR_$(SILENT))
LAPACK                    = $(LAPACK_$(HAVE_LAPACK))
HAVE_LAPACK_NATIVE_native = $(HAVE_LAPACK_NATIVE)
HAVE_LAPACK_NATIVE_       = $(HAVE_LAPACK)
LAPACK_NATIVE             = $(LAPACK_NATIVE_$(HAVE_LAPACK_NATIVE_$(UTILS_MODE)))

# Main target
.PHONY: default
default: info $(LAPACK) $(LAPACK_NATIVE)

# Info
.PHONY: info
info:
	@bash -c '$(BASH_PRETTY) ;\
        info_title "CASINO LAPACK MAKE information" ;\
        if [ "$(UTILS_MODE)" != native ] ; then\
         if [ "$(HAVE_LAPACK)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
          info_equal "F90" "$(F90) -c $(FFLAGS_libs)" ;\
          info_equal "F90-O0" "$(F90) -c $(FFLAGS0_libs)" ;\
          info_equal "AR" "$(AR)" ;\
         else\
          info_line "Using external LAPACK library." ;\
         fi ;\
        else\
         if [ "$(HAVE_LAPACK)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
          info_equal "F90" "$(F90) -c $(FFLAGS_libs)" ;\
          info_equal "F90-O0" "$(F90) -c $(FFLAGS0_libs)" ;\
          info_equal "AR" "$(AR)" ;\
         else\
          info_line "Using external LAPACK library for main binary." ;\
         fi ;\
         if [ "$(HAVE_LAPACK_NATIVE)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND_NATIVE)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND_NATIVE"\
           "$(ENVIRONMENT_COMMAND_NATIVE)" ;\
          info_equal "F90n" "$(F90_NATIVE) -c $(FFLAGS_libs_NATIVE)" ;\
          info_equal "F90n-O0" "$(F90_NATIVE) -c $(FFLAGS0_libs_NATIVE)" ;\
          info_equal "ARn" "$(AR)" ;\
         else\
          info_line "Using external LAPACK library for utilities." ;\
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
        rm -rf $(OBJDIR) $(OBJDIR_NATIVE) $(LIBDIR)/$(LAPACK_AR)\
        $(LIBDIR)/$(LAPACK_AR_NATIVE)

# Object-file list
OBJ_OPT = $(OBJDIR)/dbdsdc.o $(OBJDIR)/dbdsqr.o $(OBJDIR)/dgebak.o \
 $(OBJDIR)/dgebal.o $(OBJDIR)/dgebd2.o $(OBJDIR)/dgebrd.o $(OBJDIR)/dgeev.o \
 $(OBJDIR)/dgehd2.o $(OBJDIR)/dgehrd.o $(OBJDIR)/dgelq2.o $(OBJDIR)/dgelqf.o \
 $(OBJDIR)/dgeqr2.o $(OBJDIR)/dgeqrf.o $(OBJDIR)/dgesdd.o $(OBJDIR)/dgetf2.o \
 $(OBJDIR)/dgetrf.o $(OBJDIR)/dgetrs.o $(OBJDIR)/dhseqr.o $(OBJDIR)/disnan.o \
 $(OBJDIR)/dlabad.o \
 $(OBJDIR)/dlabrd.o $(OBJDIR)/dlacpy.o $(OBJDIR)/dladiv.o $(OBJDIR)/dlae2.o \
 $(OBJDIR)/dlaed6.o $(OBJDIR)/dlaev2.o $(OBJDIR)/dlaexc.o $(OBJDIR)/dlahqr.o \
 $(OBJDIR)/dlahr2.o $(OBJDIR)/dlaisnan.o $(OBJDIR)/dlaln2.o $(OBJDIR)/dlamrg.o \
 $(OBJDIR)/dlange.o \
 $(OBJDIR)/dlanst.o $(OBJDIR)/dlanv2.o $(OBJDIR)/dlapy2.o $(OBJDIR)/dlapy3.o \
 $(OBJDIR)/dlaqr0.o $(OBJDIR)/dlaqr1.o $(OBJDIR)/dlaqr2.o $(OBJDIR)/dlaqr3.o \
 $(OBJDIR)/dlaqr4.o $(OBJDIR)/dlaqr5.o $(OBJDIR)/dlarfb.o $(OBJDIR)/dlarf.o \
 $(OBJDIR)/dlarfg.o $(OBJDIR)/dlarft.o $(OBJDIR)/dlarfx.o $(OBJDIR)/dlartg.o \
 $(OBJDIR)/dlas2.o $(OBJDIR)/dlascl.o $(OBJDIR)/dlasd0.o $(OBJDIR)/dlasd1.o \
 $(OBJDIR)/dlasd2.o $(OBJDIR)/dlasd3.o $(OBJDIR)/dlasd4.o $(OBJDIR)/dlasd5.o \
 $(OBJDIR)/dlasd6.o $(OBJDIR)/dlasd7.o $(OBJDIR)/dlasd8.o $(OBJDIR)/dlasda.o \
 $(OBJDIR)/dlasdq.o $(OBJDIR)/dlasdt.o $(OBJDIR)/dlaset.o $(OBJDIR)/dlasq1.o \
 $(OBJDIR)/dlasq2.o $(OBJDIR)/dlasq5.o $(OBJDIR)/dlasq6.o $(OBJDIR)/dlasr.o \
 $(OBJDIR)/dlasrt.o $(OBJDIR)/dlassq.o $(OBJDIR)/dlasv2.o $(OBJDIR)/dlaswp.o \
 $(OBJDIR)/dlasyf.o \
 $(OBJDIR)/dlasy2.o $(OBJDIR)/dlazq3.o $(OBJDIR)/dlazq4.o $(OBJDIR)/dorg2r.o \
 $(OBJDIR)/dorgbr.o $(OBJDIR)/dorghr.o $(OBJDIR)/dorgl2.o $(OBJDIR)/dorglq.o \
 $(OBJDIR)/dorgqr.o $(OBJDIR)/dorm2r.o $(OBJDIR)/dormbr.o $(OBJDIR)/dorml2.o \
 $(OBJDIR)/dormlq.o $(OBJDIR)/dormqr.o $(OBJDIR)/dsterf.o $(OBJDIR)/dsyr.o \
 $(OBJDIR)/dsytf2.o $(OBJDIR)/dsytrf.o $(OBJDIR)/dsytri.o $(OBJDIR)/dtrevc.o \
 $(OBJDIR)/dtrexc.o $(OBJDIR)/ieeeck.o $(OBJDIR)/ilaenv.o $(OBJDIR)/iparmq.o \
 $(OBJDIR)/xerbla.o $(OBJDIR)/zgebak.o $(OBJDIR)/zgebal.o $(OBJDIR)/zgebd2.o \
 $(OBJDIR)/zgebrd.o $(OBJDIR)/zgeev.o $(OBJDIR)/zgehd2.o $(OBJDIR)/zgehrd.o \
 $(OBJDIR)/zgelq2.o $(OBJDIR)/zgelqf.o $(OBJDIR)/zgeqr2.o $(OBJDIR)/zgeqrf.o \
 $(OBJDIR)/zgesdd.o $(OBJDIR)/zgetf2.o $(OBJDIR)/zgetrf.o $(OBJDIR)/zgetrs.o \
 $(OBJDIR)/zheev.o $(OBJDIR)/zhegs2.o $(OBJDIR)/zhegst.o $(OBJDIR)/zhegv.o \
 $(OBJDIR)/zhetd2.o $(OBJDIR)/zhetrd.o $(OBJDIR)/zhseqr.o $(OBJDIR)/zlabrd.o \
 $(OBJDIR)/zlacgv.o $(OBJDIR)/zlacp2.o $(OBJDIR)/zlacpy.o $(OBJDIR)/zlacrm.o \
 $(OBJDIR)/zladiv.o $(OBJDIR)/zlahqr.o $(OBJDIR)/zlahr2.o $(OBJDIR)/zlange.o \
 $(OBJDIR)/zlanhe.o $(OBJDIR)/zlaqr0.o $(OBJDIR)/zlaqr1.o $(OBJDIR)/zlaqr2.o \
 $(OBJDIR)/zlaqr3.o $(OBJDIR)/zlaqr4.o $(OBJDIR)/zlaqr5.o $(OBJDIR)/zlarcm.o \
 $(OBJDIR)/zlarfb.o $(OBJDIR)/zlarf.o $(OBJDIR)/zlarfg.o $(OBJDIR)/zlarft.o \
 $(OBJDIR)/zlartg.o $(OBJDIR)/zlascl.o $(OBJDIR)/zlaset.o $(OBJDIR)/zlasr.o \
 $(OBJDIR)/zlassq.o $(OBJDIR)/zlaswp.o $(OBJDIR)/zlatrd.o $(OBJDIR)/zlatrs.o \
 $(OBJDIR)/zpotf2.o $(OBJDIR)/zpotrf.o $(OBJDIR)/zrot.o $(OBJDIR)/zsteqr.o \
 $(OBJDIR)/ztrevc.o $(OBJDIR)/ztrexc.o $(OBJDIR)/zung2l.o $(OBJDIR)/zung2r.o \
 $(OBJDIR)/zungbr.o $(OBJDIR)/zunghr.o $(OBJDIR)/zungl2.o $(OBJDIR)/zunglq.o \
 $(OBJDIR)/zungql.o $(OBJDIR)/zungqr.o $(OBJDIR)/zungtr.o $(OBJDIR)/zunm2r.o \
 $(OBJDIR)/zunmbr.o $(OBJDIR)/zunml2.o $(OBJDIR)/zunmlq.o $(OBJDIR)/zunmqr.o
OBJ_NOOPT = $(OBJDIR)/dlamch.o
OBJ_OPT_NATIVE = $(OBJDIR_NATIVE)/dbdsdc.o $(OBJDIR_NATIVE)/dbdsqr.o \
 $(OBJDIR_NATIVE)/dgebak.o $(OBJDIR_NATIVE)/dgebal.o $(OBJDIR_NATIVE)/dgebd2.o \
 $(OBJDIR_NATIVE)/dgebrd.o $(OBJDIR_NATIVE)/dgeev.o $(OBJDIR_NATIVE)/dgehd2.o \
 $(OBJDIR_NATIVE)/dgehrd.o $(OBJDIR_NATIVE)/dgelq2.o $(OBJDIR_NATIVE)/dgelqf.o \
 $(OBJDIR_NATIVE)/dgeqr2.o $(OBJDIR_NATIVE)/dgeqrf.o $(OBJDIR_NATIVE)/dgesdd.o \
 $(OBJDIR_NATIVE)/dgetf2.o $(OBJDIR_NATIVE)/dgetrf.o $(OBJDIR_NATIVE)/dgetrs.o \
 $(OBJDIR_NATIVE)/dhseqr.o $(OBJDIR_NATIVE)/disnan.o $(OBJDIR_NATIVE)/dlabad.o \
 $(OBJDIR_NATIVE)/dlabrd.o \
 $(OBJDIR_NATIVE)/dlacpy.o $(OBJDIR_NATIVE)/dladiv.o $(OBJDIR_NATIVE)/dlae2.o \
 $(OBJDIR_NATIVE)/dlaed6.o $(OBJDIR_NATIVE)/dlaev2.o $(OBJDIR_NATIVE)/dlaexc.o \
 $(OBJDIR_NATIVE)/dlahqr.o $(OBJDIR_NATIVE)/dlahr2.o \
 $(OBJDIR_NATIVE)/dlaisnan.o $(OBJDIR_NATIVE)/dlaln2.o \
 $(OBJDIR_NATIVE)/dlamrg.o $(OBJDIR_NATIVE)/dlange.o $(OBJDIR_NATIVE)/dlanst.o \
 $(OBJDIR_NATIVE)/dlanv2.o $(OBJDIR_NATIVE)/dlapy2.o $(OBJDIR_NATIVE)/dlapy3.o \
 $(OBJDIR_NATIVE)/dlaqr0.o $(OBJDIR_NATIVE)/dlaqr1.o $(OBJDIR_NATIVE)/dlaqr2.o \
 $(OBJDIR_NATIVE)/dlaqr3.o $(OBJDIR_NATIVE)/dlaqr4.o $(OBJDIR_NATIVE)/dlaqr5.o \
 $(OBJDIR_NATIVE)/dlarfb.o $(OBJDIR_NATIVE)/dlarf.o $(OBJDIR_NATIVE)/dlarfg.o \
 $(OBJDIR_NATIVE)/dlarft.o $(OBJDIR_NATIVE)/dlarfx.o $(OBJDIR_NATIVE)/dlartg.o \
 $(OBJDIR_NATIVE)/dlas2.o $(OBJDIR_NATIVE)/dlascl.o $(OBJDIR_NATIVE)/dlasd0.o \
 $(OBJDIR_NATIVE)/dlasd1.o $(OBJDIR_NATIVE)/dlasd2.o $(OBJDIR_NATIVE)/dlasd3.o \
 $(OBJDIR_NATIVE)/dlasd4.o $(OBJDIR_NATIVE)/dlasd5.o $(OBJDIR_NATIVE)/dlasd6.o \
 $(OBJDIR_NATIVE)/dlasd7.o $(OBJDIR_NATIVE)/dlasd8.o $(OBJDIR_NATIVE)/dlasda.o \
 $(OBJDIR_NATIVE)/dlasdq.o $(OBJDIR_NATIVE)/dlasdt.o $(OBJDIR_NATIVE)/dlaset.o \
 $(OBJDIR_NATIVE)/dlasq1.o $(OBJDIR_NATIVE)/dlasq2.o $(OBJDIR_NATIVE)/dlasq5.o \
 $(OBJDIR_NATIVE)/dlasq6.o $(OBJDIR_NATIVE)/dlasr.o $(OBJDIR_NATIVE)/dlasrt.o \
 $(OBJDIR_NATIVE)/dlassq.o $(OBJDIR_NATIVE)/dlasv2.o $(OBJDIR_NATIVE)/dlaswp.o \
 $(OBJDIR_NATIVE)/dlasyf.o \
 $(OBJDIR_NATIVE)/dlasy2.o $(OBJDIR_NATIVE)/dlazq3.o $(OBJDIR_NATIVE)/dlazq4.o \
 $(OBJDIR_NATIVE)/dorg2r.o $(OBJDIR_NATIVE)/dorgbr.o $(OBJDIR_NATIVE)/dorghr.o \
 $(OBJDIR_NATIVE)/dorgl2.o $(OBJDIR_NATIVE)/dorglq.o $(OBJDIR_NATIVE)/dorgqr.o \
 $(OBJDIR_NATIVE)/dorm2r.o $(OBJDIR_NATIVE)/dormbr.o $(OBJDIR_NATIVE)/dorml2.o \
 $(OBJDIR_NATIVE)/dormlq.o $(OBJDIR_NATIVE)/dormqr.o $(OBJDIR_NATIVE)/dsterf.o \
 $(OBJDIR_NATIVE)/dsyr.o   $(OBJDIR_NATIVE)/dsytf2.o $(OBJDIR_NATIVE)/dsytrf.o \
 $(OBJDIR_NATIVE)/dsytri.o \
 $(OBJDIR_NATIVE)/dtrevc.o $(OBJDIR_NATIVE)/dtrexc.o $(OBJDIR_NATIVE)/ieeeck.o \
 $(OBJDIR_NATIVE)/ilaenv.o $(OBJDIR_NATIVE)/iparmq.o $(OBJDIR_NATIVE)/xerbla.o \
 $(OBJDIR_NATIVE)/zgebak.o $(OBJDIR_NATIVE)/zgebal.o $(OBJDIR_NATIVE)/zgebd2.o \
 $(OBJDIR_NATIVE)/zgebrd.o $(OBJDIR_NATIVE)/zgeev.o $(OBJDIR_NATIVE)/zgehd2.o \
 $(OBJDIR_NATIVE)/zgehrd.o $(OBJDIR_NATIVE)/zgelq2.o $(OBJDIR_NATIVE)/zgelqf.o \
 $(OBJDIR_NATIVE)/zgeqr2.o $(OBJDIR_NATIVE)/zgeqrf.o $(OBJDIR_NATIVE)/zgesdd.o \
 $(OBJDIR_NATIVE)/zgetf2.o $(OBJDIR_NATIVE)/zgetrf.o $(OBJDIR_NATIVE)/zgetrs.o \
 $(OBJDIR_NATIVE)/zheev.o $(OBJDIR_NATIVE)/zhegs2.o $(OBJDIR_NATIVE)/zhegst.o \
 $(OBJDIR_NATIVE)/zhegv.o $(OBJDIR_NATIVE)/zhetd2.o $(OBJDIR_NATIVE)/zhetrd.o \
 $(OBJDIR_NATIVE)/zhseqr.o $(OBJDIR_NATIVE)/zlabrd.o $(OBJDIR_NATIVE)/zlacgv.o \
 $(OBJDIR_NATIVE)/zlacp2.o $(OBJDIR_NATIVE)/zlacpy.o $(OBJDIR_NATIVE)/zlacrm.o \
 $(OBJDIR_NATIVE)/zladiv.o $(OBJDIR_NATIVE)/zlahqr.o $(OBJDIR_NATIVE)/zlahr2.o \
 $(OBJDIR_NATIVE)/zlange.o $(OBJDIR_NATIVE)/zlanhe.o $(OBJDIR_NATIVE)/zlaqr0.o \
 $(OBJDIR_NATIVE)/zlaqr1.o $(OBJDIR_NATIVE)/zlaqr2.o $(OBJDIR_NATIVE)/zlaqr3.o \
 $(OBJDIR_NATIVE)/zlaqr4.o $(OBJDIR_NATIVE)/zlaqr5.o $(OBJDIR_NATIVE)/zlarcm.o \
 $(OBJDIR_NATIVE)/zlarfb.o $(OBJDIR_NATIVE)/zlarf.o $(OBJDIR_NATIVE)/zlarfg.o \
 $(OBJDIR_NATIVE)/zlarft.o $(OBJDIR_NATIVE)/zlartg.o $(OBJDIR_NATIVE)/zlascl.o \
 $(OBJDIR_NATIVE)/zlaset.o $(OBJDIR_NATIVE)/zlasr.o $(OBJDIR_NATIVE)/zlassq.o \
 $(OBJDIR_NATIVE)/zlaswp.o $(OBJDIR_NATIVE)/zlatrd.o $(OBJDIR_NATIVE)/zlatrs.o \
 $(OBJDIR_NATIVE)/zpotf2.o $(OBJDIR_NATIVE)/zpotrf.o $(OBJDIR_NATIVE)/zrot.o \
 $(OBJDIR_NATIVE)/zsteqr.o $(OBJDIR_NATIVE)/ztrevc.o $(OBJDIR_NATIVE)/ztrexc.o \
 $(OBJDIR_NATIVE)/zung2l.o $(OBJDIR_NATIVE)/zung2r.o $(OBJDIR_NATIVE)/zungbr.o \
 $(OBJDIR_NATIVE)/zunghr.o $(OBJDIR_NATIVE)/zungl2.o $(OBJDIR_NATIVE)/zunglq.o \
 $(OBJDIR_NATIVE)/zungql.o $(OBJDIR_NATIVE)/zungqr.o $(OBJDIR_NATIVE)/zungtr.o \
 $(OBJDIR_NATIVE)/zunm2r.o $(OBJDIR_NATIVE)/zunmbr.o $(OBJDIR_NATIVE)/zunml2.o \
 $(OBJDIR_NATIVE)/zunmlq.o $(OBJDIR_NATIVE)/zunmqr.o
OBJ_NOOPT_NATIVE = $(OBJDIR_NATIVE)/dlamch.o
# End object-file list

# All objects
OBJ               = $(OBJ_OPT) $(OBJ_NOOPT)
OBJ_NATIVE_       =
OBJ_NATIVE_native = $(OBJ_OPT_NATIVE) $(OBJ_NOOPT_NATIVE)
OBJ_NATIVE        = $(OBJ_NATIVE_$(UTILS_MODE))

# How to compile things
$(LAPACK): $(OBJ)
	@bash -c '$(MKOBJDIR) ; a=$(@F) ; b=$${a%.a} ;\
        $(BASH_PRETTY) ; compile_info AR lib LAPACK ;\
        rm -f $(LAPACK) ; $(AR) $(LAPACK) $(OBJ) $(REDIR)'

$(OBJ_OPT):
	@bash -c '$(MKOBJDIR) ; o=$(@F) ; b=$${o%.o} ; f=$$b.f ;\
        $(BASH_PRETTY) ; compile_info F90 lib LAPACK $$f ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(F90) $(FFLAGS_libs) -c $(LAPACK_SRCDIR)/$$f -o $@ $(REDIR)'

$(OBJ_NOOPT):
	@bash -c '$(MKOBJDIR) ; o=$(@F) ; b=$${o%.o} ; f=$$b.f ;\
        $(BASH_PRETTY) ; compile_info F90-O0 lib LAPACK $$f ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(F90) $(FFLAGS0_libs) -c $(LAPACK_SRCDIR)/$$f -o $@ $(REDIR)'

$(LAPACK_NATIVE): $(OBJ_NATIVE)
	@bash -c '$(BASH_PRETTY) ;\
        if [ "$(UTILS_MODE)" = native ] ; then \
         $(MKOBJDIR_NATIVE) ; a=$(@F) ; b=$${a%.a} ;\
         compile_info ARn lib LAPACK_NATIVE ;\
         rm -f $(LAPACK_NATIVE) ;\
         $(AR) $(LAPACK_NATIVE) $(OBJ_NATIVE) $(REDIR) ;\
        else \
         compile_info LN lib LAPACK_NATIVE ;\
         rm -f $(LAPACK_NATIVE) ;\
         ln -s $(LIBDIR)/$(LAPACK_AR) $(LAPACK_NATIVE) ;\
        fi'

$(OBJ_OPT_NATIVE):
	@bash -c '$(MKOBJDIR_NATIVE) ; o=$(@F) ; b=$${o%.o} ; f=$$b.f ;\
        $(BASH_PRETTY) ; compile_info F90n lib LAPACK_NATIVE $$f ;\
        $(ENVIRONMENT_COMMAND_NATIVE) ; cd $(OBJDIR_NATIVE) &&\
        $(F90_NATIVE) $(FFLAGS_libs_NATIVE) -c $(LAPACK_SRCDIR)/$$f\
        -o $@ $(REDIR)'

$(OBJ_NOOPT_NATIVE):
	@bash -c '$(MKOBJDIR_NATIVE) ; o=$(@F) ; b=$${o%.o} ; f=$$b.f ;\
        $(BASH_PRETTY) ; compile_info F90n-O0 lib LAPACK_NATIVE $$f ;\
        $(ENVIRONMENT_COMMAND_NATIVE) ; cd $(OBJDIR_NATIVE) &&\
        $(F90_NATIVE) $(FFLAGS0_libs_NATIVE) -c $(LAPACK_SRCDIR)/$$f\
        -o $@ $(REDIR)'

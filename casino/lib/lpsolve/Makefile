#--------------------------------------------------#
# Makefile to create a static library for LPSOLVE. #
#                                                  #
# PLR 12.2008                                      #
#--------------------------------------------------#
# Basic stuff
SHELL             = /bin/sh
LPSOLVE_SRCDIR    = $(PWD)
INCBASE           = $(LPSOLVE_SRCDIR)/../../arch/data
OBJBASE           = $(LPSOLVE_SRCDIR)/../zlib
TOOLDIR           = $(LPSOLVE_SRCDIR)/../../build-tools
LPSOLVE_AR        = liblpsolve.a
LPSOLVE_AR_NATIVE = liblpsolve_native.a
INCLUDE_COLAMD    = -I$(LPSOLVE_SRCDIR)/../colamd
SILENT            = 0
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

# System-independent variables
LIBDIR             = $(OBJBASE)/$(CASINO_ARCH)
OBJDIR             = $(LIBDIR)/lpsolve
MKOBJDIR           = [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
LPSOLVE_yes        =
LPSOLVE_no         = $(LIBDIR)/$(LPSOLVE_AR)
OBJDIR_NATIVE      = $(LIBDIR)/lpsolve_native
MKOBJDIR_NATIVE    = [ -d $(OBJDIR_NATIVE) ] || mkdir -p $(OBJDIR_NATIVE)
LPSOLVE_NATIVE_yes =
LPSOLVE_NATIVE_no  = $(LIBDIR)/$(LPSOLVE_AR_NATIVE)

# Defaults for system-dependent variables
AR                         = ar crs
REDIR_0                    =
REDIR_1                    = >& /dev/null
CFLAGS_libs                = $(CFLAGS_opt) $(CFLAGS_all)
HAVE_LPSOLVE               = no
ENVIRONMENT_COMMAND        = :
CFLAGS_libs_NATIVE         = $(CFLAGS_opt_NATIVE) $(CFLAGS_all_NATIVE)
HAVE_LPSOLVE_NATIVE        = no
ENVIRONMENT_COMMAND_NATIVE = :

# System-dependent include files
include $(INCBASE)/$(CASINO_ARCH).arch

# System-dependent variables from data in include files
REDIR                      = $(REDIR_$(SILENT))
LPSOLVE                    = $(LPSOLVE_$(HAVE_LPSOLVE))
HAVE_LPSOLVE_NATIVE_native = $(HAVE_LPSOLVE_NATIVE)
HAVE_LPSOLVE_NATIVE_       = $(HAVE_LPSOLVE)
LPSOLVE_NATIVE             = $(LPSOLVE_NATIVE_$(HAVE_LPSOLVE_NATIVE_$(UTILS_MODE)))

# Main target
.PHONY: default
default: info $(LPSOLVE) $(LPSOLVE_NATIVE)

# Info
.PHONY: info
info:
	@bash -c '$(BASH_PRETTY) ;\
        info_title "CASINO LPSOLVE MAKE information" ;\
        if [ "$(UTILS_MODE)" != native ] ; then\
         if [ "$(HAVE_LPSOLVE)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
          info_equal "CC" "$(CC) -c $(CFLAGS_libs)" ;\
          info_equal "AR" "$(AR)" ;\
         else\
          info_line "Using external LPSOLVE library." ;\
         fi ;\
        else\
         if [ "$(HAVE_LPSOLVE)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
          info_equal "CC" "$(CC) -c $(CFLAGS_libs)" ;\
          info_equal "AR" "$(AR)" ;\
         else\
          info_line "Using external LPSOLVE library for main binary." ;\
         fi ;\
         if [ "$(HAVE_LPSOLVE_NATIVE)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND_NATIVE)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND_NATIVE"\
           "$(ENVIRONMENT_COMMAND_NATIVE)" ;\
          info_equal "CCn" "$(CC_NATIVE) -c $(CFLAGS_libs_NATIVE)" ;\
          info_equal "ARn" "$(AR)" ;\
         else\
          info_line "Using external LPSOLVE library for utilities." ;\
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
        rm -rf $(OBJDIR) $(OBJDIR_NATIVE) $(LIBDIR)/$(LPSOLVE_AR)\
        $(LIBDIR)/$(LPSOLVE_AR_NATIVE)

# Object-file list
OBJ_OPT = \
$(OBJDIR)/lp_MDO.o \
$(OBJDIR)/commonlib.o \
$(OBJDIR)/mmio.o \
$(OBJDIR)/myblas.o \
$(OBJDIR)/ini.o \
$(OBJDIR)/fortify.o \
$(OBJDIR)/lp_rlp.o \
$(OBJDIR)/lp_crash.o \
$(OBJDIR)/lp_LUSOL.o \
$(OBJDIR)/lusol.o \
$(OBJDIR)/lp_Hash.o \
$(OBJDIR)/lp_lib.o \
$(OBJDIR)/lp_wlp.o \
$(OBJDIR)/lp_matrix.o \
$(OBJDIR)/lp_mipbb.o \
$(OBJDIR)/lp_MPS.o \
$(OBJDIR)/lp_params.o \
$(OBJDIR)/lp_presolve.o \
$(OBJDIR)/lp_price.o \
$(OBJDIR)/lp_pricePSE.o \
$(OBJDIR)/lp_report.o \
$(OBJDIR)/lp_scale.o \
$(OBJDIR)/lp_simplex.o \
$(OBJDIR)/lp_SOS.o \
$(OBJDIR)/lp_utils.o \
$(OBJDIR)/yacc_read.o
OBJ_OPT_NATIVE = \
$(OBJDIR_NATIVE)/lp_MDO.o \
$(OBJDIR_NATIVE)/commonlib.o \
$(OBJDIR_NATIVE)/mmio.o \
$(OBJDIR_NATIVE)/myblas.o \
$(OBJDIR_NATIVE)/ini.o \
$(OBJDIR_NATIVE)/fortify.o \
$(OBJDIR_NATIVE)/lp_rlp.o \
$(OBJDIR_NATIVE)/lp_crash.o \
$(OBJDIR_NATIVE)/lp_LUSOL.o \
$(OBJDIR_NATIVE)/lusol.o \
$(OBJDIR_NATIVE)/lp_Hash.o \
$(OBJDIR_NATIVE)/lp_lib.o \
$(OBJDIR_NATIVE)/lp_wlp.o \
$(OBJDIR_NATIVE)/lp_matrix.o \
$(OBJDIR_NATIVE)/lp_mipbb.o \
$(OBJDIR_NATIVE)/lp_MPS.o \
$(OBJDIR_NATIVE)/lp_params.o \
$(OBJDIR_NATIVE)/lp_presolve.o \
$(OBJDIR_NATIVE)/lp_price.o \
$(OBJDIR_NATIVE)/lp_pricePSE.o \
$(OBJDIR_NATIVE)/lp_report.o \
$(OBJDIR_NATIVE)/lp_scale.o \
$(OBJDIR_NATIVE)/lp_simplex.o \
$(OBJDIR_NATIVE)/lp_SOS.o \
$(OBJDIR_NATIVE)/lp_utils.o \
$(OBJDIR_NATIVE)/yacc_read.o
# End object-file list

# Handle native compilation
OBJ               = $(OBJ_OPT)
OBJ_NATIVE_       =
OBJ_NATIVE_native = $(OBJ_OPT_NATIVE)
OBJ_NATIVE        = $(OBJ_NATIVE_$(UTILS_MODE))
BUNCH_OF_DEFINES  = -DYY_NEVER_INTERACTIVE -DPARSER_LP\
 -DINVERSE_ACTIVE=INVERSE_LUSOL -DRoleIsExternalInvEngine

# How to compile things
$(LPSOLVE): $(OBJ)
	@bash -c '$(MKOBJDIR) ; a=$(@F) ; b=$${a%.a} ;\
        $(BASH_PRETTY) ; compile_info AR lib LPSOLVE ;\
        rm -f $(LPSOLVE) ; $(AR) $(LPSOLVE) $(OBJ) $(REDIR)'

$(OBJ_OPT):
	@bash -c '$(MKOBJDIR) ; o=$(@F) ; b=$${o%.o} ; f=$$b.c ;\
        $(BASH_PRETTY) ; compile_info CC lib LPSOLVE $$f ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(CC) $(CFLAGS_libs) -c -I$(LPSOLVE_SRCDIR) $(BUNCH_OF_DEFINES)\
        $(LPSOLVE_SRCDIR)/$$f $(INCLUDE_COLAMD) $(REDIR)'

$(LPSOLVE_NATIVE): $(OBJ_NATIVE)
	@bash -c '$(BASH_PRETTY) ;\
        if [ "$(UTILS_MODE)" = native ] ; then \
         $(MKOBJDIR_NATIVE) ; a=$(@F) ; b=$${a%.a} ;\
         compile_info ARn lib LPSOLVE_NATIVE ;\
         rm -f $(LPSOLVE_NATIVE) ;\
         $(AR) $(LPSOLVE_NATIVE) $(OBJ_NATIVE) $(REDIR) ;\
        else \
         compile_info LN lib LPSOLVE_NATIVE ;\
         rm -f $(LPSOLVE_NATIVE) ;\
         ln -s $(LIBDIR)/$(LPSOLVE_AR) $(LPSOLVE_NATIVE) ;\
        fi'

$(OBJ_OPT_NATIVE):
	@bash -c '$(MKOBJDIR_NATIVE) ; o=$(@F) ; b=$${o%.o} ; f=$$b.c ;\
        $(BASH_PRETTY) ; compile_info CCn lib LPSOLVE_NATIVE $$f ;\
        $(ENVIRONMENT_COMMAND_NATIVE) ; cd $(OBJDIR_NATIVE) &&\
        $(CC_NATIVE) $(CFLAGS_libs_NATIVE) -c -I$(LPSOLVE_SRCDIR)\
        $(BUNCH_OF_DEFINES) $(LPSOLVE_SRCDIR)/$$f $(INCLUDE_COLAMD) $(REDIR)'

#-------------------------------------------------#
# Makefile to create a static library for COLAMD. #
#                                                 #
# PLR 12.2013                                     #
#-------------------------------------------------#
# Basic stuff
SHELL            = /bin/sh
COLAMD_SRCDIR    = $(PWD)
INCBASE          = $(COLAMD_SRCDIR)/../../arch/data
OBJBASE          = $(COLAMD_SRCDIR)/../zlib
TOOLDIR          = $(COLAMD_SRCDIR)/../../build-tools
COLAMD_AR        = libcolamd.a
COLAMD_AR_NATIVE = libcolamd_native.a
SILENT           = 0
.SILENT:
include $(TOOLDIR)/bash-pretty.mk

# System-independent variables
LIBDIR            = $(OBJBASE)/$(CASINO_ARCH)
OBJDIR            = $(LIBDIR)/colamd
MKOBJDIR          = [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR)
COLAMD_yes        =
COLAMD_no         = $(LIBDIR)/$(COLAMD_AR)
OBJDIR_NATIVE     = $(LIBDIR)/colamd_native
MKOBJDIR_NATIVE   = [ -d $(OBJDIR_NATIVE) ] || mkdir -p $(OBJDIR_NATIVE)
COLAMD_NATIVE_yes =
COLAMD_NATIVE_no  = $(LIBDIR)/$(COLAMD_AR_NATIVE)

# Defaults for system-dependent variables
AR                         = ar crs
REDIR_0                    =
REDIR_1                    = >& /dev/null
CFLAGS_libs                = $(CFLAGS_opt) $(CFLAGS_all)
HAVE_COLAMD                = no
ENVIRONMENT_COMMAND        = :
CFLAGS_libs_NATIVE         = $(CFLAGS_opt_NATIVE) $(CFLAGS_all_NATIVE)
HAVE_COLAMD_NATIVE         = no
ENVIRONMENT_COMMAND_NATIVE = :

# System-dependent include files
include $(INCBASE)/$(CASINO_ARCH).arch

# System-dependent variables from data in include files
REDIR                     = $(REDIR_$(SILENT))
COLAMD                    = $(COLAMD_$(HAVE_COLAMD))
HAVE_COLAMD_NATIVE_native = $(HAVE_COLAMD_NATIVE)
HAVE_COLAMD_NATIVE_       = $(HAVE_COLAMD)
COLAMD_NATIVE             = $(COLAMD_NATIVE_$(HAVE_COLAMD_NATIVE_$(UTILS_MODE)))

# Main target
.PHONY: default
default: info $(COLAMD) $(COLAMD_NATIVE)

# Info
.PHONY: info
info:
	@bash -c '$(BASH_PRETTY) ;\
        info_title "CASINO COLAMD MAKE information" ;\
        if [ "$(UTILS_MODE)" != native ] ; then\
         if [ "$(HAVE_COLAMD)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
          info_equal "CC" "$(CC) -c $(CFLAGS_libs)" ;\
          info_equal "AR" "$(AR)" ;\
         else\
          info_line "Using external COLAMD library." ;\
         fi ;\
        else\
         if [ "$(HAVE_COLAMD)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND" "$(ENVIRONMENT_COMMAND)" ;\
          info_equal "CC" "$(CC) -c $(CFLAGS_libs)" ;\
          info_equal "AR" "$(AR)" ;\
         else\
          info_line "Using external COLAMD library for main binary." ;\
         fi ;\
         if [ "$(HAVE_COLAMD_NATIVE)" = no ] ; then\
          [ "$(ENVIRONMENT_COMMAND_NATIVE)" != : ]\
           && info_equal "ENVIRONMENT_COMMAND_NATIVE"\
           "$(ENVIRONMENT_COMMAND_NATIVE)" ;\
          info_equal "CCn" "$(CC_NATIVE) -c $(CFLAGS_libs_NATIVE)" ;\
          info_equal "ARn" "$(AR)" ;\
         else\
          info_line "Using external COLAMD library for utilities." ;\
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
        rm -rf $(OBJDIR) $(OBJDIR_NATIVE) $(LIBDIR)/$(COLAMD_AR)\
        $(LIBDIR)/$(COLAMD_AR_NATIVE)

# Object-file list
OBJ_OPT = \
$(OBJDIR)/colamd.o
OBJ_OPT_NATIVE = \
$(OBJDIR_NATIVE)/colamd.o
# End object-file list

# Handle native compilation
OBJ               = $(OBJ_OPT)
OBJ_NATIVE_       =
OBJ_NATIVE_native = $(OBJ_OPT_NATIVE)
OBJ_NATIVE        = $(OBJ_NATIVE_$(UTILS_MODE))

# How to compile things
$(COLAMD): $(OBJ)
	@bash -c '$(MKOBJDIR) ; a=$(@F) ; b=$${a%.a} ;\
        $(BASH_PRETTY) ; compile_info AR lib COLAMD ;\
        rm -f $(COLAMD) ; $(AR) $(COLAMD) $(OBJ) $(REDIR)'

$(OBJ_OPT):
	@bash -c '$(MKOBJDIR) ; o=$(@F) ; b=$${o%.o} ; f=$$b.c ;\
        $(BASH_PRETTY) ; compile_info CC lib COLAMD $$f ;\
        $(ENVIRONMENT_COMMAND) ; cd $(OBJDIR) &&\
        $(CC) $(CFLAGS_libs) -c -I$(COLAMD_SRCDIR) $(COLAMD_SRCDIR)/$$f\
        $(REDIR)'

$(COLAMD_NATIVE): $(OBJ_NATIVE)
	@bash -c '$(BASH_PRETTY) ;\
        if [ "$(UTILS_MODE)" = native ] ; then \
         $(MKOBJDIR_NATIVE) ; a=$(@F) ; b=$${a%.a} ;\
         compile_info ARn lib COLAMD_NATIVE ;\
         rm -f $(COLAMD_NATIVE) ;\
         $(AR) $(COLAMD_NATIVE) $(OBJ_NATIVE) $(REDIR) ;\
        else \
         compile_info LN lib COLAMD_NATIVE ;\
         rm -f $(COLAMD_NATIVE) ;\
         ln -s $(LIBDIR)/$(COLAMD_AR) $(COLAMD_NATIVE) ;\
        fi'

$(OBJ_OPT_NATIVE):
	@bash -c '$(MKOBJDIR_NATIVE) ; o=$(@F) ; b=$${o%.o} ; f=$$b.c ;\
        $(BASH_PRETTY) ; compile_info CCn lib COLAMD_NATIVE $$f ;\
        $(ENVIRONMENT_COMMAND_NATIVE) ; cd $(OBJDIR_NATIVE) &&\
        $(CC_NATIVE) $(CFLAGS_libs_NATIVE) -c -I$(COLAMD_SRCDIR)\
        $(COLAMD_SRCDIR)/$$f $(REDIR)'

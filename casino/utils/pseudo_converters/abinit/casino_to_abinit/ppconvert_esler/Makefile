#-------------------------#
# Makefile for ppconvert. #
#                         #
# PLR 11.2008             #
#-------------------------#
# Basic stuff
UTILEXE = ppconvert
SRCDIR  = $(PWD)
UTILDIR = $(SRCDIR)/../../../..
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
LDCXX                      = $(CXX)
LDCXXFLAGS_opt             = $(CXXFLAGS_opt)
LDCXXFLAGS_all             = $(CXXFLAGS_all)
ENVIRONMENT_COMMAND        = :
ENVIRONMENT_COMMAND_NATIVE = :

# System dependent include files
include $(INCBASE)/$(CASINO_ARCH).arch

# Use the same set-up as for src, unless overrides defined
UNATIVE                   = $(UNATIVE_$(UTILS_MODE))
ENVIRONMENT_COMMAND_UTILS = $(ENVIRONMENT_COMMAND$(UNATIVE))
INCLUDE_DIR_UTILS         = $(INCLUDE_DIR$(UNATIVE))
LIB_PATH_UTILS            = $(LIB_PATH)

# System-dependent variables
CXXFLAGS   = $(CXXFLAGS_opt) $(CXXFLAGS_all) $(INCLUDE_DIR_UTILS)
LDCXXFLAGS = $(LDCXXFLAGS_opt) $(LDCXXFLAGS_all) $(LIB_PATH_UTILS)
LDCXXLIBS  = $(LDCXXLIBS_opt) $(LDCXXLIBS_all)

# Main target
.PHONY: default
default: $(BINDIR)/$(UTILEXE)

# 'Clean' target
.PHONY: clean
clean:
	@rm -rf $(OBJDIR) $(BINDIR)/$(UTILEXE)

# Object-file list
OBJ = \
$(OBJDIR)/CubicSpline.o \
$(OBJDIR)/NLPPClass.o \
$(OBJDIR)/ParseCommand.o \
$(OBJDIR)/ParserClass.o \
$(OBJDIR)/XMLWriterClass2.o
# End object-file list

# Dependencies
$(OBJDIR)/CubicSpline.o: $(SRCDIR)/CubicSpline.h $(SRCDIR)/CubicSpline.cc
CubicSpline.o_name = CubicSpline

$(OBJDIR)/NLPPClass.o: $(SRCDIR)/CubicSpline.h $(SRCDIR)/ParseCommand.h \
 $(SRCDIR)/ParserClass.h $(SRCDIR)/XMLWriterClass2.h $(SRCDIR)/NLPPClass.h \
 $(SRCDIR)/NLPPClass.cc
NLPPClass.o_name = NLPPClass

$(OBJDIR)/ParseCommand.o: $(SRCDIR)/ParseCommand.cc
ParseCommand.o_name = ParseCommand

$(OBJDIR)/ParserClass.o: $(SRCDIR)/ParserClass.cc
ParserClass.o_name = ParserClass

$(OBJDIR)/XMLWriterClass2.o: $(SRCDIR)/XMLWriterClass2.cc
XMLWriterClass2.o_name = XMLWriterClass2
# End dependencies

# How to make the objects
$(OBJ):
	@[ -d $(OBJDIR) ] || mkdir -p $(OBJDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info CXX utils $(UTILEXE)\
        $($(@F)_name).cc' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ; cd $(OBJDIR)\
        && $(CXX) -c $(CXXFLAGS) $(SRCDIR)/$($(@F)_name).cc
$(BINDIR)/$(UTILEXE) : $(OBJ)
	@[ -d $(BINDIR) ] || mkdir -p $(BINDIR) ;\
        bash -c '$(BASH_PRETTY) ; compile_info LDCXX utils $(UTILEXE)' ;\
        $(ENVIRONMENT_COMMAND_UTILS) ;\
        $(LDCXX) $(LDCXXFLAGS) -o $(@) $(OBJ) $(LDCXXLIBS)

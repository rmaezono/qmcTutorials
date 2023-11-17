#----------------------------------#
# Makefile for external libraries. #
# PLR 07.2009                      #
#----------------------------------#
# Basic stuff
LIBBASE    = $(PWD)
BLASDIR    = $(LIBBASE)/BLAS
LAPACKDIR  = $(LIBBASE)/LAPACK
LPSOLVEDIR = $(LIBBASE)/lpsolve
COLAMDDIR  = $(LIBBASE)/colamd
INCBASE    = $(LIBBASE)/../arch/data
.SILENT:

.PHONY: start
start: default

# System-dependent include file
include $(INCBASE)/$(CASINO_ARCH).arch

# Info targets
.PHONY: info blas_info lapack_info lpsolve_info colamd_info
info: blas_info lapack_info lpsolve_info colamd_info
blas_info:
	@cd $(BLASDIR) ; $(MAKE) info PWD="$(BLASDIR)"
lapack_info:
	@cd $(LAPACKDIR) ; $(MAKE) info PWD="$(LAPACKDIR)"
lpsolve_info:
	@cd $(LPSOLVEDIR) ;\
        [ ! -z "$(CC)" ] && $(MAKE) info PWD="$(LPSOLVEDIR)"
colamd_info:
	@cd $(COLAMDDIR) ;\
        [ ! -z "$(CC)" ] && $(MAKE) info PWD="$(COLAMDDIR)"

# Clean targets
.PHONY: clean blas_clean lapack_clean lpsolve_clean colamd_clean
clean: blas_clean lapack_clean lpsolve_clean colamd_clean
blas_clean:
	@cd $(BLASDIR) ; $(MAKE) clean PWD="$(BLASDIR)"
lapack_clean:
	@cd $(LAPACKDIR) ; $(MAKE) clean PWD="$(LAPACKDIR)"
lpsolve_clean:
	@cd $(LPSOLVEDIR) ;\
        [ ! -z "$(CC)" ] && $(MAKE) clean PWD="$(LPSOLVEDIR)"
colamd_clean:
	@cd $(COLAMDDIR) ;\
        [ ! -z "$(CC)" ] && $(MAKE) clean PWD="$(COLAMDDIR)"

# Very-clean targets
.PHONY: vclean blas_vclean lapack_vclean lpsolve_vclean colamd_vclean
vclean: blas_vclean lapack_vclean lpsolve_vclean colamd_vclean
blas_vclean:
	@cd $(BLASDIR) ; $(MAKE) vclean PWD="$(BLASDIR)"
lapack_vclean:
	@cd $(LAPACKDIR) ; $(MAKE) vclean PWD="$(LAPACKDIR)"
lpsolve_vclean:
	@cd $(LPSOLVEDIR) ;\
        [ ! -z "$(CC)" ] && $(MAKE) vclean PWD="$(LPSOLVEDIR)"
colamd_vclean:
	@cd $(COLAMDDIR) ;\
        [ ! -z "$(CC)" ] && $(MAKE) vclean PWD="$(COLAMDDIR)"

# Build targets
.PHONY: blas lapack lpsolve colamd
default: blas lapack lpsolve colamd
blas:
	@cd $(BLASDIR) ; $(MAKE) PWD="$(BLASDIR)"
lapack:
	@cd $(LAPACKDIR) ; $(MAKE) PWD="$(LAPACKDIR)"
lpsolve:
	@cd $(LPSOLVEDIR) ;\
        [ ! -z "$(CC)" ] && $(MAKE) PWD="$(LPSOLVEDIR)"
colamd:
	@cd $(COLAMDDIR) ;\
        [ ! -z "$(CC)" ] && $(MAKE) PWD="$(COLAMDDIR)"

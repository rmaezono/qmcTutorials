#-------------------------------------------------------------------------#
# Global Makefile for CASINO.                                             #
# PLR 11.2008                                                             #
#-------------------------------------------------------------------------#
# Basic stuff
SHELL          = /bin/sh
VERSION        = opt
CASINODIR      = $(PWD)
TOOLDIR        = $(CASINODIR)/build-tools
VERSION_FILE   = $(CASINODIR)/VERSION
DIARY_FILE     = $(CASINODIR)/DIARY
GIT_REPOSITORY = $(CASINODIR)/.git
VERSION_SCRIPT = $(TOOLDIR)/update-version
include $(TOOLDIR)/bash-pretty.mk

# Prevent GNU make from cluttering the output with "make[3]: Entering directory
# `/some/random/huge/pathname/which/no/one/cares/about/but/goes/on/forever/and
# /then/a/bit/more/'".
.SILENT:

.PHONY: start
start: default

# Check CASINO_ARCH target (and for presence of git repository if required).
.PHONY: checkarch
checkarch: global_info $(VERSION_FILE)
	@if [ -z "$(CASINO_ARCH)" ] ; then\
         echo "CASINO_ARCH not set." ;\
         exit 1 ;\
        fi

# Version file
.PHONY: $(VERSION_FILE)
$(VERSION_FILE): global_info
	@cd $(CASINODIR) && $(VERSION_SCRIPT)

# Info targets
.PHONY: info global_info lib_info src_info utils_info manual_info
info: global_info lib_info src_info utils_info manual_info
global_info:
	@bash -c '$(BASH_PRETTY) ;\
        info_title "CASINO global MAKE information" ;\
        info_equal "CASINO_ARCH" "$(CASINO_ARCH)" ;\
        info_end'
lib_info: checkarch
	@cd lib ; $(MAKE) info PWD=$(CASINODIR)/lib
src_info: checkarch $(VERSION_FILE)
	@cd src ; $(MAKE) info VERSION=$(VERSION) PWD=$(CASINODIR)/src
utils_info: checkarch
	@cd utils ; $(MAKE) info PWD=$(CASINODIR)/utils
manual_info: global_info
	@[ -d manual/src ] && cd manual/src\
        && $(MAKE) info PWD=$(CASINODIR)/manual/src || :

# Clean targets
.PHONY: clean lib_clean src_clean utils_clean manual_clean
clean: lib_clean src_clean utils_clean manual_clean
lib_clean: checkarch
	@cd lib ; $(MAKE) clean PWD=$(CASINODIR)/lib
src_clean: checkarch $(VERSION_FILE)
	@cd src ; $(MAKE) all/all.clean PWD=$(CASINODIR)/src
utils_clean: checkarch
	@cd utils ; $(MAKE) clean PWD=$(CASINODIR)/utils
manual_clean: global_info
	@[ -d manual/src ] && cd manual/src\
        && $(MAKE) clean PWD=$(CASINODIR)/manual/src || :

# Very-clean targets
.PHONY: vclean lib_vclean src_vclean utils_vclean manual_vclean
vclean: lib_vclean src_vclean utils_vclean manual_vclean
lib_vclean: checkarch
	@cd lib ; $(MAKE) vclean PWD=$(CASINODIR)/lib
src_vclean: checkarch $(VERSION_FILE)
	@cd src ; $(MAKE) all/all.vclean PWD=$(CASINODIR)/src
utils_vclean: checkarch
	@cd utils ; $(MAKE) vclean PWD=$(CASINODIR)/utils
manual_vclean: global_info
	@[ -d manual/src ] && cd manual/src\
        && $(MAKE) vclean PWD=$(CASINODIR)/manual/src || :

# Build targets
.PHONY: lib src utils manual
default: lib src utils manual
lib: checkarch
	@cd lib ; $(MAKE) PWD=$(CASINODIR)/lib
src: checkarch lib $(VERSION_FILE)
	@cd src ; $(MAKE) VERSION=$(VERSION) PWD=$(CASINODIR)/src
utils: checkarch lib
	@cd utils ; $(MAKE) PWD=$(CASINODIR)/utils
manual: global_info
	@[ -d manual/src ] && cd manual/src\
        && $(MAKE) PWD=$(CASINODIR)/manual/src || :

# Extra targets for src features and versions.
.PHONY: NoFeatures opt NoFeatures/opt dev NoFeatures/dev debug\
 NoFeatures/debug prof NoFeatures/prof nofeatures nofeatures/opt\
 nofeatures/dev nofeatures/debug nofeatures/prof
NoFeatures nofeatures NoFeatures/opt nofeatures/opt opt:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=opt
NoFeatures/dev nofeatures/dev dev:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=dev
NoFeatures/debug nofeatures/debug debug:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=debug
NoFeatures/prof nofeatures/prof prof:
	@$(MAKE) USE_OPENMP=no USE_SHM=no VERSION=prof
.PHONY: Openmp Openmp/opt Openmp/dev Openmp/debug Openmp/prof\
 openmp openmp/opt openmp/dev openmp/debug openmp/prof
Openmp openmp Openmp/opt openmp/opt:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=opt
Openmp/dev openmp/dev:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=dev
Openmp/debug openmp/debug:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=debug
Openmp/prof openmp/prof:
	@$(MAKE) USE_OPENMP=yes USE_SHM=no VERSION=prof
.PHONY: Shm Shm/opt Shm/dev Shm/debug Shm/prof shm/opt shm/dev shm/debug\
 shm/prof
Shm Shm/opt shm shm/opt:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=opt
Shm/dev shm/dev:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=dev
Shm/debug shm/debug:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=debug
Shm/prof shm/prof:
	@$(MAKE) USE_OPENMP=no USE_SHM=yes VERSION=prof
.PHONY: OpenmpShm OpenmpShm/opt OpenmpShm/dev OpenmpShm/debug OpenmpShm/prof\
 openmpshm openmpshm/opt openmpshm/dev openmpshm/debug openmpshm/prof
OpenmpShm OpenmpShm/opt openmpshm openmpshm/opt:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=opt
OpenmpShm/dev openmpshm/dev:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=dev
OpenmpShm/debug openmpshm/debug:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=debug
OpenmpShm/prof openmpshm/prof:
	@$(MAKE) USE_OPENMP=yes USE_SHM=yes VERSION=prof

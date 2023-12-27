# makefile: quodas
# Kyle A. Damm (27-01-2020)

DMD := ldmd2
VERSION ?= fast
ifeq ($(VERSION), debug)
	DFLAGS := -w -g -debug
        OF := -of
endif
ifeq ($(VERSION), fast)
	DFLAGS := -w -g -O -release -inline -boundscheck=off
        OF := -of
endif
WITH_COMPLEX_NUMBERS ?= 0
ifeq ($(WITH_COMPLEX_NUMBERS), 1)
	CMPLXFLAG := -version=complex_numbers
else
	CMPLXFLAG := -version=real_numbers
endif

NTYPES := ../../ntypes
NM := ../../nm
UTIL_DIR := ../../util
include $(UTIL_DIR)/util_files.mk

FILES := main.d \
	$(NTYPES)/complex.d \
	$(NM)/number.d \
	config.d \
	block.d \
	finite_volume.d \
	flux.d \
	flow_solver.d \
	optimizer.d \
	linalg.d \
	derivative.d \
	$(UTIL_DIR)/msg_service.d

quodas: $(FILES)
	$(DMD) $(DFLAGS) $(OF)$@ $(FILES) $(CMPLXFLAG)

clean:
	- rm *.o
	- rm quodas

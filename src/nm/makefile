# makefile for the Numerical Methods package
#
# We can specify the LDC2 compiler as DMD=ldmd2 on the command-line
# when invoking this makefile.  Can also ask for gdc.
DMD ?= ldmd2

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

NTYPES_DIR := ../ntypes
NTYPES_FILES := $(NTYPES_DIR)/complex.d

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

include nm_files.mk

DEMO_PROGRAMS := ridder_demo \
	linesearch_demo \
	nelmin_demo nelmin_complex_demo \
	rungekutta_demo \
	newtoncotes_demo \
	gaussquad_demo \
	bbla_demo

TEST_PROGRAMS := bbla_test bbla_complex_test \
	bracketing_test bracketing_complex_test \
	gaussquad_test gaussquad_complex_test \
	linesearch_test linesearch_complex_test \
	newtoncotes_test newtoncotes_complex_test \
	ridder_test ridder_complex_test \
	brent_test brent_complex_test \
	secant_test secant_complex_test \
	rungekutta_test rungekutta_complex_test \
	smla_test smla_complex_test \
	rsla_test rsla_complex_test \
	luabbla_test complex_number_test \
	newton_test newton_complex_test \
	stmatrix_test \
	nelmin_test nelmin_complex_test \
	schedule_test \
	spline_test \
	splinelsq_test

LUA := ../../extern/lua-5.4.3
LIBLUA := $(LUA)/install/lib/liblua.a
LIBLUAPATH := $(LUA)/lib

ifeq ($(DMD), dmd)
    # DFLAGS := -w
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif
ifeq ($(DMD), ldmd2)
    DFLAGS := -w -O -release -inline -boundscheck=off
    OF := -of
    DLINKFLAGS := -L-L$(LIBLUAPATH) -L-ldl
    ifeq ($(PLATFORM), macosx)
	DLINKFLAGS += -L-ld_classic
    endif
endif

# DIP1008 allows throwing of exceptions in @nogc code.
# See notes in src/eilmer/makefile.
DFLAGS += -dip1008 -preview=in

# ----------------------------------------------------------------------
# Here begins the list of targets, starting with the top-level actions.
#
# The default target is test, because that is most likely your interest
# if you are invoking this makefile from within the gas models directory.
# ----------------------------------------------------------------------

test: ${TEST_PROGRAMS}
	tclsh nm-package-test.tcl

demo: ${DEMO_PROGRAMS}
	@echo "Demo codes built."

bench: rungekutta_bench_d rungekutta_bench.go rungekutta_bench_cr
	@echo "Start benchmark programs"
	./rungekutta_bench_d
	go run rungekutta_bench.go
	./rungekutta_bench_cr
	@echo "Finished benchmark programs"

clean:
	- rm -f *.o *.obj
	- rm -f ${DEMO_PROGRAMS}
	- rm -f ${TEST_PROGRAMS}
	- rm -f ./rungekutta_bench_d ./rungekutta_bench_cr
	- cd $(LUA); make clean

# Specific targets, for individual artefacts.

$(LIBLUA):
	cd $(LUA); make $(PLATFORM) local PLATFORM=$(PLATFORM)

bbla_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

bbla_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=bbla_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

bracketing_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

bracketing_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=bracketing_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

gaussquad_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

gaussquad_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=gaussquad_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

linesearch_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

linesearch_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=linesearch_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

newtoncotes_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

newtoncotes_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=newtoncotes_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

ridder_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

ridder_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=ridder_test -version=complex_number \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

secant_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

secant_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=ridder_test -version=complex_number \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

rungekutta_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

rungekutta_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=rungekutta_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

brent_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

brent_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=brent_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

schedule_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

smla_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

smla_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=smla_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

rsla_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

rsla_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=rsla_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

luabbla_test: luabbla.d bbla.d $(NTYPES_FILES) number.d $(UTIL_FILES) $(LIBLUA)
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		luabbla.d bbla.d $(NTYPES_FILES) number.d $(UTIL_FILES) $(LIBLUA) $(DLINKFLAGS)

complex_number_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=complex_number_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

newton_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

newton_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=newton_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

stmatrix_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

nelmin_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

nelmin_complex_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=nelmin_test -version=complex_numbers \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

spline_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

splinelsq_test: $(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d
	$(DMD) -of$@ -debug -g -w -dip1008 -version=$@ \
		$(NM_FILES) $(NTYPES_FILES) $(UTIL_DIR)/msg_service.d $(DLINKFLAGS)

# Some demo targets. These are left over from the early days of development.
# They should continue to work and provide a minimal example of how to get
# going with particular modules.

ridder_demo: ridder_demo.d ridder.d bracketing.d $(NTYPES_FILES)
	$(DMD) -g ridder_demo.d ridder.d bracketing.d $(NTYPES_FILES) $(DLINKFLAGS)

linesearch_demo: linesearch_demo.d linesearch.d $(NTYPES_FILES)
	$(DMD) -g linesearch_demo.d linesearch.d $(NTYPES_FILES) $(DLINKFLAGS)

nelmin_demo: nelmin_demo.d nelmin.d $(NTYPES_FILES) number.d
	$(DMD) -g nelmin_demo.d nelmin.d $(NTYPES_FILES) number.d $(DLINKFLAGS)

nelmin_complex_demo: nelmin_demo.d nelmin.d $(NTYPES_FILES) number.d
	$(DMD) -of$@ -g nelmin_demo.d nelmin.d $(NTYPES_FILES) number.d -version=complex_numbers $(DLINKFLAGS)

rungekutta_demo: rungekutta_demo.d rungekutta.d $(NTYPES_FILES)
	$(DMD) -g rungekutta_demo.d rungekutta.d $(NTYPES_FILES) $(DLINKFLAGS)

newtoncotes_demo: newtoncotes_demo.d newtoncotes.d $(NTYPES_FILES)
	$(DMD) -g newtoncotes_demo.d newtoncotes.d $(NTYPES_FILES) $(DLINKFLAGS)

gaussquad_demo: gaussquad_demo.d gaussquad.d $(NTYPES_FILES)
	$(DMD) -g gaussquad_demo.d gaussquad.d $(NTYPES_FILES) $(DLINKFLAGS)

bbla_demo: bbla_demo.d bbla.d $(NTYPES_FILES)
	$(DMD) -g -dip1008 bbla_demo.d bbla.d $(NTYPES_FILES) $(DLINKFLAGS)

# Rules for building the benchmark programs

rungekutta_bench_d: rungekutta_bench.d
	$(DMD) -of$@ -g -O -boundscheck=off rungekutta_bench.d $(DLINKFLAGS)

rungekutta_bench_cr: rungekutta_bench.cr
	crystal build -o$@ --release rungekutta_bench.cr


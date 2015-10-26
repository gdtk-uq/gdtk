# makefile for the Numerical Methods package
# Builds all of the demonstration programs

UTIL_DIR := ../util

include nm_files.mk

DEMO_PROGRAMS := ridder_demo \
	linesearch_demo \
	nelmin_demo \
	rungekutta_demo \
	newtoncotes_demo \
	gaussquad_demo \
	bbla_demo 

TEST_PROGRAMS := bbla_test \
	bracketing_test \
	gaussquad_test \
	linesearch_test \
	newtoncotes_test \
	ridder_test \
	brent_test \
	rungekutta_test

demo: ${DEMO_PROGRAMS}
	echo "Demo codes built."

test: ${TEST_PROGRAMS}
	tclsh nm-package-test.tcl

clean:
	- rm *.o
	- rm ${DEMO_PROGRAMS}
	- rm ${TEST_PROGRAMS}

ridder_demo: ridder_demo.d ridder.d bracketing.d
	dmd -g ridder_demo.d ridder.d bracketing.d

linesearch_demo: linesearch_demo.d linesearch.d
	dmd -g linesearch_demo.d linesearch.d

nelmin_demo: nelmin_demo.d nelmin.d
	dmd -g nelmin_demo.d nelmin.d

rungekutta_demo: rungekutta_demo.d rungekutta.d
	dmd -g rungekutta_demo.d rungekutta.d

newtoncotes_demo: newtoncotes_demo.d newtoncotes.d
	dmd -g newtoncotes_demo.d newtoncotes.d

gaussquad_demo: gaussquad_demo.d gaussquad.d
	dmd -g gaussquad_demo.d gaussquad.d

bbla_demo: bbla_demo.d bbla.d
	dmd -g bbla_demo.d bbla.d

bbla_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -gc -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

bracketing_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -gc -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

gaussquad_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -gc -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

linesearch_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -gc -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

newtoncotes_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -gc -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

ridder_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -gc -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

rungekutta_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -gc -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d

brent_test: $(NM_FILES) $(UTIL_DIR)/msg_service.d
	dmd -of$@ -debug -gc -w -version=$@ \
		$(NM_FILES) $(UTIL_DIR)/msg_service.d


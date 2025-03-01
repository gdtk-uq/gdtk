# object and source files needed to build ceq
# Note that we list just the D-language source files.
#
CEQ_DIR ?= .
CEQ_BUILD_DIR ?= ./build

CEQ_SRC_FILES := $(CEQ_DIR)/ceq.d

CEQ_OBJ_FILES := $(CEQ_BUILD_DIR)/ceq.o \
	$(CEQ_BUILD_DIR)/common.o \
	$(CEQ_BUILD_DIR)/linalg.o \
	$(CEQ_BUILD_DIR)/pt.o \
	$(CEQ_BUILD_DIR)/rhou.o \
	$(CEQ_BUILD_DIR)/ps.o \
	$(CEQ_BUILD_DIR)/rhot.o \
	$(CEQ_BUILD_DIR)/thermo.o

# object and source files needed to build ceq
# Note that we list just the D-language source files.
#
EQC_DIR ?= .
EQC_SRC_FILES := $(EQC_DIR)/eqc.d

EQC_OBJ_FILES := $(EQC_DIR)/eqc.o \
	$(EQC_DIR)/common.o \
	$(EQC_DIR)/linalg.o \
	$(EQC_DIR)/pt.o \
	$(EQC_DIR)/rhou.o \
	$(EQC_DIR)/ps.o \
	$(EQC_DIR)/rhot.o \
	$(EQC_DIR)/thermo.o

# object and source files needed to build ceq
#
CEQ_DIR ?= .
CEQ_FILES := $(CEQ_DIR)/ceq.d \
	$(CEQ_DIR)/ceq.o \
	$(CEQ_DIR)/common.o \
	$(CEQ_DIR)/linalg.o \
	$(CEQ_DIR)/ps.o \
	$(CEQ_DIR)/pt.o \
	$(CEQ_DIR)/rhou.o \
	$(CEQ_DIR)/thermo.o

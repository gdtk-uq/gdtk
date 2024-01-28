# makefile for chicken
#
# Select the default flavour of the build in SELECTED_TARGETS, below.
# The particular program that is compiled will be called chkn-run.
# Note that the gas model is selected at compile time.
#
# PJ 2022-09-11 build and run in the source directory
#    2022-09-19 build and install, much like our other GDTK programs
#
# To build and install CPU flavour for the ideal air gas model:
# $ make install
#
# To build and install GPU flavour for the reacting gas model:
# $ make GPU=1 AB_GAS=1 install

SELECTED_TARGETS := chkn-prep chkn-post
ifeq ($(GPU),1)
  SELECTED_TARGETS += chicken_for_cuda
  ifeq ($(FP32),1)
    SELECTED_TARGETS += chicken_for_cuda_fp32
  endif
else
  SELECTED_TARGETS += chicken_for_cpu
endif

GAS_MODEL ?= IDEAL_AIR
ifeq ($(AB_GAS),1)
  GAS_MODEL := AB_REACTING_GAS
else
  GAS_MODEL := IDEAL_AIR
endif

# The install destination.
INSTALL_DIR ?= $(HOME)/gdtkinst

# The build destination sits locally for present
BUILD_DIR := ../../build

REPO_DIR := $(shell cd ../../; pwd)
ifneq ("$(wildcard $(REPO_DIR)/.git)","")
    # Repository exists
    REVISION_STRING := $(shell git rev-parse --short HEAD)
else
    REVISION_STRING := unknown
endif

CUDA_SOURCE_FILES := main.cu config.cu simulate.cu bcs.cu block.cu cell.cu face.cu \
  flow.cu gas.cu vector3.cu number.cu

# Target for FP64 (double) numbers.
chicken_for_cuda : $(CUDA_SOURCE_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		main.cu > main_with_rev_string.cu
	nvcc -DCUDA -std c++17 -O2 --expt-relaxed-constexpr -lineinfo -o chkn-run -D$(GAS_MODEL) \
		main_with_rev_string.cu -lz -Xcompiler -fopenmp

# Target for FP32 (float) numbers.
chicken_for_cuda_fp32 : $(CUDA_SOURCE_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		main.cu > main_with_rev_string.cu
	nvcc -DCUDA -DFLOAT_NUMBERS -std c++17 -O2 --expt-relaxed-constexpr -lineinfo -o chkn-run-fp32 -D$(GAS_MODEL) \
		main_with_rev_string.cu -lz -Xcompiler -fopenmp

# Target for CPU-only flavour of the code with FP64 (double) numbers.
# -O2 is usual optimization flag
# -O0 -g for debug
chicken_for_cpu : $(CUDA_SOURCE_FILES)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
		main.cu > main_with_rev_string.cu
	nvcc -lineinfo -std c++17 -O2 --expt-relaxed-constexpr -o chkn-run -D$(GAS_MODEL) \
		main_with_rev_string.cu -lz -Xcompiler -fopenmp

# Note that profiling needs to be done via sudo
# because we need access to the hardware counters.
cuda_profile : chicken_for_cuda
	/usr/local/cuda-11.7/bin/ncu --set full --clock-control none --target-processes all \
		-o chicken -f ./chkn-run --job=$(JOB) --binary

# The program for preparing the config, grid and flow files
# for the start of a simulation.
chkn-prep : chkn_prep.py
	cp chkn_prep.py chkn-prep
	chmod +x chkn-prep

# The program for post-processing the simulation files.
chkn-post : chkn_post.py
	cp chkn_post.py chkn-post
	chmod +x chkn-post

install: $(SELECTED_TARGETS)
	- mkdir -p $(INSTALL_DIR)
	- mkdir -p $(BUILD_DIR)/bin
	- mkdir -p $(BUILD_DIR)/lib
	- mkdir -p $(BUILD_DIR)/share
	cp chkn-run $(BUILD_DIR)/bin
	- cp chkn-run-fp32 $(BUILD_DIR)/bin
	cp chkn-prep $(BUILD_DIR)/bin
	cp chkn-post $(BUILD_DIR)/bin
	cp -r ../lib/* $(BUILD_DIR)/lib
	@echo "Installing to $(INSTALL_DIR)"
	cp -r $(BUILD_DIR)/* $(INSTALL_DIR)

clean :
	- rm -f chkn-run chkn-run-fp32 chkn-prep chkn-post main_with_rev_string.cu
	- rm -rf $(BUILD_DIR)/*
	- rm -f chicken.ncu-rep

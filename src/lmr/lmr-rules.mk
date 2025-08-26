# this file gets included by makefile

MAKEDEPS :=
ifeq ($(DMD),ldc2)
	MAKEDEPS := "--makedeps="
else ifeq ($(DMD),dmd)
	MAKEDEPS := "-makedeps="
endif


SRC_DIR = ../
BIN_DIR := $(BUILD_DIR)/bin
LIB_DIR := $(BUILD_DIR)/lib
SHARE_DIR := $(BUILD_DIR)/share
DATA_DIR := $(BUILD_DIR)/data
ETC_DIR := $(BUILD_DIR)/etc


OBJ_DIR := $(BUILD_DIR)/src
REAL_OBJ_DIR := $(OBJ_DIR)/real
COMPLEX_OBJ_DIR := $(OBJ_DIR)/complex
MPI_OBJ_DIR := $(OBJ_DIR)/mpi
COMPLEX_MPI_OBJ_DIR := $(OBJ_DIR)/complex_mpi

LMR_SRCS = $(abspath $(LMR_CORE_FILES) $(LMR_LUA_FILES) \
           $(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
           $(GEOM_FILES) $(GRID_FILES) \
           $(GAS_FILES) $(EQC_SRC_FILES) $(GZIP_FILES) \
           $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
           $(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
           $(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
           $(DYAML_FILES) $(TINYENDIAN_FILES))

ABS_SRC = $(abspath $(SRC_DIR))
LMR_CMD_SRCS = $(abspath $(LMR_CMD_FILES))

REAL_OBJS := $(patsubst $(ABS_SRC)/%.d,$(REAL_OBJ_DIR)/%.o,$(LMR_SRCS))
COMPLEX_OBJS := $(patsubst $(ABS_SRC)/%.d,$(COMPLEX_OBJ_DIR)/%.o,$(LMR_SRCS))
MPI_OBJS := $(patsubst $(ABS_SRC)/%.d,$(MPI_OBJ_DIR)/%.o,$(LMR_SRCS))
COMPLEX_MPI_OBJS := $(patsubst $(ABS_SRC)/%.d,$(COMPLEX_MPI_OBJ_DIR)/%.o,$(LMR_SRCS))

REAL_CMD_OBJS := $(patsubst $(ABS_SRC)/%.d,$(REAL_OBJ_DIR)/%.o,$(LMR_CMD_SRCS))
COMPLEX_CMD_OBJS := $(patsubst $(ABS_SRC)/%.d,$(COMPLEX_OBJ_DIR)/%.o,$(LMR_CMD_SRCS))
MPI_CMD_OBJS := $(patsubst $(ABS_SRC)/%.d,$(MPI_OBJ_DIR)/%.o,$(LMR_CMD_SRCS))
COMPLEX_MPI_CMD_OBJS := $(patsubst $(ABS_SRC)/%.d,$(COMPLEX_MPI_OBJ_DIR)/%.o,$(LMR_CMD_SRCS))

ifeq ($(MPI_IMPLEMENTATION),OpenMPI)
	MPI_FILES_BUILD = $(patsubst ../%,$(OBJ_DIR)/%,$(MPI_FILES))
	MPI_BUILD_DIR = $(patsubst ../%,$(OBJ_DIR)/%,$(MPI_DIR))
endif

COMPLEX_VERSION_FLAGS := $(DVERSION)complex_numbers
MPI_VERSION_FLAGS := $(DVERSION)mpi_parallel -I$(MPI_BUILD_DIR)/source
COMPLEX_MPI_VERSION_FLAGS := $(DVERSION)complex_numbers $(DVERSION)mpi_parallel -I$(MPI_BUILD_DIR)/source

MPI_EXTRA_DEPS := $(MPI_FILES_BUILD)
COMPLEX_MPI_EXTRA_DEPS := $(MPI_FILES_BUILD)

REAL_NUMBER_TYPE := real
COMPLEX_NUMBER_TYPE := complex
MPI_NUMBER_TYPE := real
COMPLEX_MPI_NUMBER_TYPE := complex

REAL_PARALLEL_FLAVOUR := shared
COMPLEX_PARALLEL_FLAVOUR := shared
MPI_PARALLEL_FLAVOUR := $(MPI_IMPLEMENTATION)
COMPLEX_MPI_PARALLEL_FLAVOUR := $(MPI_IMPLEMENTATION)

LUA_DIR := ../../extern/lua-5.4.3
LUA_SRC_DIR = $(LUA_DIR)/src
LUA_BUILD_DIR = $(patsubst ../../%,$(BUILD_DIR)/%,$(LUA_SRC_DIR))
LUA_TARGET = $(PLATFORM)
LIBLUA = $(LUA_BUILD_DIR)/liblua.a

LPEG_SRC_DIR = $(LUA_DIR)/lpeg-1.0.2
LPEG_BUILD_DIR = $(patsubst ../../%,$(BUILD_DIR)/%,$(LPEG_SRC_DIR))
LPEG_TARGET = $(PLATFORM)
LIBLPEG = $(LPEG_BUILD_DIR)/lpeg.so

EQC_SRC_DIR := $(EQC_DIR)
EQC_BUILD_DIR = $(patsubst ../%,$(OBJ_DIR)/%,$(EQC_SRC_DIR))
LIBEQC := $(EQC_BUILD_DIR)/libeqc.a

PY_PROGRAMS_BIN := $(foreach prog,$(PY_PROGRAMS),$(BIN_DIR)/$(prog))
PROGRAMS_BIN := $(foreach prog,$(PROGRAMS) $(SUB_PROGRAMS) $(AUX_PROGRAMS),$(BIN_DIR)/$(prog))
PROGRAMS_SHARE := $(foreach prog,$(AUX_PROGRAMS_SHARE),$(SHARE_DIR)/$(prog))

CIF_DIR := $(GAS_DIR)/species-database/collision-integrals
COLLISION_INTEGRAL_FILES := $(wildcard $(CIF_DIR)/*.lua)

BUILD_DIRS := $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR) $(REAL_OBJ_DIR) $(COMPLEX_OBJ_DIR) $(MPI_OBJ_DIR) $(COMPLEX_MPI_OBJ_DIR) $(LUA_BUILD_DIR) $(LPEG_BUILD_DIR) $(EQC_BUILD_DIR) $(DATA_DIR) $(ETC_DIR) $(SHARE_DIR) $(LMR_CONFIG_BUILD_DIR)

VERSIONS = REAL MPI COMPLEX COMPLEX_MPI

default: $(PROGRAMS) $(SUB_PROGRAMS)

$(BUILD_DIRS):
	@mkdir -p $@

define GEN_COMPILE_RULE_TEMPLATE
$(1)_DEPS := $$(patsubst %.o,%.dep,$$($(1)_OBJS) \
	$$($(1)_OBJ_DIR)/lmr/main.o \
	$$($(1)_OBJ_DIR)/lmr/commands/runsim_gen.o \
	$$($(1)_OBJ_DIR)/lmr/commands/checkjacobian.o \
	$$($(1)_OBJ_DIR)/lmr/lmrconfig.o)

$$($(1)_OBJS) $$($(1)_CMD_OBJS) $$($(1)_OBJ_DIR)/lmr/main.o $$($(1)_OBJ_DIR)/lmr/commands/checkjacobian.o: $$($(1)_OBJ_DIR)/%.o: $(ABS_SRC)/%.d $$($(1)_EXTRA_DEPS) | $$($(1)_OBJ_DIR)
	@echo
	@echo "Compiling $(1): $$< -> $$@"
	@mkdir -p $$(dir $$@)
	$$(DMD) $$(FLAVOUR_FLAGS) $$(DFLAGS) -c $$(OF)$$@ \
		$$(MAKEDEPS)$$(patsubst %.o,%.dep,$$@) \
		$$(DVERSION)newton_krylov $$($(1)_VERSION_FLAGS) $$<

$$($(1)_OBJ_DIR)/lmr/commands/runsim_gen.o: $(SRC_DIR)/lmr/commands/runsim.d | $$($(1)_OBJ_DIR)
	@echo
	@echo "Generating $$(dir $$@)/runsim_gen.d"
	@mkdir -p $$(dir $$@)
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/$$($(1)_PARALLEL_FLAVOUR)/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/$$($(1)_NUMBER_TYPE)/' \
	    $$< > $$(dir $$@)runsim_gen.d
	@echo
	@echo "Compiling runsim_gen.d: $$(dir $$@)runsim_gen.d -> $$@"
	$$(DMD) $$(FLAVOUR_FLAGS) $$(DFLAGS) -c $$(OF)$$@ \
		$$(MAKEDEPS)$$(patsubst %.o,%.dep,$$@) \
		$(DVERSION)newton_krylov $(DVERSION)run_main $$($(1)_VERSION_FLAGS) $$(dir $$@)runsim_gen.d

$$($(1)_OBJ_DIR)/lmr/lmrconfig.o: $(SRC_DIR)/lmr/lmrconfig.d | $$($(1)_OBJ_DIR)
	@mkdir -p $$(dir $$@)
	@echo
	@echo "Generating $$(dir $$@)/lmrconfig.d"
	sed -e 's/PUT_REVISION_STRING_HERE/$$(REVISION_STRING)/' \
	    -e 's/PUT_FULL_REVISION_STRING_HERE/$$(FULL_REVISION_STRING)/' \
	    -e 's/PUT_REVISION_DATE_HERE/$$(REVISION_DATE)/' \
	    -e 's/PUT_COMPILER_NAME_HERE/$$(DMD)/' \
	    -e 's/PUT_BUILD_DATE_HERE/$$(BUILD_DATE)/' \
	    $$< > $$(dir $$@)lmrconfig.d
	@echo
	@echo "Compiling lmrconfig.d: $$(dir $$@)lmrconfig.d -> $$@"
	$$(DMD) $$(FLAVOUR_FLAGS) $$(DFLAGS) -c $$(OF)$$@ \
		$$(MAKEDEPS)$$(patsubst %.o,%.dep,$$@) \
		$(DVERSION)newton_krylov $$($(1)_VERSION_FLAGS) $$(dir $$@)lmrconfig.d

-include $$($(1)_DEPS)
endef

$(foreach ver,$(VERSIONS),$(eval $(call GEN_COMPILE_RULE_TEMPLATE,$(ver))))

define GEN_LIBRARY_RULE_TEMPLATE
$$(LIB$(1)): | $$($(1)_BUILD_DIR)
	ln -sf $$(abspath $$($(1)_SRC_DIR))/* $$($(1)_BUILD_DIR)/
	$$(MAKE) -C $$($(1)_BUILD_DIR) $$($(1)_TARGET)
endef

$(foreach ver,LUA LPEG EQC,$(eval $(call GEN_LIBRARY_RULE_TEMPLATE,$(ver))))

ifeq ($(MPI_IMPLEMENTATION),OpenMPI)
$(MPI_FILES_BUILD):
	cp -r $(MPI_DIR) $(MPI_BUILD_DIR)
	$(MAKE) -C $(MPI_BUILD_DIR)
endif

.PHONY: lib
lib: $(LIB_DIR)/lua_helper.lua
$(LIB_DIR)/lua_helper.lua:
	@cp -r ../lib $(BUILD_DIR)

.PHONY: share
share: $(SHARE_DIR)/diagnostics-term.gplot
$(SHARE_DIR)/diagnostics-term.gplot:
	@cp -r ./share $(BUILD_DIR)

.PHONY: lua
lua: $(LIBLUA) $(LIBLPEG) lib share | $(BIN_DIR) $(LIB_DIR)
	@cp $(LUA_BUILD_DIR)/dgd-lua $(LUA_BUILD_DIR)/dgd-luac $(BIN_DIR)/
	@cp $(LPEG_BUILD_DIR)/lpeg.so lua-modules/*.lua $(NML_LUA_MODULES) $(LIB_DIR)/
	@cp $(KINETICS_DIR)/mechanism.lua $(KINETICS_DIR)/lex_elems.lua $(KINETICS_DIR)/reaction.lua $(LIB_DIR)/

.PHONY: eqc
eqc: $(LIBEQC)

define DVERS
$(foreach _dv,$(strip $1),$(DVERSION)$(_dv))
endef

define BINARY_TEMPLATE
.PHONY: $(1)
$(1): $(BIN_DIR)/$(1)

$(BIN_DIR)/$(1): $$(strip $$($(1)_OBJS)) $$(strip $$($(1)_EXTRA_PREREQS))
	@echo
	@echo "Linking $(1): $$@"
	$$(DMD) $$(FLAVOUR_FLAGS) $$(DFLAGS) $$(OF)$$@ \
		$$(call DVERS,$$($(1)_DVERS)) $$(DVERSION)newton_krylov $$(strip $$($(1)_PRELINK)) $$^ $$(strip $$($(1)_POSTLINK)) $$(DLINKFLAGS)
endef

TARGETS := \
	lmr \
	lmr-run \
	lmrZ-run \
	lmr-mpi-run \
	lmrZ-mpi-run \
	lmr-check-jacobian \
	lmrZ-check-jacobian

lmr_OBJS := $(REAL_OBJ_DIR)/lmr/main.o $(REAL_CMD_OBJS) $(REAL_OBJS) $(REAL_OBJ_DIR)/lmr/lmrconfig.o
lmr_EXTRA_PREREQS := $(LIBEQC) $(LIBLUA)

lmr-run_OBJS := $(REAL_OBJ_DIR)/lmr/commands/runsim_gen.o $(REAL_OBJ_DIR)/lmr/commands/command.o $(REAL_OBJ_DIR)/lmr/commands/cmdhelper.o $(REAL_OBJS) $(REAL_OBJ_DIR)/lmr/lmrconfig.o
lmr-run_EXTRA_PREREQS := $(LIBLUA) $(LIBEQC)
lmr-run_DVERS := run_main

lmrZ-run_OBJS := $(COMPLEX_OBJ_DIR)/lmr/commands/runsim_gen.o $(COMPLEX_OBJS) $(COMPLEX_OBJ_DIR)/lmr/commands/command.o $(COMPLEX_OBJ_DIR)/lmr/commands/cmdhelper.o $(COMPLEX_OBJ_DIR)/lmr/lmrconfig.o
lmrZ-run_EXTRA_PREREQS := $(LIBEQC) $(LIBLUA)
lmrZ-run_DVERS := run_main complex_numbers

lmr-mpi-run_OBJS := $(MPI_OBJ_DIR)/lmr/commands/runsim_gen.o $(MPI_OBJS) $(MPI_OBJ_DIR)/lmr/commands/command.o $(MPI_OBJ_DIR)/lmr/commands/cmdhelper.o $(MPI_OBJ_DIR)/lmr/lmrconfig.o
lmr-mpi-run_EXTRA_PREREQS := $(LIBEQC) $(LIBLUA) $(MPI_FILES_BUILD)
lmr-mpi-run_DVERS := run_main mpi_parallel
lmr-mpi-run_PRELINK := -I$(BUILD_DIR)/src/extern/OpenMPI/source
lmr-mpi-run_POSTLINK := -L-lmpi

lmrZ-mpi-run_OBJS := $(COMPLEX_MPI_OBJ_DIR)/lmr/commands/runsim_gen.o $(COMPLEX_MPI_OBJS) $(COMPLEX_MPI_OBJ_DIR)/lmr/commands/command.o $(COMPLEX_MPI_OBJ_DIR)/lmr/commands/cmdhelper.o $(COMPLEX_MPI_OBJ_DIR)/lmr/lmrconfig.o
lmrZ-mpi-run_EXTRA_PREREQS := $(LIBEQC) $(LIBLUA) $(MPI_FILES_BUILD)
lmrZ-mpi-run_DVERS := run_main mpi_parallel complex_numbers
lmrZ-mpi-run_PRELINK := -I$(BUILD_DIR)/src/extern/OpenMPI/source
lmrZ-mpi-run_POSTLINK := -L-lmpi

lmr-check-jacobian_OBJS := $(REAL_OBJS) $(REAL_OBJ_DIR)/lmr/commands/checkjacobian.o $(REAL_CMD_OBJS) $(REAL_OBJ_DIR)/lmr/lmrconfig.o
lmr-check-jacobian_EXTRA_PREREQS := $(LMR_CMD)/checkjacobian.d $(LIBEQC) $(LIBLUA)

lmrZ-check-jacobian_OBJS := $(COMPLEX_OBJS) $(COMPLEX_OBJ_DIR)/lmr/commands/checkjacobian.o $(COMPLEX_CMD_OBJS) $(COMPLEX_OBJ_DIR)/lmr/lmrconfig.o
lmrZ-check-jacobian_EXTRA_PREREQS := $(LIBEQC) $(LIBLUA)
lmrZ-check-jacobian_DVERS := complex_numbers

$(foreach t,$(TARGETS),$(eval $(call BINARY_TEMPLATE,$(t))))

# if FLAVOUR=fast, then all of the build .o files will use that.
# this means that if the flavour is not debug, then the best way to additionally build lmr-debug
# is from scratch, giving it all relevant source files
.PHONY: lmr-debug
lmr-debug: $(BIN_DIR)/lmr-debug
ifneq ($(FLAVOUR), debug)
$(BIN_DIR)/lmr-debug: main.d $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(EQC_SRC_FILES) $(LIBEQC) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	$(DMD) $(DEBUG_DFLAGS) $(DFLAGS) -od=$(OBJ_DIR)/lmr $(OF)$@ \
		$(DVERSION)newton_krylov \
		main.d \
		$(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d \
		$(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
		$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
		$(GEOM_FILES) $(GRID_FILES) \
		$(GAS_FILES) $(EQC_SRC_FILES) $(GZIP_FILES) \
		$(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
		$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
		$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
		$(LIBEQC) $(LIBLUA) \
		$(DYAML_FILES) $(TINYENDIAN_FILES) \
		$(DLINKFLAGS)
else
$(BIN_DIR)/lmr-debug: $(BIN_DIR)/lmr
	@mkdir -p $(dir $@)
	cp $< $@
endif
	
$(PY_PROGRAMS_BIN): $(BIN_DIR)/lmr-%: python-programs/lmr_%.py
	@mkdir -p $(BIN_DIR)
	cp $< $@
	chmod +x $@

$(addprefix $(ETC_DIR)/,$(ETC_FILES)): $(ETC_FILES)
	@mkdir -p $(ETC_DIR)
	cp $< $@

.PHONY: lmr-complete
lmr-complete: $(SHARE_DIR)/lmr-complete.sh
$(SHARE_DIR)/lmr-complete.sh:
	echo '#/usr/bin/env bash\n' \
		'_lmr_completions()\n'\
		'{\n'\
		'COMPREPLY=($$(if [ "$$COMP_CWORD" -eq 1 ]; then\n compgen -W "compute-norms prep-gas slice-flow prep-grids prep-mapped-cells probe-flow prep-reactions limiter2vtk snapshot2vtk prep-sim residual2vtk plot-diagnostics list-species slice-solid custom-script prep-energy-exchange gradient2vtk run revision-id structured2unstructured extract-line lmrZ-check-jacobian lmr-check-jacobian" "$${COMP_WORDS[1]}"\n fi))\n'\
		'}\n'\
		'complete -F _lmr_completions lmr\n' > $@

.PHONY: prep-gas
prep-gas: $(BIN_DIR)/prep-gas
$(BIN_DIR)/prep-gas: $(DATA_DIR)/species-database.lua $(abspath $(GAS_DIR)/species_data_converter.lua $(COLLISION_INTEGRAL_FILES) $(GAS_DIR)/prep_gas.lua)
	@mkdir -p $(BIN_DIR)
	cp $(GAS_DIR)/prep_gas.lua $(BIN_DIR)/prep-gas; chmod +x $(BIN_DIR)/prep-gas
	cp $(GAS_DIR)/species_data_converter.lua $(BIN_DIR)/species-data-converter; \
		chmod +x $(BIN_DIR)/species-data-converter
	@mkdir -p $(DATA_DIR)
	cp $(COLLISION_INTEGRAL_FILES) $(DATA_DIR)/

.PHONY: species-database
species-database: $(DATA_DIR)/species-database.lua
$(DATA_DIR)/species-database.lua:
	@mkdir -p $(BUILD_DIR)/src/gas/species-database
	@ln -sf $(abspath $(abspath $(GAS_DIR)/species-database)/*) $(BUILD_DIR)/src/gas/species-database/
	$(MAKE) -C $(BUILD_DIR)/src/gas/species-database species-database.lua
	@mkdir -p $(DATA_DIR)
	@cp $(BUILD_DIR)/src/gas/species-database/species-list.txt $(DATA_DIR)/
	@cp $(BUILD_DIR)/src/gas/species-database/species-database.lua $(DATA_DIR)/

.PHONY: ugrid_partition
ugrid_partition: $(BIN_DIR)/ugrid_partition
$(BIN_DIR)/ugrid_partition: $(GRID_DIR)/ugrid_partition.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(OBJ_DIR)/lmr $(OF)$@ $< $(DLINKFLAGS)

KIN_PROGS = prep-chem chemkin2eilmer prep-kinetics
prep-chem_lua = $(KINETICS_DIR)/prep_chem.lua
chemkin2eilmer_lua = $(KINETICS_DIR)/chemkin2eilmer.lua
prep-kinetics_lua = $(KINETICS_DIR)/prep_kinetics.lua

define GEN_KINETICS_PROGRAM_RULE_TEMPLATE
.PHONY: $(1)
$(1): $$(BIN_DIR)/$(1)
$$(BIN_DIR)/$(1): $$($(1)_lua)
	@echo "Installing kinetics program: $$< -> $$@"
	@mkdir -p $$(BIN_DIR)
	cp $$< $$@; chmod +x $$@
endef

$(foreach prog,$(KIN_PROGS),$(eval $(call GEN_KINETICS_PROGRAM_RULE_TEMPLATE,$(prog))))

.PHONY: gdtk-module
gdtk-module: $(SHARE_DIR)/gdtk-module
$(SHARE_DIR)/gdtk-module:
	@mkdir -p $(SHARE_DIR)
	sed -e 's+PUT_REVISION_STRING_HERE+$(REVISION_STRING)+' \
	    -e 's+PUT_COMPILER_NAME_HERE+$(DMD)+' \
	    -e 's+PUT_INSTALL_DIR_HERE+$(INSTALL_DIR)+' \
	    -e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
	    -e 's+PUT_BUILD_DATE_HERE+$(BUILD_DATE)+' \
	    -e 's+PUT_REVISION_AGE_HERE+$(REVISION_AGE)+' \
	    -e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
	    ../eilmer/gdtk-module-template > $@

install: $(PROGRAMS_BIN) $(PY_PROGRAMS_BIN) $(PROGRAMS_SHARE) $(addprefix $(ETC_DIR)/,$(ETC_FILES)) lib share lua | $(BIN_DIR) $(LIB_DIR) $(ETC_DIR) $(DATA_DIR) $(SHARE_DIR)
	@mkdir -p $(INSTALL_DIR)
	@cp -r $(BIN_DIR) $(INSTALL_DIR)
	@cp -r $(LIB_DIR) $(INSTALL_DIR)
	@cp -r $(ETC_DIR) $(INSTALL_DIR)
	@cp -r $(DATA_DIR) $(INSTALL_DIR)
	@cp -r $(SHARE_DIR) $(INSTALL_DIR)

clean:
	- rm -rf $(BUILD_DIR)

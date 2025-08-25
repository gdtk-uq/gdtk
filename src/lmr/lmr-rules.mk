# this file gets included by makefile

# TODO: S. Imran 26-Aug-2025: simplify this stuff, there are many more variables
# than actually needed

LMR_BUILD_DIR = $(abspath $(BUILD_DIR))
LMR_OBJ_DIR := $(abspath $(LMR_BUILD_DIR)/src)
LMR_CONFIG_BUILD_DIR = $(LMR_BUILD_DIR)/config
SRC_DIR = $(abspath ../)

# CORE SOURCE FILES (everything except commands)
_LMR_CORE_SRCS_RAW = $(LMR_CORE_FILES) $(LMR_LUA_FILES) \
           $(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
           $(GEOM_FILES) $(GRID_FILES) \
           $(GAS_FILES) $(EQC_SRC_FILES) $(GZIP_FILES) \
           $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
           $(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
           $(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
           $(DYAML_FILES) $(TINYENDIAN_FILES)
_LMR_CORE_SRCS_ABS := $(abspath $(_LMR_CORE_SRCS_RAW))
LMR_CORE_SRCS := $(patsubst $(SRC_DIR)/%,%,$(_LMR_CORE_SRCS_ABS))

_LMR_CMD_SRCS_RAW = $(LMR_CMD_FILES)
_LMR_CMD_SRCS_ABS := $(abspath $(_LMR_CMD_SRCS_RAW))
LMR_CMD_SRCS := $(patsubst $(SRC_DIR)/%,%,$(_LMR_CMD_SRCS_ABS))

_ALL_SRCS_ABS := $(_LMR_CORE_SRCS_ABS) $(_LMR_CMD_SRCS_ABS)
CONTAINS_COMPLEX := $(shell grep -l "version(complex_numbers)" $(_ALL_SRCS_ABS))
CONTAINS_MPI := $(shell grep -l "version(mpi_parallel)" $(_ALL_SRCS_ABS))

LMR_BIN_DIR := $(LMR_BUILD_DIR)/bin

LMR_REAL_OBJ_DIR := $(LMR_OBJ_DIR)/real
LMR_COMPLEX_OBJ_DIR := $(LMR_OBJ_DIR)/complex
LMR_MPI_OBJ_DIR := $(LMR_OBJ_DIR)/mpi
LMR_COMPLEX_MPI_OBJ_DIR := $(LMR_OBJ_DIR)/complex_mpi

LMR_CORE_REAL_OBJS := $(patsubst %.d,$(LMR_REAL_OBJ_DIR)/%.o,$(LMR_CORE_SRCS))
LMR_CORE_COMPLEX_OBJS := $(patsubst %.d,$(LMR_COMPLEX_OBJ_DIR)/%.o,$(LMR_CORE_SRCS))
LMR_CORE_MPI_OBJS := $(patsubst %.d,$(LMR_MPI_OBJ_DIR)/%.o,$(LMR_CORE_SRCS))
LMR_CORE_COMPLEX_MPI_OBJS := $(patsubst %.d,$(LMR_COMPLEX_MPI_OBJ_DIR)/%.o,$(LMR_CORE_SRCS))

LMR_CMD_REAL_OBJS := $(patsubst %.d,$(LMR_REAL_OBJ_DIR)/%.o,$(LMR_CMD_SRCS))
LMR_CMD_COMPLEX_OBJS := $(patsubst %.d,$(LMR_COMPLEX_OBJ_DIR)/%.o,$(LMR_CMD_SRCS))
LMR_CMD_MPI_OBJS := $(patsubst %.d,$(LMR_MPI_OBJ_DIR)/%.o,$(LMR_CMD_SRCS))
LMR_CMD_COMPLEX_MPI_OBJS := $(patsubst %.d,$(LMR_COMPLEX_MPI_OBJ_DIR)/%.o,$(LMR_CMD_SRCS))

LMR_REAL_OBJS := $(LMR_CORE_REAL_OBJS) $(LMR_CMD_REAL_OBJS)
LMR_COMPLEX_OBJS := $(LMR_CORE_COMPLEX_OBJS) $(LMR_CMD_COMPLEX_OBJS)
LMR_MPI_OBJS := $(LMR_CORE_MPI_OBJS) $(LMR_CMD_MPI_OBJS)
LMR_COMPLEX_MPI_OBJS := $(LMR_CORE_COMPLEX_MPI_OBJS) $(LMR_CMD_COMPLEX_MPI_OBJS)

LUA_DIR := $(abspath ../../extern/lua-5.4.3)
LUA_BUILD_ROOT = ${BUILD_DIR}/extern/lua-5.4.3
LUA_SRC_DIR = $(LUA_DIR)/src
LPEG_SRC_DIR = $(LUA_DIR)/lpeg-1.0.2
LUA_BUILD_DIR = $(LUA_BUILD_ROOT)/src
LPEG_BUILD_DIR = $(LUA_BUILD_ROOT)/lpeg-1.0.2

LIBLUA = $(LUA_BUILD_DIR)/liblua.a
LIBLPEG = $(LPEG_BUILD_DIR)/lpeg.so

EQC_BUILD_DIR = $(LMR_OBJ_DIR)/extern/eqc
LIBEQC := $(EQC_BUILD_DIR)/libeqc.a

PY_PROGRAMS_BIN := $(foreach prog,$(PY_PROGRAMS),$(LMR_BUILD_DIR)/bin/$(prog))
LMR_PROGRAMS_BIN := $(foreach prog,$(PROGRAMS) $(SUB_PROGRAMS) $(AUX_PROGRAMS),$(LMR_BUILD_DIR)/bin/$(prog))
LMR_PROGRAMS_BIN += $(PY_PROGRAMS_BIN)
LMR_PROGRAMS_SHARE := $(foreach prog,$(AUX_PROGRAMS_SHARE),$(LMR_BUILD_DIR)/share/$(prog))

ifeq ($(MPI_IMPLEMENTATION),OpenMPI)
	MPI_FILES_BUILD = $(patsubst $(SRC_DIR)/%.d,$(LMR_OBJ_DIR)/%.d,$(abspath $(MPI_FILES)))
endif

CIF_DIR := $(GAS_DIR)/species-database/collision-integrals
COLLISION_INTEGRAL_FILES := $(CIF_DIR)/gupta_etal_1990_CI_data.lua \
	$(CIF_DIR)/wright_etal_CI_data.lua \
	$(CIF_DIR)/palmer_etal_CI_data.lua

BUILD_DIRS := $(LMR_BIN_DIR) $(LMR_REAL_OBJ_DIR) $(LMR_COMPLEX_OBJ_DIR) $(LMR_MPI_OBJ_DIR) $(LMR_COMPLEX_MPI_OBJ_DIR)

VPATH = $(SRC_DIR)

default: $(PROGRAMS) $(SUB_PROGRAMS)

$(BUILD_DIRS):
	@mkdir -p $@

$(LMR_REAL_OBJS): $(LMR_REAL_OBJ_DIR)/%.o: %.d | $(LMR_REAL_OBJ_DIR)
	@echo "Compiling real: $< -> $@"
	@mkdir -p $(dir $@)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -c -of=$@ $(DVERSION)newton_krylov $<

$(LMR_COMPLEX_OBJS): $(LMR_COMPLEX_OBJ_DIR)/%.o: %.d | $(LMR_COMPLEX_OBJ_DIR)
	@echo "Compiling complex: $< -> $@"
	@mkdir -p $(dir $@)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -c -of=$@ $(DVERSION)newton_krylov $(DVERSION)complex_numbers $<

$(LMR_MPI_OBJS): $(LMR_MPI_OBJ_DIR)/%.o: %.d $(MPI_FILES_BUILD) | $(LMR_MPI_OBJ_DIR)
	@echo "Compiling mpi: $< -> $@"
	@mkdir -p $(dir $@)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -c -of=$@ $(DVERSION)newton_krylov $(DVERSION)mpi_parallel -I$(LMR_BUILD_DIR)/src/extern/OpenMPI/source $<

$(LMR_COMPLEX_MPI_OBJS): $(LMR_COMPLEX_MPI_OBJ_DIR)/%.o: %.d $(MPI_FILES_BUILD) | $(LMR_COMPLEX_MPI_OBJ_DIR)
	@echo "Compiling complex+mpi: $< -> $@"
	@mkdir -p $(dir $@)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -c -of=$@ $(DVERSION)newton_krylov $(DVERSION)complex_numbers $(DVERSION)mpi_parallel -I$(LMR_BUILD_DIR)/src/extern/OpenMPI/source $<

$(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d: lmrconfig.d
	@mkdir -p $(dir $@)
	sed -e 's/PUT_REVISION_STRING_HERE/$(REVISION_STRING)/' \
	    -e 's/PUT_FULL_REVISION_STRING_HERE/$(FULL_REVISION_STRING)/' \
	    -e 's/PUT_REVISION_DATE_HERE/$(REVISION_DATE)/' \
	    -e 's/PUT_COMPILER_NAME_HERE/$(DMD)/' \
	    -e 's/PUT_BUILD_DATE_HERE/$(BUILD_DATE)/' \
	    $< > $@

$(LMR_CONFIG_BUILD_DIR)/runsim_shared.d: $(LMR_CMD)/runsim.d
	@mkdir -p $(dir $@)
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/shared/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/real/' \
	    $< > $@

$(LMR_CONFIG_BUILD_DIR)/runsim_shared_Z.d: $(LMR_CMD)/runsim.d
	@mkdir -p $(dir $@)
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/shared/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/complex/' \
	    $< > $@

$(LMR_CONFIG_BUILD_DIR)/runsim_mpi.d: $(LMR_CMD)/runsim.d $(MPI_FILES_BUILD)
	@mkdir -p $(dir $@)
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/$(MPI_IMPLEMENTATION)/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/real/' \
	    $< > $@

$(LMR_CONFIG_BUILD_DIR)/runsim_mpi_Z.d: $(LMR_CMD)/runsim.d $(MPI_FILES_BUILD)
	@mkdir -p $(dir $@)
	sed -e 's/PUT_PARALLEL_FLAVOUR_HERE/$(MPI_IMPLEMENTATION)/' \
	    -e 's/PUT_NUMBER_TYPE_HERE/complex/' \
	    $< > $@

$(LUA_BUILD_DIR) $(LPEG_BUILD_DIR) $(EQC_BUILD_DIR):
	mkdir -p $@

$(LIBLPEG): | $(LPEG_BUILD_DIR)
	ln -sf $(LPEG_SRC_DIR)/*.c $(LPEG_BUILD_DIR)
	ln -sf $(LPEG_SRC_DIR)/*.h $(LPEG_BUILD_DIR)
	ln -sf $(LPEG_SRC_DIR)/makefile $(LPEG_BUILD_DIR)
	$(MAKE) -C $(LPEG_BUILD_DIR) $(PLATFORM)

$(LIBLUA): | $(LUA_BUILD_DIR)
	ln -sf $(LUA_SRC_DIR)/*.c $(LUA_BUILD_DIR)
	ln -sf $(LUA_SRC_DIR)/*.h $(LUA_BUILD_DIR)
	ln -sf $(LUA_SRC_DIR)/Makefile $(LUA_BUILD_DIR)
	$(MAKE) -C $(LUA_BUILD_DIR) $(PLATFORM)

$(LIBEQC): | $(EQC_BUILD_DIR)
	ln -sf $(abspath $(EQC_DIR)/*.c) $(EQC_BUILD_DIR)
	ln -sf $(abspath $(EQC_DIR)/*.h) $(EQC_BUILD_DIR)
	ln -sf $(abspath $(EQC_DIR)/makefile) $(EQC_BUILD_DIR)
	$(MAKE) -C $(EQC_BUILD_DIR)

$(MPI_FILES_BUILD):
	cp -r $(MPI_DIR) $(LMR_BUILD_DIR)/src/extern
	$(MAKE) -C $(LMR_BUILD_DIR)/src/extern/OpenMPI


.PHONY: lua
lua: $(LIBLUA) $(LIBLPEG)
	-mkdir -p $(LMR_BUILD_DIR)/bin
	-mkdir -p $(LMR_BUILD_DIR)/lib
	@cp $(LUA_BUILD_DIR)/dgd-lua  $(LMR_BUILD_DIR)/bin
	@cp $(LUA_BUILD_DIR)/dgd-luac  $(LMR_BUILD_DIR)/bin
	@cp $(LPEG_BUILD_DIR)/lpeg.so  $(LMR_BUILD_DIR)/lib
	cp lua-modules/*.lua $(LMR_BUILD_DIR)/lib/
	cp $(NML_LUA_MODULES) $(LMR_BUILD_DIR)/lib/
	cp $(KINETICS_DIR)/mechanism.lua $(KINETICS_DIR)/lex_elems.lua $(KINETICS_DIR)/reaction.lua $(LMR_BUILD_DIR)/lib/

.PHONY: eqc
eqc: $(LIBEQC)

.PHONY: lmr
lmr: $(LMR_BUILD_DIR)/bin/lmr
$(LMR_BUILD_DIR)/bin/lmr: main.d $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_REAL_OBJS) $(LMR_CMD_REAL_OBJS) $(LIBEQC) $(LIBLUA)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(LMR_OBJ_DIR)/lmr $(OF)$@ \
		$(DVERSION)newton_krylov $^ $(DLINKFLAGS)

.PHONY: lmr-run
lmr-run: $(LMR_BUILD_DIR)/bin/lmr-run
$(LMR_BUILD_DIR)/bin/lmr-run: $(LMR_CONFIG_BUILD_DIR)/runsim_shared.d $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_REAL_OBJS) $(filter-out %/runsim.o,$(LMR_CMD_REAL_OBJS)) $(LIBLUA) $(LIBEQC)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(LMR_OBJ_DIR)/lmr $(OF)$@ \
		$(DVERSION)run_main $(DVERSION)newton_krylov \
		$^ $(DLINKFLAGS)

.PHONY: lmrZ-run
lmrZ-run: $(LMR_BUILD_DIR)/bin/lmrZ-run
$(LMR_BUILD_DIR)/bin/lmrZ-run: $(LMR_CONFIG_BUILD_DIR)/runsim_shared_Z.d $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_COMPLEX_OBJS) $(filter-out %/runsim.o,$(LMR_CMD_COMPLEX_OBJS)) $(LIBEQC) $(LIBLUA)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(LMR_OBJ_DIR)/lmr $(OF)$@ \
		$(DVERSION)run_main $(DVERSION)newton_krylov $(DVERSION)complex_numbers \
		$^ $(DLINKFLAGS)

.PHONY: lmr-mpi-run
lmr-mpi-run: $(LMR_BUILD_DIR)/bin/lmr-mpi-run
$(LMR_BUILD_DIR)/bin/lmr-mpi-run: $(LMR_CONFIG_BUILD_DIR)/runsim_mpi.d $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_MPI_OBJS) $(filter-out %/runsim.o,$(LMR_CMD_MPI_OBJS)) $(LIBEQC) $(LIBLUA) $(MPI_FILES_BUILD)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(LMR_OBJ_DIR)/lmr $(OF)$@ $(DVERSION)run_main \
		$(DVERSION)newton_krylov $(DVERSION)mpi_parallel $(MPI_LIB_DIRS_SEARCH) \
		-I$(LMR_BUILD_DIR)/src/extern/OpenMPI/source $^ -L-lmpi $(DLINKFLAGS) 

.PHONY: lmrZ-mpi-run
lmrZ-mpi-run: $(LMR_BUILD_DIR)/bin/lmrZ-mpi-run
$(LMR_BUILD_DIR)/bin/lmrZ-mpi-run: $(LMR_CONFIG_BUILD_DIR)/runsim_mpi_Z.d $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_COMPLEX_MPI_OBJS) $(filter-out %/runsim.o,$(LMR_CMD_COMPLEX_MPI_OBJS)) $(LIBEQC) $(LIBLUA) $(MPI_FILES_BUILD)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(LMR_OBJ_DIR)/lmr $(OF)$@ $(DVERSION)run_main \
		$(DVERSION)newton_krylov $(DVERSION)mpi_parallel $(MPI_LIB_DIRS_SEARCH) \
		$(DVERSION)complex_numbers \
		-I$(LMR_BUILD_DIR)/src/extern/OpenMPI/source $^ -L-lmpi $(DLINKFLAGS)

.PHONY: lmr-check-jacobian
lmr-check-jacobian: $(LMR_BUILD_DIR)/bin/lmr-check-jacobian
$(LMR_BUILD_DIR)/bin/lmr-check-jacobian: $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_REAL_OBJS) $(filter-out %/checkjacobian.o,$(LMR_CMD_REAL_OBJS)) $(LMR_CMD)/checkjacobian.d $(LIBEQC) $(LIBLUA)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(LMR_REAL_OBJ_DIR) $(OF)$@ \
		$(DVERSION)newton_krylov $^ $(DLINKFLAGS)

.PHONY: lmrZ-check-jacobian
lmrZ-check-jacobian: $(LMR_BUILD_DIR)/bin/lmrZ-check-jacobian
$(LMR_BUILD_DIR)/bin/lmrZ-check-jacobian: $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_COMPLEX_OBJS) $(filter-out %/checkjacobian.o,$(LMR_CMD_COMPLEX_OBJS)) $(LMR_CMD)/checkjacobian.d $(LIBEQC) $(LIBLUA)
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(LMR_COMPLEX_OBJ_DIR) $(OF)$@ \
		$(DVERSION)newton_krylov $(DVERSION)complex_numbers $^ $(DLINKFLAGS)

# if FLAVOUR=fast, then all of the build .o files will use that.
# this means that if the flavour is not debug, then the best way to additionally build lmr-debug
# is from scratch, giving it all relevant source files
.PHONY: lmr-debug
lmr-debug: $(LMR_BUILD_DIR)/bin/lmr-debug
ifneq ($(FLAVOUR), debug)
$(LMR_BUILD_DIR)/bin/lmr-debug: main.d $(LMR_CONFIG_BUILD_DIR)/lmrconfig_with_str_subst.d $(LMR_CORE_FILES) $(LMR_CMD_FILES) $(LMR_LUA_FILES) \
	$(LMR_BC_FILES) $(LMR_SOLID_FILES) $(LMR_EFIELD_FILES) \
	$(GEOM_FILES) $(GRID_FILES) \
	$(GAS_FILES) $(EQC_SRC_FILES) $(LIBEQC) $(LIBLUA) $(GZIP_FILES) \
	$(KINETICS_FILES) $(GAS_LUA_FILES) $(KINETICS_LUA_FILES) \
	$(NM_FILES) $(NTYPES_FILES) $(UTIL_FILES) \
	$(GASDYN_FILES) $(GASDYN_LUA_FILES) $(NM_LUA_FILES) \
	$(DYAML_FILES) $(TINYENDIAN_FILES)
	$(DMD) $(DEBUG_DFLAGS) $(DFLAGS) -od=$(LMR_OBJ_DIR)/lmr $(OF)$@ \
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
$(LMR_BUILD_DIR)/bin/lmr-debug: $(LMR_BUILD_DIR)/bin/lmr
	@mkdir -p $(dir $@)
	cp $< $@
endif
	
$(PY_PROGRAMS_BIN): $(LMR_BUILD_DIR)/bin/lmr-%: python-programs/lmr_%.py
	@mkdir -p $(LMR_BUILD_DIR)/bin
	cp $< $@
	chmod +x $@

$(addprefix $(LMR_BUILD_DIR)/etc/,$(ETC_FILES)): $(ETC_FILES)
	@mkdir -p $(LMR_BUILD_DIR)/etc
	cp $< $@

.PHONY: lmr-complete
lmr-complete: $(LMR_BUILD_DIR)/share/lmr-complete.sh
$(LMR_BUILD_DIR)/share/lmr-complete.sh:
	echo '#/usr/bin/env bash\n' \
		'_lmr_completions()\n'\
		'{\n'\
		'COMPREPLY=($$(if [ "$$COMP_CWORD" -eq 1 ]; then\n compgen -W "compute-norms prep-gas slice-flow prep-grids prep-mapped-cells probe-flow prep-reactions limiter2vtk snapshot2vtk prep-sim residual2vtk plot-diagnostics list-species slice-solid custom-script prep-energy-exchange gradient2vtk run revision-id structured2unstructured extract-line lmrZ-check-jacobian lmr-check-jacobian" "$${COMP_WORDS[1]}"\n fi))\n'\
		'}\n'\
		'complete -F _lmr_completions lmr\n' > $@

.PHONY: prep-gas
prep-gas: $(LMR_BUILD_DIR)/bin/prep-gas
$(LMR_BUILD_DIR)/bin/prep-gas: $(GAS_DIR)/prep_gas.lua $(LMR_BUILD_DIR)/data/species-database.lua \
		$(GAS_DIR)/species-database/species-list.txt $(GAS_DIR)/species_data_converter.lua $(COLLISION_INTEGRAL_FILES)
	@mkdir -p $(LMR_BUILD_DIR)/bin
	cp $(GAS_DIR)/prep_gas.lua $(LMR_BUILD_DIR)/bin/prep-gas; chmod +x $(LMR_BUILD_DIR)/bin/prep-gas
	cp $(GAS_DIR)/species_data_converter.lua $(LMR_BUILD_DIR)/bin/species-data-converter; \
		chmod +x $(LMR_BUILD_DIR)/bin/species-data-converter
	@mkdir -p $(LMR_BUILD_DIR)/data
	cp $(GAS_DIR)/species-database/species-list.txt $(LMR_BUILD_DIR)/data/
	cp $(COLLISION_INTEGRAL_FILES) $(LMR_BUILD_DIR)/data/

$(LMR_BUILD_DIR)/data/species-database.lua:
	@mkdir $(LMR_BUILD_DIR)/data
	$(MAKE) SPECIES_DATABASE_DIR=$(LMR_BUILD_DIR)/data -C $(GAS_DIR)/species-database

.PHONY: ugrid_partition
ugrid_partition: $(LMR_BUILD_DIR)/bin/ugrid_partition
$(LMR_BUILD_DIR)/bin/ugrid_partition: $(GRID_DIR)/ugrid_partition.d
	$(DMD) $(FLAVOUR_FLAGS) $(DFLAGS) -od=$(LMR_OBJ_DIR)/lmr $(OF)$@ $< $(DLINKFLAGS)
	

.PHONY: prep-chem
prep-chem: $(LMR_BUILD_DIR)/bin/prep-chem
$(LMR_BUILD_DIR)/bin/prep-chem: $(KINETICS_DIR)/prep_chem.lua lua
	- mkdir -p $(LMR_BUILD_DIR)/bin
	cp $< $@; chmod +x $@

.PHONY: chemkin2eilmer
chemkin2eilmer: $(LMR_BUILD_DIR)/bin/chemkin2eilmer
$(LMR_BUILD_DIR)/bin/chemkin2eilmer: $(KINETICS_DIR)/chemkin2eilmer.lua lua
	- mkdir -p $(LMR_BUILD_DIR)/bin
	cp $< $@; chmod +x $@

.PHONY: prep-kinetics
prep-kinetics: $(LMR_BUILD_DIR)/bin/prep-kinetics 
$(LMR_BUILD_DIR)/bin/prep-kinetics: $(KINETICS_DIR)/prep_kinetics.lua lua
	- mkdir -p $(LMR_BUILD_DIR)/bin
	cp $< $@; chmod +x $@

.PHONY: gdtk-module
gdtk-module: $(LMR_BUILD_DIR)/share/gdtk-module
$(LMR_BUILD_DIR)/share/gdtk-module:
	@mkdir -p $(LMR_BUILD_DIR)/share
	sed -e 's+PUT_REVISION_STRING_HERE+$(REVISION_STRING)+' \
	    -e 's+PUT_COMPILER_NAME_HERE+$(DMD)+' \
	    -e 's+PUT_INSTALL_DIR_HERE+$(INSTALL_DIR)+' \
	    -e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
	    -e 's+PUT_BUILD_DATE_HERE+$(BUILD_DATE)+' \
	    -e 's+PUT_REVISION_AGE_HERE+$(REVISION_AGE)+' \
	    -e 's+PUT_REPO_DIR_HERE+$(REPO_DIR)+' \
	    ../eilmer/gdtk-module-template > $@

.PHONY: lib
lib:
	@cp -r ../lib $(LMR_BUILD_DIR)

.PHONY: share
share:
	@cp -r ./share $(LMR_BUILD_DIR)

install: $(LMR_PROGRAMS_BIN) $(LMR_PROGRAMS_SHARE) $(addprefix $(LMR_BUILD_DIR)/etc/,$(ETC_FILES)) lib share
	-mkdir $(INSTALL_DIR)
	cp -r $(LMR_BUILD_DIR)/bin $(INSTALL_DIR)
	cp -r $(LMR_BUILD_DIR)/lib $(INSTALL_DIR)
	cp -r $(LMR_BUILD_DIR)/etc $(INSTALL_DIR)
	cp -r $(LMR_BUILD_DIR)/data $(INSTALL_DIR)
	cp -r $(LMR_BUILD_DIR)/share $(INSTALL_DIR)

clean:
	- rm -rf $(LMR_BUILD_DIR)

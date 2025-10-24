# makefile for pitot3
# This is a newer one which I did with some help from Peter and Rowan.
# Now installs the gas library, CEA, and PITOT3.
# Chris James (c.james4@uq.edu.au) - 12-Jun-2021

INSTALL_DIR ?= $(HOME)/gdtkinst

GDTKPY_DIR := ../gdtk-py
CEA_DIR := ../../extern/cea2

ifeq ($(wildcard $(CEA_DIR)),)
$(warning CEA not found at CEA_DIR=$(CEA_DIR))
endif

install:
	$(info Installing cea)
	- $(MAKE) -C $(CEA_DIR) install
	$(info Installing Pitot3 Python library)
	$(MAKE) pitot3 
	$(MAKE) pitot3-condition-builder 
	@echo -----------------------------------------
	@echo pitot3 should now be installed.
	@echo -----------------------------------------

clean: 
	find . -type d -name '__pycache__' -exec rm -r {} +
	find . -type d -name '*.egg-info' -exec rm -r {} +
	$(MAKE) -sC $(GDTKPY_DIR) clean
	- $(MAKE) -sC $(CEA_DIR) clean

# Priority:
#   1. VIRTUAL_ENV - Standard set by python virtual-environments
#   2. CONDA_PREFIX - Set by (ana)conda environments
#   3. INSTALL_DIR/.venv - Create our own if none already exist
CONDA_PREFIX ?= $(INSTALL_DIR)/.venv
VIRTUAL_ENV ?= $(CONDA_PREFIX)

# Override if using uv with `PIP='uv pip'`
PIP ?= $(VIRTUAL_ENV)/bin/pip
PYTHON ?= python3

# If it doesn't already exist
$(VIRTUAL_ENV):
	- mkdir -p $(INSTALL_DIR)
	$(info Creating virtual environment ...)
	$(PYTHON) -m venv $(VIRTUAL_ENV)
	$(info Updating pip)
	$(VIRTUAL_ENV)/bin/pip install --upgrade pip -q

pitot3 pitot3-condition-builder: $(VIRTUAL_ENV) gdtk-py
	$(info Installing python tools)
	$(PIP) install -e . --find-links $(INSTALL_DIR)/share/python-packages/
	- mkdir -p $(INSTALL_DIR)/bin
	- ln -sf $(VIRTUAL_ENV)/bin/$@ $(INSTALL_DIR)/bin/$@

gdtk-py:
	$(MAKE) -C $(GDTKPY_DIR) install

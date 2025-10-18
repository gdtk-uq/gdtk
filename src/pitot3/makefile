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

PYTHON ?= python3
PIP ?= $(PYTHON) -m pip

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


PY_VENV := $(INSTALL_DIR)/.venv/pitot3

.PHONY: .venv
.venv: gdtk-py
	- mkdir -p $(INSTALL_DIR)
	$(info Creating virtual environment ...)
	$(PYTHON) -m venv $(PY_VENV)
	$(info Updating pip)
	$(PY_VENV)/bin/pip install --upgrade pip -q
	$(PY_VENV)/bin/pip config set --site global.find-links $(INSTALL_DIR)/share/python-packages/
	$(info Installing python tools)
	$(PY_VENV)/bin/pip install -e . -q

pitot3 pitot3-condition-builder: .venv
	- mkdir -p $(INSTALL_DIR)/bin
	- ln -sf $(PY_VENV)/bin/$@ $(INSTALL_DIR)/bin/$@

gdtk-py:
	$(MAKE) -C $(GDTKPY_DIR) install

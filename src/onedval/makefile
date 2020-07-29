# makefile for onedval
#
# Do the bare minimum to get the python files into the installation area.
# We assume that the Python3 loadable library for Eilmer is already installed.
#
# PJ 2020-07-30

INSTALL_DIR ?= $(HOME)/dgdinst

install: onedval.py prop_avg.py cell.py
	cp onedval.py $(INSTALL_DIR)/bin/onedval
	chmod +x $(INSTALL_DIR)/bin/onedval
	cp prop_avg.py $(INSTALL_DIR)/lib/
	cp cell.py $(INSTALL_DIR)/lib/

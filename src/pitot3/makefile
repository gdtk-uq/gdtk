# makefile for pitot3
# this is based on the pitot one, which was based on the nenzfr one
# (so props to whoever did that)
#
# Chris James (c.james4@uq.edu.au) - 07-Jun-2021

INSTALL_DIR ?= $(DGD)

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk
LIBGASF := $(GAS_DIR)/libgasf.a

GASDYN_DIR := ../gasdyn
include $(GASDYN_DIR)/gasdyn_files.mk
 
install : $(INSTALL_DIR)
	@echo installing gas library
	cd $(GAS_DIR); make BUILD_DIR=$(BUILD_DIR) DMD=$(DMD) libgasf.a
	@echo Copying pitot3.py to $(INSTALL_DIR)/bin
	cp -f pitot3.py $(INSTALL_DIR)/bin
	@echo Copying pitot3_utils folder to $(INSTALL_DIR)/lib
	cp -r pitot3_utils $(INSTALL_DIR)/lib
	@echo Copying pitot3_data folder to $(INSTALL_DIR)/share
	cp -r pitot3_data $(INSTALL_DIR)/share
	@echo Making pitot3.py program executable
	chmod +x $(INSTALL_DIR)/bin/pitot3.py
	@echo -----------------------------------------
	@echo pitot3 should now be installed.
	@echo -----------------------------------------
	
clean : 
	cd $(GAS_DIR); make clean; rm -f libgas.a
	
$(INSTALL_DIR):
	- mkdir $(INSTALL_DIR)

# makefile for pitot3
# This is a newer one which I did with some help from Peter and Rowan.
# Now installs the gas library, CEA, and PITOT3.
# Chris James (c.james4@uq.edu.au) - 12-Jun-2021

INSTALL_DIR ?= $(HOME)/gdtkinst

GAS_DIR := ../gas
CEA_DIR := ../../extern/cea2

DMD ?= ldc2

install : $(INSTALL_DIR)
	@echo installing gas library
	cd $(GAS_DIR); make DMD=$(DMD) install
	@echo installing cea
	cd $(CEA_DIR); make install
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
	cd $(GAS_DIR); make clean
	cd $(CEA_DIR); make clean

$(INSTALL_DIR):
	- mkdir $(INSTALL_DIR)

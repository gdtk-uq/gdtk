INSTALL_DIR ?= $(HOME)/gdtkinst

PYTHON ?= python3
PIP ?= $(PYTHON) -m pip

GAS_DIR = ../gas
LIBGAS = $(GAS_DIR)/libgas.so

install: $(LIBGAS)
	- mkdir -p $(INSTALL_DIR)/share/python-packages
	cp $(LIBGAS) ./gdtk/libgas.so
	$(PIP) wheel --no-deps -e . -w $(INSTALL_DIR)/share/python-packages

clean:
	rm -rf build *.egg-info
	find . -type d -name '__pycache__' -exec rm -r {} +
	$(MAKE) -sC $(GAS_DIR) clean

$(LIBGAS):
	$(MAKE) -C $(GAS_DIR) libgas.so

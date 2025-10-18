INSTALL_DIR ?= $(HOME)/gdtkinst

PYTHON ?= python3
PIP ?= $(PYTHON) -m pip

GAS_DIR = ../gas
LIBGAS = $(GAS_DIR)/libgas.so

SLF_DIR = ../slf
LIBSLF = $(SLF_DIR)/libslf.so

install: $(LIBGAS) $(LIBSLF)
	- mkdir -p $(INSTALL_DIR)/share/python-packages
	cp $(LIBGAS) ./gdtk/libgas.so
	cp $(LIBSLF) ./gdtk/libslf.so
	$(PIP) wheel --no-deps -e . -w $(INSTALL_DIR)/share/python-packages

clean:
	rm -rf build *.egg-info
	rm -rf ./gdtk/*.so
	find . -type d -name '__pycache__' -exec rm -r {} +
	$(MAKE) -sC $(GAS_DIR) clean
	$(MAKE) -sC $(SLF_DIR) clean

$(LIBGAS):
	$(MAKE) -C $(GAS_DIR) libgas.so

$(LIBSLF):
	$(MAKE) -C $(SLF_DIR) libslf.so

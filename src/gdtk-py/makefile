INSTALL_DIR ?= $(HOME)/gdtkinst

PYTHON ?= python3
PIP ?= $(PYTHON) -m pip

GAS_DIR = ../gas
LIBGAS = $(GAS_DIR)/libgas.so

SLF_DIR = ../slf
LIBSLF = $(SLF_DIR)/libslf.so

wheel: build
	- mkdir -p $(INSTALL_DIR)/share/python-packages
	$(PIP) wheel --no-deps . -w $(INSTALL_DIR)/share/python-packages

build: $(LIBGAS) $(LIBSLF)
	cp $(LIBGAS) ./gdtk/libgas.so
	cp $(LIBSLF) ./gdtk/libslf.so

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

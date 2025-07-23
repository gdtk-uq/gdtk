CC = gcc
CFLAGS := -c -I. -fPIC -Wall -std=c99
CATALYST_INSTALL_DIR ?= ${HOME}/programs/catalyst

all:
	$(CC) $(CFLAGS) CatalystAdaptor.c -I $(CATALYST_INSTALL_DIR)/include/catalyst-2.0

clean: 
	rm -rf *.o

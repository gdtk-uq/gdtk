CC = gcc
CFLAGS := -I. -fPIC -Wall -std=c99 -O3
INSTALL_DIR ?= $(HOME)/eqc

all: thermo.o linalg.o common.o pt.o rhou.o ps.o rhot.o eqc.o
	$(CC) $(CFLAGS) -shared thermo.o linalg.o common.o pt.o rhou.o ps.o rhot.o eqc.o -lm -o libeqc.so
	ar r libeqc.a thermo.o linalg.o common.o pt.o rhou.o ps.o rhot.o eqc.o
	ranlib libeqc.a

thermo.o: thermo.c thermo.h
	$(CC) $(CFLAGS) -c thermo.c

linalg.o: linalg.c linalg.h
	$(CC) $(CFLAGS) -c linalg.c

common.o: common.c common.h
	$(CC) $(CFLAGS) -c common.c

pt.o: pt.c thermo.h linalg.h common.h pt.h
	$(CC) $(CFLAGS) -c pt.c

rhou.o: rhou.c thermo.h linalg.h common.h rhou.h
	$(CC) $(CFLAGS) -c rhou.c

ps.o: ps.c thermo.h linalg.h common.h ps.h
	$(CC) $(CFLAGS) -c ps.c

rhot.o: rhot.c thermo.h linalg.h common.h rhot.h
	$(CC) $(CFLAGS) -c rhot.c

eqc.o: eqc.c thermo.h pt.h rhou.h ps.h rhot.h eqc.h
	$(CC) $(CFLAGS) -c eqc.c

install: 
	mkdir -p $(INSTALL_DIR)
	cp libeqc.so $(INSTALL_DIR)
	sed -e "s+DBPATH='../thermo.inp'+DBPATH='$(INSTALL_DIR)/thermo.inp'+" \
        -e "s+LIBPATH='./libeqc.so'+LIBPATH='$(INSTALL_DIR)/libeqc.so'+" \
	    -e "s+HEADERFILE='./eqc.h'+HEADERFILE='$(INSTALL_DIR)/eqc.h'+" \
        pyeq.py > $(INSTALL_DIR)/pyeq.py
	cp 'clib.py' $(INSTALL_DIR)
	cp '../thermo.inp' $(INSTALL_DIR)
	cp 'eqc.h' $(INSTALL_DIR)
	cp -r ../tests $(INSTALL_DIR)/

moduletests: thermo.c thermo.h thermo.o linalg.c linalg.h linalg.o
	$(CC) thermo.c -D TEST -lm -o ../tests/testthermo
	$(CC) linalg.c -D TEST -lm -o ../tests/testlinalg

clean: 
	rm -rf *.o *.obj
	rm -rf libeqc.so libeqc.a

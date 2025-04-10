CEA_BUILD_DIR ?= ./build
CC = gcc
CFLAGS := -I. -fPIC -Wall -std=c99 -O3
INSTALL_DIR ?= $(HOME)/ceq
OBJS = $(CEA_BUILD_DIR)/thermo.o $(CEA_BUILD_DIR)/linalg.o $(CEA_BUILD_DIR)/common.o \
       $(CEA_BUILD_DIR)/pt.o $(CEA_BUILD_DIR)/rhou.o $(CEA_BUILD_DIR)/ps.o \
       $(CEA_BUILD_DIR)/rhot.o $(CEA_BUILD_DIR)/ceq.o

all: $(CEA_BUILD_DIR)/libceq.so $(CEA_BUILD_DIR)/libceq.a

$(CEA_BUILD_DIR)/libceq.so: $(OBJS)
	$(CC) $(CFLAGS) -shared $^ -lm -o $@

$(CEA_BUILD_DIR)/libceq.a: $(OBJS)
	ar r $@ $^
	ranlib $@

$(CEA_BUILD_DIR)/thermo.o: thermo.c thermo.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(CEA_BUILD_DIR)/linalg.o: linalg.c linalg.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(CEA_BUILD_DIR)/common.o: common.c common.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(CEA_BUILD_DIR)/pt.o: pt.c thermo.h linalg.h common.h pt.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(CEA_BUILD_DIR)/rhou.o: rhou.c thermo.h linalg.h common.h rhou.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(CEA_BUILD_DIR)/ps.o: ps.c thermo.h linalg.h common.h ps.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(CEA_BUILD_DIR)/rhot.o: rhot.c thermo.h linalg.h common.h rhot.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(CEA_BUILD_DIR)/ceq.o: ceq.c thermo.h pt.h rhou.h ps.h rhot.h ceq.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

install: 
	mkdir -p $(INSTALL_DIR)
	cp $(CEA_BUILD_DIR)/libceq.so $(INSTALL_DIR)
	sed -e "s+DBPATH='../thermo.inp'+DBPATH='$(INSTALL_DIR)/thermo.inp'+" \
        -e "s+LIBPATH='./libceq.so'+LIBPATH='$(INSTALL_DIR)/libceq.so'+" \
	    -e "s+HEADERFILE='./ceq.h'+HEADERFILE='$(INSTALL_DIR)/ceq.h'+" \
        pyeq.py > $(INSTALL_DIR)/pyeq.py
	cp clib.py $(INSTALL_DIR)
	cp ../thermo.inp $(INSTALL_DIR)
	cp ceq.h $(INSTALL_DIR)
	cp -r ../tests $(INSTALL_DIR)/

moduletests: $(CEA_BUILD_DIR)/thermo.o $(CEA_BUILD_DIR)/linalg.o
	$(CC) thermo.c -D TEST -lm -o ../tests/testthermo
	$(CC) linalg.c -D TEST -lm -o ../tests/testlinalg

clean:
	rm -rf $(CEA_BUILD_DIR)

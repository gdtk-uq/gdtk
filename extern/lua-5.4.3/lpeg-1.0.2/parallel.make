LIBNAME = lpeg
LUADIR = ../src/

# Build directory for out-of-tree builds (default "build")
LP_BUILD_DIR ?= build

COPT = -O2 -DNDEBUG
# COPT = -g

CWARNS = -Wall -Wextra -pedantic \
	-Waggregate-return \
	-Wcast-align \
	-Wcast-qual \
	-Wdisabled-optimization \
	-Wpointer-arith \
	-Wshadow \
	-Wsign-compare \
	-Wundef \
	-Wwrite-strings \
	-Wbad-function-cast \
	-Wdeclaration-after-statement \
	-Wmissing-prototypes \
	-Wnested-externs \
	-Wstrict-prototypes
# -Wunreachable-code \

CFLAGS = $(CWARNS) $(COPT) -std=c99 -I$(LUADIR) -fPIC
CC = gcc

# List of object files, now in LP_BUILD_DIR
OBJS = $(LP_BUILD_DIR)/lpvm.o $(LP_BUILD_DIR)/lpcap.o $(LP_BUILD_DIR)/lptree.o $(LP_BUILD_DIR)/lpcode.o $(LP_BUILD_DIR)/lpprint.o

# For Linux: builds shared library with Linux-specific flags.
linux:
	$(MAKE) -f parallel.make $(LP_BUILD_DIR)/lpeg.so "DLLFLAGS = -shared -fPIC"

# For Mac OS: builds shared library with Mac-specific flags.
macosx:
	$(MAKE) -f parallel.make $(LP_BUILD_DIR)/lpeg.so "DLLFLAGS = -bundle -undefined dynamic_lookup"

# Build shared library in LP_BUILD_DIR
$(LP_BUILD_DIR)/lpeg.so: $(OBJS)
	env $(CC) $(DLLFLAGS) $(OBJS) -o $(LP_BUILD_DIR)/lpeg.so

# Pattern rule to compile .c files into .o files in LP_BUILD_DIR.
$(LP_BUILD_DIR)/%.o: %.c
	@mkdir -p $(LP_BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Force objects to rebuild if the makefile changes.
$(OBJS): makefile

# If your test expects lpeg.so in the current directory, adjust this target as needed.
test: test.lua re.lua $(LP_BUILD_DIR)/lpeg.so
	./test.lua

clean:
	rm -rf $(LP_BUILD_DIR)

# Dependency rules updated for the build directory.
$(LP_BUILD_DIR)/lpcap.o: lpcap.c lpcap.h lptypes.h
$(LP_BUILD_DIR)/lpcode.o: lpcode.c lptypes.h lpcode.h lptree.h lpvm.h lpcap.h
$(LP_BUILD_DIR)/lpprint.o: lpprint.c lptypes.h lpprint.h lptree.h lpvm.h lpcap.h
$(LP_BUILD_DIR)/lptree.o: lptree.c lptypes.h lpcap.h lpcode.h lptree.h lpvm.h lpprint.h
$(LP_BUILD_DIR)/lpvm.o: lpvm.c lpcap.h lptypes.h lpvm.h lpprint.h lptree.h

# Modified Makefile for building Lua into a specified build directory.
# Build directory is specified by the make variable LUA_BUILD_DIR (default is "build")
LUA_BUILD_DIR ?= build

# == CHANGE THE SETTINGS BELOW TO SUIT YOUR ENVIRONMENT =======================

# Your platform. See PLATS for possible values.
PLAT= guess

CC= gcc -std=gnu99
CFLAGS= -O2 -Wall -Wextra -DLUA_COMPAT_5_3 $(SYSCFLAGS) $(MYCFLAGS)
LDFLAGS= $(SYSLDFLAGS) $(MYLDFLAGS)
LIBS= -lm $(SYSLIBS) $(MYLIBS)

AR= ar rcu
RANLIB= ranlib
RM= rm -f
UNAME= uname

SYSCFLAGS=
SYSLDFLAGS=
SYSLIBS=

MYCFLAGS=
MYLDFLAGS=
MYLIBS=
MYOBJS=

MAKEF = $(MAKE) -f parallel.make

# Special flags for compiler modules; -Os reduces code size.
CMCFLAGS=

# == END OF USER SETTINGS -- NO NEED TO CHANGE ANYTHING BELOW THIS LINE =======

PLATS= guess aix bsd c89 freebsd generic linux linux-readline macosx mingw posix solaris

# All files are built in the build directory (OBJ_DIR)
OBJ_DIR= $(LUA_BUILD_DIR)

# Targets (all files are prefixed with OBJ_DIR)
LUA_A= $(OBJ_DIR)/liblua.a
CORE_O= $(OBJ_DIR)/lapi.o $(OBJ_DIR)/lcode.o $(OBJ_DIR)/lctype.o $(OBJ_DIR)/ldebug.o $(OBJ_DIR)/ldo.o $(OBJ_DIR)/ldump.o $(OBJ_DIR)/lfunc.o $(OBJ_DIR)/lgc.o $(OBJ_DIR)/llex.o $(OBJ_DIR)/lmem.o $(OBJ_DIR)/lobject.o $(OBJ_DIR)/lopcodes.o $(OBJ_DIR)/lparser.o $(OBJ_DIR)/lstate.o $(OBJ_DIR)/lstring.o $(OBJ_DIR)/ltable.o $(OBJ_DIR)/ltm.o $(OBJ_DIR)/lundump.o $(OBJ_DIR)/lvm.o $(OBJ_DIR)/lzio.o

LIB_O= $(OBJ_DIR)/lauxlib.o $(OBJ_DIR)/lbaselib.o $(OBJ_DIR)/lcorolib.o $(OBJ_DIR)/ldblib.o $(OBJ_DIR)/liolib.o $(OBJ_DIR)/lmathlib.o $(OBJ_DIR)/loadlib.o $(OBJ_DIR)/loslib.o $(OBJ_DIR)/lstrlib.o $(OBJ_DIR)/ltablib.o $(OBJ_DIR)/lutf8lib.o $(OBJ_DIR)/linit.o

BASE_O= $(CORE_O) $(LIB_O) $(MYOBJS)

# RJG, 2021-11-15: Change executable name to clash with system Lua
LUA_T= $(OBJ_DIR)/dgd-lua
LUA_O= $(OBJ_DIR)/lua.o

# RJG, 2021-11-15: Change executable name to clash with system Lua
LUAC_T= $(OBJ_DIR)/dgd-luac
LUAC_O= $(OBJ_DIR)/luac.o

ALL_O= $(BASE_O) $(LUA_O) $(LUAC_O)
ALL_T= $(LUA_A) $(LUA_T) $(LUAC_T)
ALL_A= $(LUA_A)

# Pattern rule: compile any .c file into an object in the build directory.
$(OBJ_DIR)/%.o: %.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(CMCFLAGS) -c $< -o $@

# Targets start here.
default: $(PLAT)

all:	$(ALL_T)

o:	$(ALL_O)

a:	$(ALL_A)

$(LUA_A): $(BASE_O)
	$(AR) $@ $(BASE_O)
	$(RANLIB) $@

$(LUA_T): $(LUA_O) $(LUA_A)
	$(CC) -o $@ $(LDFLAGS) $(LUA_O) $(LUA_A) $(LIBS)

$(LUAC_T): $(LUAC_O) $(LUA_A)
	$(CC) -o $@ $(LDFLAGS) $(LUAC_O) $(LUA_A) $(LIBS)

test: $(LUA_T)
	./$(LUA_T) -v

clean:
	$(RM) $(ALL_T) $(ALL_O)

depend:
	@$(CC) $(CFLAGS) -MM l*.c

echo:
	@echo "PLAT= $(PLAT)"
	@echo "CC= $(CC)"
	@echo "CFLAGS= $(CFLAGS)"
	@echo "LDFLAGS= $(SYSLDFLAGS)"
	@echo "LIBS= $(LIBS)"
	@echo "AR= $(AR)"
	@echo "RANLIB= $(RANLIB)"
	@echo "RM= $(RM)"
	@echo "UNAME= $(UNAME)"
	@echo "LUA_BUILD_DIR= $(LUA_BUILD_DIR)"

# Convenience targets for popular platforms.
ALL= all

help:
	@echo "Do 'make PLATFORM' where PLATFORM is one of these:"
	@echo "   $(PLATS)"
	@echo "See doc/readme.html for complete instructions."

guess:
	@echo Guessing `$(UNAME)`
	@$(MAKEF) LUA_BUILD_DIR=$(LUA_BUILD_DIR) `$(UNAME)`

AIX aix:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) CC="xlc" CFLAGS="-O2 -DLUA_USE_POSIX -DLUA_USE_DLOPEN" SYSLIBS="-ldl" SYSLDFLAGS="-brtl -bexpall"

bsd:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) SYSCFLAGS="-DLUA_USE_POSIX -DLUA_USE_DLOPEN" SYSLIBS="-Wl,-E"

c89:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) SYSCFLAGS="-DLUA_USE_C89" CC="gcc -std=c89"
	@echo ''
	@echo '*** C89 does not guarantee 64-bit integers for Lua.'
	@echo '*** Make sure to compile all external Lua libraries'
	@echo '*** with LUA_USE_C89 to ensure consistency'
	@echo ''

FreeBSD NetBSD OpenBSD freebsd:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) SYSCFLAGS="-DLUA_USE_LINUX -DLUA_USE_READLINE -I/usr/include/edit" SYSLIBS="-Wl,-E -ledit" CC="cc"

generic: $(ALL)

Linux linux:	linux-noreadline

linux-noreadline:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) SYSCFLAGS="-DLUA_USE_LINUX -fPIC -fno-omit-frame-pointer" SYSLIBS="-Wl,-E -ldl"

linux-readline:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) SYSCFLAGS="-DLUA_USE_LINUX -DLUA_USE_READLINE" SYSLIBS="-Wl,-E -ldl -lreadline"

Darwin macos macosx:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) SYSCFLAGS="-DLUA_USE_MACOSX -DLUA_USE_READLINE" SYSLIBS="-lreadline"

mingw:
	$(MAKEF) "LUA_A=lua54.dll" "LUA_T=lua.exe" "LUA_BUILD_DIR=$(LUA_BUILD_DIR)" AR="$(CC) -shared -o" RANLIB="strip --strip-unneeded" SYSCFLAGS="-DLUA_BUILD_AS_DLL" SYSLIBS="" SYSLDFLAGS="-s" lua.exe
	$(MAKEF) "LUAC_T=luac.exe" "LUA_BUILD_DIR=$(LUA_BUILD_DIR)" luac.exe

posix:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) SYSCFLAGS="-DLUA_USE_POSIX"

SunOS solaris:
	$(MAKEF) $(ALL) LUA_BUILD_DIR=$(LUA_BUILD_DIR) SYSCFLAGS="-DLUA_USE_POSIX -DLUA_USE_DLOPEN -D_REENTRANT" SYSLIBS="-ldl"

# Compiler modules may use special flags.
$(OBJ_DIR)/llex.o: llex.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(CMCFLAGS) -c llex.c -o $(OBJ_DIR)/llex.o

$(OBJ_DIR)/lparser.o: lparser.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(CMCFLAGS) -c lparser.c -o $(OBJ_DIR)/lparser.o

$(OBJ_DIR)/lcode.o: lcode.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(CMCFLAGS) -c lcode.c -o $(OBJ_DIR)/lcode.o

# DO NOT DELETE

$(OBJ_DIR)/lapi.o: lapi.c lprefix.h lua.h luaconf.h lapi.h llimits.h lstate.h lobject.h ltm.h lzio.h lmem.h ldebug.h ldo.h lfunc.h lgc.h lstring.h ltable.h lundump.h lvm.h
$(OBJ_DIR)/lauxlib.o: lauxlib.c lprefix.h lua.h luaconf.h lauxlib.h
$(OBJ_DIR)/lbaselib.o: lbaselib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/lcode.o: lcode.c lprefix.h lua.h luaconf.h lcode.h llex.h lobject.h llimits.h lzio.h lmem.h lopcodes.h lparser.h ldebug.h lstate.h ltm.h ldo.h lgc.h lstring.h ltable.h lvm.h
$(OBJ_DIR)/lcorolib.o: lcorolib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/lctype.o: lctype.c lprefix.h lctype.h lua.h luaconf.h llimits.h
$(OBJ_DIR)/ldblib.o: ldblib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/ldebug.o: ldebug.c lprefix.h lua.h luaconf.h lapi.h llimits.h lstate.h lobject.h ltm.h lzio.h lmem.h lcode.h llex.h lopcodes.h lparser.h ldebug.h ldo.h lfunc.h lstring.h lgc.h ltable.h lvm.h
$(OBJ_DIR)/ldo.o: ldo.c lprefix.h lua.h luaconf.h lapi.h llimits.h lstate.h lobject.h ltm.h lzio.h lmem.h ldebug.h ldo.h lfunc.h lgc.h lopcodes.h lparser.h lstring.h ltable.h lundump.h lvm.h
$(OBJ_DIR)/ldump.o: ldump.c lprefix.h lua.h luaconf.h lobject.h llimits.h lstate.h ltm.h lzio.h lmem.h lundump.h
$(OBJ_DIR)/lfunc.o: lfunc.c lprefix.h lua.h luaconf.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h ldo.h lgc.h
$(OBJ_DIR)/lgc.o: lgc.c lprefix.h lua.h luaconf.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h ldo.h lfunc.h lgc.h lstring.h ltable.h
$(OBJ_DIR)/linit.o: linit.c lprefix.h lua.h luaconf.h lualib.h lauxlib.h
$(OBJ_DIR)/liolib.o: liolib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/llex.o: llex.c lprefix.h lua.h luaconf.h lctype.h llimits.h ldebug.h lstate.h lobject.h ltm.h lzio.h lmem.h ldo.h lgc.h llex.h lparser.h lstring.h ltable.h
$(OBJ_DIR)/lmathlib.o: lmathlib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/lmem.o: lmem.c lprefix.h lua.h luaconf.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h ldo.h lgc.h
$(OBJ_DIR)/loadlib.o: loadlib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/lobject.o: lobject.c lprefix.h lua.h luaconf.h lctype.h llimits.h ldebug.h lstate.h lobject.h ltm.h lzio.h lmem.h ldo.h lstring.h lgc.h lvm.h
$(OBJ_DIR)/lopcodes.o: lopcodes.c lprefix.h lopcodes.h llimits.h lua.h luaconf.h
$(OBJ_DIR)/loslib.o: loslib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/lparser.o: lparser.c lprefix.h lua.h luaconf.h lcode.h llex.h lobject.h llimits.h lzio.h lmem.h lopcodes.h lparser.h ldebug.h lstate.h ltm.h ldo.h lfunc.h lstring.h lgc.h ltable.h
$(OBJ_DIR)/lstate.o: lstate.c lprefix.h lua.h luaconf.h lapi.h llimits.h lstate.h lobject.h ltm.h lzio.h lmem.h ldebug.h ldo.h lfunc.h lgc.h llex.h lstring.h ltable.h
$(OBJ_DIR)/lstring.o: lstring.c lprefix.h lua.h luaconf.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h ldo.h lstring.h lgc.h
$(OBJ_DIR)/lstrlib.o: lstrlib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/ltable.o: ltable.c lprefix.h lua.h luaconf.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h ldo.h lgc.h lstring.h ltable.h lvm.h
$(OBJ_DIR)/ltablib.o: ltablib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/ltm.o: ltm.c lprefix.h lua.h luaconf.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h ldo.h lgc.h lstring.h ltable.h lvm.h
$(OBJ_DIR)/lua.o: lua.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/luac.o: luac.c lprefix.h lua.h luaconf.h lauxlib.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h lopcodes.h lopnames.h lundump.h
$(OBJ_DIR)/lundump.o: lundump.c lprefix.h lua.h luaconf.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h ldo.h lfunc.h lstring.h lgc.h lundump.h
$(OBJ_DIR)/lutf8lib.o: lutf8lib.c lprefix.h lua.h luaconf.h lauxlib.h lualib.h
$(OBJ_DIR)/lvm.o: lvm.c lprefix.h lua.h luaconf.h ldebug.h lstate.h lobject.h llimits.h ltm.h lzio.h lmem.h ldo.h lfunc.h lgc.h lopcodes.h lstring.h ltable.h lvm.h ljumptab.h
$(OBJ_DIR)/lzio.o: lzio.c lprefix.h lua.h luaconf.h llimits.h lmem.h lstate.h lobject.h ltm.h lzio.h

# (end of Makefile)

PLATFORM ?= linux
LUA := ../../extern/lua-5.4.3
LIBLUA := $(LUA)/install/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(LUAD_DIR)/luad/*.d \
	$(LUAD_DIR)/luad/c/*.d \
	$(LUAD_DIR)/luad/conversions/*.d

$(LIBLUA):
	cd $(LUA); make $(PLATFORM) local

lua_service : lua_service.d $(LIBLUA)
	dmd -main -unittest -g lua_service.d ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

clean:
	- rm -f *.o *.obj lua_service
	- cd $(LUA); make clean

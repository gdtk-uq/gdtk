LUA := ../../extern/lua-5.1.4
LIBLUA := $(LUA)/lib/liblua.a
LIBLUAPATH := $(LUA)/lib
LUAD_DIR := ../extern/LuaD
LUAD_FILES := $(LUAD_DIR)/luad/*.d \
	$(LUAD_DIR)/luad/c/*.d \
	$(LUAD_DIR)/luad/conversions/*.d

$(LIBLUA):
	cd $(LUA); make linux local

lua_service : lua_service.d $(LIBLUA)
	dmd -main -unittest -gc lua_service.d ${LUAD_FILES} -L-L$(LIBLUAPATH) -L-llua -L-ldl

clean:
	- rm *.o lua_service
	- cd $(LUA); make clean

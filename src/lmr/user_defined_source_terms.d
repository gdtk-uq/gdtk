/**
 * user_defined_source_terms.d
 *
 * This module handles user-defined source terms
 * that the user might specify with a Lua script.
 *
 * Authors: RG & PJ
 * Date: 2015-03-17
 */

module user_defined_source_terms;

import std.conv;
import std.stdio;
import std.string;
import util.lua;
import util.lua_service;
import lua_helper;
import gas;
import lmr.fluidfvcell;
import globalconfig;

void getUDFSourceTermsForCell(lua_State* L, FluidFVCell cell, size_t gtl,
                              double t, LocalConfig myConfig,
                              size_t blkId, size_t i, size_t j, size_t k)
{
    auto gmodel = myConfig.gmodel;
    auto cqi = myConfig.cqi;
    size_t n_species = gmodel.n_species;
    size_t n_modes = gmodel.n_modes;
    //
    // Push user function onto TOS
    lua_getglobal(L, "sourceTerms");
    // Push sim_time onto TOS
    lua_pushnumber(L, t);
    // Push cell data into an args table and onto TOS
    lua_newtable(L);
    int tblIdx = lua_gettop(L);
    pushFluidCellToTable(L, tblIdx, cell, gtl, myConfig);
    lua_pushinteger(L, blkId); lua_setfield(L, tblIdx, "blkId");
    lua_pushinteger(L, i); lua_setfield(L, tblIdx, "i");
    lua_pushinteger(L, j); lua_setfield(L, tblIdx, "j");
    lua_pushinteger(L, k); lua_setfield(L, tblIdx, "k");
    // Call sourceTerms function with (t, args)
    int number_args = 2;
    int number_results = 1;
    //
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            luaL_error(L, "error running user soure terms function: %s\n",
                       lua_tostring(L, -1));
    }
    //
    // Grab values from user-returned table at TOS
    // For any missing values, put in 0.0
    if (cqi.mass==0) cell.Qudf[cqi.mass] = getNumberFromTable(L, -1, "mass", false, 0.0);
    cell.Qudf[cqi.xMom] = getNumberFromTable(L, -1, "momentum_x", false, 0.0);
    cell.Qudf[cqi.yMom] = getNumberFromTable(L, -1, "momentum_y", false, 0.0);
    if (cqi.threeD) { cell.Qudf[cqi.zMom] = getNumberFromTable(L, -1, "momentum_z", false, 0.0); }
    cell.Qudf[cqi.totEnergy] = getNumberFromTable(L, -1, "total_energy",false, 0.0);
    version(MHD) {
        if (cqi.MHD) {
            cell.Qudf[cqi.xB] = getNumberFromTable(L, -1, "xB", false, 0.0);
            cell.Qudf[cqi.yB] = getNumberFromTable(L, -1, "yB", false, 0.0);
            cell.Qudf[cqi.zB] = getNumberFromTable(L, -1, "zB", false, 0.0);
            cell.Qudf[cqi.psi] = getNumberFromTable(L, -1, "psi", false, 0.0);
            cell.Qudf[cqi.divB] = getNumberFromTable(L, -1, "divB", false, 0.0);
        }
    }
    version(turbulence) {
        foreach(it; 0 .. myConfig.turb_model.nturb){
            string tname = myConfig.turb_model.primitive_variable_name(it);
            cell.Qudf[cqi.rhoturb+it] = getNumberFromTable(L, -1, tname, false, 0.0);
        }
    }
    version(multi_species_gas) {
        if (cqi.n_species > 1) {
            lua_getfield(L, -1, "species");
            if ( lua_istable(L, -1) ) {
                // Iterate over species by names.
                int idx = lua_gettop(L);
                lua_pushnil(L);
                while ( lua_next(L, idx) != 0 ) {
                    string key = to!string(lua_tostring(L, -2));
                    auto isp = gmodel.species_index(key);
                    if ( isp == -1 ) {
                        string errMsg = format("ERROR: In the user-defined source terms, the species name '%s'\n", key);
                        errMsg ~= "in the species table is not a valid species name.\n";
                        lua_pop(L, 1);
                        throw new LuaInputException(errMsg);
                    }
                    cell.Qudf[cqi.species+isp] = lua_tonumber(L, -1);
                    lua_pop(L, 1);
                }
                lua_pop(L, 1); // discard species table
            } else {
                lua_pop(L, 1); // discard species item first
                // For the single-species case, we just set the
                // source terms of the single-species to equal
                // that of the mass term.
                // if (n_species == 1) { There is no storage for cell.Qudf[cqi.species+0] }
                // For multi-component gases, there is really no sensible decision,
                // so leave it alone.
            }
        }
    }
    version(multi_T_gas) {
        if (n_modes > 0) {
            lua_getfield(L, -1, "energies");
            if ( !lua_isnil(L, -1) ) {
                foreach (imode; 0 .. n_modes) {
                    lua_rawgeti(L, -1, to!int(imode)+1);
                    cell.Qudf[cqi.modes+imode] = lua_tonumber(L, -1);
                    lua_pop(L, 1);
                }
            }
            lua_pop(L, 1);
        }
    }
    // Clear stack.
    lua_settop(L, 0);
} // end getUDFSourceTermsForCell()

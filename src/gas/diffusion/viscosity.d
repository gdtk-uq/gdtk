/*
 * viscosity.d
 * Interface for all viscosity models.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 */

module gas.diffusion.viscosity;

import std.string;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;
import gas.gas_model;
import gas.gas_state;
import gas.diffusion.sutherland_viscosity;
import gas.diffusion.cea_viscosity;
import gas.diffusion.chemkin_viscosity;

interface Viscosity {
    Viscosity dup() const;
    @nogc number eval(number T);
    @nogc number eval(number T, number logT);
}

interface ViscosityMixtureModel {
    ViscosityMixtureModel dup() const;
    @nogc final void update_viscosity(ref GasState Q)
    {
        Q.mu = eval(Q);
    }
    @nogc number eval(in GasState Q);
}


/**
 * Create and return a new Viscosity model.
 *
 * It is assumed that the Lua table describing the
 * Viscosity model is sitting at the top-of-stack
 * of the passed-in lua_State. At the end of the
 * function, that table is left at the top-of-stack.
 */
Viscosity createViscosityModel(lua_State *L)
{
    Viscosity vm;
    auto model = getString(L, -1, "model");

    switch ( model ) {
    case "Sutherland":
        vm = createSutherlandViscosity(L);
        break;
    case "CEA":
        vm = createCEAViscosity(L);
        break;
    case "Chemkin":
        vm = createChemkinViscosity(L);
        break;
    default:
        string errMsg = format("The requested Viscosity model '%s' is not available.", model);
        throw new Error(errMsg);
    }
    return vm;
}

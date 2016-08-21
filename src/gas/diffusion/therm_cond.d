/*
 * therm_cond.d
 * Interface for all thermal conductivity models.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 */

module gas.diffusion.therm_cond;

import std.string;
import util.lua;
import util.lua_service;
import gas.gas_model;
import gas.diffusion.sutherland_therm_cond;
import gas.diffusion.cea_therm_cond;

interface ThermalConductivity {
    ThermalConductivity dup() const;
    final void update_thermal_conductivity(GasState Q) 
    {
	for ( auto imode = 0; imode < Q.T.length; ++imode) {
	    Q.k[imode] = eval(Q, imode);
	}
    }
    double eval(in GasState Q, int imode);
}

/**
 * Create and return a new ThermalConductivity model.
 *
 * It is assumed that the Lua table describing the 
 * ThermalConductivity model is sitting at the top-of-stack
 * of the passed-in lua_State. At the end of the 
 * function, that table is left at the top-of-stack.
 */
ThermalConductivity createThermalConductivityModel(lua_State *L)
{
    ThermalConductivity tcm;
    auto model = getString(L, -1, "model");
    
    switch ( model ) {
    case "Sutherland":
	tcm = createSutherlandThermalConductivity(L);
	break;
    case "CEA":
	tcm = createCEAThermalConductivity(L);
	break;
    default:
	string errMsg = format("The requested ThermalConductivity model '%s' is not available.", model);
	throw new Error(errMsg);
    }
    
    return tcm;
}

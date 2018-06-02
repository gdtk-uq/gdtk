/*
 * therm_cond.d
 * Interface for all thermal conductivity models.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 */

module gas.diffusion.therm_cond;

import std.string;
import nm.complex;
import nm.number;
import util.lua;
import util.lua_service;
import gas.gas_model;
import gas.gas_state;
import gas.diffusion.sutherland_therm_cond;
version(complex_numbers) {
    // do nothing here
} else {
    import gas.diffusion.cea_therm_cond;
}

interface ThermalConductivity {
    ThermalConductivity dup() const;
    final void update_thermal_conductivity(GasState Q) 
    {
        Q.k = eval(Q, Q.T);
        for ( auto imode = 0; imode < Q.T_modes.length; ++imode) {
            Q.k_modes[imode] = eval(Q, Q.T_modes[imode]);
        }
    }
    number eval(ref const(GasState) Q, number T);
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
    
    version(complex_numbers) {
        // Limited number of options.
        switch ( model ) {
        case "Sutherland":
            tcm = createSutherlandThermalConductivity(L);
            break;
        default:
            string errMsg = format("The requested ThermalConductivity model '%s' is not available.", model);
            throw new Error(errMsg);
        }
    } else {
        // All options for double_numbers.
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
    }
    return tcm;
}

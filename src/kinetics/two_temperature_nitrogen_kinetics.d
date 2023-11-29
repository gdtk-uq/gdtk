/**
 * Authors: Rowan G. and Peter J.
 * Date: 2017-07-13
 *
 * A model for for vibrational relaxation kinetics in
 * single-species diatomic nitrogen.
 */

module kinetics.two_temperature_nitrogen_kinetics;

import std.math;
import std.conv;
import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;
import kinetics.thermochemical_reactor;

alias RelaxTimeFunc = @nogc number function(number, number);

final class VibRelaxNitrogen : ThermochemicalReactor {
    this(string fname, GasModel gmodel)
    {
        super(gmodel);
        _Q_eq = GasState(gmodel);
        // Read parameters from file.
        auto L = init_lua_State();
        doLuaFile(L, fname);
        lua_getglobal(L, "relaxation_time");
        string relaxTime = to!string(luaL_checkstring(L, -1));
        lua_pop(L, 1);
        lua_close(L);

        final switch(relaxTime) {
        case "Blackman":
            _relaxTimeCalc = &BlackmanRelaxationTime;
            break;
        case "Millikan-White":
            _relaxTimeCalc = &MWRelaxationTime;
            break;
        }
    } // end constructor

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        number tau = _relaxTimeCalc(Q.T, Q.p);
        // Find the total internal energy in the gas
        number uTotal = Q.u + Q.u_modes[0];
        // Find the vib energy at equilibrium with T
        _Q_eq.T = Q.T;
        _Q_eq.T_modes[0] = Q.T;
        _Q_eq.p = Q.p;
        _gmodel.update_thermo_from_pT(_Q_eq);
        number u_v_eq = _Q_eq.u_modes[0];
        number u_v = Q.u_modes[0];
        Q.u_modes[0] = u_v_eq + (u_v - u_v_eq)*exp(-tInterval/tau);
        Q.u = uTotal - Q.u_modes[0];
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);
    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        string errMsg = "eval_source_terms not implemented for two_temperature_nitrogen_kinetics.";
        throw new ThermochemicalReactorUpdateException(errMsg);
    }

private:
    GasState _Q_eq;
    RelaxTimeFunc _relaxTimeCalc;
}

@nogc
number BlackmanRelaxationTime(number T, number p)
{
    double A = 7.12e-9;
    double B = 124.07;
    number pAtm = p/P_atm;
    number tau = (A/pAtm)*exp(124.07/pow(T, 1./3.));
    return tau;
}

@nogc
number MWRelaxationTime(number T, number p)
{
    double a = 221.0;
    double b = 0.0290;
    number pAtm = p/P_atm;
    number pTau = exp(a*(pow(T, -1./3.) - b) - 18.42);
    number tau = pTau/pAtm;
    return tau;
}

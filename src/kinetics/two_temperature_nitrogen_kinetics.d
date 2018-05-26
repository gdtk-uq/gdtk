/**
 * Authors: Rowan G. and Peter J.
 * Date: 2017-07-13
 *
 * A model for for vibrational relaxation kinetics in
 * single-species diatomic nitrogen.
 */

module kinetics.two_temperature_nitrogen_kinetics;

import std.math : exp, pow;
import std.conv;

import util.lua;
import util.lua_service;
import gas;
import kinetics.thermochemical_reactor;

alias RelaxTimeFunc = double function(double, double);

final class VibRelaxNitrogen : ThermochemicalReactor {
    this(string fname, GasModel gmodel)
    {
        super(gmodel);
        _Q_eq = new GasState(gmodel);
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

    }

    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest, 
                         ref double[] params)
    {
        double tau = _relaxTimeCalc(Q.T, Q.p);
        // Find the total internal energy in the gas
        double uTotal = Q.u + Q.u_modes[0];
        // Find the vib energy at equilibrium with T
        _Q_eq.T = Q.T;
        _Q_eq.T_modes[0] = Q.T;
        _Q_eq.p = Q.p;
        _gmodel.update_thermo_from_pT(_Q_eq);
        double u_v_eq = _Q_eq.u_modes[0];
        double u_v = Q.u_modes[0];
        Q.u_modes[0] = u_v_eq + (u_v - u_v_eq)*exp(-tInterval/tau);
        Q.u = uTotal - Q.u_modes[0];
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);
    }

private:
    GasState _Q_eq;
    RelaxTimeFunc _relaxTimeCalc;
}

double BlackmanRelaxationTime(double T, double p)
{
    double A = 7.12e-9;
    double B = 124.07;
    double pAtm = p/P_atm;
    double tau = (A/pAtm)*exp(124.07/pow(T, 1./3.));
    return tau;
}

double MWRelaxationTime(double T, double p)
{
    double a = 221.0;
    double b = 0.0290;
    double pAtm = p/P_atm;
    double pTau = exp(a*(pow(T, -1./3.) - b) - 18.42);
    double tau = pTau/pAtm;
    return tau;
}

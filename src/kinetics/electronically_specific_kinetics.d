/**
 * Authors: Brad Semple
 * Date: 2018-09-04
 *
 * This module provides an integrator for the master
 * equation for an electronic state-specific system
 * 
 */

module kinetics.electronically_specific_kinetics;

import core.stdc.stdlib : exit;
import std.stdio;
import std.conv : to;
import std.string;
import std.math;
import std.algorithm;

import nm.complex;
import nm.number;
import nm.smla;
import util.lua;
import util.lua_service;
import gas;
import gas.physical_constants;
import gas.electronically_specific_gas;

import kinetics.thermochemical_reactor;
import kinetics.electronic_state_solver;

final class ElectronicallySpecificKinetics : ThermochemicalReactor {
public:

    this(GasModel gmodel)
    {
        super(gmodel);
        _numden.length = gmodel.n_species;
        _numden0.length = gmodel.n_species;

        kinetics.electronic_state_solver.Init();
    }

    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[] params)
    {
        _gmodel.massf2numden(Q, _numden0);
        //update the number density vector for time t+t_interval
        //by solving the kinetics ODE (electron temp frozen)
        
        Electronic_Solve(_numden0, _numden, Q.T_modes[0], tInterval);
        _gmodel.numden2massf(_numden, Q);
    }

private:
    number[] _numden;
    number[] _numden0; 
}

version(electronically_specific_kinetics_test) 
{
    int main() 
    {
        import util.msg_service;

        auto L = init_lua_State();
        doLuaFile(L, "../gas/sample-data/electronic_composition.lua");
        auto gm = new ElectronicallySpecificGas(L);
        auto gd = new GasState(19,19);

        gd.massf[] = 0.0;
        gd.massf[0] = 0.81403036047055; //initialises massf of NI
        gd.massf[9] = 0.185968105968037; //initialises massf of OI
        gd.massf[18] = 1.0 - (gd.massf[0] + gd.massf[9]); //tiny massf for free electron
        gd.p = 1.0e3;
        gd.T = 10000.0;
        gd.T_modes[0] = 10000.0;
        gm.update_thermo_from_pT(gd);
        

        lua_close(L);

        double _dt = 1e-10;
        double _duration = 1e-9;
        double _t = 0.0;

        auto L2 = init_lua_State();
        writeln("Parsing kinetics file.");
        writeln(gm.n_modes);
        doLuaFile(L2, "sample-input/state-specific-N2-diss.lua");

        ElectronicallySpecificKinetics psk = new ElectronicallySpecificKinetics(gm);

        writeln("Gas state BEFORE update: ", gd);

        double dtChemSuggest = -1.0;
        double dtThermSuggest = -1.0;
        number[] params;

        while (_t < _duration) {
            psk(gd, dt, dtChemSuggest, dtThermSuggest, params);
            _t+=dt;
        }
        writeln("Gas state AFTER update: ",gd);

        return 0;

    }
}
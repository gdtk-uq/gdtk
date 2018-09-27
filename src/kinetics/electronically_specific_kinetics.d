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
        _numden_input.length = gmodel.n_species - 2;
        _numden_output.length = gmodel.n_species - 2;

        kinetics.electronic_state_solver.Init();
    }

    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[maxParams] params)
    {
        foreach(int i; 0 .. _gmodel.n_species){ //give massf values for the required species
            _numden[i] = Q.massf[i];
        }
        
        _gmodel.massf2numden(Q, _numden); //Convert from mass fraction to number density

        foreach(int i;0 .. _gmodel.n_species - 2) {
            _numden_input[i] = _numden[i];
        }
        //update the number density vector for time t+t_interval
        //by solving the kinetics ODE (electron temp frozen)
        Electronic_Solve(_numden_input, _numden_output, Q.T_modes[0], tInterval);

        foreach (int i; 0 .. _gmodel.n_species - 2) {
            _numden[i] = _numden_output[i];
        } 

        _gmodel.numden2massf(_numden, Q); //Convert back to number density

        //Maybe include dissociation/recombination for N2 and O2 in a decoupled reaction here???

    }

private:
    number[] _numden; //Total set of species including static N2 and O2
    number[] _numden_input; //Input for the electronic state solver, exclusive of N2 and O2
    number[] _numden_output; //Output for the state solver
}

version(electronically_specific_kinetics_test) 
{
    int main() 
    {
        import util.msg_service;

        auto L = init_lua_State();
        doLuaFile(L, "../kinetics/sample-input/electronic_composition.lua");
        auto gm = new ElectronicallySpecificGas(L);
        auto gd = new GasState(21,1);

        gd.massf[] = 0.0;
        gd.massf[0] = 0.81403036047055; //initialises massf of NI
        gd.massf[9] = 0.185968105968037; //initialises massf of OI
        gd.massf[18] = 1.0 - (gd.massf[0] + gd.massf[9]); //tiny massf for free electron
        gd.p = 10.0;
        gd.T = 7000.0;
        gd.T_modes[0] = 10000.0;
        gm.update_thermo_from_pT(gd);
        

        lua_close(L);

        double _dt = 1e-10;
        double _duration = 1e-9;
        double _t = 0.0;

        auto L2 = init_lua_State();
        doLuaFile(L2, "sample-input/electronic_composition.lua");

        ElectronicallySpecificKinetics esk = new ElectronicallySpecificKinetics(gm);

        double dtChemSuggest = -1.0;
        double dtThermSuggest = -1.0;
        number[] params;
        while (_t < _duration) {
            esk(gd, _dt, dtChemSuggest, dtThermSuggest, params);
            _t+=_dt;
        }
        double massfsum=0.0;
        foreach(number eachmassf;gd.massf) {
            massfsum += eachmassf;
        }
        assert(approxEqual(massfsum, 1.0, 1e-2), failedUnitTest());
        return 0;

    }
}

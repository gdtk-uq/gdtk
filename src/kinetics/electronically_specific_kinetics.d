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
// TODO: Brad S. Fix state solver to cleanly compile.
//import kinetics.electronic_state_solver;

final class ElectronicallySpecificKinetics : ThermochemicalReactor {
public:

    this(GasModel gmodel)
    {
        super(gmodel);
        _numden.length = gmodel.n_species;
        _numden0.length = gmodel.n_species;
        // TODO: Brad, think about one-time initialisation of your module at this point.
        // kinetics.electronic_state_solver.Init();
    }

    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[] params) //not sure about this line
    {
        _gmodel.massf2numden(Q, _numden0);

        //update the number density vector for time t+t_interval
        //by solving the kinetics ODE (electron temp frozen)
        
        // TODO: Brad S. Reinstate when state solver is doing what you'd like.
        //Electronic_Solve(_numden0, _numden, Q.T_modes[0], tInterval);
        _gmodel.numden2massf(_numden, Q);
    }

private:
    number[] _numden;
    number[] _numden0; 
}


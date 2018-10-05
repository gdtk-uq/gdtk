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

        foreach(int i;0 .. _gmodel.n_species - 2) { //number density in state solver in #/cm^3
            _numden_input[i] = _numden[i] / 1e6;
        }
        //update the number density vector for time t+t_interval
        //by solving the kinetics ODE (electron temp frozen)
        Electronic_Solve(_numden_input, _numden_output, Q.T_modes[0], tInterval);

        foreach (int i; 0 .. _gmodel.n_species - 2) {//number density in state solver in #/cm^3
            _numden[i] = _numden_output[i] * 1e6;
        } 

        //Coefficients for Molecule Update require mol/cm^3
        foreach(int i; 0 .. _gmodel.n_species){ //convert numden to mol/cm^3
            _numden[i] = _numden[i] / (Avogadro_number*1e6);
        }
        Molecule_Update(Q, tInterval);

        foreach(int i; 0 .. _gmodel.n_species){ //convert mol/cm^3 to numden
            _numden[i] = _numden[i] * Avogadro_number*1e6;
        }

        _gmodel.numden2massf(_numden,Q);
    }

private:
    number[] _numden; //Total set of species including static N2 and O2
    number[] _numden_input; //Input for the electronic state solver, exclusive of N2 and O2
    number[] _numden_output; //Output for the state solver
    number N_sum;
    number O_sum;
    number N2_step_change;
    number O2_step_change;

    //Arhenius coefficients for N2 and O2 dissociation
    //In order: N2, N, O2, O
    //within each group: scale, A, n, E
    double[][] ArrN2_fr = [[2.5, 1.920e+17, -0.50, 113100.00],
                            [1.0, 4.150e+22, -1.50, 113100.00],
                            [1.0, 1.920e+17, -0.50, 113100.00],
                            [1.0, 1.920e+17, -0.50, 113100.00]];
    
    double[][] ArrN2_br = [[2.5, 1.090e+16, -0.50, 0.0],
                            [1.0, 2.320e+21, -1.50, 0.0],
                            [1.0, 1.090e+16, -0.50, 0.0],
                            [1.0, 1.090e+16, -0.50, 0.0]];

    double[][] ArrO2_fr = [[2.0, 3.610e+18, -1.00, 59400.00],
                            [1.0, 3.610e+18, -1.00, 59400.00],
                            [9.0, 3.610e+18, -1.00, 59400.00],
                            [25.0, 3.610e+18, -1.00, 59400.00]];

    double[][] ArrO2_br = [[2.0, 3.010e+15, -0.50, 0.0],
                            [1.0, 3.010e+15, -0.50, 0.0],
                            [9.0, 3.010e+15, -0.50, 0.0],
                            [25.0, 3.010e+15, -0.50, 0.0]];
    @nogc
    number N2_fr(GasState Q) //must only be called when in conc mol/cm^3, not massf
    {
        number rate_coef(int i){
            return ArrN2_fr[i][0]*ArrN2_fr[i][1]*(Q.T^^ArrN2_fr[i][2])*exp(-ArrN2_fr[i][3]/Q.T);
        }
        
        N_sum = 0.0;
        O_sum = 0.0;
        foreach(int i; 0 .. 8) {
            N_sum+=_numden[i];
            O_sum+=_numden[i+9];
        }
        
        return rate_coef(0)*_numden[19]*_numden[19] + rate_coef(1)*_numden[19]*N_sum + 
                rate_coef(2)*_numden[19]*_numden[20] + rate_coef(3)*_numden[19]*O_sum;
    }
    @nogc
    number N2_br(GasState Q) //must only be called when in conc mol/cm^3, not massff
    {
        number rate_coef(int i){
            return ArrN2_br[i][0]*ArrN2_br[i][1]*(Q.T^^ArrN2_br[i][2])*exp(-ArrN2_br[i][3]/Q.T);
        }
        
        N_sum = 0.0;
        O_sum = 0.0;
        foreach(int i; 0 .. 8) {
            N_sum+=_numden[i];
            O_sum+=_numden[i+9];
        }

        return rate_coef(0)*N_sum*N_sum*_numden[19] + rate_coef(1)*N_sum*N_sum*N_sum + 
                rate_coef(2)*N_sum*N_sum*_numden[20] + rate_coef(3)*N_sum*N_sum*O_sum;
    }
    @nogc
    number O2_fr(GasState Q) //must only be called when in conc mol/cm^3, not massf
    {
        number rate_coef(int i){
            return ArrO2_fr[i][0]*ArrO2_fr[i][1]*(Q.T^^ArrO2_fr[i][2])*exp(-ArrO2_fr[i][3]/Q.T);
        }
        
        N_sum = 0.0;
        O_sum = 0.0;
        foreach(int i; 0 .. 8) {
            N_sum+=_numden[i];
            O_sum+=_numden[i+9];
        }
        
        return rate_coef(0)*_numden[20]*_numden[19] + rate_coef(1)*_numden[20]*N_sum + 
                rate_coef(2)*_numden[20]*_numden[20] + rate_coef(3)*_numden[20]*O_sum;
    }
    @nogc
    number O2_br(GasState Q) //must only be called when in conc mol/cm^3, not massf
    {
        number rate_coef(int i){
            return ArrO2_br[i][0]*ArrO2_br[i][1]*(Q.T^^ArrO2_br[i][2])*exp(-ArrO2_br[i][3]/Q.T);
        }
        
        N_sum = 0.0;
        O_sum = 0.0;
        foreach(int i; 0 .. 8) {
            N_sum+=_numden[i];
            O_sum+=_numden[i+9];
        }
        
        return rate_coef(0)*O_sum*O_sum*_numden[19] + rate_coef(1)*O_sum*O_sum*N_sum + 
                rate_coef(2)*O_sum*O_sum*_numden[20] + rate_coef(3)*O_sum*O_sum*O_sum;
    }
    @nogc
    number N2_rate(GasState Q) 
    {   
        return N2_fr(Q) - N2_br(Q);
    }
    @nogc
    number O2_rate(GasState Q)
    {
        return O2_fr(Q) - O2_br(Q);
    }
    @nogc
    void Molecule_Update(GasState Q, double tInterval)
    {
        double _dt = tInterval/100.0;
        foreach(int i; 0 .. 100){
            N2_step_change = N2_rate(Q);
            O2_step_change = O2_rate(Q);
            _numden[19] += -_dt*N2_step_change;
            _numden[20] += -_dt*O2_step_change;
            _numden[0] += _dt*2*N2_step_change;
            _numden[9] += _dt*2*O2_step_change;
        }
    }
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

        number initial_uNoneq=0.0;

        // gd.massf[] = 0.0;
        // gd.massf[0] = 0.037041674288877; //initialises massf of NI
        // gd.massf[9] = 0.010577876366622; //initialises massf of OI
        // gd.massf[19] = 0.74082290750449; //N2
        // gd.massf[20] = 0.21155752733244; //O2
        // gd.massf[18] = 1.0 - (gd.massf[0] + gd.massf[9] + gd.massf[19] + gd.massf[20]); //tiny massf for free electron
        // gd.massf = [0.0313603, 0.00492971, 0.000741705, 1.06916e-06, 4.90114e-07, 
        //                 2.46998e-07, 9.58454e-08, 6.6456e-07, 6.41328e-06, 0.010005, 
        //                 0.000565079, 8.59624e-06, 2.58411e-07, 9.00322e-08, 5.80925e-08, 
        //                 3.67871e-08, 9.06483e-08, 4.16313e-07, 1.4773e-08, 0.740823, 0.211558];
        gd.massf[] = 0;
        gd.massf[18] = 1e-8;
        gd.massf[19] = 0.78;
        gd.massf[20] = 1.0 - 0.78 - 1e-8;
        gd.p = 1000000.0;
        gd.T = 7000.0;
        gd.T_modes[0] = 10000.0;
        gm.update_thermo_from_pT(gd);
        

        lua_close(L);

        double _dt = 1e-9;
        double _duration = 1e-7;
        double _t = 0.0;

        // auto L2 = init_lua_State();
        // doLuaFile(L2, "sample-input/electronic_composition.lua");   

        ElectronicallySpecificKinetics esk = new ElectronicallySpecificKinetics(gm);

        double dtChemSuggest = -1.0;
        double dtThermSuggest = -1.0;
        number[maxParams] params;
        while (_t < _duration) {
            
            esk(gd, _dt, dtChemSuggest, dtThermSuggest, params);
            //update U from change in energy modes
            gm.update_thermo_from_rhou(gd);
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

//arrhenius mol/cm^3



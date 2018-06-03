/**
 * fuel_air_mix.d
 *
 * Fuel+Air->Products reacting gas as used by JJ.
 *
 * This module deatils with the kinetics of change
 * from fuea+air to products which is inferred from
 * the turbulent mixing.
 *
 * Authors: JJ Hoste, Peter J. and Rowan G.
 * Version: 2017-04-22: initial shell copied from powers-aslam-gas.
 */

module kinetics.fuel_air_mix_kinetics;

import std.stdio : writeln;
import std.math;
import std.algorithm;
import std.conv;
import nm.complex;
import nm.number;

import gas;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;
import kinetics.reaction_mechanism;

final class MixingLimitedUpdate : ThermochemicalReactor {
    ReactionMechanism rmech; //only used if we desire to limit the reaction rate
    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        // We need to pick a number of pieces out of the gas-model file, again.
        // Although they exist in the GasModel object, they are private.
        //lua_State *L;
        auto L = init_lua_State();
        
        // get details about type of species: fuel, ox, prod or inert
        doLuaFile(L, fname); //reading the reaction file here
        _n_species = gmodel.n_species;
        lua_getglobal(L, "speciesType");
        getArrayOfStrings(L, LUA_GLOBALSINDEX, "speciesType", _species_type);
        lua_pop(L, 1); // dispose of the table
        foreach(isp;0.._n_species){
            if(_species_type[isp]=="fuel"){_n_fuel=_n_fuel+1;}
            if(_species_type[isp]=="ox"){_n_ox=_n_ox+1;}
            if(_species_type[isp]=="prod"){_n_prod=_n_prod+1;}
        }
        _n_reacting=_n_fuel+_n_ox+_n_prod;
        writeln("nb of fuel species = ",_n_fuel);
        writeln("nb of ox species = ",_n_ox);
        writeln("nb of prod species = ",_n_prod);

        // get settings for EDM model
        lua_getglobal(L, "FuelAirMix");
        // Now, pull out the numeric value parameters.
        _A_edm = getDouble(L, -1, "Aedm");
        _B_edm = getDouble(L, -1, "Bedm");
        _laminar_limit = getBool(L, -1, "lamLimit");
        lua_pop(L, 1); // dispose of the table
        //writeln(_laminar_limit);

        // if true we want to limit the reaction rate with the "no-model" one
        // in that case we initialize the ReactionMechanism object
        if (_laminar_limit) {
            lua_getglobal(L, "config");
            lua_getfield(L, -1, "tempLimits");
            lua_getfield(L, -1, "lower");
            double T_lower_limit = lua_tonumber(L, -1);
            lua_pop(L, 1);
            lua_getfield(L, -1, "upper");
            double T_upper_limit = lua_tonumber(L, -1);
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_pop(L, 1);
            lua_getglobal(L, "reaction");       
            rmech = createReactionMechanism(L, gmodel, T_lower_limit, T_upper_limit);
            lua_rawgeti(L, -1, 1); //only 1 reaction called reaction[1]
            lua_pop(L, 1);
        }

        lua_getglobal(L, "reaction");   
        lua_rawgeti(L, -1, 1); //only 1 reaction called reaction[1]
        int[] reacCoeffs; // required to set stochiometric mass ratio
        getArrayOfInts(L, -1, "reacCoeffs", reacCoeffs);
        writeln(reacCoeffs,reacCoeffs[0]);
        if (_n_reacting > 3) {
            int[] prodCoeffs;
            _nu_W.length=_n_reacting;
            getArrayOfInts(L, -1, "prodCoeffs", prodCoeffs);
            foreach(isp; 0 .. _n_fuel+_n_ox) {
                writeln(isp);
                _nu_W[isp]= 1.0 * reacCoeffs[isp]*gmodel.mol_masses[isp]; 
            }
            foreach(isp; 0 .. _n_prod) {
                _nu_W[_n_fuel+_n_ox+isp]= -1.0 * prodCoeffs[isp]*gmodel.mol_masses[isp+_n_fuel+_n_ox]; 
            }
            lua_pop(L, 1);
        } // else we don't need _nu_W in further computations
        lua_pop(L, 1); // I believe needed to get out 1 step back to all reactions
        lua_close(L);
        writeln("stochiometric mass ratio ",_sRatio);
        set_stochiometricRatio(reacCoeffs,gmodel.mol_masses);
        writeln("stochiometric mass ratio ",_sRatio);
    } // end of constructor


    void set_stochiometricRatio(int[] stochCoeff, double[] mol_masses)
    {   
        double num=0.0,denom=0.0;       
        foreach (isp; 0 .. _n_fuel) {
            denom = denom + mol_masses[isp]*stochCoeff[isp];
        }
        _nuF_WF = denom; //only needed when more than 3 reacting species/prod

        foreach (isp; 0 .. _n_ox) {
            num = num + mol_masses[_n_fuel+isp]*stochCoeff[_n_fuel+isp];
        }
        _sRatio=num/denom;
    }

    void get_massf(GasState gas, ref number[3] Ys)
    {
        // pass Ys by ref required to fill local variable
        // initialization required in case of cumulative sum below
        Ys=[to!number(0.0), to!number(0.0), to!number(0.0)];
        int i=0;
        if (_n_fuel > 1) {
            foreach (isp; 0 .. _n_fuel) {
                Ys[0]=Ys[0]+gas.massf[isp];
                i=_n_fuel-1;
            }
        } else {
            Ys[0]=gas.massf[0];
            writeln(" massf sp0 ",Ys[0]);
        }
        if (_n_ox>1) {
            foreach(isp;0.._n_ox){
                Ys[1]=Ys[1]+gas.massf[_n_fuel+isp];
            }
        } else {
            Ys[1]=gas.massf[_n_fuel];
            //writeln(" massf sp1 ",Ys[1]);
        }
        if(_n_prod>1) {
            foreach(isp;0.._n_prod){
                Ys[2]=Ys[2]+gas.massf[_n_fuel+_n_ox+isp];
                //writeln(" massf prod ",Ys[2]);
            }
        } else {
            Ys[2]=gas.massf[_n_fuel+_n_ox];
        }
    }

    number get_omegaDotF(number omega, number[3] Ys)
    {
        return -1.0 * _A_edm / _tau_edm * min(Ys[0],Ys[1]/_sRatio,_B_edm*Ys[2]/(1.0+_sRatio));
    }

    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[] params)
    {
        double local_dt_global=tInterval;
        number omega = params[0];
        _tau_edm=1.0/(_beta_star*omega);
        number _omega_dot_F = 0.0;
        number omega_dot_O = 0.0;
        number omega_dot_P = 0.0;
        number[] rates; // for laminar= "no-model" combustion
        rates.length = _n_species;
        //writeln(rmech.) 
        number[3] _Y_s;
        if(_n_reacting > 3){
            get_massf(Q,_Y_s);
            //writeln("species massf ",_Y_s);
            _omega_dot_F = get_omegaDotF(omega,_Y_s);
            //writeln("fuel RR ",_omega_dot_F);
            //writeln("nu F x W F = ",_nuF_WF);
            //writeln("nu  x W  = ",_nu_W);
            if(_laminar_limit){
                number[] conc;
                conc.length=_n_species; 
                rmech.eval_rate_constants(Q);  
                _gmodel.massf2conc(Q,conc);
                writeln("concentration ",conc);
                rmech.eval_rates(conc,rates); 
                writeln("laminar RR (1/s)= ",rates);
                number[] lamRR; //in kg/(m^3 . s)
                lamRR.length=rates.length;
                number sumLamRR=0.0;
                foreach (isp; 0 .. _n_species) {
                    lamRR[isp]= _gmodel.mol_masses[isp]/Q.rho* rates[isp];
                    sumLamRR=sumLamRR+lamRR[isp];
                }
                writeln("laminar RR (1/s) = ",lamRR);
                writeln("sum of laminar RR (1/s) = ",sumLamRR);
                writeln("Fuel RR (1/s) before lam eval= ",_omega_dot_F);
                _omega_dot_F=-1.0*(fmin(fabs(_omega_dot_F),fabs(rates[0]))); 
            }
            writeln("Fuel RR (1/s)= ",_omega_dot_F);
            foreach (isp; 0 .. _n_reacting) { //don't account for inert species
                Q.massf[isp]=Q.massf[isp] + local_dt_global * (_nu_W[isp])/ _nuF_WF* _omega_dot_F;
            }
            //number sum=0.0;
            //foreach(isp;0.._n_species){
                //sum=sum+Q.massf[isp];                                 
            //}
            //writeln("sum of species = ",sum);
        } else {
            _Y_s[0]=Q.massf[0];
            _Y_s[1]=Q.massf[1];
            _Y_s[2]=Q.massf[2];
            _omega_dot_F = get_omegaDotF(omega,_Y_s);   // units 1/s in here 
            if(_laminar_limit) {
                number[] conc;
                conc.length=_n_species; 
                rmech.eval_rate_constants(Q);  
                _gmodel.massf2conc(Q,conc);
                writeln("concentration ",conc);
                rmech.eval_rates(conc,rates); 
                writeln("laminar RR (1/s)= ",rates);
                number[] lamRR; //in kg/(m^3 . s)
                lamRR.length=rates.length;
                number sumLamRR=0.0;
                foreach (isp; 0 .. _n_species) {
                    lamRR[isp]= _gmodel.mol_masses[isp]/Q.rho* rates[isp];
                    sumLamRR=sumLamRR+lamRR[isp];
                }
                writeln("laminar RR (1/s) = ",lamRR);
                writeln("sum of laminar RR (1/s) = ",sumLamRR);
                writeln("Fuel RR (1/s) before lam eval= ",_omega_dot_F);
                _omega_dot_F=-1.0*(fmin(fabs(_omega_dot_F),fabs(rates[0]))); // units are same 
                // NOTE: when fuel/ox almost zero the laminar RR will give positive values for
                //       these reaction rates BUT the EDM one is much smaller and dominates     
                //               only reason for laminar is to avoid reaction where T too low and not enough
                //               mixing occurs.                 
            } // end if /else
            omega_dot_O = _sRatio * _omega_dot_F;
            omega_dot_P = - (1.0+_sRatio) *_omega_dot_F;
            //writeln("Fuel RR (1/s)= ",_omega_dot_F);
            //writeln(omega_dot_O);
            //writeln(omega_dot_P);
            Q.massf[0]=_Y_s[0] + local_dt_global * _omega_dot_F;
            // following is actuall still correct for laminar RR
            Q.massf[1]=_Y_s[1] + local_dt_global * _sRatio*_omega_dot_F;
            Q.massf[2]=_Y_s[2] + local_dt_global * - (1.0+_sRatio) *_omega_dot_F;               
        } //end if/ else
    } //end function opCall

private:
    // EDM model constants
    double _A_edm=4.0; 
    double _B_edm=0.5;
    number _tau_edm;
    // boolean for use of laminar reaction rate limit
    bool _laminar_limit=false;
    static immutable _beta_star=0.09;
    string[] _species_type; // fuel, ox, prod, inert

    int _n_fuel, _n_ox, _n_prod, _n_species, _n_reacting;
    double _sRatio; // mass stochiometric ratio
    // parameters used if nb of species in reaction (_n_reacting) > 3
    double _nuF_WF; // stochiometric coeff fuel x molar mass fuel
    double[] _nu_W; // the above for every species
} // end class MixingLimitedUpdate


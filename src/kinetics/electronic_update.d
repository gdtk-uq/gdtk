/** 
 * Author: Brad Semple
 * Date: 2019-01-14
 *
 * This is a restructuring of my undergraduate thesis project on the 
 * electronically specific distribution of electronc energy levels 
 * within a gas.
 * 
 * This module contains the electronically specific gas object, which
 * provides the data storage and update routines to the electronic distribution
 * within a gas species.
 * 
 * The intent of this object is to represent a single chemical species
 * (such as N, or O2+) in a gas mixture, and to break it up into 
 * multiple electronic energy levels, which are themselves groupings of 
 * individually allowed energy levels (although this is invisible to eilmer).
 * 
 * Separating the previous electronic solver allows for the better separation
 * of nitrogen and oxygen species, while also allowing for the easier expansion
 * into other chemical species for the electronic model.
 * 
 *
 */

module kinetics.electronic_update;

import std.stdio;
import std.file;
import std.string;
import std.format;
import std.conv : to;
import std.math;
import std.algorithm.iteration : sum;

import nm.bbla;
import nm.bdfLU;
import ntypes.complex;
import nm.number;
import util.msg_service;

import gas;
import gas.electronically_specific_gas;
import util.lua;
import util.lua_service;

class ElectronicUpdate
{
public:
    @property @nogc number electronic_energy() { return EnergyInNonEq(); }
    this(string rates_file, size_t[] index_location, ElectronicallySpecificGas _fullAirModel) //Index range contains: [ground energy level index, ionisation index]
    {
        auto L = init_lua_State();
        doLuaFile(L, rates_file);
        chem_name = getString(L, "gas");
        n_bins = getInt(L, "number_of_bins");
        ion_bin = getBool(L, "ionised_bin");
        n_species = n_bins;

        first_bin_index = index_location[0];
        ion_index = index_location[1];

        if (ion_bin) {
            n_species += 1;
        }
        n_reactions = 0;
        foreach(i; 0 .. n_species) {
            n_reactions += i;
        }


        lua_getglobal(L,"bins");
        foreach (isp; 0 .. n_species) {
            lua_rawgeti(L, -1, isp);
            bins ~= new ElectronicBin(L);
            lua_pop(L, 1);
        }

        foreach(i; 0 .. n_bins) {
            index_range ~= index_location[0] + i;
        }
        if (ion_bin) {
            index_range ~= index_location[1];
        }

        lua_getglobal(L, "reactions");
        foreach (isp; 0 .. n_reactions) {
            lua_rawgeti(L, -1, isp);
            reactions ~= new ElectronicReaction(L, index_range);
            lua_pop(L, 1);
        }

        rate_vector.length = n_species;
        state.length = n_species;
        guess_state.length = n_species;
        prev_state.length = n_species;
        update_vector.length = n_species;
        RHSvector.length = n_species;
        pivot.length = n_species;
        jacobian.length = n_species;
        Asolve = new Matrix!number(n_species,n_species);
        foreach (isp; 0 .. n_species) {
            jacobian[isp].length = n_species;
        }
        _AirModel = _fullAirModel;
        _numden.length = _AirModel.n_species;
    }

    @nogc final void Update(ref GasState Q, double duration)
    {   
        
        endtime = duration;
        Te = Q.T_modes[0];
        _AirModel.massf2numden(Q, _numden);

        foreach(isp; 0 .. n_species) {
            state[isp] = _numden[index_range[isp]] / 1e6; //model developed in cm^-3 number density
        }
        Ne = _numden[$-1] / 1e6;

        step_for_duration();
        foreach (isp; 0 .. n_species) {
            _numden[index_range[isp]] = state[isp] * 1e6;
        }

        _AirModel.numden2massf(_numden, Q);


    }
private:
    ElectronicallySpecificGas _AirModel;
    bool ion_bin;
    string chem_name;
    int n_bins, n_species, n_reactions;
    int newton_steps = 5;
    int[] pivot;

    size_t[] index_range;
    double dt = 1e-7;
    double t, endtime;
    size_t ion_index, first_bin_index;

    ElectronicBin[] bins;
    ElectronicReaction[] reactions;
    number[][] jacobian;
    number[] _numden, rate_vector, state, guess_state, prev_state, update_vector, RHSvector;
    Matrix!number Asolve;
    number Ne, Te;

    @nogc
    int step_for_duration()
    {
        t = 0.0;
        if (dt<endtime) {
            First_Step();
        } else {
            dt = endtime;
            First_Step();
            return 0;
        }

        while ((t+dt)<endtime) {
            Step();
        }
        if ((endtime-t) > 1e-15) {
            dt = endtime-dt;
            Step();
        }
        return 0;
    }

    @nogc
    void First_Step()
    {   
        guess_state[] = state[];

        foreach (n; 0 .. newton_steps) {

            update_rates_with_newton_guess_state();

            update_rate_vector();
            update_jacobian();

            debug {
                BDF1(update_vector, RHSvector, Asolve, pivot, dt/100.0,state,guess_state,rate_vector,jacobian);
            }
            if (ion_bin) { //Update electron number density
                Ne = (_numden[$-1] / 1e6) + (guess_state[$-1] - state[$-1]);
            }
        }

        prev_state[] = state[];
        state[] = guess_state[];

        foreach(n; 0 .. newton_steps) {

            update_rates_with_newton_guess_state();

            update_rate_vector();
            update_jacobian();
            debug {
                Var_BDF2(update_vector, RHSvector, Asolve, pivot,(99.0/100.0)*dt, dt/100.0, prev_state, state, guess_state, rate_vector,jacobian);
            }
            if (ion_bin) { //Update electron number density
                Ne = (_numden[$-1] / 1e6) + (guess_state[$-1] - state[$-1]);
            }
        }

        prev_state[] = state[];
        state[] = guess_state[];

        t+=dt;
    }

    @nogc
    void Step()
    {
        guess_state[] = state[];

        foreach (n; 0 .. newton_steps) {

            update_rates_with_newton_guess_state();

            debug {
                BDF2(update_vector, RHSvector, Asolve, pivot, dt, prev_state, state, guess_state,rate_vector,jacobian);
            }
            if (ion_bin) { //Update electron number density
                Ne = (_numden[$-1] / 1e6) + (guess_state[$-1] - state[$-1]);
            }
        }

        prev_state[] = state[];
        state[] = guess_state[];

        t+=dt;
    }

    @nogc
    void update_rates_with_newton_guess_state()
    {
        foreach(isp; 0 .. n_reactions) {
            reactions[isp].updateRates(guess_state, Te, Ne);
        }

        update_rate_vector();
        update_jacobian();
    }

    @nogc
    void update_jacobian()
    {
        foreach(isp; 0 .. n_species) {
            foreach(jsp; 0 .. n_species) {
                jacobian[isp][jsp] = 0.0;
            }
        }

        foreach(isp; 0 .. n_reactions) {
            size_t i = reactions[isp].i_index;
            size_t j = reactions[isp].j_index;
            jacobian[i][] += reactions[isp].low_dRdN[]; // dR_i / dN for each N
            jacobian[j][] -= reactions[isp].low_dRdN[]; // dR_j / dN for each N (the negative of the lower index)
        }

    }

    @nogc
    void update_rate_vector()
    {   
        foreach(isp; 0 .. n_species) {
            rate_vector[isp] = 0.0;
        }
        
        foreach (isp; 0 .. n_reactions) {
            size_t i = reactions[isp].i_index;
            size_t j = reactions[isp].j_index;
            rate_vector[i] += reactions[isp].f_rate; // forward rate for each reaction from lower species i
            rate_vector[j] -= reactions[isp].f_rate; // negative of forward rate to upper species j
        }
    }

    @nogc 
    number EnergyInNonEq()
    {
        return to!number(0);
    }
}



class ElectronicBin
{
    @property @nogc string name() const { return _name; }
    @property @nogc int binlevel() const { return _binlevel; }
    @property @nogc int lowerLevel() const {return _lowerLevel; }
    @property @nogc int upperLevel() const {return _upperLevel; }
    @property @nogc int group_degeneracy() const { return _group_degeneracy; }
    @property @nogc double electronic_energy() const {return _electronic_energy; }

    this(lua_State *L)
    {
        _name = getString(L, -1, "name");
        _binlevel = getInt(L, -1, "bin_level");
        _mol_mass = getDouble(L, -1, "M");
        _lowerLevel = getInt(L, -1, "lower_level");
        _upperLevel = getInt(L, -1, "upper_level");
        _group_degeneracy = getInt(L, -1, "group_degeneracy");
        _electronic_energy = getDouble(L, -1, "group_energy")*electron_volt_energy*Avogadro_number/_mol_mass;
    }

private:
    string _name;
    int _binlevel, _lowerLevel, _upperLevel, _group_degeneracy;
    double _electronic_energy, _mol_mass;
}

class ElectronicReaction
{
    @property @nogc number f_rate() {return _f_rate; }
    @property @nogc number[] low_dRdN() {return _low_dRdN; }
    size_t i_index;
    size_t j_index;

    this(lua_State *L, size_t[] bin_index_range)
    {
        _index_range.length = bin_index_range.length;
        _index_range[] = bin_index_range[];
        lower_bin = getInt(L, -1, "lower_bin");
        upper_bin = getInt(L, -1, "upper_bin");
        _lower_index = lower_bin - 1;
        if (upper_bin == 0) {
            _upper_index = bin_index_range.length - 1;
        } else {
            _upper_index = upper_bin-1;
        }
        A = getDouble(L, -1, "A");
        n = getDouble(L, -1, "n");
        E = getDouble(L, -1, "E");
        G1 = getDouble(L, -1, "G1");
        G2 = getDouble(L, -1, "G2");
        G3 = getDouble(L, -1, "G3");
        G4 = getDouble(L, -1, "G4");
        G5 = getDouble(L, -1, "G5");
        _low_dRdN.length = bin_index_range.length;
        foreach(isp; 0 .. _low_dRdN.length) {
            _low_dRdN[isp] = 0.0;
        }

        i_index = _lower_index;
        j_index = _upper_index;
    }

    @nogc final void updateRates(in number[] numden, number Te, number Ne)
    {   
        //Develop reaction rate coefficients from rate fitting parameters
        z = 10000.0 / Te;
        k_ij = A*((Te)^^n)*exp(-E/Te);
        k_ji = k_ij/(exp((G1/z) + G2 + G3*log(z) + G4*z + G5*(z^^2)));
        
        //Develop rates of change of the lower and upper species in the reaction
        if (upper_bin == 0) { //for an ionising reaction
            _f_rate = Ne*(Ne*k_ji*numden[_upper_index] - k_ij*numden[_lower_index]);
        } else { // for an excitation reaction
            _f_rate = Ne*(k_ji*numden[_upper_index] - k_ij*numden[_lower_index]);
        }

        //Develop analytical dR/dN for reaction for each N - many are zero.
        
        _low_dRdN[_lower_index] = -k_ij*Ne;
        if (upper_bin == 0) {
            _low_dRdN[_upper_index] = k_ji*(Ne*Ne + 2*numden[$-1]*Ne);
        } else {
            _low_dRdN[_upper_index] = k_ji*Ne;
        }
    }

private:
    int lower_bin, upper_bin;
    size_t _lower_index, _upper_index;
    size_t[] _index_range;
    double A, n, E, G1, G2, G3, G4, G5;
    number _f_rate, z, Te, Ne, k_ij, k_ji;
    number[] _low_dRdN;
}

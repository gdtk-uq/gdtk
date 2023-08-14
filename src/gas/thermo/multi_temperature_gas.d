/*
 * Multi-temperature thermodynamic module
 * 
 * The MultiTemperatureGasMixture provides a thermodynamic model
 * for a gas mixture with an arbitrary number of temperatures.
 * This will generally consist of a vibration temperature for 
 * each molecule
 * 
 * It is assumed that translation/rotation are stored in GasState.T
 * The remaining energy modes are stored in GasState.T_modes
*/

module gas.thermo.multi_temperature_gas;

import std.stdio;
import std.math;
import std.string;
import std.conv;
import std.algorithm;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;
import nm.bracketing;

import gas;
import gas.thermo.thermo_model;
import gas.thermo.energy_modes;
import gas.thermo.cea_thermo_curves;

immutable double T_REF = 298.15; // K

class MultiTemperatureGasMixture : ThermodynamicModel {
    this(lua_State *L, string[] species_names, string[] energy_mode_names)
    {
        // allocate memory
        _n_species = to!int(species_names.length);
        _n_modes = to!int(energy_mode_names.length);
        _R.length = _n_species;
        _Delta_hf.length = _n_species;
        _CpTR.length = _n_species;
        _thermo_curves.length = _n_species;
        _energy_modes.length = _n_modes;
        _energy_modes_isp.length = _n_modes;
        _reference_energies.length = _n_modes;
        _max_iterations.length = _n_modes;
        _tolerances.length = _n_modes;
        
        // set up the species data
        foreach (isp, sp_name; species_names) {
            _species_indices[sp_name] = to!int(isp);
            if (sp_name == "e-") _electron_idx = to!int(isp);
            lua_getglobal(L, "db");
            lua_getfield(L, -1, sp_name.toStringz);
            double m = getDouble(L, -1, "M");
            _R[isp] = R_universal/m;
            lua_getfield(L, -1, "thermoCoeffs");
            _thermo_curves[isp] = new CEAThermoCurve(L, _R[isp]);
            lua_pop(L, 1);
            _Delta_hf[isp] = _thermo_curves[isp].eval_h(to!number(T_REF)); //getDouble(L, -1, "Hf");
            string type = getString(L, -1, "type");
            switch (type) {
            case "electron":
                _CpTR[isp] = 0.0;
                break;
            case "atom" :
                _CpTR[isp] = (5./2.)*_R[isp];
                break;
            case "molecule":
                string molType = getString(L, -1, "molecule_type");
                _CpTR[isp] = (molType == "linear") ? (7./2.)*_R[isp] : (8./2.)*_R[isp];
                break;
            default:
                string msg = "MultiTemperatureGas: error trying to match particle type.\n";
                throw new Error(msg);
            }
            lua_pop(L, 1);
            lua_pop(L, 1);
        }

        // set up energy modes
        foreach (i_mode, mode_name; energy_mode_names) {
            string[] components; 
            lua_getglobal(L, "db");
            lua_getfield(L, -1, "modes");
            getArrayOfStrings(L, -1, mode_name, components);
            lua_pop(L, 1); // modes
            _max_iterations[i_mode] = 200;
            _tolerances[i_mode] = 1e-8;
            _energy_modes[i_mode].length = to!int(components.length);
            _energy_modes_isp[i_mode].length = to!int(components.length);
            _reference_energies[i_mode].length = to!int(components.length);
            foreach (i_comp, component; components) {
                string[] component_tokens;
                component_tokens = component.split(":");
                string species = component_tokens[0];
                int isp = _species_indices[species];
                _energy_modes_isp[i_mode][i_comp] = isp;
                if (isp == _electron_idx){
                    _electron_mode = to!int(i_mode);
                }
                string energy_type = component_tokens[1];
                lua_getfield(L, -1, species.toStringz);
                switch (energy_type) {
                    case "vib":
                        _energy_modes[i_mode][i_comp] = create_vibrational_energy_model(L, _R[isp]);
                        break;
                    case "electronic":
                        _energy_modes[i_mode][i_comp] = create_electronic_energy_model(L, _R[isp]);
                        break;
                    default:
                        throw new Error("Unknown energy type");
                }
                lua_pop(L, 1); // species data entry

                // evaluate energy at reference temperature
                _reference_energies[i_mode][i_comp] = _energy_modes[i_mode][i_comp].energy(to!number(T_REF));
            }
            lua_pop(L, 1); // db
        }
    }

    override void updateFromPU(ref GasState gs) {
        // update the thermo state from pressure and internal energy.
        // Requires the following properties to be set:
        //   + u
        //   + p
        //   + u_modes

        _update_temperature(gs);
        _update_T_modes(gs);
        _update_density(gs);
        _update_p_e(gs);
    }
    
    @nogc void updateFromPT(ref GasState gs){
        // update the thermo state from pressure and temperature.
        // Requires the following properties to be set:
        //   + T
        //   + T_modes
        //   + pressure

        _update_density(gs);
        gs.u = _mixture_transrot_energy(gs);
        _update_u_modes(gs);
        _update_p_e(gs);
    }

    @nogc void updateFromRhoU(ref GasState gs){
        // update the thermo state from rho and internal energy.
        // Requires the following properties to be set:
        //   + rho
        //   + u
        //   + u_modes

        _update_temperature(gs);
        _update_T_modes(gs);
        _update_pressure(gs);
    }

    @nogc void updateFromRhoT(ref GasState gs){
        // update the thermo state from rho and temperature.
        // Requires the following properties to be set:
        //   + rho
        //   + T
        //   + T_modes

        _update_pressure(gs);
        gs.u = _mixture_transrot_energy(gs);
        _update_u_modes(gs);
    }
    
    @nogc void updateFromRhoP(ref GasState gs){
        // Update the thermo state from rho and pressure.
        // Requires the following properties to be set:
        //   + rho
        //   + p
        //   + u_modes

        _update_T_modes(gs);
        gs.T = _transrot_temp_from_rho_p(gs);
        gs.u = _mixture_transrot_energy(gs);
        _update_p_e(gs);
    }

    @nogc void updateFromPS(ref GasState gs, number s){
        throw new Error("Not implemented");
    }

    @nogc void updateFromHS(ref GasState gs, number h, number s){
        throw new Error("Not implemented");
    }
    
    @nogc void updateSoundSpeed(ref GasState gs){
        // compute the frozen sound speed
        number gamma = dhdTConstP(gs) / dudTConstV(gs);
        gs.a = sqrt(gamma * gs.p / gs.rho);
    }

    // Methods related to computing thermo derivatives.
    @nogc number dudTConstV(in GasState gs){
        number Cv = _mixture_transrot_Cv(gs);
        foreach (imode; 0 .. _n_modes) {
            Cv += _mixture_Cv_in_mode(gs, imode);
        }
        return Cv;
    }

    @nogc number dhdTConstP(in GasState gs){
        number Cp = 0.0;
        // Note that for internal modes, Cv = Cp
        foreach (imode; 0 .. _n_modes){
            Cp += _mixture_Cv_in_mode(gs, imode);
        }
        foreach (isp; 0 .. _n_species) {
            Cp += gs.massf[isp] * _CpTR[isp];
        }
        return Cp;
    }

    @nogc number dpdrhoConstT(in GasState gs){
        number sum = 0.0;
        foreach (isp; 0 .. _n_species) {
            number T = (isp == _electron_idx) ? gs.T_modes[_electron_mode] : gs.T;
            sum += gs.massf[isp] * _R[isp] * T;
        }
        return sum;
    }

    @nogc number gasConstant(in GasState gs){
        return mass_average(gs, _R);
    }

    @nogc number gasConstant(in GasState gs, int isp){
        return to!number(_R[isp]);
    }

    @nogc number internalEnergy(in GasState gs){
        number u = _mixture_transrot_energy(gs);
        foreach (imode; 0 .. _n_modes) {
            u += _mixture_energy_in_mode(gs, imode);
        }
        return u;
    }

    @nogc number energyPerSpeciesInMode(in GasState gs, int isp, int imode){
        return _energy_per_species_in_mode(gs, isp, imode);
    }

    @nogc number enthalpy(in GasState gs){
        number u = internalEnergy(gs);
        return u + gs.p / gs.rho;
    }

    @nogc
    override void enthalpies(in GasState gs, number[] hs)
    {
        foreach(isp; 0 .. _n_species){
            hs[isp] = enthalpyPerSpecies(gs, isp);
        }
    }

    @nogc number enthalpyPerSpecies(in GasState gs, int isp){
        number h = _CpTR[isp]*(gs.T - T_REF) + _Delta_hf[isp];
        foreach (imode; 0 .. _n_modes) {
            h += _energy_per_species_in_mode(gs, isp, imode);
        }
        return h;
    }

    @nogc number enthalpyPerSpeciesInMode(in GasState gs, int isp, int imode){
        if (imode == -1) {
            return _CpTR[isp] * (gs.T - T_REF) + _Delta_hf[isp];
        }

        return _energy_per_species_in_mode(gs, isp, imode);
    }

    @nogc number entropy(in GasState gs){
        throw new Error("Not implemented");
    }

    @nogc number entropyPerSpecies(in GasState gs, int isp){
        throw new Error("Not implemented");
    }

    @nogc number cpPerSpecies(in GasState gs, int isp){
        number Cp = _CpTR[isp];

        // For internal modes, Cp = Cv
        foreach (imode; 0 .. _n_modes) {
            Cp += _Cv_per_species_in_mode(gs, isp, imode);
        }
        return Cp;
    }
    @nogc override void GibbsFreeEnergies(in GasState gs, number[] gibbs_energies)
    {
        number T = gs.T;
        number logT = log(T);
        number logp = log(gs.p/P_atm);

        foreach(isp; 0 .. _n_species) {
            number h = enthalpyPerSpecies(gs, isp);
            number s = _thermo_curves[isp].eval_s(T, logT) - _R[isp]*logp;
            gibbs_energies[isp] = h - T*s;
        }
    }

private:
    int[string] _species_indices;
    int _n_species;
    int _n_modes;
    int _electron_idx = -1;
    int _electron_mode = -1;
    double[] _R;
    number[] _Delta_hf;
    double[] _CpTR;
    InternalEnergy[][] _energy_modes;
    int[][] _energy_modes_isp; // keeps track of which species is contributing
                               // to each energy mode
    CEAThermoCurve[] _thermo_curves;
    number[][] _reference_energies; // The energy of each mode evaluated at T_REF
    double[] _tolerances;
    int[] _max_iterations;

    @nogc void _update_density(ref GasState gs){
        // update the density of the gas, assumes the following are correct:
        //   + massf
        //   + T (and T_modes)
        //   + p

        number denom = 0.0;
        foreach (isp; 0 .. _n_species) {
            number T = (_electron_idx == isp) ? gs.T_modes[_electron_mode] : gs.T;
            denom += gs.massf[isp] * _R[isp] * T;
        }
        gs.rho = gs.p/denom;
    }

    @nogc void _update_temperature(ref GasState gs) {
        // update the translational temperature, assuming the following are correct:
        //   + massf
        //   + u

        number sum_a = 0.0;
        number sum_b = 0.0;
        foreach (isp; 0 .. _n_species) {
            if (isp == _electron_idx) continue;
            sum_a += gs.massf[isp] * (_CpTR[isp]*T_REF - _Delta_hf[isp]);
            sum_b += gs.massf[isp] * (_CpTR[isp] - _R[isp]);
        }
        gs.T = (gs.u + sum_a) / sum_b;
    }

    @nogc void _update_pressure(ref GasState gs) {
        // update the pressure, assuming the following are correct:
        //   + massf
        //   + rho
        //   + T (and T_modes)

        gs.p = 0.0;
        foreach (isp; 0 .. _n_species) {
            number T = (isp == _electron_idx) ? gs.T_modes[_electron_mode] : gs.T;
            gs.p += gs.massf[isp] * gs.rho * _R[isp] * T;
        }
        _update_p_e(gs);
    }

    @nogc void _update_T_modes(ref GasState gs) {
        // update T_modes, assuming the following are correct:
        //   + u_modes
        //   + massf

        foreach (imode; 0 .. _n_modes) {
            gs.T_modes[imode] = _temperature_of_mode(gs, imode);
        }
    }

    @nogc void _update_u_modes(ref GasState gs) {
        // update u_modes, assuming the following are correct:
        //   + T_modes
        //   + massf

        foreach (imode; 0 .. _n_modes) {
            gs.u_modes[imode] = _mixture_energy_in_mode(gs, imode);
        }
    }

    @nogc void _update_p_e(ref GasState gs){
        // update the electron pressure, assuming the following are correct:
        //   + massf
        //   + T_modes

        if (_electron_idx != -1){
            gs.p_e = gs.rho*gs.massf[_electron_idx]*_R[_electron_idx]*gs.T_modes[_electron_mode];
        }
    }

    @nogc number _energy_per_species_in_mode(in GasState gs, int isp, int imode){
        if (imode == _electron_mode && isp == _electron_idx) {
            return 3./2. * _R[_electron_idx] * (gs.T_modes[_electron_mode] - T_REF);
        }

        number energy = to!number(0.0);
        InternalEnergy[] components_in_mode = _energy_modes[imode];
        number [] reference_energies = _reference_energies[imode];
        int[] modes_isp = _energy_modes_isp[imode];
        foreach (icomp, comp; components_in_mode) {
            if (modes_isp[icomp] == isp) {
                energy += comp.energy(gs.T_modes[imode]) - reference_energies[icomp];
            }
        }
        return energy;
    }

    @nogc number _mixture_energy_in_mode(in GasState gs, int imode) {
        return _mixture_energy_in_mode_at_temp(gs, gs.T_modes[imode], imode);
    }

    @nogc number _mixture_energy_in_mode_at_temp(in GasState gs, number T, int imode){
        number energy = to!number(0.0);
        InternalEnergy[] components_in_mode = _energy_modes[imode];
        int[] modes_isp = _energy_modes_isp[imode];
        number[] reference_energies = _reference_energies[imode];
        foreach (icomp, comp; components_in_mode) {
            int isp = modes_isp[icomp];
            energy += gs.massf[isp] * (comp.energy(T) - reference_energies[icomp]);
        }
        if (imode == _electron_mode) {
            energy += gs.massf[_electron_idx] * 3./2. * _R[_electron_idx] * (T - T_REF);
        }
        return energy;
    }

    @nogc number _mixture_Cv_in_mode(in GasState gs, int imode){
        return _mixture_Cv_in_mode_at_temp(gs, gs.T_modes[imode], imode);
    }

    @nogc number _mixture_Cv_in_mode_at_temp(in GasState gs, number T, int imode){
        number Cv = to!number(0.0);
        InternalEnergy[] components_in_mode = _energy_modes[imode];
        int[] modes_isp = _energy_modes_isp[imode];
        foreach (icomp, comp; components_in_mode) {
            int isp = modes_isp[icomp];
            Cv += gs.massf[isp] * comp.Cv(T);
        }
        if (imode == _electron_mode) {
            Cv += gs.massf[_electron_idx] * 3./2. * _R[_electron_idx];
        }
        return Cv;
    }

    @nogc number _Cv_per_species_in_mode(in GasState gs, int isp, int imode) {
        if (isp == _electron_idx && imode == _electron_mode) {
            return to!number(3./2. * _R[_electron_idx]);
        }

        number Cv = to!number(0.0);
        InternalEnergy[] components_in_mode = _energy_modes[imode];
        int[] modes_isp = _energy_modes_isp[imode];
        foreach (icomp, comp; components_in_mode) {
            if (modes_isp[icomp] == isp){
                Cv += comp.Cv(gs.T_modes[imode]);
            }
        }
        return Cv;
    }

    @nogc number _transrot_energy_per_species(in GasState gs, int isp) {
        if (isp == _electron_idx) return to!number(0.0);
        number h = _CpTR[isp] * (gs.T - T_REF) + _Delta_hf[isp];
        return h - _R[isp] * gs.T;
    }

    @nogc number _mixture_transrot_energy(in GasState gs) {
        number energy = to!number(0.0);
        foreach (isp; 0 .. _n_species) {
            energy += gs.massf[isp] * _transrot_energy_per_species(gs, isp);
        }
        return energy;
    }

    @nogc number _transrot_Cv_per_species(int isp){
        if (isp == _electron_idx) return to!number(0.0);
        return to!number(_CpTR[isp] - _R[isp]);
    }

    @nogc _mixture_transrot_Cv(in GasState gs) {
        number Cv = 0.0;
        foreach (isp; 0 .. _n_species) {
            Cv += gs.massf[isp] * _transrot_Cv_per_species(isp);
        }
        return Cv;
    }

    @nogc number _transrot_temp_from_rho_p(in GasState gs) {
        // We assume T_modes is set correctly
        number p_heavy = gs.p;
        if (_electron_idx != -1) {
            p_heavy -= gs.rho*gs.massf[_electron_idx]*_R[_electron_idx]*gs.T_modes[_electron_mode];
        }

        number denom = 0.0;
        foreach (isp; 0 .. _n_species) {
            if (isp == _electron_idx) continue;
            denom += gs.rho * gs.massf[isp] * _R[isp];
        }
        return p_heavy / denom;
    }

    @nogc number _temperature_of_mode(in GasState gs, int imode){
        // Use a Newton-Raphson method to solve for a modal temperature.
        // From experience, this can sometimes fail since we don't
        // necessarily have nicely behaved functions. So we stabalise 
        // it by using bisection alongside. This implementation is based
        // on `rtsafe` from "numerical recipes in C", second edition, page 366.

        immutable int max_iterations = _max_iterations[imode];
        immutable number T_MIN = to!number(10.0);
        immutable number T_MAX = to!number(200_000.0);
        immutable double tol = _tolerances[imode];
        immutable double IMAGINARY_TOL = 1.0e-30;
        immutable number target_u = gs.u_modes[imode];

        // guess that the new temperature is within +/- 5 Kelvin of gs.T_modes[imode]
        // this guess will be adjusted later by a bracketing algorithm
        number Ta = gs.T_modes[imode] - to!number(5.0);
        number Tb = gs.T_modes[imode] + to!number(5.0);


        number zero_func(number T) {
            number u = _mixture_energy_in_mode_at_temp(gs, T, imode);
            return u - target_u;
        }
        
        // try to bracket the root
        if (bracket!(zero_func, number)(Ta, Tb, T_MIN, T_MAX) == -1){
            // maybe this temperature isn't well defined, in which
            // case we should just leave the temperature alone.
            number bath_massf = to!number(0.0);
            int[] isps = _energy_modes_isp[imode];
            foreach (isp; isps) {
                bath_massf += gs.massf[isp];
            }
            if (bath_massf < 1e-15){
                return gs.T_modes[imode];
            }

            // we can't recover from this
            string msg = "Failed to bracket modal temperature";
            debug {
                msg ~= "\n";
                msg ~= format("mode: %d \n", imode);
                msg ~= format("gs: %s \n", gs);
            }
            throw new GasModelException(msg);
        }

        // the function evaluation at the bounds
        number fa, fb;
        fa = zero_func(Ta);
        fb = zero_func(Tb);
        number dT = fabs(Tb - Ta);
        number dT_old = dT;

        // use gs.T_modes[imode] as the initial guess
        number T = gs.T_modes[imode];

        // the function evaluation and derivative at the current guess
        number f = zero_func(T);
        number df = _mixture_Cv_in_mode_at_temp(gs, T, imode);

        bool converged = false;
        foreach (iter; 0 .. max_iterations) {
            bool newton_unstable = ((((T-Tb)*df-f) * ((T-Ta)*df-f)) > 0.0);
            bool newton_slow = (fabs(2.0*f) > fabs(dT_old*df));
            if (newton_unstable || newton_slow) {
                // use bisection
                dT_old = dT;
                dT = 0.5*(Tb-Ta);
                T = Ta + dT;
            }
            else {
                // use Newton
                dT_old = dT;
                dT = f/df;
                T -= dT;
            }

            // check for convergence
            version(complex_numbers) {
                if ((fabs(dT) < tol) && fabs(dT.im) < IMAGINARY_TOL) {
                    converged = true;
                    break;
                }
            }
            else {
                if (fabs(dT) < tol) {
                    converged = true;
                    break;
                }
            }

            // evaluate function for next iteration
            f = zero_func(T);
            df = _mixture_Cv_in_mode_at_temp(gs, T, imode);

            // maintain the brackets on the root
            if (f < 0.0) {
                Ta = T;
            }
            else {
                Tb = T;
            }
        }

        if (!converged) {
            string msg = "MultiTemperatureGas: Modal temperature failed to converge.\n";
            debug {
                msg ~= format("mode: %d\n", imode);
                msg ~= format("The final value was: %.16f\n", T);
                msg ~= "The supplied GasState was:\n";
                msg ~= gs.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }
        return T;
    }
}

version(multi_temperature_gas_test) {
    int main() {
        import util.msg_service;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);

        auto L = init_lua_State();
        doLuaFile(L, "sample-data/air-5sp-5T-gas-model.lua");
        string[] speciesNames;
        string[] energy_mode_names;
        getArrayOfStrings(L, "species", speciesNames);
        getArrayOfStrings(L, "energy_modes", energy_mode_names);
        auto tm = new MultiTemperatureGasMixture(L, speciesNames, energy_mode_names);
        lua_close(L);
        auto gs = GasState(5, 4);

        gs.p = 1.0e6;
        gs.T = 300.0;
        gs.T_modes = [to!number(901), to!number(302), to!number(403), to!number(704)];
        gs.massf = [to!number(0.8), to!number(0.1), to!number(0.025), to!number(0.025), to!number(0.05)];
        tm.updateFromPT(gs);
        number rho = gs.rho;
        gs.T = 450.0;
        gs.T_modes = [to!number(500), to!number(900), to!number(321), to!number(1500)];
        gs.p = 1.2e6;
        tm.updateFromRhoU(gs);

        assert(approxEqualNumbers(to!number(300.0), gs.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(901.0), gs.T_modes[0], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(302.0), gs.T_modes[1], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(403.0), gs.T_modes[2], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(704.0), gs.T_modes[3], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(1.0e6), gs.p, 1.0e-6), failedUnitTest());

        gs.T = 550.0;
        gs.T_modes = [to!number(500), to!number(900), to!number(321), to!number(700)];
        gs.rho = 0.1;
        tm.updateFromPU(gs);

        assert(approxEqualNumbers(to!number(300.0), gs.T, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(901.0), gs.T_modes[0], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(302.0), gs.T_modes[1], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(403.0), gs.T_modes[2], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(704.0), gs.T_modes[3], 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(rho), gs.rho, 1.0e-6), failedUnitTest());



        return 0;
    }
}

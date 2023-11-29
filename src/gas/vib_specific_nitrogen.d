/**
 * vib_specific_nitrogen.d
 * Authors: Rowan G., Katrina Sklavos and Peter J.
 *
 * This is a vibrationally-specific model for nitrogen
 * as descrbied in:
 *
 * D.Giordano, V.Bellucci, G.Colonna, M.Capitelli, I.Armenise and C.Bruno (1997)
 * Vibrationally Relaxing Flow of N2 past an Infinite Cylinder.
 * Journal of Thermophysics and Heat Transfer, vol 11, no 1, pp 27--35.
 *
 * The vibrational energies are associated with the thermochemical species.
 * We will have numVibLevels of these, starting with the ground state i=0.
 * Note that the vibrational quantum numbers start at 1 in the 1997 paper.
 */

module gas.vib_specific_nitrogen;

import std.algorithm.iteration;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;

import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;


class VibSpecificNitrogen: GasModel {
public:
    double _R_N2 = 296.805; // gas constant for N2
    double _M_N2 = 0.0280134; // kg/mole
    double _gamma = 7.0/5.0; // ratio of specific heats for vibrationally-frozen gas.
    double kB = Boltzmann_constant;
    int numVibLevels = 10; // default value
    double[] _vib_energy; // quantum level energies
    double[] _vib_energy_perkg;

    this(lua_State *L)
    {
        type_str = "VibSpecificNitrogen";
        // Bring table to TOS
        lua_getglobal(L, "VibSpecificNitrogen");
        if (lua_istable(L, -1)) {
            numVibLevels = getInt(L, -1, "numVibLevels");
        } else {
            // If there is not a table leave the default value for numVibLevels.
        }
        // In this model the nitrogen molecules with differing vibrational levels
        // are considered pseudo-chemical-species.
        _n_species = numVibLevels;
        _n_modes = 0;
        _species_names.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _species_names[isp] = format("N2-vib-%d", isp);
        }
        create_species_reverse_lookup();
        // The energy levels of the individual vibration levels are constant, once computed.
        _vib_energy.length = numVibLevels;
        _vib_energy_perkg.length = numVibLevels;
        foreach (i; 0 .. numVibLevels) {
            _vib_energy[i] = vib_energy(i);
            _vib_energy_perkg[i]= (Avogadro_number/_M_N2) * _vib_energy[i];
        }
    } // end constructor

    override string toString() const
    {
        char[] repr;
        repr ~= format("VibSpecificNitrogen=(numVibLevels=%d)", numVibLevels);
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q) const
    {
        Q.rho = Q.p / (_R_N2*Q.T);
        // For full internal energy, start with trans-rotational mode
        // and then add vibrational modes.
        Q.u = 2.5 * _R_N2 * Q.T;
        foreach (i; 0 .. numVibLevels) { Q.u += Q.massf[i]*_vib_energy_perkg[i]; }
    }
    override void update_thermo_from_rhou(ref GasState Q) const
    {
        // From internal energy, remove vibrational energy before computing trans-rotational temperature.
        number u = Q.u;
        foreach (i; 0 .. numVibLevels) { u -= Q.massf[i]*_vib_energy_perkg[i]; }
        Q.T = (0.4/_R_N2) * u;
        Q.p = Q.rho * _R_N2 * Q.T;
    }
    override void update_thermo_from_rhoT(ref GasState Q) const
    {
        Q.p = Q.rho * _R_N2 * Q.T;
        // Start with trans-rotational component of internal energy and add vibrational energy.
        Q.u = 2.5 * _R_N2 * Q.T;
        foreach (i; 0 .. numVibLevels) { Q.u += Q.massf[i]*_vib_energy_perkg[i]; }
    }
    override void update_thermo_from_rhop(ref GasState Q) const
    {
        Q.T = Q.p / (Q.rho * _R_N2);
        // Start with trans-rotational component of internal energy and add vibrational energy.
        Q.u = 2.5 * _R_N2 * Q.T;
        foreach (i; 0 .. numVibLevels) { Q.u += Q.massf[i]*_vib_energy_perkg[i]; }
    }
    override void update_thermo_from_ps(ref GasState Q, number s) const
    {
        throw new Error("VibSpecificNitrogen:update_thermo_from_ps NOT IMPLEMENTED.");
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s) const
    {
        throw new Error("VibSpecificNitrogen:update_thermo_from_hs NOT IMPLEMENTED.");
    }
    override void update_sound_speed(ref GasState Q) const
    {
        Q.a = sqrt(_gamma * _R_N2 * Q.T);
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        // The gas is inviscid.
        // [TODO] Add nitrogen model.
        Q.mu = 0.0;
        Q.k = 0.0;
    }
    override number dudT_const_v(in GasState Q) const
    {
        return to!number(_R_N2/(_gamma - 1.0));
    }
    override number dhdT_const_p(in GasState Q) const
    {
        throw new Error("ERROR: VibSpecificNitrogen:dhdT_const_p NOT IMPLEMENTED.");
    }
    override number dpdrho_const_T(in GasState Q) const
    {
        throw new Error("ERROR: VibSpecificNitrogen:dpdrho_const_T NOT IMPLEMENTED.");
    }
    override number gas_constant(in GasState Q) const
    {
        return to!number(_R_N2);
    }
    override number internal_energy(in GasState Q) const
    {
        return Q.u;
    }
    override number enthalpy(in GasState Q) const
    {
        return Q.u + Q.p/Q.rho;
    }
    override number entropy(in GasState Q) const
    {
        throw new Error("ERROR: VibSpecificNitrogen:entropy NOT IMPLEMENTED.");
    }

    double vib_energy(int i) const
    // Returns the quantum level energy for species i.
    {
        // The ground-state i=0 has a positive value for energy.
        // Note that the vibrational quantum numbers start at 1 in the 1997 paper.
        double w_e = 235857.0;
        double we_xe = 1432.4;
        double we_ye = -0.226;
        double h = 6.626e-34;
        double c = 2.998e8;
        double e = h*c * (w_e*(i+0.5) - we_xe*(i+0.5)^^2 + we_ye*(i+0.5)^^3);
        return e;
    }

    // Keep the following function as public
    // because the postprocessor will use it.
    number compute_Tvib(ref GasState Q, number Tguess, double tol=1.0e-9) const
    {
        // Use secant method to compute T.
        number x0 = Tguess;
        number x1 = x0 + 100.0;
        number fx0 = boltzmann_eq_species(0,x0) - Q.massf[0];
        number fx1 = boltzmann_eq_species(0,x1) - Q.massf[0];
        int max_it = 100;
        foreach(n; 0 .. max_it) {
            if (abs(fx1) < tol) {return x1;}
            number x2 = ((x0*fx1) - (x1*fx0)) / (fx1 - fx0);
            x0 = x1;
            x1 = x2;
            fx0 = fx1;
            fx1 = boltzmann_eq_species(0,x2) - Q.massf[0];
        } //end foreach
        return x1;
    } // end compute_Tvib()

    number boltzmann_eq_species(int i, number T) const
    // Returns the equilibrium population fraction for quantum-level i, given temperature.
    // i==0 is the ground state
    {
        number summ = 0;
        foreach(ej; 0 .. numVibLevels) { summ += exp(-_vib_energy[ej] / (kB*T)); }
        return (i < numVibLevels) ? exp(-_vib_energy[i]/(kB*T)) / summ : to!number(0.0);
    }
} // end class VibSpecificNitrogen

version(vib_specific_nitrogen_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/vib-specific-N2-gas.lua");
        auto gm = new VibSpecificNitrogen(L);
        lua_close(L);
        auto Q = GasState(gm.numVibLevels, 0);
        Q.p = 1.0e5;
        Q.T = 300.0;
        // Set up the species mass fractions assuming equilibrium.
        foreach (i; 0 .. gm.numVibLevels) { Q.massf[i] = gm.boltzmann_eq_species(i, Q.T); }

        double R_N2 = 296.805; // gas constant for N2
        double M_N2 = 0.0280134; // kg/mole
        double gamma = 7.0/5.0; // ratio of specific heats.

        gm.update_thermo_from_pT(Q);
        double my_rho = 1.0e5 / (R_N2 * 300.0);
        assert(isClose(Q.rho, my_rho, 1.0e-6), failedUnitTest());

        double my_u = 2.5 * R_N2 * 300.0;
        foreach (i; 0 .. gm.numVibLevels) {
            my_u += (Avogadro_number/M_N2) * gm.vib_energy(i) * Q.massf[i];
        }
        assert(isClose(Q.u, my_u, 1.0e-6), failedUnitTest());

        double Tvib = gm.compute_Tvib(Q, 400.0);
        assert(isClose(Tvib, Q.T, 1.0e-3));

        gm.update_trans_coeffs(Q);
        assert(isClose(Q.mu, 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(Q.k, 0.0, 1.0e-6), failedUnitTest());

        gm.update_sound_speed(Q);
        double my_a = sqrt(gamma * R_N2 * 300.0);
        assert(isClose(Q.a, my_a, 1.0e-6), failedUnitTest());

        return 0;
    }
}


/**
 * composite_gas.d
 * A gas model with behaviour built by composition of pieces.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2021-02-09
 */

module gas.composite_gas;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;


class CompositeGas : GasModel {
public:

    this(lua_State *L)
    {
        type_str = "CompositeGas";
        /* There are some top-level GasModel services that require us to fill in data
         * at this point.
         */
        getArrayOfStrings(L, LUA_GLOBALSINDEX, "species", _species_names);
        _n_species = cast(uint) _species_names.length;
        create_species_reverse_lookup();
        _n_modes = 0; // Single temperature gas
        if (canFind(_species_names, "e-") || canFind(_species_names, "eminus")) {
            if (!(_species_names[$-1] == "e-" || _species_names[$-1] == "eminus")) {
                throw new Error("Electrons should be last species.");
            }
            _is_plasma = true;
        }

        foreach (isp; 0 .. _n_species) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, _species_names[isp].toStringz);
            _mol_masses[isp] = getDouble(L, -1, "M");
            lua_pop(L, 1);
            lua_pop(L, 1);
        }

        ThermodynamicModel mThermo = new ThermodynamicModel(L, _species_names, _mol_masses);
        TransportPropertiesModel mTransProps = new TransportPropertiesModel(L, mThermo);
        
    }
    
    // Service methods related to thermodynamics
    // Updates to GasState
    @nogc override void update_thermo_from_pT(GasState gs)
    {
        mThermo.updateFromPT(gs);
    }
    @nogc override void update_thermo_from_rhou(GasState gs)
    {
        mThermo.updateFromRhoU(gs);
    }
    @nogc override void update_thermo_from_rhoT(GasState gs)
    {
        mThermo.updateFromRhoT(GasState gs);
    }
    @nogc override void update_thermo_from_ps(GasState gs, number s)
    {
        mThermo.updateFromPS(gs, s);
    }
    @nogc override void update_thermo_from_hs(GasState gs, number s)
    {
        mThermo.updateFromHS(gs, s);
    }
    @nogc override void update_sound_speed(GasState gs)
    {
        mThermo.updateSoundSpeed(gs);
    }
    // Return of single values
    @nogc override number dudT_const_v(in GasState gs)
    {
        retrun mThermo.dudT_const_v(gs);
    }
    @nogc override number dhdT_const_p(in GasState gs)
    {
        return mThermo.dhdT_const_p(gs);
    }
    @nogc override number dpdrho_const_T(in GasState gs)
    {
        return mThermo.dpdrho_const_T(gs);
    }
    @nogc override number gas_constant(in GasState gs)
    {
        return mThermo.gas_constant(gs);
    }
    @nogc override number internal_energy(in GasState gs)
    {
        return mThermo.internal_energy(gs);
    }
    @nogc override number enthalpy(in GasState gs)
    {
        return mThermo.enthalpy(gs);
    }
    @nogc override number enthalpy(in GasState gs, int isp)
    {
        return mThermo.enthalpy(gs, isp);
    }
    @nogc override number enthalpyPerSpeciesInMode(in GasState gs, int isp, int imode)
    {
        return mThermo.enthalpy(gs, isp, imode);
    }
    @nogc override number entropy(in GasState gs)
    {
        return mThermo.entropy(gs);
    }
    @nogc override number entropy(in GasState gs, int isp)
    {
        return mThermo.entropy(gs, isp);
    }
    @nogc override number gibbs_free_energy(GasState gs, int isp)
    {
        return mThermo.gibbs_free_energy(gs, isp);
    }

    @nogc override void update_trans_coeffs(GasState gs)
    {
        mTransProps.updateTransCoeffs(gs);
    }

private:
    ThermodynamicModel mThermo;
    TransportPropertiesModel mTransProps;
}

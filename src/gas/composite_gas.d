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
    

private:
    ThermodynamicModel mThermo;
    TransportPropertiesModel mTransProps;
}

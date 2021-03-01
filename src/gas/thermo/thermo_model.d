/*
 * thermo_model.d
 * Interface for all models which describe thermodynamic behaviour.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2021-02-09
 */

module gas.thermo.thermo_model;

import nm.number;
import gas : GasState;


interface ThermodynamicModel {
public:
    // Methods related to computing thermo state.
    @nogc void update_thermo_from_pT(GasState gs);
    @nogc void update_thermo_from_rhou(GasState gs);
    @nogc void update_thermo_from_rhoT(GasState gs);
    @nogc void update_thermo_from_rhop(GasState gs);
    @nogc void update_thermo_from_ps(GasState gs, number s);
    @nogc void update_thermo_from_hs(GasState gs, number h, number s);
    @nogc void update_sound_speed(GasState gs);
    // Methods related to computing thermo derivatives.
    @nogc number dudT_const_v(in GasState gs);
    @nogc number dhdT_const_p(in GasState gs);
    @nogc number dpdrho_const_T(in GasState gs);
    // Methods related to single properties of the mixture.
    @nogc number gas_constant(in GasState gs);
    @nogc number internal_energy(in GasState gs); 
    @nogc number enthalpy(in GasState gs);
    @nogc number enthalpy(in GasState gs, int isp);
    @nogc number enthalpy(in GasState gs, int isp, int imode);
    @nogc number entropy(in GasState gs);
    @nogc number entropy(in GasState gs, int isp);
}

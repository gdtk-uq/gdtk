/*
 * thermo_model.d
 * Interface for all models which describe thermodynamic behaviour.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2021-02-09
 */

module gas.thermo.thermo_model;

import nm.number;
import gas.gas_state : GasState;


interface ThermodynamicModel {
public:
    void updateFromPU(ref GasState gs);
    // Methods related to computing thermo state.
    @nogc void updateFromPT(ref GasState gs);
    @nogc void updateFromRhoU(ref GasState gs);
    @nogc void updateFromRhoT(ref GasState gs);
    @nogc void updateFromRhoP(ref GasState gs);
    @nogc void updateFromPS(ref GasState gs, number s);
    @nogc void updateFromHS(ref GasState gs, number h, number s);
    @nogc void updateSoundSpeed(ref GasState gs);
    // Methods related to computing thermo derivatives.
    @nogc number dudTConstV(in GasState gs);
    @nogc number dhdTConstP(in GasState gs);
    @nogc number dpdrhoConstT(in GasState gs);
    // Methods related to single properties of the mixture.
    @nogc number gasConstant(in GasState gs);
    @nogc number gasConstant(in GasState gs, int isp);
    @nogc number internalEnergy(in GasState gs);
    @nogc number energyPerSpeciesInMode(in GasState gs, int isp, int imode);
    @nogc number enthalpy(in GasState gs);
    @nogc number enthalpyPerSpecies(in GasState gs, int isp);
    @nogc number enthalpyPerSpeciesInMode(in GasState gs, int isp, int imode);
    @nogc number entropy(in GasState gs);
    @nogc number entropyPerSpecies(in GasState gs, int isp);
    @nogc number cpPerSpecies(in GasState gs, int isp);
    @nogc void enthalpies(in GasState gs, number[] hs);
    @nogc void GibbsFreeEnergies(in GasState gs, number[] gibbs_energies);
}

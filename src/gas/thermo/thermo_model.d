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
    @nogc void updateFromPT(GasState gs);
    @nogc void updateFromRhoU(GasState gs);
    @nogc void updateFromRhoT(GasState gs);
    @nogc void updateFromRhoP(GasState gs);
    @nogc void updateFromPS(GasState gs, number s);
    @nogc void updateFromHS(GasState gs, number h, number s);
    @nogc void updateSoundSpeed(GasState gs);
    // Methods related to computing thermo derivatives.
    @nogc number dudTConstV(in GasState gs);
    @nogc number dhdTConstP(in GasState gs);
    @nogc number dpdrhoConstT(in GasState gs);
    // Methods related to single properties of the mixture.
    @nogc number gasConstant(in GasState gs);
    @nogc number internalEnergy(in GasState gs); 
    @nogc number enthalpy(in GasState gs);
    @nogc number enthalpy(in GasState gs, int isp);
    @nogc number enthalpy(in GasState gs, int isp, int imode);
    @nogc number entropy(in GasState gs);
    @nogc number entropy(in GasState gs, int isp);
}

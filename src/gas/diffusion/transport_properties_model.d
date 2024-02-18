/**
 * transport_properties_model.d
 * Interface for all models that compute viscosity and thermal conductivities.
 *
 * Author: Rowan G.
 * Date: 2021-03-01
 **/

module gas.diffusion.transport_properties_model;

import nm.number;
import gas.gas_state : GasState;
import gas.gas_model : GasModel;

interface TransportPropertiesModel {
public:
    @nogc void updateTransProps(ref GasState gs, GasModel gm);
    @nogc void binaryDiffusionCoefficients(ref const(GasState) gs, ref number[][] D);
}


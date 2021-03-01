/*
 * transport_properties_model.d
 * Interface for all models that compute viscosity and thermal conductivities.
 *
 * Author: Rowan G.
 * Version: 2021-03-01
 */

module gas.diffusion.transport_properties_model;

import gas : GasState;

interface TransportPropertiesModel {
public:
    @nogc void updateTransProps(GasState gs);
}


/*
 * viscosity.d
 * Interface for all viscosity models.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-08-19 -- initial cut
 */

module gas.diffusion.viscosity;

import gas.gas_model;

interface Viscosity {
    Viscosity dup() const;
    final void update_viscosity(GasState Q)
    {
	Q.mu = eval(Q);
    }
    double eval(in GasState Q);
}

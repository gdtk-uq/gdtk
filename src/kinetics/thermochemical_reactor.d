/**
 * thermochemical_reactor.d
 * Base class for all of the kinetics package updaters.
 *
 * Author: Peter J and Rowan G.
 */

module kinetics.thermochemical_reactor;

import nm.complex;
import nm.number;
import gas;

class ThermochemicalReactorUpdateException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

// The extra parameters are intended for some special gas and kinetics models,
// such as JJ Hoste's mixing-limited fuel-air model, that need context information
// from the flow solver.
// Most gas models and reaction schemes will just ignore the params array.
immutable size_t maxParams = 10;

class ThermochemicalReactor {
public:
    this(GasModel gmodel)
    {
        // We need a reference to the original gas model object
        // to update the GasState data at a later time.
        _gmodel = gmodel;
    }

    // All the work happens when calling the concrete object
    // which updates the GasState over the (small) time, tInterval.
    //
    // The array params is there to allow extra information to be passed in.
    // For example, the mixing-limited combustion model by JJ Hoste needs
    // some information about the local flow state beyond the usual gas state.
    @nogc
    abstract void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[maxParams] params);

    // We will need to access this referenced model from the Lua functions
    // so it needs to be public.
    GasModel _gmodel;
} // end class ThermochemicalReactor

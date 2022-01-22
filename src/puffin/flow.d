// flow.d -- Part of the Puffin steady-flow calculator.
//
// PA Jacobs
// 2022-01-22
//
module flow;

import std.format;

import geom;
import gas;


class FlowState2D {
    // FlowSate2D objects are passed to our interpolation functions and flux calculator.
public:
    GasState gas;
    Vector3 vel;

    this(GasModel gm)
    {
        gas = new GasState(gm);
        vel.set(0.0, 0.0);
    }

    this(ref const(FlowState2D) other)
    {
        gas = new GasState(other.gas);
        vel.set(other.vel);
    }

    override
    string toString() const
    {
        string repr = "FlowState2D(";
        repr ~= format("gas=%s, vel=%s", gas, vel);
        repr ~= ")";
        return repr;
    }
} // end class FlowState


class CQIndex {
    // Indices into the vector of conserved quantities.
public:
    size_t n;
    size_t n_species;
    size_t n_modes;

    size_t mass;
    size_t xMom;
    size_t yMom;
    size_t totEnergy;
    size_t species;
    size_t modes;

    this(size_t nspecies, size_t nmodes)
    {
        n = 4;
        mass = 0;
        xMom = 1;
        yMom = 2;
        totEnergy = 3;
        //
        species = n; // Start of species elements.
        n_species = nspecies;
        n += nspecies;
        //
        modes = n; // Start of modes elements.
        n_modes = nmodes;
        n += nmodes;
    }

    this(ref const(CQIndex) other)
    {
        n = other.n;
        //
        mass = other.mass;
        xMom = other.xMom;
        yMom = other.yMom;
        totEnergy = other.totEnergy;
        //
        species = other.species;
        n_species = other.n_species;
        //
        modes = other.modes;
        n_modes = other.n_modes;
    }

    override
    string toString() const
    {
        string repr = "CQIndex(";
        repr ~= format("n=%d, n_species=%d, n_modes=%d", n, n_species, n_modes);
        repr ~= format(", mass=%d, xMom=%d, yMom=%d, totEnergy=%d",
                       mass, xMom, yMom, totEnergy);
        repr ~= format(", species=%d, modes=%d)", species, modes);
        return repr;
    }
} // end class CQIndex

/**
 * conservedquantities.d
 * Class for the vector of conserved quantities, for use in the CFD codes.
 *
 * Author: Peter J., Rowan G. and Kyle Damm
 * Version:
 * 2014-07-17: initial cut, to explore options.
 * 2021-05-10: Change to array storage.
 */

module conservedquantities;

import std.string;
import std.conv;
import nm.complex;
import nm.number;
import geom;
import gas;


// Underlying definition of the conserved quantities collection,
// as seen by the transient solver.

class ConservedQuantities {
public:
    number[] vec;

    this(size_t n)
    {
        vec.length = n;
    }

    this(ref const(ConservedQuantities) other)
    {
        vec.length = other.vec.length;
        foreach (i, ref e; vec) { e = other.vec[i]; }
    }

    @nogc void copy_values_from(ref const(ConservedQuantities) src)
    {
        foreach (i, ref e; vec) { e = other.vec[i]; }
    }

    @nogc void clear()
    {
        foreach (ref e; vec) { e = 0.0; }
    }

    @nogc void add(ref const(ConservedQuantities) other, double factor=1.0)
    {
        foreach (i, ref e; vec) e += other.vec[i] * factor;
    }

    @nogc void scale(double factor)
    {
        foreach (ref e; vec) { e *= factor; }
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "ConservedQuantities(vec=" ~ to!string(vec) ~ ")";
        return to!string(repr);
    }

    version(complex_numbers) {
        // When performing the complex-step Frechet derivative in the Newton-Krylov accelerator,
        // the conserved quantities accumulate imaginary components,
        // so we have to start with a clean slate, so to speak.
        @nogc void clear_imaginary_components()
        {
            foreach (ref e; vec) { e.im = 0.0; }
        }
    } // end version(complex_numbers)


    /+ Retain the old named items code, just for reference as we make the changes.
    number mass;           // density, kg/m**3
    Vector3 momentum;      // momentum/unit volume
    number total_energy;   // total energy
    version(multi_species_gas) {
        number[] massf;    // mass fractions of species
    }
    version(multi_T_gas) {
        number[] energies; // modal energies
    }
    version(MHD) {
        Vector3 B;         // magnetic field, Tesla
        number psi;        // divergence cleaning parameter for MHD
        number divB;       // divergence of the magnetic field
    }
    version(turbulence) {
        number[2] rhoturb;    // turbulent conserved
    }
    +/
} // end class ConservedQuantities


// Now that the ConservedQuantities object is a simple vector of quantities,
// this collection of indices helps us select individual elements.

class ConservedQuantitiesIndices {
public:
    bool threeD;
    size_t n;
    size_t n_species;
    size_t n_modes;
    size_t n_turb;
    size_t mass;
    size_t xMom;
    size_t yMom;
    size_t zMom;
    size_t totEnergy;
    size_t rhoturb;
    size_t xB;
    size_t yB;
    size_t zB;
    size_t psi;
    size_t divB;
    size_t species;
    size_t modes;

    this(int dimensions, size_t nturb, bool MHD, size_t nspecies, size_t nmodes) {
        mass = 0;
        xMom = 1;
        yMom = 2;
        if (dimensions == 3) {
            threeD = true;
            zMom = 3;
            totEnergy = 4;
            n = 5;
        } else {
            // Do not carry z-momentum for 2D simulations.
            threeD = false;
            totEnergy = 3;
            n = 4;
        }
        n_turb = nturb;
        if (nturb > 0) {
            rhoturb = n;
            n += nturb;
        }
        if (MHD) {
            xB = n;
            yB = n+1;
            zB = n+2;
            psi = n+3;
            divB = n+4;
            n += 5;
        }
        n_species = nspecies;
        if (nspecies > 1) {
            species = n;
            n += nspecies;
            // Note that we only carry species in the conserved-quantities vector
            // if we have a multi-species gas model.
            // A single-species gas model assumes a species fraction on 1.0
            // throughout the flow solver code.
        }
        n_modes = nmodes;
        if (nmodes > 0) {
            modes = n;
            n += nmodes;
        }
    } // end constructor

    this(const(ConservedQuantitiesIndices) other)
    {
        threeD = other.threeD;
        n = other.n;
        n_species = other.n_species;
        n_modes = other.n_modes;
        mass = other.mass;
        xMom = other.xMom;
        yMom = other.yMom;
        zMom = other.zMom;
        totEnergy = other.totEnergy;
        rhoturb = other.rhoturb;
        xB = other.xB;
        yB = other.yB;
        zB = other.zB;
        psi = other.psi;
        divB = other.divB;
        species = other.species;
        modes = other.modes;
    } // end copy constructor

} // end ConvservedQuantitiesIndices

/**
 * conservedquantities.d
 * Class for the vector of conserved quantities, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-17: initial cut, to explore options.
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

    this(int n_species, int n_modes)
    {
        version(multi_species_gas) {
            massf.length = n_species;
        } else {
            assert(n_species == 1, "only single-specied gas available");
        }
        version(multi_T_gas) {
            energies.length = n_modes;
        } else {
            assert(n_modes == 0, "single-temperature gas only");
        }
    }

    this(ref const(ConservedQuantities) other)
    {
        mass = other.mass;
        momentum = other.momentum;
        total_energy = other.total_energy;
        version(multi_species_gas) {
            massf = other.massf.dup;
        }
        version(multi_T_gas) {
            energies = other.energies.dup;
        }
        version(MHD) {
            B = other.B;
            psi = other.psi;
            divB = other.divB;
        }
        version(turbulence) {
            rhoturb = other.rhoturb.dup;
        }
    }

    @nogc
    void copy_values_from(ref const(ConservedQuantities) src)
    {
        mass = src.mass;
        momentum.set(src.momentum);
        total_energy = src.total_energy;
        version(multi_species_gas) {
            massf[] = src.massf[];
        }
        version(multi_T_gas) {
            energies[] = src.energies[];
        }
        version(MHD) {
            B.set(src.B);
            psi = src.psi;
            divB = src.divB;
        }
        version(turbulence) {
            rhoturb[] = src.rhoturb[];
        }
    }

    @nogc
    void clear()
    {
        mass = 0.0;
        momentum.clear();
        total_energy = 0.0;
        version(multi_species_gas) {
            foreach(ref mf; massf) { mf = 0.0; }
        }
        version(multi_T_gas) {
            foreach(ref e; energies) { e = 0.0; }
        }
        version(MHD) {
            B.clear();
            psi = 0.0;
            divB = 0.0;
        }
        version(turbulence) {
            foreach(ref rt; rhoturb) { rt = 0.0; }
        }
    }

    @nogc
    void add(ref const(ConservedQuantities) other, double factor=1.0)
    {
        mass += other.mass * factor;
        momentum.add(other.momentum, factor);
        total_energy += other.total_energy * factor;
        version(multi_species_gas) {
            foreach(i; 0 .. massf.length) { massf[i] += other.massf[i] * factor; }
        }
        version(multi_T_gas) {
            foreach(i; 0 .. energies.length) { energies[i] += other.energies[i] * factor; }
        }
        version(MHD) {
            B.add(other.B, factor);
            psi += other.psi * factor;
            divB += other.divB * factor;
        }
        version(turbulence) {
            foreach(i; 0 .. rhoturb.length) { rhoturb[i] += other.rhoturb[i] * factor; }
        }
    }

    @nogc
    void scale(double factor)
    {
        mass *= factor;
        momentum.scale(factor);
        total_energy *= factor;
        version(multi_species_gas) {
            foreach(ref mf; massf) { mf *= factor; }
        }
        version(multi_T_gas) {
            foreach(ref e; energies) { e *= factor; }
        }
        version(MHD) {
            B.scale(factor);
            psi *= factor;
            divB *= factor;
        }
        version(turbulence) {
            foreach(ref rt; rhoturb) { rt *= factor; }
        }
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "ConservedQuantities(";
        repr ~= "mass=" ~ to!string(mass);
        repr ~= ", momentum=" ~ to!string(momentum);
        repr ~= ", total_energy=" ~ to!string(total_energy);
        version(multi_species_gas) {
            repr ~= ", massf=" ~ to!string(massf);
        }
        version(multi_T_gas) {
            repr ~= ", energies=" ~ to!string(energies);
        }
        version(MHD) {
            repr ~= ", B=" ~ to!string(B);
            repr ~= ", psi=" ~ to!string(psi);
            repr ~= ", divB-" ~ to!string(divB);
        }
        version(turbulence) {
            repr ~= ", turb=" ~ to!string(rhoturb);
        }
        repr ~= ")";
        return to!string(repr);
    }

version(complex_numbers) {
    @nogc
    void clear_imaginary_components()
    // When performing the complex-step Frechet derivative in the Newton-Krylov accelerator,
    // the conserved quantities accumulate imaginary components, so we have to start with a clean slate, so to speak.
    {
        mass.im = 0.0;
        momentum.refx.im = 0.0;
        momentum.refy.im = 0.0;
        momentum.refz.im = 0.0;
        total_energy.im = 0.0;
        version(multi_species_gas) {
            foreach(ref mf; massf) { mf.im = 0.0; }
        }
        version(multi_T_gas) {
            foreach(ref e; energies) { e.im = 0.0; }
        }
        version(MHD) {
            B.refx.im = 0.0;
            B.refy.im = 0.0;
            B.refz.im = 0.0;
            psi.im = 0.0;
            divB.im = 0.0;
        }
        version(turbulence) {
            foreach(ref rt; rhoturb) { rt.im = 0.0; }
        }
    } // end clear_imaginary_components()
} // end version(complex)

} // end class ConservedQuantities


// When formulating the sensitivity matrices for the steady-state solvers,
// it is convenient to be able to look at the ConservedQuantities object
// as a vector of quantities that can be simply indexed.
// This collections of indices helps us do that.
// Note that this view of the ConservedQuantities vector is not quite the
// same as that seen by the explicit updates in the transient code.

struct ConservedQuantitiesIndices {
    size_t nConservedQuantities;
    size_t mass;
    size_t xMom;
    size_t yMom;
    size_t zMom;
    size_t totEnergy;
    size_t tke;
    size_t species;

    this(int dimensions, size_t nturb, size_t nmodes, size_t nspecies) {
        mass = 0;
        xMom = 1;
        yMom = 2;
        if ( dimensions == 2 ) {
            totEnergy = 3;
            nConservedQuantities = 4;
        }
        else { // 3D simulations
            zMom = 3;
            totEnergy = 4;
            nConservedQuantities = 5;
        }
        if ( nturb > 0) {
            tke = nConservedQuantities;
            nConservedQuantities += nturb;
        }
        if ( nspecies > 1) {
            species = nConservedQuantities;
            nConservedQuantities += nspecies;
        }
    }

    number get_mass(number[] cqv) { return cqv[0]; }
    void set_mass(number[] cqv, number value) { cqv[0] = value; }
    // [TODO] PJ 2021-05-08 Should we fill in more accessor functions so that we can
    // move to using simple arrays for the Vectors of conserved quantities?
    // Using simple arrays is likely to simplify much of the gasdynamic update code
    // and we are likely only to want to use these functions at the encode and
    // decode functions for the conserved quantities.
    // I think that most of the code in the ConservedQuantities class above
    // is likely to be eliminated.

} // end ConvservedQuantitiesIndices

/**
 * conservedquantities.d
 * Storage for the vector of conserved quantities, for use in the CFD codes.
 *
 * Author: Peter J., Rowan G. and Kyle Damm
 * Version:
 * 2014-07-17: initial cut, to explore options.
 * 2021-05-10: Change to array storage.
 * 2022-08-20: Bare array.
 */

module conservedquantities;

import std.string;
import std.format;
import std.conv;
import nm.complex;
import nm.number;
import geom;
import gas;


// Underlying definition of the conserved quantities collection,
// as seen by the transient solver, is just an array.
immutable alias ConservedQuantities = number[];

ConservedQuantities new_ConservedQuantities(size_t n)
{
    return new number[n];
}

@nogc void copy_values_from(ref number[] vec, ref const(number[]) src)
{
    foreach (i, ref e; vec) { e = src[i]; }
}

@nogc void clear(ref number[] vec)
{
    foreach (ref e; vec) { e = 0.0; }
}

@nogc void add(ref number[] vec, ref const(number[]) other, double factor=1.0)
{
    foreach (i, ref e; vec) e += other[i] * factor;
}

@nogc void scale(ref number[] vec, double factor)
{
    foreach (ref e; vec) { e *= factor; }
}


version(complex_numbers) {
    // When performing the complex-step Frechet derivative in the Newton-Krylov accelerator,
    // the conserved quantities accumulate imaginary components,
    // so we have to start with a clean slate, so to speak.
    @nogc void clear_imaginary_components(ref number[] vec)
    {
        foreach (ref e; vec) { e.im = 0.0; }
    }
} // end version(complex_numbers)


// Now that the ConservedQuantities object is a simple vector of quantities,
// this collection of indices helps us select individual elements.

struct ConservedQuantitiesIndices {
public:
    bool threeD;
    bool turb;
    bool MHD;
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

        bool put_mass_in_last_position = false;
        version (nk_accelerator) {
            // we will drop the mass continuity equation if we are running a multi-species calculation with
            // the steady-state solver, note that we still need an entry in the conserved quantities vector
            // for the mass (some parts of the code expect it), so we will place it in the last position
            if (nspecies > 1) { put_mass_in_last_position = true; }
        }

        if (put_mass_in_last_position) {
            xMom = 0;
            yMom = 1;
            if (dimensions == 3) {
                threeD = true;
                zMom = 2;
                totEnergy = 3;
                n = 4;
            } else {
                // Do not carry z-momentum for 2D simulations.
                threeD = false;
                totEnergy = 2;
                n = 3;
            }
            n_turb = nturb;
            if (nturb > 0) {
                turb = true;
                rhoturb = n; // Start of turbulence elements.
                n += nturb;
            } else {
                turb = false;
            }
            this.MHD = MHD;
            if (MHD) {
                xB = n;
                yB = n+1;
                zB = n+2;
                psi = n+3;
                divB = n+4;
                n += 5;
            }
            n_species = nspecies;
            species = n; // Start of species elements.
            n += nspecies;
            n_modes = nmodes;
            if (nmodes > 0) {
                modes = n; // Start of modes elements.
                n += nmodes;
            }
            // we still need the mass in the conserved quantities vector in some places of the code
            mass = n;
            n += 1;
        } else {
            // fill out the array using our standard ordering (for the transient code)
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
                turb = true;
                rhoturb = n; // Start of turbulence elements.
                n += nturb;
            } else {
                turb = false;
            }
            this.MHD = MHD;
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
                species = n; // Start of species elements.
                n += nspecies;
                // Note that we only carry species in the conserved-quantities vector
                // if we have a multi-species gas model.
                // A single-species gas model assumes a species fraction on 1.0
                // throughout the flow solver code.
            }
            n_modes = nmodes;
            if (nmodes > 0) {
                modes = n; // Start of modes elements.
                n += nmodes;
            }
        }
    } // end constructor

    this(ref const(ConservedQuantitiesIndices) other)
    {
        threeD = other.threeD;
        turb = other.turb;
        MHD = other.MHD;
        n = other.n;
        n_turb = other.n_turb;
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

    string toString() const
    {
        char[] repr;
        repr ~= "ConservedQuantitiesIndices(";
        repr ~= format("threeD=%s, turb=%s, MHD=%s", threeD, turb, MHD);
        repr ~= format(", n=%d, n_turb=%d, n_species=%d, n_modes=%d", n, n_turb, n_species, n_modes);
        repr ~= format(", mass=%d, xMom=%d, yMom=%d, zMom=%d, totEnergy=%d", mass, xMom, yMom, zMom, totEnergy);
        repr ~= format(", rhoturb=%d, xB=%d, yB=%d, zB=%d, psi=%d, divB=%d", rhoturb, xB, yB, zB, psi, divB);
        repr ~= format(", species=%d, modes=%d", species, modes);
        repr ~= ")";
        return to!string(repr);
    }
} // end struct ConservedQuantitiesIndices

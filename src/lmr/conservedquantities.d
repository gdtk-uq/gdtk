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
import ntypes.complex;
import nm.number;
import globalconfig;
import geom;
import gas;
import turbulence;

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
    string[] names;

    this(int dimensions, size_t nturb, bool MHD, size_t nspecies, size_t nmodes, ref GasModel gmodel, ref TurbulenceModel turb_model) {

        bool drop_mass_from_calculation = false;
        version (newton_krylov) {
            // Note that for a multi-species simulation, the n species conservation
            // equations and the total mass conservation equation are linearly
            // dependent. Because of this, we deduct the number of conserved
            // quantities by 1 for multi-species simulations, and remove the total
            // mass conservation equation from the system of equations we are
            // solving. For a discussion on this and other possible approaches see
            // pg. 8 of Multicomponent Flow Modelling by V. Giovangigli (1999).
            // (Comment moved here by NNG 24/04/08)
            if (nspecies > 1) { drop_mass_from_calculation = true; }
        }

        if (drop_mass_from_calculation) {
            xMom = 0; names ~= "x-mom";
            yMom = 1; names ~= "y-mom";
            if (dimensions == 3 || MHD) {
                threeD = true;
                zMom = 2; names ~= "z-mom";
                totEnergy = 3; names ~= "total-energy";
                n = 4;
            } else {
                // Do not carry z-momentum for 2D simulations (unless it's MHD).
                threeD = false;
                totEnergy = 2; names ~= "total-energy";
                n = 3;
            }
            n_turb = nturb;
            if (nturb > 0) {
                turb = true;
                rhoturb = n; // Start of turbulence elements.
                n += nturb;
                foreach (i; 0 .. nturb) names ~= turb_model.primitive_variable_name(to!int(i));
            } else {
                turb = false;
            }
            this.MHD = MHD;
            if (MHD) {
                xB = n; names ~= "x-B";
                yB = n+1; names ~= "y-B";
                zB = n+2; names ~= "z-B";
                psi = n+3; names ~= "psi";
                divB = n+4; names ~= "div-B";
                n += 5;
            }
            n_species = nspecies;
            species = n; // Start of species elements.
            n += nspecies;
            foreach (i; 0 .. nspecies) names ~= gmodel.species_name(to!int(i)).toUpper;
            n_modes = nmodes;
            if (nmodes > 0) {
                modes = n; // Start of modes elements.
                n += nmodes;
                foreach (i; 0 .. nmodes) names ~= gmodel.energy_mode_name(to!int(i));
            }
            // Set mass to a nonzero value to indicate no mass equation being solved.
            mass = 999999;
        } else {
            // fill out the array using our standard ordering (for the transient code)
            mass = 0; names ~= "mass";
            xMom = 1; names ~= "x-mom";
            yMom = 2; names ~= "y-mom";
            if (dimensions == 3 || MHD) {
                threeD = true;
                zMom = 3; names ~= "z-mom";
                totEnergy = 4; names ~= "total-energy";
                n = 5;
            } else {
                // Do not carry z-momentum for 2D simulations (unless it's MHD).
                threeD = false;
                totEnergy = 3; names ~= "total-energy";
                n = 4;
            }
            n_turb = nturb;
            if (nturb > 0) {
                turb = true;
                rhoturb = n; // Start of turbulence elements.
                n += nturb;
                foreach (i; 0 .. nturb) names ~= turb_model.primitive_variable_name(to!int(i));
            } else {
                turb = false;
            }
            this.MHD = MHD;
            if (MHD) {
                xB = n; names ~= "x-B";
                yB = n+1; names ~= "y-B";
                zB = n+2; names ~= "z-B";
                psi = n+3; names ~= "psi";
                divB = n+4; names ~= "div-B";
                n += 5;
            }
            n_species = nspecies;
            if (nspecies > 1) {
                species = n; // Start of species elements.
                n += nspecies;
                foreach (i; 0 .. nspecies) names ~= gmodel.species_name(to!int(i)).toUpper;
                // Note that we only carry species in the conserved-quantities vector
                // if we have a multi-species gas model.
                // A single-species gas model assumes a species fraction on 1.0
                // throughout the flow solver code.
            }
            n_modes = nmodes;
            if (nmodes > 0) {
                modes = n; // Start of modes elements.
                n += nmodes;
                foreach (i; 0 .. nmodes) names ~= gmodel.energy_mode_name(to!int(i));
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
	names = other.names.dup;
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

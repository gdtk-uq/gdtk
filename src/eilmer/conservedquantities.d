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

class ConservedQuantities {
public:
    number mass;         // density, kg/m**3
    Vector3 momentum;    // momentum/unit volume
    Vector3 B;           // magnetic field, Tesla
    number total_energy; // total energy
    number[] massf;      // mass fractions of species
    number[] energies;   // modal energies (mode 0 is usually transrotational)
    number psi;          // divergence cleaning parameter for MHD
    number divB;         // divergence of the magnetic field
    number tke;          // turbulent kinetic energy
    number omega;        // omega from k-omega turbulence model

    this(int n_species, int n_modes)
    {
        massf.length = n_species;
        energies.length = n_modes;
    }

    this(ref const(ConservedQuantities) other)
    {
        mass = other.mass;
        momentum = other.momentum;
        B = other.B;
        total_energy = other.total_energy;
        massf = other.massf.dup;
        energies = other.energies.dup;
        psi = other.psi;
        divB = other.divB;
        tke = other.tke;
        omega = other.omega;
    }

    @nogc
    void copy_values_from(ref const(ConservedQuantities) src)
    {
        mass = src.mass;
        momentum.set(src.momentum);
        B.set(src.B);
        total_energy = src.total_energy;
        massf[] = src.massf[];
        energies[] = src.energies[];
        psi = src.psi;
        divB = src.divB;
        tke = src.tke;
        omega = src.omega;
    }

    @nogc
    void clear()
    {
        mass = 0.0;
        momentum.clear();
        B.clear();
        total_energy = 0.0;
        foreach(ref mf; massf) { mf = 0.0; }
        foreach(ref e; energies) { e = 0.0; }
        psi = 0.0;
        divB = 0.0;
        tke = 0.0;
        omega = 0.0;
    }

    @nogc
    void add(ref const(ConservedQuantities) other, double factor=1.0)
    {
        mass += other.mass * factor;
        momentum.add(other.momentum, factor);
        B.add(other.B, factor);
        total_energy += other.total_energy * factor;
        foreach(i; 0 .. massf.length) { massf[i] += other.massf[i] * factor; }
        foreach(i; 0 .. energies.length) { energies[i] += other.energies[i] * factor; }
        psi += other.psi * factor;
        divB += other.divB * factor;
        tke += other.tke * factor;
        omega += other.omega * factor;
    }

    @nogc
    void scale(double factor)
    {
        mass *= factor;
        momentum.scale(factor);
        B.scale(factor);
        total_energy *= factor;
        foreach(ref mf; massf) { mf *= factor; }
        foreach(ref e; energies) { e *= factor; }
        psi *= factor;
        divB *= factor;
        tke *= factor;
        omega *= factor;
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "ConservedQuantities(";
        repr ~= "mass=" ~ to!string(mass);
        repr ~= ", momentum=" ~ to!string(momentum);
        repr ~= ", B=" ~ to!string(B);
        repr ~= ", total_energy=" ~ to!string(total_energy);
        repr ~= ", massf=" ~ to!string(massf);
        repr ~= ", energies=" ~ to!string(energies);
        repr ~= ", psi=" ~ to!string(psi);
        repr ~= ", divB-" ~ to!string(divB);
        repr ~= ", tke=" ~ to!string(tke);
        repr ~= ", omega=" ~ to!string(omega);
        repr ~= ")";
        return to!string(repr);
    }
} // end class ConservedQuantities

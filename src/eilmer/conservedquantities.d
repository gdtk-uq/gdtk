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
    number mass;           // density, kg/m**3
    Vector3 momentum;      // momentum/unit volume
    number total_energy;   // total energy
    version(multi_species_gas) {
        number[] massf;    // mass fractions of species
    }
    version(multi_T_gas) {
        number[] energies; // modal energies (mode 0 is usually transrotational)
    }
    version(MHD) {
        Vector3 B;         // magnetic field, Tesla
        number psi;        // divergence cleaning parameter for MHD
        number divB;       // divergence of the magnetic field
    }
    version(komega) {
        number tke;        // turbulent kinetic energy
        number omega;      // omega from k-omega turbulence model
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
        version(komega) {
            tke = other.tke;
            omega = other.omega;
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
        version(komega) {
            tke = src.tke;
            omega = src.omega;
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
        version(komega) {
            tke = 0.0;
            omega = 0.0;
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
        version(komega) {
            tke += other.tke * factor;
            omega += other.omega * factor;
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
        version(komega) {
            tke *= factor;
            omega *= factor;
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
        version(komega) {
            repr ~= ", tke=" ~ to!string(tke);
            repr ~= ", omega=" ~ to!string(omega);
        }
        repr ~= ")";
        return to!string(repr);
    }
} // end class ConservedQuantities

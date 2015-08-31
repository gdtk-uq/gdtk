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
import geom;
import gas;

class ConservedQuantities {
public:
    double mass;         // density, kg/m**3
    Vector3 momentum;    // momentum/unit volume
    Vector3 B;           // magnetic field, Tesla
    double total_energy; // total energy
    double[] massf;      // mass fractions of species
    double[] energies;   // modal energies (mode 0 is usually transrotational)
    double tke;          // turbulent kinetic energy
    double omega;        // omega from k-omega turbulence model
    // [TODO] double[] G;          // velocity dist. partial densities, kg/m**3
    // [TODO] double[] H;          // velocity dist. partial densities, (kg*s**2)/(m**5)

    this(int n_species, int n_modes)
    {
	massf.length = n_species;
	energies.length = n_modes;
    }

    this(in ConservedQuantities other)
    {
	mass = other.mass;
	momentum = other.momentum;
	B = other.B;
	total_energy = other.total_energy;
	massf = other.massf.dup;
	energies = other.energies.dup;
	tke = other.tke;
	omega = other.omega;
    }

    @nogc void copy_values_from(in ConservedQuantities src)
    {
	mass = src.mass;
	momentum.refx = src.momentum.x;
	momentum.refy = src.momentum.y;
	momentum.refz = src.momentum.z;
	B.refx = src.B.x; B.refy = src.B.y; B.refz = src.B.z;
	total_energy = src.total_energy;
	massf[] = src.massf[];
	energies[] = src.energies[];
	tke = src.tke;
	omega = src.omega;
    }

    @nogc void clear_values()
    {
	mass = 0.0;
	momentum.refx = 0.0; momentum.refy = 0.0; momentum.refz = 0.0;
	B.refx = 0.0; B.refy = 0.0; B.refz = 0.0;
	total_energy = 0.0;
	foreach(ref mf; massf) mf = 0.0;
	foreach(ref e; energies) e = 0.0;
	tke = 0.0;
	omega = 0.0;
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
	repr ~= ", tke=" ~ to!string(tke);
	repr ~= ", omega=" ~ to!string(omega);
	repr ~= ")";
	return to!string(repr);
    }
} // end class ConservedQuantities

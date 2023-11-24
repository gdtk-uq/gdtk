/**
 * Authors: Brad Semple
 *
 * Electronic species objects are used by electronic species gas model.
 * An electronic-species represents a binned collection of individually excited
 * electronic states, even if that binned state includes only one individual state.
 *
 */

module gas.electronic_species;

import std.conv : to;
import std.string : format;

import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;

import gas.gas_state;
import gas.physical_constants : electron_volt_energy, Avogadro_number, Boltzmann_constant;

class ElectronicSpecies {
    @property @nogc string name() const { return _name; }
    @property @nogc double mol_mass() const { return _mol_mass; }
    @property @nogc int lowerLevel() const {return _lowerLevel; }
    @property @nogc int upperLevel() const {return _upperLevel; }
    @property @nogc int group_degeneracy() const { return _group_degeneracy; }
    @property @nogc int dof() const { return _dof; }
    @property @nogc double electronic_energy() const {return _electronic_energy; }

    this(lua_State *L)
    {
        _name = getString(L, -1, "name");
        _lowerLevel = getInt(L, -1, "lower_level");
        _upperLevel = getInt(L, -1, "upper_level");
        _mol_mass = getDouble(L, -1, "M");
        _group_degeneracy = getInt(L, -1, "group_degeneracy");
        _dof = getInt(L,-1,"dof");
        _electronic_energy = getDouble(L, -1, "group_energy")*electron_volt_energy*Avogadro_number/_mol_mass;
    }

    // @nogc abstract number group_energy(in GasState Q) const;

private:
    string _name;
    double _mol_mass;
    int _lowerLevel;
    int _upperLevel;
    int _group_degeneracy;
    int _dof;
    double _electronic_energy;
}

// class GroupedElectronicSpecies : ElectronicSpecies {
// public:
//     this(lua_State *L)
//     {
//         super(L);
//         _group_energy = to!number(getDouble(L, -1, "group_energy"));
//         //convert to J/kg
//         _group_energy *= (electron_volt_energy*Avogadro_number)/_mol_mass;
//     }

//     override number group_energy(in GasState Q) const { return _group_energy; }

// private:
//     number _group_energy;
// }

ElectronicSpecies createElectronicSpecies(lua_State *L)
{
    return new ElectronicSpecies(L);
}

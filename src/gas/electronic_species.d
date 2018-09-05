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

import nm.complex;
import nm.number;
import util.lua;
import util.lua_service;

import gas.gas_state;
import gas.physical_constants : electron_volt_energy, Avogadro_number, Boltzmann_constant;

class ElectronicSpecies {
    @property string name() const { return _name; }
    @property int level() const { return _level; }
    @property double mol_mass() const { return _mol_mass; }
    @property int group_degeneracy() const { return _group_degeneracy; }

    this(lua_State *L)
    {
        _name = getString(L, -1, "name");
        _level = getInt(L, -1, "level");
        _mol_mass = getDouble(L, -1, "M");
        _group_degeneracy = getInt(L, -1, "group_degeneracy");
    }

    abstract number energy(in GasState Q) const;

private:
    string _name;
    double _mol_mass;
    int _level;
    int _group_degeneracy;
}

class GroupedElectronicSpecies : ElectronicSpecies {
public:
    this(lua_State *L)
    {
        super(L); //what is this?
        _group_energy = to!number(getDouble(L, -1, "group_energy"));
        //convert to J/kg
        _group_energy *= (electron_volt_energy*Avogadro_number)/_mol_mass;
    }

    override number energy(in GasState Q) const { return _group_energy; }

private:
    number _group_energy;
}

ElectronicSpecies createElectronicSpecies(lua_State *L)
{
    return new GroupedElectronicSpecies(L);
}

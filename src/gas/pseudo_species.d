/**
 * pseudo_species.d
 * Authors: Pierre Mariotto, Rowan G. and Peter J.
 *
 * Pseudo-species objects are used by the pseudo-species gas model.
 * A pseudo-species represents a collection of microstates, even if that
 * collection is only a single microstate.
 *
 */

module gas.pseudo_species;

import std.conv : to;
import std.string : format;

import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;

import gas.gas_state;
import gas.physical_constants : electron_volt_energy, Avogadro_number;


class PseudoSpecies {
    @property @nogc string name() const { return _name; }
    @property @nogc double mol_mass() const { return _mol_mass; }
    @property @nogc int DOF() const { return _dof; }
    @property @nogc int parentIdx() const { return _parentIdx; }

    this(lua_State *L)
    {
        _name = getString(L, -1, "name");
        _mol_mass = getDouble(L, -1, "M");
        _dof = getInt(L, -1, "DOF_base_mode");
        _parentIdx = getInt(L, -1, "parentIdx");
    }

    @nogc abstract number energy(in GasState Q) const;

private:
    string _name;
    double _mol_mass;
    int _dof;
    int _parentIdx;
}

class SingleStatePseudoSpecies : PseudoSpecies {
public:
    this(lua_State *L)
    {
        super(L);
        _energy = to!number(getDouble(L, -1, "energy"));
        // convert to J/kg
        _energy *= (electron_volt_energy*Avogadro_number)/_mol_mass;
    }

    override number energy(in GasState Q) const { return _energy; }

private:
    number _energy;
}

PseudoSpecies createPseudoSpecies(lua_State *L)
{
    auto type = getString(L, -1, "type");
    switch (type) {
    case "single_state":
        return new SingleStatePseudoSpecies(L);
    default:
        string msg = format("The pseudo-species type '%s' could not be created.", type);
        throw new Exception(msg);
    }
}

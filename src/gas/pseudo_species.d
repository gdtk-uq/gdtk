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

import nm.complex;
import nm.number;
import util.lua;
import util.lua_service;

import gas.gas_state;


class PseudoSpecies {
    @property string name() const { return _name; }
    @property double mol_mass() const { return _mol_mass; }
    @property int DOF() const { return _dof; }

    this(lua_State *L)
    {
        _name = getString(L, -1, "name");
        _mol_mass = getDouble(L, -1, "M");
        _dof = getInt(L, -1, "DOF_base_mode");
    }

    abstract number energy(in GasState Q) const;

private:
    string _name;
    double _mol_mass;
    int _dof;
}

class SingleStatePseudoSpecies : PseudoSpecies {
public:
    this(lua_State *L)
    {
        super(L);
        _energy = to!number(getDouble(L, -1, "energy"));
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




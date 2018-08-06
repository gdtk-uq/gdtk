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


interface PseudoSpecies {
    string name();
    number energy(in GasState Q) const;
}

class SingleStatePseudoSpecies : PseudoSpecies {
public:
    this(lua_State *L)
    {
        _name = getString(L, -1, "name");
        _energy = to!number(getDouble(L, -1, "energy"));
    }

    string name() { return _name; }

    number energy(in GasState Q) const { return _energy; }

private:
    string _name;
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




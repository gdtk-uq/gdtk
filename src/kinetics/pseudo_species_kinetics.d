/**
 * Authors: Pierre Mariotto, Rowan G. and Peter J.
 * Date: 2018-08-09
 *
 * This module provides an integrator for the master
 * equation for a state-specific system.
 *
 */

module kinetics.pseudo_species_kinetics;

import nm.complex;
import nm.number;

import gas;
import gas.pseudo_species_gas;
import kinetics.thermochemical_reactor;

final class PseudoSpeciesKinetics : ThermochemicalReactor {
public:
    this(GasModel gmodel)
    {
        super(gmodel);
        _psGmodel = cast(PseudoSpeciesGas) gmodel;
        
    }

    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[] params)
    {
        // [TODO] write update of interval tInterval
    }

private:
    
}

// --------------------------------------------------------------
// Class to provide rates of pseudo-species population changes
// --------------------------------------------------------------

class StateSpecificKineticMechanism {
public:
    @property size_t nReactions() const { return _reactions.length; }

    this(StateSpecificReaction[] reactions) 
    {
        foreach (ref reac; reactions) {
            _reactions ~= reac.dup();
        }
    }

    void evalRateConstants(GasState Q)
    {
        foreach (ref reac; _reactions) {
            reac.evalRateConstants(Q);
        }
    }

    void evalRates(in number[] conc, number[] rates)
    {
        rates[] = to!number(0.0);
        foreach (ir, ref reac; _reactions) {
            foreach (isp; reac.participants) {
                rates[isp] += reac.rateOfChange(isp);
            }
        }
    }

private:
    StateSpecificReaction[] _reactions;
}

StateSpecificReactionMechanism createSSRMechanism(lua_State *L, GasModel gmodel)
{
    StateSpecificReaction[] reactions;
    auto nReactions = to!int(lua_objlen(L, -1));
    foreach (i; 0 .. nReactions) {
        lua_rawgeti(L, -1, i);
        reactions ~= creationSSReaction(L);
        lua_pop(L, 1);
    }
    return new StateSpecificReactionMechanism(reactions);
}

// --------------------------------------------------------------
// Classes to provide behaviour dependent on reaction type.
// Reactions include:
//   + DissociationByAtom
//   + [TODO] DissociationByMolecule 
// --------------------------------------------------------------

class StateSpecificReaction {
    @property number k_f() const { return _k_f; }
    @property number k_b() const { return _k_b; }
    @property ref int[] participants() { return _participants; }

    this(lua_State *L)
    {
        lua_getfield(L, -1, "frc");
        _forward = createSSRateConstant(L);
        lua_pop(L, 1);

        lua_getfield(L, -1, "frc");
        _backward = createSSRateConstant(L);
        lua_pop(L, 1);
    }

    void evalRateConstants(GasState Q)
    {
        _forwardRateConstant.eval(Q);
        _backwardRateConstant.eval(Q);
            
    }

    final void evalRates(GasState Q)
    {
        // [TODO] Compute concentrations.
        
        evalForwardRateOfChange(conc);
        evalBackwardRateofChange(conc);
    }

    abstract number rateOfChange(int isp);

protected:
    abstract void evalForwardRateOfChange(conc);
    abstract void evalBackwardRateOfChange(conc);

private:
    int[] _participants;
    number _k_f;
    number _k_b;
    number _w_f;
    number _w_b;
}

class DissociationByAtom : StateSpecificReaction {
public:
    this(lua_State *L)
    {
        super(L);
        _moleculeIdx = getInt(L, -1, "molecule_idx");
        _atomIdx = getInt(L, -1, "atom_idx");
        _participants = [_moleculeIdx, _atomIdx];
    }

    override number rateOfChange(int isp) {
        if (isp == _moleculeIdx) {
            return _w_b - w_f_;
        }
        else if (isp == _atomIdx) {
            return 2*(_w_f - _w_b);
        }
        else {
            return to!number(0.0);
        }
    }
    
protected:
    override number evalForwardRateOfChange(in number[] conc)
    {
        // e.g. N2(v=0) + N(4s) <=> N(4s) + N(4s) + N(4s)
        //      w_f = k_f * [N2] * [N]
        _w_f = _k_f*conc[_moleculeIdx]*conc[_atomIdx];
    }
    override number evalBackwardRateOfChange(in number[] conc)
    {
        // e.g. N2(v=0) + N(4s) <=> N(4s) + N(4s) + N(4s)
        //      w_b = k_b * [N] * [N] * [N]
        _w_b = _k_b*conc[_atomIdx]*conc[_atomIdx]*[_atomIdx];
    }

private:
    int _moleculeIdx;
    int _atomIdx; 
}

StateSpecificReaction createSSReaction(lua_State *L)
{
    auto type = getString(L, -1, "type");
    switch (type) {
    case "dissociation-by-atom":
        return new DissociationByAtom(L);
    case "dissociation-by-molecule":
        throw new Error("NOT IMPLEMENTED: 'dissociation-by-molecule'");
    default:
        string msg = format("The state-specific reaction type '%s' could not be created.", type);
        throw new Exception(msg);
    }
}


// --------------------------------------------------------------
// Classes to provide various rate constant expressions.
// Rate constants include:
//   + Arrhenius
//   + [TODO] Arrhenius2
// --------------------------------------------------------------

interface StateSpecificRateConstant {
    StateSpecificRateConstant dup();
    number eval(in GasState Q);
}

class ArrheniusSSRC : StateSpecificRateConstant {
public:
    this(double A, double n, double C)
    {
        _A = A;
        _n = n;
        _C = C;
    }
    this(lua_State* L)
    {
        _A = getDouble(L, -1, "A");
        _n = getDouble(L, -1, "n");
        _C = getDouble(L, -1, "C");
    }
    ArrheniusSSRC dup()
    {
        return new ArrheniusSSRC(_A, _n, _C);
    }
    override number eval(in GasState Q)
    {
        return _A*pow(Q.T, _n)*exp(-_C/Q.T);
    }
private:
    double _A, _n, _C;
}

StateSpecificRateConstant createSSRateConstant(lua_State* L)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "Arrhenius":
        return new ArrheniusSSRC(L);
    case "Arrhenius2":
        throw new Error("NOT IMPLEMENTED: 'Arrhenius2'");
    default:
        string msg = format("The state-specific rate constant model '%s' could not be created.", model);
        throw new Exception(msg);
    }
}


/**
 * Authors: Pierre Mariotto, Rowan G. and Peter J.
 * Date: 2018-08-09
 *
 * This module provides an integrator for the master
 * equation for a state-specific system.
 *
 */

module kinetics.pseudo_species_kinetics;

import core.stdc.stdlib : exit;
import std.stdio;
import std.conv : to;
import std.string;
import std.math;
import std.algorithm;

import ntypes.complex;
import nm.number;
import nm.smla;
import util.lua;
import util.lua_service;
import gas;
import gas.pseudo_species_gas;

import kinetics.thermochemical_reactor;

extern(C) {
    // Wrapper of the fortran function `pseudosp_solve_ode`
    @nogc
    void solveODE(number *Y, size_t *neq, number *temp_init, double *dt);
}

final class PseudoSpeciesKinetics : ThermochemicalReactor {
public:

    this(GasModel gmodel)
    {
        super(gmodel);
        _psGmodel = cast(PseudoSpeciesGas) gmodel;
        _conc.length = gmodel.n_species;
        _conc0.length = gmodel.n_species;
    }

    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        _psGmodel.massf2conc(Q, _conc0);
        foreach (i; 0 .. _conc.length) {
            _conc[i] = _conc0[i];
        }

        size_t len=_conc.length;

        // Now update the concentration vector at time t+t_interval
        // by solving the kinetics ode (convection is frozen)
        // solveODE updates _conc to its new value
        // Pass pointers to the fortran subroutine
        solveODE(_conc.ptr, &len, &Q.T , &tInterval);

        _psGmodel.conc2massf(_conc, Q);

    }

    @nogc override eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        string errMsg = "eval_source_terms not implemented for pseudo_species_kinetics.";
        throw new ThermochemicalReactorUpdateException(errMsg);
    }

private:
    PseudoSpeciesGas _psGmodel;
    number[] _conc;
    number[] _conc0;
}

// --------------------------------------------------------------
// Class to provide rates of pseudo-species population changes
// --------------------------------------------------------------

// TODO The classes below should be removed because the fortran file pseudosp_solve_ode.f90
// is now in charge of updating the concentration vector.
// Nevertheless I let it there in case one decides to use it the future

class StateSpecificReactionMechanism {
public:
    @property size_t nReactions() const { return _reactions.length; }

    this(StateSpecificReaction[] reactions)
    {
        foreach (reac; reactions) {
            _reactions ~= reac.dup();
        }
    }

    void evalRateConstants(ref GasState Q)
    {
        foreach (ref reac; _reactions) {
            reac.evalRateConstants(Q);
        }
    }

    void evalRates(in number[] conc, number[] rates)
    {
        foreach (isp; 0 .. rates.length) {
            rates[isp] = to!number(0.0);
        }
        foreach (ir, ref reac; _reactions) {
            reac.evalRates(conc); // evaluate kf and kb
            // Call the method rateOfChange which compute the contribution of
            // the current reaction on the net production rate dni/dt
            reac.rateOfChange(rates); // evaluate dni/dt (RHS of the ODE)
        }
    }

    void initJacobian(SMatrix!number Jac, int nSpecies)
    {
        foreach (isp; 0 .. nSpecies) {
            bool[size_t] columnEntries;
            foreach (jsp; 0 .. nSpecies) {
                foreach (ir, reac; _reactions) {
                    if (reac.participants.canFind(isp) &&
                        reac.participants.canFind(jsp)) {
                        columnEntries[jsp] = true;
                    }
                }
            }
            size_t[] ji = columnEntries.keys.dup();
            sort(ji);
            number[] ai;
            foreach (i; 0 .. ji.length) {
                ai ~= to!number(0.0);
            }
            Jac.addRow(ai, ji);
        }
    }

    void evalJacobian(in number[] conc, SMatrix!number Jac)
    {
        foreach (ref reac; _reactions) {
            reac.evalJacobianEntries(conc);
        }

        Jac.scale(to!number(0.0)); // initialize the Jacobian at 0.
        // Now add contribution of each reaction to the Jacobian
        foreach (reac; _reactions) {
                reac.getJacobianEntry(Jac);
            }
    }

private:
    StateSpecificReaction[] _reactions;
}

StateSpecificReactionMechanism createSSRMechanism(lua_State *L, GasModel gmodel)
{
    StateSpecificReaction[] reactions;

    int nReactions = getInt(L, "number_of_reactions");

    lua_getglobal(L, "reaction");
    foreach (i; 0 .. nReactions) {
        lua_rawgeti(L, -1, i);
        reactions ~= createSSReaction(L);
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
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

    this() {}

    this(lua_State* L)
    {
        lua_getfield(L, -1, "frc");
        _forwardRateConstant = createSSRateConstant(L);
        lua_pop(L, 1);

        lua_getfield(L, -1, "brc");
        _backwardRateConstant = createSSRateConstant(L);
        lua_pop(L, 1);
    }

    abstract StateSpecificReaction dup();

    void evalRateConstants(ref GasState Q)
    {
        _k_f = _forwardRateConstant.eval(Q);
        _k_b = _backwardRateConstant.eval(Q);
    }

    final void evalRates(in number[] conc)
    {
        evalForwardRateOfChange(conc);
        evalBackwardRateOfChange(conc);
    }

    abstract void rateOfChange(ref number[] rates);

    abstract void evalJacobianEntries(in number[] conc);

    abstract void getJacobianEntry(ref SMatrix!number Jac);

protected:
    abstract void evalForwardRateOfChange(in number[] conc);
    abstract void evalBackwardRateOfChange(in number[] conc);

private:
    int[] _participants;
    StateSpecificRateConstant _forwardRateConstant;
    StateSpecificRateConstant _backwardRateConstant;
    number _k_f;
    number _k_b;
    number _w_f;
    number _w_b;
}

class DissociationByAtom : StateSpecificReaction {
public:
    this(StateSpecificRateConstant forward, StateSpecificRateConstant backward,
         int[] participants, int molecule_idx, int atom_idx)
    {
        _forwardRateConstant = forward.dup();
        _backwardRateConstant = backward.dup();
        _participants = participants.dup();
        _moleculeIdx = molecule_idx;
        _atomIdx = atom_idx;
    }
    this(lua_State* L)
    {
        super(L);
        _moleculeIdx = getInt(L, -1, "molecule_idx");
        _atomIdx = getInt(L, -1, "atom_idx");
        _participants = [_moleculeIdx, _atomIdx];
    }

    override DissociationByAtom dup()
    {
        return new DissociationByAtom(_forwardRateConstant, _backwardRateConstant, _participants,
                                      _moleculeIdx, _atomIdx);
    }

    override void rateOfChange(ref number[] rates) {
        // Store net rate of species i in rates[i]
        number w = _w_f - _w_b; // rate of progress
        rates[_moleculeIdx] = rates[_moleculeIdx] - w;
        rates[_atomIdx] = rates[_atomIdx] + 2*w;
    }

    override void evalJacobianEntries(in number[] conc)
    {
        _dwf_dmolc = _k_f*conc[_atomIdx];
        _dwf_datom = _k_f*conc[_moleculeIdx];
        _dwb_datom = 3*_k_b*conc[_atomIdx]*conc[_atomIdx];
    }

    override void getJacobianEntry(ref SMatrix!number Jac) // [TODO] change the name of this method (UpdateJacobian?)
    {
        // Contribution of reaction :
        //      molecule + atom <=> atom + atom + atom
        // on the jacobian
        Jac[_moleculeIdx,_moleculeIdx] = Jac[_moleculeIdx,_moleculeIdx] - _dwf_dmolc;
        Jac[_moleculeIdx,_atomIdx] = Jac[_moleculeIdx,_atomIdx] + _dwb_datom - _dwf_datom;
        Jac[_atomIdx,_moleculeIdx] = Jac[_atomIdx,_moleculeIdx] + 2*_dwf_dmolc;
        Jac[_atomIdx,_atomIdx] = Jac[_atomIdx,_atomIdx] + 2*(_dwf_datom - _dwb_datom);

    }

protected:
    override void evalForwardRateOfChange(in number[] conc)
    {
        // e.g. N2(v=0) + N(4s) <=> N(4s) + N(4s) + N(4s)
        //      w_f = k_f * [N2] * [N]
        _w_f = _k_f*conc[_moleculeIdx]*conc[_atomIdx];
    }
    override void evalBackwardRateOfChange(in number[] conc)
    {
        // e.g. N2(v=0) + N(4s) <=> N(4s) + N(4s) + N(4s)
        //      w_b = k_b * [N] * [N] * [N]
        _w_b = _k_b*conc[_atomIdx]*conc[_atomIdx]*conc[_atomIdx];
    }

private:
    int _moleculeIdx;
    int _atomIdx;
    number _dwf_dmolc;
    number _dwf_datom;
    number _dwb_datom;
}

class Collision2B2B : StateSpecificReaction {
  // Collision 2-bodies-2-bodies
  // r0 + r1 -> p0 + p1

public:
    this(StateSpecificRateConstant forward, StateSpecificRateConstant backward,
         int[] participants, int[] reactants_idx, int[] products_idx)
    {
        _forwardRateConstant = forward.dup();
        _backwardRateConstant = backward.dup();
        _participants = participants.dup();
        _reactants_idx = reactants_idx;
        _products_idx = products_idx;
    }
    this(lua_State* L)
    {
        super(L);
        getArrayOfInts(L, -1, "reactants_idx", _reactants_idx);
        getArrayOfInts(L, -1, "products_idx", _products_idx);
        _participants = _reactants_idx ~ _products_idx; // concatenate indices of reactants and products
    }

    override Collision2B2B dup()
    {
        return new Collision2B2B(_forwardRateConstant, _backwardRateConstant, _participants,
                                 _reactants_idx, _products_idx);
    }

    override void rateOfChange(ref number[] rates) {
        // Store net rate of species i in rates[i]
        number w = _w_f - _w_b;
        rates[_reactants_idx[0]] = rates[_reactants_idx[0]] - w;
        rates[_reactants_idx[1]] = rates[_reactants_idx[1]] - w;
        rates[_products_idx[0]] = rates[_products_idx[0]] + w;
        rates[_products_idx[1]] = rates[_products_idx[1]] + w;
    }

    override void evalJacobianEntries(in number[] conc)
    {
        // Let's call R the rate of progress of the 2-body-2-body reaction
        // nr0, nr1, np0, np1 are reactants and products density number [m^-3]
        // R = kf*nr0*nr1-kb*np0*np1 [m^3/s]
        // Let's compute the element, dR/dnr0, dR/dnr1, dR/dnp0, dR/dnp1
        _dRdnr0 = _k_f*conc[_reactants_idx[1]];
        _dRdnr1 = _k_f*conc[_reactants_idx[0]];
        _dRdnp0 = -_k_b*conc[_products_idx[1]];
        _dRdnp1 = -_k_b*conc[_products_idx[0]];
    }

    override void getJacobianEntry(ref SMatrix!number Jac)
    {

        int r0 = _reactants_idx[0];
        int r1 = _reactants_idx[1];
        int p0 = _products_idx[0];
        int p1 = _products_idx[1];

        // r0 + r1 -> p0 + p1 contributes to Jacobian[isp, jsp] with
        // (isp, jsp) = (r0,r0), (r0,r1), (r0,p0), (r0,p1),
        //              (r1,r0), (r1,r1), (r1,p0), (r1,p1)
        Jac[r0,r0] = Jac[r0,r0] - _dRdnr0;
        Jac[r0,r1] = Jac[r0,r1] - _dRdnr1;
        Jac[r0,p0] = Jac[r0,p0] - _dRdnp0;
        Jac[r0,p1] = Jac[r0,p1] - _dRdnp1;
        Jac[r1,r0] = Jac[r1,r0] - _dRdnr0;
        Jac[r1,r1] = Jac[r1,r1] - _dRdnr1;
        Jac[r1,p0] = Jac[r1,p0] - _dRdnp0;
        Jac[r1,p1] = Jac[r1,p1] - _dRdnp1;

        // ... and (isp, jsp) = (p0,r0), (p0,r1), (p0,p0), (p0,p1),
        //                      (p1,r0), (p1,r1), (p1,p0), (p1,p1)
        Jac[p0,r0] = Jac[p0,r0] + _dRdnr0;
        Jac[p0,r1] = Jac[p0,r1] + _dRdnr1;
        Jac[p0,p0] = Jac[p0,p0] + _dRdnp0;
        Jac[p0,p1] = Jac[p0,p1] + _dRdnp1;
        Jac[p1,r0] = Jac[p1,r0] + _dRdnr0;
        Jac[p1,r1] = Jac[p1,r1] + _dRdnr1;
        Jac[p1,p0] = Jac[p1,p0] + _dRdnp0;
        Jac[p1,p1] = Jac[p1,p1] + _dRdnp1;
    }

protected:
    override void evalForwardRateOfChange(in number[] conc)
    {
        // e.g. r0 + r1 <=> p0 + p1
        //      w_f = k_f * [r0] * [r1]
        _w_f = _k_f*conc[_reactants_idx[0]]*conc[_reactants_idx[1]];
    }
    override void evalBackwardRateOfChange(in number[] conc)
    {
        // e.g. r0 + r1 <=> p0 + p1
        //      w_b = k_b * [p0] * [p1]
        _w_b = _k_b*conc[_products_idx[0]]*conc[_products_idx[1]];
    }

private:
    int[] _reactants_idx;
    int[] _products_idx;
    number _dRdnr0; // derivative of rate of change with density number of r0
    number _dRdnr1; // derivative of rate of change with density number of r1
    number _dRdnp0; // derivative of rate of change with density number of p0
    number _dRdnp1; // derivative of rate of change with density number of p1
}



class Collision2B3B : StateSpecificReaction {
  // Collision 2-bodies-2-bodies
  // r0 + r1 -> p0 + p1 + p2

public:
    this(StateSpecificRateConstant forward, StateSpecificRateConstant backward,
         int[] participants, int[] reactants_idx, int[] products_idx)
    {
        _forwardRateConstant = forward.dup();
        _backwardRateConstant = backward.dup();
        _participants = participants.dup();
        _reactants_idx = reactants_idx;
        _products_idx = products_idx;
    }
    this(lua_State* L)
    {
        super(L);
        getArrayOfInts(L, -1, "reactants_idx", _reactants_idx);
        getArrayOfInts(L, -1, "products_idx", _products_idx);
        _participants = _reactants_idx ~ _products_idx; // concatenate indices of reactants and products
    }

    override Collision2B3B dup()
    {
        return new Collision2B3B(_forwardRateConstant, _backwardRateConstant, _participants,
                                 _reactants_idx, _products_idx);
    }

    override void rateOfChange(ref number[] rates) {
        // Store net rate of species i in rates[i]
        number w = _w_f - _w_b;
        rates[_reactants_idx[0]] = rates[_reactants_idx[0]] - w;
        rates[_reactants_idx[1]] = rates[_reactants_idx[1]] - w;
        rates[_products_idx[0]] = rates[_products_idx[0]] + w;
        rates[_products_idx[1]] = rates[_products_idx[1]] + w;
        rates[_products_idx[2]] = rates[_products_idx[2]] + w;
    }

    override void evalJacobianEntries(in number[] conc)
    {
        // Let's call R the rate of progress
        // of the 2-body-3-body reaction
        // nr0, nr1, np0, np1, np2 are reactants and products density number [m^-3]
        // R = kf*nr0*nr1-kb*np0*np1*np2 [m^3/s]
        // Let's compute the element, dR/dnr0, dR/dnr1, dR/dnp0, dR/dnp1, dR/dnp2
        _dRdnr0 = _k_f*conc[_reactants_idx[1]];
        _dRdnr1 = _k_f*conc[_reactants_idx[0]];
        _dRdnp0 = -_k_b*conc[_products_idx[1]]*conc[_products_idx[2]];
        _dRdnp1 = -_k_b*conc[_products_idx[0]]*conc[_products_idx[2]];
        _dRdnp2 = -_k_b*conc[_products_idx[0]]*conc[_products_idx[1]];
    }

    override void getJacobianEntry(ref SMatrix!number Jac)
    {

        int r0 = _reactants_idx[0];
        int r1 = _reactants_idx[1];
        int p0 = _products_idx[0];
        int p1 = _products_idx[1];
        int p2 = _products_idx[2];

        // r0 + r1 -> p0 + p1 + p2 contributes to Jacobian[isp, jsp] with
        // (isp, jsp) = (r0,r0), (r0,r1), (r0,p0), (r0,p1), (r0,p2),
        //              (r1,r0), (r1,r1), (r1,p0), (r1,p1), (r1,p2)
        Jac[r0,r0] = Jac[r0,r0] - _dRdnr0;
        Jac[r0,r1] = Jac[r0,r1] - _dRdnr1;
        Jac[r0,p0] = Jac[r0,p0] - _dRdnp0;
        Jac[r0,p1] = Jac[r0,p1] - _dRdnp1;
        Jac[r0,p2] = Jac[r0,p2] - _dRdnp2;
        Jac[r1,r0] = Jac[r1,r0] - _dRdnr0;
        Jac[r1,r1] = Jac[r1,r1] - _dRdnr1;
        Jac[r1,p0] = Jac[r1,p0] - _dRdnp0;
        Jac[r1,p1] = Jac[r1,p1] - _dRdnp1;
        Jac[r1,p2] = Jac[r1,p2] - _dRdnp2;

        // ... and (isp, jsp) = (p0,r0), (p0,r1), (p0,p0), (p0,p1), (p0,p2)
        //                      (p1,r0), (p1,r1), (p1,p0), (p1,p1), (p1,p2)
        //                      (p2,r0), (p2,r1), (p2,p0), (p2,p1), (p2,p2)
        Jac[p0,r0] = Jac[p0,r0] + _dRdnr0;
        Jac[p0,r1] = Jac[p0,r1] + _dRdnr1;
        Jac[p0,p0] = Jac[p0,p0] + _dRdnp0;
        Jac[p0,p1] = Jac[p0,p1] + _dRdnp1;
        Jac[p0,p2] = Jac[p0,p2] + _dRdnp2;
        Jac[p1,r0] = Jac[p1,r0] + _dRdnr0;
        Jac[p1,r1] = Jac[p1,r1] + _dRdnr1;
        Jac[p1,p0] = Jac[p1,p0] + _dRdnp0;
        Jac[p1,p1] = Jac[p1,p1] + _dRdnp1;
        Jac[p1,p2] = Jac[p1,p2] + _dRdnp2;
        Jac[p2,r0] = Jac[p2,r0] + _dRdnr0;
        Jac[p2,r1] = Jac[p2,r1] + _dRdnr1;
        Jac[p2,p0] = Jac[p2,p0] + _dRdnp0;
        Jac[p2,p1] = Jac[p2,p1] + _dRdnp1;
        Jac[p2,p2] = Jac[p2,p2] + _dRdnp2;
    }

protected:
    override void evalForwardRateOfChange(in number[] conc)
    {
        // e.g. r0 + r1 <=> p0 + p1 + p2
        //      w_f = k_f * [r0] * [r1]
        _w_f = _k_f*conc[_reactants_idx[0]]*conc[_reactants_idx[1]];
    }
    override void evalBackwardRateOfChange(in number[] conc)
    {
        // e.g. r0 + r1 <=> p0 + p1 + p2
        //      w_b = k_b * [p0] * [p1] * [p2]
        _w_b = _k_b*conc[_products_idx[0]]*conc[_products_idx[1]]*conc[_products_idx[2]];
    }

private:
    int[] _reactants_idx;
    int[] _products_idx;
    number _dRdnr0;
    number _dRdnr1;
    number _dRdnp0;
    number _dRdnp1;
    number _dRdnp2;
}


StateSpecificReaction createSSReaction(lua_State *L)
{
    auto type = getString(L, -1, "type");
    switch (type) {
    case "dissociation-by-atom":
        return new DissociationByAtom(L);
    case "collision-2B-2B":
        return new Collision2B2B(L);
    case "collision-2B-3B":
        return new Collision2B3B(L);
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

version(pseudo_species_kinetics_test) {
    int main()
    {
        import util.msg_service;

        // Set up gas model
        auto L = init_lua_State();
        doLuaFile(L, "../gas/sample-data/pseudo-species-3-components.lua");
        auto gm = new PseudoSpeciesGas(L);
        auto gd = GasState(2, 0);
        gd.massf[0] = 0.2;
        gd.massf[1] = 0.8;
        gd.p = 1.0e5;
        gd.T = 4000.0;
        gm.update_thermo_from_pT(gd);
        lua_close(L);

        // Set up kinetics object
        auto L2 = init_lua_State();
        writeln("Parsing kinetics file.");
        doLuaFile(L2, "sample-input/state-specific-N2-diss.lua");

        PseudoSpeciesKinetics psk = new PseudoSpeciesKinetics(L2, gm);

        writeln("Gas state BEFORE update: ", gd);

        double tInterval = 1.0e-6;
        double dtSuggest = -1.0;
        number[] params;

        psk(gd, 1.0e-6, dtSuggest, params);

        writeln("Gas state AFTER update: ", gd);
        // Apply thermo constraint.
        //gm.update_thermo_from_rhou(gd);

        //writeln("Gas state AFTER update: ", gd);

        return 0;

    }
}

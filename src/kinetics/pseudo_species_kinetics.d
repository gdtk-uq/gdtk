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

import nm.complex;
import nm.number;
import nm.smla;
import util.lua;
import util.lua_service;
import gas;
import gas.pseudo_species_gas;

import kinetics.thermochemical_reactor;

final class PseudoSpeciesKinetics : ThermochemicalReactor {
public:
    @property size_t nReactions() const { return _mech.nReactions(); }

    this(string fileName, GasModel gmodel)
    {
        super(gmodel);
        _psGmodel = cast(PseudoSpeciesGas) gmodel;
        auto L = init_lua_State();
        doLuaFile(L, fileName);
        _mech = createSSRMechanism(L, gmodel);
        lua_close(L);
        _conc.length = gmodel.n_species;
        _conc0.length = gmodel.n_species;
        _negFy.length = gmodel.n_species;
        _dconc.length = gmodel.n_species;
        _dconc0.length = gmodel.n_species;
        _rates.length = gmodel.n_species;
        _Jac = new SMatrix!number();
        _mech.initJacobian(_Jac, gmodel.n_species);
        _ILU0 = new SMatrix!number(_Jac);
        _gws = GMRESWorkSpace!number(gmodel.n_species, _max_iterations);
    }

    this(lua_State* L, GasModel gmodel)
    {
        super(gmodel);
        _psGmodel = cast(PseudoSpeciesGas) gmodel;
        _mech = createSSRMechanism(L, gmodel);
        _conc.length = gmodel.n_species;
        _conc0.length = gmodel.n_species;
        _negFy.length = gmodel.n_species;
        _dconc.length = gmodel.n_species;
        _dconc0.length = gmodel.n_species;
        _rates.length = gmodel.n_species;
        _Jac = new SMatrix!number();
        _mech.initJacobian(_Jac, gmodel.n_species);
        _ILU0 = new SMatrix!number(_Jac);
        _gws = GMRESWorkSpace!number(gmodel.n_species, _max_iterations);
    }

    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[] params)
    {
        _psGmodel.massf2conc(Q, _conc0);
        foreach (i; 0 .. _conc.length) {
            _conc[i] = _conc0[i];
        }
        _mech.evalRateConstants(Q);
        foreach (i; 0 .. _dconc.length) {
            _dconc[i] = to!number(0.0);
        }
        int stepCount = 0;
        foreach (step; 0 .. _max_steps) {
            stepCount++;
            // 1. Evaluate negF(Y)
            _mech.evalRates(_conc, _rates);
            foreach (i; 0 .. _rates.length) {
                _negFy[i] = -1.0*(_conc[i] - _conc0[i] - tInterval*_rates[i]);
            }
            // 2. Evaluate dfdY
            _mech.evalJacobian(_conc, _Jac);
            // 3. Form Jacobian as [I - dt*dfdy];
            _Jac.scale(to!number(-1.0*tInterval));
            // 3a. Add I on diagonal.
            foreach (isp; 0 .. _conc.length) {
                _Jac[isp,isp] = _Jac[isp,isp] + 1.0;
            }
            // 4. Compute dY
            // 4a. Compute ILU0 preconditioner.
            foreach (i; 0 .. _Jac.aa.length) {
                _ILU0.aa[i] = _Jac.aa[i];
            }
            decompILU0!number(_ILU0);
            // 4b. Use right-preconditioned GMRES to solve for dY
            //     Use previous _dconc as guess.
            foreach (i; 0 .. _dconc.length) {
                _dconc0[i] = _dconc[i];
            }
            rpcGMRES(_Jac, _ILU0, _negFy, _dconc0, _dconc,
                     _max_iterations, _GMRES_tol, _gws);
            // 5. Y_new = Y_old + dY
            foreach (i; 0 .. _conc.length) {
                _conc[i] = _conc[i] + _dconc[i];
            }
            // 6. Check for convergence.
            // [TODO] PM: We may need a more sophisticated convergence
            //            check. For the moment, I'm just checking
            //            for small changes in dY.
            double largestChange = 0.0;
            foreach (i; 0 .. _dconc.length) {
                if (fabs(_dconc[i].re) > largestChange) {
                    largestChange = fabs(_dconc[i].re);
                }
            }
            if (largestChange < _Newton_tol) {
                break;
            }
        }
        if (stepCount == _max_steps) {
            string excMsg = "Newton steps failed to converge in pseudo_species_kinetics.d\n";
            excMsg ~= format("Number of steps: %d\n", _max_steps);
            throw new ThermochemicalReactorUpdateException(excMsg);
        }
        _psGmodel.conc2massf(_conc, Q); 
    }

private:
    int _max_steps = 10;
    int _max_iterations = 30; // in matrix solve
    double _Newton_tol = 1.0e-9;
    double _GMRES_tol = 1.0e-15;
    PseudoSpeciesGas _psGmodel;
    StateSpecificReactionMechanism _mech;
    number[] _conc;
    number[] _conc0;
    number[] _negFy;
    number[] _dconc;
    number[] _dconc0;
    number[] _rates;
    SMatrix!number _Jac;
    SMatrix!number _ILU0;
    GMRESWorkSpace!number _gws;
}

// --------------------------------------------------------------
// Class to provide rates of pseudo-species population changes
// --------------------------------------------------------------

class StateSpecificReactionMechanism {
public:
    @property size_t nReactions() const { return _reactions.length; }

    this(StateSpecificReaction[] reactions) 
    {
        foreach (reac; reactions) {
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
        foreach (isp; 0 .. rates.length) {
            rates[isp] = to!number(0.0);
        }
        foreach (ir, ref reac; _reactions) {
            reac.evalRates(conc);
            foreach (isp; reac.participants) {
                rates[isp] += reac.rateOfChange(isp);
            }
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
        
        Jac.scale(to!number(0.0));
        foreach (isp; 0 .. to!int(conc.length)) {
            foreach (jsp; 0 .. to!int(conc.length)) {
                foreach (reac; _reactions) {
                    if (reac.participants.canFind(isp) && 
                        reac.participants.canFind(jsp)) {
                        Jac[isp,jsp] = Jac[isp,jsp] + reac.getJacobianEntry(isp, jsp);
                    }
                }
            }
        }
    }

private:
    StateSpecificReaction[] _reactions;
}

StateSpecificReactionMechanism createSSRMechanism(lua_State *L, GasModel gmodel)
{
    StateSpecificReaction[] reactions;
    
    int nReactions = getInt(L, LUA_GLOBALSINDEX, "number_of_reactions");

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

    void evalRateConstants(GasState Q)
    {
        _k_f = _forwardRateConstant.eval(Q);
        _k_b = _backwardRateConstant.eval(Q);
    }

    final void evalRates(in number[] conc)
    {
        evalForwardRateOfChange(conc);
        evalBackwardRateOfChange(conc);
    }

    abstract number rateOfChange(int isp);

    abstract void evalJacobianEntries(in number[] conc);

    abstract number getJacobianEntry(int isp, int jsp);

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

    override number rateOfChange(int isp) {
        if (isp == _moleculeIdx) {
            return _w_b - _w_f;
        }
        else if (isp == _atomIdx) {
            return 2*(_w_f - _w_b);
        }
        else {
            return to!number(0.0);
        }
    }

    override void evalJacobianEntries(in number[] conc)
    {
        _dwf_dmolc = _k_f*conc[_atomIdx];
        _dwf_datom = _k_f*conc[_moleculeIdx];
        _dwb_datom = 3*_k_b*conc[_atomIdx]*conc[_atomIdx];
    }
    
    override number getJacobianEntry(int isp, int jsp)
    {
        if (isp == _moleculeIdx) {
            if (jsp == _moleculeIdx) {
                return -_dwf_dmolc;
            }
            else if (jsp == _atomIdx) {
                return _dwb_datom - _dwf_datom;
            }
            else {
                return to!number(0.0);
            }
        }
        else if (isp == _atomIdx) {
            if (jsp == _moleculeIdx) {
                return 2*_dwf_dmolc;
            }
            else if (jsp == _atomIdx) {
                return 2*(_dwf_datom - _dwb_datom);
            }
            else {
                return to!number(0.0);
            }
        }
        return to!number(0.0);
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

version(pseudo_species_kinetics_test) {
    int main()
    {
        import util.msg_service;
        
        // Set up gas model
        auto L = init_lua_State();
        doLuaFile(L, "../gas/sample-data/pseudo-species-3-components.lua");
        auto gm = new PseudoSpeciesGas(L);
        auto gd = new GasState(2, 0);
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
        double dtChemSuggest = -1.0;
        double dtThermSuggest = -1.0;
        number[] params;
        
        psk(gd, 1.0e-6, dtChemSuggest, dtThermSuggest, params);
        
        writeln("Gas state AFTER update: ", gd);
        // Apply thermo constraint.
        //gm.update_thermo_from_rhou(gd);

        //writeln("Gas state AFTER update: ", gd);
        
        return 0;

    }
}

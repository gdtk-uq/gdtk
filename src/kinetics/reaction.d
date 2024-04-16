/**
 * reaction.d
 * This module holds the classes related to computing
 * rates of species change due to a single reaction.
 *
 * Author: Rowan J. Gollan
 */

module kinetics.reaction;

import std.conv;
import std.math;
import std.string;
import std.typecons;
import std.algorithm;
import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;
import util.msg_service;
import kinetics.rate_constant;
import kinetics.thermochemical_reactor : ThermochemicalReactorUpdateException;

@nogc
number compute_equilibrium_constant(GasModel gmodel, ref GasState Q, in number[] gibbs_energies,
                                    in int[] participants, in int[] nu, bool recompute_gibbs)
{
    // Note on handling rate controlling temperatures: 
    // The GasState passed in here should have all its temperatures replaced with 
    // the rate controlling temperature. So we can expect Q.T to be the rate controlling temp.
    // We do have to re-compute the gibbs free energy if we are using a different 
    // rate controlling temperature. We can simply pass Q to the gas model, and it can access
    // the rate controlling temperature through Q.T
    number dG = 0.0;
    number nuSum = 0.0;
    foreach ( isp; participants ) {
        number gibbs_energy = (recompute_gibbs) ? gmodel.gibbs_free_energy(Q, isp) : gibbs_energies[isp];
        dG += nu[isp] * gibbs_energy * gmodel.mol_masses()[isp];
        nuSum += nu[isp];
    }
    number K_p = exp(-dG/(R_universal*Q.T));
    number K_c = K_p*pow(P_atm/(R_universal*Q.T), nuSum);
    debug{
    if (K_c==0.0) {
        string msg=format("Bad equilibrium constant %e detected in reaction %s", K_c, participants);
        throw new ThermochemicalReactorUpdateException(msg);
    }
    }
    return K_c;
}

/++
 Reaction is an interface specifying the public services
 provided by an object of Reaction type.
+/
class Reaction
{
public:
    @property @nogc number k_f() const { return _k_f; }
    @property @nogc number k_b() const { return _k_b; }
    @property @nogc number K_eq() const { return _K_eq; }
    @property @nogc ref int[] participants() { return _participants; }

    this(RateConstant forward, RateConstant backward, GasModel gmodel, int rctIndex_f=-1, int rctIndex_b=-1)
    {
        _gmodel = gmodel;
        _Qw = GasState(gmodel);

        if ( forward is null ) {
            _compute_kf_then_kb = false;
            _forward = null;
            _backward = backward.dup();
        }
        else if ( backward is null ) {
            _compute_kf_then_kb = true;
            _backward = null;
            _forward = forward.dup();
        }
        else { // we have both forward and backward rate constants.
            _compute_kf_then_kb = true;
            _forward = forward.dup();
            _backward = backward.dup();
        }
        _rctIndex_f = rctIndex_f;
        _rctIndex_b = rctIndex_b;
    }

    abstract Reaction dup();

    @nogc abstract void eval_equilibrium_constant(in GasState Q, in number[] gibbs_energies, int rct_index);

    @nogc final void eval_rate_constants(in GasState Q, in number[] gibbs_energies)
    {
        if ( _compute_kf_then_kb ) {
            _k_f = eval_forward_rate_constant(Q);
            if ( _backward is null ) { // we need the equilibrium constant
                // Let's take assume that n_modes >= 1 indicates a multi-temperature
                // simulation. In which case, we need to evaluate the k_f at equilibrium
                // with the rate controlling temperature in order to determine k_b.
                if (_gmodel.n_modes >= 1) {
                    // We need to re-evaluate the forward rate constant and equilibrium
                    // constant at the backwards rate controlling temperature
                    number T = (_rctIndex_b < 0) ? Q.T : Q.T_modes[_rctIndex_b];
                    _Qw.T = T;
                    _Qw.T_modes[] = T;
                    eval_equilibrium_constant(Q, gibbs_energies, _rctIndex_b);
                    number kf_rct = eval_forward_rate_constant(_Qw);
                    _k_b = kf_rct/_K_eq;
                }
                else {
                    eval_equilibrium_constant(Q, gibbs_energies, -1);
                    _k_b = _k_f/_K_eq;
                }
            }
            else {
                _k_b = eval_backward_rate_constant(Q);
            }
        }
        else {
            _k_b = eval_backward_rate_constant(Q);
            if ( _forward is null ) { // we need the equilibrium constant
                if (_gmodel.n_modes >= 1) {
                    // re-evaluate the backward rate constsant at the 
                    // forward rate controlling temperature
                    number T = (_rctIndex_f < 0) ? Q.T : Q.T_modes[_rctIndex_f];
                    _Qw.T = T;
                    _Qw.T_modes[] = T;
                    eval_equilibrium_constant(Q, gibbs_energies, _rctIndex_f);
                    number kb_rct = eval_backward_rate_constant(Q);
                    _k_f = kb_rct * _K_eq;
                }
                else {
                    eval_equilibrium_constant(Q, gibbs_energies, -1);
                    _k_f = _k_b*_K_eq;
                }
            }
            else {
                _k_f = eval_forward_rate_constant(Q);
            }
        }
    }

    @nogc final void eval_rates(in number[] conc)
    {
        _w_f = eval_forward_rate(conc);
        _w_b = eval_backward_rate(conc);
    }

    @nogc final number forward_tickrate() const
    {
        return _w_f;
    }

    @nogc final number backward_tickrate() const
    {
        return _w_b;
    }

    @nogc abstract number production(int isp) const;
    @nogc abstract number loss(int isp) const;

protected:
    @nogc abstract number eval_forward_rate(in number[] conc);
    @nogc abstract number eval_backward_rate(in number[] conc);

private:
    RateConstant _forward, _backward;
    number _k_f, _k_b; // Storage of computed rate constants
    number _K_eq; // Equilibrium constant based on concentration
    bool _compute_kf_then_kb = true;
    number _w_f, _w_b; // Storage of computed rates of change
    GasModel _gmodel;
    GasState _Qw; // a GasState for temporary working
    int[] _participants; // Storage of indices of species that
                         // participate in this reaction
    // hacky way to store the rate controlling temperature
    // if only one of the rates is supplied
    int _rctIndex_f, _rctIndex_b;
    @nogc number eval_forward_rate_constant(in GasState Q)
    {
        return _forward.eval(Q);
    }
    @nogc number eval_backward_rate_constant(in GasState Q)
    {
        return _backward.eval(Q);
    }
}

/++
 An ElementaryReaction is a reaction whose rate law can
 be written directly from its molecularity.
 There are three forms:
   - unimolecular
   - bimolecular
   - termolecular
+/
class ElementaryReaction : Reaction
{
public:
    this(RateConstant forward, RateConstant backward, GasModel gmodel,
         int[] reac_spidx, int[] reac_coeffs, int[] prod_spidx, int[] prod_coeffs,
         size_t n_species, int rctIndex_f, int rctIndex_b)
    {
        assert(reac_spidx.length == reac_coeffs.length,
               brokenPreCondition("reac_spidx and reac_coeffs arrays not of equal length"));
        assert(prod_spidx.length == prod_coeffs.length,
               brokenPreCondition("prod_spdix and prod_coeffs arrays are not of equal length"));
        super(forward, backward, gmodel, rctIndex_f, rctIndex_b);
        bool[int] pmap;
        foreach ( i; 0 .. reac_spidx.length ) {
            _reactants ~= tuple(reac_spidx[i], reac_coeffs[i]);
            pmap[reac_spidx[i]] = true;
        }
        foreach ( i; 0 .. prod_spidx.length ) {
            _products ~= tuple(prod_spidx[i], prod_coeffs[i]);
            pmap[prod_spidx[i]] = true;
        }
        _participants = pmap.keys.dup();
        sort(_participants);
        foreach ( isp; 0 .. n_species ) {
            int nu1 = 0;
            foreach ( ref r; _reactants ) {
                if ( r[0] == isp ) {
                    nu1 = r[1];
                    break;
                }
            }
            int nu2 = 0;
            foreach ( ref p; _products ) {
                if ( p[0] == isp ) {
                    nu2 = p[1];
                    break;
                }
            }
            _nu ~= nu2 - nu1;
        }
    }
    this(RateConstant forward, RateConstant backward, GasModel gmodel, in int[] participants,
         in Tuple!(int, int)[] reactants, in Tuple!(int, int)[] products, in int[] nu)
    {
        super(forward, backward, gmodel);
        _participants = participants.dup();
        _reactants = reactants.dup();
        _products = products.dup();
        _nu = nu.dup();
    }
    override ElementaryReaction dup()
    {
        return new ElementaryReaction(_forward, _backward, _gmodel, _participants,
                                      _reactants, _products, _nu);
    }

    @nogc
    override void eval_equilibrium_constant(in GasState Q, in number[] gibbs_energies, int rct_index)
    {
        
        bool other_rct = rct_index >= 0;
        number T = (other_rct) ? Q.T_modes[rct_index] : Q.T;
        _Qw.T = T;
        _Qw.p = Q.p;
        if (_gmodel.n_modes > 0) { _Qw.T_modes[] = T; }
        _K_eq = compute_equilibrium_constant(_gmodel, _Qw, gibbs_energies, _participants, _nu, other_rct );
    }

    @nogc
    override number production(int isp) const
    {
        if ( _nu[isp] >  0 )
            return _nu[isp]*_w_f;
        else if ( _nu[isp] < 0 )
            return -_nu[isp]*_w_b;
        else
            return to!number(0.0);
    }

    @nogc
    override number loss(int isp) const
    {
        if ( _nu[isp] > 0 )
            return _nu[isp]*_w_b;
        else if ( _nu[isp] < 0 )
            return -_nu[isp]*_w_f;
        else
            return to!number(0.0);
    }
protected:
    @nogc
    override number eval_forward_rate(in number[] conc)
    {
        number val = _k_f;
        foreach ( ref r; _reactants ) {
            int isp = r[0];
            int coeff = r[1];
            val *= pow(conc[isp], coeff);
        }
        return val;
    }

    @nogc
    override number eval_backward_rate(in number[] conc)
    {
        number val = _k_b;
        foreach ( ref p; _products ) {
            int isp = p[0];
            int coeff = p[1];
            val *= pow(conc[isp], coeff);
        }
        return val;
    }

private:
    Tuple!(int, int)[] _reactants;
    Tuple!(int, int)[] _products;
    int[] _nu;
}

/++
 An AnonymousColliderReaction is a reaction that involves collision
 with a so-called anonymous partner. In reaction mechanisms, this
 is typically indicated with a partner "M".
+/
class AnonymousColliderReaction : Reaction
{
public:
    this(RateConstant forward, RateConstant backward, GasModel gmodel,
         int[] reac_spidx, int[] reac_coeffs, int[] prod_spidx, int[] prod_coeffs,
         Tuple!(int,double)[] efficiencies, size_t n_species, int rctIndex_f=-1, int rctIndex_b=-1)
    {
        assert(reac_spidx.length == reac_coeffs.length,
               brokenPreCondition("reac_spidx and reac_coeffs arrays not of equal length"));
        assert(prod_spidx.length == prod_coeffs.length,
               brokenPreCondition("prod_spdix and prod_coeffs arrays are not of equal length"));
        super(forward, backward, gmodel, rctIndex_f, rctIndex_b);
        bool[int] pmap;
        foreach ( i; 0 .. reac_spidx.length ) {
            _reactants ~= tuple(reac_spidx[i], reac_coeffs[i]);
            pmap[reac_spidx[i]] = true;
        }
        foreach ( i; 0 .. prod_spidx.length ) {
            _products ~= tuple(prod_spidx[i], prod_coeffs[i]);
            pmap[prod_spidx[i]] = true;
        }
        _participants = pmap.keys.dup();
        sort(_participants);
        foreach ( isp; 0 .. n_species ) {
            int nu1 = 0;
            foreach ( ref r; _reactants ) {
                if ( r[0] == isp ) {
                    nu1 = r[1];
                    break;
                }
            }
            int nu2 = 0;
            foreach ( ref p; _products ) {
                if ( p[0] == isp ) {
                    nu2 = p[1];
                    break;
                }
            }
            _nu ~= nu2 - nu1;
        }
        _efficiencies = efficiencies.dup();
    }
    this(RateConstant forward, RateConstant backward, GasModel gmodel, in int[] participants,
         in Tuple!(int, int)[] reactants, in Tuple!(int, int)[] products, in int[] nu, in Tuple!(int, double)[] efficiencies)
    {
        super(forward, backward, gmodel);
        _participants = participants.dup();
        _reactants = reactants.dup();
        _products = products.dup();
        _nu = nu.dup();
        _efficiencies = efficiencies.dup();
    }
    override AnonymousColliderReaction dup()
    {
        return new AnonymousColliderReaction(_forward, _backward, _gmodel, _participants,
                                             _reactants, _products, _nu, _efficiencies);
    }

    @nogc
    override void eval_equilibrium_constant(in GasState Q, in number[] gibbs_energies, int rct_index)
    {
        bool other_rct = rct_index >= 0;
        number T = (other_rct) ? Q.T_modes[rct_index] : Q.T;
        _Qw.T = T;
        _Qw.p = Q.p;
        if (_gmodel.n_modes > 0) { _Qw.T_modes[] = T; }
        _K_eq = compute_equilibrium_constant(_gmodel, _Qw, gibbs_energies, _participants, _nu, other_rct );
    }

    @nogc
    override number production(int isp) const
    {
        if ( _nu[isp] >  0 )
            return _nu[isp]*_w_f;
        else if ( _nu[isp] < 0 )
            return -_nu[isp]*_w_b;
        else
            return to!number(0.0);
    }

    @nogc
    override number loss(int isp) const
    {
        if ( _nu[isp] > 0 )
            return _nu[isp]*_w_b;
        else if ( _nu[isp] < 0 )
            return -_nu[isp]*_w_f;
        else
            return to!number(0.0);
    }
protected:
    @nogc
    override number eval_forward_rate(in number[] conc)
    {
        // eval_forward_rate() needs to be called before
        // eval_backward_rate() so that _anaonymousColliderTerm is up-to-date.
        computeAnonymousColliderTerm(conc);
        number val = _k_f*_anonymousColliderTerm;
        foreach ( ref r; _reactants ) {
            int isp = r[0];
            int coeff = r[1];
            val *= pow(conc[isp], coeff);
        }
        return val;
    }

    @nogc
    override number eval_backward_rate(in number[] conc)
    {
        number val = _k_b*_anonymousColliderTerm;
        foreach ( ref p; _products ) {
            int isp = p[0];
            double coeff = p[1];
            val *= pow(conc[isp], coeff);
        }
        return val;
    }

private:
    Tuple!(int, int)[] _reactants;
    Tuple!(int, int)[] _products;
    int[] _nu;
    Tuple!(int, double)[] _efficiencies;
    number _anonymousColliderTerm;

    @nogc
    void computeAnonymousColliderTerm(in number[] conc)
    {
        _anonymousColliderTerm = 0.0;
        foreach ( ref c; _efficiencies ) {
            _anonymousColliderTerm += c[1] * conc[c[0]];
        }
    }
}


/++
 + Creates a Reaction object from information in a Lua table.
 +
 + The table format mirrors that created by reaction.lua::transformReaction()
 + Check that source also.
 + table = {equation = "H2 + I2 <=> 2 HI",
 +          type = "elementary",
 +          frc = transform_rate_constant(fr),
 +          brc = transform_rate_constant(fr),
 +          ec = nil,
 +          reacIdx = {0, 1},
 +          reacCoeffs = {1.0, 1.0},
 +          prodIdx = {2},
 +          prodCoeffs = {2.0},
 +          pressureDependent = false,
 +          efficiencies = {}
 + }
 +/
Reaction createReaction(lua_State* L, GasModel gmodel)
{
    int n_species = gmodel.n_species;
    int rctIndex_f = -1;
    int rctIndex_b = -1;
    // We need to attempt to get table of efficiencies also
    Tuple!(int, double)[] efficiencies;
    lua_getfield(L, -1, "efficiencies");
    if ( !lua_isnil(L, -1) ) {
        lua_pushnil(L);
        while ( lua_next(L, -2) != 0 ) {
            int spIdx = to!int(lua_tointeger(L, -2));
            double efficiency = lua_tonumber(L, -1);
            efficiencies ~= tuple(spIdx, efficiency);
            lua_pop(L, 1);
        }
    }
    lua_pop(L, 1);
    // All Reactions have a forward and backward rate.
    lua_getfield(L, -1, "frc");
    auto frc = createRateConstant(L, efficiencies, gmodel);
    if (frc is null){
        rctIndex_f = getIntWithDefault(L, -1, "rctIndex", -1);
    }
    lua_pop(L, 1);
    lua_getfield(L, -1, "brc");
    auto brc = createRateConstant(L, efficiencies, gmodel);
    if (brc is null){
        rctIndex_b = getIntWithDefault(L, -1, "rctIndex", -1);
    }
    lua_pop(L, 1);

    // And most use reacIdx, reacCoeffs, prodIdx and prodCoeffs lists.
    int[] reacIdx, prodIdx;
    int[] reacCoeffs, prodCoeffs;
    getArrayOfInts(L, -1, "reacIdx", reacIdx);
    getArrayOfInts(L, -1, "reacCoeffs", reacCoeffs);
    getArrayOfInts(L, -1, "prodIdx", prodIdx);
    getArrayOfInts(L, -1, "prodCoeffs", prodCoeffs);

    // We need to specialise the creation of a Reaction
    // based on type.
    auto type = getString(L, -1, "type");
    switch (type) {
    case "elementary":
        return new ElementaryReaction(frc, brc, gmodel, reacIdx, reacCoeffs,
                                      prodIdx, prodCoeffs, n_species, rctIndex_f, rctIndex_b);
    case "anonymous_collider":
        // We need to get table of efficiencies also
        return new AnonymousColliderReaction(frc, brc, gmodel,
                                             reacIdx, reacCoeffs,
                                             prodIdx, prodCoeffs,
                                             efficiencies, n_species);
    default:
        string msg = format("The reaction type: %s is not known.", type);
        throw new Exception(msg);
    }
}


version(reaction_test) {
    int main() {
        GasModel gm = init_gas_model("sample-input/H2-I2-HI.lua");
        // Find rate of forward production for H2 + I2 reaction at 700 K.
        double[] conc = [4.54, 4.54, 0.0];
        double[] gibbs_energies = [0.0, 0.0, 0.0];
        auto rc = new ArrheniusRateConstant(1.94e14, 0.0, 20620.0, -1);
        auto gd = GasState(3, 1);
        gd.T = 700.0;
        gm.gibbs_free_energies(gd, gibbs_energies);
        auto reaction = new ElementaryReaction(rc, rc, gm, [0, 1], [1, 1],
                                               [2], [2], 3, -1, -1);
        reaction.eval_rate_constants(gd, gibbs_energies);
        reaction.eval_rates(conc);
        assert(isClose(0.0, reaction.production(0)), failedUnitTest());
        assert(isClose(0.0, reaction.production(1)), failedUnitTest());
        assert(isClose(1287.8606, reaction.production(2), 1.0e-6), failedUnitTest());
        assert(isClose(643.9303, reaction.loss(0), 1.0e-6), failedUnitTest());
        assert(isClose(643.9303, reaction.loss(1), 1.0e-6), failedUnitTest());
        assert(isClose(0.0, reaction.loss(2)), failedUnitTest());

        auto reaction2 = new AnonymousColliderReaction(rc, rc, gm, [0, 1], [1, 1],
                                                       [2], [2], [tuple(1, 1.0)], 3);
        reaction2.eval_rate_constants(gd, gibbs_energies);
        reaction2.eval_rates(conc);

        // Try a reaction with backwards rate computed from equilibrium constant.
        auto reaction3 = new ElementaryReaction(rc, null, gm, [0, 1], [1, 1],
                                                [2], [2], 3, -1, -1);
        reaction3.eval_rate_constants(gd, gibbs_energies);
        reaction3.eval_rates(conc);

        return 0;
    }
}

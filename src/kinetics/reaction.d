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
import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;
import util.msg_service;
import kinetics.rate_constant;

@nogc
number compute_equilibrium_constant(GasModel gmodel, GasState Q,
                                    in int[] participants, in double[] nu)
{
    Q.p = P_atm; // need to evaluate Gibbs energy at standard state.
    number dG = 0.0;
    number nuSum = 0.0;
    foreach ( isp; participants ) {
        dG += nu[isp] * gmodel.gibbs_free_energy(Q, isp) * gmodel.mol_masses()[isp];
        nuSum += nu[isp];
    }
    number K_p = exp(-dG/(R_universal*Q.T));
    number K_c = K_p*pow(P_atm/(R_universal*Q.T), nuSum);
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
    
    this(RateConstant forward, RateConstant backward, GasModel gmodel)
    {
        _gmodel = gmodel;
        _Qw = new GasState(gmodel);

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
    }

    abstract Reaction dup();

    @nogc abstract void eval_equilibrium_constant(in GasState Q);
        
    @nogc final void eval_rate_constants(in GasState Q)
    {
        immutable double EPS = 1.0e-16; // To prevent divide by zero
                                        // if K_eq is very very small.
        if ( _compute_kf_then_kb ) {
            _k_f = eval_forward_rate_constant(Q);
            if ( _backward is null ) { // we need the equilibrium constant
                eval_equilibrium_constant(Q);
                // Let's take assume that n_modes >= 1 indicates a multi-temperature
                // simulation. In which case, we need to evaluate the k_f at equilibrium
                // with the translational temperature in order to determine k_b.
                if (_gmodel.n_modes >= 1) {
                    _Qw.T = Q.T;
                    _Qw.T_modes[] = Q.T;
                    number kf_eq = eval_forward_rate_constant(_Qw);
                    _k_b = kf_eq/(_K_eq + EPS);
                }
                else {
                    _k_b = _k_f/(_K_eq + EPS);
                }
            }
            else {
                _k_b = eval_backward_rate_constant(Q);
            }
        }
        else {
            _k_b = eval_backward_rate_constant(Q);
            if ( _forward is null ) { // we need the equilibrium constant
                eval_equilibrium_constant(Q);
                _k_f = _k_b*_K_eq;
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
         int[] reac_spidx, double[] reac_coeffs, int[] prod_spidx, double[] prod_coeffs,
         size_t n_species)
    {
        assert(reac_spidx.length == reac_coeffs.length,
               brokenPreCondition("reac_spidx and reac_coeffs arrays not of equal length"));
        assert(prod_spidx.length == prod_coeffs.length,
               brokenPreCondition("prod_spdix and prod_coeffs arrays are not of equal length"));
        super(forward, backward, gmodel);
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
            double nu1 = 0.0;
            foreach ( ref r; _reactants ) {
                if ( r[0] == isp ) {
                    nu1 = r[1];
                    break;
                }
            }
            double nu2 = 0.0;
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
         in Tuple!(int, double)[] reactants, in Tuple!(int, double)[] products, in double[] nu)
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
    override void eval_equilibrium_constant(in GasState Q)
    {
        _Qw.T = Q.T;
	if (_gmodel.n_modes >= 1) _Qw.T_modes[] = Q.T; // equilibrium constant evaluated at thermal equilibrium 
        _K_eq = compute_equilibrium_constant(_gmodel, _Qw, _participants, _nu);
    }
    
    @nogc
    override number production(int isp) const
    {
        if ( _nu[isp] >  0.0 )
            return _nu[isp]*_w_f;
        else if ( _nu[isp] < 0.0 )
            return -_nu[isp]*_w_b;
        else
            return to!number(0.0);
    }

    @nogc
    override number loss(int isp) const
    {
        if ( _nu[isp] > 0.0 ) 
            return _nu[isp]*_w_b;
        else if ( _nu[isp] < 0.0 )
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
            double coeff = r[1].re;
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
            double coeff = p[1].re;
            val *= pow(conc[isp], coeff);
        }
        return val;
    }

private:
    Tuple!(int, double)[] _reactants;
    Tuple!(int, double)[] _products;
    double[] _nu;
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
         int[] reac_spidx, double[] reac_coeffs, int[] prod_spidx, double[] prod_coeffs,
         Tuple!(int,double)[] efficiencies, size_t n_species)
    {
        assert(reac_spidx.length == reac_coeffs.length,
               brokenPreCondition("reac_spidx and reac_coeffs arrays not of equal length"));
        assert(prod_spidx.length == prod_coeffs.length,
               brokenPreCondition("prod_spdix and prod_coeffs arrays are not of equal length"));
        super(forward, backward, gmodel);
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
            double nu1 = 0;
            foreach ( ref r; _reactants ) {
                if ( r[0] == isp ) {
                    nu1 = r[1];
                    break;
                }
            }
            double nu2 = 0;
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
         in Tuple!(int, double)[] reactants, in Tuple!(int, double)[] products, in double[] nu, in Tuple!(int, double)[] efficiencies)
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
    override void eval_equilibrium_constant(in GasState Q)
    {
        _Qw.T = Q.T;
	if (_gmodel.n_modes >= 1) _Qw.T_modes[] = Q.T; // equilibrium constant evaluated at thermal equilibrium 
        _K_eq = compute_equilibrium_constant(_gmodel, _Qw, _participants, _nu);
    }

    @nogc
    override number production(int isp) const
    {
        if ( _nu[isp] >  0.0 )
            return _nu[isp]*_w_f;
        else if ( _nu[isp] < 0.0 )
            return -_nu[isp]*_w_b;
        else
            return to!number(0.0);
    }

    @nogc
    override number loss(int isp) const
    {
        if ( _nu[isp] > 0.0 ) 
            return _nu[isp]*_w_b;
        else if ( _nu[isp] < 0.0 )
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
            double coeff = r[1];
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
    Tuple!(int, double)[] _reactants;
    Tuple!(int, double)[] _products;
    double[] _nu;
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
    lua_pop(L, 1);
    lua_getfield(L, -1, "brc");
    auto brc = createRateConstant(L, efficiencies, gmodel);
    lua_pop(L, 1);

    // And most use reacIdx, reacCoeffs, prodIdx and prodCoeffs lists.
    int[] reacIdx, prodIdx;
    double[] reacCoeffs, prodCoeffs;
    getArrayOfInts(L, -1, "reacIdx", reacIdx);
    getArrayOfDoubles(L, -1, "reacCoeffs", reacCoeffs);
    getArrayOfInts(L, -1, "prodIdx", prodIdx);
    getArrayOfDoubles(L, -1, "prodCoeffs", prodCoeffs);

    // We need to specialise the creation of a Reaction
    // based on type.
    auto type = getString(L, -1, "type");
    switch (type) {
    case "elementary":
        return new ElementaryReaction(frc, brc, gmodel, reacIdx, reacCoeffs,
                                      prodIdx, prodCoeffs, n_species);
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
        auto rc = new ArrheniusRateConstant(1.94e14, 0.0, 20620.0);
        auto gd = new GasState(3, 1);
        gd.T = 700.0;
        auto reaction = new ElementaryReaction(rc, rc, gm, [0, 1], [1.0, 1.0],
                                               [2], [2.0], 3);
        reaction.eval_rate_constants(gd);
        reaction.eval_rates(conc);
        assert(approxEqual(0.0, reaction.production(0)), failedUnitTest());
        assert(approxEqual(0.0, reaction.production(1)), failedUnitTest());
        assert(approxEqual(1287.8606, reaction.production(2), 1.0e-6), failedUnitTest());
        assert(approxEqual(643.9303, reaction.loss(0), 1.0e-6), failedUnitTest());
        assert(approxEqual(643.9303, reaction.loss(1), 1.0e-6), failedUnitTest());
        assert(approxEqual(0.0, reaction.loss(2)), failedUnitTest());

        auto reaction2 = new AnonymousColliderReaction(rc, rc, gm, [0, 1], [1., 1.],
                                                       [2], [2.], [tuple(1, 1.0)], 3);
        reaction2.eval_rate_constants(gd);
        reaction2.eval_rates(conc);
        
        // Try a reaction with backwards rate computed from equilibrium constant.
        auto reaction3 = new ElementaryReaction(rc, null, gm, [0, 1], [1., 1.],
                                                [2], [2.], 3);
        reaction3.eval_rate_constants(gd);
        reaction3.eval_rates(conc);

        return 0;
    }
}

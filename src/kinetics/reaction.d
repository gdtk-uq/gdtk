/**
 * reaction.d
 * This module holds the classes related to computing
 * rates of species change due to a single reaction.
 *
 * Author: Rowan J. Gollan
 */

module kinetics.reaction;

import std.math;
import std.string;
import std.typecons;
import std.algorithm;
import util.lua;
import util.lua_service;
import gas;
import util.msg_service;
import kinetics.rate_constant;
/++
 Reaction is an interface specifying the public services
 provided by an object of Reaction type.
+/
class Reaction
{
public:
    @property double k_f() const { return _k_f; }
    @property double k_b() const { return _k_b; }
    @property ref int[] participants() { return _participants; }
    
    this(in RateConstant forward, in RateConstant backward)
    {
	_forward = forward.dup();
	_backward = backward.dup();
    }

    abstract Reaction dup() const;
    
    final void eval_rate_constants(in GasState Q)
    {
	_k_f = eval_forward_rate_constant(Q);
	_k_b = eval_backward_rate_constant(Q);
    }

    final void eval_rates(in double[] conc)
    {
	_w_f = eval_forward_rate(conc);
	_w_b = eval_backward_rate(conc);
    }

    abstract double production(int isp) const;
    abstract double loss(int isp) const;

protected:
    abstract double eval_forward_rate(in double[] conc) const;
    abstract double eval_backward_rate(in double[] conc) const;
	    
private:
    RateConstant _forward, _backward;
    double _k_f, _k_b; // Storage of computed rate constants
    double _w_f, _w_b; // Storage of computed rates of change
    int [] _participants; // Storage of indices of species that
                          // participate in this reaction
    double eval_forward_rate_constant(in GasState Q) const
    {
	return _forward.eval(Q);
    }
    double eval_backward_rate_constant(in GasState Q) const
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
    this(in RateConstant forward, in RateConstant backward,
	 int[] reac_spidx, int[] reac_coeffs, int[] prod_spidx, int[] prod_coeffs,
	 size_t n_species)
    {
	assert(reac_spidx.length == reac_coeffs.length,
	       brokenPreCondition("reac_spidx and reac_coeffs arrays not of equal length"));
	assert(prod_spidx.length == prod_coeffs.length,
	       brokenPreCondition("prod_spdix and prod_coeffs arrays are not of equal length"));
	super(forward, backward);
	bool[int] pmap;
	foreach ( i; 0..reac_spidx.length ) {
	    _reactants ~= tuple(reac_spidx[i], reac_coeffs[i]);
	    pmap[reac_spidx[i]] = true;
	}
	foreach ( i; 0..prod_spidx.length ) {
	    _products ~= tuple(prod_spidx[i], prod_coeffs[i]);
	    pmap[prod_spidx[i]] = true;
	}
	_participants = pmap.keys.dup();
	sort(_participants);
	foreach ( isp; 0..n_species ) {
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
    this(in RateConstant forward, in RateConstant backward, in int[] participants,
	 in Tuple!(int, int)[] reactants, in Tuple!(int, int)[] products, in int[] nu)
    {
	super(forward, backward);
	_participants = participants.dup();
	_reactants = reactants.dup();
	_products = products.dup();
	_nu = nu.dup();
    }
    override ElementaryReaction dup() const
    {
	return new ElementaryReaction(_forward, _backward, _participants,
				      _reactants, _products, _nu);
    }
    override double production(int isp) const
    {
	if ( _nu[isp] >  0 )
	    return _nu[isp]*_w_f;
	else if ( _nu[isp] < 0 )
	    return -_nu[isp]*_w_b;
	else
	    return 0.0;
    }
    
    override double loss(int isp) const
    {
	if ( _nu[isp] > 0 ) 
	    return _nu[isp]*_w_b;
	else if ( _nu[isp] < 0 )
	    return -_nu[isp]*_w_f;
	else 
	    return 0.0;
    }
protected:
    override double eval_forward_rate(in double[] conc) const
    {
	double val = _k_f;
	foreach ( ref r; _reactants ) {
	    int isp = r[0];
	    int coeff = r[1];
	    val *= pow(conc[isp], coeff);
	}
	return val;
    }

    override double eval_backward_rate(in double[] conc) const
    {
	double val = _k_b;
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
 +          reacCoeffs = {1, 1},
 +          prodIdx = {2},
 +          prodCoeffs = {2},
 +          anonymousCollider = false,
 +          pressureDependent = false,
 +          efficiencies = {}
 + }
 +/
Reaction createReaction(lua_State* L, size_t n_species)
{
    // All Reactions have a forward and backward rate.
    lua_getfield(L, -1, "frc");
    auto frc = createRateConstant(L);
    lua_pop(L, 1);
    lua_getfield(L, -1, "brc");
    auto brc = createRateConstant(L);
    lua_pop(L, 1);

    // And most use reacIdx, reacCoeffs, prodIdx and prodCoeffs lists.
    int[] reacIdx, reacCoeffs, prodIdx, prodCoeffs;
    getArrayOfInts(L, -1, "reacIdx", reacIdx);
    getArrayOfInts(L, -1, "reacCoeffs", reacCoeffs);
    getArrayOfInts(L, -1, "prodIdx", prodIdx);
    getArrayOfInts(L, -1, "prodCoeffs", prodCoeffs);

    // We need to specialise the creation of a Reaction
    // based on type.
    auto type = getString(L, -1, "type");
    switch (type) {
    case "elementary":
	return new ElementaryReaction(frc, brc, reacIdx, reacCoeffs,
				      prodIdx, prodCoeffs, n_species);
    default:
	string msg = format("The reaction type: %s is not known.", type);
	throw new Exception(msg);
    }
}



unittest {
    // Find rate of forward production for H2 + I2 reaction at 700 K.
    double[] conc = [4.54, 4.54, 0.0];
    auto rc = new ArrheniusRateConstant(1.94e14, 0.0, 20620.0);
    auto gd = new GasState(3, 1);
    gd.T[0] = 700.0;
    auto reaction = new ElementaryReaction(rc, rc, [0, 1], [1, 1],
					   [2], [2], 3);
    reaction.eval_rate_constants(gd);
    reaction.eval_rates(conc);
    assert(approxEqual(0.0, reaction.production(0)), failedUnitTest());
    assert(approxEqual(0.0, reaction.production(1)), failedUnitTest());
    assert(approxEqual(1287.8606, reaction.production(2)), failedUnitTest());
    assert(approxEqual(643.9303, reaction.loss(0)), failedUnitTest());
    assert(approxEqual(643.9303, reaction.loss(1)), failedUnitTest());
    assert(approxEqual(0.0, reaction.loss(2)), failedUnitTest());
}

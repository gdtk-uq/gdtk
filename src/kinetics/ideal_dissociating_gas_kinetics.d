/**
 * ideal_dissociating_gas.d
 *
 * This is the Lighthill-Freeman ideal dissociating gas model as described in
 * M. N. Macrossan (1990)
 * Hypervelocity flow of dissociating nitrogen downstream of a blunt nose.
 * Journal of Fluid Mechanics Vol. 217, pp 167-202.
 *
 * Authors: Peter J., Katrina Sklavos and Rowan G.
 * Version: 2017-04-22: initial shell copied from powers-aslam-gas module.
 *          2017-04-11: Filled in IDG thermo details as noted in 
 *                      PJ's workbooks 2017-04-22 through 2017-05-11
 *          2017-07-27: Finished reactor function and eliminated bugs in workbook.
 */

module kinetics.ideal_dissociating_gas_kinetics;

import std.math : exp, fabs, fmax, fmin;

import gas;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;

final class UpdateIDG : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
	super(gmodel); // hang on to a reference to the gas model
	// We need to pick a number of pieces out of the gas-model file, again.
	// Although they exist in the GasModel object, they are private.
	auto L = init_lua_State();
	doLuaFile(L, fname);
	lua_getglobal(L, "IdealDissociatingGas");
	_W = getDouble(L, -1, "W"); // molecular weight, g/mole
	_Rnn = R_universal / _W * 1000; // gas constant for molecule, J/kg/K
	_T_d = getDouble(L, -1, "T_d"); // characteristic dissociation temperature, K
	_rho_d = getDouble(L, -1, "rho_d"); // characteristic density, g/cm^3
	// Rate constants follow.
	_C1 = getDouble(L, -1, "C1"); 
	_n1 = getDouble(L, -1, "n1");
	_C2 = getDouble(L, -1, "C2");
	_n2 = getDouble(L, -1, "n2");
	lua_pop(L, 1); // dispose of the table
	lua_close(L);
    }
    
    override void opCall(GasState Q, double tInterval, ref double dtSuggest,
			 ref double[] params)
    {
	// In the isolated reactor, density and internal energy remain constant.
	// Dissociation fraction and temperature may/will change.
	//
	// Initial state.
	double alpha0 = Q.massf[1];
	double T0 = Q.T;
	//
	// The IDG rate equations as functions.
	// See PJ workbook page 48, 23-Apr-2017
	double X(double alpha, double T)
	{
	    return (2.0*_C1*T^^_n1*alpha + _C2*T^^_n2*(1.0-alpha)) / _W;
	}
	double f1(double alpha, double T, double rho)
	// This is the time derivative d(alpha)/dt.
	{
	    return rho*X(alpha,T)*((1.0-alpha)*exp(-_T_d/T) - rho*alpha*alpha/_rho_d);
	}
	// The nonlinear constraint functions for the backward-Euler update.
	// See PJ's following workbook, page 4, 23-Apr-2017.
	double g1(double alpha1, double T1, double rho)
	{
	    return alpha1 - alpha0 - tInterval*f1(alpha1, T1, rho);
	}
	double f2(double alpha1, double T1, double u)
	{
	    return u - _Rnn*alpha1*_T_d - _Rnn*3.0*T1;
	}
	// Updated state, initially guessed as the same.
	double alpha1 = alpha0;
	double T1 = T0;
	foreach (i; 0 .. 30) {
	    double g1_nominal = g1(alpha1, T1, Q.rho);
	    double f2_nominal = f2(alpha1, T1, Q.u);
	    // Derivatives for the linear equations of the Newton update.
	    // See again PJ's workbook, page 4, 23-Apr-2017.
	    double df2da1 = -_Rnn*_T_d;
	    double df2dT1 = -3.0*_Rnn;
	    double del_alpha = (alpha1 > 0.5) ? -0.001: 0.001;
	    double dg1da1 = (g1(alpha1+del_alpha, T1, Q.rho) - g1_nominal)/del_alpha;
	    double del_T = fmax(0.01*T1, 10.0);
	    double dg1dT1 = (g1(alpha1, T1+del_T, Q.rho) - g1_nominal)/del_T;
	    // Solve the linear equations for the increments
	    double determinant = df2da1*dg1dT1 - dg1da1*df2dT1;
	    double dalpha1 = (-f2_nominal*dg1dT1 + g1_nominal*df2dT1)/determinant;
	    double dT1 = (-df2da1*g1_nominal + dg1da1*f2_nominal)/determinant;
	    alpha1 += dalpha1;
	    alpha1 = fmax(0.0, fmin(alpha1, 1.0));
	    T1 += dT1;
	    // Check to see if we are close enough.
	    if (fabs(dT1) < 0.01 && fabs(dalpha1) < 1.0e-6) { break; }
	} // end foreach
	Q.massf[0] = 1.0 - alpha1;
	Q.massf[1] = alpha1;
	// Since the internal energy and density in the (isolated) reactor is fixed,
	// we need to evaluate the new temperature, pressure, etc.
	_gmodel.update_thermo_from_rhou(Q);
	assert(fabs(Q.T - T1) < 0.1, "Oops, Temperature doesn't match after thermo update.");
	_gmodel.update_sound_speed(Q);
    } // end opCall

private:
    // Thermodynamic constants
    double _W; // Molecular mass in g/mole
    double _Rnn; // gas constant for molecule in J/kg/K
    double _T_d; // characteristic dissociation temperature, K
    double _rho_d; // characteristic density, g/cm^3
    // Rate constants
    double _C1, _n1, _C2, _n2;
} // end class UpdateIDG


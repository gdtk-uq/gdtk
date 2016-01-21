/* 
 * Testing area for some code

 */
module uniform_lut_demo;
import gas.uniform_lut;
import gas.gas_model;
import gas.gas_model_util;
import gas.ideal_gas;
import util.lua;
import util.lua_service;
import std.datetime;

import std.stdio;


void main() {
    writeln("Demo for uniform lut"); 
    writeln("Interpolating for the minimum energy/minimum density case");
    writeln("First row of data in the LUT");
    writeln("For cea-lut-airns.lu, this gives data: ");
    writeln("\t\t {cv_hat, cv, R_hat, cp_hat, gamma_hat, mu, k } =");
    writeln("\t\t {721.782, 742.464, 287.036, 1036.88, 1.38662, 2.7148e-05, 0.03885}");
    
    
    string file_name = "cea-lut-air.lua";    
    /*
    writeln("\n\nCreate gas model object and gas state object manually");
    lua_State* L;
    L = init_lua_State(file_name);
    auto gm1 = new UniformLUT(L);
    auto Q1 =  new GasState(gm1);
    // The following values are manually assigned to Q -they are the minimum values of 
    // emin and log(rho) in the look-up-table for air (first row of values)
    // Should lead to the following results:
    // {cv_hat, cv, R_hat, cp_hat, gamma_hat, mu, k } = 
    // {721.782, 742.464, 287.036, 1036.88, 1.38662, 2.7148e-05, 0.03885}
    Q1.e[0] = 360891;
    Q1.rho = 10.^^(-6);
    gm1.update_thermo_from_rhoe(Q1);
    gm1.update_trans_coeffs(Q1);
    writeln(Q1);
    

 
    
    writeln("Test running the model from rho, energy inputs");
    writeln("Run init_gas_model from gas_model_util.d");
    GasModel gm2 = init_gas_model(file_name);
    auto Q2 = new GasState(gm2);
    Q2.e[0] = 360891;
    Q2.rho = 10.^^(-6);
    gm2.update_thermo_from_rhoe(Q2);
    gm2.update_trans_coeffs(Q2);
    writeln(Q2);


    //writeln("\n\nTest of GasModel toString() method)"); 
    //writeln(gm2);
    
    
    writeln("\n\nTesting initialising the same state from density/temperature");
    GasModel gm3 = init_gas_model(file_name);
    auto Q3 = new GasState(gm3);
    Q3.rho = 10^^(-6.);
    Q3.T[0] = 500;
    gm3.update_thermo_from_rhoT(Q3);
    gm3.update_trans_coeffs(Q3);
    writeln(Q3);


    writeln("\n\nTesting initialising the same state from pressure/temperature");
    GasModel gm4 = init_gas_model(file_name);
    auto Q4 = new GasState(gm4);
    Q4.p = 0.143518;
    Q4.T[0] = 500;
    gm4.update_thermo_from_pT(Q4);
    gm4.update_trans_coeffs(Q4);
    writeln(Q4);
    */

    writeln("\n\n\n");
    GasModel gm5 = init_gas_model(file_name);
    auto Q5 = new GasState(gm5);
    // These values were the ones used in the version test in gas_model.d
    Q5.e = [3.0e6];
    Q5.rho = 1.0;
    gm5.update_thermo_from_rhoe(Q5);
    writeln(Q5);

    double dpdrho_const_T = gm5.dpdrho_const_T(Q5);
    writefln("For energy = %s and density = %s, dp_drho_const T is: ", Q5.e, Q5.rho);
    writeln(dpdrho_const_T);
    // double dpdrho_const_T = gm5.dpdrho_const_T(Q5);
    // writeln(dpdrho_const_T);
    // writeln("done");
    

    writeln("\n\nCompare dpdrho_const_T to ideal gas, air\n\n\n");
    GasModel gm6 = init_gas_model("sample-data/ideal-air-gas-model.lua");
    auto Q6 = new GasState(gm6);
    Q6.e = [3.0e6];
    Q6.rho = 1.0;
    gm6.update_thermo_from_rhoe(Q6);
    writeln(Q6);
    double dpdrho_const_T_ideal = gm6.dpdrho_const_T(Q6);
    writefln("For energy = %s and density = %s, dp_drho_const T is: ", Q6.e, Q6.rho);
    writeln(dpdrho_const_T_ideal);

    /*
    StopWatch sw;
    sw.start();
    int n = 10000;
    int i;
    for (i=0; i<n; ++i) {

	Q5.p = 0.143518;
	Q5.rho = 10^^(-6.);

	gm5.update_thermo_from_rhop(Q5);
	gm5.update_trans_coeffs(Q5);
    }
    sw.stop;
    long time = sw.peek().usecs;
    writeln(Q5);
    writefln("\n\nTime taken is %s microseconds", time);
    */  


}
  

/**
 * co2gas_sw.d
 * Span Wagner gas model for use in the CFD codes.
 *
 * Author: Jonathan H.
 * Version: 2015-09-17: initial cut, to explore options. Incorporated additional functions for extrapolation
 * PJ 2018-06-03
 * We make the signatures compatible with complex numbers but
 * we do not actually pass the complex values through to
 * Jonathan's functions.
 */

module gas.co2gas_sw;

import std.math;
import std.algorithm;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import util.msg_service;
import core.stdc.stdlib : exit;

import nm.complex;
import nm.number;
import nm.ridder;
import nm.bracketing;
import nm.tree_patch;
import nm.univariate_lut;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.sutherland_viscosity;
import gas.diffusion.sutherland_therm_cond;

class CO2GasSW: GasModel {
public:

    this() {
        // Default model is mostly initialized in the private data below.
        _n_species = 1;
        _n_modes = 0;
        _species_names ~= "CO2";
        _mol_masses ~= 0.04401121121333065;// value for CO2
        create_species_reverse_lookup();
        version(complex_numbers) {
            throw new Error("Do not use with complex numbers.");
        }
    }
        
    this(lua_State *L) {
        this();
        // Bring table to TOS
        lua_getglobal(L, "CO2GasSW");
        // Let's overwrite that here.
        _species_names[0] = getString(L, -1, "speciesName");
        // Now, pull out the remaining numeric value parameters.
        _mol_masses ~= getDouble(L, -1, "mMass");
        //_gamma = getDouble(L, -1, "gamma");//not used
        // Reference values for entropy
        lua_getfield(L, -1, "entropyRefValues");
        //_s1 = getDouble(L, -1, "s1");
        //_T1 = getDouble(L, -1, "T1");
        //_p1 = getDouble(L, -1, "p1");
        lua_pop(L, 1);
        // Molecular transport coefficent constants.
        lua_getfield(L, -1, "sutherlandVisc");
        _mu_ref = getDouble(L, -1, "mu_ref");
        _T_mu = getDouble(L, -1, "T_ref");
        _S_mu = getDouble(L, -1, "S");
        lua_pop(L, 1);
        lua_getfield(L, -1, "sutherlandThermCond");
        _k_ref = getDouble(L, -1, "k_ref");
        _T_k = getDouble(L, -1, "T_ref");
        _S_k = getDouble(L, -1, "S");
        lua_pop(L, 1);
        lua_getfield(L, -1, "LUTfilenames");
        string p_rhoe_filename = getString(L, -1, "p_rhoe_file");
        string a_rhoe_filename = getString(L, -1, "a_rhoe_file");
        string T_rhoe_filename = getString(L, -1, "T_rhoe_file");
        string e_rho_sat_table_filename = getString(L, -1, "e_rho_sat_file");
        string rho_sh_filename = getString(L, -1, "rho_sh_file");
        string T_sh_filename = getString(L, -1, "T_sh_file");
        lookup_hsFlag = getInt(L,-1, "lookup_hsFlag");
        lookup_rhoeFlag = getInt(L,-1, "lookup_rhoeFlag");
        lua_pop(L,1);
        if (lookup_rhoeFlag){
            P_rhoe_Tree = buildTree_fromFile(p_rhoe_filename);
            a_rhoe_Tree = buildTree_fromFile(a_rhoe_filename);
            T_rhoe_Tree = buildTree_fromFile(T_rhoe_filename);
            e_rho_sat_table = new uni_lut(e_rho_sat_table_filename);
        }
        if (lookup_hsFlag){
            rho_sh_Tree = buildTree_fromFile(rho_sh_filename);
            T_sh_Tree = buildTree_fromFile(T_sh_filename);
        }
        create_species_reverse_lookup();
    }

    this(string e_rho_sat_table_filename){
        //a special constructor used to make a gas model just with the sat-vap table
        //used to build_tree
        this.e_rho_sat_table = new uni_lut(e_rho_sat_table_filename);
        create_species_reverse_lookup();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "CO2 =(";
        repr ~= "name=\"" ~ _species_names[0] ~"\"";
        repr ~= ", Mmass=" ~ to!string(_mol_masses[0]);
        //repr ~= ", gamma=" ~ to!string(_gamma);
        //repr ~= ", s1=" ~ to!string(_s1);
        //repr ~= ", T1=" ~ to!string(_T1);
        //repr ~= ", p1=" ~ to!string(_p1);
        repr ~= ", mu_ref=" ~ to!string(_mu_ref);
        repr ~= ", T_mu=" ~ to!string(_T_mu);
        repr ~= ", S_mu=" ~ to!string(_S_mu);
        repr ~= ", k_ref=" ~ to!string(_k_ref);
        repr ~= ", T_mu= " ~ to!string(_T_k);
        repr ~= ", S_k=" ~ to!string(_S_k);
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
        debug {
            Q.rho = updateRho_PT(Q.p.re, Q.T.re);
            Q.u = updateEnergy_rhoT(Q.rho.re, Q.T.re);
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
        }
    }

    override void update_thermo_from_rhou(GasState Q) const
    {
        debug {
            if (lookup_rhoeFlag) {
                //following line assumes that both T and P trees constructed with same bounds
                double[2] uv = get_uv_rhoe(Q.rho.re, Q.u.re,
                                           T_rhoe_Tree.X_min,
                                           T_rhoe_Tree.X_max,
                                           T_rhoe_Tree.Y_min,
                                           T_rhoe_Tree.Y_max);
                Q.T = T_rhoe_Tree.search(uv[0], uv[1]).interpolateF(uv[0], uv[1]);
                Q.p = P_rhoe_Tree.search(uv[0], uv[1]).interpolateF(uv[0], uv[1]);
            } else {
                Q.T = updateTemperature_rhoe(Q.rho.re, Q.u.re);
                Q.p = updatePressure_rhoT(Q.rho.re,Q.T.re);
            }
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
        }
    }

    override void update_thermo_from_rhoT(GasState Q) const//DONE
    {
        debug {
            Q.p = updatePressure_rhoT(Q.rho.re, Q.T.re);
            Q.u = updateEnergy_rhoT(Q.rho.re, Q.T.re);
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
        }
    }

    override void update_thermo_from_rhop(GasState Q) const
    {
        debug {
            Q.T = updateT_Prho(Q.p.re, Q.rho.re);
            Q.u = updateEnergy_rhoT(Q.rho.re, Q.T.re);
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
        }
    }
    
    override void update_thermo_from_ps(GasState Q, number s) const
    {
        debug {
            Q.rho = getRho_EntropyP(s.re, Q.p.re, Q.T.re);//Q.T is modified by function
            Q.u = updateEnergy_rhoT(Q.rho.re, Q.T.re);
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
        }
    }

    override void update_thermo_from_hs(GasState Q, number h, number s) const
    {
        debug {
            if (lookup_hsFlag) {
                Q.rho = rho_sh_Tree.search(s.re,h.re).interpolateF(s.re,h.re);
                Q.T = T_sh_Tree.search(s.re,h.re).interpolateF(s.re,h.re);
            } else {
                Q.rho = getRho_sh(s.re, h.re, Q.T.re);
            }
            Q.u = updateEnergy_rhoT(Q.rho.re, Q.T.re);
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
        }
    }

    override void update_sound_speed(GasState Q) const
    {
        debug {
            if (lookup_rhoeFlag) {
                double[2] uv = get_uv_rhoe(Q.rho.re, Q.u.re,
                                           a_rhoe_Tree.X_min,
                                           a_rhoe_Tree.X_max,
                                           a_rhoe_Tree.Y_min,
                                           a_rhoe_Tree.Y_max);
                Q.a = a_rhoe_Tree.search(uv[0], uv[1]).interpolateF(uv[0], uv[1]);
            } else {
                Q.a = updateSoundSpeed_rhoT(Q.rho.re, Q.T.re);
            }
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
        }

    }
    override void update_trans_coeffs(GasState Q) const
    {
        Q.mu = sutherland_viscosity(Q.T, _T_mu, _mu_ref, _S_mu);
        Q.k = sutherland_thermal_conductivity(Q.T, _T_k, _k_ref, _S_k);
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
        throw new Exception("not implemented");
    }
    */
    override number dudT_const_v(in GasState Q) const
    {
        debug {
            return to!number(get_de_dT(Q.rho.re, Q.T.re));
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
            // return to!number(0.0);
        }
    }
    
    override number dhdT_const_p(in GasState Q) const
    {
        throw new Exception("Not implemented.");
    }
    
    override number dpdrho_const_T(in GasState Q) const
    {
        throw new Exception("Not implemented.");
    }
    
    override number gas_constant(in GasState Q) const
    {
        return to!number(_Rgas);
    }
    
    override number internal_energy(in GasState Q) const
    {
        return Q.u;
    }
    
    override number enthalpy(in GasState Q) const
    {
        debug {
            return to!number(updateEnthalpy_rhoT_original(Q.rho.re, Q.T.re));
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
            // return to!number(0.0);
        }
    }
    
    override number entropy(in GasState Q) const
    {
        debug {
            return to!number(updateEntropy_rhoT(Q.rho.re, Q.T.re));
        } else {
            assert(0, "Oops, not implemented for @nogc. PJ 2018-09-23");
            // return to!number(0.0);
        }
    }
    //------A function that re-maps the rho, T domain according to the liquid-vapour line------------
    //placed in public so it is available for building Tables
    /*
        u = 0 ---------------------------u = 0.5 -------------------------u = 1.0       v = 1.0
          |                                 |                                |
          |                                 |                                |
          |                                 |                                |
          |                                 |                                |
          |       Liquid-Vapour             |             Gas                |          rho
          |                                 |                                |
          |                                 |                                |
          |                                 |                                |
          |                                 |                                |  
        u = 0 ---------------------------u = 0.5 -------------------------u = 1.0, v = 0.0
                                                                                T
    */
    const double[2] get_uv_rhoT(double rho, double T, double rho_min, double rho_max, double T_min, double T_max)
    {
        double v = (rho - rho_min)/(rho_max - rho_min);
        double T_sat;
        if (rho > _rhoc) {
            T_sat = getT_satliq(rho);
        }
        else {
            T_sat = getT_satvap(rho);
        }
        T_sat = max(T_sat,T_min); //added to prevent mapping outside the domain if an different T_min is set
        double u;
        if (T > T_sat) {
            u = (T - T_sat)/(T_max - T_sat)*0.5 + 0.5;
        }
        else {
            u = (T - T_min)/(T_sat - T_min)*0.5;
        }
        return [u,v];
    } 

    const double[2] get_rhoT_uv(double u, double v, double rho_min, double rho_max, double T_min, double T_max)
    {
        double rho = v*rho_max + (1-v)*rho_min;
        double T_sat;
        if (rho > _rhoc) {
            T_sat = getT_satliq(rho);
        }
        else {
            T_sat = getT_satvap(rho);
        }
        T_sat = max(T_sat,T_min); //added to prevent mapping outside the domain if an different T_min is set
        double T;
        if (u > 0.5) {  
            u = 2*u - 1;
            T = u*T_max + (1 - u)*T_sat;
        }
        else {  
            u = 2*u;
            T = u*T_sat + (1 - u)*T_min;
        }
        return [rho, T];
    }

    const double[2] get_uv_rhoe(double rho, double e, double rho_min, double rho_max, double e_min, double e_max)
    {
        double v = (rho - rho_min)/(rho_max - rho_min);
        double e_sat;
        e_sat = e_rho_sat_table.quadlookup(rho);
        //e_sat = get_esat_rho(rho);
        e_sat = max(e_sat,e_min); //added to prevent mapping outside the domain if an different T_min is set
        double u;
        if (e > e_sat) {
            u = (e - e_sat)/(e_max - e_sat)*0.5 + 0.5;
        }
        else {
            u = (e - e_min)/(e_sat - e_min)*0.5;
        }
        return [u,v];
    } 

    const double[2] get_rhoe_uv(double u, double v, double rho_min, double rho_max, double e_min, double e_max)
    {
        double rho = v*rho_max + (1-v)*rho_min;
        double e_sat;
        e_sat = e_rho_sat_table.quadlookup(rho);
        //e_sat = get_esat_rho(rho);
        e_sat = max(e_sat,e_min); //added to prevent mapping outside the domain if an different e_min is set too high
        double e;
        if (u > 0.5) {  
            u = 2*u - 1;
            e = u*e_max + (1 - u)*e_sat;
        }
        else {  
            u = 2*u;
            e = u*e_sat + (1 - u)*e_min;
        }
        return [rho, e];
    }
    
    const double get_esat_rho(double rho) {
        if (rho > _rhoc) {
            return updateEnergy_rhoT(rho, getT_satliq(rho));
        }
        else { 
            return updateEnergy_rhoT(rho, getT_satvap(rho));
        }
    }

private:
    // Thermodynamic constants
        double _Rgas = 188.9241;
        double _u0 = 3.2174105e5-74841;
        double _s0 = 2.1396056e2;
        double _M = 44.01;//kg/kmol
        double _Tc = 304.1282;//K
        double _Pc = 7.3773e6;//Pc MPa
        double _rhoc = 467.6;//kg/m^304
        double _T0 = 298.15;//
        double[] _n = [0,
                         0.38856823203161e0, 0.29385475942740e1, -0.55867188534934e1, -0.76753199592477e0, 0.31729005580416e0,
                         0.54803315897767e0, 0.12279411220335e0, 0.21658961543220e1, 0.15841735109724e1, -0.23132705405503e0,
                         0.58116916431436e-1, -0.55369137205382e0, 0.48946615909422e0, -0.24275739843501e-1, 0.62494790501678e-1,
                         -0.12175860225246e0, -0.37055685270086e0, -0.16775879700426e-1, -0.11960736637987e0, -0.45619362508778e-1,
                         0.35612789270346e-1, -0.74427727132052e-2, -0.17395704902432e-2, -0.21810121289527e-1, 0.24332166559236e-1,
                         -0.37440133423463e-1, 0.14338715756878e0, -0.13491969083286e0, -0.23151225053480e-1, 0.12363125492901e-1,
                         0.21058321972940e-2, -0.33958519026368e-3, 0.55993651771592e-2, -0.30335118055646e-3, -0.21365488688320e3, 
                         0.26641569149272e5, -0.24027212204557e5, -0.28341603423999e3, 0.21247284400179e3, -0.66642276540751e0,
                         0.72608632349897e0, 0.5506866862842e-1];
        double[] _d = [0, 1,1,1,1,2,2,3,1,2,4,5,5,5,6,6,6,1,1,4,4,4,7,8,2,3,3,5,5,6,7,8,10,4,8,2,2,2,3,3];
        double[] _t = [0, 0, 0.75, 1, 2, 0.75, 2, 0.75, 1.5, 1.5, 2.5, 0, 1.5, 2, 0, 1, 2, 3, 6, 3, 6, 8, 6, 0, 7, 12, 16, 22, 24, 16, 24, 8, 2, 28, 14, 1, 0, 1, 3, 3];
        double[] _c = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6];
        double[] _alpha = [25, 25, 25, 15, 20];
        double[] _beta = [325, 300, 300, 275, 275];
        double[] _gamma = [1.16, 1.19, 1.19, 1.25, 1.22];
        double[] _epsilon =[1, 1, 1, 1, 1];
        double[] _a = [3.5, 3.5, 3.0];
        double[] _b = [0.875, 0.925, 0.875];
        double[] _beta2 = [0.3, 0.3, 0.3];
        double[] _A = [0.7, 0.7, 0.7];
        double[] _B = [0.3, 0.3, 1.0];
        double[] _C = [10.0, 10.0, 12.5];
        double[] _D = [275, 275, 275];
        double[] _a0 = [0.0, 8.37304456, -3.70454304, 2.5, 1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678];
        double[] _theta0= [0.0, 0.0, 0.0, 0.0, 3.15163, 6.11190, 6.77708, 11.32384, 27.08792];
        double[] _G = [0.0, 8.726361e3, 1.840040e2, 1.914025, -1.667825e-3, 7.305950e-7, -1.255290e-10];
        //Saturation Pressure and Density information
        double[] _a_Psat = [0, -7.0602087, 1.9391218, -1.6463597, -3.2995634];
        double[]_t_Psat = [0, 1.0, 1.5, 2.0, 4.0];
        double[]_a_rhosatvap = [0, -1.7074879, -0.82274670, -4.6008549, -10.111178, -29.742252];
        double[]_t_rhosatvap = [0, 0.340, 0.5, 1.0, 7.0/3.0, 14.0/3.0];
        double[]_a_rhosatliq = [0.0, 1.9245108, -0.62385555, -0.32731127, 0.39245142];
        double[]_t_rhosatliq = [0.0, 0.34, 0.5, 10.0/6.0, 11/6.0];
        // Molecular transport coefficent constants.
        double _mu_ref = 1.716e-5; // Pa.s
        double _T_mu = 273.0; // degrees K
        double _S_mu = 111.0; // degrees K
        double _k_ref = 0.0241; // W/(m.K) 
        double _T_k = 273.0; // degrees K
        double _S_k = 194.0; // degrees K
        
        Tree T_rhoe_Tree;
        Tree P_rhoe_Tree;
        Tree a_rhoe_Tree;
        Tree rho_sh_Tree;
        Tree T_sh_Tree;
        int lookup_hsFlag = 0;
        int lookup_rhoeFlag = 0;
        uni_lut e_rho_sat_table;
    
    //define some functions that will be available to the functions in the private area  
    const double extrapolate_rhoT_function(alias f)(double rho, double T,string extrapType = "constant")
        if (is(typeof(f(400.0, 400.0))==double))
    {
                //extrapoaltes a function mapped on rho, T
        double rho_satvap;
        double rho_satliq;
        if (T < _Tc) {
            rho_satvap = getrho_satvap(T);//be careful here about putting in supercritical fluids
            rho_satliq = getrho_satliq(T);
            if ( (rho <= rho_satliq) && (rho >= rho_satvap) ) {
                double T_sat;
                if (rho > _rhoc) {
                    T_sat = getT_satliq(rho);
                }
                else { 
                    T_sat = getT_satvap(rho);
                }
                assert(!isNaN(T_sat),format("T_sat is NaN, at rho = %f",rho));
                if ( (T_sat == _Tc) && (rho == _rhoc) ) {
                    T_sat += 0.00001; //this is to fudge it in case it lands on the critical point
                }
                double F_sat = f(rho, T_sat);
                switch (extrapType){
                case "constant":
                    return F_sat;
                case "linear":
                    double deltaT = 0.01;
                    double df_dT = (f(rho, T_sat+deltaT) - f(rho, T_sat))/deltaT;
                    return  F_sat+(T - T_sat)*df_dT;
                case "quality":
                    double v_f = 1.0/rho_satliq;
                    double v_g = 1.0/rho_satvap;
                    double v = 1/rho;
                    double qual = (v-v_f)/(v_g-v_f);
                    double f_f = f(rho_satliq, T);
                    double f_g = f(rho_satvap, T);
                    return qual*f_g + (1 - qual)*f_f;
                case "pressure":
                    //this is based on the idea that pressure is constant at a certain temperature - it goes to the satvap curve
                    return f(rho_satvap,T);
                case "soundspeed":
                    double v_f = 1.0/rho_satliq;
                    double v_g = 1.0/rho_satvap;
                    double v = 1/rho;
                    double qual = (v-v_f)/(v_g-v_f);
                    double f_f = f(rho_satliq, T);
                    double f_g = f(rho_satvap, T);
                    return pow(qual*(qual*pow(f_g,-2) + (1 - qual)*pow(f_f,-2)) + (1 - qual)*(1 - qual)*pow(f_f,-2),-0.5);
                case "zero":
                    return 0;
                default:
                    throw new Exception(format("extrapType %s is an invalid choice", extrapType));
                }
            }
        }
        return f(rho, T);
    }
    
    const double updateQuality_rhoT(double rho, double T)
    {
        double rho_satvap;
        double rho_satliq;
        if (T < _Tc) {
            rho_satvap = getrho_satvap(T);
            rho_satliq = getrho_satliq(T);
            if((rho <= rho_satliq)&&(rho >= rho_satvap)) {
                double v_f = 1.0/rho_satliq;
                double v_g = 1.0/rho_satvap;
                double v = 1/rho;
                return (v-v_f)/(v_g-v_f);
            }
        }
        return 1.0;
    }

    const double updatePressure_rhoT(double rho, double T)
    {
        string extrapType = "pressure";
        return extrapolate_rhoT_function!updatePressure_rhoT_original(rho, T,extrapType);
        
    }

    const double updatePressure_rhoT_original(double rho, double T)
    {
        //All the meanings of all the variables etc. names can be found in 
        //Span and Wagner see p37
        double Z = 1 + rho/_rhoc*get_d_alphar_d_delta(_Tc/T,rho/_rhoc);
        assert(!isNaN(Z), format("Pressure is NaN, rho: %s, T: %s", rho, T));
        return Z*rho*_Rgas*T;
    }

    const double updateEnergy_rhoT(double rho, double T)
    {
        string extrapType = "quality";
        return extrapolate_rhoT_function!updateEnergy_rhoT_original(rho, T,extrapType);
    }

    const double updateEnthalpy_rhoT(double rho, double T)
    {
        string extrapType = "quality";
        return extrapolate_rhoT_function!updateEnthalpy_rhoT_original(rho, T, extrapType);
    }

    const double updateEnthalpy_rhoT_original(double rho, double T)
    {
        return updateEnergy_rhoT_original(rho, T) + updatePressure_rhoT_original(rho, T)/rho;
    }                                                                   

    const double updateEnergy_rhoT_original(double rho, double T)
    {
        //Derivatives of Helmholtz energy w.r.t. reduced temperature give internal energy
        //See p867 of Sengers et al. EQUATIONS of STATE for FLUIDS and FLUID MIXTURES
        double d_alphar_d_tau = 0;
        double delta = rho/_rhoc;
        double tau = _Tc/T;
        for(int i = 1; i != 8; i++) {
            d_alphar_d_tau += _n[i]*_t[i]*delta^^_d[i]*tau^^(_t[i]-1);
        }
        for(int i = 8; i != 35; i++) {
            d_alphar_d_tau += _n[i]*_t[i]*delta^^_d[i]*tau^^(_t[i]-1)*exp(-delta^^_c[i]);
        }
        int j;
        for(int i = 35; i != 40; i++) {
            j = i - 35;
            d_alphar_d_tau += _n[i]*delta^^_d[i]*tau^^_t[i]*
                exp(-_alpha[j]*(delta-_epsilon[j])^^2 - _beta[j]*(tau - _gamma[j])^^2)*
                (_t[i]/tau - 2*_beta[j]*(tau-_gamma[j]));
        }
        double Psi;
        double theta;
        double Delta;
        double dPsi_dtau;
        double dDeltab_dtau;
        for (int i = 40; i != 43; i++) {
            j = i - 40;
            Psi = exp(-_C[j]*(delta - 1)^^2 - _D[j]*(tau - 1)^^2);
            theta = (1 - tau) + _A[j]*((delta-1)^^2)^^(0.5/_beta2[j]);
            Delta = theta*theta + _B[j]*((delta - 1)^^2)^^_a[j];
            dPsi_dtau = -2*_D[j]*(tau - 1)*Psi;
            dDeltab_dtau = -2*theta*_b[j]*Delta^^(_b[j] - 1);
            d_alphar_d_tau += _n[i]*delta*(dDeltab_dtau*Psi + Delta^^_b[j]*dPsi_dtau);
        }
        //----------------------------------------------------
        //Add on the ideal part of the equation
        //Coefficients adjusted to give zero for the ideal gas enthalpy at T0 = 298.15K
        //------------------------------------------------------
        double d_alpha0_d_tau = _a0[2] + _a0[3]/tau;
        for (int i = 4; i != 9; i++) {
            d_alpha0_d_tau += _a0[i]*_theta0[i]*(1.0/(1.0-exp(-_theta0[i]*tau)) -1);
        }
        return _Tc*(d_alpha0_d_tau + d_alphar_d_tau)*_Rgas; 
    }
   
    const double updateTemperature_rhoe_newtonCotes(double rho, double e, int maxIterations = 100, double Ttol = 0.1)
    {
        double T = (e +383334.16+191.86*rho)/1030.1886; // using a dodgy linear thing
        T = max(310, T); //don't go too low; // originally worked well with T = 270 - but not so well once i included the bracketing
        for (int i = 0; i != maxIterations; i++) {
            double deltaT = (updateEnergy_rhoT(rho,T)-e)/get_de_dT(rho,T);
            if (abs(deltaT) < Ttol) {
                break;
            }
            T -= deltaT;
        }
        return T;
    }
   
    const double updateTemperature_rhoe(double rho, double e,double[2] Tbracket = [80, 2000], double Ttol = 0.001)
    {
        //Uses Ridders Method
        auto zeroSpecificEnergy_T = delegate(double T)
            {
                return e - updateEnergy_rhoT(rho, T);
            };
        int max_try = 50;
        double factor = 0.05;
        double T1 = 400; double T2 = 600; 
        double T1_min = 80; double T2_max = 50000.0;
        int bracketFlag = bracket!(zeroSpecificEnergy_T,double)(T1,T2,T1_min,T2_max,max_try,factor);
        if (bracketFlag == -1) throw new Exception(format("bracketing getTemperature failed at e: %s, rho: %s, bracket: [%s, %s]", e, rho, T1, T2));
        double T = solve!(zeroSpecificEnergy_T,double)(T1,T2, Ttol);
        assert(!isNaN(T), format("T is NaN at rho: %s, e: %s", rho, e));
        return T;
    }

    const double get_de_dT(double rho, double T)
    {
        //Gets derivative of specific energy w.r.t. Temperature for evaluating T as a function of e, rho based on Newton's method
        //conveniently de_dT is c_V
        double tau = _Tc/T;
        double delta = rho/_rhoc;
        //---------------------------------------d2_alphar_d_tau2
        double d2_alphar_d_tau2 = 0;
        for(int i = 1; i != 8; i++) {
            d2_alphar_d_tau2 += _n[i]*_t[i]*(_t[i]-1)*delta^^_d[i]*tau^^(_t[i]-2);
        }
        for(int i = 8; i != 35; i++) {
            d2_alphar_d_tau2 += _n[i]*_t[i]*(_t[i]-1)*delta^^_d[i]*tau^^(_t[i]-2)*exp(-delta^^_c[i]);
        }
        
        int j;
        for(int i = 35; i != 40; i++) {
            j = i - 35;
            d2_alphar_d_tau2 += _n[i]*delta^^_d[i]*tau^^_t[i]*
                exp(-_alpha[j]*(delta-_epsilon[j])^^2 - _beta[j]*(tau - _gamma[j])^^2)*
                ((_t[i]/tau - 2*_beta[j]*(tau-_gamma[j]))^^2-_t[i]/tau/tau-2*_beta[j]);
        }       
        double Psi;
        double theta;
        double Delta;
        double dPsi_dtau;
        double dDeltab_dtau;
        double d2Psi_dtau2;
        double d2Deltab_dtau2;
        for (int i = 40; i != 43; i++) {
            j = i - 40;
            Psi = exp(-_C[j]*(delta - 1)^^2 - _D[j]*(tau - 1)^^2);
            theta = (1 - tau) + _A[j]*((delta-1)^^2)^^(0.5/_beta2[j]);
            Delta = theta*theta + _B[j]*((delta - 1)^^2)^^_a[j];
            dPsi_dtau = -2*_D[j]*(tau - 1)*Psi;
            dDeltab_dtau = -2*theta*_b[j]*Delta^^(_b[j] - 1);
            d2Psi_dtau2 = (2*_D[j]*(tau - 1)*(tau - 1) - 1)*2*_D[j]*Psi;
            d2Deltab_dtau2 = 2*_b[j]*Delta^^(_b[j]-1) + 4*theta*theta*_b[j]*(_b[j]-1)*Delta^^(_b[j]-2);
            d2_alphar_d_tau2 += _n[i]*delta*(d2Deltab_dtau2*Psi + 2*dDeltab_dtau*dPsi_dtau + Delta^^_b[j]*d2Psi_dtau2);
        }
        //----------------------------------------------------
        //Add on the ideal part of the equation
        //Copefficients adjusted to give zero for the ideal gas enthalpy at T0 = 298.15K
        //------------------------------------------------------
        double d2_alpha0_d_tau2 = -_a0[3]/tau/tau;
        for (int i = 4; i != 9; i++) {
                d2_alpha0_d_tau2 -= _a0[i]*_theta0[i]*_theta0[i]*exp(-_theta0[i]*tau)*(1.0-exp(-_theta0[i]*tau))^^-2;
        }
        return -tau*tau*(d2_alphar_d_tau2 + d2_alpha0_d_tau2)*_Rgas;
    }

    const double updateT_Prho(double P, double rho, double[2] Tbracket = [150, 2000], double tol = 0.001)
    {
        //will get Temperature based on Ridder's Method - see PJ's notes
        auto zeroPressure_T = delegate (double T)
            {
                return P - updatePressure_rhoT(rho, T);
            };
        double T = solve!(zeroPressure_T,double)(Tbracket[0], Tbracket[1], tol);
        assert(!isNaN(T), format("T is NaN at P: %s, rho: %s", P, rho));
        return T;
    }
   
    const double updateRho_PT(double P, double T, double[2] bracket = [0.0001, 2000], double tol = 1e-6)
    {
        auto zeroPressure_rho = delegate (double rho)
            {
            return P - updatePressure_rhoT(rho, T);
            };
        double rho = solve!(zeroPressure_rho,double)(bracket[0], bracket[1], tol);
        assert(!isNaN(rho), format("rho is NaN at P: %s, T: %s"));
        return rho;
    }      

    const double updateSoundSpeed_rhoT(double rho, double T)
    {
        string extrapType = "soundspeed";
        return extrapolate_rhoT_function!updateSoundSpeed_rhoT_original(rho, T, extrapType);
    }

    const double updateSoundSpeed_rhoT_original(double rho, double T)
    {
        double tau = _Tc/T;
        double delta = rho/_rhoc;
        double d_alphar_d_delta = get_d_alphar_d_delta(tau,delta);
        
        double d2_alphar_d_delta2 = 0;
        double d2_alphar_d_delta_d_tau = 0;
        for (int i = 1; i != 8; i++) {
            d2_alphar_d_delta2 += _n[i]*_d[i]*(_d[i]-1)*delta^^(_d[i]-2)*tau^^_t[i];
                
            d2_alphar_d_delta_d_tau += _n[i]*_d[i]*_t[i]*delta^^(_d[i]-1)*tau^^(_t[i]-1);
        }
        for (int i = 8; i != 35; i++) {
            d2_alphar_d_delta2 += _n[i]*exp(-delta^^_c[i])*
                (delta^^(_d[i]-2)*tau^^_t[i]*((_d[i]-_c[i]*delta^^_c[i])*(_d[i]-1.0-_c[i]*delta^^_c[i])-_c[i]^^2*delta^^_c[i]));
                
            d2_alphar_d_delta_d_tau += _n[i]*exp(-delta^^_c[i])*_t[i]*tau^^(_t[i] - 1.0)*
                delta^^(_d[i]-1)*(_d[i]-_c[i]*delta^^_c[i]);
        }
        int j;
        for (int i = 35; i != 40; i++) { 
            j = i - 35;
            d2_alphar_d_delta2 += _n[i]*tau^^_t[i]*
                exp(-_alpha[j]*(delta-_epsilon[j])^^2 - _beta[j]*(tau - _gamma[j])^^2)*
                (-2*_alpha[j]*delta^^_d[i]
                 + 4*_alpha[j]^^2*delta^^_d[i]*(delta-_epsilon[j])^^2
                 - 4*_d[i]*_alpha[j]*delta^^(_d[i]-1)*(delta-_epsilon[j])
                 + _d[i]*(_d[i]-1)*delta^^(_d[i]-2));
                
                d2_alphar_d_delta_d_tau += _n[i]*delta^^_d[i]*tau^^_t[i]*
                    exp(-_alpha[j]*(delta-_epsilon[j])^^2 - _beta[j]*(tau - _gamma[j])^^2)*
                    (_d[i]/delta - 2*_alpha[j]*(delta - _epsilon[j]))*
                    (_t[i]/tau - 2*_beta[j]*(tau - _gamma[j]));
        }
        double Psi;
        double theta;
        double Delta;
        double dPsi_ddelta;
        double dDelta_ddelta;
        double dDeltab_ddelta;
        double dPsi_dtau;
        double dDeltab_dtau;
        double d2Delta_ddelta2;
        double d2Deltab_ddelta2;
        double d2Psi_ddelta2;
        double d2Psi_ddelta_dtau;
        double d2Deltab_ddelta_dtau;
        for (int i = 40; i != 43; i++) {
            j = i - 40;
            Psi = exp(-_C[j]*(delta - 1)^^2 - _D[j]*(tau - 1)^^2);
            theta = (1 - tau) + _A[j]*((delta-1)^^2)^^(0.5/_beta2[j]);
            Delta = theta*theta + _B[j]*((delta - 1)^^2)^^_a[j];
            dPsi_ddelta = -2*_C[j]*(delta - 1)*Psi;
            dDelta_ddelta = (delta - 1)*(2*theta*_A[j]/_beta2[j]*((delta - 1)^^2)^^(0.5/_beta2[j] - 1) + 2*_B[j]*_a[j]*((delta - 1)^^2)^^(_a[j] - 1));//potential to simplify the process for d2Delta_ddelta2 -- a part was factorised out here
            dDeltab_ddelta = _b[j]*Delta^^(_b[j] - 1)*dDelta_ddelta;
            dPsi_dtau = -2*_D[j]*(tau - 1)*Psi;
            dDeltab_dtau = -2*theta*_b[j]*Delta^^(_b[j] - 1);
            d2Psi_ddelta2 = (2*_C[j]*(delta-1)^^2-1)*2*_C[j]*Psi;
            d2Deltab_ddelta_dtau = -_A[j]*_b[j]*2/_beta2[j]*Delta^^(_b[j] - 1)*(delta - 1)*((delta - 1)^^2)^^(0.5/_beta2[j] - 1)
                - 2*theta*_b[j]*(_b[j]-1)*Delta^^(_b[j] - 2)*dDelta_ddelta;
            d2Delta_ddelta2 = 2*theta*_A[j]/_beta2[j]*((delta - 1)^^2)^^(0.5/_beta2[j] - 1) + 2*_B[j]*_a[j]*((delta - 1)^^2)^^(_a[j] - 1);
            d2Delta_ddelta2 += (delta-1.0)^^2*(4*_B[j]*_a[j]*(_a[j]-1)*((delta-1)^^2)^^(_a[j]-2)
                                               +2*_A[j]^^2*_beta2[j]^^-2*(((delta-1)^^2)^^(0.5/_beta2[j]-1))^^2)
                +_A[j]*theta*4/_beta2[j]*(0.5/_beta2[j]-1)*((delta-1)^^2)^^(0.5/_beta2[j]-1);
            //expanded into the last part to prevent a 0^^negative number when delta = 1.0
            d2Deltab_ddelta2 = _b[j]*(Delta^^(_b[j] - 1)*d2Delta_ddelta2 + (_b[j] - 1)*Delta^^(_b[j]-2)*(dDelta_ddelta)^^2);
            d2Psi_ddelta_dtau = 4*_C[j]*_D[j]*(delta - 1)*(tau - 1)*Psi;
            d2_alphar_d_delta2 += _n[i]*(Delta^^_b[j]*(2*dPsi_ddelta + delta*d2Psi_ddelta2)
                                         +2*dDeltab_ddelta*(Psi + delta*dPsi_ddelta)
                                         +d2Deltab_ddelta2*delta*Psi);
            d2_alphar_d_delta_d_tau += _n[i]*(Delta^^_b[j]*(dPsi_dtau + delta*d2Psi_ddelta_dtau) 
                                              + delta*dDeltab_ddelta*dPsi_dtau 
                                              + dDeltab_dtau*(Psi + delta*dPsi_ddelta) 
                                              + d2Deltab_ddelta_dtau*delta*Psi);
        }
                
        double w2_RT =  1+2*delta*d_alphar_d_delta + delta^^2*d2_alphar_d_delta2
            -(1 + delta*d_alphar_d_delta - delta*tau*d2_alphar_d_delta_d_tau)^^2
            /(-get_de_dT(rho,T)/_Rgas);//note de_dT is equivalent to c_V
        
        return sqrt(w2_RT * _Rgas * T);
    }

    const double updateEntropy_rhoT(double rho, double T)
    {
        string extrapType = "quality";
        return extrapolate_rhoT_function!updateEntropy_rhoT_original(rho, T,extrapType);
    }

    const double updateEntropy_rhoT_original(double rho, double T) {
        double delta = rho/_rhoc;
        double tau = _Tc/T;
        double alpha0 = getalpha0(delta,tau);
        double alphaR = getalphaR(delta,tau);   
        double d_alpha0_d_tau = _a0[2] + _a0[3]/tau;
        for (int i = 4; i != 9; i++) {
            d_alpha0_d_tau += _a0[i]*_theta0[i]*(1.0/(1.0-exp(-_theta0[i]*tau)) -1);
        }
        return updateEnergy_rhoT(rho, T)/T - (alpha0 + alphaR)*_Rgas;
    }
    
    const double getalpha0(double delta, double tau)
    {
        double alpha0 = log(delta) + _a0[1] + _a0[2]*tau + _a0[3]*log(tau);//ideal part of helmholtz energy
        for (int i = 4; i != 9; i++) {
                alpha0 += _a0[i]*log(1 - exp(-_theta0[i]*tau));
        }
        return alpha0;
    }

    const double getalphaR(double delta, double tau) {
        double alphaR = 0;
        for (int i = 1; i != 8; i++) {
                alphaR += _n[i]*delta^^_d[i]*tau^^_t[i];
        }
        for (int i = 8; i != 35; i++) {
            alphaR += _n[i]*delta^^_d[i]*tau^^_t[i]*exp(-delta^^_c[i]);
        }
        int j;
        for (int i = 35; i != 40; i++) {
            j = i - 35;
            alphaR += _n[i]*delta^^_d[i]*tau^^_t[i]*
                exp(-_alpha[j]*(delta-_epsilon[j])^^2 - _beta[j]*(tau - _gamma[j])^^2);
                                                        
        }
        double Psi;
        double theta;
        double Delta;
        for (int i = 40; i != 43; i++) {
            j = i - 40;
            Psi = exp(-_C[j]*(delta - 1)^^2 - _D[j]*(tau - 1)^^2);
            theta = (1 - tau) + _A[j]*((delta-1)^^2)^^(0.5/_beta2[j]);
            Delta = theta*theta + _B[j]*((delta - 1)^^2)^^_a[j];
            alphaR += _n[i]*Delta^^_b[j]*delta*Psi;
        }
        return alphaR;
    }

    const double getT_p_rho(double p, double rho, double tol = 0.001) {
        auto zeroPressure_T = delegate (double T)
            {
                return p - updatePressure_rhoT(rho, T);
            };
        double T_max;
        if (rho < 10) T_max = 1.1*p/rho/_Rgas; //if low density, Temperature may be high -- need to extend bracket
        else T_max = 6000;
        double T = solve!(zeroPressure_T,double)(80,max(T_max,6000), tol);
        assert(!isNaN(T), format("T is NaN, P: %s, rho: %s", p, rho));
        return T;
    } 
    
    const double getEntropy_prho(double p, double rho, ref double T){
        T = getT_p_rho(p,rho);
        return updateEntropy_rhoT(rho, T);
    }

    const double getRho_EntropyP(double s, double P, ref double T){
        auto zeroEntropy_rho = delegate (double rho)
            {
                return s - getEntropy_prho(P, rho,T);
            };
        return solve!(zeroEntropy_rho,double)(0.009,1600,0.0001);
    }


    //-----------FUNCTIONS TO UPDATE FROM H, S
    const double getT_hrho(double h, double rho, double tol = 0.001) 
    {
        auto zeroEnthalpy_T = delegate (double T)
            {
                return updateEnthalpy_rhoT(rho, T) - h;
            };
        return solve!(zeroEnthalpy_T,double)(80,4000,tol);
    }
    
    const double getEntropy_hrho(double h, double rho, ref double T)
    {
        T = getT_hrho(h,rho);
        return updateEntropy_rhoT(rho,T);
    }

    const double getRho_sh(double s, double h, ref double T) {
        auto zeroEntropy_rho = delegate(double rho){
            return getEntropy_hrho(h,rho,T)-s;
        };
        return solve!(zeroEntropy_rho,double)(0.0009,1600,0.0001);
    }

    const double get_d_alphar_d_delta(double tau, double delta) {
        double d_alphar_d_delta = 0;
        for (int i = 1; i != 8; i++) {
            d_alphar_d_delta += _n[i]*_d[i]*delta^^(_d[i]-1)*tau^^_t[i];
        }
        for (int i = 8; i != 35; i++) {
            d_alphar_d_delta += _n[i]*exp(-delta^^_c[i])*(delta^^(_d[i]-1)*tau^^_t[i]*(_d[i]-_c[i]*delta^^_c[i]));
        }
        int j;
        for (int i = 35; i != 40; i++) {
            j = i - 35;
            d_alphar_d_delta += _n[i]*delta^^_d[i]*tau^^_t[i]*
                exp(-_alpha[j]*(delta-_epsilon[j])^^2 - _beta[j]*(tau - _gamma[j])^^2)*
                (_d[i]/delta-2*_alpha[j]*(delta - _epsilon[j]));
        }
        double Psi;
        double theta;
        double Delta;
        double dPsi_ddelta;
        double dDeltab_ddelta;
        for (int i = 40; i != 43; i++) {
            j = i - 40;
            Psi = exp(-_C[j]*(delta - 1)^^2 - _D[j]*(tau - 1)^^2);
                
            theta = (1 - tau) + _A[j]*((delta-1)^^2)^^(0.5/_beta2[j]);
            Delta = theta*theta + _B[j]*((delta - 1)^^2)^^_a[j];
            dPsi_ddelta = -2*_C[j]*(delta - 1)*Psi;
            dDeltab_ddelta = _b[j]*Delta^^(_b[j] - 1)*(delta - 1)*
                (2*theta*_A[j]/_beta2[j]*((delta - 1)^^2)^^(0.5/_beta2[j] - 1) + 2*_B[j]*_a[j]*((delta - 1)^^2)^^(_a[j] - 1));//fuck me;
            d_alphar_d_delta += _n[i]*(Delta^^_b[j]*(Psi + delta*dPsi_ddelta) + dDeltab_ddelta*delta*Psi);
            
        }
        return  d_alphar_d_delta;
    }

    const double getP_sat(double T)
    {
        //calculates saturated Pressure for a given temperature equation 3.13 on p16
        double summation = 0.0;
        for (int i = 1; i != 5; i++) {
            summation += _a_Psat[i]*(1.0 - T/_Tc)^^_t_Psat[i];
        }
        return exp(summation*_Tc/T)*_Pc;
    }

    const double getrho_satvap(double T)
    {
        double summation = 0;
        for (int i = 1; i != 6; i++) {
            summation += _a_rhosatvap[i]*(1.0 - T/_Tc)^^_t_rhosatvap[i];}
        return exp(summation)*_rhoc;
    }

    const double getT_satvap(double rho, double[2] bracket = [0.0, 304.1282], double tol = 1e-3) 
    {
        if (fabs(rho - _rhoc) < 0.0001) return _Tc;
        auto zeroRhosatvap_T = delegate (double T)
            {
            return rho - getrho_satvap(T);
            };
        return solve!(zeroRhosatvap_T,double)(bracket[0],bracket[1],tol);
    }

    const double getrho_satliq(double T)
    {
        double summation = 0;
        for (int i = 1; i != 5; i++) {
            summation += _a_rhosatliq[i]*(1.0 - T/_Tc)^^_t_rhosatliq[i];}
        return exp(summation)*_rhoc;
    }

    const double getT_satliq(double rho, double[2] bracket = [0.0, 304.1282], double tol = 1e-3)
    {
        if (fabs(rho - _rhoc) < 0.0001) return _Tc;
        auto zeroRhosatliq_T = delegate (double T)
            {
                return rho - getrho_satliq(T);
            };
        return solve!(zeroRhosatliq_T,double)(bracket[0],bracket[1],tol);
    }

} // end class CO2GasSW


version(co2gas_sw_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        auto gm = new CO2GasSW();
        assert(gm.species_name(0) == "CO2", failedUnitTest());
        auto gd = new GasState(gm, 7.38e6, 304.5);//point chosen close to critical point
        assert(approxEqual(gm.R(gd), 188.924, 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 1, failedUnitTest());
        assert(approxEqual(gd.p, 7.38e6, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.T, 304.5, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        assert(approxEqual(gd.rho, 340.128346, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.u, -159285.241490, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.a, 178.565447, 1.0e-4), failedUnitTest());

        gm.update_thermo_from_rhou(gd);
        gm.update_sound_speed(gd);
        assert(approxEqual(gd.p, 7.38e6, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.T, 304.5, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.a, 178.565447, 1.0e-4), failedUnitTest());

        gm.update_thermo_from_rhoT(gd);
        gm.update_sound_speed(gd);
        assert(approxEqual(gd.p, 7.38e6, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.u, -159285.241490, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.a, 178.565447, 1.0e-4), failedUnitTest());

        gm.update_thermo_from_rhop(gd);
        gm.update_sound_speed(gd);
        assert(approxEqual(gd.T, 304.5, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.u, -159285.241490, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.a, 178.565447, 1.0e-4), failedUnitTest());

        gm.update_thermo_from_ps(gd, -1183.950027);
        gm.update_sound_speed(gd);
        assert(approxEqual(gd.T, 304.5, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.rho, 340.128346, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.u, -159285.241490, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.a, 178.565447, 1.0e-4), failedUnitTest());
        gm.update_thermo_from_hs(gd, -137587.549656, -1183.950027);
        gm.update_sound_speed(gd);
        assert(approxEqual(gd.T, 304.5, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.p, 7.38e6, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.rho, 340.128346, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.u, -159285.241490, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.a, 178.565447, 1.0e-4), failedUnitTest());

        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/co2sw-gas-model.lua");
        gm = new CO2GasSW(L);
        lua_close(L);
        gd.p = 7.38e6;
        gd.T = 304.5;
        assert(approxEqual(gm.R(gd), 188.9241, 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 1, failedUnitTest());
        assert(approxEqual(gd.p, 7.38e6, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.T, 304.5, 1.0e-6), failedUnitTest());
        assert(approxEqual(gd.massf[0], 1.0, 1.0e-4), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        assert(approxEqual(gd.rho, 340.128346, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.u, -159285.241490, 1.0e-4), failedUnitTest());
        assert(approxEqual(gd.a, 178.565447, 1.0e-4), failedUnitTest());
        //looser tolerances on the updates below as they are using LuT
        gm.update_thermo_from_rhou(gd);
        gm.update_sound_speed(gd);

        assert(approxEqual(gd.p, 7.38e6, 1.0e-2), failedUnitTest());
        assert(approxEqual(gd.T, 304.5, 1.0e-2), failedUnitTest());
        assert(approxEqual(gd.a, 178.562447, 1.0e-2), failedUnitTest());
        return 0;
    }
}


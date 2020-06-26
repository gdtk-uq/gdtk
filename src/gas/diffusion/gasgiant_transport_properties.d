/**
 * Authors: Igor Segrovets and Rowan G.
 *
 */

module gas.diffusion.gasgiant_transport_properties;

import std.stdio;
import std.math;
import std.conv : to;

import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;

enum Sp {H2=0, H, Hp, He, Hep, e};

immutable double SMALL_ELECTRON_NUMBER_DENSITY = 1.0e-3;
immutable double R_u_cal = 1.987;//cal/mol/K

class GasGiantViscosity : Viscosity {
    this () {
        setGasGiantModelParameters();
        _g_22.length = _ns;
        _numden.length = _ns;
        _molef.length = _ns;
        _Delta_22.length = _ns;
 
        foreach (ref e; _g_22){
            e.length = _ns;
            foreach (ref ee; e){
                ee.length = 4;
            }
        }
        foreach (ref e; _Delta_22){
            e.length = _ns;
        }
        foreach (isp; 0 .. _ns) {       
            foreach (jsp; 0 .. isp+1){   
                if ( (isp in chargedSpecies) && (jsp in chargedSpecies) ) {continue;}   
                
                string key = speciesNames[isp] ~ ":" ~ speciesNames[jsp];
                if (!(key in G_22)){
                    key = speciesNames[jsp] ~ ":" ~ speciesNames[isp];
                }
                foreach (k, ref e; _g_22[isp][jsp]){
                    e = G_22[key][k];
                }
                _g_22[jsp][isp] = _g_22[isp][jsp];      
            }
        }
    }

    void attachGasModel(GasModel gm) {
        _gmodel = gm;
    }

    override GasGiantViscosity dup() const {
        return new GasGiantViscosity();
    }

    @nogc
    number eval(in GasState Q) {
        return viscosity(Q);    
    }
private:
    GasModel _gmodel;
    number[][][] _g_22;
    number[] _molef;
    number[] _numden; 
    number[][] _Delta_22;
    
    @nogc
    number viscosity(in GasState Q)
    {
        number T = Q.T;
        _gmodel.massf2molef(Q, _molef);
        number ne = Q.massf[Sp.e]*Q.rho / (_M[Sp.e]*1.0e-3);//numberden electron, _M mol/g->mol/kg
        ne *= 1e-6; //[1/m^3] --> [1/cm^3] 
        number sumA = 0.0;
        number sumB;
        foreach (isp; 0 .. _ns) {
            foreach (jsp; 0 .. isp+1) {
                number pi_Omega_22 = collisionIntegral(T, isp, jsp, ne);
                number _mu = _M[isp]*_M[jsp]/(_M[isp]+_M[jsp]);
                _Delta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*_mu/(to!double(PI)*R_u_cal*T))*pi_Omega_22;
                _Delta_22[jsp][isp] = _Delta_22[isp][jsp];
            }
        }                               
        foreach (isp; 0 .. _ns) {
            sumB = 0.0;
            if (_molef[isp] < SMALL_MOLE_FRACTION) continue;        
            foreach (jsp; 0 .. _ns) {
                if ((isp in chargedSpecies) && (jsp in chargedSpecies)){
                    if (ne == 0) continue; 
                }
                else {if (_molef[jsp] < SMALL_MOLE_FRACTION) continue;}
                sumB += _molef[jsp]*_Delta_22[isp][jsp];
            }
            number particleMass = (_M[isp]/Avogadro_number);
            sumA += particleMass*_molef[isp]/sumB;
        }
        number mu = sumA * (1.0e-3/1.0e-2); // convert g/(cm.s) -> kg/(m.s)
        return mu;
    }
    
    @nogc
    number collisionIntegral(number T, int isp, int jsp, number ne)
    {
        if ( (isp in chargedSpecies) && (jsp in chargedSpecies) ) { 
            return coulombicInteraction(T, isp, jsp, 2, ne);
        } else {
            return gupta(T, _g_22[isp][jsp]);
        }
    }
    
    @nogc
    number gupta(number T, number[] g)
    {
        auto x = log(T);
        auto A = g[0];
        auto B = g[1];
        auto C = g[2];
        auto D = g[3];
        auto tmpA = A*pow(x,2.0) + B*x + C;
        auto sig2Omega = D*pow(T,tmpA);
        sig2Omega = sig2Omega; 
        //PI/10 scale factor to match Bruno
        return sig2Omega;
    }
    
    @nogc
    number coulombicInteraction(number T, int isp, int jsp, int order, number ne)
    {
        auto lambda = 6.905*sqrt(T/ne); //debye shielding length cm as per source paper
        auto Tstar = 4132.5*(pow(T,(3.0/2.0))/sqrt(ne));
        double C, c, D;
        number sigma2Omega;
        //theory reproduced from Palmer and Wright 2003 DOI: 10.2514/2.6756
        if (((isp == Sp.Hep || isp == Sp.Hp)&&(jsp == Sp.Hep || jsp== Sp.Hp))||(isp==Sp.e && jsp==Sp.e)){
            // If similarly charged computed repulsive coulomb potential
            if (order == 2) {
                C = 0.157; 
                c = 0.0274;
                D = 1.235;
            } else if (order == 1) {
                C = 0.138; 
                c = 0.0106;
                D = 0.765;
            }
            auto tmpA = pow((lambda/Tstar),2.0);
            auto tmpB = D*Tstar*(1.0-C*exp(-c*Tstar));
            sigma2Omega = (5.0e+15)*tmpA*log(tmpB+1.0);
        } else {
            // else compute attractive coulomb potential
            if (order == 2){
                C = -0.146; 
                c = 0.0377;
                D = 1.262;
            } else if (order == 1){
                C = -0.476;     
                c = 0.0313;
                D = 0.784;
            }
            auto tmpA = pow((lambda/Tstar),2.0);
            auto tmpB = D*Tstar*(1.0-C*exp(-c*Tstar));
            sigma2Omega = (5.0e+15)*tmpA*log(tmpB+1.0);
        }
        return sigma2Omega;
    }
}

class GasGiantThermalConductivity : ThermalConductivity {

    override GasGiantThermalConductivity dup() const {
        return new GasGiantThermalConductivity();
    }

    //create data structures and populate
        //_g_11
    this () {
        setGasGiantModelParameters();
        _g_22.length = _ns;
        _g_11.length = _ns;
        _D_ij.length = _ns;
        _molef.length = _ns;
        _n.length = _ns;
        _Delta_11.length = _ns;
        _Delta_22.length = _ns;
        _a.length = _ns;
        _hl.length = _al.length;
        
        foreach (ref e; _g_11) {
            e.length = _ns;
            foreach (ref ee; e){
                ee.length = 4;
            }
        }
        foreach (ref e; _g_22){
            e.length = _ns;
            foreach (ref ee; e){
                ee.length = 4;
            }
        }

        foreach (ref e; _D_ij){e.length = _ns;}
        foreach (ref e; _Delta_11){e.length = _ns;}
        foreach (ref e; _Delta_22){e.length = _ns;}
        foreach (ref e; _a){e.length = _ns;}
        
        foreach (isp; 0 .. _ns) {
            foreach (jsp; 0 .. isp+1){//for lower diagonal j  
                foreach (k, ref e; _g_11[isp][jsp]){
                    if ( (isp in chargedSpecies) && (jsp in chargedSpecies) ) {
                        continue;
                    }
                    string key =speciesNames[isp] ~ ":" ~ speciesNames[jsp];
                    if (!(key in G_11)){
                        key = speciesNames[jsp] ~ ":" ~ speciesNames[isp];
                    }
                    e = G_11[key][k];
                }
                _g_11[jsp][isp] = _g_11[isp][jsp];
                //repeat for _22   
                foreach (k, ref e; _g_22[isp][jsp]){
                    if ( (isp in chargedSpecies) && (jsp in chargedSpecies) ) {
                        continue;
                    }
                    string key =speciesNames[isp] ~ ":" ~ speciesNames[jsp];
                    if (!(key in G_22)){
                        key = speciesNames[jsp] ~ ":" ~ speciesNames[isp];
                    }
                    e = G_22[key][k];
                }
                        _g_22[jsp][isp] = _g_22[isp][jsp];
            }
        }
    }//ends this()

    void attachGasModel(GasModel gm) {
        _gmodel = gm;
    }
    
    @nogc
    override number eval(ref const(GasState) Q, int imode) 
    {
        number K_fr = normalConductivity(Q, Q.T);
        number K_r = reactiveConductivity(Q, Q.T);
        number K = K_r + K_fr;
        return K;
    }
    
private:
    GasModel _gmodel;
    number[][][] _g_22;
    number[][][] _g_11;
    number[][] _D_ij; 
    number[][] _Delta_11;
    number[][] _Delta_22;
    number[][] _a;
    number[] _molef;
    number[] _n;
    number[] _hl;
    
    @nogc
    number normalConductivity(in GasState Q, number T)
    {
        //Equilibrium assumption        
        number K_tr = 0.0;
        number K_int = 0.0;
        number _tr;
        number _int;
        _gmodel.massf2molef(Q, _molef);
        number ne = eval_ne(Q);
        foreach (isp; 0 .. _ns){
            foreach (jsp; 0 .. isp+1){
                // compute Delta_22
                number pi_Omega_22 = collisionIntegral(T, isp, jsp, 2, ne);
                number _mu = _M[isp]*_M[jsp]/(_M[isp]+_M[jsp]);
                _Delta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*_mu/(to!double(PI)*R_u_cal*T))*pi_Omega_22;
                _Delta_22[jsp][isp] = _Delta_22[isp][jsp];
                // compute Delta_11
                number pi_Omega_11 = collisionIntegral(T, isp, jsp, 1, ne);
                _Delta_11[isp][jsp] = (8./3)*1.546e-20*sqrt(2.0*_mu/(to!double(PI)*R_u_cal*T))*pi_Omega_11;
                _Delta_11[jsp][isp] = _Delta_11[isp][jsp];      
                //compute a_ij
                number r_M = _M[isp]/_M[jsp]; 
                number tmpA = pow((1+r_M),2.0);
                number tmpB = (0.45-2.54*r_M); 
                _a[isp][jsp] = 1 + (1-r_M)*tmpB/tmpA;
                _a[jsp][isp] = _a[isp][jsp];
            }
        }
        
        foreach (isp; 0 .. _ns){
            _tr = 0;
            _int = 0;
            number cp_i = _gmodel.Cp(Q, isp); // for a species, J/(kg.K)
            number mol_mass_i = _M[isp]*1.0e-3; //g/mol ---> kg/mol
            number cp_i_int = (cp_i/(R_universal/mol_mass_i) - (5./2.));//[J/kg-K]/[J/mol-K][mol/kg]--> unitless
            if (cp_i_int<0) {cp_i_int=0;}   
            foreach (jsp; 0 .. _ns){        
                if ((isp in chargedSpecies) && (jsp in chargedSpecies)){
                    if (ne == 0) continue;
                }
                //K_tr
                _tr += _a[isp][jsp]*_Delta_22[isp][jsp]*_molef[jsp];
                //K_int
                _int += _molef[jsp]*_Delta_11[isp][jsp];        
            }
            K_tr += (15.0/4.0)*Boltzmann_constant*_molef[isp]/_tr; //in J/cm-K-sec 
            K_int += Boltzmann_constant*cp_i_int*_molef[isp]/_int;
        }       
        number K_f = (K_tr + K_int)*100.0; // J/cm-K-sec *100cm/m --->J/m-K-sec 
        return K_f;
    }

    @nogc
    number reactiveConductivity(in GasState Q, number T)
    {
        number ne = eval_ne(Q);//in no./cm^3
        _gmodel.massf2molef(Q, _molef);
        number h_i, h_l, i_sum, j_sum;
        number K_r = 0;
        int L;
        if ((T>=1500)&&(T<7000)){
            L = 0;
        } else if ((T>=7000)&&(T<20000)){
            L = 1;
        } else {        
            K_r = 0;
            return K_r;
        }
        h_l = 0;
        foreach (isp; 0 .. _ns){
            //compute h_l for reaction l
            h_i = _gmodel.enthalpy(Q, isp); // J/kg
            h_l += h_i*_al[L][isp]*(_M[isp]*1e-3);//convert mol mass g/mol to kg/mol
        }// J/mol units out
        _hl[L] = h_l;
        
        foreach (isp; 0 .. _ns){
            foreach (jsp; 0 .. _ns){
                // compute Delta_11
                number _mu = _M[isp]*_M[jsp]/(_M[isp]+_M[jsp]);
                number pi_Omega_11 = collisionIntegral(T, isp, jsp, 1, ne);
                _Delta_11[isp][jsp] = (8./3)*1.546e-20*sqrt(2.0*_mu/(to!double(PI)*R_u_cal*T))*pi_Omega_11;
                _Delta_11[jsp][isp] = _Delta_11[isp][jsp];                              
            }
        }
        number l_sum = 0;
        i_sum = 0;
        foreach (isp; 0 .. _ns){ 
            j_sum = 0;
            if (_molef[isp] < SMALL_MOLE_FRACTION) continue;        
            foreach (jsp; 0 .. _ns){
                if ((isp in chargedSpecies) && (jsp in chargedSpecies)){
                    if (ne == 0) continue;} else {if (_molef[jsp] < SMALL_MOLE_FRACTION) continue;} 
                
                j_sum += (_al[L][isp]*_molef[jsp] - _al[L][jsp]*_molef[isp])*_Delta_11[isp][jsp];
            }
            i_sum += j_sum*_al[L][isp]/_molef[isp]; 
        }
        if (i_sum != 0) {
            l_sum += ((_hl[L]/(R_universal*T))^^2)/(i_sum);                 
            K_r = l_sum*Boltzmann_constant;
            K_r *= 1e2; // per cm --> per m
        }
        return K_r;
    }

    @nogc
    number binaryDiffusion(in GasState Q, int k, int l)
    {
        number T = Q.T;
        number ne = eval_ne(Q); 
        foreach (isp; 0 .. _ns){
            if (_molef[isp] < SMALL_MOLE_FRACTION) continue;
            foreach (jsp; 0 .. _ns){
                if (_molef[jsp] < SMALL_MOLE_FRACTION) continue;
                if ((isp in chargedSpecies) && (jsp in chargedSpecies)){
                    if (ne < SMALL_ELECTRON_NUMBER_DENSITY) continue;
                }       
                number g = (2.628e-7);
                g *= sqrt(pow(T,3.0)*(_M[isp]+_M[jsp])/(2*_M[isp]*_M[jsp]));
                number sig2Omega_11 = collisionIntegral(T, isp, jsp, 1, ne);
                auto P_in_atm = Q.p/P_atm;
                g /= P_in_atm*sig2Omega_11;
                _D_ij[isp][jsp] = g;
            }
        }
        return _D_ij[k][l];
    }

    @nogc
    number eval_ne(in GasState Q)
    {
        number ne = Avogadro_number*Q.massf[Sp.e]*Q.rho / (_M[Sp.e]*1.0e-3); //number density of electron
        ne *= 1e-6; //1/m^3 --> 1/cm^3
        return ne;
    }
    
    @nogc
    number eval_n(in GasState Q, int isp)
    {
        number n = Avogadro_number*Q.massf[isp]*Q.rho / (_M[isp]*1.0e-3); //number density of electron
        return n; //returns in n/m^3
    }
    
    @nogc
    number collisionIntegral(number T, int isp, int jsp, int order, number ne)
    {
        if ( (isp in chargedSpecies) && (jsp in chargedSpecies) ) {
            return coulombicInteraction(T, isp, jsp, order, ne);
        } else if (order == 1){
            return gupta(T, _g_11[isp][jsp]);
        } else{
            return gupta(T, _g_22[isp][jsp]);
        }
    }
    
    @nogc
    number gupta(number T, number[] g)
    {
        auto x = log(T);
        auto A = g[0];
        auto B = g[1];
        auto C = g[2];
        auto D = g[3];
        auto tmpA = A*pow(x,2.0) + B*x + C;
        auto sig2Omega = D*pow(T,tmpA);
        sig2Omega = sig2Omega; 
        //PI/10 scale factor to match Bruno
        return sig2Omega;
    }
    
    @nogc
    number coulombicInteraction(number T, int isp, int jsp, int order, number ne)
    {
        auto lambda = 6.905*sqrt(T/ne); //debye shielding length cm
        auto Tstar = 4132.5*(pow(T,(3.0/2.0))/sqrt(ne));
        double C, c, D;
        double A, B, b;
        number sigma2Omega;
        //theory reproduced from Palmer and Wright 2003 DOI: 10.2514/2.6756
        if (((isp == Sp.Hep || isp == Sp.Hp)&&(jsp == Sp.Hep || jsp== Sp.Hp))||(isp==Sp.e && jsp==Sp.e)){
            //check if same charge and potential is repulsive
            if (order == 2){
                C = 0.157; 
                c = 0.0274;
                D = 1.235;
            } else if (order == 1){
                C = 0.138; 
                c = 0.0106;
                D = 0.765;
            } 
            auto tmpA = pow((lambda/Tstar),2.0);
            auto tmpB = D*Tstar*(1.0-C*exp(-c*Tstar));
            sigma2Omega = (5.0e+15)*tmpA*log(tmpB+1.0);                     
        } else {
            //else attractive potential
            if (order == 2){
                C = -0.146; 
                c = 0.0377;
                D = 1.262;
            } else if (order == 1){
                C = -0.476;     
                c = 0.0313;
                D = 0.784;
            }
            auto tmpA = pow((lambda/Tstar),2.0);
            auto tmpB = D*Tstar*(1.0-C*exp(-c*Tstar));
            sigma2Omega = (5.0e+15)*tmpA*log(tmpB+1.0);
        }
        return sigma2Omega;
    }
}


version(gasgiant_transport_properties_test) {
    import std.stdio;
    int main()
    {
        GasGiantViscosity ggv = new GasGiantViscosity();
        auto Q = new GasState(1, 1);
        auto f = File("gasgiant-visc-test.dat", "w");
        f.writeln("# T(K)   mu-He  mu-H2  mu-H");
        double T = 300.0;
        double dT = 100.0;
        double Tmax = 10000.0;
        /*
        while (T <= Tmax) {
            Q.T = T;
            auto mu0 = ggv.selfViscosity(Q.T, 0, 0, 0, 0.0);
            auto mu1 = ggv.selfViscosity(Q.T, 1, 1, 0, 0.0);
            auto mu2 = ggv.selfViscosity(Q.T, 2, 2, 0, 0.0);
            f.writefln("%f %.8e %.8e %.8e", T, mu0, mu1, mu2);
            T += dT;
        }
        */
        f.close();
        return 0;
    }
}




static double[6] _M;
static double[6] _m;
static double[6][3] _al;
static string[6] speciesNames = ["H2", "H", "Hp", "He", "Hep", "e"];
static double[4][string] G_22;
static double[4][string] G_11;
static int _ns = _M.length; // number species
static bool[int] chargedSpecies;
static paramsInitialised = false;

void setGasGiantModelParameters()
{
    if (paramsInitialised) return;
    else {
        paramsInitialised = true;
    }

    chargedSpecies[Sp.Hep] = true;
    chargedSpecies[Sp.Hp] = true;
    chargedSpecies[Sp.e] = true;

        
//_al holds stiochometric coefficients for reactions occuring in the gas
//each reaction can move in the forward and backwards direction
//we consider H dissociation  _al[0]
        _al[0] = [1, -2, 0, 0, 0, 0]; //H2 H H+ He Hep e
// H ionisation _al[1] 
        _al[1] = [0, 1, -1, 0, 0, -1];  //written in indep form
// He ionisation _al[2] 
        _al[2] = [0, 0, 0, -1, 1, 1]; //here both are independent form, shouldnt matter

    _M[Sp.He] = 4.002602; 
    _m[Sp.He] = (_M[Sp.He]/Avogadro_number)*10e-3;
    
    _M[Sp.H2] = 2.015882;
    _m[Sp.H2] = (_M[Sp.H2]/Avogadro_number)*10e-3;

    _M[Sp.H] = 2.015882/2.;
    _m[Sp.H] = (_M[Sp.H2]/Avogadro_number)*10e-3;

    _M[Sp.Hep] = 4.002602 - 9.1093e-29*Avogadro_number;
    _m[Sp.Hep] = (_M[Sp.He]/Avogadro_number)*10e-3;
    
    _M[Sp.Hp] = 2.015882/2. - 9.1093e-29*Avogadro_number;//g/mol
    _m[Sp.Hp] = (_M[Sp.H2]/Avogadro_number)*10e-3;

    _M[Sp.e] = 9.10938188e-31*Avogadro_number*1e3;//g/mol
    _m[Sp.e] = 9.10938188e-31;//kg

    G_22["He:He"]  = [-0.00525,0.09660,-0.803,165];
    G_22["H2:H2"]  = [-0.00546,0.09600,-0.767,221];
    G_22["H:H"]    = [-0.01990,0.41300,-3.100,75900];
    G_22["H2:H"]   = [-0.00331,0.00868,0.0265,24];
    G_22["He:H"]   = [-0.00767,0.12800,-0.953,266];
    G_22["He:H2"]  = [-0.00346,0.05240,-0.462,94.1];
    G_22["H:Hp"]   = [-0.02390,0.47000,-3.310,202000];
    G_22["H2:Hp"]  = [-0.00837,0.14900,-1.170,1430];
    G_22["He:Hp"]  = [-0.02730,0.56000,-4.180,1260000];
    G_22["H:e"]    = [-0.00190,-0.0328,0.5050,9.1];
    G_22["H2:e"]   = [-0.01950,0.38900,-2.300,372];
    G_22["He:e"]   = [-0.00649,0.12000,-0.691,19.5];
    G_22["H:Hep"]  = [-0.01620,0.28300,-1.860,1950];
    G_22["H2:Hep"] = [-0.00115,0.03910,-0.874,2810];
    G_22["He:Hep"] = [-0.00877,0.15000,-1.090,753];

    G_11["He:He"]  = [-0.00358,0.0568,-0.514,72.9];
    G_11["H2:H2"]  = [-0.00562,0.0977,-0.795,218];
    G_11["H:H"]    = [-0.0188,0.359,-2.51,11500];
    G_11["H2:H"]   = [-0.00147,-0.0279,0.227,15];
    G_11["He:H"]   = [-0.00682,0.106,-0.808,166];
    G_11["He:H2"]  = [-0.00349,0.051,-0.46,83.9];
    G_11["H:Hp"]   = [-0.000997,0.0226,-0.319,504];
    G_11["H2:Hp"]  = [-0.018,0.322,-2.17,10700];
    G_11["He:Hp"]  = [-0.0205,0.391,-2.91,76000];
    G_11["H:e"]    = [-0.00618,0.0682,-0.201,44.7];
    G_11["H2:e"]   = [-0.0227,0.464,-2.9,2710];
    G_11["He:e"]   = [-0.00824,0.168,-1.06,46.5];
    G_11["H:Hep"]  = [-0.0178,0.282,-1.64,930];
    G_11["H2:Hep"] = [0.000291,0.0154,-0.792,2860];
    G_11["He:Hep"] = [0.000824,-0.0219,0.0389,115];
}



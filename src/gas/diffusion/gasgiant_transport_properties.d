/**
 * Authors: Igor Segrovets and Rowan G.
 * Full script
 */

module gas.diffusion.gasgiant_transport_properties;

import std.stdio;
import std.string;
import std.math;
import std.conv : to;

import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;

enum Sp {H2=0, H, Hp, He, Hep, e};//species
enum pR {H2He=0, H2Hep};		  //phenomenological reaction
enum Db {Bruno=0, Palmer};		  //database
enum Ct {HeHep=0, HHp};			  //charge transfer odd l CCS

immutable double SMALL_ELECTRON_NUMBER_DENSITY = 1.0e-3;
immutable double R_u_cal = 1.987;//cal/mol/K


class GasGiantCollisionCrossSections  {
    this()
    {
        setGasGiantModelParameters();
        _g_22.length = _ns;
        _g_11.length = _ns;
        _a_22.length = _ns;
        _a_11.length = _ns;
        _c_11.length = 2;
        _c_22.length = 2;
        _p_11.length = 2;
        _p_22.length = 2; 
			
        foreach (ref e; _g_11){
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
        foreach (ref e; _a_11){
            e.length = _ns;
            foreach (ref ee; e){
                ee.length = 8;
            }
        }
        foreach (ref e; _a_22){
            e.length = _ns;
            foreach (ref ee; e){
                ee.length = 8;
            }
        }
        foreach (ref e; _c_11){
            e.length = 7;
            foreach (ref ee; e){
                ee.length = 3;
            }
        }
        foreach (ref e; _c_22){
            e.length = 7;
            foreach (ref ee; e){
                ee.length = 3;
            }
        }
        foreach (ref e; _p_11){e.length = 7;}
        foreach (ref e; _p_22){e.length = 7;}
		    
        //gupta
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
        //bruno
        foreach (isp; 0 .. _ns) {
            foreach (jsp; 0 .. isp+1){  //for lower diagonal j  
                foreach (k, ref e; _a_11[isp][jsp]){
                    if ( (isp in chargedSpecies) && (jsp in chargedSpecies) ) {
                        continue;
                    }
                    string key =speciesNames[isp] ~ ":" ~ speciesNames[jsp];
                    if (!(key in a_11)){
                        key = speciesNames[jsp] ~ ":" ~ speciesNames[isp];
                    }
                    e = a_11[key][k];
                }
                _a_11[jsp][isp] = _a_11[isp][jsp];
                //repeat for _22   
                foreach (k, ref e; _a_22[isp][jsp]){
                    if ( (isp in chargedSpecies) && (jsp in chargedSpecies) ) {
                        continue;
                    }
                    string key =speciesNames[isp] ~ ":" ~ speciesNames[jsp];
                    if (!(key in a_22)){
                        key = speciesNames[jsp] ~ ":" ~ speciesNames[isp];
                    }
                    e = a_22[key][k];
                }
                _a_22[jsp][isp] = _a_22[isp][jsp];
            }
        }
        //_c_ coeffs for phenomological potentials
        foreach (mtype; 0 .. 2){
            foreach (abeta; 0 .. 7){
                foreach (k, ref e; _c_11[mtype][abeta]){
                    int mNum = 6 - 2*mtype;
                    int _abeta = abeta + 1;
                    string key= "a" ~ to!string(_abeta) ~":m" ~ to!string(mNum);
                    e = c_11[key][k];
                }
                // repeat for _22	
                foreach (k, ref e; _c_22[mtype][abeta]){
                    int mNum = 6 - 2*mtype;
                    int _abeta = abeta + 1;
                    string key= "a" ~ to!string(_abeta) ~":m" ~ to!string(mNum);
                    e = c_22[key][k];
                }
            }
            foreach(k, ref e; _p_11[mtype]){
                number tempSum = 0;
                foreach(cj; 0 .. 3){
                    tempSum += _c_11[mtype][k][cj]*_beta[mtype]^^cj;
                }
                e = tempSum;
            }
            foreach(k, ref e; _p_22[mtype]){
                number tempSum = 0;
                foreach(cj; 0 .. 3){
                    tempSum += _c_22[mtype][k][cj]*_beta[mtype]^^cj;
                }
                e = tempSum;
            }
        }
    }//ends this()
  
    void attachGasModel(GasModel gm) {
        _gmodel = gm;
    }

    @nogc
    number eval(in GasState Q, int isp, int jsp, int order, int DB)
    {
        number ne = eval_ne(Q);
        number css = collisionIntegral(Q.T, isp, jsp, order, ne, DB);
        return css;
    }

private:
    GasModel _gmodel;
    number[][][] _g_22;
    number[][][] _g_11;
    number[][][] _a_22;
    number[][][] _a_11;
    number[][][] _c_22;
    number[][][] _c_11;
    number[][] _p_22;
    number[][] _p_11;
	

    @nogc
    number eval_ne(in GasState Q)
    {
        number ne = Avogadro_number*Q.massf[Sp.e]*Q.rho / (_M[Sp.e]*1.0e-3); //number density of electron
        ne *= 1e-6; //1/m^3 --> 1/cm^3
        return ne;
    }

    @nogc
    number collisionIntegral(number T, int isp, int jsp, int order, number ne, int database)
    {
        if ( (isp in chargedSpecies) && (jsp in chargedSpecies) ) {
            return to!number(PI)*coulombicInteraction(T, isp, jsp, order, ne, database);
        } else if (order == 1){ //order is 11
            if (database == Db.Palmer){ // database is Palmer
                return gupta(T, _g_11[isp][jsp]);
            } else { // database is Bruno
                return to!number(PI)*bruno(T, isp, jsp, 1);
            }
        }else{	//order is 22
            if (database == Db.Palmer){
                return gupta(T, _g_22[isp][jsp]);
            } else {
                return to!number(PI)*bruno(T, isp, jsp, 2);
            }
        }
    }
    
    @nogc
    number bruno(number T, int isp, int jsp, int order){
        // bruno equation handler
        if (order == 1){
            if (( (isp==Sp.e)&&(jsp !in chargedSpecies) )||( (jsp==Sp.e)&&(isp !in chargedSpecies) )){
                //If e- with neutral -> EQ19
                return electronHeavyInteraction(T, _a_11[isp][jsp]);
            } else if(( (isp==Sp.H2)&&(jsp in phenomSpecies) )||( (jsp==Sp.H2)&&(isp in phenomSpecies) )){
                //If phenomological potential -> EQ9
                if ( (isp==Sp.He)||(jsp==Sp.He) ){
                    //If m=6 uncharged reaction
                    return phenomInteraction(T, _p_11[pR.H2He], _phi[pR.H2He], _sigma[pR.H2He]);
                } else {
                    //If m=4 
                    return phenomInteraction(T, _p_11[pR.H2Hep], _phi[pR.H2Hep], _sigma[pR.H2Hep]);
                }
            } else{
                //Else must be remaining interactions
                //If He-Hep or H-Hp add resonant charge trasnfer term
                if ( ((isp==Sp.He)&&(jsp==Sp.Hep))||((jsp==Sp.He)&&(isp==Sp.Hep)) ){
                    auto Omega_el = heavyInteraction(T,_a_11[isp][jsp]);
                    auto Omega_ch = chargeTransfer(T, _d_11[Ct.HeHep]);
                    return (Omega_el^^2 + Omega_ch^^2)^^(0.5);
                } else if  ( ((isp==Sp.H)&&(jsp==Sp.Hp))||((jsp==Sp.H)&&(isp==Sp.Hp)) ){
                    auto Omega_el = heavyInteraction(T,_a_11[isp][jsp]);
                    auto Omega_ch = chargeTransfer(T, _d_11[Ct.HHp]);
                    return (Omega_el^^2 + Omega_ch^^2)^^(0.5);				
                } else {
                    return heavyInteraction(T, _a_11[isp][jsp]);
                }
            }
        } else { //repeat for order 22
            if (( (isp==Sp.e)&&(jsp !in chargedSpecies) )||( (jsp==Sp.e)&&(isp !in chargedSpecies) )){
                return electronHeavyInteraction(T, _a_22[isp][jsp]);
            } else if(( (isp==Sp.H2)&&(jsp in phenomSpecies) )||( (jsp==Sp.H2)&&(isp in phenomSpecies) )){
                if ( (isp==Sp.He)||(jsp==Sp.He) ){
                    return phenomInteraction(T, _p_22[pR.H2He], _phi[pR.H2He],_sigma[pR.H2He]);
                } else {
                    return phenomInteraction(T, _p_22[pR.H2Hep], _phi[pR.H2Hep],_sigma[pR.H2Hep]);
                }
            } else{
                return heavyInteraction(T, _a_22[isp][jsp]);
            }
        }
    }
    @nogc
    number chargeTransfer(number T, double[] d)
    {
        auto x  = log(T);
        auto d1 = d[0];
        auto d2 = d[1];
        auto d3 = d[2];
        //BRUNO EQ17
        auto sig2Omega_ch = d1 + d2*x + d3*x^^2; 
        return sig2Omega_ch;
    }
    @nogc
    number heavyInteraction(number T, number[] a)
    {
        auto x = log(T);
        auto a1 = a[0];
        auto a2 = a[1];
        auto a3 = a[2];
        auto a4 = a[3];
        auto a5 = a[4];
        auto a6 = a[5];
        auto a7 = a[6];
        auto a8 = a[7];
        //BRUNO EQ11
        auto tmpA = exp((x-a3)/a4);
        auto tmpB = exp((a3-x)/a4);
        auto tmpC = exp((x-a6)/a7);
        auto tmpD = exp((a6-x)/a7);

        auto sig2Omega = (a1+a2*x)*tmpA / (tmpA + tmpB);
        sig2Omega += a5*tmpC / (tmpC + tmpD);
        return sig2Omega;
    }
    @nogc
    number electronHeavyInteraction(number T, number[] a)//pass in _a_11 or _a_22[isp][jsp]
    {
        auto x = log(T);
        auto a1 = a[0];
        auto a2 = a[1];
        auto a3 = a[2];
        auto a4 = a[3];
        auto a5 = a[4];
        auto a6 = a[5];
        auto a7 = a[6];
        auto a8 = a[7];
        //BRUNO EQ19
        auto tmpE = pow(x,a5);
        auto tmpF = exp((x-a1)/a2);
        auto tmpG = exp(-(x-a1)/a2);
        auto tmpH = exp(-pow(((x-a7)/a8),2));

        auto sig2Omega = a3*tmpE*tmpF / (tmpF+tmpG);
        sig2Omega += a6*tmpH + a4;
        return sig2Omega;
    }
	@nogc
	number phenomInteraction(number T, number[] p, double phi, double sigma)
	{ 
            //pass in _p_11/_22  _p_xx[pR.H2He]/_p_xx[pR.H2Hep]
            auto x = log(T*kb_eV/phi);
            auto a1 = p[0];
            auto a2 = p[1];
            auto a3 = p[2];
            auto a4 = p[3];
            auto a5 = p[4];
            auto a6 = p[5];
            auto a7 = p[6];
            auto _xp1 = exp((x-a3)/a4);
            auto _xp2 = exp((a3-x)/a4);
            auto _xp3 = exp((x-a6)/a7);
            auto _xp4 = exp((a6-x)/a7);
            auto Omega = exp( ((a1+a2*x)*_xp1/(_xp1 + _xp2)) + (a5*_xp3/(_xp3 + _xp4)) );
            return sigma*Omega;
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
        return sig2Omega;
    }
    
    @nogc
    number coulombicInteraction(number T, int isp, int jsp, int order, number ne, int DB)
    {
		int ion_density_scale;
		if ( DB == Db.Bruno ){
			ion_density_scale = 2;
		} else {
			ion_density_scale = 1;
		}

		auto lambda = 6.905*sqrt(T/(ne*ion_density_scale)); //debye shielding length cm
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

class GasGiantViscosity : ViscosityMixtureModel {
//	GasGiantCollisionCrossSections _ggccs = new GasGiantCollisionCrossSections();
    this(string dbOption) {
        setGasGiantModelParameters();
        _ggccs = new GasGiantCollisionCrossSections(); // set an identifier for ggccs except it doesnt work

        _numden.length = _ns;
        _molef.length = _ns;
        _Delta_22.length = _ns;
	 
        foreach (ref e; _Delta_22){
            e.length = _ns;
        }
        _dbOption = dbOption;
    }

    void attachGasModel(GasModel gm) {
	_gmodel = gm;
    }

    override GasGiantViscosity dup() const {
        return new GasGiantViscosity(_dbOption);
    }

    @nogc
    override number eval(in GasState Q) {
        switch (_dbOption) {
        case "Bruno":
            return viscosity(Q, Db.Bruno);
        case "Palmer":
            return viscosity(Q, Db.Palmer);
        default:
            string errMsg = "GasGiantViscosity: The requested dbOption is not available.";
            throw new Error(errMsg);
        }
    }

    @nogc
    number evalwPalmer(in GasState Q){
        return viscosity(Q, Db.Palmer);
    }

private:
    string _dbOption;
    GasModel _gmodel;
    number[][][] _g_22;
    number[] _molef;
    number[] _numden; 
    number[][] _Delta_22;
    GasGiantCollisionCrossSections _ggccs;

    @nogc
    number viscosity(in GasState Q, int DB)
    {
        number T = Q.T;
        _gmodel.massf2molef(Q, _molef);
        number ne = Q.massf[Sp.e]*Q.rho / (_M[Sp.e]*1.0e-3);//numberden electron, _M mol/g->mol/kg
        ne *= 1e-6; //[1/m^3] --> [1/cm^3] 
        number sumA = 0.0;
        number sumB;
        foreach (isp; 0 .. _ns) {
            foreach (jsp; 0 .. isp+1) {
                number pi_Omega_22 = _ggccs.eval(Q, isp, jsp, 2, DB); // here call a method from the ggccs class using its identifier, but identifier not recognised.
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

}


class GasGiantThermalConductivity : ThermalConductivityMixtureModel {

    override GasGiantThermalConductivity dup() const {
        return new GasGiantThermalConductivity(_dbOption, _modelOption);
    }

    //create data structures and populate
    //_g_11
    this(string dbOption, string modelOption) {
        _dbOption = dbOption;
        _modelOption = modelOption;
        setGasGiantModelParameters();
        _ggccs = new GasGiantCollisionCrossSections();
        _molef.length = _ns;
        _n.length = _ns;
        _Delta_11.length = _ns;
        _Delta_22.length = _ns;
        _a.length = _ns;
        _hl.length = _al.length;
			
        foreach (ref e; _Delta_11){e.length = _ns;}
        foreach (ref e; _Delta_22){e.length = _ns;}
        foreach (ref e; _a){e.length = _ns;}
		    
    }//ends this()
    void attachGasModel(GasModel gm) {
        _gmodel = gm;
    }

    @nogc
    override number eval(ref const(GasState) Q, int imode)
    {
        switch (_modelOption) {
        case "frozen" :
            switch (_dbOption) {
            case "Bruno":
                return evalGuptaFrozen(Q, 0);
            case "Palmer":
                return evalGuptaFrozen_P(Q, 0);
            case "Eucken":
                return evalEuckenFrozen(Q, 0);
            default:
                string errMsg = "GasGiantThermalConductivity: The requested dbOption is not available.";
                throw new Error(errMsg);
            }
        case "reactive":
            switch (_dbOption) {
            case "Bruno":
                number K_fr = normalConductivity(Q, Q.T, Db.Bruno);
                number K_r = reactiveConductivity(Q, Q.T, Db.Bruno);
                number K = K_r + K_fr;
                return K;
            case "Palmer":
                return evalwPalmer(Q,0);
            case "Eucken":
                return evalEuckenwReactive(Q, 0);
            case "Eucken-Palmer":
                return evalEuckenwReactive_P(Q, 0);
            default:
                string errMsg = "GasGiantThermalConductivity: The requested dbOption is not available.";
                throw new Error(errMsg);
            }
        default:
            string errMsg = "GasGiantThermalConductivity: The requested modelOption '%s' is not available.";
            throw new Error(errMsg);
        }

    }

    @nogc
    number evalwPalmer(ref const(GasState) Q, int imode)
    {
        number K_fr = normalConductivity(Q, Q.T, Db.Palmer);
        number K_r = reactiveConductivity(Q, Q.T, Db.Palmer);
        //debug{writeln("palmer T=",T," K_fr = ",K_fr," K_r = ",K_r);}
        number K = K_r + K_fr;
        return K;
    }

    @nogc
    number evalEuckenwReactive(ref const(GasState) Q, int imode)
    {
        number K_fr = mixedEuckenConductivity(Q);
        number K_r = reactiveConductivity(Q, Q.T, Db.Bruno);
        number K = K_r + K_fr;
        return K;
    }

    @nogc
    number evalEuckenwReactive_P(ref const(GasState) Q, int imode)
    {
        number K_fr = mixedEuckenConductivity(Q);
        number K_r = reactiveConductivity(Q, Q.T, Db.Palmer);
        number K = K_r + K_fr;
        return K;
    }
	
    @nogc
    number evalEuckenFrozen(ref const(GasState) Q, int imode)
    {
        number K_fr = mixedEuckenConductivity(Q);
        return K_fr;
    }

    @nogc
    number evalGuptaFrozen(ref const(GasState) Q, int imode)
    {
        number K_fr = normalConductivity(Q, Q.T, Db.Bruno);
        return K_fr;
    }
	
    @nogc
    number evalGuptaFrozen_P(ref const(GasState) Q, int imode)
    {
        number K_fr = normalConductivity(Q, Q.T, Db.Palmer);
        return K_fr;
    }
private:
    string _dbOption, _modelOption;
    GasModel _gmodel; 
    number[][] _Delta_11;
    number[][] _Delta_22;
    number[][] _a; //used in normalConductivity
    number[] _molef;
    number[] _n;
    number[] _hl;
    GasGiantCollisionCrossSections _ggccs;
	
    @nogc
    number normalConductivity(in GasState Q, number T, int DB)
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
                number pi_Omega_22 = _ggccs.eval(Q, isp, jsp, 2, DB);
                number _mu = _M[isp]*_M[jsp]/(_M[isp]+_M[jsp]);
                _Delta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*_mu/(to!double(PI)*R_u_cal*T))*pi_Omega_22;
                _Delta_22[jsp][isp] = _Delta_22[isp][jsp];
                // compute Delta_11
                number pi_Omega_11 = _ggccs.eval(Q, isp, jsp, 1, DB);
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
    number reactiveConductivity(in GasState Q, number T, int DB)
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
                number pi_Omega_11 = _ggccs.eval(Q, isp, jsp, 1, DB);
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
    number mixedEuckenConductivity(in GasState Q)
    {   
        number Sc =1.32;
        number cp_i;
        number cp_mix = 0;
        number M_mix = 0;
        _gmodel.massf2molef(Q, _molef);
        
        foreach (isp; 0 .. _ns){    
            cp_i = _gmodel.Cp(Q, isp);
            if (isp == Sp.e){cp_i=0;}
            cp_mix += cp_i*_molef[isp];
            number tmp = _M[isp]*10e-3;
            M_mix += tmp*_molef[isp];
        }
        cp_mix = (cp_mix / (R_universal/M_mix))-5./2.;
        number mu = Q.mu;
        number K = ( 15./4. + Sc*cp_mix ) *R_universal*mu /M_mix; //bracket term is unitless
        return K;
    }

    @nogc
    number eval_ne(in GasState Q)
    {
        number ne = Avogadro_number*Q.massf[Sp.e]*Q.rho / (_M[Sp.e]*1.0e-3); //number density of electron
        ne *= 1e-6; //1/m^3 --> 1/cm^3
        return ne;
    }
}


version(gasgiant_transport_properties_test) {
    import std.stdio;
    int main()
    {
        GasGiantViscosity ggv = new GasGiantViscosity();
        auto Q = GasState(1, 1);
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
static string[6] speciesNames = ["H2", "H", "Hp", "He","Hep","e"];
static double[4][string] G_22;
static double[4][string] G_11;
static int _ns = _M.length; // number species
static bool[int] chargedSpecies;
static bool[int] phenomSpecies;
static double[8][string] a_22;
static double[8][string] a_11;
static double[3][string] c_11;
static double[3][string] c_22;
static double[2] _beta;
static double[2] _phi;
static double[2] _sigma;
static double kb_eV;
static double[3][2] _d_11; 

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

    phenomSpecies[Sp.He] = true;
    phenomSpecies[Sp.Hep] = true;
    // from Bruno-2010, doi:10.1063/1.3495980
    _beta[pR.H2He] =  9.29; // He-H2
    _beta[pR.H2Hep] =  9.22; // Hep-H2
    _phi[pR.H2He]  = 1.865e-3; // eV 
    _phi[pR.H2Hep]  = 177.5e-3; // eV
    // from Laricchiuta-2007, doi:10.1016/j.cplett.2007.07.097 
    _sigma[pR.H2He] = (0.8002*(9.29^^(0.049256))*3.189)^^2; // Angstrom sq
    _sigma[pR.H2Hep] = (0.7564*(9.22^^(0.064605))*2.223)^^2; // Angstrom sq

    kb_eV = 1380649.0/16021766340.0;//eV/K

    //_al holds stiochometric coefficients for reactions occuring in the gas
    //each reaction can move in the forward and backwards direction
    // H dissociation  _al[0]
    _al[0] = [1, -2, 0, 0, 0, 0]; //H2 H H+ He Hep e
    // H ionisation _al[1] 
    _al[1] = [0, 1, -1, 0, 0, -1];	//written in indep form
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
    //G_ parameters from Palmer-2013, doi:10.2514/1.A32768
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
    // c_ parameters from Laricchiuta-2007, doi:10.1016/j.cplett.2007.07.097 
    c_11["a1:m6"] = [0.7884756, -0.02438494, 0.0];
    c_11["a2:m6"] = [-0.2952759, -0.001744149, 0.0];
    c_11["a3:m6"] = [0.5020892, 0.04316985, 0.0];
    c_11["a4:m6"] = [-0.9042460, -0.04017103, 0.0];
    c_11["a5:m6"] = [-3.373058, 0.2458538, -0.004850047];
    c_11["a6:m6"] = [4.161981, 0.2202737, -0.01718010];
    c_11["a7:m6"] = [2.4625223, 0.3231308, -0.02281072];

    c_11["a1:m4"] = [0.9851755, -0.02870704, 0.0];
    c_11["a2:m4"] = [-0.4737800,-0.001370344, 0.0];
    c_11["a3:m4"] = [0.7080799, 0.004575312, 0.0];
    c_11["a4:m4"] = [-1.239439, -0.03683605, 0.0];
    c_11["a5:m4"] = [-4.638467, 0.4418904, -0.01220292];
    c_11["a6:m4"] = [3.841835, 0.3277658, -0.02660275];
    c_11["a7:m4"] = [2.317342, 0.3912768, -0.03136223];

    c_22["a1:m6"] = [0.7898524, -0.02114115, 0.0];
    c_22["a2:m6"] = [-0.2998325, -0.001243977, 0.0];
    c_22["a3:m6"] = [0.7077103, 0.03583907, 0.0];
    c_22["a4:m6"] = [-0.8946857, -0.02473947, 0.0 ];
    c_22["a5:m6"] = [-2.958969, 0.2303358, -0.005226562];
    c_22["a6:m6"] = [4.348412, 0.1920321, 0.01496557];
    c_22["a7:m6"] = [2.205440, 0.2567027, -0.01861359];

    c_22["a1:m4"] = [0.9124518, -0.02398461, 0.0];
    c_22["a2:m4"] = [-0.4697184, -0.0007809681, 0.0];
    c_22["a3:m4"] = [1.031053, 0.004069668, 0.0];
    c_22["a4:m4"] = [-1.090782, -0.02413508, 0.0];
    c_22["a5:m4"] = [-4.127243, 0.4302667, -0.01352874];
    c_22["a6:m4"] = [4.059078, 0.2597379, -0.0216995];
    c_22["a7:m4"] = [2.086906, 0.2920310, -0.02560437];
    // a_ and d_  parameters from Bruno-2010, doi:10.1063/1.3495980
    // a_ + (spare cell) = g_ makes data structure consitent between equations
    _d_11[Ct.HeHep]= [38.6185, -3.1289, 0.06341];
    _d_11[Ct.HHp]  = [63.5437, -5.0093, 0.098797];

    a_22["He:He"] = [15.9887272, -4.3159308, 8.0, -15.52249873, 14.35288253, 4.0, 4.80040319,1.0];
    a_22["H2:H2"] = [27.54387526,-1.98253166,3.88885724,-12.91940775,0.3470796,8.72131306,-0.88296275,1.0];
    a_22["H:H"]   = [22.08948804,-1.85066626,8.50932055,-7.66943974,0.77454531,9.69545318,-0.62104466,1.0];
    a_22["H2:H"]  = [7.45048892,-1.4332616,9.59201391,-1.35066206,7.15859874,9.88881724, -1.39484886,1.0];
    a_22["He:H"]  = [11.79267184,-0.91854409,9.46049069,-4.15691291,76.18192207,-2.52045886,-2.67556855,1.0];
    a_22["He:H2"] = [0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]; //    
    a_22["H:Hp"]  = [46.68783791,-0.33303803,4.1021249,-1.85454858,14.98170959,8.86285119,-1.65616736,1.0];
    a_22["H2:Hp"] = [16.29263601,-1.32780746,10.00998032,-0.70566506,17.72506447,5.16010621,-2.45568134,1.0];
    a_22["He:Hp"] = [-3.94194326,0.73466954,8.79837656,-0.85969721,76.18192207,2.86197842,-2.85970495,1.0];
    a_22["H:Hep"] = [4.47457332,-0.2270375,9.11487604,-0.84002241,13.68381322,4.46462878,-3.36880863,1.0];
    a_22["H2:Hep"]= [0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]; //
    a_22["He:Hep"]= [30.01038539,-1.71614435,5.80735212,-4.49537445,2.18797641,9.3132696,-0.88202165,1.0];
    a_22["He:e"]= [11.02594479,-2.75226422,0.25837224,-1.80065315,1.11428478,36.06175647,-41.1176569,27.16469738];
    a_22["H2:e"]  = [8.44904482,-0.85941165,0.45608703,-0.032837008,0.88418654,3.62672485,9.45483919,1.67160237];
    a_22["H:e"]= [10.3344544,-1.44880911,12.08534341,-0.018616398,0.033072315,-6.45649277,9.15932646,-2.13494419];

    a_11["He:He"] = [14.66352146,-4.63010743,8.0,-15.4027924,16.26563367,4.0,5.06175281,1.0];
    a_11["H2:H2"] = [24.0084109,-1.61027553,3.88885724,-8.89043396,0.44260972,8.88408687,-1.05402226,1.0];
    a_11["H:H"]   = [15.09506044,-1.25710008,9.57839369,-3.80371463,0.98646613,9.25705877,-0.93611707,1.0];
    a_11["H2:H"]  = [12.4906397,-1.14704753,8.76515503,-3.52921496,0.32874932,12.77040465,-3.04802967,1.0];
    a_11["He:H"]  = [10.06955597,-0.78208832,8.90551185,-4.17119119,76.18192207,-2.53113293,-2.89309888,1.0];
    a_11["He:H2"] = [0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0];  //
    a_11["H:Hp"]  = [46.68783791,-0.33303803,4.2568677,-2.03851201,14.98170958,8.59618369,-1.65616736,1.0];
    a_11["H2:Hp"] = [21.26444082,-1.66425043,9.43954208,-1.18521519,17.72506447,5.16010621,-2.45568134,1.0];
    a_11["He:Hp"] = [2.71952379,-0.04746783,8.74213997,-0.88498112,91.77715465,2.49590845,-2.82208634,1.0];
    a_11["H:Hep"] = [7.26019648,-0.45830355,8.75289388,-1.0519608,7.88022406,5.84488923,-2.83633436,1.0];
    a_11["H2:Hep"]= [0.0,1.0,1.0,1.0,0.0,1.0,1.0,1.0]; //
    a_11["He:Hep"]= [23.80276836,-1.53417198,7.14997977,-3.74339421,2.241796,8.83793442,-0.990039,1.0];
    a_11["He:e"]  = [10.59136372,-1.52697125,0.021889243,-17.3791485,1.99019765,19.13183426,-3.87345716,52.39042];
    a_11["H2:e"]  = [7.61552567,-1.31152238,1.0876703,-0.085847868,0.60307406,5.01125522,9.07403916,1.92416429];
    a_11["H:e"] = [10.35291134,-1.58301162,12.45844036,-0.2328519,0.053662857,-5.34372929,9.35561752,-2.15463427];
}

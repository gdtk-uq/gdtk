/**
 * Interface and implementation for relaxation time calculations for energy exchange mechanisms.
 *
 *  - "Conservation Equations and Physical Models for Hypersonic Air Flows in Thermal and Chemical Nonequilibrium"
 *     Peter A. Gnoffo and Roop N. Gupta and Judy L. Shinn, NASA Technical Paper 2867, 1989
 *
 *  - "Review of Chemical-Kinetic Problems for Future NASA Missions, I: Earth Entries"
 *     Chul Park, Journal of Thermophysics and Heat Transfer, Volume 7, Number 3, July-Sept 1993
 *
 *  - "Modification of chemical-kinetic parameters for 11-air species in re-entry flows"
 *     Jae Gang Kim and Sung Min Jo, International Journal of Heat and Mass Transfer, Volume 169, 2021
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.relaxation_time;

import std.string;
import std.math;
import std.conv;
import std.stdio;

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;



interface RelaxationTime {
    RelaxationTime dup();
    @nogc number eval(in GasState gs, number[] molef, number[] numden);
}

class MillikanWhiteVT : RelaxationTime {
public:
    this(int q, double a, double b)
    {
        m_q = q;
        m_a = a;
        m_b = b;
    }

    this(lua_State *L, int q)
    {
        m_q = q;
        m_a = getDouble(L, -1, "a");
        m_b = getDouble(L, -1, "b");
    }

    MillikanWhiteVT dup()
    {
        return new MillikanWhiteVT(m_q, m_a, m_b);
    }

    @nogc
    override number eval(in GasState gs, number[] molef, number[] numden)
    {
        // If bath pressure is very small, then practically no relaxation
        // occurs from collisions with particle q.
        if (molef[m_q] <= SMALL_MOLE_FRACTION)
            return to!number(-1.0);
        // Set the bath pressure at that of the partial pressure of the 'q' colliders
        // and compute in atm for use in Millikan-White expression
        number pBath = molef[m_q]*gs.p/P_atm;
        number tau = (1.0/pBath) * exp(m_a * (pow(gs.T, -1./3.) - m_b) - 18.42);
        return tau;
    }

private:
    int m_q;
    double m_a;
    double m_b;
}

class ParkHTCVT : RelaxationTime {
/*
    High temperature correction taken from Gnoffo, 1989 (equations 56-58). Originally
    from Park, 1985.

    @author: Nick Gibbons
*/
    this(int p, int q, double sigma, double mu, RelaxationTime rt, GasModel gmodel)
    {
        m_p = p;
        m_q = q;
        m_sigma = sigma;
        m_mu = mu;
        m_rt = rt.dup();
        m_gmodel = gmodel;
    }

    this(lua_State *L, int p, int q, GasModel gmodel)
    {
        m_p = p;
        m_q = q;
        m_sigma = getDouble(L, -1, "sigma");

        double mq = gmodel.mol_masses[m_q]/Avogadro_number;
        double mp = gmodel.mol_masses[m_p]/Avogadro_number;
        m_mu = mq*mp/(mq+mp); // Reduced mass of the collision pair

        lua_getfield(L, -1, "submodel");
        m_rt = createRelaxationTime(L, p, q, gmodel);
        lua_pop(L, 1);
        m_gmodel = gmodel;
    }

    ParkHTCVT dup()
    {
        return new ParkHTCVT(m_p, m_q, m_sigma, m_mu, m_rt, m_gmodel);
    }

    @nogc
    override number eval(in GasState gs, number[] molef, number[] numden)
    {
        /*
        TODO: Park's original derivation is for the relaxation of a single
        species gas and it's not entirely clear what to do with the expression
        in the multi-species case.

        I've chosen to use the hard shell collision frequency expression, with
        the number density of the colliding particle as the relevant parameter.
        This is because we want the rate *per unit mass* of the relaxing
        species, although we should probably look into a more rigorous
        derivation of this so I can be sure that this is right.
        */
        if (molef[m_q] <= SMALL_MOLE_FRACTION)
            return to!number(-1.0);

        number tau_submodel = m_rt.eval(gs, molef, numden);

        number sigma = compute_cross_section(gs);
        number cs = sqrt(8.0*Boltzmann_constant*gs.T/pi/m_mu); // Mean particle velocity
        number tau_park = 1.0/(sigma*cs*numden[m_q]);        // Notice q is the colliding particle
        return tau_submodel + tau_park;
    }

protected:
    immutable double pi = to!double(PI);
    int m_p;
    int m_q;
    double m_sigma;
    double m_mu;
    RelaxationTime m_rt;
    GasModel m_gmodel;

    @nogc
    number compute_cross_section(ref const(GasState) gs) {
        return to!number(m_sigma);
    }
}

class ParkHTC2VT : ParkHTCVT {
/*
    High temperature correction taken Park 1993 (page 387).
    @author: Nick Gibbons
*/
    this(int p, int q, double sigma, double mu, RelaxationTime rt, GasModel gmodel)
    {
        super(p, q, sigma, mu, rt, gmodel);
    }

    this(lua_State *L, int p, int q, GasModel gmodel)
    {
        super(L, p, q, gmodel);
    }

    override ParkHTC2VT dup()
    {
        return new ParkHTC2VT(m_p, m_q, m_sigma, m_mu, m_rt, m_gmodel);
    }

protected:
    @nogc
    override number compute_cross_section(ref const(GasState) gs) {
        return m_sigma*(50e3/gs.T)*(50e3/gs.T);
    }
}

class KimHTCVT : ParkHTCVT {
/*
    High temperature correction taken from Kim 2021 (Table I).
    @author: Nick Gibbons
*/
    this(int p, int q, double sigma, double exponent, double mu, RelaxationTime rt, GasModel gmodel)
    {
        m_exponent = exponent;
        super(p, q, sigma, mu, rt, gmodel);
    }

    this(lua_State *L, int p, int q, GasModel gmodel)
    {
        m_exponent = getDouble(L, -1, "exponent");
        super(L, p, q, gmodel);
    }

    override KimHTCVT dup()
    {
        return new KimHTCVT(m_p, m_q, m_sigma, m_exponent, m_mu, m_rt, m_gmodel);
    }

protected:
    double m_exponent;

    @nogc
    override number compute_cross_section(ref const(GasState) gs) {
        return m_sigma*pow(gs.T, m_exponent);
    }
}

class Schwarzentruber : RelaxationTime {
public:
    this (double m_low, double n_low, double m_high, double n_high, int q) {
        _m_low = m_low;
        _n_low = n_low;
        _m_high = m_high;
        _n_high = n_high;
        _q = q;
    }

    this(lua_State *L, int q) {
        double m_low = getDouble(L, -1, "m_low");
        double n_low = getDouble(L, -1, "n_low");
        double m_high = getDouble(L, -1, "m_high");
        double n_high = getDouble(L, -1, "n_high");
        this(m_low, n_low, m_high, n_high, q);
    }

    RelaxationTime dup() {
        return new Schwarzentruber(_m_low, _n_low, _m_high, _n_high, _q);
    }

    @nogc
    number eval(in GasState gs, number[] molef, number[] numden) {
        if (molef[_q] <= SMALL_MOLE_FRACTION) {
            return to!number(-1.0);
        }

        number Tinv = pow(gs.T, -1./3.);
        number P_bath = molef[_q] * gs.p / P_atm;
        number low = _m_low*Tinv + _n_low;
        number high = _m_high*Tinv + _n_high;
        return (exp(low) + exp(high)) / P_bath;
    }

private:
    double _m_low, _n_low, _m_high, _n_high;
    int _q;
}


/* 
   Some functions used in SSH theory for V-V and V-T relaxation times
   These service the SSH relaxation time that follow
*/
@nogc number SSH_beta(double sigma, double epsilon, double mu, double delta_E, number T){
    double tmp_a = 2*epsilon/mu;
    number tmp_b = 3*Plancks_constant*mu/(PI*PI*sigma*Boltzmann_constant*T*delta_E);
    return pow(0.5 * pow(tmp_a, 9) * pow(tmp_b, 6), 1./19);
}

@nogc number SSH_rc_star(number beta, double sigma){
    number tmp =  pow(2*beta, 1./6.) / (1. + 2./19. * beta);
    return sigma * tmp;
}

@nogc number SSH_delta_star(number beta, double sigma) {
    number tmp_a = 12. / pow(2*beta, 1./6.);
    number tmp_b = 1. + 21./19. * beta;
    return tmp_a * tmp_b / sigma;
}

@nogc number SSH_chi(number alpha_pq, number T){
    return 0.5 * pow(alpha_pq / T, 1./3.);
}

@nogc number SSH_alpha(double mu, double delta_E, number delta_star){
    double num = 16.*pow(PI, 4) * mu * pow(delta_E, 2);
    number den = delta_star*delta_star * Plancks_constant*Plancks_constant * Boltzmann_constant;
    return num / den;
}

@nogc number SSH_A(number rc_star, double sigma){
    return rc_star * rc_star / (sigma * sigma);
}

@nogc number SSH_Z_0(number delta_star, double r){
    number tmp = delta_star * r;
    return tmp + (5./2)*(1. / (tmp*tmp));
}

@nogc number SSH_Z_V(double f_m, double mu_pp, double mu_pq, number alpha_pq, 
        double theta_v, double delta_E, int i)
{
    double tmp_a = f_m / ((i+1)*PI*PI);
    double tmp_b = mu_pp / mu_pq;
    number tmp_c = alpha_pq / theta_v;
    double tmp_d = Boltzmann_constant * theta_v / delta_E;
    return tmp_a * tmp_b * tmp_c * tmp_d * tmp_d;
}

@nogc number SSH_Z_T(double delta_E, number alpha_pq, number T){
    double tmp_a = PI*PI * sqrt(3/(2*PI));
    number tmp_b = delta_E / (Boltzmann_constant * alpha_pq);
    number tmp_c = T / alpha_pq;
    number tmp_d = 3./2. * pow(tmp_c, -1./3.) - delta_E / (2*Boltzmann_constant*T);
    return tmp_a * tmp_b*tmp_b * pow(tmp_c, 1./6.) * exp(tmp_d);
}

@nogc number SSH_Z_plus (double epsilon, number chi, number T){
    number kT = Boltzmann_constant * T;
    number tmp_a = epsilon * chi / (kT);
    number tmp_b = 16./(3*PI*PI) * epsilon / kT;
    number tmp_c = -4./PI * sqrt(tmp_a) - tmp_b;
    return exp(tmp_c);
}

@nogc number collision_frequency(number n_q, double sigma, double mu_pq, number T){
    // The collision frequency we are computing is the frequency of collision 
    // between a single particle of species p and any particle of q. 
    // Standard formulations of collision frequency count the number of collisions
    // between any particle of species p with any particle of q. Hence the 
    // collision frequency we are computing is Z_(pq) / numden[p]. Thus we only
    // need the number density of species q; the number density of species p will
    // cancel. 
    //
    // Note that this also means it doesn't matter if the colliding species are like
    // or unlike -- we are only tracking a single particle of species p. The discussion
    // in Rowan's thesis on apge 89-90 about the number density is not relevant given
    // the assumptions and cancellation disscused here.
    number tmp = 2 * PI * Boltzmann_constant * T / mu_pq;
    return 2 * n_q * sigma * sigma * sqrt(tmp);
}

/*
vibration-vibration relaxation time of Schwartz, Slawsky, and Herzfeld
See section 5.3.3 of Rowan's thesis for details

Robert Watt: Ported Eilmer3 code to Eilmer4 (19/12/22)
*/
class SSH_VV : RelaxationTime {
public:
    this(lua_State *L, int p, int q){
        _theta_v_p = getDouble(L, -1, "theta_v_p");
        _theta_v_q = getDouble(L, -1, "theta_v_q");
        _mu_pq = getDouble(L, -1, "mu_pq") / Avogadro_number;
        _mu_pp = getDouble(L, -1, "mu_pp") / Avogadro_number;
        _mu_qq = getDouble(L, -1, "mu_qq") / Avogadro_number;
        _delta_E = Boltzmann_constant * (_theta_v_p - _theta_v_q);
        _sigma = getDouble(L, -1, "sigma") * 1e-10; // angstrom -> m
        _epsilon = getDouble(L, -1, "epsilon") * Boltzmann_constant;
        _p = p;
        _q = q;
        _f_m_p = getDouble(L, -1, "f_m_p");
        _f_m_q = getDouble(L, -1, "f_m_q");
        _r_eq_p = getDouble(L, -1, "r_eq_p");
        _r_eq_q = getDouble(L, -1, "r_eq_q");
    }

    this(int p, int q, double mu_pq, double mu_pp, double mu_qq, 
            double sigma, double epsilon,
            double theta_p, double theta_q, double f_m_p, double f_m_q, 
            double r_eq_p, double r_eq_q) 
    {
        _p = p;
        _q = q;
        _mu_pq = mu_pq;
        _mu_pp = mu_pp;
        _mu_qq = mu_qq;
        _sigma = sigma;
        _epsilon = epsilon;
        _theta_v_p = _theta_v_p;
        _theta_v_q = _theta_v_q;
        _delta_E = Boltzmann_constant * (_theta_v_p - _theta_v_q);
        _f_m_p = f_m_p;
        _f_m_q = f_m_q;
        _r_eq_p = r_eq_p;
        _r_eq_q = r_eq_q;
    }

    RelaxationTime dup() {
        return new SSH_VV(_p, _q, _mu_pq, _mu_pp, _mu_qq, _sigma, _epsilon, 
                          _theta_v_p, _theta_v_q, _f_m_p, _f_m_q, _r_eq_p, _r_eq_q);
    }

    @nogc number eval(in GasState gs, number[] molef, number[] numden) {
        if (molef[_p] <= SMALL_MOLE_FRACTION || molef[_q] <= SMALL_MOLE_FRACTION)
            return to!number(-1.0);
        number T = gs.T;
        number n = numden[_q];
        number Z_COLL = collision_frequency(n, _sigma, _mu_pq, T);
        number P = _transition_probability(T);
        return 1./(Z_COLL * P) * exp(-(_theta_v_p - _theta_v_q) / T);
    }

private:
    double _mu_pp, _mu_qq, _mu_pq; // (kg) reduced mass
    double _theta_v_p, _theta_v_q; // (K) characteristic vib temperature of molecules
    double _delta_E;               // (J) energy difference between characteristic 
                                   //     vibration energies
    double _sigma, _epsilon;       // Lennard-Jones parameters (m, J)
    int _p, _q;                    // partcipating species
    double _f_m_p, _f_m_q;         // mass factor
    double _r_eq_p, _r_eq_q;       // Equilibrium separation of molecule

    @nogc number _transition_probability(number T){
        number beta = SSH_beta(_sigma, _epsilon, _mu_pq, _delta_E, T);
        number rc_star = SSH_rc_star(beta, _sigma);
        number delta_star = SSH_delta_star(beta, _sigma);
        number alpha_pq = SSH_alpha(_mu_pq, _delta_E, delta_star);
        number chi = SSH_chi(alpha_pq, T);

        number A = SSH_A(rc_star, _sigma);
        number Z_0_p = SSH_Z_0(delta_star, _r_eq_p);
        number Z_0_q = SSH_Z_0(delta_star, _r_eq_q);
        number Z_V_p = SSH_Z_V(_f_m_p, _mu_pp, _mu_pq, alpha_pq, _theta_v_p, _delta_E, 2);
        number Z_V_q = SSH_Z_V(_f_m_q, _mu_qq, _mu_pq, alpha_pq, _theta_v_q, _delta_E, 1);
        number Z_T = SSH_Z_T(_delta_E, alpha_pq, T);
        number Z_plus = SSH_Z_plus(_epsilon, chi, T);

        return A/(Z_0_p*Z_0_q*Z_V_p*Z_V_q*Z_T*Z_plus);
    }
}

/* 
vibration-translation relaxation times using SSH theory
See section 5.3.3 of Rowan's thesis for details
*/
class SSH_VT : RelaxationTime {
public:
    this(lua_State *L, int p, int q){
        _theta_v_p = getDouble(L, -1, "theta_v_p");
        _mu_pq = getDouble(L, -1, "mu_pq") / Avogadro_number;
        _mu_pp = getDouble(L, -1, "mu_pp") / Avogadro_number;
        _mu_qq = getDouble(L, -1, "mu_qq") / Avogadro_number;
        _delta_E = Boltzmann_constant * _theta_v_p;
        _sigma = getDouble(L, -1, "sigma") * 1e-10; // angstrom -> m
        _epsilon = getDouble(L, -1, "epsilon") * Boltzmann_constant;
        _p = p;
        _q = q;
        _f_m_p = getDouble(L, -1, "f_m_p");
        _r_eq_p = getDouble(L, -1, "r_eq_p");
    }

    this(int p, int q, double mu_pq, double mu_pp, double mu_qq, 
            double sigma, double epsilon, 
            double theta_p, double f_m_p, double r_eq_p) 
    {
        _p = p;
        _q = q;
        _mu_pq = mu_pq;
        _mu_pp = mu_pp;
        _mu_qq = mu_qq;
        _sigma = sigma;
        _epsilon = epsilon;
        _theta_v_p = _theta_v_p;
        _delta_E = Boltzmann_constant * _theta_v_p;
        _f_m_p = f_m_p;
        _r_eq_p = r_eq_p;
    }

    RelaxationTime dup() {
        return new SSH_VT(_p, _q, _mu_pq, _mu_pp, _mu_qq, _sigma, _epsilon, 
                          _theta_v_p, _f_m_p, _r_eq_p);
    }

    @nogc number eval(in GasState gs, number [] molef, number [] numden){
        if (molef[_p] <= SMALL_MOLE_FRACTION || molef[_q] <= SMALL_MOLE_FRACTION)
            return to!number(-1.0);

        number T = gs.T;
        number n = numden[_q];
        number Z_COLL = collision_frequency(n, _sigma, _mu_pq, T);
        number P = _transition_probability(T);
        return 1. / (Z_COLL * P * (1. - exp(-_theta_v_p / T)));
    }

private:
    double _mu_pp, _mu_qq, _mu_pq; // (kg) reduced mass
    double _theta_v_p;             // (K) characteristic vib temperature of molecules
    double _delta_E;               // (J) energy difference between characteristic 
                                   //     vibration energies
    double _sigma, _epsilon;       // Lennard-Jones parameters (m, J)
    double _r_eq_p, _r_eq_q;       // Equilibrium interatomic distance
    int _p, _q;                    // partcipating species
    double _f_m_p;                 // mass factor

    @nogc number _transition_probability(number T) {
        number beta = SSH_beta(_sigma, _epsilon, _mu_pq, _delta_E, T);
        number rc_star = SSH_rc_star(beta, _sigma);
        number delta_star = SSH_delta_star(beta, _sigma);
        number alpha_pq = SSH_alpha(_mu_pq, _delta_E, delta_star);
        number chi = SSH_chi(alpha_pq, T);

        number A = SSH_A(rc_star, _sigma);
        number Z_0_p = SSH_Z_0(delta_star, _r_eq_p);
        number Z_V_p = SSH_Z_V(_f_m_p, _mu_pp, _mu_pq, alpha_pq, _theta_v_p, _delta_E, 0);
        number Z_T = SSH_Z_T(_delta_E, alpha_pq, T);
        number Z_plus = SSH_Z_plus(_epsilon, chi, T);
        return A/(Z_0_p*Z_V_p*Z_T*Z_plus);
    }
}

class CandlerEV : RelaxationTime {
    // Free electron/vibration relaxation times. This uses the curve fit form given by Candler
    // (Eq. 2.7.14 of his thesis) for N2.

    this(double a_low, double b_low, double c_low, double a_high, double b_high, double c_high)
    {
        m_a_low = a_low; m_b_low = b_low; m_c_low = c_low;
        m_a_high = a_high; m_b_high = b_high; m_c_high = c_high;
    }

    this(lua_State *L)
    {
        m_a_low = getDouble(L, -1, "a_low");
        m_b_low = getDouble(L, -1, "b_low");
        m_c_low = getDouble(L, -1, "c_low");
        m_a_high = getDouble(L, -1, "a_high");
        m_b_high = getDouble(L, -1, "b_high");
        m_c_high = getDouble(L, -1, "c_high");
    }

    CandlerEV dup()
    {
        return new CandlerEV(m_a_low, m_b_low, m_c_low, m_a_high, m_b_high, m_c_high);
    }


    @nogc
    number eval(in GasState gs, number[] molef, number[] numden)
    {
        if (molef[$-1] < SMALL_MOLE_FRACTION)
            return to!number(-1.0);
        // assume free electron temperature is within the last temperature
        number Te = gs.T_modes[$-1];
        number log_Te = log10(Te);
        number p_e = gs.p_e;
        number log_pe_tau;
        
        // the curve fit fit by Candler doesn't have a continuous derivative at
        // 7000 K. So we blend it near 7000K to make sure the gradient is 
        // continuous.
        if (Te <= 6950.0) {
            log_pe_tau = _eval_low_temp(log_Te);
        }
        else if (Te >= 7050.0) {
            log_pe_tau = _eval_high_temp(log_Te);
        }
        else {
            number dT = (Te - 6950.0) / 100.0;
            number w0 = cos(PI/2 * dT.re) * cos(PI/2 * dT.re);
            number w1 = 1 - w0;
            log_pe_tau = _eval_low_temp(log_Te)*w0 + _eval_high_temp(log_Te)*w1;
        }
        // Divide by electron pressure converted to atmospheres
        return pow(10.0, log_pe_tau)/(p_e/101325.0);
    }

    @nogc
    number _eval_low_temp(number log_Te) {
        return m_a_low * log_Te * log_Te + m_b_low * log_Te + m_c_low;
    }

    @nogc
    number _eval_high_temp(number log_Te) {
        return m_a_high * log_Te * log_Te + m_b_high * log_Te + m_c_high;
    }

private:
    double m_a_low, m_b_low, m_c_low, m_a_high, m_b_high, m_c_high;
}

RelaxationTime createRelaxationTime(lua_State *L, int p, int q, GasModel gmodel)
{
    auto model = getString(L, -1, "model");
    switch (model) {
    case "Millikan-White":
	return new MillikanWhiteVT(L, q);
    case "Schwarzentruber":
        return new Schwarzentruber(L, q);
    case "ParkHTC":
	return new ParkHTCVT(L, p, q, gmodel);
    case "ParkHTC2":
	return new ParkHTC2VT(L, p, q, gmodel);
    case "KimHTC":
	return new KimHTCVT(L, p, q, gmodel);
    case "SSH_VT":
        return new SSH_VT(L, p, q);
    case "SSH-VV":
        return new SSH_VV(L, p, q);
    case "Candler":
        return new CandlerEV(L);
    default:
	string msg = format("The relaxation time model: %s is not known.", model);
	throw new Error(msg);
    }
}


/**
 * Interface and implementation for energy exchange mechanisms.
 *
 * References:
 *  - "Conservation Equations and Physical Models for Hypersonic Air Flows in Thermal and Chemical Nonequilibrium"
 *     Peter A. Gnoffo and Roop N. Gupta and Judy L. Shinn, NASA Technical Paper 2867, 1989
 *
 *  - "Effects of hydrogen impurities on shock structure and stability in ionizing monatomic gases. Part 1. Argon"
 *     I. I. Glass and W. S. Lio, Journal of Fluid Mechanics vol. 85, part 1, pp 55-77, 1978
*
 *  - "Theory and Validation of a Physically Consistent Couplied Vibration-Chemistry-Vibration Model"
 *    O. Knab and H. H. Fruhauf and E. W. Masserschmid, Journal of Thermophysics and Heat Transfer,
 *    Volume 9, Number 2, April-June 1995
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.energy_exchange_mechanism;

import std.format;
import std.math;
import std.stdio;
import std.conv;

import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;
import gas.thermo.energy_modes;

import kinetics.relaxation_time;
import kinetics.exchange_cross_section;
import kinetics.reaction_mechanism;
import kinetics.exchange_chemistry_coupling;

class EnergyExchangeMechanism {
    @property @nogc number tau() const { return m_tau; }
    @property @nogc int mode_p() const { return m_mode_p; }
    @property @nogc int mode_q() const { return m_mode_q; }
    @nogc void evalRelaxationTime(in GasState gs, number[] molef, number[] numden)
    {
        m_tau = mRT.eval(gs, molef, numden);
    }
    @nogc abstract number rate(in GasState gs, in GasState gseq, number[] molef, number[] numden, in ReactionMechanism rmech);

private:
    number m_tau;
    RelaxationTime mRT;
    GasModel mGmodel;
    int m_mode_p, m_mode_q;
}

class LandauTeller : EnergyExchangeMechanism {
public:

    this(lua_State *L, int mode_p, int mode_q, GasModel gmodel)
    {
        string pspecies = getString(L, -1, "p");
        string qspecies = getString(L, -1, "q");
        m_p = gmodel.species_index(pspecies);
        m_q = gmodel.species_index(qspecies);
        m_mode_p = mode_p;
        m_mode_q = mode_q;
        mGmodel = gmodel;
        lua_getfield(L, -1, "relaxation_time");
        mRT = createRelaxationTime(L, m_p, m_q, gmodel);
        lua_pop(L, 1);
    }

    this(int p, int q, int mode_p, int mode_q, RelaxationTime RT, GasModel gmodel)
    {
        m_p = p;
        m_q = q;
        m_mode_p = mode_p;
        m_mode_q = mode_q;
        mRT = RT.dup();
        mGmodel = gmodel;
    }

    @nogc
    override number rate(in GasState gs, in GasState gsEq, number[] molef, number[] numden, in ReactionMechanism rMech)
    {
        // If bath pressure is very small, then practically no relaxation
        // occurs from collisions with particle q.
        if (m_tau < 0.0) // signals very small mole fraction of q
            return to!number(0.0);
        number evStar = mGmodel.energyPerSpeciesInMode(gsEq, m_p, m_mode_p);
        number ev = mGmodel.energyPerSpeciesInMode(gs, m_p, m_mode_p);
        // NOTE 1. tau has already been weighted by colliding mole fractions.
        //         This is taken care of as part of by using bath pressure in
        //         calculation of relaxation time.
        // NOTE 2. massf scaling is applied here to convert J/s/kg-of-species-ip to J/s/kg-of-mixture
        return gs.massf[m_p] * (evStar - ev)/m_tau;
    }

private:
    int m_p;
    int m_q;
}

class ElectronExchangeET : EnergyExchangeMechanism {
    /*
        This rate expression accounts energy transfer between the electron
        thermal mode and the heavy particle thermal mode via elastic
        collisions, in flows where there is a separate electron temperature.

        Notes: A potential problem with this is that the reactor is updating u and uv,
        rather than T and Tv, so this routine isn't seeing the increments done during
        the RKF step.
        @author: Nick Gibbons
    */
    this(lua_State *L, int mode, GasModel gmodel)
    {
        string pspecies = getString(L, -1, "p");
        string qspecies = getString(L, -1, "q");
        m_p = gmodel.species_index(pspecies);
        m_q = gmodel.species_index(qspecies);
        m_e = m_p;
        m_mode_p = mode;
        m_mode_q = -1;
        mGmodel = gmodel;
        lua_getfield(L, -1, "exchange_cross_section");
        mCS = createExchangeCrossSection(L, m_e, m_q, m_mode_p);
        lua_pop(L, 1);

        m_me = gmodel.mol_masses[m_e]/Avogadro_number;
        m_mq = gmodel.mol_masses[m_q]/Avogadro_number;
        check_species_indices(m_e, m_q, gmodel);
    }

    this(int p, int q, int mode, ExchangeCrossSection CS, GasModel gmodel)
    {
        m_p = p;
        m_q = q;
        m_e = p;
        m_mode_p = mode;
        m_mode_q = -1;
        mCS = CS.dup();
        mGmodel = gmodel;

        m_me = gmodel.mol_masses[m_e]/Avogadro_number;
        m_mq = gmodel.mol_masses[m_q]/Avogadro_number;
        check_species_indices(m_e, m_q, gmodel);
    }

    @nogc
    override number rate(in GasState gs, in GasState gsEq, number[] molef, number[] numden, in ReactionMechanism rMech)
    {
    /*
        Electron Exchange rate from Glass and Liu, 1978, equation (11)

        Notes: A potential problem with this is that the reactor is updating u and uv,
        rather than T and Tv, so this routine isn't seeing the increments done during
        the RKF step.
    */
        immutable double pi = to!double(PI);
        immutable double six_times_sqrt2 = 6.0*sqrt(2.0);

        number ne = numden[m_e];
        number nq = numden[m_q];
        number Eq = Boltzmann_constant*gs.T;
        number Ee = Boltzmann_constant*gs.T_modes[m_mode_p];
        number Qeq = mCS(gs, numden); // Exchange collision cross section

        number rate = six_times_sqrt2*ne*nq*sqrt(m_me*Ee/pi)*Qeq/m_mq*(Eq - Ee);
        return rate/gs.rho; // Rate per unit mass of mixture
    }

    @nogc
    override void evalRelaxationTime(in GasState gs, number[] molef, number[] numden)
    {
        return;
    }


private:
    int m_p;
    int m_q;
    int m_e;
    ExchangeCrossSection mCS;
    GasModel mGmodel;
    double m_me, m_mq;

    @nogc final void check_species_indices(int e, int p, GasModel gmodel) {
        if (gmodel.species_name(e)!="e-") {
            throw new Error("Relaxer species in ElectronExchangeET must be an electron!");
        }
        if (gmodel.species_name(p)=="e-") {
            throw new Error("Collider species in ElectronExchangeET should not be an electron!");
        }
    }
}

class MarroneTreanorCV : EnergyExchangeMechanism {
    /*
        Model for the vibrational/chemistry coupling from Knab et al. 1995,
        following the formulation of Marrone and Treanor, 1963

        @author: Nick Gibbons
    */
    this(lua_State *L, int mode, GasModel gmodel)
    {
        m_mode_p = mode;
        m_mode_q = -1;
        mGmodel = gmodel;

        mReactionIdx = getInt(L, -1, "reaction_index");
        mSpeciesIdx  = gmodel.species_index(getString(L, -1, "p"));
        lua_getfield(L, -1, "coupling_model");
        mECC = createExchangeChemistryCoupling(L, gmodel, mode, mSpeciesIdx);
        lua_pop(L, 1);
    }

    this(int mode, GasModel gmodel, int reactionidx, int speciesidx, ExchangeChemistryCoupling ECC)
    {
        m_mode_p = mode;
        m_mode_q = -1;
        mGmodel = gmodel;
        mReactionIdx = reactionidx;
        mSpeciesIdx = speciesidx;
        mECC = ECC.dup();
    }

    @nogc
    override number rate(in GasState gs, in GasState gsEq, number[] molef, number[] numden, in ReactionMechanism rMech)
    {
    /*
        Based on Knab, 1995 equation 14. The reaction object handles the sign for us, so
        the production rate is always the amount of mSpeciesIdx appearing.
    */

    // Using Knab's formulation, this is mol/m3/s * J/mol -> J/m3/s
    number rate = rMech.production_rate(mReactionIdx, mSpeciesIdx)*mECC.Gappear(gs)
                - rMech.loss_rate(mReactionIdx, mSpeciesIdx)*mECC.Gvanish(gs);

    // Convert from J/m3/s (energy density of a specific oscillator) 
    // to J/kg/s (total vibrational energy of per unit mass of mixture)
    return rate/gs.rho;
    }

    @nogc
    override void evalRelaxationTime(in GasState gs, number[] molef, number[] numden)
    {
        // TODO: Maybe precompute some things here
        return;
    }

private:
    int mReactionIdx, mSpeciesIdx;
    GasModel mGmodel;
    ExchangeChemistryCoupling mECC;
}

class ThivetVV : EnergyExchangeMechanism {
    /* 
        Model for vibration exchanging with other vibration modes
    */
    this(lua_State *L, int mode_p, int mode_q, GasModel gmodel){
        mGmodel = gmodel; 
        string pspecies = getString(L, -1, "p");
        string qspecies = getString(L, -1, "q");
        m_p = gmodel.species_index(pspecies);
        m_q = gmodel.species_index(qspecies);
        m_mode_p = mode_p;
        m_mode_q = mode_q;

        lua_getfield(L, -1, "relaxation_time");
        mRT = createRelaxationTime(L, m_p, m_q, gmodel);
        lua_pop(L, 1);

        _theta_v_p = getDouble(L, -1, "theta_v_p");
        _theta_v_q = getDouble(L, -1, "theta_v_q");
        double theta_D_p = getDouble(L, -1, "theta_D_p");
        double theta_D_q = getDouble(L, -1, "theta_D_q");
        double R_p = getDouble(L, -1, "R_p");
        double R_q = getDouble(L, -1, "R_q");
        _e_hat_q = new HarmonicOscillator(_theta_v_q, R_q);
        _e_bar_p = new TruncatedHarmonicOscillator(_theta_v_p, theta_D_p, R_p);
        _e_bar_q = new TruncatedHarmonicOscillator(_theta_v_q, theta_D_q, R_q);
        _gm = gmodel;
    }

    @nogc override number rate(in GasState gs, in GasState gseq, number[] molef, 
                                    number[] numden, in ReactionMechanism rmech){
        // The relaxation time was set to -1 if there wasn't enough of either
        // species present
        if (m_tau < to!number(0.0)) {
            return to!number(0.0);
        }

        number e_p = _e_bar_p.energy(gs.T_modes[m_mode_p]);
        number e_q = _e_bar_q.energy(gs.T_modes[m_mode_q]);
        number e_hat_q = _e_hat_q.energy(gs.T);
        number e_bar_p = _e_bar_p.energy(gs.T);
        number e_bar_q = _e_bar_q.energy(gs.T);

        number tmp_a = e_q / e_hat_q * (e_bar_p - e_p);
        number tmp_b = e_p / e_hat_q * (e_bar_q - e_q);
        number exp_p = 1. - exp(-_theta_v_p / gs.T);
        number exp_q = 1. - exp(-_theta_v_q / gs.T);

        // multiply by massf[p] to convert J/kg-of-species-p to J/kg-of-mixture
        // The scaling by molef[q] is taken care of by using numden[q] in the
        // relaxation time instead of the total number density.
        return gs.massf[m_p] / m_tau * (exp_p / exp_q * tmp_a - tmp_b);

    }

private:
    int m_p, m_q;           // species indecies of participating species
    double _theta_v_p, _theta_v_q;
    InternalEnergy _e_hat_q;
    InternalEnergy _e_bar_p, _e_bar_q;
    GasModel _gm;
}

EnergyExchangeMechanism createEnergyExchangeMechanism(lua_State *L, int mode_p, int mode_q, GasModel gmodel)
{
    auto rateModel = getString(L, -1, "rate");
    switch (rateModel) {
    case "Landau-Teller":
        return new LandauTeller(L, mode_p, mode_q, gmodel);
    case "ElectronExchange":
        return new ElectronExchangeET(L, mode_p, gmodel);
    case "Marrone-Treanor":
        return new MarroneTreanorCV(L, mode_p, gmodel);
    case "Thivet-SSH":
        return new ThivetVV(L, mode_p, mode_q, gmodel); 
    default:
        string msg = format("The EE mechanism rate model: %s is not known.", rateModel);
        throw new Error(msg);
    }
}

/**
 * gas_model.d
 *
 * Contents: The gas model file has a number of parts.
 *   1. The GasModel base class for specifying how specific gas models should behave.
 *   2. Utility functions to transform mass-fraction and mole-fraction data arrays.
 *   3. Fill-in functions for gas model classes that don't implement some of
 *      the functions declared in the base class.
 * Note that GasState class, which specifies the storage arrangement for the data
 * defining a gas state is in a separate module.
 *
 * Authors: Peter J. and Rowan G.
 * Versions:
 *   2014-06-22, first cut, exploring the options.
 *   2015--2016, lots of experiments
 *   2017-01-06, introduce ChemicalReactor base class
 *   2017-11-10, move ThermochemicalReactor to its own module in kinetics package
 *   2018-06-02, adapted to complex numbers for Kyle
 */


//----------------------------------------------------------------------------------------
// PART 1. The GasModel class
//----------------------------------------------------------------------------------------

module gas.gas_model;

import std.conv;
import std.math;
import std.stdio;
import std.string;
import std.algorithm.searching : canFind;
import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import util.msg_service;
import gas.gas_state;
import gas.physical_constants;

immutable double SMALL_MOLE_FRACTION = 1.0e-15;
immutable double MIN_MASS_FRACTION = 1.0e-30;
immutable double MIN_MOLES = 1.0e-30;
immutable double T_MIN = 10.0;
immutable double T_MAX = 100000.0;
immutable double MASSF_ERROR_TOL = 1.0e-6;

immutable string[] equilibriumEnergyModeNames = ["equilibrium",
                                                 "transrotational",
                                                 "translational"];

class GasModelException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

class GasModel {
public:
    string type_str = ""; // To be replaced when constructing a specific model.
    @nogc @property bool is_plasma() const { return _is_plasma; }
    @nogc @property uint n_species() const { return _n_species; }
    @nogc @property uint n_heavy() const {
        uint nheavy = _n_species;
        if (_is_plasma) {
            if (_n_species <= 2) {
                throw new GasModelException("Not enough species for electrons to be separate.");
            }
            nheavy -= 1;
        }
        return nheavy;
    }
    @nogc @property uint n_modes() const { return _n_modes; }
    @nogc @property ref double[] mol_masses() { return _mol_masses; }
    @nogc @property ref double[] LJ_sigmas() { return _LJ_sigmas; }
    @nogc @property ref double[] LJ_epsilons() { return _LJ_epsilons; }
    @nogc @property ref double[] Le() { return _Lewis_numbers; }
    @nogc @property ref double[] electronic_energy() { return _electronic_energy; }
    @nogc @property ref double[] charge() { return _charge; }
    @nogc final string species_name(size_t i) const { return _species_names[i]; }
    final int species_index(string spName) const { return _species_indices.get(spName, -1); }
    @nogc final string energy_mode_name(int i) const { return _energy_mode_names[i]; }
    final int energy_mode_index(string modeName) const {
        if (equilibriumEnergyModeNames.canFind(modeName)) {
            return -1;
        }
        return _energy_mode_indices.get(modeName, -99);
    }

    void create_species_reverse_lookup()
    {
        _species_indices.clear;
        foreach (int isp; 0 .. _n_species) {
            _species_indices[_species_names[isp]] = isp;
        }
    }

    void create_energy_mode_reverse_lookup()
    {
        _energy_mode_indices.clear;
        foreach (int imode; 0 .. _n_modes) {
            _energy_mode_indices[_energy_mode_names[imode]] = imode;
        }
    }

    // Methods to be overridden.
    //
    // Although the following methods are not intended to alter their
    // GasModel object in a formal sense, they are not marked const.
    // The reason for this non-const-ness is that some GasModel classes
    // have private workspace that needs to be alterable.
    void update_thermo_from_pu(ref GasState Q) {
	return;
    }
    @nogc abstract void update_thermo_from_pT(ref GasState Q);
    @nogc abstract void update_thermo_from_rhou(ref GasState Q);
    @nogc abstract void update_thermo_from_rhoT(ref GasState Q);
    @nogc abstract void update_thermo_from_rhop(ref GasState Q);
    @nogc abstract void update_thermo_from_ps(ref GasState Q, number s);
    @nogc abstract void update_thermo_from_hs(ref GasState Q, number h, number s);
    @nogc abstract void update_sound_speed(ref GasState Q);
    @nogc abstract void update_trans_coeffs(ref GasState Q);
    // const void update_diff_coeffs(ref GasState Q) {}

    // Methods to be overridden.
    @nogc abstract number dudT_const_v(in GasState Q);
    @nogc abstract number dhdT_const_p(in GasState Q);
    @nogc abstract number dpdrho_const_T(in GasState Q);
    @nogc abstract number gas_constant(in GasState Q);
    @nogc number gas_constant(in GasState Q, int isp)
    {
        // For the single-species gases, provide a default implementation
        // but we need to be careful to override this for multi-component gases.
        return gas_constant(Q);
    }
    @nogc abstract number internal_energy(in GasState Q); // u+sum(u_modes)
    @nogc number internal_energy(in GasState Q, int isp)
    {
        // For the single-species gases, provide a default implementation
        // but we need to be careful to override this for multi-component gases.
        return internal_energy(Q);
    }
    @nogc number energyPerSpeciesInMode(in GasState Q, int isp, int imode)
    {
        // For single-species, single temperature gas a default implementation is:
        return internal_energy(Q);
    }
    //
    @nogc abstract number enthalpy(in GasState Q);
    @nogc number enthalpy(in GasState Q, int isp)
    {
        // For the single-species gases, provide a default implementation
        // but we need to be careful to override this for multi-component gases.
        return enthalpy(Q);
    }
    @nogc void enthalpies(in GasState Q, number[] hs)
    {
        foreach(isp; 0 .. _n_species){
            hs[isp] = enthalpy(Q, isp);
        }
    }
    @nogc number enthalpyPerSpeciesInMode(in GasState Q, int isp, int imode)
    {
        // For a single-species, one-temperature gas, provide a default implementation.
        // This should be overridden for multi-temperature gases.
        return to!number(0.0);
    }
    @nogc abstract number entropy(in GasState Q);
    @nogc number entropy(in GasState Q, int isp)
    {
        // For the single-species gases, provide a default implementation
        // but we need to be careful to override this for multi-component gases.
        return entropy(Q);
    }
    //
    @nogc number gibbs_free_energy(ref GasState Q, int isp)
    {
        number h = enthalpy(Q, isp);
        number s = entropy(Q, isp);
        number g = h - Q.T*s;
        return g;
     }
    @nogc void gibbs_free_energies(ref GasState Q, number[] gibbs_energies)
    {
        foreach(isp; 0 .. _n_species){
            number h = enthalpy(Q, isp);
            number s = entropy(Q, isp);
            gibbs_energies[isp] = h - Q.T*s;
        }
    }
    //
    @nogc final number Cv(in GasState Q) { return dudT_const_v(Q); }
    @nogc final number Cp(in GasState Q) { return dhdT_const_p(Q); }
    @nogc number Cp(in GasState Q, int isp)
    {
        // For the single-species gases, provide a default implementation
        // but we need to be careful to override this for multi-component gases.
        return Cp(Q);
    }
    @nogc final number R(in GasState Q)  { return gas_constant(Q); }
    @nogc final number gamma(in GasState Q) { return Cp(Q)/Cv(Q); }
    @nogc final number Prandtl(in GasState Q) { return Cp(Q)*Q.mu/Q.k; }
    @nogc
    final number molecular_mass(ref const(GasState) Q) const
    in {
        debug { assert(Q.massf.length == _mol_masses.length,
                       brokenPreCondition("Inconsistent array lengths.")); }
    }
    do {
        return mixture_molecular_mass(Q.massf, _mol_masses);
    }
    //
    @nogc
    final void massf2molef(ref const(GasState) Q, number[] molef) const
    in {
        debug { assert(Q.massf.length == molef.length,
                       brokenPreCondition("Inconsistent array lengths.")); }
    }
    do {
        gas.gas_model.massf2molef(Q.massf, _mol_masses, molef);
    }
    //
    @nogc
    final void molef2massf(const(number[]) molef, ref GasState Q) const
    in {
        debug { assert(Q.massf.length == molef.length,
                       brokenPreCondition("Inconsistent array lengths.")); }
    }
    do {
        gas.gas_model.molef2massf(molef, _mol_masses, Q.massf);
    }
    //
    @nogc
    final void massf2conc(ref GasState Q, number[] conc) const
    in {
        debug {
            assert(Q.massf.length == conc.length,
                   brokenPreCondition("Inconsistent array lengths."));
        }
    }
    do {
        if (isNaN(Q.rho) || (Q.rho <= 0.0)) {
            throw new GasModelException("Invalid density.");
        }
        foreach ( i; 0.._n_species ) {
            conc[i] = Q.massf[i]*Q.rho / _mol_masses[i];
            if ( conc[i] < MIN_MOLES ) conc[i] = 0.0;
        }
    }
    //
    @nogc
    final void conc2massf(const(number[]) conc, ref GasState Q) const
    in {
        debug {
            assert(Q.massf.length == conc.length,
                   brokenPreCondition("Inconsistent array lengths."));
        }
    }
    do {
        if (isNaN(Q.rho) || (Q.rho <= 0.0)) {
            throw new GasModelException("Invalid density.");
        }
        foreach ( i; 0.._n_species ) {
            Q.massf[i] = conc[i]*_mol_masses[i] / Q.rho;
            if ( Q.massf[i] < MIN_MASS_FRACTION ) Q.massf[i] = 0.0;
        }
    }
    //
    @nogc
    final void rates2source(number[] rates, number[] source) const
    do {
        foreach ( i; 0.._n_species ) {
            source[i] = rates[i]*_mol_masses[i];
        }
    }
    //
    @nogc
    final void massf2numden(ref const(GasState) Q, number[] numden) const
    in {
        debug { assert(Q.massf.length == numden.length,
                       brokenPreCondition("Inconsistent array lengths.")); }
    }
    do {
        foreach ( i; 0.._n_species ) {
            numden[i] = Avogadro_number*Q.massf[i]*Q.rho / _mol_masses[i];
        }
    }
    //
    @nogc
    final void numden2massf(const(number[]) numden, ref GasState Q) const
    in {
        debug { assert(Q.massf.length == numden.length,
                       brokenPreCondition("Inconsistent array lengths.")); }
    }
    do {
        foreach ( i; 0.._n_species ) {
            Q.massf[i] = numden[i]*_mol_masses[i] / (Avogadro_number*Q.rho);
            if ( Q.massf[i] < MIN_MASS_FRACTION ) Q.massf[i] = 0.0;
        }
    }
    //
    @nogc
    void balance_charge(ref GasState Q)
    {
        // The flow solver may call upon this function for plasmas,
        // when it wants the electron fraction updated so that it
        // matches the ion fraction(s).
        //
        // If relevant, the derived GasModel should know what to do.
        //
        if (_is_plasma) {
            // ... such that we should never arrive here.
            throw new Error("Gas model does not set the electron mass-fraction.");
        }
    }
    @nogc
    void binary_diffusion_coefficients(ref const(GasState) Q, ref number[][] D)
    {
        // Calling this (optional) method without overriding it should be an error
        // But the report to the user should depend on what they are trying to do
        if (_n_species<=1) {
            throw new Error("Species diffusion invalid for single species gas.");
        } else {
            throw new Error("Gas model has no binary diffusion implementation.");
        }
    }
    @nogc
    void minimum_mixture_energy(ref GasState Q)
    {
        // NOTE: This method changes the state in Q and puts the minimum
        //       energy values in u and u_modes.
        // Assumption: perfect gas behaviour.
        // Here we set pressure to an arbitrary value. This assumes that
        // pressure has no part in the caloric equation of state.
        // In other words, we assume no real gas effects.
        Q.p = P_atm;
        Q.T = T_MIN;
        foreach (ref T; Q.T_modes) T = T_MIN;
        update_thermo_from_pT(Q);
    }

protected:
    // Default to non-plasma gas model, where all species are treated alike.
    // The quasi-neutral plasma model assumes that the last species is the electron,
    // and that the number of electrons balances the numbers of ions.
    bool _is_plasma = false;

    // The following data need to be properly initialized by the derived class.
    uint _n_species;
    uint _n_modes;
    string[] _species_names;
    int[string] _species_indices;
    string[] _energy_mode_names;
    int[string] _energy_mode_indices;
    double[] _mol_masses; // kg/mol
    double[] _LJ_sigmas;
    double[] _LJ_epsilons;
    double[] _Lewis_numbers;
    double[] _electronic_energy;
    double[] _charge;
} // end class GasModel


//----------------------------------------------------------------------------------------
// PART 2. Utility functions to transform mass-fraction and mole-fraction data arrays.
//----------------------------------------------------------------------------------------

@nogc
void scale_mass_fractions(ref number[] massf, double tolerance=0.0,
                          double assert_error_tolerance=0.1)
{
    pragma(inline, true);
    auto my_nsp = massf.length;
    if (my_nsp == 1) {
        // Single species, always expect massf[0]==1.0, so we can take a short-cut.
        if (fabs(massf[0] - 1.0) > assert_error_tolerance) {
            throw new GasModelException("Single species mass fraction far from 1.0.");
        }
        massf[0] = 1.0;
    } else {
        // Multiple species, do the full job.
        number massf_sum = 0.0;
        foreach(isp; 0 .. my_nsp) {
            massf[isp] = massf[isp] >= 0.0 ? massf[isp] : to!number(0.0);
            massf_sum += massf[isp];
        }
        if (fabs(massf_sum - 1.0) > assert_error_tolerance) {
            string msg = "Sum of species mass fractions far from 1.0";
            debug {
                msg ~= format(": fabs(massf_sum - 1.0) = %s \n", fabs(massf_sum - 1.0));
                msg ~= format("  assert_error_tolerance = %s \n", assert_error_tolerance);
                msg ~= format("  tolerance = %s \n", tolerance);
                msg ~= "  massf = [";
                foreach (mf; massf) {msg ~= format(" %e", mf); }
                msg ~= "]\n";
            }
            throw new GasModelException(msg);
        }
        if ( fabs(massf_sum - 1.0) > tolerance ) {
            foreach(isp; 0 .. my_nsp) massf[isp] /= massf_sum;
        }
    }
    return;
} // end scale_mass_fractions()

@nogc pure number mass_average(in GasState Q, in number[] phi)
in {
    assert(Q.massf.length == phi.length, "Inconsistent array lengths.");
}
do {
    number result = 0.0;
    foreach ( i; 0..Q.massf.length ) result += Q.massf[i] * phi[i];
    return result;
}
version(complex_numbers) {
    // We want to retain the flavour with double numbers passed in.
    @nogc pure number mass_average(in GasState Q, in double[] phi)
        in {
            assert(Q.massf.length == phi.length, "Inconsistent array lengths.");
        }
    do {
        number result = 0.0;
        foreach ( i; 0..Q.massf.length ) result += Q.massf[i] * phi[i];
        return result;
    }
}

@nogc pure number mole_average(in number[] molef, in number[] phi)
in {
    assert(molef.length == phi.length, "Inconsistent array lengths.");
}
do {
    number result = 0.0;
    foreach ( i; 0..molef.length ) result += molef[i] * phi[i];
    return result;
}
version(complex_numbers) {
    // We want to retain the flavour with double numbers.
    @nogc pure number mole_average(in number[] molef, in double[] phi)
        in {
            assert(molef.length == phi.length, "Inconsistent array lengths.");
        }
    do {
        number result = 0.0;
        foreach ( i; 0..molef.length ) result += molef[i] * phi[i];
        return result;
    }
}

@nogc pure number mixture_molecular_mass(in number[] massf, in double[] mol_masses)
in {
    assert(massf.length == mol_masses.length, "Inconsistent array lengths.");
}
do {
    number M_inv = 0.0;
    foreach (i; 0 .. massf.length) M_inv += massf[i] / mol_masses[i];
    return 1.0/M_inv;
}

@nogc void massf2molef(in number[] massf, in double[] mol_masses, number[] molef)
in {
    assert(massf.length == mol_masses.length, "Inconsistent array lengths.");
    assert(massf.length == molef.length, "Inconsistent array lengths.");
}
do {
    number mixMolMass = mixture_molecular_mass(massf, mol_masses);
    foreach ( i; 0..massf.length ) molef[i] = massf[i] * mixMolMass / mol_masses[i];
}

@nogc void molef2massf(in number[] molef, in double[] mol_masses, number[] massf)
in {
    assert(massf.length == mol_masses.length, "Inconsistent array lengths.");
    assert(massf.length == molef.length, "Inconsistent array lengths.");
}
do {
    number mixMolMass = mole_average(molef, mol_masses);
    foreach ( i; 0..massf.length ) massf[i] = molef[i] * mol_masses[i] / mixMolMass;
}


//----------------------------------------------------------------------------------------
// PART 3. Fill-in functions for gas models that don't define all functions
//         specified in the base class GasModel
//----------------------------------------------------------------------------------------

/* The following functions:
   update_thermo_state_pT(), update_thermo_state_rhoT(), update_thermo_state_rhop()
   are for updating the thermo state from when  the gas model does not  have a method
   for doing so with those variables, but does have a defined method for
   update_thermo_from_rhou(). (e.g. in the UniformLUT class)
   A guess is made for rho & e and that guess is iterated using the  Newton-Raphson
   method.
   A GasModel object is a function parameter so that the update method from rho,e for
   that gas model can be called.
   The funcions:
   update_thermo_state_ps(), and update_thermo_state_hs() actually iterate on the update
   update_thermo_from_pT(), as  called from the gas model (though the p,T update can be
   defined as the function defined here that itself iterates on the update method for
   rho,e)
*/

immutable MAX_RELATIVE_STEP = 0.1;
immutable MAX_STEPS = 50;

@nogc void update_thermo_state_pT(GasModel gmodel, ref GasState Q)
{
    number drho, rho_old, rho_new, e_old, e_new, de;
    number drho_sign, de_sign;
    number Cv_eff, R_eff, T_old;
    number fp_old, fT_old, fp_new, fT_new;
    number dfp_drho, dfT_drho, dfp_de, dfT_de, det;
    int converged, count;

    number p_given = Q.p;
    number T_given = Q.T;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.

    number fT_tol = 1.0e-6 * T_given;
    number fp_tol = 1.0e-6 * p_given;
    number fp_tol_fail = 0.02 * p_given;
    number fT_tol_fail = 0.02 * T_given;

    Q.rho = 1.0; // kg/m**3
    Q.u = 2.0e5; // J/kg
    // Get an idea of the gas properties by calling the original
    // equation of state with some starting values for density
    // and internal energy.
    gmodel.update_thermo_from_rhou(Q);
    T_old = Q.T;
    R_eff = Q.p / (Q.rho * T_old);
    de = 0.01 * Q.u;
    Q.u += de;
    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Getting gas properties, failed in update_thermo_state_pT.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }
    Cv_eff = de / (Q.T - T_old);
    // Now, get a better guess for the appropriate density and
    // internal energy.
    e_old = Q.u + (T_given - Q.T) * Cv_eff;
    rho_old = p_given / (R_eff * T_given);
    //
    // Evaluate the error functions using this guess.
    Q.rho = rho_old;
    Q.u = e_old;
    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess at iteration 2 failed in update_thermo_state_pT.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }
    fp_old = p_given - Q.p;
    fT_old = T_given - Q.T;
    //
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fp_old) < fp_tol) && (fabs(fT_old) < fT_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
        // Perturb first dimension to get derivatives.
        rho_new = rho_old * 1.001;
        e_new = e_old;
        Q.rho = rho_new;
        Q.u = e_new;
        try { gmodel.update_thermo_from_rhou(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration failed in update_thermo_state_pT.";
            debug {
                msg ~= format("\nIteration %s failed at call A in %s\n", count, __FUNCTION__);
                msg ~= format("Exception message from update_thermo_from_rhou() was:\n\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }
        fp_new = p_given - Q.p;
        fT_new = T_given - Q.T;
        dfp_drho = (fp_new - fp_old) / (rho_new - rho_old);
        dfT_drho = (fT_new - fT_old) / (rho_new - rho_old);
        // Perturb other dimension to get derivatives.
        rho_new = rho_old;
        e_new = e_old * 1.001;
        Q.rho = rho_new;
        Q.u = e_new;

        try { gmodel.update_thermo_from_rhou(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration failed in update_thermo_state_pT.";
            debug {
                msg ~= format("\nIteration %s failed at call B in %", count, __FUNCTION__);
                msg ~= format("Exception message from update_thermo_from_rhou() was:\n\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }

        fp_new = p_given - Q.p;
        fT_new = T_given - Q.T;
        dfp_de = (fp_new - fp_old) / (e_new - e_old);
        dfT_de = (fT_new - fT_old) / (e_new - e_old);
        det = dfp_drho * dfT_de - dfT_drho * dfp_de;
        if( fabs(det) < 1.0e-12 ) {
            string msg = "Iteration failed in update_thermo_state_pT, nearly zero determinant.";
            debug { msg ~= format("    Nearly zero determinant, det = ", det); }
            throw new GasModelException(msg);
        }
        drho = (-dfT_de * fp_old + dfp_de * fT_old) / det;
        de = (dfT_drho * fp_old - dfp_drho * fT_old) / det;
        if( fabs(drho) > MAX_RELATIVE_STEP * rho_old ) {
            // move a little toward the goal
            drho_sign = (drho > 0.0 ? 1.0 : -1.0);
            drho = drho_sign * MAX_RELATIVE_STEP * rho_old;
        }
        if( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
            // move a little toward the goal
            de_sign = (de > 0.0 ? 1.0 : -1.0);
            de = de_sign * MAX_RELATIVE_STEP * e_old;
        }
        rho_old += drho;
        e_old += de;
        // Make sure of consistent thermo state.
        Q.rho = rho_old;
        Q.u = e_old;
        try { gmodel.update_thermo_from_rhou(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration failed in update_thermo_state_pT.";
            debug {
                msg ~= format("\nIteration %d\n", count);
                msg ~= format("Exception message from update_thermo_from_rhou() was:\n\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }
        // Prepare for next iteration.
        fp_old = p_given - Q.p;
        fT_old = T_given - Q.T;
        converged = (fabs(fp_old) < fp_tol) && (fabs(fT_old) < fT_tol);
        ++count;
    } // end while

    if (count >= MAX_STEPS) {
        string msg = "Iterations did not converge in update_thermo_state_pT.";
        debug {
            msg ~= format("\n    fp_old = %g, fT_old = %g\n", fp_old, fT_old);
            msg ~= format("    p_given = %.10s, T_given, %.5s\n", p_given, T_given);
            msg ~= "  Supplied Q:" ~ Q.toString();
        }
        throw new GasModelException(msg);
    }

    if( (fabs(fp_old) > fp_tol_fail) || (fabs(fT_old) > fT_tol_fail) ) {
        string msg = "Iterations failed badly in update_thermo_state_pT.";
        debug {
            msg ~= format("\n    p_given = %.10s, T_given, %.5s\n", p_given, T_given);
            msg ~= format("    fp_old = %g, fT_old = %g\n", fp_old, fT_old);
            msg ~= "  Supplied Q:" ~ Q.toString();
        }
        throw new GasModelException(msg);
    }
} // end update_thermo_state_pT()

@nogc void update_thermo_state_rhoT(GasModel gmodel, ref GasState Q)
{
    // This method can be applied to single-species models only
    number e_old, e_new, de, tmp, de_sign;
    number Cv_eff, T_old;
    number dfT_de, fT_old, fT_new;
    int converged, count;

    number rho_given = Q.rho;
    number T_given = Q.T;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    number fT_tol = 1.0e-6 * T_given;
    number fT_tol_fail = 0.02 * T_given;

    // Get an idea of the gas properties by calling the original
    // equation of state with some dummy values for density
    // and internal energy.

    Q.rho = rho_given; // kg/m**3
    Q.u = 2.0e5; // J/kg
    gmodel.update_thermo_from_rhou(Q);

    T_old = Q.T;
    de = 0.01 * Q.u;
    Q.u += de;

    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess at iteration 0 failed in update_thermo_state_rhoT.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }

    Cv_eff = de / (Q.T - T_old);
    // Now, get a better guess for the appropriate density and internal energy.
    e_old = Q.u + (T_given - Q.T) * Cv_eff;
    // Evaluate state variables using this guess.
    Q.rho = rho_given;
    Q.u = e_old;

    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess at iteration 1 failed in update_thermo_state_rhoT.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }

    fT_old = T_given - Q.T;
    // Perturb to get derivative.
    e_new = e_old * 1.001;
    Q.rho = rho_given;
    Q.u = e_new;

    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess at iteration 2 failed in update_thermo_state_rhoT.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }

    fT_new = T_given - Q.T;
    dfT_de = (fT_new - fT_old) / (e_new - e_old);

    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fT_new) < fabs(fT_old) ) {
        tmp = fT_new; fT_new = fT_old; fT_old = tmp;
        tmp = e_new; e_new = e_old; e_old = tmp;
    }
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fT_old) < fT_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
        de = -fT_old / dfT_de;
        if ( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
            // move a little toward the goal
            de_sign = (de > 0.0 ? 1.0 : -1.0);
            de = de_sign * MAX_RELATIVE_STEP * fabs(e_old);
        }
        e_new = e_old + de;
        Q.rho = rho_given;
        Q.u = e_new;
        try { gmodel.update_thermo_from_rhou(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration failed in update_thermo_state_rhoT.";
            debug {
                msg ~= format("\nIteration %d", count);
                msg ~= format("Exception message from update_thermo_from_rhou() was:\n\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }
        fT_new = T_given - Q.T;
        dfT_de = (fT_new - fT_old) / (e_new - e_old);
        // Prepare for the next iteration.
        ++count;
        fT_old = fT_new;
        e_old = e_new;
        converged = fabs(fT_old) < fT_tol;
    }   // end while
    // Ensure that we have the current data for all EOS variables.
    Q.rho = rho_given;
    Q.u = e_old;

    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Function update_thermo_state_rhoT failed after finishing iterations.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }

     if (count >= MAX_STEPS) {
        string msg = "Iterations did not converge in update_thermo_state_rhoT.";
        debug {
            msg ~= format("\n    fT_old = %g\n", fT_old);
            msg ~= format("    rho_given = %.5s, T_given, %.5s\n", rho_given, T_given);
            msg ~= "  Supplied Q:" ~ Q.toString;
        }
        throw new GasModelException(msg);
    }
    if (fabs(fT_old) > fT_tol_fail) {
        string msg = "Iterations failed badly in update_thermo_state_rhoT.";
        debug {
            msg ~= format("    Iterations failed badly.\n");
            msg ~= format("    rho_given = %.5s, T_given, %.5s\n", rho_given, T_given);
            msg ~= "  Supplied Q:" ~ Q.toString();
        }
        throw new GasModelException(msg);
    }
} // end update_thermo_state_rhoT()

@nogc void update_thermo_state_rhop(GasModel gmodel, ref GasState Q)
{
    number e_old, e_new, de, dedp, tmp, de_sign;
    number p_old;
    number dfp_de, fp_old, fp_new;
    int converged, count;

    number rho_given = Q.rho;
    number p_given = Q.p;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    number fp_tol = 1.0e-6 * p_given;
    number fp_tol_fail = 0.02 * p_given;

    // Get an idea of the gas properties by calling the original
    // equation of state with some dummy values for density
    // and internal energy.
    Q.rho = rho_given; // kg/m**3
    Q.u = 2.0e5; // J/kg
    gmodel.update_thermo_from_rhou(Q);
    p_old = Q.p;
    de = 0.01 * Q.u;
    Q.u += de;

    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess at iteration 0 failed in update_thermo_state_rhop.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }

    dedp = de / (Q.p - p_old);
    // Now, get a better guess for the appropriate internal energy.
    e_old = Q.u + (p_given - Q.p) * dedp;
    //     printf( "Initial guess e_old= %g dedp= %g\n", e_old, dedp );
    // Evaluate state variables using this guess.
    Q.rho = rho_given;
    Q.u = e_old;

    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess at iteration 1 failed in update_thermo_state_rhop.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }

    fp_old = p_given - Q.p;
    // Perturb to get derivative.
    e_new = e_old * 1.001;
    Q.rho = rho_given;
    Q.u = e_new;

    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess at iteration 2 failed in update_thermo_state_rhop.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }

    fp_new = p_given - Q.p;
    dfp_de = (fp_new - fp_old) / (e_new - e_old);

    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fp_new) < fabs(fp_old) ) {
        tmp = fp_new; fp_new = fp_old; fp_old = tmp;
        tmp = e_new; e_new = e_old; e_old = tmp;
    }
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fp_old) < fp_tol);
    count = 0;
    while (!converged && count < MAX_STEPS) {
        de = -fp_old / dfp_de;
        if ( fabs(de) > MAX_RELATIVE_STEP * e_old ) {
            // move a little toward the goal
            de_sign = (de > 0.0 ? 1.0 : -1.0);
            de = de_sign * MAX_RELATIVE_STEP * fabs(e_old);
        }
        e_new = e_old + de;
        Q.rho = rho_given;
        Q.u = e_new;

        try { gmodel.update_thermo_from_rhou(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration failed in update_thermo_state_rhop.";
            debug {
                msg ~= format("\nIteration %d", count);
                msg ~= format("Exception message from update_thermo_from_rhou() was:\n\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }

        fp_new = p_given - Q.p;
        dfp_de = (fp_new - fp_old) / (e_new - e_old);
        // Prepare for next iteration.
        ++count;
        fp_old = fp_new;
        e_old = e_new;
        converged = fabs(fp_old) < fp_tol;
    }   // end while
    // Ensure that we have the current data for all EOS variables.
    Q.rho = rho_given;
    Q.u = e_old;

    try { gmodel.update_thermo_from_rhou(Q); }
    catch (Exception caughtException) {
        string msg = "Function update_thermo_state_rhop failed after finishing iterations.";
        debug {
            msg ~= format("\nException message from update_thermo_from_rhou() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }

      if (count >= MAX_STEPS) {
        string msg = "Iterations did not converge in update_thermo_state_rhop.";
        debug {
            msg ~= format("\n    fp_old = %g, e_old = %g\n", fp_old, e_old);
            msg ~= format("    rho_given = %.5s, p_given, %.8s\n", rho_given, p_given);
            msg ~= "  Supplied Q:" ~ Q.toString;
        }
        throw new GasModelException(msg);
    }

    if (fabs(fp_old) > fp_tol_fail) {
        string msg = "Iterations failed badly in update_thermo_state_rhop.";
        debug {
            msg ~= format("\n    rho_given = %.5s, T_given, %.8s\n", rho_given, p_given);
            msg ~= "  Supplied Q:" ~ Q.toString();
        }
        throw new GasModelException(msg);
    }
} // end update_thermo_state_rhop()

@nogc void update_thermo_state_ps(GasModel gmodel, ref GasState Q, number s)
{
    number T_old, T_new, dT, tmp, dT_sign;
    number dfs_dT, fs_old, fs_new;
    int converged, count;

    number s_given = s;
    number p_given = Q.p;

    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    number fs_tol = 1.0e-6 * s_given;
    number fs_tol_fail = 0.02 * s_given;
    //
    // Fill in the thermo state assuming that T is a good guess.
    // Note that pressure and mass fractions must be correctly set but,
    // maybe T is not already set.
    if (isNaN(Q.T)) {
        // Since a guess for T is required, we have to set a value.
        Q.T = 1000.0;
        // [TODO] PJ 2023-08-07
        // We set an arbitrary number for the moment but we should do better,
        // maybe setting up a list of candidates and selecting the best.
    }
    T_old = Q.T;
    try { gmodel.update_thermo_from_pT(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess for T: failed in update_thermo_state_ps.";
        debug {
            msg ~= format("\nException message from update_thermo_from_pT() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }
    number s_old = gmodel.entropy(Q);
    fs_old = s_given - s_old;
    //
    // Perturb T to get a derivative estimate
    T_new = T_old * 1.001;
    Q.T = T_new;
    try { gmodel.update_thermo_from_pT(Q); }
    catch (Exception caughtException) {
        string msg = "Starting guess perturbed T: failed in update_thermo_state_ps.";
        debug {
            msg ~= format("\nException message from update_thermo_from_pT() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }
    number s_new = gmodel.entropy(Q);
    fs_new = s_given - s_new;
    dfs_dT = (fs_new - fs_old)/(T_new - T_old);
    // At the start of iteration, we want *_old to be the best guess.
    if ( fabs(fs_new) < fabs(fs_old) ) {
        tmp = fs_new; fs_new = fs_old; fs_old = tmp;
        tmp = s_new; s_new = s_old; s_old = tmp;
    }
    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fs_old) < fs_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
        dT = -fs_old / dfs_dT;
        if ( fabs(dT) > MAX_RELATIVE_STEP * T_old ) {
            // move a little toward the goal
            dT_sign = (dT > 0.0 ? 1.0 : -1.0);
            dT = dT_sign * MAX_RELATIVE_STEP * fabs(T_old);
        }
        T_new = T_old + dT;
        Q.T = T_new;
        try { gmodel.update_thermo_from_pT(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration failed in update_thermo_state_ps.";
            debug {
                msg ~= format("\nIteration %d", count);
                msg ~= format("\nException message from update_thermo_from_pT() was:\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }
        s_new = gmodel.entropy(Q);
        fs_new = s_given - s_new;
        dfs_dT = (fs_new - fs_old) / (T_new - T_old);
        // Prepare for next iteration.
        ++count;
        fs_old = fs_new;
        T_old = T_new;
        converged = (fabs(fs_old) < fs_tol);
    }   // end while
    //
    // Ensure that we have the current data for all EOS variables.
    Q.T = T_old;
    try { gmodel.update_thermo_from_pT(Q); }
    catch (Exception caughtException) {
        string msg = "Function update_thermo_state_ps failed after finishing iterations.";
        debug {
            msg ~= format("\nException message from update_thermo_from_pT() was:\n\n");
            msg ~= to!string(caughtException);
        }
        throw new GasModelException(msg);
    }
    if ( count >= MAX_STEPS ) {
        string msg = "Iterations did not converge in update_thermo_state_ps.";
        debug {
            msg ~= format("\n    fs_old = %g  count=%d MAX_STEPS=%d\n", fs_old, count, MAX_STEPS);
            msg ~= format("    p_given = %.8s, s_given, %.5s\n", p_given, s_given);
            msg ~= "  Supplied Q:" ~ Q.toString;
        }
        throw new GasModelException(msg);
    }

    if ( fabs(fs_old) > fs_tol_fail ) {
        string msg = "Iterations failed badly in update_thermo_state_ps.";
        debug {
            msg ~= format("\n    p_given = %.8s, s_given, %.5s\n", p_given, s_given);
            msg ~= "  Supplied Q:" ~ Q.toString();
        }
        throw new GasModelException(msg);
    }
} // end update_thermo_state_ps()

@nogc void update_thermo_state_hs(GasModel gmodel, ref GasState Q, number h, number s)
{
    number dp, p_old, p_new, T_old, T_new, dT;
    number dp_sign, dT_sign;
    number fh_old, fs_old, fh_new, fs_new;
    number dfh_dp, dfs_dp, dfh_dT, dfs_dT, det;
    int converged, count;

    number h_given = h;
    number s_given = s;
    // When using single-sided finite-differences on the
    // curve-fit EOS functions, we really cannot expect
    // much more than 0.1% tolerance here.
    // However, we want a tighter tolerance so that the starting values
    // don't get shifted noticeably.
    number fh_tol = 1.0e-6 * h_given;
    number fs_tol = 1.0e-6 * s_given;
    number fh_tol_fail = 0.02 * h_given;
    number fs_tol_fail = 0.02 * s_given;

    // Use current gas state as guess
    p_old = Q.p;
    T_old = Q.T;
    number h_new = gmodel.enthalpy(Q);
    number s_new = gmodel.entropy(Q);
    fh_old = h_given - h_new;
    fs_old = s_given - s_new;

    // Update the guess using Newton iterations
    // with the partial derivatives being estimated
    // via finite differences.
    converged = (fabs(fh_old) < fh_tol) && (fabs(fs_old) < fs_tol);
    count = 0;
    while ( !converged && count < MAX_STEPS ) {
        // Perturb first dimension to get derivatives.
        p_new = p_old * 1.001;
        T_new = T_old;
        Q.p = p_new;
        Q.T = T_new;
        try { gmodel.update_thermo_from_pT(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration at call A failed in update_thermo_state_hs.";
            debug {
                msg ~= format("\nIteration %d", count);
                msg ~= format("Exception message from update_thermo_from_pT() was:\n\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }
        h_new = gmodel.enthalpy(Q);
        s_new = gmodel.entropy(Q);
        fh_new = h_given - h_new;
        fs_new = s_given - s_new;
        dfh_dp = (fh_new - fh_old) / (p_new - p_old);
        dfs_dp = (fs_new - fs_old) / (p_new - p_old);
        // Perturb other dimension to get derivatives.
        p_new = p_old;
        T_new = T_old * 1.001;
        Q.p = p_new;
        Q.T = T_new;
        try { gmodel.update_thermo_from_pT(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration at call B failed in update_thermo_state_hs.";
            debug {
                msg ~= format("Iteration %d", count);
                msg ~= format("Exception message from update_thermo_from_pT() was:\n\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }
        h_new = gmodel.enthalpy(Q);
        s_new = gmodel.entropy(Q);
        fh_new = h_given - h_new;
        fs_new = s_given - s_new;
        dfh_dT = (fh_new - fh_old) / (T_new - T_old);
        dfs_dT = (fs_new - fs_old) / (T_new - T_old);

        det = dfh_dp * dfs_dT - dfs_dp * dfh_dT;

        if( fabs(det) < 1.0e-12 ) {
            string msg = "Nearly zero determinant in update_thermo_state_hs.";
            debug { msg ~= format("\n    det = ", det); }
            throw new GasModelException(msg);
        }
        dp = (-dfs_dT * fh_old + dfh_dT * fs_old) / det;
        dT = (dfs_dp * fh_old - dfh_dp * fs_old) / det;
        if( fabs(dp) > MAX_RELATIVE_STEP * p_old ) {
            // move a little toward the goal
            dp_sign = (dp > 0.0 ? 1.0 : -1.0);
            dp = dp_sign * MAX_RELATIVE_STEP * p_old;
        }
        if( fabs(dT) > MAX_RELATIVE_STEP * T_old ) {
            // move a little toward the goal
            dT_sign = (dT > 0.0 ? 1.0 : -1.0);
            dT = dT_sign * MAX_RELATIVE_STEP * T_old;
        }
        p_old += dp;
        T_old += dT;
        // Make sure of consistent thermo state.
        Q.p = p_old;
        Q.T = T_old;
        try { gmodel.update_thermo_from_pT(Q); }
        catch (Exception caughtException) {
            string msg = "Iteration failed at call C in update_thermo_state_hs.";
            debug {
                msg ~= format("\nIteration %d", count);
                msg ~= format("Exception message from update_thermo_from_pT() was:\n\n");
                msg ~= to!string(caughtException);
            }
            throw new GasModelException(msg);
        }
        h_new = gmodel.enthalpy(Q);
        s_new = gmodel.entropy(Q);
        // Prepare for next iteration.
        fh_old = h_given - h_new;
        fs_old = s_given - s_new;
        converged = (fabs(fh_old) < fh_tol) && (fabs(fs_old) < fs_tol);
        ++count;
    } // end while

    if ( count >= MAX_STEPS ) {
        string msg = "Iterations did not converge in update_thermo_state_hs.";
        debug {
            msg ~= format("\n    fh_old = %g, fs_old = %g\n", fh_old, fs_old);
            msg ~= format("    h_given = %.10s, h_given, %.5s\n", h_given, s_given);
            msg ~= "  Supplied Q:" ~ Q.toString();
        }
        throw new GasModelException(msg);
    }

    if( (fabs(fh_old) > fh_tol_fail) || (fabs(fs_old) > fs_tol_fail) ) {
        string msg = "Iterations failed badly in update_thermo_state_hs.";
        debug {
            msg ~= format("\n    h_given = %.10s, h_given, %.5s\n", h_given, s_given);
            msg ~= "  Supplied Q:" ~ Q.toString();
        }
        throw new GasModelException(msg);
    }
} // end update_thermo_state_hs()

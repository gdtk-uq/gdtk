/**
 * fluxcalc.d
 * Convective-Flux calculators, for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-07-23: initial cut, to explore options.
 */

module fluxcalc;

import std.math;
import std.stdio;
import std.conv;
import geom;
import gas;
import flowstate;
import conservedquantities;
import fvcore;
import fvinterface;
import globalconfig;


void compute_interface_flux(ref FlowState Lft, ref FlowState Rght, ref FVInterface IFace,
			    GasModel gmodel, double omegaz=0.0)
/** \brief Compute the inviscid fluxes (in 2D) across the cell interfaces.
 *
 * This is the top-level function that calls the previously selected
 * flux calculator.
 * Much of the detailed work is delegated.
 *
 * \param Lft    : reference to the LEFT flow state
 * \param Rght   : reference to the RIGHT flow state
 * \param IFace  : reference to the interface where the fluxes are to be stored
 *
 * Note that the FlowState objects are tampered with, but should be put back right
 * by the end of the function. [TODO] check this assertion.
 */
{
    // Transform to interface frame of reference.
    Lft.vel.refx -= IFace.gvel.x; Lft.vel.refy -= IFace.gvel.y; Lft.vel.refz -= IFace.gvel.z;
    Rght.vel.refx -= IFace.gvel.x; Rght.vel.refy -= IFace.gvel.y; Rght.vel.refz -= IFace.gvel.z;
    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    Lft.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    Rght.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    // Also transform the magnetic field
    if ( GlobalConfig.MHD ) {
	Lft.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
        Rght.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    }
    // Compute the fluxes in the local frame of the interface.
    final switch (GlobalConfig.flux_calculator) {
    case FluxCalculator.efm:
        efmflx(Lft, Rght, IFace, gmodel);
	break;
    case FluxCalculator.ausmdv:
        ausmdv(Lft, Rght, IFace);
	break;
    case FluxCalculator.adaptive:
        adaptive_flux(Lft, Rght, IFace, gmodel);
	break;
    case FluxCalculator.ausm_plus_up:
        ausm_plus_up(Lft, Rght, IFace);
	break;
    case FluxCalculator.hlle:
        hlle(Lft, Rght, IFace, gmodel);
	break;
    } // end switch
    ConservedQuantities F = IFace.F;
    if (omegaz != 0.0) {
	// Rotating frame.
	double x = IFace.pos.x;
	double y = IFace.pos.y;
	double rsq = x*x + y*y;
	// The conserved quantity is rothalpy,
	// so we need to take -(u**2)/2 off the total energy flux.
	// Note that rotating frame velocity u = omegaz * r.
	F.total_energy -= F.mass * 0.5*omegaz*omegaz*rsq;
    }
    // Transform fluxes back from interface frame of reference to local frame of reference.
    // Flux of Total Energy
    double v_sqr = IFace.gvel.x*IFace.gvel.x + IFace.gvel.y*IFace.gvel.y + IFace.gvel.z*IFace.gvel.z; 
    F.total_energy += 0.5 * F.mass * v_sqr + F.momentum.dot(IFace.gvel);
    // Flux of momentum
    F.momentum += F.mass * IFace.gvel;
    // 
    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    // Also, transform the interface (grid) velocity
    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    // and transform the magnetic field
    if (GlobalConfig.MHD) {
	F.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    }
    return;
} // end compute_interface_flux()

@nogc
void set_flux_vector_in_local_frame(ref ConservedQuantities F, ref FlowState fs)
{
    double rho = fs.gas.rho;
    double un = fs.vel.x;
    double vt1 = fs.vel.y;
    double vt2 = fs.vel.z;
    double p = fs.gas.p;
    double e = 0.0; foreach(elem; fs.gas.e) e += elem;
    double ke = 0.5 * (un*un + vt1*vt1 + vt2*vt2); // Kinetic energy per unit volume.
    
    // Mass flux (mass / unit time / unit area)
    F.mass = rho * un; // The mass flux is relative to the moving interface.
    // Flux of normal momentum
    F.momentum.refx = F.mass * un + p;
    // Flux of tangential momentum
    F.momentum.refy = F.mass * vt1;
    F.momentum.refz = F.mass * vt2;
    // Flux of Total Energy
    F.total_energy = F.mass * (e + ke) + p * un;
    F.tke = F.mass * fs.tke;  // turbulence kinetic energy
    F.omega = F.mass * fs.omega;  // pseudo vorticity
    // Species mass flux
    for ( size_t isp = 0; isp < F.massf.length; ++isp ) {
	F.massf[isp] = F.mass * fs.gas.massf[isp];
    }
    // Individual energies.
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    for ( size_t imode = 1; imode < F.energies.length; ++imode ) {
	F.energies[imode] = F.mass * fs.gas.e[imode];
    }
} // end set_flux_vector_in_local_frame()

@nogc
void set_flux_vector_in_global_frame(ref FVInterface IFace, ref FlowState fs, 
				     double omegaz=0.0)
{
    ConservedQuantities F = IFace.F;
    // Record velocity to restore fs at end.
    double vx = fs.vel.x; double vy = fs.vel.y; double vz = fs.vel.z; 
    // Transform to interface frame of reference.
    // Beware: fs.vel is changed here and restored below.
    fs.vel.refx -= IFace.gvel.x; fs.vel.refy -= IFace.gvel.y; fs.vel.refz -= IFace.gvel.z; 
    IFace.gvel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    fs.vel.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    // also transform the magnetic field
    if ( GlobalConfig.MHD ) {
	fs.B.transform_to_local_frame(IFace.n, IFace.t1, IFace.t2);
    }
    set_flux_vector_in_local_frame(IFace.F, fs);
    if ( omegaz != 0.0 ) {
	// Rotating frame.
	double x = IFace.pos.x;
	double y = IFace.pos.y;
	double rsq = x*x + y*y;
	// The conserved quantity is rothalpy,
	// so we need to take -(u**2)/2 off the total energy flux.
	// Note that rotating frame velocity u = omegaz * r.
	F.total_energy -= F.mass * 0.5*omegaz*omegaz*rsq;
    }

    // Transform fluxes back from interface frame of reference to local frame of reference.
    /* Flux of Total Energy */
    double v_sqr = IFace.gvel.x*IFace.gvel.x + IFace.gvel.y*IFace.gvel.y + IFace.gvel.z*IFace.gvel.z;
    F.total_energy += 0.5 * F.mass * v_sqr + F.momentum.dot(IFace.gvel);
    /* Flux of momentum */
    F.momentum += F.mass * IFace.gvel;

    // Rotate momentum fluxes back to the global frame of reference.
    F.momentum.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    // also transform the interface (grid) velocity
    IFace.gvel.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);	  
    // also transform the magnetic field
    if ( GlobalConfig.MHD ) {
	F.B.transform_to_global_frame(IFace.n, IFace.t1, IFace.t2);
    }
    fs.vel.refx = vx; fs.vel.refy = vy; fs.vel.refz = vz; // restore fs.vel
} // end set_flux_vector_in_global_frame()

@nogc
void ausmdv(in FlowState Lft, in FlowState Rght, ref FVInterface IFace)
// Wada and Liou's flux calculator.
// 
// Implemented from details in their AIAA paper 
// with hints from Ian Johnston.
// Y. Wada and M. -S. Liou (1994)
// A flux splitting scheme with high-resolution and robustness for discontinuities.
// AIAA-94-0083.
{
    // Unpack the flow-state vectors for either side of the interface.
    // Store in work vectors, those quantities that will be neede later.
    double rL = Lft.gas.rho;
    double pL = Lft.gas.p;
    double pLrL = pL / rL;
    double uL = Lft.vel.x;
    double vL = Lft.vel.y;
    double wL = Lft.vel.z;
    double eL = 0.0; foreach(elem; Lft.gas.e) eL += elem;
    double aL = Lft.gas.a;
    double keL = 0.5 * (uL * uL + vL * vL + wL * wL);
    double HL = eL + pLrL + keL;

    double rR = Rght.gas.rho;
    double pR = Rght.gas.p;
    double pRrR = pR / rR;
    double uR = Rght.vel.x;
    double vR = Rght.vel.y;
    double wR = Rght.vel.z;
    double eR = 0.0; foreach(elem; Rght.gas.e) eR += elem;
    double aR = Rght.gas.a;
    double keR = 0.5 * (uR * uR + vR * vR + wR * wR);
    double HR = eR + pRrR + keR;

    /*
     * This is the main part of the flux calculator.
     */
    /*
     * Weighting parameters (eqn 32) for velocity splitting.
     */
    double alphaL = 2.0 * pLrL / (pLrL + pRrR);
    double alphaR = 2.0 * pRrR / (pLrL + pRrR);
    /*
     * Common sound speed (eqn 33) and Mach numbers.
     */
    double am = fmax(aL, aR);
    double ML = uL / am;
    double MR = uR / am;
    /*
     * Left state: 
     * pressure splitting (eqn 34) 
     * and velocity splitting (eqn 30)
     */
    double pLplus, uLplus;
    double duL = 0.5 * (uL + fabs(uL));
    if (fabs(ML) <= 1.0) {
	pLplus = pL * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) * 0.25;
	uLplus = alphaL * ((uL + am) * (uL + am) / (4.0 * am) - duL) + duL;
    } else {
	pLplus = pL * duL / uL;
	uLplus = duL;
    }
    /*
     * Right state: 
     * pressure splitting (eqn 34) 
     * and velocity splitting (eqn 31)
     */
    double pRminus, uRminus;
    double duR = 0.5 * (uR - fabs(uR));
    if (fabs(MR) <= 1.0) {
	pRminus = pR * (MR - 1.0) * (MR - 1.0) * (2.0 + MR) * 0.25;
	uRminus =
	    alphaR * (-(uR - am) * (uR - am) / (4.0 * am) - duR) +
	    duR;
    } else {
	pRminus = pR * duR / uR;
	uRminus = duR;
    }
    /*
     * Mass Flux (eqn 29)
     */
    // The mass flux is relative to the moving interface.
    double ru_half = uLplus * rL + uRminus * rR;
    /*
     * Pressure flux (eqn 34)
     */
    double p_half = pLplus + pRminus;
    /*
     * Momentum flux: normal direction
     *
     * Compute blending parameter s (eqn 37),
     * the momentum flux for AUSMV (eqn 21) and AUSMD (eqn 21)
     * and blend (eqn 36).
     */
    double dp = pL - pR;
    const double K_SWITCH = 10.0;
    dp = K_SWITCH * fabs(dp) / fmin(pL, pR);
    double s = 0.5 * fmin(1.0, dp);

    double ru2_AUSMV = uLplus * rL * uL + uRminus * rR * uR;
    double ru2_AUSMD = 0.5 * (ru_half * (uL + uR) - fabs(ru_half) * (uR - uL));

    double ru2_half = (0.5 + s) * ru2_AUSMV + (0.5 - s) * ru2_AUSMD;

    /*
     * Assemble components of the flux vector.
     */
    ConservedQuantities F = IFace.F;
    F.mass = ru_half;
    F.momentum.refx = ru2_half + p_half;
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	F.momentum.refy = ru_half * vL;
	F.momentum.refz = ru_half * wL;
	F.total_energy = ru_half * HL;
	F.tke = ru_half * Lft.tke;
	F.omega = ru_half * Lft.omega;
    } else {
	/* Wind is blowing from the right */
	F.momentum.refy = ru_half * vR;
	F.momentum.refz = ru_half * wR;
	F.total_energy = ru_half * HR;
	F.tke = ru_half * Rght.tke;
	F.omega = ru_half * Rght.omega;
    }
    /* 
     * Species mass fluxes 
     */
    size_t nsp = F.massf.length;
    for (size_t isp = 0; isp < nsp; ++isp) {
	if (ru_half >= 0.0) {
	    /* Wind is blowing from the left */
	    F.massf[isp] = ru_half * Lft.gas.massf[isp];
	} else {
	    /* Wind is blowing from the right */
	    F.massf[isp] = ru_half * Rght.gas.massf[isp];
	}
    } /* isp loop */
    /* 
     * Individual energies 
     */
    size_t nmodes = F.energies.length;
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	for (size_t imode = 1; imode < nmodes; ++imode) {
	    F.energies[imode] = ru_half * Lft.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += ru_half * Lft.gas.p_e / Lft.gas.rho;
    } else {
	/* Wind is blowing from the right */
	for (size_t imode = 1; imode < nmodes; ++imode) {
	    F.energies[imode] = ru_half * Rght.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += ru_half * Rght.gas.p_e / Rght.gas.rho;
    }

    /*
     * Apply entropy fix (section 3.5 in Wada and Liou's paper)
     */
    const double C_EFIX = 0.125;
    bool caseA = ((uL - aL) < 0.0) && ((uR - aR) > 0.0);
    bool caseB = ((uL + aL) < 0.0) && ((uR + aR) > 0.0);

    double d_ua = 0.0;
    if (caseA && !caseB)
	d_ua = C_EFIX * ((uR - aR) - (uL - aL));
    if (caseB && !caseA)
	d_ua = C_EFIX * ((uR + aR) - (uL + aL));

    if (d_ua != 0.0) {
	F.mass -= d_ua * (rR - rL);
	F.momentum.refx -= d_ua * (rR * uR - rL * uL);
	F.momentum.refy -= d_ua * (rR * vR - rL * vL);
	F.momentum.refz -= d_ua * (rR * wR - rL * wL);
	F.total_energy -= d_ua * (rR * HR - rL * HL);
	F.tke -= d_ua * (rR * Rght.tke - rL * Lft.tke);
	F.omega -= d_ua * (rR * Rght.omega - rL * Lft.omega);
    }   /* end of entropy fix (d_ua != 0) */

    for ( int isp = 0; isp < nsp; ++isp ) {
	if (d_ua != 0.0) F.massf[isp] -= d_ua * (rR * Rght.gas.massf[isp] - 
						 rL * Lft.gas.massf[isp]);
    }

    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    for ( int imode = 1; imode < nmodes; ++imode ) {
	if (d_ua != 0.0)
	    F.energies[imode] -= d_ua * (rR * Rght.gas.e[imode] - rL * Lft.gas.e[imode]);
    }
} // end ausmdv()


void efmflx(in FlowState Lft, in FlowState Rght, ref FVInterface IFace, GasModel gmodel)
/** \brief Compute the fluxes across an interface using
 * the Equilibrium Flux	Method of Macrossan & Pullin
 *
 * \param Lft    : IN     : array of Left flow states
 *     (with velocities in local frame of reference)
 * \param Rght   : IN     : array of Right flow state
 * \param IF     : IN/OUT : array of interface flux data structures
 *
 * \verbatim
 * interface data -- contains...
 *     flux of mass across the interface (kg/s/m**2)
 *     flux of normal momentum
 *     flux of tangential momentum
 *     flux of energy
 *     array of species fluxes 
 *     vibrational energies
 *     free-electron energy
 * \endverbatim 
 */
{
    /*
     * Local variable names reflect the names used in the original
     * FORTRAN code by MNM.
     */
    const double PHI = 1.0;
    double vnL, vpL, vnR, vpR, vqL, vqR;
    double rtL, cmpL, rtR, cmpR;
    double hvsqL, hvsqR;
    double wL, wR, dL, dR;
    double rhoL, rhoR, presL, presR, tR, tL;
    double eL, eR, hL, hR;
    double snL, snR, exL, exR, efL, efR;
    double fmsL, fmsR;
    double cv, cp, con, gam, Rgas;
    double cvL, cvR, RgasL, RgasR;
    double rLsqrt, rRsqrt, alpha;
    int statusf;

    /*
     * Calculate Constants
     */
    /* dtwspi = 1.0 / (2.0 * sqrt ( 3.14159265359 )); */
    const double dtwspi = 0.282094792;

    /*
     * Unpack Left flow state
     */
    rhoL = Lft.gas.rho;
    presL = Lft.gas.p;
    eL = 0.0; foreach(elem; Lft.gas.e) eL += elem;
    hL = eL + presL/rhoL;
    tL = Lft.gas.T[0];
    vnL = Lft.vel.x;
    vpL = Lft.vel.y;
    vqL = Lft.vel.z;

    /*
     * Unpack Right flow state
     */
    rhoR = Rght.gas.rho;
    presR = Rght.gas.p;
    eR = 0.0; foreach(elem; Rght.gas.e) eR += elem;
    hR = eR + presR/rhoR;
    tR = Rght.gas.T[0];
    vnR = Rght.vel.x;
    vpR = Rght.vel.y;
    vqR = Rght.vel.z;

    /* Derive the gas "constants" from the local conditions. */
    cvL = gmodel.Cv(Lft.gas);
    RgasL = presL / (rhoL * tL);
    cvR = gmodel.Cv(Rght.gas);
    RgasR = presR / (rhoR * tR);

    rLsqrt = sqrt(rhoL);
    rRsqrt = sqrt(rhoR);
    alpha = rLsqrt / (rLsqrt + rRsqrt);

    cv = alpha * cvL + (1.0 - alpha) * cvR;
    Rgas = alpha * RgasL + (1.0 - alpha) * RgasR;

    cp = cv + Rgas;
    gam = cp / cv;

    /*
     * Start EFM calculation proper.
     */
    con = 0.5 * (gam + 1.0) / (gam - 1.0);

    rtL = Rgas * tL;
    cmpL = sqrt(2.0 * rtL);
    hvsqL = 0.5 * (vnL * vnL + vpL * vpL + vqL * vqL);

    snL = vnL / (PHI * cmpL);
    exxef(snL, exL, efL);

    wL = 0.5 * (1.0 + efL);
    dL = exL * dtwspi;

    rtR = presR / rhoR;
    cmpR = sqrt(2.0 * rtR);
    hvsqR = 0.5 * (vnR * vnR + vpR * vpR + vqR * vqR);

    snR = vnR / (PHI * cmpR);
    exxef(snR, exR, efR);

    wR = 0.5 * (1.0 - efR);
    dR = -exR * dtwspi;

    /*
     * combine the fluxes
     */
    fmsL = (wL * rhoL * vnL) + (dL * cmpL * rhoL);
    fmsR = (wR * rhoR * vnR) + (dR * cmpR * rhoR);

    ConservedQuantities F = IFace.F;
    F.mass = fmsL + fmsR;

    F.momentum.refx = fmsL * vnL + fmsR * vnR + wL * presL + wR * presR;
    F.momentum.refy = fmsL * vpL + fmsR * vpR;
    F.momentum.refz = fmsL * vqL + fmsR * vqR;

    F.total_energy = (wL * rhoL * vnL) * (hvsqL + hL) +
	(wR * rhoR * vnR) * (hvsqR + hR) +
	(dL * cmpL * rhoL) * (hvsqL + con * rtL) +
	(dR * cmpR * rhoR) * (hvsqR + con * rtR);

    if (F.mass > 0.0) {
	F.tke = F.mass * Lft.tke;
	F.omega = F.mass * Lft.omega;
    } else {
	F.tke = F.mass * Rght.tke;
	F.omega = F.mass * Rght.omega;
    }

    /* 
     * Species mass flux.
     * Presently, this is implemented by assuming that
     * the wind is blowing one way or the other and then
     * picking the appropriate side for the species fractions.
     * Such an approach may not be fully compatible with the
     * EFM approach where there can be fluxes from both sides.
     */
    for (size_t isp = 0; isp < F.massf.length; ++isp ) {
	if (F.mass > 0.0)
	    F.massf[isp] = (F.mass) * Lft.gas.massf[isp];
	else
	    F.massf[isp] = (F.mass) * Rght.gas.massf[isp];
    }   /* isp loop */

    // Individual energies.
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    if (F.mass > 0.0) {
	for (size_t imode = 1; imode < F.energies.length; ++imode) {
	    F.energies[imode] = (F.mass) * Lft.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	// F.energies[$-1] += (F.mass) * Lft.gas.p_e / Lft.gas.rho; [TODO]
    } else {
	for (size_t imode = 1; imode < F.energies.length; ++imode) {
	    F.energies[imode] = F.mass * Rght.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	// F.energies[nmodes-1] += F.mass * Rght.gas.p_e / Rght.gas.rho; [TODO]
    }
} // end efmflx()

@nogc
void exxef(double sn, ref double exx, ref double ef)
/** \brief Compute exp(-x**2) and erf(x) with a polynomial approximation.
 *
 * \param sn   : IN  : x
 * \param &exx : OUT : exp(x**2)
 * \param &ef  : OUT : erf(x)  error function
 */
{
    double snsq, ef1, y;

    const double P = 0.327591100;
    const double A1 = 0.254829592;
    const double A2 = -0.284496736;
    const double A3 = 1.421413741;
    const double A4 = -1.453152027;
    const double A5 = 1.061405429;
    const double LIMIT = 5.0;
    const double EXLIM = 0.138879e-10;
    const double EFLIM = 1.0;

    //#   define DSIGN(val,sgn) ( (sgn >= 0.0)? fabs(val): -fabs(val) )

    if (fabs(sn) > LIMIT) {
        exx = EXLIM;
        ef1 = EFLIM;
    } else {
        snsq = sn * sn;
        exx = exp(-snsq);
        y = 1.0 / (1.0 + P * fabs(sn));
        ef1 = 1.0 - y * (A1 + y * (A2 + y * (A3 + y * (A4 + A5 * y)))) * exx;
    }
    ef = copysign(ef1, sn);
} // end exxef


void adaptive_flux(in FlowState Lft, in FlowState Rght, ref FVInterface IFace, GasModel gmodel)
// This adaptive flux calculator uses uses the Equilibrium Flux Method.
// near shocks and AUSMDV away from shocks, however, we really don't want
// EFM to be used across interfaces with strong shear.
// EFM should still be able to do it's work as we really needed it for the
// situations where the shock is closely aligned with the grid.
// In that situation, we don't expect a stong shear at the interface.
//
// The actual work is passed off to the original flux calculation functions.
{
    double sound_speed = 0.5 * (Lft.gas.a + Rght.gas.a);
    double shear_y = fabs(Lft.vel.y - Rght.vel.y) / sound_speed;
    double shear_z = fabs(Lft.vel.z - Rght.vel.z) / sound_speed;
    bool shear_is_small = fmax(shear_y, shear_z) <= GlobalConfig.shear_tolerance;
    if ( (Lft.S == 1 || Rght.S == 1) && shear_is_small ) {
	efmflx(Lft, Rght, IFace, gmodel);
    } else {
	ausmdv(Lft, Rght, IFace);
    }
} // end adaptive_flux()

@nogc
void ausm_plus_up(in FlowState Lft, in FlowState Rght, ref FVInterface IFace)
// Liou's 2006 AUSM+up flux calculator
//
// A new version of the AUSM-family schemes, based 
// on the low Mach number asymptotic analysis.
// Ironically, this flux calculator causes simulations
// initialised with 0.0 m/s velocities to crash.
//
// RJG -- 26-Apr-2013
// Added a (+ EPSILON) to help with any divide by zero problems.
// That being said, I'm not sure this helps with the 
// crashes at zero velocity because it would seem that the flow
// of control would pass through a different branch for these cases.
//
// M. -S. Liou (2006)
// A sequel to AUSM, Part II: AUSM+-up for all speeds
// Journal of Computational Physics, Vol 214, pp 137-170
//
// This code: W. Y. K. Chan & P. A. Jacobs
{
    // Some helper functions
    @nogc double M1plus(double M) { return 0.5*(M + fabs(M)); }
    @nogc double M1minus(double M) { return 0.5*(M - fabs(M)); }
    @nogc double M2plus(double M) { return 0.25*(M + 1.0)*(M + 1.0); }
    @nogc double M2minus(double M) { return -0.25*(M - 1.0)*(M - 1.0); } 
    @nogc double M4plus(double M, double beta) {
	if ( fabs(M) >= 1.0 ) {
	    return M1plus(M);
	} else {
	    double M2p = M2plus(M);
	    double M2m = M2minus(M);
	    return M2p*(1.0 - 16.0*beta*M2m);
	}
    }
    @nogc double M4minus(double M, double beta) {
	if ( fabs(M) >= 1.0 ) {
	    return M1minus(M);
	} else {
	    double M2p = M2plus(M);
	    double M2m = M2minus(M);
	    return M2m*(1.0 + 16.0*beta*M2p);
	}
    }
    @nogc double P5plus(double M, double alpha) {
	if ( fabs(M) >= 1.0 ) {
	    return (1.0/M)*M1plus(M);
	} else {
	    double M2p = M2plus(M);
	    double M2m = M2minus(M);
	    return M2p*((2.0 - M) - 16.0*alpha*M*M2m);
	}
    }
    @nogc double P5minus(double M, double alpha) {
	if ( fabs(M) >= 1.0 ) {
	    return (1.0/M)*M1minus(M);
	} else {
	    double M2p = M2plus(M);
	    double M2m = M2minus(M);
	    return M2m*((-2.0 - M) + 16.0*alpha*M*M2p);
	}
    }
    // Unpack the flow-state vectors for either side of the interface.
    // Store in work vectors, those quantities that will be neede later.
    double rL = Lft.gas.rho;
    double pL = Lft.gas.p;
    double uL = Lft.vel.x;
    double vL = Lft.vel.y;
    double wL = Lft.vel.z;
    double eL = 0.0; foreach(elem; Lft.gas.e) eL += elem;
    double aL = Lft.gas.a;
    double keL = 0.5 * (uL * uL + vL * vL + wL * wL);
    double HL = eL + pL/rL + keL;

    double rR = Rght.gas.rho;
    double pR = Rght.gas.p;
    double uR = Rght.vel.x;
    double vR = Rght.vel.y;
    double wR = Rght.vel.z;
    double eR = 0.0; foreach(elem; Rght.gas.e) eR += elem;
    double aR = Rght.gas.a;
    double keR = 0.5 * (uR * uR + vR * vR + wR * wR);
    double HR = eR + pR/rR + keR;

    /*
     * This is the main part of the flux calculator.
     */
    /*
     * Interface sound speed (eqns 28 & 30). 
     * An approximation is used instead of these equations as
     * suggested by Liou in his paper (see line below eqn 69).
     */
    double a_half = 0.5 * (aR + aL);
    /*
     * Left and right state Mach numbers (eqn 69).
     */
    double ML = uL / a_half;
    double MR = uR / a_half;
    /*
     * Mean local Mach number (eqn 70).
     */
    double MbarSq = (uL*uL + uR*uR) / (2.0 * a_half *a_half);
    /*
     * Reference Mach number (eqn 71).
     */
    double M0Sq = fmin(1.0, fmax(MbarSq, GlobalConfig.M_inf));
    /*
     * Some additional parameters.
     */
    double fa = sqrt(M0Sq) * (2.0 - sqrt(M0Sq));   // eqn 72
    double alpha = 0.1875 * (-4.0 + 5 * fa * fa);  // eqn 76
    double beta = 0.125;                           // eqn 76

    /*
     * Left state: 
     * M4plus(ML)
     * P5plus(ML)
     */
    double M4plus_ML = M4plus(ML, beta);
    double P5plus_ML = P5plus(ML, alpha);
    
    /*
     * Right state: 
     * M4minus(MR)
     * P5minus(MR)
     */
    double M4minus_MR = M4minus(MR, beta);
    double P5minus_MR = P5minus(MR, alpha);

    /*
     * Pressure diffusion modification for 
     * mass flux (eqn 73) and pressure flux (eqn 75).
     */
    const double KP = 0.25;
    const double KU = 0.75;
    const double SIGMA = 1.0;
    double r_half = 0.5*(rL + rR);
    double Mp = -KP / fa * fmax((1.0 - SIGMA * MbarSq), 0.0) * (pR - pL) / (r_half*a_half*a_half);
    double Pu = -KU * P5plus_ML * P5minus_MR * (rL + rR) * fa * a_half * (uR - uL);
    /*
     * Mass Flux (eqns 73 & 74).
     */
    double M_half = M4plus_ML + M4minus_MR + Mp;
    double ru_half = a_half * M_half;
    if ( M_half > 0.0 ) {
       ru_half *= rL;
    } else {
       ru_half *= rR;
    }
    /*
     * Pressure flux (eqn 75).
     */
    double p_half = P5plus_ML*pL + P5minus_MR*pR + Pu;

    /*
     * Momentum flux: normal direction
     */
    double ru2_half;
    if ( ru_half >= 0.0 ) {
	ru2_half = ru_half * uL;
    } else {
	ru2_half = ru_half * uR;
    }
    /*
     * Assemble components of the flux vector.
     */
    ConservedQuantities F = IFace.F;
    F.mass = ru_half;
    F.momentum.refx = ru2_half + p_half;
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	F.momentum.refy = ru_half * vL;
	F.momentum.refz = ru_half * wL;
	F.total_energy = ru_half * HL;
	F.tke = ru_half * Lft.tke;
	F.omega = ru_half * Lft.omega;
    } else {
	/* Wind is blowing from the right */
	F.momentum.refy = ru_half * vR;
	F.momentum.refz = ru_half * wR;
	F.total_energy = ru_half * HR;
	F.tke = ru_half * Rght.tke;
	F.omega = ru_half * Rght.omega;
    }
    
    /* 
     * Species mass fluxes 
     */
    size_t nsp = F.massf.length;
    for (size_t isp = 0; isp < nsp; ++isp) {
	if (ru_half >= 0.0) {
	    /* Wind is blowing from the left */
	    F.massf[isp] = ru_half * Lft.gas.massf[isp];
	} else {
	    /* Wind is blowing from the right */
	    F.massf[isp] = ru_half * Rght.gas.massf[isp];
	}
    } /* isp loop */
    /* 
     * Individual energies 
     */
    size_t nmodes = F.energies.length;
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	for (size_t imode = 1; imode < nmodes; ++imode) {
	    F.energies[imode] = ru_half * Lft.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += ru_half * Lft.gas.p_e / Lft.gas.rho;
    } else {
	/* Wind is blowing from the right */
	for (size_t imode = 1; imode < nmodes; ++imode) {
	    F.energies[imode] = ru_half * Rght.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	F.energies[nmodes-1] += ru_half * Rght.gas.p_e / Rght.gas.rho;
    }
} // end ausm_plus_up()


void hlle(in FlowState Lft, in FlowState Rght, ref FVInterface IFace, GasModel gmodel)
// HLLE fluxes for MHD.
// From V. Wheatley Matlab implementation
// Author D. M. Bond
// Port to D by PJ, 2014-07-24
{
    @nogc double SAFESQRT(double x) { return (fabs(x)>1.0e-14) ? sqrt(x) : 0.0; }
    // Unpack the flow-state vectors for either side of the interface.
    // Store in work vectors, those quantities that will be neede later.
    double rL = Lft.gas.rho;
    double pL = Lft.gas.p;
    double uL = Lft.vel.x;
    double vL = Lft.vel.y;
    double wL = Lft.vel.z;
    double BxL = Lft.B.x;
    double ByL = Lft.B.y;
    double BzL = Lft.B.z;
    double rR = Rght.gas.rho;
    double pR = Rght.gas.p;
    double uR = Rght.vel.x;
    double vR = Rght.vel.y;
    double wR = Rght.vel.z;
    double BxR = Rght.B.x;
    double ByR = Rght.B.y;
    double BzR = Rght.B.z;

    // Derive the gas "constants" from the local conditions.
    double cvL = gmodel.Cv(Lft.gas);
    double RgasL = gmodel.R(Lft.gas);
    double cvR = gmodel.Cv(Rght.gas);
    double RgasR = gmodel.R(Rght.gas);
    double rLsqrt = sqrt(rL);
    double rRsqrt = sqrt(rR);
    double alpha = rLsqrt / (rLsqrt + rRsqrt);
    double cv = alpha * cvL + (1.0 - alpha) * cvR;
    double Rgas = alpha * RgasL + (1.0 - alpha) * RgasR;
    double cp = cv + Rgas;
    double gam = cp / cv;

    // Compute Roe Average State (currently simple average)
    double rho = 0.5*(rL+rR);
    double p   = 0.5*(pL+pR);
    double u   = 0.5*(uL+uR);
    //v   = 0.5*(vL+vR);
    //w   = 0.5*(wL+wR);
    double Bx  = 0.5*(BxL+BxR);
    double By  = 0.5*(ByL+ByR);
    double Bz  = 0.5*(BzL+BzR);

    // Compute Eigenvalues of Roe Matrix
    //u2=u*u;
    //v2=v*v;
    //w2=w*w;
    //uu=u2+v2+w2;
    double a2 = gam*p/rho;
    double Bx2 = Bx*Bx;
    double Bt2 = By*By + Bz*Bz;
    double BB = Bx2 + Bt2;
    double ca2 = Bx2/rho;
    double alf = a2+BB/rho;
    double als = SAFESQRT(alf*alf-4.0*a2*ca2);
    double cf2 = 0.5*(alf+als);
    double cf = sqrt(cf2);
    double wp = u+cf;
    double wm = u-cf;

    // Compute the Jump in Conserved Variables between L and R
    double BxL2 = BxL*BxL;
    double BtL2 = ByL*ByL + BzL*BzL;
    double BBL = BxL2 + BtL2;
    double ptL = pL + 0.5*BBL;
    double uL2 = uL*uL;
    double uuL = uL2 + vL*vL + wL*wL;
    double aL2 = gam*pL/rL;
    double caL2 = BxL2/rL;
    double alfL = aL2+BBL/rL;
    double alsL = SAFESQRT(alfL*alfL-4.0*aL2*caL2);
    double cfL2 = 0.5*(alfL+alsL);
    double cfL = sqrt(cfL2);
    //wpL = uL+cfL;
    double wmL = uL-cfL;

    double BxR2 = BxR*BxR;
    double BtR2 = ByR*ByR + BzR*BzR;
    double BBR = BxR2 + BtR2;
    double ptR = pR + 0.5*BBR;
    double uR2 = uR*uR;
    double uuR = uR2 + vR*vR + wR*wR;
    double aR2 = gam*pR/rR;
    double caR2 = BxR2/rR;
    double alfR = aR2+BBR/rR;
    double alsR = SAFESQRT(alfR*alfR-4.0*aR2*caR2);
    double cfR2 = 0.5*(alfR+alsR);
    double cfR = sqrt(cfR2);
    double wpR = uR+cfR;
    //wmR = uR-cfR;

    double[8] dU;
    dU[0] = rR - rL;
    dU[1] = rR*uR - rL*uL;
    dU[2] = rR*vR - rL*vL;
    dU[3] = rR*wR - rL*wL;
    dU[4] = BxR - BxL;
    dU[5] = ByR - ByL;
    dU[6] = BzR - BzL;
    dU[7] = (pR - pL)/(gam-1.0) + 0.5*(rR*uuR+BBR) - 0.5*(rL*uuL+BBL);

    double bl = fmin(wmL, wm);
    double br = fmax(wpR, wp);
    double blm = fmin(bl, 0.0);
    double brp = fmax(br, 0.0);

    double fmassL = rL*uL;
    double fmassR = rR*uR;

    double fmomxL = rL*uL2 - BxL2 + ptL;
    double fmomxR = rR*uR2 - BxR2 + ptR;

    double fmomyL = rL*uL*vL - BxL*ByL;
    double fmomyR = rR*uR*vR - BxR*ByR;

    double fmomzL = rL*uL*wL - BxL*BzL;
    double fmomzR = rR*uR*wR - BxR*BzR;

    double fBxL = 0.0;
    double fBxR = 0.0;

    double fByL = uL*ByL - vL*BxL;
    double fByR = uR*ByR - vR*BxR;

    double fBzL = uL*BzL - wL*BxL;
    double fBzR = uR*BzR - wR*BxR;

    double fenergyL = (pL/(gam-1.0)+0.5*(rL*uuL+BBL)+ptL)*uL - (uL*BxL+vL*ByL+wL*BzL)*BxL;
    double fenergyR = (pR/(gam-1.0)+0.5*(rR*uuR+BBR)+ptR)*uR - (uR*BxR+vR*ByR+wR*BzR)*BxR;

    double iden = 1.0/(brp - blm);
    double fac1 = brp*blm;

    ConservedQuantities F = IFace.F;

    F.mass = (brp*fmassL   - blm*fmassR   + fac1*dU[0])*iden;

    F.momentum.refx = (brp*fmomxL   - blm*fmomxR   + fac1*dU[1])*iden;
    F.momentum.refy = (brp*fmomyL   - blm*fmomyR   + fac1*dU[2])*iden;
    F.momentum.refz = (brp*fmomzL   - blm*fmomzR   + fac1*dU[3])*iden;

    F.B.refx = (brp*fBxL     - blm*fBxR     + fac1*dU[4])*iden;
    F.B.refy = (brp*fByL     - blm*fByR     + fac1*dU[5])*iden;
    F.B.refz = (brp*fBzL     - blm*fBzR     + fac1*dU[6])*iden;

    F.total_energy = (brp*fenergyL - blm*fenergyR + fac1*dU[7])*iden;

    /*
     * Species mass fluxes
     */
    for (size_t isp = 0; isp < F.massf.length; ++isp ) {
	if (F.mass >= 0.0) {
	    /* Wind is blowing from the left */
	    F.massf[isp] = F.mass * Lft.gas.massf[isp];
	} else {
	    /* Wind is blowing from the right */
	    F.massf[isp] = F.mass * Rght.gas.massf[isp];
	}
    } /* isp loop */
    /*
     * Individual energies
     */
    // NOTE: renergies[0] is never used so skipping (DFP 10/12/09)
    if (F.mass >= 0.0) {
	/* Wind is blowing from the left */
	for (size_t imode = 1; imode < F.energies.length; ++imode ) {
	    F.energies[imode] = F.mass * Lft.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	// F.energies[$-1] += F.mass * Lft.gas.p_e / Lft.gas.rho; [TODO]
    } else {
	/* Wind is blowing from the right */
	for (size_t imode = 1; imode < F.energies.length; ++imode ) {
	    F.energies[imode] = F.mass * Rght.gas.e[imode];
	}
	// NOTE: - the following relies on the free-electron mode being the last mode
	//       - for single temp models F_renergies isn't used
	//       - for multitemp modes with no free-electrons p_e is zero
	// Add electron pressure work term onto final energy mode
	// F.energies[$-1] += F.mass * Rght.gas.p_e / Rght.gas.rho; [TODO]
    }
} // end hlle()

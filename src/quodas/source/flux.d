// flux.d
module flux;
import std.conv;
import std.math;
import finite_volume;
import number;
import complexify;


void compute_flux(FVInterface f, number gamma, string flux_calc) {
        
    if (flux_calc == "van_leer") { van_leer(f, gamma); }
    else if (flux_calc == "ausmdv") { ausmdv(f, gamma); }
    else {
        string msg = text("flux_calc is not valid. Select either 'aumsdv' or 'van_leer'.");
        throw new Error(msg);
    }
    
} // end compute_flux()

void van_leer(FVInterface f, number gamma) {
    /*
      Van Leer flux calculator
      Ref.:
          Riemann Solvers and Numerical Methods for Fluid Dynamics, 
          3rd Edition, Toro, 2009, pg 277.
     
      note: lft  == plus
            rght == minus
    */
    number F_lft, M_lft, a_lft;
    number F_rght, M_rght, a_rght;
    number lft_rho = f.left_fs.rho;
    number rght_rho = f.right_fs.rho;
    number lft_u = f.left_fs.vel;
    number rght_u = f.right_fs.vel;
    number lft_p = f.left_fs.p;
    number rght_p = f.right_fs.p;

    // left state
    a_lft = sqrt( (gamma*lft_p)/lft_rho );
    M_lft = lft_u/a_lft;

    // right state
    a_rght = sqrt( (gamma*rght_p)/rght_rho );
    M_rght = rght_u/a_rght;

    number M = 0.5*(M_lft+M_rght);
    if (M >= 1.0) {
        f.Fmass = lft_rho * a_lft * M_lft;
        f.Fmom = lft_rho * a_lft*a_lft*(M_lft*M_lft+1.0/gamma);
        f.Fenergy = lft_rho * a_lft*a_lft*a_lft * M_lft * (0.5*M_lft*M_lft + 1.0/(gamma-1.0));
    } else if (M <= -1.0) {
	f.Fmass = rght_rho * a_rght * M_rght;
	f.Fmom = rght_rho * a_rght*a_rght*(M_rght*M_rght+1.0/gamma);
	f.Fenergy = rght_rho * a_rght*a_rght*a_rght * M_rght * (0.5*M_rght*M_rght + 1.0/(gamma-1.0));
    } else { 
	F_lft = 0.25*lft_rho*a_lft*(1.0+M_lft)*(1.0+M_lft);
	F_rght = -0.25*rght_rho*a_rght*(1.0-M_rght)*(1.0-M_rght);
	f.Fmass = F_lft+F_rght;

	F_lft = 0.25*lft_rho*a_lft*(1.0+M_lft)*(1.0+M_lft) * ( 2.0*a_lft/gamma * ( (gamma-1.0)/2.0 * M_lft + 1.0) );
	F_rght = -0.25*rght_rho*a_rght*(1.0-M_rght)*(1.0-M_rght) * ( 2.0*a_rght/gamma * ( (gamma-1.0)/2.0 * M_rght - 1.0) );
	f.Fmom = F_lft+F_rght;

	F_lft = 0.25*lft_rho*a_lft*(1.0+M_lft)*(1.0+M_lft) *
	    ((2.0*a_lft*a_lft)/(gamma*gamma - 1.0) * ((gamma-1.0)/2.0 * M_lft + 1.0) * ((gamma-1.0)/2.0 * M_lft + 1.0));
	F_rght = -0.25*rght_rho*a_rght*(1.0-M_rght)*(1.0-M_rght) *
	    ((2.0*a_rght*a_rght)/(gamma*gamma - 1.0) * ((gamma-1.0)/2.0 * M_rght - 1.0)*((gamma-1.0)/2.0 * M_rght - 1.0));
	f.Fenergy = F_lft+F_rght;
    }
} // end van_leer()

void ausmdv(FVInterface f, number gamma) {
    /* 
       Wada and Liou's flux calculator.
       Ref.:
           A flux splitting scheme with high-resolution and robustness for discontinuities,
           Y. Wada and M. -S. Liou, 1994, AIAA-94-0083.
       
       note: implementation taken verbatim from Eilmer4 compressible gas flow code
             with permission from Rowan J. Gollan & Peter A. Jacobs
             http://cfcfd.mechmining.uq.edu.au/eilmer/
    */
    
    number K_SWITCH = 10.0; //10.0;
    number C_EFIX = 0.125; //0.125;
    number rL;
    number rR;
    number pL;
    number pR;
    number uL;
    number uR;
    number aL;
    number aR;
    number HL;
    number HR;
    number pLrL;
    number pRrR;
    number ML;
    number MR;
    number eL;
    number eR;
    number keL;
    number keR;
    number alphaL;
    number alphaR;
    number am;
    number pLplus;
    number pRminus;
    number uLplus;
    number uRminus;
    number duL;
    number duR;
    number p_half;
    number ru_half;
    number ru2_half;
    number dp;
    number s;
    number ru2_AUSMV;
    number ru2_AUSMD;
    int caseA, caseB;
    number d_ua;
    /*
     * Unpack the flow-state vectors for either side of the interface.
     * Store in work vectors, those quantities that will be neede later.
     */
    rL = f.left_fs.rho;
    pL = f.left_fs.p;
    pLrL = pL / rL;
    uL = f.left_fs.vel;
    eL = pL/((gamma-1.0)*rL); 
    aL = sqrt((pL*gamma)/rL);
    keL = 0.5 * (uL * uL);
    HL = eL + pLrL + keL;

    rR = f.right_fs.rho;
    pR = f.right_fs.p;
    pRrR = pR / rR;
    uR = f.right_fs.vel;
    eR = pR/((gamma-1.0)*rR); 
    aR = sqrt((pR*gamma)/rR);
    keR = 0.5 * (uR * uR);
    HR = eR + pRrR + keR;
    /*
     * This is the main part of the flux calculator.
     */
    /*
     * Weighting parameters (eqn 32) for velocity splitting.
     */
    alphaL = 2.0 * pLrL / (pLrL + pRrR);
    alphaR = 2.0 * pRrR / (pLrL + pRrR);
    /*
     * Common sound speed (eqn 33) and Mach numbers.
     */
    am = fmax(aL, aR);
    ML = uL / am;
    MR = uR / am;
    /*
     * Left state: 
     * pressure splitting (eqn 34) 
     * and velocity splitting (eqn 30)
     */
    duL = 0.5 * (uL + fabs(uL));
    if (fabs(ML) <= 1.0) {
	pLplus = pL * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) * 0.25;
	uLplus =
	    alphaL * ((uL + am) * (uL + am) / (4.0 * am) - duL) +
	    duL;
    } else {
	pLplus = pL * duL / uL;
	uLplus = duL;
    }
    /*
     * Right state: 
     * pressure splitting (eqn 34) 
     * and velocity splitting (eqn 31)
     */
    duR = 0.5 * (uR - fabs(uR));
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
    ru_half = uLplus * rL + uRminus * rR;
    /*
     * Pressure flux (eqn 34)
     */
    p_half = pLplus + pRminus;
    /*
     * Momentum flux: normal direction
     *
     * Compute blending parameter s (eqn 37),
     * the momentum flux for AUSMV (eqn 21) and AUSMD (eqn 21)
     * and blend (eqn 36).
     */
    dp = pL - pR;
    dp = K_SWITCH * fabs(dp) / fmin(pL, pR);
    s = 0.5 * fmin(1.0, dp);

    ru2_AUSMV = uLplus * rL * uL + uRminus * rR * uR;
    ru2_AUSMD = 0.5 * (ru_half * (uL + uR) -
		       fabs(ru_half) * (uR - uL));

    ru2_half = (0.5 + s) * ru2_AUSMV + (0.5 - s) * ru2_AUSMD;

    /*
     * Assemble components of the flux vector.
     */
    f.Fmass = ru_half;
    f.Fmom = ru2_half + p_half;
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	f.Fenergy = ru_half * HL;
    } else {
	/* Wind is blowing from the right */
	f.Fenergy = ru_half * HR;
    }
    /*
     * Apply entropy fix (section 3.5 in Wada and Liou's paper)
     */
    caseA = ((uL - aL) < 0.0) && ((uR - aR) > 0.0);
    caseB = ((uL + aL) < 0.0) && ((uR + aR) > 0.0);

    d_ua = 0.0;
    if (caseA && !caseB)
	d_ua = C_EFIX * ((uR - aR) - (uL - aL));
    if (caseB && !caseA)
	d_ua = C_EFIX * ((uR + aR) - (uL + aL));

    if (d_ua != 0.0) {
	f.Fmass -= d_ua * (rR - rL);
	f.Fmom -= d_ua * (rR * uR - rL * uL);
	f.Fenergy -= d_ua * (rR * HR - rL * HL);
    }   /* end of entropy fix (d_ua != 0) */
} // end ausmdv()

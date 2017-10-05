/++ 
    1D (with area change) Euler-Adjoint solver for a supersonic nozzle
    author: Kyle Damm, 13/07/2017
++/
import std.stdio;
import std.math;
import std.string;
import std.file;
import std.conv;
import std.algorithm;
import std.exception;

//-----------------------------------------------------
// Set global simulations parameters
//-----------------------------------------------------

immutable double length = 10.0;            // length of nozzle in metres
immutable size_t ncells = 100;             // number of cells
immutable size_t nghost = 4;              // number of ghost cells
immutable size_t ninterfaces = ncells+1;  // number of interfaces

//-----------------------------------------------------
// Classes and helper functions
//-----------------------------------------------------

struct flowstate
{
    double p;     // cell pressure
    double rho;   // cell density
    double u;     // cell velocity
}

class fvcell {
public:
    size_t id;        // id number
    double xpos;      // position
    double dx;        // width
    double vol;       // volume
    double p;         // pressure
    double rho;       // density
    double u;         // velocity
    double ru;        // momentum
    double rE;        // energy
    
    this(size_t id, double xpos, double dx) {
	this.id = id;
	this.xpos = xpos;
	this.dx = dx;
    }
}

class fvinterface {
public:
    size_t id;           // id number
    double xpos;         // position
    double area;         // area
    double mass;         // mass flux
    double momentum;     // momentum flux
    double energy;       // energy flux
    flowstate left_fs;   // left flowstate
    flowstate right_fs;  // right flowstate
    
    this(size_t id, double xpos) {
	this.id = id;
	this.xpos = xpos;
    }
}

void first_order_interpolation(fvcell L0, fvcell R0, fvinterface f) {
    // First order interpolation -- copies data from cell center to interface left, and right flowstates.

    // left flow state
    f.left_fs.p = L0.p; f.left_fs.rho = L0.rho; f.left_fs.u = L0.u;

    // right flow state
    f.right_fs.p = R0.p; f.right_fs.rho = R0.rho; f.right_fs.u = R0.u;

} // end first_order_interpolation

void second_order_interpolation_with_van_albada_limiting(fvcell L0, fvcell L1, fvcell R0, fvcell R1, fvinterface f) {
    // Second order interpolation with Van albada limiting (as outlined in Johnston's thesis [1999]).

    double delta_minus; double delta_plus; double S; double eps = 1.0e-12; double k = 0.0;

    // left fs state
    // pressure
    delta_minus = (L0.p - L1.p)/(0.5*(L0.dx+L1.dx));
    delta_plus = (R0.p - L0.p)/(0.5*(R0.dx+L0.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.left_fs.p = L0.p + (L0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);

    // density
    delta_minus = (L0.rho - L1.rho)/(0.5*(L0.dx+L1.dx));
    delta_plus = (R0.rho - L0.rho)/(0.5*(R0.dx+L0.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.left_fs.rho = L0.rho+ (L0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);

    // velocity
    delta_minus = (L0.u - L1.u)/(0.5*(L0.dx+L1.dx));
    delta_plus = (R0.u - L0.u)/(0.5*(R0.dx+L0.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.left_fs.u = L0.u + (L0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);

    // right flow state
    // pressure
    delta_minus = (R0.p - L0.p)/(0.5*(R0.dx+L0.dx));
    delta_plus = (R1.p - R0.p)/(0.5*(R0.dx+R1.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.right_fs.p = R0.p + (R0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);

    // density
    delta_minus = (R0.rho - L0.rho)/(0.5*(R0.dx+L0.dx));
    delta_plus = (R1.rho - R0.rho)/(0.5*(R0.dx+R1.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.right_fs.rho = R0.rho + (R0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);

    // velocity
    delta_minus = (R0.u - L0.u)/(0.5*(R0.dx+L0.dx));
    delta_plus = (R1.u - R0.u)/(0.5*(R0.dx+R1.dx));
    S = (2.0*delta_plus*delta_minus+eps)/(delta_plus*delta_plus+delta_minus*delta_minus+eps);
    f.right_fs.u = R0.u + (R0.dx*S)/4.0 * ((1.0-S*k)*delta_minus+(1.0+S*k)*delta_plus);

} // end second_order_interpolation_with_van_albada_limiting

void van_leer_flux_calculator(fvinterface f, double gamma) {
    // note lft == plus, rght == minus

    double F_lft; double M_lft; double a_lft;
    double F_rght; double M_rght; double a_rght;
    double lft_rho = f.left_fs.rho; double rght_rho = f.right_fs.rho;
    double lft_u = f.left_fs.u; double rght_u = f.right_fs.u;
    double lft_p = f.left_fs.p; double rght_p = f.right_fs.p;

    // left state
    a_lft = sqrt( (gamma*lft_p)/lft_rho );
    M_lft = lft_u/a_lft;

    // right state
    a_rght = sqrt( (gamma*rght_p)/rght_rho );
    M_rght = rght_u/a_rght;

    double M = 0.5*(M_lft+M_rght); // average Mach number
        
    if (M >= 1.0) {
	// mass flux
	f.mass = lft_rho * a_lft * M_lft;
	// momentum flux
	f.momentum = lft_rho * a_lft*a_lft*(M_lft*M_lft+1.0/gamma);
	// energy flux
	f.energy = lft_rho * a_lft*a_lft*a_lft * M_lft * (0.5*M_lft*M_lft + 1.0/(gamma-1.0));
    }
    else if (M <= -1.0) {
	// mass flux
	f.mass = rght_rho * a_rght * M_rght;
	// momentum flux
	f.momentum = rght_rho * a_rght*a_rght*(M_rght*M_rght+1.0/gamma);
	// energy flux
	f.energy = rght_rho * a_rght*a_rght*a_rght * M_rght * (0.5*M_rght*M_rght + 1.0/(gamma-1.0));
    }
    else { 
	// mass flux
	F_lft = 0.25*lft_rho*a_lft*pow( (1.0+M_lft), 2.0 );
	F_rght = -0.25*rght_rho*a_rght*pow( (1.0-M_rght), 2.0 );
	f.mass = F_lft+F_rght;
	// momentum flux
	F_lft = 0.25*lft_rho*a_lft*pow( (1.0+M_lft), 2.0 ) * ( 2.0*a_lft/gamma * ( (gamma-1.0)/2.0 * M_lft + 1.0) );
	F_rght = -0.25*rght_rho*a_rght*pow( (1.0-M_rght), 2.0 ) * ( 2.0*a_rght/gamma * ( (gamma-1.0)/2.0 * M_rght - 1.0) );
	f.momentum = F_lft+F_rght;
	// energy flux
	F_lft = 0.25*lft_rho*a_lft*pow( (1.0+M_lft), 2.0 ) *
	    ((2.0*a_lft*a_lft)/(gamma*gamma - 1.0) * pow( (gamma-1.0)/2.0 * M_lft + 1.0, 2.0 ));
	F_rght = -0.25*rght_rho*a_rght*pow( (1.0-M_rght), 2.0 ) *
	    ((2.0*a_rght*a_rght)/(gamma*gamma - 1.0) * pow( (gamma-1.0)/2.0 * M_rght - 1.0, 2.0 ));
	f.energy = F_lft+F_rght;
    }
} // end van_leer_flux_calculator

void ausm_plus_up_flux_calculator(fvinterface f, double M_inf, double gamma)
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
    double M1plus(double M) { return 0.5*(M + fabs(M)); }
    double M1minus(double M) { return 0.5*(M - fabs(M)); }
    double M2plus(double M) { return 0.25*(M + 1.0)*(M + 1.0); }
    double M2minus(double M) { return -0.25*(M - 1.0)*(M - 1.0); } 
    double M4plus(double M, double beta) {
	if ( fabs(M) >= 1.0 ) {
	    return M1plus(M);
	} else {
	    double M2p = M2plus(M);
	    double M2m = M2minus(M);
	    return M2p*(1.0 - 16.0*beta*M2m);
	}
    }
    double M4minus(double M, double beta) {
	if ( fabs(M) >= 1.0 ) {
	    return M1minus(M);
	} else {
	    double M2p = M2plus(M);
	    double M2m = M2minus(M);
	    return M2m*(1.0 + 16.0*beta*M2p);
	}
    }
    double P5plus(double M, double alpha) {
	if ( fabs(M) >= 1.0 ) {
	    return (1.0/M)*M1plus(M);
	} else {
	    double M2p = M2plus(M);
	    double M2m = M2minus(M);
	    return M2p*((2.0 - M) - 16.0*alpha*M*M2m);
	}
    }
    double P5minus(double M, double alpha) {
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
    double rL = f.left_fs.rho;
    double pL = f.left_fs.p;
    double pLrL = pL / rL;
    double uL = f.left_fs.u;
    double eL = pL/((gamma-1.0)*rL);
    double aL = sqrt((pL*gamma)/rL);
    double keL = 0.5 * (uL * uL);
    double HL = eL + pLrL + keL;
    //
    double rR = f.right_fs.rho;
    double pR = f.right_fs.p;
    double pRrR = pR / rR;
    double uR = f.right_fs.u;
    double eR = pR/((gamma-1.0)*rR);
    double aR = sqrt((pR*gamma)/rR);
    double keR = 0.5 * (uR * uR);
    double HR = eR + pRrR + keR;

    //
    // This is the main part of the flux calculator.
    //
    // Interface sound speed (eqns 28 & 30). 
    // An approximation is used instead of these equations as
    // suggested by Liou in his paper (see line below eqn 69).
    double a_half = 0.5 * (aR + aL);
    // Left and right state Mach numbers (eqn 69).
    double ML = uL / a_half;
    double MR = uR / a_half;
    // Mean local Mach number (eqn 70).
    double MbarSq = (uL*uL + uR*uR) / (2.0 * a_half *a_half);
    // Reference Mach number (eqn 71).
    double M0Sq = fmin(1.0, fmax(MbarSq, M_inf));
     // Some additional parameters.
    double fa = sqrt(M0Sq) * (2.0 - sqrt(M0Sq));   // eqn 72
    double alpha = 0.1875 * (-4.0 + 5 * fa * fa);  // eqn 76
    double beta = 0.125;                           // eqn 76
    // Left state: 
    // M4plus(ML)
    // P5plus(ML)
    double M4plus_ML = M4plus(ML, beta);
    double P5plus_ML = P5plus(ML, alpha);
    // Right state: 
    // M4minus(MR)
    // P5minus(MR)
    double M4minus_MR = M4minus(MR, beta);
    double P5minus_MR = P5minus(MR, alpha);
    // Pressure diffusion modification for 
    // mass flux (eqn 73) and pressure flux (eqn 75).
    const double KP = 0.25;
    const double KU = 0.75;
    const double SIGMA = 1.0;
    double r_half = 0.5*(rL + rR);
    double Mp = -KP / fa * fmax((1.0 - SIGMA * MbarSq), 0.0) * (pR - pL) / (r_half*a_half*a_half);
    double Pu = -KU * P5plus_ML * P5minus_MR * (rL + rR) * fa * a_half * (uR - uL);
    // Mass Flux (eqns 73 & 74).
    double M_half = M4plus_ML + M4minus_MR + Mp;
    double ru_half = a_half * M_half;
    if ( M_half > 0.0 ) {
       ru_half *= rL;
    } else {
       ru_half *= rR;
    }
    // Pressure flux (eqn 75).
    double p_half = P5plus_ML*pL + P5minus_MR*pR + Pu;
    // Momentum flux: normal direction
    double ru2_half;
    if (ru_half >= 0.0) {
	ru2_half = ru_half * uL;
    } else {
	ru2_half = ru_half * uR;
    }
    // Assemble components of the flux vector.
    f.mass = ru_half;
    if (ru_half >= 0.0) {
	// Wind is blowing from the left.
	f.momentum = ru2_half+p_half;
	f.energy = ru_half * HL;
    } else {
	// Wind is blowing from the right.
	f.momentum = ru2_half+p_half;
	f.energy = ru_half * HR;
    }
} // end ausm_plus_up()


void ausmdv_flux_calculator(fvinterface f, double gamma) {
    double K_SWITCH = 10.0;
    double C_EFIX = 0.125;
    double rL, rR;
    double pL, pR;
    double uL, uR;
    double aL, aR;
    double HL, HR;
    double pLrL, pRrR;
    double ML, MR;
    double eL, eR;
    double keL, keR;
    double alphaL, alphaR, am;
    double pLplus, pRminus;
    double uLplus, uRminus;
    double duL, duR;

    double p_half, ru_half, ru2_half;
    double dp, s, ru2_AUSMV, ru2_AUSMD;

    int caseA, caseB;
    double d_ua;
    
    /*
     * Unpack the flow-state vectors for either side of the interface.
     * Store in work vectors, those quantities that will be neede later.
     */
    rL = f.left_fs.rho;
    pL = f.left_fs.p;
    pLrL = pL / rL;
    uL = f.left_fs.u;
    eL = pL/((gamma-1.0)*rL); 
    aL = sqrt((pL*gamma)/rL);
    keL = 0.5 * (uL * uL);
    HL = eL + pLrL + keL;

    rR = f.right_fs.rho;
    pR = f.right_fs.p;
    pRrR = pR / rR;
    uR = f.right_fs.u;
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
    f.mass = ru_half;
    f.momentum = ru2_half + p_half;
    if (ru_half >= 0.0) {
	/* Wind is blowing from the left */
	f.energy = ru_half * HL;
    } else {
	/* Wind is blowing from the right */
	f.energy = ru_half * HR;
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
	f.mass -= d_ua * (rR - rL);
	f.momentum -= d_ua * (rR * uR - rL * uL);
	f.energy -= d_ua * (rR * HR - rL * HL);
    }   /* end of entropy fix (d_ua != 0) */
} // end ausmdv_flux_calculator

void encode_conserved_variables(fvcell cell, double gamma) {
    double e = cell.p / (cell.rho*(gamma - 1.0));
    double ke = 0.5*cell.u*cell.u;
    cell.ru = cell.rho * cell.u;
    cell.rE = cell.rho*(e + ke);

} // end encode_conserved_variables

void decode_conserved_variables(fvcell cell, double gamma) {
    cell.u = cell.ru/cell.rho;
    double ke = 0.5*cell.u*cell.u;
    double e = cell.rE/cell.rho - ke;
    cell.p = cell.rho*e*(gamma-1.0);
} // end decode_conserved_variables

void matrixMult(ref double[3*ncells][3*ncells] A, ref double[3*ncells][3*ncells] B, ref double[3*ncells][3*ncells] C) {
    for (int i = 0; i < 3*ncells; i++) {
        for (int j = 0; j < 3*ncells; j++) {
            C[i][j] = 0;
            for (int k = 0; k < 3*ncells; k++) {
                C[i][j] += A[i][k]*B[k][j];
	    }
	}
    }
} // end matrixMult

void linearSolve(ref double[3*ncells+1][3*ncells] c, double very_small_value=1.0e-16) {
    // solves Ax=b using Gauss Jordan eilimination,
    // where c = [A|b] output is [I|x] where x is the solution vector.
    
    void swapRows(ref double[3*ncells+1][3*ncells] c, size_t i1, size_t i2) {
	swap(c[i1], c[i2]);
    }
    
    foreach(j; 0 .. c.length) {
	// Select pivot.
	size_t p = j;
	foreach(i; j+1 .. c.length) {
	    if ( abs(c[i][j]) > abs(c[p][j]) ) p = i;
	}
	if ( abs(c[p][j]) < very_small_value ) {
	    throw new Exception("matrix is essentially singular");
	}
	if ( p != j ) swapRows(c, p,j);
	// Scale row j to get unity on the diagonal.
	double cjj = c[j][j];
	foreach(col; 0 .. c[0].length) c[j][col] /= cjj;
	// Do the elimination to get zeros in all off diagonal values in column j.
	foreach(i; 0 .. c.length) {
	    if ( i == j ) continue;
	    double cij = c[i][j];
	    foreach(col; 0 .. c.length) c[i][col] -= cij * c[j][col]; 
	}
    } // end foreach j
} // end linearSolve

void matrixInv(ref double[3*ncells][3*ncells] matrix, ref double[3*ncells][3*ncells] inverse) {
    // gauss jordan elimination of [A|I] where the output is [I|B] where B is the inverse of A.

    int N = 3*ncells;
    static double[2*3*ncells][3*ncells] c;
    foreach(i; 0 .. N) {
	foreach(j; 0.. N) {
	    c[i][j] = matrix[i][j];
	}
    }
    foreach(i; 0 .. N) {
	foreach(j; N .. 2*N) {
	    if (N+i == j) {
		c[i][j] = 1.0;
	    } else {
		c[i][j] = 0.0;
	    }
	}
    }
    foreach(j; 0 .. N) {
	    // Select pivot.
	    size_t p = j;
	    foreach(i; j+1 .. N) {
		if ( abs(c[i][j]) > abs(c[p][j]) ) p = i;
	    }
	    //if (abs(c[p][j]) <= very_small_value) return -1; // singular
	    if ( p != j ) { // Swap rows
		foreach(col; 0 .. 2*N) {
		    double tmp = c[p][col]; c[p][col] = c[j][col]; c[j][col] = tmp;
		}
	    }
	    // Scale row j to get unity on the diagonal.
	    double cjj = c[j][j];
	    foreach(col; 0 .. 2*N) c[j][col] /= cjj;
	    // Do the elimination to get zeros in all off diagonal values in column j.
	    foreach(i; 0 .. N) {
		if ( i == j ) continue;
		double cij = c[i][j];
		foreach(col; 0 .. 2*N) c[i][col] -= cij * c[j][col]; 
	    }
    } // end foreach j
    foreach(i; 0 .. N) {
	foreach(j; N .. 2*N) {
	    inverse[i][j-N] = c[i][j];
	}
    }
} // end matrixInv

void solve(ref fvcell[ncells+nghost] global_cells,
	   ref fvinterface[ninterfaces] global_interfaces,
	   ref fvinterface[ninterfaces] perturbed_interfaces,
	   double dt, double sim_time, double final_time, size_t step, size_t max_step, string flux_calc, string simulation_type,
	   size_t interpolation_order, size_t outflow_bc_order, ref double[3*ncells][3*ncells] JV,
	   ref double[3*ncells][3*ncells] JU, ref double[3*ncells][3*ncells] transform,
	   ref double[3*ncells] R, ref double[3*ncells] dJdV, ref double[3*ncells] dJdU, ref double[3*ncells][3*ncells] JUT,
	   ref double[3*ncells][3*ncells] invJUT, ref double[3*ncells][3*ncells] JVT,
	   ref double[3*ncells][3*ncells] invJVT, ref double[3*ncells] psi,
	   double eps, ref double J, ref double dLdB, ref double dLdC, ref double dLdD, double dx, double gamma,
	   double b, double c, double d, double yo, double scale, string solver, string adjoint_form) {
    //-----------------------------------------------------
    // Explicit Flow Solver
    //-----------------------------------------------------

    // set initial conditions
    double Mach = 1.5;
    double p_inflow = 101.325e3; // Pa
    double rho_inflow = 1.0;   // kg/m^3
    double a_inflow = sqrt((p_inflow*gamma)/rho_inflow);
    double u_inflow = Mach * a_inflow;  // m/s
    
    foreach(i; 0 .. ncells) {
	global_cells[i+2].p = p_inflow;
	global_cells[i+2].rho = rho_inflow;
	global_cells[i+2].u = 0.6 * a_inflow;
    }
    
    
    // begin by computing conserved quantities at cell centers
    foreach(i; 0 .. ncells+nghost) {
	encode_conserved_variables(global_cells[i], gamma);	
    }
    
    while (sim_time <= final_time && step < max_step) {
	//writeln("==================================");
	//writeln("STEP: ", step, " TIME: ", sim_time);
	//writeln("==================================");
	
	//-----------------------------------------------------
	// Apply Boundary Conditions
	//-----------------------------------------------------
	
	if (simulation_type == "sod") {
	    // Sod's shocktube bcs
	    
	    // left boundary -- wall
	    global_cells[0].rho = global_cells[2].rho;
	    global_cells[1].rho = global_cells[2].rho;
	    global_cells[0].u = -global_cells[2].u;
	    global_cells[1].u = -global_cells[2].u;
	    global_cells[0].p = global_cells[2].p;
	    global_cells[1].p = global_cells[2].p;
	    // right boundary -- wall
	    global_cells[ncells+2].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+3].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+2].u = -global_cells[ncells+1].u;
	    global_cells[ncells+3].u = -global_cells[ncells+1].u;
	    global_cells[ncells+2].p = global_cells[ncells+1].p;
	    global_cells[ncells+3].p = global_cells[ncells+1].p;
	}
	else if (simulation_type == "nozzle") {
	    // Nozzle bcs
	    
	    // left boundary -- inflow
	    global_cells[0].rho = global_cells[0].rho;
	    global_cells[1].rho = global_cells[1].rho;
	    global_cells[0].u = global_cells[0].u;
	    global_cells[1].u = global_cells[1].u;
	    global_cells[0].p = global_cells[0].p;
	    global_cells[1].p = global_cells[1].p;
	    
	    if (outflow_bc_order == 1) {
		// right boundary -- outflow
		
		// first order
		global_cells[ncells+2].rho = global_cells[ncells+1].rho;
		global_cells[ncells+3].rho = global_cells[ncells+1].rho;
		global_cells[ncells+2].u = global_cells[ncells+1].u;
		global_cells[ncells+3].u = global_cells[ncells+1].u;
		global_cells[ncells+2].p = 2.5*global_cells[0].p;
		global_cells[ncells+3].p = 2.5*global_cells[0].p;
		//global_cells[ncells+2].p = global_cells[ncells+1].p;
		//global_cells[ncells+3].p = global_cells[ncells+1].p;
	    }
	    else {
		// second order
		double drhodx; double dudx; double dpdx;
		drhodx = (global_cells[ncells+1].rho - global_cells[ncells].rho)/dx;
		dudx = (global_cells[ncells+1].u - global_cells[ncells].u)/dx;
		dpdx = (global_cells[ncells+1].p - global_cells[ncells].p)/dx;
		
		global_cells[ncells+2].rho = global_cells[ncells+1].rho+drhodx*dx;
		global_cells[ncells+3].rho = global_cells[ncells+1].rho+drhodx*2.0*dx;
		global_cells[ncells+2].u = global_cells[ncells+1].u+dudx*dx;
		global_cells[ncells+3].u = global_cells[ncells+1].u+dudx*2.0*dx;
		global_cells[ncells+2].p = global_cells[ncells+1].p+dpdx*dx;
		global_cells[ncells+3].p = global_cells[ncells+1].p+dpdx*2.0*dx;
	    }
	}
	
	//-----------------------------------------------------
	// Interpolate cell centered values to interface
	//-----------------------------------------------------
	
	foreach(i; 0 .. ninterfaces) {
	    if (interpolation_order == 1){first_order_interpolation(global_cells[i+1],
								    global_cells[i+2],
								    global_interfaces[i]);}
	    else { second_order_interpolation_with_van_albada_limiting(global_cells[i+1],
								       global_cells[i],
								       global_cells[i+2],
								       global_cells[i+3],
								       global_interfaces[i]); }
	}
	
	//-----------------------------------------------------
	// Compute Interface Flux
	//-----------------------------------------------------
	
	foreach(i; 0 .. ninterfaces) {
	    if (flux_calc == "van_leer") van_leer_flux_calculator(global_interfaces[i], gamma);
	    else if (flux_calc == "ausmdv") ausmdv_flux_calculator(global_interfaces[i], gamma);
	    else if (flux_calc == "ausm_plus_up") {
		double M_inf = global_cells[0].u/(sqrt((global_cells[0].p*gamma)/global_cells[0].rho));
		ausm_plus_up_flux_calculator(global_interfaces[i], M_inf, gamma);
	    }
	}
	
	if (simulation_type == "nozzle") {
	    // for the nozzle simulation we must overwrite the inflow interface flux
	    
	    fvcell cell = global_cells[0];
	    double e = cell.p / (cell.rho*(gamma - 1.0));
	    double ke = 0.5*cell.u*cell.u;
	    
	    global_interfaces[0].mass = cell.rho*cell.u;
	    global_interfaces[0].momentum = cell.p+cell.rho*cell.u*cell.u;
	    global_interfaces[0].energy = (cell.rho*e + cell.rho*ke +cell.p)*cell.u;
	}
	
	//-----------------------------------------------------
	// Integrate flux and update cells
	//-----------------------------------------------------
	
	foreach(i; 0 .. ncells) {
	    fvcell cell = global_cells[i+2];
	    fvinterface fin = global_interfaces[i];
	    fvinterface fout = global_interfaces[i+1];
	    double delta_rho; double delta_momentum; double delta_energy;
		
	    // update mass
	    delta_rho = dt/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	    cell.rho = delta_rho + cell.rho;
	    
	    // update momentum
	    delta_momentum = dt/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area+cell.p*(fout.area-fin.area));
	    cell.ru = cell.ru + delta_momentum;
	    
	    // update total energy
	    delta_energy = dt/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	    cell.rE = cell.rE + delta_energy;
	    
	    decode_conserved_variables(global_cells[i+2], gamma);
	    //writef("%f, %f, %f \n", cell.rho, cell.u, cell.p);
	}
	
	// update time and step
	sim_time += dt;
	step += 1;
    } // end while loop

    if (solver == "simulation") return;
    
    //-----------------------------------------------------
    // Construct Jacobian
    //-----------------------------------------------------
    
    //-----------------------------------------------------
    // Apply Boundary Conditions
    //-----------------------------------------------------
    
    if (simulation_type == "sod") {
	// Sod's shocktube bcs
	
	// left boundary -- wall
	global_cells[0].rho = global_cells[2].rho; global_cells[1].rho = global_cells[2].rho;
	global_cells[0].u = -global_cells[2].u; global_cells[1].u = -global_cells[2].u;
	global_cells[0].p = global_cells[2].p; global_cells[1].p = global_cells[2].p;
	// right boundary -- wall
	global_cells[ncells+2].rho = global_cells[ncells+1].rho; global_cells[ncells+3].rho = global_cells[ncells+1].rho;
	global_cells[ncells+2].u = -global_cells[ncells+1].u; global_cells[ncells+3].u = -global_cells[ncells+1].u;
	global_cells[ncells+2].p = global_cells[ncells+1].p; global_cells[ncells+3].p = global_cells[ncells+1].p;
    }
    else if (simulation_type == "nozzle") {
	// Nozzle bcs
	
	// left boundary -- inflow
	global_cells[0].rho = global_cells[0].rho; global_cells[1].rho = global_cells[1].rho;
	global_cells[0].u = global_cells[0].u; global_cells[1].u = global_cells[1].u;
	global_cells[0].p = global_cells[0].p; global_cells[1].p = global_cells[1].p;
	
	if (outflow_bc_order == 1) {
	    // right boundary -- outflow
	    
	    // first order
	    global_cells[ncells+2].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+3].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+2].u = global_cells[ncells+1].u;
	    global_cells[ncells+3].u = global_cells[ncells+1].u;
	    global_cells[ncells+2].p = 2.5*global_cells[0].p;
	    global_cells[ncells+3].p = 2.5*global_cells[0].p;
	    //global_cells[ncells+2].p = global_cells[ncells+1].p;
	    //global_cells[ncells+3].p = global_cells[ncells+1].p;
	}
	else {
	    // second order
	    double drhodx; double dudx; double dpdx;
	    drhodx = (global_cells[ncells+1].rho - global_cells[ncells].rho)/dx;
	    dudx = (global_cells[ncells+1].u - global_cells[ncells].u)/dx;
	    dpdx = (global_cells[ncells+1].p - global_cells[ncells].p)/dx;
	    
	    global_cells[ncells+2].rho = global_cells[ncells+1].rho+drhodx*dx;
	    global_cells[ncells+3].rho = global_cells[ncells+1].rho+drhodx*2.0*dx;
	    global_cells[ncells+2].u = global_cells[ncells+1].u+dudx*dx;
	    global_cells[ncells+3].u = global_cells[ncells+1].u+dudx*2.0*dx;
	    global_cells[ncells+2].p = 2.5*global_cells[0].p+dpdx*dx;
	    global_cells[ncells+3].p = 2.5*global_cells[0].p+dpdx*2.0*dx;
	    //global_cells[ncells+2].p = global_cells[ncells+1].p;
	    //global_cells[ncells+3].p = global_cells[ncells+1].p;
	}
    }
    
    //-----------------------------------------------------
    // Interpolate cell centered values to interface
    //-----------------------------------------------------
    
    foreach(i; 0 .. ninterfaces) {
	if (interpolation_order == 1){first_order_interpolation(global_cells[i+1],
								global_cells[i+2],
								global_interfaces[i]);}
	else { second_order_interpolation_with_van_albada_limiting(global_cells[i+1],
								   global_cells[i],
								   global_cells[i+2],
								   global_cells[i+3],
								   global_interfaces[i]); }
    }
    
    //-----------------------------------------------------
    // Compute Interface Flux
    //-----------------------------------------------------
    
    foreach(i; 0 .. ninterfaces) {
	if (flux_calc == "van_leer") van_leer_flux_calculator(global_interfaces[i], gamma);
	else if (flux_calc == "ausmdv") ausmdv_flux_calculator(global_interfaces[i], gamma);
	else if (flux_calc == "ausm_plus_up") {
	    double M_inf = global_cells[0].u/(sqrt((global_cells[0].p*gamma)/global_cells[0].rho));
	    ausm_plus_up_flux_calculator(global_interfaces[i], M_inf, gamma);
	}
    }
    
    
    if (simulation_type == "nozzle") {
	// for the nozzle simulation we must overwrite the inflow interface flux
	
	fvcell cell = global_cells[0];
	double e = cell.p / (cell.rho*(gamma - 1.0));
	double ke = 0.5*cell.u*cell.u;
	
	global_interfaces[0].mass = cell.rho*cell.u;
	global_interfaces[0].momentum = cell.p+cell.rho*cell.u*cell.u;
	global_interfaces[0].energy = (cell.rho*e + cell.rho*ke +cell.p)*cell.u;
    }
    
    
    // Numerical Jacobian
    //-----------------------------------------------------
    // Construct R.H.S. residual vector (R)
    //-----------------------------------------------------
    
    foreach(i; 0 .. ncells) {
	fvinterface fin = global_interfaces[i];
	fvinterface fout = global_interfaces[i+1];
	fvcell cell = global_cells[i+2];
	// mass flux
	R[i*3] = -1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	// momentum flux
	R[i*3+1] = -1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	// total energy flux
	R[i*3+2] = -1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
    }
    //--------------------------------------------------------
    // Construct flow Jacobian w.r.t. primitive variables (JV)
    //--------------------------------------------------------
    
    foreach(i; 0 .. ncells) {
	// save current cells original flowstate
	
	double orig_rho = global_cells[i+2].rho;
	double orig_u = global_cells[i+2].u;
	double orig_p = global_cells[i+2].p;
	
	//--------------------------------
	// perturb density in current cell
	//--------------------------------
	
	global_cells[i+2].rho += (global_cells[i+2].rho*eps+eps);
	
	//-----------------------------------------------------
	// Apply Boundary Conditions
	//-----------------------------------------------------
	
	if (simulation_type == "sod") {
	    // Sod's shocktube bcs
	    
	    // left boundary -- wall
	    global_cells[0].rho = global_cells[2].rho;
	    global_cells[1].rho = global_cells[2].rho;
	    global_cells[0].u = -global_cells[2].u;
	    global_cells[1].u = -global_cells[2].u;
	    global_cells[0].p = global_cells[2].p;
	    global_cells[1].p = global_cells[2].p;
	    // right boundary -- wall
	    global_cells[ncells+2].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+3].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+2].u = -global_cells[ncells+1].u;
	    global_cells[ncells+3].u = -global_cells[ncells+1].u;
	    global_cells[ncells+2].p = global_cells[ncells+1].p;
	    global_cells[ncells+3].p = global_cells[ncells+1].p;
	}
	else if (simulation_type == "nozzle") {
	    // Nozzle bcs
	    
	    // left boundary -- inflow
	    global_cells[0].rho = global_cells[0].rho;
	    global_cells[1].rho = global_cells[1].rho;
	    global_cells[0].u = global_cells[0].u;
	    global_cells[1].u = global_cells[1].u;
	    global_cells[0].p = global_cells[0].p;
	    global_cells[1].p = global_cells[1].p;
	    
	    if (outflow_bc_order == 1) {
		// right boundary -- outflow
		
		// first order
		global_cells[ncells+2].rho = global_cells[ncells+1].rho;
		global_cells[ncells+3].rho = global_cells[ncells+1].rho;
		global_cells[ncells+2].u = global_cells[ncells+1].u;
		global_cells[ncells+3].u = global_cells[ncells+1].u;
		global_cells[ncells+2].p = 2.5*global_cells[0].p;
		global_cells[ncells+3].p = 2.5*global_cells[0].p;
	    }
	    else {
		// second order
		double drhodx; double dudx; double dpdx;
		drhodx = (global_cells[ncells+1].rho - global_cells[ncells].rho)/dx;
		dudx = (global_cells[ncells+1].u - global_cells[ncells].u)/dx;
		dpdx = (global_cells[ncells+1].p - global_cells[ncells].p)/dx;
		
		global_cells[ncells+2].rho = global_cells[ncells+1].rho+drhodx*dx;
		global_cells[ncells+3].rho = global_cells[ncells+1].rho+drhodx*2.0*dx;
		global_cells[ncells+2].u = global_cells[ncells+1].u+dudx*dx;
		global_cells[ncells+3].u = global_cells[ncells+1].u+dudx*2.0*dx;
		global_cells[ncells+2].p = global_cells[ncells+1].p+dpdx*dx;
		global_cells[ncells+3].p = global_cells[ncells+1].p+dpdx*2.0*dx;
	    }		    
	}
	
	//-----------------------------------------------------
	// Interpolate cell centered values to interface
	//-----------------------------------------------------
	
	foreach(j; 0 .. ninterfaces) {
	    if (interpolation_order == 1){first_order_interpolation(global_cells[j+1],
								    global_cells[j+2],
								    perturbed_interfaces[j]);}
	    else { second_order_interpolation_with_van_albada_limiting(global_cells[j+1],
								       global_cells[j],
								       global_cells[j+2],
								       global_cells[j+3],
								       perturbed_interfaces[j]); }
	}
	
	//-----------------------------------------------------
	// Compute Interface Flux
	//-----------------------------------------------------
	
	foreach(j; 0 .. ninterfaces) {
	    if (flux_calc == "van_leer") van_leer_flux_calculator(perturbed_interfaces[j], gamma);
	    else if (flux_calc == "ausmdv") ausmdv_flux_calculator(perturbed_interfaces[j], gamma);
	    else if (flux_calc == "ausm_plus_up") {
		double M_inf = global_cells[0].u/(sqrt((global_cells[0].p*gamma)/global_cells[0].rho));
		ausm_plus_up_flux_calculator(perturbed_interfaces[j], M_inf, gamma);
	    }
	}
	
	if (simulation_type == "nozzle") {
	    // for the nozzle simulation we must overwrite the inflow interface flux
	    
	    fvcell cell = global_cells[0];
	    double e = cell.p / (cell.rho*(gamma - 1.0));
	    double ke = 0.5*cell.u*cell.u;
	    
	    perturbed_interfaces[0].mass = cell.rho*cell.u;
	    perturbed_interfaces[0].momentum = cell.p+cell.rho*cell.u*cell.u;
	    perturbed_interfaces[0].energy = (cell.rho*e + cell.rho*ke +cell.p)*cell.u;
	}
	
	//-----------------------------------------------------
	// Fill column of Jacobian via Frechet derivative
	//-----------------------------------------------------
	
	foreach(j; 0 .. ncells) {
	    double resd0;
	    double resd1;
	    fvinterface fin = global_interfaces[j];
	    fvinterface fout = global_interfaces[j+1];
	    fvinterface pfin = perturbed_interfaces[j];
	    fvinterface pfout = perturbed_interfaces[j+1];
	    fvcell cell = global_cells[j+2];
	    
	    // mass flux
	    resd0 = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.mass*pfin.area - pfout.mass*pfout.area);
	    //if(i==0 && j==ncells-1) writeln(resd0, ", ", resd1, ", ",
	    //fin.mass, ", ", pfin.mass, ", ",
	    //pfout.mass, ", ", fout.mass);
	    JV[j*3][i*3] = (resd1-resd0)/(orig_rho*eps+eps);
	    
	    // fill transform matrix
	    if (i==j) transform[j*3][i*3] = 1.0;
	    else transform[j*3][i*3] = 0.0;
	    
	    // momentum flux
	    resd0 = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	    resd1 = 1.0/cell.vol * (pfin.momentum*pfin.area-pfout.momentum*pfout.area+ cell.p*(pfout.area-pfin.area));
	    JV[j*3+1][i*3] = (resd1-resd0)/(orig_rho*eps+eps);
	    
	    // fill transform matrix
	    if (i==j) transform[j*3+1][i*3] = -orig_u/orig_rho;
	    else transform[j*3+1][i*3] = 0.0;
	    
	    // total energy flux
	    resd0 = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.energy*pfin.area - pfout.energy*pfout.area);
	    JV[j*3+2][i*3] = (resd1-resd0)/(orig_rho*eps+eps);

	    //writeln("coords: ", j*3+2, ", ", i*3);
	    //writef("dFe/drho = %.16f", (pfout.energy-fout.energy)/(orig_rho*eps+eps));
	    //writef(", perturbed flux = %.16f", pfout.energy);
	    //writef(", unperturbed flux = %.16f \n", fout.energy);
	    	    
	    // fill transform matrix
	    if (i==j) transform[j*3+2][i*3] = 0.5*(gamma-1.0)*orig_u*orig_u;
	    else transform[j*3+2][i*3] = 0.0;
	}
	
	//--------------------------------
	// restore density in current cell
	//--------------------------------
	
	global_cells[i+2].rho = orig_rho;
	
	//--------------------------------
	// perturb velocity in current cell
	//--------------------------------
	
	global_cells[i+2].u += (global_cells[i+2].u*eps+eps);
	
	//-----------------------------------------------------
	// Apply Boundary Conditions
	//-----------------------------------------------------
	
	if (simulation_type == "sod") {
	    // Sod's shocktube bcs
	    
	    // left boundary -- wall
	    global_cells[0].rho = global_cells[2].rho;
	    global_cells[1].rho = global_cells[2].rho;
	    global_cells[0].u = -global_cells[2].u;
	    global_cells[1].u = -global_cells[2].u;
	    global_cells[0].p = global_cells[2].p;
	    global_cells[1].p = global_cells[2].p;
	    // right boundary -- wall
	    global_cells[ncells+2].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+3].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+2].u = -global_cells[ncells+1].u;
	    global_cells[ncells+3].u = -global_cells[ncells+1].u;
	    global_cells[ncells+2].p = global_cells[ncells+1].p;
	    global_cells[ncells+3].p = global_cells[ncells+1].p;
	}
	else if (simulation_type == "nozzle") {
	    // Nozzle bcs
	    
	    // left boundary -- inflow
	    global_cells[0].rho = global_cells[0].rho;
	    global_cells[1].rho = global_cells[1].rho;
	    global_cells[0].u = global_cells[0].u;
	    global_cells[1].u = global_cells[1].u;
	    global_cells[0].p = global_cells[0].p;
	    global_cells[1].p = global_cells[1].p;
	    
	    if (outflow_bc_order == 1) {
		// right boundary -- outflow
		
		// first order
		global_cells[ncells+2].rho = global_cells[ncells+1].rho;
		global_cells[ncells+3].rho = global_cells[ncells+1].rho;
		global_cells[ncells+2].u = global_cells[ncells+1].u;
		global_cells[ncells+3].u = global_cells[ncells+1].u;
		global_cells[ncells+2].p = 2.5*global_cells[0].p;
		global_cells[ncells+3].p = 2.5*global_cells[0].p;
	    }
	    else {
		// second order
		double drhodx; double dudx; double dpdx;
		drhodx = (global_cells[ncells+1].rho - global_cells[ncells].rho)/dx;
		dudx = (global_cells[ncells+1].u - global_cells[ncells].u)/dx;
		dpdx = (global_cells[ncells+1].p - global_cells[ncells].p)/dx;
		
		global_cells[ncells+2].rho = global_cells[ncells+1].rho+drhodx*dx;
		global_cells[ncells+3].rho = global_cells[ncells+1].rho+drhodx*2.0*dx;
		global_cells[ncells+2].u = global_cells[ncells+1].u+dudx*dx;
		global_cells[ncells+3].u = global_cells[ncells+1].u+dudx*2.0*dx;
		global_cells[ncells+2].p = global_cells[ncells+1].p+dpdx*dx;
		global_cells[ncells+3].p = global_cells[ncells+1].p+dpdx*2.0*dx;
	    }
	}
	
	//-----------------------------------------------------
	// Interpolate cell centered values to interface
	//-----------------------------------------------------
	
	foreach(j; 0 .. ninterfaces) {
	    if (interpolation_order == 1){first_order_interpolation(global_cells[j+1],
								    global_cells[j+2],
								    perturbed_interfaces[j]);}
	    else { second_order_interpolation_with_van_albada_limiting(global_cells[j+1],
								       global_cells[j],
								       global_cells[j+2],
								       global_cells[j+3],
								       perturbed_interfaces[j]); }
	}
	
	//-----------------------------------------------------
	// Compute Interface Flux
	//-----------------------------------------------------
	
	foreach(j; 0 .. ninterfaces) {
	    if (flux_calc == "van_leer") van_leer_flux_calculator(perturbed_interfaces[j], gamma);
	    else if (flux_calc == "ausmdv") ausmdv_flux_calculator(perturbed_interfaces[j], gamma);
	    else if (flux_calc == "ausm_plus_up") {
		double M_inf = global_cells[0].u/(sqrt((global_cells[0].p*gamma)/global_cells[0].rho));
		ausm_plus_up_flux_calculator(perturbed_interfaces[j], M_inf, gamma);
	    }
	}
	
	if (simulation_type == "nozzle") {
	    // for the nozzle simulation we must overwrite the inflow interface flux
	    
	    fvcell cell = global_cells[0];
	    double e = cell.p / (cell.rho*(gamma - 1.0));
	    double ke = 0.5*cell.u*cell.u;
	    
	    perturbed_interfaces[0].mass = cell.rho*cell.u;
	    perturbed_interfaces[0].momentum = cell.p+cell.rho*cell.u*cell.u;
	    perturbed_interfaces[0].energy = (cell.rho*e + cell.rho*ke +cell.p)*cell.u;
	}
	
	
	//-----------------------------------------------------
	// Fill column of Jacobian via Frechet derivative
	//-----------------------------------------------------
	
	foreach(j; 0 .. ncells) {
	    double resd0;
	    double resd1;
	    fvinterface fin = global_interfaces[j];
	    fvinterface fout = global_interfaces[j+1];
	    fvinterface pfin = perturbed_interfaces[j];
	    fvinterface pfout = perturbed_interfaces[j+1];
	    fvcell cell = global_cells[j+2];
	    
	    // mass flux
	    resd0 = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.mass*pfin.area - pfout.mass*pfout.area);
	    JV[j*3][i*3+1] = (resd1-resd0)/(orig_u*eps+eps);
	    
	    // fill transform matrix
	    if (i==j) transform[j*3][i*3+1] = 0.0;
	    else transform[j*3][i*3+1] = 0.0;
	    
	    // momentum flux
	    resd0 = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
	    resd1 = 1.0/cell.vol * (pfin.momentum*pfin.area-pfout.momentum*pfout.area+cell.p*(pfout.area-pfin.area));
	    JV[j*3+1][i*3+1] = (resd1-resd0)/(orig_u*eps+eps);
	    
	    // fill transform matrix
	    if (i==j) transform[j*3+1][i*3+1] = 1.0/orig_rho;
	    else transform[j*3+1][i*3+1] = 0.0;
	    
	    // total energy flux
	    resd0 = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.energy*pfin.area - pfout.energy*pfout.area);
	    JV[j*3+2][i*3+1] = (resd1-resd0)/(orig_u*eps+eps);

	    //writeln("coords = ", j*3+2, ", ", i*3);
	    //writef("dFe/du = %.16f", (pfout.energy-fout.energy)/(orig_u*eps+eps));
	    //writef(", perturbed flux = %.16f ", pfout.energy);
	    //writef(", unperturbed flux = %.16f \n", fout.energy);
	    // fill transform matrix
	    if (i==j) transform[j*3+2][i*3+1] = -(gamma-1.0)*orig_u;
	    else transform[j*3+2][i*3+1] = 0.0;
	}
	
	//--------------------------------
	// restore velocity in current cell
	//--------------------------------
	
	global_cells[i+2].u = orig_u;
	
	//--------------------------------
	// perturb pressure in current cell
	//--------------------------------
	
	global_cells[i+2].p += (global_cells[i+2].p*eps+eps);
	
	//-----------------------------------------------------
	// Apply Boundary Conditions
	//-----------------------------------------------------
	
	if (simulation_type == "sod") {
	    // Sod's shocktube bcs
	    
	    // left boundary -- wall
	    global_cells[0].rho = global_cells[2].rho;
	    global_cells[1].rho = global_cells[2].rho;
	    global_cells[0].u = -global_cells[2].u;
	    global_cells[1].u = -global_cells[2].u;
	    global_cells[0].p = global_cells[2].p;
	    global_cells[1].p = global_cells[2].p;
	    // right boundary -- wall
	    global_cells[ncells+2].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+3].rho = global_cells[ncells+1].rho;
	    global_cells[ncells+2].u = -global_cells[ncells+1].u;
	    global_cells[ncells+3].u = -global_cells[ncells+1].u;
	    global_cells[ncells+2].p = global_cells[ncells+1].p;
	    global_cells[ncells+3].p = global_cells[ncells+1].p;
	}
	else if (simulation_type == "nozzle") {
	    // Nozzle bcs
	    
	    // left boundary -- inflow
	    global_cells[0].rho = global_cells[0].rho;
	    global_cells[1].rho = global_cells[1].rho;
	    global_cells[0].u = global_cells[0].u;
	    global_cells[1].u = global_cells[1].u;
	    global_cells[0].p = global_cells[0].p;
	    global_cells[1].p = global_cells[1].p;
	    
	    if (outflow_bc_order == 1) {
		// right boundary -- outflow
		
		// first order
		global_cells[ncells+2].rho = global_cells[ncells+1].rho;
		global_cells[ncells+3].rho = global_cells[ncells+1].rho;
		global_cells[ncells+2].u = global_cells[ncells+1].u;
		global_cells[ncells+3].u = global_cells[ncells+1].u;
		global_cells[ncells+2].p = 2.5*global_cells[0].p;
		global_cells[ncells+3].p = 2.5*global_cells[0].p;
	    }
	    else {
		// second order
		double drhodx; double dudx; double dpdx;
		drhodx = (global_cells[ncells+1].rho - global_cells[ncells].rho)/dx;
		dudx = (global_cells[ncells+1].u - global_cells[ncells].u)/dx;
		dpdx = (global_cells[ncells+1].p - global_cells[ncells].p)/dx;
		
		global_cells[ncells+2].rho = global_cells[ncells+1].rho+drhodx*dx;
		global_cells[ncells+3].rho = global_cells[ncells+1].rho+drhodx*2.0*dx;
		global_cells[ncells+2].u = global_cells[ncells+1].u+dudx*dx;
		global_cells[ncells+3].u = global_cells[ncells+1].u+dudx*2.0*dx;
		global_cells[ncells+2].p = global_cells[ncells+1].p+dpdx*dx;
		global_cells[ncells+3].p = global_cells[ncells+1].p+dpdx*2.0*dx;
	    }
	}
	
	//-----------------------------------------------------
	// Interpolate cell centered values to interface
	//-----------------------------------------------------
	
	foreach(j; 0 .. ninterfaces) {
	    if (interpolation_order == 1){first_order_interpolation(global_cells[j+1],
								    global_cells[j+2],
								    perturbed_interfaces[j]);}
	    else { second_order_interpolation_with_van_albada_limiting(global_cells[j+1],
								       global_cells[j],
								       global_cells[j+2],
								       global_cells[j+3],
								       perturbed_interfaces[j]); }
	}
	
	//-----------------------------------------------------
	// Compute Interface Flux
	//-----------------------------------------------------
	
	foreach(j; 0 .. ninterfaces) {
	    if (flux_calc == "van_leer") van_leer_flux_calculator(perturbed_interfaces[j], gamma);
	    else if (flux_calc == "ausmdv") ausmdv_flux_calculator(perturbed_interfaces[j], gamma);
	    else if (flux_calc == "ausm_plus_up") {
		double M_inf = global_cells[0].u/(sqrt((global_cells[0].p*gamma)/global_cells[0].rho));
		ausm_plus_up_flux_calculator(perturbed_interfaces[j], M_inf, gamma);
	    }
	}
	
	if (simulation_type == "nozzle") {
	    // for the nozzle simulation we must overwrite the inflow interface flux
	    
	    fvcell cell = global_cells[0];
	    double e = cell.p / (cell.rho*(gamma - 1.0));
	    double ke = 0.5*cell.u*cell.u;
	    
	    perturbed_interfaces[0].mass = cell.rho*cell.u;
	    perturbed_interfaces[0].momentum = cell.p+cell.rho*cell.u*cell.u;
	    perturbed_interfaces[0].energy = (cell.rho*e + cell.rho*ke +cell.p)*cell.u;
	}
	
	//-----------------------------------------------------
	// Fill column of Jacobian via Frechet derivative
	//-----------------------------------------------------
	
	foreach(j; 0 .. ncells) {
	    double resd0;
	    double resd1;
	    fvinterface fin = global_interfaces[j];
	    fvinterface fout = global_interfaces[j+1];
	    fvinterface pfin = perturbed_interfaces[j];
	    fvinterface pfout = perturbed_interfaces[j+1];
	    fvcell cell = global_cells[j+2];
	    
	    // mass flux
	    resd0 = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.mass*pfin.area - pfout.mass*fout.area);
	    JV[j*3][i*3+2] = (resd1-resd0)/(orig_p*eps+eps);
	    
	    // fill transform matrix
	    if (i==j) transform[j*3][i*3+2] = 0.0;
	    else transform[j*3][i*3+2] = 0.0;
	    
	    // momentum flux
	    if (i == j) resd0 = 1.0/cell.vol *(fin.momentum*fin.area-fout.momentum*fout.area+orig_p*(fout.area-fin.area));
	    else resd0 = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area+cell.p*(fout.area-fin.area));
	    resd1 = 1.0/cell.vol * (pfin.momentum*pfin.area-pfout.momentum*pfout.area+cell.p*(pfout.area-pfin.area));
	    JV[j*3+1][i*3+2] = (resd1-resd0)/(orig_p*eps+eps);
	    
	    // fill transform matrix
	    if (i==j) transform[j*3+1][i*3+2] = 0.0;
	    else transform[j*3+1][i*3+2] = 0.0;
	    
	    // total energy flux
	    resd0 = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
	    resd1 = 1.0/cell.vol * (pfin.energy*pfin.area - pfout.energy*pfout.area);
	    JV[j*3+2][i*3+2] = (resd1-resd0)/(orig_p*eps+eps);
	    
	    // fill transform matrix
	    if (i==j) transform[j*3+2][i*3+2] = gamma-1.0;
	    else transform[j*3+2][i*3+2] = 0.0;
	}
	
	//--------------------------------
	// restore pressure in current cell
	//--------------------------------
	
	global_cells[i+2].p = orig_p;
	
    }
    
    /*
    // Transform matrix
    foreach(i; 0 .. ncells) {
    // save current cells original flowstate
    double orig_rho = global_cells[i+2].rho;
    double orig_u = global_cells[i+2].u;
    double orig_p = global_cells[i+2].p;
    
    foreach(j; 0 .. ncells) {
    // fill transform matrix
    if (i==j) transform[j*3][i*3] = 1.0;
    else transform[j*3][i*3] = 0.0;
    
    // fill transform matrix
    if (i==j) transform[j*3+1][i*3] = -orig_u/orig_rho;
    else transform[j*3+1][i*3] = 0.0;
    
    // fill transform matrix
    if (i==j) transform[j*3+2][i*3] = 0.5*(gamma-1.0)*orig_u*orig_u;
    else transform[j*3+2][i*3] = 0.0;
    }
    
    foreach(j; 0 .. ncells) {
    // fill transform matrix
    if (i==j) transform[j*3][i*3+1] = 0.0;
    else transform[j*3][i*3+1] = 0.0;
    
    // fill transform matrix
    if (i==j) transform[j*3+1][i*3+1] = 1.0/orig_rho;
    else transform[j*3+1][i*3+1] = 0.0;
    
    // fill transform matrix
    if (i==j) transform[j*3+2][i*3+1] = -(gamma-1.0)*orig_u;
    else transform[j*3+2][i*3+1] = 0.0;
    }
    
    foreach(j; 0 .. ncells) {
    // fill transform matrix
    if (i==j) transform[j*3][i*3+2] = 0.0;
    else transform[j*3][i*3+2] = 0.0;
    
    // fill transform matrix
    if (i==j) transform[j*3+1][i*3+2] = 0.0;
    else transform[j*3+1][i*3+2] = 0.0;
    
    // fill transform matrix
    if (i==j) transform[j*3+2][i*3+2] = gamma-1.0;
    else transform[j*3+2][i*3+2] = 0.0;
    }
    }
    
    // Analytical Jacobian
    foreach(i; 0..ncells) {
    
    //--------------------------------
    // Density derivatives
    //--------------------------------
    
    foreach(j; 0 .. ncells) {
    fvinterface fin = global_interfaces[j];
    fvinterface fout = global_interfaces[j+1];
    fvcell lftcell = global_cells[j+1];
    fvcell cell = global_cells[j+2];
    
    double in_dFmassdrho; double in_dFmomdrho; double in_dFedrho;
    double out_dFmassdrho; double out_dFmomdrho; double out_dFedrho;
    
    // first determine left interface Mach number
    // left state
    double aL = sqrt( (gamma*fin.left_fs.p)/fin.left_fs.rho );
    double ML = fin.left_fs.u/aL;
    
    // right state
    double aR = sqrt( (gamma*fin.right_fs.p)/fin.right_fs.rho );
    double MR = fin.right_fs.u/aR;
    
    double Min = 0.5*(ML+MR); // average Mach number
    
    // next determine right interface Mach number
    // left state
    aL = sqrt( (gamma*fout.left_fs.p)/fout.left_fs.rho );
    ML = fout.left_fs.u/aL;
    
    // right state
    aR = sqrt( (gamma*fout.right_fs.p)/fout.right_fs.rho );
    MR = fout.right_fs.u/aR;
    
    double Mout = 0.5*(ML+MR); // average Mach number
    
    if (i == j-1) {
    if (Min >= 1.0) {
    double uL = fin.left_fs.u;
    in_dFmassdrho = uL;
    in_dFmomdrho = uL*uL;
    in_dFedrho = 0.5*uL*uL*uL;
    }
    else if (Min <= -1.0) {	
    in_dFmassdrho = 0.0;
    in_dFmomdrho = 0.0;
    in_dFedrho = 0.0;
    }
    else {
    double rhoL = fin.left_fs.rho; double rhoR = fin.right_fs.rho;
    double uL = fin.left_fs.u; double uR = fin.right_fs.u;
    double pL = fin.left_fs.p; double pR = fin.right_fs.p;
    double g = gamma;
    
    in_dFmassdrho = 0.25*uL*(uL/sqrt(g*pL/rhoL)+1)+0.125*sqrt(g*pL/rhoL)*pow((uL/sqrt(g*pL/rhoL) + 1),2);
    in_dFmomdrho = 0.25*pL*uL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)/(rhoL*sqrt(g*pL/rhoL))
    + 0.5*pL*uL*(uL/sqrt(g*pL/rhoL) + 1)*(uL*(g/2.0- 0.5)/sqrt(g*pL/rhoL) + 1)/(rhoL*sqrt(g*pL/rhoL)); 
    in_dFedrho = 0.5*g*pL*uL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/(rhoL*(pow(g,2) - 1))
    + 0.5*g*pL*uL*(uL/sqrt(g*pL/rhoL)+1)*pow((uL*(g/2.0-0.5)/sqrt(g*pL/rhoL)+1),2)/(rhoL*(pow(g,2)-1))
    - 0.25*g*pL*sqrt(g*pL/rhoL)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *pow((uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1),2)/(rhoL*(pow(g,2) - 1));
    }
    // mass flux
    JV[j*3][i*3] = -1.0/cell.vol * (in_dFmassdrho*fin.area);
    
    // momentum flux
    JV[j*3+1][i*3] = -1.0/cell.vol * (in_dFmomdrho*fin.area);
    
    // energy flux
    JV[j*3+2][i*3] = -1.0/cell.vol * (in_dFedrho*fin.area);
    }
    
    else if (i == j+1) {
    if (Mout >= 1.0) {
    out_dFmassdrho = 0.0;
    out_dFmomdrho = 0.0;
    out_dFedrho = 0.0;
    }
    else if (Mout <= -1.0) {
    double uR = fout.right_fs.u;
    out_dFmassdrho = uR;
    out_dFmomdrho = uR*uR;
    out_dFedrho = 0.5*uR*uR*uR;
    }
    else {
    double rhoL = fout.left_fs.rho; double rhoR = fout.right_fs.rho;
    double uL = fout.left_fs.u; double uR = fout.right_fs.u;
    double pL = fout.left_fs.p; double pR = fout.right_fs.p;
    double g = gamma;
    
    out_dFmassdrho =0.25*uR*(-uR/sqrt(g*pR/rhoR)+1)-0.125*sqrt(g*pR/rhoR)*pow((-uR/sqrt(g*pR/rhoR) + 1),2);
    out_dFmomdrho = -0.25*pR*uR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)/(rhoR*sqrt(g*pR/rhoR))
    + 0.5*pR*uR*(-uR/sqrt(g*pR/rhoR) + 1)*(uR*(g/2.0-0.5)/sqrt(g*pR/rhoR) - 1)/(rhoR*sqrt(g*pR/rhoR));
    out_dFedrho = -0.5*g*pR*uR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/(rhoR*(pow(g,2) - 1))
    + 0.5*g*pR*uR*(-uR/sqrt(g*pR/rhoR) + 1)
    *pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(rhoR*(pow(g,2) - 1))
    + 0.25*g*pR*sqrt(g*pR/rhoR)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(rhoR*(pow(g,2) - 1));
    }
    // mass flux
    JV[j*3][i*3] = -1.0/cell.vol * (-out_dFmassdrho*fout.area);
    
    // momentum flux
    JV[j*3+1][i*3] = -1.0/cell.vol * (-out_dFmomdrho*fout.area);
    
    // energy flux
    JV[j*3+2][i*3] = -1.0/cell.vol * (-out_dFedrho*fout.area);
    }
    
    else if (i == j) {
    if (Min >= 1.0) {
    in_dFmassdrho = 0.0;
    in_dFmomdrho = 0.0;
    in_dFedrho = 0.0;
    }
    else if (Min <= -1.0) {
    double uR = fin.right_fs.u;
    in_dFmassdrho = uR;
    in_dFmomdrho = uR*uR;
    in_dFedrho = 0.5*uR*uR*uR;
    }
    else {
    double rhoL = fin.left_fs.rho; double rhoR = fin.right_fs.rho;
    double uL = fin.left_fs.u; double uR = fin.right_fs.u;
    double pL = fin.left_fs.p; double pR = fin.right_fs.p;
    double g = gamma;
    
    in_dFmassdrho = 0.25*uR*(-uR/sqrt(g*pR/rhoR)+1)-0.125*sqrt(g*pR/rhoR)*pow((-uR/sqrt(g*pR/rhoR) + 1),2);
    in_dFmomdrho = -0.25*pR*uR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)/(rhoR*sqrt(g*pR/rhoR))
    + 0.5*pR*uR*(-uR/sqrt(g*pR/rhoR) + 1)*(uR*(g/2.0-0.5)/sqrt(g*pR/rhoR) - 1)/(rhoR*sqrt(g*pR/rhoR));
    in_dFedrho = -0.5*g*pR*uR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/(rhoR*(pow(g,2) - 1))
    + 0.5*g*pR*uR*(-uR/sqrt(g*pR/rhoR) + 1)
    *pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(rhoR*(pow(g,2) - 1))
    + 0.25*g*pR*sqrt(g*pR/rhoR)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(rhoR*(pow(g,2) - 1));
    }
    
    if (Mout >= 1.0) {
    double uL = fout.left_fs.u;
    out_dFmassdrho = uL;
    out_dFmomdrho = uL*uL;
    out_dFedrho = 0.5*uL*uL*uL;
    }
    else if (Mout <= -1.0) {	
    out_dFmassdrho = 0.0;
    out_dFmomdrho = 0.0;
    out_dFedrho = 0.0;
    }
    else {
    double rhoL = fout.left_fs.rho; double rhoR = fout.right_fs.rho;
    double uL = fout.left_fs.u; double uR = fout.right_fs.u;
    double pL = fout.left_fs.p; double pR = fout.right_fs.p;
    double g = gamma;
    
    out_dFmassdrho = 0.25*uL*(uL/sqrt(g*pL/rhoL)+1)+0.125*sqrt(g*pL/rhoL)*pow((uL/sqrt(g*pL/rhoL) + 1),2);
    out_dFmomdrho = 0.25*pL*uL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)/(rhoL*sqrt(g*pL/rhoL))
    + 0.5*pL*uL*(uL/sqrt(g*pL/rhoL) + 1)*(uL*(g/2.0-0.5)/sqrt(g*pL/rhoL) + 1)/(rhoL*sqrt(g*pL/rhoL)); 
    out_dFedrho = 0.5*g*pL*uL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/(rhoL*(pow(g,2) - 1))
    + 0.5*g*pL*uL*(uL/sqrt(g*pL/rhoL)+1)*pow((uL*(g/2.0-0.5)/sqrt(g*pL/rhoL)+1),2)/(rhoL*(pow(g,2)-1))
    - 0.25*g*pL*sqrt(g*pL/rhoL)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *pow((uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1),2)/(rhoL*(pow(g,2) - 1));
    }
    
    // mass flux
    JV[j*3][i*3] = -1.0/cell.vol * (in_dFmassdrho*fin.area-out_dFmassdrho*fout.area);
    
    // momentum flux
    JV[j*3+1][i*3] = -1.0/cell.vol * (in_dFmomdrho*fin.area-out_dFmomdrho*fout.area);
    
    // energy flux
    JV[j*3+2][i*3] = -1.0/cell.vol * (in_dFedrho*fin.area-out_dFedrho*fout.area);
    }
    
    else {
    // mass flux
    JV[j*3][i*3] = 0.0;
    
    // momentum flux
    JV[j*3+1][i*3] = 0.0;
    
    // energy flux
    JV[j*3+2][i*3] = 0.0;
    }
    }
    //--------------------------------
    // velocity derivatives
    //--------------------------------
    foreach(j; 0 .. ncells) {
    fvinterface fin = global_interfaces[j];
    fvinterface fout = global_interfaces[j+1];
    fvcell lftcell = global_cells[j+1];
    fvcell cell = global_cells[j+2];
    
    double in_dFmassdu; double in_dFmomdu; double in_dFedu;
    double out_dFmassdu; double out_dFmomdu; double out_dFedu;
    
    // first determine left interface Mach number
    // left state
    double aL = sqrt( (gamma*fin.left_fs.p)/fin.left_fs.rho );
    double ML = fin.left_fs.u/aL;
    
    // right state
    double aR = sqrt( (gamma*fin.right_fs.p)/fin.right_fs.rho );
    double MR = fin.right_fs.u/aR;
    
    double Min = 0.5*(ML+MR); // average Mach number
    
    // next determine right interface Mach number
    // left state
    aL = sqrt( (gamma*fout.left_fs.p)/fout.left_fs.rho );
    ML = fout.left_fs.u/aL;
    
    // right state
    aR = sqrt( (gamma*fout.right_fs.p)/fout.right_fs.rho );
    MR = fout.right_fs.u/aR;
    
    double Mout = 0.5*(ML+MR); // average Mach number
    
    if (i == j-1) {
    if (Min >= 1.0) {
    double rhoL = fin.left_fs.rho;
    double uL = fin.left_fs.u;
    double pL = fin.left_fs.p;
    
    in_dFmassdu = rhoL;
    in_dFmomdu = 2*rhoL*uL;
    in_dFedu = (3.0/2.0) * rhoL*uL*uL + (pL*gamma)/(gamma-1);
    }
    else if (Min <= -1.0) {	
    in_dFmassdu = 0.0;
    in_dFmomdu = 0.0;
    in_dFedu = 0.0;
    }
    else {
    double rhoL = fin.left_fs.rho; double rhoR = fin.right_fs.rho;
    double uL = fin.left_fs.u; double uR = fin.right_fs.u;
    double pL = fin.left_fs.p; double pR = fin.right_fs.p;
    double g = gamma;
    
    in_dFmassdu = 0.5*rhoL*(uL/sqrt(g*pL/rhoL) + 1);
    in_dFmomdu = 0.5*pL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)/sqrt(g*pL/rhoL)
    + 1.0*pL*(uL/sqrt(g*pL/rhoL) + 1)*(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/sqrt(g*pL/rhoL);
    in_dFedu = 1.0*g*pL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/(pow(g,2) - 1)
    + 1.0*g*pL*(uL/sqrt(g*pL/rhoL) + 1)*pow((uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1),2)/(pow(g,2) - 1);
    }
    // mass flux
    JV[j*3][i*3+1] = -1.0/cell.vol * (in_dFmassdu*fin.area);
    
    // momentum flux
    JV[j*3+1][i*3+1] = -1.0/cell.vol * (in_dFmomdu*fin.area);
    
    // energy flux
    JV[j*3+2][i*3+1] = -1.0/cell.vol * (in_dFedu*fin.area);
    }
    
    else if (i == j+1) {
    if (Mout >= 1.0) {
    out_dFmassdu = 0.0;
    out_dFmomdu = 0.0;
    out_dFedu = 0.0;
    }
    else if (Mout <= -1.0) {
    double rhoR = fout.right_fs.rho;
    double uR = fout.right_fs.u;
    double pR = fout.right_fs.p;
    
    out_dFmassdu = rhoR;
    out_dFmomdu = 2.0*rhoR*uR;
    out_dFedu = (3.0/2.0) *rhoR*uR*uR + (pR*gamma)/(gamma-1);
    }
    else {
    double rhoL = fout.left_fs.rho; double rhoR = fout.right_fs.rho;
    double uL = fout.left_fs.u; double uR = fout.right_fs.u;
    double pL = fout.left_fs.p; double pR = fout.right_fs.p;
    double g = gamma;
    
    out_dFmassdu = 0.5*rhoR*(-uR/sqrt(g*pR/rhoR) + 1);
    out_dFmomdu = -0.5*pR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)/sqrt(g*pR/rhoR)
    + 1.0*pR*(-uR/sqrt(g*pR/rhoR) + 1)*(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/sqrt(g*pR/rhoR);
    out_dFedu = -1.0*g*pR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/(pow(g,2) - 1)
    + 1.0*g*pR*(-uR/sqrt(g*pR/rhoR) + 1)*pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(pow(g,2) - 1);
    }
    // mass flux
    JV[j*3][i*3+1] = -1.0/cell.vol * (-out_dFmassdu*fout.area);
    
    // momentum flux
    JV[j*3+1][i*3+1] = -1.0/cell.vol * (-out_dFmomdu*fout.area);
    
    // energy flux
    JV[j*3+2][i*3+1] = -1.0/cell.vol * (-out_dFedu*fout.area);
    }
    
    else if (i == j) {
    if (Min >= 1.0) {
    in_dFmassdu = 0.0;
    in_dFmomdu = 0.0;
    in_dFedu = 0.0;
    }
    else if (Min <= -1.0) {
    double rhoR = fin.right_fs.rho;
    double uR = fin.right_fs.u;
    double pR = fin.right_fs.p;
    
    in_dFmassdu = rhoR;
    in_dFmomdu = 2.0*rhoR*uR;
    in_dFedu = (3.0/2.0) *rhoR*uR*uR + (pR*gamma)/(gamma-1);
    }
    else {
    double rhoL = fin.left_fs.rho; double rhoR = fin.right_fs.rho;
    double uL = fin.left_fs.u; double uR = fin.right_fs.u;
    double pL = fin.left_fs.p; double pR = fin.right_fs.p;
    double g = gamma;
    
    in_dFmassdu = 0.5*rhoR*(-uR/sqrt(g*pR/rhoR) + 1);
    in_dFmomdu = -0.5*pR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)/sqrt(g*pR/rhoR)
    + 1.0*pR*(-uR/sqrt(g*pR/rhoR) + 1)*(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/sqrt(g*pR/rhoR);
    in_dFedu = -1.0*g*pR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/(pow(g,2) - 1)
    + 1.0*g*pR*(-uR/sqrt(g*pR/rhoR) + 1)*pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(pow(g,2) - 1);
    }
    
    if (Mout >= 1.0) {
    double rhoL = fout.left_fs.rho;
    double uL = fout.left_fs.u;
    double pL = fout.left_fs.p;
    
    out_dFmassdu = rhoL;
    out_dFmomdu = 2*rhoL*uL;
    out_dFedu = (3.0/2.0) *rhoL*uL*uL + (pL*gamma)/(gamma-1);
    }
    else if (Mout <= -1.0) {	
    out_dFmassdu = 0.0;
    out_dFmomdu = 0.0;
    out_dFedu = 0.0;
    }
    else {
    double rhoL = fout.left_fs.rho; double rhoR = fout.right_fs.rho;
    double uL = fout.left_fs.u; double uR = fout.right_fs.u;
    double pL = fout.left_fs.p; double pR = fout.right_fs.p;
    double g = gamma;
    
    out_dFmassdu = 0.5*rhoL*(uL/sqrt(g*pL/rhoL) + 1);
    out_dFmomdu = 0.5*pL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)/sqrt(g*pL/rhoL)
    + 1.0*pL*(uL/sqrt(g*pL/rhoL) + 1)*(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/sqrt(g*pL/rhoL);
    out_dFedu = 1.0*g*pL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/(pow(g,2) - 1)
    + 1.0*g*pL*(uL/sqrt(g*pL/rhoL) + 1)*pow((uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1),2)/(pow(g,2) - 1);
    }
    
    // mass flux
    JV[j*3][i*3+1] = -1.0/cell.vol * (in_dFmassdu*fin.area-out_dFmassdu*fout.area);
    
    // momentum flux
    JV[j*3+1][i*3+1] = -1.0/cell.vol * (in_dFmomdu*fin.area-out_dFmomdu*fout.area);
    
    // energy flux
    JV[j*3+2][i*3+1] = -1.0/cell.vol * (in_dFedu*fin.area-out_dFedu*fout.area);
    }
    
    else {
    // mass flux
    JV[j*3][i*3+1] = 0.0;
    
    // momentum flux
    JV[j*3+1][i*3+1] = 0.0;
    
    // energy flux
    JV[j*3+2][i*3+1] = 0.0;
    }
    }
    //--------------------------------
    // pressure derivatives
    //--------------------------------
    foreach(j; 0 .. ncells) {
    fvinterface fin = global_interfaces[j];
    fvinterface fout = global_interfaces[j+1];
    fvcell lftcell = global_cells[j+1];
    fvcell cell = global_cells[j+2];
    
    double in_dFmassdp; double in_dFmomdp; double in_dFedp;
    double out_dFmassdp; double out_dFmomdp; double out_dFedp;
    
    // first determine left interface Mach number
    // left state
    double aL = sqrt( (gamma*fin.left_fs.p)/fin.left_fs.rho );
    double ML = fin.left_fs.u/aL;
    
    // right state
    double aR = sqrt( (gamma*fin.right_fs.p)/fin.right_fs.rho );
    double MR = fin.right_fs.u/aR;
    
    double Min = 0.5*(ML+MR); // average Mach number
    
    // next determine right interface Mach number
    // left state
    aL = sqrt( (gamma*fout.left_fs.p)/fout.left_fs.rho );
    ML = fout.left_fs.u/aL;
    
    // right state
    aR = sqrt( (gamma*fout.right_fs.p)/fout.right_fs.rho );
    MR = fout.right_fs.u/aR;
    
    double Mout = 0.5*(ML+MR); // average Mach number
    
    if (i == j-1) {
    if (Min >= 1.0) {
    double rhoL = fin.left_fs.rho;
    double uL = fin.left_fs.u;
    double pL = fin.left_fs.p;
    
    in_dFmassdp = 0.0;
    in_dFmomdp = 1.0;
    in_dFedp = (uL*gamma)/(gamma-1);
    }
    else if (Min <= -1.0) {	
    in_dFmassdp = 0.0;
    in_dFmomdp = 0.0;
    in_dFedp = 0.0;
    }
    else {
    double rhoL = fin.left_fs.rho; double rhoR = fin.right_fs.rho;
    double uL = fin.left_fs.u; double uR = fin.right_fs.u;
    double pL = fin.left_fs.p; double pR = fin.right_fs.p;
    double g = gamma;
    
    in_dFmassdp = -0.25*rhoL*uL*(uL/sqrt(g*pL/rhoL) + 1)/pL
    + 0.125*rhoL*sqrt(g*pL/rhoL)*pow((uL/sqrt(g*pL/rhoL) + 1),2)/pL;
    in_dFmomdp = -0.25*uL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)/sqrt(g*pL/rhoL)
    - 0.5*uL*(uL/sqrt(g*pL/rhoL) + 1)*(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/sqrt(g*pL/rhoL)
    + 0.5*pow((uL/sqrt(g*pL/rhoL) + 1),2)*(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1);
    in_dFedp = -0.5*g*uL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/(pow(g,2) - 1) - 0.5*g*uL*(uL/sqrt(g*pL/rhoL) + 1)
    *pow((uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1),2)/(pow(g,2) - 1)
    + 0.75*g*sqrt(g*pL/rhoL)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *pow((uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1),2)/(pow(g,2) - 1);
    }
    // mass flux
    JV[j*3][i*3+2] = -1.0/cell.vol * (in_dFmassdp*fin.area);
    
    // momentum flux
    JV[j*3+1][i*3+2] = -1.0/cell.vol * (in_dFmomdp*fin.area);
    
    // energy flux
    JV[j*3+2][i*3+2] = -1.0/cell.vol * (in_dFedp*fin.area);
    }
	
    else if (i == j+1) {
    if (Mout >= 1.0) {
    out_dFmassdp = 0.0;
    out_dFmomdp = 0.0;
    out_dFedp = 0.0;
    }
    else if (Mout <= -1.0) {
    double rhoR = fout.right_fs.rho;
    double uR = fout.right_fs.u;
    double pR = fout.right_fs.p;
    
    out_dFmassdp = 0.0;
    out_dFmomdp = 1.0;
    out_dFedp = (uR*gamma)/(gamma-1);
    }
    else {
    double rhoL = fout.left_fs.rho; double rhoR = fout.right_fs.rho;
    double uL = fout.left_fs.u; double uR = fout.right_fs.u;
    double pL = fout.left_fs.p; double pR = fout.right_fs.p;
    double g = gamma;
    
    out_dFmassdp = -0.25*rhoR*uR*(-uR/sqrt(g*pR/rhoR) + 1)/pR
    - 0.125*rhoR*sqrt(g*pR/rhoR)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)/pR;
    out_dFmomdp = 0.25*uR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)/sqrt(g*pR/rhoR)
    - 0.5*uR*(-uR/sqrt(g*pR/rhoR) + 1)*(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/sqrt(g*pR/rhoR)
    - 0.5*pow((-uR/sqrt(g*pR/rhoR) + 1),2)*(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1);
    out_dFedp = 0.5*g*uR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/(pow(g,2) - 1)
    - 0.5*g*uR*(-uR/sqrt(g*pR/rhoR) + 1)*pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(pow(g,2) - 1)
    - 0.75*g*sqrt(g*pR/rhoR)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(pow(g,2) - 1);
    }
    // mass flux
    JV[j*3][i*3+2] = -1.0/cell.vol * (-out_dFmassdp*fout.area);
    
    // momentum flux
    JV[j*3+1][i*3+2] = -1.0/cell.vol * (-out_dFmomdp*fout.area);
    
    // energy flux
    JV[j*3+2][i*3+2] = -1.0/cell.vol * (-out_dFedp*fout.area);
    }
	
    else if (i == j) {
    if (Min >= 1.0) {
    in_dFmassdp = 0.0;
    in_dFmomdp = 0.0;
    in_dFedp = 0.0;
    }
    else if (Min <= -1.0) {
    double rhoR = fin.right_fs.rho;
    double uR = fin.right_fs.u;
    double pR = fin.right_fs.p;
    
    in_dFmassdp = 0.0;
    in_dFmomdp = 1.0;
    in_dFedp = (uR*gamma)/(gamma-1);
    }
    else {
    double rhoL = fin.left_fs.rho; double rhoR = fin.right_fs.rho;
    double uL = fin.left_fs.u; double uR = fin.right_fs.u;
    double pL = fin.left_fs.p; double pR = fin.right_fs.p;
    double g = gamma;
	
    in_dFmassdp = -0.25*rhoR*uR*(-uR/sqrt(g*pR/rhoR) + 1)/pR
    - 0.125*rhoR*sqrt(g*pR/rhoR)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)/pR;
    in_dFmomdp = 0.25*uR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)/sqrt(g*pR/rhoR)
    - 0.5*uR*(-uR/sqrt(g*pR/rhoR) + 1)*(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/sqrt(g*pR/rhoR)
    - 0.5*pow((-uR/sqrt(g*pR/rhoR) + 1),2)*(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1);
    in_dFedp = 0.5*g*uR*(g/2.0 - 0.5)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *(uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1)/(pow(g,2) - 1)
    - 0.5*g*uR*(-uR/sqrt(g*pR/rhoR) + 1)*pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(pow(g,2) - 1)
    - 0.75*g*sqrt(g*pR/rhoR)*pow((-uR/sqrt(g*pR/rhoR) + 1),2)
    *pow((uR*(g/2.0 - 0.5)/sqrt(g*pR/rhoR) - 1),2)/(pow(g,2) - 1);
    }
    
    if (Mout >= 1.0) {
    double rhoL = fout.left_fs.rho;
    double uL = fout.left_fs.u;
    double pL = fout.left_fs.p;
    
    out_dFmassdp = 0.0;
    out_dFmomdp = 1.0;
    out_dFedp = (uL*gamma)/(gamma-1);
    }
    else if (Mout <= -1.0) {	
    out_dFmassdp = 0.0;
    out_dFmomdp = 0.0;
    out_dFedp = 0.0;
    }
    else {
    double rhoL = fout.left_fs.rho; double rhoR = fout.right_fs.rho;
    double uL = fout.left_fs.u; double uR = fout.right_fs.u;
    double pL = fout.left_fs.p; double pR = fout.right_fs.p;
    double g = gamma;
    
    out_dFmassdp = -0.25*rhoL*uL*(uL/sqrt(g*pL/rhoL) + 1)/pL
    + 0.125*rhoL*sqrt(g*pL/rhoL)*pow((uL/sqrt(g*pL/rhoL) + 1),2)/pL;
    out_dFmomdp = -0.25*uL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)/sqrt(g*pL/rhoL)
    - 0.5*uL*(uL/sqrt(g*pL/rhoL) + 1)*(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/sqrt(g*pL/rhoL)
    + 0.5*pow((uL/sqrt(g*pL/rhoL) + 1),2)*(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1);
    out_dFedp = -0.5*g*uL*(g/2.0 - 0.5)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *(uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1)/(pow(g,2) - 1) - 0.5*g*uL*(uL/sqrt(g*pL/rhoL) + 1)
    *pow((uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1),2)/(pow(g,2) - 1)
    + 0.75*g*sqrt(g*pL/rhoL)*pow((uL/sqrt(g*pL/rhoL) + 1),2)
    *pow((uL*(g/2.0 - 0.5)/sqrt(g*pL/rhoL) + 1),2)/(pow(g,2) - 1);
    }
    
    // mass flux
    JV[j*3][i*3+2] = -1.0/cell.vol * (in_dFmassdp*fin.area-out_dFmassdp*fout.area);
    
    // momentum flux
    JV[j*3+1][i*3+2] = -1.0/cell.vol * (in_dFmomdp*fin.area-out_dFmomdp*fout.area +1.0*(fout.area-fin.area));
    
    // energy flux
    JV[j*3+2][i*3+2] = -1.0/cell.vol * (in_dFedp*fin.area-out_dFedp*fout.area);
    }
    
    else {
    // mass flux
    JV[j*3][i*3+2] = 0.0;
    
    // momentum flux
    JV[j*3+1][i*3+2] = 0.0;
    
    // energy flux
    JV[j*3+2][i*3+2] = 0.0;
    }
    }
    } // end analytical Jacobian
    */
    //writeln(JV);
    //-----------------------------------------------------
    // Transform JV to JU
    //-----------------------------------------------------
    //writeln("J = ", JV);
    //writeln(JV.length);
    matrixMult(JV, transform, JU); 
    
    //--------------------------------------------------------
    // Adjoint Solver -- [JU]^T * psi = -dJdU
    //--------------------------------------------------------
    
    // cost function is defined as:
    // J(Q) = 0.5*integral[0->l] (p-p*)^2 dx
    
    // target pressure distribution saved in file target.dat
    double[ncells] p_target;
    auto file = File("target_output.dat", "r");
    foreach(i; 0 .. ncells) {
	auto lineContent = file.readln().strip();
	auto tokens = lineContent.split();
	p_target[i] = to!double(tokens[2]);
    }
    
    //--------------------------------------------------------
    // Transpose JU
    //--------------------------------------------------------
    
    // transpose
    foreach (i; 0 .. JU.length) {
	foreach (j; 0 .. JU.length) {
		JUT[i][j] = JU[j][i];
		JVT[i][j] = JV[j][i];
	}
    }
    //writeln("--------------------------------------------------");
    //writeln("JVT = ", JVT);
    //--------------------------------------------------------
    // Invert [JU]^T
    //--------------------------------------------------------
    
    matrixInv(JUT, invJUT);
    matrixInv(JVT, invJVT);
    //writeln(invJVT);
    //--------------------------------------------------------
    // Construct dJdV vector
    //--------------------------------------------------------
	
    // first calculate un-perturbed (J)
    //double J = 0.0;
    J = 0.0;
    foreach (i; 0..ncells) {
	J += 0.5*(global_cells[i+2].p - p_target[i])*(global_cells[i+2].p - p_target[i]); //*global_cells[i+2].dx;
    }
    //writeln(J, ", ", b, ", ", c, ", ", d, ", ", p_target[50], ", ", global_cells[52].p, ", ", global_cells[52].dx);
    //writef("%.16f, %f, %f, %f \n", J, b, c, d);
    //writef("J = %.18f \n", J);
    
    /*
    // Numerically form dJdV using Frechet Derivatives
    
    double J_perturb;	    
    
    // now loop through cells and perturb primitive variables
    foreach (i; 0 .. ncells) {
    // store original cell flowstate
    double orig_rho = global_cells[i+2].rho;
    double orig_u = global_cells[i+2].u;
    double orig_p = global_cells[i+2].p;
    
    // perturb density
    global_cells[i+2].rho += (global_cells[i+2].rho*eps+eps);
    J_perturb = 0.0;
    foreach (j; 0..ncells) {
    J_perturb += 0.5*(global_cells[j+2].p - p_target[j])*(global_cells[j+2].p - p_target[j])*global_cells[j+2].dx;
    }
    dJdV[i*3] = (J_perturb-J)/(orig_rho*eps+eps);
    global_cells[i+2].rho = orig_rho;
    
    // perturb velocity
    global_cells[i+2].u += (global_cells[i+2].u*eps+eps);
    J_perturb = 0.0;
    foreach (j; 0..ncells) {
    J_perturb += 0.5*(global_cells[j+2].p - p_target[j])*(global_cells[j+2].p - p_target[j])*global_cells[j+2].dx;
    }
    dJdV[i*3+1] = (J_perturb-J)/(orig_u*eps+eps);
    global_cells[i+2].u = orig_u;
    
    // perturb pressure
    global_cells[i+2].p += (global_cells[i+2].p*eps+eps);
    J_perturb = 0.0;
    foreach (j; 0..ncells) {
    J_perturb += 0.5*(global_cells[j+2].p - p_target[j])*(global_cells[j+2].p - p_target[j])*global_cells[j+2].dx;
    }
    dJdV[i*3+2] = (J_perturb-J)/(orig_p*eps+eps);
    global_cells[i+2].p = orig_p;
    }
    */
    
    // Analytically form dJdV by hand differentiation
    foreach (i; 0..ncells) {
	foreach (j; 0..ncells) {
	    dJdV[i*3] = 0.0;
	    dJdV[i*3+1] = 0.0;
	    dJdV[i*3+2] = 0.5*(2.0*global_cells[i+2].p-2.0*p_target[i]); //global_cells[i+2].dx
	}
    }
    //writeln("-------------------------------------------------------------------");
    //writeln("dJdV = ", dJdV);
    /*
    //--------------------------------------------------------
    // Transform dJdV to dJdU
    //--------------------------------------------------------
	
    foreach (i; 0 .. transform.length) {
	dJdU[i] = 0.0;
	foreach (j; 0 .. transform.length) {
	    dJdU[i] += dJdV[j] * transform[j][i];
	}
    }
    */
    
    
    // We can also construct dJdU directly by hand differentiating J w.r.t. the conserved variables
    foreach (i; 0..ncells) {
	double U1 = global_cells[i+2].rho; double U2 = global_cells[i+2].ru; double U3 = global_cells[i+2].rE;
	double delx = global_cells[i+2].dx; double pstar = p_target[i];
	dJdU[i*3] = delx/2.0 * (U3*(gamma-1.0)*(gamma-1.0)*U2*U2/(U1*U1)
				-U2*U2*U2*U2*(gamma-1.0)*(gamma-1.0)/(2.0*U1*U1*U1)
				-U2*U2/(U1*U1)*(gamma-1.0)*pstar);
	dJdU[i*3+1] = delx/2.0 * (-2.0*U3*(gamma-1.0)*(gamma-1.0)*U2/U1
				  +U2*U2*U2*(gamma-1.0)*(gamma-1.0)/(U1*U1)
				  +2.0*U2/U1*(gamma-1.0)*pstar);
	dJdU[i*3+2] = delx/2.0 * (2.0*U3*(gamma-1.0)*(gamma-1.0)
				  -(gamma-1.0)*(gamma-1.0)*U2*U2/(U1)
				  -2.0*(gamma-1.0)*pstar);
    }
    
    //--------------------------------------------------------
    // Compute inv[JU]^T * -dJdQ
    //--------------------------------------------------------
    
    // We can use the inverse of JU directly
    foreach (i; 0 .. ncells) {
	double psi_rho = 0.0; double psi_ru = 0.0; double psi_rE = 0.0;
	foreach (j; 0 .. invJUT.length) {
	    if (adjoint_form == "conservative") {
		psi_rho += invJUT[i*3][j] * -dJdU[j];
		psi_ru += invJUT[i*3+1][j] * -dJdU[j];
		psi_rE += invJUT[i*3+2][j] * -dJdU[j];
	    }
	    else if (adjoint_form == "primitive") {
		psi_rho += invJVT[i*3][j] * -dJdV[j];
		psi_ru += invJVT[i*3+1][j] * -dJdV[j];
		psi_rE += invJVT[i*3+2][j] * -dJdV[j];
	    }
	    else writeln("-------------------- unknown adjoint form --------------------");
	}
	psi[i*3] = psi_rho;
	psi[i*3+1] = psi_ru;
	psi[i*3+2] = psi_rE;
    }
    //writeln(psi);

    double[3*ncells] temp;
    foreach (i; 0 .. ncells) {
	double temp1 = 0.0; double temp2 = 0.0; double temp3 = 0.0;
	foreach (j; 0 .. invJUT.length) {
	    temp1 += JVT[i*3][j] * psi[j];
	    //writeln(temp1, ", ", JVT[i*3][j], ", ", psi[j]);
	    temp2 += JVT[i*3+1][j] * psi[j];
	    temp3 += JVT[i*3+2][j] * psi[j];
	}
	temp[i*3] = temp1;
	temp[i*3+1] = temp2;
	temp[i*3+2] = temp3;
    }
    //writeln(temp);

    //writeln(invJVT);
    
    /*
    // We can also solve the system proposed by many researchers iteratively
    // {1.0/dt * [I] + [dRdU]^T} * dpsi_n = -dJdU - [dRdU]^T * psi_n
    
    static double[3*ncells][3*ncells] LHS;
    static double[3*ncells][3*ncells] invLHS;
    static double[3*ncells] psiJUT;
    static double[3*ncells] RHS;
    
    // initial psi guess
    foreach (i; 0 .. psi.length) {
    psi[i] = 1.0;
    }
    
    // construct L.H.S.
    foreach (i; 0 .. JUT.length) {
    foreach (j; 0 .. JUT.length) {
    if (i == j) LHS[i][j] = JUT[i][j] + 1.0/dt;
    else LHS[i][j] = JUT[i][j];
    }
    }
    
    // invert L.H.S.
    matrixInv(LHS, invLHS);
    
    size_t psi_step = 0;
    while (psi_step < 10000) {
    
    // Construct R.H.S.
    foreach (i; 0 .. JUT.length) {
    psiJUT[i] = 0.0;
    foreach (j; 0 .. JUT.length) {
    psiJUT[i] += JUT[i][j] * psi[j];
    }
    }
    
    foreach (i; 0 .. JUT.length) {
    RHS[i] = -dJdU[i] - psiJUT[i];
    }
    
    // perform Matrix/vector multiplication
    foreach (i; 0 .. ncells) {
    double psi_rho = 0.0; double psi_ru = 0.0; double psi_rE = 0.0;
    foreach (j; 0 .. R.length) {
    psi_rho += invLHS[i*3][j] * RHS[j];
    psi_ru += invLHS[i*3+1][j] * RHS[j];
    psi_rE += invLHS[i*3+2][j] * RHS[j];
    }
    psi[i*3] += psi_rho;
    psi[i*3+1] += psi_ru;
    psi[i*3+2] += psi_rE;	
    }
    psi_step += 1;
    }
    */
    
    //-----------------------------------------------------
    // Construct dR/dD note: dRdD = dRdA * dAdD
    //-----------------------------------------------------
    static double[3][3*ncells] dRdD;
    static double[ninterfaces][3*ncells] dRdA;
    static double[3][ninterfaces] dAdD;
    
    /*
    // Numerically constrcut dRdA using Frechet derivative
    // first compute current residual vector
    foreach(i; 0 .. ncells) {
    fvinterface fin = global_interfaces[i];
    fvinterface fout = global_interfaces[i+1];
    fvcell cell = global_cells[i+2];
    // mass flux
    R[i*3] = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
    // momentum flux
    R[i*3+1] = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
    // total energy flux
    R[i*3+2] = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
    }
    
    foreach (i; 0..ninterfaces) {
    // Perturb the areas
    double orig_area = perturbed_interfaces[i].area;
    perturbed_interfaces[i].area += perturbed_interfaces[i].area*eps+eps;
    
    // set perturbed cell volume ---------------------------
    foreach(j; 0 .. ncells) {
    global_cells[j+2].vol = 0.5*global_cells[j+2].dx*(perturbed_interfaces[j].area+perturbed_interfaces[j+1].area);
    }
    
    // now construct dR/dA using Frechet derivative
    foreach(j; 0 .. ncells) {
    double R_perturbed;
    fvinterface fin = perturbed_interfaces[j];
    fvinterface fout = perturbed_interfaces[j+1];
    fvcell cell = global_cells[j+2];
    
    // mass flux
    R_perturbed = 1.0/cell.vol * (fin.mass*fin.area - fout.mass*fout.area);
    dRdA[i][j*3] = (R_perturbed-R[j*3])/(orig_area*eps+eps);
    
    // momentum flux
    R_perturbed = 1.0/cell.vol * (fin.momentum*fin.area - fout.momentum*fout.area + cell.p*(fout.area-fin.area));
    dRdA[i][j*3+1] = (R_perturbed-R[j*3+1])/(orig_area*eps+eps);
    
    // total energy flux
    R_perturbed = 1.0/cell.vol * (fin.energy*fin.area - fout.energy*fout.area);
    dRdA[i][j*3+2] = (R_perturbed-R[j*3+2])/(orig_area*eps+eps);
    }
    perturbed_interfaces[i].area = orig_area;
    }
    */
    
    /*
    // Analytically construct dAdD by hand differentiation
    foreach (i; 0..ninterfaces) {
    dAdD[i] = global_interfaces[i].xpos;
    }
    */
    /*
    // Numerically construct dAdD by Frechet derivatives
    //perturb b
    double b_orig = b;
    b += b*eps + eps;
    a = yo - b*tanh(-d);
    foreach(i; 0 .. ninterfaces) {
    double radius = a + b*tanh(c*global_interfaces[i].xpos -d);
    perturbed_interfaces[i].area = PI*radius*radius;
    }
     foreach(i; 0..ninterfaces) {
    dAdD[i][0] = (perturbed_interfaces[i].area - global_interfaces[i].area)/(b*eps + eps);
    }
    b = b_orig;
    a = yo - b*tanh(-d);
    */
    // Analytically construct dRdA by hand differentiation
    foreach (i; 0..ncells) {
	foreach(j; 0 .. ninterfaces) {
	    fvinterface fin = global_interfaces[i];
	    fvinterface fout = global_interfaces[i+1];
	    fvcell cell = global_cells[i+2];
	    
	    // mass flux
	    if (i == j) dRdA[i*3][j] = -1.0* -fin.mass/cell.vol;
	    else if (j == i+1) dRdA[i*3][j] = -1.0* fout.mass/cell.vol;
	    else dRdA[i*3][j] = 0.0;
	    
	    // momentum flux
	    if (i == j) dRdA[i*3+1][j] = -1.0* -(fin.momentum-cell.p)/cell.vol;
	    else if (j == i+1) dRdA[i*3+1][j] = -1.0* -(-fout.momentum+cell.p)/cell.vol;
	    else dRdA[i*3+1][j] = 0.0;
	    
	    // total energy flux
	    if (i == j) dRdA[i*3+2][j] = -1.0* -fin.energy/cell.vol;
	    else if (j == i+1) dRdA[i*3+2][j] = -1.0* fout.energy/cell.vol;
	    else dRdA[i*3+2][j] = 0.0;
	}
    }
    foreach(i; 0..ninterfaces) {
	fvinterface F = global_interfaces[i];
	dAdD[i][0] = 2.0*tanh(d/scale) + 2.0*tanh((c*F.xpos - d)/scale);
	dAdD[i][1] = 2.0*b*F.xpos*(-pow(tanh((c*F.xpos - d)/scale),2) + 1)/scale;
	dAdD[i][2] = 2.0*b*(-pow(tanh(d/scale),2) + 1)/scale + 2.0*b*(pow(tanh((c*F.xpos - d)/scale),2) - 1)/scale;
    }
    // dRdD = dRdA*dAdD
    foreach(i; 0..ncells*3 ) {
	dRdD[i][0] = 0.0;
	dRdD[i][1] = 0.0;
	dRdD[i][2] = 0.0;
	foreach(j; 0..ninterfaces ) {
	    dRdD[i][0] += dRdA[i][j]*dAdD[j][0];
	    dRdD[i][1] += dRdA[i][j]*dAdD[j][1];
	    dRdD[i][2] += dRdA[i][j]*dAdD[j][2];
	}
    }
    
    //-----------------------------------------------------
    // Compute dLdD = dRdD * psi
    //-----------------------------------------------------
    dLdB = 0.0;
    dLdC = 0.0;
    dLdD = 0.0;
    foreach (j; 0 .. psi.length) {
	dLdB += dRdD[j][0] * psi[j];
	dLdC += dRdD[j][1] * psi[j];
	dLdD += dRdD[j][2] * psi[j];
    }
}

//-----------------------------------------------------
// Main Body
//-----------------------------------------------------

void main() {
    //-----------------------------------------------------
    // Select Simulation type
    //-----------------------------------------------------
    string simulation_type = "nozzle"; // options: sod, nozzle
    
    //-----------------------------------------------------
    // Construct Mesh
    //-----------------------------------------------------

    fvcell[ncells+nghost] global_cells;
    double dx = length/ncells; // cell width

    // construct ghost cells ---------------------

    global_cells[0] = new fvcell(10000001, -1.5*dx, dx);
    global_cells[1] = new fvcell(10000002, -0.5*dx, dx);
    global_cells[ncells+2] = new fvcell(10000003, length+0.5*dx, dx);
    global_cells[ncells+3] = new fvcell(10000004, length+1.5*dx, dx);

    // construct inner domain cells --------------------

    double temp_pos = 0.5*dx;
    foreach(i; 0 .. ncells) {
	global_cells[i+2] = new fvcell(i, temp_pos, dx);
	temp_pos += dx ;
    }

    // construct interfaces ----------------------

    fvinterface[ninterfaces] global_interfaces;
    fvinterface[ninterfaces] perturbed_interfaces; // for implicit solver
    temp_pos = 0.0;
    foreach(i; 0 .. ninterfaces) {
	global_interfaces[i] = new fvinterface(i, temp_pos);
	perturbed_interfaces[i] = new fvinterface(i, temp_pos);
	temp_pos += dx;
    }

    // set interface areas -----------------------
    double inlet_area; double exit_area;
    double yo; double a; double b; double c; double d; double scale;
    if (simulation_type == "sod") {
	// shock tube geometry (no area variaton)
	foreach(i; 0 .. ninterfaces) {
	    global_interfaces[i].area = 1.0;
	    perturbed_interfaces[i].area = 1.0;
	}
    }
    else if (simulation_type == "nozzle") {
	// linear nozzle, with design variable a
	// target double b = 0.347, c = 0.8, d = 4.0, yo = 1.05
	//yo = 1.05;
	//b = 0.32, c = 1.0, d = 3.8;
	// current target: b = 0.05, c = 0.85, d = 4.0;
	b = 0.07, c = 0.8, d = 3.8;
	//b = 0.05, c = 0.85, d = 4.0;
	yo = 0.105;
	scale = 1.5;
	a = yo - b*tanh(-d/scale);
	foreach(i; 0 .. ninterfaces) {
	    double radius = a + b*tanh((c*global_interfaces[i].xpos -d)/scale);
	    global_interfaces[i].area =  2.0*radius;
	    perturbed_interfaces[i].area = 2.0*radius;
	}
    }
    // set cell volume ---------------------------
    
    foreach(i; 0 .. ncells) {
	global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
    }

    //-----------------------------------------------------
    // Set air properties
    //-----------------------------------------------------
    
    double gamma = 1.4;
   
    //-----------------------------------------------------
    // Initial Condiions
    //-----------------------------------------------------
    if (simulation_type == "sod" ) {
	//-----------------------------------------------------
	// Initial Conditions -- Sod's Shocktube
	//-----------------------------------------------------
	// left fill state
	double p_init = 1.0e05; // Pa
	double rho_init = 1.0; // kg/m^3
	double u_init = 0.0; // m/s
	foreach(i; 0 .. to!int(0.5*ncells)) {
	    global_cells[i+2].p = p_init;
	    global_cells[i+2].rho = rho_init;
	    global_cells[i+2].u = u_init;
	}
	
	// right fill state
	p_init = 1.0e04; // Pa
	rho_init = 0.125; // kg/m^3
	u_init = 0.0; // m/s
	foreach(i; to!int(0.5*ncells) .. ncells) {
	    global_cells[i+2].p = p_init;
	    global_cells[i+2].rho = rho_init;
	    global_cells[i+2].u = u_init;
	}
    }
    else if (simulation_type == "nozzle") {
	//-----------------------------------------------------
	// Initial and Inflow Conditions -- Nozzle
	//-----------------------------------------------------

	// use the analytical solution as initial guess (based on inflow conditions) 
	/*
	double p_init; double rho_init; double u_init;
      	foreach(i; 0 .. ncells) {
	    double cA = 0.5*(global_interfaces[i].area+global_interfaces[i+1].area);
	    double AStar = 5.925925926e-03;
	    double rhoStar = 3.19933979;
	    double pStar = 413351.3941;
	    double tolerance = 0.001;
	    double delta = 100.0;
	    double M = 2.0;  // Mach number
	    double dM = 0.00001;
	    while (abs(delta) > tolerance) {
		M += dM;
		delta = (cA/AStar)*(cA/AStar) -
		    1.0/(M*M)*pow(( (2.0/(gamma+1)) * (1.0 + (gamma-1.0)/2.0 * (M*M))), (gamma+1.0)/(gamma-1.0));
	    }
	    p_init = pStar *
		pow((1.0 + (gamma-1.0)/2.0 * M*M), -gamma/(gamma-1.0)) / pow((1.0 + (gamma-1.0)/2.0), -gamma/(gamma-1.0));
	    rho_init = rhoStar *
		pow((1.0 + (gamma-1.0)/2.0 * M*M), -1.0/(gamma-1.0)) / pow((1.0 + (gamma-1.0)/2.0), -1.0/(gamma-1.0));
	    u_init = M * sqrt((p_init*gamma)/rho_init);
	    // overriding the anlytical condition with the inflow conditions
	    global_cells[i+2].p = p_init;
	    global_cells[i+2].rho = rho_init;
	    global_cells[i+2].u = u_init;
	    }
	*/

	// set inflow conditions
	double Mach = 1.5;
	double p_inflow = 101.325e3; // Pa
	double rho_inflow = 1.0;   // kg/m^3
	double a_inflow = sqrt((p_inflow*gamma)/rho_inflow);
	double u_inflow = Mach * a_inflow;  // m/s

	foreach(i; 0 .. ncells) {
	    global_cells[i+2].p = p_inflow;
	    global_cells[i+2].rho = rho_inflow;
	    global_cells[i+2].u = 0.6 * a_inflow;
	}

	// Nozzle inflow state ----------------------------------------------
	global_cells[0].p = p_inflow; global_cells[1].p = p_inflow;
	global_cells[0].rho = rho_inflow; global_cells[1].rho = rho_inflow;
	global_cells[0].u = u_inflow; global_cells[1].u = u_inflow;

	global_cells[ncells+2].p = global_cells[ncells+1].p; global_cells[ncells+3].p = global_cells[ncells+1].p;
	global_cells[ncells+2].rho =  global_cells[ncells+1].rho; global_cells[ncells+3].rho =  global_cells[ncells+1].rho;
	global_cells[ncells+2].u =  global_cells[ncells+1].u; global_cells[ncells+3].u =  global_cells[ncells+1].u;	
    }
    // ---------------------------------------------------------------------------------- //
    // Now we're ready to perform some finite volume calculations on the cells
    // ---------------------------------------------------------------------------------- //

    //-----------------------------------------------------
    // Set some more  simulations parameters
    //-----------------------------------------------------
    double dt = 1.0e-06;                                // time step size, s
    double sim_time = 0.0;                              // current simulation time, s
    double final_time;                                  // final simulation time, s
    if (simulation_type == "sod") final_time = 6.0e-04; 
    else if (simulation_type == "nozzle") final_time = 2.0; 
    string flux_calc = "ausmdv"; // van_leer, ausmdv, ausm_plus_up
    size_t step = 0;
    size_t opt_iter = 0;
    size_t max_step = 100000000;
    size_t interpolation_order = 1;
    size_t outflow_bc_order = 1;
    string solver = "verification"; // options: simulation, optimisation, verification

    double dLdB_old = 0.0; double dLdC_old = 0.0; double dLdD_old = 0.0;
    double dLdB = 1e10; double dLdC = 1e10; double dLdD = 1e10;
    double J = 0.0;
    double Jref = 0.0;
    // implicit solver settings and arrays  ------------------------
    string adjoint_form = "primitive";            // options: conservative form, primtive form
    static double[3*ncells][3*ncells] JV;        // flow Jacobian w.r.t. primitive variables
    static double[3*ncells][3*ncells] JU;        // flow Jacobian w.r.t. conserved variables
    static double[3*ncells][3*ncells] transform; // transform from JV to JU
    static double[3*ncells] R;                   // R.H.S. resiudal vector
    
    // adjoint solver arrays ---------------------------------------
    static double[3*ncells] dJdV;                // sensitivty of cost function (J)
                                                 // w.r.t. primitive flow varaibles (V)
    static double[3*ncells] dJdU;                // sensitivty of cost function (J)
                                                 // w.r.t. conserved flow varaibles (U)
    static double[3*ncells][3*ncells] JUT;       // JU transpose
    static double[3*ncells][3*ncells] invJUT;    // JU transpose inverse
    static double[3*ncells][3*ncells] JVT;       // JU transpose
    static double[3*ncells][3*ncells] invJVT;    // JU transpose inverse
    static double[3*ncells] psi;                 // adjoint variables
	
    double eps = 1.0e-8;                         // finite difference epsilon
    double tol = 1.0e-05;

    double[3] pk;
    double[3] pk_old;
    if (solver == "simulation") {
	solve(global_cells, global_interfaces, perturbed_interfaces,
	      dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
	      JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD, dx,
	      gamma, b, c, d, yo, scale, solver, adjoint_form);
    }
    else if(solver == "verification") {
	writeln("-------------------------------------");
	writeln("Compute gradient via adjoint method");
	writeln("-------------------------------------");
	solve(global_cells, global_interfaces, perturbed_interfaces,
	      dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
	      JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD, dx,
	      gamma, b, c, d, yo, scale, solver, adjoint_form);
	double adjoint_dLdB = dLdB; double adjoint_dLdC = dLdC; double adjoint_dLdD = dLdD;
	writeln("-------------------------------------");
	writeln("Compute gradient via finite difference method");
	writeln("-------------------------------------");
	// perturb b
	writeln("-- perturb b --");
	double fd_dLdB;
	double pertb = b + 1.0e-06;
	// update interface areas ---------------------------
	a = yo - pertb*tanh(-d/scale);
	foreach(i; 0 .. ninterfaces) {
	    double radius = a + pertb*tanh((c*global_interfaces[i].xpos -d)/scale);
	    global_interfaces[i].area =  2.0*radius;
 	    perturbed_interfaces[i].area = 2.0*radius;
	}
	// update cell volumes ---------------------------
	foreach(i; 0 .. ncells) {
	    global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
	}
	double Jminus; double Jplus;
	solve(global_cells, global_interfaces, perturbed_interfaces,
	      dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
	      JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD, dx,
	      gamma, pertb, c, d, yo, scale, solver, adjoint_form);
	Jplus = J;
        pertb = b - 1.0e-06;
	// update interface areas ---------------------------
	a = yo - pertb*tanh(-d/scale);
	foreach(i; 0 .. ninterfaces) {
	    double radius = a + pertb*tanh((c*global_interfaces[i].xpos -d)/scale);
	    global_interfaces[i].area =  2.0*radius;
	    perturbed_interfaces[i].area = 2.0*radius;
	}
	// update cell volumes ---------------------------
	foreach(i; 0 .. ncells) {
	    global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
	}
	solve(global_cells, global_interfaces, perturbed_interfaces,
	      dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
	      JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD, dx,
	      gamma, pertb, c, d, yo, scale, solver, adjoint_form);
	Jminus = J;
	fd_dLdB = (Jplus - Jminus)/(2.0 * 1.0e-06);
	// perturb c
	writeln("-- perturb c --");
	double fd_dLdC;
	double pertc = c + 1.0e-06;
	// update interface areas ---------------------------
	a = yo - b*tanh(-d/scale);
	foreach(i; 0 .. ninterfaces) {
	    double radius = a + b*tanh((pertc*global_interfaces[i].xpos -d)/scale);
	    global_interfaces[i].area =  2.0*radius;
	    perturbed_interfaces[i].area = 2.0*radius;
	}
	// update cell volumes ---------------------------
	foreach(i; 0 .. ncells) {
	    global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
	}
	solve(global_cells, global_interfaces, perturbed_interfaces,
	      dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
	      JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD, dx,
	      gamma, b, pertc, d, yo, scale, solver, adjoint_form);
	Jplus = J;
	pertc = c - 1.0e-06;
	// update interface areas ---------------------------
	a = yo - b*tanh(-d/scale);
	foreach(i; 0 .. ninterfaces) {
	    double radius = a + b*tanh((pertc*global_interfaces[i].xpos -d)/scale);
	    global_interfaces[i].area =  2.0*radius;
	    perturbed_interfaces[i].area = 2.0*radius;
	}
	// update cell volumes ---------------------------
	foreach(i; 0 .. ncells) {
	    global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
	}
	solve(global_cells, global_interfaces, perturbed_interfaces,
	      dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
	      JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD, dx,
	      gamma, b, pertc, d, yo, scale, solver, adjoint_form);
	Jminus = J;
	fd_dLdC = (Jplus - Jminus)/(2.0 * 1.0e-06);
	// perturb d
	writeln("-- perturb d --");
	double fd_dLdD;
	double pertd = d + 1.0e-06;
	// update interface areas ---------------------------
	a = yo - b*tanh(-pertd/scale);
	foreach(i; 0 .. ninterfaces) {
	    double radius = a + b*tanh((c*global_interfaces[i].xpos -pertd)/scale);
	    global_interfaces[i].area =  2.0*radius;
	    perturbed_interfaces[i].area = 2.0*radius;
	}
	// update cell volumes ---------------------------
	foreach(i; 0 .. ncells) {
	    global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
	}
	solve(global_cells, global_interfaces, perturbed_interfaces,
	      dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
	      JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD,
	      dx, gamma, b, c, pertd, yo, scale, solver, adjoint_form);
	Jplus = J;
	pertd = d - 1.0e-06;
	// update interface areas ---------------------------
	a = yo - b*tanh(-pertd/scale);
	foreach(i; 0 .. ninterfaces) {
	    double radius = a + b*tanh((c*global_interfaces[i].xpos -pertd)/scale);
	    global_interfaces[i].area =  2.0*radius;
	    perturbed_interfaces[i].area = 2.0*radius;
	}
	// update cell volumes ---------------------------
	foreach(i; 0 .. ncells) {
	    global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
	}
	solve(global_cells, global_interfaces, perturbed_interfaces,
	      dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
	      JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD,
	      dx, gamma, b, c, pertd, yo, scale, solver, adjoint_form);
	Jminus = J;
	fd_dLdD = (Jplus - Jminus)/(2.0 * 1.0e-06);
	writeln("-------------------------------------");
	writeln("Results");
	writeln("-------------------------------------");
	writeln("fd_b = ", fd_dLdB);
	writeln("adjoint_b = ", adjoint_dLdB);
	writeln("% error b = ", (adjoint_dLdB-fd_dLdB)/adjoint_dLdB * 100.0);
	writeln("fd_c = ", fd_dLdC);
	writeln("adjoint_c = ", adjoint_dLdC);
	writeln("% error c = ", (adjoint_dLdC-fd_dLdC)/adjoint_dLdC * 100.0);
	writeln("fd_d = ", fd_dLdD);
	writeln("adjoint_d = ", adjoint_dLdD);
	writeln("% error d = ", (adjoint_dLdD-fd_dLdD)/adjoint_dLdD * 100.0);
    }
    else if (solver == "optimisation") {
	writeln(b, ", ", c, ", ", d);
	while (opt_iter < 1e6) {
	    //-----------------------------------------------------
	    // conjugate gradient update
	    //-----------------------------------------------------
	    	    	    
	    // 1. compute flow solution and gradients for current geometry
	    solve(global_cells, global_interfaces, perturbed_interfaces,
		  dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
		  JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD,
		  dx, gamma, b, c, d, yo, scale, solver, adjoint_form);

	    // if on the first step, then set the reference J
	    if (opt_iter == 0) Jref = J;

	    // if the gradients and cost function are small enough then we are at a minimum
	    if (J/Jref < tol) break; // && abs(dLdB) < 1.0 && abs(dLdC) < 1.0 && abs(dLdD) < 1.0) break;
	    
	    // 2. compute search direction
	    double beta;
	    if (dLdB_old == 0.0 && dLdC_old == 0.0 && dLdD_old == 0.0) {
		double norm = sqrt(dLdB*dLdB + dLdC*dLdC + dLdD*dLdD);
		pk[0] = -dLdB/norm;
		pk[1] = -dLdC/norm;
		pk[2] = -dLdD/norm;
	    }
	    else if( (opt_iter/70.0 - trunc(opt_iter/70.0)) == 0.0) {
		double norm = sqrt(dLdB*dLdB + dLdC*dLdC + dLdD*dLdD);
		pk[0] = -dLdB/norm;
		pk[1] = -dLdC/norm;
		pk[2] = -dLdD/norm;
	    }
	    else {
		double norm = sqrt(dLdB*dLdB + dLdC*dLdC + dLdD*dLdD);
		double beta_numer = dLdB*dLdB + dLdC*dLdC + dLdD*dLdD;
		double beta_denom = dLdB_old*dLdB_old + dLdC_old*dLdC_old + dLdD_old*dLdD_old;
		beta = beta_numer/beta_denom;
		pk[0] = -dLdB/norm + beta*pk_old[0]; 
		pk[1] = -dLdC/norm + beta*pk_old[1]; 
		pk[2] = -dLdD/norm + beta*pk_old[2]; 
	    }
	    // 3. perform line search to find step size ak
	    double rho = 0.5;
	    double mu = 1.0e-4;
	    double ak = 1.0;
	    double bb = b;
	    double cc = c;
	    double dd = d;
	    double Jk = J;
	    double dLdBB = dLdB;
	    double dLdCC = dLdC;
	    double dLdDD = dLdD;
	    
	    b = b + ak*pk[0];
	    c = c + ak*pk[1];
	    d = d + ak*pk[2];

	    // update interface areas ---------------------------
	    a = yo - b*tanh(-d/scale);
	    foreach(i; 0 .. ninterfaces) {
		double radius = a + b*tanh((c*global_interfaces[i].xpos -d)/scale);
		global_interfaces[i].area =  2.0*radius;
		perturbed_interfaces[i].area = 2.0*radius;
	    }
	    // update cell volumes ---------------------------
	    foreach(i; 0 .. ncells) {
		global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
	    }
	    
	    solve(global_cells, global_interfaces, perturbed_interfaces,
		  dt, sim_time, final_time, step, max_step, flux_calc, simulation_type, interpolation_order, outflow_bc_order,
		  JV, JU, transform, R, dJdV, dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD,
		  dx, gamma, b, c, d, yo, scale, solver, adjoint_form);

	    double Jk1 = J;
	    while (Jk1 > Jk + mu*ak*(dLdBB*pk[0]+dLdCC*pk[1]+dLdDD*pk[2])) {
		    ak = ak*rho;
		    b = bb + ak*pk[0];
		    c = cc + ak*pk[1];
		    d = dd + ak*pk[2];

		    // update interface areas ---------------------------
		    a = yo - b*tanh(-d/scale);
		    foreach(i; 0 .. ninterfaces) {
			double radius = a + b*tanh((c*global_interfaces[i].xpos -d)/scale);
			global_interfaces[i].area =  2.0*radius;
			perturbed_interfaces[i].area = 2.0*radius;
		    }
		    // update cell volumes ---------------------------
		    foreach(i; 0 .. ncells) {
			global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
		    }
		    
		    solve(global_cells, global_interfaces, perturbed_interfaces,
			  dt, sim_time, final_time, step, max_step, flux_calc, simulation_type,
			  interpolation_order, outflow_bc_order, JV, JU, transform, R, dJdV,
			  dJdU, JUT, invJUT, JVT, invJVT, psi, eps, J, dLdB, dLdC, dLdD, dx,
			  gamma, b, c, d, yo, scale, solver, adjoint_form);
		    Jk1 = J;
		    //writeln("ak = ", ak, ", J_k = ", Jk, ", J_k+1 = ", Jk1,  ", b = ", b, " c = ", c, ", d = ", d);
	    }

	    // update interface areas ---------------------------
	    a = yo - b*tanh(-d/scale);
	    foreach(i; 0 .. ninterfaces) {
		double radius = a + b*tanh((c*global_interfaces[i].xpos -d)/scale);
		global_interfaces[i].area =  2.0*radius;
		perturbed_interfaces[i].area = 2.0*radius;
	    }
	    // update cell volumes ---------------------------
	    foreach(i; 0 .. ncells) {
		global_cells[i+2].vol = 0.5*global_cells[i+2].dx*(global_interfaces[i].area+global_interfaces[i+1].area);
	    }

	    dLdB_old = dLdB; 
	    dLdC_old = dLdC;
	    dLdD_old = dLdD; 
	    pk_old[0] = pk[0];
	    pk_old[1] = pk[1];
	    pk_old[2] = pk[2];
	    opt_iter += 1;
	    writeln("--------------------------------------------------------------------------");
	    writeln(opt_iter, ". J/Jref = ", J/Jref, ", b = ", b, ", c = ", c, ", d = ", d);
	    writeln("gradients: ", ", dLdB = ", dLdB,  ", dLdC = ", dLdC,  ", dLdD = ", dLdD);
	    writeln("--------------------------------------------------------------------------");
	}
    writeln(b, ", ", c, ", ", d);
    }

    writeln("--------------------------------------------------------------------------");
    writeln("FINAL RESULTS:", " b = ", b, ", c = ", c, ", d = ", d);
    writeln("gradients: ", ", dLdB = ", dLdB,  ", dLdC = ", dLdC,  ", dLdD = ", dLdD);
    writeln("--------------------------------------------------------------------------");

    //--------------------------------------------------------
    // Write properties out for plotting at end of simulation
    //--------------------------------------------------------
    
    foreach(i; 0 .. ncells) {
	fvcell cell = global_cells[i+2];
	auto writer = format("%f %f %f %f %f %f %f \n", cell.xpos, cell.rho, cell.p, cell.u, psi[i*3], psi[i*3+1], psi[i*3+2]);
	append("output.dat", writer);
    }
    foreach(i; 0 .. ninterfaces) {
	fvinterface f = global_interfaces[i];
	auto writer = format("%f %f \n", f.xpos, f.area/2.0);
	append("contour.dat", writer);
    }   
}

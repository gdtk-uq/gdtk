/**
 * solid_loose_coupling_update.d
 *
 * Author: Nigel Chong and Rowan G.
 * Date: 2017-04-18
 */

module solid_loose_coupling_update;

import globalconfig;
import globaldata;
import solidblock; 
import solidfvcell;  
import std.stdio; 
import std.parallelism; 
import nm.bbla;
import std.math;
import solidprops; 
import std.datetime; 
import solid_udf_source_terms;


Matrix eval_dedts(Matrix eip1, int ftl, double sim_time)
// Evaulates the energy derivatives of all the cells in the solid block and puts them into an array (Matrix object)
// Evaluates derivative when cell energy is "eip1" 
 
// Returns: Matrix of doubles
{ 
    auto n = eip1.nrows;
    Matrix ret = zeros(n, 1); // array that will be returned (empty - to be populated)
    foreach (sblk; parallel(solidBlocks, 1)) {
	if (!sblk.active) continue;
	sblk.applyPreSpatialDerivAction(sim_time, ftl);
	sblk.clearSources(); 
	sblk.computeSpatialDerivatives(ftl); 
	sblk.computeFluxes();
    }
    if (GlobalConfig.apply_bcs_in_parallel) {
	foreach (sblk; parallel(solidBlocks, 1)) {
	    if (sblk.active) { sblk.applyPostFluxAction(sim_time, ftl); }
	}
    } else {
	foreach (sblk; solidBlocks) {
	    if (sblk.active) { sblk.applyPostFluxAction(sim_time, ftl); }
	}
    }
    int i = 0;
    foreach (sblk; parallel(solidBlocks, 1)) {
	foreach (scell; sblk.activeCells) { 
	    if (GlobalConfig.udfSolidSourceTerms) {
		addUDFSourceTermsToSolidCell(sblk.myL, scell, sim_time);
	    }
	    scell.e[ftl] = eip1[i, 0];
	    scell.T = updateTemperature(sblk.sp, scell.e[ftl]);
	    scell.timeDerivatives(ftl, GlobalConfig.dimensions);
	    ret[i, 0] = scell.dedt[ftl];
	    i += 1;
	}
    }
    return ret;
}


Matrix e_vec(int ftl)
// Puts energies from desired ftl into an array (Matrix object) 

// Returns: Matrix of doubles
{
    double[] test1;
    foreach (sblk; solidBlocks){ 
	foreach(scell; sblk.activeCells){
	    test1 ~= scell.e[ftl];
	}
    } 
    auto len = test1.length; 
    Matrix ret = zeros(len, 1); 
    int count = 0;
    foreach (entry; test1){ 
	ret[count, 0] = entry; 
	count += 1;
    }
    return ret;
} 

Matrix F(Matrix eip1, Matrix ei, double dt_global, double sim_time) 
// Function for which root needs to be found (for Newton's method)
// Re-arranged implicit euler equation 

// Returns: Matrix of doubles
{
    return eip1 - ei - dt_global*eval_dedts(eip1, 1, sim_time); // lets get rid of dedt_vec... 
}

  
Matrix copy_vec(Matrix mat)  
// Copies an original vector to a new vector (i.e. leaves original untouched) 

// Returns: Matrix of doubles
{
    auto n = mat.nrows; 
    auto cop = zeros(n, 1); 
    foreach(i; 0 .. n){
	cop[i, 0] = mat[i, 0];
    }
    return cop; 
} 

double norm(Matrix vec)
// Calculates the euclidean norm of a vector 

// Returns: double
{
    auto n = vec.nrows;
    double total = 0; 
    foreach(i;0 .. n){	    
	total += pow(vec[i,0],2);
    }
    double ret = pow(total, 0.5);
    return ret; 
} 

Matrix Jac(Matrix Yip1, Matrix Yi, double dt, double sim_time, double eps = 1e-2)
// Jacobian of F 

// Returns: Matrix of doubles
{ 
    auto n = Yip1.nrows;
    auto ret = eye(n);
    Matrix Fout;
    Matrix Fin = F(Yip1, Yi, dt, sim_time);
    Matrix cvec;
    foreach(i;0 .. n){ 
	cvec = copy_vec(Yip1);
	cvec[i, 0] += eps; 
	Fout = F(cvec, Yi, dt, sim_time); 
	foreach(j; 0 .. n){
	    ret[j, i] = (Fout[j, 0] - Fin[j, 0])/(eps);		      
	}
    }
    return ret;
}


Matrix GMRES(Matrix A, Matrix b, Matrix x0)
// Executes GMRES_step, can be used to repeat GMRES steps checking for multiple iterations, otherwise relatively useless 

// Returns: Matrix of doubles
{    
    auto m = 10;
    auto xm = GMRES_step(A, b, x0, m);
    
    return xm;
}


Matrix GMRES_step(Matrix A, Matrix b, Matrix x0, int m)
// Takes a single GMRES step using Arnoldi 
// Minimisation problem is solved using Least Squares via Gauss Jordan (see function below) 

// Returns: Matrix of doubles
{ 
    auto r0 = b - dot(A, x0);
    double beta = norm(r0);  
    auto v = vstack([r0])/beta; 
    auto H = zeros(m+1, m); 
    auto e1 = zeros(m+1, 1); 
    e1[0, 0] = 1; 
    auto b1 = beta*e1;
    Matrix vi; 
    Matrix wj; 
    Matrix temp;
    double hjp1j;
    Matrix H1;
    Matrix b2; 
    Matrix y; 
    double hij;
    double tol = 1e-3;
    foreach(j; 0 .. m){
	if(j==0){
	    wj = dot(A, v);
	    hij = dot(flatten(wj), v)[0, 0]; 
	    H[0, 0] = hij; 
	    wj = wj - hij*v;  
	    wj = flatten(wj);
	}
	else{
	    wj = flatten(dot(A, v.sliceDup(0, v.nrows, j, j+1)));
	    foreach(i; 0 .. j){ 
		vi = v.sliceDup(0, v.nrows, i, i+1); 
		hij = dot(wj, vi)[0, 0];
		H[i, j] = hij;  
		wj = vertit(wj);
		wj = wj - hij*vi;   
		wj = flatten(wj);
	    } 
	}
	wj = vertit(wj);
	hjp1j = norm(wj);
	H[j+1, j] = hjp1j;    
	if(fabs(hjp1j) < tol){   
	    H1 = H.sliceDup(0, j+1, 0, j+1); 
	    b2 = b1.sliceDup(0, j+1, 0, b1.ncols);  
	    return x0 + dot(v, lsq_gmres(H1, b2)); 
	}
	else{
	    if(j==0){
		v = hstack([v, wj/hjp1j]);
	    }
	    else{
		v = hstack([v, wj/hjp1j]); 
	    }
	}
    }
    return x0 + dot(v.sliceDup(0, v.nrows, 0, v.ncols-1), lsq_gmres(H, b1)); 
}

Matrix flatten(Matrix vec)
// E.g. Takes Matrix[[1], [2], [3]] and returns Matrix[[1, 2, 3]]

// Returns: Matrix of doubles
{
    auto len = vec.nrows;
    Matrix ret = zeros(1, len);
    foreach(i; 0 .. len){
	ret[0, i] = vec[i, 0];
    }
    return ret;
} 

Matrix vertit(Matrix vec)
// E.g. Takes Matrix[[1, 2, 3]] and returns Matrix [[1], [2], [3]]

// Returns: Matrix of doubles
{
    auto len = vec.ncols; 
    Matrix ret = zeros(len, 1); 
    foreach(i; 0 .. len){
	ret[i, 0] = vec[0, i];
    }
    return ret; 
}
    

Matrix lsq_gmres(Matrix A, Matrix b)
// Solves the minimisation problem using least squares and direct solution (Gauss Jordan elimination)

// Returns: Matrix of doubles
{
    auto Ht = transpose(A); 
    auto c = hstack([dot(Ht, A), dot(Ht, b)]);
    gaussJordanElimination(c);  
    auto nr = c.nrows;
    return arr_2_vert_vec(c.getColumn(nr));
}

Matrix arr_2_vert_vec(double[] arr) 
// Takes an array (i.e. D equivalent of a python list) and converts it into a vertical Matrix 
// E.g. [1, 2, 3] --> Matrix[[1], [2], [3]]

// Returns: Matrix of doubles
{
    auto len = arr.length;
    auto ret = zeros(len, 1);
    foreach(i; 0 .. len){
	ret[i, 0] = arr[i];
    }
    return ret; 
}


void post(Matrix eip1, Matrix dei)
// Does the post processing of an implicit method step 

// Returns: none
{
    int i = 0;
    foreach (sblk; solidBlocks){ 
	foreach(scell; sblk.activeCells){
	    double eicell = eip1[i, 0];
	    scell.e[0] = eicell;
	    scell.de_prev = dei[i, 0];
	    scell.T = updateTemperature(sblk.sp, eicell); //not necesary? 
	    i += 1; 
	}
    }
}

void solid_domains_backward_euler_update(double sim_time, double dt_global) 
// Executes implicit method
{  
    writeln("== Begin Step ==");
    auto t0 = Clock.currTime();
    int n = 0;
    foreach (sblk; parallel(solidBlocks, 1)) {
	foreach (scell; sblk.activeCells) {
	    if (sim_time-dt_global == 0){ 
		// initialise e0 values
		scell.e[0] = updateEnergy(sblk.sp, scell.T); 
	    }
	    // make a guess for e[1]
	    scell.e[1] = scell.e[0] + scell.de_prev; 
	    n += 1;
	}
    }
    double eps = 1e-2;
    double omega = 0.01;
    
    auto ei = e_vec(0); 
    // retrieving guess
    auto eip1 = e_vec(0);
    // Begin prepping for newton loop
    auto Xk = copy_vec(eip1);
    auto Xkp1 = Xk + eps;
    auto Fk = F(Xk, ei, dt_global, sim_time); 
    auto Fkp1 = F(Xkp1, ei, dt_global, sim_time); // done so that fabs(normfk - normfkp1) satisfied
    int count = 0;
    Matrix Jk;
    Matrix mFk; // e.g. minusFk --> mFk
    Matrix dX;   
    while(fabs(norm(Fk) - norm(Fkp1)) > eps){ 
	count += 1; 
	if (count != 1){
	    Xk = Xkp1; 
	    Fk = Fkp1;  
	}
	Jk = Jac(Xk, ei, dt_global, sim_time);
	mFk = -1*Fk;
	dX = GMRES(Jk, mFk, Xk);
	
	Xkp1 = Xk + dX; 
	Fkp1 = F(Xkp1, ei, dt_global, sim_time); 
	
	// Break out if stuck or taking too many computations
	if (count == 10){
	    break;
	}

    }						
    eip1 = Xkp1;
    Matrix dei = eip1 - ei;
    Fk = Fkp1;
    post(eip1, dei); // post processing incl. writing solns for e, de and updating T.
    writeln("== End Step ==");
}


// main.d: a steady laminar flamelet calculator by NNG.

module slf;

import std.stdio;
import std.math;
import std.mathspecial;
import std.format;
import std.algorithm;
import std.string;
import gas;
import gas.physical_constants;
import kinetics;
import nm.smla;
import nm.bbla;
import nm.number;
import nm.complex;

/* Code borrowed from zero_rk:
   Approximate inverse error function implemented from:
   "A handy approximation for the error function and its inverse"
   by Sergei Winitzki
*/

@nogc
double erfc_inv(double q) {
    if(q <= 0.0) return float.infinity;
    if(q >= 2.0) return -float.infinity;
    double p = 1 - q;
    double sgn = p < 0 ? -1 : 1;

    double logpt = log((1 - p)*(1 + p));
    double pt1 = 2/(PI*0.147) + 0.5*logpt;
    double pt2 = 1/(0.147) * logpt;

    return(sgn*sqrt(-pt1 + sqrt(pt1*pt1 - pt2)));
}

@nogc
void second_derivs_from_cent_diffs(ref const Parameters pm, number[] U, number[] U2nd){
/*
    Use a central difference stencil to get second order derivatives,
    assuming fixed value end conditions

*/
    size_t N = pm.N;
    size_t neq = pm.neq;

    size_t lft, ctr, rgt;
    ctr = 0*neq;
    rgt = 1*neq;
    foreach(ii; 0 .. neq) U2nd[ctr + ii] = (U[rgt + ii] - 2.0*U[ctr + ii] + pm.U0[ii])/pm.dZ/pm.dZ;

    foreach(i; 1 .. N-1) {
        lft = (i-1)*neq;
        ctr = i*neq;
        rgt = (i+1)*neq;

        foreach(ii; 0 .. neq) U2nd[ctr+ii] = (U[rgt + ii] - 2.0*U[ctr+ii] + U[lft + ii])/pm.dZ/pm.dZ;
    }

    
    lft = (N-2)*neq;
    ctr = (N-1)*neq;
    foreach(ii; 0 .. neq) U2nd[ctr + ii] = (pm.U1[ii] - 2.0*U[ctr + ii] + U[lft + ii])/pm.dZ/pm.dZ;
    return;
}

@nogc
void compute_residual(GasModel gm, ThermochemicalReactor reactor, GasState gs, number[] omegaMi, ref const Parameters pm, number[] U, number[] U2nd, number[] R, bool v=false){
/*
    The residual or right hand side is the time derivatives of the equations.
    See Lapointe et al. equations (1) and (2)
*/
    size_t N = pm.N;
    size_t neq = pm.neq;
    size_t nsp = pm.nsp;

    foreach(i; 0 .. N){
        double arg = erfc_inv(2.0*pm.Z[i]);
        double chi = 0.5*exp(-2.0*arg*arg);
        size_t idx = i*neq;
        bool verbose = v &&(i==15);

        gs.T = U[idx+nsp];
        gs.p = pm.p;
        foreach(j, Yj; U[idx .. idx+nsp]) gs.massf[j] = Yj;
        gm.update_thermo_from_pT(gs);
        //if (isNaN(gs.rho) || (gs.rho <= 0.0)) {
        //    throw new Exception(format("Invalid density. Gasstate is %s", gs));
        //}
        reactor.eval_source_terms(gm, gs, omegaMi);
        number cp = gm.Cp(gs);

        R[idx+nsp] = chi/2.0*U2nd[idx+nsp];
        number asdf = 0.0;
        for(int isp=0; isp<nsp; isp++){
            double Mi = gm.mol_masses[isp];
            number hi = gm.enthalpy(gs, isp);

            R[idx+isp] = chi/2.0*U2nd[idx+isp] + omegaMi[isp]/gs.rho;
            debug { if (verbose) writefln("   chi/2.0*U2nd[idx+isp] %e U2nd[idx+isp] %e omegaMi[isp]/gs.rho %e", chi/2.0*U2nd[idx+isp], U2nd[idx+isp], omegaMi[isp]/gs.rho); }
            R[idx+nsp] -= 1.0/gs.rho/cp*omegaMi[isp]*hi;
            asdf -= 1.0/gs.rho/cp*omegaMi[isp]*hi;
        }
        debug{
            if (verbose) writefln("   T= chi/2.0*U2nd[idx+nsp] %e chi %e U2nd %e", chi/2.0*U2nd[idx+nsp], chi, U2nd[idx+nsp]);
            if (verbose) writefln("Computing residual for cell %d Z %e T %e", i, pm.Z[i], gs.T);
            if (verbose) writefln(" Y: %s ", gs.massf);
            if (verbose) writefln("   asdf= %e", asdf);
        }
    }
}

//void compute_jacobian(GasModel gm, ThermochemicalReactor reactor, double[] U0, double[] U1, double dZ, double p, size_t N, size_t neq, size_t nsp, number[] Up, double[] Z, number[] U, number[] U2nd, number[] R, number[] Rp, Matrix!double J){
///*
//    Fill out a sparse matrix with the derivatives of the governing equations
//    computed using real-valued finite differences.
//*/
//    J._data[] = 0.0;
//    foreach(i; 0 .. neq*N) Up[i] = U[i];
//    double eps = 1e-9;
//
//    // we're computing dRi/dUj, with j being the column and i the row index
//    foreach(j; 0 .. neq){
//        // We can perturb every every third cell and compute the residuals in one go.
//        // This is different to how Eilmer does it, where we can compute the residuals
//        // on a subset of cells which are known to be affected by a given perturbation.
//        foreach(loop; 0 .. 3){
//            for (size_t cell=loop; cell<N; cell+=3){
//                size_t idx = cell*neq + j;
//                Up[idx] += eps;
//            }
//
//            second_derivs_from_cent_diffs(dZ, N, neq, nsp, U0, U1, Up, U2nd);
//            compute_residual(gm, reactor, p, N, neq, nsp, Z, Up, U2nd, Rp);
//
//            for (size_t cell=loop; cell<N; cell+=3){
//                size_t lft = (cell-1)*neq;
//                size_t ctr = cell*neq;
//                size_t rgt = (cell+1)*neq;
//                size_t col = ctr + j; 
//
//                if (cell>0){ // Do left cell
//                    foreach(i; 0 .. neq){
//                        double dRdU = (Rp[lft+i] - R[lft+i])/eps;
//                        size_t row = lft + i;
//                        J[row, col] = dRdU;
//                    }
//                }
//
//                // do the centre cell
//                foreach(i; 0 .. neq){
//                    double dRdU = (Rp[ctr+i] - R[ctr+i])/eps;
//                    size_t row = ctr + i;
//                    J[row, col] = dRdU;
//                }
//
//                if (cell<N-1) { // Do right cell
//                    foreach(i; 0 .. neq){
//                        double dRdU = (Rp[rgt+i] - R[rgt+i])/eps;
//                        size_t row = rgt + i;
//                        J[row, col] = dRdU;
//                    }
//                }
//                // finally, undo the perturbation
//                size_t idx = cell*neq + j;
//                Up[idx] -= eps;
//            }
//        }
//    }
//    return;
//}

void compute_jacobian2(GasModel gm, ThermochemicalReactor reactor, GasState gs, ref const Parameters pm, number[] omegaMi, number[] Up, number[] U, number[] U2nd, number[] R, number[] Rp, Matrix!double J){
/*
    Fill out a dense matrix with the derivatives of the governing equations
    computed using real-valued finite differences.
*/
    size_t n = pm.n;

    J._data[] = 0.0;
    foreach(i; 0 .. n) Up[i] = U[i];
    double eps = 1e-16;

    // we're computing dRi/dUj, with j being the column and i the row index
    foreach(j; 0 .. n){
        Up[j].im = eps;

        second_derivs_from_cent_diffs(pm, Up, U2nd);
        compute_residual(gm, reactor, gs, omegaMi, pm, Up, U2nd, Rp);

        foreach(i; 0 .. n){
           double dRdU = (Rp[i].im)/eps;
           J[i, j] = dRdU;
        }
        Up[j].im = 0.0;
        foreach(i; 0 .. n) Rp[i].im = 0.0;
    }
    return;
}

struct Parameters {
    size_t nsp;
    size_t neq;
    size_t N;
    size_t n;

    double p;
    double dZ;
    double T0;
    double T1;

    double[] Z;
    double[] Y0;
    double[] Y1;
    double[] U0;
    double[] U1;
}

void write_solution_to_file(ref const Parameters pm, number[] U, string filename){

    File outfile = File(filename, "wb");
    size_t[1] ibuff; double[1] dbuff; // buffer arrays

    ibuff[0] = pm.nsp;  outfile.rawWrite(ibuff);
    ibuff[0] = pm.neq;  outfile.rawWrite(ibuff);
    ibuff[0] = pm.N;    outfile.rawWrite(ibuff);
    ibuff[0] = pm.n;    outfile.rawWrite(ibuff);

    dbuff[0] = pm.p;    outfile.rawWrite(dbuff);
    dbuff[0] = pm.dZ;   outfile.rawWrite(dbuff);
    dbuff[0] = pm.T0;   outfile.rawWrite(dbuff);
    dbuff[0] = pm.T1;   outfile.rawWrite(dbuff);

    foreach(i; 0 .. pm.N)  { dbuff[0] = pm.Z[i];  outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.nsp){ dbuff[0] = pm.Y0[i]; outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.nsp){ dbuff[0] = pm.Y1[i]; outfile.rawWrite(dbuff); }
    foreach(i; 0 .. pm.n)  { dbuff[0] = U[i].re;     outfile.rawWrite(dbuff); }
    outfile.close();
    return;
}

void read_solution_from_file(ref const Parameters pm, number[] U, string filename){

    File infile = File(filename, "rb");
    size_t[1] ibuff; double[1] dbuff; // buffer arrays

    infile.rawRead(ibuff); if (pm.nsp != ibuff[0]) throw new Error(format("Sim nsp %d does not match file %d", pm.nsp, ibuff[0]));
    infile.rawRead(ibuff); if (pm.neq != ibuff[0]) throw new Error(format("Sim neq %d does not match file %d", pm.neq, ibuff[0]));
    infile.rawRead(ibuff); if (pm.N   != ibuff[0]) throw new Error(format("Sim N   %d does not match file %d", pm.N,   ibuff[0]));
    infile.rawRead(ibuff); if (pm.n   != ibuff[0]) throw new Error(format("Sim n   %d does not match file %d", pm.n,   ibuff[0]));

    // Let's not bother checking these. It's kind of whatever.
    infile.rawRead(dbuff); //if (pm.p   != dbuff[0]) throw new Error(format("Sim p   %e does not match file %e", pm.p,   dbuff[0]));
    infile.rawRead(dbuff); //if (pm.dZ  != dbuff[0]) throw new Error(format("Sim dZ  %e does not match file %e", pm.dZ,  dbuff[0]));
    infile.rawRead(dbuff); //if (pm.T0  != dbuff[0]) throw new Error(format("Sim T0  %e does not match file %e", pm.T0,  dbuff[0]));
    infile.rawRead(dbuff); //if (pm.T1  != dbuff[0]) throw new Error(format("Sim T1  %e does not match file %e", pm.T1,  dbuff[0]));

    // At this stage we actually don't want to override the Parameters struct
    foreach(i; 0 .. pm.N)  { infile.rawRead(dbuff); /* pm.Z[i]  = dbuff[0]; */}
    foreach(i; 0 .. pm.nsp){ infile.rawRead(dbuff); /* pm.Y0[i] = dbuff[0]; */}
    foreach(i; 0 .. pm.nsp){ infile.rawRead(dbuff); /* pm.Y1[i] = dbuff[0]; */}
    foreach(i; 0 .. pm.n)  { infile.rawRead(dbuff); U[i].re = dbuff[0]; }
    infile.close();
    return;
}

void gaussian_initial_condition(ref const Parameters pm, number[] U){
    double sigma = 0.1;

    foreach(i; 0 .. pm.N){
        size_t idx = i*pm.neq;
        double factor = 0.5*tanh(2*6.0*pm.Z[i] - 6.0) + 0.5;
        foreach(isp; 0 .. pm.nsp) U[idx+isp] = (1.0 - factor)*pm.Y0[isp] + factor*pm.Y1[isp];
        U[idx+pm.nsp].re = 1500.0*exp(-(pm.Z[i]-0.5)*(pm.Z[i]-0.5)/2.0/sigma/sigma) + 300.0;
        version(complex_numbers) { U[idx+pm.nsp].im = 0.0; }
    }
}

void RK4_explicit_time_increment(number[] U, number[] U2nd, number[] R,
                                 number[] Ua, number[] Ub, number[] Uc,
                                 number[] Ra, number[] Rb, number[] Rc,
                                 number[] omegaMi, GasState gs,
                                 number[] dU,
                                 double dt, GasModel gm, ThermochemicalReactor reactor, Parameters pm, bool verbose){

    size_t n = pm.n;
    second_derivs_from_cent_diffs(pm, U, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, U, U2nd, R, verbose);

    foreach(i; 0 .. n) Ua[i] = U[i] + dt/2.0*R[i];
    second_derivs_from_cent_diffs(pm, Ua, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, Ua, U2nd, Ra, false);

    foreach(i; 0 .. n) Ub[i] = U[i] + dt/2.0*Ra[i];
    second_derivs_from_cent_diffs(pm, Ub, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, Ub, U2nd, Rb, false);

    foreach(i; 0 .. n) Uc[i] = U[i] + dt*Rb[i];
    second_derivs_from_cent_diffs(pm, Uc, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, Uc, U2nd, Rc, false);

    foreach(i; 0 .. n) dU[i] = dt/6.0*(R[i] + 2.0*Ra[i] + 2.0*Rb[i] + Rc[i]);
    return;
}

void Euler_implicit_time_increment(number[] U, number[] U2nd, number[] R, number[] Up, number[] Rp, double[] Rr,
                                   Matrix!double LHS, Matrix!double J, 
                                   number[] omegaMi, GasState gs,
                                   number[] dU,
                                   double dt, GasModel gm, ThermochemicalReactor reactor, Parameters pm, bool verbose){

    second_derivs_from_cent_diffs(pm, U, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, U, U2nd, R, verbose);

    size_t n = pm.n;
    compute_jacobian2(gm, reactor, gs, pm, omegaMi, Up, U, U2nd, R,  Rp, J);
    foreach(i; 0 .. n){
        foreach(j; 0 .. n){
            LHS[i,j] = (i==j) ? 1.0/dt : 0.0;
            LHS[i,j] -= J[i,j];
        }
    }

    foreach(i, Ri; R) Rr[i] = Ri.re;
    auto perm = decomp!double(LHS);
    auto x = new Matrix!double(Rr);
    solve!double(LHS, x, perm);

    // x is solved in place to be delta U
    foreach(i; 0 .. n) dU[i] = x[i,0];
    return;
}

int main(string[] args)
{
    int exitFlag = 0; // Presume OK in the beginning.
    string name = "flame";
    if (args.length>1) name = args[1];
    
    GasModel gm = init_gas_model("gm.lua");
    ThermochemicalReactor reactor = init_thermochemical_reactor(gm, "rr.lua", "");
    GasState gs = GasState(gm);
    Parameters pm = Parameters();

    pm.nsp = gm.n_species;
    pm.neq = pm.nsp+1;

    pm.p = 75e3;
    pm.N = 32;
    pm.n = pm.N*pm.neq;
    pm.Z.length = pm.N;
    foreach(i; 1 .. pm.N+1) pm.Z[i-1] = i/(pm.N+1.0);
    pm.dZ = 1.0/(pm.N+1.0);

    // Boundary Conditions
    pm.T0 = 300.0; pm.T1 = 300.0;
    pm.Y0.length = pm.nsp; foreach(isp; 0 .. pm.nsp) pm.Y0[isp] = 0.0;
    pm.Y1.length = pm.nsp; foreach(isp; 0 .. pm.nsp) pm.Y1[isp] = 0.0;
    
    pm.Y0[gm.species_index("N2")] = 0.767;
    pm.Y0[gm.species_index("O2")] = 0.233;

    pm.Y1[gm.species_index("N2")] = 0.88;
    pm.Y1[gm.species_index("H2")] = 0.12;

    pm.U0.length = pm.neq; foreach(isp; 0 .. pm.nsp) pm.U0[isp] = pm.Y0[isp]; pm.U0[pm.nsp] = pm.T0;
    pm.U1.length = pm.neq; foreach(isp; 0 .. pm.nsp) pm.U1[isp] = pm.Y1[isp]; pm.U1[pm.nsp] = pm.T1;

    size_t neq = pm.neq;
    size_t N = pm.N;
    size_t n = pm.n;
    size_t nsp = pm.nsp;

    // Initial Guess
    number[] U,Up,dU;
    double[] Rr;
    number[] Ua,Ub,Uc,Ra,Rb,Rc;
    number[] R,Rp;
    number[] U2nd;

    U.length = neq*N;
    Up.length = neq*N;
    dU.length = neq*N;
    U2nd.length = neq*N;
    R.length = (neq)*N;
    Rp.length = (neq)*N;
    Ua.length = (neq)*N;
    Ub.length = (neq)*N;
    Uc.length = (neq)*N;
    Ra.length = (neq)*N;
    Rb.length = (neq)*N;
    Rc.length = (neq)*N;
    Rr.length = (neq)*N;
    number[] omegaMi; omegaMi.length = pm.nsp;

    // Set up Jacobian. We have a tridiagonal block matrix
    //auto J = new SMatrix!double();
    //J.aa.length = (neq*neq)*3*N - 2*(neq*neq);
    //foreach(i; 0 .. N){
    //    foreach(jj; 0 .. neq){
    //        // The ia array holds where in ja a certain row starts
    //        J.ia ~= J.ja.length; 

    //        // ia holds which column each entry would be in in a real matrix
    //        // we start with the left side block matrix, assuming this isn't
    //        // i==0, which has no neighbour on that side.
    //        if (i>0) {
    //            foreach(ii; 0 .. neq) J.ja ~= (i-1)*neq + ii;
    //        }
    //        // this cell block matrix entries
    //        foreach(ii; 0 .. neq) J.ja ~= i*neq + ii;
    //        // right cell block matrix entries
    //        if (i<N-1) {
    //            foreach(ii; 0 .. neq) J.ja ~= (i+1)*neq + ii;
    //        }
    //        // Then we do all that again for the next row down.
    //    }
    //}
    //writefln("Gotta add one more element to the array: J.ja.length %d J.aa.length %d", J.ja.length, J.aa.length);
    //J.ia ~= J.ja.length;

    size_t extent = neq*N;
    auto J = new Matrix!double(extent, extent);
    auto LHS = new Matrix!double(extent, extent);

    gaussian_initial_condition(pm, U);

    second_derivs_from_cent_diffs(pm, U, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, U, U2nd, R, false);

    double GR0 = 0.0; foreach(Ri; R) GR0 += Ri.re*Ri.re;
    GR0 = sqrt(GR0);


    immutable int maxiters = 10000;
    double dt = 1e-6;
    foreach(iter; 0 .. maxiters) {
        // Let's try computing some derivatives
        //compute_jacobian2(gm, reactor,  U0, U1, dZ, p, N, neq, nsp, Up, Z, U, U2nd, R,  Rp, J);

        // This might not be the best way of doing this but who knows.
        //decompILU0(J);
        //solve(J, R);
        //foreach(i, Ri; R) Rc[i] = Ri.re;
        //auto perm = decomp!double(J);
        //auto x = new Matrix!double(Rc);
        //solve!double(J, x, perm);

        //double global_relaxation_factor = 1.0;
        //double maxdY = -1.0;
        //foreach(cell; 0 .. N){
        //    foreach(j; 0 .. nsp){
        //        double dY = x[neq*cell+j,0].re;
        //        maxdY = fmax(maxdY, fabs(dY));
        //        double relaxation_factor = fmin(1.0, 1e-1/(fabs(dY)+1e-16));
        //        global_relaxation_factor = fmin(global_relaxation_factor, relaxation_factor);
        //    }
        //}
        //writefln("relax factor is %e from dYmax %e", global_relaxation_factor, maxdY);

        //// R is solved in place to be -delta U
        //foreach(i; 0 .. n) U[i] -= global_relaxation_factor*x[i,0];
        //foreach(cell; 0 .. N){
        //    foreach(j; 0 .. nsp){
        //        U[cell*neq + j] = fmin(fmax(U[cell*neq + j], 0.0), 1.0) ;
        //    }
        //}
        bool verbose = false;
        //if (iter%1000==0) verbose=true;

        if (iter==20){
            GR0 = 0.0; foreach(Ri; R) GR0 += Ri.re*Ri.re;
            GR0 = sqrt(GR0);
        }
        //RK4_explicit_time_increment(U,  U2nd,  R,  Ua,  Ub,  Uc, Ra,  Rb,  Rc, omegaMi, gs, dU, dt, gm, reactor, pm, verbose);
        Euler_implicit_time_increment(U,  U2nd,  R,  Up,  Rp, Rr, LHS, J, omegaMi, gs, dU, dt, gm, reactor, pm, verbose);

        foreach(i; 0 .. n) U[i] += dU[i];

        double GR = 0.0; foreach(Ri; R) GR += Ri.re*Ri.re;
        GR = sqrt(GR);
        double GRR = GR/GR0;

        if (iter%10==0){ 
            writefln("iter %d GR0 %e GR %e GRR %e", iter, GR0, GR, GRR);
        }
        if (GRR<1e-6) break;
        //if (iter==maxiters-1) throw new Error(format("Convergence failed after %s iterations, GRR was %e", maxiters, GRR));
        if (iter==maxiters-1) writefln("Warning: Simulation timed out at %d iterations with GRR %e", iter, GRR);
    }
    write_solution_to_file(pm, U, format("%s.sol", name));
    //read_solution_from_file(pm, Ua, format("%s.sol", name));
    writefln("Done!");


    return exitFlag;
} // end main()

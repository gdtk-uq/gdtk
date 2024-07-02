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

/* This function borrowed from zero_rk:
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

void multiply_tridigonal_block_matrix(ref const Parameters pm, Matrix!(double)[3][] LHS, Matrix!(double)[] U, Matrix!(double)[] R){
    size_t neq = pm.neq;
    size_t N   = pm.N;


    R[0] = LHS[0][1].dot(U[0]) + LHS[0][2].dot(U[0+1]);
    foreach(cell; 1 .. N-1){
        R[cell] = LHS[cell][0].dot(U[cell-1]) + LHS[cell][1].dot(U[cell]) + LHS[cell][2].dot(U[cell+1]);
    }
    R[N-1] = LHS[N-1][0].dot(U[N-2]) + LHS[N-1][1].dot(U[N-1]);
    return;

}

void solve_tridigonal_block_matrix(ref const Parameters pm, ref TridiagonalSolveWorkspace tdws, Matrix!(double)[3][] LHS, Matrix!(double)[] U, Matrix!(double)[] R){
/*
    Solve a sparse tridiagonal block matrix problem. We assume the data is packed as follows:
    [[##, B0, C0],
     [A1, B1, C1],
         ...
     [AN-1, BN-1, ##]]

*/

    size_t neq = pm.neq;
    size_t N   = pm.N;
    auto Ai = tdws.Ai;
    auto Bi = tdws.Bi;
    auto Ci = tdws.Ci;
    auto Ui = tdws.Ui;
    auto Ri = tdws.Ri;

    // First eliminate the A blocks using block row multiply and adds
    foreach(step; 0 .. N-1){
        Bi._data[] = LHS[step][1]._data[];
        Ci._data[] = LHS[step][2]._data[];
        Ri._data[] = R[step]._data[];

        auto perm = decomp!double(Bi);
        solve!double(Bi, Ri, perm);
        solve!double(Bi, Ci, perm);
        R[step]._data[] = Ri._data[];
        LHS[step][2]._data[] = Ci._data[];
        LHS[step][1].eye();

        Ai._data[] = LHS[step+1][0]._data[];
        LHS[step+1][0].zeros();
        LHS[step+1][1] -= Ai.dot(Ci);
        R[step+1]      -= Ai.dot(Ri);
    }

    // Now do a backward substituion on the remaining blocks to solve for U
    size_t end = N-1;
    Bi._data[] = LHS[end][1]._data[];
    Ri._data[] = R[end]._data[];

    auto perm = decomp!double(Bi);
    solve!double(Bi, Ri, perm);
    U[end]._data[] = Ri._data[];

    // Danger Looping down an unsigned integer is a bad idea...
    for (size_t step=N-2; step>=0; step--){
        Bi._data[] = LHS[step][1]._data[];
        Ci._data[] = LHS[step][2]._data[];
        Ri._data[] = R[step]._data[];
        Ui._data[] = U[step+1]._data[];

        Ri -= Ci.dot(Ui);

        perm = decomp!double(Bi);
        solve!double(Bi, Ri, perm);
        U[step]._data[] = Ri._data[];
        if (step==0) break;
    }

    return;
}



void compute_sparse_jacobian(GasModel gm, ThermochemicalReactor reactor, GasState gs, ref const Parameters pm, number[] omegaMi, number[] Up, number[] U, number[] U2nd, number[] R, number[] Rp, Matrix!(double)[3][] J){
/*
    Fill out a dense matrix with the derivatives of the governing equations
    computed using complex-valued finite differences.
*/
    size_t N = pm.N;
    size_t n = pm.n;
    size_t neq = pm.neq;

    foreach(i; 0 .. n) Up[i] = U[i];
    double eps = 1e-16;

    // we're computing dRi/dUj, with j being the column and i the row index
    foreach(j; 0 .. neq){
        // We can perturb every every third cell and compute the residuals in one go.
        // This is different to how Eilmer does it, where we can compute the residuals
        // on a subset of cells which are known to be affected by a given perturbation.
        foreach(loop; 0 .. 3){
            for (size_t cell=loop; cell<N; cell+=3){
                size_t idx = cell*neq + j;
                Up[idx].im = eps;
            }

            second_derivs_from_cent_diffs(pm, Up, U2nd);
            compute_residual(gm, reactor, gs, omegaMi, pm, Up, U2nd, Rp);

            for (size_t cell=loop; cell<N; cell+=3){
                size_t lft = (cell-1);
                size_t ctr = cell;
                size_t rgt = (cell+1);

                if (cell>0){ // Do left cell
                    foreach(i; 0 .. neq){
                        double dRdU = (Rp[lft*neq+i].im)/eps;
                        J[lft][2][i, j] = dRdU;
                    }
                }

                // do the centre cell
                foreach(i; 0 .. neq){
                    double dRdU = (Rp[ctr*neq+i].im)/eps;
                    J[ctr][1][i, j] = dRdU;
                }

                if (cell<N-1) { // Do right cell
                    foreach(i; 0 .. neq){
                        double dRdU = (Rp[rgt*neq+i].im)/eps;
                        J[rgt][0][i, j] = dRdU;
                    }
                }
                // finally, undo the perturbation
                size_t idx = cell*neq + j;
                Up[idx].im = 0.0;
            }
        }
    }
    return;
}

void compute_jacobian(GasModel gm, ThermochemicalReactor reactor, GasState gs, ref const Parameters pm, number[] omegaMi, number[] Up, number[] U, number[] U2nd, number[] R, number[] Rp, Matrix!double J){
/*
    Fill out a dense matrix with the derivatives of the governing equations
    computed using complex-valued finite differences.
*/
    size_t N = pm.N;
    size_t n = pm.n;
    size_t neq = pm.neq;

    J._data[] = 0.0;
    foreach(i; 0 .. n) Up[i] = U[i];
    double eps = 1e-16;

    // we're computing dRi/dUj, with j being the column and i the row index
    foreach(j; 0 .. neq){
        // We can perturb every every third cell and compute the residuals in one go.
        // This is different to how Eilmer does it, where we can compute the residuals
        // on a subset of cells which are known to be affected by a given perturbation.
        foreach(loop; 0 .. 3){
            for (size_t cell=loop; cell<N; cell+=3){
                size_t idx = cell*neq + j;
                Up[idx].im = eps;
            }

            second_derivs_from_cent_diffs(pm, Up, U2nd);
            compute_residual(gm, reactor, gs, omegaMi, pm, Up, U2nd, Rp);

            for (size_t cell=loop; cell<N; cell+=3){
                size_t lft = (cell-1)*neq;
                size_t ctr = cell*neq;
                size_t rgt = (cell+1)*neq;
                size_t col = ctr + j;

                if (cell>0){ // Do left cell
                    foreach(i; 0 .. neq){
                        double dRdU = (Rp[lft+i].im)/eps;
                        size_t row = lft + i;
                        J[row, col] = dRdU;
                    }
                }

                // do the centre cell
                foreach(i; 0 .. neq){
                    double dRdU = (Rp[ctr+i].im)/eps;
                    size_t row = ctr + i;
                    J[row, col] = dRdU;
                }

                if (cell<N-1) { // Do right cell
                    foreach(i; 0 .. neq){
                        double dRdU = (Rp[rgt+i].im)/eps;
                        size_t row = rgt + i;
                        J[row, col] = dRdU;
                    }
                }
                // finally, undo the perturbation
                size_t idx = cell*neq + j;
                Up[idx].im = 0.0;
            }
        }
    }
    return;
}

void compute_jacobian2(GasModel gm, ThermochemicalReactor reactor, GasState gs, ref const Parameters pm, number[] omegaMi, number[] Up, number[] U, number[] U2nd, number[] R, number[] Rp, Matrix!double J){
/*
    Fill out a dense matrix with the derivatives of the governing equations
    computed using complex finite differences.
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

struct TridiagonalSolveWorkspace{
    Matrix!double Ai;
    Matrix!double Bi;
    Matrix!double Ci;
    Matrix!double Ui;
    Matrix!double Ri;

    this(size_t neq){
        Ai = new Matrix!double(neq, neq);
        Bi = new Matrix!double(neq, neq);
        Ci = new Matrix!double(neq, neq);
        Ui  = new Matrix!double(neq);
        Ri  = new Matrix!double(neq);
    }
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

void Euler_implicit_time_increment(number[] U, number[] U2nd, number[] R, number[] Up, number[] Rp,
                                   Matrix!(double)[3][] J2, Matrix!(double)[] U2, Matrix!(double)[] R2,
                                   ref TridiagonalSolveWorkspace tdws,
                                   number[] omegaMi, GasState gs,
                                   number[] dU,
                                   double dt, GasModel gm, ThermochemicalReactor reactor, Parameters pm, bool verbose){


    size_t neq = pm.neq;
    size_t N = pm.N;
    size_t n = pm.n;


    second_derivs_from_cent_diffs(pm, U, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, U, U2nd, R, false);

    compute_sparse_jacobian(gm, reactor, gs, pm, omegaMi, Up, U, U2nd, R,  Rp, J2);

    // Get put into [I/dt - J] U = R form, for the Euler/Newton update
    foreach(cell; 0 .. N){
        J2[cell][0]._data[] *= -1.0; 
        J2[cell][1]._data[] *= -1.0; 
        J2[cell][2]._data[] *= -1.0; 
        foreach(i; 0 .. neq) J2[cell][1][i, i] += 1.0/dt;
        foreach(i; 0 .. neq) R2[cell][i, 0] = R[cell*neq + i].re;
    }

    solve_tridigonal_block_matrix(pm, tdws, J2, U2, R2);

    foreach(cell; 0 .. N){
        foreach(j; 0 .. neq){
            dU[cell*neq +  j] = U2[cell][j, 0];
        }
    }

    return;
}

double compute_line_search(GasModel gm, ThermochemicalReactor reactor, GasState gs, ref const Parameters pm, number[] omegaMi, number[] dU, number[] Up, number[] U, number[] U2nd, number[] R, number[] Rp, double GRold){
/*
    Limit the size of the update to ensure that we get a reduction in the residual.

*/
    
    size_t n = pm.n;
    double lambda = 1.0;

    while (lambda>1e-6){
        foreach(i; 0 .. n) Up[i] = U[i] + lambda*dU[i];
        second_derivs_from_cent_diffs(pm, Up, U2nd);
        compute_residual(gm, reactor, gs, omegaMi, pm, Up, U2nd, Rp, false);

        double GR = 0.0; foreach(Ri; Rp) GR += Ri.re*Ri.re;
        GR = sqrt(GR);

        if (GR>GRold) {
            lambda *= 0.5;
        } else {
            return lambda;
        }
    }
    throw new Error("Line search failed to find an acceptable relaxation factor");
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
    pm.N = 48;
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
    number[] omegaMi; omegaMi.length = pm.nsp;

    // Set up Jacobian. We have a tridiagonal block matrix
    Matrix!(double)[3][] J2;
    Matrix!(double)[] U2;
    Matrix!(double)[] R2;

    foreach(cell; 0 .. N) {
        auto Jc0 = new Matrix!double(neq, neq);
        auto Jc1 = new Matrix!double(neq, neq);
        auto Jc2 = new Matrix!double(neq, neq);
        J2 ~= [Jc0, Jc1, Jc2];
        auto Unew = new Matrix!double(neq);
        U2 ~= Unew;
        auto Rnew = new Matrix!double(neq);
        R2 ~= Rnew;
    }

    TridiagonalSolveWorkspace tdws = TridiagonalSolveWorkspace(neq);

    gaussian_initial_condition(pm, U);
    second_derivs_from_cent_diffs(pm, U, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, U, U2nd, R, false);

    double GRmax = 0.0; foreach(Ri; R) GRmax += Ri.re*Ri.re;
    GRmax = sqrt(GRmax);


    immutable int maxiters = 4000;
    double dt = 5e-7;
    bool verbose = false;
    double GRRold = 1.0;
    double tau = 4e-2;
    foreach(iter; 0 .. maxiters) {

        //RK4_explicit_time_increment(U,  U2nd,  R,  Ua,  Ub,  Uc, Ra,  Rb,  Rc, omegaMi, gs, dU, dt, gm, reactor, pm, verbose);
        Euler_implicit_time_increment(U,  U2nd,  R,  Up,  Rp, J2, U2, R2, tdws, omegaMi, gs, dU, dt, gm, reactor, pm, verbose);

        foreach(i; 0 .. n) U[i] += dU[i];

        double GR = 0.0; foreach(Ri; R) GR += Ri.re*Ri.re;
        GR = sqrt(GR);
        if (iter<200) GRmax = fmax(GR, GRmax);
        double GRR = GR/GRmax;

        if (iter%10==0){
            writefln("iter %d dt %e GRmax %e GR %e GRR %e", iter, dt, GRmax, GR, GRR);
        }
        if (GRR<tau){
            double dtnew = dt*pow(GRRold/GRR, 1.6);
            dtnew = fmin(dtnew, 1.2*dt);
            dtnew = fmax(dtnew, 0.5*dt);
            dtnew = fmin(dtnew, 1.0);
            dt = dtnew;
        }
        GRRold = GRR;

        if (GRR<1e-6) break;
        if (iter==maxiters-1) writefln("Warning: Simulation timed out at %d iterations with GRR %e", iter, GRR);
    }
    write_solution_to_file(pm, U, format("%s.sol", name));
    //read_solution_from_file(pm, Ua, format("%s.sol", name));
    writefln("Done!");


    return exitFlag;
} // end main()

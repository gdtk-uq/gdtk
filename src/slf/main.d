// main.d: a steady laminar flamelet calculator by NNG.

module slf;

import std.stdio;
import std.math;
import std.mathspecial;
import std.format;
import std.string;
import std.conv;
import std.datetime.stopwatch : StopWatch;

// slf specific modules
import io;
import linalg;
import misc;

import gas;
import gas.physical_constants;
import kinetics;
import nm.bbla;
import nm.number;
import nm.complex;


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
        number chi = pm.D*exp(-2.0*arg*arg);
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
                                   double dt, GasModel gm, ThermochemicalReactor reactor, ref const Parameters pm, bool verbose){


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

pure @nogc
double compute_global_residual(in number[] R) {
    double GR = 0.0;
    foreach(Ri; R) {
        GR += Ri.re*Ri.re;
    }
    GR = sqrt(GR);
    return GR;
}

pure @nogc
double update_dt(double GRRold, double GRR, double dt){
    double dtnew = dt*pow(GRRold/GRR, 1.6);
    dtnew = fmin(dtnew, 1.2*dt);
    dtnew = fmax(dtnew, 0.5*dt);
    dtnew = fmin(dtnew, 1.0);
    return dtnew;
}

double[] get_derivatives_from_adjoint(number[] U, number[] U2nd, number[] R, number[] Up, number[] Rp,
                                   Matrix!(double)[3][] J2, Matrix!(double)[] U2, Matrix!(double)[] R2,
                                   ref TridiagonalSolveWorkspace tdws,
                                   number[] omegaMi, GasState gs,
                                   number[] dU,
                                   double dt, GasModel gm, ThermochemicalReactor reactor, ref Parameters pm, bool verbose){

    size_t neq = pm.neq;
    size_t N = pm.N;
    size_t n = pm.n;

    second_derivs_from_cent_diffs(pm, U, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, U, U2nd, R, false);

    compute_sparse_jacobian(gm, reactor, gs, pm, omegaMi, Up, U, U2nd, R,  Rp, J2);

    // compute the RHS vector, which is the partial derivate of R with x constant
    // and p changed
    double eta = 1e-12;
    pm.p.im = eta;
    second_derivs_from_cent_diffs(pm, U, U2nd);
    compute_residual(gm, reactor, gs, omegaMi, pm, U, U2nd, Rp, false);


    foreach(cell; 0 .. N){
        J2[cell][0]._data[] *= -1.0; 
        J2[cell][1]._data[] *= -1.0; 
        J2[cell][2]._data[] *= -1.0; 
        foreach(i; 0 .. neq) R2[cell][i, 0] = (Rp[cell*neq + i].im)/(eta);
    }
    solve_tridigonal_block_matrix(pm, tdws, J2, U2, R2);

    double[] dUdp; dUdp.length = pm.n;
    foreach(cell; 0 .. N){
        foreach(j; 0 .. neq){
            dUdp[cell*neq +  j] = U2[cell][j, 0];
        }
    }
    pm.p.im = 0.0;
    return dUdp;
}

int main(string[] args)
{
    int exitFlag = 0;
    string name = "slf";
    if (args.length>1) name = args[1];

    Config config = read_config_from_file(format("%s.yaml", name));

    GasModel gm = init_gas_model(config.gas_file_name);
    ThermochemicalReactor reactor = init_thermochemical_reactor(gm, config.reaction_file_name, "");
    GasState gs = GasState(gm);

    Parameters pm = Parameters(config, gm);
    StopWatch sw;


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

    immutable int maxiters = 8000;
    immutable double targetGRR = 1e-12;
    double dt = 5e-7;
    bool verbose = false;
    double GRRold = 1.0;
    double GRold = 1e99;
    double tau = 1e-1;
    string[] log;

    sw.start();
    foreach(iter; 0 .. maxiters) {

        //RK4_explicit_time_increment(U,  U2nd,  R,  Ua,  Ub,  Uc, Ra,  Rb,  Rc, omegaMi, gs, dU, dt, gm, reactor, pm, verbose);
        Euler_implicit_time_increment(U,  U2nd,  R,  Up,  Rp, J2, U2, R2, tdws, omegaMi, gs, dU, dt, gm, reactor, pm, verbose);

        foreach(i; 0 .. n) U[i] += dU[i];

        double GR = compute_global_residual(R);
        if (iter<150) GRmax = fmax(GR, GRmax);
        double GRR = GR/GRmax;

        string output = format("iter %d dt %e GRmax %e GR %e GRR %e", iter, dt, GRmax, GR, GRR);
        log ~= output;
        if (iter%50==0){
            writeln(output);
        }
        if ((GRR<tau) && (iter>150)){
            dt = update_dt(GRRold, GRR, dt);
        }
        GRRold = GRR;
        GRold = GR;

        if (GRR<targetGRR) {
            writefln("Reached target residual %e (%e), stopping...", targetGRR, GRR);
            break;
        }
        if (iter==maxiters-1) writefln("Warning: Simulation timed out at %d iterations with GRR %e", iter, GRR);
    }
    sw.stop();
    long wall_clock_elapsed = sw.peek.total!"seconds";

    write_solution_to_file(pm, U, format("%s.sol", name));
    write_log_to_file(log, format("%s.log", name));

    //double[] dUdp;
    //dUdp = get_derivatives_from_adjoint(U,  U2nd,  R,  Up,  Rp, J2, U2, R2, tdws, omegaMi, gs, dU, dt, gm, reactor, pm, verbose);
    //write_solution_to_file(pm, dUdp, format("%s-derivatives.sol", name));
    writefln("Done simulation in %d seconds.", wall_clock_elapsed);


    return exitFlag;
} // end main()

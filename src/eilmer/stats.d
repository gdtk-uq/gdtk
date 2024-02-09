/**
 * stats.d
 * Code for computing flow statistics.
 *
 * Author: Nick Gibbons
 * Version: 2024-02-02
 */

module stats;

import std.string;
import std.format;
import std.array;
import std.stdio;
import std.math;
import std.file;

import geom;
import flowstate;
import gas;
import turbulence;


struct StatInstance {
    double rho, velx, vely, velz, p, a, mu, k, mu_t, k_t, S, u, T;
    double[] turb;
    double[] massf;

    this(size_t nturb, size_t nspecies){
        turb.length = nturb;
        massf.length = nspecies;
    }

    void reset() {
        rho = 0.0;
        velx = 0.0;
        vely = 0.0;
        velz = 0.0;
        p = 0.0;
        a = 0.0;
        mu = 0.0;
        k = 0.0;
        mu_t = 0.0;
        k_t = 0.0;
        S = 0.0;
        u = 0.0;
        T = 0.0;
        turb[] = 0.0;
        massf[] = 0.0;
    }

    void plus_equals(double weight, FlowState* fs){
        rho   += weight * fs.gas.rho.re; 
        velx  += weight * fs.vel.x.re; 
        vely  += weight * fs.vel.y.re;
        velz  += weight * fs.vel.z.re;
        p     += weight * fs.gas.p.re;
        a     += weight * fs.gas.a.re;
        mu    += weight * fs.gas.mu.re;
        k     += weight * fs.gas.k.re; 
        mu_t  += weight * fs.mu_t;
        k_t   += weight * fs.k_t.re;
        S     += weight * fs.S.re;
        u     += weight * fs.gas.u;
        T     += weight * fs.gas.T;
        version(turbulence) { foreach(it; 0 .. turb.length) turb[it] += weight*fs.turb[it].re; }
        version(multi_species_gas) { foreach(isp; 0 .. massf.length) massf[isp] += weight*fs.gas.massf[isp].re; }
    }
}


class FlowStats {
    // Compute and store statistical moments of the flow
    immutable double zmax =  10e-3;
    immutable double zmin = -10e-3;
    immutable size_t NTIME = 100;
    size_t NSPACE;
    double[] z;
    double dz, filter_width, b_squared;
    size_t tidx=0;

    double[] weights;
    double[NTIME] times;
    StatInstance[][NTIME] data;

    this(size_t nspecies, size_t nturb, size_t NSPACE) {
        this.NSPACE=NSPACE;
        dz = (zmax - zmin)/NSPACE;
        filter_width = 2.0*dz;
        b_squared = filter_width*filter_width;

        z.length = NSPACE;
        foreach(j; 0 .. NSPACE) {
            z[j] = zmin + dz/2.0 + j*dz;
        }

        weights.length = NSPACE;

        foreach(t; 0 .. NTIME) {
            foreach(i; 0 .. NSPACE) {
                data[t] ~= StatInstance(nspecies, nturb);
            }
        }

        weights[] = 0.0;
        reset_buffers();
    }

    void init_map_arrays_for_block(Vector3[] positions, ref size_t[][] cell_to_bin_map, ref double[][] bin_weights){
    /*
        We don't want to have to compute the weights and stuff every time, so each
        block should store a precomputed array of them.

        Maybe this should be called with FluidBLock? Or cell data? Good questions.
    */

        size_t[30] nbins;
        nbins[] = 0;
        bin_weights.length = positions.length;
        cell_to_bin_map.length = positions.length;

        immutable double threshold = 1e-4;

        foreach(i; 0 .. positions.length){
            foreach(j, zj; z){
                double dist = (positions[i].y - zj);
                double arg = -1.0*(dist*dist)/2.0/b_squared;
                double weight = exp(arg);

                // Many of the cells will be quite far from the bin,
                // so discard any that have weights below a threshold
                if (weight < threshold) continue;

                // otherwise let's save both the weight and record of which bins
                // this cell cares about.
                cell_to_bin_map[i] ~= j;
                bin_weights[i]     ~= weight;

                // Also, the actual weight array in bin-space does change with time,
                // so we can compute that here and keep it saved.
                weights[j] += weight;
            }
            nbins[cell_to_bin_map[i].length] += 1;
        }
        //foreach(i, j; nbins){
        //    writefln("%d cells have %d bins", j, i);
        //}
    }

    void update(double time, size_t[][] cell_to_bin_map, double[][] bin_weights, FlowState[] fss){
        assert(tidx<=NTIME);
        size_t[] ncheck;
        ncheck.length = NSPACE;
        ncheck[] = 0;

        times[tidx] = time;
        foreach(i; 0 .. cell_to_bin_map.length){
            foreach(jj, bin; cell_to_bin_map[i]){
                double weight = bin_weights[i][jj];
                data[tidx][bin].plus_equals(weight, &(fss[i]));
                ncheck[bin] += 1;
            }
        }

        //foreach(bin, ncells; ncheck) {
        //    writefln("Bin %d received contributions from %d cells", bin, ncells);
        //}
    }

    void increment_time_index(){
        tidx += 1;
    } 

    void reset_time_index(){
        tidx = 0;
    } 

    void dump_stats_to_file(TurbulenceModel tm, GasModel gm){
        // First check if the file exists and write the header, and then the positions
        // of the averaging kernels

        immutable string filename = "stats.txt";
        if (!exists(filename)) {
            auto statsfile = File(filename, "w");
            statsfile.write("rho velx vely velz p a mu k mu_t k_t S u T");
            version(turbulence) {
                foreach(it; 0 .. tm.nturb){
                    statsfile.writef(" %s", tm.primitive_variable_name(it));
                }
            }
            version(multi_species_gas) {
                foreach (isp; 0 .. gm.n_species) {
                    statsfile.writef(" Y_%s", gm.species_name(isp));
                }
            }
            statsfile.write("\n");
            statsfile.close();
        }

        // IN MPI we do a reduce here and only the master program writes.
        // We mnight have to rethink this if using this code for full 
        // grid stats. I guess in that case we can do block local writes
        // to separate files.

        auto f = File(filename, "a");
        foreach(t; 0 .. tidx){
            f.writef("t= %.18e n= %d\n", times[t], NSPACE);
            foreach(j; 0 .. NSPACE){
                StatInstance* s = &(data[t][j]);
                double weight = weights[j];
                f.writef("%.18e",  s.rho.re/weight);
                f.writef(" %.18e", s.velx.re/weight);
                f.writef(" %.18e", s.vely.re/weight);
                f.writef(" %.18e", s.velz.re/weight);
                f.writef(" %.18e", s.p.re/weight);
                f.writef(" %.18e", s.a.re/weight);
                f.writef(" %.18e", s.mu.re/weight);
                f.writef(" %.18e", s.k.re/weight);
                f.writef(" %.18e", s.mu_t.re/weight);
                f.writef(" %.18e", s.k_t.re/weight);
                f.writef(" %.18e", s.S.re/weight);
                f.writef(" %.18e", s.u.re/weight);
                f.writef(" %.18e", s.T.re/weight);
                version(turbulence) {
                    foreach(turb; s.turb){
                        f.writef(" %.18e", turb.re/weight);
                    }
                }
                version(multi_species_gas) {
                    foreach (Y; s.massf) { f.writef(" %.18e", Y.re/weight); }
                }
                f.write("\n");
            }
        }
        f.close();
    }

    void reset_buffers(){
        foreach(t; 0 .. NTIME) {
            times[t] = 0.0;
            foreach(j; 0 .. NSPACE) {
                data[t][j].reset();
            }
        }
    }
}
        


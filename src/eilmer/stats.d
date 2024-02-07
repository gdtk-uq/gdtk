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

    double[][NTIME] weights;
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
        foreach(t; 0 .. NTIME) {
            weights[t].length = NSPACE;
            foreach(i; 0 .. NSPACE) {
                data[t] ~= StatInstance(nspecies, nturb);
            }
        }
        reset_buffers();
    }

    void update(double time, Vector3[] positions, FlowState[] fss){
        assert(tidx<=NTIME);

        foreach(i; 0 .. positions.length){
            // TDB: Actually we want to store these weights and indices for any cells 
            // whose weights are above a threshold, which means ragged arrays stored 
            // in the relevant blocks.
            times[tidx] = time;
            foreach(j, zj; z){
                double dist = (positions[i].y - zj);
                double arg = -1.0*(dist*dist)/2.0/b_squared;
                double weight = exp(arg);
                weights[tidx][j] += weight;
                data[tidx][j].plus_equals(weight, &(fss[i]));
            }
        }
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
            statsfile.write("rho velx vely velz p a mu k mu_t k_t S");
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
                double weight = weights[t][j];
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
                weights[t][j] = 0.0;
                data[t][j].reset();
            }
        }
    }
}
        


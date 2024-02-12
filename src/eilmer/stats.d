/**
 * stats.d
 * Code for computing flow statistics.
 *
 * Author: Nick Gibbons
 * Version: 2024-02-02
 */

module stats;

import std.format;
import std.array;
import std.stdio;
import std.conv;
import std.math;
import std.file;

import geom;
import flowstate;
import gas;
import turbulence;
version(mpi_parallel) {
    import mpi;
}


struct StatsData {
    size_t ntimes, nspace, nvars, nturb, nspecies;
    int nbuffer;
    double[] buffer;

    this(size_t ntimes, size_t nspace, size_t nturb, size_t nspecies){
        this.ntimes = ntimes;
        this.nspace = nspace;
        this.nvars = 13+nturb+nspecies;
        this.nturb = nturb;
        this.nspecies=nspecies;
        buffer.length = ntimes*nspace*nvars;
        nbuffer = to!int(buffer.length);
    }

    void reset() {
        buffer[] = 0.0;
    }

    void plus_equals(size_t srt, size_t end, double weight, FlowState* fs){
        double[] section = buffer[srt .. end];

        section[0]  += weight * fs.gas.rho.re;
        section[1]  += weight * fs.vel.x.re;
        section[2]  += weight * fs.vel.y.re;
        section[3]  += weight * fs.vel.z.re;
        section[4]  += weight * fs.gas.p.re;
        section[5]  += weight * fs.gas.a.re;
        section[6]  += weight * fs.gas.mu.re;
        section[7]  += weight * fs.gas.k.re;
        section[8]  += weight * fs.mu_t;
        section[9]  += weight * fs.k_t.re;
        section[10] += weight * fs.S.re;
        section[11] += weight * fs.gas.u;
        section[12] += weight * fs.gas.T;
        version(turbulence) { foreach(it; 0 .. nturb) section[13+it] += weight*fs.turb[it].re; }
        version(multi_species_gas) {
            foreach(isp; 0 .. nspecies) section[13+nturb+isp] += weight*fs.gas.massf[isp].re;
        }
    }
}


class FlowStats {
    // Compute and store statistical moments of the flow
    immutable double zmax =  10e-3;
    immutable double zmin = -10e-3;
    immutable size_t NTIME = 100;
    size_t NSPACE, NVARS;
    double[] z;
    double dz, filter_width, b_squared;
    size_t tidx=0;

    double[] weights;
    double[NTIME] times;
    StatsData data;

    this(size_t nturb, size_t nspecies, size_t NSPACE) {
        this.NVARS = 13+nspecies+nturb;
        this.NSPACE=NSPACE;
        dz = (zmax - zmin)/NSPACE;
        filter_width = 2.0*dz;
        b_squared = filter_width*filter_width;

        z.length = NSPACE;
        foreach(j; 0 .. NSPACE) {
            z[j] = zmin + dz/2.0 + j*dz;
        }

        weights.length = NSPACE;
        data = StatsData(NTIME, NSPACE, nturb, nspecies);

        weights[] = 0.0;
        data.reset();
    }

    void init_map_arrays_for_block(Vector3[] positions, ref size_t[][] cell_to_bin_map, ref double[][] bin_weights){
    /*
        We don't want to have to compute the weights and stuff every time, so each
        block should store a precomputed array of them.

        Maybe this should be called with FluidBLock? Or cell data? Good questions.
    */

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
        }
    }

    void sum_weights_over_mpi(bool is_master_task){
        // In MPI only the master task needs to know the total weights
        version(mpi_parallel){
            int n = to!int(weights.length);

            if (is_master_task) {
                MPI_Reduce(MPI_IN_PLACE, weights.ptr, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            } else {
                MPI_Reduce(weights.ptr, weights.ptr, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
        }
    }

    void update(double time, size_t[][] cell_to_bin_map, double[][] bin_weights, FlowState[] fss){
        assert(tidx<=NTIME);

        times[tidx] = time;
        foreach(i; 0 .. cell_to_bin_map.length){
            foreach(jj, bin; cell_to_bin_map[i]){
                size_t srt = tidx*NSPACE*NVARS + bin*NVARS;
                size_t end = tidx*NSPACE*NVARS + bin*NVARS + NVARS;
                double weight = bin_weights[i][jj];
                data.plus_equals(srt, end, weight, &(fss[i]));
            }
        }
    }

    void increment_time_index(){
        tidx += 1;
    } 

    void reset_time_index(){
        tidx = 0;
    } 

    void dump_stats_to_file(bool is_master_task, TurbulenceModel tm, GasModel gm){
        // First check if the file exists and write the header, and then the positions
        // of the averaging kernels

        immutable string filename = "stats.txt";
        if (!exists(filename) && is_master_task) {
            writef(" Creating file stats.txt...");
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
        version(mpi_parallel){
            if (is_master_task) {
                MPI_Reduce(MPI_IN_PLACE, data.buffer.ptr, data.nbuffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            } else {
                MPI_Reduce(data.buffer.ptr, data.buffer.ptr, data.nbuffer, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
        }

        if (is_master_task) {
            writef(" flushing %d lines of data...", tidx);
            auto f = File(filename, "a");
            foreach(t; 0 .. tidx){
                f.writef("t= %.18e n= %d\n", times[t], NSPACE);
                foreach(j; 0 .. NSPACE){
                    double weight = weights[j];
                    foreach(k; 0 .. NVARS){
                        size_t index = t*(NSPACE*NVARS) + j*NVARS + k;
                        if (k==0) {
                            f.writef("%.18e",  data.buffer[index]/weight);
                        } else {
                            f.writef(" %.18e",  data.buffer[index]/weight);
                        }
                    }
                    f.write("\n");
                }
            }
            f.close();
            writefln(" Complete.");
        }
        return;
    }

    void reset_buffers(){
        data.reset();
        times[] = 0.0;
    }
}

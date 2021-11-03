/**
 * field.d
 * Machinery for solving electromagnetic fields through the fluid domain
 *
 * Author: Nick Gibbons
 * Version: 2021-01-20: Prototyping
 */

module field;

import std.conv;
import std.algorithm;
import std.format;
import std.stdio;
import std.math;
import nm.complex;
import nm.number;

import fvcell;
import fvinterface;
import fluidblock;
import geom;
import globalconfig;
import fieldbc;
import fieldgmres;
import fieldconductivity;
import fieldexchange;

class ElectricField {
    this(const FluidBlock[] localFluidBlocks) {
        N = 0;
        foreach(block; localFluidBlocks){
            block_offsets ~= N;
            N += to!int(block.cells.length);
        }
        A.length = N*nbands;
        Ai.length = N*nbands;
        b.length = N;
        max_iter = 100; // eventually you need to do something about this.
        phi.length = N;
        phi0.length = N;

        string conductivity_model_name = "test";
        conductivity = create_conductivity_model(conductivity_model_name); // TODO: Config variable?

        // I don't want random bits of the field module hanging off the boundary conditions.
        // Doing it this way is bad encapsulation, but it makes sure that other people only break my code
        // rather than the other way around.
        field_bcs.length = localFluidBlocks.length;
        foreach(i, block; localFluidBlocks){
            field_bcs[i].length = block.bc.length;
            foreach(j, bc; block.bc){
                field_bcs[i][j] = create_field_bc(bc.field_bc_name, bc, block_offsets, conductivity_model_name, N);
            }
        }

        version(mpi_parallel){
            exchanger = new Exchanger(GlobalConfig.mpi_rank_for_local_task, GlobalConfig.mpi_size, MPISharedField.instances);
            gmres = new GMResFieldSolver(exchanger);
        } else {
            gmres = new GMResFieldSolver();
        }

        return;
    }

    void solve_efield(FluidBlock[] localFluidBlocks) {
        A[] = 0.0;
        b[] = 0.0;
        Ai[] = -1;
    
        FVCell other;
        foreach(blkid, block; localFluidBlocks){
            foreach(cell; block.cells){
                int k = cell.id + block_offsets[blkid];
                Ai[nbands*k + 2] = k;

                foreach(io, face; cell.iface){
                    int iio = (io>1) ? to!int(io+1) : to!int(io); // -> [0,1,3,4] since 2 is the entry for "cell" 
                    double sign = cell.outsign[io];
        
                    if (face.is_on_boundary) {
                        int Aio;
                        double Akk, Ako, bk;
                        auto field_bc = field_bcs[blkid][face.bc_id];
                        field_bc(sign, face, cell, Akk, bk, Ako, Aio);
                        Ai[k*nbands + iio] = Aio;
                        A[k*nbands + iio] += Ako;
                        A[k*nbands + 2] += Akk;
                        b[k] += bk;
                    } else {
                        other = face.left_cell;
                        if (other==cell) other = face.right_cell;

                        double S = face.length.re;
                        double sigmaF = conductivity(face.fs, face.pos);
                        double d = distance_between(other.pos[0], cell.pos[0]);
                        Ai[k*nbands + iio] = other.id + block_offsets[blkid];
                        A[k*nbands + iio]+=  1.0*S/d*sigmaF;
                        A[k*nbands + 2] += -1.0*S/d*sigmaF;
                    }
                }
            }
        }
    
        // The Jacobi Preconditioner, a very simple scheme. Good for diagonally dominant matrices
        if (precondition){
            foreach(blkid, block; localFluidBlocks){
                foreach(cell; block.cells){
                    int k = cell.id + block_offsets[blkid];
                    double Akk = A[k*nbands + 2];
                    b[k] /= Akk;
                    foreach(iio; 0 .. nbands) A[k*nbands + iio] /= Akk;
                }
            }
        }
    
        phi0[] = 0.0;
        gmres.solve(N, nbands, A, Ai, b, phi0, phi, max_iter);

        // Unpack the solution into the "electric_potential" members stored in the cells
        size_t i = 0;
        foreach(block; localFluidBlocks){
            foreach(cell; block.cells){
                cell.electric_potential = phi[i];
                i += 1;
            }
        }
    }
private:
    immutable int nbands = 5; // 5 for a 2D structured grid
    immutable bool precondition = true;

    GMResFieldSolver gmres;
    ConductivityModel conductivity;
    FieldBC[][] field_bcs;
    int[] block_offsets; // FIXME: Badness with the block id's not matching their order in local fluid blocks
    version(mpi_parallel){
        Exchanger exchanger;
    }

    double[] A;
    int[] Ai;
    double[] b;
    double[] phi,phi0;
    int max_iter;
    int N;
}

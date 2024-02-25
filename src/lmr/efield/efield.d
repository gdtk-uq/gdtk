/**
 * Machinery for solving electromagnetic fields through the fluid domain
 *
 * Author: Nick Gibbons
 * Version: 2021-01-20: Prototyping
 */

module efield;

import std.conv;
import std.algorithm;
import std.format;
import std.stdio;
import std.math;
import ntypes.complex;
import nm.number;

import lmr.fluidfvcell;
import fvinterface;
import fluidblock;
import geom;
import globalconfig;
import efieldbc;
import efieldgmres;
import efieldconductivity;
import efieldexchange;
import efieldderivatives;
version(mpi_parallel){
    import mpi;
}

immutable uint ZNG_interior = 0b0000;
immutable uint ZNG_north    = 0b0001;
immutable uint ZNG_east     = 0b0010;
immutable uint ZNG_south    = 0b0100;
immutable uint ZNG_west     = 0b1000;
immutable uint[4] ZNG_types = [ZNG_north, ZNG_east, ZNG_south, ZNG_west];

class ElectricField {
    this(const FluidBlock[] localFluidBlocks, const string field_conductivity_model) {
        N = 0;
        foreach(block; localFluidBlocks){
            block_offsets ~= N;
            N += to!int(block.cells.length);
        }
        A.length = N*nbands;
        Ai.length = N*nbands;
        b.length = N;
        max_iter = N; // eventually you need to do something about this.
        phi.length = N;
        phi0.length = N;
        auto gmodel = GlobalConfig.gmodel_master;
        conductivity = create_conductivity_model(field_conductivity_model, gmodel);

        // I don't want random bits of the field module hanging off the boundary conditions.
        // Doing it this way is bad encapsulation, but it makes sure that other people only break my code
        // rather than the other way around.
        field_bcs.length = localFluidBlocks.length;
        foreach(i, block; localFluidBlocks){
            field_bcs[i].length = block.bc.length;
            foreach(j, bc; block.bc){
                field_bcs[i][j] = create_field_bc(bc.field_bc, bc, block_offsets, field_conductivity_model, N);
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

    void solve_efield(FluidBlock[] localFluidBlocks, bool verbose) {
        A[] = 0.0;
        b[] = 0.0;
        Ai[] = -1;

        FluidFVCell other;
        foreach(blkid, block; localFluidBlocks){
            auto gmodel = block.myConfig.gmodel;
            foreach(cell; block.cells){

                // Set up the sparse matrix indexes. This actually only needs to be done once, technically
                int k = cell.id + block_offsets[blkid];
                Ai[nbands*k + 2] = k;

                // We need the distances to the unknown phi locations, which may be faces or cells
                double[4] dx, dy, nx, ny;
                foreach(io, face; cell.iface){
                    double sign = cell.outsign[io];
                    Vector3 pos;
                    int iio = (io>1) ? to!int(io+1) : to!int(io); // -> [0,1,3,4] since 2 is the entry for "cell"

                    if (face.is_on_boundary) {
                        auto field_bc = field_bcs[blkid][face.bc_id];
                        pos = field_bc.other_pos(face);
                        Ai[k*nbands + iio] = field_bc.other_id(face);
                    } else {
                        other = face.left_cell;
                        if (other==cell) other = face.right_cell;
                        pos = other.pos[0];
                        Ai[k*nbands + iio] = other.id + block_offsets[blkid];
                    }

                    nx[io] = sign*face.n.x.re;
                    ny[io] = sign*face.n.y.re;
                    dx[io] = pos.x.re - cell.pos[0].x.re;
                    dy[io] = pos.y.re - cell.pos[0].y.re;
                }

                double dxN = dx[0]; double dyN = dy[0];
                double dxE = dx[1]; double dyE = dy[1];
                double dxS = dx[2]; double dyS = dy[2];
                double dxW = dx[3]; double dyW = dy[3];

                double nxN = nx[0]; double nyN = ny[0];
                double nxE = nx[1]; double nyE = ny[1];
                double nxS = nx[2]; double nyS = ny[2];
                double nxW = nx[3]; double nyW = ny[3];

                double[4] fdx, fdy, fdxx, fdyy;
                double _Ix, _Iy, _Ixx, _Iyy, D;

                // Figure out what kind of cell we are, since ones with ZeroNormalGradient
                // have different equations for the derivatives...
				uint celltype = ZNG_interior;
                foreach(io, face; cell.iface){
					if (face.is_on_boundary) {
						auto field_bc = field_bcs[blkid][face.bc_id];
						if ((cast(ZeroNormalGradient) field_bc) !is null) {
                            celltype = celltype | ZNG_types[io];
						}
                    }
				}

                switch (celltype) {
                case ZNG_north:
                    D = mixin(ZGN_D);
                    fdx[0] = mixin(ZGN_Nx);
                    fdx[1] = mixin(ZGN_Ex);
                    fdx[2] = mixin(ZGN_Sx);
                    fdx[3] = mixin(ZGN_Wx);
                    _Ix    = mixin(ZGN_Ix);

                    fdy[0] = mixin(ZGN_Ny);
                    fdy[1] = mixin(ZGN_Ey);
                    fdy[2] = mixin(ZGN_Sy);
                    fdy[3] = mixin(ZGN_Wy);
                    _Iy    = mixin(ZGN_Iy);

                    break;
                case ZNG_east:
                    D = mixin(ZGE_D);
                    fdx[0] = mixin(ZGE_Nx);
                    fdx[1] = mixin(ZGE_Ex);
                    fdx[2] = mixin(ZGE_Sx);
                    fdx[3] = mixin(ZGE_Wx);
                    _Ix    = mixin(ZGE_Ix);

                    fdy[0] = mixin(ZGE_Ny);
                    fdy[1] = mixin(ZGE_Ey);
                    fdy[2] = mixin(ZGE_Sy);
                    fdy[3] = mixin(ZGE_Wy);
                    _Iy    = mixin(ZGE_Iy);

                    break;
                case ZNG_south:
                    D = mixin(ZGS_D);
                    fdx[0] = mixin(ZGS_Nx);
                    fdx[1] = mixin(ZGS_Ex);
                    fdx[2] = mixin(ZGS_Sx);
                    fdx[3] = mixin(ZGS_Wx);
                    _Ix    = mixin(ZGS_Ix);

                    fdy[0] = mixin(ZGS_Ny);
                    fdy[1] = mixin(ZGS_Ey);
                    fdy[2] = mixin(ZGS_Sy);
                    fdy[3] = mixin(ZGS_Wy);
                    _Iy    = mixin(ZGS_Iy);

                    break;
                case ZNG_west:
                    D = mixin(ZGW_D);
                    fdx[0] = mixin(ZGW_Nx);
                    fdx[1] = mixin(ZGW_Ex);
                    fdx[2] = mixin(ZGW_Sx);
                    fdx[3] = mixin(ZGW_Wx);
                    _Ix    = mixin(ZGW_Ix);

                    fdy[0] = mixin(ZGW_Ny);
                    fdy[1] = mixin(ZGW_Ey);
                    fdy[2] = mixin(ZGW_Sy);
                    fdy[3] = mixin(ZGW_Wy);
                    _Iy    = mixin(ZGW_Iy);

                    break;
                case ZNG_interior:
                    D = mixin(R_D);
                    fdx[0] = mixin(R_Nx);
                    fdx[1] = mixin(R_Ex);
                    fdx[2] = mixin(R_Sx);
                    fdx[3] = mixin(R_Wx);
                    _Ix    = mixin(R_Ix);

                    fdy[0] = mixin(R_Ny);
                    fdy[1] = mixin(R_Ey);
                    fdy[2] = mixin(R_Sy);
                    fdy[3] = mixin(R_Wy);
                    _Iy    = mixin(R_Iy);

                    break;
                default:
                    string errMsg = format("An invalid ZNGtype '%s' was requested.", celltype);
                    throw new Error(errMsg);
                }

                foreach(io, face; cell.iface){
                    int iio = (io>1) ? to!int(io+1) : to!int(io); // -> [0,1,3,4] since 2 is 
                    face.fs.gas.sigma = conductivity(face.fs.gas, face.pos, gmodel); // TODO: Redundant work.
                    double sign = cell.outsign[io];
                    double S = face.length.re;
                    double sigmaF = face.fs.gas.sigma.re;
                    double dxF = face.pos.x.re - cell.pos[0].x.re;
                    double dyF = face.pos.y.re - cell.pos[0].y.re;
                    double nxF = sign*face.n.x.re;
                    double nyF = sign*face.n.y.re;
                    double emag = sqrt(dx[io]*dx[io] + dy[io]*dy[io]);
                    double ehatx = dx[io]/emag;
                    double ehaty = dy[io]/emag;

                    // Hybrid method
                    double facx = nxF - ehatx*ehatx*nxF - ehatx*ehaty*nyF;
                    double facy = nyF - ehaty*ehatx*nxF - ehaty*ehaty*nyF;
                    double fac = (ehatx*nxF + ehaty*nyF)/emag;

                    // Finite difference stencil
                    //    double facx = nxF;
                    //    double facy = nyF;
                    //    double fac = 0.0;
                    //// Direct normal gradient only
                    //    double facx = 0.0;
                    //    double facy = 0.0;
                    //    double fac = (ehatx*nxF + ehaty*nyF)/emag;

                    // PART ONE: Using the hybrid stencil, we first have a direct contribution to the
                    // flux, based on the cell and the other cell on the other side of the face
                    if (face.is_on_boundary) {
                        auto field_bc = field_bcs[blkid][face.bc_id];

                        // For a ZNG boundary, we want to use the stencil gradients for face i's
                        // contribution to the fluxes. It's ugly, but works for the moment
                        if ((cast(ZeroNormalGradient) field_bc) !is null){
                            facx = nxF;
                            facy = nyF;
                            fac = 0.0;
                        } else {
                            A[k*nbands + 2]  += field_bc.lhs_direct_component(fac, face);
                            A[k*nbands + iio]+= field_bc.lhs_other_component(fac, face);
                            b[k]             -= field_bc.rhs_direct_component(sign, fac, face);
                        }
                    } else {
                        A[k*nbands + 2] +=  -1.0*S*fac*sigmaF;
                        A[k*nbands + iio]+=  1.0*S*fac*sigmaF;
                    }

                    // PART TWO: The other part of the gradient comes from a finite difference stencil,
                    // which has components from all of the nearby cells, and cell k:
                    A[k*nbands + 2] +=  S/D*sigmaF*(facx*(_Ix) + facy*(_Iy));

                    // Each jface makes a contribution to the flux through "face"
                    foreach(jo, jface; cell.iface){
                        int jjo = (jo>1) ? to!int(jo+1) : to!int(jo); // -> [0,1,3,4] since 2 is the entry for "cell"
                        if (jface.is_on_boundary) {
                            auto field_bc = field_bcs[blkid][jface.bc_id];
                            A[k*nbands + jjo] += S*sigmaF*field_bc.lhs_stencil_component(D, facx, facy, fdx[jo], fdy[jo], jface);
                            b[k]              -= S*sigmaF*field_bc.rhs_stencil_component(D, facx, facy, fdx[jo], fdy[jo], jface);
                        } else {
                            A[k*nbands + jjo] += S/D*sigmaF*(facx*(fdx[jo])
                                                           + facy*(fdy[jo]));
                        }
                    }
                }
            }
        }

        // The Jacobi Preconditioner, a very simple scheme. Good for diagonally dominant matrices
        if (precondition){
            foreach(blkid, block; localFluidBlocks){
                foreach(cell; block.cells){
                    int k = cell.id + block_offsets[blkid];
                    //writefln(" A[i,:]=[%e,%e,%e,%e,%e] b=[%e]", A[k*nbands+0], A[k*nbands+1], A[k*nbands+2], A[k*nbands+3], A[k*nbands+4], b[k]);
                    double Akk = A[k*nbands + 2];
                    b[k] /= Akk;
                    foreach(iio; 0 .. nbands) A[k*nbands + iio] /= Akk;
                }
            }
        }

        phi0[] = 0.0;
        gmres.solve(N, nbands, A, Ai, b, phi0, phi, max_iter, verbose);

        // Unpack the solution into the "electric_potential" members stored in the cells
        size_t i = 0;
        foreach(block; localFluidBlocks){
            foreach(cell; block.cells){
                cell.electric_potential = phi[i];
                i += 1;
            }
        }

        // We also have enough info to set the ghost cell potentials too
        version(mpi_parallel){ exchanger.update_buffers(phi);}

        foreach(blkid, block; localFluidBlocks){
            foreach(j, bc; block.bc){
                auto field_bc = field_bcs[blkid][j];
                if (!field_bc.isShared) continue;

                foreach(fidx, f; bc.faces){
                    // Idx is the real cell's location in the phi array
                    int idx = field_bc.other_id(f);

                    double phiidx;
                    version(mpi_parallel) {
                        phiidx = exchanger.external_cell_buffer[idx-phi.length];
                    } else {
                        phiidx = phi[idx];
                    }

                    // Figure out which side the ghost cell is on and set 
                    if (bc.outsigns[fidx] == 1) {
                        f.right_cell.electric_potential = phiidx;
                    } else {
                        f.left_cell.electric_potential = phiidx;
                    }
                }
            }
        }
        return;
    }

    void compute_electric_field_vector(FluidBlock[] localFluidBlocks) {
    /*
        With the electric potential field solved for, compute its gradients
        and store the electric field vector.

        Notes:
         - TODO: MPI

        @author: Nick Gibbons
    */

        FluidFVCell other;
        foreach(blkid, block; localFluidBlocks){
            foreach(cell; block.cells){
                double[4] dx, dy, nx, ny, phis;

                foreach(io, face; cell.iface){
                    Vector3 pos;
                    double phi;
                    double sign = cell.outsign[io];
                    if (face.is_on_boundary) {
                        auto field_bc = field_bcs[blkid][face.bc_id];
                        phi = field_bc.phif(face);
                        pos = field_bc.other_pos(face);
                    } else {
                        other = face.left_cell;
                        if (other==cell) other = face.right_cell;
                        pos = other.pos[0];
                        phi = other.electric_potential;
                    }
                    nx[io] = sign*face.n.x.re;
                    ny[io] = sign*face.n.y.re;
                    dx[io] = pos.x.re - cell.pos[0].x.re;
                    dy[io] = pos.y.re - cell.pos[0].y.re;
                    phis[io] = phi;
                }

                double dxN = dx[0]; double dyN = dy[0]; double nxN = nx[0]; double nyN = ny[0];
                double dxE = dx[1]; double dyE = dy[1]; double nxE = nx[1]; double nyE = ny[1];
                double dxS = dx[2]; double dyS = dy[2]; double nxS = nx[2]; double nyS = ny[2];
                double dxW = dx[3]; double dyW = dy[3]; double nxW = nx[3]; double nyW = ny[3];

                double[4] fdx, fdy;
                double _Ix, _Iy, D;

                // Figure out what kind of cell we are, since ones with ZeroNormalGradient
                // have different equations for the derivatives...
				uint celltype = ZNG_interior;
                foreach(io, face; cell.iface){
					if (face.is_on_boundary) {
						auto field_bc = field_bcs[blkid][face.bc_id];
						if ((cast(ZeroNormalGradient) field_bc) !is null) {
                            celltype = celltype | ZNG_types[io];
						}
                    }
				}

                switch (celltype) {
                case ZNG_north:
                    D = mixin(ZGN_D);
                    fdx[0] = mixin(ZGN_Nx);
                    fdx[1] = mixin(ZGN_Ex);
                    fdx[2] = mixin(ZGN_Sx);
                    fdx[3] = mixin(ZGN_Wx);
                    _Ix    = mixin(ZGN_Ix);

                    fdy[0] = mixin(ZGN_Ny);
                    fdy[1] = mixin(ZGN_Ey);
                    fdy[2] = mixin(ZGN_Sy);
                    fdy[3] = mixin(ZGN_Wy);
                    _Iy    = mixin(ZGN_Iy);

                    break;
                case ZNG_east:
                    D = mixin(ZGE_D);
                    fdx[0] = mixin(ZGE_Nx);
                    fdx[1] = mixin(ZGE_Ex);
                    fdx[2] = mixin(ZGE_Sx);
                    fdx[3] = mixin(ZGE_Wx);
                    _Ix    = mixin(ZGE_Ix);

                    fdy[0] = mixin(ZGE_Ny);
                    fdy[1] = mixin(ZGE_Ey);
                    fdy[2] = mixin(ZGE_Sy);
                    fdy[3] = mixin(ZGE_Wy);
                    _Iy    = mixin(ZGE_Iy);

                    break;
                case ZNG_south:
                    D = mixin(ZGS_D);
                    fdx[0] = mixin(ZGS_Nx);
                    fdx[1] = mixin(ZGS_Ex);
                    fdx[2] = mixin(ZGS_Sx);
                    fdx[3] = mixin(ZGS_Wx);
                    _Ix    = mixin(ZGS_Ix);

                    fdy[0] = mixin(ZGS_Ny);
                    fdy[1] = mixin(ZGS_Ey);
                    fdy[2] = mixin(ZGS_Sy);
                    fdy[3] = mixin(ZGS_Wy);
                    _Iy    = mixin(ZGS_Iy);

                    break;
                case ZNG_west:
                    D = mixin(ZGW_D);
                    fdx[0] = mixin(ZGW_Nx);
                    fdx[1] = mixin(ZGW_Ex);
                    fdx[2] = mixin(ZGW_Sx);
                    fdx[3] = mixin(ZGW_Wx);
                    _Ix    = mixin(ZGW_Ix);

                    fdy[0] = mixin(ZGW_Ny);
                    fdy[1] = mixin(ZGW_Ey);
                    fdy[2] = mixin(ZGW_Sy);
                    fdy[3] = mixin(ZGW_Wy);
                    _Iy    = mixin(ZGW_Iy);

                    break;
                case ZNG_interior:
                    D = mixin(R_D);
                    fdx[0] = mixin(R_Nx);
                    fdx[1] = mixin(R_Ex);
                    fdx[2] = mixin(R_Sx);
                    fdx[3] = mixin(R_Wx);
                    _Ix    = mixin(R_Ix);

                    fdy[0] = mixin(R_Ny);
                    fdy[1] = mixin(R_Ey);
                    fdy[2] = mixin(R_Sy);
                    fdy[3] = mixin(R_Wy);
                    _Iy    = mixin(R_Iy);

                    break;
                default:
                    string errMsg = format("An invalid ZNGtype '%s' was requested.", celltype);
                    throw new Error(errMsg);
                }

                double Ex = (_Ix*cell.electric_potential + fdx[0]*phis[0] + fdx[1]*phis[1] + fdx[2]*phis[2] + fdx[3]*phis[3])/D;
                double Ey = (_Iy*cell.electric_potential + fdy[0]*phis[0] + fdy[1]*phis[1] + fdy[2]*phis[2] + fdy[3]*phis[3])/D;
                cell.electric_field[0] = Ex;
                cell.electric_field[1] = Ey;
            }
        }
    }

    void compute_boundary_current(FluidBlock[] localFluidBlocks, ref double current_in, ref double current_out) {
    /*
        Loop over the boundaries of the domain and compute the total electrical current flowing in and out.
        We put the contributions of each face into different buckets depending on their sign, negative means
        current flow in and positive out.

        Notes:
         - TODO: MPI
         - Caution: We assume that the field and conductivity are already set

        @author: Nick Gibbons
    */

        double Iin = 0.0;
        double Iout = 0.0;

        foreach(blkid, block; localFluidBlocks){
            foreach(cell; block.cells){
                foreach(io, face; cell.iface){
                    if (!face.is_on_boundary) continue;

                    auto field_bc = field_bcs[blkid][face.bc_id];
                    double phif = field_bc.phif(face);
                    if (phif==0.0) continue; // FIXME: We don't really need this
                    if (field_bc.isShared) continue;

                    double sign = cell.outsign[io];
                    double S = face.length.re;
                    double sigmaF = face.fs.gas.sigma.re;
                    double dxF = face.pos.x.re - cell.pos[0].x.re;
                    double dyF = face.pos.y.re - cell.pos[0].y.re;
                    double nxF = sign*face.n.x.re;
                    double nyF = sign*face.n.y.re;
                    double emag = sqrt(dxF*dxF + dyF*dyF);
                    double ehatx = dxF/emag;
                    double ehaty = dyF/emag;

                    // Hybrid method
                    double facx = nxF - ehatx*ehatx*nxF - ehatx*ehaty*nyF;
                    double facy = nyF - ehaty*ehatx*nxF - ehaty*ehaty*nyF;
                    double fac = (ehatx*nxF + ehaty*nyF)/emag;

                    double phigrad_dot_n = nxF*cell.electric_field[0] + nyF*cell.electric_field[1];

                    //double phigrad_dot_n = facx*cell.electric_field[0] + 
                    //                       facy*cell.electric_field[1] + 
                    //                       fac*(phif - cell.electric_potential);

                    double I = phigrad_dot_n*S*sigmaF;
                    if (I<0.0) {
                        Iin -= I;
                    } else if (I>0.0) {
                        Iout += I;
                    }
                }
            }
        }
        version(mpi_parallel){
            MPI_Allreduce(MPI_IN_PLACE, &Iin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Iout, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        current_in = Iin;
        current_out = Iout;
	}

    void compute_boundary_current_old(FluidBlock[] localFluidBlocks) {
    /*
        Loop over the boundaries of the domain and compute the total electrical current flowing in and out.
        We put the contributions of each face into different buckets depending on their sign, negative means
        current flow in and positive out.

        Notes:
         - TODO: MPI

        @author: Nick Gibbons
    */
        double Iin = 0.0;
        double Iout = 0.0;

        writeln("Called field.compute_boundary_current() ...");
        foreach(blkid, block; localFluidBlocks){
            auto gmodel = block.myConfig.gmodel;
            foreach(cell; block.cells){
                foreach(io, face; cell.iface){
                    if (face.is_on_boundary) {
                        face.fs.gas.sigma = conductivity(face.fs.gas, face.pos, gmodel);
                        double sign = cell.outsign[io];
                        auto field_bc = field_bcs[blkid][face.bc_id];
                        double I = field_bc.compute_current(sign, face, cell);
                        if (I<0.0) {
                            Iin -= I;
                        } else if (I>0.0) {
                            Iout += I;
                        }
                    }
                }
            }
        }
        writefln("    Current in:  %f (A/m)", Iin);
        writefln("    Current out: %f (A/m)", Iout);
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

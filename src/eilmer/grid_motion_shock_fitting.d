// grid_motion_shock_fitting.d
// Module for implementing shock-fitting with a moving grid in Eilmer4.
//
// 2015-Nov Kyle D. original implmentation (moving grid, shock fitting).
// 2019 Lachlan Whyborn multiple blocks and MPI
// 2021-02-15 PJ complete rework to have a more "global" view of the shock.

module grid_motion_shock_fitting;

import std.stdio;
import std.conv;
import std.math;
import std.algorithm;
import nm.complex;
import nm.number;
version(mpi_parallel) {
    import mpi;
    import mpi.util;
}

import globalconfig;
import globaldata;
import flowstate;
import fvcore;
import fvvertex;
import fvinterface;
import fvcell;
import onedinterp;
import bc;
import fluidblock;
import fluidblockarray;
import sfluidblock;
import geom;
import grid_motion;
import bc;
version(mpi_parallel) {
    import bc.ghost_cell_effect.full_face_copy : MPI_Wait_a_while, make_mpi_tag;
}


void compute_vtx_velocities_for_sf(FBArray fba)
{
    // The shock-fitting is coordinated by the FBArray which has a global view
    // of the shock boundary.  The boundary may be composed of several blocks
    // which, for a shared memory simulation, will all be in local memory.
    //
    // The general plan is to do the calculation in phases:
    // 1. Use the conservation equations to compute wave speed estimates
    //    at all faces on the shock boundary.
    // 2. Compute vertex velocities for shock-boundary locations, aligning those
    //    velocities with the "rails" along which we would like the shock
    //    vertices to move.
    // 3. Propagate velocities across all vertices in the blocks, decreasing
    //    the magnitude of the velocity as the east-boundary of the FBArray
    //    is approached.
    //
    // In an MPI context, it may be that the blocks are in different MPI tasks.
    // The design is that all tasks keep a global copy of the FBArray data
    // and use their local FluidBlocks to fill in the appropriate sections
    // of the global arrays.  Messages are then passed between the MPI tasks
    // to synchronize the content of the global arrays. The devil is in the details.
    //
    // Start by computing wave speeds and rail directions at all west-most faces.
    //
    int blkId = fba.blockArray[0][0][0];
    auto blk = cast(SFluidBlock) globalBlocks[blkId];
    auto bc = cast(BFE_ConstFlux) blk.bc[Face.west].postConvFluxAction[0];
    if (!bc) { throw new Error("Did not find an appropriate boundary-face effect."); }
    auto inflow = bc.fstate;
    bool allow_reconstruction = GlobalConfig.shock_fitting_allow_flow_reconstruction;
    // Work across all west-most blocks in the array, storing the wave speeds at face centres.
    foreach (jb; 0 .. fba.njb) {
        int j0 = 0; if (jb > 0) { foreach(j; 0 .. jb) { j0 += fba.njcs[j]; } }
        foreach (kb; 0 .. fba.nkb) {
            int k0 = 0; if (kb > 0) { foreach(k; 0 .. kb) { k0 += fba.nkcs[k]; } }
            blkId = fba.blockArray[0][jb][kb];
            if (canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                blk = cast(SFluidBlock) globalBlocks[blkId];
                FlowState Rght = new FlowState(blk.myConfig.gmodel);
                foreach (k; 0 .. blk.nkc) {
                    foreach (j; 0 .. blk.njc) {
                        auto f = blk.get_ifi(0,j,k);
                        if (allow_reconstruction) {
                            // Reconstruct the flow state just behind the shock from
                            // the flow states in the first two cells.
                            // Note that the actual order of reconstruction will be
                            // determined by GlocalConfig.interpolation_order.
                            FVCell cR0 = blk.get_cell(0,j,k);
                            FVCell cR1 = blk.get_cell(1,j,k);
                            Rght.copy_values_from(cR0.fs);
                            blk.one_d.interp_l0r2(f, cR0, cR1, cR0.iLength, cR1.iLength, Rght);
                            fba.face_ws[j0+j][k0+k] = wave_speed(inflow, Rght, f.n);
                        } else {
                            // Using the first cell-centre state for R0 is first-order.
                            fba.face_ws[j0+j][k0+k] = wave_speed(inflow, blk.get_cell(0,j,k).fs, f.n);
                        }
                        fba.face_pos[j0+j][k0+k] = f.pos;
                    }
                }
                // We need the vertex positions to do upwinding along the shock boundary.
                foreach (k; 0 .. blk.nkv) {
                    foreach (j; 0 .. blk.njv) {
                        fba.vtx_pos[j0+j][k0+k] = blk.get_vtx(0,j,k).pos[0];
                    }
                }
            } // end if canFind
        } // end foreach kb
    } // end foreach jb
    //
    version(mpi_parallel) {
        // MPI-parallel flavour.
        // MPI synchronization of fba.face_ws, fba.face_pos and fba.vtx_pos
        foreach (jb; 0 .. fba.njb) {
            int j0 = 0; if (jb > 0) { foreach(j; 0 .. jb) { j0 += fba.njcs[j]; } }
            foreach (kb; 0 .. fba.nkb) {
                int k0 = 0; if (kb > 0) { foreach(k; 0 .. kb) { k0 += fba.nkcs[k]; } }
                blkId = fba.blockArray[0][jb][kb];
                blk = cast(SFluidBlock) globalBlocks[blkId];
                int src_task = GlobalConfig.mpi_rank_for_block[blkId];
                int items;
                //
                // Broadcast the face velocities.
                if (canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    // Local task owns the block so we pack its data into the buffer.
                    assert(src_task == GlobalConfig.mpi_rank_for_local_task,
                           "Oops, source task should be local MPI task.");
                    items = 0;
                    foreach (k; 0 .. blk.nkc) {
                        foreach (j; 0 .. blk.njc) {
                            fba.buffer[items++] = fba.face_ws[j0+j][k0+k].re;
                        }
                    }
                } else {
                    // Local task does not own this block,
                    // but we need to know how many items are broadcast.
                    items = to!int(blk.nkc * blk.njc);
                }
                MPI_Bcast(fba.buffer.ptr, items, MPI_DOUBLE, src_task, fba.mpicomm);
                if (!canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    assert(src_task != GlobalConfig.mpi_rank_for_local_task,
                           "Oops, source task should not be local task.");
                    items = 0;
                    foreach (k; 0 .. blk.nkc) {
                        foreach (j; 0 .. blk.njc) {
                            fba.face_ws[j0+j][k0+k] = fba.buffer[items++];
                        }
                    }
                }
                //
                // Broadcast the face positions.
                if (canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    assert(src_task == GlobalConfig.mpi_rank_for_local_task,
                           "Oops, source task should be local MPI task.");
                    items = 0;
                    foreach (k; 0 .. blk.nkc) {
                        foreach (j; 0 .. blk.njc) {
                            fba.buffer[items++] = fba.face_pos[j0+j][k0+k].x.re;
                            fba.buffer[items++] = fba.face_pos[j0+j][k0+k].y.re;
                            fba.buffer[items++] = fba.face_pos[j0+j][k0+k].z.re;
                        }
                    }
                } else {
                    items = to!int(blk.nkc * blk.njc * 3);
                }
                MPI_Bcast(fba.buffer.ptr, items, MPI_DOUBLE, src_task, fba.mpicomm);
                if (!canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    assert(src_task != GlobalConfig.mpi_rank_for_local_task,
                           "Oops, source task should not be local task.");
                    items = 0;
                    foreach (k; 0 .. blk.nkc) {
                        foreach (j; 0 .. blk.njc) {
                            fba.face_pos[j0+j][k0+k].set(fba.buffer[items++],
                                                         fba.buffer[items++],
                                                         fba.buffer[items++]);
                        }
                    }
                }
                //
                // Broadcast the vertex positions.
                if (canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    assert(src_task == GlobalConfig.mpi_rank_for_local_task,
                           "Oops, source task should be local MPI task.");
                    items = 0;
                    foreach (k; 0 .. blk.nkv) {
                        foreach (j; 0 .. blk.njv) {
                            fba.buffer[items++] = fba.vtx_pos[j0+j][k0+k].x.re;
                            fba.buffer[items++] = fba.vtx_pos[j0+j][k0+k].y.re;
                            fba.buffer[items++] = fba.vtx_pos[j0+j][k0+k].z.re;
                        }
                    }
                } else {
                    items = to!int(blk.nkv * blk.njv * 3);
                }
                MPI_Bcast(fba.buffer.ptr, items, MPI_DOUBLE, src_task, fba.mpicomm);
                if (!canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    assert(src_task != GlobalConfig.mpi_rank_for_local_task,
                           "Oops, source task should not be local task.");
                    items = 0;
                    foreach (k; 0 .. blk.nkv) {
                        foreach (j; 0 .. blk.njv) {
                            fba.vtx_pos[j0+j][k0+k].set(fba.buffer[items++],
                                                        fba.buffer[items++],
                                                        fba.buffer[items++]);
                        }
                    }
                } // end if !canFind
            } // end foreach kb
        } // end foreach jb
    } // end version !mpi_parallel
    //
    // Compute rail directions for the vertices at the boundary.
    // The rails are presently (2021-02-16) fixed so there is lots of waste,
    // however, we expect to make the rails curve at a later date and so will
    // need to do these calculations based on the current boundary positions.
    //
    foreach (k; 0 .. fba.nkv) {
        foreach (j; 0 .. fba.njv) {
            Vector3 d = fba.p_east[j][k]; d -= fba.p_west[j][k]; d.normalize();
            fba.vtx_dir[j][k].set(d);
        }
    }
    //
    // Compute the shock-boundary vertex velocities with an upwind-weighting
    // of the face-centre velocities, as described by Ian.
    // Note that we want the vertex velocities that are aligned with the rails
    // so we start with the rail-direction vector and scale that
    // with the wave-speed estimate to get the vertex velocity.
    //
    if (GlobalConfig.dimensions == 2) {
        // The shock boundary is a line.
        int k = 0;
        // First vertex has only one face, so just use that velocity.
        Vector3 v = fba.vtx_dir[0][k]; v.scale(fba.face_ws[0][k]); fba.vtx_vel[0][k].set(v);
        // Do Ian's upwind weighting for vertices between faces.
        // I think that Ian used the post-shock flow properties but,
        // for the moment, use the free-stream properties in the Mach number weights.
        // Across the shock the tangential velocity will be unchanged,
        // however, we use the free-stream sound speed because it is readily
        // available whereas Ian used the post-shock sound speed.
        // Presumably, our Mach numbers will be higher once the shock has
        // fused with the boundary.
        foreach (j; 1 .. fba.njv-1) {
            Vector3 tA = fba.vtx_pos[j][k]; tA -= fba.face_pos[j-1][k]; tA.normalize();
            number MA = dot(inflow.vel, tA) / inflow.gas.a;
            number wA = Mach_weighting(MA);
            Vector3 tB = fba.vtx_pos[j][k]; tB -= fba.face_pos[j][k]; tB.normalize();
            number MB = dot(inflow.vel, tB) / inflow.gas.a;
            number wB = Mach_weighting(MB);
            number ws;
            if (fabs(wA+wB) > 0.0) {
                ws = (wA*fba.face_ws[j-1][k] + wB*fba.face_ws[j][k])/(wA+wB);
            } else {
                ws = 0.5*fba.face_ws[j-1][k] + 0.5*fba.face_ws[j][k];
            }
            v = fba.vtx_dir[j][k]; v.scale(ws); fba.vtx_vel[j][k].set(v);
        }
        v = fba.vtx_dir[$-1][k]; v.scale(fba.face_ws[$-1][k]); fba.vtx_vel[$-1][k].set(v);
    } else {
        // The shock boundary is a surface.
        throw new Error("Shock-velocity upwinding is not yet implemented in 3D.");
    }
    //
    // Propagate the shock-boundary velocities to the vertices
    // for all blocks in the array.
    //
    foreach (jb; 0 .. fba.njb) {
        int j0 = 0; if (jb > 0) { foreach(j; 0 .. jb) { j0 += fba.njcs[j]; } }
        foreach (kb; 0 .. fba.nkb) {
            int k0 = 0; if (kb > 0) { foreach(k; 0 .. kb) { k0 += fba.nkcs[k]; } }
            foreach (ib; 0 .. fba.nib) {
                blkId = fba.blockArray[ib][jb][kb];
                if (canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    blk = cast(SFluidBlock) globalBlocks[blkId];
                    int i0 = 0; if (ib > 0) { foreach(i; 0 .. ib) { i0 += fba.nics[i]; } }
                    foreach (k; 0 .. blk.nkv) {
                        foreach (j; 0 .. blk.njv) {
                            auto bndry_vel = fba.vtx_vel[j0+j][k0+k];
                            bndry_vel.scale(blk.myConfig.shock_fitting_scale_factor);
                            foreach (i; 0 .. blk.niv) {
                                auto vtx_vel = blk.get_vtx(i,j,k).vel;
                                vtx_vel[0].set(bndry_vel);
                                vtx_vel[0].scale(fba.velocity_weights[i0+i][j0+j][k0+k]);
                                // Note that we set only the first element of the velocity array.
                            }
                        }
                    }
                } // end if canFind
            } // end foreach ib
        } // end foreach kb
    } // end foreach jb
} // end compute_vtx_velocities_for_sf()


@nogc
number wave_speed(const(FlowState) L0, const(FlowState) R0, const(Vector3) n)
{
    // Compute wave speed at the mid-point of the boundary face
    // using the approach described in Ian Johnston's thesis.
    // See PJ workbook pages 32-35, 2019-11-09.
    //
    //  west boundary
    //       |
    //       +------+------+
    //   L0  |  R0  |  R1  |
    //       +------+------+
    //     n-->
    //
    // L0 is the free-stream
    // R0 is post-shock, presumably.
    // The face normal, n, is pointing into the domain (i.e. into cell R0).
    //
    number veln = dot(L0.vel, n);
    // Ian's shock detector looks for a significan density difference, equation 4.23
    immutable double kappa = 0.2;
    if ((R0.gas.rho - L0.gas.rho) > kappa*R0.gas.rho) {
        // Estimate shock-wave speed from conservation equations.
        // Conservation of mass, equation 4.5 in Ian's thesis.
        number ws1 = (L0.gas.rho*dot(L0.vel,n) - R0.gas.rho*dot(R0.vel,n)) / (L0.gas.rho - R0.gas.rho);
        // Conservation of momentum, equation 4.6 in Ian's thesis.
        double pRpL = R0.gas.p.re - L0.gas.p.re;
        number ws2 = veln - sgn(pRpL)/L0.gas.rho * sqrt(fabs(pRpL / (1.0/L0.gas.rho - 1.0/R0.gas.rho)));
        immutable double alpha = 0.5;
        return alpha*ws1 + (1.0-alpha)*ws2;
    } else {
        // Estimate shock-wave speed using local sound speed.
        return veln - L0.gas.a;
    }
} // end wave_speed()

@nogc
number Mach_weighting(number M)
{
    // Weighting function defined on page 78 of Ian's thesis.
    // M is mach number in direction pointing toward the vertex.
    number Mp1 = M + 1.0;
    return (M <= 1.0) ? (0.125*(Mp1^^2 + Mp1*fabs(Mp1))) : M;
}

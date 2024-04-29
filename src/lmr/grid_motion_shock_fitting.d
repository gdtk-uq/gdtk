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
import ntypes.complex;
import nm.number;
version(mpi_parallel) {
    import mpi;
}

import globalconfig;
import globaldata;
import flowstate;
import fvvertex;
import fvinterface;
import lmr.fluidfvcell;
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
    // In an MPI context, it may be that the blocks are in different MPI tasks.
    // The design is that all tasks keep a global copy of the FBArray data
    // and use their local FluidBlocks to fill in the appropriate sections
    // of the global arrays.  Messages are then passed between the MPI tasks
    // to synchronize the content of the global arrays. The devil is in the details.
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
    bool allow_reconstruction = GlobalConfig.shock_fitting_allow_flow_reconstruction;
    double filter_scale = GlobalConfig.shock_fitting_filter_velocity_scale;
    bool assume_symmetry = GlobalConfig.shock_fitting_assume_symmetry_at_first_point;
    int blkId = fba.blockArray[0][0][0];
    auto blk = cast(SFluidBlock) globalBlocks[blkId];
    auto bc = cast(BFE_ConstFlux) blk.bc[Face.west].postConvFluxAction[0];
    if (!bc) { throw new Error("Did not find an appropriate boundary-face effect."); }
    auto nominal_inflow = bc.fstate;
    auto inflow = bc.fstate.dup(); // Start with a copy that we will partially overwrite.
    SourceFlow sourceFlow;
    if (bc.r > 0.0) {
        if (GlobalConfig.dimensions == 2 && !GlobalConfig.axisymmetric) {
            throw new Error("In a 2D shock-fitted flow, conical inflow is only for axisymmetric flows.");
        }
        sourceFlow = bc.sflow;
    }
    //
    // Start by computing wave speeds and rail directions at all west-most faces.
    //
    // Work across all west-most blocks in the array, storing the wave speeds at face centres.
    foreach (jb; 0 .. fba.njb) {
        int j0 = 0; if (jb > 0) { foreach(j; 0 .. jb) { j0 += fba.njcs[j]; } }
        foreach (kb; 0 .. fba.nkb) {
            int k0 = 0; if (kb > 0) { foreach(k; 0 .. kb) { k0 += fba.nkcs[k]; } }
            blkId = fba.blockArray[0][jb][kb];
            if (canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                blk = cast(SFluidBlock) globalBlocks[blkId];
                FlowState Rght = FlowState(blk.myConfig.gmodel, blk.myConfig.turb_model.nturb);
                foreach (k; 0 .. blk.nkc) {
                    foreach (j; 0 .. blk.njc) {
                        auto f = blk.get_ifi(0,j,k);
                        inflow.gas.p = nominal_inflow.gas.p;
                        inflow.gas.rho = nominal_inflow.gas.rho;
                        inflow.gas.u = nominal_inflow.gas.u;
                        inflow.vel.x = nominal_inflow.vel.x;
                        inflow.vel.y = nominal_inflow.vel.y;
                        inflow.vel.z = nominal_inflow.vel.z;
                        if (bc.r > 0.0) {
                            // We want to adjust the inflow velocities to be conical.
                            double dx = f.pos.x.re - bc.x0;
                            double dy = f.pos.y.re - bc.y0;
                            double dz = f.pos.z.re - bc.z0;
                            double hypot = sqrt(dx*dx + dy*dy + dz*dz);
                            double[4] deltas = sourceFlow.get_rho_v_p_u_increments(hypot-bc.r);
                            inflow.gas.rho += deltas[0];
                            double v = nominal_inflow.vel.x.re + deltas[1];
                            inflow.vel.x = v * dx/hypot;
                            inflow.vel.y = v * dy/hypot;
                            inflow.vel.z = v * dz/hypot;
                            inflow.gas.p += deltas[2];
                            inflow.gas.u += deltas[3];
                        }
                        if (allow_reconstruction) {
                            // Reconstruct the flow state just behind the shock from
                            // the flow states in the first two cells.
                            // Note that the actual order of reconstruction will be
                            // determined by GlobalConfig.interpolation_order.
                            FluidFVCell cR0 = blk.get_cell(0,j,k);
                            FluidFVCell cR1 = blk.get_cell(1,j,k);
                            Rght.copy_values_from(cR0.fs);
                            // blk.one_d.interp_l0r2(f, cR0, cR1, cR0.iLength, cR1.iLength, Rght);
                            blk.one_d.interp_l0r2(f, Rght, Rght); // note that Rght is used as a dummy variable here
                            fba.face_ws[j0+j][k0+k] = wave_speed(inflow, Rght, f.n);
                        } else {
                            // Using the first cell-centre state for R0 is first-order.
                            fba.face_ws[j0+j][k0+k] = wave_speed(inflow, *(blk.get_cell(0,j,k).fs), f.n);
                        }
                        fba.face_pos[j0+j][k0+k] = f.pos;
                        fba.face_a[j0+j][k0+k] = blk.get_cell(0,j,k).fs.gas.a;
                    }
                }
                // We need the vertex positions to do upwinding along the shock boundary
                // and also to compute the counter-kink velocities.
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
                // Broadcast the face velocities and post-shock sound-speeds.
                if (canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    // Local MPI task owns the block so we pack its data into the buffer
                    // for communication to all other MPI tasks.
                    assert(src_task == GlobalConfig.mpi_rank_for_local_task,
                           "Oops, source task should be local MPI task.");
                    items = 0;
                    foreach (k; 0 .. blk.nkc) {
                        foreach (j; 0 .. blk.njc) {
                            fba.buffer[items++] = fba.face_ws[j0+j][k0+k].re;
                            fba.buffer[items++] = fba.face_a[j0+j][k0+k].re;
                        }
                    }
                } else {
                    // Local task does not own this block,
                    // but we need to know how many items are broadcast.
                    items = to!int(blk.nkc * blk.njc * 2);
                }
                MPI_Bcast(fba.buffer.ptr, items, MPI_DOUBLE, src_task, fba.mpicomm);
                if (!canFind(GlobalConfig.localFluidBlockIds, blkId)) {
                    // The local MPI task does not own this block so get the data
                    // from the buffer and copy it into the local fba object.
                    assert(src_task != GlobalConfig.mpi_rank_for_local_task,
                           "Oops, source task should not be local task.");
                    items = 0;
                    foreach (k; 0 .. blk.nkc) {
                        foreach (j; 0 .. blk.njc) {
                            fba.face_ws[j0+j][k0+k] = fba.buffer[items++];
                            fba.face_a[j0+j][k0+k] = fba.buffer[items++];
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
    // The rails point "in", toward the solid body that is presumed to be
    // at the east face of the superblock.
    //
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
    // of the face-centre velocities, as described by Ian Johnston in his thesis.
    // Note that we want the vertex velocities that are aligned with the rails
    // so we start with the rail-direction vector and scale that
    // with the wave-speed estimate to get the vertex velocity.
    //
    if (GlobalConfig.dimensions == 2) {
        // The shock boundary is a line in 2D.
        int k = 0;
        // First vertex has only one face, so just use that velocity.
        Vector3 v = fba.vtx_dir[0][k]; v.scale(fba.face_ws[0][k]); fba.vtx_vel[0][k].set(v);
        // Do Ian's upwind weighting for vertices between faces.
        // I think that Ian used the post-shock flow properties.
        // Across the shock the tangential velocity will be unchanged, so we use that.
        foreach (j; 1 .. fba.njv-1) {
            Vector3 tA = fba.vtx_pos[j][k]; tA -= fba.face_pos[j-1][k]; tA.normalize();
            number MA = dot(nominal_inflow.vel, tA) / fba.face_a[j-1][k];
            number wA = Mach_weighting(MA);
            Vector3 tB = fba.vtx_pos[j][k]; tB -= fba.face_pos[j][k]; tB.normalize();
            number MB = dot(nominal_inflow.vel, tB) / fba.face_a[j][k];
            number wB = Mach_weighting(MB);
            number ws;
            if (fabs(wA+wB) > 0.0) {
                ws = (wA*fba.face_ws[j-1][k] + wB*fba.face_ws[j][k])/(wA+wB);
            } else {
                // Would not expect to have both weightings zero but, just in case.
                ws = 0.5*fba.face_ws[j-1][k] + 0.5*fba.face_ws[j][k];
            }
            v = fba.vtx_dir[j][k]; v.scale(ws); fba.vtx_vel[j][k].set(v);
        }
        v = fba.vtx_dir[$-1][k]; v.scale(fba.face_ws[$-1][k]); fba.vtx_vel[$-1][k].set(v);
        if (filter_scale > 0.0) {
            // Work along the shock and compute additional velocity components
            // to counter kinks in the shock.
            // We work with a stencil of 4 vertices at a time.
            // For details, see PJ's workbook page 37, 2021-08-09.
            foreach (j; 0 .. fba.njv-3) {
                Vector3 p0 = fba.vtx_pos[j][k];
                Vector3 p1 = fba.vtx_pos[j+1][k];
                Vector3 p2 = fba.vtx_pos[j+2][k];
                Vector3 p3 = fba.vtx_pos[j+3][k];
                Vector3 pmid = p1; pmid.scale(0.5); pmid.add(p2, 0.5);
                number L = distance_between(p1, p2);
                //
                // Compute displacements of the two interior points.
                Vector3 dp1 = pmid; dp1.scale(2.0/3.0); dp1.add(p0, 1.0/3.0); dp1 -= p1;
                Vector3 dp2 = pmid; dp2.scale(2.0/3.0); dp2.add(p3, 1.0/3.0); dp2 -= p2;
                if (dot(dp1, dp2) < 0.0) {
                    // Compute counter-kink velocity only for a corrugation,
                    // not in a situation that would just flatten the shock.
                    //
                    // Start with rail direction and scale with local sound speed.
                    Vector3 v_inc = fba.vtx_dir[j+1][k];
                    v_inc.scale(dot(v_inc, dp1)*filter_scale*inflow.gas.a/L);
                    fba.vtx_vel[j+1][k].add(v_inc);
                    //
                    v_inc = fba.vtx_dir[j+2][k];
                    v_inc.scale(dot(v_inc, dp2)*filter_scale*inflow.gas.a/L);
                    fba.vtx_vel[j+2][k].add(v_inc);
                }
            }
            if (assume_symmetry) {
                // Treat the symmetry vertex.
                // Try to fit a quadratic curve to the first few points on the shock
                // and use the x-location of that curve as the ideal point for the shock.
                // For details, see PJ's workbook pages 34,35, 2021-08-05.
                size_t N = 4;
                number x0 = fba.vtx_pos[0][k].x; number y0 = fba.vtx_pos[0][k].y;
                number Sx = 0.0; number Sxy2 = 0.0; number Sy2 = 0.0; number Sy4 = 0.0;
                foreach (j; 1 .. N+1) {
                    number x = fba.vtx_pos[j][k].x;
                    number y = fba.vtx_pos[j][k].y;
                    number y2 = (y-y0)^^2;
                    Sx += x; Sxy2 += x * y2; Sy2 += y2; Sy4 += y2^^2;
                }
                number denom = Sy2^^2 - N*Sy4;
                if (fabs(denom) > 0.0) {
                    number xstar = (Sxy2*Sy2 - Sx*Sy4)/denom;
                    number L = distance_between(fba.vtx_pos[0][k], fba.vtx_pos[1][k]);
                    Vector3 v_inc = fba.vtx_dir[0][k];
                    v_inc.scale((xstar - x0)*filter_scale*inflow.gas.a/L);
                    fba.vtx_vel[0][k].add(v_inc);
                } // else, we cannot solve for xstar so leave the velocity untouched.
            }
        }
    } else {
        // The shock boundary is a surface in 3D.
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

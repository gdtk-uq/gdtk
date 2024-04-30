/**
 * Eilmer4 boundary interface loads writing functions.
 *
 * Author: Kyle Damm
 * First code: 2017-03-17
 * Edits by Will Landsberg and Tim Cullen
 *
 * 2018-01-20 Some code clean up by PJ;
 * [TODO] more is needed for MPI context.
 *
 * History:
 *   2024-03-11 RJG update for eilmer 5
 */

module lmr.loads;

import std.stdio;
import std.string;
import std.file;
import std.array;
import std.algorithm;
import std.format;
import std.parallelism;
import std.math;
import std.conv;
import std.typecons : Tuple;
import std.json;
import ntypes.complex;
import nm.number;

import lmrconfig;
import json_helper;
import globalconfig;
import globaldata;
import fileutil;
import solidfvcell;
import solidblock;
import flowstate;
import flowgradients;
import geom;
import fluidblock;
import sfluidblock: SFluidBlock;
import ufluidblock: UFluidBlock;
import fvinterface;
import bc;
import mass_diffusion;
import conservedquantities;

version(mpi_parallel) {
    import mpi;
}

/**
 * Write loads when in timemarching mode.
 *
 * Authors: KAD, PAJ and RJG
 * Date: 2024-03-11
 */
void writeLoadsToFile(double sim_time, int current_loads_tindx) {
    foreach (blk; localFluidBlocks) {
        if (blk.active) {
            final switch (blk.grid_type) {
            case Grid_t.unstructured_grid:
                auto ublk = cast(UFluidBlock) blk;
                assert(ublk !is null, "Oops, this should be a UFluidBlock object.");
                apply_unstructured_grid(ublk, sim_time, current_loads_tindx);
                break;
            case Grid_t.structured_grid:
                auto sblk = cast(SFluidBlock) blk;
                assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                apply_structured_grid(sblk, sim_time, current_loads_tindx);
            } // end final switch
        } // if (blk.active)
    } // end foreach
} // end write_boundary_loads_to_file()


void apply_unstructured_grid(UFluidBlock blk, double sim_time, int current_loads_tindx)
{
    foreach (bndary; blk.bc) {
        if (canFind(GlobalConfig.group_names_for_loads, bndary.group)) {
            string fname = generate_boundary_load_file(blk.id, bndary.which_boundary, current_loads_tindx, sim_time, bndary.group);
            foreach (i, iface; bndary.faces) {
                // cell width normal to surface
                number w = (bndary.outsigns[i] == 1) ? iface.left_cell.L_min : iface.right_cell.L_min;
                computeAndStoreLoads(iface, bndary.outsigns[i], w, fname);
            }
        }
    }
} // end apply_unstructured_grid()

void apply_structured_grid(SFluidBlock blk, double sim_time, int current_loads_tindx) {
    foreach (bndary; blk.bc) {
        if (canFind(GlobalConfig.group_names_for_loads, bndary.group)) {
            string fname = generate_boundary_load_file(blk.id, bndary.which_boundary, current_loads_tindx, sim_time, bndary.group);
            // For structured blocks, all cell faces along a boundary point out or
            // all faces point in, so just use a constant for the outsign value.
            final switch (bndary.which_boundary) {
            case Face.north:
                foreach (k; 0 .. blk.nkc) {
                    foreach (i; 0 .. blk.nic) {
                        auto cell = blk.get_cell(i,blk.njc-1,k);
                        auto f = cell.iface[Face.north];
                        number cellWidthNormalToSurface = cell.jLength;
                        computeAndStoreLoads(f, 1, cellWidthNormalToSurface, fname);
                    }
                }
                break;
            case Face.east:
                foreach (k; 0 .. blk.nkc) {
                    foreach (j; 0 .. blk.njc) {
                        auto cell = blk.get_cell(blk.nic-1,j,k);
                        auto f = cell.iface[Face.east];
                        number cellWidthNormalToSurface = cell.iLength;
                        computeAndStoreLoads(f, 1, cellWidthNormalToSurface, fname);
                    }
                }
                break;
            case Face.south:
                foreach (k; 0 .. blk.nkc) {
                    foreach (i; 0 .. blk.nic) {
                        auto cell = blk.get_cell(i,0,k);
                        auto f = cell.iface[Face.south];
                        number cellWidthNormalToSurface = cell.jLength;
                        computeAndStoreLoads(f, -1, cellWidthNormalToSurface, fname);
                    }
                }
                break;
            case Face.west:
                foreach (k; 0 .. blk.nkc) {
                    foreach (j; 0 .. blk.njc) {
                        auto cell = blk.get_cell(0,j,k);
                        auto f = cell.iface[Face.west];
                        number cellWidthNormalToSurface = cell.iLength;
                        computeAndStoreLoads(f, -1, cellWidthNormalToSurface, fname);
                    }
                }
                break;
            case Face.top:
                foreach (i; 0 .. blk.nic) {
                    foreach (j; 0 .. blk.njc) {
                        auto cell = blk.get_cell(i,j,blk.nkc-1);
                        auto f = cell.iface[Face.top];
                        number cellWidthNormalToSurface = cell.kLength;
                        computeAndStoreLoads(f, 1, cellWidthNormalToSurface, fname);
                    }
                }
                break;
            case Face.bottom:
                foreach (i; 0 .. blk.nic) {
                    foreach (j; 0 .. blk.njc) {
                        auto cell = blk.get_cell(i,j,0);
                        auto f = cell.iface[Face.bottom];
                        number cellWidthNormalToSurface = cell.kLength;
                        computeAndStoreLoads(f, -1, cellWidthNormalToSurface, fname);
                    }
                }
                break;
            } // end switch which_boundary
        } // end if
    } // end foreach bndary
} // end apply_structured_grid()

string generate_boundary_load_file(int blkid, int bcid, int current_loads_tindx, double sim_time, string group)
{
    // generate data file -- naming format tindx_groupname.dat
    string fname = loadsFilename(current_loads_tindx, blkid, bcid, group);
    // create file
    auto f = File(fname, "w");
    f.writeln("pos.x pos.y pos.z n.x n.y n.z area cellWidthNormalToSurface outsign"~
              " p rho T vel.x vel.y vel.z mu a"~
              " Re y+"~
              " tau_wall.x tau_wall.y tau_wall.z"~
              " q_total q_cond q_diff");
    f.close();
    return fname;
} // end generate_boundary_load_file()

int count_written_loads()
{
    int nPrevLoadsWritten;
    string fname = lmrCfg.loadsDir ~ "/" ~ lmrCfg.loadsPrefix ~ "-metadata";
    if (exists(fname)) {
        auto lines = readText(fname).splitLines();
        if (lines.length == 1) {
            // looks like we found a single line.
            nPrevLoadsWritten = 0;
        } else {
            // assume we have a valid line to work with
            auto finalLine = lines[$-1];
            nPrevLoadsWritten = to!int(finalLine.split[0]) + 1;
        }
    } else {
        nPrevLoadsWritten = 0;
    }

    return nPrevLoadsWritten;
}

void init_current_loads_indx_dir(int current_loads_tindx)
{
    string dirName = lmrCfg.loadsDir ~ "/" ~ format(lmrCfg.snapshotIdxFmt, current_loads_tindx);
    ensure_directory_is_present(dirName);
}

void wait_for_current_indx_dir(int current_loads_tindx)
{
    string dirName = lmrCfg.loadsDir ~ "/" ~ format(lmrCfg.snapshotIdxFmt, current_loads_tindx);
    wait_for_directory_to_be_present(dirName);
}

void init_loads_metadata_file()
{
    string fname = lmrCfg.loadsDir ~ "/" ~ lmrCfg.loadsPrefix ~ "-metadata";
    std.file.write(fname, "loads_index    sim_time    step\n");
}

void update_loads_metadata_file(double sim_time, int step, int current_loads_tindx)
{
    string fname = lmrCfg.loadsDir ~ "/" ~ lmrCfg.loadsPrefix ~ "-metadata";
    std.file.append(fname, format("%04d %.18e %d\n", current_loads_tindx, sim_time, step));
}

void computeAndStoreLoads(FVInterface iface, int outsign, number cellWidthNormalToSurface, string fname)
{
    auto gmodel = GlobalConfig.gmodel_master;
    auto cqi = GlobalConfig.cqi;
    FlowState* fs = iface.fs;
    FlowGradients* grad = iface.grad;
    // iface orientation
    number nx = iface.n.x; number ny = iface.n.y; number nz = iface.n.z;
    // iface properties
    number mu_wall = fs.gas.mu; // laminar viscosity
    number P = fs.gas.p;
    number T_wall = fs.gas.T;
    number a_wall = fs.gas.a;
    number rho_wall = fs.gas.rho;
    number u_wall = fs.vel.x;
    number v_wall = fs.vel.y;
    number w_wall = fs.vel.z;
    number Re_wall, nu_wall, u_star, y_plus;
    number q_total, q_cond, q_diff;
    number tau_wall_x, tau_wall_y, tau_wall_z;
    if (GlobalConfig.viscous) {
        iface.F.clear();
        iface.average_turbulent_transprops();
        iface.viscous_flux_calc();
        // compute heat load
        q_cond = iface.q_conduction;
        q_diff = iface.q_diffusion;
        q_total = q_cond + q_diff;
        // compute stress tensor at interface in global reference frame
        iface.F.clear();
        iface.viscous_flux_calc();
        tau_wall_x = iface.F[cqi.xMom];
        tau_wall_y = iface.F[cqi.yMom];
        tau_wall_z = (cqi.threeD) ? iface.F[cqi.zMom] : to!number(0.0);
        number tau_wall = sqrt(tau_wall_x^^2 + tau_wall_y^^2 + tau_wall_z^^2);
        // compute y+
        nu_wall = mu_wall / rho_wall;
        u_star = sqrt(tau_wall / rho_wall);
        y_plus = u_star*(cellWidthNormalToSurface/2.0) / nu_wall;
        // compute cell Reynolds number
        Re_wall = rho_wall * a_wall * cellWidthNormalToSurface / mu_wall;
    } else {
        // For an inviscid simulation, we have only pressure.
        tau_wall_x = 0.0;
        tau_wall_y = 0.0;
        tau_wall_z = 0.0;
        q_total = 0.0;
        q_cond = 0.0;
        q_diff = 0.0;
        Re_wall = 0.0;
        y_plus = 0.0;
    }
    // store in file
    auto writer = format("%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %d %20.16e " ~
                         "%20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e %20.16e " ~
                         "%20.16e %20.16e %20.16e %20.16e %20.16e\n",
                         iface.pos.x.re, iface.pos.y.re, iface.pos.z.re, nx.re, ny.re, nz.re,
                         iface.area[0].re, cellWidthNormalToSurface.re, outsign,
                         P.re, rho_wall.re, T_wall.re, u_wall.re, v_wall.re, w_wall.re,
                         mu_wall.re, a_wall.re, Re_wall.re, y_plus.re,
                         tau_wall_x.re, tau_wall_y.re, tau_wall_z.re,
                         q_total.re, q_cond.re, q_diff.re);
    std.file.append(fname, writer);
}




struct RunTimeLoads {
    Tuple!(size_t,size_t)[] surfaces;
    Vector3 momentCentre;
    Vector3 resultantForce = Vector3(0.0, 0.0, 0.0);
    Vector3 resultantMoment = Vector3(0.0, 0.0, 0.0);
}

double[] groupedLoads;

void initRunTimeLoads(JSONValue jsonData)
{
    int nLoadsGroups = getJSONint(jsonData, "ngroups", 0);
    groupedLoads.length = nLoadsGroups*6; // 6 components: F.x, F.y, F.z, M.x, M.y, M.z
    runTimeLoads.length = nLoadsGroups;
    foreach (grpIdx; 0 .. nLoadsGroups) {
        auto jsonDataForGroup = jsonData["group-" ~ to!string(grpIdx)];
        string groupName = getJSONstring(jsonDataForGroup, "groupLabel", "");
        runTimeLoadsByName[groupName] = grpIdx;
        double x = getJSONdouble(jsonDataForGroup, "momentCtr_x", 0.0);
        double y = getJSONdouble(jsonDataForGroup, "momentCtr_y", 0.0);
        double z = getJSONdouble(jsonDataForGroup, "momentCtr_z", 0.0);
        runTimeLoads[grpIdx].momentCentre = Vector3(x, y, z);
        // Now look over all blocks in this process and configure the surfaces list
        foreach (iblk, blk; localFluidBlocks) {
            foreach (ibndry, bndry; blk.bc) {
                if (groupName == bndry.group) {
                    runTimeLoads[grpIdx].surfaces ~= Tuple!(size_t,size_t)(iblk, ibndry);
                }
            }
        }
    }

} // end initRunTimeLoads()

//@nogc
void computeRunTimeLoads()
{
    auto cqi = GlobalConfig.cqi;
    // Make sure each process has finished getting its viscous fluxes in place.
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    // On each process, tally up forces and moments by group
    foreach (grpIdx, ref group; runTimeLoads) {
        group.resultantForce.clear();
        group.resultantMoment.clear();
        foreach (blkNbndry; group.surfaces) {
            auto iblk = blkNbndry[0];
            auto ibndry = blkNbndry[1];
            foreach (iface, face; localFluidBlocks[iblk].bc[ibndry].faces) {
                int outsign = localFluidBlocks[iblk].bc[ibndry].outsigns[iface];
                number area = face.area[0] * outsign;
                number pdA = face.fs.gas.p * area;
                Vector3 pressureForce = Vector3(pdA*face.n.x, pdA*face.n.y, pdA*face.n.z);
                Vector3 momentArm; momentArm.set(face.pos); momentArm -= group.momentCentre;
                Vector3 momentContrib; cross(momentContrib, momentArm, pressureForce);
                group.resultantForce += pressureForce;
                group.resultantMoment += momentContrib;
                if (GlobalConfig.viscous) {
                    // The shear stresses have been calculated and stored in the flux
                    // for momentum just before calling this function.
                    Vector3 shearForce = Vector3(face.F[cqi.xMom]*area,
                                                 face.F[cqi.yMom]*area,
                                                 (cqi.threeD) ? face.F[cqi.zMom]*area : to!number(0.0));
                    cross(momentContrib, momentArm, shearForce);
                    group.resultantForce += shearForce;
                    group.resultantMoment += momentContrib;
                }
            }
        }
    }

    // The remainder is only required for parallel work.
    // We do the following:
    // 1. update the groupedLoads array locally on each process
    // 2. reduce the array with a sum across all processes
    // 3. make that sum available on all processes.
    // 4. place the array values back in our convenient data structure
    version(mpi_parallel) {
        version(nk_accelerator) {
            throw new Error("computeRunTimeLoads(): not available in e4-nk-dist.");
        }
        else {
            foreach (grpIdx, group; runTimeLoads) {
                // Place local values in groupedLoads on each process
                groupedLoads[grpIdx*6+0] = group.resultantForce.x.re;
                groupedLoads[grpIdx*6+1] = group.resultantForce.y.re;
                groupedLoads[grpIdx*6+2] = group.resultantForce.z.re;
                groupedLoads[grpIdx*6+3] = group.resultantMoment.x.re;
                groupedLoads[grpIdx*6+4] = group.resultantMoment.y.re;
                groupedLoads[grpIdx*6+5] = group.resultantMoment.z.re;
            }
            // Make sure each process has updated its loads calculation complete
            // before going on.
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, groupedLoads.ptr, to!int(groupedLoads.length), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            foreach (grpIdx, ref group; runTimeLoads) {
                group.resultantForce.x = groupedLoads[grpIdx*6+0];
                group.resultantForce.y = groupedLoads[grpIdx*6+1];
                group.resultantForce.z = groupedLoads[grpIdx*6+2];
                group.resultantMoment.x = groupedLoads[grpIdx*6+3];
                group.resultantMoment.y = groupedLoads[grpIdx*6+4];
                group.resultantMoment.z = groupedLoads[grpIdx*6+5];
            }
        } // END version(!nk_accelerator)
    } // end version(mpi_parallel)
}

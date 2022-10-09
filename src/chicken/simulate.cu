// simulate.cu
// Include file for chicken, high-level simulation functions.
//
// PJ 2022-09-09

#ifndef SIMULATE_INCLUDED
#define SIMULATE_INCLUDED

#include <cmath>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <limits>

#include "number.cu"
#include "vector3.cu"
#include "gas.cu"
#include "flow.cu"
#include "vertex.cu"
#include "face.cu"
#include "cell.cu"
#include "block.cu"
#include "config.cu"

using namespace std;

namespace SimState {
    number dt = 0.0;
    int step = 0;
    number t = 0.0;
    number t_plot = 0.0;
    int next_plot_indx = 1;
    int steps_since_last_plot = 0;
};

vector<Block> fluidBlocks;

__host__
void initialize_simulation(int tindx_start)
{
    char nameBuf[256];
    filesystem::path pth{Config::job};
    if (!filesystem::exists(pth) || !filesystem::is_directory(pth)) {
        throw runtime_error("Job directory is not present in current directory.");
    }
    read_config_file(Config::job + "/config.json");
    // Read initial grids and flow data
    for (int blk_id=0; blk_id < Config::nFluidBlocks; ++blk_id) {
        BConfig& blk_config = Config::blk_configs[blk_id];
        int i = blk_config.i; int j = blk_config.j; int k = blk_config.k;
        if (blk_id != Config::blk_ids[i][j][k]) {
            throw runtime_error("blk_id=" + to_string(blk_id) + " is inconsistent: i=" +
                                to_string(i) + " j=" + to_string(j) + " k=" + to_string(k));
        }
        Block blk;
        blk.configure(Config::nics[i], Config::njcs[j], Config::nkcs[k], blk_config.active);
        sprintf(nameBuf, "/grid/grid-%04d-%04d-%04d.gz", i, j, k);
        string fileName = Config::job + string(nameBuf);
        blk.readGrid(fileName);
        sprintf(nameBuf, "/flow/t%04d/flow-%04d-%04d-%04d.zip", tindx_start, i, j, k);
        fileName = Config::job + string(nameBuf);
        blk.readFlow(fileName);
        blk.computeGeometry();
        if (blk_config.active) { blk.encodeConserved(0); }
        fluidBlocks.push_back(blk);
    }
    //
    // Set up the simulation control parameters.
    SimState::t = 0.0;
    // Read times.data file to determine starting time.
    ifstream timesFile(Config::job+"/times.data", ifstream::binary);
    if (timesFile.good()) {
        string line;
        while (getline(timesFile, line)) {
            if (line.empty()) continue;
            if (line.find("#") < string::npos) continue; // Skip the comment.
            stringstream ss(line);
            int tindx; number tme;
            ss >> tindx >> tme;
            if (tindx == tindx_start) {
                SimState::t = tme;
                SimState::next_plot_indx = tindx + 1;
            }
        }
    }
    SimState::t_plot = SimState::t + Config::dt_plot_schedule.get_value(SimState::t);
    SimState::steps_since_last_plot = 0;
    return;
} // initialize_simulation()

__host__
void write_flow_data(int tindx, number tme)
{
    cout << "Write flow data at tindx=" << tindx << " time=" << tme << endl;
    //
    char nameBuf[256];
    sprintf(nameBuf, "%s/flow/t%04d", Config::job.c_str(), tindx);
    string flowDir = string(nameBuf);
    if (!filesystem::exists(flowDir)) { filesystem::create_directories(flowDir); }
    for (int blk_id=0; blk_id < Config::nFluidBlocks; ++blk_id) {
        BConfig& blk_config = Config::blk_configs[blk_id];
        int i = blk_config.i; int j = blk_config.j; int k = blk_config.k;
        sprintf(nameBuf, "%s/flow-%04d-%04d-%04d.zip", flowDir.c_str(), i, j, k);
        string fileName = string(nameBuf);
        fluidBlocks[blk_id].writeFlow(fileName);
    }
    // Update the times file.
    ofstream timesFile(Config::job+"/times.data", ofstream::binary|ofstream::app);
    timesFile << tindx << " " << tme << endl;
    return;
} // end write_flow_data()

// Repetitive boundary condition code is hidden here.
#include "bcs.cu"

__host__
void apply_boundary_conditions()
// Since the boundary-condition code needs a view of all blocks and
// most of the coperations are switching between code to copy specific data,
// we expect the CPU to apply the boundary conditions more effectively than the GPU.
// Measurements might tell us otherwise.
{
    for (int iblk=0; iblk < Config::nFluidBlocks; iblk++) {
        BConfig& blk_config = Config::blk_configs[iblk];
        if (blk_config.active) {
            for (int ibc=0; ibc < 6; ibc++) {
                switch (blk_config.bcCodes[ibc]) {
                case BCCode::wall_with_slip: bc_wall_with_slip(iblk, ibc); break;
                case BCCode::wall_no_slip: bc_wall_no_slip(iblk, ibc); break;
                case BCCode::exchange: bc_exchange(iblk, ibc); break;
                case BCCode::inflow: bc_inflow(iblk, ibc, Config::flow_states[blk_config.bc_fs[ibc]]); break;
                case BCCode::outflow: bc_outflow(iblk, ibc); break;
                default:
                    throw runtime_error("Invalid bcCode: "+to_string(blk_config.bcCodes[ibc]));
                }
            } // end for ibc
        } // end if active
    } // end for iblk
} // end apply_boundary_conditions()

__host__
void march_in_time()
{
    cout << "march_in_time() start" << endl;
    SimState::dt = Config::dt_init;
    SimState::step = 0;
    //
    while (SimState::step < Config::max_step && SimState::t < Config::max_time) {
        //
        // Occasionally determine allowable time step.
        if (SimState::step > 0 && (SimState::step % Config::cfl_count)==0) {
            number smallest_dt = numeric_limits<number>::max();
            number cfl = Config::cfl_schedule.get_value(SimState::t);
            for (Block& blk : fluidBlocks) {
                if (blk.active) { smallest_dt = fmin(smallest_dt, blk.estimate_allowed_dt(cfl)); }
            }
            SimState::dt = smallest_dt;
        }
        // Attempt a step, stage 1.
        apply_boundary_conditions();
        //
        int bad_cell_count = 0;
        for (Block& blk : fluidBlocks) {
            if (blk.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_count += blk.update_stage_1(SimState::dt);
            }
        }
        if (bad_cell_count == 0) {
            // After a successful step, copy the conserved data back to level 0.
            for (Block& blk : fluidBlocks) {
                if (blk.active) { blk.copy_conserved_data(1, 0); }
            }
        } else {
            throw runtime_error("Bad cell count: "+to_string(bad_cell_count));
        }
        //
        SimState::t += SimState::dt;
        SimState::step += 1;
        SimState::steps_since_last_plot += 1;
        //
        if (SimState::step > 0 && (SimState::step % Config::print_count)==0) {
            cout << "Step=" << SimState::step << " t=" << SimState::t
                 << " dt=" << SimState::dt << " cfl=" << Config::cfl_schedule.get_value(SimState::t)
                 << endl;
        }
        //
        // Occasionally dump the flow data for making plots.
        if (SimState::t >= SimState::t_plot) {
            write_flow_data(SimState::next_plot_indx, SimState::t);
            SimState::steps_since_last_plot = 0;
            SimState::next_plot_indx += 1;
            SimState::t_plot = SimState::t + Config::dt_plot_schedule.get_value(SimState::t);
        }
    } // end while loop
    cout << "march_in_time() end" << endl;
    return;
} // end march_in_time()

__host__
void finalize_simulation()
{
    if (SimState::steps_since_last_plot > 0) {
        write_flow_data(SimState::next_plot_indx, SimState::t);
    }
    return;
} // end finalize_simulation()

#endif

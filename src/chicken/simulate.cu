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
#include <chrono>

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
    auto clock_start = chrono::system_clock::now();
    read_config_file(Config::job + "/config.json");
    // Read initial grids and flow data
    for (int blk_id=0; blk_id < Config::nFluidBlocks; ++blk_id) {
        BConfig& cfg = Config::blk_configs[blk_id];
        int i = cfg.i; int j = cfg.j; int k = cfg.k;
        if (blk_id != Config::blk_ids[i][j][k]) {
            throw runtime_error("blk_id=" + to_string(blk_id) + " is inconsistent: i=" +
                                to_string(i) + " j=" + to_string(j) + " k=" + to_string(k));
        }
        cfg.fill_in_dimensions(Config::nics[i], Config::njcs[j], Config::nkcs[k]);
        Block blk; blk.configure(cfg);
        sprintf(nameBuf, "/grid/grid-%04d-%04d-%04d.gz", i, j, k);
        string fileName = Config::job + string(nameBuf);
        blk.readGrid(cfg, fileName);
        sprintf(nameBuf, "/flow/t%04d/flow-%04d-%04d-%04d.zip", tindx_start, i, j, k);
        fileName = Config::job + string(nameBuf);
        blk.readFlow(cfg, fileName);
        blk.computeGeometry(cfg);
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
    auto clock_now = chrono::system_clock::now();
    auto clock_ms = chrono::duration_cast<chrono::milliseconds>(clock_now - clock_start);
    cout << "initialize_simulation() finished in " << clock_ms.count() << "ms" << endl;
    return;
} // initialize_simulation()


__host__
void write_flow_data(int tindx, number tme)
{
    cout << "Write flow data at tindx=" << tindx
         << " time=" << scientific << setprecision(3) << tme << endl;
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
        fluidBlocks[blk_id].writeFlow(blk_config, fileName);
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
        if (!blk_config.active) continue;
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
    } // end for iblk
} // end apply_boundary_conditions()


__host__
void march_in_time_using_cpu_only()
// Variant of the main simulation function which uses only the CPU.
// We retain this function as a reasonably-easy-to-read reference code,
// while be build the GPU variant.
{
    if (Config::verbosity > 0) cout << "march_in_time_using_cpu_only() start" << endl;
    auto clock_start = chrono::system_clock::now();
    SimState::dt = Config::dt_init;
    SimState::step = 0;
    for (int ib=0; ib < Config::nFluidBlocks; ib++) {
        BConfig& cfg = Config::blk_configs[ib];
        if (cfg.active) fluidBlocks[ib].encodeConserved(cfg, 0);
    }
    //
    while (SimState::step < Config::max_step && SimState::t < Config::max_time) {
        //
        // Occasionally determine allowable time step.
        if (SimState::step > 0 && (SimState::step % Config::cfl_count)==0) {
            number smallest_dt = numeric_limits<number>::max();
            number cfl = Config::cfl_schedule.get_value(SimState::t);
            for (int ib=0; ib < Config::nFluidBlocks; ib++) {
                BConfig& cfg = Config::blk_configs[ib];
                Block& blk = fluidBlocks[ib];
                if (cfg.active) smallest_dt = fmin(smallest_dt, blk.estimate_allowed_dt(cfg, cfl));
            }
            SimState::dt = smallest_dt;
        }
        //
        // Gas-dynamic update over three stages with TVD-RK3 weights.
        int bad_cell_count = 0;
        // Stage 1.
        // number t = SimState::t; // Only needed if we have time-dependent source terms or BCs.
        apply_boundary_conditions();
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = Config::blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_count += blk.update_stage_1(cfg, SimState::dt);
            }
        }
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 1 bad cell count: "+to_string(bad_cell_count));
        }
        // Stage 2
        // t = SimState::t + 0.5*SimState::dt;
        apply_boundary_conditions();
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = Config::blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_count += blk.update_stage_2(cfg, SimState::dt);
            }
        }
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 2 bad cell count: "+to_string(bad_cell_count));
        }
        // Stage 3
        // t = SimState::t + SimState::dt;
        apply_boundary_conditions();
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = Config::blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_count += blk.update_stage_3(cfg, SimState::dt);
            }
        }
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 3 bad cell count: "+to_string(bad_cell_count));
        }
        // After a successful gasdynamic update, copy the conserved data back to level 0.
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = Config::blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) blk.copy_conserved_data(cfg, 1, 0);
        }
        //
        SimState::t += SimState::dt;
        SimState::step += 1;
        SimState::steps_since_last_plot += 1;
        //
        // Occasionally write the current step and time to the console.
        if (SimState::step > 0 && (SimState::step % Config::print_count)==0) {
            auto clock_now = chrono::system_clock::now();
            auto clock_ms = chrono::duration_cast<chrono::milliseconds>(clock_now - clock_start);
            double wall_clock_elapsed = clock_ms.count()/1000.0;
            double wall_clock_per_step = wall_clock_elapsed / SimState::step;
            double WCtFT = (Config::max_time - SimState::t) / SimState::dt * wall_clock_per_step;
            double WCtMS = (Config::max_step - SimState::step) * wall_clock_per_step;
            cout << "Step=" << SimState::step
                 << " t=" << scientific << setprecision(3) << SimState::t
                 << " dt=" << SimState::dt
                 << " cfl=" << fixed << Config::cfl_schedule.get_value(SimState::t)
                 << " WC=" << wall_clock_elapsed << "s"
                 << " WCtFT=" << WCtFT << "s"
                 << " WCtMS=" << WCtMS << "s"
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
    if (Config::verbosity > 0) cout << "march_in_time_using_cpu_only() end" << endl;
    return;
} // end march_in_time_using_cpu_only()


__host__
void march_in_time_using_gpu()
// Variant of the main simulation function where we may offload work to the GPU.
{
    if (Config::verbosity > 0) cout << "march_in_time_using_gpu() start" << endl;
    auto clock_start = chrono::system_clock::now();
    SimState::dt = Config::dt_init;
    SimState::step = 0;
    for (int ib=0; ib < Config::nFluidBlocks; ib++) {
        BConfig& cfg = Config::blk_configs[ib];
        if (cfg.active) fluidBlocks[ib].encodeConserved(cfg, 0);
    }
    //
    while (SimState::step < Config::max_step && SimState::t < Config::max_time) {
        //
        // Occasionally determine allowable time step.
        if (SimState::step > 0 && (SimState::step % Config::cfl_count)==0) {
            number smallest_dt = numeric_limits<number>::max();
            number cfl = Config::cfl_schedule.get_value(SimState::t);
            for (int ib=0; ib < Config::nFluidBlocks; ib++) {
                BConfig& cfg = Config::blk_configs[ib];
                Block& blk = fluidBlocks[ib];
                if (cfg.active) smallest_dt = fmin(smallest_dt, blk.estimate_allowed_dt(cfg, cfl));
            }
            SimState::dt = smallest_dt;
        }
        //
        // Gas-dynamic update over three stages with TVD-RK3 weights.
        int bad_cell_count = 0;
        // Stage 1.
        // number t = SimState::t; // Only needed if we have time-dependent source terms or BCs.
        apply_boundary_conditions();
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = Config::blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_count += blk.update_stage_1(cfg, SimState::dt);
            }
        }
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 1 bad cell count: "+to_string(bad_cell_count));
        }
        // Stage 2
        // t = SimState::t + 0.5*SimState::dt;
        apply_boundary_conditions();
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = Config::blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_count += blk.update_stage_2(cfg, SimState::dt);
            }
        }
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 2 bad cell count: "+to_string(bad_cell_count));
        }
        // Stage 3
        // t = SimState::t + SimState::dt;
        apply_boundary_conditions();
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = Config::blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_count += blk.update_stage_3(cfg, SimState::dt);
            }
        }
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 3 bad cell count: "+to_string(bad_cell_count));
        }
        // After a successful gasdynamic update, copy the conserved data back to level 0.
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = Config::blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) blk.copy_conserved_data(cfg, 1, 0);
        }
        //
        SimState::t += SimState::dt;
        SimState::step += 1;
        SimState::steps_since_last_plot += 1;
        //
        // Occasionally write the current step and time to the console.
        if (SimState::step > 0 && (SimState::step % Config::print_count)==0) {
            auto clock_now = chrono::system_clock::now();
            auto clock_ms = chrono::duration_cast<chrono::milliseconds>(clock_now - clock_start);
            double wall_clock_elapsed = clock_ms.count()/1000.0;
            double wall_clock_per_step = wall_clock_elapsed / SimState::step;
            double WCtFT = (Config::max_time - SimState::t) / SimState::dt * wall_clock_per_step;
            double WCtMS = (Config::max_step - SimState::step) * wall_clock_per_step;
            cout << "Step=" << SimState::step
                 << " t=" << scientific << setprecision(3) << SimState::t
                 << " dt=" << SimState::dt
                 << " cfl=" << fixed << Config::cfl_schedule.get_value(SimState::t)
                 << " WC=" << wall_clock_elapsed << "s"
                 << " WCtFT=" << WCtFT << "s"
                 << " WCtMS=" << WCtMS << "s"
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
    if (Config::verbosity > 0) cout << "march_in_time_using_gpu() end" << endl;
    return;
} // end march_in_time_using_gpu()


__host__
void finalize_simulation()
{
    if (SimState::steps_since_last_plot > 0) {
        write_flow_data(SimState::next_plot_indx, SimState::t);
    }
    for (Block& blk : fluidBlocks) {
        blk.releaseMemory();
    }
    return;
} // end finalize_simulation()

#endif

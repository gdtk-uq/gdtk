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
#include <omp.h>

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
Block* fluidBlocks_on_gpu;
vector<BConfig> blk_configs;
BConfig* blk_configs_on_gpu;

__host__
void initialize_simulation(int tindx_start)
{
    char nameBuf[256];
    filesystem::path pth{Config::job};
    if (!filesystem::exists(pth) || !filesystem::is_directory(pth)) {
        throw runtime_error("Job directory is not present in current directory.");
    }
    auto clock_start = chrono::system_clock::now();
    blk_configs = read_config_file(Config::job + "/config.json");
    // Read initial grids and flow data
    size_t bytes_allocated = 0;
    size_t cells_in_simulation = 0;
    for (int blk_id=0; blk_id < Config::nFluidBlocks; ++blk_id) {
        BConfig& cfg = blk_configs[blk_id];
        bytes_allocated += sizeof(BConfig);
        int i = cfg.i; int j = cfg.j; int k = cfg.k;
        if (blk_id != Config::blk_ids[i][j][k]) {
            throw runtime_error("blk_id=" + to_string(blk_id) + " is inconsistent: i=" +
                                to_string(i) + " j=" + to_string(j) + " k=" + to_string(k));
        }
        cfg.fill_in_dimensions(Config::nics[i], Config::njcs[j], Config::nkcs[k]);
        Block blk;
        bytes_allocated += sizeof(Block) + blk.configure(cfg);
        sprintf(nameBuf, "/grid/grid-%04d-%04d-%04d.gz", i, j, k);
        string fileName = Config::job + string(nameBuf);
        blk.readGrid(cfg, fileName);
        sprintf(nameBuf, "/flow/t%04d/flow-%04d-%04d-%04d.zip", tindx_start, i, j, k);
        fileName = Config::job + string(nameBuf);
        blk.readFlow(cfg, fileName);
        blk.computeGeometry(cfg);
        fluidBlocks.push_back(blk);
        if (cfg.active) cells_in_simulation += cfg.nic*cfg.njc*cfg.nkc;
    }
    cout << "Cells in simulation: " << cells_in_simulation << endl;
    cout << "Bytes allocated on CPU: " << fixed << setprecision(3) << bytes_allocated/1.0e6 << "MB" << endl;
#ifdef CUDA
    // We need to put a copy of the block and config data onto the GPU.
    int nbytes = blk_configs.size()*sizeof(BConfig);
    auto status = cudaMalloc(&blk_configs_on_gpu, nbytes);
    if (status) {
        cerr << cudaGetErrorString(cudaGetLastError()) << endl;
        throw runtime_error("Could not allocate blk_configs on gpu.");
    }
    status = cudaMemcpy(blk_configs_on_gpu, blk_configs.data(), nbytes, cudaMemcpyHostToDevice);
    if (status) {
        cerr << cudaGetErrorString(cudaGetLastError()) << endl;
        throw runtime_error("Could not copy blk_configs to gpu.");
    }
    //
    nbytes = fluidBlocks.size()*sizeof(Block);
    status = cudaMalloc(&fluidBlocks_on_gpu, nbytes);
    if (status) {
        cerr << cudaGetErrorString(cudaGetLastError()) << endl;;
        throw runtime_error("Could not allocate fluidBlocks on gpu.");
    }
    status = cudaMemcpy(fluidBlocks_on_gpu, fluidBlocks.data(), nbytes, cudaMemcpyHostToDevice);
    if (status) {
        cerr << cudaGetErrorString(cudaGetLastError()) << endl;;
        throw runtime_error("Could not copy fluidBlocks to gpu.");
    }
#endif
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
        BConfig& blk_config = blk_configs[blk_id];
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
    #pragma omp parallel for
    for (int iblk=0; iblk < Config::nFluidBlocks; iblk++) {
        BConfig& blk_config = blk_configs[iblk];
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
    // A couple of global arrays to regulate the simulation.
    vector<number> allowed_dts; allowed_dts.resize(Config::nFluidBlocks);
    vector<int> bad_cell_counts; bad_cell_counts.resize(Config::nFluidBlocks);
    //
    #pragma omp parallel for
    for (int ib=0; ib < Config::nFluidBlocks; ib++) {
        BConfig& cfg = blk_configs[ib];
        if (cfg.active) fluidBlocks[ib].encodeConserved(cfg, 0);
    }
    //
    while (SimState::step < Config::max_step && SimState::t < Config::max_time) {
        //
        // Occasionally determine allowable time step.
        if (SimState::step > 0 && (SimState::step % Config::cfl_count)==0) {
            number smallest_dt = numeric_limits<number>::max();
            for (auto& adt : allowed_dts) adt = smallest_dt;
            number cfl = Config::cfl_schedule.get_value(SimState::t);
            #pragma omp parallel for
            for (int ib=0; ib < Config::nFluidBlocks; ib++) {
                BConfig& cfg = blk_configs[ib];
                Block& blk = fluidBlocks[ib];
                if (cfg.active) allowed_dts[ib] = blk.estimate_allowed_dt(cfg, cfl);
            }
            for (auto adt : allowed_dts) smallest_dt = fmin(smallest_dt, adt);
            SimState::dt = smallest_dt;
        }
        //
        // Gas-dynamic update over three stages with TVD-RK3 weights.
        int bad_cell_count = 0;
        // Stage 1.
        // number t = SimState::t; // Only needed if we have time-dependent source terms or BCs.
        apply_boundary_conditions();
        for (auto& bcc : bad_cell_counts) bcc = 0;
        #pragma omp parallel for
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_counts[ib] = blk.update_stage_1(cfg, SimState::dt);
            }
        }
        for (auto bcc : bad_cell_counts) bad_cell_count += bcc;
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 1 bad cell count: "+to_string(bad_cell_count));
        }
        // Stage 2
        // t = SimState::t + 0.5*SimState::dt;
        apply_boundary_conditions();
        for (auto& bcc : bad_cell_counts) bcc = 0;
        #pragma omp parallel for
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_counts[ib] = blk.update_stage_2(cfg, SimState::dt);
            }
        }
        for (auto bcc : bad_cell_counts) bad_cell_count += bcc;
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 2 bad cell count: "+to_string(bad_cell_count));
        }
        // Stage 3
        // t = SimState::t + SimState::dt;
        apply_boundary_conditions();
        for (auto& bcc : bad_cell_counts) bcc = 0;
        #pragma omp parallel for
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (cfg.active) {
                blk.calculate_fluxes(Config::x_order);
                bad_cell_counts[ib] = blk.update_stage_3(cfg, SimState::dt);
            }
        }
        for (auto bcc : bad_cell_counts) bad_cell_count += bcc;
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 3 bad cell count: "+to_string(bad_cell_count));
        }
        // After a successful gasdynamic update, copy the conserved data back to level 0.
        #pragma omp parallel for
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
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
        BConfig& cfg = blk_configs[ib];
        if (!cfg.active) continue;
        auto& blk = fluidBlocks[ib];
        // Transfer block data, including the initial flow states, to the GPU and encode.
        int nbytes = blk.cells.size()*sizeof(FVCell);
        auto status = cudaMemcpy(blk.cells_on_gpu, blk.cells.data(), nbytes, cudaMemcpyHostToDevice);
        if (status) {
            cerr << cudaGetErrorString(cudaGetLastError()) << endl;
            throw runtime_error("Could not copy blk.cells to gpu.");
        }
        // No need to send conserved quantities and their time-derivatives
        // but we do want to send faces and vertices.
        nbytes = blk.iFaces.size()*sizeof(FVFace);
        status = cudaMemcpy(blk.iFaces_on_gpu, blk.iFaces.data(), nbytes, cudaMemcpyHostToDevice);
        if (status) {
            cerr << cudaGetErrorString(cudaGetLastError()) << endl;
            throw runtime_error("Could not copy blk.iFaces to gpu.");
        }
        nbytes = blk.jFaces.size()*sizeof(FVFace);
        status = cudaMemcpy(blk.jFaces_on_gpu, blk.jFaces.data(), nbytes, cudaMemcpyHostToDevice);
        if (status) {
            cerr << cudaGetErrorString(cudaGetLastError()) << endl;
            throw runtime_error("Could not copy blk.jFaces to gpu.");
        }
        nbytes = blk.kFaces.size()*sizeof(FVFace);
        status = cudaMemcpy(blk.kFaces_on_gpu, blk.kFaces.data(), nbytes, cudaMemcpyHostToDevice);
        if (status) {
            cerr << cudaGetErrorString(cudaGetLastError()) << endl;
            throw runtime_error("Could not copy blk.kFaces to gpu.");
        }
        nbytes = blk.vertices.size()*sizeof(Vector3);
        status = cudaMemcpy(blk.vertices_on_gpu, blk.vertices.data(), nbytes, cudaMemcpyHostToDevice);
        if (status) {
            cerr << cudaGetErrorString(cudaGetLastError()) << endl;
            throw runtime_error("Could not copy blk.vertices to gpu.");
        }
        // Now, do the encode of flow states to conserved quantities.
        Block& blk_on_gpu = fluidBlocks_on_gpu[ib];
        BConfig& cfg_on_gpu = blk_configs_on_gpu[ib];
        int nGPUblocks = cfg.nGPUblocks_for_cells;
        int nGPUthreads = cfg.threads_per_GPUblock;
        encodeConserved_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, 0);
        auto cudaError = cudaGetLastError();
        if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
    }
    //
    // A couple of global variables for keeping an eye on the simulation process.
    int bad_cell_count = 0;
    int* bad_cell_count_on_gpu;
    auto status = cudaMalloc(&bad_cell_count_on_gpu, sizeof(int));
    if (status) throw runtime_error("Could not allocate bad_cell_count_on_gpu.");
    //
    long long int* smallest_dt_picos_on_gpu;
    status = cudaMalloc(&smallest_dt_picos_on_gpu, sizeof(long long int));
    if (status) throw runtime_error("Could not allocate smallest_dt_picos_on_gpu.");
    //
    // Main stepping loop.
    //
    while (SimState::step < Config::max_step && SimState::t < Config::max_time) {
        //
        // Occasionally determine allowable time step.
        if (SimState::step > 0 && (SimState::step % Config::cfl_count)==0) {
            long long int smallest_dt_picos = numeric_limits<long long int>::max();
            status = cudaMemcpy(smallest_dt_picos_on_gpu, &smallest_dt_picos,
                sizeof(long long int), cudaMemcpyHostToDevice);
            if (status) throw runtime_error("Could not copy smallest_dt_picos to gpu.");
            number cfl = Config::cfl_schedule.get_value(SimState::t);
            for (int ib=0; ib < Config::nFluidBlocks; ib++) {
                BConfig& cfg = blk_configs[ib];
                if (!cfg.active) continue;
                Block& blk = fluidBlocks[ib];
                Block& blk_on_gpu = fluidBlocks_on_gpu[ib];
                BConfig& cfg_on_gpu = blk_configs_on_gpu[ib];
                int nGPUblocks = cfg.nGPUblocks_for_cells;
                int nGPUthreads = cfg.threads_per_GPUblock;
                estimate_allowed_dt_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu,
                                                                       cfl, smallest_dt_picos_on_gpu);
                auto cudaError = cudaGetLastError();
                if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            }
            status = cudaMemcpy(&smallest_dt_picos, smallest_dt_picos_on_gpu,
                sizeof(long long int), cudaMemcpyDeviceToHost);
            if (status) throw runtime_error("Could not copy smallest_dt_picos from gpu to host cpu.");
            SimState::dt = smallest_dt_picos * 1.0e-12;
        }
        //
        // Gas-dynamic update over three stages with TVD-RK3 weights.
        bad_cell_count = 0;
        status = cudaMemcpy(bad_cell_count_on_gpu, &bad_cell_count, sizeof(int), cudaMemcpyHostToDevice);
        if (status) throw runtime_error("Could not copy bad_cell_count to gpu.");
        //
        // Stage 1.
        // number t = SimState::t; // Only needed if we have time-dependent source terms or BCs.
        apply_boundary_conditions();
        // [TODO] PJ 2022-10-15
        // Boundary-conditions are done on the host CPU, affecting only the ghost-cell data.
        // If we continue to do this, we should copy just the ghost cell data onto the GPU,
        // however, we could do the boundary conditions on the GPU is all of the blocks fit.
        // That would eliminate lots of copying back and forth.
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            int nbytes = blk.cells.size()*sizeof(FVCell);
            auto status = cudaMemcpy(blk.cells_on_gpu, blk.cells.data(), nbytes, cudaMemcpyHostToDevice);
            if (status) {
                cerr << cudaGetErrorString(cudaGetLastError()) << endl;
                throw runtime_error("Could not copy blk.cells to gpu.");
            }
        }
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            // Do the stage-1 update on the GPU.
            Block& blk_on_gpu = fluidBlocks_on_gpu[ib];
            BConfig& cfg_on_gpu = blk_configs_on_gpu[ib];
            int nGPUblocks = cfg.nGPUblocks_for_iFaces;
            int nGPUthreads = cfg.threads_per_GPUblock;
            calculate_iFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            auto cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_jFaces;
            nGPUthreads = cfg.threads_per_GPUblock;
            calculate_jFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_kFaces;
            nGPUthreads = cfg.threads_per_GPUblock;
            calculate_kFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_cells;
            nGPUthreads = cfg.threads_per_GPUblock;
            update_stage_1_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu,
                                                              SimState::dt, bad_cell_count_on_gpu);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
        }
        status = cudaMemcpy(&bad_cell_count, bad_cell_count_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
        if (status) throw runtime_error("Could not copy bad_cell_count from gpu to host cpu.");
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 1 bad cell count: "+to_string(bad_cell_count));
        }
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            int nbytes = blk.cells.size()*sizeof(FVCell);
            auto status = cudaMemcpy(blk.cells.data(), blk.cells_on_gpu, nbytes, cudaMemcpyDeviceToHost);
            if (status) {
                cerr << cudaGetErrorString(cudaGetLastError()) << endl;
                throw runtime_error("Could not copy blk.cells from gpu to cpu.");
            }
        }
        //
        // Stage 2
        // t = SimState::t + 0.5*SimState::dt;
        apply_boundary_conditions();
        //
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            int nbytes = blk.cells.size()*sizeof(FVCell);
            auto status = cudaMemcpy(blk.cells_on_gpu, blk.cells.data(), nbytes, cudaMemcpyHostToDevice);
            if (status) {
                cerr << cudaGetErrorString(cudaGetLastError()) << endl;
                throw runtime_error("Could not copy blk.cells to gpu.");
            }
        }
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            // Do the stage-1 update on the GPU.
            Block& blk_on_gpu = fluidBlocks_on_gpu[ib];
            BConfig& cfg_on_gpu = blk_configs_on_gpu[ib];
            int nGPUblocks = cfg.nGPUblocks_for_iFaces;
            int nGPUthreads = cfg.threads_per_GPUblock;
            calculate_iFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            auto cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_jFaces;
            nGPUthreads = cfg.threads_per_GPUblock;
            calculate_jFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_kFaces;
            nGPUthreads = cfg.threads_per_GPUblock;
            calculate_kFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_cells;
            nGPUthreads = cfg.threads_per_GPUblock;
            update_stage_2_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu,
                                                              SimState::dt, bad_cell_count_on_gpu);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
        }
        status = cudaMemcpy(&bad_cell_count, bad_cell_count_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
        if (status) throw runtime_error("Could not copy bad_cell_count from gpu to host cpu.");
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 2 bad cell count: "+to_string(bad_cell_count));
        }
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            int nbytes = blk.cells.size()*sizeof(FVCell);
            auto status = cudaMemcpy(blk.cells.data(), blk.cells_on_gpu, nbytes, cudaMemcpyDeviceToHost);
            if (status) {
                cerr << cudaGetErrorString(cudaGetLastError()) << endl;
                throw runtime_error("Could not copy blk.cells from gpu to cpu.");
            }
        }
        //
        // Stage 3
        // t = SimState::t + SimState::dt;
        apply_boundary_conditions();
        //
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            int nbytes = blk.cells.size()*sizeof(FVCell);
            auto status = cudaMemcpy(blk.cells_on_gpu, blk.cells.data(), nbytes, cudaMemcpyHostToDevice);
            if (status) {
                cerr << cudaGetErrorString(cudaGetLastError()) << endl;
                throw runtime_error("Could not copy blk.cells to gpu.");
            }
        }
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            // Do the stage-1 update on the GPU.
            Block& blk_on_gpu = fluidBlocks_on_gpu[ib];
            BConfig& cfg_on_gpu = blk_configs_on_gpu[ib];
            int nGPUblocks = cfg.nGPUblocks_for_iFaces;
            int nGPUthreads = cfg.threads_per_GPUblock;
            calculate_iFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            auto cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_jFaces;
            nGPUthreads = cfg.threads_per_GPUblock;
            calculate_jFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_kFaces;
            nGPUthreads = cfg.threads_per_GPUblock;
            calculate_kFace_fluxes_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, Config::x_order);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
            //
            nGPUblocks = cfg.nGPUblocks_for_cells;
            nGPUthreads = cfg.threads_per_GPUblock;
            update_stage_3_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu,
                                                              SimState::dt, bad_cell_count_on_gpu);
            cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
        }
        status = cudaMemcpy(&bad_cell_count, bad_cell_count_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
        if (status) throw runtime_error("Could not copy bad_cell_count from gpu to host cpu.");
        if (bad_cell_count > 0) {
            throw runtime_error("Stage 3 bad cell count: "+to_string(bad_cell_count));
        }
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            int nbytes = blk.cells.size()*sizeof(FVCell);
            auto status = cudaMemcpy(blk.cells.data(), blk.cells_on_gpu, nbytes, cudaMemcpyDeviceToHost);
            if (status) {
                cerr << cudaGetErrorString(cudaGetLastError()) << endl;
                throw runtime_error("Could not copy blk.cells from gpu to cpu.");
            }
        }
        // After a successful gasdynamic update, copy the conserved data back to level 0.
        for (int ib=0; ib < Config::nFluidBlocks; ib++) {
            BConfig& cfg = blk_configs[ib];
            Block& blk = fluidBlocks[ib];
            if (!cfg.active) continue;
            Block& blk_on_gpu = fluidBlocks_on_gpu[ib];
            BConfig& cfg_on_gpu = blk_configs_on_gpu[ib];
            int nGPUblocks = cfg.nGPUblocks_for_cells;
            int nGPUthreads = cfg.threads_per_GPUblock;
            copy_conserved_data_on_gpu<<<nGPUblocks,nGPUthreads>>>(blk_on_gpu, cfg_on_gpu, 1, 0);
            auto cudaError = cudaGetLastError();
            if (cudaError) throw runtime_error(cudaGetErrorString(cudaError));
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

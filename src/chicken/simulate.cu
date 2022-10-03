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
#include <string>
#include <filesystem>
#include <limits>

#include "number.cu"
#include "vector3.cu"
#include "gas.cu"
#include "flow.cu"
#include "vertex.cu"
#include "face.cu"
#include "flux.cu"
#include "cell.cu"
#include "block.cu"
#include "config.cu"

using namespace std;

struct SimState {
    number dt;
    int step;
    int max_step;
    number t_final;
    number dt_plot;
};

vector<Block*> fluidBlocks;

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
    for (int k=0; k < Config::nkb; ++k) {
        for (int j=0; j < Config::njb; ++j) {
            for (int i=0; i < Config::nib; ++i) {
                if (Config::blk_ids[i][j][k] >= 0) {
                    // Only defined blocks in the array will have a non-zero id.
                    Block* blk_ptr = new Block{};
                    int blk_id = Config::blk_ids[i][j][k];
                    blk_ptr->configure(Config::nics[i], Config::njcs[j], Config::nkcs[k]);
                    sprintf(nameBuf, "/grid/grid-%04d-%04d-%04d.gz", i, j, k);
                    string fileName = Config::job + string(nameBuf);
                    blk_ptr->readGrid(fileName);
                    sprintf(nameBuf, "/flow/t%04d/flow-%04d-%04d-%04d.zip", tindx_start, i, j, k);
                    fileName = Config::job + string(nameBuf);
                    blk_ptr->readFlow(fileName);
                    blk_ptr->computeGeometry();
                    blk_ptr->encodeConserved(0);
                    cout << "Sample cell data: " << blk_ptr->cells[blk_ptr->activeCellIndex(0,0,0)].toString() << endl;
                    cout << "Sample iFace data: " << blk_ptr->iFaces[blk_ptr->iFaceIndex(0,0,0)].toString() << endl;
                    cout << "Sample jFace data: " << blk_ptr->jFaces[blk_ptr->jFaceIndex(0,0,0)].toString() << endl;
                    cout << "Sample kFace data: " << blk_ptr->kFaces[blk_ptr->kFaceIndex(0,0,0)].toString() << endl;
                    fluidBlocks.push_back(blk_ptr);
                    if (blk_id+1 != fluidBlocks.size()) {
                        throw runtime_error("Inconsistent blk_id and position in fluidBlocks array.");
                    }
                }
            }
        }
    }
    if (fluidBlocks.size() != Config::nFluidBlocks) {
        throw runtime_error("Inconsistent number of blocks: "+
                            to_string(fluidBlocks.size())+" "+to_string(Config::nFluidBlocks));
    }
    return;
} // initialize_simulation()

__host__
void write_flow_data(int tindx)
{
    char nameBuf[256];
    sprintf(nameBuf, "%s/flow/t%04d", Config::job.c_str(), tindx);
    string flowDir = string(nameBuf);
    if (!filesystem::exists(flowDir)) { filesystem::create_directories(flowDir); }
    for (int k=0; k < Config::nkb; ++k) {
        for (int j=0; j < Config::njb; ++j) {
            for (int i=0; i < Config::nib; ++i) {
                if (Config::blk_ids[i][j][k] >= 0) {
                    // Only defined blocks in the array will have a non-zero id.
                    int blk_id = Config::blk_ids[i][j][k];
                    Block* blk_ptr = fluidBlocks[blk_id];
                    sprintf(nameBuf, "%s/flow-%04d-%04d-%04d.zip", flowDir.c_str(), i, j, k);
                    string fileName = string(nameBuf);
                    blk_ptr->writeFlow(fileName);
                }
            }
        }
    }
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
        auto* blk_config = &(Config::blk_configs[iblk]);
        auto* blk_ptr = fluidBlocks[iblk];
        for (int ibc=0; ibc < 6; ibc++) {
            switch (blk_config->bcCodes[ibc]) {
            case BCCode::wall_with_slip: bc_wall_with_slip(blk_ptr, ibc); break;
            case BCCode::wall_no_slip: bc_wall_no_slip(blk_ptr, ibc); break;
            case BCCode::exchange: bc_exchange(iblk, ibc); break;
            case BCCode::inflow: bc_inflow(blk_ptr, ibc, Config::flow_states[blk_config->bc_fs[ibc]]); break;
            case BCCode::outflow: bc_outflow(blk_ptr, ibc); break;
            default:
                throw runtime_error("Invalid bcCode: "+to_string(blk_config->bcCodes[ibc]));
            }
        } // end for ibc
    } // end for iblk
} // end apply_boundary_conditions()

__host__
void gasdynamic_update(number dt)
{
    apply_boundary_conditions();
    // update_stage_1 for all blocks
    //
} // end gasdynamic_update()

__host__
void march_in_time()
{
    // Occasionally determine allowable time step.
    number dt = 1.0e-6; // [TODO] get from config
    for (auto* blk_ptr : fluidBlocks) {
        dt = fmin(dt, blk_ptr->estimate_allowed_dt(0.5));
    }
    gasdynamic_update(dt);
    // Call gasdynamic_update an number of times.
    for (auto* blk_ptr : fluidBlocks) {
        int bad_cell = blk_ptr->decodeConserved(0);
    }
    return;
}

__host__
void finalize_simulation()
{
    // Exercise the writing of flow data, even we have done no calculations.
    write_flow_data(1);
    return;
}

#endif

// simulate.cu
// Include file for chicken, high-level simulation functions.
//
// PJ 2022-09-09

#ifndef SIMULATE_INCLUDED
#define SIMULATE_INCLUDED

//#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>
#include <filesystem>

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

void do_something(); // left over from the CUDA workshop experiment

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
                    Block* blkptr = new Block{};
                    int blk_id = Config::blk_ids[i][j][k];
                    blkptr->configure(Config::nics[i], Config::njcs[j], Config::nkcs[k],
                                      Config::blk_configs[blk_id].bcCodes);
                    sprintf(nameBuf, "/grid/grid-%04d-%04d-%04d.gz", i, j, k);
                    string fileName = Config::job + string(nameBuf);
                    blkptr->readGrid(fileName);
                    sprintf(nameBuf, "/flow/t%04d/flow-%04d-%04d-%04d.zip", tindx_start, i, j, k);
                    fileName = Config::job + string(nameBuf);
                    blkptr->readFlow(fileName);
                    blkptr->computeGeometry();
                    blkptr->encodeConserved(0);
                    cout << "Sample cell data: " << blkptr->cells[blkptr->activeCellIndex(0,0,0)].toString() << endl;
                    cout << "Sample iFace data: " << blkptr->iFaces[blkptr->iFaceIndex(0,0,0)].toString() << endl;
                    cout << "Sample jFace data: " << blkptr->jFaces[blkptr->jFaceIndex(0,0,0)].toString() << endl;
                    cout << "Sample kFace data: " << blkptr->kFaces[blkptr->kFaceIndex(0,0,0)].toString() << endl;
                    fluidBlocks.push_back(blkptr);
                    if (blk_id+1 != fluidBlocks.size()) {
                        throw runtime_error("Inconsistent blk_id and fluidBlocks array.");
                    }
                }
            }
        }
    }
    do_something();
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
                    Block* blkptr = fluidBlocks[blk_id];
                    sprintf(nameBuf, "%s/flow-%04d-%04d-%04d.zip", flowDir.c_str(), i, j, k);
                    string fileName = string(nameBuf);
                    blkptr->writeFlow(fileName);
                }
            }
        }
    }
    return;
} // end write_flow_data()

__host__
void march_in_time()
{
    for (auto* blkptr : fluidBlocks) {
        int bad_cell = blkptr->decodeConserved(0);
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

//---------------------------------------------------------------------------
//
// Bits left over from the CUDA workshop experiment.
// Initial hack adapts the vector addition example from the CUDA workshop
// to look a bit closer to our Puffin CFD code.
//
void host_process(vector<FlowState>& fss)
{
    for (auto& fs : fss) {
        auto& gas = fs.gas;
        auto& vel = fs.vel;
        number v2 = vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;
        number v = sqrt(v2);
        number M = v/gas.a;
        number g = GasModel::g;
        number t1 = 1.0f + 0.5f*(g-1.0)*M*M;
        // Compute stagnation condition.
        number p_total = gas.p * pow(t1, (g/(g-1.0)));
        number T_total = gas.T * t1;
        gas.p = p_total;
        gas.T = T_total;
        gas.update_from_pT();
        vel = {0.0, 0.0, 0.0};
    }
    cout << "inside host_process: fss[0]= " << fss[0].toString() << endl;
}

__global__ void device_process(FlowState* fss, int N)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    //
    if (idx < N) {
        auto& fs = fss[idx];
        auto gas = fs.gas;
        auto vel = fs.vel;
        number v2 = vel.x*vel.x + vel.y*vel.y + vel.z*vel.z;
        number v = sqrt(v2);
        number M = v/gas.a;
        number g = GasModel::g;
        number t1 = 1.0f + 0.5f*(g-1.0f)*M*M;
        // Compute stagnation condition.
        number p_total = gas.p * pow(t1, (g/(g-1.0f)));
        number T_total = gas.T * t1;
        gas.p = p_total;
        gas.T = T_total;
        gas.update_from_pT();
        vel = {0.0f, 0.0f, 0.0f};
        fs.gas = gas;
        fs.vel = vel;
    }
}

void print_sample(vector<FlowState> fss)
{
    for (int idx=0; idx < 3; idx++) {
        auto& fs = fss[idx];
        cout << "fs= " << fs.toString() << endl;
    }
    cout << "..." << endl;
    int N = fss.size();
    for (int idx=N-3; idx < N; idx++) {
        auto& fs = fss[idx];
        cout << "fs=" << fs.toString() << endl;
    }
}

void do_something()
{
    // Host data is in a standard C++ vector.
    vector<FlowState> fss_h;
    const int N = 32*512;
    for (int idx=0; idx < N; idx++) {
        auto gas = GasState{0.0, 0.0, 100.0e3, 300.0, 0.0};
        gas.update_from_pT();
        auto vel = Vector3{1000.0, 99.0, 0.0};
        fss_h.push_back(FlowState{gas, vel});
    }
    #ifdef CUDA
    if (!filesystem::exists(filesystem::status("/proc/driver/nvidia"))) {
        throw runtime_error("Cannot find NVIDIA driver in /proc/driver.");
    }
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    cout << "Found " << nDevices << " CUDA devices." << endl;
    if (nDevices > 0) {
        cout << "We have a CUDA device, so use it." << endl;
        // Pointer to device arrays.
        FlowState* fss_d;
        int sze = N * sizeof(FlowState);
        cudaMalloc(&fss_d, sze);
        cudaMemcpy(fss_d, fss_h.data(), sze, cudaMemcpyHostToDevice);
        //
        const int threads_per_block = 128;
        const int nblocks = N/threads_per_block;
        device_process<<<nblocks,threads_per_block>>>(fss_d, N);
        cout << cudaGetErrorString(cudaGetLastError()) << endl;
        //
        cudaMemcpy(fss_h.data(), fss_d, sze, cudaMemcpyDeviceToHost);
        cudaFree(fss_d);
    } else {
        cout << "Fall back to CPU-only processing." << endl;
        host_process(fss_h);
    }
    #else
    host_process(fss_h);
    #endif
    print_sample(fss_h);
    fss_h.resize(0);
    return;
}

#endif

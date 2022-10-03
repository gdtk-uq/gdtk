// left-over-code.cu
// CUDA bootcamp left-overs.
// Put here so that code might be cribbed.


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

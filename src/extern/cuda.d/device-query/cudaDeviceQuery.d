// author: KD
// CUDA device query

import std.stdio;
import std.string;
import cuda_d.cuda;
import cuda_d.cuda_runtime_api;
import cuda_d.cublas_api;

void main() {
    cuInit(0); // this has to be 0
    CUdevprop prop;
    int count;

    CUresult cuDeviceGetCountReturnState = cuDeviceGetCount( &count );
    if (cuDeviceGetCountReturnState == 0) writeln("Number of CUDA Devices: ", count);
    else writeln("cuDeviceGetCountReturnState = ", cuDeviceGetCountReturnState);
   
    for (int dev=0; dev < count; dev++) {
   
        cuDeviceGetProperties( &prop, dev );
	CUresult cuDeviceGetPropertiesReturnState = cuDeviceGetProperties( &prop, dev );
        if (cuDeviceGetPropertiesReturnState == 0) writeln("Successfully collected properties for Device ", dev); 
      	
	writef( "--- General Information for device %d ---\n", dev );
	char[64] name;
	int len = 64;		
	cuDeviceGetName (&name, len, dev);
	writeln( "Name: ", name );

	int major; int minor;
	cuDeviceComputeCapability (&major, &minor, dev);
	writef( "Compute capability: %d.%d \n", major, minor);
	
	writef( "Clock rate: %d kilohertz\n", prop.clockRate );
        
        writef( "--- Memory Information for device %d ---\n", dev );
	size_t bytes;
        cuDeviceTotalMem_v2 (&bytes, dev);
        writef( "Total global mem: %d\n", bytes );
        
        writef( "Total constant Mem: %d\n", prop.totalConstantMemory );
        writef( "Max mem pitch: %d\n", prop.memPitch );
        writef( "Texture Alignment: %d\n", prop.textureAlign );
        writef( "--- MP Information for device %d ---\n", dev );

       	writef( "Shared mem per mp: %d\n", prop.sharedMemPerBlock );
        writef( "Registers per mp: %d\n", prop.regsPerBlock );
        
        writef( "Threads in warp: %d\n", prop.SIMDWidth );
        
        writef( "Max threads per block: %d\n", prop.maxThreadsPerBlock );
        writef( "Max thread dimensions: (%d, %d, %d)\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2] );
        writef( "Max grid dimensions: (%d, %d, %d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2] );
	writef("\n");
    }
}

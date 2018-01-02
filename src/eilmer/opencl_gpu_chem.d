// Authors: KD, RJG & PJ
//

import core.stdc.stdlib;
import core.stdc.stdio;
import std.file;
import std.stdio;
import std.string;
import std.conv;

import fvcore: FlowSolverException;
import globalconfig;
import globaldata;
import fvcell;
import cl;
import gas;
import kinetics;
 
immutable string openCLProgName = "alpha_qss_kernel.cl";
immutable int MAX_SOURCE_SIZE = 0x100000;

struct OpenCLContainer {
    cl_int status;
    cl_uint numPlatforms;
    cl_platform_id* platforms;
    cl_uint numDevices;
    cl_device_id* devices;
    cl_context context;
    cl_command_queue cmdQueue;

    cl_mem bufkf, bufkb, bufMM;
    cl_mem bufY, bufYp, bufYpo, bufYc;
    cl_mem bufh;
    cl_mem bufalpha, bufq, bufp, bufqp, bufpp;
    cl_mem bufq_bar, bufp_bar, bufalpha_bar;
    cl_mem bufcond, bufsigma, bufdebugging;

    cl_program program;
    cl_kernel kernel;

    size_t nreac_datasize;
    size_t datasize;
    size_t ncell;
    size_t debugnum;

    this(size_t ncells, size_t nreacElements, size_t nspElements, size_t debugNum)
    {
        nreac_datasize = double.sizeof*nreacElements;
        datasize = double.sizeof*nspElements;
        this.ncell = ncells;
        debugnum = debugNum;
    }

    ~this()
    {
        clReleaseCommandQueue(cmdQueue);
        clReleaseMemObject(bufkf);
        clReleaseMemObject(bufkb);
        clReleaseMemObject(bufMM);
        clReleaseMemObject(bufY);
        clReleaseMemObject(bufYp);
        clReleaseMemObject(bufYpo);
        clReleaseMemObject(bufYc);
        clReleaseMemObject(bufh);
        clReleaseMemObject(bufalpha);
        clReleaseMemObject(bufq);
        clReleaseMemObject(bufp);
        clReleaseMemObject(bufqp);
        clReleaseMemObject(bufpp);
        clReleaseMemObject(bufq_bar);
        clReleaseMemObject(bufp_bar);
        clReleaseMemObject(bufalpha_bar);
        clReleaseMemObject(bufdebugging);
        clReleaseContext(context);
        clReleaseKernel(kernel);
        clReleaseProgram(program);
    }
};

class GPUChem {
public:
    this()
    {
        writeln("GPUChem:this()");
        _gmodel = GlobalConfig.gmodel_master;
        auto myChemUpdate = cast(ChemistryUpdate) gasBlocks[0].myConfig.thermochemUpdate;
        if (myChemUpdate !is null) { 
            _rmech = myChemUpdate.rmech.dup();
        } else {
            throw new Exception("Opps, incorrect ThermochemicalReactor.");
        }
        size_t nsp = _gmodel.n_species();
        size_t nreac = _rmech.n_reactions();
        _conc.length = nsp;

        storeReferencesToCells();
        _ncell = _cells.length;
        _nreacElements = _ncell*nreac;
        _nspElements = _ncell*nsp;
        int numDebug = 0;

        writeln("Init opencl container.");
        _oc = OpenCLContainer(_ncell, _nreacElements, _nspElements, numDebug);
        
        writeln("Init work arrays...");
        initWorkArrays();

        writeln("Configure opencl container");
        configureOpenCLContainer();

    }

    ~this()
    {
        free(_kf);
        free(_kb);
        free(_Y);
        free(_h);
    }

    /* We can only call this function after the blocks have
     * been configured in the main routines.
     */
    void storeReferencesToCells()
    {
        foreach (blk; gasBlocks) {
            foreach (cell; blk.cells) {
                _cells ~= cell;
            }
        }
    }

    void initWorkArrays()
    {
        _kf = cast(double*)malloc(_oc.nreac_datasize);
        _kb = cast(double*)malloc(_oc.nreac_datasize);
        _Y = cast(double*)malloc(_oc.datasize);
        _h = cast(double*)malloc(_oc.ncell*double.sizeof);
    }

    void configureOpenCLContainer()
    {
        //Find available devices----------------------------------------------------------------------------
        //retrieve number of platforms
        _oc.status = clGetPlatformIDs(0, null, &(_oc.numPlatforms));
        
        //allocate enough space for each platform
        _oc.platforms = cast(cl_platform_id*) malloc(_oc.numPlatforms*cl_platform_id.sizeof);
        
        //fill in platforms
        _oc.status = clGetPlatformIDs(_oc.numPlatforms, _oc.platforms, null);

        //retrieve number of devices, note for GPU change CPU to GPU
        _oc.status = clGetDeviceIDs(_oc.platforms[0], CL_DEVICE_TYPE_GPU, 0, null, &(_oc.numDevices));

        //allocate enough space for each device
        _oc.devices = cast(cl_device_id*)malloc(_oc.numDevices*cl_device_id.sizeof);

        //fill in the devices, for for GPU change CPU to GPU
        _oc.status=clGetDeviceIDs(_oc.platforms[0], CL_DEVICE_TYPE_GPU, _oc.numDevices, _oc.devices, null);

        //D. Create Context and Command Queue for each device for use---------------------------------------------
        _oc.context = clCreateContext(null, _oc.numDevices, _oc.devices, null, null, &(_oc.status));
        _oc.cmdQueue = clCreateCommandQueue(_oc.context, _oc.devices[0], 0, &(_oc.status));

        //E. Create buffers which will contain data on device-----------------------------------------------------
        _oc.bufkf = clCreateBuffer(_oc.context, CL_MEM_READ_ONLY, _oc.nreac_datasize, null, &(_oc.status));
        _oc.bufkb = clCreateBuffer(_oc.context, CL_MEM_READ_ONLY, _oc.nreac_datasize, null, &(_oc.status));
        _oc.bufMM = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.nreac_datasize, null, &(_oc.status));
        _oc.bufY = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufYp = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufYpo = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufYc = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufh = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.ncell*double.sizeof, null, &(_oc.status));
        _oc.bufalpha = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufq = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufp = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufqp = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufpp = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufq_bar = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufp_bar = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufalpha_bar = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.datasize, null, &(_oc.status));
        _oc.bufdebugging = clCreateBuffer(_oc.context, CL_MEM_READ_WRITE, _oc.debugnum*double.sizeof, null, &(_oc.status));

        //create a program with source code
        string fname = openCLProgName;
        FILE *fp = fopen(fname.toStringz, "r");
        if ( fp == null ) {
            throw new FlowSolverException("ERROR: Failed to open gpu chemistry kernel.");
        }
        char* programSource = cast(char*) malloc(MAX_SOURCE_SIZE);
        size_t sourceSize = fread(programSource, 1, MAX_SOURCE_SIZE, fp);
        fclose(fp);

        _oc.program = clCreateProgramWithSource(_oc.context, 1, cast(const char**)&programSource,
                                                cast(const size_t *)&sourceSize, &(_oc.status));
        //build (compile) the program for the device
        _oc.status = clBuildProgram(_oc.program, _oc.numDevices, _oc.devices, null, null, null);

        //create the kernel
        _oc.kernel = clCreateKernel(_oc.program, "new_qss", &(_oc.status));

        //Let's the user know some interesting facts ahout the particular GPU doing the calculations
        writeln("Information about chemistry accelerator device (GPU):");
        size_t  maxWorkgroupSize;
        clGetDeviceInfo(_oc.devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, maxWorkgroupSize.sizeof, &maxWorkgroupSize, null);
        writeln("maxWorkgroupSize= ", maxWorkgroupSize);

        cl_uint  maxComputeUnit;
        clGetDeviceInfo(_oc.devices[0], CL_DEVICE_MAX_COMPUTE_UNITS, maxComputeUnit.sizeof, &maxComputeUnit, null);
        writeln("maxComputeUnit= ", maxComputeUnit);

        cl_ulong  GlobeMemSize;
        clGetDeviceInfo(_oc.devices[0], CL_DEVICE_GLOBAL_MEM_SIZE, GlobeMemSize.sizeof, &GlobeMemSize, null);
        writeln("GlobeMemSize= ", GlobeMemSize);

        cl_char[1024] name;
        clGetDeviceInfo(_oc.devices[0], CL_DEVICE_NAME, name.sizeof, &name, null);
        //writeln("name= ", name);
    
        cl_char[1024]  vendor;
        clGetDeviceInfo(_oc.devices[0], CL_DEVICE_VENDOR, vendor.sizeof, &vendor, null);
        //writeln("vendor= ", vendor);
    }

    void thermochemical_increment(double dt_flow)
    {
        size_t nsp = _gmodel.n_species();
        size_t nreac = _rmech.n_reactions();

        prepareDataForGPU();
        writeDataToGPUBuffers();
        prepareKernelArgsForGPU(dt_flow);

        // Define an index space for work items
        size_t globalWorkSize = _ncell;

        // Execute the kernel
        _oc.status = clEnqueueNDRangeKernel(_oc.cmdQueue, _oc.kernel, 1, null, &globalWorkSize, null, 0 , null, null);
        // Read data from GPU
        clEnqueueReadBuffer(_oc.cmdQueue, _oc.bufY, CL_TRUE, 0, _oc.datasize, cast(void *)_Y, 0, null, null);
        clEnqueueReadBuffer(_oc.cmdQueue, _oc.bufh, CL_TRUE, 0, _oc.ncell*double.sizeof, cast(void *)_h, 0, null, null);

        // Now use the new concentrations to update the rest of the gas state
        foreach ( i, cell; _cells ) {
            foreach ( isp; 0 .. nsp ) _conc[isp] = _Y[i+isp*_ncell];
            _gmodel.conc2massf(_conc, cell.fs.gas);
            cell.dt_chem = _h[i];
            _gmodel.update_thermo_from_rhou(cell.fs.gas);
            if ( GlobalConfig.viscous ) _gmodel.update_trans_coeffs(cell.fs.gas);
            foreach ( isp; 0 .. nsp )
                cell.U[0].massf[isp] = cell.fs.gas.rho*cell.fs.gas.massf[isp];
        }
    }

    void prepareDataForGPU()
    {
        size_t nsp = _gmodel.n_species();
        size_t nreac = _rmech.n_reactions();

        foreach (i, cell; _cells) {
            // 1. Compute concentrations and store
            _gmodel.massf2conc(cell.fs.gas, _conc);
            foreach ( isp; 0 .. nsp ) _Y[i+isp*_ncell] = _conc[isp];

            // 2. Compute kf and kb and store
            _rmech.eval_rate_constants(cell.fs.gas);
            foreach ( int ir; 0 .. to!int(nreac) ) {
                _kf[i+ir*_ncell] = _rmech.k_f(ir);
                _kb[i+ir*_ncell] = _rmech.k_b(ir);
            }

            // 3. Finally, populate stepsize (h) array
            if ( cell.dt_chem <= 0.0 )
                _h[i] = _rmech.estimateStepSize(_conc);
            else
                _h[i] = cell.dt_chem;
        }
    }

    void writeDataToGPUBuffers()
    {
        _oc.status = clEnqueueWriteBuffer(_oc.cmdQueue, _oc.bufkf, CL_FALSE, 0, _oc.nreac_datasize, cast(void *)_kf, 0, null, null);
        _oc.status = clEnqueueWriteBuffer(_oc.cmdQueue, _oc.bufkb, CL_FALSE, 0, _oc.nreac_datasize, cast(void *)_kb, 0, null, null);
        _oc.status = clEnqueueWriteBuffer(_oc.cmdQueue, _oc.bufY, CL_FALSE, 0, _oc.datasize, cast(void *)_Y, 0, null, null);
        _oc.status = clEnqueueWriteBuffer(_oc.cmdQueue, _oc.bufh, CL_FALSE, 0, _oc.ncell*double.sizeof, cast(void*)_h, 0, null, null);
    }
    
    void prepareKernelArgsForGPU(double dt_global)
    {
        size_t nsp = _gmodel.n_species();
        size_t nreac = _rmech.n_reactions();
        _oc.status = clSetKernelArg(_oc.kernel, 0, cl_mem.sizeof, cast(void *)&(_oc.bufkf));
        _oc.status = clSetKernelArg(_oc.kernel, 1, cl_mem.sizeof, cast(void *)&(_oc.bufkb));
        _oc.status = clSetKernelArg(_oc.kernel, 2, cl_mem.sizeof, cast(void *)&(_oc.bufMM));
        _oc.status = clSetKernelArg(_oc.kernel, 3, cl_mem.sizeof, cast(void *)&(_oc.bufY));
        _oc.status = clSetKernelArg(_oc.kernel, 4, cl_mem.sizeof, cast(void *)&(_oc.bufYc));
        _oc.status = clSetKernelArg(_oc.kernel, 5, cl_mem.sizeof, cast(void *)&(_oc.bufYp));
        _oc.status = clSetKernelArg(_oc.kernel, 6, cl_mem.sizeof, cast(void *)&(_oc.bufh));
        _oc.status = clSetKernelArg(_oc.kernel, 7, cl_mem.sizeof, cast(void *)&(_oc.bufalpha));
        _oc.status = clSetKernelArg(_oc.kernel, 8, cl_mem.sizeof, cast(void *)&(_oc.bufq));
        _oc.status = clSetKernelArg(_oc.kernel, 9, cl_mem.sizeof, cast(void *)&(_oc.bufp));
        _oc.status = clSetKernelArg(_oc.kernel, 10, cl_mem.sizeof, cast(void *)&(_oc.bufqp));
        _oc.status = clSetKernelArg(_oc.kernel, 11, cl_mem.sizeof, cast(void *)&(_oc.bufpp));
        _oc.status = clSetKernelArg(_oc.kernel, 12, cl_mem.sizeof, cast(void *)&(_oc.bufq_bar));
        _oc.status = clSetKernelArg(_oc.kernel, 13, cl_mem.sizeof, cast(void *)&(_oc.bufp_bar));
        _oc.status = clSetKernelArg(_oc.kernel, 14, cl_mem.sizeof, cast(void *)&(_oc.bufalpha_bar));
        _oc.status = clSetKernelArg(_oc.kernel, 15, cl_mem.sizeof, cast(void *)&(_oc.bufYpo));
        _oc.status = clSetKernelArg(_oc.kernel, 16, cl_int.sizeof, cast(void *)&nsp);
        _oc.status = clSetKernelArg(_oc.kernel, 17, cl_int.sizeof, cast(void *)&nreac);
        _oc.status = clSetKernelArg(_oc.kernel, 18, cl_int.sizeof, cast(void *)&_ncell);
        _oc.status = clSetKernelArg(_oc.kernel, 19, cl_double.sizeof, cast(void *)&dt_global);
        _oc.status = clSetKernelArg(_oc.kernel, 20, cl_mem.sizeof, cast(void *)&(_oc.bufdebugging));
    }

private:
    FVCell[] _cells;
    GasModel _gmodel;
    ReactionMechanism _rmech;
    size_t _ncell;
    size_t _nreacElements;
    size_t _nspElements;
    double* _kf;
    double* _kb;
    double* _Y;
    double* _h;
    double[] _conc;

    OpenCLContainer _oc;
}

void initGPUChem()
{

    GlobalConfig.gpuChem = new GPUChem();
}

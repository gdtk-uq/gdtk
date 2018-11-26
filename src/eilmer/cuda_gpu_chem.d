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
import gas;
import kinetics.reaction_mechanism;
import kinetics.chemistry_update;

import cuda_d.cuda;
import cuda_d.cuda_runtime_api;
import cuda_d.cublas_api;

extern (C) void kernel_launcher(double *kf, double *kb, double *MM, double *Y, double *Yp, double *Yc,
                                double *h, double *alpha, double *q, double *p, double *qp, double *pp,
                                double *q_bar, double *p_bar, double *alpha_bar, double *Ypo, size_t nspec,
                                size_t numreac, size_t numcell, double dt_global, double *debugging,
                                size_t blkDimx, size_t blkDimy, size_t grdDimx, size_t grdDimy);

struct CudaContainer {
    double *bufkf; double *bufkb; double *bufMM;
    double *bufY; double *bufYp; double *bufYpo; double *bufYc;
    double *bufh;
    double *bufalpha; double *bufq; double *bufp; double *bufqp; double *bufpp;
    double *bufq_bar; double *bufp_bar; double *bufalpha_bar;
    double *bufcond; double *bufsigma; double *bufdebugging;

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
        cudaFree(bufkf);
        cudaFree(bufkb);
        cudaFree(bufMM);
        cudaFree(bufY);
        cudaFree(bufYp);
        cudaFree(bufYpo);
        cudaFree(bufYc);
        cudaFree(bufh);
        cudaFree(bufalpha);
        cudaFree(bufq);
        cudaFree(bufp);
        cudaFree(bufqp);
        cudaFree(bufpp);
        cudaFree(bufq_bar);
        cudaFree(bufp_bar);
        cudaFree(bufalpha_bar);
        cudaFree(bufdebugging);
    }
};

class GPUChem {
public:
    this()
    {
        writeln("GPUChem:this()");
        _gmodel = GlobalConfig.gmodel_master;
        auto myChemUpdate = cast(ChemistryUpdate) localFluidBlocks[0].myConfig.thermochemUpdate;
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

        writeln("Init cuda container.");
        _cc = CudaContainer(_ncell, _nreacElements, _nspElements, numDebug);
        
        writeln("Init work arrays...");
        initWorkArrays();

        writeln("Configure cuda container");
        configureCudaContainer();

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
        foreach (blk; localFluidBlocks) {
            foreach (cell; blk.cells) {
                _cells ~= cell;
            }
        }
    }

    void initWorkArrays()
    {
        _kf = cast(double*)malloc(_cc.nreac_datasize);
        _kb = cast(double*)malloc(_cc.nreac_datasize);
        _Y = cast(double*)malloc(_cc.datasize);
        _h = cast(double*)malloc(_cc.ncell*double.sizeof);
    }

    void configureCudaContainer()
    {
        //E. Create buffers which will contain data on device-----------------------------------------------------
        cudaMalloc( cast(void**)&_cc.bufkf, _cc.nreac_datasize);
        cudaMalloc( cast(void**)&_cc.bufkb, _cc.nreac_datasize);
        cudaMalloc( cast(void**)&_cc.bufMM, _cc.nreac_datasize);
        cudaMalloc( cast(void**)&_cc.bufY, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufYp, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufYpo, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufYc, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufh, _cc.ncell*double.sizeof);
        cudaMalloc( cast(void**)&_cc.bufalpha, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufq, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufp, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufqp, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufpp, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufq_bar, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufp_bar, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufalpha_bar, _cc.datasize);
        cudaMalloc( cast(void**)&_cc.bufdebugging, _cc.debugnum*double.sizeof);
 
    }

    void thermochemical_increment(double dt_flow)
    {
        size_t nsp = _gmodel.n_species();
        size_t nreac = _rmech.n_reactions();

        prepareDataForGPU();
        writeDataToGPUBuffers();
        //prepareKernelArgsForGPU(dt_flow);

        // Define an index space for work items
        size_t blkDimx = 1;
        size_t blkDimy = 1;
        size_t grdDimx = _ncell;
        size_t grdDimy = 1;

        // Execute the kernel
        kernel_launcher(_cc.bufkf, _cc.bufkb, _cc.bufMM, _cc.bufY, _cc.bufYp, _cc.bufYc, _cc.bufh, _cc.bufalpha,
                        _cc.bufq, _cc.bufp, _cc.bufqp, _cc.bufpp, _cc.bufq_bar, _cc.bufp_bar, _cc.bufalpha_bar,
                        _cc.bufYpo, nsp, nreac, _ncell, dt_flow, _cc.bufdebugging, blkDimx,
                        blkDimy, grdDimx, grdDimy);
        
        // Read data from GPU
        cudaMemcpy( cast(void*)_Y, cast(void*)_cc.bufY, _cc.datasize, cudaMemcpyKind.cudaMemcpyDeviceToHost );
        cudaMemcpy( cast(void*)_h, cast(void*)_cc.bufh, _cc.ncell*double.sizeof, cudaMemcpyKind.cudaMemcpyDeviceToHost );
        
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
        cudaMemcpy( cast(void*)_cc.bufkf, cast(void*)_kf, _cc.nreac_datasize, cudaMemcpyKind.cudaMemcpyHostToDevice );
        cudaMemcpy( cast(void*)_cc.bufkb, cast(void*)_kb, _cc.nreac_datasize, cudaMemcpyKind.cudaMemcpyHostToDevice );
        cudaMemcpy( cast(void*)_cc.bufY, cast(void*)_Y, _cc.datasize, cudaMemcpyKind.cudaMemcpyHostToDevice );
        cudaMemcpy( cast(void*)_cc.bufh, cast(void*)_h, _cc.ncell*double.sizeof, cudaMemcpyKind.cudaMemcpyHostToDevice );
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

    CudaContainer _cc;
}

void initGPUChem()
{

    GlobalConfig.gpuChem = new GPUChem();
}

/*
  Densified datastructures for the fluid solver.
 
  Author: Nick Gibbons
  Version: 2024-04-22

*/

module lmr.coredata;

import geom;
import nm.number;

import lmr.conservedquantities;
import lmr.flowgradients;
import lmr.flowstate;
import lmr.lsqinterp;
import lmr.globalconfig : LocalConfig;
import lmr.onedinterp : L2R2InterpData, L3R3InterpData;

struct FluidCellData{
    size_t[] all_cell_idxs;
    size_t[] nfaces;
    size_t[][] c2f;
    size_t[][] c2v;
    double[] dt_local;
    int[][] outsigns;
    bool[] data_is_bad;
    size_t[][] halo_cell_ids;
    size_t[][] halo_face_ids;
    number[] areas;
    number[] volumes;
    Vector3[] positions;
    number[3][] lengths;
    number[] wall_distances;
    bool[] in_turbulent_zone;
    size_t[][] cell_cloud_indices;
    Vector3[][] face_distances;
    FlowState[] flowstates;
    FlowGradients[] gradients;
    ConservedQuantities U0, U1, U2, U3, U4;
    ConservedQuantities dUdt0, dUdt1, dUdt2, dUdt3, dUdt4;
    ConservedQuantities source_terms;
    WLSQGradWorkspace[] workspaces;
    LSQInterpWorkspace[] lsqws;
    LSQInterpGradients[] lsqgradients;
    version(newton_krylov) {
        double[] cell_jacobians;
        FlowGradients[] saved_gradients;
        LSQInterpGradients[] saved_lsqgradients;
        ConservedQuantities saved_source_terms;
        bool[] doNotPerturb;
    }

    void allocate_minimal(size_t ncells, size_t nghost, LocalConfig myConfig){
    /*
        FluidFVCells expect a certain some things in the FVCellData to be built
        in order for their constructors to work properly. This routine
        allocates only those things; principly for creating ghost fluid cells
        in the solid_gas_full_face_copy code.

        TODO: Move the rest of the allocation to this file from fluidblock.d
        @author: Nick N. Gibbons
    */
        auto gmodel = myConfig.gmodel;
        size_t neq = myConfig.cqi.n;
        size_t nturb = myConfig.turb_model.nturb;
        size_t nftl = myConfig.n_flow_time_levels;
        size_t ngtl = myConfig.n_grid_time_levels;

        areas.length     = ngtl*(ncells + nghost);
        volumes.length   = ngtl*(ncells + nghost);
        positions.length = ngtl*(ncells + nghost);
        lengths.length   = ncells + nghost;

        U0.length = (ncells + nghost)*neq*nftl;
        if (nftl>1) U1.length = (ncells + nghost)*neq*nftl;
        if (nftl>2) U2.length = (ncells + nghost)*neq*nftl;
        if (nftl>3) U3.length = (ncells + nghost)*neq*nftl;
        if (nftl>4) U4.length = (ncells + nghost)*neq*nftl;
        dUdt0.length = (ncells + nghost)*neq*nftl;
        if (nftl>1) dUdt1.length = (ncells + nghost)*neq*nftl;
        if (nftl>2) dUdt2.length = (ncells + nghost)*neq*nftl;
        if (nftl>3) dUdt3.length = (ncells + nghost)*neq*nftl;
        if (nftl>4) dUdt4.length = (ncells + nghost)*neq*nftl;
        source_terms.length = (ncells + nghost)*neq;

        flowstates.reserve(ncells + nghost);
        gradients.reserve(ncells + nghost);
        workspaces.reserve(ncells + nghost);
        foreach (n; 0 .. ncells+nghost) flowstates ~= FlowState(gmodel, nturb);
        foreach (n; 0 .. ncells+nghost) gradients ~= FlowGradients(myConfig);
        foreach (n; 0 .. ncells+nghost) workspaces ~= WLSQGradWorkspace();
    }
}


struct LR {size_t left,right;}
struct LLLRRR {size_t L2,L1,L0,R0,R1,R2;}

struct FVInterfaceData{
    size_t[] all_face_idxs;
    LR[] f2c;
    Vector3[] dL;
    Vector3[] dR;
    number[] areas;
    Vector3[] normals;
    Vector3[] tangents1;
    Vector3[] tangents2;
    Vector3[] positions;
    Vector3[] grid_velocities;
    FlowState[] flowstates;
    LLLRRR[] stencil_idxs;
    L2R2InterpData[] l2r2_interp_data;
    L3R3InterpData[] l3r3_interp_data;
    FlowGradients[] gradients;
    WLSQGradWorkspace[] workspaces;
    ConservedQuantities fluxes;
    bool[] left_interior_only;
    bool[] right_interior_only;
}

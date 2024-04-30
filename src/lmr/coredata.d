/*
  Densified datastructures for the fluid solver.
 
  Author: Nick Gibbons
  Version: 2024-04-22

*/

module lmr.coredata;

import nm.number;
import geom;
import flowstate;
import flowgradients;
import conservedquantities;
import lsqinterp;
import onedinterp : L2R2InterpData, L3R3InterpData;

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

/*
  Densified datastructures for the fluid solver.
 
  Author: Nick Gibbons
  Version: 2024-04-22

*/

module lmr.coredata;

import nm.number;
import geom;
//import gas;
//import kinetics;
import flowstate;
import flowgradients;
import conservedquantities;
//import globalconfig;
import lsqinterp;

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

// face.chpl
// Face functions (such as flux calculators) for Chicken-in-Chapel.
//
// PJ 2024-01-30

module Face {
  use Math;
  use Geom;
  use Gas;
  use Flow;

  // Boundary condition codes, to decide what to do for the ghost cells.
  // Periodic boundary conditions should just work if we wrap the index
  // in each direction.
  // There's not enough information here to have arbitrary block connections.
  enum BCCode {
    wall_with_slip,
    wall_no_slip_adiabatic,
    wall_no_slip_fixed_T,
    exchange,
    inflow,
    outflow,
    inflow_function
  };

  var BCName = [
    BCCode.wall_with_slip => "wall_with_slip",
    BCCode.wall_no_slip_adiabatic => "wall_no_slip_adiabatic",
    BCCode.wall_no_slip_fixed_T => "wall_no_slip_fixed_T",
    BCCode.exchange => "exchange",
    BCCode.inflow => "inflow",
    BCCode.outflow => "outflow",
    BCCode.inflow_function => "inflow_function"
  ];

  proc BCCode_from_name(name: string): BCCode {
    select name {
      when "wall_with_slip" do return BCCode.wall_with_slip;
      when "wall_no_slip_adiabatic" do return BCCode.wall_no_slip_adiabatic;
      when "wall_no_slip_fixed_T" do return BCCode.wall_no_slip_fixed_T;
      when "wall_no_slip" do return BCCode.wall_no_slip_adiabatic; // alias
      when "exchange" do return BCCode.exchange;
      when "inflow" do return BCCode.inflow;
      when "outflow" do return BCCode.outflow;
      when "inflow_function" do return BCCode.inflow_function;
      otherwise return BCCode.wall_with_slip;
    }
  }

  enum BCFunc {
    nofun,
    supersonic_vortex,
    laminar_boundary_layer,
    manufactured_solution
  };

  var BCFuncName = [
    BCFunc.nofun => "none",
    BCFunc.supersonic_vortex => "supersonic_vortex",
    BCFunc.laminar_boundary_layer => "laminar_boundary_layer",
    BCFunc.manufactured_solution => "manufactured_solution"
  ];

  proc BCFunc_from_name(name: string): BCFunc {
    select name {
      when "nofun" do return BCFunc.nofun;
      when "none" do return BCFunc.nofun; // alias from original Chicken
      when "supersonic_vortex" do return BCFunc.supersonic_vortex;
      when "laminar_boundary_layer" do return BCFunc.laminar_boundary_layer;
      when "manufactured_solution" do return BCFunc.manufactured_solution;
      otherwise return BCFunc.nofun;
    }
  }

  // Interpolation functions that will be used in the fluc calculator.

  // A smooth slope limiter.
  proc van_albada_limit1(a: real, b: real): real {
    const eps = 1.0e-12;
    var s = (a*b + Math.abs(a*b) + eps)/(a*a + b*b + eps);
    return s;
  }

  // Reconstruct values, qL,qR, at the middle interface
  // for a stencil of 4 cell-centred values.
  // Assume equal cell widths.
  proc interp_l2r2_scalar(qL1: real, qL0: real, qR0: real, qR1: real,
                          out qL: real, out qR: real) {
    // Set up differences and limiter values.
    const delLminus = (qL0 - qL1);
    const del = (qR0 - qL0);
    const delRplus = (qR1 - qR0);
    const sL = van_albada_limit1(delLminus, del);
    const sR = van_albada_limit1(del, delRplus);
    // The actual high-order reconstruction, possibly limited.
    qL = qL0 + sL * 0.125 * (3.0*del + delLminus);
    qR = qR0 - sR * 0.125 * (delRplus + 3.0*del);
  }

  // Flux calculations are done in the context of a face on a cell.

  param cloud_ncmax = 2;
  param cloud_nfmax = 8;
  param cloud_nmax = cloud_ncmax + cloud_nfmax;

  record FVFace {
    var pos: Vector3; // midpoint position in space
    var area: real;
    var n: Vector3;  // unit normal
    var t1: Vector3; // unit tangent 1
    var t2: Vector3 ; // unit tangent 2
    var F: [Flow.DCon]real; // flux vector for conserved quantities
    // We will keep connections to the pieces composing the face
    // as indices into global arrays.
    var vtx: [0..#4]int = [-1,-1,-1,-1];
    var left_cells: [0..#2]int = [-1,-1];
    var right_cells: [0..#2]int = [-1,-1];
    // To apply boundary conditions at the face, we need to carry some extra information.
    // Not all faces are boundary faces, so we start with dummy values.
    var bcId: int = -1;
    var bcCode: BCCode = BCCode.wall_with_slip;
    var other_blkId: int = -1;
    var other_cells: [0..#2]int = [-1,-1];
    var inflowId: int = -1;
    var TWall: real = 300.0;
    var bcFun: BCFunc = BCFunc.nofun;
    // For the gradient calculations that form part of the viscous fluxes
    // we keep lists of faces and cells that form a cloud of points around
    // this face-centre.
    var cells_in_cloud: [0..#cloud_ncmax]int = [-1,-1];
    var faces_in_cloud: [0..#cloud_nfmax]int = [-1,-1,-1,-1,-1,-1,-1,-1];
    var cloud_nc: int = 0;
    var cloud_nf: int = 0;
    // Prepared least-squares solution for cloud of cell- and face-FlowStates.
    var wx, wy, wz: [0..#cloud_nmax]real;
    // We also need the FlowState at this face-centre.  It will be set during
    // the convective-flux calculation or by the boundary-condition code for a wall.
    var fs: FlowState;

    // Specific convective-flux calculators here...

    // Compute the face's flux vector from left and right flow states.
    // Wada and Liou's flux calculator, implemented from details in their AIAA paper,
    // with hints from Ian Johnston.
    // Y. Wada and M. -S. Liou (1994)
    // A flux splitting scheme with high-resolution and robustness for discontinuities.
    // AIAA-94-0083.
    proc ref ausmdv(const ref fsL: FlowState, const ref fsR: FlowState) {
      var velL: Vector3 = fsL.vel;
      var velR: Vector3 = fsR.vel;
      velL.transformToLocalFrame(n, t1, t2);
      velR.transformToLocalFrame(n, t1, t2);
      //
      var rhoL = fsL.gs.rho;
      var pL = fsL.gs.p;
      var pLrL = pL/rhoL;
      var velxL = velL.x;
      var velyL = velL.y;
      var velzL = velL.z;
      var uL = fsL.gs.e;
      var aL = fsL.gs.a;
      var keL = 0.5*(velxL*velxL + velyL*velyL + velzL*velzL);
      var HL = uL + pLrL + keL;
      //
      var rhoR = fsR.gs.rho;
      var pR = fsR.gs.p;
      var pRrR = pR/rhoR;
      var velxR = velR.x;
      var velyR = velR.y;
      var velzR = velR.z;
      var uR = fsR.gs.e;
      var aR = fsR.gs.a;
      var keR = 0.5*(velxR*velxR + velyR*velyR + velzR*velzR);
      var HR = uR + pR/rhoR + keR;
      //
      // This is the main part of the flux calculator.
      //
      // Weighting parameters (eqn 32) for velocity splitting.
      var alphaL = 2.0*pLrL/(pLrL+pRrR);
      var alphaR = 2.0*pRrR/(pLrL+pRrR);
      // Common sound speed (eqn 33) and Mach numbers.
      var am = max(aL, aR);
      var ML = velxL/am;
      var MR = velxR/am;
      // Left state:
      // pressure splitting (eqn 34)
      // and velocity splitting (eqn 30)
      var pLplus, velxLplus: real;
      var dvelxL = 0.5*(velxL+abs(velxL));
      if abs(ML) <= 1.0 {
        pLplus = pL*(ML+1.0)*(ML+1.0)*(2.0-ML)*0.25;
        velxLplus = alphaL*((velxL+am)*(velxL+am)/(4.0*am) - dvelxL) + dvelxL;
      } else {
        pLplus = pL * dvelxL / velxL;
        velxLplus = dvelxL;
      }
      // Right state:
      // pressure splitting (eqn 34)
      // and velocity splitting (eqn 31)
      var pRminus, velxRminus: real;
      var dvelxR = 0.5*(velxR-abs(velxR));
      if abs(MR) <= 1.0 {
        pRminus = pR*(MR-1.0)*(MR-1.0)*(2.0+MR)*0.25;
        velxRminus = alphaR*(-(velxR-am)*(velxR-am)/(4.0*am) - dvelxR) + dvelxR;
      } else {
        pRminus = pR * dvelxR / velxR;
        velxRminus = dvelxR;
      }
      // The mass flux. (eqn 29)
      var massL = velxLplus*rhoL;
      var massR = velxRminus*rhoR;
      var mass_half = massL+massR;
      // Pressure flux (eqn 34)
      var p_half = pLplus + pRminus;
      // Momentum flux: normal direction
      // Compute blending parameter s (eqn 37),
      // the momentum flux for AUSMV (eqn 21) and AUSMD (eqn 21)
      // and blend (eqn 36).
      var dp = pL - pR;
      const K_SWITCH = 10.0;
      dp = K_SWITCH * abs(dp) / min(pL, pR);
      var s = 0.5 * min(1.0, dp);
      var rvel2_AUSMV = massL*velxL + massR*velxR;
      var rvel2_AUSMD = 0.5*(mass_half*(velxL+velxR) - abs(mass_half)*(velxR-velxL));
      var rvel2_half = (0.5+s)*rvel2_AUSMV + (0.5-s)*rvel2_AUSMD;
      //
      // Assemble components of the flux vector (eqn 36).
      var momentum: Vector3;
      if mass_half >= zero {
        forall i in DS do F[mass+i] = mass_half * fsL.gs.massf[i];
        momentum.set(rvel2_half+p_half, mass_half*fsL.vel.y, mass_half*fsL.vel.z);
        F[totEnergy] = mass_half*HL;
      } else {
        forall i in DS do F[mass+i] = mass_half * fsR.gs.massf[i];
        momentum.set(rvel2_half+p_half, mass_half*fsR.vel.y, mass_half*fsR.vel.z);
        F[totEnergy] = mass_half*HR;
      }
      momentum.transformToGlobalFrame(n, t1, t2);
      F[xMom] = momentum.x;
      F[yMom] = momentum.y;
      F[zMom] = momentum.z;
      return;
    } // end ausmdv()

    // [TODO] other flux calculators

    // And one generic flux calculation function.

    proc ref calculate_convective_flux(
      const ref fsL1: FlowState, const ref fsL0: FlowState,
      const ref fsR0: FlowState, const ref fsR1: FlowState,
      const x_order: int, const gm: GasModel)
    {
      // First-order reconstruction is just a copy from the nearest cell centre.
      var fsL: FlowState = fsL0;
      var fsR: FlowState = fsR0;
      if (x_order > 1) {
        // We will interpolate only some GasState properties...
        interp_l2r2_scalar(fsL1.gs.rho, fsL0.gs.rho, fsR0.gs.rho,
                           fsR1.gs.rho, fsL.gs.rho, fsR.gs.rho);
        interp_l2r2_scalar(fsL1.gs.e, fsL0.gs.e, fsR0.gs.e,
                           fsR1.gs.e, fsL.gs.e, fsR.gs.e);
        forall i in DS {
          interp_l2r2_scalar(fsL1.gs.massf[i], fsL0.gs.massf[i], fsR0.gs.massf[i],
                             fsR1.gs.massf[i], fsL.gs.massf[i], fsR.gs.massf[i]);
        }
        // and make the rest consistent.
        gm.update_from_rhoe(fsL.gs);
        gm.update_from_rhoe(fsR.gs);
        // Velocity components.
        interp_l2r2_scalar(fsL1.vel.x, fsL0.vel.x, fsR0.vel.x,
                           fsR1.vel.x, fsL.vel.x, fsR.vel.x);
        interp_l2r2_scalar(fsL1.vel.y, fsL0.vel.y, fsR0.vel.y,
                           fsR1.vel.y, fsL.vel.y, fsR.vel.y);
        interp_l2r2_scalar(fsL1.vel.z, fsL0.vel.z, fsR0.vel.z,
                           fsR1.vel.z, fsL.vel.z, fsR.vel.z);
      }
      // Use the reconstructed values near the face in a simple flux calculator.
      ausmdv(fsL, fsR);
      // For later use in gradient calculations for viscous fluxes.
      fs.setAsAverage(fsL,fsR);
    } // end calculate_convective_flux()

  } // end record FVFace

} // end module Face

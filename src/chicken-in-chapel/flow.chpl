// flow.chpl
// FlowState and ConservedQuantities for Chicken-in-Chapel.
//
// PJ 2024-01-29

module Flow {
  use Geom;
  use Gas;

  // For sizing and indexing into the Flux vectors of conserved quantities.
  param mass: int = 0;               // starting index for species densities
  param xMom: int = mass + nSpecies; // index for x-momentum
  param yMom: int = xMom + 1;
  param zMom: int = yMom + 1;
  param totEnergy: int = zMom + 1;
  param nCon: int = nSpecies + 4;    // total number of elements
  const DCon = {0..nCon-1};

  record FlowState {
    var gs: GasState;
    var vel: Vector3;

    proc setAsAverage(a: FlowState, b: FlowState) {
      gs.setAsAverage(a.gs, b.gs);
      vel.setAsAverage(a.vel, b.vel);
    }

    proc ref encodeConserved(ref U: [DCon]real) {
      var rho = gs.rho;
      forall i in DS do U[mass+i] = rho * gs.massf[i];
      U[xMom] = rho * vel.x;
      U[yMom] = rho * vel.y;
      U[zMom] = rho * vel.z;
      U[totEnergy] = rho * (gs.e + 0.5*vel.dot(vel));
    } // end encodeConserved()

    proc ref decodeConserved(const ref U: [DCon]real, gm: GasModel) {
      // [TODO] Check for NaNs.
      var rho = 0.0;
      for i in DS do rho += U[mass+i];
      // [TODO] or, maybe, just check rho greater then zero.
      var dinv = 1.0/rho;
      vel.set(U[xMom]*dinv, U[yMom]*dinv, U[zMom]*dinv);
      var e = U[totEnergy] * dinv;
      var ke = 0.5*vel.dot(vel);
      e -= ke;
      // Put data into GasState.
      gs.rho = rho;
      gs.e = e;
      for i in DS do gs.massf[i] = U[mass+i] * dinv;
      gm.update_from_rhoe(gs);
    } // end decodeConserved()

  } // end record FlowState

} // end module Flow

// flow.chpl
// FlowState and ConservedQuantities for Chicken-in-Chapel.
//
// PJ 2024-01-29

module Flow {
  use Geom;
  use Gas;

  class FlowError: Error { }

  // For sizing and indexing into the Flux vectors of conserved quantities.
  // Note that speciesCap is set at compile time and the number of species
  // in the gas model needs to be at or below this capacity.
  param mass: int = 0;                 // starting index for species densities
  param xMom: int = mass + speciesCap; // index for x-momentum
  param yMom: int = xMom + 1;
  param zMom: int = yMom + 1;
  param totEnergy: int = zMom + 1;
  param nCon: int = speciesCap + 4;    // total number of elements
  const DCon = {0..nCon-1};

  record FlowState {
    var gs: GasState;
    var vel: Vector3;

    proc set(const ref gs: GasState, const ref vel: Vector3) {
      this.gs = gs;
      this.vel = vel;
    }

    proc set(const ref fs: FlowState) {
      this.gs = fs.gs;
      this.vel = vel;
    }

    proc setAsAverage(const ref a: FlowState, const ref b: FlowState) {
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

    proc ref decodeConserved(const ref U: [DCon]real, const ref gm: GasModel) throws {
      var rho = 0.0;
      for i in DS do rho += U[mass+i];
      if !(rho > 0.0) {
        writeln("Density not positive.");
        throw new owned FlowError();
      }
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

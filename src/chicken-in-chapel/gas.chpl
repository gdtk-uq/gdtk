// gas.chpl
// Gas models and states for Chicken-in-Chapel.
//
// PJ 2024-01-28: Just enough to implement single-temperature models.
//
// [TODO]: consider initializing the gas model objects explicitly,
// reading their parameter values from a JSON string.
//

module Gas {
  use Math;

  // We decide the number of chemical species at compile time.
  config param nSpecies: int = 1;
  const DS = {0..nSpecies-1};


  record GasState {
    var rho: real;        // density, kg/m^3
    var massf: [DS]real;  // species mass-fractions
    var e: real;          // internal energy, J/kg
    var p: real;          // pressure, Pa
    var T: real;          // Trans-rotational temperature, degrees K
    var a: real;          // sound speed, m/s

    proc setAsAverage(ref gs0: GasState, ref gs1: GasState) {
      rho = 0.5*(gs0.rho + gs1.rho);
      forall i in DS do massf[i] = 0.5*(gs0.massf[i] + gs1.massf[i]);
      e = 0.5*(gs0.e + gs1.e);
      p = 0.5*(gs0.p + gs1.p);
      T = 0.5*(gs0.T + gs1.T);
      a = 0.5*(gs0.a + gs1.a);
      // Note, to get a consistent GasState, we need to call
      // either update_from_pT() or update_from_rhoe() on return.
    }
  }


  class GasModel {
    proc update_from_pT(ref gs: GasState) {
      // need to override
    }

    proc update_from_rhoe(ref gs: GasState) {
      // need to override
    }

    proc trans_coeffs(ref gs: GasState): (real, real) {
      // need to override
    }

    proc update_chemistry(ref gs: GasState, dt: real) {
      // need to override
    }
  } // end class GasModel


  class IdealGas: GasModel {
    const name: string = "IdealGas";
    //
    const g: real = 1.4;   // ratio of specific heats
    const R: real = 287.1; // gas constant, J/kg/K
    const Cv: real = R / (g-1.0);
    const Cp: real = Cv + R;
    //
    // Sutherland viscosity model, air.
    const mu_ref = 1.716e-5; // Pa.s
    const T_ref = 273.0;     // K
    const S = 111.0;         // K
    //
    const Prandtl: real = 0.72;

    override proc update_from_pT(ref gs: GasState) {
      gs.e = gs.T*Cv;
      gs.rho = gs.p / (gs.T * R);
      gs.a = sqrt(g * R * gs.T);
    }

    override proc update_from_rhoe(ref gs: GasState) {
      gs.T = gs.e/Cv;
      gs.p = gs.rho * R* gs.T;
      gs.a = sqrt(g * R * gs.T);
    }

    override proc trans_coeffs(ref gs: GasState): (real, real) {
      var mu = mu_ref*sqrt(gs.T/T_ref)*(gs.T/T_ref)*(T_ref + S)/(gs.T + S);
      var k = Cp*mu/Prandtl;
      return (mu, k);
    }

    override proc update_chemistry(ref gs: GasState, dt: real) {
      // Do nothing to change mass fractions.
      this.update_from_rhoe(gs);
    }
  } // end class IdealGas


  class ABReactingGas: GasModel {
    // Reacting gas, species A to species B (with YB=massf[0])
    const name: string = "ABReactingGas";

    const g: real = 6.0/5.0;      // ratio of specific heats
    const R: real = 287.0;        // gas constant, J/kg/K
    const q: real = 300_000.0;    // heat of reaction, J/kg
    const alpha: real = 1000.0;   // reaction rate, 1/s
    const Ti: real = 362.58;      // ignition temperature, K
    const Cv: real = R / (g-1.0);
    const Cp: real = Cv + R;
    //
    // Sutherland viscosity model, air.
    const mu_ref = 1.716e-5; // Pa.s
    const T_ref = 273.0;     // K
    const S = 111.0;         // K
    //
    const Prandtl: real = 0.72;

    override proc update_from_pT(ref gs: GasState) {
      const YB = gs.massf[0];
      gs.e = gs.T*Cv - YB*q;
      gs.rho = gs.p / (gs.T * R);
      gs.a = sqrt(g * R * gs.T);
    }

    override proc update_from_rhoe(ref gs: GasState) {
      const YB = gs.massf[0];
      gs.T = (gs.e + YB*q)/Cv;
      gs.p = gs.rho * R* gs.T;
      gs.a = sqrt(g * R * gs.T);
    }

    override proc trans_coeffs(ref gs: GasState): (real, real) {
      var mu = mu_ref*sqrt(gs.T/T_ref)*(gs.T/T_ref)*(T_ref + S)/(gs.T + S);
      var k = Cp*mu/Prandtl;
      return (mu, k);
    }

    override proc update_chemistry(ref gs: GasState, dt: real) {
      var YB = gs.massf[0];
      const YA = 1.0 - YB;
      if gs.T > Ti {
        YB = 1.0 - YA*exp(-alpha*dt);
      }
      gs.massf[0] = YB;
      this.update_from_rhoe(gs);
    }
  } // end class ABReactingGas

} // end module Gas

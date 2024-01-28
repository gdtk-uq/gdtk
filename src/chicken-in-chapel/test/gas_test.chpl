// gas_test.chpl
// Try the gas models and states for Chicken-in-Chapel.
//
// Build with:
// $ chpl gas_test.chpl ../gas.chpl
//
// Run with:
// $ ./gas_test
//
// PJ 2024-01-28: Set up the initial code using classes for GasModels
//

use Gas;

proc main() {
  var gm1 = new IdealGas();
  writeln("gm1=", gm1);
  var gs1 = new GasState(p=100.0e3, T=300.0);
  gm1.update_from_pT(gs1);
  writeln("gs1=", gs1);
  var (mu, k) = gm1.trans_coeffs(gs1);
  writeln("mu=", mu, " k=", k);

  var gm2 = new ABReactingGas(Ti=362.58);  // explicitly set something
  writeln("gm2=", gm2);
  var gs2 = new GasState(p=100.0e3, T=400.0, massf=[0.0,]);
  gm2.update_from_pT(gs2);
  writeln("before chemistry update, gs2=", gs2);
  gm2.update_chemistry(gs2, 1.0e-3);
  writeln("after chemistry_update, gs2=", gs2);
}

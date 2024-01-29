// flow_test.chpl
// Try the FlowState module for Chicken-in-Chapel.
//
// Build with:
// $ chpl flow_test.chpl ../vector3.chpl ../gas.chpl ../flow.chpl
//
// Run with:
// $ ./flow_test
//
// PJ 2024-01-29
//

use Geom;
use Gas;
use Flow;

proc main() {
  var gm1 = new IdealGas();
  writeln("gm1=", gm1);
  // We assume that speciesCap == 2
  var gs1 = new GasState(p=100.0e3, T=300.0, massf=[1.0,0.0]);
  gm1.update_from_pT(gs1);
  writeln("gs1=", gs1);
  var (mu, k) = gm1.trans_coeffs(gs1);
  writeln("mu=", mu, " k=", k);
  var v1 = new Vector3(1.0, 2.0, 3.0);
  var fs1 = new FlowState(gs=gs1, vel=v1);
  writeln("fs1=", fs1);
  var U1: [DCon]real;
  fs1.encodeConserved(U1);
  writeln("U1=", U1);
  var fs2 = new FlowState();
  fs2.decodeConserved(U1, gm1:GasModel);
  writeln("fs2=", fs2);
  // gm1.checkMassFractions(fs2.gs); // oops, dead value because ownership transferred above
}

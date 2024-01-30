// face_test.chpl
// Try the Face module for Chicken-in-Chapel.
//
// Build with:
// $ chpl face_test.chpl ../face.chpl ../flow.chpl ../gas.chpl ../vector3.chpl
//
// Run with:
// $ ./face_test
//
// PJ 2024-01-30
//

use Geom;
use Gas;
use Flow;
use Face;

proc main() {
  var bc1 = BCCode_from_name("exchange");
  writeln("bc1=", bc1, " name=", BCName[bc1]);
  var bcf1 = BCFunc_from_name("none");
  writeln("bcf1=", bcf1, " name=", BCFuncName[bcf1]);

  var a = 1.0; var b = 1.2; var s = van_albada_limit1(a, b);
  writeln("a=", a, " b=", b, " van_albada_limit1=", s);
  a = 1.0; b = -1.0; s = van_albada_limit1(a, b);
  writeln("a=", a, " b=", b, " van_albada_limit1=", s);

  var qL, qR;
  interp_l2r2_scalar(1.0, 2.0, 3.0, 4.0, qL, qR);
  writeln("1 2 3 4 --> qL=", qL, " qR=", qR);
  interp_l2r2_scalar(1.0, 2.0, 2.0, 1.0, qL, qR);
  writeln("1 2 2 1 --> qL=", qL, " qR=", qR);

  var gm1 = new IdealGas();
  writeln("gm1=", gm1);
  // We assume that speciesCap == 2
  var gs1 = new GasState(p=100.0e3, T=300.0, massf=[1.0,0.0]);
  gm1.update_from_pT(gs1);
  var v1 = new Vector3(1000.0, 10.0, 1.0);
  var fs1 = new FlowState(gs=gs1, vel=v1);
  writeln("fs1=", fs1);

  var f: FVFace;
  f.n.set(1.0,0.0,0.0); f.t1.set(0.0,1.0,0.0); f.t2.set(0.0,0.0,1.0);
  f.calculate_convective_flux(fs1, fs1, fs1, fs1, 2, gm1: borrowed GasModel);
  writeln("FVFace f=", f);
}

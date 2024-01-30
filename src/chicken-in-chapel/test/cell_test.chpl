// cell_test.chpl
// Try the Cell module for Chicken-in-Chapel.
//
// Build with:
// $ chpl cell_test.chpl ../cell.chpl ../face.chpl ../flow.chpl ../gas.chpl ../vector3.chpl
//
// Run with:
// $ ./cell_test
//
// PJ 2024-01-30
//

use Geom;
use Gas;
use Flow;
use Face;
use Cell;

proc main() {
  writeln("Test Cell...");
  var f1 = FaceId_from_name("kplus");
  writeln("f1=", f1, " name=", FaceIdName[f1]);
  var s1 = SourceTerms_from_name("none");
  writeln("s1=", s1, " name=", SourceTermsName[s1]);
  for i in DS do writeln(speciesIOName(i));
  
}


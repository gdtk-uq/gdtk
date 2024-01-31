// rsla_test.chpl
// Try out our simple matrix calculations.
//
// Build with:
// $ chpl rsla_test.chpl ../rsla.chpl
//
// Run with:
// $ ./rsla_test
//
// Peter J. 2024-01-31

use Rsla;

proc main() {
  var a, ainv: [D2x2]real;
  var b, x, y: [D2x1] real;
  a[0,0]=0.0; a[0,1]=2.0; a[1,0]=2.0; a[1,1]=2.0;
  x[0]=-0.5; x[1]=1.0;
  MVMult(a, x, b);
  var sflag = MInverse(a, ainv);
  writeln("sflag=", sflag);
  writeln("ainv=", ainv);
  MVMult(ainv, b, y);
  writeln("y=", y);
  //
  var a3, a3inv: [D3x3]real;
  var b3, x3, y3: [D3x1] real;
  a3[0,0]=0.0; a3[0,1]=2.0; a3[0,2]=0.0;
  a3[1,0]=2.0; a3[1,1]=2.0; a3[1,2]=3.0;
  a3[2,0]=4.0; a3[2,1]=-3.0; a3[2,2]=0.0;
  x3[0]=-0.5; x3[1]=1.0; x3[2]=1.0/3;
  MVMult(a3, x3, b3);
  sflag = MInverse(a3, a3inv);
  writeln("sflag=", sflag);
  writeln("a3inv=", a3inv);
  MVMult(a3inv, b3, y3);
  writeln("y3=", y3);
  //
  var a4, a4inv: [{0..#4,0..#4}]real;
  var b4, x4, y4: [{0..#4}] real;
  a4[0,0]=0.0; a4[0,1]=2.0; a4[0,2]=0.0; a4[0,3]=1.0;
  a4[1,0]=2.0; a4[1,1]=2.0; a4[1,2]=3.0; a4[1,3]=2.0;
  a4[2,0]=4.0; a4[2,1]=-3.0; a4[2,2]=0.0; a4[2,3]=1.0;
  a4[3,0]=6.0; a4[3,1]=1.0; a4[3,2]=-6.0; a4[3,3]=-5.0;
  x4[0]=-0.5; x4[1]=1.0; x4[2]=1.0/3; x4[3]=-2.0;
  MVMult(a4, x4, b4);
  writeln("b4=", b4);
  sflag = MInverse(a4, a4inv);
  writeln("sflag=", sflag);
  writeln("a4inv=", a4inv);
  MVMult(a4inv, b4, y4);
  writeln("y4=", y4);
}

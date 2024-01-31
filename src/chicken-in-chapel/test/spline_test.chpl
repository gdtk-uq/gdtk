// spline_test.chpl
// Try out our spline class.
//
// Build with:
// $ chpl spline_test.chpl ../spline.chpl ../rsla.chpl
//
// Run with:
// $ ./spline_test
//
// Peter J. 2024-02-01

use Rsla;
use Spline;

proc runge(x: real): real { return 1.0/(1.0+25.0*x*x); }

proc main() {
  writeln("Test of the Spline class.");
  const ns = 5; const Dsample = {0..#5};
  const x0 = -1.0; const x1 = 1.0;
  const dx = (x1-x0)/(ns-1);
  var x_sample, y_sample: [Dsample]real;
  foreach i in Dsample {
    const xx = x0 + dx*i;
    x_sample[i] = xx; y_sample[i] = runge(xx);
  }
  var s: CubicSpline; s.set(x_sample, y_sample);
  writeln("s=", s);
}


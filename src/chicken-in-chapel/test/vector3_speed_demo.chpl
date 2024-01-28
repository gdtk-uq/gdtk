// vector3_speed_demo.chpl
// A pretend expensive calculation using records
// as we might like to use them in our CFD code.
//
// Build with:
// $ chpl --fast vector3_speed_demo.chpl ../vector3.chpl
//
// Run example on Lenovo T490:
// $ ./vector3_speed_demo
// Begin...
// pi_estimate: 3.14185
// error: 0.000258546
// elapsed_time: 14.043 ms
// Done.
//
// Peter J. 2024-01-22

use Geom;
use Random;
use Math;
use Time;

proc main() {
  writeln("Begin...");

  const n : int = 10_000_000;
  const D = {1..n};
  var aa : [D]Vector3;
  var rs : randomStream(real);
  for v in aa {
    v.x = rs.getNext(); v.y = rs.getNext(); v.z = 0.0;
  }

  var inside_count : [D]int;
  var sw : stopwatch;
  sw.start();
  forall (v, i) in zip(aa, inside_count) {
    if dot(v,v) < 1.0 {
      i = 1;
    }
  }
  sw.stop();
  var sum = + reduce inside_count;
  var elapsed_time = sw.elapsed();
  var pi_estimate = (4.0 * sum) / n;
  writeln("pi_estimate: ", pi_estimate);
  writeln("error: ", pi_estimate-Math.pi);
  writeln("elapsed_time: ", elapsed_time*1000, " ms");

  writeln("Done.");
}


// vector3_speed_gpu.chpl
// A pretend expensive calculation using records
// as we might like to use them in our CFD code.
//
// Build with:
// $ chpl --fast vector3_speed_gpu.chpl vector3.chpl
//
// Run example on Dell Optiplex 990 with NVIDIA GTX 1650:
// peterj@helmholtz ~/work/play/chapel/geom $ ./vector3_speed_gpu 
// Begin...
// Fill arrays.
// Start computation.
// Tally points inside.
// pi_estimate: 3.14138
// error: -0.000208254
// elapsed_time: 3.005 ms
// number of kernel launches=3
// Done.
//
// Peter J. 2024-01-22
//          2024-01-24 check kernel launches

use Geom;
use Random;
use Math;
use Time;
use GpuDiagnostics;

proc startCountingKernelLaunches() {
  resetGpuDiagnostics();
  startGpuDiagnostics();
}

proc numKernelLaunches() {
  stopGpuDiagnostics();
  var result = + reduce getGpuDiagnostics().kernel_launch;
  startCountingKernelLaunches();
  return result;
}

proc main() {
  writeln("Begin...");
  startCountingKernelLaunches();

  on here.gpus[0] {
    const n : int = 10_000_000;
    const D = {1..n};
    var aa : [D]Vector3;
    var rs : randomStream(real);
    writeln("Fill arrays.");
    for v in aa {
      v.x = rs.getNext(); v.y = rs.getNext(); v.z = 0.0;
    }

    writeln("Start computation.");
    var inside_count : [D]int;
    var sw : stopwatch;
    sw.start();
    forall (v, i) in zip(aa, inside_count) {
      if dot(v,v) < 1.0 {
        i = 1;
      }
    }
    sw.stop();
    writeln("Tally points inside.");
    var sum = + reduce inside_count;
    var elapsed_time = sw.elapsed();
    var pi_estimate = (4.0 * sum) / n;
    writeln("pi_estimate: ", pi_estimate);
    writeln("error: ", pi_estimate-Math.pi);
    writeln("elapsed_time: ", elapsed_time*1000, " ms");
  }

  writeln("number of kernel launches=", numKernelLaunches());
  writeln("Done.");
}


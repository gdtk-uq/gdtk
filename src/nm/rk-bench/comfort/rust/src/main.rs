// rk-bench/main.rs
// Try out the Runge-Kutta ODE stepper as a benchmark program.
//
// Run with the command:
// $ cargo run
//
// Author: Peter J.
// Version: 2014-Jun-15, Dlang version adapted from the mech2700 class example
//          2014-Jul-09, preallocate work arrays and pass them in.
//          2018-May-26, work with double or complex numbers
//          2018-May-30, accept the type of the dependent variables as a parameter
//          2022-May-20, Build as a single-source-file program.
//          2022-May-23, Go version
//          2022-May-24, Crystal version
//          2025-Aug-20, NArray version moves away from preallocated arrays
//                       for neatness but at the cost of calculation speed.
//          2025-Aug-22, Rust ndarray version adapted from Crystal version.

use ndarray::{array, Array1};

// Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
//
// Params:
//    f is a callable function that returns the derivative of y wrt t
//       The signature of this function is dydt = f(t, y) where
//       t is a float value, y is an array of number values,
//       dydt is the array with the computed derivatives.
//    t0: is the starting value of the independent variable
//    h: the requested step size
//    y0: an array of starting values for the dependent variables
//        It is assumed that the y-elements are indexed 0 .. n-1
//        where n = y0.size
//
// Returns:
//    t0+h: the final value of the dependent variable
//    y1: an array of final values of the dependent variables
//    err: estimates of the errors in the values of y1
//
fn rkf45_step(
    f : fn(f64, Array1<f64>) -> Array1<f64>,
    t0 : f64,
    h : f64,
    y0 : Array1<f64>
) -> (f64, Array1<f64>, Array1<f64>) {
    // Build up the sample point information as per the text book descriptions.
    let k1 = f(t0, y0.clone());
    let k2 = f(t0 + h/4.0, y0.clone() + 0.25*h*k1.clone());
    let k3 = f(t0 + 3.0*h/8.0, y0.clone() + 3.0*h*k1.clone()/32.0 + 9.0*h*k2.clone()/32.0);
    let k4 = f(t0 + 12.0*h/13.0, y0.clone() + 1932.0*h*k1.clone()/2197.0 - 7200.0*h*k2.clone()/2197.0 +
               7296.0*h*k3.clone()/2197.0);
    let k5 = f(t0 + h, y0.clone() + 439.0*h*k1.clone()/216.0 - 8.0*h*k2.clone() + 3680.0*h*k3.clone()/513.0 -
               845.0*h*k4.clone()/4104.0);
    let k6 = f(t0 + h/2.0, y0.clone() - 8.0*h*k1.clone()/27.0 + 2.0*h*k2.clone() - 3544.0*h*k3.clone()/2565.0 +
               1859.0*h*k4.clone()/4104.0 - 11.0*h*k5.clone()/40.0);
    // Now, do the integration as a weighting of the sampled data.
    let y1 = y0 + 16.0*h*k1.clone()/135.0 + 6656.0*h*k3.clone()/12825.0 +
        28561.0*h*k4.clone()/56430.0 - 9.0*h*k5.clone()/50.0 + 2.0*h*k6.clone()/55.0;
    let err = (h*k1/360.0 - 128.0*h*k3/4275.0 - 2197.0*h*k4/75240.0 +
               h*k5/50.0 + 2.0*h*k6/55.0).mapv(f64::abs);
    (t0+h, y1, err)
}

// Third-order system with a simple analytic solution.
// Adapted from section 11.3 in Cheney and Kincaid, 6th ed.
// Except for the zero-based indexing, the notation is
// chosen to match that in the text.

use libm::exp;

fn solution1(t: f64) -> [f64; 3] {
    let x = exp(-3.0*t)/6.0*(6.0-50.0*exp(t)+10.0*exp(2.0*t)+34.0*exp(3.0*t));
    let y = exp(-3.0*t)/6.0*(12.0-125.0*exp(t)+40.0*exp(2.0*t)+73.0*exp(3.0*t));
    let z = exp(-3.0*t)/6.0*(14.0-200.0*exp(t)+70.0*exp(2.0*t)+116.0*exp(3.0*t));
    [x, y, z]
}

use std::time::Instant;

fn main() {
    println!("Begin demonstration of ODE stepper (Rust)...");

    fn test_system1 (_t : f64, x : Array1<f64>) -> Array1<f64> {
        let dxdt0 =  -8.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 12.0;
        let dxdt1 = -17.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 29.0;
        let dxdt2 = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0;
        array![dxdt0, dxdt1, dxdt2]
    }

    let mut t0 = 0.0;
    let nstep = 10000;
    let h = 1.0/(nstep as f64);
    let mut x0 = array![0.0, 0.0, 0.0];
    let start = Instant::now();
    for _ in 0..nstep {
        let (t1, x1, _err) = rkf45_step(test_system1, t0, h, x0);
        t0 = t1; x0 = x1;
    }
    let elapsed = start.elapsed();
    println!("  elapsed_time= {:?}", elapsed);
    println!("  x           = {}", x0);
    println!("  exact       = {:?}", solution1(t0));
    println!("Done.");
}

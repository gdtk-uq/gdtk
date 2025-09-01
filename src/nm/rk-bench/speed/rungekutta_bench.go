/** rungekutta_demo.go
 *
 * Try out the Runge-Kutta ODE stepper as a benchmark program.
 *
 * Author: Peter J.
 * Version: 2014-Jun-15, adapted from the mech2700 class example
 *          2014-Jul-09, preallocate work arrays and pass them in.
 *          2018-May-26, work with double or complex numbers
 *          2018-May-30, accept the type of the dependent variables as a parameter
 *          2022-May-20, Build as a single-source-file program.
 *          2022-May-23, Go version
 */

package main

import (
    "fmt"
    "math"
    "time"
)

/**
 * Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
 *
 * Params:
 *     f is a callable function that returns the derivative of y wrt t
 *        The signature of this function is f(t, y, dydt) where
 *        t is a float value, y is an array of number values
 *        and dydt is the array to hold the computed derivatives.
 *     t0: is the starting value of the independent variable
 *     h: the requested step size
 *     y0: an array of starting values for the dependent variables
 *         It is assumed that the y-elements are indexed 0 .. n-1
 *         where n = y0.length
 *     y1: an array of final values of the dependent variables
 *     err: estimates of the errors in the values of y1
 *
 * Returns:
 *     the final value of the dependent variable
 */

func rkf45_step(
    f func(float64, []float64, []float64),
    t0 float64, h float64,
    y0 []float64, y1 []float64, err []float64,
    work_arrays [7][]float64) float64 {
    n := len(y0)
    // Assuming a system of equations, we need arrays for the intermediate data.
    k1 := work_arrays[1]
    k2 := work_arrays[2]
    k3 := work_arrays[3]
    k4 := work_arrays[4]
    k5 := work_arrays[5]
    k6 := work_arrays[6]
    ytmp := work_arrays[0]
    // Build up the sample point information as per the text book descriptions.
    // We assign the result of intermediate array expressions to ytmp
    // because that's needed for D.
    f(t0, y0, k1)
    for j := 0; j < n; j++ { ytmp[j] = y0[j] + 0.25*h*k1[j] }
    f(t0 + h/4.0, ytmp, k2)
    for j := 0; j < n; j++ { ytmp[j] = y0[j] + 3.0*h*k1[j]/32.0 + 9.0*h*k2[j]/32.0 }
    f(t0 + 3.0*h/8.0, ytmp, k3)
    for j := 0; j < n; j++ { ytmp[j] = y0[j] + 1932.0*h*k1[j]/2197.0 - 7200.0*h*k2[j]/2197.0 + 7296.0*h*k3[j]/2197.0 }
    f(t0 + 12.0*h/13.0, ytmp, k4)
    for j := 0; j < n; j++ { ytmp[j] = y0[j] + 439.0*h*k1[j]/216.0 - 8.0*h*k2[j] + 3680.0*h*k3[j]/513.0 - 845.0*h*k4[j]/4104.0 }
    f(t0 + h, ytmp, k5)
    for j := 0; j < n; j++ {
        ytmp[j] = y0[j] - 8.0*h*k1[j]/27.0 + 2.0*h*k2[j] -
            3544.0*h*k3[j]/2565.0 + 1859.0*h*k4[j]/4104.0 - 11.0*h*k5[j]/40.0
    }
    f(t0 + h/2.0, ytmp, k6)
    // Now, do the integration as a weighting of the sampled data.
    for j := 0; j < n; j++ {
        y1[j] = y0[j] + 16.0*h*k1[j]/135.0 + 6656.0*h*k3[j]/12825.0 +
            28561.0*h*k4[j]/56430.0 - 9.0*h*k5[j]/50.0 + 2.0*h*k6[j]/55.0
        err[j] = h*k1[j]/360.0 - 128.0*h*k3[j]/4275.0 - 2197.0*h*k4[j]/75240.0 + h*k5[j]/50.0 + 2.0*h*k6[j]/55.0
        err[j] = math.Abs(err[j])
    }
    return t0 + h;
} // end rkf45_step()


/** Test system 1
 * Third-order system with a simple analytic solution.
 * Adapted from section 11.3 in Cheney and Kincaid, 6th ed.
 * Except for the zero-based indexing, the notation is
 * chosen to match that in the text.
 */
func testSystem1(t float64, x []float64, dxdt []float64) {
    dxdt[0] =  -8.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 12.0
    dxdt[1] = -17.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 29.0
    dxdt[2] = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0
}

func solution1(t float64) []float64 {
    x := math.Exp(-3.0*t)/6.0*(6.0-50.0*math.Exp(t)+10.0*math.Exp(2.0*t)+34.0*math.Exp(3.0*t))
    y := math.Exp(-3.0*t)/6.0*(12.0-125.0*math.Exp(t)+40.0*math.Exp(2.0*t)+73.0*math.Exp(3.0*t));
    z := math.Exp(-3.0*t)/6.0*(14.0-200.0*math.Exp(t)+70.0*math.Exp(2.0*t)+116.0*math.Exp(3.0*t));
    return []float64{x, y, z}
}

func main() {
    fmt.Println("Begin demonstration of ODE stepper (Golang)...")

    x1 := make([]float64, 3)
    err := make([]float64, 3)
    // errsum := make([]float64, 3)
    var work [7][]float64
    for i := 0; i < 7; i++ { work[i] = make([]float64, 3) }

    t0 := 0.0
    t1 := 0.0
    Nstep := 10000
    h := 1.0/float64(Nstep)
    x0 := []float64{0.0, 0.0, 0.0}
    start_time := time.Now()
    for i := 0; i < Nstep; i++ {
        t1 = rkf45_step(testSystem1, t0, h, x0, x1, err, work)
        for j := 0; j < 3; j++ { x0[j] = x1[j] }
        t0 = t1
    }
    elapsed_time := time.Now().Sub(start_time).Microseconds()
    fmt.Println("  elapsed_time=", elapsed_time, " us")
    fmt.Println("  x1 = ", x1)
    fmt.Println("  exact = ", solution1(t1))

    fmt.Println("Done.")
}

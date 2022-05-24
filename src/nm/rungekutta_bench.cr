# rungekutta_demo.cr
# Try out the Runge-Kutta ODE stepper as a benchmark program.
#
# Run with the command:
# $ crystal run rungekutta_bench.cr
#
# Author: Peter J.
# Version: 2014-Jun-15, adapted from the mech2700 class example
#          2014-Jul-09, preallocate work arrays and pass them in.
#          2018-May-26, work with double or complex numbers
#          2018-May-30, accept the type of the dependent variables as a parameter
#          2022-May-20, Build as a single-source-file program.
#          2022-May-23, Go version
#          2022-May-24, Crystal version

require "math"

# Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
#
# Params:
#     f is a callable function that returns the derivative of y wrt t
#        The signature of this function is f(t, y, dydt) where
#        t is a float value, y is an array of number values,
#        dydt is the array with the computed derivatives.
#     t0: is the starting value of the independent variable
#     h: the requested step size
#     y0: an array of starting values for the dependent variables
#         It is assumed that the y-elements are indexed 0 .. n-1
#         where n = y0.length
#     y1: an array of final values of the dependent variables
#     err: estimates of the errors in the values of y1
#
# Returns:
#     the final value of the dependent variable
#
def rkf45_step(
     f : Proc(Float64, Array(Float64), Array(Float64), Nil),
     t0 : Float64, h : Float64,
     y0 : Array(Float64), y1 : Array(Float64), err : Array(Float64),
     work_arrays : Array(Array(Float64))
   )
  n = y0.size
  # Assuming a system of equations, we need arrays for the intermediate data.
  ytmp, k1, k2, k3, k4, k5, k6 = work_arrays
  # Build up the sample point information as per the text book descriptions.
  # We assign the result of intermediate array expressions to ytmp
  # because that's needed for D.
  f.call(t0, y0, k1)
  (0...n).each do |j| ytmp[j] = y0[j] + 0.25*h*k1[j] end
  f.call(t0 + h/4.0, ytmp, k2)
  (0...n).each do |j| ytmp[j] = y0[j] + 3.0*h*k1[j]/32.0 + 9.0*h*k2[j]/32.0 end
  f.call(t0 + 3.0*h/8.0, ytmp, k3)
  (0...n).each do |j|
    ytmp[j] = y0[j] + 1932.0*h*k1[j]/2197.0 - 7200.0*h*k2[j]/2197.0 + 7296.0*h*k3[j]/2197.0
  end
  f.call(t0 + 12.0*h/13.0, ytmp, k4)
  (0...n).each do |j|
    ytmp[j] = y0[j] + 439.0*h*k1[j]/216.0 - 8.0*h*k2[j] + 3680.0*h*k3[j]/513.0 - 845.0*h*k4[j]/4104.0
  end
  f.call(t0 + h, ytmp, k5)
  (0...n).each do |j|
    ytmp[j] = y0[j] - 8.0*h*k1[j]/27.0 + 2.0*h*k2[j] -
              3544.0*h*k3[j]/2565.0 + 1859.0*h*k4[j]/4104.0 - 11.0*h*k5[j]/40.0
  end
  f.call(t0 + h/2.0, ytmp, k6)
  # Now, do the integration as a weighting of the sampled data.
  (0...n).each do |j|
    y1[j] = y0[j] + 16.0*h*k1[j]/135.0 + 6656.0*h*k3[j]/12825.0 +
            28561.0*h*k4[j]/56430.0 - 9.0*h*k5[j]/50.0 + 2.0*h*k6[j]/55.0
    err[j] = h*k1[j]/360.0 - 128.0*h*k3[j]/4275.0 - 2197.0*h*k4[j]/75240.0 +
             h*k5[j]/50.0 + 2.0*h*k6[j]/55.0
    err[j] = (err[j]).abs
  end
  t0 + h
end # rkf45_step()


# Third-order system with a simple analytic solution.
# Adapted from section 11.3 in Cheney and Kincaid, 6th ed.
# Except for the zero-based indexing, the notation is
# chosen to match that in the text.

def solution1(t : Float64)
  x = Math.exp(-3.0*t)/6.0*(6.0-50.0*Math.exp(t)+10.0*Math.exp(2.0*t)+34.0*Math.exp(3.0*t))
  y = Math.exp(-3.0*t)/6.0*(12.0-125.0*Math.exp(t)+40.0*Math.exp(2.0*t)+73.0*Math.exp(3.0*t));
  z = Math.exp(-3.0*t)/6.0*(14.0-200.0*Math.exp(t)+70.0*Math.exp(2.0*t)+116.0*Math.exp(3.0*t));
  [x, y, z]
end

def main()
  puts "Begin demonstration of ODE stepper (Crystal)..."

  testSystem1 = -> (t : Float64, x : Array(Float64), dxdt : Array(Float64)) do
    dxdt[0] =  -8.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 12.0
    dxdt[1] = -17.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 29.0
    dxdt[2] = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0
    Nil
  end

  x1 = Array(Float64).new(3, 0.0)
  err = Array(Float64).new(3, 0.0)
  # errsum = Array(Float64, 3, 0.0)
  work = Array(Array(Float64)).new
  (0...7).each do work << Array(Float64).new(3, 0.0) end

  t0 = 0.0
  t1 = 0.0
  nstep = 10000
  h = 1.0/nstep
  x0 = [0.0, 0.0, 0.0] of Float64
  elapsed_time = Time.measure do
    (0...nstep).each do
      t1 = rkf45_step(testSystem1, t0, h, x0, x1, err, work)
      (0...3).each do |j| x0[j] = x1[j] end
      t0 = t1
    end
  end
  puts "  elapsed_time= #{elapsed_time}"
  puts "  x1 = #{x1}"
  puts "  exact = #{solution1(t1)}"

  puts "Done."
end

main

# rungekutta_bench.cr
# Try out the Runge-Kutta ODE stepper as a benchmark program.
#
# Run with the command:
# $ crystal run rungekutta_bench.cr
#
# Author: Peter J.
# Version: 2014-Jun-15, Dlang version adapted from the mech2700 class example
#          2014-Jul-09, preallocate work arrays and pass them in.
#          2018-May-26, work with double or complex numbers
#          2018-May-30, accept the type of the dependent variables as a parameter
#          2022-May-20, Build as a single-source-file program.
#          2022-May-23, Go version
#          2022-May-24, Crystal version
#          2025-Aug-20, NArray version moves away from preallocated arrays
#                       for neatness but at the cost of calculation speed.

require "math"
require "./numarray.cr"

# Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
#
# Params:
#     f is a callable function that returns the derivative of y wrt t
#        The signature of this function is dydt = f(t, y) where
#        t is a float value, y is an array of number values,
#        dydt is the array with the computed derivatives.
#     t0: is the starting value of the independent variable
#     h: the requested step size
#     y0: an array of starting values for the dependent variables
#         It is assumed that the y-elements are indexed 0 .. n-1
#         where n = y0.size
#
# Returns:
#     t0+h: the final value of the dependent variable
#     y1: an array of final values of the dependent variables
#     err: estimates of the errors in the values of y1
#
def rkf45_step(
     f : Proc(Float64, NArray(Float64), NArray(Float64)),
     t0 : Float64,
     h : Float64,
     y0 : NArray(Float64)
  )
  # Build up the sample point information as per the text book descriptions.
  k1 = f.call(t0, y0.clone)
  k2 = f.call(t0 + h/4.0, y0 + 0.25*h*k1)
  k3 = f.call(t0 + 3.0*h/8.0, y0 + 3.0*h*k1/32.0 + 9.0*h*k2/32.0)
  k4 = f.call(t0 + 12.0*h/13.0, y0 + 1932.0*h*k1/2197.0 - 7200.0*h*k2/2197.0 +
                                7296.0*h*k3/2197.0)
  k5 = f.call(t0 + h, y0 + 439.0*h*k1/216.0 - 8.0*h*k2 + 3680.0*h*k3/513.0 -
                      845.0*h*k4/4104.0)
  k6 = f.call(t0 + h/2.0, y0 - 8.0*h*k1/27.0 + 2.0*h*k2 - 3544.0*h*k3/2565.0 +
                          1859.0*h*k4/4104.0 - 11.0*h*k5/40.0)
  # Now, do the integration as a weighting of the sampled data.
  y1 = y0 + 16.0*h*k1/135.0 + 6656.0*h*k3/12825.0 +
       28561.0*h*k4/56430.0 - 9.0*h*k5/50.0 + 2.0*h*k6/55.0
  err = (h*k1/360.0 - 128.0*h*k3/4275.0 - 2197.0*h*k4/75240.0 +
         h*k5/50.0 + 2.0*h*k6/55.0).abs!
  return t0+h, y1, err
end


# Third-order system with a simple analytic solution.
# Adapted from section 11.3 in Cheney and Kincaid, 6th ed.
# Except for the zero-based indexing, the notation is
# chosen to match that in the text.

include Math

def solution1(t : Float64)
  x = exp(-3.0*t)/6.0*(6.0-50.0*exp(t)+10.0*exp(2.0*t)+34.0*exp(3.0*t))
  y = exp(-3.0*t)/6.0*(12.0-125.0*exp(t)+40.0*exp(2.0*t)+73.0*exp(3.0*t))
  z = exp(-3.0*t)/6.0*(14.0-200.0*exp(t)+70.0*exp(2.0*t)+116.0*exp(3.0*t))
  [x, y, z]
end

def main()
  puts "Begin demonstration of ODE stepper (Crystal)..."

  testSystem1 = -> (t : Float64, x : NArray(Float64)) do
    dxdt0 =  -8.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 12.0
    dxdt1 = -17.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 29.0
    dxdt2 = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0
    NArray.new([dxdt0, dxdt1, dxdt2])
  end

  t0 = 0.0
  t1 = 0.0
  nstep = 10000
  h = 1.0/nstep
  x0 = NArray.new([0.0, 0.0, 0.0])
  elapsed_time = Time.measure do
    nstep.times do
      t1, x1, err = rkf45_step(testSystem1, t0, h, x0)
      t0 = t1; x0 = x1
    end
  end
  puts "  elapsed_time= #{elapsed_time}"
  puts "  x     = #{x0}"
  puts "  exact = #{solution1(t0)}"
  puts "Done."
end

main

# zero_solvers.rb
=begin
Solve nonlinear functions of a single variable.

Author: Rowan J Gollan and Peter J.

Versions:
  06-Dec-2004
  08-May-2011: Dan's bisection_method added by PJ
  16-Apr-2012: PJ, make more efficient by not evaluating f redundantly
    Also, make the code more compact (so that it fits in the editor window).
  29-Dec-2019: PJ, Python3 port.  Make better use of exceptions.
  10-Feb-2020: PJ, Ruby port
=end

class ZeroSolversError < ::StandardError
end

module ZeroSolvers

  def secant(f, x0, x1, tol=1.0e-11, limits=[], max_iterations=1000, tf=false)
    # The iterative secant method for zero-finding in one-dimension.
    #
    # f: user-defined function f(x)
    # x0: first guess
    # x1: second guess, presumably close to x0
    # tol: stopping tolerance for f(x)=0
    # max_iterations: to stop the iterations running forever, just in case...
    # tf: boolean flag to turn on printing of intermediate states
    #
    # Returns: x such that f(x)=0
    #
    # We're going to arrange x0 as the oldest (furtherest) point
    # and x1 and the closer-to-the-solution point.
    # x2, when we compute it, will be the newest sample point.
    f0 = f.call(x0); f1 = f.call(x1)
    if f0.abs < f1.abs then
      x0, f0, x1, f1 = x1, f1, x0, f0
    end
    max_iterations.times do |i|
      begin
        x2 = x1 - f1 * (x0 - x1) / (f0 - f1)
      rescue ZeroDivisionError
        raise ZeroSolversError, 'Cannot proceed with zero slope.'
      end
      if limits != [] then
        x2 = [limits[0], x2].max
        x2 = [limits[1], x2].min
      end
      f2 = f.call(x2)
      if tf then
        puts "  %d \t  %f \t %f \t %f \t %e" % [i+1, x0, x1, x2, f2]
      end
      x0, f0, x1, f1 = x1, f1, x2, f2
      if f2.abs < tol then
        return x2
      end
    end
    raise ZeroSolversError, 'Did not converge after %d iterations.' % max_iterations
  end

  def bisection(f, bx, ux, tol=1.0e-6)
    # The iterative bisection method for zero-finding in one-dimension.
    #
    # f: user-defined function f(x)
    # bx: bottom-limit of bracket
    # ux: upper-limit of bracket
    # tol: stopping tolerance on bracket size
    #
    # Returns: x such that f(x)=0
    while (ux-bx).abs > tol do
      midpoint = 0.5*(bx+ux)
      if f.call(bx) * f.call(midpoint) > 0 then
        bx = midpoint
      else
        ux = midpoint
      end
    end
    return 0.5*(bx+ux)
  end

end # module ZeroSolvers

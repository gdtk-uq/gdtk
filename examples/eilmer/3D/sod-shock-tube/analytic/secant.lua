-- secant.lua
-- PAJ 2018-04-27
-- Adapted from the Python version in zero_solvers.py
--

--[=[
The iterative secant method for zero-finding in one-dimension.

   f: user-defined function f(x)
   x0: first guess
   x1: second guess, presumably close to x0
   tol: stopping tolerance for f(x)=0
   limits: array of the min and max values allowed for x
   max_iterations: to stop the iterations running forever, just in case...

   returns: x such that f(x)=0 or the string 'FAIL'
--]=]
function secant(f, x0, x1, tol, limits, max_iterations)
   tol = tol or 1.0e-11
   limits = limits or {}
   max_iterations = max_iterations or 1000
   -- We're going to arrange x0 as the oldest (furtherest) point
   -- and x1 and the closer-to-the-solution point.
   -- x2, when we compute it, will be the newest sample point.
   local f0 = f(x0); local f1 = f(x1)
   if math.abs(f0) < math.abs(f1) then
      x0, f0, x1, f1 = x1, f1, x0, f0
   end
   for i=1,max_iterations do
      x2 = x1 - f1 * (x0 - x1) / (f0 - f1)
      if #limits == 2 then
         x2 = math.min(limits[2], math.max(limits[1], x2))
      end
      f2 = f(x2)
      if TEST_MODULE then
         print(string.format('  %d \t  %f \t %f \t %f \t %e', i+1,x0,x1,x2,f2))
      end
      x0, f0, x1, f1 = x1, f1, x2, f2
      if math.abs(f2) < tol then return x2 end
   end
   if TEST_MODULE then
      print('The secant method did not converge after ', i+1, ' iterations')
   end
   return 'FAIL'
end

--[=[
print("Begin secant_method test...")
TEST_MODULE = true
function ftest(x) return (math.sin(x)) end
print("one solution at x=", secant(ftest, 6.0, 6.5))
print("expected x=", (2.0*math.pi))
print("Done.")
--]=]

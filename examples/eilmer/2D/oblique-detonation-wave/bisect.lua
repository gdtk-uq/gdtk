-- bisect.lua
-- A little brute-force bracketing function solver.
-- Peter J, 2017-01-06

function solve(f, a, b, tol)
   tol = tol or 1.0e-9
   local fa = f(a)
   if math.abs(fa) < tol then return a end
   local fb = f(b)
   if math.abs(fb) < tol then return b end
   assert(fa*fb < 0.0, "crap bracket")
   for i=1, 50 do
      local c = 0.5*(a+b)
      local fc = f(c)
      if math.abs(fc) < tol then return c end
      if fa*fc < 0.0 then
	 b = c; fb = fc
      else
	 a = c; fa = fc
      end
   end
   -- surely 50 iterations has cut the interval enough
   return c
end

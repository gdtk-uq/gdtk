-- secant_test.lua

local secant = require "secant"

print("Begin secant_method test...")
TEST_MODULE = true
function ftest(x) return math.sin(x) end -- normal behaviour
function ftest2(x) return 1.0e-12 end -- early success
function ftest3(x) return x^2 + 1.0 end -- no zero to be found
xsol, err = secant.solve(ftest, 6.0, 6.5)
if err then
   print("error=", err)
else
   print("one solution at x=", xsol)
end
print("expected x=", (2.0*math.pi))
print("Done.")

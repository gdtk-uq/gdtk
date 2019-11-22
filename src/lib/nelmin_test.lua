--[[
nelmin_test.lua
Serves as a test and demonstration.

peterj@helmholtz ~/dgd/src/lib $ gas-calc nelmin_test.lua
test 1: simple quadratic with zero at (1, 1, ...)
actual f(1.000, 1.000, 1.000, 1.000)=1.642e-10
convergence-flag = 	false
number-of-fn-evaluations = 	365
number-of-restarts = 	3
test 2: Example 3.3 from Olsson and Nelson
expect f(0.811, -0.585)=-67.1
actual f(0.811, -0.585)=-67.1
convergence-flag = 	true
number-of-fn-evaluations = 	187
number-of-restarts = 	0
test 3: Example 3.5 in Olsson and Nelson, nonlinear least-squares
expect f(1.801, -1.842, -0.463, -1.205)=0.0009
actual f(1.801, -1.842, -0.463, -1.205)=0.0009
convergence-flag = 	true
number-of-fn-evaluations = 	618
number-of-restarts = 	0
]]

local nelmin = require "nelmin"

print("test 1: simple quadratic with zero at (1, 1, ...)")
local testfun1 = function(x)
   local sum = 0.0
   for i,elem in ipairs(x) do
      sum = sum + (elem-1)*(elem-1)
   end
   return sum
end
local x = {0, 0, 0, 0}
local x_best, f_best, converged, nfe, nrestarts = nelmin.minimise(testfun1, x)
print(string.format("actual f(%.3f, %.3f, %.3f, %.3f)=%.3e",
                    x_best[1], x_best[2], x_best[3], x_best[4], f_best))
print("convergence-flag = ", converged)
print("number-of-fn-evaluations = ", nfe)
print("number-of-restarts = ", nrestarts)

print("test 2: Example 3.3 from Olsson and Nelson")
print("expect f(0.811, -0.585)=-67.1")
local testfun2 = function(x)
   if ((x[1]*x[1] + x[2]*x[2]) > 1.0) then
      return 1.0e38
   else
      local yp = 53.69+7.26*x[1]-10.33*x[2]+7.22*x[1]*x[1]+6.43*x[2]*x[2]+11.36*x[1]*x[2]
      local ys = 82.17-1.01*x[1]-8.61*x[2]+1.40*x[1]*x[1]-8.76*x[2]*x[2]-7.20*x[1]*x[2]
      return -yp + math.abs(ys-87.8)
   end
end
x = {0, 0}
local dx = {0.5, 0.5}
x_best, f_best, converged, nfe, nrestarts = nelmin.minimise(testfun2, x, dx)
print(string.format("actual f(%.3f, %.3f)=%.1f", x_best[1], x_best[2], f_best))
print("convergence-flag = ", converged)
print("number-of-fn-evaluations = ", nfe)
print("number-of-restarts = ", nrestarts)

print("test 3: Example 3.5 in Olsson and Nelson, nonlinear least-squares")
print("expect f(1.801, -1.842, -0.463, -1.205)=0.0009")
local testfun3 = function(z)
   local x = {0.25, 0.50, 1.00, 1.70, 2.00, 4.00}
   local y = {0.25, 0.40, 0.60, 0.58, 0.54, 0.27}
   local a1, a2, alpha1, alpha2 = z[1], z[2], z[3], z[4]
   local sum_residuals = 0.0
   for i=1,6 do
      local t = x[i]
      local eta = a1*math.exp(alpha1*t) + a2*math.exp(alpha2*t)
      local r = y[i] - eta
      sum_residuals = sum_residuals + r*r
   end
   return sum_residuals
end
x = {1, 1, -0.5, -2.5}
local dx = {0.1, 0.1, 0.1, 0.1}
x_best, f_best, converged, nfe, nrestarts = nelmin.minimise(testfun3, x, dx, 1.0e-9, 800)
print(string.format("actual f(%.3f, %.3f, %.3f, %.3f)=%.4f",
                    x_best[1], x_best[2], x_best[3], x_best[4], f_best))
print("convergence-flag = ", converged)
print("number-of-fn-evaluations = ", nfe)
print("number-of-restarts = ", nrestarts)

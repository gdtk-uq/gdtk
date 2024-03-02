-- analytic.lua
-- Analytic solution functions for the oblique detonation wave
-- of
-- JM Powers and TD Aslam (2006)
-- Exact solution for multidimensional compressible reactive flow
-- for verifying numerical algorithms.
-- AIAA Journal Vol. 44 No. 2 pages 337-344
--
-- This script by Peter J. and Rowan G.
-- 2017-01-06
--
-- Parameters from table 1 (set as global variables)
beta = math.pi/4 -- radians
sinB = math.sin(beta)
cosB = math.cos(beta)
T1 = 300.0 -- K
M1 = 3.0 
rho1 = 1.0 -- kg/m^3
R = 287.0 -- J/kg/K
alpha = 1000.0 -- 1/s
gamma = 6.0/5
q = 300000.0 -- J/kg
p1 = rho1*R*T1
a1 = math.sqrt(gamma*R*T1)
u1 = M1*a1
-- print("u1=", u1, "p1=", p1)

function calculate_X(lmbda)
   local MsinBeta2 = (M1*sinB)^2
   local a1 = (1/((gamma+1)*M1*sinB)) * (a1/alpha)
   local a2 = 1+gamma*MsinBeta2
   local a3 = MsinBeta2 - 1
   local a4 = 2*MsinBeta2/((MsinBeta2-1)^2) *
      (gamma^2-1)/gamma * q/(R*T1)
   local oneMinusA4L = 1-a4*lmbda
   local oneMinusA4 = 1-a4
   local t1 = 2*a3*(math.sqrt(oneMinusA4L) - 1)
   local t2 = math.pow(1/(1-lmbda), a2)
   local t3 = 1 - math.sqrt(oneMinusA4L/oneMinusA4)
   local t4 = 1 + math.sqrt(1/oneMinusA4)
   local t5 = 1 + math.sqrt(oneMinusA4L/oneMinusA4)
   local t6 = 1 - math.sqrt(1/oneMinusA4)
   local X = a1*(t1 + math.log(t2*math.pow((t3*t4)/(t5*t6),
				  a3*math.sqrt(oneMinusA4))))
   return X
end

function calculate_rho(lmbda)
   local MsinBeta2 = (M1*sinB)^2
   local t1 = rho1*(gamma+1)*MsinBeta2
   local t2 = 1+gamma*MsinBeta2
   local t3 = t2^2
   local t4 = (gamma+1)*MsinBeta2
   local t5 = (gamma-1)/gamma * 2*lmbda*q/(R*T1)
   local t6 = (gamma-1)*MsinBeta2
   local rho = t1/(t2 - math.sqrt(t3 - t4*(2+t5+t6)))
   return rho
end

function calculate_U(lmbda, rho)
   local U = rho1*u1*sinB/rho
   return U
end

function calculate_T(lmbda, rho)
   local t1 = p1/(rho*R)
   local t2 = (rho1*u1*sinB)^2 / (rho*R)
   local t3 = 1/rho1 - 1/rho
   local T = t1 + t2*t3
   return T
end

function calculate_p(lmbda, rho)
   local t2 = (rho1*u1*sinB)^2
   local t3 = 1/rho1 - 1/rho
   local p = p1 + t2*t3
   return p
end

function calculate_Yw(lmbda)
   local Yw = (u1*cosB/alpha) * math.log(1/(1-lmbda))
   return Yw
end

--[[
local lmbda = 0.1
local rho = calculate_rho(lmbda)
print("lmbda=", lmbda, "X=", calculate_X(lmbda), "rho=", rho)
print("U=", calculate_U(lmbda, rho), "T=", calculate_T(lmbda, rho))
print("p=", calculate_p(lmbda, rho), "Yw=", calculate_Yw(lmbda))
--]]

function transform_xy_to_XY(x, y)
   local X = x*sinB - y*cosB
   local Y = x*cosB + y*sinB
   return X, Y
end

function transform_XY_to_xy(X, Y)
   local x = X*sinB + Y*cosB
   local y = Y*sinB - X*cosB
   return x, y
end

function transform_UV_to_uv(U, V)
   local u = U*sinB + V*cosB
   local v = V*sinB - U*cosB
   return u, v
end

dofile("bisect.lua")

function find_lmbda_from_x(x)
   function errx(lmbda)
      local X = calculate_X(lmbda)
      local Yw = calculate_Yw(lmbda)
      local xg, yg = transform_XY_to_xy(X, Yw)
      return x - xg
   end
   local lmbda = solve(errx, 0.0, 0.999)
   return lmbda
end

--[[
local xp = 0.797749
local lmbda = find_lmbda_from_x(xp)
print("xp=", xp, "lmbda=", lmbda)
--]]

function create_wall_function(xmin, xmax)
   function wall_function(t)
      local lmbda = find_lmbda_from_x(t*(xmax-xmin))
      local X = calculate_X(lmbda)
      local Yw = calculate_Yw(lmbda)
      local x, y = transform_XY_to_xy(X, Yw)
      return {x=x+xmin, y=y, z=0}
   end
   return wall_function
end

--[[
myWallFn = create_wall_function(0.0, 1.0)
for _,t in ipairs{0.0, 0.066116, 0.141020} do
   local point = myWallFn(t)
   print("x=", point.x, "y=", point.y, "z=", point.z)
end
--]]

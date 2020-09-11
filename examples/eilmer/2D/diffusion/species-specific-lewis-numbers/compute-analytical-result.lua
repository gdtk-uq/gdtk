-- Author: Rowan J. Gollan
-- Date: 2017-07-02
--
-- This script is used to compute the time-dependent
-- diffusion of two species at the same temperature
-- and pressure that are initially separated by a
-- membrane.
--
-- The exact solution is presented in the paper by
-- Sanchez et al.
--
-- Reference:
-- Sanchez, Vera, Cinan (2006)
-- Exact Solutions for Transient Mixing of Two Gases of
-- Different Densities.
-- Physics of Fluids, 18(7)
--
-- Process this script using gas-calc:
--
-- > gas-calc compute-analytical-result.lua
--

outfile = 'exact-profile.dat'

gm = GasModel:new{'thermally-perfect-N2-O2-mix.lua'}
Q = GasState:new{gm}

t = 1.0e-6
D = 1.645306207462e-05
p = 100e3
T = 273.2
Q.p = p
Q.T = T
Q.massf = {N2=1, O2=0}
gm:updateThermoFromPT(Q)
rhoA = Q.rho
-- Now switch gas composition
Q.massf = {N2=0, O2=1}
gm:updateThermoFromPT(Q)
rhoB = Q.rho

-- plotting ranges
xMin = -2.0e-5
xMax = 2.0e-5
dx = 0.25e-6

-- Error function implementation found at:
-- http://hewgill.com/picomath/lua/erf.lua.html

function erf(x)
    -- constants
    local a1 =  0.254829592
    local a2 = -0.284496736
    local a3 =  1.421413741
    local a4 = -1.453152027
    local a5 =  1.061405429
    local p  =  0.3275911

    -- Save the sign of x
    local sign = 1
    if x < 0 then
        sign = -1
    end
    local x1 = math.abs(x)

    -- A&S formula 7.1.26
    local t = 1.0/(1.0 + p*x1)
    local y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x1*x1)

    return sign*y
end

function calcRho(r, t)
   arg = r / (2.0*math.sqrt(D*t))
   rho = 0.5 * (1 - erf(arg)) * (rhoA - rhoB) + rhoB
   return rho
end

function calcFA(rho)
   return (rhoA * (rhoB - rho))/(rho * (rhoB - rhoA))
end

function main()
   f = assert(io.open(outfile, 'w'))
   f:write("#   x(m)    rho    f_N2    f_O2\n")
   for x=xMin,xMax,dx do
      rho = calcRho(x, t)
      f_N2 = calcFA(rho)
      f_O2 = 1.0 - f_N2
      f:write(string.format("%.8e  %.8f  %.8f  %.8f\n",
			    x, rho, f_N2, f_O2))
   end
   f:close()
end

main()
   




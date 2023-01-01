-- billig.lua
-- Fred Billig's correlations for hypersonic shock-wave shapes.
--
-- These are a direct implementation of equations 5.36, 5.37 and 5.38
-- from J.D. Anderson's text Hypersonic and High Temperature Gas Dynamics
--
-- PJ 2016-10-16 adapted from Python module
--
-- RJG 2019-03-06 Re-packaged as module.
--
-- NOTE on use as a module.
-- ------------------------
-- This module will not work stand-alone since it is dependent on having
-- idealgasflow functions registered in the namespace. If you are working
-- with the dgd code collection, this module WILL work with the custom-built
-- version of lua provided during install: dgd-lua
--

local billig = {}

function billig.delta_over_R(M, axisymmetric)
   -- Calculates the normalised stand-off distance.
   -- Assume cylindrical nose
   local d_R = 0.386 * math.exp(4.67/(M*M))
   if axisymmetric then
      -- Spherical nose
      d_R = 0.143 * math.exp(3.24/(M*M))
   end
   return d_R
end

function billig.Rc_over_R(M, axisymmetric)
   -- Calculates the normalised radius of curvature of the shock.
   -- Assume cylindrical nose
   local Rc_R = 1.386 * math.exp(1.8/math.pow(M-1, 0.75))
   if axisymmetric then
      -- Spherical nose
      Rc_R = 1.143 * math.exp(0.54/math.pow(M-1, 1.2))
   end
   return Rc_R
end

function billig.x_from_y(y, M, theta, axisymmetric, R_nose)
   -- Determine the x-coordinate of a point on the shock wave.
   --
   -- y: y-coordinate of the point on the shock wave
   -- M: free-stream Mach number
   -- theta: angle of the downstream surface (radians)
   --    with respect to free-stream direction.
   --    For example, a blunted plate or cylinder will have
   --    a surface angle of zero.
   -- axisymmetric: flag
   --    false or nil : cylinder-wedge
   --    true: sphere-cone
   -- R_nose: radius of the forebody (either cylinder or sphere)
   --
   -- It is assumed that, for the ideal gas, gamma=1.4.
   -- That's the only value relevant to the data used for
   -- Billig's correlations.
   --
   R_nose = R_nose or 1.0
   theta = theta or 0.0
   local Rc = R_nose * billig.Rc_over_R(M, axisymmetric)
   local d = R_nose * billig.delta_over_R(M, axisymmetric)
   local beta
   if theta == 0 then
      beta = math.asin(1.0/M)
   else
      if axisymmetric then
         -- Use shock angle on a cone
         beta = idealgasflow.beta_cone2(M, theta)
      else
         -- Use shock angle on a wedge
         beta = idealgasflow.beta_obl(M, theta)
      end
   end
   local tan_beta = math.tan(beta)
   local cot_beta = 1.0/tan_beta
   local x = R_nose + d - Rc * cot_beta^2 * (math.sqrt(1+(y*tan_beta/Rc)^2) - 1)
   return x
end

function billig.y_from_x(x, M, theta, axisymmetric, R_nose)
   -- Determine the y-coordinate of a point on the shock wave.
   --
   -- x: x-coordinate of the point on the shock wave
   -- M: free-stream Mach number
   -- theta: angle of the downstream surface (radians)
   --        with respect to free-stream direction
   -- axisymmetric: flag
   --    false or nil : cylinder-wedge
   --    true: sphere-cone
   -- R_nose: radius of the forebody (either cylinder or sphere)
   --
   R_nose = R_nose or 1.0
   theta = theta or 0.0
   local Rc = R_nose * Rc_over_R(M, axi)
   local d = R_nose * delta_over_R(M, axi)
   local beta
   if theta == 0 then
      beta = math.asin(1.0/theta)
   else
      if axisymmetric then
         -- Use shock angle on a cone
         beta = idealgasflow.beta_cone2(M, theta)
      else
         -- Assume use of shock angle on a wedge
         beta = idealgasflow.beta_obl(M, theta)
      end
   end
   local tan_beta = math.tan(beta)
   local cot_beta = 1.0/tan_beta
   local tmpA = (x - R_nose - d)/(-Rc * cot_beta^2) + 1
   local y = math.sqrt( ((tmpA^2 - 1) * Rc^2) / (tan_beta^2) )
   return y
end

return billig

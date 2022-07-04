-- billig_patch.lua
-- A convenience function for bluff-body simulations.
-- PJ 2019-05-25
--    2022-06-29 introduce quadrant==3 option
--    2022-07-04, Allow different x and y scales

local billig = require 'billig'
local billig_patch = {}

function billig_patch.make_patch(t)
   -- Construct a surface patch for use in a bluff-body simulation
   -- using arguments found in the supplied table.
   -- For the quadrant==2 case the patch is above the x-axis.
   --
   --                      ++d[#d]
   --                    +   |
   --                  +     | North
   --        West    +       |
   --               +        |
   --              +      ---b
   --             +      /
   --            +      /
   --            +     |
   --          d[1]----a     c  ----> x
   --            South
   --
   -- For the quandrant==3 case the patch is below the x-axis.
   --
   --            North
   --          d[1]----a     c  ----> x
   --            +     |
   --            +      \
   --             +      \
   --              +      ---b
   --        West   +        |
   --                +       | South
   --                  +     |
   --                    +   |
   --                      ++d[#d]
   --
   if not type(t) == "table" then
      error("Expected a single table containing named arguments.", 2)
   end
   -- Arguments without default values.
   -- Free-stream Mach number.
   local Minf = t.Minf
   -- Radius of cylinder or sphere.
   local R = t.R
   -- Arguments with default values.
   local quadrant = t.quadrant or 2
   -- Position of centre of body.
   local xc = t.xc or 0.0
   local yc = t.yc or 0.0
   -- Scale for accommodating thermochemical variation
   -- away from ideal low-Temperature air.
   -- Set the scales using the most specific information available,
   -- eventually defaulting to a scale of 1.0 if nothing is specified.
   local x_scale = nil
   local y_scale = nil
   if t.x_scale then x_scale = t.x_scale end
   if t.y_scale then y_scale = t.y_scale end
   if t.scale then
      x_scale = x_scale or t.scale
      y_scale = y_scale or t.scale
   end
   x_scale = x_scale or 1.0
   y_scale = y_scale or 1.0
   -- Assume 2D cylinder (axisymmetric=true for a sphere).
   local axisymmetric = t.axisymmetric or false
   -- Angle of aft-body with respect to freestream direction.
   local theta = t.theta or 0.0
   --
   local a = Vector3:new{x=xc-R, y=yc}
   local b = Vector3:new{x=xc, y=yc+R}
   if quadrant == 3 then b = Vector3:new{x=xc, y=yc-R} end
   local c = Vector3:new{x=xc, y=yc}
   local body
   if quadrant == 2 then
      body = Arc:new{p0=a, p1=b, centre=c}
   else
      body = Arc:new{p0=b, p1=a, centre=c}
   end
   --
   -- In order to have a grid that fits reasonably close the the shock,
   -- use Billig's shock shape correlation to generate
   -- a few sample points along the expected shock position.
   print("Points on Billig's correlation.")
   local xys = {}
   for i,y in ipairs({0.0, 0.2, 0.4, 0.6, 1.0, 1.4, 1.6, 2.0, 2.37}) do
      -- y is normalized for R=1
      local x = billig.x_from_y(y*R, Minf, theta, axisymmetric, R)
      xys[#xys+1] = {x=x, y=y*R}  -- a new coordinate pair
      -- print("x=", x, "y=", y*R)
   end
   -- Scale the Billig distances, depending on the expected behaviour
   -- relative to the gamma=1.4 ideal gas.
   local d = {} -- will use a list to keep the nodes for the shock boundary
   local shock
   local xaxis
   local outlet
   local patch
   if quadrant == 2 then
      for i, xy in ipairs(xys) do
         -- the outer boundary should be a little further than the shock itself
         d[#d+1] = Vector3:new{x=-x_scale*xy.x+xc, y=y_scale*xy.y+yc}
      end
      shock = ArcLengthParameterizedPath:new{underlying_path=Spline:new{points=d}}
      --
      xaxis = Line:new{p0=d[1], p1=a} -- shock to nose of body
      outlet = Line:new{p0=d[#d], p1=b} -- shock to last-point of body
      --
      patch = CoonsPatch:new{south=xaxis, north=outlet, west=shock, east=body}
   else
      -- For quadrant 3, the shock path progresses fro the outer point to the axis.
      for i=1,#xys do
         -- the outer boundary should be a little further than the shock itself
         local xy = xys[#xys-i+1]
         d[#d+1] = Vector3:new{x=-x_scale*xy.x+xc, y=-y_scale*xy.y+yc}
      end
      shock = ArcLengthParameterizedPath:new{underlying_path=Spline:new{points=d}}
      --
      xaxis = Line:new{p0=d[#d], p1=a} -- shock to nose of body
      outlet = Line:new{p0=d[1], p1=b} -- shock to last-point of body
      --
      patch = CoonsPatch:new{south=outlet, north=xaxis, west=shock, east=body}
   end
   return {patch=patch, points={a=a, b=b, c=c, d=d}}
end

return billig_patch

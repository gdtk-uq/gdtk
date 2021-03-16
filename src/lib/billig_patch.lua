-- billig_patch.lua
-- A convenience function for bluff-body simulations.
-- PJ 2019-05-25

module(..., package.seeall)

require 'billig'

function make_patch(t)
   -- Construct a surface patch for use in a bluff-body simulation
   -- using arguments found in the supplied table.
   if not type(t) == "table" then
      error("Expected a single table containing named arguments.", 2)
   end
   -- Arguments without default values.
   -- Free-stream Mach number.
   local Minf = t.Minf
   -- Radius of cylinder or sphere.
   local R = t.R
   -- Arguments with default values.
   -- Position of centre of body.
   local xc = t.xc or 0.0
   local yc = t.yc or 0.0
   -- Scale for accommodating thermochemical variation
   -- away from ideal low-Temperature air.
   local scale = t.scale or 1.0
   -- Assume 2D cylinder (axisymmetric=true for a sphere).
   local axisymmetric = t.axisymmetric or false
   local theta = theta or 0.0
   --
   local a = Vector3:new{x=xc-R, y=yc}
   local b = Vector3:new{x=xc, y=yc+R}
   local c = Vector3:new{x=xc, y=yc}
   local body = Arc:new{p0=a, p1=b, centre=c}
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
   for i, xy in ipairs(xys) do
      -- the outer boundary should be a little further than the shock itself
      d[#d+1] = Vector3:new{x=-scale*xy.x+xc, y=scale*xy.y+yc}
   end
   -- print("front of grid: d[1]=", d[1])
   local shock = ArcLengthParameterizedPath:new{underlying_path=Spline:new{points=d}}
   --
   local xaxis = Line:new{p0=d[1], p1=a} -- first-point of shock to nose of body
   local outlet = Line:new{p0=d[#d], p1=b}  -- top-point of shock to top of body
   --
   local patch = CoonsPatch:new{south=xaxis, north=outlet, west=shock, east=body}
   return {patch=patch, points={a=a, b=b, c=c, d=d}}
end

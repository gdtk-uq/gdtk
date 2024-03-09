print("Mohammadian's convex-ramp experiment.")
-- Peter J., Dan Potter and Rowan Gollan, Eilmer3 flavour, 15-Aug-2013
-- PJ, Rowan and Kyle, Eilmer5 flavour, 2024-03-09

config.dimensions = 2

-- Mohammadian used the inch as his length scale.
m_per_inch = 0.0254
mm = 1.0e-3 -- metres per mm

function ramp(t)
   -- Parametric definition of ramp profile in xy-plane.
   -- Here, we join the initial straight 18-degree ramp to the polynomial.
   local alpha = 18.0*math.pi/180.0 -- angle of initial straight section
   local sin18 = math.sin(alpha)
   local cos18 = math.cos(alpha)
   local tan18 = math.tan(alpha)
   local x_join_inch = 3.0
   local y_join_inch = x_join_inch * tan18
   local L1 = x_join_inch/cos18 -- length of initial straight section
   local L2 = 4.14677 -- length of fairing (computed via maxima)
   local t2 = (L1+L2) * t
   local x_inch, y_inch
   if t2 < L1 then
      x_inch = t2 * cos18
      y_inch = t2 * sin18
   else
      s = (t2 - L1)/L2 * 4.0577
      g = 0.0026 * math.pow(s,4) - 0.0211 * math.pow(s,3)
      x_inch = x_join_inch + s * cos18 - g * sin18
      y_inch = y_join_inch + s * sin18 + g * cos18
   end
   return {x=x_inch*m_per_inch, y=y_inch*m_per_inch, z=0.0}
end

-- Leading edge of ramp
local a0 = Vector3:new(ramp(0.0)); local a1 = a0 + Vector3:new{x=0.0, y=5*mm}
local b0 = Vector3:new(ramp(1.0)); local b1 = b0 + Vector3:new{x=-10*mm, y=40*mm}
-- For the final straight section, angle continues at final angle of transition.
local x_length = 10*m_per_inch - b0.x
local beta = -1.90*math.pi/180.0
-- end of model
local c0 = Vector3:new{x=b0.x+x_length, y=b0.y+x_length*math.tan(beta)}
local c1 = Vector3:new{x=c0.x, y=b1.y}
--
wedge = makePatch{north=Line:new{p0=a1, p1=b1},
                  east=Line:new{p0=b0, p1=b1},
                  south=LuaFnPath:new{luaFnName="ramp"},
                  west=Line:new{p0=a0, p1=a1}}
tail = CoonsPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}
--
-- Mesh the patches, with particular discretisation.
local rcfx = RobertsFunction:new{end0=true, end1=false, beta=1.2}
local rcfy = RobertsFunction:new{end0=true, end1=false, beta=1.1}
local factor = 1  -- was 2
local ni0 = math.floor(240*factor)
local ni1 = math.floor(60*factor)
local nj0 = math.floor(40*factor)
--
-- We split the patches into roughly equal blocks so that
-- we make good use of our multicore machines.
registerFluidGridArray{
   grid=StructuredGrid:new{psurface=wedge, niv=nj0+1, njv=nj0+1},
   nib=10, njb=2,
   fsTag="initial",
   bcTags={west="inflow", north="inflow", south="noslipwall"}
}
registerFluidGridArray{
   grid=StructuredGrid:new{psurface=tail, niv=ni1+1, njv=nj0+1},
   nib=4, njb=2,
   fsTag="initial",
   bcTags={east="outflow", north="inflow", south="noslipwall"}
}
identifyGridConnections()

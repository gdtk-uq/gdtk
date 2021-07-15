-- fjl-gmsh.lua
-- Generate a file containing Gmsh loops for the fluidic oscillator -- longer body
-- PJ, 2019-08-01
-- $ e4shared --custom-script --script-file=fjl-gmsh.lua
--
f = assert(io.open("fjl-gmsh.txt", "w"))

pTag, cTag, lTag = 0, 0, 0 -- points, curves, loops
outerPath = "M0.025,-0.010;V0.010;H0.0;V0.001;".. -- first vertical segment is outflow
   "C-0.001,0.001 -0.002,0.003 -0.004,0.003;"..
   "C-0.006,0.003 -0.006,0.0016 -0.008,0.0016;"..
   "C-0.010,0.0016 -0.010,0.003 -0.012,0.003;"..
   "H-0.014;"..
   "C-0.0156,0.003 -0.017,0.0018 -0.017,0.0;"..
   "C-0.017,-0.0018 -0.0156,-0.003 -0.014,-0.003;"..
   "H-0.012;"..
   "C-0.010,-0.003 -0.010,-0.0016 -0.008,-0.0016;"..
   "C-0.006,-0.0016 -0.006,-0.003 -0.004,-0.003;"..
   "C-0.002,-0.003 -0.001,-0.001 0.0,-0.001;"..
   "V-0.010;H0.025;Z"
svgPth1 = SVGPath:new{path=outerPath}
pTag, cTag, lTag, gmshStr = svgPth1:toGmshString(pTag, cTag, lTag,
                                                 "SLIP_WALL", 5.0e-5)
f:write(gmshStr)

injectorPath = "M-0.0126,-0.0006;V0.0006;".. -- first vertical segment is inflow
   "H-0.012;V0.0016;H-0.014;C-0.015,0.0016 -0.0156,0.001 -0.0156,0.0;"..
   "C-0.0156,-0.001 -0.015,-0.0016 -0.014,-0.0016;H-0.012;V-0.0006;H-0.0126;Z"
svgPth2 = SVGPath:new{path=injectorPath}
pTag, cTag, lTag, gmshStr = svgPth2:toGmshString(pTag, cTag, lTag,
                                                 "SLIP_WALL", 5.0e-5)
f:write(gmshStr)

bbStr = [[// Make sure that the boundaries appear in the grid file.
Physical Curve("INFLOW", 91) = {17}; // injector inflow line
Physical Curve("OUTFLOW", 92) = {1}; // outflow line for whole domain
Plane Surface(1) = {1, 2};
Physical Surface("my surface", 1) = {1};
]]
f:write(bbStr)
f:close()

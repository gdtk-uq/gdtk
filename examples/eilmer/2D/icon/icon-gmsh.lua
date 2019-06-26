-- icon-gmsh.lua
-- Generate a file containing Gmsh loops for the Eilmer iconic "E".
-- PJ, 2019-06-19
-- $ e4shared --custom-script --script-file=icon-gmsh.lua
--
f = assert(io.open("icon-gmsh.txt", "w"))
bbStr = '// bounding box\n'..
   'len = 0.05;\n'..
   'Point(1) = {-2.0, -2.0, 0.0, len};\n'..
   'Point(2) = {7.0, -2.0, 0.0, len};\n'..
   'Point(3) = {7.0, 7.0, 0.0, len};\n'..
   'Point(4) = {-2.0, 7.0, 0.0, len};\n'..
   'Line(1) = {1, 2}; // south\n'..
   'Line(2) = {2, 3}; // east\n'..
   'Line(3) = {3, 4}; // north\n'..
   'Line(4) = {4, 1}; // west\n'..
   'Curve Loop(1) = {4, 1, 2, 3};\n'..
   'Physical Curve("SLIPWALL", 5) = {1, 3}; // north and south\n'..
   'Physical Curve("INFLOW", 6) = {2}; // east\n'..
   'Physical Curve("OUTFLOW", 7) = {4}; // west\n'
f:write(bbStr)
--
labels = {-- "bounding box",
          "TOP_BAR", "MIDDLE_BAR", "LOWER_BAR"}
svgStrings = {
   -- "M-2.0,-2.0;h9.0;v9.0;h-9.0;Z",
   "M1.0,4.0;h2.5;q1.25,0.0 1.5,1.0;h-2.5;q-1.25,0.0 -1.5,-1.0;Z",
   "M0.5,2.0;h1.5;q1.25,0.0 1.5,1.0;h-1.5;q-1.25,0.0 -1.5,-1.0;Z",
   "M0.0,0.0;h2.5;q1.25,0.0 1.5,1.0;h-2.5;q-1.25,0.0 -1.5,-1.0;Z"
}
pTag, cTag, lTag = 4, 4, 1
for i,svgStr in ipairs(svgStrings) do
   svgPth = SVGPath:new{path=svgStr}
   pTag, cTag, lTag, gmshStr = svgPth:toGmshString(pTag, cTag, lTag, labels[i], 0.05)
   f:write(gmshStr)
end
--
f:write('Plane Surface(1) = {1, 2, 3, 4};\n')
f:write('Physical Surface("my surface", 1) = {1};\n')
f:close()

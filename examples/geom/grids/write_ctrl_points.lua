function writeCtrlPts2VTK(pts, fname)
   local niv = #pts
   local njv = #(pts[1])
   local f = assert(io.open(fname, "w"))
   f:write("# vtk DataFile Version 2.0\n")
   f:write("\n")
   f:write("ASCII\n")
   f:write("\n")
   f:write("DATASET STRUCTURED_GRID\n")
   f:write(string.format("DIMENSIONS %d %d 1\n", niv, njv))
   f:write(string.format("POINTS %d float\n", (niv * njv)))
   for j=1,njv do
      for i=1,niv do
         local p = pts[i][j]
         f:write(string.format("%.18e %.18e %.18e\n", p.x, p.y, 0.0))
      end
   end
   f:close()
end

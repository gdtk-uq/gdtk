-- Author: Rowan J. Gollan
-- Date: 2019-04-30

-- We know the reference shock angle is 48.96
beta_ref = math.rad(48.96)

-- Read in original grid
grid = UnstructuredGrid:new{filename="cone20.su2", fmt="su2text"}

-- Build a list of all cells on shock
p0 = Vector3:new{x=0.2, y=0.0}
p1 = Vector3:new{x=1.0, y=(1.0-0.2)*math.tan(beta_ref)}
shock = Line:new{p0=p0, p1=p1}

-- a line just before shock
p0b = Vector3:new{p0}
p0b.x = p0b.x - 0.01
p1b = Vector3:new{p1}
p1b.x = p1b.x - 0.01
bShock = Line:new{p0=p0b, p1=p1b}

-- a line just after shock
p0a = Vector3:new{p0}
p0a.x = p0b.x + 0.01
p1a = Vector3:new{p1}
p1a.x = p1a.x + 0.01
aShock = Line:new{p0=p0a, p1=p1a}

cellsToRefine = {}
npts = 100
dt = 1.0/(npts-1)
for i=1,npts do
   t = (i-1)*dt

   p = shock(t)
   cellId = grid:find_nearest_cell_centre{x=p.x, y=p.y}
   cellsToRefine[cellId] = true

   p = bShock(t)
   cellId = grid:find_nearest_cell_centre{x=p.x, y=p.y}
   cellsToRefine[cellId] = true

   p = aShock(t)
   cellId = grid:find_nearest_cell_centre{x=p.x, y=p.y}
   cellsToRefine[cellId] = true
end

noCells = 0
for _,__ in pairs(cellsToRefine) do noCells = noCells + 1 end

print("no. of cells that intersect with shock: ", noCells)

-- Now actually split cells.
for id,_ in pairs(cellsToRefine) do
   grid:splitTriangleCell(id)
end

-- And write new grid to disk
grid:write_to_su2_file("cone20-refined-near-shock.su2")



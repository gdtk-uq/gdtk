
-- Grid Geometry Specs   --            bo   ___---* co
-- ri = 0.1524           --        _-- *---/      |    -- Notes: - ai is always at the origin
-- L = 1.295 - ri
dofile('conditions.lua')
thetai = math.rad(9.0)   --       /               |              - ci and co are at x=L
ro = 0.25                --      /    _-*---------* ci           - oo is downstream of oi by diffo
diffo = 0.08             --     /    /  bi           
thetao = math.rad(25.0)  --     |    |
L = L+ri                 --  ao *----*ai *oi *oo

oi = Vector3:new{x=ri, y=0.0}
ai = oi + ri*Vector3:new{x=-1.0, y=0.0}
bi = oi + ri*Vector3:new{x=-math.sin(thetai), y=math.cos(thetai)}
dxci = L - bi.x
dyci = dxci*math.tan(thetai)
ci = Vector3:new{x=bi.x + dxci, y=bi.y+dyci}
aibi = Arc:new{p0=ai, p1=bi, centre=oi}
bici = Line:new{p0=bi, p1=ci}
aici = Polyline:new{segments={aibi, bici}}


aoco = Spline2:new{filename="initial_shock_shape.dat"}
ao = aoco(0.0)
co = aoco(1.0)

aoai = Line:new{p0=ao, p1=ai}
coci = Line:new{p0=co, p1=ci}

niv = 20+1; njv = 30+1;
surface = CoonsPatch:new{north=coci, south=aoai, east=aici, west=aoco}

cluster_east = GeometricFunction:new{
   a = 1/njv/3, r=1.1, N=njv, reverse=false
}
cluster_west = GeometricFunction:new{
   a = 1/njv/3, r=1.1, N=njv, reverse=false
}
clusterlist = {
   north=cluster_north, 
   south=cluster_south,
   east=cluster_east, 
   west=cluster_west
}

grid = StructuredGrid:new{
   psurface=surface,
   niv=niv,
   njv=njv,
   cfList=clusterlist
}

registerFluidGridArray{
   grid=grid,
   fsTag="initial",
   bcTags = {
      north = "outflow",
      west = "inflow",
      east = "wall",
      south = "symm"
   },
   nib=2, njb=4,
   shock_fitting = false,
}

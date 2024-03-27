print("Under-expanded Jet experiment by Orescanin and Austin.")
-- PJ 2021-03-27 Eilmer4 example
--    2024-03-27 Eilmer5 port
--
config.dimensions = 2
config.axisymmetric = true

-- Geometry of the flow domain.
D = 0.010 -- diameter of nozzle throat, metres
H = 0.100 -- height of flow domain
LX = 0.200 -- downstream edge of flow domain
a0 = {x=0.0, y=0.0}; a1 = {x=0.0, y=D/2}; a2 = {x=0.0, y=0.025}
a3 = {x=-0.005, y=0.040}; a4 = {x=-0.005, y=H}
b0 = {x=-0.005, y=0.0}; b1 = {x=-0.005, y=D/2}
b2 = {x=-0.010, y=0.010}; b3 = {x=-0.010, y=0.0275}
b4 = {x=-0.025, y=0.040}; b5 = {x=-0.025, y=H}
c0 = {x=-0.048, y=0.0}; c1 = {x=-0.048, y=D/2}
c2 = {x=-0.048, y=0.010}; c3 = {x=-0.048, y=0.0275}
d0 = {x=0.055, y=0.0}; d1 = {x=0.055, y=0.005}; d2 = {x=0.055, y=0.025}
d3 = {x=0.055, y=0.040}; d4 = {x=0.055, y=H}
e0 = {x=LX, y=0.0}; e1 = {x=LX, y=0.005}; e2 = {x=LX, y=0.025}
e3 = {x=LX, y=0.040}; e4 = {x=LX, y=H}

surfs = {}
surfs[0] = CoonsPatch:new{p00=b0, p10=a0, p11=a1, p01=b1}
surfs[1] = CoonsPatch:new{p00=c0, p10=b0, p11=b1, p01=c1}
surfs[2] = CoonsPatch:new{p00=c1, p10=b1, p11=b2, p01=c2}
surfs[3] = CoonsPatch:new{p00=c2, p10=b2, p11=b3, p01=c3}

surfs[4] = CoonsPatch:new{p00=a0, p10=d0, p11=d1, p01=a1}
surfs[5] = CoonsPatch:new{p00=a1, p10=d1, p11=d2, p01=a2}
surfs[6] = CoonsPatch:new{p00=a2, p10=d2, p11=d3, p01=a3}
surfs[7] = CoonsPatch:new{p00=a3, p10=d3, p11=d4, p01=a4}
surfs[8] = CoonsPatch:new{p00=b4, p10=a3, p11=a4, p01=b5}

surfs[9] = CoonsPatch:new{p00=d0, p10=e0, p11=e1, p01=d1}
surfs[10] = CoonsPatch:new{p00=d1, p10=e1, p11=e2, p01=d2}
surfs[11] = CoonsPatch:new{p00=d2, p10=e2, p11=e3, p01=d3}
surfs[12] = CoonsPatch:new{p00=d3, p10=e3, p11=e4, p01=d4}

-- Mesh the patches, with particular discretisation.
dx = 0.25e-3
ncb = math.floor(0.048/dx)
nba = math.floor(0.005/dx)
print("ncb=", ncb, "nba=", nba)
n01 = math.floor(0.005/dx)
n12 = math.floor(0.005/dx)
n23=math.floor(0.0225/dx)
print("n01=", n01, "n12=", n12, "n23=", n23)

grids = {}
grids[0] = StructuredGrid:new{psurface=surfs[0], niv=nba+1, njv=n01+1}
grids[1] = StructuredGrid:new{psurface=surfs[1], niv=ncb+1, njv=n01+1}
grids[2] = StructuredGrid:new{psurface=surfs[2], niv=ncb+1, njv=n12+1}
grids[3] = StructuredGrid:new{psurface=surfs[3], niv=ncb+1, njv=n23+1}

nad = math.floor(0.055/dx)
nn12 = math.floor(0.020/dx)
nn23 = math.floor(0.015/dx)
nn34 = math.floor((H-0.040)/4/dx)
cfy = RobertsFunction:new{end0=true, end1=false, beta=1.1}
grids[4] = StructuredGrid:new{psurface=surfs[4], niv=nad+1, njv=n01+1}
grids[5] = StructuredGrid:new{psurface=surfs[5], niv=nad+1, njv=nn12+1}
grids[6] = StructuredGrid:new{psurface=surfs[6], niv=nad+1, njv=nn23+1}
grids[7] = StructuredGrid:new{psurface=surfs[7], niv=nad+1, njv=nn34+1,
                              cfList={west=cfy, east=cfy}}

grids[8] = StructuredGrid:new{psurface=surfs[8], niv=nn12+1, njv=nn34+1,
                              cfList={west=cfy, east=cfy}}

nde = nad
cfx = RobertsFunction:new{end0=true, end1=false, beta=1.1}
grids[9] = StructuredGrid:new{psurface=surfs[9], niv=nad+1, njv=n01+1,
                              cfList={south=cfx, north=cfx}}
grids[10] = StructuredGrid:new{psurface=surfs[10], niv=nad+1, njv=nn12+1,
                               cfList={south=cfx, north=cfx}}
grids[11] = StructuredGrid:new{psurface=surfs[11], niv=nad+1, njv=nn23+1,
                               cfList={south=cfx, north=cfx}}
grids[12] = StructuredGrid:new{psurface=surfs[12], niv=nad+1, njv=nn34+1,
                               cfList={south=cfx, north=cfx,
                                       west=cfy, east=cfy}}

-- Define the flow-solution blocks and stitch them together with BCs.
registerFluidGrid{grid=grids[0], fsTag="reservoir"}
registerFluidGridArray{grid=grids[1], nib=2, njb=1, fsTag="reservoir"}
registerFluidGridArray{grid=grids[2], nib=2, njb=1, fsTag="reservoir"}
registerFluidGridArray{grid=grids[3], nib=2, njb=1, fsTag="reservoir"}

registerFluidGridArray{grid=grids[4], nib=3, njb=1, fsTag="tank"}
registerFluidGridArray{grid=grids[5], nib=3, njb=1, fsTag="tank"}
registerFluidGridArray{grid=grids[6], nib=3, njb=1, fsTag="tank"}
registerFluidGridArray{grid=grids[7], nib=3, njb=1, fsTag="tank", bcTags={north="ambient"}}

registerFluidGrid{grid=grids[8], fsTag="tank", bcTags={north="ambient"}}

registerFluidGridArray{grid=grids[9], nib=3, njb=1, fsTag="tank", bcTags={east="outflow"}}
registerFluidGridArray{grid=grids[10], nib=3, njb=1, fsTag="tank", bcTags={east="outflow"}}
registerFluidGridArray{grid=grids[11], nib=3, njb=1, fsTag="tank", bcTags={east="outflow"}}
registerFluidGridArray{grid=grids[12], nib=3, njb=1, fsTag="tank", bcTags={north="ambient", east="outflow"}}

identifyGridConnections()

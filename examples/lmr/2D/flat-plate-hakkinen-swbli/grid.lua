-- grid.lua
print("generate grids for a simulation of Hakkinen et al's 1959 experiment.")
-- Peter J. and Rowan G. 2024-03-06
--
config.dimensions = 2
factor = 1 -- Scale discretization with this value, nominally 4 for good resolution.
--
--   y
--   ^   a1---b1---c1---d1   Shock generator
--   |   |    |     |    |
--   |   | 0  |  1  | 2  |   patches
--   |   |    |     |    |
--   0   a0---b0---c0---d0   Flat plate with boundary layer
--
--             0---> x
mm = 1.0e-3 -- metres per mm
-- Leading edge of shock generator and inlet to the flow domain.
L1 = 10.0*mm; H1 = 37.36*mm
a0 = Vector3:new{x=-L1, y=0.0}
a1 = a0+Vector3:new{x=0.0,y=H1}
-- Angle of inviscid shock generator.
alpha = 3.09*math.pi/180.0
tan_alpha = math.tan(alpha)
-- Start of flat plate with boundary layer.
b0 = Vector3:new{x=0.0, y=0.0}
b1 = b0+Vector3:new{x=0.0,y=H1-L1*tan_alpha}
-- End of shock generator is only part way long the plate.
L3 = 67*mm
c0 = Vector3:new{x=L3, y=0.0}
c1 = c0+Vector3:new{x=0.0,y=H1-(L1+L3)*tan_alpha}
-- End of plate, and of the whole flow domain.
L2 = 90.0*mm
d0 = Vector3:new{x=L2, y=0.0}
d1 = d0+Vector3:new{x=0.0,y=H1}
-- Now, define the three patches.
patch0 = CoonsPatch:new{p00=a0, p10=b0, p11=b1, p01=a1}
patch1 = CoonsPatch:new{p00=b0, p10=c0, p11=c1, p01=b1}
patch2 = CoonsPatch:new{p00=c0, p10=d0, p11=d1, p01=c1}
--
-- Discretization of the flow domain.
--
-- We want to cluster the cells toward the surface of the flat plate.
-- where the boundary layer will be developing.
rcf = RobertsFunction:new{end0=true,end1=true,beta=1.1}
ni0 = math.floor(20*factor); nj0 = math.floor(80*factor)
grid0 = StructuredGrid:new{psurface=patch0, niv=ni0+1, njv=nj0+1,
                           cfList={east=rcf,west=rcf}}
grid1 = StructuredGrid:new{psurface=patch1, niv=7*ni0+1, njv=nj0+1,
                           cfList={east=rcf,west=rcf}}
grid2 = StructuredGrid:new{psurface=patch2, niv=2*ni0+1, njv=nj0+1,
                           cfList={east=rcf,west=rcf}}
--
-- Build the flow blocks and attach boundary tags.
--
registerFluidGridArray{grid=grid0, fsTag='initial', nib=1, njb=2,
                       bcTags={west='inflow'}}
registerFluidGridArray{grid=grid1, fsTag='initial', nib=7, njb=2,
                       bcTags={south='noslipwall'}}
registerFluidGridArray{grid=grid2, fsTag='initial', nib=2, njb=2,
                       bcTags={south='noslipwall', east='outflow'}}
identifyGridConnections()

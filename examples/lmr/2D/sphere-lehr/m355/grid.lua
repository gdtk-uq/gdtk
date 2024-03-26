-- Rowan J. Gollan and Peter J: 2024-03-26, Ported to Eilmer5
--
-- This script is used to setup a simlulation of Lehr's
-- hemispherical projectile fired into a detonable gas.
--
-- Reference:
--   Lehr, H. (1972)
--   Acta Astronautica, 17, pp.589--597
--
print("Lehr experiment M=3.55 -- set up grid")
R = 7.5e-3 -- nose radius, metres
config.axisymmetric = true

a = {x=0.0, y=0.0}
b = {x=-R, y=0.0}
c = {x=0.0, y=R}
d = {{x=-1.5*R,y=0.0}, {x=-1.5*R,y=R}, {x=-R,y=2*R}, {x=0.0,y=3*R}}

psurf = CoonsPatch:new{
   north=Line:new{p0=d[#d], p1=c}, east=Arc:new{p0=b, p1=c, centre=a},
   south=Line:new{p0=d[1], p1=b}, west=Bezier:new{points=d}
}
ni = 128; nj = 256
-- Shock-fitting is coordinated across four blocks.
registerFluidGridArray{
   grid=StructuredGrid:new{psurface=psurf, niv=ni+1, njv=nj+1},
   nib=1, njb=4,
   fsTag="initial",
   shock_fitting=true,
   bcTags={west="inflow_sf", north="outflow"}
}

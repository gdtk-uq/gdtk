-- grid.lua
print("Riggins drag reduction calculation -- set up grid")
dofile("parameters.lua")

a = {x=R, y=0.0}
b = {x=0.0, y=0.0}
c = {x=R, y=R}
d = {{x=-h,y=0.0}, {x=-h,y=R}, {x=-R,y=2.5*R}, {x=R,y=3*R}}
e = {x=xmax, y=R}
f = {x=xmax, y=d[4].y}

psurf0 = CoonsPatch:new{
   north=Line:new{p0=d[#d], p1=c}, east=Arc:new{p0=b, p1=c, centre=a},
   south=Line:new{p0=d[1], p1=b}, west=Bezier:new{points=d}
}
-- The north edge of psurf0 corresponds to west edge of psurf1
psurf1 = CoonsPatch:new{p00=c, p10=e, p11=f, p01=d[4]}

ni0 = 256; nj0 = 256
registerFluidGridArray{
   grid=StructuredGrid:new{psurface=psurf0, niv=ni0+1, njv=nj0+1},
   nib=2, njb=2,
   fsTag="initial",
   bcTags={west="inflow", east="body"}
}
registerFluidGridArray{
   grid=StructuredGrid:new{psurface=psurf1, niv=math.floor(ni0/2)+1, njv=ni0+1},
   nib=1, njb=2,
   fsTag="initial",
   bcTags={north="inflow", south="body", east="outflow"}
}
identifyGridConnections()

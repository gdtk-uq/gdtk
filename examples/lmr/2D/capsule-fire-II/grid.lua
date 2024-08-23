-- Original notes:
-- 
-- -- Simulation of Fire II Flight Experiments
-- -- From doi.org/10.2514/6.2007-605
-- -- @author: Nick N. Gibbons 
-- -- 2020-11-30
--
-- Updates: 2024-08-23
--          RJG, ported for eilmer v5
--               changed to use control-point patch
--

config.dimensions = 2

Ri = 0.9347             
ri = 0.0102
A  = 0.3358
thetai = math.asin((A-ri)/(Ri-ri))
L = 0.05
diffo = 0.07
Ro = Ri + diffo
ro = ri + diffo
psio = math.rad(5.0)

oi = Vector3:new{x=Ri, y=0.0}
ai = oi + Ri*Vector3:new{x=-1.0, y=0.0}
bi = oi + Ri*Vector3:new{x=-math.cos(thetai), y=math.sin(thetai)}
pi = oi + (Ri-ri)*Vector3:new{x=-math.cos(thetai), y=math.sin(thetai)}
ci = pi + ri*Vector3:new{x=0.0, y=1.0}
di = ci + L*Vector3:new{x=1.0, y=0.0}

aibi = Arc:new{p0=ai, p1=bi, centre=oi}
bici = Arc:new{p0=bi, p1=ci, centre=pi}
cidi = Line:new{p0=ci, p1=di}
aidi = ArcLengthParameterizedPath:new{underlying_path=Polyline:new{segments={aibi, bici, cidi}}}

thetao = 1.5*thetai
oo = Vector3:new{x=Ri, y=0.0}
ao = oo + Ro*Vector3:new{x=-1.0, y=0.0}
do_= oo + Ro*Vector3:new{x=-math.cos(thetao), y=math.sin(thetao)}

aodo = Arc:new{p0=ao, p1=do_, centre=oo}
aoai = Line:new{p0=ao, p1=ai}
dodi = Line:new{p0=do_, p1=di}

ncpi = 4
ncpj = 13
surface = ControlPointPatch:new{north=dodi, south=aoai, east=aidi, west=aodo,
                                    ncpi=ncpi, ncpj=ncpj, guide_patch="channel-e2w"}

-- Adjust some control points to obtain better quality in shoulder region
surface:setCtrlPt(1, 11, {x=0.075, y=0.46})
surface:setCtrlPt(2, 11, {x=0.096, y=0.375})
surface:setCtrlPt(2, 10, {x=0.055, y=0.36})

cluster_east = GaussianFunction:new{m=0.001, s=0.7, ratio=0.5}
cluster_west = GaussianFunction:new{m=0.001, s=0.5, ratio=0.5}
clusterlist = {north=none, south=none, east=cluster_east, west=cluster_west}

niv =32+1; njv =128+1;
grid = StructuredGrid:new{psurface=surface, niv=niv, njv=njv, cfList=clusterlist}

registerFluidGridArray{
   grid=grid,
   nib=2, njb=8,
   fsTag='initial',
   shock_fitting=true,
   bcTags={west='inflow_sf', north='outflow'}
}






# double-ramp.py
# A simple 3D simulation of flow over a double ramp.
# Christine 2022-11-04 adapted from the Eilmer 2D sharp cone example.
# 
#  ^j, y
#  |
# z1-a1--b1--c1--d1   (H)
# |  |   |   |   |
# |  |   |   |   |
# |  |   |   | 2 |
# |  |   | 1 c0__d0
# |  | 0 |   /
# |  |   |  /
# |  |   | /
# |  |   b0
# |  |  /    
# |  | / 
# z0-a0-----------L   -->i,x
#
#
mm = 0.001  # metres
b0_x = 92.08*mm
b0_y = 42.94*mm
c0_x = 153.69*mm
c0_y = 130.925*mm
z0_x = -150.0*mm
a_x = -20*mm
Box_H = 350.0*mm
Box_L = 250*mm
H = Vector3(0.0, 0.0, 0.01)
z0 = Vector3(z0_x, 0.0, 0.0); 
z1 = Vector3(z0_x, Box_H, 0.0);
a0 = Vector3(a_x, 0.0, 0.0); 
a1 = Vector3(a_x, Box_H, 0.0);
b0 = Vector3(b0_x, b0_y, 0.0); b1 = Vector3(b0_x, Box_H, 0.0); 
c0 = Vector3(c0_x, c0_y, 0.0);c1 = Vector3(c0_x, Box_H, 0.0);
d0 = Vector3(Box_L, c0_y, 0.0);d1 = Vector3(Box_L, Box_H, 0.0);
vol0 = TFIVolume(p000=z0, p100=a0, p110=a1, p010=z1,
                 p001=z0+H, p101=a0+H, p111=a1+H, p011=z1+H)
grd0 = StructuredGrid(pvolume=vol0, niv=94, njv=98, nkv=3)
vol1 = TFIVolume(p000=a0, p100=b0, p110=b1, p010=a1,
                 p001=a0+H, p101=b0+H, p111=b1+H, p011=a1+H)
grd1 = StructuredGrid(pvolume=vol1, niv=58, njv=98, nkv=3)
vol2 = TFIVolume(p000=b0, p100=c0, p110=c1, p010=b1,
                 p001=b0+H, p101=c0+H, p111=c1+H, p011=b1+H)
grd2 = StructuredGrid(pvolume=vol2, niv=58, njv=98, nkv=3)
vol3 = TFIVolume(p000=c0, p100=d0, p110=d1, p010=c1,
                 p001=c0+H, p101=d0+H, p111=d1+H, p011=c1+H)
grd3 = StructuredGrid(pvolume=vol3, niv=58, njv=98, nkv=3)
#
initial = FlowState(p=5595.0, T=304.0)
inflow = FlowState(p=95.84e3, T=1103.0, velx=1000.0)
#
b0 = FluidBlock(i=0, grid=grd0, initialState=initial,                	    			bcs={'iminus':InflowBC(inflow),'iplus':ExchangeBC(),
		'jminus':WallWithSlipBC(),'jplus':OutflowBC()})
b1 = FluidBlock(i=1, grid=grd1, initialState=initial,                	    			bcs={'iminus':ExchangeBC(),'iplus':ExchangeBC(),
		'jplus':OutflowBC()})
b2 = FluidBlock(i=2, grid=grd2, initialState=initial,
                bcs={'iminus':ExchangeBC(),'iplus':ExchangeBC(),
                'jplus':OutflowBC()})
b3 = FluidBlock(i=3, grid=grd3, initialState=initial,
                bcs={'iminus':ExchangeBC(),'iplus':OutflowBC(),
                'jplus':OutflowBC()})
#
config.max_time = 5.0e-3
config.max_step = 50000
add_cfl_value(0.0, 0.5)
add_dt_plot(0.0, 1.0e-5)


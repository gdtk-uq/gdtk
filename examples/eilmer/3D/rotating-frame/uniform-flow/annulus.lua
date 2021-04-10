-- annulus.lua
-- Test case for 3D flow in a rotating frame.
--
-- Peter J. & Jens Kunze 2020-08-27

config.dimensions = 3
Ro = 0.100 -- outer radius of annulus, metres
Ri = 0.080 -- inner radius
L = 0.100  -- (axial) length of annulus

-- Free-stream properties
T_inf = 300.0   -- degrees K
p_inf = 100.0e3 -- Pa
V_inf = 400.0   -- m/s supersonic, just
config.title = "Rotating frame with uniform velocity flowing through it."
print(config.title)

nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
inflow = FlowState:new{p=p_inf, T=T_inf, velz=V_inf}

a = Vector3:new{x=Ro, y=0, z=0}; b = Vector3:new{x=Ri, y=0, z=0}
c = Vector3:new{x=0, y=0, z=0}
d = Vector3:new{x=0, y=Ro, z=0}; e = Vector3:new{x=0, y=Ri, z=0}
f = Vector3:new{x=-Ro, y=0, z=0}; g = Vector3:new{x=-Ri, y=0, z=0}
h = Vector3:new{x=0, y=-Ro, z=0}; i = Vector3:new{x=0, y=-Ri, z=0}

bottom_face0 = CoonsPatch:new{south=Arc:new{p0=e, p1=b, centre=c},
                              north=Arc:new{p0=d, p1=a, centre=c},
                              west=Line:new{p0=e, p1=d},
                              east=Line:new{p0=b, p1=a}}
e2 = e + Vector3:new{x=0, y=0, z=L}
vol0 = SweptSurfaceVolume:new{face0123=bottom_face0, edge04=Line:new{p0=e, p1=e2}}
grid0 = StructuredGrid:new{pvolume=vol0, niv=21, njv=5, nkv=11}
blk0 = FBArray:new{grid=grid0, initialState=inflow, omegaz=10.0,
		       bcList={bottom=InFlowBC_Supersonic:new{flowState=inflow},
			       top=OutFlowBC_Simple:new{}},
                       nkb=2}

bottom_face1 = CoonsPatch:new{south=Arc:new{p0=g, p1=e, centre=c},
                              north=Arc:new{p0=f, p1=d, centre=c},
                              west=Line:new{p0=g, p1=f},
                              east=Line:new{p0=e, p1=d}}
g2 = g + Vector3:new{x=0, y=0, z=L}
vol1 = SweptSurfaceVolume:new{face0123=bottom_face1, edge04=Line:new{p0=g, p1=g2}}
grid1 = StructuredGrid:new{pvolume=vol1, niv=21, njv=5, nkv=11}
blk1 = FBArray:new{grid=grid1, initialState=inflow, omegaz=10.0,
		       bcList={bottom=InFlowBC_Supersonic:new{flowState=inflow},
			       top=OutFlowBC_Simple:new{}},
                       nkb=2}

bottom_face2 = CoonsPatch:new{south=Arc:new{p0=i, p1=g, centre=c},
                              north=Arc:new{p0=h, p1=f, centre=c},
                              west=Line:new{p0=i, p1=h},
                              east=Line:new{p0=g, p1=f}}
i2 = i + Vector3:new{x=0, y=0, z=L}
vol2 = SweptSurfaceVolume:new{face0123=bottom_face2, edge04=Line:new{p0=i, p1=i2}}
grid2 = StructuredGrid:new{pvolume=vol2, niv=21, njv=5, nkv=11}
blk2 = FBArray:new{grid=grid2, initialState=inflow, omegaz=10.0,
		       bcList={bottom=InFlowBC_Supersonic:new{flowState=inflow},
			       top=OutFlowBC_Simple:new{}},
                       nkb=2}

bottom_face3 = CoonsPatch:new{south=Arc:new{p0=b, p1=i, centre=c},
                              north=Arc:new{p0=a, p1=h, centre=c},
                              west=Line:new{p0=b, p1=a},
                              east=Line:new{p0=i, p1=h}}
b2 = b + Vector3:new{x=0, y=0, z=L}
vol3 = SweptSurfaceVolume:new{face0123=bottom_face3, edge04=Line:new{p0=b, p1=b2}}
grid3 = StructuredGrid:new{pvolume=vol3, niv=21, njv=5, nkv=11}
blk3 = FBArray:new{grid=grid3, initialState=inflow, omegaz=10.0,
		       bcList={bottom=InFlowBC_Supersonic:new{flowState=inflow},
			       top=OutFlowBC_Simple:new{}},
                       nkb=2}

identifyBlockConnections()

config.max_time = 1.0e-3
config.max_step = 4000
config.dt_init = 1.0e-6

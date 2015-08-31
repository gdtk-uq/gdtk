title = "Method of Manufactured Solutions, coupled domains case."
print(title)
config.dimensions = 2
config.title = title

-- Parse number of cells in 'i' direction from file.
file = io.open("case.txt", "r")
ncells = tonumber(file:read("*line"))
file:close()

L = 1.0
R = 287.0
k_g = 10000.0
k_s = 10*k_g

---- Analytical solution ------
local cos = math.cos
local sin = math.sin
local pi = math.pi

rho0=1.0; rhox=0.1; rhoy=0.15; rhoxy=0.08; arhox=0.75; arhoy=1.0; arhoxy=1.25;
u0 = 1.0; ux = 0.1; uy = u0; uxy = -ux; aux = 5.0/3; auy = -1.0; auxy = aux;
v0 = 0.9; vx = -0.02; vy = -v0; vxy = -vx; avx = 1.5; avy = 0.5; avxy = avx;	
T0 = 350; Tx = -10.0; Ty = 25.0; aTx = 1.5; aTy = 1.0; Ti = 350.0; aTx2 = 0.75;

function rho(x, y)
   return rho0 + rhox*sin(arhox*pi*x/L) + rhoy*cos(arhoy*pi*y/L) + rhoxy*cos(arhoxy*pi*x*y/(L*L));
end

function u(x, y)
   return u0 + ux*cos(aux*pi*x/L) + uy*cos(auy*pi*y/L) + uxy*cos(auxy*pi*x*y/(L*L))
end

function v(x, y)
   return  v0 + vx*cos(avx*pi*x/L) + vy*sin(avy*pi*y/L) + vxy*cos(avxy*pi*x*y/(L*L))
end

function T(x, y)
   return T0 + Tx*cos(aTx*pi*x/L) + Ty*cos(aTx2*pi*x/L)*sin(aTy*pi*y/L)
end

function Ts(x, y)
   return T0 + Tx*cos(aTx*pi*x/L) +  Ty*(k_g/k_s)*cos(aTx2*pi*x/L)*sin(aTy*pi*y/L)
end
----- end: Analytical solution -----

setGasModel('very-viscous-air.lua')

--[=[
-- Uncomment these functions if the analytical solution
-- is desired as the initial condition
function gasFillFn(x, y, z)
   T_g = T(x, y)
   rho_g = rho(x, y)
   p_g = rho_g*R*T_g
   return FlowState:new{p=p_g, T=T_g, velx=u(x, y), vely=v(x, y)}
end

function solidFillFn(x, y, z)
   return Ts(x, y)
end
--]=]

initT = 350.0
initial = FlowState:new{p=1.0e5, T=initT, velx=0.0, vely=0.0}

a = Vector3:new{0.0, 0.0}
b = Vector3:new{L, 0.0}
c = Vector3:new{0.0, L}
d = Vector3:new{L, L}
e = Vector3:new{0.0, 1.5*L}
f = Vector3:new{L, 1.5*L}

nx0 = ncells
ny0 = ncells
ny1 = math.floor(ny0/2)

grid0 = StructuredGrid:new{psurface=makePatch{Line:new{c,d}, Line:new{b,d}, Line:new{a,b}, Line:new{a,c}},
			  niv=nx0+1, njv=ny0+1}
grid1 = StructuredGrid:new{psurface=makePatch{Line:new{e,f}, Line:new{d,f}, Line:new{c,d}, Line:new{c, e}},
			   niv=nx0+1, njv=ny1+1}

blk0 = SBlock:new{grid=grid0, fillCondition=initial, label="blk0"}

blk0.bcList[north] = AdjacentToSolidBC:new{otherBlock=0,
					   otherFace="south",
					   orientation=-1}
blk0.bcList[east] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
blk0.bcList[south] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
blk0.bcList[west] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}

blk1 = SSolidBlock:new{grid=grid1, initTemperature=initT,
		       properties={rho=10000, k=100000, Cp=100}}

-- Set boundary conditions
blk1.bcList[north] = SolidUserDefinedBC:new{fileName='udf-bc.lua'}
blk1.bcList[east] = SolidUserDefinedBC:new{fileName='udf-bc.lua'}
blk1.bcList[south] = SolidAdjacentToGasBC:new{otherBlock=0,
					      otherFace="north",
					      orientation=-1}
blk1.bcList[west] = SolidUserDefinedBC:new{fileName='udf-bc.lua'}


config.interpolation_order = 2
config.gasdynamic_update_scheme = "euler"
config.flux_calculator = 'ausm_plus_up'
config.viscous = true
config.max_time = 4.0 -- seconds
config.max_step = 10000000
config.dt_init = 5.0e-6
config.fixed_time_step = true
config.dt_plot = config.max_time/40
--config.apply_limiter = false
--config.extrema_clipping = false

config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.udf_solid_source_terms = true
config.udf_solid_source_terms_file = 'udf-solid-source-terms.lua'

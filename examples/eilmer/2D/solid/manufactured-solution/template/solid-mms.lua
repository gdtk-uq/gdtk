-- A Method of Manufactured Solutions case to exercise the solid domain solver
--
-- Author: Rowan J. Gollan
-- Date: Aug-2015
--
-- General notes:
-- + The solution used here in the solid domain is the same
--   as the solution we used in the coupled case. What differs
--   is that we explicitly set the boundary condition on the
--   SOUTH boundary. There is no connection to the gas domain.
--
-- + Speaking of the gas domain, note that it is only a dummy
--   flow domain. This is because the eilmer4 implementation 
--   requires at least one flow block to be present. The grid
--   is only 2x2 cells and the flow solution in that block is
--   not meaningful.
--

title = "Method of Manufactured Solutions, solid only."
print(title)
config.dimensions = 2
config.title = title

-- Parse "case.txt" to find the number of cells in 'i' direction.
file = io.open("case.txt", "r")
ncells = tonumber(file:read("*line"))
file:close()

L = 1.0
R = 287.0
k_g = 10000.0
k_s = 10*k_g

-- Note the analytical solution is only here if one desires 
-- an initial condition based on using the analytical solution.
-- It is commented out since we use an initial fill temperature
-- of 350.0 K everywhere.
--[==[
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
--]==]

setGasModel('very-viscous-air.lua')
initT = 350.0
initial = FlowState:new{p=1.0e5, T=initT, velx=0.0, vely=0.0}

--[==[
-- Analytical fill functions, if required 
function gasFillFn(x, y, z)
   T_g = T(x, y)
   rho_g = rho(x, y)
   p_g = rho_g*R*T_g
   return FlowState:new{p=p_g, T=T_g, velx=u(x, y), vely=v(x, y)}
end

function solidFillFn(x, y, z)
   return Ts(x, y)
end
--]==]

-- ^^^ end of analytical fill condition.

nx = ncells
ny = math.floor(ncells/2)

c = Vector3:new{0.0, L}
d = Vector3:new{L, L}
e = Vector3:new{0.0, 1.5*L}
f = Vector3:new{L, 1.5*L}

-- Yes grid0 appears to lie on top of grid1
-- grid0 is a dummy grid with 2x2 cells to support the dummy flow domain.
grid0 = StructuredGrid:new{psurface=makePatch{Line:new{e,f}, Line:new{d,f}, Line:new{c,d}, Line:new{c, e}},
			   niv=3, njv=3}
grid1 = StructuredGrid:new{psurface=makePatch{Line:new{e,f}, Line:new{d,f}, Line:new{c,d}, Line:new{c, e}},
			   niv=nx+1, njv=ny+1}

blk0 = SBlock:new{grid=grid0, fillCondition=initial, label="blk0"}

blk0.bcList[north] = BoundaryCondition:new{
   preReconAction = { UserDefinedGhostCell:new{fileName='udf-bc.lua'} },
   preSpatialDerivAction = { UserDefinedInterface:new{fileName='udf-bc.lua'},
			     UpdateThermoTransCoeffs:new()
   }
}
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
blk1.bcList[south] = SolidUserDefinedBC:new{fileName='udf-bc.lua'}
blk1.bcList[west] = SolidUserDefinedBC:new{fileName='udf-bc.lua'}

config.interpolation_order = 2
config.gasdynamic_update_scheme = "euler"
config.flux_calculator = 'ausm_plus_up'
config.viscous = true
config.max_time = 4.0 -- seconds
config.max_step = 10000000
config.dt_init = 5.0e-6
config.fixed_time_step = true
config.dt_plot = config.max_time/20
config.udf_solid_source_terms = true
config.udf_solid_source_terms_file = 'udf-solid-source-terms.lua'

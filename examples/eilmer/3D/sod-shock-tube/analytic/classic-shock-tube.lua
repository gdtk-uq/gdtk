-- classic-shock-tube.lua
-- Run with a command like:
-- $ e4shared --custom-post --script-file=classic-shock-tube.lua
-- Note that you need the gas-model file in the current directory.
-- You may copy it from the adjacent sg/ folder.
-- PJ, 2018-04-27
print("Compute the flow conditions expected in the Sod shock tube.")
--
print("shock-tube fill conditions with air driving air")
gm = GasModel:new{'ideal-air-gas-model.lua'}
states = {}
for i=1,4 do states[i] = GasState:new{gm} end
states[1].p = 1.0e4; states[1].T = 278.8 -- driven tube, initial
states[2].p = 1.0e4; states[2].T = 278.8 -- intermediate, post-shock
states[3].p = 1.0e4; states[3].T = 278.8 -- intermediate, post-expansion
states[4].p = 1.0e5; states[4].T = 348.4 -- driver tube, initial
for i,state in ipairs(states) do
   gm:updateThermoFromPT(state)
   gm:updateSoundSpeed(state)
end
print("state 1:"); printValues(states[1])
--
-- For the unsteady expansion of the driver gas, regulation of the amount
-- of expansion is determined by the shock-processed test gas.
-- Across the contact surface between these gases, the pressure and velocity
-- have to match so we set up some trials of various pressures and check 
-- that velocities match.
function error_in_velocity(p3p4)
   -- Compute the velocity mismatch for a given pressure ratio across the expansion.
   -- Across the expansion, we get a test-gas velocity, V3g.
   local p3 = p3p4*states[4].p
   local V3g
   states[3], V3g = gasflow.finite_wave_dp(states[4], 0.0, 'cplus', p3)
   -- Across the contact surface.
   local p2 = p3
   print("current guess for p3 and p2=", p2)
   local V1s, V2, V2g = gasflow.normal_shock_p2p1(states[1], p2/states[1].p)
   return (V3g - V2g)/V3g
end
--
dofile("secant.lua")
p3p4 = secant(error_in_velocity, 0.1, 0.11, 1.0e-3)
print("From secant solve: p3/p4=", p3p4)
print("Expanded driver gas:")
p3 = p3p4*states[4].p
states[3], V3g = gasflow.finite_wave_dp(states[4], 0.0, 'cplus', p3)
print("V3g=", V3g)
print("state 3:"); printValues(states[3])
print("Shock-processed test gas:")
V1s, V2, V2g = gasflow.normal_shock_p2p1(states[1], p3/states[1].p)
states[2], V2, V2g = gasflow.normal_shock(states[1], V1s)
print("V1s=", V1s, "V2g=", V2g)
print("state 2:"); printValues(states[2])
assert(math.abs(V2g - V3g)/V3g < 1.0e-3, "mismatch in velocities")
--
-- Make a record for plotting against the Eilmer3 simulation data.
-- We reconstruct the expected data along a tube 0.0 <= x <= 1.0
-- at t=100us, where the diaphragm is at x=0.5.
x_centre = 0.5 -- metres
t = 600.0e-6 -- seconds
f = io.open('analytic.data', 'w')
f:write('# 1:x(m)  2:rho(kg/m**3) 3:p(Pa) 4:T(K) 5:V(m/s)\n')
print('Left end')
x = 0.0
f:write(string.format('%g %g %g %g %g\n',
                      x, states[4].rho, states[4].p, states[4].T, 0.0))
print('Upstream head of the unsteady expansion.')
x = x_centre - states[4].a * t
f:write(string.format('%g %g %g %g %g\n',
                      x, states[4].rho, states[4].p, states[4].T, 0.0))
print('The unsteady expansion in n steps.')
n = 100
dp = (states[3].p - states[4].p) / n
state = GasState:new{gm}
copyValues(states[4], state)
s = gm:entropy(states[4])
V = 0.0
p = states[4].p
for i=1,n do
   rhoa = state.rho * state.a
   dV = -dp / rhoa
   V = V + dV
   p = p + dp
   state.p = p
   gm:updateThermoFromPS(state, s)
   gm:updateSoundSpeed(state)
   x = x_centre + t * (V - state.a)
   f:write(string.format('%g %g %g %g %g\n',
                         x, state.rho, state.p, state.T, V))
end
print('Downstream tail of expansion.')
x = x_centre + t * (V3g - states[3].a)
f:write(string.format('%g %g %g %g %g\n',
                      x, states[3].rho, states[3].p, states[3].T, V3g))
print('Contact surface.')
x = x_centre + t * V3g
f:write(string.format('%g %g %g %g %g\n',
                      x, states[3].rho, states[3].p, states[3].T, V3g))
x = x_centre + t * V2g  -- should not have moved
f:write(string.format('%g %g %g %g %g\n',
                      x, states[2].rho, states[2].p, states[2].T, V2g))
print('Shock front')
x = x_centre + t * V1s  -- should not have moved
f:write(string.format('%g %g %g %g %g\n',
                      x, states[2].rho, states[2].p, states[2].T, V2g))
f:write(string.format('%g %g %g %g %g\n',
                      x, states[1].rho, states[1].p, states[1].T, 0.0))
print('Right end')
x = 1.0
f:write(string.format('%g %g %g %g %g\n',
                      x, states[1].rho, states[1].p, states[1].T, 0.0))
f:close()

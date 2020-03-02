local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi
local max = math.max
local sqrt = math.sqrt

function gasFillFn(x, y, z)
   <insert-expressions-here>
   return FlowState:new{p=p, T=T, velx=velx, vely=vely, velz=velz, tke=tke, omega=omega, mu_t=mu_t, k_t=k_t}
end

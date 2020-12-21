local sin = math.sin
local cos = math.cos
local exp = math.exp
local pi = math.pi
local max = math.max
local sqrt = math.sqrt
local tanh = math.tanh

function gasFillFn(x, y, z)
   $expressions
   return FlowState:new{p=p, T=T, velx=velx, vely=vely, velz=velz, mu_t=mu_t, k_t=k_t, nuhat=nuhat}
end

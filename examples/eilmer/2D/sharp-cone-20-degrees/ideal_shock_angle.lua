-- ideal_shock_angle.lua
-- Invoke with the command line:
-- $ e4shared --custom-post --script-file=ideal_shock_angle.lua
V1=1000.0; p1=95.84e3; T1=1103.0; theta=math.rad(20.0)
beta = idealgasflow.beta_cone(V1, p1, T1, theta)
print("beta=", math.deg(beta), "degrees")

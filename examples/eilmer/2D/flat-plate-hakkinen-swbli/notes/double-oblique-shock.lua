-- double-oblique-shock.lua
-- Estimate pressure rise across a reflected oblique shock.
-- PJ, 01-May-2013, 2016-11-01 for Lua version
-- $ e4shared --custom-post --script-file=double-oblique-shock.lua
--
print("Begin...")
M1 = 2.0
p1 = 1.0
g = 1.4
print("First shock:")
delta1 = 3.09 * math.pi/180.0
beta1 = idealgasflow.beta_obl(M1,delta1,g)
p2 = idealgasflow.p2_p1_obl(M1,beta1,g)
M2 = idealgasflow.M2_obl(M1,beta1,delta1,g)
print("   beta1=", beta1, "p2=", p2, "M2=", M2)
--
print("Reflected shock:")
delta2 = delta1
beta2 = idealgasflow.beta_obl(M2,delta2,g)
p3 = p2 * idealgasflow.p2_p1_obl(M2,beta2,g)
M3 = idealgasflow.M2_obl(M2,beta2,delta2,g)
print("   beta2=", beta2, "p3=", p3, "M3=", M3)
print("Done.")

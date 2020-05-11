model = "IdealGasAB"

-- Set some reference constants
-- We'll use these in our Eilmer input script also
R = 287.0 -- J/(kg.k)
T_ref = 300.0 -- K
rho_ref = 1.0 -- J/(kg.m^3)
S_K = 7.1247 -- m/s (in Kotov case)
S_cj = 2090.6

IdealGasAB = {
  type = "Yee-Kotov",
  R = R,
  gamma = 1.4,
  q = 25*R*T_ref,
  K_0 = 16418*S_cj/S_K,
  T_ign = 25*T_ref
}

model = 'VeryViscousAir'

Rgas = 287.0
g = 1.4
Prandtl = 1.0
Cv = Rgas/(g-1)
Cp = g * Cv
mu = 10.0
k = mu*Cp/Prandtl

VeryViscousAir = {
  k = k,
  mu = mu
}

ncellsList = [8,16,32]
fluxCalc = 'ausmdv'
derivCalc = 'least_squares'
derivLcn = 'cells'
blocking = 'multi'
threading = 'multi'
xOrder = 2

turbulence_model = 'spalart_allmaras_dwall_const'
norm = ['rho', 'p', 'T', 'vel.x', 'vel.y', 'vel.z', 'mu_t', 'k_t', 'nuhat']

explicit = 'false'
exe_path = "$DGD/bin/"
template_dir = "./template"

# constants
L = 1.0

Rgas = 287.0
gamma = 1.4
Cv = Rgas/(gamma-1)
Cp = gamma*Cv
Prandtl = 1.0
PrT = 0.89

mu = 10.0
k = Cp*mu/Prandtl

wall_distance = 1000.0

only_analysis = False
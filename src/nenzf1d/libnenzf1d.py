"""
Python wrapper for nenzf1d.

@author: Nick Gibbons
"""

from cffi import FFI
from numpy import linspace , array, int32, zeros

ffi = FFI()
ffi.cdef("""
    int cwrap_run();
    typedef struct PlainConfig PlainConfig;
    struct PlainConfig
    {
        int n_gm1_filename;
        char* gm1_filename;
        int n_gm2_filename;
        char* gm2_filename;
        int n_reactions_filename;
        char* reactions_filename;
        int n_reactions_filename2;
        char* reactions_filename2;
        int nn_species;
        int* n_species;
        char** species;
        int n_molef;
        double* molef;
        double T1;
        double p1;
        double Vs;
        double pe;
        double meq_throat;
        double ar;
        double pp_ps;
        double C;
        int n_xi;
        double* xi;
        int n_di;
        double* di;
        double x_end;
        double t_final;
        double t_inc;
        double t_inc_factor;
        double t_inc_max;
    };
    void pass_a_struct(PlainConfig cfg);
""")
lib = ffi.dlopen('libnenzf1d.so')

T1 = 300.0;
molef = array([1.0, 0.0, 0.0])
gm1_filename = "air-5sp-gm.lua"

species =  ['N2', 'O2', 'N', 'O', 'NO']
molef = {'N2': 0.79, 'O2': 0.21}
gm1_filename = "air-5sp-eq.lua"
gm2_filename = "air-5sp-1T.lua"
reactions_filename = "air-5sp-1T-reactions.lua"

T1 = 300        # K
p1 = 192.0e3    # Pa
Vs = 2268.0     # m/s
pe = 45.8e6     # Pa
ar = 271.16     # Mach 8 nozzle
pp_ps = 7.01e-3 # Check to make sure this is being used...
C = 0.96         # pPitot/(rho*v^2)

# Define the expanding part of the nozzle as a schedule of diameters with position.
# Values are sampled from M8_COORDS.txt file.
xi = array([0.0000, 5.007e-3, 1.038e-2, 1.998e-2, 5.084e-2, 0.10097, 0.20272, 0.40123, 0.60631, 0.80419, 1.110])
di = array([0.0164, 0.01676, 0.01840, 0.02330, 0.04332, 0.07457, 0.12397, 0.18691,0.22705, 0.25263, 0.27006])

molefa = zeros(len(species))
for k,v in molef.items(): molefa[species.index(k)] = v
molef = molefa

config = ffi.new("struct PlainConfig *")

# When working with FFI, it's very important that your objects created with ffi.new don't get
# garbage collected while the underlying C code is working on them.

config.n_gm1_filename = len(gm1_filename)
gm1_filename_p = ffi.new("char[]", bytes(gm1_filename, encoding="utf-8"))
config.gm1_filename = gm1_filename_p

config.n_gm2_filename = len(gm2_filename)
gm2_filename_p = ffi.new("char[]", bytes(gm2_filename, encoding="utf-8"))
config.gm2_filename = gm2_filename_p

config.n_reactions_filename = len(reactions_filename)
reactions_filename_p = ffi.new("char[]", bytes(reactions_filename, encoding="utf-8"))
config.reactions_filename = reactions_filename_p

species_p = [ffi.new("char[]", bytes(spname, encoding="utf-8")) for spname in species]
n_species = array([len(spname) for spname in species], dtype=int32)
n_species_p = ffi.from_buffer("int[]", n_species)
nn_species = n_species.size
species_pp = ffi.new("char *[]", species_p)
config.nn_species = nn_species
config.n_species = n_species_p
config.species = species_pp

config.n_molef = molef.size;
molef_p = ffi.from_buffer("double[]", molef)
config.molef = molef_p

config.T1 = T1;
config.p1 = p1;
config.Vs = Vs;
config.pe = pe;
config.ar = ar;
config.pp_ps = pp_ps;
config.C = C;

config.n_xi = xi.size;
xi_p = ffi.from_buffer("double[]", xi)
config.xi = xi_p

config.n_di = di.size;
di_p = ffi.from_buffer("double[]", di)
config.di = di_p

lib.pass_a_struct(config[0])


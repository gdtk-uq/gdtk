"""
Python wrapper for nenzf1d.

@author: Nick Gibbons
"""

from pathlib import Path
from dataclasses import dataclass

from cffi import FFI
from numpy import array, zeros, ndarray

ffi = FFI()
ffi.cdef("""
    int cwrap_init();
    typedef struct PlainConfig PlainConfig;
    struct PlainConfig
    {
        char* gm1_filename;
        char* gm2_filename;
        char* reactions_filename;
        char* reactions_filename2;
        int n_species;
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
    typedef struct PlainResult PlainResult;
    struct PlainResult
    {
        double x;
        double area_ratio;
        double velocity;
        double Mach_number;
        double p_pitot;
        double rayleigh_pitot;
        double pressure, density, temperature, viscosity;
        int n_T_modes;
        double* T_modes;
        int n_massf;
        double* massf;
        double massflux_rel_err, enthalpy_rel_err, pitot_rel_err;
    };
    PlainResult run(int verbosityLevel, PlainConfig cfg);
""")

so = ffi.dlopen(str(Path(__file__).with_name("libnenzf1d.so")))
so.cwrap_init()


@dataclass
class Nenzf1d:
    species: list[str]
    molef: dict[str, float]
    gas_model_1: str
    gas_model_2: str
    reactions: str
    T1: float
    p1: float
    Vs: float
    xi: ndarray
    di: ndarray
    reaction_file2: str = ""
    pe: float = 0.0
    meq_throat: float = 1.0
    ar: float = 1.0
    pp_ps: float = 0.0
    C: float = 1.0
    x_end: float = 0.0
    t_final: float = 1.0e-2
    t_inc: float = 1.0e-10
    t_inc_factor: float = 1.0001
    t_inc_max: float = 1.0e-7

    def run(self, verbosityLevel=1):
        config = ffi.new("struct PlainConfig *")

        # When working with FFI, it's very important that your objects created with ffi.new don't get
        # garbage collected while the underlying C code is working on them.
        # That's the reason there's so many definitions of thing_p = ffi.new(...) below.

        # FIXME: Do optional parameters...
        gm1_filename_p = ffi.new("char[]", bytes(self.gas_model_1, encoding="utf-8"))
        config.gm1_filename = gm1_filename_p

        gm2_filename_p = ffi.new("char[]", bytes(self.gas_model_2, encoding="utf-8"))
        config.gm2_filename = gm2_filename_p

        reactions_filename_p = ffi.new(
            "char[]", bytes(self.reactions, encoding="utf-8")
        )
        config.reactions_filename = reactions_filename_p

        species_p = [
            ffi.new("char[]", bytes(spname, encoding="utf-8"))
            for spname in self.species
        ]
        n_species = len(species_p)
        species_pp = ffi.new("char *[]", species_p)
        config.n_species = n_species
        config.species = species_pp

        config.n_molef = self.molef.size
        molef_p = ffi.from_buffer("double[]", self.molef)
        config.molef = molef_p

        config.T1 = self.T1
        config.p1 = self.p1
        config.Vs = self.Vs
        config.pe = self.pe
        config.ar = self.ar
        config.pp_ps = self.pp_ps
        config.C = self.C
        config.n_xi = self.xi.size
        xi_p = ffi.from_buffer("double[]", self.xi)
        config.xi = xi_p

        config.n_di = self.di.size
        di_p = ffi.from_buffer("double[]", self.di)
        config.di = di_p

        result = so.run(verbosityLevel, config[0])
        self.x = result.x
        self.area_ratio = result.area_ratio
        self.velocity = result.velocity
        self.Mach_number = result.Mach_number
        self.p_pitot = result.p_pitot
        self.rayleigh_pitot = result.rayleigh_pitot
        self.pressure = result.pressure
        self.density = result.density
        self.temperature = result.temperature
        self.viscosity = result.viscosity
        if result.n_T_modes > 0:
            self.T_modes = ffi.unpack(result.T_modes, result.n_T_modes)
        else:
            self.T_modes = []
        if result.n_massf > 0:
            self.massf = ffi.unpack(result.massf, result.n_massf)
        else:
            self.massf = []
        self.massflux_rel_err = result.massflux_rel_err
        self.enthalpy_rel_err = result.enthalpy_rel_err
        self.pitot_rel_err = result.pitot_rel_err

    def print_exit_condition(self):
        print("Exit condition:")
        print("  x           {:g} m".format(self.x))
        print("  area-ratio  {:g}".format(self.area_ratio))
        print("  velocity    {:g} km/s".format(self.velocity / 1000.0))
        print("  Mach        {:g}".format(self.Mach_number))
        print("  p_pitot     {:g} kPa (C.rho.V^2)".format(self.p_pitot / 1000.0))
        if self.rayleigh_pitot > 0.0:
            print(
                "  p_pitot     {:g} kPa (Rayleigh-Pitot, frozen)".format(
                    self.rayleigh_pitot / 1000.0
                )
            )
        print("  pressure    {:g} kPa".format(self.pressure / 1000.0))
        print("  density     {:g} kg/m^3".format(self.density))
        print("  temperature {:g} K".format(self.temperature))
        for species, Y in zip(self.species, self.massf):
            print("  massf[{}]{}{:g}".format(species, " " * (5 - len(species)), Y))
        print("  viscosity   {:g} Pa.s".format(self.viscosity))
        print("Expansion error-indicators:")
        print("  relerr-mass {:g}".format(self.massflux_rel_err))
        print("  relerr-H    {:g}".format(self.enthalpy_rel_err))
        if self.rayleigh_pitot > 0.0:
            print("  relerr-pitot {:g}".format(self.pitot_rel_err))


if __name__ == "__main__":
    # Conditions matching T4 shot 11742
    T1 = 300.0
    molef = array([1.0, 0.0, 0.0])

    species = ["N2", "O2", "N", "O", "NO"]
    molef = {"N2": 0.79, "O2": 0.21}
    gas_model_1 = "air-5sp-eq.lua"
    gas_model_2 = "air-5sp-1T.lua"
    reactions = "air-5sp-1T-reactions.lua"

    T1 = 300  # K
    p1 = 192.0e3  # Pa
    Vs = 2268.0  # m/s
    pe = 45.8e6  # Pa
    ar = 271.16  # Mach 8 nozzle
    pp_ps = 7.01e-3  # Check to make sure this is being used...
    C = 0.96  # pPitot/(rho*v^2)

    # Define the expanding part of the nozzle as a schedule of diameters with position.
    # Values are sampled from M8_COORDS.txt file.
    xi = array(
        [
            0.0000,
            5.007e-3,
            1.038e-2,
            1.998e-2,
            5.084e-2,
            0.10097,
            0.20272,
            0.40123,
            0.60631,
            0.80419,
            1.110,
        ]
    )
    di = array(
        [
            0.0164,
            0.01676,
            0.01840,
            0.02330,
            0.04332,
            0.07457,
            0.12397,
            0.18691,
            0.22705,
            0.25263,
            0.27006,
        ]
    )

    molefa = zeros(len(species))
    for k, v in molef.items():
        molefa[species.index(k)] = v
    molef = molefa

    sim = Nenzf1d(
        species,
        molef,
        gas_model_1,
        gas_model_2,
        reactions,
        T1,
        p1,
        Vs,
        xi,
        di,
        pe=pe,
        ar=ar,
        pp_ps=pp_ps,
        C=C,
    )
    sim.run()
    sim.print_exit_condition()

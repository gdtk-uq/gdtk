/**
 * physical_constants.d
 *
 * We store some oft used physical constants
 * in one place for use with the gas models.
 *
 * Author: Rowan G. and Peter J.
 */

module gas.physical_constants;

///Universal gas constant (S.I. units)
immutable double R_universal = 8.31451;   // J/(mol.K) -- Tipler (1991)
immutable double R_universal_cal = 1.987; // cal/(mol.K)
// One atmosphere, in Pascals
immutable double P_atm = 101.325e3;                  // Pa
immutable double Boltzmann_constant = 1.380658e-23;  // J/K -- Tipler (1991)
immutable double Avogadro_number = 6.02214e23;
immutable double electron_volt_energy = 1.60219e-19;     // J
immutable double elementary_charge = 1.602176634e-19;    // C (Exact as of 2019 SI unit redefinition)
immutable double vacuum_permittivity = 8.8541878128e-12; // A^2.s^4/kg/m^3
immutable double Plancks_constant = 6.62607015e-34;      // J/Hz
immutable double speed_of_light = 299792458.0;           // m/s

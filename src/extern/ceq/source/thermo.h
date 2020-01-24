#ifndef thermo_h
#define thermo_h

extern const double Ru;
double compute_Cp0_R(double Tin, double* lewis);
double compute_H0_RT(double Tin, double* lewis);
double compute_S0_R(double Tin, double* lewis);
double compute_G0_RT(double T, double* lewis);

#endif

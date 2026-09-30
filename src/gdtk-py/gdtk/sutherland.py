"""
sutherland.py -- Sutherland form of viscosity and thermal conductivity for a few gases.

Species available: Air, N2, O2, H2, CO2, CO, Ar.

Reference for coefficients: White FM (2006) Viscous Fluid Flow, pp.28

.. Versions: Wilson (original)
             PJ (refactor) 29-Oct-2010
             Sphinx doc-comments 15-Jan-2013
             Port to GDTK and Python3
"""

mu_ref = {'Air':1.716e-5, 'N2':1.663e-5, 'O2':1.919e-5, 'H2':8.411e-6,
          'CO2':1.370e-5, 'CO':1.657e-5, 'Ar':1.716e-5}
T_ref = {'Air':273.0, 'N2':273.0, 'O2':273.0, 'H2':273.0,
         'CO2':273.0, 'CO':273.0, 'Ar':273.0}
S_mu = {'Air':111.0, 'N2':107.0, 'O2':139.0, 'H2':97.0,
     'CO2':222.0, 'CO':136.0, 'Ar':144.0}
k_ref = {'Air':0.0241, 'N2':0.0242, 'O2':0.0244, 'H2':0.1680,
         'CO2':0.0146, 'CO':0.0232, 'Ar':0.0163}
S_k = {'Air':194.0, 'N2':150.0, 'O2':240.0, 'H2':120.0,
       'CO2':1800.0, 'CO':180.0, 'Ar':170.0}

def sutherland(T, T_ref, S):
    """
    Sutherland's expression relating the quantity, at temperature T,
    to the value at the reference temperature.

    :param T: temperature in degrees K
    :param T_ref: reference temperature (degrees K)
    :param S: Sutherland constant (degrees K)
    :returns: ratio of property at temperature T to reference value.
    """
    return (T / T_ref)**1.5 * (T_ref + S)/(T + S)

def mu(T, gname):
    """
    Coefficient of dynamic viscosity.

    :param T: temperature in degrees K
    :param gname: string specifying gas name (for dictionary look-up)
    :returns: coefficient of dynamic viscosity, in Pa.s.
    """
    return mu_ref[gname] * sutherland(T, T_ref[gname], S_mu[gname]);

def k(T, gname):
    """
    Coefficient of thermal conductivity.

    :param T: temperature in degrees K
    :param gname: string specifying gas name (for dictionary look-up)
    :returns: thermal conductivity, in W/(m.K)
    """
    return k_ref[gname] * sutherland(T, T_ref[gname], S_k[gname]);

if __name__ == '__main__':
    print("Sutherland expression for evaluating transport coefficients.")
    gas = 'Air'; T = 300.0
    print("For gas=", gas, "T=", T, ",we compute mu=", mu(T,gas), "k=", k(T,gas))
    print("Done.")

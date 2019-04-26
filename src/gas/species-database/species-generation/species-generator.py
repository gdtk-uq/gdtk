"""
Revision 2
Author: Oliver Street
Date:   15/04/2019

Title:     Chemkin-Species to Eilmer-Species 
Author:    Samuel Dillon
Date:      05/02/2019


Rewrote code, to allow for more ready conversion
of multiple species + streamlined code.
Added the generation of Chemkin style viscosity and thermal
conductivity calculations, for use in a polynomial of form:
    
ln(X) = A + B*(ln(T)) + C*(ln(T))^2 + D*(ln(T))^3

Result:
    X = eta_k (viscosity), lambda_k (conductivity)
Input:
    T = Temperature over range given by thermo.dat
    A, B, C, D = Least squares fitted coefficients


species-generator
 ├ species-generator.py
 ├ tran.dat
 ├ therm.dat
 ├ species-files
 |   └ *species*.lua
 └ OMEGA-tables
     ├ OMEGA11_table.txt
     └ OMEGA22_table.txt

----------------------------------------------
              INSTRUCTIONS:                  |
Create files name therm.dat and tran.dat
which contain Chemkin style thermo and transport
data, respectively. Define the source in the
section at the top of the code, otherwise
it will be unreferenced. It will look into 
the species-files subdirectory, either making 
a new species file or appending any missing 
information to the end. Any files which were
not listed in either the thermo or tran files
will be purged from the species-file directory
after completion, leaving only those of interest
behind.
----------------------------------------------
  
WARNING:
DOES NOT INCLUDE: Lewis Numbers, CEA Thermo coeffs, CEA Visc, or CEA Therm Cond values
ONLY CONFIGURED FOR: C, H and O atomic structures
  
EXAMPLE: therm.dat format - standard Chemkin format:
H2                TPIS78H   2   00   00   00G   200.000  3500.000   1000.00    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01 8.46810200E+03    4
OH                RUS 78O   1H   1   00   00G   200.000  3500.000  1000.000    1
 3.09288767E+00 5.48429716E-04 1.26505228E-07-8.79461556E-11 1.17412376E-14    2
 3.85865700E+03 4.47669610E+00 3.99201543E+00-2.40131752E-03 4.61793841E-06    3
-3.88113333E-09 1.36411470E-12 3.61508056E+03-1.03925458E-01 8.81310600E+03    4

EXAMPLE: tran.dat format - standard Chemkin format:
H2                 1    38.000     2.920     0.000     0.790   280.000          
H2O                2   572.400     2.605     1.844     0.000     4.000          

See Chemkin thermo and transport manuals for meaning of values



Last Change Date: 19/02/2019
Revised by Oliver: 15/04/2019
"""

'''
===============================================================================
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
===============================================================================
'''

import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.optimize as opt

'''
IMPORTANT: Define source of information here:
'''
#==============================================================================
reference = ''
#==============================================================================

#if reference == '': reference = 'Source not specified'

# Do you want to plot the data?
plot = False


'''
===============================================================================
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
===============================================================================
'''
def bilinear_interp(dat, T_k, d_k):
    '''
    Bilinear interpolation. Used for generating the OMEGA11 and OMEGA22 values,
    as described by the tables in 'Transport Properties of Polar Gases 1961' by
    Monchick and Mason, as used in the Chemkin Transport Theory Manual.
    Input:
        dat = Raw OMEGA22 table as presented by Monchick and Mason. | axb array
        T_k = Reduced temperature. | scalar
        d_k = Reduced dipole moment. | scalar
    Output:
        OMEGA = OMEGA (11 or 22) value for parameters given
    '''
    # Reduced temperature table index
    T_header = np.concatenate((np.arange(0.1, 1.0, 0.1), 
                           np.arange(1.0, 2.0, 0.2), 
                           np.arange(2.0, 4.0, 0.5), 
                           np.arange(4.0, 10.0, 1.0), 
                           np.arange(10.0, 20.0, 2.0), 
                           np.arange(20.0, 40.0, 5.0), 
                           np.array([40.0, 50.0, 75.0, 100.0])))
    
    # Reduced dipole moment table index
    d_header = np.array([-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.25,
                         0.5, 0.75, 1.0, 1.5, 2.0, 2.5])  

    # Nearest, lower header index to provided temperature and dipole
    T_i = 0; d_i = 0

    # If temperature, dipole is outside bounds, set to nearest in bounds, else
    # nearest lower value
    for i in range(len(T_header)-1):
        if T_k <= min(T_header): T_i = 0
        elif T_k >= max(T_header): T_i = len(T_header)
        elif T_k >= T_header[i] and T_k < T_header[i+1]: T_i = i
    
    for i in range(len(d_header)-1):
        if d_k <= min(d_header): d_i = 0
        elif d_k >= max(d_header): d_i = len(d_header)
        elif d_k >= d_header[i] and d_k < d_header[i+1]: d_i = i
    
    # Set up indices to do linear interpolation over. Originally setup this way
    # Separate from previous step so a quad/cubic interpolation could be used.
    # Accuracy increase is likely minor, and much more difficult to implement.
    T_ind = np.array([T_i, T_i+1]).astype(np.uint8)
    d_ind = np.array([d_i, d_i+1]).astype(np.uint8)
    
    # If relevant points are near bounds, shift so extrapolation is done
    if T_i == len(T_header)-1: T_ind -= 1
    if d_i == len(d_header)-1: d_ind -= 1
    if T_i == len(T_header): T_ind -= 2
    if d_i == len(d_header): d_ind -= 2
    
    # Get actual header data instead of indices
    T_lin = [T_header[T_ind[0]], T_header[T_ind[1]]]
    d_lin = [d_header[d_ind[0]], d_header[d_ind[1]]]
    
    # Get local relevant data bounding OMEGA value of interest. i.e.
    # * = point of interest
    #    dat00          dat01
    #                 
    #                *
    #    dat10          dat11
    
    dat_4 = [[ dat[ T_ind[0] ][ d_ind[0] ], dat[ T_ind[0] ][ d_ind[1] ] ],
             [ dat[ T_ind[1] ][ d_ind[0] ], dat[ T_ind[1] ][ d_ind[1] ] ]]
    
    # Do the bilinear interpolation
    OMEGA = (T_lin[1] - T_k)  / (T_lin[1] - T_lin[0]) *                    \
              ((d_lin[1] - d_k) / (d_lin[1] - d_lin[0]) * (dat_4[0][0]) +  \
               (d_k - d_lin[0]) / (d_lin[1] - d_lin[0]) * (dat_4[0][1])) + \
              (T_k - T_lin[0])  / (T_lin[1] - T_lin[0]) *                  \
              ((d_lin[1] - d_k) / (d_lin[1] - d_lin[0]) * (dat_4[1][0]) +  \
               (d_k - d_lin[0]) / (d_lin[1] - d_lin[0]) * (dat_4[1][1]))
    
    return OMEGA 
'''
===============================================================================
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
===============================================================================
'''
def CEA_curves(T_lower, T_upper, lower_coeffs, upper_coeffs):
    Cp_on_R = []; H_on_RT = []; S_on_R = []
    T = [T_lower, T_upper]; C = [lower_coeffs, upper_coeffs]
    for i in [0,1]:
        Cp_on_R += list(C[i][0] + C[i][1]*T[i] + C[i][2]*T[i]**2 + C[i][3]*T[i]**3 + C[i][4]*T[i]**4)
        H_on_RT += list(C[i][0] + C[i][1]*T[i]/2 + C[i][2]*T[i]**2/3 + C[i][3]*T[i]**3/4 + C[i][4]*T[i]**4/5 + C[i][5]/T[i])
        S_on_R  += list(C[i][0]*np.log(T[i]) + C[i][1]*T[i] + C[i][2]*T[i]**2/2 + C[i][3]*T[i]**3/3 + C[i][4]*T[i]**4/4 + C[i][6])  
    
    return np.array(Cp_on_R), np.array(H_on_RT), np.array(S_on_R)
'''
===============================================================================
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
===============================================================================
'''
def F(epsilon_on_kb, T): 
    '''
    Rotational relaxation function, as described in Chemkin manual.
    Input:
        epsilon_on_kb = Lennard-Jones Potential Well depth, from data file
        T             = Temperature value, K
    '''
    F_T = 1 + (np.pi**1.5 / 2)*(epsilon_on_kb/T)**0.5 + \
            (np.pi**2 / 4 + 2)*(epsilon_on_kb/T) + \
            (np.pi**1.5)*(epsilon_on_kb/T)**1.5
    return F_T

'''
===============================================================================
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
===============================================================================
'''
if not os.path.exists('species-files'):
    os.makedirs('species-files')
    
if not os.path.exists(os.path.join('OMEGA-tables','OMEGA11_table.txt')):
    raise ValueError('OMEGA11_table.txt is missing from directory OMEGA-tables.'+
                     ' Either find it or copy from Monchick and Mason.')
if not os.path.exists(os.path.join('OMEGA-tables','OMEGA22_table.txt')):
    raise ValueError('OMEGA22_table.txt is missing from directory OMEGA-tables.'+
                     ' Either find it or copy from Monchick and Mason.')

print('--- Reading in OMEGA tables')
# Read in OMEGA11 and OMEGA22 tables, source: Monchick and Mason
with open(os.path.join('OMEGA-tables','OMEGA11_table.txt'), 'rt') as OMEGA11_file:
    OMEGA11_table = OMEGA11_file.read().split('\n')[:-1]
    for i in range(len(OMEGA11_table)): OMEGA11_table[i] = OMEGA11_table[i].split(' ')
    OMEGA11_table = np.array([np.array(x, dtype=np.float64) for x in OMEGA11_table])
with open(os.path.join('OMEGA-tables','OMEGA22_table.txt'), 'rt') as OMEGA22_file:
    OMEGA22_table = OMEGA22_file.read().split('\n')[:-1]
    for i in range(len(OMEGA22_table)): OMEGA22_table[i] = OMEGA22_table[i].split(' ')
    OMEGA22_table = np.array([np.array(x, dtype=np.float64) for x in OMEGA22_table])

# Need to do checks for all needed files

# Atomic mass information, taken from periodic table, kg/mol
m_atom = {'H': 1.0079e-3,
          'C': 12.0107e-3,
          'O': 15.9994e-3,
          'N': 14.0067e-3,
          'AR': 39.948e-3,
          'HE': 4.002602e-3}

thermolines = []; transportlines = []

'''
Read Thermo file
'''
print('--- Reading thermo data')
try:
    with open('therm.dat', 'rt') as thermfile:
        thermolines = thermfile.read().split('\n')
except:
    raise ValueError('therm.dat does not exist in the current directory')
# Remove blank rows from thermo data, especially leading or trailing.
for i in reversed(range(len(thermolines))):
    if thermolines[i] == '': del thermolines[i]
    
if thermolines[-1][:9] == 'ENDOFDATA': del thermolines[-1] 
if thermolines[0][:6] == 'THERMO': del thermolines[:2]

# 4 thermo rows per species, if not than the data file is wrong.
if len(thermolines) % 4 != 0:
    raise ValueError('therm.dat does not have enough 4 rows per species, and may'+
                     ' not contain all valid data. Check your files, and remove headers.')     
    
# Get number of species in thermo data file
species_count = len(thermolines)//4 
thermo_species_list = []; thermo_data = {}

for i in range(species_count):
    thermo_species_list.append(thermolines[4*i].split(' ')[0])
    thermo_data[thermo_species_list[i]] = thermolines[4*i:4*(i+1)]

'''
Read Transport file
'''
print('--- Reading transport data')
try:
    with open('tran.dat', 'rt') as transport_file:
        transportlines = transport_file.read().split('\n')
except:
    raise ValueError('tran.dat does not exist in the current directory') 
        
for i in reversed(range(len(transportlines))):
    if transportlines[i] == '': del transportlines[i]

trans_species_list = []; trans_data = {}

for i in range(len(transportlines)):
    trans_species_list.append(transportlines[i].split(' ')[0])
    trans_data[trans_species_list[i]] = transportlines[i]

'''
Find any differences in species lists
'''
trans_set = set(trans_species_list)
thermo_set = set(thermo_species_list)

missing_from_thermo = trans_set - thermo_set
missing_from_trans = thermo_set - trans_set
common_species = thermo_set - missing_from_thermo - missing_from_trans

# Give species which will be used, missing species.
print(
'==============================================\n'+      
'Generating files for:          {}\n'.format(list(common_species))+
'Missing data from thermo.txt:  {}\n'.format(list(missing_from_thermo))+
'Missing data from trandat.txt: {}\n'.format(list(missing_from_trans))+
'==============================================\n')  

file_list = []

'''
===============================================================================
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
===============================================================================
'''

# Generate useful data from raw data
for species_name in common_species:
    '''
    Thermo Data
    '''
    
    print('--- Working on: {}'.format(species_name))
    
    # Get charge
    if species_name[-1] == '-': charge = -1
    elif species_name[-1] == '+': charge = 1
    else: charge = 0
    
    thermolines_i = thermo_data[species_name]
    header = thermolines_i[0] # Get header data
    atom_count = {'C': 0,
                  'H': 0,
                  'O': 0,
                  'AR': 0,
                  'N': 0,
                  'HE': 0} # Composition of atoms, from header
    
    for i in [24, 29, 34, 39]: # Locations of atom identifiers in header
        if header[i] != ' ' and header[i] != '0': # Check if identifier is valid i.e. nothing in >29
            atom_count[header[i:i+2].replace(' ','')] = int(header[i+4:i+5]) # Get atom count
    
    m_spec = 0 # molar mass of species, kg/mol
    
    # Get molecular mass using atom counts and atom weights
    for a in list(atom_count.keys()): m_spec += m_atom[a] * atom_count[a]
    
    lowertemp = float(header[48:54])
    uppertemp = float(header[57:63])
    try: midtemp = float(header[65:74])
    except: 
        midtemp = 1000
        print('------ WARNING: Middle temperature could not be read. Defaulting to 1000.0K')
    print('------ Temperature Range (K): {}, {}, {}'.format(lowertemp, midtemp, uppertemp))
    l_coeff = []; h_coeff = [] # Thermo constants for low and high temperatures
    
    for i in range(4): # Change strange exponent symbols
        thermolines_i[i] = thermolines_i[i].replace("d", "E").replace("D", "E")
    
    # Get data from 1st row
    for j in range(5):   h_coeff.append(float(thermolines_i[1][15*j:15*(j+1)].replace(" ", "")))
        
    # Get data from 2nd row
    for j in range(2):   h_coeff.append(float(thermolines_i[2][15*j:15*(j+1)].replace(" ", "")))
    for j in range(2,5): l_coeff.append(float(thermolines_i[2][15*j:15*(j+1)].replace(" ", "")))
        
    # Get data from 3rd row
    for j in range(4):   l_coeff.append(float(thermolines_i[3][15*j:15*(j+1)].replace(" ", "")))

    '''
    Transport Data
    '''
    transline_i = trans_data[species_name].split(' ')
    for i in reversed(range(len(transline_i))):
        if transline_i[i] == '': del transline_i[i]
    
    geometry      = transline_i[1] # 0, 1 or 2 -> mono, linear, non-linear
    epsilon_on_kb = float(transline_i[2])            # Kelvin
    
    sigma_raw     = transline_i[3]
    sigma         = float(sigma_raw)*1E-10           # Convert from Angstrom
    
    mu_raw        = transline_i[4]
    mu            = float(mu_raw)*3.1623E-25         # Convert from Debye
    
    alpha_raw     = transline_i[5]
    alpha         = float(alpha_raw)*(1E-10)**3      # Convert from Angstrom^3
    
    Z_rot_298     = float(transline_i[6])            # Unitless
    
    T_disc = 4 # Temperature discretisation.
    T_l = np.arange(lowertemp,midtemp,T_disc) # Lower temperature range
    T_h = np.arange(midtemp,uppertemp+1,T_disc) # Upper Temperature range
    T = np.append(T_l, T_h) # Total Temperature Range
    
    # Get Cp_on_R, from thermo coefficients
#    Cp_on_R = np.append(l_coeff[0] + l_coeff[1]*T_l + l_coeff[2]*T_l**2 + 
#                        l_coeff[3]*T_l**3 + l_coeff[4]*T_l**4,
#                        h_coeff[0] + h_coeff[1]*T_h + h_coeff[2]*T_h**2 + 
#                        h_coeff[3]*T_h**3 + h_coeff[4]*T_h**4)
    
    Cp_on_R, H_on_RT, S_on_R = CEA_curves(T_l, T_h, l_coeff, h_coeff)
    
    Cp_on_R_300 = l_coeff[0] + l_coeff[1]*300 + l_coeff[2]*300**2 + \
                    l_coeff[3]*300**3 + l_coeff[4]*300**4
    
    # Define some necessary values for later calculation
    kb = 1.38064852E-23        # Boltzmann Constant  m^2 kg / s^2 K
    R = 8.3144598              # Universal Gas Constant    J/K.mol
    avo = 6.0221409e+23        # Avogadros number
    m = m_spec / avo           # Molecular mass of species, kg per molecule
    Cv = Cp_on_R*R - R         # Specific heat w/r volume   J/K.mol
    P = 1                      # Unit pressure. rho*D_kk mean pressures cancel
    rho = P*m_spec / (R*T)     # 'Density' of species
    

        
    # Reduced Temperature and Dipole Moment
    T_red = T/epsilon_on_kb
    d_red = 0.5 * mu**2 * kb / (epsilon_on_kb * sigma**3)
    
    print('------ Performing bilinear interpolation of OMEGA tables')
    # Perform linear interpolations of the tables to get OMEGA values
    OMEGA11 = np.zeros(len(T)); OMEGA22 = np.zeros(len(T))
    for i in range(len(T_red)):
        OMEGA11[i] = bilinear_interp(OMEGA11_table, T_red[i], d_red)
        OMEGA22[i] = bilinear_interp(OMEGA22_table, T_red[i], d_red)
    
    # Calculate the viscosity term from method given by the Chemkin transport
    # manual. viscosity = eta
    # This term is used for a least squares fit to get coefficients
    eta = (5/16) * np.sqrt(np.pi*m*kb*T) / (np.pi*sigma**2 * OMEGA22)
    
    # Get rotational relaxation terms
    Z_rot = Z_rot_298 * F(epsilon_on_kb, 298)/F(epsilon_on_kb, T)
    
    # Get self diffusion term, as described by the manual
    D_kk = (3/8)*np.sqrt(np.pi* kb**3 * T**3 /m)/(P*np.pi*sigma**2 *OMEGA11)

    # Check molecule geometry for components
    # trans, rot, vib = translational, rotational, vibrational
    if geometry == '0':
        Cv_trans = R*3/2
        Cv_rot = 0
        Cv_vib = 0
    elif geometry == '1':
        Cv_trans = R*3/2
        Cv_rot = R
        Cv_vib = Cv - 5/2*R
    elif geometry == '2':
        Cv_trans = R*3/2
        Cv_rot = R*3/2
        Cv_vib = (Cv - 3*R)
    else:
        raise ValueError('Invalid geometry value for species: {}'
                         .format(species_name)) 
        
    # Factors
    A = 5/2 - rho*D_kk/eta
    B = Z_rot + (2/np.pi) * ((5/3)*(Cv_rot/R) + rho*D_kk/eta)
    
    # Contributions
    f_trans = (5/2) * (1 - (2/np.pi) * (Cv_rot/Cv_trans) * (A/B))
    f_rot = (rho * D_kk / eta) * (1 + (2/np.pi) * (A/B))
    f_vib = (rho * D_kk) / eta
    
    # lam = lambda, thermal conductivity
    lam = (eta/m_spec)*(f_trans*Cv_trans + f_rot*Cv_rot + f_vib*Cv_vib)
    
    # gamma = ratio of specific heats at 300 K
    gamma = 1/(1-(1/Cp_on_R_300))

    '''
    Some demonstrative plots
    '''
    if plot == True:
        print('------ Plotting data')
        fig, [ax_1, ax_2, ax_3] = plt.subplots(figsize=(15,5), nrows=1, ncols=3)
        
        ax_1.plot(T, eta, color='b', ls='-', label='CHEMKIN Viscosity')
        ax_2.plot(T, lam, color='r', ls='-', label='CHEMKIN Conductivity')
        
        
        
        ax_3.plot(T, Cp_on_R, color='r', ls='-', label='Grimech Cp/R')
        ax_3.plot(T, H_on_RT, color='g', ls='-', label='Grimech H/RT')
        ax_3.plot(T, S_on_R, color='b', ls='-', label='Grimech S/R')
        
        ax_2.set_title('{} Comparison of Conductivities'.format(species_name))
        ax_2.set_ylabel(r'$\frac{W}{m . K}$', fontsize=14)
        
        ax_1.set_title('{} Comparison of Viscosities'.format(species_name))
        ax_1.set_ylabel(r'Pa S', fontsize=14)
        
        ax_3.set_title('{} Cp/R'.format(species_name))
        
        ax_1.set_xlabel('K'); ax_2.set_xlabel('K'); ax_3.set_xlabel('K')
    
    '''
    Least Squares fit to get Viscosity and Thermal Conductivity coefficients
    '''
    print('------ Performing curve fit for viscosity, thermal conductivity')
    # Curve fit polynomial
    def func(x, A, B, C, D):
        return A + B*(x) + C*(x)**2 + D*(x)**3

    eta_params, _ = opt.curve_fit(func, np.log(T), np.log(eta), p0=[1,1,1,1])
    lam_params, _ = opt.curve_fit(func, np.log(T), np.log(lam), p0=[1,1,1,1])
    
    if plot == True:
        A, B, C, D = eta_params
        fitted = np.exp(func(np.log(T), A, B, C, D))
        ax_1.plot(T, fitted, color='k', ls='--', lw=1, label='Curvefit Viscosity')
        
        A, B, C, D = lam_params
        fitted = np.exp(func(np.log(T), A, B, C, D))
        ax_2.plot(T, fitted, color='k', ls='--', lw=1, label='Curvefit Conductivity')
    
        ax_1.legend(); ax_2.legend(); ax_3.legend()

    '''
    Write Data to file
    '''
    
    if species_name[-1] == '*': 
        formatted_name = "['"+species_name+"']"
        file_name = species_name[:-1]+'_star'
    else:
        formatted_name = "."+species_name
        file_name = species_name
    
    formatted_name = formatted_name.replace('HE','He').replace('AR','Ar')
    file_name = file_name.replace('HE','He').replace('AR','Ar')
    
    
    print('------ Writing data to file: {}.lua\n'.format(file_name))
    
    base_txt = [
    "--Auto-Generated File" ,
    "db{} = {{}}".format(formatted_name) ,
    "db{}.atomicConstituents = {{C={},H={},O={},Ar={},He={},N={}}}"
        .format(formatted_name, atom_count['C'], atom_count['H'], 
                atom_count['O'],  atom_count['AR'], 
                atom_count['HE'], atom_count['N']),
    "db{}.charge = {}".format(formatted_name, charge) ,
    "db{}.M = {{".format(formatted_name) ,
    "    value = {:.6f},".format(m_spec) ,
    "    units = 'kg/mol'," ,
    "   description = 'molecular mass'," ,
    "   reference = 'Periodic table'" ,
    "}" ,
    "db{}.gamma = {{".format(formatted_name) ,
    "   value = {:.6f},".format(gamma) ,
    "   units = 'non-dimensional'," ,
    "   description = 'ratio of specific heats at 300.0K'," ,
    "   reference = 'evaluated using Cp/R from Chemkin-II coefficients'" ,
    "}" ,
    "db{}.sigma = {{".format(formatted_name) ,
    "   value = {},".format(sigma_raw) ,
    "   units = 'Angstrom'," ,
    "   description = 'Lennard-Jones potential distance'," ,
    "   reference = ' '" ,
    "}" ,
    "db{}.epsilon = {{".format(formatted_name) ,
    "   value = {:.6f},".format(epsilon_on_kb) ,
    "   units = 'K'," ,
    "   description = 'Lennard-Jones potential well depth.'," ,
    "   reference = ' '" ,
    "}"]
    grimech_txt = [
    "db{}.grimechThermoCoeffs = {{".format(formatted_name) ,
    "   notes = ' '," ,
    "   nsegments = 2, " ,
    "   segment0 ={" ,
    "      T_lower = {:.3f},".format(lowertemp) ,
    "      T_upper = {:.3f},".format(midtemp) ,
    "      coeffs = {" ,
    "          {:.12e},".format(0.0) ,
    "          {:.12e},".format(0.0) ,
    "          {:.12e},".format(l_coeff[0]) ,
    "          {:.12e},".format(l_coeff[1]) ,
    "          {:.12e},".format(l_coeff[2]) ,
    "          {:.12e},".format(l_coeff[3]) ,
    "          {:.12e},".format(l_coeff[4]) ,
    "          {:.12e},".format(l_coeff[5]) ,
    "          {:.12e},".format(l_coeff[6]) ,
    "      }" ,
    "   }," ,
    "   segment1 = {" ,
    "      T_lower = {:.3f},".format(midtemp) ,
    "      T_upper = {:.3f},".format(uppertemp) ,
    "      coeffs = {" ,
    "          {:.12e},".format(0.0) ,
    "          {:.12e},".format(0.0) ,
    "          {:.12e},".format(h_coeff[0]) ,
    "          {:.12e},".format(h_coeff[1]) ,
    "          {:.12e},".format(h_coeff[2]) ,
    "          {:.12e},".format(h_coeff[3]) ,
    "          {:.12e},".format(h_coeff[4]) ,
    "          {:.12e},".format(h_coeff[5]) ,
    "          {:.12e},".format(h_coeff[6]) ,
    "      }" ,
    "   }" ,
    "}"]
    chemkin_visc_txt = [
    "db{}.chemkinViscosity = {{".format(formatted_name) ,
    "   notes = 'Generated by species-generator.py'," ,
    "   nsegments = 1, " ,
    "   segment0 ={" ,
    "      T_lower = {:.3f},".format(lowertemp) ,
    "      T_upper = {:.3f},".format(uppertemp) ,
    "      A = {:.12e},".format(eta_params[0]) ,
    "      B = {:.12e},".format(eta_params[1]) ,
    "      C = {:.12e},".format(eta_params[2]) ,
    "      D = {:.12e},".format(eta_params[3]) ,
    "   }" ,
    "}"]
    chemkin_therm_txt = [
    "db{}.chemkinThermCond = {{".format(formatted_name) ,
    "   notes = 'Generated by species-generator.py'," ,
    "   nsegments = 1, " ,
    "   segment0 ={" ,
    "      T_lower = {:.3f},".format(lowertemp) ,
    "      T_upper = {:.3f},".format(uppertemp) ,
    "      A = {:.12e},".format(lam_params[0]) ,
    "      B = {:.12e},".format(lam_params[1]) ,
    "      C = {:.12e},".format(lam_params[2]) ,
    "      D = {:.12e},".format(lam_params[3]) ,
    "   }" ,
    "}"]
    
    if not os.path.exists(os.path.join('species-files','{}.lua'.format(file_name))):
        species_file = open(os.path.join('species-files','{}.lua').format(file_name), "w")
        txt = base_txt + grimech_txt + chemkin_visc_txt + chemkin_therm_txt
        for line in txt:
            species_file.write(line+'\n')
    else:
        existing_txt = open(os.path.join('species-files','{}.lua').format(file_name), "r").read()
        txt = ['']
        if existing_txt.find('atomicConstituents') == -1: txt += base_txt
        if existing_txt.find('grimechThermoCoeffs') == -1: txt += grimech_txt
        if existing_txt.find('chemkinViscosity') == -1: txt += chemkin_visc_txt
        if existing_txt.find('chemkinThermCond') == -1: txt += chemkin_therm_txt
            
        species_file = open(os.path.join('species-files','{}.lua').format(file_name), "a")
        for line in txt:
            species_file.write(line+'\n')

    species_file.close()
    file_list.append(file_name+'.lua')
    
file_all = os.listdir('species-files')
for file in file_all:
    if file not in file_list:
        os.remove(os.path.join('species-files',file))

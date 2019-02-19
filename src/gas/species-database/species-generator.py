# -*- coding: utf-8 -*-
"""
Title:     Chemkin-Species to Eilmer-Species 
Author:    Samuel Dillon
Date:      05/02/2019
"""
"""
----------------------------------------------
              INSTRUCTIONS:                  |
Simply change the species variable and       |
the 'numb' list to your liking. The numb     |
list contain the formula numbers. Then, press|
run and a file will be produced containing   |
the relevant information.                    |
----------------------------------------------
              Extra Notes:
+ If you want to add more elements go 
  to the molecular weights .py file 
  and add them.
+ Make sure in your working directory
  are the thermo files, transport files
  and the molecular weights file.
"""
"""
Last Change Date: 19/02/2019
"""
import numpy as np

import re
#import species_list
#Molecular Weights kg/mol
hydrogen = 1.0079e-3
carbon = 12.0107e-3
oxygen = 15.9994e-3

#            Change these only
# -------------------------------------------
species = 'C2H4'
numb = [2, 4, 0] #[carbon, hydrogen, oxygen]
# -------------------------------------------

mweight = float("%.5g" % ((numb[0]*carbon)+(numb[1]*\
                          hydrogen)+(numb[2]*oxygen)))
thermo= open("thermo.txt","rt")
thermolines = []
transportlines = []
with open ('thermo.txt', 'rt') as thermfile:
    for line in thermfile:
        thermolines.append(line)


        
i=0
while i < len(thermolines):
    data=None
    low_constants = None
    #row = thermolines[i].find(species)
    if species[-1] == "-":
        row = re.search(r'\b(%s)\B' % species, thermolines[i])
    else:
        row = re.search(r'\b(%s)\b' % species, thermolines[i])
    if row == None:
        i = i+1
    else:
        data = [0, 0, 0, 0]
        data[0] = thermolines[i]
        data[1] = thermolines[i+1]
        data[2] = thermolines[i+2]
        data[3] = thermolines[i+3]
        lowertemp = data[0][48]+data[0][49]+data[0][50]+data[0][51]+data[0][52]+data[0][53]+data[0][54]
        uppertemp = data[0][57]+data[0][58]+data[0][59]+data[0][60]+data[0][61]+data[0][62]+data[0][63]
        low_constants = []
        high_constants = []
        i=0
        while i < len(data[1]):
            if i+14 > len(data[1]):
                break
            ahigh = data[1][i]+data[1][i+1]+data[1][i+2]+data[1][i+3]+data[1][i+4]+data[1][i+5]+data[1][i+6] \
                +data[1][i+7] + data[1][i+8]+data[1][i+9]+data[1][i+10]+data[1][i+11]+data[1][i+12]+data[1][i+13]+ \
                data[1][i+14] 
            high_constants.append(ahigh)
            i+=15
        i=0
        while i < 30:
            ahigh = data[2][i]+data[2][i+1]+data[2][i+2]+data[2][i+3]+data[2][i+4]+data[2][i+5]+data[2][i+6] \
                +data[2][i+7] + data[2][i+8]+data[2][i+9]+data[2][i+10]+data[2][i+11]+data[2][i+12]+data[2][i+13]+ \
                data[2][i+14] 
            high_constants.append(ahigh)
            i+=15
        i=30
        while i < len(data[2]):
            if i+14 > len(data[2]):
                break
            alow = data[2][i]+data[2][i+1]+data[2][i+2]+data[2][i+3]+data[2][i+4]+data[2][i+5]+data[2][i+6] \
                +data[2][i+7] + data[2][i+8]+data[2][i+9]+data[2][i+10]+data[2][i+11]+data[2][i+12]+data[2][i+13]+ \
                data[2][i+14] 
            low_constants.append(alow)
            i+=15      
        i=0
        while i < len(data[3]):
            if i+14 > len(data[3]):
                break
            alow = data[3][i]+data[3][i+1]+data[3][i+2]+data[3][i+3]+data[3][i+4]+data[3][i+5]+data[3][i+6] \
                +data[3][i+7] + data[3][i+8]+data[3][i+9]+data[3][i+10]+data[3][i+11]+data[3][i+12]+data[3][i+13]+ \
                data[3][i+14] 
            low_constants.append(alow)
            i+=15
        break
    
if data == None:
    raise ValueError('Could not Match this species')

with open ('trandat.txt', 'rt') as transport_file:
    for line in transport_file:
        transportlines.append(line)
i=0
while i < len(transportlines):
    data=None
    if species[-1] == "-":
        row1 = re.search(r'\b(%s)\B' % species, transportlines[i])
    else:
        row1 = re.search(r'\b(%s)\b' % species, transportlines[i])
    if row1 == None:
        i = i+1
    else:
        data = transportlines[i]
        sigma = data[34]+data[35] + data[36]+data[37] + data[38] + data[39] + data[40]
        epsilon = data[23]+data[24]+data[25] + data[26]+data[27] + data[28] + data[29] + data[30]
        break
print(species[-1])
temp_const = np.array([1, 300, 90000, 27000000, 8100000000, 0, 0])
CponR = 0.0
maxhigh = len(high_constants)
maxlow = len(low_constants)
for i in range(0, 6):
    low_constants[i] = low_constants[i].replace(" ", "") # Removes whitespace 
    low_constants[i] = low_constants[i].replace("d", "E")
    sumterm = float(low_constants[i])*temp_const[i]
    CponR = CponR + sumterm 

for i in range(0, maxlow):
    low_constants[i] = low_constants[i].replace(" ", "") # Removes whitespace
    if low_constants[i] != "":
        low_constants[i] = low_constants[i].replace("d", "E")
        low_constants[i] = float(low_constants[i])  

for i in range(0, maxhigh):
    high_constants[i] = high_constants[i].replace(" ", "") # Removes whitespace 
    if high_constants[i] != "":
        high_constants[i] = high_constants[i].replace("d", "E")
        high_constants[i] = float(high_constants[i])

gamma = float("%.5g" % (1-(1/CponR))**(-1))

  
    

values = (species, species, numb[0], numb[1], numb[2], \
          species, \
          species,\
          mweight,\
          species,\
          gamma,\
          species,\
          sigma,\
          species,\
          epsilon,\
          species,\
          float(lowertemp),\
          0.0, \
          0.0, \
          low_constants[0],\
          low_constants[1],\
          low_constants[2],\
          low_constants[3],\
          low_constants[4],\
          low_constants[5],\
          low_constants[6],\
          float(uppertemp),\
          0.0,\
          0.0,\
          high_constants[0],\
          high_constants[1],\
          high_constants[2],\
          high_constants[3],\
          high_constants[4],\
          high_constants[5],\
          high_constants[6],\
          )

# ---------------------------------------------------------------- 
txt = "--Auto-Generated File\n\
db.%s = {}\
\ndb.%s.atomicConstituents = {C=%s,H=%s,O=%s,}\
\ndb.%s.charge = 0\
\ndb.%s.M = {\
   \n\t value = %s,\
   \n\t units = 'kg/mol',\
   \n\tdescription = 'molecular mass',\
   \n\treference = 'Periodic table'\
\n\t}\
\ndb.%s.gamma = {\
   \n\tvalue = %s,\
   \n\tunits = 'non-dimensional',\
   \n\tdescription = 'ratio of specific heats at 300.0K',\
   \n\treference = 'evaluated using Cp/R from Chemkin-II coefficients'\
\n\t}\
\ndb.%s.sigma = {\
   \n\tvalue = %s,\
   \n\tunits = 'Angstrom',\
   \n\tdescription = 'Lennard-Jones potential distance',\
   \n\treference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'\
\n\t}\
\ndb.%s.epsilon = {\
   \n\tvalue = %s,\
   \n\tunits = 'K',\
   \n\tdescription = 'Lennard-Jones potential well depth.',\
   \n\treference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'\
\n\t}\
\ndb.%s.grimechThermoCoeffs = {\
   \n\tnotes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',\
   \n\tnsegments = 2, \
   \n\tsegment0 ={\
      \n\tT_lower = %.1f,\
      \n\tT_upper = 1000.0,\
      \n\tcoeffs = {\
         \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
      \n\t\t}\
   \n\t},\
   \n\tsegment1 = {\
      \n\tT_lower = 1000.0,\
      \n\tT_upper = %.1f,\
      \n\tcoeffs = {\
         \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
      \n\t\t}\
   \n\t}\
\n\
}\
\n\n" % values
                  
# ----------------------------------------------------------------                  
txt1 = "--Auto-Generated File\n\
db[\"%s\"] = {}\ndb[\"%s\"].atomicConstituents = {C=%s,H=%s,O=%s,}\
\ndb[\"%s\"].charge = 0\
\ndb[\"%s\"].M = {\
   \n\tvalue = %s,\
   \n\tunits = 'kg/mol',\
   \n\tdescription = 'molecular mass',\
   \n\treference = 'Periodic table'\
\n\t}\
\ndb[\"%s\"].gamma = {\
   \n\tvalue = %s,\
   \n\tunits = 'non-dimensional',\
   \n\tdescription = 'ratio of specific heats at 300.0K',\
   \n\treference = 'evaluated using Cp/R from Chemkin-II coefficients'\
\n\t}\
\ndb[\"%s\"].sigma = {\
   \n\tvalue = %s,\
   \n\tunits = 'Angstrom',\
   \n\tdescription = 'Lennard-Jones potential distance',\
   \n\treference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'\
\n\t}\
\ndb[\"%s\"].epsilon = {\
   \n\tvalue = %s,\
   \n\tunits = 'K',\
   \n\tdescription = 'Lennard-Jones potential well depth.',\
   \n\treference = ' A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang)- transport file.'\
\n\t}\
\ndb[\"%s\"].grimechThermoCoeffs = {\
   \n\tnotes = 'data from A DETAILED KINETIC MODEL OF C2- AND C3- FUEL COMBUSTION (Hai Wang) thermo file.',\
   \n\tnsegments = 2, \
   \n\tsegment0 ={\
      \n\tT_lower = %.1f,\
      \n\tT_upper = 1000.0,\
      \n\tcoeffs = {\
         \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
      \n\t\t}\
   \n\t},\
   \n\tsegment1 = {\
      \n\tT_lower = 1000.0,\
      \n\tT_upper = %.1f,\
      \n\tcoeffs = {\
         \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
         \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
          \n\t\t% 14.12e,\
      \n\t\t}\
   \n\t}\
\n\
}\
\n\n" % values
# ---------------------------------------------------------------- 
                  
species_file = open('%s.lua' % species, "w")
if species.find('-') <= 0:
    species_file.write(txt)
else: 
    species_file.write(txt1)
    
species_file.close()

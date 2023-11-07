"""
A program to combine PITOT3 gas models into new ones...

This is just a script for now, but I figured that I should get it in the repository...

Chris James (c.james4@uq.edu.au) - 07/11/23

"""

VERSION_STRING = '07-Nov-2023'

from pitot3_utils.pitot3_classes import eilmer4_CEAGas_input_file_reader, eilmer4_CEAGas_input_file_creator
from datetime import datetime

gas_model_1 = 'giant-planet-h2-80-he'
gas_model_2 = 'air23species'

gas_model_filename_1 = f'cea-{gas_model_1}-gas-model.lua'
gas_model_filename_2 =  f'cea-{gas_model_2}-gas-model.lua'

# this is by moles
percentage_gas_model_1 = 0.9
percentage_gas_model_2 = 1 - percentage_gas_model_1

combined_mixtureName = f'{percentage_gas_model_1*100:.0f}%-{gas_model_1}-{percentage_gas_model_2*100:.0f}%-{gas_model_2}'

output_filename = f'cea-{combined_mixtureName}-gas-model.lua'

print("The new gas model will be called:")
print(output_filename)

# now we open the two files and pull stuff out of them to combine them

mixtureName_1, speciesList_1, reactants_1, inputUnits_1, withIons_1, trace_1 = eilmer4_CEAGas_input_file_reader(gas_model_filename_1)

mixtureName_2, speciesList_2, reactants_2, inputUnits_2, withIons_2, trace_2 = eilmer4_CEAGas_input_file_reader(gas_model_filename_2)

# start with the species from the first gas model and then add any new species which are not in the species list for the first one

speciesList = speciesList_1

for species in speciesList_2:
    if species not in speciesList:
        speciesList.append(species)

# we add the reactants for each gas multiplied by their factor
# remembering that we need to take into account that the same gases might be in each gas model...

reactants = {}

for reactant in reactants_1:
    reactants[reactant] = reactants_1[reactant]*percentage_gas_model_1

for reactant in reactants_2:
    # first we need to check that it isn't already in the dictionary...
    if reactant not in reactants:
        reactants[reactant] = reactants_2[reactant] * percentage_gas_model_2
    else:
        #we add it to the existing value...
        reactants[reactant] = reactants[reactant] + reactants_2[reactant] * percentage_gas_model_2

#print(reactants)

if inputUnits_1 != inputUnits_2:
    print("Both of the gas models do not have the same inputUnits. This is not going to work...")
    raise Exception
else:
    # just set it to the first one if they're both the same
    inputUnits = inputUnits_1

if withIons_1 or withIons_2:
    print("At least one of the gas models have ions, so ions will be used in the calculation.")
    withIons = True
else:
    withIons = False

trace = min(trace_1, trace_2)

print(f"Using a trace value of {trace} as it was the smaller of the two values.")

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

header = f"""-- CEA Gas model made automatically by PITOT3 on the {dt_string}
-- This gas model is a combination of two gas models:
-- {percentage_gas_model_1*100:.0f}% of {gas_model_filename_1}
-- {percentage_gas_model_2*100:.0f}% of {gas_model_filename_2}"""

eilmer4_CEAGas_input_file_creator(output_filename, combined_mixtureName, speciesList, reactants,
                                      inputUnits, withIons, trace, header)








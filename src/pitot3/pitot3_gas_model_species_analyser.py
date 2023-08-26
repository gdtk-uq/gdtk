"""
This is a program that goes through the folder with all of the PITOT3 preset gas models and makes a list of
all of the species used with CEAGAS models in the code which we will later user to get their molecular weights.

Chris James (c.james4@uq.edu.au) - 23/11/22

"""

import os
from datetime import date

pitot3_preset_gas_models_folder = os.path.expanduser("~") + '/gdtkinst/share/pitot3_data/preset_gas_models'

CEA_thermo_inp_location = os.path.expanduser("~") + '/gdtkinst/share/cea-cases/thermo.inp'

species_output_filename = os.path.expanduser("~") + '/gdtkinst/share/pitot3_data/PITOT3_species_molecular_weights.yaml'

total_species_list = []

with os.scandir(pitot3_preset_gas_models_folder) as folder_contents:

    for file in folder_contents: #technically there could be folders too, but don't worry about that...
        # first check that it is a .lua file

        filename = file.name

        if '.lua' in filename: # open the file and we work through it
            print(filename)

            # we are after two lines, one with 'model' in it which also has "CEAGas" in it and one with speciesList in it.
            # so we loop through looking for those...

            is_a_CEAGas_gas_model = False # start by assuming it is false...

            with open(f"{pitot3_preset_gas_models_folder}/{filename}") as file:
                for line in file:
                    if 'model' in line and 'CEAGas' in line: # it is a CEAGas gas model...
                        is_a_CEAGas_gas_model = True
                    elif 'speciesList' in line: # this is the line with the species...
                        split_line = line.strip().split('=') # split it into before and after the variable name, after the variable name will be the specieslist

                        raw_species_list_string = split_line[1].strip()

                        # we need to split the input into a list...
                        split_species_list = raw_species_list_string.split(',')

                        speciesList = []

                        characters_to_replace = ['{', '}', "'", '"',' ']

                        for species in split_species_list:
                            for character in characters_to_replace:
                                species = species.replace(character, '')

                            if species: # to remove empty values
                                speciesList.append(species)

                        print(speciesList)

                if is_a_CEAGas_gas_model and speciesList:

                    # add any species which isn't already in the species list

                    for species in speciesList:
                        if species not in total_species_list:
                            total_species_list.append(species)

# now sort out final species list and then print it...
total_species_list.sort()

print(f"There are {len(total_species_list)} total species in the code.")
print("The species are:")
print(total_species_list)

# now we need to loop through and get their molecular weights from the CEA thermo.inp file...

species_MW_dict = {}

# this is a bit brute force but we will run through the whole CEA thermo.inp file and check if every
# particle is in our total species list and if so, we grab its MW and store it.

with open(CEA_thermo_inp_location) as CEA_thermo_inp_file:

    have_found_species = False  # use this to find species

    for line in CEA_thermo_inp_file:

        if line[-1] != '!' and 'thermo' not in line: # ! seems to be the comment character and thermo is the first line
            # if the line doesn't begin with a space, it is an element and it the first character will be a letter too... (we check both)
            if line [0] != ' ' and line[0].isalpha():
                species_name = line.split(' ')[0]

                #print(species_name)
                have_found_species = True # in this case the MW will be on the next line...
            elif have_found_species: # we found a species on the line before and now want is info
                # the MW ENDS at index 65 and has a changing length, so we grab the line up to that value, and then it will be the last value
                #print(line[0:65])

                split_line = line[0:65].split(' ')
                molecular_weight = float(split_line[-1])

                #print(molecular_weight)
                have_found_species = False

                if species_name in total_species_list:
                    species_MW_dict[species_name] = molecular_weight

print('-'*60)
print("The list of species found in CEA are:")
print(total_species_list)

print("Their molecular weights are:")
print(species_MW_dict)

print(f"We have found {len(species_MW_dict)} total species in CEA.")

print('-'*60)
if len(total_species_list) == len(species_MW_dict.keys()):
    print("The anount of species in the code and the amount of species found from CEA agree.")
else:
    print("WARNING: The anount of species in the code and the amount of species found from CEA DO NOT agree. There may be an issue")

# now output to file...

with open(species_output_filename, mode = 'w') as species_output_file:
    header = f'# PITOT3 species molecular weight list created {date.today().strftime("%d-%b-%Y")}. Units are g/mol'

    species_output_file.write(header + '\n')

    for species in species_MW_dict:

        molecular_weight = species_MW_dict[species]

        output_line = f'{species} : {molecular_weight}'

        species_output_file.write(output_line + '\n')

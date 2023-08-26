"""
pitot3_classes.py

Just a file to store all of the classes for PITOT3,
to make the main file a bit cleaner.

Chris James (c.james4@uq.edu.au) 01/01/21

"""

import sys, os, math
import yaml

from gdtk.gas import GasModel, GasState, GasFlow
from gdtk.ideal_gas_flow import p0_p
from gdtk.numeric.zero_solvers import secant

def eilmer4_CEAGas_input_file_creator(output_filename, mixtureName, speciesList, reactants,
                                      inputUnits, withIons, trace = 1.0e-6,
                                      header = '-- CEA Gas model made automatically by PITOT3'):
    """Just a function to make an input file for the CEAGas object..."""


    if '.lua' not in output_filename:
        output_filename += '.lua'

    with open(output_filename, 'w') as gas_file:
        gas_file.write(header + '\n')
        gas_file.write('\n')

        gas_file.write('model = "CEAGas"' + '\n')
        gas_file.write('\n')

        gas_file.write('CEAGas = {' + '\n')
        gas_file.write('  mixtureName = "{0}",'.format(mixtureName) + '\n')
        gas_file.write('  speciesList = {0},'.format(speciesList).replace('[','{').replace(']','}') + '\n')
        gas_file.write('  reactants = {0},'.format(reactants).replace(':',' =').replace("'",'') + '\n')
        gas_file.write("  inputUnits = '{0}',".format(inputUnits) + '\n')
        if withIons:
            gas_file.write('  withIons = true,' + '\n')
        else:
            gas_file.write('  withIons = false,' + '\n')
        gas_file.write('  trace = {0}'.format(trace) + '\n')
        gas_file.write('}')

    return

def eilmer4_CEAGas_gmodel_without_ions_creator(gmodel_filename):

    """
    This is a function which takes a link to a CEA Gas file with ions and makes a version of it without ions
    in the current folder. It then returns the filename so the gmodel can be used in the program.

    This function is needed as sometimes PITOT3 fails in expansions at low temperatures where ions asre not needed, anyway,
    this gmodel without ions will hopefully allow these cases to be solved.

    :param gmodel_filename:
    :return:
    """

    # first we load the file line by line and pull out what we need
    # we are making the assumption that what we have is a proper good quality functioning CEA Gas gas model...

    variables_to_look_for_list = ['mixtureName', 'speciesList', 'reactants', 'inputUnits', 'trace']

    # model = "CEAGas"
    #
    # CEAGas = {
    #   mixtureName = "driver_gas",
    #   speciesList = {'He', 'Ar'},
    #   reactants = {He = 0.8, Ar = 0.2},
    #   inputUnits = 'moles',
    #   withIons = false,
    #   trace = 1e-06
    # }

    with open(gmodel_filename) as gmodel_file:
        for line in gmodel_file:

            if line[0:2] != '--' and '=' in line: # -- is the lua comment character and any line we need will have = in it...

                split_line = line.strip().split('=')

                variable_name = split_line[0].strip()
                variable = split_line[1].strip()

                if variable_name in variables_to_look_for_list:
                    # we just have a case for each variable as that seemed easiest to me...
                    if variable_name == 'mixtureName':
                        # we add without ions to the end of it

                        characters_to_replace = ["'", '"', ',']

                        for character in characters_to_replace:
                            variable = variable.replace(character, '')

                        mixtureName = variable + '-without-ions'
                    if variable_name == 'speciesList':
                        # we need to split the input into a list...
                        split_species_list = variable.split(',')

                        speciesList = []

                        characters_to_replace = ['{', '}', "'", '"']

                        for species in split_species_list:
                            for character in characters_to_replace:
                                species = species.replace(character, '')

                            # any species with a + or a - in it is a charged particle, remove them!

                            if '+' not in species and '-' not in species and species:
                                speciesList.append(species)

                    if variable_name == 'reactants':
                        # this will already be kind of split for us by the original = split... so let us recombine that
                        # and then keep going

                        variable = '='.join(split_line[1:])

                        split_reactants_list = variable.split(',')

                        characters_to_replace = ['{', '}', "'", '"', ' ', ',']

                        reactants = {}

                        for reactant in split_reactants_list:
                            for character in characters_to_replace:
                                reactant = reactant.replace(character, '')

                            split_reactant = reactant.split('=')

                            if len(split_reactant) == 2:
                                reactants[split_reactant[0]] = float(split_reactant[1])

                    if variable_name == 'inputUnits':

                        characters_to_replace = ["'", '"', ',']

                        for character in characters_to_replace:
                            variable = variable.replace(character, '')

                        inputUnits = variable

                    if variable_name == 'trace':

                        trace = float(variable)

    # now the final variable is withIons which will be False
    withIons = False

    # now we just need to give ourselves a filename for the new gmodel
    # and make the file...

    gmodel_justfilename = gmodel_filename.split('/')[-1] # remove the folders etc.

    # the pitot3 format is kind of cea-gmodel-name-gas-model.lua
    # assuming that lets remove anything like that...
    parts_to_remove = ['cea-', '-gas-model', '.lua']

    gmodel_name = gmodel_justfilename

    for part in parts_to_remove:
        gmodel_name = gmodel_name.replace(part, '')

    # now we reconstitute it in the PITOT3 style...

    gmodel_without_ions_filename = f'cea-{gmodel_name}-without-ions-gas-model.lua'


    header = '-- modified gas model made with PITOT3 without ions for low temperature operation'

    eilmer4_CEAGas_input_file_creator(gmodel_without_ions_filename, mixtureName, speciesList, reactants, inputUnits,
                                      withIons, trace, header)

    return gmodel_without_ions_filename

def eilmer4_IdealGas_gas_model_creator(gas_state):
    """
    This function is to make a dummy ideal gas gas model for PITOT3 to do the frozen normal shock calculation
    over the test model with.

    This gas model should not be used for other things! As it will have dummy entropy, viscocity and thermal conductivity
    values which are just based on air.

    :param gas_state:
    :return:
    """

    ideal_gas_gmodel_filename = 'PITOT3_ideal_gas_test_section_gmodel.lua'

    # we need the specific heat ratio (gamma) and molecular mass
    # from the CEA gas model to put into the ideal gas gas model.

    gamma = gas_state.gamma
    molecular_mass = gas_state.ceaSavedData['Mmass']

    with open(ideal_gas_gmodel_filename, mode = 'w') as ideal_gas_gmodel_file:

        gas_model_string = \
f"""-- this is a dummy ideal gas gas model made by PITOT3 purely 
-- to perform a frozen shock over the test model.
-- Do not use it for anything else! As the entropy, viscosity and 
-- thermal conductivity values are just dummy values for air.

model = "IdealGas"

IdealGas = {{
    speciesName = 'dummy PITOT3 test section species',
    mMass = {molecular_mass},
    gamma = {gamma},
    entropyRefValues = {{
        s1 = 0.0,
        T1 = 298.15,
        p1 = 101.325e3
    }},
    viscosity = {{
    model = 'Sutherland',
    mu_ref = 1.716e-5,
    T_ref = 273.0,
    S = 111.0,
    }},
    thermCondModel = {{
    model = 'Sutherland',
    T_ref = 273.0,
    k_ref = 0.0241,
    S = 194.0
    }}
}}

"""

        ideal_gas_gmodel_file.write(gas_model_string)

    return ideal_gas_gmodel_filename


def finite_wave_dp_wrapper(state1, v1, characteristic, p2, state2, gas_flow, steps=100,
                           gmodel_without_ions = None, cutoff_temp_for_no_ions = 5000.0, number_of_calculations = 10):

    """
    Similar to how in PITOT3 I had a "normal shock wrapper" to deal with cases where the normal shock just didn't work
    this is my way of trying to do this in PITOT3. It breaks the finite wave dp up into a set of blocks, so if any of them
    bail out, we can try to turn ions off and keep going (if we're below the characteristic temp...)


    :param state1:
    :param v1:
    :param characteristic:
    :param p2:
    :param state2:
    :param steps:
    :param gmodel_without_ions:
    :param cutoff_temp:
    :return:
    """

    import numpy as np

    p1 = state1.p

    pressure_list = np.geomspace(p1, p2, number_of_calculations + 1)

    substeps = int(steps/number_of_calculations)

    # let us make a new local versions of states 1 and 2 for this local calculation...

    state1_local = GasState(state1.gmodel)

    state1_local.p = state1.p
    state1_local.T = state1.T

    # if we're copying a thermally perfect gas, we also need to port over the mass fractions, so we'll do that here...
    if state1.gmodel.type_str == 'ThermallyPerfectGas':
        state1_local.massf = state1.massf

    state1_local.update_thermo_from_pT()
    state1_local.update_sound_speed()

    state2_local = GasState(state1.gmodel)

    # if we're copying a thermally perfect gas, we also need to port over the mass fractions, so we'll do that here...
    if state1.gmodel.type_str == 'ThermallyPerfectGas':
        state2_local.massf = state1.massf

    # now we loop through all of the pressures
    for pressure in pressure_list:

        if pressure != p1: # skip that first value...

            try:

                v2g = gas_flow.finite_wave_dp(state1_local, v1, characteristic, pressure, state2_local, steps=substeps)

                have_performed_calculation_without_ions = False # didn't need it

            except Exception as e:
                print(e)
                print("Unsteady expansion calculation failed, probably due to a CEA error.")
                
                if state1_local.T < cutoff_temp_for_no_ions and gmodel_without_ions:
                    print(f"We are below the set cutoff temperature of {cutoff_temp_for_no_ions} K so we are going to try performing this part of the calculation without ions.")

                    state1_local_without_ions = GasState(gmodel_without_ions)

                    state2_local_without_ions = GasState(gmodel_without_ions)

                    state1_local_without_ions.p = state1_local.p
                    state1_local_without_ions.T = state1_local.T

                    state1_local_without_ions.update_thermo_from_pT()
                    state1_local_without_ions.update_sound_speed()

                    gas_flow_without_ions = GasFlow(gmodel_without_ions)

                    v2g = gas_flow_without_ions.finite_wave_dp(state1_local_without_ions, v1, characteristic, pressure,
                                                               state2_local_without_ions, steps=substeps)

                    have_performed_calculation_without_ions = True

                    print("Calculation was successful. Turning ions back on.")

                    # turning ions back on for next calculation...

                    # set the gas state here without ions, in case it fails here
                    state2_local.gmodel = gmodel_without_ions

                    state2_local.p = state2_local_without_ions.p
                    state2_local.T = state2_local_without_ions.T

                    state2_local.update_thermo_from_pT()
                    state2_local.update_sound_speed()

                    # then turn the ions back on here...

                    state2_local.gmodel = state1.gmodel
                else:
                    raise Exception("pitot3_classes.finite_wave_dp_wrapper: Unsteady expansion failed but the temperature is too high to justify turning ions off.")

            # now we need to make our state2 state 1 for the next step and v2g v1 as well

            v1 = v2g

            if have_performed_calculation_without_ions:
                # we turn ions off when we set the state here in case it would fail otherwise...

                state1_local.gmodel = gmodel_without_ions

            state1_local.p = state2_local.p
            state1_local.T = state2_local.T

            state1_local.update_thermo_from_pT()
            state1_local.update_sound_speed()

            if have_performed_calculation_without_ions:
                state1_local.gmodel = state1.gmodel

    # we may even need to turn off ions here too...

    # all of these statements will be from the last iteration here, so that is fine...
    if have_performed_calculation_without_ions:
            # we turn ions off when we set the state here in case it would fail otherwise...

            state2.gmodel = gmodel_without_ions

    state2.p = state2_local.p
    state2.T = state2_local.T

    # if we're copying a thermally perfect gas, we also need to port over the mass fractions, so we'll do that here...
    if state2_local.gmodel.type_str == 'ThermallyPerfectGas':
        state2.massf = state2_local.massf

    state2.update_thermo_from_pT()
    state2.update_sound_speed()

    if have_performed_calculation_without_ions:
        state2.gmodel = state1.gmodel

    return v2g

def finite_wave_dv_wrapper(state1, v1, characteristic, v2_target, state2, gas_flow, steps=100, t_min=200.0,
                           gmodel_without_ions=None, cutoff_temp_for_no_ions=5000.0, number_of_calculations=10):
    """
    Similar to how in PITOT3 I had a "normal shock wrapper" to deal with cases where the normal shock just didn't work
    this is my way of trying to do this in PITOT3. It breaks the finite wave dv up into a set of blocks, so if any of them
    bail out, we can try to turn ions off and keep going (if we're below the characteristic temp...)


    :param state1:
    :param v1:
    :param characteristic:
    :param v2_target:
    :param state2:
    :param steps:
    :param gmodel_without_ions:
    :param cutoff_temp:
    :return:
    """

    import numpy as np

    # I used geomspace above as the pressure changes logarthimically, velocity is more linear so I have just used linspace here...
    velocity_list = np.linspace(v1, v2_target, number_of_calculations + 1)

    substeps = int(steps / number_of_calculations)

    # let us make a new local versions of states 1 and 2 for this local calculation...

    state1_local = GasState(state1.gmodel)

    state1_local.p = state1.p
    state1_local.T = state1.T

    state1_local.update_thermo_from_pT()
    state1_local.update_sound_speed()

    state2_local = GasState(state1.gmodel)

    # now we loop through all of the pressures
    for velocity in velocity_list:
        if velocity != v1:  # skip that first value...

            try:

                v2g = gas_flow.finite_wave_dv(state1_local, v1, characteristic, velocity, state2_local, steps=substeps, t_min= t_min)

                have_performed_calculation_without_ions = False # didn't need it

            except Exception as e:
                print(e)
                print("Unsteady expansion calculation failed, probably due to a CEA error.")

                if state1_local.T < cutoff_temp_for_no_ions and gmodel_without_ions:
                    print(
                        f"We are below the set cutoff temperature of {cutoff_temp_for_no_ions} K so we are going to try performing this part of the calculation without ions.")

                    state1_local_without_ions = GasState(gmodel_without_ions)

                    state2_local_without_ions = GasState(gmodel_without_ions)

                    state1_local_without_ions.p = state1_local.p
                    state1_local_without_ions.T = state1_local.T

                    state1_local_without_ions.update_thermo_from_pT()
                    state1_local_without_ions.update_sound_speed()

                    gas_flow_without_ions = GasFlow(gmodel_without_ions)

                    v2g = gas_flow_without_ions.finite_wave_dv(state1_local_without_ions, v1, characteristic, velocity,
                                                               state2_local_without_ions, steps=substeps, t_min=t_min)

                    have_performed_calculation_without_ions = True

                    print("Calculation was successful. Turning ions back on")

                    # turning ions back on for next calculation...

                    # set the gas state here without ions, in case it fails here
                    state2_local.gmodel = gmodel_without_ions

                    state2_local.p = state2_local_without_ions.p
                    state2_local.T = state2_local_without_ions.T

                    state2_local.update_thermo_from_pT()
                    state2_local.update_sound_speed()

                    # then turn the ions back on here...

                    state2_local.gmodel = state1.gmodel

                else:
                    raise Exception(
                        "pitot3_classes.finite_wave_dp_wrapper: Unsteady expansion failed but the temperature is too high to justify turning ions off.")

            # now we need to make our state2 state 1 for the next step and v2g v1 as well

            v1 = v2g

            if have_performed_calculation_without_ions:
                # we turn ions off when we set the state here in case it would fail otherwise...

                state1_local.gmodel = gmodel_without_ions

            state1_local.p = state2_local.p
            state1_local.T = state2_local.T

            state1_local.update_thermo_from_pT()
            state1_local.update_sound_speed()

            if have_performed_calculation_without_ions:
                state1_local.gmodel = state1.gmodel

    # we may even need to turn off ions here too...

    # all of these statements will be from the last iteration here, so that is fine...
    if have_performed_calculation_without_ions:
            # we turn ions off when we set the state here in case it would fail otherwise...

            state2.gmodel = gmodel_without_ions

    state2.p = state2_local.p
    state2.T = state2_local.T

    state2.update_thermo_from_pT()
    state2.update_sound_speed()

    if have_performed_calculation_without_ions:
        state2.gmodel = state1.gmodel

    return v2g

def pitot3_input_yaml_file_creator(config_dict, output_filename):
    """
    Function that takes a PITOT3 config dict and makes a related yaml file.

    I am building this for a condition builder, mainly, to get the simulations ready to go...

    :param config_dict:
    :return:
    """

    print(f'Exporting config dictionary to the file {output_filename}.')

    with open(output_filename, 'w') as output_file:
        # yaml.dump did not deal with the Nonetypes correctly, so I just did this manually here...

        output_file.write("# PITOT3 input file created automatically by the function pitot3_input_yaml_file_creator \n \n")

        for key in config_dict.keys():
            item = config_dict[key]

            # I don't think a proper string (with quotations) is needed for yaml files
            if isinstance(item, str):
                output_file.write(f"{key} : '{item}' \n")
            else:
                output_file.write(f"{key} : {item} \n")

    return

def pitot3_states_dict_json_output_file_creator(states_dict, output_filename):

    """
    Function that takes the PITOT3 result states dictionary and makes a dictionary output version of each state
    and outputs it to a json file.

    :param states_dictionary:
    :param output_filename:
    :return:
    """

    import json

    print('-'*60)
    print(f'Exporting facility states data to the json file {output_filename}.')
    print('-'*60)

    states_dict_output_dict = {} # bad name, I know...

    for state_name in states_dict:

        facility_state = states_dict[state_name]

        if state_name in ['s2','s7','s8','s10e']: # just calculate the transport coeffcients and unit Re for these important states... (as I had some issues with it sometimes...
            states_dict_output_dict[state_name] = facility_state.get_dictionary_output(add_trans_coeffs_and_unit_Re = True)
        else:
            states_dict_output_dict[state_name] = facility_state.get_dictionary_output()

    with open(output_filename, "w") as output_file:
        json.dump(states_dict_output_dict, output_file)

    return

def pitot3_states_dict_json_output_file_loader(json_output_filename):
    """

    Just a simple function to load the json output file back into Python.

    :param output_filename:
    :return:
    """

    import json

    with open(json_output_filename, "r") as json_output_file:
        states_dict_output_dict = json.load(json_output_file)

    return states_dict_output_dict

def pitot3_pickle_output_file_creator(dict_of_objects, output_filename):

    """
    A function to pickle the PITOT3 objects, mainly for the condition builder intermediate results.

    :param dict_of_objects:
    :param output_filename:
    :return:
    """

    import pickle

    with open(output_filename, "wb") as output_file:
        pickle.dump(dict_of_objects, output_file, pickle.HIGHEST_PROTOCOL)

    return

def pitot3_pickle_output_file_loader(pickle_output_filename):
    """

    Just a simple function to load the pickle output file back into Python.

    :param pickle_output_filename:
    :return:
    """

    import pickle

    with open(pickle_output_filename, "rb") as pickle_output_file:
        dict_of_objects = pickle.load(pickle_output_file)

    return dict_of_objects

def pitot3_single_line_output_file_creator(config_data, object_dict, states_dict, json_output_states_dict = None, output_to_file = True,
                                           condition_builder_output = False, return_title_and_result_lists = False):

    """
    This is a function to allow us to make the PITOT3 output a single line csv (with a lot of columns)
    handy for loading simply, sometimes, and will be useful for creating output for the condition builder.

    The json output states dict can be used for running this function instead of the states dict itself
    (This is handy when using pickled objects where the original gas model may not be available

    :param config_data:
    :param object_dict:
    :param states_dict:
    :return:
    """

    # we just make two lists which we print to file on two lines when we're done...
    title_list = []
    value_list = []

    # first we need to know the type of facility and whether we have a secondary driver or a nozzle
    # or not as that changes some decisions later on...

    facility_type = config_data['facility_type']
    if 'secondary_driver' in object_dict:
        secondary_driver_flag = True
    else:
        secondary_driver_flag = False
    if 'nozzle' in object_dict:
        nozzle_flag = True
    else:
        nozzle_flag = False

    # it will also help to know what the test section state is, so we should get that too...

    test_section = object_dict['test_section']
    freestream_state = test_section.get_entrance_state()
    freestream_state_name = freestream_state.get_state_name()

    # start by adding some stuff from the config dictionary...

    config_dict_initial_variables = []

    if condition_builder_output:
        # we add the test number to the start too...
        config_dict_initial_variables += ['test_number']

    config_dict_initial_variables += ['driver_condition']

    if nozzle_flag:
        config_dict_initial_variables += ['area_ratio']

    if secondary_driver_flag:
        config_dict_initial_variables += ['secondary_driver_gas_gas_model', 'secondary_driver_gas_name']

    config_dict_initial_variables += ['test_gas_gas_model', 'test_gas_name']

    if facility_type == 'expansion_tube':
        config_dict_initial_variables += ['accelerator_gas_gas_model', 'accelerator_gas_name']

    if secondary_driver_flag: config_dict_initial_variables += ['psd1']

    config_dict_initial_variables += ['p1']

    if facility_type == 'expansion_tube': config_dict_initial_variables += ['p5']

    for variable in config_dict_initial_variables:
        title_list.append(variable)
        value_list.append(config_data[variable])

    # now add the shock speeds...

    if secondary_driver_flag:
        title_list.append('vsd1')
        value_list.append(object_dict['secondary_driver'].get_shock_speed())

    title_list.append('vs1')
    value_list.append(object_dict['shock_tube'].get_shock_speed())

    if facility_type == 'expansion_tube':
        title_list.append('vs2')
        value_list.append(object_dict['acceleration_tube'].get_shock_speed())

    # we start by just printing the test section values at the start, as they are useful to be able to find easily!

    freestream_variables_to_print_list = ['Ht', 'h', 'Ue', 'pitot_p','total_p']

    if states_dict:
        freestream_facility_state = states_dict[freestream_state_name]
        freestream_facilty_state_dictionary_output = freestream_facility_state.get_dictionary_output()

    else:
        freestream_facilty_state_dictionary_output = json_output_states_dict[freestream_state_name]

    for variable in freestream_variables_to_print_list:

        title_list.append(variable)
        value_list.append(freestream_facilty_state_dictionary_output[variable])

    # now we go through states and print some values, with some special values printed out when it recognises certain states...

    states_to_print_list = []

    if secondary_driver_flag:
        states_to_print_list += ['sd2', 'sd3']

    states_to_print_list += ['s2', 's3']

    if facility_type == 'expansion_tube':
        states_to_print_list += ['s7', 's6']

    if nozzle_flag:
        states_to_print_list += ['s8']

    states_to_print_list += ['s10f', 's10e']

    if 'cone_half_angle_degrees' in config_data:
        states_to_print_list += ['s10c']

    if 'wedge_angle_degrees' in config_data:
        states_to_print_list += ['s10w']

    generic_variables_to_print_list = ['p', 'T', 'rho', 'v', 'M', 'a']

    # dictionary to store when we want a variable to print more things...

    extra_variables_to_print_dict = {}
    extra_variables_to_print_dict['s2'] = ['gamma', 'R', 'Ht', 'h']
    extra_variables_to_print_dict[freestream_state_name] = ['gamma', 'R', 'unit_Re']
    extra_variables_to_print_dict['s10e'] = ['unit_Re']

    for state_name in states_to_print_list:

        state_name_without_the_s = state_name[1:]

        if state_name in extra_variables_to_print_dict:
            variables_to_print_list = generic_variables_to_print_list + extra_variables_to_print_dict[state_name]
        else:
            variables_to_print_list = generic_variables_to_print_list

        if 'unit_Re' in variables_to_print_list:
            add_trans_coeffs_and_unit_Re = True
        else:
            add_trans_coeffs_and_unit_Re = False

        # we need a special case here for if optional states (such as s10w) fail
        # we just set everything to None...
        if states_dict and state_name in states_dict:
            facility_state = states_dict[state_name]
            facilty_state_dictionary_output = facility_state.get_dictionary_output(add_trans_coeffs_and_unit_Re=add_trans_coeffs_and_unit_Re)

        elif json_output_states_dict and state_name in json_output_states_dict:
            facilty_state_dictionary_output = json_output_states_dict[state_name]
        else:
            # we just make a fake dictionary output for the state where everything is None
            facilty_state_dictionary_output = {}

            for variable in variables_to_print_list:
                facilty_state_dictionary_output[variable] = 'None'

        for variable in variables_to_print_list:
            if state_name not in ['sd2', 'sd3']:
                variable_name = f"{variable}{state_name_without_the_s}"
            else:
                variable_name = f"{variable}{state_name}"

            title_list.append(variable_name)
            value_list.append(facilty_state_dictionary_output[variable])

    # finally, we output the composition of selected states...

    states_to_print_composition_list = ['s2']

    if facility_type == 'expansion_tube':
        states_to_print_composition_list += ['s7']

    if nozzle_flag:
        states_to_print_composition_list += ['s8']

    states_to_print_composition_list += ['s10e']

    # now get the composition we want to do:

    outputUnits = config_data['outputUnits']

    for state_name in states_to_print_composition_list:

        if states_dict:
            facility_state = states_dict[state_name]
            facilty_state_dictionary_output = facility_state.get_dictionary_output()

            gmodel_type = facility_state.get_gas_state_gmodel_type()

        else:
            facilty_state_dictionary_output = json_output_states_dict[state_name]

            gmodel_type = facilty_state_dictionary_output['gmodel_type']

        if gmodel_type == 'CEAGas':

            if outputUnits == 'moles':
                composition = facilty_state_dictionary_output['composition_moles']
                composition_type_variable = 'X'

            elif outputUnits == 'massf':
                composition = facilty_state_dictionary_output['composition_massf']
                composition_type_variable = 'c'

            for species in composition.keys():

                variable_name = f"{state_name} {composition_type_variable}{species}"

                title_list.append(variable_name)
                value_list.append(composition[species])

    # now we output the results to a file

    title_line = ''

    for variable_name in title_list:
        if variable_name != title_list[-1]:
            title_line += f'{variable_name},'
        else:
            title_line += f'{variable_name}'

    result_line = ''

    for value in value_list:
        if value != value_list[-1]:
            result_line += f'{value},'
        else:
            result_line += f'{value}'

    if output_to_file:

        output_filename = config_data['output_filename'] + '_one_line_output.csv'

        with open(output_filename, 'w') as output_file:
            output_file.write(title_line + '\n')
            output_file.write(result_line)

    if not return_title_and_result_lists:
        return title_line, result_line
    else:
        return title_line, result_line, title_list, value_list

def pitot3_remote_run_creator(config_dict, run_folder, pitot_3_input_file_filename):
    """
    Function to make a folder with a given name and put a PITOT3 input file in it which can be ran later.

    :param config_dict:
    :param run_folder:
    :return:
    """

    import os

    print('-'*60)
    print(f"Making a folder called {run_folder} and putting a PITOT3.yaml input file called {pitot_3_input_file_filename} in it.")
    print('-'*60)

    if not os.path.exists(run_folder):
        os.mkdir(run_folder)
    else:
        print("Just a warning that the specified run folder already exists, so it isn't being created now.")

    if pitot_3_input_file_filename[-5:] != '.yaml':
        pitot_3_input_file_filename += '.yaml'

    pitot3_input_yaml_file_creator(config_dict, run_folder + '/' + pitot_3_input_file_filename)

    return

class Facility(object):
    """
    Class to load in and store test facility information for PITOT3.
    """

    def __init__(self, cfg):
        """

        :param cfg: Python dictionary which contains configuration information which is loaded into the class,
                            this might have come from some kind of config file (which has already been loaded),
                            but could also just be generated in Python.

        """

        # -------------------- initial setup ----------------------------

        # I am going to keep this as a private variable here, as I worry that something may be updated elsewhere in the class
        # and not here, which would lead to some confusion if people then want to pull variables out of here...
        # (I also may just let the variable lapse in future versions...)
        self.__cfg = cfg

        self.facility_name = cfg['facility_name']
        self.facility_type = cfg['facility_type']
        # could just set these to False if they aren't in the file...
        self.secondary_driver = cfg['secondary_driver']
        self.nozzle = cfg['nozzle']

        self.driver_conditions_folder = cfg['driver_conditions_folder']

        # all facilities will have shock tubes, but only some will have secondary drivers, acceleration tubes, and nozzles
        # so we check for this
        # (will add error trapping later on...)

        #-------------------- shock tube ----------------------------
        if 'shock_tube_length' in cfg:
            self.shock_tube_length = float(cfg['shock_tube_length'])
        else:
            self.shock_tube_length = None

        if 'shock_tube_diameter' in cfg:
            self.shock_tube_diameter = float(cfg['shock_tube_diameter'])
        else:
            self.shock_tube_diameter = None
        if 'shock_tube_sensors' in cfg:
            self.shock_tube_sensors  = cfg['shock_tube_sensors']
        else:
            self.shock_tube_sensors = None

        #-------------------- acceleration tube ----------------------------
        if self.facility_type == 'expansion_tube':
            if 'acceleration_tube_length' in cfg:
                self.acceleration_tube_length = float(cfg['acceleration_tube_length'])
            else:
                self.acceleration_tube_length = None

            if 'acceleration_tube_diameter' in cfg:
                self.acceleration_tube_diameter = float(cfg['acceleration_tube_diameter'])
            else:
                self.acceleration_tube_diameter = None
            if 'acceleration_tube_sensors' in cfg:
                self.acceleration_tube_sensors = cfg['acceleration_tube_sensors']
            else:
                self.acceleration_tube_sensors = None

        #-------------------- secondary driver tube ----------------------------
        if self.secondary_driver :
            if 'secondary_driver_tube_length' in cfg:
                self.secondary_driver_tube_length = float(cfg['secondary_driver_tube_length'])
            else:
                self.secondary_driver_tube_length = None

            if 'secondary_driver_tube_diameter' in cfg:
                self.secondary_driver_tube_diameter = float(cfg['secondary_driver_tube_diameter'])
            else:
                self.secondary_driver_tube_diameter = None
            if 'secondary_driver_tube_sensors' in cfg:
                self.secondary_driver_tube_sensors = cfg['secondary_driver_tube_sensors']
            else:
                self.secondary_driver_tube_sensors = None

        #---------------------------- nozzle ----------------------------
        if self.nozzle:
            self.nozzle_geometric_area_ratio = cfg['nozzle_geometric_area_ratio']

        #---------------------------- sensor locations ----------------------------
        if 'sensor_locations' in cfg:
            self.sensor_locations = cfg['sensor_locations']
        else:
            self.sensor_locations = None

        return

    def get_facility_name(self):
        """
        Return the facility name
        """

        return self.facility_name

    def get_facility_type(self):
        """
        Return the facility type
        """

        return self.facility_type

    def get_secondary_driver_flag(self):
        """
        Return the secondary driver flag
        :return:
        """

        return self.secondary_driver

    def get_nozzle_flag(self):
        """
        Return the nozzle flag
        :return:
        """

        return self.nozzle

    def get_secondary_driver_length_and_diameter(self):
        """
        Return the secondary driver length and diameter
        :return:
        """
        return self.secondary_driver_tube_length, self.secondary_driver_tube_diameter

    def get_shock_tube_length_and_diameter(self):
        """
        Return the shock tube length and diameter
        :return:
        """
        return self.shock_tube_length, self.shock_tube_diameter

    def get_acceleration_tube_length_and_diameter(self):
        """
        Return the acceleration tube length and diameter
        :return:
        """

        # TODO: probably needs some error trapping or something like that for when it isn't an expansion tube...
        return self.acceleration_tube_length, self.acceleration_tube_diameter

    def get_sensor_locations(self):
        """
        Return the sensor locations dictionary, if it exists
        """

        if self.sensor_locations:
            return self.sensor_locations
        else:
            print("There are no sensor locations loaded for this facility. None will be returned")
            return None

    def get_nozzle_geometric_area_ratio(self):
        """
        Return the nozzle geometric area ratio
        :return:
        """

        if self.nozzle and self.nozzle_geometric_area_ratio:
            return self.nozzle_geometric_area_ratio
        else:
            print("Either this configuration does not have a nozzle or no nozzle geometric area ratio has been provided. None will be returned.")
            return None

    def get_driver_conditions_folder(self):
        """
        Return driver conditions folder
        """

        return self.driver_conditions_folder

class Driver(object):
    """
    Class to store and calculate facility driver information.
    """

    def __init__(self, cfg, p_0 = None, T_0 = None, preset_gas_models_folder = None, outputUnits = 'massf', species_MW_dict = None, D_shock_tube = None):
        """

        :param cfg: Python dictionary which contains configuration information which is loaded into the class,
                            this might have come from some kind of config file (which has already been loaded),
                            but could also just be generated in Python.
        :param p_0: standard atmospheric pressure, in Pa, used to make a reference state for the driver.
        :param T_0: standard atmospheric temperature, in K, used to make a reference state for the driver.

        """

        # -------------------- initial setup ----------------------------

        # I am going to keep this as a private variable here, as I worry that something may be updated elsewhere in the class
        # and not here, which would lead to some confusion if people then want to pull variables out of here...
        # (I also may just let the variable lapse in future versions...)
        self.__cfg = cfg

        self.driver_condition_name = cfg['driver_condition_name']
        self.driver_condition_type = cfg['driver_condition_type']

        # -------------------- gas composition ----------------------------

        self.driver_gas_model = cfg['driver_gas_model']

        # the user specified CEA model for now...
        if self.driver_gas_model == 'CEAGas':

            self.driver_fill_composition = cfg['driver_fill_composition']
            self.driver_speciesList = cfg['driver_speciesList']
            self.driver_inputUnits = cfg['driver_inputUnits']
            self.driver_withIons = cfg['driver_withIons']

            # now build a gas model file and attach it to this object...

            eilmer4_CEAGas_input_file_creator('PITOT3_cea_driver_condition', 'driver_gas', self.driver_speciesList, self.driver_fill_composition,
                                              self.driver_inputUnits, self.driver_withIons)

            driver_gmodel_location = 'PITOT3_cea_driver_condition.lua'

        elif self.driver_gas_model == 'CEAGas-preset':

            self.driver_gas_model = 'CEAGas'
            self.driver_fill_gas_name = cfg['driver_fill_gas_name']

            driver_gmodel_location = '{0}/cea-{1}-gas-model.lua'.format(preset_gas_models_folder, self.driver_fill_gas_name)

        elif self.driver_gas_model == 'ideal-preset':

            self.driver_gas_model = 'ideal'
            self.driver_fill_gas_name = cfg['driver_fill_gas_name']

            driver_gmodel_location = '{0}/ideal-{1}-gas-model.lua'.format(preset_gas_models_folder, self.driver_fill_gas_name)

        elif self.driver_gas_model == 'thermally-perfect-preset':
            self.driver_gas_model = 'thermally-perfect'
            self.driver_fill_gas_name = cfg['driver_fill_gas_name']

            # we need the composition for the thermally perfect gas as the gas model file does not store it (just the species)
            self.driver_fill_composition = cfg['driver_fill_composition']
            self.driver_inputUnits = cfg['driver_inputUnits']

            driver_gmodel_location = '{0}/thermally-perfect-{1}-gas-model.lua'.format(preset_gas_models_folder,
                                                                                      self.driver_fill_gas_name)

            # the gas state input requires mass fractions, so if we inputted them great, if not, we need to calculate them...
            if self.driver_inputUnits == 'massf':
                self.driver_fill_composition_massf = self.driver_fill_composition
            else:
                # we need the molecular masses to get from mole fractions to mass fractions
                # we will need to make a gas model to do that...

                gmodel = GasModel(os.path.expandvars(driver_gmodel_location))

                self.driver_speciesList = gmodel.species_names

                molecular_mass_list = gmodel.mol_masses

                molecular_mass_dict = {}

                for species, molecular_mass in zip (self.driver_speciesList, molecular_mass_list):
                    molecular_mass_dict[species] = molecular_mass

                # now we need the total molecular mass for our calculation
                total_molecular_mass = 0.0

                for species in self.driver_speciesList:
                    total_molecular_mass += molecular_mass_dict[species]*self.driver_fill_composition[species]

                # which we can now use to turn out mole fractions into mass fractions:

                self.driver_fill_composition_massf = {}

                for species in self.driver_speciesList:
                    self.driver_fill_composition_massf[species] = self.driver_fill_composition[species]*(molecular_mass_dict[species]/ total_molecular_mass)

                # TO DO: add making sure that the mass fractions add up to 1...

        elif self.driver_gas_model == 'custom' and cfg['driver_fill_gas_filename']:

            self.driver_fill_gas_filename = cfg['driver_fill_gas_filename']

            driver_gmodel_location = self.driver_fill_gas_filename

        self.gmodel = GasModel(os.path.expandvars(driver_gmodel_location))

        # now what we have to do now depends on the type of driver condition
        # we just have empirical for now, which is easy, as they just specify p4 and T4...

        # when do i want to make this stuff facility gas objects?

        if self.driver_condition_type == 'empirical':
            # we just need to take p4 and T4 out of the file and then make a gas object...
            self.p4 = float(cfg['p4'])
            self.T4 = float(cfg['T4'])

        elif self.driver_condition_type == 'cold-driven-ideal':
            # this is similar to empirical above, but in this case the driver
            # is assumed to be at room temperature so PITOT's ambient temp can be used
            # and it doesn't need to be specified...
            self.p4 = float(cfg['p4'])
            self.T4 = T_0

        elif self.driver_condition_type in ['isentropic-compression-p4', 'isentropic-compression-compression-ratio']:
            # do isentropic compression from the fill state to p4
            self.driver_p = float(cfg['driver_p'])
            if 'driver_T' in cfg:
                self.driver_T = float(cfg['driver_T'])
            else:
                print(f"Setting driver temperature to the default value of {T:.2f} K")
                self.driver_T = T_0

            # we call this state 4i
            state4i = GasState(self.gmodel)
            state4i.p = self.driver_p
            state4i.T = self.driver_T

            if self.driver_gas_model == 'thermally-perfect':
                state4i.massf = self.driver_fill_composition_massf

            state4i.update_thermo_from_pT()
            state4i.update_sound_speed()

            gamma = state4i.gamma

            # we assume that the state 4 driver is stationary
            v4i = 0.0

            # make a reference gas state here...
            reference_gas_state = GasState(self.gmodel)
            reference_gas_state.p = p_0
            reference_gas_state.T = T_0

            if self.driver_gas_model == 'thermally-perfect':
                reference_gas_state.massf = self.driver_fill_composition_massf

            reference_gas_state.update_thermo_from_pT()
            reference_gas_state.update_sound_speed()

            # now make our facility driver object...
            self.state4i = Facility_State('s4i', state4i, v4i, reference_gas_state=reference_gas_state,
                                          species_MW_dict=species_MW_dict, outputUnits=outputUnits)

            if self.driver_condition_type == 'isentropic-compression-p4':
                self.p4 = float(cfg['p4'])

                print (f"Performing isentropic compression from the driver fill condition to {self.p4/1.0e6:.2f} MPa.")

                self.T4 = state4i.T * (self.p4 / state4i.p) ** (1.0 - (1.0 / gamma))  # K

            elif self.driver_condition_type == 'isentropic-compression-compression-ratio':

                self.compression_ratio = cfg['compression_ratio']

                print (f"Performing isentropic compression from driver fill condition over compression ratio of {self.compression_ratio}.")
                pressure_ratio =  self.compression_ratio**gamma #pressure ratio is compression ratio to the power of gamma

                self.p4 = state4i.p*pressure_ratio

            self.T4 = state4i.T * (self.p4 / state4i.p) ** (1.0 - (1.0 / gamma))  # K

        state4 = GasState(self.gmodel)
        state4.p = self.p4
        state4.T = self.T4

        if self.driver_gas_model == 'thermally-perfect':
            state4.massf = self.driver_fill_composition_massf

        state4.update_thermo_from_pT()
        state4.update_sound_speed()

        # we assume that the state 4 driver is stationary
        v4 = 0.0

        # I have had to add some of the room temperature only gas state stuff to here for detonation driver
        # configurations, as I was having some issues setting the reference gas state with H2/O2/He gas compositions:

        if self.driver_gas_model == 'CEAGas' and hasattr(self, 'driver_fill_gas_name') and self.driver_fill_gas_name:
            fill_gmodel_location = '{0}/cea-{1}-gas-model.lua'.format(preset_gas_models_folder, self.driver_fill_gas_name)
            # we create this link and then use if it exists...
            fill_room_temperature_only_gmodel_location = '{0}/cea-{1}-room-temperature-only-gas-model.lua'.format(preset_gas_models_folder, self.driver_fill_gas_name)
        elif self.driver_gas_model == 'custom' and self.driver_fill_gas_filename:
            fill_gmodel_location = self.driver_fill_gas_filename
            # we make it the gas name with -room-temperature-only-gas-model added... (in the same folder)
            fill_room_temperature_only_gmodel_location = fill_gmodel_location.split('.')[0] + '-room-temperature-only-gas-model.' + fill_gmodel_location.split('.')[-1]
        else:
            fill_room_temperature_only_gmodel_location = None

        if fill_room_temperature_only_gmodel_location and os.path.exists(os.path.expandvars(fill_room_temperature_only_gmodel_location)):

            # if there is a room temperature only object we "trick" the fill state gas model to use that when it sets the gas
            # state and then replace the correct gas state after. This seemed to be the best / cleanest way to do it.

            room_temperature_only_gmodel = GasModel(os.path.expandvars(fill_room_temperature_only_gmodel_location))

            reference_gas_state = GasState(self.gmodel)

            reference_gas_state.gmodel = room_temperature_only_gmodel

            reference_gas_state.p = p_0
            reference_gas_state.T = T_0

            reference_gas_state.update_thermo_from_pT()
            reference_gas_state.update_sound_speed()

            reference_gas_state.gmodel = self.gmodel

        else:

            # we don't do anything abnormal...

            room_temperature_only_gmodel = None

            # make a reference gas state here...
            reference_gas_state = GasState(self.gmodel)
            reference_gas_state.p = p_0
            reference_gas_state.T = T_0

            if self.driver_gas_model == 'thermally-perfect':
                reference_gas_state.massf = self.driver_fill_composition_massf

            reference_gas_state.update_thermo_from_pT()
            reference_gas_state.update_sound_speed()

        # now make our facility driver object...
        self.state4 = Facility_State('s4', state4, v4,
                                     reference_gas_state=reference_gas_state, room_temperature_only_gmodel = room_temperature_only_gmodel,
                                     outputUnits=outputUnits, species_MW_dict=species_MW_dict)

        #-------------------- M throat ---------------------------------------------------------

        if 'M_throat' in cfg:

            self.M_throat = cfg['M_throat']

            print(f'M_throat = {self.M_throat}')

        elif 'M_throat' not in cfg and 'D_throat' in cfg and D_shock_tube:
            # we can calculate M_throat ourselves using the throat diameter and the shock tube diameter

            from gdtk.ideal_gas_flow import A_Astar
            from scipy.optimize import newton

            print("Calculating M_throat from the provided D_throat and D_shock_tube values:")

            D_throat = cfg['D_throat']
            A_throat = (math.pi/4.0)*D_throat**2.0

            A_shock_tube = (math.pi/4.0)*D_shock_tube**2.0

            print(f'D_throat = {D_throat} m ({D_throat*1000.0} mm)')
            print(f'D_shock_tube = {D_shock_tube} m ({D_shock_tube*1000.0} mm)')

            print(f'A_throat = {A_throat:.2e} m**2 ({A_throat*1.0e6:.2f} mm**2)')
            print(f'A_shock_tube = {A_shock_tube:.2e} m**2 ({A_shock_tube*1.0e6:.2f} mm**2)')

            A_shock_tube_over_A_throat = A_shock_tube/A_throat

            print(f"A_shock_tube / A_throat = {A_shock_tube_over_A_throat:.2f}")

            M_throat_eqn = lambda M_throat : A_shock_tube_over_A_throat - A_Astar(M_throat, g=self.state4.get_gas_state().gamma)

            self.M_throat = newton(M_throat_eqn, 1.1)

            print(f"Calculated M_throat = {self.M_throat:.2f}")

        if self.M_throat > 0.0:
            self.driver_exit_state_name = '3s'
        else:
            self.driver_exit_state_name = '4'

        #--------------------------------------------------------------------------------------

        if self.driver_exit_state_name == '3s':
            driver_gas_flow_object = GasFlow(self.gmodel)

            state3s = GasState(self.gmodel)

            v3s = driver_gas_flow_object.expand_from_stagnation(self.state4.get_gas_state(),
                                                                1.0 / p0_p(self.M_throat, self.state4.get_gas_state().gamma),
                                                                state3s)

            # expand to Mach is having issues, so we will just try expand_from_stagnation which we used in the original version of PITOT...
            # if self.M_throat > 1.0:
            #     # we do it in two steps here, once to Mach 1 and then once beyond it....
            #
            #     state3s_intermediate = GasState(self.gmodel)
            #     v3s_intermediate = driver_gas_flow_object.expand_to_mach(self.state4.get_gas_state(), 1.0, state3s_intermediate)
            #     print(state3s_intermediate)
            #     print(v3s_intermediate)
            #     print(self.M_throat)
            #     v3s = driver_gas_flow_object.expand_to_mach(state3s_intermediate, self.M_throat, state3s)
            #     print(state3s)
            #     print(v3s)
            #     print (v3s/state3s.a)
            # else:
            #     v3s = driver_gas_flow_object.expand_to_mach(self.state4.get_gas_state(), self.M_throat, state3s)

            # now make our facility driver object...
            # user the same reference gas state as above...
            self.state3s = Facility_State('s3s', state3s, v3s,
                                          reference_gas_state=reference_gas_state,
                                          outputUnits=outputUnits, species_MW_dict=species_MW_dict, )

        return

    def get_driver_condition_name(self):
        """
        Return the facility name
        """

        return self.driver_condition_name

    def get_driver_condition_type(self):
        """
        Return the driver condition type
        """

        return self.driver_condition_type

    def get_exit_state_name(self):
        """
        Return the driver exit state string.
        :return:
        """

        return self.driver_exit_state_name

    def get_M_throat(self):
        """
        Return the M_throat value
        :return:
        """

        return self.M_throat

    def get_exit_state(self):
        """
        Return the exit state
        """

        # if we have a state 3s, that is it, otherwise, it is state4
        if hasattr(self, 'state3s'):
            return self.state3s
        else:
            return self.state4

    def get_driver_rupture_state(self):
        """
        Return the driver rupture condition (state 4, by convention).
        This returns the gas state, the velocity, and the Mach number.
        """

        return self.state4

    def get_driver_throat_state(self):
        """
        Return the driver throat condition if it exists
        """

        if hasattr(self, 'state3s'):
            return self.state3s
        else:
            print("This driver condition does not have a state 3s")

    def get_facility_states(self):
        """
        This function returns a list of the unique facility states (i.e. it shouldn't return anything
        which would be returned using the function with the same name on another object) which this object contains,
        in an order which is relevant for the PITOT3 text output.

        This was mainly made for the output, but may have other uses.

        """

        facility_states = []

        if hasattr(self, 'state4i'):
            facility_states.append(self.state4i)

        facility_states.append(self.state4)

        if hasattr(self, 'state3s'):
            facility_states.append(self.state3s)

        return facility_states

class Diaphragm(object):
    """
    Class for diaphragm objects in PITOT3. Mainly for when we are considering non-ideal diaphragm models,
    but will keep it here even when working with ideal diaphragms so that we have it around for those conditions later on.
    """

    def __init__(self, cfg):
        """

        :param cfg: Python dictionary which contains configuration information which is loaded into the class,
                            this might have come from some kind of config file (which has already been loaded),
                            but could also just be generated in Python.

        """

        # -------------------- initial setup ----------------------------

        # I am going to keep this as a private variable here, as I worry that something may be updated elsewhere in the class
        # and not here, which would lead to some confusion if people then want to pull variables out of here...
        # (I also may just let the variable lapse in future versions...)
        self.__cfg = cfg

        self.diaphragm_name = cfg['diaphragm_name']
        self.diaphragm_type = cfg['diaphragm_type']

        # may need to decide at some point if we just need labels here or real objects...
        self.diaphragm_entrance_state_name = cfg['diaphragm_entrance_state_name']
        self.diaphragm_entrance_state = cfg['diaphragm_entrance_state']

        if self.diaphragm_type in ['reflected_shock', 'velocity_loss_factor_and_reflected_shock']:
            # grab our input Mach number
            if 'Mr_input' not in cfg:
                print("Setting the reflected shock at the diaphragm to 'maximum' as the user has not specified a value.")
                self.Mr_input = 'maximum'
            else:
                self.Mr_input = cfg['Mr_input']

        if self.diaphragm_type in ['velocity_loss_factor', 'velocity_loss_factor_and_reflected_shock']:
            self.velocity_loss_factor = cfg['velocity_loss_factor']

        if self.diaphragm_type == 'ideal':
            # ideal just means that the diaphragm does nothing...
            # but more complex models, we would need to do calculations in here...
            self.diaphragm_exit_state_name = self.diaphragm_entrance_state_name
            self.diaphragm_exit_state = self.diaphragm_entrance_state

        elif self.diaphragm_type == 'reflected_shock':

            self.perform_reflected_shock_at_diaphragm(entrance_state=self.diaphragm_entrance_state)

            self.diaphragm_exit_state_name = self.diaphragm_entrance_state_name + 'r'
            self.diaphragm_exit_state = self.reflected_shocked_state

        elif self.diaphragm_type == 'velocity_loss_factor':

            self.perform_velocity_loss_factor_at_diaphragm()

            self.diaphragm_exit_state_name = self.diaphragm_entrance_state_name + 'l'
            self.diaphragm_exit_state = self.velocity_loss_factor_state

        elif self.diaphragm_type == 'velocity_loss_factor_and_reflected_shock':

            self.perform_velocity_loss_factor_at_diaphragm()

            self.perform_reflected_shock_at_diaphragm(entrance_state=self.velocity_loss_factor_state)

            self.diaphragm_exit_state_name = self.diaphragm_entrance_state_name + 'r'
            self.diaphragm_exit_state = self.reflected_shocked_state

        return

    def perform_reflected_shock_at_diaphragm(self, entrance_state):
        """
        Function to perform a reflected shock at the diaphragm.
        :return:
        """

        print('-' * 60)
        print(f"Doing reflected shock at the {self.diaphragm_name} which the user has asked for.")

        # make a gas state to work with...

        entrance_gas_gas_state = entrance_state.get_gas_state()
        entrance_gas_v = entrance_state.get_v()

        entrance_gas_gmodel = entrance_gas_gas_state.gmodel
        entrance_gas_gas_flow = GasFlow(entrance_gas_gmodel)

        reflected_shocked_gas_gas_state = GasState(entrance_gas_gmodel)

        # if the value is maximum, we just have to do a full reflected shock...
        if self.Mr_input == 'maximum':
            print("Performing the maximum strength reflected shock.")
            self.vr = entrance_gas_gas_flow.reflected_shock(entrance_gas_gas_state, entrance_gas_v,
                                                            reflected_shocked_gas_gas_state)
            self.Mr = (
                                  entrance_gas_v + self.vr) / entrance_gas_gas_state.a  # normally this would be V2 - Vr, but it's plus here as Vr has been left positive

            reflected_shocked_gas_v = 0.0
        else:
            print(f"Performing a reflected shock with a user selected Mach number of {self.Mr_input}")
            # we have a numerical value, so things are a bit harder...
            self.Mr = self.Mr_input
            self.vr = self.Mr * entrance_gas_gas_state.a - entrance_gas_v

            v2r, v2rg = entrance_gas_gas_flow.normal_shock(entrance_gas_gas_state, self.vr + entrance_gas_v,
                                                         reflected_shocked_gas_gas_state)

            reflected_shocked_gas_v = entrance_gas_v - v2rg

        if entrance_state.reference_gas_state:
            reference_gas_state = entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None

        self.reflected_shocked_state = Facility_State(entrance_state.get_state_name() + 'r',
                                                      reflected_shocked_gas_gas_state, reflected_shocked_gas_v,
                                                      reference_gas_state=reference_gas_state,
                                                      outputUnits=entrance_state.outputUnits,
                                                      species_MW_dict=entrance_state.species_MW_dict)

        print(f'Mr = {self.Mr:.2f}, vr = {self.vr:.2f}')
        print(self.reflected_shocked_state)

        if self.reflected_shocked_state.get_gas_state_gmodel_type() == 'CEAGas':
            print(self.reflected_shocked_state.get_reduced_composition_two_line_output_string())

        return

    def perform_velocity_loss_factor_at_diaphragm(self):

        print('-' * 60)
        print(f"Multiplying the velocity at the {self.diaphragm_name} by a velocity_loss_factor of {self.velocity_loss_factor}")

        entrance_gas_gas_state = self.diaphragm_entrance_state.get_gas_state()
        entrance_gas_v = self.diaphragm_entrance_state.get_v()

        if self.diaphragm_entrance_state.reference_gas_state:
            reference_gas_state = self.diaphragm_entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None

        self.velocity_loss_factor_state = Facility_State(self.diaphragm_entrance_state_name + 'l',
                                                         entrance_gas_gas_state,
                                                         entrance_gas_v * self.velocity_loss_factor,
                                                         reference_gas_state=reference_gas_state,
                                                         outputUnits=self.diaphragm_entrance_state.outputUnits,
                                                         species_MW_dict=self.diaphragm_entrance_state.species_MW_dict)
        print(self.velocity_loss_factor_state)

        return



    def get_diaphragm_name(self):
        """
        Return the diaphragm name.
        """

        return self.diaphragm_name

    def get_diaphragm_type(self):
        """
        Return the diaphragm type.
        """

        return self.diaphragm_type

    def get_entrance_state_name(self):
        """
        Return the name of the diaphragm entrance state.
        """

        return self.diaphragm_entrance_state_name

    def get_entrance_state(self):
        """
        Return the diaphragm entrance state.
        """

        return self.diaphragm_entrance_state

    def get_exit_state_name(self):
        """
        Return the name of the diaphragm exit state.
        """

        return self.diaphragm_exit_state_name

    def get_exit_state(self):
        """
        Return the diaphragm exit state.
        """

        return self.diaphragm_exit_state

    def get_vr_Mr(self):
        """
        Return the vr and Mr values if they have been calculated
        :return:
        """

        if self.diaphragm_type in ['reflected_shock', 'velocity_loss_factor_and_reflected_shock']:
            return self.vr, self.Mr
        else:
            return None
    def get_facility_states(self):
        """
        This function returns a list of the unique facility states (i.e. it shouldn't return anything
        which would be returned using the function with the same name on another object) which this object contains,
        in an order which is relevant for the PITOT3 text output.

        This was mainly made for the output, but may have other uses.

        """

        if self.diaphragm_type == 'ideal':
            # nothing unique to return in this case...
            return []
        elif self.diaphragm_type == 'reflected_shock':
            return [self.reflected_shocked_state]
        elif self.diaphragm_type == 'velocity_loss_factor':
            return [self.velocity_loss_factor_state]
        elif self.diaphragm_type == 'velocity_loss_factor_and_reflected_shock':
            return [self.velocity_loss_factor_state, self.reflected_shocked_state]

class Facility_State(object):
    """
    This is the class for storing the 'facility states' in PITOT3 i.e. the numbers on an x-t diaphragm of an impulse facility.
    I have avoided calling it 'Facility_Gas_State' as it is not actually a gas state, as it stores velocity information
    (and any other important information which we decide we need) for particular states as they move around the code.

    In the old version of PITOT this was in the information stored in the states, V, and M dictionaries, which here I have
    decided to build into one simple object, as we generally always want them together, so I created this very simple wrapper for them.

    As it has the velocity, I have also added the ability for it to calculate the pitot and total condition!

    """

    def __init__(self, state_name, gas_state, v, reference_gas_state = None, room_temperature_only_gmodel = None,
                 outputUnits = 'massf', species_MW_dict = None, make_gmodel_without_ions = None):
        """

        :param state_name: just a string labelled the state, probably just '4', '3s', '1' etc.
        :param gas_state: eilmer gas object for the state
        :param v: velocity related to the object. Mach number will be found from this and the gas state too.
        :param reference_gas_state: a reference gas state which we can use to get the stagnation enthalpy of the object...

        """

        self.state_name = state_name
        self.gas_state = gas_state
        self.v = v
        self.M = self.v/self.gas_state.a

        self.reference_gas_state = reference_gas_state

        self.gas_state_gmodel = self.gas_state.gmodel

        # we save this here for when we reload the objects and may not have the gasmodel itself handy...
        self.gas_state_gmodel_typestr = self.gas_state_gmodel.type_str

        # we make a no-ions version of the gmodel if needed.
        # (we do this for any CEA gas model with ions
        # as they can fail at low temperatures, and sometimes
        # we want to be able to keep going then...)

        if self.get_gas_state_gmodel_type() == 'CEAGas':
            # this assumes the gas state has been created properly at some point...
            # (I added the ability to have a flag here as there was a bug where if the gas state was
            # made with the ions turned off... then the next gas state wouldn't make this object...)
            if 'e-' in self.gas_state.ceaSavedData['massf'] or make_gmodel_without_ions:
                gmodel_without_ions_filename = eilmer4_CEAGas_gmodel_without_ions_creator(self.gas_state_gmodel.file_name)

                self.gas_state_gmodel_without_ions = GasModel(gmodel_without_ions_filename)

            else:
                self.gas_state_gmodel_without_ions = None
        else:
            self.gas_state_gmodel_without_ions = None

        self.room_temperature_only_gmodel = room_temperature_only_gmodel

        self.outputUnits = outputUnits

        if self.get_gas_state_gmodel_type() == 'CEAGas' and self.outputUnits == 'moles' and not species_MW_dict:
            raise Exception("pitot3_classes.Facility_State: Cannot select OutputUnits of 'moles' without providing a species molecular weight dictionary.")

        # in the end, doing it this way did not work when I did calculations without ions,
        # so I'll just store the whole species MW dict with every facility state... it isn't a big deal.
        #if species_MW_dict and self.get_gas_state_gmodel_type() == 'CEAGas':
        #    # make a species MW dict with just the species in the gas...
        #    self.species_MW_dict = {}
        #
        #    for species in self.gas_state.ceaSavedData['massf'].keys():
        #        self.species_MW_dict[species] = species_MW_dict[species]

        if self.get_gas_state_gmodel_type() == 'CEAGas':
            self.species_MW_dict = species_MW_dict
        else:
            self.species_MW_dict = None

        return

    def __str__(self):
        text = "Facility_State(state {0}: p = {1:.2f} Pa, T = {2:.2f} K, gam = {3:.2f}, R = {4:.2f} J/kg K, V = {5:.2f} m/s, M = {6:.2f})".format(self.state_name,
                                                                                                               self.gas_state.p,
                                                                                                               self.gas_state.T,
                                                                                                               self.gas_state.gamma,
                                                                                                               self.gas_state.R,
                                                                                                               self.v, self.M)

        return text

    def get_state_name(self):
        """
        Return just the state name.
        """

        return self.state_name

    def get_gas_state(self):
        """
        Return just the gas state.
        """

        return self.gas_state

    def get_gas_state_gmodel(self):
        """
        Return just the gas model of the gas state.
        """

        return self.gas_state_gmodel

    def get_gas_state_gmodel_without_ions(self):
        """
        Return the gas model without ions if we have it...
        """

        return self.gas_state_gmodel_without_ions

    def get_room_temperature_only_gmodel(self):
        """
        Return the gas model without ions if we have it...
        """

        return self.room_temperature_only_gmodel

    def get_gas_state_gmodel_type(self):
        """
        Return just the type of the gas model's gas state
        """

        return self.gas_state_gmodel_typestr

    def get_v(self):
        """
        Return just the velocity.
        """

        return self.v

    def get_M(self):
        """
        Return just the Mach number.
        """

        return self.M

    def get_gas_state_and_v(self):
        """
        Return the gas state and the velocity.
        """

        return self.gas_state, self.v

    def calculate_pitot_condition(self):
        """
        Function to calculate the pitot condition for the facility state, which can then be returned using the related function.

        Remember that this state doesn't need the velocity object as it has been brought to rest.
        :return:
        """

        # start by making a gas flow object which we can use...
        gas_flow_object = GasFlow(self.gas_state.gmodel)

        # TODO: add stuff to check that this actually worked.

        try:
            pitot_state = GasState(self.gas_state.gmodel)
            gas_flow_object.pitot_condition(self.gas_state, self.v, pitot_state)
            self.pitot_state = pitot_state
        except Exception as err:
            print("Error: {0}".format(err))
            print(
                "pitot3_classes.Facility_State.calculate_pitot_condition() Failed to calculate pitot condition for state {0}".format(
                    self.state_name))
            print("self.pitot_state will be set to None.")
            self.pitot_state = None

        return

    def get_pitot_condition(self):
        """
        Returns the pitot condition. Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'pitot_state') and self.pitot_state:
            return self.pitot_state
        elif hasattr(self, 'pitot_state') and not self.pitot_state:
            # means that it failed to calculate the Pitot state
            print("Cannot return the pitot condition as PITOT3 failed to calculate it in the past.")
            print("Will return None.")
            return None
        else:
            self.calculate_pitot_condition()
            if not self.pitot_state:
                print("Cannot return the pitot condition as PITOT3 failed to calculate it.")
                print("Will return None.")
                return None
            else:
                return self.pitot_state

    def get_pitot_pressure(self):
        """
        Returns the pitot pressure (the pressure for the pitot condition). Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'pitot_state') and self.pitot_state:
            return self.pitot_state.p
        elif hasattr(self, 'pitot_state') and not self.pitot_state:
            # means that it failed to calculate the Pitot state
            print("Cannot return the pitot pressure as PITOT3 failed to calculate the pitot condition in the past.")
            print("Will return None.")
            return None
        else:
            self.calculate_pitot_condition()
            if not self.pitot_state:
                print("Cannot return the pitot pressure as PITOT3 failed to calculate it.")
                print("Will return None.")
                return None
            else:
                return self.pitot_state.p

    def calculate_total_condition(self):
        """
        Function to calculate the total condition for the facility state, which can then be returned using the related function.

        Remember that this state doesn't need the velocity object as it has been brought to rest.
        :return:
        """

        # start by making a gas flow object which we can use...
        gas_flow_object = GasFlow(self.gas_state.gmodel)

        # TODO: add stuff to check that this actually worked.

        try:
            total_state = GasState(self.gas_state.gmodel)
            gas_flow_object.total_condition(self.gas_state, self.v, total_state)
            self.total_state = total_state
        except Exception as err:
            print("Error: {0}".format(err))
            print("pitot3_classes.Facility_State.calculate_total_condition() Failed to calculate total condition for state {0}".format(self.state_name))
            print("self.total_state will be set to None.")
            self.total_state = None

        return

    def get_total_condition(self):
        """
        Returns the total condition. Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'total_state') and self.total_state:
            return self.total_state
        elif hasattr(self, 'total_state') and not self.total_state:
            # means that it failed to calculate the total state
            print("Cannot return the total condition as PITOT3 failed to calculate the total condition in the past.")
            print("Will return None.")
            return None
        else:
            self.calculate_total_condition()
            if not self.total_state:
                print("Cannot return the total condition as PITOT3 failed to calculate it.")
                print("Will return None.")
                return None
            else:
                return self.total_state

    def get_total_pressure(self):
        """
        Returns the total pressure (the pressure of the total condition). Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'total_state') and self.total_state:
            return self.total_state.p
        elif hasattr(self, 'total_state') and not self.total_state:
            # means that it failed to calculate the total state
            print("Cannot return the total pressure as PITOT3 failed to calculate the total condition in the past.")
            print("Will return None.")
            return None
        else:
            self.calculate_total_condition()
            if not self.total_state:
                print("Cannot return the total pressure as PITOT3 failed to calculate the total condition.")
                print("Will return None.")
                return None
            else:
                return self.total_state.p

    def get_total_temperature(self):
        """
        Returns the total temperature (the temperature of the total condition). Will calculate it if it doesn't exist yet.
        :return:
        """

        if hasattr(self, 'total_state') and self.total_state:
            return self.total_state.T
        elif hasattr(self, 'total_state') and not self.total_state:
            # means that it failed to calculate the total state
            print("Cannot return the total temperature as PITOT3 failed to calculate the total condition in the past.")
            print("Will return None.")
            return None
        else:
            self.calculate_total_condition()
            if not self.total_state:
                print("Cannot return the total temperature as PITOT3 failed to calculate the total condition.")
                print("Will return None.")
                return None
            else:
                return self.total_state.T

    def set_reference_gas_state(self, reference_gas_state):
        """
        Function to either add or update the reference gas state on this object.

        """

        self.reference_gas_state = reference_gas_state

        return

    def get_reference_gas_state(self):
        """
        Returns the reference gas state.

        """

        if self.reference_gas_state:
            return self.reference_gas_state
        else:
            print("This FacilityState does not have a reference gas state. None will be returned.")
            return None

        return

    def calculate_total_enthalpy(self):
        """
        If a reference gas state is available, this function will calculate the total enthalpy.

        :return:
        """

        if not self.reference_gas_state:
            print("This FacilityState does not have a reference gas state which is required to calculate total enthalpy.")

        else:
            # we need get the total state to get the total enthalpy
            # if we use this function it will return it if we already have it, or calculate it if we don't
            total_state = self.get_total_condition()

            if total_state:

                reference_gas_state = self.get_reference_gas_state()

                # total enthalpy is the sensible enthalpy of the total condition - the sensible enthalpy of the reference state
                total_enthalpy = total_state.enthalpy - reference_gas_state.enthalpy

                self.total_enthalpy = total_enthalpy
            else:
                print("Total state could not be calculated, so self.total_enthalpy will be set to None.")
                self.total_enthalpy = None

        return

    def get_total_enthalpy(self):
        """
        Function to return the total enthalpy if we have it or can get it.
        """

        # TO DO: we could also just add variables like total_enthalpy to the class when we make it?
        if hasattr(self, 'total_enthalpy'):
            return self.total_enthalpy
        elif self.reference_gas_state:
            # just calculate the enthalpy and then return it...
            self.calculate_total_enthalpy()
            if not self.total_enthalpy:
                print("Cannot return the total enthalpy as PITOT3 failed to calculate it.")
                print("Will return None.")
                return None
            else:
                return self.total_enthalpy
        else:
            print("This FacilityState does not have a reference gas state which is required to calculate enthalpy.")
            print("None will be returned.")
            return None

    def calculate_sensible_enthalpy(self):
        """
        If a reference gas state is available, this function will calculate the sensible enthalpy.

        :return:
        """

        if not self.reference_gas_state:
            print("This FacilityState does not have a reference gas state which is required to calculate sensible enthalpy.")

        else:
            reference_gas_state = self.get_reference_gas_state()

            # then sensible enthalpy is just the sensible enthalpy of the gas state minus the reference state
            gas_state = self.get_gas_state()
            sensible_enthalpy = gas_state.enthalpy - reference_gas_state.enthalpy

            self.sensible_enthalpy = sensible_enthalpy

        return

    def get_sensible_enthalpy(self):
        """
        Function to return the sensible enthalpy if we have it or we can get it.
        """

        # TO DO: we could also just add variables like total_enthalpy to the class when we make it?
        if hasattr(self, 'sensible_enthalpy'):
            return self.sensible_enthalpy
        elif self.reference_gas_state:
            # just calculate the enthalpy and then return it...
            self.calculate_sensible_enthalpy()
            return self.sensible_enthalpy
        else:
            print("This FacilityState does not have a reference gas state which is required to calculate enthalpy.")
            print("None will be returned.")
            return None

    def calculate_flight_equivalent_velocity(self):
        """
        Function to calculate the flight equivalent velocity
        :return:
        """

        if not self.reference_gas_state:
            print("This FacilityState does not have a reference gas state which is required to calculate flight equivalent velocity.")

        self.calculate_total_enthalpy()

        if self.total_enthalpy:
            try:
                self.flight_equivalent_velocity = math.sqrt(2.0 * self.total_enthalpy)
            except Exception as e:
                print(e)
                print("Failed to calculate flight equivalent velocity.")
                print("self.flight_equivalent_velocity will be set to None.")
                self.flight_equivalent_velocity = None
        else:
            print("Total enthalpy could not be calculated, so the flight equivalent velocity could not be calculated either.")
            print("self.flight_equivalent_velocity will be set to None.")
            self.flight_equivalent_velocity = None

        return

    def get_flight_equivalent_velocity(self):
        """
        Function to return the flight equivalent velocity if we have it or can get it.
        """

        if hasattr(self, 'flight_equivalent_velocity'):
            return self.flight_equivalent_velocity
        elif self.reference_gas_state:
            # just calculate the flight equivalent velocity and return it.
            self.calculate_flight_equivalent_velocity()
            if not self.flight_equivalent_velocity:
                print("Cannot return the flight equivalent velocity as PITOT3 failed to calculate it.")
                print("Will return None.")
                return None
            else:
                return self.flight_equivalent_velocity
        else:
            print("This FacilityState does not have a reference gas state which is required to calculate the flight equivalent velocity.")
            print("None will be returned.")
            return None

    def get_gamma_and_R_string(self):
        """
        Returns a string with the processed gamma and R for printing. I found that I was printing them a lot together,
        so seemed a like a good idea.

        """

        gas_state = self.get_gas_state()

        return "gam = {0:.2f}, R = {1:.2f} J/kg K".format(gas_state.gamma, gas_state.R)

    def get_molecular_mass(self):
        """
        Just a function to get the molecular mass as I always find the stuff hidden in the CEA output hard to find...
        (if it is the CEA gas object the molecular mass on the gas object is empty...

        :return:
        """

        if self.get_gas_state_gmodel_type() == 'CEAGas':

            return self.get_gas_state().ceaSavedData['Mmass']
        else:
            return self.get_gas_state().molecular_mass

    def get_species_massf_dict(self):
        """
        If the gas is a CEAGas, this function will return the CEA mass fractions dictionary.

        :return:
        """

        if self.get_gas_state_gmodel_type() == 'CEAGas':

            return self.get_gas_state().ceaSavedData['massf']

        else:
            print("GasModel is not a CEAGas, so this function isn't useful. Will return None.")
            return None

    def get_species_moles_dict(self):
        """
        This is like the function above, but it converts the massf dictionary to mole fractions using the
        molecular weights of all of the species which we should have in self.species_MW_dict

        :return:
        """

        if self.get_gas_state_gmodel_type() == 'CEAGas':

            # start by getting the reduced mass fractions

            species_massf_dict = self.get_species_massf_dict()

            species_moles_dict = {}

            gas_total_MW = self.get_molecular_mass()*1000.0 # to get into g/mol

            for species in species_massf_dict:
                mass_fraction = species_massf_dict[species]
                species_MW = self.species_MW_dict[species]

                mole_fraction = mass_fraction*(gas_total_MW/species_MW)

                species_moles_dict[species] = mole_fraction

            return species_moles_dict

    def get_reduced_species_massf_dict(self):
        """
        If the gas is a CEAGas, it will go through the gas dictionary and remove any values which are empty,
        which is very easy for printing things like fill states which probably only have 1-3 gases out of the full set of
        high temperature possibilities.

        Obviously not useful if the GasModel isn't a CEAGas...

        """

        if self.get_gas_state_gmodel_type() == 'CEAGas':
            # a fill state will not have many different species, so we should just the species which actually exist
            species_massf_dict = self.get_species_massf_dict()
            reduced_species_massf_dict = {}
            for species in species_massf_dict.keys():
                if species_massf_dict[species] > 0.0:
                    reduced_species_massf_dict[species] = species_massf_dict[species]

            return reduced_species_massf_dict

        else:
            print("GasModel is not a CEAGas, so this function isn't useful. Will return None.")
            return None

    def get_reduced_species_moles_dict(self):
        """
        This is like the function above, but for mole fractions.

        :return:
        """

        if self.get_gas_state_gmodel_type() == 'CEAGas':

            # # start by getting the reduced mass fractions
            #
            # reduced_species_massf_dict = self.get_reduced_species_massf_dict()
            #
            # reduced_species_moles_dict = {}
            #
            # gas_total_MW = self.get_molecular_mass()*1000.0 # to get into g/mol
            #
            # for species in reduced_species_massf_dict:
            #     mass_fraction = reduced_species_massf_dict[species]
            #     species_MW = self.species_MW_dict[species]
            #
            #     mole_fraction = mass_fraction*(gas_total_MW/species_MW)
            #
            #     reduced_species_moles_dict[species] = mole_fraction
            #
            # return reduced_species_moles_dict

            if self.get_gas_state_gmodel_type() == 'CEAGas':
                # a fill state will not have many different species, so we should just the species which actually exist
                species_moles_dict = self.get_species_moles_dict()
                reduced_species_moles_dict = {}
                for species in species_moles_dict.keys():
                    if species_moles_dict[species] > 0.0:
                        reduced_species_moles_dict[species] = species_moles_dict[species]

                return reduced_species_moles_dict

        else:
            print("GasModel is not a CEAGas, so this function isn't useful. Will return None.")
            return None

    def get_reduced_species_moles_dict_for_printing(self, number_of_digits = 3):
        """

        This takes the output from the function above and turns it into strings of particular length to make it look better printing.

        I actually then just make it a complete string so I can control what it looks like better...

        :return:
        """

        original_reduced_species_moles_dict = self.get_reduced_species_moles_dict()

        str_reduced_species_moles_dict = {}

        for species in original_reduced_species_moles_dict:
            mole_fraction = original_reduced_species_moles_dict[species]

            if mole_fraction > 0.001:
                mole_fraction_string = f'{mole_fraction:.{number_of_digits}}'
            else:
                mole_fraction_string = f'{mole_fraction:.{number_of_digits}e}'

            str_reduced_species_moles_dict[species] = mole_fraction_string

        # now I turn it into a string for better output...

        output_string = ''

        for species in str_reduced_species_moles_dict:
            if species == list(str_reduced_species_moles_dict.keys())[0]:
                output_string += f"{{'{species}': {str_reduced_species_moles_dict[species]}, "
            elif species == list(str_reduced_species_moles_dict.keys())[-1]:
                output_string += f"'{species}': {str_reduced_species_moles_dict[species]}}}"
            else:
                output_string += f"'{species}': {str_reduced_species_moles_dict[species]}, "

        return output_string

    def get_reduced_composition(self):
        """
        Depending on the outputUnits selected, this function will run the appropriate function above
        and return the dictionary.

        :return:

        """

        if self.outputUnits == 'massf':
            return self.get_reduced_species_massf_dict()
        elif self.outputUnits == 'moles':
            return self.get_reduced_species_moles_dict()

    def get_reduced_composition_for_printing(self):
        """
        Depending on the outputUnits selected, this function will run the appropriate function above
        and return the dictionary.

        :return:

        """

        if self.outputUnits == 'massf':
            return self.get_reduced_species_massf_dict()
        elif self.outputUnits == 'moles':
            return self.get_reduced_species_moles_dict_for_printing()

    def get_reduced_composition_with_units(self):
        """
        Depending on the outputUnits selected, this function will run the appropriate function above
        and return the dictionary with the units too.

        :return:

        """

        return self.get_reduced_composition(), self.outputUnits

    def get_reduced_composition_single_line_output_string(self):
        """
        This is kind of like the two functions above, but it returns a string for printing.

        :return:

        """

        return f'{self.get_reduced_composition_for_printing()} (by {self.outputUnits})'

    def get_reduced_composition_two_line_output_string(self):
        """
        This is kind of like the two functions above, but it returns a string for printing.
        This is the longer two line version.

        :return:

        """

        return f'species in {self.get_state_name()} at equilibrium (by {self.outputUnits}): \n{self.get_reduced_composition_for_printing()}'

    def get_mu(self):
        """
        Function to return the dynamic viscosity mu, as it seems that one needs to the update the transport
        coefficients for it to exist at all...

        Mu is in Pa.s
        :return:
        """

        gas_state = self.get_gas_state()

        if gas_state.mu == 0.0: # this seems to be the default value before it is set. We only set it if needed, as we may be trying to load a pickle object and not want to use CEA again...
            gas_state.update_trans_coeffs()

        return gas_state.mu

    def get_unit_Reynolds_number(self):
        """
        Function to return the unit Reynolds number.
        :return:
        """

        # this is easy, just (rho*v)/mu

        gas_state = self.get_gas_state()
        v = self.get_v()
        mu = self.get_mu()

        return (gas_state.rho*v)/mu

    def get_dictionary_output(self, add_trans_coeffs_and_unit_Re = False):
        """
        This function puts everything useful in the state into a dictionary for when we might want to return it in this form...

        do_not_add_transport_coefficients input is used for the state10f state which has bogus transport coefficients...

        :return:
        """

        state_name = self.state_name

        gas_state = self.get_gas_state()

        gmodel_type = self.get_gas_state_gmodel_type()

        molecular_mass = self.get_molecular_mass()

        v = self.get_v()
        M = self.get_M()

        Ht = self.get_total_enthalpy()
        Ue = self.get_flight_equivalent_velocity()

        total_p = self.get_total_pressure()
        total_T = self.get_total_temperature()
        pitot_p = self.get_pitot_pressure()

        output_dict = {'state_name':state_name, 'gmodel_type':gmodel_type,
                       'rho':gas_state.rho, 'p':gas_state.p, 'T':gas_state.T, 'a':gas_state.a,
                       'h':gas_state.enthalpy, 'u':gas_state.u, 's':gas_state.entropy,
                       'gamma':gas_state.gamma, 'R':gas_state.R, 'Cp':gas_state.Cp, 'molecular_mass':molecular_mass,
                       'v':v, 'M':M, 'Ht':Ht, 'Ue':Ue, 'total_p':total_p, 'total_T':total_T, 'pitot_p':pitot_p}

        if add_trans_coeffs_and_unit_Re:
            output_dict['unit_Re'] = self.get_unit_Reynolds_number()
            output_dict['k'] = gas_state.k
            output_dict['mu'] = gas_state.mu

        if gmodel_type == 'CEAGas':
            output_dict['composition_massf'] = self.get_species_massf_dict()

            if self.outputUnits == 'moles':
                output_dict['composition_moles'] = self.get_species_moles_dict()

        return output_dict

class Tube(object):
    """
    This is the class for 'tube' sections in PITOT3, which covers the shock tube, secondary driver tube, and acceleration tube.
    It takes an inlet state which is driving the shock, a fill state, and will calculate the expanding test gas state,
    the shocked fill gas state, and the shock speed.

    """

    # TO DO: these variables could be in a better order, as things like the tube length and diameters are not actually
    # necessary for a lot of things...

    def __init__(self, tube_name, tube_length, tube_diameter, fill_pressure, fill_temperature, fill_gas_model, fill_gas_name, fill_gas_filename,
                 fill_state_name, shocked_fill_state_name, entrance_state_name, entrance_state, unsteadily_expanded_entrance_state_name,
                 expand_to, expansion_factor,
                 preset_gas_models_folder, unsteady_expansion_steps, vs_guess_1, vs_guess_2, vs_limits, vs_tolerance,
                 outputUnits = 'massf', species_MW_dict = None):

        self.tube_name = tube_name

        # give the shock a name here too...
        if self.tube_name == 'secondary_driver':
            self.shock_name = 'vsd'
        elif self.tube_name == 'shock_tube':
            self.shock_name = 'vs1'
        elif self.tube_name == 'acceleration_tube':
            self.shock_name = 'vs2'
        else:
            # just call it vs...
            self.shock_name = 'vs'

        self.tube_length = tube_length
        self.tube_diameter = tube_diameter

        self.fill_pressure = fill_pressure
        self.fill_temperature = fill_temperature

        self.fill_gas_model = fill_gas_model
        self.fill_gas_name = fill_gas_name
        self.fill_gas_filename = fill_gas_filename
        self.fill_state_name = fill_state_name

        if self.fill_gas_model == 'CEAGas' and self.fill_gas_name:
            fill_gmodel_location = '{0}/cea-{1}-gas-model.lua'.format(preset_gas_models_folder, self.fill_gas_name)
            # we create this link and then use if it exists...
            fill_room_temperature_only_gmodel_location = '{0}/cea-{1}-room-temperature-only-gas-model.lua'.format(preset_gas_models_folder, self.fill_gas_name)
        elif self.fill_gas_model == 'custom' and self.fill_gas_filename:
            fill_gmodel_location = self.fill_gas_filename
            # we make it the gas name with -room-temperature-only-gas-model added... (in the same folder)
            fill_room_temperature_only_gmodel_location = fill_gmodel_location.split('.')[0] + '-room-temperature-only-gas-model.' + fill_gmodel_location.split('.')[-1]
        else:
            fill_room_temperature_only_gmodel_location = None

        fill_gmodel = GasModel(os.path.expandvars(fill_gmodel_location))

        if fill_room_temperature_only_gmodel_location and os.path.exists(os.path.expandvars(fill_room_temperature_only_gmodel_location)):

            # if there is a room temperature only object we "trick" the fill state gas model to use that when it sets the gas
            # state and then replace the correct gas state after. This seemed to be the best / cleanest way to do it.

            fill_room_temperature_only_gmodel = GasModel(os.path.expandvars(fill_room_temperature_only_gmodel_location))

            fill_state_gas_object = GasState(fill_gmodel)

            fill_state_gas_object.gmodel = fill_room_temperature_only_gmodel

            fill_state_gas_object.p = self.fill_pressure
            fill_state_gas_object.T = self.fill_temperature

            fill_state_gas_object.update_thermo_from_pT()
            fill_state_gas_object.update_sound_speed()

            fill_state_gas_object.gmodel = fill_gmodel

        else:

            fill_room_temperature_only_gmodel = None

            fill_state_gas_object = GasState(fill_gmodel)

            fill_state_gas_object.p = self.fill_pressure
            fill_state_gas_object.T = self.fill_temperature

            fill_state_gas_object.update_thermo_from_pT()
            fill_state_gas_object.update_sound_speed()

        # fill state isn't moving...
        fill_state_v = 0.0 #m/s

        # now make the related facility state object:
        # the fill state can be its own reference state...
        self.fill_state = Facility_State(self.fill_state_name, fill_state_gas_object, fill_state_v,
                                         reference_gas_state=fill_state_gas_object, room_temperature_only_gmodel = fill_room_temperature_only_gmodel,
                                         outputUnits = outputUnits, species_MW_dict = species_MW_dict)

        self.shocked_fill_state_name = shocked_fill_state_name

        self.entrance_state_name = entrance_state_name
        self.entrance_state = entrance_state

        # I thought that we may as well have two names for these objects, as we want to follow the entrance / exit state
        # convention of the rest of the code, but also have a sensible internal name in this class...
        self.unsteadily_expanding_state_name = self.entrance_state_name
        self.unsteadily_expanding_state = self.entrance_state

        self.unsteadily_expanded_entrance_state_name = unsteadily_expanded_entrance_state_name

        self.expand_to = expand_to

        self.unsteady_expansion_steps = unsteady_expansion_steps

        self.vs_guess_1 = vs_guess_1
        self.vs_guess_2 = vs_guess_2
        self.vs_limits = vs_limits
        self.vs_tolerance = vs_tolerance

        return

    def calculate_shock_speed(self):
        """
        Function to do the shock speed calculation for a given fill state and unsteadily expanding driver state.

        """

        # the function below sets it up, code is below...

        def error_in_velocity_function(vs, fill_state = self.fill_state, unsteadily_expanding_state = self.unsteadily_expanding_state,
                                       steps=self.unsteady_expansion_steps, shock_name = self.shock_name):
            """Compute the velocity mismatch for a given shock speed."""

            print ('-' * 60)
            print (f"Current guess for {shock_name} = {vs:.2f} m/s")

            # the shocked state has the same gmodel as the fill_state
            # which we can use to get the shocked state and its gas flow object...
            fill_state_gmodel = fill_state.get_gas_state().gmodel

            shocked_gas_state = GasState(fill_state_gmodel)
            fill_state_gas_flow = GasFlow(fill_state_gmodel)

            # do the shock

            v2, v2g = fill_state_gas_flow.normal_shock(fill_state.get_gas_state(), vs, shocked_gas_state)

            # now set up and do the unsteady expansion

            # Across the contact surface, p3 == p2 (i.e. the post-shock pressure is teh same as the unsteadily expanded pressure)
            p3 = shocked_gas_state.p

            unsteadily_expanding_state_gmodel = unsteadily_expanding_state.get_gas_state().gmodel

            unsteadily_expanded_gas_state = GasState(unsteadily_expanding_state_gmodel)
            unsteadily_expanding_state_gas_flow = GasFlow(unsteadily_expanding_state_gmodel)

            v3g = finite_wave_dp_wrapper(unsteadily_expanding_state.get_gas_state(), unsteadily_expanding_state.get_v(),
                                         'cplus', p3, unsteadily_expanded_gas_state, unsteadily_expanding_state_gas_flow,
                                         steps = steps, gmodel_without_ions=unsteadily_expanding_state.get_gas_state_gmodel_without_ions())

            # shocked state label number will be the fill state label + 1, the unsteadily expanded state will be that number + 2
            # this is easy for shock tube and acceleration tube as we have s1 or s5, but harder for secondary driver where it is sd1
            # so we'll do something different based on the length...
            fill_state_name = fill_state.get_state_name()

            if len(fill_state_name) == 2:
                shocked_state_label_number = int(fill_state_name[1]) + 1
                unsteadily_expanded_state_label_number = int(fill_state_name[1]) + 2
            elif len(fill_state_name) == 3:
                shocked_state_label_number = 'sd' + str(int(fill_state_name[2]) + 1)
                unsteadily_expanded_state_label_number = 'sd' + str(int(fill_state_name[2]) + 2)

            print(f"Current p{shocked_state_label_number} = {shocked_gas_state.p:.2f} Pa, current p{unsteadily_expanded_state_label_number} = {unsteadily_expanded_gas_state.p:.2f} Pa.")
            print(f"Current v{shocked_state_label_number}g = {v2g:.2f} m/s, current v{unsteadily_expanded_state_label_number}g = {v3g:.2f} m/s.")
            if abs((v2g - v3g) / v2g) > 0.001:
                print(f"Current (v{shocked_state_label_number}g - v{unsteadily_expanded_state_label_number}g) / v{shocked_state_label_number}g = {(v2g - v3g) / v2g:.6f}.")
            else:
                print(f"Current (v{shocked_state_label_number}g - v{unsteadily_expanded_state_label_number}g) / v{shocked_state_label_number}g = {(v2g - v3g) / v2g:.3e}.")

            return (v2g - v3g) / v2g

        print('-'*60)
        print(f"Calculating {self.tube_name} shock speed ({self.shock_name}).")
        print(f"{self.tube_name} fill state is:")
        print(self.fill_state)
        print("Unsteadily expanding entry state is:")
        print(self.unsteadily_expanding_state)

        # calculate the shock speed using a secant solver
        self.vs = secant(error_in_velocity_function, self.vs_guess_1, self.vs_guess_2, limits = self.vs_limits, tol=self.vs_tolerance)

        self.Ms = self.vs / self.fill_state.get_gas_state().a

        print ('-' * 60)
        print (f"From secant solve: {self.shock_name} = {self.vs:.2f} m/s")

        return

    def calculate_shock_speed_and_related_states(self):
        """
        This just runs the calculate shock speed and state calculating functions.
        :return:
        """

        self.calculate_shock_speed()
        self.calculate_shocked_and_unsteadily_expanded_states()

    def set_shock_speed(self, vs):
        """
        Function to set the shock speed. Basically for use with the experimental mode where the shock speed is forced.

        :param vs: shock speed in m/s

        """

        print ('-'*60)
        print(f"Setting {self.tube_name} shock speed ({ self.shock_name}) to {vs:.2f} m/s")

        self.vs = vs
        self.Ms = self.vs / self.fill_state.get_gas_state().a

        return

    def set_shock_speed_and_calculate_related_states(self, vs):
        """
        Function to set the shock speed and then calculate the related states...
        :param vs:
        :return:
        """

        self.set_shock_speed(vs)
        self.calculate_shocked_and_unsteadily_expanded_states()

        return

    def calculate_shocked_and_unsteadily_expanded_states(self):

        """
        After the shock speed is added (either calculated or set by the user), this allows the shocked and unsteadily
        expanded states to be found.

        :return:
        """

        if hasattr(self, 'vs'):
            print ('-'*60)
            print(f"Now that {self.shock_name} is known, finding conditions at states {self.shocked_fill_state_name} and {self.unsteadily_expanded_entrance_state_name}.")

            # now that we have the shock speed, we can get the final versions of the post-shock and unsteadily expanded states
            # (we don't use the values from the secant solver as we may want to change what we expand the unsteadily expanded state to)

            fill_state_gmodel = self.fill_state.get_gas_state().gmodel

            shocked_gas_state = GasState(fill_state_gmodel)
            fill_state_gas_flow = GasFlow(fill_state_gmodel)

            v2, v2g = fill_state_gas_flow.normal_shock(self.fill_state.get_gas_state(), self.vs, shocked_gas_state)

            # for the shocked state we can just use the fill state as the reference state...
            # we need to carry around the room temperature only gas model for if we need it for CO2 gases in the nozzle expansion...
            if self.fill_state.get_room_temperature_only_gmodel():
                self.shocked_state = Facility_State(self.shocked_fill_state_name, shocked_gas_state, v2g,
                                                    reference_gas_state=self.fill_state.get_gas_state(),
                                                    room_temperature_only_gmodel=self.fill_state.get_room_temperature_only_gmodel(),
                                                    outputUnits = self.fill_state.outputUnits, species_MW_dict = self.fill_state.species_MW_dict)
            else:
                self.shocked_state = Facility_State(self.shocked_fill_state_name, shocked_gas_state, v2g,
                                                    reference_gas_state=self.fill_state.get_gas_state(),
                                                    outputUnits = self.fill_state.outputUnits, species_MW_dict = self.fill_state.species_MW_dict)

            # Across the contact surface, p3 == p2 (i.e. the post-shock pressure is teh same as the unsteadily expanded pressure)
            p3 = shocked_gas_state.p

            # TO DO: do I keep this as entrance state here or return it to the unsteadily expanding state like above in the function?

            unsteadily_expanding_state_gmodel = self.unsteadily_expanding_state.get_gas_state().gmodel
            unsteadily_expanding_state_gas_flow = GasFlow(unsteadily_expanding_state_gmodel)

            unsteadily_expanded_gas_state = GasState(unsteadily_expanding_state_gmodel)

            # TO DO: need to get finite wave dv working so we can get the expansion factor stuff in here...
            # need something like below...
            #State 7 is being expanded to V6 (8050.30994903) multiplied by an expansion factor of 1.0.

            if self.expand_to == 'flow_behind_shock':
                # just do the same calculation as the secant solver, this is the ideal case

                v3g = finite_wave_dp_wrapper(self.unsteadily_expanding_state.get_gas_state(),
                                             self.unsteadily_expanding_state.get_v(),
                                             'cplus', p3, unsteadily_expanded_gas_state,
                                             unsteadily_expanding_state_gas_flow,
                                             steps=self.unsteady_expansion_steps,
                                             gmodel_without_ions=self.unsteadily_expanding_state.get_gas_state_gmodel_without_ions())

            elif self.expand_to == 'shock_speed':
                # we use finite wave_dv and expand to the shock speed instead...
                v3g = finite_wave_dv_wrapper(self.unsteadily_expanding_state.get_gas_state(),
                                             self.unsteadily_expanding_state.get_v(),
                                             'cplus', self.vs, unsteadily_expanded_gas_state,
                                             unsteadily_expanding_state_gas_flow,
                                             steps=self.unsteady_expansion_steps,
                                             gmodel_without_ions=self.unsteadily_expanding_state.get_gas_state_gmodel_without_ions())

            # for this state, if the entrance state had a reference gas state, we can grab it as the reference state...

            if self.unsteadily_expanding_state.reference_gas_state:
                reference_gas_state = self.unsteadily_expanding_state.get_reference_gas_state()
            else:
                reference_gas_state = None

            # if we had the gas model without ions, we need to ensure that the new facility state also will make a gas model
            # without ions, as it will ahve been set without ions, so the FacilityState constructor will miss it...

            if self.unsteadily_expanding_state.get_gas_state_gmodel_without_ions():
                make_gmodel_without_ions = True
            else:
                make_gmodel_without_ions = None

            # we need to carry around the room temperature only gas model for if we need it for CO2 gases in the nozzle expansion...
            if self.unsteadily_expanding_state.get_room_temperature_only_gmodel():
                self.unsteadily_expanded_state = Facility_State(self.unsteadily_expanded_entrance_state_name, unsteadily_expanded_gas_state, v3g,
                                                                reference_gas_state=reference_gas_state,
                                                                room_temperature_only_gmodel=self.unsteadily_expanding_state.get_room_temperature_only_gmodel(),
                                                                outputUnits = self.unsteadily_expanding_state.outputUnits, species_MW_dict = self.unsteadily_expanding_state.species_MW_dict,
                                                                make_gmodel_without_ions = make_gmodel_without_ions)
            else:
                self.unsteadily_expanded_state = Facility_State(self.unsteadily_expanded_entrance_state_name,
                                                                unsteadily_expanded_gas_state, v3g,
                                                                reference_gas_state=reference_gas_state,
                                                                outputUnits = self.unsteadily_expanding_state.outputUnits, species_MW_dict = self.unsteadily_expanding_state.species_MW_dict,
                                                                make_gmodel_without_ions = make_gmodel_without_ions)

            for facility_state in [self.shocked_state, self.unsteadily_expanded_state]:
                print('-'*60)
                print(facility_state)
                if facility_state.get_gas_state_gmodel_type() == 'CEAGas':
                    print(facility_state.get_reduced_composition_two_line_output_string())

                # add the stagnation enthalpy, if we can:
                if facility_state.reference_gas_state:
                    total_enthalpy = facility_state.get_total_enthalpy()

                    if total_enthalpy:
                        print (f"The total enthalpy (Ht) at state {facility_state.get_state_name()} is {total_enthalpy/1.0e6:.2f} MJ/kg.")
                    else:
                        print(f"Was unable to calculate the total enthalpy at state {facility_state.get_state_name()} so it cannot be printed")
            return

        else:
            # TO DO: should this be an exception?
            print("Shock speed must be either calculated or specified for this function to be used.")
            return

    def calculate_rst_stagnation_conditions(self):
        """
        This function will perform reflected shocks on the post-shock and unsteadily expanded states (if they exist)
        to get the stagnation conditions for an RST.

        I decided to put this here for an RST, as that stuff is all happening in the shock tube, whereas in an expansion tube,
        a reflected shock at the end of the shock tube is caused by the diaphragm, so is in the diaphragm object...

        TO DO: need to give custom names for these states if we want it...

        :return:
        """

        if hasattr(self, 'shocked_state') and hasattr(self, 'unsteadily_expanded_state'):

            # do the stagnated fill state, which is the main gas state which we want
            shocked_gas_gas_state = self.get_shocked_state().get_gas_state()
            shocked_gas_v = self.get_shocked_state().get_v()
            shocked_state_gmodel = shocked_gas_gas_state.gmodel

            shocked_state_gas_flow = GasFlow(shocked_state_gmodel)
            stagnated_fill_gas_state = GasState(shocked_state_gmodel) # remember taht the shocked gas is the fill state...

            self.vr = shocked_state_gas_flow.reflected_shock(shocked_gas_gas_state, shocked_gas_v, stagnated_fill_gas_state)

            self.Mr = (shocked_gas_v + self.vr)/ shocked_gas_gas_state.a #normally this would be V2 - Vr, but it's plus here as Vr has been left positive

            self.stagnated_fill_gas = Facility_State('s5', stagnated_fill_gas_state, 0.0,
                                                     reference_gas_state=self.get_shocked_state().get_reference_gas_state())

            # then do the same thing for the unsteadily expanded driver gas as this is important for RSTs as well...

            unsteadily_expanded_gas_state = self.get_unsteadily_expanded_state().get_gas_state()
            unsteadily_expanded_state_v = self.get_unsteadily_expanded_state().get_v()
            unsteadily_expanded_state_gmodel = unsteadily_expanded_gas_state.gmodel

            unsteadily_expanded_state_gas_flow = GasFlow(unsteadily_expanded_state_gmodel)
            stagnated_unsteadily_expanded_gas_state = GasState(unsteadily_expanded_state_gmodel)

            self.vrd = unsteadily_expanded_state_gas_flow.reflected_shock(unsteadily_expanded_gas_state, unsteadily_expanded_state_v,
                                                                          stagnated_unsteadily_expanded_gas_state)

            self.Mrd = (unsteadily_expanded_state_v + self.vrd) / unsteadily_expanded_gas_state.a  # normally this would be V3 - Vr, but it's plus here as Vr has been left positive

            self.stagnated_unsteadily_expanding_gas = Facility_State('s5d', stagnated_unsteadily_expanded_gas_state, 0.0,
                                                                     reference_gas_state=self.get_unsteadily_expanded_state().get_reference_gas_state())

            for facility_state in [self.stagnated_fill_gas, self.stagnated_unsteadily_expanding_gas]:
                print('-'*60)
                print(facility_state)
                if facility_state.get_gas_state_gmodel_type() == 'CEAGas':
                    print(facility_state.get_reduced_composition_two_line_output_string())

                # add the stagnation enthalpy, if we can:
                if facility_state.reference_gas_state:
                    total_enthalpy = facility_state.get_total_enthalpy()

                    print (f"The total enthalpy (Ht) at state {facility_state.get_state_name()} is {total_enthalpy/1.0e6:.2f} MJ/kg.")

        else:
            print("Shocked and unsteadily expanded states must have been calculated for this function to be used.")
            return

    def get_tube_name(self):
        """
        Return the tube name

        :return:
        """

        return self.tube_name

    def get_shock_speed(self):
        """
        Return the shock speed
        """

        if hasattr(self, 'vs'):
            return self.vs
        else:
            print("Shock speed has not yet been calculated.")

    def get_shock_Mach_number(self):
        """
        Return the shock Mach number
        """

        if hasattr(self, 'Ms'):
            return self.Ms
        else:
            print("Shock speed has not yet been calculated.")

    def get_fill_gas_model(self):
        """
        Return the fill gas model

        """

        return self.fill_gas_model

    def get_fill_gas_name(self):
        """
        Retrun the fill gas name
        :return:
        """

        return self.fill_gas_name

    def get_fill_gas_filename(self):
        """
        Return the fill gas filename
        :return:
        """

        return self.fill_gas_filename

    def get_fill_state(self):
        """
        Return the fill state

        """

        return self.fill_state

    def get_shocked_state(self):
        """
        Return the shocked state
        :return:
        """
        if hasattr(self, 'shocked_state'):
            return self.shocked_state
        else:
            print("Shocked state has not yet been calculated.")

    def get_unsteadily_expanded_state(self):
        """
        Return the unsteadily expanded state
        :return:
        """
        if hasattr(self, 'get_unsteadily_expanded_state'):
            return self.unsteadily_expanded_state
        else:
            print("The unsteadily expanded state has not yet been calculated.")

    def get_exit_state_name(self):
        """
        Get the exit state name. This is the shocked state for the shock and secondary driver tubes
        and the unsteadily expanded state for the acceleration tube.
        :return:
        """
        if self.tube_name == 'acceleration_tube':
            if hasattr(self, 'unsteadily_expanded_state'):
                return self.unsteadily_expanded_state.get_state_name()
            else:
                print("The exit state (the unsteadily expanded state) has not yet been calculated.")
                return None

        else:
            if hasattr(self, 'shocked_state'):
                if hasattr(self, 'stagnated_fill_gas'):
                    # it is a reflected shock tunnel so the exit state is actually the stagnated fill state...
                    return self.stagnated_fill_gas.get_state_name()
                else:
                    # just a normal shock tube, so the shocked state continues on...
                    return self.shocked_state.get_state_name()
            else:
                print("The exit state (the shocked state) has not yet been calculated.")
                return None

    def get_exit_state(self):
        """
        Get the exit state. This is the shocked state for the shock and secondary driver tubes
        and the unsteadily expanded state for the acceleration tube.
        :return:
        """
        if self.tube_name == 'acceleration_tube':
            if hasattr(self, 'unsteadily_expanded_state'):
                return self.unsteadily_expanded_state
            else:
                print("The exit state (the unsteadily expanded state) has not yet been calculated.")
                return None
        else:
            if hasattr(self, 'shocked_state'):
                if hasattr(self, 'stagnated_fill_gas'):
                    # it is a reflected shock tunnel so the exit state is actually the stagnated fill state...
                    return self.stagnated_fill_gas
                else:
                    # just a normal shock tube, so the shocked state continues on...
                    return self.shocked_state
            else:
                print("The exit state (the shocked state) has not yet been calculated.")
                return None

    def get_unsteadily_expanded_state(self):
        """
        Return the unsteadily expanded state
        :return:
        """
        if hasattr(self, 'unsteadily_expanded_state'):
            return self.unsteadily_expanded_state
        else:
            print("Unsteadily expanded state has not yet been calculated.")

    def get_facility_states(self):
        """
        This function returns a list of the unique facility states (i.e. it shouldn't return anything
        which would be returned using the function with the same name on another object) which this object contains,
        in an order which is relevant for the PITOT3 text output.

        This was mainly made for the output, but may have other uses.

        """

        facility_states = [self.fill_state, self.shocked_state, self.unsteadily_expanded_state]

        # then add the reflected shock states if we have them...

        if hasattr(self, 'stagnated_fill_gas') and hasattr(self, 'stagnated_unsteadily_expanding_gas'):
            facility_states += [self.stagnated_fill_gas, self.stagnated_unsteadily_expanding_gas]

        return facility_states

    def get_tube_length_and_diameter(self):
        """
        Return the tube length and diameter if they exist
        :return:
        """

        if self.tube_length and self.tube_diameter:
            return self.tube_length, self.tube_diameter
        else:
            print("This function requires the tube length and diameter to exist. None will be returned.")
            return None

    def get_tube_length(self):
        """
        Return the tube length if it exists
        :return:
        """

        if self.tube_length:
            return self.tube_length
        else:
            print("This function requires the tube length to exist. None will be returned.")
            return None

    def get_tube_diameter(self):
        """
        Return the tube diameter if it exists
        :return:
        """

        if self.tube_diameter:
            return self.tube_diameter
        else:
            print("This function requires the tube diameter to exist. None will be returned.")
            return None

class Nozzle(object):
    """
    This is the class for nozzles in PITOT3. It currently only works for supersonic nozzles as I'm still unsure if I'll
    make it do the expansion to the throat Mach number like PITOT used to do for RSTs, or if I'll put that elsewhere...

    """

    def __init__(self, entrance_state_name, entrance_state, exit_state_name, area_ratio, nozzle_expansion_tolerance,
                 facility_type = None, expansion_tube_nozzle_expansion_minimum_p2_over_p1 = None,
                 maximum_temp_for_room_temperature_only_gmodel = 1100.0, cutoff_temp_for_no_ions = 5000.0):
        """

        :param entrance_state_name:
        :param entrance_state:
        :param exit_state_name:
        :param area_ratio:
        """

        self.entrance_state_name = entrance_state_name
        self.entrance_state = entrance_state

        self.exit_state_name = exit_state_name

        self.area_ratio = area_ratio

        self.facility_type = facility_type

        # TO DO: maybe add some comments about what it is doing...

        print('-'*60)
        print(f"Starting steady expansion through the nozzle using an area ratio of {self.area_ratio}")

        # make the gas model object

        entrance_state_gmodel = self.entrance_state.get_gas_state().gmodel
        entrance_state_gas_flow_object = GasFlow(entrance_state_gmodel)

        exit_gas_state = GasState(entrance_state_gmodel)

        if facility_type == 'reflected_shock_tunnel':
            print("Due to the fact that this is a reflected shock simulation, we must first expand to the throat condition (state 6).")

            from gdtk.ideal_gas_flow import p0_p

            state6 = GasState(entrance_state_gmodel)

            entrance_gas_state = self.entrance_state.get_gas_state()

            # setting it to 1.01 as we want it to be above 1 so it supersonic!
            # (the p0/p is ideal gas).
            # We also have some code to catch up if it ends up subsonic below too
            v6 = entrance_state_gas_flow_object.expand_from_stagnation(entrance_gas_state,
                                                                       1.0 / p0_p(1.01, entrance_gas_state.gamma),
                                                                       state6)

            M6 = v6 / state6.a

            print (f"M6 expected is 1.0 and M6 found is {M6:.4f}.")

            if M6 < 1.0 and M6 >= 0.98:
                print ("M6 is just below 1.0 so it is being set to 1.0 so the code can continue.")
                # this must be above 1 to do the supersonic steady_flow_with_area_change...
                import copy
                old_velocity = copy.copy(v6)

                while M6 < 1.0:
                    v6 += 0.001
                    M6 = v6 / state6.a
                print (f"The velocity has also been slightly adjusted from {old_velocity:.4f} m/s to {v6:.4f} m/s.")
            elif M6 < 0.98:
                raise Exception("pitot3_classes.Nozzle: M throat is too small, so there must be an issue here...")

            # if the entrance state has a reference gas state, we can grab that as the exit state will have the same one.
            if self.entrance_state.reference_gas_state:
                reference_gas_state = self.entrance_state.get_reference_gas_state()
            else:
                reference_gas_state = None
            self.state6 = Facility_State('s6', state6, v6,
                                         reference_gas_state=reference_gas_state,
                                         outputUnits=self.entrance_state.outputUnits, species_MW_dict=self.entrance_state.species_MW_dict)

            v_exit = entrance_state_gas_flow_object.steady_flow_with_area_change(self.state6.get_gas_state(),
                                                                                 self.state6.get_v(),
                                                                                 self.area_ratio, exit_gas_state,
                                                                                 tol=nozzle_expansion_tolerance)


        elif facility_type == 'expansion_tube' and expansion_tube_nozzle_expansion_minimum_p2_over_p1:

            # for expansion tube cases, I reduced the default p2p1_min to 0.01 as most expansion tube nozzles do not have large area ratios.
            # 0.01 was a value which should work with large nozzles such as the X3 Mach 12 nozzle
            # (this now stored in the default config so the user can change it if they want).
            # the main case that fails is CO2 cases, everything else seems to be able to cope with it.

            try:

                v_exit = entrance_state_gas_flow_object.steady_flow_with_area_change(self.entrance_state.get_gas_state(),
                                                                                     self.entrance_state.get_v(),
                                                                                     self.area_ratio, exit_gas_state,
                                                                                     tol=nozzle_expansion_tolerance,
                                                                                     p2p1_min=expansion_tube_nozzle_expansion_minimum_p2_over_p1)
            except Exception as e:
                print(e)
                print("Nozzle expansion has failed.")

                # we have a case for using the room temperature only gs model for CO2
                # and turning off ionisation for cases with ionisaton...

                if self.entrance_state.get_room_temperature_only_gmodel() and self.entrance_state.get_gas_state().T < maximum_temp_for_room_temperature_only_gmodel:
                    print(f"Our gas state has a room temperature only gas model and we are below the maximum temperature for using that model of {maximum_temp_for_room_temperature_only_gmodel}.")
                    print('So we are going to try doing the nozzle expansion with that model.')

                    original_gmodel = self.entrance_state.get_gas_state_gmodel()

                    room_temperature_only_gmodel = self.entrance_state.get_room_temperature_only_gmodel()
                    room_temperature_only_gas_flow = GasFlow(room_temperature_only_gmodel)

                    entrance_gas_state = self.entrance_state.get_gas_state()

                    entrance_gas_state.gmodel = room_temperature_only_gmodel

                    exit_gas_state = GasState(room_temperature_only_gmodel)

                    v_exit = room_temperature_only_gas_flow.steady_flow_with_area_change(
                        entrance_gas_state,
                        self.entrance_state.get_v(),
                        self.area_ratio, exit_gas_state,
                        tol=nozzle_expansion_tolerance,
                        p2p1_min=expansion_tube_nozzle_expansion_minimum_p2_over_p1)

                    entrance_gas_state.gmodel = original_gmodel
                    exit_gas_state.gmodel = original_gmodel

                elif self.entrance_state.get_gas_state_gmodel_without_ions() and self.entrance_state.get_gas_state().T < cutoff_temp_for_no_ions:
                    print(
                        f"We are below the set cutoff temperature of {cutoff_temp_for_no_ions} K so we are going to try performing the nozzle expansion without ions.")

                    original_gmodel = self.entrance_state.get_gas_state_gmodel()

                    gmodel_without_ions = self.entrance_state.get_gas_state_gmodel_without_ions()
                    gas_flow_without_ions = GasFlow(gmodel_without_ions )

                    entrance_gas_state = self.entrance_state.get_gas_state()

                    entrance_gas_state.gmodel = gmodel_without_ions

                    exit_gas_state = GasState(gmodel_without_ions)

                    v_exit = gas_flow_without_ions.steady_flow_with_area_change(
                        entrance_gas_state,
                        self.entrance_state.get_v(),
                        self.area_ratio, exit_gas_state,
                        tol=nozzle_expansion_tolerance,
                        p2p1_min=expansion_tube_nozzle_expansion_minimum_p2_over_p1)

                    entrance_gas_state.gmodel = original_gmodel
                    exit_gas_state.gmodel = original_gmodel

        else:
            v_exit = entrance_state_gas_flow_object.steady_flow_with_area_change(self.entrance_state.get_gas_state(), self.entrance_state.get_v(),
                                                                                 self.area_ratio, exit_gas_state,
                                                                                 tol=nozzle_expansion_tolerance)

        # if the entrance state has a reference gas state, we can grab that as the exit state will have the same one.
        if self.entrance_state.reference_gas_state:
            reference_gas_state = self.entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None
        self.exit_state = Facility_State(self.exit_state_name, exit_gas_state, v_exit,
                                         reference_gas_state=reference_gas_state,
                                         outputUnits=self.entrance_state.outputUnits, species_MW_dict=self.entrance_state.species_MW_dict)

        print (self.exit_state)
        if self.exit_state.get_gas_state_gmodel_type() == 'CEAGas':
            print(self.exit_state.get_reduced_composition_two_line_output_string())

        return

    def get_area_ratio(self):
        """
        Return the area ratio
        :return:
        """

        return self.area_ratio

    def get_entrance_state_name(self):
        """
        Return the name of the nozzle entrance state.
        """

        return self.entrance_state_name

    def get_entrance_state(self):
        """
        Return the nozzle entrance state.
        """

        return self.entrance_state

    def get_exit_state_name(self):
        """
        Return the name of the nozzle exit state.
        """

        return self.exit_state_name

    def get_exit_state(self):
        """
        Return the nozzle exit state.
        """

        return self.exit_state

    def get_facility_states(self):
        """
        This function returns a list of the unique facility states (i.e. it shouldn't return anything
        which would be returned using the function with the same name on another object) which this object contains,
        in an order which is relevant for the PITOT3 text output.

        This was mainly made for the output, but may have other uses.

        """

        if self.facility_type == 'reflected_shock_tunnel':
            return [self.state6, self.exit_state]
        else:
            return [self.exit_state]

class Test_Section(object):
    """
    This is the class for the test section in PITOT3. Just holds the test section state, but can then also be used
    to do calculations over the model.

    """

    def __init__(self, entrance_state_name, entrance_state, test_section_post_shock_state_name):
        """

        :param entrance_state_name:
        :param entrance_state:
        """

        self.entrance_state_name = entrance_state_name
        self.entrance_state = entrance_state

        # remember that the test section state is just the entrance state here...
        self.test_section_state_name = self.entrance_state_name
        self.test_section_state = self.entrance_state

        self.test_section_post_shock_state_name = test_section_post_shock_state_name

        return

    def get_entrance_state_name(self):
        """
        Return the name of the test section entrance state.
        """

        return self.entrance_state_name

    def get_entrance_state(self):
        """
        Return the test section entrance state.
        """
        return self.entrance_state


    def get_test_section_state_name(self):
        """
        Return the name of the test section state
        """

        return self.test_section_state_name


    def get_test_section_state(self):
        """
        Return the test section state.
        """
        return self.test_section_state

    def calculate_post_normal_shock_state(self):
        """
        Calculate the post-normal shock state over the test model.

        :return:
        """

        if self.test_section_state.get_gas_state_gmodel_type() == 'CEAGas':

            print('-'*60)
            print("Starting frozen normal shock calculation over the test model.")

            ideal_gas_gmodel_filename = eilmer4_IdealGas_gas_model_creator(self.test_section_state.get_gas_state())

            test_section_state_ideal_gas_gmodel = GasModel(ideal_gas_gmodel_filename)
            test_section_state_ideal_gas_gas_flow_object = GasFlow(test_section_state_ideal_gas_gmodel)

            test_section_ideal_gas_gas_state = GasState(test_section_state_ideal_gas_gmodel)

            test_section_ideal_gas_gas_state.p = self.test_section_state.get_gas_state().p
            test_section_ideal_gas_gas_state.T = self.test_section_state.get_gas_state().T

            test_section_ideal_gas_gas_state.update_thermo_from_pT()
            test_section_ideal_gas_gas_state.update_sound_speed()

            post_normal_shock_ideal_gas_gas_state = GasState(test_section_state_ideal_gas_gmodel)

            v10f, v10g = test_section_state_ideal_gas_gas_flow_object.normal_shock(test_section_ideal_gas_gas_state,
                                                                                   self.test_section_state.get_v(),
                                                                                   post_normal_shock_ideal_gas_gas_state)

            #v10 here instead of v10g as the model is stationary in the lab frame...
            self.post_normal_shock_ideal_gas_state = Facility_State(f'{self.test_section_post_shock_state_name}f',
                                                                    post_normal_shock_ideal_gas_gas_state, v10f)

            print (self.post_normal_shock_ideal_gas_state)

        print('-'*60)
        print("Starting equilibrium normal shock calculation over the test model.")

        test_section_state_gmodel = self.test_section_state.get_gas_state().gmodel
        test_section_state_gas_flow_object = GasFlow(test_section_state_gmodel)

        post_normal_shock_gas_state = GasState(test_section_state_gmodel)

        v10, v10g = test_section_state_gas_flow_object.normal_shock(self.test_section_state.get_gas_state(),
                                                                    self.test_section_state.get_v(),
                                                                    post_normal_shock_gas_state)

        # if the entrance state has a reference gas state, we can grab that as the post-shock state will have the same one.
        if self.entrance_state.reference_gas_state:
            reference_gas_state = self.entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None

        #v10 here instead of v10g as the model is stationary in the lab frame...
        self.post_normal_shock_state = Facility_State(f'{self.test_section_post_shock_state_name}e',
                                                      post_normal_shock_gas_state, v10,
                                                      reference_gas_state=reference_gas_state,
                                                      outputUnits=self.entrance_state.outputUnits,
                                                      species_MW_dict=self.entrance_state.species_MW_dict)

        print (self.post_normal_shock_state)

        if self.post_normal_shock_state.get_gas_state_gmodel_type() == 'CEAGas':
            print(self.post_normal_shock_state.get_reduced_composition_two_line_output_string())

    def get_post_normal_shock_state(self):
        """
        Return the test section post-normal shock state. (and calculate it if it doesn't exist already.)
        """

        if hasattr(self, 'post_normal_shock_state'):
            return self.post_normal_shock_state
        else:
            # calculate it and then return it
            self.calculate_post_normal_shock_state()
            return self.post_normal_shock_state

    def get_post_normal_shock_ideal_gas_state(self):
        """
        Return the test section post-normal shock state. (and calculate it if it doesn't exist already.)
        """

        if hasattr(self, 'post_normal_shock_ideal_gas_state'):
            return self.post_normal_shock_ideal_gas_state
        else:
            if self.test_section_state.get_gas_state_gmodel_type() == 'CEAGas':
                # calculate it and then return it
                self.calculate_post_normal_shock_state()
                return self.post_normal_shock_ideal_gas_state
            else:
                print("The gas model is not CEAGas, so there will not be a post-normal shock ideal gas state. Returning None.")
                return None

    def calculate_post_conical_shock_state(self, cone_half_angle_degrees):
        """
        Calculate the post-conical shock state for a given cone half angle (in degrees)

        :param cone_half_angle_degrees:
        :return:
        """

        self.cone_half_angle_degrees = cone_half_angle_degrees
        cone_half_angle = math.radians(self.cone_half_angle_degrees)

        print('-' * 60)
        print(f"Starting equilibrium conical shock calculation with a cone half angle of {self.cone_half_angle_degrees} degrees.")

        test_section_state_gmodel = self.test_section_state.get_gas_state().gmodel
        test_section_state_gas_flow_object = GasFlow(test_section_state_gmodel)

        post_conical_shock_gas_state = GasState(test_section_state_gmodel)

        print("Test section freestream state is:")
        print(self.test_section_state)

        # start by getting the shock angle (beta)
        beta = test_section_state_gas_flow_object.beta_cone(self.test_section_state.get_gas_state(),
                                                            self.test_section_state.get_v(),
                                                            cone_half_angle)

        self.conical_shock_half_angle_degrees = math.degrees(beta)

        print(f"Shock angle over the cone is {self.conical_shock_half_angle_degrees} degrees")

        cone_half_angle_calculated, v10c = test_section_state_gas_flow_object.theta_cone(self.test_section_state.get_gas_state(),
                                                                                   self.test_section_state.get_v(),
                                                                                   beta,
                                                                                   post_conical_shock_gas_state)

        cone_half_angle_calculated_degrees = math.degrees(cone_half_angle_calculated)

        print (f"The calculated cone half-angle should be the same as the specified one: {cone_half_angle_calculated_degrees} deg = {self.cone_half_angle_degrees} deg")

        # if the entrance state has a reference gas state, we can grab that as the post-shock state will have the same one.
        if self.entrance_state.reference_gas_state:
            reference_gas_state = self.entrance_state.get_reference_gas_state()
        else:
            reference_gas_state = None

        self.post_conical_shock_state = Facility_State(f'{self.test_section_post_shock_state_name}c',
                                                       post_conical_shock_gas_state, v10c,
                                                       reference_gas_state=reference_gas_state,
                                                       outputUnits=self.entrance_state.outputUnits,
                                                       species_MW_dict=self.entrance_state.species_MW_dict)

        print(self.post_conical_shock_state)

        if self.post_conical_shock_state.get_gas_state_gmodel_type() == 'CEAGas':
            print(self.post_conical_shock_state.get_reduced_composition_two_line_output_string())

        return

    def calculate_post_wedge_shock_state(self, wedge_angle_degrees):
        """
        Calculate the post-wedge shock state for a given wedge angle (in degrees)

        :param wedge_angle_degrees:
        :return:
        """

        self.wedge_angle_degrees = wedge_angle_degrees
        wedge_angle = math.radians(self.wedge_angle_degrees)

        print('-' * 60)
        print(f"Starting equilibrium wedge shock calculation with a wedge angle of {self.wedge_angle_degrees} degrees.")

        print("Test section freestream state is:")
        print(self.test_section_state)

        test_section_state_gmodel = self.test_section_state.get_gas_state().gmodel
        test_section_state_gas_flow_object = GasFlow(test_section_state_gmodel)

        post_wedge_shock_gas_state = GasState(test_section_state_gmodel)

        try:

            # start by getting the shock angle (beta)
            beta = test_section_state_gas_flow_object.beta_oblique(self.test_section_state.get_gas_state(),
                                                                   self.test_section_state.get_v(),
                                                                   wedge_angle)

            self.wedge_shock_angle_degrees = math.degrees(beta)

            print(f"Shock angle over the wedge is {self.wedge_shock_angle_degrees} degrees")

            wedge_angle_calculated, v10w = test_section_state_gas_flow_object.theta_oblique(self.test_section_state.get_gas_state(),
                                                                                            self.test_section_state.get_v(),
                                                                                            beta,
                                                                                            post_wedge_shock_gas_state)

            wedge_angle_calculated_degrees = math.degrees(wedge_angle_calculated)

            # TO DO: could add the check on this angle here which PITOT had...

            print (f"The calculated wedge angle should be the same as the specified one: {wedge_angle_calculated_degrees} deg = {self.wedge_angle_degrees} deg")

            # if the entrance state has a reference gas state, we can grab that as the post-shock state will have the same one.
            if self.entrance_state.reference_gas_state:
                reference_gas_state = self.entrance_state.get_reference_gas_state()
            else:
                reference_gas_state = None

            self.post_wedge_shock_state = Facility_State(f'{self.test_section_post_shock_state_name}w',
                                                         post_wedge_shock_gas_state, v10w,
                                                         reference_gas_state=reference_gas_state,
                                                         outputUnits=self.entrance_state.outputUnits,
                                                         species_MW_dict=self.entrance_state.species_MW_dict)

            print(self.post_wedge_shock_state)

            if self.post_wedge_shock_state.get_gas_state_gmodel_type() == 'CEAGas':
                print(self.post_wedge_shock_state.get_reduced_composition_two_line_output_string())

        except Exception as e:
            print("Wedge shock calculation has failed:")
            print(e)
            print("Result will not be added to the output.")

            self.post_wedge_shock_state = None

        return

    def get_facility_states(self):
        """
        This function returns a list of the unique facility states (i.e. it shouldn't return anything
        which would be returned using the function with the same name on another object) which this object contains,
        in an order which is relevant for the PITOT3 text output.

        This was mainly made for the output, but may have other uses.

        """

        if self.test_section_state.get_gas_state_gmodel_type() == 'CEAGas':
            facility_state_list = [self.post_normal_shock_ideal_gas_state, self.post_normal_shock_state]
        else:
            facility_state_list = [self.post_normal_shock_state]

        if hasattr(self, 'post_conical_shock_state'):
            facility_state_list += [self.post_conical_shock_state]

        if hasattr(self, 'post_wedge_shock_state') and self.post_wedge_shock_state:
            facility_state_list += [self.post_wedge_shock_state]

        return facility_state_list

def state_output_for_final_output(facility_state):
    """
    This function makes the output line for each state in the final PITOT3 output.

    It just collects all of the variables and then returns them...

    """

    gas_state = facility_state.get_gas_state()

    state_name = facility_state.get_state_name()

    p = gas_state.p
    T = gas_state.T
    a = gas_state.a
    v = facility_state.get_v()
    M = facility_state.get_M()
    rho = gas_state.rho
    # don't try to get the pitot and total state if the gas isn't moving, as some CEA mixtures don't like that.
    if v > 0.0:
        pitot_p = facility_state.get_pitot_pressure()
        p0 = facility_state.get_total_pressure()
    else:
        pitot_p = p
        p0 = p

    # we get these out of order as Ht is just h if the gas isn't moving so we can use h if needed
    if facility_state.reference_gas_state:
        h = facility_state.get_sensible_enthalpy()
    else:
        h = '-'

    if facility_state.reference_gas_state and v > 0.0:
        Ht = facility_state.get_total_enthalpy()
    elif facility_state.reference_gas_state and v == 0.0:
        Ht = h
    else:
        Ht = '-'

    # now we go through and create the output line...
    # with some small changes for different variables in different situations...

    output_line = ''

    output_line += "{0:<6}".format(state_name)
    if p < 1.0e6:
        output_line += "{0:<10.7}".format(p)
    else:
        output_line += "{0:<10.3e}".format(p)
    if T < 1000.0:
        output_line += "{0:<8.2f}".format(T)
    elif T >= 10000.0:
        output_line += "{0:<8.0f}".format(T)
    else:
        output_line += "{0:<8.1f}".format(T)
    output_line += "{0:<6.0f}".format(a)
    if v < 10000.0:
        output_line += "{0:<8.1f}".format(v)
    else:
        output_line += "{0:<8.0f}".format(v)
    output_line += "{0:<6.2f}".format(M)
    if rho < 0.1:
        output_line += "{0:<9.2e}".format(rho)
    else:
        output_line += "{0:<9.5f}".format(rho)
    if pitot_p == None:
        # it failed...
        output_line += "{0:<8}".format("failed")
    elif pitot_p/1000.0 < 10000.0:
        output_line += "{0:<8.1f}".format(pitot_p/1000.0) # to get kPa
    else:
        output_line += "{0:<8.0f}".format(pitot_p/1000.0)  # to get kPa
    if p0 == None:
        # it failed...
        output_line += "{0:<7}".format("failed")
    elif p0/1.0e6 < 1000.0:
        output_line += "{0:<7.2f}".format(p0/1.0e6) # to get MPa
    elif 1000.0 <= p0/1.0e6 < 10000.0:
        output_line += "{0:<7.1f}".format(p0 / 1.0e6)  # to get MPa
    else:
        output_line += "{0:<7.0f}".format(p0 / 1.0e6)  # to get MPa
    if Ht == '-':
        output_line += "{0:<6}".format(Ht)
    else:
        if Ht == None:
            # it failed...
            output_line += "{0:<7}".format("failed")
        elif Ht/1.0e6 < 100.0:
            output_line += "{0:<7.2f}".format(Ht/1.0e6) # to get MJ/kg
        else:
            output_line += "{0:<7.1f}".format(Ht/1.0e6)  # to get MJ/kg
    if h == '-':
        output_line += "{0:<6}".format(h)
    else:
        if h == None:
            output_line += "{0:<5}".format("failed")
        else:
            output_line += "{0:<5.2f}".format(h / 1.0e6)  # to get MJ/kg

    return output_line

def expansion_tube_test_time_calculator(acceleration_tube_tube_object):
    """
    Function which takes the acceleration tube tube object and returns the basic test time if it is possible to get it.

    This is the test time between when the contact surface separating the accelerator and unsteady expanded test gases
    arrives at teh end of the tube and the unsteady expansion which is following it.

    In the future, I will try to add the other form of test time termination which are reflected u+a waves...

    What I have called the "basic test time" is T_usx in:
    A. Paull & R. J. Stalker, "Test Flow Disturbances in an Expansion Tube", J. Fluid Mech. (1992), vol. 245, pp. 493-521 (p499).

    :param acceleration_tube_tube_object:
    :return:
    """

    if acceleration_tube_tube_object.get_tube_length():
        acceleration_tube_length = acceleration_tube_tube_object.get_tube_length()
    else:
        print("This function needs the acceleration tube length to be available. None will be returned.")
        return None

    # technically they could also not have the different states in the acceleration tube object too...

    unsteadily_expanded_test_gas_state = acceleration_tube_tube_object.get_unsteadily_expanded_state()

    unsteadily_expanded_test_gas_velocity = unsteadily_expanded_test_gas_state.get_v()
    unsteadily_expanded_test_gas_sound_speed = unsteadily_expanded_test_gas_state.get_gas_state().a

    # shorter version of variables so that the function below makes a bit more sense...
    x_a = acceleration_tube_length
    v7 = unsteadily_expanded_test_gas_velocity
    a7 = unsteadily_expanded_test_gas_sound_speed

    # this is close to the Paull and Stalker equation, except we have state 7 as the test gas (whereas they used state 5,
    # similar to a reflected shock tunnel) and we have used v instead of u for velocity
    T_usx = (x_a*a7)/(v7*(v7 - a7))

    return T_usx

def pitot3_results_output(config_data, gas_path, object_dict, generate_output_files = True):
    """
    Function which takes the config data, the gas path and the object_dict (with everything defined by name)
    and outputs the results of the program

    :param config_data:
    :param gas_path:
    :param object_dict:
    :return:
    """

    # print some setup stuff, and then go through the states in the gas path...

    print('-' * 60)
    print("Test completed. Printing output now.")
    print('-' * 60)

    # if we have a file to output to, we will go through all of the outputs once to the screen, and then a second time
    # to the file, by changing the location of the print function in the outputs...
    output_list = [sys.stdout]

    # if we have an output file we add it to the list
    # I would normally open files with with, but this is a standard operation so can probably deal with not having it...
    # / the code will be a bit more cleaner that way...

    if 'output_filename' in config_data and generate_output_files:
        if '.txt' not in config_data['output_filename']:
            output_file = open(config_data['output_filename'] + '.txt', "w")
        else:
            output_file = open(config_data['output_filename'], "w")
        output_list.append(output_file)

    for output_stream in output_list:
        # words starting with vowels causing issues below...
        if config_data['facility_type'] in ['expansion_tube']:
            print(f"PITOT3 Version {config_data['VERSION_STRING']} doing an {config_data['facility_type']} calculation.",
                  file=output_stream)
        else:
            print(f"PITOT3 Version {config_data['VERSION_STRING']} doing a {config_data['facility_type']} calculation.",
                  file=output_stream)
        print(f"Calculation mode is '{config_data['mode']}'.", file=output_stream)
        if config_data['facility']:
            print(f"Facility is '{config_data['facility']}'.", file=output_stream)

        driver = object_dict['driver']
        state4 = driver.get_driver_rupture_state()
        if config_data['driver_condition'] != 'custom':
            print(f"Driver condition is '{config_data['driver_condition']}'. Driver gas model is {state4.get_gas_state_gmodel_type()}.",
                  file=output_stream)
        else:
            print(f"Using custom driver condition from the file {config_data['driver_condition_filename']}.",
                  file=output_stream)
            print(f"Driver gas model is {state4.get_gas_state_gmodel_type()}.",
                  file=output_stream)

        if state4.get_gas_state_gmodel_type() == 'CEAGas':
            print(f"Driver gas composition is {state4.get_reduced_composition_single_line_output_string()} ({state4.get_gamma_and_R_string()}).",
                  file=output_stream)
        else:
            print(f"Driver gas {state4.get_gamma_and_R_string()}.", file=output_stream)
        if 'nozzle' in object_dict:
            print(f"Nozzle area ratio is {object_dict['nozzle'].get_area_ratio()}.", file=output_stream)

        if 'secondary_driver' in object_dict:
            secondary_driver = object_dict['secondary_driver']

            secondary_driver_fill_state = secondary_driver.get_fill_state()

            if secondary_driver.get_fill_gas_model() != 'custom':
                print(f"Secondary driver gas ({secondary_driver_fill_state.get_state_name()}) is {secondary_driver.get_fill_gas_name()}. Secondary driver gas model is { secondary_driver_fill_state.get_gas_state_gmodel_type()}.",
                      file=output_stream)
            else:
                print(f'Using custom secondary driver gas from the file {secondary_driver.get_fill_gas_filename()}.',
                      file=output_stream)
                print(f"Secondary driver gas model is {secondary_driver_fill_state.get_gas_state_gmodel_type()}.",
                    file=output_stream)

            if secondary_driver_fill_state.get_gas_state_gmodel_type() == 'CEAGas':
                print(f"Secondary driver gas composition is {secondary_driver_fill_state.get_reduced_composition_single_line_output_string()} ({secondary_driver_fill_state.get_gamma_and_R_string()}).",
                      file=output_stream)
            else:
                print(f"Secondary driver gas {secondary_driver_fill_state.get_gamma_and_R_string()}.",
                      file=output_stream)

        # we need the test gas here, which is the shock tube fill state...

        shock_tube = object_dict['shock_tube']

        shock_tube_fill_state = shock_tube.get_fill_state()

        if shock_tube.get_fill_gas_model() != 'custom':
            print(f"Test gas ({shock_tube_fill_state.get_state_name()}) is {shock_tube.get_fill_gas_name()}. Test gas gas model is {shock_tube_fill_state.get_gas_state_gmodel_type()}.",
                  file=output_stream)
        else:
            print(f'Using custom test gas from the file {shock_tube.get_fill_gas_filename()}.',
                  file=output_stream)
            print(f"Test gas gas model is {shock_tube_fill_state.get_gas_state_gmodel_type()}.",
                  file=output_stream)

        if shock_tube_fill_state.get_gas_state_gmodel_type() == 'CEAGas':

            print(f"Test gas composition is {shock_tube_fill_state.get_reduced_composition_single_line_output_string()} ({shock_tube_fill_state.get_gamma_and_R_string()})",
                  file=output_stream)
        else:
            print(f"Test gas {shock_tube_fill_state.get_gamma_and_R_string()}",
                  file=output_stream)

        # if we have an acceleration tube, we are an expansion tube...
        if 'acceleration_tube' in object_dict:
            acceleration_tube = object_dict['acceleration_tube']

            acceleration_tube_fill_state = acceleration_tube.get_fill_state()

            if acceleration_tube.get_fill_gas_model() != 'custom':
                print(f"Accelerator gas ({acceleration_tube_fill_state.get_state_name()}) is {acceleration_tube.get_fill_gas_name()}. Accelerator gas gas model is {acceleration_tube_fill_state.get_gas_state_gmodel_type()}.",
                      file=output_stream)
            else:
                print(f'Using custom accelerator gas from the file {acceleration_tube.get_fill_gas_filename()}.',
                      file=output_stream)
                print(f"Accelerator gas gas model is {acceleration_tube_fill_state.get_gas_state_gmodel_type()}.",
                    file=output_stream)

            if acceleration_tube_fill_state.get_gas_state_gmodel_type() == 'CEAGas':
                print(f"Accelerator gas composition is {acceleration_tube_fill_state.get_reduced_composition_single_line_output_string()} ({acceleration_tube_fill_state.get_gamma_and_R_string()}).",
                      file=output_stream)
            else:
                print(f"Accelerator gas {acceleration_tube_fill_state.get_gamma_and_R_string()}.",
                      file=output_stream)

        if 'secondary_driver' in locals():
            print(f'vsd = {secondary_driver.get_shock_speed():.2f} m/s, Msd = {secondary_driver.get_shock_Mach_number():.2f}',
                  file=output_stream)

        shock_speed_output = f'vs1 = {shock_tube.get_shock_speed():.2f} m/s, Ms1 = {shock_tube.get_shock_Mach_number():.2f}'

        if 'acceleration_tube' in locals():
            shock_speed_output += f', vs2 = {acceleration_tube.get_shock_speed():.2f} m/s, Ms2 = {acceleration_tube.get_shock_Mach_number():.2f}'

        print(shock_speed_output, file=output_stream)

        # doing the reflected shock at the diaphragm stuff a bit hard to semi-automate it...

        if 'secondary_driver' in locals():
            tube_list = ['secondary_driver','shock_tube']
            diaphragm_list = ['secondary_diaphragm', 'tertiary_diaphragm']
        else:
            tube_list = ['shock_tube']
            diaphragm_list = ['secondary_diaphragm']

        for tube, diaphragm in zip(tube_list, diaphragm_list):
            if object_dict[diaphragm].get_diaphragm_type() in ['reflected_shock', 'velocity_loss_factor_and_reflected_shock']:

                tube_name = object_dict[tube].get_tube_name()
                # this tries to just get 'sd', or 'st' which obviously assumes that the user has kept a similar form...
                tube_name_reduced = f"{tube_name.split('_')[0][0]}{tube_name.split('_')[1][0]}"
                vr, Mr = object_dict[diaphragm].get_vr_Mr()

                print(f"NOTE: a user specified reflected shock was done at the end of the {tube_name}.",
                      file=output_stream)

                print(f"vr-{tube_name_reduced} = {vr:.2f} m/s, Mr-{tube_name_reduced} = {Mr:.2f}",
                      file=output_stream)

        #---------------------------------------------------------------------------------------------------------------
        # doing the main output below...

        key = "{0:6}{1:10}{2:8}{3:6}{4:8}{5:6}{6:9}{7:8}{8:7}{9:7}{10:5}".format("state", "p", "T", "a", "v", "M",
                                                                                 "rho", "pitot_p", "stgn_p", "Ht", "h")
        print(key, file=output_stream)

        units = "{0:6}{1:10}{2:8}{3:6}{4:8}{5:6}{6:9}{7:8}{8:7}{9:7}{10:5}".format("", "Pa", "K", "m/s", "m/s", "",
                                                                                   "kg/m^3", "kPa", "MPa", "MJ/kg",
                                                                                   "MJ/kg", )
        print(units, file=output_stream)

        # now we build up the gas states for the output...
        states_list = []

        for object in gas_path:
            states_list += object.get_facility_states()

        for facility_state in states_list:
            output_line = state_output_for_final_output(facility_state)
            print(output_line, file=output_stream)

        #---------------------------------------------------------------------------------------------------------------
        # now do some final extra outputs at the bottom...
        # start by pulling out our freestream and post-shock test section states which we will need later...
        test_section = object_dict['test_section']

        freestream_state = test_section.get_entrance_state()
        test_section_post_normal_shock_state = test_section.get_post_normal_shock_state()

        freestream_total_temperature = freestream_state.get_total_temperature()

        if not freestream_total_temperature:
            # the nozzle expansion is isenthalpic, so if the freestream total enthalpy failed, we can try the nozzle inlet state...
            if 'nozzle' in object_dict:
                print("The freestream total temperature was not available so we will attempt to get it from the nozzle inlet state.")
                nozzle = object_dict['nozzle']
                nozzle_inlet_state = nozzle.get_entrance_state()

                freestream_total_temperature = nozzle_inlet_state.get_total_temperature()

        freestream_flight_equivalent_velocity = freestream_state.get_flight_equivalent_velocity()

        if not freestream_flight_equivalent_velocity :
            # the nozzle expansion is isenthalpic, so if the freestream total enthalpy failed, we can try the nozzle inlet state...
            if 'nozzle' in object_dict:
                print("The freestream flight equivalent velocity was not available so we will attempt to get it from the nozzle inlet state.")
                nozzle = object_dict['nozzle']
                nozzle_inlet_state = nozzle.get_entrance_state()

                freestream_flight_equivalent_velocity = nozzle_inlet_state.get_flight_equivalent_velocity()

                # if this works we should set the freestream state Ht and Ue to this value too...
                # (I will find a way to make this fix up the table output too...

                freestream_state.total_enthalpy = nozzle_inlet_state.get_total_enthalpy()
                freestream_state.flight_equivalent_velocity = freestream_flight_equivalent_velocity

        if freestream_total_temperature:
            print(f"The freestream ({freestream_state.get_state_name()}) total temperature (Tt) is {freestream_total_temperature:.2f} K.",
                  file=output_stream)

        if freestream_flight_equivalent_velocity:
            print(f"The freestream ({freestream_state.get_state_name()}) flight equivalent velocity (Ue) is {freestream_flight_equivalent_velocity:.2f} m/s.",
                  file=output_stream)

        if 'acceleration_tube' in locals():
            if acceleration_tube.tube_length:
                basic_test_time = expansion_tube_test_time_calculator(acceleration_tube)
                print(f"Basic test time = {basic_test_time * 1.0e6:.2f} microseconds.", file=output_stream)

        if hasattr(test_section, 'post_wedge_shock_state') and test_section.post_wedge_shock_state:
            print(f"Wedge angle was {test_section.wedge_angle_degrees:.2f} degrees. Wedge shock angle (beta) was found to be {test_section.wedge_shock_angle_degrees:.2f} degrees.",
                  file=output_stream)
        elif hasattr(test_section, 'post_wedge_shock_state') and not test_section.post_wedge_shock_state:
            print(f"Calculation with a Wedge angle of {test_section.wedge_angle_degrees:.2f} degrees was attempted but it failed.",
                  file=output_stream)

        print(f"Freestream ({freestream_state.get_state_name()}) unit Reynolds number is {freestream_state.get_unit_Reynolds_number():.2f} /m (related mu is {freestream_state.get_mu():.2e} Pa.s).",
              file=output_stream)

        print(f"Post normal shock equilibrium ({test_section_post_normal_shock_state.get_state_name()}) unit Reynolds number is {test_section_post_normal_shock_state.get_unit_Reynolds_number():.2f} /m (related mu is {test_section_post_normal_shock_state.get_mu():.2e} Pa.s).",
              file=output_stream)

        if freestream_state.get_gas_state_gmodel_type()== 'CEAGas':
            print(f"Species in the freestream state ({freestream_state.get_state_name()}) at equilibrium (by {freestream_state.outputUnits}):",
                  file=output_stream)
            print(freestream_state.get_reduced_species_moles_dict_for_printing(), file=output_stream)

        if test_section.get_post_normal_shock_state().get_gas_state_gmodel_type() == 'CEAGas':
            print(f"Species in the shock layer at equilibrium ({test_section.get_post_normal_shock_state().get_state_name()}) (by {test_section.get_post_normal_shock_state().outputUnits}):",
                  file=output_stream)
            print(test_section.get_post_normal_shock_state().get_reduced_species_moles_dict_for_printing(), file=output_stream)

    # some extra stuff at the bottom here, we make a states dict so we can return it later on...

    states_dict = {}

    for state in states_list:
        state_name = state.get_state_name()

        states_dict[state_name] = state

    if generate_output_files:
        # and we output a cut down version of the states to a json file...

        pitot3_states_dict_json_output_file_creator(states_dict, config_data['output_filename'] + '.json')

        # and the full result to a pickle...

        dict_of_objects = {'config_data':config_data, 'gas_path':gas_path,
                           'object_dict':object_dict, 'states_dict':states_dict}

        pitot3_pickle_output_file_creator(dict_of_objects, config_data['output_filename'] + '.pickle')

        # and we output a one line csv of the output too...
        pitot3_single_line_output_file_creator(config_data, object_dict, states_dict)

    return states_dict


def cleanup_function(cleanup_generated_gas_models = False):
    """Function to clean up temporary files created during the running of the program."""

    import os

    print ("-" * 60)
    print ("Removing temporary files created during the running of the program.")
    print ("-" * 60)

    files_to_remove_list = ['thermo.inp', 'thermo.out', 'thermo.lib', 'tmp.inp', 'tmp.out', 'trans.inp', 'trans.out', 'trans.lib']

    if cleanup_generated_gas_models:
        files_to_remove_list += ['PITOT3_ideal_gas_test_section_gmodel.lua', 'PITOT3_cea_driver_condition.lua']

        # we also need to scan the folder for any gas models without ions...

        files_and_folders_in_directory = os.listdir(os.getcwd())

        for file_or_folder in files_and_folders_in_directory:
            if 'without-ions-gas-model.lua' in file_or_folder:
                files_to_remove_list.append(file_or_folder)

    for file in files_to_remove_list:
        if os.path.isfile(file):
            os.remove(file)

    return
